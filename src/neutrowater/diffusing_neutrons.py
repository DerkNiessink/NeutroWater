import numpy as np
from numpy.random import random
import pandas as pd
from typing import Sequence, Any
import multiprocessing
from tqdm.contrib.concurrent import process_map
from dataclasses import dataclass
from importlib import resources

from neutrowater.models.neutrons import Neutrons, Neutron, Vector
from neutrowater.models.tank import Tank
from neutrowater.models.collisions import Collision, random_direction
from neutrowater.process.data_processor import (
    TotalProcessor,
    SpectrumProcessor,
    AbsorptionProcessor,
)
from neutrowater.process.angular_processor import AngularProcessor
from neutrowater.models.maxwell_boltzmann import MaxwellBoltzmann
from neutrowater import data


@dataclass
class Parameters:
    """
    Parameters:
    - nNeutrons (int): Number of neutrons to simulate
    - molecule_structure (Sequence): number of atoms in the molecule (Default:
        (2, 1) for H20)
    - nuclei_masses (Sequence): Masses of the nuclei in the molecule, should be
        in the same order as the molecule_structure. (Default: (1, 16) for H20)
    - radius_tank (float): Radius of the tank in meters. Default is 1.
    - height_tank (float): Height of the tank in meters. Default is 1.
    - position_tank (np.ndarray): Position of the tank. Default is [0, 0, 0].
    - temperature (float): Temperature [K] of the medium. Default is 293.
    """

    nNeutrons: int
    molecule_structure: Sequence = (2, 1)
    nuclei_masses: Sequence = (1, 16)
    radius_tank: float = 1
    height_tank: float = 1
    position_tank: tuple = (0.0, 0.0, 0.0)
    temperature: float = 293

    def __post_init__(self):
        files = resources.files(data)

        self.total_data: Sequence[pd.DataFrame] = [
            pd.read_csv(files / "h_cross_t.txt", sep=r"\s+"),
            pd.read_csv(files / "o_cross_t.txt", sep=r"\s+"),
        ]
        self.scattering_data: Sequence[pd.DataFrame] = [
            pd.read_csv(files / "h_cross_s.txt", sep=r"\s+"),
            pd.read_csv(files / "o_cross_s.txt", sep=r"\s+"),
        ]
        self.absorption_data: Sequence[pd.DataFrame] = [
            pd.read_csv(files / "h_cross_a.txt", sep=r"\s+"),
            pd.read_csv(files / "o_cross_a.txt", sep=r"\s+"),
        ]
        self.angular_data: Sequence[pd.DataFrame] = [
            pd.read_csv(files / "H_angular.txt", sep=r"\s+"),
            pd.read_csv(files / "O_angular.txt", sep=r","),
        ]
        self.spectrum_data: pd.DataFrame = pd.read_csv(
            files / "neutron_spectrum_normalized.txt", sep=","
        )


class DiffusingNeutrons:
    """
    Class that simulates multiple neutrons from diffusing in a medium.

    Args:
        - p (Parameters): Parameters for the simulation.
    """

    def __init__(
        self,
        p: Parameters,
    ):
        self.kT = (
            1.380649e-23 * p.temperature * 6.24150907 * 10**18
        )  # [m^2 kg / (s^2 K) * K * (eV/J) = J*(eV/J) = eV]
        self.nCollisions = 0
        self.mol_struc = p.molecule_structure
        self.tank = Tank(p.radius_tank, p.height_tank, p.position_tank)

        # For interpolating the total cross section data and computing the
        # mean-free-path.
        self.total_processor = TotalProcessor(p.total_data)

        # For interpolating the scattering and absorption cross section and
        # computing the aborption ratio.
        self.absorption_processor = AbsorptionProcessor(
            p.scattering_data, p.absorption_data
        )

        # Sample from the interpolated energy spectrum
        spectrum_processor = SpectrumProcessor(p.spectrum_data)
        initial_energies = spectrum_processor.sample(num_samples=p.nNeutrons)

        # Convert the energies MeV -> eV
        initial_energies = [energy * 10**6 for energy in initial_energies]

        # All neutrons start at the origin
        initial_positions = [np.array([0, 0, 0]) for _ in range(p.nNeutrons)]
        self.neutrons = Neutrons(initial_energies, initial_positions)

        # For handling the collisions with the atomic nuclei, hydrogen and oxygen.
        self.nuclei_masses = p.nuclei_masses

        # Maxwell-Boltzmann distribution for the thermal energy
        self.mw = MaxwellBoltzmann(T=p.temperature)

        # For handling the angular distribution of the scattering.
        self.angular_processor = AngularProcessor(p.angular_data, p.nuclei_masses)

    def diffuse(self, nCollisions: int):
        """
        Let the neutrons diffuse in the medium.

        Args:
            - nCollisions (int): number of times each neutron collides with an
            atomic nucleus.
        """
        self.nCollisions += nCollisions
        # Using not all cores but 2 less than the total number of cores, tends
        # to be faster
        num_processes = multiprocessing.cpu_count() - 2

        # Split the neutrons into chunks
        chunk_size = (len(self.neutrons) + num_processes - 1) // num_processes
        chunks = [
            self.neutrons[i : i + chunk_size]
            for i in range(0, len(self.neutrons), chunk_size)
        ]

        # Diffuse the chuncks of neutrons in parallel
        with multiprocessing.Pool(processes=num_processes) as pool:
            results = list(
                process_map(self._diffuse_chunk, chunks, [nCollisions] * len(chunks)),
            )

        # Update the neutrons list with the results from the parallel processess
        self.neutrons = [neutron for chunk in results for neutron in chunk]

    def _diffuse_chunk(self, chunk: Neutrons, nCollisions: int):
        """
        Diffuse a chunk of all neutrons in the medium.

        Args:
            - chunk (Sequence[Neutron]): chunk of neutrons to diffuse
            - nCollisions (int): number of times each neutron collides with an atomic
            nucleus.

        Returns: the neutrons after diffusing in the medium.
        """
        np.random.seed()  # Set a new seed for each process
        return [self._diffuse_neutron(neutron, nCollisions) for neutron in chunk]

    def _diffuse_neutron(self, neutron: Neutron, nCollisions: int) -> Neutron:
        """
        Diffuse a single neutron in the medium.

        Args:
            - neutron (Neutron): Neutron to diffuse
            - nCollisions (int): number of times the neutron collides with an atomic
            nucleus.

        Returns: the neutron after diffusing in the medium.
        """
        neutron.direction = random_direction()

        for _ in range(nCollisions):

            if not self.tank.inside(neutron.position):
                break

            # Sample the distance the neutron travels from an exponential distribution
            # with mean free path as the scale parameter.
            l = -self.total_processor.get_mfp(neutron.energy) * np.log(random())
            neutron.travel(l)

            # Determine the nucleus the neutron collides with
            mass = (
                self.nuclei_masses[0]
                if random() < self.total_processor.get_ratio(neutron.energy)
                else self.nuclei_masses[1]
            )

            collision = Collision(
                initial_E=neutron.energy,
                initial_direction=neutron.direction,
                mass=mass,
                scattering_cosine=self.angular_processor.get_CM_cosines(
                    mass, neutron.energy, 1
                )[0],
                absorption=self._absorbed(neutron, self.nuclei_masses.index(mass)),
                thermal=neutron.energy < 10 * self.kT,
            )

            # Update the neutron's energy and direction
            neutron.collide(collision.energy_loss_frac, collision.scattering_direction)

            # If the neutron is absorbed, break the loop
            if collision.energy_loss_frac == 1:
                break

        return neutron

    def _absorbed(self, neutron: Neutron, index: int) -> bool:
        """
        Check if the neutron is absorbed by the nucleus.

        Args:
            - neutron (Neutron): Neutron to check.
            - index (int): index of the nucleus in the molecule.

        Returns: True if the neutron is absorbed, False otherwise.
        """
        if (
            random()
            < self.absorption_processor.get_absorption_rates(neutron.energy)[index]
        ):
            return True
        else:
            return False
