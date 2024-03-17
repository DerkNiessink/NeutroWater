import numpy as np
import pandas as pd
from typing import Sequence, Any
import multiprocessing
from tqdm.contrib.concurrent import process_map
from dataclasses import dataclass

from neutrons.models.neutrons import Neutrons, Neutron, Vector
from neutrons.models.tank import Tank
from neutrons.models.collisions import Collisions
from neutrons.process.data_processor import (
    TotalProcessor,
    SpectrumProcessor,
    AbsorptionProcessor,
)
from neutrons.process.maxwell_boltzmann import MaxwellBoltzmann


@dataclass
class Parameters:
    """
    Parameters:
    - total_data (Sequence): Total cross section data of the form of
        a pandas DataFrame with two columns containing: energy [eV] and total
        cross section [barns].
    - scattering_data (Sequence): Scattering cross section data of
        the form of a pandas DataFrame with two columns containing:
        energy [eV] and scattering cross section [barns].
    - absorption_data (Sequence): Absorption cross section data of the form
        of a pandas DataFrame with two columns containing: energy [eV] and
        absorption cross section [barns].
    - spectrum_Data (pd.DataFrame): Data for the probability distribution of the
        initial energy spectrum from the neutron source wich should have two columns
        containing: energy [eV] and probability.
    - nNeutrons (int): Number of neutrons to simulate
    - molecule_structure (Sequence): number of atoms in the molecule (Default:
        (2, 1) for H20)
    - nuclei_masses (Sequence): Masses of the nuclei in the molecule, should be
        in the same order as the molecule_structure. (Default: (1, 16) for H20)
    - radius_tank (float): Radius of the tank in meters. Default is 1.
    - height_tank (float): Height of the tank in meters. Default is 1.
    - position_tank (np.ndarray): Position of the tank. Default is [0, 0, 0].
    - xi (float): Logarithmic reduction of neutron energy per collision. Default
        is 0.920 (H20).
    - temperature (float): Temperature [K] of the medium.
    """

    total_data: Sequence[pd.DataFrame]
    scattering_data: Sequence[pd.DataFrame]
    absorption_data: Sequence[pd.DataFrame]
    spectrum_data: pd.DataFrame
    nNeutrons: int
    molecule_structure: Sequence = (2, 1)
    nuclei_masses: Sequence = (1, 16)
    radius_tank: float = 1
    height_tank: float = 1
    position_tank: Vector = np.array([0.0, 0.0, 0.0])
    xi: float = 0.920
    temperature: float = 293


class DiffusingNeutrons:
    """
    Class that simulates multiple neutrons from diffusing in a medium.

    Args:
        p (Parameters): Parameters for the simulation.
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
        self.tank = Tank(p.radius_tank, p.height_tank, p.position_tank, p.xi)

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
        self.collisions = Collisions(masses=p.nuclei_masses)
        self.nuclei_masses = p.nuclei_masses

        # Maxwell-Boltzmann distribution for the thermal energy
        self.mw = MaxwellBoltzmann(T=p.temperature)

    def _random_direction(self) -> Vector:
        """
        Sample a random 3D direction.

        Returns a np.ndarray of N 3D (np.ndarray) vectors.
        """
        vec = np.random.normal(size=3)
        return vec / np.linalg.norm(vec, axis=-1)

    def diffuse(self, nCollisions: int):
        """
        Let the neutrons diffuse in the medium.

        Args:
            nCollisions (int): number of times each neutron collides with an
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
            chunk (Sequence[Neutron]): chunk of neutrons to diffuse
            nCollisions (int): number of times each neutron collides with an atomic
            nucleus.
        """
        np.random.seed()  # Set a new seed for each process
        return [self._diffuse_neutron(neutron, nCollisions) for neutron in chunk]

    def _diffuse_neutron(self, neutron: Neutron, nCollisions: int):
        """
        Diffuse a single neutron in the medium.

        Args:
            neutron (Neutron): Neutron to diffuse
            nCollisions (int): number of times the neutron collides with an atomic
            nucleus.
        """
        direction = self._random_direction()

        for _ in range(nCollisions):

            if not self.tank.inside(neutron.position):
                break

            neutron.travel(self.total_processor.get_mfp(neutron.energy), direction)

            collision_func = (
                self._handle_thermal_collision
                if neutron.energy < 10 * self.kT
                else self._handle_collision
            )

            # Determine the nucleus the neutron collides with and handle the collision
            # accordingly.
            direction, energy_loss_frac = (
                collision_func(neutron, self.nuclei_masses[0])
                if np.random.random() < self.total_processor.get_ratio(neutron.energy)
                else collision_func(neutron, self.nuclei_masses[1])
            )

            neutron.collide(energy_loss_frac)
            # If the neutron is absorbed, break the loop
            if energy_loss_frac == 1:
                break

        return neutron

    def _handle_collision(self, neutron: Neutron, mass: float) -> tuple[Vector, float]:
        """
        Handle a collisions of a neutron with a nucleus.

        neutron (Neutron): Neutron to handle.
        index (int): index of the nucleus in the molecule.
        mass (float): mass of the nucleus.

        Returns a tuple with the new direction and the fraction of energy lost.
        """
        return (
            (np.zeros(3), 1)
            if self._absorbed(neutron, self.nuclei_masses.index(mass))
            else (
                self._get_direction(self.collisions.theta(mass)),
                self.collisions.energy_loss_frac(mass),
            )
        )

    def _handle_thermal_collision(
        self, neutron: Neutron, mass: float
    ) -> tuple[Vector, float]:
        """
        Handle a collisions of a thermal neutron with a nucleus.

        neutron (Neutron): Neutron to handle.
        mass (float): mass of the nucleus.

        Returns a tuple with the new direction and the fraction of energy lost.
        """
        return (
            (np.zeros(3), 1)
            if self._absorbed(neutron, self.nuclei_masses.index(mass))
            else (
                self._random_direction(),
                self.mw.thermal_energy() / neutron.energy,
            )
        )

    def _get_direction(self, theta: float) -> Vector:
        """
        Get a direction in cartesian coordinates.

        theta (float): angle in radians.
        """
        phi = np.random.random() * 2 * np.pi
        return np.array(
            [
                np.sin(phi) * np.cos(theta),
                np.sin(phi) * np.sin(theta),
                np.cos(phi),
            ]
        )

    def _absorbed(self, neutron: Neutron, index: int) -> bool:
        """
        Check if the neutron is absorbed by the nucleus.

        neutron (Neutron): Neutron to check.
        index (int): index of the nucleus in the molecule.
        """
        if (
            np.random.random()
            < self.absorption_processor.get_absorption_rates(neutron.energy)[index]
        ):
            return True
        else:
            return False
