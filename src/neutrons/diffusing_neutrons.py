import numpy as np
import pandas as pd
from typing import Sequence

from neutrons.neutrons import Neutrons
from neutrons.tank import Tank
from neutrons.data_processor import CrossSectionProcessor, SpectrumProcessor


class DiffusingNeutrons:
    """
    Class that simulates multiple neutrons from diffusing in a medium.
    """

    def __init__(
        self,
        cross_section_data: Sequence[pd.DataFrame],
        spectrum_data: pd.DataFrame,
        nNeutrons: int,
        molecule_structure: Sequence = (2, 1),
        radius_tank: float = 1,
        height_tank: float = 1,
        position_tank: np.ndarray[np.float64] = np.array([0.0, 0.0, 0.0]),
        xi: float = 0.920,
        temperature: float = 293,
    ):
        """
        Parameters:
        - cross_section_data (Sequence): Cross section data of the form of a pandas
        DataFrame with two columns containing: energy [eV] and cross section [barns]
        - spectrum_Data (pd.DataFrame): Data for the probability distribution of the
        initial energy spectrum from the neutron source wich should have two columns
        containing: energy [eV] and probability.
        - nNeutrons (int): Number of neutrons to simulate
        - molecule_structure (Sequence): number of atoms in the molecule (Default:
        (2, 1) for H20)
        - radius_tank (float): Radius of the tank in meters. Default is 1.
        - height_tank (float): Height of the tank in meters. Default is 1.
        - position_tank (np.ndarray): Position of the tank. Default is [0, 0, 0].
        - xi (float): Logarithmic reduction of neutron energy per collision. Default
        is 0.920 (H20).
        - temperature (float): Temperature [K] of the medium.
        """
        self.kT = 8.617333262145e-5 * temperature
        self.energy_loss_frac = 1 / np.exp(xi)
        self.mol_struc = molecule_structure
        self.tank = Tank(radius_tank, height_tank, position_tank, xi)

        # For interpolating the cross section data and computing the mean-free-path.
        self.cross_processor = CrossSectionProcessor(cross_section_data)

        # Sample from the interpolated energy spectrum
        spectrum_processor = SpectrumProcessor(spectrum_data)
        initial_energies = spectrum_processor.sample(num_samples=nNeutrons)
        # Convert the energies MeV -> eV
        initial_energies = [energy * 10**6 for energy in initial_energies]

        # All neutrons start at the origin
        initial_positions = [np.array([0, 0, 0]) for _ in range(nNeutrons)]
        self.neutrons = Neutrons(initial_energies, initial_positions)

    def _random_directions(self, N: int) -> np.ndarray[np.ndarray[np.float64]]:
        """
        Sample random 3D directions.

        Args:
            N (int): number of vectors to generate.

        Returns a np.ndarray of N  3D (np.ndarray) vectors.
        """
        vecs = np.random.normal(size=(N, 3))
        mags = np.linalg.norm(vecs, axis=-1)
        return vecs / mags[..., np.newaxis]

    def _thermal_energy(self):
        """
        Sample an energy (in eV) from the Maxwell-Boltzmann distribution at a given temperature.

        Args:
            temperature (float): Temperature in Kelvin.

        Returns a float sampled energy in eV.
        """
        v = np.random.normal(loc=np.sqrt(2 * self.kT), size=1)
        return 0.5 * v**2

    def diffuse(self, nCollisions: int):
        """
        Let the neutrons diffuse in the medium

        Args:
            nCollisions (int): number of times each neutron collides with an atomic nucleus.
        """
        for neutron in self.neutrons:
            directions = self._random_directions(nCollisions)
            for dir in directions:
                if not self.tank.inside(neutron.positions[-1]):
                    break
                neutron.travel(self.cross_processor.get_mfp(neutron.energies[-1]), dir)
                neutron.collide(self.energy_loss_frac)

    def get_positions(self) -> list:
        """
        Get the positions of the neutrons.

        Returns a list where each element is a list of length nNeutrons:
        [[np.ndarray, np.ndarray, ...], [np.ndarray, np.ndarray, ...], ...]
        """
        return [neutron.positions for neutron in self.neutrons]

    def get_energies(self) -> list:
        """
        Get the energies of the neutrons.

        Returns a list where each element is a list of length nNeutrons:
        [[float, float, ...], [float, float, ...], ...]
        """
        return [neutron.energies for neutron in self.neutrons]

    def get_number_escaped(self) -> int:
        """
        Get the number of neutrons that escaped the tank.

        Returns an int.
        """
        return sum(
            [
                1
                for neutron in self.neutrons
                if not self.tank.inside(neutron.positions[-1])
            ]
        )
