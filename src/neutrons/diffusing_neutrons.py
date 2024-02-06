import numpy as np
import pandas as pd
from typing import Sequence

from neutrons.neutrons import Neutrons
from neutrons.tank import Tank
from neutrons.data_processor import DataProcessor


class DiffusingNeutrons:
    """Class that simulates multiple neutrons diffusing in a medium"""

    def __init__(
        self,
        data: Sequence[pd.DataFrame],
        initial_positions: Sequence[np.ndarray[np.float64]],
        initial_energies: Sequence[float],
        molecule_structure: Sequence = (2, 1),
        radius_tank: float = 1,
        height_tank: float = 1,
        position_tank: np.ndarray[np.float64] = np.array([0.0, 0.0, 0.0]),
        xi: float = 0.920,
    ):
        """
        Parameters:
        - data (Sequence): cross section data of of the form of a pandas DataFrame
        with columns "energy(eV)" and "sigma_t(b)"
        - initial_positions (Sequence): Initial positions of the neutrons.
        - initial_energies (Sequence): Initial energies of the neutrons.
        - molecule_structure (Sequence): number of atoms in the molecule (Default:
        (2, 1) for H20)
        - radius_tank (float): Radius of the tank in meters. Default is 1.
        - height_tank (float): Height of the tank in meters. Default is 1.
        - position_tank (np.ndarray): Position of the tank. Default is [0, 0, 0].
        - xi (float): Logarithmic reduction of neutron energy per collision. Default
        is 0.920 (H20).
        """
        self.energy_loss_frac = 1 / np.exp(xi)
        self.mol_struc = molecule_structure
        self.data_processor = DataProcessor(data)
        self.neutrons = Neutrons(initial_energies, initial_positions)
        self.tank = Tank(radius_tank, height_tank, position_tank, xi)
        self.nNeutrons = len(self.neutrons)

    def _random_directions(self, N: int) -> np.ndarray[np.ndarray[np.float64]]:
        """
        Sample random 3D directions.

        N (int): number of vectors to generate.

        Returns a np.ndarray of N  3D (np.ndarray) vectors.
        """
        vecs = np.random.normal(size=(N, 3))
        mags = np.linalg.norm(vecs, axis=-1)
        return vecs / mags[..., np.newaxis]

    def diffuse(self, nCollisions: int):
        """
        Let the neutrons diffuse in the medium

        nCollisions (int): number of times each neutron collides with an atomic nucleus.
        """
        for neutron in self.neutrons:
            directions = self._random_directions(nCollisions)
            for dir in directions:
                if not self.tank.inside(neutron.positions[-1]):
                    break
                neutron.travel(self.data_processor.get_mfp(neutron.energies[-1]), dir)
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
