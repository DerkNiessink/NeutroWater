import numpy as np
from typing import Sequence

from neutrons.neutrons import Neutrons
from neutrons.tank import Tank


class DiffusingNeutrons:
    """Class that simulates multiple neutrons diffusing in a medium"""

    def __init__(
        self,
        mean_free_paths: Sequence[float],
        energies: Sequence[float],
        initial_positions: Sequence[np.ndarray[np.float64]],
        initial_energies: Sequence[float],
        radius_tank: float = 1,
        height_tank: float = 1,
        position_tank: np.ndarray[np.float64] = np.array([0.0, 0.0, 0.0]),
        xi: float = 0.920,
    ):
        """
        Parameters:
        - mean_free_paths (list): Mean-free-paths in water as a function of neutron energy
        in meters.
        - energies (list): Neutron energies corresponding to the mean-free-paths.
        - initial_positions (list): Initial positions of the neutrons.
        - initial_energies (list): Initial energies of the neutrons.
        - radius_tank (float): Radius of the tank in meters. Default is 1.
        - height_tank (float): Height of the tank in meters. Default is 1.
        - position_tank (np.ndarray): Position of the tank. Default is [0, 0, 0].
        - xi (float): Logarithmic reduction of neutron energy per collision. Default is
        0.920 (water).
        """
        self.energy_loss_frac = 1 / np.exp(xi)
        self.mean_free_paths = mean_free_paths
        self.energies = energies
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

    def _mf_lookup(self, energy: float):
        """
        Get the closest value in the mean-free-path data, that corresponds
        to the current energy of the neutron.
        """
        return self.mean_free_paths[idx_of_closest(self.energies, energy)]

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
                neutron.travel(self._mf_lookup(neutron.energies[-1]), dir)
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


def idx_of_closest(lst: list, K: float):
    lst = np.asarray(lst)
    return (np.abs(lst - K)).argmin()
