import numpy as np
from typing import Sequence


class Neutron:
    """
    Class that holds the parameters of the Neutron, providing methods for
    traveling and colliding, which adjust the relevant class parameters.
    """

    def __init__(self, initial_energy: float, initial_position: np.ndarray[np.float64]):
        self.positions: list = [initial_position]
        self.energies: list = [initial_energy]

    def travel(self, distance: float, direction: np.ndarray[np.float64]):
        """
        Update the position by letting the neutron travel in a given direction
        for a given distance.
        """
        self.positions.append(self.positions[-1] + distance * direction)

    def collide(self, energy_loss_frac: float):
        """
        Update the energy by letting the neutron collide resulting in a given
        energy loss.
        """
        self.energies.append(self.energies[-1] * energy_loss_frac)


class Neutrons:
    """Class that holds multiple neutrons  and their parameters."""

    def __init__(
        self, energies: Sequence[float], positions: Sequence[np.ndarray[np.float64]]
    ):
        """
        energies (Sequence[float]): initial energies of the neutrons, should be
            same length as positions.
        positions (Sequence[np.ndarray[np.float64]]): initial positions of the
            neutrons, should be same length as energies.
        """
        self.energies = energies
        self.positions = positions
        self._init_neutrons()

    def _init_neutrons(self):
        """
        Initialize the neutrons with an initial energy and position
        """
        self.neutrons = [
            Neutron(energy, position)
            for energy, position in zip(self.energies, self.positions)
        ]

    def __iter__(self) -> iter:
        """
        Allows for iterating over the neutrons in the class.
        """
        return iter(self.neutrons)

    def __len__(self) -> int:
        """
        Returns the number of neutrons in the class.
        """
        return len(self.neutrons)

    def __getitem__(self, index: int) -> Neutron:
        """
        Returns the neutron at a given index.
        """
        return self.neutrons[index]
