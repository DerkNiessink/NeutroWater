import numpy as np
from typing import Sequence, overload
from collections.abc import Iterator
import typing


Vector = np.ndarray[typing.Any, np.dtype[np.float64]]


class Neutron:
    """
    Class that holds the parameters of the Neutron, providing methods for
    traveling and colliding, which adjust the relevant class parameters.
    """

    def __init__(self, initial_energy: float, initial_position: Vector):
        self.positions: list = [initial_position]
        self.energies: list = [initial_energy]

    @property
    def position(self) -> Vector:
        return self.positions[-1]

    @property
    def energy(self) -> float:
        return self.energies[-1]

    @position.setter
    def position(self, vector: Vector) -> None:
        self.positions.append(vector)

    @energy.setter
    def energy(self, value: float) -> None:
        self.energies.append(value)

    def travel(self, distance: float):
        """
        Update the position by letting the neutron travel in a given direction
        for a given distance.
        """
        self.position = self.position + distance * self.direction

    def collide(self, energy_loss_frac: float, direction: Vector):
        """
        Update the energy by letting the neutron collide resulting in a given
        energy loss.
        """
        self.direction = direction
        self.energy = self.energy * energy_loss_frac


class Neutrons:
    """Class that holds multiple neutrons and their parameters."""

    def __init__(
        self,
        energies: Sequence[float],
        positions: Sequence[Vector],
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

    def __iter__(self) -> Iterator[Neutron]:
        """
        Allows for iterating over the neutrons in the class.
        """
        return iter(self.neutrons)

    def __len__(self) -> int:
        """
        Returns the number of neutrons in the class.
        """
        return len(self.neutrons)

    @overload
    def __getitem__(self, index: int) -> Neutron: ...

    @overload
    def __getitem__(self, index: slice) -> list[Neutron]: ...

    def __getitem__(self, index: int | slice) -> Neutron | list[Neutron]:
        """
        Returns the neutron at a given index.
        """
        return self.neutrons[index]
