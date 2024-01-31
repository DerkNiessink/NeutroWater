import numpy as np


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
