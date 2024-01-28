import numpy as np


class Neutron:
    """
    Class that holds the parameters of the Neutron, providing methods for
    traveling and colliding, which adjust the relevant class parameters.
    """

    def __init__(self, initial_energy: float, initial_position: float):
        self.initial_energy: float = initial_energy
        self.initial_position: np.ndarray[np.float64] = np.array([0.0, 0.0, 0.0])
        self.position: np.ndarray[np.float64] = np.copy(initial_position)
        self.energy: float = initial_energy
        self.positions: list = [initial_position]
        self.energies: list = [initial_energy]

    def travel(self, distance: float, direction: np.ndarray[np.float64]):
        """
        Update the position by letting the neutron travel in a given direction
        for a given distance.
        """
        self.position += distance * direction
        self.positions.append(np.copy(self.position))

    def collide(self, energy_loss_frac: float):
        """
        Update the energy by letting the neutron collide resulting in a given
        energy loss.
        """
        self.energy *= energy_loss_frac
        self.energies.append(self.energy)


class DiffusingNeutrons:
    """Class that simulates multiple neutrons diffusing in a medium"""

    def __init__(
        self,
        nNeutrons: int,
        nCollisions: int,
        mean_free_paths: list,
        energies: list,
        xi: float = 0.920,
    ):
        """
        nNeutrons (int): number of neutrons to simulate.
        nCollisions (int): number of times each neutron collides with an atomic nucleus.
        mean_free_paths (list) / energies (list): mean-free-paths in water as a
            function of neutron energy.
        xi (float): logarithmic reduction of neutron energy per collision.
        """
        self.nNeutrons = nNeutrons
        self.nCollisions = nCollisions
        self.energy_loss_frac = 1 / np.exp(xi)
        self.mean_free_paths = mean_free_paths
        self.energies = energies
        self.neutrons = []
        self._init_neutrons()

    def _init_neutrons(self):
        """
        Initialize the neutrons with an initial energy and position
        """
        for _ in range(self.nNeutrons):
            self.neutrons.append(Neutron(5 * 10**6, np.array([0.0, 0.0, 0.0])))

    def _random_directions(self) -> np.ndarray[np.ndarray[np.float64]]:
        """
        Sample random 3D directions

        Eeturns a np.ndarray of 3D (np.ndarray) vectors
        """
        vec = np.random.randn(self.nCollisions, 3)
        return vec / np.linalg.norm(vec, axis=0)

    def _mf_lookup(self, energy):
        """
        Get the closest value in the mean-free-path data, that corresponds
        to the current energy of the neutron.
        """
        return self.mean_free_paths[idx_of_closest(self.energies, energy)]

    def diffuse(self):
        """Let the neutrons diffuse in the medium"""
        for neutron in self.neutrons:
            directions = self._random_directions()
            for dir in directions:
                neutron.travel(self._mf_lookup(neutron.energy), dir)
                neutron.collide(self.energy_loss_frac)


def idx_of_closest(lst: list, K: float):
    lst = np.asarray(lst)
    return (np.abs(lst - K)).argmin()
