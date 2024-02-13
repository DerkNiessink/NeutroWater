import numpy as np

from neutrons.diffusing_neutrons import DiffusingNeutrons


class Measurer:
    def __init__(self, sim: DiffusingNeutrons):
        self.neutrons = sim.neutrons
        self.sim = sim

    def positions(self) -> list:
        """
        Get the positions of the neutrons.

        Returns a list where each element is a list of length nNeutrons:
        [[np.ndarray, np.ndarray, ...], [np.ndarray, np.ndarray, ...], ...]
        """
        return [neutron.positions for neutron in self.neutrons]

    def energies(self) -> list:
        """
        Get the energies of the neutrons.

        Returns a list where each element is a list of length nNeutrons:
        [[float, float, ...], [float, float, ...], ...]
        """
        return [neutron.energies for neutron in self.neutrons]

    def number_escaped(self) -> int:
        """
        Get the number of neutrons that escaped the tank.

        Returns an int.
        """
        return sum(
            [
                1
                for neutron in self.neutrons
                if not self.sim.tank.inside(neutron.positions[-1])
            ],
            start=0,
        )

    def density(self, r, dr):
        """
        Get the density of the medium at a given radius r.

        Args:
            r (float): radius at which to compute the density.
            dr (float): thickness of the shell.

        Returns a float.
        """
        n = 0
        for neutron in self.neutrons:
            if (
                np.linalg.norm(neutron.positions[-1]) < r + dr
                and np.linalg.norm(neutron.positions[-1]) > r - dr
            ):
                n += 1

        return n / (4 * np.pi * r**2)

    def thermalize_positions(self, E=0.2) -> list:
        """
        Get the positions of the neutrons where their energy is less than E for
        the first time.

        Args:
            E (float): energy threshold.
        """

        thermalize_positions = []
        for energies, positions in zip(self.energies(), self.positions()):
            # Get the first element of the list of energies that is less than 0.2
            res = list(filter((lambda val: val < E), energies))
            if len(res) > 0:
                index = energies.index(res[0])
                thermalize_positions.append(positions[index])

        return thermalize_positions

    def thermalize_distances(self, E=0.2) -> list:
        """
        Get the distance of the neutrons where their energy is less than E for
        the first time.

        Args:
            E (float): energy threshold.
        """

        thermalize_positions = self.thermalize_positions(E)
        return [np.linalg.norm(pos) for pos in thermalize_positions]
