import numpy as np

from neutrons.diffusing_neutrons import DiffusingNeutrons, Vector


class Measurer:
    """Class to measure the properties of the neutrons after the simulation."""

    def __init__(self, sim: DiffusingNeutrons):
        self.neutrons = sim.neutrons
        self.sim = sim

    def positions(self) -> list[list[Vector]]:
        """
        Get the positions of the neutrons.

        Returns a list where each element is a list of length nNeutrons:
        [[np.ndarray, np.ndarray, ...], [np.ndarray, np.ndarray, ...], ...]
        """
        return [neutron.positions for neutron in self.neutrons]

    def energies(self) -> list[list[float]]:
        """
        Get the energies of the neutrons.

        Returns a list where each element is a list of length nNeutrons:
        [[float, float, ...], [float, float, ...], ...]
        """
        return [neutron.energies for neutron in self.neutrons]

    def number_total(self) -> int:
        """
        Get the total number of neutrons.

        Returns an int.
        """
        return len(self.neutrons)

    def number_escaped(self) -> int:
        """
        Get the number of neutrons that escaped the tank.

        Returns an int.
        """
        return sum(
            [
                1
                for neutron in self.neutrons
                if not self.sim.tank.inside(neutron.position)
            ],
            start=0,
        )

    def density(self, r: float, dr: float):
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
                np.linalg.norm(neutron.position) < r + dr
                and np.linalg.norm(neutron.position) > r - dr
            ):
                n += 1

        return n / (4 * np.pi * r**2)

    def number_thermal(self, E=0.2):
        """
        Get the number of neutrons that are thermalized.

        Returns an int.
        """
        n = 0
        for energies in self.energies():
            if energies[-1] < E:
                n += 1
        return n

    def thermalize_positions(self) -> list:
        """
        Get the positions of the neutrons where their energy is less than E for
        the first time.

        Args:
            E (float): energy threshold.
        """

        thermalize_positions = []
        for energies, positions in zip(self.energies(), self.positions()):
            # Get the first element of the list of energies that is less than 10kT
            res = list(filter((lambda val: val < (10 * self.sim.kT)), energies))
            if len(res) > 0:
                index = energies.index(res[0])
                thermalize_positions.append(positions[index])

        return thermalize_positions

    def thermalize_distances(self) -> list:
        """
        Get the distance of the neutrons where their energy is less than E for
        the first time.

        Args:
            E (float): energy threshold.
        """
        return [np.linalg.norm(pos) for pos in self.thermalize_positions()]

    def number_absorbed(self) -> int:
        """
        Get the number of neutrons that were absorbed.

        Returns an int.
        """
        n = 0
        for neutron in self.neutrons:
            if len(neutron.positions) < self.sim.nCollisions and self.sim.tank.inside(
                neutron.position
            ):
                n += 1
        return n

    def absorbed_positions(self) -> list[Vector]:
        """
        Get the positions of the neutrons that were absorbed.

        Returns a list of np.ndarray.
        """
        return [
            neutron.positions[-1]
            for neutron in self.neutrons
            if len(neutron.positions) < self.sim.nCollisions
            and self.sim.tank.inside(neutron.position)
        ]

    def absorbed_distances(self) -> list[np.float64]:
        """
        Get the distances of the neutrons that were absorbed.

        Returns a list of floats.
        """
        return [np.linalg.norm(pos) for pos in self.absorbed_positions()]

    def flux(self, r: float) -> float:
        """
        Get the flux of the neutrons.

        Returns a float.
        """
        count = 0
        for positions in self.positions():
            positions = np.array(positions)
            norms = np.linalg.norm(positions, axis=1)
            # Count the number of neutrons that have been in the shell
            count += np.sum((norms[:-1] < r) & (r < norms[1:])) + np.sum(
                (norms[1:] < r) & (r < norms[:-1])
            )
        return count / (4 * np.pi * r**2)

    def energy_spectrum(self, r: float) -> list[float]:
        """
        Get the energy spectrum of the neutrons at a given radius r.

        Args:
            r (float): radius at which to compute the energy spectrum.

        Returns a list of floats.
        """
        result_energies = []
        for positions, energies in zip(self.positions(), self.energies()):
            positions = np.array(positions)
            energies = np.array(energies)
            norms = np.linalg.norm(positions, axis=1)
            # Create boolean arrays for the conditions
            cond1 = np.concatenate([(norms[:-1] < r) & (r < norms[1:]), [False]])
            cond2 = np.concatenate([[False], (norms[1:] < r) & (r < norms[:-1])])
            # Get the energies of the neutrons that have been in the shell
            result_energies += list(energies[cond1]) + list(energies[cond2])
        return result_energies
