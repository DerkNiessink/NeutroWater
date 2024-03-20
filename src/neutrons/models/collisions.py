import numpy as np
from typing import Sequence


class Collisions:
    """
    Class to simulate the collisions of neutrons with atomic nuclei.
    """

    def __init__(self, masses: Sequence[float]):
        """
        Args:
            - masses: list of atomic masses of the nuclei in the medium.
        """
        self.masses = masses
        self.min_energies = [((m - 1) / (m + 1)) ** 2 for m in masses]
        self.constants = [np.sqrt((1 + m) ** 2 / (4 * m)) for m in masses]

    def energy_loss_frac(self, mass: float):
        """
        Get the fraction of energy lost in a collision with a nucleus of mass `mass`.

        Args:
            - mass (float): atomic mass of the nucleus.

        Returns: fraction of energy lost.
        """
        return np.random.uniform(self.min_energies[self.masses.index(mass)], 1)

    def theta(self, mass: float) -> float:
        """
        Get the scattering angle of a neutron with a nucleus of mass `mass`.

        Args:
            - mass (float): atomic mass of the nucleus.

        Returns: scattering angle in radians.
        """
        Esc = self.energy_loss_frac(mass)
        E = Esc if mass == 1 else 1 - Esc
        return np.arccos(self.constants[self.masses.index(mass)] * np.sqrt(E))
