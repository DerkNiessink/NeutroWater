import numpy as np


# J -> eV
J = 6.24150907e18  # eV

# Boltzmann constant
kB = 1.3806452e-23


class MaxwellBoltzmann:
    """
    Class to sample from the Maxwell-Boltzmann distribution.
    """

    def __init__(self, m=1.67493e-27, T=293):
        self.m = m
        self.T = T

    def distribution(self, v):
        """
        Maxwell-Boltzmann distribution for velocities.

        Args:
            v (float): velocity
        """
        return (
            (self.m / (2 * np.pi * kB * self.T)) ** 1.5
            * 4
            * np.pi
            * v**2
            * np.exp(-self.m * v**2 / (2 * kB * self.T))
        )

    def thermal_velocity(self) -> float:
        """
        Samples a thermal velocity of the particle.
        """
        return self._sample()[0]

    def thermal_energy(self) -> float:
        """
        Sample a thermal energy of the particle.
        """
        return 0.5 * self.m * self.thermal_velocity() ** 2 * J

    def _sample(
        self,
        x_min=0,
        x_max=10000,
        y_max=0.0004,
        num_samples=1,
    ) -> list:
        """
        Sample from a function using Monte Carlo sampling.

        Args:
            - f (Callable): function to sample from.
            - x_min (float)
            - x_max (float)
            - y_max (float)
            - num_samples (int)
        """
        sampled_points = []
        while len(sampled_points) < num_samples:
            x = np.random.uniform(x_min, x_max)
            y = np.random.uniform(0, y_max)

            px = self.distribution(x)

            if y <= px:
                sampled_points.append(x)

        return sampled_points
