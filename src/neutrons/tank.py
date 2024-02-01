import numpy as np
from dataclasses import dataclass


@dataclass
class Tank:
    """
    Class that holds the parameters of the tank, providing a method for
    checking if a position is inside the tank.

    radius (float): radius of the tank.
    height (float): height of the tank.
    xi (float): logarithmic reduction of neutron energy per collision
                (depends on the medium in the tank).
    position (np.ndarray): position of the tank, default [0, 0, 0].
    """

    radius: float
    height: float
    position: float
    xi: float

    def __post_init__(self):
        """
        Calculate the volume and energy loss fraction.
        """
        self.volume = np.pi * self.radius**2 * self.height
        self.energy_loss_frac = 1 / np.exp(self.xi)

    def inside(self, position: np.ndarray) -> bool:
        """
        Check if a position is inside the tank.

        position (np.ndarray): position to check.
        """
        return (
            np.linalg.norm(position[:2]) < self.radius
            and abs(position[2]) < self.height / 2
        )
