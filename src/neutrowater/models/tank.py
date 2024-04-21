import numpy as np
from dataclasses import dataclass


@dataclass
class Tank:
    """
    Class that holds the parameters of the tank, providing a method for
    checking if a position is inside the tank.

    Args:
        - radius (float): radius of the tank.
        - height (float): height of the tank.
        - position (tuple): position of the tank, default (0, 0, 0).
    """

    radius: float
    height: float
    position: tuple = (0.0, 0.0, 0.0)

    def __post_init__(self):
        """
        Calculate the volume and energy loss fraction.
        """
        self.volume = np.pi * self.radius**2 * self.height

    def inside(self, position: np.ndarray) -> bool:
        """
        Check if a position is inside the tank.

        Args:
            - position (np.ndarray): position to check.

        Returns: True if the position is inside the tank, False otherwise.
        """
        return bool(
            np.linalg.norm(position[:2] - self.position[:2]) < self.radius
            and abs(position[2] - self.position[2]) < self.height / 2
        )
