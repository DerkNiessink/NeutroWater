import numpy as np
from dataclasses import dataclass

from neutrowater.models.neutrons import Vector
from neutrowater.models.maxwell_boltzmann import MaxwellBoltzmann


@dataclass
class Collision:
    """
    Class to simulate the collisions of neutrons with atomic nuclei.

    Args:D
        - initial_E (float): initial energy of the neutron.
        - initial_direction (Vector): initial direction of the neutron in the labframe.
        - mass (float): mass of the nucleus.
        - scattering_cosine (float): cosine of the scattering angle in the CM frame.
        - absorption (bool): flag to indicate if the collision is an absorption.
        - thermal (bool): flag to indicate if the collision is thermal.

    """

    initial_E: float
    initial_direction: Vector
    mass: float
    scattering_cosine: float
    absorption: bool
    thermal: bool

    def __post_init__(self):
        """
        Define the initial conditions of the collision.
        """

        self.mw = MaxwellBoltzmann()

        # initial velocity of neutron in the labframe
        self.v_n = self.initial_direction * np.sqrt(self.initial_E)

        # Velocity of the CM frame
        v_cm = self.v_n / (1 + self.mass)

        # Initial velocity of neutron in the CM frame
        V_n = self.v_n - v_cm

        # Resulting velocity of the neutron in the labframe
        self.v_n_new = v_cm + self._direction_CM() * np.linalg.norm(V_n)

    def _direction_CM(self) -> Vector:
        """
        Get the post-collision scattering direction of the neutron in the CM frame.

        Returns: a np.ndarray of the post-collision scattering direction.
        """

        # Scattering angle in the CM frame
        theta = np.arccos(self.scattering_cosine)

        # Random azimuthal angle
        phi = np.random.uniform(0, 2 * np.pi)

        # Compute direction in the CM frame
        u = np.sin(theta) * np.cos(phi)
        v = np.sin(theta) * np.sin(phi)
        w = np.cos(theta)

        direction = np.array([u, v, w])

        # Normalize the direction for numerical stability
        return direction / np.linalg.norm(direction)

    @property
    def scattering_direction(self) -> Vector:
        """
        Get the post-collision scattering direction of the neutron.

        Returns: a np.ndarray of the post-collision scattering direction.
        """
        if self.thermal:
            return random_direction()
        else:
            return self.v_n_new / np.linalg.norm(self.v_n_new)

    @property
    def energy_loss_frac(self) -> float:
        """
        Get the fraction of energy lost in a collision for the current collision.

        Returns: fraction of energy lost.
        """
        if self.absorption:
            return 1.0
        elif self.thermal:
            return self.mw.thermal_energy() / self.initial_E
        else:
            return float(
                np.linalg.norm(self.v_n_new) ** 2 / np.linalg.norm(self.v_n) ** 2
            )


def random_direction() -> Vector:
    """
    Sample a random 3D direction.

    Returns: a np.ndarray of N 3D (np.ndarray) vectors.
    """
    vec = np.random.normal(size=3)
    return vec / np.linalg.norm(vec, axis=-1)
