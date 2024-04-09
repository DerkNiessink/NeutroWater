import numpy as np
from dataclasses import dataclass

from neutrons.models.neutrons import Vector
from neutrons.models.maxwell_boltzmann import MaxwellBoltzmann


@dataclass
class Collision:
    """
    Class to simulate the collisions of neutrons with atomic nuclei.

    Args:
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
        self.v_cm = self.v_n / (1 + self.mass)

        # Initial velocity of neutron in the CM frame
        self.V_n = self.v_n * (1 - 1 / (1 + self.mass))

        # Scattering direction of the neutron in the CM frame
        self.scattering_direction = (
            random_direction() if self.thermal else self._compute_scattering_direction()
        )

    def _compute_energy_loss_frac(self) -> float:
        """
        Get the fraction of energy lost in a collision for a given scattering direction
        and a nucleus of a given mass.

        Returns: fraction of energy lost.
        """

        # New velocity of neutron in the labframe after collision
        v_n_new = np.linalg.norm(self.V_n) * self.scattering_direction + self.v_cm

        return float(np.linalg.norm(v_n_new) ** 2 / np.linalg.norm(self.v_n) ** 2)

    def _compute_scattering_direction(self) -> Vector:
        """
        Get the post-collision scattering direction of the neutron in the CM frame.

        Returns: a np.ndarray of the post-collision scattering direction.
        """

        # Scattering cosine in the CM frame
        mu = self.scattering_cosine

        # Random azimuthal angle
        phi = np.random.uniform(0, 2 * np.pi)

        # Components of the initial velocity of neutron in the CM frame
        u, v, w = self.V_n / np.linalg.norm(self.V_n)

        # Relate the pre and post neutron directions in the CM frame
        u = mu * u + np.sqrt(1 - mu**2) * (
            u * w * np.cos(phi) - v * np.sin(phi)
        ) / np.sqrt(1 - w**2)
        v = mu * v + np.sqrt(1 - mu**2) * (
            v * w * np.cos(phi) + u * np.sin(phi)
        ) / np.sqrt(1 - w**2)
        w = mu * w - np.sqrt(1 - mu**2) * np.sqrt(1 - w**2) * np.cos(phi)

        direction = np.array([u, v, w])

        return direction / np.linalg.norm(direction)

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
            return self._compute_energy_loss_frac()


def random_direction() -> Vector:
    """
    Sample a random 3D direction.

    Returns: a np.ndarray of N 3D (np.ndarray) vectors.
    """
    vec = np.random.normal(size=3)
    return vec / np.linalg.norm(vec, axis=-1)
