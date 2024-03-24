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
        - initial_direction (Vector): initial direction of the neutron.
        - mass (float): mass of the nucleus.
        - scattering_cosine (float): cosine of the scattering angle.
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
        self.mw = MaxwellBoltzmann()
        if self.thermal:
            self.scattering_direction = random_direction()
        else:
            self.scattering_direction = self._compute_scattering_direction()

    def _compute_energy_loss_frac(self) -> float:
        """
        Get the fraction of energy lost in a collision for a given scattering direction
        and a nucleus of a given mass.
        """

        # Initial velocity of neutron in the labframe
        v_n = self.initial_direction * np.sqrt(2 * self.initial_E)

        # Velocity of the center of mass
        v_cm = v_n / (1 + self.mass)

        # New velocity of neutron in the labframe after collision
        V_n = np.linalg.norm(v_n - v_cm) * self.scattering_direction + v_cm

        return float(np.linalg.norm(V_n) ** 2 / np.linalg.norm(v_n) ** 2)

    def _compute_scattering_direction(self) -> Vector:
        """
        Get the scattering direction of the neutron after a collision.
        """

        # Compute the scattering cosine in the lab frame
        mu = (1 + self.mass * self.scattering_cosine) / np.sqrt(
            self.mass**2 + 2 * self.mass * self.scattering_cosine + 1
        )
        phi = np.random.uniform(0, 2 * np.pi)
        u, v, w = self.initial_direction

        # Relate the pre and post neutron directions
        u = mu * u + np.sqrt(1 - mu**2) * (
            u * w * np.cos(phi) - v * np.sin(phi)
        ) / np.sqrt(1 - w**2)
        v = mu * v + np.sqrt(1 - mu**2) * (
            v * w * np.cos(phi) + u * np.sin(phi)
        ) / np.sqrt(1 - w**2)
        w = mu * w - np.sqrt(1 - mu**2) * np.cos(phi) * np.sqrt(1 - w**2)
        return np.array([u, v, w])

    @property
    def energy_loss_frac(self) -> float:
        """
        Get the fraction of energy lost in a collision for the current collision.
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
