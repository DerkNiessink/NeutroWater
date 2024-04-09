import pandas as pd
import numpy as np
from typing import Sequence


class AngularProcessor:
    """
    Class that holds the angular distribution data for neutrons in water and
    provides a method for sampling the angular distribution for a given energy.
    """

    def __init__(self, data: Sequence[pd.DataFrame], masses: Sequence[float]):
        """
        Args:
            - data (Sequence[pd.DataFrame]): angular distribution data of the form
            of a sequence of pandas DataFrames with the first column energy and the
            rest of the columns the coefficients of the Legendre series.
            - masses (Sequence[float]): list of atomic masses of the nuclei in the medium.
        """

        self.data = data
        self.bounds = (-1, 1)
        self.masses = masses

    def _legendre_series(self, data_coeffs):
        """
        Evaluate the coefficients of the Legendre series. The coefficients are
        multiplied by the factor (2 * (i + 1) + 1) / 2, where i is the index of
        the coefficient.

        Args:
            - coeffs (list): coefficients of the Legendre series in the data.

        Returns: list of coefficients of the Legendre series.
        """
        coeffs = [coef * (2 * (i + 1) + 1) / 2 for i, coef in enumerate(data_coeffs)]
        # Add the 0th coefficient.
        coeffs.insert(0, 0.5)
        return coeffs

    def _distribution_index(self, E: float, j: int):
        """
        Find the index of the angular distribution for a given energy.

        Args:
            - E (float): energy of the neutron in eV.
            - j (int): index of the mass in the list

        Returns: index of the distribution.
        """

        # Find the index of the energy in the data closest to the given energy
        energies = self.data[j].iloc[:, 0].to_list()
        for i, val in enumerate(energies):
            if val > E:
                break

        # Use interpolation factor to decide which index to use
        f = (E - energies[i - 1]) / (energies[i] - energies[i - 1])

        if np.random.random() > f:
            return i - 1
        else:
            return i

    def get_CM_cosines(
        self,
        mass: float,
        E: float,
        num_samples: int,
    ) -> list:
        """
        Sample the spectrum using Monte Carlo sampling.

        Args:
            - mass (float): mass of the nucleus.
            - E (float): energy of the neutron in eV.
            - num_samples (int): number of samples to take.

        Returns: list of sampled cosines of the scattering angle.
        """
        j = self.masses.index(mass)

        # Find the index of the distribution.
        index = self._distribution_index(E, j)

        # Get the coefficients of the Legendre series with the index.
        coeffs = self._legendre_series(self.data[j].iloc[index, :].tolist()[1:])

        # Sample the cosines of the scattering angle.
        sampled_points = []
        while len(sampled_points) < num_samples:
            x = np.random.uniform(*self.bounds)
            y = np.random.uniform(0, 1)

            px = np.polynomial.legendre.legval(x, coeffs)

            if y <= px:
                sampled_points.append(x)

        return sampled_points
