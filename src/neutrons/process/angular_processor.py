import pandas as pd
import numpy as np
from typing import Sequence


class AngularProcessor:

    def __init__(self, data: Sequence[pd.DataFrame]):
        self.data = data
        self.bounds = (-1, 1)

    def _legendre_series(self, coeffs):
        coeffs = [coef * (2 * (i + 1) + 1) / 2 for i, coef in enumerate(coeffs)]
        coeffs.insert(0, 0.5)
        return coeffs

    def _distribution_index(self, E: float, j: int):
        """
        Find the index of the angular distribution for a given energy.
        """
        energies = self.data[j].iloc[:, 0].to_list()
        for i, val in enumerate(energies):
            if val > E:
                break

        f = (E - energies[i - 1]) / (energies[i] - energies[i - 1])

        if np.random.random() < f:
            return i - 1
        else:
            return i

    def get_CM_cosines(
        self,
        A,
        E,
        num_samples: int,
    ) -> list:
        """
        Sample the spectrum using Monte Carlo sampling.

        Parameters:
        - num_samples (int): number of samples to take.
        """
        j = 0 if A == 1 else 1
        index = self._distribution_index(E, j)
        sampled_points = []
        coeffs = self._legendre_series(self.data[j].iloc[index, :].tolist()[1:])

        while len(sampled_points) < num_samples:
            x = np.random.uniform(*self.bounds)
            y = np.random.uniform(0, 1)

            px = np.polynomial.legendre.legval(x, coeffs)

            if y <= px:
                sampled_points.append(x)

        return sampled_points
