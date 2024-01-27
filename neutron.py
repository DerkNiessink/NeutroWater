import numpy as np
import pandas as pd


def read_mf_data() -> pd.DataFrame:
    """
    Read the mean-free-path data.

    Returns a dataframe containing the energies and corresponding mean-free-
    paths.
    """
    return pd.read_csv("data/mean_free_path.csv")


def closest_value(
    data1: pd.DataFrame, data2: pd.DataFrame, val: float
) -> tuple[float, float]:
    """
    Find the index of a value in data1 closest to val, and return
    the corresponding value in data2.
    """
    df_sort = data1.iloc[(data1 - val).abs().argsort()[:1]]
    i = df_sort.index.tolist()[0]
    return data2[i]


class Neutron:
    """Class that simulates a Neutron travelling through water"""

    def __init__(self, N: int, xi=0.920, E0=5 * 10**6):
        """
        N (float): number of collisions.
        xi (float): mean logarithmic energy reduction, default is for water.
        E0 (float): initial energy of the neutron in eV, default 5 MeV.
        """
        self.N = N
        self.E = E0
        self.mfs = []

        # Sample N random directions.
        self._sample_directions()
        # Average fraction of energy left after each collision.
        self.frac = 1 / np.exp(xi)
        # mean-free-path dataset
        self.data = read_mf_data()

    def walk(self, r0: np.ndarray):
        """
        Let the neutron walk from a given initial position.

        r0 (np.ndarray): initial position of the neutron.
        """
        r = r0
        self.positions = [np.copy(r0)]
        for i in range(self.N):
            # let the neutron walk for the mean free distance in a random
            # direction
            r += self._mf_distance() * self.directions[i]
            # lower the energy for every collision
            self.E = self.E * self.frac
            self.positions.append(np.copy(r))

    def _sample_directions(self):
        """Sample N random directions"""
        vec = np.random.randn(self.N, 3)
        self.directions = vec / np.linalg.norm(vec, axis=0)

    def _mf_distance(self):
        """
        Get the closest value in the mean-free-path data, that corresponds
        to the current energy of the neutron.
        """
        mf = closest_value(self.data["Energy(eV)"], self.data["lambda(m)"], self.E)
        self.mfs.append(mf)
        return mf
