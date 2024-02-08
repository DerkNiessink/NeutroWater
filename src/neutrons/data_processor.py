import pandas as pd
import numpy as np
from scipy.interpolate import Akima1DInterpolator
from typing import Sequence


class DataProcessor:
    """
    Class that holds the cross section data for neutrons in water and provides a
    method for calculating the mean free path for a given energy. The data is of
    the form of a pandas DataFrame with columns "energy(eV)" and "sigma_t(b)".
    """

    def __init__(self, data: Sequence[pd.DataFrame]) -> None:
        """
        interpolate the data using scipy.interpolate.Akima1DInterpolator.

        Parameters:
        - data (Sequence[pd.DataFrame]): cross section data.
        """
        self.interpolaters = [self.interpolate(d) for d in data]

    @classmethod
    def interpolate(
        self, data: pd.DataFrame, column1="energy(eV)", column2="sigma_t(b)"
    ) -> Akima1DInterpolator:
        """
        Preprocess the data to be used for the mean free path calculations.

        Parameters:
            - data (pd.DataFrame): cross section data of the form of a pandas
            DataFrame with columns "energy(eV)" and "sigma_t(b)"
        """
        data = data.drop_duplicates(subset=column1)
        xp = data[column1].values
        fp = data[column2].values
        x_log = np.log(xp)
        y_log = np.log(fp)
        return Akima1DInterpolator(x_log, y_log)

    @classmethod
    def cross_section(self, f: Akima1DInterpolator, energy: float) -> float:
        """
        Get the cross section for a given energy in m^2.

        Parameters:
        - f (Akima1DInterpolator): interpolater for the cross section data.
        - energy (float): energy of the neutron in eV.
        """
        return np.exp(f(np.log(energy))) * 10 ** (-28)  # convert barns -> m^2

    def get_mfp(
        self,
        energy: float,
        n: float = 3.3327677048150236e28,
        nMolecules: Sequence = (2, 1),
    ) -> float:
        """
        Get the mean free path in meters for a given energy.

        Parameters:
        - energy (float): energy of the neutron in eV.
        - n (float): number of atoms per volume in m^-3 (Default is for H20).
        - nMolecules (Sequence): number of atoms in the molecule (Default:
        (2, 1) for H20)
        """
        cross_sections = [self.cross_section(f, energy) for f in self.interpolaters]
        return 1 / (
            sum([nMolecules[i] * n * cs for i, cs in enumerate(cross_sections)])
        )
