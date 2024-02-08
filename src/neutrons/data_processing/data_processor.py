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

    def __init__(self, data: Sequence[pd.DataFrame], log: bool = True) -> None:
        """
        interpolate the data using scipy.interpolate.Akima1DInterpolator.

        Parameters:
        - data (Sequence[pd.DataFrame]): cross section data.
        """
        self.interpolaters = [self.interpolate(d, log) for d in data]

    @classmethod
    def interpolate(self, data: pd.DataFrame, log: bool = True) -> Akima1DInterpolator:
        """
        Preprocess the data to be used for the mean free path calculations.

        Parameters:
            - data (pd.DataFrame): cross section data of the form of a pandas
            DataFrame with columns "energy(eV)" and "sigma_t(b)"
        """
        data = data.drop_duplicates(subset=data.columns.values.tolist()[0])
        x = data.iloc[:, 0]
        y = data.iloc[:, 1]
        if log:
            x = np.log(x)
            y = np.log(y)
        return Akima1DInterpolator(x, y)
