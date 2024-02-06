import pandas as pd
import numpy as np
from scipy.interpolate import Akima1DInterpolator


class NeutronData:
    def __init__(self, O_data: pd.DataFrame, H_data: pd.DataFrame):
        self.H_interpolater = self.interpolate(H_data)
        self.O_interpolater = self.interpolate(O_data)

    def interpolate(self, data: pd.DataFrame):
        """
        Preprocess the data to be used for the mean free path calculations.
        """
        data = data.drop_duplicates(subset="energy(eV)")
        xp = data["energy(eV)"].values
        fp = data["sigma_t(b)"].values
        x_log = np.log(xp)
        y_log = np.log(fp)
        return Akima1DInterpolator(x_log, y_log)

    def get_mfp(self, energy: float) -> float:
        """
        Get the mean free path for a given energy.
        """
        pass
