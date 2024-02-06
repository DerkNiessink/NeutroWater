import pandas as pd
import numpy as np
from scipy.interpolate import Akima1DInterpolator


class DataProcessor:
    """
    Class that holds the cross section data for neutrons in water and provides a
    method for calculating the mean free path for a given energy.
    """

    # constants for H2O
    rho = 997  # kg/m^3
    M = 0.01801528  # kg/mol
    N_A = 6.02214076 * 10**23  # mol^-1
    n = rho * N_A / M  # kg/m^3 * 1/mol * mol/kg = m^-3

    def __init__(self, O_data: pd.DataFrame, H_data: pd.DataFrame):
        """
        Convert the cross section data from barns to m^2 and interpolate the data.
        """
        self.H_interpolater = self.interpolate(H_data)
        self.O_interpolater = self.interpolate(O_data)

    def interpolate(self, data: pd.DataFrame) -> Akima1DInterpolator:
        """
        Preprocess the data to be used for the mean free path calculations.
        """
        data = data.drop_duplicates(subset="energy(eV)")
        xp = data["energy(eV)"].values
        fp = data["sigma_t(b)"].values
        x_log = np.log(xp)
        y_log = np.log(fp)
        return Akima1DInterpolator(x_log, y_log)

    def cross_section(self, f: Akima1DInterpolator, energy: float) -> float:
        """
        Get the cross section for a given energy in m^2.
        """
        return np.exp(f(np.log(energy))) * 10 ** (-28)  # convert barns -> m^2

    def get_mfp(self, energy: float) -> float:
        """
        Get the mean free path for a given energy.
        """

        cross_section_H = self.cross_section(self.H_interpolater, energy)
        cross_section_O = self.cross_section(self.O_interpolater, energy)

        return 1 / (self.n * (2 * cross_section_H + cross_section_O))
