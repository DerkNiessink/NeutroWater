import pandas as pd
import numpy as np
from scipy.interpolate import Akima1DInterpolator
from typing import Sequence


class DataProcessor:
    """
    Class that holds data and provides a method for interpolating. The data is
    of the form of a sequence of pandas DataFrames with columns x and y.

    Attributes:
    - interpolaters (Sequence[Akima1DInterpolator]): interpolaters for the data.
    """

    def __init__(self, data: Sequence[pd.DataFrame], log: bool) -> None:
        """
        interpolate the data using scipy.interpolate.Akima1DInterpolator.

        Parameters:
        - data (Sequence[pd.DataFrame]): cross section data of the form of a
        sequence of pandas DataFrames with columns x and y.
        - log (bool): if True, the data is log-log transformed.
        """
        self.interpolaters = [self.interpolate(d, log) for d in data]

    @classmethod
    def interpolate(self, data: pd.DataFrame, log: bool = False) -> Akima1DInterpolator:
        """
        Preprocess the data to be used for the mean free path calculations.

        Parameters:
            - data (pd.DataFrame): cross section data of the form of a pandas
            DataFrame with columns x and y.
            - log (bool): if True, the data is log-log transformed.
        """
        data = data.drop_duplicates(subset=data.columns.values.tolist()[0])
        x = data.iloc[:, 0]
        y = data.iloc[:, 1]
        if log:
            x = np.log(x)
            y = np.log(y)
        return Akima1DInterpolator(x, y)


class CrossSectionProcessor(DataProcessor):
    """
    Class that holds the cross section data for neutrons in water and provides a
    method for calculating the mean free path for a given energy. The data is of
    the form of a sequence of pandas DataFrames with columns x and y.

    Attributes:
    - interpolaters (Sequence[Akima1DInterpolator]): interpolaters for the data.
    """

    def __init__(self, data: Sequence[pd.DataFrame]):
        """
        Parameters:
        - data (Sequence[pd.DataFrame]): cross section data of the form of a
        sequence of pandas DataFrames with columns x and y.
        """
        super().__init__(data=data, log=True)

    def cross_section(
        self,
        energy: float,
        f: Akima1DInterpolator = None,
    ) -> float:
        """
        Get the cross section for a given energy in m^2.

        Parameters:
        - energy (float): energy of the neutron in eV.
        - f (Akima1DInterpolator): interpolater for the cross section data.
        """
        f = self.interpolaters[0] if f is None else f
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
        cross_sections = [self.cross_section(energy, f) for f in self.interpolaters]
        return 1 / (
            sum([nMolecules[i] * n * cs for i, cs in enumerate(cross_sections)])
        )


class SpectrumProcessor(DataProcessor):
    """
    Class that holds the spectrum data for neutrons and provides a method for
    sampling the spectrum.

    Attributes:
    - interpolaters (Sequence[Akima1DInterpolator]): interpolaters for the data.
    """

    def __init__(self, data: pd.DataFrame):
        """
        Parameters:
        - data (pd.DataFrame): spectrum data of the form of a pandas DataFrame with
        columns x and y.
        - log (bool): if True, the data is log-log transformed.
        """
        x = data.iloc[:, 0]
        self.bounds = (min(x), max(x))
        super().__init__(data=[data], log=False)

    def sample(
        self,
        num_samples: int,
    ) -> list:
        """
        Sample the spectrum using Monte Carlo sampling.

        Parameters:
        - num_samples (int): number of samples to take.
        """
        f = self.interpolaters[0]
        sampled_points = []
        while len(sampled_points) < num_samples:
            x = np.random.uniform(*self.bounds)
            y = np.random.uniform(0, 1)

            px = f(x)

            if y <= px:
                sampled_points.append(x)

        return sampled_points
