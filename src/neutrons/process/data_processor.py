import pandas as pd
import numpy as np
from scipy.interpolate import Akima1DInterpolator
from typing import Sequence
import itertools


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

        Args:
            - data (Sequence[pd.DataFrame]): cross section data of the form of a
            sequence of pandas DataFrames with columns x and y.
            - log (bool): if True, the data is log-log transformed.
        """
        self.interpolaters = [self._interpolate(d, log) for d in data]

    @staticmethod
    def _interpolate(data: pd.DataFrame, log: bool = False) -> Akima1DInterpolator:
        """
        Preprocess the data to be used for the mean free path calculations.

        Args:
            - data (pd.DataFrame): cross section data of the form of a pandas
            DataFrame with columns x and y.
            - log (bool): if True, the data is log-log transformed.
        """
        data = data.drop_duplicates(subset=data.columns.values.tolist()[0])
        data = data.loc[(data != 0.0).all(axis=1)]

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
        Args:
            - data (Sequence[pd.DataFrame]): cross section data of the form of a
            sequence of pandas DataFrames with columns x and y.
        """
        super().__init__(data=data, log=True)

    def cross_section(
        self,
        energy: float,
        f: Akima1DInterpolator | None = None,
    ) -> float:
        """
        Get the cross section for a given energy in m^2.

        Args:
            - energy (float): energy of the neutron in eV.
            - f (Akima1DInterpolator): interpolater for the cross section data.
        """
        f = self.interpolaters[0] if f is None else f
        return np.exp(f(np.log(energy))) * 10 ** (-28)  # convert barns -> m^2


class TotalProcessor(CrossSectionProcessor):
    """
    Class that holds the total cross section data for neutrons and provides a
    method for calculating the mean free path for a given energy.
    """

    def get_mfp(
        self,
        energy: float,
        n: float = 3.3327677048150236e28,
        nMolecules: Sequence = (2, 1),
    ) -> float:
        """
        Get the mean free path in meters for a given energy.

        Args:
            - energy (float): energy of the neutron in eV.
            - n (float): number of atoms per volume in m^-3 (Default is for H20).
            - nMolecules (Sequence): number of atoms in the molecule (Default:
            (2, 1) for H20)
        """
        cross_sections = [self.cross_section(energy, f) for f in self.interpolaters]
        return 1 / (
            sum([nMolecules[i] * n * cs for i, cs in enumerate(cross_sections)])
        )

    def get_ratio(
        self,
        energy: float,
        nMolecules: Sequence = (2, 1),
    ):
        """
        Get the ratio of the total cross section of the first atom to the total
        cross section of the molecule.

        Args:
            - energy (float): energy of the neutron in eV.
            - n (float): number of atoms per volume in m^-3 (Default is for H20).
            - nMolecules (Sequence): number of atoms in the molecule (Default:
            (2, 1) for H20)
        """
        cross_sections = [self.cross_section(energy, f) for f in self.interpolaters]
        return (
            cross_sections[0]
            * nMolecules[0]
            / sum(
                [nMolecules[i] * cross_sections[i] for i in range(len(cross_sections))]
            )
        )


class AbsorptionProcessor(CrossSectionProcessor):
    """
    Class that holds the absorption cross section data for neutrons and provides
    a method for calculating the absorption rate for a given energy.
    """

    def __init__(
        self,
        scattering_data: Sequence[pd.DataFrame],
        absorption_data: Sequence[pd.DataFrame],
    ):
        """
        Args:
            - scattering_data (Sequence[pd.DataFrame]): scattering cross section data
            of the form of a sequence of pandas DataFrames with columns x and y
            - absorption_data (Sequence[pd.DataFrame]): absorption cross section data
            of the form of a sequence of pandas DataFrames with columns x and y.

        Scattering and absorption data should be in the same order of atoms and
        of the same length.
        """
        # List with the scattering and absorption data interleaved.
        data = list(itertools.chain(*zip(scattering_data, absorption_data)))
        super().__init__(data)

    def get_total_absorption_rate(
        self,
        energy,
        nMolecules: Sequence = (2, 1),
    ) -> float:
        """
        Get the absorption rate for a given energy.

        Args:
            - energy (float): energy of the neutron in eV.
            - n (float): number of atoms per volume in m^-3 (Default is for H20).
            - nMolecules (Sequence): number of atoms in the molecule (Default:
            (2, 1) for H20)

        Returns: absorption rate
        """
        # list with the scattering and absorption data interleaved.
        cross_sections = [self.cross_section(energy, f) for f in self.interpolaters]
        return 1 / sum(
            [
                nMolecules[i // 2] * (cross_sections[i] / cross_sections[i + 1])
                for i in range(0, len(cross_sections), 2)
            ]
        )

    def get_absorption_rates(self, energy: float) -> list:
        """
        Get the absorption rates for a given energy.

        Args:
            - energy (float): energy of the neutron in eV.

        Returns: list of absorption rates for each atom.
        """
        cross_sections = [self.cross_section(energy, f) for f in self.interpolaters]
        return [
            cross_sections[i + 1] / cross_sections[i]
            for i in range(0, len(cross_sections), 2)
        ]


class SpectrumProcessor(DataProcessor):
    """
    Class that holds the spectrum data for neutrons and provides a method for
    sampling the spectrum.

    Attributes:
        - interpolaters (Sequence[Akima1DInterpolator]): interpolaters for the data.
    """

    def __init__(self, data: pd.DataFrame):
        """
        Args:
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

        Args:
            - num_samples (int): number of samples to take.

        Returns: list of sampled energies.
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
