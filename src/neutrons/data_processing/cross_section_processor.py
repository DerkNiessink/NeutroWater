import numpy as np
from scipy.interpolate import Akima1DInterpolator
from typing import Sequence

from neutrons.data_processing.data_processor import DataProcessor


class CrossSectionProcessor(DataProcessor):

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
