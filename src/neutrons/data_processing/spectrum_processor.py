import numpy as np
from scipy.interpolate import Akima1DInterpolator
from typing import Sequence

from neutrons.data_processing.data_processor import DataProcessor


class SpectrumProcessor(DataProcessor):

    def __init__(self, data, log):
        x = data.iloc[:, 0]
        self.bounds = (min(x), max(x))
        super().__init__(data, log)

    @classmethod
    def rejection_sampler(self, f: Akima1DInterpolator, num_samples: int) -> list:
        sampled_points = []
        while len(sampled_points) < num_samples:
            x = np.random.uniform(*self.bounds)
            y = np.random.uniform(0, 1)

            px = f(x)

            if y <= px:
                sampled_points.append(x)

        return sampled_points
