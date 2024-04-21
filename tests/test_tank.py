from neutrowater.models.tank import Tank
import numpy as np
import pytest


class TestTank:
    @pytest.mark.parametrize(
        "radius, height, position, xi",
        [
            (1, 1, np.array([0, 0, 0]), 0.920),
            (2, 2, np.array([0, 5, 0]), 0.920),
            (1, 1, np.array([1, 1, 1]), 0.920),
            (1, 2, np.array([0, 0, 0]), 0.5),
        ],
    )
    def test_init(self, radius, height, position, xi):
        tank = Tank(radius, height, position, xi)
        assert tank.radius == radius
        assert tank.height == height
        assert np.all(tank.position == position)
        assert tank.xi == xi
        assert tank.volume == np.pi * radius**2 * height
        assert tank.energy_loss_frac == 1 / np.exp(xi)

    @pytest.mark.parametrize(
        "test_position, expected",
        [
            (np.array([0, 0, 0]), True),
            (np.array([0.5, 0.5, 0.5]), True),
            (np.array([1, 1, 1]), False),
            (np.array([0, 0, 2]), False),
            (np.array([0, 0, -2]), False),
        ],
    )
    def test_inside(self, test_position, expected):
        tank = Tank(1, 2, np.array([0, 0, 0]), 0.920)
        assert tank.inside(test_position) == expected
