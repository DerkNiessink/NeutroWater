import pytest
import numpy as np

from neutrons.diffusing_neutrons import DiffusingNeutrons


@pytest.fixture
def get_diffusing_neutron():
    return DiffusingNeutrons(
        [1, 2, 3],
        [1, 2, 3],
        [np.array([0, 0, 0])],
        [1],
        1,
        1,
        np.array([0, 0, 0]),
        0.920,
    )


@pytest.fixture
def get_diffusing_neutron_2():
    return DiffusingNeutrons(
        [5, 3, 3],
        [1, 2, 8],
        [np.array([1, 2, 2])],
        [20],
        5,
        3,
        np.array([0, 5, 0]),
        0.2,
    )


class TestDiffusingNeutrons:
    @pytest.mark.parametrize(
        "mean_free_paths, energies, initial_positions, initial_energies, radius_tank, height_tank, position_tank, xi",
        [
            (
                [1, 2, 3],
                [1, 2, 3],
                [np.array([0, 0, 0])],
                [1],
                1,
                1,
                np.array([0, 0, 0]),
                0.920,
            ),
            (
                [1, 2, 3],
                [1, 2, 3],
                [np.array([0, 0, 0])],
                [1],
                1,
                1,
                np.array([0, 0, 0]),
                0.920,
            ),
            (
                [1, 2, 3],
                [1, 2, 3],
                [np.array([0, 0, 0])],
                [1],
                1,
                1,
                np.array([0, 0, 0]),
                0.920,
            ),
        ],
    )
    def test_init(
        self,
        mean_free_paths,
        energies,
        initial_positions,
        initial_energies,
        radius_tank,
        height_tank,
        position_tank,
        xi,
    ):
        diffusing_neutrons = DiffusingNeutrons(
            mean_free_paths,
            energies,
            initial_positions,
            initial_energies,
            radius_tank,
            height_tank,
            position_tank,
            xi,
        )
        assert diffusing_neutrons.energy_loss_frac == 1 / np.exp(xi)
        assert diffusing_neutrons.mean_free_paths == mean_free_paths
        assert diffusing_neutrons.energies == energies
        assert diffusing_neutrons.nNeutrons == len(initial_positions)
        assert diffusing_neutrons.tank.radius == radius_tank
        assert diffusing_neutrons.tank.height == height_tank
        assert np.all(diffusing_neutrons.tank.position == position_tank)
        assert diffusing_neutrons.tank.xi == xi
        assert diffusing_neutrons.tank.volume == np.pi * radius_tank**2 * height_tank
        assert diffusing_neutrons.tank.energy_loss_frac == 1 / np.exp(xi)
        assert diffusing_neutrons.neutrons.energies == initial_energies
        assert np.all(diffusing_neutrons.neutrons.positions == initial_positions)

    def test_random_directions(self):
        diffusing_neutrons = DiffusingNeutrons(
            [1, 2, 3],
            [1, 2, 3],
            [np.array([0, 0, 0])],
            [1],
            1,
            1,
            np.array([0, 0, 0]),
            0.920,
        )
        N = 10
        directions = diffusing_neutrons._random_directions(N)
        assert len(directions) == N
        for direction in directions:
            assert np.linalg.norm(direction) == pytest.approx(1, abs=1e-6)

    @pytest.mark.parametrize(
        "energy, expected",
        [
            (5, 2),
            (4, 2),
            (2, 3),
        ],
    )
    def test_mf_lookup(self, energy, expected):
        diffusing_neutrons = DiffusingNeutrons(
            [1, 2, 3],
            [3, 4, 2],
            [np.array([0, 0, 0])],
            [1],
            1,
            1,
            np.array([0, 0, 0]),
            0.920,
        )
        assert diffusing_neutrons._mf_lookup(energy) == expected
