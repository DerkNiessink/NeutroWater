from neutrons.models.neutrons import Neutron, Neutrons

import pytest
import numpy as np


class TestNeutron:
    @pytest.mark.parametrize(
        "initial_energy, initial_position",
        [
            (1, np.array([0, 0, 0])),
            (2, np.array([1, 1, 1])),
            (3, np.array([2, 2, 2])),
        ],
    )
    def test_init(self, initial_energy, initial_position):
        self.neutron = Neutron(initial_energy, initial_position)
        assert self.neutron.energies == [initial_energy]
        assert np.all(self.neutron.positions == [initial_position])

    @pytest.mark.parametrize(
        "distance, direction, expected",
        [
            (1, np.array([1, 0, 0]), np.array([1, 0, 0])),
            (2, np.array([0, 1, 0]), np.array([0, 2, 0])),
            (3, np.array([0, 0, 1]), np.array([0, 0, 3])),
        ],
    )
    def test_travel(self, distance, direction, expected):
        self.neutron = Neutron(1, np.array([0, 0, 0]))
        self.neutron.travel(distance, direction)
        assert np.all(self.neutron.positions[-1] == expected)

    @pytest.mark.parametrize(
        "energy_loss_frac, expected",
        [
            (0.5, 1),
            (0.75, 1.5),
            (0.25, 0.5),
        ],
    )
    def test_collide(self, energy_loss_frac, expected):
        self.neutron = Neutron(2, np.array([0, 0, 0]))
        self.neutron.collide(energy_loss_frac)
        assert self.neutron.energies[-1] == expected


@pytest.mark.parametrize(
    "energies, positions",
    [
        (
            [1, 2, 3],
            [np.array([0, 0, 0]), np.array([1, 1, 1]), np.array([2, 2, 2])],
        ),
        ([1, 2], [np.array([0, 0, 0]), np.array([1, 1, 1])]),
    ],
)
class TestNeutrons:
    def test_init(self, energies, positions):
        self.neutrons = Neutrons(energies, positions)
        assert np.all(self.neutrons.energies == energies)
        assert np.all(self.neutrons.positions == positions)
        assert len(self.neutrons) == len(energies)
        for i, neutron in enumerate(self.neutrons.neutrons):
            assert neutron.energies == [energies[i]]
            assert np.all(neutron.positions == [positions[i]])

    def test_iter(self, energies, positions):
        self.neutrons = Neutrons(energies, positions)
        for i, neutron in enumerate(self.neutrons):
            assert neutron.energies == [energies[i]]
            assert np.all(neutron.positions == [positions[i]])
