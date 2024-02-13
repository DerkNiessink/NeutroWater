import pytest
import numpy as np
import pandas as pd

from neutrons.diffusing_neutrons import DiffusingNeutrons
from neutrons.post_measure import Measurer

test_data = [
    {
        "class_object": DiffusingNeutrons(
            [
                pd.read_csv("data/h_cross_t.txt", sep="\s+"),
                pd.read_csv("data/o_cross_t.txt", sep="\s+"),
            ],
            pd.read_csv("data/neutron_spectrum_normalized.txt", sep=","),
            3,
            (2, 1),
            50,
            50,
            np.array([0, 0, 0]),
            0.920,
        ),
        "energy_expected_dict": {6.0: 3, 2.0: 3, 0.1: 5.0},
        "expected_collisions": 16,
    },
    {
        "class_object": DiffusingNeutrons(
            [
                pd.read_csv("data/h_cross_t.txt", sep="\s+"),
                pd.read_csv("data/o_cross_t.txt", sep="\s+"),
            ],
            pd.read_csv("data/neutron_spectrum_normalized.txt", sep=","),
            3,
            (2, 1),
            50,
            50,
            np.array([0, 0, 0]),
            0.725,
        ),
        "energy_expected_dict": {2.0: 2.0, 1.6: 2.0, 8.0: 3.0},
        "expected_collisions": 20,
    },
]


class TestDiffusingNeutrons:

    # Only use the DiffusingNeutrons objects
    @pytest.mark.parametrize(
        "diffusing_neutrons",
        [test_data[0]["class_object"], test_data[1]["class_object"]],
    )
    def test_random_directions(self, diffusing_neutrons: DiffusingNeutrons) -> None:
        """
        Test the random directions of the neutrons, which should have a norm of 1.
        """
        N = 10
        directions = diffusing_neutrons._random_directions(N)
        assert len(directions) == N
        for direction in directions:
            assert np.linalg.norm(direction) == pytest.approx(1, abs=1e-6)

    # Use only the first DiffusingNeutrons object.
    @pytest.mark.parametrize(
        "diffusing_neutrons, expected_collisions",
        [
            (test_data[0]["class_object"], test_data[0]["expected_collisions"]),
            (test_data[1]["class_object"], test_data[1]["expected_collisions"]),
        ],
    )
    def test_diffuse(
        self, diffusing_neutrons: DiffusingNeutrons, expected_collisions: int
    ) -> None:
        """
        Test the diffuse method and the final energy of the neutrons, which
        should be 1 eV for both tests (H20 and D as moderators).
        """
        # Set the initial energy of the neutrons to 2 MeV
        for neutron in diffusing_neutrons.neutrons:
            neutron.energies[0] = 2 * 10**6

        diffusing_neutrons.diffuse(nCollisions=expected_collisions - 5)
        diffusing_neutrons.diffuse(nCollisions=5)
        measure = Measurer(diffusing_neutrons)
        assert measure.energies()[0][-1] == pytest.approx(1, abs=0.2)

    def test_get_number_escaped(self) -> None:
        """
        Test the get_number_escaped method, in the case all neutrons should escape
        from the tank.
        """
        neutrons = DiffusingNeutrons(
            [
                pd.read_csv("data/h_cross_t.txt", sep="\s+"),
                pd.read_csv("data/o_cross_t.txt", sep="\s+"),
            ],
            pd.read_csv("data/neutron_spectrum_normalized.txt", sep=","),
            3,
            (2, 1),
            0.1,
            0.1,
            np.array([0, 0, 0]),
            0.920,
        )
        neutrons.neutrons[0].positions[0] = np.array([0.1, 0, 0])
        neutrons.neutrons[1].positions[0] = np.array([1.01, 0, 0])
        neutrons.neutrons[2].positions[0] = np.array([0.5, 0, 0])

        neutrons.diffuse(nCollisions=30)
        measure = Measurer(neutrons)
        assert measure.number_escaped() == 3
