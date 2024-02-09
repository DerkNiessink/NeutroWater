import pandas as pd
import pytest

from neutrons.data_processor import CrossSectionProcessor


test_data = {
    "O_data": pd.read_csv("data/o_cross_t.txt", sep="\s+"),
    "H_data": pd.read_csv("data/h_cross_t.txt", sep="\s+"),
}


@pytest.mark.parametrize(
    "O_data, H_data",
    [
        (test_data["O_data"], test_data["H_data"]),
    ],
)
class TestDataProcessor:

    # Test data from ../data/o_cross_t.txt and converted to m^2
    @pytest.mark.parametrize(
        "energy, expected",
        [
            (1.360720e07, 1.612130e00 * 10**-28),
            (2.000000e08, 4.002810e-01 * 10**-28),
            (6.121900e06, 9.904218e-01 * 10**-28),
        ],
    )
    def test_interpolate(self, energy, expected, O_data, H_data):
        data_processor = CrossSectionProcessor([H_data, O_data])
        f = data_processor.interpolate(O_data, log=True)
        cross_section = data_processor.cross_section(energy, f)
        assert cross_section == pytest.approx(expected, rel=1e-6)

    # Test data from ../data/mean_free_path.csv
    @pytest.mark.parametrize(
        "energy, expected",
        [
            (6.444825e-05, 0.00031210004621352416),
            (0.0003413832, 0.0007147400949399097),
            (0.0253, 0.0045985095787233355),
        ],
    )
    def test_get_mfp(self, energy, expected, O_data, H_data):
        data_processor = CrossSectionProcessor([H_data, O_data])
        mfp = data_processor.get_mfp(energy)
        assert mfp == pytest.approx(expected, rel=1e-3)
