import numpy as np
import pytest

from pysisyphus.cos import multi_state


@pytest.fixture
def data(this_dir):
    """Data from tables S2 to S5 from SI of 10.1063/5.0021923"""
    data = np.load(this_dir / "multi_state_inp.npz")
    return data


@pytest.mark.parametrize(
    "key, ref_inds",
    (
        ("S1", [9]),
        ("S2", [4]),
        ("S3", [4]),
        ("S4", [5]),
        ("S5", [3, 9]),
    ),
)
def test_select(key, ref_inds, data):
    key_data = data[f"{key}"]
    all_energies = key_data[:, [1, 2]]
    inds = multi_state.determine_meisc_images(all_energies)
    print(key, all_energies.shape, inds)
    np.testing.assert_allclose(inds, ref_inds)
