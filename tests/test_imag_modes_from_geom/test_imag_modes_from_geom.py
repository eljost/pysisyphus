import numpy as np
import pytest

from pysisyphus.helpers import geom_loader, imag_modes_from_geom


def test_imag_modes_from_geom(this_dir):
    # Obtained at the xtb-GFN2 level of theory
    H = np.loadtxt(this_dir / "calculated_final_cart_hessian")
    geom = geom_loader("lib:hcn_iso_ts_opt_xtb.xyz")
    geom._hessian = H

    imag_modes = imag_modes_from_geom(geom)
    assert len(imag_modes) == 1

    imag_mode = imag_modes[0]
    assert imag_mode.nu == pytest.approx(-1426.207289152)

    with open("imag_mode.trj", "w") as handle:
        handle.write(imag_mode.trj_str)
