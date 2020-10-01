import pytest

from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("pyscf")
def test_run_geom_section_union():
    ref_energy = -160.7433945831947

    def run_assert(run_dict, prim_len):
        res = run_from_dict(run_dict)
        geom = res.calced_geoms[0]
        assert geom.energy == pytest.approx(ref_energy)
        int_ = geom.internal
        assert len(int_.prim_indices_set) == prim_len

        # import numpy as np
        # H = geom.hessian
        # w, v = np.linalg.eigh(H)
        # first = 5
        # print("eigvals", w[:first])

    def get_run_dict():
        return {
            "calc": {
                "type": "pyscf",
                "basis": "sto3g",
            },
            "geom": {
                "type": "redund",
                "fn": "lib:test_union_ts_001.xyz",
            },
        }

    run_assert(get_run_dict(), prim_len=56)

    # Use union of internal coordinates from the two geometries
    # in the .trj file.
    run_dict = get_run_dict()
    run_dict["geom"]["union"] = "lib:test_union_rx_001.trj"
    run_assert(run_dict, prim_len=67)
