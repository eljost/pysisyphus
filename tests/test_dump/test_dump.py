import h5py

from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("pyscf")
def test_hcn_neb():
    run_dict = {
        "preopt": {
            "max_cycles": 3,
        },
        "interpol": {
            "type": "idpp",
            "between": 3,
        },
        "cos": {
            "type": "neb",
        },
        "opt": {
            "type": "qm",
            "align": True,
        },
        "tsopt": {
            "type": "rsirfo",
        },
        "calc": {
            "type": "pyscf",
            "pal": 1,
            "basis": "321g",
        },
        "geom": {
            "type": "cart",
            "fn": ["lib:hcn.xyz", "lib:hcn_iso_ts.xyz", "lib:nhc.xyz"],
        },
    }
    results = run_from_dict(run_dict)

    assert results.cos_opt.is_converged
    assert results.cos_opt.cur_cycle == 18
    assert results.ts_opt.is_converged
    assert results.ts_opt.cur_cycle == 1

    with h5py.File("optimization.h5", "r") as handle:
        groups = list(handle.keys())
    ref_groups = ("first_pre", "last_pre", "opt", "tsopt")
    assert set(groups) == set(ref_groups)
