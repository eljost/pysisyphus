from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("pyscf")
def test_diels_alder_growing_string():
    run_dict = {
        "preopt": {
            "max_cycles": 5,
        },
        "cos": {
            "type": "gs",
            "fix_ends": True,
            "max_nodes": 8,
        },
        "opt": {
            "type": "string",
            "stop_in_when_full": 0,
        },
        "tsopt": {
            "type": "rsirfo",
            "do_hess": True,
            "hessian_recalc": 5,
        },
        "irc": {
            "type": "eulerpc",
            "opt_ends": "fragments",
            "corr_func": "scipy",
            "rms_grad_thresh": 2.5e-3,
        },
        "calc": {
            "type": "pyscf",
            "pal": 2,
            "basis": "321g",
        },
        "xyz": "lib:diels_alder_interpolated.trj",
        "coord_type": "dlc",
    }
    result = run_from_dict(run_dict)
