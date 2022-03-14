import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("xtb")
def test_opt_coord_type(this_dir):
    preopt_ct, tsopt_ct, endopt_ct = ("dlc", "redund", "tric")
    run_dict = {
        "geom": {
            "type": "cart",
            "fn": str(this_dir / "test_geoms.trj"),
        },
        "calc": {
            "type": "xtb",
            "quiet": True,
            "pal": 2,
        },
        "preopt": {
            "geom": {
                "type": preopt_ct,
                "freeze_atoms": [0],
            },
            "max_cycles": 1,
        },
        "interpol": {
            "type": "redund",
            "between": 1,
        },
        "cos": {
            "type": "neb",
        },
        "opt": {
            "type": "sd",
            "max_cycles": 1,
        },
        "tsopt": {
            "geom": {
                "type": tsopt_ct,
                "freeze_atoms": [0],
            },
            "max_cycles": 1,
        },
        "irc": {
            "max_cycles": 1,
        },
        "endopt": {
            "max_cycles": 1,
            "geom": {
                "type": endopt_ct,
                "freeze_atoms": [0],
            },
        },
    }
    results = run_from_dict(run_dict)

    assert results.ts_geom.coord_type == tsopt_ct
    assert all([geom.coord_type == endopt_ct for geom in results.end_geoms])
