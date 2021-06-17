import pytest

from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@pytest.mark.skip
@using("xtb")
def test_oniom_opt_big():
    run_dict = {
        "geom": {
            "type": "redund",
            # frag.py raffinose.xyz 21 14 33 10
            "fn": "lib:birkholz/raffinose.xyz",
        },
        "calc": {
            "type": "oniom",
            "calcs": {
                "real": {
                    "type": "xtb",
                    "gfn": 0,
                    "pal": 6,
                    "quiet": True,
                },
                "high": {
                    "type": "xtb",
                    "gfn": 2,
                    "pal": 6,
                    "quiet": True,
                },
            },
            "models": {
                "high": {
                    "inds": "11..21,33..36,50..57",
                    "calc": "high",
                },
            },
        },
        "opt": {
            "thresh": "gau_loose",
        },
    }
    results = run_from_dict(run_dict)

    assert results.opt_geom.energy == pytest.approx(-110.75489874)


@using("xtb")
@pytest.mark.parametrize(
    "opt_dict",
    [
        {
            "type": "rfo",
        },
        {
            "type": "lbfgs",
            "line_search": True,
        },
        {
            "type": "oniom",
        },
    ],
    ids=["rfo_ref", "lbfgs_ref", "oniom"],
)
def test_oniom_opt_small(opt_dict):
    opt_dict.update({
        "thresh": "gau",
        "step": "full",
        # "step": "high",
        "dump": True,
    })
    run_dict = {
        "geom": {"type": "cart", "fn": "lib:acetaldehyd_oniom.xyz"},
        "calc": {
            "type": "oniom",
            "calcs": {
                "real": {
                    "type": "xtb",
                    "gfn": 0,
                    "pal": 6,
                    "quiet": True,
                },
                "high": {
                    "type": "xtb",
                    "gfn": 2,
                    "pal": 6,
                    "quiet": True,
                },
            },
            "models": {
                "high": {
                    "inds": [4, 5, 6],
                    "calc": "high",
                },
            },
        },
        "opt": opt_dict,
    }
    results = run_from_dict(run_dict)

    assert results.opt_geom.energy == pytest.approx(-10.419331913)


@pytest.mark.skip
def test_dmp():
    run_dict = {
        # "geom": {"type": "cart", "fn": "inp.xyz"},
        "geom": {"type": "redund", "fn": "inp.xyz"},
        "calc": {
            "type": "oniom",
            "calcs": {
                "real": {
                    "type": "pyscf",
                    "basis": "sto3g",
                    "pal": 2,
                    "verbose": 0,
                },
                "high": {
                    "type": "pyscf",
                    "basis": "321g",
                    "pal": 2,
                    "verbose": 0,
                },
            },
            "models": {
                "high": {
                    "inds": [15, 14, 1, 0],
                    "calc": "high",
                },
            },
        },
        # "opt": {
            # "type": "oniom",
            # # "thresh": "gau_loose",
            # "thresh": "gau",
        # },
        "opt": {"thresh": "gau"},
    }
    results = run_from_dict(run_dict)

    # assert results.opt_geom.energy == pytest.approx(-10.419331913)
