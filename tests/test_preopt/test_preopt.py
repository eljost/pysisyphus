import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("pyscf")
@pytest.mark.parametrize("coord_type, cos_type", [("dlc", "gs"), ("cart", "neb")])
def test_preopt(coord_type, cos_type, this_dir):
    run_dict = {
        "geom": {
            "type": coord_type,
            "fn": "lib:xtb_rx/00_c2no2.trj",
        },
        "calc": {
            "type": "pyscf",
            "basis": "sto3g",
            "verbose": 0,
            "pal": 2,
        },
        "preopt": {
            "max_cycles": 3,
            "geom": {
                "type": "tric",
            },
        },
        "cos": {
            "type": cos_type,
        },
        "opt": {
            "type": "sd",
            "max_cycles": 1,
            "align": False,
        },
    }
    if cos_type == "neb":
        run_dict["interpol"] = {
            "type": "redund",
            "between": 1,
        }
    results = run_from_dict(run_dict)

    # first_en = -258.2759597602544
    # assert results.preopt_first_geom.energy == pytest.approx(first_en)
    # assert results.cos.images[0].energy == pytest.approx(first_en)

    # last_en = -258.42817036587434
    # assert results.preopt_last_geom.energy == pytest.approx(last_en)
    # assert results.cos.images[-1].energy == pytest.approx(last_en)

    first_ref = geom_loader(this_dir / "ref_first_preopt.xyz")
    assert first_ref.rmsd(results.cos.images[0]) <= 1e-4
    last_ref = geom_loader(this_dir / "ref_last_preopt.xyz")
    assert last_ref.rmsd(results.cos.images[-1]) <= 1e-4
