import pytest
import yaml

from pysisyphus.calculators import OpenMolcas
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("openmolcas")
def test_openmolcas_s1_opt(this_dir):
    geom = geom_loader(this_dir / "trans_butadien.xyz", coord_type="redund")
    inporb_fn = this_dir / "butadien_vdzp.RasOrb"

    kwargs = {
        "basis": "ano-rcc-vdzp",
        "inporb": inporb_fn,
        "charge": 0,
        "mult": 1,
        "rasscf": {
            "ciroot": "5 5 1",
            "mdrlxroot": 2,
        },
    }
    calc = OpenMolcas(**kwargs)
    geom.set_calculator(calc)
    opt_kwargs = {
        "dump": True,
        "thresh": "gau",
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 6
    assert geom.energy == pytest.approx(-154.8758862347)


@using("openmolcas")
def test_cmspdft(this_dir):
    with open(this_dir / "cms_pdft_opt.yaml") as handle:
        run_dict = yaml.safe_load(handle)
    # Overwrite hardcoded path in yaml with path to test directory
    run_dict["calc"]["inporb"] = this_dir / "cms_pdft.RasOrb"
    results = run_from_dict(run_dict)
    assert results.opt_geom.energy == pytest.approx(-437.1389)
