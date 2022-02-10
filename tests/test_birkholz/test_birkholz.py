import pytest

from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader
# from pysisyphus.helpers_pure import filter_fixture_store
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


@pytest.mark.benchmark
@using("xtb")
@pytest.mark.parametrize(
    "xyz_fn, charge, mult", [
    ("artemisin.xyz", 0, 1),
    ("aspartame.xyz", 0, 1),
    ("avobenzone.xyz", 0, 1),
    ("azadirachtin.xyz", 0, 1),
    ("bisphenol_a.xyz", 0, 1),
    ("cetirizine.xyz", 0, 1),
    ("codeine.xyz", 0, 1),
    ("diisobutyl_phthalate.xyz", 0, 1),
    ("easc.xyz", 0, 1),
    ("estradiol.xyz", 0, 1),
    ("inosine.xyz", 1, 1),
    ("maltose.xyz", 0, 1),
    ("mg_porphin.xyz", 0, 1),
    ("ochratoxin_a.xyz", 0, 1),
    ("penicillin_v.xyz", 0, 1),
    ("raffinose.xyz", 0, 1),
    ("sphingomyelin.xyz", 0, 1),
    ("tamoxifen.xyz", 0, 1),
    ("vitamin_c.xyz", 0, 1),
    ("zn_edta.xyz", -2, 1),
    ]
)
def test_birkholz_set(xyz_fn, charge, mult, results_bag):
    print(f"@{xyz_fn}")
    coord_type = "redund"
    # coord_type = "cart"
    geom = geom_loader(f"lib:birkholz/{xyz_fn}", coord_type=coord_type)
    calc = XTB(charge=charge, mult=mult, pal=4)
    geom.set_calculator(calc)

    opt_kwargs_base = {
        "max_cycles": 150,
        "thresh": "baker",
        "overachieve_factor": 2,
        "dump": True,
        "adapt_step_func": True,
    }

    opt_kwargs = opt_kwargs_base.copy()
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    results_bag.cycles = opt.cur_cycle + 1
    results_bag.is_converged = opt.is_converged

    assert opt.is_converged


@pytest.mark.benchmark
# @filter_fixture_store("test_birkholz_set")
def test_birkholz_set_synthesis(fixture_store):
    for i, fix in enumerate(fixture_store):
        print(i, fix)

    tot_cycles = 0
    converged = 0
    bags = fixture_store["results_bag"]
    for k, v in bags.items():
        if not k.startswith("test_birkholz_set"):
            continue
        print(k)
        try:
            tot_cycles += v["cycles"]
            converged += 1 if v["is_converged"] else 0
            for kk, vv in v.items():
                print("\t", kk, vv)
        except KeyError:
            print("\tFailed!")
    print(f"Total cycles: {tot_cycles}")
    print(f"Converged: {converged}/{len(bags)}")
