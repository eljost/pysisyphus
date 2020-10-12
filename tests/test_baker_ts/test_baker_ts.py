import itertools as it

import pytest

from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import do_final_hessian, get_baker_ts_geoms_flat, geom_loader
from pysisyphus.intcoords.augment_bonds import augment_bonds
from pysisyphus.testing import using_pyscf, using_gaussian16
from pysisyphus.tsoptimizers import *


def get_geoms():
    fails = (
        "15_hocl.xyz",  # SCF goes completely nuts
    )
    works = (
        "01_hcn.xyz",
        "02_hcch.xyz",
        "03_h2co.xyz",
        "04_ch3o.xyz",
        "05_cyclopropyl.xyz",
        "06_bicyclobutane.xyz",
        "07_bicyclobutane.xyz",
        "08_formyloxyethyl.xyz",
        "09_parentdieslalder.xyz",
        "12_ethane_h2_abstraction.xyz",
        "13_hf_abstraction.xyz",
        "14_vinyl_alcohol.xyz",
        "16_h2po4_anion.xyz",
        "17_claisen.xyz",
        "18_silyene_insertion.xyz",
        "19_hnccs.xyz",
        # "20_hconh3_cation.xyz",
        "21_acrolein_rot.xyz",
        # "22_hconhoh.xyz",
        "23_hcn_h2.xyz",
        "24_h2cnh.xyz",
        "25_hcnh2.xyz",
    )
    alpha_negative = ()
    no_imag = (
        "10_tetrazine.xyz",
        "11_trans_butadiene.xyz",
    )
    only = (
        # "18_silyene_insertion.xyz",
    )
    use = (
        # fails,
        works,
        # alpha_negative,
        # no_imag,
        # only,
        )
    use_names = list(it.chain(*use))
    geom_data = get_baker_ts_geoms_flat(
        coord_type="redund",
        coord_kwargs={"rebuild": True,},
    )
    return [_ for _ in geom_data if _[0] in use_names]


@using_pyscf
@pytest.mark.parametrize(
    "name, geom, charge, mult, ref_energy",
    get_geoms()
)
def test_baker_tsopt(name, geom, charge, mult, ref_energy, results_bag):
    calc_kwargs = {
        "charge": charge,
        "mult": mult,
        "pal": 1,
    }

    print(f"@Running {name}")
    # geom.set_calculator(Gaussian16(route="HF/3-21G", **calc_kwargs))
    geom.set_calculator(PySCF(basis="321g", **calc_kwargs))
    if geom.coord_type != "cart":
        geom = augment_bonds(geom)

    opt_kwargs = {
        "thresh": "baker",
        "max_cycles": 50,
        "trust_radius": 0.3,
        "trust_max": 0.3,
        # "min_line_search": True,
        # "max_line_search": True,
    }
    # opt = RSPRFOptimizer(geom, **opt_kwargs)
    opt = RSIRFOptimizer(geom, **opt_kwargs)
    # opt = TRIM(geom, **opt_kwargs)
    opt.run()

    # Without symmetry restriction this lower lying TS will be obtained.
    if name.startswith("22_"):
        ref_energy = -242.25695787

    results_bag.cycles = opt.cur_cycle + 1
    results_bag.is_converged = opt.is_converged
    results_bag.energy = geom.energy
    results_bag.ref_energy = ref_energy

    print(f"\t@Converged: {opt.is_converged}, {opt.cur_cycle+1} cycles")

    assert geom.energy == pytest.approx(ref_energy)
    print("\t@Energies match!")

    return opt.cur_cycle+1


def test_baker_tsopt_synthesis(fixture_store):
    for i, fix in enumerate(fixture_store):
        print(i, fix)

    tot_cycles = 0
    converged = 0
    bags = fixture_store["results_bag"]
    for k, v in bags.items():
        print(k)
        try:
            tot_cycles += v["cycles"]
            energy_matches = v["energy"] == pytest.approx(v["ref_energy"])
            converged += 1 if v["is_converged"] and energy_matches else 0
            for kk, vv in v.items():
                print("\t", kk, vv)
        except KeyError:
            print("\tFailed!")
    print(f"Total cycles: {tot_cycles}")
    print(f"Converged: {converged}/{len(bags)}")


@using_pyscf
@pytest.mark.parametrize(
    "define_prims, proj, ref_cycle", [
        (None, True, 14),
        pytest.param(None, False, 12),
        pytest.param([[1, 5], [0, 4], [4, 10], [5, 11], [13, 1], [12, 0]], False, 12),
        pytest.param([[1, 5], [0, 4], [13, 1], [12, 0]], False, 10),
    ]
)
def test_diels_alder_ts(define_prims, ref_cycle, proj):
    """
        https://onlinelibrary.wiley.com/doi/epdf/10.1002/jcc.21494
    """

    coord_kwargs = None
    augment = True

    if define_prims:
        coord_kwargs = {
            "define_prims": define_prims,
        }
        augment = False

    geom = geom_loader("lib:baker_ts/09_parentdieslalder.xyz",
                       coord_type="redund",
                       coord_kwargs=coord_kwargs,
    )

    calc_kwargs = {
        "charge": 0,
        "mult": 1,
        "pal": 4,
    }
    # geom.set_calculator(Gaussian16(route="HF/3-21G", **calc_kwargs))
    geom.set_calculator(PySCF(basis="321g", **calc_kwargs))
    if augment:
        geom = augment_bonds(geom, proj=proj)

    opt_kwargs = {
        "thresh": "baker",
        "max_cycles": 50,
        "trust_radius": 0.3,
        "trust_max": 0.3,
        "hessian_recalc": 5,
        "dump": True,
        "overachieve_factor": 2,
    }
    opt = RSIRFOptimizer(geom, **opt_kwargs)
    opt.run()

    print(f"\t@Converged: {opt.is_converged}, {opt.cur_cycle+1} cycles")

    ref_energy = -231.60320857
    assert geom.energy == pytest.approx(ref_energy)
    print("\t@Energies match!")
    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
