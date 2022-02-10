import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
# from pysisyphus.helpers_pure import filter_fixture_store
from pysisyphus.intcoords.augment_bonds import augment_bonds
from pysisyphus.testing import using
from pysisyphus.tsoptimizers import *


def calc_getter(charge, mult):
    return PySCF(basis="321g", pal=1, verbose=0, charge=charge, mult=mult)


BakerTSBm = Benchmark(
    "baker_ts",
    coord_type="redund",
    calc_getter=calc_getter,
    exclude=(1, 9, 10, 14, 21, 23),
)


@pytest.mark.benchmark
@using("pyscf")
@pytest.mark.parametrize("fn, geom, ref_energy", BakerTSBm)
def test_baker_tsopt(fn, geom, ref_energy, results_bag):
    geom = augment_bonds(geom)

    opt_kwargs = {
        "thresh": "baker",
        "max_cycles": 50,
        "trust_radius": 0.1,
        "trust_max": 0.3,
        "min_line_search": True,
        "max_line_search": True,
    }
    opt = RSPRFOptimizer(geom, **opt_kwargs)
    # opt = RSIRFOptimizer(geom, **opt_kwargs)
    # opt = TRIM(geom, **opt_kwargs)
    opt.run()

    results_bag.cycles = opt.cur_cycle + 1
    results_bag.is_converged = opt.is_converged
    results_bag.energy = geom.energy
    results_bag.ref_energy = ref_energy

    assert geom.energy == pytest.approx(ref_energy)


@pytest.mark.benchmark
# @filter_fixture_store("test_baker_tsopt")
def test_baker_tsopt_synthesis(fixture_store):
    for i, fix in enumerate(fixture_store):
        print(i, fix)

    tot_cycles = 0
    converged = 0
    bags = fixture_store["results_bag"]
    for k, v in bags.items():
        if not k.startswith("test_baker_tsopt"):
            continue
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


@using("pyscf")
@pytest.mark.parametrize(
    "proj, ref_cycle",
    [
        pytest.param(True, 14, marks=pytest.mark.skip_ci),
        pytest.param(False, 12, marks=pytest.mark.skip_ci),
    ],
)
def test_diels_alder_ts(ref_cycle, proj):
    """
    https://onlinelibrary.wiley.com/doi/epdf/10.1002/jcc.21494
    """

    geom = geom_loader(
        "lib:baker_ts/09_parentdieslalder.xyz",
        coord_type="redund",
    )

    calc_kwargs = {
        "charge": 0,
        "mult": 1,
        "pal": 4,
    }
    geom.set_calculator(PySCF(basis="321g", **calc_kwargs))
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
