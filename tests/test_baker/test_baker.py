import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using_pyscf


def calc_getter(charge, mult):
    return PySCF(basis="sto3g", pal=4, charge=charge, mult=mult)


BakerBm = Benchmark("baker", coord_type="redund", calc_getter=calc_getter)


@using_pyscf
@pytest.mark.parametrize("fn, geom, ref_energy", BakerBm)
def test_baker_gs_opt(fn, geom, ref_energy, results_bag):
    opt_kwargs = {
        "thresh": "baker",
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    results_bag.cycles = opt.cur_cycle + 1
    results_bag.is_converged = opt.is_converged
    results_bag.energy = geom.energy
    results_bag.ref_energy = ref_energy

    assert geom.energy == pytest.approx(ref_energy)


def test_baker_synthesis(fixture_store):
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
