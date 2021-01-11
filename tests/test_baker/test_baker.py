import os

import numpy as np
import pytest

from pysisyphus.benchmark_sets import get_baker_fns, get_baker_ref_energies
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using_pyscf


@using_pyscf
@pytest.mark.parametrize(
    # Drop prefix
    "fn", get_baker_fns()[1]
)
def test_baker_gs_opt(fn, results_bag):
    geom = geom_loader(f"lib:baker/{fn}", coord_type="redund")
    ref_energy = get_baker_ref_energies()[fn]
    opt_kwargs = {
        "thresh": "baker",
        "adapt_step_func": False,
    }
    print(f"@Running {fn}")
    pal = min(os.cpu_count(), 4)
    geom.set_calculator(PySCF(basis="sto3g", pal=pal, verbose=0))
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    results_bag.cycles = opt.cur_cycle + 1
    results_bag.is_converged = opt.is_converged
    results_bag.energy = geom.energy
    results_bag.ref_energy = ref_energy

    assert np.allclose(geom.energy, ref_energy)

    return opt.cur_cycle + 1, opt.is_converged


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
