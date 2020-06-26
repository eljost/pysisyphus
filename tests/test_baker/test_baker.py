import os

import numpy as np
import pytest

from pysisyphus.helpers import get_baker_geoms
from pysisyphus.color import red, green
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.testing import using_pyscf



@using_pyscf
@pytest.mark.parametrize(
    "name, geom, ref_energy",
    [(name, geom, ref_energy) for name, (geom, ref_energy)
     in get_baker_geoms(coord_type="redund").items()]
)
def test_baker_gs_opt(name, geom, ref_energy, results_bag):
    opt_kwargs = {
        "thresh": "baker",
    }
    print(f"@Running {name}")
    pal = min(os.cpu_count(), 4)
    geom.set_calculator(PySCF(basis="sto3g", pal=pal))
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
        tot_cycles += v["cycles"]
        converged += 1 if v["is_converged"] else 0
        for kk, vv in v.items():
            print("\t", kk, vv)
    print(f"Total cycles: {tot_cycles}")
    print(f"Converged: {converged}/{len(bags)}")
