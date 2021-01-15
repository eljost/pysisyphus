from pathlib import Path
import pickle
import os
import shutil

from natsort import natsorted
import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators import Dimer
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.testing import using


BakerTSBm = Benchmark(
    "baker_ts",
    coord_type="cart",
    # exclude=(10, ),
    # inv_exclude=True,
)


@using("pyscf")
@pytest.mark.parametrize("fn, geom, charge, mult, ref_energy", BakerTSBm.geom_iter)
def test_baker_ts_dimer(fn, geom, charge, mult, ref_energy, results_bag, this_dir):
    with open(this_dir / "Ns", "rb") as handle:
        N_INITS = pickle.load(handle)

    id_ = fn[:2]
    calc_kwargs = {
        "charge": charge,
        "mult": mult,
        "pal": 2,
        "verbose": 0,
        "base_name": Path(fn).stem,
    }
    calc = PySCF("321g", **calc_kwargs)

    N_raw = N_INITS[id_]
    dimer_kwargs = {
        "rotation_method": "fourier",
        "calculator": calc,
        "N_raw": N_INITS[id_],
        "length": 0.0189,
        "rotation_tol": 5,
        "rotation_disable_pos_curv": True,
        "trans_force_f_perp": True,
    }
    dimer = Dimer(**dimer_kwargs)
    geom.set_calculator(dimer)

    opt_kwargs = {
        "thresh": "baker",
        "precon": True,
        "max_step_element": 0.25,
        "max_cycles": 50,
        "c_stab": 0.103,
        "dump": True,
    }
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    shutil.copy("final_geometry.xyz", f"{id_}_final_geometry.xyz")

    energies_match = geom.energy == pytest.approx(ref_energy)

    results_bag.cycles = opt.cur_cycle + 1
    results_bag.is_converged = opt.is_converged
    results_bag.energies_match = energies_match
    results_bag.force_evals = dimer.force_evals

    assert opt.is_converged
    assert energies_match

    print(f"@{fn} converged using {dimer.force_evals} force evaluations")


def test_baker_ts_dimer_synthesis(fixture_store):
    tot_cycles = 0
    converged = 0
    tot_force_evals = 0
    bags = fixture_store["results_bag"]
    for k, v in bags.items():
        print(k)
        try:
            tot_cycles += v["cycles"]
            converged += 1 if v["energies_match"] else 0
            tot_force_evals += v["force_evals"]
            for kk, vv in v.items():
                print("\t", kk, vv)
        except KeyError:
            print("\tFailed!")
    print(f"Total cycles: {tot_cycles}")
    print(f"Total force evaluations: {tot_force_evals}")
    print(f"Converged: {converged}/{len(bags)}")
