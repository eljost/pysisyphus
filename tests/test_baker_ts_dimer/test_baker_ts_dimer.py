from pathlib import Path
import pickle
import shutil

import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators import Dimer
from pysisyphus.calculators.PySCF import PySCF
# from pysisyphus.helpers_pure import filter_fixture_store
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.testing import using


BakerTSBm = Benchmark(
    "baker_ts",
    coord_type="cart",
    exclude=(10, ),
    # inv_exclude=True,
)


@pytest.mark.benchmark
@using("pyscf")
@pytest.mark.parametrize("fn, geom, charge, mult, ref_energy", BakerTSBm.geom_iter)
def test_baker_ts_dimer(fn, geom, charge, mult, ref_energy, results_bag, this_dir):
    init_orients = "Ns"
    # init_orients = "Ns_start"
    with open(this_dir / init_orients, "rb") as handle:
        N_INITS = pickle.load(handle)

    id_ = fn[:2]
    calc_kwargs = {
        "charge": charge,
        "mult": mult,
        "base_name": Path(fn).stem,
    }
    calc = PySCF("321g", pal=2, verbose=0, **calc_kwargs)

    dimer_kwargs = {
        "rotation_method": "fourier",
        "calculator": calc,
        "N_raw": N_INITS[id_],
        "length": 0.0189,
        "rotation_tol": 5,
        "rotation_disable_pos_curv": True,
        "trans_force_f_perp": True,
        # "rotation_thresh": 0.002,  # Shang 2010
        # "rotation_thresh": 0.001,
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
        # "max_force_only": True,  # Shang 2010
        # "rms_force": 0.0013,  # Shang 2010
    }
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    shutil.copy("final_geometry.xyz", f"{id_}_final_geometry.xyz")

    energies_match = geom.energy == pytest.approx(ref_energy)

    results_bag.cycles = opt.cur_cycle + 1
    results_bag.is_converged = opt.is_converged
    results_bag.energies_match = energies_match
    results_bag.dimer_force_evals = dimer.force_evals

    assert opt.is_converged
    assert energies_match

    opt_cycs = opt.cur_cycle + 1
    print(f"@{fn} converged: {opt_cycs} optimization cycles, "
          f"{dimer.force_evals} force evaluations."
    )


@pytest.mark.benchmark
# @filter_fixture_store("test_baker_ts_dimer")
def test_baker_ts_dimer_synthesis(fixture_store):
    converged = 0
    tot_cycles = 0
    tot_dimer_force_evals = 0
    tot_cycles_failed = 0
    tot_dimer_force_evals_failed = 0

    bags = fixture_store["results_bag"]
    for k, v in bags.items():
        if not k.startswith("test_baker_ts_dimer"):
            continue
        print(k)
        try:
            energies_match = v["energies_match"]
            converged += 1 if energies_match else 0
            cycles = v["cycles"]
            force_evals = v["dimer_force_evals"]
            if energies_match:
                tot_cycles += cycles
                tot_dimer_force_evals += force_evals
            else:
                tot_cycles_failed += cycles
                tot_dimer_force_evals_failed += force_evals
            for kk, vv in v.items():
                print("\t", kk, vv)
        except KeyError:
            print("\tFailed!")

    print("### Converged")
    print(f"\tTotal cycles: {tot_cycles}")
    print(f"\tTotal dimer force evaluations: {tot_dimer_force_evals}")
    print(f"\tConverged: {converged}/{len(bags)}")

    print("### Failed")
    print(f"\tTotal cycles: {tot_cycles_failed}")
    print(f"\tTotal dimer force evaluations:  {tot_dimer_force_evals_failed}")
