#!/usr/bin/env python3

from pathlib import Path
import os
import sys
import time

import cloudpickle
from natsort import natsorted
import numpy as np
import pytest

from pysisyphus.calculators import Dimer, Gaussian16
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_xyz_file, geom_from_library, \
                               get_baker_ts_geoms_flat, do_final_hessian
from pysisyphus.tsoptimizers.dimer import dimer_method
from pysisyphus.tsoptimizers.dimerv2 import dimer_method as dimer_method_v2
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.testing import using


def make_N_init_dict():
    THIS_DIR = Path(os.path.abspath(os.path.dirname(__file__)))
    xyz_path = THIS_DIR / "../../xyz_files/baker_ts"
    xyzs = natsorted(xyz_path.glob("*.xyz"))
    N_dict = dict()
    for guess, initial in [xyzs[2*i:2*i+2] for i in range(25)]:
        assert "downhill" in initial.stem
        assert guess.stem[:2] == initial.stem[:2]
        guess_geom = geom_from_xyz_file(guess)
        initial_geom = geom_from_xyz_file(initial)
        N_init = guess_geom.coords - initial_geom.coords
        N_dict[guess.name] = N_init
    return N_dict


# def get_N_10_11_15_dict():
    # pickle_path = Path("10_11_N_init.pickle")
    # with open(pickle_path, "rb") as handle:
        # N_inits = cloudpickle.load(handle)
    # return N_inits


# def run():
    # start = time.time()
    # dimv2_kwargs = {
        # "max_step": 0.25,
        # "R": 0.0189,
        # "max_cycles": 1,
        # "rot_kwargs": {
            # "max_cycles": 15,
            # "alpha": 0.05,
        # }
    # }

    # N_init_dict = make_N_init_dict()
    # N_init_10_11_15 = get_N_10_11_15_dict()
    # N_init_dict.update(N_init_10_11_15)
    # results_list = list()
    # # for f, g in zip(xyzs, geoms):
    # for num in BAKER_DICT.keys():
        # xyz_fn, charge, mult, ref_en = BAKER_DICT[num]
        # geom = geom_from_library(f"baker_ts/{xyz_fn}")
        # f = Path(xyz_fn)
        # print(f.stem)
        # dim_kwargs = dimer_kwargs.copy()
        # dim_kwargs["N_init"] = N_init_dict[num]
        # results = run_geom(f.stem, geom, charge, mult, dimer_kwargs=dim_kwargs)
        # # v2_kwargs = dimv2_kwargs.copy()
        # # v2_kwargs["N"]= N_init_dict[num]
        # # results = run_geomv2(f.stem, geom, charge, mult, dimer_kwargs=v2_kwargs)
        # g0 = results.geom0
        # ts_en = g0.energy
        # en_diff = ref_en - ts_en
        # results_list.append((f.name, ts_en, en_diff, results.converged, results.force_evals))
        # print("Results so far:")
        # for f_, ts_en, den, conved, evals_ in results_list:
            # print(f"\t{f_}: {ts_en:.6f} ({den:.6f}), {conved}, {evals_}")
        # print()
        # print()
        # sys.stdout.flush()
    # results_list.append(dimer_kwargs)

    # dimer_kwargs["time"] = time.time()
    # hash_ = hash(frozenset(dimer_kwargs.items()))
    # pickle_fn = f"results_{hash_}.pickle"
    # with open(pickle_fn, "wb") as handle:
        # cloudpickle.dump(results_list, handle)
    # print(f"Save pickled results to {pickle_fn}.")
    # end = time.time()
    # duration = end - start
    # print(f"Whole run took {duration:.0f} seconds.")


# def run_geomv2(stem, geom, charge, mult, dimer_kwargs):
    # # Gaussian 16
    # calc_kwargs = {
        # "route": "HF/3-21G",
        # "pal": 4,
        # "mem": 1000,
        # "charge": charge,
        # "mult": mult,
    # }
    # def calc_getter():
        # return Gaussian16(**calc_kwargs)

    # geom.set_calculator(calc_getter())

    # dimer_kwargs["calc_getter"] = calc_getter
    # coords = dimer_method_v2(geom, **dimer_kwargs)
    # ts_xyz = geom.as_xyz()
    # ts_fn = f"{stem}_dimer_ts.xyz"
    # with open(ts_fn, "w") as handle:
        # handle.write(ts_xyz)
    # print(f"Wrote dimer result to {ts_fn}")


@using("pyscf")
# @using("gaussian16")
@pytest.mark.parametrize(
    "name, geom, charge, mult, ref_energy",
    [_ for _ in get_baker_ts_geoms_flat()
     if _[0][:2] not in ("10", "11", "15", "17", "20", "18", "05")
     # if _[0][:2] in ("10", "11", "15")  # no cart. imag. freqs.
     # if _[0][:2] in ("17", "20", "18", "05")  # just fails :)
    ]
)
def test_baker_ts_dimer(name, geom, charge, mult, ref_energy):
    print(f"@Got geom '{name}'")
    # Load initial dimers
    N_init_dict = make_N_init_dict()

    calc_kwargs = {
        "charge": charge,
        "mult": mult,
        "pal": 2,
        "base_name": Path(name).stem,
    }
    def calc_getter():
        return PySCF(basis="321g", **calc_kwargs)

    # def calc_getter():
        # return Gaussian16(route="HF/3-21G", **calc_kwargs)
    geom.set_calculator(calc_getter())

    dimer_kwargs = {
        "max_step": 0.25,
        # 1e-2 Angstroem in bohr
        "dR_base": 0.0189,
        "rot_opt": "lbfgs",
        # "trans_opt": "mb",
        "trans_opt": "lbfgs",
        # "trans_memory": 10, # bad idea
        "angle_tol": 5,
        "f_thresh": 1e-3,
        "max_cycles": 50,
        "f_tran_mod": True,
        # "rot_type": "direct",
        "multiple_translations": True,
    }
    dimer_kwargs["N_init"] = N_init_dict[name]
    geoms = (geom, )
    results = dimer_method(geoms, calc_getter, **dimer_kwargs)

    same_energy = geom.energy == pytest.approx(ref_energy)
    print(f"@Same energy: {str(same_energy): >5}, {name}, "
          f"{results.force_evals} force evaluations")
    if not same_energy:
        do_final_hessian(geom)

    # This way pytest prints the actual values... instead of just the boolean
    assert geom.energy == pytest.approx(ref_energy)


@using("pyscf")
@pytest.mark.parametrize(
    "name, geom, charge, mult, ref_energy",
    [_ for _ in get_baker_ts_geoms_flat()
     if _[0][:2] not in ("10", "11", "15", "17", "20", "18", "05")
    ]
)
def test_baker_dimer_new(name, geom, charge, mult, ref_energy):
    print(f"@Got geom '{name}'")
    # Load initial dimers
    N_init_dict = make_N_init_dict()

    calc_kwargs = {
        "charge": charge,
        "mult": mult,
        "pal": 2,
        "base_name": Path(name).stem,
    }
    calc = PySCF("321g", **calc_kwargs)

    dimer_kwargs = {
        "rotation_method": "fourier",
        "calculator": calc,
        "N_init": N_init_dict[name],
        "length": 0.0189,
        "rotation_tol": 5,
    }
    dimer = Dimer(**dimer_kwargs)
    geom.set_calculator(dimer)

    opt_kwargs = {
        "precon": True,
        "max_step_element": 0.25,
        "max_cycles": 50,
        "thresh": "baker",
        "c_stab": 0.103,
    }
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(ref_energy)

    print(f"@{name} converged using {dimer.force_evals} force evaluations")
