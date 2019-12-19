#!/usr/bin/env python3

from pathlib import Path
import os
import sys
import time

import cloudpickle
from natsort import natsorted
import numpy as np

from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.calculators.XTB import XTB
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_xyz_file, geom_from_library
from pysisyphus.tsoptimizers.dimer import dimer_method
from pysisyphus.tsoptimizers.dimerv2 import dimer_method as dimer_method_v2


def make_N_init_dict():
    xyz_path = Path("/scratch/programme/pysisyphus/xyz_files/baker_ts")
    xyzs = natsorted(xyz_path.glob("*.xyz"))
    N_dict = dict()
    for guess, initial in [xyzs[2*i:2*i+2] for i in range(25)]:
        assert "downhill" in initial.stem
        assert guess.stem[:2] == initial.stem[:2]
        num = initial.stem[:2]
        guess_geom = geom_from_xyz_file(guess)
        initial_geom = geom_from_xyz_file(initial)
        N_init = guess_geom.coords - initial_geom.coords
        N_dict[num] = N_init
    return N_dict


def get_N_10_11_15_dict():
    pickle_path = Path("10_11_N_init.pickle")
    with open(pickle_path, "rb") as handle:
        N_inits = cloudpickle.load(handle)
    return N_inits


def run():
    start = time.time()
    np.random.seed(20180325)
    xyz_path = Path("/scratch/Code/parsezmat/")
    xyzs = natsorted(xyz_path.glob("*.xyz"))
    geoms = [geom_from_xyz_file(fn) for fn in xyzs]

    BAKER_DICT = {
        "01": ("01_hcn.xyz", 0, 1, -92.24604),
        "02": ("02_hcch.xyz", 0, 1, -76.29343),
        "03": ("03_h2co.xyz", 0, 1, -113.05003),
        "04": ("04_ch3o.xyz", 0, 2, -113.69365),
        "05": ("05_cyclopropyl.xyz", 0, 2, -115.72100),
        "06": ("06_bicyclobutane.xyz", 0, 1, -153.90494),
        "07": ("07_bicyclobutane.xyz", 0, 1, -153.89754),
        "08": ("08_formyloxyethyl.xyz", 0, 2, -264.64757),
        "09": ("09_parentdieslalder.xyz", 0, 1, -231.60321),
        # "10": ("10_tetrazine.xyz", 0, 1, -292.81026),
        # "11": ("11_trans_butadiene.xyz", 0, 1, -154.05046),
        "12": ("12_ethane_h2_abstraction.xyz", 0, 1, -78.54323),
        "13": ("13_hf_abstraction.xyz", 0, 1, -176.98453),
        "14": ("14_vinyl_alcohol.xyz", 0, 1, -151.91310),
        # "15": ("15_hocl.xyz", 0, 1, -569.89752),
        "16": ("16_h2po4_anion.xyz", -1, 1, -637.92388),
        "17": ("17_claisen.xyz", 0, 1, -267.23859),
        "18": ("18_silyene_insertion.xyz", 0, 1, -367.20778),
        "19": ("19_hnccs.xyz", 0, 1, -525.43040),
        "20": ("20_hconh3_cation.xyz", +1, 1, -168.24752),
        "21": ("21_acrolein_rot.xyz", 0, 1, -189.67574),
        "22": ("22_hconhoh.xyz", 0, 1, -242.25529),
        "23": ("23_hcn_h2.xyz", 0, 1, -93.31114),
        "24": ("24_h2cnh.xyz", 0, 1, -93.33296),
        "25": ("25_hcnh2.xyz", 0, 1, -93.28172),
    }

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

    dimv2_kwargs = {
        "max_step": 0.25,
        "R": 0.0189,
        "max_cycles": 1,
        "rot_kwargs": {
            "max_cycles": 15,
            "alpha": 0.05,
        }
    }

    N_init_dict = make_N_init_dict()
    N_init_10_11_15 = get_N_10_11_15_dict()
    N_init_dict.update(N_init_10_11_15)
    results_list = list()
    # for f, g in zip(xyzs, geoms):
    for num in BAKER_DICT.keys():
        xyz_fn, charge, mult, ref_en = BAKER_DICT[num]
        geom = geom_from_library(f"baker_ts/{xyz_fn}")
        f = Path(xyz_fn)
        print(f.stem)
        dim_kwargs = dimer_kwargs.copy()
        dim_kwargs["N_init"] = N_init_dict[num]
        results = run_geom(f.stem, geom, charge, mult, dimer_kwargs=dim_kwargs)
        # v2_kwargs = dimv2_kwargs.copy()
        # v2_kwargs["N"]= N_init_dict[num]
        # results = run_geomv2(f.stem, geom, charge, mult, dimer_kwargs=v2_kwargs)
        g0 = results.geom0
        ts_en = g0.energy
        en_diff = ref_en - ts_en
        results_list.append((f.name, ts_en, en_diff, results.converged, results.force_evals))
        print("Results so far:")
        for f_, ts_en, den, conved, evals_ in results_list:
            print(f"\t{f_}: {ts_en:.6f} ({den:.6f}), {conved}, {evals_}")
        print()
        print()
        sys.stdout.flush()
    results_list.append(dimer_kwargs)

    dimer_kwargs["time"] = time.time()
    hash_ = hash(frozenset(dimer_kwargs.items()))
    pickle_fn = f"results_{hash_}.pickle"
    with open(pickle_fn, "wb") as handle:
        cloudpickle.dump(results_list, handle)
    print(f"Save pickled results to {pickle_fn}.")
    end = time.time()
    duration = end - start
    print(f"Whole run took {duration:.0f} seconds.")


def run_geom(stem, geom, charge, mult, dimer_kwargs=None):
    # Gaussian 16
    calc_kwargs = {
        "route": "HF/3-21G",
        "pal": 4,
        "mem": 1000,
        "charge": charge,
        "mult": mult,
    }
    def calc_getter():
        return Gaussian16(**calc_kwargs)

    # # XTB
    # calc_kwargs = {
        # "pal": 4,
        # "charge": charge,
        # "mult": mult,
    # }
    # def calc_getter():
        # return XTB(**calc_kwargs)

    geom.set_calculator(calc_getter())
    geoms = [geom, ]

    if dimer_kwargs is None:
        dimer_kwargs = {
            "max_step": 0.5,
            # 1e-2 Angstroem
            "dR_base": 0.0189,
            "rot_opt": "lbfgs",
            "angle_tol": 5,
            "f_thresh": 1e-3,
            "max_cycles": 100,
        }
    results = dimer_method(geoms, calc_getter, **dimer_kwargs)
    geom0 = results.geom0
    ts_xyz = geom0.as_xyz()
    ts_fn = f"{stem}_dimer_ts.xyz"
    with open(ts_fn, "w") as handle:
        handle.write(ts_xyz)
    print(f"Wrote dimer result to {ts_fn}")
    return results


def run_geomv2(stem, geom, charge, mult, dimer_kwargs):
    # Gaussian 16
    calc_kwargs = {
        "route": "HF/3-21G",
        "pal": 4,
        "mem": 1000,
        "charge": charge,
        "mult": mult,
    }
    def calc_getter():
        return Gaussian16(**calc_kwargs)

    geom.set_calculator(calc_getter())

    dimer_kwargs["calc_getter"] = calc_getter
    coords = dimer_method_v2(geom, **dimer_kwargs)
    ts_xyz = geom.as_xyz()
    ts_fn = f"{stem}_dimer_ts.xyz"
    with open(ts_fn, "w") as handle:
        handle.write(ts_xyz)
    print(f"Wrote dimer result to {ts_fn}")


if __name__ == "__main__":
    run()
