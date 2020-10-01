#!/usr/bin/env python3

import itertools as it
from pathlib import Path
from pprint import pprint
import shutil
import time

import numpy as np
import pandas as pd

from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.color import red, green
from pysisyphus.helpers import get_baker_ts_geoms, do_final_hessian, \
                               geom_from_library, get_baker_ts_geoms_flat, \
                               geom_loader
from pysisyphus.intcoords.augment_bonds import augment_bonds
from pysisyphus.tsoptimizers import *


def print_summary(converged, failed, cycles, ran, runid):
    ran_ = f"{ran+1:02d}"
    print(f"converged: {converged:02d}/{ran_}")
    print(f"   failed: {failed:d}")
    print(f"   cycles: {cycles}")
    print(f"      run: {runid}")


def run_baker_ts_opts(geoms, meta, coord_type="cart", thresh="baker", runid=0):
    """From 10.1002/(SICI)1096-987X(199605)17:7<888::AID-JCC12>3.0.CO;2-7"""
    start = time.time()

    converged = 0
    failed = 0
    cycles = 0
    opt_kwargs = {
        "thresh": thresh,
        # "max_cycles": 150,
        "max_cycles": 100,
        # "max_cycles": 50,
        "dump": True,
        "trust_radius": 0.3,
        "trust_max": 0.3,
        # "max_micro_cycles": 1,
    }
    results = dict()
    for i, (name, geom) in enumerate(geoms.items()):
        print(f"@Running {name}")
        charge, mult, ref_energy = meta[name]
        calc_kwargs = {
            "charge": charge,
            "mult": mult,
            "pal": 4,
        }
        geom.set_calculator(Gaussian16(route="HF/3-21G", **calc_kwargs))
        geom = augment_bonds(geom)
        # geom.set_calculator(PySCF(basis="321g", **calc_kwargs))

        # opt = RSPRFOptimizer(geom, **opt_kwargs)
        opt = RSIRFOptimizer(geom, **opt_kwargs)
        # opt = RSIRFOptimizer(geom, **opt_kwargs)
        # opt = TRIM(geom, **opt_kwargs)
        opt.run()
        if opt.is_converged:
            converged += 1
        else:
            failed += 1
        cycles += opt.cur_cycle + 1
        energies_match = np.allclose(geom.energy, ref_energy)
        try:
            assert np.allclose(geom.energy, ref_energy)
            # Backup TS if optimization succeeded
            # ts_xyz_fn = Path(name).stem + "_opt_ts.xyz"
            # out_path = Path("/scratch/programme/pysisyphus/xyz_files/baker_ts_opt/")
            print(green(f"\t@Energies MATCH for {name}! ({geom.energy:.6f}, {ref_energy:.6f})"))
            # with open(out_path / ts_xyz_fn, "w") as handle:
                # handle.write(geom.as_xyz())
        except AssertionError as err:
            print(red(f"\t@Calculated energy {geom.energy:.6f} and reference "
                      f"energy {ref_energy:.6f} DON'T MATCH'."))
        print()
        print_summary(converged & energies_match, failed, cycles, i, runid)
        print()
        results[name] = (opt.cur_cycle + 1, opt.is_converged)
        pprint(results)
        print()
        # do_final_hessian(geom, False)
        # print()

    end = time.time()
    duration = end - start

    print(f"  runtime: {duration:.1f} s")
    print_summary(converged, failed, cycles, i, runid)
    return results, duration, cycles


@pytest.mark.benchmark
@using_gaussian16
def _test_baker_ts_optimizations():
    coord_type = "redund"
    # coord_type = "dlc"
    # coord_type = "cart"
    thresh = "baker"
    runs = 1

    all_results = list()
    durations = list()
    all_cycles = list()
    for i in range(runs):
        geoms, meta = get_baker_ts_geoms(coord_type=coord_type)
        # only = "01_hcn.xyz"
        # only = "24_h2cnh.xyz"
        # only = "15_hocl.xyz"
        # only = "02_hcch.xyz"
        # geoms = {
            # only: geoms[only],
        # }

        fails = (
            "09_parentdieslalder.xyz",
            "12_ethane_h2_abstraction.xyz",
            "22_hconhoh.xyz",
            "17_claisen.xyz",
            "15_hocl.xyz",
        )
        works = (
            "05_cyclopropyl.xyz",
            "08_formyloxyethyl.xyz",
            "14_vinyl_alcohol.xyz",
            "16_h2po4_anion.xyz",
            "18_silyene_insertion.xyz",
            "04_ch3o.xyz",
            "06_bicyclobutane.xyz",
            "07_bicyclobutane.xyz",
            "23_hcn_h2.xyz",
            "01_hcn.xyz",
            "25_hcnh2.xyz",
        )
        math_error_but_works = (
            # [..]/intcoords/derivatives.py", line 640, in d2q_d
            # x99 = 1/sqrt(x93)
            #   ValueError: math domain error
            # ZeroDivison Fix
            "20_hconh3_cation.xyz",
            "24_h2cnh.xyz",
            "13_hf_abstraction.xyz",
            "19_hnccs.xyz",
            "21_acrolein_rot.xyz",
            "03_h2co.xyz",
        )
        alpha_negative = (
            "02_hcch.xyz",
        )
        no_imag = (
            "10_tetrazine.xyz",
            "11_trans_butadiene.xyz",
        )
        only = (
            "18_silyene_insertion.xyz",
            # "21_acrolein_rot.xyz",
            # "22_hconhoh.xyz",
        )
        use = (
            # fails,
            works,
            math_error_but_works,
            # alpha_negative,
            # no_imag,
            # only,
        )
        geoms = {key: geoms[key] for key in it.chain(*use)}

        # geoms = {"05_cyclopropyl.xyz": geoms["05_cyclopropyl.xyz"]}

        results, duration, cycles = run_baker_ts_opts(
                                        geoms,
                                        meta,
                                        coord_type,
                                        thresh,
                                        runid=i
        )
        all_results.append(results)
        durations.append(duration)
        all_cycles.append(cycles)
        print(f"@Run {i}, {cycles} cycles")
        print(f"@All cycles: {all_cycles}")
        print(f"@This runtime: {duration:.1f} s")
        print(f"@Total runtime: {sum(durations):.1f} s")
        print(f"@")
    return

    names = list(results.keys())
    cycles = {
        name: [result[name][0] for result in all_results] for name in names
    }
    df = pd.DataFrame.from_dict(cycles)

    df_name = f"cycles_{coord_type}_{runs}_runs_{thresh}.pickle"
    df.to_pickle(df_name)
    print(f"Pickled dataframe to {df_name}")
    print(f"{runs} runs took {sum(durations):.1f} seconds.")


if __name__ == "__main__":
    _test_baker_ts_optimizations()
