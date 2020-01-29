#!/usr/bin/env python3

from pprint import pprint
import time

import numpy as np
import pandas as pd
import pytest

from pysisyphus.helpers import get_baker_geoms
from pysisyphus.color import red, green
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
# from pysisyphus.calculators.XTB import XTB
# from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.testing import using_pyscf


def print_summary(converged, failed, cycles, ran, runid):
    ran_ = f"{ran+1:02d}"
    print(f"converged: {converged:02d}/{ran_}")
    print(f"   failed: {failed:02d}/{ran_}")
    print(f"   cycles: {cycles}")
    print(f"      run: {runid}")


def run_baker_opts(baker_geoms, coord_type="cart", thresh="gau_tight",
                   gediis=False, gdiis=False, poly=False, runid=0):
    """From https://onlinelibrary.wiley.com/doi/epdf/10.1002/jcc.540140910"""
    start = time.time()

    converged = 0
    failed = 0
    cycles = 0
    opt_kwargs = {
        "thresh": thresh,
        "max_cycles": 150,
        "line_search": poly,
        "gediis": gediis,
        "gdiis": gdiis,
        # "overachieve_factor": 2,
        # "dump": True,
    }
    results = dict()
    for i, (name, (geom, refen)) in enumerate(baker_geoms.items()):
        print(f"@Running {name}")
        # geom.set_calculator(XTB(pal=4))
        # geom.set_calculator(Gaussian16(route="HF/STO-3G", pal=4))
        geom.set_calculator(PySCF(basis="sto3g", pal=2))
        opt = RFOptimizer(geom, **opt_kwargs)
        opt.run()
        if opt.is_converged:
            converged += 1
        else:
            failed += 1
        cycles += opt.cur_cycle + 1
        try:
            assert geom.energy == pytest.approx(refen)
            print(green(f"@Energies match for {name}! ({geom.energy:.6f}, {refen:.6f})"))
        except KeyError:
            print(red(f"@Couldn't find reference energy for {name}"))
        except AssertionError as err:
            print(red(f"@Calculated energy {geom.energy:.6f} and reference "
                      f"energy {refen:.6f} don't match."))
        print()
        print_summary(converged, failed, cycles, i, runid)
        print()
        results[name] = (opt.cur_cycle + 1, opt.is_converged)
        pprint(results)
        print()

    end = time.time()
    duration = end - start

    print(f"  runtime: {duration:.1f} s")
    print_summary(converged, failed, cycles, i, runid)
    return results, duration, cycles


def run_baker_minimum_optimizations():
    coord_type = "redund"
    # coord_type = "cart"
    gediis = False
    gdiis = True
    poly = True
    # thresh = "gau_tight"
    thresh = "baker"
    runs = 1


    all_results = list()
    durations = list()
    all_cycles = list()
    for i in range(runs):
        baker_geoms = get_baker_geoms(coord_type=coord_type)

        results, duration, cycles = run_baker_opts(
                                        baker_geoms,
                                        coord_type,
                                        thresh,
                                        gediis,
                                        gdiis,
                                        poly,
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
        break
    return

    names = list(results.keys())
    cycles = {
        name: [result[name][0] for result in all_results] for name in names
    }
    df = pd.DataFrame.from_dict(cycles)

    ged = "_gediis" if gediis else ""
    gd = "_gdiis" if gdiis else ""
    pod = "_poly" if poly else ""
    df_name = f"cycles_{coord_type}_{runs}_runs_{thresh}{ged}{gd}{pod}.pickle"
    df.to_pickle(df_name)
    print(f"Pickled dataframe to {df_name}")
    print(f"{runs} runs took {sum(durations):.1f} seconds.")


@using_pyscf
@pytest.mark.parametrize(
    "name, geom, ref_energy",
    [(name, geom, ref_energy) for name, (geom, ref_energy)
     in get_baker_geoms(coord_type="redund").items()]
)
def test_baker_gs_opt(name, geom, ref_energy):
    opt_kwargs = {
        "thresh": "baker",
    }
    print(f"@Running {name}")
    geom.set_calculator(PySCF(basis="sto3g", pal=2))
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()
    assert np.allclose(geom.energy, ref_energy)


if __name__ == "__main__":
    run_baker_minimum_optimizations()
