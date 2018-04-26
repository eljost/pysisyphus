#!/usr/bin/env python3

import os
import pathlib
from pathlib import Path

#import numpy as np

from pysisyphus.helpers import geom_from_library, geom_from_xyz_file
from pysisyphus.calculators.Turbomole import Turbomole
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
from pysisyphus.optimizers.SteepestDescent import SteepestDescent

THIS_DIR = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))


def check():
    import re
    import matplotlib.pyplot as plt
    from natsort import natsorted
    import numpy as np
    p = Path(".")
    all_ens = list()
    for log in natsorted(p.glob("calculator*.out")):
        with open(log) as handle:
            text = handle.read()
        ens = [float(e) for e in 
               re.findall("Total energy:\s*([\d\.\-]+)", text)]
        all_ens.append(ens)
    arr = np.array(all_ens, dtype=float)
    fig, ax = plt.subplots()
    for i, row in enumerate(all_ens):
        xs = np.full_like(row, i)
        ax.plot(xs, row, "+")
    plt.show()


def test_butadiene_track_opt():
    in_path = THIS_DIR / "butadiene"
    geom = geom_from_xyz_file(in_path / "butadiene_hf_sto3g.xyz")
    turbo = Turbomole(in_path, track=True, wfo_basis="def2-svp")
    geom.set_calculator(turbo)

    #fn = "/scratch/programme/pysisyphus/tests/test_turbo_butadien_td_opt/wfo_backup.out"
    #with open(fn) as handle:
    #    stdout = handle.read()
    #wfo = turbo.wfow
    #a = wfo.parse_wfoverlap_out(stdout)
    #print(a)
    #print(a.reshape(-1, 6))

    opt_kwargs = {
        "max_cycles": 10,
        "dump": True,
    }
    opt = ConjugateGradient(geom, **opt_kwargs)
    opt = SteepestDescent(geom, **opt_kwargs)
    opt.run()


def test_butadiene_twist_track_opt():
    in_path = THIS_DIR / "butadiene_twist"
    geom = geom_from_xyz_file(in_path / "02_buta_twist.xyz")
    turbo = Turbomole(in_path, track=True, wfo_basis="def2-svp")
    geom.set_calculator(turbo)

    opt_kwargs = {
        "max_cycles": 3,
        "dump": True,
    }
    opt = SteepestDescent(geom, **opt_kwargs)
    opt.run()


if __name__ == "__main__":
    #test_butadiene_track_opt()
    test_butadiene_twist_track_opt()
    check()
