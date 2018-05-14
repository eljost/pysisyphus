#!/usr/bin/env python3

import os
import pathlib
from pathlib import Path

import cloudpickle
from natsort import natsorted
import numpy as np

from pysisyphus.helpers import geom_from_library, geom_from_xyz_file, geoms_from_trj
from pysisyphus.calculators.Turbomole import Turbomole
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
from pysisyphus.optimizers.SteepestDescent import SteepestDescent


THIS_DIR = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
np.set_printoptions(suppress=True, precision=4)


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


def test_butadiene_twist_track_cc2_opt():
    in_path = THIS_DIR / "butadiene_twist_cc2"
    geom = geom_from_xyz_file(in_path / "02_buta_twist.xyz")
    turbo = Turbomole(in_path, track=True, wfo_basis="def2-svp")
    geom.set_calculator(turbo)

    opt_kwargs = {
        "max_cycles": 3,
        "dump": True,
    }
    opt = SteepestDescent(geom, **opt_kwargs)
    opt.run()


def test_wfo_ref():
    in_path = THIS_DIR / "butadiene_twist"
    geom = geom_from_xyz_file(in_path / "02_buta_twist.xyz")
    turbo = Turbomole(in_path, track=True, wfo_basis="def2-svp")
    geom.set_calculator(turbo)

    opt_kwargs = {
        "max_cycles": 2,
        "dump": True,
    }
    opt = SteepestDescent(geom, **opt_kwargs)
    opt.run()

    #wfo_ref = "wfo.ref"
    #with open(wfo_ref) as handle:
    #    text = handle.read()

    #wfow = turbo.wfow
    #for key in wfow.matrix_types:
    #    mat = wfow.parse_wfoverlap_out(text, key)
    #    print(mat)


def test_wfo_compare():
    in_path = THIS_DIR / "butadiene_twist"
    geom = geom_from_xyz_file(in_path / "02_buta_twist.xyz")
    turbo = Turbomole(in_path, track=True, wfo_basis="def2-svp")
    geom.set_calculator(turbo)

    opt_kwargs = {
        "max_cycles": 1,
        "dump": True,
    }
    opt = SteepestDescent(geom, **opt_kwargs)
    opt.run()
    wfow = turbo.wfow
    wfow.compare(wfow)


def test_wfo_compare_butadiene_cc2_sto3g():
    in_path = THIS_DIR / "butadiene_cc2_sto3g"
    geom = geom_from_xyz_file(in_path / "buta.xyz")
    turbo = Turbomole(in_path, track=True, wfo_basis="sto3g")
    geom.set_calculator(turbo)

    # opt_kwargs = {
        # "max_cycles": 1,
        # "dump": True,
    # }
    # opt = SteepestDescent(geom, **opt_kwargs)
    # opt.run()
    turbo.run_calculation(geom.atoms, geom.coords)
    wfow = turbo.wfow
    wfow.compare(wfow)

def test_wfo_compare_butadiene_cc2():
    in_path = THIS_DIR / "butadiene_cc2"
    geom = geom_from_xyz_file(in_path / "buta.xyz")
    turbo = Turbomole(in_path, track=True, wfo_basis="def2-svp")
    geom.set_calculator(turbo)

    opt_kwargs = {
        "max_cycles": 1,
        "dump": True,
    }
    opt = SteepestDescent(geom, **opt_kwargs)
    opt.run()
    wfow = turbo.wfow
    wfow.compare(wfow)

def test_wfo_compare_neon():
    in_path = THIS_DIR / "neon"
    geom = geom_from_xyz_file(in_path / "neon.xyz")
    turbo = Turbomole(in_path, track=True, wfo_basis="def2-svp")
    geom.set_calculator(turbo)

    opt_kwargs = {
        "max_cycles": 1,
        "dump": True,
    }
    opt = SteepestDescent(geom, **opt_kwargs)
    opt.run()
    wfow = turbo.wfow
    wfow.compare(wfow)


def test_wfo_compare_neon_dimer():
    in_path = THIS_DIR / "neon_dimer"
    geom = geom_from_xyz_file(in_path / "neon_dimer.xyz")
    turbo = Turbomole(in_path, track=True, wfo_basis="def2-svp")
    geom.set_calculator(turbo)

    opt_kwargs = {
        "max_cycles": 5,
        "dump": True,
        "track": True,
        # "convergence": {
            # "max_force_thresh": 2.5e-8,
        # }
    }
    opt = SteepestDescent(geom, **opt_kwargs)
    #import pdb; pdb.set_trace()
    opt.run()
    # wfow = turbo.wfow
    # wfow.compare(wfow)


def test_wfo_compare_sto3g():
    in_path = THIS_DIR / "butadiene_twist_sto3g"
    geom = geom_from_xyz_file(in_path / "02_buta_twist.xyz")
    turbo = Turbomole(in_path, track=True, wfo_basis="sto-3g")
    geom.set_calculator(turbo)

    opt_kwargs = {
        "max_cycles": 1,
        "dump": True,
    }
    opt = SteepestDescent(geom, **opt_kwargs)
    opt.run()
    wfow = turbo.wfow
    wfow.compare(wfow)


def test_diabatize():
    geoms = geoms_from_trj(THIS_DIR / "ma_proton_transfer/interpolated.trj")
    in_path = THIS_DIR / "ma_turbo"
    calc_kwargs = {
        "track": True,
        "wfo_basis": "sto-3g",
    }
    turbos = list()
    wfos = list()
    for i, geom in enumerate(geoms):
        pickle_fn = f"wfo_pickle_{i}"
        turbo = Turbomole(in_path, calc_number=i, **calc_kwargs)
        geom.set_calculator(turbo)
        forces = geom.forces
        print(f"cycle {i}")
        turbos.append(turbo)
        wfos.append(turbo.wfow)
        with open(pickle_fn, "wb") as handle:
            cloudpickle.dump(turbo.wfow, handle)


def test_diabatize():
    geoms = geoms_from_trj(THIS_DIR / "ma_proton_transfer/interpolated.trj")[:4]
    #in_path = THIS_DIR / "ma_turbo"
    in_path = THIS_DIR / "ma_turbo_no_exopt"
    calc_kwargs = {
        "track": True,
        "wfo_basis": "sto-3g",
    }
    
    #geoms = geoms_from_trj(THIS_DIR / "biaryl_trj/biaryl_first_13.trj")[:4]
    #in_path = THIS_DIR / "biaryl"
    #calc_kwargs = {
    #    "track": True,
    #    "wfo_basis": "def2-svp",
    #}

    turbos = list()
    wfos = list()
    mos = None
    np.set_printoptions(precision=2, suppress=True)
    # for i, geom in enumerate(geoms):
        # pickle_fn = f"wfo_pickle_{i}"
        # turbo = Turbomole(in_path, calc_number=i, **calc_kwargs)
        # if i > 0:
            # turbo.mos = mos
        # none_ = turbo.get_tddft(geom.atoms, geom.coords)
        # mos = turbo.mos
        # print(f"cycle {i}")
        # turbos.append(turbo)
        # wfos.append(turbo.wfow)
        # with open(pickle_fn, "wb") as handle:
            # cloudpickle.dump(turbo.wfow, handle)
        # if i > 0:
            # wfo1 = wfos[i-1]
            # inds = wfo1.compare(turbo.wfow)
            # print(inds)

    mos = natsorted(THIS_DIR.glob("*.mos"))
    td_vecs = natsorted(THIS_DIR.glob("*.ciss_a"))
    for i, (geom, mos, td_vec) in enumerate(zip(geoms, mos, td_vecs)):
        print("cycle", i, geom, mos, td_vec)
        turbo = Turbomole(in_path, calc_number=i, **calc_kwargs)
        turbo.mos = mos
        turbo.td_vec_fn = td_vec
        turbo.check_for_root_flip(geom.atoms, geom.coords)
        wfos.append(turbo.wfow)
        if i > 0:
            wfo1 = wfos[i-1]
            inds = wfo1.compare(turbo.wfow)


def diabatize_pickled():
    wfos = list()
    for fn in natsorted(THIS_DIR.glob("wfo_pickle_*")):
        print(fn)
        with open(fn, "rb") as handle:
            wfo = cloudpickle.load(handle)
        print(wfo)
        wfos.append(wfo)
    print("Found {len(wfos)} WFOWrapper pickles.")

    diabats = list()
    for i, wfo1 in enumerate(wfos[:-1]):
        wfo2 = wfos[i+1]
        max_ovlp_inds = wfo1.compare(wfo2)
        diabats.append(max_ovlp_inds)
    """
    dia_arr: shape (no. of images-1, no. of states)
     The N-th row contains the roots with maximum overlap between the
    N-th and (N+1)-th WFOWrapper objects.
     The index i of item j in the k-th row gives the root in WFO1 while
    j gives the root with maximum overlap in WFO2.
     Transposing the whole array makes this easier to understand. Now
     every rows hold the overlaps for one adiabatic state.
    """
    dia_arr = np.array(diabats)
    print(dia_arr)
    print()
    print(dia_arr.T)


# def diabatize_pickled():
    # fn = "wfo_pickle_0"
    # with open(fn, "rb") as handle:
        # wfo = cloudpickle.load(handle)
    # print(wfo)
    # inds = wfo.compare(wfo)
    # print(inds)
    # fn = "wfo_0.001.out"
    # with open(fn) as handle:
        # text = handle.read()
    # wfo.parse_wfoverlap_out(text)


if __name__ == "__main__":
    #test_butadiene_track_opt()
    #test_butadiene_twist_track_opt()
    #test_butadiene_twist_track_cc2_opt()
    #check()
    #test_wfo_ref()
    #test_wfo_compare()
    #print()
    # test_wfo_compare_neon()
    #test_wfo_compare_neon_dimer()
    #test_wfo_compare_butadiene_cc2()
    test_wfo_compare_butadiene_cc2_sto3g()
    #test_wfo_compare_sto3g()
    #test_diabatize()
    #diabatize_pickled()
