#!/usr/bin/env python3

import argparse
from enum import Enum
import os
from pathlib import Path
from pprint import pprint
import shutil
import subprocess
import sys
import time

from distributed import Client, LocalCluster
import matplotlib.pyplot as plt
import numpy as np
import psutil
import yaml

from pysisyphus.calculators import OpenMolcas
from pysisyphus.constants import AU2EV
from pysisyphus.color import red
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging


init_logging()


MODE = Enum("MODE", ("DIRECTIONS", "STANDLONE"))


def timef():
    return time.strftime("%x %X")


def get_calc_getter(calc_kwargs):
    _calc_kwargs = calc_kwargs.copy()
    assert "charge" in _calc_kwargs
    assert "mult" in _calc_kwargs

    def calc_getter(calc_number, inporb):
        calc_kwargs = _calc_kwargs.copy()
        calc_kwargs.update(
            {
                "calc_number": calc_number,
                "inporb": inporb,
                # "rassi": {
                # "mees": None,
                # },
            }
        )
        calc = OpenMolcas(**calc_kwargs)
        return calc

    return calc_getter


def run_mopicgen(index, calc):
    try:
        mop_path = Path(f"{index:03d}")
        os.mkdir(mop_path)
    except FileExistsError:
        pass

    def subproc_run(cmd):
        subprocess.run(cmd, cwd=mop_path, capture_output=True)

    molden = calc.inporb.with_suffix(".rasscf.molden")
    mop_cmd = ("mopicgen", str(molden), "--fracmos", "--hydrogen")
    subproc_run(mop_cmd)

    run_path = str((mop_path / "run.sh").absolute())
    run_cmd = ("bash", run_path)
    subproc_run(run_cmd)
    shutil.copy(mop_path / "montage.0.jpg", f"montage.{index:03d}.jpg")


def run_calculation(geom, index, calc):
    key = f"geom_{index:03d}"
    print(f"{key} ... started at {timef()}")
    print(f"{key} ... inporb={calc.inporb}")
    sys.stdout.flush()
    start = time.time()
    results = calc.get_energy(geom.atoms, geom.cart_coords)
    dur = time.time() - start
    print(red(f"{key} ... finished at {timef()}"))
    print(f"{key} ... calculation took {dur: >8.2f} s")
    run_mopicgen(index, calc)
    print(f"{key} ... ran mopicgen")
    return results


def run_direction(geoms, indices, inporb, calc_kwargs):
    ngeoms = len(indices)
    print(f"Running direction for {ngeoms} geometries with {indices=}")
    calc_getter = get_calc_getter(calc_kwargs)
    prev_inporb = inporb
    all_results = list()
    for index in indices:
        calc = calc_getter(calc_number=index, inporb=prev_inporb)
        geom = geoms[index]
        results = run_calculation(geom, index, calc)
        all_results.append(results)
        prev_inporb = calc.inporb

    return all_results


def run_directions(geoms, start_index: int, inporb, calc_kwargs: dict, client=None):
    ngeoms = len(geoms)
    end_index = ngeoms - 1
    assert 0 <= start_index < ngeoms

    directions = list()
    if start_index == 0:
        directions.append(list(range(ngeoms)))
    elif start_index == end_index:
        directions.append(list(range(ngeoms))[::-1])
    else:
        directions.extend(
            [
                list(range(0, start_index))[::-1],
                list(range(start_index, ngeoms)),
            ]
        )

    def wrapper(indices):
        return run_direction(geoms, indices, inporb, calc_kwargs)

    # Carry out calculations
    #
    # in parallel
    if client is not None:
        futures = client.map(wrapper, directions)
        results = client.gather(futures)
    # or in serial
    else:
        results = list()
        for indices in directions:
            results.append(run_direction(geoms, indices, inporb, calc_kwargs))

    all_results = [None] * ngeoms
    for res, indices in zip(results, directions):
        for i, index in enumerate(indices):
            all_results[index] = res[i]

    return all_results


def run_calculations(geoms, indices, inporbs, calc_kwargs, client):
    calc_getter = get_calc_getter(calc_kwargs)
    all_results = list()

    def wrapper(i):
        inporb = inporbs[i]
        geom = geoms[i]
        calc = calc_getter(calc_number=i, inporb=inporb)
        results = run_calculation(geom, i, calc)
        return results

    if client is not None:
        futures = client.map(wrapper, indices)
        all_results = client.gather(futures)
    else:
        all_results = list(map(wrapper, indices))
    return all_results


def dump_all_energies(all_results, cwd):
    # Dump all energies
    all_energies = [res["all_energies"] for res in all_results]
    all_energies = np.array(all_energies)
    np.save(cwd / "all_energies.npy", all_energies)

    indices = np.arange(all_energies.shape[0])

    # Plot PEC
    all_energies -= all_energies.min()
    all_energies *= AU2EV
    fig, ax = plt.subplots()
    for i, state in enumerate(all_energies.T):
        label = f"root {i}"
        ax.plot(indices, state, "o-", label=label)
    ax.set_xlabel("Index")
    ax.set_ylabel("ΔE / eV")
    ax.legend()
    fig.tight_layout()
    pdf_fn = cwd / "all_energies.pdf"
    fig.savefig(pdf_fn)
    print(f"Dumped PEC to '{pdf_fn}'")


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("yaml")
    parser.add_argument("--cluster", action="store_true")
    parser.add_argument("--nthreads", type=int, default=psutil.cpu_count(logical=False))
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])
    yaml_fn = args.yaml
    cluster = args.cluster
    nthreads = args.nthreads

    loader = yaml.SafeLoader
    with open(yaml_fn) as handle:
        yaml_str = handle.read()
    run_dict = yaml.load(yaml_str, Loader=loader)

    os.environ["MOLCAS"] = run_dict["molcas"]

    trj_fn = run_dict["geom"]["fn"]
    geoms = geom_loader(trj_fn)
    ngeoms = len(geoms)

    calc_kwargs = run_dict["calc"]
    calc_type = calc_kwargs.pop("type")
    assert calc_type == "openmolcas"
    print("calc_kwargs:")
    pprint(calc_kwargs)
    print()

    start_index = run_dict.get("start_index", None)
    inporb_glob = run_dict.get("inporb_glob", None)

    cwd = Path(".")
    calc_kwargs["out_dir"] = cwd / "qm_calcs"

    if inporb_glob is not None:
        inporbs = list(sorted(Path(".").glob(inporb_glob)))
        assert ngeoms == len(inporbs)
        indices = np.arange(ngeoms)
    elif start_index is not None:
        # TODO: Give inporb in calc section or outside, next to start_index?
        inporb = calc_kwargs.pop("inporb")
        inporbs = [inporb]
        indices = [start_index]
    else:
        raise Exception("Either start_index or inporb_glob must be present!")

    mode = (
        MODE.DIRECTIONS if (len(indices) == 1 and len(inporbs) == 1) else MODE.STANDLONE
    )

    if cluster:
        n_workers = 1
        # TODO: threads depend on run_mode
        threads_per_worker = 2 if mode == MODE.DIRECTIONS else nthreads
        cluster_kwargs = {
            "threads_per_worker": threads_per_worker,
            "n_workers": n_workers,
        }
        cluster = LocalCluster(**cluster_kwargs)
        link = cluster.dashboard_link
        print(
            f"Started cluster with {n_workers=} and {threads_per_worker=}. "
            f"Dashboard link is {link}"
        )
        client = Client(cluster)
    else:
        client = None

    start = time.time()
    # Direction mode, where active space orbitals are propgagated along a coordinate
    if mode == MODE.DIRECTIONS:
        all_results = run_directions(geoms, indices[0], inporbs[0], calc_kwargs, client)
    # Parallel mode, where each geometry already has active space orbitals and
    # calculations can run in parallel.
    elif mode == MODE.STANDLONE:
        all_results = run_calculations(geoms, indices, inporbs, calc_kwargs, client)
    else:
        raise Exception("How did I end up here?")
    dur = time.time() - start
    print(f"Calculations took {dur:.2f} s")

    dump_all_energies(all_results, cwd)


if __name__ == "__main__":
    run()
