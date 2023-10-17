#!/usr/bin/env python3

import argparse
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
import yaml

from pysisyphus.calculators import OpenMolcas
from pysisyphus.constants import AU2EV
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging


init_logging()


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


def run_direction(geoms, indices, inporb, calc_kwargs):
    ngeoms = len(indices)
    print(f"Running direction for {ngeoms} geometries with {indices=}")
    calc_getter = get_calc_getter(calc_kwargs)
    prev_inporb = inporb
    all_results = list()
    for index in indices:
        geom = geoms[index]
        key = f"geom_{index:03d}"
        calc = calc_getter(calc_number=index, inporb=prev_inporb)
        print(f"{key} ... started at {timef()}")
        print(f"{key} ... inporb={prev_inporb}")
        sys.stdout.flush()
        start = time.time()
        results = calc.get_energy(geom.atoms, geom.cart_coords)
        dur = time.time() - start
        print(f"{key} ... finished at {timef()}")
        print(f"{key} ... calculation took {dur: >8.2f} s")
        run_mopicgen(index, calc)
        print(f"{key} ... ran mopicgen")

        all_results.append(results)
        prev_inporb = calc.inporb

    return all_results


def run_calculations(geoms, start_index: int, calc_kwargs: dict, inporb, client=None):
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
    # or in parallel
    else:
        results = list()
        for indices in directions:
            results.append(run_direction(geoms, indices, inporb, calc_kwargs))

    all_results = [None] * ngeoms
    for res, indices in zip(results, directions):
        for i, index in enumerate(indices):
            all_results[index] = res[i]

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
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])
    yaml_fn = args.yaml
    cluster = args.cluster

    loader = yaml.SafeLoader
    with open(yaml_fn) as handle:
        yaml_str = handle.read()
    run_dict = yaml.load(yaml_str, Loader=loader)

    os.environ["MOLCAS"] = run_dict["molcas"]

    trj_fn = run_dict["geom"]["fn"]
    geoms = geom_loader(trj_fn)

    calc_kwargs = run_dict["calc"]
    calc_type = calc_kwargs.pop("type")
    assert calc_type == "openmolcas"
    print("calc_kwargs:")
    pprint(calc_kwargs)
    print()

    inporb = calc_kwargs.pop("inporb")
    start_index = run_dict["start_index"]

    cwd = Path(".")
    calc_kwargs["out_dir"] = cwd / "qm_calcs"

    if cluster:
        cluster_kwargs = {
            "threads_per_worker": 2,
            "n_workers": 1,
        }
        cluster = LocalCluster(**cluster_kwargs)
        link = cluster.dashboard_link
        print(f"Started cluster. Dashboard link is {link}")
        client = Client(cluster)
    else:
        client = None

    start = time.time()
    all_results = run_calculations(geoms, start_index, calc_kwargs, inporb, client)
    dur = time.time() - start
    print(f"Calculations took {dur:.2f} s")
    dump_all_energies(all_results, cwd)


if __name__ == "__main__":
    run()
