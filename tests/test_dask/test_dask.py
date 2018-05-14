#!/usr/bin/env python

import argparse
import sys

from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_xyz_file
from pysisyphus.cos.NEB import NEB

from distributed import LocalCluster
import numpy as np

np.set_printoptions(suppress=True, precision=4)

def run_distributed(dask_cluster):
    atoms = ("H", "H")
    geoms = list()
    for i in range(7):
        bond_length = 0.8+i*0.2
        print(f"{i}: {bond_length:.02f}")
        coords = np.array((0., 0., 0., 0., 0., bond_length))
        geom = Geometry(atoms, coords)
        kwargs = {
            "keywords": "BP86 def2-TZVP",
            "charge": 0,
            "mult": 1,
            "calc_number": i,
            "blocks": "%tddft nroots 2 iroot 1 end",
            "track": True,
        }
        orca = ORCA(**kwargs)
        geom.set_calculator(orca)
        geoms.append(geom)

    neb_kwargs = {
        "dask_cluster": dask_cluster,
    }
    neb = NEB(geoms, **neb_kwargs)
    forces = neb.forces
    for f in forces.reshape(-1, 6):
        print(f, f"{np.linalg.norm(f):.2f}")


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("scheduler")
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])
    run_distributed(args.scheduler)


if __name__ == "__main__":
    run()
