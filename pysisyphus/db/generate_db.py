import argparse
import os
import sys

import numpy as np

from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.helpers import geom_loader, do_final_hessian
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.db import LEVELS, MOLECULES, GUESS_DIR, LEVEL_DIR, OPT_FLAG, THIS_DIR
import pysisyphus.db.helpers as dbhelpers


OPT_KWARGS = {
    "thresh": "gau_tight",
    "max_cycles": 100,
}


def generate_db():
    geoms = [geom_loader(GUESS_DIR / mol.fn, coord_type="redund") for mol in MOLECULES]

    for level in LEVELS:
        level_name, calc_cls, calc_kwargs = level
        try:
            out_dir = LEVEL_DIR / level_name
            os.mkdir(out_dir)
        except FileExistsError:
            pass

        for geom, mol in zip(geoms, MOLECULES):
            print(highlight_text(f"Running '{mol.name}' AT LEVEL '{level_name}'"))

            out_fn = dbhelpers.get_path(mol.name, level_name)
            try:
                with open(out_fn, "r") as handle:
                    _ = handle.readline()
                    line2 = handle.readline()
                if OPT_FLAG in line2:
                    print("Already optimized! Skipping!")
                    continue
            except FileNotFoundError:
                pass

            geom = geom.copy()
            geom.comment = level_name
            calc = calc_cls(**calc_kwargs)
            geom.set_calculator(calc)
            opt = RFOptimizer(geom, **OPT_KWARGS)
            opt.run()

            assert opt.is_converged
            hess_res = do_final_hessian(geom, save_hessian=False)
            neg_eigvals = hess_res.neg_eigvals
            if neg_eigvals.size > 0:
                assert all(
                    np.abs(neg_eigvals) <= 2e-5
                ), f"Calculation at level '{level_name}' failed for {mol.name}"

            geom.standard_orientation()
            geom.comment += f", {OPT_FLAG}"

            with open(out_fn, "w") as handle:
                handle.write(geom.as_xyz())
            print()


def clean():
    fns = "final_geometry.xyz cur_out optimizer.log internal_coords.log".split()
    fns += list(THIS_DIR.glob("calculator_000.*"))
    for fn in fns:
        try:
            os.remove(THIS_DIR / fn)
        except FileNotFoundError:
            print(f"Could not delete {fn}")


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--generate", action="store_true", help="Run all necessary optimizations."
    )
    parser.add_argument("--clean", action="store_true", help="Clean directory")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    if args.generate:
        generate_db()
        clean()


if __name__ == "__main__":
    run()
