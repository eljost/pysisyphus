import os
from pathlib import Path
import shutil

import numpy as np

from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.helpers import geom_loader, do_final_hessian, highlight_text
from pysisyphus.db import LEVELS, MOLECULES


OPT_KWARGS = {
    "thresh": "gau_tight",
    "max_cycles": 100,
}
GUESS_DIR = Path("./guess")
LEVEL_DIR = Path("./levels")


def generate_db():
    geoms = [geom_loader(GUESS_DIR / mol.fn, coord_type="redund") for mol in MOLECULES]

    for level in LEVELS:
        level_name, calc_cls, calc_kwargs = level
        try:
            out_dir = LEVEL_DIR / level_name
            os.mkdir(out_dir)
        except FileExistsError:
            shutil.rmtree(out_dir)
            os.mkdir(out_dir)

        for geom, mol in zip(geoms, MOLECULES):
            geom = geom.copy()
            print(highlight_text(f"Running '{mol.name}' AT LEVEL '{level_name}'"))
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

            out_fn = f"{mol.name}.xyz"
            with open(out_dir / out_fn, "w") as handle:
                handle.write(geom.as_xyz())
            print()


def clean():
    fns = "final_geometry.xyz cur_out optimizer.log internal_coords.log".split()
    for fn in fns:
        try:
            os.remove(fn)
        except FileNotFoundError:
            print(f"Could not delete {fn}")


if __name__ == "__main__":
    try:
        generate_db()
        clean()
    except AssertionError as err:
        raise (err)
