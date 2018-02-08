#!/usr/bin/env python3

import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from pysisyphus.constants import BOHR2ANG


class Calculator:
    logger = logging.getLogger("calculator")

    def __init__(self, calc_number=0, charge=0, mult=1,
                 base_name="calculator", last_calc_cycle=None):
        self.charge = int(charge)
        self.mult = int(mult)
        # Index of the image this calculator belongs too in
        # in a ChainOfStates calculation.
        self.calc_number = calc_number
        self.base_name = base_name
        self.name = f"{base_name}_{calc_number}"

        # Extensions of the files to keep after running a calculation.
        # Usually overridden in derived classes.
        self.to_keep = ()
        # How many calculations were already run
        self.calc_counter = 0
        # Handle restarts
        if last_calc_cycle:
            self.calc_counter = int(last_calc_cycle)+1
            self.reattach(int(last_calc_cycle))
            self.log(f"set {self.calc_counter} for this calculation")

        self.inp_fn = "calc.inp"
        self.out_fn = "calc.out"
        # When this is set the run() method will use this path
        # instead of creating a new one.
        # Currently this is only used with the Turbomole calculator.
        self.path_already_prepared = None

    def reattach(self, last_calc_cycle):
        raise Exception("Not implemented!")

    def log(self, message):
        self.logger.debug(f"{self.name}_cyc_{self.calc_counter:03d}, "
                          + message)

    def get_energy(self, atoms, coords):
        raise Exception("Not implemented!")

    def get_forces(self, atoms, coords):
        raise Exception("Not implemented!")

    def get_hessian(self, atoms, coords):
        raise Exception("Not implemented!")

    def make_fn(self, ext, counter=None, abspath=False):
        if not counter:
            counter = self.calc_counter
        fn = f"{self.name}.{counter:03d}.{ext}"
        if abspath:
            fn = os.path.abspath(fn)
        return fn

    def prepare_path(self, use_in_run=False):
        """When use_in_run is True run() will not create a new path
        but reuse this path instead."""
        prefix = f"{self.name}_{self.calc_counter:03d}_"
        path = Path(tempfile.mkdtemp(prefix=prefix))
        if use_in_run:
            self.path_already_prepared = path
        return path

    def prepare(self, inp):
        if not self.path_already_prepared:
            path = self.prepare_path()
        else:
            path = self.path_already_prepared

        # Calculators like Turbomole got no input.
        if inp:
            inp_path = path / self.inp_fn
            with open(inp_path, "w") as handle:
                handle.write(inp)

        return path

    def prepare_coords(self, atoms, coords):
        """Convert Bohr to Angstrom."""
        coords = coords.reshape(-1, 3) * BOHR2ANG
        coords = "\n".join(
                ["{} {:10.08f} {:10.08f} {:10.08f}".format(a, *c) for a, c in zip(atoms, coords)]
        )
        return coords

    def run_after(self, path):
        pass

    def run(self, inp, calc, add_args=None, env=None, shell=False, hold=False):
        path = self.prepare(inp)
        self.log(f"running in {path}")
        args = [self.base_cmd, self.inp_fn]
        if add_args:
            args.extend(add_args)
        if not env:
            env = os.environ.copy()
        with open(path / self.out_fn, "w") as handle:
            result = subprocess.Popen(args, cwd=path,
                                      stdout=handle, stderr=subprocess.PIPE,
                                      env=env, shell=shell)
            result.wait()
        try:
            self.run_after(path)
            results = self.parser_funcs[calc](path)
            self.keep(path)
        except Exception as err:
            print(err)
            print()
            print("Crashed input:")
            print(inp)
            backup_dir = Path(os.getcwd()) / f"crashed_{self.name}"
            if backup_dir.exists():
                shutil.rmtree(backup_dir)
            shutil.copytree(path, backup_dir)
            sys.exit()
        finally:
            if not hold:
                self.clean(path)
                self.calc_counter += 1

        self.path_already_prepared = None
        return results

    def keep(self, path):
        kept_fns = dict()
        for ext in self.to_keep:
            pattern = f"*{ext}"
            globbed = list(path.glob(pattern))
            assert(len(globbed) <= 1), (f"Expected at most one file ending with "
                f"{pattern} in {path}. Found {len(globbed)} instead!"
            )
            if len(globbed) == 0:
                continue
            old_fn = globbed[0]
            new_fn = os.path.abspath(self.make_fn(ext))
            shutil.copy(old_fn, new_fn)
            kept_fns[ext] = new_fn
        return kept_fns

    def clean(self, path):
        shutil.rmtree(path)
        self.log(f"cleaned {path}")

    def reattach(self, calc_counter):
        self.calc_counter = calc_counter


if __name__ == "__main__":
    calc = Calculator()
    calc.base_cmd = "sleep"
    calc.run("dummy_inp", "dummy_calc_type")
