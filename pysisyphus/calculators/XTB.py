#!/usr/bin/env python3

import glob
import os
import re

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.calculators.parser import parse_turbo_gradient
from pysisyphus.constants import BOHR2ANG
from pysisyphus.calculators.parser import parse_turbo_gradient
from pysisyphus.helpers import geom_from_xyz_file
from pysisyphus.xyzloader import make_xyz_str


class XTB(Calculator):

    conf_key = "xtb"

    def __init__(self, gbsa="", gfn=2, acc=1.0, **kwargs):
        super(XTB, self).__init__(**kwargs)

        self.gbsa = gbsa
        self.gfn = gfn
        self.acc = acc

        valid_gfns = (1, 2)
        assert self.gfn in valid_gfns, "Invalid gfn argument. " \
            f"Allowed arguments are: {', '.join(valid_gfns)}!"
        self.uhf = self.mult - 1

        self.inp_fn = "xtb.xyz"
        self.out_fn = "xtb.out"
        self.to_keep = ("out:xtb.out", "grad", "xtbopt.xyz", "g98.out")

        self.parser_funcs = {
            "grad": self.parse_gradient,
            "hess": self.parse_hessian,
            "opt": self.parse_opt,
            "noparse": lambda path: None,
        }

        self.base_cmd = self.get_cmd("cmd")

    def reattach(self, last_calc_cycle):
        pass

    def prepare_coords(self, atoms, coords):
        coords = coords * BOHR2ANG
        return make_xyz_str(atoms, coords.reshape((-1, 3)))

    def prepare_input(self, atoms, coords, calc_type):
        return None

    def prepare_add_args(self):
        add_args = f"--gfn {self.gfn} --chrg {self.charge} --uhf {self.uhf} --acc {self.acc}".split()
        # Use solvent model if specified
        if self.gbsa:
            gbsa = f"--gbsa {self.gbsa}".split()
            add_args = add_args + gbsa
        return add_args

    def get_pal_env(self):
        env_copy = os.environ.copy()
        env_copy["OMP_NUM_THREADS"] = str(self.pal)
        env_copy["MKL_NUM_THREADS"] = str(self.pal)
        env_copy["OMP_STACKSIZE"] = "1000m"

        return env_copy

    def get_energy(self, atoms, coords):
        results = self.get_forces(atoms, coords)
        del results["forces"]
        return results

    def get_forces(self, atoms, coords):
        inp = self.prepare_coords(atoms, coords)
        add_args = self.prepare_add_args() + ["--grad"]
        self.log(f"Executing {self.base_cmd} {add_args}")
        kwargs = {
            "calc": "grad",
            "add_args": add_args,
            "env": self.get_pal_env(),
        }
        results = self.run(inp, **kwargs)
        return results

    def get_hessian(self, atoms, coords):
        inp = self.prepare_coords(atoms, coords)
        add_args = self.prepare_add_args() + ["--hess"]
        self.log(f"Executing {self.base_cmd} {add_args}")
        kwargs = {
            "calc": "hess",
            "add_args": add_args,
            "env": self.get_pal_env(),
        }
        results = self.run(inp, **kwargs)
        return results

    def run_calculation(self, atoms, coords):
        inp = self.prepare_coords(atoms, coords)
        kwargs = {
                "calc": "noparse",
                "env": self.get_pal_env(),
        }
        results = self.run(inp, **kwargs)
        return results

    def run_opt(self, atoms, coords, keep=True):
        inp = self.prepare_coords(atoms, coords)
        add_args = self.prepare_add_args() + ["--opt", "tight"]
        self.log(f"Executing {self.base_cmd} {add_args}")
        kwargs = {
            "calc": "opt",
            "add_args": add_args,
            "env": self.get_pal_env(),
            "keep": keep,
        }
        opt_geom = self.run(inp, **kwargs)
        return opt_geom

    def parse_opt(self, path):
        xtbopt = path / "xtbopt.xyz"
        if not xtbopt.exists():
            self.log(f"{self.calc_number:03d} failed")
            return None
        opt_geom = geom_from_xyz_file(xtbopt)
        opt_geom.energy = self.parse_energy(path)
        return opt_geom

    def parse_energy(self, path):
        with open(path / self.out_fn) as handle:
            text = handle.read()
        energy_re = "total E\s*:\s*([-\d\.]+)"
        energy = float(re.search(energy_re, text)[1])
        return energy

    def parse_gradient(self, path):
        return parse_turbo_gradient(path)

    def parse_hessian(self, path):
        with open(path / "hessian") as handle:
            text = handle.read()
        hessian = np.array(text.split()[1:], dtype=float)
        coord_num = int(hessian.size**0.5)
        hessian = hessian.reshape(coord_num, coord_num)
        energy = self.parse_energy(path)
        results = {
            "energy": energy,
            "hessian": hessian,
        }
        return results

    def __str__(self):
        return "XTB calculator"


if __name__ == "__main__":
    from pathlib import Path
    path = Path("/scratch/projekte/phosphor_fprakt/07_09_neb/02_irc/tmpo")
    xtb = XTB()
    xtb.parse_hessian(path)
