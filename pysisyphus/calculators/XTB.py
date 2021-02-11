from collections import namedtuple
import io
import os
import re
import shutil
import textwrap

import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.calculators.parser import parse_turbo_gradient
from pysisyphus.constants import BOHR2ANG, BOHRPERFS2AU
from pysisyphus.helpers import geom_from_xyz_file, geoms_from_trj
from pysisyphus.xyzloader import make_xyz_str


OptResult = namedtuple("OptResult", "opt_geom opt_log")


class XTB(Calculator):

    conf_key = "xtb"

    def __init__(
        self,
        gbsa="",
        gfn=2,
        acc=1.0,
        topo=None,
        topo_update=None,
        mem=1000,
        quiet=False,
        **kwargs,
    ):
        super().__init__(**kwargs)

        self.gbsa = gbsa
        self.gfn = gfn
        self.acc = acc
        self.topo = topo
        self.topo_update = topo_update
        self.mem = mem
        self.quiet = quiet

        self.topo_used = 0
        valid_gfns = (0, 1, 2, "ff")
        assert self.gfn in valid_gfns, \
            f"Invalid gfn argument. Allowed arguments are: {', '.join(valid_gfns)}!"
        self.uhf = self.mult - 1

        self.inp_fn = "xtb.xyz"
        self.out_fn = "xtb.out"
        self.to_keep = (
            "out:xtb.out",
            "gradient",
            "xtbopt.xyz",
            "g98.out",
            "xtb.trj",
        )
        if self.quiet:
            self.to_keep = ()

        self.parser_funcs = {
            "grad": self.parse_gradient,
            "hess": self.parse_hessian,
            "opt": self.parse_opt,
            "md": self.parse_md,
            "topo": self.parse_topo,
            "noparse": lambda path: None,
        }

        self.base_cmd = self.get_cmd("cmd")

    def reattach(self, last_calc_cycle):
        pass

    def prepare_coords(self, atoms, coords):
        coords = coords * BOHR2ANG
        return make_xyz_str(atoms, coords.reshape((-1, 3)))

    def prepare_input(self, atoms, coords, calc_type):
        # Check if the topology has to be recreated/updated
        if (
            self.topo_used > 0
            and self.topo_update
            and (self.topo_used % self.topo_update == 0)
        ):
            results = self.run_topo(atoms, coords)
            self.topo = results["topo"]
            self.log(f"Updated topology! Saved to '{self.topo}'.")
        if self.topo:
            path = self.prepare_path(use_in_run=True)
            shutil.copy(self.topo, path / "gfnff_topo")
            self.log(f"Using toplogy given in {self.topo}.")
            self.topo_used += 1

    def prepare_add_args(self):
        add_args = f"--chrg {self.charge} --uhf {self.uhf} --acc {self.acc}".split()
        # Use solvent model if specified
        if self.gbsa:
            gbsa = f"--gbsa {self.gbsa}".split()
            add_args = add_args + gbsa
        # Select parametrization
        gfn = ["--gfnff"] if self.gfn == "ff" else f"--gfn {self.gfn}".split()
        add_args = add_args + gfn
        return add_args

    def get_pal_env(self):
        env_copy = os.environ.copy()
        env_copy["OMP_NUM_THREADS"] = str(self.pal)
        env_copy["MKL_NUM_THREADS"] = str(self.pal)
        # Per thread
        env_copy["OMP_STACKSIZE"] = f"{self.mem}M"

        return env_copy

    def get_energy(self, atoms, coords, prepare_kwargs=None):
        if prepare_kwargs is None:
            prepare_kwargs = {}  # lgtm [py/unused-local-variable]
        results = self.get_forces(atoms, coords)
        del results["forces"]
        return results

    def get_forces(self, atoms, coords, prepare_kwargs=None):
        if prepare_kwargs is None:
            prepare_kwargs = {}  # lgtm [py/unused-local-variable]
        self.prepare_input(atoms, coords, "forces")
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

    def get_hessian(self, atoms, coords, prepare_kwargs=None):
        if prepare_kwargs is None:
            prepare_kwargs = {}  # lgtm [py/unused-local-variable]
        self.prepare_input(atoms, coords, "hessian")
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
        self.prepare_input(atoms, coords, "calculation")
        inp = self.prepare_coords(atoms, coords)
        kwargs = {
            "calc": "noparse",
            "env": self.get_pal_env(),
        }
        results = self.run(inp, **kwargs)
        return results

    def run_topo(self, atoms, coords):
        inp = self.prepare_coords(atoms, coords)
        kwargs = {
            "calc": "topo",
            "cmd": [self.base_cmd, "topo"],
            "env": self.get_pal_env(),
        }
        results = self.run(inp, **kwargs)
        return results

    def parse_topo(self, path):
        fn = "gfnff_topo"
        topo = path / fn
        target = self.make_fn(fn)
        shutil.copy(topo, target)

        return {
            "topo": target,
        }

    def get_mdrestart_str(self, coords, velocities):
        """coords and velocities have to given in au!"""
        vals = np.concatenate((coords, velocities), axis=1)

        with io.StringIO() as io_stream:
            np.savetxt(io_stream, vals, fmt="% .14e")
            mdrestart = io_stream.getvalue()
        # What does the -1.0 mean?
        mdrestart = "-1.0\n" + mdrestart.replace("e", "D")
        mdrestart = textwrap.indent(mdrestart, " ")
        return mdrestart

    def write_mdrestart(self, path, mdrestart_str):
        with open(path / "mdrestart", "wb") as handle:
            handle.write(mdrestart_str.encode("ascii"))

    def run_md(self, atoms, coords, t, dt, velocities=None, dump=1):
        """Expecting t and dt in fs, even though xtb wants t in ps!"""

        restart = "false"
        path = self.prepare_path(use_in_run=True)
        if velocities is not None:
            coords3d = coords.reshape(-1, 3)
            velocities3d = velocities.reshape(-1, 3) * BOHRPERFS2AU
            assert (
                coords3d.shape == velocities3d.shape
            ), "Shape of coordinates and velocities doesn't match!"
            mdrestart_str = self.get_mdrestart_str(coords3d, velocities3d)
            self.write_mdrestart(path, mdrestart_str)
            restart = "true"
        md_str = textwrap.dedent(
            """
        $md
            hmass=1
            dump={dump}  # fs
            nvt=false
            restart={restart}
            time={time}  # ps
            shake=0
            step={step}  # fs
            velo=false
        $end"""
        )
        t_ps = t / 1000
        md_str_fmt = md_str.format(restart=restart, time=t_ps, step=dt, dump=dump)
        with open(path / "xcontrol", "w") as handle:
            handle.write(md_str_fmt)
        inp = self.prepare_coords(atoms, coords)

        add_args = self.prepare_add_args() + ["--input", "xcontrol", "--md"]
        self.log(f"Executing {self.base_cmd} {add_args}")
        kwargs = {
            "calc": "md",
            "add_args": add_args,
            "env": self.get_pal_env(),
            "keep": True,
        }
        geoms = self.run(inp, **kwargs)
        return geoms

    def parse_md(self, path):
        assert (path / "xtbmdok").exists(), "File xtbmdok does not exist!"
        geoms = geoms_from_trj(path / "xtb.trj")
        return geoms

    def run_opt(self, atoms, coords, keep=True, keep_log=False):
        inp = self.prepare_coords(atoms, coords)
        add_args = self.prepare_add_args() + ["--opt", "tight"]
        self.log(f"Executing {self.base_cmd} {add_args}")
        kwargs = {
            "calc": "opt",
            "add_args": add_args,
            "env": self.get_pal_env(),
            "keep": keep,
            "parser_kwargs": {"keep_log": keep_log},
        }
        opt_result = self.run(inp, **kwargs)
        return opt_result

    def parse_opt(self, path, keep_log=False):
        xtbopt = path / "xtbopt.xyz"
        if not xtbopt.exists():
            self.log(f"{self.calc_number:03d} failed")
            return None
        opt_geom = geom_from_xyz_file(xtbopt)
        opt_geom.energy = self.parse_energy(path)

        opt_log = None
        if keep_log:
            opt_log = geoms_from_trj(path / "xtbopt.log")

        opt_result = OptResult(opt_geom=opt_geom, opt_log=opt_log)
        return opt_result

    def parse_energy(self, path):
        with open(path / self.out_fn) as handle:
            text = handle.read()
        energy_re = "TOTAL ENERGY\s*([-\d\.]+) Eh"
        energy = float(re.search(energy_re, text)[1])
        return energy

    def parse_gradient(self, path):
        return parse_turbo_gradient(path)

    def parse_hessian(self, path):
        with open(path / "hessian") as handle:
            text = handle.read()
        hessian = np.array(text.split()[1:], dtype=float)
        coord_num = int(hessian.size ** 0.5)
        hessian = hessian.reshape(coord_num, coord_num)
        energy = self.parse_energy(path)
        results = {
            "energy": energy,
            "hessian": hessian,
        }
        return results

    def __str__(self):
        return "XTB calculator"
