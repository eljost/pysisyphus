#!/usr/bin/env python

"""
    PYTHONPATH="." luigi-deps-tree --module test_luigi ReactionPath --id- 1
    PYTHONPATH="." luigi --module test_luigi ReactionPath --local-scheduler --id- 1
    python test_luigi.py inp/input.xyz 53 9 -1 11 54 -1 9 11 1
"""

import argparse
from collections import deque
import os
from pathlib import Path
import shutil
import sys
import tempfile

import luigi
import numpy as np

from pysisyphus.calculators import XTB
from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.helpers import geom_loader, do_final_hessian
from pysisyphus.intcoords.helpers import get_weighted_bond_mode_getter
from pysisyphus.cos.GrowingNT import GrowingNT
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.irc import EulerPC


class Params(luigi.Config):
    id_ = luigi.IntParameter()
    base = luigi.Parameter(default="out")
    step = luigi.IntParameter(default=0)

    @property
    def key(self):
        return f"{self.id_:04d}_{self.step:02d}"

    @property
    def out_dir(self):
        return Path(f"{self.base}_{self.id_:04d}")

    @property
    def step_str(self):
        return f"{self.step:02d}"

    def get_path(self, fn):
        out_dir = self.out_dir
        if not out_dir.exists():
            os.mkdir(out_dir)
        return out_dir / f"{self.step_str}_{fn}"

    def get_prefix(self, prefix):
        return f"{self.step_str}_{prefix}"

    def get_calc(self):
        return XTB(
            pal=6, quiet=True, retry_etemp=1000, gbsa="methanol", out_dir=self.get_path("qm_calcs/")
        )

    def backup_from_dir(self, dir_, fn, dest_fn=None):
        if dest_fn is None:
            dest_fn = self.get_path(fn)
        shutil.copy(Path(dir_) / fn, dest_fn)


class Minimum(Params, luigi.Task):
    def output(self):
        return luigi.LocalTarget(self.get_path("input.xyz"))

    def run(self):
        geom = geoms[self.id_]

        with self.output().open("w") as handle:
            handle.write(geom.as_xyz())


class GNT(Params, luigi.Task):
    def output(self):
        return luigi.LocalTarget(self.get_path("gnt_ts.xyz"))

    def requires(self):
        return Minimum(self.id_, self.base, self.step)

    def run(self):
        geom = geom_loader(self.input().path)
        geom.set_calculator(self.get_calc())
        bonds = get_weighted_bond_mode(geom.atoms, geom.coords3d)

        print(f"@@@ {self.key}: Using bonds: {bonds}")

        gnt_kwargs = {
            "step_len": 0.1,
            "bonds": [bonds[self.id_]],
            "stop_after_ts": True,
            "rms_thresh": 0.003,
        }
        gnt = GrowingNT(geom, **gnt_kwargs)
        with tempfile.TemporaryDirectory() as tmp_dir:
            opt_kwargs = {
                "max_cycles": 1_000,
                "dump": True,
                "out_dir": tmp_dir,
                "prefix": "gnt",
            }
            opt = PreconLBFGS(gnt, **opt_kwargs)
            opt.run()
            self.backup_from_dir(tmp_dir, "gnt_optimization.trj")
        with self.output().open("w") as handle:
            handle.write(gnt.ts_images[0].as_xyz())


class TSOpt(Params, luigi.Task):
    def output(self):
        return (
            luigi.LocalTarget(self.get_path("ts_opt.xyz")),
            luigi.LocalTarget(self.get_path("ts_hessian.h5")),
        )

    def requires(self):
        return GNT(self.id_, self.base, self.step)

    def run(self):
        geom = geom_loader(self.input().path, coord_type="redund")
        geom.set_calculator(self.get_calc())

        with tempfile.TemporaryDirectory() as tmp_dir:
            tsopt = RSPRFOptimizer(
                # tsopt = RSIRFOptimizer(
                geom,
                max_cycles=75,
                out_dir=tmp_dir,
                prefix="ts",
                dump=True,
                hessian_recalc=5,
                trust_max=0.3,
                overachieve_factor=3.0,
                # max_line_search=True,
                # min_line_search=True,
            )
            tsopt.run()
            self.backup_from_dir(tmp_dir, "ts_optimization.trj")
            assert tsopt.is_converged

            xyz_out, hess_out = self.output()
            with xyz_out.open("w") as handle:
                handle.write(geom.as_xyz())
            do_final_hessian(geom, out_dir=tmp_dir)
            self.backup_from_dir(tmp_dir, "final_hessian.h5", hess_out.path)


class IRC(Params, luigi.Task):
    def output(self):
        return (
            luigi.LocalTarget(self.get_path("irc_first.xyz")),
            luigi.LocalTarget(self.get_path("irc_last.xyz")),
        )

    def requires(self):
        return TSOpt(self.id_, self.base, self.step)

    def run(self):
        ts_xyz, ts_hess = self.input()
        geom = geom_loader(ts_xyz.path)
        geom.set_calculator(self.get_calc())

        with tempfile.TemporaryDirectory() as tmp_dir:
            irc = EulerPC(
                geom,
                out_dir=tmp_dir,
                hessian_init=ts_hess.path,
                hessian_recalc=10,
                rms_grad_thresh=0.0005,
            )
            irc.run()
            self.backup_from_dir(tmp_dir, "finished_irc.trj")
        first_geom = geom.copy()
        first_geom.coords = irc.all_coords[0]
        last_geom = geom.copy()
        last_geom.coords = irc.all_coords[-1]
        for geom_, target in zip((first_geom, last_geom), self.output()):
            with target.open("w") as handle:
                handle.write(geom_.as_xyz())


class EndOpt(Params, luigi.Task):
    pref_map = {
        "first": 0,
        "last": -1,
    }
    prefix = luigi.Parameter()

    @property
    def pref_(self):
        return f"{self.prefix}_endopt"

    def output(self):
        return luigi.LocalTarget(self.get_path(f"{self.pref_}.xyz"))

    def requires(self):
        return IRC(self.id_, self.base, self.step)

    def run(self):
        ind = self.pref_map[self.prefix]
        geom = geom_loader(self.input()[ind].path, coord_type="tric")
        geom.set_calculator(self.get_calc())
        with tempfile.TemporaryDirectory() as tmp_dir:
            opt = RFOptimizer(
                geom,
                dump=True,
                overachieve_factor=3.0,
                out_dir=tmp_dir,
            )
            opt.run()
            with self.output().open("w") as handle:
                handle.write(geom.as_xyz())


class FirstEndOpt(EndOpt):
    prefix = "first"


class LastEndOpt(EndOpt):
    prefix = "last"


class ReactionPath(Params, luigi.Task):
    def requires(self):
        return (
            Minimum(self.id_, self.base, self.step),
            FirstEndOpt(self.id_, self.base, self.step),
            TSOpt(self.id_, self.base, self.step),
            LastEndOpt(self.id_, self.base, self.step),
        )

    def run(self):
        if self.step >= 1:
            self.complete = lambda: True

        minimum, first_endopt, tsopt, last_endopt = self.input()
        ref_geom, first_geom, last_geom = [
            geom_loader(target.path) for target in (minimum, first_endopt, last_endopt)
        ]
        first_rmsd = ref_geom.rmsd(first_geom)
        last_rmsd = ref_geom.rmsd(last_geom)
        print("first", first_rmsd, "last", last_rmsd)

        calc = self.get_calc()
        atoms = ref_geom.atoms

        def get_energy(geom):
            return calc.get_energy(atoms, geom.cart_coords)["energy"]

        ts_geom = geom_loader(tsopt[0].path)
        ts_energy = get_energy(ts_geom)
        first_energy = get_energy(first_geom)
        last_energy = get_energy(last_geom)
        energies = np.array((first_energy, ts_energy, last_energy))
        energies -= energies.min()
        energies *= AU2KJPERMOL
        comments = [f"{en:.2f} kJ mol⁻¹" for en in energies]
        first_kj, ts_kj, last_kj = energies
        first_ts = ts_kj - first_kj
        last_ts = ts_kj - last_kj
        print(f"@@@ {self.key}: TS - first: {first_ts:.2f} kJ mol⁻¹")
        print(f"@@@ {self.key}: TS - last: {last_ts:.2f} kJ mol⁻¹")

        with open(self.get_path("first_ts_last.trj"), "w") as handle:
            handle.write(
                "\n".join(
                    [
                        geom.as_xyz(comment=comment)
                        for geom, comment in zip(
                            (first_geom, ts_geom, last_geom), comments
                        )
                    ]
                )
            )

        # pick geom with higher rmsd

        next_ind = 1 if last_rmsd > first_rmsd else 0
        next_key = ("first", "last")[next_ind]
        next_geom = (first_geom, last_geom)[next_ind]
        # check if bond constraints are satisfied at next_geom
        bonds = get_weighted_bond_mode(next_geom.atoms, next_geom.coords3d)
        print(f"@@@ {self.key}: Trying to continue with {next_key}_geom")
        print(f"@@@ {self.key}: bonds for next step {bonds}\n@@@")
        if bonds:
            new_step = self.step + 1
            if new_step == 1:
                with open(self.out_dir / f"{new_step:02d}_input.xyz", "w") as handle:
                    handle.write(next_geom.as_xyz())
                yield ReactionPath(self.id_, self.base, new_step)
        else:
            self.complete = lambda: True


class ReactionPaths(luigi.WrapperTask):
    def requires(self):
        for id_, _ in enumerate(geom):
            yield ReactionPath(id_=id_)


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("fn")
    parser.add_argument("bonds", nargs="+", type=int)
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    fn = args.fn
    bonds = args.bonds

    global geoms
    geoms = list(geom_loader(fn, iterable=True))
    global bonds
    bonds = np.load(bonds)

    luigi.build((ReactionPaths(),), local_scheduler=True)


if __name__ == "__main__":
    run()
