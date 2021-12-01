#!/usr/bin/env python

import argparse
import os
from pathlib import Path
import shutil
import sys
import tempfile

import luigi
from luigi.tools.deps_tree import print_tree
import psutil
import yaml

from pysisyphus.calculators import ORCA5, XTB
from pysisyphus.drivers.pka import direct_cycle
from pysisyphus.helpers import geom_loader, do_final_hessian
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


class Params(luigi.Config):
    id_ = luigi.IntParameter()
    name = luigi.Parameter()
    h_ind = luigi.IntParameter()
    is_base = luigi.BoolParameter(default=False)
    charge = luigi.IntParameter(default=0)
    base = luigi.Parameter(default="out")

    @property
    def key(self):
        return f"{self.id_:04d}_{self.name}"

    @property
    def out_dir(self):
        return Path(f"{self.base}_{self.key}")

    def get_path(self, fn):
        out_dir = self.out_dir
        is_base_str = "_base" if self.is_base else ""
        fn = f"{self.name}{is_base_str}_{fn}"
        if not out_dir.exists():
            os.mkdir(out_dir)
        return out_dir / fn

    def backup_from_dir(self, dir_, fn, dest_fn=None):
        if dest_fn is None:
            dest_fn = self.get_path(fn)
        shutil.copy(Path(dir_) / fn, dest_fn)

    def calc_charge(self):
        charge = self.charge + (-1 if self.is_base else 0)
        return charge

    @property
    def qm_out_dir(self):
        qm_out_dir = self.get_path("qm_calcs")
        return qm_out_dir


class InputGeometry(Params, luigi.Task):
    def output(self):
        return luigi.LocalTarget(self.get_path("input.xyz"))

    def requires(self):
        if self.is_base:
            return Minimization(self.id_, self.name, self.h_ind, is_base=False)
        else:
            return None

    def run(self):
        # Derive initial geometry of the base from the optimized acid
        if self.is_base:
            acid_geom = geom_loader(self.input()[0].path)
            # Be sure that it is actually an H atom.
            assert acid_geom.atoms[self.h_ind].lower() == "h"
            geom = acid_geom.get_subgeom_without((self.h_ind,))
        else:
            geom = geom_queue[self.id_]

        with self.output().open("w") as handle:
            handle.write(geom.as_xyz())


class Minimization(Params, luigi.Task):
    def output(self):
        return (
            luigi.LocalTarget(self.get_path("opt.xyz")),
            luigi.LocalTarget(self.get_path("opt_hessian.h5")),
        )

    def requires(self):
        return InputGeometry(self.id_, self.name, self.h_ind, self.is_base)

    def run(self):
        geom = geom_loader(self.input().path, coord_type="redund")
        geom.set_calculator(
            get_calc(charge=self.calc_charge(), out_dir=self.qm_out_dir)
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            opt = RFOptimizer(
                geom,
                dump=True,
                overachieve_factor=4.0,
                out_dir=tmp_dir,
                thresh="gau",
            )
            opt.run()
            #
            xyz_out, hess_out = self.output()
            with xyz_out.open("w") as handle:
                handle.write(geom.as_xyz())
            do_final_hessian(geom, out_dir=tmp_dir)
            self.backup_from_dir(tmp_dir, "final_hessian.h5", hess_out.path)


class SolvEnergy(Params, luigi.Task):
    def output(self):
        return luigi.LocalTarget(self.get_path("solv_energy"))

    def requires(self):
        return Minimization(self.id_, self.name, self.h_ind, self.is_base)

    def run(self):
        geom = geom_loader(self.input()[0].path)
        geom.set_calculator(
            get_solv_calc(charge=self.calc_charge(), out_dir=self.qm_out_dir)
        )
        solv_energy = geom.energy
        with self.output().open("w") as handle:
            handle.write(str(solv_energy))


class DirectCycle(Params, luigi.Task):
    def output(self):
        return luigi.LocalTarget(self.get_path("pka"))

    def requires(self):
        return (
            Minimization(self.id_, self.name, self.h_ind, is_base=False),
            Minimization(self.id_, self.name, self.h_ind, is_base=True),
            SolvEnergy(self.id_, self.name, self.h_ind, is_base=False),
            SolvEnergy(self.id_, self.name, self.h_ind, is_base=True),
        )

    def run(self):
        acid_inp, base_inp, acid_solv_fn, base_solv_fn = self.input()
        _, acid_h5 = [inp.path for inp in acid_inp]
        _, base_h5 = [inp.path for inp in base_inp]

        # Load solvated electronic energy from files
        with open(acid_solv_fn.path) as handle:
            acid_solv_en = float(handle.read())
        with open(base_solv_fn.path) as handle:
            base_solv_en = float(handle.read())
        pKa = direct_cycle(acid_h5, base_h5, acid_solv_en, base_solv_en)
        print(f"@@@ pka={pKa:.4f}")
        with self.output().open("w") as handle:
            handle.write(str(pKa))


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("yaml")
    return parser.parse_args(args) 

def run():
    args = parse_args(sys.argv[1:])

    with open(args.yaml) as handle:
        run_dict = yaml.load(handle.read(), Loader=yaml.SafeLoader)

    inputs = list()
    for acid, acid_dict in run_dict["acids"].items():
        fn = acid_dict["fn"]
        h_ind = acid_dict["h_ind"]
        assert Path(fn).exists(), \
            f"File '{fn}' does not exist!"
        geom = geom_loader(fn)
        assert geom.atoms[h_ind].lower() == "h", \
            f"Atom at index {h_ind} in '{fn}' is not a hydrogen atom!"
        print(f"Checked {acid}.")
        inputs.append((fn, h_ind))

    global geom_queue
    geom_queue = [geom_loader(fn) for fn, _ in inputs]

    pal = psutil.cpu_count(logical=False)
    calc_cls = ORCA5
    rdc = run_dict["calc"]
    gas_calc_cls = rdc.pop("type")
    gas_kwargs = rdc#run_dict["calc"]
    rdsc  = run_dict["solv_calc"]
    solv_calc_cls = rdsc.pop("type")
    solv_kwargs = rdsc
    assert gas_calc_cls == solv_calc_cls
    calc_dict = {
        "orca5": ORCA5,
        "xtb": XTB,
    }
    calc_cls = calc_dict[gas_calc_cls]

    cycles = list()
    for id_, (fn, h_ind) in enumerate(inputs):
        name = Path(fn).stem
        cycles.append(DirectCycle(id_=id_, name=name, h_ind=h_ind))

    global get_calc

    def get_calc(charge, out_dir):
        return calc_cls(
            pal=pal, charge=charge, out_dir=out_dir, base_name="gas", **gas_kwargs
        )

    global get_solv_calc

    def get_solv_calc(charge, out_dir):
        return calc_cls(
            pal=pal,
            charge=charge,
            out_dir=out_dir,
            base_name="solv",
            **solv_kwargs,
        )

    print(print_tree(cycles[0], indent="", last=True))

    luigi.build(
        # (Minimization(id_=0, name="formicacid", h_ind=4, is_base=False),)
        # (Minimization(id_=0, name="formicacid", h_ind=4, is_base=True), )
        cycles,
        local_scheduler=True,
    )


if __name__ == "__main__":
    run()
