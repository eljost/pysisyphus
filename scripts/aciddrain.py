#!/usr/bin/env python

import argparse
import os
from pathlib import Path
import shutil
import sys
import tempfile

import luigi
import psutil
import yaml
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

from pysisyphus.calculators import ORCA5, XTB
from pysisyphus.drivers.pka import direct_cycle, G_aq_from_h5_hessian
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.helpers import geom_loader, do_final_hessian
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


class Params(luigi.Config):
    id_ = luigi.IntParameter()
    name = luigi.Parameter()
    h_ind = luigi.IntParameter()
    is_base = luigi.BoolParameter(default=False)
    charge = luigi.IntParameter()
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
            return Minimization(
                self.id_, self.name, self.h_ind, is_base=False, charge=self.charge
            )
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


class PreMinimization(Params, luigi.Task):
    def output(self):
        return luigi.LocalTarget(self.get_path("preopt.xyz"))

    def requires(self):
        return InputGeometry(self.id_, self.name, self.h_ind, self.is_base, self.charge)

    def run(self):
        print(highlight_text(f"Pre-Minimization {self.key}"))
        geom = geom_loader(self.input().path, coord_type="redund")
        geom.set_calculator(
            get_xtb_calc(charge=self.calc_charge(), out_dir=self.qm_out_dir)
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            opt = RFOptimizer(
                geom,
                dump=True,
                overachieve_factor=2.0,
                thresh="gau_loose",
                out_dir=tmp_dir,
                max_cycles=500,
            )
            opt.run()
            assert opt.is_converged
            with self.output().open("w") as handle:
                handle.write(geom.as_xyz())


class Minimization(Params, luigi.Task):
    def output(self):
        return (
            luigi.LocalTarget(self.get_path("opt.xyz")),
            luigi.LocalTarget(self.get_path("opt_hessian.h5")),
        )

    def requires(self):
        # Derive the initial base geometry from the optimized acid geometry.
        # Maybe some cycles could be saved when the base is also pre-
        # optimized, but things could also go wrong.
        if self.is_base:
            return InputGeometry(
                self.id_, self.name, self.h_ind, self.is_base, self.charge
            )
        # Only preoptimize initial acid geometry, not the base.
        else:
            return PreMinimization(
                self.id_, self.name, self.h_ind, self.is_base, self.charge
            )

    def run(self):
        geom = geom_loader(self.input().path, coord_type="redund")
        geom.set_calculator(
            get_calc(charge=self.calc_charge(), out_dir=self.qm_out_dir)
        )
        final_hess_fn = "final_hessian.h5"

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir = Path(tmp_dir)
            opt_kwargs_ = {
                "dump": True,
                "thresh": "gau",
                "out_dir": tmp_dir,
                "max_cycles": 250,
            }
            # Iterate until no imaginary frequencies are present
            for i in range(5):
                base_infix = " base" if self.is_base else ""
                print(highlight_text(f"Minimization {self.key}{base_infix}, cycle {i}"))
                opt_kwargs = opt_kwargs_.copy()
                if i > 0:
                    opt_kwargs.update(
                        {
                            "hessian_init": tmp_dir / final_hess_fn,
                            # Disable overachieve factor to avoid immediate convergence
                            # in the first cycle
                            "overachieve_factor": 0.0,
                        }
                    )

                opt = RFOptimizer(
                    geom,
                    **opt_kwargs,
                )
                opt.run()
                assert opt.is_converged
                #
                xyz_out, hess_out = self.output()
                with xyz_out.open("w") as handle:
                    handle.write(geom.as_xyz())
                hess_result = do_final_hessian(geom, out_dir=tmp_dir)
                if len(hess_result.neg_eigvals) == 0:
                    self.backup_from_dir(tmp_dir, final_hess_fn, hess_out.path)
                    break
                else:
                    suffix = f".{i:02d}"
                    self.backup_from_dir(tmp_dir, final_hess_fn, hess_out.path + suffix)
                    self.backup_from_dir(
                        tmp_dir, "final_geometry.xyz", xyz_out.path + suffix
                    )
                    # Delete everything in tmp_dir, besides final_hess_fn, as it will be
                    # reused.
                    files = [
                        f
                        for f in tmp_dir.glob("./*")
                        if f.is_file() and f.name != final_hess_fn
                    ]
                    for f in files:
                        os.remove(f)
            else:
                raise Exception("Minimization failed!")


class SolvEnergy(Params, luigi.Task):
    def output(self):
        return luigi.LocalTarget(self.get_path("solv_energy"))

    def requires(self):
        return Minimization(self.id_, self.name, self.h_ind, self.is_base, self.charge)

    def run(self):
        geom = geom_loader(self.input()[0].path)
        geom.set_calculator(
            get_solv_calc(charge=self.calc_charge(), out_dir=self.qm_out_dir)
        )
        solv_energy = geom.energy
        with self.output().open("w") as handle:
            handle.write(str(solv_energy))


class DirectCycle(Params, luigi.Task):
    pka_exp = luigi.FloatParameter()

    def output(self):
        return luigi.LocalTarget(self.get_path("summary.yaml"))

    def requires(self):
        return (
            Minimization(
                self.id_, self.name, self.h_ind, is_base=False, charge=self.charge
            ),
            Minimization(
                self.id_, self.name, self.h_ind, is_base=True, charge=self.charge
            ),
            SolvEnergy(
                self.id_, self.name, self.h_ind, is_base=False, charge=self.charge
            ),
            SolvEnergy(
                self.id_, self.name, self.h_ind, is_base=True, charge=self.charge
            ),
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
        G_aq_H = -6.28 + (-265.9)  # corresponds to the default
        pka = direct_cycle(acid_h5, base_h5, acid_solv_en, base_solv_en, G_aq_H=G_aq_H)
        print(f"@@@ {self.name}: pKa={pka:.4f}")

        G_acid_aq = G_aq_from_h5_hessian(acid_h5, acid_solv_en)
        G_base_aq = G_aq_from_h5_hessian(base_h5, base_solv_en)
        G_diss_aq = G_base_aq - G_acid_aq

        results = {
            "name": self.name,
            "h_ind": self.h_ind,
            "G_acid_aq": G_acid_aq,
            "G_base_aq": G_base_aq,
            "G_diss_aq": G_diss_aq,
            "pka_calc": pka,
            "pka_exp": self.pka_exp,
        }
        with self.output().open("w") as handle:
            yaml.dump(results, handle)


class DirectCycler(luigi.Task):
    yaml_inp = luigi.Parameter()

    def requires(self):
        with open(self.yaml_inp) as handle:
            run_dict = yaml.load(handle.read(), Loader=yaml.SafeLoader)

        for id_, (acid, acid_dict) in enumerate(run_dict["acids"].items()):
            fn = acid_dict["fn"]
            h_ind = acid_dict["h_ind"]
            charge = acid_dict.get("charge", 0)
            name = Path(fn).stem
            pka_exp = acid_dict["pka_exp"]
            yield DirectCycle(
                id_=id_, name=name, h_ind=h_ind, charge=charge, pka_exp=pka_exp
            )

    def output(self):
        return luigi.LocalTarget("summary.yaml")

    def run(self):
        summary = {}
        for dc in self.input():
            with dc.open() as handle:
                results = yaml.load(handle, Loader=yaml.SafeLoader)

            name = results["name"]
            pka_calc = results["pka_calc"]
            pka_exp = results["pka_exp"]
            h_ind = results["h_ind"]
            summary[name] = {
                "h_ind": h_ind,
                "pka_exp": pka_exp,
                "pka_calc": pka_calc,
            }
            pka_diff = pka_calc - pka_exp
            print(
                f"@@@\t{dc.path} pkas: calc={pka_calc:.2f}, exp={pka_exp:.2}, Δ={pka_diff:.2f}"
            )

        with self.output().open("w") as handle:
            yaml.dump(summary, handle)


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("yaml")
    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument("--lfer", action="store_true")
    action_group.add_argument("--pka", action="store_true")
    return parser.parse_args(args)


def run_pka(yaml_inp):
    with open(yaml_inp) as handle:
        run_dict = yaml.load(handle.read(), Loader=yaml.SafeLoader)

    inputs = list()
    for acid, acid_dict in run_dict["acids"].items():
        fn = acid_dict["fn"]
        h_ind = acid_dict["h_ind"]
        assert Path(fn).exists(), f"File '{fn}' does not exist!"
        geom = geom_loader(fn)
        assert (
            geom.atoms[h_ind].lower() == "h"
        ), f"Atom at index {h_ind} in '{fn}' is not a hydrogen atom!"
        print(f"Checked {acid}.")
        inputs.append((fn, h_ind))

    global geom_queue
    geom_queue = [geom_loader(fn) for fn, _ in inputs]

    # Calculator setup
    pal = psutil.cpu_count(logical=False)
    calc_cls = ORCA5
    rdc = run_dict["calc"]
    gas_calc_cls = rdc.pop("type")
    gas_kwargs = rdc
    rdsc = run_dict["solv_calc"]
    solv_calc_cls = rdsc.pop("type")
    solv_kwargs = rdsc
    assert gas_calc_cls == solv_calc_cls
    calc_dict = {
        "orca5": ORCA5,
        "xtb": XTB,
    }
    calc_cls = calc_dict[gas_calc_cls]

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

    global get_xtb_calc

    def get_xtb_calc(charge, out_dir):
        return XTB(pal=pal, charge=charge, out_dir=out_dir, base_name="xtb")

    luigi.build(
        # (Minimization(id_=0, name="formicacid", h_ind=4, is_base=False),)
        # (Minimization(id_=0, name="formicacid", h_ind=4, is_base=True), )
        (DirectCycler(yaml_inp),),
        local_scheduler=True,
    )


def run_lfer(yaml_inp):
    with open(yaml_inp) as handle:
        run_dict = yaml.load(handle.read(), Loader=yaml.SafeLoader)
    names = list()
    exps = list()
    calcs = list()
    for name, values in run_dict.items():
        names.append(name)
        exps.append(values["pka_exp"])
        calcs.append(values["pka_calc"])
    min_ = min(min(exps), min(calcs)) - 1
    max_ = max(max(exps), max(calcs)) + 1
    lims = (min_, max_)

    _, (ax0, ax1) = plt.subplots(ncols=2)
    ax0.scatter(exps, calcs, s=30)
    ax0.set_xlabel("pKa experimental")
    ax0.set_ylabel("pKa calcualted")
    ax0.set_title("experimental vs. calculated")

    # Linear regression
    res = linregress(calcs, exps)
    print("pka_calc,corr = m * pka_calc + n")
    print(f"m={res.slope}")
    print(f"n={res.intercept}")
    print(f"f(x)={res.slope:.6f}*x + {res.intercept:.6f}")

    def linfit(xs):
        return res.slope * xs + res.intercept

    fmt = " >6.2f"
    for name, exp_, calc in zip(names, exps, calcs):
        calc_corr = linfit(calc)
        print(f"{name: >32s}: exp={exp_:{fmt}} calc={calc:{fmt}} corr={calc_corr:{fmt}}")

    calc_corr = linfit(np.array(calcs))
    ax1.scatter(exps, calc_corr, s=30)
    ax1.set_xlabel("pKa calculated, LFER corrected")
    ax1.set_ylabel("pKa experimental")
    ax1.set_title("LFER")

    for ax in (ax0, ax1):
        ax.plot(lims, lims, c="k", ls="--")
        ax.set_xlim(lims)
        ax.set_ylim(lims)

    plt.tight_layout()
    plt.show()


def run():
    args = parse_args(sys.argv[1:])

    if args.pka:
        run_pka(args.yaml)
    elif args.lfer:
        run_lfer(args.yaml)


if __name__ == "__main__":
    run()
