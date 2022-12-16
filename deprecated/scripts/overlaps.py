#!/usr/bin/env python3

"""Example YAML input:

geom:
 fn: lib:h2o2_hf_321g_opt.xyz
calc1:
 type: orca5
 keywords: hf sto-3g
 blocks: "%tddft nroots 2 iroot 1 end"
 pal: 2
calc2:
 type: orca5
 keywords: hf sto-3g
 blocks: "%tddft nroots 2 iroot 1 end"
 pal: 2
# Either wf|tden
ovlp_type: wf

"""

import argparse
from pprint import pprint
import sys
import time

import numpy as np
import yaml

from pysisyphus.calculators import ORCA, ORCA5, Gaussian16
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging


init_logging()
np.set_printoptions(suppress=True, precision=6)


def parse_args(args):

    parser = argparse.ArgumentParser()

    parser.add_argument("yaml")
    parser.add_argument("--ovlp-fn", dest="ovlp_fn", default="ovlp_mat.dat")
    parser.add_argument("--skip-calcs", dest="do_calc", action="store_false")
    parser.add_argument("--h5-fns", dest="h5_fns", nargs=2, default=None)
    parser.add_argument("--conf-thresh", dest="conf_thresh", default=0.001, type=float)

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    with open(args.yaml) as handle:
        run_dict = yaml.load(handle.read(), Loader=yaml.SafeLoader)
    pprint(run_dict)
    print()

    geom = geom_loader(run_dict["geom"]["fn"])

    CALCS = {"orca": ORCA, "gaussian16": Gaussian16, "orca5": ORCA5}

    ovlp_type = run_dict["ovlp_type"]

    def get_calc(key):
        calc_kwargs = run_dict[key]
        calc_kwargs["ovlp_type"] = ovlp_type
        calc_key = calc_kwargs.pop("type")
        dump_fn = f"overlap_data_{key}.h5"
        calc = CALCS[calc_key](**calc_kwargs, base_name=key, dump_fn=dump_fn)
        assert calc.root, "No 'root' set on calculator. Please specify an initial root."
        return calc

    calc1 = get_calc("calc1")
    calc2 = get_calc("calc2")

    calc_args = (geom.atoms, geom.coords)

    def calc_es(calc):
        print(f"Calculating ES for {calc} ... ", end="")
        start = time.time()
        calc.get_energy(*calc_args)
        dur = time.time() - start
        print(f"finished in {dur:.1f} s.")
        calc.store_overlap_data(*calc_args)
        calc.dump_overlap_data()

    if args.do_calc:
        calc_es(calc1)
        calc_es(calc2)
    else:
        try:
            h5_fn1, h5_fn2 = args.h5_fns
        except TypeError:
            h5_fn1 = calc1.dump_fn
            h5_fn2 = calc2.dump_fn
        print(f"Taking overlap_data from '{h5_fn1}' and '{h5_fn2}'.")
        calc1 = calc1.from_overlap_data(h5_fn1, set_wfow=True)
        calc2 = calc2.from_overlap_data(h5_fn2)

    conf_thresh = args.conf_thresh
    calc1.conf_thresh = conf_thresh
    if ovlp_type == "wf":
        calc1.wfow.conf_thresh = conf_thresh

    ao_ovlp = calc1.get_sao_from_mo_coeffs(calc1.mo_coeff_list[-1])
    print("Recreate S_AO from MO coeffs at calc1")
    ovlp_funcs = {
        "tden": "tden_overlap_with_calculator",
        "wf": "wf_overlap_with_calculator",
    }
    ovlp_func = ovlp_funcs[ovlp_type]
    print(f"Calculating {ovlp_type} overlaps")
    ovlp_mat = getattr(calc1, ovlp_func)(calc2, ao_ovlp=ao_ovlp)
    if ovlp_type == "wf":
        ovlp_mat = ovlp_mat[0]
    print("Rows along states of calc1, columns along states of calc2")
    print(ovlp_mat)
    ovlp_fn = f"{ovlp_type}_{args.ovlp_fn}"
    np.savetxt(ovlp_fn, ovlp_mat)
    print(f"Dumped overlap matrix to '{ovlp_fn}'.")


if __name__ == "__main__":
    run()
