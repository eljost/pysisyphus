#!/usr/bin/env python3

import itertools as it
import logging
import os
from pathlib import Path
import re

from natsort import natsorted
import numpy as np

from pysisyphus.calculators.Gaussian09 import Gaussian09
from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.calculators.Turbomole import Turbomole
from pysisyphus.helpers import (geom_from_xyz_file, index_array_from_overlaps,
                                np_print
)


class Overlapper:
    orca_exts = ("out", "gbw", "cis")
    gaussian_exts = ("log", "fchk", "dump_635r")

    logger = logging.getLogger("overlapper")

    def __init__(self, path, ovlp_with="previous", prev_n=0,
                 calc_key=None, calc_kwargs=None):
        self.path = Path(path)
        self.ovlp_with = ovlp_with
        assert ovlp_with in ("previous", "first")
        self.calc_key = calc_key
        self.calc_kwargs = calc_kwargs
        self.calc_kwargs["out_dir"] = path

        mobj = re.match("previous(\d+)", self.ovlp_with)
        self.prev_n = prev_n
        assert self.prev_n >= 0

        self.setter_dict = {
            "gaussian09": self.set_g16_files,
            "gaussian16": self.set_g16_files,
            "orca": self.set_orca_files,
            "turbomole": self.set_turbo_files,
        }
        self.files_from_dir_dict = {
            "orca": self.set_orca_files_from_dir,
            "gaussian09": self.set_gaussian16_files_from_dir,
            "gaussian16": self.set_gaussian16_files_from_dir,
        }

    def log(self, message, lvl="info"):
        func = getattr(self.logger, lvl)
        func(message)

    def keyfunc(self, element):
        regex = "_(\d+)\.(\d+)\."
        mobj = re.search(regex, element)
        return tuple([int(num) for num in mobj.groups()])

    def discover_files(self, path):
        image_str = "image_"
        calc_str = "calculator_"
        files = [str(f) for f in path.glob(image_str + "*")]
        if len(files) > 1:
            base_str = image_str
        else:
            files = [str(f) for f in path.glob(calc_str + "*")]
            base_str = calc_str
        print(f"Found {len(files)} files starting with '{base_str}'. "
              f"I assume that base string is '{base_str}'.")
        files = sorted(files, key=self.keyfunc)
        files_dict = dict()
        for key, elements in it.groupby(files, self.keyfunc):
            files_dict[key] = list(elements)
        return files_dict

    def discover_files_in_dir(self, path, exts):
        files_list = list()
        for ext in exts:
            glob = f"*{ext}"
            fns = [_ for _ in path.glob(glob)]
            assert len(fns) == 1, f"Searched for *.{ext} and was expecting " \
                                  f"one file but found {len(fns)} files " \
                                  f"instead: {fns}"
            fn = str(fns[0])
            files_list.append(fn)
        return files_list

    def discover_geometries(self, path):
        xyz_fns = natsorted(path.glob("*.xyz"))
        geoms = [geom_from_xyz_file(xyz_fn) for xyz_fn in xyz_fns]
        self.restore_calculators(geoms)

        return geoms

    def file_by_ext(self, iterable, ext):
        matches = [f for f in iterable if f.endswith(ext)]
        if len(matches) == 0:
            raise Exception(f"Couldn't file with extension '{ext}'!")
        assert len(matches) == 1
        return matches[0]

    def set_files_on_calculator(self, geom, files_dict, calc_class, exts,
                                calc_number, cycle_number=0):
        key = (calc_number, cycle_number)
        files = files_dict[key]
        calc = calc_class(calc_number=calc_number, **self.calc_kwargs)
        geom.set_calculator(calc)
        print(f"Setting files on calculator_{calc_number:03d}:")
        for ext in exts:
            file_ext = self.file_by_ext(files, ext)
            setattr(calc, ext, file_ext)
            print(f"\t{file_ext}")

    def set_files_on_calculators(self, geoms, files_dict, calc_class,
                                exts):
        for i, geom in enumerate(geoms):
            calc_number, cycle_number = i, 0
            self.set_files_on_calculator(geom, files_dict, calc_class, exts,
                                         calc_number, cycle_number)

    def set_files_from_dir(self, geom, path, calc_number):
        func = self.files_from_dir_dict[self.calc_key]
        func(geom, path, calc_number)

    def set_orca_files_from_dir(self, geom, path, calc_number):
        exts = self.orca_exts
        files_list = self.discover_files_in_dir(path, exts)
        files_dict = {
            (calc_number, 0): files_list,
        }
        self.set_files_on_calculator(geom, files_dict, ORCA, exts, calc_number)
        geom.calculator.store_overlap_data(geom.atoms, geom.coords)

    def set_orca_files(self, geoms, files_dict):
        self.set_files_on_calculators(geoms, files_dict, ORCA, self.orca_exts)
        for geom in geoms:
            geom.calculator.store_overlap_data(geom.atoms, geom.coords)

    def set_gaussian16_files_from_dir(self, geom, path, calc_number):
        exts = self.gaussian_exts

        files_list = self.discover_files_in_dir(path, exts)
        log_file, *files_list = files_list
        files_dict = {
            (calc_number, 0): files_list,
        }
        exts_without_log = exts[1:]
        assert "log" not in exts_without_log
        self.set_files_on_calculator(geom, files_dict, Gaussian16,
                                     exts_without_log,
                                     calc_number
        )
        log_path = Path(log_file)
        nmos, roots = geom.calculator.parse_log(log_path)
        calc = geom.calculator
        calc.nmos = nmos
        calc.roots = roots
        calc.store_overlap_data(geom.atoms, geom.coords)

    def set_g16_files(self, geoms, files_dict):
        exts = ("fchk", "dump_635r")
        self.set_files_on_calculators(geoms, files_dict, Gaussian16, exts)

        first_log = Path(self.file_by_ext(files_dict[(0, 0)], ".log"))
        nmos, roots = geoms[0].calculator.parse_log(first_log)
        for geom in geoms:
            calc = geom.calculator
            calc.nmos = nmos
            calc.roots = roots
            calc.store_overlap_data(geom.atoms, geom.coords)

    def set_turbo_files(self, geoms, files_dict):
        exts = ("mos", "ciss_a", "out", "control")
        self.set_files_on_calculators(geoms, files_dict, Turbomole, exts)

        for geom in geoms:
            calc = geom.calculator
            if hasattr(calc, "ciss_a"):
                calc.td_vec_fn = calc.ciss_a
            elif hasattr(calc, "ccres"):
                calc.td_vec_fn = calc.ccres
            calc.store_overlap_data(geom.atoms, geom.coords)
    
    def restore_calculators(self, geoms):
        files_dict = self.discover_files(self.path)
        unique_calculators = set([calc_num for calc_num, cycle_num in files_dict])
        # assert len(unique_calculators) <= len(geoms), ("Number of discovered "
            # f"unique calculators ({len(unique_calculators)}) is bigger than the "
            # f"number of discovered geometries ({len(geoms)})."
        # )
        print(f"Found {len(unique_calculators)} unique calculators.")
        print(f"Found {len(geoms)} geometries.")
        calc_num = min(len(unique_calculators), len(geoms))
        setter_func = self.setter_dict[self.calc_key]
        setter_func(geoms[:calc_num], files_dict)
        print(f"Restored {calc_num} calculators.")
        return calc_num

    def similar_overlaps(self, overlaps_for_state, ovlp_thresh=.1, diff_thresh=.2):
        """Return True if overlaps for a state are similar."""
        # Find overlaps above ovlp_thresh
        above_inds = np.where(np.abs(overlaps_for_state) > ovlp_thresh)[0]
        # Unambiguous assignment. There is a one to one correspondence between
        # the states.
        if len(above_inds) == 1:
            return False
        # Given the a full row containing overlaps this may evaluate to True if
        # something went wrong and the overlaps ARE that small. Given only a subset
        # of a full row, e.g. when only the first N states are considered this may
        # evaluate to True if the index of the current state lies below N. E.g. if we
        # check state 6, but got only the overlaps from state 1 to 5.
        elif len(above_inds) == 0:
            return False

        above_thresh = np.abs(overlaps_for_state[above_inds])
        max_ovlp_ind = above_thresh.argmax()
        max_ovlp = above_thresh[max_ovlp_ind]
        without_max = np.delete(above_thresh, max_ovlp_ind)
        # Consider the differences between the maximum overlap and the smaller ones.
        diffs = np.abs(max_ovlp - without_max)
        # Return True if any difference is below the threshold
        return any(diffs < diff_thresh)

    def get_ovlp_func(self, ovlp_type, double_mol=False, recursive=False,
                      consider_first=None):
        def wf_ovlp(calc1, calc2, ao_ovlp):
            ovlp_mats = calc1.wfow.overlaps_with(calc2.wfow, ao_ovlp=ao_ovlp)
            ovlp_mat = ovlp_mats[0]
            return ovlp_mat

        def tden_ovlp(calc1, calc2, ao_ovlp):
            return calc1.tdens_overlap_with_calculator(calc2,
                                                       ao_ovlp=ao_ovlp)
        ovlp_dict = {
            "wf": wf_ovlp,
            "tden": tden_ovlp,
        }
        valid_ovlps = "/".join([str(k) for k in ovlp_dict.keys()])
        assert ovlp_type in ovlp_dict.keys(), \
            f"Invalid ovlp_type! Valid keys are {valid_ovlps}."

        ovlp_func_ = ovlp_dict[ovlp_type]

        def ovlp_func(geoms, i, j, depth=2, ao_ovlp=None):
            ith_geom = geoms[i]
            jth_geom = geoms[j]
            ith_calc = geoms[i].calculator
            jth_calc = geoms[j].calculator
            icn = ith_calc.calc_number
            jcn = jth_calc.calc_number
            if double_mol:
                true_ovlp_mat_fn = f"ao_ovlp_true_{icn:03d}_{jcn:03d}"
                try:
                    ao_ovlp = np.loadtxt(true_ovlp_mat_fn)
                    self.logger.info(f"Using true AO overlaps from {true_ovlp_mat_fn}.")
                except:
                    self.logger.info("Doing double molecule calculation to get "
                                     "AO overlaps."
                    )
                    ao_ovlp = jth_geom.calc_double_ao_overlap(ith_geom)
                    np.savetxt(f"ao_ovlp_true_{icn:03d}_{jcn:03d}", ao_ovlp)
            self.log(f"Calculationg overlaps for steps {icn:03d} and {jcn:03d}.")
            ovlp_mat = ovlp_func_(ith_calc, jth_calc, ao_ovlp)

            self.log(ovlp_mat)

            ovlp_mat_fn = f"{ovlp_type}_ovlp_mat_{icn:03d}_{jcn:03d}"
            np.savetxt(self.path / ovlp_mat_fn, ovlp_mat)

            similar = any(
                [self.similar_overlaps(per_state)
                 for per_state in ovlp_mat[:,:consider_first]]
            )
            if similar:
                self.log( "Some entries of the overlap matrix between steps "
                         f"{icn:03d} and {jcn:03d} are very similar!")
            if recursive and similar and (i > 0) and depth > 0:
                self.log(f"Comparing {icn-1:03d} and {jcn:03d} now, "
                         f"because steps {icn:03d} and {jcn:03d} were "
                          "too similar."
                )
                return ovlp_func(geoms, i-1, j, depth-1)
            return ovlp_mat

        return ovlp_func

    @np_print
    def overlaps_for_geoms(self, geoms, ovlp_type="wf", double_mol=False,
                           recursive=False, consider_first=None, skip=0):
        # if skip > 0 and recursive:
            # raise Exception("recursive = True and skip > 0 can't be used "
                            # "together."
            # )
        ovlp_func = self.get_ovlp_func(ovlp_type, double_mol, recursive,
                                       consider_first)

        if double_mol:
            assert hasattr(geoms[0].calculator, "run_double_mol_calculation"), \
                   "Double molecule calculation not implemented for " \
                  f"{self.calc_key}."

        self.log(f"Doing {ovlp_type.upper()}-overlaps.")

        inds_list = list()
        ovlp_mats = list()
        is_similar = lambda ovlp_mat: any([self.similar_overlaps(per_state)
                                           for per_state in ovlp_mat[:,:consider_first]]
        )
        for i in range(len(geoms)-1):
            # We can be sure that i is always a valid index.
            j = i+(1+skip)
            if self.ovlp_with == "first":
                i = 0
            elif self.prev_n:
                i = max(i - self.prev_n, 0)
            if j >= len(geoms):
                break
            ovlp_mat = ovlp_func(geoms, i,  j)
            ovlp_mats.append(ovlp_mat)
            index_array = index_array_from_overlaps(ovlp_mat)
            inds_list.append(index_array)
            self.log(index_array)
        inds_arr = np.array(inds_list)
        ovlp_mats = np.array(ovlp_mats)
        max_ovlp_inds_fn = f"{ovlp_type}_max_ovlp_inds"
        ovlp_mats_fn = f"{ovlp_type}_ovlp_mats"
        np.savetxt(self.path / max_ovlp_inds_fn, inds_arr, fmt="%i")
        np.save(ovlp_mats_fn, ovlp_mats)
        self.log("")
        self.log("Max overlap indices.")
        for i, row in enumerate(inds_arr):
            self.log(f"Step {i:02d}: {row}")


if __name__ == "__main__":
    # path = Path("/scratch/programme/pysisyphus/tests_staging/test_diabatizer/cb3_def2svp")
    path = Path("/scratch/programme/pysisyphus/tests_staging/test_diabatizer/cb3_def2svp/first_five")
    calc_kwargs = {
        "keywords": "CAM-B3LYP def2-SVP RIJCOSX D3BJ TightSCF",
        "blocks": "%tddft nroots 5 tda false end %maxcore 1000",
        "track": True,
        "pal": 4,
    }
    calc_key = "orca"
    ovl = Overlapper(path, calc_key, calc_kwargs)
    geoms = ovl.discover_geometries(path)
    # files_dict = dia.discover_files(path)
    # dia.restore_calculators("orca")
    ovl.overlaps(geoms)
