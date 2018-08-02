#!/usr/bin/env python3

import itertools as it
import os
from pathlib import Path
import re

from natsort import natsorted
import numpy as np

from pysisyphus.calculators.Gaussian09 import Gaussian09
from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.calculators.Turbomole import Turbomole
from pysisyphus.helpers import geom_from_xyz_file


class Overlapper:
    orca_exts = ("out", "gbw", "cis")

    def __init__(self, path, calc_key=None, calc_kwargs=None):
        self.path = Path(path)
        self.calc_key = calc_key
        self.calc_kwargs = calc_kwargs
        self.calc_kwargs["out_dir"] = path

        self.setter_dict = {
            "g09": self.set_g16_files,
            "g16": self.set_g16_files,
            "orca": self.set_orca_files,
            "turbo": self.set_turbo_files,
        }
        self.files_from_dir_dict = {
            "orca": self.set_orca_files_from_dir,
        }

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

    def disover_files_in_dir(self, path, exts):
        files_list = list()
        for ext in exts:
            glob = f"*{ext}"
            fns = [_ for _ in path.glob(glob)]
            assert len(fns) == 1
            fn = str(fns[0])
            files_list.append(fn)
        return files_list

    def discover_geometries(self, path):
        xyz_fns = natsorted(path.glob("*.xyz"))
        geoms = [geom_from_xyz_file(xyz_fn) for xyz_fn in xyz_fns]
        self.restore_calculators(geoms, self.calc_key)

        return geoms

    def file_by_ext(self, iterable, ext):
        matches = [f for f in iterable if f.endswith(ext)]
        if len(matches) == 0:
            raise Excep
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
        files_list = self.disover_files_in_dir(path, exts)
        files_dict = {
            (calc_number, 0): files_list,
        }
        self.set_files_on_calculator(geom, files_dict, ORCA, exts, calc_number)
        geom.calculator.store_overlap_data(geom.atoms, geom.coords)

    def set_orca_files(self, geoms, files_dict):
        self.set_files_on_calculators(geoms, files_dict, ORCA, self.orca_exts)
        for geom in geoms:
            geom.calculator.store_overlap_data(geom.atoms, geom.coords)


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
        exts = ("mos", "ciss_a", "out")
        self.set_files_on_calculators(geoms, files_dict, Turbomole, exts)
    
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

    def overlaps_for_geoms(self, geoms, ovlp_type="wf", double_mol=False):
        def wf_ovlp(calc1, calc2, ao_ovlp):
            ovlp_mats = calc1.wfow.overlaps_with(calc2.wfow, ao_ovlp=ao_ovlp)
            ovlp_mat = ovlp_mats[0]
            return ovlp_mat

        def tden_ovlp(calc1, calc2, ao_ovlp=None):
            return calc1.tdens_overlap_with_calculator(calc2,
                                                       ao_ovlp=ao_ovlp)

        ovlp_dict = {
            "wf": wf_ovlp,
            "tden": tden_ovlp,
        }
        try:
            ovlp_func = ovlp_dict[ovlp_type]
        except KeyError as err:
            valid_ovlps = "/".join([str(k) for k in ovlp_dict.keys()])
            print(f"Please supply a valid ovlp_type ({valid_ovlps})! Exiting!")
            return
        if double_mol:
            assert hasattr(geoms[0].calculator, "run_double_mol_calculation"), \
                   "Double molecule calculation not implemented for " \
                  f"{self.calc_key}."

        print(f"Doing {ovlp_type.upper()}-overlaps.")

        inds_list = list()
        for i, geom in enumerate(geoms[1:]):
            ao_ovlp = None
            if double_mol:
                prev_geom = geoms[i]
                ao_ovlp = geom.calc_double_ao_overlap(prev_geom)
                np.savetxt(f"double_ovlp_{i:03d}", ao_ovlp)

            prev_calc = geoms[i].calculator
            cur_calc = geom.calculator
            print(f"Step {i:02d}, comparing Calculator {prev_calc.calc_number:02d} "
                  f"with Calculator {cur_calc.calc_number:02d}:"
            )
            ovlp_mat = ovlp_func(prev_calc, cur_calc, ao_ovlp)
            print(ovlp_mat)
            # Determine indices of maximum overlap
            index_array = prev_calc.index_array_from_overlaps(ovlp_mat)
            inds_list.append(index_array)
            ovlp_mat_fn = f"{ovlp_type}_ovlp_mat_{i:03d}"
            np.savetxt(self.path / ovlp_mat_fn, ovlp_mat)
            print(index_array)
        inds_arr = np.array(inds_list)
        max_ovlp_inds_fn = f"{ovlp_type}_max_ovlp_inds"
        np.savetxt(self.path / max_ovlp_inds_fn, inds_arr, fmt="%i")
        print()
        print("Max overlap indices.")
        print(inds_arr)


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
    geoms = ovl.discover_geometries(dia.path)
    # files_dict = dia.discover_files(path)
    # dia.restore_calculators("orca")
    ovl.overlaps(geoms)
