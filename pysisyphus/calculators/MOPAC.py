#!/usr/bin/env python3

from math import sqrt
import re
import textwrap

import numpy as np

from pysisyphus.constants import BOHR2ANG, AU2EV, AU2KCALMOL
from pysisyphus.calculators.Calculator import Calculator


np.set_printoptions(suppress=True, precision=4, linewidth=120)


class MOPAC(Calculator):
    """http://openmopac.net/manual/"""

    conf_key = "mopac"

    MULT_STRS = {
        1: "SINGLET",
        2: "DOUBLET",
        3: "TRIPLET",
        4: "QUARTET",
        5: "QUINTET",
        6: "SEXTET",
        7: "SEPTET",
        8: "OCTET",
    }

    CALC_TYPES = {
        "energy": "1SCF",
        "gradient": "1SCF GRADIENTS",
        "hessian": "DFORCE FORCE LET",
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.uhf = "UHF" if self.mult != 1 else ""

        _ = "mopac"
        self.inp_fn = f"{_}.mop"
        self.out_fn = f"{_}.out"
        self.aux_fn = f"{_}.aux"
        self.to_keep = ("mop", "out", "arc", "aux")

        self.parser_funcs = {
            "energy": self.parse_energy,
            "grad": self.parse_grad,
            "hessian": self.parse_hessian,
        }

        self.base_cmd = self.get_cmd("cmd")

        """
        1SCF: Do only SCF
        AUX: Creates a checkpoint file
        NOREO: Dont reorient geometry
        """

        self.inp = textwrap.dedent("""
        NOSYM PM7 {mult} CHARGE={charge} {calc_type} {uhf} THREADS={pal} AUX(6,PRECISION=9) NOREOR

 
        {coord_str}
        """).strip()

    def prepare_coords(self, atoms, coords, opt=False):
        coords = coords.reshape(-1, 3) * BOHR2ANG
        # Optimization flag for coordinate
        of = 1 if opt else 0
        coord_str = "\n".join(
                [f"{a} {c[0]: 10.08f} {of} {c[1]: 10.08f} {of} {c[2]: 10.08f} {of}"
                 for a, c in zip(atoms, coords)]
        )
        return coord_str

    def prepare_input(self, atoms, coords, calc_type, opt=False):
        coord_str = self.prepare_coords(atoms, coords, opt)

        inp = self.inp.format(
                charge=self.charge,
                mult=self.MULT_STRS[self.mult],
                uhf=self.uhf,
                calc_type=self.CALC_TYPES[calc_type],
                coord_str=coord_str,
                pal=self.pal,
                # mem=self.mem,
        )
        return inp

    def get_energy(self, atoms, coords):
        calc_type = "energy"
        inp = self.prepare_input(atoms, coords, calc_type)
        # with open("inp.mop", "w") as handle:
            # handle.write(inp)
        # import sys; sys.exit()
        results = self.run(inp, calc="energy")
        return results

    def get_forces(self, atoms, coords):
        calc_type = "gradient"
        inp = self.prepare_input(atoms, coords, calc_type, opt=True)
        results = self.run(inp, calc="grad")
        return results

    def get_hessian(self, atoms, coords):
        calc_type = "hessian"
        inp = self.prepare_input(atoms, coords, calc_type, opt=True)
        results = self.run(inp, calc="hessian")
        return results

    def read_aux(self, path):
        with open(path / self.aux_fn) as handle:
            text = handle.read()
        return text

    def parse_energy(self, path):
        # with open(path / self.out_fn) as handle:
            # text = handle.read()
        # energy_re = "TOTAL ENERGY\s+=\s+([\-\.\d]+) EV"

        text = self.read_aux(path)
        energy_re = "HEAT_OF_FORMATION:KCAL/MOL=([\d\-D+\.]+)"
        mobj = re.search(energy_re, text)
        energy = float(mobj[1].replace("D", "E")) / AU2KCALMOL

        result = {
            "energy": energy,
        }
        return result


    def parse_grad(self, path):
        text = self.read_aux(path)
        grad_re = "GRADIENTS:KCAL.+$\s+(.+)$"
        mobj = re.search(grad_re, text, re.MULTILINE)
        # Gradients are given in kcal*mol/angstrom
        gradients = np.array(mobj[1].split(), dtype=float)
        # Convert to hartree/bohr
        gradients /= AU2KCALMOL / BOHR2ANG

        forces = -gradients
        result = {
            "forces": forces,
        }
        result.update(self.parse_energy(path))
        return result

    def parse_hessian(self, path):
        text = self.read_aux(path)

        mode_re = "NORMAL_MODES\[\d+\]=\s+([\s\.\-\d]+)\s+#Warning"
        normal_modes = re.search(mode_re, text)[1].strip().split()
        normal_modes = np.array(normal_modes, dtype=float)
        coord_num = int(sqrt(normal_modes.size))
        assert normal_modes.size % coord_num == 0
        normal_modes = normal_modes.reshape(-1, coord_num)

        hess_re = " #  Lower half triangle only\s+([\s\.\-\d]+)\s+NORMAL_MODE"
        tril_hess = re.search(hess_re, text)[1].strip().split()
        tril_hess = np.array(tril_hess, dtype=float)
        hessian = np.zeros((coord_num, coord_num))
        # TODO: get coord_num from somewhere else?!
        tril_indices = np.tril_indices(coord_num)
        hessian[tril_indices] = tril_hess

        triu_indices = np.triu_indices(coord_num, k=1)
        hessian[triu_indices] = hessian.T[triu_indices]
        # Hessian is given in millidyne/angstrom and massweighted

        result = {
            "hessian": hessian,
        }
        result.update(self.parse_energy(path))
        return result

    def __str__(self):
        return f"MOPAC({self.name})"
