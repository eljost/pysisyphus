import re
import textwrap

import numpy as np

from pysisyphus.constants import BOHR2ANG, AU2KCALPERMOL
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.helpers_pure import file_or_str


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

    METHODS = [
        m.lower()
        for m in "AM1 PM3 PM6 PM6-DH2 PM6-D3 PM6-DH+ PM6-DH2 PM6-DH2X "  # lgtm [py/non-iterable-in-for-loop]
        "PM6-D3H4 PM6-D3H4X PM7 PM7-TS".split()
    ]

    def __init__(self, method="PM7", **kwargs):
        super().__init__(**kwargs)

        self.method = method
        assert (
            self.method.lower() in self.METHODS
        ), f"Invalid method={self.method}! Supported methods are ({self.METHODS})"

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

        self.base_cmd = self.get_cmd()

        """
        1SCF: Do only SCF
        AUX: Creates a checkpoint file
        NOREO: Dont reorient geometry
        """

        self.inp = textwrap.dedent(
            """
        NOSYM {method} {mult} CHARGE={charge} {calc_type} {uhf} THREADS={pal} AUX(6,PRECISION=9) NOREOR

 
        {coord_str}
        """
        ).strip()

        self.log(f"Created MOPAC calculator using the '{self.method}' method.")

    def prepare_coords(self, atoms, coords, opt=False):
        coords = coords.reshape(-1, 3) * BOHR2ANG
        # Optimization flag for coordinate
        of = 1 if opt else 0
        coord_str = "\n".join(
            [
                f"{a} {c[0]: 10.08f} {of} {c[1]: 10.08f} {of} {c[2]: 10.08f} {of}"
                for a, c in zip(atoms, coords)
            ]
        )
        return coord_str

    def prepare_input(self, atoms, coords, calc_type, opt=False):
        coord_str = self.prepare_coords(atoms, coords, opt)

        inp = self.inp.format(
            method=self.method,
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
        return self.parse_energy_from_aux(self.read_aux(path))

    @staticmethod
    @file_or_str(".aux", method=False)
    def parse_energy_from_aux(text):
        energy_re = r"HEAT_OF_FORMATION:KCAL/MOL=([\d\-D+\.]+)"
        mobj = re.search(energy_re, text)
        energy = float(mobj[1].replace("D", "E")) / AU2KCALPERMOL

        result = {
            "energy": energy,
        }
        return result

    def parse_grad(self, path):
        text = self.read_aux(path)
        # grad_re = "GRADIENTS:KCAL.+$\s+(.+)$"
        # mobj = re.search(grad_re, text, re.MULTILINE)
        grad_re = r"GRADIENTS:KCAL/MOL/ANGSTROM\[\d+]=\s+(.+)\s+OVERLAP_MATRIX"
        mobj = re.search(grad_re, text, re.DOTALL)
        # Gradients are given in kcal*mol/angstrom
        gradients = np.array(mobj[1].split(), dtype=float)
        # Convert to hartree/bohr
        gradients = gradients / AU2KCALPERMOL / BOHR2ANG

        forces = -gradients
        result = {
            "forces": forces,
        }
        result.update(self.parse_energy(path))
        return result

    def parse_hessian(self, path):
        return self.parse_hessian_from_aux(self.read_aux(path))

    @staticmethod
    @file_or_str(".aux", method=False)
    def parse_hessian_from_aux(text):
        # Parse employed masses, as the given hessian is mass-weighted
        # and we have to un-weigh it.
        mass_re = re.compile(
            r"ISOTOPIC_MASSES\[(\d+)\]=\s*(.+?)ROTAT_CONSTS", re.DOTALL
        )
        # mobj = re.search(mass_re, text, re.MULTILINE)
        mass_mobj = mass_re.search(text)
        masses = np.array(mass_mobj[2].strip().split(), dtype=float)
        # This matrix is used to un-weigh the hessian
        M = np.diag(np.sqrt(np.repeat(masses, 3)))
        # For N atoms we expect 3N cartesian coordinates
        coord_num = masses.size * 3

        hess_re = r" #  Lower half triangle only\s+([\s\.\-\d]+)\s+NORMAL_MODE"
        tril_hess = re.search(hess_re, text)[1].strip().split()
        tril_hess = np.array(tril_hess, dtype=float)
        assert tril_hess.size == sum(range(coord_num + 1))
        hessian_m = np.zeros((coord_num, coord_num))
        tril_indices = np.tril_indices(coord_num)
        hessian_m[tril_indices] = tril_hess

        triu_indices = np.triu_indices(coord_num, k=1)
        hessian_m[triu_indices] = hessian_m.T[triu_indices]

        # Hessian is given in mdyn/(Å*amu).
        # In a first step we have to unweigh the hessian using the matrix
        # built from the parsed masses.
        hessian = M @ hessian_m @ M
        # Then we have to convert mdyn/Å to Hartree/Bohr²
        #     mdyn/Å = 100 kg/s²
        #     Hartree/Bohr² ~  1556.8931 kg/s²
        #
        #     1 mydn/Å * (100 / 1556.8931 Hartree/Bohr² * Å/mydn) = 0.06423 Hartree/Bohr²
        hessian *= 0.06423

        energy = MOPAC.parse_energy_from_aux(text)["energy"]

        result = {
            "energy": energy,
            "hessian": hessian,
        }
        return result

    def __str__(self):
        return f"MOPAC({self.name})"
