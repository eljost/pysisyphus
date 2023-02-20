import re
import shutil
import textwrap

import numpy as np

from pysisyphus.calculators.Calculator import Calculator

class CFOUR(Calculator):

    conf_key = "cfour"

    def __init__(
            self,
            cfour_input,
            keep_molden=True,
            **kwargs,
    ):
        super().__init__(**kwargs)

        self.cfour_input = cfour_input

        self.inp_fn = 'ZMAT'
        self.out_fn = 'out.log'
        self.to_keep = ("out.log", "density:__den.dat")
        self.initden = None
        if keep_molden:
            self.to_keep = self.to_keep + ("MOLDEN*",)

        self.base_cmd = self.get_cmd("cmd")
        self.parser_funcs = {
            "energy": self.parse_energy,
            "grad": self.parse_gradient,
        }

        self.float_regex = r"([-]?\d+\.\d+)"  ## CFOUR doesn't use scientific notation for final energy or gradient.

    def prepare(self, inp):
        path = super().prepare(inp)
        if self.initden:
            shutil.copy(self.initden,f"{path}/initden.dat")
        return path

    def keep(self, path):
        kept_fns = super().keep(path)
        try:
            self.initden = kept_fns["density"]
        except KeyError:
            self.log("den.dat not found!")
            return

    def prepare_input(self, atoms, coords, calc_type):
        xyz_string = self.prepare_coords(atoms, coords, angstrom=False)
        cfour_keyword_string = '\n'.join(f"{key.upper()}={str(value).upper()}" for key, value in self.cfour_input.items())
        grad_string = "DERIV_LEVEL=1\n" if calc_type == "grad" else ""

        ### Note: CFOUR abhors blank lines between keywords. Make sure no extra whitespace is added between keyword lines.
        ### First dedent, then format.
        input_str = textwrap.dedent("""\
        CFOUR calculation
        {xyz_string}

        *CFOUR(UNITS=BOHR
        COORDINATES=CARTESIAN
        {grad_string}{cfour_keyword_string})

        """).format(xyz_string=xyz_string, grad_string=grad_string, cfour_keyword_string=cfour_keyword_string)

        return input_str

    def get_energy(self, atoms, coords):
        return self.run_calculation(atoms, coords, "energy")

    def get_forces(self, atoms, coords):
        return self.run_calculation(atoms, coords, "grad")

    def parse_energy(self, path):
        energy_fn = path / self.out_fn
        with open(energy_fn) as handle:
            text = handle.read()

        regex = "\s*The final electronic energy is\s*" + self.float_regex
        mobj = re.search(regex, text, re.DOTALL)
        energy = float(mobj.groups()[0])

        return energy

    def parse_gradient(self, path):
        ## Adapted from OpenMolcas calculator
        results = {}
        gradient_fn = path / self.out_fn
        with open(gradient_fn) as handle:
            text = handle.read()

        # Search for the block containing the gradient table
        regex = "gradient from JOBARC(.+)--executable xjoda finished"
        floats = [self.float_regex for i in range(3)]
        line_regex = r"^\s*" + r"\s*".join(floats) + r"\s*$"  ## Nothing but floats on the gradient lines

        mobj = re.search(regex, text, re.DOTALL)
        gradient = list()
        for line in mobj.groups()[0].split("\n"):
            # Now look for the lines containing the gradient
            mobj = re.match(line_regex, line.strip())
            if not mobj:
                continue
            gradient.append(mobj.groups())
        gradient = np.array(gradient, dtype=float)
        gradient_rotated = self.rotate_gradient(text, gradient).flatten()

        energy = self.parse_energy(path)
        results["energy"] = energy
        results["forces"] = -gradient_rotated

        return results

    def rotate_gradient(self, text, gradient):
        cfour_coords_3d = self.read_geom(text)
        pysis_coords_3d = np.reshape(self.input_coords, (-1,3))

        cfour_centroid = self.calc_centroid(cfour_coords_3d)
        pysis_centroid = self.calc_centroid(pysis_coords_3d)

        rot_matrix = self.calc_rot_matrix(cfour_coords_3d, pysis_coords_3d, cfour_centroid, pysis_centroid)
        t = self.calc_translation(rot_matrix, cfour_centroid, pysis_centroid)

        return (gradient @ rot_matrix.T) + t

    def read_geom(self, text):
        regex = r"Coordinates used in calculation \(QCOMP\)(.+)Interatomic distance matrix"
        floats =  [self.float_regex for i in range(3)]
        line_regex = r"\s*[A-Z]+\s*\d+\s*" + r"\s*".join(floats) ## Element symbol, atomic number, then x, y, z in bohr

        mobj = re.search(regex, text, re.DOTALL)
        geom = list()
        for line in mobj.groups()[0].split("\n"):
            mobj = re.match(line_regex, line.strip())
            if not mobj:
                continue
            geom.append(mobj.groups())
        geom_3d = np.array(geom, dtype=float)

        return geom_3d

    def calc_centroid(self, coords_3d):
        centroid = np.sum(coords_3d, axis=0)/coords_3d.shape[0]
        return centroid

    def calc_rot_matrix(self, cfour_coords_3d, pysis_coords_3d, cfour_centroid, pysis_centroid):
        H = (cfour_coords_3d - cfour_centroid).T @ (pysis_coords_3d - pysis_centroid)
        assert H.shape == (3,3)
        U, S, V = np.linalg.svd(H)
        R = V @ U.T

        ## Cover corner case
        if np.linalg.det(R) < 0:
            U, S, V = np.linalg.svd(R)
            V[:,2] = V[:,2]*-1
            R = V @ U.T

        return R

    def calc_translation(self, rotation_matrix, cfour_centroid, pysis_centroid):
        return pysis_centroid - rotation_matrix @ cfour_centroid

    def run_calculation(self, atoms, coords, calc_type):
        self.input_coords = coords  ## For use later to rotate CFOUR gradient to the pysisyphus frame
        inp = self.prepare_input(atoms, coords, calc_type)
        results = self.run(inp, calc=calc_type)
        return results
