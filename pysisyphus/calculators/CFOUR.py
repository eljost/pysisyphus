import re
import shutil
import textwrap

import numpy as np

from pysisyphus.calculators.Calculator import Calculator

cfour_float_regex = r"([-]?\d+\.\d+)" # CFOUR doesn't use scientific notation for final energy or gradient.

def parse_cfour_energy(out_fn):
    with open(out_fn) as handle:
        text = handle.read()

    regex = "\s*The final electronic energy is\s*" + cfour_float_regex
    mobj = re.search(regex, text, re.DOTALL)
    energy = float(mobj.groups()[0])

    return {
        "energy": energy
    }

def parse_cfour_gradient(grd_fn, out_fn, coords_pysis_frame_3d):
    results = {}
    gradient = np.loadtxt(grd_fn, skiprows=1)
    natoms = int(gradient.shape[0] / 2)
    gradient = gradient[natoms:, 1:]

    coords_comp_frame_3d = read_cfour_geom(out_fn)
    gradient_rotated = rotate_gradient(gradient, coords_pysis_frame_3d, coords_comp_frame_3d).flatten()

    energy = parse_cfour_energy(out_fn)["energy"]
    results["energy"] = energy
    results["forces"] = -gradient_rotated

    return results

def read_cfour_geom(out_fn):
    with open(out_fn) as handle:
        text = handle.read()
    
    regex = r"Coordinates used in calculation \(QCOMP\)(.+)Interatomic distance matrix"
    floats =  [cfour_float_regex] * 3
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

def rotate_gradient(gradient, coords_pysis_frame_3d, coords_comp_frame_3d):
    # Following http://nghiaho.com/?page_id=671
    # Permalink: https://web.archive.org/web/20230128004529/http://nghiaho.com/?page_id=671

    comp_frame_centroid = calc_centroid(coords_comp_frame_3d)
    pysis_frame_centroid = calc_centroid(coords_pysis_frame_3d)

    rot_matrix = calc_rot_matrix(coords_comp_frame_3d, coords_pysis_frame_3d, comp_frame_centroid, pysis_frame_centroid)

    return (gradient @ rot_matrix.T)

def calc_centroid(coords_3d):
    centroid = np.sum(coords_3d, axis=0)/coords_3d.shape[0]
    return centroid

def calc_rot_matrix(source_coords_3d, target_coords_3d, source_centroid, target_centroid):
    H = (source_coords_3d - source_centroid).T @ (target_coords_3d - target_centroid)
    assert H.shape == (3,3)
    U, S, V = np.linalg.svd(H)
    R = V @ U.T

    # Cover corner case
    if np.linalg.det(R) < 0:
        U, S, V = np.linalg.svd(R)
        V[:,2] = V[:,2]*-1
        R = V @ U.T

    return R

class CFOUR(Calculator):

    conf_key = "cfour"

    def __init__(
            self,
            cfour_input,
            wavefunction_dump=True,
            initden_file=None,
            **kwargs,
    ):
        """CFOUR calculator.

        Wrapper handling CFOUR ground state energy and gradient calculations.

        Parameters
        ----------
        cfour_input : dict
            CFOUR keywords and values. Note: "on" must be encapsulated in quotes to avoid being translated to True by YAML.
        keep_molden : bool, optional
            Whether or not to keep ground state SCF orbitals for each geometry step.
        initden_file: str, optional
            Path to an input initden file for use as a guess SCF density.
        """
        super().__init__(**kwargs)

        self.cfour_input = cfour_input

        self.inp_fn = 'ZMAT'
        self.out_fn = 'out.log'
        self.gradient_fn = 'GRD'
        self.to_keep = ("out.log", "density:__den.dat")
        self.initden = None
        self.wavefunction_dump = wavefunction_dump

        if self.wavefunction_dump:
            self.to_keep = self.to_keep + ("MOLDEN*",)
        if initden_file:
            self.initden = initden_file

        self.base_cmd = self.get_cmd("cmd")
        self.parser_funcs = {
            "energy": self.parse_energy,
            "grad": self.parse_gradient,
        }

        # Convert pysisyphus keywords to CFOUR keywords
        # Explicitly not adding support for managing CFOUR parallelism - 
        # user-configured cmd script in .pysisyphusrc should set either OMP_NUM_THREADS 
        # or the environment variable they defined at compile-time for MPI parallelism
        cfour_input['MEM_SIZE'] = self.mem
        cfour_input['MEM_UNIT'] = 'MB'
        cfour_input['CHARGE'] = self.charge
        cfour_input['MULT'] = self.mult

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

    def prepare_input(self, atoms, coords, calc_type):
        xyz_string = self.prepare_coords(atoms, coords)
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

    def prepare_coords(self, atoms, coords):
        """Get 3d coords in Bohr

        Reshape internal 1d coords to 3d.

        Parameters
        ----------
        atoms : iterable
            Atom descriptors (element symbols).
        coords: np.array, 1d
            1D-array holding coordinates in Bohr.

        Returns
        -------
        coords: np.array, 3d
            3D-array holding coordinates in Bohr.
        """
        coords = coords.reshape(-1, 3)
        coords = "\n".join(
            [
                "{} {:10.08f} {:10.08f} {:10.08f}".format(a, *c)
                for a, c in zip(atoms, coords)
            ]
        )
        return coords

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self.run_calculation(atoms, coords, "energy", **prepare_kwargs)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self.run_calculation(atoms, coords, "grad", **prepare_kwargs)

    def parse_energy(self, path):
        return parse_cfour_energy(path / self.out_fn)

    def parse_gradient(self, path):
        return parse_cfour_gradient(path / self.gradient_fn, path / self.out_fn, self.pysis_frame_coords.reshape((-1,3)))

    def run_calculation(self, atoms, coords, calc_type, **prepare_kwargs):
        self.pysis_frame_coords = coords  # For use later to rotate CFOUR gradient to the pysisyphus frame
        inp = self.prepare_input(atoms, coords, calc_type, **prepare_kwargs)
        results = self.run(inp, calc=calc_type)
        return results

    def __str__(self):
        return f"CFOUR({self.name})"