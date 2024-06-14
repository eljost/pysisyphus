from dataclasses import dataclass
from math import sqrt
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Optional
import warnings

from jinja2 import Template
import numpy as np
import pyparsing as pp

# from pysisyphus.calculators.cosmo_data import COSMO_RADII
from pysisyphus.calculators.OverlapCalculator import (
    GroundStateContext,
    OverlapCalculator,
)
from pysisyphus.calculators.parser import (
    parse_turbo_gradient,
    parse_turbo_ccre0_ascii,
    parse_turbo_mos,
    # parse_turbo_exstates,
    parse_turbo_exstates_re as parse_turbo_exstates,
)
from pysisyphus.helpers_pure import file_or_str, get_random_path
from pysisyphus.wavefunction.excited_states import (
    make_density_matrices_for_root,
    norm_ci_coeffs,
)


@dataclass
class ExSpectrumRoot:
    root: int
    sym: str
    exc_energy: float
    osc_vel: float
    osc_len: float


@file_or_str(".exspectrum")
def parse_exspectrum(text: str) -> list[ExSpectrumRoot]:
    """Parse root data from exspectrum file."""
    ex_roots = list()
    for line in text.strip().split("\n"):
        line = line.strip()
        is_comment = line.startswith("#")
        is_empty = line == ""
        if is_comment or is_empty:
            continue
        # root, sym, energy_au, energy_ev, energy_cm⁻¹, energy_nm, osc_vel, osc_len
        root, sym, exc_energy, *_, osc_vel, osc_len = line.split()
        root = int(root)
        exc_energy = float(exc_energy)
        osc_vel = float(osc_vel)
        osc_len = float(osc_len)
        ex_root = ExSpectrumRoot(
            root=root,
            sym=sym,
            exc_energy=exc_energy,
            osc_vel=osc_vel,
            osc_len=osc_len,
        )
        ex_roots.append(ex_root)
    return ex_roots


def index_strs_for_atoms(atoms):
    indices = dict()
    for i, atom in enumerate(atoms, 1):
        indices.setdefault(atom.lower(), list()).append(i)

    pairs = dict()
    for key, values in indices.items():
        pairs[key] = list()
        i_end = len(values) - 1
        start = None
        for i, v in enumerate(values):
            if start is None:
                start = v
            if (i == i_end) or (v + 1 != values[i + 1]):
                end = None if (start == v) else v
                pairs[key].append((start, end))
                start = None

    atom_strs = dict()
    for key, prs in pairs.items():
        tmp = list()
        for start, end in prs:
            if end is None:
                itm = str(start)
            else:
                itm = f"{start}-{end}"
            tmp.append(itm)
        atom_str = f"{key}  " + ",".join(tmp)
        atom_strs[key] = atom_str

    return atom_strs


def get_cosmo_data_groups(atoms, epsilon, rsolv=None, refind=None, dcosmo_rs=None):
    cosmo_dgs = {
        "cosmo": {
            "epsilon": epsilon,
        },
        "cosmo_out": "file=out.ccf",
        "cosmo_data": "file=cosmo_transfer.tmp",
    }
    if rsolv is not None:
        cosmo_dgs["cosmo"]["rsolv"] = rsolv
    if refind is not None:
        cosmo_dgs["cosmo"]["refind"] = refind
    # Transform atom list into TURBOMOLE-style compressed indices
    #   H, H, H, O, O, H, O
    # is transformed to
    #   h 1-3, 6
    #   o 4-5, 7
    # etc.
    index_strs = index_strs_for_atoms(atoms)
    cosmo_atom_strs = list()
    cosmo_atoms = dict()
    # This does not yet work; so we stick with the default radii.
    # for key, index_str in index_strs.items():
    # TODO: reenable COSMO_RADII import
    # radius = COSMO_RADII[key].radius
    # cosmo_atoms[index_str] = f"\nradius={radius:.6f}"
    # cosmo_dgs["cosmo_atoms"] = cosmo_atoms

    if dcosmo_rs:
        cosmo_dgs["dcosmo_rs"] = f"file={dcosmo_rs}"
    return cosmo_dgs


DATA_GROUP_TPL = Template(
    """
{% for dg in data_groups %}
${{ dg }}
    {% for key, value in data_groups[dg].items() %}
    {{ key }}{% if value %} ={{ value }}{% endif %}

    {% endfor %}
{% endfor %}
        """,
    trim_blocks=True,
    lstrip_blocks=True,
)
SIMPLE_INPUT_TPL = Template(
    """
$atoms
  basis={{ basis }}
$coord file=coord
$grad    file=gradient
$symmetry c1
$eht charge={{ charge }} unpaired={{ unpaired }}
{{ data_groups }}
$end
"""
)


def render_data_groups(raw_data_groups):
    data_groups = dict()
    for dg, kws in raw_data_groups.items():
        # datagroup w/o keywords
        if kws is None:
            data_groups[dg] = dict()
        # datagroup w/ one keyword on same line. Here we just merge the name and the keyword,
        # so they appear on the same line.
        elif type(kws) in (str, int, float):
            data_groups[f"{dg} {kws}"] = dict()
        # Otherwise a dict is expected
        elif type(kws) is dict:
            data_groups[dg] = kws
        else:
            raise Exception(f"Can't handle input '{dg}': '{kws}'!")
    return DATA_GROUP_TPL.render(data_groups=data_groups)


def control_from_simple_input(simple_inp, charge, mult, cosmo_kwargs=None):
    """Create control file from 'simple input'.

    See examples/opt/26_turbomole_simple_input/ for an example."""
    unpaired = mult - 1
    simple_inp = simple_inp.copy()
    basis = simple_inp.pop("basis")
    # Add memory statement, if not already present
    if "maxcor" not in simple_inp.keys():
        simple_inp["maxcor"] = "500 MiB per_core"

    data_groups = render_data_groups(simple_inp)

    rendered = SIMPLE_INPUT_TPL.render(
        basis=basis,
        charge=charge,
        unpaired=unpaired,
        data_groups=data_groups,
    )
    return rendered


def parse_frozen_nmos(text: str) -> tuple[list[tuple[int, int], tuple[int, int]], bool]:
    """Determine number of occ. & and virt. orbitals used in ES calculations."""

    frozen_re = re.compile(
        r"number of non-frozen orbitals\s+:\s+(?P<nmos>\d+)"
        r"\s+number of non-frozen occupied orbitals :\s+(?P<nocc>\d+)"
    )
    matches = frozen_re.findall(text)
    mo_nums = list()
    for nmos, nocc in matches:
        nmos = int(nmos)
        nocc = int(nocc)
        nvirt = nmos - nocc
        mo_nums.append((nocc, nvirt))

    restricted = len(mo_nums) == 1
    if restricted:
        mo_nums.append((0, 0))
    return mo_nums, restricted


@file_or_str(
    "ciss_a",
    "cist_a",
    "ucis_a",
    "sing_a",
    "trip_a",
    "unrs_a",
    "dipl_a",
    exact=True,
    add_exts=True,
)
def parse_ci_coeffs(
    text: str,
    nocc_a: int,
    nvirt_a: int,
    nocc_b: int,
    nvirt_b: int,
    restricted_same_ab=False,
):
    """Parse CI coefficients from escf/egrad calculations."""

    states_data = Turbomole.parse_td_vectors(text)
    expect_a = nocc_a * nvirt_a
    expect_b = nocc_b * nvirt_b
    shape_a = (nocc_a, nvirt_a)
    shape_b = (nocc_b, nvirt_b)

    Xa = np.zeros((len(states_data), *shape_a))
    Ya = np.zeros((len(states_data), *shape_a))
    Xb = np.zeros((len(states_data), *shape_b))
    Yb = np.zeros((len(states_data), *shape_b))

    # Whether a Y vector is present
    with_deexc = len(states_data[0]["vector"]) == 2 * (expect_a + expect_b)
    if with_deexc:
        XpYa = np.empty(shape_a)
        XmYa = np.empty(shape_a)
        XpYb = np.empty(shape_b)
        XmYb = np.empty(shape_b)

    for i, state_data in enumerate(states_data):
        coeffs = np.array(state_data["vector"])
        # TD-DFT/TD-HF
        if with_deexc:
            start = 0
            # X + Y
            XpYa[:] = coeffs[start : start + expect_a].reshape(shape_a)
            start += expect_a
            XpYb = coeffs[start : start + expect_b].reshape(shape_b)
            start += expect_b
            # X - Y
            XmYa[:] = coeffs[start : start + expect_a].reshape(shape_a)
            start += expect_a
            XmYb = coeffs[start : start + expect_b].reshape(shape_b)
            start += expect_b
            # Recover X and Y vectors from X + Y and X - Y
            Xa[i] = (XpYa + XmYa) / 2.0
            Ya[i] = XpYa - Xa[i]
            Xb[i] = (XpYb + XmYb) / 2.0
            Yb[i] = XpYb - Xb[i]
        # TDA/CIS
        else:
            Xa[i] = coeffs[:expect_a].reshape(shape_a)
            Xb[i] = coeffs[expect_a:].reshape(shape_b)
    if restricted_same_ab and Xb.size == 0:
        Xb = Xa.copy()
        Yb = Ya.copy()
    return Xa, Ya, Xb, Yb


def get_density_matrices_for_root(
    log_fn, vec_fn, root, rlx_vec_fn=None, Ca=None, Cb=None
):
    with open(log_fn) as handle:
        text = handle.read()

    # Number of non-frozen/active molecular orbitals
    ((nfocc_a, nfvirt_a), (nfocc_b, nfvirt_b)), restricted = parse_frozen_nmos(text)

    # Transition density
    Xa, Ya, Xb, Yb = parse_ci_coeffs(vec_fn, nfocc_a, nfvirt_a, nfocc_b, nfvirt_b)
    # Relaxed density correction; there is never a deexciation part
    if rlx_vec_fn:
        ov_corr_a, _, ov_corr_b, _ = parse_ci_coeffs(
            rlx_vec_fn, nfocc_a, nfvirt_a, nfocc_b, nfvirt_b
        )
        # We expect only one root in dipl_a
        assert ov_corr_a.shape[0] == 1
        ov_corr_a = ov_corr_a[0]
        assert ov_corr_b.shape[0] == 1
        ov_corr_b = ov_corr_b[0]
    else:
        ov_corr_a = None
        ov_corr_b = None

    return make_density_matrices_for_root(
        root - 1, restricted, Xa, Ya, Xb, Yb, ov_corr_a, ov_corr_b, Ca, Cb
    )


class TurbomoleGroundStateContext(GroundStateContext):
    def __enter__(self):
        super().__enter__()
        self.energy_cmd_bak = self.calc.energy_cmd
        self.calc.energy_cmd = self.calc.scf_cmd

    def __exit__(self, exc_type, exc_value, exc_traceback):
        super().__exit__(exc_type, exc_value, exc_traceback)
        self.calc.energy_cmd = self.energy_cmd_bak


# def get_overlap_data_from_base_name(base: str | Path, vec_fn: str | Path):
def get_overlap_data_from_base_name(base: str | Path, ci_suffix: str):
    base = Path(base)

    with open(base.with_suffix(".turbomole.out")) as handle:
        text = handle.read()
    ((nfocc_a, nfvirt_a), (nfocc_b, nfvirt_b)), _ = parse_frozen_nmos(text)

    """
    # At first it seemed clever to just look for the available files and pick the
    # first one to be found, but of course this is a stupid idea. When multiple
    # calculations are run in the same directory restricted and unrestricted
    # calculations can become mixed ...
    vec_exts = ("ciss_a", "ucis_a", "sing_a", "unrs_a")
    if vec_fn is None:
        for ext in vec_exts:
            if (vec_fn := base.with_suffix(f".{ext}")).exists():
                break
        else:
            raise Exception("Couldn't find any files containing CI coefficients")
    else:
        vec_fn = Path(vec_fn)
    """
    vec_fn = base.with_suffix(ci_suffix)

    # CI coefficients
    Xa, Ya, Xb, Yb = parse_ci_coeffs(
        vec_fn, nfocc_a, nfvirt_a, nfocc_b, nfvirt_b, restricted_same_ab=True
    )
    Xa, Ya, Xb, Yb = norm_ci_coeffs(Xa, Ya, Xb, Yb)

    mos_fn = base.with_suffix(".mos")
    if mos_fn.exists():
        with open(mos_fn) as handle:
            text = handle.read()
        Ca = parse_turbo_mos(text)
        Cb = Ca.copy()
    else:
        with open(base.with_suffix(".alpha")) as handle:
            text = handle.read()
        Ca = parse_turbo_mos(text)
        with open(base.with_suffix(".beta")) as handle:
            text = handle.read()
        Cb = parse_turbo_mos(text)

    return Ca, Cb, Xa, Ya, Xb, Yb


"""
def parse_shells(text, key: Literal["closed", "alpha", "beta"]):
    regex = re.compile(rf"{key} shells\s+(\w)\s*(\d+)-(\d+)")
    mobj = regex.search(text)
    from_ = int(mobj.groups(1))
    to_ = int(mobj.groups(2))
    assert from_ == 1
    return to_
"""


class Turbomole(OverlapCalculator):
    conf_key = "turbomole"
    _set_plans = (
        "out",
        "control",
        "alpha",
        "beta",
        "ccres",
        "exspectrum",
        "exstates",
        "mos",
        "mwfn_wf",
        ("ciss_a", "td_vec_fn"),
        ("cist_a", "td_vec_fn"),
        ("sing_a", "td_vec_fn"),
        ("trip_a", "td_vec_fn"),
        ("ucis_a", "td_vec_fn"),
        ("unrs_a", "td_vec_fn"),
        ("dipl_a", "rlx_vec_fn"),
    )

    def __init__(
        self,
        control_path=None,
        numfreq=False,
        simple_input=None,
        double_mol_path=None,
        cosmo_kwargs=None,
        wavefunction_dump=True,
        **kwargs,
    ):
        super().__init__(**kwargs)

        assert (control_path is not None) or (
            simple_input is not None
        ), "Please either provide a prepared 'control_path' or 'simple_input'!"

        # Handle simple input
        if simple_input:
            control_path = (self.out_dir / get_random_path("control_path")).absolute()
            self.log(
                "Set 'control_path' to '{control_path}'. Creating 'control' from simple input in it."
            )
            try:
                control_path.mkdir()
            except FileExistsError:
                # Clean directory if already exists
                for fn in control_path.iterdir():
                    fn.unlink()
            control_str = control_from_simple_input(
                simple_input, charge=self.charge, mult=self.mult
            )
            with open(control_path / "control", "w") as handle:
                handle.write(control_str)
        # End of simple input handling

        # Set provided control_path or use the one generated for simple_input
        self.control_path = Path(control_path).absolute()
        self.numfreq = numfreq

        self.double_mol_path = double_mol_path
        if self.double_mol_path:
            self.double_mol_path = Path(self.double_mol_path)
        # Check if the overlap matrix will be printed and assert
        # that no SCF iterations are done.
        if self.double_mol_path:
            with open(self.double_mol_path / "control") as handle:
                text = handle.read()
            assert re.search(r"\$intsdebug\s*sao", text) and re.search(
                r"\$scfiterlimit\s*0", text
            ), ("Please set " "$intsdebug sao and $scfiterlimit 0 !")
        self.cosmo_kwargs = cosmo_kwargs
        if self.cosmo_kwargs:
            assert (
                "epsilon" in self.cosmo_kwargs
            ), "If 'cosmo_kwargs' is given 'epsilon' must be specified!"
        self.wavefunction_dump = wavefunction_dump

        self.to_keep = (
            "control",
            "mos",
            "alpha",
            "beta",
            "out",
            "ciss_a",
            "cist_a",
            "ucis_a",
            "sing_a",
            "trip_a",
            "unrs_a",
            "dipl_a",
            "gradient",
            "__ccre*",
            "exstates",
            "coord",
            "mwfn_wf:wavefunction.molden",
            "input.xyz",
            "pc_gradients",
            "nprhessian",
            "exspectrum",
        )

        self.parser_funcs = {
            "energy": self.parse_energy,
            "force": self.parse_force,
            "hessian": self.parse_hessian,
            "double_mol": self.parse_double_mol,
            "noparse": lambda path: None,
        }

        # Turbomole uses the 'control' file implicitly
        self.inp_fn = ""
        self.out_fn = "turbomole.out"
        # MO coefficient files
        self.mos = None
        self.alpha = None
        self.beta = None

        # Read control file
        with open(self.control_path / "control") as handle:
            text = handle.read().lower()
        # Prepare base_cmd
        scf_cmd = "dscf"
        second_cmd = "grad"
        # Check for RI
        self.ri = ("$rij" in text) or ("$rik" in text)
        if self.ri:
            scf_cmd = "ridft"
            second_cmd = "rdgrad"
            self.log("Found RI calculation.")

        self.uhf = ("$uhf" in text) or (self.mult > 1)

        # TODO: drop this block below
        """
        try:
            # It seems as they changed whats written in the control file in version 7.7
            # self.set_occ_and_mo_nums(text)
            pass
        except TypeError:
            warnings.warn(
                "Parsing of occupied and virtual MO numbers failed! Disabling "
                "excited state tracking!"
            )
            self.track = False
        """

        exopt_present = "$exopt" in text
        geoopt_present = "$geoopt" in text
        # At most one of these flags (exopt/geoopt) can be present
        assert not (exopt_present and geoopt_present), (
            "Found $exopt and $ricc2 in the control file! $exopt is used "
            "for TD-DFT/TDA gradients whereas $ricc2 with 'geoopt ...' "
            "leads to ricc2 gradients. Please delete one of the keywords!"
        )

        self.td = False
        self.td_vec_fn = None
        self.ricc2 = False
        self.ricc2_opt = False
        # Check for excited state calculation
        if exopt_present:
            exopt_re = r"\$exopt\s*(\d+)"
            # Only use root from $exopt if no explicit root was given
            if self.root is None:
                try:
                    self.root = int(re.search(exopt_re, text)[1])
                except TypeError:
                    self.log(
                        f"Couldn't parse root from $exopt and no explicit root was given!"
                    )
            second_cmd = "egrad"
            self.prepare_td(text)
            self.td = True
        elif "$soes" in text:
            second_cmd = "escf"
            self.td = True
            self.prepare_td(text)
            if not exopt_present:
                warnings.warn(
                    "ES calculations w/ escf requested, but '$exopt' isn't present "
                    "in the control file. In gradient/force calculations Turbomole "
                    "will default to the highest ES, which may lead to suprising results. "
                    "If gradients/forces for certain ES are desired please add a valid "
                    "$exopt datagroup to the control file!"
                )
        elif ("$ricc2" in text) and ("$excitations" in text):
            self.ricc2 = True
            self.ricc2_opt = "geoopt" in text
            second_cmd = "ricc2"
            self.prepare_td(text)
            self.root = self.get_ricc2_root(text)
            try:
                frozen_mos = int(re.search(r"implicit core=\s*(\d+)", text)[1])
            except TypeError:
                frozen_mos = 0
            self.frozen_mos = frozen_mos
            self.log(f"Found {self.frozen_mos} frozen orbitals.")
        if self.track:
            assert (
                self.td or self.ricc2
            ), "track=True can only be used in connection with excited state calculations."
        # Right now this can't handle a root flip from some excited state
        # to the ground state ... Then we would need grad/rdgrad again,
        # instead of egrad.
        self.scf_cmd = scf_cmd
        self.second_cmd = second_cmd
        self.numforce_cmd = f"{self.second_cmd}; NumForce -central -d 0.005"
        if self.ri:
            self.numforce_cmd += " -ri"

        # Setup several cmds, depending on the calc type
        def get_cmd(cmd):
            return ";".join((self.scf_cmd, cmd))

        if self.td:
            self.energy_cmd = get_cmd("escf")
            self.forces_cmd = get_cmd("egrad")
            self.hessian_cmd = "not_yet_implemented"
        elif self.ricc2:
            ricc2_cmd = get_cmd("ricc2")
            self.energy_cmd = ricc2_cmd
            self.forces_cmd = ricc2_cmd
            self.hessian_cmd = "not_yet_implemented"
        else:
            self.energy_cmd = self.scf_cmd
            self.forces_cmd = get_cmd(second_cmd)
            self.hessian_cmd = (
                get_cmd(self.numforce_cmd) if self.numfreq else get_cmd("aoforce")
            )
        self.log("Prepared commands:")
        self.log("\tEnergy cmd: " + self.energy_cmd)
        self.log("\tForces cmd: " + self.forces_cmd)
        self.log("\tHessian cmd: " + self.hessian_cmd)

        if (self.td or self.ricc2) and (self.root is None):
            warnings.warn(
                "No root set! Either include '$exopt' for TDA/TDDFT or "
                "'geoopt' for ricc2 in the control or supply a value for 'root'! "
            )

    def set_occ_and_mo_nums(self, text):
        # Determine number of basis functions
        nbf_re = r"nbf\(AO\)=(\d+)"
        nbf = int(re.search(nbf_re, text)[1])

        self.occ_mos = None
        self.virt_mos = None

        # Determine number of occupied orbitals
        if not self.uhf:
            # Sometimes Turbomole does NOT use a range, but only a single integer,
            # e.g., for hydrogen (H_2).
            occ_re = r"closed shells\s+(\w)\s*(\d+)-?(\d*)"
            occ_mobj = re.search(occ_re, text)
            occ_mos = occ_mobj[3]
            if occ_mos == "":
                occ_mos = occ_mobj[2]
            self.occ_mos = int(occ_mos)

            self.log(f"Found {self.occ_mos} occupied MOs.")
            # Number of spherical basis functions. May be different from CAO
            # Determine number of virtual orbitals
            self.virt_mos = nbf - self.occ_mos
        else:
            alpha_re = r"alpha shells\s+(\w)\s*\d+-(\d+)"
            alpha_mos = int(re.search(alpha_re, text)[2])
            self.log(f"Found {alpha_mos} occupied alpha MOs.")

            beta_re = r"beta shells\s+(\w)\s*\d+-(\d+)"
            beta_mos = int(re.search(beta_re, text)[2])
            self.log(f"Found {beta_mos} occupied beta MOs.")

    def get_ricc2_root(self, text):
        regex = r"geoopt.+?state=\((.+?)\)"
        mobj = re.search(regex, text)
        if not mobj:
            root = None
        elif mobj[1] == "x":
            root = 0
        else:
            assert mobj[1].startswith("a "), "symmetry isn't supported!"
            root = int(mobj[1][-1])
        return root

    def prepare_td(self, text):
        self.log("Preparing for excited state (gradient) calculations")
        self.td_vec_fn = None
        self.ci_coeffs = None
        self.mo_inds = None

    def prepare_point_charges(self, point_charges):
        """$point_charges
        <x> <y> <z> <q>
        """
        lines = [f"{x:.12} {y:.12f} {z:.12f} {q:.12f}" for x, y, z, q in point_charges]

        return "$point_charges\n\t" + "\n".join(lines)

    def prepare_input(self, atoms, coords, calc_type, point_charges=None):
        """To rectify this we have to construct the basecmd
        dynamically and construct it ad hoc. We could set a RI flag
        in the beginning and select the correct scf binary here from
        it. Then we select the following binary on demand, e.g. aoforce
        or rdgrad or egrad etc."""

        valid_calc_types = ("energy", "force", "double_mol", "noparse", "hessian")
        if calc_type not in valid_calc_types:
            raise Exception(
                f"Invalid calc_type '{calc_type}'! Supported "
                f"calc_types are '{valid_calc_types}'."
            )

        path = self.prepare_path(use_in_run=True)
        if calc_type == "double_mol":
            copy_from = self.double_mol_path
        else:
            copy_from = self.control_path
        # Copy everything from the reference control_dir into this path
        # Use self.control_path for all calculations except the double
        # molecule calculation.
        """Maybe we shouldn't copy everything because it may give convergence
        problems? Right now we use the initial MO guess generated in the
        reference path for all images along the path."""
        for glob in copy_from.glob("./*"):
            shutil.copy(glob, path)
        xyz_str = self.prepare_xyz_string(atoms, coords)
        with open(path / "input.xyz", "w") as handle:
            handle.write(xyz_str)
        # Write coordinates
        coord_str = self.prepare_turbo_coords(atoms, coords)
        coord_fn = path / "coord"
        with open(coord_fn, "w") as handle:
            handle.write(coord_str)
        # Copy MO coefficients from previous cycle with this calculator
        # if present.
        if self.mos:
            shutil.copy(self.mos, path / "mos")
            self.log(f"Using {self.mos} as MO guess.")
        elif self.alpha and self.beta:
            shutil.copy(self.alpha, path / "alpha")
            shutil.copy(self.beta, path / "beta")
            self.log(f"Using {self.alpha} and {self.beta} as MO guesses.")
        if self.td_vec_fn:
            # The suffix contains the true name with a leading
            # dot, that we drop.
            td_vec_fn = self.td_vec_fn.suffix[1:]
            shutil.copy(self.td_vec_fn, path / td_vec_fn)
            self.log(f"Using '{self.td_vec_fn}' as escf guess.")

        # Set memory
        self.sub_control(r"\$maxcor.+", f"$maxcor {self.mem} MiB per_core")

        if self.cosmo_kwargs is not None:
            cosmo_data_groups = get_cosmo_data_groups(**self.cosmo_kwargs, atoms=atoms)
            cosmo_rendered = render_data_groups(cosmo_data_groups)
            self.sub_control(r"\$end", cosmo_rendered + "\n$end")

        root_log_msg = f"with current root information: {self.root}"
        if self.root and self.td:
            repl = f"$exopt {self.root}"
            self.sub_control(r"\$exopt\s*(\d+)", f"$exopt {self.root}", root_log_msg)
            self.log(f"Using '{repl}'")

        if self.root and self.ricc2:
            repl = f"state=(a {self.root})"
            self.sub_control(
                r"state=\(a\s+(?P<state>\d+)\)", f"state=(a {self.root})", root_log_msg
            )
            self.log(f"Using '{repl}' for geoopt.")

        if point_charges is not None:
            charge_num = len(point_charges)
            pc_str = self.prepare_point_charges(point_charges)
            self.sub_control(
                r"\$end", pc_str + "\n$end", f"appended {charge_num} point charges"
            )
            # Activate calculation of gradients on point charges
            self.sub_control(r"\$drvopt", "$drvopt\npoint charges\n")
            # Write point charge gradients to file
            self.sub_control(
                r"\$end", "$point_charge_gradients file=pc_gradients\n$end"
            )

        if calc_type == "hessian":
            self.append_control("$noproj\n$nprhessian file=nprhessian")

    def sub_control(self, pattern, repl, log_msg="", **kwargs):
        path = self.path_already_prepared
        assert path
        self.log(f"Updating control file in '{path}' {log_msg}")
        control_path = path / "control"
        with open(control_path) as handle:
            text = handle.read()
        # text = re.sub(pattern, repl, text, **kwargs)
        text, nsubs = re.subn(pattern, repl, text, **kwargs)
        with open(control_path, "w") as handle:
            handle.write(text)

    def append_control(self, to_append, log_msg="", **kwargs):
        self.sub_control(r"\$end", f"{to_append}\n$end", log_msg, **kwargs)

    def get_pal_env(self):
        env_copy = os.environ.copy()
        env_copy["PARA_ARCH"] = "SMP"
        env_copy["PARNODES"] = str(self.pal)
        env_copy["SMPCPUS"] = str(self.pal)

        return env_copy

    def store_and_track(self, results, func, atoms, coords, **prepare_kwargs):
        # Increase counter manually as it wasn't increased in Calculator.run()
        self.calc_counter += 1
        if self.track:
            prev_run_path = self.last_run_path
            self.store_overlap_data(atoms, coords)
            # Redo the calculation with the updated root
            if self.track_root():
                results = func(atoms, coords, **prepare_kwargs)
            self.last_run_path = prev_run_path
        try:
            shutil.rmtree(self.last_run_path)
        except FileNotFoundError:
            self.log(f"'{self.last_run_path}' has already been deleted!")
        results["all_energies"] = self.parse_all_energies()
        return results

    def get_energy(self, atoms, coords, **prepare_kwargs):
        self.prepare_input(atoms, coords, "energy", **prepare_kwargs)
        kwargs = {
            "calc": "energy",
            "shell": True,
            "hold": self.track,
            "env": self.get_pal_env(),
            "cmd": self.energy_cmd,
        }
        results = self.run(None, **kwargs)
        results = self.store_and_track(
            results, self.get_energy, atoms, coords, **prepare_kwargs
        )
        return results

    def get_all_energies(self, atoms, coords, **prepare_kwargs):
        results = self.get_energy(atoms, coords, **prepare_kwargs)

        with open(self.out) as handle:
            text = handle.read()

        ((nfocc_a, nfvirt_a), (nfocc_b, nfvirt_b)), restricted = parse_frozen_nmos(text)

        results["td_1tdms"] = parse_ci_coeffs(
            self.td_vec_fn,
            nfocc_a,
            nfvirt_a,
            nfocc_b,
            nfvirt_b,
            restricted_same_ab=True,
        )
        return results

    def get_forces(self, atoms, coords, cmd=None, **prepare_kwargs):
        self.prepare_input(atoms, coords, "force", **prepare_kwargs)

        if cmd is None:
            cmd = self.forces_cmd

        kwargs = {
            "calc": "force",
            "shell": True,  # To allow chained commands like 'ridft; rdgrad'
            "hold": self.track,  # Keep the files for WFOverlap
            "env": self.get_pal_env(),
            "cmd": cmd,
        }
        # Use inp=None because we don't use any dedicated input besides
        # the previously prepared control file and the current coords.
        results = self.run(None, **kwargs)
        results = self.store_and_track(
            results, self.get_forces, atoms, coords, **prepare_kwargs
        )
        return results

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        if self.td or self.ricc2:
            raise Exception("ricc2 or TD-DFT/TDA hessian not yet supported!")

        self.prepare_input(atoms, coords, "hessian", **prepare_kwargs)
        kwargs = {
            "calc": "hessian",
            "shell": True,  # To allow chained commands like 'ridft; rdgrad'
            "hold": self.track,
            "env": self.get_pal_env(),
            "cmd": self.hessian_cmd,
        }
        results = self.run(None, **kwargs)
        results = self.store_and_track(
            results, self.get_hessian, atoms, coords, **prepare_kwargs
        )
        return results

    def get_stored_wavefunction(self, **kwargs):
        return self.load_wavefunction_from_file(self.mwfn_wf, **kwargs)

    def get_wavefunction(self, atoms, coords, **prepare_kwargs):
        with TurbomoleGroundStateContext(self):
            results = self.get_energy(atoms, coords, **prepare_kwargs)
            results["wavefunction"] = self.get_stored_wavefunction()
        return results

    def get_relaxed_density(self, atoms, coords, root, **prepare_kwargs):
        root_bak = self.root
        self.root = root
        results = self.get_forces(atoms, coords, **prepare_kwargs)
        self.root = root_bak
        if self.uhf:
            with open(self.alpha) as handle:
                text = handle.read()
            Ca = parse_turbo_mos(text)
            with open(self.beta) as handle:
                text = handle.read()
            Cb = parse_turbo_mos(text)
        else:
            with open(self.mos) as handle:
                text = handle.read()
            Ca = parse_turbo_mos(text)
            Cb = Ca.copy()
        Pa, Pb = get_density_matrices_for_root(
            self.out, self.td_vec_fn, root, rlx_vec_fn=self.rlx_vec_fn, Ca=Ca, Cb=Cb
        )
        density = np.stack((Pa, Pb), axis=0)
        results["density"] = density
        return results

    def run_calculation(self, atoms, coords, **prepare_kwargs):
        return self.get_energy(atoms, coords, **prepare_kwargs)

    def run_double_mol_calculation(self, atoms, coords1, coords2):
        if not self.double_mol_path:
            self.log(
                "Skipping double molecule calculations as double mol "
                "path is not specified.!"
            )
            return None
        self.log("Running double molecule calculation")
        double_atoms = atoms + atoms
        double_coords = np.hstack((coords1, coords2))
        self.prepare_input(double_atoms, double_coords, "double_mol")
        kwargs = {
            "calc": "double_mol",
            "shell": True,
            "keep": False,
            "hold": True,
            "cmd": self.scf_cmd,
            "env": self.get_pal_env(),
        }
        results = self.run(None, **kwargs)
        return results

    def parse_double_mol(self, path):
        """Parse a double molecule overlap matrix from Turbomole output
        to be used with WFOWrapper."""
        with open(path / self.out_fn) as handle:
            text = handle.read()
        regex = r"OVERLAP\(SAO\)\s+-+([\d\.E\-\s*\+]+)\s+-+"
        ovlp_str = re.search(regex, text)[1]
        ovlp = np.array(ovlp_str.strip().split(), dtype=np.float64)
        mo_num = self.occ_mos + self.virt_mos
        double_mo_num = 2 * mo_num
        full_ovlp = np.zeros((double_mo_num, double_mo_num))
        full_ovlp[np.tril_indices(double_mo_num)] = ovlp
        double_mol_S = full_ovlp[mo_num:, :mo_num]
        return double_mol_S

    def parse_mos(self):
        pass

    def parse_energy(self, path):
        with open(path / self.out_fn) as handle:
            text = handle.read()
        en_regex = re.compile(r"Total energy\s*:?\s*=?\s*([\d\-\.]+)", re.IGNORECASE)
        tot_ens = en_regex.findall(text)

        # Only modify energy when self.root is set; otherwise stick with the GS energy.
        if self.td and self.root is not None:
            # Drop ground state energy that is repeated. That is why we don't subtract
            # 1 from self.root.
            tot_en = tot_ens[1:][self.root]
        elif self.ricc2 and self.ricc2_opt:
            results = parse_turbo_gradient(path)
            tot_en = results["energy"]
        elif self.ricc2 and not self.ricc2_opt:
            raise Exception("Implement me!")
        else:
            tot_en = tot_ens[0]

        tot_en = float(tot_en)
        return {
            "energy": tot_en,
        }

    @staticmethod
    @file_or_str(".sing_a", ".ciss_a")
    def parse_tddft_tden(text):
        eigval_re = re.compile(r"(\d+)\s+eigenvalue\s+=\s+([\d\.\-D\+]+)")
        eigvals = eigval_re.findall(text)
        state_inds, exc_ens = zip(*eigvals)
        exc_ens = [exc_en.replace("D", "E") for exc_en in exc_ens]
        return np.array(exc_ens, dtype=float)

    def parse_all_energies(self):
        # Parse eigenvectors from escf/egrad calculation
        gs_energy = self.parse_gs_energy()
        # if self.root and self.second_cmd in ("escf", "egrad"):
        if self.second_cmd in ("escf", "egrad") and hasattr(self, "exspectrum"):
            # I don't know why, but sometimes sing_a contains wrong excitation energies...
            # exc_energies = Turbomole.parse_tddft_tden(self.td_vec_fn)
            roots = parse_exspectrum(self.exspectrum)
            exc_energies = np.array([root.exc_energy for root in roots])
        # Parse eigenvectors from ricc2 calculation
        elif self.second_cmd == "ricc2":
            with open(self.exstates) as handle:
                exstates_text = handle.read()
            exc_energies_by_model = parse_turbo_exstates(exstates_text)
            # Drop CCS and take energies from whatever model was used
            exc_energies = [
                (model, exc_ens)
                for model, exc_ens in exc_energies_by_model.items()
                if model != "CCS"
            ]
            assert len(exc_energies) == 1
            _, exc_energies = exc_energies[0]

        else:
            exc_energies = np.array(list())

        if exc_energies.size == 0:
            all_energies = np.array((gs_energy,))
        else:
            all_energies = np.full(len(exc_energies) + 1, gs_energy)
            all_energies[1:] += exc_energies

        return all_energies

    def parse_ci_coeffs(self):
        if self.second_cmd in ("escf", "egrad"):
            with open(self.td_vec_fn) as handle:
                text = handle.read()
            ci_coeffs = self.parse_td_vectors(text)
            ci_coeffs = [cc["vector"] for cc in ci_coeffs]
        # Parse eigenvectors from ricc2 calculation
        elif self.second_cmd == "ricc2":
            ci_coeffs = [self.parse_cc2_vectors(ccre) for ccre in self.ccres]

        ci_coeffs = np.array(ci_coeffs)
        states = ci_coeffs.shape[0]
        tden_size = self.occ_mos * self.virt_mos
        if ci_coeffs.shape[1] == (2 * tden_size):
            self.log("TDDFT calculation with X and Y vectors present. ")
            XpY = ci_coeffs[:, :tden_size]
            XmY = ci_coeffs[:, tden_size:]
            X = (XpY + XmY) / 2.0
            Y = XpY - X
        else:
            X = ci_coeffs
            Y = np.zeros_like(X)
        tden_shape = (states, self.occ_mos, self.virt_mos)
        X = X.reshape(tden_shape)
        Y = Y.reshape(tden_shape)

        return X, Y

    def parse_force(self, path):
        results = parse_turbo_gradient(path)
        return results

    def parse_hessian(self, path, fn=None):
        if fn is None:
            if self.numfreq:
                fn = path / "numforce" / "nprhessian"
            else:
                fn = path / "nprhessian"

        with open(fn) as handle:
            text = handle.read()

        split = text.strip().split()
        assert split[0] == "$nprhessian"
        assert split[-1] == "$end"

        def is_float(str_):
            return "." in str_

        hess_items = [item for item in split if is_float(item)]
        coord_num = int(sqrt(len(hess_items)))
        assert coord_num**2 == len(hess_items)
        hessian = np.array(hess_items, dtype=float).reshape(-1, coord_num)

        energy = self.parse_energy(path)["energy"]

        results = {
            "energy": energy,
            "hessian": hessian,
        }
        return results

    @staticmethod
    def parse_td_vectors(text):
        """For TDA calculations only the X vector is present in the
        ciss_a/etc. file. In TDDFT calculations there are twise as
        much items compared with TDA. The first half corresponds to
        (X+Y) and the second half to (X-Y). X can be calculated as
        X = ((X+Y)+(X-Y))/2. Y is then given as Y = (X+Y)-X. The
        normalization can then by checked as
            np.concatenate((X, Y)).dot(np.concatenate((X, -Y)))
        and should be 1."""

        def to_float(s, loc, toks):
            match = toks[0].replace("D", "E")
            return float(match)

        float_ = pp.Word(pp.nums + ".-D+").setParseAction(to_float)
        integer = pp.Word(pp.nums).setParseAction(lambda t: int(t[0]))
        float_chrs = pp.nums + "D.+-"
        float_20 = pp.Word(float_chrs, exact=20).setParseAction(to_float)

        word = pp.Word(pp.alphanums)
        line_word = word().setWhitespaceChars(" ")
        line = pp.Group(pp.Optional(word) + line_word[...])("title")
        title = pp.Literal("$title") + pp.Optional(line)
        symmetry = pp.Literal("$symmetry") + pp.Word(pp.alphanums).setResultsName(
            "symmetry"
        )
        tensor_dim = pp.Literal("$tensor space dimension") + integer.setResultsName(
            "tensor_dim"
        )
        scfinstab = pp.Literal("$scfinstab") + pp.Word(pp.alphanums).setResultsName(
            "scfinstab"
        )
        subspace_dim = pp.Literal(
            "$current subspace dimension"
        ) + integer.setResultsName("subspace_dim")
        converged = pp.Literal("$current iteration converged")
        eigenpairs = pp.Literal("$eigenpairs")
        eigenpair = pp.Group(
            integer.setResultsName("state")
            + pp.Literal("eigenvalue =")
            + float_.setResultsName("eigenvalue")
            + pp.Group(pp.OneOrMore(float_20)).setResultsName("vector")
        )
        end = pp.Literal("$end")

        parser = (
            title
            + symmetry
            + tensor_dim
            + scfinstab
            + subspace_dim
            + converged
            + eigenpairs
            + pp.OneOrMore(eigenpair).setResultsName("eigenpairs")
            + end
        )
        result = parser.parseString(text)
        states = result["subspace_dim"]
        eigenpairs = result["eigenpairs"]
        eigenpair_list = [eigenpairs[i].asDict() for i in range(states)]
        return eigenpair_list

    def parse_cc2_vectors(self, ccre):
        with open(ccre) as handle:
            text = handle.read()
        coeffs = parse_turbo_ccre0_ascii(text)
        coeffs = coeffs.reshape(-1, self.virt_mos)

        eigenpairs_full = np.zeros((self.occ_mos, self.virt_mos))
        eigenpairs_full[self.frozen_mos :, :] = coeffs
        """
        from_inds, to_inds = np.where(np.abs(eigenpairs_full) > 0.1)
        for i, (from_, to_) in enumerate(zip(from_inds, to_inds)):
            sq = eigenpairs_full[from_, to_]**2
            print(f"{from_+1:02d} -> {to_+self.occ_mos+1:02d}: {sq:.2%}")
        print()
        """

        return eigenpairs_full

    def parse_gs_energy(self):
        """Several places are possible:
        $subenergy from control file
        total energy from turbomole.out
        Final MP2 energy from turbomole.out with ADC(2)
        Final CC2 energy from turbomole.out with CC(2)
        """
        float_re = r"([\d\-\.E]+)"
        regexs = [
            # CC2 ground state energy
            ("out", r"Final CC2 energy\s*:\s*" + float_re, 0),
            # ADC(2) ground state energy
            ("out", r"Final MP2 energy\s*:\s*" + float_re, 0),
            ("control", r"\$subenergy.*$\s*" + float_re, re.MULTILINE),
            # DSCF ground state energy
            ("out", r"total energy\s*=\s*" + float_re, 0),
            # From egrad when a rootflip occured. Then only the excited
            # state calculation will be redone and the ground state calculation
            # won't be present in the out-file.
            ("out", r"Ground state\s*?Total energy:\s+" + float_re, re.MULTILINE),
        ]
        for file_attr, regex, flag in regexs:
            regex_ = re.compile(regex, flags=flag)
            with open(getattr(self, file_attr)) as handle:
                text = handle.read()
            mobj = regex_.search(text)
            try:
                gs_energy = float(mobj[1])
                self.log(
                    f"Parsed ground state energy from '{file_attr}' using "
                    f"regex '{regex[:11]}'."
                )
                return gs_energy
            except TypeError:
                continue
        raise Exception("Couldn't parse ground state energy!")

    """
    def prepare_overlap_data(self, path):
        all_energies = self.parse_all_energies()
        X, Y = self.parse_ci_coeffs()
        self.log(f"Reading MO coefficients from '{self.mos}'.")
        with open(self.mos) as handle:
            text = handle.read()
        C = parse_turbo_mos(text)
        self.log(f"Reading electronic energies from '{self.out}'.")
        return C, X, Y, all_energies
    """

    def prepare_overlap_data(self, path):
        base = self.out.parent / self.out.stem
        Ca, Cb, Xa, Ya, Xb, Yb = get_overlap_data_from_base_name(
            base, ci_suffix=self.td_vec_fn.suffix
        )
        all_energies = self.parse_all_energies()
        return Ca, Xa, Ya, Cb, Xb, Yb, all_energies

    def run_after(self, path):
        # Convert binary CCRE0 files to ASCII for easier parsing
        for ccre in path.glob("CCRE0-*"):
            cmd = f"ricctools -dump {ccre.name}".split()
            result = subprocess.Popen(
                cmd, cwd=path, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            result.wait()
            # ricctools seem to crash sometimes, even though the respective ASCII
            # ccre-files are generated.
            subprocess.run(
                "actual -r".split(),
                cwd=path,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

        if self.wavefunction_dump or self.td:
            self.make_molden(path)
        # With ricc2 we probably have a frozen core that we have to disable
        # temporarily before creating the molden file. Afterwards we restore
        # the original control file with the frozen core.
        elif self.ricc2:
            # Backup original control file
            ctrl_backup = path / "control.backup"
            shutil.copy(path / "control", ctrl_backup)
            # We have to remove line with implicit core in the control file
            with open(path / "control") as handle:
                text = handle.read()
            lines = text.split("\n")
            lines = [l for l in lines if "implicit core" not in l]
            with open(path / "control", "w") as handle:
                handle.write("\n".join(lines))
            self.make_molden(path)
            # Restore control backup
            shutil.copy(ctrl_backup, path / "control")

    @staticmethod
    def make_molden(path):
        cmd = "tm2molden norm".split()
        fn = "wavefunction.molden"
        stdin = f"""{fn}

        """
        res = subprocess.Popen(
            cmd,
            cwd=path,
            universal_newlines=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = res.communicate(stdin)
        res.terminate()

    def get_chkfiles(self):
        if self.uhf:
            chkfiles = {
                "alpha": self.alpha,
                "beta": self.beta,
            }
        else:
            chkfiles = {
                "mos": self.mos,
            }
        return chkfiles

    def set_chkfiles(self, chkfiles):
        try:
            if self.uhf:
                alpha = chkfiles["alpha"]
                beta = chkfiles["beta"]
                self.alpha = alpha
                self.beta = beta
                self.log(f"Set chkfile '{alpha}' and '{beta}' as alpha and beta.")
            else:
                mos = chkfiles["mos"]
                self.mos = mos
                self.log(f"Set chkfile '{mos}' as mos.")
        except KeyError:
            self.log("Found no chkfile information in chkfiles!")

    def __str__(self):
        return "Turbomole calculator"
