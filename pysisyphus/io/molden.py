import re

import numpy as np
import pyparsing as pp
from scipy.special import gamma

from pysisyphus.elem_data import ATOMIC_NUMBERS, nuc_charges_for_atoms
from pysisyphus.Geometry import Geometry
from pysisyphus.io.xyz import parse_xyz
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.wavefunction import (
    get_l,
    MoldenShells,
    ORCAMoldenShells,
    Shell,
    Wavefunction,
)
from pysisyphus.wavefunction.helpers import BFType


@file_or_str(".molden", ".input")
def parse_molden_atoms_coords(molden_str):
    _, geometries = re.split(r"\[GEOMETRIES\] \(XYZ\)", molden_str)
    atoms_coords, comments = parse_xyz(geometries, with_comment=True)

    return atoms_coords, comments


def geoms_from_molden(fn, **kwargs):
    atoms_coords, comments = parse_molden_atoms_coords(fn)
    geoms = [
        Geometry(atoms, coords.flatten(), comment=comment, **kwargs)
        for (atoms, coords), comment in zip(atoms_coords, comments)
    ]
    return geoms


def parse_mo(mo_lines):
    # sym, ene_key, energy, spin_key, spin, occup_key, occup, *rest = mo_lines
    sym, ene_line, spin_line, occup_line, *rest = mo_lines
    ene_key, energy = ene_line.split("=")
    spin_key, spin = spin_line.split("=")
    occup_key, occup = occup_line.split("=")
    assert (
        (ene_key.strip() == "ene")
        and (spin_key.strip() == "spin")
        and (occup_key.strip() == "occup")
    )
    prev_ind = 0
    coeffs = list()
    for line in rest:
        ind, coeff = line.split()
        ind = int(ind)
        assert ind == prev_ind + 1
        prev_ind = ind
        coeffs.append(coeff)
    return {
        "sym": sym.strip(),
        "ene": float(energy),
        "spin": spin.strip(),
        "occup": float(occup),
        "coeffs": np.array(coeffs),
    }


@file_or_str(".molden", ".input")
def parse_molden(text, with_mos=True):
    int_ = pp.common.integer
    real = pp.common.real
    sci_real = pp.common.sci_real

    def get_section(name):
        return pp.CaselessLiteral(f"[{name}]")

    molden_format = get_section("Molden Format")
    n_atoms = get_section("N_Atoms") + int_
    atoms_header = get_section("Atoms")
    title = get_section("Title") + pp.Optional(
        pp.ZeroOrMore(~atoms_header + pp.Word(pp.printables))
    )

    alpha_re = re.compile("([a-zA-Z]+)")

    def drop_number(toks):
        """OpenMolcas writes atoms like H1, C2, ...
        We don't want the number."""
        return alpha_re.match(toks[0]).group(1)

    # [Atoms]
    # element_name number atomic_number x y z
    unit = pp.one_of(("AU Angs (AU) (Angs)"), caseless=True)
    atom_symbol = pp.Word(pp.alphanums).set_parse_action(drop_number)
    xyz = pp.Group(sci_real + sci_real + sci_real)
    atom_line = pp.Group(
        atom_symbol.set_results_name("symbol")
        + int_.set_results_name("number")
        + int_.set_results_name("atomic_number")
        + xyz.set_results_name("xyz")
    )
    atoms = (
        atoms_header
        + pp.Optional(unit).set_results_name("unit")
        + pp.OneOrMore(atom_line).set_results_name("atoms")
    )

    # [Charge], written by OpenMOLCAS
    charge = get_section("Charge")
    charge_kind = pp.one_of(("(Mulliken)",), caseless=True)
    charges = charge + charge_kind + pp.OneOrMore(sci_real)

    def pgto_exps_coeffs(s, loc, toks):
        data = toks.asList()
        exps, coeffs = np.array(data, dtype=float).reshape(-1, 2).T
        pr = pp.ParseResults.from_dict(
            {
                "exponents": exps,
                "coeffs": coeffs,
            }
        )
        # Somehow this function will result in exponents and coeffs being both
        # stored in lists of length 1.
        return pr

    # [GTO]
    ang_mom = pp.one_of("s p sp d f g h i", caseless=True)
    pgto = sci_real + sci_real
    shell = pp.Group(
        ang_mom.set_results_name("ang_mom")
        + int_.set_results_name("contr_depth")
        # Usually 1.0, can be left out
        + pp.Optional(real.set_whitespace_chars(" ").suppress() + pp.LineEnd())
        + pp.OneOrMore(pgto).set_parse_action(pgto_exps_coeffs)
    )
    atom_gtos = pp.Group(
        int_.set_results_name("number")  # Center
        + pp.Optional(pp.Literal("0"))
        + pp.OneOrMore(shell).set_results_name("shells")
    )
    gto = (
        pp.CaselessLiteral("[GTO]")
        + pp.Optional(unit)
        + pp.Group(pp.OneOrMore(atom_gtos)).set_results_name("gto")
    )

    # [5D] [7F] [9G] etc.
    ang_mom_flags = pp.ZeroOrMore(pp.one_of("[5D] [7F] [9G]", caseless=True))

    # Pyparsing code to parse MOs ... slow. Regex-based code is much faster.
    # [MO]
    # sym_label = pp.Word(pp.printables)
    # spin = pp.one_of("Alpha Beta", caseless=True)
    # mo_coeff = pp.Suppress(int_) + real
    # mo = pp.Group(
    # pp.CaselessLiteral("Sym=").suppress()
    # + sym_label.set_results_name("sym")
    # + pp.CaselessLiteral("Ene=").suppress()
    # + sci_real.set_results_name("ene")
    # + pp.CaselessLiteral("Spin=").suppress()
    # + spin.set_results_name("spin")
    # + pp.CaselessLiteral("Occup=").suppress()
    # + real.set_results_name("occup")
    # + pp.Group(pp.OneOrMore(mo_coeff)).set_results_name("coeffs")
    # )
    # mos = pp.CaselessLiteral("[MO]") + pp.OneOrMore(mo).set_results_name("mos")
    mos = pp.CaselessLiteral("[MO]") + pp.OneOrMore(
        pp.Regex(r"[^\[]+")  # Match everything until the next [ is encountered
    ).set_results_name("mos")

    # Actual parser
    parser = molden_format + pp.Each(
        [
            pp.Optional(expr)
            for expr in (n_atoms, title, atoms, charges, gto, ang_mom_flags)
        ]
    )
    if with_mos:
        parser += mos

    res = parser.parse_string(text)
    as_dict = res.asDict()
    if with_mos:
        mos_raw = as_dict["mos"][0]

        by_mo = [bm.strip().split("\n") for bm in mos_raw.lower().strip().split("sym=")]
        assert by_mo[0] == [""], by_mo[0]
        by_mo = by_mo[1:]

        as_dict["mos"] = [parse_mo(mo_lines) for mo_lines in by_mo]
    return as_dict


def get_xtb_nuc_charges(atoms, as_ecp_electrons=False):
    """Modified nuclear charges w/o core electrons.

    Adapated from Multiwfn 3.8
    """

    org_nuc_charges = nuc_charges_for_atoms(atoms)
    ecp_electrons = np.zeros_like(org_nuc_charges, dtype=int)
    for i, Z in enumerate(org_nuc_charges):
        # H, He
        if 1 <= Z <= 2:
            ecpel = 0
        # Li - Ne
        elif 3 <= Z <= 10:
            ecpel = 2
        # Na - Ar
        elif 11 <= Z <= 18:
            ecpel = 10
        # K - Cu
        elif 19 <= Z <= 29:
            ecpel = 18
        # Zn - Kr
        elif 30 <= Z <= 36:
            ecpel = 28
        # Rb - Ag
        elif 37 <= Z <= 47:
            ecpel = 36
        # Cd - Xe
        elif 48 <= Z <= 54:
            ecpel = 46
        # Cs - Ba
        elif 55 <= Z <= 56:
            ecpel = 54
        # La - Lu, crazy
        elif 57 <= Z <= 71:
            ecpel = Z - 3  # Resulting charge is always 3
        # Hf - Au
        elif 72 <= Z <= 79:
            ecpel = 68
        # Hg - Rn
        elif 80 <= Z <= 86:
            ecpel = 78
        ecp_electrons[i] = ecpel
    if as_ecp_electrons:
        result = ecp_electrons
    else:
        result = org_nuc_charges - ecp_electrons
    return result


def parse_molden_atoms(data):
    atoms = list()
    coords = list()
    nuc_charges = list()
    for atom in data["atoms"]:
        atoms.append(atom["symbol"])
        coords.extend(atom["xyz"])
        Z = atom["atomic_number"]
        nuc_charges.append(Z)
    atoms = tuple(atoms)
    coords = np.array(coords, dtype=float).flatten()
    nuc_charges = np.array(nuc_charges, dtype=int)
    return atoms, coords, nuc_charges


@file_or_str(".molden", ".input")
def shells_from_molden(text):
    data = parse_molden(text, with_mos=False)

    atoms, coords, _ = parse_molden_atoms(data)
    atomic_numbers = [ATOMIC_NUMBERS[atom.lower()] for atom in atoms]
    coords3d = coords.reshape(-1, 3)

    _shells = list()
    for atom_gtos, center, atomic_num in zip(data["gto"], coords3d, atomic_numbers):
        center_ind = atom_gtos["number"] - 1
        for shell in atom_gtos["shells"]:
            L = shell["ang_mom"]
            exps = shell["exponents"][0]
            coeffs = shell["coeffs"][0]
            shell = Shell(
                L=L,
                center=center,
                coeffs=coeffs,
                exps=exps,
                center_ind=center_ind,
                atomic_num=atomic_num,
            )
            _shells.append(shell)

    shells = MoldenShells(_shells)
    return shells


def radial_integral(l, exponent):
    """
    Integrates
        (r r**l * exp(-exponent * r**2))**2 dr from r=0 to r=oo
    as described in the SI of the JANPA paper [1] (see top of page 8,
    second integral in the square root.

    In my opinion, the integrals lacks a factor 'r'. Below, some sympy code
    can be found to solve this integral (including 1*r).

    import sympy as sym
    r, z = sym.symbols("r z", positive=True)
    l = sym.symbols("l", integer=True, positive=True)
    sym.integrate((r * r**l * sym.exp(-z*r**2))**2, (r, 0, sym.oo))

    The 'solved' integral on page 8 is correct again.

     ∞
     ⌠
     ⎮              2
     ⎮  2  2⋅l  -2⋅r ⋅z
     ⎮ r ⋅r   ⋅ℯ        dr = (2*z)**(-l - 1/2)*gamma(l + 3/2)/(4*z)
     ⌡
     0
    """
    return (2 * exponent) ** (-l - 1 / 2) * gamma(l + 3 / 2) / (4 * exponent)


@file_or_str(".molden", ".input")
def shells_from_orca_molden(text):
    molden_shells = shells_from_molden(text)

    dividers = {
        0: 1,
        1: 1,
        2: 3,
        3: 15,
        4: 35,
    }

    def fix_contr_coeffs(l, coeffs, exponents):
        """Fix contraction coefficients. Based on equations found in the SI
        of the JANPA paper [1]."""
        l = get_l(l)
        divider = dividers[l]
        rad_ints = radial_integral(l, exponents)
        norms2 = rad_ints * 4 * np.pi / (2 * l + 1) / divider
        norms = np.sqrt(norms2)
        normed_coeffs = coeffs * norms
        return normed_coeffs

    _shells = list()
    for shell in molden_shells.shells:
        L, center, _, exps, *_ = shell.as_tuple()
        fixed_coeffs = fix_contr_coeffs(L, shell.coeffs_org, exps)
        fixed_shell = Shell(
            L=L,
            center=center,
            coeffs=fixed_coeffs,
            exps=exps,
            center_ind=shell.center_ind,
            atomic_num=shell.atomic_num,
        )
        _shells.append(fixed_shell)
    shells = ORCAMoldenShells(_shells)
    return shells


@file_or_str(".molden", ".input")
def wavefunction_from_molden(
    text,
    charge=None,
    orca_contr_renorm=False,
    xtb_nuc_charges: bool = False,
    **wf_kwargs,
):
    """Construct Wavefunction object from .molden file.

    orca_contr_renorm fixes contraction coefficients, as required for
    ORCA and XTB.
    """
    data = parse_molden(text)
    atoms, coords, nuc_charges = parse_molden_atoms(data)
    nuc_charge = sum(nuc_charges)
    # XTB only describes valence electrons. The missing electrons will be
    # treated as ECP electrons.
    if xtb_nuc_charges:
        ecp_electrons = get_xtb_nuc_charges(atoms, as_ecp_electrons=True)
        wf_kwargs["ecp_electrons"] = ecp_electrons
        nuc_charge = nuc_charge - sum(ecp_electrons)

    spins = list()
    occ_a = 0.0
    occ_b = 0.0
    occ_a_unpaired = 0.0
    occ_b_unpaired = 0.0
    Ca = list()
    Cb = list()
    for mo in data["mos"]:
        spin = mo["spin"].lower()
        occ = mo["occup"]
        assert occ >= 0.0
        is_unpaired = 0.0 < occ < 2.0
        spins.append(spin)
        coeffs = mo["coeffs"]
        if spin == "alpha":
            occ_a += occ
            Ca.append(coeffs)
            if is_unpaired:
                occ_a_unpaired += occ
        elif spin == "beta":
            occ_b += occ
            Cb.append(coeffs)
            if is_unpaired:
                occ_b_unpaired += occ
        else:
            raise Exception(f"Spin can only be 'alpha' or 'beta', but got '{spin}'!")
    assert occ_a.is_integer
    assert occ_b.is_integer
    occ_a = int(occ_a)
    occ_b = int(occ_b)
    occ_a_unpaired = int(occ_a_unpaired)
    occ_b_unpaired = int(occ_b_unpaired)
    unpaired_electrons = occ_a_unpaired or occ_b_unpaired
    unrestricted = set(spins) == {"Alpha", "Beta"} or unpaired_electrons
    restricted = not unrestricted

    # MOs must be in columns
    Ca = np.array(Ca, dtype=float).T
    Cb = np.array(Cb, dtype=float).T
    # Restricted calculation
    if restricted:
        Cb = Ca.copy()
        occ_a = occ_b = occ_a // 2
    # The second case is hit w/ XTB molden files and restriced-open-shell calculatons
    elif unrestricted and Cb.size == 0:
        assert occ_a >= occ_b
        common_electrons = (occ_a - occ_a_unpaired) // 2
        Cb = Ca.copy()
        occ_a -= common_electrons
        occ_b = common_electrons
    C = np.stack((Ca, Cb))

    molden_charge = nuc_charge - (occ_a + occ_b)
    if charge is None:
        charge = molden_charge
    # Multiplicity = (2S + 1), S = number of unpaired elecs. * 0.5
    mult = (occ_a - occ_b) + 1

    if orca_contr_renorm:
        shells_func = shells_from_orca_molden
    else:
        shells_func = shells_from_molden

    shells = shells_func(text)

    return Wavefunction(
        atoms=atoms,
        coords=coords,
        charge=charge,
        mult=mult,
        unrestricted=unrestricted,
        occ=(occ_a, occ_b),
        C=C,
        bf_type=BFType.PURE_SPHERICAL,
        shells=shells,
        **wf_kwargs,
    )
