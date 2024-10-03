import numpy as np
import pyparsing as pp

from pysisyphus.elem_data import ATOMIC_NUMBERS, nuc_charges_for_atoms
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.wavefunction import get_l, Shell, AOMixShells, Wavefunction
from pysisyphus.wavefunction.helpers import BFType

AOMIX_EXTS = (".in", ".aomix")


@file_or_str(*AOMIX_EXTS)
def parse_aomix(text):
    real = pp.common.real
    sreal = pp.common.sci_real
    int_ = pp.common.integer
    header = pp.CaselessLiteral("[AOMix Format]")
    title = pp.CaselessKeyword("[Title]")
    energy = pp.Suppress(pp.CaselessLiteral("[SCF Energy / Hartree]")) + sreal
    xyz = pp.Group(sreal + sreal + sreal)
    atom_line = pp.Group(
        pp.Word(pp.alphas).set_results_name("atom")
        + int_.set_results_name("id")
        + int_.set_results_name("atom_num")
        + xyz.set_results_name("xyz")
    )
    atoms = (
        pp.Suppress(pp.CaselessLiteral("[Atoms]"))
        + pp.Optional(pp.CaselessLiteral("AU")).set_results_name("unit")
        + pp.OneOrMore(atom_line).set_results_name("atom_lines")
    )
    prim_gto = pp.Group(sreal + sreal)
    contr_gto = pp.Group(
        pp.Word(pp.alphas).set_results_name("l")
        + int_.set_results_name("contr_depth")
        + real
        + pp.OneOrMore(prim_gto).set_results_name("prim_gtos")
    )
    center_contr_gto = pp.Group(
        int_.set_results_name("center")
        + int_
        + pp.OneOrMore(contr_gto).set_results_name("contr_gtos")
    )
    gtos = pp.CaselessLiteral("[GTO]") + pp.OneOrMore(
        center_contr_gto
    ).set_results_name("center_gtos")
    ao_label = pp.Word(pp.alphanums)
    mo_line = pp.Suppress(int_ + int_ + pp.Word(pp.alphas) + int_ + ao_label) + sreal
    sym_label = int_ + pp.one_of(("a", "b"), caseless=True)
    mo = pp.Group(
        pp.CaselessLiteral("Sym=")
        + sym_label.set_results_name("sym_label")
        + pp.CaselessLiteral("Ene=")
        + sreal.set_results_name("ene")
        + pp.CaselessLiteral("Spin=")
        + pp.one_of(("Alpha", "Beta"), caseless=True).set_results_name("spin")
        + pp.CaselessLiteral("Occup=")
        + sreal.set_results_name("occup")
        + pp.Group(pp.OneOrMore(mo_line)).set_results_name("coeffs")
    )
    mos = pp.CaselessLiteral("[MO]") + pp.OneOrMore(mo).set_results_name("mos")

    parser = (
        header
        + pp.Optional(title).set_results_name("title")
        + pp.Optional(energy).set_results_name("energy")
        + atoms
        + gtos
        + mos
    )
    result = parser.parseString(text)
    as_dict = result.asDict()

    return as_dict


def shells_from_aomix_dict(aomix_dict, **kwargs):
    atom_lines = aomix_dict["atom_lines"]
    center_gtos = aomix_dict["center_gtos"]
    _shells = list()
    for atom_line, cgtos in zip(atom_lines, center_gtos):
        atomic_num = ATOMIC_NUMBERS[atom_line["atom"].lower()]
        center_ind = cgtos["center"]
        assert atom_line["id"] == center_ind
        center_ind -= 1  # Start from 0 internally
        center = np.array(atom_line["xyz"])
        for cgto in cgtos["contr_gtos"]:
            exps, coeffs = zip(*cgto["prim_gtos"])
            exps = np.array(exps)
            coeffs = np.array(coeffs)
            L = get_l(cgto["l"])
            shell = Shell(
                L=L,
                atomic_num=atomic_num,
                center=center,
                coeffs=coeffs,
                exps=exps,
                center_ind=center_ind,
            )
            _shells.append(shell)
    shells = AOMixShells(_shells, **kwargs)
    return shells


def wavefunction_from_aomix_dict(aomix_dict, **kwargs):
    kwargs = kwargs.copy()
    shell_kwargs = kwargs.pop("shell_kwargs", {})
    shells = shells_from_aomix_dict(aomix_dict, **shell_kwargs)

    atoms = list()
    coords3d = list()
    for atom_line in aomix_dict["atom_lines"]:
        atoms.append(atom_line["atom"].lower())
        coords3d.append(atom_line["xyz"])
    atoms = tuple(atoms)
    coords = np.array(coords3d).flatten()

    occ_a = 0
    occ_b = 0
    unrestricted = False
    C_a = list()
    C_b = list()
    for mo in aomix_dict["mos"]:
        spin = mo["spin"].lower()
        occup = mo["occup"]
        coeffs = mo["coeffs"]
        if spin == "alpha":
            occ_a += occup
            C_a.append(coeffs)
        elif spin == "beta":
            occ_b += occup
            C_b.append(coeffs)
        else:
            raise Exception(f"Unknown {spin=}!")

        if not unrestricted:
            unrestricted = spin == "beta"
    occ_a = int(occ_a)
    occ_b = int(occ_b)
    if not unrestricted:
        C_b = C_a.copy()
        occ_a = occ_b = occ_a // 2
        mult = 1
    else:
        mult = (occ_a - occ_b) + 1

    nuc_charge = nuc_charges_for_atoms(atoms).sum()
    charge = nuc_charge - occ_a - occ_b
    C = np.array((np.transpose(C_a), np.transpose(C_b)))

    _, nao, _ = C.shape
    bf_type = {
        shells.cart_size: BFType.CARTESIAN,
        shells.sph_size: BFType.PURE_SPHERICAL,
    }[nao]

    wf_kwargs = {
        "atoms": atoms,
        "coords": coords,
        "charge": charge,
        "mult": mult,
        "unrestricted": unrestricted,
        "occ": (occ_a, occ_b),
        "C": C,
        "bf_type": bf_type,
        "shells": shells,
        **kwargs,
    }

    # _wf_kwargs.update(wf_kwargs)
    return Wavefunction(**wf_kwargs)


@file_or_str(*AOMIX_EXTS)
def shells_from_aomix(text):
    aomix_dict = parse_aomix(text)
    return shells_from_aomix_dict(aomix_dict)


@file_or_str(*AOMIX_EXTS)
def wavefunction_from_aomix(text, **kwargs):
    aomix_dict = parse_aomix(text)
    return wavefunction_from_aomix_dict(aomix_dict, **kwargs)
