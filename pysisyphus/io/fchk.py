import re

import numpy as np

from pysisyphus.elem_data import INV_ATOMIC_NUMBERS
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.wavefunction.shells import Shell, FCHKShells
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.helpers import BFType


@file_or_str(".fchk")
def parse_fchk(text):
    header_re = r"\s+".join(
        (
            r"(?P<keyword>[\w\-\s/\(\)=\*]+)",
            r"(?P<kind>[HIRCL])",
            r"(?P<length>N=)?",
            r"(?P<num>[\d\-\+\.E]+)",
        )
    )
    header_re = re.compile(header_re)

    cur_keyword = None
    lines = text.strip().split("\n")
    # title1, titl2 are not used
    _, _, *lines = lines

    conv_funcs = {
        "R": lambda line: [float(_) for _ in line.split()],  # real
        "I": lambda line: [int(_) for _ in line.split()],  # integer
        "C": lambda line: [line],  # character
        "H": lambda line: [bool(_) for _ in line.split()],  # logical
        "L": lambda line: [bool(_) for _ in line.split()],  # logical
    }

    data = {}
    for line in lines:
        line = line.strip()
        if mobj := header_re.match(line):
            cur_keyword = mobj.group("keyword").strip()
            conv_func = conv_funcs[mobj.group("kind")]
            if mobj.group("length"):
                data[cur_keyword] = list()
            else:
                data[cur_keyword] = conv_func(mobj.group("num"))[0]
            continue
        data[cur_keyword].extend(conv_func(line))
    return data


@file_or_str(".fchk")
def shells_from_fchk(text):
    data = parse_fchk(text)

    def to_arr(key, dtype):
        return np.array(data[key], dtype=dtype)

    """
    These coordinates may not coincide with the basis functions?! Nonetheless,
    basis functions are usually atom-centered. But just to be sure we use the
    entries from 'Coordinates of each shell'.
    # coords3d = to_arr("Current cartesian coordinates", float).reshape(-1, 3)
    """
    atomic_nums = to_arr("Atomic numbers", int)

    centers = to_arr("Shell to atom map", int) - 1
    ang_moms = to_arr("Shell types", int)
    prims_per_shell = to_arr("Number of primitives per shell", int)

    prim_exponents = to_arr("Primitive exponents", float)
    prim_coeffs = to_arr("Contraction coefficients", float)
    try:
        p_prim_coeffs = to_arr("P(S=P) Contraction coefficients", float)
    except KeyError:
        p_prim_coeffs = np.zeros_like(prim_coeffs)
    shell_coords3d = to_arr("Coordinates of each shell", float).reshape(-1, 3)

    _shells = list()
    prim_ind = 0
    for i, (center_ind, L, prim_num) in enumerate(
        zip(centers, ang_moms, prims_per_shell)
    ):
        # 0: s, 1: p, -1: sp, 2: 6d, -2: 5d, 3: 10f, -3: 7f, etc...
        assert L <= 1, "Only spherical basis functions are supported."
        # Convert SP-shell with L = -1 to a s-shell with L = 0.
        # The p-shell is also created later.
        L, is_sp = (0, True) if (L == -1) else (L, False)

        coeffs = list()
        p_coeffs = list()
        exps = list()
        for _ in range(prim_num):
            coeffs.append(prim_coeffs[prim_ind])
            p_coeffs.append(p_prim_coeffs[prim_ind])
            exps.append(prim_exponents[prim_ind])
            prim_ind += 1
        atomic_num = atomic_nums[center_ind]
        center = shell_coords3d[i]

        def append_shell(coeffs, L=L):
            shell = Shell(
                L=abs(L),
                atomic_num=atomic_num,
                center=center,
                coeffs=coeffs,
                exps=exps,
                center_ind=center_ind,
            )
            _shells.append(shell)

        # Create shell with actual L.
        # If we deal with a SP shell we also create the P shell below.
        append_shell(coeffs)
        if is_sp:
            # Explicitly create p-shell
            append_shell(p_coeffs, L=1)
    shells = FCHKShells(_shells)
    return shells


def atoms_from_data(data):
    atomic_numbers = data["Atomic numbers"]
    atoms = tuple([INV_ATOMIC_NUMBERS[Z] for Z in atomic_numbers])
    return atoms


@file_or_str(".fchk")
def wavefunction_from_fchk(text):
    data = parse_fchk(text)

    charge = data["Charge"]
    mult = data["Multiplicity"]
    atoms = atoms_from_data(data)
    coords = data["Current cartesian coordinates"]
    unrestricted = "Beta Orbital Energies" in data

    def mos_for_spin(spin):
        mo_ens = data[f"{spin} Orbital Energies"]
        mo_num = len(mo_ens)
        mo_coeffs = np.array(data[f"{spin} MO coefficients"])
        mo_coeffs = mo_coeffs.reshape(mo_num, mo_num)
        return mo_coeffs

    C_a = mos_for_spin("Alpha").T
    occ_a = data["Number of alpha electrons"]
    if unrestricted:
        C_b = mos_for_spin("Beta").T
        occ_b = data["Number of beta electrons"]
    else:
        C_b = C_a.copy()
        occ_b = occ_a
    C = np.stack((C_a, C_b))
    shells = shells_from_fchk(text)
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
    )


@file_or_str(".fchk")
def geom_from_fchk(text, **geom_kwargs):
    data = parse_fchk(text)
    atoms = atoms_from_data(data)

    try:
        coords = data["Opt point       1 Geometries"]
    except KeyError:
        coords = data["Current cartesian coordinates"]
    coords = np.array(coords)
    coords = coords.reshape(-1, 3 * len(atoms))

    geoms = [Geometry(atoms, coords_, **geom_kwargs) for coords_ in coords]
    if len(geoms) == 1:
        geoms = geoms[0]
    return geoms
