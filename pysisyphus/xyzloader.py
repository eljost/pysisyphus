import re

import numpy as np

from pysisyphus.constants import BOHR2ANG


def make_xyz_str(atoms, coords, comment=""):
    assert len(atoms) == len(coords)

    coord_fmt = "{: 03.8f}"
    line_fmt = "{:>3s} " + " ".join(
        [
            coord_fmt,
        ]
        * 3
    )

    body = [line_fmt.format(a, *xyz) for a, xyz in zip(atoms, coords)]
    body = "\n".join(body)

    return f"{len(atoms)}\n{comment}\n{body}"


def make_trj_str(atoms, coords_list, comments=None):
    if comments is None:
        comments = ["" for _ in coords_list]
    xyz_strings = [
        make_xyz_str(atoms, coords, comment)
        for coords, comment in zip(coords_list, comments)
    ]
    return "\n".join(xyz_strings)


def coords_to_trj(trj_fn, atoms, coords_list, comments=None):
    coords = np.array(coords_list)
    coords = coords.reshape(-1, len(atoms), 3) * BOHR2ANG
    trj_str = make_trj_str(atoms, coords, comments)
    with open(trj_fn, "w") as handle:
        handle.write(trj_str)
    return trj_fn


def make_trj_str_from_geoms(geoms, comments=None, energy_comments=False):
    atoms = geoms[0].atoms
    coords_list = [geom.coords3d * BOHR2ANG for geom in geoms]
    if energy_comments and comments is None:
        comments = [str(geom._energy) for geom in geoms]
    elif comments is not None:
        assert len(comments) == len(geoms)

    return make_trj_str(atoms, coords_list, comments)


def write_geoms_to_trj(geoms, fn, comments=None):
    trj_str = make_trj_str_from_geoms(geoms, comments)
    with open(fn, "w") as handle:
        handle.write(trj_str)


def split_xyz_str(xyz_str):
    """Example:

    xyz:
     1

     X -1.2 1.4 0.0
     1

     X 2.0 4.0 0.0

    """
    float_ = r"([\+\d\-\.]+)"
    header_re = re.compile(r"(\d+)")
    coord_re = re.compile(fr"[a-zA-Z]+\s+{float_}\s+{float_}\s+{float_}")

    lines = [l.strip() for l in xyz_str.strip().split("\n")]

    lines_remaining = len(lines)
    cur_line = 0
    valid_xyz_strs = list()
    while lines_remaining:
        header_mobj = header_re.match(lines[cur_line])
        expect_lines = int(header_mobj[1])
        slice_ = slice(
            cur_line + 2, cur_line + 2 + expect_lines
        )  # lgtm [py/hash-unhashable-value]
        check_lines = lines[slice_]
        assert len(check_lines) == expect_lines
        assert all([coord_re.match(line.strip()) for line in check_lines])

        valid_xyz_strs.append(str(expect_lines) + "\n\n" + "\n".join(check_lines))

        lines_read = expect_lines + 2
        lines_remaining -= lines_read
        cur_line += lines_read
    atoms_coords = [parse_xyz_str(_, with_comment=False) for _ in valid_xyz_strs]
    return atoms_coords


def parse_xyz_str(xyz_str, with_comment):
    """Parse a xyz string.

    Paramters
    ---------
    xyz_str : str
        The contents of a .xyz file.
    with_comment : bool
        Return comment line if True.

    Returns
    -------
    atoms : list
        List of length N (N = number of atoms) holding the
        element symbols.
    coords: np.array
        An array of shape (N, 3) holding the xyz coordinates.
    comment_line : str, optional
        Comment line if with_comment argument was True.
    """

    xyz_lines = xyz_str.strip().split("\n")
    comment_line = xyz_lines[1]

    # Only consider the first four items on a line
    atom_num = int(xyz_lines[0])
    atoms_present = len(xyz_lines) - 2
    assert len(xyz_lines) == atom_num + 2, \
        f"Expected {atom_num} atoms, but found only {atoms_present}!"
    atoms_coords = [
        line.strip().split()[:4] for line in xyz_str.strip().split("\n")[2:]
    ]
    atoms, coords = zip(*[(a, c) for a, *c in atoms_coords])
    coords = np.array(coords, dtype=float)
    if with_comment:
        return atoms, coords, comment_line
    else:
        return atoms, coords


def parse_xyz_file(xyz_fn, with_comment=False):
    with open(xyz_fn) as handle:
        xyz_str = handle.read()

    return parse_xyz_str(xyz_str, with_comment)


def parse_trj_file(trj_fn, with_comments=False):
    with open(trj_fn) as handle:
        trj_str = handle.read()

    return parse_trj_str(trj_str, with_comments)


def parse_trj_str(trj_str, with_comments=False):
    trj_lines = trj_str.strip().split("\n")
    number_of_atoms = int(trj_lines[0].strip())
    xyz_lines = number_of_atoms + 2
    # Split the trj file in evenly sized strings
    xyz_strs = [
        "\n".join(trj_lines[i : i + xyz_lines])
        for i in range(0, len(trj_lines), xyz_lines)
    ]
    xyzs = [parse_xyz_str(xyz_str, with_comments) for xyz_str in xyz_strs]

    assert len(xyzs) == (len(trj_lines) / xyz_lines)
    return xyzs
