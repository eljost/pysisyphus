import numpy as np

def make_xyz_str(atoms, coords, comment=""):
    assert(len(atoms) == len(coords))

    coord_fmt = "{: 03.8f}"
    line_fmt = "{:>3s} " + " ".join([coord_fmt, ]*3)

    body = [line_fmt.format(a, *xyz)
            for a, xyz
            in zip(atoms, coords)]
    body = "\n".join(body)
 
    return "{}\n{}\n{}".format(len(atoms), comment, body)


def make_trj_str(atoms, coords_list):
    xyz_strings = [make_xyz_str(atoms, coords) for coords in coords_list]
    return "\n".join(xyz_strings)


def parse_xyz_str(xyz_str):
    """Parse a xyz string.

    Paramters
    ---------
    xyz_str : str
        The contents of a .xyz file.

    Returns
    -------
    atoms : list
        List of length N (N = number of atoms) holding the
        element symbols.
    coords: np.array
        An array of shape (N, 3) holding the xyz coordinates.
    """

    # Only consider the first four items on a line
    atoms_coords = [line.strip().split()[:4]
                    for line in xyz_str.strip().split("\n")[2:]
    ]
    atoms, coords = zip(*[(a, c) for a, *c in atoms_coords])
    coords = np.array(coords, dtype=np.float)
    return atoms, coords


def parse_xyz_file(xyz_fn):
    with open(xyz_fn) as handle:
        xyz_str = handle.read()

    return parse_xyz_str(xyz_str)


def parse_trj_file(trj_fn):
    with open(trj_fn) as handle:
        trj_str = handle.read()

    return parse_trj_str(trj_str)


def parse_trj_str(trj_str):
    trj_lines = trj_str.strip().split("\n")
    number_of_atoms = int(trj_lines[0].strip())
    xyz_lines = number_of_atoms + 2
    # Split the trj file in evenly sized strings
    xyz_strs = ["\n".join(trj_lines[i:i+xyz_lines])
                for i in range(0, len(trj_lines), xyz_lines)]
    xyzs = [parse_xyz_str(xyz_str) for xyz_str in xyz_strs]

    assert(len(xyzs) == (len(trj_lines) / xyz_lines))
    return xyzs
