from collections import namedtuple

import numpy as np

from pysisyphus.constants import ANG2BOHR as ANG2BOHR
from pysisyphus.Geometry import Geometry


ZLine = namedtuple(
    "ZLine", "atom rind r aind a dind d", defaults=(None, None, None, None, None, None)
)


def geom_from_zmat(
    zmat,
    atoms=None,
    coords3d=None,
    geom=None,
    start_at=None,
    drop_dummy=True,
    **geom_kwargs
):
    """Adapted from https://github.com/robashaw/geomConvert by Robert Shaw."""

    if isinstance(zmat, str):
        zmat = zmat_from_str(zmat)

    zmat_atoms = [zline.atom for zline in zmat]
    # Extend supplied geometry by zmat
    if geom is not None:
        atoms = geom.atoms
        coords3d = geom.coords3d
        start_at = len(geom.atoms)

    if atoms is not None:
        atoms = list(atoms) + zmat_atoms
    else:
        atoms = zmat_atoms

    # Grow coordindate array and assign old coordinates
    if coords3d is not None:
        _coords3d = np.zeros((len(coords3d) + len(zmat), 3))
        _coords3d[: len(coords3d)] = coords3d
        coords3d = _coords3d
    else:
        coords3d = np.zeros((len(zmat), 3), dtype=float)

    if atoms or coords3d:
        assert len(coords3d) == len(atoms)

    if start_at is None:
        start_at = 0

    for i, zline in enumerate(zmat, start_at):
        assert all(
            [
                (ind is None) or (ind >= 0)
                for ind in (zline.rind, zline.aind, zline.dind)
            ]
        ), "Found invalid atom index. Atom indices start with 1, not 0!"

        r = zline.r
        # First atom is placed at the origin
        if i == 0:
            continue
        # Bond along x-axis
        elif i == 1:
            coords3d[i, 0] = r
        # Angle in xy-plane from polar coordinates
        elif i == 2:
            r"""
            M       P <- add
             \     /
              u   v
               \ /
                O
            """
            theta = np.deg2rad(zline.a)
            # Center
            O = coords3d[zline.rind]
            # Bond, pointing away from O to M
            u = coords3d[zline.aind] - O
            # Direction of u along x axis (left/right)
            sign = np.sign(u[0])
            # Polar coordinates
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            # Translate from center with correct orientation
            coords3d[i] = O + (sign * x, sign * y, 0.0)
        # Dihedral in xyz-space from spherical coordinates
        else:
            theta, phi = np.deg2rad((zline.a, zline.d))

            sin_theta = np.sin(theta)
            cos_theta = np.cos(theta)
            sin_phi = np.sin(phi)
            cos_phi = np.cos(phi)

            x = r * cos_theta
            y = r * sin_theta * cos_phi
            z = r * sin_theta * sin_phi

            r"""
            M <- add      N
             \           /
              u         v
               \       /
                O--w--P
            """

            O = coords3d[zline.rind]
            P = coords3d[zline.aind]
            N = coords3d[zline.dind]

            # Local axis system
            v_ = P - N
            v = v_ / np.linalg.norm(v_)
            w_ = O - P
            w = w_ / np.linalg.norm(w_)
            a = np.cross(v, w)
            a /= np.linalg.norm(a)
            b = np.cross(a, w)
            b /= np.linalg.norm(b)
            coords3d[i] = O - w * x + b * y + a * z

    if drop_dummy:
        atoms_ = list()
        coords3d_ = list()
        for atom, xyz in zip(atoms, coords3d):
            if atom.lower() == "x":
                continue
            atoms_.append(atom)
            coords3d_.append(xyz)
        atoms = atoms_
        coords3d = np.array(coords3d_)

    geom = Geometry(atoms, coords3d, **geom_kwargs)
    return geom


def geom_from_zmat_str(text, coord_type="cart", coord_kwargs=None):
    zmat = zmat_from_str(text)
    return geom_from_zmat(zmat, coord_type=coord_type, coord_kwargs=coord_kwargs)


def zmat_from_str(text):
    def dec(str_):
        return int(str_) - 1

    def to_bohr(str_):
        return float(str_) * ANG2BOHR

    def convert(items):
        funcs = (str, dec, to_bohr, dec, float, dec, float)
        return [f(item) for f, item in zip(funcs, items)]

    zmat = list()
    for line in text.strip().split("\n"):
        line = line.strip()
        if line.startswith("#"):
            continue
        zmat.append(ZLine(*convert(line.split())))
    return zmat


def zmat_from_fn(fn):
    with open(fn) as handle:
        text = handle.read()
    return zmat_from_str(text)


def geom_from_zmat_fn(fn, **geom_kwargs):
    zmat = zmat_from_fn(fn)
    geom = geom_from_zmat(zmat, **geom_kwargs)
    return geom
