import numpy as np

from pysisyphus.elem_data import COVALENT_RADII as CR
from pysisyphus.intcoords.exceptions import DifferentPrimitivesException
from pysisyphus.intcoords.RedundantCoords import RedundantCoords
from pysisyphus.intcoords.setup import get_bond_sets, BOND_FACTOR
from pysisyphus.intcoords.Stretch import Stretch


def get_tangent(prims1, prims2, dihedral_inds, normalize=False):
    """Normalized tangent between primitive internal coordinates.

    Tangent pointing from prims1 to prims2 in primitive
    internal coordinates, taking into account the periodicity of
    dihedral angles.

    Parameters
    ----------
    prims1 : np.array
        1d-array of primitive internal coordinates in the order
        (stretches, bends, dihedrals).
    prims2 : np.array
        See prims1.
    dihedral_inds : list of int
        Dihedral indices in prims1 and prims2.

    Returns
    -------
    tangent : np.array
        1d array containing the normalized tangent pointing from prims1 to prims2.
    """
    try:
        tangent = prims2 - prims1
    except ValueError:
        raise DifferentPrimitivesException
    diheds = tangent[dihedral_inds].copy()
    diheds_plus = diheds.copy() + 2 * np.pi
    diheds_minus = diheds.copy() - 2 * np.pi
    bigger = np.abs(diheds) > np.abs(diheds_plus)
    diheds[bigger] = diheds_plus[bigger]
    bigger = np.abs(diheds) > np.abs(diheds_minus)
    diheds[bigger] = diheds_minus[bigger]
    tangent[dihedral_inds] = diheds
    if normalize:
        tangent /= np.linalg.norm(tangent)
    return tangent


def get_step(geom, coords):
    assert len(geom.coords) == len(coords)

    if geom.coord_type == "cart":
        diff = geom.coords - coords
    elif geom.coord_type in ("redund", "dlc"):
        diff = -get_tangent(
            geom.internal.prim_coords, coords, geom.internal.dihedral_inds
        )
    else:
        raise Exception("Invalid coord_type!")

    # Convert to DLC
    if geom.coord_type == "dlc":
        diff = geom.internal.U.T.dot(diff)

    return diff


def merge_coordinate_definitions(geom1, geom2):
    typed_prims1 = geom1.internal.typed_prims
    typed_prims2 = geom2.internal.typed_prims
    union = list(set(typed_prims1) | set(typed_prims2))
    union.sort()
    # Check if internal coordinates that are only present in one of the two
    # geometries are valid in the other. If not, we omit these primitives.
    redundant = RedundantCoords(geom1.atoms, geom1.cart_coords, typed_prims=union)
    valid_typed_prims = redundant.typed_prims
    return valid_typed_prims


def form_coordinate_union(geom1, geom2):
    def not_cartesian(geom):
        return geom.coord_type != "cart"

    assert not_cartesian(geom1) and not_cartesian(geom2)
    # The first call yields all primitives from geom1 that are also valid at geom2.
    typed_prims1 = merge_coordinate_definitions(geom1, geom2)
    # The second call yields all primitives from geom2 that are also valid at geom1.
    typed_prims2 = merge_coordinate_definitions(geom2, geom1)

    intersection = list(set(typed_prims1) & set(typed_prims2))
    intersection.sort()
    return intersection


def get_weighted_bond_mode(weighted_bonds, coords3d, remove_translation=True):
    bond_mode = np.zeros_like(coords3d.flatten())
    for *indices, weight in weighted_bonds:
        _, grad = Stretch._calculate(coords3d, indices, gradient=True)
        """
        The gradient gives us the direction into which the bond increases, but
        we want that positive weights correspond to bond formation (distance
        decrease) and negative weights to bond breaking (distance increase),
        so we reverse the sign of the weight.
        """
        bond_mode += -weight * grad

    if remove_translation:
        bm3d = bond_mode.reshape(-1, 3)
        bond_mode = (bm3d - bm3d.mean(axis=0)[None, :]).flatten()

    bond_mode /= np.linalg.norm(bond_mode)
    return bond_mode


def get_weighted_bond_mode_getter(
    target_weighted_bonds, bond_factor=1.2, fractional=False
):
    """Create input for intcoords.helpers.get_weighted_bond_mode.

    Compared to the rest of pysisyphus this method uses a slightly lowered
    bond factor, so it is more strict regarding what is considered a bond
    and what not."""

    def func(atoms, coords3d):
        bond_sets = [
            set(bond)
            for bond in get_bond_sets(atoms, coords3d, bond_factor=bond_factor).tolist()
        ]
        bonds = list()
        frac_weights = list()
        for tbond in target_weighted_bonds:
            *inds, weight = tbond
            assert weight != 0.0
            tset = set(inds)
            """
            Skip target bond if weight is positive (bond is to be formed) and
            bond is already present. Similarly, skip when the bond is to be broken
            and it is not present.
            """
            if (weight > 0 and tset in bond_sets) or (
                weight < 0 and tset not in bond_sets
            ):
                continue
            bonds.append(tbond)

            if fractional:
                from_, to_ = inds
                ref_len = CR[atoms[from_].lower()] + CR[atoms[to_].lower()]
                act_len = np.linalg.norm(coords3d[from_] - coords3d[to_])
                fw = (act_len / ref_len) ** weight
                frac_weights.append(fw)

        if fractional:
            frac_weights = np.array(frac_weights)
            frac_weights /= frac_weights.sum()
            frac_weights = frac_weights ** 2
            frac_weights /= frac_weights.max()
            bonds = [(*inds, fw) for ((*inds, _), fw) in zip(bonds, frac_weights)]
        return bonds

    return func


def get_bond_difference(geom1, geom2, bond_factor=BOND_FACTOR):
    """Return formed and broken bonds when going from geom1 to geom2."""

    assert geom1.atoms == geom2.atoms

    def as_sets(geom):
        bonds = get_bond_sets(geom.atoms, geom.coords3d, bond_factor=BOND_FACTOR)
        return set([frozenset(b) for b in bonds.tolist()])

    def as_lists(bonds):
        return [list(b) for b in bonds]

    bonds1 = as_sets(geom1)
    bonds2 = as_sets(geom2)
    formed = as_lists(bonds2 - bonds1)
    broken = as_lists(bonds1 - bonds2)
    return formed, broken


def verbose_bond_difference(formed, broken, key1, key2, atoms=None):
    def atom_bonds(bonds):
        if len(bonds) == 0:
            ab = ""
        elif atoms is not None:
            ab = [f"{atoms[from_]}{from_}-{atoms[to_]}{to_}" for from_, to_ in bonds]
        else:
            ab = bonds
        return ab

    def pos(bonds):
        """Plural or singular?"""
        return ("", "is") if (len(bonds) == 1) else ("s", "are")

    sf, verbf = pos(formed)
    sb, verbb = pos(broken)

    # return (
    # f"{key1}->{key2}: "
    # f"{len(formed)} bond{sf} formed {atom_bonds(formed)}, "
    # f"{len(broken)} bond{sb} broken {atom_bonds(broken)}."
    # )
    return (
        f"{len(formed)} bond{sf} formed {atom_bonds(formed)}",
        f"{len(broken)} bond{sb} broken {atom_bonds(broken)}",
    )


def get_bond_differences_verbose(
    geom1, geom2, bond_factor=BOND_FACTOR, key1="geom1", key2="geom2"
):
    formed, broken = get_bond_difference(geom1, geom2, BOND_FACTOR)
    fverb, bverb = verbose_bond_difference(formed, broken, key1, key2, atoms=geom1.atoms)
    return formed, broken, fverb, bverb
