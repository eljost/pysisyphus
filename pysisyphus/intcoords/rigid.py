import functools
import itertools as it
from typing import List, Optional, Set, Tuple, Union
import warnings

import numpy as np
from numpy.typing import NDArray
import networkx as nx
from scipy.spatial.transform import Rotation as spRot

from pysisyphus.Geometry import Geometry
from pysisyphus.intcoords.PrimTypes import PrimTypes, normalize_prim_input
from pysisyphus.intcoords.setup_fast import find_bonds
from pysisyphus.linalg import cross3, norm3


ListOrTuple = Union[List, Tuple]


@functools.singledispatch
def to_graph(atoms: ListOrTuple, coords3d: NDArray[float]) -> nx.Graph:
    """Construct nx.Graph for given atoms and 3c Cartesian coordinates."""
    bonds = find_bonds(atoms, coords3d)
    bonded_set = set(it.chain(*[bond.tolist() for bond in bonds]))
    full_set = set(range(len(atoms)))
    G = nx.from_edgelist(bonds)
    for node in full_set - bonded_set:
        G.add_node(node)
    return G


@to_graph.register
def _(geom: Geometry):
    return to_graph(geom.atoms, geom.coords3d)


def get_fragments(G: nx.Graph, del_bonds: Optional[List] = None) -> List[Set]:
    """Determine fragments in graph w/ possible bond deletion."""
    G_cut = G.copy()
    if del_bonds is None:
        del_bonds = list()

    for del_bond in del_bonds:
        u, v = del_bond
        try:
            G_cut.remove_edge(u, v)
        except nx.exception.NetworkXError:
            warnings.warn(f"'{del_bond}' bond is not present!")
    comps = nx.connected_components(G_cut)
    return list(comps)


def get_connected_groups(frags, anchors):
    """Determine groups that are connected to given anchors.

    Given n anchors, at most n groups are returned; possibly less."""
    frag_sets = [set(frag) for frag in frags]
    groups = list()
    for anchor in anchors:
        for frag in frag_sets:
            if anchor in frag:
                groups.append(list(frag))
                break
    return groups


def vec_rot(c3d, vec, group1, group2, rad):
    vec = vec / np.linalg.norm(vec)
    deg = np.rad2deg(rad)
    deg2 = deg / 2.0  # Rotate both groups in opposite directions by deg/2

    if (len(group1) == 1) or (len(group2) == 1):
        warnings.warn("Found 1-atom fragment.")
        # Double degrees; otherwise the other group will only do a half-rotation.
        deg *= 2

    def rot_group(group, pos=True):
        factor = 1.0 if pos else -1.0
        return spRot.from_rotvec(factor * deg2 * vec, degrees=True).apply(c3d[group])

    c3d[group1] = rot_group(group1)
    c3d[group2] = rot_group(group2, False)
    return c3d


def vec_trans(c3d, vec, group1, group2, dist):
    vec = vec / np.linalg.norm(vec)

    def trans_group(group, pos=True):
        factor = 1.0 if pos else -1.0
        return c3d[group] + factor * vec[None, :] * dist

    c3d[group1] = trans_group(group1)
    c3d[group2] = trans_group(group2, False)
    return c3d


def get_vec_trans_step_func(
    atoms: ListOrTuple,
    coords3d: NDArray[float],
    indices: List,
):
    assert len(indices) == 2
    bond = indices
    G = to_graph(atoms, coords3d)
    frags = get_fragments(G, (bond,))
    if len(frags) == 1:
        warnings.warn("Found only one fragment. Moving only bonded atoms.")
        raise Exception("Not yet implemented.")
    gr1, gr2 = get_connected_groups(frags, bond)

    from_, to_ = bond
    bond_vec = coords3d[from_] - coords3d[to_]
    bond_vec /= np.linalg.norm(bond_vec)

    coords3d = coords3d.copy()

    def stepper(step):
        return vec_trans(coords3d, bond_vec, gr1, gr2, step).copy()

    return stepper


def get_rot_vec(coords3d, indices):
    # Torsion; use central bond vector.
    if len(indices) == 2:
        from_, to_ = indices
        vec = coords3d[from_] - coords3d[to_]
    elif len(indices) == 4:
        _, from_, to_, _ = indices
        vec = coords3d[from_] - coords3d[to_]
    # Bend; use vector normal to bend plane.
    elif len(indices) == 3:
        m, o, n = indices
        u = coords3d[m] - coords3d[o]
        v = coords3d[n] - coords3d[o]
        vec = cross3(u, v)
        vec /= norm3(vec)
    else:
        raise Exception("Not implemented!")
    vec /= np.linalg.norm(vec)
    return vec


def get_vec_rot_step_func(
    atoms: ListOrTuple,
    coords3d: NDArray[float],
    indices: List,
):
    # Indices for rotation vector from bond. All groups connected to the bonded
    # atoms will be rotated.
    if len(indices) == 2:
        del_bonds = [
            indices,
        ]
        anchors = indices
    # Indices for rotation from vector normal to the plane of a bend.
    # Two bonds will be deleted to determine the groups.
    elif len(indices) == 3:
        del_bonds = [indices[:2], indices[1:]]
        anchors = indices[0], indices[2]
    # Indices for rotation vector from central bond in a torsion. Only the atoms
    # directly connected to the terminal atoms of the torsion will be rotated.
    elif len(indices) == 4:
        m, o, p, n = indices
        del_bonds = [
            [m, o],
            [n, p],
            [m, n],
        ]
        anchors = [m, n]
    else:
        raise Exception("'indices' must be an iterable of 3 or 4 integers!")
    anchors = list(anchors)

    # Determine fragments ...
    G = to_graph(atoms, coords3d)
    frags = get_fragments(G, del_bonds)
    # ... and do some basic sanity checking.
    nfrags = len(frags)
    single_atom_frags = [i for i, frag in enumerate(frags) if len(frag) == 1]
    if nfrags == 1:
        warnings.warn("Only one fragment. System is invariant to rotation.")
    elif len(single_atom_frags) == nfrags:
        warnings.warn("Only 1-atom fragments. They are invariant to rotation.")

    # Determine groups that will be rotated
    gr1, gr2 = get_connected_groups(frags, anchors)
    # Determine the rotation vector.
    rot_vec = get_rot_vec(coords3d, indices)
    # Don't operate on the original coordinates, but on a copy of it.
    coords3d = coords3d.copy()

    def stepper(rad):
        return vec_rot(coords3d, rot_vec, gr1, gr2, rad).copy()

    return stepper


def get_stepper(atoms, coords3d, typed_prim, step_size, nsteps):
    # TODO: Also return value of internal coordinate from stepper_func
    #       This is probably not always possible, e.g., when rotating around a bond
    #       or around a vector orthogonal to a plane.
    # TODO: Don't rotate atoms in cycles, similar to GaussView.
    # TODO: Handle corner cases.

    typed_prim = normalize_prim_input(typed_prim)[0]
    pt, *indices = typed_prim
    step_func_getters = {
        PrimTypes.BOND: get_vec_trans_step_func,  # Translation of bonded groups
        PrimTypes._ROT_BOND: get_vec_rot_step_func,  # Rotation of all bonded groups
        PrimTypes.BEND: get_vec_rot_step_func,  # Rotation normal to bend plane
        # Rotation of groups bonded to terminal atoms of torsion
        PrimTypes.PROPER_DIHEDRAL: get_vec_rot_step_func,
    }
    step_func = step_func_getters[pt](atoms, coords3d, indices)
    cum_step_size = 0.0  # Cumulated step size
    step_sizes = list()

    def stepper():
        def take_step(step_size):
            nonlocal cum_step_size

            cum_step_size += step_size
            step_sizes.append(step_size)
            return cum_step_size, step_func(step_size)

        # Yield the original coordinates in the first cycle.
        yield take_step(0.0)

        for _ in range(nsteps):
            yield take_step(step_size)

    return stepper


@functools.singledispatch
def rot_around_bond(
    atoms: ListOrTuple, coords3d: NDArray[float], bond: List, rad: float
):
    return get_vec_rot_step_func(atoms, coords3d, indices=bond)(rad)


@rot_around_bond.register
def _(geom: Geometry, bond, rad):
    return rot_around_bond(geom.atoms, geom.coords3d, bond, rad)
