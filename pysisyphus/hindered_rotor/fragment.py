import functools

import networkx as nx
import numpy as np
import scipy as sp

from pysisyphus.Geometry import Geometry
from pysisyphus.intcoords.PrimTypes import PrimTypes as PT


def fragment_geom(
    geom: Geometry, rm_bond: tuple[int, int]
) -> tuple[list[int], list[int]]:
    """Split geometry comprising N fragments into N+1 fragments by bond deletion.

    This function is probably not completely general, but should support common cases
    like an initially bonded molecule that is fragmented by removal of a bond and
    non-bonded van-der-Waals complexes.

    The method will fail for circular systems where removal of one bond breaks the ring.
    """
    geom_bond = geom.copy(
        coord_type="redund",
        coord_kwargs={
            "bonds_only": True,
            # Here, we explicitly define the bond to be removed to support geometrie,
            # where the bond is not actually present, as in non-covalently bound strutures.
            "define_prims": [(PT.BOND, *rm_bond)],
        },
    )
    # Only continue with actual bonds; drop all interfragment bonds etc.
    act_bonds = [
        indices for pt, *indices in geom_bond.internal.typed_prims if pt == PT.BOND
    ]
    # Create a graph from the bond topology
    G = nx.from_edgelist(act_bonds)
    ncomponents_org = len(list(nx.connected_components(G)))

    # Remove the selected bond to fragment the graph
    G.remove_edge(*rm_bond)
    components_frag = list(nx.connected_components(G))
    ncomponents_frag = len(components_frag)
    # Check, that removal of the bond actually fragmented the system and increased the number
    # of components by one.
    assert ncomponents_frag == ncomponents_org + 1

    sg_left, sg_right = [G.subgraph(c).copy() for c in components_frag]
    nodes_left = list(map(int, sg_left))
    nodes_right = list(map(int, sg_right))
    return nodes_left, nodes_right


@functools.singledispatch
def rotate_fragment_around_bond(
    coords3d: np.ndarray, bond: tuple[int, int], fragment: list[int], deg: float
) -> np.ndarray:
    """Rotate fragment bonded to the 2nd atom in bond around said bond."""
    # Determine unit vector pointing along bond
    from_, to_ = bond
    bond_vec = coords3d[from_] - coords3d[to_]
    bond_vec /= np.linalg.norm(bond_vec)

    # Extract fragment coordiantes
    frag_coords3d = coords3d[fragment].copy()

    # Degrees of rotation are encoded in the length of the vector
    Rot = sp.spatial.transform.Rotation.from_rotvec(deg * bond_vec, degrees=True)
    frag_coords3d_rot = Rot.apply(frag_coords3d)

    # Construct updated coordinates, that contain the rotated fragment
    coords3d_rot = coords3d.copy()
    coords3d_rot[fragment] = frag_coords3d_rot
    return coords3d_rot


@rotate_fragment_around_bond.register
def _(geom: Geometry, bond: tuple[int, int], deg: float) -> np.ndarray:
    # Note the [1] at the end ... we drop the left frag and only continue with
    # the right fragment.
    right_frag = fragment_geom(geom, bond)[1]
    return rotate_fragment_around_bond(geom.coords3d, bond, right_frag, deg)
