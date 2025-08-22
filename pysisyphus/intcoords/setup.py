# [1] http://link.aip.org/link/doi/10.1063/1.1515483
#     The efficient optimization of molecular geometries using redundant internal
#     coordinates
#     Bakken, Helgaker, 2002
# [2] https://doi.org/10.1063/1.479510
#     Geometry optimization in generalized natural internal coordinates
#     von Arnim, Ahlrichs, 1999
# [3] https://doi.org/10.1002/jcc.21494
#     The choice of internal coordinates in complex chemical systems
#     Nemeth, Challacombe, Van Veenendaal, 2010

from collections import namedtuple
import itertools as it
from typing import FrozenSet, Optional, Set
import math

import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans

from pysisyphus.constants import BOHR2ANG
from pysisyphus.config import BEND_MIN_DEG, DIHED_MAX_DEG
from pysisyphus.helpers_pure import log, sort_by_central, merge_sets
from pysisyphus.elem_data import VDW_RADII, COVALENT_RADII as CR
from pysisyphus.intcoords import Stretch, Bend, LinearBend, Torsion
from pysisyphus.intcoords.setup_fast import find_bonds as find_bonds_fast
from pysisyphus.intcoords.PrimTypes import (
    Bonds as PTBonds,
    PrimTypes,
    PrimMap,
    Rotations,
)
from pysisyphus.intcoords.valid import bend_valid, dihedral_valid


BOND_FACTOR = 1.3


def get_pair_covalent_radii(atoms):
    atoms = [a.lower() for a in atoms]
    cov_radii = np.array([CR[a.lower()] for a in atoms])
    pair_cov_radii = np.array([r1 + r2 for r1, r2 in it.combinations(cov_radii, 2)])
    return pair_cov_radii


def get_bond_mat(geom, bond_factor=BOND_FACTOR):
    cdm = pdist(geom.coords3d)
    pair_cov_radii = get_pair_covalent_radii(geom.atoms)
    bond_mat = squareform(cdm <= (pair_cov_radii * bond_factor))
    return bond_mat


def get_bond_sets(
    atoms, coords3d, bond_factor=BOND_FACTOR, return_cdm=False, return_cbm=False
):
    """I'm sorry, but this function does not return sets, but an int ndarray."""
    cdm = pdist(coords3d)
    # Generate indices corresponding to the atom pairs in the
    # condensed distance matrix cdm.
    atom_inds = list(it.combinations(range(len(coords3d)), 2))
    atom_inds = np.array(atom_inds, dtype=int)
    scaled_cr_sums = bond_factor * get_pair_covalent_radii(atoms)
    # condensed bond matrix
    cbm = cdm <= scaled_cr_sums
    bond_inds = atom_inds[cbm]
    if not return_cbm and not return_cdm:
        return bond_inds
    add_returns = tuple(
        [mat for flag, mat in ((return_cdm, cdm), (return_cbm, cbm)) if flag]
    )
    return (bond_inds,) + add_returns


def get_fragments(
    atoms,
    coords,
    bond_inds=None,
    bond_factor=BOND_FACTOR,
    ignore_atom_inds=None,
    ignore_bonds=None,
    with_unconnected_atoms=False,
):
    """This misses unconnected single atoms w/ 'with_unconnected_atoms=False'!"""
    coords3d = coords.reshape(-1, 3)
    if ignore_atom_inds is None:
        ignore_atom_inds = list()
    if ignore_bonds is None:
        ignore_bonds = set()
    ignore_bonds = set([frozenset(ib) for ib in ignore_bonds])
    ignore_atom_inds = set(ignore_atom_inds)

    if bond_inds is None:
        # Bond indices without interfragment bonds and/or hydrogen bonds
        bond_inds = get_bond_sets(atoms, coords3d, bond_factor=bond_factor)

    bond_ind_sets = [frozenset(bi) for bi in bond_inds]
    bond_ind_sets = [
        bi
        for bi in bond_ind_sets
        if (not bi & ignore_atom_inds) and (bi not in ignore_bonds)
    ]
    fragments = merge_sets(bond_ind_sets)

    # Also add single-atom fragments for unconnected atoms that don't participate
    # in any bonds.
    if with_unconnected_atoms:
        all_atoms = set([i for i, _ in enumerate(atoms)]) - ignore_atom_inds
        atoms_in_fragments = set(it.chain(*fragments))
        unconnected_atoms = all_atoms - atoms_in_fragments
        unconnected_fragments = [frozenset((atom,)) for atom in unconnected_atoms]
        fragments = fragments + unconnected_fragments
    return fragments


def connect_fragments(
    cdm, fragments, max_aux=3.78, aux_factor=BOND_FACTOR, logger=None
):
    """Determine the smallest interfragment bond for a list
    of fragments and a condensed distance matrix. For more than a few fragments
    this function performs poorly, as each fragment is connected to all reamaining
    fragments, leading to an explosion of bonds, bends and dihedrals."""
    if len(fragments) > 1:
        log(
            logger,
            f"Detected {len(fragments)} fragments. Generating interfragment bonds.",
        )
    dist_mat = squareform(cdm)
    interfrag_inds = list()
    aux_interfrag_inds = list()
    for frag1, frag2 in it.combinations(fragments, 2):
        log(logger, f"\tConnecting {len(frag1)} atom and {len(frag2)} atom fragment")
        inds = [(i1, i2) for i1, i2 in it.product(frag1, frag2)]
        distances = np.array([dist_mat[ind] for ind in inds])

        # Determine minimum distance bond
        min_ind = distances.argmin()
        min_dist = distances[min_ind]
        interfrag_bond = tuple(inds[min_ind])
        interfrag_inds.append(interfrag_bond)
        log(logger, f"\tMinimum distance bond: {interfrag_bond}, {min_dist:.4f} au")

        # Determine auxiliary interfragment bonds that are either below max_aux
        # (default 2 Å, ≈ 3.78 au), or less than aux_factor (default 1.3) times the
        # minimum interfragment distance.
        below_max_aux = [
            ind for ind in inds if (dist_mat[ind] < max_aux) and (ind != interfrag_bond)
        ]
        if below_max_aux:
            log(
                logger,
                f"\tAux. interfrag bonds below {max_aux*BOHR2ANG:.2f} Å:\n"
                + "\n".join(
                    [f"\t\t{ind}: {dist_mat[ind]:.4f} au" for ind in below_max_aux]
                ),
            )
        scaled_min_dist = aux_factor * min_dist
        above_min_dist = [
            ind
            for ind in inds
            if (dist_mat[ind] < scaled_min_dist)
            and (ind != interfrag_bond)
            and (ind not in below_max_aux)
        ]
        if above_min_dist:
            log(
                logger,
                f"\tAux. interfrag bonds below {aux_factor:.2f} * min_dist:\n"
                + "\n".join(
                    [f"\t\t{ind}: {dist_mat[ind]:.4f} au" for ind in above_min_dist]
                ),
            )
        aux_interfrag_inds.extend(below_max_aux)
        aux_interfrag_inds.extend(above_min_dist)
    # Or as Philipp proposed: two loops over the fragments and only
    # generate interfragment distances. So we get a full matrix with
    # the original indices but only the required distances.
    return interfrag_inds, aux_interfrag_inds


def connect_fragments_kmeans(
    cdm,
    fragments,
    atoms,
    aux_below_thresh=3.7807,  # 2 Å
    aux_add_dist=2.8356,  # 1.5 Å
    aux_keep=5,
    aux_no_hh=True,
    min_dist_thresh=5.0,
    min_scale=1.2,
    logger=None,
):
    """Generate (auxiliary) interfragment bonds.

    In the a first step, minimum distance interfragment bonds (IFBs) are determined between
    all possible fragment pairs. Similarly, possible auxiliary IFBs are determined.
    Candidates for auxiliary IFBs are:
        IFB <= aux_below_thresh, default 2 Å
        IFB <= (minimum distance IFB + aux_add_dist), default 1.5 Å
    By default, only the first aux_keep (default = 5) auxiliary IFBs are kept.

    Connecting all fragments can lead to bonds between very distant atoms. If more than
    two fragments are present we cluster the minimum distance IFB distances using KMeans,
    to determine a reasonable length for valid IFBs. We start out with two clusters and
    increase the number of cluster until the center of one cluster is around the scaled
    global minimum distance between the fragments. The center of this cluster is then used
    as a cutoff vor valid IFBs.

    After pruning all possible IFBs we can determine the fragment pairs, that are actually
    connected. This information is then used to also prune possible interfragment bonds.
    Only auxiliary IFBs between fragments that are actually connected via IFBs are kept.
    """
    # atoms = [atom.lower() for atom in atoms]
    if len(fragments) > 1:
        log(
            logger,
            f"Detected {len(fragments)} fragments. Generating interfragment bonds.",
        )
    dist_mat = squareform(cdm)

    frag_pairs = list()
    interfrag_inds = list()
    interfrag_dists = list()
    aux_dict = dict()
    for i, j in it.combinations(range(len(fragments)), 2):
        frag1 = fragments[i]
        frag2 = fragments[j]
        log(logger, f"\tConnecting {len(frag1)}-atom and {len(frag2)}-atom fragments.")
        # Pairs of possible interfragment bonds
        inds = np.array([(i1, i2) for i1, i2 in it.product(frag1, frag2)], dtype=int)
        distances = np.array([dist_mat[k, l] for k, l in inds])

        frag_pairs.append((i, j))
        # Determine minimum distance bond
        min_ind = distances.argmin()
        min_dist = distances[min_ind]
        interfrag_inds.append(inds[min_ind])
        interfrag_dists.append(min_dist)
        """
        But also consider other bonds with reasonable lengths as auxiliary interfrag bonds.
        Contrary to [1] we don't use a scaled minimum interfragment distances (MID),
        but just add a fixed value to the MID. Scaling a possibly big MID would
        probably include too many coordinates.
        """
        # if aux_no_hh:
        # is_hh = np.array(
        # [(atoms[k] == "h") and (atoms[k] == atoms[l]) for k, l in inds]
        # )
        # else:
        # is_hh = np.zeros(len(inds))

        aux_mask = np.logical_or(
            distances <= aux_below_thresh,
            distances <= (min_dist + aux_add_dist),
        )
        # aux_mask = distances <= 1.1 * min_dist
        # aux_mask = np.logical_and(aux_mask, ~is_hh)
        # Don't include interfragment bond
        aux_mask[min_ind] = False
        aux_inds = inds[aux_mask]

        if aux_keep >= 0:
            aux_dists = distances[aux_mask]
            aux_keep_inds = aux_dists.argsort()[:aux_keep]
            aux_inds = aux_inds[aux_keep_inds]

        aux_dict[(i, j)] = aux_inds
    frag_pairs = np.array(frag_pairs, dtype=int)

    """
    The code below tries to determine a reasonable distance, to define
    interfragment bonds. We don't want to use all interfragment bonds defined
    above, as this would connect all fragments to each other, resulting in an
    explosion of the number of internal coordinates.

    Here we try to cluster the interfragment distances, until one cluster center
    comes close to the scaled minimum interfragment distance. We then use this
    distance to filter all interfragment bonds.
    """
    if len(fragments) > 2:
        dists = np.reshape(interfrag_dists, (-1, 1))
        min_dist = dists.min()

        for n_clusters in range(2, 10):
            kmeans = KMeans(n_clusters=n_clusters)
            _ = kmeans.fit_predict(dists)
            min_center = kmeans.cluster_centers_.min()

            if min_center <= (min_scale * min_dist):
                break
        else:
            raise Exception("Not handled!")

        interfrag_inds = np.array(interfrag_inds)
        mask = dists <= min_scale * min_center
        # We also use interfragment bonds that still have a reasonable length.
        if (min_dist_thresh is not None) and (min_dist_thresh > 0.0):
            mask |= dists <= min_dist_thresh
        mask = mask.flatten()
        interfrag_inds = interfrag_inds[mask]
        conn_frags_mask = mask
    else:
        conn_frags_mask = np.ones(len(frag_pairs), dtype=bool)

    # Only keep auxiliary interfragment bonds between actually connected fragments
    actually_connected_frags = frag_pairs[conn_frags_mask]
    aux_interfrag_inds = list(
        it.chain(*[aux_dict[(i, j)] for i, j in actually_connected_frags])
    )
    aux_interfrag_inds = np.array(aux_interfrag_inds, dtype=int)
    return interfrag_inds, aux_interfrag_inds


def connect_fragments_ahlrichs(
    cdm,
    fragments,
    atoms,
    min_dist_scale=1.1,
    scale=1.2,
    avoid_h=False,
    avoid_hh=True,
    logger=None,
    max_dist=3.0 / BOHR2ANG,
):
    """Ahlrich/von Arnim connection scheme.

    See II. A in [2] and p. 2082 in [3]."""
    atoms = [atom.lower() for atom in atoms]
    if len(fragments) > 1:
        log(
            logger,
            f"Detected {len(fragments)} fragments. Generating interfragment bonds.",
        )
    dist_mat = squareform(cdm)

    # Tracks already connected fragments
    connected_fragments: Set[FrozenSet] = set()
    # Will contain interfragment atom pairs.
    interfrag_inds = list()
    h_sums = list()
    all_fragment_inds = [i for i, _ in enumerate(fragments)]
    # In the beginning, all fragments are unconnected
    unconnected_fragment_inds = all_fragment_inds.copy()
    h_inds = set([i for i, atom in enumerate(atoms) if atom == "h"])
    h_mask = np.array([1 if atom == "h" else 0 for atom in atoms])

    def get_pairs(unconnected, all_fragments):
        candidates = set(
            [
                frozenset((i, j))
                for i, j in it.product(unconnected, all_fragments)
                if i != j
            ]
        )
        return candidates - connected_fragments

    cycle = 0
    more_than_one_frag = len(fragments) > 1
    # Skip loop when only one fragment is present.
    # Loop until all fragments are somehow connected.
    while more_than_one_frag and len(unconnected_fragment_inds) > 0:
        candidate_pairs = get_pairs(unconnected_fragment_inds, all_fragment_inds)
        for i, j in candidate_pairs:
            frag1 = fragments[i]
            frag2 = fragments[j]
            # Pairs of possible interfragment bonds
            inds = np.array(
                [(i1, i2) for i1, i2 in it.product(frag1, frag2)], dtype=int
            )
            distances = np.array([dist_mat[k, l] for k, l in inds])

            # Determine minimum distance bond
            min_ind = distances.argmin()
            min_dist = distances[min_ind]

            # If we want to avoid interfragment bonds involving hydrogen we check
            # if 'min_ind' and the associated 'min_dist' must be updated to involve
            # a bond without hydrogen.
            # Alternatively, 'no_hh' can also be set to True; this only avoids
            # interfragment bonds between two hydrogens.
            if avoid_h:
                sort_inds = np.argsort(distances, kind="stable")
                for k, dist in zip(sort_inds, distances[sort_inds]):
                    if set(inds[k]) & h_inds:
                        continue
                    if dist >= 1.5 * min_dist:
                        break
                    min_ind = k
                    min_dist = distances[k]
                    break

            # If the minimum distances is below the current allowed maximum distance
            # we connect the fragments and try to define additional interfragment bonds
            # below or equal to 'min_dist + offset'.
            if min_dist <= max_dist:
                # Scaling of bigger 'min_dist' values by a fixed factor can lead
                # to definition of many additional bonds. We restrict the offset
                # to at most 1 Å.
                offset = min((min_dist_scale - 1.0) * min_dist, 1.0 / BOHR2ANG)
                mask = distances <= (min_dist + offset)
                cur_interfrag_inds = inds[mask]
                cur_h_sums = h_mask[inds[mask]].sum(axis=1)
                no_hh = cur_h_sums < 2
                # If fragment pairs is connected by any no-HH bond we drop all
                # HH-interfragment bonds.
                if avoid_hh and any(no_hh):
                    cur_interfrag_inds = cur_interfrag_inds[no_hh]
                    cur_h_sums = cur_h_sums[no_hh]
                interfrag_inds.extend(cur_interfrag_inds)
                h_sums.extend(cur_h_sums)
                # Indicate that the just connected fragments don't have to be
                # connected anymore. In the current cycle of the while loop additional
                # bonds to the just connected fragments can still be defined.
                unconnected_fragment_inds = [
                    k for k in unconnected_fragment_inds if k not in (i, j)
                ]
                connected_fragments.add(frozenset((i, j)))
                log(
                    logger,
                    f"\tMacro cycle {cycle}, connected fragments {i} ({len(frag1)} atom) "
                    f"and fragment {j} ({len(frag2)} atoms). Min. dist is "
                    f"{min_dist:.4f} au; current allowed max. dist. is {max_dist:.4f} au.",
                )
        # If there are still unconnected fragments present we allow longer
        # interfragment bonds.
        max_dist *= scale
        cycle += 1
        assert max_dist != math.inf, "Max. dist grew to infinity. Something went wrong!"

    interfrag_inds = np.array(interfrag_inds, dtype=int)
    return interfrag_inds, list()


def get_hydrogen_bond_inds(atoms, coords3d, bond_inds, logger=None):
    tmp_sets = [frozenset(bi) for bi in bond_inds]
    hydrogen_inds = [i for i, a in enumerate(atoms) if a.lower() == "h"]
    x_inds = [i for i, a in enumerate(atoms) if a.lower() in "n o f p s cl".split()]
    hydrogen_bond_inds = list()
    bond_sets = [set(map(int, bond)) for bond in bond_inds]
    for h_ind, x_ind in it.product(hydrogen_inds, x_inds):
        as_set = set((h_ind, x_ind))
        if as_set not in tmp_sets:
            continue
        # Check if distance of H to another electronegative atom Y is
        # greater than (1.1 * sum of their covalent radii) but smaller than
        # (0.9 * sum of their van der Waals radii). If the
        # angle X-H-Y is greater than 90° a hydrogen bond is asigned.
        y_inds = set(x_inds) - set((x_ind,))
        for y_ind in y_inds:
            y_atom = atoms[y_ind].lower()
            cov_rad_sum = CR["h"] + CR[y_atom]
            distance = Stretch._calculate(coords3d, (h_ind, y_ind))
            vdw_rad_sum = VDW_RADII["h"] + VDW_RADII[y_atom]
            angle = Bend._calculate(coords3d, (x_ind, h_ind, y_ind))
            if (
                (1.1 * cov_rad_sum < distance < 0.9 * vdw_rad_sum)
                and (angle > np.pi / 2)
                # Avoid adding hydrogen bonds that are already in bond_inds.
                # This can happen when the user manually defined a bond, that
                # is actually a hydrogen bond.
                and {h_ind, y_ind} not in bond_sets
            ):
                hydrogen_bond_inds.append((h_ind, y_ind))
                log(
                    logger,
                    f"Detected hydrogen bond between atoms {h_ind} "
                    f"({atoms[h_ind]}) and {y_ind} ({atoms[y_ind]})",
                )

    return hydrogen_bond_inds


def get_hydrogen_bond_inds_v2(atoms, coords3d, bond_inds, logger=None):
    def to_set(iterable):
        return {frozenset(_) for _ in iterable}

    atoms_lower = [atom.lower() for atom in atoms]
    # Determine Hydrogen indices
    org_h_inds = {i for i, atom in enumerate(atoms_lower) if atom == "h"}

    org_bond_sets = to_set(bond_inds)
    # Determine Hydrogen bonding partners in original bond set
    org_h_partners = dict()
    for org_bond in org_bond_sets:
        if h_set := org_bond & org_h_inds:
            (x_ind,) = org_bond - h_set
            (h_ind,) = h_set
            org_h_partners.setdefault(h_ind, list()).append(x_ind)

    hx = {"h", "n", "o", "f", "p", "s", "cl"}
    # Determine indices of potential hydrogen bond acceptors and hydrogen to
    # carry out a search for bonds with a bigger bond factor.
    hx_inds = [i for i, atom in enumerate(atoms) if atom.lower() in hx]
    hx_atoms = [atoms[i] for i in hx_inds]
    hx_map = {j: i for j, i in enumerate(hx_inds)}
    # See pysisyphus.elemdata.HBOND_FACTORS
    hx_coords3d = coords3d[hx_inds]
    # Search bonds with KDTree
    h_bonds = find_bonds_fast(hx_atoms, hx_coords3d, bond_factor=2.3)
    h_bonds = to_set(h_bonds)
    # Map back to original indices. Until now the indices were only in the basis
    # of the reduced number of atoms.
    h_bonds = {frozenset((hx_map[i], hx_map[j])) for i, j in h_bonds}
    # Drop h_bonds that were already defined as "normal" bonds
    h_bonds = h_bonds - org_bond_sets

    hydrogen_bond_inds = list()
    # h_bonds can still contain XX and HH bonds
    for h_bond in h_bonds:
        hi = h_bond & org_h_inds
        # Skip XX and HH. We are only interest in 'h_bond' with one hydrogen.
        if len(hi) in (0, 2):
            continue
        (h_ind,) = hi
        (y_partner,) = h_bond - hi
        v = coords3d[y_partner] - coords3d[h_ind]
        for x_partner in org_h_partners[h_ind]:
            u = coords3d[x_partner] - coords3d[h_ind]
            # > 90°
            if u.dot(v) < 0.0:
                hydrogen_bond_inds.append((h_ind, y_partner))
    return hydrogen_bond_inds


def get_bend_inds(coords3d, bond_inds, min_deg, max_deg, logger=None):
    bond_sets = {frozenset(bi) for bi in bond_inds}

    bend_inds = list()
    for bond_set1, bond_set2 in it.combinations(bond_sets, 2):
        union = bond_set1 | bond_set2
        if len(union) == 3:
            indices, _ = sort_by_central(bond_set1, bond_set2)
            if not bend_valid(coords3d, indices, min_deg, max_deg):
                log(logger, f"Bend {indices} is not valid!")
                continue
            bend_inds.append(indices)

    return bend_inds


def get_linear_bend_inds(coords3d, cbm, bends, min_deg=175, max_bonds=4, logger=None):
    linear_bends = list()
    complements = list()

    if min_deg is None:
        return linear_bends, complements

    bm = squareform(cbm)
    for bend in bends:
        deg = np.rad2deg(Bend._calculate(coords3d, bend))
        bonds = sum(bm[bend[1]])
        if (deg >= min_deg) and (bonds <= max_bonds):
            log(
                logger,
                f"Bend {bend}={deg:.1f}° is (close to) linear. "
                "Creating linear bend & complement.",
            )
            linear_bends.append(bend)
            complements.append(bend)
    return linear_bends, complements


def get_dihedral_inds(coords3d, bond_inds, bend_inds, max_deg, logger=None):
    max_rad = np.deg2rad(max_deg)
    bond_dict = dict()
    for from_, to_ in bond_inds:
        bond_dict.setdefault(from_, list()).append(to_)
        bond_dict.setdefault(to_, list()).append(from_)
    proper_dihedral_inds = list()
    improper_candidates = list()
    improper_dihedral_inds = list()

    def log_dihed_skip(inds):
        log(
            logger,
            f"Skipping generation of dihedral {inds} "
            "as some of the the atoms are (close too) linear.",
        )

    def set_dihedral_index(dihedral_ind, proper=True):
        dihed = tuple(dihedral_ind)
        check_in = proper_dihedral_inds if proper else improper_dihedral_inds
        # Check if this dihedral is already present
        if (dihed in check_in) or (dihed[::-1] in check_in):
            return
        # Assure that the angles are below 175° (3.054326 rad)
        if not dihedral_valid(coords3d, dihedral_ind, deg_thresh=max_deg):
            log_dihed_skip(dihedral_ind)
            return
        if proper:
            proper_dihedral_inds.append(dihed)
        else:
            improper_dihedral_inds.append(dihed)

    for bond, bend in it.product(bond_inds, bend_inds):
        central = bend[1]
        bend_set = set(bend)
        bond_set = set(bond)
        # Check if the two sets share one common atom. If not continue.
        intersect = bend_set & bond_set
        if len(intersect) != 1:
            continue

        # TODO: check collinearity of bond and bend.

        # When the common atom between bond and bend is a terminal, and not a central atom
        # in the bend we create a proper dihedral. Improper dihedrals are only created
        # when no proper dihedrals have been found.
        if central not in bond_set:
            # The new terminal atom in the dihedral is the one, that doesn' intersect.
            terminal = tuple(bond_set - intersect)[0]
            intersecting_atom = tuple(intersect)[0]
            bend_terminal = tuple(bend_set - {central} - intersect)[0]

            bend_rad = Bend._calculate(coords3d, bend)
            # Bend atoms are nearly collinear. Check if we can skip the central bend atom
            # and use an atom that is conneced to the terminal atom of the bend or bond.
            if bend_rad >= max_rad:
                bend_terminal_bonds = set(bond_dict[bend_terminal]) - bend_set
                bond_terminal_bonds = set(bond_dict[terminal]) - bond_set
                set_dihedrals = [
                    (terminal, intersecting_atom, bend_terminal, betb)
                    for betb in bend_terminal_bonds
                ] + [
                    (bend_terminal, intersecting_atom, terminal, botb)
                    for botb in bond_terminal_bonds
                ]
                # Hardcoded for now ... look ahead to next shell of atoms
                if not any(
                    [
                        dihedral_valid(coords3d, inds, deg_thresh=max_deg)
                        for inds in set_dihedrals
                    ]
                ):
                    set_dihedrals = []
                    for betb in bend_terminal_bonds:
                        bend_terminal_bonds_v2 = (
                            set(bond_dict[betb]) - bend_set - bond_set
                        )
                        set_dihedrals = [
                            (terminal, intersecting_atom, betb, betb_v2)
                            for betb_v2 in bend_terminal_bonds_v2
                        ]
                    for botb in bond_terminal_bonds:
                        bond_terminal_bonds_v2 = (
                            set(bond_dict[botb]) - bend_set - bond_set
                        )
                        set_dihedrals = [
                            (bend_terminal, intersecting_atom, botb, botb_v2)
                            for botb_v2 in bond_terminal_bonds_v2
                        ]
            elif intersecting_atom == bend[0]:
                set_dihedrals = [[terminal] + list(bend)]
            else:
                set_dihedrals = [list(bend) + [terminal]]
            [set_dihedral_index(dihed) for dihed in set_dihedrals]
        # If the common atom is the central atom we try to form an out
        # of plane bend / improper torsion. They may be created later on.
        else:
            fourth_atom = list(bond_set - intersect)
            dihedral_ind = list(bend) + fourth_atom
            # This way dihedrals may be generated that contain linear
            # atoms and these would be undefinied. So we check for this.
            if dihedral_valid(coords3d, dihedral_ind, deg_thresh=max_deg):
                improper_candidates.append(dihedral_ind)
            else:
                log_dihed_skip(dihedral_ind)

    # Now try to create the remaining improper dihedrals.
    if (len(coords3d) >= 4) and (len(proper_dihedral_inds) == 0):
        log(
            logger,
            "Could not define any proper dihedrals! Generating improper dihedrals!",
        )
        for improp in improper_candidates:
            set_dihedral_index(improp, proper=False)
        log(
            logger,
            "Permutational symmetry not considerd in generation of "
            "improper dihedrals.",
        )

    return proper_dihedral_inds, improper_dihedral_inds


def sort_by_prim_type(to_sort=None):
    if to_sort is None:
        to_sort = list()

    by_prim_type = [[], [], []]
    for item in to_sort:
        len_ = len(item)
        # len -> index
        #   2 ->     0 (bond)
        #   3 ->     1 (bend)
        #   4 ->     2 (torsion)
        by_prim_type[len_ - 2].append(tuple(item))
    return by_prim_type


CoordInfo = namedtuple(
    "CoordInfo",
    "bonds hydrogen_bonds interfrag_bonds aux_interfrag_bonds "
    "bends linear_bends linear_bend_complements "
    # "dihedrals typed_prims fragments cdm cbm".split(),
    "proper_dihedrals improper_dihedrals "
    "translation_inds rotation_inds cartesian_inds "
    "typed_prims fragments".split(),
)


def setup_redundant(
    atoms,
    coords3d,
    factor=BOND_FACTOR,
    define_prims=None,
    min_deg=BEND_MIN_DEG,
    dihed_max_deg=DIHED_MAX_DEG,
    lb_min_deg=None,
    lb_max_bonds=4,
    min_weight=None,
    tric=False,
    hybrid=False,
    interfrag_hbonds=True,
    hbond_angles=False,
    freeze_atoms=None,
    define_for=None,
    internals_with_frozen=False,
    rm_for_frag: Optional[set] = None,
    logger=None,
):
    if define_prims is None:
        define_prims = list()
    if freeze_atoms is None:
        freeze_atoms = list()
    if define_for is None:
        define_for = list()
    if rm_for_frag is None:
        rm_for_frag = set()

    log(
        logger,
        f"Detecting primitive internals for {len(atoms)} atoms.\n"
        f"Excluding {len(freeze_atoms)} frozen atoms from the internal coordinate setup.",
    )

    # Mask array. By default all atomes are used to generate internal coordinates.
    use_atoms = np.ones_like(atoms, dtype=bool)
    # Only use atoms in 'define_for' to generate internal coordinates
    if define_for:
        use_atoms[:] = False  # Disable/mask all others
        use_atoms[define_for] = True
    # If not explicitly enabled, don't form internal coordinates containing frozen atoms.
    # With 'internals_with_frozen', the bonds will be filtered for bonds, containing
    # at most one frozen atom.
    elif not internals_with_frozen:
        use_atoms[freeze_atoms] = False
    freeze_atom_set = set(freeze_atoms)
    atoms = [atom for mobile, atom in zip(use_atoms, atoms) if mobile]
    coords3d = coords3d[use_atoms]

    # Maps (different) indices of mobile atoms back to their original indices
    freeze_map = {
        sub_ind: org_ind for sub_ind, org_ind in enumerate(np.where(use_atoms)[0])
    }
    mobile_org_inds = set(freeze_map.values())

    def keep_coord(prim_cls, prim_inds):
        return (
            True
            if (min_weight is None)
            else (prim_cls._weight(atoms, coords3d, prim_inds, 0.12) >= min_weight)
        )

    def keep_coords(prims, prim_cls):
        return [prim for prim in prims if keep_coord(prim_cls, prim)]

    # Bonds
    bonds, cdm, cbm = get_bond_sets(
        atoms,
        coords3d,
        bond_factor=factor,
        return_cdm=True,
        return_cbm=True,
    )
    # Determine which bonds were added via define_prims. These will be removed
    # from the list and are already added here to the automatically defined bonds.
    defined_bonds = [(pt, *inds) for pt, *inds in define_prims if pt == PrimTypes.BOND]
    if defined_bonds:
        # Update define_prims and remove the bonds
        define_prims = [tp for tp in define_prims if tp not in defined_bonds]
        defined_bond_inds = [[from_, to_] for _, from_, to_ in defined_bonds]
        bonds = np.concatenate((bonds, defined_bond_inds), axis=0)
    if internals_with_frozen:
        bonds = [bond for bond in bonds if len(set(bond) & freeze_atom_set) <= 1]
    bonds = [tuple(bond) for bond in bonds]
    bonds = keep_coords(bonds, Stretch)
    bonds = [bond for bond in bonds if rm_for_frag.isdisjoint(set(bond))]

    # Fragments
    fragments = merge_sets(bonds) + [
        frozenset((rmed_atom,)) for rmed_atom in rm_for_frag
    ]
    # Check for unbonded single atoms and create fragments for them.
    bonded_set = set(tuple(np.ravel(bonds)))
    # Set of unbonded, single atoms
    unbonded_set = set(range(len(atoms))) - bonded_set - freeze_atom_set
    log(
        logger,
        f"Merging bonded atoms yielded {len(fragments)} fragment(s) and "
        f"{len(unbonded_set)} atoms.",
    )
    # Create an additional single atom set for all unbonded single atoms
    fragments.extend([frozenset((atom,)) for atom in unbonded_set])

    interfrag_bonds = list()
    aux_interfrag_bonds = list()
    translation_inds = list()
    rotation_inds = list()
    # With translational & rotational internal coordinates (TRIC) we don't need
    # interfragment coordinates.
    if tric:
        translation_inds = [list(fragment) for fragment in fragments]
        # Exclude rotational coordinates for atomic species (1 atom)
        rotation_inds = [list(fragment) for fragment in fragments if len(fragment) > 1]
    # Without TRIC we have to somehow connect all fragments.
    else:
        # interfrag_bonds, aux_interfrag_bonds = connect_fragments_kmeans(
        interfrag_bonds, aux_interfrag_bonds = connect_fragments_ahlrichs(
            cdm, fragments, atoms, logger=logger
        )

    # Hydrogen bonds
    assert interfrag_hbonds, "Disabling interfrag_hbonds is not yet supported!"
    hydrogen_bonds = get_hydrogen_bond_inds(atoms, coords3d, bonds, logger=logger)
    hydrogen_set = [frozenset(bond) for bond in hydrogen_bonds]

    def remove_h_bonds(bond_list):
        return [bond for bond in bond_list if set(bond) not in hydrogen_set]

    # Remove newly obtained hydrogen bonds from other lists
    interfrag_bonds = remove_h_bonds(interfrag_bonds)
    aux_interfrag_bonds = remove_h_bonds(aux_interfrag_bonds)
    bonds = remove_h_bonds(bonds)
    aux_bonds = list()  # Not defined by default
    # Don't use auxilary interfragment bonds for bend detection
    bonds_for_bends = [
        bonds,
    ]
    # If we use regular redundant internals (not TRIC) we define interfragment
    # bends.
    if not tric:
        bonds_for_bends += [interfrag_bonds]
        if hbond_angles:
            bonds_for_bends += [hydrogen_bonds]
    bonds_for_bends = set([frozenset(bond) for bond in it.chain(*bonds_for_bends)])

    # Bends
    bends = get_bend_inds(
        coords3d,
        bonds_for_bends,
        min_deg=min_deg,
        # All bends will be checked, for being (close to) linear and will be removed from
        # bend_inds, if needed. Thats why we keep 180° here.
        max_deg=180.0,
        logger=logger,
    )
    bends = keep_coords(bends, Bend)

    # Linear Bends and orthogonal complements
    linear_bends, linear_bend_complements = get_linear_bend_inds(
        coords3d,
        cbm,
        bends,
        min_deg=lb_min_deg,
        max_bonds=lb_max_bonds,
        logger=logger,
    )
    # Remove linear bends from bends
    bends = [bend for bend in bends if bend not in linear_bends]
    linear_bends = keep_coords(linear_bends, LinearBend)
    linear_bend_complements = keep_coords(linear_bend_complements, LinearBend)

    # Dihedrals
    bends_for_dihedrals = bends + linear_bends
    proper_dihedrals, improper_dihedrals = get_dihedral_inds(
        coords3d,
        bonds_for_bends,
        bends_for_dihedrals,
        max_deg=dihed_max_deg,
        logger=logger,
    )
    proper_dihedrals = keep_coords(proper_dihedrals, Torsion)
    improper_dihedrals = keep_coords(improper_dihedrals, Torsion)
    # Improper dihedrals are disabled for now in TRIC
    if tric:
        improper_dihedrals = []

    cartesian_inds = []
    if hybrid:
        cartesian_inds = [i for i, _ in enumerate(atoms)]

    """
    When additional primitives are given in 'define_prims', we want to append
    them to the correct lists, that may contain already some primitive internals
    detected by our algorithms. Here we define a map between the PrimTypes and the
    present lists.
    """
    defined_cartesians = list()
    defined_translations = list()
    defined_rotations = list()

    define_map = {
        PrimTypes.BOND: "bonds",
        PrimTypes.AUX_BOND: "aux_bonds",
        PrimTypes.HYDROGEN_BOND: "hydrogen_bonds",
        PrimTypes.INTERFRAG_BOND: "interfrag_bonds",
        PrimTypes.AUX_INTERFRAG_BOND: "aux_interfrag_bonds",
        PrimTypes.BEND: "bends",
        PrimTypes.LINEAR_BEND: "linear_bends",
        PrimTypes.LINEAR_BEND_COMPLEMENT: "linear_bend_complements",
        PrimTypes.PROPER_DIHEDRAL: "proper_dihedrals",
        PrimTypes.IMPROPER_DIHEDRAL: "improper_dihedrals",
    }
    unmapped_typed_prims = list()
    for type_, *indices in define_prims:
        try:
            key = define_map[type_]
        except KeyError:
            try:
                key = define_map[PrimTypes(type_)]
            except KeyError:
                """
                With the current approach, some primitives in 'define_prims' can't
                be mapped to their respective lists, e.g., given CARTESIAN_X/Y/Z.
                Currently, every item in 'cartesian_inds' is expanded to three
                Cartesians (X/Y/Z). When only the X-component is to be defined, adding
                only the atom index to the list would result in all three components
                to be generated.

                So instead of adding them to their respective lists we keep them in
                'unmapped_typed_prims' and use them later as is.
                """
                unmapped_typed_prims.append((type_, *indices))
                continue
        locals()[key].append(tuple(indices))

    def make_tp(prim_type, *indices):
        """Map possibly modified indices to their original indices.

        With frozen atoms, the indices used to set up internal coordinates do
        not correspond to the actual indices. Here we map them back.
        """
        try:
            org_indices = [freeze_map[ind] for ind in indices]
            """The given 'indices' may not be present in freeze_map. This can happen,
            when coordinates are rebuilt, frozen atoms are excluded from using them
            in the coordinate definition and a previous set of coordinates ALREADY
            using the original indices is supplied."""
        except KeyError as error:
            """In such a case we check if the given coordinates are already fully
            defined in terms of original indices. If so, we use them as is."""
            if set(indices) < mobile_org_inds:
                org_indices = indices
            else:
                raise error
        # Convert from numpy.int64 to normal int
        org_indices = tuple(map(int, org_indices))
        return (prim_type, *org_indices)

    # Shortcut for PrimTypes Enum
    pt = PrimTypes
    # Create actual typed prims with the desired indices
    typed_prims = (
        # Bonds, two indices
        [make_tp(pt.BOND, *bond) for bond in bonds]
        + [make_tp(pt.AUX_BOND, *abond) for abond in aux_bonds]
        + [make_tp(pt.HYDROGEN_BOND, *hbond) for hbond in hydrogen_bonds]
        + [make_tp(pt.INTERFRAG_BOND, *ifbond) for ifbond in interfrag_bonds]
        + [make_tp(pt.AUX_INTERFRAG_BOND, *aifbond) for aifbond in aux_interfrag_bonds]
        # Bends, three indices
        + [make_tp(pt.BEND, *bend) for bend in bends]
        + [make_tp(pt.LINEAR_BEND, *lbend) for lbend in linear_bends]
        + [
            make_tp(pt.LINEAR_BEND_COMPLEMENT, *lbendc)
            for lbendc in linear_bend_complements
        ]
        # Dihedral, four indices
        + [make_tp(pt.PROPER_DIHEDRAL, *pdihedral) for pdihedral in proper_dihedrals]
        + [
            make_tp(pt.IMPROPER_DIHEDRAL, *idihedral)
            for idihedral in improper_dihedrals
        ]
        + [make_tp(pt.CARTESIAN_X, cind) for cind in cartesian_inds]
        + [make_tp(pt.CARTESIAN_Y, cind) for cind in cartesian_inds]
        + [make_tp(pt.CARTESIAN_Z, cind) for cind in cartesian_inds]
        + [make_tp(pt.CARTESIAN_Z, cind) for cind in cartesian_inds]
    )

    # Translational and rotational coordinates result in 3 different coordinates each
    for frag in translation_inds:
        typed_prims += [
            make_tp(pt.TRANSLATION_X, *frag),
            make_tp(pt.TRANSLATION_Y, *frag),
            make_tp(pt.TRANSLATION_Z, *frag),
        ]
    for frag in rotation_inds:
        typed_prims += [
            make_tp(pt.ROTATION_A, *frag),
            make_tp(pt.ROTATION_B, *frag),
            make_tp(pt.ROTATION_C, *frag),
        ]
    typed_prims += unmapped_typed_prims
    # Drop duplicated typed_prims
    typed_prims = tuple(dict.fromkeys(typed_prims))

    coord_info = CoordInfo(
        bonds=bonds,
        hydrogen_bonds=hydrogen_bonds,
        interfrag_bonds=interfrag_bonds,
        aux_interfrag_bonds=aux_interfrag_bonds,
        bends=bends,
        linear_bends=linear_bends,
        linear_bend_complements=linear_bend_complements,
        proper_dihedrals=proper_dihedrals,
        improper_dihedrals=improper_dihedrals,
        translation_inds=translation_inds,
        rotation_inds=rotation_inds,
        cartesian_inds=cartesian_inds,
        typed_prims=typed_prims,
        fragments=fragments,
    )
    return coord_info


def setup_redundant_from_geom(geom, *args, **kwargs):
    return setup_redundant(geom.atoms, geom.coords3d, *args, **kwargs)


def get_primitives(coords3d, typed_prims, logger=None):
    primitives = list()
    for type_, *indices in typed_prims:
        cls = PrimMap[type_]
        cls_kwargs = {"indices": indices}
        if type_ in Rotations:
            cls_kwargs["ref_coords3d"] = coords3d
        primitives.append(cls(**cls_kwargs))

    msg = (
        "Defined primitives\n"
        + "\n".join(
            [f"\t{i:03d}: {str(p.indices): >14}" for i, p in enumerate(primitives)]
        )
        + "\n"
    )
    log(logger, msg)
    return primitives
