# [1] https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.23910
#     Birkholz, Schlegel, 2015

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.drivers import relaxed_scan
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.intcoords.setup import get_bond_mat
from pysisyphus.intcoords.PrimTypes import normalize_prim_inputs


def bond_order(r, r0, b=2):
    """Bond order for given bond length and reference length.

    Eq. (3) in [1]."""
    return max(0, (b * (r0 / r) - 1) / (b - 1))


def bond_orders(coords3d, bond_indices, r0s):
    """List of bond orders."""
    bos = list()
    for r0, (a, b) in zip(r0s, bond_indices):
        r = np.linalg.norm(coords3d[a] - coords3d[b])
        bo = bond_order(r, r0)
        bos.append(bo)
    return np.array(bos)


def get_r0s(geom, bond_indices):
    """Reference bond lengths as sum of covalent radii."""
    crs = geom.covalent_radii
    a, b = bond_indices.T
    return crs[a] + crs[b]


def bond_orders_for_geom(geom, bond_indices):
    """Wrapper for bond_orders for simple use with Geometry."""
    r0s = get_r0s(geom, bond_indices)
    return bond_orders(geom.coords3d, bond_indices, r0s)


def length_for_bond_order(bo, r0, b=2):
    """Return bond length for given bond order and reference length.

    Eq. (3) in [1]."""
    return b / ((b - 1) * bo + 1) * r0


def birkholz_interpolation(geoms, calc_getter, recreate=True):
    assert len(geoms) >= 2
    start = geoms[0]
    end = geoms[-1]
    geoms = (start, end)
    atoms = start.atoms
    print(f"Coordinates of 'start' geometry\n{start.as_xyz()})")
    print(f"Coordinates of 'end' geometry\n{end.as_xyz()})\n")

    print(highlight_text("Bonding analysis"))
    # Bond matrices
    bm_start, bm_end = [get_bond_mat(geom) for geom in geoms]
    # Determine formed/broken bonds. The lower triangular part of the bond
    # matrices has to be zeroed, to avoid double counting. Changes in bonding
    # are determined with logical XOR.
    #
    # Start  End    XOR    Comment
    # -----  ---    ---    -------
    # True   True   False  Bond is present at both geometries.
    # False  False  False  Bond is absent at both geometries.
    # True   False  True   Bond is broken when going from start to end.
    # False  True   True   Bond is formed when going from start to end.
    bm_diff = np.triu(bm_start) ^ np.triu(bm_end)
    bond_indices = np.stack(np.nonzero(bm_diff), axis=1)
    bonding_changes = bm_diff.sum()
    print(f"{bonding_changes} bonds are formed/broken when going from start to end.\n")

    # Determine bond orders at start and end
    all_bos = [bond_orders_for_geom(geom, bond_indices) for geom in geoms]
    # Reference bond lengths
    r0s = get_r0s(start, bond_indices)

    # Print summary
    for title, geom, bos in zip(("Start", "End"), geoms, all_bos):
        print(highlight_text(title, level=1))
        for bo_, (a, b) in zip(bos, bond_indices):
            aa = atoms[a]
            ba = atoms[b]
            print(f"Bond ({aa}{a: >4},{ba}{b: >4}): BO={bo_:.2f}")
        print(
            f"mean(BOs)={bos.mean():.2f}, max(BOs)={bos.max():.2f}, min(BOs)={bos.min():.2f}"
        )
        print()

    # Guess goal bond lenghts from mean bond order
    start_bos, end_bos = all_bos
    mean_bos = (start_bos + end_bos) / 2
    # Determine goal bond lenghts that yield the mean bond orders
    goal_lengths = length_for_bond_order(mean_bos, r0s)

    # Constrain formed/broken bonds
    constrain_prims = normalize_prim_inputs([["BOND", *bond] for bond in bond_indices])

    # Start relaxed scans towards TS from 'start' and 'end'
    print(highlight_text("Relaxed scan towards TS, from 'start'"))
    start_guess, *_ = relaxed_scan(
        start, calc_getter, constrain_prims, target_values=goal_lengths, title="start"
    )
    print()
    start_guess.dump_xyz("start_guess")

    print(highlight_text("Relaxed scan towards TS, from 'end'"))
    end_guess, *_ = relaxed_scan(
        end, calc_getter, constrain_prims, target_values=goal_lengths, title="end"
    )
    print()
    end_guess.dump_xyz("end_guess")

    # Compare energies and take the one with lower energy
    start_en = start_guess.energy
    end_en = end_guess.energy
    print(f"Energy of 'start'-guess: {start_en:.6f} au")
    print(f"Energy of   'end'-guess: {end_en:.6f} au")
    start_lower = start_en <= end_en
    lower = "'start'-guess" if start_lower else "'end'-guess"
    print(f"{lower} has lower energy.")
    ts_guess = start_guess if start_lower else end_guess

    if recreate:
        print("Recreated Geometry object.")
        ts_guess = Geometry(ts_guess.atoms, ts_guess.cart_coords, coord_type="redund")
        ts_guess.set_calculator(calc_getter())
    print()
    print(f"Coordinates of TS guess\n{ts_guess.as_xyz()})\n")
    ts_guess_fn = "ts_guess.xyz"
    ts_guess.dump_xyz(ts_guess_fn)
    print(f"Dumped TS-guess to '{ts_guess_fn}'.")

    return ts_guess
