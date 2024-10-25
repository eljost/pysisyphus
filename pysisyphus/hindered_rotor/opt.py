import functools
from pathlib import Path

import numpy as np

from pysisyphus.drivers import opt
from pysisyphus.hindered_rotor import fragment as hr_fragment
from pysisyphus.Geometry import Geometry


# Define a custom exception for an empty store
class EmptyStoreException(Exception):
    pass


def find_closest_key(rad_store, target_radian):
    # Raise an exception if the dictionary is empty
    if not rad_store:
        raise EmptyStoreException("The radian dictionary is empty.")

    # Normalize the target radian value to the range [0, 2*pi]
    target_radian = target_radian % (2 * np.pi)

    # Extract keys and calculate distances
    target_radians = np.array(list(rad_store.keys()))
    distances = np.abs(target_radians - target_radian)
    distances = np.minimum(distances, 2 * np.pi - distances)

    # Report calculated distances
    print(f"Target radian: {target_radian}")
    # print(f"Calculated distances: {distances}")

    # Select the minimum distance
    closest_index = np.argmin(distances)
    closest_key = target_radians[closest_index]
    print(f"Closest key: {closest_key}")

    return closest_key


def opt_closure(geom, torsion, calc_getter, opt_kwargs=None, out_dir=Path(".")):
    _opt_kwargs = {
        "type": "rfo",
        "thresh": "gau",
    }
    if opt_kwargs is None:
        opt_kwargs = dict()
    _opt_kwargs.update(opt_kwargs)
    opt_kwargs = _opt_kwargs
    opt_cls = opt_kwargs.pop("type")

    atoms = geom.atoms
    _geom = geom.copy(coord_type="redund")
    typed_prims = _geom.internal.typed_prims
    bond = torsion[1:3]
    tp_constraint = ("PROPER_DIHEDRAL", *torsion)
    tp_bond = ("BOND", *bond)
    calc_number = 0

    # This dict will store optimized coordinates and associated energies
    rad_store = dict()

    _, fragment = hr_fragment.fragment_geom(geom, bond)

    def energy_getter(rad):
        nonlocal calc_number

        rad_key = f"{rad:5.3f}_rad"

        # It is beneficial to start a new optimization from the geometry, that
        # is closest to the desired radian value. Below, we try to look up the
        # closest geometry.
        try:
            rad_trial = find_closest_key(rad_store, rad)
            _, c3d_trial = rad_store[rad_trial]
        except EmptyStoreException:
            c3d_trial = geom.coords3d.copy()
            rad_trial = 0.0

        # When we found a geometry that is closer to the desired radian value
        # we have to correct the radian for rotation by the radian value of the
        # stored geometry.
        target_deg = np.rad2deg(rad - rad_trial)
        c3d_rot = hr_fragment.rotate_fragment_around_bond(
            c3d_trial, bond, fragment, target_deg
        )
        opt_geom = Geometry(
            atoms,
            c3d_rot,
            coord_type="redund",
            coord_kwargs={
                # TODO: Using the same/trying to use the same typed_prims is probably
                # neither required nor sensible. Check if this can or should be disabled.
                "typed_prims": typed_prims,
                #
                "constrain_prims": [tp_constraint],
                # Especially for systems that are initially not bonded it is beneficial
                # to also define the central bond of the torsion, so additional coordinates
                # are defined.
                "define_prims": [tp_bond, tp_constraint],
            },
        )
        # Set appropriate calc_number to distinguish the log
        calc_getter_wrapped = functools.partial(calc_getter, calc_number=calc_number)
        title = f"Rad {rad:5.3f} Optimization"
        opt_result = opt.run_opt(
            opt_geom, calc_getter_wrapped, opt_cls, opt_kwargs, title=title
        )
        # Continue with (potentially) replaced optimized geometry
        opt_geom = opt_result.geom
        # Dump optimized geometry
        opt_geom.dump_xyz(out_dir / f"{calc_number:03d}_opt_{rad_key}.xyz")
        coords3d_opt = opt_geom.coords3d
        energy = opt_geom.energy
        # Store optimize geometry so it can later be utilized to start new optimizations
        # with similar radian.
        rad_store[rad] = (energy, coords3d_opt.copy())
        calc_number += 1
        return energy

    return energy_getter, rad_store
