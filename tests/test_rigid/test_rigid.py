import numpy as np
import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.rigid import rot_around_bond, get_stepper


d2r = np.deg2rad


@pytest.fixture
def geom():
    return geom_loader("lib:methane_iminium_cation.xyz")


@pytest.mark.parametrize(
    "bond, out_fn", (
        ([0, 1], "rot.xyz"),
        ([0, 4], "rot_term.xyz"),
        ([4, 5], "no_bond.xyz"),
        ([0, 5], "one_frag.xyz"),
    )
)
def test_rot_vec(bond, out_fn, this_dir, geom):
    rad = d2r(90.0)
    geom.coords3d = rot_around_bond(geom, bond, rad)
    geom.dump_xyz(this_dir / out_fn)


@pytest.mark.parametrize(
    "key, typed_prim, step_size, nsteps",
    (
        ("trans", ("BOND", 0, 1), 0.333, 3),
        ("bend_rot", ("BEND", 4, 0, 1), d2r(10.0), 3),
        ("torsion", ("TORSION", 4, 0, 1, 3), d2r(10.0), 5),
        ("rot_bond", ("ROT_BOND", 0, 1), d2r(10.0), 5),
    ),
)
def test_stepper(key, typed_prim, step_size, nsteps, geom, this_dir):
    stepper = get_stepper(geom.atoms, geom.coords3d, typed_prim, step_size, nsteps)

    xyzs = [geom.as_xyz(cart_coords=c3d) for _, c3d in stepper()]
    fn = f"{key}_stepper.trj"
    with open(this_dir / fn, "w") as handle:
        handle.write("\n".join(xyzs))
    assert len(xyzs) == nsteps + 1


@pytest.mark.parametrize(
    "key, typed_prim, step_size, nsteps",
    (
        ("torsion", ("TORSION", 2, 0, 1, 3), d2r(10.0), 9),
        ("rot_bond", ("ROT_BOND", 0, 1), d2r(10.0), 9),
    ),
)
def test_cycle_stepper(key, typed_prim, step_size, nsteps, this_dir):
    fn = "lib:cycle_geom.xyz"
    geom = geom_loader(fn)
    stepper = get_stepper(geom.atoms, geom.coords3d, typed_prim, step_size, nsteps)

    xyzs = [geom.as_xyz(cart_coords=c3d) for _, c3d in stepper()]
    fn = f"cycle_{key}_stepper.trj"
    with open(this_dir / fn, "w") as handle:
        handle.write("\n".join(xyzs))
    assert len(xyzs) == nsteps + 1
