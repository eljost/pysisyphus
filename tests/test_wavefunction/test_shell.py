from pathlib import Path

import numpy as np
import pytest

from pysisyphus.config import WF_LIB_DIR
from pysisyphus.io.cube import Cube, write_cube
from pysisyphus.wavefunction import Wavefunction


@pytest.mark.parametrize(
    "fn",
    (WF_LIB_DIR / "orca_ch4_def2svp.json",),
)
def test_cgto_eval(fn):
    wf = Wavefunction.from_orca_json(fn)
    n_points = 5
    xyz = np.random.rand(n_points, 3)
    shells = wf.shells
    ref_vals = np.array([shells.eval_single(r) for r in xyz])
    vals = shells.eval(xyz)
    np.testing.assert_allclose(vals, ref_vals)


def get_grid(coords3d, num=10, offset=3.0):
    minx, miny, minz = coords3d.min(axis=0) - offset
    maxx, maxy, maxz = coords3d.max(axis=0) + offset
    # minx = miny = minz = -6.0
    # maxx = maxy = maxz =  6.0
    X, Y, Z = np.mgrid[
        minx : maxx : num * 1j,
        miny : maxz : num * 1j,
        minz : maxz : num * 1j,
    ]
    xyz = np.stack((X.flatten(), Y.flatten(), Z.flatten()), axis=1)
    spacing = np.array((maxx - minx, maxy - miny, maxz - minz)) / (num - 1)
    return xyz, spacing


@pytest.mark.parametrize(
    "fn",
    (
        # WF_LIB_DIR / "orca_ch4_def2svp.json",
        # WF_LIB_DIR / "orca_h2o_def2svp.json",
        "bz/00_benzene.json",
    ),
)
def test_grid(fn):
    wf = Wavefunction.from_orca_json(fn)
    Pa, Pb = wf.P
    num = 50
    xyz, spacing = get_grid(wf.coords3d, num=num)
    vals = wf.shells.eval(xyz, sph=True)
    rho_a = 2 * np.einsum("uv,iu,iv->i", Pa, vals, vals)
    # rho_b = np.einsum("uv,iu,iv->i", Pb, vals, vals)

    from skimage import measure
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    iso_val = 0.1
    verts, faces, *_ = measure.marching_cubes(
        rho_a.reshape(num, num, num), iso_val, spacing=spacing
    )
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], cmap="Spectral", lw=1)
    ax.set_zlim(-4, 4)
    plt.show()


@pytest.mark.parametrize(
    "fn",
    (
        "bz/00_benzene.json",
        "h/00_h.json",
    ),
)
def test_density_cube(fn):
    stem = Path(fn).stem
    cube_fn = f"{stem}.cube"
    wf = Wavefunction.from_orca_json(fn)
    Pa, Pb = wf.P
    num = 50
    xyz, spacing = get_grid(wf.coords3d, num=num)
    axes = np.diag(spacing)
    vals = wf.shells.eval(xyz, sph=True)
    rho_a = np.einsum("uv,iu,iv->i", Pa, vals, vals)
    rho_b = np.einsum("uv,iu,iv->i", Pb, vals, vals)
    rho = rho_a + rho_b
    dV = np.product(spacing)
    print(f"rho.sum() * dV = {rho.sum() * dV:.8f}")
    origin = xyz[0]
    vol_data = rho_a.reshape(num, num, num)
    cube = write_cube(
        atoms=wf.atoms,
        coords3d=wf.coords3d,
        vol_data=vol_data,
        origin=origin,
        axes=axes,
    )
    # with open("benzene.cube", "w") as handle:
    with open(cube_fn, "w") as handle:
        handle.write(cube)
    print("wrote cube_fn", cube_fn)


def get_grid2(origin, npoints, axes):
    minx, miny, minz = origin
    nx, ny, nz = npoints
    # Only support diagonal axes matrix for now
    # Reconstruct matrix from diagonal
    dx, dy, dz = spacing = np.diag(axes)
    np.testing.assert_allclose(np.diag(spacing), axes)

    maxx, maxy, maxz = origin + spacing * (npoints - 1)
    X, Y, Z = np.mgrid[
        minx : maxx : nx * 1j,
        miny : maxz : ny * 1j,
        minz : maxz : nz * 1j,
    ]
    xyz = np.stack((X.flatten(), Y.flatten(), Z.flatten()), axis=1)
    return xyz, spacing


@pytest.mark.parametrize(
    "fn, ref_fn",
    (
        # ("bz/00_benzene.json", "bz/density.cub"),
        ("h/00_h.json", "h/density.cub"),
        ("h/00_h.molden.input", "h/density.cub"),
        ("h/01_h_1p.json", "h/1p_density.cub"),
        ("ch2/00_ch2.json", "ch2/density.cub"),
        ("ch2/00_ch2.molden.input", "ch2/density.cub"),
        ("c/00_c.molden.input", "c/00_density.cub"),
        ("c/01_c.molden.input", "c/01_density.cub"),
    ),
)
def test_mwfn_comp(fn, ref_fn):
    ref_cube = Cube.from_file(ref_fn)
    grid, spacing = get_grid2(ref_cube.origin, ref_cube.npoints, ref_cube.axes)

    stem = Path(fn).stem
    cube_fn = f"{stem}.cube"
    if fn.endswith(".json"):
        wf = Wavefunction.from_orca_json(fn)
    elif fn.endswith(".molden.input"):
        wf = Wavefunction.from_orca_molden(fn)
    else:
        raise Exception("Whoot?")

    Pa, Pb = wf.P
    vals = wf.shells.eval(grid, sph=True)
    rho_a = np.einsum("uv,iu,iv->i", Pa, vals, vals)
    rho_b = np.einsum("uv,iu,iv->i", Pb, vals, vals)
    rho = rho_a + rho_b
    dV = np.product(spacing)
    print(f"rho.sum() * dV = {rho.sum() * dV:.8f}")
    origin = ref_cube.origin
    vol_data = rho.reshape(*ref_cube.npoints)
    cube = write_cube(
        atoms=wf.atoms,
        coords3d=wf.coords3d,
        vol_data=vol_data,
        origin=origin,
        axes=ref_cube.axes,
    )
    with open(cube_fn, "w") as handle:
        handle.write(cube)
    print("wrote cube_fn", cube_fn)
    np.testing.assert_allclose(vol_data, ref_cube.vol_data)
