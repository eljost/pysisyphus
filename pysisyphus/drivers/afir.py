from functools import reduce
from pathlib import Path

import numpy as np

from pysisyphus.calculators.AFIR import AFIR, AFIRPath
from pysisyphus.config import OUT_DIR_DEFAULT
from pysisyphus.cos.NEB import NEB
from pysisyphus.drivers.opt import run_opt
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import pick_image_inds
from pysisyphus.intcoords.helpers import get_bond_difference


def generate_random_union(geoms, offset=1.0, copy=True):
    assert 2 <= len(geoms) <= 6
    # Center, rotate and displace from origin acoording to approximate radius
    # and an offset.
    # Displace along +x, -x, +y, -y, +z, -z. Affords at max 6 fragments.
    axis_inds = (0, 0, 1, 1, 2, 2)
    randomized = list()
    for i, geom in enumerate(geoms):
        if copy:
            geom = geom.copy()
        geom.center()
        geom.rotate()
        step_size = geom.approximate_radius() + offset
        step = np.zeros(3)
        axis = axis_inds[i]
        # Alternate between negative and positive direction along x/y/z
        step[axis] = (-1) ** i * step_size
        geom.coords3d += step[None, :]
        randomized.append(geom)
    union = reduce(lambda geom1, geom2: geom1 + geom2, randomized)
    return union


def relax_afir_path(atoms, cart_coords, calc_getter, images=15, out_dir=None):
    image_inds = pick_image_inds(cart_coords, images=images)
    images = [Geometry(atoms, cart_coords[i]) for i in image_inds]

    # Relax last image
    opt_kwargs = {
        "dump": True,
        "prefix": "last",
        "out_dir": out_dir,
    }
    last_image = images[-1]
    last_image_backup = last_image.copy()
    run_opt(last_image, calc_getter, opt_key="rfo", opt_kwargs=opt_kwargs)
    _, broken = get_bond_difference(last_image, last_image_backup)
    if broken:
        return

    cos_kwargs = {}
    cos = NEB(images, **cos_kwargs)
    cos_opt_kwargs = {
        "align": True,
        "dump": True,
        "max_cycles": 30,
        "out_dir": out_dir,
    }
    run_opt(cos, calc_getter, opt_key="lbfgs", opt_kwargs=cos_opt_kwargs)


def run_mc_afir_paths(geoms, calc_getter, num=5, gamma=None, t=None, T=298.15):
    assert gamma or t, "Either parameter gamme of time t must be given!"

    # Set up list of fragments
    i = 0
    fragments = list()
    for geom in geoms:
        atom_num = len(geom.atoms)
        fragments.append(np.arange(atom_num) + i)
        i += atom_num

    unions = [generate_random_union(geoms, copy=True) for _ in range(num)]
    # with open("unions.trj", "w") as handle:
    # handle.write("\n".join([union.as_xyz() for union in unions]))

    afir_paths = list()
    for i, union in enumerate(unions):
        out_dir = Path(f"out_{i:03d}")

        union_backup = union.copy()
        # Preoptimize union, so we start from an actual minimum
        union_opt_kwargs = {
            "prefix": "union",
            "out_dir": out_dir,
        }
        run_opt(union, calc_getter, opt_key="rfo", opt_kwargs=union_opt_kwargs)

        afir_path = multicomponent_afir(
            union, calc_getter, fragments, gamma=gamma, out_dir=out_dir
        )
        afir_paths.append(afir_path)
        formed, broken = get_bond_difference(union, union_backup)
        if formed or broken:
            relax_afir_path(
                union.atoms, afir_path.cart_coords, calc_getter, out_dir=out_dir
            )


def multicomponent_afir(geom, calc_getter, fragments, gamma, out_dir=None):
    actual_calc = calc_getter(out_dir=out_dir / OUT_DIR_DEFAULT)

    def afir_calc_getter():
        afir_calc = AFIR(
            actual_calc, fragment_indices=fragments, gamma=gamma, out_dir=out_dir
        )
        return afir_calc

    opt_kwargs = {
        "dump": True,
        "out_dir": out_dir,
        "prefix": "afir",
        "max_cycles": 200,
    }
    opt_result = run_opt(geom, afir_calc_getter, opt_key="rfo", opt_kwargs=opt_kwargs)
    opt = opt_result.opt

    afir_path = AFIRPath(
        atoms=geom.atoms,
        cart_coords=np.array(opt.cart_coords),
        energies=np.array(opt.true_energies),
        forces=np.array(opt.true_forces),
        opt_is_converged=opt.is_converged,
    )

    return afir_path


def analyze_afir_path(energies):
    energies = np.array(energies)
    # min_ind = energies.argmin()
    max_ind = energies.argmax()

    local_minima = list()
    local_maxima = list()
    stationary_points = list()
    for i, en in enumerate(energies[1:-1], 1):
        prev_en = energies[i - 1]
        next_en = energies[i + 1]
        if is_minimum := prev_en > en < next_en:
            local_minima.append(i)
        if is_ts := prev_en < en > next_en:
            local_maxima.append(i)
        if is_minimum or is_ts:
            stationary_points.append(i)

    if max_ind in stationary_points:
        ts_ind = np.where(stationary_points == max_ind)[0][0]
        prev_min_ind = stationary_points[ts_ind - 1]
        next_min_ind = stationary_points[ts_ind + 1]
        sp_inds = [prev_min_ind, max_ind, next_min_ind]
    else:
        sp_inds = list()
    return sp_inds
