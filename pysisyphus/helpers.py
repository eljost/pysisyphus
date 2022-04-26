from collections import namedtuple
import getpass
import itertools as it
import logging
from math import log
import os
from pathlib import Path
import re
import sys
import time

import numpy as np
import scipy as sp
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from pysisyphus.config import p_DEFAULT, T_DEFAULT, LIB_DIR
from pysisyphus.constants import ANG2BOHR, AU2KJPERMOL
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import (
    eigval_to_wavenumber,
    report_isotopes,
    highlight_text,
)
from pysisyphus.io import (
    geom_from_cjson,
    geom_from_crd,
    geom_from_mol2,
    geom_from_pdb,
    save_hessian as save_h5_hessian,
    geom_from_zmat_fn,
    geoms_from_inline_xyz,
    geom_from_pubchem_name,
)
from pysisyphus.thermo import (
    can_thermoanalysis,
    print_thermoanalysis,
)
from pysisyphus.xyzloader import parse_xyz_file, parse_trj_file, make_trj_str


def geom_from_xyz_file(xyz_fn, coord_type="cart", **coord_kwargs):
    kwargs = {
        "coord_type": coord_type,
    }
    kwargs.update(coord_kwargs)
    xyz_fn = str(xyz_fn)
    atoms, coords, comment = parse_xyz_file(xyz_fn, with_comment=True)
    coords *= ANG2BOHR
    geom = Geometry(
        atoms,
        coords.flatten(),
        comment=comment,
        **kwargs,
    )
    return geom


def geoms_from_trj(trj_fn, first=None, coord_type="cart", **coord_kwargs):
    trj_fn = str(trj_fn)
    kwargs = {
        "coord_type": coord_type,
    }
    kwargs.update(coord_kwargs)
    atoms_coords_comments = parse_trj_file(trj_fn, with_comments=True)[:first]
    geoms = [
        Geometry(atoms, coords.flatten() * ANG2BOHR, comment=comment, **kwargs)
        for atoms, coords, comment in atoms_coords_comments
    ]
    return geoms


def geom_loader(fn, coord_type="cart", iterable=False, **coord_kwargs):
    """After introducing the pubchem functionality I don't like this
    function anymore :) Too complicated."""
    fn = str(fn)
    org_fn = fn

    split_ = re.split(r"\[(-?\d+)\]$", fn)
    fn = split_.pop(0)
    if split_:
        index = int(split_.pop(0))
    else:
        index = None
    ext = "" if "\n" in fn else Path(fn).suffix

    funcs = {
        ".mol2": geom_from_mol2,
        ".crd": geom_from_crd,
        ".xyz": geom_from_xyz_file,
        ".trj": geoms_from_trj,
        ".pdb": geom_from_pdb,
        ".cjson": geom_from_cjson,
        ".zmat": geom_from_zmat_fn,
        "": geoms_from_inline_xyz,
    }
    assert ext in funcs, f"Unknown filetype for '{fn}'!"
    func = funcs[ext]

    if fn.startswith("lib:"):
        fn = str(LIB_DIR / fn[4:])
    elif fn.startswith("pubchem:"):
        func = geom_from_pubchem_name
        fn = fn[8:]

    kwargs = {
        "coord_type": coord_type,
    }
    kwargs.update(coord_kwargs)
    geom = func(fn, **kwargs)

    if index is not None:
        geom = geom[index]

    if iterable and org_fn.startswith("pubchem:"):
        geom = (geom,)
    if iterable and (ext in (".trj", "")) and index is None:
        geom = tuple(geom)
    elif iterable:
        geom = (geom,)

    return geom


def align_geoms(geoms):
    # http://nghiaho.com/?page_id=671#comment-559906
    first_geom = geoms[0]
    coords3d = first_geom.coords3d
    centroid = coords3d.mean(axis=0)
    last_centered = coords3d - centroid
    first_geom.coords3d = last_centered
    atoms_per_image = len(first_geom.atoms)

    # Don't rotate the first image, so just add identitiy matrices
    # for every atom.
    rot_mats = [np.eye(3)] * atoms_per_image
    for i, geom in enumerate(geoms[1:], 1):
        coords3d = geom.coords3d
        centroid = coords3d.mean(axis=0)
        # Center next image
        centered = coords3d - centroid
        tmp_mat = centered.T.dot(last_centered)
        U, W, Vt = np.linalg.svd(tmp_mat)
        rot_mat = U.dot(Vt)
        # Avoid reflections
        if np.linalg.det(rot_mat) < 0:
            U[:, -1] *= -1
            rot_mat = U.dot(Vt)
        # Rotate the coords
        rotated3d = centered.dot(rot_mat)
        geom.coords3d = rotated3d
        last_centered = rotated3d
        # Append one rotation matrix per atom
        rot_mats.extend([rot_mat] * atoms_per_image)
    return rot_mats


def procrustes(geometry):
    # http://nghiaho.com/?page_id=671#comment-559906
    image0 = geometry.images[0]
    coords3d = image0.coords3d
    centroid = coords3d.mean(axis=0)
    last_centered = coords3d - centroid
    geometry.set_coords_at(0, last_centered.flatten())
    atoms_per_image = len(image0.atoms)

    # Don't rotate the first image, so just add identitiy matrices
    # for every atom.
    rot_mats = [np.eye(3)] * atoms_per_image
    for i, image in enumerate(geometry.images[1:], 1):
        coords3d = image.coords3d
        centroid = coords3d.mean(axis=0)
        # Center next image
        centered = coords3d - centroid
        tmp_mat = centered.T.dot(last_centered)
        U, W, Vt = np.linalg.svd(tmp_mat)
        rot_mat = U.dot(Vt)
        # Avoid reflections
        if np.linalg.det(rot_mat) < 0:
            U[:, -1] *= -1
            rot_mat = U.dot(Vt)
        # Rotate the coords
        rotated3d = centered.dot(rot_mat)
        geometry.set_coords_at(i, rotated3d.flatten())
        last_centered = rotated3d
        # Append one rotation matrix per atom
        rot_mats.extend([rot_mat] * atoms_per_image)
    return rot_mats


def align_coords(coords_list):
    coords_list = np.array(coords_list)
    coord_num = len(coords_list)
    aligned_coords = np.empty_like(coords_list).reshape(coord_num, -1, 3)

    coords0 = coords_list[0]
    coords0_3d = coords0.reshape(-1, 3)
    centroid = coords0_3d.mean(axis=0)
    prev_centered = coords0_3d - centroid
    aligned_coords[0] = prev_centered

    for i, coords in enumerate(coords_list[1:], 1):
        coords3d = coords.reshape(-1, 3)
        centroid = coords3d.mean(axis=0)
        # Center next image
        centered = coords3d - centroid
        tmp_mat = centered.T.dot(prev_centered)
        U, W, Vt = np.linalg.svd(tmp_mat)
        rot_mat = U.dot(Vt)
        # Avoid reflections
        if np.linalg.det(rot_mat) < 0:
            U[:, -1] *= -1
            rot_mat = U.dot(Vt)
        # Rotate the coords
        rotated3d = centered.dot(rot_mat)
        aligned_coords[i] = rotated3d
        prev_centered = rotated3d
    aligned_coords.reshape(coord_num, -1)
    return aligned_coords


def fit_rigid(geometry, vectors=(), vector_lists=(), hessian=None):
    rotated_vector_lists = list()
    rotated_hessian = None

    rot_mats = procrustes(geometry)
    G = sp.linalg.block_diag(*rot_mats)
    rotated_vectors = [vec.dot(G) for vec in vectors]
    for vl in vector_lists:
        rvl = [vec.dot(G) for vec in vl]
        rotated_vector_lists.append(rvl)

    if hessian is not None:
        # rotated_hessian = G.dot(hessian).dot(G.T)
        # rotated_hessian = G.T.dot(hessian).dot(G)
        rotated_hessian = G * hessian * G.T
    return rotated_vectors, rotated_vector_lists, rotated_hessian


def slugify_worker(dask_worker):
    slug = re.sub("tcp://", "host_", dask_worker)
    slug = re.sub(r"\.", "_", slug)
    slug = re.sub(":", "-", slug)
    return slug


def match_geoms(ref_geom, geom_to_match, hydrogen=False):
    """
    See
        [1] 10.1021/ci400534h
        [2] 10.1021/acs.jcim.6b00516
    """

    logging.warning(
        "helpers.match_geoms is deprecated!"
        "Use stocastic.align.match_geom_atoms instead!"
    )

    assert len(ref_geom.atoms) == len(geom_to_match.atoms), "Atom numbers don't match!"

    ref_coords, _ = ref_geom.coords_by_type
    coords_to_match, inds_to_match = geom_to_match.coords_by_type
    atoms = ref_coords.keys()
    for atom in atoms:
        # Only match hydrogens if explicitly requested
        if atom == "H" and not hydrogen:
            continue
        print("atom", atom)
        ref_coords_for_atom = ref_coords[atom]
        coords_to_match_for_atom = coords_to_match[atom]
        # Pairwise distances between two collections
        # Atoms of ref_geom are along the rows, atoms of geom_to_match
        # along the columns.
        cd = cdist(ref_coords_for_atom, coords_to_match_for_atom)
        print(cd)
        # Hungarian method, row_inds are returned already sorted.
        row_inds, col_inds = linear_sum_assignment(cd)
        print("col_inds", col_inds)
        old_inds = inds_to_match[atom]
        new_inds = old_inds[col_inds]
        print("old_inds", old_inds)
        print("new_inds", new_inds)
        new_coords_for_atom = coords_to_match_for_atom[new_inds]
        # print(ref_coords_for_atom)
        # print(new_coords_for_atom)
        # print(ref_coords_for_atom-new_coords_for_atom)
        # Update coordinates
        print("old coords")
        c3d = geom_to_match.coords3d
        print(c3d)
        # Modify the coordinates directly
        c3d[old_inds] = new_coords_for_atom
        # coords_to_match[atom] = coords_to_match_for_atom[new_inds]


def check_for_end_sign(check_user=True, cwd="."):
    signs = (
        "stop",
        "converged",
        "exit",
    )
    sign_found = False
    cur_user = getpass.getuser()
    cwd = Path(cwd)

    def sign_owner(path):
        if not check_user:
            return True
        else:
            return path.owner() == cur_user

    for sign in signs:
        sign_path = cwd / sign
        if sign_path.exists() and sign_owner(sign_path):
            print(f"Found sign '{sign}'. Ending run.")
            os.remove(sign_path)
            sign_found = sign

            if sign == "exit":
                sys.exit()
    return sign_found


def index_array_from_overlaps(overlaps, axis=1):
    """It is assumed that the overlaps between two points with indices
    i and j with (j > i) are computed and that i changes along the first
    axis (axis=0) and j changes along the second axis (axis=1).

    So the first row of the overlap matrix (overlaps[0]) should contain
    the overlaps between state 0 at index i and all states at index j.

    argmax along axis 1 returns the indices of the most overlapping states
    at index j with the states at index i, given by the item index in the
    indices array. E.g.:
        [0 1 3 2] indicates a root flip in a system with four states when
        going from index i to index j. Root 2 at i became root 3 at j and
        vice versa.
    """
    # indices = np.argmax(overlaps**2, axis=1)
    indices = np.argmax(np.abs(overlaps), axis=1)
    return indices


def np_print(func, precision=2, suppress=True, linewidth=120):
    def wrapped(*args, **kwargs):
        org_print_dict = dict(np.get_printoptions())
        np.set_printoptions(suppress=suppress, precision=precision, linewidth=linewidth)
        result = func(*args, **kwargs)
        np.set_printoptions(**org_print_dict)
        return result

    return wrapped


def confirm_input(message):
    full_message = message + " (yes/no)\n"
    inp = input(full_message)
    return inp == "yes"


def get_geom_getter(ref_geom, calc_setter):
    def geom_from_coords(coords):
        new_geom = ref_geom.copy()
        new_geom.coords = coords
        new_geom.set_calculator(calc_setter())
        return new_geom

    return geom_from_coords


def get_coords_diffs(coords, align=False, normalize=True):
    if align:
        coords = align_coords(coords)
    cds = [
        0,
    ]
    for i in range(len(coords) - 1):
        diff = np.linalg.norm(coords[i + 1] - coords[i])
        cds.append(diff)
    cds = np.cumsum(cds)
    if normalize:
        cds /= cds.max()
    return cds


def pick_image_inds(cart_coords, images: int):
    """Pick approx. evenly distributed images from given Cartesian coordinates."""
    cds = get_coords_diffs(cart_coords, align=True)
    target_cds = iter(np.linspace(0, 1, num=images))
    target_cd = next(target_cds)
    image_inds = list()
    for i, cd in enumerate(cds):
        if len(image_inds) == (images - 1):
            break
        if cd >= target_cd:
            image_inds.append(i)
            target_cd = next(target_cds)

    image_inds.append(len(cart_coords) - 1)
    return image_inds


def shake_coords(coords, scale=0.1, seed=None):
    if seed:
        np.random.seed(seed)
    offset = np.random.normal(scale=scale, size=coords.size)
    return coords + offset


def rms(arr):
    return np.sqrt(np.mean(arr**2))


def norm_max_rms(arr):
    norm = np.linalg.norm(arr)
    max_ = np.max(np.abs(arr))
    rms_ = rms(arr)
    return norm, max_, rms_


def complete_fragments(atoms, fragments):
    lengths = [len(frag) for frag in fragments]

    frag_atoms = sum(lengths)

    all_inds = set(range(len(atoms)))
    frag_inds = set(it.chain(*fragments))
    rest_inds = all_inds - frag_inds

    assert len(frag_inds) + len(rest_inds) == len(atoms)
    assert frag_inds & rest_inds == set()

    if rest_inds:
        fragments.append(tuple(rest_inds))
    fragments = tuple([tuple(frag) for frag in fragments])

    return fragments


FinalHessianResult = namedtuple(
    "FinalHessianResult",
    "neg_eigvals eigvals nus imag_fns thermo",
)


def do_final_hessian(
    geom,
    save_hessian=True,
    write_imag_modes=False,
    is_ts=False,
    prefix="",
    T=T_DEFAULT,
    p=p_DEFAULT,
    ev_thresh=-1e-6,
    print_thermo=False,
    out_dir=None,
):
    if out_dir is None:
        out_dir = "."
    out_dir = Path(out_dir)

    print(highlight_text("Hessian at final geometry", level=1))
    print()

    report_isotopes(geom, "the_frequencies")

    start = time.time()
    print("... started Hessian calculation")
    hessian = geom.cart_hessian
    dur = time.time() - start
    print(f"... calculation took {dur/60:.2f} min")
    print("... mass-weighing cartesian hessian")
    mw_hessian = geom.mass_weigh_hessian(hessian)
    print("... doing Eckart-projection")
    proj_hessian = geom.eckart_projection(mw_hessian)
    eigvals, _ = np.linalg.eigh(proj_hessian)

    neg_inds = eigvals < ev_thresh
    neg_eigvals = eigvals[neg_inds]
    neg_num = sum(neg_inds)
    eigval_str = np.array2string(eigvals[:10], precision=4)
    print()
    print("First 10 eigenvalues", eigval_str)
    if neg_num > 0:
        wavenumbers = eigval_to_wavenumber(neg_eigvals)
        wavenum_str = np.array2string(wavenumbers, precision=2)
        print("Imaginary frequencies:", wavenum_str, "cm⁻¹")

    if prefix:
        prefix = f"{prefix}_"

    # Dump HDF Hessian
    if save_hessian:
        final_h5_hessian_fn = prefix + "final_hessian.h5"
        save_h5_hessian(out_dir / final_h5_hessian_fn, geom)
        print(f"Wrote Hessian data HD5 file '{final_h5_hessian_fn}'.")

    imag_fns = list()
    if write_imag_modes:
        imag_modes = imag_modes_from_geom(geom)
        for i, imag_mode in enumerate(imag_modes):
            trj_fn = out_dir / (prefix + f"imaginary_mode_{i:03d}.trj")
            imag_fns.append(trj_fn)
            with open(trj_fn, "w") as handle:
                handle.write(imag_mode.trj_str)
            print(
                f"Wrote imaginary mode with ṽ={imag_mode.nu: >10.2f} cm⁻¹ to '{trj_fn}'"
            )
        print()

    thermo = None
    if can_thermoanalysis:
        thermo = geom.get_thermoanalysis(geom, T=T, p=p)
        if print_thermo:
            print_thermoanalysis(thermo, geom=geom, is_ts=is_ts)

    res = FinalHessianResult(
        neg_eigvals=neg_eigvals,
        eigvals=eigvals,
        nus=eigval_to_wavenumber(eigvals),
        imag_fns=imag_fns,
        thermo=thermo,
    )
    return res


def print_barrier(ref_energy, comp_energy, ref_str, comp_str):
    barrier = (ref_energy - comp_energy) * AU2KJPERMOL
    print(f"Barrier between {ref_str} and {comp_str}: {barrier:.1f} kJ mol⁻¹")
    return barrier


def get_tangent_trj_str(atoms, coords, tangent, comment=None, points=10, displ=None):
    if displ is None:
        # Linear equation. Will give displ~3 for 30 atoms and
        # displ ~ 1 for 3 atoms.
        # displ = 2/27 * len(atoms) + 0.78

        # Logarithmic function f(x) = a*log(x) + b
        # f(3) = ~1 and (f30) = ~2 with a = 0.43429 and b = 0.52288
        # I guess this works better, because only some atoms move, even in bigger
        # systems and the log function converges against a certain value, whereas
        # the linear function just keeps growing.
        displ = 0.43429 * log(len(atoms)) + 0.52288
    step_sizes = np.linspace(-displ, displ, 2 * points + 1)
    steps = step_sizes[:, None] * tangent
    trj_coords = coords[None, :] + steps
    trj_coords = trj_coords.reshape(step_sizes.size, -1, 3) / ANG2BOHR

    comments = None
    if comment:
        comments = [comment] * step_sizes.size
    trj_str = make_trj_str(atoms, trj_coords, comments=comments)

    return trj_str


def imag_modes_from_geom(geom, freq_thresh=-10, points=10, displ=None):
    NormalMode = namedtuple("NormalMode", "nu mode trj_str")

    # We don't want to do start any calculation here, so we directly access
    # the attribute underlying the geom.hessian property.
    nus, _, _, cart_displs = geom.get_normal_modes(geom._hessian)
    below_thresh = nus < freq_thresh

    imag_modes = list()
    for nu, eigvec in zip(nus[below_thresh], cart_displs[:, below_thresh].T):
        comment = f"{nu:.2f} cm⁻¹"
        trj_str = get_tangent_trj_str(
            geom.atoms,
            geom.cart_coords,
            eigvec,
            comment=comment,
            points=points,
            displ=displ,
        )
        imag_modes.append(
            NormalMode(
                nu=nu,
                mode=eigvec,
                trj_str=trj_str,
            )
        )

    return imag_modes
