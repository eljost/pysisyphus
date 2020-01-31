#!/usr/bin/env python3

from collections import namedtuple
import itertools as it
import logging
import os
from pathlib import Path
import re

import numpy as np
import scipy as sp
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from pysisyphus.constants import ANG2BOHR, AU2J, AMU2KG, BOHR2M
from pysisyphus.Geometry import Geometry
from pysisyphus.xyzloader import parse_xyz_file, parse_trj_file


THIS_DIR = Path(os.path.abspath(os.path.dirname(__file__)))


def geom_from_xyz_file(xyz_fn, **kwargs):
    atoms, coords, comment = parse_xyz_file(xyz_fn, with_comment=True)
    coords *= ANG2BOHR
    geom = Geometry(atoms, coords.flatten(), comment=comment, **kwargs)
    return geom


def geoms_from_trj(trj_fn, first=None, **kwargs):
    atoms_coords_comments = parse_trj_file(trj_fn, with_comments=True)[:first]
    geoms = [Geometry(atoms, coords.flatten()*ANG2BOHR, comment=comment, **kwargs)
             for atoms, coords, comment in atoms_coords_comments
    ]
    return geoms


def geom_from_library(xyz_fn, **kwargs):
    xyz_dir = THIS_DIR / "../xyz_files/"
    xyz_fn = xyz_dir / xyz_fn
    if xyz_fn.suffix == ".xyz":
        return geom_from_xyz_file(xyz_fn, **kwargs)
    elif xyz_fn.suffix == ".trj":
        return geoms_from_trj(xyz_fn, **kwargs)
    else:
        raise Exception("Unknown filetype!")


def get_baker_geoms(**kwargs):
    baker_path = THIS_DIR / "../xyz_files/baker"
    xyz_fns = baker_path.glob("*.xyz")
    geoms = {
        xyz_fn.name: geom_from_xyz_file(xyz_fn, **kwargs) for xyz_fn in xyz_fns
    }
    # del geoms["acetylene.xyz"]
    # From 10.1002/jcc.540140910
    sto3g_energies = {
        "water.xyz": -74.96590,
        "ammonia.xyz": -55.45542,
        "ethane.xyz": -78.30618,
        "acetylene.xyz": -75.85625,
        "allene.xyz": -114.42172,
        "hydroxysulphane.xyz": -468.12592,
        "benzene.xyz": -227.89136,
        "methylamine.xyz": -94.01617,
        "ethanol.xyz": -152.13267,
        "acetone.xyz": -189.53603,
        "disilylether.xyz": -648.58003,
        "135trisilacyclohexane.xyz": -976.13242,
        "benzaldehyde.xyz": -339.12084,
        "13difluorobenzene.xyz": -422.81106,
        "135trifluorobenzene.xyz": -520.27052,
        "neopentane.xyz": -194.04677,
        "furan.xyz": -225.75126,
        "naphthalene.xyz": -378.68685,
        "15difluoronaphthalene.xyz": -573.60633,
        "2hydroxybicyclopentane.xyz": -265.46482,
        "achtar10.xyz": -356.28265,
        "acanil01.xyz": -432.03012,
        "benzidine.xyz": -563.27798,
        "pterin.xyz": -569.84884,
        "difuropyrazine.xyz": -556.71910,
        "mesityloxide.xyz": -304.05919,
        "histidine.xyz": -538.54910,
        "dimethylpentane.xyz": -271.20088,
        "caffeine.xyz": -667.73565,
        "menthone.xyz": -458.44639,
    }
    # return geoms, sto3g_energies
    # Join both dicts
    baker_dict = {
        name: (geom, sto3g_energies[name]) for name, geom in geoms.items()
    }
    return baker_dict


def get_baker_ts_geoms(**kwargs):
    baker_ts_path = THIS_DIR / "../xyz_files/baker_ts"
    xyz_fns = baker_ts_path.glob("*.xyz")
    geoms = {
        xyz_fn.name: geom_from_xyz_file(xyz_fn, **kwargs) for xyz_fn in xyz_fns
        if not ("downhill" in xyz_fn.name)
    }
    # From 10.1002/(SICI)1096-987X(199605)17:7<888::AID-JCC12>3.0.CO;2-7
    meta_data = {
        "01_hcn.xyz": (0, 1, -92.24604),
        "02_hcch.xyz": (0, 1, -76.29343),
        "03_h2co.xyz": (0, 1, -113.05003),
        "04_ch3o.xyz": (0, 2, -113.69365),
        "05_cyclopropyl.xyz": (0, 2, -115.72100),
        "06_bicyclobutane.xyz": (0, 1, -153.90494),
        "07_bicyclobutane.xyz": (0, 1, -153.89754),
        "08_formyloxyethyl.xyz": (0, 2, -264.64757),
        "09_parentdieslalder.xyz": (0, 1, -231.60321),
        # 10 and 11 don't have any imaginary frequencies at the given
        # geometry, so they may be skipped.
        "10_tetrazine.xyz": (0, 1, -292.81026),
        "11_trans_butadiene.xyz": (0, 1, -154.05046),
        "12_ethane_h2_abstraction.xyz": (0, 1, -78.54323),
        "13_hf_abstraction.xyz": (0, 1, -176.98453),
        "14_vinyl_alcohol.xyz": (0, 1, -151.91310),
        # 15 does not have an imaginary mode in cartesian coordinates
        "15_hocl.xyz": (0, 1, -596.87865),
        "16_h2po4_anion.xyz": (-1, 1, -637.92388),
        "17_claisen.xyz": (0, 1, -267.23859),
        "18_silyene_insertion.xyz": (0, 1, -367.20778),
        "19_hnccs.xyz": (0, 1, -525.43040),
        # The energy given in the paper (-168.24752 au) is the correct one
        # if one forms the central (0,1) bond (0-based indexing). If this
        # bond is missing, as it is if we autogenerate with bond-factor=1.3
        # then a TS with -168.241392 will be found.
        # For now we will use the original value from the paper.
        "20_hconh3_cation.xyz": (1, 1, -168.24752),
        "21_acrolein_rot.xyz": (0, 1, -189.67574),
        "22_hconhoh.xyz": (0, 1, -242.25529),
        "23_hcn_h2.xyz": (0, 1, -93.31114),
        "24_h2cnh.xyz": (0, 1, -93.33296),
        "25_hcnh2.xyz": (0, 1, -93.28172),
    }
    return geoms, meta_data


def get_baker_ts_geoms_flat(**kwargs):
     geoms, meta_data = get_baker_ts_geoms(**kwargs)
     # Key: (geometry, charge, mult, ref_energy)
     return [(mol,) + (geom, ) + meta_data[mol] for mol, geom in geoms.items()]


def get_baker_opt_ts_geoms(**kwargs):
    meta_data = {
        "01_hcn_opt_ts.xyz": (0, 1),
        # "02_hcch_opt_ts.xyz": (0, 1),
        "03_h2co_opt_ts.xyz": (0, 1),
        "04_ch3o_opt_ts.xyz": (0, 2),
        "05_cyclopropyl_opt_ts.xyz": (0, 2),
        "06_bicyclobutane_opt_ts.xyz": (0, 1),
        "07_bicyclobutane_opt_ts.xyz": (0, 1),
        "08_formyloxyethyl_opt_ts.xyz": (0, 2),
        # "09_parentdieslalder_opt_ts.xyz": (0, 1),
        # "10_tetrazine_opt_ts.xyz": (0, 1),
        # "11_trans_butadiene_opt_ts.xyz": (0, 1),
        # "12_ethane_h2_abstraction_opt_ts.xyz": (0, 1),
        "13_hf_abstraction_opt_ts.xyz": (0, 1),
        "14_vinyl_alcohol_opt_ts.xyz": (0, 1),
        # "15_hocl_opt_ts.xyz": (0, 1),
        "16_h2po4_anion_opt_ts.xyz": (-1, 1),
        # "17_claisen_opt_ts.xyz": (0, 1),
        "18_silyene_insertion_opt_ts.xyz": (0, 1),
        "19_hnccs_opt_ts.xyz": (0, 1),
        "20_hconh3_cation_opt_ts.xyz": (1, 1),
        "21_acrolein_rot_opt_ts.xyz": (0, 1),
        # "22_hconhoh_opt_ts.xyz": (0, 1),
        "23_hcn_h2_opt_ts.xyz": (0, 1),
        "24_h2cnh_opt_ts.xyz": (0, 1),
        "25_hcnh2_opt_ts.xyz": (0, 1),
    }
    baker_ts_path = THIS_DIR / "../xyz_files/baker_opt_ts"
    geoms = {
        xyz_fn: geom_from_xyz_file(baker_ts_path / xyz_fn, **kwargs)
        for xyz_fn in meta_data.keys()
    }

    return geoms, meta_data


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
    rot_mats = [np.eye(3)]*atoms_per_image
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
        rot_mats.extend([rot_mat]*atoms_per_image)
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
    rot_mats = [np.eye(3)]*atoms_per_image
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
        rot_mats.extend([rot_mat]*atoms_per_image)
    return rot_mats


def fit_rigid(geometry, vectors=(), vector_lists=(), hessian=None):
    rotated_vector_lists = list()
    rotated_hessian = None

    rot_mats = procrustes(geometry)
    G = sp.linalg.block_diag(*rot_mats)
    rotated_vectors = [vec.dot(G) for vec in vectors]
    for vl in vector_lists:
        rvl = [vec.dot(G) for vec in vl]
        rotated_vector_lists.append(rvl)

    if hessian:
        #rotated_hessian = G.dot(hessian).dot(G.T)
        rotated_hessian = G.T.dot(hessian).dot(G)
        #rotated_hessian = G*hessian*G.T
    return rotated_vectors, rotated_vector_lists, rotated_hessian


def chunks(l, n):
    """Yield successive n-sized chunks from l.
    https://stackoverflow.com/a/312464
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]


def slugify_worker(dask_worker):
    slug = re.sub("tcp://", "host_", dask_worker)
    slug = re.sub("\.", "_", slug)
    slug = re.sub(":", "-", slug)
    return slug


def match_geoms(ref_geom, geom_to_match, hydrogen=False):
    """
    See
        [1] 10.1021/ci400534h
        [2] 10.1021/acs.jcim.6b00516
    """

    logging.warning("helpers.match_geoms is deprecated!"
                    "Use stocastic.align.match_geom_atoms instead!")

    assert len(ref_geom.atoms) == len(geom_to_match.atoms), \
        "Atom numbers don't match!"

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


def check_for_stop_sign():
    stop_signs = ("stop", "STOP")
    stop_sign_found = False

    for ss in stop_signs:
        if os.path.exists(ss):
            print("Found stop sign. Stopping run.")
            os.remove(ss)
            stop_sign_found = True
    return stop_sign_found


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
        np.set_printoptions(suppress=suppress,
                            precision=precision,
                            linewidth=linewidth)
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


def get_coords_diffs(coords):
    cds = [0, ]
    for i in range(len(coords)-1):
        diff = np.linalg.norm(coords[i+1]-coords[i])
        cds.append(diff)
    cds = np.cumsum(cds)
    cds /= cds.max()
    return cds


def shake_coords(coords, scale=0.1, seed=None):
    if seed:
        np.random.seed(seed)
    offset = np.random.normal(scale=scale, size=coords.size)
    return coords + offset


def highlight_text(text, width=80):
    full_length = len(text) + 4
    pad_len = width - full_length
    pad_len = (pad_len - (pad_len % 2)) // 2
    pad = " " * pad_len
    full_row = "#" * full_length
    highlight = f"""{pad}{full_row}\n{pad}# {text.upper()} #\n{pad}{full_row}"""
    return highlight


def rms(arr):
    return np.sqrt(np.mean(arr**2))


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


def eigval_to_wavenumber(ev):
    conv = AU2J/(AMU2KG*BOHR2M**2)

    return np.sign(ev) * np.sqrt(np.abs(ev)*conv)/(2*np.pi*3e10)


FinalHessianResult = namedtuple("FinalHessianResult",
                                "neg_eigvals"
)


def do_final_hessian(geom, save_hessian=True):
    print("Calculating hessian at final geometry.")

    # TODO: Add cartesian_hessian property to Geometry to avoid
    # accessing a "private" attribute.
    hessian = geom.cart_hessian
    print("... mass-weighing cartesian hessian")
    mw_hessian = geom.mass_weigh_hessian(hessian)
    print("... doing eckart-projection")
    proj_hessian = geom.eckart_projection(mw_hessian)
    eigvals, eigvecs = np.linalg.eigh(proj_hessian)
    ev_thresh = -1e-6

    neg_inds = eigvals < ev_thresh
    neg_eigvals = eigvals[neg_inds]
    neg_num = sum(neg_inds)
    eigval_str = np.array2string(eigvals[:10], precision=4)
    print()
    print("First 10 eigenvalues", eigval_str)
    # print(f"Self found {neg_num} eigenvalue(s) < {ev_thresh}.")
    if neg_num > 0:
        wavenumbers = eigval_to_wavenumber(neg_eigvals)
        wavenum_str = np.array2string(wavenumbers, precision=2)
        print("Imaginary frequencies:", wavenum_str, "cm⁻¹")

    if save_hessian:
        final_hessian_fn = "calculated_final_cart_hessian"
        np.savetxt(final_hessian_fn, hessian)
        print()
        print(f"Wrote final (not mass-weighted) hessian to '{final_hessian_fn}'.")

    res = FinalHessianResult(neg_eigvals=neg_eigvals)
    return res
