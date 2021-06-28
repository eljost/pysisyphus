# [1] https://doi.org/10.1002/jcc.26495
#     Habershon, 2021

"""
prp a901cdfacc579eb63b193cbc9043212e8b57746f
pysis 340ab6105ac4156f0613b4d0e8f080d9f195530c
do_trans accidentally disabled in transtorque
"""

from functools import reduce
import itertools as it

import numpy as np

from pysisyphus.calculators import (
    HardSphere,
    TransTorque,
    AtomAtomTransTorque,
    Composite,
)
from pysisyphus.constants import BOHR2ANG
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import align_coords
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.init_logging import init_logging
from pysisyphus.intcoords.setup import get_fragments, get_bond_sets
from pysisyphus.xyzloader import coords_to_trj, make_xyz_str

init_logging()


class SteepestDescent:
    def __init__(
        self,
        geom,
        max_cycles=1000,
        max_step=0.05,
        rms_force=0.05,
        rms_force_only=True,
        prefix=None,
        dump=False,
        print_mod=25,
    ):
        self.geom = geom
        self.max_cycles = max_cycles
        self.max_step = max_step
        self.rms_force = rms_force
        self.rms_force_only = rms_force_only
        self.prefix = prefix
        self.dump = dump
        self.print_mod = print_mod

        self.all_coords = np.zeros((max_cycles, self.geom.coords.size))

    def run(self):
        coords = self.geom.coords.copy()

        to_dump = []

        for i in range(self.max_cycles):
            self.all_coords[i] = coords.copy()
            if self.dump and (i % 100) == 0:
                to_dump.append(self.geom.as_xyz(cart_coords=coords))
            results = self.geom.get_energy_and_forces_at(coords)
            forces = results["forces"]
            norm = np.linalg.norm(forces)
            rms = np.sqrt(np.mean(forces ** 2))
            if rms <= self.rms_force:
                print(f"Converged in cycle {i}. Breaking.")
                break

            if i > 0:
                beta = forces.dot(forces) / self.prev_forces.dot(self.prev_forces)
                step = forces + beta * self.prev_step
            else:
                step = forces.copy()
            # step = forces.copy()

            step *= min(self.max_step / np.abs(step).max(), 1)
            if i % self.print_mod == 0:
                print(
                    f"{i:03d}: |forces|={norm: >12.6f} "
                    f"rms(forces)={np.sqrt(np.mean(forces**2)): >12.6f} "
                    f"|step|={np.linalg.norm(step): >12.6f}"
                )
            coords += step

            self.prev_step = step
            self.prev_forces = forces
        self.geom.coords = coords
        self.all_coords = self.all_coords[: i + 1]

        if to_dump:
            with open("optimization.trj", "w") as handle:
                handle.write("\n".join(to_dump))


def get_fragments_and_bonds(geoms):
    if isinstance(geoms, Geometry) or len(geoms) == 1:
        geom = geoms
        atoms = geom.atoms
        coords3d = geom.coords3d
        bonds = [frozenset(bond) for bond in get_bond_sets(atoms, coords3d)]
        fragments = get_fragments(atoms, coords3d.flatten(), bond_inds=bonds)
        frag_inds = list(it.chain(*fragments))
        if len(frag_inds) != len(atoms):
            all_inds = list(range(len(atoms)))
            missing_inds = set(all_inds) - set(frag_inds)
            for mi in missing_inds:
                fragments.append(frozenset((mi,)))

        frag_bonds = [
            list(filter(lambda bond: bond <= frag, bonds)) for frag in fragments
        ]
        # frag_atoms = [[a for i, a in enumerate(atoms) if i in frag] for frag in fragments]

        # Assert that we do not have any interfragment bonds
        assert reduce((lambda x, y: x + len(y)), frag_bonds, 0) == len(bonds)
        union_geom = geom.copy(coord_type="cart")
    else:
        # Form union, determine consistent new indices for all atoms and calculate bonds
        raise Exception()

    # return fragments, frag_bonds, set(bonds), frag_atoms
    return fragments, frag_bonds, set(bonds), union_geom


def get_rot_mat(coords3d_1, coords3d_2, center=False):
    coords3d_1 = coords3d_1.copy().reshape(-1, 3)
    coords3d_2 = coords3d_2.copy().reshape(-1, 3)

    def _center(coords3d):
        return coords3d - coords3d.mean(axis=0)

    if center:
        coords3d_1 = _center(coords3d_1)
        coords3d_2 = _center(coords3d_2)

    tmp_mat = coords3d_1.T.dot(coords3d_2)
    U, W, Vt = np.linalg.svd(tmp_mat)
    rot_mat = U.dot(Vt)
    # Avoid reflections
    if np.linalg.det(rot_mat) < 0:
        U[:, -1] *= -1
        rot_mat = U.dot(Vt)
    return rot_mat


def get_steps_to_active_atom_mean(
    frag_lists, iter_frag_lists, ind_dict, coords3d, skip=True
):
    frag_num = len(frag_lists)
    steps = np.zeros((frag_num, 3))
    for m, frag_m in enumerate(frag_lists):
        step_m = np.zeros(3)
        for n, _ in enumerate(iter_frag_lists):
            if skip and m == n:
                continue
            active_inds = ind_dict[(n, m)]
            if len(active_inds) == 0:
                continue
            step_m += coords3d[active_inds].mean(axis=0)
        step_m /= frag_num
        steps[m] = step_m
    return steps


def report_frags(rgeom, pgeom, rfrags, pfrags, rbond_diff, pbond_diff):
    for name, geom in (("Reactant(s)", rgeom), ("Product(s)", pgeom)):
        print(f"{name}: {geom}\n\n{geom.as_xyz()}\n")

    def get_frag_atoms(geom, frag):
        atoms = geom.atoms
        return [atoms[i] for i in frag]

    for name, geom, frags in (("reactant", rgeom, rfrags), ("product", pgeom, pfrags)):
        print(f"{len(frags)} fragment(s) in {name} image\n")
        for frag in frags:
            frag_atoms = get_frag_atoms(geom, frag)
            frag_coords = geom.coords3d[list(frag)]
            frag_xyz = make_xyz_str(frag_atoms, frag_coords * BOHR2ANG)
            print(frag_xyz + "\n")

    def print_bonds(geom, bonds):
        for from_, to_ in bonds:
            from_atom, to_atom = [geom.atoms[i] for i in (from_, to_)]
            print(f"\t({from_: >3d}{from_atom} - {to_: >3d}{to_atom})")

    print("Bonds broken in reactant image:")
    print_bonds(rgeom, rbond_diff)
    print()
    print("Bonds formed in product image:")
    print_bonds(pgeom, pbond_diff)
    print()


def report_mats(name, mats):
    for (m, n), indices in mats.items():
        print(f"{name}({m}, {n}): {indices}")
    print()


def center_fragments(frag_list, geom):
    c3d = geom.coords3d
    for frag in frag_list:
        mean = c3d[frag].mean(axis=0)
        c3d[frag] -= mean[None, :]


def get_which_frag(frags):
    which_frag = dict()
    for frag_ind, frag in enumerate(frags):
        which_frag.update({atom_ind: frag_ind for atom_ind in frag})
    return which_frag


def form_A(frags, which_frag, formed_bonds):
    """Construct the A-matrices.

    AR[(m, n)] (AP[(m, n)]) contains the subset of atoms in Rm (Pm) that forms
    bonds with Rn (Pn).
    """
    A = dict()
    for m, n in formed_bonds:
        key = (which_frag[m], which_frag[n])
        A.setdefault(key, list()).append(m)
        A.setdefault(key[::-1], list()).append(n)
    return A


CONFIG = {
    "s2_hs_kappa": 1.0,
    "s4_hs_kappa": 50.0,
    "s4_v_kappa": 1.0,
    "s4_w_kappa": 1.0,
    "s5_v_kappa": 1.0,
    "s5_w_kappa": 3.0,
    "s5_hs_kappa": 10.0,
    "s5_z_kappa": 2.0,
    "s5_trans": True,
    "s5_rms_force": 0.01,
}


def precon_pos_rot(reactants, products, prefix=None, config=CONFIG):
    c = config

    if prefix is None:
        prefix = ""

    def make_fn(fn):
        return prefix + fn

    rfrags, rfrag_bonds, rbonds, runion = get_fragments_and_bonds(reactants)
    pfrags, pfrag_bonds, pbonds, punion = get_fragments_and_bonds(products)

    pbond_diff = pbonds - rbonds  # Present in product(s)
    rbond_diff = rbonds - pbonds  # Present in reactant(s)
    involved_atoms = set(tuple(it.chain(*pbond_diff)))
    involved_atoms |= set(tuple(it.chain(*rbond_diff)))

    which_rfrag = get_which_frag(rfrags)
    which_pfrag = get_which_frag(pfrags)

    rfrag_lists = [list(frag) for frag in rfrags]
    pfrag_lists = [list(frag) for frag in pfrags]

    report_frags(runion, punion, rfrags, pfrags, rbond_diff, pbond_diff)

    def form_C(m_frags, n_frags):
        """Construct the C-matrices.

        Returns a dict with (m, n) keys, containing the respective
        unions of rectant fragment n and product fragment m.
        """
        C = dict()
        for m, m_frag in enumerate(m_frags):
            for n, n_frag in enumerate(n_frags):
                C[(m, n)] = list(m_frag & n_frag)
        return C

    CR = form_C(rfrags, pfrags)
    assert len(set(it.chain(*CR.values()))) == len(runion.atoms)
    CP = {(n, m): union for (m, n), union in CR.items()}
    print("CR(m, n), subset of atoms in molecule Rn which are in Pm after reaction.")
    report_mats("CR", CR)
    print("CP(m, n), subset of atoms in molecule Pn which are in Rm before reaction.")
    report_mats("CP", CP)

    def form_B(C):
        """Construct the B-matrices.

        Returns a dict with (m, n) keys, containing the respective
        subsets of C[(m, n)] that acutally participate in bond-breaking/forming.
        """
        B = dict()
        for (m, n), union in C.items():
            key = (m, n)
            B.setdefault(key, set())
            B[key] |= set(union) & involved_atoms
        for k, v in B.items():
            B[k] = list(v)
        return B

    BR = form_B(CR)
    BP = form_B(CP)
    print(
        "BR(m, n), subset of atoms in CRnm actually involved in bond forming/breaking."
    )
    report_mats("BR", BR)
    print(
        "BP(m, n), subset of atoms in CPnm actually involved in bond forming/breaking."
    )
    report_mats("BP", BP)

    AR = form_A(rfrags, which_rfrag, pbond_diff)
    AP = form_A(pfrags, which_pfrag, rbond_diff)
    print("AR(m, n), subset of atoms in Rm that form bonds to atoms in Rn.")
    report_mats("AR", AR)
    print(
        "AP(m, n), subset of atoms in Pm which had bonds with Pn (formerly bonded in R)."
    )
    report_mats("AP", AP)

    def form_G(A):
        G = dict()
        for (m, n), inds in A.items():
            G.setdefault(m, set())
            G[m] |= set(inds)

        for k, v in G.items():
            G[k] = list(v)
            assert len(v) > 0
        return G

    GR = form_G(AR)
    # GP = form_G(AP)
    print(f"GR: {GR}")
    # print(f"GP: {GP}")

    # Initial, centered, coordinates and 5 stages
    r_coords = np.zeros((6, runion.coords.size))
    p_coords = np.zeros((6, punion.coords.size))

    def backup_coords(stage):
        assert 0 <= stage < 6
        r_coords[stage] = runion.coords.copy()
        p_coords[stage] = punion.coords.copy()

    """
    STAGE 1
    Initial positioning of reactant and product molecules
    """

    # Center fragments at their geometric average
    center_fragments(rfrag_lists, runion)
    center_fragments(pfrag_lists, punion)

    backup_coords(0)

    # Translate reactant molecules
    alphas = get_steps_to_active_atom_mean(
        rfrag_lists, rfrag_lists, AR, runion.coords3d
    )
    for rfrag, alpha in zip(rfrag_lists, alphas):
        runion.coords3d[rfrag] += alpha

    # Translate product molecules
    betas = get_steps_to_active_atom_mean(
        pfrag_lists, rfrag_lists, BR, punion.coords3d, skip=False
    )
    sigmas = get_steps_to_active_atom_mean(
        pfrag_lists, rfrag_lists, CR, punion.coords3d, skip=False
    )
    bs_half = (betas + sigmas) / 2
    for pfrag, bsh in zip(pfrag_lists, bs_half):
        punion.coords3d[pfrag] += bsh

    backup_coords(1)
    print()

    """
    STAGE 2
    Intra-image Inter-molecular Hard-Sphere forces
    """

    print(highlight_text("Stage 2, Hard-Sphere Forces"))

    s2_hs_kappa = c["s2_hs_kappa"]

    def hardsphere_sd_opt(geom, frag_lists, title):
        print(highlight_text(title, level=1))
        calc = HardSphere(geom, frag_lists, kappa=s2_hs_kappa)
        geom.set_calculator(calc)
        opt_kwargs = {
            "max_cycles": 1000,
            "max_step": 0.5,
            "rms_force": 0.05,
        }
        opt = SteepestDescent(geom, **opt_kwargs)
        opt.run()

    hardsphere_sd_opt(runion, rfrag_lists, "Reactants")
    hardsphere_sd_opt(punion, pfrag_lists, "Products")

    backup_coords(2)
    print()

    """
    STAGE 3
    Initial orientation of molecules
    """

    print(highlight_text("Stage 3, Initial Orientation"))

    # Rotate R fragments
    if len(rfrag_lists) > 1:
        alphas = get_steps_to_active_atom_mean(rfrag_lists, rfrag_lists, AR, runion.coords3d)
        gammas = np.zeros_like(alphas)
        for m, rfrag in enumerate(rfrag_lists):
            Gm = GR[m]
            gammas[m] = runion.coords3d[Gm].mean(axis=0)
        r_means = np.array([runion.coords3d[frag].mean(axis=0) for frag in rfrag_lists])

        for m, rfrag in enumerate(rfrag_lists):
            gm = r_means[m]
            rot_mat = get_rot_mat(gammas[m] - gm, alphas[m] - gm)
            rot_coords = (runion.coords3d[rfrag] - gm).dot(rot_mat)
            runion.coords3d[rfrag] = rot_coords + gm - rot_coords.mean(axis=0)

    Ns = [0] * len(pfrag_lists)
    for (m, n), CPmn in CP.items():
        Ns[m] += len(CPmn)

    # Rotate P fragments
    for m, pfrag in enumerate(pfrag_lists):
        pc3d = punion.coords3d[pfrag]
        gm = pc3d.mean(axis=0)
        r0Pm = pc3d - gm[None, :]
        mu_Pm = np.zeros_like(r0Pm)
        N = Ns[m]
        for n, rfrag in enumerate(rfrag_lists):
            CPmn = CP[(m, n)]
            RPmRn = get_rot_mat(
                punion.coords3d[CPmn], runion.coords3d[CPmn], center=True
            )
            print(f"m={m}, n={n}, len(CPmn)={len(CPmn)}")
            # Eq. (A2) in [1]
            r0Pmn = np.einsum("ij,jk->ki", RPmRn, r0Pm.T)
            mu_Pm += len(CPmn) ** 2 / N * r0Pmn
        rot_mat = get_rot_mat(r0Pm, mu_Pm, center=True)
        rot_coords = r0Pm.dot(rot_mat)
        punion.coords3d[pfrag] = rot_coords + gm - rot_coords.mean(axis=0)

    backup_coords(3)
    print()

    """
    STAGE 4
    Alignment of reactive atoms

    This stage involves three forces: hard-sphere forces and two kinds
    of average translational (^t) and rotational (^r) forces (v and w,
    (A3) - (A5) in [1]).

    v^t and v^r arise from atoms in A^Rnm and A^Rmn, that is atoms that
    participate in bond forming/breaking in R. The translational force
    is usually attractive, which is counteracted by the repulsive hard-sphere
    forces.
    """

    print(highlight_text("Stage 4, Alignment Of Reactive Atoms"))

    def composite_sd_opt(geom, keys_calcs, title, rms_force=0.05):
        print(highlight_text(title, level=1))
        final = " + ".join([k for k in keys_calcs.keys()])
        calc = Composite(final, keys_calcs=keys_calcs)
        geom.set_calculator(calc)
        opt_kwargs = {
            "max_step": 0.05,
            "max_cycles": 2000,
            "rms_force": rms_force,
        }
        opt = SteepestDescent(geom, **opt_kwargs)
        opt.run()

    def get_vr_trans_torque(kappa=1.0, do_trans=True):
        return TransTorque(
            rfrag_lists, rfrag_lists, AR, AR, kappa=kappa, do_trans=do_trans
        )

    def r_weight_func(m, n, a, b):
        """As required for (A5) in [1]."""
        return 1 if a in BR[(m, n)] else 0.5

    def get_wr_trans_torque(kappa=1.0, do_trans=True):
        return TransTorque(
            rfrag_lists,
            pfrag_lists,
            CR,
            CP,
            weight_func=r_weight_func,
            skip=False,
            b_coords3d=punion.coords3d,
            kappa=kappa,
            do_trans=do_trans,
        )

    def get_vp_trans_torque(kappa=1.0, do_trans=True):
        return TransTorque(
            pfrag_lists, pfrag_lists, AP, AP, kappa=kappa, do_trans=do_trans
        )

    def p_weight_func(m, n, a, b):
        """As required for (A5) in [1]."""
        return 1 if a in BP[(m, n)] else 0.5

    def get_wp_trans_torque(kappa=1.0, do_trans=True):
        return TransTorque(
            pfrag_lists,
            rfrag_lists,
            CP,
            CR,
            weight_func=p_weight_func,
            skip=False,
            b_coords3d=runion.coords3d,
            kappa=kappa,
            do_trans=do_trans,
        )

    s4_hs_kappa = c["s4_hs_kappa"]
    s4_v_kappa = c["s4_v_kappa"]
    s4_w_kappa = c["s4_w_kappa"]

    vr_trans_torque = get_vr_trans_torque(kappa=s4_v_kappa)
    wr_trans_torque = get_wr_trans_torque(kappa=s4_w_kappa)
    r_keys_calcs = {
        "hardsphere": HardSphere(runion, rfrag_lists, kappa=s4_hs_kappa),
        "v": vr_trans_torque,
        "w": wr_trans_torque,
    }
    composite_sd_opt(runion, r_keys_calcs, "Reactants")

    vp_trans_torque = get_vp_trans_torque(kappa=s4_v_kappa)
    wp_trans_torque = get_wp_trans_torque(kappa=s4_w_kappa)
    p_keys_calcs = {
        "hardsphere": HardSphere(punion, pfrag_lists, kappa=s4_hs_kappa),
        "v": vp_trans_torque,
        "w": wp_trans_torque,
    }
    composite_sd_opt(punion, p_keys_calcs, "Products")

    backup_coords(4)
    print()

    """
    STAGE 5
    Refinement of atomic positions using further hard-sphere forces.
    """

    print(highlight_text("Stage 5, Refinement"))

    s5_v_kappa = c["s5_v_kappa"]
    s5_w_kappa = c["s5_w_kappa"]
    s5_hs_kappa = c["s5_hs_kappa"]
    s5_z_kappa = c["s5_z_kappa"]
    s5_trans = c["s5_trans"]
    s5_rms_force = c["s5_rms_force"]

    vr_trans_torque = get_vr_trans_torque(kappa=s5_v_kappa, do_trans=s5_trans)
    wr_trans_torque = get_wr_trans_torque(kappa=s5_w_kappa, do_trans=s5_trans)
    zr_aa_trans_torque = AtomAtomTransTorque(runion, rfrag_lists, AR, kappa=s5_z_kappa)
    r_keys_calcs = {
        "v": vr_trans_torque,
        "w": wr_trans_torque,
        "hardsphere": HardSphere(runion, rfrag_lists, kappa=s5_hs_kappa),
        "z": zr_aa_trans_torque,
    }
    composite_sd_opt(runion, r_keys_calcs, "Reactants", rms_force=s5_rms_force)

    vp_trans_torque = get_vp_trans_torque(kappa=s5_v_kappa, do_trans=s5_trans)
    wp_trans_torque = get_wp_trans_torque(kappa=s5_w_kappa, do_trans=s5_trans)
    zp_aa_trans_torque = AtomAtomTransTorque(punion, pfrag_lists, AP, kappa=s5_z_kappa)
    p_keys_calcs = {
        "v": vp_trans_torque,
        "w": wp_trans_torque,
        "hardsphere": HardSphere(punion, pfrag_lists, kappa=s5_hs_kappa),
        "z": zp_aa_trans_torque,
    }
    composite_sd_opt(punion, p_keys_calcs, "Products", rms_force=s5_rms_force)

    backup_coords(5)
    print()

    with open(make_fn("s5.trj"), "w") as handle:
        handle.write("\n".join([geom.as_xyz() for geom in (runion, punion)]))

    def dump_stages(fn, atoms, coords_list):
        align_coords(coords_list)
        comments = [f"Stage {i}" for i in range(coords_list.shape[0])]
        fn = make_fn(fn)
        coords_to_trj(fn, atoms, coords_list, comments=comments)

    dump_stages("r_coords.trj", runion.atoms, r_coords)
    dump_stages("p_coords.trj", punion.atoms, p_coords)

    runion.set_calculator(None)
    punion.set_calculator(None)
    return runion, punion
