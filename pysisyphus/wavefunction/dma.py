# [1] https://doi.org/10.1080/00268970110089432
#     Distributed multipole analysis Methods and applications
#     Stone, Alderton, 1985
# [2] https://doi.org/10.1016/0009-2614(81)85452-8
#     Distributed multipole analysis, or how to describe a molecular
#     charge distribution
#     Stone, 1981
# [3] https://doi.org/10.1063/1.3430523
#     Convergence of the multipole expansion for 1,2 Coulomb interactions:
#     The modified multipole shifting algorithm
#     Solano, Pendas, Francisco, Blanco, Popelier, 2010
# [4] https://doi.org/10.1063/5.0076630
#     Multi-center decomposition of molecular densities:
#     A mathematical perspective
#     Benda, Cances, Ehrlacher, Stamm, 2022


import itertools as it
import textwrap
import time
from typing import Optional, Tuple

import jinja2
import numpy as np
from scipy.special import binom

from pysisyphus.constants import ANG2BOHR, BOHR2ANG
from pysisyphus.elem_data import ATOMIC_NUMBERS, nuc_charges_for_atoms
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.numint import get_mol_grid, get_gdma_atomic_grid
from pysisyphus.wavefunction import gdma_int
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.ints.multipole3d_sph import multipole3d_sph
from pysisyphus.wavefunction.helpers import lm_iter
from pysisyphus.wavefunction.solid_harmonics import Rlm
from pysisyphus.xyzloader import make_xyz_str


_LE_MAX = 2  # Up to quadrupoles
_TWO_LE_MAX = 2 * _LE_MAX
_SQRT2 = 1 / np.sqrt(2.0)
_DIST_THRESH = 1e-6
_SQRT3 = np.sqrt(3.0)


def get_index_maps(L_max: int = _LE_MAX) -> dict[tuple[int, int], int]:
    """Map between (l, m) keys and 1d indices in the multipole arrays.

    The ordering here must match the ordering in the integral functions.
    """

    index_map = dict()
    index = 0
    for l in range(L_max + 1):
        for m in range(-l, l + 1):
            key = (l, m)
            index_map[key] = index
            index += 1
    return index_map


# Map between 2d (l, m) and 1d-indices.
_INDEX_MAP = get_index_maps()


def get_C(lmax: int = _TWO_LE_MAX) -> dict[tuple[int, int], complex]:
    """Factors for complex->real transformation of multipoles.

    Defined below A9 in [3].
    """
    C = dict()
    for ma in range(-lmax, lmax + 1):
        for mb in range(-lmax, lmax + 1):
            key = (ma, mb)
            if key == (0, 0):
                value = 1.0
            elif abs(ma) != abs(mb):
                value = 0.0
            elif ma > 0 and mb > 0 and ma == mb:
                value = (-1) ** ma * _SQRT2
            elif ma > 0 and (-ma == mb):
                value = _SQRT2
            elif ma < 0 and (-ma == mb):
                value = -1j * (-1) ** mb * _SQRT2
            elif ma < 0 and mb < 0 and ma == mb:
                value = 1j * _SQRT2
            else:
                raise Exception("Safeguard")
            C[key] = value
    return C


C = get_C()
# Complex conjugate values of C
CCONJ = {k: np.conj(v) for k, v in C.items()}


def get_binom_lmkq(L_max: int = _LE_MAX) -> np.ndarray:
    """Bionmial-coefficient prefactor required for translating multipoles.

    See eq. (2.7) in [1] and eq. (11) in [2]. Returns sqrt of the product
    of two binomial coefficients."""
    lmkq = np.zeros((L_max + 1, 2 * L_max + 1, L_max + 1, 2 * L_max + 1))

    for l in range(L_max + 1):
        for m in range(-l, l + 1):
            for k in range(l + 1):
                for q in range(-k, k + 1):
                    binom1 = binom(l + m, k + q)
                    binom2 = binom(l - m, k - q)
                    lmkq[l, m, k, q] = binom1 * binom2
    lmkq = np.sqrt(lmkq)
    return lmkq


def get_W(xyz: np.ndarray, lmkq: np.ndarray, L_max: int = _LE_MAX) -> np.ndarray:
    """Multipole-shift matrix W.

    See eq. (A11) in [3]."""
    num = (np.arange(L_max + 1) * 2 + 1).sum()
    # Complex type is not strictly needed, as the matrix is always real. But adding
    # to a real matrix raises a warning.
    W = np.zeros((num, num), dtype=complex)
    for (l, m), (ldash, mdash) in it.product(lm_iter(L_max), lm_iter(L_max)):
        lm_ind = _INDEX_MAP[(l, m)]
        lmdash_ind = _INDEX_MAP[(ldash, mdash)]
        for m1 in range(-l, l + 1):
            Cmm1 = C[(m, m1)]
            if Cmm1 == 0.0:
                continue
            for m2 in range(-ldash, ldash + 1):
                Cmdashm2 = CCONJ[(mdash, m2)]
                if Cmdashm2 == 0.0:
                    continue
                # prefact = LMKQ[l, m1, ldash, m2]
                prefact = lmkq[l, m1, ldash, m2]
                if prefact == 0.0:
                    continue
                for m3 in range(ldash - l, l - ldash + 1):
                    Cprod = CCONJ[(m3, m1 - m2)]
                    if Cprod == 0.0:
                        continue
                    Cprod *= Cmm1 * Cmdashm2
                    # Index order is consistent with eq. (A11) in [3] (W_{l'm', lm}).
                    W[lmdash_ind, lm_ind] += prefact * Cprod * Rlm(l - ldash, m3, *xyz)
    return np.real(W)


def move_multipoles(
    L_max: int,
    lmkq: np.ndarray,
    multipoles: np.ndarray,
    P: np.ndarray,
    S: np.ndarray,
) -> np.ndarray:
    """Move multipoles from P to S using matrix multiplication.

    Follows the approach outlined in [3].
    """
    W = get_W(P - S, lmkq, L_max)
    # Eq. (7) in [3]
    return multipoles @ W


def move_multipoles_gen(
    L_max: int,
    lmkq: np.ndarray,
    Q: np.ndarray,
    P: np.ndarray,
    S: np.ndarray,
) -> np.ndarray:
    """Move multipoles from P to S using matrix multiplication.

    Follows the approach outlined in [3]. The code below is generated
    with 'gen_m2m.py' in the 'scripts' subdirectory. It is basically
    the unrolled matrix multiplication.

    The code below works up to quadrupoles.
    """
    x, y, z = P - S
    # Eq. (7) in [3], unrolled and piped through CSE
    x0 = y * Q[0]
    x1 = z * Q[0]
    x2 = 1.73205080756888 * x
    x3 = 1.73205080756888 * y
    x4 = 1.73205080756888 * x0
    x5 = 1.73205080756888 * Q[2]
    x6 = 1.73205080756888 * z
    x7 = x**2
    x8 = y**2

    return np.array(
        (
            Q[0],
            x0 + Q[1],
            x1 + Q[2],
            x * Q[0] + Q[3],
            x * x4 + x2 * Q[1] + x3 * Q[3] + Q[4],
            x4 * z + x5 * y + x6 * Q[1] + Q[5],
            -x * Q[3]
            - y * Q[1]
            + 2.0 * z * Q[2]
            - (0.5 * x7 + 0.5 * x8 - z**2) * Q[0]
            + Q[6],
            1.73205080756888 * x * x1 + x * x5 + x6 * Q[3] + Q[7],
            x2 * Q[3] - x3 * Q[1] + 0.866025403784438 * (x7 - x8) * Q[0] + Q[8],
        )
    )


def closest_sites_and_weights(
    org_coords: np.ndarray,
    sites: np.ndarray,
    site_radii: np.ndarray,
    thresh: float = _DIST_THRESH,
):
    """Determine closest expansion sites and associated weights.

    Implements Stone's closest-neigbours-take-it-all-strategy.
    Alternatively, Vigne-Maeder could be implemented (eq. (78) in [4])."""
    # Determine distances to expansion sites
    dists = np.linalg.norm(sites - org_coords[None, :], axis=1)
    # Take possibly different site radii into account
    dists /= site_radii

    # Determine closest site
    closest_index = dists.argmin()
    # Look for other sites that are also/equally close
    delta_dists = dists - dists[closest_index]
    delta_dist_mask = np.abs(delta_dists) <= thresh
    # Indices of close neighbours
    neighbours = np.arange(len(delta_dists))[delta_dist_mask]
    # Distribute the multipoles equally among closest sites
    nneighs = len(neighbours)
    weights = np.full(nneighs, 1 / nneighs)
    return neighbours, weights


def accumulate_multipoles(dx, dy, dz, prefact, multipoles):
    dx2 = dx**2
    dy2 = dy**2
    dz2 = dz**2
    # Charges
    multipoles[0] += prefact
    # Dipole moments
    multipoles[1] += prefact * dy
    multipoles[2] += prefact * dz
    multipoles[3] += prefact * dx
    # Quadrupole moments
    multipoles[4] += prefact * _SQRT3 * dx * dy
    multipoles[5] += prefact * _SQRT3 * dy * dz
    multipoles[6] += prefact * (2 * dz2 - dy2 - dx2) / 2.0
    multipoles[7] += prefact * _SQRT3 * dx * dz
    multipoles[8] += prefact * _SQRT3 * (-dy2 + dx2) / 2.0


def get_redistribution_func(sites, site_radii, distributed_multipoles):
    # Prefactors required for redistributing multipoles
    lmkq = get_binom_lmkq()

    def redistribute_multipoles(nat_coords, multipoles):
        neighbours, weights = closest_sites_and_weights(nat_coords, sites, site_radii)
        for neigh_index, weight in zip(neighbours, weights):
            S = sites[neigh_index]
            # Don't move if natural- and expansion-center coincide
            if np.linalg.norm(S - nat_coords) <= _DIST_THRESH:
                moved_multipoles = multipoles
            else:
                moved_multipoles = move_multipoles_gen(
                    _LE_MAX,
                    lmkq,  # Binomial coefficient prefactors
                    weight * multipoles,
                    nat_coords,  # Original/natural center
                    S,  # Target center/DMA site
                )
            # Sum contributions of the different basis function pairs
            distributed_multipoles[neigh_index] += moved_multipoles

    return redistribute_multipoles


def dma(
    wavefunction: Wavefunction,
    sites: np.ndarray,
    site_radii: np.ndarray,
    exp_cutoff: float = -36.0,
    dens_thresh: float = 1e-9,
    switch: float = 4.0,
    addons: bool = True,
):
    assert sites.ndim == 2
    assert site_radii.ndim == 1
    assert len(sites) == len(site_radii)
    assert switch >= 0.0, f"Switch must be a float >= 0.0 but got '{switch}'!"

    # Number of multipoles that will be calculated
    nmultipoles = (2 * np.arange(_LE_MAX + 1) + 1).sum()
    # Number of expansion sites
    nsites = len(sites)
    # Array that will hold the final distributed multipoles
    distributed_multipoles = np.zeros((nsites, nmultipoles))

    # Nuclear coordinates
    coords3d = wavefunction.coords.reshape(-1, 3)
    # Sum alpha and beta densities into total density
    P = np.sum(wavefunction.P, axis=0)

    shells = wavefunction.shells
    reorder_c2s_coeffs = shells.reorder_c2s_coeffs

    # Get function that redistributes multipoles onto the expansion sites
    redistribute_multipoles = get_redistribution_func(
        sites, site_radii, distributed_multipoles
    )

    header_tpl = jinja2.Template(
        textwrap.dedent(
            """
    {{ title }}

    {{ wavefunction }}
    Expansion Sites:
    {{ sites }}
    Switch at Exponent-Sum: {{ switch }}

    """
        )
    )
    header_rendered = header_tpl.render(
        title=highlight_text("Distributed Multipole Analysis"),
        wavefunction=wavefunction,
        sites=make_xyz_str([""] * nsites, sites * BOHR2ANG),
        switch=switch,
    )
    print(header_rendered)

    dma_dur = time.time()
    # Take nuclear charges into account and distribute them onto expansion sites
    nuc_charges = nuc_charges_for_atoms(wavefunction.atoms)
    for i, nuc_charge in enumerate(nuc_charges):
        nuc_multipoles = np.zeros(nmultipoles)
        nuc_multipoles[0] = nuc_charge
        nuc_coords = coords3d[i]
        redistribute_multipoles(nuc_coords, nuc_multipoles)

    # Loop over unique shell pairs
    for i, shell_a in enumerate(shells):
        La, A, das, axs, a_ind, a_size = shell_a.as_tuple()
        a_slice_cart = slice(a_ind, a_ind + a_size)
        a_slice_sph = shell_a.sph_slice
        a_c2s = reorder_c2s_coeffs[a_slice_sph, a_slice_cart]
        for j, shell_b in enumerate(shells[i:], i):
            Lb, B, dbs, bxs, b_ind, b_size = shell_b.as_tuple()
            b_slice_sph = shell_b.sph_slice
            # Slice of the density matrix
            P_slice = P[a_slice_sph, b_slice_sph]
            # Every off-diagonal element of the density matrix contributes two times.
            sym_factor = 2.0 if (i != j) else 1.0

            # Screen out whole shell pairs
            #
            # The default dens_thresh of 1e-9 is very conservative and will
            # probably never trigger for smaller molecules. Anything above 1e-9,
            # e.g. 1e-8 already leads to quadrupole errors for methane (HF/def2-svp)
            # on the order of 1e-8.
            dens_max = max(abs(P_slice.min()), P_slice.max())
            dens_abs = sym_factor * dens_max
            if dens_abs <= dens_thresh:
                continue

            AB = A - B
            AB_dist2 = AB.dot(AB)
            b_slice_cart = slice(b_ind, b_ind + b_size)
            b_c2s = reorder_c2s_coeffs[b_slice_sph, b_slice_cart]
            # Pick appropriate integral function
            func = multipole3d_sph[(La, Lb)]

            # Loop over primitive pairs.
            for da, ax in zip(das, axs):
                ax_AB_dist = ax * AB_dist2
                for db, bx in zip(dbs, bxs):
                    px = ax + bx

                    # Do numerical integration otherwise
                    if px <= switch:
                        continue

                    # exp-argument of gaussian overlap
                    exp_arg = (bx * ax_AB_dist) / px
                    # When exp_arg in 'exp(-exp_arg)' is very small we can skip this
                    # pair of primitives. The default thresh is -36 (exp(-36) = 2.32e-16),
                    # so the screening is very conservative. Using 'exp_cutoff = 10' results
                    # in errors up to 1e-6 in the test cases.
                    if -exp_arg <= exp_cutoff:
                        continue

                    # Natural center of Gaussian primitive pair
                    P_prim = (ax * A + bx * B) / (ax + bx)

                    # Calculate multipole integrals.
                    # The minus is taken into account here, as well as the symmetry factor.
                    prim_multipoles = -sym_factor * func(ax, da, A, bx, db, B)
                    # Convert to spherical basis functions & reorder, so the basis function
                    # order matches the density matrix.
                    prim_multipoles = np.einsum(
                        "ij,mjk,kl->mil",
                        a_c2s,
                        prim_multipoles,
                        b_c2s.T,
                        optimize="greedy",
                    )

                    # Contract integrals with density matrix elements.
                    contr_multipoles = prim_multipoles * P_slice[None, :, :]
                    contr_multipoles = contr_multipoles.sum(axis=(1, 2))
                    # Redistribute onto expansion sites
                    redistribute_multipoles(P_prim, contr_multipoles)
            # End Loop over primitive pairs
    dma_dur = time.time() - dma_dur
    print(f"DMA part took {dma_dur:.4f} s")

    # Handle numerical integration/do GDMA
    if switch > 0.0:
        # Defaul to the numba implementation
        diffuse_density_func = gdma_int.get_diffuse_density_numba

        # Try to import fater Fortran routines, if enabled & available
        if addons:
            try:
                from pysisyphus_addons.gdma.gdma_int import (
                    # Overwrite previously assigned numba function
                    get_diffuse_density as diffuse_density_func,
                )
            except ModuleNotFoundError:
                pass

        # Set up molecular grid for numerical integration ...
        grid_dur = time.time()
        # TODO: make atom_grid_kwargs adjustable
        mol_grid = get_mol_grid(wavefunction, atom_grid_kwargs={"kind": "g5"})
        # mol_grid = get_mol_grid(wavefunction, grid_func=get_gdma_atomic_grid)
        grid_dur = time.time() - grid_dur
        print(f"Grid generation took {grid_dur:.4f} s")
        print(f"There are {mol_grid.npoints} grid points")

        # ... and accumulate the density on the grid
        numint_dur = time.time()
        rho = diffuse_density_func(wavefunction, mol_grid, switch)
        numint_dur = time.time() - numint_dur
        print(f"Numerical integration took {numint_dur:.4f} s")

        for i, (grid_center, xyz, ww) in enumerate(mol_grid.grid_iter()):
            np.testing.assert_allclose(grid_center, sites[i])
            slice_ = mol_grid.slices[i]
            multipoles = np.zeros(nmultipoles)
            # Take minus sign into account here
            ww_rho = -ww * rho[slice_]
            for xyz_, wrho_ in zip(xyz[::-1], ww_rho[::-1]):
                dist = xyz_ - grid_center
                accumulate_multipoles(*dist, wrho_, multipoles)

            redistribute_multipoles(grid_center, multipoles)
    return distributed_multipoles


def get_radii(atoms: Tuple[str]) -> np.ndarray:
    """Default site radii from DMA. 0.325 Å for Hydrogen, 0.65 Å else."""
    radii = [0.325 if atom.lower() == "h" else 0.65 for atom in atoms]
    radii = np.array(radii) / BOHR2ANG
    return radii


def run_dma(
    wf: Wavefunction,
    sites: Optional[np.ndarray] = None,
    site_labels: list[str] = None,
    site_radii: Optional[np.ndarray] = None,
    **kwargs,
):
    if sites is not None:
        assert (
            site_radii is not None
        ), "When sites are given their radii must also be provided!"

    if sites is None:
        dma_sites = wf.coords.reshape(-1, 3).copy()
        site_labels = wf.atoms
    else:
        dma_sites = np.array(sites).reshape(-1, 3)
    if site_labels is None:
        site_labels = ()

    # TODO: handling radius assignment when sites are given ...
    if site_radii is None:
        site_radii = get_radii(wf.atoms)

    assert len(site_radii) == len(dma_sites)

    start = time.time()
    distributed_multipoles = dma(wf, dma_sites, site_radii, **kwargs)
    dur = time.time() - start
    print(f"Total analysis took {dur:.4f} s\n")

    cfmt = " >10.6f"
    mfmt = " >10.6f"
    for i, (lbl, dma_site) in enumerate(zip(site_labels, dma_sites)):
        x, y, z = dma_site
        radius = site_radii[i]
        print(
            f"{lbl.capitalize()}: "
            f"({x=:{cfmt}},  {y=:{cfmt}}, {z=:{cfmt}}, {radius=:{cfmt}}) au"
        )
        site_multipoles = distributed_multipoles[i]
        index = 0
        for l in range(_LE_MAX + 1):
            size = 2 * l + 1
            l_multipoles = site_multipoles[index : index + size]
            l_norm = np.linalg.norm(l_multipoles)
            multi_strs = [
                f"{sm:{mfmt}}" for sm in site_multipoles[index : index + size]
            ]
            print(f"\t|Q_{l}|={l_norm:{mfmt}}: {', '.join(multi_strs)}")
            index += size
        print()
    return distributed_multipoles


ORIENT_TPL = jinja2.Template(
    """
title {{ title }}
default rank {{ rank }}

units bohr

types
{% for atom, charge in atoms_charges -%}
{{ atom }} Z {{ charge }}
{% endfor -%}
end

molecule {{ molecule_name }} at 0.0 0.0 0.0
{%- for atom, (x, y, z), atom_type, bonds, multipoles in molecule %}
{{ atom }} {{ ffmt(x) }} {{ ffmt(y) }} {{ ffmt(z) }} type {{ atom_type }} rank {{ rank }}
{% for rank in range(rank+1) -%}
    {% for rmp in multipoles[rank**2:(rank+1)**2] %}{{ ffmt(rmp) }} {% endfor %}
{% endfor -%}
{% endfor -%}
end

edit {{ molecule_name }}
bonds auto
end

Display potential
  Tracing off
  Grid
    name {{ molecule_name }}-vdw2
    Molecule {{ molecule_name }}
    Radii scale 2
    Step 0.25 B
  End
  Map vdw2
    Molecule {{ molecule_name }}
    Potential
  End

  Background RGBX e8e8e8
  Colour-map
    0    210  0.25  1
    6    240  0.75 1
    12   300  1.0  0
    18   360  0.75 1
    24   390  0.25  1
  End

  Show vdw2
    contour eV -{{ eV_lim}} 0.0 {{ eV_lim }}
    Colour-scale min -{{ eV_lim }} max {{ eV_lim }} eV
    ticks -{{ eV_lim }} 0.0 {{ eV_lim }}
    Ball-and-stick
  End
End

Comment "Display closed"

Finish"""
)


def to_orient(
    atoms,
    coords3d,
    multipoles,
    rank=2,
    molecule_name="Molecule1",
    title="Orient input from pysisyphus",
    eV_lim=0.5,
):
    atoms = [atom.lower() for atom in atoms]
    charges = [ATOMIC_NUMBERS[atom] for atom in atoms]
    atoms_charges = zip(atoms, charges)

    molecule = list()
    assert rank == 2
    order = [0, 2, 3, 1, 6, 7, 5, 8, 4]
    for i, atom in enumerate(atoms):
        atom_coords = coords3d[i]
        bonds = ()  # Currently unused
        atom_multipoles = multipoles[i][order]
        # Resort multipoles!
        atom_type = atom
        molecule.append((atom, atom_coords, atom_type, bonds, atom_multipoles))

    def ffmt(num):
        return f"{num: >14.8f}"

    rendered = ORIENT_TPL.render(
        title=title,
        rank=rank,
        atoms_charges=atoms_charges,
        molecule_name=molecule_name,
        molecule=molecule,
        eV_lim=eV_lim,
        ffmt=ffmt,
    )
    return rendered
