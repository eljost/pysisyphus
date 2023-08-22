# [1] https://doi.org/10.26434/chemrxiv-2022-253hc-v2
#     Identifying the Marcus dimension of electron transfer from
#     ab initio calculations
#     Šrut, Lear, Krewald, 2022 on chemrxiv
# [2] https://doi.org/10.1039/D3SC01402A
#     Identifying the Marcus dimension of electron transfer from
#     ab initio calculations
#     Šrut, Lear, Krewald, 2023, actual published version

import datetime
from enum import Enum
from pathlib import Path
import sys
import time
from typing import Callable, List, Tuple

from distributed import Client
import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.dynamics import get_wigner_sampler
from pysisyphus.helpers import get_tangent_trj_str, get_fragment_xyzs, highlight_text
from pysisyphus.helpers_pure import approx_float, eigval_to_wavenumber, rms
from pysisyphus.io import geom_from_hessian
from pysisyphus.linalg import are_collinear
from pysisyphus.TablePrinter import TablePrinter
from pysisyphus.wavefunction.pop_analysis import iao_charges_from_wf
from pysisyphus.xyzloader import make_xyz_str


_CORR_THRESH = 0.10


########################################
# Plot & report expansion coefficients #
########################################


def create_marcus_coefficient_plot(nqs, coeffs, batch=None, ncalcs=None):
    fig, ax = plt.subplots()
    ax.bar(np.arange(nqs) + 6, coeffs)
    ax.set_xlabel("Mode")
    ax.set_ylabel("Coefficient")
    title_frags = list()
    if batch is not None:
        title_frags.append(f"Batch {batch}")
    if ncalcs is not None:
        title_frags.append(f"{ncalcs} calculations")
    title = ", ".join(title_frags)
    fig.suptitle(title)
    return fig, ax


def report_marcus_coefficients(coeffs, nus, drop_first, nhighest=10):
    breakpoint()
    highest_bs = np.abs(coeffs).argsort()[-nhighest:][::-1]
    coeff_table = TablePrinter(
        header=("internal", "0-based", "1-based", "coeff. b", "ṽ / cm⁻¹"),
        col_fmts=("int", "int", "int", "float", "float_short"),
    )
    print("10 highest (absolute) Normal mode expansion coefficients:")
    coeff_table.print_header()
    # for i, hb in enumerate(highest_bs):
    for hb in highest_bs:
        coeff_table.print_row(
            (hb, hb + drop_first, hb + drop_first + 1, coeffs[hb], nus[hb])
        )


def report_marcus_coefficients(coeffs, nus, drop_first, nhighest=10):
    sort_inds = np.abs(coeffs).argsort()
    highest_inds = sort_inds[-nhighest:][::-1]
    coeff_table = TablePrinter(
        header=("internal", "0-based", "1-based", "coeff. b", "ṽ / cm⁻¹"),
        col_fmts=("int", "int", "int", "float", "float_short"),
    )
    print("10 highest (absolute) Normal mode expansion coefficients:")
    coeff_table.print_header()
    for hi in highest_inds:
        coeff_table.print_row(
            (hi, hi + drop_first, hi + drop_first + 1, coeffs[hi], nus[hi])
        )


##############################
# Plot & report correlations #
##############################


def create_correlation_plots(
    qs, nus, corrs, depos, drop_first, batch, out_dir=".", corr_thresh=_CORR_THRESH
):
    out_dir = Path(out_dir)
    fig, ax = plt.subplots()
    # Normal coordinate limits to get sensible y-axis limits
    deposmin = depos.min()
    deposmax = depos.max()
    deposfit = np.linspace(deposmin, deposmax, 2)
    qmin = qs.min()
    qmax = qs.max()
    print(f"Correlation plots with ρ >= {corr_thresh:.6f}:")
    # One plot per normal coordinate
    nplots = 0
    for i, (q, nu) in enumerate(zip(qs.T, nus)):
        corr = corrs[i]
        if abs(corr) < corr_thresh:
            continue
        print(f"\t... processing mode {i:03d}, ρ={corr: >12.6f}")
        ax.scatter(depos, q)
        # Linear fit
        coef = np.polynomial.polynomial.polyfit(depos, q, deg=1, full=False)
        poly = np.polynomial.polynomial.Polynomial(coef, domain=(deposmin, deposmax))
        ax.plot(deposfit, poly(deposfit), c="red")
        ax.set_xlim(deposmin, deposmax)
        ax.set_ylim(qmin, qmax)
        ax.set_xlabel("epos")
        ax.set_ylabel(f"$q_{{{i+drop_first}}}$")
        title = rf"Batch {batch:02d}: Mode {i:03d}, {nu: >8.2f} cm⁻¹, $\rho$ = {corr:> 8.4f}"
        ax.set_title(title)
        fn = f"correlation_batch_{batch:02d}_mode_{i:03d}.png"
        fig.savefig(out_dir / fn)
        ax.cla()
        nplots += 1
    print(f"Created {nplots} correlation plots.\n")


def report_correlations(nus, corrs, drop_first, corr_thresh=_CORR_THRESH):
    corrs_above = np.where(np.abs(corrs) >= corr_thresh)[0]
    corrs_above_sorted = np.abs(corrs[corrs_above]).argsort()[::-1]
    nabove = len(corrs_above)
    print(f"Found {nabove} modes with correlation coefficients >= {corr_thresh:.4f}.\n")

    corr_table = TablePrinter(
        header=("internal", "0-based", "1-based", "corr. ρ", "ṽ / cm⁻¹"),
        col_fmts=("int", "int", "int", "float", "float_short"),
    )
    corr_table.print_header()
    for cas in corrs_above_sorted:
        ca = corrs_above[cas]
        corr_table.print_row(
            (ca, ca + drop_first, ca + drop_first + 1, corrs[ca], nus[ca])
        )
    print()


def get_marcus_dim(eigvecs, masses, qs, depos):
    """Determine Marcus dimension from normal coordinates and properties.

    Parameters
    ----------
    eigevecs
        2d array with shape (3N, nmodes) holding normal modes, with N being the
        number of atoms.
    masses
        1d array with shape (N, ) holding atomic masses.
    qs
        2d array with shape (nsamples, nmodes) holding the normal coordinates
        for each sample.
    depos
        1d array with shape (nsamples, ) holding the differences in electron
        positions from the equilbrium geometries.

    Returns
    -------
    corrs
        Pearson correlation coefficients.
    coeffs
        Normal mode expansion coefficients.
    marcus_dim
        Normalized Marcus dimension in unweighted Cartesiand coordinates
        of shape (3*natoms, ).
    """
    sqrt_masses_rep = np.repeat(np.sqrt(masses), 3)

    # Number of normal modes
    nmodes = eigvecs.shape[1]

    corrs = np.zeros(nmodes)
    for i, q in enumerate(qs.T):
        corr = np.corrcoef(depos, q)[0, 1]
        corrs[i] = corr

    # Do least-squares.
    print("Solving least squares 'qs . coeffs = depos' for 'coeffs'")
    # Residuals are an empty array
    coeffs, _, rank, singvals = np.linalg.lstsq(qs, depos, rcond=None)

    # Construct Marcus dimension as linear combination of normal modes
    # according to coefficients b from least-squares fit.
    marcus_dim = np.einsum("i,ki->k", coeffs, eigvecs)
    # Remove mass-weighting
    marcus_dim /= sqrt_masses_rep
    # and renormalize vector.
    marcus_dim /= np.linalg.norm(marcus_dim)
    return corrs, coeffs, marcus_dim


Property = Enum("Property", ["EPOS", "EEXC"])


def get_displaced_coordinates(wigner_sampler, eigvecs, coords_eq, M):
    """Displaced (normal) coordinates from Wigner sampling."""
    # Displaced coordinates, drop velocities
    displ_coords, _ = wigner_sampler()
    displ_coords = displ_coords.flatten()
    # Calculate mass-weighted displacements from equilibrium geometry and
    # transform to normal coordinates.
    displ = (displ_coords - coords_eq).flatten()
    displ_mw = M @ displ
    # Transform mass-weighted displacements to normal coordinates
    norm_coords = eigvecs.T @ displ_mw
    return displ_coords, norm_coords


def epos_from_wf(wf, fragments):
    """Fragment mulliken charges using intrinsic atomic orbital."""
    pop_ana = iao_charges_from_wf(wf)

    assert approx_float(pop_ana.tot_charge, wf.charge)

    def sum_spin_pop(spin_pop):
        per_frags = [spin_pop[frag] for frag in fragments]
        frag_sums = [pf.sum() for pf in per_frags]
        weights = np.arange(1, len(frag_sums) + 1)
        # Centroid of spin density; basically eq. (3) from [2].
        epos = (weights * frag_sums).sum()
        return frag_sums, epos

    _, tot_epos = sum_spin_pop(pop_ana.spin_pop)  # alpha + beta
    _, alpha_epos = sum_spin_pop(pop_ana.alpha_spin_pop)  # alpha only
    return tot_epos, alpha_epos


def epos_property(geom, fragments, eq_property):
    wf = geom.calculator.get_stored_wavefunction()
    # Only use alpha part
    # TODO: make this toggable?
    _, alpha_epos = epos_from_wf(wf, fragments)
    property = alpha_epos - eq_property
    return property


def en_exc_property(geom, fragments, eq_property):
    gs_energy, *es_energies = geom.all_energies
    assert len(es_energies) > 0
    property = es_energies[0] - gs_energy
    return property


def batched_marcus_dim(
    h5_fn: str,
    calc_getter: Callable,
    fragments: List[Tuple[int]],
    T: float = 300,
    property=Property.EPOS,
    batch_size: int = 25,
    max_batches: int = 20,
    rms_thresh: float = 0.005,
    correlations: bool = False,
    corr_thresh: float = _CORR_THRESH,
    scheduler=None,
    out_dir=".",
):
    assert T >= 0.0
    assert batch_size > 0
    assert max_batches > 0
    assert rms_thresh > 0.0
    max_ncalcs = batch_size * max_batches

    property_funcs = {
        property.EPOS: epos_property,
        property.EEXC: en_exc_property,
    }
    property_func = property_funcs[property]

    out_dir = Path(out_dir)

    # Load geometry from pysisyphus HDF5 Hessian
    geom, h5_attrs = geom_from_hessian(h5_fn, with_attrs=True)
    try:
        charge = h5_attrs["charge"]
        mult = h5_attrs["mult"]
    except KeyError:
        charge = None
        mult = None

    if (charge is not None) and (mult is not None):
        print(f"       Charge: {charge}")
        print(f" Multiplicity: {mult}")
    print(f"  Temperature: {T:.2f} K")
    print(f"   Batch size: {batch_size}")
    print(f" Max. batches: {max_batches}")
    print(f"rms threshold: {rms_thresh}")
    print(f"     Property: {property}")

    print(f"Doing at most {batch_size*max_batches} calculations")

    # Dump fragments
    fragments_trj = "\n".join(
        get_fragment_xyzs(geom, fragments, with_geom=True, with_dummies=True)
    )
    fragment_trj_fn = "fragments.trj"
    with open(fragment_trj_fn, "w") as handle:
        handle.write(fragments_trj)
    print(f"Wrote fragments to '{fragment_trj_fn}'")

    coords_eq = geom.cart_coords
    # The function will probably never be applied to linear molecules,
    # but why not check this.
    linear = are_collinear(geom.coords3d)
    drop_first = 5 if linear else 6
    eigvals, eigvecs = np.linalg.eigh(
        geom.eckart_projection(geom.mw_hessian, full=True)
    )
    nus = eigval_to_wavenumber(eigvals)
    eigvecs = eigvecs[:, drop_first:]
    nus = nus[drop_first:]
    nmodes = len(nus)
    geom.set_calculator(calc_getter(0, base_name="eq"))

    # Calculate wavefunction at equilibrium geometry
    print("Starting calculation of equilibrium geometry")
    geom.all_energies
    print("Finished calculation of equilibrium geometry")
    property_eq = property_func(geom, fragments, 0.0)

    # Function that create displaced geometries by drawing from a Wigner distribution
    wigner_sampler = get_wigner_sampler(h5_fn, temperature=T)

    # Arrays holding normal coordinates and properties
    all_norm_coords = np.zeros((max_ncalcs, nmodes))
    all_properties = np.zeros(max_ncalcs)

    to_save = {
        # Property key
        "property": str(property),
        # Equilibrium geometry
        "cart_coords_eq": geom.cart_coords,
        "property_eq": property_eq,
        # Samples
        "normal_coordinates": all_norm_coords,
        "properties": all_properties,
    }

    masses = geom.masses
    sqrt_masses_rep = np.repeat(np.sqrt(masses), 3)
    M = np.diag(1 / sqrt_masses_rep)
    prev_coeffs = None
    prev_marcus_dim = None
    rms_coeffs = None
    rms_marcus_dim = None

    def calculate_property(i):
        displ_geom = geom.copy()
        displ_geom.coords = batch_displ_coords[i - start_ind]
        displ_geom.set_calculator(calc_getter(i))
        displ_geom.all_energies
        property = property_func(displ_geom, fragments, property_eq)
        return property

    if scheduler is not None:
        client = Client(scheduler)
    else:
        client = None

    print()
    for batch in range(max_batches):
        batch_str = f"batch_{batch:02d}"
        start_ind = batch * batch_size
        end_ind = start_ind + batch_size
        now = datetime.datetime.now()
        print(highlight_text(f"Batch {batch}") + "\n")
        print(f"Starting at index {start_ind} on {now.strftime('%c')}")
        print(f"Doing calculations with indices {start_ind} to {end_ind-1}.")

        batch_displ_coords = np.zeros((batch_size, coords_eq.size))
        # Get and store displaced coordinates and normal coordinates
        for i in range(start_ind, end_ind):
            displ_coords, norm_coords = get_displaced_coordinates(
                wigner_sampler, eigvecs, coords_eq, M
            )
            batch_displ_coords[i - start_ind] = displ_coords
            all_norm_coords[i] = norm_coords

        # Loop over cacluclations in current batch
        batch_start = time.time()
        batch_range = range(start_ind, end_ind)
        if client is not None:
            futures = client.map(calculate_property, batch_range)
            batch_properties = client.gather(futures)
            all_properties[start_ind:end_ind] = batch_properties
        else:
            for i in batch_range:
                if i % 5 == 0:
                    print(f"\t{i}")
                    sys.stdout.flush()
                all_properties[i] = calculate_property(i)
        # End loop over calculations in one batch

        # Dump results
        np.savez("marcus_dim_results.npz", **to_save)

        batch_dur = time.time() - batch_start
        calc_dur = batch_dur / batch_size
        print(f"Did {batch_size} calculations.")
        print(f"Full batch took {batch_dur:.2f} s, {calc_dur:.2f} s / calculation")
        print(f"Total number of calculations done until now: {end_ind}")

        # Actually calculate Marcus dimension using least-squares
        corrs, coeffs, marcus_dim = get_marcus_dim(
            eigvecs, masses, all_norm_coords[:end_ind], all_properties[:end_ind]
        )

        # XYZ representation of Marcus dimension and animation
        marcus_dim_xyz_str = make_xyz_str(geom.atoms, marcus_dim.reshape(-1, 3))
        print(f"Current Marcus dimension:\n{marcus_dim_xyz_str}\n")
        marcus_dim_trj = get_tangent_trj_str(
            geom.atoms, geom.cart_coords, marcus_dim, comment="Marcus dimension"
        )
        # Dump animated Marcus dimension
        with open(f"marcus_dim_{batch_str}.trj", "w") as handle:
            handle.write(marcus_dim_trj)

        # Report expansion coefficients and rms values
        report_marcus_coefficients(coeffs, nus, drop_first)
        print()
        if prev_coeffs is not None:
            rms_coeffs = rms(prev_coeffs - coeffs)
            rms_marcus_dim = rms(prev_marcus_dim - marcus_dim)
            print(
                f"rms(Δcoeffs)={rms_coeffs:.6f}, rms(ΔMarcus dim)={rms_marcus_dim:.6f}"
            )

        # Plot of expansion coefficients
        fig, _ = create_marcus_coefficient_plot(
            nmodes, coeffs, batch=batch, ncalcs=end_ind
        )
        plot_fn = out_dir / f"marcus_coeffs_{batch_str}.png"
        fig.savefig(plot_fn)
        print(f"Saved expansion coefficent plot to '{plot_fn}'.\n")

        # Correlations between normal modes and properties
        if correlations:
            report_correlations(nus, corrs, drop_first, corr_thresh=corr_thresh)
            create_correlation_plots(
                all_norm_coords[:end_ind],
                nus,
                corrs,
                all_properties[:end_ind],
                drop_first,
                batch=batch,
                out_dir=out_dir,
                corr_thresh=corr_thresh,
            )

        if (rms_marcus_dim is not None) and (rms_marcus_dim <= rms_thresh):
            print(
                f"\nCalculation of Marcus dimension converged after {end_ind} calculations!"
            )
            break
        prev_coeffs = coeffs
        prev_marcus_dim = marcus_dim

        print()
    return marcus_dim
