# [1] https://doi.org/10.26434/chemrxiv-2022-253hc-v2
#     Identifying the Marcus dimension of electron transfer from
#     ab initio calculations
#     Šrut, Lear, Krewald, 2022 on chemrxiv
# [2] https://doi.org/10.1039/D3SC01402A
#     Identifying the Marcus dimension of electron transfer from
#     ab initio calculations
#     Šrut, Lear, Krewald, 2023, actual published version

import datetime
import functools
from pathlib import Path
import sys
import time
from typing import Callable, List, Optional, Tuple

from distributed import Client
import matplotlib.pyplot as plt
import numpy as np

from pysisyphus import plot_helpers
from pysisyphus.constants import BOHR2ANG, C
from pysisyphus.exceptions import (
    CalculationFailedException,
    RunAfterCalculationFailedException,
)
from pysisyphus.marcus_dim.config import RMS_THRESH, FIT_RESULTS_FN
import pysisyphus.marcus_dim.types as mdtypes
from pysisyphus.dynamics.wigner import get_wigner_sampler, plot_normal_coords
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import (
    get_tangent_trj_str,
    get_fragment_xyzs,
    check_for_end_sign,
)
from pysisyphus.helpers_pure import (
    approx_float,
    eigval_to_wavenumber,
    highlight_text,
    rms,
)
from pysisyphus.TablePrinter import TablePrinter
from pysisyphus.wavefunction.pop_analysis import (
    iao_charges,
    mulliken_charges,
)
from pysisyphus.xyzloader import make_xyz_str


_CORR_THRESH = 0.10


########################################
# Plot & report expansion coefficients #
########################################


def create_marcus_coefficient_plot(nqs, coeffs, batch=None, ncalcs=None):
    fig, ax = plot_helpers.fig_ax_from_canvas_agg()
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
    sort_inds = np.abs(coeffs).argsort()
    highest_inds = sort_inds[-nhighest:][::-1]
    coeff_table = TablePrinter(
        header=("internal", "0-based", "1-based", "coeff. b", "ṽ / cm⁻¹"),
        col_fmts=("int_short", "int_short", "int_short", "float", "float_short"),
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
    fig, ax = plot_helpers.fig_ax_from_canvas_agg()
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
        ax.set_xlabel("Δepos")
        ax.set_ylabel(f"$q_{{{i + drop_first}}}$")
        title = rf"Batch {batch:02d}: Mode {i:03d}, {nu: >8.2f} cm⁻¹, $\rho$ = {corr:> 8.4f}"
        ax.set_title(title)
        fn = f"correlation_batch_{batch:02d}_mode_{i:03d}.svg"
        fig.savefig(out_dir / fn)
        ax.cla()
        nplots += 1
    print(f"Created {nplots} correlation plots.\n")
    plt.close(fig)


def report_correlations(nus, corrs, drop_first, corr_thresh=_CORR_THRESH):
    corrs_above = np.where(np.abs(corrs) >= corr_thresh)[0]
    corrs_above_sorted = np.abs(corrs[corrs_above]).argsort()[::-1]
    nabove = len(corrs_above)
    print(f"Found {nabove} modes with correlation coefficients >= {corr_thresh:.4f}.\n")

    corr_table = TablePrinter(
        header=("internal", "0-based", "1-based", "corr. ρ", "ṽ / cm⁻¹"),
        col_fmts=("int_short", "int_short", "int_short", "float", "float_short"),
    )
    corr_table.print_header()
    for cas in corrs_above_sorted:
        ca = corrs_above[cas]
        corr_table.print_row(
            (ca, ca + drop_first, ca + drop_first + 1, corrs[ca], nus[ca])
        )
    print()


def get_marcus_dim(cart_displs, qs, properties):
    """Determine Marcus dimension from normal coordinates and properties.

    Parameters
    ----------
    cart_displs
        2d array with shape (3N, nmodes) holding normal modes, with N being the
        number of atoms.
    masses
        1d array with shape (N, ) holding atomic masses.
    qs
        2d array with shape (nsamples, nmodes) holding the normal coordinates
        for each sample.
    properties
        1d array with shape (nsamples, ) holding the calculated properties, e.g,.
        electron position (centroid of spin density) or exciation energy.

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

    # Only use non-NaN values
    mask = ~np.isnan(properties)
    qs = qs[mask]
    properties = properties[mask]

    # Number of normal modes
    nmodes = cart_displs.shape[1]

    corrs = np.zeros(nmodes)
    for i, q in enumerate(qs.T):
        corr = np.corrcoef(properties, q)[0, 1]
        corrs[i] = corr

    # Do least-squares.
    print("Solving least squares 'qs . coeffs = depos' for 'coeffs'")
    # coeffs, residual, rank, singvals = np.linalg.lstsq(...)
    coeffs, *_ = np.linalg.lstsq(qs, properties, rcond=None)

    # Construct Marcus dimension as linear combination of normal modes
    # according to coefficients b from least-squares fit ...
    marcus_dim = np.einsum("i,ki->k", coeffs, cart_displs)
    # ... and renormalize vector.
    marcus_dim /= np.linalg.norm(marcus_dim)

    # Marcus dimension in normal coordinates.
    # Normalize coefficients for convenience
    marcus_dim_q = coeffs / np.linalg.norm(coeffs)

    return corrs, coeffs, marcus_dim, marcus_dim_q


def get_displaced_coordinates(wigner_sampler, cart_displs, coords_eq):
    """Displaced (normal) coordinates from Wigner sampling."""
    # Sample from wigner distribution
    sample = wigner_sampler()
    displ_coords = sample.coords3d.flatten()
    # Calculate displacements from equilibrium geometry and transform to normal coordinates.
    norm_coords = cart_displs.T @ (displ_coords - coords_eq).flatten()
    return displ_coords, norm_coords


POP_FUNCS = {
    mdtypes.PopKind.IAO: iao_charges,
    mdtypes.PopKind.MULLIKEN: mulliken_charges,
}


def epos_from_wf(wf, fragments, pop_kind: mdtypes.PopKind = mdtypes.PopKind.IAO):
    """Fragment mulliken charges using intrinsic atomic orbital."""
    pop_func = POP_FUNCS[pop_kind]
    pop_ana = pop_func(wf)

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


def epos_property(geom, fragments, pop_kind):
    wf = geom.calculator.get_stored_wavefunction()
    # Only use alpha part
    # TODO: make this toggable?
    _, alpha_epos = epos_from_wf(wf, fragments, pop_kind=pop_kind)
    property = alpha_epos
    return property


def epos_property_iao(geom, fragments):
    return epos_property(geom, fragments, pop_kind=mdtypes.PopKind.IAO)


def epos_property_mulliken(geom, fragments):
    return epos_property(geom, fragments, pop_kind=mdtypes.PopKind.MULLIKEN)


def en_exc_property(geom, fragments):
    gs_energy, *es_energies = geom.all_energies
    assert len(es_energies) > 0
    property = es_energies[0] - gs_energy
    return property


def calculate_property(
    i, coords, geom, calc_getter, pal, fragments, prop_eq, property_func
):
    displ_geom = geom.copy()
    displ_geom.coords = coords
    displ_geom.set_calculator(calc_getter(pal=pal, calc_number=i, base_name="displ"))
    prop = np.nan
    try:
        # TODO: move actual calculation into property functions?!
        # displ_geom.all_energies
        # TODO: do a different calculation for DEXC
        displ_geom.calc_wavefunction()
        prop = property_func(displ_geom, fragments)
        # Subtract property at equilibrium geometry
        # If we would not subtract the equilibrium property we would have to pad
        # the normal coordinate matrix with an additional column.
        prop = prop - prop_eq
    except CalculationFailedException as err:
        print(err)
        print(f"Calculation {i:03d} failed!")
    except RunAfterCalculationFailedException as err:
        print(err)
        print(f"Postprocessing of calculation {i:03d} failed!")
    return prop


def get_mass(marcus_dim: np.ndarray, masses: np.ndarray) -> float:
    """Mass of fictious particle moving along Marcus dimension.

    Eq. (2) in SI of [2].

    Parameters
    ----------
    marcus_dim
        Array of shape (3N, ) or (N, 3) with N being the number of atoms.
    masses
        Array of shape (N, ) containing the atom masses in atomic mass units.

    Returns
    -------
    mass
        Mass of fictious particle moving along Marcus dimension.
    """
    marcus_dim3d = marcus_dim.reshape(-1, 3)
    natoms, _ = marcus_dim3d.shape
    assert natoms == len(masses)
    return (marcus_dim3d**2 * masses[:, None]).sum()


def get_wavenumber(
    marcus_dim: np.ndarray, nus: np.ndarray, eigvecs: np.ndarray, masses: np.ndarray
) -> float:
    """Wavenumber of Marcus dimension.

    Eq. (3) in SI of [2].

    The angular frequency is obtained from the square root of the quotient between
    force constants and reduced mass.

        ω = sqrt(k/μ)

    The 'normal' frequency ν is obtained as

        ω_i = 2π ν_i = sqrt(k_i/μ_i)
        ν_i = 1/2π * sqrt(k_i/μ_i)
        ν_i = sqrt(k_i/(4π²μ_i))          (1)

    Force constant k_i is

        k_i = 4π²ν_i²μ_i                  (2)

    The general ideo of Eq. (3) in the SI of [2] seems to be to calculate the
    force constants of all normal modes in the numerator (Ω² C^T M C; compare to
    ν² μ in (2)), then contract these values with the Marcus dimension in normal
    coordinates. This number is then divided by the mass of the Marcus dimension.
    The square root of this fraction yields a frequency with unit 1/s!

    Parameters
    ----------
    marcus_dim
        Array of shape (3N, ) or (N, 3) with N being the number of atoms.
    nus
        Array of wavenumbers of shape (N, ).
    eigvecs
        Array of eigenvectors of the projected mass-weighted Hessian with shape (3N, 3N).
    masses
        Array of shape (N, ) containing the atom masses in atomic mass units.

    Returns
    -------
    nu_marcus
        Wavenumber along Marcus dimension in cm⁻¹.
    """
    masses_rep = np.repeat(masses, 3)
    a = eigvecs.T @ (np.sqrt(masses_rep) * marcus_dim)
    a = a / np.linalg.norm(a)

    nus_m = nus * 100  # Wavenumber in 1/m
    freqs_s = C * nus_m
    # Matrix multiplications with diagonal matrices are replaced by
    # element-wise multiplications.
    force_constant = (
        a[None, :]
        @ ((freqs_s**2)[:, None] * eigvecs.T)
        @ (masses_rep[:, None] * eigvecs)
        @ a
    )
    mass = get_mass(marcus_dim, masses)
    freq_marcus = np.sqrt(force_constant / mass)
    nu_marcus = float(freq_marcus / C / 100)
    return nu_marcus


def get_wavenumber_from_coeffs(coeffs: np.ndarray, nus: np.ndarray) -> float:
    """Wavenumber in cm⁻¹ of Marcus dimension from fit coefficients."""
    coeffs2 = coeffs**2
    coeff_norm = (coeffs2).sum()
    marcus_nu = (coeffs2 * nus).sum() / coeff_norm
    return marcus_nu


def dump_marcus_dim_xyz_trj(
    atoms, coords3d_eq: np.ndarray, marcus_dim: np.ndarray, out_path: Path
):
    # XYZ representation of Marcus dimension and animation
    marcus_dim_xyz = make_xyz_str(atoms, BOHR2ANG * marcus_dim.reshape(-1, 3))
    xyz_path = out_path.with_suffix(".xyz")
    with open(xyz_path, "w") as handle:
        handle.write(marcus_dim_xyz)

    print(f"Current Marcus dimension:\n{marcus_dim_xyz}\n")
    marcus_dim_trj = get_tangent_trj_str(
        atoms, coords3d_eq.flatten(), marcus_dim, comment="Marcus dimension"
    )
    # Dump animated Marcus dimension
    trj_path = out_path.with_suffix(".trj")
    with open(trj_path, "w") as handle:
        handle.write(marcus_dim_trj)
    return marcus_dim_xyz, marcus_dim_trj


def setup_calc_client(geom, scheduler):
    if scheduler is not None:
        client = Client(scheduler)
        pal = 1
        sched_info = client.scheduler_info()
        n_workers = len(sched_info["workers"])
        calc_msg = f"Running {n_workers} calculations in parallel with {pal=} each."
    else:
        client = None
        pal = geom.calculator.pal
        calc_msg = f"Running calculations in serial with {pal=}."
    return client, pal, calc_msg


def run_batch_calcs(batch_size, start_ind, sampler, client, calc_func):
    # Get and store displaced coordinates and normal coordinates
    batch_displ_coords = list()
    batch_norm_coords = list()
    for i in range(batch_size):
        displ_coords, norm_coords = sampler()
        batch_displ_coords.append(displ_coords)
        batch_norm_coords.append(norm_coords)

    # Loop over cacluclations in current batch
    batch_range = range(start_ind, start_ind + batch_size)
    if client is not None:
        futures = client.map(calc_func, batch_range, batch_displ_coords, pure=False)
        batch_properties = client.gather(futures)
    else:
        batch_properties = list()
        for i, displ_coords in enumerate(batch_displ_coords, start_ind):
            if i % 5 == 0:
                print(f"\t{i}")
                sys.stdout.flush()
            prop = calc_func(i, displ_coords)
            batch_properties.append(prop)
    return batch_displ_coords, batch_norm_coords, batch_properties


def fit_marcus_dim(
    geom: Geometry,
    calc_getter: Callable,
    fragments: List[Tuple[int]],
    T: float,
    property=mdtypes.Property.EPOS_IAO,
    batch_size: int = 25,
    max_batches: int = 20,
    rms_thresh: float = RMS_THRESH,
    correlations: bool = False,
    corr_thresh: float = _CORR_THRESH,
    scheduler=None,
    seed: Optional[int] = None,
    out_dir=".",
):
    assert T >= 0.0
    assert batch_size > 0
    assert max_batches > 0
    assert rms_thresh > 0.0
    max_ncalcs = batch_size * max_batches

    property_funcs = {
        property.EPOS_IAO: epos_property_iao,
        property.EPOS_MULLIKEN: epos_property_mulliken,
        property.EEXC: en_exc_property,
    }
    property_func = property_funcs[property]

    out_dir = Path(out_dir)

    charge = geom.calculator.charge
    mult = geom.calculator.mult

    if (charge is not None) and (mult is not None):
        print(f"       Charge: {charge}")
        print(f" Multiplicity: {mult}")
    print(f"  Temperature: {T:.2f} K")
    print(f"   Batch size: {batch_size}")
    print(f"  Max batches: {max_batches}")
    print(f"rms threshold: {rms_thresh}")
    print(f"     Property: {property}")

    print(f"Doing at most {batch_size * max_batches} calculations")

    # Dump fragments
    fragments_trj = "\n".join(
        get_fragment_xyzs(geom, fragments, with_geom=True, with_dummies=True)
    )
    fragment_trj_fn = "fragments.trj"
    with open(fragment_trj_fn, "w") as handle:
        handle.write(fragments_trj)
    print(f"Wrote fragments to '{fragment_trj_fn}'")

    coords_eq = geom.cart_coords
    nus, eigvals, eigvecs, cart_displs = geom.get_normal_modes()
    nmodes = len(nus)
    drop_first = 3 * len(geom.atoms) - nmodes

    # Calculate wavefunction at equilibrium geometry
    print("Starting calculation at equilibrium geometry")
    # TODO: figure out why this is even needed here ... probably because the calculations aren't done
    # inside the property function(s).
    gs_energy_eq, *_ = geom.all_energies
    # TODO: using 'geom.calc_wavefunction()' leads to a pickling-error later on
    prop_eq = property_func(geom, fragments)
    print("Finished calculation at equilibrium geometry")

    # Function that create displaced geometries by drawing from a Wigner distribution
    wigner_sampler, seed = get_wigner_sampler(geom, temperature=T, seed=seed)
    sampler_wrapped = functools.partial(
        get_displaced_coordinates,
        wigner_sampler=wigner_sampler,
        cart_displs=cart_displs,
        coords_eq=coords_eq,
    )
    print(f"Seed for Wigner sampling: {seed}")

    # Arrays holding displace Cartesian coordinates, normal coordinates and properties
    all_displ_coords = np.zeros((max_ncalcs, coords_eq.size))
    all_norm_coords = np.zeros((max_ncalcs, nmodes))
    all_properties = np.zeros(max_ncalcs)

    masses = geom.masses

    to_save = {
        # Property key
        "property": str(property),
        # Equilibrium geometry
        "cart_coords_eq": geom.cart_coords,
        "property_eq": prop_eq,
        "gs_energy_eq": gs_energy_eq,
        # Samples
        "hessian": geom.cart_hessian,
        "mw_hessian": geom.mw_hessian,
        "linear": geom.is_linear,
        "wigner_seed": seed,
        "cart_coords": all_displ_coords,
        "normal_coordinates": all_norm_coords,
        "properties": all_properties,
        "masses": masses,
        "eigvals": eigvals,
        "eigvecs": eigvecs,
        "cart_displs": cart_displs,
        "nus": nus,
        # Additional keys will be added/updated throughout the run.
        # "marcus_dim": np.zeros_like(geom.cart_coords),
        # "coeffs": np.zeros_like(all_norm_coords),
        "rms_thresh": rms_thresh,
        "max_batches": max_batches,
    }
    results_fn = out_dir / FIT_RESULTS_FN

    client, pal, calc_msg = setup_calc_client(geom, scheduler)
    print(calc_msg)

    calculate_property_wrapped = functools.partial(
        calculate_property,
        geom=geom,
        calc_getter=calc_getter,
        pal=pal,
        fragments=fragments,
        prop_eq=prop_eq,
        property_func=property_func,
    )
    print()
    sys.stdout.flush()

    # Start calculation batches
    prev_coeffs = None
    prev_marcus_dim = None
    rms_coeffs = None
    rms_marcus_dim = None
    for batch in range(max_batches):
        batch_str = f"batch_{batch:02d}"
        start_ind = batch * batch_size
        end_ind = start_ind + batch_size
        now = datetime.datetime.now()
        print(highlight_text(f"Batch {batch}") + "\n")
        print(f"Starting at index {start_ind} on {now.strftime('%c')}")
        print(f"Doing calculations with indices {start_ind} to {end_ind - 1}.")
        sys.stdout.flush()

        # Run the actual calculatiions
        batch_dur = time.time()
        batch_displ_coords, batch_norm_coords, batch_properties = run_batch_calcs(
            batch_size,
            start_ind,
            sampler_wrapped,
            client,
            calculate_property_wrapped,
        )
        batch_dur = time.time() - batch_dur
        all_displ_coords[start_ind:end_ind] = batch_displ_coords
        all_norm_coords[start_ind:end_ind] = batch_norm_coords
        all_properties[start_ind:end_ind] = batch_properties

        calc_dur = batch_dur / batch_size
        print(f"Did {batch_size} calculations.")
        print(f"Full batch took {batch_dur:.2f} s, {calc_dur:.2f} s / calculation")
        print(f"Total number of calculations done until now: {end_ind}")
        sys.stdout.flush()

        # Actually calculate Marcus dimension using least-squares
        corrs, coeffs, marcus_dim, marcus_dim_q = get_marcus_dim(
            cart_displs,
            all_norm_coords[:end_ind],
            all_properties[:end_ind],
        )
        # Report expansion coefficients and rms values
        report_marcus_coefficients(coeffs, nus, drop_first)
        # Dump Marcus dimension and .trj file
        marcus_dim_trj_fn = out_dir / f"marcus_dim_{batch_str}.trj"
        # Suppress the 2nd argument (marcus_dim_trj)
        marcus_dim_xyz, _ = dump_marcus_dim_xyz_trj(
            geom.atoms, geom.coords3d, marcus_dim, marcus_dim_trj_fn
        )

        # Property change along Marcus dimension from fitted coefficients
        prop_change_along_marcus_dim = marcus_dim_q.dot(coeffs)

        # Mass and wavenumber of Marcus dimension mode
        mass_marcus = get_mass(marcus_dim, masses)
        # nu_marcus = get_wavenumber(marcus_dim, nus, eigvecs, masses)
        nu_marcus = get_wavenumber_from_coeffs(coeffs, nus)
        print(f"Mass along Marcus dimension: {mass_marcus:.6f} amu")
        print(f"Wavenumber of Marcus dimension: {nu_marcus:.6f} cm⁻¹")

        # Dump results
        to_save["coeffs"] = coeffs
        to_save["marcus_dim"] = marcus_dim
        to_save["marcus_dim_q"] = marcus_dim_q
        to_save["mass_marcus"] = mass_marcus
        to_save["nu_marcus"] = nu_marcus
        to_save["prop_change_along_marcus_dim"] = prop_change_along_marcus_dim
        to_save["end_ind"] = end_ind

        sys.stdout.flush()
        print()
        if prev_coeffs is not None:
            rms_coeffs = rms(prev_coeffs - coeffs)
            rms_marcus_dim = rms(prev_marcus_dim - marcus_dim)
            print(
                f"rms(Δcoeffs)={rms_coeffs:.6f}, rms(ΔMarcus dim)={rms_marcus_dim:.6f}"
            )
            rms_converged = rms_marcus_dim <= rms_thresh
            # Will result in scalar entries
            to_save["rms_coeffs"] = rms_coeffs
            to_save["rms_marcus_dim"] = rms_marcus_dim
            to_save["rms_converged"] = rms_converged
        else:
            rms_converged = False
        # Save the actual data here
        np.savez(results_fn, **to_save)

        # Plot of expansion coefficients
        fig, _ = create_marcus_coefficient_plot(
            nmodes, coeffs, batch=batch, ncalcs=end_ind
        )
        plot_fn = out_dir / f"marcus_coeffs_{batch_str}.svg"
        fig.savefig(plot_fn)
        print(f"Saved expansion coefficent plot to '{plot_fn}'.\n")
        plt.close(fig)

        # Plot summary of utilized normal coordinates
        fig, _ = plot_normal_coords(all_norm_coords[:end_ind])
        fig.savefig(out_dir / "normal_coords.svg", dpi=200)
        plt.close(fig)

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
        sys.stdout.flush()

        # Convergence check
        if (rms_marcus_dim is not None) and rms_converged:
            print(
                f"\nCalculation of Marcus dimension converged after {end_ind} calculations!"
            )
            break

        sign = check_for_end_sign()
        if sign in ("stop", "converged"):
            print(f"Found '{sign}' sign. Breaking.")
            break
        prev_coeffs = coeffs
        prev_marcus_dim = marcus_dim
        print()
    else:
        print("Reached maximum number of requested batches!")
    # End of loop

    # Write final Marcus dimension
    with open(out_dir / "marcus_dim.xyz", "w") as handle:
        handle.write(marcus_dim_xyz)

    # Report failed calculations
    calc_indices = np.arange(end_ind)
    failed_mask = np.isnan(all_properties[:end_ind])
    for i in calc_indices[failed_mask]:
        print(f"Calculation {i:03d} failed!")
    print(f"{failed_mask.sum()}/{end_ind} calculations failed.")
    print()
    print(f"Final Marcus dimension:\n{marcus_dim_xyz}\n")
    print(f"Mass along final Marcus dimension: {mass_marcus:.6f} amu")
    print(f"Wavenumber of final Marcus dimension: {nu_marcus:.6f} cm⁻¹")

    return to_save
