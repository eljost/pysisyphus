import math
from pathlib import Path
import sys
import warnings

from distributed import Client
import matplotlib.pyplot as plt
import numpy as np


from pysisyphus.constants import AU2EV
from pysisyphus.exceptions import (
    CalculationFailedException,
    RunAfterCalculationFailedException,
)
from pysisyphus.marcus_dim.config import (
    MAX_BOND_CHANGE,
    SCAN_RESULTS_FN,
    SCAN_STEP_LENGTH,
)
from pysisyphus.marcus_dim.fit import get_marcus_dim
import pysisyphus.marcus_dim.types as mdtypes


def scan_coords_for_max_bond_change(
    B: np.ndarray,
    cart_displs: np.ndarray,
    max_bond_change: float = MAX_BOND_CHANGE,
    step_length: float = 0.01,
):
    assert max_bond_change > 0.0
    assert 0.0 < step_length <= max_bond_change
    cart_dipls = cart_displs.copy()
    cart_displs /= np.linalg.norm(cart_displs)
    cart_displs = cart_displs.flatten()
    step = cart_displs * step_length
    steps = list()
    cur_step = np.zeros_like(step)
    while True:
        cur_step += step
        cur_int_step = B @ cur_step
        max_int_step = np.abs(cur_int_step).max()
        steps.append(cur_step.copy())
        if max_int_step >= max_bond_change:
            break
    steps = np.array(steps)
    factors = np.arange(len(steps), dtype=float) + 1.0
    factors *= step_length
    return factors, steps


def get_classII_scan_info(
    marcus_dim: np.ndarray,
    normal_coords: np.ndarray,
    properties: np.ndarray,
    eigvecs: np.ndarray,
    masses: np.ndarray,
    max_bond_change: float = MAX_BOND_CHANGE,
):
    """Determe directions & lengths for Robin-Day class II systems.

    In class II systems we start at one of two minima and the scan along the Marcus
    dimension will be asymmetric. One direction of the scan will pass a barrier,
    towards the second minimum. This (sub)-scan will be longer. The scan in the opposite
    direction can be shorter and no big changes regarding electron position are expected.
    """
    assert max_bond_change > 0.0
    _, coeffs, marcus_dim, marcus_dim_q = get_marcus_dim(
        eigvecs, masses, normal_coords, properties
    )

    prop_change_along_marcus_dim = (marcus_dim_q * coeffs).sum()
    # When prop_change_along_marcus_dim is positive, then a positive displacement along
    # the Marcus dimension will increase the excitation energy and move away from the barrier.
    # If the change is negative, then displacing against the Marcus dimension will bring
    # us towards the barrier.
    to_lower_eexc_factor = -1 if prop_change_along_marcus_dim > 0 else 1

    """
    # Alternatively one could also project the Marcus dimension in normal coordinates
    # onto the calculation that yielded the biggest property change.
    max_change_ind = np.abs(props).argmax()
    max_change = props[max_change_ind]
    max_change_q = normal_coords[max_change_ind]
    ovlp = (
        max_change_q.dot(marcus_dim_q)
        / np.linalg.norm(max_change_q)
        / np.linalg.norm(marcus_dim_q)
    )
    """

    # In the direction towards lower excitation energies we have to pass the barrier
    # in the scan, so we scan longer in this direction.
    marcus_dim_to_lower_eexc = marcus_dim * to_lower_eexc_factor
    marcus_dim_to_higher_eexc = -marcus_dim_to_lower_eexc
    max_change_to_lower_eexc = max_bond_change
    max_change_to_higher_eexc = max_bond_change / 2
    return (
        marcus_dim_to_lower_eexc,
        max_change_to_lower_eexc,
        marcus_dim_to_higher_eexc,
        max_change_to_higher_eexc,
    )


def get_scan_factors_and_steps(
    B,
    normal_coords,
    properties,
    eigvecs,
    masses,
    marcus_dim,
    rd_class: mdtypes.RobinDay,
    max_bond_change: float = MAX_BOND_CHANGE,
    step_length: float = SCAN_STEP_LENGTH,
):
    # TODO: switch to ENUM and also derive property from RDCLASS
    if rd_class == mdtypes.RobinDay.CLASS2:
        # Class II
        (
            md_lower,
            max_change_lower,
            md_higher,
            max_change_higher,
        ) = get_classII_scan_info(
            marcus_dim,
            normal_coords,
            properties,
            eigvecs,
            masses,
            max_bond_change=max_bond_change,
        )

        scan_factors_lower, scan_steps_lower = scan_coords_for_max_bond_change(
            B,
            md_lower,
            max_bond_change=max_change_lower,
            step_length=step_length,
        )
        scan_factors_higher, scan_steps_higher = scan_coords_for_max_bond_change(
            B,
            md_higher,
            max_bond_change=max_change_higher,
            step_length=step_length,
        )
        scan_factors_left = scan_factors_lower
        scan_steps_left = scan_steps_lower
        scan_factors_right = scan_factors_higher
        scan_steps_right = scan_steps_higher
    elif rd_class == mdtypes.RobinDay.CLASS3:
        max_change = max_bond_change / 2.0
        scan_factors, scan_steps = scan_coords_for_max_bond_change(
            B,
            marcus_dim,
            max_bond_change=max_change,
            step_length=step_length,
        )
        scan_factors_left = scan_factors
        scan_steps_left = scan_steps
        scan_factors_right = -scan_factors
        scan_steps_right = -scan_steps
    else:
        raise Exception(f"Unknown Robin-Day class: '{rd_class}'!")

    scan_factors = np.concatenate(
        (
            scan_factors_left[::-1],
            [0.0],
            scan_factors_right,
        )
    )
    zero_step = np.zeros_like(marcus_dim)
    scan_steps = np.concatenate(
        (scan_steps_left[::-1], zero_step[None, :], scan_steps_right), axis=0
    )
    return scan_factors, scan_steps


def scan_dir(
    x0,
    direction,
    prefact,
    get_property,
    step_size=0.05,
    add_steps=10,
    max_steps=500,
    min_steps=10,
    grad_thresh=1e-2,
):
    assert step_size > 0.0, f"{step_size=} must be positive!"
    assert add_steps >= 0, f"{add_steps=} must must be positive!"
    assert max_steps > 0, f"{max_steps=} must be positive!"
    assert min_steps >= 0, f"{min_steps=} must be positive!"
    assert grad_thresh > 0.0, f"{grad_thresh} must be positive!"

    step = prefact * step_size * direction
    stop_in = add_steps
    xcur = x0 + step

    converged = False
    grad = np.nan
    prop_prev = None
    abs_grad_prev = None
    grad_decreased_already = False

    all_factors = prefact * (np.arange(max_steps) + 1)
    all_coords = np.empty((max_steps, *x0.shape))
    all_energies = np.empty((max_steps, 2))
    all_properties = np.empty(max_steps)
    for i in range(max_steps):
        factor = all_factors[i]
        all_coords[i] = xcur
        # Calculate & store property
        try:
            energies, prop = get_property(factor, xcur)
        except (CalculationFailedException, RunAfterCalculationFailedException):
            print("Calculation failed! Breaking.")
            break
        all_properties[i] = prop
        all_energies[i] = energies

        # Determine gradient from finite differences
        if prop_prev is not None:
            grad = (prop - prop_prev) / step_size
            abs_grad = abs(grad)
        else:
            abs_grad = None

        if abs_grad_prev is not None:
            grad_decreased = abs_grad < abs_grad_prev
            grad_decreased_already = grad_decreased_already or grad_decreased
        else:
            grad_decreased = None
            grad_decreased_already = False
        print(f"{i=:03d}, property={prop: >12.4f}, {grad=: >12.4f}")
        sys.stdout.flush()

        did_min_steps = i >= min_steps - 1
        # Break when the gradient already decreased once and increased
        # unexpectedly aftwards. But do at least 'min_steps' steps.
        if (
            did_min_steps
            and grad_decreased_already
            and (grad_decreased is not None)
            and not grad_decreased
        ):
            print("Unexpected increase of abs(grad(property))! Breaking")
            break

        # Check gradient convergence; this check is skipped once convergence is indicated.
        if not converged and (
            # Bad idea to force to continue when converged?!
            converged := did_min_steps
            and np.abs(grad) <= grad_thresh
        ):
            print("Converged!")

        # If requested, we carry out additional steps, if requested.
        if add_steps and converged:
            stop_in -= 1

        # Break directly if converged and we don't want to do any additional steps.
        if converged and add_steps == 0:
            break
        elif add_steps and stop_in < 0:
            print("Did additional steps")
            break

        # Update variables
        xcur = xcur + step
        prop_prev = prop
        if abs_grad is not None:
            abs_grad_prev = abs_grad
    else:
        print("Reached maximum number of cycles in scan.")

    all_energies = np.array(all_energies)
    # Truncate arrays and drop empty part. This will also drop the last calculation
    # that lead to the break from the loop.
    end_ind = i
    return (
        all_factors[:end_ind] * step_size,
        all_coords[:end_ind],
        all_energies[:end_ind],
        all_properties[:end_ind],
    )


def scan(
    coords_init,
    direction,
    get_properties,
    pos_min_steps=None,
    neg_min_steps=None,
    out_dir=".",
    **kwargs,
):
    out_dir = Path(out_dir)
    dir_norm = np.linalg.norm(direction)
    if not math.isclose(dir_norm, 1.0):
        warnings.warn(f"norm(direction)={dir_norm:.6f} is not 1.0! Renormalizing.")
        direction = direction / dir_norm
    # Carry out calculation on initial geometry.
    ens0, prop0 = get_properties(0.0, coords_init)

    print("Positive direction")
    pos_kwargs = kwargs.copy()
    if pos_min_steps is not None:
        assert pos_min_steps >= 0
        pos_kwargs["min_steps"] = pos_min_steps
    pos_facts, pos_coords, pos_ens, pos_props = scan_dir(
        coords_init, direction, 1.0, get_properties, **pos_kwargs
    )
    print()

    print("Negative direction")
    neg_kwargs = kwargs.copy()
    if neg_min_steps is not None:
        assert neg_min_steps >= 0
        neg_kwargs["min_steps"] = neg_min_steps
    neg_facts, neg_coords, neg_ens, neg_props = scan_dir(
        coords_init, direction, -1.0, get_properties, **neg_kwargs
    )

    def concat(neg, init, pos):
        return np.concatenate((neg[::-1], [init], pos))

    all_facts = concat(neg_facts, 0, pos_facts)
    all_coords = concat(neg_coords, coords_init, pos_coords)
    all_energies = concat(neg_ens, ens0, pos_ens)
    # When we consider the difference w.r.t. initial geometry then
    # the property is always 0.0.
    all_properties = concat(neg_props, prop0, pos_props)

    to_save = {
        "factors": all_facts,
        "coords": all_coords,
        "energies": all_energies,
        "properties": all_properties,
        "scan_done": True,
    }
    scan_results_fn = out_dir / SCAN_RESULTS_FN
    np.savez(scan_results_fn, **to_save)

    return mdtypes.MarcusDimScanResult(
        factors=all_facts,
        coords=all_coords,
        energies=all_energies,
        properties=all_properties,
    )


def scan_parallel(
    coords_init,
    get_properties,
    marcus_dim,
    normal_coords,
    properties,
    eigvecs,
    masses,
    B,
    rd_class,
    max_bond_change: float = MAX_BOND_CHANGE,
    step_length: float = SCAN_STEP_LENGTH,
    scheduler=None,
    out_dir=".",
):
    out_dir = Path(out_dir)
    marcus_dim_norm = np.linalg.norm(marcus_dim)
    if not math.isclose(marcus_dim_norm, 1.0):
        warnings.warn(
            f"norm(direction)={marcus_dim_norm:.6f} is not 1.0! Renormalizing."
        )
        marcus_dim = marcus_dim / marcus_dim_norm
    # Carry out calculation on initial geometry.
    scan_factors, scan_steps = get_scan_factors_and_steps(
        B,
        normal_coords,
        properties,
        eigvecs,
        masses,
        marcus_dim,
        rd_class,
        max_bond_change=max_bond_change,
        step_length=step_length,
    )
    scan_coords = coords_init[None, :] + scan_steps
    ncalcs = scan_factors.size
    scan_calc_numbers = np.arange(ncalcs)
    print(f"There are {ncalcs} calculations in the scan.")
    print(f"Max bond length change: {max_bond_change:>8.4f} au")
    print(f"           Step length: {step_length:>8.4f} au")

    client = None
    if scheduler is not None:
        client = Client(scheduler)

    if client is not None:
        # Distribute calculations
        futures = client.map(
            get_properties,
            # scan_factors,
            scan_calc_numbers,
            scan_coords,
            pure=False,
        )
        scan_energies, scan_properties = zip(*client.gather(futures))
    else:
        scan_energies = list()
        scan_properties = list()
        # for factor, coords in zip(scan_factors, scan_coords):
        # ens, props = get_properties(factor, coords)
        for calc_number, coords in zip(scan_calc_numbers, scan_coords):
            ens, props = get_properties(calc_number, coords)
            scan_energies.append(ens)
            scan_properties.append(props)

    # Drop Nans from the arrays
    def isnan(obj):
        try:
            result = math.isnan(obj)
        except TypeError:
            result = False
        return result

    mask = [i for i, en in enumerate(scan_energies) if not isnan(en)]
    scan_energies = np.array([scan_energies[i] for i in mask])
    scan_properties = np.array([scan_properties[i] for i in mask])

    scan_factors = scan_factors[mask]
    scan_coords = scan_coords[mask]

    to_save = {
        "factors": scan_factors,
        "coords": scan_coords,
        "energies": scan_energies,
        "properties": scan_properties,
        "scan_done": True,
    }
    scan_results_fn = out_dir / SCAN_RESULTS_FN
    np.savez(scan_results_fn, **to_save)

    return mdtypes.MarcusDimScanResult(
        factors=scan_factors,
        coords=scan_coords,
        energies=scan_energies,
        properties=scan_properties,
    )


def plot_scan(factors, energies, properties, dummy_scan=False):
    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
    xlim = factors[[0, -1]]
    energies = energies - energies.min()
    energies *= AU2EV
    ax0.plot(factors, energies, "o-")
    ax0.set_ylabel("dE / eV")
    ax1.plot(factors, properties, "o-", label="props")
    ax1.set_ylabel("e$^-$ position")
    ax1.legend()
    ax1.set_xlabel("Marcus dimension / $a_0$")
    for ax in (ax0, ax1):
        ax.axvline(0.0, c="k", ls="--")
        ax.set_xlim(xlim)
    title = "Scan along Marcus dimension"
    if dummy_scan:
        title += " (dummy values!)"
    fig.suptitle(title)
    fig.tight_layout()
    return fig, (ax0, ax1)
