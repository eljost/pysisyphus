import math
from pathlib import Path
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np


from pysisyphus.constants import AU2EV
from pysisyphus.marcus_dim.config import SCAN_RESULTS_FN


def scan_dir(
    x0,
    direction,
    get_property,
    step_size=0.05,
    add_steps=10,
    max_steps=500,
    min_steps=5,
    grad_thresh=1e-2,
):
    assert step_size > 0.0, f"{step_size=} must be positive!"
    assert add_steps >= 0, f"{add_steps=} must must be positive!"
    assert max_steps > 0, f"{max_steps=} must be positive!"
    assert min_steps >= 0, f"{min_steps=} must be positive!"
    assert grad_thresh > 0.0, f"{grad_thresh} must be positive!"

    step = step_size * direction
    stop_in = add_steps
    xcur = x0 + step

    converged = False
    grad = np.nan
    prop_prev = None
    abs_grad_prev = None
    grad_decreased_already = False

    all_factors = np.arange(max_steps) + 1
    all_coords = np.empty((max_steps, *x0.shape))
    all_energies = np.empty((max_steps, 2))
    all_properties = np.empty(max_steps)
    for i in range(max_steps):
        factor = all_factors[i]
        all_coords[i] = xcur
        # Calculate & store property
        energies, prop = get_property(factor, xcur)
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

        # Break when the gradient already decreased once and increased
        # unexpectedly aftwards. But do at least 'min_steps' steps.
        if (
            (i >= min_steps)
            and grad_decreased_already
            and (grad_decreased is not None)
            and not grad_decreased
        ):
            print("Unexpected increase of abs(grad(property))! Breaking")
            break

        # Check gradient convergence; this check is skipped once convergence is indicated.
        if not converged and (converged := np.abs(grad) <= grad_thresh):
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


def scan(coords_init, direction, get_properties, out_dir=".", **kwargs):
    out_dir = Path(out_dir)
    dir_norm = np.linalg.norm(direction)
    if not math.isclose(dir_norm, 1.0):
        warnings.warn(f"norm(direction)={dir_norm:.6f} is not 1.0! Renormalizing.")
        direction = direction / dir_norm
    # Carry out calculation on initial geometry.
    ens0, prop0 = get_properties(0.0, coords_init)

    def get_property_changes(factor, xi):
        """Get property changes w.r.t. initial geometry.

        TODO: rename this because prop0 isn't substracted ..."""
        ens, prop = get_properties(factor, xi)
        return ens, prop  # - prop0

    print("Positive direction")
    pos_dir = direction
    pos_facts, pos_coords, pos_ens, pos_props = scan_dir(
        coords_init, pos_dir, get_property_changes, **kwargs
    )
    print()

    print("Negative direction")
    neg_dir = -1.0 * direction
    neg_facts, neg_coords, neg_ens, neg_props = scan_dir(
        coords_init, neg_dir, get_property_changes, **kwargs
    )

    def concat(neg, init, pos):
        return np.concatenate((neg[::-1], [init], pos))

    all_facts = concat(-neg_facts, 0.0, pos_facts)
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
        "scan_converged": True,
    }
    scan_results_fn = out_dir / SCAN_RESULTS_FN
    np.savez(scan_results_fn, **to_save)

    return all_facts, all_coords, all_energies, all_properties


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
