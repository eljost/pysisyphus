import functools
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import align_coords
from pysisyphus.hindered_rotor import (
    opt as hr_opt,
    torsion_gpr,
    TorsionGPRResult,
    RotorInfo,
)
from pysisyphus import xyzloader


def run(
    geom: Geometry,
    rotor_info: RotorInfo,
    calc_getter=None,
    single_point_calc_getter=None,
    energy_getter=None,
    opt_kwargs: Optional[dict] = None,
    npoints: int = 721,
    max_cycles: int = 75,
    partfunc_thresh: float = 1e-4,
    std_thresh: float = 1e-4,
    temperature: float = 298.15,
    plot: bool = True,
    out_dir: str | Path = ".",
):
    # Logical XOR; we need one of them but not both.
    assert bool(calc_getter) != bool(
        energy_getter
    ), "Either 'calc_getter' or 'energy_getter' must be provided but not both!"
    if opt_kwargs is None:
        opt_kwargs = {}
    out_dir = Path(out_dir)
    if not out_dir.exists():
        out_dir.mkdir()

    print(rotor_info.render_report())

    # When calc_getter is given we probably have to do actual constrained optimizations,
    # so we prepare the energy_getter.
    if calc_getter is not None:
        energy_getter = hr_opt.opt_closure(
            geom,
            rotor_info,
            calc_getter,
            opt_kwargs,
            single_point_calc_getter=single_point_calc_getter,
            out_dir=out_dir,
        )
    rad_store = energy_getter.rad_store

    # The step is later required for finite differences
    grid = np.linspace(torsion_gpr.PERIOD_LOW, torsion_gpr.PERIOD_HIGH, num=npoints)
    part_callback = functools.partial(callback, plot=plot, out_dir=out_dir)

    mass = rotor_info.imom_left
    gpr_status = torsion_gpr.run_gpr(
        grid.reshape(-1, 1),
        energy_getter,
        mass,
        callback=part_callback,
        max_cycles=max_cycles,
        partfunc_thresh=partfunc_thresh,
        std_thresh=std_thresh,
        temperature=temperature,
    )
    rads = list(rad_store.keys())
    xyzs = list()
    en_fmt = " >20.8f"
    # Final alignment of coordinates onto the initial geometry at Δrad = 0.0
    coords3d_aligned = list()
    for ind in np.argsort(rads):
        key = rads[ind]
        *_, c3d = rad_store[key]
        coords3d_aligned.append(c3d)
    coords3d_aligned = align_coords(coords3d_aligned)

    for ind in np.argsort(rads):
        key = rads[ind]
        en, sp_en, _ = rad_store[key]
        comment = f"{en:{en_fmt}}"
        if sp_en is not np.nan:
            comment += f", sp_energy={sp_en:{en_fmt}}"
        xyz = xyzloader.make_xyz_str_au(
            rotor_info.atoms, coords3d_aligned[ind], comment=comment
        )
        xyzs.append(xyz)
    trj = "\n".join(xyzs)
    with open(out_dir / "torsion_scan.trj", "w") as handle:
        handle.write(trj)

    result = TorsionGPRResult.TorsionGPRResult(
        rotor_info=rotor_info,
        grid=grid,
        energies=gpr_status.energies,
        std=gpr_status.std,
        grid_train=gpr_status.x_train,
        energies_train=gpr_status.y_train,
        eigvals=gpr_status.eigvals_absolute,
        eigvecs=gpr_status.eigvecs,
        temperature=temperature,
        weights=gpr_status.weights,
    )
    # TODO:
    # - report scan results similar to relaxed scans
    # - add support for known symmetry ... modify periodicity?!
    # - report numbers from Calculator.run_call_counts
    result.report()
    return result


def plot_summary(gpr_status, boltzmann_thresh=0.95):
    grid = gpr_status.grid.flatten()

    energies = gpr_status.energies
    en_min = energies.min()
    energies = (energies - en_min) * AU2KJPERMOL
    y_train = (gpr_status.y_train - en_min) * AU2KJPERMOL
    std_pred = gpr_status.std * AU2KJPERMOL

    fig = plt.figure(layout="constrained", figsize=(10, 5))
    gs = GridSpec(3, 2, figure=fig)
    ax = fig.add_subplot(gs[:2, 0])
    ax_acq = fig.add_subplot(gs[2, 0])
    ax_numerov = fig.add_subplot(gs[:, 1])
    ax.plot(grid, energies, c="orange", label="Fit")
    ax.axhline(energies.max(), c="k", ls=":")
    ax.scatter(gpr_status.x_train, y_train, s=50, c="red", label="Samples", zorder=5)
    # Next point
    x_next = gpr_status.x_next
    next_lbl = f"next trial at x={x_next: >8.4f}"
    ax.scatter(
        x_next, energies[gpr_status.ind_next], s=50, marker="D", label="next", zorder=5
    )

    # Wavefunction and energy levels
    weights = gpr_status.weights
    weights_cum = np.cumsum(weights)
    nstates = np.argmax(weights_cum > boltzmann_thresh) + 1
    weight_tot = weights[:nstates].sum()
    eigvals = gpr_status.eigvals_absolute
    eigvals = (eigvals - en_min) * AU2KJPERMOL
    y_range = energies.max() - energies.min()
    scale = y_range / nstates / 2.0
    ax_numerov.plot(grid, energies, c="orange")
    for j, wj in enumerate(eigvals[:nstates]):
        ax_numerov.axhline(wj, c="k", alpha=0.5, ls="--")
        # ax_numerov.plot(grid[:-1].flatten(), wj + scale * eigvecs[:, j], alpha=0.5)
        ax_numerov.plot(
            # no abs because we have real wavefunctions
            grid[:-1],
            wj + scale * gpr_status.eigvecs[:, j] ** 2,
            alpha=0.5,
            c="k",
        )
    ax_numerov.barh(eigvals[:nstates], weights[:nstates], label="Boltzmann weight")
    ax_numerov.set_title(
        f"First {nstates} state(s) with Σp_Boltzmann={weight_tot:.2f} "
        f"at {gpr_status.temperature:.2f} K"
    )
    ax_numerov.legend(loc="upper right")

    for ax_ in (ax, ax_numerov):
        ax_.fill_between(
            grid,
            energies - std_pred,
            energies + std_pred,
            alpha=0.1,
            color="grey",
            label="± std",
        )
        if energies.max() > 0.0:
            ax_.set_ylim(0, energies.max() * 1.25)
        ax_.set_ylabel("ΔE / kJ mol⁻¹")

    ax.legend(
        loc="upper center",
        ncols=5,
        prop={
            "size": 6,
        },
    )
    ax.set_title(f"Macro cycle {gpr_status.cycle:03d}, {next_lbl}")

    acq_iter = (
        ("Fleck", gpr_status.acq_func, "D"),
        ("Std", std_pred, "P"),
    )

    # Acquisition function
    for lbl, acq_func, marker in acq_iter:
        amax = acq_func.argmax()
        acq_norm = acq_func / acq_func.max()
        ax_acq.set_yscale("log")
        (acq_line,) = ax_acq.plot(grid, acq_norm, label=lbl)
        color = acq_line.get_color()
        ax_acq.scatter(
            grid[amax],
            acq_norm[amax],
            s=50,
            marker=marker,
            c=color,
            label=f"{lbl} next",
        )
    ax_acq.set_title("log(normalized acquisition function)")
    ax_acq.set_xlabel("ΔTorsion / rad")
    ax_acq.legend(
        loc="lower center",
        ncols=2,
        prop={
            "size": 6,
        },
    )
    for ax_ in (ax, ax_numerov, ax_acq):
        ax_.set_xlabel("ΔTorsion / rad")
        ax_.set_xlim(grid[0], grid[-1])
    return fig


def callback(gpr_status, plot=False, out_dir=Path(".")):
    gpr_status.dump_potential(out_dir)

    if plot:
        fig = plot_summary(gpr_status)
        out_fn = f"cycle_{gpr_status.cycle:03d}.png"
        fig.savefig(out_dir / out_fn, dpi=200)
        plt.close()
