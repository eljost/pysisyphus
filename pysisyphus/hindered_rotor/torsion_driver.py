import dataclasses
import functools
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

from pysisyphus.constants import AU2KJPERMOL, AMU2AU, AU2SEC, AU2NU, AU2EV
from pysisyphus.finite_diffs import periodic_fd_2_8
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.hindered_rotor import (
    fragment as hr_fragment,
    inertmom,
    opt as hr_opt,
    torsion_gpr,
)
from pysisyphus.partfuncs import partfuncs as pf
from pysisyphus.TablePrinter import TablePrinter


@dataclasses.dataclass
class TorsionGPRResult:
    atoms: tuple[str, ...]
    coords3d: np.ndarray
    inertmom_left: float
    fragment_left: list[int]
    inertmom_right: float
    fragment_right: list[int]
    grid: np.ndarray
    energies: np.ndarray
    std: np.ndarray
    grid_train: np.ndarray
    energies_train: np.ndarray
    eigvals: np.ndarray
    eigvecs: np.ndarray
    temperature: float
    # Boltzmann weights
    weights: np.ndarray
    # TODO: Store optimized geometries?!

    @property
    def dx(self):
        return abs(self.grid[1] - self.grid[0])

    def calc_hr_partfunc(self, temperature: Optional[float] = None):
        if temperature is None:
            temperature = self.temperature
        eigvals = self.eigvals - self.eigvals.min()
        return pf.sos_partfunc(eigvals, temperature)

    def calc_cancel_partfunc(self, temperature: Optional[float] = None):
        if temperature is None:
            temperature = self.temperature
        # TODO: make calculation of force constant more robust; currently
        # index 0 is hardcoded, but this is not necessarily the minimum.
        # If an index != 0 is supposed to be used we would also have to update
        # the stored moments of inertia, as they are currently calculated for
        # the initial geometry (with index 0) only.
        #
        # Idea: add index as argument; if not provided use actual minimum and
        # recalculate moments of inertia?
        energies_cut = self.energies[:-1]
        force_constant = periodic_fd_2_8(0, energies_cut, self.dx)
        ho_freq = np.sqrt(force_constant / self.inertmom_left) / (2 * np.pi)
        ho_freq_si = ho_freq / AU2SEC
        return pf.harmonic_quantum_partfunc(ho_freq_si, temperature)

    def calc_corr_factor(self, temperature: Optional[float] = None):
        if temperature is None:
            temperature = self.temperature
        hr_partfunc = self.calc_hr_partfunc(temperature)
        cancel_partfunc = self.calc_cancel_partfunc(temperature)
        return hr_partfunc / cancel_partfunc

    @property
    def B(self):
        return 1 / (2 * self.inertmom_left)

    @property
    def BeV(self):
        return self.B * AU2EV

    def __post_init__(self):
        # Calculation of partition function correction factor (HR / HO)
        #
        # Currently, we evaluate the partition function at the initial geometry.
        # TODO: make evaluation more flexible; Maybe pick the global minimum of the scan
        # and don't always stay at the initial geometry.
        # TODO: move partfunc stuff into separate function?!
        self.hr_partfunc = self.calc_hr_partfunc()
        self.cancel_partfunc = self.calc_cancel_partfunc()
        self.corr_factor = self.calc_corr_factor()

    def report(self, boltzmann_thresh: float = 0.9999):
        print(highlight_text("Hindered Rotor Scan Summary"))
        print(f"Rotor moment of inertia: {self.inertmom_left: >14.8e} au")
        print(f"           B = ħ²/(2·I): {self.BeV*1e3:>14.8e} meV")
        print(f"1d hindered rotor partition function: {self.hr_partfunc: >14.6f}")
        print(f"        HO cancel partition function: {self.cancel_partfunc: >14.6f}")
        print(f"           Correction factor (HR/HO): {self.corr_factor: >14.6f}")
        print(f"Reporting states until Σp_Boltz. >= {boltzmann_thresh}")
        print(f"                         Temperature: {self.temperature: >14.4f} K")
        header = ("# State", "E/Eh", "E/kJ mol⁻¹", "ΔE/cm⁻¹", "p_Boltz.", "Σp_Boltz.")
        col_fmts = ("{:3d}", "float", "{: 12.6f}", "{: >8.2f}", "float", "float")
        table = TablePrinter(header, col_fmts, width=12, sub_underline=False)

        weights_cum = np.cumsum(self.weights)
        nstates = np.argmax(weights_cum > boltzmann_thresh) + 1

        eigvals = self.eigvals - self.eigvals.min()
        eigvals_kJ = eigvals * AU2KJPERMOL
        eigvals_nu = eigvals * AU2NU
        nrows = len(eigvals)

        table.print_header()
        table.print_rows(
            (range(nrows), eigvals, eigvals_kJ, eigvals_nu, self.weights, weights_cum),
            first_n=nstates,
        )


def run(
    geom: Geometry,
    indices: list[int],
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
    assert (calc_getter is not None) or (
        energy_getter is not None
    ), "Either 'calc_getter' or 'energy_getter' must be provided!"
    if opt_kwargs is None:
        opt_kwargs = {}
    out_dir = Path(out_dir)
    if not out_dir.exists():
        out_dir.mkdir()

    bond = indices[1:3]

    fragment_left, fragment_right = hr_fragment.fragment_geom(geom, bond)
    imom_left, imom_right = inertmom.get_top_moment_of_inertia(
        geom.coords3d, geom.masses, fragment_left, bond, m=2, n=3
    )
    imom_left *= AMU2AU
    imom_right *= AMU2AU
    mass = imom_left
    if calc_getter is not None:
        energy_getter = hr_opt.opt_closure(
            geom,
            indices,
            calc_getter,
            opt_kwargs,
            single_point_calc_getter=single_point_calc_getter,
            out_dir=out_dir,
        )
    rad_store = energy_getter.rad_store

    # The step is later required for finite differences
    grid, dx = np.linspace(
        torsion_gpr.PERIOD_LOW, torsion_gpr.PERIOD_HIGH, num=npoints, retstep=True
    )
    part_callback = functools.partial(callback, plot=plot, out_dir=out_dir)

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
    for ind in np.argsort(rads):
        key = rads[ind]
        en, sp_en, c3d = rad_store[key]
        comment = f"{en:{en_fmt}}"
        if sp_en is not np.nan:
            comment += f", sp_energy={sp_en:{en_fmt}}"
        xyz = geom.as_xyz(cart_coords=c3d, comment=comment)
        xyzs.append(xyz)
    trj = "\n".join(xyzs)
    with open(out_dir / "torsion_scan.trj", "w") as handle:
        handle.write(trj)

    result = TorsionGPRResult(
        geom.atoms,
        geom.coords3d,
        grid=grid,
        energies=gpr_status.energies,
        std=gpr_status.std,
        inertmom_left=imom_left,
        fragment_left=fragment_left,
        inertmom_right=imom_right,
        fragment_right=fragment_right,
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
