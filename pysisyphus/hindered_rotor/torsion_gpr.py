import dataclasses
from pathlib import Path

import numpy as np
import optuna
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process.kernels import ExpSineSquared, ConstantKernel
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.utils._testing import ignore_warnings

from pysisyphus.constants import KBAU, AU2KJPERMOL


# TODO: do this in a function and revert afterwards?!
optuna.logging.set_verbosity(optuna.logging.WARNING)


# Periodicity bounds
PERIOD_BOUNDS = (0.0, 2 * np.pi)
PERIOD_LOW, PERIOD_HIGH = PERIOD_BOUNDS
# Length scale bound
LS_BOUNDS = (0.1, 2 * np.pi)
# Constant value bound
CONST_BOUNDS = (1e-6, 1e3)


def get_gpr(length_scale, constant_value) -> GaussianProcessRegressor:
    expsine_kwargs = {
        "length_scale": length_scale,
        "length_scale_bounds": LS_BOUNDS,
        "periodicity": 2 * np.pi,
        "periodicity_bounds": "fixed",
    }
    expsine = ExpSineSquared(**expsine_kwargs)
    const_kwargs = {
        "constant_value": constant_value,
        "constant_value_bounds": CONST_BOUNDS,
    }
    const = ConstantKernel(**const_kwargs)
    kernel = const * expsine
    gpr_kwargs = {
        "n_restarts_optimizer": 0,  # Default; 1 opt. cycle
        "normalize_y": True,
    }
    gpr = GaussianProcessRegressor(kernel=kernel, **gpr_kwargs)
    return gpr


def fleck_acquisition_func(
    temperature: float, en_range: float, mean: np.ndarray, variance: np.ndarray
) -> np.ndarray:
    mean = mean - mean.min()
    beta = en_range / (KBAU * temperature)
    p = 1 / np.sqrt(2 * np.pi * variance)
    term1 = beta * np.exp(-beta * mean)
    term2 = p * np.log(p)
    fleck_func = term1 - term2

    # Shift minimum to 0.0
    # TODO: do shifting outside of this function, because this is only required for plotting
    fleck_func -= fleck_func.min()

    return fleck_func


@dataclasses.dataclass
class GPRStatus:
    cycle: int
    grid: np.ndarray
    # Predicted energies and standard deviations
    energies: np.ndarray
    std: np.ndarray
    x_train: np.ndarray
    y_train: np.ndarray
    ind_next: int
    acq_func: np.ndarray
    temperature: float

    @property
    def x_next(self):
        return self.grid[self.ind_next][0]

    def dump_potential(self, out_dir=Path(".")):
        fn = out_dir / "potential_pred.dat"
        np.savetxt(
            fn,
            np.stack((self.grid.flatten(), self.energies, self.std), axis=1),
            header="Grid/rad pot/Eh std/Eh",
        )
        return fn


def run_gpr(
    grid,
    wrapped,
    callback=None,
    max_cycles=50,
    en_thresh=1e-5,
    en_range=50 / AU2KJPERMOL,
    temperature=298.15,
) -> GPRStatus:
    x_train = list()
    y_train = list()

    # Calculate the initial geometry first at 0.0.
    x_next = PERIOD_LOW

    @ignore_warnings(category=ConvergenceWarning)
    def objective(trial):
        length_scale = trial.suggest_float("length_scale", *LS_BOUNDS)
        constant_value = trial.suggest_float("constant_value", *CONST_BOUNDS)

        gpr = get_gpr(length_scale, constant_value)
        gpr.fit(np.array(x_train).reshape(-1, 1), y_train)
        _, std_pred = gpr.predict(np.array(grid).reshape(-1, 1), return_std=True)
        return std_pred.sum()

    prev_best_params = None
    # Start of macro loop
    for cur_cycle in range(max_cycles):
        x_train.append(x_next)
        # Evaluate new point
        y_train.append(wrapped(x_next))

        # print(f"@ There are {len(x_train)} training points in macro cycle {cur_cycle}.")

        study = optuna.create_study(direction="minimize")
        if prev_best_params is not None:
            study.enqueue_trial(prev_best_params)
        study.optimize(objective, n_trials=25)

        prev_best_params = study.best_params

        gpr = get_gpr(**study.best_params)
        gpr.fit(np.array(x_train).reshape(-1, 1), y_train)
        y_pred, std_pred = gpr.predict(grid, return_std=True)
        print(
            f"@@@ Cycle {cur_cycle:03d}, max. predicted standard deviation: "
            f"{std_pred.max(): >12.6e}"
        )

        # Update the selected energy range after some points have been calculated.
        if len(x_train) > 5:
            en_range = max(y_train) - min(y_train)
        # Add new trial point to enlarge the training set using the Fleck
        # acquisition function.
        fleck_func = fleck_acquisition_func(
            temperature,
            en_range,
            mean=y_pred,
            variance=std_pred**2,
        )
        max_ind = int(fleck_func.argmax())

        gpr_status = GPRStatus(
            cycle=cur_cycle,
            grid=grid,
            energies=y_pred,
            std=std_pred,
            x_train=np.array(x_train),
            y_train=np.array(y_train),
            ind_next=max_ind,
            acq_func=fleck_func,
            temperature=temperature,
        )

        if callback is not None:
            callback(gpr_status)

        # Check for convergence
        max_std = std_pred.max()
        # TODO: one could also check for the convergence of the ground state
        # energy/energies of some wavefunction(s).
        if max_std <= en_thresh:
            print(f"Converged in macro cycle {cur_cycle}!")
            break

        # Pick next trial point
        x_next = float(grid[max_ind, 0])
    # End of macro loop
    return gpr_status
