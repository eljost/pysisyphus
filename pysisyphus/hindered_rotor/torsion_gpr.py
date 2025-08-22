import dataclasses
from pathlib import Path

import numpy as np
import optuna
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process.kernels import ExpSineSquared, ConstantKernel
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.utils._testing import ignore_warnings
from sklearn.model_selection import cross_val_score
from sklearn.metrics import root_mean_squared_error

from pysisyphus import numerov
from pysisyphus.constants import KBAU
from pysisyphus.drivers import boltzmann
from pysisyphus.helpers import check_for_end_sign
from pysisyphus.partfuncs import partfuncs as pf

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
    mass: float
    temperature: float
    eigvals_absolute: np.ndarray
    eigvecs: np.ndarray
    partfunc: float
    weights: np.ndarray

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


def numerov_callback(grid: np.ndarray, energies: np.ndarray, mass: float):
    # Drop last data point because our potential is periodic
    grid_cut = grid[:-1].flatten()
    energies_cut = energies[:-1]
    en_min = energies_cut.min()
    energies_cut = energies_cut - en_min

    def energy_getter(i, x):
        return energies_cut[i]

    # Run Numerov method to obtain wavefunctions and associated eigenvalues
    eigvals, eigvecs = numerov.run(grid_cut, energy_getter, mass, periodic=True)
    eigvals_absolute = eigvals + en_min
    return eigvals_absolute, eigvecs


@ignore_warnings(category=ConvergenceWarning)
def objective(trial, x_train, y_arr, start_validation: int = 10):
    length_scale = trial.suggest_float("length_scale", *LS_BOUNDS)
    constant_value = trial.suggest_float("constant_value", *CONST_BOUNDS)

    gpr = get_gpr(length_scale, constant_value)
    # After some training data is available we do cross validation,
    # to actually converge against the true potential energy curve.
    if len(x_train) >= start_validation:
        score = cross_val_score(
            gpr,
            x_train,
            y_arr,
            cv=5,
            scoring="neg_mean_squared_error",
        )
        error = -score.mean()
    # In the beginning we use the negative of the log marginal likelihood
    # as error metric.
    else:
        gpr.fit(x_train, y_arr)
        error = -gpr.log_marginal_likelihood(gpr.kernel_.theta)
    return error


def run_gpr(
    grid,
    wrapped,
    mass: float,
    callback=None,
    max_cycles: int = 50,
    partfunc_thresh: float = 1e-4,
    std_thresh: float = 1e-4,
    temperature: float = 298.15,
) -> GPRStatus:
    assert max_cycles > 0
    en_range = KBAU * temperature
    x_train = list()
    y_train = list()

    # Calculate the initial geometry first at 0.0.
    x_next = PERIOD_LOW

    prev_partfunc = float("nan")
    prev_best_params = None
    # Start of macro loop
    for cur_cycle in range(max_cycles):
        x_train.append(x_next)
        # Evaluate new point
        y_train.append(wrapped(x_next))

        x_arr = np.reshape(x_train, (-1, 1))
        y_arr = np.array(y_train)
        # Shift energies into range [0, 1] by subtracting the minimum energy
        # and dividing by the maximum energy.
        y_min = y_arr.min()
        y_arr = y_arr - y_min
        y_max = y_arr.max()
        if cur_cycle > 0:
            # TODO: enable additional check to verify that the energy range is
            # sensible?!
            # assert y_max > 1e-8
            y_arr = y_arr / y_max

        # print(f"@ There are {len(x_train)} training points in macro cycle {cur_cycle}.")

        study = optuna.create_study(direction="minimize")
        # TODO: check if the two lines below are needed
        if prev_best_params is not None:
            study.enqueue_trial(prev_best_params)
        # Optimize hyperparameters w/ optuna
        study.optimize(
            lambda trial: objective(trial, x_arr, y_arr),
            n_trials=25,
        )
        # Update best parameters for next cycle
        prev_best_params = study.best_params

        # Predict some data
        gpr = get_gpr(**study.best_params)
        gpr.fit(x_arr, y_arr)
        y_pred_norm, std_pred_norm = gpr.predict(grid, return_std=True)
        # Also predict exactly at the training points so we can later calculate the
        # root mean square error between predictions and actual training points.
        y_train_pred = gpr.predict(np.reshape(x_train, (-1, 1)))

        # Recover actual energies and standard deviation
        y_pred = y_min + (y_pred_norm * y_max)
        y_train_pred = y_min + (y_train_pred * y_max)
        if cur_cycle > 0:
            std_pred = std_pred_norm * y_max
        else:
            std_pred = std_pred_norm
        rms_err = root_mean_squared_error(y_train, y_train_pred)

        # Update the selected energy range after some points have been calculated.
        if len(x_train) > 5:
            en_range = max(y_train) - min(y_train)
        if y_max > 0.0:
            en_range_norm = en_range / y_max
        else:
            en_range_norm = KBAU * temperature
        # Add new trial point to enlarge the training set using the Fleck
        # acquisition function.
        fleck_func = fleck_acquisition_func(
            temperature,
            en_range_norm,
            mean=y_pred_norm,
            variance=std_pred_norm**2,
        )
        # Integer index that determines the next grid point to be calculated to
        # enlarge the training set.
        max_ind = int(fleck_func.argmax())
        x_next = float(grid[max_ind, 0])

        eigvals_absolute, eigvecs = numerov_callback(grid, y_pred, mass)
        eigvals_shift = eigvals_absolute - eigvals_absolute.min()
        weights = boltzmann.boltzmann_weights(eigvals_shift, temperature)
        partfunc = pf.sos_partfunc(eigvals_shift, temperature)
        partfunc_diff = abs(partfunc - prev_partfunc)
        prev_partfunc = partfunc
        print(
            f"@@@ Cycle {cur_cycle:03d}, max(std_pred):{std_pred.max():8.4e}, "
            f"rms(error):{rms_err:8.4e} pf={partfunc: >8.4f}, |Δpf|={partfunc_diff: >8.4e}"
        )

        gpr_status = GPRStatus(
            cycle=cur_cycle,
            grid=grid,
            energies=y_pred,
            std=std_pred,
            x_train=np.array(x_train),
            y_train=np.array(y_train),
            ind_next=max_ind,
            acq_func=fleck_func,
            mass=mass,
            temperature=temperature,
            eigvals_absolute=eigvals_absolute,
            eigvecs=eigvecs,
            partfunc=partfunc,
            weights=weights,
        )
        if callback is not None:
            callback(gpr_status)

        # Indicate convergence when the predicted standard deviation and the change
        # of the partition function are sufficiently small.
        if (partfunc_diff <= partfunc_thresh) and (std_pred.max() <= std_thresh):
            print(f"Converged in macro cycle {cur_cycle}!")
            break

        # Check if the operator want's to manually stop the calculation
        sign = check_for_end_sign(add_signs="converged_hr")
        if sign == "converged_hr":
            print("Operator indicated convergence!")
            break

    # End of macro loop
    return gpr_status
