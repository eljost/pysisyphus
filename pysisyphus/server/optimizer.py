import types

import numpy as np

from pysisyphus.optimizers.poly_fit import poly_line_search


def get_step_func(key="sd", alpha=0.1, gdiis=False, line_search=False):
    alpha0 = alpha
    print("alpha0", alpha0, "gdiis", gdiis, "line_search", line_search)

    def sd_step_func(self, coords, energy, grad):
        grad_norm = np.linalg.norm(grad)
        direction = -grad / grad_norm

        alpha = min(grad_norm, alpha0)
        step = alpha * direction

        return step

    prev_grad = None  # lgtm [py/unused-local-variable]
    prev_p = None  # lgtm [py/unused-local-variable]
    prev_energy = None  # lgtm [py/unused-local-variable] 
    prev_step = None  # lgtm [py/unused-local-variable]
    def cg_step_func(self, coords, energy, grad):
        nonlocal prev_grad
        nonlocal prev_p
        nonlocal prev_energy
        nonlocal prev_step

        grad_norm = np.linalg.norm(grad)
        alpha = min(grad_norm, alpha0)

        # Steepst descent in the first cycle
        if prev_grad is None:
            prev_grad = grad
            prev_p = -grad
            prev_energy = energy

            direction = -grad / grad_norm
            step = alpha * direction
            prev_step = step
            return step

        if line_search:
            fit_step, fit_grad, fit_energy = poly_line_search(energy, prev_energy,
                                                              grad, prev_grad,
                                                              prev_step, coords)
        # Polak-Ribierre
        beta = grad.dot(grad-prev_grad) / prev_grad.dot(prev_grad)
        p = -grad + beta*prev_p

        if line_search and fit_step is not None:
            step = fit_step
            grad = fit_grad
            energy = fit_energy
        else:
            step = np.linalg.norm(prev_step) * p / np.linalg.norm(p)

        prev_grad = grad
        prev_energy = energy
        prev_p = p
        prev_step = step

        return step

    step_funcs = {
        "sd": sd_step_func,
        "cg": cg_step_func,
    }
    return step_funcs[key]


class OptState:

    def __init__(self, key="sd", alpha=0.3, gdiis=False, line_search=True):
        self.key = key
        self.alpha = alpha
        self.gdiis = gdiis
        self.line_search = line_search

        self.coords = list()
        self.grads = list()
        self.energies = list()
        self.steps = list()

        step_func_kwargs = {
            "key": self.key,
            "alpha": self.alpha,
            "gdiis": self.gdiis,
            "line_search": self.line_search,
        }
        step_func = get_step_func(**step_func_kwargs)
        # Register step_func
        self.step_func = types.MethodType(step_func, self)

        print("Initialized OptState")

    def step(self, coords, energy, grad):
        self.coords.append(coords)
        self.energies.append(energy)
        self.grads.append(grad)

        step = self.step_func(coords, energy, grad)

        return step, "SUCCESS"
