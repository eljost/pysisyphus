#!/usr/bin/env python3

import time

import numpy as np
from scipy.optimize import minimize

from pysisyphus.helpers import check_for_stop_sign
from pysisyphus.optimizers.Optimizer import Optimizer


class StopOptException(Exception):
    pass


class SciPyOptimizer(Optimizer):

    def __init__(self, geometry, method="l-bfgs-b", **kwargs):
        super(SciPyOptimizer, self).__init__(geometry, **kwargs)

        if self.align:
            print("Ignoring align in SciPyOptimizer.")
        self.method = method
        self.options = {
            "disp": True,
            "maxiter": self.max_cycles,
            "gtol": 1e-3,
        }

    def callback(self, xk):
        self.cur_cycle += 1
        forces = self.geometry.forces
        step = self.coords[-1] - xk

        self.steps.append(step)

        if self.is_cos:
            self.tangents.append(self.geometry.get_tangents())

        self.check_convergence()

        if self.dump:
            self.write_cycle_to_file()

        if self.is_zts:
            self.geometry.reparametrize()

        if check_for_stop_sign():
            raise StopOptException()

    def fun(self, coords):
        start_time = time.time()

        self.coords.append(self.geometry.coords)
        self.geometry.coords = coords
        forces = self.geometry.forces
        self.forces.append(forces)
        self.energies.append(self.geometry.energy)
        forces_rms = np.sqrt(np.mean(np.square(forces)))

        end_time = time.time()
        elapsed_seconds = end_time - start_time
        self.cycle_times.append(elapsed_seconds)

        # gradient = -forces
        return forces_rms, -forces

    def run(self):
        # self.print_header()
        x0 = self.geometry.coords
        try:
            self.opt_res = minimize(self.fun, x0, jac=True, method=self.method,
                               callback=self.callback, options=self.options)
            if self.opt_res.success:
                print("Converged!")
            else:
                print("Didn't converge.")
        except StopOptException:
            self.log("Caught StopOptException. Stopping.")
