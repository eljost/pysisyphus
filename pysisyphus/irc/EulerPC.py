#!/usr/bin/env python3


# [1  ] https://aip.scitation.org/doi/pdf/10.1063/1.3514202?class=pdf
#       Original EulerPC
#       Hratchian, Schlegel, 2010
# [2  ] https://aip.scitation.org/doi/pdf/10.1063/1.1724823?class=pdf
#       Original HPC
#       Hratchian, Schlegel, 2004
# [3  ] https://pubs.rsc.org/en/content/articlepdf/2017/cp/c7cp03722h
#       EulerPC re-implementation
#       Meisner, Kästner, 2017
# [3.1] http://www.rsc.org/suppdata/c7/cp/c7cp03722h/c7cp03722h1.pdf
#       Corresponding SI

from pprint import pprint


import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.Geometry import Geometry
# from pysisyphus.optimizers.hessian_updates import bfgs_update
from pysisyphus.irc.DWI import DWI
from pysisyphus.irc.IRC import IRC
from pysisyphus.optimizers.hessian_updates import bfgs_update


class EulerPC(IRC):

    def prepare(self, *args, **kwargs):
        super().prepare(*args, **kwargs)

        # Initialize the distance weighted interpolator with the data
        # from the initial displacement.
        self.dwi = DWI()
        mw_grad = self.mw_gradient
        energy = self.energy
        # TODO: Hessian update
        # dH, _ = bfgs_update(self.ts_hessian, dx, dg)
        mw_hessian = self.mw_hessian
        self.dwi.update(self.mw_coords, energy, mw_grad, mw_hessian)

    def step(self):
        ##################
        # PREDICTOR STEP #
        ##################

        mw_grad = self.mw_gradient
        mw_grad_norm = np.linalg.norm(mw_grad)

        # Simple euler integration
        euler_step_length = self.step_length / 250
        # TODO: Avoid recalculation of hessian
        mw_hessian = self.mw_hessian

        def taylor_gradient(step):
            """Return gradient from Taylor expansion of energy to 2nd order."""
            return mw_grad + mw_hessian @ step

        # Create a copy of the inital coordinates for the determination
        # of the actual step size in the predictor Euler integration.
        init_mw_coords = self.mw_coords.copy()

        # These variables will hold the coordinates and gradients along
        # the Euler integration and will be updated frequently.
        euler_mw_coords = self.mw_coords.copy()
        euler_mw_grad = mw_grad.copy()

        m_sqrt = np.sqrt(self.geometry.masses_rep)
        for i in range(500):
            # Calculate step length in non-mass-weighted coordinates
            cur_length = np.linalg.norm((euler_mw_coords - init_mw_coords) / m_sqrt)

            # Check if we achieved the desired step length
            if cur_length > self.step_length:
                self.log(f"Predictor Euler integration converged with Δs={cur_length:.4f}!")
                break
            step_ = euler_step_length * -euler_mw_grad / np.linalg.norm(euler_mw_grad)
            euler_mw_coords += step_
            # Determine actual step by comparing the current and the initial coordinates
            euler_step = euler_mw_coords - init_mw_coords
            euler_mw_grad = taylor_gradient(euler_step)
        else:
            # Assume convergence when predictor Euler integration does not converge.
            self.mw_coords = euler_mw_coords
            self.converged = True
            return

        # Calculate energy and gradient at new predictor geometry. These
        # results will be added to the DWI for use in the corrector step.
        self.mw_coords = euler_mw_coords
        mw_grad = self.mw_gradient
        energy = self.energy
        mw_hessian = self.mw_hessian
        self.dwi.update(self.mw_coords, energy, mw_grad, mw_hessian)

        ##################
        # CORRECTOR STEP #
        ##################

        all_coords = list()
        richardson = dict()
        errors = list()
        for k in range(10):
            points = 10*(2**k) + 1
            corr_step_length  = self.step_length / (points - 1)
            cur_coords = init_mw_coords.copy()
            k_coords = list()
            length = 0
            while True:
                k_coords.append(cur_coords.copy())
                if length >= self.step_length:
                    self.log(f"mBS integration with k={k} and "
                             f"step_length={corr_step_length:.6f} "
                             f"converged with total step length={length:.6f}")
                    break
                energy, gradient = self.dwi.interpolate(cur_coords, gradient=True)
                cur_coords += corr_step_length * -gradient/np.linalg.norm(gradient)
                length += corr_step_length
                # Check for oscillation
                try:
                    prev_coords = k_coords[-2]
                    osc_norm = np.linalg.norm(cur_coords - prev_coords)
                    if osc_norm <= corr_step_length:
                        print("Detected oscillation. Breaking!")
                        # TODO: handle this by restarting everyhting with a smaller stepsize.
                        # Check 10.1039/c7cp03722h SI
                        assert False, "This case is not yet handled"
                        break
                except IndexError:
                    pass
            richardson[(k, 0)] = cur_coords

            # Refine using Richardson extrapolation
            # Set additional values using Richard extrapolation
            for j in range(1, k+1):
                # print(f"k={k},j={j}")
                richardson[(k, j)] = ((2**j) * richardson[(k, j-1)] - richardson[(k-1, j-1)]) \
                                     / (2**j-1)
            if k > 0:
                # Error estimate according to Numerical Recipes Eq. (17.3.9)
                # error = np.linalg.norm(richardson[(k,k)] - richardson[(k,k-1)])
                # RMS error
                error = np.sqrt(np.mean((richardson[(k,k)] - richardson[(k,k-1)])**2))
                errors.append(error)
                if error <= 1e-5:
                    self.log(f"Extrapolation converged (error={error:.4e})!")
                    break
            all_coords.append(np.array(k_coords))
        else:
            raise Exception("Richardson did not converge!")
        
        self.mw_coords = richardson[(k,k)]
