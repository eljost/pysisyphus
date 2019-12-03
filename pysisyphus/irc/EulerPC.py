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


import numpy as np

from pysisyphus.irc.DWI import DWI
from pysisyphus.irc.IRC import IRC
from pysisyphus.optimizers.hessian_updates import bfgs_update, bofill_update


class EulerPC(IRC):

    def __init__(self, *args, hessian_recalc=None, hessian_update="bofill", **kwargs):
        super().__init__(*args, **kwargs)

        self.hessian_recalc = hessian_recalc
        self.hessian_update = {
            "bfgs": bfgs_update,
            "bofill": bofill_update,
        }
        self.hessian_update_func = self.hessian_update[hessian_update]

    def prepare(self, *args, **kwargs):
        super().prepare(*args, **kwargs)

        # Initialize the distance weighted interpolator with the data
        # from the initial displacement.
        self.dwi = DWI()
        mw_grad = self.mw_gradient
        energy = self.energy
        self.mw_H = self.geometry.mass_weigh_hessian(self.ts_hessian)

        dx = self.mw_coords - self.ts_mw_coords
        dg = mw_grad - self.ts_mw_gradient
        dH, key = self.hessian_update_func(self.ts_hessian, dx, dg)
        self.log(f"Did {key} hessian update.")
        self.mw_H += dH
        self.dwi.update(self.mw_coords, energy, mw_grad, self.mw_H.copy())

    def step(self):
        ##################
        # PREDICTOR STEP #
        ##################

        mw_grad = self.mw_gradient
        mw_grad_norm = np.linalg.norm(mw_grad)

        # Simple euler integration
        euler_step_length = self.step_length / 250

        # TODO: Avoid recalculation of hessian
        # mw_hessian = self.mw_hessian

        # def taylor_gradient(step):
            # """Return gradient from Taylor expansion of energy to 2nd order."""
            # return mw_grad + mw_hessian @ step
        def taylor_gradient(step):
            """Return gradient from Taylor expansion of energy to 2nd order."""
            return mw_grad + self.mw_H @ step

        # Create a copy of the inital coordinates for the determination
        # of the actual step size in the predictor Euler integration.
        init_mw_coords = self.mw_coords.copy()
        m_sqrt = np.sqrt(self.geometry.masses_rep)
        def get_integration_length(cur_mw_coords):
            """Returns length of integration path done in mass-weighted coordinates
            in un-mass-weighted coordinates."""
            return np.linalg.norm((cur_mw_coords - init_mw_coords) / m_sqrt)

        # These variables will hold the coordinates and gradients along
        # the Euler integration and will be updated frequently.
        euler_mw_coords = self.mw_coords.copy()
        euler_mw_grad = mw_grad.copy()

        for i in range(500):
            # Calculate step length in non-mass-weighted coordinates
            cur_length = get_integration_length(euler_mw_coords)

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

        dx = self.mw_coords - self.irc_mw_coords[-1]
        dg = mw_grad - self.irc_mw_gradients[-1]
        if self.hessian_recalc and (self.cur_step % self.hessian_recalc == 0):
            self.mw_H = self.mw_hessian
            self.log("Calculated excact hessian")
        else:
            dH, key = self.hessian_update_func(self.mw_H, dx, dg)
            self.mw_H += dH
            self.log(f"Did {key} hessian update.")
        self.dwi.update(self.mw_coords.copy(), energy, mw_grad, self.mw_H.copy())

        ##################
        # CORRECTOR STEP #
        ##################

        all_coords = list()
        richardson = dict()
        errors = list()
        self.log("Starting mBS integration")
        for k in range(10):
            points = 10*(2**k) + 1
            corr_step_length  = self.step_length / (points - 1)
            cur_coords = init_mw_coords.copy()
            k_coords = list()
            cur_length = 0

            # Only use given points
            # for point in range(points):
                # k_coords.append(cur_coords.copy())
                # energy, gradient = self.dwi.interpolate(cur_coords, gradient=True)
                # cur_coords += corr_step_length * -gradient/np.linalg.norm(gradient)
            # cur_length = get_integration_length(cur_coords)
            # self.log(f"\tk={k:02d} points={points: >4d} "
                     # f"step_length={corr_step_length:.4f} Δs={cur_length:.4f}")

            # Integrate until the desired spacing is reached
            while True:
                k_coords.append(cur_coords.copy())
                if cur_length >= self.step_length:
                    self.log(f"\tk={k:02d} points={points: >4d} "
                             f"step_length={corr_step_length:.4f} Δs={cur_length:.4f}")
                    break
                energy, gradient = self.dwi.interpolate(cur_coords, gradient=True)
                cur_coords += corr_step_length * -gradient/np.linalg.norm(gradient)
                cur_length += corr_step_length
                # cur_length = get_integration_length(cur_coords)

                # Check for oscillation
                try:
                    prev_coords = k_coords[-2]
                    osc_norm = np.linalg.norm(cur_coords - prev_coords)
                    # TODO: handle this by restarting everyhting with a smaller stepsize.
                    # Check 10.1039/c7cp03722h SI
                    if osc_norm <= corr_step_length:
                        self.log("Detected oscillation. This is not handled right now!")
                        # assert False, "This case is not yet handled"
                        # break
                except IndexError:
                    pass
            richardson[(k, 0)] = cur_coords

            # Refine using Richardson extrapolation
            # Set additional values using Richard extrapolation
            for j in range(1, k+1):
                richardson[(k, j)] = ((2**j) * richardson[(k, j-1)] - richardson[(k-1, j-1)]) \
                                     / (2**j-1)
            if k > 0:
                # Error estimate according to Numerical Recipes Eq. (17.3.9)
                # error = np.linalg.norm(richardson[(k,k)] - richardson[(k,k-1)])

                # RMS error
                error = np.sqrt(np.mean((richardson[(k,k)] - richardson[(k,k-1)])**2))
                errors.append(error)
                if error <= 1e-5:
                    self.log(f"mBS integration converged (error={error:.4e})!")
                    break
            all_coords.append(np.array(k_coords))
        else:
            raise Exception("Richardson did not converge!")
        
        self.mw_coords = richardson[(k,k)]
