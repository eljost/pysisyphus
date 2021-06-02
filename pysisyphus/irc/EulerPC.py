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
# [4  ] https://aip.scitation.org/doi/pdf/10.1063/1.3593456?class=pdf
#       Hratchian, Frisch, 2011
# [6  ] https://aip.scitation.org/doi/10.1063/1.3593456<Paste>
#       Hratchian, Frisch
# 	Further improvements for DWI; not implemented

import time

import numpy as np
from scipy.integrate import solve_ivp

from pysisyphus.helpers import rms
from pysisyphus.io.hessian import save_hessian
from pysisyphus.irc.DWI import DWI
from pysisyphus.irc.IRC import IRC
from pysisyphus.optimizers.hessian_updates import bfgs_update, bofill_update


class EulerPC(IRC):
    def __init__(
        self,
        *args,
        hessian_recalc=None,
        hessian_update="bofill",
        hessian_init="calc",
        max_pred_steps=500,
        loose_cycles=3,
        dump_dwi=False,
        scipy_method=None,
        corr_func="mbs",
        **kwargs,
    ):
        super().__init__(*args, hessian_init=hessian_init, **kwargs)

        self.hessian_recalc = hessian_recalc
        self.hessian_update = {
            "bfgs": bfgs_update,
            "bofill": bofill_update,
        }
        self.hessian_update_func = self.hessian_update[hessian_update]
        self.max_pred_steps = int(max_pred_steps)
        self.loose_cycles = loose_cycles
        self.dump_dwi = dump_dwi

        self.scipy_method = scipy_method
        corr_funcs = {
            "mbs": self.corrector_step,
            # "scipy_v1": self.scipy_corrector_step_v1,
            "scipy": self.scipy_corrector_step,
        }
        if (self.scipy_method is not None) and corr_func != "scipy":
            self.log(
                f"scipy_method={scipy_method} given, but corr_func={corr_func}. "
                "Setting corr_func=scipy to use scipy integrator."
            )
            corr_func = "scipy"
        self.corr_func = corr_funcs[corr_func]

    def prepare(self, *args, **kwargs):
        super().prepare(*args, **kwargs)

        # Initialize the distance weighted interpolator with the data
        # from the initial displacement.
        self.dwi = DWI()
        mw_grad = self.mw_gradient
        energy = self.energy

        # Store starting information for distances weighted interpolation
        self.dwi.update(self.mw_coords, energy, mw_grad, self.mw_hessian.copy())

        if self.downhill:
            return

        # Do a first Hessian update with the information between the TS
        # and the initially displaced geometry.
        dx = self.mw_coords - self.ts_mw_coords
        dg = mw_grad - self.ts_mw_gradient
        dH, key = self.hessian_update_func(self.mw_hessian, dx, dg)
        self.log(f"Did {key} hessian update.")
        self.mw_hessian += dH

    def get_integration_length_func(self, init_mw_coords):
        def get_integration_length(cur_mw_coords):
            """Returns length of integration path done in mass-weighted coordinates
            in un-mass-weighted coordinates."""
            return np.linalg.norm((cur_mw_coords - init_mw_coords) / self.m_sqrt)

        return get_integration_length

    def step(self):
        ##################
        # PREDICTOR STEP #
        ##################

        mw_grad = self.mw_gradient
        energy = self.energy

        if self.cur_cycle > 0:
            if self.hessian_recalc and (self.cur_cycle % self.hessian_recalc == 0):
                self.mw_hessian = self.geometry.mw_hessian
                h5_fn = f"hess_calc_irc_{self.direction}_cyc{self.cur_cycle}.h5"
                save_hessian(h5_fn, self.geometry)
                self.log("Calculated excact hessian")
            else:
                dx = self.mw_coords - self.irc_mw_coords[-2]
                dg = mw_grad - self.irc_mw_gradients[-2]
                dH, key = self.hessian_update_func(self.mw_hessian, dx, dg)
                self.mw_hessian += dH
                self.log(f"Did {key} hessian update before predictor step.")
            self.dwi.update(
                self.mw_coords.copy(), energy, mw_grad, self.mw_hessian.copy()
            )

        # Create a copy of the inital coordinates for the determination
        # of the actual step size in the predictor Euler integration.
        init_mw_coords = self.mw_coords.copy()

        get_integration_length = self.get_integration_length_func(init_mw_coords)

        # Calculate predictor Euler-integration step length. See get_conv_fact
        # method definition for a comment on this.
        conv_fact = self.get_conv_fact(mw_grad)
        euler_step_length = self.step_length / (self.max_pred_steps / conv_fact)

        def taylor_gradient(step):
            """Return gradient from Taylor expansion of energy to 2nd order."""
            return mw_grad + self.mw_hessian @ step

        # These variables will hold the coordinates and gradients along
        # the Euler integration and will be updated frequently.
        euler_mw_coords = self.mw_coords.copy()
        euler_mw_grad = mw_grad.copy()
        self.log(
            f"Predictor-Euler-integration with Δs={euler_step_length:.6f} "
            f"for up to {self.max_pred_steps} steps\n     #  |step|  d|step|"
        )
        prev_cur_length = 0.0
        for i in range(self.max_pred_steps):
            # Calculate step length in non-mass-weighted coordinates
            cur_length = get_integration_length(euler_mw_coords)
            if i % 50 == 0:
                diff = cur_length - prev_cur_length
                self.log(f"\t{i:03d}: {cur_length:.4f} Δ={diff:.4f}")
                prev_cur_length = cur_length

            # Check if we achieved the desired step length.
            if cur_length >= self.step_length:
                self.log(
                    "Predictor-Euler integration converged with "
                    f"Δs={cur_length:.4f} (desired Δs={self.step_length:.4f}) "
                    f"after {i+1} steps!"
                )
                break
            step_ = euler_step_length * -euler_mw_grad / np.linalg.norm(euler_mw_grad)
            euler_mw_coords += step_
            # Determine actual step by comparing the current and the initial coordinates
            euler_step = euler_mw_coords - init_mw_coords
            euler_mw_grad = taylor_gradient(euler_step)
        else:
            self.log(
                f"Predictor-Euler integration did not converge in {i+1} "
                f"steps. Δs={cur_length:.4f}."
            )

            # Check if we are already sufficiently converged. If so signal
            # convergence.
            self.mw_coords = euler_mw_coords

            # Use rms of gradient from taylor expansion for convergence check.
            euler_grad = self.unweight_vec(euler_mw_grad)
            rms_grad = rms(euler_grad)

            # Or check true gradient? But this would need an additional calculation,
            # so I disabled it for now.
            # rms_grad = rms(self.gradient)

            if self.cur_cycle < self.loose_cycles:
                self.log(
                    f"Current cycle {self.cur_cycle} is still in 'loose' mode.\n"
                    "Continuing IRC integration even though predictor integration "
                    f"did not succeed.\n{self.loose_cycles - self.cur_cycle - 1} loose "
                    "cycles remaining."
                )
            # elif rms_grad <= 5*self.rms_grad_thresh:
            elif rms_grad <= self.rms_grad_thresh:
                self.log("Sufficient convergence achieved on rms(grad)")
                self.converged = True
                return
        self.log("")

        # Calculate energy and gradient at new predicted geometry. Update the
        # hessian accordingly. These results will be added to the DWI for use
        # in the corrector step.
        self.mw_coords = euler_mw_coords
        self.log("Calculating energy and gradient at predictor step geometry.")
        mw_grad = self.mw_gradient
        energy = self.energy

        # Hessian update
        dx = self.mw_coords - self.irc_mw_coords[-1]
        dg = mw_grad - self.irc_mw_gradients[-1]
        dH, key = self.hessian_update_func(self.mw_hessian, dx, dg)
        self.mw_hessian += dH
        self.log(f"Did {key} hessian update after predictor step.\n")
        self.dwi.update(self.mw_coords.copy(), energy, mw_grad, self.mw_hessian.copy())
        if self.dump_dwi:
            self.dwi.dump(
                f"dwi_{self.cur_direction}_{self.cur_cycle:0{self.cycle_places}d}.h5"
            )

        corrected_mw_coords = self.corr_func(init_mw_coords, self.step_length, self.dwi)
        self.mw_coords = corrected_mw_coords
        corr_step_length = get_integration_length(self.mw_coords)
        self.log(f"Corrected unweighted step length: {corr_step_length:.6f}")

    def corrector_step(self, init_mw_coords, step_length, dwi):
        self.log("Corrector step using mBS integration")

        get_integration_length = self.get_integration_length_func(init_mw_coords)

        errors = list()
        richardson = dict()
        for k in range(15):
            points = 20 * (2 ** k)
            corr_step_length = step_length / (points - 1)
            cur_coords = init_mw_coords.copy()
            k_coords = list()
            cur_length = 0

            # Integrate until the desired spacing is reached
            while True:
                k_coords.append(cur_coords.copy())
                if abs(step_length - cur_length) < 0.5 * corr_step_length:
                    self.log(
                        f"\tk={k:02d} points={points: >4d} "
                        f"step_length={corr_step_length:.4f} Δs={cur_length:.4f}"
                    )
                    break

                energy, gradient = dwi.interpolate(cur_coords, gradient=True)
                cur_coords += corr_step_length * -gradient / np.linalg.norm(gradient)
                # cur_length += corr_step_length
                cur_length = get_integration_length(cur_coords)

                # Check for oscillation
                try:
                    prev_coords = k_coords[-2]
                    osc_norm = np.linalg.norm(cur_coords - prev_coords)
                    # TODO: Handle this by restarting everything with a smaller stepsize?
                    # Check 10.1039/c7cp03722h SI
                    if osc_norm <= corr_step_length:
                        self.log(
                            "\tDetected oscillation in Corrector-Euler "
                            f"integration for k={k:02d} and {points} points.\n"
                            "\tAborting corrector integration!"
                        )
                        return prev_coords
                except IndexError:
                    pass
            richardson[(k, 0)] = cur_coords

            # Refine using Richardson extrapolation
            # Set additional values using Richard extrapolation
            for j in range(1, k + 1):
                richardson[(k, j)] = (
                    (2 ** j) * richardson[(k, j - 1)] - richardson[(k - 1, j - 1)]
                ) / (2 ** j - 1)
            # Can only be done after the second successful integration
            if k > 0:
                # Error estimate according to Numerical Recipes Eq. (17.3.9).
                # We compare the last two entries/columns in the current row.
                # RMS error
                error = np.sqrt(
                    np.mean((richardson[(k, k)] - richardson[(k, k - 1)]) ** 2)
                )
                errors.append(error)
                if error <= 1e-5:
                    self.log(f"mBS integration converged (error={error:.4e})!")
                    break
        else:
            raise Exception("Richardson did not converge!")

        self.log(
            f"Returning corrected mass-weighted coordinates from richardson[({k},{k})]"
        )
        return richardson[(k, k)]

    def scipy_corrector_step(self, init_mw_coords, step_length, dwi):
        """Solve IRC equation dx/ds = -g/|g| on DWI PES in mass-weighted
        coordinates. Integration done until self.step_length in unweighted
        coordinates is achieved."""

        # Integrator
        if self.scipy_method is None:
            method = "Radau"
        else:
            method = self.scipy_method

        # Jacobian/hessian
        jac = None
        if method in "Radau BDF".split():
            # jac = dwi.hessians[0]
            jac = self.mw_hessian
        else:
            self.log(
                f"The chosen method {method} is probably a bad choice "
                "for integrating IRCs, as they it is not well suited "
                "for stiff problems. See SciPy documentation of "
                "solve_ivp() for more information."
            )

        self.log(f"Corrector step using SciPy and {method} integration")

        _, init_mw_grad = dwi.interpolate(init_mw_coords, gradient=True)
        conv_fact = self.get_conv_fact(init_mw_grad, min_fact=1.0)

        def fun(t, cur_mw_coords):
            energy, gradient = dwi.interpolate(cur_mw_coords, gradient=True)
            return -gradient / np.linalg.norm(gradient)

        t_span = (0, self.step_length * conv_fact)
        self.log(f"\tt_span={t_span}")
        int_start = time.time()
        ode_res = solve_ivp(
            fun=fun,
            t_span=t_span,
            y0=init_mw_coords,
            method=method,
            jac=jac,
        )
        int_end = time.time()
        int_duration = int_end - int_start

        attrs = "message nfev njev nlu sol status success".split()
        for attr in attrs:
            self.log(f"\t{attr}: {getattr(ode_res, attr)}")
        self.log(f"Corrector integration took {int_duration:.4} s")

        # Integration reached end of t_span
        # try:
        # assert ode_res.status == 0
        # except AssertionError:
        # import pdb; pdb.set_trace()
        corrected_mw_coords = ode_res.y[:, -1]
        return corrected_mw_coords

    # def scipy_corrector_step_v1(self, init_mw_coords, step_length, dwi):
    # get_integration_length = self.get_integration_length_func(init_mw_coords)

    # # Termination event
    # def integration_length_reached_event(t, cur_mw_coords):
    # cur_length = get_integration_length(cur_mw_coords)
    # return self.step_length - cur_length
    # integration_length_reached_event.terminal = True

    # def fun(t, cur_mw_coords):
    # energy, gradient = dwi.interpolate(cur_mw_coords, gradient=True)
    # return -gradient

    # ode_res = solve_ivp(
    # fun=fun,
    # t_span=(0, 10),
    # y0=init_mw_coords,
    # method="RK45",
    # events=integration_length_reached_event,
    # # first_step=self.step_length / 50,
    # # max_step=self.step_length / 50,
    # # dense_output=True,
    # )

    # # Expect termination event
    # # try:
    # # assert ode_res.status == 1
    # # except AssertionError:
    # # import pdb; pdb.set_trace()
    # import pdb; pdb.set_trace()
    # corrected_mw_coords = ode_res.y[:,-1]
    # return corrected_mw_coords
