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

    # def __init__(self, *args, **kwargs):
        # # initial_displacement is already called in the constructor
        # # of IRC. As we override initial_displacement and make use of
        # # the DWI object there we already define it here.
        # self.dwi = DWI()

        # super().__init__(*args, **kwargs)

    def prepare(self, *args, **kwargs):
        super().prepare(*args, **kwargs)

        # Initialize the distance weighted interpolator with the data
        # from the initial displacement.
        self.dwi = DWI()
        mw_grad = self.mw_gradient
        energy = self.energy
        mw_hessian = self.mw_hessian
        self.dwi.update(self.mw_coords, energy, mw_grad, mw_hessian)

    # def initial_displacement(self):
        # step = super().initial_displacement()

        # TODO: Hessian update
        # dH, _ = bfgs_update(self.ts_hessian, dx, dg)

        # Initialize the distance weighted interpolator with the data
        # from the initial displacement.
        # mw_grad = self.mw_gradient
        # energy = self.energy
        # mw_hessian = self.mw_hessian
        # self.dwi.update(self.mw_coords, energy, mw_grad, mw_hessian)

        # return step

    def step(self):
        mw_grad = self.mw_gradient
        mw_grad_norm = np.linalg.norm(mw_grad)

        # Simple euler integration
        euler_step_length = self.step_length / 250
        mw_hessian = self.mw_hessian
        # print("mw_hessian")
        # print(mw_hessian)
        # print("init coords", self.mw_coords)
        # print("init grad", mw_grad)


        taylor_mw_grad = mw_grad.copy()
        tnmwg = np.linalg.norm(taylor_mw_grad)
        def taylor_grad(step):
            nmwhs = np.linalg.norm(mw_hessian @ step)
            # print(f"\tnorm(g)={tnmwg:.4f}, norm(H@s)={np.linalg.norm(mw_hessian @ step):.4f}")
            return taylor_mw_grad + mw_hessian @ step

        init_mw_coords = self.mw_coords.copy()
        euler_mw_coords = self.mw_coords
        m_sqrt = np.sqrt(self.geometry.masses_rep)
        # print("m_sqrt", m_sqrt)

        # calc = self.geometry.calculator
        # ax = calc.ax
        # ax.scatter(*(init_mw_coords / m_sqrt).T[:2], s=20, label="Init. displ")
        # e_coords = list()
        for i in range(500):
            # e_coords.append(euler_mw_coords.copy())
            # Calculate step length in non-mass-weighted coordinates
            cur_length = np.linalg.norm((euler_mw_coords - init_mw_coords) / m_sqrt)
            if i % 50 == 0:
                print(f"{i:02d}: cur_length={cur_length:.6f}")
                dir_ = -mw_grad / np.linalg.norm(mw_grad)
                print(f"dir: {dir_}")

            # Check if we achieved the desired step length up to some threshold
            # if (abs(cur_length - self.step_length) < euler_step_length) \
               # or cur_length > self.step_length:
            if cur_length > self.step_length:
                print(f"Predictor Euler integration converged with Δs={cur_length:.4f}!")
                break
            step_ = euler_step_length * -mw_grad / np.linalg.norm(mw_grad)
            euler_mw_coords += step_
            step = euler_mw_coords - init_mw_coords
            mw_grad = taylor_grad(step)

            # step_ = euler_step_length * -grad_ / np.linalg.norm(grad_)
        else:
            # raise Exception("Handle case where Euler integration does not converge. "
                            # "Kästner paper suggests that IRC is probably converged."
            # )
            self.mw_coords = euler_mw_coords
            return

        # for i in range(500):
            # e_coords.append(euler_mw_coords.copy())
            # # Calculate step length in non-mass-weighted coordinates
            # cur_length = np.linalg.norm((euler_mw_coords - init_mw_coords) / m_sqrt)
            # if i % 50 == 0:
                # print(f"{i:02d}: cur_length={cur_length:.6f}")

            # # Check if we achieved the desired step length up to some threshold
            # if (abs(cur_length - self.step_length) < euler_step_length) \
               # or cur_length > self.step_length:
                # print(f"Predictor Euler integration converged with Δs={cur_length:.4f}!")
                # break
            # step_ = euler_step_length * -mw_grad / np.linalg.norm(mw_grad)
            # euler_mw_coords += step_
            # step = euler_mw_coords - init_mw_coords
            # mw_grad = taylor_grad(step)

            # # step_ = euler_step_length * -grad_ / np.linalg.norm(grad_)
        # else:
            # raise Exception("Handle case where Euler integration does not converge. "
                            # "Kästner paper suggests that IRC is probably converged."
            # )
        self.mw_coords = euler_mw_coords
        # predictor_step_length = cur_length
        predictor_step_length = np.linalg.norm(euler_mw_coords - init_mw_coords)

        # e_coords = np.array(e_coords) / m_sqrt
        # ax.plot(*e_coords.T[:2], "r-")
        # ax.legend()
        # plt.show()
        # import pdb; pdb.set_trace()

        # Calculate energy and gradient at new predictor geometry.
        # These results will be added to the DWI for use in the corrector
        # step.
        mw_grad = self.mw_gradient
        energy = self.energy
        mw_hessian = self.mw_hessian
        self.dwi.update(self.mw_coords, energy, mw_grad, mw_hessian)

        all_coords = list()
        richardson = dict()
        errors = list()
        for k in range(10):
            points = 10*(2**k) + 1
            corr_step_length  = predictor_step_length / (points - 1)
            # print("corr_step_length", corr_step_length)
            cur_coords = init_mw_coords.copy()
            k_coords = list()
            length = 0
            while True:
                k_coords.append(cur_coords.copy())
                if length >= predictor_step_length:
                    # print(f"Converged! length={length:.8f}, length-step={length-corr_step_length:.8f}")
                    print(f"Converged! length={length:.8f}")
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
                print(f"k={k},j={j}")
                richardson[(k, j)] = ((2**j) * richardson[(k, j-1)] - richardson[(k-1, j-1)]) \
                                     / (2**j-1)
            if k > 0:
                # RMS of coordinates
                # error = np.sqrt(np.mean((richardson[(k, k)] - richardson[(k-1, k-1)])**2))
                # Error estimate according to Numerical Recipes Eq. (17.3.9)
                error = np.linalg.norm(richardson[(k,k)] - richardson[(k,k-1)])
                # print(f"\terror={error:.8e}, error_={error_:.8e}")
                # print(f"\terror={error:.8e}, err={err:.8e}")
                print(f"\terror={error:.8e}")
                errors.append(error)
                # if error <= 1e-6:
                if error <= 1e-4:
                    break
            all_coords.append(np.array(k_coords))
        else:
            raise Exception("Richardson did not converge!")
        
        self.mw_coords = richardson[(k,k)]

        print("Richardson table")
        # pprint(richardson)
        print()
        print("Errors")
        erros = np.array(errors)
        print(errors)

        # ax.plot(*coords.T[:2], "r-")

        self.irc_mw_coords.append(self.mw_coords)
        print("######COORDS")
        print(self.irc_mw_coords)
        print("######COORDS")
        # calc = self.geometry.calculator
        # calc.plot()
        # ax = calc.ax
        # cs = np.array(self.irc_mw_coords) / m_sqrt
        # import pdb; pdb.set_trace()
        # ax.plot(*cs.T[:2], "ro-")
        # plt.show()
