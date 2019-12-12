#!/usr/bin/env python3


import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.XTB import XTB
from pysisyphus.helpers import geom_from_xyz_file
from pysisyphus.optimizers.hessian_updates import bfgs_update, flowchart_update



def rms(arr):
    return np.sqrt(np.mean(arr**2))

def max_(arr):
    return np.abs(arr.max())


def rfo(gradient, H, trust=None):
    H_aug = np.bmat(
                ((H, gradient[:,None]),
                 (gradient[None,:], [[0]]))
    )
    eigvals, eigvecs = np.linalg.eigh((H_aug+H_aug.T)/2)
    aug_step = np.asarray(eigvecs[:,0]).flatten()
    lambda_ = aug_step[-1]
    step = aug_step[:-1] / lambda_
    if trust:
        norm = np.linalg.norm(step)
        if np.linalg.norm(step) > trust:
            step = step / norm * trust
    return step


def update_trust_radius(trust_radius, coeff, last_step_norm):
    # Nocedal, Numerical optimization Chapter 4, Algorithm 4.1
    if coeff < 0.25:
        # trust_radius = max(trust_radius/4, 0.001)
        trust_radius = max(trust_radius/4, 0.1)#0.001)
    # Only increase trust radius if last step norm was at least 80% of it
    # See [5], Appendix, step size and direction control
    # elif coeff > 0.75 and (last_step_norm >= .8*trust_radius):
    elif coeff > 0.75 and (abs(last_step_norm-trust_radius) < 1e-6):
        trust_radius = min(trust_radius*1.414, 2)
    # elif coeff > 0.75 and (last_step_norm >= .9*trust_radius):
        # trust_radius = min(trust_radius*1.414, 2)
    return trust_radius


def run():
    # geom = AnaPot.get_geom((-0.366, 2.03, 0))
    # H = geom.mw_hessian
    # geom = geom_from_xyz_file("azetidine_guess.xyz", coord_type="redund")
    geom = geom_from_xyz_file("azetidine_guess.xyz")
    geom.set_calculator(XTB())
    # H = geom.mw_hessian
    # H = geom.hessian
    H = geom.get_initial_hessian()
    import pdb; pdb.set_trace()
    M = geom.mm_sqrt_inv
    trust = 0.8

    max_cycles = 75
    steps = list()
    gradients = list()
    coords = list()
    energies = list()
    pred_changes = list()

    trj = open("opt.trj", "w")
    for i in range(max_cycles):
        coords.append(geom.coords.copy())
        trj.write(geom.as_xyz()+"\n")
        g = geom.gradient
        gradients.append(g)
        energies.append(geom.energy)

        if i > 0:
            last_step_norm = np.linalg.norm(steps[-1])
            pred = pred_changes[-1]
            actual = energies[-1] - energies[-2]
            # predicted will be defined when i > 0
            coeff = predicted / actual  # noqa: F821
            trust = update_trust_radius(trust, coeff, last_step_norm)
            # Hess update
            dg = gradients[-1] - gradients[-2]
            dx = steps[-1]
            # dH, _ = bfgs_update(H, dx, dg)
            dH, _ = flowchart_update(H, dx, dg)
            H = H + dH
            # H = geom.hessian

        # Convert gradient to normal mode gradient
        # cart_step = rfo(g, H, trust=trust)

        # eigvals, eigvecsT = np.linalg.eigh(H)
        # gq = eigvecsT.T.dot(g)
        # H_diag = np.diag(eigvals)
        # # Convert gradient to normal mode gradient
        # dq = rfo(gq, H_diag)#, trust=trust)
        # cart_step = eigvecsT.dot(dq)
        # norm = np.linalg.norm(cart_step)
        # if norm > trust:
            # cart_step = cart_step / norm * trust

        mwH = M.dot(H).dot(M)
        vm, vemT = np.linalg.eigh(mwH)
        S = M.dot(vemT)
        gqm = S.T.dot(g)
        Hm_diag = np.diag(vm)
        dqm = rfo(gqm, Hm_diag)
        cart_step = S.dot(dqm)
        norm = np.linalg.norm(cart_step)
        if norm > trust:
            cart_step = cart_step / norm * trust

        step_norm = np.linalg.norm(cart_step)
        # print(f"norm(step)={step_norm:.6f}")
        rms_g = rms(g)
        max_g = max_(g)
        rms_s = rms(cart_step)
        max_s = max_(cart_step)
        norm_g = np.linalg.norm(g)
        # print(f"Cycle {i:02d}: "
              # f"\tmax(grad)={max_g:.6f} rms(grad)={rms_g:.6f} "
              # f"max(step)={max_s:.6f} rms(step)={rms_s:.6f}")
        print(f"{i:02d}: {max_g:.6f} {rms_g:.6f} {max_s:.6f} {rms_s:.6f} {norm_g:.6f} {geom.energy:.6f} {trust:.6f}")

        converged = (max_g <= 4.5e-4) and (rms_g <= 3e-4) #\
                    # and (max_s <= 1.8e-3) and (rms_s <= 1.2e-3)
        if converged:
            print("Converged!")
            break


        steps.append(cart_step)
        new_coords = geom.coords + cart_step
        geom.coords = new_coords

        predicted = cart_step.dot(g) + 0.5 * cart_step.dot(H).dot(cart_step)
        pred_changes.append(predicted)

    print(f"Energy: {geom.energy:.8f}")
    pot = geom.calculator
    # pot.plot()
    # coords_arr = np.array(coords)
    # pot.ax.plot(coords_arr[:,0], coords_arr[:,1], "o-")
    # plt.show()
    trj.close()
    


if __name__ == "__main__":
    run()
