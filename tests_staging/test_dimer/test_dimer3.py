#!/usr/bin/env python3

# [1] https://pubs.rsc.org/en/content/articlepdf/2002/cp/b108658h
# [2] http://www.acme.byu.edu/wp-content/uploads/2014/09/Krylov1.pdf
# [3] https://www.youtube.com/watch?v=0t7WJybTmFg
# [4] Improved mode following
#     https://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b01216

from collections import namedtuple

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.spatial.distance import pdist
from scipy.optimize import minimize

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.AnaPot import AnaPot

np.set_printoptions(suppress=True, precision=4)


RotResult = namedtuple("RotResult",
                       "geom1 geom2 N C coords",
)
TransResult = namedtuple("TransResult",
                         "geom0",
)


def update_dimer_ends(geom0, geom1, geom2, N, R):
    coords1 = geom0.coords + R*N
    geom1.coords = coords1
    coords2 = geom0.coords - R*N
    geom2.coords = coords2


def get_dimer_ends(geom0, N, R, calc_getter):
    dummy_coords = np.zeros_like(geom0.coords)
    geom1 = Geometry(geom0.atoms, dummy_coords)
    geom2 = Geometry(geom0.atoms, dummy_coords)
    update_dimer_ends(geom0, geom1, geom2, N, R)
    geom1.set_calculator(calc_getter())

    return geom1, geom2


def get_frot(f0, f1, N):
    frot = 2*(f1 - f0) - 2*(f1 - f0).dot(N)*N
    return frot


class DIISError(Exception):
    pass


def diis(error_vecs, coords, forces):
    cycles = len(error_vecs)

    last_cr = None
    for use_last in range(2, cycles+1):
        A = np.zeros((use_last, use_last))
        # Start with the latest point and add previous points one by one
        err_vecs = np.array(error_vecs[::-1])[:use_last]
        # Scale error vector so that norm of smallest error vector is 1
        err_norms = np.linalg.norm(err_vecs, axis=1)
        scale_factor = 1 / err_norms.min()
        err_vecs  *= scale_factor
        err_norms = np.linalg.norm(err_vecs, axis=1)

        for i, e1 in enumerate(err_vecs):
            for j in range(i, len(err_vecs)):
                e2 = err_vecs[j]
                A[i, j] = e1.dot(e2)
                A[j, i] = A[i, j]
        cr = np.linalg.solve(A, np.ones(A.shape[0]))
        if any(np.abs(cr) > 1e8):
            break
        last_cr = cr
    if last_cr is None:
        raise DIISError("DIIS failed!")
    used_last = len(last_cr)
    cs = last_cr / np.sum(last_cr)
    last_coords = coords[::-1][:used_last]
    last_error_vecs = error_vecs[::-1][:used_last]
    last_forces = forces[::-1][:used_last]

    # Form linear combinations
    coords = np.sum(cs[:,None]*last_coords, axis=0)
    error = np.sum(cs[:,None]*last_error_vecs, axis=0)
    force = np.sum(cs[:,None]*last_forces, axis=0)
    return error, coords, force


def get_rot_optimizer(alpha=0.05):
    cycle = 0
    error_vecs = list()
    all_coords = list()
    all_forces = list()
    def get_rot_step(frot, coords1):
        nonlocal cycle

        error_vecs.append(frot.copy())
        all_coords.append(coords1.copy())
        all_forces.append(frot.copy())

        # Plain steepest descent step
        step = alpha*frot
        # Try DIIS from the second iteration onwards
        if cycle > 0:
            try:
                error, coords1_, frot_ = diis(error_vecs, all_coords, all_forces)
                # Inter-/extrapolated coords + step from inter-/extrapolated
                # rotational forces.
                new_diis_coords = coords1_ + alpha*frot_
                # Determine actual step as difference between the current coordinates
                # and the DIIS coordinates.
                step = new_diis_coords - coords1
            except DIISError:
                pass
        cycle += 1
        return step
    return get_rot_step


def update_rotated_endpoints(geom0, geom1, geom2, geom1_step, R):
    tmp_coords = geom1.coords + geom1_step
    geom1.coords = tmp_coords

    N = geom1.coords - geom0.coords
    N /= np.linalg.norm(N)
    # Reconstrain dimer onto hypersphere
    x1 = geom0.coords + R*N
    x2 = geom0.coords - R*N
    geom1.coords = x1
    geom2.coords = x2
    return geom1, geom2, N


def curvature(f0, f1, N, R):
    C = 2*(f0-f1).dot(N)/(2*R)
    return C


def rotate_dimer(geom0, geom1, geom2, N, R):
    rot_optimizer = get_rot_optimizer()
    f0 = geom0.forces

    coords = list()
    for i in range(15):
        f1 = geom1.forces
        frot = get_frot(f0, f1, N)
        norm = np.linalg.norm(frot)
        C = curvature(f0, f1, N, R)
        print(f"Cycle {i:02d}: norm(forces)={norm:.06f}, C={C: .06f}")
        if norm < 1e-3:
            print("Converged")
            break
        step = rot_optimizer(frot, geom1.coords)
        geom1, geom2, N = update_rotated_endpoints(geom0, geom1, geom2, step, R)
        coords.append(
            (geom1.coords.copy(), geom0.coords.copy(), geom2.coords.copy())
        )
    rot_result = RotResult(
        geom1=geom1,
        geom2=geom2,
        N=N,
        C=C,
        coords=coords,
    )
    return rot_result


def translate_dimer(geom0, N, C):
    f0 = geom0.forces

    f_parallel = f0.dot(N)*N
    if C > 0:
        f_eff = -f_parallel
    else:
        f_eff = f0 - 2*f_parallel
    alpha = 0.5
    step = alpha * f_eff
    new_coords0 = geom0.coords + alpha*step
    geom0.coords = new_coords0
    trans_result = TransResult(
                    geom0=geom0,
    )
    return trans_result


def run_dimer():
    x0 = (-0.5767, 1.6810, 0)
    geom0 = AnaPot.get_geom(x0)

    N = np.array((-0.9, 0.43, 0))
    N /= np.linalg.norm(N)
    R = 0.125

    geom1, geom2 = get_dimer_ends(geom0, N, R, AnaPot)
    coords = list()
    for i in range(35):
        f0 = geom0.forces
        norm_f0 = np.linalg.norm(f0)
        print(f"{i:0d} {norm_f0:.6f}")
        if norm_f0 < 1e-3:
            print("Converged!")
            break
        # Rotation
        rot_result = rotate_dimer(geom0, geom1, geom2, N, R)
        geom1, geom2, N, C, _ = rot_result
        # Translation
        trans_result = translate_dimer(geom0, N, C)
        update_dimer_ends(geom0, geom1, geom2, N, R)
        # Backup for plotting
        cs = (geom1.coords, geom0.coords, geom2.coords)
        coords.append(cs)

    # geom0.calculator.plot_eigenvalue_structure(grid=100)
    coords = np.array(coords)
    pot = AnaPot()
    pot.plot()
    ax = pot.ax
    for i, rot_cycle in enumerate(coords):
        ax.plot(*rot_cycle.T[:2], "o-", label=f"Cycle {i:02d}")
    ax.legend()
    plt.show()


def lanczos(geom, dx=1e-4, dl=1e-2, r_prev=None, max_cycles=10):
    if r_prev is None:
        r_prev = np.random.rand(*geom.coords.shape)
    beta_prev = np.linalg.norm(r_prev)
    q_prev = np.zeros_like(r_prev)

    alphas = list()
    betas = list()
    w_mins = list()
    for i in range(max_cycles):
        q = r_prev / beta_prev
        # Finite difference
        grad_l = geom.gradient
        # Perturbed geometry
        coords_k = geom.coords + dx*q
        results_k = geom.get_energy_and_forces_at(coords_k)
        grad_k = -results_k["forces"]
        # u = H.dot(q)
        # Finite difference alternative
        u = (grad_k - grad_l) / dx
        r = u - beta_prev*q_prev
        alpha = q.dot(r)
        r = r - alpha*q
        beta = np.linalg.norm(r)

        alphas.append(alpha)
        betas.append(beta)
        size = len(alphas)
        T = np.zeros((size, size))
        diag_inds = np.diag_indices(size)
        T[diag_inds] = alphas
        if len(alphas) > 1:
            for i, b in enumerate(betas[:-1], 1):
                j = i-1
                T[i,j] = b
                T[j,i] = b
        beta_prev = beta
        q_prev = q
        r_prev = r
        w, v = np.linalg.eig(T)
        min_ind = w.argmin()
        w_min = w[min_ind]
        w_mins.append(w_min)
        print(f"w_min={w_min: .6f}")
        try:
            w_diff = abs(w_min - w_mins[-2])
            if w_diff < dl:
                print("Converged")
                break
        except IndexError:
            pass
    w_min = w_mins[-1]
    v_min = v[:,min_ind]
    return w_min, v_min


def run_lanczos():
    x0 = (-0.5767, 1.6810, 0)
    geom0 = AnaPot.get_geom(x0)
    w_min, v_min = lanczos(geom0, dx=1e-5)
    H = geom0.hessian
    w, v = np.linalg.eigh(H)
    w_min_ref = w[0]
    print(f"w_min_ref={w_min_ref:.6f}")


def hypersphere_displacement(size, samples=50, radius=0.3):
    # https://www.sciencedirect.com/science/article/pii/S0378437116001230
    coords = np.random.rand(samples, size) - 0.5
    norms = np.linalg.norm(coords, axis=1)
    # coords = radius * coords / norms[:,None]
    coords = radius * coords / norms[:,None]
    # print("norms", np.linalg.norm(coords, axis=1))

    cdm = pdist(coords)

    def energy(coords, shape=coords.shape):
        _ = coords.reshape(shape)
        norms = np.linalg.norm(_, axis=1)
        cdm = pdist(coords.reshape(shape))
        return np.sum(1/cdm) + np.sum((norms-radius)**8)
    res = minimize(energy, coords)
    print(res)
    x = res.x.reshape(coords.shape)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(*coords.T)
    ax.scatter(*x.T, color="r")

    xnorms = np.linalg.norm(x, axis=1)
    print("xnorms", xnorms)

    plt.show()




def run():
    # run_dimer()
    run_lanczos()
    # hypersphere_displacement(3)


if __name__ == "__main__":
    run()
