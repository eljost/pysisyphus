# [1] https://pubs.rsc.org/en/content/articlepdf/2002/cp/b108658h
#     Farkas, 2001, GDIIS
# [2] https://aip.scitation.org/doi/abs/10.1063/1.4878944
#     Schaefer, 2014, DIIS + Dimer

from collections import namedtuple
from functools import partial

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import rms
from pysisyphus.optimizers.closures import lbfgs_closure_ as lbfgs_closure


RotResult = namedtuple("RotResult",
                       "geom1 geom2 N C coords",
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
    # We don't need a calculator on geom2

    return geom1, geom2


def get_f_rot(f0, f1, N):
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

        for i, e1 in enumerate(err_vecs):
            for j in range(i, len(err_vecs)):
                e2 = err_vecs[j]
                A[i, j] = e1.dot(e2)
                A[j, i] = A[i, j]
        det = np.linalg.norm(A)
        print(f"det(A)={det:.6f}")
        cr = np.linalg.solve(A, np.ones(A.shape[0]))
        if any(np.abs(cr) > 1e8):
            break
        last_cr = cr
    if last_cr is None:
        raise DIISError("DIIS failed!")
    used_last = len(last_cr)
    cs = last_cr / np.sum(last_cr)
    print(f"DIIS with {used_last} vectors. Coeffs: ", cs)
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


def rotate_dimer(geom0, geom1, geom2, N, R, f_thresh=2.5e-3, max_cycles=10,
                 alpha=0.05):
    rot_optimizer = get_rot_optimizer(alpha=alpha)
    f0 = geom0.forces

    coords = list()
    for i in range(max_cycles):
        f1 = geom1.forces
        f_rot = get_f_rot(f0, f1, N)
        f_rot_rms = np.linalg.norm(f_rot)
        C = curvature(f0, f1, N, R)
        print(f"Cycle {i:02d}: rms(f_rot)={f_rot_rms:.06f}, C={C: .06f}")
        if f_rot_rms < f_thresh:
            print("Converged")
            break
        step = rot_optimizer(f_rot, geom1.coords)
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


def get_f_trans(geom0, N, C):
    f0 = geom0.forces

    f_parallel = f0.dot(N)*N
    if C > 0:
        f_trans = -f_parallel
    else:
        f_trans = f0 - 2*f_parallel

    return f_trans


def restrict_step_length(_, step, max_length):
    norm = np.linalg.norm(step)
    if norm > max_length:
        direction = step / norm
        step = direction * max_length
    return step


def dimer_method(geom0, N, R, calc_getter, max_cycles=50, f_thresh=1e-3,
                 max_step=0.3,
                 rot_kwargs=None, trans_kwargs=None):
    if rot_kwargs is None:
        rot_kwargs = {}
    if trans_kwargs is None:
        trans_kwargs = {}

    geom1, geom2 = get_dimer_ends(geom0, N, R, calc_getter)
    coords = [(geom1.coords, geom0.coords, geom2.coords), ]
    restrict_step = partial(restrict_step_length, max_length=max_step)
    trans_optimizer = lbfgs_closure(lambda _, *args: get_f_trans(*args),
                                    restrict_step=restrict_step)
    for i in range(max_cycles):
        f0 = geom0.forces
        f0_rms = rms(f0)
        print(f"{i:0d} rms(f0)={f0_rms:.6f}")
        if f0_rms < f_thresh:
            print("Converged!")
            break

        # Rotation
        rot_result = rotate_dimer(geom0, geom1, geom2, N, R, **rot_kwargs)
        geom1, geom2, N, C, _ = rot_result

        # Translation
        step, f_trans = trans_optimizer(geom0.coords, geom0, N, C)
        new_coords0 = geom0.coords + step
        geom0.coords = new_coords0

        # Steepet descent
        # f_trans = get_f_trans(geom0, N, C)
        # alpha = 0.5
        # step = alpha * f_trans
        # new_coords0 = geom0.coords + alpha*step
        # geom0.coords = new_coords0

        update_dimer_ends(geom0, geom1, geom2, N, R)

        coords.append(
            (geom1.coords, geom0.coords, geom2.coords)
        )
    # TODO: return dimerresults
    return coords
