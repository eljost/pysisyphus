import numpy as np


# [1]  https://doi.org/10.1063/1.1809574


def lanczos(coords, grad_getter, dx=5e-3, dl=1e-2, guess=None, max_cycles=25,
            reortho=True, logger=None):
    """Lanczos method to determine smallest eigenvalue & -vector.

    See [1] for description of algorithm.
    """

    def log(msg):
        if logger is not None:
            logger.debug(msg)

    log("Lanczos Algorithm")
    r_prev = guess
    if r_prev is None:
        r_prev = np.random.rand(coords.size)
    beta_prev = np.linalg.norm(r_prev)
    q_prev = np.zeros_like(r_prev)

    alphas = list()
    betas = list()
    w_mins = list()
    Q = np.zeros((coords.size, max_cycles))
    qs = list()
    # Gradient at current coordinates; does not change.
    grad_l = grad_getter(coords)
    for i in range(max_cycles):
        # Normalize q
        q = r_prev / beta_prev
        Q[:,i] = q
        qs.append(q.copy())
        # Approximate action of Hessian on q (u = Hq) by finite differences.
        #
        # Gradient at perturbed coordinates
        grad_k = grad_getter(coords + dx*q)
        u = (grad_k - grad_l) / dx
        # Residue
        r = u - beta_prev*q_prev
        alpha = q.dot(r)
        r -= alpha*q
        # Reorthogonalization of r against the present Lanczos vectors in Q
        if reortho:
            r -= Q.dot(Q.T.dot(r))
        beta = np.linalg.norm(r)

        alphas.append(alpha)
        betas.append(beta)
        size = len(alphas)
        # Construct tri-diagonal matrix T
        T = np.zeros((size, size))
        diag_inds = np.diag_indices(size)
        T[diag_inds] = alphas
        if len(alphas) > 1:
            for j, b in enumerate(betas[:-1], 1):
                k = j-1
                T[j,k] = b
                T[k,j] = b

        # Values for next cycle
        beta_prev = beta
        q_prev = q
        r_prev = r
        # Diagonalize T
        w, v = np.linalg.eigh(T)
        w_min = w[0]
        log(f"Cycle {i: >3d}: w_min={w_min: .6f}")

        # Check eigenvalue convergence
        if (i > 0) and (abs((w_min - w_mins[-1])/w_mins[-1]) < dl):
            log("Converged")
            break

        w_mins.append(w_min)
    v_min = v[:,0]

    # Form eigenvector from linear combination of Lanczos vectors
    eigenvector = (v_min[:,None] * qs).sum(axis=0)
    eigenvector /= np.linalg.norm(eigenvector)
    return w_min, eigenvector


def geom_lanczos(geom, *args, **kwargs):
    """Wraps Lanczos algorithm for use with Geometry objects."""
    coords = geom.coords.copy()

    def grad_getter(coords):
        results = geom.get_energy_and_forces_at(coords)
        return -results["forces"]

    return lanczos(coords, grad_getter, *args, **kwargs)
