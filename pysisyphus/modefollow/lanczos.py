import numpy as np


def lanczos(coords, grad_getter, dx=5e-3, dl=1e-2, r_prev=None, max_cycles=25,
            logger=None):
    """Lanczos method to determine smallest eigenvalue & -vector."""

    def log(msg):
        if logger is not None:
            logger.debug(msg)

    log("Lanczos Algorithm")
    if r_prev is None:
        r_prev = np.random.rand(*coords.shape)
    beta_prev = np.linalg.norm(r_prev)
    q_prev = np.zeros_like(r_prev)

    alphas = list()
    betas = list()
    w_mins = list()
    qs = list()
    # Gradient at current coordinates; does not change.
    grad_l = grad_getter(coords)
    for i in range(max_cycles):
        log("Cycle {i: >3d}")
        # Normalize q
        q = r_prev / beta_prev
        qs.append(q.copy())
        # Approximate action of Hessian on q (u = Hq) by finite differences.
        #
        # Gradient at perturbed coordinates
        grad_k = grad_getter(coords + dx*q)
        u = (grad_k - grad_l) / dx
        # Residue
        r = u - beta_prev*q_prev
        alpha = q.dot(r)
        r = r - alpha*q
        beta = np.linalg.norm(r)

        alphas.append(alpha)
        betas.append(beta)
        size = len(alphas)
        # Construct tri-diagonal matrix T
        T = np.zeros((size, size))
        diag_inds = np.diag_indices(size)
        T[diag_inds] = alphas
        if len(alphas) > 1:
            for i, b in enumerate(betas[:-1], 1):
                j = i-1
                T[i,j] = b
                T[j,i] = b
        # Set values for next cycle
        beta_prev = beta
        q_prev = q
        r_prev = r
        # Diagonaliez T
        w, v = np.linalg.eigh(T)
        w_min = w[0]
        w_mins.append(w_min)
        print(f"w_min={w_min: .6f}")
        try:
            w_diff = abs(w_min - w_mins[-2])
            if w_diff < dl:
                print("Converged")
                break
        except IndexError:
            pass
    v_min = v[:,0]

    # Form eigenvector from linear combination of Lanczos vectors
    eigenvector = (v_min[:,None] * qs).sum(axis=0)
    eigenvector /= np.linalg.norm(eigenvector)
    return w_min, eigenvector


def geom_lanczos(geom, *args, **kwargs):
    coords = geom.coords.copy()
    def grad_getter(coords):
        results = geom.get_energy_and_forces_at(coords)
        return -results["forces"]
    return lanczos(coords, grad_getter, *args, **kwargs)


