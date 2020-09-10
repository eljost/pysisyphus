import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.helpers import geom_loader


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


def test_anapot():
    x0 = (-0.5767, 1.6810, 0)
    geom0 = AnaPot.get_geom(x0)
    w_min, v_min = lanczos(geom0, dx=1e-5)
    H = geom0.hessian
    w, v = np.linalg.eigh(H)
    w_min_ref = w[0]
    print(f"w_min_ref={w_min_ref:.6f}")


def test_hcn():
    geom = geom_loader("lib:hcn_iso_hf_sto3g_ts_opt.xyz")
    from pysisyphus.calculators.PySCF import PySCF
    calc = PySCF(pal=2, basis="sto3g")
    geom.set_calculator(calc)
    r_prev = np.array((-0.4, -0.1, -0.0, 0.0, 0.5, 0.0, 0.4, -0.2, 0.0))
    w_min, v_min = lanczos(geom, dx=1e-5, r_prev=r_prev)

    H = geom.hessian
    w, v = np.linalg.eigh(H)
    w_ref = w[0]
    v_ref = v[:,0]
    import pdb; pdb.set_trace()
    # print(
