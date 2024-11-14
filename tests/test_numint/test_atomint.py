import time

import numpy as np
import pytest

from pysisyphus.numint import get_atomic_grid
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.density_numba import eval_density


@pytest.mark.parametrize("kind", ("g1", "g2", "g3", "g4", "g5", "g6", "g7"))
@pytest.mark.parametrize(
    "fn, N_ref",
    (
        ("lib:00_h.json", 1.0),
        ("lib:01_he.json", 2.0),
        ("lib:02_ne.json", 10.0),
        ("lib:03_c.json", 6.0),
    ),
)
def test_atomint(kind, fn, N_ref):
    wf = Wavefunction.from_file(fn)
    atom = wf.atoms[0]
    origin = wf.coords3d[0]
    xyz, ww = get_atomic_grid(atom, origin, kind=kind)
    shellstructs = wf.shells.as_numba_shellstructs()
    P = wf.P_tot
    precontr = wf.shells.cart2sph_coeffs.T @ wf.shells.P_sph.T
    rho = np.zeros(len(xyz))
    start = time.time()
    eval_density(shellstructs, xyz, P, precontr, rho)
    dur = time.time() - start
    N = (rho * ww).sum()
    print(
        f"@ Atom={atom.capitalize(): >3s}: {xyz.shape[0]: >10d} points, {N=:.8f}, "
        f"density took {dur:.4e} s"
    )
    assert N == pytest.approx(N_ref, abs=5e-5)
