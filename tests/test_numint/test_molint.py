import time

import numpy as np
import pytest

from pysisyphus.numint import get_mol_grid
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.density_numba import eval_density


@pytest.mark.parametrize(
    "fn, N_ref",
    (
        ("lib:02_ne.json", 10.0),
        ("lib:04_h2.json", 2.0),
        ("lib:orca_h2o_def2svp.json", 10.0),
        ("lib:orca_benzene_quad.json", 42.0),
        ("lib:orca_ch4_qzvpp.json", 10.0),
    ),
)
def test_molint(fn, N_ref):
    wf = Wavefunction.from_file(fn)
    # Get molecular grid
    mol_grid = get_mol_grid(wf.atoms, wf.coords3d)
    xyz = mol_grid.xyz

    # Evaluate density at grid points
    rho = np.zeros(len(xyz))
    shellstructs = wf.shells.as_numba_shellstructs()
    precontr = wf.shells.cart2sph_coeffs.T @ wf.shells.P_sph.T
    start = time.time()
    eval_density(shellstructs, xyz, wf.P_tot, precontr, rho)
    dens_dur = time.time() - start

    N = (rho * mol_grid.weights).sum()
    print(f"{xyz.shape=}, N_integrated={N:.6f}, {dens_dur=:.4f} s")
    assert N == pytest.approx(N_ref, abs=3e-4)
