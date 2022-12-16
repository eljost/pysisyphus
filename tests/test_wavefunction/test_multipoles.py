import numpy as np
import pytest

from pysisyphus.calculators.ORCA import parse_orca_cis
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.excited_states import (
    norm_ci_coeffs,
    get_state_to_state_transition_density,
)

np.set_printoptions(suppress=True, precision=5)


@pytest.mark.parametrize(
    "base, unrestricted", (
        ("00_pyrrole_rhf", False),
        ("01_pyrrole_uhf", True),
        )
)
def test_transition_multipole_moments(base, unrestricted, this_dir):
    base_fn = this_dir / "pyrrole" / base
    json_fn = base_fn.with_suffix(".json")
    cis_fn = base_fn.with_suffix(".cis")
    wf = Wavefunction.from_orca_json(json_fn)
    Xa, Ya, Xb, Yb = parse_orca_cis(cis_fn)
    Xa, Ya = norm_ci_coeffs(Xa, Ya)

    gs_es_args = [Xa, ]
    # Reference values
    if unrestricted:
        Xa, Ya, Xb, Yb = norm_ci_coeffs(Xa, Yb, Xb, Yb)
        gs_es_args = [Xa, Xb]
        tdms_ref = np.array(
            (
                (0.00000, -0.00000, 0.00000),
                (0.00000, -0.00000, 0.00000),
                (0.00000, -0.00000, 0.00001),
            )
        )
        es2es_tdms_ref = np.array(
            (
                (-0.00000, 0.00000, -0.09534),
                (0.00000, 0.00000, 0.00000),
                (-0.00000, 0.00000, 0.00000),
            )
        )
    else:
        tdms_ref = np.array(
            (
                (0.00000, -0.00000, 0.00001),
                (0.00000, -0.00000, -0.05814),
                (-0.00079, 0.01290, -0.00000),
            )
        )
        es2es_tdms_ref = np.array(
            (
                (0.54050, 0.03649, -0.00000),
                (0.00000, -0.00000, 0.00000),
                (-0.00000, 0.00000, 0.07424),
            )
        )
    
    # Calculate and check GS-ES transition moments
    tdms = wf.get_transition_dipole_moment(*gs_es_args)
    np.testing.assert_allclose(tdms, tdms_ref, atol=1e-5)

    # Calculate and check ES-ES transition moments
    nstates = Xa.shape[0]
    estdens_a = list()
    estdens_b = list()
    # Loop over excited state pairs
    for i in range(nstates):
        for j in range(i + 1, nstates):
            Xai = Xa[i]
            Xaj = Xa[j]
            estdena = get_state_to_state_transition_density(Xai, Xaj)
            estdens_a.append(estdena)
            if unrestricted:
                Xbi = Xb[i]
                Xbj = Xb[j]
                estdenb = get_state_to_state_transition_density(Xbi, Xbj)
                estdens_b.append(estdenb)
    estdens_a = np.array(estdens_a)
    estdens_b = np.array(estdens_b)
    args = [
        estdens_a,
    ]
    if unrestricted:
        args.append(estdens_b)
    es_tdms = wf.get_transition_dipole_moment(*args, full=True)
    np.testing.assert_allclose(es_tdms, es2es_tdms_ref, atol=1e-5)
