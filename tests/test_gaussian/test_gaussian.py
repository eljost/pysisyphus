import numpy as np
import pytest

from pysisyphus.calculators import Gaussian16
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


@pytest.mark.parametrize("route", ("symmetry", "nosymm", "force", "opt", "freq", "irc"))
def test_invalid_keywords(route):
    with pytest.raises(AssertionError):
        Gaussian16(route)


@using("gaussian16")
def test_all_energies():
    geom = Gaussian16.geom_from_fn("lib:h2o.xyz", route="hf def2svp td=(nstates=2)")
    all_energies = geom.all_energies
    len(all_energies) == 3
    assert all_energies[2] == pytest.approx(-75.5569539)


@using("gaussian16")
@pytest.mark.parametrize(
    "root, ref_energy, ref_dpm", (
    (0, -98.185807, (0.0, 0.0, 0.7423065)),
    # Funnily enough the HF cation has a degenerate GS, so D0 and D1 have
    # identical energies and dipole moments. Higher multipoles differ though.
    (1, -98.185807, (0.0, 0.0, 0.7423065)),
    (4, -97.547199, (0.0, 0.0, 2.0825352)),
    )
)
def test_relaxed_density_hf_cation(root, ref_energy, ref_dpm):
    calc_kwargs = {
        "charge": 1,
        "mult": 2,
        "route": "hf sto-3g cis(nstates=4)",
        "root": root,
        "base_name": f"root_{root}",
    }

    geom = geom_loader("lib:hf_hf_sto3g_opt.xyz")
    calc = Gaussian16(**calc_kwargs)
    geom.set_calculator(calc)

    result = geom.calc_relaxed_density(root)
    dens = result["density"]
    wf = geom.wavefunction
    dens_key = wf.set_relaxed_density(root, dens)
    # Gaussian seems to put the origin for the multipole calculations at (0.0, 0.0, 0.0)
    with wf.current_density(dens_key):
        dpm = wf.get_dipole_moment(origin=np.zeros(3))

    energy = result["energy"]
    assert energy == pytest.approx(ref_energy)

    atol = 1e-5
    np.testing.assert_allclose(dpm, ref_dpm, atol=atol)
