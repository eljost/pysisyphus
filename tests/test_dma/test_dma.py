import numpy as np
import pytest

from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.dma import run_dma


@pytest.mark.parametrize(
    "wf_fn, ref_multis, sites",
    (
        # Match DMA exactly
        ("09_o_atom_def2_tzvp.fchk", "09_o_atom_ref.npy", None),
        ("00_h2.fchk", "00_h2_ref.npy", None),
        ("05_h_porb.fchk", "05_h_porb_ref.npy", None),
        (
            "15_o_double_cation.fchk",
            "15_o_double_cation_ref.npy",
            (0.5, 0.5, 0.5),
        ),
        ("18_h2_small.fchk", "18_h2_small_ref.npy", None),
        ("17_h2_prims.fchk", "17_h2_prims_ref.npy", None),
        ("07_oh_minus.fchk", "07_oh_minus_ref.npy", None),
        ("24_oh_minus_tzvp.fchk", "24_oh_minus_tzvp_ref.npy", None),
        (
            "24_oh_minus_tzvp.fchk",
            "24_oh_minus_tzvp_ref_shift.npy",
            (0.0, 0.0, 0.0),
        ),
        ("13_h2o_sto3g.fchk", "13_h2o_sto3g_ref.npy", None),
        # Without GDMA and only DMA (shift 0.0) DMA produces an unsymmetric result
        # for ethene.
        (
            "06_ethene_from_orca_geom.fchk",
            "06_ethene_from_orca_geom_ref.npy",
            None,
        ),
        ("16_methane_def2svp.fchk", "16_methane_def2svp_ref.npy", None),
        ("23_he.fchk", "23_he_ref.npy", None),
    ),
)
def test_dma(wf_fn, ref_multis, sites, this_dir):
    data_dir = this_dir / "data"
    wf = Wavefunction.from_file(data_dir / wf_fn)
    print(wf_fn, "\n", wf)

    distributed_multipoles = run_dma(wf, sites, switch=0.0)
    if ref_multis is not None:
        ref_multis = data_dir / ref_multis
        if not ref_multis.exists():
            np.save(ref_multis, distributed_multipoles)
            print(f"Dumped multipoles to '{ref_multis}'")
        multis_ref = np.load(ref_multis)
        np.testing.assert_allclose(distributed_multipoles, multis_ref, atol=1e-14)
        print("Multipoles match!")
    else:
        print("Skipped check!")


@pytest.mark.parametrize(
    "wf_fn, ref_multis",
    (("24_oh_minus_tzvp.fchk", "24_oh_minus_tzvp_switch40_ref.npy"),),
)
def test_gdma(wf_fn, ref_multis, this_dir):
    data_dir = this_dir / "data"
    wf = Wavefunction.from_file(data_dir / wf_fn)
    print(wf_fn, "\n", wf)

    distributed_multipoles = run_dma(wf, switch=4.0)
    if ref_multis is not None:
        ref_multis = data_dir / ref_multis
        if not ref_multis.exists():
            np.save(ref_multis, distributed_multipoles)
            print(f"Dumped multipoles to '{ref_multis}'")
        multis_ref = np.load(ref_multis)
        np.testing.assert_allclose(distributed_multipoles, multis_ref, atol=1e-14)
        print("Multipoles match!")
