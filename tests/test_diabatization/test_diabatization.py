import numpy as np
import pytest
import yaml

from pysisyphus.diabatization.driver_yaml import dq_diabatization_from_run_dict
from pysisyphus.diabatization import driver as dia_driver
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


def test_d_diabatization(this_dir):
    yaml_fn = this_dir / "guanine_indole_dia.yaml"
    with open(yaml_fn) as handle:
        run_dict = yaml.load(handle, Loader=yaml.SafeLoader)
    dia_res = dq_diabatization_from_run_dict(run_dict)
    assert dia_res.is_converged
    assert dia_res.P == pytest.approx(172.4953, abs=1e-3)


@pytest.mark.parametrize(
    "dia_kinds",
    (
        pytest.param(
            dia_driver.DiaKind.EDMISTON_RUEDENBERG, marks=using("pysisyphus_addons")
        ),
        (dia_driver.DiaKind.BOYS),
    ),
)
def test_diabatization_driver(dia_kinds, this_dir):
    """Example from http://dx.doi.org/10.1063/1.3148777 (BeNa2 cation)."""
    inp_dir = this_dir / "00_bena2_cat"
    out_dir = this_dir / "00_bena2_cat_out"
    base = inp_dir / "00_bena2_dcat.log"
    base_name = base.stem
    wf, all_ens, Xa, Ya, Xb, Yb = dia_driver.parse_orca(base)
    states = [1, 2]

    dia_inp = dia_driver.DiaInput(
        wf=wf,
        all_ens=all_ens,
        Xa=Xa,
        Ya=Ya,
        Xb=Xb,
        Yb=Yb,
        states=states,
        base_name=base_name,
    )
    # Update states by calculating overlaps if requested
    dia_results, cube_fns = dia_driver.run_dia(
        dia_inp,
        dia_kinds=dia_kinds,
        cube_kinds=dia_driver.CubeKind.NONE,
        out_dir=out_dir,
        force=True,
        h5_fn=out_dir / "00_bena2_cat_dia_result.h5",
    )
    # Only one key will be present
    key = list(dia_results.keys())[0]
    dia_res = dia_results[key]
    dia_ens = dia_res.dia_ens
    np.testing.assert_allclose(dia_ens, [4.50496744, 4.50496744], atol=1e-8)
    assert dia_res.couplings[(0, 1)] == pytest.approx(0.2793929, abs=1e-7)


@using("pyscf")
def test_pyscf_diabatization_driver(this_dir):
    """Example for diabatization of a PySCF HF/TDA calculation using the new driver.

    This test replicates the BeNa^{2+} example from  section III on p. 234102-8
    in 10.1063/1.3148777. The level of theory is HF/6-31G* and CIS for the excited
    states. A diabatization of the first two excited singlet states is carried out.

    """
    from pysisyphus.calculators.PySCF import PySCF

    geom = geom_loader(this_dir / "00_bena2_cat" / "00_bena2_dcat.xyz")
    calc = PySCF(
        charge=2,
        mult=1,
        basis="631Gs",
        method="tdahf",
        nroots=5,
        verbose=4,
    )
    geom.set_calculator(calc)

    # Accessing geom.wavefunction would start a plain GS calculation.
    # In this example, we request a GS + ES calculation via geom.all_energies,
    # and then create a Wavefunction object from the (stored) mf-object in the
    # PySCF calculator.
    #
    # wf = geom.wavefunction
    #

    # Ground state energy and excited state energies
    all_ens = geom.all_energies
    # Pysisyphus wavefunction object
    wf = calc.get_stored_wavefunction()
    # Unrelaxed transition density matrices
    Xa, Ya, Xb, Yb = geom.td_1tdms
    # Do diabatization of the first two excited states
    states = [1, 2]

    base_name = "00_bena2_cat_pyscf"
    out_dir = this_dir / f"{base_name}_out"
    dia_inp = dia_driver.DiaInput(
        wf=wf,
        all_ens=all_ens,
        Xa=Xa,
        Ya=Ya,
        Xb=Xb,
        Yb=Yb,
        states=states,
        base_name=base_name,
    )
    # Update states by calculating overlaps if requested
    dia_results, cube_fns = dia_driver.run_dia(
        dia_inp,
        dia_kinds=dia_driver.DiaKind.BOYS,
        cube_kinds=dia_driver.CubeKind.NONE,
        out_dir=out_dir,
        force=True,
        h5_fn=out_dir / f"{base_name}_dia_result.h5",
    )
    # Only one key will be present
    key = list(dia_results.keys())[0]
    dia_res = dia_results[key]
    dia_ens = dia_res.dia_ens

    # Check diabatic exciation energies and diabatic couplings against reference
    np.testing.assert_allclose(dia_ens, [4.50497, 4.50497], atol=1e-5)
    assert dia_res.couplings[(0, 1)] == pytest.approx(0.27938, abs=1e-5)
