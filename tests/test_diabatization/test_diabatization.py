import numpy as np
import pytest
import yaml

from pysisyphus.drivers import dq_diabatization_from_run_dict
from pysisyphus.diabatization import driver as dia_driver
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
    )
    # Only one key will be present
    key = list(dia_results.keys())[0]
    dia_res = dia_results[key]
    dia_ens = dia_res.dia_ens
    np.testing.assert_allclose(dia_ens, [4.50496744, 4.50496744], atol=1e-8)
    assert dia_res.couplings[(0, 1)] == pytest.approx(0.2793929, abs=1e-7)
