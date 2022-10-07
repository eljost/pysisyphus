import pytest

from pysisyphus.drivers import dq_diabatization_from_yaml


def test_d_diabatization(this_dir):
    dia_res = dq_diabatization_from_yaml(this_dir / "guanine_indole_dia.yaml")
    assert dia_res.is_converged
    assert dia_res.P == pytest.approx(172.4953, abs=1e-3)
