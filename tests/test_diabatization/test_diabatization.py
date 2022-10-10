import pytest
import yaml

from pysisyphus.drivers import dq_diabatization_from_run_dict


def test_d_diabatization(this_dir):
    yaml_fn = this_dir / "guanine_indole_dia.yaml"
    with open(yaml_fn) as handle:
        run_dict = yaml.load(handle, Loader=yaml.SafeLoader)
    dia_res = dq_diabatization_from_run_dict(run_dict)
    assert dia_res.is_converged
    assert dia_res.P == pytest.approx(172.4953, abs=1e-3)
