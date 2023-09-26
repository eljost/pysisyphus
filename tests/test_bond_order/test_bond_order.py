import numpy as np
import pytest

from pysisyphus.config import WF_LIB_DIR
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.bond_order import (
    improved_mayer_bond_order,
    report_bond_orders,
)


def test_mayer():
    wf = Wavefunction.from_file(WF_LIB_DIR / "indole_bp86_def2svp.bson")
    BOs = improved_mayer_bond_order(wf)
    report_bond_orders(BOs)

    np.testing.assert_allclose(BOs, BOs.T)
    assert BOs[0, 6] == pytest.approx(1.1905, abs=1e-4)
    assert BOs[7, 8] == pytest.approx(1.3327, abs=1e-4)
