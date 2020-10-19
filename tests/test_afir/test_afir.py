import numpy as np
import pytest

from pysisyphus.calculators.AFIR import AFIR
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.testing import using


@using("pyscf")
def test_afir():
    """Example (R1) from
        https://aip.scitation.org/doi/pdf/10.1063/1.3457903?class=pdf

    See Fig. 2 and Fig. 4
    """

    geom = geom_loader("lib:ohch3f_anion_cs.xyz")
    # OH group is a fragment
    fragment_indices = [(5, 6), ]
    calc = PySCF(basis="6-31g*", xc="b3lyp", pal=2, charge=-1)
    gamma = 100
    afir = AFIR(calc, fragment_indices, gamma)
    geom.set_calculator(afir)

    opt = RFOptimizer(geom, dump=True, trust_max=.3)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 47

    # Broken C-Cl bond
    c3d = geom.coords3d
    assert np.linalg.norm(c3d[0]-c3d[4]) == pytest.approx(4.805665, abs=1e-4)
    # Formed O-C bond
    assert np.linalg.norm(c3d[0]-c3d[5]) == pytest.approx(2.674330, abs=1e-4)
