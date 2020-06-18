import h5py
import numpy as np
import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.io.hessian import save_hessian
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


@pytest.fixture
def geom():
    fn = "lib:h2o.xyz"
    geom = geom_loader(fn, coord_type="redund")
    calc = PySCF(basis="def2svp", xc="bp86", pal=2)
    geom.set_calculator(calc)
    return geom


@using("pyscf")
def test_save_hessian(geom):
    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()

    hess_fn = "_h2o_hessian.h5"
    save_hessian(hess_fn, geom)

    with h5py.File(hess_fn, "r") as handle:
        hessian = handle["hessian"][:]

    np.testing.assert_allclose(hessian, geom.cart_hessian)
