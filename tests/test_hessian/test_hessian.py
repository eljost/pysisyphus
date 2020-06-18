import h5py
import numpy as np

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.io.hessian import save_hessian
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


@using("pyscf")
def test_save_hessian():
    # fn, hess_fn = "lib:hcfclbr_xtbopt.xyz", "chfclbr.h5"
    fn, hess_fn = "lib:h2o.xyz", "h2o_hessian.h5"
    geom = geom_loader(fn, coord_type="redund")

    calc = PySCF(basis="def2svp", xc="bp86", pal=2)
    geom.set_calculator(calc)

    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()

    save_hessian(hess_fn, geom)

    with h5py.File(hess_fn, "r") as handle:
        hessian = handle["hessian"][:]

    np.testing.assert_allclose(hessian, geom.cart_hessian)
