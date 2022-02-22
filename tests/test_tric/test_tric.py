import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using
from pysisyphus.calculators.XTB import XTB
from pysisyphus.calculators.PySCF import PySCF


@pytest.mark.skip_ci
@using("xtb")
@pytest.mark.parametrize("recalc_B", [True, False])
def test_tric_solvated_sulfate(recalc_B):
    coord_kwargs = {
        # Whether the B-matrix is repeatedly recalculated in the internal
        # Cartesian backtransformation
        "recalc_B": recalc_B,
    }
    geom = geom_loader(
        "lib:SulfateInSolution.xyz", coord_type="tric", coord_kwargs=coord_kwargs
    )
    geom.set_calculator(XTB(pal=4, charge=-2, gfn=2, quiet=True))
    opt = RFOptimizer(
        geom,
        hessian_init="simple",
        max_cycles=500,
        trust_radius=0.3,
        adapt_step_func=True,
        # dump=True,
    )
    opt.run()

    assert opt.is_converged
    assert geom.energy <= -167.54


@using("pyscf")
def test_tric_two_waters():
    geom = geom_loader("lib:water_dimer_oniomopt_test.pdb", coord_type="tric")
    # for i, tp in enumerate(geom.internal.typed_prims):
    # print(i, tp)
    calc = PySCF(basis="sto3g", pal=2, verbose=0)
    geom.set_calculator(calc)
    opt = RFOptimizer(geom, dump=True)
    opt.run()

    assert opt.is_converged
