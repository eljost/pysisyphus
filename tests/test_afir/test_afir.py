import numpy as np
import pytest

from pysisyphus.calculators.AFIR import AFIR
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.calculators import XTB
from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.linalg import finite_difference_hessian
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


init_logging()


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs, ref_cycle, ccl_dist, oc_dist",
    [
        pytest.param(
            PySCF,
            {"basis": "6-31g*", "xc": "b3lyp", "pal": 2},
            28,
            4.794052,
            2.677647,
            marks=using("pyscf"),
        ),
        pytest.param(XTB, {}, 27, 5.244300, 2.6294451, marks=using("xtb")),
    ],
)
def test_ohch3f_anion(calc_cls, calc_kwargs, ref_cycle, ccl_dist, oc_dist):
    """Example (R1) from
        https://aip.scitation.org/doi/pdf/10.1063/1.3457903?class=pdf

    See Fig. 3 and Fig. 4
    """

    geom = geom_loader("lib:ohch3f_anion_cs.xyz")
    # OH group is a fragment
    fragment_indices = [
        (5, 6),
    ]
    gamma = 100 / AU2KJPERMOL
    calc = calc_cls(charge=-1, **calc_kwargs)
    afir = AFIR(calc, fragment_indices, gamma, ignore_hydrogen=True)
    geom.set_calculator(afir)

    opt = RFOptimizer(geom, dump=True, trust_max=0.3)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle

    # Broken C-Cl bond
    c3d = geom.coords3d
    assert np.linalg.norm(c3d[0] - c3d[4]) == pytest.approx(ccl_dist, abs=1e-4)
    # Formed O-C bond
    assert np.linalg.norm(c3d[0] - c3d[5]) == pytest.approx(oc_dist, abs=1e-4)


@using("xtb")
def test_three_frag_afir():
    geom = geom_loader("lib:afir3test.xyz", coord_type="redund")
    fragment_indices = [
        (0, 1, 2),
        (3, 4, 5, 6),
    ]
    calc = XTB()
    gamma = 150 / AU2KJPERMOL
    afir = AFIR(calc, fragment_indices, gamma, ignore_hydrogen=False)
    geom.set_calculator(afir)

    opt = RFOptimizer(geom, dump=True, overachieve_factor=2)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 34
    assert geom.energy == pytest.approx(-22.575014)

    c3d = geom.coords3d
    assert np.linalg.norm(c3d[3] - c3d[9]) == pytest.approx(2.6235122)
    assert np.linalg.norm(c3d[2] - c3d[0]) == pytest.approx(3.8615080)


@using("pyscf")
def test_afir_hessian():
    atoms = ("N", "N")
    coords = np.zeros((2, 3))
    coords[1, 2] = 3.0
    geom = Geometry(atoms, coords)

    calc = PySCF(basis="sto3g")
    afir_calc = AFIR(
        calc,
        fragment_indices=[[0]],
        gamma=100 / AU2KJPERMOL,
    )
    geom.set_calculator(afir_calc)

    afir_hess = geom.hessian

    def grad_func(coords):
        return -geom.get_energy_and_forces_at(coords)["forces"]

    fd_hessian = finite_difference_hessian(geom.coords, grad_func, acc=2)
    np.set_printoptions(suppress=True, precision=4)
    # Locally it works with atol=1e-4, but GH action seem to require 4e-4
    np.testing.assert_allclose(afir_hess, fd_hessian, atol=4e-4)
