import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords import (
    Bend,
    LinearBend,
    Stretch,
    Bend,
    Torsion,
)
from pysisyphus.intcoords.PrimTypes import PrimTypes as PT
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


def compare_internals(xyz_fn):
    geom = geom_loader(xyz_fn, coord_type="redund")

    cls = {
        2: Stretch,
        3: Bend,
        4: Torsion,
    }
    for pi in geom.internal._prim_internals:
        inds = pi.inds
        obj = cls[len(inds)](inds)
        print(obj)
        v, g = obj.calculate(geom.coords3d, gradient=True)
        assert v == pytest.approx(pi.val)
        np.testing.assert_allclose(g, pi.grad)


def test_compare_h2o2():
    xyz_fn = "lib:h2o2_hf_321g_opt.xyz"
    compare_internals(xyz_fn)


def test_compare_biaryl():
    xyz_fn = "lib:biaryl_bare_pm6_splined_hei.xyz"
    compare_internals(xyz_fn)


@using("pyscf")
def test_allene_opt():
    geom = geom_loader("lib:08_allene.xyz", coord_type="redund")

    calc = PySCF(basis="321g", pal=1)
    geom.set_calculator(calc)
    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()

    assert opt.is_converged
    # assert opt.cur_cycle == 7
    assert geom.energy == pytest.approx(-115.21991342)


def test_hydrogen_bonds_fragments():
    geom = geom_loader("lib:hydrogen_bond_fragments_test.xyz", coord_type="redund")
    int_ = geom.internal
    hydrogen_bonds = sum(
        [1 for pt, *indices in geom.internal.typed_prims if pt == PT.HYDROGEN_BOND]
    )
    assert hydrogen_bonds == 1
    assert len(int_.fragments) == 3
