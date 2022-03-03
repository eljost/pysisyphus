import numpy as np
import pytest

from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.helpers import geom_loader
from pysisyphus.calculators import XTB
from pysisyphus.testing import using


@using("xtb")
@pytest.mark.parametrize(
    "coord_type",
    (
        "cartesian",
        "mwcartesian",
    ),
)
def test_frozen_cartesian(coord_type):
    freeze_atoms = [0, 1, 2, 3]
    geom = geom_loader(
        "lib:benzene.xyz", freeze_atoms=freeze_atoms, coord_type=coord_type
    )
    coords3d_0 = geom.coords3d.copy()
    dof = 3 * len(geom.atoms) - 3 * len(freeze_atoms)
    assert geom.coords.size == 3 * len(geom.atoms) - 3 * len(freeze_atoms)

    geom.set_calculator(XTB(pal=2))
    opt = LBFGS(geom, thresh="gau", dump=True)
    opt.run()

    np.testing.assert_allclose(geom.coords3d[freeze_atoms], coords3d_0[freeze_atoms])
    H = geom.hessian
    assert H.shape == (dof, dof)
