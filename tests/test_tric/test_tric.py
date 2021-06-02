import itertools as it

import numpy as np
import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.PrimTypes import PrimTypes as PT
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.calculators.XTB import XTB
from pysisyphus.testing import using


def trics_for_frag(inds):
    return (
        (PT.TRANSLATION_X, inds),
        (PT.TRANSLATION_Y, inds),
        (PT.TRANSLATION_Z, inds),
        (PT.ROTATION_A, inds),
        (PT.ROTATION_B, inds),
        (PT.ROTATION_C, inds),
    )


@using("xtb")
def test_two_waters(this_dir):
    inds = ((0, 1, 2), (3, 4, 5))
    print(inds)
    tp0 = trics_for_frag(inds[0])
    tp1 = trics_for_frag(inds[1])
    tps = tp0 + tp1
    print(tps)

    geom = geom_loader(
        this_dir / "2h2o.xyz",
        coord_type="redund",
        coord_kwargs={"typed_prims": tps},
    )
    # geom.jmol()
    calc = XTB(pal=2)
    geom.set_calculator(calc)

    opt_kwargs = {
        "dump": True,
        "max_cycles": 5,
    }
    opt = LBFGS(geom, **opt_kwargs)
    opt.run()


@using("xtb")
def test_21_waters(this_dir):
    waters = 21
    inds = np.arange(waters*3).reshape(-1, 3)
    tps = list(it.chain(*[trics_for_frag(i) for i in inds]))

    geom = geom_loader(
        this_dir / "output.pdb",
        coord_type="redund",
        coord_kwargs={"typed_prims": tps},
    )
    # geom.jmol()
    calc = XTB(gfn="ff", pal=2)
    geom.set_calculator(calc)

    opt_kwargs = {
        "dump": True,
        "max_cycles": 50,
    }
    opt = LBFGS(geom, **opt_kwargs)
    opt.run()

# from geometric.internal import Rotator as geomRot
# from pysisyphus.intcoords.Rotation import RotationA as pysisRot
# from pysisyphus.linalg import get_rot_mat

# @pytest.fixture
# def prep_comp_h2o():
    # ref = geom_loader("start.xyz")
    # rot = geom_loader("rot.xyz")

    # inds = (3, 4, 5)
    # return (ref.coords3d, rot.coords3d, inds)


# @pytest.fixture
# def prep_comp_biaryl():
    # # ref = geom_loader("biaryl.xyz")  # 30
    # # ref = geom_loader("C60-Ih.xyz")  # 60
    # ref = geom_loader("ptphen.xyz")  # 95
    # inds = list(range(len(ref.atoms)))
    # ref3d = ref.coords3d

    # # Create rotated coordinates
    # abc = 0.3, 0.2, 0.1
    # R = get_rot_mat(abc)
    # rot3d = R.dot(ref3d.T).T
    # return ref3d, rot3d, inds


# def test_this(prep_comp_biaryl):
    # ref3d, rot3d, inds = prep_comp_biaryl
    # # gr = geomRot(inds, ref3d)
    # pr = pysisRot(inds, ref_coords3d=ref3d)
    # val, grad = pr.calculate(rot3d, gradient=True)
    # # print(val, grad)


# def test_geom(prep_comp_biaryl):
    # ref3d, rot3d, inds = prep_comp_biaryl
    # gr = geomRot(inds, ref3d)
    # val = gr.value(rot3d)
    # grad = gr.derivative(rot3d)
    # # print(val, grad)
