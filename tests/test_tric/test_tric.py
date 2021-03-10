import itertools as it

import numpy as np

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
