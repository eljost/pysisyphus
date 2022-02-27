from pathlib import Path

import pytest

from pysisyphus.calculators import EGO, XTB, ORCA
from pysisyphus.helpers import geom_loader, align_geoms
from pysisyphus.init_logging import init_logging
from pysisyphus.testing import using
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


init_logging()


@pytest.mark.skip
@using("xtb")
@pytest.mark.parametrize(
    "fn, cycles", (
        ("gbutene.xyz", 6),
        ("butene.xyz", 6),
        ("tbutene.xyz", 6),
        ("pyranose.xyz", 9)
    )
)
def test_ego(fn, cycles):
    def get_act_calc():
        return ORCA("b3lyp 6-31G*", pal=6)
        # return XTB(pal=6, quiet=True)

    stem = Path(fn).stem
    org_calc = get_act_calc()
    org_geom = geom_loader(fn, coord_type="redund")
    org_geom.set_calculator(org_calc)
    org_opt = RFOptimizer(org_geom, prefix="org", dump=True, thresh="gau")
    org_opt.run()

    cycles = 6
    sgo_geom = org_geom
    egos = list()
    sgos = list()
    for i in range(cycles):
        ego_geom = sgo_geom.copy(coord_type="cart")
        ego_geom.standard_orientation()
        act_calc = get_act_calc()
        calc = EGO(act_calc, ego_geom)
        ego_geom.set_calculator(calc)
        ego_opt = RFOptimizer(ego_geom, max_cycles=150, prefix=f"{i}_ego", dump=True)
        ego_opt.run()
        egos.append(ego_geom.copy())

        sgo_geom = ego_geom.copy(coord_type="redund")
        sgo_geom.set_calculator(act_calc)
        sgo_opt = RFOptimizer(sgo_geom, prefix=f"{i}_sgo", dump=True, thresh="gau_tight")
        sgo_opt.run()
        sgos.append(sgo_geom.copy())

    geoms = [org_geom, ]
    comments = ["org", ]
    for i in range(cycles):
        geoms.append(egos[i])
        comments.append(f"EGO {i}")
        geoms.append(sgos[i])
        comments.append(f"SGO {i}")
    align_geoms(geoms)
    trj = "\n".join([geom.as_xyz(comment=comment) for geom, comment in zip(geoms, comments)])
    with open(f"{stem}_sego.trj", "w") as handle:
        handle.write(trj)

    with open(f"{stem}_sgos.trj", "w") as handle:
        handle.write(
            "\n".join([geom.as_xyz() for geom in sgos])
        )
