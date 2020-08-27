
import numpy as np
import pytest

from pysisyphus.calculators import ORCA
from pysisyphus.cos.NEB import NEB
from pysisyphus.helpers import geom_from_library
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.run import run_tsopt_from_cos
from pysisyphus.testing import using


# @pytest.mark.parametrize(
    # "tsopt_key, tsopt_kwargs", [
        # # ("dimer", ),
        # ("rsprfo", {}),
    # ]
# )
# def test_run_tsopt_from_cos(tsopt_key, tsopt_kwargs):
    # from pysisyphus.calculators.AnaPot import AnaPot
    # from pysisyphus.interpolate.helpers import interpolate_all
    # initial = AnaPot.get_geom((-1.05274, 1.02776, 0))
    # final = AnaPot.get_geom((1.94101, 3.85427, 0))
    # geoms = (initial, final)


    # images = interpolate_all((initial, final), 10)
    # for image in images:
        # image.set_calculator(AnaPot())

    # cos = NEB(images)
    # opt_kwargs = {
        # # "thresh": "gau_loose",
        # "max_cycles": 23,
    # }
    # opt = SteepestDescent(cos, **opt_kwargs)
    # # import pdb; pdb.set_trace()
    # opt.run()

    # run_tsopt_from_cos(cos, tsopt_key, tsopt_kwargs)


@pytest.mark.skip
@pytest.mark.parametrize(
    "tsopt_key, tsopt_kwargs", [
        # ("dimer", ),
        pytest.param("rsprfo", {}, marks=using("orca")),
    ]
)
def test_run_tsopt_from_cos(tsopt_key, tsopt_kwargs):
    geoms = geom_from_library("ala_dipeptide_iso_b3lyp_631gd_10_images.trj")

    orca_kwargs = {
        "keywords": "b3lyp 6-31G* rijcosx",
        "pal": 2,
        "mem": 2000,
        "charge": 0,
        "mult": 1,
    }

    calc_number = 0
    def calc_getter():
        nonlocal calc_number
        calc = ORCA(**orca_kwargs, calc_number=calc_number)
        calc_number += 1
        return calc

    for i, geom in enumerate(geoms):
        geom.set_calculator(calc_getter())

    cos = NEB(geoms)
    # We don't actually need to run an optimization step.
    # When calling 'get_hei_index()' in the run_tsopt... function energy
    # calculations will be done.

    # opt_kwargs = {
        # "max_cycles": 1,
    # }
    # opt = SteepestDescent(cos, **opt_kwargs)
    # opt.run()

    tsopt_kwargs_ = {
        "trust_max": 0.3,
        "thresh": "gau_loose",
        "overachieve_factor": 1,
        # "hessian_recalc": 10,
        # "do_hess": False,
    }
    tsopt_kwargs_.update(tsopt_kwargs)

    run_tsopt_from_cos(cos, tsopt_key, tsopt_kwargs_, calc_getter)
