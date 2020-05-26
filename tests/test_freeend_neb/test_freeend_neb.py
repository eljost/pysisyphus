from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.FreeEndNEBPot import FreeEndNEBPot
from pysisyphus.cos.FreeEndNEB import FreeEndNEB
from pysisyphus.optimizers.SteepestDescent import SteepestDescent


def test_freeend_anapot():
    calc = AnaPot()
    geoms = calc.get_path(15)
    geoms = geoms[3:11]

    fe_neb = FreeEndNEB(geoms)
    opt = SteepestDescent(fe_neb)
    opt.run()

    # calc.anim_opt(opt, show=True)

    assert opt.is_converged
    assert opt.cur_cycle == 31


def test_freeend_freeendpot():
    calc = FreeEndNEBPot()
    geoms = calc.get_path(50)[:35]

    kwargs = {
        "k_min": 5,
        "k_max": 10,
    }
    fe_neb = FreeEndNEB(geoms, **kwargs)
    opt = SteepestDescent(fe_neb, max_cycles=100)
    opt.run()

    # calc.anim_opt(opt, show=True)

    assert opt.is_converged
    assert opt.cur_cycle == 93
