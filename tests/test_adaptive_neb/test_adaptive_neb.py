from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.cos.AdaptiveNEB import AdaptiveNEB
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


def test_anapot_aneb():
    calc = AnaPot()
    all_geoms = calc.get_path(10)
    aneb = AdaptiveNEB(all_geoms)

    opt_kwargs = {
        "rms_force": 3e-3,
    }
    opt = ConjugateGradient(aneb, **opt_kwargs)
    opt.run()

    # ap = calc.anim_opt(opt, show=True)

    assert opt.is_converged
    assert opt.cur_cycle == 21


@using("pyscf")
def test_hcn_aneb():
    run_dict = {
        "preopt": {
            "max_cycles": 3,
        },
        "interpol": {
            "type": "idpp",
            "between": 1,
        },
        "cos": {
            "type": "aneb",
        },
        "opt": {
            "type": "qm",
            "max_cycles": 25,
            "align": True,
        },
        "calc": {
            "type": "pyscf",
            "pal": 2,
            "basis": "321g",
        },
        "geom": {
            "type": "cart",
            "fn": ["lib:hcn.xyz", "lib:hcn_iso_ts.xyz", "lib:nhc.xyz"],
        }
    }
    results = run_from_dict(run_dict)

    assert results.cos_opt.is_converged
    assert results.cos_opt.cur_cycle == 18
    assert results.cos.level == 3
