from pprint import pprint

import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.DimerMethod import DimerMethod
from pysisyphus.tsoptimizers.dimer import dimer_method


def test_dimer():
    calc_getter = lambda: AnaPot()
    geom_ref = AnaPot.get_geom((-0.2, 1.1, 0))
    geom = geom_ref.copy()
    geom.set_calculator(calc_getter())

    N_init = np.array((0.3, 0.7, 0.))

    # Reference
    ref_kwargs = {
        "max_cycles": 1,
        "ana_2dpot": True,
        "calc_getter": calc_getter,
        "f_tran_mod": False,
        "N_init": N_init,
    }
    res = dimer_method([geom_ref, ], **ref_kwargs)
    last_cycle = res.dimer_cycles[-1]
    pprint(last_cycle)

    f0 = last_cycle.f0
    f_tran = last_cycle.f_tran

    # New implementation
    dimer_kwargs = {
        "calculator": AnaPot(),
        "geometry": geom,
        "N_init": N_init,
    }
    dimer = DimerMethod(**dimer_kwargs)
    geom.set_calculator(dimer)

    dim_forces = geom.forces
    # import pdb; pdb.set_trace()
    np.testing.assert_allclose(dim_forces, f_tran)

    from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
    opt = PreconLBFGS(geom, precon=False, line_search=None)
    opt.run()

    AnaPot().plot_opt(opt)


    # calc = geom.calculator
    # calc.plot(show=True)
