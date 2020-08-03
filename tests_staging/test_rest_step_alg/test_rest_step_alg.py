import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.optimizers.RSA import RSA


def test_rsa():
    calc_cls = AnaPot
    start = (0.667, 1.609, 0.)
    ref_coords = (1.941, 3.8543, 0.)
    geom = calc_cls.get_geom(start)

    print("@Using", calc_cls)

    opt_kwargs = {
        "thresh": "gau_tight",
        "hessian_init": "calc",
    }
    opt = RSA(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged

    ref_coords = np.array(ref_coords)
    diff = ref_coords - geom.coords
    diff_norm = np.linalg.norm(diff)
    print(f"@\tnorm(diff)={diff_norm:.8f}")
    assert diff_norm < 6e-5

    # import matplotlib.pyplot as plt
    # calc = geom.calculator
    # calc.plot()
    # coords = np.array(opt.coords)
    # ax = calc.ax
    # ax.plot(*coords.T[:2], "ro-")
    # plt.show()
