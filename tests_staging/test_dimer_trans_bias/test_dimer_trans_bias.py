import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.Dimer import Dimer
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_library
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.testing import using


init_logging()


def test_bias_translation():
    # coords = (-1.05, 1.02, 0)
    coords = (-1.009, 1.024, 0)
    geom = AnaPot.get_geom(coords)
    N_raw = np.array((0.9603, 0.2789, 0.))

    calc = geom.calculator
    # calc.plot_eigenvalue_structure()
    dimer_kwargs = {
        "calculator": calc,
        "N_raw": N_raw,
        "rotation_disable": True,
        # "rotation_method": "direct",
        # "bias_rotation": True,
        "bias_translation": True,
        # "bias_gaussian_dot": 0.5,
    }
    dimer = Dimer(**dimer_kwargs)
    geom.set_calculator(dimer)


    import matplotlib.pyplot as plt
    def plot(dimer, coords0, title=None, step_norms=None):
        if step_norms is None:
            step_norms = list()
            step_norms_max = 0.5
        else:
            step_norms_max = max(step_norms)

        N = dimer.N
        # steps = np.linspace(-0.025, 0.125, 30)
        steps = np.linspace(-0.05, max(step_norms_max, 0.1), 50)
        stepsN = steps[:,None]*N
        Ncoords = coords0 + steps[:,None] * N
        results = [dimer.get_forces(geom.atoms, coords) for coords in Ncoords]
        ens, forces = zip(*[(r["energy"], r["forces"]) for r in results])
        forces = np.array(forces)
        norms = np.linalg.norm(forces, axis=1)
        step1_norm = np.linalg.norm(np.array((-0.759042, 1.09659, 0.) - coords0))
        # fig, (ax0, ax1) = plt.subplots(nrows=2)
        fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
        ax0.plot(steps, ens)
        ax0.set_title("energy")
        ax1.plot(steps, norms, "o-")
        ax1.set_title("norm(forces)")
        # ax2.plot(*forces.T[:2], "o-")
        # for i, f in enumerate(forces):
            # ax2.annotate(i, f[:2])

        # ax2.quiver(Ncoords.[::2]
        # ax2.quiver(Ncoords[:,0], Ncoords[:,1], forces[:,0], forces[:,1])
        ax2.quiver(stepsN[:,0], stepsN[:,1], forces[:,0], forces[:,1], scale=10)
        ax2.axvline(0.0529)

        for ax in (ax0, ax1):
            ax.axvline(0, color="r", ls="--")

        for i, x in enumerate(step_norms):
            for ax in (ax0, ax1):
                ax.axvline(x, color="k", ls="--")
                mi, ma = ax.get_ylim()
                y = (mi+ma)/2
                ax.annotate(i, (x, y))
        if title:
            fig.suptitle(title)
        return fig

    coords0 = geom.coords.copy()
    # f_ub = plot(dimer, coords0, "unbiased")

    g0 = dimer.add_gaussian(geom.atoms, geom.coords, dimer.N)
    print(g0)
    # f_b = plot(dimer, coords0, "biased")
    plt.show()
    # return

    opt_kwargs = {
        "precon": False,
        # "max_step_element": 0.06,
        # "max_step_element": 0.03,
        "max_step_element": 0.0529,
        # "max_cycles": 9,
        "max_cycles": 10,
        "dump": True,
    }
    # opt = PreconLBFGS(geom, **opt_kwargs)
    from pysisyphus.optimizers.PreconSteepestDescent import PreconSteepestDescent
    opt = PreconSteepestDescent(geom, **opt_kwargs)
    opt.run()


    step_norms = np.linalg.norm([c - coords0 for c in opt.coords], axis=1)
    f_b = plot(dimer, coords0, "biased", step_norms)
    print(step_norms)
    for i, step in enumerate(opt.steps):
        print(i, step)
    # calc.plot_opt(opt)
    plt.show()
