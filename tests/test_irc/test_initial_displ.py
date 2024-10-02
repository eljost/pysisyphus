import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
import pysisyphus.calculators.Pots1d as pots1d
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.irc.initial_displ import cubic_displ_for_geom
from pysisyphus.irc import EulerPC
from pysisyphus.testing import using
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer


@pytest.fixture
def anapot_ts():
    ts_coords = (0.61173113, 1.49297317, 0.0)
    geom = AnaPot.get_geom(ts_coords)
    return geom


def test_cubic_displ(anapot_ts):
    step_plus, step_minus, _ = cubic_displ_for_geom(anapot_ts)
    ref_step_plus = [-0.016441, -0.011105, 0.0]
    np.testing.assert_allclose(step_plus, ref_step_plus, atol=1e-6)


@pytest.mark.parametrize(
    "displ",
    [
        "energy",
        "energy_cubic",
    ],
)
def test_irc_cubic_displ(displ, anapot_ts):
    irc_kwargs = {
        "displ": displ,
        "displ_energy": 1,
    }
    irc = EulerPC(anapot_ts, **irc_kwargs)
    irc.run()

    # anapot_ts.calculator.plot_irc(irc, show=True)


@using("pyscf")
@pytest.mark.parametrize("displ", ("energy", "length", "energy_cubic"))
def test_hcn_initial_displ(displ):
    geom = geom_loader("lib:hcn_iso_hf_sto3g_ts_opt.xyz")
    geom.set_calculator(PySCF(pal=2, basis="sto3g"))
    irc = EulerPC(geom, displ=displ, max_cycles=1)
    irc.run()


def print_cubic1dpot_ts_coord():
    geom = pots1d.Cubic1dPot().get_saddles(i=0)
    opt = RSPRFOptimizer(geom, thresh="gau_vtight", hessian_recalc=1)
    opt.run()
    print("opt x coord", geom.coords[0])

    irc = EulerPC(geom, displ="energy_cubic", max_cycles=5, displ_energy=0.5)
    irc.run()


@pytest.mark.parametrize("pot_cls", (pots1d.Cubic1dPot, pots1d.Quadratic1dPot))
@pytest.mark.parametrize("energy_diff_ref", (-0.001, -0.0005, -0.0001))
def test_cubic_displ_exact(pot_cls, energy_diff_ref):
    # The two quadratic and cubic 1d potentials can be handled nearly exactly.
    # It is not fully exact, as the third derivatives are estimated via finite
    # differences.
    geom = pot_cls().get_saddles(i=0)

    step_plus, step_minus, _ = cubic_displ_for_geom(geom, dE=energy_diff_ref)

    energy = geom.energy

    for step in (step_plus, step_minus):
        energy_tmp = geom.get_energy_at(geom.coords + step)
        energy_diff = energy_tmp - energy
        # There is a nice agreement between actual and desired energy lowering
        assert energy_diff == pytest.approx(energy_diff_ref, abs=1e-10)
        # print(f"{energy_diff: >12.6f} au, {energy_diff_ref: >12.6f} au")


def plot_quadratic_cubic_displacement():
    import mymplrc
    import matplotlib.pyplot as plt

    from pysisyphus.irc.initial_displ import get_curv_vec

    ts_coords = (0.61173113, 1.49297317, 0.0)
    geom = AnaPot.get_geom(ts_coords)

    step_plus, step_minus, res = cubic_displ_for_geom(geom)

    H = geom.mw_hessian
    w, v = np.linalg.eigh(H)
    w0 = w[0]
    v0 = v[:, 0]
    v1 = get_curv_vec(H, res.G_vec, v0, w0)

    def step(ds):
        return ds * v0 + ds**2 * v1 / 2

    lim = 1.5
    x0 = geom.mw_coords
    quad_steps = np.array([-lim * v0, lim * v0]) + x0
    third_steps = np.array([step(ds) for ds in np.linspace(-lim, lim)]) + x0

    calc = geom.calculator
    calc.plot(figsize=(8, 4.5))
    fig = calc.fig
    ax = calc.ax

    ax.scatter(*ts_coords[:2], s=50, c="red", label="TS", zorder=5)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        "$V(x, y) = 4 + 4.5x - 4y + x^2 + 2y^2-2xy + x^4 - 2x^2 y$",
    )

    i = 0

    def save():
        nonlocal i
        ax.legend()
        fig.savefig(f"initial_displ_{i:03d}.png")
        i += 1

    fig.tight_layout()
    # fig.savefig("init_displ_000.pdf")
    save()
    ax.plot(*quad_steps.T[:2], lw=3, label="2nd derivs.")
    save()
    ax.plot(*third_steps.T[:2], lw=3, label="3rd derivs.")
    save()

    # plt.show()


if __name__ == "__main__":
    plot_quadratic_cubic_displacement()
