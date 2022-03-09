#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.calculators import Gaussian16, Turbomole
from pysisyphus.constants import BOHR2ANG
from pysisyphus.helpers import geom_loader
from pysisyphus.irc import *
from pysisyphus.testing import using


@pytest.fixture
def this_dir(request):
    return Path(request.module.__file__).parents[0]


def assert_anapot_irc(irc):
    fc = irc.all_coords[0]
    bc = irc.all_coords[-1]
    forward_ref = np.array((-1.0527, 1.0278,  0.))
    backward_ref = np.array((1.941, 3.8543, 0.))
    forward_diff = np.linalg.norm(fc - forward_ref)
    backward_diff = np.linalg.norm(bc - backward_ref)
    assert forward_diff == pytest.approx(0.05, abs=0.1)
    assert backward_diff == pytest.approx(0.05, abs=0.1)


def plot_irc(irc, title=None):
    geom = irc.geometry
    calc = geom.calculator
    levels = np.linspace(-3, 4, 120)
    calc.plot()
    ax = calc.ax
    ax.plot(*irc.all_coords.T[:2], "ro-")
    if title:
        ax.set_title(title)
    plt.show()


@pytest.mark.parametrize(
    "irc_cls, mod_kwargs, ref", [
        (DampedVelocityVerlet, {"v0": 0.1, "max_cycles": 400,}, None),
        (Euler, {"step_length": 0.05,}, None),
        (EulerPC, {}, None),
        (GonzalezSchlegel, {}, None),
        (IMKMod, {}, None),
        (RK4, {}, None),
        (LQA, {}, None),
    ]
)
def test_anapot_irc(irc_cls, mod_kwargs, ref):
    geom = AnaPot().get_geom((0.61173, 1.49297, 0.))

    kwargs = {
        "step_length": 0.1,
        "rms_grad_thresh": 1e-2,
    }
    kwargs.update(**mod_kwargs)

    irc = irc_cls(geom, **kwargs)
    irc.run()

    # geom.calculator.plot_irc(irc, show=True)
    assert_anapot_irc(irc)


@pytest.mark.parametrize(
    "step_length", [
        (0.1),
        (0.2),
        (0.3),
        (0.4),
    ]
)
def test_imk(step_length):
    geom = AnaPot().get_geom((0.61173, 1.49297, 0.))

    irc_kwargs = {
        "step_length": step_length,
        "rms_grad_thresh": 1e-2,
        "corr_first": True,
        "corr_first_energy": True,
        "corr_bisec_size": 0.0025,
        "corr_bisec_energy": True,
    }

    irc = IMKMod(geom, **irc_kwargs)
    irc.run()

    # plot_irc(irc, irc.__class__.__name__)
    assert_anapot_irc(irc)


@pytest.mark.parametrize(
    "calc_cls, kwargs_", [
        pytest.param(PySCF,
            {"basis": "321g", },
            marks=using("pyscf")
        ),
        pytest.param(Gaussian16,
            {"route": "HF/3-21G"},
            marks=using("gaussian16")
        ),
        pytest.param(Turbomole,
            {"control_path": "./hf_abstr_control_path", "pal": 1},
            marks=using("turbomole")
        ),
    ]
)
def test_hf_abstraction_dvv(calc_cls, kwargs_, this_dir):
    geom = geom_loader("lib:hfabstraction_hf321g_displ_forward.xyz")

    calc_kwargs = {
        "pal": 2,
    }
    calc_kwargs.update(kwargs_)

    if "control_path" in calc_kwargs:
        calc_kwargs["control_path"] = this_dir / calc_kwargs["control_path"]

    print("Using", calc_cls)
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    irc_kwargs = {
        "dt0": 0.5,
        "v0": 0.04,
        "downhill": True,
        "max_cycles": 150,
    }
    dvv = DampedVelocityVerlet(geom, **irc_kwargs)
    dvv.run()

    c3d = geom.coords3d * BOHR2ANG
    def bond(i,j): return np.linalg.norm(c3d[i]-c3d[j])

    assert bond(2, 7) == pytest.approx(0.93, abs=0.01)
    assert bond(4, 7) == pytest.approx(2.42, abs=0.01)
    assert bond(2, 0) == pytest.approx(2.23, abs=0.01)


@using("pyscf")
@pytest.mark.parametrize(
    "irc_cls, irc_kwargs, fw_cycle, bw_cycle",
    [
        (EulerPC, {"hessian_recalc": 10, "dump_dwi": False,}, 30, 38),
        # (EulerPC, {"hessian_recalc": 10, "corr_func": "scipy",}, 19, 23),
    ]
)
def test_hcn_irc(irc_cls, irc_kwargs, fw_cycle, bw_cycle):
    geom = geom_loader("lib:hcn_iso_hf_sto3g_ts_opt.xyz")

    calc = PySCF(
            basis="sto3g",
    )
    geom.set_calculator(calc)

    irc = irc_cls(geom, **irc_kwargs, rms_grad_thresh=1e-4)
    irc.run()

    # approx. +- 0.5 kJ/mol
    ref_energies = [pytest.approx(en) for en in (-91.6444238, -91.67520895)]
    assert irc.forward_energies[0] in ref_energies
    assert irc.backward_energies[-1] in ref_energies


@pytest.mark.parametrize(
    "scipy_method",
    [
        (None),
        ("RK45"),
        ("DOP853"),
    ]
)
def test_eulerpc_scipy(scipy_method):
    geom = AnaPot().get_geom((0.61173, 1.49297, 0.))

    kwargs = {
        "step_length": 0.2,
        "rms_grad_thresh": 1e-2,
        "corr_func": "scipy",
        "scipy_method": scipy_method,
    }

    irc = EulerPC(geom, **kwargs)
    irc.run()

    # plot_irc(irc, irc.__class__.__name__)
    assert_anapot_irc(irc)


@using("pyscf")
@pytest.mark.parametrize(
    "hessian_init, ref_cycle", [
        ("calc", 28),
        pytest.param("fischer", 0, marks=pytest.mark.xfail),
        pytest.param("unit", 0, marks=pytest.mark.xfail),
        ("lindh", 28),
        ("simple", 28),
        ("swart", 28),
    ]
)
def test_downhill_irc_model_hessian(hessian_init, ref_cycle):
    geom = geom_loader("lib:hcn_downhill_model_hessian.xyz")

    calc = PySCF(basis="sto3g", pal=2)
    geom.set_calculator(calc)

    irc_kwargs = {
        "hessian_init": hessian_init,
        "rms_grad_thresh": 5e-3,
        "downhill": True,
    }

    irc = EulerPC(geom, **irc_kwargs)
    irc.run()

    assert irc.downhill_energies[-1] == pytest.approx(-91.67517096968325)
    assert irc.downhill_cycle == ref_cycle


# @pytest.mark.skip()
@pytest.mark.parametrize(
    "step_length", [
        0.1,
        0.2,
        0.3,
        # 0.4  # requires hessian_recalc=1
    ]
)
def test_mb_gs2(step_length):
    calc = MullerBrownPot()
    geom = calc.get_saddles(i=0)

    irc_kwargs = {
        "step_length": step_length,
        "line_search": False,
        # "hessian_recalc": 1,
    }
    irc = GonzalezSchlegel(geom, **irc_kwargs)
    irc.run()
    # calc.plot_irc(irc, show=True, title=f"length {step_length:.2f}")

    assert irc.forward_is_converged
    assert irc.backward_is_converged


@using("pyscf")
@pytest.mark.parametrize(
    "step_length", [
        0.1,
        0.2,
        0.3,
        # 0.4,  # sometimes fails in the CI
    ]
)
def test_hcn_iso_gs2(step_length):
    geom = geom_loader("lib:hcn_iso_hf_sto3g_ts_opt.xyz")
    calc = PySCF(basis="sto3g", verbose=0)
    geom.set_calculator(calc)
    irc_kwargs = {
        "step_length": step_length,
        "displ_energy": 0.0005,
    }
    irc = GonzalezSchlegel(geom, **irc_kwargs)
    irc.run()

    assert irc.forward_is_converged
    assert irc.backward_is_converged


@pytest.mark.parametrize(
    "step_length", [
        0.1,
        0.2,
        # 0.3,
        # 0.4,
    ]
)
def test_mb_eulerpc(step_length):
    calc = MullerBrownPot()
    geom = calc.get_saddles(0)

    irc_kwargs = {
        "step_length": step_length,
        # Using Scipy here takes forever...
        # "corr_func": "scipy",
        # "scipy_method": "BDF",
    }
    irc = EulerPC(geom, **irc_kwargs)
    irc.run()
    # calc.plot_irc(irc, show=True, title=f"length {step_length:.2f}")

    forward_coords = irc.all_coords[0]
    backward_coords = irc.all_coords[-1]
    assert np.linalg.norm(forward_coords - (-0.558, 1.441, 0.0)) <= 2e-2
    assert np.linalg.norm(backward_coords - (-0.050, 0.466, 0.0)) <= 5e-3
