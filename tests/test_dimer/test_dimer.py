import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.Dimer import Dimer
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.testing import using


init_logging()


@pytest.mark.parametrize(
    "rotation_method, ref_cycle",
    [
        ("direct", 9),
        ("fourier", 9),
    ]
)
def test_dimer(rotation_method, ref_cycle):
    coords = (-0.2, 1.1, 0)
    geom = Geometry(("X", ), coords)
    N_raw = np.array((0.83, 0.27, 0.))

    # New implementation
    dimer_kwargs = {
        "rotation_method": rotation_method,
        "calculator": AnaPot(),
        "N_raw": N_raw,
        "rotation_remove_trans": False,
    }
    dimer = Dimer(**dimer_kwargs)
    geom.set_calculator(dimer)

    opt_kwargs = {
        "precon": False,
        "line_search": None,
        "max_step_element": 0.25,
        "thresh": "gau_tight",
        "max_cycles": 15,
    }
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()
    # AnaPot().plot_opt(opt)

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
    assert geom.energy == pytest.approx(2.80910484)

    # AnaPot().plot_opt(opt)


@using("pyscf")
@pytest.mark.parametrize(
    "bonds, ref_cycle",
    [
        (None, 9),
        ([[1, 2, -1], [2, 0, 1]], 9),
    ]
)
def test_dimer_hcn(bonds, ref_cycle):
    geom = geom_loader("lib:baker_ts/01_hcn.xyz")
    ref_energy = -92.24604
    N_raw = " 0.5858  0.      0.0543 " \
             "-0.7697 -0.      0.061 " \
             "0.2027  0.     -0.1295".split()
    if bonds is not None:
        N_raw = None

    calc = PySCF("321g", pal=2)

    dimer_kwargs = {
        "rotation_method": "fourier",
        "calculator": calc,
        "N_raw": N_raw,
        "bonds": bonds,
    }
    dimer = Dimer(**dimer_kwargs)
    geom.set_calculator(dimer)

    opt_kwargs = {
        "precon": True,
        "max_step_element": 0.25,
        "max_cycles": 15,
    }
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
    assert geom.energy == pytest.approx(ref_energy)


@pytest.mark.parametrize(
    "bonds",
    [
        [(0, 4, 1)],
        [(0, 4, -1)],
    ]
)
def test_N_raw(bonds):
    geom = geom_loader("lib:baker_ts/08_formyloxyethyl.xyz")

    dimer_kwargs = {
        "calculator": None,
        "bonds": bonds,
    }

    dimer = Dimer(**dimer_kwargs)
    dimer.set_N_raw(geom.coords)
    N = dimer.N.reshape(-1, 3)

    from_, to_, weight = bonds[0]
    ref_row = np.array((0.66128738, 0.24330643, -0.05916908))
    np.testing.assert_allclose(N[from_], weight * ref_row)
    np.testing.assert_allclose(N[to_], -weight * ref_row)


def test_bias_rotation():
    geom = geom_loader("lib:baker_ts/01_hcn.xyz")
    bonds = ((1, 2, -1), (2, 0, 1))
    calc = PySCF("321g", pal=2)

    dimer_kwargs = {
        "rotation_method": "fourier",
        "calculator": calc,
        "bonds": bonds,
        # "rotation_thresh": 2e-4,
        "bias_rotation": True,
    }
    dimer = Dimer(**dimer_kwargs)
    geom.set_calculator(dimer)

    f = geom.forces
    print(f)


def test_add_gaussian():
    geom = AnaPot.get_geom((-0.2, 1.1, 0))
    N_raw = np.array((0.3, 0.7, 0.))

    calc = geom.calculator
    dimer_kwargs = {
        "calculator": calc,
        "N_raw": N_raw,
        "rotation_remove_trans": False,
    }
    dimer = Dimer(**dimer_kwargs)
    geom.set_calculator(dimer)

    g0 = dimer.add_gaussian(geom.atoms, geom.coords, dimer.N)

    assert g0.height == pytest.approx(0.1708984)


@using("pyscf")
@pytest.mark.parametrize(
    "rotation_remove_trans", [
        True,
        False,
    ]
)
def test_remove_translation(rotation_remove_trans):

    name = "01_hcn.xyz with translations removed" if rotation_remove_trans \
           else "01_hcn.xyz without translations removed"
    ref_energy = -92.24604
    geom = geom_loader("lib:baker_ts/01_hcn.xyz")
    downhill_geom = geom_loader("lib:baker_ts/01_hcn_downhill.xyz")
    N_raw = geom.coords - downhill_geom.coords

    calc = PySCF("321g", pal=2)

    dimer_kwargs = {
        "rotation_method": "fourier",
        "calculator": calc,
        "N_raw": N_raw,
        "length": 0.0189,
        "rotation_tol": 5,
        "rotation_remove_trans": rotation_remove_trans,
    }
    dimer = Dimer(**dimer_kwargs)
    geom.set_calculator(dimer)

    opt_kwargs = {
        "precon": True,
        "max_step_element": 0.25,
        "max_cycles": 50,
        "thresh": "baker",
        "c_stab": 0.103,
    }
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(ref_energy)

    print(f"@{name} converged using {dimer.force_evals} force evaluations")

