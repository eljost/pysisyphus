# Tested against Turbomole 7.4.1

import numpy as np
import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators import Turbomole
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.testing import using


@pytest.fixture
def geom():
    return geom_loader("lib:h2o_bp86_def2svp_opt.xyz")


@using("turbomole")
def test_turbomole_hessian(geom, this_dir):
    turbo_kwargs = {
        "control_path": this_dir / "./control_path_dft_gs",
    }
    calc = Turbomole(**turbo_kwargs)
    geom.set_calculator(calc)

    mw_H = geom.mw_hessian
    proj_H = geom.eckart_projection(mw_H)
    w, v = np.linalg.eigh(proj_H)
    nus = eigval_to_wavenumber(w)
    print("nus / cm⁻¹:", nus)

    # Turbomole reference values
    # turbo_ref_nus = np.array((1607.81, 3684.62, 3783.64))
    ref_nus = np.array((1607.984768, 3684.472251, 3783.356437))
    np.testing.assert_allclose(nus[-3:], ref_nus, atol=1e-2)


def test_turbomole_big_hessian_parsing(this_dir):
    geom = geom_loader("lib:ch4_12.xyz")
    control_path = this_dir / "control_path_big_hess"
    turbo_kwargs = {
        "control_path": control_path,
    }
    calc = Turbomole(**turbo_kwargs)
    geom.set_calculator(calc)
    results = calc.parse_hessian(path=control_path)

    hessian = results["hessian"]
    assert hessian.size == 180 ** 2


@using("turbomole")
@pytest.mark.parametrize(
    "control_path, ref_energy",
    [
        # Ground state
        ("./control_path_dft_gs", -76.36357867674),
        # Excited state
        ("./control_path_dft_es1", -76.0926146085),
        # ricc2
        ("./control_path_ricc2", -75.8716368247),
    ],
)
def test_h2o_energy(control_path, ref_energy, geom, this_dir):
    turbo_kwargs = {
        "control_path": this_dir / control_path,
    }
    calc = Turbomole(**turbo_kwargs)
    geom.set_calculator(calc)

    energy = geom.energy

    assert energy == pytest.approx(ref_energy)


@using("turbomole")
@pytest.mark.parametrize(
    "control_path, ref_energy, ref_force_norm",
    [
        # Ground state
        ("./control_path_dft_gs", -76.36357867674, 1.30342385e-5),
        # Excited state gradient, TDDFT
        ("./control_path_dft_es1", -76.0926146085, 0.16006233),
        # Excited state gradient, ricc2
        ("./control_path_ricc2", -75.8716368247, 0.15925937),
    ],
)
def test_h2o_forces(control_path, ref_energy, ref_force_norm, geom, this_dir):
    turbo_kwargs = {
        "control_path": this_dir / control_path,
    }
    calc = Turbomole(**turbo_kwargs)
    geom.set_calculator(calc)

    forces = geom.forces
    energy = geom.energy

    norm = np.linalg.norm(forces)

    assert norm == pytest.approx(ref_force_norm, abs=1e-4)
    assert energy == pytest.approx(ref_energy)


@pytest.mark.skip
@using("turbomole")
def test_turbomole_cos(this_dir):
    def calc_getter(charge, mult):
        calc_kwargs = {
            "charge": charge,
            "mult": mult,
            "control_path": this_dir / "control_cos",
            "pal": 2,
        }
        return Turbomole(**calc_kwargs)

    def gs_calc_getter():
        return calc_getter(charge=0, mult=1)

    bench = Benchmark("xtb_rx", calc_getter=calc_getter)
    geoms = bench.get_geoms(11, set_calculator=True)
    for i, geom in enumerate(geoms):
        en = geom.energy
        print(f"{i:02d}: {en:.6f} au")
    start, _, end = geoms
    images = (start, end)
    cos_kwargs = {
        "calc_getter": gs_calc_getter,
        "max_nodes": 9,
        "climb": True,
    }
    cos = GrowingString(images, calc_getter=gs_calc_getter, max_nodes=9, climb=True)
    opt_kwargs = {
        "rms_force": 0.002,
        "rms_force_only": True,
        "dump": True,
    }
    opt = StringOptimizer(cos, **opt_kwargs)
    opt.run()
    assert opt.is_converged

    ens = [image.energy for image in cos.images]
    assert max(ens) == pytest.approx(-178.8393291373)
    assert opt.cur_cycle == 9
