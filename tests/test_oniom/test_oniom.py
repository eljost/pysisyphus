import numpy as np
import pytest

from pysisyphus.calculators.ONIOMv2 import ONIOM
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import do_final_hessian, geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


init_logging()


@using("gaussian16")
def test_energy():
    geom = geom_loader("lib:alkyl17_sto3g_opt.xyz")

    real = set(range(len(geom.atoms)))
    medmin = set((0, 1, 2, 3, 4, 5, 6, 46, 47, 48, 49, 50, 51, 52))
    med = list(real - medmin)
    h1 = list(range(13, 22))
    h2 = list(range(31, 40))

    calcs = {
        "real": {
            "route": "HF/STO-3G",
        },
        "medium": {
            "route": "HF/3-21G",
        },
        "high1": {
            "route": "HF/6-31G",
        },
        "high2": {
            "route": "HF/6-311G",
        },
    }
    for key, calc in calcs.items():
        calc["type"] = "g16"
        calc["pal"] = 2
        calc["mult"] = 1
        calc["charge"] = 0

    models = {
        "med": {
            # "inds": (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
            "inds": med,
            "calc": "medium",
        },
        "h1": {
            # "inds": (4, 5, 6),
            "inds": h1,
            "calc": "high1",
        },
        "h2": {
            # "inds": (10, 11, 12),
            "inds": h2,
            "calc": "high2",
        },
    }

    layers = ["med", ["h1", "h2"]]

    oniom = ONIOM(calcs, models, geom, layers)

    assert oniom.layer_num == 3

    geom.set_calculator(oniom)

    assert geom.energy == pytest.approx(-661.3512410069466)


@pytest.mark.parametrize(
    "calcs, ref_energy, ref_force_norm",
    [
        # From https://doi.org/10.1016/S0166-1280(98)00475-8
        pytest.param(
            {
                "real": {"type": "g16", "route": "hf sto-3g"},
                "high": {"type": "g16", "route": "b3lyp d95v"},
            },
            -153.07432042299052,
            0.03768246934785125,
            marks=using("gaussian16"),
        ),
        # The following two tests should yield identical results
        pytest.param(
            {
                "real": {"type": "g16", "route": "hf sto-3g"},
                "high": {"type": "g16", "route": "b3lyp 3-21g"},
            },
            -152.4529060634755,
            0.018462670668992546,
            marks=using("gaussian16"),
        ),
        pytest.param(
            {
                "real": {"type": "pyscf", "basis": "sto3g"},
                "high": {"type": "pyscf", "xc": "b3lypg", "basis": "321g"},
            },
            -152.4529060634755,
            0.01839279960703439,
            marks=using("pyscf"),
        ),
    ],
)
def test_gradient(calcs, ref_energy, ref_force_norm):
    geom = geom_loader("lib:acetaldehyd_oniom.xyz", coord_type="redund")

    high = [4, 5, 6]

    for key, calc in calcs.items():
        calc["pal"] = 2
        calc["mult"] = 1
        calc["charge"] = 0

    models = {
        "high": {
            "inds": high,
            "calc": "high",
        },
    }

    # layers = ["high"]
    # No layers specified
    layers = None

    oniom = ONIOM(calcs, models, geom, layers)

    assert oniom.layer_num == 2

    geom.set_calculator(oniom)

    # Calculate forces
    forces = geom.forces
    energy = geom.energy

    assert np.linalg.norm(forces) == pytest.approx(ref_force_norm)
    assert energy == pytest.approx(ref_energy)


@pytest.mark.parametrize(
    "calc_key, embedding, ref_energy, ref_force_norm",
    [
        # No embedding
        pytest.param("g16", None, -582.392035, 0.085568849, marks=using("gaussian16")),
        pytest.param("pyscf", None, -582.392035, 0.078387703, marks=using("pyscf")),
        # Electronic embedding
        pytest.param(
            "g16",
            "electronic",
            -582.3997769406087,
            0.08582761,
            marks=using("gaussian16"),
        ),
        pytest.param(
            "pyscf", "electronic", -582.3997769406087, 0.07861744, marks=using("pyscf")
        ),
    ],
)
def test_electronic_embedding(calc_key, embedding, ref_energy, ref_force_norm):
    geom = geom_loader(
        "lib:oniom_ee_model_system.xyz",
        coord_type="redund",
        coord_kwargs={
            "hbond_angles": True,
        },
    )

    all_ = set(range(len(geom.atoms)))
    high = list(sorted(all_ - set((21, 20, 19, 15, 14, 13))))

    calcs_dict = {
        "g16": (
            {
                "real": {"type": "g16", "route": "hf sto-3g"},
                "high": {"type": "g16", "route": "hf 3-21g"},
            }
        ),
        "pyscf": (
            {
                "real": {
                    "type": "pyscf",
                    "basis": "sto3g",
                },
                "high": {"type": "pyscf", "basis": "321g"},
            }
        ),
    }
    calcs = calcs_dict[calc_key]

    for key, calc in calcs.items():
        calc["pal"] = 2
        calc["mult"] = 1
        calc["charge"] = 0

    models = {
        "high": {
            "inds": high,
            "calc": "high",
        },
    }

    oniom_kwargs = {
        "embedding": embedding,
    }
    oniom = ONIOM(calcs, models, geom, **oniom_kwargs)
    geom.set_calculator(oniom)

    # Calculate forces and energy
    forces = geom.forces
    energy = geom.energy

    assert energy == pytest.approx(ref_energy)
    assert np.linalg.norm(forces) == pytest.approx(ref_force_norm)


@using("gaussian16")
def test_oniom_13_coupling():
    geom = geom_loader("lib:oniom_13_coupling_example.xyz")

    calcs = {
        "real": {
            "route": "HF/STO-3G",
        },
        "mid": {
            "route": "HF/3-21G",
        },
        "high": {
            "route": "PM6",
        },
    }
    for key, calc in calcs.items():
        calc["type"] = "g16"
        calc["pal"] = 2
        calc["mult"] = 1
        calc["charge"] = 0

    # "real" must not be defined; will be automatically set
    models = {
        "mid": {
            "inds": range(4, 14),
            "calc": "mid",
        },
        "high": {
            "inds": range(4, 10),
            "calc": "high",
        },
    }

    # "real" must not be defined; will be automatically set
    layers = ["mid", "high"]

    oniom = ONIOM(calcs, models, geom, layers)

    assert oniom.layer_num == 3

    geom.set_calculator(oniom)

    assert geom.energy == pytest.approx(-77.420587)


@using("gaussian16")
def test_acetaldehyde_opt():
    """
    From https://doi.org/10.1016/S0166-1280(98)00475-8
    """
    geom = geom_loader("lib:acetaldehyd_oniom.xyz", coord_type="redund")

    high = [4, 5, 6]

    calcs = {
        "real": {"type": "g16", "route": "hf sto-3g"},
        "high": {"type": "g16", "route": "b3lyp d95v"},
    }

    for key, calc in calcs.items():
        calc["pal"] = 2
        calc["mult"] = 1
        calc["charge"] = 0

    models = {
        "high": {
            "inds": high,
            "calc": "high",
        },
    }

    oniom = ONIOM(calcs, models, geom, layers=None)
    geom.set_calculator(oniom)

    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()
    assert geom.energy == pytest.approx(-153.07526171)

    hess_result = do_final_hessian(geom, save_hessian=False)
    nus = hess_result.nus
    print("wavenumbers / cm⁻¹:", nus)

    assert nus[-1] == pytest.approx(3759.5872)
    assert nus[-3] == pytest.approx(3560.9944)


@using("gaussian16")
def test_yaml_oniom():

    run_dict = {
        "geom": {
            "type": "redund",
            "fn": "lib:acetaldehyd_oniom.xyz",
        },
        "calc": {
            "type": "oniom",
            "calcs": {
                "real": {
                    "type": "g16",
                    "route": "hf sto-3g",
                    "pal": 2,
                },
                "high": {
                    "type": "g16",
                    "route": "b3lyp d95v",
                    "pal": 2,
                },
            },
            "models": {
                "high": {
                    "inds": [4, 5, 6],
                    "calc": "high",
                }
            },
        },
        "opt": {
            "thresh": "gau_tight",
        },
    }
    res = run_from_dict(run_dict)

    opt = res.opt
    assert opt.is_converged
    assert opt.cur_cycle == 7
    assert res.opt_geom.energy == pytest.approx(-153.07526171)


@pytest.mark.skip_ci
@using("pyscf")
def test_oniom3():
    run_dict = {
        "geom": {
            "type": "redund",
            "fn": "lib:oniom3alkyl.pdb",
        },
        "calc": {
            "type": "oniom",
            "calcs": {
                "real": {
                    "type": "pyscf",
                    "basis": "sto3g",
                    "pal": 2,
                },
                "mid": {
                    "type": "pyscf",
                    "basis": "321g",
                    "pal": 2,
                },
                "high": {
                    "type": "pyscf",
                    # "basis": "sto3g",
                    "basis": "431g",
                    "pal": 2,
                },
            },
            "models": {
                "high": {
                    "inds": list(range(7, 16)),
                    "calc": "high",
                },
                "mid": {
                    "inds": list(range(4, 19)),
                    "calc": "mid",
                },
            },
        },
        "opt": {
            "thresh": "gau_tight",
        },
    }
    res = run_from_dict(run_dict)
    print()

    opt = res.opt
    assert opt.is_converged
    assert opt.cur_cycle == 7

    geom = res.opt_geom
    res = do_final_hessian(geom, save_hessian=False)
    nus = res.nus
    np.testing.assert_allclose(nus[[-1, -5]], (3750.1537948, 3566.366994), atol=1e-2)


@pytest.mark.skip
@using("pyscf")
def test_oniomopt_water_dimer():
    run_dict = {
        "xyz": "lib:water_dimer_oniomopt_test.pdb",
        "coord_type": "cart",
        "calc": {
            # "type": "pyscf",
            # "basis": "sto-3g",
            "type": "oniom",
            "calcs": {
                "real": {
                    "type": "pyscf",
                    "basis": "sto3g",
                    "pal": 2,
                },
                "high": {
                    "type": "pyscf",
                    "basis": "321g",
                    "pal": 2,
                },
            },
            "models": {
                "high": {
                    "inds": [0, 1, 2],
                    "calc": "high",
                },
            },
        },
        "opt": {
            "micro_cycles": [3, 1],
            "type": "oniom",
            # "rms_force": 0.0025,
            "rms_force": 0.005,
            "max_cycles": 10,
        },
    }
    res = run_from_dict(run_dict)
    print()

    # opt = res.opt
    # assert opt.is_converged
    # assert opt.cur_cycle == 13


@pytest.mark.skip
@using("pyscf")
def test_oniom_microiters():
    run_dict = {
        "xyz": "lib:oniom_microiters_test.pdb",
        # "xyz": "lib:acetaldehyd_oniom.xyz",
        "coord_type": "cart",
        "calc": {
            # "type": "pyscf",
            # "basis": "sto-3g",
            "type": "oniom",
            "calcs": {
                "real": {
                    "type": "pyscf",
                    "basis": "sto3g",
                    "pal": 2,
                },
                "high": {
                    "type": "pyscf",
                    "basis": "321g",
                    "pal": 2,
                },
            },
            "models": {
                "high": {
                    "inds": [10, 11, 12, 13, 14],
                    # "inds": [4, 5, 6],
                    "calc": "high",
                },
            },
        },
        "opt": {
            "micro_cycles": [3, 1],
            "type": "oniom",
            # "rms_force": 0.0025,
            # "rms_force": 0.005,
            # "max_cycles": 7,
            # "max_cycles": 3,
        },
    }
    res = run_from_dict(run_dict)
    print()

    # opt = res.opt
    # assert opt.is_converged
    # assert opt.cur_cycle == 13


@using("orca")
@pytest.mark.parametrize(
    "fn, mult, high_inds, ref_energy",
    [
        (
            "lib:subst_effect/toluene_b3lypG_631g.xyz",
            1,
            (0, 7, 8, 9),
            -271.470945476921,
        ),
        (
            "lib:subst_effect/toluene_minus_H_b3lypG_631g.xyz",
            2,
            (0, 7, 8),
            -270.824806805671,
        ),
    ],
)
def test_composite_oniom(fn, mult, high_inds, ref_energy):
    charge = 0
    pal = 6
    run_dict = {
        "geom": {
            "type": "cart",
            "fn": fn,
        },
        "calc": {
            "type": "oniom",
            "calcs": {
                "real": {
                    "type": "orca",
                    "keywords": "b3lyp_G 6-31G(d) tightscf",
                    "charge": charge,
                    "mult": mult,
                    "pal": pal,
                },
                "high": {
                    "type": "composite",
                    "calcs": {
                        "ccsdt": {
                            "type": "orca",
                            "keywords": "ccsd(t) 6-31G(d) tightscf",
                        },
                        "mp2_high": {
                            "type": "orca",
                            "keywords": "mp2 6-311+G(2df,2p) tightscf",
                        },
                        "mp2_low": {
                            "type": "orca",
                            "keywords": "mp2 6-31G(d) tightscf",
                        },
                    },
                    "charge": charge,
                    "mult": mult,
                    "pal": pal,
                    "final": "ccsdt + mp2_high - mp2_low",
                },
            },
            "models": {
                "high": {
                    "inds": high_inds,
                    "calc": "high",
                },
            },
        },
    }
    results = run_from_dict(run_dict)
    geom = results.calced_geoms[0]
    en = geom._energy
    print(f"{fn}: {en:.6f} au")
    assert geom._energy == pytest.approx(ref_energy)


@pytest.fixture
def pyscf_acetaldehyd_getter():
    geom = geom_loader("lib:acetaldehyd_oniom.xyz", coord_type="redund")

    calcs = {
        "high": {
            "type": "pyscf",
            "basis": "321g",
            "verbose": 0,
        },
        "real": {
            "type": "pyscf",
            "basis": "sto3g",
            "verbose": 0,
        },
    }
    models = {
        "high": {
            "inds": [4, 5, 6],
            "calc": "high",
        },
    }

    def get_oniom(**kwargs):
        oniom = ONIOM(calcs, models, geom, layers=None, **kwargs)
        geom.set_calculator(oniom)
        return geom

    return get_oniom


@using("pyscf")
@pytest.mark.parametrize(
    "embedding, ref_energy",
    [
        ("", -151.8130757),
        ("electronic", -151.817564),
        ("electronic_rc", -151.814822),
        ("electronic_rcd", -151.817018),
    ],
)
def test_oniom_ee_charge_distribution(embedding, ref_energy, pyscf_acetaldehyd_getter):
    geom = pyscf_acetaldehyd_getter(embedding=embedding)
    en = geom.energy
    assert en == pytest.approx(ref_energy)


@using("pyscf")
def test_layer_calc(pyscf_acetaldehyd_getter):
    geom = pyscf_acetaldehyd_getter()
    calc = geom.calculator

    real_calc = PySCF(basis="sto3g")
    args = geom.atoms, geom.cart_coords
    #  Reference energies
    real_low_en = real_calc.get_energy(*args)["energy"]
    oniom_en = geom.energy
    ref_energies = (real_low_en, oniom_en)
    # Reference forces
    real_low_forces = real_calc.get_forces(*args)["forces"]
    oniom_forces = geom.cart_forces
    ref_forces = (real_low_forces, oniom_forces)
    # Reference hessians
    real_low_hessian = real_calc.get_hessian(*args)["hessian"]
    oniom_hessian = geom.cart_hessian
    ref_hessians = (real_low_hessian, oniom_hessian)
    energies = list()
    all_forces = list()
    all_hessians = list()
    for i, _ in enumerate(calc.layers):
        lcalc = calc.get_layer_calc(i)
        energy = lcalc.get_energy(*args)["energy"]
        energies.append(energy)
        forces = lcalc.get_forces(*args)["forces"]
        all_forces.append(forces)
        hessian = lcalc.get_hessian(*args)["hessian"]
        all_hessians.append(hessian)
    np.testing.assert_allclose(energies, ref_energies)
    # When using XTB instea of PySCF we don't have to increase the
    # tolerance. There seem to be some numerical instabilities in PySCF.
    np.testing.assert_allclose(all_forces, ref_forces, atol=1e-7)
    np.testing.assert_allclose(all_hessians, ref_hessians, atol=2e-7)
