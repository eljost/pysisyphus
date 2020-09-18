import numpy as np
import pytest

from pysisyphus.calculators.ONIOMv2 import ONIOM
from pysisyphus.helpers import do_final_hessian, geom_loader
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using_gaussian16, using_pyscf


init_logging()


@using_gaussian16
def test_energy():
    geom = geom_loader("lib:alkyl17_sto3g_opt.xyz")

    real = set(range(len(geom.atoms)))
    medmin = set((0,1,2,3,4,5,6, 46,47,48,49,50,51,52))
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
        "med" : {
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
        }
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
            {"real": {"type": "g16", "route": "hf sto-3g"},
             "high": {"type": "g16", "route": "b3lyp d95v"}},
            -153.07432042299052 , 0.03768246934785125,
            marks=using_gaussian16,
        ),
        # The following two tests should yield identical results
        pytest.param(
            {"real": {"type": "g16", "route": "hf sto-3g"},
             "high": {"type": "g16", "route": "b3lyp 3-21g"}},
            -152.4529060634755 , 0.018462670668992546,
            marks=using_gaussian16,
        ),
        pytest.param(
            {"real": {"type": "pyscf", "basis": "sto3g"},
             "high": {"type": "pyscf", "xc": "b3lypg", "basis": "321g"}},
            -152.4529060634755, 0.01839279960703439,
            marks=using_pyscf,
        ),
])
def test_gradient(calcs, ref_energy, ref_force_norm):
    geom = geom_loader("lib:acetaldehyd_oniom.xyz", coord_type="redund")

    real = set(range(len(geom.atoms)))
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
        pytest.param("g16", None,   -582.392035, 0.085568849,
                     marks=using_gaussian16),
        pytest.param("pyscf", None, -582.392035, 0.088564339,
                     marks=using_pyscf),

        # Electronic embedding
        pytest.param("g16", "electronic",   -582.3997769406087, 0.08582761,
                     marks=using_gaussian16),
        pytest.param("pyscf", "electronic", -582.3997769406087, 0.08878196,
                     marks=using_pyscf),
])
def test_electronic_embedding(calc_key, embedding, ref_energy, ref_force_norm):
    geom = geom_loader("lib:oniom_ee_model_system.xyz", coord_type="redund")

    all_ = set(range(len(geom.atoms)))
    high = list(sorted(all_ - set((21, 20, 19, 15, 14, 13))))

    calcs_dict = {
        "g16": ({"real": {"type": "g16", "route": "hf sto-3g"},
                 "high": {"type": "g16", "route": "hf 3-21g"},
        }),
        "pyscf": ({"real": {"type": "pyscf", "basis": "sto3g",},
                   "high": {"type": "pyscf", "basis": "321g"},
        }),
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


@using_gaussian16
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
        }
    }

    # "real" must not be defined; will be automatically set
    layers = ["mid", "high"]

    oniom = ONIOM(calcs, models, geom, layers)

    assert oniom.layer_num == 3

    geom.set_calculator(oniom)

    assert geom.energy == pytest.approx(-77.420587)


@using_gaussian16
def test_acetaldehyde_opt():
    """
        From https://doi.org/10.1016/S0166-1280(98)00475-8
    """
    geom = geom_loader("lib:acetaldehyd_oniom.xyz", coord_type="redund")

    real = set(range(len(geom.atoms)))
    high = [4, 5, 6]

    calcs = {
        "real": {"type": "g16", "route": "hf sto-3g"},
         "high": {"type": "g16", "route": "b3lyp d95v"}
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


@using_gaussian16
def test_yaml_oniom():

    run_dict = {
        "xyz": "lib:acetaldehyd_oniom.xyz",
        "coord_type": "redund",
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
            }
        },
        "opt": {
            "thresh": "gau_tight",
        }
    }
    res = run_from_dict(run_dict)

    opt = res.opt
    assert opt.is_converged
    assert opt.cur_cycle == 7
    assert res.opt_geom.energy == pytest.approx(-153.07526171)


@using_pyscf
def test_oniom3():
    run_dict = {
        "xyz": "lib:oniom3alkyl.pdb",
        "coord_type": "redund",
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
            }
        },
        "opt": {
            "thresh": "gau_tight",
        },
    }
    res = run_from_dict(run_dict)
    print()

    opt = res.opt
    assert opt.is_converged
    assert opt.cur_cycle == 6

    geom = res.opt_geom
    res = do_final_hessian(geom, save_hessian=False)
    nus = res.nus
    assert nus[-1] == pytest.approx(3747.54937)
    assert nus[-5] == pytest.approx(3563.89449)


@pytest.mark.skip
@using_pyscf
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
            }
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
@using_pyscf
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
            }
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
