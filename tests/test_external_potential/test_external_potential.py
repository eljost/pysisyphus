import numpy as np
import pytest

from pysisyphus.constants import KB, BOHR2ANG
from pysisyphus.calculators import ExternalPotential, XTB
from pysisyphus.calculators.LennardJones import LennardJones
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_library
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.testing import using


init_logging()


def test_lj_external_potential():
    atoms = ("X", )
    coords = (0., 0., 0.)
    geom = Geometry(atoms, coords)

    calc = LennardJones()
    geom.set_calculator(calc)

    ref_energy = geom.energy
    potentials = [
        {
            "type": "logfermi",
            "beta": 6,
            "T": 300,
            "radius": 8,
        },
    ]
    ext_calc = ExternalPotential(calc, potentials=potentials)
    geom.set_calculator(ext_calc)
    ext_energy = geom.energy
    assert ext_energy == pytest.approx(ref_energy)

    displ_coords = (9., 0., 0)
    geom.coords = displ_coords
    ext_energy = geom.energy
    assert ext_energy == pytest.approx(0.00570261)

    forces = geom.forces
    ref_forces = np.array((-0.0056861, 0., 0.))
    np.testing.assert_allclose(forces, ref_forces, atol=1e-6)


def test_geometry_sphere_radius():
    geom = geom_from_library("benzene.xyz")
    radius = geom.get_sphere_radius(offset=4)

    assert radius == pytest.approx(8.7142660)


@using("xtb")
def test_h2o_xtb_opt():
    xtb_kwargs = {
        "pal": 2,
    }
    opt_kwargs = {
        "max_cycles": 500,
    }

    ref_geom = geom_from_library("h2o_30_sphere.xyz")
    ref_com = ref_geom.center_of_mass
    ref_radius = 7 / BOHR2ANG
    print("Reference center of mass", ref_com)
    # ref_geom.set_calculator(XTB(**xtb_kwargs))
    # ref_opt = PreconLBFGS(ref_geom, **opt_kwargs)
    # ref_opt.run()
    # assert ref_opt.is_converged
    # return

    geom = geom_from_library("h2o_30_sphere_translated_2_2_2.xyz")
    com = geom.center_of_mass
    print("Translated center of mass", com)
    org_diff = ref_com - com
    org_diff_norm = np.linalg.norm(org_diff)
    print(f"norm(COM-differences) {org_diff_norm:.6f}")

    # Actual calculator
    xtb_calc = XTB(**xtb_kwargs)

    # External potential
    # Jmol commands
    #   draw sphere diameter 14. {-0.7 0.85 -0.2} translucent
    #   delete $sphere
    ext_kwargs = {
        "potentials": [
            {
                "type": "logfermi",
                "beta": 6,
                "T": 300,
                "radius": ref_radius,
                "origin": ref_com,
            },
        ],
    }
    calc = ExternalPotential(xtb_calc, **ext_kwargs)
    geom.set_calculator(calc)
    print("External potentials")
    for pot in calc.potentials:
        print("\t", pot)

    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    com_opt = geom.center_of_mass
    print()
    print("Original translated center of mass", com)
    print()
    print("Origin of external potential", ref_com)
    print("Optimized center of mass", com_opt)
    opt_diff = ref_com - com_opt
    opt_diff_norm = np.linalg.norm(opt_diff)
    print("Difference to origin of external potential", opt_diff)
    print(f"norm(COM-differences) {opt_diff_norm:.6f}")

    unbiased_energy = xtb_calc.get_energy(geom.atoms, geom.coords)["energy"]
    print(f"Unbiased energy: {unbiased_energy:.6f}")

    assert opt.is_converged
    assert opt_diff_norm < org_diff_norm
    assert opt_diff_norm < 4
