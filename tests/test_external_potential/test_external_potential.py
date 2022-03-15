import numpy as np
import pytest

from pysisyphus.constants import BOHR2ANG
from pysisyphus.calculators import ExternalPotential, XTB
from pysisyphus.calculators.ExternalPotential import HarmonicSphere, Restraint
from pysisyphus.calculators.LennardJones import LennardJones
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


init_logging()


def test_lj_external_potential():
    atoms = ("X",)
    coords = (0.0, 0.0, 0.0)
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

    displ_coords = (9.0, 0.0, 0)
    geom.coords = displ_coords
    ext_energy = geom.energy
    assert ext_energy == pytest.approx(0.00570261)

    forces = geom.forces
    ref_forces = np.array((-0.0056861, 0.0, 0.0))
    np.testing.assert_allclose(forces, ref_forces, atol=1e-6)


def test_geometry_sphere_radius():
    geom = geom_loader("lib:benzene.xyz")
    radius = geom.get_sphere_radius(offset=4)

    assert radius == pytest.approx(8.7142660)


@pytest.mark.skip
@using("xtb")
def test_h2o_xtb_opt():
    xtb_kwargs = {
        "pal": 2,
    }
    opt_kwargs = {
        "max_cycles": 500,
    }

    ref_geom = geom_loader("lib:h2o_30_sphere.xyz")
    ref_com = ref_geom.center_of_mass
    ref_radius = 7 / BOHR2ANG
    print("Reference center of mass", ref_com)
    # ref_geom.set_calculator(XTB(**xtb_kwargs))
    # ref_opt = PreconLBFGS(ref_geom, **opt_kwargs)
    # ref_opt.run()
    # assert ref_opt.is_converged
    # return

    geom = geom_loader("lib:h2o_30_sphere_translated_2_2_2.xyz")
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


@pytest.mark.skip
def test_lj_external_potential_opt():
    np.random.seed(20182503)

    radius = 9
    atom_num = 50
    coords3d = (np.random.rand(atom_num, 3) - 0.5) * 2 * radius
    atoms = ("Ar",) * atom_num

    geom = Geometry(atoms, coords3d.flatten())

    lj_calc = LennardJones()
    geom.set_calculator(lj_calc)

    opt_kwargs = {
        "max_cycles": 250,
        "precon_getter_update": 50,
        "c_stab": 0.5,
    }
    # opt = QuickMin(geom, **opt_kwargs)
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()


@pytest.mark.parametrize(
    "coords, ref_energy, ref_norm_forces",
    [
        ((0.0, 0.0, 0.0), 0.0, 0.0),
        ((1.0, 0.0, 0.0), 0.0, 0.0),
        ((2.0, 0.0, 0.0), 4.0, 4.0),
        ((2.0, 2.0, 0.0), 8.0, 2 * 8.0 ** 0.5),
        ((4.0, 0.0, 0.0), 16.0, 8.0),
    ],
)
def test_harmonic_sphere(coords, ref_energy, ref_norm_forces):
    atoms = ("X",)
    coords = np.array(coords)
    geom = Geometry(atoms, coords)

    potentials = [
        {
            "type": "harmonic_sphere",
            "k": 1,
            "radius": 1,
        },
    ]
    calc = ExternalPotential(potentials=potentials)
    geom.set_calculator(calc)
    energy = geom.energy
    assert energy == pytest.approx(ref_energy)

    forces = geom.forces
    norm_forces = np.linalg.norm(forces)
    assert norm_forces == pytest.approx(ref_norm_forces)


@pytest.mark.parametrize(
    "coords, ref_pressure",
    [
        ((0.0, 0.0, 0.0), 0.0),
        ((1.0, 0.0, 0.0), 0.0),
        ((2.0, 0.0, 0.0), 0.31830989),
        ((2 * 2.0, 0.0, 0.0), 2 * 0.31830989),
    ],
)
def test_harmonic_sphere_pressure(coords, ref_pressure):
    c3d = np.reshape(coords, (-1, 3))
    calc = HarmonicSphere(k=1, radius=1.0)

    assert calc.surface_area == pytest.approx(4 * np.pi)

    p = calc.instant_pressure(c3d)
    assert p == pytest.approx(ref_pressure)


# def test_h2o_restraint():
# geom = geom_loader("lib:h2o.xyz")
# print("atoms", geom.atoms)
# from pysisyphus.constants import ANG2BOHR
# restraints = [
# [["BOND", 0, 1], 10000, 2.0*ANG2BOHR],
# # [["BOND", 1, 2], 10, 2.5],
# # [["BOND", 0, 2], 2, 1.8],
# # [["BOND", 2, 3], 0.1, 1.6],
# # [["BEND", 1, 0, 2], 1.0, np.deg2rad(140)],
# ]
# # restraints = []
# calc = XTB(quiet=True)
# R = Restraint(restraints, geom)
# en, grad = R.calc(geom.coords3d, gradient=True)

# potentials = [
# {
# "type": "restraint",
# "restraints": restraints,
# "geom": geom,
# },
# ]
# ext = ExternalPotential(calc, potentials=potentials)
# geom.set_calculator(ext)
# # from pysisyphus.optimizers.RFOptimizer import RFOptimizer
# # opt = RFOptimizer(geom, dump=True)
# from pysisyphus.optimizers.SteepestDescent import SteepestDescent
# opt = SteepestDescent(geom, dump=True, max_step=0.05, max_cycles=1000)
# opt.run()


@using("pyscf")
@pytest.mark.parametrize("ref_val", np.linspace(0.6, 1.5, num=5))
def test_h2_restraint(ref_val):
    geom = geom_loader("lib:h2.xyz", coord_type="redund")
    restraints = [
        [["BOND", 0, 1], 10, ref_val],
    ]
    calc = PySCF(basis="sto3g", verbose=0)
    R = Restraint(restraints, geom)
    en, grad = R.calc(geom.coords3d, gradient=True)

    potentials = [
        {
            "type": "restraint",
            "restraints": restraints,
            "geom": geom,
        },
    ]
    ext = ExternalPotential(calc, potentials=potentials)
    geom.set_calculator(ext)
    opt = RFOptimizer(geom, dump=True)
    opt.run()
    val = geom.coords[0]
    print(f"@@@ {ref_val:.8f} {geom.coords[0]:.8f})")
    assert val == pytest.approx(ref_val, abs=0.08)


@using("pyscf")
@pytest.mark.parametrize("ref_val", np.linspace(-np.pi / 2, np.pi / 2, 5))
def test_torsion_restraint(ref_val):
    geom = geom_loader("lib:h2o2_hf_321g_opt.xyz", coord_type="redund")
    restraints = [
        [["TORSION", 1, 3, 2, 0], 0.5, ref_val],
    ]
    calc = PySCF(basis="sto3g", verbose=0)
    R = Restraint(restraints, geom)
    en, grad = R.calc(geom.coords3d, gradient=True)

    potentials = [
        {
            "type": "restraint",
            "restraints": restraints,
            "geom": geom,
        },
    ]
    ext = ExternalPotential(calc, potentials=potentials)
    geom.set_calculator(ext)
    opt = RFOptimizer(geom, dump=True)
    opt.run()
    val = geom.coords[-1]
    val = np.rad2deg(val)
    ref_val = np.rad2deg(ref_val)
    print(f"@@@ ref={ref_val:.4f} val={val:.4f} Δ={val-ref_val:.4f})")
    # assert val == pytest.approx(ref_val, abs=0.08)
