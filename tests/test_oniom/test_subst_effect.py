import numpy as np
import pytest

from pysisyphus.constants import AU2KCALMOL
from pysisyphus.calculators import ORCA, Composite
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.testing import using
from pysisyphus.run import run_from_dict

try:
    from thermoanalysis.QCData import QCData
    from thermoanalysis.thermo import thermochemistry
except ModuleNotFoundError:
    pass


def hydrogen_bde(geom, mult, diss_geom, diss_mult, calc, thermo_calc):
    H_en = -0.5

    def get_H(geom, mult):
        calc.mult = mult
        thermo_calc.mult = mult
        calc.reset()
        thermo_calc.reset()
        energy = calc.get_energy(geom.atoms, geom.coords)["energy"]
        hessian = thermo_calc.get_hessian(geom.atoms, geom.coords)["hessian"]

        mw_hessian = geom.mass_weigh_hessian(hessian)
        proj_hessian = geom.eckart_projection(mw_hessian)
        eigvals, eigvecs = np.linalg.eigh(proj_hessian)
        vibfreqs = eigval_to_wavenumber(eigvals)

        thermo_dict = {
            "masses": geom.masses,
            "vibfreqs": vibfreqs,
            "coords3d": geom.coords3d,
            "energy": energy,
            "mult": mult,
        }

        qcd = QCData(thermo_dict)
        thermo = thermochemistry(qcd, temperature=298.15)
        H = thermo.H
        return H

    H = get_H(geom, mult)
    diss_H = get_H(diss_geom, diss_mult)

    print("\tH", H)
    print("\tdiss_H", diss_H)
    bde = diss_H + H_en - H
    print("\tBDE", bde)

    return bde


def test_calc_bde():
    thermo_calc = ORCA("B3lyp_G 6-31G tightscf", pal=6)

    print()
    model = geom_loader("ch4_model.xyz")
    diss_model = geom_loader("01_ch3.xyz")

    real = geom_loader("lib:toluene_b3lypG_631g.xyz")
    diss_real = geom_loader("lib:toluene_minus_H_b3lypG_631g.xyz")
    # high_calc = ORCA("BP86 def2-sv(p)", pal=6)
    high_calc = Composite(**{"from_dict": {
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
    "charge": 0,
    "mult": 1,
    "pal": 6,
    "final": "ccsdt + mp2_high - mp2_low",
    }
    )
    high_calc = Composite(**{"from_dict": {
        "quick": {
            "type": "orca",
            "keywords": "pm3",
        },
    },
    "charge": 0,
    "mult": 1,
    "pal": 1,
    "final": "quick",
    }
    )
    low_calc = ORCA("mp2 6-31G(d)", pal=6) 

    low_real = hydrogen_bde(real, 1, diss_real, 2, low_calc, thermo_calc)
    print("low_real", low_real)
    low_model = hydrogen_bde(model, 1, diss_model, 2, low_calc, thermo_calc)
    print("low_model", low_model)
    S_low = low_real - low_model
    print("S_low", S_low)

    high_real = hydrogen_bde(real, 1, diss_real, 2, high_calc, thermo_calc)
    print("high_real", high_real)
    high_model = hydrogen_bde(model, 1, diss_model, 2, high_calc, thermo_calc)
    print("high_model", high_model)
    S_high = high_real - high_model
    print("S_high", S_high)

    dS = S_low - S_high
    dS_kcal = dS * 630
    print(f"dS={dS:.6f} au, dS_kcal={dS_kcal:.2f} kcal/mol")


def test_oniom_s_value():
    geom = geom_loader("lib:toluene_b3lypG_631g.xyz")

    from pysisyphus.calculators import ORCA

    def get_calc():
        return ORCA("rhf 6-31G(d)", pal=6)

    calc = get_calc()

    geom.set_calculator(calc)
    parent_calc = calc

    from pysisyphus.calculators.ONIOMv2 import Model

    mdl = Model(
        "model", "high", calc, "real", "low", parent_calc, (0, 7, 8, 9), list(range(15))
    )
    atoms = geom.atoms
    coords3d = geom.coords3d
    mdl.create_links(atoms, coords3d)
    ha, hc = mdl.capped_atoms_coords(atoms, coords3d)
    from pysisyphus.Geometry import Geometry

    high = Geometry(ha, hc)

    high.set_calculator(get_calc())

    model_en = high.energy
    print("model en", model_en)

    real_en = geom.energy
    print("real en", real_en)

    print("methane")
    ch4_g2ms, *_ = g2ms(high, mos=9, do_opt=False, do_freq=False, pal=6)
    print()
    print("toluene")
    tol_g2ms, *_ = g2ms(geom, mos=25, do_opt=False, do_freq=False, pal=6)
    s_high = tol_g2ms - ch4_g2ms
    print(f"S_high: {s_high:.6f}")
    from pysisyphus.constants import AU2KCALMOL

    print(f"S_high: {s_high*AU2KCALMOL:.6f}")
    # s_high = -230.693594  # with hlc == -3.8 kcal/mol

    low_calc = ORCA("mp2 6-31G(d) tightscf", pal=6)
    ch4_low = low_calc.get_energy(high.atoms, high.cart_coords)["energy"]
    low_calc = ORCA("mp2 6-31G(d) tightscf", pal=6)
    tol_low = low_calc.get_energy(geom.atoms, geom.cart_coords)["energy"]
    s_low = tol_low - ch4_low
    print(f"S_low: {s_low:.6f} au")
    dS = s_high - s_low
    print(f"ΔS={dS:.6f} au")
    dS_kcal = dS * AU2KCALMOL
    print(f"ΔS={dS_kcal:.6f} kcal/mol")


@using("thermoanalysis")
def g2ms(geom, mos, do_opt=False, do_freq=True, pal=6, charge=0, mult=1):
    cmp_ = {
        "charge": charge,
        "mult": mult,
        "pal": pal,
    }
    opt_calc = ORCA("b3lyp_G 6-31G(d) tightscf", **cmp_)
    if do_opt:
        geom.set_calculator(opt_calc)
        opt_kwargs = {
            "thresh": "gau_tight",
            "dump": True,
        }
        opt = RFOptimizer(geom, **opt_kwargs)
        opt.run()
        assert opt.is_converged

    ccsdt = ORCA("ccsd(t) 6-31G(d) tightscf", **cmp_)
    mp2_high = ORCA("mp2 6-311+G(2df,2p) tightscf", **cmp_)
    mp2_low = ORCA("mp2 6-31G(d) tightscf", **cmp_)

    atoms = geom.atoms
    coords = geom.cart_coords

    def get_energy(calc):
        return calc.get_energy(atoms, coords)["energy"]

    ccsdt_en = get_energy(ccsdt)
    mp2_high_en = get_energy(mp2_high)
    mp2_low_en = get_energy(mp2_low)
    print(f"ccsdt: {ccsdt_en:.6f} au")
    print(f"mp2/high: {mp2_high_en:.6f} au")
    print(f"mp2/low: {mp2_low_en:.6f} au")

    diff = ccsdt_en + mp2_high_en - mp2_low_en
    print(f"diff: {diff:.6f} au")
    hlc = -3.8 / AU2KCALMOL
    print(f"mos*hlc: {mos*hlc:.6f} au")

    g2ms = diff + mos * hlc
    print(f"g2ms: {g2ms:.6f} au")

    zpe = None
    g2ms_zpe = None
    if do_freq:
        hessian = opt_calc.get_hessian(geom.atoms, geom.cart_coords)["hessian"]
        mw_hessian = geom.mass_weigh_hessian(hessian)
        proj_hessian = geom.eckart_projection(mw_hessian)
        eigvals, eigvecs = np.linalg.eigh(proj_hessian)
        vibfreqs = eigval_to_wavenumber(eigvals)

        thermo_dict = {
            "masses": geom.masses,
            "vibfreqs": vibfreqs,
            "coords3d": geom.coords3d,
            "energy": g2ms,
            "mult": mult,
        }

        qcd = QCData(thermo_dict)
        thermo = thermochemistry(qcd, temperature=298.15)
        zpe = thermo.ZPE
        g2ms_zpe = g2ms + zpe

    return g2ms, zpe, g2ms_zpe


def test_g2ms():
    geom = geom_loader("h2_g2ms_opt.xyz")
    h2_g2ms, zpe, h2_g2ms_zpe = g2ms(geom, mos=1, do_opt=False, pal=1)

    h_geom = geom_loader("h.xyz")
    calc = ORCA("hf 6-311G(d,p) tightscf", mult=2)
    en = calc.get_energy(h_geom.atoms, h_geom.coords)["energy"]

    d = h2_g2ms_zpe - 2 * en
    d_kcal = d * AU2KCALMOL

    print(f"atomization energy: {d:.6f} au")
    print(f"atomization energy: {d_kcal:.6f} kcal_mol")

    ref = 103.3

    diff = d_kcal + ref
    print(f"Δ={diff:.2f} kcal/mol")


def test_toluene():
    tol = geom_loader("lib:toluene_b3lypG_631g.xyz")
    mH = geom_loader("lib:toluene_minus_H_b3lypG_631g.xyz")

    print("Toluene")
    tol_g2ms, *_ = g2ms(
        tol, mos=0, do_opt=False, do_freq=False, pal=6, charge=0, mult=1
    )

    print("Toluene minus H")
    mH_g2ms, *_ = g2ms(mH, mos=0, do_opt=False, do_freq=False, pal=6, charge=0, mult=2)


@using("orca")
@pytest.mark.parametrize(
    "fn, mult, high_inds, ref_energy",
    [
        ("lib:toluene_b3lypG_631g.xyz", 1, (0, 7, 8, 9), -271.470945476921),
        ("lib:toluene_minus_H_b3lypG_631g.xyz", 2, (0, 7, 8), -270.824806805671),
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
                    "from_dict": {
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


def test_deltas():
    from pysisyphus.calculators import Composite
    from pysisyphus.Geometry import Geometry
    from pysisyphus.helpers import geom_loader
    calc = Composite(**{"from_dict": {
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
    "charge": 0,
    "mult": 1,
    "pal": 6,
    "final": "ccsdt + mp2_high - mp2_low",
    }
    )
    parent_calc = ORCA("HF 6-31G(d) tightscf", pal=6)

    geom = geom_loader("lib:toluene_b3lypG_631g.xyz")
    # real_en = calc.get_energy(geom.atoms, geom.coords)["energy"]
    # real_en = -271.023663500836
    # print("real en", real_en)

    from pysisyphus.calculators.ONIOMv2 import Model
    mdl = Model(
        "model", "high", calc, "real", "low", parent_calc, (0, 7, 8, 9), list(range(15))
    )
    mdl.create_links(geom.atoms, geom.coords)
    # dS = mdl.get_delta_S(geom.atoms, geom.coords)
    model_atoms, model_coords = mdl.capped_atoms_coords(geom.atoms, geom.coords)
    model = Geometry(model_atoms, model_coords)
    with open("ch4_model.xyz", "w") as handle:
        handle.write(model.as_xyz())
    model.jmol()
    # model_en = calc.get_energy(model.atoms, model.coords)["energy"]
    # print("model en", model_en)
    # dS = real_en - model_en
    # dS_kcal = dS * 630
    # print(f"dS={dS:.6f} au, dS={dS_kcal:.2f} kcal/mol")
