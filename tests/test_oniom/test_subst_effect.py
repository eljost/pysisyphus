import numpy as np
import pytest

from pysisyphus.calculators import ORCA, Composite
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.testing import using

try:
    from thermoanalysis.QCData import QCData
    from thermoanalysis.thermo import thermochemistry
except ModuleNotFoundError:
    pass


def calc_H(geom, calc_getter, thermo_calc_getter, charge, mult):
    calc = calc_getter(charge=charge, mult=mult)
    thermo_calc = thermo_calc_getter(charge=charge, mult=mult)

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


def hydrogen_bde(geom, mult, diss_geom, diss_mult, calc_getter, thermo_calc_getter):
    H_en = -0.5
    charge = 0

    H = calc_H(geom, calc_getter, thermo_calc_getter, charge, mult)
    diss_H = calc_H(diss_geom, calc_getter, thermo_calc_getter, charge, diss_mult)

    print("\tH", H)
    print("\tdiss_H", diss_H)
    bde = diss_H + H_en - H
    print("\tBDE", bde)

    return bde


@using("thermoanalysis")
@using("orca")
@pytest.mark.parametrize(
    "low_keywords", [
    "mp2 6-31G(d)",
    "hf 6-31G(d)",
    ]
)
def test_calc_bde(low_keywords):
    def thermo_calc_getter(charge=0, mult=1):
        return ORCA("B3lyp_G 6-31G tightscf", charge=charge, mult=mult, pal=6)

    def high_calc_getter(charge=0, mult=1):
        return Composite(
            **{
                "from_dict": {
                    "ccsdt": {
                        "type": "orca",
                        "keywords": "ccsd(t) 6-31G(d) tightscf",
                    },
                    "mp2_high": {
                        "type": "orca",
                        "keywords": "mp2 6-311+G(2df,2p) tightscf rohf",
                    },
                    "mp2_low": {
                        "type": "orca",
                        "keywords": "mp2 6-31G(d) tightscf rohf",
                    },
                },
                "charge": charge,
                "mult": mult,
                "pal": 6,
                "final": "ccsdt + mp2_high - mp2_low",
            }
        )

    def low_calc_getter(charge=0, mult=1):
        return ORCA(low_keywords, charge=charge, mult=mult, pal=6)

    pref = "lib:subst_effect/"
    model = geom_loader(pref + "ch4_model.xyz")
    diss_model = geom_loader(pref + "ch3_b3lypG_631g.xyz")
    real = geom_loader(pref + "toluene_b3lypG_631g.xyz")
    diss_real = geom_loader(pref + "toluene_minus_H_b3lypG_631g.xyz")

    low_real = hydrogen_bde(real, 1, diss_real, 2, low_calc_getter, thermo_calc_getter)
    print("low_real", low_real)
    low_model = hydrogen_bde(model, 1, diss_model, 2, low_calc_getter, thermo_calc_getter)
    print("low_model", low_model)
    S_low = low_real - low_model
    print("S_low", S_low)

    # high_real = hydrogen_bde(real, 1, diss_real, 2, high_calc_getter, thermo_calc_getter)
    # print("high_real", high_real)
    # high_model = hydrogen_bde(model, 1, diss_model, 2, high_calc_getter, thermo_calc_getter)
    # print("high_model", high_model)
    # S_high = high_real - high_model
    S_high = -0.0195031236259027
    print("S_high", S_high)

    dS = S_high - S_low
    dS_kcal = dS * 630
    print(f"@@@{low_keywords}")
    print(f"@@@ dS={dS:.6f} au, dS_kcal={dS_kcal:.2f} kcal/mol")
