from collections import namedtuple

import h5py
import numpy as np
from scipy.constants import N_A as NA
from thermoanalysis.QCData import QCData
from thermoanalysis.thermo import thermochemistry

from pysisyphus.calculators import Gaussian16, ORCA
from pysisyphus.constants import AU2J
from pysisyphus.helpers import geom_loader, highlight_text, eigval_to_wavenumber
from pysisyphus.io.hessian import save_hessian
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using

"""
See
    https://pubs.acs.org/doi/pdf/10.1021/jp025853n
    https://www.sciencedirect.com/science/article/pii/S0009261404012977
    Grimme review
    https://pubs.acs.org/doi/pdf/10.1021/acs.jpca.0c05052
"""

# Farady constant
F = 96_485.3329  # As mol⁻¹

RedoxResult = namedtuple("RedoxResult",
                         "energy_gas hessian_gas energy_solv dG_solv geom_gas charge mult"
)



def redox_result_to_qcdata(redox_result):
    geom = redox_result.geom_gas
    H = geom.eckart_projection(geom.mass_weigh_hessian(redox_result.hessian_gas))
    w, v = np.linalg.eigh(H)
    nus = eigval_to_wavenumber(w)

    thermo_dict = {
        "masses": geom.masses,
        "vibfreqs": nus,
        "coords3d": geom.coords3d,
        "energy": redox_result.energy_gas,
        "mult": redox_result.mult,
    }
    qcd = QCData(thermo_dict)
    return qcd


def calculate_potential(result_1, result_2, T=298.15):
    qcd_1 = redox_result_to_qcdata(result_1)
    qcd_2 = redox_result_to_qcdata(result_2)

    thermo_1 = thermochemistry(qcd_1, temperature=T)
    thermo_2 = thermochemistry(qcd_2, temperature=T)

    DG_solv_neut = result_1.dG_solv
    DG_solv_oxd = result_2.dG_solv
    print(f"ΔG_solv_neut = {DG_solv_neut: >8.6f} E_h")
    print(f"ΔG_solv_oxd  = {DG_solv_oxd: >8.6f} E_h")

    DG_ox_g = thermo_2.G - thermo_1.G
    print(f"ΔG_oxd_gas  = {DG_ox_g: >8.6f} E_h")
    DG_ox_solv = DG_ox_g + DG_solv_oxd - DG_solv_neut
    print(f"ΔG_oxd_solv  = {DG_ox_solv: >8.6f} E_h")

    n = abs(result_1.charge - result_2.charge)
    print(f"Δelectron(s) = {n}")
    # To Joule mol⁻¹
    redox_pot = DG_ox_solv * AU2J * NA / (n * F)
    print(f"E^0/+ = {redox_pot:.2f} V")

    return redox_pot


def run_calculations(geom, charge, mult, calc_getter_gas, calc_getter_solv, opt=False):
    def get_name(base): return f"{base}_{charge}_{mult}"

    print(highlight_text(f"Charge={charge}, Mult={mult}"))

    gas_name = get_name("gas")
    calc_kwargs = {
        "charge": charge,
        "mult": mult,
        "base_name": gas_name,
    }

    opt_kwargs = {
        "thresh": "gau",
        "h5_group_name": f"opt_{charge}_{mult}",
        "dump": True,
        "prefix": f"{gas_name}_",
    }

    # Gas phase calculation, optimization and frequency
    gas_calc = calc_getter_gas(calc_kwargs)
    geom.set_calculator(gas_calc)
    if opt:
        opt = RFOptimizer(geom, **opt_kwargs)
        opt.run()
        assert opt.is_converged
        print(highlight_text("Gas phase optimization finished!", level=1))
    Hg = geom.cart_hessian
    energyg = geom.energy
    save_hessian(f"{get_name('gas')}.h5", geom,
                 cart_hessian=Hg, energy=energyg, mult=mult)
    print(highlight_text("Gas phase hessian finished!", level=1))

    # Solvent calculation, frequency
    solv_name = get_name("solv")
    solv_kwargs = calc_kwargs.copy()
    solv_kwargs["base_name"] = solv_name
    solv_calc = calc_getter_solv(solv_kwargs)
    solv_geom = geom.copy()
    solv_geom.set_calculator(solv_calc)
    energys = solv_geom.energy
    with open(f"{solv_name}.energy", "w") as handle:
        handle.write(str(energys))
    print(highlight_text("Solvated energy finished!", level=1))

    dG_solv = energys - energyg

    res = RedoxResult(
            energy_gas=energyg,
            hessian_gas=Hg,
            energy_solv=energys,
            dG_solv=dG_solv,
            geom_gas=geom,
            charge=charge,
            mult=mult,
    )
    return res


def test_redox(inp_1, inp_2, calc_cls, gas_kwargs, solv_kwargs):
    geom_1 = geom_loader(inp_1["xyz"], coord_type="redund")
    charge_1 = inp_1["charge"]
    mult_1 = inp_1["mult"]
    opt_1 = inp_1["opt"]

    geom_2 = geom_loader(inp_2["xyz"], coord_type="redund")
    charge_2= inp_2["charge"]
    mult_2 = inp_2["mult"]
    opt_2 = inp_2["opt"]


    gas_calc_number = 0
    def calc_getter_gas(calc_kwargs):
        nonlocal gas_calc_number
        calc = calc_cls(**gas_kwargs, calc_number=gas_calc_number, **calc_kwargs)
        gas_calc_number += 1
        return calc

    solv_calc_number = 0
    def calc_getter_solv(calc_kwargs):
        nonlocal solv_calc_number
        calc = calc_cls(**solv_kwargs, calc_number=gas_calc_number, **calc_kwargs)
        solv_calc_number += 1
        return calc

    result_1 = run_calculations(geom_1, charge_1, mult_1,
                                calc_getter_gas, calc_getter_solv, opt=opt_1)
    print()
    result_2 = run_calculations(geom_2, charge_2, mult_2,
                                calc_getter_gas, calc_getter_solv, opt=opt_2)

    calculate_potential(result_1, result_2, T=298.15)


@using("orca")
def test_nitrobenzene():
    inp_1 = {
        "xyz": "lib:redox/nitrobenzene.xyz",
        "charge": 0,
        "mult": 1,
        "opt": False,
    }
    inp_2 = {
        "xyz": "lib:redox/nitrobenzene_anion.xyz",
        "charge": -1,
        "mult": 2,
        "opt": False,
    }

    orca_gas = {
        "keywords": "uks b3lyp 6-31G** rijcosx autoaux",
        "pal": 4,
    }
    orca_solv = orca_gas.copy()
    orca_solv["keywords"] += " cpcm(acetonitrile)"
    orca_kwargs = {
        "calc_cls": ORCA,
        "gas_kwargs": orca_gas,
        "solv_kwargs": orca_solv,
    }
    test_redox(inp_1, inp_2, **orca_kwargs)
    import pdb; pdb.set_trace()
    pass


if __name__ == "__main__":
    ##############
    # Gaussian16 #
    ##############

    g16_gas = {
        "route": "upm6",
        "pal": 4,
    }
    g16_solv = g16_gas.copy()
    g16_solv["route"] += " scrf(pcm,solvent=dimethylsulfoxide)"
    g16_kwargs = {
        "calc_cls": Gaussian16,
        "gas_kwargs": g16_gas,
        "solv_kwargs": g16_solv,
    }

    ########
    # ORCA #
    ########

    orca_gas = {
        "keywords": "uks bp86 def2-svp",
        "pal": 4,
    }
    orca_solv = orca_gas.copy()
    orca_solv["keywords"] += " cpcm(dmso)"
    orca_kwargs = {
        "calc_cls": ORCA,
        "gas_kwargs": orca_gas,
        "solv_kwargs": orca_solv,
    }

    ############
    # FERROCEN #
    ############

    fn = "lib:redox/ferrocen_pm6.xyz"
    inp_1 = {
        "xyz": fn,
        "charge": 0,
        "mult": 1,
        "opt": False,
    }
    inp_2 = {
        "xyz": fn,
        "charge": 1,
        "mult": 2,
        "opt": True,
    }

    test_redox(inp_1, inp_2, **g16_kwargs)
    # test_redox(fn, **orca_kwargs)
