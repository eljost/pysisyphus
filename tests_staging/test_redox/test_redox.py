from collections import namedtuple

import h5py
from scipy.constants import N_A as NA
from thermoanalysis.QCData import QCData
from thermoanalysis.thermo import thermochemistry

from pysisyphus.calculators import Gaussian16, ORCA
from pysisyphus.constants import AU2J
from pysisyphus.helpers import geom_loader, highlight_text
from pysisyphus.io.hessian import save_hessian
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

"""
See
    https://pubs.acs.org/doi/pdf/10.1021/jp025853n
"""


RedoxResult = namedtuple("RedoxResult",
                         "energy_gas hessian_gas energy_solv geom_gas"
)


def run(geom, charge, mult, calc_getter_gas, calc_getter_solv):
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

    res = RedoxResult(
            energy_gas=energyg,
            hessian_gas=Hg,
            energy_solv=energys,
            geom_gas=geom,
    )
    return res


def test_redox(fn, calc_cls, gas_kwargs, solv_kwargs):
    ref_geom = geom_loader(fn, coord_type="redund")

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

    ref_charge_mult = (0, 1)
    oxd_charge_mult = (1, 2)

    ref_res = run(ref_geom, *ref_charge_mult, calc_getter_gas, calc_getter_solv)
    print()

    oxd_geom = ref_res.geom_gas.copy()
    oxd_res = run(oxd_geom, *oxd_charge_mult, calc_getter_gas, calc_getter_solv)


def test_thermo(n=1, T=298.15):
    # Neutral
    gas01 = "gas_0_1.h5"
    solv01 = "solv_0_1.energy"
    # Oxidized
    gas12 = "gas_1_2.h5"
    solv12 = "solv_1_2.energy"

    def process(h5):
        qc = QCData(h5)
        thermo = thermochemistry(qc, temperature=T)
        return qc, thermo

    # Neutral
    qcg01, thermog01 = process(gas01)
    with open(solv01) as handle:
        energy_solv01 = float(handle.read())

    # Oxidized
    qcg12, thermog12 = process(gas12)
    with open(solv12) as handle:
        energy_solv12 = float(handle.read())

    # Farady constant
    F = 96_485.3329  # As mol⁻¹

    def DG_solv(solv_energy, h5_gas):
        with h5py.File(h5_gas, "r") as handle:
            gas_energy = handle.attrs["energy"]
        DG_solv = solv_energy - gas_energy
        return DG_solv

    DG_solv_neut = DG_solv(energy_solv01, gas01)
    DG_solv_oxd = DG_solv(energy_solv12, gas12)
    print(f"ΔG_solv_neut = {DG_solv_neut: >8.6f} E_h")
    print(f"ΔG_solv_oxd  = {DG_solv_oxd: >8.6f} E_h")

    DG_ox_g = thermog12.G - thermog01.G
    print(f"ΔG_oxd_gas  = {DG_ox_g: >8.6f} E_h")
    DG_ox_solv = DG_ox_g + DG_solv_oxd - DG_solv_neut
    print(f"ΔG_oxd_solv  = {DG_ox_solv: >8.6f} E_h")

    print(f"Δelectron(s) = {n}")
    # To Joule mol⁻¹
    redox_pot = DG_ox_solv * AU2J * NA / (n * F)
    print(f"E^0/+ = {redox_pot:.2f} V")

    return redox_pot


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

    fn = "fe0_freq.pdb"
    # test_redox(fn, **g16_kwargs)
    # test_redox(fn, **orca_kwargs)
    test_thermo()
