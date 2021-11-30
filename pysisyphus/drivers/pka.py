from math import log

from pysisyphus.constants import AU2KCALMOL, KBAU
from pysisyphus.io import geom_from_hessian

def direct_cycle(
    acid_h5,
    base_h5,
    acid_solv_en,
    base_solv_en,
    G_aq_H=None,
    G_gas_H=-6.28,
    dG_solv_H=-265.9,
    T=298.15,
):
    def G_aq_from_hessian(h5, solv_en):
        geom = geom_from_hessian(h5)
        thermo = geom.get_thermoanalysis(T=T)
        dG_solv = solv_en - thermo.U_el
        G_aq = thermo.G + dG_solv
        return G_aq

    G_aq_acid = G_aq_from_hessian(acid_h5, acid_solv_en)
    G_aq_base = G_aq_from_hessian(base_h5, base_solv_en)

    if G_aq_H is None:
        G_aq_H = (G_gas_H + dG_solv_H) / AU2KCALMOL

    dG_aq = G_aq_H + G_aq_base - G_aq_acid
    corr = KBAU * T * log(24.46)
    dG_aq_corr = dG_aq + corr

    pKa = dG_aq_corr / (KBAU * T * log(10))
    return pKa
