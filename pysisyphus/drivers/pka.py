from math import log

from pysisyphus.config import p_DEFAULT, T_DEFAULT
from pysisyphus.constants import AU2KCALPERMOL, KBAU
from pysisyphus.helpers_pure import standard_state_corr
from pysisyphus import logger
from pysisyphus.io import geom_from_hessian


def G_aq_from_h5_hessian(h5_hessian, solv_en, T=T_DEFAULT, p=p_DEFAULT):
    geom = geom_from_hessian(h5_hessian)
    thermo = geom.get_thermoanalysis(T=T, p=p)
    dG_solv = solv_en - thermo.U_el
    G_aq = float(thermo.G + dG_solv)
    return G_aq


def direct_cycle(
    acid_h5,
    base_h5,
    acid_solv_en,
    base_solv_en,
    G_aq_H=None,
    G_gas_H=-6.28,
    dG_solv_H=-265.9,
    T=T_DEFAULT,
    p=p_DEFAULT,
):

    G_aq_acid = G_aq_from_h5_hessian(acid_h5, acid_solv_en, T=T, p=p)
    G_aq_base = G_aq_from_h5_hessian(base_h5, base_solv_en, T=T, p=p)

    if G_aq_H is None:
        G_aq_H = G_gas_H + dG_solv_H
    G_aq_H /= AU2KCALPERMOL

    dG_aq = G_aq_H + G_aq_base - G_aq_acid
    dG_aq_corr = dG_aq + standard_state_corr(T=T, p=p)

    pKa = dG_aq_corr / (KBAU * T * log(10))

    logger.debug(
        f"pKa calculation using direct cycle:\n"
        f"\tAcid HDF5: {acid_h5}\n\tBase HDF5: {base_h5}\n"
        f"\tT={T:.4f}, p={p:.1f} Pa\n"
        f"\tG_aq_acid={G_aq_acid:.6f} au\n"
        f"\tG_aq_base={G_aq_base:.6f} au\n"
        f"\tG_aq_H={G_aq_H:.6f} au\n"
        f"\tpKa={pKa:.6f}"
    )

    return pKa
