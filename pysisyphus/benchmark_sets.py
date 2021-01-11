def get_s22_fns():
    """https://doi.org/10.1039/B600027D"""
    fns = (
        "00_adenine_thymine_wc.xyz",
        "01_adenine_thymine_stack.xyz",
        "02_ammonia_dimer.xyz",
        "03_water_dimer.xyz",
        "04_methane_dimer.xyz",
        "05_ethene_dimer.xyz",
        "06_ethene_ethine.xyz",
        "07_formic_acid_dimer.xyz",
        "08_formamide_dimer.xyz",
        "09_benzene_water.xyz",
        "10_benzene_ammonia.xyz",
        "11_benzene_methane.xyz",
        "12_benzene_dimer_c2v.xyz",
        "13_benzene_dimer_c2h.xyz",
        "14_indole_benzene_t-shape.xyz",
        "15_indole_benzene_stack.xyz",
        "16_pyrazine_dimer.xyz",
        "17_2-pyridoxine_2-aminopyridine.xyz",
        "18_phenol_dimer.xyz",
        "19_uracil_dimer_stack.xyz",
        "20_uracil_dimer_hb.xyz",
        "21_benzene_hcn.xyz",
    )
    prefix = "lib:s22/"
    return prefix, fns


def get_baker_fns():
    """"""
    fns = (
        "00_water.xyz",
        "01_ammonia.xyz",
        "02_ethane.xyz",
        "03_acetylene.xyz",
        "04_allene.xyz",
        "05_hydroxysulphane.xyz",
        "06_benzene.xyz",
        "07_methylamine.xyz",
        "08_ethanol.xyz",
        "09_acetone.xyz",
        "10_disilylether.xyz",
        "11_135trisilacyclohexane.xyz",
        "12_benzaldehyde.xyz",
        "13_13difluorobenzene.xyz",
        "14_135trifluorobenzene.xyz",
        "15_neopentane.xyz",
        "16_furan.xyz",
        "17_naphthalene.xyz",
        "18_15difluoronaphthalene.xyz",
        "19_2hydroxybicyclopentane.xyz",
        "20_achtar10.xyz",
        "21_acanil01.xyz",
        "22_benzidine.xyz",
        "23_pterin.xyz",
        "24_difuropyrazine.xyz",
        "25_mesityloxide.xyz",
        "26_histidine.xyz",
        "27_dimethylpentane.xyz",
        "28_caffeine.xyz",
        "29_menthone.xyz",
    )
    prefix = "lib:baker/"
    return prefix, fns


def get_baker_ref_energies():
    """HF/STO-3G"""
    ref_energies = {
        "00_water.xyz": -74.96590,
        "01_ammonia.xyz": -55.45542,
        "02_ethane.xyz": -78.30618,
        "03_acetylene.xyz": -75.85625,
        "04_allene.xyz": -114.42172,
        "05_hydroxysulphane.xyz": -468.12592,
        "06_benzene.xyz": -227.89136,
        "07_methylamine.xyz": -94.01617,
        "08_ethanol.xyz": -152.13267,
        "09_acetone.xyz": -189.53603,
        "10_disilylether.xyz": -648.58003,
        "11_135trisilacyclohexane.xyz": -976.13242,
        "12_benzaldehyde.xyz": -339.12084,
        "13_13difluorobenzene.xyz": -422.81106,
        "14_135trifluorobenzene.xyz": -520.27052,
        "15_neopentane.xyz": -194.04677,
        "16_furan.xyz": -225.75126,
        "17_naphthalene.xyz": -378.68685,
        "18_15difluoronaphthalene.xyz": -573.60633,
        "19_2hydroxybicyclopentane.xyz": -265.46482,
        "20_achtar10.xyz": -356.28265,
        "21_acanil01.xyz": -432.03012,
        "22_benzidine.xyz": -563.27798,
        "23_pterin.xyz": -569.84884,
        "24_difuropyrazine.xyz": -556.71910,
        "25_mesityloxide.xyz": -304.05919,
        "26_histidine.xyz": -538.54910,
        "27_dimethylpentane.xyz": -271.20088,
        "28_caffeine.xyz": -667.73565,
        "29_menthone.xyz": -458.44639,
    }
    return ref_energies
