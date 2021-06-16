def get_baker_data():
    """10.1002/jcc.540140910

    HF/STO-3G
    """
    data = (
        ("00_water.xyz", 0, 1, -74.96590),
        ("01_ammonia.xyz", 0, 1, -55.45542),
        ("02_ethane.xyz", 0, 1, -78.30618),
        ("03_acetylene.xyz", 0, 1, -75.85625),
        ("04_allene.xyz", 0, 1, -114.42172),
        ("05_hydroxysulphane.xyz", 0, 1, -468.12592),
        ("06_benzene.xyz", 0, 1, -227.89136),
        ("07_methylamine.xyz", 0, 1, -94.01617),
        ("08_ethanol.xyz", 0, 1, -152.13267),
        ("09_acetone.xyz", 0, 1, -189.53603),
        ("10_disilylether.xyz", 0, 1, -648.58003),
        ("11_135trisilacyclohexane.xyz", 0, 1, -976.13242),
        ("12_benzaldehyde.xyz", 0, 1, -339.12084),
        ("13_13difluorobenzene.xyz", 0, 1, -422.81106),
        ("14_135trifluorobenzene.xyz", 0, 1, -520.27052),
        ("15_neopentane.xyz", 0, 1, -194.04677),
        ("16_furan.xyz", 0, 1, -225.75126),
        ("17_naphthalene.xyz", 0, 1, -378.68685),
        ("18_15difluoronaphthalene.xyz", 0, 1, -573.60633),
        ("19_2hydroxybicyclopentane.xyz", 0, 1, -265.46482),
        ("20_achtar10.xyz", 0, 1, -356.28265),
        ("21_acanil01.xyz", 0, 1, -432.03012),
        ("22_benzidine.xyz", 0, 1, -563.27798),
        ("23_pterin.xyz", 0, 1, -569.84884),
        ("24_difuropyrazine.xyz", 0, 1, -556.71910),
        ("25_mesityloxide.xyz", 0, 1, -304.05919),
        ("26_histidine.xyz", 0, 1, -538.54910),
        ("27_dimethylpentane.xyz", 0, 1, -271.20088),
        ("28_caffeine.xyz", 0, 1, -667.73565),
        ("29_menthone.xyz", 0, 1, -458.44639),
    )
    prefix = "lib:baker/"
    return prefix, data


def get_baker_ts_data():
    """10.1002/(SICI)1096-987X(199605)17:7<888::AID-JCC12>3.0.CO;2-7

    HF/3-21G
    """
    data = (
        ("01_hcn.xyz", 0, 1, -92.24604),
        ("02_hcch.xyz", 0, 1, -76.29343),
        ("03_h2co.xyz", 0, 1, -113.05003),
        ("04_ch3o.xyz", 0, 2, -113.69365),
        ("05_cyclopropyl.xyz", 0, 2, -115.72100),
        ("06_bicyclobutane.xyz", 0, 1, -153.90494),
        ("07_bicyclobutane.xyz", 0, 1, -153.89754),
        ("08_formyloxyethyl.xyz", 0, 2, -264.64757),
        ("09_parentdieslalder.xyz", 0, 1, -231.60321),
        # 10 and 11 don't have any imaginary frequencies at the given
        # geometry, so they may be skipped. Until now (jan. 2021) I've
        # never seen anybody mention, how they treated these cases ...
        # It's not discussed in the 2002 Bakken paper and also not
        # in the 1998 Besalu/Bofill paper. And of course not in the
        # original 1996 Baker paper ... any advice on these two cases
        # is greatly appreciated.
        ("10_tetrazine.xyz", 0, 1, -292.81026),
        ("11_trans_butadiene.xyz", 0, 1, -154.05046),
        ("12_ethane_h2_abstraction.xyz", 0, 1, -78.54323),
        ("13_hf_abstraction.xyz", 0, 1, -176.98453),
        ("14_vinyl_alcohol.xyz", 0, 1, -151.91310),
        # 15 does not have an imaginary mode in cartesian coordinates
        ("15_hocl.xyz", 0, 1, -569.897524),
        ("16_h2po4_anion.xyz", -1, 1, -637.92388),
        ("17_claisen.xyz", 0, 1, -267.23859),
        ("18_silyene_insertion.xyz", 0, 1, -367.20778),
        ("19_hnccs.xyz", 0, 1, -525.43040),
        # The energy given in the paper (-168.24752 au) is the correct one
        # if one forms the central (0,1) bond (0-based indexing). If this
        # bond is missing, as it is if we autogenerate with bond-factor=1.3
        # then a TS with -168.241392 will be found.
        # For now we will use the original value from the paper.
        ("20_hconh3_cation.xyz", 1, 1, -168.24752),
        ("21_acrolein_rot.xyz", 0, 1, -189.67574),
        # The published energy -242.25529 corresponds to a planar TS. Without
        # symmetry restrictions the planar TS relaxes to -242.25695785.
        # As our algorithms obtain the unconstrained TS, we will use the
        # updated energy.
        # ("22_hconhoh.xyz", 0, 1, -242.25529),
        ("22_hconhoh.xyz", 0, 1, -242.256958),
        ("23_hcn_h2.xyz", 0, 1, -93.31114),
        ("24_h2cnh.xyz", 0, 1, -93.33296),
        ("25_hcnh2.xyz", 0, 1, -93.28172),
    )
    prefix = "lib:baker_ts/"
    return prefix, data


def get_s22_data():
    """
    https://doi.org/10.1039/B600027D
    """
    data = (
        ("00_adenine_thymine_wc.xyz", 0, 1, None),
        ("01_adenine_thymine_stack.xyz", 0, 1, None),
        ("02_ammonia_dimer.xyz", 0, 1, None),
        ("03_water_dimer.xyz", 0, 1, None),
        ("04_methane_dimer.xyz", 0, 1, None),
        ("05_ethene_dimer.xyz", 0, 1, None),
        ("06_ethene_ethine.xyz", 0, 1, None),
        ("07_formic_acid_dimer.xyz", 0, 1, None),
        ("08_formamide_dimer.xyz", 0, 1, None),
        ("09_benzene_water.xyz", 0, 1, None),
        ("10_benzene_ammonia.xyz", 0, 1, None),
        ("11_benzene_methane.xyz", 0, 1, None),
        ("12_benzene_dimer_c2v.xyz", 0, 1, None),
        ("13_benzene_dimer_c2h.xyz", 0, 1, None),
        ("14_indole_benzene_t-shape.xyz", 0, 1, None),
        ("15_indole_benzene_stack.xyz", 0, 1, None),
        ("16_pyrazine_dimer.xyz", 0, 1, None),
        ("17_2-pyridoxine_2-aminopyridine.xyz", 0, 1, None),
        ("18_phenol_dimer.xyz", 0, 1, None),
        ("19_uracil_dimer_stack.xyz", 0, 1, None),
        ("20_uracil_dimer_hb.xyz", 0, 1, None),
        ("21_benzene_hcn.xyz", 0, 1, None),
    )
    prefix = "lib:s22/"
    return prefix, data


def get_zimmerman_data():
    """
    https://dx.doi.org/10.1021/ct400319w
    """
    size = 105
    prefix = "lib:zimmerman/"
    data = list()
    for i in range(size):
        fn = f"case_{i:03d}.trj"
        data.append((fn, 0, 1, None))
    data = tuple(data)
    return prefix, data


def get_zimmerman_xtb_data():
    """
    Reoptimization of
        https://dx.doi.org/10.1021/ct400319w
    at the gfn2-xtb level of theory.

    Includes set 1 (first 72 entries), minus (0-based) ids
        7, 8, 13, 14, 15, 34, 39
    """
    size = 65
    prefix = "lib:zimmerman_xtb/"
    data = list()
    for i in range(size):
        fn = f"{i:02d}_zm_xtb.trj"
        data.append((fn, 0, 1, None))
    data = tuple(data)
    return prefix, data


def get_birkholz_rx_data():
    """
    https://doi.org/10.1002/jcc.23910
    """
    data = (
        ("00_c2no2.trj", 0, 1, None),
        ("01_c5ht.trj", 0, 1, None),
        ("02_hcn.trj", 0, 1, None),
        ("03_cope.trj", 0, 1, None),
        ("04_cpht.trj", 0, 1, None),
        ("05_cycbut.trj", 0, 1, None),
        ("06_dacp2.trj", 0, 1, None),
        ("07_dacp_eth.trj", 0, 1, None),
        ("08_dfcp.trj", 0, 1, None),
        ("09_ene.trj", 0, 1, None),
        ("10_grignard.trj", 0, 1, None),
        ("11_h2co.trj", 0, 1, None),
        ("12_hf_eth.trj", 0, 1, None),
        ("13_hydro.trj", 0, 1, None),
        ("14_meoh.trj", 0, 1, None),
        ("15_oxirane.trj", -1, 1, None),
        ("16_oxycope.trj", 0, 1, None),
        ("17_silane.trj", 0, 1, None),
        ("18_sn2.trj", -1, 1, None),
        ("19_sulfolene.trj", 0, 1, None),
    )
    prefix = "lib:birkholz_rx/"
    return prefix, data


def get_xtb_rx_data():
    data = (
        ("00_c2no2.trj", 0, 1, None),
        ("01_c5ht.trj", 0, 1, None),
        ("02_hcn.trj", 0, 1, None),
        ("03_cope.trj", 0, 1, None),
        ("04_cpht.trj", 0, 1, None),
        ("05_cycbut.trj", 0, 1, None),
        ("06_dacp2.trj", 0, 1, None),
        ("07_dacp_eth.trj", 0, 1, None),
        ("08_ene.trj", 0, 1, None),
        ("09_grignard.trj", 0, 1, None),
        ("10_h2co.trj", 0, 1, None),
        ("11_hf_eth.trj", 0, 1, None),
        ("12_hydro.trj", 0, 1, None),
        ("13_meoh.trj", 0, 1, None),
        ("14_oxirane.trj", -1, 1, None),
        ("15_oxycope.trj", 0, 1, None),
        ("16_silane.trj", 0, 1, None),
        ("17_sulfolene.trj", 0, 1, None),
        ("18_mobh35_14.trj", 0, 1, None),
        ("19_mobh35_30.trj", 0, 1, None),
    )
    prefix = "lib:xtb_rx/"
    return prefix, data


def get_precon_pos_rot_data():
    """
    https://doi.org/10.1002/jcc.23910
    """
    data = (
        ("00_c2no2.trj", 0, 1, None),
        ("06_dacp2.trj", 0, 1, None),
        ("07_dacp_eth.trj", 0, 1, None),
        ("08_dfcp.trj", 0, 1, None),
        ("09_ene.trj", 0, 1, None),
        ("11_h2co.trj", 0, 1, None),
        ("12_hf_eth.trj", 0, 1, None),
        ("14_meoh.trj", 0, 1, None),
        ("18_sn2.trj", -1, 1, None),
        ("19_sulfolene.trj", 0, 1, None),
    )
    prefix = "lib:birkholz_rx/"
    return prefix, data
