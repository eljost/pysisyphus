from time import time

import numpy as np

from pysisyphus.config import T_DEFAULT, p_DEFAULT
from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.helpers_pure import highlight_text, standard_state_corr
from pysisyphus.thermo import print_thermoanalysis


def do_endopt_ts_barriers(
    ts_geom,
    left_geoms,
    right_geoms=None,
    left_fns=None,
    right_fns=None,
    do_thermo=False,
    T=T_DEFAULT,
    p=p_DEFAULT,
    solv_calc_getter=None,
    do_standard_state_corr=False,
):
    print(highlight_text("Barrier heights after end optimization(s)"))
    print()

    # Only allow standard state correction when solvent calculator is specified.
    do_ssc = do_standard_state_corr and (solv_calc_getter is not None)
    ssc = standard_state_corr(T=T, p=p) if do_ssc else 0.0

    if right_geoms is None:
        right_geoms = []

    if left_fns is None:
        left_fns = ["" for geom in left_geoms]
    if right_fns is None:
        right_fns = ["" for geom in right_geoms]

    ts_atom_num = len(ts_geom.atoms)

    def drop_total_geom(geoms):
        atom_nums = [len(geom.atoms) for geom in geoms]
        tot_atom_num = sum(atom_nums)
        total_geom = None
        if tot_atom_num > ts_atom_num:
            total_ind = atom_nums.index(max(atom_nums))
            total_geom = geoms[total_ind]
            geoms = [geom for i, geom in enumerate(geoms) if i != total_ind]
        return geoms, total_geom

    left_geoms, left_total_geom = drop_total_geom(left_geoms)
    right_geoms, right_total_geom = drop_total_geom(right_geoms)

    def tot_atom_num(geoms):
        return sum([len(geom.atoms) for geom in geoms])

    assert tot_atom_num(left_geoms) == ts_atom_num
    right_tot_atom_num = tot_atom_num(right_geoms)
    assert (right_tot_atom_num == ts_atom_num) or (right_tot_atom_num == 0)

    # Gibbs free energies
    if do_thermo:
        en_key = "free energy"

        def get_thermo(geom, title):
            thermo = geom.get_thermoanalysis(geom, T=T)
            print_thermoanalysis(thermo, geom=geom, level=1, title=title)
            print()
            return thermo

        ts_thermo = get_thermo(ts_geom, "TS")
        ts_energy = ts_thermo.G
        left_thermos = [get_thermo(geom, fn) for geom, fn in zip(left_geoms, left_fns)]
        left_energies = [thermo.G for thermo in left_thermos]
        right_thermos = [
            get_thermo(geom, fn) for geom, fn in zip(right_geoms, right_fns)
        ]
        right_energies = [thermo.G for thermo in right_thermos]
    # Electronic energies only
    else:
        en_key = "energy"

        ts_energy = ts_geom.energy
        left_energies = [geom.energy for geom in left_geoms]
        right_energies = [geom.energy for geom in right_geoms]
    print()

    def get_solv_correction(geom):
        atom_num = len(geom.atoms)
        print(f"Solvated energy calculation for {atom_num} atoms ...", end="")
        start = time()
        solv_energy = solv_calc_getter().get_energy(geom.atoms, geom.cart_coords)[
            "energy"
        ]
        dur = time() - start
        print(f" finished in {dur:.1f} s.")
        # Add standard state correction. ssc will be 0.0 if, disabled.
        return (solv_energy - geom.energy) + ssc

    if solv_calc_getter:
        ts_solv_corr = get_solv_correction(ts_geom)
        left_solv_corrs = [get_solv_correction(geom) for geom in left_geoms]
        right_solv_corrs = [get_solv_correction(geom) for geom in right_geoms]
        left_solv_corr = sum(left_solv_corrs)
        right_solv_corr = sum(right_solv_corrs)
        print()
    else:
        ts_solv_corr = left_solv_corr = right_solv_corr = 0.0

    ts_energy += ts_solv_corr
    left_energy = sum(left_energies) + left_solv_corr
    right_energy = sum(right_energies) + right_solv_corr
    energies = [left_energy, ts_energy]
    # TS is always only 1 geometry
    geom_nums = [len(left_geoms), 1]

    fns = ["Left", "TS"]
    if right_geoms:
        energies.append(right_energy)
        geom_nums.append(len(right_geoms))
        fns.append("Right")
    max_len = max(len(s) for s in fns)

    energies = np.array(energies)
    energies -= energies.min()
    min_ind = energies.argmin()
    energies *= AU2KJPERMOL

    print(highlight_text("Barriers", level=1))
    print()
    print(f"Temperature: {T:.2f} K")
    print(f"Pressure: {p:.1f} Pa\n")

    print("Corrections:")
    print(f"Change of standard-state: {do_ssc}, {ssc*AU2KJPERMOL:.2} kJ mol⁻¹")
    print(f"                 Solvent: {solv_calc_getter is not None}")
    print(f"         Thermochemistry: {do_thermo}")
    print()


    def print_geoms_fns(geoms, fns):
        for i, (geom, fn) in enumerate(zip(geoms, fns)):
            print(f"\t{i}: {fn} ({geom}, {len(geom.atoms)} atoms)")

    def get_geom_key(geoms): return "geometry" if len(geoms) == 1 else "geometries"
    print(f"Left {get_geom_key(left_geoms)}:")
    print_geoms_fns(left_geoms, left_fns)
    print("TS geometry:")
    print_geoms_fns((ts_geom, ), ("", ))
    if right_geoms:
        print(f"Right {get_geom_key(right_geoms)}:")
        print_geoms_fns(right_geoms, right_fns)
    print()

    print(f"Minimum {en_key} of {energies[min_ind]} kJ mol⁻¹ at '{fns[min_ind]}'.")
    print()
    for fn, en, gnum in zip(fns, energies, geom_nums):
        geom_str = "geometry" if gnum == 1 else "geometries"
        is_sum = " " if gnum == 1 else "Σ"
        print(f"\t{is_sum}{fn:>{max_len}s}: {en:>8.2f} kJ mol⁻¹ ({gnum} {geom_str})")
    print()
