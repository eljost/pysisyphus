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
    do_ssc = do_standard_state_corr and solv_calc_getter
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
        if tot_atom_num > ts_atom_num:
            total_ind = atom_nums.index(max(atom_nums))
            geoms = [geom for i, geom in enumerate(geoms) if i != total_ind]
        return geoms

    right_geoms = drop_total_geom(right_geoms)
    left_geoms = drop_total_geom(left_geoms)

    def tot_atom_num(geoms):
        return sum([len(geom.atoms) for geom in geoms])

    assert tot_atom_num(left_geoms) == ts_atom_num
    right_tot_atom_num = tot_atom_num(right_geoms)
    assert (right_tot_atom_num == ts_atom_num) or (right_tot_atom_num == 0)

    # Gibbs free energies
    if do_thermo:
        print(f"Including thermochemical corrections at {T:.2f} K.\n")
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
        print("Thermochemical corrections are NOT included!")
        en_key = "energy"

        ts_energy = ts_geom.energy
        left_energies = [geom.energy for geom in left_geoms]
        right_energies = [geom.energy for geom in right_geoms]
    print()

    def get_solv_correction(geom):
        solv_energy = solv_calc_getter().get_energy(geom.atoms, geom.cart_coords)[
            "energy"
        ]
        # Add standard state correction. ssc will be 0.0 if, disabled.
        return (solv_energy - geom.energy) + ssc

    if solv_calc_getter:
        ts_solv_corr = get_solv_correction(ts_geom)
        left_solv_corrs = [get_solv_correction(geom) for geom in left_geoms]
        right_solv_corrs = [get_solv_correction(geom) for geom in right_geoms]
        left_solv_corr = sum(left_solv_corrs)
        right_solv_corr = sum(right_solv_corrs)
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
    if do_ssc:
        print(
            f"Standard-state correction (gas->solution): {ssc*AU2KJPERMOL:.2} kJ mol⁻¹"
        )
    if solv_calc_getter is not None:
        print("Including solvent correction (dG_solv = E_solv - E_gas)")
    print()


    def print_geoms_fns(geoms, fns):
        for i, (geom, fn) in enumerate(zip(geoms, fns)):
            print(f"\t{i}: {fn} ({geom}, {len(geom.atoms)} atoms)")

    print("Left geometries:")
    print_geoms_fns(left_geoms, left_fns)
    if right_geoms:
        print("Right geometries:")
        print_geoms_fns(right_geoms, right_fns)
    print()

    print(f"Minimum {en_key} of {energies[min_ind]} kJ mol⁻¹ at '{fns[min_ind]}'.")
    print()
    for fn, en, gnum in zip(fns, energies, geom_nums):
        geom_str = "geometry" if gnum == 1 else "geometries"
        is_sum = " " if gnum == 1 else "Σ"
        print(f"\t{is_sum}{fn:>{max_len}s}: {en:>8.2f} kJ mol⁻¹ ({gnum} {geom_str})")
    print()
