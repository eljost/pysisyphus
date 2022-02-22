from pathlib import Path
from time import time

import numpy as np

from pysisyphus.config import T_DEFAULT, p_DEFAULT
from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.helpers_pure import highlight_text, standard_state_corr
from pysisyphus.TablePrinter import TablePrinter
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
    calc_getter=None,
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
    ts_fn = "TS"

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

    left_atom_num = tot_atom_num(left_geoms)
    assert left_atom_num == ts_atom_num, (
        f"Atom number mismatch between TS ({ts_atom_num} atoms) and "
        f"left side ({left_atom_num} atoms)! Aborting barrier calculation"
    )
    right_tot_atom_num = tot_atom_num(right_geoms)
    assert (right_tot_atom_num == ts_atom_num) or (right_tot_atom_num == 0)

    def get_energy(geom, base_name):
        try:
            geom.energy
        except AttributeError:
            try:
                geom.set_calculator(calc_getter(base_name=base_name))
            except TypeError:
                raise Exception(
                    f"Energy isn't set at '{base_name}' geometry "
                    "and no 'calc_getter' was supplied!\nPlease either set "
                    "the desired quantities (energy/(Hessian)) or supply "
                    "a 'calc_getter' to calculate them."
                )
        return geom.energy

    def get_energies(geoms, base_name):
        return [
            get_energy(geom, f"{base_name}_{i:02d}") for i, geom in enumerate(geoms)
        ]

    en_key = "energy"
    ts_energy = get_energy(ts_geom, "ts")
    left_energies = get_energies(left_geoms, "left")
    right_energies = get_energies(right_geoms, "right")

    def zeros(geoms):
        return np.zeros(len(geoms))

    # Gibbs free energies
    if do_thermo:
        en_key = "free energy"

        def get_thermo(geom, title, is_ts=False):
            thermo = geom.get_thermoanalysis(geom, T=T, p=p)
            print_thermoanalysis(thermo, geom=geom, is_ts=is_ts, level=1, title=title)
            print()
            return thermo

        ts_thermo = get_thermo(ts_geom, "TS", is_ts=True)
        ts_dG = ts_thermo.dG
        left_thermos = [get_thermo(geom, fn) for geom, fn in zip(left_geoms, left_fns)]
        left_dGs = [thermo.dG for thermo in left_thermos]
        right_thermos = [
            get_thermo(geom, fn) for geom, fn in zip(right_geoms, right_fns)
        ]
        right_dGs = [thermo.dG for thermo in right_thermos]
    else:
        ts_dG = 0.0
        left_dGs = zeros(left_geoms)
        right_dGs = zeros(right_geoms)

    def get_solv_correction(geom, fn=None, name=None):
        fn = str(fn)
        if fn is not None:
            infix = f"'{Path(fn).name: >40s}'"
        else:
            atom_num = len(geom.atoms)
            infix = f"{atom_num} atoms"

        print(f"Solvated calculation for {infix} ...", end="")
        start = time()
        solv_calc_kwargs = {}
        if name is not None:
            solv_calc_kwargs["base_name"] = f"solv_{name}"
        solv_energy = solv_calc_getter(**solv_calc_kwargs).get_energy(
            geom.atoms, geom.cart_coords
        )["energy"]
        dur = time() - start
        print(f" finished in {dur:.1f} s.")
        # Add standard state correction. ssc will be 0.0 if, disabled.
        return solv_energy, (solv_energy - geom.energy) + ssc

    def get_solv_corrections(geoms, fns, name):
        return zip(
            *[get_solv_correction(geom, fn, name) for geom, fn in zip(geoms, fns)]
        )

    if solv_calc_getter is not None:
        print(highlight_text("Calculations with solvent model", level=1) + "\n")
        ts_solv_energy, ts_solv_corr = get_solv_correction(ts_geom, ts_fn, "ts")
        left_solv_energies, left_solv_corrs = get_solv_corrections(
            left_geoms, left_fns, "left"
        )
        right_solv_energies, right_solv_corrs = get_solv_corrections(
            right_geoms, right_fns, "right"
        )
        print()
    # Use zeros, so we can add the term later.
    else:
        ts_solv_corr = 0.0
        left_solv_corrs = zeros(left_geoms)
        right_solv_corrs = zeros(right_geoms)

    # Calculate total contributions for both sides
    left_dG = sum(left_dGs)
    right_dG = sum(right_dGs)
    left_solv_corr = sum(left_solv_corrs)
    right_solv_corr = sum(right_solv_corrs)

    ts_energy_corr = ts_energy + ts_solv_corr + ts_dG
    left_energy_corr = sum(left_energies) + left_dG + left_solv_corr
    right_energy_corr = sum(right_energies) + right_dG + right_solv_corr
    energies_corr = [left_energy_corr, ts_energy_corr]
    # TS is always only 1 geometry
    geom_nums = [len(left_geoms), 1]

    fns = ["Left", "TS"]
    if right_geoms:
        energies_corr.append(right_energy_corr)
        geom_nums.append(len(right_geoms))
        fns.append("Right")
    max_len = max(len(s) for s in fns)

    def print_table(table, header, width=15):
        col_fmts = ["str"] + (len(header) - 1) * ["float"]
        tp = TablePrinter(header, col_fmts, width=width, sub_underline=False)
        tp.print_header()
        for fn, row in zip(all_fns, table.T):
            fn = str(fn)  # fns may be PosixPath etc.
            fn_cut = fn[:6] + "..." + fn[-6:] if (len(fn) > width) else fn
            tp.print_row((fn_cut, *row))
        print()

    all_fns = left_fns + [ts_fn] + right_fns
    # Electronic energies
    all_energies = np.concatenate((left_energies, [ts_energy], right_energies))
    to_stack = [all_energies]
    gas_titles = ["File", "E_el"]
    if do_thermo:
        # Thermochemical corrections
        all_dGs = np.concatenate((left_dGs, [ts_dG], right_dGs))
        # Free Gibbs energies in the gas phase
        all_Gs_gas = all_energies + all_dGs
        to_stack += [all_dGs, all_Gs_gas]
        gas_titles += ["dG_gas", "G_gas"]

    gas_table = np.stack(to_stack)
    print("All tabulated quantities are given in au.\n")
    print_table(gas_table, gas_titles)
    full_titles = gas_titles

    if solv_calc_getter is not None:
        all_solv_energies = np.concatenate(
            (left_solv_energies, [ts_solv_energy], right_solv_energies)
        )
        all_solv_corrs = np.concatenate(
            (left_solv_corrs, [ts_solv_corr], right_solv_corrs)
        )
        to_stack = [all_energies, all_solv_energies, all_solv_corrs]
        solv_titles = ["File", "E_el", "E_solv", "dG_solv"]
        if do_thermo:
            all_Gs_sol = all_Gs_gas + all_solv_corrs
            to_stack += [all_Gs_sol]
            solv_titles += ["G_sol"]
        solv_table = np.stack(to_stack)
        if do_ssc:
            print("dG_solv includes a correction for change of standard state.\n")
        print_table(solv_table, solv_titles)
        full_titles += solv_titles[2:]

    energies_corr = np.array(energies_corr)
    energies_corr -= energies_corr.min()
    min_ind = energies_corr.argmin()
    energies_corr *= AU2KJPERMOL

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
            fn_name = Path(fn).name
            print(f"\t{i}: {fn_name} ({geom}, {len(geom.atoms)} atoms)")

    def get_geom_key(geoms):
        return "geometry" if len(geoms) == 1 else "geometries"

    print(f"Left {get_geom_key(left_geoms)}:")
    print_geoms_fns(left_geoms, left_fns)
    print("TS geometry:")
    print_geoms_fns((ts_geom,), ("",))
    if right_geoms:
        print(f"Right {get_geom_key(right_geoms)}:")
        print_geoms_fns(right_geoms, right_fns)
    print()

    print(f"Minimum {en_key} of {energies_corr[min_ind]} kJ mol⁻¹ at '{fns[min_ind]}'.")
    print()
    for fn, en, gnum in zip(fns, energies_corr, geom_nums):
        geom_str = "geometry" if gnum == 1 else "geometries"
        is_sum = " " if gnum == 1 else "Σ"
        print(f"\t{is_sum}{fn:>{max_len}s}: {en:>8.2f} kJ mol⁻¹ ({gnum} {geom_str})")
    print()
