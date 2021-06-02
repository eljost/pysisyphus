import numpy as np

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.thermo import get_thermoanalysis, print_thermoanalysis


def do_endopt_ts_barriers(
    ts_geom,
    left_geoms,
    right_geoms=None,
    left_fns=None,
    right_fns=None,
    do_thermo=False,
    T=298.15,
):
    print(highlight_text("Barrier heights after end optimization(s)"))
    print()

    if right_geoms is None:
        right_geoms = []
        right_fns = []

    # Gibbs free energies
    if do_thermo:
        print(f"Including thermochemical corrections at {T:.2f} K.\n")
        en_key = "free energy"

        def get_thermo(geom, title):
            thermo = get_thermoanalysis(geom, T=T)
            print_thermoanalysis(thermo, geom=geom, level=1, title=title)
            print()
            return thermo

        ts_thermo = get_thermo(ts_geom, "TS")
        ts_energy = ts_thermo.G
        left_thermos = [
            get_thermo(geom, fn) for geom, fn in zip(left_geoms, left_fns)
        ]
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

    left_energy = sum(left_energies)
    right_energy = sum(right_energies)
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


