import pytest

from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using
from pysisyphus.run import run_from_dict, do_endopt_ts_barriers


@using("xtb")
def test_run_dict_endopt_barriers():
    run_dict = {
        "geom": {
            "type": "redund",
            "fn": "lib:hcn_iso_ts_opt_xtb.xyz",
        },
        "tsopt": {
            "type": "rsirfo",
        },
        "irc": {
            "type": "eulerpc",
            "rms_grad_thresh": 0.015,
        },
        "endopt": {
            "thresh": "gau_loose",
        },
        "calc": {
            "type": "xtb",
            "pal": 1,
        },
    }
    results = run_from_dict(run_dict)

    assert results.ts_opt.is_converged
    assert results.ts_opt.cur_cycle == 0
    assert results.ts_geom.energy == pytest.approx(-5.38737424)


import numpy as np

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.helpers import highlight_text
from pysisyphus.thermo import get_thermoanalysis, print_thermoanalysis


def do_endopt_ts_barriers(
    ts_geom,
    left_geoms,
    right_geoms=None,
    do_thermo=False,
    T=298.15,
):
    print(highlight_text("Barrier heights after end optimizations"))

    if right_geoms is None:
        right_geoms = []

    # Gibbs free energies
    if do_thermo:
        print(f"\nIncluding thermochemical corrections at {T:.2f} K.\n")
        en_key = "free energy"

        def get_thermo(geom, title):
            thermo = get_thermoanalysis(geom, T=T)
            print_thermoanalysis(thermo, geom=geom, level=1, title=title)
            print()
            return thermo

        ts_thermo = get_thermo(ts_geom, "TS")
        ts_energy = ts_thermo.G
        left_thermos = [
            get_thermo(geom, f"forward_{i:03d}") for i, geom in enumerate(left_geoms)
        ]
        left_energies = [thermo.G for thermo in left_thermos]
        right_thermos = [
            get_thermo(geom, f"backward_{i:03d}") for i, geom in enumerate(right_geoms)
        ]
        right_energies = [thermo.G for thermo in right_thermos]
    # Electronic energies only
    else:
        print("Thermochemical corrections are NOT included!")
        en_key = "energy"

        ts_energy = ts_geom.energy
        left_energies = [geom.energy for geom in left_geoms]
        right_energies = [geom.energy for geom in right_geoms]

    left_energy = sum(left_energies)
    right_energy = sum(right_energies)
    energies = [left_energy, ts_energy]
    # TS is always only 1 geometry
    geom_nums = [len(left_geoms), 1]
    if right_geoms:
        energies.append(right_energy)
        geom_nums.append(len(right_geoms))
        fns = ("Forward", "TS", "Backward")
    else:
        fns = ("Downhill", "TS")
    max_len = max(len(s) for s in fns)
    energies = np.array(energies)
    energies -= energies.min()
    min_ind = energies.argmin()
    energies *= AU2KJPERMOL

    print(highlight_text("Barriers", level=1))
    print()
    print(f"Minimum {en_key} of {energies[min_ind]} kJ mol⁻¹ at '{fns[min_ind]}'.")
    print()
    for fn, en, gnum in zip(fns, energies, geom_nums):
        geom_str = "geometry" if gnum == 1 else "geometries"
        is_sum = "" if gnum == 1 else "Σ"
        print(f"\t{is_sum}{fn:>{max_len}s}: {en:>8.2f} kJ mol⁻¹ ({gnum} {geom_str})")
    print()


@using("thermoanalysis")
@using("xtb")
@pytest.mark.parametrize(
    "do_thermo, downhill",
    [(True, True), (True, False), (False, True), (False, False)],
)
def test_do_endopt_ts_barriers(this_dir, do_thermo, downhill):
    def get_geom(fn):
        geom = geom_loader(this_dir / fn)
        geom.set_calculator(XTB(pal=2))
        return geom

    ts_geom = get_geom("hcn_ts_opt.xyz")
    fw_fn = "hcn_forward_end_opt.xyz"
    bw_fn = "hcn_backward_end_opt.xyz"
    forward_geom = get_geom(fw_fn)
    backward_geom = get_geom(bw_fn)
    end_fns = (fw_fn, bw_fn)
    left_geoms = [forward_geom]
    if downhill:
        right_geoms = None
    else:
        right_geoms = [backward_geom]
    # do_endopt_ts_barriers(end_geoms, end_fns, ts_geom, do_thermo=do_thermo)
    print("@right", right_geoms)
    do_endopt_ts_barriers(ts_geom, left_geoms, right_geoms, do_thermo=do_thermo)
    # modify run_endopt in run.py in a way that it returns forwad_endopt
    # and backward_endopt separately
