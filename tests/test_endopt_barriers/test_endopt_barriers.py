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


def do_endopt_ts_barriers(end_geoms, end_fns, ts_geom, do_thermo=False, T=298.15):
    import numpy as np

    from pysisyphus.constants import AU2KJPERMOL
    from pysisyphus.helpers import highlight_text
    from pysisyphus.thermo import get_thermoanalysis, print_thermoanalysis

    print(highlight_text("Barrier heights after end optimizations"))

    forward_geoms, backward_geoms = end_geoms

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
        forward_thermos = [
            get_thermo(geom, f"forward_{i:03d}") for i, geom in enumerate(forward_geoms)
        ]
        forward_energies = [thermo.G for thermo in forward_thermos]
        backward_thermos = [
            get_thermo(geom, f"backward_{i:03d}")
            for i, geom in enumerate(backward_geoms)
        ]
        backward_energies = [thermo.G for thermo in backward_thermos]
    # Electronic energies only
    else:
        print("Thermochemical corrections are NOT included!")
        en_key = "energy"

        ts_energy = ts_geom.energy
        forward_geom, backward_geom = end_geoms
        forward_energies = [geom.energy for geom in forward_geoms]
        backward_energies = [geom.energy for geom in backward_geoms]

    forward_energy = sum(forward_energies)
    backward_energy = sum(backward_energies)
    energies = np.array((forward_energy, ts_energy, backward_energy))
    energies -= energies.min()
    min_ind = energies.argmin()
    energies *= AU2KJPERMOL

    forward_fn, backward_fn = "Forward", "Backward"
    fns = (forward_fn, "TS", backward_fn)
    max_len = max(len(s) for s in fns)

    print(highlight_text("Barriers", level=1))
    print()
    print(f"Minimum {en_key} of {energies[min_ind]} kJ mol⁻¹ at '{fns[min_ind]}'.")
    print()
    for fn, en in zip(fns, energies):
        print(f"\t{fn:>{max_len}s}: {en:>8.2f} kJ mol⁻¹")
        "Σ"
    print()


@using("thermoanalysis")
@using("xtb")
@pytest.mark.parametrize(
    "do_thermo",
    [
        True,
        False,
    ],
)
def test_do_endopt_ts_barriers(this_dir, do_thermo):
    def get_geom(fn):
        geom = geom_loader(this_dir / fn)
        geom.set_calculator(XTB(pal=2))
        return geom

    ts_geom = get_geom("hcn_ts_opt.xyz")
    fw_fn = "hcn_forward_end_opt.xyz"
    bw_fn = "hcn_backward_end_opt.xyz"
    forward_geom = get_geom(fw_fn)
    backward_geom = get_geom(bw_fn)
    end_geoms = [forward_geom], [backward_geom]
    end_fns = (fw_fn, bw_fn)
    do_endopt_ts_barriers(end_geoms, end_fns, ts_geom, do_thermo=do_thermo)
    # modify run_endopt in run.py in a way that it returns forwad_endopt
    # and backward_endopt separately
