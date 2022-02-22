from pathlib import Path
import shutil
import tempfile

import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.helpers import align_geoms
# from pysisyphus.helpers_pure import filter_fixture_store
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using
from pysisyphus.xyzloader import write_geoms_to_trj


Bh = Benchmark(
    "xtb_rx",
    # only=(18, 19),
)


@pytest.mark.benchmark
@using("xtb")
@pytest.mark.parametrize("fn, geoms, charge, mult, ref_energy", Bh.geom_iter)
def test_xtb_rx(fn, geoms, charge, mult, ref_energy, results_bag):
    start, ts_ref_org, end = geoms
    id_ = fn[:2]

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        inp_ts = str(tmp_path / "ts_input.trj")
        with open(inp_ts, "w") as handle:
            handle.write(ts_ref_org.as_xyz())
        ts_run_dict = {
            "geom": {
                "type": "redund",
                "fn": inp_ts,
            },
            "calc": {
                "type": "xtb",
                "pal": 6,
                "mem": 750,
                "charge": charge,
                "mult": mult,
                "quiet": True,
            },
            "tsopt": {
                "type": "rsirfo",
                "hessian_recalc": 1,
                "trust_max": 0.3,
                "thresh": "gau",
                "do_hess": True,
            },
        }
        ts_results = run_from_dict(ts_run_dict)
        # Reference values
        ts_ref = ts_results.ts_geom
        ts_ref_energy = ts_ref.energy
        ts_ref_imag = ts_ref.get_imag_frequencies()[0]

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        inp_trj = str(tmp_path / "gs_inputs.trj")
        write_geoms_to_trj((start, end), inp_trj)

        run_dict = {
            "geom": {
                "type": "dlc",
                "fn": inp_trj,
            },
            "calc": {
                "type": "xtb",
                "pal": 1,
                "mem": 750,
                "charge": charge,
                "mult": mult,
            },
            "preopt": {
                "max_cycles": 5,
            },
            "cos": {
                "type": "gs",
                "climb": True,
                "climb_rms": 0.0075,
                "reset_dlc": True,
            },
            "opt": {
                "type": "string",
                "max_step": 0.2,
                "rms_force": 0.005,
                "rms_force_only": True,
            },
            "tsopt": {
                "type": "rsirfo",
                "thresh": "gau",
                "trust_max": 0.5,
                "do_hess": True,
            },
        }
        if id_ == "02":
            run_dict["geom"]["type"] = "cart"
        elif id_ == "19":
            run_dict["opt"]["rms_force"] = 0.003

        results = run_from_dict(run_dict)
        ts_geom = results.ts_geom
        ts_energy = ts_geom.energy
        ts_imag = ts_geom.get_imag_frequencies()[0]

        opt = results.opt
        ts_opt = results.ts_opt

        rmsd = ts_ref.rmsd(ts_geom)
        diff = ts_ref_energy - ts_energy
        cmt = "Ref" if diff < 0.0 else " TS"
        rmsd_fmt = " >12.4f"
        print(f"RMSD: {rmsd:{rmsd_fmt}}")
        print(f" TS energy: {ts_energy:.6f}")
        print(f"Ref energy: {ts_ref_energy:.6f}")
        print(f"      Diff: {diff:.6f}")
        print(
            f"@@@{id_} COMPARE@@@: rmsd={rmsd:{rmsd_fmt}}, ΔE= {diff: .6f} {cmt} is lower, "
            # f"Ref: {ts_ref_imag: >8.1f}, TS: {ts_imag: >8.1f} cm⁻¹"
            f"fn={fn[:10]}, cycs: opt={opt.cur_cycle+1: >2d}, tsopt={ts_opt.cur_cycle+1: >2d}"
        )

        results_bag.opt_converged = opt.is_converged
        results_bag.opt_cycles = opt.cur_cycle + 1
        # results_bag.tsopt_converged = ts_opt.is_converged
        results_bag.tsopt_cycles = ts_opt.cur_cycle + 1
        results_bag.rmsd = rmsd

        assert results.ts_opt.is_converged
        shutil.copy("ts_opt.xyz", f"{id_}_ts_opt.xyz")

        # Dump TS geoms
        ts_ref_org.comment = "TS ref org"
        ts_ref.comment = "TS ref opt"
        ts_geom.comment = "TS opt from cos"
        ts_geoms = (ts_ref_org, ts_ref, ts_geom)
        align_geoms(ts_geoms)
        ts_fns = f"{id_}_ts_geoms.trj"
        write_geoms_to_trj(ts_geoms, ts_fns)


@pytest.mark.benchmark
@using("xtb")
# @filter_fixture_store("test_xtb_rx")
def test_xtb_rx_synthesis(fixture_store):
    for i, fix in enumerate(fixture_store):
        print(i, fix)

    tot_opt_cycles = 0
    opt_converged = 0
    tot_tsopt_cycles = 0
    # tsopt_converged = 0
    bags = fixture_store["results_bag"]
    for k, v in bags.items():
        if not k.startswith("test_xtb_rx"):
            continue
        print(k)
        opt_converged += 1 if v["opt_converged"] else 0
        for kk, vv in v.items():
            print("\t", kk, vv)
        if opt_converged:
            tot_opt_cycles += v["opt_cycles"]
            tot_tsopt_cycles += v["tsopt_cycles"]
    print(f"Opts converged: {opt_converged}/{len(bags)}")
    print(f"Total opt cycles: {tot_opt_cycles}")
    print(f"Total tsopt cycles: {tot_tsopt_cycles}")
