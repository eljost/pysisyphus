from pathlib import Path
import shutil
import tempfile

import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.helpers import align_geoms
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using
from pysisyphus.xyzloader import write_geoms_to_trj


Bh = Benchmark(
    "birkholz_rx",
    # exclude=list(range(14)),
    # 03 cope; GSM can't handle this
    # 16 cope; GSM can't handle this
    # 18 has no TS at the GFN2-XTB level of theory
    exclude=(3, 16, 18),
)


@pytest.mark.benchmark
@using("xtb")
@pytest.mark.parametrize("fn, geoms, charge, mult, ref_energy", Bh.geom_iter)
def test_birkholz_rx_gsm(fn, geoms, charge, mult, ref_energy, results_bag):
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
                "type": "cart" if id_ == "02" else "dlc",
                "fn": inp_trj,
            },
            "calc": {
                "type": "xtb",
                "pal": 6,
                "mem": 750,
                "charge": charge,
                "mult": mult,
            },
            "preopt": {
                "max_cycles": 5,
            },
            "cos": {
                "type": "gs",
                # "fix_ends": True,
                # "max_nodes": 11,
                # "reparam_check": "rms",
                # "perp_thresh": 0.075,
                "climb": True,
                "climb_rms": 0.01,
                # "climb_lanczos": True,
                # "climb_lanczos_rms": 0.0075,
                # "reset_dlc": True,
            },
            "opt": {
                "type": "string",
                "max_step": 0.2,
                # "lbfgs_when_full": True,
                # "max_step": 0.25,
                # "keep_last": 10,
                "rms_force": 0.005,
                "rms_force_only": True,
                # "double_damp": True,
            },
            "tsopt": {
                "type": "rsirfo",
                "do_hess": True,
                "thresh": "gau",
                "trust_max": 0.3,
                "max_cycles": 100,
            },
        }

        results = run_from_dict(run_dict)
        ts_geom = results.ts_geom
        ts_energy = ts_geom.energy
        ts_imag = ts_geom.get_imag_frequencies()[0]

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
            f"Ref: {ts_ref_imag: >8.1f}, TS: {ts_imag: >8.1f} cm⁻¹"
        )

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
