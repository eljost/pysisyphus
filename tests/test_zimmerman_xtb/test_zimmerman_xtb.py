from pathlib import Path
import shutil
import tempfile

import pytest

# from pysisyphus.helpers import geom_loader, align_geoms
from pysisyphus.benchmarks import Benchmark
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using
from pysisyphus.xyzloader import write_geoms_to_trj


ZmXTB = Benchmark("zimmerman_xtb")


@pytest.mark.skip
@using("xtb")
@pytest.mark.parametrize("fn, geoms, ref_energy", ZmXTB)
def test_zimmerman_xtb_gsm(fn, geoms, ref_energy, results_bag):
    start, ts_ref, end = geoms
    id_ = fn[:2]

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
                "pal": 6,
                "quiet": True,
            },
            "cos": {
                "type": "gs",
                "fix_ends": True,
                "max_nodes": 11,
                "reparam_check": "rms",
                # "perp_thresh": 0.075,
                "climb": True,
                "climb_rms": 0.0075,
                "climb_lanczos": True,
                "climb_lanczos_rms": 0.0075,
                # "reset_dlc": True,
            },
            "opt": {
                "type": "string",
                "lbfgs_when_full": True,
                "max_step": 0.25,
                "keep_last": 10,
                "rms_force": 0.0035,
                "rms_force_only": True,
                "double_damp": True,
            },
            "tsopt": {
                "type": "rsirfo",
                "do_hess": True,
                "hessian_recalc": 3,
                # "thresh": "gau",
                "trust_max": 0.3,
                "max_cycles": 100,
                "root": 0,
            },
        }

        results = run_from_dict(run_dict)
        ts_geom = results.ts_geom
        ts_energy = ts_geom.energy
        ts_imag = ts_geom.get_imag_frequencies()[0]

        # Reference values
        ts_ref_results = ts_geom.calculator.get_hessian(ts_ref.atoms, ts_ref.cart_coords)
        ts_ref_energy = ts_ref_results["energy"]
        ts_ref._hessian = ts_ref_results["hessian"]
        ts_ref_imag = ts_ref.get_imag_frequencies()[0]

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
