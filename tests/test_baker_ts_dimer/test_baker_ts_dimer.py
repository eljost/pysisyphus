from pathlib import Path
import os

from natsort import natsorted
import numpy as np
import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators import Dimer
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.testing import using


def make_N_init_dict():
    THIS_DIR = Path(os.path.abspath(os.path.dirname(__file__)))
    xyz_path = THIS_DIR / "../../xyz_files/baker_ts"
    xyzs = natsorted(xyz_path.glob("*.xyz"))
    N_dict = dict()
    for guess, initial in [xyzs[2 * i : 2 * i + 2] for i in range(25)]:
        assert "downhill" in initial.stem
        assert guess.stem[:2] == initial.stem[:2]
        guess_geom = geom_from_xyz_file(guess)
        initial_geom = geom_from_xyz_file(initial)
        N_init = guess_geom.coords - initial_geom.coords
        N_dict[guess.name] = N_init
    return N_dict


# Create reasonable initial dimer orientations
N_INITS = make_N_init_dict()
BakerTSBm = Benchmark(
    "baker_ts", coord_type="cart", exclude=(9, 10, 14)
)


@using("pyscf")
@pytest.mark.parametrize("fn, geom, charge, mult, ref_energy", BakerTSBm.geom_iter)
def test_baker_ts_dimer(fn, geom, charge, mult, ref_energy):
    calc_kwargs = {
        "charge": charge,
        "mult": mult,
        "pal": 2,
        "verbose": 0,
        "base_name": Path(fn).stem,
    }
    calc = PySCF("321g", **calc_kwargs)

    dimer_kwargs = {
        "rotation_method": "fourier",
        "calculator": calc,
        "N_raw": N_INITS[fn],
        "length": 0.0189,
        "rotation_tol": 5,
    }
    dimer = Dimer(**dimer_kwargs)
    geom.set_calculator(dimer)

    opt_kwargs = {
        "thresh": "baker",
        "precon": True,
        "max_step_element": 0.25,
        "max_cycles": 50,
        "c_stab": 0.103,
    }
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(ref_energy)

    print(f"@{fn} converged using {dimer.force_evals} force evaluations")
