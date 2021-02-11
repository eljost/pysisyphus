from pathlib import Path
import os

from natsort import natsorted
import numpy as np
import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators import Dimer, Gaussian16
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_xyz_file, geom_from_library, do_final_hessian
from pysisyphus.tsoptimizers.dimer import dimer_method
from pysisyphus.tsoptimizers.dimerv2 import dimer_method as dimer_method_v2
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



BakerTSBm = Benchmark(
    "baker_ts", coord_type="cart", exclude=(4, 9, 10, 14, 16, 17, 19)
)


@using("pyscf")
@pytest.mark.parametrize("fn, geom, charge, mult, ref_energy", BakerTSBm.geom_iter)
def test_baker_ts_dimer(fn, geom, charge, mult, ref_energy, results_bag):
    # Load initial dimers directions
    N_init_dict = make_N_init_dict()

    calc_kwargs = {
        "charge": charge,
        "mult": mult,
        "pal": 2,
        "base_name": Path(fn).stem,
    }

    def calc_getter():
        return PySCF(basis="321g", verbose=0, **calc_kwargs)

    geom.set_calculator(calc_getter())

    dimer_kwargs = {
        "max_step": 0.25,
        # 1e-2 Angstroem in bohr
        "dR_base": 0.0189,
        "rot_opt": "lbfgs",
        # "trans_opt": "mb",
        "trans_opt": "lbfgs",
        # "trans_memory": 10, # bad idea
        "angle_tol": 5,
        "f_thresh": 1e-3,
        "max_cycles": 50,
        "f_tran_mod": True,
        # "rot_type": "direct",
        "multiple_translations": True,
    }
    dimer_kwargs["N_init"] = N_init_dict[fn]
    geoms = (geom,)
    results = dimer_method(geoms, calc_getter, **dimer_kwargs)

    same_energy = geom.energy == pytest.approx(ref_energy)
    print(
        f"@Same energy: {str(same_energy): >5}, {fn}, "
        f"{results.force_evals} force evaluations"
    )
    if not same_energy:
        do_final_hessian(geom)

    # This way pytest prints the actual values... instead of just the boolean
    assert geom.energy == pytest.approx(ref_energy)
