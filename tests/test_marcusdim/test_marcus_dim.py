from pathlib import Path
import tempfile

import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.constants import AU2EV
from pysisyphus.io import geom_from_hessian
from pysisyphus.init_logging import init_logging
from pysisyphus.marcus_dim import run_marcus_dim
import pysisyphus.marcus_dim.types as mdtypes
from pysisyphus.testing import using

_1EV = 1.0 / AU2EV
init_logging()


class AllEnsXTB(XTB):
    def get_all_energies(self, atoms, coords, **prepare_kwargs):
        en_res = self.get_energy(atoms, coords, **prepare_kwargs)
        energy = en_res["energy"]
        dummy_en = energy + _1EV
        result = {
            "energy": energy,
            "all_energies": np.array((energy, dummy_en)),
        }
        return result

    def get_wavefunction(self, atoms, coords, **prepare_kwargs):
        result = self.get_energy(atoms, coords, **prepare_kwargs)
        result["wavefunction"] = self.get_stored_wavefunction()
        return result


@using("xtb")
@pytest.mark.parametrize("cluster", (True, False))
def test_run_marcus_dim(cluster, this_dir):
    use_tmpdir = True
    # use_tmpdir = False
    if use_tmpdir:
        out_dir = tempfile.mkdtemp()
    else:
        out_dir = f"cluster_{str(cluster).lower()}"
    out_dir = Path(out_dir)

    pal = 8

    def calc_getter(with_td=False, **calc_kwargs):
        _calc_kwargs = {
            "charge": -1,
            "mult": 2,
            "pal": pal,
            "mem": 2000,
            "wavefunction_dump": True,
            "out_dir": out_dir / "qm_calcs",
        }
        if with_td:
            pass
        _calc_kwargs.update(calc_kwargs)
        calc = AllEnsXTB(**_calc_kwargs)
        return calc

    h5_fn = "pdnb_xtb_final_hessian.h5"
    calc = calc_getter(base_name="eq", pal=pal)
    geom = geom_from_hessian(h5_fn, calculator=calc)

    fragments = [
        [0, 1, 2],  # O N O
        [3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        [13, 14, 15],  # O N O
    ]
    T = 300.0
    seed = 20182020
    calc_kwargs = {
        "fragments": fragments,
        "temperature": T,
        "property": mdtypes.Property.EPOS_IAO,
        "seed": seed,
    }
    batch_kwargs = {
        "rms_thresh": 0.005,
        "batch_size": 32,
        "max_batches": 20,
    }
    marcus_dim_data, scan_data, model_data = run_marcus_dim(
        geom,
        calc_getter=calc_getter,
        calc_kwargs=calc_kwargs,
        batch_kwargs=batch_kwargs,
        rd_class=mdtypes.RobinDay.CLASS3,
        cluster=cluster,
        out_dir=out_dir,
    )
    if use_tmpdir:
        print(f"Calculations ran in '{out_dir}'")

    nu_marcus = marcus_dim_data["nu_marcus"]
    mass_marcus = marcus_dim_data["mass_marcus"]
    model_a = model_data["a"]
    coupling_a = model_a.coupling
    np.testing.assert_allclose(
        (nu_marcus, mass_marcus, coupling_a),
        (1173.12, 15.42, 0.018),
        atol=1e-2,
    )
    # assert coupling_a == pytest.approx(0.018, abs=1e-2)
