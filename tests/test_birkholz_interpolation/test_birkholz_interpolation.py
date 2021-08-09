import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.drivers import birkholz_interpolation
from pysisyphus.helpers import geom_loader, do_final_hessian
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.tsoptimizers import RSPRFOptimizer
from pysisyphus.testing import using


@using("pyscf")
@pytest.mark.parametrize(
    "fn, ref_energy",
    [
        ("lib:birkholz_rx/02_hcn_original.trj", -91.56485102),
        # ("lib:birkholz_rx/03_cope.trj", -230.06026151),
    ],
)
def test_birkholz_interpolation(fn, ref_energy):
    geoms = geom_loader(fn)

    def calc_getter():
        return PySCF(basis="sto3g", verbose=0)

    ts_guess = birkholz_interpolation(geoms, calc_getter)

    print(highlight_text("TS-Optimization"))
    tsopt_kwargs = {
        "dump": True,
        "thresh": "gau",
        "trust_max": 0.3,
    }
    tsopt = RSPRFOptimizer(ts_guess, **tsopt_kwargs)
    tsopt.run()

    do_final_hessian(ts_guess, write_imag_modes=True)
    assert ts_guess.energy == pytest.approx(ref_energy)
