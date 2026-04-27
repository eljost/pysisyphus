"""Test the N_hessian='calc' option that seeds N_raw from an analytic Hessian.

Kept in its own file to avoid importing optional heavy deps (pyscf) from the
sibling test_dimer.py module-level imports.
"""

import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.Dimer import Dimer
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS


def test_dimer_n_hessian_calc():
    """N_hessian='calc' seeds N_raw from the wrapped calculator's analytic Hessian.

    Starts near AnaPot's saddle point at (0.6117, 1.4930, 0) so the analytic
    Hessian has a clean imaginary mode. The dimer rotation converges from this
    analytic seed without any precomputed N_raw input.
    """
    coords = (0.55, 1.45, 0)
    geom = Geometry(("X",), coords)
    dimer = Dimer(
        calculator=AnaPot(),
        N_hessian="calc",
        rotation_remove_trans=False,
    )
    geom.set_calculator(dimer)

    opt = PreconLBFGS(
        geom,
        precon=False,
        line_search=None,
        max_step_element=0.25,
        thresh="gau_tight",
        max_cycles=15,
    )
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(2.80910484)
