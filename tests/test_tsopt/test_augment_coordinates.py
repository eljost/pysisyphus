#!/usr/bin/env python3

import pytest

from pysisyphus.calculators import Gaussian16
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.augment_bonds import augment_bonds
from pysisyphus.testing import using
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer


@using("pyscf")
@pytest.mark.parametrize(
    "augment, ref_cycle",
    [
        (True, 6),
        # 57 cycles are needed when 45 < VALID_BEND_DEG < 170
        # (False, 57),
        (False, 7),
    ],
)
def test_augment_coordinates_silyl(augment, ref_cycle):
    geom = geom_loader("lib:baker_ts/18_silyene_insertion.xyz", coord_type="redund")

    opt_kwargs = {
        "thresh": "baker",
        "max_cycles": 100,
        "dump": True,
        "trust_radius": 0.3,
        "trust_max": 0.3,
        "augment_bonds": augment,
    }
    calc_kwargs = {
        "charge": 0,
        "mult": 1,
        "pal": 4,
    }
    calc = PySCF(basis="321g", **calc_kwargs)

    geom.set_calculator(calc)

    opt = RSIRFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle

    ref_en = -367.20778
    assert opt.geometry.energy == pytest.approx(ref_en)


@using("gaussian16")
@pytest.mark.parametrize(
    "augment, ref_cycle",
    [
        (True, 28),
        (False, 40),
    ],
)
def test_augment_biaryl_bare(augment, ref_cycle):
    geom = geom_loader("lib:biaryl_bare_pm6_splined_hei.xyz", coord_type="redund")

    calc = Gaussian16("PM6", pal=4)
    geom.set_calculator(calc)

    opt_kwargs = {
        "thresh": "gau_tight",
        "augment_bonds": augment,
    }
    opt = RSIRFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
