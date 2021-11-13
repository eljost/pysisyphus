import numpy as np
import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


@using("pyscf")
def test_mulliken_charges():
    geom = geom_loader("lib:oniom_ee_model_system.xyz", coord_type="redund")

    calc_kwargs = {
        "basis": "sto3g",
    }
    calc = PySCF(**calc_kwargs)
    geom.set_calculator(calc)

    energy = geom.energy
    assert energy == pytest.approx(-5.789815189924002E+02)

    charges = calc.parse_charges()
    chrg_ref = """
        -4.84551410E-02  9.03447614E-02  8.75522347E-02 -4.66142813E-01  2.09615873E-01
      2.03308060E-01  3.18063759E-01 -3.52693455E-01  2.59066455E-01 -2.90419166E-01
     -4.29794447E-01  2.21955145E-01  1.77292603E-01 -3.98003060E-01  2.17833663E-01
      2.08840824E-01 -4.30825317E-01  1.78166353E-01  2.25524078E-01 -4.10146906E-01
      2.03058796E-01  2.25857702E-01"""
    ref_charges = np.array(chrg_ref.strip().split(), dtype=float)

    np.testing.assert_allclose(charges, ref_charges, atol=1e-5)
