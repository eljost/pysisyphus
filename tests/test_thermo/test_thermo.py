import h5py
import pytest

from pysisyphus.calculators import ORCA
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using
from pysisyphus.thermo import get_thermoanalysis, print_thermoanalysis


@using("thermoanalysis")
def test_thermoanalysis(this_dir):
    hess_fn = this_dir / "h2o_hessian.h5"
    with h5py.File(hess_fn, "r") as handle:
        masses = handle["masses"][:]
        vibfreqs = handle["vibfreqs"][:]
        coords3d = handle["coords3d"][:]
        energy = handle.attrs["energy"]
        mult = handle.attrs["mult"]

    thermo_dict = {
        "masses": masses,
        "vibfreqs": vibfreqs,
        "coords3d": coords3d,
        "energy": energy,
        "mult": mult,
    }

    from thermoanalysis.QCData import QCData
    from thermoanalysis.thermo import thermochemistry

    qcd = QCData(thermo_dict)
    thermo = thermochemistry(qcd, temperature=298.15)

    assert thermo.M == pytest.approx(18.01528)
    assert thermo.dG == pytest.approx(0.002267160)


@pytest.fixture
def hcn_geom():
    """Optimized at HF/STO-3G"""
    geom = geom_loader("lib:hcn_sto3g_freq_ref.xyz")
    return geom


@using("pyscf")
@using("thermoanalysis")
def test_get_thermoanalysis(hcn_geom):
    hcn_geom.set_calculator(PySCF(basis="sto3g", verbose=0, pal=2))
    thermo = get_thermoanalysis(hcn_geom)
    print_thermoanalysis(thermo)

    assert thermo.dG == pytest.approx(-0.006280, abs=1e-5)


@using("orca")
@using("thermoanalysis")
def test_hcn_thermo(hcn_geom):
    hcn_geom.set_calculator(ORCA(keywords="HF sto-3g"))
    thermo = get_thermoanalysis(hcn_geom)
    print_thermoanalysis(thermo)

    assert thermo.dG == pytest.approx(-0.006280, abs=1e-5)
