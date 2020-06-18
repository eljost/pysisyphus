import h5py
import pytest

from pysisyphus.testing import using


@using("thermoanalysis")
def test_thermoanalysis():
    hess_fn = "h2o_hessian.h5"
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
