import pytest

from pysisyphus.calculators import EnergyMin, XTB
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.cos.NEB import NEB
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.testing import using
from pysisyphus.interpolate import interpolate
from pysisyphus.optimizers.LBFGS import LBFGS


init_logging()


@pytest.mark.parametrize(
    "calc_cls, ref_energy",
    [
        pytest.param(XTB, -6.27036401645, marks=using("xtb")),
        pytest.param(
            lambda **kwargs: PySCF(basis="sto3g", **kwargs),
            -77.07214828,
            marks=using("pyscf"),
        ),
    ],
)
@using("xtb")
def test_energy_min_calc(calc_cls, ref_energy):
    geom = geom_loader("lib:ethene.xyz")
    calc1 = calc_cls(mult=1)
    calc2 = calc_cls(mult=3, calc_number=1)
    calc = EnergyMin(calc1, calc2)
    geom.set_calculator(calc)
    energy = geom.energy

    assert energy == pytest.approx(ref_energy)


@pytest.mark.skip
def test_energy_min_cos():
    image1 = geom_loader("lib:ethene.xyz")
    image2 = geom_loader("lib:ethene_rot.xyz")

    images = interpolate(image1, image2, between=9, kind="redund")
    images = [image.copy(coord_type="cart") for image in images]
    print(images)

    def get_calculator():
        calc1 = PySCF(basis="631g*", xc="b3lyp", mult=1, pal=2)
        calc2 = PySCF(basis="631g*", xc="b3lyp", mult=3, pal=2)
        calc = EnergyMin(calc1, calc2)
        return calc

    for image in images:
        image.set_calculator(get_calculator())

    cos = NEB(images)
    opt = LBFGS(cos, dump=True, max_cycles=10)
    opt.run()
