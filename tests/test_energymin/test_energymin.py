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
        pytest.param(XTB, -6.27121126, marks=using("xtb")),
        pytest.param(
            lambda **kwargs: PySCF(basis="sto3g", **kwargs),
            -77.07357762,
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


# @pytest.mark.skip
@using("pyscf")
@pytest.mark.parametrize("energy_min_mix, ref_cycle", ((True, 7), (False, 3)))
def test_energy_min_cos(energy_min_mix, ref_cycle):
    fn1 = "lib:ethene_b3lyp_631gd.xyz"
    fn2 = "lib:ethene_rot_b3lyp_631gd.xyz"
    image1 = geom_loader(fn1)
    image2 = geom_loader(fn2)

    def calc_getter(mult):
        return PySCF(basis="631g*", xc="b3lyp", mult=mult, pal=8)

    between = 10
    images = interpolate(image1, image2, between=between, kind="redund")
    images = [image.copy(coord_type="cart") for image in images]

    def get_calculator():
        calc1 = calc_getter(mult=1)
        calc2 = calc_getter(mult=3)
        calc = EnergyMin(calc1, calc2)
        return calc

    for image in images:
        image.set_calculator(get_calculator())

    h5_group_name = "mix_opt" if energy_min_mix else "opt"
    cos = NEB(images, progress=True, energy_min_mix=energy_min_mix)
    opt = LBFGS(cos, dump=True, max_step=0.1, align=True, h5_group_name=h5_group_name)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
