import pytest

from pysisyphus.constants import AU2KJPERMOL, C
from pysisyphus.drivers import (
    eyring_rate,
    harmonic_tst_rate,
    bell_corr,
    eckart_corr,
    eckart_corr_brown,
    wigner_corr,
)
from pysisyphus.drivers.rates import get_rates_for_geoms, render_rx_rates
from pysisyphus.io import geom_from_hessian
from pysisyphus.testing import using


T = 298.15
BARRIER = 104.5 / AU2KJPERMOL  # about 25 kcal mol⁻¹


def test_eyring():
    rate = eyring_rate(BARRIER, temperature=T)
    assert rate == pytest.approx(3.059526e-06)


@pytest.mark.parametrize(
    "rs_part_func, ts_part_func, ref_rate", ((1.0, 1.0, 3.059526e-6),)
)
def test_harmonic_tst(rs_part_func, ts_part_func, ref_rate):
    rate = harmonic_tst_rate(BARRIER, T, rs_part_func, ts_part_func)
    assert rate == pytest.approx(ref_rate)


def test_wigner_corr():
    imag_frequency = -1.6e13
    kappa = wigner_corr(temperature=298.15, imag_frequency=imag_frequency)
    assert kappa == pytest.approx(1.276379)


@pytest.mark.parametrize(
    "nu, T, kappa_ref",
    (
        (572, 272, 1.52),
        (742, 420, 1.33),
        (1000, 225, -57.44),
    ),
)
def test_bell_corr(nu, T, kappa_ref):
    imag_frequency = nu * C * 100
    kappa = bell_corr(temperature=T, imag_frequency=imag_frequency)
    assert kappa == pytest.approx(kappa_ref, abs=1e-2)


@pytest.mark.parametrize(
    "fw_barrier_height, bw_barrier_height, nu, ref_rate",
    (
        (7.62, 15, 1000, 2.43885387),
        (15, 7.62, 1000, 2.43885387),
    ),
)
def test_eckart_corr(fw_barrier_height, bw_barrier_height, nu, ref_rate):
    temperature = 298.15
    imag_frequency = nu * C * 100  # Convert from wavenumbers (1/cm) to (1/s)

    fw_barrier_height /= AU2KJPERMOL  # from kJ mol⁻¹ to au
    bw_barrier_height /= AU2KJPERMOL  # from kJ mol⁻¹ to au

    # Implementation using scipy.integrate.quad
    kappa = eckart_corr(
        fw_barrier_height, bw_barrier_height, temperature, imag_frequency
    )
    # Implementation by Brown using 6-point Gaussian quadrature
    assert kappa == pytest.approx(ref_rate)

    kappa_brown = eckart_corr_brown(
        fw_barrier_height, bw_barrier_height, temperature, imag_frequency
    )
    # Gives slightly different results
    assert kappa_brown == pytest.approx(ref_rate, abs=5e-3)


@using("thermoanalysis")
@pytest.mark.parametrize(
    "with_product",
    (True, False),
)
def test_rx_rates(with_product, this_dir):
    reactant = geom_from_hessian(this_dir / "peroxyl_r1_h5s/reactant.h5")
    product = geom_from_hessian(this_dir / "peroxyl_r1_h5s/product.h5")
    ts_geom = geom_from_hessian(this_dir / "peroxyl_r1_h5s/ts.h5")
    reactant_geoms = (reactant,)
    if with_product:
        product_geoms = (product,)
    else:
        product_geoms = tuple()

    rx_rates = get_rates_for_geoms(T, reactant_geoms, ts_geom, product_geoms)
    for rxr in rx_rates:
        print(render_rx_rates(rxr))
