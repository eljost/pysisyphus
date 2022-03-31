import pytest

from pysisyphus.constants import AU2KJPERMOL, C
from pysisyphus.drivers import (
    eyring_rate,
    harmonic_tst_rate,
    wigner_corr,
    eckart_corr,
    eckart_corr_brown,
)


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
