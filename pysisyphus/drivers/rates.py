# [1] https://doi.org/10.1002/ange.201914943
#     Heavy-Atom Tunneling in Organic Reactions
#     Castro, Karney, 2020
# [2] https://doi.org/10.1002/wcms.1165
#     Theory and simulation of atom tunneling in chemical reactions
#     Kästner, 2014
# [3] https://doi.org/10.1002/qua.25686
#     Eyringpy
#     Dzib, Merino et al, 2018
# [4] https://pubs.acs.org/doi/pdf/10.1021/j100809a040
#     TUNNELLING CORRECTIONS FOR UNSYMMETRICAL ECKART POTENTIAL ENERGY BARRIERS
#     Johnston, Heicklen, 1961
# [5] https://dx.doi.org/10.6028%2Fjres.086.014
#     A Method of Calculating Tunneling Corrections For Eckart Potential Barriers
#     Brown, 1981
# [6] https://doi.org/10.1021/jp062775l
#     Theoretical Determination of the Rate Constant for OH Hydrogen Abstraction
#     from Toluene
#     Uc, Alvarez-Idaboy, Galano, García-Cruz, Vivier-Bunge, 2006

from dataclasses import dataclass
from math import sin
from typing import Literal, Optional, Sequence
import warnings

import jinja2
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import scipy.integrate as integrate

try:
    from thermoanalysis.QCData import QCData
    from thermoanalysis.thermo import thermochemistry

    can_rates = True
except (ImportError, ModuleNotFoundError):
    can_rates = False

from pysisyphus.constants import (
    AU2KJPERMOL,
    AU2SEC,
    BOHR2ANG,
    C,
    KB,
    KBAU,
    PLANCK,
    PLANCKAU,
    R,
)
from pysisyphus.io import geom_from_hessian


@dataclass
class ReactionRates:
    from_: str
    barrier: float  # in E_h
    barrier_si: float  # in kJ mol⁻¹
    temperature: float  # in K
    imag_wavenumber: float  # in cm⁻¹
    imag_frequency: float  # in s⁻¹
    rate_eyring: float  # in s⁻¹
    kappa_eyring: float
    unit: str
    molecularity: int
    rate_wigner: Optional[float] = None  # in s⁻¹
    kappa_wigner: Optional[float] = None
    rate_bell: Optional[float] = None  # in s⁻¹
    kappa_bell: Optional[float] = None
    rate_eckart: Optional[float] = None  # in s⁻¹
    kappa_eckart: Optional[float] = None
    rate_eckart_brown: Optional[float] = None  # in s⁻¹
    kappa_eckart_brown: Optional[float] = None


RX_RATES_TPL = """
             From: {{ "{: >18s}".format(rxr.from_) }} to TS
          Barrier: {{ "{: >18.6f}".format(rxr.barrier) }} E_h ( {{ rxr.barrier_si|f2 }} kJ mol⁻¹)
      Temperature: {{ "{: >18.2f}".format(rxr.temperature) }} K

Transition vector:
       Wavenumber: {{ "{: >18.2f}".format(rxr.imag_wavenumber) }} cm⁻¹
        Frequency: {{ "{: >18.6e}".format(rxr.imag_frequency|abs) }} s⁻¹

 Type        kappa       rate / {{ "%13s"|format(unit) }}       comment
 -------  ------------ ----------------------- ------------------
 {{ rate_line("Eyring", rxr.kappa_eyring, rxr.rate_eyring) }}
{%- if rxr.rate_wigner %}
 {{ rate_line("Wigner", rxr.kappa_wigner, rxr.rate_wigner, "1d tunnel corr.") }}
{%- endif %}
{%- if rxr.rate_bell %}
 {{ rate_line("Bell", rxr.kappa_bell, rxr.rate_bell, "1d tunnel corr.") }}
{%- endif %}
{%- if rxr.rate_eckart %}
 {{ rate_line("Eckart", rxr.kappa_eckart, rxr.rate_eckart, "1d tunnel corr.") }}
{%- endif %}
{%- if rxr.rate_eckart_brown %}
 {{ rate_line("Eckart'", rxr.kappa_eckart_brown, rxr.rate_eckart_brown, "1d tunnel corr.") }}
{%- endif %}
 -------  ------------ ----------------------- ------------------
"""


def render_rx_rates(rx_rates: ReactionRates) -> str:
    def rate_line(name, kappa, rate, comment=""):
        return f"{name: <12s} {kappa: >8.4f} {rate: >23.8e} {comment: >18s}"

    env = jinja2.Environment()
    env.filters["f2"] = lambda _: f"{_:.2f}"
    env.filters["f4"] = lambda _: f"{_:.4f}"
    env.filters["f6"] = lambda _: f"{_:.6f}"
    env.filters["e8"] = lambda _: f"{_: >16.8e}"

    tpl = env.from_string(RX_RATES_TPL)

    rendered = tpl.render(rxr=rx_rates, rate_line=rate_line, unit=rx_rates.unit)
    return rendered


# def eyring_rate(
# barrier_height: float,
# temperature: npt.ArrayLike,
# trans_coeff: float = 1.0,
# ) -> npt.NDArray[np.float64]:
# """Rate constant in 1/s from the Eyring equation.

# See https://pubs.acs.org/doi/10.1021/acs.organomet.8b00456
# "Reaction Rates and Barriers" on p. 3234 and eq. (8).

# Parameters
# ----------
# barrier_height
# Barrier height (energy, enthalpy, gibbs energy, ...) in Hartree.
# temperature
# Temperature in Kelvin.
# trans_coeff
# Unitless transmission coefficient, e.g., obtained from Wigner or Eckart
# correction.

# Returns
# -------
# rate
# Eyring reaction rate in 1/s.
# """
# temperature = np.array(temperature, dtype=float)
# prefactor = trans_coeff * KB * temperature / PLANCK
# rate = prefactor * np.exp(-barrier_height / (KBAU * temperature))
# return rate


def eyring_rate(
    barrier_height: npt.ArrayLike,
    temperature: npt.ArrayLike,
    molecularity: int,
    pressure: float,
    # trans_coeff: float = 1.0,
    # degen: int = 1,
) -> npt.NDArray[np.float64]:
    """Eyring equation rate constant. Unit depends on the molecularity of the
    reaction. It will be s⁻¹ for unimolecular reactions and cm³ mol⁻¹ s⁻¹ for
    bimolecular reactions and so on ((cm³ mol⁻¹ s⁻¹)**(molecularity-1).

    TODO: add citation


    Parameters
    ----------
    barrier_height
        Barrier height (energy, enthalpy, gibbs energy, ...) in Hartree.
    temperature
        Temperature in Kelvin.
    pressure
        Pressure in Pascal.
    trans_coeff
        Unitless transmission coefficient, e.g., obtained from Wigner or Eckart
        correction.

    Returns
    -------
    rate
        Eyring reaction rate in 1/s.
    """
    temperature = np.array(temperature, dtype=float)
    # unit(prefactor1) = 1 / s
    prefactor = KB * temperature / PLANCK
    # unit(prefactor2) = Depends on the moleculary: (cm³ / mol)**molecularity
    prefactor2 = (R * temperature * 1e6 / pressure) ** (molecularity - 1)
    rate = prefactor * prefactor2 * np.exp(-barrier_height / (KBAU * temperature))
    return rate


def wigner_corr(
    temperature: float,
    imag_frequency: float,
) -> float:
    """Tunneling correction according to Wigner.

    See https://doi.org/10.1002/qua.25686 eq. (12) and
    https://doi.org/10.1007/978-3-642-59033-7_9 for the original publication.

    Parameters
    ----------
    temperature
        Temperature in Kelvin.
    imag_frequency
        Imaginary frequency in 1/s.

    Returns
    -------
    kappa
        Unitless tunneling correction according to Wigner.
    """
    kappa = 1 + 1 / 24 * (PLANCK * abs(imag_frequency) / (KB * temperature)) ** 2
    return kappa


def bell_corr(
    temperature: float,
    imag_frequency: float,
) -> float:
    """Tunneling correction according to Bell.

    See https://onlinelibrary.wiley.com/doi/10.1002/anie.201708489 eq. (1) and
    eq. (2).

    Parameters
    ----------
    temperature
        Temperature in Kelvin.
    imag_frequency
        Imaginary frequency in 1/s.

    Returns
    -------
    kappa
        Unitless tunneling correction according to Bell. Negative kappas are
        meaningless.
    """
    u = PLANCK * abs(imag_frequency) / (KB * temperature)
    kappa = (u / 2) / sin(u / 2)
    return kappa


class EckartNotApplicableError(Exception):
    pass


def eckart_corr(
    fw_barrier_height: float,
    bw_barrier_height: float,
    temperature: float,
    imag_frequency: float,
) -> float:
    """Tunneling correction according to Eckart.

    See [3], [4] and [5]. The correction should be independent of the order
    fw_barrier_height/bw_barrier_height.

    Parameters
    ----------
    fw_barrier_height
        Barrier height in forward direction in Hartree.
    bw_barrier_height
        Barrier height in backward direction in Hartree.
    temperature
        Temperature in Kelvin.
    imag_frequency
        Frequency in 1/s of the imaginary mode at the TS.

    Returns
    -------
    kappa
        Unitless tunneling correction according to Eckart.
    """
    warnings.warn("'eckart_corr' seems to fail for small absolute imag_wavenumbers!")
    kBT = KBAU * temperature  # Hartree
    two_pi = 2 * np.pi

    imag_frequency = abs(imag_frequency) * AU2SEC  # Convert from 1/s to 1/au
    hnu = PLANCKAU * imag_frequency
    u = hnu / kBT  # unitless
    v1 = fw_barrier_height / kBT
    v2 = bw_barrier_height / kBT
    alpha1, alpha2 = [
        two_pi * barrier / hnu for barrier in (fw_barrier_height, bw_barrier_height)
    ]
    quot = 1 / (1 / np.sqrt(alpha1) + 1 / np.sqrt(alpha2))
    d = 1 / two_pi * np.sqrt(np.abs(4 * alpha1 * alpha2 - np.pi**2))
    cosh_d = np.cosh(two_pi * np.abs(d))

    def eps(E):
        return (E - fw_barrier_height) / kBT

    def ai(E, vi):
        epsilon = eps(E)
        return np.sqrt(2 * (epsilon + vi) / (np.pi * u)) * quot

    def _a1(E):
        return ai(E, vi=v1)

    def _a2(E):
        return ai(E, vi=v2)

    def eckart_probability(E):
        a1 = _a1(E)
        a2 = _a2(E)
        plus = np.cosh(two_pi * (a1 + a2))
        minus = np.cosh(two_pi * (a1 - a2))
        return (plus - minus) / (plus + cosh_d)

    def eckart_func(E):
        epsilon = eps(E)
        P = eckart_probability(E)
        return P * np.exp(-epsilon)

    # Integration limits
    E_low = max(0, fw_barrier_height - bw_barrier_height)
    E_max = 1  # 1 Hartree max
    kappa, _ = integrate.quad(eckart_func, E_low, E_max, limit=250)
    # Make kappa unitless again by dividing by kBT, as we integrated over the energy.
    kappa /= kBT
    return kappa


def tunl(alph1: float, alph2: float, U: float, strict: bool = False):
    """Eckart correction factor for rate constants according to Brown.

    Python adaption of the TUNL subroutine in 4. Appendix of [5].

    Parameters
    ----------
    alph1
        Unitless barrier height descriptor. 2π V1 / (h nu*); see (2) in [5].
    alph2
        Unitless barrier heigth descriptor. 2π V2 / (h nu*); see (2) in [5].
    U
        Unitless curvature descriptor. h nu* / kT; see (2) in [5].
    strict
        If enabled, arguments are bound checked. Will raise AssertionError if
        they are out of bonds. TUNL was found to yield accurate results when
        the arguments are within bounds.

    Returns
    -------
    G
        Unitless tunneling correction according to Eckart.
    """
    if strict:
        try:
            alphs = np.array((alph1, alph2))
            assert (0.5 <= alphs).all() and (alphs <= 20).all() and 0 < U <= 16
            if (alphs >= 8).all():
                assert U <= 12
            if (alphs >= 16).all():
                assert U <= 10
        except AssertionError:
            raise EckartNotApplicableError()

    # Quadrature arguments
    x = np.array((-0.9324695, -0.6612094, -0.2386192, 0.2386192, 0.6612094, 0.9324695))
    w = np.array((0.1713245, 0.3607616, 0.4679139, 0.4679139, 0.3607616, 0.1713245))

    pi = np.pi
    upi2 = U / (2 * pi)

    C = 0.125 * pi * U * (1 / np.sqrt(alph1) + 1.0 / np.sqrt(alph2)) ** 2
    v1 = upi2 * alph1
    v2 = upi2 * alph2
    D = 4 * alph1 * alph2 - pi**2
    DF = np.cosh(np.sqrt(D)) if D >= 0 else np.cos(np.sqrt(-D))

    ez = -v1 if (v2 >= v1) else -v2
    eb = min(
        C * (np.log(2.0 * (1.0 + DF) / 0.014) / (2 * pi)) ** 2 - (v1 + v2) / 2, 3.2
    )
    em = (eb - ez) / 2
    ep = (eb + ez) / 2

    # Original loop
    E = em * x + ep
    A1 = pi * np.sqrt((E + v1) / C)
    A2 = pi * np.sqrt((E + v2) / C)
    fm = np.cosh(A1 - A2)
    fp = np.cosh(A1 + A2)
    G = (w * np.exp(-E) * (fp - fm) / (fp + DF)).sum()
    # Outside the loop
    G = em * G + np.exp(-eb)
    return G


def eckart_corr_brown(
    fw_barrier_height: float,
    bw_barrier_height: float,
    temperature: float,
    imag_frequency: float,
    strict: bool = True,
) -> float:
    """Tunneling correction according to Eckart.

    Wrapper for the TUNL subroutine as given in the appendix of [5].

    Parameters
    ----------
    fw_barrier_height
        Barrier height in forward direction in Hartree.
    bw_barrier_height
        Barrier height in backward direction in Hartree.
    temperature
        Temperature in Kelvin.
    imag_frequency
        Frequency in 1/s of the imaginary mode at the TS.
    strict
        Enable bound-checking in tunl-function.

    Returns
    -------
    kappa
        Unitless tunneling correction according to Eckart.
    """
    kBT = KBAU * temperature
    hnu = PLANCKAU * imag_frequency * AU2SEC
    u = hnu / kBT

    def alpha(barr):
        return 2 * np.pi * barr / hnu

    alpha1 = alpha(fw_barrier_height)
    alpha2 = alpha(bw_barrier_height)
    kappa = tunl(alpha1, alpha2, u, strict=strict)
    return kappa


def get_rate_unit(molecularity: int):
    if molecularity == 1:
        unit = "s⁻¹"
    elif molecularity == 2:
        unit = "cm³ mol⁻¹ s⁻¹"
    else:
        unit = f"(cm³ mol⁻¹)**{molecularity-1} s⁻¹"
    return unit


def get_rates(temperature, reactant_thermos, ts_thermo, *, degen, product_thermos=None):
    """Driver for reaction rate calculation."""
    fw_molecularity = len(reactant_thermos)
    fw_unit = get_rate_unit(fw_molecularity)
    assert fw_molecularity >= 1
    pressure = ts_thermo.p

    imag_wavenumber = ts_thermo.org_wavenumbers[0]
    assert (
        imag_wavenumber < 0.0
    ), f"There is no imaginary mode! Lowest wavenumber is {imag_wavenumber: >10.4f} cm⁻¹!"
    imag_frequency = (
        imag_wavenumber * C * 100
    )  # Convert from wavenumbers (1/cm) to (1/s)
    kappa_wigner = wigner_corr(temperature, imag_frequency)
    kappa_bell = bell_corr(temperature, imag_frequency)

    Gs_reactant = [thermo.G for thermo in reactant_thermos]
    G_reactant = sum(Gs_reactant)
    G_TS = ts_thermo.G
    fw_barrier_height = G_TS - G_reactant
    fw_rate_eyring = degen * eyring_rate(
        fw_barrier_height, temperature, fw_molecularity, pressure=pressure
    )

    if product_thermos:
        Gs_product = [thermo.G for thermo in product_thermos]
        G_product = sum(Gs_product)
        bw_barrier_height = G_TS - G_product
        # bw_rate_eyring = degen * eyring_rate(bw_barrier_height, temperature, bw_molecularity, pressure)
        kappa_eckart = eckart_corr(
            fw_barrier_height, bw_barrier_height, temperature, imag_frequency
        )
        try:
            kappa_eckart_brown = eckart_corr_brown(
                fw_barrier_height,
                bw_barrier_height,
                temperature,
                abs(imag_frequency),
                strict=False,
            )
        except EckartNotApplicableError:
            kappa_eckart_brown = None
    else:
        kappa_eckart = None
        kappa_eckart_brown = None

    def make_rx_rates(from_, barrier, rate_eyring, kappa_eyring=1.0):
        kwargs = {}
        if kappa_wigner and not np.isnan(kappa_wigner):
            kwargs.update(
                {
                    "kappa_wigner": kappa_wigner,
                    "rate_wigner": kappa_wigner * rate_eyring,
                }
            )
        if kappa_bell and kappa_bell > 0.0:
            kwargs.update(
                {
                    "kappa_bell": kappa_bell,
                    "rate_bell": kappa_bell * rate_eyring,
                }
            )
        if kappa_eckart and not np.isnan(kappa_eckart):
            kwargs.update(
                {
                    "kappa_eckart": kappa_eckart,
                    "rate_eckart": kappa_eckart * rate_eyring,
                }
            )
        if kappa_eckart_brown and not np.isnan(kappa_eckart_brown):
            kwargs.update(
                {
                    "kappa_eckart_brown": kappa_eckart_brown,
                    "rate_eckart_brown": kappa_eckart_brown * rate_eyring,
                }
            )
        rx_rates = ReactionRates(
            from_=from_,
            barrier=barrier,
            barrier_si=barrier * AU2KJPERMOL,
            temperature=temperature,
            imag_wavenumber=imag_wavenumber,
            imag_frequency=imag_frequency,
            rate_eyring=rate_eyring,
            kappa_eyring=kappa_eyring,
            molecularity=fw_molecularity,
            unit=fw_unit,
            **kwargs,
        )
        return rx_rates

    reactant_rx_rates = make_rx_rates("Reactant(s)", fw_barrier_height, fw_rate_eyring)
    return reactant_rx_rates
    # rx_rates = [
    # reactant_rx_rates,
    # ]
    # if product_thermos:
    # product_rx_rates = make_rx_rates(
    # "Product(s)", bw_barrier_height, bw_rate_eyring
    # )
    # rx_rates.append(product_rx_rates)
    # return rx_rates


def get_rates_for_geoms(temperature, reactant_geoms, ts_geom, degen, product_geoms):
    def get_thermos(geoms):
        return [geom.get_thermoanalysis(T=temperature) for geom in geoms]

    reactant_thermos = get_thermos(reactant_geoms)
    ts_thermo = get_thermos((ts_geom,))[0]
    product_thermos = get_thermos(product_geoms)
    return get_rates(
        temperature,
        reactant_thermos,
        ts_thermo,
        degen=degen,
        product_thermos=product_thermos,
    )


@dataclass
class RateInput:
    atoms: tuple[str, ...]
    coords3d: np.ndarray
    masses: np.ndarray
    mult: int
    energy: float
    vibfreqs: np.ndarray

    # def to_qcdata(self) -> QCData:
    def to_qcdata(self):
        coords3d_ang = self.coords3d * BOHR2ANG
        kwargs = {
            "masses": self.masses,
            "wavenumbers": self.vibfreqs,
            "coords3d": coords3d_ang,
            "scf_energy": self.energy,
            "mult": self.mult,
        }
        qcdata = QCData(kwargs)
        return qcdata

    @staticmethod
    def from_h5_hessian(
        h5_fn,
        alt_energy: Optional[float] = None,
    ):
        geom, attrs = geom_from_hessian(h5_fn, with_attrs=True)
        # This could also be directly loaded/read from the Hessian
        vibfreqs, *_ = geom.get_normal_modes()
        energy = alt_energy if alt_energy is not None else attrs["energy"]
        mult = attrs["mult"]
        return RateInput(
            atoms=geom.atoms,
            coords3d=geom.coords3d,
            masses=geom.masses,
            mult=mult,
            energy=energy,
            vibfreqs=vibfreqs,
        )


def get_rates_for_rate_inputs(
    temperatures: Sequence[float],
    reactant_inps: list[RateInput],
    ts_inp: RateInput,
    degen: int,
    product_inps: Optional[list[RateInput]] = None,
    thermo_kwargs: Optional[dict] = None,
):
    """

    Parameters
    ----------
    temperatures
        1d array of temperatures for which the rates will be calculated.
    degen
        Positive integer. Degeneracy factor of the path.
    vib_scale
        Scale factor for the vibrational wavenumbers.
    zpe_scale
        Scale factor for the zero point energy.
    """

    if not can_rates:
        raise ModuleNotFoundError(
            "get_rates_for_rate_inputs() requires the thermoanalysis packages "
            "(https://github.com/eljost/thermoanalysis)"
        )

    assert min(temperatures) > 0
    assert degen >= 1
    if thermo_kwargs is None:
        thermo_kwargs = dict()
    if product_inps is None:
        product_inps = list()

    reactant_qcds = [reactant.to_qcdata() for reactant in reactant_inps]
    ts_qcd = ts_inp.to_qcdata()
    product_qcds = [product.to_qcdata() for product in product_inps]

    all_rx_rates = list()
    for temperature in temperatures:
        reactant_thermos = [
            thermochemistry(qcd, temperature=temperature, **thermo_kwargs.copy())
            for qcd in reactant_qcds
        ]
        ts_thermo = thermochemistry(
            ts_qcd, temperature=temperature, **thermo_kwargs.copy()
        )
        product_thermos = [
            thermochemistry(qcd, temperature=temperature, **thermo_kwargs.copy())
            for qcd in product_qcds
        ]
        rx_rates = get_rates(
            temperature,
            reactant_thermos,
            ts_thermo,
            degen=degen,
            product_thermos=product_thermos,
        )
        all_rx_rates.append(rx_rates)
    return all_rx_rates


def plot_rates(
    all_rx_rates,
    kind: Literal["eyring", "wigner", "bell", "eckart", "eckart_brown"],
    show=False,
):
    rxr0 = all_rx_rates[0]
    unit = rxr0.unit

    # Sort from low to high temperatures
    all_rx_rates = sorted(all_rx_rates, key=lambda rxr: rxr.temperature)
    temperatures = np.array([rxr.temperature for rxr in all_rx_rates])
    attr = f"rate_{kind}"
    rates = np.array([getattr(rxr, attr) for rxr in all_rx_rates])

    _temps = 1000 / temperatures

    fig, ax = plt.subplots()

    ax.plot(_temps, rates, "o-", label=kind)
    ax.legend()
    ax.set_yscale("log")
    ax.set_xlim(_temps[-1], _temps[0])
    ax.set_xlabel("1000 T⁻¹ / K")
    # Label depends on molecularity
    ax.set_ylabel(f"k / {unit}")
    ax.axhline(rates[0], c="k", ls=":")
    ax.axhline(rates[-1], c="k", ls=":")
    ax.annotate(f"{rates[0]: >8.2e} {unit}", (_temps[-1], rates[0] * 1.1))
    ax.annotate(f"{rates[-1]: >8.2e} {unit}", (_temps[-1], rates[-1] * 1.1))
    fig.tight_layout()
    if show:
        plt.show()
    return fig, ax
