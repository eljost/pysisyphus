import numpy as np
import pytest

from pysisyphus.drivers import rates


@pytest.mark.parametrize(
    "ts_h5_fn, ts_alt_energy, degen, name, rate_wigner0",
    (
        ("02_ts3_final_hessian.h5", -347.057742, 3, "TS3", 1006290420.93),
        ("03_ts2_final_hessian.h5", -347.064226, 2, "TS2", 32609879744.87),
    ),
)
def test_ethyl_propionate_rates(
    ts_h5_fn, ts_alt_energy, degen, name, rate_wigner0, this_dir
):
    """From https://doi.org/10.1016/j.fuel.2024.131492

    # Gaussian 16 M06-2X/6-311++G(d, p) w/ CBS-QB3 single points
    """
    base = this_dir / "ethyl_propionate"

    Hri = rates.RateInput.from_h5_hessian(
        base / "00_h_final_hessian.h5", alt_energy=-0.499818
    )
    EPri = rates.RateInput.from_h5_hessian(
        base / "01_ep_final_hessian.h5", alt_energy=-346.578342
    )
    TSri = rates.RateInput.from_h5_hessian(base / ts_h5_fn, alt_energy=ts_alt_energy)

    temperatures = np.linspace(500, 2500, num=21)
    thermo_kwargs = {
        "kind": "qrrho",
        "scale_factor": 0.983,
        "zpe_scale_factor": 0.97,
    }
    all_rx_rates = rates.get_rates_for_rate_inputs(
        temperatures, [Hri, EPri], TSri, degen=degen, thermo_kwargs=thermo_kwargs
    )
    fig, ax = rates.plot_rates(all_rx_rates, "wigner", show=False)
    # fig.savefig(this_dir / f"rates_{name}.png")
    rxr0 = all_rx_rates[0]
    assert rxr0.rate_wigner == pytest.approx(rate_wigner0)
