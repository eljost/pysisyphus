import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.interpolate.Interpolator import Interpolator
from pysisyphus.interpolate.Redund import Redund


@pytest.mark.parametrize(
    "interpol_cls",
    [
        Interpolator,
        # LST,
        # IDPP,
        Redund,
    ],
)
def test_oxygen_extrapolate(interpol_cls):
    between = 3
    extrapolate = 3
    geoms = geom_loader("lib:oxygen_extrpol.trj")
    interpolator = interpol_cls(
        geoms, between=between, extrapolate=extrapolate, align=True
    )
    geoms = interpolator.interpolate_all()
    interpolator.all_geoms_to_trj("interpolated.trj")

    assert len(geoms) == (between + 2 + 2 * extrapolate)


@pytest.mark.parametrize(
    "interpol_cls, label",
    (
        (Interpolator, "cart"),
        (Redund, "redund"),
    ),
)
def test_diels_alder_extrapolate(interpol_cls, label):
    between = 5
    extrapolate = 4
    geoms = geom_loader("lib:diels_alder_extrapol.trj")
    interpolator = interpol_cls(
        geoms, between=between, extrapolate=extrapolate, extrapolate_damp=0.5,
        align=False
    )
    geoms = interpolator.interpolate_all()
    assert len(geoms) == (between + 2 + 2 * extrapolate)
    fn = f"diels_alder_interpolated_{label}.trj"
    interpolator.all_geoms_to_trj(fn)
    # plot_da(fn, label)
