from pysisyphus.interpolate import IDPP, LST, Interpolator, Redund


INTERPOLATE = {
    "idpp": IDPP.IDPP,
    "lst": LST.LST,
    "redund": Redund.Redund,
    "linear": Interpolator.Interpolator,
}


def interpolate(initial_geom, final_geom, between, kind="linear",
                only_between=False, align=False, interpol_kwargs=None):
    if interpol_kwargs is None:
        interpol_kwargs = {}
    """ only_between: return only the interpolated images."""
    interpolate_class = INTERPOLATE[kind]
    interpolator = interpolate_class((initial_geom, final_geom), between,
                                     align=align, **interpol_kwargs)
    if only_between:
        return interpolator.interpolate(initial_geom, final_geom)
    else:
        return interpolator.interpolate_all()


def interpolate_all(geoms, between, kind="linear", align=False):
    interpolate_class = INTERPOLATE[kind]
    interpolator = interpolate_class(geoms, between, align=align)

    return interpolator.interpolate_all()
