from pysisyphus.calculators.AnaPotBase import AnaPotBase


class SerpentinePot(AnaPotBase):
    def __init__(self):
        # NOTE: the whole potential is scaled further below; see the 'scale'
        # argument in the initializer call of the parent class
        V_str = (
            "0.4*exp(-0.05*x**2 - 0.05*y**2) "
            "- 0.5*exp(-0.025*x**2 - 0.025*y**2) "
            "- 0.1*exp(-2.5*(0.2*x - 1)**2 - 25.0*(0.2*y - 1)**2) "
            "- 0.1*exp(-0.625*(0.2*x + 1)**2 - 25.0*(0.2*y + 1)**2) "
            "+ 0.25*exp(-1.2*(0.25*x - 1)**2 - 4.0*(0.5*y - 1)**2) "
            "+ 0.25*exp(-1.2*(0.25*x + 1)**2 - 4.0*(0.5*y + 1)**2)"
        )
        lim = 8.0
        xlim = (-lim, lim)
        ylim = (-lim, lim)
        levels = 30
        minima = (
            (-5.0, -5.0, 0.0),
            (5.0, 5.0, 0.0),
        )
        saddles = ((0.0, 0.0, 0.0),)

        super().__init__(
            V_str=V_str,
            xlim=xlim,
            ylim=ylim,
            levels=levels,
            minima=minima,
            saddles=saddles,
            scale=0.25,
        )

    def __str__(self):
        return "SerpentinePot calculator"
