# [1] http://dx.doi.org/10.1142/S0217732319502080
#     Exact solutions of a quartic potential
#     Dong, Sun, Aoki, Chen, Dong, 2019
# [2] https://doi.org/10.1007/s10910-022-01328-9
#     Exact solutions of an asymmetric double well potential
#     Sun, Dong, Bezerra, Dong, 20222

from collections.abc import Callable
from typing import Optional

import numpy as np
import numpy.typing as npt
import sympy as sym

from pysisyphus.calculators.Calculator import Calculator


def get_sym_a_b(en_barr=None, sep=None):
    """Solve for parameters a and b of double well potential.

    f(x, a, b) = a*x**2 + b*x**4

    Parameters
    ----------
    en_barr
        Potential at x = 0.0, barrier height.
    sep
        Separation of minima. Minima will be located at ±(sep/2.0).
    Returns
    -------
    a
        Double well potential parameter a. Scales quadratic term.
    b
        Double well potential parameter b. Scales quartic term.
    """
    a = sym.symbols("a", negative=True, real=True)
    b = sym.symbols("b", positive=True, real=True)
    if en_barr is None:
        en_barr = sym.symbols("en_barr", positive=True, real=True)
    if sep is None:
        sep = sym.symbols("sep", positive=True, real=True)

    expr1 = -(a**2) / (4 * b) + en_barr
    aabs = sym.Abs(a)
    expr2 = sym.sqrt(aabs) / sym.sqrt(2 * b) - sep / 2.0
    results = sym.solve([expr1, expr2], a, b)
    assert len(results) >= 1
    a0, b0 = results[0]
    return a0, b0


def get_sym_double_well(en_barr: float, sep: float, use_sym=False) -> Callable:
    """1D Symmetric double well potential, centered at 0.0.

    Parameters
    ----------
    en_barr
        Potential at x = 0.0, barrier height.
    sep
        Separation of minima. Minima will be located at
        ±(sep/2.0).

    Returns
    -------
    inner
        Function to evaluate the potential. Takes one float.
    """
    if use_sym:
        a, b = get_sym_a_b(en_barr, sep)
    else:
        a = -8.0 * en_barr / sep**2
        b = 16.0 * en_barr / sep**4

    def inner(x):
        x2 = x**2
        return a * x2 + b * x2**2 + en_barr

    return inner


def get_quadratic(en_min: float, x1: float, y1: float) -> Callable:
    """1D quadratic potential f(x) = a * x**2 + en_min

    Parameters
    ----------
    en_min
        Function value at x = 0.0; vertical shift.
    x1
        Potential argument x1. Used to calculate parameter a.
    y1
        Potential value at x1. Used to calculate parameter a.

    Returns
    -------
    innner
         Function to evaluate the potential. Takes one float.
    """

    a = (y1 - en_min) / x1**2

    def inner(x):
        return a * x**2 + en_min

    return inner


def get_pot_func(en_barr, en_above_barr, en_above_gs_min, sep):
    """Get function to evaluate 2 PEC (double well GS and quadratic ES).

    Parameters
    ----------
    en_barr
        Barrier height in the GS at x = 0.0
    en_above_barr
        Energy difference between potential value in the ES at x = 0.0 and
        the barrier height.
    en_above_gs_min
        Energy in the ES above the GS minimum.

    Returns
    -------
    innner
         Function to evaluate the potential. Takes one float,
         returns array of two floats.
    """
    gs = get_sym_double_well(en_barr, sep)
    gs_min = -sep / 2.0
    es = get_quadratic(en_barr + en_above_barr, gs_min, en_above_gs_min)

    def inner(x):
        return np.array((gs(x), es(x)), dtype=float)

    return inner


def logistic(
    x: npt.ArrayLike,
    L: float = 1.0,
    k: float = 1.0,
    x0: float = 0.0,
    factor: float = 1.0,
) -> np.ndarray:
    """Logistic function.

    Parameters
    ----------
    x
        Actual function argument.
    L
        Supremum value of the function.
    k
        Logistic growth rate/steepness.
    x0
        Function midpoint.
    factor
        Scaling factor for x. With factor = -1.0 the function
        can be mirrored at the y-axis.

    Returns
    -------
    y
        Value of the logistic function with the given parameters
        at the given point(s).
    """
    return L / (1 + np.exp(-k * (factor * (x - x0))))


def get_epos_property(ref_coords, flip=False):
    x0 = 0.0
    factor = -1.0 if flip else 1.0
    k = 50  # Steepness

    def inner(geom, *args):
        x = (geom.cart_coords - ref_coords)[0]
        return logistic(x, k=k, factor=factor, x0=x0)

    return inner


class RobinDayClass2(Calculator):
    """1D model of Robin-Day Class II mixed-valence system with 2 states.

    4500 +------------------------------------------------+
         |                                                |
         |                                                |
    4000 |                                                |
         |                                                |
         |                                                |
    3500 |**                                            **|
         |  ***                                      ***  |
         |     **                                  **     |
    3000 |       ***                            ***       |
         |       |  ****                    ****          |
         |       |      ********    ********              |
    2500 |       |              ****                      |
         |       |                |                       |
         |       |                |                       |
    2000 |       |                |                       |
         |       |                |                       |
         | en_above_gs_min   en_above_barr                |
    1500 |       |                |                       |
         |*      |<-----sep/2---->|                       |
         |*      |                |                      *|
    1000 | *     |             ******                    *|
         | *     |          ***   |  ***                * |
         |  *    |        **      |     **              * |
     500 |  *    |      **     en_barr    **           *  |
         |   *   |     *          |         *         *   |
         |    *  |   ***          |          ***     *    |
       0 +------------------------------------------------+
       -0.3    -0.2    -0.1       0      0.1     0.2     0.3


    Parameters
    ----------
    en_barr
        Barrier height in the GS at x = 0.0
    en_above_barr
        Energy difference between potential value in the ES at x = 0.0 and
        the barrier height.
    en_above_gs_min
        Energy in the ES above the GS minimum.
    ref_coords
        Reference coordinates that will be substracted from the provided
        coordinates when get_energy() is called. This allows this calculator
        to be used with an actual Geometry.
    flip
        Whether the logistic function mimicing the electron-position along
        the curve is mirrored at the line passing through the midpoint.

    Returns
    -------
    calculator
         Pysisyphus calculator. Atom argument is ignored and only the FIRST item
         of the provided coordinates is utilized as argument to the analytical
         potentials.
    """

    def __init__(
        self,
        en_barr: float,
        en_above_barr: float,
        en_above_gs_min: float,
        sep: float,
        ref_coords: Optional[np.ndarray] = None,
        flip_epos=False,
        **kwargs,
    ):
        assert en_barr > 0.0
        assert en_above_barr > 0.0
        assert en_above_gs_min > (en_barr + en_above_barr)
        assert sep > 0.0

        super().__init__(**kwargs)
        self.pot_func = get_pot_func(en_barr, en_above_barr, en_above_gs_min, sep)
        try:
            ref_coords = ref_coords.copy()
        except AttributeError:
            self.log("Didn't set any reference coordinates.")
        self.ref_coords = ref_coords
        self.flip_epos = flip_epos

        self.epos_property = get_epos_property(
            self.ref_coords.copy(), flip=self.flip_epos
        )

    def get_energy(self, atoms, coords):
        if self.ref_coords is not None:
            coords = coords - self.ref_coords
        x = coords[0]
        all_energies = self.pot_func(x)
        results = {
            "energy": all_energies[0],
            "all_energies": all_energies,
        }
        return results


def do_plots(show=False):
    en_barr = 1000
    en_above_barr = 1500
    en_above_gs_min = 3000
    sep = 0.4
    calc = RobinDayClass2(en_barr, en_above_barr, en_above_gs_min, sep)
    print(calc)

    num = 50
    xs = np.linspace(-0.3, 0.3, num=num)
    all_coords = np.zeros((num, 3))
    all_coords[:, 0] = xs

    all_energies = np.empty((xs.size, 2))
    for i, coords in enumerate(all_coords):
        res = calc.get_energy(None, coords)
        all_energies[i] = res["all_energies"]
    gs, es = all_energies.T

    try:
        import termplotlib as tpl

        tfig = tpl.figure()
        width = 60
        height = 30
        ylim = 0, 4500
        tfig.plot(xs, es, label="ES", width=width, height=height, ylim=ylim)
        tfig.plot(xs, gs, label="GS", width=width, height=height, ylim=ylim)
        if show:
            tfig.show()
    except ModuleNotFoundError:
        print("'termplotlib' is missing (python -m pip install termplotlib)!")

    import matplotlib.pyplot as plt

    _, ax = plt.subplots()
    ax.plot(xs, gs, label="GS")
    ax.plot(xs, es, label="ES")
    ax.set_ylim(0, 8000)
    ax.legend()
    if show:
        plt.show()


if __name__ == "__main__":
    do_plots(show=False)
