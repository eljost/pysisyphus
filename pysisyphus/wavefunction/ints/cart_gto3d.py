import numpy


_L_MAX = 4


def cart_gto3d_0(ax, da, Xa, Ya, Za):
    """3D Cartesian s-Gaussian shell.

    Exponent a, centered at A, evaluated at (Xa, Ya, Za) + A.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1,), dtype=float)

    # 1 item(s)
    result[0] = numpy.sum(
        0.71270547035499
        * da
        * numpy.sqrt(ax**1.5)
        * numpy.exp(-ax * (Xa**2 + Ya**2 + Za**2))
    )
    return result


def cart_gto3d_1(ax, da, Xa, Ya, Za):
    """3D Cartesian p-Gaussian shell.

    Exponent a, centered at A, evaluated at (Xa, Ya, Za) + A.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3,), dtype=float)

    x0 = (
        1.42541094070998
        * da
        * numpy.sqrt(ax**2.5)
        * numpy.exp(-ax * (Xa**2 + Ya**2 + Za**2))
    )

    # 3 item(s)
    result[0] = numpy.sum(Xa * x0)
    result[1] = numpy.sum(Ya * x0)
    result[2] = numpy.sum(Za * x0)
    return result


def cart_gto3d_2(ax, da, Xa, Ya, Za):
    """3D Cartesian d-Gaussian shell.

    Exponent a, centered at A, evaluated at (Xa, Ya, Za) + A.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6,), dtype=float)

    x0 = Xa**2
    x1 = Ya**2
    x2 = Za**2
    x3 = 0.423777208123758 * da * numpy.sqrt(ax**3.5) * numpy.exp(-ax * (x0 + x1 + x2))
    x4 = 3.88393417365859 * x3
    x5 = 6.72717132202972 * Xa * x3

    # 6 item(s)
    result[0] = numpy.sum(x0 * x4)
    result[1] = numpy.sum(Ya * x5)
    result[2] = numpy.sum(Za * x5)
    result[3] = numpy.sum(x1 * x4)
    result[4] = numpy.sum(6.72717132202972 * Ya * Za * x3)
    result[5] = numpy.sum(x2 * x4)
    return result


def cart_gto3d_3(ax, da, Xa, Ya, Za):
    """3D Cartesian f-Gaussian shell.

    Exponent a, centered at A, evaluated at (Xa, Ya, Za) + A.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10,), dtype=float)

    x0 = Xa**2
    x1 = Ya**2
    x2 = Za**2
    x3 = 0.423777208123758 * da * numpy.sqrt(ax**4.5) * numpy.exp(-ax * (x0 + x1 + x2))
    x4 = 3.47389633297403 * x3
    x5 = 7.76786834731718 * x3
    x6 = x0 * x5
    x7 = Xa * x5

    # 10 item(s)
    result[0] = numpy.sum(Xa**3 * x4)
    result[1] = numpy.sum(Ya * x6)
    result[2] = numpy.sum(Za * x6)
    result[3] = numpy.sum(x1 * x7)
    result[4] = numpy.sum(13.4543426440594 * Xa * Ya * Za * x3)
    result[5] = numpy.sum(x2 * x7)
    result[6] = numpy.sum(Ya**3 * x4)
    result[7] = numpy.sum(Za * x1 * x5)
    result[8] = numpy.sum(Ya * x2 * x5)
    result[9] = numpy.sum(Za**3 * x4)
    return result


def cart_gto3d_4(ax, da, Xa, Ya, Za):
    """3D Cartesian g-Gaussian shell.

    Exponent a, centered at A, evaluated at (Xa, Ya, Za) + A.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15,), dtype=float)

    x0 = Xa**2
    x1 = Ya**2
    x2 = Za**2
    x3 = 0.423777208123758 * da * numpy.sqrt(ax**5.5) * numpy.exp(-ax * (x0 + x1 + x2))
    x4 = 2.62601879356243 * x3
    x5 = 6.94779266594806 * x3
    x6 = Xa**3 * x5
    x7 = x0 * x3
    x8 = 8.96956176270629 * x7
    x9 = 15.5357366946344 * Za
    x10 = Ya**3
    x11 = Xa * x5
    x12 = x1 * x3
    x13 = Za**3

    # 15 item(s)
    result[0] = numpy.sum(Xa**4 * x4)
    result[1] = numpy.sum(Ya * x6)
    result[2] = numpy.sum(Za * x6)
    result[3] = numpy.sum(x1 * x8)
    result[4] = numpy.sum(Ya * x7 * x9)
    result[5] = numpy.sum(x2 * x8)
    result[6] = numpy.sum(x10 * x11)
    result[7] = numpy.sum(Xa * x12 * x9)
    result[8] = numpy.sum(15.5357366946344 * Xa * Ya * x2 * x3)
    result[9] = numpy.sum(x11 * x13)
    result[10] = numpy.sum(Ya**4 * x4)
    result[11] = numpy.sum(Za * x10 * x5)
    result[12] = numpy.sum(8.96956176270629 * x12 * x2)
    result[13] = numpy.sum(Ya * x13 * x5)
    result[14] = numpy.sum(Za**4 * x4)
    return result


cart_gto3d = {
    (0,): cart_gto3d_0,
    (1,): cart_gto3d_1,
    (2,): cart_gto3d_2,
    (3,): cart_gto3d_3,
    (4,): cart_gto3d_4,
}
