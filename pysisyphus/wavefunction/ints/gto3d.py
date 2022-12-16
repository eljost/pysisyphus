import numpy


def cart_gto3d_0(a, Xa, Ya, Za):
    """3D Cartesian s-Gaussian shell.

    Exponent a, centered at A, evaluated at (Xa, Ya, Za) + A.

    Generated code; DO NOT modify by hand!"""

    # 1 item(s)
    return numpy.array([numpy.exp(-a * (Xa**2 + Ya**2 + Za**2))])


def cart_gto3d_1(a, Xa, Ya, Za):
    """3D Cartesian p-Gaussian shell.

    Exponent a, centered at A, evaluated at (Xa, Ya, Za) + A.

    Generated code; DO NOT modify by hand!"""

    x0 = numpy.exp(-a * (Xa**2 + Ya**2 + Za**2))

    # 3 item(s)
    return numpy.array([Xa * x0, Ya * x0, Za * x0])


def cart_gto3d_2(a, Xa, Ya, Za):
    """3D Cartesian d-Gaussian shell.

    Exponent a, centered at A, evaluated at (Xa, Ya, Za) + A.

    Generated code; DO NOT modify by hand!"""

    x0 = Xa**2
    x1 = Ya**2
    x2 = Za**2
    x3 = numpy.exp(-a * (x0 + x1 + x2))
    x4 = Xa * x3

    # 6 item(s)
    return numpy.array([x0 * x3, Ya * x4, Za * x4, x1 * x3, Ya * Za * x3, x2 * x3])


def cart_gto3d_3(a, Xa, Ya, Za):
    """3D Cartesian f-Gaussian shell.

    Exponent a, centered at A, evaluated at (Xa, Ya, Za) + A.

    Generated code; DO NOT modify by hand!"""

    x0 = Xa**2
    x1 = Ya**2
    x2 = Za**2
    x3 = numpy.exp(-a * (x0 + x1 + x2))
    x4 = x0 * x3
    x5 = Xa * x3

    # 10 item(s)
    return numpy.array(
        [
            Xa**3 * x3,
            Ya * x4,
            Za * x4,
            x1 * x5,
            Ya * Za * x5,
            x2 * x5,
            Ya**3 * x3,
            Za * x1 * x3,
            Ya * x2 * x3,
            Za**3 * x3,
        ]
    )


def cart_gto3d_4(a, Xa, Ya, Za):
    """3D Cartesian g-Gaussian shell.

    Exponent a, centered at A, evaluated at (Xa, Ya, Za) + A.

    Generated code; DO NOT modify by hand!"""

    x0 = Xa**2
    x1 = Ya**2
    x2 = Za**2
    x3 = numpy.exp(-a * (x0 + x1 + x2))
    x4 = Xa**3 * x3
    x5 = x0 * x3
    x6 = Ya**3
    x7 = Xa * x3
    x8 = Za**3

    # 15 item(s)
    return numpy.array(
        [
            Xa**4 * x3,
            Ya * x4,
            Za * x4,
            x1 * x5,
            Ya * Za * x5,
            x2 * x5,
            x6 * x7,
            Za * x1 * x7,
            Ya * x2 * x7,
            x7 * x8,
            Ya**4 * x3,
            Za * x3 * x6,
            x1 * x2 * x3,
            Ya * x3 * x8,
            Za**4 * x3,
        ]
    )
