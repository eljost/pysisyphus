import numpy


def cart_gto3d_0(a, A, R):
    """3D Cartesian s-Gaussian shell.

    Exponent a, centered at A, evaluated at R.

    Generated code; DO NOT modify by hand!"""

    # 1 item(s)
    return numpy.array(
        [numpy.exp(-a * ((A[0] - R[0]) ** 2 + (A[1] - R[1]) ** 2 + (A[2] - R[2]) ** 2))]
    )


def cart_gto3d_1(a, A, R):
    """3D Cartesian p-Gaussian shell.

    Exponent a, centered at A, evaluated at R.

    Generated code; DO NOT modify by hand!"""

    x0 = A[0] - R[0]
    x1 = A[1] - R[1]
    x2 = A[2] - R[2]
    x3 = numpy.exp(-a * (x0 ** 2 + x1 ** 2 + x2 ** 2))

    # 3 item(s)
    return numpy.array([-x0 * x3, -x1 * x3, -x2 * x3])


def cart_gto3d_2(a, A, R):
    """3D Cartesian d-Gaussian shell.

    Exponent a, centered at A, evaluated at R.

    Generated code; DO NOT modify by hand!"""

    x0 = A[0] - R[0]
    x1 = -x0
    x2 = A[1] - R[1]
    x3 = A[2] - R[2]
    x4 = numpy.exp(-a * (x0 ** 2 + x2 ** 2 + x3 ** 2))
    x5 = -x2
    x6 = x1 * x4
    x7 = -x3

    # 6 item(s)
    return numpy.array(
        [x1 ** 2 * x4, x5 * x6, x6 * x7, x4 * x5 ** 2, x4 * x5 * x7, x4 * x7 ** 2]
    )


def cart_gto3d_3(a, A, R):
    """3D Cartesian f-Gaussian shell.

    Exponent a, centered at A, evaluated at R.

    Generated code; DO NOT modify by hand!"""

    x0 = A[0] - R[0]
    x1 = -x0
    x2 = A[1] - R[1]
    x3 = A[2] - R[2]
    x4 = numpy.exp(-a * (x0 ** 2 + x2 ** 2 + x3 ** 2))
    x5 = -x2
    x6 = x1 ** 2 * x4
    x7 = -x3
    x8 = x5 ** 2
    x9 = x1 * x4
    x10 = x7 ** 2

    # 10 item(s)
    return numpy.array(
        [
            x1 ** 3 * x4,
            x5 * x6,
            x6 * x7,
            x8 * x9,
            x5 * x7 * x9,
            x10 * x9,
            x4 * x5 ** 3,
            x4 * x7 * x8,
            x10 * x4 * x5,
            x4 * x7 ** 3,
        ]
    )


def cart_gto3d_4(a, A, R):
    """3D Cartesian g-Gaussian shell.

    Exponent a, centered at A, evaluated at R.

    Generated code; DO NOT modify by hand!"""

    x0 = A[0] - R[0]
    x1 = -x0
    x2 = A[1] - R[1]
    x3 = A[2] - R[2]
    x4 = numpy.exp(-a * (x0 ** 2 + x2 ** 2 + x3 ** 2))
    x5 = -x2
    x6 = x1 ** 3 * x4
    x7 = -x3
    x8 = x5 ** 2
    x9 = x1 ** 2 * x4
    x10 = x7 ** 2
    x11 = x5 ** 3
    x12 = x1 * x4
    x13 = x7 ** 3

    # 15 item(s)
    return numpy.array(
        [
            x1 ** 4 * x4,
            x5 * x6,
            x6 * x7,
            x8 * x9,
            x5 * x7 * x9,
            x10 * x9,
            x11 * x12,
            x12 * x7 * x8,
            x10 * x12 * x5,
            x12 * x13,
            x4 * x5 ** 4,
            x11 * x4 * x7,
            x10 * x4 * x8,
            x13 * x4 * x5,
            x4 * x7 ** 4,
        ]
    )
