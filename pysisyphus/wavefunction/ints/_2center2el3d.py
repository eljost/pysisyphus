import numpy

from pysisyphus.wavefunction.ints.boys import boys


_L_MAX = 5


def _2center2el3d_00(ax, da, A, bx, db, B):
    """Cartesian (s|s) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 1), dtype=float)

    x0 = ax + bx

    # 1 item(s)
    result[0, 0] = numpy.sum(
        17.7715317526335
        * da
        * db
        * x0 ** (-0.5)
        * numpy.sqrt(ax**1.5)
        * numpy.sqrt(bx**1.5)
        * boys(
            0,
            ax * bx * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2) / x0,
        )
        / (ax * bx)
    )
    return result


def _2center2el3d_01(ax, da, A, bx, db, B):
    """Cartesian (s|p) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 3), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = (
        35.5430635052669
        * da
        * db
        * x0 ** (-0.5)
        * numpy.sqrt(ax**1.5)
        * numpy.sqrt(bx**2.5)
        * boys(
            1,
            ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2),
        )
        / (ax * bx)
    )

    # 3 item(s)
    result[0, 0] = numpy.sum(-x2 * (-x1 * (ax * A[0] + bx * B[0]) + B[0]))
    result[0, 1] = numpy.sum(-x2 * (-x1 * (ax * A[1] + bx * B[1]) + B[1]))
    result[0, 2] = numpy.sum(-x2 * (-x1 * (ax * A[2] + bx * B[2]) + B[2]))
    return result


def _2center2el3d_02(ax, da, A, bx, db, B):
    """Cartesian (s|d) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 6), dtype=float)

    x0 = bx ** (-1.0)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = ax * bx * x2 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x4 = 17.4934183276249
    x5 = 2.0 * x0 * x4
    x6 = ax ** (-1.0)
    x7 = x1 ** (-0.5)
    x8 = x0**2 * x4 * x6 * x7 * boys(0, x3) - 0.5 * x0 * x1 ** (-1.5) * x5 * boys(1, x3)
    x9 = x2 * (ax * A[0] + bx * B[0]) - B[0]
    x10 = boys(2, x3)
    x11 = x6 * x7
    x12 = x10 * x11 * x5
    x13 = da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**3.5)
    x14 = 1.17305816955079 * x13
    x15 = x2 * (ax * A[1] + bx * B[1]) - B[1]
    x16 = 3.14159265358979
    x17 = 22.6274169979695 * x0 * x10 * x11 * x13 * x16 * x9
    x18 = x2 * (ax * A[2] + bx * B[2]) - B[2]

    # 6 item(s)
    result[0, 0] = numpy.sum(x14 * (x12 * x9**2 + x8))
    result[0, 1] = numpy.sum(x15 * x17)
    result[0, 2] = numpy.sum(x17 * x18)
    result[0, 3] = numpy.sum(x14 * (x12 * x15**2 + x8))
    result[0, 4] = numpy.sum(22.6274169979695 * x0 * x10 * x11 * x13 * x15 * x16 * x18)
    result[0, 5] = numpy.sum(x14 * (x12 * x18**2 + x8))
    return result


def _2center2el3d_03(ax, da, A, bx, db, B):
    """Cartesian (s|f) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 10), dtype=float)

    x0 = bx ** (-1.0)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = x2 * (ax * A[0] + bx * B[0]) - B[0]
    x4 = ax * bx * x2 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x5 = 17.4934183276249
    x6 = 2.0 * x0 * x5
    x7 = x1 ** (-1.5) * x6 * boys(2, x4)
    x8 = ax ** (-1.0)
    x9 = x1 ** (-0.5)
    x10 = boys(1, x4)
    x11 = 0.5 * x0 * (2.0 * x0 * x10 * x5 * x8 * x9 - x7)
    x12 = boys(3, x4)
    x13 = x8 * x9
    x14 = x12 * x13 * x6
    x15 = x14 * x3**2
    x16 = da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**4.5)
    x17 = 0.179587122125167 * x16
    x18 = 5.84237394672177 * x17
    x19 = x2 * (ax * A[1] + bx * B[1]) - B[1]
    x20 = x0 * x19 * (2.0 * x0 * x10 * x5 * x8 * x9 - x7)
    x21 = 13.0639452948436 * x17
    x22 = x2 * (ax * A[2] + bx * B[2]) - B[2]
    x23 = x0 * x22 * (2.0 * x0 * x10 * x5 * x8 * x9 - x7)
    x24 = 0.5 * x23
    x25 = x14 * x19**2
    x26 = x11 + x25
    x27 = x21 * x3
    x28 = x11 + x14 * x22**2

    # 10 item(s)
    result[0, 0] = numpy.sum(
        x18 * x3 * (x0 * (2.0 * x0 * x10 * x5 * x8 * x9 - x7) + x11 + x15)
    )
    result[0, 1] = numpy.sum(0.5 * x21 * (2.0 * x15 * x19 + x20))
    result[0, 2] = numpy.sum(x21 * (x15 * x22 + x24))
    result[0, 3] = numpy.sum(x26 * x27)
    result[0, 4] = numpy.sum(142.172254021068 * x0 * x12 * x13 * x16 * x19 * x22 * x3)
    result[0, 5] = numpy.sum(x27 * x28)
    result[0, 6] = numpy.sum(x18 * (x19 * x26 + x20))
    result[0, 7] = numpy.sum(x21 * (x22 * x25 + x24))
    result[0, 8] = numpy.sum(x19 * x21 * x28)
    result[0, 9] = numpy.sum(x18 * (x22 * x28 + x23))
    return result


def _2center2el3d_04(ax, da, A, bx, db, B):
    """Cartesian (s|g) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 15), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = x1 * (ax * A[0] + bx * B[0]) - B[0]
    x3 = bx ** (-1.0)
    x4 = ax * x1
    x5 = bx * x4 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x6 = boys(3, x5)
    x7 = 17.4934183276249
    x8 = 2.0 * x3 * x7
    x9 = x0 ** (-1.5) * x8
    x10 = x6 * x9
    x11 = x10 * x2
    x12 = ax ** (-1.0)
    x13 = x0 ** (-0.5)
    x14 = boys(2, x5)
    x15 = 0.5 * x3
    x16 = x15 * (-x10 + 2.0 * x12 * x13 * x14 * x3 * x7)
    x17 = boys(4, x5)
    x18 = x2**2
    x19 = x12 * x13 * x8
    x20 = x18 * x19
    x21 = x17 * x20
    x22 = boys(1, x5)
    x23 = x15 * (2.0 * x12 * x13 * x3 * x7 * boys(0, x5) - x22 * x9)
    x24 = x15 * (2.0 * x12 * x13 * x22 * x3 * x7 - x14 * x9)
    x25 = x14 * x19
    x26 = 1.5 * x3
    x27 = 0.179587122125167 * da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**5.5)
    x28 = 4.41641957979107 * x27
    x29 = x1 * (ax * A[1] + bx * B[1]) - B[1]
    x30 = x10 * x29
    x31 = x3 * (2.0 * x12 * x13 * x14 * x29 * x3 * x7 - x30)
    x32 = x21 * x29
    x33 = 11.6847478934435 * x27
    x34 = x1 * (ax * A[2] + bx * B[2]) - B[2]
    x35 = x3 * x34 * (-x10 + 2.0 * x12 * x13 * x14 * x3 * x7)
    x36 = 0.5 * x35
    x37 = x29**2
    x38 = x19 * x37
    x39 = x17 * x38
    x40 = x16 + x39
    x41 = x23 + x25 * x37 - x4 * (x24 + x38 * x6)
    x42 = 15.084944665313 * x27
    x43 = x3 * x34 * (2.0 * x12 * x13 * x14 * x29 * x3 * x7 - x30)
    x44 = 26.1278905896872 * x27
    x45 = x34**2
    x46 = x19 * x45
    x47 = x16 + x17 * x46
    x48 = x23 + x25 * x45 - x4 * (x24 + x46 * x6)
    x49 = x15 * x48
    x50 = x29 * x40 + x31
    x51 = x2 * x33
    x52 = x34 * x39 + x36
    x53 = x2 * x44
    x54 = x34 * x47 + x35

    # 15 item(s)
    result[0, 0] = numpy.sum(
        x28
        * (
            x2 * (x2 * (x16 + x21) - x3 * (x11 - 2.0 * x12 * x13 * x14 * x2 * x3 * x7))
            + x26 * (x18 * x25 + x23 - x4 * (x20 * x6 + x24))
        )
    )
    result[0, 1] = numpy.sum(
        0.5
        * x33
        * (
            x2 * (x31 + 2.0 * x32)
            - 2.0 * x29 * x3 * (x11 - 2.0 * x12 * x13 * x14 * x2 * x3 * x7)
        )
    )
    result[0, 2] = numpy.sum(
        x33
        * (
            x2 * (x21 * x34 + x36)
            - x3 * x34 * (x11 - 2.0 * x12 * x13 * x14 * x2 * x3 * x7)
        )
    )
    result[0, 3] = numpy.sum(x42 * (x15 * x41 + x18 * x40))
    result[0, 4] = numpy.sum(0.5 * x44 * (2.0 * x32 * x34 + x43))
    result[0, 5] = numpy.sum(x42 * (x18 * x47 + x49))
    result[0, 6] = numpy.sum(x50 * x51)
    result[0, 7] = numpy.sum(x52 * x53)
    result[0, 8] = numpy.sum(x29 * x47 * x53)
    result[0, 9] = numpy.sum(x51 * x54)
    result[0, 10] = numpy.sum(x28 * (x26 * x41 + x29 * x50))
    result[0, 11] = numpy.sum(x33 * (x29 * x52 + x43))
    result[0, 12] = numpy.sum(x42 * (x37 * x47 + x49))
    result[0, 13] = numpy.sum(x29 * x33 * x54)
    result[0, 14] = numpy.sum(x28 * (x26 * x48 + x34 * x54))
    return result


def _2center2el3d_05(ax, da, A, bx, db, B):
    """Cartesian (s|h) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 21), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = x1 * (ax * A[0] + bx * B[0]) - B[0]
    x3 = bx ** (-1.0)
    x4 = ax * x1
    x5 = bx * x4 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x6 = boys(4, x5)
    x7 = 17.4934183276249
    x8 = 2.0 * x3
    x9 = x7 * x8
    x10 = x0 ** (-1.5) * x9
    x11 = x10 * x6
    x12 = x11 * x2
    x13 = ax ** (-1.0)
    x14 = x0 ** (-0.5)
    x15 = boys(3, x5)
    x16 = 0.5 * x3
    x17 = x16 * (-x11 + 2.0 * x13 * x14 * x15 * x3 * x7)
    x18 = boys(5, x5)
    x19 = x2**2
    x20 = x13 * x14 * x9
    x21 = x19 * x20
    x22 = x18 * x21
    x23 = x10 * x15
    x24 = boys(2, x5)
    x25 = x16 * (2.0 * x13 * x14 * x24 * x3 * x7 - x23)
    x26 = x21 * x6
    x27 = x25 + x26
    x28 = x10 * x24
    x29 = boys(1, x5)
    x30 = x16 * (2.0 * x13 * x14 * x29 * x3 * x7 - x28)
    x31 = x15 * x20
    x32 = x19 * x31
    x33 = x30 + x32
    x34 = 1.5 * x3
    x35 = 0.179587122125167 * da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**6.5)
    x36 = 2.94427971986071 * x35
    x37 = x1 * (ax * A[1] + bx * B[1]) - B[1]
    x38 = x12 * x37
    x39 = x11 * x37
    x40 = x3 * (2.0 * x13 * x14 * x15 * x3 * x37 * x7 - x39)
    x41 = x22 * x37
    x42 = x3 * x37 * (2.0 * x13 * x14 * x29 * x3 * x7 - x28)
    x43 = x3 * x37 * (2.0 * x13 * x14 * x24 * x3 * x7 - x23)
    x44 = 8.83283915958214 * x35
    x45 = x1 * (ax * A[2] + bx * B[2]) - B[2]
    x46 = x3 * x45 * (-x11 + 2.0 * x13 * x14 * x15 * x3 * x7)
    x47 = 0.5 * x46
    x48 = x3 * x45 * (2.0 * x13 * x14 * x29 * x3 * x7 - x28)
    x49 = 0.5 * x48
    x50 = x3 * x45 * (2.0 * x13 * x14 * x24 * x3 * x7 - x23)
    x51 = 0.5 * x50
    x52 = x37**2
    x53 = x30 + x31 * x52
    x54 = x20 * x52
    x55 = x54 * x6
    x56 = x25 + x55
    x57 = x4 * x56
    x58 = x18 * x54
    x59 = x17 + x58
    x60 = x53 - x57
    x61 = 13.4923846833851 * x35
    x62 = x3 * x45 * (2.0 * x13 * x14 * x15 * x3 * x37 * x7 - x39)
    x63 = 23.3694957868871 * x35
    x64 = x45**2
    x65 = x30 + x31 * x64
    x66 = x20 * x64
    x67 = x25 + x6 * x66
    x68 = x4 * x67
    x69 = x17 + x18 * x66
    x70 = x19 * x69
    x71 = x65 - x68
    x72 = x16 * x71
    x73 = x37 * x59 + x40
    x74 = x37 * x53 - x4 * (x37 * x56 + x43) + x42
    x75 = 13.4923846833851 * x35
    x76 = x45 * x58 + x47
    x77 = x31 * x45 * x52 - x4 * (x45 * x55 + x51) + x49
    x78 = 30.169889330626 * x35
    x79 = x3 * x37 * (x65 - x68)
    x80 = x45 * x69 + x46
    x81 = -x4 * (x45 * x67 + x50) + x45 * x65 + x48
    x82 = x16 * x81
    x83 = x34 * x60 + x37 * x73
    x84 = x2 * x44
    x85 = x37 * x76 + x62
    x86 = x2 * x63
    x87 = x52 * x69 + x72
    x88 = x34 * x71 + x45 * x80

    # 21 item(s)
    result[0, 0] = numpy.sum(
        x2
        * x36
        * (
            x2 * (x2 * (x17 + x22) - x3 * (x12 - 2.0 * x13 * x14 * x15 * x2 * x3 * x7))
            - x34 * (x27 * x4 - x33)
            + x8
            * (
                x3 * (2.0 * x13 * x14 * x29 * x3 * x7 - x28)
                + x33
                - x4 * (x27 + x3 * (2.0 * x13 * x14 * x24 * x3 * x7 - x23))
            )
        )
    )
    result[0, 1] = numpy.sum(
        0.5
        * x44
        * (
            x2
            * (
                x2 * (x40 + 2.0 * x41)
                + 2.0 * x3 * (2.0 * x13 * x14 * x15 * x2 * x3 * x37 * x7 - x38)
            )
            + x34 * (2.0 * x32 * x37 - x4 * (2.0 * x26 * x37 + x43) + x42)
        )
    )
    result[0, 2] = numpy.sum(
        x44
        * (
            x2
            * (
                x2 * (x22 * x45 + x47)
                - x3 * x45 * (x12 - 2.0 * x13 * x14 * x15 * x2 * x3 * x7)
            )
            + x34 * (x32 * x45 - x4 * (x26 * x45 + x51) + x49)
        )
    )
    result[0, 3] = numpy.sum(x2 * x61 * (x16 * x60 + x19 * x59 + x3 * (x53 - x57)))
    result[0, 4] = numpy.sum(
        0.5
        * x63
        * (
            x2 * (2.0 * x41 * x45 + x62)
            + 2.0 * x3 * x45 * (2.0 * x13 * x14 * x15 * x2 * x3 * x37 * x7 - x38)
        )
    )
    result[0, 5] = numpy.sum(x2 * x61 * (x3 * (x65 - x68) + x70 + x72))
    result[0, 6] = numpy.sum(x75 * (x16 * x74 + x19 * x73))
    result[0, 7] = numpy.sum(x78 * (x16 * x77 + x19 * x76))
    result[0, 8] = numpy.sum(0.5 * x78 * (2.0 * x37 * x70 + x79))
    result[0, 9] = numpy.sum(x75 * (x19 * x80 + x82))
    result[0, 10] = numpy.sum(x83 * x84)
    result[0, 11] = numpy.sum(x85 * x86)
    result[0, 12] = numpy.sum(x2 * x78 * x87)
    result[0, 13] = numpy.sum(x37 * x80 * x86)
    result[0, 14] = numpy.sum(x84 * x88)
    result[0, 15] = numpy.sum(x36 * (x37 * x83 + x74 * x8))
    result[0, 16] = numpy.sum(x44 * (x34 * x77 + x37 * x85))
    result[0, 17] = numpy.sum(x61 * (x37 * x87 + x79))
    result[0, 18] = numpy.sum(x75 * (x52 * x80 + x82))
    result[0, 19] = numpy.sum(x37 * x44 * x88)
    result[0, 20] = numpy.sum(x36 * (x45 * x88 + x8 * x81))
    return result


def _2center2el3d_10(ax, da, A, bx, db, B):
    """Cartesian (p|s) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 1), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = (
        35.5430635052669
        * da
        * db
        * x0 ** (-0.5)
        * numpy.sqrt(ax**2.5)
        * numpy.sqrt(bx**1.5)
        * boys(
            1,
            ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2),
        )
        / (ax * bx)
    )

    # 3 item(s)
    result[0, 0] = numpy.sum(-x2 * (-x1 * (ax * A[0] + bx * B[0]) + A[0]))
    result[1, 0] = numpy.sum(-x2 * (-x1 * (ax * A[1] + bx * B[1]) + A[1]))
    result[2, 0] = numpy.sum(-x2 * (-x1 * (ax * A[2] + bx * B[2]) + A[2]))
    return result


def _2center2el3d_11(ax, da, A, bx, db, B):
    """Cartesian (p|p) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 3), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x3 = x0 ** (-0.5) / (ax * bx)
    x4 = 34.9868366552497 * x3
    x5 = 0.5 * x4 * boys(1, x2) / (ax + bx)
    x6 = -x1 * (ax * A[0] + bx * B[0])
    x7 = -x6 - A[0]
    x8 = -x6 - B[0]
    x9 = boys(2, x2)
    x10 = x4 * x9
    x11 = da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**2.5)
    x12 = 2.03179634989571 * x11
    x13 = -x1 * (ax * A[1] + bx * B[1])
    x14 = -x13 - B[1]
    x15 = 71.0861270105339 * x11 * x3 * x9
    x16 = x15 * x7
    x17 = -x1 * (ax * A[2] + bx * B[2])
    x18 = -x17 - B[2]
    x19 = -x13 - A[1]
    x20 = x15 * x19
    x21 = -x17 - A[2]
    x22 = x15 * x21

    # 9 item(s)
    result[0, 0] = numpy.sum(x12 * (x10 * x7 * x8 + x5))
    result[0, 1] = numpy.sum(x14 * x16)
    result[0, 2] = numpy.sum(x16 * x18)
    result[1, 0] = numpy.sum(x20 * x8)
    result[1, 1] = numpy.sum(x12 * (x10 * x14 * x19 + x5))
    result[1, 2] = numpy.sum(x18 * x20)
    result[2, 0] = numpy.sum(x22 * x8)
    result[2, 1] = numpy.sum(x14 * x22)
    result[2, 2] = numpy.sum(x12 * (x10 * x18 * x21 + x5))
    return result


def _2center2el3d_12(ax, da, A, bx, db, B):
    """Cartesian (p|d) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 6), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = bx ** (-1.0)
    x5 = ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x6 = boys(2, x5)
    x7 = 17.4934183276249
    x8 = 2.0 * x4 * x7
    x9 = x6 * x8
    x10 = ax ** (-1.0)
    x11 = x0 ** (-0.5)
    x12 = -0.5 * x0 ** (-1.5) * x4 * x9 + x10 * x11 * x4**2 * x7 * boys(1, x5)
    x13 = -x2 - B[0]
    x14 = boys(3, x5)
    x15 = x10 * x11
    x16 = x14 * x15 * x8
    x17 = x12 + x13**2 * x16
    x18 = 0.5 * x15 / (ax + bx)
    x19 = x13 * x18
    x20 = 4.0 * x4 * x6 * x7
    x21 = da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**3.5)
    x22 = 0.179587122125167 * x21
    x23 = 13.0639452948436 * x22
    x24 = -x1 * (ax * A[1] + bx * B[1])
    x25 = -x24 - B[1]
    x26 = x18 * x9
    x27 = x25 * x26
    x28 = x16 * x25
    x29 = x13 * x3
    x30 = 22.6274169979695 * x22
    x31 = -x1 * (ax * A[2] + bx * B[2])
    x32 = -x31 - B[2]
    x33 = x26 * x32
    x34 = x16 * x32
    x35 = x12 + x16 * x25**2
    x36 = x23 * x3
    x37 = 3.14159265358979
    x38 = 45.2548339959391 * x14 * x15 * x21 * x32 * x37 * x4
    x39 = x12 + x16 * x32**2
    x40 = -x24 - A[1]
    x41 = x23 * x40
    x42 = x19 * x9
    x43 = x28 * x40
    x44 = x18 * x20
    x45 = -x31 - A[2]
    x46 = x23 * x45
    x47 = x13 * x45

    # 18 item(s)
    result[0, 0] = numpy.sum(x23 * (x17 * x3 + x19 * x20))
    result[0, 1] = numpy.sum(x30 * (x27 + x28 * x29))
    result[0, 2] = numpy.sum(x30 * (x29 * x34 + x33))
    result[0, 3] = numpy.sum(x35 * x36)
    result[0, 4] = numpy.sum(x25 * x3 * x38)
    result[0, 5] = numpy.sum(x36 * x39)
    result[1, 0] = numpy.sum(x17 * x41)
    result[1, 1] = numpy.sum(x30 * (x13 * x43 + x42))
    result[1, 2] = numpy.sum(x13 * x38 * x40)
    result[1, 3] = numpy.sum(x23 * (x25 * x44 + x35 * x40))
    result[1, 4] = numpy.sum(x30 * (x32 * x43 + x33))
    result[1, 5] = numpy.sum(x39 * x41)
    result[2, 0] = numpy.sum(x17 * x46)
    result[2, 1] = numpy.sum(45.2548339959391 * x14 * x15 * x21 * x25 * x37 * x4 * x47)
    result[2, 2] = numpy.sum(x30 * (x34 * x47 + x42))
    result[2, 3] = numpy.sum(x35 * x46)
    result[2, 4] = numpy.sum(x30 * (x27 + x28 * x32 * x45))
    result[2, 5] = numpy.sum(x23 * (x32 * x44 + x39 * x45))
    return result


def _2center2el3d_13(ax, da, A, bx, db, B):
    """Cartesian (p|f) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 10), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = bx ** (-1.0)
    x5 = -x2 - B[0]
    x6 = ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x7 = boys(3, x6)
    x8 = 17.4934183276249
    x9 = 2.0 * x4 * x8
    x10 = x0 ** (-1.5) * x9
    x11 = x10 * x7
    x12 = ax ** (-1.0)
    x13 = x0 ** (-0.5)
    x14 = boys(2, x6)
    x15 = 0.5 * x4
    x16 = x15 * (-x11 + 2.0 * x12 * x13 * x14 * x4 * x8)
    x17 = boys(4, x6)
    x18 = x12 * x13
    x19 = x18 * x9
    x20 = x19 * x5**2
    x21 = x17 * x20
    x22 = x5 * (x16 + x21 - x4 * (x11 - 2.0 * x12 * x13 * x14 * x4 * x8))
    x23 = 0.5 / (ax + bx)
    x24 = x15 * (-x10 * x14 + 2.0 * x12 * x13 * x4 * x8 * boys(1, x6))
    x25 = x23 * (x20 * x7 + x24)
    x26 = 0.179587122125167 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**4.5)
    x27 = 11.6847478934435 * x26
    x28 = -x1 * (ax * A[1] + bx * B[1])
    x29 = -x28 - B[1]
    x30 = x29 * x4 * (-x11 + 2.0 * x12 * x13 * x14 * x4 * x8)
    x31 = x21 * x29 + 0.5 * x30
    x32 = x23 * x29 * x7
    x33 = x32 * x5
    x34 = 4.0 * x18 * x4 * x8
    x35 = x33 * x34
    x36 = 26.1278905896872 * x26
    x37 = -x1 * (ax * A[2] + bx * B[2])
    x38 = -x37 - B[2]
    x39 = x38 * x4 * (-x11 + 2.0 * x12 * x13 * x14 * x4 * x8)
    x40 = 0.5 * x39
    x41 = x21 * x38 + x40
    x42 = x34 * x38
    x43 = x23 * x5 * x7
    x44 = x42 * x43
    x45 = x19 * x29**2
    x46 = x23 * (x24 + x45 * x7)
    x47 = x17 * x45
    x48 = x16 + x47
    x49 = x3 * x5
    x50 = x19 * x38
    x51 = x17 * x29 * x50
    x52 = 45.2548339959391 * x26
    x53 = x19 * x38**2
    x54 = x23 * (x24 + x53 * x7)
    x55 = x16 + x17 * x53
    x56 = x29 * x48 + x30
    x57 = x27 * x3
    x58 = x38 * x47 + x40
    x59 = x3 * x36
    x60 = x29 * x55
    x61 = x38 * x55 + x39
    x62 = -x28 - A[1]
    x63 = x27 * x62
    x64 = x5 * x62
    x65 = x32 * x42
    x66 = -x37 - A[2]
    x67 = x27 * x66
    x68 = x5 * x66

    # 30 item(s)
    result[0, 0] = numpy.sum(x27 * (x22 * x3 + 3.0 * x25))
    result[0, 1] = numpy.sum(x36 * (x3 * x31 + x35))
    result[0, 2] = numpy.sum(x36 * (x3 * x41 + x44))
    result[0, 3] = numpy.sum(x36 * (x46 + x48 * x49))
    result[0, 4] = numpy.sum(x52 * (x32 * x50 + x49 * x51))
    result[0, 5] = numpy.sum(x36 * (x49 * x55 + x54))
    result[0, 6] = numpy.sum(x56 * x57)
    result[0, 7] = numpy.sum(x58 * x59)
    result[0, 8] = numpy.sum(x59 * x60)
    result[0, 9] = numpy.sum(x57 * x61)
    result[1, 0] = numpy.sum(x22 * x63)
    result[1, 1] = numpy.sum(x36 * (x25 + x31 * x62))
    result[1, 2] = numpy.sum(x36 * x41 * x62)
    result[1, 3] = numpy.sum(x36 * (x35 + x48 * x64))
    result[1, 4] = numpy.sum(x52 * (x43 * x50 + x51 * x64))
    result[1, 5] = numpy.sum(x36 * x55 * x64)
    result[1, 6] = numpy.sum(x27 * (3.0 * x46 + x56 * x62))
    result[1, 7] = numpy.sum(x36 * (x58 * x62 + x65))
    result[1, 8] = numpy.sum(x36 * (x54 + x60 * x62))
    result[1, 9] = numpy.sum(x61 * x63)
    result[2, 0] = numpy.sum(x22 * x67)
    result[2, 1] = numpy.sum(x31 * x36 * x66)
    result[2, 2] = numpy.sum(x36 * (x25 + x41 * x66))
    result[2, 3] = numpy.sum(x36 * x48 * x68)
    result[2, 4] = numpy.sum(x52 * (x19 * x33 + x51 * x68))
    result[2, 5] = numpy.sum(x36 * (x44 + x55 * x68))
    result[2, 6] = numpy.sum(x56 * x67)
    result[2, 7] = numpy.sum(x36 * (x46 + x58 * x66))
    result[2, 8] = numpy.sum(x36 * (x60 * x66 + x65))
    result[2, 9] = numpy.sum(x27 * (3.0 * x54 + x61 * x66))
    return result


def _2center2el3d_14(ax, da, A, bx, db, B):
    """Cartesian (p|g) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 15), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = -x2 - B[0]
    x5 = bx ** (-1.0)
    x6 = ax * x1
    x7 = bx * x6 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x8 = boys(4, x7)
    x9 = 17.4934183276249
    x10 = 2.0 * x5 * x9
    x11 = x0 ** (-1.5) * x10
    x12 = x11 * x8
    x13 = x12 * x4
    x14 = ax ** (-1.0)
    x15 = x0 ** (-0.5)
    x16 = boys(3, x7)
    x17 = 0.5 * x5
    x18 = x17 * (-x12 + 2.0 * x14 * x15 * x16 * x5 * x9)
    x19 = boys(5, x7)
    x20 = x4**2
    x21 = x14 * x15
    x22 = x10 * x21
    x23 = x20 * x22
    x24 = x19 * x23
    x25 = boys(2, x7)
    x26 = x17 * (-x11 * x25 + 2.0 * x14 * x15 * x5 * x9 * boys(1, x7))
    x27 = x11 * x16
    x28 = x17 * (2.0 * x14 * x15 * x25 * x5 * x9 - x27)
    x29 = x23 * x8
    x30 = x28 + x29
    x31 = x16 * x22
    x32 = 1.5 * x5
    x33 = x32 * (x20 * x31 + x26 - x30 * x6) + x4 * (
        x4 * (x18 + x24) - x5 * (x13 - 2.0 * x14 * x15 * x16 * x4 * x5 * x9)
    )
    x34 = 0.5 / (ax + bx)
    x35 = x34 * x4 * (x30 + x5 * (2.0 * x14 * x15 * x25 * x5 * x9 - x27))
    x36 = 0.179587122125167 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**5.5)
    x37 = 8.83283915958214 * x36
    x38 = -x1 * (ax * A[1] + bx * B[1])
    x39 = -x38 - B[1]
    x40 = x12 * x39
    x41 = x5 * (2.0 * x14 * x15 * x16 * x39 * x5 * x9 - x40)
    x42 = x24 * x39
    x43 = -x39 * x5 * (x13 - 2.0 * x14 * x15 * x16 * x4 * x5 * x9) + 0.5 * x4 * (
        x41 + 2.0 * x42
    )
    x44 = x39 * x5 * (2.0 * x14 * x15 * x25 * x5 * x9 - x27)
    x45 = 0.5 * x34 * (2.0 * x29 * x39 + x44)
    x46 = 23.3694957868871 * x36
    x47 = -x1 * (ax * A[2] + bx * B[2])
    x48 = -x47 - B[2]
    x49 = x48 * x5 * (-x12 + 2.0 * x14 * x15 * x16 * x5 * x9)
    x50 = 0.5 * x49
    x51 = x4 * (x24 * x48 + x50) - x48 * x5 * (x13 - 2.0 * x14 * x15 * x16 * x4 * x5 * x9)
    x52 = x48 * x5 * (2.0 * x14 * x15 * x25 * x5 * x9 - x27)
    x53 = 0.5 * x52
    x54 = x34 * (x29 * x48 + x53)
    x55 = x39**2
    x56 = x22 * x55
    x57 = x19 * x56
    x58 = x18 + x57
    x59 = x56 * x8
    x60 = x28 + x59
    x61 = x26 + x31 * x55 - x6 * x60
    x62 = x17 * x61 + x20 * x58
    x63 = x34 * x4
    x64 = x60 * x63
    x65 = 30.169889330626 * x36
    x66 = x48 * x5 * (2.0 * x14 * x15 * x16 * x39 * x5 * x9 - x40)
    x67 = x42 * x48 + 0.5 * x66
    x68 = 4.0 * x21 * x39 * x48 * x5 * x63 * x8 * x9
    x69 = 52.2557811793745 * x36
    x70 = x48**2
    x71 = x22 * x70
    x72 = x18 + x19 * x71
    x73 = x28 + x71 * x8
    x74 = x26 + x31 * x70 - x6 * x73
    x75 = x17 * x74
    x76 = x20 * x72 + x75
    x77 = x34 * x73
    x78 = x4 * x77
    x79 = x34 * (x39 * x60 + x44)
    x80 = x39 * x58 + x41
    x81 = x3 * x4
    x82 = x34 * (x48 * x59 + x53)
    x83 = x48 * x57 + x50
    x84 = x39 * x77
    x85 = x39 * x72
    x86 = x34 * (x48 * x73 + x52)
    x87 = x48 * x72 + x49
    x88 = x32 * x61 + x39 * x80
    x89 = x3 * x37
    x90 = x39 * x83 + x66
    x91 = x3 * x46
    x92 = x55 * x72 + x75
    x93 = x39 * x87
    x94 = x32 * x74 + x48 * x87
    x95 = -x38 - A[1]
    x96 = x37 * x95
    x97 = x4 * x95
    x98 = -x47 - A[2]
    x99 = x37 * x98
    x100 = x4 * x98

    # 45 item(s)
    result[0, 0] = numpy.sum(x37 * (x3 * x33 + 4.0 * x35))
    result[0, 1] = numpy.sum(x46 * (x3 * x43 + 3.0 * x45))
    result[0, 2] = numpy.sum(x46 * (x3 * x51 + 3.0 * x54))
    result[0, 3] = numpy.sum(x65 * (x3 * x62 + 2.0 * x64))
    result[0, 4] = numpy.sum(x69 * (x3 * x67 + x68))
    result[0, 5] = numpy.sum(x65 * (x3 * x76 + 2.0 * x78))
    result[0, 6] = numpy.sum(x46 * (x79 + x80 * x81))
    result[0, 7] = numpy.sum(x69 * (x81 * x83 + x82))
    result[0, 8] = numpy.sum(x69 * (x81 * x85 + x84))
    result[0, 9] = numpy.sum(x46 * (x81 * x87 + x86))
    result[0, 10] = numpy.sum(x88 * x89)
    result[0, 11] = numpy.sum(x90 * x91)
    result[0, 12] = numpy.sum(x3 * x65 * x92)
    result[0, 13] = numpy.sum(x91 * x93)
    result[0, 14] = numpy.sum(x89 * x94)
    result[1, 0] = numpy.sum(x33 * x96)
    result[1, 1] = numpy.sum(x46 * (x35 + x43 * x95))
    result[1, 2] = numpy.sum(x46 * x51 * x95)
    result[1, 3] = numpy.sum(x65 * (2.0 * x45 + x62 * x95))
    result[1, 4] = numpy.sum(x69 * (x54 + x67 * x95))
    result[1, 5] = numpy.sum(x65 * x76 * x95)
    result[1, 6] = numpy.sum(x46 * (3.0 * x64 + x80 * x97))
    result[1, 7] = numpy.sum(x69 * (x68 + x83 * x97))
    result[1, 8] = numpy.sum(x69 * (x78 + x85 * x97))
    result[1, 9] = numpy.sum(x46 * x87 * x97)
    result[1, 10] = numpy.sum(x37 * (4.0 * x79 + x88 * x95))
    result[1, 11] = numpy.sum(x46 * (3.0 * x82 + x90 * x95))
    result[1, 12] = numpy.sum(x65 * (2.0 * x84 + x92 * x95))
    result[1, 13] = numpy.sum(x46 * (x86 + x93 * x95))
    result[1, 14] = numpy.sum(x94 * x96)
    result[2, 0] = numpy.sum(x33 * x99)
    result[2, 1] = numpy.sum(x43 * x46 * x98)
    result[2, 2] = numpy.sum(x46 * (x35 + x51 * x98))
    result[2, 3] = numpy.sum(x62 * x65 * x98)
    result[2, 4] = numpy.sum(x69 * (x45 + x67 * x98))
    result[2, 5] = numpy.sum(x65 * (2.0 * x54 + x76 * x98))
    result[2, 6] = numpy.sum(x100 * x46 * x80)
    result[2, 7] = numpy.sum(x69 * (x100 * x83 + x64))
    result[2, 8] = numpy.sum(x69 * (x100 * x85 + x68))
    result[2, 9] = numpy.sum(x46 * (x100 * x87 + 3.0 * x78))
    result[2, 10] = numpy.sum(x88 * x99)
    result[2, 11] = numpy.sum(x46 * (x79 + x90 * x98))
    result[2, 12] = numpy.sum(x65 * (2.0 * x82 + x92 * x98))
    result[2, 13] = numpy.sum(x46 * (3.0 * x84 + x93 * x98))
    result[2, 14] = numpy.sum(x37 * (4.0 * x86 + x94 * x98))
    return result


def _2center2el3d_15(ax, da, A, bx, db, B):
    """Cartesian (p|h) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 21), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = -x2 - B[0]
    x5 = bx ** (-1.0)
    x6 = ax * x1
    x7 = bx * x6 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x8 = boys(5, x7)
    x9 = 17.4934183276249
    x10 = 2.0 * x5
    x11 = x10 * x9
    x12 = x0 ** (-1.5) * x11
    x13 = x12 * x8
    x14 = x13 * x4
    x15 = ax ** (-1.0)
    x16 = x0 ** (-0.5)
    x17 = boys(4, x7)
    x18 = 0.5 * x5
    x19 = x18 * (-x13 + 2.0 * x15 * x16 * x17 * x5 * x9)
    x20 = boys(6, x7)
    x21 = x4**2
    x22 = x11 * x15 * x16
    x23 = x21 * x22
    x24 = x20 * x23
    x25 = x12 * x17
    x26 = boys(3, x7)
    x27 = x18 * (2.0 * x15 * x16 * x26 * x5 * x9 - x25)
    x28 = x23 * x8
    x29 = x27 + x28
    x30 = x12 * x26
    x31 = boys(2, x7)
    x32 = x18 * (2.0 * x15 * x16 * x31 * x5 * x9 - x30)
    x33 = x17 * x22
    x34 = x21 * x33
    x35 = x32 + x34
    x36 = 1.5 * x5
    x37 = x25 * x4
    x38 = x29 * x4 + x5 * (2.0 * x15 * x16 * x26 * x4 * x5 * x9 - x37)
    x39 = x10 * (
        x35 * x4 - x38 * x6 + x4 * x5 * (2.0 * x15 * x16 * x31 * x5 * x9 - x30)
    ) - x4 * (
        x36 * (x29 * x6 - x35)
        - x4 * (x4 * (x19 + x24) - x5 * (x14 - 2.0 * x15 * x16 * x17 * x4 * x5 * x9))
    )
    x40 = 0.5 / (ax + bx)
    x41 = x18 * (-x12 * x31 + 2.0 * x15 * x16 * x5 * x9 * boys(1, x7))
    x42 = x22 * x26
    x43 = x40 * (x36 * (x21 * x42 - x35 * x6 + x41) + x38 * x4)
    x44 = 0.179587122125167 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**6.5)
    x45 = 5.88855943972142 * x44
    x46 = -x1 * (ax * A[1] + bx * B[1])
    x47 = -x46 - B[1]
    x48 = x14 * x47
    x49 = x13 * x47
    x50 = x5 * (2.0 * x15 * x16 * x17 * x47 * x5 * x9 - x49)
    x51 = x24 * x47
    x52 = x47 * x5 * (2.0 * x15 * x16 * x31 * x5 * x9 - x30)
    x53 = x25 * x47
    x54 = x5 * (2.0 * x15 * x16 * x26 * x47 * x5 * x9 - x53)
    x55 = x28 * x47
    x56 = 0.5 * x54 + x55
    x57 = 0.5 * x36 * (2.0 * x34 * x47 + x52 - 2.0 * x56 * x6) + 0.5 * x4 * (
        x4 * (x50 + 2.0 * x51)
        + 2.0 * x5 * (2.0 * x15 * x16 * x17 * x4 * x47 * x5 * x9 - x48)
    )
    x58 = x40 * (x4 * x56 + x47 * x5 * (2.0 * x15 * x16 * x26 * x4 * x5 * x9 - x37))
    x59 = 17.6656783191643 * x44
    x60 = -x1 * (ax * A[2] + bx * B[2])
    x61 = -x60 - B[2]
    x62 = x5 * x61 * (-x13 + 2.0 * x15 * x16 * x17 * x5 * x9)
    x63 = 0.5 * x62
    x64 = x5 * x61 * (2.0 * x15 * x16 * x31 * x5 * x9 - x30)
    x65 = 0.5 * x64
    x66 = x5 * x61 * (2.0 * x15 * x16 * x26 * x5 * x9 - x25)
    x67 = 0.5 * x66
    x68 = x28 * x61 + x67
    x69 = x36 * (x34 * x61 - x6 * x68 + x65) + x4 * (
        x4 * (x24 * x61 + x63) - x5 * x61 * (x14 - 2.0 * x15 * x16 * x17 * x4 * x5 * x9)
    )
    x70 = x40 * (x4 * x68 + x5 * x61 * (2.0 * x15 * x16 * x26 * x4 * x5 * x9 - x37))
    x71 = x47**2
    x72 = x32 + x33 * x71
    x73 = x22 * x71
    x74 = x73 * x8
    x75 = x27 + x74
    x76 = x6 * x75
    x77 = x20 * x73
    x78 = x19 + x77
    x79 = x72 - x76
    x80 = x4 * (x18 * x79 + x21 * x78 + x5 * (x72 - x76))
    x81 = x41 + x42 * x71 - x6 * x72
    x82 = x40 * (x18 * x81 + x21 * x75)
    x83 = 3.0 * x82
    x84 = 26.9847693667702 * x44
    x85 = x5 * x61 * (2.0 * x15 * x16 * x17 * x47 * x5 * x9 - x49)
    x86 = 0.5 * x4 * (2.0 * x51 * x61 + x85) + x5 * x61 * (
        2.0 * x15 * x16 * x17 * x4 * x47 * x5 * x9 - x48
    )
    x87 = x5 * x61 * (2.0 * x15 * x16 * x26 * x47 * x5 * x9 - x53)
    x88 = 0.5 * x40 * (2.0 * x55 * x61 + x87)
    x89 = 46.7389915737742 * x44
    x90 = x61**2
    x91 = x32 + x33 * x90
    x92 = x22 * x90
    x93 = x27 + x8 * x92
    x94 = x6 * x93
    x95 = x19 + x20 * x92
    x96 = x21 * x95
    x97 = x91 - x94
    x98 = x18 * x97
    x99 = x4 * (x5 * (x91 - x94) + x96 + x98)
    x100 = x41 + x42 * x90 - x6 * x91
    x101 = x100 * x18
    x102 = x40 * (x101 + x21 * x93)
    x103 = 3.0 * x102
    x104 = x47 * x78 + x50
    x105 = x47 * x75 + x54
    x106 = -x105 * x6 + x47 * x72 + x52
    x107 = x104 * x21 + x106 * x18
    x108 = x4 * x40
    x109 = x105 * x108
    x110 = 26.9847693667702 * x44
    x111 = x61 * x77 + x63
    x112 = x61 * x74 + x67
    x113 = -x112 * x6 + x33 * x61 * x71 + x65
    x114 = x111 * x21 + x113 * x18
    x115 = x108 * x112
    x116 = 2.0 * x115
    x117 = 60.3397786612521 * x44
    x118 = x47 * x5 * (x91 - x94)
    x119 = 0.5 * x118 + x47 * x96
    x120 = x108 * x47 * x93
    x121 = 2.0 * x120
    x122 = x61 * x95 + x62
    x123 = x61 * x93 + x66
    x124 = -x123 * x6 + x61 * x91 + x64
    x125 = x124 * x18
    x126 = x122 * x21 + x125
    x127 = x123 * x40
    x128 = x127 * x4
    x129 = x40 * (x105 * x47 + x36 * x81)
    x130 = x104 * x47 + x36 * x79
    x131 = x3 * x4
    x132 = x40 * (x112 * x47 + x87)
    x133 = x111 * x47 + x85
    x134 = x40 * (x101 + x71 * x93)
    x135 = x71 * x95 + x98
    x136 = x127 * x47
    x137 = x122 * x47
    x138 = x40 * (x100 * x36 + x123 * x61)
    x139 = x122 * x61 + x36 * x97
    x140 = x10 * x106 + x130 * x47
    x141 = x3 * x45
    x142 = x113 * x36 + x133 * x47
    x143 = x3 * x59
    x144 = x118 + x135 * x47
    x145 = x122 * x71 + x125
    x146 = x139 * x47
    x147 = x10 * x124 + x139 * x61
    x148 = -x46 - A[1]
    x149 = x148 * x45
    x150 = 2.0 * x88
    x151 = x148 * x4
    x152 = 3.0 * x134
    x153 = -x60 - A[2]
    x154 = x153 * x45
    x155 = x153 * x4

    # 63 item(s)
    result[0, 0] = numpy.sum(x45 * (x3 * x39 + 5.0 * x43))
    result[0, 1] = numpy.sum(x59 * (x3 * x57 + 4.0 * x58))
    result[0, 2] = numpy.sum(x59 * (x3 * x69 + 4.0 * x70))
    result[0, 3] = numpy.sum(x84 * (x3 * x80 + x83))
    result[0, 4] = numpy.sum(x89 * (x3 * x86 + 3.0 * x88))
    result[0, 5] = numpy.sum(x84 * (x103 + x3 * x99))
    result[0, 6] = numpy.sum(x110 * (x107 * x3 + 2.0 * x109))
    result[0, 7] = numpy.sum(x117 * (x114 * x3 + x116))
    result[0, 8] = numpy.sum(x117 * (x119 * x3 + x121))
    result[0, 9] = numpy.sum(x110 * (x126 * x3 + 2.0 * x128))
    result[0, 10] = numpy.sum(x59 * (x129 + x130 * x131))
    result[0, 11] = numpy.sum(x89 * (x131 * x133 + x132))
    result[0, 12] = numpy.sum(x117 * (x131 * x135 + x134))
    result[0, 13] = numpy.sum(x89 * (x131 * x137 + x136))
    result[0, 14] = numpy.sum(x59 * (x131 * x139 + x138))
    result[0, 15] = numpy.sum(x140 * x141)
    result[0, 16] = numpy.sum(x142 * x143)
    result[0, 17] = numpy.sum(x144 * x3 * x84)
    result[0, 18] = numpy.sum(x110 * x145 * x3)
    result[0, 19] = numpy.sum(x143 * x146)
    result[0, 20] = numpy.sum(x141 * x147)
    result[1, 0] = numpy.sum(x149 * x39)
    result[1, 1] = numpy.sum(x59 * (x148 * x57 + x43))
    result[1, 2] = numpy.sum(x148 * x59 * x69)
    result[1, 3] = numpy.sum(x84 * (x148 * x80 + 2.0 * x58))
    result[1, 4] = numpy.sum(x89 * (x148 * x86 + x70))
    result[1, 5] = numpy.sum(x148 * x84 * x99)
    result[1, 6] = numpy.sum(x110 * (x107 * x148 + x83))
    result[1, 7] = numpy.sum(x117 * (x114 * x148 + x150))
    result[1, 8] = numpy.sum(x117 * (x102 + x119 * x148))
    result[1, 9] = numpy.sum(x110 * x126 * x148)
    result[1, 10] = numpy.sum(x59 * (4.0 * x109 + x130 * x151))
    result[1, 11] = numpy.sum(x89 * (3.0 * x115 + x133 * x151))
    result[1, 12] = numpy.sum(x117 * (x121 + x135 * x151))
    result[1, 13] = numpy.sum(x89 * (x128 + x137 * x151))
    result[1, 14] = numpy.sum(x139 * x151 * x59)
    result[1, 15] = numpy.sum(x45 * (5.0 * x129 + x140 * x148))
    result[1, 16] = numpy.sum(x59 * (4.0 * x132 + x142 * x148))
    result[1, 17] = numpy.sum(x84 * (x144 * x148 + x152))
    result[1, 18] = numpy.sum(x110 * (2.0 * x136 + x145 * x148))
    result[1, 19] = numpy.sum(x59 * (x138 + x146 * x148))
    result[1, 20] = numpy.sum(x147 * x149)
    result[2, 0] = numpy.sum(x154 * x39)
    result[2, 1] = numpy.sum(x153 * x57 * x59)
    result[2, 2] = numpy.sum(x59 * (x153 * x69 + x43))
    result[2, 3] = numpy.sum(x153 * x80 * x84)
    result[2, 4] = numpy.sum(x89 * (x153 * x86 + x58))
    result[2, 5] = numpy.sum(x84 * (x153 * x99 + 2.0 * x70))
    result[2, 6] = numpy.sum(x107 * x110 * x153)
    result[2, 7] = numpy.sum(x117 * (x114 * x153 + x82))
    result[2, 8] = numpy.sum(x117 * (x119 * x153 + x150))
    result[2, 9] = numpy.sum(x110 * (x103 + x126 * x153))
    result[2, 10] = numpy.sum(x130 * x155 * x59)
    result[2, 11] = numpy.sum(x89 * (x109 + x133 * x155))
    result[2, 12] = numpy.sum(x117 * (x116 + x135 * x155))
    result[2, 13] = numpy.sum(x89 * (3.0 * x120 + x137 * x155))
    result[2, 14] = numpy.sum(x59 * (4.0 * x128 + x139 * x155))
    result[2, 15] = numpy.sum(x140 * x154)
    result[2, 16] = numpy.sum(x59 * (x129 + x142 * x153))
    result[2, 17] = numpy.sum(x84 * (2.0 * x132 + x144 * x153))
    result[2, 18] = numpy.sum(x110 * (x145 * x153 + x152))
    result[2, 19] = numpy.sum(x59 * (4.0 * x136 + x146 * x153))
    result[2, 20] = numpy.sum(x45 * (5.0 * x138 + x147 * x153))
    return result


def _2center2el3d_20(ax, da, A, bx, db, B):
    """Cartesian (d|s) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 1), dtype=float)

    x0 = ax ** (-1.0)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = ax * bx * x2 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x4 = 17.4934183276249
    x5 = 2.0 * x0 * x4
    x6 = bx ** (-1.0)
    x7 = x1 ** (-0.5)
    x8 = x0**2 * x4 * x6 * x7 * boys(0, x3) - 0.5 * x0 * x1 ** (-1.5) * x5 * boys(1, x3)
    x9 = x2 * (ax * A[0] + bx * B[0]) - A[0]
    x10 = boys(2, x3)
    x11 = x6 * x7
    x12 = x10 * x11 * x5
    x13 = da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**1.5)
    x14 = 1.17305816955079 * x13
    x15 = x2 * (ax * A[1] + bx * B[1]) - A[1]
    x16 = 3.14159265358979
    x17 = 22.6274169979695 * x0 * x10 * x11 * x13 * x16 * x9
    x18 = x2 * (ax * A[2] + bx * B[2]) - A[2]

    # 6 item(s)
    result[0, 0] = numpy.sum(x14 * (x12 * x9**2 + x8))
    result[1, 0] = numpy.sum(x15 * x17)
    result[2, 0] = numpy.sum(x17 * x18)
    result[3, 0] = numpy.sum(x14 * (x12 * x15**2 + x8))
    result[4, 0] = numpy.sum(22.6274169979695 * x0 * x10 * x11 * x13 * x15 * x16 * x18)
    result[5, 0] = numpy.sum(x14 * (x12 * x18**2 + x8))
    return result


def _2center2el3d_21(ax, da, A, bx, db, B):
    """Cartesian (d|p) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 3), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x5 = 17.4934183276249
    x6 = 2.0 * x5
    x7 = x6 * boys(2, x4)
    x8 = ax ** (-1.0)
    x9 = bx ** (-1.0)
    x10 = x0 ** (-0.5)
    x11 = x10 * x9
    x12 = x11 * x8
    x13 = 0.5 * x12 * x7 / (ax + bx)
    x14 = -x2 - B[0]
    x15 = x14 * x8
    x16 = x11 * x15
    x17 = boys(3, x4)
    x18 = x17 * x6
    x19 = x16 * x18
    x20 = x19 * x3
    x21 = x0 ** (-1.5) * x7
    x22 = boys(1, x4)
    x23 = 0.5 * x8
    x24 = x23 * (2.0 * x10 * x14 * x22 * x5 * x8 * x9 - x15 * x21)
    x25 = da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**2.5)
    x26 = 0.179587122125167 * x25
    x27 = 13.0639452948436 * x26
    x28 = -x1 * (ax * A[1] + bx * B[1])
    x29 = -x28 - B[1]
    x30 = x21 * x8
    x31 = x23 * x29 * (2.0 * x10 * x22 * x5 * x8 * x9 - x30)
    x32 = x12 * x18
    x33 = x3**2 * x32
    x34 = -x1 * (ax * A[2] + bx * B[2])
    x35 = -x34 - B[2]
    x36 = x23 * x35 * (2.0 * x10 * x22 * x5 * x8 * x9 - x30)
    x37 = -x28 - A[1]
    x38 = x13 * x37
    x39 = 22.6274169979695 * x26
    x40 = x29 * x32
    x41 = x37 * x40
    x42 = x13 + x41
    x43 = x3 * x39
    x44 = 3.14159265358979
    x45 = 45.2548339959391 * x17 * x25 * x37 * x44
    x46 = x12 * x3
    x47 = -x34 - A[2]
    x48 = x13 * x47
    x49 = x32 * x35
    x50 = x13 + x47 * x49
    x51 = x37**2
    x52 = x47**2

    # 18 item(s)
    result[0, 0] = numpy.sum(x27 * (x13 * x3 + x24 + x3 * (x13 + x20)))
    result[0, 1] = numpy.sum(x27 * (x29 * x33 + x31))
    result[0, 2] = numpy.sum(x27 * (x33 * x35 + x36))
    result[1, 0] = numpy.sum(x39 * (x20 * x37 + x38))
    result[1, 1] = numpy.sum(x42 * x43)
    result[1, 2] = numpy.sum(x35 * x45 * x46)
    result[2, 0] = numpy.sum(x39 * (x20 * x47 + x48))
    result[2, 1] = numpy.sum(45.2548339959391 * x17 * x25 * x29 * x44 * x46 * x47)
    result[2, 2] = numpy.sum(x43 * x50)
    result[3, 0] = numpy.sum(x27 * (x19 * x51 + x24))
    result[3, 1] = numpy.sum(x27 * (x31 + x37 * x42 + x38))
    result[3, 2] = numpy.sum(x27 * (x36 + x49 * x51))
    result[4, 0] = numpy.sum(x16 * x45 * x47)
    result[4, 1] = numpy.sum(x39 * (x41 * x47 + x48))
    result[4, 2] = numpy.sum(x37 * x39 * x50)
    result[5, 0] = numpy.sum(x27 * (x19 * x52 + x24))
    result[5, 1] = numpy.sum(x27 * (x31 + x40 * x52))
    result[5, 2] = numpy.sum(x27 * (x36 + x47 * x50 + x48))
    return result


def _2center2el3d_22(ax, da, A, bx, db, B):
    """Cartesian (d|d) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 6), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = bx * x1
    x5 = ax * x4 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x6 = boys(3, x5)
    x7 = x0 ** (-1.5)
    x8 = bx ** (-1.0)
    x9 = 17.4934183276249
    x10 = x8 * x9
    x11 = 2.0 * x10
    x12 = x11 * x7
    x13 = ax ** (-1.0)
    x14 = x0 ** (-0.5)
    x15 = boys(2, x5)
    x16 = 0.5 * x8
    x17 = x16 * (-x12 * x6 + 2.0 * x13 * x14 * x15 * x8 * x9)
    x18 = boys(4, x5)
    x19 = -x2 - B[0]
    x20 = x19**2
    x21 = x13 * x14
    x22 = x11 * x21
    x23 = x20 * x22
    x24 = x17 + x18 * x23
    x25 = x24 * x3
    x26 = x19 * x6
    x27 = 0.5 / (ax + bx)
    x28 = x10 * x21
    x29 = 4.0 * x27 * x28
    x30 = x26 * x29
    x31 = boys(1, x5)
    x32 = x16 * (-x12 * x31 + 2.0 * x13 * x14 * x8 * x9 * boys(0, x5))
    x33 = x16 * (-x12 * x15 + 2.0 * x13 * x14 * x31 * x8 * x9)
    x34 = x15 * x22
    x35 = 0.5 * x13
    x36 = x35 * (x20 * x34 + x32 - x4 * (x23 * x6 + x33))
    x37 = 2.0 * x27
    x38 = x28 * x37
    x39 = x15 * x38
    x40 = x22 * x3
    x41 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**3.5)
    x42 = 15.084944665313 * x41
    x43 = -x1 * (ax * A[1] + bx * B[1])
    x44 = -x43 - B[1]
    x45 = x44 * x6
    x46 = x38 * x45
    x47 = x18 * x19
    x48 = x40 * x47
    x49 = x44 * x48
    x50 = 2.0 * x13 * x26 * x7 * x9
    x51 = x35 * x44 * (2.0 * x13 * x14 * x15 * x19 * x8 * x9 - x50)
    x52 = 26.1278905896872 * x41
    x53 = -x1 * (ax * A[2] + bx * B[2])
    x54 = -x53 - B[2]
    x55 = x54 * x6
    x56 = x38 * x55
    x57 = x48 * x54
    x58 = x35 * x54 * (2.0 * x13 * x14 * x15 * x19 * x8 * x9 - x50)
    x59 = x3**2
    x60 = x44**2
    x61 = x22 * x60
    x62 = x17 + x18 * x61
    x63 = x35 * (x32 + x34 * x60 - x4 * (x33 + x6 * x61))
    x64 = 2.0 * x13 * x35 * x54 * x9 * (x14 * x15 * x44 * x8 - x45 * x7)
    x65 = x18 * x44 * x54
    x66 = x54**2
    x67 = x22 * x66
    x68 = x17 + x18 * x67
    x69 = x35 * (x32 + x34 * x66 - x4 * (x33 + x6 * x67))
    x70 = -x43 - A[1]
    x71 = x22 * x70
    x72 = x27 * (x39 + x45 * x71)
    x73 = x26 * x38
    x74 = x44 * x47
    x75 = x71 * x74
    x76 = x73 + x75
    x77 = 45.2548339959391 * x41
    x78 = x56 * x70
    x79 = x62 * x70
    x80 = x29 * x45
    x81 = x79 + x80
    x82 = x3 * x52
    x83 = x56 + x65 * x71
    x84 = x3 * x77
    x85 = -x53 - A[2]
    x86 = x46 * x85
    x87 = x22 * x85
    x88 = x27 * (x39 + x55 * x87)
    x89 = x47 * x54
    x90 = x73 + x87 * x89
    x91 = x46 + x65 * x87
    x92 = x29 * x55 + x68 * x85
    x93 = x70**2
    x94 = x52 * x70
    x95 = x73 * x85
    x96 = x85**2

    # 36 item(s)
    result[0, 0] = numpy.sum(x42 * (x3 * (x25 + x30) + x36 + x37 * (x26 * x40 + x39)))
    result[0, 1] = numpy.sum(x52 * (x3 * x46 + x3 * (x46 + x49) + x51))
    result[0, 2] = numpy.sum(x52 * (x3 * x56 + x3 * (x56 + x57) + x58))
    result[0, 3] = numpy.sum(x42 * (x59 * x62 + x63))
    result[0, 4] = numpy.sum(x52 * (x22 * x59 * x65 + x64))
    result[0, 5] = numpy.sum(x42 * (x59 * x68 + x69))
    result[1, 0] = numpy.sum(x52 * x70 * (x25 + x30))
    result[1, 1] = numpy.sum(x77 * (x3 * x76 + x72))
    result[1, 2] = numpy.sum(x77 * (x57 * x70 + x78))
    result[1, 3] = numpy.sum(x81 * x82)
    result[1, 4] = numpy.sum(x83 * x84)
    result[1, 5] = numpy.sum(x68 * x70 * x82)
    result[2, 0] = numpy.sum(x52 * x85 * (x25 + x30))
    result[2, 1] = numpy.sum(x77 * (x49 * x85 + x86))
    result[2, 2] = numpy.sum(x77 * (x3 * x90 + x88))
    result[2, 3] = numpy.sum(x62 * x82 * x85)
    result[2, 4] = numpy.sum(x84 * x91)
    result[2, 5] = numpy.sum(x82 * x92)
    result[3, 0] = numpy.sum(x42 * (x24 * x93 + x36))
    result[3, 1] = numpy.sum(x52 * (x51 + x70 * x73 + x70 * x76))
    result[3, 2] = numpy.sum(x52 * (x22 * x89 * x93 + x58))
    result[3, 3] = numpy.sum(x42 * (x63 + x70 * x81 + 2.0 * x72))
    result[3, 4] = numpy.sum(x52 * (x64 + x70 * x83 + x78))
    result[3, 5] = numpy.sum(x42 * (x68 * x93 + x69))
    result[4, 0] = numpy.sum(x24 * x85 * x94)
    result[4, 1] = numpy.sum(x77 * (x75 * x85 + x95))
    result[4, 2] = numpy.sum(x70 * x77 * x90)
    result[4, 3] = numpy.sum(x52 * x85 * (x79 + x80))
    result[4, 4] = numpy.sum(x77 * (x70 * x91 + x88))
    result[4, 5] = numpy.sum(x92 * x94)
    result[5, 0] = numpy.sum(x42 * (x24 * x96 + x36))
    result[5, 1] = numpy.sum(x52 * (x22 * x74 * x96 + x51))
    result[5, 2] = numpy.sum(x52 * (x58 + x85 * x90 + x95))
    result[5, 3] = numpy.sum(x42 * (x62 * x96 + x63))
    result[5, 4] = numpy.sum(x52 * (x64 + x85 * x91 + x86))
    result[5, 5] = numpy.sum(x42 * (x69 + x85 * x92 + 2.0 * x88))
    return result


def _2center2el3d_23(ax, da, A, bx, db, B):
    """Cartesian (d|f) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 10), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = bx ** (-1.0)
    x5 = -x2 - B[0]
    x6 = bx * x1
    x7 = ax * x6 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x8 = boys(4, x7)
    x9 = x0 ** (-1.5)
    x10 = 17.4934183276249
    x11 = x10 * x4
    x12 = 2.0 * x11
    x13 = x12 * x9
    x14 = x13 * x8
    x15 = ax ** (-1.0)
    x16 = x0 ** (-0.5)
    x17 = boys(3, x7)
    x18 = 0.5 * x4
    x19 = x18 * (2.0 * x10 * x15 * x16 * x17 * x4 - x14)
    x20 = boys(5, x7)
    x21 = x5**2
    x22 = x15 * x16
    x23 = x12 * x22
    x24 = x21 * x23
    x25 = x20 * x24
    x26 = x5 * (x19 + x25 + x4 * (2.0 * x10 * x15 * x16 * x17 * x4 - x14))
    x27 = x26 * x3
    x28 = 0.5 / (ax + bx)
    x29 = x13 * x17
    x30 = boys(2, x7)
    x31 = x18 * (2.0 * x10 * x15 * x16 * x30 * x4 - x29)
    x32 = x24 * x8
    x33 = x31 + x32
    x34 = x28 * x33
    x35 = 3.0 * x34
    x36 = x13 * x30
    x37 = boys(1, x7)
    x38 = x18 * (2.0 * x10 * x15 * x16 * x37 * x4 - x36)
    x39 = x17 * x23
    x40 = x21 * x39
    x41 = 0.5 * x15
    x42 = (
        x41
        * x5
        * (
            x38
            + x4 * (2.0 * x10 * x15 * x16 * x37 * x4 - x36)
            + x40
            - x6 * (x33 + x4 * (2.0 * x10 * x15 * x16 * x30 * x4 - x29))
        )
    )
    x43 = x11 * x22
    x44 = x17 * x43
    x45 = x44 * x5
    x46 = 4.0 * x28
    x47 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**4.5)
    x48 = 13.4923846833851 * x47
    x49 = -x1 * (ax * A[1] + bx * B[1])
    x50 = -x49 - B[1]
    x51 = x4 * x50 * (2.0 * x10 * x15 * x16 * x17 * x4 - x14)
    x52 = x25 * x50 + 0.5 * x51
    x53 = x3 * x52
    x54 = x50 * x8
    x55 = x5 * x54
    x56 = x43 * x46
    x57 = x55 * x56
    x58 = x4 * x50 * (2.0 * x10 * x15 * x16 * x37 * x4 - x36)
    x59 = x4 * x50 * (2.0 * x10 * x15 * x16 * x30 * x4 - x29)
    x60 = 0.5 * x41 * (2.0 * x40 * x50 + x58 - x6 * (2.0 * x32 * x50 + x59))
    x61 = 2.0 * x28
    x62 = x44 * x61
    x63 = x50 * x62
    x64 = x3 * x5
    x65 = x23 * x64
    x66 = 30.169889330626 * x47
    x67 = -x1 * (ax * A[2] + bx * B[2])
    x68 = -x67 - B[2]
    x69 = x4 * x68 * (2.0 * x10 * x15 * x16 * x17 * x4 - x14)
    x70 = 0.5 * x69
    x71 = x25 * x68 + x70
    x72 = x3 * x71
    x73 = x68 * x8
    x74 = x5 * x73
    x75 = x56 * x74
    x76 = x4 * x68 * (2.0 * x10 * x15 * x16 * x37 * x4 - x36)
    x77 = 0.5 * x76
    x78 = x4 * x68 * (2.0 * x10 * x15 * x16 * x30 * x4 - x29)
    x79 = 0.5 * x78
    x80 = x41 * (x40 * x68 - x6 * (x32 * x68 + x79) + x77)
    x81 = x62 * x68
    x82 = x50**2
    x83 = x23 * x82
    x84 = x8 * x83
    x85 = x31 + x84
    x86 = x28 * x85
    x87 = x20 * x83
    x88 = x19 + x87
    x89 = x64 * x88
    x90 = x38 + x39 * x82
    x91 = x5 * x6
    x92 = x41 * (x5 * x90 - x85 * x91)
    x93 = x54 * x68
    x94 = x43 * x61
    x95 = x93 * x94
    x96 = x20 * x50 * x68
    x97 = 2.0 * x10 * x15 * x41 * x68 * (x16 * x17 * x4 * x5 * x50 - x55 * x9)
    x98 = 52.2557811793745 * x47
    x99 = x68**2
    x100 = x23 * x99
    x101 = x100 * x8 + x31
    x102 = x101 * x28
    x103 = x100 * x20 + x19
    x104 = x103 * x64
    x105 = x38 + x39 * x99
    x106 = x41 * (-x101 * x91 + x105 * x5)
    x107 = x3**2
    x108 = x50 * x88 + x51
    x109 = x41 * (x50 * x90 + x58 - x6 * (x50 * x85 + x59))
    x110 = x68 * x87 + x70
    x111 = x41 * (x39 * x68 * x82 - x6 * (x68 * x84 + x79) + x77)
    x112 = x41 * x50 * (-x101 * x6 + x105)
    x113 = x103 * x50
    x114 = x103 * x68 + x69
    x115 = x41 * (x105 * x68 - x6 * (x101 * x68 + x78) + x76)
    x116 = -x49 - A[1]
    x117 = 23.3694957868871 * x47
    x118 = x116 * x52
    x119 = x118 + x34
    x120 = x45 * x61
    x121 = x116 * x23
    x122 = x61 * (x120 + x121 * x55)
    x123 = x44 * x46
    x124 = x28 * (x116 * x85 + x123 * x50)
    x125 = x5 * x88
    x126 = x116 * x125
    x127 = x126 + x57
    x128 = x28 * (x121 * x93 + x81)
    x129 = x74 * x94
    x130 = x5 * x96
    x131 = x121 * x130 + x129
    x132 = 90.5096679918781 * x47
    x133 = x102 * x116
    x134 = x108 * x116
    x135 = 3.0 * x86
    x136 = x134 + x135
    x137 = x117 * x3
    x138 = x56 * x93
    x139 = x110 * x116 + x138
    x140 = x3 * x98
    x141 = x102 + x113 * x116
    x142 = -x67 - A[2]
    x143 = x142 * x57
    x144 = x142 * x71 + x34
    x145 = x142 * x23
    x146 = x28 * (x120 + x145 * x74)
    x147 = 2.0 * x146
    x148 = x142 * x86
    x149 = x28 * (x145 * x93 + x63)
    x150 = x55 * x94
    x151 = x130 * x145 + x150
    x152 = x28 * (x101 * x142 + x123 * x68)
    x153 = x103 * x5
    x154 = x142 * x153 + x75
    x155 = x110 * x142 + x86
    x156 = x113 * x142 + x138
    x157 = 3.0 * x102 + x114 * x142
    x158 = x116**2
    x159 = x116 * x117
    x160 = x142 * x34
    x161 = x116 * x98
    x162 = 2.0 * x149
    x163 = x142**2

    # 60 item(s)
    result[0, 0] = numpy.sum(
        x48 * (3.0 * x28 * (x3 * x33 + x45 * x46) + x3 * (x27 + x35) + x42)
    )
    result[0, 1] = numpy.sum(x66 * (x3 * (x53 + x57) + x60 + x61 * (x54 * x65 + x63)))
    result[0, 2] = numpy.sum(x66 * (x3 * (x72 + x75) + x61 * (x65 * x73 + x81) + x80))
    result[0, 3] = numpy.sum(x66 * (x3 * x86 + x3 * (x86 + x89) + x92))
    result[0, 4] = numpy.sum(x98 * (x3 * x95 + x3 * (x65 * x96 + x95) + x97))
    result[0, 5] = numpy.sum(x66 * (x102 * x3 + x106 + x3 * (x102 + x104)))
    result[0, 6] = numpy.sum(x48 * (x107 * x108 + x109))
    result[0, 7] = numpy.sum(x66 * (x107 * x110 + x111))
    result[0, 8] = numpy.sum(x66 * (x107 * x113 + x112))
    result[0, 9] = numpy.sum(x48 * (x107 * x114 + x115))
    result[1, 0] = numpy.sum(x116 * x117 * (x27 + x35))
    result[1, 1] = numpy.sum(x98 * (x119 * x3 + x122))
    result[1, 2] = numpy.sum(x116 * x98 * (x72 + x75))
    result[1, 3] = numpy.sum(x98 * (x124 + x127 * x3))
    result[1, 4] = numpy.sum(x132 * (x128 + x131 * x3))
    result[1, 5] = numpy.sum(x98 * (x104 * x116 + x133))
    result[1, 6] = numpy.sum(x136 * x137)
    result[1, 7] = numpy.sum(x139 * x140)
    result[1, 8] = numpy.sum(x140 * x141)
    result[1, 9] = numpy.sum(x114 * x116 * x137)
    result[2, 0] = numpy.sum(x117 * x142 * (x27 + x35))
    result[2, 1] = numpy.sum(x98 * (x142 * x53 + x143))
    result[2, 2] = numpy.sum(x98 * (x144 * x3 + x147))
    result[2, 3] = numpy.sum(x98 * (x142 * x89 + x148))
    result[2, 4] = numpy.sum(x132 * (x149 + x151 * x3))
    result[2, 5] = numpy.sum(x98 * (x152 + x154 * x3))
    result[2, 6] = numpy.sum(x108 * x137 * x142)
    result[2, 7] = numpy.sum(x140 * x155)
    result[2, 8] = numpy.sum(x140 * x156)
    result[2, 9] = numpy.sum(x137 * x157)
    result[3, 0] = numpy.sum(x48 * (x158 * x26 + x42))
    result[3, 1] = numpy.sum(x66 * (x116 * x119 + x116 * x34 + x60))
    result[3, 2] = numpy.sum(x66 * (x158 * x71 + x80))
    result[3, 3] = numpy.sum(x66 * (x116 * x127 + x122 + x92))
    result[3, 4] = numpy.sum(x98 * (x116 * x129 + x116 * x131 + x97))
    result[3, 5] = numpy.sum(x66 * (x106 + x153 * x158))
    result[3, 6] = numpy.sum(x48 * (x109 + x116 * x136 + 3.0 * x124))
    result[3, 7] = numpy.sum(x66 * (x111 + x116 * x139 + 2.0 * x128))
    result[3, 8] = numpy.sum(x66 * (x112 + x116 * x141 + x133))
    result[3, 9] = numpy.sum(x48 * (x114 * x158 + x115))
    result[4, 0] = numpy.sum(x142 * x159 * x26)
    result[4, 1] = numpy.sum(x98 * (x118 * x142 + x160))
    result[4, 2] = numpy.sum(x144 * x161)
    result[4, 3] = numpy.sum(x98 * (x126 * x142 + x143))
    result[4, 4] = numpy.sum(x132 * (x116 * x151 + x146))
    result[4, 5] = numpy.sum(x154 * x161)
    result[4, 6] = numpy.sum(x117 * x142 * (x134 + x135))
    result[4, 7] = numpy.sum(x98 * (x116 * x155 + x162))
    result[4, 8] = numpy.sum(x98 * (x116 * x156 + x152))
    result[4, 9] = numpy.sum(x157 * x159)
    result[5, 0] = numpy.sum(x48 * (x163 * x26 + x42))
    result[5, 1] = numpy.sum(x66 * (x163 * x52 + x60))
    result[5, 2] = numpy.sum(x66 * (x142 * x144 + x160 + x80))
    result[5, 3] = numpy.sum(x66 * (x125 * x163 + x92))
    result[5, 4] = numpy.sum(x98 * (x142 * x150 + x142 * x151 + x97))
    result[5, 5] = numpy.sum(x66 * (x106 + x142 * x154 + x147))
    result[5, 6] = numpy.sum(x48 * (x108 * x163 + x109))
    result[5, 7] = numpy.sum(x66 * (x111 + x142 * x155 + x148))
    result[5, 8] = numpy.sum(x66 * (x112 + x142 * x156 + x162))
    result[5, 9] = numpy.sum(x48 * (x115 + x142 * x157 + 3.0 * x152))
    return result


def _2center2el3d_24(ax, da, A, bx, db, B):
    """Cartesian (d|g) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 15), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = -x2 - B[0]
    x5 = bx ** (-1.0)
    x6 = ax * x1
    x7 = bx * x6 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x8 = boys(5, x7)
    x9 = 17.4934183276249
    x10 = x5 * x9
    x11 = 2.0 * x10
    x12 = x0 ** (-1.5) * x11
    x13 = x12 * x8
    x14 = x13 * x4
    x15 = ax ** (-1.0)
    x16 = x0 ** (-0.5)
    x17 = boys(4, x7)
    x18 = 0.5 * x5
    x19 = x18 * (-x13 + 2.0 * x15 * x16 * x17 * x5 * x9)
    x20 = boys(6, x7)
    x21 = x4**2
    x22 = x15 * x16
    x23 = x11 * x22
    x24 = x21 * x23
    x25 = x20 * x24
    x26 = x12 * x17
    x27 = boys(3, x7)
    x28 = x18 * (2.0 * x15 * x16 * x27 * x5 * x9 - x26)
    x29 = x24 * x8
    x30 = x28 + x29
    x31 = x12 * x27
    x32 = boys(2, x7)
    x33 = x18 * (2.0 * x15 * x16 * x32 * x5 * x9 - x31)
    x34 = x17 * x23
    x35 = x21 * x34
    x36 = x33 + x35
    x37 = 1.5 * x5
    x38 = -x37 * (x30 * x6 - x36) + x4 * (
        x4 * (x19 + x25) - x5 * (x14 - 2.0 * x15 * x16 * x17 * x4 * x5 * x9)
    )
    x39 = x3 * x38
    x40 = 0.5 / (ax + bx)
    x41 = x26 * x4
    x42 = x30 * x4 + x5 * (2.0 * x15 * x16 * x27 * x4 * x5 * x9 - x41)
    x43 = x40 * x42
    x44 = 4.0 * x43
    x45 = x31 * x4
    x46 = boys(1, x7)
    x47 = x18 * (-x12 * x46 + 2.0 * x15 * x16 * x5 * x9 * boys(0, x7))
    x48 = x18 * (-x12 * x32 + 2.0 * x15 * x16 * x46 * x5 * x9)
    x49 = x23 * x27
    x50 = x21 * x49 + x48
    x51 = x23 * x32
    x52 = bx * x1
    x53 = 0.5 * x15
    x54 = x53 * (
        x37 * (x21 * x51 + x47 - x50 * x6)
        + x4 * (x36 * x4 + x5 * (2.0 * x15 * x16 * x32 * x4 * x5 * x9 - x45))
        + x52 * (x37 * (x36 * x6 - x50) - x4 * x42)
    )
    x55 = x36 * x40
    x56 = 4.0 * x40
    x57 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**5.5)
    x58 = 10.1992841329868 * x57
    x59 = -x1 * (ax * A[1] + bx * B[1])
    x60 = -x59 - B[1]
    x61 = x4 * x60
    x62 = x13 * x60
    x63 = x5 * (2.0 * x15 * x16 * x17 * x5 * x60 * x9 - x62)
    x64 = x25 * x60
    x65 = 0.5 * x4 * (x63 + 2.0 * x64) - x5 * (
        x13 * x61 - 2.0 * x15 * x16 * x17 * x4 * x5 * x60 * x9
    )
    x66 = x3 * x65
    x67 = x26 * x60
    x68 = x5 * (2.0 * x15 * x16 * x27 * x5 * x60 * x9 - x67)
    x69 = x29 * x60
    x70 = 0.5 * x68 + x69
    x71 = x40 * x70
    x72 = 3.0 * x71
    x73 = x31 * x60
    x74 = x5 * (2.0 * x15 * x16 * x32 * x5 * x60 * x9 - x73)
    x75 = x35 * x60
    x76 = (
        0.5
        * x53
        * (
            x4 * (x74 + 2.0 * x75)
            + 2.0 * x5 * (2.0 * x15 * x16 * x32 * x4 * x5 * x60 * x9 - x31 * x61)
            - 2.0
            * x52
            * (x4 * x70 + x5 * (2.0 * x15 * x16 * x27 * x4 * x5 * x60 * x9 - x26 * x61))
        )
    )
    x77 = x10 * x17 * x22
    x78 = x56 * x77
    x79 = x61 * x78
    x80 = 3.0 * x40
    x81 = 26.9847693667702 * x57
    x82 = -x1 * (ax * A[2] + bx * B[2])
    x83 = -x82 - B[2]
    x84 = x5 * x83 * (-x13 + 2.0 * x15 * x16 * x17 * x5 * x9)
    x85 = 0.5 * x84
    x86 = x4 * (x25 * x83 + x85) - x5 * x83 * (x14 - 2.0 * x15 * x16 * x17 * x4 * x5 * x9)
    x87 = x3 * x86
    x88 = x5 * x83 * (2.0 * x15 * x16 * x27 * x5 * x9 - x26)
    x89 = 0.5 * x88
    x90 = x29 * x83 + x89
    x91 = x40 * x90
    x92 = 3.0 * x91
    x93 = x5 * x83 * (2.0 * x15 * x16 * x32 * x5 * x9 - x31)
    x94 = 0.5 * x93
    x95 = x53 * (
        x4 * (x35 * x83 + x94)
        + x5 * x83 * (2.0 * x15 * x16 * x32 * x4 * x5 * x9 - x45)
        - x52 * (x4 * x90 + x5 * x83 * (2.0 * x15 * x16 * x27 * x4 * x5 * x9 - x41))
    )
    x96 = x4 * x83
    x97 = x78 * x96
    x98 = x60**2
    x99 = x23 * x98
    x100 = x20 * x99
    x101 = x100 + x19
    x102 = x8 * x99
    x103 = x102 + x28
    x104 = x33 + x34 * x98
    x105 = -x103 * x6 + x104
    x106 = x101 * x21 + x105 * x18
    x107 = x106 * x3
    x108 = 2.0 * x40
    x109 = x103 * x4
    x110 = x108 * x109
    x111 = x48 + x49 * x98
    x112 = -x111 * x6 + x47 + x51 * x98
    x113 = -x104 * x6 + x111
    x114 = x53 * (x104 * x21 + x112 * x18 - x52 * (x103 * x21 + x113 * x18))
    x115 = x104 * x40
    x116 = x3 * x4
    x117 = 34.8371874529163 * x57
    x118 = x5 * x83 * (2.0 * x15 * x16 * x17 * x5 * x60 * x9 - x62)
    x119 = 0.5 * x118 + x64 * x83
    x120 = x8 * x83
    x121 = x120 * x61
    x122 = x10 * x121 * x22 * x56
    x123 = x5 * x83 * (2.0 * x15 * x16 * x32 * x5 * x60 * x9 - x73)
    x124 = x5 * x83 * (2.0 * x15 * x16 * x27 * x5 * x60 * x9 - x67)
    x125 = 0.5 * x53 * (x123 - x52 * (x124 + 2.0 * x69 * x83) + 2.0 * x75 * x83)
    x126 = x108 * x77
    x127 = x60 * x83
    x128 = x116 * x60
    x129 = 60.3397786612521 * x57
    x130 = x83**2
    x131 = x130 * x23
    x132 = x131 * x20 + x19
    x133 = x131 * x8 + x28
    x134 = x130 * x34 + x33
    x135 = -x133 * x6 + x134
    x136 = x135 * x18
    x137 = x132 * x21 + x136
    x138 = x137 * x3
    x139 = x133 * x4
    x140 = x130 * x49 + x48
    x141 = x130 * x51 - x140 * x6 + x47
    x142 = x141 * x18
    x143 = -x134 * x6 + x140
    x144 = x143 * x18
    x145 = x53 * (x134 * x21 + x142 - x52 * (x133 * x21 + x144))
    x146 = x134 * x40
    x147 = x103 * x60 + x68
    x148 = x147 * x40
    x149 = x101 * x60 + x63
    x150 = x116 * x149
    x151 = x104 * x60 + x74
    x152 = x4 * x52
    x153 = x53 * (-x147 * x152 + x151 * x4)
    x154 = x102 * x83 + x89
    x155 = x154 * x40
    x156 = x100 * x83 + x85
    x157 = x34 * x83 * x98 + x94
    x158 = x53 * (-x152 * x154 + x157 * x4)
    x159 = x133 * x60
    x160 = x159 * x40
    x161 = x53 * (x134 * x4 * x60 - x152 * x159)
    x162 = x133 * x83 + x88
    x163 = x162 * x40
    x164 = x132 * x83 + x84
    x165 = x116 * x164
    x166 = x134 * x83 + x93
    x167 = x53 * (-x152 * x162 + x166 * x4)
    x168 = x3**2
    x169 = x105 * x37 + x149 * x60
    x170 = x53 * (x112 * x37 + x151 * x60 - x52 * (x113 * x37 + x147 * x60))
    x171 = x118 + x156 * x60
    x172 = x53 * (x123 + x157 * x60 - x52 * (x124 + x154 * x60))
    x173 = x132 * x98 + x136
    x174 = x53 * (x134 * x98 + x142 - x52 * (x133 * x98 + x144))
    x175 = x53 * x60 * (-x162 * x52 + x166)
    x176 = x164 * x60
    x177 = x135 * x37 + x164 * x83
    x178 = x53 * (x141 * x37 + x166 * x83 - x52 * (x143 * x37 + x162 * x83))
    x179 = -x59 - A[1]
    x180 = 17.6656783191643 * x57
    x181 = x179 * x65
    x182 = x181 + x43
    x183 = x179 * x70 + x55
    x184 = 46.7389915737742 * x57
    x185 = x106 * x179
    x186 = 2.0 * x71
    x187 = x185 + x186
    x188 = x179 * x4
    x189 = x103 * x188 + x79
    x190 = 60.3397786612521 * x57
    x191 = x119 * x179 + x91
    x192 = x121 * x23
    x193 = x108 * (x126 * x96 + x179 * x192)
    x194 = 104.511562358749 * x57
    x195 = x133 * x188
    x196 = x40 * (3.0 * x115 + x147 * x179)
    x197 = x149 * x188
    x198 = x109 * x80
    x199 = x197 + x198
    x200 = x127 * x78
    x201 = x40 * (x154 * x179 + x200)
    x202 = x122 + x156 * x188
    x203 = x40 * (x146 + x159 * x179)
    x204 = x132 * x61
    x205 = x139 * x40 + x179 * x204
    x206 = x163 * x179
    x207 = x169 * x179
    x208 = 4.0 * x148
    x209 = x207 + x208
    x210 = x180 * x3
    x211 = 3.0 * x155 + x171 * x179
    x212 = x184 * x3
    x213 = x108 * x159 + x173 * x179
    x214 = x190 * x3
    x215 = x163 + x176 * x179
    x216 = -x82 - A[2]
    x217 = x216 * x86 + x43
    x218 = x40 * (x216 * x90 + x55)
    x219 = x119 * x216 + x71
    x220 = x108 * (x126 * x61 + x192 * x216)
    x221 = x137 * x216 + 2.0 * x91
    x222 = x40 * (x139 * x216 + x97)
    x223 = x148 * x216
    x224 = x40 * (x115 + x154 * x216)
    x225 = x109 * x40
    x226 = x216 * x4
    x227 = x156 * x226 + x225
    x228 = x40 * (x159 * x216 + x200)
    x229 = x122 + x204 * x216
    x230 = x40 * (3.0 * x146 + x162 * x216)
    x231 = x139 * x80 + x164 * x226
    x232 = x148 + x171 * x216
    x233 = 2.0 * x155 + x173 * x216
    x234 = x159 * x80 + x176 * x216
    x235 = 4.0 * x163 + x177 * x216
    x236 = x179**2
    x237 = x179 * x180
    x238 = x216 * x43
    x239 = x179 * x184
    x240 = x216**2

    # 90 item(s)
    result[0, 0] = numpy.sum(
        x58 * (x3 * (x39 + x44) + x54 + x56 * (x3 * x42 + 3.0 * x55))
    )
    result[0, 1] = numpy.sum(x81 * (x3 * (x66 + x72) + x76 + x80 * (x3 * x70 + x79)))
    result[0, 2] = numpy.sum(x81 * (x3 * (x87 + x92) + x80 * (x3 * x90 + x97) + x95))
    result[0, 3] = numpy.sum(
        x117 * (x108 * (x103 * x116 + x115) + x114 + x3 * (x107 + x110))
    )
    result[0, 4] = numpy.sum(
        x129 * (x108 * (x120 * x128 * x23 + x126 * x127) + x125 + x3 * (x119 * x3 + x122))
    )
    result[0, 5] = numpy.sum(
        x117 * (x108 * (x116 * x133 + x146) + x145 + x3 * (x108 * x139 + x138))
    )
    result[0, 6] = numpy.sum(x81 * (x148 * x3 + x153 + x3 * (x148 + x150)))
    result[0, 7] = numpy.sum(x129 * (x155 * x3 + x158 + x3 * (x116 * x156 + x155)))
    result[0, 8] = numpy.sum(x129 * (x160 * x3 + x161 + x3 * (x128 * x132 + x160)))
    result[0, 9] = numpy.sum(x81 * (x163 * x3 + x167 + x3 * (x163 + x165)))
    result[0, 10] = numpy.sum(x58 * (x168 * x169 + x170))
    result[0, 11] = numpy.sum(x81 * (x168 * x171 + x172))
    result[0, 12] = numpy.sum(x117 * (x168 * x173 + x174))
    result[0, 13] = numpy.sum(x81 * (x168 * x176 + x175))
    result[0, 14] = numpy.sum(x58 * (x168 * x177 + x178))
    result[1, 0] = numpy.sum(x179 * x180 * (x39 + x44))
    result[1, 1] = numpy.sum(x184 * (x182 * x3 + x183 * x80))
    result[1, 2] = numpy.sum(x179 * x184 * (x87 + x92))
    result[1, 3] = numpy.sum(x190 * (x108 * x189 + x187 * x3))
    result[1, 4] = numpy.sum(x194 * (x191 * x3 + x193))
    result[1, 5] = numpy.sum(x190 * (x108 * x195 + x138 * x179))
    result[1, 6] = numpy.sum(x184 * (x196 + x199 * x3))
    result[1, 7] = numpy.sum(x194 * (x201 + x202 * x3))
    result[1, 8] = numpy.sum(x194 * (x203 + x205 * x3))
    result[1, 9] = numpy.sum(x184 * (x165 * x179 + x206))
    result[1, 10] = numpy.sum(x209 * x210)
    result[1, 11] = numpy.sum(x211 * x212)
    result[1, 12] = numpy.sum(x213 * x214)
    result[1, 13] = numpy.sum(x212 * x215)
    result[1, 14] = numpy.sum(x177 * x179 * x210)
    result[2, 0] = numpy.sum(x180 * x216 * (x39 + x44))
    result[2, 1] = numpy.sum(x184 * x216 * (x66 + x72))
    result[2, 2] = numpy.sum(x184 * (x217 * x3 + 3.0 * x218))
    result[2, 3] = numpy.sum(x190 * x216 * (x107 + x110))
    result[2, 4] = numpy.sum(x194 * (x219 * x3 + x220))
    result[2, 5] = numpy.sum(x190 * (x221 * x3 + 2.0 * x222))
    result[2, 6] = numpy.sum(x184 * (x150 * x216 + x223))
    result[2, 7] = numpy.sum(x194 * (x224 + x227 * x3))
    result[2, 8] = numpy.sum(x194 * (x228 + x229 * x3))
    result[2, 9] = numpy.sum(x184 * (x230 + x231 * x3))
    result[2, 10] = numpy.sum(x169 * x210 * x216)
    result[2, 11] = numpy.sum(x212 * x232)
    result[2, 12] = numpy.sum(x214 * x233)
    result[2, 13] = numpy.sum(x212 * x234)
    result[2, 14] = numpy.sum(x210 * x235)
    result[3, 0] = numpy.sum(x58 * (x236 * x38 + x54))
    result[3, 1] = numpy.sum(x81 * (x179 * x182 + x179 * x43 + x76))
    result[3, 2] = numpy.sum(x81 * (x236 * x86 + x95))
    result[3, 3] = numpy.sum(x117 * (x108 * x183 + x114 + x179 * x187))
    result[3, 4] = numpy.sum(x129 * (x125 + x179 * x191 + x179 * x91))
    result[3, 5] = numpy.sum(x117 * (x137 * x236 + x145))
    result[3, 6] = numpy.sum(x81 * (x153 + x179 * x199 + x189 * x80))
    result[3, 7] = numpy.sum(x129 * (x158 + x179 * x202 + x193))
    result[3, 8] = numpy.sum(x129 * (x161 + x179 * x205 + x195 * x40))
    result[3, 9] = numpy.sum(x81 * (x164 * x236 * x4 + x167))
    result[3, 10] = numpy.sum(x58 * (x170 + x179 * x209 + 4.0 * x196))
    result[3, 11] = numpy.sum(x81 * (x172 + x179 * x211 + 3.0 * x201))
    result[3, 12] = numpy.sum(x117 * (x174 + x179 * x213 + 2.0 * x203))
    result[3, 13] = numpy.sum(x81 * (x175 + x179 * x215 + x206))
    result[3, 14] = numpy.sum(x58 * (x177 * x236 + x178))
    result[4, 0] = numpy.sum(x216 * x237 * x38)
    result[4, 1] = numpy.sum(x184 * (x181 * x216 + x238))
    result[4, 2] = numpy.sum(x217 * x239)
    result[4, 3] = numpy.sum(x190 * x216 * (x185 + x186))
    result[4, 4] = numpy.sum(x194 * (x179 * x219 + x218))
    result[4, 5] = numpy.sum(x179 * x190 * x221)
    result[4, 6] = numpy.sum(x184 * x216 * (x197 + x198))
    result[4, 7] = numpy.sum(x194 * (x179 * x227 + x220))
    result[4, 8] = numpy.sum(x194 * (x179 * x229 + x222))
    result[4, 9] = numpy.sum(x231 * x239)
    result[4, 10] = numpy.sum(x180 * x216 * (x207 + x208))
    result[4, 11] = numpy.sum(x184 * (x179 * x232 + 3.0 * x224))
    result[4, 12] = numpy.sum(x190 * (x179 * x233 + 2.0 * x228))
    result[4, 13] = numpy.sum(x184 * (x179 * x234 + x230))
    result[4, 14] = numpy.sum(x235 * x237)
    result[5, 0] = numpy.sum(x58 * (x240 * x38 + x54))
    result[5, 1] = numpy.sum(x81 * (x240 * x65 + x76))
    result[5, 2] = numpy.sum(x81 * (x216 * x217 + x238 + x95))
    result[5, 3] = numpy.sum(x117 * (x106 * x240 + x114))
    result[5, 4] = numpy.sum(x129 * (x125 + x216 * x219 + x216 * x71))
    result[5, 5] = numpy.sum(x117 * (x145 + x216 * x221 + 2.0 * x218))
    result[5, 6] = numpy.sum(x81 * (x149 * x240 * x4 + x153))
    result[5, 7] = numpy.sum(x129 * (x158 + x216 * x225 + x216 * x227))
    result[5, 8] = numpy.sum(x129 * (x161 + x216 * x229 + x220))
    result[5, 9] = numpy.sum(x81 * (x167 + x216 * x231 + 3.0 * x222))
    result[5, 10] = numpy.sum(x58 * (x169 * x240 + x170))
    result[5, 11] = numpy.sum(x81 * (x172 + x216 * x232 + x223))
    result[5, 12] = numpy.sum(x117 * (x174 + x216 * x233 + 2.0 * x224))
    result[5, 13] = numpy.sum(x81 * (x175 + x216 * x234 + 3.0 * x228))
    result[5, 14] = numpy.sum(x58 * (x178 + x216 * x235 + 4.0 * x230))
    return result


def _2center2el3d_25(ax, da, A, bx, db, B):
    """Cartesian (d|h) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 21), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = -x2 - B[0]
    x5 = bx ** (-1.0)
    x6 = ax * x1
    x7 = bx * x6 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x8 = boys(6, x7)
    x9 = 17.4934183276249
    x10 = 2.0 * x5
    x11 = x10 * x9
    x12 = x0 ** (-1.5) * x11
    x13 = x12 * x8
    x14 = x13 * x4
    x15 = ax ** (-1.0)
    x16 = x0 ** (-0.5)
    x17 = boys(5, x7)
    x18 = 0.5 * x5
    x19 = x18 * (-x13 + 2.0 * x15 * x16 * x17 * x5 * x9)
    x20 = boys(7, x7)
    x21 = x4**2
    x22 = x11 * x15 * x16
    x23 = x21 * x22
    x24 = x20 * x23
    x25 = x12 * x17
    x26 = boys(4, x7)
    x27 = x18 * (2.0 * x15 * x16 * x26 * x5 * x9 - x25)
    x28 = x23 * x8
    x29 = x27 + x28
    x30 = x12 * x26
    x31 = boys(3, x7)
    x32 = x18 * (2.0 * x15 * x16 * x31 * x5 * x9 - x30)
    x33 = x17 * x22
    x34 = x21 * x33
    x35 = x32 + x34
    x36 = 1.5 * x5
    x37 = x25 * x4
    x38 = x29 * x4 + x5 * (2.0 * x15 * x16 * x26 * x4 * x5 * x9 - x37)
    x39 = x30 * x4
    x40 = x35 * x4 + x5 * (2.0 * x15 * x16 * x31 * x4 * x5 * x9 - x39)
    x41 = -x10 * (x38 * x6 - x40) - x4 * (
        x36 * (x29 * x6 - x35)
        - x4 * (x4 * (x19 + x24) - x5 * (x14 - 2.0 * x15 * x16 * x17 * x4 * x5 * x9))
    )
    x42 = x3 * x41
    x43 = 0.5 / (ax + bx)
    x44 = x12 * x31
    x45 = boys(2, x7)
    x46 = x18 * (2.0 * x15 * x16 * x45 * x5 * x9 - x44)
    x47 = x22 * x26
    x48 = x21 * x47
    x49 = x46 + x48
    x50 = -x36 * (x35 * x6 - x49) + x38 * x4
    x51 = x43 * x50
    x52 = 5.0 * x51
    x53 = x12 * x45
    x54 = boys(1, x7)
    x55 = x18 * (2.0 * x15 * x16 * x5 * x54 * x9 - x53)
    x56 = x22 * x31
    x57 = x21 * x56
    x58 = x55 + x57
    x59 = x4 * (x49 + x5 * (2.0 * x15 * x16 * x45 * x5 * x9 - x44))
    x60 = bx * x1
    x61 = 0.5 * x15
    x62 = x61 * (
        x10 * (x4 * x5 * (2.0 * x15 * x16 * x5 * x54 * x9 - x53) + x4 * x58 - x59 * x6)
        - x4 * (x36 * (x49 * x6 - x58) - x4 * x40)
        + x60 * (x10 * (x40 * x6 - x59) - x4 * x50)
    )
    x63 = x40 * x43
    x64 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**6.5)
    x65 = 6.79952275532455 * x64
    x66 = -x1 * (ax * A[1] + bx * B[1])
    x67 = -x66 - B[1]
    x68 = x4 * x67
    x69 = x13 * x67
    x70 = x5 * (2.0 * x15 * x16 * x17 * x5 * x67 * x9 - x69)
    x71 = x24 * x67
    x72 = x25 * x67
    x73 = x5 * (2.0 * x15 * x16 * x26 * x5 * x67 * x9 - x72)
    x74 = x28 * x67
    x75 = 0.5 * x73 + x74
    x76 = x30 * x67
    x77 = x5 * (2.0 * x15 * x16 * x31 * x5 * x67 * x9 - x76)
    x78 = x34 * x67
    x79 = 0.5 * x77 + x78
    x80 = -x36 * (x6 * x75 - x79) + 0.5 * x4 * (
        x4 * (x70 + 2.0 * x71)
        - 2.0 * x5 * (x13 * x68 - 2.0 * x15 * x16 * x17 * x4 * x5 * x67 * x9)
    )
    x81 = x3 * x80
    x82 = x4 * x75 + x5 * (2.0 * x15 * x16 * x26 * x4 * x5 * x67 * x9 - x25 * x68)
    x83 = x43 * x82
    x84 = 4.0 * x83
    x85 = x5 * x67 * (2.0 * x15 * x16 * x5 * x54 * x9 - x53)
    x86 = x5 * x67 * (2.0 * x15 * x16 * x45 * x5 * x9 - x44)
    x87 = x48 * x67 + 0.5 * x86
    x88 = (
        0.5
        * x61
        * (
            x36 * (2.0 * x57 * x67 - 2.0 * x6 * x87 + x85)
            + 2.0
            * x4
            * (x4 * x79 + x5 * (2.0 * x15 * x16 * x31 * x4 * x5 * x67 * x9 - x30 * x68))
            + 2.0 * x60 * (x36 * (x6 * x79 - x87) - x4 * x82)
        )
    )
    x89 = x43 * x79
    x90 = 4.0 * x43
    x91 = 20.3985682659737 * x64
    x92 = -x1 * (ax * A[2] + bx * B[2])
    x93 = -x92 - B[2]
    x94 = x13 * x93
    x95 = x5 * (2.0 * x15 * x16 * x17 * x5 * x9 * x93 - x94)
    x96 = 0.5 * x95
    x97 = x25 * x93
    x98 = x5 * (2.0 * x15 * x16 * x26 * x5 * x9 * x93 - x97)
    x99 = 0.5 * x98
    x100 = x28 * x93 + x99
    x101 = x30 * x93
    x102 = x5 * (-x101 + 2.0 * x15 * x16 * x31 * x5 * x9 * x93)
    x103 = 0.5 * x102
    x104 = x103 + x34 * x93
    x105 = -x36 * (x100 * x6 - x104) + x4 * (
        x4 * (x24 * x93 + x96) - x5 * x93 * (x14 - 2.0 * x15 * x16 * x17 * x4 * x5 * x9)
    )
    x106 = x105 * x3
    x107 = x100 * x4 + x5 * x93 * (2.0 * x15 * x16 * x26 * x4 * x5 * x9 - x37)
    x108 = x107 * x43
    x109 = 4.0 * x108
    x110 = x5 * x93 * (2.0 * x15 * x16 * x5 * x54 * x9 - x53)
    x111 = 0.5 * x110
    x112 = x5 * x93 * (2.0 * x15 * x16 * x45 * x5 * x9 - x44)
    x113 = 0.5 * x112
    x114 = x113 + x48 * x93
    x115 = x61 * (
        x36 * (x111 - x114 * x6 + x57 * x93)
        + x4 * (x104 * x4 + x5 * x93 * (2.0 * x15 * x16 * x31 * x4 * x5 * x9 - x39))
        - x60 * (x107 * x4 - x36 * (x104 * x6 - x114))
    )
    x116 = x104 * x43
    x117 = x67**2
    x118 = x117 * x33 + x32
    x119 = x117 * x22
    x120 = x119 * x8
    x121 = x120 + x27
    x122 = x121 * x6
    x123 = x119 * x20
    x124 = x123 + x19
    x125 = x118 - x122
    x126 = x4 * (x124 * x21 + x125 * x18 + x5 * (x118 - x122))
    x127 = x126 * x3
    x128 = x117 * x47 + x46
    x129 = -x118 * x6 + x128
    x130 = x121 * x21 + x129 * x18
    x131 = x130 * x43
    x132 = 3.0 * x131
    x133 = x117 * x56 + x55
    x134 = -x128 * x6 + x133
    x135 = x118 * x4
    x136 = -x61 * (
        x4 * x5 * (x128 * x6 - x133)
        - x4 * (x118 * x21 + x134 * x18)
        + x60 * (x130 * x4 + x5 * (x128 * x4 - x135 * x6))
    )
    x137 = 2.0 * x43
    x138 = 3.0 * x43
    x139 = 31.1593277158494 * x64
    x140 = x5 * x93 * (2.0 * x15 * x16 * x17 * x5 * x67 * x9 - x69)
    x141 = 0.5 * x4 * (x140 + 2.0 * x71 * x93) + x5 * (
        2.0 * x15 * x16 * x17 * x4 * x5 * x67 * x9 * x93 - x68 * x94
    )
    x142 = x5 * x93 * (2.0 * x15 * x16 * x26 * x5 * x67 * x9 - x72)
    x143 = 0.5 * x142 + x74 * x93
    x144 = x5 * x93 * (2.0 * x15 * x16 * x31 * x5 * x67 * x9 - x76)
    x145 = (
        -0.5
        * x61
        * (
            -x4 * (x144 + 2.0 * x78 * x93)
            + 2.0 * x5 * (x101 * x68 - 2.0 * x15 * x16 * x31 * x4 * x5 * x67 * x9 * x93)
            + 2.0
            * x60
            * (
                x143 * x4
                + x5 * (2.0 * x15 * x16 * x26 * x4 * x5 * x67 * x9 * x93 - x68 * x97)
            )
        )
    )
    x146 = x15 * x16 * x17 * x5 * x68 * x9 * x90 * x93
    x147 = 53.9695387335403 * x64
    x148 = x93**2
    x149 = x148 * x33 + x32
    x150 = x148 * x22
    x151 = x150 * x8 + x27
    x152 = x151 * x6
    x153 = x150 * x20 + x19
    x154 = x153 * x21
    x155 = x149 - x152
    x156 = x155 * x18
    x157 = x4 * (x154 + x156 + x5 * (x149 - x152))
    x158 = x157 * x3
    x159 = x151 * x21
    x160 = x148 * x47 + x46
    x161 = -x149 * x6 + x160
    x162 = x161 * x18
    x163 = x159 + x162
    x164 = x163 * x43
    x165 = 3.0 * x164
    x166 = x148 * x56 + x55
    x167 = x149 * x21
    x168 = -x160 * x6 + x166
    x169 = x168 * x18
    x170 = x149 * x4
    x171 = -x61 * (
        x4 * x5 * (x160 * x6 - x166)
        - x4 * (x167 + x169)
        + x60 * (x163 * x4 + x5 * (x160 * x4 - x170 * x6))
    )
    x172 = x124 * x67 + x70
    x173 = x121 * x67 + x73
    x174 = x118 * x67 + x77
    x175 = -x173 * x6 + x174
    x176 = x172 * x21 + x175 * x18
    x177 = x176 * x3
    x178 = x173 * x4
    x179 = x128 * x67 + x86
    x180 = x133 * x67 - x179 * x6 + x85
    x181 = -x174 * x6 + x179
    x182 = x61 * (x174 * x21 + x18 * x180 - x60 * (x173 * x21 + x18 * x181))
    x183 = x174 * x43
    x184 = x3 * x4
    x185 = 31.1593277158495 * x64
    x186 = x123 * x93 + x96
    x187 = x120 * x93 + x99
    x188 = x103 + x117 * x33 * x93
    x189 = -x187 * x6 + x188
    x190 = x18 * x189 + x186 * x21
    x191 = x187 * x4
    x192 = x137 * x191
    x193 = x113 + x117 * x47 * x93
    x194 = x111 + x117 * x56 * x93 - x193 * x6
    x195 = -x188 * x6 + x193
    x196 = x61 * (x18 * x194 + x188 * x21 - x60 * (x18 * x195 + x187 * x21))
    x197 = x188 * x43
    x198 = 69.6743749058326 * x64
    x199 = x5 * x67 * (x149 - x152)
    x200 = x154 * x67 + 0.5 * x199
    x201 = x151 * x68
    x202 = x137 * x201
    x203 = x5 * x67 * (-x160 * x6 + x166)
    x204 = x149 * x67
    x205 = x5 * (x160 * x67 - x204 * x6)
    x206 = 0.5 * x61 * (2.0 * x167 * x67 + x203 - x60 * (2.0 * x159 * x67 + x205))
    x207 = x184 * x67
    x208 = x153 * x93 + x95
    x209 = x151 * x93 + x98
    x210 = x102 + x149 * x93
    x211 = -x209 * x6 + x210
    x212 = x18 * x211
    x213 = x208 * x21 + x212
    x214 = x213 * x3
    x215 = x209 * x4
    x216 = x112 + x160 * x93
    x217 = x110 + x166 * x93 - x216 * x6
    x218 = x18 * x217
    x219 = -x210 * x6 + x216
    x220 = x18 * x219
    x221 = x61 * (x21 * x210 + x218 - x60 * (x209 * x21 + x220))
    x222 = x210 * x43
    x223 = x129 * x36 + x173 * x67
    x224 = x223 * x43
    x225 = x125 * x36 + x172 * x67
    x226 = x184 * x225
    x227 = x134 * x36 + x174 * x67
    x228 = x4 * x60
    x229 = x61 * (-x223 * x228 + x227 * x4)
    x230 = x142 + x187 * x67
    x231 = x230 * x43
    x232 = x140 + x186 * x67
    x233 = x144 + x188 * x67
    x234 = x61 * (-x228 * x230 + x233 * x4)
    x235 = x117 * x151 + x162
    x236 = x235 * x43
    x237 = x117 * x153 + x156
    x238 = x117 * x149 + x169
    x239 = x61 * (-x228 * x235 + x238 * x4)
    x240 = x209 * x67
    x241 = x240 * x43
    x242 = x61 * (x210 * x4 * x67 - x228 * x240)
    x243 = x161 * x36 + x209 * x93
    x244 = x243 * x43
    x245 = x155 * x36 + x208 * x93
    x246 = x184 * x245
    x247 = x168 * x36 + x210 * x93
    x248 = x61 * (-x228 * x243 + x247 * x4)
    x249 = x3**2
    x250 = x10 * x175 + x225 * x67
    x251 = x61 * (x10 * x180 + x227 * x67 - x60 * (x10 * x181 + x223 * x67))
    x252 = x189 * x36 + x232 * x67
    x253 = x61 * (x194 * x36 + x233 * x67 - x60 * (x195 * x36 + x230 * x67))
    x254 = x199 + x237 * x67
    x255 = x61 * (x203 + x238 * x67 - x60 * (x205 + x235 * x67))
    x256 = x117 * x208 + x212
    x257 = x61 * (x117 * x210 + x218 - x60 * (x117 * x209 + x220))
    x258 = x61 * x67 * (-x243 * x60 + x247)
    x259 = x245 * x67
    x260 = x10 * x211 + x245 * x93
    x261 = x61 * (x10 * x217 + x247 * x93 - x60 * (x10 * x219 + x243 * x93))
    x262 = -x66 - A[1]
    x263 = 11.7771188794428 * x64
    x264 = x262 * x80
    x265 = x264 + x51
    x266 = x262 * x82 + x63
    x267 = 35.3313566383285 * x64
    x268 = x126 * x262
    x269 = 2.0 * x83
    x270 = x268 + x269
    x271 = x138 * (x130 * x262 + 2.0 * x89)
    x272 = x108 + x141 * x262
    x273 = x116 + x143 * x262
    x274 = 93.4779831475484 * x64
    x275 = x176 * x262
    x276 = x132 + x275
    x277 = x262 * x4
    x278 = x135 * x138 + x173 * x277
    x279 = 53.9695387335404 * x64
    x280 = x137 * x143
    x281 = x190 * x262 + x280
    x282 = x146 + x187 * x277
    x283 = 120.679557322504 * x64
    x284 = x164 + x200 * x262
    x285 = x262 * x68
    x286 = x137 * (x151 * x285 + x170 * x43)
    x287 = x209 * x277
    x288 = x43 * (4.0 * x183 + x223 * x262)
    x289 = x225 * x277
    x290 = x178 * x90 + x289
    x291 = x43 * (3.0 * x197 + x230 * x262)
    x292 = x138 * x191 + x232 * x277
    x293 = x43 * (x137 * x204 + x235 * x262)
    x294 = x202 + x237 * x277
    x295 = x43 * (x222 + x240 * x262)
    x296 = x208 * x285 + x215 * x43
    x297 = x244 * x262
    x298 = x250 * x262
    x299 = 5.0 * x224
    x300 = x298 + x299
    x301 = x263 * x3
    x302 = 4.0 * x231 + x252 * x262
    x303 = x267 * x3
    x304 = 3.0 * x236
    x305 = x254 * x262 + x304
    x306 = x147 * x3
    x307 = x137 * x240 + x256 * x262
    x308 = x279 * x3
    x309 = x244 + x259 * x262
    x310 = -x92 - A[2]
    x311 = x105 * x310 + x51
    x312 = x43 * (x107 * x310 + x63)
    x313 = x132 * x310
    x314 = x141 * x310 + x83
    x315 = x143 * x310 + x89
    x316 = 2.0 * x108 + x157 * x310
    x317 = x43 * (2.0 * x116 + x163 * x310)
    x318 = 3.0 * x317
    x319 = x310 * x4
    x320 = x173 * x319
    x321 = x131 + x190 * x310
    x322 = x135 * x43 + x187 * x319
    x323 = x137 * x322
    x324 = x200 * x310 + x280
    x325 = x310 * x68
    x326 = x146 + x151 * x325
    x327 = x137 * x326
    x328 = x165 + x213 * x310
    x329 = x43 * (x138 * x170 + x215 * x310)
    x330 = x224 * x310
    x331 = x43 * (x183 + x230 * x310)
    x332 = x178 * x43 + x232 * x319
    x333 = x43 * (2.0 * x197 + x235 * x310)
    x334 = x192 + x237 * x319
    x335 = x43 * (x138 * x204 + x240 * x310)
    x336 = x138 * x201 + x208 * x325
    x337 = x43 * (4.0 * x222 + x243 * x310)
    x338 = x215 * x90 + x245 * x319
    x339 = x224 + x252 * x310
    x340 = 2.0 * x231 + x254 * x310
    x341 = x256 * x310 + x304
    x342 = x240 * x90 + x259 * x310
    x343 = 5.0 * x244 + x260 * x310
    x344 = x262**2
    x345 = x262 * x263
    x346 = x310 * x51
    x347 = x262 * x267
    x348 = x137 * x315
    x349 = 3.0 * x333
    x350 = x310**2

    # 126 item(s)
    result[0, 0] = numpy.sum(
        x65 * (x3 * (x42 + x52) + 5.0 * x43 * (x3 * x50 + 4.0 * x63) + x62)
    )
    result[0, 1] = numpy.sum(
        x91 * (x3 * (x81 + x84) + x88 + x90 * (x3 * x82 + 3.0 * x89))
    )
    result[0, 2] = numpy.sum(
        x91 * (x115 + x3 * (x106 + x109) + x90 * (x107 * x3 + 3.0 * x116))
    )
    result[0, 3] = numpy.sum(
        x139 * (x136 + x138 * (x130 * x3 + x135 * x137) + x3 * (x127 + x132))
    )
    result[0, 4] = numpy.sum(
        x147 * (x138 * (x143 * x3 + x146) + x145 + x3 * (x138 * x143 + x141 * x3))
    )
    result[0, 5] = numpy.sum(
        x139 * (x138 * (x137 * x170 + x163 * x3) + x171 + x3 * (x158 + x165))
    )
    result[0, 6] = numpy.sum(
        x185 * (x137 * (x173 * x184 + x183) + x182 + x3 * (x137 * x178 + x177))
    )
    result[0, 7] = numpy.sum(
        x198 * (x137 * (x184 * x187 + x197) + x196 + x3 * (x190 * x3 + x192))
    )
    result[0, 8] = numpy.sum(
        x198 * (x137 * (x151 * x207 + x204 * x43) + x206 + x3 * (x200 * x3 + x202))
    )
    result[0, 9] = numpy.sum(
        x185 * (x137 * (x184 * x209 + x222) + x221 + x3 * (x137 * x215 + x214))
    )
    result[0, 10] = numpy.sum(x91 * (x224 * x3 + x229 + x3 * (x224 + x226)))
    result[0, 11] = numpy.sum(x147 * (x231 * x3 + x234 + x3 * (x184 * x232 + x231)))
    result[0, 12] = numpy.sum(x198 * (x236 * x3 + x239 + x3 * (x184 * x237 + x236)))
    result[0, 13] = numpy.sum(x147 * (x241 * x3 + x242 + x3 * (x207 * x208 + x241)))
    result[0, 14] = numpy.sum(x91 * (x244 * x3 + x248 + x3 * (x244 + x246)))
    result[0, 15] = numpy.sum(x65 * (x249 * x250 + x251))
    result[0, 16] = numpy.sum(x91 * (x249 * x252 + x253))
    result[0, 17] = numpy.sum(x139 * (x249 * x254 + x255))
    result[0, 18] = numpy.sum(x185 * (x249 * x256 + x257))
    result[0, 19] = numpy.sum(x91 * (x249 * x259 + x258))
    result[0, 20] = numpy.sum(x65 * (x249 * x260 + x261))
    result[1, 0] = numpy.sum(x262 * x263 * (x42 + x52))
    result[1, 1] = numpy.sum(x267 * (x265 * x3 + x266 * x90))
    result[1, 2] = numpy.sum(x262 * x267 * (x106 + x109))
    result[1, 3] = numpy.sum(x147 * (x270 * x3 + x271))
    result[1, 4] = numpy.sum(x274 * (x138 * x273 + x272 * x3))
    result[1, 5] = numpy.sum(x147 * x262 * (x158 + x165))
    result[1, 6] = numpy.sum(x279 * (x137 * x278 + x276 * x3))
    result[1, 7] = numpy.sum(x283 * (x137 * x282 + x281 * x3))
    result[1, 8] = numpy.sum(x283 * (x284 * x3 + x286))
    result[1, 9] = numpy.sum(x279 * (x137 * x287 + x214 * x262))
    result[1, 10] = numpy.sum(x267 * (x288 + x290 * x3))
    result[1, 11] = numpy.sum(x274 * (x291 + x292 * x3))
    result[1, 12] = numpy.sum(x283 * (x293 + x294 * x3))
    result[1, 13] = numpy.sum(x274 * (x295 + x296 * x3))
    result[1, 14] = numpy.sum(x267 * (x246 * x262 + x297))
    result[1, 15] = numpy.sum(x300 * x301)
    result[1, 16] = numpy.sum(x302 * x303)
    result[1, 17] = numpy.sum(x305 * x306)
    result[1, 18] = numpy.sum(x307 * x308)
    result[1, 19] = numpy.sum(x303 * x309)
    result[1, 20] = numpy.sum(x260 * x262 * x301)
    result[2, 0] = numpy.sum(x263 * x310 * (x42 + x52))
    result[2, 1] = numpy.sum(x267 * x310 * (x81 + x84))
    result[2, 2] = numpy.sum(x267 * (x3 * x311 + 4.0 * x312))
    result[2, 3] = numpy.sum(x147 * (x127 * x310 + x313))
    result[2, 4] = numpy.sum(x274 * (x138 * x315 + x3 * x314))
    result[2, 5] = numpy.sum(x147 * (x3 * x316 + x318))
    result[2, 6] = numpy.sum(x279 * (x137 * x320 + x177 * x310))
    result[2, 7] = numpy.sum(x283 * (x3 * x321 + x323))
    result[2, 8] = numpy.sum(x283 * (x3 * x324 + x327))
    result[2, 9] = numpy.sum(x279 * (x3 * x328 + 2.0 * x329))
    result[2, 10] = numpy.sum(x267 * (x226 * x310 + x330))
    result[2, 11] = numpy.sum(x274 * (x3 * x332 + x331))
    result[2, 12] = numpy.sum(x283 * (x3 * x334 + x333))
    result[2, 13] = numpy.sum(x274 * (x3 * x336 + x335))
    result[2, 14] = numpy.sum(x267 * (x3 * x338 + x337))
    result[2, 15] = numpy.sum(x250 * x301 * x310)
    result[2, 16] = numpy.sum(x303 * x339)
    result[2, 17] = numpy.sum(x306 * x340)
    result[2, 18] = numpy.sum(x308 * x341)
    result[2, 19] = numpy.sum(x303 * x342)
    result[2, 20] = numpy.sum(x301 * x343)
    result[3, 0] = numpy.sum(x65 * (x344 * x41 + x62))
    result[3, 1] = numpy.sum(x91 * (x262 * x265 + x262 * x51 + x88))
    result[3, 2] = numpy.sum(x91 * (x105 * x344 + x115))
    result[3, 3] = numpy.sum(x139 * (x136 + x137 * x266 + x262 * x270))
    result[3, 4] = numpy.sum(x147 * (x108 * x262 + x145 + x262 * x272))
    result[3, 5] = numpy.sum(x139 * (x157 * x344 + x171))
    result[3, 6] = numpy.sum(x185 * (x182 + x262 * x276 + x271))
    result[3, 7] = numpy.sum(x198 * (x137 * x273 + x196 + x262 * x281))
    result[3, 8] = numpy.sum(x198 * (x164 * x262 + x206 + x262 * x284))
    result[3, 9] = numpy.sum(x185 * (x213 * x344 + x221))
    result[3, 10] = numpy.sum(x91 * (x229 + x262 * x290 + x278 * x90))
    result[3, 11] = numpy.sum(x147 * (x138 * x282 + x234 + x262 * x292))
    result[3, 12] = numpy.sum(x198 * (x239 + x262 * x294 + x286))
    result[3, 13] = numpy.sum(x147 * (x242 + x262 * x296 + x287 * x43))
    result[3, 14] = numpy.sum(x91 * (x245 * x344 * x4 + x248))
    result[3, 15] = numpy.sum(x65 * (x251 + x262 * x300 + 5.0 * x288))
    result[3, 16] = numpy.sum(x91 * (x253 + x262 * x302 + 4.0 * x291))
    result[3, 17] = numpy.sum(x139 * (x255 + x262 * x305 + 3.0 * x293))
    result[3, 18] = numpy.sum(x185 * (x257 + x262 * x307 + 2.0 * x295))
    result[3, 19] = numpy.sum(x91 * (x258 + x262 * x309 + x297))
    result[3, 20] = numpy.sum(x65 * (x260 * x344 + x261))
    result[4, 0] = numpy.sum(x310 * x345 * x41)
    result[4, 1] = numpy.sum(x267 * (x264 * x310 + x346))
    result[4, 2] = numpy.sum(x311 * x347)
    result[4, 3] = numpy.sum(x147 * x310 * (x268 + x269))
    result[4, 4] = numpy.sum(x274 * (x262 * x314 + x312))
    result[4, 5] = numpy.sum(x147 * x262 * x316)
    result[4, 6] = numpy.sum(x279 * (x275 * x310 + x313))
    result[4, 7] = numpy.sum(x283 * (x262 * x321 + x348))
    result[4, 8] = numpy.sum(x283 * (x262 * x324 + x317))
    result[4, 9] = numpy.sum(x262 * x279 * x328)
    result[4, 10] = numpy.sum(x267 * (x289 * x310 + x320 * x90))
    result[4, 11] = numpy.sum(x274 * (x138 * x322 + x262 * x332))
    result[4, 12] = numpy.sum(x283 * (x262 * x334 + x327))
    result[4, 13] = numpy.sum(x274 * (x262 * x336 + x329))
    result[4, 14] = numpy.sum(x338 * x347)
    result[4, 15] = numpy.sum(x263 * x310 * (x298 + x299))
    result[4, 16] = numpy.sum(x267 * (x262 * x339 + 4.0 * x331))
    result[4, 17] = numpy.sum(x147 * (x262 * x340 + x349))
    result[4, 18] = numpy.sum(x279 * (x262 * x341 + 2.0 * x335))
    result[4, 19] = numpy.sum(x267 * (x262 * x342 + x337))
    result[4, 20] = numpy.sum(x343 * x345)
    result[5, 0] = numpy.sum(x65 * (x350 * x41 + x62))
    result[5, 1] = numpy.sum(x91 * (x350 * x80 + x88))
    result[5, 2] = numpy.sum(x91 * (x115 + x310 * x311 + x346))
    result[5, 3] = numpy.sum(x139 * (x126 * x350 + x136))
    result[5, 4] = numpy.sum(x147 * (x145 + x310 * x314 + x310 * x83))
    result[5, 5] = numpy.sum(x139 * (x171 + x310 * x316 + 2.0 * x312))
    result[5, 6] = numpy.sum(x185 * (x176 * x350 + x182))
    result[5, 7] = numpy.sum(x198 * (x131 * x310 + x196 + x310 * x321))
    result[5, 8] = numpy.sum(x198 * (x206 + x310 * x324 + x348))
    result[5, 9] = numpy.sum(x185 * (x221 + x310 * x328 + x318))
    result[5, 10] = numpy.sum(x91 * (x225 * x350 * x4 + x229))
    result[5, 11] = numpy.sum(x147 * (x234 + x310 * x332 + x320 * x43))
    result[5, 12] = numpy.sum(x198 * (x239 + x310 * x334 + x323))
    result[5, 13] = numpy.sum(x147 * (x138 * x326 + x242 + x310 * x336))
    result[5, 14] = numpy.sum(x91 * (x248 + x310 * x338 + 4.0 * x329))
    result[5, 15] = numpy.sum(x65 * (x250 * x350 + x251))
    result[5, 16] = numpy.sum(x91 * (x253 + x310 * x339 + x330))
    result[5, 17] = numpy.sum(x139 * (x255 + x310 * x340 + 2.0 * x331))
    result[5, 18] = numpy.sum(x185 * (x257 + x310 * x341 + x349))
    result[5, 19] = numpy.sum(x91 * (x258 + x310 * x342 + 4.0 * x335))
    result[5, 20] = numpy.sum(x65 * (x261 + x310 * x343 + 5.0 * x337))
    return result


def _2center2el3d_30(ax, da, A, bx, db, B):
    """Cartesian (f|s) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 1), dtype=float)

    x0 = ax ** (-1.0)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = x2 * (ax * A[0] + bx * B[0]) - A[0]
    x4 = ax * bx * x2 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x5 = 17.4934183276249
    x6 = 2.0 * x0 * x5
    x7 = x1 ** (-1.5) * x6 * boys(2, x4)
    x8 = bx ** (-1.0)
    x9 = x1 ** (-0.5)
    x10 = boys(1, x4)
    x11 = 0.5 * x0 * (2.0 * x0 * x10 * x5 * x8 * x9 - x7)
    x12 = boys(3, x4)
    x13 = x8 * x9
    x14 = x12 * x13 * x6
    x15 = x14 * x3**2
    x16 = da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**1.5)
    x17 = 0.179587122125167 * x16
    x18 = 5.84237394672177 * x17
    x19 = x2 * (ax * A[1] + bx * B[1]) - A[1]
    x20 = x0 * x19 * (2.0 * x0 * x10 * x5 * x8 * x9 - x7)
    x21 = 13.0639452948436 * x17
    x22 = x2 * (ax * A[2] + bx * B[2]) - A[2]
    x23 = x0 * x22 * (2.0 * x0 * x10 * x5 * x8 * x9 - x7)
    x24 = 0.5 * x23
    x25 = x14 * x19**2
    x26 = x11 + x25
    x27 = x21 * x3
    x28 = x11 + x14 * x22**2

    # 10 item(s)
    result[0, 0] = numpy.sum(
        x18 * x3 * (x0 * (2.0 * x0 * x10 * x5 * x8 * x9 - x7) + x11 + x15)
    )
    result[1, 0] = numpy.sum(0.5 * x21 * (2.0 * x15 * x19 + x20))
    result[2, 0] = numpy.sum(x21 * (x15 * x22 + x24))
    result[3, 0] = numpy.sum(x26 * x27)
    result[4, 0] = numpy.sum(142.172254021068 * x0 * x12 * x13 * x16 * x19 * x22 * x3)
    result[5, 0] = numpy.sum(x27 * x28)
    result[6, 0] = numpy.sum(x18 * (x19 * x26 + x20))
    result[7, 0] = numpy.sum(x21 * (x22 * x25 + x24))
    result[8, 0] = numpy.sum(x19 * x21 * x28)
    result[9, 0] = numpy.sum(x18 * (x22 * x28 + x23))
    return result


def _2center2el3d_31(ax, da, A, bx, db, B):
    """Cartesian (f|p) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 3), dtype=float)

    x0 = ax ** (-1.0)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = bx * x2
    x4 = ax * x3 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x5 = boys(2, x4)
    x6 = 0.5 / (ax + bx)
    x7 = bx ** (-1.0)
    x8 = x1 ** (-0.5)
    x9 = 17.4934183276249
    x10 = 2.0 * x0 * x9
    x11 = x10 * x7 * x8
    x12 = x11 * x6
    x13 = x12 * x5
    x14 = boys(3, x4)
    x15 = x11 * x14
    x16 = -x2 * (ax * A[0] + bx * B[0])
    x17 = -x16 - A[0]
    x18 = -x16 - B[0]
    x19 = x17 * x18
    x20 = boys(1, x4)
    x21 = x11 * x20 * x6
    x22 = x11 * x5
    x23 = x1 ** (-1.5) * x10
    x24 = 0.5 * x0
    x25 = x24 * (2.0 * x0 * x20 * x7 * x8 * x9 - x23 * x5)
    x26 = x17**2
    x27 = x12 * x14
    x28 = x11 * boys(4, x4)
    x29 = x19 * x28
    x30 = x14 * x23
    x31 = x18 * x30
    x32 = x24 * (2.0 * x0 * x18 * x5 * x7 * x8 * x9 - x31)
    x33 = x17 * x27
    x34 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**2.5)
    x35 = 11.6847478934435 * x34
    x36 = -x2 * (ax * A[1] + bx * B[1])
    x37 = -x36 - B[1]
    x38 = x30 * x37
    x39 = x24 * (2.0 * x0 * x37 * x5 * x7 * x8 * x9 - x38)
    x40 = x26 * x28
    x41 = x37 * x40
    x42 = -x2 * (ax * A[2] + bx * B[2])
    x43 = -x42 - B[2]
    x44 = x30 * x43
    x45 = x24 * (2.0 * x0 * x43 * x5 * x7 * x8 * x9 - x44)
    x46 = x40 * x43
    x47 = -x36 - A[1]
    x48 = x27 * x47
    x49 = x29 * x47
    x50 = x0 * x47 * (2.0 * x0 * x18 * x5 * x7 * x8 * x9 - x31)
    x51 = 26.1278905896872 * x34
    x52 = x37 * x47
    x53 = x28 * x52
    x54 = x27 + x53
    x55 = x0 * (x21 + x22 * x37 * x47 - x3 * (x13 + x15 * x52))
    x56 = x0 * x47 * (2.0 * x0 * x43 * x5 * x7 * x8 * x9 - x44)
    x57 = -x42 - A[2]
    x58 = x27 * x57
    x59 = x0 * x57 * (2.0 * x0 * x18 * x5 * x7 * x8 * x9 - x31)
    x60 = 0.5 * x59
    x61 = x0 * x57 * (2.0 * x0 * x37 * x5 * x7 * x8 * x9 - x38)
    x62 = 0.5 * x61
    x63 = x43 * x57
    x64 = x27 + x28 * x63
    x65 = x0 * (x21 + x22 * x43 * x57 - x3 * (x13 + x15 * x63))
    x66 = 0.5 * x65
    x67 = x47**2
    x68 = x6 * (x15 * x67 + x25)
    x69 = x28 * x67
    x70 = x18 * x69
    x71 = x32 + x70
    x72 = x39 + x47 * x54 + x48
    x73 = x17 * x51
    x74 = x43 * x69 + x45
    x75 = x48 * x57
    x76 = 45.2548339959391 * x34
    x77 = x53 * x57 + x58
    x78 = x17 * x76
    x79 = x57**2
    x80 = x6 * (x15 * x79 + x25)
    x81 = x28 * x79
    x82 = x18 * x81 + x32
    x83 = x37 * x81 + x39
    x84 = x45 + x57 * x64 + x58
    x85 = x47 * x51

    # 30 item(s)
    result[0, 0] = numpy.sum(
        x35
        * (
            x0 * (x17 * x18 * x22 + x21 - x3 * (x13 + x15 * x19))
            + x17 * (x17 * (x27 + x29) + x32 + x33)
            + x6 * (x15 * x26 + x25)
        )
    )
    result[0, 1] = numpy.sum(
        x17 * x35 * (x0 * (2.0 * x0 * x37 * x5 * x7 * x8 * x9 - x38) + x39 + x41)
    )
    result[0, 2] = numpy.sum(
        x17 * x35 * (x0 * (2.0 * x0 * x43 * x5 * x7 * x8 * x9 - x44) + x45 + x46)
    )
    result[1, 0] = numpy.sum(
        0.5 * x51 * (2.0 * x17 * (x48 + x49) + 2.0 * x33 * x47 + x50)
    )
    result[1, 1] = numpy.sum(0.5 * x51 * (2.0 * x26 * x54 + x55))
    result[1, 2] = numpy.sum(0.5 * x51 * (2.0 * x46 * x47 + x56))
    result[2, 0] = numpy.sum(x51 * (x17 * (x29 * x57 + x58) + x33 * x57 + x60))
    result[2, 1] = numpy.sum(x51 * (x41 * x57 + x62))
    result[2, 2] = numpy.sum(x51 * (x26 * x64 + x66))
    result[3, 0] = numpy.sum(x51 * (x17 * x71 + x68))
    result[3, 1] = numpy.sum(x72 * x73)
    result[3, 2] = numpy.sum(x73 * x74)
    result[4, 0] = numpy.sum(x76 * (x49 * x57 + x75))
    result[4, 1] = numpy.sum(x77 * x78)
    result[4, 2] = numpy.sum(x47 * x64 * x78)
    result[5, 0] = numpy.sum(x51 * (x17 * x82 + x80))
    result[5, 1] = numpy.sum(x73 * x83)
    result[5, 2] = numpy.sum(x73 * x84)
    result[6, 0] = numpy.sum(x35 * (x47 * x71 + x50))
    result[6, 1] = numpy.sum(x35 * (x47 * x72 + x55 + x68))
    result[6, 2] = numpy.sum(x35 * (x47 * x74 + x56))
    result[7, 0] = numpy.sum(x51 * (x57 * x70 + x60))
    result[7, 1] = numpy.sum(x51 * (x47 * x77 + x62 + x75))
    result[7, 2] = numpy.sum(x51 * (x64 * x67 + x66))
    result[8, 0] = numpy.sum(x82 * x85)
    result[8, 1] = numpy.sum(x51 * (x47 * x83 + x80))
    result[8, 2] = numpy.sum(x84 * x85)
    result[9, 0] = numpy.sum(x35 * (x57 * x82 + x59))
    result[9, 1] = numpy.sum(x35 * (x57 * x83 + x61))
    result[9, 2] = numpy.sum(x35 * (x57 * x84 + x65 + x80))
    return result


def _2center2el3d_32(ax, da, A, bx, db, B):
    """Cartesian (f|d) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 6), dtype=float)

    x0 = ax ** (-1.0)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = -x2 * (ax * A[0] + bx * B[0])
    x4 = -x3 - A[0]
    x5 = bx * x2
    x6 = ax * x5 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x7 = boys(2, x6)
    x8 = x1 ** (-1.5)
    x9 = bx ** (-1.0)
    x10 = 17.4934183276249
    x11 = x10 * x9
    x12 = 2.0 * x11
    x13 = x12 * x8
    x14 = x1 ** (-0.5)
    x15 = 0.5 * x9
    x16 = x15 * (2.0 * x0 * x10 * x14 * x9 * boys(1, x6) - x13 * x7)
    x17 = -x3 - B[0]
    x18 = x17**2
    x19 = boys(3, x6)
    x20 = x0 * x14
    x21 = x12 * x20
    x22 = x19 * x21
    x23 = x16 + x18 * x22
    x24 = x15 * (2.0 * x0 * x10 * x14 * x7 * x9 - x13 * x19)
    x25 = boys(4, x6)
    x26 = x18 * x21
    x27 = x24 + x25 * x26
    x28 = x17 * x19
    x29 = 0.5 / (ax + bx)
    x30 = x11 * x20
    x31 = 4.0 * x29 * x30
    x32 = x17 * x7
    x33 = x15 * (2.0 * x0 * x10 * x14 * x19 * x9 - x13 * x25)
    x34 = boys(5, x6)
    x35 = x26 * x34 + x33
    x36 = x35 * x4
    x37 = x17 * x25
    x38 = x31 * x37
    x39 = x27 * x5
    x40 = 0.5 * x0
    x41 = x40 * (x23 - x39)
    x42 = 2.0 * x29
    x43 = x30 * x42
    x44 = x19 * x43
    x45 = x21 * x37
    x46 = x4 * x45
    x47 = x44 + x46
    x48 = 2.0 * x0 * x10 * x8
    x49 = x40 * (2.0 * x0 * x10 * x14 * x17 * x7 * x9 - x28 * x48)
    x50 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**3.5)
    x51 = 13.4923846833851 * x50
    x52 = -x2 * (ax * A[1] + bx * B[1])
    x53 = -x52 - B[1]
    x54 = x44 * x53
    x55 = x43 * x7
    x56 = x53 * x55
    x57 = x17 * x22
    x58 = x53 * x57
    x59 = x19 * x48
    x60 = x40 * x53 * (2.0 * x0 * x10 * x14 * x7 * x9 - x59)
    x61 = x25 * x53
    x62 = x4**2
    x63 = x21 * x62
    x64 = x43 * x61
    x65 = x21 * x34
    x66 = x17 * x65
    x67 = x4 * x66
    x68 = x53 * x67
    x69 = x37 * x48
    x70 = x53 * x69
    x71 = x40 * (2.0 * x0 * x10 * x14 * x17 * x19 * x53 * x9 - x70)
    x72 = x4 * x64
    x73 = 23.3694957868871 * x50
    x74 = -x2 * (ax * A[2] + bx * B[2])
    x75 = -x74 - B[2]
    x76 = x44 * x75
    x77 = x55 * x75
    x78 = x57 * x75
    x79 = x40 * x75 * (2.0 * x0 * x10 * x14 * x7 * x9 - x59)
    x80 = x25 * x75
    x81 = x43 * x80
    x82 = x67 * x75
    x83 = x69 * x75
    x84 = x40 * (2.0 * x0 * x10 * x14 * x17 * x19 * x75 * x9 - x83)
    x85 = x4 * x81
    x86 = x53**2
    x87 = x16 + x22 * x86
    x88 = x21 * x86
    x89 = x24 + x25 * x88
    x90 = x5 * x89
    x91 = x33 + x34 * x88
    x92 = x62 * x91
    x93 = x40 * (x87 - x90)
    x94 = x61 * x75
    x95 = x48 * x94
    x96 = x40 * (2.0 * x0 * x10 * x14 * x19 * x53 * x75 * x9 - x95)
    x97 = x53 * x75
    x98 = x75**2
    x99 = x16 + x22 * x98
    x100 = x21 * x98
    x101 = x100 * x25 + x24
    x102 = x101 * x5
    x103 = x100 * x34 + x33
    x104 = x103 * x62
    x105 = x40 * (-x102 + x99)
    x106 = -x52 - A[1]
    x107 = x106 * x36
    x108 = x106 * x38
    x109 = x0 * x106 * (x23 - x39)
    x110 = x106 * x44
    x111 = 30.169889330626 * x50
    x112 = x21 * x61
    x113 = x106 * x112
    x114 = x113 + x44
    x115 = x114 * x29
    x116 = x37 * x43
    x117 = x106 * x53
    x118 = x117 * x66
    x119 = x116 + x118
    x120 = x28 * x43
    x121 = x32 * x43
    x122 = x0 * (x106 * x58 + x121 - x5 * (x117 * x45 + x120))
    x123 = 52.2557811793745 * x50
    x124 = x106 * x81
    x125 = x0 * x106 * (2.0 * x0 * x10 * x14 * x17 * x19 * x75 * x9 - x83)
    x126 = x106 * x91
    x127 = x31 * x61
    x128 = x126 + x127
    x129 = x31 * x53
    x130 = x0 * (x106 * x87 + x129 * x7 - x5 * (x106 * x89 + x129 * x19))
    x131 = x65 * x97
    x132 = x106 * x131 + x81
    x133 = x21 * x94
    x134 = x22 * x97
    x135 = x0 * (x106 * x134 - x5 * (x106 * x133 + x76) + x77)
    x136 = x0 * x106 * (-x102 + x99)
    x137 = -x74 - A[2]
    x138 = x0 * x137 * (x23 - x39)
    x139 = 0.5 * x138
    x140 = x137 * x44
    x141 = x137 * x64
    x142 = x0 * x137 * (2.0 * x0 * x10 * x14 * x17 * x19 * x53 * x9 - x70)
    x143 = 0.5 * x142
    x144 = x21 * x80
    x145 = x137 * x144 + x44
    x146 = x145 * x29
    x147 = x137 * x75
    x148 = x116 + x147 * x66
    x149 = x148 * x4
    x150 = x0 * (x121 + x137 * x78 - x5 * (x120 + x147 * x45))
    x151 = 0.5 * x150
    x152 = x0 * x137 * (x87 - x90)
    x153 = 0.5 * x152
    x154 = x131 * x137 + x64
    x155 = x0 * (x134 * x137 - x5 * (x133 * x137 + x54) + x56)
    x156 = 0.5 * x155
    x157 = x103 * x137 + x31 * x80
    x158 = x31 * x75
    x159 = x0 * (x137 * x99 + x158 * x7 - x5 * (x101 * x137 + x158 * x19))
    x160 = 0.5 * x159
    x161 = x106**2
    x162 = x161 * x35
    x163 = x162 + x41
    x164 = x29 * (x161 * x45 + x49)
    x165 = x29 * (x106 * x114 + x110 + x60)
    x166 = x106 * x116
    x167 = x106 * x119 + x166 + x71
    x168 = x29 * (x144 * x161 + x79)
    x169 = x161 * x66 * x75 + x84
    x170 = x106 * x128 + 2.0 * x115 + x93
    x171 = x111 * x4
    x172 = x106 * x132 + x124 + x96
    x173 = x123 * x4
    x174 = x103 * x161 + x105
    x175 = x29 * (x113 * x137 + x140)
    x176 = x116 * x137
    x177 = x118 * x137 + x176
    x178 = 90.5096679918781 * x50
    x179 = x106 * x146
    x180 = x137 * (x126 + x127)
    x181 = x106 * x154 + x146
    x182 = x137**2
    x183 = x182 * x35 + x41
    x184 = x29 * (x182 * x45 + x49)
    x185 = x29 * (x112 * x182 + x60)
    x186 = x182 * x53 * x66 + x71
    x187 = x29 * (x137 * x145 + x140 + x79)
    x188 = x137 * x148 + x176 + x84
    x189 = x182 * x91 + x93
    x190 = x137 * x154 + x141 + x96
    x191 = x105 + x137 * x157 + 2.0 * x146
    x192 = x106 * x111

    # 60 item(s)
    result[0, 0] = numpy.sum(
        x51
        * (
            x0 * (x23 * x4 + x31 * x32 - x5 * (x27 * x4 + x28 * x31))
            + x4 * (x4 * (x36 + x38) + x41 + x42 * x47)
            + x42 * (x4 * x44 + x4 * x47 + x49)
        )
    )
    result[0, 1] = numpy.sum(
        x73
        * (
            x0 * (x4 * x58 - x5 * (x46 * x53 + x54) + x56)
            + x29 * (x60 + x61 * x63)
            + x4 * (x4 * (x64 + x68) + x71 + x72)
        )
    )
    result[0, 2] = numpy.sum(
        x73
        * (
            x0 * (x4 * x78 - x5 * (x46 * x75 + x76) + x77)
            + x29 * (x63 * x80 + x79)
            + x4 * (x4 * (x81 + x82) + x84 + x85)
        )
    )
    result[0, 3] = numpy.sum(x4 * x51 * (x0 * (x87 - x90) + x92 + x93))
    result[0, 4] = numpy.sum(
        x4
        * x73
        * (
            x0 * (2.0 * x0 * x10 * x14 * x19 * x53 * x75 * x9 - x95)
            + x34 * x63 * x97
            + x96
        )
    )
    result[0, 5] = numpy.sum(x4 * x51 * (-x0 * (x102 - x99) + x104 + x105))
    result[1, 0] = numpy.sum(
        0.5 * x111 * (x109 + 2.0 * x4 * (x107 + x108) + 2.0 * x42 * (x106 * x46 + x110))
    )
    result[1, 1] = numpy.sum(0.5 * x123 * (4.0 * x115 * x4 + 2.0 * x119 * x4**2 + x122))
    result[1, 2] = numpy.sum(
        0.5 * x123 * (2.0 * x106 * x85 + x125 + 2.0 * x4 * (x106 * x82 + x124))
    )
    result[1, 3] = numpy.sum(0.5 * x111 * (2.0 * x128 * x62 + x130))
    result[1, 4] = numpy.sum(0.5 * x123 * (2.0 * x132 * x62 + x135))
    result[1, 5] = numpy.sum(0.5 * x111 * (2.0 * x104 * x106 + x136))
    result[2, 0] = numpy.sum(
        x111 * (x137 * x4 * (x36 + x38) + x139 + x42 * (x137 * x46 + x140))
    )
    result[2, 1] = numpy.sum(x123 * (x137 * x72 + x143 + x4 * (x137 * x68 + x141)))
    result[2, 2] = numpy.sum(x123 * (x146 * x4 + x151 + x4 * (x146 + x149)))
    result[2, 3] = numpy.sum(x111 * (x137 * x92 + x153))
    result[2, 4] = numpy.sum(x123 * (x154 * x62 + x156))
    result[2, 5] = numpy.sum(x111 * (x157 * x62 + x160))
    result[3, 0] = numpy.sum(x111 * (x163 * x4 + 2.0 * x164))
    result[3, 1] = numpy.sum(x123 * (x165 + x167 * x4))
    result[3, 2] = numpy.sum(x123 * (x168 + x169 * x4))
    result[3, 3] = numpy.sum(x170 * x171)
    result[3, 4] = numpy.sum(x172 * x173)
    result[3, 5] = numpy.sum(x171 * x174)
    result[4, 0] = numpy.sum(x123 * x137 * (x107 + x108))
    result[4, 1] = numpy.sum(x178 * (x175 + x177 * x4))
    result[4, 2] = numpy.sum(x178 * (x106 * x149 + x179))
    result[4, 3] = numpy.sum(x173 * x180)
    result[4, 4] = numpy.sum(x178 * x181 * x4)
    result[4, 5] = numpy.sum(x106 * x157 * x173)
    result[5, 0] = numpy.sum(x111 * (x183 * x4 + 2.0 * x184))
    result[5, 1] = numpy.sum(x123 * (x185 + x186 * x4))
    result[5, 2] = numpy.sum(x123 * (x187 + x188 * x4))
    result[5, 3] = numpy.sum(x171 * x189)
    result[5, 4] = numpy.sum(x173 * x190)
    result[5, 5] = numpy.sum(x171 * x191)
    result[6, 0] = numpy.sum(x51 * (x106 * x163 + x109))
    result[6, 1] = numpy.sum(x73 * (x106 * x167 + x122 + x164))
    result[6, 2] = numpy.sum(x73 * (x106 * x169 + x125))
    result[6, 3] = numpy.sum(x51 * (x106 * x170 + x130 + 2.0 * x165))
    result[6, 4] = numpy.sum(x73 * (x106 * x172 + x135 + x168))
    result[6, 5] = numpy.sum(x51 * (x106 * x174 + x136))
    result[7, 0] = numpy.sum(x111 * (x137 * x162 + x139))
    result[7, 1] = numpy.sum(x123 * (x106 * x177 + x137 * x166 + x143))
    result[7, 2] = numpy.sum(x123 * (x148 * x161 + x151))
    result[7, 3] = numpy.sum(x111 * (x106 * x180 + x153 + 2.0 * x175))
    result[7, 4] = numpy.sum(x123 * (x106 * x181 + x156 + x179))
    result[7, 5] = numpy.sum(x111 * (x157 * x161 + x160))
    result[8, 0] = numpy.sum(x183 * x192)
    result[8, 1] = numpy.sum(x123 * (x106 * x186 + x184))
    result[8, 2] = numpy.sum(x106 * x123 * x188)
    result[8, 3] = numpy.sum(x111 * (x106 * x189 + 2.0 * x185))
    result[8, 4] = numpy.sum(x123 * (x106 * x190 + x187))
    result[8, 5] = numpy.sum(x191 * x192)
    result[9, 0] = numpy.sum(x51 * (x137 * x183 + x138))
    result[9, 1] = numpy.sum(x73 * (x137 * x186 + x142))
    result[9, 2] = numpy.sum(x73 * (x137 * x188 + x150 + x184))
    result[9, 3] = numpy.sum(x51 * (x137 * x189 + x152))
    result[9, 4] = numpy.sum(x73 * (x137 * x190 + x155 + x185))
    result[9, 5] = numpy.sum(x51 * (x137 * x191 + x159 + 2.0 * x187))
    return result


def _2center2el3d_33(ax, da, A, bx, db, B):
    """Cartesian (f|f) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 10), dtype=float)

    x0 = ax ** (-1.0)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = -x2 * (ax * A[0] + bx * B[0])
    x4 = -x3 - A[0]
    x5 = bx ** (-1.0)
    x6 = -x3 - B[0]
    x7 = bx * x2
    x8 = ax * x7 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x9 = boys(3, x8)
    x10 = x1 ** (-1.5)
    x11 = 17.4934183276249
    x12 = x11 * x5
    x13 = 2.0 * x12
    x14 = x10 * x13
    x15 = x14 * x9
    x16 = x1 ** (-0.5)
    x17 = boys(2, x8)
    x18 = 0.5 * x5
    x19 = x18 * (2.0 * x0 * x11 * x16 * x17 * x5 - x15)
    x20 = x6**2
    x21 = boys(4, x8)
    x22 = x0 * x16
    x23 = x13 * x22
    x24 = x21 * x23
    x25 = x20 * x24
    x26 = x19 + x25
    x27 = x6 * (x26 + x5 * (2.0 * x0 * x11 * x16 * x17 * x5 - x15))
    x28 = 0.5 / (ax + bx)
    x29 = x18 * (2.0 * x0 * x11 * x16 * x5 * boys(1, x8) - x14 * x17)
    x30 = x23 * x9
    x31 = x20 * x30 + x29
    x32 = x28 * x31
    x33 = x14 * x21
    x34 = x18 * (2.0 * x0 * x11 * x16 * x5 * x9 - x33)
    x35 = boys(5, x8)
    x36 = x20 * x23
    x37 = x35 * x36
    x38 = x34 + x37
    x39 = x6 * (x38 + x5 * (2.0 * x0 * x11 * x16 * x5 * x9 - x33))
    x40 = x26 * x28
    x41 = x14 * x35
    x42 = x18 * (2.0 * x0 * x11 * x16 * x21 * x5 - x41)
    x43 = boys(6, x8)
    x44 = x36 * x43
    x45 = x6 * (x42 + x44 + x5 * (2.0 * x0 * x11 * x16 * x21 * x5 - x41))
    x46 = x4 * x45
    x47 = x28 * x38
    x48 = 3.0 * x47
    x49 = x39 * x7
    x50 = 0.5 * x0
    x51 = x50 * (x27 - x49)
    x52 = x38 * x4
    x53 = x21 * x6
    x54 = x12 * x22
    x55 = 4.0 * x28 * x54
    x56 = x53 * x55
    x57 = x52 + x56
    x58 = 3.0 * x28
    x59 = x50 * (-x26 * x7 + x31)
    x60 = 2.0 * x28
    x61 = x54 * x60
    x62 = x61 * x9
    x63 = x4 * x6
    x64 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**4.5)
    x65 = 12.0679557322504 * x64
    x66 = -x2 * (ax * A[1] + bx * B[1])
    x67 = -x66 - B[1]
    x68 = x5 * x67 * (2.0 * x0 * x11 * x16 * x17 * x5 - x15)
    x69 = x25 * x67 + 0.5 * x68
    x70 = x5 * x67 * (2.0 * x0 * x11 * x16 * x5 * x9 - x33)
    x71 = x37 * x67 + 0.5 * x70
    x72 = x53 * x67
    x73 = x55 * x72
    x74 = x6 * x67
    x75 = x55 * x9
    x76 = x74 * x75
    x77 = x5 * x67 * (2.0 * x0 * x11 * x16 * x21 * x5 - x41)
    x78 = x44 * x67 + 0.5 * x77
    x79 = x4 * x78
    x80 = x35 * x74
    x81 = x55 * x80
    x82 = x7 * x71
    x83 = x50 * (x69 - x82)
    x84 = x21 * x67
    x85 = x61 * x84
    x86 = x35 * x67
    x87 = x23 * x63
    x88 = x86 * x87
    x89 = x85 + x88
    x90 = 2.0 * x0 * x10 * x11
    x91 = x50 * (2.0 * x0 * x11 * x16 * x5 * x6 * x67 * x9 - x72 * x90)
    x92 = 26.9847693667702 * x64
    x93 = -x2 * (ax * A[2] + bx * B[2])
    x94 = -x93 - B[2]
    x95 = x5 * x94 * (2.0 * x0 * x11 * x16 * x17 * x5 - x15)
    x96 = 0.5 * x95
    x97 = x25 * x94 + x96
    x98 = x5 * x94 * (2.0 * x0 * x11 * x16 * x5 * x9 - x33)
    x99 = 0.5 * x98
    x100 = x37 * x94 + x99
    x101 = x56 * x94
    x102 = x6 * x94
    x103 = x102 * x75
    x104 = x5 * x94 * (2.0 * x0 * x11 * x16 * x21 * x5 - x41)
    x105 = 0.5 * x104
    x106 = x105 + x44 * x94
    x107 = x106 * x4
    x108 = x102 * x35
    x109 = x108 * x55
    x110 = x100 * x7
    x111 = x50 * (-x110 + x97)
    x112 = x21 * x94
    x113 = x112 * x61
    x114 = x87 * x94
    x115 = x114 * x35
    x116 = x113 + x115
    x117 = x90 * x94
    x118 = x50 * (2.0 * x0 * x11 * x16 * x5 * x6 * x9 * x94 - x117 * x53)
    x119 = x67**2
    x120 = x119 * x30 + x29
    x121 = x120 * x28
    x122 = x119 * x24 + x19
    x123 = x122 * x6
    x124 = x122 * x28
    x125 = x119 * x23
    x126 = x125 * x35
    x127 = x126 + x34
    x128 = x4**2
    x129 = x50 * (x120 - x122 * x7)
    x130 = x127 * x28
    x131 = x125 * x43
    x132 = x131 + x42
    x133 = x132 * x63
    x134 = x6 * x7
    x135 = x127 * x134
    x136 = x50 * (x122 * x6 - x135)
    x137 = x130 * x4
    x138 = x67 * x94
    x139 = x24 * x67
    x140 = x50 * (2.0 * x0 * x11 * x16 * x5 * x67 * x9 * x94 - x117 * x84)
    x141 = x86 * x94
    x142 = x141 * x23
    x143 = x141 * x61
    x144 = x50 * (2.0 * x0 * x11 * x16 * x21 * x5 * x6 * x67 * x94 - x117 * x80)
    x145 = 46.7389915737742 * x64
    x146 = x94**2
    x147 = x146 * x30 + x29
    x148 = x147 * x28
    x149 = x146 * x24 + x19
    x150 = x149 * x6
    x151 = x149 * x28
    x152 = x146 * x23
    x153 = x152 * x35 + x34
    x154 = x50 * (x147 - x149 * x7)
    x155 = x153 * x28
    x156 = x152 * x43 + x42
    x157 = x156 * x63
    x158 = x134 * x153
    x159 = x50 * (x149 * x6 - x158)
    x160 = x155 * x4
    x161 = x122 * x67 + x68
    x162 = x127 * x67 + x70
    x163 = x162 * x7
    x164 = x132 * x67 + x77
    x165 = x128 * x164
    x166 = x50 * (x161 - x163)
    x167 = x24 * x94
    x168 = x119 * x167 + x96
    x169 = x126 * x94 + x99
    x170 = x169 * x7
    x171 = x105 + x131 * x94
    x172 = x50 * (x168 - x170)
    x173 = x153 * x67
    x174 = x173 * x7
    x175 = x50 * (x149 * x67 - x174)
    x176 = x156 * x67
    x177 = x149 * x94 + x95
    x178 = x153 * x94 + x98
    x179 = x178 * x7
    x180 = x104 + x156 * x94
    x181 = x128 * x180
    x182 = x50 * (x177 - x179)
    x183 = -x66 - A[1]
    x184 = x183 * x46
    x185 = x183 * x48
    x186 = x0 * x183 * (x27 - x49)
    x187 = x183 * x78
    x188 = x187 + x47
    x189 = x53 * x61
    x190 = x23 * x80
    x191 = x183 * x190
    x192 = x189 + x191
    x193 = x192 * x60
    x194 = x0 * (x183 * x69 + x32 - x7 * (x183 * x71 + x40))
    x195 = x28 * (x139 * x183 + x62)
    x196 = 60.3397786612521 * x64
    x197 = x0 * x183 * (-x110 + x97)
    x198 = x113 * x183
    x199 = x127 * x183
    x200 = x55 * x84
    x201 = x199 + x200
    x202 = x201 * x28
    x203 = x132 * x6
    x204 = x183 * x203
    x205 = x204 + x81
    x206 = x0 * (x123 * x183 - x7 * (x199 * x6 + x73) + x76)
    x207 = x113 + x142 * x183
    x208 = x207 * x28
    x209 = x108 * x61
    x210 = x23 * x43 * x74 * x94
    x211 = x183 * x210 + x209
    x212 = x138 * x24 * x6
    x213 = x0 * (x102 * x62 + x183 * x212 - x7 * x94 * (x189 + x191))
    x214 = 104.511562358749 * x64
    x215 = x155 * x183
    x216 = x0 * x183 * (x149 * x6 - x158)
    x217 = x164 * x183
    x218 = 3.0 * x130
    x219 = x217 + x218
    x220 = x0 * (3.0 * x121 + x161 * x183 - x7 * (3.0 * x124 + x162 * x183))
    x221 = x141 * x55
    x222 = x171 * x183 + x221
    x223 = x200 * x94
    x224 = x138 * x75
    x225 = x0 * (x168 * x183 + x224 - x7 * (x169 * x183 + x223))
    x226 = x155 + x176 * x183
    x227 = x149 * x67
    x228 = x0 * (x148 + x183 * x227 - x7 * (x151 + x173 * x183))
    x229 = x0 * x183 * (x177 - x179)
    x230 = -x93 - A[2]
    x231 = x0 * x230 * (x27 - x49)
    x232 = 0.5 * x231
    x233 = x230 * x81
    x234 = x0 * x230 * (x69 - x82)
    x235 = 0.5 * x234
    x236 = x230 * x85
    x237 = x106 * x230 + x47
    x238 = x237 * x4
    x239 = x108 * x23
    x240 = x189 + x230 * x239
    x241 = x240 * x28
    x242 = 2.0 * x241
    x243 = x0 * (x230 * x97 + x32 - x7 * (x100 * x230 + x40))
    x244 = 0.5 * x243
    x245 = x28 * (x167 * x230 + x62)
    x246 = x130 * x230
    x247 = x0 * x230 * (x122 * x6 - x135)
    x248 = 0.5 * x247
    x249 = x142 * x230 + x85
    x250 = x249 * x28
    x251 = x61 * x80
    x252 = x210 * x230 + x251
    x253 = x0 * (x212 * x230 + x62 * x74 - x7 * (x190 * x230 * x94 + x61 * x72))
    x254 = 0.5 * x253
    x255 = x153 * x230
    x256 = x112 * x55 + x255
    x257 = x256 * x28
    x258 = x156 * x6
    x259 = x109 + x230 * x258
    x260 = x259 * x4
    x261 = x0 * (x103 + x150 * x230 - x7 * (x101 + x255 * x6))
    x262 = 0.5 * x261
    x263 = x0 * x230 * (x161 - x163)
    x264 = 0.5 * x263
    x265 = x130 + x171 * x230
    x266 = x0 * (x121 + x168 * x230 - x7 * (x124 + x169 * x230))
    x267 = 0.5 * x266
    x268 = x176 * x230 + x221
    x269 = x0 * (x224 + x227 * x230 - x7 * (x223 + x255 * x67))
    x270 = 0.5 * x269
    x271 = 3.0 * x155 + x180 * x230
    x272 = x0 * (3.0 * x148 + x177 * x230 - x7 * (3.0 * x151 + x178 * x230))
    x273 = 0.5 * x272
    x274 = x183**2
    x275 = x274 * x45
    x276 = x275 + x51
    x277 = x28 * (x274 * x38 + x59)
    x278 = x183 * x47
    x279 = x183 * x188 + x278 + x83
    x280 = x60 * (x183 * x189 + x183 * x192 + x91)
    x281 = x106 * x274 + x111
    x282 = x28 * (x118 + x239 * x274)
    x283 = x28 * (x129 + x183 * x201 + 2.0 * x195)
    x284 = x136 + x183 * x205 + x193
    x285 = x28 * (x140 + x183 * x207 + x198)
    x286 = x144 + x183 * x209 + x183 * x211
    x287 = x28 * (x153 * x274 + x154)
    x288 = x159 + x258 * x274
    x289 = x166 + x183 * x219 + 3.0 * x202
    x290 = x4 * x92
    x291 = x172 + x183 * x222 + 2.0 * x208
    x292 = x196 * x4
    x293 = x175 + x183 * x226 + x215
    x294 = x180 * x274 + x182
    x295 = x230 * x47
    x296 = x187 * x230 + x295
    x297 = x189 * x230
    x298 = x60 * (x191 * x230 + x297)
    x299 = x230 * x28 * (x199 + x200)
    x300 = x204 * x230 + x233
    x301 = x28 * (x183 * x249 + x245)
    x302 = x183 * x252 + x241
    x303 = x183 * x257
    x304 = x230 * (x217 + x218)
    x305 = x145 * x4
    x306 = 2.0 * x250
    x307 = x183 * x265 + x306
    x308 = x214 * x4
    x309 = x183 * x268 + x257
    x310 = x230**2
    x311 = x310 * x45 + x51
    x312 = x28 * (x310 * x38 + x59)
    x313 = x310 * x78 + x83
    x314 = x28 * (x190 * x310 + x91)
    x315 = 2.0 * x314
    x316 = x111 + x230 * x237 + x295
    x317 = x28 * (x118 + x230 * x240 + x297)
    x318 = 2.0 * x317
    x319 = x28 * (x127 * x310 + x129)
    x320 = x136 + x203 * x310
    x321 = x28 * (x140 + x230 * x249 + x236)
    x322 = x144 + x230 * x251 + x230 * x252
    x323 = x28 * (x154 + x230 * x256 + 2.0 * x245)
    x324 = x159 + x230 * x259 + x242
    x325 = x164 * x310 + x166
    x326 = x172 + x230 * x265 + x246
    x327 = x175 + x230 * x268 + x306
    x328 = x182 + x230 * x271 + 3.0 * x257
    x329 = x183 * x92
    x330 = x183 * x196
    x331 = 2.0 * x321

    # 100 item(s)
    result[0, 0] = numpy.sum(
        x65
        * (
            x0 * (x27 * x4 + 3.0 * x32 - x7 * (x39 * x4 + 3.0 * x40))
            + x4 * (x4 * (x46 + x48) + x51 + x57 * x58)
            + x58 * (x4 * x57 + x59 + x60 * (x24 * x63 + x62))
        )
    )
    result[0, 1] = numpy.sum(
        x92
        * (
            x0 * (x4 * x69 - x7 * (x4 * x71 + x73) + x76)
            + x4 * (x4 * (x79 + x81) + x60 * x89 + x83)
            + x60 * (x4 * x85 + x4 * x89 + x91)
        )
    )
    result[0, 2] = numpy.sum(
        x92
        * (
            x0 * (x103 + x4 * x97 - x7 * (x100 * x4 + x101))
            + x4 * (x111 + x116 * x60 + x4 * (x107 + x109))
            + x60 * (x113 * x4 + x116 * x4 + x118)
        )
    )
    result[0, 3] = numpy.sum(
        x92
        * (
            x0 * (x121 + x123 * x4 - x7 * (x124 + x127 * x63))
            + x28 * (x127 * x128 + x129)
            + x4 * (x136 + x137 + x4 * (x130 + x133))
        )
    )
    result[0, 4] = numpy.sum(
        x145
        * (
            x0 * (x138 * x62 + x139 * x63 * x94 - x7 * x94 * (x85 + x88))
            + x28 * (x128 * x142 + x140)
            + x4 * (x143 * x4 + x144 + x4 * (x114 * x43 * x67 + x143))
        )
    )
    result[0, 5] = numpy.sum(
        x92
        * (
            x0 * (x148 + x150 * x4 - x7 * (x151 + x153 * x63))
            + x28 * (x128 * x153 + x154)
            + x4 * (x159 + x160 + x4 * (x155 + x157))
        )
    )
    result[0, 6] = numpy.sum(x4 * x65 * (x0 * (x161 - x163) + x165 + x166))
    result[0, 7] = numpy.sum(x4 * x92 * (x0 * (x168 - x170) + x128 * x171 + x172))
    result[0, 8] = numpy.sum(x4 * x92 * (x0 * (x149 * x67 - x174) + x128 * x176 + x175))
    result[0, 9] = numpy.sum(x4 * x65 * (x0 * (x177 - x179) + x181 + x182))
    result[1, 0] = numpy.sum(
        0.5 * x92 * (2.0 * x183 * x58 * (x52 + x56) + x186 + 2.0 * x4 * (x184 + x185))
    )
    result[1, 1] = numpy.sum(
        0.5
        * x196
        * (x194 + 2.0 * x4 * (x188 * x4 + x193) + 2.0 * x60 * (x192 * x4 + x195))
    )
    result[1, 2] = numpy.sum(
        0.5
        * x196
        * (2.0 * x183 * x4 * (x107 + x109) + x197 + 2.0 * x60 * (x115 * x183 + x198))
    )
    result[1, 3] = numpy.sum(0.5 * x196 * (4.0 * x202 * x4 + 2.0 * x205 * x4**2 + x206))
    result[1, 4] = numpy.sum(0.5 * x214 * (4.0 * x208 * x4 + 2.0 * x211 * x4**2 + x213))
    result[1, 5] = numpy.sum(
        0.5 * x196 * (2.0 * x160 * x183 + x216 + 2.0 * x4 * (x157 * x183 + x215))
    )
    result[1, 6] = numpy.sum(0.5 * x92 * (2.0 * x128 * x219 + x220))
    result[1, 7] = numpy.sum(0.5 * x196 * (2.0 * x128 * x222 + x225))
    result[1, 8] = numpy.sum(0.5 * x196 * (2.0 * x128 * x226 + x228))
    result[1, 9] = numpy.sum(0.5 * x92 * (2.0 * x181 * x183 + x229))
    result[2, 0] = numpy.sum(
        x92 * (x230 * x4 * (x46 + x48) + x230 * x58 * (x52 + x56) + x232)
    )
    result[2, 1] = numpy.sum(
        x196 * (x235 + x4 * (x230 * x79 + x233) + x60 * (x230 * x88 + x236))
    )
    result[2, 2] = numpy.sum(
        x196 * (x244 + x4 * (x238 + x242) + x60 * (x240 * x4 + x245))
    )
    result[2, 3] = numpy.sum(x196 * (x137 * x230 + x248 + x4 * (x133 * x230 + x246)))
    result[2, 4] = numpy.sum(x214 * (x250 * x4 + x254 + x4 * (x250 + x252 * x4)))
    result[2, 5] = numpy.sum(x196 * (x257 * x4 + x262 + x4 * (x257 + x260)))
    result[2, 6] = numpy.sum(x92 * (x165 * x230 + x264))
    result[2, 7] = numpy.sum(x196 * (x128 * x265 + x267))
    result[2, 8] = numpy.sum(x196 * (x128 * x268 + x270))
    result[2, 9] = numpy.sum(x92 * (x128 * x271 + x273))
    result[3, 0] = numpy.sum(x92 * (x276 * x4 + 3.0 * x277))
    result[3, 1] = numpy.sum(x196 * (x279 * x4 + x280))
    result[3, 2] = numpy.sum(x196 * (x281 * x4 + 2.0 * x282))
    result[3, 3] = numpy.sum(x196 * (x283 + x284 * x4))
    result[3, 4] = numpy.sum(x214 * (x285 + x286 * x4))
    result[3, 5] = numpy.sum(x196 * (x287 + x288 * x4))
    result[3, 6] = numpy.sum(x289 * x290)
    result[3, 7] = numpy.sum(x291 * x292)
    result[3, 8] = numpy.sum(x292 * x293)
    result[3, 9] = numpy.sum(x290 * x294)
    result[4, 0] = numpy.sum(x145 * x230 * (x184 + x185))
    result[4, 1] = numpy.sum(x214 * (x296 * x4 + x298))
    result[4, 2] = numpy.sum(x183 * x214 * (x238 + x242))
    result[4, 3] = numpy.sum(x214 * (x299 + x300 * x4))
    result[4, 4] = numpy.sum(181.019335983756 * x64 * (x301 + x302 * x4))
    result[4, 5] = numpy.sum(x214 * (x183 * x260 + x303))
    result[4, 6] = numpy.sum(x304 * x305)
    result[4, 7] = numpy.sum(x307 * x308)
    result[4, 8] = numpy.sum(x308 * x309)
    result[4, 9] = numpy.sum(x183 * x271 * x305)
    result[5, 0] = numpy.sum(x92 * (x311 * x4 + 3.0 * x312))
    result[5, 1] = numpy.sum(x196 * (x313 * x4 + x315))
    result[5, 2] = numpy.sum(x196 * (x316 * x4 + x318))
    result[5, 3] = numpy.sum(x196 * (x319 + x320 * x4))
    result[5, 4] = numpy.sum(x214 * (x321 + x322 * x4))
    result[5, 5] = numpy.sum(x196 * (x323 + x324 * x4))
    result[5, 6] = numpy.sum(x290 * x325)
    result[5, 7] = numpy.sum(x292 * x326)
    result[5, 8] = numpy.sum(x292 * x327)
    result[5, 9] = numpy.sum(x290 * x328)
    result[6, 0] = numpy.sum(x65 * (x183 * x276 + x186))
    result[6, 1] = numpy.sum(x92 * (x183 * x279 + x194 + x277))
    result[6, 2] = numpy.sum(x92 * (x183 * x281 + x197))
    result[6, 3] = numpy.sum(x92 * (x183 * x284 + x206 + x280))
    result[6, 4] = numpy.sum(x145 * (x183 * x286 + x213 + x282))
    result[6, 5] = numpy.sum(x92 * (x183 * x288 + x216))
    result[6, 6] = numpy.sum(x65 * (x183 * x289 + x220 + 3.0 * x283))
    result[6, 7] = numpy.sum(x92 * (x183 * x291 + x225 + 2.0 * x285))
    result[6, 8] = numpy.sum(x92 * (x183 * x293 + x228 + x287))
    result[6, 9] = numpy.sum(x65 * (x183 * x294 + x229))
    result[7, 0] = numpy.sum(x92 * (x230 * x275 + x232))
    result[7, 1] = numpy.sum(x196 * (x183 * x296 + x230 * x278 + x235))
    result[7, 2] = numpy.sum(x196 * (x237 * x274 + x244))
    result[7, 3] = numpy.sum(x196 * (x183 * x300 + x248 + x298))
    result[7, 4] = numpy.sum(x214 * (x183 * x241 + x183 * x302 + x254))
    result[7, 5] = numpy.sum(x196 * (x259 * x274 + x262))
    result[7, 6] = numpy.sum(x92 * (x183 * x304 + x264 + 3.0 * x299))
    result[7, 7] = numpy.sum(x196 * (x183 * x307 + x267 + 2.0 * x301))
    result[7, 8] = numpy.sum(x196 * (x183 * x309 + x270 + x303))
    result[7, 9] = numpy.sum(x92 * (x271 * x274 + x273))
    result[8, 0] = numpy.sum(x311 * x329)
    result[8, 1] = numpy.sum(x196 * (x183 * x313 + x312))
    result[8, 2] = numpy.sum(x316 * x330)
    result[8, 3] = numpy.sum(x196 * (x183 * x320 + x315))
    result[8, 4] = numpy.sum(x214 * (x183 * x322 + x317))
    result[8, 5] = numpy.sum(x324 * x330)
    result[8, 6] = numpy.sum(x92 * (x183 * x325 + 3.0 * x319))
    result[8, 7] = numpy.sum(x196 * (x183 * x326 + x331))
    result[8, 8] = numpy.sum(x196 * (x183 * x327 + x323))
    result[8, 9] = numpy.sum(x328 * x329)
    result[9, 0] = numpy.sum(x65 * (x230 * x311 + x231))
    result[9, 1] = numpy.sum(x92 * (x230 * x313 + x234))
    result[9, 2] = numpy.sum(x92 * (x230 * x316 + x243 + x312))
    result[9, 3] = numpy.sum(x92 * (x230 * x320 + x247))
    result[9, 4] = numpy.sum(x145 * (x230 * x322 + x253 + x314))
    result[9, 5] = numpy.sum(x92 * (x230 * x324 + x261 + x318))
    result[9, 6] = numpy.sum(x65 * (x230 * x325 + x263))
    result[9, 7] = numpy.sum(x92 * (x230 * x326 + x266 + x319))
    result[9, 8] = numpy.sum(x92 * (x230 * x327 + x269 + x331))
    result[9, 9] = numpy.sum(x65 * (x230 * x328 + x272 + 3.0 * x323))
    return result


def _2center2el3d_34(ax, da, A, bx, db, B):
    """Cartesian (f|g) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 15), dtype=float)

    x0 = ax ** (-1.0)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = -x2 * (ax * A[0] + bx * B[0])
    x4 = -x3 - A[0]
    x5 = -x3 - B[0]
    x6 = bx ** (-1.0)
    x7 = ax * x2
    x8 = bx * x7 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x9 = boys(4, x8)
    x10 = x1 ** (-1.5)
    x11 = 17.4934183276249
    x12 = x11 * x6
    x13 = 2.0 * x12
    x14 = x10 * x13
    x15 = x14 * x9
    x16 = x15 * x5
    x17 = x1 ** (-0.5)
    x18 = boys(3, x8)
    x19 = 0.5 * x6
    x20 = x19 * (2.0 * x0 * x11 * x17 * x18 * x6 - x15)
    x21 = x5**2
    x22 = boys(5, x8)
    x23 = x0 * x17
    x24 = x13 * x23
    x25 = x22 * x24
    x26 = x21 * x25
    x27 = x20 + x26
    x28 = x27 * x5 + x6 * (2.0 * x0 * x11 * x17 * x18 * x5 * x6 - x16)
    x29 = boys(2, x8)
    x30 = x19 * (2.0 * x0 * x11 * x17 * x6 * boys(1, x8) - x14 * x29)
    x31 = x14 * x18
    x32 = x19 * (2.0 * x0 * x11 * x17 * x29 * x6 - x31)
    x33 = x24 * x9
    x34 = x21 * x33
    x35 = x32 + x34
    x36 = x18 * x24
    x37 = 1.5 * x6
    x38 = x28 * x5 + x37 * (x21 * x36 + x30 - x35 * x7)
    x39 = 0.5 / (ax + bx)
    x40 = x5 * (x35 + x6 * (2.0 * x0 * x11 * x17 * x29 * x6 - x31))
    x41 = x39 * x40
    x42 = x14 * x22
    x43 = x42 * x5
    x44 = x19 * (2.0 * x0 * x11 * x17 * x6 * x9 - x42)
    x45 = boys(6, x8)
    x46 = x21 * x24
    x47 = x45 * x46
    x48 = x44 + x47
    x49 = x48 * x5 + x6 * (2.0 * x0 * x11 * x17 * x5 * x6 * x9 - x43)
    x50 = -x37 * (x27 * x7 - x35) + x49 * x5
    x51 = x28 * x39
    x52 = bx * x2
    x53 = x14 * x45
    x54 = x5 * x53
    x55 = x19 * (2.0 * x0 * x11 * x17 * x22 * x6 - x53)
    x56 = boys(7, x8)
    x57 = x46 * x56
    x58 = x37 * (x27 - x48 * x7) + x5 * (
        x5 * (x55 + x57) + x6 * (2.0 * x0 * x11 * x17 * x22 * x5 * x6 - x54)
    )
    x59 = x4 * x58
    x60 = x39 * x49
    x61 = 4.0 * x60
    x62 = x50 * x52
    x63 = 0.5 * x0
    x64 = x63 * (x38 - x62)
    x65 = x4 * x49
    x66 = x27 * x39
    x67 = 3.0 * x66
    x68 = x65 + x67
    x69 = 4.0 * x39
    x70 = x63 * (-x28 * x52 + x40)
    x71 = x12 * x23
    x72 = x71 * x9
    x73 = x5 * x72
    x74 = 3.0 * x39
    x75 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**5.5)
    x76 = 9.12251705727742 * x75
    x77 = -x2 * (ax * A[1] + bx * B[1])
    x78 = -x77 - B[1]
    x79 = x5 * x78
    x80 = x15 * x78
    x81 = x6 * (2.0 * x0 * x11 * x17 * x18 * x6 * x78 - x80)
    x82 = x26 * x78
    x83 = 0.5 * x81 + x82
    x84 = x5 * x83 + x6 * (2.0 * x0 * x11 * x17 * x18 * x5 * x6 * x78 - x15 * x79)
    x85 = x6 * x78 * (2.0 * x0 * x11 * x17 * x29 * x6 - x31)
    x86 = x34 * x78 + 0.5 * x85
    x87 = x39 * x86
    x88 = x42 * x78
    x89 = x6 * (2.0 * x0 * x11 * x17 * x6 * x78 * x9 - x88)
    x90 = x47 * x78
    x91 = 0.5 * x89 + x90
    x92 = x5 * x91 + x6 * (2.0 * x0 * x11 * x17 * x5 * x6 * x78 * x9 - x42 * x79)
    x93 = x39 * x83
    x94 = x53 * x78
    x95 = x6 * (2.0 * x0 * x11 * x17 * x22 * x6 * x78 - x94)
    x96 = x57 * x78
    x97 = 0.5 * x5 * (x95 + 2.0 * x96) + x6 * (
        2.0 * x0 * x11 * x17 * x22 * x5 * x6 * x78 - x53 * x79
    )
    x98 = x4 * x97
    x99 = x39 * x91
    x100 = 3.0 * x99
    x101 = x52 * x92
    x102 = x63 * (-x101 + x84)
    x103 = x4 * x91
    x104 = x22 * x79
    x105 = x69 * x71
    x106 = x104 * x105
    x107 = x103 + x106
    x108 = x63 * (-x52 * x83 + x86)
    x109 = 2.0 * x39
    x110 = x109 * x72
    x111 = x110 * x78
    x112 = x4 * x5
    x113 = x25 * x78
    x114 = 24.1359114645008 * x75
    x115 = -x2 * (ax * A[2] + bx * B[2])
    x116 = -x115 - B[2]
    x117 = x116 * x6 * (2.0 * x0 * x11 * x17 * x18 * x6 - x15)
    x118 = 0.5 * x117
    x119 = x116 * x26 + x118
    x120 = x116 * x6 * (2.0 * x0 * x11 * x17 * x18 * x5 * x6 - x16) + x119 * x5
    x121 = x116 * x6 * (2.0 * x0 * x11 * x17 * x29 * x6 - x31)
    x122 = 0.5 * x121
    x123 = x116 * x34 + x122
    x124 = x123 * x39
    x125 = x116 * x6 * (2.0 * x0 * x11 * x17 * x6 * x9 - x42)
    x126 = 0.5 * x125
    x127 = x116 * x47 + x126
    x128 = x116 * x6 * (2.0 * x0 * x11 * x17 * x5 * x6 * x9 - x43) + x127 * x5
    x129 = x119 * x39
    x130 = x116 * x6 * (2.0 * x0 * x11 * x17 * x22 * x6 - x53)
    x131 = 0.5 * x130
    x132 = x116 * x6 * (2.0 * x0 * x11 * x17 * x22 * x5 * x6 - x54) + x5 * (
        x116 * x57 + x131
    )
    x133 = x132 * x4
    x134 = x127 * x39
    x135 = 3.0 * x134
    x136 = x128 * x52
    x137 = x63 * (x120 - x136)
    x138 = x127 * x4
    x139 = x116 * x22
    x140 = x139 * x5
    x141 = x105 * x140
    x142 = x138 + x141
    x143 = x63 * (-x119 * x52 + x123)
    x144 = x110 * x116
    x145 = x116 * x25
    x146 = x78**2
    x147 = x146 * x25 + x20
    x148 = x146 * x33 + x32
    x149 = x146 * x36 - x148 * x7 + x30
    x150 = x147 * x21 + x149 * x19
    x151 = x146 * x24
    x152 = x151 * x45
    x153 = x152 + x44
    x154 = -x147 * x7 + x148
    x155 = x153 * x21 + x154 * x19
    x156 = x147 * x39
    x157 = x156 * x5
    x158 = x148 * x5
    x159 = x151 * x56
    x160 = x159 + x55
    x161 = x147 - x153 * x7
    x162 = x160 * x21 + x161 * x19
    x163 = x162 * x4
    x164 = x153 * x5
    x165 = x109 * x164
    x166 = x155 * x52
    x167 = x63 * (x150 - x166)
    x168 = x112 * x153
    x169 = x156 + x168
    x170 = x5 * x52
    x171 = x63 * (-x147 * x170 + x148 * x5)
    x172 = 31.1593277158494 * x75
    x173 = x116 * x6 * (2.0 * x0 * x11 * x17 * x18 * x6 * x78 - x80)
    x174 = x116 * x82 + 0.5 * x173
    x175 = x116 * x6 * (2.0 * x0 * x11 * x17 * x6 * x78 * x9 - x88)
    x176 = x116 * x90 + 0.5 * x175
    x177 = x104 * x116
    x178 = x105 * x177
    x179 = x69 * x72
    x180 = x116 * x179
    x181 = x180 * x79
    x182 = x116 * x6 * (2.0 * x0 * x11 * x17 * x22 * x6 * x78 - x94)
    x183 = x116 * x96 + 0.5 * x182
    x184 = x116 * x45
    x185 = x105 * x184 * x79
    x186 = x63 * (x174 - x176 * x52)
    x187 = x109 * x71
    x188 = x139 * x78
    x189 = x187 * x188
    x190 = x184 * x24 * x78
    x191 = x112 * x190 + x189
    x192 = 2.0 * x0 * x11 * x63 * (-x10 * x177 + x116 * x17 * x5 * x6 * x78 * x9)
    x193 = 53.9695387335403 * x75
    x194 = x116**2
    x195 = x194 * x25 + x20
    x196 = x194 * x33 + x32
    x197 = x194 * x36 - x196 * x7 + x30
    x198 = x19 * x197
    x199 = x195 * x21 + x198
    x200 = x194 * x24
    x201 = x200 * x45 + x44
    x202 = -x195 * x7 + x196
    x203 = x19 * x202
    x204 = x201 * x21 + x203
    x205 = x195 * x39
    x206 = x205 * x5
    x207 = x196 * x5
    x208 = x200 * x56 + x55
    x209 = x195 - x201 * x7
    x210 = x19 * x209
    x211 = x208 * x21 + x210
    x212 = x211 * x4
    x213 = x201 * x5
    x214 = x204 * x52
    x215 = x63 * (x199 - x214)
    x216 = x112 * x201
    x217 = x205 + x216
    x218 = x195 * x5
    x219 = x63 * (x196 * x5 - x218 * x52)
    x220 = x148 * x78 + x85
    x221 = x220 * x39
    x222 = x147 * x78 + x81
    x223 = x222 * x5
    x224 = x222 * x39
    x225 = x153 * x78 + x89
    x226 = x4**2
    x227 = x63 * (x220 - x222 * x52)
    x228 = x225 * x39
    x229 = x160 * x78 + x95
    x230 = x112 * x229
    x231 = x63 * (-x170 * x225 + x222 * x5)
    x232 = x228 * x4
    x233 = x116 * x146 * x33 + x122
    x234 = x233 * x39
    x235 = x118 + x145 * x146
    x236 = x235 * x5
    x237 = x235 * x39
    x238 = x116 * x152 + x126
    x239 = x63 * (x233 - x235 * x52)
    x240 = x238 * x39
    x241 = x116 * x159 + x131
    x242 = x63 * (-x170 * x238 + x235 * x5)
    x243 = x196 * x78
    x244 = x205 * x78
    x245 = x195 * x78
    x246 = x63 * (x196 * x78 - x245 * x52)
    x247 = x201 * x78
    x248 = x247 * x39
    x249 = x208 * x78
    x250 = x63 * (-x170 * x247 + x195 * x5 * x78)
    x251 = x116 * x196 + x121
    x252 = x251 * x39
    x253 = x116 * x195 + x117
    x254 = x253 * x5
    x255 = x253 * x39
    x256 = x116 * x201 + x125
    x257 = x63 * (x251 - x253 * x52)
    x258 = x256 * x39
    x259 = x116 * x208 + x130
    x260 = x112 * x259
    x261 = x63 * (-x170 * x256 + x253 * x5)
    x262 = x258 * x4
    x263 = x149 * x37 + x222 * x78
    x264 = x154 * x37 + x225 * x78
    x265 = x264 * x52
    x266 = x161 * x37 + x229 * x78
    x267 = x226 * x266
    x268 = x63 * (x263 - x265)
    x269 = x173 + x235 * x78
    x270 = x175 + x238 * x78
    x271 = x270 * x52
    x272 = x182 + x241 * x78
    x273 = x63 * (x269 - x271)
    x274 = x146 * x195 + x198
    x275 = x146 * x201 + x203
    x276 = x275 * x52
    x277 = x146 * x208 + x210
    x278 = x63 * (x274 - x276)
    x279 = x256 * x78
    x280 = x279 * x52
    x281 = x63 * (x253 * x78 - x280)
    x282 = x259 * x78
    x283 = x116 * x253 + x197 * x37
    x284 = x116 * x256 + x202 * x37
    x285 = x284 * x52
    x286 = x116 * x259 + x209 * x37
    x287 = x226 * x286
    x288 = x63 * (x283 - x285)
    x289 = -x77 - A[1]
    x290 = x289 * x59
    x291 = x289 * x61
    x292 = x0 * x289 * (x38 - x62)
    x293 = 20.3985682659737 * x75
    x294 = x289 * x97
    x295 = x294 + x60
    x296 = x289 * x91
    x297 = x296 + x66
    x298 = x0 * (x289 * x84 + x41 - x52 * (x289 * x92 + x51))
    x299 = x109 * x73
    x300 = x289 * x5
    x301 = x109 * (x113 * x300 + x299)
    x302 = x0 * x289 * (x120 - x136)
    x303 = x139 * x300
    x304 = x162 * x289
    x305 = 2.0 * x99
    x306 = x304 + x305
    x307 = x153 * x300
    x308 = x106 + x307
    x309 = x0 * (x150 * x289 - x52 * (x155 * x289 + 2.0 * x93) + 2.0 * x87)
    x310 = x39 * (x147 * x289 + x179 * x78)
    x311 = 69.6743749058326 * x75
    x312 = x134 + x183 * x289
    x313 = x140 * x187 + x190 * x300
    x314 = x109 * x313
    x315 = x0 * (x124 + x174 * x289 - x52 * (x129 + x176 * x289))
    x316 = x113 * x116
    x317 = x39 * (x144 + x289 * x316)
    x318 = 120.679557322504 * x75
    x319 = x201 * x300
    x320 = x0 * x289 * (x199 - x214)
    x321 = x205 * x289
    x322 = x225 * x289
    x323 = 3.0 * x156
    x324 = x322 + x323
    x325 = x324 * x39
    x326 = x229 * x300
    x327 = x164 * x74
    x328 = x326 + x327
    x329 = x0 * (x158 * x74 + x223 * x289 - x5 * x52 * (x322 + x323))
    x330 = x238 * x289
    x331 = x105 * x188
    x332 = x330 + x331
    x333 = x332 * x39
    x334 = x185 + x241 * x300
    x335 = x0 * (x181 + x236 * x289 - x52 * (x178 + x330 * x5))
    x336 = x205 + x247 * x289
    x337 = x336 * x39
    x338 = x213 * x39 + x249 * x300
    x339 = x0 * (x207 * x39 + x218 * x289 * x78 - x52 * (x206 + x247 * x300))
    x340 = x258 * x289
    x341 = x0 * (x253 * x289 * x5 - x256 * x300 * x52)
    x342 = x266 * x289
    x343 = 4.0 * x228
    x344 = x342 + x343
    x345 = x0 * (4.0 * x221 + x263 * x289 - x52 * (4.0 * x224 + x264 * x289))
    x346 = 3.0 * x240 + x272 * x289
    x347 = x0 * (3.0 * x234 + x269 * x289 - x52 * (3.0 * x237 + x270 * x289))
    x348 = x109 * x247 + x277 * x289
    x349 = x0 * (x109 * x243 + x274 * x289 - x52 * (2.0 * x244 + x275 * x289))
    x350 = x258 + x282 * x289
    x351 = x253 * x78
    x352 = x0 * (x252 + x289 * x351 - x52 * (x255 + x279 * x289))
    x353 = x0 * x289 * (x283 - x285)
    x354 = -x115 - A[2]
    x355 = x0 * x354 * (x38 - x62)
    x356 = 0.5 * x355
    x357 = x0 * x354 * (-x101 + x84)
    x358 = 0.5 * x357
    x359 = x354 * x5
    x360 = x22 * x359 * x78
    x361 = x105 * x360
    x362 = x132 * x354 + x60
    x363 = x362 * x4
    x364 = x127 * x354 + x66
    x365 = x364 * x39
    x366 = 3.0 * x365
    x367 = x0 * (x120 * x354 + x41 - x52 * (x128 * x354 + x51))
    x368 = 0.5 * x367
    x369 = x39 * (x145 * x359 + x299)
    x370 = 2.0 * x369
    x371 = x0 * x354 * (x150 - x166)
    x372 = 0.5 * x371
    x373 = x156 * x354
    x374 = x183 * x354 + x99
    x375 = x104 * x187 + x190 * x359
    x376 = x109 * x375
    x377 = x0 * (x174 * x354 - x52 * (x176 * x354 + x93) + x87)
    x378 = 0.5 * x377
    x379 = x39 * (x111 + x316 * x354)
    x380 = 2.0 * x134 + x211 * x354
    x381 = x380 * x4
    x382 = x141 + x213 * x354
    x383 = x382 * x39
    x384 = 2.0 * x383
    x385 = x0 * (2.0 * x124 + x199 * x354 - x52 * (2.0 * x129 + x204 * x354))
    x386 = 0.5 * x385
    x387 = x195 * x354
    x388 = x39 * (x180 + x387)
    x389 = x228 * x354
    x390 = x0 * (x222 * x354 * x5 - x225 * x359 * x52)
    x391 = 0.5 * x390
    x392 = x238 * x354
    x393 = x156 + x392
    x394 = x39 * x393
    x395 = x164 * x39
    x396 = x241 * x359 + x395
    x397 = x0 * (x158 * x39 + x236 * x354 - x52 * (x157 + x392 * x5))
    x398 = 0.5 * x397
    x399 = x247 * x354 + x331
    x400 = x39 * x399
    x401 = x185 + x249 * x359
    x402 = x0 * (x181 + x387 * x79 - x52 * (x178 + x247 * x359))
    x403 = 0.5 * x402
    x404 = x256 * x354
    x405 = 3.0 * x205 + x404
    x406 = x39 * x405
    x407 = x213 * x74 + x259 * x359
    x408 = x4 * x407
    x409 = x0 * (x207 * x74 + x254 * x354 - x52 * (3.0 * x206 + x404 * x5))
    x410 = 0.5 * x409
    x411 = x0 * x354 * (x263 - x265)
    x412 = 0.5 * x411
    x413 = x228 + x272 * x354
    x414 = x0 * (x221 + x269 * x354 - x52 * (x224 + x270 * x354))
    x415 = 0.5 * x414
    x416 = 2.0 * x240 + x277 * x354
    x417 = x0 * (2.0 * x234 + x274 * x354 - x52 * (2.0 * x237 + x275 * x354))
    x418 = 0.5 * x417
    x419 = x247 * x74 + x282 * x354
    x420 = x0 * (x243 * x74 + x351 * x354 - x52 * (3.0 * x244 + x404 * x78))
    x421 = 0.5 * x420
    x422 = 4.0 * x258 + x286 * x354
    x423 = x0 * (4.0 * x252 + x283 * x354 - x52 * (4.0 * x255 + x284 * x354))
    x424 = 0.5 * x423
    x425 = x289**2
    x426 = x425 * x58
    x427 = x426 + x64
    x428 = x39 * (x425 * x49 + x70)
    x429 = x289 * x60
    x430 = x102 + x289 * x295 + x429
    x431 = x108 + x289 * x297 + x289 * x66
    x432 = x132 * x425 + x137
    x433 = x39 * (x127 * x425 + x143)
    x434 = x109 * x297 + x167 + x289 * x306
    x435 = x171 + x289 * x308 + x301
    x436 = x134 * x289 + x186 + x289 * x312
    x437 = x109 * (x187 * x303 + x192 + x289 * x313)
    x438 = x211 * x425 + x215
    x439 = x39 * (x213 * x425 + x219)
    x440 = x39 * (x227 + x289 * x324 + 3.0 * x310)
    x441 = x231 + x289 * x328 + x308 * x74
    x442 = x39 * (x239 + x289 * x332 + 2.0 * x317)
    x443 = x242 + x289 * x334 + x314
    x444 = x39 * (x246 + x289 * x336 + x321)
    x445 = x250 + x289 * x338 + x319 * x39
    x446 = x39 * (x256 * x425 + x257)
    x447 = x259 * x425 * x5 + x261
    x448 = x268 + x289 * x344 + 4.0 * x325
    x449 = x293 * x4
    x450 = x273 + x289 * x346 + 3.0 * x333
    x451 = x193 * x4
    x452 = x278 + x289 * x348 + 2.0 * x337
    x453 = x311 * x4
    x454 = x281 + x289 * x350 + x340
    x455 = x286 * x425 + x288
    x456 = 35.3313566383285 * x75
    x457 = x354 * x60
    x458 = x294 * x354 + x457
    x459 = x354 * x66
    x460 = x296 * x354 + x459
    x461 = 93.4779831475484 * x75
    x462 = x354 * (x304 + x305)
    x463 = x307 * x354 + x361
    x464 = 120.679557322504 * x75
    x465 = x289 * x374 + x365
    x466 = x109 * (x289 * x375 + x369)
    x467 = 209.023124717498 * x75
    x468 = x354 * x39 * (x322 + x323)
    x469 = x354 * (x326 + x327)
    x470 = 2.0 * x379
    x471 = x39 * (x289 * x393 + x470)
    x472 = x289 * x396 + x376
    x473 = x39 * (x289 * x399 + x388)
    x474 = x289 * x401 + x383
    x475 = x289 * x406
    x476 = x354 * (x342 + x343)
    x477 = x4 * x456
    x478 = x289 * x413 + 3.0 * x394
    x479 = x4 * x461
    x480 = x289 * x416 + 2.0 * x400
    x481 = x289 * x419 + x406
    x482 = x354**2
    x483 = x482 * x58 + x64
    x484 = x39 * (x482 * x49 + x70)
    x485 = x102 + x482 * x97
    x486 = x39 * (x108 + x482 * x91)
    x487 = x137 + x354 * x362 + x457
    x488 = x39 * (x143 + x354 * x364 + x459)
    x489 = x162 * x482 + x167
    x490 = x39 * (x164 * x482 + x171)
    x491 = x186 + x354 * x374 + x354 * x99
    x492 = x109 * (x187 * x360 + x192 + x354 * x375)
    x493 = x215 + x354 * x380 + 2.0 * x365
    x494 = x39 * (x219 + x354 * x382 + x370)
    x495 = x39 * (x225 * x482 + x227)
    x496 = x229 * x482 * x5 + x231
    x497 = x39 * (x239 + x354 * x393 + x373)
    x498 = x242 + x354 * x395 + x354 * x396
    x499 = x39 * (x246 + x354 * x399 + x470)
    x500 = x250 + x354 * x401 + x376
    x501 = x39 * (x257 + x354 * x405 + 3.0 * x388)
    x502 = x261 + x354 * x407 + 3.0 * x383
    x503 = x266 * x482 + x268
    x504 = x273 + x354 * x413 + x389
    x505 = x278 + x354 * x416 + 2.0 * x394
    x506 = x281 + x354 * x419 + 3.0 * x400
    x507 = x288 + x354 * x422 + 4.0 * x406
    x508 = x289 * x293
    x509 = x193 * x289

    # 150 item(s)
    result[0, 0] = numpy.sum(
        x76
        * (
            x0 * (x38 * x4 + 4.0 * x41 - x52 * (x4 * x50 + 4.0 * x51))
            + x4 * (x4 * (x59 + x61) + x64 + x68 * x69)
            + x69 * (x4 * x68 + x70 + x74 * (x27 * x4 + x69 * x73))
        )
    )
    result[0, 1] = numpy.sum(
        x114
        * (
            x0 * (x4 * x84 - x52 * (x4 * x92 + 3.0 * x93) + 3.0 * x87)
            + x4 * (x102 + x107 * x74 + x4 * (x100 + x98))
            + x74 * (x107 * x4 + x108 + x109 * (x111 + x112 * x113))
        )
    )
    result[0, 2] = numpy.sum(
        x114
        * (
            x0 * (x120 * x4 + 3.0 * x124 - x52 * (x128 * x4 + 3.0 * x129))
            + x4 * (x137 + x142 * x74 + x4 * (x133 + x135))
            + x74 * (x109 * (x112 * x145 + x144) + x142 * x4 + x143)
        )
    )
    result[0, 3] = numpy.sum(
        x172
        * (
            x0 * (x109 * x158 + x150 * x4 - x52 * (x155 * x4 + 2.0 * x157))
            + x109 * (x156 * x4 + x169 * x4 + x171)
            + x4 * (x109 * x169 + x167 + x4 * (x163 + x165))
        )
    )
    result[0, 4] = numpy.sum(
        x193
        * (
            x0 * (x174 * x4 + x181 - x52 * (x176 * x4 + x178))
            + x109 * (x189 * x4 + x191 * x4 + x192)
            + x4 * (x109 * x191 + x186 + x4 * (x183 * x4 + x185))
        )
    )
    result[0, 5] = numpy.sum(
        x172
        * (
            x0 * (x109 * x207 + x199 * x4 - x52 * (x204 * x4 + 2.0 * x206))
            + x109 * (x205 * x4 + x217 * x4 + x219)
            + x4 * (x109 * x217 + x215 + x4 * (x109 * x213 + x212))
        )
    )
    result[0, 6] = numpy.sum(
        x114
        * (
            x0 * (x221 + x223 * x4 - x52 * (x112 * x225 + x224))
            + x39 * (x225 * x226 + x227)
            + x4 * (x231 + x232 + x4 * (x228 + x230))
        )
    )
    result[0, 7] = numpy.sum(
        x193
        * (
            x0 * (x234 + x236 * x4 - x52 * (x112 * x238 + x237))
            + x39 * (x226 * x238 + x239)
            + x4 * (x240 * x4 + x242 + x4 * (x112 * x241 + x240))
        )
    )
    result[0, 8] = numpy.sum(
        x193
        * (
            x0 * (x112 * x245 + x243 * x39 - x52 * (x216 * x78 + x244))
            + x39 * (x226 * x247 + x246)
            + x4 * (x248 * x4 + x250 + x4 * (x112 * x249 + x248))
        )
    )
    result[0, 9] = numpy.sum(
        x114
        * (
            x0 * (x252 + x254 * x4 - x52 * (x112 * x256 + x255))
            + x39 * (x226 * x256 + x257)
            + x4 * (x261 + x262 + x4 * (x258 + x260))
        )
    )
    result[0, 10] = numpy.sum(x4 * x76 * (x0 * (x263 - x265) + x267 + x268))
    result[0, 11] = numpy.sum(x114 * x4 * (x0 * (x269 - x271) + x226 * x272 + x273))
    result[0, 12] = numpy.sum(x172 * x4 * (x0 * (x274 - x276) + x226 * x277 + x278))
    result[0, 13] = numpy.sum(x114 * x4 * (x0 * (x253 * x78 - x280) + x226 * x282 + x281))
    result[0, 14] = numpy.sum(x4 * x76 * (x0 * (x283 - x285) + x287 + x288))
    result[1, 0] = numpy.sum(
        0.5 * x293 * (2.0 * x289 * x69 * (x65 + x67) + x292 + 2.0 * x4 * (x290 + x291))
    )
    result[1, 1] = numpy.sum(
        0.5
        * x193
        * (x298 + 2.0 * x4 * (x295 * x4 + x297 * x74) + 2.0 * x74 * (x297 * x4 + x301))
    )
    result[1, 2] = numpy.sum(
        0.5
        * x193
        * (
            2.0 * x289 * x4 * (x133 + x135)
            + x302
            + 2.0 * x74 * (x105 * x303 + x138 * x289)
        )
    )
    result[1, 3] = numpy.sum(
        0.5
        * x311
        * (2.0 * x109 * (x308 * x4 + x310) + x309 + 2.0 * x4 * (x109 * x308 + x306 * x4))
    )
    result[1, 4] = numpy.sum(
        0.5
        * x318
        * (2.0 * x109 * (x313 * x4 + x317) + x315 + 2.0 * x4 * (x312 * x4 + x314))
    )
    result[1, 5] = numpy.sum(
        0.5
        * x311
        * (
            2.0 * x109 * (x216 * x289 + x321)
            + x320
            + 2.0 * x4 * (x109 * x319 + x212 * x289)
        )
    )
    result[1, 6] = numpy.sum(0.5 * x193 * (4.0 * x325 * x4 + 2.0 * x328 * x4**2 + x329))
    result[1, 7] = numpy.sum(0.5 * x318 * (4.0 * x333 * x4 + 2.0 * x334 * x4**2 + x335))
    result[1, 8] = numpy.sum(0.5 * x318 * (4.0 * x337 * x4 + 2.0 * x338 * x4**2 + x339))
    result[1, 9] = numpy.sum(
        0.5 * x193 * (2.0 * x262 * x289 + x341 + 2.0 * x4 * (x260 * x289 + x340))
    )
    result[1, 10] = numpy.sum(0.5 * x293 * (2.0 * x226 * x344 + x345))
    result[1, 11] = numpy.sum(0.5 * x193 * (2.0 * x226 * x346 + x347))
    result[1, 12] = numpy.sum(0.5 * x311 * (2.0 * x226 * x348 + x349))
    result[1, 13] = numpy.sum(0.5 * x193 * (2.0 * x226 * x350 + x352))
    result[1, 14] = numpy.sum(0.5 * x293 * (2.0 * x287 * x289 + x353))
    result[2, 0] = numpy.sum(
        x293 * (x354 * x4 * (x59 + x61) + x354 * x69 * (x65 + x67) + x356)
    )
    result[2, 1] = numpy.sum(
        x193 * (x354 * x4 * (x100 + x98) + x358 + x74 * (x103 * x354 + x361))
    )
    result[2, 2] = numpy.sum(
        x193 * (x368 + x4 * (x363 + x366) + x74 * (x364 * x4 + x370))
    )
    result[2, 3] = numpy.sum(
        x311 * (x109 * (x168 * x354 + x373) + x354 * x4 * (x163 + x165) + x372)
    )
    result[2, 4] = numpy.sum(
        x318 * (x109 * (x375 * x4 + x379) + x378 + x4 * (x374 * x4 + x376))
    )
    result[2, 5] = numpy.sum(
        x311 * (x109 * (x382 * x4 + x388) + x386 + x4 * (x381 + x384))
    )
    result[2, 6] = numpy.sum(x193 * (x232 * x354 + x391 + x4 * (x230 * x354 + x389)))
    result[2, 7] = numpy.sum(x318 * (x394 * x4 + x398 + x4 * (x394 + x396 * x4)))
    result[2, 8] = numpy.sum(x318 * (x4 * x400 + x4 * (x4 * x401 + x400) + x403))
    result[2, 9] = numpy.sum(x193 * (x4 * x406 + x4 * (x406 + x408) + x410))
    result[2, 10] = numpy.sum(x293 * (x267 * x354 + x412))
    result[2, 11] = numpy.sum(x193 * (x226 * x413 + x415))
    result[2, 12] = numpy.sum(x311 * (x226 * x416 + x418))
    result[2, 13] = numpy.sum(x193 * (x226 * x419 + x421))
    result[2, 14] = numpy.sum(x293 * (x226 * x422 + x424))
    result[3, 0] = numpy.sum(x293 * (x4 * x427 + 4.0 * x428))
    result[3, 1] = numpy.sum(x193 * (x4 * x430 + x431 * x74))
    result[3, 2] = numpy.sum(x193 * (x4 * x432 + 3.0 * x433))
    result[3, 3] = numpy.sum(x311 * (x109 * x435 + x4 * x434))
    result[3, 4] = numpy.sum(x318 * (x4 * x436 + x437))
    result[3, 5] = numpy.sum(x311 * (x4 * x438 + 2.0 * x439))
    result[3, 6] = numpy.sum(x193 * (x4 * x441 + x440))
    result[3, 7] = numpy.sum(x318 * (x4 * x443 + x442))
    result[3, 8] = numpy.sum(x318 * (x4 * x445 + x444))
    result[3, 9] = numpy.sum(x193 * (x4 * x447 + x446))
    result[3, 10] = numpy.sum(x448 * x449)
    result[3, 11] = numpy.sum(x450 * x451)
    result[3, 12] = numpy.sum(x452 * x453)
    result[3, 13] = numpy.sum(x451 * x454)
    result[3, 14] = numpy.sum(x449 * x455)
    result[4, 0] = numpy.sum(x354 * x456 * (x290 + x291))
    result[4, 1] = numpy.sum(x461 * (x4 * x458 + x460 * x74))
    result[4, 2] = numpy.sum(x289 * x461 * (x363 + x366))
    result[4, 3] = numpy.sum(x464 * (x109 * x463 + x4 * x462))
    result[4, 4] = numpy.sum(x467 * (x4 * x465 + x466))
    result[4, 5] = numpy.sum(x289 * x464 * (x381 + x384))
    result[4, 6] = numpy.sum(x461 * (x4 * x469 + x468))
    result[4, 7] = numpy.sum(x467 * (x4 * x472 + x471))
    result[4, 8] = numpy.sum(x467 * (x4 * x474 + x473))
    result[4, 9] = numpy.sum(x461 * (x289 * x408 + x475))
    result[4, 10] = numpy.sum(x476 * x477)
    result[4, 11] = numpy.sum(x478 * x479)
    result[4, 12] = numpy.sum(x4 * x464 * x480)
    result[4, 13] = numpy.sum(x479 * x481)
    result[4, 14] = numpy.sum(x289 * x422 * x477)
    result[5, 0] = numpy.sum(x293 * (x4 * x483 + 4.0 * x484))
    result[5, 1] = numpy.sum(x193 * (x4 * x485 + 3.0 * x486))
    result[5, 2] = numpy.sum(x193 * (x4 * x487 + 3.0 * x488))
    result[5, 3] = numpy.sum(x311 * (x4 * x489 + 2.0 * x490))
    result[5, 4] = numpy.sum(x318 * (x4 * x491 + x492))
    result[5, 5] = numpy.sum(x311 * (x4 * x493 + 2.0 * x494))
    result[5, 6] = numpy.sum(x193 * (x4 * x496 + x495))
    result[5, 7] = numpy.sum(x318 * (x4 * x498 + x497))
    result[5, 8] = numpy.sum(x318 * (x4 * x500 + x499))
    result[5, 9] = numpy.sum(x193 * (x4 * x502 + x501))
    result[5, 10] = numpy.sum(x449 * x503)
    result[5, 11] = numpy.sum(x451 * x504)
    result[5, 12] = numpy.sum(x453 * x505)
    result[5, 13] = numpy.sum(x451 * x506)
    result[5, 14] = numpy.sum(x449 * x507)
    result[6, 0] = numpy.sum(x76 * (x289 * x427 + x292))
    result[6, 1] = numpy.sum(x114 * (x289 * x430 + x298 + x428))
    result[6, 2] = numpy.sum(x114 * (x289 * x432 + x302))
    result[6, 3] = numpy.sum(x172 * (x109 * x431 + x289 * x434 + x309))
    result[6, 4] = numpy.sum(x193 * (x289 * x436 + x315 + x433))
    result[6, 5] = numpy.sum(x172 * (x289 * x438 + x320))
    result[6, 6] = numpy.sum(x114 * (x289 * x441 + x329 + x435 * x74))
    result[6, 7] = numpy.sum(x193 * (x289 * x443 + x335 + x437))
    result[6, 8] = numpy.sum(x193 * (x289 * x445 + x339 + x439))
    result[6, 9] = numpy.sum(x114 * (x289 * x447 + x341))
    result[6, 10] = numpy.sum(x76 * (x289 * x448 + x345 + 4.0 * x440))
    result[6, 11] = numpy.sum(x114 * (x289 * x450 + x347 + 3.0 * x442))
    result[6, 12] = numpy.sum(x172 * (x289 * x452 + x349 + 2.0 * x444))
    result[6, 13] = numpy.sum(x114 * (x289 * x454 + x352 + x446))
    result[6, 14] = numpy.sum(x76 * (x289 * x455 + x353))
    result[7, 0] = numpy.sum(x293 * (x354 * x426 + x356))
    result[7, 1] = numpy.sum(x193 * (x289 * x458 + x354 * x429 + x358))
    result[7, 2] = numpy.sum(x193 * (x362 * x425 + x368))
    result[7, 3] = numpy.sum(x311 * (x109 * x460 + x289 * x462 + x372))
    result[7, 4] = numpy.sum(x318 * (x289 * x365 + x289 * x465 + x378))
    result[7, 5] = numpy.sum(x311 * (x380 * x425 + x386))
    result[7, 6] = numpy.sum(x193 * (x289 * x469 + x391 + x463 * x74))
    result[7, 7] = numpy.sum(x318 * (x289 * x472 + x398 + x466))
    result[7, 8] = numpy.sum(x318 * (x289 * x383 + x289 * x474 + x403))
    result[7, 9] = numpy.sum(x193 * (x407 * x425 + x410))
    result[7, 10] = numpy.sum(x293 * (x289 * x476 + x412 + 4.0 * x468))
    result[7, 11] = numpy.sum(x193 * (x289 * x478 + x415 + 3.0 * x471))
    result[7, 12] = numpy.sum(x311 * (x289 * x480 + x418 + 2.0 * x473))
    result[7, 13] = numpy.sum(x193 * (x289 * x481 + x421 + x475))
    result[7, 14] = numpy.sum(x293 * (x422 * x425 + x424))
    result[8, 0] = numpy.sum(x483 * x508)
    result[8, 1] = numpy.sum(x193 * (x289 * x485 + x484))
    result[8, 2] = numpy.sum(x487 * x509)
    result[8, 3] = numpy.sum(x311 * (x289 * x489 + 2.0 * x486))
    result[8, 4] = numpy.sum(x318 * (x289 * x491 + x488))
    result[8, 5] = numpy.sum(x289 * x311 * x493)
    result[8, 6] = numpy.sum(x193 * (x289 * x496 + 3.0 * x490))
    result[8, 7] = numpy.sum(x318 * (x289 * x498 + x492))
    result[8, 8] = numpy.sum(x318 * (x289 * x500 + x494))
    result[8, 9] = numpy.sum(x502 * x509)
    result[8, 10] = numpy.sum(x293 * (x289 * x503 + 4.0 * x495))
    result[8, 11] = numpy.sum(x193 * (x289 * x504 + 3.0 * x497))
    result[8, 12] = numpy.sum(x311 * (x289 * x505 + 2.0 * x499))
    result[8, 13] = numpy.sum(x193 * (x289 * x506 + x501))
    result[8, 14] = numpy.sum(x507 * x508)
    result[9, 0] = numpy.sum(x76 * (x354 * x483 + x355))
    result[9, 1] = numpy.sum(x114 * (x354 * x485 + x357))
    result[9, 2] = numpy.sum(x114 * (x354 * x487 + x367 + x484))
    result[9, 3] = numpy.sum(x172 * (x354 * x489 + x371))
    result[9, 4] = numpy.sum(x193 * (x354 * x491 + x377 + x486))
    result[9, 5] = numpy.sum(x172 * (x354 * x493 + x385 + 2.0 * x488))
    result[9, 6] = numpy.sum(x114 * (x354 * x496 + x390))
    result[9, 7] = numpy.sum(x193 * (x354 * x498 + x397 + x490))
    result[9, 8] = numpy.sum(x193 * (x354 * x500 + x402 + x492))
    result[9, 9] = numpy.sum(x114 * (x354 * x502 + x409 + 3.0 * x494))
    result[9, 10] = numpy.sum(x76 * (x354 * x503 + x411))
    result[9, 11] = numpy.sum(x114 * (x354 * x504 + x414 + x495))
    result[9, 12] = numpy.sum(x172 * (x354 * x505 + x417 + 2.0 * x497))
    result[9, 13] = numpy.sum(x114 * (x354 * x506 + x420 + 3.0 * x499))
    result[9, 14] = numpy.sum(x76 * (x354 * x507 + x423 + 4.0 * x501))
    return result


def _2center2el3d_35(ax, da, A, bx, db, B):
    """Cartesian (f|h) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 21), dtype=float)

    x0 = ax ** (-1.0)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = -x2 * (ax * A[0] + bx * B[0])
    x4 = -x3 - A[0]
    x5 = -x3 - B[0]
    x6 = bx ** (-1.0)
    x7 = ax * x2
    x8 = bx * x7 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x9 = boys(5, x8)
    x10 = 17.4934183276249
    x11 = 2.0 * x6
    x12 = x10 * x11
    x13 = x1 ** (-1.5) * x12
    x14 = x13 * x9
    x15 = x14 * x5
    x16 = x1 ** (-0.5)
    x17 = boys(4, x8)
    x18 = 0.5 * x6
    x19 = x18 * (2.0 * x0 * x10 * x16 * x17 * x6 - x14)
    x20 = x5**2
    x21 = boys(6, x8)
    x22 = x0 * x16
    x23 = x12 * x22
    x24 = x21 * x23
    x25 = x20 * x24
    x26 = x19 + x25
    x27 = x26 * x5 + x6 * (2.0 * x0 * x10 * x16 * x17 * x5 * x6 - x15)
    x28 = x13 * x17
    x29 = boys(3, x8)
    x30 = x18 * (2.0 * x0 * x10 * x16 * x29 * x6 - x28)
    x31 = x23 * x9
    x32 = x20 * x31
    x33 = x30 + x32
    x34 = x13 * x29
    x35 = boys(2, x8)
    x36 = x18 * (2.0 * x0 * x10 * x16 * x35 * x6 - x34)
    x37 = x17 * x23
    x38 = x20 * x37
    x39 = x36 + x38
    x40 = 1.5 * x6
    x41 = x27 * x5 - x40 * (x33 * x7 - x39)
    x42 = x28 * x5
    x43 = x33 * x5 + x6 * (2.0 * x0 * x10 * x16 * x29 * x5 * x6 - x42)
    x44 = (
        x11 * (x39 * x5 - x43 * x7 + x5 * x6 * (2.0 * x0 * x10 * x16 * x35 * x6 - x34))
        + x41 * x5
    )
    x45 = 0.5 / (ax + bx)
    x46 = x18 * (2.0 * x0 * x10 * x16 * x6 * boys(1, x8) - x13 * x35)
    x47 = x23 * x29
    x48 = x40 * (x20 * x47 - x39 * x7 + x46) + x43 * x5
    x49 = x45 * x48
    x50 = x13 * x21
    x51 = x5 * x50
    x52 = x18 * (2.0 * x0 * x10 * x16 * x6 * x9 - x50)
    x53 = boys(7, x8)
    x54 = x20 * x23
    x55 = x53 * x54
    x56 = x52 + x55
    x57 = x5 * x56 + x6 * (2.0 * x0 * x10 * x16 * x5 * x6 * x9 - x51)
    x58 = -x40 * (x26 * x7 - x33) + x5 * x57
    x59 = -x11 * (x27 * x7 - x43) + x5 * x58
    x60 = x41 * x45
    x61 = bx * x2
    x62 = x13 * x53
    x63 = x5 * x62
    x64 = x18 * (2.0 * x0 * x10 * x16 * x21 * x6 - x62)
    x65 = boys(8, x8)
    x66 = x54 * x65
    x67 = x11 * (x27 - x57 * x7) + x5 * (
        x40 * (x26 - x56 * x7)
        + x5 * (x5 * (x64 + x66) + x6 * (2.0 * x0 * x10 * x16 * x21 * x5 * x6 - x63))
    )
    x68 = x4 * x67
    x69 = x45 * x58
    x70 = 5.0 * x69
    x71 = x59 * x61
    x72 = 0.5 * x0
    x73 = x72 * (x44 - x71)
    x74 = x4 * x58
    x75 = x27 * x45
    x76 = 4.0 * x75
    x77 = x74 + x76
    x78 = 5.0 * x45
    x79 = x72 * (-x41 * x61 + x48)
    x80 = x33 * x45
    x81 = 4.0 * x45
    x82 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**6.5)
    x83 = 6.08167803818495 * x82
    x84 = -x2 * (ax * A[1] + bx * B[1])
    x85 = -x84 - B[1]
    x86 = x5 * x85
    x87 = x14 * x85
    x88 = x6 * (2.0 * x0 * x10 * x16 * x17 * x6 * x85 - x87)
    x89 = x25 * x85
    x90 = 0.5 * x88 + x89
    x91 = x5 * x90 + x6 * (2.0 * x0 * x10 * x16 * x17 * x5 * x6 * x85 - x14 * x86)
    x92 = x6 * x85 * (2.0 * x0 * x10 * x16 * x35 * x6 - x34)
    x93 = x28 * x85
    x94 = x6 * (2.0 * x0 * x10 * x16 * x29 * x6 * x85 - x93)
    x95 = x32 * x85
    x96 = 0.5 * x94 + x95
    x97 = x38 * x40 * x85 - x40 * x7 * x96 + 0.5 * x40 * x92 + x5 * x91
    x98 = x5 * x96 + x6 * (2.0 * x0 * x10 * x16 * x29 * x5 * x6 * x85 - x28 * x86)
    x99 = x45 * x98
    x100 = x50 * x85
    x101 = x6 * (2.0 * x0 * x10 * x16 * x6 * x85 * x9 - x100)
    x102 = x55 * x85
    x103 = 0.5 * x101 + x102
    x104 = x103 * x5 + x6 * (2.0 * x0 * x10 * x16 * x5 * x6 * x85 * x9 - x50 * x86)
    x105 = x104 * x5 - x40 * (x7 * x90 - x96)
    x106 = x45 * x91
    x107 = x62 * x85
    x108 = x6 * (2.0 * x0 * x10 * x16 * x21 * x6 * x85 - x107)
    x109 = x66 * x85
    x110 = -x40 * (x103 * x7 - x90) + 0.5 * x5 * (
        x5 * (x108 + 2.0 * x109)
        + 2.0 * x6 * (2.0 * x0 * x10 * x16 * x21 * x5 * x6 * x85 - x62 * x86)
    )
    x111 = x110 * x4
    x112 = x104 * x45
    x113 = 4.0 * x112
    x114 = x105 * x61
    x115 = x72 * (-x114 + x97)
    x116 = x104 * x4
    x117 = x45 * x90
    x118 = 3.0 * x117
    x119 = x116 + x118
    x120 = x72 * (-x61 * x91 + x98)
    x121 = x10 * x22 * x6 * x81 * x9
    x122 = x121 * x86
    x123 = 3.0 * x45
    x124 = 18.2450341145548 * x82
    x125 = -x2 * (ax * A[2] + bx * B[2])
    x126 = -x125 - B[2]
    x127 = x126 * x14
    x128 = x6 * (2.0 * x0 * x10 * x126 * x16 * x17 * x6 - x127)
    x129 = 0.5 * x128
    x130 = x126 * x25 + x129
    x131 = x126 * x6 * (2.0 * x0 * x10 * x16 * x17 * x5 * x6 - x15) + x130 * x5
    x132 = x126 * x6 * (2.0 * x0 * x10 * x16 * x35 * x6 - x34)
    x133 = 0.5 * x132
    x134 = x126 * x6 * (2.0 * x0 * x10 * x16 * x29 * x6 - x28)
    x135 = 0.5 * x134
    x136 = x126 * x32 + x135
    x137 = x131 * x5 + x40 * (x126 * x38 + x133 - x136 * x7)
    x138 = x126 * x6 * (2.0 * x0 * x10 * x16 * x29 * x5 * x6 - x42) + x136 * x5
    x139 = x138 * x45
    x140 = x126 * x50
    x141 = x6 * (2.0 * x0 * x10 * x126 * x16 * x6 * x9 - x140)
    x142 = 0.5 * x141
    x143 = x126 * x55 + x142
    x144 = x126 * x6 * (2.0 * x0 * x10 * x16 * x5 * x6 * x9 - x51) + x143 * x5
    x145 = x144 * x5 - x40 * (x130 * x7 - x136)
    x146 = x131 * x45
    x147 = x126 * x62
    x148 = x6 * (2.0 * x0 * x10 * x126 * x16 * x21 * x6 - x147)
    x149 = 0.5 * x148
    x150 = x40 * (x130 - x143 * x7) + x5 * (
        x126 * x6 * (2.0 * x0 * x10 * x16 * x21 * x5 * x6 - x63)
        + x5 * (x126 * x66 + x149)
    )
    x151 = x150 * x4
    x152 = x144 * x45
    x153 = 4.0 * x152
    x154 = x145 * x61
    x155 = x72 * (x137 - x154)
    x156 = x144 * x4
    x157 = x130 * x45
    x158 = 3.0 * x157
    x159 = x156 + x158
    x160 = x72 * (-x131 * x61 + x138)
    x161 = x121 * x126
    x162 = x161 * x5
    x163 = x85**2
    x164 = x163 * x37 + x36
    x165 = x163 * x31 + x30
    x166 = x163 * x24 + x19
    x167 = x164 - x165 * x7
    x168 = x166 * x20 + x167 * x18
    x169 = x5 * (x168 + x6 * (x164 - x165 * x7))
    x170 = x163 * x47 - x164 * x7 + x46
    x171 = x165 * x20 + x170 * x18
    x172 = x171 * x45
    x173 = 3.0 * x172
    x174 = x166 * x5
    x175 = x163 * x23
    x176 = x175 * x53
    x177 = x176 + x52
    x178 = x165 - x166 * x7
    x179 = x177 * x20 + x178 * x18
    x180 = x179 * x5 + x6 * (x165 * x5 - x174 * x7)
    x181 = x168 * x45
    x182 = 3.0 * x181
    x183 = x177 * x7
    x184 = x175 * x65
    x185 = x184 + x64
    x186 = x166 - x183
    x187 = x5 * (x18 * x186 + x185 * x20 + x6 * (x166 - x183))
    x188 = x187 * x4
    x189 = x179 * x45
    x190 = 3.0 * x189
    x191 = x180 * x61
    x192 = x72 * (x169 - x191)
    x193 = x179 * x4
    x194 = 2.0 * x45
    x195 = x174 * x194
    x196 = x193 + x195
    x197 = x72 * (-x168 * x61 + x171)
    x198 = x165 * x45
    x199 = 27.869749962333 * x82
    x200 = x126 * x6 * (2.0 * x0 * x10 * x16 * x17 * x6 * x85 - x87)
    x201 = x126 * x89 + 0.5 * x200
    x202 = x201 * x5 + x6 * (
        2.0 * x0 * x10 * x126 * x16 * x17 * x5 * x6 * x85 - x127 * x86
    )
    x203 = x126 * x6 * (2.0 * x0 * x10 * x16 * x29 * x6 * x85 - x93)
    x204 = x126 * x95 + 0.5 * x203
    x205 = x126 * x6 * (2.0 * x0 * x10 * x16 * x6 * x85 * x9 - x100)
    x206 = x102 * x126 + 0.5 * x205
    x207 = x206 * x5 + x6 * (
        2.0 * x0 * x10 * x126 * x16 * x5 * x6 * x85 * x9 - x140 * x86
    )
    x208 = x126 * x6 * (2.0 * x0 * x10 * x16 * x21 * x6 * x85 - x107)
    x209 = 0.5 * x5 * (2.0 * x109 * x126 + x208) + x6 * (
        2.0 * x0 * x10 * x126 * x16 * x21 * x5 * x6 * x85 - x147 * x86
    )
    x210 = x72 * (x202 - x207 * x61)
    x211 = x10 * x126 * x21 * x22 * x6 * x81 * x86
    x212 = x206 * x4 + x211
    x213 = x72 * (-x201 * x61 + x204)
    x214 = x4 * x5
    x215 = x214 * x85
    x216 = x126 * x24
    x217 = 48.2718229290016 * x82
    x218 = x126**2
    x219 = x218 * x37 + x36
    x220 = x218 * x31 + x30
    x221 = x220 * x5
    x222 = x19 + x218 * x24
    x223 = x20 * x222
    x224 = x219 - x220 * x7
    x225 = x18 * x224
    x226 = x223 + x225
    x227 = x226 * x5 + x6 * (x219 * x5 - x221 * x7)
    x228 = x218 * x47 - x219 * x7 + x46
    x229 = x18 * x228
    x230 = x20 * x220 + x229
    x231 = x230 * x45
    x232 = 3.0 * x231
    x233 = x222 * x5
    x234 = x218 * x23
    x235 = x234 * x53 + x52
    x236 = x20 * x235
    x237 = x220 - x222 * x7
    x238 = x18 * x237
    x239 = x236 + x238
    x240 = x239 * x5 + x6 * (x220 * x5 - x233 * x7)
    x241 = x226 * x45
    x242 = 3.0 * x241
    x243 = x235 * x7
    x244 = x234 * x65 + x64
    x245 = x20 * x244
    x246 = x222 - x243
    x247 = x18 * x246
    x248 = x5 * (x245 + x247 + x6 * (x222 - x243))
    x249 = x248 * x4
    x250 = x239 * x45
    x251 = 3.0 * x250
    x252 = x240 * x61
    x253 = x72 * (x227 - x252)
    x254 = x239 * x4
    x255 = x194 * x233
    x256 = x254 + x255
    x257 = x72 * (-x226 * x61 + x230)
    x258 = x220 * x45
    x259 = x166 * x85 + x88
    x260 = x165 * x85 + x94
    x261 = x164 * x85 - x260 * x7 + x92
    x262 = x18 * x261 + x20 * x259
    x263 = x101 + x177 * x85
    x264 = -x259 * x7 + x260
    x265 = x18 * x264 + x20 * x263
    x266 = x259 * x45
    x267 = x266 * x5
    x268 = x260 * x5
    x269 = x108 + x185 * x85
    x270 = x259 - x263 * x7
    x271 = x18 * x270 + x20 * x269
    x272 = x271 * x4
    x273 = x263 * x5
    x274 = x265 * x61
    x275 = x72 * (x262 - x274)
    x276 = x214 * x263
    x277 = x266 + x276
    x278 = x5 * x61
    x279 = x72 * (-x259 * x278 + x260 * x5)
    x280 = 27.869749962333 * x82
    x281 = x129 + x163 * x216
    x282 = x126 * x163 * x31 + x135
    x283 = x126 * x163 * x37 + x133 - x282 * x7
    x284 = x18 * x283 + x20 * x281
    x285 = x126 * x176 + x142
    x286 = -x281 * x7 + x282
    x287 = x18 * x286 + x20 * x285
    x288 = x281 * x45
    x289 = 2.0 * x288
    x290 = x289 * x5
    x291 = x282 * x5
    x292 = x194 * x291
    x293 = x126 * x184 + x149
    x294 = x281 - x285 * x7
    x295 = x18 * x294 + x20 * x293
    x296 = x285 * x5
    x297 = x194 * x296
    x298 = x72 * (x284 - x287 * x61)
    x299 = x214 * x285 + x288
    x300 = x72 * (-x278 * x281 + x282 * x5)
    x301 = 62.3186554316989 * x82
    x302 = x6 * x85 * (x219 - x220 * x7)
    x303 = x223 * x85 + 0.5 * x302
    x304 = x222 * x85
    x305 = x6 * (x220 * x85 - x304 * x7)
    x306 = x236 * x85 + 0.5 * x305
    x307 = x255 * x85
    x308 = x221 * x85
    x309 = x194 * x308
    x310 = x6 * x85 * (x222 - x243)
    x311 = x245 * x85 + 0.5 * x310
    x312 = x235 * x86
    x313 = x194 * x312
    x314 = x72 * (x303 - x306 * x61)
    x315 = x304 * x45
    x316 = x215 * x235 + x315
    x317 = x61 * x85
    x318 = x72 * (x220 * x5 * x85 - x233 * x317)
    x319 = x126 * x222 + x128
    x320 = x126 * x220 + x134
    x321 = x126 * x219 + x132 - x320 * x7
    x322 = x18 * x321
    x323 = x20 * x319 + x322
    x324 = x126 * x235 + x141
    x325 = -x319 * x7 + x320
    x326 = x18 * x325
    x327 = x20 * x324 + x326
    x328 = x319 * x45
    x329 = x328 * x5
    x330 = x320 * x5
    x331 = x126 * x244 + x148
    x332 = x319 - x324 * x7
    x333 = x18 * x332
    x334 = x20 * x331 + x333
    x335 = x334 * x4
    x336 = x324 * x5
    x337 = x327 * x61
    x338 = x72 * (x323 - x337)
    x339 = x214 * x324
    x340 = x328 + x339
    x341 = x319 * x5
    x342 = x72 * (x320 * x5 - x341 * x61)
    x343 = x170 * x40 + x260 * x85
    x344 = x343 * x45
    x345 = x167 * x40 + x259 * x85
    x346 = x345 * x5
    x347 = x345 * x45
    x348 = x178 * x40 + x263 * x85
    x349 = x4**2
    x350 = x72 * (x343 - x345 * x61)
    x351 = x348 * x45
    x352 = x186 * x40 + x269 * x85
    x353 = x214 * x352
    x354 = x72 * (-x278 * x348 + x345 * x5)
    x355 = x351 * x4
    x356 = x203 + x282 * x85
    x357 = x356 * x45
    x358 = x200 + x281 * x85
    x359 = x358 * x5
    x360 = x358 * x45
    x361 = x205 + x285 * x85
    x362 = x72 * (x356 - x358 * x61)
    x363 = x361 * x45
    x364 = x208 + x293 * x85
    x365 = x72 * (-x278 * x361 + x358 * x5)
    x366 = x163 * x220 + x229
    x367 = x366 * x45
    x368 = x163 * x222 + x225
    x369 = x368 * x5
    x370 = x368 * x45
    x371 = x163 * x235 + x238
    x372 = x72 * (x366 - x368 * x61)
    x373 = x371 * x45
    x374 = x163 * x244 + x247
    x375 = x72 * (-x278 * x371 + x368 * x5)
    x376 = x320 * x85
    x377 = x328 * x85
    x378 = x72 * (-x317 * x319 + x320 * x85)
    x379 = x324 * x85
    x380 = x379 * x45
    x381 = x72 * (-x278 * x379 + x319 * x5 * x85)
    x382 = x126 * x320 + x228 * x40
    x383 = x382 * x45
    x384 = x126 * x319 + x224 * x40
    x385 = x384 * x5
    x386 = x384 * x45
    x387 = x126 * x324 + x237 * x40
    x388 = x72 * (x382 - x384 * x61)
    x389 = x387 * x45
    x390 = x126 * x331 + x246 * x40
    x391 = x214 * x390
    x392 = x72 * (-x278 * x387 + x384 * x5)
    x393 = x389 * x4
    x394 = x11 * x261 + x345 * x85
    x395 = x11 * x264 + x348 * x85
    x396 = x395 * x61
    x397 = x11 * x270 + x352 * x85
    x398 = x349 * x397
    x399 = x72 * (x394 - x396)
    x400 = x283 * x40 + x358 * x85
    x401 = x286 * x40 + x361 * x85
    x402 = x401 * x61
    x403 = x294 * x40 + x364 * x85
    x404 = x72 * (x400 - x402)
    x405 = x302 + x368 * x85
    x406 = x305 + x371 * x85
    x407 = x406 * x61
    x408 = x310 + x374 * x85
    x409 = x72 * (x405 - x407)
    x410 = x163 * x319 + x322
    x411 = x163 * x324 + x326
    x412 = x411 * x61
    x413 = x163 * x331 + x333
    x414 = x72 * (x410 - x412)
    x415 = x387 * x85
    x416 = x415 * x61
    x417 = x72 * (x384 * x85 - x416)
    x418 = x390 * x85
    x419 = x11 * x321 + x126 * x384
    x420 = x11 * x325 + x126 * x387
    x421 = x420 * x61
    x422 = x11 * x332 + x126 * x390
    x423 = x349 * x422
    x424 = x72 * (x419 - x421)
    x425 = -x84 - A[1]
    x426 = x425 * x68
    x427 = x425 * x70
    x428 = x0 * x425 * (x44 - x71)
    x429 = 13.5990455106491 * x82
    x430 = x110 * x425
    x431 = x430 + x69
    x432 = x104 * x425
    x433 = x432 + x75
    x434 = x0 * (x425 * x97 + x49 - x61 * (x105 * x425 + x60))
    x435 = x425 * x90 + x80
    x436 = 40.7971365319473 * x82
    x437 = x0 * x425 * (x137 - x154)
    x438 = x187 * x425
    x439 = 2.0 * x112
    x440 = x438 + x439
    x441 = x179 * x425
    x442 = 2.0 * x117
    x443 = x441 + x442
    x444 = x123 * x443
    x445 = x0 * (x169 * x425 - x61 * (2.0 * x106 + x180 * x425) + 2.0 * x99)
    x446 = x122 + x174 * x425
    x447 = x152 + x209 * x425
    x448 = x157 + x206 * x425
    x449 = x0 * (x139 + x202 * x425 - x61 * (x146 + x207 * x425))
    x450 = x425 * x5
    x451 = x450 * x85
    x452 = x194 * (x126 * x31 * x45 * x5 + x216 * x451)
    x453 = 107.939077467081 * x82
    x454 = x0 * x425 * (x227 - x252)
    x455 = x271 * x425
    x456 = x190 + x455
    x457 = x263 * x450
    x458 = x123 * x174
    x459 = x457 + x458
    x460 = x0 * (x173 + x262 * x425 - x61 * (x182 + x265 * x425))
    x461 = x45 * (3.0 * x198 + x259 * x425)
    x462 = 62.3186554316989 * x82
    x463 = x194 * x206
    x464 = x295 * x425 + x463
    x465 = x211 + x285 * x450
    x466 = x194 * x204
    x467 = x194 * x201
    x468 = x0 * (x284 * x425 + x466 - x61 * (x287 * x425 + x467))
    x469 = x161 * x85
    x470 = x45 * (x281 * x425 + x469)
    x471 = 139.348749811665 * x82
    x472 = x250 + x311 * x425
    x473 = x233 * x45
    x474 = x235 * x451 + x473
    x475 = x194 * x474
    x476 = x0 * (x231 + x303 * x425 - x61 * (x241 + x306 * x425))
    x477 = x45 * (x258 + x304 * x425)
    x478 = x324 * x450
    x479 = x0 * x425 * (x323 - x337)
    x480 = x328 * x425
    x481 = x348 * x425
    x482 = 4.0 * x266
    x483 = x481 + x482
    x484 = x45 * x483
    x485 = x352 * x450
    x486 = x273 * x81 + x485
    x487 = x0 * (x268 * x81 + x346 * x425 - x5 * x61 * (x481 + x482))
    x488 = x361 * x425
    x489 = 3.0 * x288
    x490 = x488 + x489
    x491 = x45 * x490
    x492 = x123 * x296 + x364 * x450
    x493 = x0 * (x123 * x291 + x359 * x425 - x5 * x61 * (x488 + x489))
    x494 = x371 * x425
    x495 = x194 * x304 + x494
    x496 = x45 * x495
    x497 = x313 + x374 * x450
    x498 = x0 * (x309 + x369 * x425 - x61 * (x307 + x494 * x5))
    x499 = x328 + x379 * x425
    x500 = x45 * x499
    x501 = x331 * x451 + x336 * x45
    x502 = x0 * (x330 * x45 + x341 * x425 * x85 - x61 * (x329 + x379 * x450))
    x503 = x389 * x425
    x504 = x0 * (x384 * x425 * x5 - x387 * x450 * x61)
    x505 = x397 * x425
    x506 = 5.0 * x351
    x507 = x505 + x506
    x508 = x0 * (5.0 * x344 + x394 * x425 - x61 * (5.0 * x347 + x395 * x425))
    x509 = 4.0 * x363 + x403 * x425
    x510 = x0 * (4.0 * x357 + x400 * x425 - x61 * (4.0 * x360 + x401 * x425))
    x511 = 3.0 * x373
    x512 = x408 * x425 + x511
    x513 = 3.0 * x367
    x514 = 3.0 * x370
    x515 = x0 * (x405 * x425 + x513 - x61 * (x406 * x425 + x514))
    x516 = x194 * x379 + x413 * x425
    x517 = x0 * (x194 * x376 + x410 * x425 - x61 * (2.0 * x377 + x411 * x425))
    x518 = x389 + x418 * x425
    x519 = x384 * x85
    x520 = x0 * (x383 + x425 * x519 - x61 * (x386 + x415 * x425))
    x521 = x0 * x425 * (x419 - x421)
    x522 = -x125 - A[2]
    x523 = x0 * x522 * (x44 - x71)
    x524 = 0.5 * x523
    x525 = x0 * x522 * (-x114 + x97)
    x526 = 0.5 * x525
    x527 = x150 * x522 + x69
    x528 = x4 * x527
    x529 = x144 * x522 + x75
    x530 = x45 * x529
    x531 = 4.0 * x530
    x532 = x0 * (x137 * x522 + x49 - x61 * (x145 * x522 + x60))
    x533 = 0.5 * x532
    x534 = x45 * (x130 * x522 + x80)
    x535 = x190 * x522
    x536 = x0 * x522 * (x169 - x191)
    x537 = 0.5 * x536
    x538 = x112 + x209 * x522
    x539 = x117 + x206 * x522
    x540 = x0 * (x202 * x522 - x61 * (x106 + x207 * x522) + x99)
    x541 = 0.5 * x540
    x542 = x5 * x522
    x543 = x542 * x85
    x544 = x194 * (x216 * x543 + x31 * x45 * x86)
    x545 = 2.0 * x152 + x248 * x522
    x546 = x4 * x545
    x547 = 2.0 * x157 + x239 * x522
    x548 = x45 * x547
    x549 = 3.0 * x548
    x550 = x0 * (2.0 * x139 + x227 * x522 - x61 * (2.0 * x146 + x240 * x522))
    x551 = 0.5 * x550
    x552 = x45 * (x162 + x233 * x522)
    x553 = x263 * x542
    x554 = x0 * x522 * (x262 - x274)
    x555 = 0.5 * x554
    x556 = x266 * x522
    x557 = x189 + x295 * x522
    x558 = x174 * x45
    x559 = x285 * x542 + x558
    x560 = x194 * x559
    x561 = x0 * (x172 + x284 * x522 - x61 * (x181 + x287 * x522))
    x562 = 0.5 * x561
    x563 = x45 * (x198 + x281 * x522)
    x564 = x311 * x522 + x463
    x565 = x211 + x235 * x543
    x566 = x194 * x565
    x567 = x0 * (x303 * x522 + x466 - x61 * (x306 * x522 + x467))
    x568 = 0.5 * x567
    x569 = x45 * (x304 * x522 + x469)
    x570 = x251 + x334 * x522
    x571 = x4 * x570
    x572 = x123 * x233
    x573 = x336 * x522 + x572
    x574 = x45 * x573
    x575 = 2.0 * x574
    x576 = x0 * (x232 + x323 * x522 - x61 * (x242 + x327 * x522))
    x577 = 0.5 * x576
    x578 = x319 * x522
    x579 = x45 * (3.0 * x258 + x578)
    x580 = x351 * x522
    x581 = x0 * (x345 * x5 * x522 - x348 * x542 * x61)
    x582 = 0.5 * x581
    x583 = x361 * x522
    x584 = x266 + x583
    x585 = x45 * x584
    x586 = x273 * x45 + x364 * x542
    x587 = x0 * (x268 * x45 + x359 * x522 - x61 * (x267 + x5 * x583))
    x588 = 0.5 * x587
    x589 = x371 * x522
    x590 = x289 + x589
    x591 = x45 * x590
    x592 = x297 + x374 * x542
    x593 = x0 * (x292 + x369 * x522 - x61 * (x290 + x5 * x589))
    x594 = 0.5 * x593
    x595 = x123 * x304 + x379 * x522
    x596 = x45 * x595
    x597 = x123 * x312 + x331 * x543
    x598 = x0 * (x123 * x308 + x578 * x86 - x61 * (x379 * x542 + x572 * x85))
    x599 = 0.5 * x598
    x600 = x387 * x522
    x601 = 4.0 * x328 + x600
    x602 = x45 * x601
    x603 = x336 * x81 + x390 * x542
    x604 = x4 * x603
    x605 = x0 * (x330 * x81 + x385 * x522 - x61 * (4.0 * x329 + x5 * x600))
    x606 = 0.5 * x605
    x607 = x0 * x522 * (x394 - x396)
    x608 = 0.5 * x607
    x609 = x351 + x403 * x522
    x610 = x0 * (x344 + x400 * x522 - x61 * (x347 + x401 * x522))
    x611 = 0.5 * x610
    x612 = 2.0 * x363 + x408 * x522
    x613 = x0 * (2.0 * x357 + x405 * x522 - x61 * (2.0 * x360 + x406 * x522))
    x614 = 0.5 * x613
    x615 = x413 * x522 + x511
    x616 = x0 * (x410 * x522 + x513 - x61 * (x411 * x522 + x514))
    x617 = 0.5 * x616
    x618 = x379 * x81 + x418 * x522
    x619 = x0 * (x376 * x81 + x519 * x522 - x61 * (4.0 * x377 + x600 * x85))
    x620 = 0.5 * x619
    x621 = 5.0 * x389 + x422 * x522
    x622 = x0 * (5.0 * x383 + x419 * x522 - x61 * (5.0 * x386 + x420 * x522))
    x623 = 0.5 * x622
    x624 = x425**2
    x625 = x624 * x67
    x626 = x625 + x73
    x627 = x45 * (x58 * x624 + x79)
    x628 = x425 * x69
    x629 = x115 + x425 * x431 + x628
    x630 = x120 + x425 * x433 + x425 * x75
    x631 = x150 * x624 + x155
    x632 = x45 * (x144 * x624 + x160)
    x633 = x192 + x194 * x433 + x425 * x440
    x634 = x123 * (x194 * x435 + x197 + x425 * x443)
    x635 = x152 * x425 + x210 + x425 * x447
    x636 = x157 * x425 + x213 + x425 * x448
    x637 = x248 * x624 + x253
    x638 = x45 * (x239 * x624 + x257)
    x639 = x275 + x425 * x456 + x444
    x640 = x123 * x446 + x279 + x425 * x459
    x641 = x194 * x448 + x298 + x425 * x464
    x642 = x300 + x425 * x465 + x452
    x643 = x250 * x425 + x314 + x425 * x472
    x644 = x194 * (x318 + x425 * x473 + x425 * x474)
    x645 = x334 * x624 + x338
    x646 = x45 * (x336 * x624 + x342)
    x647 = x45 * (x350 + x425 * x483 + 4.0 * x461)
    x648 = x354 + x425 * x486 + x459 * x81
    x649 = x45 * (x362 + x425 * x490 + 3.0 * x470)
    x650 = x123 * x465 + x365 + x425 * x492
    x651 = x45 * (x372 + x425 * x495 + 2.0 * x477)
    x652 = x375 + x425 * x497 + x475
    x653 = x45 * (x378 + x425 * x499 + x480)
    x654 = x381 + x425 * x501 + x45 * x478
    x655 = x45 * (x387 * x624 + x388)
    x656 = x390 * x5 * x624 + x392
    x657 = x399 + x425 * x507 + 5.0 * x484
    x658 = x4 * x429
    x659 = x404 + x425 * x509 + 4.0 * x491
    x660 = x4 * x436
    x661 = x409 + x425 * x512 + 3.0 * x496
    x662 = x301 * x4
    x663 = x414 + x425 * x516 + 2.0 * x500
    x664 = x4 * x462
    x665 = x417 + x425 * x518 + x503
    x666 = x422 * x624 + x424
    x667 = 23.5542377588857 * x82
    x668 = x522 * x69
    x669 = x430 * x522 + x668
    x670 = x522 * x75
    x671 = x432 * x522 + x670
    x672 = 70.6627132766571 * x82
    x673 = x522 * (x438 + x439)
    x674 = x123 * x522 * (x441 + x442)
    x675 = x425 * x538 + x530
    x676 = x425 * x539 + x534
    x677 = 186.955966295097 * x82
    x678 = x455 * x522 + x535
    x679 = x522 * (x457 + x458)
    x680 = 107.939077467081 * x82
    x681 = x194 * x539
    x682 = x425 * x557 + x681
    x683 = x425 * x559 + x544
    x684 = 241.359114645008 * x82
    x685 = x425 * x564 + x548
    x686 = x194 * (x425 * x565 + x552)
    x687 = x45 * x522 * (x481 + x482)
    x688 = x485 * x522 + x553 * x81
    x689 = x45 * (x425 * x584 + 3.0 * x563)
    x690 = x123 * x559 + x425 * x586
    x691 = x45 * (x425 * x590 + 2.0 * x569)
    x692 = x425 * x592 + x566
    x693 = x45 * (x425 * x595 + x579)
    x694 = x425 * x597 + x574
    x695 = x425 * x602
    x696 = x522 * (x505 + x506)
    x697 = x4 * x667
    x698 = x425 * x609 + 4.0 * x585
    x699 = x4 * x672
    x700 = 3.0 * x591
    x701 = x425 * x612 + x700
    x702 = x425 * x615 + 2.0 * x596
    x703 = x425 * x618 + x602
    x704 = x522**2
    x705 = x67 * x704 + x73
    x706 = x45 * (x58 * x704 + x79)
    x707 = x110 * x704 + x115
    x708 = x45 * (x104 * x704 + x120)
    x709 = x155 + x522 * x527 + x668
    x710 = x45 * (x160 + x522 * x529 + x670)
    x711 = x187 * x704 + x192
    x712 = x45 * (x179 * x704 + x197)
    x713 = 3.0 * x712
    x714 = x112 * x522 + x210 + x522 * x538
    x715 = x117 * x522 + x213 + x522 * x539
    x716 = x253 + x522 * x545 + 2.0 * x530
    x717 = x45 * (x257 + x522 * x547 + 2.0 * x534)
    x718 = 3.0 * x717
    x719 = x271 * x704 + x275
    x720 = x45 * (x273 * x704 + x279)
    x721 = x189 * x522 + x298 + x522 * x557
    x722 = x300 + x522 * x558 + x522 * x559
    x723 = x194 * x722
    x724 = x314 + x522 * x564 + x681
    x725 = x318 + x522 * x565 + x544
    x726 = x194 * x725
    x727 = x338 + x522 * x570 + x549
    x728 = x45 * (x342 + x522 * x573 + 3.0 * x552)
    x729 = x45 * (x348 * x704 + x350)
    x730 = x352 * x5 * x704 + x354
    x731 = x45 * (x362 + x522 * x584 + x556)
    x732 = x365 + x45 * x553 + x522 * x586
    x733 = x45 * (x372 + x522 * x590 + 2.0 * x563)
    x734 = x375 + x522 * x592 + x560
    x735 = x45 * (x378 + x522 * x595 + 3.0 * x569)
    x736 = x123 * x565 + x381 + x522 * x597
    x737 = x45 * (x388 + x522 * x601 + 4.0 * x579)
    x738 = x392 + x522 * x603 + 4.0 * x574
    x739 = x397 * x704 + x399
    x740 = x404 + x522 * x609 + x580
    x741 = x409 + x522 * x612 + 2.0 * x585
    x742 = x414 + x522 * x615 + x700
    x743 = x417 + x522 * x618 + 4.0 * x596
    x744 = x424 + x522 * x621 + 5.0 * x602
    x745 = x425 * x429
    x746 = x425 * x436
    x747 = x194 * x715
    x748 = 3.0 * x733

    # 210 item(s)
    result[0, 0] = numpy.sum(
        x83
        * (
            x0 * (x4 * x44 + 5.0 * x49 - x61 * (x4 * x59 + 5.0 * x60))
            + x4 * (x4 * (x68 + x70) + x73 + x77 * x78)
            + x78 * (x4 * x77 + x79 + x81 * (x27 * x4 + 3.0 * x80))
        )
    )
    result[0, 1] = numpy.sum(
        x124
        * (
            x0 * (x4 * x97 - x61 * (x105 * x4 + 4.0 * x106) + 4.0 * x99)
            + x4 * (x115 + x119 * x81 + x4 * (x111 + x113))
            + x81 * (x119 * x4 + x120 + x123 * (x122 + x4 * x90))
        )
    )
    result[0, 2] = numpy.sum(
        x124
        * (
            x0 * (x137 * x4 + 4.0 * x139 - x61 * (x145 * x4 + 4.0 * x146))
            + x4 * (x155 + x159 * x81 + x4 * (x151 + x153))
            + x81 * (x123 * (x130 * x4 + x162) + x159 * x4 + x160)
        )
    )
    result[0, 3] = numpy.sum(
        x199
        * (
            x0 * (x169 * x4 + x173 - x61 * (x180 * x4 + x182))
            + x123 * (x194 * (x174 * x4 + x198) + x196 * x4 + x197)
            + x4 * (x123 * x196 + x192 + x4 * (x188 + x190))
        )
    )
    result[0, 4] = numpy.sum(
        x217
        * (
            x0 * (x123 * x204 + x202 * x4 - x61 * (x123 * x201 + x207 * x4))
            + x123 * (x194 * (x126 * x31 * x45 * x85 + x215 * x216) + x212 * x4 + x213)
            + x4 * (x123 * x212 + x210 + x4 * (x123 * x206 + x209 * x4))
        )
    )
    result[0, 5] = numpy.sum(
        x199
        * (
            x0 * (x227 * x4 + x232 - x61 * (x240 * x4 + x242))
            + x123 * (x194 * (x233 * x4 + x258) + x256 * x4 + x257)
            + x4 * (x123 * x256 + x253 + x4 * (x249 + x251))
        )
    )
    result[0, 6] = numpy.sum(
        x280
        * (
            x0 * (x194 * x268 + x262 * x4 - x61 * (x265 * x4 + 2.0 * x267))
            + x194 * (x266 * x4 + x277 * x4 + x279)
            + x4 * (x194 * x277 + x275 + x4 * (x194 * x273 + x272))
        )
    )
    result[0, 7] = numpy.sum(
        x301
        * (
            x0 * (x284 * x4 + x292 - x61 * (x287 * x4 + x290))
            + x194 * (x288 * x4 + x299 * x4 + x300)
            + x4 * (x194 * x299 + x298 + x4 * (x295 * x4 + x297))
        )
    )
    result[0, 8] = numpy.sum(
        x301
        * (
            x0 * (x303 * x4 + x309 - x61 * (x306 * x4 + x307))
            + x194 * (x315 * x4 + x316 * x4 + x318)
            + x4 * (x194 * x316 + x314 + x4 * (x311 * x4 + x313))
        )
    )
    result[0, 9] = numpy.sum(
        x280
        * (
            x0 * (x194 * x330 + x323 * x4 - x61 * (x327 * x4 + 2.0 * x329))
            + x194 * (x328 * x4 + x340 * x4 + x342)
            + x4 * (x194 * x340 + x338 + x4 * (x194 * x336 + x335))
        )
    )
    result[0, 10] = numpy.sum(
        x124
        * (
            x0 * (x344 + x346 * x4 - x61 * (x214 * x348 + x347))
            + x4 * (x354 + x355 + x4 * (x351 + x353))
            + x45 * (x348 * x349 + x350)
        )
    )
    result[0, 11] = numpy.sum(
        x217
        * (
            x0 * (x357 + x359 * x4 - x61 * (x214 * x361 + x360))
            + x4 * (x363 * x4 + x365 + x4 * (x214 * x364 + x363))
            + x45 * (x349 * x361 + x362)
        )
    )
    result[0, 12] = numpy.sum(
        x301
        * (
            x0 * (x367 + x369 * x4 - x61 * (x214 * x371 + x370))
            + x4 * (x373 * x4 + x375 + x4 * (x214 * x374 + x373))
            + x45 * (x349 * x371 + x372)
        )
    )
    result[0, 13] = numpy.sum(
        x217
        * (
            x0 * (x215 * x319 + x376 * x45 - x61 * (x339 * x85 + x377))
            + x4 * (x380 * x4 + x381 + x4 * (x215 * x331 + x380))
            + x45 * (x349 * x379 + x378)
        )
    )
    result[0, 14] = numpy.sum(
        x124
        * (
            x0 * (x383 + x385 * x4 - x61 * (x214 * x387 + x386))
            + x4 * (x392 + x393 + x4 * (x389 + x391))
            + x45 * (x349 * x387 + x388)
        )
    )
    result[0, 15] = numpy.sum(x4 * x83 * (x0 * (x394 - x396) + x398 + x399))
    result[0, 16] = numpy.sum(x124 * x4 * (x0 * (x400 - x402) + x349 * x403 + x404))
    result[0, 17] = numpy.sum(x199 * x4 * (x0 * (x405 - x407) + x349 * x408 + x409))
    result[0, 18] = numpy.sum(x280 * x4 * (x0 * (x410 - x412) + x349 * x413 + x414))
    result[0, 19] = numpy.sum(x124 * x4 * (x0 * (x384 * x85 - x416) + x349 * x418 + x417))
    result[0, 20] = numpy.sum(x4 * x83 * (x0 * (x419 - x421) + x423 + x424))
    result[1, 0] = numpy.sum(
        0.5 * x429 * (2.0 * x4 * (x426 + x427) + 2.0 * x425 * x78 * (x74 + x76) + x428)
    )
    result[1, 1] = numpy.sum(
        0.5
        * x436
        * (
            2.0 * x4 * (x4 * x431 + x433 * x81)
            + x434
            + 2.0 * x81 * (x123 * x435 + x4 * x433)
        )
    )
    result[1, 2] = numpy.sum(
        0.5
        * x436
        * (2.0 * x4 * x425 * (x151 + x153) + 2.0 * x425 * x81 * (x156 + x158) + x437)
    )
    result[1, 3] = numpy.sum(
        0.5
        * x301
        * (2.0 * x123 * (x194 * x446 + x4 * x443) + 2.0 * x4 * (x4 * x440 + x444) + x445)
    )
    result[1, 4] = numpy.sum(
        0.5
        * x453
        * (2.0 * x123 * (x4 * x448 + x452) + 2.0 * x4 * (x123 * x448 + x4 * x447) + x449)
    )
    result[1, 5] = numpy.sum(
        0.5
        * x301
        * (2.0 * x123 * x425 * (x254 + x255) + 2.0 * x4 * x425 * (x249 + x251) + x454)
    )
    result[1, 6] = numpy.sum(
        0.5
        * x462
        * (2.0 * x194 * (x4 * x459 + x461) + 2.0 * x4 * (x194 * x459 + x4 * x456) + x460)
    )
    result[1, 7] = numpy.sum(
        0.5
        * x471
        * (2.0 * x194 * (x4 * x465 + x470) + 2.0 * x4 * (x194 * x465 + x4 * x464) + x468)
    )
    result[1, 8] = numpy.sum(
        0.5
        * x471
        * (2.0 * x194 * (x4 * x474 + x477) + 2.0 * x4 * (x4 * x472 + x475) + x476)
    )
    result[1, 9] = numpy.sum(
        0.5
        * x462
        * (
            2.0 * x194 * (x339 * x425 + x480)
            + 2.0 * x4 * (x194 * x478 + x335 * x425)
            + x479
        )
    )
    result[1, 10] = numpy.sum(
        0.5 * x436 * (2.0 * x4**2 * x486 + 4.0 * x4 * x484 + x487)
    )
    result[1, 11] = numpy.sum(
        0.5 * x453 * (2.0 * x4**2 * x492 + 4.0 * x4 * x491 + x493)
    )
    result[1, 12] = numpy.sum(
        0.5 * x471 * (2.0 * x4**2 * x497 + 4.0 * x4 * x496 + x498)
    )
    result[1, 13] = numpy.sum(
        0.5 * x453 * (2.0 * x4**2 * x501 + 4.0 * x4 * x500 + x502)
    )
    result[1, 14] = numpy.sum(
        0.5 * x436 * (2.0 * x393 * x425 + 2.0 * x4 * (x391 * x425 + x503) + x504)
    )
    result[1, 15] = numpy.sum(0.5 * x429 * (2.0 * x349 * x507 + x508))
    result[1, 16] = numpy.sum(0.5 * x436 * (2.0 * x349 * x509 + x510))
    result[1, 17] = numpy.sum(0.5 * x301 * (2.0 * x349 * x512 + x515))
    result[1, 18] = numpy.sum(0.5 * x462 * (2.0 * x349 * x516 + x517))
    result[1, 19] = numpy.sum(0.5 * x436 * (2.0 * x349 * x518 + x520))
    result[1, 20] = numpy.sum(0.5 * x429 * (2.0 * x423 * x425 + x521))
    result[2, 0] = numpy.sum(
        x429 * (x4 * x522 * (x68 + x70) + x522 * x78 * (x74 + x76) + x524)
    )
    result[2, 1] = numpy.sum(
        x436 * (x4 * x522 * (x111 + x113) + x522 * x81 * (x116 + x118) + x526)
    )
    result[2, 2] = numpy.sum(
        x436 * (x4 * (x528 + x531) + x533 + x81 * (x4 * x529 + 3.0 * x534))
    )
    result[2, 3] = numpy.sum(
        x301 * (x123 * x522 * (x193 + x195) + x4 * (x188 * x522 + x535) + x537)
    )
    result[2, 4] = numpy.sum(
        x453 * (x123 * (x4 * x539 + x544) + x4 * (x123 * x539 + x4 * x538) + x541)
    )
    result[2, 5] = numpy.sum(
        x301 * (x123 * (x4 * x547 + 2.0 * x552) + x4 * (x546 + x549) + x551)
    )
    result[2, 6] = numpy.sum(
        x462 * (x194 * (x276 * x522 + x556) + x4 * (x194 * x553 + x272 * x522) + x555)
    )
    result[2, 7] = numpy.sum(
        x471 * (x194 * (x4 * x559 + x563) + x4 * (x4 * x557 + x560) + x562)
    )
    result[2, 8] = numpy.sum(
        x471 * (x194 * (x4 * x565 + x569) + x4 * (x4 * x564 + x566) + x568)
    )
    result[2, 9] = numpy.sum(
        x462 * (x194 * (x4 * x573 + x579) + x4 * (x571 + x575) + x577)
    )
    result[2, 10] = numpy.sum(x436 * (x355 * x522 + x4 * (x353 * x522 + x580) + x582))
    result[2, 11] = numpy.sum(x453 * (x4 * x585 + x4 * (x4 * x586 + x585) + x588))
    result[2, 12] = numpy.sum(x471 * (x4 * x591 + x4 * (x4 * x592 + x591) + x594))
    result[2, 13] = numpy.sum(x453 * (x4 * x596 + x4 * (x4 * x597 + x596) + x599))
    result[2, 14] = numpy.sum(x436 * (x4 * x602 + x4 * (x602 + x604) + x606))
    result[2, 15] = numpy.sum(x429 * (x398 * x522 + x608))
    result[2, 16] = numpy.sum(x436 * (x349 * x609 + x611))
    result[2, 17] = numpy.sum(x301 * (x349 * x612 + x614))
    result[2, 18] = numpy.sum(x462 * (x349 * x615 + x617))
    result[2, 19] = numpy.sum(x436 * (x349 * x618 + x620))
    result[2, 20] = numpy.sum(x429 * (x349 * x621 + x623))
    result[3, 0] = numpy.sum(x429 * (x4 * x626 + 5.0 * x627))
    result[3, 1] = numpy.sum(x436 * (x4 * x629 + x630 * x81))
    result[3, 2] = numpy.sum(x436 * (x4 * x631 + 4.0 * x632))
    result[3, 3] = numpy.sum(x301 * (x4 * x633 + x634))
    result[3, 4] = numpy.sum(x453 * (x123 * x636 + x4 * x635))
    result[3, 5] = numpy.sum(x301 * (x4 * x637 + 3.0 * x638))
    result[3, 6] = numpy.sum(x462 * (x194 * x640 + x4 * x639))
    result[3, 7] = numpy.sum(x471 * (x194 * x642 + x4 * x641))
    result[3, 8] = numpy.sum(x471 * (x4 * x643 + x644))
    result[3, 9] = numpy.sum(x462 * (x4 * x645 + 2.0 * x646))
    result[3, 10] = numpy.sum(x436 * (x4 * x648 + x647))
    result[3, 11] = numpy.sum(x453 * (x4 * x650 + x649))
    result[3, 12] = numpy.sum(x471 * (x4 * x652 + x651))
    result[3, 13] = numpy.sum(x453 * (x4 * x654 + x653))
    result[3, 14] = numpy.sum(x436 * (x4 * x656 + x655))
    result[3, 15] = numpy.sum(x657 * x658)
    result[3, 16] = numpy.sum(x659 * x660)
    result[3, 17] = numpy.sum(x661 * x662)
    result[3, 18] = numpy.sum(x663 * x664)
    result[3, 19] = numpy.sum(x660 * x665)
    result[3, 20] = numpy.sum(x658 * x666)
    result[4, 0] = numpy.sum(x522 * x667 * (x426 + x427))
    result[4, 1] = numpy.sum(x672 * (x4 * x669 + x671 * x81))
    result[4, 2] = numpy.sum(x425 * x672 * (x528 + x531))
    result[4, 3] = numpy.sum(x453 * (x4 * x673 + x674))
    result[4, 4] = numpy.sum(x677 * (x123 * x676 + x4 * x675))
    result[4, 5] = numpy.sum(x425 * x453 * (x546 + x549))
    result[4, 6] = numpy.sum(x680 * (x194 * x679 + x4 * x678))
    result[4, 7] = numpy.sum(x684 * (x194 * x683 + x4 * x682))
    result[4, 8] = numpy.sum(x684 * (x4 * x685 + x686))
    result[4, 9] = numpy.sum(x425 * x680 * (x571 + x575))
    result[4, 10] = numpy.sum(x672 * (x4 * x688 + x687))
    result[4, 11] = numpy.sum(x677 * (x4 * x690 + x689))
    result[4, 12] = numpy.sum(x684 * (x4 * x692 + x691))
    result[4, 13] = numpy.sum(x677 * (x4 * x694 + x693))
    result[4, 14] = numpy.sum(x672 * (x425 * x604 + x695))
    result[4, 15] = numpy.sum(x696 * x697)
    result[4, 16] = numpy.sum(x698 * x699)
    result[4, 17] = numpy.sum(x4 * x453 * x701)
    result[4, 18] = numpy.sum(x4 * x680 * x702)
    result[4, 19] = numpy.sum(x699 * x703)
    result[4, 20] = numpy.sum(x425 * x621 * x697)
    result[5, 0] = numpy.sum(x429 * (x4 * x705 + 5.0 * x706))
    result[5, 1] = numpy.sum(x436 * (x4 * x707 + 4.0 * x708))
    result[5, 2] = numpy.sum(x436 * (x4 * x709 + 4.0 * x710))
    result[5, 3] = numpy.sum(x301 * (x4 * x711 + x713))
    result[5, 4] = numpy.sum(x453 * (x123 * x715 + x4 * x714))
    result[5, 5] = numpy.sum(x301 * (x4 * x716 + x718))
    result[5, 6] = numpy.sum(x462 * (x4 * x719 + 2.0 * x720))
    result[5, 7] = numpy.sum(x471 * (x4 * x721 + x723))
    result[5, 8] = numpy.sum(x471 * (x4 * x724 + x726))
    result[5, 9] = numpy.sum(x462 * (x4 * x727 + 2.0 * x728))
    result[5, 10] = numpy.sum(x436 * (x4 * x730 + x729))
    result[5, 11] = numpy.sum(x453 * (x4 * x732 + x731))
    result[5, 12] = numpy.sum(x471 * (x4 * x734 + x733))
    result[5, 13] = numpy.sum(x453 * (x4 * x736 + x735))
    result[5, 14] = numpy.sum(x436 * (x4 * x738 + x737))
    result[5, 15] = numpy.sum(x658 * x739)
    result[5, 16] = numpy.sum(x660 * x740)
    result[5, 17] = numpy.sum(x662 * x741)
    result[5, 18] = numpy.sum(x664 * x742)
    result[5, 19] = numpy.sum(x660 * x743)
    result[5, 20] = numpy.sum(x658 * x744)
    result[6, 0] = numpy.sum(x83 * (x425 * x626 + x428))
    result[6, 1] = numpy.sum(x124 * (x425 * x629 + x434 + x627))
    result[6, 2] = numpy.sum(x124 * (x425 * x631 + x437))
    result[6, 3] = numpy.sum(x199 * (x194 * x630 + x425 * x633 + x445))
    result[6, 4] = numpy.sum(x217 * (x425 * x635 + x449 + x632))
    result[6, 5] = numpy.sum(x199 * (x425 * x637 + x454))
    result[6, 6] = numpy.sum(x280 * (x425 * x639 + x460 + x634))
    result[6, 7] = numpy.sum(x301 * (x194 * x636 + x425 * x641 + x468))
    result[6, 8] = numpy.sum(x301 * (x425 * x643 + x476 + x638))
    result[6, 9] = numpy.sum(x280 * (x425 * x645 + x479))
    result[6, 10] = numpy.sum(x124 * (x425 * x648 + x487 + x640 * x81))
    result[6, 11] = numpy.sum(x217 * (x123 * x642 + x425 * x650 + x493))
    result[6, 12] = numpy.sum(x301 * (x425 * x652 + x498 + x644))
    result[6, 13] = numpy.sum(x217 * (x425 * x654 + x502 + x646))
    result[6, 14] = numpy.sum(x124 * (x425 * x656 + x504))
    result[6, 15] = numpy.sum(x83 * (x425 * x657 + x508 + 5.0 * x647))
    result[6, 16] = numpy.sum(x124 * (x425 * x659 + x510 + 4.0 * x649))
    result[6, 17] = numpy.sum(x199 * (x425 * x661 + x515 + 3.0 * x651))
    result[6, 18] = numpy.sum(x280 * (x425 * x663 + x517 + 2.0 * x653))
    result[6, 19] = numpy.sum(x124 * (x425 * x665 + x520 + x655))
    result[6, 20] = numpy.sum(x83 * (x425 * x666 + x521))
    result[7, 0] = numpy.sum(x429 * (x522 * x625 + x524))
    result[7, 1] = numpy.sum(x436 * (x425 * x669 + x522 * x628 + x526))
    result[7, 2] = numpy.sum(x436 * (x527 * x624 + x533))
    result[7, 3] = numpy.sum(x301 * (x194 * x671 + x425 * x673 + x537))
    result[7, 4] = numpy.sum(x453 * (x425 * x530 + x425 * x675 + x541))
    result[7, 5] = numpy.sum(x301 * (x545 * x624 + x551))
    result[7, 6] = numpy.sum(x462 * (x425 * x678 + x555 + x674))
    result[7, 7] = numpy.sum(x471 * (x194 * x676 + x425 * x682 + x562))
    result[7, 8] = numpy.sum(x471 * (x425 * x548 + x425 * x685 + x568))
    result[7, 9] = numpy.sum(x462 * (x570 * x624 + x577))
    result[7, 10] = numpy.sum(x436 * (x425 * x688 + x582 + x679 * x81))
    result[7, 11] = numpy.sum(x453 * (x123 * x683 + x425 * x690 + x588))
    result[7, 12] = numpy.sum(x471 * (x425 * x692 + x594 + x686))
    result[7, 13] = numpy.sum(x453 * (x425 * x574 + x425 * x694 + x599))
    result[7, 14] = numpy.sum(x436 * (x603 * x624 + x606))
    result[7, 15] = numpy.sum(x429 * (x425 * x696 + x608 + 5.0 * x687))
    result[7, 16] = numpy.sum(x436 * (x425 * x698 + x611 + 4.0 * x689))
    result[7, 17] = numpy.sum(x301 * (x425 * x701 + x614 + 3.0 * x691))
    result[7, 18] = numpy.sum(x462 * (x425 * x702 + x617 + 2.0 * x693))
    result[7, 19] = numpy.sum(x436 * (x425 * x703 + x620 + x695))
    result[7, 20] = numpy.sum(x429 * (x621 * x624 + x623))
    result[8, 0] = numpy.sum(x705 * x745)
    result[8, 1] = numpy.sum(x436 * (x425 * x707 + x706))
    result[8, 2] = numpy.sum(x709 * x746)
    result[8, 3] = numpy.sum(x301 * (x425 * x711 + 2.0 * x708))
    result[8, 4] = numpy.sum(x453 * (x425 * x714 + x710))
    result[8, 5] = numpy.sum(x301 * x425 * x716)
    result[8, 6] = numpy.sum(x462 * (x425 * x719 + x713))
    result[8, 7] = numpy.sum(x471 * (x425 * x721 + x747))
    result[8, 8] = numpy.sum(x471 * (x425 * x724 + x717))
    result[8, 9] = numpy.sum(x425 * x462 * x727)
    result[8, 10] = numpy.sum(x436 * (x425 * x730 + 4.0 * x720))
    result[8, 11] = numpy.sum(x453 * (x123 * x722 + x425 * x732))
    result[8, 12] = numpy.sum(x471 * (x425 * x734 + x726))
    result[8, 13] = numpy.sum(x453 * (x425 * x736 + x728))
    result[8, 14] = numpy.sum(x738 * x746)
    result[8, 15] = numpy.sum(x429 * (x425 * x739 + 5.0 * x729))
    result[8, 16] = numpy.sum(x436 * (x425 * x740 + 4.0 * x731))
    result[8, 17] = numpy.sum(x301 * (x425 * x741 + x748))
    result[8, 18] = numpy.sum(x462 * (x425 * x742 + 2.0 * x735))
    result[8, 19] = numpy.sum(x436 * (x425 * x743 + x737))
    result[8, 20] = numpy.sum(x744 * x745)
    result[9, 0] = numpy.sum(x83 * (x522 * x705 + x523))
    result[9, 1] = numpy.sum(x124 * (x522 * x707 + x525))
    result[9, 2] = numpy.sum(x124 * (x522 * x709 + x532 + x706))
    result[9, 3] = numpy.sum(x199 * (x522 * x711 + x536))
    result[9, 4] = numpy.sum(x217 * (x522 * x714 + x540 + x708))
    result[9, 5] = numpy.sum(x199 * (x522 * x716 + x550 + 2.0 * x710))
    result[9, 6] = numpy.sum(x280 * (x522 * x719 + x554))
    result[9, 7] = numpy.sum(x301 * (x522 * x721 + x561 + x712))
    result[9, 8] = numpy.sum(x301 * (x522 * x724 + x567 + x747))
    result[9, 9] = numpy.sum(x280 * (x522 * x727 + x576 + x718))
    result[9, 10] = numpy.sum(x124 * (x522 * x730 + x581))
    result[9, 11] = numpy.sum(x217 * (x522 * x732 + x587 + x720))
    result[9, 12] = numpy.sum(x301 * (x522 * x734 + x593 + x723))
    result[9, 13] = numpy.sum(x217 * (x123 * x725 + x522 * x736 + x598))
    result[9, 14] = numpy.sum(x124 * (x522 * x738 + x605 + 4.0 * x728))
    result[9, 15] = numpy.sum(x83 * (x522 * x739 + x607))
    result[9, 16] = numpy.sum(x124 * (x522 * x740 + x610 + x729))
    result[9, 17] = numpy.sum(x199 * (x522 * x741 + x613 + 2.0 * x731))
    result[9, 18] = numpy.sum(x280 * (x522 * x742 + x616 + x748))
    result[9, 19] = numpy.sum(x124 * (x522 * x743 + x619 + 4.0 * x735))
    result[9, 20] = numpy.sum(x83 * (x522 * x744 + x622 + 5.0 * x737))
    return result


def _2center2el3d_40(ax, da, A, bx, db, B):
    """Cartesian (g|s) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 1), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = x1 * (ax * A[0] + bx * B[0]) - A[0]
    x3 = ax ** (-1.0)
    x4 = bx * x1
    x5 = ax * x4 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x6 = boys(3, x5)
    x7 = 17.4934183276249
    x8 = 2.0 * x3 * x7
    x9 = x0 ** (-1.5) * x8
    x10 = x6 * x9
    x11 = x10 * x2
    x12 = bx ** (-1.0)
    x13 = x0 ** (-0.5)
    x14 = boys(2, x5)
    x15 = 0.5 * x3
    x16 = x15 * (-x10 + 2.0 * x12 * x13 * x14 * x3 * x7)
    x17 = boys(4, x5)
    x18 = x2**2
    x19 = x12 * x13 * x8
    x20 = x18 * x19
    x21 = x17 * x20
    x22 = boys(1, x5)
    x23 = x15 * (2.0 * x12 * x13 * x3 * x7 * boys(0, x5) - x22 * x9)
    x24 = x15 * (2.0 * x12 * x13 * x22 * x3 * x7 - x14 * x9)
    x25 = x14 * x19
    x26 = 1.5 * x3
    x27 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**1.5)
    x28 = 4.41641957979107 * x27
    x29 = x1 * (ax * A[1] + bx * B[1]) - A[1]
    x30 = x10 * x29
    x31 = x3 * (2.0 * x12 * x13 * x14 * x29 * x3 * x7 - x30)
    x32 = x21 * x29
    x33 = 11.6847478934435 * x27
    x34 = x1 * (ax * A[2] + bx * B[2]) - A[2]
    x35 = x3 * x34 * (-x10 + 2.0 * x12 * x13 * x14 * x3 * x7)
    x36 = 0.5 * x35
    x37 = x29**2
    x38 = x19 * x37
    x39 = x17 * x38
    x40 = x16 + x39
    x41 = x23 + x25 * x37 - x4 * (x24 + x38 * x6)
    x42 = 15.084944665313 * x27
    x43 = x3 * x34 * (2.0 * x12 * x13 * x14 * x29 * x3 * x7 - x30)
    x44 = 26.1278905896872 * x27
    x45 = x34**2
    x46 = x19 * x45
    x47 = x16 + x17 * x46
    x48 = x23 + x25 * x45 - x4 * (x24 + x46 * x6)
    x49 = x15 * x48
    x50 = x29 * x40 + x31
    x51 = x2 * x33
    x52 = x34 * x39 + x36
    x53 = x2 * x44
    x54 = x34 * x47 + x35

    # 15 item(s)
    result[0, 0] = numpy.sum(
        x28
        * (
            x2 * (x2 * (x16 + x21) - x3 * (x11 - 2.0 * x12 * x13 * x14 * x2 * x3 * x7))
            + x26 * (x18 * x25 + x23 - x4 * (x20 * x6 + x24))
        )
    )
    result[1, 0] = numpy.sum(
        0.5
        * x33
        * (
            x2 * (x31 + 2.0 * x32)
            - 2.0 * x29 * x3 * (x11 - 2.0 * x12 * x13 * x14 * x2 * x3 * x7)
        )
    )
    result[2, 0] = numpy.sum(
        x33
        * (
            x2 * (x21 * x34 + x36)
            - x3 * x34 * (x11 - 2.0 * x12 * x13 * x14 * x2 * x3 * x7)
        )
    )
    result[3, 0] = numpy.sum(x42 * (x15 * x41 + x18 * x40))
    result[4, 0] = numpy.sum(0.5 * x44 * (2.0 * x32 * x34 + x43))
    result[5, 0] = numpy.sum(x42 * (x18 * x47 + x49))
    result[6, 0] = numpy.sum(x50 * x51)
    result[7, 0] = numpy.sum(x52 * x53)
    result[8, 0] = numpy.sum(x29 * x47 * x53)
    result[9, 0] = numpy.sum(x51 * x54)
    result[10, 0] = numpy.sum(x28 * (x26 * x41 + x29 * x50))
    result[11, 0] = numpy.sum(x33 * (x29 * x52 + x43))
    result[12, 0] = numpy.sum(x42 * (x37 * x47 + x49))
    result[13, 0] = numpy.sum(x29 * x33 * x54)
    result[14, 0] = numpy.sum(x28 * (x26 * x48 + x34 * x54))
    return result


def _2center2el3d_41(ax, da, A, bx, db, B):
    """Cartesian (g|p) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 3), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = ax ** (-1.0)
    x2 = ax + bx
    x3 = x2 ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - A[0]
    x6 = bx * x3
    x7 = ax * x6 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x8 = boys(3, x7)
    x9 = 17.4934183276249
    x10 = 2.0 * x1 * x9
    x11 = x10 * x2 ** (-1.5)
    x12 = x11 * x8
    x13 = bx ** (-1.0)
    x14 = x2 ** (-0.5)
    x15 = boys(2, x7)
    x16 = 0.5 * x1
    x17 = x16 * (2.0 * x1 * x13 * x14 * x15 * x9 - x12)
    x18 = x5**2
    x19 = boys(4, x7)
    x20 = x10 * x13 * x14
    x21 = x19 * x20
    x22 = x18 * x21
    x23 = x17 + x22
    x24 = x20 * x8
    x25 = x0 * x24
    x26 = -x4 - B[0]
    x27 = x26 * x5
    x28 = x21 * x27
    x29 = x25 + x28
    x30 = x15 * x20
    x31 = x0 * x30
    x32 = x24 * x26
    x33 = x32 * x5
    x34 = x31 + x33
    x35 = x0 * x21
    x36 = x20 * boys(5, x7)
    x37 = x27 * x36
    x38 = x11 * x26
    x39 = x19 * x38
    x40 = x16 * (2.0 * x1 * x13 * x14 * x26 * x8 * x9 - x39)
    x41 = x35 * x5
    x42 = boys(1, x7)
    x43 = x16 * (2.0 * x1 * x13 * x14 * x26 * x42 * x9 - x15 * x38)
    x44 = x16 * x26 * (2.0 * x1 * x13 * x14 * x15 * x9 - x12)
    x45 = 1.5 * x1
    x46 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**2.5)
    x47 = 8.83283915958214 * x46
    x48 = -x3 * (ax * A[1] + bx * B[1])
    x49 = -x48 - B[1]
    x50 = x11 * x49
    x51 = x19 * x50
    x52 = x5 * x51
    x53 = x16 * (2.0 * x1 * x13 * x14 * x49 * x8 * x9 - x51)
    x54 = x18 * x36
    x55 = x49 * x54
    x56 = x16 * (2.0 * x1 * x13 * x14 * x42 * x49 * x9 - x15 * x50)
    x57 = x16 * x49 * (2.0 * x1 * x13 * x14 * x15 * x9 - x12)
    x58 = x24 * x49
    x59 = -x3 * (ax * A[2] + bx * B[2])
    x60 = -x59 - B[2]
    x61 = x11 * x60
    x62 = x19 * x61
    x63 = x5 * x62
    x64 = x16 * (2.0 * x1 * x13 * x14 * x60 * x8 * x9 - x62)
    x65 = x54 * x60
    x66 = x16 * (2.0 * x1 * x13 * x14 * x42 * x60 * x9 - x15 * x61)
    x67 = x16 * x60 * (2.0 * x1 * x13 * x14 * x15 * x9 - x12)
    x68 = x24 * x60
    x69 = -x48 - A[1]
    x70 = x25 * x69
    x71 = x31 * x69
    x72 = x1 * x69 * (2.0 * x1 * x13 * x14 * x15 * x9 - x12)
    x73 = x35 * x69
    x74 = x37 * x69
    x75 = x39 * x69
    x76 = x1 * (2.0 * x1 * x13 * x14 * x26 * x69 * x8 * x9 - x75)
    x77 = x41 * x69
    x78 = 23.3694957868871 * x46
    x79 = x58 * x69
    x80 = x31 + x79
    x81 = x49 * x69
    x82 = x21 * x81
    x83 = x25 + x82
    x84 = x6 * x83
    x85 = x36 * x81
    x86 = x35 + x85
    x87 = x1 * (x80 - x84)
    x88 = x1 * x69 * (2.0 * x1 * x13 * x14 * x60 * x8 * x9 - x62)
    x89 = -x59 - A[2]
    x90 = x25 * x89
    x91 = x31 * x89
    x92 = x1 * x89 * (2.0 * x1 * x13 * x14 * x15 * x9 - x12)
    x93 = 0.5 * x92
    x94 = x35 * x89
    x95 = x1 * x89 * (2.0 * x1 * x13 * x14 * x26 * x8 * x9 - x39)
    x96 = 0.5 * x95
    x97 = x1 * x89 * (2.0 * x1 * x13 * x14 * x49 * x8 * x9 - x51)
    x98 = 0.5 * x97
    x99 = x31 + x68 * x89
    x100 = x60 * x89
    x101 = x100 * x21 + x25
    x102 = x101 * x6
    x103 = x100 * x36 + x35
    x104 = x103 * x18
    x105 = x1 * (-x102 + x99)
    x106 = 0.5 * x105
    x107 = x69**2
    x108 = x107 * x21
    x109 = x108 + x17
    x110 = x0 * x109
    x111 = x107 * x36
    x112 = x111 * x26
    x113 = x112 + x40
    x114 = x107 * x32 + x43 - x6 * (x108 * x26 + x44)
    x115 = 30.169889330626 * x46
    x116 = x53 + x69 * x86 + x73
    x117 = x56 - x6 * (x57 + x69 * x83 + x70) + x69 * x80 + x71
    x118 = x111 * x60 + x64
    x119 = x107 * x68 - x6 * (x108 * x60 + x67) + x66
    x120 = x73 * x89
    x121 = x1 * x89 * (2.0 * x1 * x13 * x14 * x26 * x69 * x8 * x9 - x75)
    x122 = 52.2557811793745 * x46
    x123 = x85 * x89 + x94
    x124 = x1 * (-x6 * (x82 * x89 + x90) + x79 * x89 + x91)
    x125 = x1 * x69 * (-x102 + x99)
    x126 = x89**2
    x127 = x126 * x21
    x128 = x127 + x17
    x129 = x0 * x128
    x130 = x126 * x36
    x131 = x130 * x26 + x40
    x132 = x131 * x5
    x133 = x126 * x32 + x43 - x6 * (x127 * x26 + x44)
    x134 = x133 * x16
    x135 = x130 * x49 + x53
    x136 = x126 * x58 + x56 - x6 * (x127 * x49 + x57)
    x137 = x136 * x16
    x138 = x103 * x89 + x64 + x94
    x139 = -x6 * (x101 * x89 + x67 + x90) + x66 + x89 * x99 + x91
    x140 = x139 * x16
    x141 = x0 * (x109 * x69 + x72)
    x142 = x113 * x69 + x76
    x143 = x110 + x116 * x69 + x87
    x144 = x5 * x78
    x145 = x118 * x69 + x88
    x146 = x0 * (x108 * x89 + x93)
    x147 = x112 * x89 + x96
    x148 = x120 + x123 * x69 + x98
    x149 = x122 * x5
    x150 = x103 * x107 + x106
    x151 = x129 * x69
    x152 = x129 + x135 * x69
    x153 = x0 * (x128 * x89 + x92)
    x154 = x131 * x89 + x95
    x155 = x135 * x89 + x97
    x156 = x105 + x129 + x138 * x89
    x157 = x69 * x78

    # 45 item(s)
    result[0, 0] = numpy.sum(
        x47
        * (
            x0 * x5 * (x1 * (2.0 * x1 * x13 * x14 * x15 * x9 - x12) + x23)
            + x45 * (x0 * x30 * x5 + x34 * x5 + x43 - x6 * (x25 * x5 + x29 * x5 + x44))
            + x5
            * (x0 * x23 - x1 * (x29 * x6 - x34) + x5 * (x40 + x41 + x5 * (x35 + x37)))
        )
    )
    result[0, 1] = numpy.sum(
        x47
        * (
            x45 * (x18 * x58 + x56 - x6 * (x22 * x49 + x57))
            + x5
            * (x1 * (2.0 * x1 * x13 * x14 * x49 * x5 * x8 * x9 - x52) + x5 * (x53 + x55))
        )
    )
    result[0, 2] = numpy.sum(
        x47
        * (
            x45 * (x18 * x68 - x6 * (x22 * x60 + x67) + x66)
            + x5
            * (x1 * (2.0 * x1 * x13 * x14 * x5 * x60 * x8 * x9 - x63) + x5 * (x64 + x65))
        )
    )
    result[1, 0] = numpy.sum(
        0.5
        * x78
        * (
            x0 * (2.0 * x22 * x69 + x72)
            + 2.0 * x1 * (x33 * x69 - x6 * (x28 * x69 + x70) + x71)
            + x5 * (2.0 * x5 * (x73 + x74) + x76 + 2.0 * x77)
        )
    )
    result[1, 1] = numpy.sum(
        0.5 * x5 * x78 * (2.0 * x1 * (x80 - x84) + 2.0 * x18 * x86 + x87)
    )
    result[1, 2] = numpy.sum(
        0.5
        * x78
        * (
            2.0 * x1 * x69 * (2.0 * x1 * x13 * x14 * x5 * x60 * x8 * x9 - x63)
            + x5 * (2.0 * x65 * x69 + x88)
        )
    )
    result[2, 0] = numpy.sum(
        x78
        * (
            x0 * (x22 * x89 + x93)
            + x1 * (x33 * x89 - x6 * (x28 * x89 + x90) + x91)
            + x5 * (x41 * x89 + x5 * (x37 * x89 + x94) + x96)
        )
    )
    result[2, 1] = numpy.sum(
        x78
        * (
            x1 * x89 * (2.0 * x1 * x13 * x14 * x49 * x5 * x8 * x9 - x52)
            + x5 * (x55 * x89 + x98)
        )
    )
    result[2, 2] = numpy.sum(x5 * x78 * (-x1 * (x102 - x99) + x104 + x106))
    result[3, 0] = numpy.sum(x115 * (x110 * x5 + x114 * x16 + x5 * (x110 + x113 * x5)))
    result[3, 1] = numpy.sum(x115 * (x116 * x18 + x117 * x16))
    result[3, 2] = numpy.sum(x115 * (x118 * x18 + x119 * x16))
    result[4, 0] = numpy.sum(
        0.5 * x122 * (x121 + 2.0 * x5 * (x120 + x74 * x89) + 2.0 * x77 * x89)
    )
    result[4, 1] = numpy.sum(0.5 * x122 * (2.0 * x123 * x18 + x124))
    result[4, 2] = numpy.sum(0.5 * x122 * (2.0 * x104 * x69 + x125))
    result[5, 0] = numpy.sum(x115 * (x129 * x5 + x134 + x5 * (x129 + x132)))
    result[5, 1] = numpy.sum(x115 * (x135 * x18 + x137))
    result[5, 2] = numpy.sum(x115 * (x138 * x18 + x140))
    result[6, 0] = numpy.sum(x78 * (x141 + x142 * x5))
    result[6, 1] = numpy.sum(x143 * x144)
    result[6, 2] = numpy.sum(x144 * x145)
    result[7, 0] = numpy.sum(x122 * (x146 + x147 * x5))
    result[7, 1] = numpy.sum(x148 * x149)
    result[7, 2] = numpy.sum(x149 * x150)
    result[8, 0] = numpy.sum(x122 * (x132 * x69 + x151))
    result[8, 1] = numpy.sum(x149 * x152)
    result[8, 2] = numpy.sum(x138 * x149 * x69)
    result[9, 0] = numpy.sum(x78 * (x153 + x154 * x5))
    result[9, 1] = numpy.sum(x144 * x155)
    result[9, 2] = numpy.sum(x144 * x156)
    result[10, 0] = numpy.sum(x47 * (x114 * x45 + x142 * x69))
    result[10, 1] = numpy.sum(x47 * (x117 * x45 + x141 + x143 * x69))
    result[10, 2] = numpy.sum(x47 * (x119 * x45 + x145 * x69))
    result[11, 0] = numpy.sum(x78 * (x121 + x147 * x69))
    result[11, 1] = numpy.sum(x78 * (x124 + x146 + x148 * x69))
    result[11, 2] = numpy.sum(x78 * (x125 + x150 * x69))
    result[12, 0] = numpy.sum(x115 * (x107 * x131 + x134))
    result[12, 1] = numpy.sum(x115 * (x137 + x151 + x152 * x69))
    result[12, 2] = numpy.sum(x115 * (x107 * x138 + x140))
    result[13, 0] = numpy.sum(x154 * x157)
    result[13, 1] = numpy.sum(x78 * (x153 + x155 * x69))
    result[13, 2] = numpy.sum(x156 * x157)
    result[14, 0] = numpy.sum(x47 * (x133 * x45 + x154 * x89))
    result[14, 1] = numpy.sum(x47 * (x136 * x45 + x155 * x89))
    result[14, 2] = numpy.sum(x47 * (x139 * x45 + x153 + x156 * x89))
    return result


def _2center2el3d_42(ax, da, A, bx, db, B):
    """Cartesian (g|d) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 6), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = ax ** (-1.0)
    x5 = bx * x1
    x6 = ax * x5 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x7 = boys(4, x6)
    x8 = x0 ** (-1.5)
    x9 = bx ** (-1.0)
    x10 = 17.4934183276249
    x11 = x10 * x9
    x12 = 2.0 * x11
    x13 = x12 * x8
    x14 = x0 ** (-0.5)
    x15 = boys(3, x6)
    x16 = 0.5 * x9
    x17 = x16 * (2.0 * x10 * x14 * x15 * x4 * x9 - x13 * x7)
    x18 = boys(5, x6)
    x19 = -x2 - B[0]
    x20 = x19**2
    x21 = x12 * x14
    x22 = x21 * x4
    x23 = x20 * x22
    x24 = x17 + x18 * x23
    x25 = x24 * x3
    x26 = x4 * x7
    x27 = 0.5 / (ax + bx)
    x28 = x11 * x14
    x29 = 4.0 * x27 * x28
    x30 = x19 * x29
    x31 = x26 * x30
    x32 = x25 + x31
    x33 = boys(2, x6)
    x34 = -2.0 * x10 * x14 * x33 * x4 * x9
    x35 = -x16 * (x13 * x15 + x34)
    x36 = x22 * x7
    x37 = x20 * x36 + x35
    x38 = x3 * x37
    x39 = x15 * x4
    x40 = x30 * x39
    x41 = x38 + x40
    x42 = x16 * (2.0 * x10 * x14 * x4 * x7 * x9 - x13 * x18)
    x43 = boys(6, x6)
    x44 = x23 * x43 + x42
    x45 = x3 * x44
    x46 = x18 * x4
    x47 = x30 * x46
    x48 = x24 * x5
    x49 = 0.5 * x4
    x50 = x49 * (x37 - x48)
    x51 = 2.0 * x27
    x52 = x28 * x51
    x53 = x26 * x52
    x54 = x19 * x3
    x55 = x18 * x22
    x56 = x54 * x55
    x57 = x53 + x56
    x58 = 2.0 * x10 * x8
    x59 = x26 * x58
    x60 = x19 * x59
    x61 = x49 * (2.0 * x10 * x14 * x15 * x19 * x4 * x9 - x60)
    x62 = x3 * x53
    x63 = x3 * x57 + x61 + x62
    x64 = x39 * x52
    x65 = x36 * x54 + x64
    x66 = x4 * x52
    x67 = x33 * x66
    x68 = x21 * x39
    x69 = x19 * x3 * x68 + x67
    x70 = x39 * x58
    x71 = -x49 * (x34 + x70)
    x72 = x3**2
    x73 = x36 * x72
    x74 = boys(1, x6)
    x75 = x16 * (2.0 * x10 * x14 * x4 * x9 * boys(0, x6) - x13 * x74)
    x76 = x16 * (2.0 * x10 * x14 * x4 * x74 * x9 - x13 * x33)
    x77 = x20 * x68 + x76
    x78 = x22 * x33
    x79 = x49 * (x20 * x78 - x5 * x77 + x75)
    x80 = x49 * (-x37 * x5 + x77)
    x81 = 1.5 * x4
    x82 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**3.5)
    x83 = 10.1992841329868 * x82
    x84 = -x1 * (ax * A[1] + bx * B[1])
    x85 = -x84 - B[1]
    x86 = x59 * x85
    x87 = x49 * (2.0 * x10 * x14 * x15 * x4 * x85 * x9 - x86)
    x88 = x55 * x72
    x89 = x85 * x88
    x90 = x87 + x89
    x91 = x53 * x85
    x92 = x19 * x85
    x93 = x3 * x92
    x94 = x55 * x93
    x95 = x91 + x94
    x96 = x64 * x85
    x97 = x36 * x92
    x98 = x3 * x97
    x99 = x96 + x98
    x100 = x18 * x66
    x101 = x100 * x85
    x102 = x22 * x43
    x103 = x102 * x93
    x104 = x46 * x58
    x105 = x104 * x92
    x106 = x49 * (2.0 * x10 * x14 * x19 * x4 * x7 * x85 * x9 - x105)
    x107 = x101 * x3
    x108 = x49 * (2.0 * x10 * x14 * x19 * x33 * x4 * x85 * x9 - x70 * x92)
    x109 = x49 * x85 * (2.0 * x10 * x14 * x15 * x19 * x4 * x9 - x60)
    x110 = 17.6656783191643 * x82
    x111 = -x1 * (ax * A[2] + bx * B[2])
    x112 = -x111 - B[2]
    x113 = x112 * x59
    x114 = x49 * (2.0 * x10 * x112 * x14 * x15 * x4 * x9 - x113)
    x115 = x112 * x88
    x116 = x114 + x115
    x117 = x112 * x53
    x118 = x112 * x19
    x119 = x118 * x3
    x120 = x119 * x55
    x121 = x117 + x120
    x122 = x112 * x64
    x123 = x118 * x36
    x124 = x123 * x3
    x125 = x122 + x124
    x126 = x100 * x112
    x127 = x102 * x119
    x128 = x104 * x118
    x129 = x49 * (2.0 * x10 * x112 * x14 * x19 * x4 * x7 * x9 - x128)
    x130 = x126 * x3
    x131 = x49 * (2.0 * x10 * x112 * x14 * x19 * x33 * x4 * x9 - x118 * x70)
    x132 = x112 * x49 * (2.0 * x10 * x14 * x15 * x19 * x4 * x9 - x60)
    x133 = x85**2
    x134 = x133 * x36 + x35
    x135 = x133 * x22
    x136 = x135 * x18 + x17
    x137 = x136 * x5
    x138 = x137 * x3
    x139 = x135 * x43 + x42
    x140 = x139 * x72
    x141 = x49 * (x134 - x137)
    x142 = x133 * x68 + x76
    x143 = x49 * (x133 * x78 - x142 * x5 + x75)
    x144 = x49 * (-x134 * x5 + x142)
    x145 = x112 * x85
    x146 = x104 * x145
    x147 = x49 * (2.0 * x10 * x112 * x14 * x4 * x7 * x85 * x9 - x146)
    x148 = x102 * x145
    x149 = x49 * (2.0 * x10 * x112 * x14 * x33 * x4 * x85 * x9 - x145 * x70)
    x150 = x112 * x49 * (2.0 * x10 * x14 * x15 * x4 * x85 * x9 - x86)
    x151 = x112**2
    x152 = x151 * x36 + x35
    x153 = x151 * x22
    x154 = x153 * x18 + x17
    x155 = x154 * x5
    x156 = x155 * x3
    x157 = x153 * x43 + x42
    x158 = x157 * x72
    x159 = x49 * (x152 - x155)
    x160 = x151 * x68 + x76
    x161 = x49 * (x151 * x78 - x160 * x5 + x75)
    x162 = x49 * (-x152 * x5 + x160)
    x163 = -x84 - A[1]
    x164 = x163 * x45
    x165 = x163 * x47
    x166 = x163 * x48
    x167 = x4 * (x163 * x37 - x166)
    x168 = x163 * x53
    x169 = x163 * x56
    x170 = x168 + x169
    x171 = x163 * x4 * (2.0 * x10 * x14 * x15 * x19 * x4 * x9 - x60)
    x172 = 26.9847693667702 * x82
    x173 = x163 * x68 * x85 + x67
    x174 = x173 * x27
    x175 = x19 * x64
    x176 = x163 * x97
    x177 = x175 + x176
    x178 = x163 * x85
    x179 = x178 * x36 + x64
    x180 = x179 * x27
    x181 = x19 * x53
    x182 = x163 * x92
    x183 = x182 * x55
    x184 = x181 + x183
    x185 = x178 * x55
    x186 = x185 + x53
    x187 = x4 * (x173 - x179 * x5)
    x188 = x186 * x27
    x189 = x100 * x19
    x190 = x102 * x182
    x191 = x189 + x190
    x192 = x4 * (x177 - x184 * x5)
    x193 = 46.7389915737742 * x82
    x194 = x117 * x163
    x195 = x122 * x163
    x196 = x163 * x4 * (2.0 * x10 * x112 * x14 * x15 * x4 * x9 - x113)
    x197 = x126 * x163
    x198 = x163 * x4 * (2.0 * x10 * x112 * x14 * x19 * x4 * x7 * x9 - x128)
    x199 = x134 * x163
    x200 = x29 * x85
    x201 = x200 * x39
    x202 = x199 + x201
    x203 = x136 * x163
    x204 = x200 * x26
    x205 = x203 + x204
    x206 = x205 * x5
    x207 = x139 * x163
    x208 = x200 * x46
    x209 = x207 + x208
    x210 = x4 * (x202 - x206)
    x211 = x145 * x36
    x212 = x122 + x163 * x211
    x213 = x145 * x55
    x214 = x117 + x163 * x213
    x215 = x214 * x5
    x216 = x126 + x148 * x163
    x217 = x4 * (x212 - x215)
    x218 = x163 * x4 * (x152 - x155)
    x219 = -x111 - A[2]
    x220 = x219 * x4 * (x37 - x48)
    x221 = 0.5 * x220
    x222 = x219 * x53
    x223 = x219 * x56 + x222
    x224 = x219 * x4 * (2.0 * x10 * x14 * x15 * x19 * x4 * x9 - x60)
    x225 = 0.5 * x224
    x226 = x219 * x91
    x227 = x219 * x96
    x228 = x219 * x4 * (2.0 * x10 * x14 * x15 * x4 * x85 * x9 - x86)
    x229 = 0.5 * x228
    x230 = x101 * x219
    x231 = x219 * x4 * (2.0 * x10 * x14 * x19 * x4 * x7 * x85 * x9 - x105)
    x232 = 0.5 * x231
    x233 = x112 * x219 * x68 + x67
    x234 = x233 * x27
    x235 = x123 * x219 + x175
    x236 = x112 * x219
    x237 = x236 * x36 + x64
    x238 = x237 * x27
    x239 = x118 * x219
    x240 = x181 + x239 * x55
    x241 = x236 * x55 + x53
    x242 = x4 * (x233 - x237 * x5)
    x243 = 0.5 * x242
    x244 = x241 * x27
    x245 = x102 * x239 + x189
    x246 = x245 * x3
    x247 = x240 * x5
    x248 = x4 * (x235 - x247)
    x249 = 0.5 * x248
    x250 = x244 * x3
    x251 = x219 * x4 * (x134 - x137)
    x252 = 0.5 * x251
    x253 = x211 * x219 + x96
    x254 = x213 * x219 + x91
    x255 = x254 * x5
    x256 = x101 + x148 * x219
    x257 = x4 * (x253 - x255)
    x258 = 0.5 * x257
    x259 = x112 * x29
    x260 = x152 * x219 + x259 * x39
    x261 = x154 * x219 + x259 * x26
    x262 = x261 * x5
    x263 = x157 * x219 + x259 * x46
    x264 = x263 * x72
    x265 = x4 * (x260 - x262)
    x266 = 0.5 * x265
    x267 = x163**2
    x268 = x267 * x44
    x269 = x268 + x50
    x270 = x267 * x55
    x271 = x19 * x270
    x272 = x271 + x61
    x273 = x27 * x272
    x274 = x267 * x37 - x5 * (x24 * x267 + x80) + x79
    x275 = x267 * x36
    x276 = x27 * (x275 + x71)
    x277 = 34.8371874529163 * x82
    x278 = x163 * x186 + x168 + x87
    x279 = x27 * x278
    x280 = x163 * x189
    x281 = x106 + x163 * x191 + x280
    x282 = x108 + x163 * x175 + x163 * x177 - x5 * (x109 + x163 * x184 + x168 * x19)
    x283 = 60.3397786612521 * x82
    x284 = x112 * x270 + x114
    x285 = x27 * x284
    x286 = x102 * x118 * x267 + x129
    x287 = x118 * x275 + x131 - x5 * (x118 * x270 + x132)
    x288 = x141 + x163 * x209 + 2.0 * x188
    x289 = x143 + x163 * x202 + 2.0 * x174 - x5 * (x144 + x163 * x205 + 2.0 * x180)
    x290 = x147 + x163 * x216 + x197
    x291 = x149 + x163 * x212 + x195 - x5 * (x150 + x163 * x214 + x194)
    x292 = x157 * x267 + x159
    x293 = x152 * x267 + x161 - x5 * (x154 * x267 + x162)
    x294 = x219 * x4 * (x163 * x37 - x166)
    x295 = x168 * x219
    x296 = 60.3397786612521 * x82
    x297 = x185 * x219 + x222
    x298 = x27 * x297
    x299 = x189 * x219
    x300 = x190 * x219 + x299
    x301 = x181 * x219
    x302 = x175 * x219
    x303 = x4 * (x176 * x219 + x302 - x5 * (x183 * x219 + x301))
    x304 = 104.511562358749 * x82
    x305 = x163 * x244
    x306 = x163 * x4 * (x235 - x247)
    x307 = x219 * (x207 + x208)
    x308 = x219 * x4 * (x199 + x201 - x5 * (x203 + x204))
    x309 = x163 * x256 + x244
    x310 = x4 * (x163 * x253 + x234 - x5 * (x163 * x254 + x238))
    x311 = x163 * x4 * (x260 - x262)
    x312 = x219**2
    x313 = x312 * x44 + x50
    x314 = x3 * x313
    x315 = x312 * x55
    x316 = x19 * x315 + x61
    x317 = x27 * x316
    x318 = 2.0 * x317
    x319 = x312 * x37 - x5 * (x24 * x312 + x80) + x79
    x320 = x319 * x49
    x321 = x312 * x36
    x322 = x27 * (x321 + x71)
    x323 = x315 * x85 + x87
    x324 = x27 * x323
    x325 = x102 * x312 * x92 + x106
    x326 = x108 + x321 * x92 - x5 * (x109 + x315 * x92)
    x327 = x326 * x49
    x328 = x114 + x219 * x241 + x222
    x329 = x27 * x328
    x330 = x129 + x219 * x245 + x299
    x331 = x3 * x330
    x332 = x131 + x219 * x235 + x302 - x5 * (x132 + x219 * x240 + x301)
    x333 = x332 * x49
    x334 = x139 * x312 + x141
    x335 = x134 * x312 + x143 - x5 * (x136 * x312 + x144)
    x336 = x335 * x49
    x337 = x147 + x219 * x256 + x230
    x338 = x149 + x219 * x253 + x227 - x5 * (x150 + x219 * x254 + x226)
    x339 = x338 * x49
    x340 = x159 + x219 * x263 + 2.0 * x244
    x341 = x161 + x219 * x260 + 2.0 * x234 - x5 * (x162 + x219 * x261 + 2.0 * x238)
    x342 = x341 * x49
    x343 = x163 * x269 + x167
    x344 = x27 * (x163 * x272 + x171)
    x345 = x27 * (x163 * x278 + x187 + x276)
    x346 = x163 * x281 + x192 + x273
    x347 = x27 * (x163 * x284 + x196)
    x348 = x163 * x286 + x198
    x349 = x163 * x288 + x210 + 2.0 * x279
    x350 = x172 * x3
    x351 = x163 * x290 + x217 + x285
    x352 = x193 * x3
    x353 = x163 * x292 + x218
    x354 = x219 * x268 + x221
    x355 = x27 * (x219 * x271 + x225)
    x356 = x27 * (x163 * x297 + x229 + x295)
    x357 = x163 * x300 + x219 * x280 + x232
    x358 = x27 * (x241 * x267 + x243)
    x359 = x245 * x267 + x249
    x360 = x163 * x307 + x252 + 2.0 * x298
    x361 = x296 * x3
    x362 = x163 * x309 + x258 + x305
    x363 = x3 * x304
    x364 = x263 * x267 + x266
    x365 = x27 * (x163 * x323 + x322)
    x366 = x163 * x325 + x317
    x367 = x163 * x329
    x368 = x163 * x334 + 2.0 * x324
    x369 = x163 * x337 + x329
    x370 = x219 * x313 + x220
    x371 = x27 * (x219 * x316 + x224)
    x372 = x27 * (x219 * x323 + x228)
    x373 = x219 * x325 + x231
    x374 = x27 * (x219 * x328 + x242 + x322)
    x375 = x219 * x330 + x248 + x317
    x376 = x219 * x334 + x251
    x377 = x219 * x337 + x257 + x324
    x378 = x219 * x340 + x265 + 2.0 * x329
    x379 = x163 * x172

    # 90 item(s)
    result[0, 0] = numpy.sum(
        x83
        * (
            x3
            * (
                x3 * (x3 * (x45 + x47) + x50 + x51 * x57)
                - x4 * (x32 * x5 - x41)
                + x51 * x63
            )
            + x51 * (x27 * (x71 + x73) + x3 * x63 - x4 * (x5 * x65 - x69))
            + x81 * (x3 * x41 - x5 * (x3 * x32 + x51 * x65 + x80) + x51 * x69 + x79)
        )
    )
    result[0, 1] = numpy.sum(
        x110
        * (
            x27 * x3 * (x4 * (2.0 * x10 * x14 * x15 * x4 * x85 * x9 - x86) + x90)
            + x3
            * (
                x27 * x90
                + x3 * (x106 + x107 + x3 * (x101 + x103))
                - x4 * (x5 * x95 - x99)
            )
            + x81 * (x108 + x3 * x96 + x3 * x99 - x5 * (x109 + x3 * x95 + x62 * x85))
        )
    )
    result[0, 2] = numpy.sum(
        x110
        * (
            x27 * x3 * (x116 + x4 * (2.0 * x10 * x112 * x14 * x15 * x4 * x9 - x113))
            + x3
            * (
                x116 * x27
                + x3 * (x129 + x130 + x3 * (x126 + x127))
                - x4 * (x121 * x5 - x125)
            )
            + x81 * (x122 * x3 + x125 * x3 + x131 - x5 * (x112 * x62 + x121 * x3 + x132))
        )
    )
    result[0, 3] = numpy.sum(
        x83
        * (
            x3 * (x3 * (x140 + x141) + x4 * (x134 * x3 - x138))
            + x81 * (x134 * x72 + x143 - x5 * (x136 * x72 + x144))
        )
    )
    result[0, 4] = numpy.sum(
        x110
        * (
            x3**2
            * (
                x147
                + x148 * x72
                + x4 * (2.0 * x10 * x112 * x14 * x4 * x7 * x85 * x9 - x146)
            )
            + x81 * (x145 * x73 + x149 - x5 * (x145 * x88 + x150))
        )
    )
    result[0, 5] = numpy.sum(
        x83
        * (
            x3 * (x3 * (x158 + x159) + x4 * (x152 * x3 - x156))
            + x81 * (x152 * x72 + x161 - x5 * (x154 * x72 + x162))
        )
    )
    result[1, 0] = numpy.sum(
        0.5
        * x172
        * (
            2.0 * x163 * x4 * (x38 + x40 - x5 * (x25 + x31))
            + x3 * (x167 + 2.0 * x170 * x51 + 2.0 * x3 * (x164 + x165))
            + x51 * (2.0 * x163 * x62 + 2.0 * x170 * x3 + x171)
        )
    )
    result[1, 1] = numpy.sum(
        0.5
        * x193
        * (
            x27 * (2.0 * x186 * x72 + x187)
            + x3 * (2.0 * x188 * x3 + x192 + 2.0 * x3 * (x188 + x191 * x3))
            + 2.0 * x4 * (x174 + x177 * x3 - x5 * (x180 + x184 * x3))
        )
    )
    result[1, 2] = numpy.sum(
        0.5
        * x193
        * (
            x27 * (2.0 * x115 * x163 + x196)
            + x3 * (2.0 * x130 * x163 + x198 + 2.0 * x3 * (x127 * x163 + x197))
            + 2.0 * x4 * (x124 * x163 + x195 - x5 * (x120 * x163 + x194))
        )
    )
    result[1, 3] = numpy.sum(
        0.5 * x172 * x3 * (2.0 * x209 * x72 + x210 + 2.0 * x4 * (x202 - x206))
    )
    result[1, 4] = numpy.sum(
        0.5 * x193 * x3 * (2.0 * x216 * x72 + x217 + 2.0 * x4 * (x212 - x215))
    )
    result[1, 5] = numpy.sum(
        0.5
        * x172
        * (2.0 * x163 * x4 * (x152 * x3 - x156) + x3 * (2.0 * x158 * x163 + x218))
    )
    result[2, 0] = numpy.sum(
        x172
        * (
            x219 * x4 * (x38 + x40 - x5 * (x25 + x31))
            + x3 * (x219 * x3 * (x45 + x47) + x221 + x223 * x51)
            + x51 * (x219 * x62 + x223 * x3 + x225)
        )
    )
    result[2, 1] = numpy.sum(
        x193
        * (
            x27 * (x219 * x89 + x229)
            + x3 * (x107 * x219 + x232 + x3 * (x103 * x219 + x230))
            + x4 * (x219 * x98 + x227 - x5 * (x219 * x94 + x226))
        )
    )
    result[2, 2] = numpy.sum(
        x193
        * (
            x27 * (x241 * x72 + x243)
            + x3 * (x249 + x250 + x3 * (x244 + x246))
            + x4 * (x234 + x235 * x3 - x5 * (x238 + x240 * x3))
        )
    )
    result[2, 3] = numpy.sum(
        x172 * (x219 * x4 * (x134 * x3 - x138) + x3 * (x140 * x219 + x252))
    )
    result[2, 4] = numpy.sum(x193 * x3 * (x256 * x72 + x258 + x4 * (x253 - x255)))
    result[2, 5] = numpy.sum(x172 * x3 * (x264 + x266 + x4 * (x260 - x262)))
    result[3, 0] = numpy.sum(
        x277 * (x274 * x49 + x3 * (x269 * x3 + 2.0 * x273) + x51 * (x272 * x3 + x276))
    )
    result[3, 1] = numpy.sum(x283 * (x279 * x3 + x282 * x49 + x3 * (x279 + x281 * x3)))
    result[3, 2] = numpy.sum(x283 * (x285 * x3 + x287 * x49 + x3 * (x285 + x286 * x3)))
    result[3, 3] = numpy.sum(x277 * (x288 * x72 + x289 * x49))
    result[3, 4] = numpy.sum(x283 * (x290 * x72 + x291 * x49))
    result[3, 5] = numpy.sum(x277 * (x292 * x72 + x293 * x49))
    result[4, 0] = numpy.sum(
        0.5
        * x296
        * (2.0 * x219 * x3 * (x164 + x165) + x294 + 2.0 * x51 * (x169 * x219 + x295))
    )
    result[4, 1] = numpy.sum(0.5 * x304 * (4.0 * x298 * x3 + 2.0 * x3**2 * x300 + x303))
    result[4, 2] = numpy.sum(
        0.5 * x304 * (2.0 * x163 * x250 + 2.0 * x3 * (x163 * x246 + x305) + x306)
    )
    result[4, 3] = numpy.sum(0.5 * x296 * (2.0 * x307 * x72 + x308))
    result[4, 4] = numpy.sum(0.5 * x304 * (2.0 * x309 * x72 + x310))
    result[4, 5] = numpy.sum(0.5 * x296 * (2.0 * x163 * x264 + x311))
    result[5, 0] = numpy.sum(
        x277 * (x3 * (x314 + x318) + x320 + x51 * (x3 * x316 + x322))
    )
    result[5, 1] = numpy.sum(x283 * (x3 * x324 + x3 * (x3 * x325 + x324) + x327))
    result[5, 2] = numpy.sum(x283 * (x3 * x329 + x3 * (x329 + x331) + x333))
    result[5, 3] = numpy.sum(x277 * (x334 * x72 + x336))
    result[5, 4] = numpy.sum(x283 * (x337 * x72 + x339))
    result[5, 5] = numpy.sum(x277 * (x340 * x72 + x342))
    result[6, 0] = numpy.sum(x172 * (x3 * x343 + 2.0 * x344))
    result[6, 1] = numpy.sum(x193 * (x3 * x346 + x345))
    result[6, 2] = numpy.sum(x193 * (x3 * x348 + x347))
    result[6, 3] = numpy.sum(x349 * x350)
    result[6, 4] = numpy.sum(x351 * x352)
    result[6, 5] = numpy.sum(x350 * x353)
    result[7, 0] = numpy.sum(x296 * (x3 * x354 + 2.0 * x355))
    result[7, 1] = numpy.sum(x304 * (x3 * x357 + x356))
    result[7, 2] = numpy.sum(x304 * (x3 * x359 + x358))
    result[7, 3] = numpy.sum(x360 * x361)
    result[7, 4] = numpy.sum(x362 * x363)
    result[7, 5] = numpy.sum(x361 * x364)
    result[8, 0] = numpy.sum(x163 * x296 * (x314 + x318))
    result[8, 1] = numpy.sum(x304 * (x3 * x366 + x365))
    result[8, 2] = numpy.sum(x304 * (x163 * x331 + x367))
    result[8, 3] = numpy.sum(x361 * x368)
    result[8, 4] = numpy.sum(x363 * x369)
    result[8, 5] = numpy.sum(x163 * x340 * x361)
    result[9, 0] = numpy.sum(x172 * (x3 * x370 + 2.0 * x371))
    result[9, 1] = numpy.sum(x193 * (x3 * x373 + x372))
    result[9, 2] = numpy.sum(x193 * (x3 * x375 + x374))
    result[9, 3] = numpy.sum(x350 * x376)
    result[9, 4] = numpy.sum(x352 * x377)
    result[9, 5] = numpy.sum(x350 * x378)
    result[10, 0] = numpy.sum(x83 * (x163 * x343 + x274 * x81))
    result[10, 1] = numpy.sum(x110 * (x163 * x346 + x282 * x81 + x344))
    result[10, 2] = numpy.sum(x110 * (x163 * x348 + x287 * x81))
    result[10, 3] = numpy.sum(x83 * (x163 * x349 + x289 * x81 + 2.0 * x345))
    result[10, 4] = numpy.sum(x110 * (x163 * x351 + x291 * x81 + x347))
    result[10, 5] = numpy.sum(x83 * (x163 * x353 + x293 * x81))
    result[11, 0] = numpy.sum(x172 * (x163 * x354 + x294))
    result[11, 1] = numpy.sum(x193 * (x163 * x357 + x303 + x355))
    result[11, 2] = numpy.sum(x193 * (x163 * x359 + x306))
    result[11, 3] = numpy.sum(x172 * (x163 * x360 + x308 + 2.0 * x356))
    result[11, 4] = numpy.sum(x193 * (x163 * x362 + x310 + x358))
    result[11, 5] = numpy.sum(x172 * (x163 * x364 + x311))
    result[12, 0] = numpy.sum(x277 * (x267 * x313 + x320))
    result[12, 1] = numpy.sum(x283 * (x163 * x317 + x163 * x366 + x327))
    result[12, 2] = numpy.sum(x283 * (x267 * x330 + x333))
    result[12, 3] = numpy.sum(x277 * (x163 * x368 + x336 + 2.0 * x365))
    result[12, 4] = numpy.sum(x283 * (x163 * x369 + x339 + x367))
    result[12, 5] = numpy.sum(x277 * (x267 * x340 + x342))
    result[13, 0] = numpy.sum(x370 * x379)
    result[13, 1] = numpy.sum(x193 * (x163 * x373 + x371))
    result[13, 2] = numpy.sum(x163 * x193 * x375)
    result[13, 3] = numpy.sum(x172 * (x163 * x376 + 2.0 * x372))
    result[13, 4] = numpy.sum(x193 * (x163 * x377 + x374))
    result[13, 5] = numpy.sum(x378 * x379)
    result[14, 0] = numpy.sum(x83 * (x219 * x370 + x319 * x81))
    result[14, 1] = numpy.sum(x110 * (x219 * x373 + x326 * x81))
    result[14, 2] = numpy.sum(x110 * (x219 * x375 + x332 * x81 + x371))
    result[14, 3] = numpy.sum(x83 * (x219 * x376 + x335 * x81))
    result[14, 4] = numpy.sum(x110 * (x219 * x377 + x338 * x81 + x372))
    result[14, 5] = numpy.sum(x83 * (x219 * x378 + x341 * x81 + 2.0 * x374))
    return result


def _2center2el3d_43(ax, da, A, bx, db, B):
    """Cartesian (g|f) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 10), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = ax ** (-1.0)
    x5 = bx ** (-1.0)
    x6 = -x2 - B[0]
    x7 = bx * x1
    x8 = ax * x7 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x9 = boys(5, x8)
    x10 = x0 ** (-1.5)
    x11 = 17.4934183276249
    x12 = x11 * x5
    x13 = 2.0 * x12
    x14 = x10 * x13
    x15 = x14 * x9
    x16 = x0 ** (-0.5)
    x17 = boys(4, x8)
    x18 = 0.5 * x5
    x19 = x18 * (2.0 * x11 * x16 * x17 * x4 * x5 - x15)
    x20 = boys(6, x8)
    x21 = x6**2
    x22 = x16 * x4
    x23 = x13 * x22
    x24 = x21 * x23
    x25 = x20 * x24
    x26 = x19 + x25
    x27 = x6 * (x26 + x5 * (2.0 * x11 * x16 * x17 * x4 * x5 - x15))
    x28 = x27 * x3
    x29 = 0.5 / (ax + bx)
    x30 = x14 * x17
    x31 = boys(3, x8)
    x32 = x18 * (2.0 * x11 * x16 * x31 * x4 * x5 - x30)
    x33 = x23 * x9
    x34 = x21 * x33
    x35 = x32 + x34
    x36 = x29 * x35
    x37 = 3.0 * x36
    x38 = x28 + x37
    x39 = -2.0 * x11 * x16 * x31 * x4 * x5 * x6
    x40 = x35 * x6 - x5 * (x30 * x6 + x39)
    x41 = x3 * x40
    x42 = x14 * x31
    x43 = boys(2, x8)
    x44 = x18 * (2.0 * x11 * x16 * x4 * x43 * x5 - x42)
    x45 = x17 * x23
    x46 = x21 * x45
    x47 = x44 + x46
    x48 = x29 * x47
    x49 = 3.0 * x48
    x50 = x41 + x49
    x51 = x14 * x20
    x52 = x18 * (2.0 * x11 * x16 * x4 * x5 * x9 - x51)
    x53 = boys(7, x8)
    x54 = x24 * x53
    x55 = x6 * (x5 * (2.0 * x11 * x16 * x4 * x5 * x9 - x51) + x52 + x54)
    x56 = x3 * x55
    x57 = x26 * x29
    x58 = 3.0 * x57
    x59 = x27 * x7
    x60 = 0.5 * x4
    x61 = x60 * (x40 - x59)
    x62 = x26 * x3
    x63 = x6 * x9
    x64 = x12 * x22
    x65 = 4.0 * x29 * x64
    x66 = x63 * x65
    x67 = x62 + x66
    x68 = 3.0 * x29
    x69 = x35 * x7
    x70 = x60 * (x47 - x69)
    x71 = 2.0 * x29
    x72 = x64 * x71
    x73 = x17 * x72
    x74 = x3 * x6
    x75 = x33 * x74
    x76 = x73 + x75
    x77 = x3 * x67 + x70 + x71 * x76
    x78 = x17 * x6
    x79 = x65 * x78
    x80 = x3 * x35 + x79
    x81 = x31 * x6
    x82 = x3 * x47 + x65 * x81
    x83 = 2.0 * x10 * x11 * x4
    x84 = x78 * x83
    x85 = -x60 * (x39 + x84)
    x86 = x3 * x73
    x87 = x14 * x43
    x88 = boys(1, x8)
    x89 = x18 * (2.0 * x11 * x16 * x4 * x5 * x88 - x87)
    x90 = x23 * x31
    x91 = x21 * x90
    x92 = x6 * (x47 + x5 * (2.0 * x11 * x16 * x4 * x43 * x5 - x42))
    x93 = x60 * (
        x5 * x6 * (2.0 * x11 * x16 * x4 * x5 * x88 - x87) + x6 * (x89 + x91) - x7 * x92
    )
    x94 = x60 * (-x40 * x7 + x92)
    x95 = 1.5 * x4
    x96 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**4.5)
    x97 = 9.12251705727742 * x96
    x98 = -x1 * (ax * A[1] + bx * B[1])
    x99 = -x98 - B[1]
    x100 = x5 * x99 * (2.0 * x11 * x16 * x17 * x4 * x5 - x15)
    x101 = 0.5 * x100 + x25 * x99
    x102 = x101 * x3
    x103 = x66 * x99
    x104 = x102 + x103
    x105 = -2.0 * x11 * x16 * x31 * x4 * x5 * x99
    x106 = -x5 * (x105 + x30 * x99)
    x107 = 0.5 * x106 + x34 * x99
    x108 = x107 * x3
    x109 = x79 * x99
    x110 = x108 + x109
    x111 = x5 * x99 * (2.0 * x11 * x16 * x4 * x5 * x9 - x51)
    x112 = 0.5 * x111 + x54 * x99
    x113 = x112 * x3
    x114 = x20 * x99
    x115 = x114 * x6
    x116 = x115 * x65
    x117 = x101 * x7
    x118 = x60 * (x107 - x117)
    x119 = x9 * x99
    x120 = x119 * x72
    x121 = x23 * x74
    x122 = x114 * x121
    x123 = x120 + x122
    x124 = x63 * x83
    x125 = x124 * x99
    x126 = x60 * (2.0 * x11 * x16 * x17 * x4 * x5 * x6 * x99 - x125)
    x127 = x120 * x3
    x128 = x123 * x3 + x126 + x127
    x129 = x73 * x99
    x130 = x33 * x99
    x131 = x130 * x74
    x132 = x129 + x131
    x133 = x31 * x72
    x134 = x133 * x99
    x135 = x45 * x99
    x136 = x134 + x135 * x74
    x137 = x17 * x83
    x138 = -x60 * (x105 + x137 * x99)
    x139 = x3**2
    x140 = x5 * x99 * (2.0 * x11 * x16 * x4 * x5 * x88 - x87)
    x141 = x5 * x99 * (2.0 * x11 * x16 * x4 * x43 * x5 - x42)
    x142 = 0.5 * x141 + x46 * x99
    x143 = 0.5 * x60 * (x140 - 2.0 * x142 * x7 + 2.0 * x91 * x99)
    x144 = x60 * (-x107 * x7 + x142)
    x145 = 20.3985682659737 * x96
    x146 = -x1 * (ax * A[2] + bx * B[2])
    x147 = -x146 - B[2]
    x148 = x147 * x5 * (2.0 * x11 * x16 * x17 * x4 * x5 - x15)
    x149 = 0.5 * x148
    x150 = x147 * x25 + x149
    x151 = x150 * x3
    x152 = x147 * x66
    x153 = x151 + x152
    x154 = -2.0 * x11 * x147 * x16 * x31 * x4 * x5
    x155 = -x5 * (x147 * x30 + x154)
    x156 = 0.5 * x155
    x157 = x147 * x34 + x156
    x158 = x157 * x3
    x159 = x147 * x79
    x160 = x158 + x159
    x161 = x147 * x5 * (2.0 * x11 * x16 * x4 * x5 * x9 - x51)
    x162 = 0.5 * x161
    x163 = x147 * x54 + x162
    x164 = x163 * x3
    x165 = x147 * x65
    x166 = x20 * x6
    x167 = x165 * x166
    x168 = x150 * x7
    x169 = x60 * (x157 - x168)
    x170 = x147 * x9
    x171 = x170 * x72
    x172 = x121 * x147 * x20
    x173 = x171 + x172
    x174 = x124 * x147
    x175 = x60 * (2.0 * x11 * x147 * x16 * x17 * x4 * x5 * x6 - x174)
    x176 = x171 * x3
    x177 = x173 * x3 + x175 + x176
    x178 = x147 * x73
    x179 = x147 * x33
    x180 = x178 + x179 * x74
    x181 = x133 * x147
    x182 = x147 * x45
    x183 = x181 + x182 * x74
    x184 = -x60 * (x137 * x147 + x154)
    x185 = x147 * x5 * (2.0 * x11 * x16 * x4 * x5 * x88 - x87)
    x186 = 0.5 * x185
    x187 = x147 * x5 * (2.0 * x11 * x16 * x4 * x43 * x5 - x42)
    x188 = 0.5 * x187
    x189 = x147 * x46 + x188
    x190 = x60 * (x147 * x91 + x186 - x189 * x7)
    x191 = x60 * (-x157 * x7 + x189)
    x192 = x99**2
    x193 = x192 * x45 + x44
    x194 = x192 * x33 + x32
    x195 = x194 * x7
    x196 = x192 * x23
    x197 = x196 * x20
    x198 = x19 + x197
    x199 = x139 * x198
    x200 = x60 * (x193 - x195)
    x201 = x199 + x200
    x202 = x194 * x29
    x203 = x198 * x74
    x204 = x202 + x203
    x205 = x193 * x29
    x206 = x194 * x6
    x207 = x206 * x3
    x208 = x205 + x207
    x209 = x198 * x29
    x210 = x196 * x53
    x211 = x210 + x52
    x212 = x211 * x74
    x213 = x6 * x7
    x214 = x198 * x213
    x215 = x60 * (x194 * x6 - x214)
    x216 = x209 * x3
    x217 = x192 * x90 + x89
    x218 = x6 * x60 * (-x193 * x7 + x217)
    x219 = x60 * (x193 * x6 - x206 * x7)
    x220 = x147 * x83
    x221 = x119 * x220
    x222 = x60 * (2.0 * x11 * x147 * x16 * x17 * x4 * x5 * x99 - x221)
    x223 = x114 * x147
    x224 = x223 * x23
    x225 = x139 * x224 + x222
    x226 = x147 * (x120 + x122)
    x227 = x147 * (x129 + x131)
    x228 = x223 * x72
    x229 = x147 * x99
    x230 = x229 * x53
    x231 = x60 * (2.0 * x11 * x147 * x16 * x4 * x5 * x6 * x9 * x99 - x115 * x220)
    x232 = x60 * (2.0 * x11 * x147 * x16 * x31 * x4 * x5 * x6 * x99 - x229 * x84)
    x233 = x147 * x60 * (2.0 * x11 * x16 * x17 * x4 * x5 * x6 * x99 - x125)
    x234 = 35.3313566383285 * x96
    x235 = x147**2
    x236 = x235 * x45 + x44
    x237 = x235 * x33 + x32
    x238 = x237 * x7
    x239 = x23 * x235
    x240 = x19 + x20 * x239
    x241 = x139 * x240
    x242 = x60 * (x236 - x238)
    x243 = x241 + x242
    x244 = x237 * x29
    x245 = x240 * x74
    x246 = x244 + x245
    x247 = x236 * x29
    x248 = x237 * x6
    x249 = x248 * x3
    x250 = x247 + x249
    x251 = x240 * x29
    x252 = x239 * x53 + x52
    x253 = x252 * x74
    x254 = x213 * x240
    x255 = x60 * (x237 * x6 - x254)
    x256 = x251 * x3
    x257 = x235 * x90 + x89
    x258 = x6 * x60 * (-x236 * x7 + x257)
    x259 = x60 * (x236 * x6 - x248 * x7)
    x260 = x106 + x194 * x99
    x261 = x100 + x198 * x99
    x262 = x261 * x7
    x263 = x262 * x3
    x264 = x111 + x211 * x99
    x265 = x139 * x264
    x266 = x60 * (x260 - x262)
    x267 = x141 + x193 * x99
    x268 = x60 * (x140 + x217 * x99 - x267 * x7)
    x269 = x60 * (-x260 * x7 + x267)
    x270 = x156 + x179 * x192
    x271 = x147 * x197 + x149
    x272 = x271 * x7
    x273 = x147 * x210 + x162
    x274 = x60 * (x270 - x272)
    x275 = x182 * x192 + x188
    x276 = x60 * (x147 * x192 * x90 + x186 - x275 * x7)
    x277 = x60 * (-x270 * x7 + x275)
    x278 = x240 * x99
    x279 = x278 * x7
    x280 = x60 * (x237 * x99 - x279)
    x281 = x252 * x99
    x282 = x60 * x99 * (-x236 * x7 + x257)
    x283 = x237 * x99
    x284 = x60 * (x236 * x99 - x283 * x7)
    x285 = x147 * x237 + x155
    x286 = x147 * x240 + x148
    x287 = x286 * x7
    x288 = x287 * x3
    x289 = x147 * x252 + x161
    x290 = x139 * x289
    x291 = x60 * (x285 - x287)
    x292 = x147 * x236 + x187
    x293 = x60 * (x147 * x257 + x185 - x292 * x7)
    x294 = x60 * (-x285 * x7 + x292)
    x295 = -x98 - A[1]
    x296 = x295 * x47
    x297 = x295 * x56
    x298 = x295 * x58
    x299 = x295 * x59
    x300 = x4 * (x295 * x40 - x299)
    x301 = x295 * x62
    x302 = x295 * x66
    x303 = x301 + x302
    x304 = x295 * x4 * (x47 - x69)
    x305 = x295 * x73
    x306 = 24.1359114645008 * x96
    x307 = x107 * x295
    x308 = x307 + x48
    x309 = x72 * x81
    x310 = x45 * x6
    x311 = x295 * x310 * x99 + x309
    x312 = x311 * x71
    x313 = x101 * x295
    x314 = x313 + x36
    x315 = x72 * x78
    x316 = x33 * x6
    x317 = x295 * x316 * x99 + x315
    x318 = x317 * x71
    x319 = x112 * x295
    x320 = x319 + x57
    x321 = x63 * x72
    x322 = x115 * x23
    x323 = x295 * x322
    x324 = x321 + x323
    x325 = x324 * x71
    x326 = x4 * (x308 - x314 * x7)
    x327 = x130 * x295
    x328 = x327 + x73
    x329 = x29 * x328
    x330 = x3 * x324 + x329
    x331 = x4 * (x311 - x317 * x7)
    x332 = 53.9695387335403 * x96
    x333 = x295 * x4 * (x157 - x168)
    x334 = x171 * x295
    x335 = x172 * x295 + x334
    x336 = x295 * x4 * (2.0 * x11 * x147 * x16 * x17 * x4 * x5 * x6 - x174)
    x337 = x65 * x99
    x338 = x193 * x295 + x31 * x337
    x339 = x29 * x338
    x340 = x206 * x295
    x341 = x109 + x340
    x342 = x17 * x337 + x194 * x295
    x343 = x29 * x342
    x344 = x198 * x295
    x345 = x344 * x6
    x346 = x103 + x345
    x347 = x119 * x65
    x348 = x344 + x347
    x349 = x4 * (x338 - x342 * x7)
    x350 = x29 * x348
    x351 = x211 * x6
    x352 = x295 * x351
    x353 = x116 + x352
    x354 = x4 * (x341 - x346 * x7)
    x355 = x135 * x147
    x356 = x181 + x295 * x355
    x357 = x29 * x356
    x358 = x147 * x315
    x359 = x229 * x316
    x360 = x295 * x359 + x358
    x361 = x147 * x327 + x178
    x362 = x29 * x361
    x363 = x147 * x321
    x364 = x147 * x323 + x363
    x365 = x171 + x224 * x295
    x366 = x4 * (x356 - x361 * x7)
    x367 = x29 * x365
    x368 = x147 * x166
    x369 = x368 * x72
    x370 = x23 * x230 * x6
    x371 = x295 * x370 + x369
    x372 = x4 * (x360 - x364 * x7)
    x373 = 93.4779831475484 * x96
    x374 = x247 * x295
    x375 = x244 * x295
    x376 = x295 * x4 * (x236 - x238)
    x377 = x251 * x295
    x378 = x295 * x4 * (x237 * x6 - x254)
    x379 = x260 * x295
    x380 = 3.0 * x205
    x381 = x379 + x380
    x382 = x261 * x295
    x383 = 3.0 * x202
    x384 = x382 + x383
    x385 = x384 * x7
    x386 = x264 * x295
    x387 = 3.0 * x209
    x388 = x386 + x387
    x389 = x4 * (x381 - x385)
    x390 = x17 * x229 * x65
    x391 = x270 * x295 + x390
    x392 = x147 * x347
    x393 = x271 * x295 + x392
    x394 = x393 * x7
    x395 = x114 * x165
    x396 = x273 * x295 + x395
    x397 = x4 * (x391 - x394)
    x398 = x247 + x283 * x295
    x399 = x244 + x278 * x295
    x400 = x399 * x7
    x401 = x251 + x281 * x295
    x402 = x4 * (x398 - x400)
    x403 = x295 * x4 * (x285 - x287)
    x404 = -x146 - A[2]
    x405 = x4 * x404 * (x40 - x59)
    x406 = 0.5 * x405
    x407 = x404 * (x62 + x66)
    x408 = x4 * x404 * (x47 - x69)
    x409 = 0.5 * x408
    x410 = x404 * x73
    x411 = x103 * x404
    x412 = x109 * x404
    x413 = x116 * x404
    x414 = x4 * x404 * (x107 - x117)
    x415 = 0.5 * x414
    x416 = x120 * x404
    x417 = x122 * x404 + x416
    x418 = x4 * x404 * (2.0 * x11 * x16 * x17 * x4 * x5 * x6 * x99 - x125)
    x419 = 0.5 * x418
    x420 = x157 * x404 + x48
    x421 = x147 * x310 * x404 + x309
    x422 = x29 * x421
    x423 = 2.0 * x422
    x424 = x150 * x404 + x36
    x425 = x147 * x404
    x426 = x315 + x316 * x425
    x427 = x29 * x426
    x428 = 2.0 * x427
    x429 = x163 * x404 + x57
    x430 = x3 * x429
    x431 = x23 * x368
    x432 = x321 + x404 * x431
    x433 = x29 * x432
    x434 = 2.0 * x433
    x435 = x424 * x7
    x436 = x4 * (x420 - x435)
    x437 = 0.5 * x436
    x438 = x179 * x404 + x73
    x439 = x29 * x438
    x440 = x3 * x432
    x441 = x439 + x440
    x442 = x4 * (x421 - x426 * x7)
    x443 = 0.5 * x442
    x444 = x205 * x404
    x445 = x202 * x404
    x446 = x4 * x404 * (x193 - x195)
    x447 = 0.5 * x446
    x448 = x209 * x404
    x449 = x4 * x404 * (x194 * x6 - x214)
    x450 = 0.5 * x449
    x451 = x134 + x355 * x404
    x452 = x29 * x451
    x453 = x315 * x99
    x454 = x359 * x404 + x453
    x455 = x129 + x130 * x425
    x456 = x29 * x455
    x457 = x321 * x99
    x458 = x322 * x425 + x457
    x459 = x120 + x224 * x404
    x460 = x4 * (x451 - x455 * x7)
    x461 = 0.5 * x460
    x462 = x29 * x459
    x463 = x115 * x72
    x464 = x370 * x404 + x463
    x465 = x4 * (x454 - x458 * x7)
    x466 = 0.5 * x465
    x467 = x165 * x31 + x236 * x404
    x468 = x29 * x467
    x469 = x159 + x248 * x404
    x470 = x165 * x17 + x237 * x404
    x471 = x29 * x470
    x472 = x240 * x404
    x473 = x152 + x472 * x6
    x474 = x170 * x65 + x472
    x475 = x4 * (x467 - x470 * x7)
    x476 = 0.5 * x475
    x477 = x29 * x474
    x478 = x252 * x6
    x479 = x167 + x404 * x478
    x480 = x3 * x479
    x481 = x473 * x7
    x482 = x4 * (x469 - x481)
    x483 = 0.5 * x482
    x484 = x3 * x477
    x485 = x4 * x404 * (x260 - x262)
    x486 = 0.5 * x485
    x487 = x205 + x270 * x404
    x488 = x202 + x271 * x404
    x489 = x488 * x7
    x490 = x209 + x273 * x404
    x491 = x4 * (x487 - x489)
    x492 = 0.5 * x491
    x493 = x283 * x404 + x390
    x494 = x392 + x472 * x99
    x495 = x494 * x7
    x496 = x281 * x404 + x395
    x497 = x4 * (x493 - x495)
    x498 = 0.5 * x497
    x499 = 3.0 * x247 + x285 * x404
    x500 = 3.0 * x244 + x286 * x404
    x501 = x500 * x7
    x502 = 3.0 * x251 + x289 * x404
    x503 = x139 * x502
    x504 = x4 * (x499 - x501)
    x505 = 0.5 * x504
    x506 = x295**2
    x507 = x506 * x55
    x508 = x507 + x61
    x509 = x26 * x506
    x510 = x509 + x70
    x511 = x29 * x510
    x512 = x40 * x506 - x7 * (x27 * x506 + x94) + x93
    x513 = x29 * (x316 * x506 + x85)
    x514 = 31.1593277158494 * x96
    x515 = x295 * x57
    x516 = x118 + x295 * x320 + x515
    x517 = x295 * x321
    x518 = x126 + x295 * x324 + x517
    x519 = x518 * x71
    x520 = x143 + x29 * x296 + x295 * x308 - x7 * (x144 + x295 * x314 + x295 * x36)
    x521 = x29 * (x138 + x295 * x328 + x305)
    x522 = 69.6743749058326 * x96
    x523 = x163 * x506 + x169
    x524 = x175 + x431 * x506
    x525 = x29 * x524
    x526 = x157 * x506 + x190 - x7 * (x150 * x506 + x191)
    x527 = x29 * (x179 * x506 + x184)
    x528 = x200 + x295 * x348 + 2.0 * x329
    x529 = x29 * x528
    x530 = x215 + x295 * x353 + x325
    x531 = x218 + x295 * x341 + x312 - x7 * (x219 + x295 * x346 + x318)
    x532 = x222 + x295 * x365 + x334
    x533 = x29 * x532
    x534 = x231 + x295 * x369 + x295 * x371
    x535 = x232 + x295 * x358 + x295 * x360 - x7 * (x233 + x295 * x363 + x295 * x364)
    x536 = 120.679557322504 * x96
    x537 = x240 * x506
    x538 = x242 + x537
    x539 = x29 * x538
    x540 = x255 + x478 * x506
    x541 = x248 * x506 + x258 - x7 * (x259 + x537 * x6)
    x542 = x266 + x295 * x388 + 3.0 * x350
    x543 = x268 + x295 * x381 + 3.0 * x339 - x7 * (x269 + x295 * x384 + 3.0 * x343)
    x544 = x274 + x295 * x396 + 2.0 * x367
    x545 = x276 + x295 * x391 + 2.0 * x357 - x7 * (x277 + x295 * x393 + 2.0 * x362)
    x546 = x280 + x295 * x401 + x377
    x547 = x282 + x295 * x398 + x374 - x7 * (x284 + x295 * x399 + x375)
    x548 = x289 * x506 + x291
    x549 = x285 * x506 + x293 - x7 * (x286 * x506 + x294)
    x550 = x4 * x404 * (x295 * x40 - x299)
    x551 = x404 * x57
    x552 = x319 * x404 + x551
    x553 = x321 * x404
    x554 = x323 * x404 + x553
    x555 = x554 * x71
    x556 = x404 * x48
    x557 = x36 * x404
    x558 = x4 * (x307 * x404 + x556 - x7 * (x313 * x404 + x557))
    x559 = x29 * (x327 * x404 + x410)
    x560 = 120.679557322504 * x96
    x561 = x295 * x4 * (x420 - x435)
    x562 = x295 * x439
    x563 = x404 * (x344 + x347)
    x564 = x29 * x563
    x565 = x352 * x404 + x413
    x566 = x4 * (x340 * x404 + x412 - x7 * (x345 * x404 + x411))
    x567 = x295 * x459 + x439
    x568 = x29 * x567
    x569 = x295 * x464 + x433
    x570 = x4 * (x295 * x454 + x422 - x7 * (x295 * x458 + x427))
    x571 = 209.023124717498 * x96
    x572 = x295 * x477
    x573 = x295 * x4 * (x469 - x481)
    x574 = x404 * (x386 + x387)
    x575 = x4 * x404 * (x379 + x380 - x7 * (x382 + x383))
    x576 = 2.0 * x462
    x577 = x295 * x490 + x576
    x578 = 2.0 * x452
    x579 = 2.0 * x456
    x580 = x4 * (x295 * x487 + x578 - x7 * (x295 * x488 + x579))
    x581 = x295 * x496 + x477
    x582 = x4 * (x295 * x493 + x468 - x7 * (x295 * x494 + x471))
    x583 = x295 * x4 * (x499 - x501)
    x584 = x404**2
    x585 = x55 * x584 + x61
    x586 = x3 * x585
    x587 = x26 * x584 + x70
    x588 = x29 * x587
    x589 = 3.0 * x588
    x590 = x40 * x584 - x7 * (x27 * x584 + x94) + x93
    x591 = x590 * x60
    x592 = x29 * (x316 * x584 + x85)
    x593 = x112 * x584 + x118
    x594 = x126 + x322 * x584
    x595 = x29 * x594
    x596 = 2.0 * x595
    x597 = x107 * x584 + x143 - x7 * (x101 * x584 + x144)
    x598 = x597 * x60
    x599 = x29 * (x130 * x584 + x138)
    x600 = x169 + x404 * x429 + x551
    x601 = x3 * x600
    x602 = x175 + x404 * x432 + x553
    x603 = x29 * x602
    x604 = 2.0 * x603
    x605 = x190 + x404 * x420 + x556 - x7 * (x191 + x404 * x424 + x557)
    x606 = x60 * x605
    x607 = x29 * (x184 + x404 * x438 + x410)
    x608 = x198 * x584
    x609 = x200 + x608
    x610 = x29 * x609
    x611 = x215 + x351 * x584
    x612 = x206 * x584 + x218 - x7 * (x219 + x6 * x608)
    x613 = x60 * x612
    x614 = x222 + x404 * x459 + x416
    x615 = x29 * x614
    x616 = x231 + x404 * x463 + x404 * x464
    x617 = x232 + x404 * x453 + x404 * x454 - x7 * (x233 + x404 * x457 + x404 * x458)
    x618 = x60 * x617
    x619 = x242 + x404 * x474 + 2.0 * x439
    x620 = x29 * x619
    x621 = x255 + x404 * x479 + x434
    x622 = x3 * x621
    x623 = x258 + x404 * x469 + x423 - x7 * (x259 + x404 * x473 + x428)
    x624 = x60 * x623
    x625 = x264 * x584 + x266
    x626 = x260 * x584 + x268 - x7 * (x261 * x584 + x269)
    x627 = x60 * x626
    x628 = x274 + x404 * x490 + x448
    x629 = x276 + x404 * x487 + x444 - x7 * (x277 + x404 * x488 + x445)
    x630 = x60 * x629
    x631 = x280 + x404 * x496 + x576
    x632 = x282 + x404 * x493 + x578 - x7 * (x284 + x404 * x494 + x579)
    x633 = x60 * x632
    x634 = x291 + x404 * x502 + 3.0 * x477
    x635 = x293 + x404 * x499 + 3.0 * x468 - x7 * (x294 + x404 * x500 + 3.0 * x471)
    x636 = x60 * x635
    x637 = x295 * x508 + x300
    x638 = x29 * (x295 * x510 + x304)
    x639 = x295 * x516 + x326 + x511
    x640 = x71 * (x295 * x518 + x331 + x513)
    x641 = x295 * x523 + x333
    x642 = x29 * (x295 * x524 + x336)
    x643 = x29 * (x295 * x528 + x349 + 2.0 * x521)
    x644 = x295 * x530 + x354 + x519
    x645 = x29 * (x295 * x532 + x366 + x527)
    x646 = x295 * x534 + x372 + x525
    x647 = x29 * (x295 * x538 + x376)
    x648 = x295 * x540 + x378
    x649 = x295 * x542 + x389 + 3.0 * x529
    x650 = x3 * x306
    x651 = x295 * x544 + x397 + 2.0 * x533
    x652 = x3 * x332
    x653 = x295 * x546 + x402 + x539
    x654 = x295 * x548 + x403
    x655 = x404 * x507 + x406
    x656 = x29 * (x404 * x509 + x409)
    x657 = x295 * x552 + x404 * x515 + x415
    x658 = x71 * (x295 * x554 + x404 * x517 + x419)
    x659 = x429 * x506 + x437
    x660 = x29 * (x432 * x506 + x443)
    x661 = x29 * (x295 * x563 + x447 + 2.0 * x559)
    x662 = x295 * x565 + x450 + x555
    x663 = x29 * (x295 * x567 + x461 + x562)
    x664 = x295 * x433 + x295 * x569 + x466
    x665 = x29 * (x474 * x506 + x476)
    x666 = x479 * x506 + x483
    x667 = x295 * x574 + x486 + 3.0 * x564
    x668 = x295 * x577 + x492 + 2.0 * x568
    x669 = x3 * x560
    x670 = x295 * x581 + x498 + x572
    x671 = x502 * x506 + x505
    x672 = x295 * x593 + x588
    x673 = x71 * (x295 * x594 + x592)
    x674 = x29 * (x295 * x609 + 2.0 * x599)
    x675 = x295 * x611 + x596
    x676 = x29 * (x295 * x614 + x607)
    x677 = x295 * x616 + x603
    x678 = x295 * x620
    x679 = x295 * x625 + 3.0 * x610
    x680 = 2.0 * x615
    x681 = x295 * x628 + x680
    x682 = x295 * x631 + x620
    x683 = x404 * x585 + x405
    x684 = x29 * (x404 * x587 + x408)
    x685 = x404 * x593 + x414
    x686 = x29 * (x404 * x594 + x418)
    x687 = 2.0 * x686
    x688 = x404 * x600 + x436 + x588
    x689 = x29 * (x404 * x602 + x442 + x592)
    x690 = 2.0 * x689
    x691 = x29 * (x404 * x609 + x446)
    x692 = x404 * x611 + x449
    x693 = x29 * (x404 * x614 + x460 + x599)
    x694 = x404 * x616 + x465 + x595
    x695 = x29 * (x404 * x619 + x475 + 2.0 * x607)
    x696 = x404 * x621 + x482 + x604
    x697 = x404 * x625 + x485
    x698 = x404 * x628 + x491 + x610
    x699 = x404 * x631 + x497 + x680
    x700 = x404 * x634 + x504 + 3.0 * x620
    x701 = x295 * x306
    x702 = x295 * x332
    x703 = 2.0 * x693

    # 150 item(s)
    result[0, 0] = numpy.sum(
        x97
        * (
            x3
            * (
                x3 * (x3 * (x56 + x58) + x61 + x67 * x68)
                - x4 * (x38 * x7 - x50)
                + x68 * x77
            )
            + x68 * (x3 * x77 - x4 * (x7 * x80 - x82) + x71 * (x3 * x76 + x85 + x86))
            + x95 * (x3 * x50 + x68 * x82 - x7 * (x3 * x38 + x68 * x80 + x94) + x93)
        )
    )
    result[0, 1] = numpy.sum(
        x145
        * (
            x3
            * (
                x128 * x71
                + x3 * (x118 + x123 * x71 + x3 * (x113 + x116))
                - x4 * (x104 * x7 - x110)
            )
            + x71 * (x128 * x3 + x29 * (x130 * x139 + x138) - x4 * (x132 * x7 - x136))
            + x95 * (x110 * x3 + x136 * x71 + x143 - x7 * (x104 * x3 + x132 * x71 + x144))
        )
    )
    result[0, 2] = numpy.sum(
        x145
        * (
            x3
            * (
                x177 * x71
                + x3 * (x169 + x173 * x71 + x3 * (x164 + x167))
                - x4 * (x153 * x7 - x160)
            )
            + x71 * (x177 * x3 + x29 * (x139 * x179 + x184) - x4 * (x180 * x7 - x183))
            + x95 * (x160 * x3 + x183 * x71 + x190 - x7 * (x153 * x3 + x180 * x71 + x191))
        )
    )
    result[0, 3] = numpy.sum(
        x145
        * (
            x29 * x3 * (x201 + x4 * (x193 - x195))
            + x3
            * (
                x201 * x29
                + x3 * (x215 + x216 + x3 * (x209 + x212))
                - x4 * (x204 * x7 - x208)
            )
            + x95
            * (x193 * x29 * x3 + x208 * x3 + x218 - x7 * (x202 * x3 + x204 * x3 + x219))
        )
    )
    result[0, 4] = numpy.sum(
        x234
        * (
            x29 * x3 * (x225 + x4 * (2.0 * x11 * x147 * x16 * x17 * x4 * x5 * x99 - x221))
            + x3
            * (
                x225 * x29
                + x3 * (x228 * x3 + x231 + x3 * (x121 * x230 + x228))
                - x4 * (x226 * x7 - x227)
            )
            + x95
            * (x227 * x3 + x229 * x86 + x232 - x7 * (x127 * x147 + x226 * x3 + x233))
        )
    )
    result[0, 5] = numpy.sum(
        x145
        * (
            x29 * x3 * (x243 + x4 * (x236 - x238))
            + x3
            * (
                x243 * x29
                + x3 * (x255 + x256 + x3 * (x251 + x253))
                - x4 * (x246 * x7 - x250)
            )
            + x95
            * (x236 * x29 * x3 + x250 * x3 + x258 - x7 * (x244 * x3 + x246 * x3 + x259))
        )
    )
    result[0, 6] = numpy.sum(
        x97
        * (
            x3 * (x3 * (x265 + x266) + x4 * (x260 * x3 - x263))
            + x95 * (x139 * x260 + x268 - x7 * (x139 * x261 + x269))
        )
    )
    result[0, 7] = numpy.sum(
        x145
        * (
            x3**2 * (x139 * x273 + x274 + x4 * (x270 - x272))
            + x95 * (x139 * x270 + x276 - x7 * (x139 * x271 + x277))
        )
    )
    result[0, 8] = numpy.sum(
        x145
        * (
            x3**2 * (x139 * x281 + x280 + x4 * (x237 * x99 - x279))
            + x95 * (x139 * x283 + x282 - x7 * (x241 * x99 + x284))
        )
    )
    result[0, 9] = numpy.sum(
        x97
        * (
            x3 * (x3 * (x290 + x291) + x4 * (x285 * x3 - x288))
            + x95 * (x139 * x285 + x293 - x7 * (x139 * x286 + x294))
        )
    )
    result[1, 0] = numpy.sum(
        0.5
        * x306
        * (
            x3 * (2.0 * x3 * (x297 + x298) + x300 + 2.0 * x303 * x68)
            + 2.0 * x4 * (x295 * x41 - x295 * x7 * (x28 + x37) + x296 * x68)
            + x68 * (2.0 * x3 * x303 + x304 + 2.0 * x71 * (x295 * x75 + x305))
        )
    )
    result[1, 1] = numpy.sum(
        0.5
        * x332
        * (
            x3 * (2.0 * x3 * (x3 * x320 + x325) + x326 + 2.0 * x330 * x71)
            + 2.0 * x4 * (x3 * x308 + x312 - x7 * (x3 * x314 + x318))
            + x71 * (2.0 * x3 * x329 + 2.0 * x3 * x330 + x331)
        )
    )
    result[1, 2] = numpy.sum(
        0.5
        * x332
        * (
            2.0 * x295 * x4 * (x158 + x159 - x7 * (x151 + x152))
            + x3 * (2.0 * x295 * x3 * (x164 + x167) + x333 + 2.0 * x335 * x71)
            + x71 * (2.0 * x176 * x295 + 2.0 * x3 * x335 + x336)
        )
    )
    result[1, 3] = numpy.sum(
        0.5
        * x332
        * (
            x29 * (2.0 * x139 * x348 + x349)
            + x3 * (2.0 * x3 * x350 + 2.0 * x3 * (x3 * x353 + x350) + x354)
            + 2.0 * x4 * (x3 * x341 + x339 - x7 * (x3 * x346 + x343))
        )
    )
    result[1, 4] = numpy.sum(
        0.5
        * x373
        * (
            x29 * (2.0 * x139 * x365 + x366)
            + x3 * (2.0 * x3 * x367 + 2.0 * x3 * (x3 * x371 + x367) + x372)
            + 2.0 * x4 * (x3 * x360 + x357 - x7 * (x3 * x364 + x362))
        )
    )
    result[1, 5] = numpy.sum(
        0.5
        * x332
        * (
            x29 * (2.0 * x241 * x295 + x376)
            + x3 * (2.0 * x256 * x295 + 2.0 * x3 * (x253 * x295 + x377) + x378)
            + 2.0 * x4 * (x249 * x295 + x374 - x7 * (x245 * x295 + x375))
        )
    )
    result[1, 6] = numpy.sum(
        0.5 * x3 * x306 * (2.0 * x139 * x388 + x389 + 2.0 * x4 * (x381 - x385))
    )
    result[1, 7] = numpy.sum(
        0.5 * x3 * x332 * (2.0 * x139 * x396 + x397 + 2.0 * x4 * (x391 - x394))
    )
    result[1, 8] = numpy.sum(
        0.5 * x3 * x332 * (2.0 * x139 * x401 + 2.0 * x4 * (x398 - x400) + x402)
    )
    result[1, 9] = numpy.sum(
        0.5
        * x306
        * (2.0 * x295 * x4 * (x285 * x3 - x288) + x3 * (2.0 * x290 * x295 + x403))
    )
    result[2, 0] = numpy.sum(
        x306
        * (
            x3 * (x3 * x404 * (x56 + x58) + x406 + x407 * x68)
            + x4 * x404 * (x41 + x49 - x7 * (x28 + x37))
            + x68 * (x3 * x407 + x409 + x71 * (x404 * x75 + x410))
        )
    )
    result[2, 1] = numpy.sum(
        x332
        * (
            x3 * (x3 * (x113 * x404 + x413) + x415 + x417 * x71)
            + x4 * (x108 * x404 + x412 - x7 * (x102 * x404 + x411))
            + x71 * (x127 * x404 + x3 * x417 + x419)
        )
    )
    result[2, 2] = numpy.sum(
        x332
        * (
            x3 * (x3 * (x430 + x434) + x437 + x441 * x71)
            + x4 * (x3 * x420 + x423 - x7 * (x3 * x424 + x428))
            + x71 * (x3 * x439 + x3 * x441 + x443)
        )
    )
    result[2, 3] = numpy.sum(
        x332
        * (
            x29 * (x199 * x404 + x447)
            + x3 * (x216 * x404 + x3 * (x212 * x404 + x448) + x450)
            + x4 * (x207 * x404 + x444 - x7 * (x203 * x404 + x445))
        )
    )
    result[2, 4] = numpy.sum(
        x373
        * (
            x29 * (x139 * x459 + x461)
            + x3 * (x3 * x462 + x3 * (x3 * x464 + x462) + x466)
            + x4 * (x3 * x454 + x452 - x7 * (x3 * x458 + x456))
        )
    )
    result[2, 5] = numpy.sum(
        x332
        * (
            x29 * (x139 * x474 + x476)
            + x3 * (x3 * (x477 + x480) + x483 + x484)
            + x4 * (x3 * x469 + x468 - x7 * (x3 * x473 + x471))
        )
    )
    result[2, 6] = numpy.sum(
        x306 * (x3 * (x265 * x404 + x486) + x4 * x404 * (x260 * x3 - x263))
    )
    result[2, 7] = numpy.sum(x3 * x332 * (x139 * x490 + x4 * (x487 - x489) + x492))
    result[2, 8] = numpy.sum(x3 * x332 * (x139 * x496 + x4 * (x493 - x495) + x498))
    result[2, 9] = numpy.sum(x3 * x306 * (x4 * (x499 - x501) + x503 + x505))
    result[3, 0] = numpy.sum(
        x514
        * (x3 * (x3 * x508 + 3.0 * x511) + x512 * x60 + x68 * (x3 * x510 + 2.0 * x513))
    )
    result[3, 1] = numpy.sum(
        x522 * (x3 * (x3 * x516 + x519) + x520 * x60 + x71 * (x3 * x518 + x521))
    )
    result[3, 2] = numpy.sum(
        x522 * (x3 * (x3 * x523 + 2.0 * x525) + x526 * x60 + x71 * (x3 * x524 + x527))
    )
    result[3, 3] = numpy.sum(x522 * (x3 * x529 + x3 * (x3 * x530 + x529) + x531 * x60))
    result[3, 4] = numpy.sum(x536 * (x3 * x533 + x3 * (x3 * x534 + x533) + x535 * x60))
    result[3, 5] = numpy.sum(x522 * (x3 * x539 + x3 * (x3 * x540 + x539) + x541 * x60))
    result[3, 6] = numpy.sum(x514 * (x139 * x542 + x543 * x60))
    result[3, 7] = numpy.sum(x522 * (x139 * x544 + x545 * x60))
    result[3, 8] = numpy.sum(x522 * (x139 * x546 + x547 * x60))
    result[3, 9] = numpy.sum(x514 * (x139 * x548 + x549 * x60))
    result[4, 0] = numpy.sum(
        0.5
        * x332
        * (2.0 * x3 * x404 * (x297 + x298) + 2.0 * x404 * x68 * (x301 + x302) + x550)
    )
    result[4, 1] = numpy.sum(
        0.5
        * x560
        * (2.0 * x3 * (x3 * x552 + x555) + x558 + 2.0 * x71 * (x3 * x554 + x559))
    )
    result[4, 2] = numpy.sum(
        0.5
        * x560
        * (2.0 * x295 * x3 * (x430 + x434) + x561 + 2.0 * x71 * (x295 * x440 + x562))
    )
    result[4, 3] = numpy.sum(0.5 * x560 * (2.0 * x3**2 * x565 + 4.0 * x3 * x564 + x566))
    result[4, 4] = numpy.sum(0.5 * x571 * (2.0 * x3**2 * x569 + 4.0 * x3 * x568 + x570))
    result[4, 5] = numpy.sum(
        0.5 * x560 * (2.0 * x295 * x484 + 2.0 * x3 * (x295 * x480 + x572) + x573)
    )
    result[4, 6] = numpy.sum(0.5 * x332 * (2.0 * x139 * x574 + x575))
    result[4, 7] = numpy.sum(0.5 * x560 * (2.0 * x139 * x577 + x580))
    result[4, 8] = numpy.sum(0.5 * x560 * (2.0 * x139 * x581 + x582))
    result[4, 9] = numpy.sum(0.5 * x332 * (2.0 * x295 * x503 + x583))
    result[5, 0] = numpy.sum(
        x514 * (x3 * (x586 + x589) + x591 + x68 * (x3 * x587 + 2.0 * x592))
    )
    result[5, 1] = numpy.sum(
        x522 * (x3 * (x3 * x593 + x596) + x598 + x71 * (x3 * x594 + x599))
    )
    result[5, 2] = numpy.sum(
        x522 * (x3 * (x601 + x604) + x606 + x71 * (x3 * x602 + x607))
    )
    result[5, 3] = numpy.sum(x522 * (x3 * x610 + x3 * (x3 * x611 + x610) + x613))
    result[5, 4] = numpy.sum(x536 * (x3 * x615 + x3 * (x3 * x616 + x615) + x618))
    result[5, 5] = numpy.sum(x522 * (x3 * x620 + x3 * (x620 + x622) + x624))
    result[5, 6] = numpy.sum(x514 * (x139 * x625 + x627))
    result[5, 7] = numpy.sum(x522 * (x139 * x628 + x630))
    result[5, 8] = numpy.sum(x522 * (x139 * x631 + x633))
    result[5, 9] = numpy.sum(x514 * (x139 * x634 + x636))
    result[6, 0] = numpy.sum(x306 * (x3 * x637 + 3.0 * x638))
    result[6, 1] = numpy.sum(x332 * (x3 * x639 + x640))
    result[6, 2] = numpy.sum(x332 * (x3 * x641 + 2.0 * x642))
    result[6, 3] = numpy.sum(x332 * (x3 * x644 + x643))
    result[6, 4] = numpy.sum(x373 * (x3 * x646 + x645))
    result[6, 5] = numpy.sum(x332 * (x3 * x648 + x647))
    result[6, 6] = numpy.sum(x649 * x650)
    result[6, 7] = numpy.sum(x651 * x652)
    result[6, 8] = numpy.sum(x652 * x653)
    result[6, 9] = numpy.sum(x650 * x654)
    result[7, 0] = numpy.sum(x332 * (x3 * x655 + 3.0 * x656))
    result[7, 1] = numpy.sum(x560 * (x3 * x657 + x658))
    result[7, 2] = numpy.sum(x560 * (x3 * x659 + 2.0 * x660))
    result[7, 3] = numpy.sum(x560 * (x3 * x662 + x661))
    result[7, 4] = numpy.sum(x571 * (x3 * x664 + x663))
    result[7, 5] = numpy.sum(x560 * (x3 * x666 + x665))
    result[7, 6] = numpy.sum(x652 * x667)
    result[7, 7] = numpy.sum(x668 * x669)
    result[7, 8] = numpy.sum(x669 * x670)
    result[7, 9] = numpy.sum(x652 * x671)
    result[8, 0] = numpy.sum(x295 * x332 * (x586 + x589))
    result[8, 1] = numpy.sum(x560 * (x3 * x672 + x673))
    result[8, 2] = numpy.sum(x295 * x560 * (x601 + x604))
    result[8, 3] = numpy.sum(x560 * (x3 * x675 + x674))
    result[8, 4] = numpy.sum(x571 * (x3 * x677 + x676))
    result[8, 5] = numpy.sum(x560 * (x295 * x622 + x678))
    result[8, 6] = numpy.sum(x652 * x679)
    result[8, 7] = numpy.sum(x669 * x681)
    result[8, 8] = numpy.sum(x669 * x682)
    result[8, 9] = numpy.sum(x295 * x634 * x652)
    result[9, 0] = numpy.sum(x306 * (x3 * x683 + 3.0 * x684))
    result[9, 1] = numpy.sum(x332 * (x3 * x685 + x687))
    result[9, 2] = numpy.sum(x332 * (x3 * x688 + x690))
    result[9, 3] = numpy.sum(x332 * (x3 * x692 + x691))
    result[9, 4] = numpy.sum(x373 * (x3 * x694 + x693))
    result[9, 5] = numpy.sum(x332 * (x3 * x696 + x695))
    result[9, 6] = numpy.sum(x650 * x697)
    result[9, 7] = numpy.sum(x652 * x698)
    result[9, 8] = numpy.sum(x652 * x699)
    result[9, 9] = numpy.sum(x650 * x700)
    result[10, 0] = numpy.sum(x97 * (x295 * x637 + x512 * x95))
    result[10, 1] = numpy.sum(x145 * (x295 * x639 + x520 * x95 + x638))
    result[10, 2] = numpy.sum(x145 * (x295 * x641 + x526 * x95))
    result[10, 3] = numpy.sum(x145 * (x295 * x644 + x531 * x95 + x640))
    result[10, 4] = numpy.sum(x234 * (x295 * x646 + x535 * x95 + x642))
    result[10, 5] = numpy.sum(x145 * (x295 * x648 + x541 * x95))
    result[10, 6] = numpy.sum(x97 * (x295 * x649 + x543 * x95 + 3.0 * x643))
    result[10, 7] = numpy.sum(x145 * (x295 * x651 + x545 * x95 + 2.0 * x645))
    result[10, 8] = numpy.sum(x145 * (x295 * x653 + x547 * x95 + x647))
    result[10, 9] = numpy.sum(x97 * (x295 * x654 + x549 * x95))
    result[11, 0] = numpy.sum(x306 * (x295 * x655 + x550))
    result[11, 1] = numpy.sum(x332 * (x295 * x657 + x558 + x656))
    result[11, 2] = numpy.sum(x332 * (x295 * x659 + x561))
    result[11, 3] = numpy.sum(x332 * (x295 * x662 + x566 + x658))
    result[11, 4] = numpy.sum(x373 * (x295 * x664 + x570 + x660))
    result[11, 5] = numpy.sum(x332 * (x295 * x666 + x573))
    result[11, 6] = numpy.sum(x306 * (x295 * x667 + x575 + 3.0 * x661))
    result[11, 7] = numpy.sum(x332 * (x295 * x668 + x580 + 2.0 * x663))
    result[11, 8] = numpy.sum(x332 * (x295 * x670 + x582 + x665))
    result[11, 9] = numpy.sum(x306 * (x295 * x671 + x583))
    result[12, 0] = numpy.sum(x514 * (x506 * x585 + x591))
    result[12, 1] = numpy.sum(x522 * (x295 * x588 + x295 * x672 + x598))
    result[12, 2] = numpy.sum(x522 * (x506 * x600 + x606))
    result[12, 3] = numpy.sum(x522 * (x295 * x675 + x613 + x673))
    result[12, 4] = numpy.sum(x536 * (x295 * x603 + x295 * x677 + x618))
    result[12, 5] = numpy.sum(x522 * (x506 * x621 + x624))
    result[12, 6] = numpy.sum(x514 * (x295 * x679 + x627 + 3.0 * x674))
    result[12, 7] = numpy.sum(x522 * (x295 * x681 + x630 + 2.0 * x676))
    result[12, 8] = numpy.sum(x522 * (x295 * x682 + x633 + x678))
    result[12, 9] = numpy.sum(x514 * (x506 * x634 + x636))
    result[13, 0] = numpy.sum(x683 * x701)
    result[13, 1] = numpy.sum(x332 * (x295 * x685 + x684))
    result[13, 2] = numpy.sum(x688 * x702)
    result[13, 3] = numpy.sum(x332 * (x295 * x692 + x687))
    result[13, 4] = numpy.sum(x373 * (x295 * x694 + x689))
    result[13, 5] = numpy.sum(x696 * x702)
    result[13, 6] = numpy.sum(x306 * (x295 * x697 + 3.0 * x691))
    result[13, 7] = numpy.sum(x332 * (x295 * x698 + x703))
    result[13, 8] = numpy.sum(x332 * (x295 * x699 + x695))
    result[13, 9] = numpy.sum(x700 * x701)
    result[14, 0] = numpy.sum(x97 * (x404 * x683 + x590 * x95))
    result[14, 1] = numpy.sum(x145 * (x404 * x685 + x597 * x95))
    result[14, 2] = numpy.sum(x145 * (x404 * x688 + x605 * x95 + x684))
    result[14, 3] = numpy.sum(x145 * (x404 * x692 + x612 * x95))
    result[14, 4] = numpy.sum(x234 * (x404 * x694 + x617 * x95 + x686))
    result[14, 5] = numpy.sum(x145 * (x404 * x696 + x623 * x95 + x690))
    result[14, 6] = numpy.sum(x97 * (x404 * x697 + x626 * x95))
    result[14, 7] = numpy.sum(x145 * (x404 * x698 + x629 * x95 + x691))
    result[14, 8] = numpy.sum(x145 * (x404 * x699 + x632 * x95 + x703))
    result[14, 9] = numpy.sum(x97 * (x404 * x700 + x635 * x95 + 3.0 * x695))
    return result


def _2center2el3d_44(ax, da, A, bx, db, B):
    """Cartesian (g|g) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 15), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = ax ** (-1.0)
    x5 = -x2 - B[0]
    x6 = bx ** (-1.0)
    x7 = ax * x1
    x8 = bx * x7 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x9 = boys(6, x8)
    x10 = x0 ** (-1.5)
    x11 = 17.4934183276249
    x12 = x11 * x6
    x13 = 2.0 * x12
    x14 = x10 * x13
    x15 = x14 * x9
    x16 = x15 * x5
    x17 = x0 ** (-0.5)
    x18 = boys(5, x8)
    x19 = 0.5 * x6
    x20 = x19 * (2.0 * x11 * x17 * x18 * x4 * x6 - x15)
    x21 = boys(7, x8)
    x22 = x5**2
    x23 = x17 * x4
    x24 = x13 * x23
    x25 = x22 * x24
    x26 = x21 * x25
    x27 = x20 + x26
    x28 = x27 * x5 + x6 * (2.0 * x11 * x17 * x18 * x4 * x5 * x6 - x16)
    x29 = x14 * x18
    x30 = boys(4, x8)
    x31 = x19 * (2.0 * x11 * x17 * x30 * x4 * x6 - x29)
    x32 = x24 * x9
    x33 = x22 * x32
    x34 = x31 + x33
    x35 = x14 * x30
    x36 = boys(3, x8)
    x37 = x19 * (2.0 * x11 * x17 * x36 * x4 * x6 - x35)
    x38 = x18 * x24
    x39 = x22 * x38
    x40 = x37 + x39
    x41 = 1.5 * x6
    x42 = x28 * x5 - x41 * (x34 * x7 - x40)
    x43 = x3 * x42
    x44 = 0.5 / (ax + bx)
    x45 = x29 * x5
    x46 = x34 * x5 + x6 * (2.0 * x11 * x17 * x30 * x4 * x5 * x6 - x45)
    x47 = x44 * x46
    x48 = 4.0 * x47
    x49 = x43 + x48
    x50 = bx * x1
    x51 = x14 * x36
    x52 = boys(2, x8)
    x53 = x19 * (2.0 * x11 * x17 * x4 * x52 * x6 - x51)
    x54 = x24 * x30
    x55 = x22 * x54
    x56 = x53 + x55
    x57 = -x41 * (x40 * x7 - x56) + x46 * x5
    x58 = x3 * x57
    x59 = x35 * x5
    x60 = x40 * x5 + x6 * (2.0 * x11 * x17 * x36 * x4 * x5 * x6 - x59)
    x61 = x44 * x60
    x62 = 4.0 * x61
    x63 = x58 + x62
    x64 = x14 * x21
    x65 = x5 * x64
    x66 = x19 * (2.0 * x11 * x17 * x4 * x6 * x9 - x64)
    x67 = boys(8, x8)
    x68 = x25 * x67
    x69 = -x41 * (x27 * x7 - x34) + x5 * (
        x5 * (x66 + x68) + x6 * (2.0 * x11 * x17 * x4 * x5 * x6 * x9 - x65)
    )
    x70 = x3 * x69
    x71 = x28 * x44
    x72 = 4.0 * x71
    x73 = x42 * x50
    x74 = 0.5 * x4
    x75 = x74 * (x57 - x73)
    x76 = x28 * x3
    x77 = x34 * x44
    x78 = 3.0 * x77
    x79 = x76 + x78
    x80 = 4.0 * x44
    x81 = x46 * x50
    x82 = x74 * (x60 - x81)
    x83 = x3 * x34
    x84 = x18 * x5
    x85 = x12 * x23
    x86 = x80 * x85
    x87 = x84 * x86
    x88 = x83 + x87
    x89 = 3.0 * x44
    x90 = x3 * x79 + x82 + x88 * x89
    x91 = x40 * x44
    x92 = x3 * x46 + 3.0 * x91
    x93 = x44 * x56
    x94 = x3 * x60 + 3.0 * x93
    x95 = x74 * (-x40 * x50 + x56)
    x96 = 2.0 * x44
    x97 = x85 * x96
    x98 = x30 * x97
    x99 = x3 * x5
    x100 = x5 * x51
    x101 = boys(1, x8)
    x102 = x19 * (-x101 * x14 + 2.0 * x11 * x17 * x4 * x6 * boys(0, x8))
    x103 = x19 * (2.0 * x101 * x11 * x17 * x4 * x6 - x14 * x52)
    x104 = x24 * x36
    x105 = x103 + x104 * x22
    x106 = x24 * x52
    x107 = x41 * (x105 - x56 * x7) + x5 * x60
    x108 = x74 * (
        -x107 * x50
        + x41 * (x102 - x105 * x7 + x106 * x22)
        + x5 * (x5 * x56 - x6 * (x100 - 2.0 * x11 * x17 * x4 * x5 * x52 * x6))
    )
    x109 = x74 * (x107 - x50 * x57)
    x110 = 1.5 * x4
    x111 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**5.5)
    x112 = 6.89597470414309 * x111
    x113 = -x1 * (ax * A[1] + bx * B[1])
    x114 = -x113 - B[1]
    x115 = x114 * x5
    x116 = x114 * x15
    x117 = x6 * (2.0 * x11 * x114 * x17 * x18 * x4 * x6 - x116)
    x118 = x114 * x26
    x119 = 0.5 * x117 + x118
    x120 = x119 * x5 + x6 * (2.0 * x11 * x114 * x17 * x18 * x4 * x5 * x6 - x115 * x15)
    x121 = x120 * x3
    x122 = x114 * x29
    x123 = x6 * (2.0 * x11 * x114 * x17 * x30 * x4 * x6 - x122)
    x124 = x114 * x33
    x125 = 0.5 * x123 + x124
    x126 = x125 * x44
    x127 = 3.0 * x126
    x128 = x121 + x127
    x129 = -2.0 * x11 * x114 * x17 * x30 * x4 * x5 * x6
    x130 = x125 * x5 - x6 * (x115 * x29 + x129)
    x131 = x130 * x3
    x132 = x114 * x35
    x133 = x6 * (2.0 * x11 * x114 * x17 * x36 * x4 * x6 - x132)
    x134 = x114 * x39
    x135 = 0.5 * x133 + x134
    x136 = x135 * x44
    x137 = x131 + 3.0 * x136
    x138 = x114 * x64
    x139 = x6 * (2.0 * x11 * x114 * x17 * x4 * x6 * x9 - x138)
    x140 = x114 * x68
    x141 = 0.5 * x5 * (x139 + 2.0 * x140) + x6 * (
        2.0 * x11 * x114 * x17 * x4 * x5 * x6 * x9 - x115 * x64
    )
    x142 = x141 * x3
    x143 = x119 * x44
    x144 = 3.0 * x143
    x145 = x120 * x50
    x146 = x74 * (x130 - x145)
    x147 = x119 * x3
    x148 = x115 * x9
    x149 = x148 * x86
    x150 = x147 + x149
    x151 = x125 * x50
    x152 = x74 * (x135 - x151)
    x153 = x114 * x18
    x154 = x153 * x97
    x155 = x114 * x99
    x156 = x155 * x32
    x157 = x154 + x156
    x158 = x150 * x3 + x152 + x157 * x96
    x159 = x115 * x18
    x160 = x159 * x86
    x161 = x125 * x3 + x160
    x162 = x30 * x86
    x163 = x115 * x162
    x164 = x135 * x3 + x163
    x165 = 2.0 * x10 * x11 * x4
    x166 = -x74 * (x129 + x159 * x165)
    x167 = x114 * x51
    x168 = x6 * (2.0 * x11 * x114 * x17 * x4 * x52 * x6 - x167)
    x169 = x114 * x55
    x170 = x135 * x5 + x6 * (2.0 * x11 * x114 * x17 * x36 * x4 * x5 * x6 - x115 * x35)
    x171 = (
        0.5
        * x74
        * (
            -2.0 * x170 * x50
            + x5 * (x168 + 2.0 * x169)
            + 2.0 * x6 * (2.0 * x11 * x114 * x17 * x4 * x5 * x52 * x6 - x115 * x51)
        )
    )
    x172 = x74 * (-x130 * x50 + x170)
    x173 = 18.2450341145548 * x111
    x174 = -x1 * (ax * A[2] + bx * B[2])
    x175 = -x174 - B[2]
    x176 = x175 * x6 * (2.0 * x11 * x17 * x18 * x4 * x6 - x15)
    x177 = 0.5 * x176
    x178 = x175 * x26 + x177
    x179 = x175 * x6 * (2.0 * x11 * x17 * x18 * x4 * x5 * x6 - x16) + x178 * x5
    x180 = x179 * x3
    x181 = x175 * x6 * (2.0 * x11 * x17 * x30 * x4 * x6 - x29)
    x182 = 0.5 * x181
    x183 = x175 * x33 + x182
    x184 = x183 * x44
    x185 = 3.0 * x184
    x186 = x180 + x185
    x187 = -2.0 * x11 * x17 * x175 * x30 * x4 * x5 * x6
    x188 = x183 * x5 - x6 * (x175 * x45 + x187)
    x189 = x188 * x3
    x190 = x175 * x6 * (2.0 * x11 * x17 * x36 * x4 * x6 - x35)
    x191 = 0.5 * x190
    x192 = x175 * x39 + x191
    x193 = x192 * x44
    x194 = x189 + 3.0 * x193
    x195 = x175 * x6 * (2.0 * x11 * x17 * x4 * x6 * x9 - x64)
    x196 = 0.5 * x195
    x197 = x175 * x6 * (2.0 * x11 * x17 * x4 * x5 * x6 * x9 - x65) + x5 * (
        x175 * x68 + x196
    )
    x198 = x197 * x3
    x199 = x178 * x44
    x200 = 3.0 * x199
    x201 = x179 * x50
    x202 = x74 * (x188 - x201)
    x203 = x178 * x3
    x204 = x175 * x5
    x205 = x204 * x9
    x206 = x205 * x86
    x207 = x203 + x206
    x208 = x183 * x50
    x209 = x74 * (x192 - x208)
    x210 = x175 * x18
    x211 = x210 * x97
    x212 = x175 * x32
    x213 = x212 * x99
    x214 = x211 + x213
    x215 = x207 * x3 + x209 + x214 * x96
    x216 = x175 * x87
    x217 = x183 * x3 + x216
    x218 = x162 * x204
    x219 = x192 * x3 + x218
    x220 = x165 * x175
    x221 = -x74 * (x187 + x220 * x84)
    x222 = x175 * x6 * (2.0 * x11 * x17 * x4 * x52 * x6 - x51)
    x223 = 0.5 * x222
    x224 = x175 * x6 * (2.0 * x11 * x17 * x36 * x4 * x5 * x6 - x59) + x192 * x5
    x225 = -x74 * (
        x175 * x6 * (x100 - 2.0 * x11 * x17 * x4 * x5 * x52 * x6)
        + x224 * x50
        - x5 * (x175 * x55 + x223)
    )
    x226 = x74 * (-x188 * x50 + x224)
    x227 = x114**2
    x228 = x227 * x24
    x229 = x21 * x228
    x230 = x20 + x229
    x231 = x227 * x32 + x31
    x232 = x227 * x38 + x37
    x233 = -x231 * x7 + x232
    x234 = x19 * x233 + x22 * x230
    x235 = x234 * x3
    x236 = x231 * x44
    x237 = x236 * x5
    x238 = x235 + 2.0 * x237
    x239 = x227 * x54 + x53
    x240 = -x232 * x7 + x239
    x241 = x19 * x240 + x22 * x231
    x242 = x241 * x3
    x243 = x232 * x5
    x244 = x242 + x243 * x96
    x245 = x228 * x67
    x246 = x245 + x66
    x247 = -x230 * x7 + x231
    x248 = x19 * x247 + x22 * x246
    x249 = x248 * x3
    x250 = x230 * x5
    x251 = x250 * x96
    x252 = x234 * x50
    x253 = x74 * (x241 - x252)
    x254 = x230 * x99
    x255 = x236 + x254
    x256 = x5 * x50
    x257 = x74 * (-x231 * x256 + x232 * x5)
    x258 = x236 * x3
    x259 = x255 * x3 + x257 + x258
    x260 = x232 * x44
    x261 = x231 * x99 + x260
    x262 = x239 * x44
    x263 = x243 * x3 + x262
    x264 = x3**2
    x265 = x74 * (-x232 * x50 + x239)
    x266 = x103 + x104 * x227
    x267 = x102 + x106 * x227 - x266 * x7
    x268 = -x239 * x7 + x266
    x269 = x19 * x268 + x22 * x232
    x270 = x74 * (x19 * x267 + x22 * x239 - x269 * x50)
    x271 = x74 * (-x241 * x50 + x269)
    x272 = 23.5542377588857 * x111
    x273 = x175 * x6 * (2.0 * x11 * x114 * x17 * x18 * x4 * x6 - x116)
    x274 = x118 * x175 + 0.5 * x273
    x275 = x149 * x175
    x276 = x274 * x3 + x275
    x277 = -2.0 * x11 * x114 * x17 * x175 * x30 * x4 * x6
    x278 = -x6 * (x122 * x175 + x277)
    x279 = x124 * x175 + 0.5 * x278
    x280 = x160 * x175
    x281 = x279 * x3 + x280
    x282 = x175 * x6 * (2.0 * x11 * x114 * x17 * x4 * x6 * x9 - x138)
    x283 = x140 * x175 + 0.5 * x282
    x284 = x175 * x21
    x285 = x115 * x284
    x286 = x285 * x86
    x287 = x74 * (-x274 * x50 + x279)
    x288 = x114 * x175
    x289 = x288 * x9
    x290 = x289 * x97
    x291 = x155 * x24 * x284 + x290
    x292 = x74 * (2.0 * x11 * x114 * x17 * x175 * x18 * x4 * x5 * x6 - x148 * x220)
    x293 = x290 * x3 + x291 * x3 + x292
    x294 = x154 * x175 + x155 * x212
    x295 = x175 * x38
    x296 = x155 * x295 + x288 * x98
    x297 = -x74 * (x153 * x220 + x277)
    x298 = x114 * x175 * x32
    x299 = x175 * x6 * (2.0 * x11 * x114 * x17 * x4 * x52 * x6 - x167)
    x300 = x175 * x6 * (2.0 * x11 * x114 * x17 * x36 * x4 * x6 - x132)
    x301 = x134 * x175 + 0.5 * x300
    x302 = 0.5 * x74 * (2.0 * x169 * x175 + x299 - 2.0 * x301 * x50)
    x303 = x74 * (-x279 * x50 + x301)
    x304 = 40.7971365319473 * x111
    x305 = x175**2
    x306 = x24 * x305
    x307 = x20 + x21 * x306
    x308 = x305 * x32 + x31
    x309 = x305 * x38 + x37
    x310 = -x308 * x7 + x309
    x311 = x19 * x310
    x312 = x22 * x307 + x311
    x313 = x3 * x312
    x314 = x308 * x44
    x315 = x314 * x5
    x316 = x313 + 2.0 * x315
    x317 = x305 * x54 + x53
    x318 = -x309 * x7 + x317
    x319 = x19 * x318
    x320 = x22 * x308 + x319
    x321 = x3 * x320
    x322 = x309 * x5
    x323 = x321 + x322 * x96
    x324 = x306 * x67 + x66
    x325 = -x307 * x7 + x308
    x326 = x19 * x325
    x327 = x22 * x324 + x326
    x328 = x3 * x327
    x329 = x307 * x5
    x330 = x312 * x50
    x331 = x74 * (x320 - x330)
    x332 = x307 * x99
    x333 = x314 + x332
    x334 = x74 * (-x256 * x308 + x309 * x5)
    x335 = x3 * x314
    x336 = x3 * x333 + x334 + x335
    x337 = x309 * x44
    x338 = x308 * x99
    x339 = x337 + x338
    x340 = x317 * x44
    x341 = x3 * x322 + x340
    x342 = x74 * (-x309 * x50 + x317)
    x343 = x103 + x104 * x305
    x344 = x102 + x106 * x305 - x343 * x7
    x345 = x19 * x344
    x346 = -x317 * x7 + x343
    x347 = x19 * x346
    x348 = x22 * x309 + x347
    x349 = x74 * (x22 * x317 + x345 - x348 * x50)
    x350 = x74 * (-x320 * x50 + x348)
    x351 = x114 * x232 + x133
    x352 = x114 * x231 + x123
    x353 = x352 * x50
    x354 = x114 * x230 + x117
    x355 = x264 * x354
    x356 = x74 * (x351 - x353)
    x357 = x355 + x356
    x358 = x352 * x44
    x359 = x354 * x99
    x360 = x358 + x359
    x361 = x351 * x44
    x362 = x352 * x5
    x363 = x3 * x362
    x364 = x361 + x363
    x365 = x354 * x44
    x366 = x114 * x246 + x139
    x367 = x366 * x99
    x368 = x74 * (-x256 * x354 + x352 * x5)
    x369 = x3 * x365
    x370 = x114 * x239 + x168
    x371 = x5 * x74 * (-x351 * x50 + x370)
    x372 = x74 * (x351 * x5 - x362 * x50)
    x373 = x191 + x227 * x295
    x374 = x182 + x212 * x227
    x375 = x374 * x50
    x376 = x175 * x229 + x177
    x377 = x74 * (x373 - x375)
    x378 = x264 * x376 + x377
    x379 = x374 * x44
    x380 = x376 * x99 + x379
    x381 = x373 * x44
    x382 = x374 * x5
    x383 = x3 * x382 + x381
    x384 = x376 * x44
    x385 = x175 * x245 + x196
    x386 = x74 * (-x256 * x376 + x374 * x5)
    x387 = x175 * x227 * x54 + x223
    x388 = x5 * x74 * (-x373 * x50 + x387)
    x389 = x74 * (x373 * x5 - x382 * x50)
    x390 = x114 * x50
    x391 = x308 * x390
    x392 = x74 * (x114 * x309 - x391)
    x393 = x114 * x307
    x394 = x264 * x393 + x392
    x395 = x114 * x314
    x396 = x114 * x332 + x395
    x397 = x114 * x337
    x398 = x114 * x338 + x397
    x399 = x393 * x44
    x400 = x74 * (x114 * x308 * x5 - x256 * x393)
    x401 = x114 * x74 * (x317 * x5 - x322 * x50)
    x402 = x115 * x308
    x403 = x74 * (x114 * x309 * x5 - x402 * x50)
    x404 = x175 * x309 + x190
    x405 = x175 * x308 + x181
    x406 = x405 * x50
    x407 = x175 * x307 + x176
    x408 = x264 * x407
    x409 = x74 * (x404 - x406)
    x410 = x408 + x409
    x411 = x405 * x44
    x412 = x407 * x99
    x413 = x411 + x412
    x414 = x404 * x44
    x415 = x405 * x5
    x416 = x3 * x415
    x417 = x414 + x416
    x418 = x407 * x44
    x419 = x175 * x324 + x195
    x420 = x419 * x99
    x421 = x74 * (-x256 * x407 + x405 * x5)
    x422 = x3 * x418
    x423 = x175 * x317 + x222
    x424 = x5 * x74 * (-x404 * x50 + x423)
    x425 = x74 * (x404 * x5 - x415 * x50)
    x426 = x114 * x352 + x240 * x41
    x427 = x114 * x354 + x233 * x41
    x428 = x427 * x50
    x429 = x3 * x428
    x430 = x114 * x366 + x247 * x41
    x431 = x264 * x430
    x432 = x74 * (x426 - x428)
    x433 = x114 * x351 + x268 * x41
    x434 = x74 * (x114 * x370 + x267 * x41 - x433 * x50)
    x435 = x74 * (-x426 * x50 + x433)
    x436 = x114 * x374 + x278
    x437 = x114 * x376 + x273
    x438 = x437 * x50
    x439 = x114 * x385 + x282
    x440 = x74 * (x436 - x438)
    x441 = x114 * x373 + x300
    x442 = x74 * (x114 * x387 + x299 - x441 * x50)
    x443 = x74 * (-x436 * x50 + x441)
    x444 = x227 * x308 + x319
    x445 = x227 * x307 + x311
    x446 = x445 * x50
    x447 = x227 * x324 + x326
    x448 = x74 * (x444 - x446)
    x449 = x227 * x309 + x347
    x450 = x74 * (x227 * x317 + x345 - x449 * x50)
    x451 = x74 * (-x444 * x50 + x449)
    x452 = x390 * x407
    x453 = x74 * (x114 * x405 - x452)
    x454 = x114 * x419
    x455 = x114 * x74 * (-x404 * x50 + x423)
    x456 = x114 * x405
    x457 = x74 * (x114 * x404 - x456 * x50)
    x458 = x175 * x405 + x318 * x41
    x459 = x175 * x407 + x310 * x41
    x460 = x459 * x50
    x461 = x3 * x460
    x462 = x175 * x419 + x325 * x41
    x463 = x264 * x462
    x464 = x74 * (x458 - x460)
    x465 = x175 * x404 + x346 * x41
    x466 = x74 * (x175 * x423 + x344 * x41 - x465 * x50)
    x467 = x74 * (-x458 * x50 + x465)
    x468 = -x113 - A[1]
    x469 = x468 * x60
    x470 = x468 * x70
    x471 = x468 * x72
    x472 = x468 * x73
    x473 = x4 * (x468 * x57 - x472)
    x474 = x468 * x76
    x475 = x468 * x78
    x476 = x474 + x475
    x477 = x4 * x468 * (x60 - x81)
    x478 = x468 * x5
    x479 = x18 * x86
    x480 = x130 * x468
    x481 = x480 + x61
    x482 = x135 * x468 + x93
    x483 = x120 * x468
    x484 = x47 + x483
    x485 = x125 * x468 + x91
    x486 = x141 * x468
    x487 = x486 + x71
    x488 = x119 * x468
    x489 = x488 + x77
    x490 = x4 * (x481 - x484 * x50)
    x491 = x84 * x97
    x492 = x115 * x32
    x493 = x468 * x492
    x494 = x491 + x493
    x495 = x494 * x96
    x496 = x3 * x489 + x495
    x497 = x4 * (x482 - x485 * x50)
    x498 = x114 * x468
    x499 = x44 * (x38 * x498 + x98)
    x500 = 48.2718229290016 * x111
    x501 = x192 * x468
    x502 = x4 * x468 * (x188 - x201)
    x503 = x175 * x478 * x9
    x504 = x203 * x468 + x503 * x86
    x505 = x4 * x468 * (x192 - x208)
    x506 = x211 * x468
    x507 = x241 * x468
    x508 = 2.0 * x136 + x507
    x509 = x163 + x243 * x468
    x510 = x234 * x468
    x511 = 2.0 * x126
    x512 = x510 + x511
    x513 = x231 * x468
    x514 = x160 + x5 * x513
    x515 = x248 * x468
    x516 = 2.0 * x143
    x517 = x515 + x516
    x518 = x230 * x478
    x519 = x149 + x518
    x520 = x4 * (-x50 * x512 + x508)
    x521 = x153 * x86
    x522 = x513 + x521
    x523 = x44 * x522
    x524 = x3 * x519 + x523
    x525 = x4 * (-x50 * x514 + x509)
    x526 = 62.3186554316989 * x111
    x527 = x193 + x279 * x468
    x528 = x115 * x295 * x468 + x204 * x98
    x529 = x528 * x96
    x530 = x184 + x274 * x468
    x531 = x115 * x212
    x532 = x175 * x491 + x468 * x531
    x533 = x532 * x96
    x534 = x199 + x283 * x468
    x535 = x24 * x285
    x536 = x205 * x97 + x468 * x535
    x537 = x536 * x96
    x538 = x4 * (-x50 * x530 + x527)
    x539 = x211 + x212 * x498
    x540 = x44 * x539
    x541 = x3 * x536 + x540
    x542 = x4 * (-x50 * x532 + x528)
    x543 = 107.939077467081 * x111
    x544 = x314 * x478
    x545 = x322 * x468
    x546 = x307 * x478
    x547 = x4 * x468 * (x320 - x330)
    x548 = x314 * x468
    x549 = x332 * x468 + x548
    x550 = x478 * x50
    x551 = x4 * (-x308 * x550 + x309 * x468 * x5)
    x552 = 3.0 * x262 + x351 * x468
    x553 = x44 * x552
    x554 = x362 * x468
    x555 = x243 * x89 + x554
    x556 = 3.0 * x260 + x352 * x468
    x557 = x44 * x556
    x558 = x354 * x468
    x559 = x5 * x558
    x560 = 3.0 * x236
    x561 = x5 * x560 + x559
    x562 = x558 + x560
    x563 = x4 * (-x50 * x556 + x552)
    x564 = x44 * x562
    x565 = x366 * x478
    x566 = x250 * x89
    x567 = x565 + x566
    x568 = x4 * (-x50 * x561 + x555)
    x569 = x162 * x288
    x570 = x373 * x468 + x569
    x571 = x44 * x570
    x572 = x280 + x382 * x468
    x573 = x175 * x521
    x574 = x374 * x468 + x573
    x575 = x44 * x574
    x576 = x376 * x468
    x577 = x275 + x5 * x576
    x578 = x289 * x86
    x579 = x576 + x578
    x580 = x4 * (-x50 * x574 + x570)
    x581 = x44 * x579
    x582 = x286 + x385 * x478
    x583 = x4 * (-x50 * x577 + x572)
    x584 = x114 * x309
    x585 = x340 + x468 * x584
    x586 = x44 * x585
    x587 = x322 * x44 + x402 * x468
    x588 = x308 * x498 + x337
    x589 = x44 * x588
    x590 = x315 + x393 * x478
    x591 = x314 + x393 * x468
    x592 = x4 * (-x50 * x588 + x585)
    x593 = x44 * x591
    x594 = x115 * x324
    x595 = x329 * x44 + x468 * x594
    x596 = x4 * (-x50 * x590 + x587)
    x597 = x414 * x468
    x598 = x411 * x468
    x599 = x4 * x468 * (x404 - x406)
    x600 = x418 * x468
    x601 = x4 * (x405 * x468 * x5 - x407 * x550)
    x602 = x426 * x468
    x603 = 4.0 * x361
    x604 = x602 + x603
    x605 = x427 * x468
    x606 = 4.0 * x358
    x607 = x605 + x606
    x608 = x50 * x607
    x609 = x430 * x468
    x610 = 4.0 * x365
    x611 = x609 + x610
    x612 = x4 * (x604 - x608)
    x613 = 3.0 * x381 + x436 * x468
    x614 = 3.0 * x379 + x437 * x468
    x615 = x50 * x614
    x616 = 3.0 * x384 + x439 * x468
    x617 = x4 * (x613 - x615)
    x618 = 2.0 * x397 + x444 * x468
    x619 = 2.0 * x395 + x445 * x468
    x620 = x50 * x619
    x621 = x393 * x96 + x447 * x468
    x622 = x4 * (x618 - x620)
    x623 = x414 + x456 * x468
    x624 = x407 * x498 + x411
    x625 = x50 * x624
    x626 = x418 + x454 * x468
    x627 = x4 * (x623 - x625)
    x628 = x4 * x468 * (x458 - x460)
    x629 = -x174 - A[2]
    x630 = x4 * x629 * (x57 - x73)
    x631 = 0.5 * x630
    x632 = x629 * (x76 + x78)
    x633 = x4 * x629 * (x60 - x81)
    x634 = 0.5 * x633
    x635 = x5 * x629
    x636 = x135 * x629
    x637 = x4 * x629 * (x130 - x145)
    x638 = 0.5 * x637
    x639 = x149 * x629
    x640 = x147 * x629 + x639
    x641 = x4 * x629 * (x135 - x151)
    x642 = 0.5 * x641
    x643 = x154 * x629
    x644 = x188 * x629 + x61
    x645 = x192 * x629 + x93
    x646 = x44 * x645
    x647 = x179 * x629 + x47
    x648 = x183 * x629 + x91
    x649 = x44 * x648
    x650 = x197 * x629 + x71
    x651 = x3 * x650
    x652 = x178 * x629 + x77
    x653 = x44 * x652
    x654 = 3.0 * x653
    x655 = x50 * x647
    x656 = x4 * (x644 - x655)
    x657 = 0.5 * x656
    x658 = x3 * x652
    x659 = x212 * x635 + x491
    x660 = x44 * x659
    x661 = 2.0 * x660
    x662 = x658 + x661
    x663 = x4 * (-x50 * x648 + x645)
    x664 = 0.5 * x663
    x665 = x295 * x629
    x666 = x44 * (x665 + x98)
    x667 = x236 * x629
    x668 = x5 * x667
    x669 = x243 * x629
    x670 = x4 * x629 * (x241 - x252)
    x671 = 0.5 * x670
    x672 = x254 * x629 + x667
    x673 = x50 * x635
    x674 = x4 * (-x231 * x673 + x232 * x5 * x629)
    x675 = 0.5 * x674
    x676 = x136 + x279 * x629
    x677 = x115 * (x665 + x98)
    x678 = x677 * x96
    x679 = x126 + x274 * x629
    x680 = x159 * x97 + x531 * x629
    x681 = x680 * x96
    x682 = x143 + x283 * x629
    x683 = x148 * x97
    x684 = x535 * x629 + x683
    x685 = x684 * x96
    x686 = x4 * (-x50 * x679 + x676)
    x687 = 0.5 * x686
    x688 = x154 + x298 * x629
    x689 = x44 * x688
    x690 = x3 * x684 + x689
    x691 = x4 * (-x50 * x680 + x677)
    x692 = 0.5 * x691
    x693 = 2.0 * x193 + x320 * x629
    x694 = x218 + x322 * x629
    x695 = x44 * x694
    x696 = 2.0 * x184 + x312 * x629
    x697 = x308 * x629
    x698 = x216 + x5 * x697
    x699 = x44 * x698
    x700 = 2.0 * x199 + x327 * x629
    x701 = x3 * x700
    x702 = x206 + x329 * x629
    x703 = x44 * x702
    x704 = 2.0 * x703
    x705 = x50 * x696
    x706 = x4 * (x693 - x705)
    x707 = 0.5 * x706
    x708 = x210 * x86 + x697
    x709 = x44 * x708
    x710 = x3 * x702
    x711 = x709 + x710
    x712 = x4 * (-x50 * x698 + x694)
    x713 = 0.5 * x712
    x714 = x361 * x629
    x715 = x358 * x629
    x716 = x4 * x629 * (x351 - x353)
    x717 = 0.5 * x716
    x718 = x365 * x629
    x719 = x4 * (x352 * x5 * x629 - x354 * x673)
    x720 = 0.5 * x719
    x721 = x262 + x373 * x629
    x722 = x44 * x721
    x723 = x243 * x44 + x382 * x629
    x724 = x260 + x374 * x629
    x725 = x44 * x724
    x726 = x376 * x629
    x727 = x237 + x5 * x726
    x728 = x236 + x726
    x729 = x4 * (-x50 * x724 + x721)
    x730 = 0.5 * x729
    x731 = x44 * x728
    x732 = x250 * x44
    x733 = x385 * x635 + x732
    x734 = x4 * (-x50 * x727 + x723)
    x735 = 0.5 * x734
    x736 = x569 + x584 * x629
    x737 = x44 * x736
    x738 = x115 * x697 + x280
    x739 = x114 * x697 + x573
    x740 = x44 * x739
    x741 = x275 + x393 * x635
    x742 = x393 * x629 + x578
    x743 = x4 * (-x50 * x739 + x736)
    x744 = 0.5 * x743
    x745 = x44 * x742
    x746 = x286 + x594 * x629
    x747 = x4 * (-x50 * x741 + x738)
    x748 = 0.5 * x747
    x749 = 3.0 * x340 + x404 * x629
    x750 = x44 * x749
    x751 = x322 * x89 + x415 * x629
    x752 = 3.0 * x337 + x405 * x629
    x753 = x44 * x752
    x754 = x407 * x629
    x755 = 3.0 * x315 + x5 * x754
    x756 = 3.0 * x314 + x754
    x757 = x4 * (-x50 * x752 + x749)
    x758 = 0.5 * x757
    x759 = x44 * x756
    x760 = x329 * x89 + x419 * x635
    x761 = x3 * x760
    x762 = x50 * x755
    x763 = x4 * (x751 - x762)
    x764 = 0.5 * x763
    x765 = x3 * x759
    x766 = x4 * x629 * (x426 - x428)
    x767 = 0.5 * x766
    x768 = x361 + x436 * x629
    x769 = x358 + x437 * x629
    x770 = x50 * x769
    x771 = x365 + x439 * x629
    x772 = x4 * (x768 - x770)
    x773 = 0.5 * x772
    x774 = 2.0 * x381 + x444 * x629
    x775 = 2.0 * x379 + x445 * x629
    x776 = x50 * x775
    x777 = 2.0 * x384 + x447 * x629
    x778 = x4 * (x774 - x776)
    x779 = 0.5 * x778
    x780 = 3.0 * x397 + x456 * x629
    x781 = x114 * x754 + 3.0 * x395
    x782 = x50 * x781
    x783 = x393 * x89 + x454 * x629
    x784 = x4 * (x780 - x782)
    x785 = 0.5 * x784
    x786 = 4.0 * x414 + x458 * x629
    x787 = 4.0 * x411 + x459 * x629
    x788 = x50 * x787
    x789 = 4.0 * x418 + x462 * x629
    x790 = x264 * x789
    x791 = x4 * (x786 - x788)
    x792 = 0.5 * x791
    x793 = x468**2
    x794 = x69 * x793
    x795 = x75 + x794
    x796 = x28 * x793
    x797 = x796 + x82
    x798 = x44 * x797
    x799 = x108 - x50 * (x109 + x42 * x793) + x57 * x793
    x800 = x44 * (x34 * x793 + x95)
    x801 = x468 * x71
    x802 = x146 + x468 * x487 + x801
    x803 = x468 * x77
    x804 = x152 + x468 * x489 + x803
    x805 = x171 + x44 * x469 + x468 * x481 - x50 * (x172 + x468 * x47 + x468 * x484)
    x806 = x18 * x97
    x807 = x96 * (x166 + x468 * x494 + x478 * x806)
    x808 = x197 * x793 + x202
    x809 = x178 * x793 + x209
    x810 = x44 * x809
    x811 = x188 * x793 + x225 - x50 * (x179 * x793 + x226)
    x812 = x5 * x793
    x813 = x44 * (x212 * x812 + x221)
    x814 = x253 + x468 * x517 + x489 * x96
    x815 = x257 + x468 * x519 + x495
    x816 = x270 + x468 * x508 + x482 * x96 - x50 * (x271 + x468 * x512 + x485 * x96)
    x817 = x44 * (x265 + x468 * x522 + 2.0 * x499)
    x818 = 80.4530382150027 * x111
    x819 = x199 * x468 + x287 + x468 * x534
    x820 = x292 + x468 * x536 + x503 * x97
    x821 = x820 * x96
    x822 = x302 + x44 * x501 + x468 * x527 - x50 * (x184 * x468 + x303 + x468 * x530)
    x823 = x44 * (x297 + x468 * x539 + x506)
    x824 = 139.348749811665 * x111
    x825 = x327 * x793 + x331
    x826 = x329 * x793 + x334
    x827 = x44 * x826
    x828 = x320 * x793 + x349 - x50 * (x312 * x793 + x350)
    x829 = x44 * (x308 * x793 + x342)
    x830 = x356 + x468 * x562 + 3.0 * x523
    x831 = x44 * x830
    x832 = x368 + x468 * x567 + x519 * x89
    x833 = x371 + x468 * x555 - x50 * (x372 + x468 * x561 + x514 * x89) + x509 * x89
    x834 = x377 + x468 * x579 + 2.0 * x540
    x835 = x44 * x834
    x836 = x386 + x468 * x582 + x537
    x837 = x388 + x468 * x572 - x50 * (x389 + x468 * x577 + x533) + x529
    x838 = x392 + x468 * x591 + x548
    x839 = x44 * x838
    x840 = x400 + x44 * x546 + x468 * x595
    x841 = x401 + x44 * x545 + x468 * x587 - x50 * (x403 + x468 * x590 + x544)
    x842 = x407 * x793
    x843 = x409 + x842
    x844 = x44 * x843
    x845 = x419 * x812 + x421
    x846 = x415 * x793 + x424 - x50 * (x425 + x5 * x842)
    x847 = x432 + x468 * x611 + 4.0 * x564
    x848 = x434 + x468 * x604 - x50 * (x435 + x468 * x607 + 4.0 * x557) + 4.0 * x553
    x849 = x440 + x468 * x616 + 3.0 * x581
    x850 = x442 + x468 * x613 - x50 * (x443 + x468 * x614 + 3.0 * x575) + 3.0 * x571
    x851 = x448 + x468 * x621 + 2.0 * x593
    x852 = x450 + x468 * x618 - x50 * (x451 + x468 * x619 + 2.0 * x589) + 2.0 * x586
    x853 = x453 + x468 * x626 + x600
    x854 = x455 + x468 * x623 - x50 * (x457 + x468 * x624 + x598) + x597
    x855 = x462 * x793 + x464
    x856 = x458 * x793 + x466 - x50 * (x459 * x793 + x467)
    x857 = x4 * x629 * (x468 * x57 - x472)
    x858 = x629 * x71
    x859 = x486 * x629 + x858
    x860 = x629 * x77
    x861 = x488 * x629 + x860
    x862 = x61 * x629
    x863 = x47 * x629
    x864 = x4 * (x480 * x629 - x50 * (x483 * x629 + x863) + x862)
    x865 = x635 * x806
    x866 = x96 * (x493 * x629 + x865)
    x867 = x4 * x468 * (x644 - x655)
    x868 = x629 * (x515 + x516)
    x869 = x518 * x629 + x639
    x870 = x4 * (-x50 * x629 * (x510 + x511) + x507 * x629 + x636 * x96)
    x871 = x44 * x629 * (x513 + x521)
    x872 = x468 * x682 + x653
    x873 = x468 * x684 + x660
    x874 = x873 * x96
    x875 = x4 * (x468 * x676 - x50 * (x468 * x679 + x649) + x646)
    x876 = x44 * (x468 * x688 + x666)
    x877 = 241.359114645008 * x111
    x878 = x4 * x468 * (x693 - x705)
    x879 = x468 * x709
    x880 = x629 * (x558 + x560)
    x881 = x44 * x880
    x882 = x629 * (x565 + x566)
    x883 = x4 * (-x50 * (x559 * x629 + x560 * x635) + x554 * x629 + x669 * x89)
    x884 = 2.0 * x689
    x885 = x468 * x728 + x884
    x886 = x44 * x885
    x887 = x468 * x733 + x685
    x888 = x4 * (x468 * x723 - x50 * (x468 * x727 + x681) + x678)
    x889 = x468 * x742 + x709
    x890 = x44 * x889
    x891 = x468 * x746 + x703
    x892 = x4 * (x468 * x738 - x50 * (x468 * x741 + x699) + x695)
    x893 = x468 * x759
    x894 = x4 * x468 * (x751 - x762)
    x895 = x629 * (x609 + x610)
    x896 = x4 * x629 * (-x50 * (x605 + x606) + x602 + x603)
    x897 = x468 * x771 + 3.0 * x731
    x898 = x4 * (x468 * x768 - x50 * (x468 * x769 + 3.0 * x725) + 3.0 * x722)
    x899 = x468 * x777 + 2.0 * x745
    x900 = x4 * (x468 * x774 - x50 * (x468 * x775 + 2.0 * x740) + 2.0 * x737)
    x901 = x468 * x783 + x759
    x902 = x4 * (x468 * x780 - x50 * (x468 * x781 + x753) + x750)
    x903 = x4 * x468 * (x786 - x788)
    x904 = x629**2
    x905 = x69 * x904 + x75
    x906 = x3 * x905
    x907 = x28 * x904 + x82
    x908 = x44 * x907
    x909 = 4.0 * x908
    x910 = x108 - x50 * (x109 + x42 * x904) + x57 * x904
    x911 = x74 * x910
    x912 = x44 * (x34 * x904 + x95)
    x913 = x141 * x904 + x146
    x914 = x119 * x904 + x152
    x915 = x44 * x914
    x916 = x130 * x904 + x171 - x50 * (x120 * x904 + x172)
    x917 = x74 * x916
    x918 = x44 * (x166 + x492 * x904)
    x919 = 2.0 * x918
    x920 = x202 + x629 * x650 + x858
    x921 = x3 * x920
    x922 = x209 + x629 * x652 + x860
    x923 = x44 * x922
    x924 = 3.0 * x923
    x925 = x225 - x50 * (x226 + x629 * x647 + x863) + x629 * x644 + x862
    x926 = x74 * x925
    x927 = x44 * (x221 + x629 * x659 + x865)
    x928 = 2.0 * x927
    x929 = x248 * x904 + x253
    x930 = x250 * x904 + x257
    x931 = x44 * x930
    x932 = x241 * x904 + x270 - x50 * (x234 * x904 + x271)
    x933 = x74 * x932
    x934 = x44 * (x231 * x904 + x265)
    x935 = x143 * x629 + x287 + x629 * x682
    x936 = x292 + x629 * x683 + x629 * x684
    x937 = x936 * x96
    x938 = x302 + x44 * x636 - x50 * (x126 * x629 + x303 + x629 * x679) + x629 * x676
    x939 = x74 * x938
    x940 = x44 * (x297 + x629 * x688 + x643)
    x941 = x331 + x629 * x700 + 2.0 * x653
    x942 = x3 * x941
    x943 = x334 + x629 * x702 + x661
    x944 = x44 * x943
    x945 = 2.0 * x944
    x946 = x349 - x50 * (x350 + x629 * x696 + 2.0 * x649) + x629 * x693 + 2.0 * x646
    x947 = x74 * x946
    x948 = x44 * (x342 + x629 * x708 + 2.0 * x666)
    x949 = x354 * x904
    x950 = x356 + x949
    x951 = x44 * x950
    x952 = x366 * x5 * x904 + x368
    x953 = x362 * x904 + x371 - x50 * (x372 + x5 * x949)
    x954 = x74 * x953
    x955 = x377 + x629 * x728 + x667
    x956 = x44 * x955
    x957 = x386 + x629 * x732 + x629 * x733
    x958 = x388 + x44 * x669 - x50 * (x389 + x629 * x727 + x668) + x629 * x723
    x959 = x74 * x958
    x960 = x392 + x629 * x742 + x884
    x961 = x44 * x960
    x962 = x400 + x629 * x746 + x685
    x963 = x401 - x50 * (x403 + x629 * x741 + x681) + x629 * x738 + x678
    x964 = x74 * x963
    x965 = x409 + x629 * x756 + 3.0 * x709
    x966 = x44 * x965
    x967 = x421 + x629 * x760 + 3.0 * x703
    x968 = x3 * x967
    x969 = x424 - x50 * (x425 + x629 * x755 + 3.0 * x699) + x629 * x751 + 3.0 * x695
    x970 = x74 * x969
    x971 = x430 * x904 + x432
    x972 = x426 * x904 + x434 - x50 * (x427 * x904 + x435)
    x973 = x74 * x972
    x974 = x440 + x629 * x771 + x718
    x975 = x442 - x50 * (x443 + x629 * x769 + x715) + x629 * x768 + x714
    x976 = x74 * x975
    x977 = x448 + x629 * x777 + 2.0 * x731
    x978 = x450 - x50 * (x451 + x629 * x775 + 2.0 * x725) + x629 * x774 + 2.0 * x722
    x979 = x74 * x978
    x980 = x453 + x629 * x783 + 3.0 * x745
    x981 = x455 - x50 * (x457 + x629 * x781 + 3.0 * x740) + x629 * x780 + 3.0 * x737
    x982 = x74 * x981
    x983 = x464 + x629 * x789 + 4.0 * x759
    x984 = x466 - x50 * (x467 + x629 * x787 + 4.0 * x753) + x629 * x786 + 4.0 * x750
    x985 = x74 * x984
    x986 = x468 * x795 + x473
    x987 = x44 * (x468 * x797 + x477)
    x988 = x468 * x802 + x490 + x798
    x989 = x468 * x804 + x497 + x800
    x990 = x468 * x808 + x502
    x991 = x44 * (x468 * x809 + x505)
    x992 = x468 * x814 + x520 + x804 * x96
    x993 = x468 * x815 + x525 + x807
    x994 = x468 * x819 + x538 + x810
    x995 = x96 * (x468 * x820 + x542 + x813)
    x996 = x468 * x825 + x547
    x997 = x44 * (x468 * x826 + x551)
    x998 = x44 * (x468 * x830 + x563 + 3.0 * x817)
    x999 = x468 * x832 + x568 + x815 * x89
    x1000 = x44 * (x468 * x834 + x580 + 2.0 * x823)
    x1001 = x468 * x836 + x583 + x821
    x1002 = x44 * (x468 * x838 + x592 + x829)
    x1003 = x468 * x840 + x596 + x827
    x1004 = x44 * (x468 * x843 + x599)
    x1005 = x468 * x845 + x601
    x1006 = x468 * x847 + x612 + 4.0 * x831
    x1007 = x173 * x3
    x1008 = x468 * x849 + x617 + 3.0 * x835
    x1009 = x3 * x500
    x1010 = x468 * x851 + x622 + 2.0 * x839
    x1011 = x3 * x526
    x1012 = x468 * x853 + x627 + x844
    x1013 = x468 * x855 + x628
    x1014 = x629 * x794 + x631
    x1015 = x44 * (x629 * x796 + x634)
    x1016 = x468 * x859 + x629 * x801 + x638
    x1017 = x468 * x861 + x629 * x803 + x642
    x1018 = x650 * x793 + x657
    x1019 = x44 * (x652 * x793 + x664)
    x1020 = x468 * x868 + x671 + x861 * x96
    x1021 = x468 * x869 + x675 + x866
    x1022 = x468 * x653 + x468 * x872 + x687
    x1023 = x96 * (x468 * x660 + x468 * x873 + x692)
    x1024 = x700 * x793 + x707
    x1025 = x44 * (x702 * x793 + x713)
    x1026 = x44 * (x468 * x880 + x717 + 3.0 * x871)
    x1027 = x468 * x882 + x720 + x869 * x89
    x1028 = x44 * (x468 * x885 + x730 + 2.0 * x876)
    x1029 = x468 * x887 + x735 + x874
    x1030 = x44 * (x468 * x889 + x744 + x879)
    x1031 = x468 * x703 + x468 * x891 + x748
    x1032 = x44 * (x756 * x793 + x758)
    x1033 = x760 * x793 + x764
    x1034 = x468 * x895 + x767 + 4.0 * x881
    x1035 = x3 * x304
    x1036 = x468 * x897 + x773 + 3.0 * x886
    x1037 = x3 * x543
    x1038 = x468 * x899 + x779 + 2.0 * x890
    x1039 = x3 * x824
    x1040 = x468 * x901 + x785 + x893
    x1041 = x789 * x793 + x792
    x1042 = x468 * x913 + x908
    x1043 = x468 * x914 + x912
    x1044 = x468 * x929 + 2.0 * x915
    x1045 = x468 * x930 + x919
    x1046 = x468 * x935 + x923
    x1047 = x96 * (x468 * x936 + x927)
    x1048 = x44 * (x468 * x950 + 3.0 * x934)
    x1049 = x468 * x952 + 3.0 * x931
    x1050 = 2.0 * x940
    x1051 = x44 * (x1050 + x468 * x955)
    x1052 = x468 * x957 + x937
    x1053 = x44 * (x468 * x960 + x948)
    x1054 = x468 * x962 + x944
    x1055 = x468 * x966
    x1056 = x468 * x971 + 4.0 * x951
    x1057 = x468 * x974 + 3.0 * x956
    x1058 = x468 * x977 + 2.0 * x961
    x1059 = x468 * x980 + x966
    x1060 = x629 * x905 + x630
    x1061 = x44 * (x629 * x907 + x633)
    x1062 = x629 * x913 + x637
    x1063 = x44 * (x629 * x914 + x641)
    x1064 = x629 * x920 + x656 + x908
    x1065 = x44 * (x629 * x922 + x663 + x912)
    x1066 = x629 * x929 + x670
    x1067 = x44 * (x629 * x930 + x674)
    x1068 = x629 * x935 + x686 + x915
    x1069 = x96 * (x629 * x936 + x691 + x918)
    x1070 = x629 * x941 + x706 + 2.0 * x923
    x1071 = x44 * (x629 * x943 + x712 + x928)
    x1072 = x44 * (x629 * x950 + x716)
    x1073 = x629 * x952 + x719
    x1074 = x44 * (x629 * x955 + x729 + x934)
    x1075 = x629 * x957 + x734 + x931
    x1076 = x44 * (x1050 + x629 * x960 + x743)
    x1077 = x629 * x962 + x747 + x937
    x1078 = x44 * (x629 * x965 + x757 + 3.0 * x948)
    x1079 = x629 * x967 + x763 + 3.0 * x944
    x1080 = x629 * x971 + x766
    x1081 = x629 * x974 + x772 + x951
    x1082 = x629 * x977 + x778 + 2.0 * x956
    x1083 = x629 * x980 + x784 + 3.0 * x961
    x1084 = x629 * x983 + x791 + 4.0 * x966
    x1085 = x173 * x468
    x1086 = x468 * x500

    # 225 item(s)
    result[0, 0] = numpy.sum(
        x112
        * (
            x110 * (x108 + x3 * x63 - x50 * (x109 + x3 * x49 + x80 * x92) + x80 * x94)
            + x3
            * (
                x3 * (x3 * (x70 + x72) + x75 + x79 * x80)
                - x4 * (x49 * x50 - x63)
                + x80 * x90
            )
            + x80
            * (
                x3 * x90
                - x4 * (x50 * x92 - x94)
                + x89 * (x3 * x88 + x95 + x96 * (x38 * x99 + x98))
            )
        )
    )
    result[0, 1] = numpy.sum(
        x173
        * (
            x110 * (x137 * x3 + x164 * x89 + x171 - x50 * (x128 * x3 + x161 * x89 + x172))
            + x3
            * (
                x158 * x89
                + x3 * (x146 + x150 * x89 + x3 * (x142 + x144))
                - x4 * (x128 * x50 - x137)
            )
            + x89
            * (
                x158 * x3
                - x4 * (x161 * x50 - x164)
                + x96 * (x154 * x3 + x157 * x3 + x166)
            )
        )
    )
    result[0, 2] = numpy.sum(
        x173
        * (
            x110 * (x194 * x3 + x219 * x89 + x225 - x50 * (x186 * x3 + x217 * x89 + x226))
            + x3
            * (
                x215 * x89
                + x3 * (x202 + x207 * x89 + x3 * (x198 + x200))
                - x4 * (x186 * x50 - x194)
            )
            + x89
            * (
                x215 * x3
                - x4 * (x217 * x50 - x219)
                + x96 * (x211 * x3 + x214 * x3 + x221)
            )
        )
    )
    result[0, 3] = numpy.sum(
        x272
        * (
            x110 * (x244 * x3 + x263 * x96 + x270 - x50 * (x238 * x3 + x261 * x96 + x271))
            + x3
            * (
                x259 * x96
                + x3 * (x253 + x255 * x96 + x3 * (x249 + x251))
                - x4 * (x238 * x50 - x244)
            )
            + x96 * (x259 * x3 - x4 * (x261 * x50 - x263) + x44 * (x231 * x264 + x265))
        )
    )
    result[0, 4] = numpy.sum(
        x304
        * (
            x110 * (x281 * x3 + x296 * x96 + x302 - x50 * (x276 * x3 + x294 * x96 + x303))
            + x3
            * (
                x293 * x96
                + x3 * (x287 + x291 * x96 + x3 * (x283 * x3 + x286))
                - x4 * (x276 * x50 - x281)
            )
            + x96 * (x293 * x3 - x4 * (x294 * x50 - x296) + x44 * (x264 * x298 + x297))
        )
    )
    result[0, 5] = numpy.sum(
        x272
        * (
            x110 * (x3 * x323 + x341 * x96 + x349 - x50 * (x3 * x316 + x339 * x96 + x350))
            + x3
            * (
                x3 * (x3 * (x328 + x329 * x96) + x331 + x333 * x96)
                + x336 * x96
                - x4 * (x316 * x50 - x323)
            )
            + x96 * (x3 * x336 - x4 * (x339 * x50 - x341) + x44 * (x264 * x308 + x342))
        )
    )
    result[0, 6] = numpy.sum(
        x173
        * (
            x110
            * (x3 * x351 * x44 + x3 * x364 + x371 - x50 * (x3 * x358 + x3 * x360 + x372))
            + x3 * x44 * (x357 + x4 * (x351 - x353))
            + x3
            * (
                x3 * (x3 * (x365 + x367) + x368 + x369)
                + x357 * x44
                - x4 * (x360 * x50 - x364)
            )
        )
    )
    result[0, 7] = numpy.sum(
        x304
        * (
            x110
            * (x3 * x373 * x44 + x3 * x383 + x388 - x50 * (x3 * x379 + x3 * x380 + x389))
            + x3 * x44 * (x378 + x4 * (x373 - x375))
            + x3
            * (
                x3 * (x3 * x384 + x3 * (x384 + x385 * x99) + x386)
                + x378 * x44
                - x4 * (x380 * x50 - x383)
            )
        )
    )
    result[0, 8] = numpy.sum(
        x304
        * (
            x110 * (x3 * x397 + x3 * x398 + x401 - x50 * (x114 * x335 + x3 * x396 + x403))
            + x3 * x44 * (x394 + x4 * (x114 * x309 - x391))
            + x3
            * (
                x3 * (x3 * x399 + x3 * (x155 * x324 + x399) + x400)
                + x394 * x44
                - x4 * (x396 * x50 - x398)
            )
        )
    )
    result[0, 9] = numpy.sum(
        x173
        * (
            x110
            * (x3 * x404 * x44 + x3 * x417 + x424 - x50 * (x3 * x411 + x3 * x413 + x425))
            + x3 * x44 * (x4 * (x404 - x406) + x410)
            + x3
            * (
                x3 * (x3 * (x418 + x420) + x421 + x422)
                - x4 * (x413 * x50 - x417)
                + x410 * x44
            )
        )
    )
    result[0, 10] = numpy.sum(
        x112
        * (
            x110 * (x264 * x426 + x434 - x50 * (x264 * x427 + x435))
            + x3 * (x3 * (x431 + x432) + x4 * (x3 * x426 - x429))
        )
    )
    result[0, 11] = numpy.sum(
        x173
        * (
            x110 * (x264 * x436 + x442 - x50 * (x264 * x437 + x443))
            + x3**2 * (x264 * x439 + x4 * (x436 - x438) + x440)
        )
    )
    result[0, 12] = numpy.sum(
        x272
        * (
            x110 * (x264 * x444 + x450 - x50 * (x264 * x445 + x451))
            + x3**2 * (x264 * x447 + x4 * (x444 - x446) + x448)
        )
    )
    result[0, 13] = numpy.sum(
        x173
        * (
            x110 * (x264 * x456 + x455 - x50 * (x114 * x408 + x457))
            + x3**2 * (x264 * x454 + x4 * (x114 * x405 - x452) + x453)
        )
    )
    result[0, 14] = numpy.sum(
        x112
        * (
            x110 * (x264 * x458 + x466 - x50 * (x264 * x459 + x467))
            + x3 * (x3 * (x463 + x464) + x4 * (x3 * x458 - x461))
        )
    )
    result[1, 0] = numpy.sum(
        0.5
        * x173
        * (
            x3 * (2.0 * x3 * (x470 + x471) + x473 + 2.0 * x476 * x80)
            + 2.0 * x4 * (-x468 * x50 * (x43 + x48) + x468 * x58 + x469 * x80)
            + x80 * (2.0 * x3 * x476 + x477 + 2.0 * x89 * (x468 * x83 + x478 * x479))
        )
    )
    result[1, 1] = numpy.sum(
        0.5
        * x500
        * (
            x3 * (2.0 * x3 * (x3 * x487 + x489 * x89) + x490 + 2.0 * x496 * x89)
            + 2.0 * x4 * (x3 * x481 + x482 * x89 - x50 * (x3 * x484 + x485 * x89))
            + x89 * (2.0 * x3 * x496 + x497 + 2.0 * x96 * (x3 * x494 + x499))
        )
    )
    result[1, 2] = numpy.sum(
        0.5
        * x500
        * (
            x3 * (2.0 * x3 * x468 * (x198 + x200) + x502 + 2.0 * x504 * x89)
            + 2.0 * x4 * (x189 * x468 - x468 * x50 * (x180 + x185) + x501 * x89)
            + x89 * (2.0 * x3 * x504 + x505 + 2.0 * x96 * (x213 * x468 + x506))
        )
    )
    result[1, 3] = numpy.sum(
        0.5
        * x526
        * (
            x3 * (2.0 * x3 * (x3 * x517 + x519 * x96) + x520 + 2.0 * x524 * x96)
            + 2.0 * x4 * (x3 * x508 - x50 * (x3 * x512 + x514 * x96) + x509 * x96)
            + x96 * (2.0 * x3 * x523 + 2.0 * x3 * x524 + x525)
        )
    )
    result[1, 4] = numpy.sum(
        0.5
        * x543
        * (
            x3 * (2.0 * x3 * (x3 * x534 + x537) + x538 + 2.0 * x541 * x96)
            + 2.0 * x4 * (x3 * x527 - x50 * (x3 * x530 + x533) + x529)
            + x96 * (2.0 * x3 * x540 + 2.0 * x3 * x541 + x542)
        )
    )
    result[1, 5] = numpy.sum(
        0.5
        * x526
        * (
            x3 * (2.0 * x3 * (x328 * x468 + x546 * x96) + x547 + 2.0 * x549 * x96)
            + 2.0 * x4 * (x321 * x468 - x50 * (x313 * x468 + 2.0 * x544) + x545 * x96)
            + x96 * (2.0 * x3 * x549 + 2.0 * x335 * x468 + x551)
        )
    )
    result[1, 6] = numpy.sum(
        0.5
        * x500
        * (
            x3 * (2.0 * x3 * x564 + 2.0 * x3 * (x3 * x567 + x564) + x568)
            + 2.0 * x4 * (x3 * x555 - x50 * (x3 * x561 + x557) + x553)
            + x44 * (2.0 * x264 * x562 + x563)
        )
    )
    result[1, 7] = numpy.sum(
        0.5
        * x543
        * (
            x3 * (2.0 * x3 * x581 + 2.0 * x3 * (x3 * x582 + x581) + x583)
            + 2.0 * x4 * (x3 * x572 - x50 * (x3 * x577 + x575) + x571)
            + x44 * (2.0 * x264 * x579 + x580)
        )
    )
    result[1, 8] = numpy.sum(
        0.5
        * x543
        * (
            x3 * (2.0 * x3 * x593 + 2.0 * x3 * (x3 * x595 + x593) + x596)
            + 2.0 * x4 * (x3 * x587 - x50 * (x3 * x590 + x589) + x586)
            + x44 * (2.0 * x264 * x591 + x592)
        )
    )
    result[1, 9] = numpy.sum(
        0.5
        * x500
        * (
            x3 * (2.0 * x3 * (x420 * x468 + x600) + 2.0 * x422 * x468 + x601)
            + 2.0 * x4 * (x416 * x468 - x50 * (x412 * x468 + x598) + x597)
            + x44 * (2.0 * x408 * x468 + x599)
        )
    )
    result[1, 10] = numpy.sum(
        0.5 * x173 * x3 * (2.0 * x264 * x611 + 2.0 * x4 * (x604 - x608) + x612)
    )
    result[1, 11] = numpy.sum(
        0.5 * x3 * x500 * (2.0 * x264 * x616 + 2.0 * x4 * (x613 - x615) + x617)
    )
    result[1, 12] = numpy.sum(
        0.5 * x3 * x526 * (2.0 * x264 * x621 + 2.0 * x4 * (x618 - x620) + x622)
    )
    result[1, 13] = numpy.sum(
        0.5 * x3 * x500 * (2.0 * x264 * x626 + 2.0 * x4 * (x623 - x625) + x627)
    )
    result[1, 14] = numpy.sum(
        0.5
        * x173
        * (x3 * (2.0 * x463 * x468 + x628) + 2.0 * x4 * x468 * (x3 * x458 - x461))
    )
    result[2, 0] = numpy.sum(
        x173
        * (
            x3 * (x3 * x629 * (x70 + x72) + x631 + x632 * x80)
            + x4 * x629 * (-x50 * (x43 + x48) + x58 + x62)
            + x80 * (x3 * x632 + x634 + x89 * (x479 * x635 + x629 * x83))
        )
    )
    result[2, 1] = numpy.sum(
        x500
        * (
            x3 * (x3 * x629 * (x142 + x144) + x638 + x640 * x89)
            + x4 * (x131 * x629 - x50 * x629 * (x121 + x127) + x636 * x89)
            + x89 * (x3 * x640 + x642 + x96 * (x156 * x629 + x643))
        )
    )
    result[2, 2] = numpy.sum(
        x500
        * (
            x3 * (x3 * (x651 + x654) + x657 + x662 * x89)
            + x4 * (x3 * x644 - x50 * (x3 * x647 + 3.0 * x649) + 3.0 * x646)
            + x89 * (x3 * x662 + x664 + x96 * (x3 * x659 + x666))
        )
    )
    result[2, 3] = numpy.sum(
        x526
        * (
            x3 * (x3 * x629 * (x249 + x251) + x671 + x672 * x96)
            + x4 * (x242 * x629 - x50 * (x235 * x629 + 2.0 * x668) + x669 * x96)
            + x96 * (x258 * x629 + x3 * x672 + x675)
        )
    )
    result[2, 4] = numpy.sum(
        x543
        * (
            x3 * (x3 * (x3 * x682 + x685) + x687 + x690 * x96)
            + x4 * (x3 * x676 - x50 * (x3 * x679 + x681) + x678)
            + x96 * (x3 * x689 + x3 * x690 + x692)
        )
    )
    result[2, 5] = numpy.sum(
        x526
        * (
            x3 * (x3 * (x701 + x704) + x707 + x711 * x96)
            + x4 * (x3 * x693 - x50 * (x3 * x696 + 2.0 * x699) + 2.0 * x695)
            + x96 * (x3 * x709 + x3 * x711 + x713)
        )
    )
    result[2, 6] = numpy.sum(
        x500
        * (
            x3 * (x3 * (x367 * x629 + x718) + x369 * x629 + x720)
            + x4 * (x363 * x629 - x50 * (x359 * x629 + x715) + x714)
            + x44 * (x355 * x629 + x717)
        )
    )
    result[2, 7] = numpy.sum(
        x543
        * (
            x3 * (x3 * x731 + x3 * (x3 * x733 + x731) + x735)
            + x4 * (x3 * x723 - x50 * (x3 * x727 + x725) + x722)
            + x44 * (x264 * x728 + x730)
        )
    )
    result[2, 8] = numpy.sum(
        x543
        * (
            x3 * (x3 * x745 + x3 * (x3 * x746 + x745) + x748)
            + x4 * (x3 * x738 - x50 * (x3 * x741 + x740) + x737)
            + x44 * (x264 * x742 + x744)
        )
    )
    result[2, 9] = numpy.sum(
        x500
        * (
            x3 * (x3 * (x759 + x761) + x764 + x765)
            + x4 * (x3 * x751 - x50 * (x3 * x755 + x753) + x750)
            + x44 * (x264 * x756 + x758)
        )
    )
    result[2, 10] = numpy.sum(
        x173 * (x3 * (x431 * x629 + x767) + x4 * x629 * (x3 * x426 - x429))
    )
    result[2, 11] = numpy.sum(x3 * x500 * (x264 * x771 + x4 * (x768 - x770) + x773))
    result[2, 12] = numpy.sum(x3 * x526 * (x264 * x777 + x4 * (x774 - x776) + x779))
    result[2, 13] = numpy.sum(x3 * x500 * (x264 * x783 + x4 * (x780 - x782) + x785))
    result[2, 14] = numpy.sum(x173 * x3 * (x4 * (x786 - x788) + x790 + x792))
    result[3, 0] = numpy.sum(
        x272
        * (x3 * (x3 * x795 + 4.0 * x798) + x74 * x799 + x80 * (x3 * x797 + 3.0 * x800))
    )
    result[3, 1] = numpy.sum(
        x526 * (x3 * (x3 * x802 + x804 * x89) + x74 * x805 + x89 * (x3 * x804 + x807))
    )
    result[3, 2] = numpy.sum(
        x526
        * (x3 * (x3 * x808 + 3.0 * x810) + x74 * x811 + x89 * (x3 * x809 + 2.0 * x813))
    )
    result[3, 3] = numpy.sum(
        x818 * (x3 * (x3 * x814 + x815 * x96) + x74 * x816 + x96 * (x3 * x815 + x817))
    )
    result[3, 4] = numpy.sum(
        x824 * (x3 * (x3 * x819 + x821) + x74 * x822 + x96 * (x3 * x820 + x823))
    )
    result[3, 5] = numpy.sum(
        x818 * (x3 * (x3 * x825 + 2.0 * x827) + x74 * x828 + x96 * (x3 * x826 + x829))
    )
    result[3, 6] = numpy.sum(x526 * (x3 * x831 + x3 * (x3 * x832 + x831) + x74 * x833))
    result[3, 7] = numpy.sum(x824 * (x3 * x835 + x3 * (x3 * x836 + x835) + x74 * x837))
    result[3, 8] = numpy.sum(x824 * (x3 * x839 + x3 * (x3 * x840 + x839) + x74 * x841))
    result[3, 9] = numpy.sum(x526 * (x3 * x844 + x3 * (x3 * x845 + x844) + x74 * x846))
    result[3, 10] = numpy.sum(x272 * (x264 * x847 + x74 * x848))
    result[3, 11] = numpy.sum(x526 * (x264 * x849 + x74 * x850))
    result[3, 12] = numpy.sum(x818 * (x264 * x851 + x74 * x852))
    result[3, 13] = numpy.sum(x526 * (x264 * x853 + x74 * x854))
    result[3, 14] = numpy.sum(x272 * (x264 * x855 + x74 * x856))
    result[4, 0] = numpy.sum(
        0.5
        * x304
        * (2.0 * x3 * x629 * (x470 + x471) + 2.0 * x629 * x80 * (x474 + x475) + x857)
    )
    result[4, 1] = numpy.sum(
        0.5
        * x543
        * (2.0 * x3 * (x3 * x859 + x861 * x89) + x864 + 2.0 * x89 * (x3 * x861 + x866))
    )
    result[4, 2] = numpy.sum(
        0.5
        * x543
        * (2.0 * x3 * x468 * (x651 + x654) + 2.0 * x468 * x89 * (x658 + x661) + x867)
    )
    result[4, 3] = numpy.sum(
        0.5
        * x824
        * (2.0 * x3 * (x3 * x868 + x869 * x96) + x870 + 2.0 * x96 * (x3 * x869 + x871))
    )
    result[4, 4] = numpy.sum(
        0.5
        * x877
        * (2.0 * x3 * (x3 * x872 + x874) + x875 + 2.0 * x96 * (x3 * x873 + x876))
    )
    result[4, 5] = numpy.sum(
        0.5
        * x824
        * (2.0 * x3 * x468 * (x701 + x704) + x878 + 2.0 * x96 * (x468 * x710 + x879))
    )
    result[4, 6] = numpy.sum(0.5 * x543 * (2.0 * x3**2 * x882 + 4.0 * x3 * x881 + x883))
    result[4, 7] = numpy.sum(0.5 * x877 * (2.0 * x3**2 * x887 + 4.0 * x3 * x886 + x888))
    result[4, 8] = numpy.sum(0.5 * x877 * (2.0 * x3**2 * x891 + 4.0 * x3 * x890 + x892))
    result[4, 9] = numpy.sum(
        0.5 * x543 * (2.0 * x3 * (x468 * x761 + x893) + 2.0 * x468 * x765 + x894)
    )
    result[4, 10] = numpy.sum(0.5 * x304 * (2.0 * x264 * x895 + x896))
    result[4, 11] = numpy.sum(0.5 * x543 * (2.0 * x264 * x897 + x898))
    result[4, 12] = numpy.sum(0.5 * x824 * (2.0 * x264 * x899 + x900))
    result[4, 13] = numpy.sum(0.5 * x543 * (2.0 * x264 * x901 + x902))
    result[4, 14] = numpy.sum(0.5 * x304 * (2.0 * x468 * x790 + x903))
    result[5, 0] = numpy.sum(
        x272 * (x3 * (x906 + x909) + x80 * (x3 * x907 + 3.0 * x912) + x911)
    )
    result[5, 1] = numpy.sum(
        x526 * (x3 * (x3 * x913 + 3.0 * x915) + x89 * (x3 * x914 + x919) + x917)
    )
    result[5, 2] = numpy.sum(
        x526 * (x3 * (x921 + x924) + x89 * (x3 * x922 + x928) + x926)
    )
    result[5, 3] = numpy.sum(
        x818 * (x3 * (x3 * x929 + 2.0 * x931) + x933 + x96 * (x3 * x930 + x934))
    )
    result[5, 4] = numpy.sum(
        x824 * (x3 * (x3 * x935 + x937) + x939 + x96 * (x3 * x936 + x940))
    )
    result[5, 5] = numpy.sum(
        x818 * (x3 * (x942 + x945) + x947 + x96 * (x3 * x943 + x948))
    )
    result[5, 6] = numpy.sum(x526 * (x3 * x951 + x3 * (x3 * x952 + x951) + x954))
    result[5, 7] = numpy.sum(x824 * (x3 * x956 + x3 * (x3 * x957 + x956) + x959))
    result[5, 8] = numpy.sum(x824 * (x3 * x961 + x3 * (x3 * x962 + x961) + x964))
    result[5, 9] = numpy.sum(x526 * (x3 * x966 + x3 * (x966 + x968) + x970))
    result[5, 10] = numpy.sum(x272 * (x264 * x971 + x973))
    result[5, 11] = numpy.sum(x526 * (x264 * x974 + x976))
    result[5, 12] = numpy.sum(x818 * (x264 * x977 + x979))
    result[5, 13] = numpy.sum(x526 * (x264 * x980 + x982))
    result[5, 14] = numpy.sum(x272 * (x264 * x983 + x985))
    result[6, 0] = numpy.sum(x173 * (x3 * x986 + 4.0 * x987))
    result[6, 1] = numpy.sum(x500 * (x3 * x988 + x89 * x989))
    result[6, 2] = numpy.sum(x500 * (x3 * x990 + 3.0 * x991))
    result[6, 3] = numpy.sum(x526 * (x3 * x992 + x96 * x993))
    result[6, 4] = numpy.sum(x543 * (x3 * x994 + x995))
    result[6, 5] = numpy.sum(x526 * (x3 * x996 + 2.0 * x997))
    result[6, 6] = numpy.sum(x500 * (x3 * x999 + x998))
    result[6, 7] = numpy.sum(x543 * (x1000 + x1001 * x3))
    result[6, 8] = numpy.sum(x543 * (x1002 + x1003 * x3))
    result[6, 9] = numpy.sum(x500 * (x1004 + x1005 * x3))
    result[6, 10] = numpy.sum(x1006 * x1007)
    result[6, 11] = numpy.sum(x1008 * x1009)
    result[6, 12] = numpy.sum(x1010 * x1011)
    result[6, 13] = numpy.sum(x1009 * x1012)
    result[6, 14] = numpy.sum(x1007 * x1013)
    result[7, 0] = numpy.sum(x304 * (x1014 * x3 + 4.0 * x1015))
    result[7, 1] = numpy.sum(x543 * (x1016 * x3 + x1017 * x89))
    result[7, 2] = numpy.sum(x543 * (x1018 * x3 + 3.0 * x1019))
    result[7, 3] = numpy.sum(x824 * (x1020 * x3 + x1021 * x96))
    result[7, 4] = numpy.sum(x877 * (x1022 * x3 + x1023))
    result[7, 5] = numpy.sum(x824 * (x1024 * x3 + 2.0 * x1025))
    result[7, 6] = numpy.sum(x543 * (x1026 + x1027 * x3))
    result[7, 7] = numpy.sum(x877 * (x1028 + x1029 * x3))
    result[7, 8] = numpy.sum(x877 * (x1030 + x1031 * x3))
    result[7, 9] = numpy.sum(x543 * (x1032 + x1033 * x3))
    result[7, 10] = numpy.sum(x1034 * x1035)
    result[7, 11] = numpy.sum(x1036 * x1037)
    result[7, 12] = numpy.sum(x1038 * x1039)
    result[7, 13] = numpy.sum(x1037 * x1040)
    result[7, 14] = numpy.sum(x1035 * x1041)
    result[8, 0] = numpy.sum(x304 * x468 * (x906 + x909))
    result[8, 1] = numpy.sum(x543 * (x1042 * x3 + x1043 * x89))
    result[8, 2] = numpy.sum(x468 * x543 * (x921 + x924))
    result[8, 3] = numpy.sum(x824 * (x1044 * x3 + x1045 * x96))
    result[8, 4] = numpy.sum(x877 * (x1046 * x3 + x1047))
    result[8, 5] = numpy.sum(x468 * x824 * (x942 + x945))
    result[8, 6] = numpy.sum(x543 * (x1048 + x1049 * x3))
    result[8, 7] = numpy.sum(x877 * (x1051 + x1052 * x3))
    result[8, 8] = numpy.sum(x877 * (x1053 + x1054 * x3))
    result[8, 9] = numpy.sum(x543 * (x1055 + x468 * x968))
    result[8, 10] = numpy.sum(x1035 * x1056)
    result[8, 11] = numpy.sum(x1037 * x1057)
    result[8, 12] = numpy.sum(x1039 * x1058)
    result[8, 13] = numpy.sum(x1037 * x1059)
    result[8, 14] = numpy.sum(x1035 * x468 * x983)
    result[9, 0] = numpy.sum(x173 * (x1060 * x3 + 4.0 * x1061))
    result[9, 1] = numpy.sum(x500 * (x1062 * x3 + 3.0 * x1063))
    result[9, 2] = numpy.sum(x500 * (x1064 * x3 + 3.0 * x1065))
    result[9, 3] = numpy.sum(x526 * (x1066 * x3 + 2.0 * x1067))
    result[9, 4] = numpy.sum(x543 * (x1068 * x3 + x1069))
    result[9, 5] = numpy.sum(x526 * (x1070 * x3 + 2.0 * x1071))
    result[9, 6] = numpy.sum(x500 * (x1072 + x1073 * x3))
    result[9, 7] = numpy.sum(x543 * (x1074 + x1075 * x3))
    result[9, 8] = numpy.sum(x543 * (x1076 + x1077 * x3))
    result[9, 9] = numpy.sum(x500 * (x1078 + x1079 * x3))
    result[9, 10] = numpy.sum(x1007 * x1080)
    result[9, 11] = numpy.sum(x1009 * x1081)
    result[9, 12] = numpy.sum(x1011 * x1082)
    result[9, 13] = numpy.sum(x1009 * x1083)
    result[9, 14] = numpy.sum(x1007 * x1084)
    result[10, 0] = numpy.sum(x112 * (x110 * x799 + x468 * x986))
    result[10, 1] = numpy.sum(x173 * (x110 * x805 + x468 * x988 + x987))
    result[10, 2] = numpy.sum(x173 * (x110 * x811 + x468 * x990))
    result[10, 3] = numpy.sum(x272 * (x110 * x816 + x468 * x992 + x96 * x989))
    result[10, 4] = numpy.sum(x304 * (x110 * x822 + x468 * x994 + x991))
    result[10, 5] = numpy.sum(x272 * (x110 * x828 + x468 * x996))
    result[10, 6] = numpy.sum(x173 * (x110 * x833 + x468 * x999 + x89 * x993))
    result[10, 7] = numpy.sum(x304 * (x1001 * x468 + x110 * x837 + x995))
    result[10, 8] = numpy.sum(x304 * (x1003 * x468 + x110 * x841 + x997))
    result[10, 9] = numpy.sum(x173 * (x1005 * x468 + x110 * x846))
    result[10, 10] = numpy.sum(x112 * (x1006 * x468 + x110 * x848 + 4.0 * x998))
    result[10, 11] = numpy.sum(x173 * (3.0 * x1000 + x1008 * x468 + x110 * x850))
    result[10, 12] = numpy.sum(x272 * (2.0 * x1002 + x1010 * x468 + x110 * x852))
    result[10, 13] = numpy.sum(x173 * (x1004 + x1012 * x468 + x110 * x854))
    result[10, 14] = numpy.sum(x112 * (x1013 * x468 + x110 * x856))
    result[11, 0] = numpy.sum(x173 * (x1014 * x468 + x857))
    result[11, 1] = numpy.sum(x500 * (x1015 + x1016 * x468 + x864))
    result[11, 2] = numpy.sum(x500 * (x1018 * x468 + x867))
    result[11, 3] = numpy.sum(x526 * (x1017 * x96 + x1020 * x468 + x870))
    result[11, 4] = numpy.sum(x543 * (x1019 + x1022 * x468 + x875))
    result[11, 5] = numpy.sum(x526 * (x1024 * x468 + x878))
    result[11, 6] = numpy.sum(x500 * (x1021 * x89 + x1027 * x468 + x883))
    result[11, 7] = numpy.sum(x543 * (x1023 + x1029 * x468 + x888))
    result[11, 8] = numpy.sum(x543 * (x1025 + x1031 * x468 + x892))
    result[11, 9] = numpy.sum(x500 * (x1033 * x468 + x894))
    result[11, 10] = numpy.sum(x173 * (4.0 * x1026 + x1034 * x468 + x896))
    result[11, 11] = numpy.sum(x500 * (3.0 * x1028 + x1036 * x468 + x898))
    result[11, 12] = numpy.sum(x526 * (2.0 * x1030 + x1038 * x468 + x900))
    result[11, 13] = numpy.sum(x500 * (x1032 + x1040 * x468 + x902))
    result[11, 14] = numpy.sum(x173 * (x1041 * x468 + x903))
    result[12, 0] = numpy.sum(x272 * (x793 * x905 + x911))
    result[12, 1] = numpy.sum(x526 * (x1042 * x468 + x468 * x908 + x917))
    result[12, 2] = numpy.sum(x526 * (x793 * x920 + x926))
    result[12, 3] = numpy.sum(x818 * (x1043 * x96 + x1044 * x468 + x933))
    result[12, 4] = numpy.sum(x824 * (x1046 * x468 + x468 * x923 + x939))
    result[12, 5] = numpy.sum(x818 * (x793 * x941 + x947))
    result[12, 6] = numpy.sum(x526 * (x1045 * x89 + x1049 * x468 + x954))
    result[12, 7] = numpy.sum(x824 * (x1047 + x1052 * x468 + x959))
    result[12, 8] = numpy.sum(x824 * (x1054 * x468 + x468 * x944 + x964))
    result[12, 9] = numpy.sum(x526 * (x793 * x967 + x970))
    result[12, 10] = numpy.sum(x272 * (4.0 * x1048 + x1056 * x468 + x973))
    result[12, 11] = numpy.sum(x526 * (3.0 * x1051 + x1057 * x468 + x976))
    result[12, 12] = numpy.sum(x818 * (2.0 * x1053 + x1058 * x468 + x979))
    result[12, 13] = numpy.sum(x526 * (x1055 + x1059 * x468 + x982))
    result[12, 14] = numpy.sum(x272 * (x793 * x983 + x985))
    result[13, 0] = numpy.sum(x1060 * x1085)
    result[13, 1] = numpy.sum(x500 * (x1061 + x1062 * x468))
    result[13, 2] = numpy.sum(x1064 * x1086)
    result[13, 3] = numpy.sum(x526 * (2.0 * x1063 + x1066 * x468))
    result[13, 4] = numpy.sum(x543 * (x1065 + x1068 * x468))
    result[13, 5] = numpy.sum(x1070 * x468 * x526)
    result[13, 6] = numpy.sum(x500 * (3.0 * x1067 + x1073 * x468))
    result[13, 7] = numpy.sum(x543 * (x1069 + x1075 * x468))
    result[13, 8] = numpy.sum(x543 * (x1071 + x1077 * x468))
    result[13, 9] = numpy.sum(x1079 * x1086)
    result[13, 10] = numpy.sum(x173 * (4.0 * x1072 + x1080 * x468))
    result[13, 11] = numpy.sum(x500 * (3.0 * x1074 + x1081 * x468))
    result[13, 12] = numpy.sum(x526 * (2.0 * x1076 + x1082 * x468))
    result[13, 13] = numpy.sum(x500 * (x1078 + x1083 * x468))
    result[13, 14] = numpy.sum(x1084 * x1085)
    result[14, 0] = numpy.sum(x112 * (x1060 * x629 + x110 * x910))
    result[14, 1] = numpy.sum(x173 * (x1062 * x629 + x110 * x916))
    result[14, 2] = numpy.sum(x173 * (x1061 + x1064 * x629 + x110 * x925))
    result[14, 3] = numpy.sum(x272 * (x1066 * x629 + x110 * x932))
    result[14, 4] = numpy.sum(x304 * (x1063 + x1068 * x629 + x110 * x938))
    result[14, 5] = numpy.sum(x272 * (2.0 * x1065 + x1070 * x629 + x110 * x946))
    result[14, 6] = numpy.sum(x173 * (x1073 * x629 + x110 * x953))
    result[14, 7] = numpy.sum(x304 * (x1067 + x1075 * x629 + x110 * x958))
    result[14, 8] = numpy.sum(x304 * (x1069 + x1077 * x629 + x110 * x963))
    result[14, 9] = numpy.sum(x173 * (3.0 * x1071 + x1079 * x629 + x110 * x969))
    result[14, 10] = numpy.sum(x112 * (x1080 * x629 + x110 * x972))
    result[14, 11] = numpy.sum(x173 * (x1072 + x1081 * x629 + x110 * x975))
    result[14, 12] = numpy.sum(x272 * (2.0 * x1074 + x1082 * x629 + x110 * x978))
    result[14, 13] = numpy.sum(x173 * (3.0 * x1076 + x1083 * x629 + x110 * x981))
    result[14, 14] = numpy.sum(x112 * (4.0 * x1078 + x1084 * x629 + x110 * x984))
    return result


def _2center2el3d_45(ax, da, A, bx, db, B):
    """Cartesian (g|h) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 21), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = ax ** (-1.0)
    x5 = -x2 - B[0]
    x6 = bx ** (-1.0)
    x7 = ax * x1
    x8 = bx * x7 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x9 = boys(7, x8)
    x10 = x0 ** (-1.5)
    x11 = 17.4934183276249
    x12 = 2.0 * x6
    x13 = x11 * x12
    x14 = x10 * x13
    x15 = x14 * x9
    x16 = x15 * x5
    x17 = x0 ** (-0.5)
    x18 = boys(6, x8)
    x19 = 0.5 * x6
    x20 = x19 * (2.0 * x11 * x17 * x18 * x4 * x6 - x15)
    x21 = boys(8, x8)
    x22 = x5**2
    x23 = x17 * x4
    x24 = x13 * x23
    x25 = x22 * x24
    x26 = x21 * x25
    x27 = x20 + x26
    x28 = x27 * x5 + x6 * (2.0 * x11 * x17 * x18 * x4 * x5 * x6 - x16)
    x29 = x14 * x18
    x30 = boys(5, x8)
    x31 = x19 * (2.0 * x11 * x17 * x30 * x4 * x6 - x29)
    x32 = x24 * x9
    x33 = x22 * x32
    x34 = x31 + x33
    x35 = x14 * x30
    x36 = boys(4, x8)
    x37 = x19 * (2.0 * x11 * x17 * x36 * x4 * x6 - x35)
    x38 = x18 * x24
    x39 = x22 * x38
    x40 = x37 + x39
    x41 = 1.5 * x6
    x42 = x28 * x5 - x41 * (x34 * x7 - x40)
    x43 = x29 * x5
    x44 = x34 * x5 + x6 * (2.0 * x11 * x17 * x30 * x4 * x5 * x6 - x43)
    x45 = x35 * x5
    x46 = x40 * x5 + x6 * (2.0 * x11 * x17 * x36 * x4 * x5 * x6 - x45)
    x47 = -x12 * (x44 * x7 - x46) + x42 * x5
    x48 = x3 * x47
    x49 = 0.5 / (ax + bx)
    x50 = x14 * x36
    x51 = boys(3, x8)
    x52 = x19 * (2.0 * x11 * x17 * x4 * x51 * x6 - x50)
    x53 = x24 * x30
    x54 = x22 * x53
    x55 = x52 + x54
    x56 = -x41 * (x40 * x7 - x55) + x44 * x5
    x57 = x49 * x56
    x58 = 5.0 * x57
    x59 = x48 + x58
    x60 = bx * x1
    x61 = x5 * x50
    x62 = x5 * x55 + x6 * (2.0 * x11 * x17 * x4 * x5 * x51 * x6 - x61)
    x63 = -x12 * (x46 * x7 - x62) + x5 * x56
    x64 = x3 * x63
    x65 = x14 * x51
    x66 = boys(2, x8)
    x67 = x19 * (2.0 * x11 * x17 * x4 * x6 * x66 - x65)
    x68 = x24 * x36
    x69 = x22 * x68
    x70 = x67 + x69
    x71 = -x41 * (x55 * x7 - x70) + x46 * x5
    x72 = x49 * x71
    x73 = 5.0 * x72
    x74 = x64 + x73
    x75 = x14 * x21
    x76 = x5 * x75
    x77 = x19 * (2.0 * x11 * x17 * x4 * x6 * x9 - x75)
    x78 = boys(9, x8)
    x79 = x25 * x78
    x80 = -x12 * (x28 * x7 - x44) - x5 * (
        x41 * (x27 * x7 - x34)
        - x5 * (x5 * (x77 + x79) + x6 * (2.0 * x11 * x17 * x4 * x5 * x6 * x9 - x76))
    )
    x81 = x3 * x80
    x82 = x42 * x49
    x83 = 5.0 * x82
    x84 = x47 * x60
    x85 = 0.5 * x4
    x86 = x85 * (x63 - x84)
    x87 = x3 * x42
    x88 = x44 * x49
    x89 = 4.0 * x88
    x90 = x87 + x89
    x91 = 5.0 * x49
    x92 = x56 * x60
    x93 = x85 * (x71 - x92)
    x94 = x3 * x44
    x95 = x40 * x49
    x96 = 3.0 * x95
    x97 = x94 + x96
    x98 = 4.0 * x49
    x99 = x3 * x90 + x93 + x97 * x98
    x100 = x46 * x49
    x101 = 4.0 * x100 + x3 * x56
    x102 = x49 * x62
    x103 = 4.0 * x102 + x3 * x71
    x104 = x85 * (-x46 * x60 + x62)
    x105 = x11 * x23 * x6 * x98
    x106 = x105 * x30
    x107 = 3.0 * x49
    x108 = x14 * x66
    x109 = boys(1, x8)
    x110 = x19 * (-x108 + 2.0 * x109 * x11 * x17 * x4 * x6)
    x111 = x24 * x51
    x112 = x111 * x22
    x113 = x110 + x112
    x114 = x5 * (x6 * (2.0 * x11 * x17 * x4 * x6 * x66 - x65) + x70)
    x115 = x12 * (x114 - x62 * x7) + x5 * x71
    x116 = x85 * (
        -x115 * x60
        - x12
        * (-x113 * x5 + x114 * x7 + x5 * x6 * (x108 - 2.0 * x109 * x11 * x17 * x4 * x6))
        + x5 * (x41 * (x113 - x7 * x70) + x5 * x62)
    )
    x117 = x85 * (x115 - x60 * x63)
    x118 = 1.5 * x4
    x119 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**6.5)
    x120 = 4.59731646942873 * x119
    x121 = -x1 * (ax * A[1] + bx * B[1])
    x122 = -x121 - B[1]
    x123 = x122 * x5
    x124 = x122 * x15
    x125 = x6 * (2.0 * x11 * x122 * x17 * x18 * x4 * x6 - x124)
    x126 = x122 * x26
    x127 = 0.5 * x125 + x126
    x128 = x127 * x5 + x6 * (2.0 * x11 * x122 * x17 * x18 * x4 * x5 * x6 - x123 * x15)
    x129 = x122 * x29
    x130 = x6 * (2.0 * x11 * x122 * x17 * x30 * x4 * x6 - x129)
    x131 = x122 * x33
    x132 = 0.5 * x130 + x131
    x133 = x122 * x35
    x134 = x6 * (2.0 * x11 * x122 * x17 * x36 * x4 * x6 - x133)
    x135 = x122 * x39
    x136 = 0.5 * x134 + x135
    x137 = x128 * x5 - x41 * (x132 * x7 - x136)
    x138 = x137 * x3
    x139 = x132 * x5 + x6 * (2.0 * x11 * x122 * x17 * x30 * x4 * x5 * x6 - x123 * x29)
    x140 = x139 * x49
    x141 = 4.0 * x140
    x142 = x138 + x141
    x143 = x122 * x50
    x144 = x6 * (2.0 * x11 * x122 * x17 * x4 * x51 * x6 - x143)
    x145 = x122 * x54
    x146 = 0.5 * x144 + x145
    x147 = x139 * x5 - x41 * (x136 * x7 - x146)
    x148 = x147 * x3
    x149 = x136 * x5 + x6 * (2.0 * x11 * x122 * x17 * x36 * x4 * x5 * x6 - x123 * x35)
    x150 = x149 * x49
    x151 = x148 + 4.0 * x150
    x152 = x122 * x75
    x153 = x6 * (2.0 * x11 * x122 * x17 * x4 * x6 * x9 - x152)
    x154 = x122 * x79
    x155 = -x41 * (x127 * x7 - x132) + 0.5 * x5 * (
        x5 * (x153 + 2.0 * x154)
        + 2.0 * x6 * (2.0 * x11 * x122 * x17 * x4 * x5 * x6 * x9 - x123 * x75)
    )
    x156 = x155 * x3
    x157 = x128 * x49
    x158 = 4.0 * x157
    x159 = x137 * x60
    x160 = x85 * (x147 - x159)
    x161 = x128 * x3
    x162 = x132 * x49
    x163 = 3.0 * x162
    x164 = x161 + x163
    x165 = x139 * x60
    x166 = x85 * (x149 - x165)
    x167 = x132 * x3
    x168 = x123 * x18
    x169 = x105 * x168
    x170 = x167 + x169
    x171 = x107 * x170 + x164 * x3 + x166
    x172 = x136 * x49
    x173 = x139 * x3 + 3.0 * x172
    x174 = x146 * x49
    x175 = x149 * x3 + 3.0 * x174
    x176 = x85 * (-x136 * x60 + x146)
    x177 = x122 * x49 * x53
    x178 = x3 * x5
    x179 = x122 * x178
    x180 = 2.0 * x49
    x181 = x122 * x6 * (-x108 + 2.0 * x109 * x11 * x17 * x4 * x6)
    x182 = x122 * x6 * (2.0 * x11 * x17 * x4 * x6 * x66 - x65)
    x183 = x122 * x69 + 0.5 * x182
    x184 = x149 * x5 - x41 * (x146 * x7 - x183)
    x185 = (
        0.5
        * x85
        * (
            -2.0 * x184 * x60
            + x41 * (2.0 * x112 * x122 + x181 - 2.0 * x183 * x7)
            + 2.0
            * x5
            * (
                x146 * x5
                + x6 * (2.0 * x11 * x122 * x17 * x4 * x5 * x51 * x6 - x123 * x50)
            )
        )
    )
    x186 = x85 * (-x147 * x60 + x184)
    x187 = 13.7919494082862 * x119
    x188 = -x1 * (ax * A[2] + bx * B[2])
    x189 = -x188 - B[2]
    x190 = x15 * x189
    x191 = x6 * (2.0 * x11 * x17 * x18 * x189 * x4 * x6 - x190)
    x192 = 0.5 * x191
    x193 = x189 * x26 + x192
    x194 = x189 * x6 * (2.0 * x11 * x17 * x18 * x4 * x5 * x6 - x16) + x193 * x5
    x195 = x189 * x29
    x196 = x6 * (2.0 * x11 * x17 * x189 * x30 * x4 * x6 - x195)
    x197 = 0.5 * x196
    x198 = x189 * x33 + x197
    x199 = x189 * x35
    x200 = x6 * (2.0 * x11 * x17 * x189 * x36 * x4 * x6 - x199)
    x201 = 0.5 * x200
    x202 = x189 * x39 + x201
    x203 = x194 * x5 - x41 * (x198 * x7 - x202)
    x204 = x203 * x3
    x205 = x189 * x6 * (2.0 * x11 * x17 * x30 * x4 * x5 * x6 - x43) + x198 * x5
    x206 = x205 * x49
    x207 = 4.0 * x206
    x208 = x204 + x207
    x209 = x189 * x50
    x210 = x6 * (2.0 * x11 * x17 * x189 * x4 * x51 * x6 - x209)
    x211 = 0.5 * x210
    x212 = x189 * x54 + x211
    x213 = x205 * x5 - x41 * (x202 * x7 - x212)
    x214 = x213 * x3
    x215 = x189 * x6 * (2.0 * x11 * x17 * x36 * x4 * x5 * x6 - x45) + x202 * x5
    x216 = x215 * x49
    x217 = x214 + 4.0 * x216
    x218 = x189 * x75
    x219 = x6 * (2.0 * x11 * x17 * x189 * x4 * x6 * x9 - x218)
    x220 = 0.5 * x219
    x221 = -x41 * (x193 * x7 - x198) + x5 * (
        x189 * x6 * (2.0 * x11 * x17 * x4 * x5 * x6 * x9 - x76) + x5 * (x189 * x79 + x220)
    )
    x222 = x221 * x3
    x223 = x194 * x49
    x224 = 4.0 * x223
    x225 = x203 * x60
    x226 = x85 * (x213 - x225)
    x227 = x194 * x3
    x228 = x198 * x49
    x229 = 3.0 * x228
    x230 = x227 + x229
    x231 = x205 * x60
    x232 = x85 * (x215 - x231)
    x233 = x198 * x3
    x234 = x105 * x189
    x235 = x18 * x234
    x236 = x235 * x5
    x237 = x233 + x236
    x238 = x107 * x237 + x230 * x3 + x232
    x239 = x202 * x49
    x240 = x205 * x3 + 3.0 * x239
    x241 = x212 * x49
    x242 = x215 * x3 + 3.0 * x241
    x243 = x85 * (-x202 * x60 + x212)
    x244 = x189 * x53
    x245 = x244 * x49
    x246 = x189 * x38
    x247 = x189 * x6 * (-x108 + 2.0 * x109 * x11 * x17 * x4 * x6)
    x248 = 0.5 * x247
    x249 = x189 * x6 * (2.0 * x11 * x17 * x4 * x6 * x66 - x65)
    x250 = 0.5 * x249
    x251 = x189 * x69 + x250
    x252 = x215 * x5 - x41 * (x212 * x7 - x251)
    x253 = x85 * (
        -x252 * x60
        + x41 * (x112 * x189 + x248 - x251 * x7)
        + x5 * (x189 * x6 * (2.0 * x11 * x17 * x4 * x5 * x51 * x6 - x61) + x212 * x5)
    )
    x254 = x85 * (-x213 * x60 + x252)
    x255 = x122**2
    x256 = x255 * x38 + x37
    x257 = x255 * x32 + x31
    x258 = x257 * x5
    x259 = x24 * x255
    x260 = x21 * x259
    x261 = x20 + x260
    x262 = x256 - x257 * x7
    x263 = x19 * x262 + x22 * x261
    x264 = x263 * x5 + x6 * (x256 * x5 - x258 * x7)
    x265 = x264 * x3
    x266 = x255 * x53 + x52
    x267 = -x256 * x7 + x266
    x268 = x19 * x267 + x22 * x257
    x269 = x268 * x49
    x270 = 3.0 * x269
    x271 = x265 + x270
    x272 = -x266 * x5
    x273 = x256 * x5
    x274 = x268 * x5 - x6 * (x272 + x273 * x7)
    x275 = x274 * x3
    x276 = x255 * x68 + x67
    x277 = -x266 * x7 + x276
    x278 = x19 * x277 + x22 * x256
    x279 = x278 * x49
    x280 = 3.0 * x279
    x281 = x275 + x280
    x282 = x261 * x7
    x283 = x259 * x78
    x284 = x283 + x77
    x285 = x257 - x282
    x286 = x5 * (x19 * x285 + x22 * x284 + x6 * (x257 - x282))
    x287 = x286 * x3
    x288 = x263 * x49
    x289 = 3.0 * x288
    x290 = x264 * x60
    x291 = x85 * (x274 - x290)
    x292 = x263 * x3
    x293 = x180 * x258
    x294 = x292 + x293
    x295 = x268 * x60
    x296 = x85 * (x278 - x295)
    x297 = x256 * x49
    x298 = x258 * x3
    x299 = x297 + x298
    x300 = x180 * x299 + x294 * x3 + x296
    x301 = x180 * x273 + x268 * x3
    x302 = x266 * x5
    x303 = x180 * x302 + x278 * x3
    x304 = -x85 * (x272 + x273 * x60)
    x305 = x110 + x111 * x255
    x306 = -x276 * x7 + x305
    x307 = x278 * x5 + x6 * (x276 * x5 - x302 * x7)
    x308 = -x85 * (
        x307 * x60 + x5 * x6 * (x276 * x7 - x305) - x5 * (x19 * x306 + x22 * x266)
    )
    x309 = x85 * (-x274 * x60 + x307)
    x310 = 21.0675507148243 * x119
    x311 = x189 * x6 * (2.0 * x11 * x122 * x17 * x18 * x4 * x6 - x124)
    x312 = x126 * x189 + 0.5 * x311
    x313 = x312 * x5 + x6 * (
        2.0 * x11 * x122 * x17 * x18 * x189 * x4 * x5 * x6 - x123 * x190
    )
    x314 = x189 * x6 * (2.0 * x11 * x122 * x17 * x30 * x4 * x6 - x129)
    x315 = x131 * x189 + 0.5 * x314
    x316 = x107 * x315 + x3 * x313
    x317 = -2.0 * x11 * x122 * x17 * x189 * x30 * x4 * x5 * x6
    x318 = x315 * x5 - x6 * (x123 * x195 + x317)
    x319 = x189 * x6 * (2.0 * x11 * x122 * x17 * x36 * x4 * x6 - x133)
    x320 = x135 * x189 + 0.5 * x319
    x321 = x107 * x320 + x3 * x318
    x322 = x189 * x6 * (2.0 * x11 * x122 * x17 * x4 * x6 * x9 - x152)
    x323 = 0.5 * x5 * (2.0 * x154 * x189 + x322) + x6 * (
        2.0 * x11 * x122 * x17 * x189 * x4 * x5 * x6 * x9 - x123 * x218
    )
    x324 = x85 * (-x313 * x60 + x318)
    x325 = x123 * x234 * x9
    x326 = x3 * x312 + x325
    x327 = x85 * (-x315 * x60 + x320)
    x328 = x122 * x189 * x38
    x329 = x328 * x49
    x330 = x189 * x32
    x331 = x179 * x330 + x329
    x332 = x180 * x331 + x3 * x326 + x327
    x333 = x168 * x189
    x334 = x105 * x333
    x335 = x3 * x315 + x334
    x336 = x106 * x189
    x337 = x123 * x336
    x338 = x3 * x320 + x337
    x339 = -x85 * (2.0 * x10 * x11 * x333 * x4 + x317)
    x340 = x189 * x6 * (2.0 * x11 * x122 * x17 * x4 * x51 * x6 - x143)
    x341 = x320 * x5 + x6 * (
        2.0 * x11 * x122 * x17 * x189 * x36 * x4 * x5 * x6 - x123 * x199
    )
    x342 = (
        0.5
        * x85
        * (
            -2.0 * x341 * x60
            + x5 * (2.0 * x145 * x189 + x340)
            + 2.0
            * x6
            * (2.0 * x11 * x122 * x17 * x189 * x4 * x5 * x51 * x6 - x123 * x209)
        )
    )
    x343 = x85 * (-x318 * x60 + x341)
    x344 = 36.4900682291097 * x119
    x345 = x189**2
    x346 = x345 * x38 + x37
    x347 = x31 + x32 * x345
    x348 = x347 * x5
    x349 = x24 * x345
    x350 = x20 + x21 * x349
    x351 = x22 * x350
    x352 = x346 - x347 * x7
    x353 = x19 * x352
    x354 = x351 + x353
    x355 = x354 * x5 + x6 * (x346 * x5 - x348 * x7)
    x356 = x3 * x355
    x357 = x22 * x347
    x358 = x345 * x53 + x52
    x359 = -x346 * x7 + x358
    x360 = x19 * x359
    x361 = x357 + x360
    x362 = x361 * x49
    x363 = 3.0 * x362
    x364 = x356 + x363
    x365 = -x358 * x5
    x366 = x346 * x5
    x367 = x361 * x5 - x6 * (x365 + x366 * x7)
    x368 = x3 * x367
    x369 = x22 * x346
    x370 = x345 * x68 + x67
    x371 = -x358 * x7 + x370
    x372 = x19 * x371
    x373 = x369 + x372
    x374 = x373 * x49
    x375 = 3.0 * x374
    x376 = x368 + x375
    x377 = x350 * x7
    x378 = x349 * x78 + x77
    x379 = x22 * x378
    x380 = x347 - x377
    x381 = x19 * x380
    x382 = x5 * (x379 + x381 + x6 * (x347 - x377))
    x383 = x3 * x382
    x384 = x354 * x49
    x385 = 3.0 * x384
    x386 = x355 * x60
    x387 = x85 * (x367 - x386)
    x388 = x3 * x354
    x389 = x180 * x348
    x390 = x388 + x389
    x391 = x361 * x60
    x392 = x85 * (x373 - x391)
    x393 = x346 * x49
    x394 = x3 * x348
    x395 = x393 + x394
    x396 = x180 * x395 + x3 * x390 + x392
    x397 = x180 * x366 + x3 * x361
    x398 = x358 * x5
    x399 = x180 * x398 + x3 * x373
    x400 = -x85 * (x365 + x366 * x60)
    x401 = x110 + x111 * x345
    x402 = x22 * x358
    x403 = -x370 * x7 + x401
    x404 = x19 * x403
    x405 = x373 * x5 + x6 * (x370 * x5 - x398 * x7)
    x406 = -x85 * (x405 * x60 + x5 * x6 * (x370 * x7 - x401) - x5 * (x402 + x404))
    x407 = x85 * (-x367 * x60 + x405)
    x408 = x122 * x261 + x125
    x409 = x122 * x257 + x130
    x410 = x122 * x256 + x134
    x411 = -x409 * x7 + x410
    x412 = x19 * x411 + x22 * x408
    x413 = x3 * x412
    x414 = x409 * x49
    x415 = x414 * x5
    x416 = x413 + 2.0 * x415
    x417 = x122 * x266 + x144
    x418 = -x410 * x7 + x417
    x419 = x19 * x418 + x22 * x409
    x420 = x3 * x419
    x421 = x410 * x5
    x422 = x180 * x421 + x420
    x423 = x122 * x284 + x153
    x424 = -x408 * x7 + x409
    x425 = x19 * x424 + x22 * x423
    x426 = x3 * x425
    x427 = x408 * x5
    x428 = x412 * x60
    x429 = x85 * (x419 - x428)
    x430 = x178 * x408
    x431 = x414 + x430
    x432 = x5 * x60
    x433 = x85 * (-x409 * x432 + x410 * x5)
    x434 = x3 * x414
    x435 = x3 * x431 + x433 + x434
    x436 = x410 * x49
    x437 = x178 * x409 + x436
    x438 = x417 * x49
    x439 = x3 * x421 + x438
    x440 = x3**2
    x441 = x85 * (-x410 * x60 + x417)
    x442 = x122 * x276 + x182
    x443 = x122 * x305 + x181 - x442 * x7
    x444 = -x417 * x7 + x442
    x445 = x19 * x444 + x22 * x410
    x446 = x85 * (x19 * x443 + x22 * x417 - x445 * x60)
    x447 = x85 * (-x419 * x60 + x445)
    x448 = 21.0675507148243 * x119
    x449 = x189 * x260 + x192
    x450 = x197 + x255 * x330
    x451 = x201 + x246 * x255
    x452 = -x450 * x7 + x451
    x453 = x19 * x452 + x22 * x449
    x454 = x450 * x49
    x455 = 2.0 * x454
    x456 = x455 * x5
    x457 = x3 * x453 + x456
    x458 = x211 + x244 * x255
    x459 = -x451 * x7 + x458
    x460 = x19 * x459 + x22 * x450
    x461 = x451 * x5
    x462 = x180 * x461
    x463 = x3 * x460 + x462
    x464 = x189 * x283 + x220
    x465 = -x449 * x7 + x450
    x466 = x19 * x465 + x22 * x464
    x467 = x449 * x5
    x468 = x180 * x467
    x469 = x85 * (-x453 * x60 + x460)
    x470 = x178 * x449 + x454
    x471 = x85 * (-x432 * x450 + x451 * x5)
    x472 = x3 * x454 + x3 * x470 + x471
    x473 = x451 * x49
    x474 = x178 * x450 + x473
    x475 = x458 * x49
    x476 = x3 * x461 + x475
    x477 = x85 * (-x451 * x60 + x458)
    x478 = x189 * x255 * x68 + x250
    x479 = x111 * x189 * x255 + x248 - x478 * x7
    x480 = -x458 * x7 + x478
    x481 = x19 * x480 + x22 * x451
    x482 = x85 * (x19 * x479 + x22 * x458 - x481 * x60)
    x483 = x85 * (-x460 * x60 + x481)
    x484 = 47.1084755177714 * x119
    x485 = x122 * x347
    x486 = x6 * (x122 * x346 - x485 * x7)
    x487 = x122 * x351 + 0.5 * x486
    x488 = x122 * x389
    x489 = x3 * x487 + x488
    x490 = -x122 * x358
    x491 = x122 * x346
    x492 = -x6 * (x490 + x491 * x7)
    x493 = x122 * x357 + 0.5 * x492
    x494 = x122 * x366
    x495 = x180 * x494
    x496 = x3 * x493 + x495
    x497 = x122 * x6 * (x347 - x377)
    x498 = x122 * x379 + 0.5 * x497
    x499 = x123 * x350
    x500 = x180 * x499
    x501 = x85 * (-x487 * x60 + x493)
    x502 = x485 * x49
    x503 = x179 * x350 + x502
    x504 = x122 * x60
    x505 = x85 * (x122 * x346 * x5 - x348 * x504)
    x506 = x3 * x502 + x3 * x503 + x505
    x507 = x122 * x393
    x508 = x122 * x394 + x507
    x509 = x122 * x358
    x510 = x3 * x494 + x49 * x509
    x511 = -x85 * (x490 + x491 * x60)
    x512 = x122 * x6 * (-x370 * x7 + x401)
    x513 = x6 * (x122 * x370 - x509 * x7)
    x514 = x122 * x369 + 0.5 * x513
    x515 = 0.5 * x85 * (2.0 * x122 * x402 + x512 - 2.0 * x514 * x60)
    x516 = x85 * (-x493 * x60 + x514)
    x517 = x189 * x350 + x191
    x518 = x189 * x347 + x196
    x519 = x189 * x346 + x200
    x520 = -x518 * x7 + x519
    x521 = x19 * x520
    x522 = x22 * x517 + x521
    x523 = x3 * x522
    x524 = x49 * x518
    x525 = x5 * x524
    x526 = x523 + 2.0 * x525
    x527 = x189 * x358 + x210
    x528 = -x519 * x7 + x527
    x529 = x19 * x528
    x530 = x22 * x518 + x529
    x531 = x3 * x530
    x532 = x5 * x519
    x533 = x180 * x532 + x531
    x534 = x189 * x378 + x219
    x535 = -x517 * x7 + x518
    x536 = x19 * x535
    x537 = x22 * x534 + x536
    x538 = x3 * x537
    x539 = x5 * x517
    x540 = x522 * x60
    x541 = x85 * (x530 - x540)
    x542 = x178 * x517
    x543 = x524 + x542
    x544 = x85 * (-x432 * x518 + x5 * x519)
    x545 = x3 * x524
    x546 = x3 * x543 + x544 + x545
    x547 = x49 * x519
    x548 = x178 * x518
    x549 = x547 + x548
    x550 = x49 * x527
    x551 = x3 * x532 + x550
    x552 = x85 * (-x519 * x60 + x527)
    x553 = x189 * x370 + x249
    x554 = x189 * x401 + x247 - x553 * x7
    x555 = x19 * x554
    x556 = -x527 * x7 + x553
    x557 = x19 * x556
    x558 = x22 * x519 + x557
    x559 = x85 * (x22 * x527 + x555 - x558 * x60)
    x560 = x85 * (-x530 * x60 + x558)
    x561 = x122 * x410 + x277 * x41
    x562 = x122 * x409 + x267 * x41
    x563 = x562 * x60
    x564 = x122 * x408 + x262 * x41
    x565 = x440 * x564
    x566 = x85 * (x561 - x563)
    x567 = x565 + x566
    x568 = x49 * x562
    x569 = x178 * x564
    x570 = x568 + x569
    x571 = x49 * x561
    x572 = x5 * x562
    x573 = x3 * x572
    x574 = x571 + x573
    x575 = x49 * x564
    x576 = x122 * x423 + x285 * x41
    x577 = x178 * x576
    x578 = x85 * (-x432 * x564 + x5 * x562)
    x579 = x3 * x575
    x580 = x122 * x417 + x306 * x41
    x581 = x5 * x85 * (-x561 * x60 + x580)
    x582 = x85 * (x5 * x561 - x572 * x60)
    x583 = x122 * x451 + x319
    x584 = x122 * x450 + x314
    x585 = x584 * x60
    x586 = x122 * x449 + x311
    x587 = x85 * (x583 - x585)
    x588 = x440 * x586 + x587
    x589 = x49 * x584
    x590 = x178 * x586 + x589
    x591 = x49 * x583
    x592 = x5 * x584
    x593 = x3 * x592 + x591
    x594 = x49 * x586
    x595 = x122 * x464 + x322
    x596 = x85 * (-x432 * x586 + x5 * x584)
    x597 = x122 * x458 + x340
    x598 = x5 * x85 * (-x583 * x60 + x597)
    x599 = x85 * (x5 * x583 - x592 * x60)
    x600 = x255 * x346 + x372
    x601 = x255 * x347 + x360
    x602 = x60 * x601
    x603 = x255 * x350 + x353
    x604 = x85 * (x600 - x602)
    x605 = x440 * x603 + x604
    x606 = x49 * x601
    x607 = x178 * x603 + x606
    x608 = x49 * x600
    x609 = x5 * x601
    x610 = x3 * x609 + x608
    x611 = x49 * x603
    x612 = x255 * x378 + x381
    x613 = x85 * (-x432 * x603 + x5 * x601)
    x614 = x255 * x358 + x404
    x615 = x5 * x85 * (-x60 * x600 + x614)
    x616 = x85 * (x5 * x600 - x60 * x609)
    x617 = x504 * x518
    x618 = x85 * (x122 * x519 - x617)
    x619 = x122 * x517
    x620 = x440 * x619 + x618
    x621 = x122 * x524
    x622 = x122 * x542 + x621
    x623 = x122 * x547
    x624 = x122 * x548 + x623
    x625 = x49 * x619
    x626 = x85 * (x122 * x5 * x518 - x432 * x619)
    x627 = x122 * x85 * (x5 * x527 - x532 * x60)
    x628 = x123 * x518
    x629 = x85 * (x122 * x5 * x519 - x60 * x628)
    x630 = x189 * x519 + x371 * x41
    x631 = x189 * x518 + x359 * x41
    x632 = x60 * x631
    x633 = x189 * x517 + x352 * x41
    x634 = x440 * x633
    x635 = x85 * (x630 - x632)
    x636 = x634 + x635
    x637 = x49 * x631
    x638 = x178 * x633
    x639 = x637 + x638
    x640 = x49 * x630
    x641 = x5 * x631
    x642 = x3 * x641
    x643 = x640 + x642
    x644 = x49 * x633
    x645 = x189 * x534 + x380 * x41
    x646 = x178 * x645
    x647 = x85 * (-x432 * x633 + x5 * x631)
    x648 = x3 * x644
    x649 = x189 * x527 + x403 * x41
    x650 = x5 * x85 * (-x60 * x630 + x649)
    x651 = x85 * (x5 * x630 - x60 * x641)
    x652 = x12 * x418 + x122 * x562
    x653 = x12 * x411 + x122 * x564
    x654 = x60 * x653
    x655 = x3 * x654
    x656 = x12 * x424 + x122 * x576
    x657 = x440 * x656
    x658 = x85 * (x652 - x654)
    x659 = x12 * x444 + x122 * x561
    x660 = x85 * (x12 * x443 + x122 * x580 - x60 * x659)
    x661 = x85 * (-x60 * x652 + x659)
    x662 = x122 * x584 + x41 * x459
    x663 = x122 * x586 + x41 * x452
    x664 = x60 * x663
    x665 = x122 * x595 + x41 * x465
    x666 = x85 * (x662 - x664)
    x667 = x122 * x583 + x41 * x480
    x668 = x85 * (x122 * x597 + x41 * x479 - x60 * x667)
    x669 = x85 * (-x60 * x662 + x667)
    x670 = x122 * x601 + x492
    x671 = x122 * x603 + x486
    x672 = x60 * x671
    x673 = x122 * x612 + x497
    x674 = x85 * (x670 - x672)
    x675 = x122 * x600 + x513
    x676 = x85 * (x122 * x614 + x512 - x60 * x675)
    x677 = x85 * (-x60 * x670 + x675)
    x678 = x255 * x518 + x529
    x679 = x255 * x517 + x521
    x680 = x60 * x679
    x681 = x255 * x534 + x536
    x682 = x85 * (x678 - x680)
    x683 = x255 * x519 + x557
    x684 = x85 * (x255 * x527 + x555 - x60 * x683)
    x685 = x85 * (-x60 * x678 + x683)
    x686 = x504 * x633
    x687 = x85 * (x122 * x631 - x686)
    x688 = x122 * x645
    x689 = x122 * x85 * (-x60 * x630 + x649)
    x690 = x122 * x631
    x691 = x85 * (x122 * x630 - x60 * x690)
    x692 = x12 * x528 + x189 * x631
    x693 = x12 * x520 + x189 * x633
    x694 = x60 * x693
    x695 = x3 * x694
    x696 = x12 * x535 + x189 * x645
    x697 = x440 * x696
    x698 = x85 * (x692 - x694)
    x699 = x12 * x556 + x189 * x630
    x700 = x85 * (x12 * x554 + x189 * x649 - x60 * x699)
    x701 = x85 * (-x60 * x692 + x699)
    x702 = -x121 - A[1]
    x703 = x702 * x71
    x704 = x702 * x81
    x705 = x702 * x83
    x706 = x702 * x84
    x707 = x4 * (x63 * x702 - x706)
    x708 = x702 * x87
    x709 = x702 * x89
    x710 = x708 + x709
    x711 = x4 * x702 * (x71 - x92)
    x712 = 12.1633560763699 * x119
    x713 = x147 * x702
    x714 = x713 + x72
    x715 = x102 + x149 * x702
    x716 = x137 * x702
    x717 = x57 + x716
    x718 = x100 + x139 * x702
    x719 = x155 * x702
    x720 = x719 + x82
    x721 = x128 * x702
    x722 = x721 + x88
    x723 = x4 * (-x60 * x717 + x714)
    x724 = x132 * x702
    x725 = x724 + x95
    x726 = x107 * x725 + x3 * x722
    x727 = x4 * (-x60 * x718 + x715)
    x728 = x49 * x5 * x53
    x729 = x123 * x702
    x730 = x180 * (x38 * x729 + x728)
    x731 = x215 * x702
    x732 = x4 * x702 * (x213 - x225)
    x733 = x702 * (x227 + x229)
    x734 = x4 * x702 * (x215 - x231)
    x735 = x5 * x702
    x736 = x274 * x702
    x737 = 2.0 * x150 + x736
    x738 = 2.0 * x174 + x278 * x702
    x739 = x107 * x738
    x740 = x264 * x702
    x741 = 2.0 * x140
    x742 = x740 + x741
    x743 = 2.0 * x172 + x268 * x702
    x744 = x107 * x743
    x745 = x286 * x702
    x746 = 2.0 * x157
    x747 = x745 + x746
    x748 = x263 * x702
    x749 = 2.0 * x162
    x750 = x748 + x749
    x751 = x107 * x750
    x752 = x4 * (-x60 * x742 + x737)
    x753 = x258 * x702
    x754 = x169 + x753
    x755 = x180 * x754 + x3 * x750
    x756 = x4 * (-x60 * x743 + x738)
    x757 = x49 * (x106 * x122 + x256 * x702)
    x758 = 55.7394999246661 * x119
    x759 = x216 + x318 * x702
    x760 = x241 + x320 * x702
    x761 = x206 + x313 * x702
    x762 = x239 + x315 * x702
    x763 = x223 + x323 * x702
    x764 = x228 + x312 * x702
    x765 = x4 * (-x60 * x761 + x759)
    x766 = x189 * x38 * x49 * x5 + x330 * x729
    x767 = x180 * x766
    x768 = x3 * x764 + x767
    x769 = x4 * (-x60 * x762 + x760)
    x770 = x122 * x702
    x771 = x49 * (x245 + x246 * x770)
    x772 = 96.5436458580033 * x119
    x773 = x373 * x702
    x774 = x4 * x702 * (x367 - x386)
    x775 = x702 * (x388 + x389)
    x776 = x4 * x702 * (x373 - x391)
    x777 = x393 * x702
    x778 = x419 * x702
    x779 = x280 + x778
    x780 = x107 * x302 + x421 * x702
    x781 = x412 * x702
    x782 = x270 + x781
    x783 = x409 * x702
    x784 = x107 * x273 + x5 * x783
    x785 = x425 * x702
    x786 = x289 + x785
    x787 = x408 * x735
    x788 = x107 * x258
    x789 = x787 + x788
    x790 = x4 * (-x60 * x782 + x779)
    x791 = 3.0 * x297
    x792 = x783 + x791
    x793 = x49 * x792
    x794 = x3 * x789 + x793
    x795 = x4 * (-x60 * x784 + x780)
    x796 = 55.7394999246661 * x119
    x797 = x180 * x320
    x798 = x460 * x702 + x797
    x799 = x337 + x461 * x702
    x800 = x180 * x315
    x801 = x453 * x702 + x800
    x802 = x450 * x702
    x803 = x334 + x5 * x802
    x804 = x180 * x312
    x805 = x466 * x702 + x804
    x806 = x325 + x449 * x735
    x807 = x4 * (-x60 * x801 + x798)
    x808 = x122 * x235
    x809 = x802 + x808
    x810 = x49 * x809
    x811 = x3 * x806 + x810
    x812 = x4 * (-x60 * x803 + x799)
    x813 = 124.637310863398 * x119
    x814 = x374 + x493 * x702
    x815 = x398 * x49 + x494 * x702
    x816 = x180 * x815
    x817 = x362 + x487 * x702
    x818 = x348 * x770 + x366 * x49
    x819 = x180 * x818
    x820 = x384 + x498 * x702
    x821 = x348 * x49
    x822 = x350 * x729 + x821
    x823 = x180 * x822
    x824 = x4 * (-x60 * x817 + x814)
    x825 = x393 + x485 * x702
    x826 = x49 * x825
    x827 = x3 * x822 + x826
    x828 = x4 * (-x60 * x818 + x815)
    x829 = x524 * x735
    x830 = x532 * x702
    x831 = x517 * x735
    x832 = x4 * x702 * (x530 - x540)
    x833 = x524 * x702
    x834 = x542 * x702 + x833
    x835 = x60 * x735
    x836 = x4 * (x5 * x519 * x702 - x518 * x835)
    x837 = 4.0 * x438 + x561 * x702
    x838 = x49 * x837
    x839 = x572 * x702
    x840 = x421 * x98 + x839
    x841 = 4.0 * x436 + x562 * x702
    x842 = x49 * x841
    x843 = x564 * x702
    x844 = x5 * x843
    x845 = 4.0 * x414
    x846 = x5 * x845 + x844
    x847 = x843 + x845
    x848 = x4 * (-x60 * x841 + x837)
    x849 = x49 * x847
    x850 = x576 * x735
    x851 = x427 * x98 + x850
    x852 = x4 * (-x60 * x846 + x840)
    x853 = 3.0 * x475 + x583 * x702
    x854 = x49 * x853
    x855 = x107 * x461 + x592 * x702
    x856 = 3.0 * x473 + x584 * x702
    x857 = x49 * x856
    x858 = x586 * x702
    x859 = 3.0 * x454
    x860 = x5 * (x858 + x859)
    x861 = x858 + x859
    x862 = x4 * (-x60 * x856 + x853)
    x863 = x49 * x861
    x864 = x107 * x467 + x595 * x735
    x865 = x4 * (-x60 * x860 + x855)
    x866 = x180 * x509 + x600 * x702
    x867 = x49 * x866
    x868 = x495 + x609 * x702
    x869 = 2.0 * x507 + x601 * x702
    x870 = x49 * x869
    x871 = x603 * x702
    x872 = x488 + x5 * x871
    x873 = x180 * x485 + x871
    x874 = x4 * (-x60 * x869 + x866)
    x875 = x49 * x873
    x876 = x500 + x612 * x735
    x877 = x4 * (-x60 * x872 + x868)
    x878 = x122 * x519
    x879 = x550 + x702 * x878
    x880 = x49 * x879
    x881 = x49 * x532 + x628 * x702
    x882 = x518 * x770 + x547
    x883 = x49 * x882
    x884 = x525 + x619 * x735
    x885 = x524 + x619 * x702
    x886 = x4 * (-x60 * x882 + x879)
    x887 = x49 * x885
    x888 = x49 * x539 + x534 * x729
    x889 = x4 * (-x60 * x884 + x881)
    x890 = x640 * x702
    x891 = x637 * x702
    x892 = x4 * x702 * (x630 - x632)
    x893 = x644 * x702
    x894 = x4 * (x5 * x631 * x702 - x633 * x835)
    x895 = x652 * x702
    x896 = 5.0 * x571
    x897 = x895 + x896
    x898 = x653 * x702
    x899 = 5.0 * x568
    x900 = x898 + x899
    x901 = x60 * x900
    x902 = x656 * x702
    x903 = 5.0 * x575
    x904 = x902 + x903
    x905 = x4 * (x897 - x901)
    x906 = 4.0 * x591 + x662 * x702
    x907 = 4.0 * x589 + x663 * x702
    x908 = x60 * x907
    x909 = 4.0 * x594 + x665 * x702
    x910 = x4 * (x906 - x908)
    x911 = 3.0 * x608
    x912 = x670 * x702 + x911
    x913 = 3.0 * x606
    x914 = x671 * x702 + x913
    x915 = x60 * x914
    x916 = 3.0 * x611
    x917 = x673 * x702 + x916
    x918 = x4 * (x912 - x915)
    x919 = 2.0 * x623 + x678 * x702
    x920 = 2.0 * x621 + x679 * x702
    x921 = x60 * x920
    x922 = x180 * x619 + x681 * x702
    x923 = x4 * (x919 - x921)
    x924 = x640 + x690 * x702
    x925 = x633 * x770 + x637
    x926 = x60 * x925
    x927 = x644 + x688 * x702
    x928 = x4 * (x924 - x926)
    x929 = x4 * x702 * (x692 - x694)
    x930 = -x188 - A[2]
    x931 = x4 * x930 * (x63 - x84)
    x932 = 0.5 * x931
    x933 = x930 * (x87 + x89)
    x934 = x4 * x930 * (x71 - x92)
    x935 = 0.5 * x934
    x936 = x149 * x930
    x937 = x4 * x930 * (x147 - x159)
    x938 = 0.5 * x937
    x939 = x930 * (x161 + x163)
    x940 = x4 * x930 * (x149 - x165)
    x941 = 0.5 * x940
    x942 = x123 * x930
    x943 = x105 * x18 * x942
    x944 = x213 * x930 + x72
    x945 = x102 + x215 * x930
    x946 = x49 * x945
    x947 = x203 * x930 + x57
    x948 = x100 + x205 * x930
    x949 = x49 * x948
    x950 = x221 * x930 + x82
    x951 = x3 * x950
    x952 = x194 * x930 + x88
    x953 = x49 * x952
    x954 = 4.0 * x953
    x955 = x60 * x947
    x956 = x4 * (x944 - x955)
    x957 = 0.5 * x956
    x958 = x3 * x952
    x959 = x198 * x930 + x95
    x960 = x49 * x959
    x961 = 3.0 * x960
    x962 = x958 + x961
    x963 = x4 * (-x60 * x948 + x945)
    x964 = 0.5 * x963
    x965 = x5 * x930
    x966 = x49 * (x246 * x965 + x728)
    x967 = 2.0 * x966
    x968 = x270 * x930
    x969 = x278 * x930
    x970 = x107 * x969
    x971 = x289 * x930
    x972 = x4 * x930 * (x274 - x290)
    x973 = 0.5 * x972
    x974 = x930 * (x292 + x293)
    x975 = x4 * x930 * (x278 - x295)
    x976 = 0.5 * x975
    x977 = x297 * x930
    x978 = x150 + x318 * x930
    x979 = x174 + x320 * x930
    x980 = x140 + x313 * x930
    x981 = x172 + x315 * x930
    x982 = x157 + x323 * x930
    x983 = x162 + x312 * x930
    x984 = x4 * (-x60 * x980 + x978)
    x985 = 0.5 * x984
    x986 = x123 * x38 * x49 + x330 * x942
    x987 = x180 * x986
    x988 = x3 * x983 + x987
    x989 = x4 * (-x60 * x981 + x979)
    x990 = 0.5 * x989
    x991 = x49 * (x177 + x328 * x930)
    x992 = 2.0 * x216 + x367 * x930
    x993 = 2.0 * x241 + x373 * x930
    x994 = x49 * x993
    x995 = 3.0 * x994
    x996 = 2.0 * x206 + x355 * x930
    x997 = 2.0 * x239 + x361 * x930
    x998 = x49 * x997
    x999 = 3.0 * x998
    x1000 = 2.0 * x223 + x382 * x930
    x1001 = x1000 * x3
    x1002 = 2.0 * x228 + x354 * x930
    x1003 = x1002 * x49
    x1004 = 3.0 * x1003
    x1005 = x60 * x996
    x1006 = x4 * (-x1005 + x992)
    x1007 = 0.5 * x1006
    x1008 = x1002 * x3
    x1009 = x348 * x930
    x1010 = x1009 + x236
    x1011 = x1010 * x49
    x1012 = 2.0 * x1011
    x1013 = x1008 + x1012
    x1014 = x4 * (-x60 * x997 + x993)
    x1015 = 0.5 * x1014
    x1016 = x49 * (x336 + x346 * x930)
    x1017 = x414 * x930
    x1018 = x1017 * x5
    x1019 = x421 * x930
    x1020 = x408 * x965
    x1021 = x4 * x930 * (x419 - x428)
    x1022 = 0.5 * x1021
    x1023 = x1017 + x430 * x930
    x1024 = x60 * x965
    x1025 = x4 * (-x1024 * x409 + x410 * x5 * x930)
    x1026 = 0.5 * x1025
    x1027 = x279 + x460 * x930
    x1028 = x302 * x49 + x461 * x930
    x1029 = x1028 * x180
    x1030 = x269 + x453 * x930
    x1031 = x450 * x930
    x1032 = x1031 * x5 + x273 * x49
    x1033 = x1032 * x180
    x1034 = x288 + x466 * x930
    x1035 = x258 * x49
    x1036 = x1035 + x449 * x965
    x1037 = x1036 * x180
    x1038 = x4 * (x1027 - x1030 * x60)
    x1039 = 0.5 * x1038
    x1040 = x1031 + x297
    x1041 = x1040 * x49
    x1042 = x1036 * x3 + x1041
    x1043 = x4 * (x1028 - x1032 * x60)
    x1044 = 0.5 * x1043
    x1045 = x493 * x930 + x797
    x1046 = x337 + x494 * x930
    x1047 = x1046 * x180
    x1048 = x487 * x930 + x800
    x1049 = x1009 * x122 + x334
    x1050 = x1049 * x180
    x1051 = x498 * x930 + x804
    x1052 = x325 + x350 * x942
    x1053 = x1052 * x180
    x1054 = x4 * (x1045 - x1048 * x60)
    x1055 = 0.5 * x1054
    x1056 = x485 * x930 + x808
    x1057 = x1056 * x49
    x1058 = x1052 * x3 + x1057
    x1059 = x4 * (x1046 - x1049 * x60)
    x1060 = 0.5 * x1059
    x1061 = x375 + x530 * x930
    x1062 = x107 * x398 + x532 * x930
    x1063 = x1062 * x49
    x1064 = x363 + x522 * x930
    x1065 = x518 * x930
    x1066 = x1065 * x5 + x107 * x366
    x1067 = x1066 * x49
    x1068 = x385 + x537 * x930
    x1069 = x1068 * x3
    x1070 = x107 * x348
    x1071 = x1070 + x539 * x930
    x1072 = x1071 * x49
    x1073 = 2.0 * x1072
    x1074 = x1064 * x60
    x1075 = x4 * (x1061 - x1074)
    x1076 = 0.5 * x1075
    x1077 = x1065 + 3.0 * x393
    x1078 = x1077 * x49
    x1079 = x1071 * x3
    x1080 = x1078 + x1079
    x1081 = x4 * (x1062 - x1066 * x60)
    x1082 = 0.5 * x1081
    x1083 = x571 * x930
    x1084 = x568 * x930
    x1085 = x4 * x930 * (x561 - x563)
    x1086 = 0.5 * x1085
    x1087 = x575 * x930
    x1088 = x4 * (-x1024 * x564 + x5 * x562 * x930)
    x1089 = 0.5 * x1088
    x1090 = x438 + x583 * x930
    x1091 = x1090 * x49
    x1092 = x421 * x49 + x592 * x930
    x1093 = x436 + x584 * x930
    x1094 = x1093 * x49
    x1095 = x586 * x930
    x1096 = x1095 * x5 + x415
    x1097 = x1095 + x414
    x1098 = x4 * (x1090 - x1093 * x60)
    x1099 = 0.5 * x1098
    x1100 = x1097 * x49
    x1101 = x427 * x49 + x595 * x965
    x1102 = x4 * (x1092 - x1096 * x60)
    x1103 = 0.5 * x1102
    x1104 = 2.0 * x475 + x600 * x930
    x1105 = x1104 * x49
    x1106 = x462 + x609 * x930
    x1107 = 2.0 * x473 + x601 * x930
    x1108 = x1107 * x49
    x1109 = x603 * x930
    x1110 = x1109 * x5 + x456
    x1111 = x1109 + x455
    x1112 = x4 * (x1104 - x1107 * x60)
    x1113 = 0.5 * x1112
    x1114 = x1111 * x49
    x1115 = x468 + x612 * x965
    x1116 = x4 * (x1106 - x1110 * x60)
    x1117 = 0.5 * x1116
    x1118 = x107 * x509 + x878 * x930
    x1119 = x1118 * x49
    x1120 = x1065 * x123 + x107 * x494
    x1121 = x1065 * x122 + 3.0 * x507
    x1122 = x1121 * x49
    x1123 = x1070 * x122 + x619 * x965
    x1124 = x107 * x485 + x619 * x930
    x1125 = x4 * (x1118 - x1121 * x60)
    x1126 = 0.5 * x1125
    x1127 = x1124 * x49
    x1128 = x107 * x499 + x534 * x942
    x1129 = x4 * (x1120 - x1123 * x60)
    x1130 = 0.5 * x1129
    x1131 = 4.0 * x550 + x630 * x930
    x1132 = x1131 * x49
    x1133 = x532 * x98 + x641 * x930
    x1134 = 4.0 * x547 + x631 * x930
    x1135 = x1134 * x49
    x1136 = x633 * x930
    x1137 = x1136 * x5 + 4.0 * x525
    x1138 = x1136 + 4.0 * x524
    x1139 = x4 * (x1131 - x1134 * x60)
    x1140 = 0.5 * x1139
    x1141 = x1138 * x49
    x1142 = x539 * x98 + x645 * x965
    x1143 = x1142 * x3
    x1144 = x1137 * x60
    x1145 = x4 * (x1133 - x1144)
    x1146 = 0.5 * x1145
    x1147 = x1141 * x3
    x1148 = x4 * x930 * (x652 - x654)
    x1149 = 0.5 * x1148
    x1150 = x571 + x662 * x930
    x1151 = x568 + x663 * x930
    x1152 = x1151 * x60
    x1153 = x575 + x665 * x930
    x1154 = x4 * (x1150 - x1152)
    x1155 = 0.5 * x1154
    x1156 = 2.0 * x591 + x670 * x930
    x1157 = 2.0 * x589 + x671 * x930
    x1158 = x1157 * x60
    x1159 = 2.0 * x594 + x673 * x930
    x1160 = x4 * (x1156 - x1158)
    x1161 = 0.5 * x1160
    x1162 = x678 * x930 + x911
    x1163 = x679 * x930 + x913
    x1164 = x1163 * x60
    x1165 = x681 * x930 + x916
    x1166 = x4 * (x1162 - x1164)
    x1167 = 0.5 * x1166
    x1168 = 4.0 * x623 + x690 * x930
    x1169 = x1136 * x122 + 4.0 * x621
    x1170 = x1169 * x60
    x1171 = x619 * x98 + x688 * x930
    x1172 = x4 * (x1168 - x1170)
    x1173 = 0.5 * x1172
    x1174 = 5.0 * x640 + x692 * x930
    x1175 = 5.0 * x637 + x693 * x930
    x1176 = x1175 * x60
    x1177 = 5.0 * x644 + x696 * x930
    x1178 = x1177 * x440
    x1179 = x4 * (x1174 - x1176)
    x1180 = 0.5 * x1179
    x1181 = x702**2
    x1182 = x1181 * x80
    x1183 = x1182 + x86
    x1184 = x1181 * x42
    x1185 = x1184 + x93
    x1186 = x1185 * x49
    x1187 = x116 + x1181 * x63 - x60 * (x117 + x1181 * x47)
    x1188 = x49 * (x104 + x1181 * x44)
    x1189 = 15.7028251725905 * x119
    x1190 = x702 * x82
    x1191 = x1190 + x160 + x702 * x720
    x1192 = x702 * x88
    x1193 = x1192 + x166 + x702 * x722
    x1194 = x185 + x49 * x703 - x60 * (x186 + x57 * x702 + x702 * x717) + x702 * x714
    x1195 = x176 + x702 * x725 + x702 * x95
    x1196 = x1181 * x221 + x226
    x1197 = x1181 * x194 + x232
    x1198 = x1197 * x49
    x1199 = x1181 * x213 + x253 - x60 * (x1181 * x203 + x254)
    x1200 = x49 * (x1181 * x198 + x243)
    x1201 = x180 * x722 + x291 + x702 * x747
    x1202 = x180 * x725 + x296 + x702 * x750
    x1203 = x107 * x1202
    x1204 = x180 * x715 + x308 - x60 * (x180 * x718 + x309 + x702 * x742) + x702 * x737
    x1205 = x304 + x702 * x754 + x730
    x1206 = 71.9593849780538 * x119
    x1207 = x223 * x702 + x324 + x702 * x763
    x1208 = x228 * x702 + x327 + x702 * x764
    x1209 = x342 + x49 * x731 - x60 * (x206 * x702 + x343 + x702 * x761) + x702 * x759
    x1210 = x180 * (x246 * x49 * x735 + x339 + x702 * x766)
    x1211 = x1181 * x382 + x387
    x1212 = x1181 * x354 + x392
    x1213 = x1212 * x49
    x1214 = x1181 * x367 + x406 - x60 * (x1181 * x355 + x407)
    x1215 = x49 * (x1181 * x348 + x400)
    x1216 = x429 + x702 * x786 + x751
    x1217 = x107 * x754 + x433 + x702 * x789
    x1218 = x446 - x60 * (x447 + x702 * x782 + x744) + x702 * x779 + x739
    x1219 = x49 * (x441 + x702 * x792 + 3.0 * x757)
    x1220 = 71.9593849780538 * x119
    x1221 = x180 * x764 + x469 + x702 * x805
    x1222 = x471 + x702 * x806 + x767
    x1223 = x180 * x760 + x482 - x60 * (x180 * x762 + x483 + x702 * x801) + x702 * x798
    x1224 = x49 * (x477 + x702 * x809 + 2.0 * x771)
    x1225 = 160.906076430005 * x119
    x1226 = x384 * x702 + x501 + x702 * x820
    x1227 = x505 + x702 * x821 + x702 * x822
    x1228 = x1227 * x180
    x1229 = x49 * x773 + x515 - x60 * (x362 * x702 + x516 + x702 * x817) + x702 * x814
    x1230 = x49 * (x511 + x702 * x825 + x777)
    x1231 = x1181 * x537 + x541
    x1232 = x1181 * x539 + x544
    x1233 = x1232 * x49
    x1234 = x1181 * x530 + x559 - x60 * (x1181 * x522 + x560)
    x1235 = x49 * (x1181 * x518 + x552)
    x1236 = x566 + x702 * x847 + 4.0 * x793
    x1237 = x1236 * x49
    x1238 = x578 + x702 * x851 + x789 * x98
    x1239 = x581 - x60 * (x582 + x702 * x846 + x784 * x98) + x702 * x840 + x780 * x98
    x1240 = x587 + x702 * x861 + 3.0 * x810
    x1241 = x1240 * x49
    x1242 = x107 * x806 + x596 + x702 * x864
    x1243 = x107 * x799 + x598 - x60 * (x107 * x803 + x599 + x702 * x860) + x702 * x855
    x1244 = x604 + x702 * x873 + 2.0 * x826
    x1245 = x1244 * x49
    x1246 = x613 + x702 * x876 + x823
    x1247 = -x60 * (x616 + x702 * x872 + x819) + x615 + x702 * x868 + x816
    x1248 = x618 + x702 * x885 + x833
    x1249 = x1248 * x49
    x1250 = x49 * x831 + x626 + x702 * x888
    x1251 = x49 * x830 - x60 * (x629 + x702 * x884 + x829) + x627 + x702 * x881
    x1252 = x1181 * x633
    x1253 = x1252 + x635
    x1254 = x1253 * x49
    x1255 = x1181 * x5 * x645 + x647
    x1256 = x1181 * x641 - x60 * (x1252 * x5 + x651) + x650
    x1257 = x658 + x702 * x904 + 5.0 * x849
    x1258 = -x60 * (x661 + x702 * x900 + 5.0 * x842) + x660 + x702 * x897 + 5.0 * x838
    x1259 = x666 + x702 * x909 + 4.0 * x863
    x1260 = -x60 * (x669 + x702 * x907 + 4.0 * x857) + x668 + x702 * x906 + 4.0 * x854
    x1261 = x674 + x702 * x917 + 3.0 * x875
    x1262 = -x60 * (x677 + x702 * x914 + 3.0 * x870) + x676 + x702 * x912 + 3.0 * x867
    x1263 = x682 + x702 * x922 + 2.0 * x887
    x1264 = -x60 * (x685 + x702 * x920 + 2.0 * x883) + x684 + x702 * x919 + 2.0 * x880
    x1265 = x687 + x702 * x927 + x893
    x1266 = -x60 * (x691 + x702 * x925 + x891) + x689 + x702 * x924 + x890
    x1267 = x1181 * x696 + x698
    x1268 = x1181 * x692 - x60 * (x1181 * x693 + x701) + x700
    x1269 = x4 * x930 * (x63 * x702 - x706)
    x1270 = 27.1980910212982 * x119
    x1271 = x82 * x930
    x1272 = x1271 + x719 * x930
    x1273 = x88 * x930
    x1274 = x1273 + x721 * x930
    x1275 = x72 * x930
    x1276 = x57 * x930
    x1277 = x4 * (x1275 - x60 * (x1276 + x716 * x930) + x713 * x930)
    x1278 = x930 * x95
    x1279 = x1278 + x724 * x930
    x1280 = 81.5942730638946 * x119
    x1281 = x4 * x702 * (x944 - x955)
    x1282 = x930 * (x745 + x746)
    x1283 = x930 * (x748 + x749)
    x1284 = x107 * x1283
    x1285 = x4 * (x180 * x936 - x60 * x930 * (x740 + x741) + x736 * x930)
    x1286 = x753 * x930 + x943
    x1287 = x702 * x982 + x953
    x1288 = x702 * x983 + x960
    x1289 = x4 * (-x60 * (x702 * x980 + x949) + x702 * x978 + x946)
    x1290 = x180 * (x702 * x986 + x966)
    x1291 = 215.878154934161 * x119
    x1292 = x4 * x702 * (-x1005 + x992)
    x1293 = x785 * x930 + x971
    x1294 = x930 * (x787 + x788)
    x1295 = x4 * (-x60 * (x781 * x930 + x968) + x778 * x930 + x970)
    x1296 = x49 * x930 * (x783 + x791)
    x1297 = 124.637310863398 * x119
    x1298 = x180 * x983
    x1299 = x1034 * x702 + x1298
    x1300 = x1036 * x702 + x987
    x1301 = x180 * x979
    x1302 = x180 * x981
    x1303 = x4 * (x1027 * x702 + x1301 - x60 * (x1030 * x702 + x1302))
    x1304 = 2.0 * x991
    x1305 = x49 * (x1040 * x702 + x1304)
    x1306 = 278.69749962333 * x119
    x1307 = x1003 + x1051 * x702
    x1308 = x1011 + x1052 * x702
    x1309 = x1308 * x180
    x1310 = x4 * (x1045 * x702 - x60 * (x1048 * x702 + x998) + x994)
    x1311 = x49 * (x1016 + x1056 * x702)
    x1312 = x4 * x702 * (x1061 - x1074)
    x1313 = x1078 * x702
    x1314 = x930 * (x843 + x845)
    x1315 = x1314 * x49
    x1316 = x1020 * x98 + x850 * x930
    x1317 = x4 * (x1019 * x98 - x60 * (x844 * x930 + x845 * x965) + x839 * x930)
    x1318 = 3.0 * x1041 + x1097 * x702
    x1319 = x1318 * x49
    x1320 = x1036 * x107 + x1101 * x702
    x1321 = x4 * (x1028 * x107 + x1092 * x702 - x60 * (x1032 * x107 + x1096 * x702))
    x1322 = 2.0 * x1057 + x1111 * x702
    x1323 = x1322 * x49
    x1324 = x1053 + x1115 * x702
    x1325 = x4 * (x1047 + x1106 * x702 - x60 * (x1050 + x1110 * x702))
    x1326 = x1078 + x1124 * x702
    x1327 = x1326 * x49
    x1328 = x1072 + x1128 * x702
    x1329 = x4 * (x1063 + x1120 * x702 - x60 * (x1067 + x1123 * x702))
    x1330 = x1141 * x702
    x1331 = x4 * x702 * (x1133 - x1144)
    x1332 = x930 * (x902 + x903)
    x1333 = x4 * x930 * (-x60 * (x898 + x899) + x895 + x896)
    x1334 = 4.0 * x1100 + x1153 * x702
    x1335 = x4 * (4.0 * x1091 + x1150 * x702 - x60 * (4.0 * x1094 + x1151 * x702))
    x1336 = 3.0 * x1114
    x1337 = x1159 * x702 + x1336
    x1338 = 3.0 * x1105
    x1339 = 3.0 * x1108
    x1340 = x4 * (x1156 * x702 + x1338 - x60 * (x1157 * x702 + x1339))
    x1341 = 2.0 * x1127 + x1165 * x702
    x1342 = x4 * (2.0 * x1119 + x1162 * x702 - x60 * (2.0 * x1122 + x1163 * x702))
    x1343 = x1141 + x1171 * x702
    x1344 = x4 * (x1132 + x1168 * x702 - x60 * (x1135 + x1169 * x702))
    x1345 = x4 * x702 * (x1174 - x1176)
    x1346 = x930**2
    x1347 = x1346 * x80 + x86
    x1348 = x1347 * x3
    x1349 = x1346 * x42 + x93
    x1350 = x1349 * x49
    x1351 = 5.0 * x1350
    x1352 = x116 + x1346 * x63 - x60 * (x117 + x1346 * x47)
    x1353 = x1352 * x85
    x1354 = x49 * (x104 + x1346 * x44)
    x1355 = x1346 * x155 + x160
    x1356 = x128 * x1346 + x166
    x1357 = x1356 * x49
    x1358 = x1346 * x147 + x185 - x60 * (x1346 * x137 + x186)
    x1359 = x1358 * x85
    x1360 = x49 * (x132 * x1346 + x176)
    x1361 = x1271 + x226 + x930 * x950
    x1362 = x1361 * x3
    x1363 = x1273 + x232 + x930 * x952
    x1364 = x1363 * x49
    x1365 = 4.0 * x1364
    x1366 = x1275 + x253 - x60 * (x1276 + x254 + x930 * x947) + x930 * x944
    x1367 = x1366 * x85
    x1368 = x49 * (x1278 + x243 + x930 * x959)
    x1369 = x1346 * x286 + x291
    x1370 = x1346 * x263 + x296
    x1371 = x1370 * x49
    x1372 = 3.0 * x1371
    x1373 = x1346 * x274 + x308 - x60 * (x1346 * x264 + x309)
    x1374 = x1373 * x85
    x1375 = x49 * (x1346 * x258 + x304)
    x1376 = x157 * x930 + x324 + x930 * x982
    x1377 = x162 * x930 + x327 + x930 * x983
    x1378 = x342 + x49 * x936 - x60 * (x140 * x930 + x343 + x930 * x980) + x930 * x978
    x1379 = x1378 * x85
    x1380 = x180 * (x339 + x38 * x49 * x942 + x930 * x986)
    x1381 = x1000 * x930 + x387 + 2.0 * x953
    x1382 = x1381 * x3
    x1383 = x1002 * x930 + x392 + 2.0 * x960
    x1384 = x1383 * x49
    x1385 = 3.0 * x1384
    x1386 = x406 - x60 * (x407 + x930 * x996 + 2.0 * x949) + x930 * x992 + 2.0 * x946
    x1387 = x1386 * x85
    x1388 = x49 * (x1010 * x930 + x400 + x967)
    x1389 = x1346 * x425 + x429
    x1390 = x1346 * x427 + x433
    x1391 = x1390 * x49
    x1392 = x1346 * x419 + x446 - x60 * (x1346 * x412 + x447)
    x1393 = x1392 * x85
    x1394 = x49 * (x1346 * x409 + x441)
    x1395 = x1034 * x930 + x288 * x930 + x469
    x1396 = x1035 * x930 + x1036 * x930 + x471
    x1397 = x1396 * x180
    x1398 = x1027 * x930 + x482 + x49 * x969 - x60 * (x1030 * x930 + x269 * x930 + x483)
    x1399 = x1398 * x85
    x1400 = x49 * (x1040 * x930 + x477 + x977)
    x1401 = x1051 * x930 + x1298 + x501
    x1402 = x1052 * x930 + x505 + x987
    x1403 = x1402 * x180
    x1404 = x1045 * x930 + x1301 + x515 - x60 * (x1048 * x930 + x1302 + x516)
    x1405 = x1404 * x85
    x1406 = x49 * (x1056 * x930 + x1304 + x511)
    x1407 = x1004 + x1068 * x930 + x541
    x1408 = x1407 * x3
    x1409 = 3.0 * x1011 + x1071 * x930 + x544
    x1410 = x1409 * x49
    x1411 = 2.0 * x1410
    x1412 = x1061 * x930 + x559 - x60 * (x1064 * x930 + x560 + x999) + x995
    x1413 = x1412 * x85
    x1414 = x49 * (3.0 * x1016 + x1077 * x930 + x552)
    x1415 = x1346 * x564
    x1416 = x1415 + x566
    x1417 = x1416 * x49
    x1418 = x1346 * x5 * x576 + x578
    x1419 = x1346 * x572 + x581 - x60 * (x1415 * x5 + x582)
    x1420 = x1419 * x85
    x1421 = x1017 + x1097 * x930 + x587
    x1422 = x1421 * x49
    x1423 = x1020 * x49 + x1101 * x930 + x596
    x1424 = x1019 * x49 + x1092 * x930 + x598 - x60 * (x1018 + x1096 * x930 + x599)
    x1425 = x1424 * x85
    x1426 = 2.0 * x1041 + x1111 * x930 + x604
    x1427 = x1426 * x49
    x1428 = x1037 + x1115 * x930 + x613
    x1429 = x1029 + x1106 * x930 - x60 * (x1033 + x1110 * x930 + x616) + x615
    x1430 = x1429 * x85
    x1431 = 3.0 * x1057 + x1124 * x930 + x618
    x1432 = x1431 * x49
    x1433 = x1052 * x107 + x1128 * x930 + x626
    x1434 = (
        x1046 * x107 + x1120 * x930 - x60 * (x1049 * x107 + x1123 * x930 + x629) + x627
    )
    x1435 = x1434 * x85
    x1436 = 4.0 * x1078 + x1138 * x930 + x635
    x1437 = x1436 * x49
    x1438 = 4.0 * x1072 + x1142 * x930 + x647
    x1439 = x1438 * x3
    x1440 = 4.0 * x1063 + x1133 * x930 - x60 * (4.0 * x1067 + x1137 * x930 + x651) + x650
    x1441 = x1440 * x85
    x1442 = x1346 * x656 + x658
    x1443 = x1346 * x652 - x60 * (x1346 * x653 + x661) + x660
    x1444 = x1443 * x85
    x1445 = x1087 + x1153 * x930 + x666
    x1446 = x1083 + x1150 * x930 - x60 * (x1084 + x1151 * x930 + x669) + x668
    x1447 = x1446 * x85
    x1448 = 2.0 * x1100 + x1159 * x930 + x674
    x1449 = 2.0 * x1091 + x1156 * x930 - x60 * (2.0 * x1094 + x1157 * x930 + x677) + x676
    x1450 = x1449 * x85
    x1451 = x1165 * x930 + x1336 + x682
    x1452 = x1162 * x930 + x1338 - x60 * (x1163 * x930 + x1339 + x685) + x684
    x1453 = x1452 * x85
    x1454 = 4.0 * x1127 + x1171 * x930 + x687
    x1455 = 4.0 * x1119 + x1168 * x930 - x60 * (4.0 * x1122 + x1169 * x930 + x691) + x689
    x1456 = x1455 * x85
    x1457 = 5.0 * x1141 + x1177 * x930 + x698
    x1458 = 5.0 * x1132 + x1174 * x930 - x60 * (5.0 * x1135 + x1175 * x930 + x701) + x700
    x1459 = x1458 * x85
    x1460 = x1183 * x702 + x707
    x1461 = x49 * (x1185 * x702 + x711)
    x1462 = x1186 + x1191 * x702 + x723
    x1463 = x1188 + x1193 * x702 + x727
    x1464 = x1196 * x702 + x732
    x1465 = x49 * (x1197 * x702 + x734)
    x1466 = x1193 * x180 + x1201 * x702 + x752
    x1467 = x107 * (x1195 * x180 + x1202 * x702 + x756)
    x1468 = x1198 + x1207 * x702 + x765
    x1469 = x1200 + x1208 * x702 + x769
    x1470 = x1211 * x702 + x774
    x1471 = x49 * (x1212 * x702 + x776)
    x1472 = x1203 + x1216 * x702 + x790
    x1473 = x107 * x1205 + x1217 * x702 + x795
    x1474 = x1208 * x180 + x1221 * x702 + x807
    x1475 = x1210 + x1222 * x702 + x812
    x1476 = x1213 + x1226 * x702 + x824
    x1477 = x180 * (x1215 + x1227 * x702 + x828)
    x1478 = x1231 * x702 + x832
    x1479 = x49 * (x1232 * x702 + x836)
    x1480 = x49 * (4.0 * x1219 + x1236 * x702 + x848)
    x1481 = x1217 * x98 + x1238 * x702 + x852
    x1482 = x49 * (3.0 * x1224 + x1240 * x702 + x862)
    x1483 = x107 * x1222 + x1242 * x702 + x865
    x1484 = x49 * (2.0 * x1230 + x1244 * x702 + x874)
    x1485 = x1228 + x1246 * x702 + x877
    x1486 = x49 * (x1235 + x1248 * x702 + x886)
    x1487 = x1233 + x1250 * x702 + x889
    x1488 = x49 * (x1253 * x702 + x892)
    x1489 = x1255 * x702 + x894
    x1490 = 5.0 * x1237 + x1257 * x702 + x905
    x1491 = x3 * x712
    x1492 = 4.0 * x1241 + x1259 * x702 + x910
    x1493 = x3 * x344
    x1494 = 3.0 * x1245 + x1261 * x702 + x918
    x1495 = x3 * x758
    x1496 = 2.0 * x1249 + x1263 * x702 + x923
    x1497 = x3 * x796
    x1498 = x1254 + x1265 * x702 + x928
    x1499 = x1267 * x702 + x929
    x1500 = x1182 * x930 + x932
    x1501 = x49 * (x1184 * x930 + x935)
    x1502 = x1190 * x930 + x1272 * x702 + x938
    x1503 = x1192 * x930 + x1274 * x702 + x941
    x1504 = x1181 * x950 + x957
    x1505 = x49 * (x1181 * x952 + x964)
    x1506 = x1274 * x180 + x1282 * x702 + x973
    x1507 = x107 * (x1279 * x180 + x1283 * x702 + x976)
    x1508 = x1287 * x702 + x702 * x953 + x985
    x1509 = x1288 * x702 + x702 * x960 + x990
    x1510 = x1000 * x1181 + x1007
    x1511 = x49 * (x1002 * x1181 + x1015)
    x1512 = x1022 + x1284 + x1293 * x702
    x1513 = x1026 + x107 * x1286 + x1294 * x702
    x1514 = x1039 + x1288 * x180 + x1299 * x702
    x1515 = x1044 + x1290 + x1300 * x702
    x1516 = x1003 * x702 + x1055 + x1307 * x702
    x1517 = x180 * (x1011 * x702 + x1060 + x1308 * x702)
    x1518 = x1068 * x1181 + x1076
    x1519 = x49 * (x1071 * x1181 + x1082)
    x1520 = x49 * (x1086 + 4.0 * x1296 + x1314 * x702)
    x1521 = x1089 + x1294 * x98 + x1316 * x702
    x1522 = x49 * (x1099 + 3.0 * x1305 + x1318 * x702)
    x1523 = x107 * x1300 + x1103 + x1320 * x702
    x1524 = x49 * (x1113 + 2.0 * x1311 + x1322 * x702)
    x1525 = x1117 + x1309 + x1324 * x702
    x1526 = x49 * (x1126 + x1313 + x1326 * x702)
    x1527 = x1072 * x702 + x1130 + x1328 * x702
    x1528 = x49 * (x1138 * x1181 + x1140)
    x1529 = x1142 * x1181 + x1146
    x1530 = x1149 + 5.0 * x1315 + x1332 * x702
    x1531 = x1270 * x3
    x1532 = x1155 + 4.0 * x1319 + x1334 * x702
    x1533 = x1280 * x3
    x1534 = x1161 + 3.0 * x1323 + x1337 * x702
    x1535 = x3 * x813
    x1536 = x1167 + 2.0 * x1327 + x1341 * x702
    x1537 = x1297 * x3
    x1538 = x1173 + x1330 + x1343 * x702
    x1539 = x1177 * x1181 + x1180
    x1540 = x1350 + x1355 * x702
    x1541 = x1354 + x1356 * x702
    x1542 = 2.0 * x1357 + x1369 * x702
    x1543 = x107 * (2.0 * x1360 + x1370 * x702)
    x1544 = x1364 + x1376 * x702
    x1545 = x1368 + x1377 * x702
    x1546 = x1372 + x1389 * x702
    x1547 = 3.0 * x1375 + x1390 * x702
    x1548 = x1377 * x180
    x1549 = x1395 * x702 + x1548
    x1550 = x1380 + x1396 * x702
    x1551 = x1384 + x1401 * x702
    x1552 = x180 * (x1388 + x1402 * x702)
    x1553 = x49 * (4.0 * x1394 + x1416 * x702)
    x1554 = 4.0 * x1391 + x1418 * x702
    x1555 = x49 * (3.0 * x1400 + x1421 * x702)
    x1556 = x107 * x1396 + x1423 * x702
    x1557 = x49 * (2.0 * x1406 + x1426 * x702)
    x1558 = x1403 + x1428 * x702
    x1559 = x49 * (x1414 + x1431 * x702)
    x1560 = x1410 + x1433 * x702
    x1561 = x1437 * x702
    x1562 = 5.0 * x1417 + x1442 * x702
    x1563 = 4.0 * x1422 + x1445 * x702
    x1564 = 3.0 * x1427
    x1565 = x1448 * x702 + x1564
    x1566 = 2.0 * x1432 + x1451 * x702
    x1567 = x1437 + x1454 * x702
    x1568 = x1347 * x930 + x931
    x1569 = x49 * (x1349 * x930 + x934)
    x1570 = x1355 * x930 + x937
    x1571 = x49 * (x1356 * x930 + x940)
    x1572 = x1350 + x1361 * x930 + x956
    x1573 = x49 * (x1354 + x1363 * x930 + x963)
    x1574 = x1369 * x930 + x972
    x1575 = x49 * (x1370 * x930 + x975)
    x1576 = 3.0 * x1575
    x1577 = x1357 + x1376 * x930 + x984
    x1578 = x1360 + x1377 * x930 + x989
    x1579 = x1006 + 2.0 * x1364 + x1381 * x930
    x1580 = x49 * (x1014 + 2.0 * x1368 + x1383 * x930)
    x1581 = 3.0 * x1580
    x1582 = x1021 + x1389 * x930
    x1583 = x49 * (x1025 + x1390 * x930)
    x1584 = x1038 + x1371 + x1395 * x930
    x1585 = x1043 + x1375 + x1396 * x930
    x1586 = x1585 * x180
    x1587 = x1054 + x1401 * x930 + x1548
    x1588 = x1059 + x1380 + x1402 * x930
    x1589 = x1588 * x180
    x1590 = x1075 + x1385 + x1407 * x930
    x1591 = x49 * (x1081 + 3.0 * x1388 + x1409 * x930)
    x1592 = x49 * (x1085 + x1416 * x930)
    x1593 = x1088 + x1418 * x930
    x1594 = x49 * (x1098 + x1394 + x1421 * x930)
    x1595 = x1102 + x1391 + x1423 * x930
    x1596 = x49 * (x1112 + 2.0 * x1400 + x1426 * x930)
    x1597 = x1116 + x1397 + x1428 * x930
    x1598 = x49 * (x1125 + 3.0 * x1406 + x1431 * x930)
    x1599 = x107 * x1402 + x1129 + x1433 * x930
    x1600 = x49 * (x1139 + 4.0 * x1414 + x1436 * x930)
    x1601 = x1145 + 4.0 * x1410 + x1438 * x930
    x1602 = x1148 + x1442 * x930
    x1603 = x1154 + x1417 + x1445 * x930
    x1604 = x1160 + 2.0 * x1422 + x1448 * x930
    x1605 = x1166 + x1451 * x930 + x1564
    x1606 = x1172 + 4.0 * x1432 + x1454 * x930
    x1607 = x1179 + 5.0 * x1437 + x1457 * x930
    x1608 = x702 * x712
    x1609 = x344 * x702
    x1610 = x1578 * x180
    x1611 = 3.0 * x1596

    # 315 item(s)
    result[0, 0] = numpy.sum(
        x120
        * (
            x118 * (x103 * x91 + x116 + x3 * x74 - x60 * (x101 * x91 + x117 + x3 * x59))
            + x3
            * (
                x3 * (x3 * (x81 + x83) + x86 + x90 * x91)
                - x4 * (x59 * x60 - x74)
                + x91 * x99
            )
            + x91
            * (
                x3 * x99
                - x4 * (x101 * x60 - x103)
                + x98 * (x104 + x107 * (x106 * x5 + x3 * x40) + x3 * x97)
            )
        )
    )
    result[0, 1] = numpy.sum(
        x187
        * (
            x118 * (x151 * x3 + x175 * x98 + x185 - x60 * (x142 * x3 + x173 * x98 + x186))
            + x3
            * (
                x171 * x98
                + x3 * (x160 + x164 * x98 + x3 * (x156 + x158))
                - x4 * (x142 * x60 - x151)
            )
            + x98
            * (
                x107 * (x170 * x3 + x176 + x180 * (x177 + x179 * x38))
                + x171 * x3
                - x4 * (x173 * x60 - x175)
            )
        )
    )
    result[0, 2] = numpy.sum(
        x187
        * (
            x118 * (x217 * x3 + x242 * x98 + x253 - x60 * (x208 * x3 + x240 * x98 + x254))
            + x3
            * (
                x238 * x98
                + x3 * (x226 + x230 * x98 + x3 * (x222 + x224))
                - x4 * (x208 * x60 - x217)
            )
            + x98
            * (
                x107 * (x180 * (x178 * x246 + x245) + x237 * x3 + x243)
                + x238 * x3
                - x4 * (x240 * x60 - x242)
            )
        )
    )
    result[0, 3] = numpy.sum(
        x310
        * (
            x107
            * (
                x180 * (x297 * x3 + x299 * x3 + x304)
                + x3 * x300
                - x4 * (x301 * x60 - x303)
            )
            + x118
            * (x107 * x303 + x281 * x3 + x308 - x60 * (x107 * x301 + x271 * x3 + x309))
            + x3
            * (
                x107 * x300
                + x3 * (x107 * x294 + x291 + x3 * (x287 + x289))
                - x4 * (x271 * x60 - x281)
            )
        )
    )
    result[0, 4] = numpy.sum(
        x344
        * (
            x107
            * (
                x180 * (x3 * x329 + x3 * x331 + x339)
                + x3 * x332
                - x4 * (x335 * x60 - x338)
            )
            + x118
            * (x107 * x338 + x3 * x321 + x342 - x60 * (x107 * x335 + x3 * x316 + x343))
            + x3
            * (
                x107 * x332
                + x3 * (x107 * x326 + x3 * (x107 * x312 + x3 * x323) + x324)
                - x4 * (x316 * x60 - x321)
            )
        )
    )
    result[0, 5] = numpy.sum(
        x310
        * (
            x107
            * (
                x180 * (x3 * x393 + x3 * x395 + x400)
                + x3 * x396
                - x4 * (x397 * x60 - x399)
            )
            + x118
            * (x107 * x399 + x3 * x376 + x406 - x60 * (x107 * x397 + x3 * x364 + x407))
            + x3
            * (
                x107 * x396
                + x3 * (x107 * x390 + x3 * (x383 + x385) + x387)
                - x4 * (x364 * x60 - x376)
            )
        )
    )
    result[0, 6] = numpy.sum(
        x448
        * (
            x118
            * (x180 * x439 + x3 * x422 + x446 - x60 * (x180 * x437 + x3 * x416 + x447))
            + x180 * (x3 * x435 - x4 * (x437 * x60 - x439) + x49 * (x409 * x440 + x441))
            + x3
            * (
                x180 * x435
                + x3 * (x180 * x431 + x3 * (x180 * x427 + x426) + x429)
                - x4 * (x416 * x60 - x422)
            )
        )
    )
    result[0, 7] = numpy.sum(
        x484
        * (
            x118
            * (x180 * x476 + x3 * x463 + x482 - x60 * (x180 * x474 + x3 * x457 + x483))
            + x180 * (x3 * x472 - x4 * (x474 * x60 - x476) + x49 * (x440 * x450 + x477))
            + x3
            * (
                x180 * x472
                + x3 * (x180 * x470 + x3 * (x3 * x466 + x468) + x469)
                - x4 * (x457 * x60 - x463)
            )
        )
    )
    result[0, 8] = numpy.sum(
        x484
        * (
            x118
            * (x180 * x510 + x3 * x496 + x515 - x60 * (x180 * x508 + x3 * x489 + x516))
            + x180 * (x3 * x506 - x4 * (x508 * x60 - x510) + x49 * (x440 * x485 + x511))
            + x3
            * (
                x180 * x506
                + x3 * (x180 * x503 + x3 * (x3 * x498 + x500) + x501)
                - x4 * (x489 * x60 - x496)
            )
        )
    )
    result[0, 9] = numpy.sum(
        x448
        * (
            x118
            * (x180 * x551 + x3 * x533 + x559 - x60 * (x180 * x549 + x3 * x526 + x560))
            + x180 * (x3 * x546 - x4 * (x549 * x60 - x551) + x49 * (x440 * x518 + x552))
            + x3
            * (
                x180 * x546
                + x3 * (x180 * x543 + x3 * (x180 * x539 + x538) + x541)
                - x4 * (x526 * x60 - x533)
            )
        )
    )
    result[0, 10] = numpy.sum(
        x187
        * (
            x118
            * (x3 * x49 * x561 + x3 * x574 + x581 - x60 * (x3 * x568 + x3 * x570 + x582))
            + x3 * x49 * (x4 * (x561 - x563) + x567)
            + x3
            * (
                x3 * (x3 * (x575 + x577) + x578 + x579)
                - x4 * (x570 * x60 - x574)
                + x49 * x567
            )
        )
    )
    result[0, 11] = numpy.sum(
        x344
        * (
            x118
            * (x3 * x49 * x583 + x3 * x593 + x598 - x60 * (x3 * x589 + x3 * x590 + x599))
            + x3 * x49 * (x4 * (x583 - x585) + x588)
            + x3
            * (
                x3 * (x3 * x594 + x3 * (x178 * x595 + x594) + x596)
                - x4 * (x590 * x60 - x593)
                + x49 * x588
            )
        )
    )
    result[0, 12] = numpy.sum(
        x484
        * (
            x118
            * (x3 * x49 * x600 + x3 * x610 - x60 * (x3 * x606 + x3 * x607 + x616) + x615)
            + x3 * x49 * (x4 * (x600 - x602) + x605)
            + x3
            * (
                x3 * (x3 * x611 + x3 * (x178 * x612 + x611) + x613)
                - x4 * (x60 * x607 - x610)
                + x49 * x605
            )
        )
    )
    result[0, 13] = numpy.sum(
        x344
        * (
            x118 * (x3 * x623 + x3 * x624 - x60 * (x122 * x545 + x3 * x622 + x629) + x627)
            + x3 * x49 * (x4 * (x122 * x519 - x617) + x620)
            + x3
            * (
                x3 * (x3 * x625 + x3 * (x179 * x534 + x625) + x626)
                - x4 * (x60 * x622 - x624)
                + x49 * x620
            )
        )
    )
    result[0, 14] = numpy.sum(
        x187
        * (
            x118
            * (x3 * x49 * x630 + x3 * x643 - x60 * (x3 * x637 + x3 * x639 + x651) + x650)
            + x3 * x49 * (x4 * (x630 - x632) + x636)
            + x3
            * (
                x3 * (x3 * (x644 + x646) + x647 + x648)
                - x4 * (x60 * x639 - x643)
                + x49 * x636
            )
        )
    )
    result[0, 15] = numpy.sum(
        x120
        * (
            x118 * (x440 * x652 - x60 * (x440 * x653 + x661) + x660)
            + x3 * (x3 * (x657 + x658) + x4 * (x3 * x652 - x655))
        )
    )
    result[0, 16] = numpy.sum(
        x187
        * (
            x118 * (x440 * x662 - x60 * (x440 * x663 + x669) + x668)
            + x3**2 * (x4 * (x662 - x664) + x440 * x665 + x666)
        )
    )
    result[0, 17] = numpy.sum(
        x310
        * (
            x118 * (x440 * x670 - x60 * (x440 * x671 + x677) + x676)
            + x3**2 * (x4 * (x670 - x672) + x440 * x673 + x674)
        )
    )
    result[0, 18] = numpy.sum(
        x448
        * (
            x118 * (x440 * x678 - x60 * (x440 * x679 + x685) + x684)
            + x3**2 * (x4 * (x678 - x680) + x440 * x681 + x682)
        )
    )
    result[0, 19] = numpy.sum(
        x187
        * (
            x118 * (x440 * x690 - x60 * (x122 * x634 + x691) + x689)
            + x3**2 * (x4 * (x122 * x631 - x686) + x440 * x688 + x687)
        )
    )
    result[0, 20] = numpy.sum(
        x120
        * (
            x118 * (x440 * x692 - x60 * (x440 * x693 + x701) + x700)
            + x3 * (x3 * (x697 + x698) + x4 * (x3 * x692 - x695))
        )
    )
    result[1, 0] = numpy.sum(
        0.5
        * x712
        * (
            x3 * (2.0 * x3 * (x704 + x705) + x707 + 2.0 * x710 * x91)
            + 2.0 * x4 * (-x60 * x702 * (x48 + x58) + x64 * x702 + x703 * x91)
            + x91 * (2.0 * x3 * x710 + 2.0 * x702 * x98 * (x94 + x96) + x711)
        )
    )
    result[1, 1] = numpy.sum(
        0.5
        * x344
        * (
            x3 * (2.0 * x3 * (x3 * x720 + x722 * x98) + x723 + 2.0 * x726 * x98)
            + 2.0 * x4 * (x3 * x714 - x60 * (x3 * x717 + x718 * x98) + x715 * x98)
            + x98 * (2.0 * x107 * (x3 * x725 + x730) + 2.0 * x3 * x726 + x727)
        )
    )
    result[1, 2] = numpy.sum(
        0.5
        * x344
        * (
            x3 * (2.0 * x3 * x702 * (x222 + x224) + x732 + 2.0 * x733 * x98)
            + 2.0 * x4 * (x214 * x702 - x60 * x702 * (x204 + x207) + x731 * x98)
            + x98 * (2.0 * x107 * (x233 * x702 + x235 * x735) + 2.0 * x3 * x733 + x734)
        )
    )
    result[1, 3] = numpy.sum(
        0.5
        * x758
        * (
            x107 * (2.0 * x180 * (x3 * x754 + x757) + 2.0 * x3 * x755 + x756)
            + x3 * (2.0 * x107 * x755 + 2.0 * x3 * (x3 * x747 + x751) + x752)
            + 2.0 * x4 * (x3 * x737 - x60 * (x3 * x742 + x744) + x739)
        )
    )
    result[1, 4] = numpy.sum(
        0.5
        * x772
        * (
            x107 * (2.0 * x180 * (x3 * x766 + x771) + 2.0 * x3 * x768 + x769)
            + x3 * (2.0 * x107 * x768 + 2.0 * x3 * (x107 * x764 + x3 * x763) + x765)
            + 2.0 * x4 * (x107 * x760 + x3 * x759 - x60 * (x107 * x762 + x3 * x761))
        )
    )
    result[1, 5] = numpy.sum(
        0.5
        * x758
        * (
            x107 * (2.0 * x180 * (x394 * x702 + x777) + 2.0 * x3 * x775 + x776)
            + x3 * (2.0 * x107 * x775 + 2.0 * x3 * x702 * (x383 + x385) + x774)
            + 2.0 * x4 * (x107 * x773 + x368 * x702 - x60 * x702 * (x356 + x363))
        )
    )
    result[1, 6] = numpy.sum(
        0.5
        * x796
        * (
            x180 * (2.0 * x3 * x793 + 2.0 * x3 * x794 + x795)
            + x3 * (2.0 * x180 * x794 + 2.0 * x3 * (x180 * x789 + x3 * x786) + x790)
            + 2.0 * x4 * (x180 * x780 + x3 * x779 - x60 * (x180 * x784 + x3 * x782))
        )
    )
    result[1, 7] = numpy.sum(
        0.5
        * x813
        * (
            x180 * (2.0 * x3 * x810 + 2.0 * x3 * x811 + x812)
            + x3 * (2.0 * x180 * x811 + 2.0 * x3 * (x180 * x806 + x3 * x805) + x807)
            + 2.0 * x4 * (x180 * x799 + x3 * x798 - x60 * (x180 * x803 + x3 * x801))
        )
    )
    result[1, 8] = numpy.sum(
        0.5
        * x813
        * (
            x180 * (2.0 * x3 * x826 + 2.0 * x3 * x827 + x828)
            + x3 * (2.0 * x180 * x827 + 2.0 * x3 * (x3 * x820 + x823) + x824)
            + 2.0 * x4 * (x3 * x814 - x60 * (x3 * x817 + x819) + x816)
        )
    )
    result[1, 9] = numpy.sum(
        0.5
        * x796
        * (
            x180 * (2.0 * x3 * x834 + 2.0 * x545 * x702 + x836)
            + x3 * (2.0 * x180 * x834 + 2.0 * x3 * (x180 * x831 + x538 * x702) + x832)
            + 2.0 * x4 * (x180 * x830 + x531 * x702 - x60 * (x523 * x702 + 2.0 * x829))
        )
    )
    result[1, 10] = numpy.sum(
        0.5
        * x344
        * (
            x3 * (2.0 * x3 * x849 + 2.0 * x3 * (x3 * x851 + x849) + x852)
            + 2.0 * x4 * (x3 * x840 - x60 * (x3 * x846 + x842) + x838)
            + x49 * (2.0 * x440 * x847 + x848)
        )
    )
    result[1, 11] = numpy.sum(
        0.5
        * x772
        * (
            x3 * (2.0 * x3 * x863 + 2.0 * x3 * (x3 * x864 + x863) + x865)
            + 2.0 * x4 * (x3 * x855 - x60 * (x3 * x860 + x857) + x854)
            + x49 * (2.0 * x440 * x861 + x862)
        )
    )
    result[1, 12] = numpy.sum(
        0.5
        * x813
        * (
            x3 * (2.0 * x3 * x875 + 2.0 * x3 * (x3 * x876 + x875) + x877)
            + 2.0 * x4 * (x3 * x868 - x60 * (x3 * x872 + x870) + x867)
            + x49 * (2.0 * x440 * x873 + x874)
        )
    )
    result[1, 13] = numpy.sum(
        0.5
        * x772
        * (
            x3 * (2.0 * x3 * x887 + 2.0 * x3 * (x3 * x888 + x887) + x889)
            + 2.0 * x4 * (x3 * x881 - x60 * (x3 * x884 + x883) + x880)
            + x49 * (2.0 * x440 * x885 + x886)
        )
    )
    result[1, 14] = numpy.sum(
        0.5
        * x344
        * (
            x3 * (2.0 * x3 * (x646 * x702 + x893) + 2.0 * x648 * x702 + x894)
            + 2.0 * x4 * (-x60 * (x638 * x702 + x891) + x642 * x702 + x890)
            + x49 * (2.0 * x634 * x702 + x892)
        )
    )
    result[1, 15] = numpy.sum(
        0.5 * x3 * x712 * (2.0 * x4 * (x897 - x901) + 2.0 * x440 * x904 + x905)
    )
    result[1, 16] = numpy.sum(
        0.5 * x3 * x344 * (2.0 * x4 * (x906 - x908) + 2.0 * x440 * x909 + x910)
    )
    result[1, 17] = numpy.sum(
        0.5 * x3 * x758 * (2.0 * x4 * (x912 - x915) + 2.0 * x440 * x917 + x918)
    )
    result[1, 18] = numpy.sum(
        0.5 * x3 * x796 * (2.0 * x4 * (x919 - x921) + 2.0 * x440 * x922 + x923)
    )
    result[1, 19] = numpy.sum(
        0.5 * x3 * x344 * (2.0 * x4 * (x924 - x926) + 2.0 * x440 * x927 + x928)
    )
    result[1, 20] = numpy.sum(
        0.5
        * x712
        * (x3 * (2.0 * x697 * x702 + x929) + 2.0 * x4 * x702 * (x3 * x692 - x695))
    )
    result[2, 0] = numpy.sum(
        x712
        * (
            x3 * (x3 * x930 * (x81 + x83) + x91 * x933 + x932)
            + x4 * x930 * (-x60 * (x48 + x58) + x64 + x73)
            + x91 * (x3 * x933 + x930 * x98 * (x94 + x96) + x935)
        )
    )
    result[2, 1] = numpy.sum(
        x344
        * (
            x3 * (x3 * x930 * (x156 + x158) + x938 + x939 * x98)
            + x4 * (x148 * x930 - x60 * x930 * (x138 + x141) + x936 * x98)
            + x98 * (x107 * (x167 * x930 + x943) + x3 * x939 + x941)
        )
    )
    result[2, 2] = numpy.sum(
        x344
        * (
            x3 * (x3 * (x951 + x954) + x957 + x962 * x98)
            + x4 * (x3 * x944 - x60 * (x3 * x947 + 4.0 * x949) + 4.0 * x946)
            + x98 * (x107 * (x3 * x959 + x967) + x3 * x962 + x964)
        )
    )
    result[2, 3] = numpy.sum(
        x758
        * (
            x107 * (x180 * (x298 * x930 + x977) + x3 * x974 + x976)
            + x3 * (x107 * x974 + x3 * (x287 * x930 + x971) + x973)
            + x4 * (x275 * x930 - x60 * (x265 * x930 + x968) + x970)
        )
    )
    result[2, 4] = numpy.sum(
        x772
        * (
            x107 * (x180 * (x3 * x986 + x991) + x3 * x988 + x990)
            + x3 * (x107 * x988 + x3 * (x107 * x983 + x3 * x982) + x985)
            + x4 * (x107 * x979 + x3 * x978 - x60 * (x107 * x981 + x3 * x980))
        )
    )
    result[2, 5] = numpy.sum(
        x758
        * (
            x107 * (x1013 * x3 + x1015 + x180 * (x1010 * x3 + x1016))
            + x3 * (x1007 + x1013 * x107 + x3 * (x1001 + x1004))
            + x4 * (x3 * x992 - x60 * (x3 * x996 + x999) + x995)
        )
    )
    result[2, 6] = numpy.sum(
        x796
        * (
            x180 * (x1023 * x3 + x1026 + x434 * x930)
            + x3 * (x1022 + x1023 * x180 + x3 * (x1020 * x180 + x426 * x930))
            + x4 * (x1019 * x180 + x420 * x930 - x60 * (2.0 * x1018 + x413 * x930))
        )
    )
    result[2, 7] = numpy.sum(
        x813
        * (
            x180 * (x1041 * x3 + x1042 * x3 + x1044)
            + x3 * (x1039 + x1042 * x180 + x3 * (x1034 * x3 + x1037))
            + x4 * (x1027 * x3 + x1029 - x60 * (x1030 * x3 + x1033))
        )
    )
    result[2, 8] = numpy.sum(
        x813
        * (
            x180 * (x1057 * x3 + x1058 * x3 + x1060)
            + x3 * (x1055 + x1058 * x180 + x3 * (x1051 * x3 + x1053))
            + x4 * (x1045 * x3 + x1047 - x60 * (x1048 * x3 + x1050))
        )
    )
    result[2, 9] = numpy.sum(
        x796
        * (
            x180 * (x1078 * x3 + x1080 * x3 + x1082)
            + x3 * (x1076 + x1080 * x180 + x3 * (x1069 + x1073))
            + x4 * (x1061 * x3 + 2.0 * x1063 - x60 * (x1064 * x3 + 2.0 * x1067))
        )
    )
    result[2, 10] = numpy.sum(
        x344
        * (
            x3 * (x1089 + x3 * (x1087 + x577 * x930) + x579 * x930)
            + x4 * (x1083 + x573 * x930 - x60 * (x1084 + x569 * x930))
            + x49 * (x1086 + x565 * x930)
        )
    )
    result[2, 11] = numpy.sum(
        x772
        * (
            x3 * (x1100 * x3 + x1103 + x3 * (x1100 + x1101 * x3))
            + x4 * (x1091 + x1092 * x3 - x60 * (x1094 + x1096 * x3))
            + x49 * (x1097 * x440 + x1099)
        )
    )
    result[2, 12] = numpy.sum(
        x813
        * (
            x3 * (x1114 * x3 + x1117 + x3 * (x1114 + x1115 * x3))
            + x4 * (x1105 + x1106 * x3 - x60 * (x1108 + x1110 * x3))
            + x49 * (x1111 * x440 + x1113)
        )
    )
    result[2, 13] = numpy.sum(
        x772
        * (
            x3 * (x1127 * x3 + x1130 + x3 * (x1127 + x1128 * x3))
            + x4 * (x1119 + x1120 * x3 - x60 * (x1122 + x1123 * x3))
            + x49 * (x1124 * x440 + x1126)
        )
    )
    result[2, 14] = numpy.sum(
        x344
        * (
            x3 * (x1146 + x1147 + x3 * (x1141 + x1143))
            + x4 * (x1132 + x1133 * x3 - x60 * (x1135 + x1137 * x3))
            + x49 * (x1138 * x440 + x1140)
        )
    )
    result[2, 15] = numpy.sum(
        x712 * (x3 * (x1149 + x657 * x930) + x4 * x930 * (x3 * x652 - x655))
    )
    result[2, 16] = numpy.sum(x3 * x344 * (x1153 * x440 + x1155 + x4 * (x1150 - x1152)))
    result[2, 17] = numpy.sum(x3 * x758 * (x1159 * x440 + x1161 + x4 * (x1156 - x1158)))
    result[2, 18] = numpy.sum(x3 * x796 * (x1165 * x440 + x1167 + x4 * (x1162 - x1164)))
    result[2, 19] = numpy.sum(x3 * x344 * (x1171 * x440 + x1173 + x4 * (x1168 - x1170)))
    result[2, 20] = numpy.sum(x3 * x712 * (x1178 + x1180 + x4 * (x1174 - x1176)))
    result[3, 0] = numpy.sum(
        x1189
        * (
            x1187 * x85
            + x3 * (x1183 * x3 + 5.0 * x1186)
            + x91 * (x1185 * x3 + 4.0 * x1188)
        )
    )
    result[3, 1] = numpy.sum(
        x484
        * (
            x1194 * x85
            + x3 * (x1191 * x3 + x1193 * x98)
            + x98 * (x107 * x1195 + x1193 * x3)
        )
    )
    result[3, 2] = numpy.sum(
        x484
        * (
            x1199 * x85
            + x3 * (x1196 * x3 + 4.0 * x1198)
            + x98 * (x1197 * x3 + 3.0 * x1200)
        )
    )
    result[3, 3] = numpy.sum(
        x1206
        * (x107 * (x1202 * x3 + x1205 * x180) + x1204 * x85 + x3 * (x1201 * x3 + x1203))
    )
    result[3, 4] = numpy.sum(
        x813
        * (x107 * (x1208 * x3 + x1210) + x1209 * x85 + x3 * (x107 * x1208 + x1207 * x3))
    )
    result[3, 5] = numpy.sum(
        x1206
        * (
            x107 * (x1212 * x3 + 2.0 * x1215)
            + x1214 * x85
            + x3 * (x1211 * x3 + 3.0 * x1213)
        )
    )
    result[3, 6] = numpy.sum(
        x1220
        * (x1218 * x85 + x180 * (x1217 * x3 + x1219) + x3 * (x1216 * x3 + x1217 * x180))
    )
    result[3, 7] = numpy.sum(
        x1225
        * (x1223 * x85 + x180 * (x1222 * x3 + x1224) + x3 * (x1221 * x3 + x1222 * x180))
    )
    result[3, 8] = numpy.sum(
        x1225 * (x1229 * x85 + x180 * (x1227 * x3 + x1230) + x3 * (x1226 * x3 + x1228))
    )
    result[3, 9] = numpy.sum(
        x1220
        * (x1234 * x85 + x180 * (x1232 * x3 + x1235) + x3 * (x1231 * x3 + 2.0 * x1233))
    )
    result[3, 10] = numpy.sum(
        x484 * (x1237 * x3 + x1239 * x85 + x3 * (x1237 + x1238 * x3))
    )
    result[3, 11] = numpy.sum(
        x813 * (x1241 * x3 + x1243 * x85 + x3 * (x1241 + x1242 * x3))
    )
    result[3, 12] = numpy.sum(
        x1225 * (x1245 * x3 + x1247 * x85 + x3 * (x1245 + x1246 * x3))
    )
    result[3, 13] = numpy.sum(
        x813 * (x1249 * x3 + x1251 * x85 + x3 * (x1249 + x1250 * x3))
    )
    result[3, 14] = numpy.sum(
        x484 * (x1254 * x3 + x1256 * x85 + x3 * (x1254 + x1255 * x3))
    )
    result[3, 15] = numpy.sum(x1189 * (x1257 * x440 + x1258 * x85))
    result[3, 16] = numpy.sum(x484 * (x1259 * x440 + x1260 * x85))
    result[3, 17] = numpy.sum(x1206 * (x1261 * x440 + x1262 * x85))
    result[3, 18] = numpy.sum(x1220 * (x1263 * x440 + x1264 * x85))
    result[3, 19] = numpy.sum(x484 * (x1265 * x440 + x1266 * x85))
    result[3, 20] = numpy.sum(x1189 * (x1267 * x440 + x1268 * x85))
    result[4, 0] = numpy.sum(
        0.5
        * x1270
        * (x1269 + 2.0 * x3 * x930 * (x704 + x705) + 2.0 * x91 * x930 * (x708 + x709))
    )
    result[4, 1] = numpy.sum(
        0.5
        * x1280
        * (
            x1277
            + 2.0 * x3 * (x1272 * x3 + x1274 * x98)
            + 2.0 * x98 * (x107 * x1279 + x1274 * x3)
        )
    )
    result[4, 2] = numpy.sum(
        0.5
        * x1280
        * (x1281 + 2.0 * x3 * x702 * (x951 + x954) + 2.0 * x702 * x98 * (x958 + x961))
    )
    result[4, 3] = numpy.sum(
        0.5
        * x813
        * (
            2.0 * x107 * (x1283 * x3 + x1286 * x180)
            + x1285
            + 2.0 * x3 * (x1282 * x3 + x1284)
        )
    )
    result[4, 4] = numpy.sum(
        0.5
        * x1291
        * (
            2.0 * x107 * (x1288 * x3 + x1290)
            + x1289
            + 2.0 * x3 * (x107 * x1288 + x1287 * x3)
        )
    )
    result[4, 5] = numpy.sum(
        0.5
        * x813
        * (
            2.0 * x107 * x702 * (x1008 + x1012)
            + x1292
            + 2.0 * x3 * x702 * (x1001 + x1004)
        )
    )
    result[4, 6] = numpy.sum(
        0.5
        * x1297
        * (
            x1295
            + 2.0 * x180 * (x1294 * x3 + x1296)
            + 2.0 * x3 * (x1293 * x3 + x1294 * x180)
        )
    )
    result[4, 7] = numpy.sum(
        0.5
        * x1306
        * (
            x1303
            + 2.0 * x180 * (x1300 * x3 + x1305)
            + 2.0 * x3 * (x1299 * x3 + x1300 * x180)
        )
    )
    result[4, 8] = numpy.sum(
        0.5
        * x1306
        * (x1310 + 2.0 * x180 * (x1308 * x3 + x1311) + 2.0 * x3 * (x1307 * x3 + x1309))
    )
    result[4, 9] = numpy.sum(
        0.5
        * x1297
        * (
            x1312
            + 2.0 * x180 * (x1079 * x702 + x1313)
            + 2.0 * x3 * x702 * (x1069 + x1073)
        )
    )
    result[4, 10] = numpy.sum(
        0.5 * x1280 * (4.0 * x1315 * x3 + 2.0 * x1316 * x3**2 + x1317)
    )
    result[4, 11] = numpy.sum(
        0.5 * x1291 * (4.0 * x1319 * x3 + 2.0 * x1320 * x3**2 + x1321)
    )
    result[4, 12] = numpy.sum(
        0.5 * x1306 * (4.0 * x1323 * x3 + 2.0 * x1324 * x3**2 + x1325)
    )
    result[4, 13] = numpy.sum(
        0.5 * x1291 * (4.0 * x1327 * x3 + 2.0 * x1328 * x3**2 + x1329)
    )
    result[4, 14] = numpy.sum(
        0.5 * x1280 * (2.0 * x1147 * x702 + x1331 + 2.0 * x3 * (x1143 * x702 + x1330))
    )
    result[4, 15] = numpy.sum(0.5 * x1270 * (2.0 * x1332 * x440 + x1333))
    result[4, 16] = numpy.sum(0.5 * x1280 * (2.0 * x1334 * x440 + x1335))
    result[4, 17] = numpy.sum(0.5 * x813 * (2.0 * x1337 * x440 + x1340))
    result[4, 18] = numpy.sum(0.5 * x1297 * (2.0 * x1341 * x440 + x1342))
    result[4, 19] = numpy.sum(0.5 * x1280 * (2.0 * x1343 * x440 + x1344))
    result[4, 20] = numpy.sum(0.5 * x1270 * (2.0 * x1178 * x702 + x1345))
    result[5, 0] = numpy.sum(
        x1189 * (x1353 + x3 * (x1348 + x1351) + x91 * (x1349 * x3 + 4.0 * x1354))
    )
    result[5, 1] = numpy.sum(
        x484
        * (x1359 + x3 * (x1355 * x3 + 4.0 * x1357) + x98 * (x1356 * x3 + 3.0 * x1360))
    )
    result[5, 2] = numpy.sum(
        x484 * (x1367 + x3 * (x1362 + x1365) + x98 * (x1363 * x3 + 3.0 * x1368))
    )
    result[5, 3] = numpy.sum(
        x1206 * (x107 * (x1370 * x3 + 2.0 * x1375) + x1374 + x3 * (x1369 * x3 + x1372))
    )
    result[5, 4] = numpy.sum(
        x813 * (x107 * (x1377 * x3 + x1380) + x1379 + x3 * (x107 * x1377 + x1376 * x3))
    )
    result[5, 5] = numpy.sum(
        x1206 * (x107 * (x1383 * x3 + 2.0 * x1388) + x1387 + x3 * (x1382 + x1385))
    )
    result[5, 6] = numpy.sum(
        x1220 * (x1393 + x180 * (x1390 * x3 + x1394) + x3 * (x1389 * x3 + 2.0 * x1391))
    )
    result[5, 7] = numpy.sum(
        x1225 * (x1399 + x180 * (x1396 * x3 + x1400) + x3 * (x1395 * x3 + x1397))
    )
    result[5, 8] = numpy.sum(
        x1225 * (x1405 + x180 * (x1402 * x3 + x1406) + x3 * (x1401 * x3 + x1403))
    )
    result[5, 9] = numpy.sum(
        x1220 * (x1413 + x180 * (x1409 * x3 + x1414) + x3 * (x1408 + x1411))
    )
    result[5, 10] = numpy.sum(x484 * (x1417 * x3 + x1420 + x3 * (x1417 + x1418 * x3)))
    result[5, 11] = numpy.sum(x813 * (x1422 * x3 + x1425 + x3 * (x1422 + x1423 * x3)))
    result[5, 12] = numpy.sum(x1225 * (x1427 * x3 + x1430 + x3 * (x1427 + x1428 * x3)))
    result[5, 13] = numpy.sum(x813 * (x1432 * x3 + x1435 + x3 * (x1432 + x1433 * x3)))
    result[5, 14] = numpy.sum(x484 * (x1437 * x3 + x1441 + x3 * (x1437 + x1439)))
    result[5, 15] = numpy.sum(x1189 * (x1442 * x440 + x1444))
    result[5, 16] = numpy.sum(x484 * (x1445 * x440 + x1447))
    result[5, 17] = numpy.sum(x1206 * (x1448 * x440 + x1450))
    result[5, 18] = numpy.sum(x1220 * (x1451 * x440 + x1453))
    result[5, 19] = numpy.sum(x484 * (x1454 * x440 + x1456))
    result[5, 20] = numpy.sum(x1189 * (x1457 * x440 + x1459))
    result[6, 0] = numpy.sum(x712 * (x1460 * x3 + 5.0 * x1461))
    result[6, 1] = numpy.sum(x344 * (x1462 * x3 + x1463 * x98))
    result[6, 2] = numpy.sum(x344 * (x1464 * x3 + 4.0 * x1465))
    result[6, 3] = numpy.sum(x758 * (x1466 * x3 + x1467))
    result[6, 4] = numpy.sum(x772 * (x107 * x1469 + x1468 * x3))
    result[6, 5] = numpy.sum(x758 * (x1470 * x3 + 3.0 * x1471))
    result[6, 6] = numpy.sum(x796 * (x1472 * x3 + x1473 * x180))
    result[6, 7] = numpy.sum(x813 * (x1474 * x3 + x1475 * x180))
    result[6, 8] = numpy.sum(x813 * (x1476 * x3 + x1477))
    result[6, 9] = numpy.sum(x796 * (x1478 * x3 + 2.0 * x1479))
    result[6, 10] = numpy.sum(x344 * (x1480 + x1481 * x3))
    result[6, 11] = numpy.sum(x772 * (x1482 + x1483 * x3))
    result[6, 12] = numpy.sum(x813 * (x1484 + x1485 * x3))
    result[6, 13] = numpy.sum(x772 * (x1486 + x1487 * x3))
    result[6, 14] = numpy.sum(x344 * (x1488 + x1489 * x3))
    result[6, 15] = numpy.sum(x1490 * x1491)
    result[6, 16] = numpy.sum(x1492 * x1493)
    result[6, 17] = numpy.sum(x1494 * x1495)
    result[6, 18] = numpy.sum(x1496 * x1497)
    result[6, 19] = numpy.sum(x1493 * x1498)
    result[6, 20] = numpy.sum(x1491 * x1499)
    result[7, 0] = numpy.sum(x1270 * (x1500 * x3 + 5.0 * x1501))
    result[7, 1] = numpy.sum(x1280 * (x1502 * x3 + x1503 * x98))
    result[7, 2] = numpy.sum(x1280 * (x1504 * x3 + 4.0 * x1505))
    result[7, 3] = numpy.sum(x813 * (x1506 * x3 + x1507))
    result[7, 4] = numpy.sum(x1291 * (x107 * x1509 + x1508 * x3))
    result[7, 5] = numpy.sum(x813 * (x1510 * x3 + 3.0 * x1511))
    result[7, 6] = numpy.sum(x1297 * (x1512 * x3 + x1513 * x180))
    result[7, 7] = numpy.sum(x1306 * (x1514 * x3 + x1515 * x180))
    result[7, 8] = numpy.sum(x1306 * (x1516 * x3 + x1517))
    result[7, 9] = numpy.sum(x1297 * (x1518 * x3 + 2.0 * x1519))
    result[7, 10] = numpy.sum(x1280 * (x1520 + x1521 * x3))
    result[7, 11] = numpy.sum(x1291 * (x1522 + x1523 * x3))
    result[7, 12] = numpy.sum(x1306 * (x1524 + x1525 * x3))
    result[7, 13] = numpy.sum(x1291 * (x1526 + x1527 * x3))
    result[7, 14] = numpy.sum(x1280 * (x1528 + x1529 * x3))
    result[7, 15] = numpy.sum(x1530 * x1531)
    result[7, 16] = numpy.sum(x1532 * x1533)
    result[7, 17] = numpy.sum(x1534 * x1535)
    result[7, 18] = numpy.sum(x1536 * x1537)
    result[7, 19] = numpy.sum(x1533 * x1538)
    result[7, 20] = numpy.sum(x1531 * x1539)
    result[8, 0] = numpy.sum(x1270 * x702 * (x1348 + x1351))
    result[8, 1] = numpy.sum(x1280 * (x1540 * x3 + x1541 * x98))
    result[8, 2] = numpy.sum(x1280 * x702 * (x1362 + x1365))
    result[8, 3] = numpy.sum(x813 * (x1542 * x3 + x1543))
    result[8, 4] = numpy.sum(x1291 * (x107 * x1545 + x1544 * x3))
    result[8, 5] = numpy.sum(x702 * x813 * (x1382 + x1385))
    result[8, 6] = numpy.sum(x1297 * (x1546 * x3 + x1547 * x180))
    result[8, 7] = numpy.sum(x1306 * (x1549 * x3 + x1550 * x180))
    result[8, 8] = numpy.sum(x1306 * (x1551 * x3 + x1552))
    result[8, 9] = numpy.sum(x1297 * x702 * (x1408 + x1411))
    result[8, 10] = numpy.sum(x1280 * (x1553 + x1554 * x3))
    result[8, 11] = numpy.sum(x1291 * (x1555 + x1556 * x3))
    result[8, 12] = numpy.sum(x1306 * (x1557 + x1558 * x3))
    result[8, 13] = numpy.sum(x1291 * (x1559 + x1560 * x3))
    result[8, 14] = numpy.sum(x1280 * (x1439 * x702 + x1561))
    result[8, 15] = numpy.sum(x1531 * x1562)
    result[8, 16] = numpy.sum(x1533 * x1563)
    result[8, 17] = numpy.sum(x1535 * x1565)
    result[8, 18] = numpy.sum(x1537 * x1566)
    result[8, 19] = numpy.sum(x1533 * x1567)
    result[8, 20] = numpy.sum(x1457 * x1531 * x702)
    result[9, 0] = numpy.sum(x712 * (x1568 * x3 + 5.0 * x1569))
    result[9, 1] = numpy.sum(x344 * (x1570 * x3 + 4.0 * x1571))
    result[9, 2] = numpy.sum(x344 * (x1572 * x3 + 4.0 * x1573))
    result[9, 3] = numpy.sum(x758 * (x1574 * x3 + x1576))
    result[9, 4] = numpy.sum(x772 * (x107 * x1578 + x1577 * x3))
    result[9, 5] = numpy.sum(x758 * (x1579 * x3 + x1581))
    result[9, 6] = numpy.sum(x796 * (x1582 * x3 + 2.0 * x1583))
    result[9, 7] = numpy.sum(x813 * (x1584 * x3 + x1586))
    result[9, 8] = numpy.sum(x813 * (x1587 * x3 + x1589))
    result[9, 9] = numpy.sum(x796 * (x1590 * x3 + 2.0 * x1591))
    result[9, 10] = numpy.sum(x344 * (x1592 + x1593 * x3))
    result[9, 11] = numpy.sum(x772 * (x1594 + x1595 * x3))
    result[9, 12] = numpy.sum(x813 * (x1596 + x1597 * x3))
    result[9, 13] = numpy.sum(x772 * (x1598 + x1599 * x3))
    result[9, 14] = numpy.sum(x344 * (x1600 + x1601 * x3))
    result[9, 15] = numpy.sum(x1491 * x1602)
    result[9, 16] = numpy.sum(x1493 * x1603)
    result[9, 17] = numpy.sum(x1495 * x1604)
    result[9, 18] = numpy.sum(x1497 * x1605)
    result[9, 19] = numpy.sum(x1493 * x1606)
    result[9, 20] = numpy.sum(x1491 * x1607)
    result[10, 0] = numpy.sum(x120 * (x118 * x1187 + x1460 * x702))
    result[10, 1] = numpy.sum(x187 * (x118 * x1194 + x1461 + x1462 * x702))
    result[10, 2] = numpy.sum(x187 * (x118 * x1199 + x1464 * x702))
    result[10, 3] = numpy.sum(x310 * (x118 * x1204 + x1463 * x180 + x1466 * x702))
    result[10, 4] = numpy.sum(x344 * (x118 * x1209 + x1465 + x1468 * x702))
    result[10, 5] = numpy.sum(x310 * (x118 * x1214 + x1470 * x702))
    result[10, 6] = numpy.sum(x448 * (x118 * x1218 + x1467 + x1472 * x702))
    result[10, 7] = numpy.sum(x484 * (x118 * x1223 + x1469 * x180 + x1474 * x702))
    result[10, 8] = numpy.sum(x484 * (x118 * x1229 + x1471 + x1476 * x702))
    result[10, 9] = numpy.sum(x448 * (x118 * x1234 + x1478 * x702))
    result[10, 10] = numpy.sum(x187 * (x118 * x1239 + x1473 * x98 + x1481 * x702))
    result[10, 11] = numpy.sum(x344 * (x107 * x1475 + x118 * x1243 + x1483 * x702))
    result[10, 12] = numpy.sum(x484 * (x118 * x1247 + x1477 + x1485 * x702))
    result[10, 13] = numpy.sum(x344 * (x118 * x1251 + x1479 + x1487 * x702))
    result[10, 14] = numpy.sum(x187 * (x118 * x1256 + x1489 * x702))
    result[10, 15] = numpy.sum(x120 * (x118 * x1258 + 5.0 * x1480 + x1490 * x702))
    result[10, 16] = numpy.sum(x187 * (x118 * x1260 + 4.0 * x1482 + x1492 * x702))
    result[10, 17] = numpy.sum(x310 * (x118 * x1262 + 3.0 * x1484 + x1494 * x702))
    result[10, 18] = numpy.sum(x448 * (x118 * x1264 + 2.0 * x1486 + x1496 * x702))
    result[10, 19] = numpy.sum(x187 * (x118 * x1266 + x1488 + x1498 * x702))
    result[10, 20] = numpy.sum(x120 * (x118 * x1268 + x1499 * x702))
    result[11, 0] = numpy.sum(x712 * (x1269 + x1500 * x702))
    result[11, 1] = numpy.sum(x344 * (x1277 + x1501 + x1502 * x702))
    result[11, 2] = numpy.sum(x344 * (x1281 + x1504 * x702))
    result[11, 3] = numpy.sum(x758 * (x1285 + x1503 * x180 + x1506 * x702))
    result[11, 4] = numpy.sum(x772 * (x1289 + x1505 + x1508 * x702))
    result[11, 5] = numpy.sum(x758 * (x1292 + x1510 * x702))
    result[11, 6] = numpy.sum(x796 * (x1295 + x1507 + x1512 * x702))
    result[11, 7] = numpy.sum(x813 * (x1303 + x1509 * x180 + x1514 * x702))
    result[11, 8] = numpy.sum(x813 * (x1310 + x1511 + x1516 * x702))
    result[11, 9] = numpy.sum(x796 * (x1312 + x1518 * x702))
    result[11, 10] = numpy.sum(x344 * (x1317 + x1513 * x98 + x1521 * x702))
    result[11, 11] = numpy.sum(x772 * (x107 * x1515 + x1321 + x1523 * x702))
    result[11, 12] = numpy.sum(x813 * (x1325 + x1517 + x1525 * x702))
    result[11, 13] = numpy.sum(x772 * (x1329 + x1519 + x1527 * x702))
    result[11, 14] = numpy.sum(x344 * (x1331 + x1529 * x702))
    result[11, 15] = numpy.sum(x712 * (x1333 + 5.0 * x1520 + x1530 * x702))
    result[11, 16] = numpy.sum(x344 * (x1335 + 4.0 * x1522 + x1532 * x702))
    result[11, 17] = numpy.sum(x758 * (x1340 + 3.0 * x1524 + x1534 * x702))
    result[11, 18] = numpy.sum(x796 * (x1342 + 2.0 * x1526 + x1536 * x702))
    result[11, 19] = numpy.sum(x344 * (x1344 + x1528 + x1538 * x702))
    result[11, 20] = numpy.sum(x712 * (x1345 + x1539 * x702))
    result[12, 0] = numpy.sum(x1189 * (x1181 * x1347 + x1353))
    result[12, 1] = numpy.sum(x484 * (x1350 * x702 + x1359 + x1540 * x702))
    result[12, 2] = numpy.sum(x484 * (x1181 * x1361 + x1367))
    result[12, 3] = numpy.sum(x1206 * (x1374 + x1541 * x180 + x1542 * x702))
    result[12, 4] = numpy.sum(x813 * (x1364 * x702 + x1379 + x1544 * x702))
    result[12, 5] = numpy.sum(x1206 * (x1181 * x1381 + x1387))
    result[12, 6] = numpy.sum(x1220 * (x1393 + x1543 + x1546 * x702))
    result[12, 7] = numpy.sum(x1225 * (x1399 + x1545 * x180 + x1549 * x702))
    result[12, 8] = numpy.sum(x1225 * (x1384 * x702 + x1405 + x1551 * x702))
    result[12, 9] = numpy.sum(x1220 * (x1181 * x1407 + x1413))
    result[12, 10] = numpy.sum(x484 * (x1420 + x1547 * x98 + x1554 * x702))
    result[12, 11] = numpy.sum(x813 * (x107 * x1550 + x1425 + x1556 * x702))
    result[12, 12] = numpy.sum(x1225 * (x1430 + x1552 + x1558 * x702))
    result[12, 13] = numpy.sum(x813 * (x1410 * x702 + x1435 + x1560 * x702))
    result[12, 14] = numpy.sum(x484 * (x1181 * x1438 + x1441))
    result[12, 15] = numpy.sum(x1189 * (x1444 + 5.0 * x1553 + x1562 * x702))
    result[12, 16] = numpy.sum(x484 * (x1447 + 4.0 * x1555 + x1563 * x702))
    result[12, 17] = numpy.sum(x1206 * (x1450 + 3.0 * x1557 + x1565 * x702))
    result[12, 18] = numpy.sum(x1220 * (x1453 + 2.0 * x1559 + x1566 * x702))
    result[12, 19] = numpy.sum(x484 * (x1456 + x1561 + x1567 * x702))
    result[12, 20] = numpy.sum(x1189 * (x1181 * x1457 + x1459))
    result[13, 0] = numpy.sum(x1568 * x1608)
    result[13, 1] = numpy.sum(x344 * (x1569 + x1570 * x702))
    result[13, 2] = numpy.sum(x1572 * x1609)
    result[13, 3] = numpy.sum(x758 * (2.0 * x1571 + x1574 * x702))
    result[13, 4] = numpy.sum(x772 * (x1573 + x1577 * x702))
    result[13, 5] = numpy.sum(x1579 * x702 * x758)
    result[13, 6] = numpy.sum(x796 * (x1576 + x1582 * x702))
    result[13, 7] = numpy.sum(x813 * (x1584 * x702 + x1610))
    result[13, 8] = numpy.sum(x813 * (x1580 + x1587 * x702))
    result[13, 9] = numpy.sum(x1590 * x702 * x796)
    result[13, 10] = numpy.sum(x344 * (4.0 * x1583 + x1593 * x702))
    result[13, 11] = numpy.sum(x772 * (x107 * x1585 + x1595 * x702))
    result[13, 12] = numpy.sum(x813 * (x1589 + x1597 * x702))
    result[13, 13] = numpy.sum(x772 * (x1591 + x1599 * x702))
    result[13, 14] = numpy.sum(x1601 * x1609)
    result[13, 15] = numpy.sum(x712 * (5.0 * x1592 + x1602 * x702))
    result[13, 16] = numpy.sum(x344 * (4.0 * x1594 + x1603 * x702))
    result[13, 17] = numpy.sum(x758 * (x1604 * x702 + x1611))
    result[13, 18] = numpy.sum(x796 * (2.0 * x1598 + x1605 * x702))
    result[13, 19] = numpy.sum(x344 * (x1600 + x1606 * x702))
    result[13, 20] = numpy.sum(x1607 * x1608)
    result[14, 0] = numpy.sum(x120 * (x118 * x1352 + x1568 * x930))
    result[14, 1] = numpy.sum(x187 * (x118 * x1358 + x1570 * x930))
    result[14, 2] = numpy.sum(x187 * (x118 * x1366 + x1569 + x1572 * x930))
    result[14, 3] = numpy.sum(x310 * (x118 * x1373 + x1574 * x930))
    result[14, 4] = numpy.sum(x344 * (x118 * x1378 + x1571 + x1577 * x930))
    result[14, 5] = numpy.sum(x310 * (x118 * x1386 + 2.0 * x1573 + x1579 * x930))
    result[14, 6] = numpy.sum(x448 * (x118 * x1392 + x1582 * x930))
    result[14, 7] = numpy.sum(x484 * (x118 * x1398 + x1575 + x1584 * x930))
    result[14, 8] = numpy.sum(x484 * (x118 * x1404 + x1587 * x930 + x1610))
    result[14, 9] = numpy.sum(x448 * (x118 * x1412 + x1581 + x1590 * x930))
    result[14, 10] = numpy.sum(x187 * (x118 * x1419 + x1593 * x930))
    result[14, 11] = numpy.sum(x344 * (x118 * x1424 + x1583 + x1595 * x930))
    result[14, 12] = numpy.sum(x484 * (x118 * x1429 + x1586 + x1597 * x930))
    result[14, 13] = numpy.sum(x344 * (x107 * x1588 + x118 * x1434 + x1599 * x930))
    result[14, 14] = numpy.sum(x187 * (x118 * x1440 + 4.0 * x1591 + x1601 * x930))
    result[14, 15] = numpy.sum(x120 * (x118 * x1443 + x1602 * x930))
    result[14, 16] = numpy.sum(x187 * (x118 * x1446 + x1592 + x1603 * x930))
    result[14, 17] = numpy.sum(x310 * (x118 * x1449 + 2.0 * x1594 + x1604 * x930))
    result[14, 18] = numpy.sum(x448 * (x118 * x1452 + x1605 * x930 + x1611))
    result[14, 19] = numpy.sum(x187 * (x118 * x1455 + 4.0 * x1598 + x1606 * x930))
    result[14, 20] = numpy.sum(x120 * (x118 * x1458 + 5.0 * x1600 + x1607 * x930))
    return result


def _2center2el3d_50(ax, da, A, bx, db, B):
    """Cartesian (h|s) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((21, 1), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = x1 * (ax * A[0] + bx * B[0]) - A[0]
    x3 = ax ** (-1.0)
    x4 = bx * x1
    x5 = ax * x4 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x6 = boys(4, x5)
    x7 = 17.4934183276249
    x8 = 2.0 * x3
    x9 = x7 * x8
    x10 = x0 ** (-1.5) * x9
    x11 = x10 * x6
    x12 = x11 * x2
    x13 = bx ** (-1.0)
    x14 = x0 ** (-0.5)
    x15 = boys(3, x5)
    x16 = 0.5 * x3
    x17 = x16 * (-x11 + 2.0 * x13 * x14 * x15 * x3 * x7)
    x18 = boys(5, x5)
    x19 = x2**2
    x20 = x13 * x14 * x9
    x21 = x19 * x20
    x22 = x18 * x21
    x23 = x10 * x15
    x24 = boys(2, x5)
    x25 = x16 * (2.0 * x13 * x14 * x24 * x3 * x7 - x23)
    x26 = x21 * x6
    x27 = x25 + x26
    x28 = x10 * x24
    x29 = boys(1, x5)
    x30 = x16 * (2.0 * x13 * x14 * x29 * x3 * x7 - x28)
    x31 = x15 * x20
    x32 = x19 * x31
    x33 = x30 + x32
    x34 = 1.5 * x3
    x35 = 0.179587122125167 * da * db * numpy.sqrt(ax**6.5) * numpy.sqrt(bx**1.5)
    x36 = 2.94427971986071 * x35
    x37 = x1 * (ax * A[1] + bx * B[1]) - A[1]
    x38 = x12 * x37
    x39 = x11 * x37
    x40 = x3 * (2.0 * x13 * x14 * x15 * x3 * x37 * x7 - x39)
    x41 = x22 * x37
    x42 = x3 * x37 * (2.0 * x13 * x14 * x29 * x3 * x7 - x28)
    x43 = x3 * x37 * (2.0 * x13 * x14 * x24 * x3 * x7 - x23)
    x44 = 8.83283915958214 * x35
    x45 = x1 * (ax * A[2] + bx * B[2]) - A[2]
    x46 = x3 * x45 * (-x11 + 2.0 * x13 * x14 * x15 * x3 * x7)
    x47 = 0.5 * x46
    x48 = x3 * x45 * (2.0 * x13 * x14 * x29 * x3 * x7 - x28)
    x49 = 0.5 * x48
    x50 = x3 * x45 * (2.0 * x13 * x14 * x24 * x3 * x7 - x23)
    x51 = 0.5 * x50
    x52 = x37**2
    x53 = x30 + x31 * x52
    x54 = x20 * x52
    x55 = x54 * x6
    x56 = x25 + x55
    x57 = x4 * x56
    x58 = x18 * x54
    x59 = x17 + x58
    x60 = x53 - x57
    x61 = 13.4923846833851 * x35
    x62 = x3 * x45 * (2.0 * x13 * x14 * x15 * x3 * x37 * x7 - x39)
    x63 = 23.3694957868871 * x35
    x64 = x45**2
    x65 = x30 + x31 * x64
    x66 = x20 * x64
    x67 = x25 + x6 * x66
    x68 = x4 * x67
    x69 = x17 + x18 * x66
    x70 = x19 * x69
    x71 = x65 - x68
    x72 = x16 * x71
    x73 = x37 * x59 + x40
    x74 = x37 * x53 - x4 * (x37 * x56 + x43) + x42
    x75 = 13.4923846833851 * x35
    x76 = x45 * x58 + x47
    x77 = x31 * x45 * x52 - x4 * (x45 * x55 + x51) + x49
    x78 = 30.169889330626 * x35
    x79 = x3 * x37 * (x65 - x68)
    x80 = x45 * x69 + x46
    x81 = -x4 * (x45 * x67 + x50) + x45 * x65 + x48
    x82 = x16 * x81
    x83 = x34 * x60 + x37 * x73
    x84 = x2 * x44
    x85 = x37 * x76 + x62
    x86 = x2 * x63
    x87 = x52 * x69 + x72
    x88 = x34 * x71 + x45 * x80

    # 21 item(s)
    result[0, 0] = numpy.sum(
        x2
        * x36
        * (
            x2 * (x2 * (x17 + x22) - x3 * (x12 - 2.0 * x13 * x14 * x15 * x2 * x3 * x7))
            - x34 * (x27 * x4 - x33)
            + x8
            * (
                x3 * (2.0 * x13 * x14 * x29 * x3 * x7 - x28)
                + x33
                - x4 * (x27 + x3 * (2.0 * x13 * x14 * x24 * x3 * x7 - x23))
            )
        )
    )
    result[1, 0] = numpy.sum(
        0.5
        * x44
        * (
            x2
            * (
                x2 * (x40 + 2.0 * x41)
                + 2.0 * x3 * (2.0 * x13 * x14 * x15 * x2 * x3 * x37 * x7 - x38)
            )
            + x34 * (2.0 * x32 * x37 - x4 * (2.0 * x26 * x37 + x43) + x42)
        )
    )
    result[2, 0] = numpy.sum(
        x44
        * (
            x2
            * (
                x2 * (x22 * x45 + x47)
                - x3 * x45 * (x12 - 2.0 * x13 * x14 * x15 * x2 * x3 * x7)
            )
            + x34 * (x32 * x45 - x4 * (x26 * x45 + x51) + x49)
        )
    )
    result[3, 0] = numpy.sum(x2 * x61 * (x16 * x60 + x19 * x59 + x3 * (x53 - x57)))
    result[4, 0] = numpy.sum(
        0.5
        * x63
        * (
            x2 * (2.0 * x41 * x45 + x62)
            + 2.0 * x3 * x45 * (2.0 * x13 * x14 * x15 * x2 * x3 * x37 * x7 - x38)
        )
    )
    result[5, 0] = numpy.sum(x2 * x61 * (x3 * (x65 - x68) + x70 + x72))
    result[6, 0] = numpy.sum(x75 * (x16 * x74 + x19 * x73))
    result[7, 0] = numpy.sum(x78 * (x16 * x77 + x19 * x76))
    result[8, 0] = numpy.sum(0.5 * x78 * (2.0 * x37 * x70 + x79))
    result[9, 0] = numpy.sum(x75 * (x19 * x80 + x82))
    result[10, 0] = numpy.sum(x83 * x84)
    result[11, 0] = numpy.sum(x85 * x86)
    result[12, 0] = numpy.sum(x2 * x78 * x87)
    result[13, 0] = numpy.sum(x37 * x80 * x86)
    result[14, 0] = numpy.sum(x84 * x88)
    result[15, 0] = numpy.sum(x36 * (x37 * x83 + x74 * x8))
    result[16, 0] = numpy.sum(x44 * (x34 * x77 + x37 * x85))
    result[17, 0] = numpy.sum(x61 * (x37 * x87 + x79))
    result[18, 0] = numpy.sum(x75 * (x52 * x80 + x82))
    result[19, 0] = numpy.sum(x37 * x44 * x88)
    result[20, 0] = numpy.sum(x36 * (x45 * x88 + x8 * x81))
    return result


def _2center2el3d_51(ax, da, A, bx, db, B):
    """Cartesian (h|p) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((21, 3), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = -x2 * (ax * A[0] + bx * B[0])
    x4 = -x3 - A[0]
    x5 = ax ** (-1.0)
    x6 = bx * x2
    x7 = ax * x6 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x8 = boys(4, x7)
    x9 = 17.4934183276249
    x10 = 2.0 * x5
    x11 = x10 * x9
    x12 = x1 ** (-1.5) * x11
    x13 = x12 * x8
    x14 = x13 * x4
    x15 = bx ** (-1.0)
    x16 = x1 ** (-0.5)
    x17 = boys(3, x7)
    x18 = 0.5 * x5
    x19 = x18 * (-x13 + 2.0 * x15 * x16 * x17 * x5 * x9)
    x20 = boys(5, x7)
    x21 = x4**2
    x22 = x11 * x15 * x16
    x23 = x21 * x22
    x24 = x20 * x23
    x25 = x19 + x24
    x26 = x25 * x4 - x5 * (x14 - 2.0 * x15 * x16 * x17 * x4 * x5 * x9)
    x27 = x12 * x17
    x28 = boys(2, x7)
    x29 = x18 * (2.0 * x15 * x16 * x28 * x5 * x9 - x27)
    x30 = x23 * x8
    x31 = x29 + x30
    x32 = boys(1, x7)
    x33 = x18 * (-x12 * x28 + 2.0 * x15 * x16 * x32 * x5 * x9)
    x34 = x17 * x22
    x35 = x21 * x34 + x33
    x36 = 1.5 * x5
    x37 = x0 * x22
    x38 = x37 * x8
    x39 = x20 * x22
    x40 = -x3 - B[0]
    x41 = x4 * x40
    x42 = x39 * x41
    x43 = x38 + x42
    x44 = x0 * x34
    x45 = x22 * x8
    x46 = x40 * x45
    x47 = x4 * x46
    x48 = x44 + x47
    x49 = x20 * x37
    x50 = boys(6, x7)
    x51 = x22 * x50
    x52 = x41 * x51
    x53 = x12 * x20
    x54 = x40 * x53
    x55 = x18 * (2.0 * x15 * x16 * x40 * x5 * x8 * x9 - x54)
    x56 = x4 * x49
    x57 = x13 * x40
    x58 = x18 * (2.0 * x15 * x16 * x17 * x40 * x5 * x9 - x57)
    x59 = x38 * x4
    x60 = x4 * x43 + x58 + x59
    x61 = x27 * x40
    x62 = x18 * (2.0 * x15 * x16 * x28 * x40 * x5 * x9 - x61)
    x63 = x34 * x4
    x64 = x0 * x63
    x65 = x4 * x48 + x62 + x64
    x66 = x22 * x28
    x67 = x0 * x66
    x68 = x40 * x63 + x67
    x69 = x0 * x22 * x32
    x70 = 0.179587122125167 * da * db * numpy.sqrt(ax**6.5) * numpy.sqrt(bx**2.5)
    x71 = 5.88855943972142 * x70
    x72 = -x2 * (ax * A[1] + bx * B[1])
    x73 = -x72 - B[1]
    x74 = x53 * x73
    x75 = x4 * x74
    x76 = x18 * (2.0 * x15 * x16 * x5 * x73 * x8 * x9 - x74)
    x77 = x23 * x50
    x78 = x73 * x77
    x79 = x13 * x73
    x80 = x18 * (2.0 * x15 * x16 * x17 * x5 * x73 * x9 - x79)
    x81 = x24 * x73
    x82 = x80 + x81
    x83 = x27 * x73
    x84 = x18 * (2.0 * x15 * x16 * x28 * x5 * x73 * x9 - x83)
    x85 = x30 * x73
    x86 = x84 + x85
    x87 = -x2 * (ax * A[2] + bx * B[2])
    x88 = -x87 - B[2]
    x89 = x53 * x88
    x90 = x4 * x89
    x91 = x18 * (2.0 * x15 * x16 * x5 * x8 * x88 * x9 - x89)
    x92 = x77 * x88
    x93 = x13 * x88
    x94 = x18 * (2.0 * x15 * x16 * x17 * x5 * x88 * x9 - x93)
    x95 = x24 * x88
    x96 = x94 + x95
    x97 = x27 * x88
    x98 = x18 * (2.0 * x15 * x16 * x28 * x5 * x88 * x9 - x97)
    x99 = x30 * x88
    x100 = x98 + x99
    x101 = -x72 - A[1]
    x102 = x101 * x13
    x103 = x5 * (2.0 * x101 * x15 * x16 * x17 * x5 * x9 - x102)
    x104 = x101 * x24
    x105 = 0.5 * x103 + x104
    x106 = x101 * x38
    x107 = x101 * x42
    x108 = x106 + x107
    x109 = x101 * x44
    x110 = x101 * x47
    x111 = x109 + x110
    x112 = x101 * x49
    x113 = x101 * x52
    x114 = x101 * x54
    x115 = x5 * (2.0 * x101 * x15 * x16 * x40 * x5 * x8 * x9 - x114)
    x116 = x101 * x56
    x117 = x101 * x5 * (2.0 * x15 * x16 * x28 * x40 * x5 * x9 - x61)
    x118 = x101 * x5 * (2.0 * x15 * x16 * x17 * x40 * x5 * x9 - x57)
    x119 = 17.6656783191643 * x70
    x120 = x45 * x73
    x121 = x101 * x120
    x122 = x121 + x44
    x123 = x101 * x73
    x124 = x123 * x39
    x125 = x124 + x38
    x126 = x125 * x6
    x127 = x123 * x51
    x128 = x127 + x49
    x129 = x5 * (x122 - x126)
    x130 = x101 * x34 * x73 + x67
    x131 = x5 * (x101 * x66 * x73 - x130 * x6 + x69)
    x132 = x5 * (-x122 * x6 + x130)
    x133 = x101 * x5 * (2.0 * x15 * x16 * x5 * x8 * x88 * x9 - x89)
    x134 = x101 * x5 * (2.0 * x15 * x16 * x28 * x5 * x88 * x9 - x97)
    x135 = x101 * x5 * (2.0 * x15 * x16 * x17 * x5 * x88 * x9 - x93)
    x136 = -x87 - A[2]
    x137 = x136 * x5 * (-x13 + 2.0 * x15 * x16 * x17 * x5 * x9)
    x138 = 0.5 * x137
    x139 = x136 * x24 + x138
    x140 = x136 * x38
    x141 = x136 * x42 + x140
    x142 = x136 * x44
    x143 = x136 * x47 + x142
    x144 = x136 * x49
    x145 = x136 * x5 * (2.0 * x15 * x16 * x40 * x5 * x8 * x9 - x54)
    x146 = 0.5 * x145
    x147 = x136 * x5 * (2.0 * x15 * x16 * x28 * x40 * x5 * x9 - x61)
    x148 = 0.5 * x147
    x149 = x136 * x5 * (2.0 * x15 * x16 * x17 * x40 * x5 * x9 - x57)
    x150 = 0.5 * x149
    x151 = x136 * x5 * (2.0 * x15 * x16 * x5 * x73 * x8 * x9 - x74)
    x152 = 0.5 * x151
    x153 = x136 * x5 * (2.0 * x15 * x16 * x28 * x5 * x73 * x9 - x83)
    x154 = 0.5 * x153
    x155 = x136 * x5 * (2.0 * x15 * x16 * x17 * x5 * x73 * x9 - x79)
    x156 = 0.5 * x155
    x157 = x45 * x88
    x158 = x136 * x157 + x44
    x159 = x136 * x88
    x160 = x159 * x39 + x38
    x161 = x160 * x6
    x162 = x161 * x4
    x163 = x159 * x51 + x49
    x164 = x163 * x21
    x165 = x5 * (x158 - x161)
    x166 = 0.5 * x165
    x167 = x136 * x34 * x88 + x67
    x168 = x5 * (x136 * x66 * x88 - x167 * x6 + x69)
    x169 = 0.5 * x168
    x170 = x5 * (-x158 * x6 + x167)
    x171 = 0.5 * x170
    x172 = x101**2
    x173 = x172 * x34 + x33
    x174 = x0 * x173
    x175 = x172 * x46 + x62
    x176 = x172 * x45 + x29
    x177 = x0 * x176
    x178 = x172 * x39
    x179 = x178 * x40
    x180 = x179 + x58
    x181 = x178 + x19
    x182 = x173 - x176 * x6
    x183 = x0 * x181
    x184 = x172 * x51
    x185 = x184 * x40
    x186 = x185 + x55
    x187 = x175 - x180 * x6
    x188 = 26.9847693667702 * x70
    x189 = x101 * x122 + x109 + x84
    x190 = x101 * x125 + x106 + x80
    x191 = x190 * x6
    x192 = x101 * x128 + x112 + x76
    x193 = x189 - x191
    x194 = x157 * x172 + x98
    x195 = x178 * x88 + x94
    x196 = x195 * x6
    x197 = x184 * x88 + x91
    x198 = x194 - x196
    x199 = x106 * x136
    x200 = x109 * x136
    x201 = x136 * x5 * (2.0 * x101 * x15 * x16 * x17 * x5 * x9 - x102)
    x202 = x112 * x136
    x203 = x136 * x5 * (2.0 * x101 * x15 * x16 * x40 * x5 * x8 * x9 - x114)
    x204 = 46.7389915737742 * x70
    x205 = x121 * x136 + x142
    x206 = x124 * x136 + x140
    x207 = x206 * x6
    x208 = x127 * x136 + x144
    x209 = x5 * (x205 - x207)
    x210 = x101 * x5 * (x158 - x161)
    x211 = x136**2
    x212 = x211 * x34 + x33
    x213 = x0 * x212
    x214 = x211 * x46 + x62
    x215 = x211 * x45 + x29
    x216 = x0 * x215
    x217 = x211 * x39
    x218 = x217 * x40 + x58
    x219 = x19 + x217
    x220 = x212 - x215 * x6
    x221 = x18 * x220
    x222 = x0 * x219
    x223 = x211 * x51
    x224 = x223 * x40 + x55
    x225 = x224 * x4
    x226 = x218 * x6
    x227 = x214 - x226
    x228 = x18 * x227
    x229 = x222 * x4
    x230 = x120 * x211 + x84
    x231 = x217 * x73 + x80
    x232 = x231 * x6
    x233 = x223 * x73 + x76
    x234 = x230 - x232
    x235 = x18 * x234
    x236 = x136 * x158 + x142 + x98
    x237 = x136 * x160 + x140 + x94
    x238 = x237 * x6
    x239 = x136 * x163 + x144 + x91
    x240 = x21 * x239
    x241 = x236 - x238
    x242 = x18 * x241
    x243 = x101 * x181 + x103
    x244 = x0 * x243
    x245 = x101 * x186 + x115
    x246 = x101 * x175 + x117 - x6 * (x101 * x180 + x118)
    x247 = 26.9847693667702 * x70
    x248 = x101 * x192 + x129 + x183
    x249 = x101 * x189 + x131 + x174 - x6 * (x101 * x190 + x132 + x177)
    x250 = x101 * x197 + x133
    x251 = x101 * x194 + x134 - x6 * (x101 * x195 + x135)
    x252 = x136 * x178 + x138
    x253 = x0 * x252
    x254 = x136 * x185 + x146
    x255 = x136 * x172 * x46 + x148 - x6 * (x136 * x179 + x150)
    x256 = 60.3397786612521 * x70
    x257 = x101 * x208 + x152 + x202
    x258 = x101 * x205 + x154 + x200 - x6 * (x101 * x206 + x156 + x199)
    x259 = x163 * x172 + x166
    x260 = x158 * x172 + x169 - x6 * (x160 * x172 + x171)
    x261 = x101 * x222
    x262 = x101 * x5 * (x214 - x226)
    x263 = x101 * x233 + x222
    x264 = x5 * (x101 * x230 + x213 - x6 * (x101 * x231 + x216))
    x265 = x101 * x5 * (x236 - x238)
    x266 = x136 * x219 + x137
    x267 = x0 * x266
    x268 = x136 * x224 + x145
    x269 = x268 * x4
    x270 = x136 * x214 + x147 - x6 * (x136 * x218 + x149)
    x271 = x18 * x270
    x272 = x136 * x233 + x151
    x273 = x136 * x230 + x153 - x6 * (x136 * x231 + x155)
    x274 = x18 * x273
    x275 = x136 * x239 + x165 + x222
    x276 = x136 * x236 + x168 + x213 - x6 * (x136 * x237 + x170 + x216)
    x277 = x18 * x276
    x278 = x0 * (x101 * x243 + x182 * x36)
    x279 = x101 * x245 + x187 * x36
    x280 = x101 * x248 + x193 * x36 + x244
    x281 = x119 * x4
    x282 = x101 * x250 + x198 * x36
    x283 = x0 * (x101 * x252 + x201)
    x284 = x101 * x254 + x203
    x285 = x101 * x257 + x209 + x253
    x286 = x204 * x4
    x287 = x101 * x259 + x210
    x288 = x0 * (x172 * x219 + x221)
    x289 = x172 * x224 + x228
    x290 = x101 * x263 + x235 + x261
    x291 = x256 * x4
    x292 = x172 * x239 + x242
    x293 = x101 * x267
    x294 = x101 * x272 + x267
    x295 = x0 * (x136 * x266 + x220 * x36)
    x296 = x136 * x268 + x227 * x36
    x297 = x136 * x272 + x234 * x36
    x298 = x136 * x275 + x241 * x36 + x267
    x299 = x101 * x119

    # 63 item(s)
    result[0, 0] = numpy.sum(
        x71
        * (
            x0 * (x26 * x4 - x36 * (x31 * x6 - x35))
            + x10
            * (
                x0 * x35
                + x4 * x65
                + x5 * (x4 * x40 * x66 - x6 * x68 + x69)
                - x6 * (x0 * x31 + x4 * x60 - x5 * (x48 * x6 - x68))
            )
            + x4
            * (
                x0 * x26
                - x36 * (x6 * x60 - x65)
                + x4
                * (x0 * x25 + x4 * (x4 * (x49 + x52) + x55 + x56) - x5 * (x43 * x6 - x48))
            )
        )
    )
    result[0, 1] = numpy.sum(
        x71
        * (
            x10
            * (
                x4 * x5 * (2.0 * x15 * x16 * x28 * x5 * x73 * x9 - x83)
                + x4 * x86
                - x6
                * (x4 * x82 - x5 * x73 * (x14 - 2.0 * x15 * x16 * x17 * x4 * x5 * x9))
            )
            - x4
            * (
                x36 * (x6 * x82 - x86)
                - x4
                * (
                    x4 * (x76 + x78)
                    + x5 * (2.0 * x15 * x16 * x4 * x5 * x73 * x8 * x9 - x75)
                )
            )
        )
    )
    result[0, 2] = numpy.sum(
        x71
        * (
            x10
            * (
                x100 * x4
                + x4 * x5 * (2.0 * x15 * x16 * x28 * x5 * x88 * x9 - x97)
                - x6
                * (x4 * x96 - x5 * x88 * (x14 - 2.0 * x15 * x16 * x17 * x4 * x5 * x9))
            )
            + x4
            * (
                x36 * (x100 - x6 * x96)
                + x4
                * (
                    x4 * (x91 + x92)
                    + x5 * (2.0 * x15 * x16 * x4 * x5 * x8 * x88 * x9 - x90)
                )
            )
        )
    )
    result[1, 0] = numpy.sum(
        0.5
        * x119
        * (
            2.0
            * x0
            * (-x101 * x5 * (x14 - 2.0 * x15 * x16 * x17 * x4 * x5 * x9) + x105 * x4)
            + x36
            * (
                2.0 * x101 * x64
                + 2.0 * x111 * x4
                + x117
                - x6 * (2.0 * x101 * x59 + 2.0 * x108 * x4 + x118)
            )
            + x4
            * (
                2.0 * x0 * x105
                + x4 * (x115 + 2.0 * x116 + 2.0 * x4 * (x112 + x113))
                - 2.0 * x5 * (x108 * x6 - x111)
            )
        )
    )
    result[1, 1] = numpy.sum(
        0.5
        * x119
        * (
            x36 * (2.0 * x122 * x21 + x131 - x6 * (2.0 * x125 * x21 + x132))
            + x4**2 * (2.0 * x128 * x21 + x129 + 2.0 * x5 * (x122 - x126))
        )
    )
    result[1, 2] = numpy.sum(
        0.5
        * x119
        * (
            x36 * (2.0 * x101 * x99 + x134 - x6 * (2.0 * x101 * x95 + x135))
            + x4
            * (
                2.0 * x101 * x5 * (2.0 * x15 * x16 * x4 * x5 * x8 * x88 * x9 - x90)
                + x4 * (2.0 * x101 * x92 + x133)
            )
        )
    )
    result[2, 0] = numpy.sum(
        x119
        * (
            x0 * (-x136 * x5 * (x14 - 2.0 * x15 * x16 * x17 * x4 * x5 * x9) + x139 * x4)
            + x36 * (x136 * x64 + x143 * x4 + x148 - x6 * (x136 * x59 + x141 * x4 + x150))
            + x4
            * (
                x0 * x139
                + x4 * (x136 * x56 + x146 + x4 * (x136 * x52 + x144))
                - x5 * (x141 * x6 - x143)
            )
        )
    )
    result[2, 1] = numpy.sum(
        x119
        * (
            x36 * (x136 * x85 + x154 - x6 * (x136 * x81 + x156))
            + x4
            * (
                x136 * x5 * (2.0 * x15 * x16 * x4 * x5 * x73 * x8 * x9 - x75)
                + x4 * (x136 * x78 + x152)
            )
        )
    )
    result[2, 2] = numpy.sum(
        x119
        * (
            x36 * (x158 * x21 + x169 - x6 * (x160 * x21 + x171))
            + x4 * (x4 * (x164 + x166) + x5 * (x158 * x4 - x162))
        )
    )
    result[3, 0] = numpy.sum(
        x188
        * (
            x0 * (x18 * x182 + x181 * x21)
            + x4 * (x18 * x187 + x183 * x4 + x4 * (x183 + x186 * x4))
            + x5 * (x174 + x175 * x4 - x6 * (x177 + x180 * x4))
        )
    )
    result[3, 1] = numpy.sum(x188 * x4 * (x18 * x193 + x192 * x21 + x5 * (x189 - x191)))
    result[3, 2] = numpy.sum(x188 * x4 * (x18 * x198 + x197 * x21 + x5 * (x194 - x196)))
    result[4, 0] = numpy.sum(
        0.5
        * x204
        * (
            x0 * (2.0 * x104 * x136 + x201)
            + x4 * (2.0 * x116 * x136 + x203 + 2.0 * x4 * (x113 * x136 + x202))
            + 2.0 * x5 * (x110 * x136 + x200 - x6 * (x107 * x136 + x199))
        )
    )
    result[4, 1] = numpy.sum(
        0.5 * x204 * x4 * (2.0 * x208 * x21 + x209 + 2.0 * x5 * (x205 - x207))
    )
    result[4, 2] = numpy.sum(
        0.5
        * x204
        * (2.0 * x101 * x5 * (x158 * x4 - x162) + x4 * (2.0 * x101 * x164 + x210))
    )
    result[5, 0] = numpy.sum(
        x188
        * (
            x0 * (x21 * x219 + x221)
            + x4 * (x228 + x229 + x4 * (x222 + x225))
            + x5 * (x213 + x214 * x4 - x6 * (x216 + x218 * x4))
        )
    )
    result[5, 1] = numpy.sum(x188 * x4 * (x21 * x233 + x235 + x5 * (x230 - x232)))
    result[5, 2] = numpy.sum(x188 * x4 * (x240 + x242 + x5 * (x236 - x238)))
    result[6, 0] = numpy.sum(x247 * (x18 * x246 + x244 * x4 + x4 * (x244 + x245 * x4)))
    result[6, 1] = numpy.sum(x247 * (x18 * x249 + x21 * x248))
    result[6, 2] = numpy.sum(x247 * (x18 * x251 + x21 * x250))
    result[7, 0] = numpy.sum(x256 * (x18 * x255 + x253 * x4 + x4 * (x253 + x254 * x4)))
    result[7, 1] = numpy.sum(x256 * (x18 * x258 + x21 * x257))
    result[7, 2] = numpy.sum(x256 * (x18 * x260 + x21 * x259))
    result[8, 0] = numpy.sum(
        0.5 * x256 * (2.0 * x101 * x229 + x262 + 2.0 * x4 * (x101 * x225 + x261))
    )
    result[8, 1] = numpy.sum(0.5 * x256 * (2.0 * x21 * x263 + x264))
    result[8, 2] = numpy.sum(0.5 * x256 * (2.0 * x101 * x240 + x265))
    result[9, 0] = numpy.sum(x247 * (x267 * x4 + x271 + x4 * (x267 + x269)))
    result[9, 1] = numpy.sum(x247 * (x21 * x272 + x274))
    result[9, 2] = numpy.sum(x247 * (x21 * x275 + x277))
    result[10, 0] = numpy.sum(x119 * (x278 + x279 * x4))
    result[10, 1] = numpy.sum(x280 * x281)
    result[10, 2] = numpy.sum(x281 * x282)
    result[11, 0] = numpy.sum(x204 * (x283 + x284 * x4))
    result[11, 1] = numpy.sum(x285 * x286)
    result[11, 2] = numpy.sum(x286 * x287)
    result[12, 0] = numpy.sum(x256 * (x288 + x289 * x4))
    result[12, 1] = numpy.sum(x290 * x291)
    result[12, 2] = numpy.sum(x291 * x292)
    result[13, 0] = numpy.sum(x204 * (x101 * x269 + x293))
    result[13, 1] = numpy.sum(x286 * x294)
    result[13, 2] = numpy.sum(x101 * x275 * x286)
    result[14, 0] = numpy.sum(x119 * (x295 + x296 * x4))
    result[14, 1] = numpy.sum(x281 * x297)
    result[14, 2] = numpy.sum(x281 * x298)
    result[15, 0] = numpy.sum(x71 * (x10 * x246 + x101 * x279))
    result[15, 1] = numpy.sum(x71 * (x10 * x249 + x101 * x280 + x278))
    result[15, 2] = numpy.sum(x71 * (x10 * x251 + x101 * x282))
    result[16, 0] = numpy.sum(x119 * (x101 * x284 + x255 * x36))
    result[16, 1] = numpy.sum(x119 * (x101 * x285 + x258 * x36 + x283))
    result[16, 2] = numpy.sum(x119 * (x101 * x287 + x260 * x36))
    result[17, 0] = numpy.sum(x188 * (x101 * x289 + x262))
    result[17, 1] = numpy.sum(x188 * (x101 * x290 + x264 + x288))
    result[17, 2] = numpy.sum(x188 * (x101 * x292 + x265))
    result[18, 0] = numpy.sum(x247 * (x172 * x268 + x271))
    result[18, 1] = numpy.sum(x247 * (x101 * x294 + x274 + x293))
    result[18, 2] = numpy.sum(x247 * (x172 * x275 + x277))
    result[19, 0] = numpy.sum(x296 * x299)
    result[19, 1] = numpy.sum(x119 * (x101 * x297 + x295))
    result[19, 2] = numpy.sum(x298 * x299)
    result[20, 0] = numpy.sum(x71 * (x10 * x270 + x136 * x296))
    result[20, 1] = numpy.sum(x71 * (x10 * x273 + x136 * x297))
    result[20, 2] = numpy.sum(x71 * (x10 * x276 + x136 * x298 + x295))
    return result


def _2center2el3d_52(ax, da, A, bx, db, B):
    """Cartesian (h|d) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((21, 6), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = ax ** (-1.0)
    x5 = bx * x1
    x6 = ax * x5 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x7 = boys(5, x6)
    x8 = x0 ** (-1.5)
    x9 = bx ** (-1.0)
    x10 = 17.4934183276249
    x11 = x10 * x9
    x12 = 2.0 * x11
    x13 = x12 * x8
    x14 = x0 ** (-0.5)
    x15 = boys(4, x6)
    x16 = 0.5 * x9
    x17 = x16 * (2.0 * x10 * x14 * x15 * x4 * x9 - x13 * x7)
    x18 = boys(6, x6)
    x19 = -x2 - B[0]
    x20 = x19**2
    x21 = 2.0 * x4
    x22 = x11 * x14
    x23 = x21 * x22
    x24 = x20 * x23
    x25 = x17 + x18 * x24
    x26 = x25 * x3
    x27 = x19 * x7
    x28 = 0.5 / (ax + bx)
    x29 = x22 * x4
    x30 = 4.0 * x28 * x29
    x31 = x27 * x30
    x32 = x26 + x31
    x33 = x15 * x8
    x34 = boys(3, x6)
    x35 = -2.0 * x10 * x14 * x34 * x4 * x9
    x36 = -x16 * (x12 * x33 + x35)
    x37 = x23 * x7
    x38 = x20 * x37 + x36
    x39 = x3 * x38
    x40 = x19 * x30
    x41 = x15 * x40
    x42 = x39 + x41
    x43 = x16 * (2.0 * x10 * x14 * x4 * x7 * x9 - x13 * x18)
    x44 = boys(7, x6)
    x45 = x24 * x44 + x43
    x46 = x3 * x45
    x47 = x18 * x19
    x48 = x30 * x47
    x49 = x25 * x5
    x50 = 0.5 * x4
    x51 = x50 * (x38 - x49)
    x52 = 2.0 * x28
    x53 = x29 * x52
    x54 = x53 * x7
    x55 = x23 * x47
    x56 = x3 * x55
    x57 = x54 + x56
    x58 = x10 * x21
    x59 = x58 * x8
    x60 = x27 * x59
    x61 = x50 * (2.0 * x10 * x14 * x15 * x19 * x4 * x9 - x60)
    x62 = x3 * x54
    x63 = x3 * x57 + x61 + x62
    x64 = x15 * x53
    x65 = x19 * x37
    x66 = x3 * x65
    x67 = x64 + x66
    x68 = x34 * x53
    x69 = x15 * x23
    x70 = x19 * x69
    x71 = x3 * x70
    x72 = x68 + x71
    x73 = x33 * x58
    x74 = -x50 * (x35 + x73)
    x75 = x3**2
    x76 = x37 * x75
    x77 = x74 + x76
    x78 = x28 * x77 + x3 * x63 - x4 * (x5 * x67 - x72)
    x79 = boys(2, x6)
    x80 = x16 * (2.0 * x10 * x14 * x4 * x79 * x9 - x13 * x34)
    x81 = x20 * x69 + x80
    x82 = x50 * (-x38 * x5 + x81)
    x83 = x3 * x32 + x52 * x67 + x82
    x84 = x16 * (2.0 * x10 * x14 * x4 * x9 * boys(1, x6) - x13 * x79)
    x85 = x23 * x34
    x86 = x20 * x85 + x84
    x87 = x50 * (-x5 * x81 + x86)
    x88 = x3 * x42 + x52 * x72 + x87
    x89 = 1.5 * x4
    x90 = x3 * x81 + x34 * x40
    x91 = x34 * x59
    x92 = x19 * x50 * (2.0 * x10 * x14 * x4 * x79 * x9 - x91)
    x93 = x3 * x68 + x3 * x72 + x92
    x94 = x19 * x73
    x95 = x50 * (2.0 * x10 * x14 * x19 * x34 * x4 * x9 - x94)
    x96 = x3 * x64
    x97 = x3 * x67 + x95 + x96
    x98 = x3 * x73
    x99 = 0.179587122125167 * da * db * numpy.sqrt(ax**6.5) * numpy.sqrt(bx**3.5)
    x100 = 6.79952275532455 * x99
    x101 = -x1 * (ax * A[1] + bx * B[1])
    x102 = -x101 - B[1]
    x103 = x59 * x7
    x104 = x102 * x103
    x105 = x104 * x3
    x106 = x50 * (2.0 * x10 * x102 * x14 * x15 * x4 * x9 - x104)
    x107 = x102 * x18
    x108 = x23 * x75
    x109 = x107 * x108
    x110 = x106 + x109
    x111 = x110 * x3 + x4 * (2.0 * x10 * x102 * x14 * x15 * x3 * x4 * x9 - x105)
    x112 = x102 * x73
    x113 = x50 * (2.0 * x10 * x102 * x14 * x34 * x4 * x9 - x112)
    x114 = x102 * x76
    x115 = x113 + x114
    x116 = x102 * x50 * (2.0 * x10 * x14 * x4 * x79 * x9 - x91)
    x117 = x102 * x69
    x118 = x116 + x117 * x75
    x119 = x102 * x54
    x120 = x102 * x56
    x121 = x119 + x120
    x122 = x102 * x64
    x123 = x102 * x66
    x124 = x122 + x123
    x125 = x107 * x53
    x126 = x23 * x44
    x127 = x126 * x19
    x128 = x127 * x3
    x129 = x102 * x128
    x130 = x47 * x59
    x131 = x102 * x130
    x132 = x50 * (2.0 * x10 * x102 * x14 * x19 * x4 * x7 * x9 - x131)
    x133 = x125 * x3
    x134 = x102 * x60
    x135 = x50 * (2.0 * x10 * x102 * x14 * x15 * x19 * x4 * x9 - x134)
    x136 = x102 * x62
    x137 = x121 * x3 + x135 + x136
    x138 = x102 * x94
    x139 = x50 * (2.0 * x10 * x102 * x14 * x19 * x34 * x4 * x9 - x138)
    x140 = x102 * x96
    x141 = x124 * x3 + x139 + x140
    x142 = x102 * x68
    x143 = x102 * x71 + x142
    x144 = x53 * x79
    x145 = x102 * x144
    x146 = x19 * x85
    x147 = x102 * x146
    x148 = 11.7771188794428 * x99
    x149 = -x1 * (ax * A[2] + bx * B[2])
    x150 = -x149 - B[2]
    x151 = x103 * x150
    x152 = x151 * x3
    x153 = x50 * (2.0 * x10 * x14 * x15 * x150 * x4 * x9 - x151)
    x154 = x150 * x18
    x155 = x108 * x154
    x156 = x153 + x155
    x157 = x156 * x3 + x4 * (2.0 * x10 * x14 * x15 * x150 * x3 * x4 * x9 - x152)
    x158 = x150 * x50 * (2.0 * x10 * x14 * x34 * x4 * x9 - x73)
    x159 = x150 * x76 + x158
    x160 = x150 * x50 * (2.0 * x10 * x14 * x4 * x79 * x9 - x91)
    x161 = x150 * x69
    x162 = x160 + x161 * x75
    x163 = x150 * x54
    x164 = x150 * x56
    x165 = x163 + x164
    x166 = x150 * x64
    x167 = x150 * x66
    x168 = x166 + x167
    x169 = x154 * x53
    x170 = x128 * x150
    x171 = x130 * x150
    x172 = x50 * (2.0 * x10 * x14 * x150 * x19 * x4 * x7 * x9 - x171)
    x173 = x169 * x3
    x174 = x150 * x60
    x175 = x50 * (2.0 * x10 * x14 * x15 * x150 * x19 * x4 * x9 - x174)
    x176 = x150 * x62
    x177 = x165 * x3 + x175 + x176
    x178 = x150 * x94
    x179 = x50 * (2.0 * x10 * x14 * x150 * x19 * x34 * x4 * x9 - x178)
    x180 = x150 * x96
    x181 = x168 * x3 + x179 + x180
    x182 = x150 * x68
    x183 = x150 * x71 + x182
    x184 = x144 * x150
    x185 = x146 * x150
    x186 = x102**2
    x187 = x186 * x37 + x36
    x188 = x186 * x23
    x189 = x17 + x18 * x188
    x190 = x189 * x5
    x191 = x190 * x3
    x192 = x188 * x44 + x43
    x193 = x192 * x75
    x194 = x50 * (x187 - x190)
    x195 = x189 * x75
    x196 = x186 * x69 + x80
    x197 = x50 * (-x187 * x5 + x196)
    x198 = x195 + x197
    x199 = x187 * x75
    x200 = x186 * x85 + x84
    x201 = x50 * (-x196 * x5 + x200)
    x202 = x199 + x201
    x203 = x107 * x150
    x204 = x203 * x59
    x205 = x50 * (2.0 * x10 * x102 * x14 * x150 * x4 * x7 * x9 - x204)
    x206 = x102 * x150
    x207 = x150 * x50 * (2.0 * x10 * x102 * x14 * x15 * x4 * x9 - x104)
    x208 = x108 * x203 + x207
    x209 = x150 * x50 * (2.0 * x10 * x102 * x14 * x34 * x4 * x9 - x112)
    x210 = x114 * x150 + x209
    x211 = x150**2
    x212 = x211 * x37 + x36
    x213 = x211 * x23
    x214 = x17 + x18 * x213
    x215 = x214 * x5
    x216 = x215 * x3
    x217 = x213 * x44 + x43
    x218 = x217 * x75
    x219 = x50 * (x212 - x215)
    x220 = x214 * x75
    x221 = x211 * x69 + x80
    x222 = x50 * (-x212 * x5 + x221)
    x223 = x220 + x222
    x224 = x212 * x75
    x225 = x211 * x85 + x84
    x226 = x50 * (-x221 * x5 + x225)
    x227 = x224 + x226
    x228 = -x101 - A[1]
    x229 = x228 * x26
    x230 = x228 * x31
    x231 = x229 + x230
    x232 = x228 * x39
    x233 = x228 * x41
    x234 = x232 + x233
    x235 = x228 * x46
    x236 = x228 * x48
    x237 = x228 * x49
    x238 = x4 * (x228 * x38 - x237)
    x239 = x228 * x54
    x240 = x228 * x56
    x241 = x239 + x240
    x242 = x228 * x60
    x243 = x4 * (2.0 * x10 * x14 * x15 * x19 * x228 * x4 * x9 - x242)
    x244 = x228 * x62
    x245 = x241 * x3 + 0.5 * x243 + x244
    x246 = x228 * x64
    x247 = x228 * x66 + x246
    x248 = x228 * x68
    x249 = x228 * x71 + x248
    x250 = x228 * x4 * (2.0 * x10 * x14 * x34 * x4 * x9 - x73)
    x251 = x228 * x4 * (-x5 * x81 + x86)
    x252 = x228 * x4 * (-x38 * x5 + x81)
    x253 = 20.3985682659737 * x99
    x254 = x117 * x228
    x255 = x254 + x68
    x256 = x102 * x228
    x257 = x256 * x37
    x258 = x257 + x64
    x259 = x258 * x5
    x260 = x107 * x23
    x261 = x228 * x260
    x262 = x261 + x54
    x263 = x4 * (x255 - x259)
    x264 = x262 * x75 + 0.5 * x263
    x265 = x258 * x28
    x266 = x27 * x53
    x267 = x256 * x55
    x268 = x266 + x267
    x269 = x265 + x268 * x3
    x270 = x255 * x28
    x271 = x19 * x64
    x272 = x102 * x228 * x65
    x273 = x271 + x272
    x274 = x270 + x273 * x3
    x275 = x262 * x28
    x276 = x47 * x53
    x277 = x127 * x256
    x278 = x276 + x277
    x279 = x4 * (-x268 * x5 + x273)
    x280 = x19 * x68
    x281 = x102 * x228 * x70 + x280
    x282 = x144 * x19
    x283 = x4 * (x147 * x228 - x281 * x5 + x282)
    x284 = x4 * (-x273 * x5 + x281)
    x285 = 35.3313566383285 * x99
    x286 = x228 * x4 * (2.0 * x10 * x14 * x15 * x150 * x4 * x9 - x151)
    x287 = x155 * x228 + 0.5 * x286
    x288 = x163 * x228
    x289 = x164 * x228 + x288
    x290 = x166 * x228
    x291 = x167 * x228 + x290
    x292 = x169 * x228
    x293 = x228 * x4 * (2.0 * x10 * x14 * x150 * x19 * x4 * x7 * x9 - x171)
    x294 = x228 * x4 * (2.0 * x10 * x14 * x150 * x19 * x34 * x4 * x9 - x178)
    x295 = x228 * x4 * (2.0 * x10 * x14 * x15 * x150 * x19 * x4 * x9 - x174)
    x296 = x187 * x228
    x297 = x102 * x30
    x298 = x15 * x297
    x299 = x296 + x298
    x300 = x189 * x228
    x301 = x297 * x7
    x302 = x300 + x301
    x303 = x302 * x5
    x304 = x192 * x228
    x305 = x107 * x30
    x306 = x304 + x305
    x307 = x4 * (x299 - x303)
    x308 = x196 * x228 + x297 * x34
    x309 = x4 * (x200 * x228 + x297 * x79 - x308 * x5)
    x310 = x4 * (-x299 * x5 + x308)
    x311 = x206 * x37
    x312 = x166 + x228 * x311
    x313 = x203 * x23
    x314 = x163 + x228 * x313
    x315 = x314 * x5
    x316 = x126 * x206
    x317 = x169 + x228 * x316
    x318 = x4 * (x312 - x315)
    x319 = x117 * x150
    x320 = x182 + x228 * x319
    x321 = x102 * x150 * x85
    x322 = x4 * (x184 + x228 * x321 - x320 * x5)
    x323 = x4 * (-x312 * x5 + x320)
    x324 = x228 * x4 * (x212 - x215)
    x325 = x228 * x4 * (-x221 * x5 + x225)
    x326 = x228 * x4 * (-x212 * x5 + x221)
    x327 = -x149 - A[2]
    x328 = x327 * (x26 + x31)
    x329 = x327 * (x39 + x41)
    x330 = x327 * x4 * (x38 - x49)
    x331 = 0.5 * x330
    x332 = x327 * x54
    x333 = x327 * x56 + x332
    x334 = x327 * x4 * (2.0 * x10 * x14 * x15 * x19 * x4 * x9 - x60)
    x335 = 0.5 * x334
    x336 = x3 * x333 + x327 * x62 + x335
    x337 = x327 * x64
    x338 = x327 * x66 + x337
    x339 = x327 * x68
    x340 = x327 * x71 + x339
    x341 = x327 * x4 * (2.0 * x10 * x14 * x34 * x4 * x9 - x73)
    x342 = 0.5 * x341
    x343 = x327 * x4 * (-x5 * x81 + x86)
    x344 = 0.5 * x343
    x345 = x327 * x38
    x346 = x4 * (x327 * x81 - x345 * x5)
    x347 = 0.5 * x346
    x348 = x327 * x4 * (2.0 * x10 * x102 * x14 * x15 * x4 * x9 - x104)
    x349 = 0.5 * x348
    x350 = x109 * x327 + x349
    x351 = x119 * x327
    x352 = x120 * x327 + x351
    x353 = x122 * x327
    x354 = x123 * x327 + x353
    x355 = x125 * x327
    x356 = x327 * x4 * (2.0 * x10 * x102 * x14 * x19 * x4 * x7 * x9 - x131)
    x357 = 0.5 * x356
    x358 = x327 * x4 * (2.0 * x10 * x102 * x14 * x19 * x34 * x4 * x9 - x138)
    x359 = 0.5 * x358
    x360 = x327 * x4 * (2.0 * x10 * x102 * x14 * x15 * x19 * x4 * x9 - x134)
    x361 = 0.5 * x360
    x362 = x161 * x327 + x68
    x363 = x150 * x327
    x364 = x363 * x37 + x64
    x365 = x364 * x5
    x366 = x154 * x23
    x367 = x327 * x366 + x54
    x368 = x367 * x75
    x369 = x4 * (x362 - x365)
    x370 = 0.5 * x369
    x371 = x368 + x370
    x372 = x28 * x364
    x373 = x266 + x363 * x55
    x374 = x3 * x373
    x375 = x372 + x374
    x376 = x28 * x362
    x377 = x150 * x327 * x65 + x271
    x378 = x3 * x377
    x379 = x376 + x378
    x380 = x28 * x367
    x381 = x127 * x363 + x276
    x382 = x3 * x381
    x383 = x373 * x5
    x384 = x4 * (x377 - x383)
    x385 = 0.5 * x384
    x386 = x3 * x380
    x387 = x150 * x327 * x70 + x280
    x388 = x4 * (x185 * x327 + x282 - x387 * x5)
    x389 = 0.5 * x388
    x390 = x4 * (-x377 * x5 + x387)
    x391 = 0.5 * x390
    x392 = x327 * x4 * (x187 - x190)
    x393 = 0.5 * x392
    x394 = x327 * x4 * (-x196 * x5 + x200)
    x395 = 0.5 * x394
    x396 = x327 * x4 * (-x187 * x5 + x196)
    x397 = 0.5 * x396
    x398 = x122 + x311 * x327
    x399 = x119 + x313 * x327
    x400 = x399 * x5
    x401 = x125 + x316 * x327
    x402 = x4 * (x398 - x400)
    x403 = 0.5 * x402
    x404 = x142 + x319 * x327
    x405 = x4 * (x145 + x321 * x327 - x404 * x5)
    x406 = 0.5 * x405
    x407 = x4 * (-x398 * x5 + x404)
    x408 = 0.5 * x407
    x409 = x150 * x30
    x410 = x15 * x409 + x212 * x327
    x411 = x214 * x327 + x409 * x7
    x412 = x411 * x5
    x413 = x3 * x412
    x414 = x154 * x30 + x217 * x327
    x415 = x414 * x75
    x416 = x4 * (x410 - x412)
    x417 = 0.5 * x416
    x418 = x221 * x327 + x34 * x409
    x419 = x4 * (x225 * x327 + x409 * x79 - x418 * x5)
    x420 = 0.5 * x419
    x421 = x4 * (-x410 * x5 + x418)
    x422 = 0.5 * x421
    x423 = x228**2
    x424 = x38 * x423 + x87
    x425 = x423 * x70 + x92
    x426 = x28 * x425
    x427 = x25 * x423
    x428 = x427 + x82
    x429 = x37 * x423
    x430 = x19 * x429
    x431 = x430 + x95
    x432 = x28 * x431
    x433 = x423 * x45
    x434 = x433 + x51
    x435 = x423 * x55
    x436 = x435 + x61
    x437 = x28 * x436
    x438 = x424 - x428 * x5
    x439 = x429 + x74
    x440 = x28 * x439
    x441 = x3 * x436 + x440
    x442 = x425 - x431 * x5
    x443 = 31.1593277158494 * x99
    x444 = x116 + x228 * x255 + x248
    x445 = x28 * x444
    x446 = x19 * x246
    x447 = x139 + x228 * x273 + x446
    x448 = x113 + x228 * x258 + x246
    x449 = x28 * x448
    x450 = x228 * x266
    x451 = x135 + x228 * x268 + x450
    x452 = x106 + x228 * x262 + x239
    x453 = x444 - x448 * x5
    x454 = x28 * x452
    x455 = x228 * x276
    x456 = x132 + x228 * x278 + x455
    x457 = x447 - x451 * x5
    x458 = 53.9695387335403 * x99
    x459 = x160 + x161 * x423
    x460 = x28 * x459
    x461 = x150 * x430 + x179
    x462 = x150 * x429 + x158
    x463 = x28 * x462
    x464 = x150 * x435 + x175
    x465 = x153 + x366 * x423
    x466 = x459 - x462 * x5
    x467 = x28 * x465
    x468 = x127 * x150 * x423 + x172
    x469 = x461 - x464 * x5
    x470 = x201 + x228 * x299 + 2.0 * x270
    x471 = x197 + x228 * x302 + 2.0 * x265
    x472 = x471 * x5
    x473 = x194 + x228 * x306 + 2.0 * x275
    x474 = x470 - x472
    x475 = x209 + x228 * x312 + x290
    x476 = x207 + x228 * x314 + x288
    x477 = x476 * x5
    x478 = x205 + x228 * x317 + x292
    x479 = x475 - x477
    x480 = x212 * x423 + x226
    x481 = x214 * x423 + x222
    x482 = x481 * x5
    x483 = x217 * x423 + x219
    x484 = x480 - x482
    x485 = x327 * x4 * (x228 * x38 - x237)
    x486 = x239 * x327
    x487 = x240 * x327 + x486
    x488 = x327 * x4 * (2.0 * x10 * x14 * x15 * x19 * x228 * x4 * x9 - x242)
    x489 = x254 * x327 + x339
    x490 = x28 * x489
    x491 = x271 * x327
    x492 = x272 * x327 + x491
    x493 = x257 * x327 + x337
    x494 = x28 * x493
    x495 = x266 * x327
    x496 = x267 * x327 + x495
    x497 = x261 * x327 + x332
    x498 = x4 * (x489 - x493 * x5)
    x499 = x28 * x497
    x500 = x276 * x327
    x501 = x277 * x327 + x500
    x502 = x4 * (x492 - x496 * x5)
    x503 = 93.4779831475484 * x99
    x504 = x228 * x376
    x505 = x228 * x372
    x506 = x228 * x4 * (x362 - x365)
    x507 = x228 * x380
    x508 = x228 * x4 * (x377 - x383)
    x509 = x327 * (x296 + x298)
    x510 = x327 * (x300 + x301)
    x511 = x5 * x510
    x512 = x327 * (x304 + x305)
    x513 = x4 * (x509 - x511)
    x514 = x228 * x398 + x376
    x515 = x228 * x399 + x372
    x516 = x5 * x515
    x517 = x228 * x401 + x380
    x518 = x4 * (x514 - x516)
    x519 = x228 * x4 * (x410 - x412)
    x520 = x327**2
    x521 = x38 * x520 + x87
    x522 = x520 * x70 + x92
    x523 = x28 * x522
    x524 = x25 * x520 + x82
    x525 = x37 * x520
    x526 = x19 * x525
    x527 = x526 + x95
    x528 = x28 * x527
    x529 = x45 * x520 + x51
    x530 = x3 * x529
    x531 = x520 * x55
    x532 = x531 + x61
    x533 = x28 * x532
    x534 = 2.0 * x533
    x535 = x5 * x524
    x536 = x521 - x535
    x537 = x50 * x536
    x538 = x525 + x74
    x539 = x28 * x538
    x540 = x3 * x532
    x541 = x539 + x540
    x542 = -x5 * x527 + x522
    x543 = x50 * x542
    x544 = x116 + x117 * x520
    x545 = x28 * x544
    x546 = x102 * x526 + x139
    x547 = x102 * x525 + x113
    x548 = x28 * x547
    x549 = x102 * x531 + x135
    x550 = x106 + x260 * x520
    x551 = -x5 * x547 + x544
    x552 = x50 * x551
    x553 = x28 * x550
    x554 = x102 * x127 * x520 + x132
    x555 = -x5 * x549 + x546
    x556 = x50 * x555
    x557 = x160 + x327 * x362 + x339
    x558 = x28 * x557
    x559 = x179 + x327 * x377 + x491
    x560 = x158 + x327 * x364 + x337
    x561 = x28 * x560
    x562 = x175 + x327 * x373 + x495
    x563 = x153 + x327 * x367 + x332
    x564 = -x5 * x560 + x557
    x565 = x50 * x564
    x566 = x28 * x563
    x567 = x172 + x327 * x381 + x500
    x568 = x3 * x567
    x569 = x5 * x562
    x570 = x559 - x569
    x571 = x50 * x570
    x572 = x3 * x566
    x573 = x187 * x520 + x201
    x574 = x189 * x520 + x197
    x575 = x5 * x574
    x576 = x192 * x520 + x194
    x577 = x573 - x575
    x578 = x50 * x577
    x579 = x209 + x327 * x398 + x353
    x580 = x207 + x327 * x399 + x351
    x581 = x5 * x580
    x582 = x205 + x327 * x401 + x355
    x583 = x579 - x581
    x584 = x50 * x583
    x585 = x226 + x327 * x410 + 2.0 * x376
    x586 = x222 + x327 * x411 + 2.0 * x372
    x587 = x5 * x586
    x588 = x219 + x327 * x414 + 2.0 * x380
    x589 = x588 * x75
    x590 = x585 - x587
    x591 = x50 * x590
    x592 = x228 * x434 + x238
    x593 = x228 * x436 + x243
    x594 = x28 * x593
    x595 = x228 * x424 + x251 - x5 * (x228 * x428 + x252)
    x596 = x28 * (x228 * x439 + x250)
    x597 = 31.1593277158495 * x99
    x598 = x228 * x452 + x263 + x440
    x599 = x28 * x598
    x600 = x228 * x456 + x279 + x437
    x601 = x228 * x447 + x283 + x426 - x5 * (x228 * x451 + x284 + x432)
    x602 = 53.9695387335404 * x99
    x603 = x228 * x465 + x286
    x604 = x28 * x603
    x605 = x228 * x468 + x293
    x606 = x228 * x461 + x294 - x5 * (x228 * x464 + x295)
    x607 = x228 * x473 + x307 + 2.0 * x454
    x608 = x228 * x470 + x309 + 2.0 * x445 - x5 * (x228 * x471 + x310 + 2.0 * x449)
    x609 = x228 * x478 + x318 + x467
    x610 = x228 * x475 + x322 + x460 - x5 * (x228 * x476 + x323 + x463)
    x611 = x228 * x483 + x324
    x612 = x228 * x480 + x325 - x5 * (x228 * x481 + x326)
    x613 = x327 * x433 + x331
    x614 = x327 * x435 + x335
    x615 = x28 * x614
    x616 = x344 + x345 * x423 - x5 * (x327 * x427 + x347)
    x617 = x28 * (x327 * x429 + x342)
    x618 = 69.6743749058326 * x99
    x619 = x228 * x497 + x349 + x486
    x620 = x28 * x619
    x621 = x228 * x501 + x327 * x455 + x357
    x622 = x228 * x492 + x327 * x446 + x359 - x5 * (x228 * x496 + x327 * x450 + x361)
    x623 = 120.679557322504 * x99
    x624 = x367 * x423 + x370
    x625 = x28 * x624
    x626 = x381 * x423 + x385
    x627 = x377 * x423 + x389 - x5 * (x373 * x423 + x391)
    x628 = x228 * x512 + x393 + 2.0 * x499
    x629 = x228 * x509 + x395 + 2.0 * x490 - x5 * (x228 * x510 + x397 + 2.0 * x494)
    x630 = x228 * x517 + x403 + x507
    x631 = x228 * x514 + x406 - x5 * (x228 * x515 + x408 + x505) + x504
    x632 = x414 * x423 + x417
    x633 = x410 * x423 + x420 - x5 * (x411 * x423 + x422)
    x634 = x228 * x4 * (x521 - x535)
    x635 = x228 * x539
    x636 = x228 * x550 + x539
    x637 = x28 * x636
    x638 = x228 * x554 + x533
    x639 = x4 * (x228 * x546 - x5 * (x228 * x549 + x528) + x523)
    x640 = x228 * x566
    x641 = x228 * x4 * (x559 - x569)
    x642 = x228 * x576 + 2.0 * x553
    x643 = x4 * (x228 * x573 - x5 * (x228 * x574 + 2.0 * x548) + 2.0 * x545)
    x644 = x228 * x582 + x566
    x645 = x4 * (x228 * x579 - x5 * (x228 * x580 + x561) + x558)
    x646 = x228 * x4 * (x585 - x587)
    x647 = x327 * x529 + x330
    x648 = x3 * x647
    x649 = x327 * x532 + x334
    x650 = x28 * x649
    x651 = 2.0 * x650
    x652 = x327 * x521 + x343 - x5 * (x327 * x524 + x346)
    x653 = x50 * x652
    x654 = x28 * (x327 * x538 + x341)
    x655 = x327 * x550 + x348
    x656 = x28 * x655
    x657 = x327 * x554 + x356
    x658 = x327 * x546 + x358 - x5 * (x327 * x549 + x360)
    x659 = x50 * x658
    x660 = x327 * x563 + x369 + x539
    x661 = x28 * x660
    x662 = x327 * x567 + x384 + x533
    x663 = x3 * x662
    x664 = x327 * x559 + x388 - x5 * (x327 * x562 + x390 + x528) + x523
    x665 = x50 * x664
    x666 = x327 * x576 + x392
    x667 = x327 * x573 + x394 - x5 * (x327 * x574 + x396)
    x668 = x50 * x667
    x669 = x327 * x582 + x402 + x553
    x670 = x327 * x579 + x405 - x5 * (x327 * x580 + x407 + x548) + x545
    x671 = x50 * x670
    x672 = x327 * x588 + x416 + 2.0 * x566
    x673 = x327 * x585 + x419 - x5 * (x327 * x586 + x421 + 2.0 * x561) + 2.0 * x558
    x674 = x50 * x673
    x675 = x228 * x592 + x438 * x89
    x676 = x28 * (x228 * x593 + x442 * x89)
    x677 = x28 * (x228 * x598 + x453 * x89 + x596)
    x678 = x228 * x600 + x457 * x89 + x594
    x679 = x28 * (x228 * x603 + x466 * x89)
    x680 = x228 * x605 + x469 * x89
    x681 = x228 * x607 + x474 * x89 + 2.0 * x599
    x682 = x253 * x3
    x683 = x228 * x609 + x479 * x89 + x604
    x684 = x285 * x3
    x685 = x228 * x611 + x484 * x89
    x686 = x228 * x613 + x485
    x687 = x28 * (x228 * x614 + x488)
    x688 = x28 * (x228 * x619 + x498 + x617)
    x689 = x228 * x621 + x502 + x615
    x690 = x28 * (x228 * x624 + x506)
    x691 = x228 * x626 + x508
    x692 = x228 * x628 + x513 + 2.0 * x620
    x693 = x3 * x458
    x694 = x228 * x630 + x518 + x625
    x695 = x3 * x503
    x696 = x228 * x632 + x519
    x697 = x423 * x529 + x537
    x698 = x28 * (x423 * x532 + x543)
    x699 = x28 * (x228 * x636 + x552 + x635)
    x700 = x228 * x533 + x228 * x638 + x556
    x701 = x28 * (x423 * x563 + x565)
    x702 = x423 * x567 + x571
    x703 = x228 * x642 + x578 + 2.0 * x637
    x704 = x3 * x618
    x705 = x228 * x644 + x584 + x640
    x706 = x423 * x588 + x591
    x707 = x28 * (x228 * x655 + x654)
    x708 = x228 * x657 + x650
    x709 = x228 * x661
    x710 = x228 * x666 + 2.0 * x656
    x711 = x228 * x669 + x661
    x712 = x327 * x647 + x536 * x89
    x713 = x28 * (x327 * x649 + x542 * x89)
    x714 = x28 * (x327 * x655 + x551 * x89)
    x715 = x327 * x657 + x555 * x89
    x716 = x28 * (x327 * x660 + x564 * x89 + x654)
    x717 = x327 * x662 + x570 * x89 + x650
    x718 = x327 * x666 + x577 * x89
    x719 = x327 * x669 + x583 * x89 + x656
    x720 = x327 * x672 + x590 * x89 + 2.0 * x661
    x721 = x228 * x253

    # 126 item(s)
    result[0, 0] = numpy.sum(
        x100
        * (
            x21
            * (
                x3 * x88
                + x4 * (x3 * x86 + x40 * x79 - x5 * x90)
                - x5 * (x3 * x83 - x4 * (x42 * x5 - x90) + x52 * x97)
                + x52 * x93
            )
            + x3
            * (
                x3
                * (
                    x3 * (x3 * (x46 + x48) + x51 + x52 * x57)
                    - x4 * (x32 * x5 - x42)
                    + x52 * x63
                )
                + x52 * x78
                - x89 * (x5 * x83 - x88)
            )
            + x52
            * (
                x28 * (x3 * x77 + x4 * (2.0 * x10 * x14 * x3 * x34 * x4 * x9 - x98))
                + x3 * x78
                - x89 * (x5 * x97 - x93)
            )
        )
    )
    result[0, 1] = numpy.sum(
        x148
        * (
            x21
            * (
                x118 * x28
                + x141 * x3
                + x4 * (-x143 * x5 + x145 + x147 * x3)
                - x5 * (x115 * x28 + x137 * x3 - x4 * (x124 * x5 - x143))
            )
            + x28 * (x111 * x3 - x89 * (x115 * x5 - x118))
            + x3
            * (
                x111 * x28
                + x3
                * (
                    x110 * x28
                    + x3 * (x132 + x133 + x3 * (x125 + x129))
                    - x4 * (x121 * x5 - x124)
                )
                - x89 * (x137 * x5 - x141)
            )
        )
    )
    result[0, 2] = numpy.sum(
        x148
        * (
            x21
            * (
                x162 * x28
                + x181 * x3
                + x4 * (-x183 * x5 + x184 + x185 * x3)
                - x5 * (x159 * x28 + x177 * x3 - x4 * (x168 * x5 - x183))
            )
            + x28 * (x157 * x3 - x89 * (x159 * x5 - x162))
            + x3
            * (
                x157 * x28
                + x3
                * (
                    x156 * x28
                    + x3 * (x172 + x173 + x3 * (x169 + x170))
                    - x4 * (x165 * x5 - x168)
                )
                - x89 * (x177 * x5 - x181)
            )
        )
    )
    result[0, 3] = numpy.sum(
        x100
        * x3
        * (
            -x21
            * (-x202 + x4 * (x196 * x5 - x200) + x5 * (x198 - x4 * (x187 * x5 - x196)))
            + x3 * (x3 * (x193 + x194) + x4 * (x187 * x3 - x191))
            - x89 * (x198 * x5 - x202)
        )
    )
    result[0, 4] = numpy.sum(
        x148
        * (
            x21
            * (
                x210 * x3
                + x4 * (2.0 * x10 * x102 * x14 * x150 * x3 * x34 * x4 * x9 - x206 * x98)
                - x5
                * (
                    x150 * x4 * (2.0 * x10 * x102 * x14 * x15 * x3 * x4 * x9 - x105)
                    + x208 * x3
                )
            )
            + x3
            * (
                x3**2
                * (
                    x108 * x206 * x44
                    + x205
                    + x4 * (2.0 * x10 * x102 * x14 * x150 * x4 * x7 * x9 - x204)
                )
                - x89 * (x208 * x5 - x210)
            )
        )
    )
    result[0, 5] = numpy.sum(
        x100
        * x3
        * (
            -x21
            * (-x227 + x4 * (x221 * x5 - x225) + x5 * (x223 - x4 * (x212 * x5 - x221)))
            + x3 * (x3 * (x218 + x219) + x4 * (x212 * x3 - x216))
            - x89 * (x223 * x5 - x227)
        )
    )
    result[1, 0] = numpy.sum(
        0.5
        * x253
        * (
            x3
            * (
                2.0 * x245 * x52
                + x3 * (x238 + 2.0 * x241 * x52 + 2.0 * x3 * (x235 + x236))
                - 2.0 * x4 * (x231 * x5 - x234)
            )
            + x52
            * (
                2.0 * x245 * x3
                + x28 * (2.0 * x228 * x76 + x250)
                - 2.0 * x4 * (x247 * x5 - x249)
            )
            + x89
            * (
                2.0 * x234 * x3
                + 2.0 * x249 * x52
                + x251
                - x5 * (2.0 * x231 * x3 + 2.0 * x247 * x52 + x252)
            )
        )
    )
    result[1, 1] = numpy.sum(
        0.5
        * x285
        * (
            2.0 * x28 * x3 * (x264 + x4 * (x255 - x259))
            + x3
            * (
                2.0 * x264 * x28
                + x3 * (2.0 * x275 * x3 + x279 + 2.0 * x3 * (x275 + x278 * x3))
                - 2.0 * x4 * (x269 * x5 - x274)
            )
            + x89
            * (
                2.0 * x255 * x28 * x3
                + 2.0 * x274 * x3
                + x283
                - x5 * (2.0 * x265 * x3 + 2.0 * x269 * x3 + x284)
            )
        )
    )
    result[1, 2] = numpy.sum(
        0.5
        * x285
        * (
            2.0
            * x28
            * (
                x228 * x4 * (2.0 * x10 * x14 * x15 * x150 * x3 * x4 * x9 - x152)
                + x287 * x3
            )
            + x3
            * (
                2.0 * x28 * x287
                + x3 * (2.0 * x173 * x228 + x293 + 2.0 * x3 * (x170 * x228 + x292))
                - 2.0 * x4 * (x289 * x5 - x291)
            )
            + x89
            * (
                2.0 * x180 * x228
                + 2.0 * x291 * x3
                + x294
                - x5 * (2.0 * x176 * x228 + 2.0 * x289 * x3 + x295)
            )
        )
    )
    result[1, 3] = numpy.sum(
        0.5
        * x253
        * (
            x3**2 * (2.0 * x306 * x75 + x307 + 2.0 * x4 * (x299 - x303))
            + x89 * (2.0 * x299 * x75 + x309 - x5 * (2.0 * x302 * x75 + x310))
        )
    )
    result[1, 4] = numpy.sum(
        0.5
        * x285
        * (
            x3**2 * (2.0 * x317 * x75 + x318 + 2.0 * x4 * (x312 - x315))
            + x89 * (2.0 * x312 * x75 + x322 - x5 * (2.0 * x314 * x75 + x323))
        )
    )
    result[1, 5] = numpy.sum(
        0.5
        * x253
        * (
            x3 * (2.0 * x228 * x4 * (x212 * x3 - x216) + x3 * (2.0 * x218 * x228 + x324))
            + x89 * (2.0 * x224 * x228 + x325 - x5 * (2.0 * x220 * x228 + x326))
        )
    )
    result[2, 0] = numpy.sum(
        x253
        * (
            x3
            * (
                x3 * (x3 * x327 * (x46 + x48) + x331 + x333 * x52)
                + x336 * x52
                - x4 * (x328 * x5 - x329)
            )
            + x52 * (x28 * (x327 * x76 + x342) + x3 * x336 - x4 * (x338 * x5 - x340))
            + x89 * (x3 * x329 + x340 * x52 + x344 - x5 * (x3 * x328 + x338 * x52 + x347))
        )
    )
    result[2, 1] = numpy.sum(
        x285
        * (
            x28
            * (
                x3 * x350
                + x327 * x4 * (2.0 * x10 * x102 * x14 * x15 * x3 * x4 * x9 - x105)
            )
            + x3
            * (
                x28 * x350
                + x3 * (x133 * x327 + x3 * (x129 * x327 + x355) + x357)
                - x4 * (x352 * x5 - x354)
            )
            + x89
            * (x140 * x327 + x3 * x354 + x359 - x5 * (x136 * x327 + x3 * x352 + x361))
        )
    )
    result[2, 2] = numpy.sum(
        x285
        * (
            x28 * x3 * (x371 + x4 * (x362 - x365))
            + x3
            * (
                x28 * x371
                + x3 * (x3 * (x380 + x382) + x385 + x386)
                - x4 * (x375 * x5 - x379)
            )
            + x89
            * (x28 * x3 * x362 + x3 * x379 + x389 - x5 * (x3 * x372 + x3 * x375 + x391))
        )
    )
    result[2, 3] = numpy.sum(
        x253
        * (
            x3 * (x3 * (x193 * x327 + x393) + x327 * x4 * (x187 * x3 - x191))
            + x89 * (x199 * x327 + x395 - x5 * (x195 * x327 + x397))
        )
    )
    result[2, 4] = numpy.sum(
        x285
        * (
            x3**2 * (x4 * (x398 - x400) + x401 * x75 + x403)
            + x89 * (x398 * x75 + x406 - x5 * (x399 * x75 + x408))
        )
    )
    result[2, 5] = numpy.sum(
        x253
        * (
            x3 * (x3 * (x415 + x417) + x4 * (x3 * x410 - x413))
            + x89 * (x410 * x75 + x420 - x5 * (x411 * x75 + x422))
        )
    )
    result[3, 0] = numpy.sum(
        x443
        * (
            x3 * (x3 * (x3 * x434 + 2.0 * x437) + x438 * x50 + x441 * x52)
            + x4 * (x3 * x424 + 2.0 * x426 - x5 * (x3 * x428 + 2.0 * x432))
            + x52 * (x3 * x440 + x3 * x441 + x442 * x50)
        )
    )
    result[3, 1] = numpy.sum(
        x458
        * (
            x28 * (x452 * x75 + x453 * x50)
            + x3 * (x3 * x454 + x3 * (x3 * x456 + x454) + x457 * x50)
            + x4 * (x3 * x447 + x445 - x5 * (x3 * x451 + x449))
        )
    )
    result[3, 2] = numpy.sum(
        x458
        * (
            x28 * (x465 * x75 + x466 * x50)
            + x3 * (x3 * x467 + x3 * (x3 * x468 + x467) + x469 * x50)
            + x4 * (x3 * x461 + x460 - x5 * (x3 * x464 + x463))
        )
    )
    result[3, 3] = numpy.sum(x3 * x443 * (x4 * (x470 - x472) + x473 * x75 + x474 * x50))
    result[3, 4] = numpy.sum(x3 * x458 * (x4 * (x475 - x477) + x478 * x75 + x479 * x50))
    result[3, 5] = numpy.sum(x3 * x443 * (x4 * (x480 - x482) + x483 * x75 + x484 * x50))
    result[4, 0] = numpy.sum(
        0.5
        * x458
        * (
            x3 * (2.0 * x3 * x327 * (x235 + x236) + x485 + 2.0 * x487 * x52)
            + 2.0 * x327 * x4 * (x232 + x233 - x5 * (x229 + x230))
            + x52 * (2.0 * x244 * x327 + 2.0 * x3 * x487 + x488)
        )
    )
    result[4, 1] = numpy.sum(
        0.5
        * x503
        * (
            x28 * (2.0 * x497 * x75 + x498)
            + x3 * (2.0 * x3 * x499 + 2.0 * x3 * (x3 * x501 + x499) + x502)
            + 2.0 * x4 * (x3 * x492 + x490 - x5 * (x3 * x496 + x494))
        )
    )
    result[4, 2] = numpy.sum(
        0.5
        * x503
        * (
            x28 * (2.0 * x228 * x368 + x506)
            + x3 * (2.0 * x228 * x386 + 2.0 * x3 * (x228 * x382 + x507) + x508)
            + 2.0 * x4 * (x228 * x378 - x5 * (x228 * x374 + x505) + x504)
        )
    )
    result[4, 3] = numpy.sum(
        0.5 * x3 * x458 * (2.0 * x4 * (x509 - x511) + 2.0 * x512 * x75 + x513)
    )
    result[4, 4] = numpy.sum(
        0.5 * x3 * x503 * (2.0 * x4 * (x514 - x516) + 2.0 * x517 * x75 + x518)
    )
    result[4, 5] = numpy.sum(
        0.5
        * x458
        * (2.0 * x228 * x4 * (x3 * x410 - x413) + x3 * (2.0 * x228 * x415 + x519))
    )
    result[5, 0] = numpy.sum(
        x443
        * (
            x3 * (x3 * (x530 + x534) + x52 * x541 + x537)
            + x4 * (x3 * x521 - x5 * (x3 * x524 + 2.0 * x528) + 2.0 * x523)
            + x52 * (x3 * x539 + x3 * x541 + x543)
        )
    )
    result[5, 1] = numpy.sum(
        x458
        * (
            x28 * (x550 * x75 + x552)
            + x3 * (x3 * x553 + x3 * (x3 * x554 + x553) + x556)
            + x4 * (x3 * x546 - x5 * (x3 * x549 + x548) + x545)
        )
    )
    result[5, 2] = numpy.sum(
        x458
        * (
            x28 * (x563 * x75 + x565)
            + x3 * (x3 * (x566 + x568) + x571 + x572)
            + x4 * (x3 * x559 - x5 * (x3 * x562 + x561) + x558)
        )
    )
    result[5, 3] = numpy.sum(x3 * x443 * (x4 * (x573 - x575) + x576 * x75 + x578))
    result[5, 4] = numpy.sum(x3 * x458 * (x4 * (x579 - x581) + x582 * x75 + x584))
    result[5, 5] = numpy.sum(x3 * x443 * (x4 * (x585 - x587) + x589 + x591))
    result[6, 0] = numpy.sum(
        x597 * (x3 * (x3 * x592 + 2.0 * x594) + x50 * x595 + x52 * (x3 * x593 + x596))
    )
    result[6, 1] = numpy.sum(x602 * (x3 * x599 + x3 * (x3 * x600 + x599) + x50 * x601))
    result[6, 2] = numpy.sum(x602 * (x3 * x604 + x3 * (x3 * x605 + x604) + x50 * x606))
    result[6, 3] = numpy.sum(x597 * (x50 * x608 + x607 * x75))
    result[6, 4] = numpy.sum(x602 * (x50 * x610 + x609 * x75))
    result[6, 5] = numpy.sum(x597 * (x50 * x612 + x611 * x75))
    result[7, 0] = numpy.sum(
        x618 * (x3 * (x3 * x613 + 2.0 * x615) + x50 * x616 + x52 * (x3 * x614 + x617))
    )
    result[7, 1] = numpy.sum(x623 * (x3 * x620 + x3 * (x3 * x621 + x620) + x50 * x622))
    result[7, 2] = numpy.sum(x623 * (x3 * x625 + x3 * (x3 * x626 + x625) + x50 * x627))
    result[7, 3] = numpy.sum(x618 * (x50 * x629 + x628 * x75))
    result[7, 4] = numpy.sum(x623 * (x50 * x631 + x630 * x75))
    result[7, 5] = numpy.sum(x618 * (x50 * x633 + x632 * x75))
    result[8, 0] = numpy.sum(
        0.5
        * x618
        * (2.0 * x228 * x3 * (x530 + x534) + 2.0 * x52 * (x228 * x540 + x635) + x634)
    )
    result[8, 1] = numpy.sum(0.5 * x623 * (2.0 * x3**2 * x638 + 4.0 * x3 * x637 + x639))
    result[8, 2] = numpy.sum(
        0.5 * x623 * (2.0 * x228 * x572 + 2.0 * x3 * (x228 * x568 + x640) + x641)
    )
    result[8, 3] = numpy.sum(0.5 * x618 * (2.0 * x642 * x75 + x643))
    result[8, 4] = numpy.sum(0.5 * x623 * (2.0 * x644 * x75 + x645))
    result[8, 5] = numpy.sum(0.5 * x618 * (2.0 * x228 * x589 + x646))
    result[9, 0] = numpy.sum(
        x597 * (x3 * (x648 + x651) + x52 * (x3 * x649 + x654) + x653)
    )
    result[9, 1] = numpy.sum(x602 * (x3 * x656 + x3 * (x3 * x657 + x656) + x659))
    result[9, 2] = numpy.sum(x602 * (x3 * x661 + x3 * (x661 + x663) + x665))
    result[9, 3] = numpy.sum(x597 * (x666 * x75 + x668))
    result[9, 4] = numpy.sum(x602 * (x669 * x75 + x671))
    result[9, 5] = numpy.sum(x597 * (x672 * x75 + x674))
    result[10, 0] = numpy.sum(x253 * (x3 * x675 + 2.0 * x676))
    result[10, 1] = numpy.sum(x285 * (x3 * x678 + x677))
    result[10, 2] = numpy.sum(x285 * (x3 * x680 + x679))
    result[10, 3] = numpy.sum(x681 * x682)
    result[10, 4] = numpy.sum(x683 * x684)
    result[10, 5] = numpy.sum(x682 * x685)
    result[11, 0] = numpy.sum(x458 * (x3 * x686 + 2.0 * x687))
    result[11, 1] = numpy.sum(x503 * (x3 * x689 + x688))
    result[11, 2] = numpy.sum(x503 * (x3 * x691 + x690))
    result[11, 3] = numpy.sum(x692 * x693)
    result[11, 4] = numpy.sum(x694 * x695)
    result[11, 5] = numpy.sum(x693 * x696)
    result[12, 0] = numpy.sum(x618 * (x3 * x697 + 2.0 * x698))
    result[12, 1] = numpy.sum(x623 * (x3 * x700 + x699))
    result[12, 2] = numpy.sum(x623 * (x3 * x702 + x701))
    result[12, 3] = numpy.sum(x703 * x704)
    result[12, 4] = numpy.sum(x3 * x623 * x705)
    result[12, 5] = numpy.sum(x704 * x706)
    result[13, 0] = numpy.sum(x228 * x458 * (x648 + x651))
    result[13, 1] = numpy.sum(x503 * (x3 * x708 + x707))
    result[13, 2] = numpy.sum(x503 * (x228 * x663 + x709))
    result[13, 3] = numpy.sum(x693 * x710)
    result[13, 4] = numpy.sum(x695 * x711)
    result[13, 5] = numpy.sum(x228 * x672 * x693)
    result[14, 0] = numpy.sum(x253 * (x3 * x712 + 2.0 * x713))
    result[14, 1] = numpy.sum(x285 * (x3 * x715 + x714))
    result[14, 2] = numpy.sum(x285 * (x3 * x717 + x716))
    result[14, 3] = numpy.sum(x682 * x718)
    result[14, 4] = numpy.sum(x684 * x719)
    result[14, 5] = numpy.sum(x682 * x720)
    result[15, 0] = numpy.sum(x100 * (x21 * x595 + x228 * x675))
    result[15, 1] = numpy.sum(x148 * (x21 * x601 + x228 * x678 + x676))
    result[15, 2] = numpy.sum(x148 * (x21 * x606 + x228 * x680))
    result[15, 3] = numpy.sum(x100 * (x21 * x608 + x228 * x681 + 2.0 * x677))
    result[15, 4] = numpy.sum(x148 * (x21 * x610 + x228 * x683 + x679))
    result[15, 5] = numpy.sum(x100 * (x21 * x612 + x228 * x685))
    result[16, 0] = numpy.sum(x253 * (x228 * x686 + x616 * x89))
    result[16, 1] = numpy.sum(x285 * (x228 * x689 + x622 * x89 + x687))
    result[16, 2] = numpy.sum(x285 * (x228 * x691 + x627 * x89))
    result[16, 3] = numpy.sum(x253 * (x228 * x692 + x629 * x89 + 2.0 * x688))
    result[16, 4] = numpy.sum(x285 * (x228 * x694 + x631 * x89 + x690))
    result[16, 5] = numpy.sum(x253 * (x228 * x696 + x633 * x89))
    result[17, 0] = numpy.sum(x443 * (x228 * x697 + x634))
    result[17, 1] = numpy.sum(x458 * (x228 * x700 + x639 + x698))
    result[17, 2] = numpy.sum(x458 * (x228 * x702 + x641))
    result[17, 3] = numpy.sum(x443 * (x228 * x703 + x643 + 2.0 * x699))
    result[17, 4] = numpy.sum(x458 * (x228 * x705 + x645 + x701))
    result[17, 5] = numpy.sum(x443 * (x228 * x706 + x646))
    result[18, 0] = numpy.sum(x597 * (x423 * x647 + x653))
    result[18, 1] = numpy.sum(x602 * (x228 * x650 + x228 * x708 + x659))
    result[18, 2] = numpy.sum(x602 * (x423 * x662 + x665))
    result[18, 3] = numpy.sum(x597 * (x228 * x710 + x668 + 2.0 * x707))
    result[18, 4] = numpy.sum(x602 * (x228 * x711 + x671 + x709))
    result[18, 5] = numpy.sum(x597 * (x423 * x672 + x674))
    result[19, 0] = numpy.sum(x712 * x721)
    result[19, 1] = numpy.sum(x285 * (x228 * x715 + x713))
    result[19, 2] = numpy.sum(x228 * x285 * x717)
    result[19, 3] = numpy.sum(x253 * (x228 * x718 + 2.0 * x714))
    result[19, 4] = numpy.sum(x285 * (x228 * x719 + x716))
    result[19, 5] = numpy.sum(x720 * x721)
    result[20, 0] = numpy.sum(x100 * (x21 * x652 + x327 * x712))
    result[20, 1] = numpy.sum(x148 * (x21 * x658 + x327 * x715))
    result[20, 2] = numpy.sum(x148 * (x21 * x664 + x327 * x717 + x713))
    result[20, 3] = numpy.sum(x100 * (x21 * x667 + x327 * x718))
    result[20, 4] = numpy.sum(x148 * (x21 * x670 + x327 * x719 + x714))
    result[20, 5] = numpy.sum(x100 * (x21 * x673 + x327 * x720 + 2.0 * x716))
    return result


def _2center2el3d_53(ax, da, A, bx, db, B):
    """Cartesian (h|f) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((21, 10), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = ax ** (-1.0)
    x5 = bx ** (-1.0)
    x6 = -x2 - B[0]
    x7 = bx * x1
    x8 = ax * x7 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x9 = boys(6, x8)
    x10 = 17.4934183276249
    x11 = x0 ** (-1.5) * x10
    x12 = 2.0 * x11
    x13 = x12 * x5
    x14 = x13 * x9
    x15 = x0 ** (-0.5)
    x16 = boys(5, x8)
    x17 = 0.5 * x5
    x18 = x17 * (2.0 * x10 * x15 * x16 * x4 * x5 - x14)
    x19 = boys(7, x8)
    x20 = x6**2
    x21 = 2.0 * x4
    x22 = x10 * x15
    x23 = x22 * x5
    x24 = x21 * x23
    x25 = x20 * x24
    x26 = x19 * x25
    x27 = x18 + x26
    x28 = x6 * (x27 + x5 * (2.0 * x10 * x15 * x16 * x4 * x5 - x14))
    x29 = x28 * x3
    x30 = 0.5 / (ax + bx)
    x31 = x16 * x5
    x32 = x12 * x31
    x33 = boys(4, x8)
    x34 = x17 * (2.0 * x10 * x15 * x33 * x4 * x5 - x32)
    x35 = x24 * x9
    x36 = x20 * x35
    x37 = x34 + x36
    x38 = x30 * x37
    x39 = 3.0 * x38
    x40 = x29 + x39
    x41 = -2.0 * x10 * x15 * x33 * x4 * x5 * x6
    x42 = x37 * x6 - x5 * (x32 * x6 + x41)
    x43 = x3 * x42
    x44 = x13 * x33
    x45 = boys(3, x8)
    x46 = -2.0 * x10 * x15 * x4 * x45 * x5
    x47 = -x17 * (x44 + x46)
    x48 = x22 * x31
    x49 = x21 * x48
    x50 = x20 * x49
    x51 = x47 + x50
    x52 = x30 * x51
    x53 = 3.0 * x52
    x54 = x43 + x53
    x55 = x13 * x19
    x56 = x17 * (2.0 * x10 * x15 * x4 * x5 * x9 - x55)
    x57 = boys(8, x8)
    x58 = x25 * x57
    x59 = x6 * (x5 * (2.0 * x10 * x15 * x4 * x5 * x9 - x55) + x56 + x58)
    x60 = x3 * x59
    x61 = x27 * x30
    x62 = 3.0 * x61
    x63 = x28 * x7
    x64 = 0.5 * x4
    x65 = x64 * (x42 - x63)
    x66 = x27 * x3
    x67 = x6 * x9
    x68 = 4.0 * x30 * x4
    x69 = x23 * x68
    x70 = x67 * x69
    x71 = x66 + x70
    x72 = 3.0 * x30
    x73 = x37 * x7
    x74 = x64 * (x51 - x73)
    x75 = 2.0 * x30
    x76 = x4 * x75
    x77 = x48 * x76
    x78 = x3 * x6
    x79 = x35 * x78
    x80 = x77 + x79
    x81 = x3 * x71 + x74 + x75 * x80
    x82 = x3 * x37
    x83 = x48 * x68
    x84 = x6 * x83
    x85 = x82 + x84
    x86 = x3 * x51
    x87 = x33 * x69
    x88 = x6 * x87
    x89 = x86 + x88
    x90 = x11 * x21
    x91 = x16 * x90
    x92 = x6 * x91
    x93 = -x64 * (x41 + x92)
    x94 = x3 * x77
    x95 = x3 * x80 + x93 + x94
    x96 = x3 * x81 - x4 * (x7 * x85 - x89) + x75 * x95
    x97 = x6 * (x5 * (2.0 * x10 * x15 * x4 * x45 * x5 - x44) + x51)
    x98 = x64 * (-x42 * x7 + x97)
    x99 = x3 * x40 + x72 * x85 + x98
    x100 = x13 * x45
    x101 = boys(2, x8)
    x102 = x17 * (2.0 * x10 * x101 * x15 * x4 * x5 - x100)
    x103 = x21 * x33
    x104 = x103 * x23
    x105 = x104 * x20
    x106 = x102 + x105
    x107 = x6 * (x106 + x5 * (2.0 * x10 * x101 * x15 * x4 * x5 - x100))
    x108 = x64 * (x107 - x7 * x97)
    x109 = x108 + x3 * x54 + x72 * x89
    x110 = 1.5 * x4
    x111 = x17 * (2.0 * x10 * x15 * x4 * x5 * boys(1, x8) - x101 * x13)
    x112 = x24 * x45
    x113 = x111 + x112 * x20
    x114 = x113 * x30
    x115 = x106 * x30
    x116 = 3.0 * x115 + x3 * x97
    x117 = x64 * (-x106 * x7 + x113)
    x118 = x23 * x76
    x119 = x118 * x45
    x120 = x104 * x78 + x119
    x121 = x117 + x120 * x75 + x3 * x89
    x122 = x64 * (x106 - x51 * x7)
    x123 = x118 * x33
    x124 = x123 + x49 * x78
    x125 = x122 + x124 * x75 + x3 * x85
    x126 = x103 * x11
    x127 = -x64 * (x126 + x46)
    x128 = x3**2
    x129 = x128 * x49
    x130 = 0.179587122125167 * da * db * numpy.sqrt(ax**6.5) * numpy.sqrt(bx**4.5)
    x131 = 6.08167803818495 * x130
    x132 = -x1 * (ax * A[1] + bx * B[1])
    x133 = -x132 - B[1]
    x134 = x133 * x5 * (2.0 * x10 * x15 * x16 * x4 * x5 - x14)
    x135 = x133 * x26 + 0.5 * x134
    x136 = x135 * x3
    x137 = x133 * x6
    x138 = x69 * x9
    x139 = x137 * x138
    x140 = x136 + x139
    x141 = -2.0 * x10 * x133 * x15 * x33 * x4 * x5
    x142 = -x5 * (x133 * x32 + x141)
    x143 = x133 * x36 + 0.5 * x142
    x144 = x143 * x3
    x145 = x137 * x83
    x146 = x144 + x145
    x147 = x133 * x5 * (2.0 * x10 * x15 * x4 * x5 * x9 - x55)
    x148 = x133 * x58 + 0.5 * x147
    x149 = x148 * x3
    x150 = x137 * x19
    x151 = x150 * x69
    x152 = x135 * x7
    x153 = x64 * (x143 - x152)
    x154 = x118 * x9
    x155 = x133 * x154
    x156 = x19 * x24
    x157 = x156 * x78
    x158 = x133 * x157
    x159 = x155 + x158
    x160 = x9 * x90
    x161 = x137 * x160
    x162 = x64 * (2.0 * x10 * x133 * x15 * x16 * x4 * x5 * x6 - x161)
    x163 = x155 * x3
    x164 = x159 * x3 + x162 + x163
    x165 = x133 * x77
    x166 = x133 * x35
    x167 = x166 * x78
    x168 = x165 + x167
    x169 = x123 * x133
    x170 = x133 * x49
    x171 = x170 * x78
    x172 = x169 + x171
    x173 = x133 * x91
    x174 = -x64 * (x141 + x173)
    x175 = x128 * x166
    x176 = x174 + x175
    x177 = x164 * x3 + x176 * x30 - x4 * (x168 * x7 - x172)
    x178 = x133 * x5 * (2.0 * x10 * x15 * x4 * x45 * x5 - x44)
    x179 = x133 * x50 + 0.5 * x178
    x180 = x64 * (-x143 * x7 + x179)
    x181 = x140 * x3 + x168 * x75 + x180
    x182 = x133 * x5 * (2.0 * x10 * x101 * x15 * x4 * x5 - x100)
    x183 = x105 * x133 + 0.5 * x182
    x184 = x64 * (-x179 * x7 + x183)
    x185 = x146 * x3 + x172 * x75 + x184
    x186 = x137 * x87
    x187 = x179 * x3 + x186
    x188 = x45 * x69
    x189 = x137 * x188
    x190 = x64 * (2.0 * x10 * x133 * x15 * x4 * x45 * x5 * x6 - x126 * x137)
    x191 = x169 * x3 + x172 * x3 + x190
    x192 = x133 * x92
    x193 = x64 * (2.0 * x10 * x133 * x15 * x33 * x4 * x5 * x6 - x192)
    x194 = x133 * x94 + x168 * x3 + x193
    x195 = 13.5990455106491 * x130
    x196 = -x1 * (ax * A[2] + bx * B[2])
    x197 = -x196 - B[2]
    x198 = x197 * x5 * (2.0 * x10 * x15 * x16 * x4 * x5 - x14)
    x199 = 0.5 * x198
    x200 = x197 * x26 + x199
    x201 = x200 * x3
    x202 = x197 * x6
    x203 = x138 * x202
    x204 = x201 + x203
    x205 = -2.0 * x10 * x15 * x197 * x33 * x4 * x5
    x206 = -x5 * (x197 * x32 + x205)
    x207 = 0.5 * x206
    x208 = x197 * x36 + x207
    x209 = x208 * x3
    x210 = x202 * x83
    x211 = x209 + x210
    x212 = x197 * x5 * (2.0 * x10 * x15 * x4 * x5 * x9 - x55)
    x213 = 0.5 * x212
    x214 = x197 * x58 + x213
    x215 = x214 * x3
    x216 = x19 * x202
    x217 = x216 * x69
    x218 = x200 * x7
    x219 = x64 * (x208 - x218)
    x220 = x154 * x197
    x221 = x157 * x197
    x222 = x220 + x221
    x223 = x160 * x202
    x224 = x64 * (2.0 * x10 * x15 * x16 * x197 * x4 * x5 * x6 - x223)
    x225 = x220 * x3
    x226 = x222 * x3 + x224 + x225
    x227 = x197 * x77
    x228 = x197 * x35
    x229 = x228 * x78
    x230 = x227 + x229
    x231 = x123 * x197
    x232 = x197 * x49
    x233 = x232 * x78
    x234 = x231 + x233
    x235 = x197 * x91
    x236 = -x64 * (x205 + x235)
    x237 = x128 * x228
    x238 = x236 + x237
    x239 = x226 * x3 + x238 * x30 - x4 * (x230 * x7 - x234)
    x240 = x197 * x5 * (2.0 * x10 * x15 * x4 * x45 * x5 - x44)
    x241 = 0.5 * x240
    x242 = x197 * x50 + x241
    x243 = x64 * (-x208 * x7 + x242)
    x244 = x204 * x3 + x230 * x75 + x243
    x245 = x197 * x5 * (2.0 * x10 * x101 * x15 * x4 * x5 - x100)
    x246 = 0.5 * x245
    x247 = x105 * x197 + x246
    x248 = x64 * (-x242 * x7 + x247)
    x249 = x211 * x3 + x234 * x75 + x248
    x250 = x202 * x87
    x251 = x242 * x3 + x250
    x252 = x188 * x202
    x253 = x64 * (2.0 * x10 * x15 * x197 * x4 * x45 * x5 * x6 - x126 * x202)
    x254 = x231 * x3 + x234 * x3 + x253
    x255 = x197 * x64 * (2.0 * x10 * x15 * x33 * x4 * x5 * x6 - x92)
    x256 = x197 * x94 + x230 * x3 + x255
    x257 = x133**2
    x258 = x257 * x49 + x47
    x259 = x257 * x35 + x34
    x260 = x259 * x7
    x261 = x260 * x3
    x262 = x24 * x257
    x263 = x19 * x262
    x264 = x18 + x263
    x265 = x128 * x264
    x266 = x64 * (x258 - x260)
    x267 = x265 + x266
    x268 = x267 * x3 + x4 * (x258 * x3 - x261)
    x269 = x102 + x104 * x257
    x270 = x64 * (-x258 * x7 + x269)
    x271 = x128 * x259 + x270
    x272 = x111 + x112 * x257
    x273 = x64 * (-x269 * x7 + x272)
    x274 = x128 * x258 + x273
    x275 = x259 * x30
    x276 = x264 * x78
    x277 = x275 + x276
    x278 = x258 * x30
    x279 = x259 * x6
    x280 = x279 * x3
    x281 = x278 + x280
    x282 = x264 * x30
    x283 = x262 * x57
    x284 = x283 + x56
    x285 = x284 * x78
    x286 = x6 * x7
    x287 = x264 * x286
    x288 = x64 * (x259 * x6 - x287)
    x289 = x282 * x3
    x290 = x64 * (x258 * x6 - x279 * x7)
    x291 = x275 * x3
    x292 = x277 * x3 + x290 + x291
    x293 = x258 * x6
    x294 = x64 * (x269 * x6 - x293 * x7)
    x295 = x258 * x3
    x296 = x295 * x30
    x297 = x281 * x3 + x294 + x296
    x298 = x272 * x30
    x299 = x269 * x6
    x300 = x269 * x30
    x301 = x295 * x6 + x300
    x302 = x133 * x197
    x303 = x160 * x302
    x304 = x64 * (2.0 * x10 * x133 * x15 * x16 * x197 * x4 * x5 - x303)
    x305 = x156 * x302
    x306 = x128 * x305 + x304
    x307 = x3 * (x306 + x4 * (2.0 * x10 * x133 * x15 * x16 * x197 * x4 * x5 - x303))
    x308 = x197 * x64 * (2.0 * x10 * x133 * x15 * x33 * x4 * x5 - x173)
    x309 = x175 * x197 + x308
    x310 = x64 * (2.0 * x10 * x133 * x15 * x197 * x4 * x45 * x5 - x126 * x302)
    x311 = x129 * x302 + x310
    x312 = x154 * x302
    x313 = x157 * x302 + x312
    x314 = x167 * x197 + x302 * x77
    x315 = x19 * x302
    x316 = x118 * x315
    x317 = x24 * x57
    x318 = x150 * x197
    x319 = x64 * (2.0 * x10 * x133 * x15 * x197 * x4 * x5 * x6 * x9 - x318 * x90)
    x320 = x197 * x64 * (2.0 * x10 * x133 * x15 * x16 * x4 * x5 * x6 - x161)
    x321 = x3 * x312 + x3 * x313 + x320
    x322 = x197 * x64 * (2.0 * x10 * x133 * x15 * x33 * x4 * x5 * x6 - x192)
    x323 = x3 * x314 + x302 * x94 + x322
    x324 = x170 * x197
    x325 = x123 * x302 + x324 * x78
    x326 = x104 * x133
    x327 = 23.5542377588857 * x130
    x328 = x197**2
    x329 = x328 * x49 + x47
    x330 = x328 * x35 + x34
    x331 = x330 * x7
    x332 = x3 * x331
    x333 = x24 * x328
    x334 = x18 + x19 * x333
    x335 = x128 * x334
    x336 = x64 * (x329 - x331)
    x337 = x335 + x336
    x338 = x3 * x337 + x4 * (x3 * x329 - x332)
    x339 = x128 * x330
    x340 = x102 + x104 * x328
    x341 = x64 * (-x329 * x7 + x340)
    x342 = x339 + x341
    x343 = x111 + x112 * x328
    x344 = x64 * (-x340 * x7 + x343)
    x345 = x128 * x329 + x344
    x346 = x30 * x330
    x347 = x334 * x78
    x348 = x346 + x347
    x349 = x30 * x329
    x350 = x330 * x6
    x351 = x3 * x350
    x352 = x349 + x351
    x353 = x30 * x334
    x354 = x333 * x57 + x56
    x355 = x354 * x78
    x356 = x286 * x334
    x357 = x64 * (x330 * x6 - x356)
    x358 = x3 * x353
    x359 = x64 * (x329 * x6 - x350 * x7)
    x360 = x3 * x346
    x361 = x3 * x348 + x359 + x360
    x362 = x329 * x6
    x363 = x64 * (x340 * x6 - x362 * x7)
    x364 = x3 * x329
    x365 = x30 * x364
    x366 = x3 * x352 + x363 + x365
    x367 = x30 * x343
    x368 = x340 * x6
    x369 = x30 * x340
    x370 = x364 * x6 + x369
    x371 = x133 * x259 + x142
    x372 = x133 * x264 + x134
    x373 = x372 * x7
    x374 = x3 * x373
    x375 = x133 * x284 + x147
    x376 = x128 * x375
    x377 = x64 * (x371 - x373)
    x378 = x128 * x372
    x379 = x133 * x258 + x178
    x380 = x64 * (-x371 * x7 + x379)
    x381 = x378 + x380
    x382 = x128 * x371
    x383 = x133 * x269 + x182
    x384 = x64 * (-x379 * x7 + x383)
    x385 = x382 + x384
    x386 = x207 + x228 * x257
    x387 = x197 * x263 + x199
    x388 = x387 * x7
    x389 = x197 * x283 + x213
    x390 = x64 * (x386 - x388)
    x391 = x232 * x257 + x241
    x392 = x64 * (-x386 * x7 + x391)
    x393 = x128 * x387 + x392
    x394 = x104 * x197
    x395 = x246 + x257 * x394
    x396 = x64 * (-x391 * x7 + x395)
    x397 = x128 * x386 + x396
    x398 = x133 * x334
    x399 = x398 * x7
    x400 = x64 * (x133 * x330 - x399)
    x401 = x133 * x354
    x402 = x133 * x330
    x403 = x64 * (x133 * x329 - x402 * x7)
    x404 = x133 * x335 + x403
    x405 = x133 * x329
    x406 = x64 * (x133 * x340 - x405 * x7)
    x407 = x133 * x339 + x406
    x408 = x197 * x330 + x206
    x409 = x197 * x334 + x198
    x410 = x409 * x7
    x411 = x3 * x410
    x412 = x197 * x354 + x212
    x413 = x128 * x412
    x414 = x64 * (x408 - x410)
    x415 = x128 * x409
    x416 = x197 * x329 + x240
    x417 = x64 * (-x408 * x7 + x416)
    x418 = x415 + x417
    x419 = x128 * x408
    x420 = x197 * x340 + x245
    x421 = x64 * (-x416 * x7 + x420)
    x422 = x419 + x421
    x423 = -x132 - A[1]
    x424 = x29 * x423
    x425 = x39 * x423
    x426 = x424 + x425
    x427 = x423 * x43
    x428 = x423 * x51
    x429 = x427 + x428 * x72
    x430 = x423 * x60
    x431 = x423 * x62
    x432 = x423 * x63
    x433 = x4 * (x42 * x423 - x432)
    x434 = x423 * x66
    x435 = x423 * x70
    x436 = x434 + x435
    x437 = x423 * x73
    x438 = x4 * (x423 * x51 - x437)
    x439 = x423 * x77
    x440 = x423 * x79
    x441 = x439 + x440
    x442 = x3 * x436 + 0.5 * x438 + x441 * x75
    x443 = x423 * (x82 + x84)
    x444 = x423 * (x86 + x88)
    x445 = x4 * x423 * (2.0 * x10 * x15 * x33 * x4 * x5 * x6 - x92)
    x446 = x4 * x423 * (x107 - x7 * x97)
    x447 = x4 * x423 * (-x42 * x7 + x97)
    x448 = 18.2450341145548 * x130
    x449 = x135 * x423
    x450 = x38 + x449
    x451 = x6 * x77
    x452 = x35 * x6
    x453 = x133 * x423 * x452
    x454 = x451 + x453
    x455 = x454 * x75
    x456 = x3 * x450 + x455
    x457 = x143 * x423
    x458 = x457 + x52
    x459 = x123 * x6
    x460 = x49 * x6
    x461 = x133 * x460
    x462 = x423 * x461
    x463 = x459 + x462
    x464 = x463 * x75
    x465 = x3 * x458 + x464
    x466 = x148 * x423
    x467 = x466 + x61
    x468 = x118 * x67
    x469 = x150 * x24
    x470 = x423 * x469
    x471 = x468 + x470
    x472 = x471 * x75
    x473 = x4 * (-x450 * x7 + x458)
    x474 = x166 * x423
    x475 = x474 + x77
    x476 = x30 * x475
    x477 = x3 * x471 + x476
    x478 = x4 * (-x454 * x7 + x463)
    x479 = x3 * x476 + x3 * x477 + 0.5 * x478
    x480 = x123 + x170 * x423
    x481 = x30 * x480
    x482 = x3 * x454 + x481
    x483 = x119 + x326 * x423
    x484 = x30 * x483
    x485 = x3 * x463 + x484
    x486 = x4 * (-x480 * x7 + x483)
    x487 = x115 + x179 * x423
    x488 = x4 * (x114 + x183 * x423 - x487 * x7)
    x489 = x4 * (-x458 * x7 + x487)
    x490 = 40.7971365319473 * x130
    x491 = x423 * (x201 + x203)
    x492 = x423 * (x209 + x210)
    x493 = x4 * x423 * (x208 - x218)
    x494 = x220 * x423
    x495 = x221 * x423 + x494
    x496 = x4 * x423 * (2.0 * x10 * x15 * x16 * x197 * x4 * x5 * x6 - x223)
    x497 = x225 * x423 + x3 * x495 + 0.5 * x496
    x498 = x227 * x423
    x499 = x229 * x423 + x498
    x500 = x231 * x423
    x501 = x233 * x423 + x500
    x502 = x4 * x423 * (2.0 * x10 * x15 * x197 * x33 * x4 * x5 - x235)
    x503 = x4 * x423 * (-x242 * x7 + x247)
    x504 = x4 * x423 * (-x208 * x7 + x242)
    x505 = x258 * x423
    x506 = x133 * x87
    x507 = x505 + x506
    x508 = x259 * x423
    x509 = x133 * x83
    x510 = x508 + x509
    x511 = x510 * x7
    x512 = x264 * x423
    x513 = x133 * x138
    x514 = x512 + x513
    x515 = x4 * (x507 - x511)
    x516 = x128 * x514 + 0.5 * x515
    x517 = x30 * x510
    x518 = x512 * x6
    x519 = x139 + x518
    x520 = x3 * x519 + x517
    x521 = x30 * x507
    x522 = x279 * x423
    x523 = x145 + x522
    x524 = x3 * x523 + x521
    x525 = x30 * x514
    x526 = x284 * x6
    x527 = x423 * x526
    x528 = x151 + x527
    x529 = x4 * (-x519 * x7 + x523)
    x530 = x186 + x293 * x423
    x531 = x4 * (x189 + x299 * x423 - x530 * x7)
    x532 = x4 * (-x523 * x7 + x530)
    x533 = x231 + x324 * x423
    x534 = x197 * x474 + x227
    x535 = x534 * x7
    x536 = x220 + x305 * x423
    x537 = x4 * (x533 - x535)
    x538 = x128 * x536 + 0.5 * x537
    x539 = x30 * x534
    x540 = x154 * x202
    x541 = x24 * x318
    x542 = x423 * x541 + x540
    x543 = x3 * x542 + x539
    x544 = x30 * x533
    x545 = x302 * x452
    x546 = x202 * x77 + x423 * x545
    x547 = x3 * x546 + x544
    x548 = x30 * x536
    x549 = x118 * x216
    x550 = x137 * x317
    x551 = x197 * x423 * x550 + x549
    x552 = x4 * (-x542 * x7 + x546)
    x553 = x302 * x460
    x554 = x123 * x202 + x423 * x553
    x555 = x104 * x6
    x556 = x4 * (x119 * x202 + x302 * x423 * x555 - x554 * x7)
    x557 = x4 * (-x546 * x7 + x554)
    x558 = 70.6627132766571 * x130
    x559 = x4 * x423 * (x329 - x331)
    x560 = x335 * x423 + 0.5 * x559
    x561 = x346 * x423
    x562 = x347 * x423 + x561
    x563 = x349 * x423
    x564 = x351 * x423 + x563
    x565 = x353 * x423
    x566 = x4 * x423 * (x330 * x6 - x356)
    x567 = x4 * x423 * (x340 * x6 - x362 * x7)
    x568 = x4 * x423 * (x329 * x6 - x350 * x7)
    x569 = x371 * x423
    x570 = 3.0 * x278
    x571 = x569 + x570
    x572 = x372 * x423
    x573 = 3.0 * x275
    x574 = x572 + x573
    x575 = x574 * x7
    x576 = x375 * x423
    x577 = 3.0 * x282
    x578 = x576 + x577
    x579 = x4 * (x571 - x575)
    x580 = 3.0 * x300 + x379 * x423
    x581 = x4 * (3.0 * x298 + x383 * x423 - x580 * x7)
    x582 = x4 * (-x571 * x7 + x580)
    x583 = x302 * x83
    x584 = x386 * x423 + x583
    x585 = x138 * x302
    x586 = x387 * x423 + x585
    x587 = x586 * x7
    x588 = x315 * x69
    x589 = x389 * x423 + x588
    x590 = x4 * (x584 - x587)
    x591 = x302 * x87
    x592 = x391 * x423 + x591
    x593 = x188 * x302
    x594 = x4 * (x395 * x423 - x592 * x7 + x593)
    x595 = x4 * (-x584 * x7 + x592)
    x596 = x349 + x402 * x423
    x597 = x346 + x398 * x423
    x598 = x597 * x7
    x599 = x353 + x401 * x423
    x600 = x4 * (x596 - x598)
    x601 = x133 * x340
    x602 = x369 + x405 * x423
    x603 = x4 * (x367 + x423 * x601 - x602 * x7)
    x604 = x4 * (-x596 * x7 + x602)
    x605 = x4 * x423 * (x408 - x410)
    x606 = x4 * x423 * (-x416 * x7 + x420)
    x607 = x4 * x423 * (-x408 * x7 + x416)
    x608 = -x196 - A[2]
    x609 = x608 * (x29 + x39)
    x610 = x608 * (x43 + x53)
    x611 = x4 * x608 * (x42 - x63)
    x612 = 0.5 * x611
    x613 = x608 * (x66 + x70)
    x614 = x4 * x608 * (x51 - x73)
    x615 = 0.5 * x614
    x616 = x608 * x77
    x617 = x608 * x79 + x616
    x618 = x3 * x613 + x615 + x617 * x75
    x619 = x608 * (x82 + x84)
    x620 = x608 * (x86 + x88)
    x621 = x4 * x608 * (2.0 * x10 * x15 * x33 * x4 * x5 * x6 - x92)
    x622 = 0.5 * x621
    x623 = x4 * x608 * (x107 - x7 * x97)
    x624 = 0.5 * x623
    x625 = x42 * x608
    x626 = x4 * (x608 * x97 - x625 * x7)
    x627 = 0.5 * x626
    x628 = x139 * x608
    x629 = x136 * x608 + x628
    x630 = x145 * x608
    x631 = x144 * x608 + x630
    x632 = x151 * x608
    x633 = x4 * x608 * (x143 - x152)
    x634 = 0.5 * x633
    x635 = x155 * x608
    x636 = x158 * x608 + x635
    x637 = x4 * x608 * (2.0 * x10 * x133 * x15 * x16 * x4 * x5 * x6 - x161)
    x638 = 0.5 * x637
    x639 = x163 * x608 + x3 * x636 + x638
    x640 = x165 * x608
    x641 = x167 * x608 + x640
    x642 = x169 * x608
    x643 = x171 * x608 + x642
    x644 = x4 * x608 * (2.0 * x10 * x133 * x15 * x33 * x4 * x5 - x173)
    x645 = 0.5 * x644
    x646 = x4 * x608 * (-x179 * x7 + x183)
    x647 = 0.5 * x646
    x648 = x4 * x608 * (-x143 * x7 + x179)
    x649 = 0.5 * x648
    x650 = x200 * x608 + x38
    x651 = x3 * x650
    x652 = x197 * x608
    x653 = x451 + x452 * x652
    x654 = x30 * x653
    x655 = 2.0 * x654
    x656 = x651 + x655
    x657 = x208 * x608 + x52
    x658 = x3 * x657
    x659 = x197 * x460
    x660 = x459 + x608 * x659
    x661 = x30 * x660
    x662 = 2.0 * x661
    x663 = x658 + x662
    x664 = x214 * x608 + x61
    x665 = x3 * x664
    x666 = x156 * x202
    x667 = x468 + x608 * x666
    x668 = x30 * x667
    x669 = 2.0 * x668
    x670 = x650 * x7
    x671 = x4 * (x657 - x670)
    x672 = 0.5 * x671
    x673 = x228 * x608 + x77
    x674 = x30 * x673
    x675 = x3 * x667
    x676 = x674 + x675
    x677 = x653 * x7
    x678 = x4 * (x660 - x677)
    x679 = 0.5 * x678
    x680 = x3 * x674
    x681 = x3 * x676 + x679 + x680
    x682 = x123 + x232 * x608
    x683 = x30 * x682
    x684 = x3 * x653 + x683
    x685 = x119 + x394 * x608
    x686 = x30 * x685
    x687 = x3 * x660 + x686
    x688 = x4 * (-x682 * x7 + x685)
    x689 = 0.5 * x688
    x690 = x115 + x242 * x608
    x691 = x4 * (x114 + x247 * x608 - x690 * x7)
    x692 = 0.5 * x691
    x693 = x4 * (-x657 * x7 + x690)
    x694 = 0.5 * x693
    x695 = x4 * x608 * (x258 - x260)
    x696 = 0.5 * x695
    x697 = x265 * x608 + x696
    x698 = x275 * x608
    x699 = x276 * x608 + x698
    x700 = x278 * x608
    x701 = x280 * x608 + x700
    x702 = x282 * x608
    x703 = x4 * x608 * (x259 * x6 - x287)
    x704 = 0.5 * x703
    x705 = x4 * x608 * (x269 * x6 - x293 * x7)
    x706 = 0.5 * x705
    x707 = x4 * x608 * (x258 * x6 - x279 * x7)
    x708 = 0.5 * x707
    x709 = x169 + x324 * x608
    x710 = x165 + x166 * x652
    x711 = x7 * x710
    x712 = x155 + x305 * x608
    x713 = x4 * (x709 - x711)
    x714 = 0.5 * x713
    x715 = x128 * x712 + x714
    x716 = x30 * x710
    x717 = x137 * x154
    x718 = x541 * x608 + x717
    x719 = x3 * x718 + x716
    x720 = x30 * x709
    x721 = x137 * x77 + x545 * x608
    x722 = x3 * x721 + x720
    x723 = x30 * x712
    x724 = x118 * x150
    x725 = x550 * x652 + x724
    x726 = x4 * (-x7 * x718 + x721)
    x727 = 0.5 * x726
    x728 = x123 * x137 + x553 * x608
    x729 = x4 * (x119 * x137 + x302 * x555 * x608 - x7 * x728)
    x730 = 0.5 * x729
    x731 = x4 * (-x7 * x721 + x728)
    x732 = 0.5 * x731
    x733 = x197 * x87 + x329 * x608
    x734 = x197 * x83 + x330 * x608
    x735 = x7 * x734
    x736 = x334 * x608
    x737 = x138 * x197 + x736
    x738 = x128 * x737
    x739 = x4 * (x733 - x735)
    x740 = 0.5 * x739
    x741 = x738 + x740
    x742 = x30 * x734
    x743 = x203 + x6 * x736
    x744 = x3 * x743
    x745 = x742 + x744
    x746 = x30 * x733
    x747 = x210 + x350 * x608
    x748 = x3 * x747
    x749 = x746 + x748
    x750 = x30 * x737
    x751 = x354 * x6
    x752 = x217 + x608 * x751
    x753 = x3 * x752
    x754 = x7 * x743
    x755 = x4 * (x747 - x754)
    x756 = 0.5 * x755
    x757 = x3 * x750
    x758 = x250 + x362 * x608
    x759 = x4 * (x252 + x368 * x608 - x7 * x758)
    x760 = 0.5 * x759
    x761 = x4 * (-x7 * x747 + x758)
    x762 = 0.5 * x761
    x763 = x4 * x608 * (x371 - x373)
    x764 = 0.5 * x763
    x765 = x4 * x608 * (-x379 * x7 + x383)
    x766 = 0.5 * x765
    x767 = x4 * x608 * (-x371 * x7 + x379)
    x768 = 0.5 * x767
    x769 = x278 + x386 * x608
    x770 = x275 + x387 * x608
    x771 = x7 * x770
    x772 = x282 + x389 * x608
    x773 = x4 * (x769 - x771)
    x774 = 0.5 * x773
    x775 = x300 + x391 * x608
    x776 = x4 * (x298 + x395 * x608 - x7 * x775)
    x777 = 0.5 * x776
    x778 = x4 * (-x7 * x769 + x775)
    x779 = 0.5 * x778
    x780 = x402 * x608 + x583
    x781 = x133 * x736 + x585
    x782 = x7 * x781
    x783 = x401 * x608 + x588
    x784 = x4 * (x780 - x782)
    x785 = 0.5 * x784
    x786 = x405 * x608 + x591
    x787 = x4 * (x593 + x601 * x608 - x7 * x786)
    x788 = 0.5 * x787
    x789 = x4 * (-x7 * x780 + x786)
    x790 = 0.5 * x789
    x791 = 3.0 * x349 + x408 * x608
    x792 = 3.0 * x346 + x409 * x608
    x793 = x7 * x792
    x794 = x3 * x793
    x795 = 3.0 * x353 + x412 * x608
    x796 = x128 * x795
    x797 = x4 * (x791 - x793)
    x798 = 0.5 * x797
    x799 = 3.0 * x369 + x416 * x608
    x800 = x4 * (3.0 * x367 + x420 * x608 - x7 * x799)
    x801 = 0.5 * x800
    x802 = x4 * (-x7 * x791 + x799)
    x803 = 0.5 * x802
    x804 = x423**2
    x805 = x108 + x42 * x804
    x806 = x117 + x51 * x804
    x807 = x30 * x806
    x808 = x28 * x804
    x809 = x808 + x98
    x810 = x122 + x37 * x804
    x811 = x30 * x810
    x812 = x59 * x804
    x813 = x65 + x812
    x814 = x27 * x804
    x815 = x74 + x814
    x816 = x30 * x815
    x817 = -x7 * x809 + x805
    x818 = x452 * x804
    x819 = x818 + x93
    x820 = x30 * x819
    x821 = x3 * x815 + 2.0 * x820
    x822 = -x7 * x810 + x806
    x823 = x30 * (x127 + x49 * x804)
    x824 = 27.869749962333 * x130
    x825 = x30 * x428
    x826 = x184 + x423 * x458 + x825
    x827 = x190 + x423 * x459 + x423 * x463
    x828 = x75 * x827
    x829 = x38 * x423
    x830 = x180 + x423 * x450 + x829
    x831 = x193 + x423 * x454 + x439 * x6
    x832 = x75 * x831
    x833 = x423 * x61
    x834 = x153 + x423 * x467 + x833
    x835 = x423 * x468
    x836 = x162 + x423 * x471 + x835
    x837 = x75 * x836
    x838 = -x7 * x830 + x826
    x839 = x174 + x423 * x475 + x439
    x840 = x30 * x839
    x841 = x3 * x836 + x840
    x842 = -x7 * x831 + x827
    x843 = 62.3186554316989 * x130
    x844 = x208 * x804 + x248
    x845 = x253 + x659 * x804
    x846 = x30 * x845
    x847 = x200 * x804 + x243
    x848 = x197 * x818 + x255
    x849 = x30 * x848
    x850 = x214 * x804 + x219
    x851 = x224 + x666 * x804
    x852 = x30 * x851
    x853 = -x7 * x847 + x844
    x854 = x228 * x804 + x236
    x855 = x30 * x854
    x856 = x3 * x851 + x855
    x857 = -x7 * x848 + x845
    x858 = x273 + x423 * x507 + 2.0 * x484
    x859 = x30 * x858
    x860 = x294 + x423 * x523 + x464
    x861 = x270 + x423 * x510 + 2.0 * x481
    x862 = x30 * x861
    x863 = x290 + x423 * x519 + x455
    x864 = x266 + x423 * x514 + 2.0 * x476
    x865 = -x7 * x861 + x858
    x866 = x30 * x864
    x867 = x288 + x423 * x528 + x472
    x868 = -x7 * x863 + x860
    x869 = x310 + x423 * x533 + x500
    x870 = x30 * x869
    x871 = x202 * x439 + x322 + x423 * x546
    x872 = x308 + x423 * x534 + x498
    x873 = x30 * x872
    x874 = x320 + x423 * x540 + x423 * x542
    x875 = x304 + x423 * x536 + x494
    x876 = -x7 * x872 + x869
    x877 = x30 * x875
    x878 = x319 + x423 * x549 + x423 * x551
    x879 = -x7 * x874 + x871
    x880 = 107.939077467081 * x130
    x881 = x329 * x804 + x344
    x882 = x30 * x881
    x883 = x350 * x804 + x363
    x884 = x330 * x804 + x341
    x885 = x30 * x884
    x886 = x334 * x804
    x887 = x359 + x6 * x886
    x888 = x336 + x886
    x889 = -x7 * x884 + x881
    x890 = x30 * x888
    x891 = x357 + x751 * x804
    x892 = -x7 * x887 + x883
    x893 = x384 + x423 * x571 + 3.0 * x521
    x894 = x380 + x423 * x574 + 3.0 * x517
    x895 = x7 * x894
    x896 = x377 + x423 * x578 + 3.0 * x525
    x897 = x893 - x895
    x898 = x396 + x423 * x584 + 2.0 * x544
    x899 = x392 + x423 * x586 + 2.0 * x539
    x900 = x7 * x899
    x901 = x390 + x423 * x589 + 2.0 * x548
    x902 = x898 - x900
    x903 = x406 + x423 * x596 + x563
    x904 = x403 + x423 * x597 + x561
    x905 = x7 * x904
    x906 = x400 + x423 * x599 + x565
    x907 = x903 - x905
    x908 = x408 * x804 + x421
    x909 = x409 * x804 + x417
    x910 = x7 * x909
    x911 = x412 * x804 + x414
    x912 = x908 - x910
    x913 = x4 * x608 * (x42 * x423 - x432)
    x914 = x608 * (x434 + x435)
    x915 = x4 * x608 * (x423 * x51 - x437)
    x916 = x439 * x608
    x917 = 48.2718229290016 * x130
    x918 = x52 * x608
    x919 = x457 * x608 + x918
    x920 = x459 * x608
    x921 = x462 * x608 + x920
    x922 = x75 * x921
    x923 = x38 * x608
    x924 = x449 * x608 + x923
    x925 = x451 * x608
    x926 = x453 * x608 + x925
    x927 = x75 * x926
    x928 = x608 * x61
    x929 = x466 * x608 + x928
    x930 = x468 * x608
    x931 = x470 * x608 + x930
    x932 = x75 * x931
    x933 = x4 * (-x7 * x924 + x919)
    x934 = x474 * x608 + x616
    x935 = x30 * x934
    x936 = x3 * x931 + x935
    x937 = x4 * (-x7 * x926 + x921)
    x938 = x423 * x660
    x939 = x4 * x423 * (x657 - x670)
    x940 = x423 * x674
    x941 = x423 * x675 + x940
    x942 = x4 * x423 * (x660 - x677)
    x943 = x608 * (x505 + x506)
    x944 = x30 * x943
    x945 = x522 * x608 + x630
    x946 = x608 * (x508 + x509)
    x947 = x30 * x946
    x948 = x518 * x608 + x628
    x949 = x608 * (x512 + x513)
    x950 = x4 * (-x7 * x946 + x943)
    x951 = x30 * x949
    x952 = x527 * x608 + x632
    x953 = x4 * (-x7 * x948 + x945)
    x954 = x423 * x709 + x686
    x955 = x30 * x954
    x956 = x423 * x721 + x661
    x957 = x423 * x710 + x683
    x958 = x30 * x957
    x959 = x423 * x718 + x654
    x960 = x423 * x712 + x674
    x961 = x4 * (-x7 * x957 + x954)
    x962 = x30 * x960
    x963 = x423 * x725 + x668
    x964 = x4 * (-x7 * x959 + x956)
    x965 = 186.955966295097 * x130
    x966 = x423 * x746
    x967 = x423 * x742
    x968 = x4 * x423 * (x733 - x735)
    x969 = x423 * x750
    x970 = x4 * x423 * (x747 - x754)
    x971 = x608 * (x569 + x570)
    x972 = x608 * (x572 + x573)
    x973 = x7 * x972
    x974 = x608 * (x576 + x577)
    x975 = x4 * (x971 - x973)
    x976 = 2.0 * x720
    x977 = x423 * x769 + x976
    x978 = 2.0 * x716
    x979 = x423 * x770 + x978
    x980 = x7 * x979
    x981 = 2.0 * x723
    x982 = x423 * x772 + x981
    x983 = x4 * (x977 - x980)
    x984 = x423 * x780 + x746
    x985 = x423 * x781 + x742
    x986 = x7 * x985
    x987 = x423 * x783 + x750
    x988 = x4 * (x984 - x986)
    x989 = x4 * x423 * (x791 - x793)
    x990 = x608**2
    x991 = x108 + x42 * x990
    x992 = x117 + x51 * x990
    x993 = x30 * x992
    x994 = x28 * x990 + x98
    x995 = x122 + x37 * x990
    x996 = x30 * x995
    x997 = x59 * x990 + x65
    x998 = x3 * x997
    x999 = x27 * x990 + x74
    x1000 = x30 * x999
    x1001 = 3.0 * x1000
    x1002 = x7 * x994
    x1003 = -x1002 + x991
    x1004 = x1003 * x64
    x1005 = x3 * x999
    x1006 = x452 * x990
    x1007 = x1006 + x93
    x1008 = x1007 * x30
    x1009 = 2.0 * x1008
    x1010 = x1005 + x1009
    x1011 = -x7 * x995 + x992
    x1012 = x1011 * x64
    x1013 = x30 * (x127 + x49 * x990)
    x1014 = x143 * x990 + x184
    x1015 = x190 + x461 * x990
    x1016 = x1015 * x30
    x1017 = 2.0 * x1016
    x1018 = x135 * x990 + x180
    x1019 = x1006 * x133 + x193
    x1020 = x1019 * x30
    x1021 = 2.0 * x1020
    x1022 = x148 * x990 + x153
    x1023 = x162 + x469 * x990
    x1024 = x1023 * x30
    x1025 = 2.0 * x1024
    x1026 = x1014 - x1018 * x7
    x1027 = x1026 * x64
    x1028 = x166 * x990 + x174
    x1029 = x1028 * x30
    x1030 = x1023 * x3 + x1029
    x1031 = x1015 - x1019 * x7
    x1032 = x1031 * x64
    x1033 = x248 + x608 * x657 + x918
    x1034 = x253 + x608 * x660 + x920
    x1035 = x1034 * x30
    x1036 = 2.0 * x1035
    x1037 = x243 + x608 * x650 + x923
    x1038 = x255 + x608 * x653 + x925
    x1039 = x1038 * x30
    x1040 = 2.0 * x1039
    x1041 = x219 + x608 * x664 + x928
    x1042 = x1041 * x3
    x1043 = x224 + x608 * x667 + x930
    x1044 = x1043 * x30
    x1045 = 2.0 * x1044
    x1046 = x1037 * x7
    x1047 = x1033 - x1046
    x1048 = x1047 * x64
    x1049 = x236 + x608 * x673 + x616
    x1050 = x1049 * x30
    x1051 = x1043 * x3
    x1052 = x1050 + x1051
    x1053 = x1034 - x1038 * x7
    x1054 = x1053 * x64
    x1055 = x258 * x990 + x273
    x1056 = x1055 * x30
    x1057 = x279 * x990 + x294
    x1058 = x259 * x990 + x270
    x1059 = x1058 * x30
    x1060 = x264 * x990
    x1061 = x1060 * x6 + x290
    x1062 = x1060 + x266
    x1063 = x1055 - x1058 * x7
    x1064 = x1063 * x64
    x1065 = x1062 * x30
    x1066 = x288 + x526 * x990
    x1067 = x1057 - x1061 * x7
    x1068 = x1067 * x64
    x1069 = x310 + x608 * x709 + x642
    x1070 = x1069 * x30
    x1071 = x137 * x616 + x322 + x608 * x721
    x1072 = x308 + x608 * x710 + x640
    x1073 = x1072 * x30
    x1074 = x320 + x608 * x717 + x608 * x718
    x1075 = x304 + x608 * x712 + x635
    x1076 = x1069 - x1072 * x7
    x1077 = x1076 * x64
    x1078 = x1075 * x30
    x1079 = x319 + x608 * x724 + x608 * x725
    x1080 = x1071 - x1074 * x7
    x1081 = x1080 * x64
    x1082 = x344 + x608 * x733 + 2.0 * x686
    x1083 = x1082 * x30
    x1084 = x363 + x608 * x747 + x662
    x1085 = x341 + x608 * x734 + 2.0 * x683
    x1086 = x1085 * x30
    x1087 = x359 + x608 * x743 + x655
    x1088 = x336 + x608 * x737 + 2.0 * x674
    x1089 = x1082 - x1085 * x7
    x1090 = x1089 * x64
    x1091 = x1088 * x30
    x1092 = x357 + x608 * x752 + x669
    x1093 = x1092 * x3
    x1094 = x1087 * x7
    x1095 = x1084 - x1094
    x1096 = x1095 * x64
    x1097 = x1091 * x3
    x1098 = x371 * x990 + x384
    x1099 = x372 * x990 + x380
    x1100 = x1099 * x7
    x1101 = x375 * x990 + x377
    x1102 = x1098 - x1100
    x1103 = x1102 * x64
    x1104 = x396 + x608 * x769 + x700
    x1105 = x392 + x608 * x770 + x698
    x1106 = x1105 * x7
    x1107 = x390 + x608 * x772 + x702
    x1108 = x1104 - x1106
    x1109 = x1108 * x64
    x1110 = x406 + x608 * x780 + x976
    x1111 = x403 + x608 * x781 + x978
    x1112 = x1111 * x7
    x1113 = x400 + x608 * x783 + x981
    x1114 = x1110 - x1112
    x1115 = x1114 * x64
    x1116 = x421 + x608 * x791 + 3.0 * x746
    x1117 = x417 + x608 * x792 + 3.0 * x742
    x1118 = x1117 * x7
    x1119 = x414 + x608 * x795 + 3.0 * x750
    x1120 = x1119 * x128
    x1121 = x1116 - x1118
    x1122 = x1121 * x64
    x1123 = x423 * x813 + x433
    x1124 = x423 * x815 + x438
    x1125 = x1124 * x30
    x1126 = x423 * x805 + x446 - x7 * (x423 * x809 + x447)
    x1127 = x30 * (x423 * x819 + x445)
    x1128 = 27.869749962333 * x130
    x1129 = x423 * x834 + x473 + x816
    x1130 = x423 * x836 + x478 + x820
    x1131 = x1130 * x75
    x1132 = x423 * x826 + x488 - x7 * (x423 * x830 + x489 + x811) + x807
    x1133 = x30 * (x423 * x839 + x486 + x823)
    x1134 = 62.3186554316989 * x130
    x1135 = x423 * x850 + x493
    x1136 = x423 * x851 + x496
    x1137 = x1136 * x30
    x1138 = x423 * x844 + x503 - x7 * (x423 * x847 + x504)
    x1139 = x30 * (x423 * x854 + x502)
    x1140 = x423 * x864 + x515 + 2.0 * x840
    x1141 = x1140 * x30
    x1142 = x423 * x867 + x529 + x837
    x1143 = x423 * x860 + x531 - x7 * (x423 * x863 + x532 + x832) + x828
    x1144 = x423 * x875 + x537 + x855
    x1145 = x1144 * x30
    x1146 = x423 * x878 + x552 + x852
    x1147 = x423 * x871 + x556 - x7 * (x423 * x874 + x557 + x849) + x846
    x1148 = 107.939077467081 * x130
    x1149 = x423 * x888 + x559
    x1150 = x1149 * x30
    x1151 = x423 * x891 + x566
    x1152 = x423 * x883 + x567 - x7 * (x423 * x887 + x568)
    x1153 = x423 * x896 + x579 + 3.0 * x866
    x1154 = x423 * x893 + x581 - x7 * (x423 * x894 + x582 + 3.0 * x862) + 3.0 * x859
    x1155 = x423 * x901 + x590 + 2.0 * x877
    x1156 = x423 * x898 + x594 - x7 * (x423 * x899 + x595 + 2.0 * x873) + 2.0 * x870
    x1157 = x423 * x906 + x600 + x890
    x1158 = x423 * x903 + x603 - x7 * (x423 * x904 + x604 + x885) + x882
    x1159 = x423 * x911 + x605
    x1160 = x423 * x908 + x606 - x7 * (x423 * x909 + x607)
    x1161 = x608 * x812 + x612
    x1162 = x608 * x814 + x615
    x1163 = x1162 * x30
    x1164 = x624 + x625 * x804 - x7 * (x608 * x808 + x627)
    x1165 = x30 * (x608 * x818 + x622)
    x1166 = x423 * x929 + x608 * x833 + x634
    x1167 = x423 * x931 + x608 * x835 + x638
    x1168 = x1167 * x75
    x1169 = x423 * x919 + x608 * x825 + x647 - x7 * (x423 * x924 + x608 * x829 + x649)
    x1170 = x30 * (x423 * x934 + x645 + x916)
    x1171 = 139.348749811665 * x130
    x1172 = x664 * x804 + x672
    x1173 = x667 * x804 + x679
    x1174 = x1173 * x30
    x1175 = x657 * x804 + x692 - x7 * (x650 * x804 + x694)
    x1176 = x30 * (x673 * x804 + x689)
    x1177 = x423 * x949 + x696 + 2.0 * x935
    x1178 = x1177 * x30
    x1179 = x423 * x952 + x704 + x932
    x1180 = x423 * x945 - x7 * (x423 * x948 + x708 + x927) + x706 + x922
    x1181 = x423 * x960 + x714 + x940
    x1182 = x1181 * x30
    x1183 = x423 * x668 + x423 * x963 + x727
    x1184 = x30 * x938 + x423 * x956 - x7 * (x423 * x654 + x423 * x959 + x732) + x730
    x1185 = 241.359114645008 * x130
    x1186 = x737 * x804 + x740
    x1187 = x1186 * x30
    x1188 = x752 * x804 + x756
    x1189 = -x7 * (x743 * x804 + x762) + x747 * x804 + x760
    x1190 = x423 * x974 + x764 + 3.0 * x951
    x1191 = x423 * x971 - x7 * (x423 * x972 + x768 + 3.0 * x947) + x766 + 3.0 * x944
    x1192 = x423 * x982 + x774 + 2.0 * x962
    x1193 = x423 * x977 - x7 * (x423 * x979 + x779 + 2.0 * x958) + x777 + 2.0 * x955
    x1194 = x423 * x987 + x785 + x969
    x1195 = x423 * x984 - x7 * (x423 * x985 + x790 + x967) + x788 + x966
    x1196 = x795 * x804 + x798
    x1197 = -x7 * (x792 * x804 + x803) + x791 * x804 + x801
    x1198 = x4 * x423 * (-x1002 + x991)
    x1199 = x1000 + x1022 * x423
    x1200 = x1008 + x1023 * x423
    x1201 = x1200 * x75
    x1202 = x4 * (x1014 * x423 - x7 * (x1018 * x423 + x996) + x993)
    x1203 = x30 * (x1013 + x1028 * x423)
    x1204 = x4 * x423 * (x1033 - x1046)
    x1205 = x1050 * x423
    x1206 = 2.0 * x1029 + x1062 * x423
    x1207 = x1206 * x30
    x1208 = x1025 + x1066 * x423
    x1209 = x4 * (x1017 + x1057 * x423 - x7 * (x1021 + x1061 * x423))
    x1210 = x1050 + x1075 * x423
    x1211 = x1210 * x30
    x1212 = x1044 + x1079 * x423
    x1213 = x4 * (x1035 + x1071 * x423 - x7 * (x1039 + x1074 * x423))
    x1214 = x1091 * x423
    x1215 = x4 * x423 * (x1084 - x1094)
    x1216 = 3.0 * x1065 + x1101 * x423
    x1217 = x4 * (3.0 * x1056 + x1098 * x423 - x7 * (3.0 * x1059 + x1099 * x423))
    x1218 = 2.0 * x1078
    x1219 = x1107 * x423 + x1218
    x1220 = 2.0 * x1070
    x1221 = 2.0 * x1073
    x1222 = x4 * (x1104 * x423 + x1220 - x7 * (x1105 * x423 + x1221))
    x1223 = x1091 + x1113 * x423
    x1224 = x4 * (x1083 + x1110 * x423 - x7 * (x1086 + x1111 * x423))
    x1225 = x4 * x423 * (x1116 - x1118)
    x1226 = x608 * x997 + x611
    x1227 = x1226 * x3
    x1228 = x608 * x999 + x614
    x1229 = x1228 * x30
    x1230 = 3.0 * x1229
    x1231 = x608 * x991 + x623 - x7 * (x608 * x994 + x626)
    x1232 = x1231 * x64
    x1233 = x30 * (x1007 * x608 + x621)
    x1234 = x1022 * x608 + x633
    x1235 = x1023 * x608 + x637
    x1236 = x1235 * x30
    x1237 = 2.0 * x1236
    x1238 = x1014 * x608 + x646 - x7 * (x1018 * x608 + x648)
    x1239 = x1238 * x64
    x1240 = x30 * (x1028 * x608 + x644)
    x1241 = x1000 + x1041 * x608 + x671
    x1242 = x1241 * x3
    x1243 = x1008 + x1043 * x608 + x678
    x1244 = x1243 * x30
    x1245 = 2.0 * x1244
    x1246 = x1033 * x608 + x691 - x7 * (x1037 * x608 + x693 + x996) + x993
    x1247 = x1246 * x64
    x1248 = x30 * (x1013 + x1049 * x608 + x688)
    x1249 = x1062 * x608 + x695
    x1250 = x1249 * x30
    x1251 = x1066 * x608 + x703
    x1252 = x1057 * x608 - x7 * (x1061 * x608 + x707) + x705
    x1253 = x1252 * x64
    x1254 = x1029 + x1075 * x608 + x713
    x1255 = x1254 * x30
    x1256 = x1024 + x1079 * x608 + x726
    x1257 = x1016 + x1071 * x608 - x7 * (x1020 + x1074 * x608 + x731) + x729
    x1258 = x1257 * x64
    x1259 = 2.0 * x1050 + x1088 * x608 + x739
    x1260 = x1259 * x30
    x1261 = x1045 + x1092 * x608 + x755
    x1262 = x1261 * x3
    x1263 = x1036 + x1084 * x608 - x7 * (x1040 + x1087 * x608 + x761) + x759
    x1264 = x1263 * x64
    x1265 = x1101 * x608 + x763
    x1266 = x1098 * x608 - x7 * (x1099 * x608 + x767) + x765
    x1267 = x1266 * x64
    x1268 = x1065 + x1107 * x608 + x773
    x1269 = x1056 + x1104 * x608 - x7 * (x1059 + x1105 * x608 + x778) + x776
    x1270 = x1269 * x64
    x1271 = x1113 * x608 + x1218 + x784
    x1272 = x1110 * x608 + x1220 - x7 * (x1111 * x608 + x1221 + x789) + x787
    x1273 = x1272 * x64
    x1274 = 3.0 * x1091 + x1119 * x608 + x797
    x1275 = 3.0 * x1083 + x1116 * x608 - x7 * (3.0 * x1086 + x1117 * x608 + x802) + x800
    x1276 = x1275 * x64
    x1277 = x110 * x817 + x1123 * x423
    x1278 = x30 * (x110 * x822 + x1124 * x423)
    x1279 = x110 * x838 + x1125 + x1129 * x423
    x1280 = x75 * (x110 * x842 + x1127 + x1130 * x423)
    x1281 = x110 * x853 + x1135 * x423
    x1282 = x30 * (x110 * x857 + x1136 * x423)
    x1283 = x30 * (x110 * x865 + 2.0 * x1133 + x1140 * x423)
    x1284 = x110 * x868 + x1131 + x1142 * x423
    x1285 = x30 * (x110 * x876 + x1139 + x1144 * x423)
    x1286 = x110 * x879 + x1137 + x1146 * x423
    x1287 = x30 * (x110 * x889 + x1149 * x423)
    x1288 = x110 * x892 + x1151 * x423
    x1289 = x110 * x897 + 3.0 * x1141 + x1153 * x423
    x1290 = x3 * x448
    x1291 = x110 * x902 + 2.0 * x1145 + x1155 * x423
    x1292 = x3 * x490
    x1293 = x110 * x907 + x1150 + x1157 * x423
    x1294 = x110 * x912 + x1159 * x423
    x1295 = x1161 * x423 + x913
    x1296 = x30 * (x1162 * x423 + x915)
    x1297 = x1163 + x1166 * x423 + x933
    x1298 = x75 * (x1165 + x1167 * x423 + x937)
    x1299 = x1172 * x423 + x939
    x1300 = x30 * (x1173 * x423 + x942)
    x1301 = x30 * (2.0 * x1170 + x1177 * x423 + x950)
    x1302 = x1168 + x1179 * x423 + x953
    x1303 = x30 * (x1176 + x1181 * x423 + x961)
    x1304 = x1174 + x1183 * x423 + x964
    x1305 = x30 * (x1186 * x423 + x968)
    x1306 = x1188 * x423 + x970
    x1307 = 3.0 * x1178 + x1190 * x423 + x975
    x1308 = x3 * x917
    x1309 = 2.0 * x1182 + x1192 * x423 + x983
    x1310 = x3 * x880
    x1311 = x1187 + x1194 * x423 + x988
    x1312 = x1196 * x423 + x989
    x1313 = x1004 + x804 * x997
    x1314 = x30 * (x1012 + x804 * x999)
    x1315 = x1000 * x423 + x1027 + x1199 * x423
    x1316 = x75 * (x1008 * x423 + x1032 + x1200 * x423)
    x1317 = x1041 * x804 + x1048
    x1318 = x30 * (x1043 * x804 + x1054)
    x1319 = x30 * (x1064 + 2.0 * x1203 + x1206 * x423)
    x1320 = x1068 + x1201 + x1208 * x423
    x1321 = x30 * (x1077 + x1205 + x1210 * x423)
    x1322 = x1044 * x423 + x1081 + x1212 * x423
    x1323 = x30 * (x1088 * x804 + x1090)
    x1324 = x1092 * x804 + x1096
    x1325 = x1103 + 3.0 * x1207 + x1216 * x423
    x1326 = x3 * x843
    x1327 = x1109 + 2.0 * x1211 + x1219 * x423
    x1328 = x1171 * x3
    x1329 = x1115 + x1214 + x1223 * x423
    x1330 = x1119 * x804 + x1122
    x1331 = x1229 + x1234 * x423
    x1332 = x75 * (x1233 + x1235 * x423)
    x1333 = x30 * (2.0 * x1240 + x1249 * x423)
    x1334 = x1237 + x1251 * x423
    x1335 = x30 * (x1248 + x1254 * x423)
    x1336 = x1244 + x1256 * x423
    x1337 = x1260 * x423
    x1338 = 3.0 * x1250 + x1265 * x423
    x1339 = 2.0 * x1255
    x1340 = x1268 * x423 + x1339
    x1341 = x1260 + x1271 * x423
    x1342 = x1003 * x110 + x1226 * x608
    x1343 = x30 * (x1011 * x110 + x1228 * x608)
    x1344 = x1026 * x110 + x1234 * x608
    x1345 = x30 * (x1031 * x110 + x1235 * x608)
    x1346 = 2.0 * x1345
    x1347 = x1047 * x110 + x1229 + x1241 * x608
    x1348 = x30 * (x1053 * x110 + x1233 + x1243 * x608)
    x1349 = 2.0 * x1348
    x1350 = x30 * (x1063 * x110 + x1249 * x608)
    x1351 = x1067 * x110 + x1251 * x608
    x1352 = x30 * (x1076 * x110 + x1240 + x1254 * x608)
    x1353 = x1080 * x110 + x1236 + x1256 * x608
    x1354 = x30 * (x1089 * x110 + 2.0 * x1248 + x1259 * x608)
    x1355 = x1095 * x110 + x1245 + x1261 * x608
    x1356 = x110 * x1102 + x1265 * x608
    x1357 = x110 * x1108 + x1250 + x1268 * x608
    x1358 = x110 * x1114 + x1271 * x608 + x1339
    x1359 = x110 * x1121 + 3.0 * x1260 + x1274 * x608
    x1360 = x423 * x448
    x1361 = x423 * x490
    x1362 = 2.0 * x1352

    # 210 item(s)
    result[0, 0] = numpy.sum(
        x131
        * (
            x21
            * (
                x109 * x3
                + x121 * x72
                + x4 * (x107 * x3 + 3.0 * x114 - x116 * x7)
                - x7 * (x125 * x72 + x3 * x99 + x4 * (x116 - x54 * x7))
            )
            + x3
            * (
                x110 * (x109 - x7 * x99)
                + x3
                * (
                    x3 * (x3 * (x60 + x62) + x65 + x71 * x72)
                    - x4 * (x40 * x7 - x54)
                    + x72 * x81
                )
                + x72 * x96
            )
            + x72
            * (
                x110 * (x121 - x125 * x7)
                + x3 * x96
                + x75 * (x3 * x95 + x30 * (x127 + x129) + x4 * (x120 - x124 * x7))
            )
        )
    )
    result[0, 1] = numpy.sum(
        x195
        * (
            x21
            * (
                x185 * x3
                + x191 * x75
                + x4 * (x183 * x3 - x187 * x7 + x189)
                - x7 * (x181 * x3 + x194 * x75 - x4 * (x146 * x7 - x187))
            )
            + x3
            * (
                -x110 * (x181 * x7 - x185)
                + x177 * x75
                + x3
                * (
                    x164 * x75
                    + x3 * (x153 + x159 * x75 + x3 * (x149 + x151))
                    - x4 * (x140 * x7 - x146)
                )
            )
            + x75
            * (
                x110 * (x191 - x194 * x7)
                + x177 * x3
                + x3 * x30 * (x176 + x4 * (2.0 * x10 * x133 * x15 * x33 * x4 * x5 - x173))
            )
        )
    )
    result[0, 2] = numpy.sum(
        x195
        * (
            x21
            * (
                x249 * x3
                + x254 * x75
                + x4 * (x247 * x3 - x251 * x7 + x252)
                - x7 * (x244 * x3 + x256 * x75 - x4 * (x211 * x7 - x251))
            )
            + x3
            * (
                -x110 * (x244 * x7 - x249)
                + x239 * x75
                + x3
                * (
                    x226 * x75
                    + x3 * (x219 + x222 * x75 + x3 * (x215 + x217))
                    - x4 * (x204 * x7 - x211)
                )
            )
            + x75
            * (
                x110 * (x254 - x256 * x7)
                + x239 * x3
                + x3 * x30 * (x238 + x4 * (2.0 * x10 * x15 * x197 * x33 * x4 * x5 - x235))
            )
        )
    )
    result[0, 3] = numpy.sum(
        x195
        * (
            x21
            * (
                x274 * x30
                + x297 * x3
                + x4 * (x298 + x299 * x3 - x301 * x7)
                - x7 * (x271 * x30 + x292 * x3 - x4 * (x281 * x7 - x301))
            )
            + x3
            * (
                -x110 * (x292 * x7 - x297)
                + x268 * x30
                + x3
                * (
                    x267 * x30
                    + x3 * (x288 + x289 + x3 * (x282 + x285))
                    - x4 * (x277 * x7 - x281)
                )
            )
            - x30 * (x110 * (x271 * x7 - x274) - x268 * x3)
        )
    )
    result[0, 4] = numpy.sum(
        x327
        * (
            x21
            * (
                x3 * x323
                + x30 * x311
                + x4 * (x119 * x302 + x197 * x326 * x78 - x325 * x7)
                - x7 * (x3 * x321 + x30 * x309 - x4 * (x314 * x7 - x325))
            )
            + x3
            * (
                -x110 * (x321 * x7 - x323)
                + x3
                * (
                    x3 * (x3 * x316 + x3 * (x302 * x317 * x78 + x316) + x319)
                    + x30 * x306
                    - x4 * (x313 * x7 - x314)
                )
                + x30 * x307
            )
            - x30 * (x110 * (x309 * x7 - x311) - x3 * x307)
        )
    )
    result[0, 5] = numpy.sum(
        x195
        * (
            x21
            * (
                x3 * x366
                + x30 * x345
                + x4 * (x3 * x368 + x367 - x370 * x7)
                - x7 * (x3 * x361 + x30 * x342 - x4 * (x352 * x7 - x370))
            )
            + x3
            * (
                -x110 * (x361 * x7 - x366)
                + x3
                * (
                    x3 * (x3 * (x353 + x355) + x357 + x358)
                    + x30 * x337
                    - x4 * (x348 * x7 - x352)
                )
                + x30 * x338
            )
            - x30 * (x110 * (x342 * x7 - x345) - x3 * x338)
        )
    )
    result[0, 6] = numpy.sum(
        x131
        * x3
        * (
            -x110 * (x381 * x7 - x385)
            - x21
            * (-x385 + x4 * (x379 * x7 - x383) + x7 * (x381 - x4 * (x371 * x7 - x379)))
            + x3 * (x3 * (x376 + x377) + x4 * (x3 * x371 - x374))
        )
    )
    result[0, 7] = numpy.sum(
        x195
        * x3
        * (
            -x110 * (x393 * x7 - x397)
            - x21
            * (-x397 + x4 * (x391 * x7 - x395) + x7 * (x393 - x4 * (x386 * x7 - x391)))
            + x3**2 * (x128 * x389 + x390 + x4 * (x386 - x388))
        )
    )
    result[0, 8] = numpy.sum(
        x195
        * (
            x21
            * (
                x133 * x4 * (x3 * x340 - x364 * x7)
                + x3 * x407
                - x3 * x7 * (x4 * (x133 * x329 - x402 * x7) + x404)
            )
            - x3
            * (
                x110 * (x404 * x7 - x407)
                - x3**2 * (x128 * x401 + x4 * (x133 * x330 - x399) + x400)
            )
        )
    )
    result[0, 9] = numpy.sum(
        x131
        * x3
        * (
            -x110 * (x418 * x7 - x422)
            - x21
            * (x4 * (x416 * x7 - x420) - x422 + x7 * (-x4 * (x408 * x7 - x416) + x418))
            + x3 * (x3 * (x413 + x414) + x4 * (x3 * x408 - x411))
        )
    )
    result[1, 0] = numpy.sum(
        0.5
        * x448
        * (
            x110
            * (
                2.0 * x3 * x429
                + 2.0 * x444 * x72
                + x446
                - x7 * (2.0 * x3 * x426 + 2.0 * x443 * x72 + x447)
            )
            + x3
            * (
                x3 * (2.0 * x3 * (x430 + x431) + x433 + 2.0 * x436 * x72)
                - 2.0 * x4 * (x426 * x7 - x429)
                + 2.0 * x442 * x72
            )
            + x72
            * (
                2.0 * x3 * x442
                - 2.0 * x4 * (x443 * x7 - x444)
                + x75 * (2.0 * x3 * x441 + 2.0 * x423 * x94 + x445)
            )
        )
    )
    result[1, 1] = numpy.sum(
        0.5
        * x490
        * (
            x110
            * (
                2.0 * x3 * x465
                + 2.0 * x485 * x75
                + x488
                - x7 * (2.0 * x3 * x456 + 2.0 * x482 * x75 + x489)
            )
            + x3
            * (
                x3 * (2.0 * x3 * (x3 * x467 + x472) + x473 + 2.0 * x477 * x75)
                - 2.0 * x4 * (x456 * x7 - x465)
                + 2.0 * x479 * x75
            )
            + x75
            * (
                2.0 * x3 * x479
                + x30 * (2.0 * x128 * x475 + x486)
                - 2.0 * x4 * (x482 * x7 - x485)
            )
        )
    )
    result[1, 2] = numpy.sum(
        0.5
        * x490
        * (
            x110
            * (
                2.0 * x3 * x492
                + 2.0 * x501 * x75
                + x503
                - x7 * (2.0 * x3 * x491 + 2.0 * x499 * x75 + x504)
            )
            + x3
            * (
                x3 * (2.0 * x3 * x423 * (x215 + x217) + x493 + 2.0 * x495 * x75)
                - 2.0 * x4 * (x491 * x7 - x492)
                + 2.0 * x497 * x75
            )
            + x75
            * (
                2.0 * x3 * x497
                + x30 * (2.0 * x237 * x423 + x502)
                - 2.0 * x4 * (x499 * x7 - x501)
            )
        )
    )
    result[1, 3] = numpy.sum(
        0.5
        * x490
        * (
            x110
            * (
                2.0 * x3 * x30 * x507
                + 2.0 * x3 * x524
                + x531
                - x7 * (2.0 * x3 * x517 + 2.0 * x3 * x520 + x532)
            )
            + 2.0 * x3 * x30 * (x4 * (x507 - x511) + x516)
            + x3
            * (
                x3 * (2.0 * x3 * x525 + 2.0 * x3 * (x3 * x528 + x525) + x529)
                + 2.0 * x30 * x516
                - 2.0 * x4 * (x520 * x7 - x524)
            )
        )
    )
    result[1, 4] = numpy.sum(
        0.5
        * x558
        * (
            x110
            * (
                2.0 * x3 * x30 * x533
                + 2.0 * x3 * x547
                + x556
                - x7 * (2.0 * x3 * x539 + 2.0 * x3 * x543 + x557)
            )
            + 2.0 * x3 * x30 * (x4 * (x533 - x535) + x538)
            + x3
            * (
                x3 * (2.0 * x3 * x548 + 2.0 * x3 * (x3 * x551 + x548) + x552)
                + 2.0 * x30 * x538
                - 2.0 * x4 * (x543 * x7 - x547)
            )
        )
    )
    result[1, 5] = numpy.sum(
        0.5
        * x490
        * (
            x110
            * (
                2.0 * x3 * x564
                + 2.0 * x365 * x423
                + x567
                - x7 * (2.0 * x3 * x562 + 2.0 * x360 * x423 + x568)
            )
            + x3
            * (
                x3 * (2.0 * x3 * (x355 * x423 + x565) + 2.0 * x358 * x423 + x566)
                + 2.0 * x30 * x560
                - 2.0 * x4 * (x562 * x7 - x564)
            )
            + 2.0 * x30 * (x3 * x560 + x4 * x423 * (x3 * x329 - x332))
        )
    )
    result[1, 6] = numpy.sum(
        0.5
        * x448
        * (
            x110 * (2.0 * x128 * x571 + x581 - x7 * (2.0 * x128 * x574 + x582))
            + x3**2 * (2.0 * x128 * x578 + 2.0 * x4 * (x571 - x575) + x579)
        )
    )
    result[1, 7] = numpy.sum(
        0.5
        * x490
        * (
            x110 * (2.0 * x128 * x584 + x594 - x7 * (2.0 * x128 * x586 + x595))
            + x3**2 * (2.0 * x128 * x589 + 2.0 * x4 * (x584 - x587) + x590)
        )
    )
    result[1, 8] = numpy.sum(
        0.5
        * x490
        * (
            x110 * (2.0 * x128 * x596 + x603 - x7 * (2.0 * x128 * x597 + x604))
            + x3**2 * (2.0 * x128 * x599 + 2.0 * x4 * (x596 - x598) + x600)
        )
    )
    result[1, 9] = numpy.sum(
        0.5
        * x448
        * (
            x110 * (2.0 * x419 * x423 + x606 - x7 * (2.0 * x415 * x423 + x607))
            + x3
            * (x3 * (2.0 * x413 * x423 + x605) + 2.0 * x4 * x423 * (x3 * x408 - x411))
        )
    )
    result[2, 0] = numpy.sum(
        x448
        * (
            x110 * (x3 * x610 + x620 * x72 + x624 - x7 * (x3 * x609 + x619 * x72 + x627))
            + x3
            * (
                x3 * (x3 * x608 * (x60 + x62) + x612 + x613 * x72)
                - x4 * (x609 * x7 - x610)
                + x618 * x72
            )
            + x72
            * (
                x3 * x618
                - x4 * (x619 * x7 - x620)
                + x75 * (x3 * x617 + x608 * x94 + x622)
            )
        )
    )
    result[2, 1] = numpy.sum(
        x490
        * (
            x110 * (x3 * x631 + x643 * x75 + x647 - x7 * (x3 * x629 + x641 * x75 + x649))
            + x3
            * (
                x3 * (x3 * (x149 * x608 + x632) + x634 + x636 * x75)
                - x4 * (x629 * x7 - x631)
                + x639 * x75
            )
            + x75 * (x3 * x639 + x30 * (x175 * x608 + x645) - x4 * (x641 * x7 - x643))
        )
    )
    result[2, 2] = numpy.sum(
        x490
        * (
            x110 * (x3 * x663 + x687 * x75 + x692 - x7 * (x3 * x656 + x684 * x75 + x694))
            + x3
            * (
                x3 * (x3 * (x665 + x669) + x672 + x676 * x75)
                - x4 * (x656 * x7 - x663)
                + x681 * x75
            )
            + x75 * (x3 * x681 + x30 * (x128 * x673 + x689) - x4 * (x684 * x7 - x687))
        )
    )
    result[2, 3] = numpy.sum(
        x490
        * (
            x110
            * (x296 * x608 + x3 * x701 - x7 * (x291 * x608 + x3 * x699 + x708) + x706)
            + x3
            * (
                x3 * (x289 * x608 + x3 * (x285 * x608 + x702) + x704)
                + x30 * x697
                - x4 * (x699 * x7 - x701)
            )
            + x30 * (x3 * x697 + x4 * x608 * (x258 * x3 - x261))
        )
    )
    result[2, 4] = numpy.sum(
        x558
        * (
            x110
            * (x3 * x30 * x709 + x3 * x722 - x7 * (x3 * x716 + x3 * x719 + x732) + x730)
            + x3 * x30 * (x4 * (x709 - x711) + x715)
            + x3
            * (
                x3 * (x3 * x723 + x3 * (x3 * x725 + x723) + x727)
                + x30 * x715
                - x4 * (x7 * x719 - x722)
            )
        )
    )
    result[2, 5] = numpy.sum(
        x490
        * (
            x110
            * (x3 * x30 * x733 + x3 * x749 - x7 * (x3 * x742 + x3 * x745 + x762) + x760)
            + x3 * x30 * (x4 * (x733 - x735) + x741)
            + x3
            * (
                x3 * (x3 * (x750 + x753) + x756 + x757)
                + x30 * x741
                - x4 * (x7 * x745 - x749)
            )
        )
    )
    result[2, 6] = numpy.sum(
        x448
        * (
            x110 * (x382 * x608 - x7 * (x378 * x608 + x768) + x766)
            + x3 * (x3 * (x376 * x608 + x764) + x4 * x608 * (x3 * x371 - x374))
        )
    )
    result[2, 7] = numpy.sum(
        x490
        * (
            x110 * (x128 * x769 - x7 * (x128 * x770 + x779) + x777)
            + x3**2 * (x128 * x772 + x4 * (x769 - x771) + x774)
        )
    )
    result[2, 8] = numpy.sum(
        x490
        * (
            x110 * (x128 * x780 - x7 * (x128 * x781 + x790) + x788)
            + x3**2 * (x128 * x783 + x4 * (x780 - x782) + x785)
        )
    )
    result[2, 9] = numpy.sum(
        x448
        * (
            x110 * (x128 * x791 - x7 * (x128 * x792 + x803) + x801)
            + x3 * (x3 * (x796 + x798) + x4 * (x3 * x791 - x794))
        )
    )
    result[3, 0] = numpy.sum(
        x824
        * (
            x3 * (x3 * (x3 * x813 + 3.0 * x816) + x64 * x817 + x72 * x821)
            + x4 * (x3 * x805 - x7 * (x3 * x809 + 3.0 * x811) + 3.0 * x807)
            + x72 * (x3 * x821 + x64 * x822 + x75 * (x3 * x819 + x823))
        )
    )
    result[3, 1] = numpy.sum(
        x843
        * (
            x3 * (x3 * (x3 * x834 + x837) + x64 * x838 + x75 * x841)
            + x4 * (x3 * x826 - x7 * (x3 * x830 + x832) + x828)
            + x75 * (x3 * x840 + x3 * x841 + x64 * x842)
        )
    )
    result[3, 2] = numpy.sum(
        x843
        * (
            x3 * (x3 * (x3 * x850 + 2.0 * x852) + x64 * x853 + x75 * x856)
            + x4 * (x3 * x844 - x7 * (x3 * x847 + 2.0 * x849) + 2.0 * x846)
            + x75 * (x3 * x855 + x3 * x856 + x64 * x857)
        )
    )
    result[3, 3] = numpy.sum(
        x843
        * (
            x3 * (x3 * x866 + x3 * (x3 * x867 + x866) + x64 * x868)
            + x30 * (x128 * x864 + x64 * x865)
            + x4 * (x3 * x860 - x7 * (x3 * x863 + x862) + x859)
        )
    )
    result[3, 4] = numpy.sum(
        x880
        * (
            x3 * (x3 * x877 + x3 * (x3 * x878 + x877) + x64 * x879)
            + x30 * (x128 * x875 + x64 * x876)
            + x4 * (x3 * x871 - x7 * (x3 * x874 + x873) + x870)
        )
    )
    result[3, 5] = numpy.sum(
        x843
        * (
            x3 * (x3 * x890 + x3 * (x3 * x891 + x890) + x64 * x892)
            + x30 * (x128 * x888 + x64 * x889)
            + x4 * (x3 * x883 - x7 * (x3 * x887 + x885) + x882)
        )
    )
    result[3, 6] = numpy.sum(x3 * x824 * (x128 * x896 + x4 * (x893 - x895) + x64 * x897))
    result[3, 7] = numpy.sum(x3 * x843 * (x128 * x901 + x4 * (x898 - x900) + x64 * x902))
    result[3, 8] = numpy.sum(x3 * x843 * (x128 * x906 + x4 * (x903 - x905) + x64 * x907))
    result[3, 9] = numpy.sum(x3 * x824 * (x128 * x911 + x4 * (x908 - x910) + x64 * x912))
    result[4, 0] = numpy.sum(
        0.5
        * x917
        * (
            x3 * (2.0 * x3 * x608 * (x430 + x431) + 2.0 * x72 * x914 + x913)
            + 2.0 * x4 * x608 * (x427 + x428 * x72 - x7 * (x424 + x425))
            + x72 * (2.0 * x3 * x914 + 2.0 * x75 * (x440 * x608 + x916) + x915)
        )
    )
    result[4, 1] = numpy.sum(
        0.5
        * x880
        * (
            x3 * (2.0 * x3 * (x3 * x929 + x932) + 2.0 * x75 * x936 + x933)
            + 2.0 * x4 * (x3 * x919 - x7 * (x3 * x924 + x927) + x922)
            + x75 * (2.0 * x3 * x935 + 2.0 * x3 * x936 + x937)
        )
    )
    result[4, 2] = numpy.sum(
        0.5
        * x880
        * (
            x3 * (2.0 * x3 * x423 * (x665 + x669) + 2.0 * x75 * x941 + x939)
            + 2.0 * x4 * (x423 * x658 - x423 * x7 * (x651 + x655) + x75 * x938)
            + x75 * (2.0 * x3 * x941 + 2.0 * x423 * x680 + x942)
        )
    )
    result[4, 3] = numpy.sum(
        0.5
        * x880
        * (
            x3 * (2.0 * x3 * x951 + 2.0 * x3 * (x3 * x952 + x951) + x953)
            + x30 * (2.0 * x128 * x949 + x950)
            + 2.0 * x4 * (x3 * x945 - x7 * (x3 * x948 + x947) + x944)
        )
    )
    result[4, 4] = numpy.sum(
        0.5
        * x965
        * (
            x3 * (2.0 * x3 * x962 + 2.0 * x3 * (x3 * x963 + x962) + x964)
            + x30 * (2.0 * x128 * x960 + x961)
            + 2.0 * x4 * (x3 * x956 - x7 * (x3 * x959 + x958) + x955)
        )
    )
    result[4, 5] = numpy.sum(
        0.5
        * x880
        * (
            x3 * (2.0 * x3 * (x423 * x753 + x969) + 2.0 * x423 * x757 + x970)
            + x30 * (2.0 * x423 * x738 + x968)
            + 2.0 * x4 * (x423 * x748 - x7 * (x423 * x744 + x967) + x966)
        )
    )
    result[4, 6] = numpy.sum(
        0.5 * x3 * x917 * (2.0 * x128 * x974 + 2.0 * x4 * (x971 - x973) + x975)
    )
    result[4, 7] = numpy.sum(
        0.5 * x3 * x880 * (2.0 * x128 * x982 + 2.0 * x4 * (x977 - x980) + x983)
    )
    result[4, 8] = numpy.sum(
        0.5 * x3 * x880 * (2.0 * x128 * x987 + 2.0 * x4 * (x984 - x986) + x988)
    )
    result[4, 9] = numpy.sum(
        0.5
        * x917
        * (x3 * (2.0 * x423 * x796 + x989) + 2.0 * x4 * x423 * (x3 * x791 - x794))
    )
    result[5, 0] = numpy.sum(
        x824
        * (
            x3 * (x1004 + x1010 * x72 + x3 * (x1001 + x998))
            + x4 * (x3 * x991 - x7 * (x3 * x994 + 3.0 * x996) + 3.0 * x993)
            + x72 * (x1010 * x3 + x1012 + x75 * (x1007 * x3 + x1013))
        )
    )
    result[5, 1] = numpy.sum(
        x843
        * (
            x3 * (x1027 + x1030 * x75 + x3 * (x1022 * x3 + x1025))
            + x4 * (x1014 * x3 + x1017 - x7 * (x1018 * x3 + x1021))
            + x75 * (x1029 * x3 + x1030 * x3 + x1032)
        )
    )
    result[5, 2] = numpy.sum(
        x843
        * (
            x3 * (x1048 + x1052 * x75 + x3 * (x1042 + x1045))
            + x4 * (x1033 * x3 + x1036 - x7 * (x1037 * x3 + x1040))
            + x75 * (x1050 * x3 + x1052 * x3 + x1054)
        )
    )
    result[5, 3] = numpy.sum(
        x843
        * (
            x3 * (x1065 * x3 + x1068 + x3 * (x1065 + x1066 * x3))
            + x30 * (x1062 * x128 + x1064)
            + x4 * (x1056 + x1057 * x3 - x7 * (x1059 + x1061 * x3))
        )
    )
    result[5, 4] = numpy.sum(
        x880
        * (
            x3 * (x1078 * x3 + x1081 + x3 * (x1078 + x1079 * x3))
            + x30 * (x1075 * x128 + x1077)
            + x4 * (x1070 + x1071 * x3 - x7 * (x1073 + x1074 * x3))
        )
    )
    result[5, 5] = numpy.sum(
        x843
        * (
            x3 * (x1096 + x1097 + x3 * (x1091 + x1093))
            + x30 * (x1088 * x128 + x1090)
            + x4 * (x1083 + x1084 * x3 - x7 * (x1086 + x1087 * x3))
        )
    )
    result[5, 6] = numpy.sum(x3 * x824 * (x1101 * x128 + x1103 + x4 * (x1098 - x1100)))
    result[5, 7] = numpy.sum(x3 * x843 * (x1107 * x128 + x1109 + x4 * (x1104 - x1106)))
    result[5, 8] = numpy.sum(x3 * x843 * (x1113 * x128 + x1115 + x4 * (x1110 - x1112)))
    result[5, 9] = numpy.sum(x3 * x824 * (x1120 + x1122 + x4 * (x1116 - x1118)))
    result[6, 0] = numpy.sum(
        x1128
        * (
            x1126 * x64
            + x3 * (x1123 * x3 + 3.0 * x1125)
            + x72 * (x1124 * x3 + 2.0 * x1127)
        )
    )
    result[6, 1] = numpy.sum(
        x1134 * (x1132 * x64 + x3 * (x1129 * x3 + x1131) + x75 * (x1130 * x3 + x1133))
    )
    result[6, 2] = numpy.sum(
        x1134
        * (x1138 * x64 + x3 * (x1135 * x3 + 2.0 * x1137) + x75 * (x1136 * x3 + x1139))
    )
    result[6, 3] = numpy.sum(
        x1134 * (x1141 * x3 + x1143 * x64 + x3 * (x1141 + x1142 * x3))
    )
    result[6, 4] = numpy.sum(
        x1148 * (x1145 * x3 + x1147 * x64 + x3 * (x1145 + x1146 * x3))
    )
    result[6, 5] = numpy.sum(
        x1134 * (x1150 * x3 + x1152 * x64 + x3 * (x1150 + x1151 * x3))
    )
    result[6, 6] = numpy.sum(x1128 * (x1153 * x128 + x1154 * x64))
    result[6, 7] = numpy.sum(x1134 * (x1155 * x128 + x1156 * x64))
    result[6, 8] = numpy.sum(x1134 * (x1157 * x128 + x1158 * x64))
    result[6, 9] = numpy.sum(x1128 * (x1159 * x128 + x1160 * x64))
    result[7, 0] = numpy.sum(
        x843
        * (
            x1164 * x64
            + x3 * (x1161 * x3 + 3.0 * x1163)
            + x72 * (x1162 * x3 + 2.0 * x1165)
        )
    )
    result[7, 1] = numpy.sum(
        x1171 * (x1169 * x64 + x3 * (x1166 * x3 + x1168) + x75 * (x1167 * x3 + x1170))
    )
    result[7, 2] = numpy.sum(
        x1171
        * (x1175 * x64 + x3 * (x1172 * x3 + 2.0 * x1174) + x75 * (x1173 * x3 + x1176))
    )
    result[7, 3] = numpy.sum(
        x1171 * (x1178 * x3 + x1180 * x64 + x3 * (x1178 + x1179 * x3))
    )
    result[7, 4] = numpy.sum(
        x1185 * (x1182 * x3 + x1184 * x64 + x3 * (x1182 + x1183 * x3))
    )
    result[7, 5] = numpy.sum(
        x1171 * (x1187 * x3 + x1189 * x64 + x3 * (x1187 + x1188 * x3))
    )
    result[7, 6] = numpy.sum(x843 * (x1190 * x128 + x1191 * x64))
    result[7, 7] = numpy.sum(x1171 * (x1192 * x128 + x1193 * x64))
    result[7, 8] = numpy.sum(x1171 * (x1194 * x128 + x1195 * x64))
    result[7, 9] = numpy.sum(x843 * (x1196 * x128 + x1197 * x64))
    result[8, 0] = numpy.sum(
        0.5
        * x843
        * (x1198 + 2.0 * x3 * x423 * (x1001 + x998) + 2.0 * x423 * x72 * (x1005 + x1009))
    )
    result[8, 1] = numpy.sum(
        0.5
        * x1171
        * (x1202 + 2.0 * x3 * (x1199 * x3 + x1201) + 2.0 * x75 * (x1200 * x3 + x1203))
    )
    result[8, 2] = numpy.sum(
        0.5
        * x1171
        * (x1204 + 2.0 * x3 * x423 * (x1042 + x1045) + 2.0 * x75 * (x1051 * x423 + x1205))
    )
    result[8, 3] = numpy.sum(
        0.5 * x1171 * (4.0 * x1207 * x3 + 2.0 * x1208 * x3**2 + x1209)
    )
    result[8, 4] = numpy.sum(
        0.5 * x1185 * (4.0 * x1211 * x3 + 2.0 * x1212 * x3**2 + x1213)
    )
    result[8, 5] = numpy.sum(
        0.5 * x1171 * (2.0 * x1097 * x423 + x1215 + 2.0 * x3 * (x1093 * x423 + x1214))
    )
    result[8, 6] = numpy.sum(0.5 * x843 * (2.0 * x1216 * x128 + x1217))
    result[8, 7] = numpy.sum(0.5 * x1171 * (2.0 * x1219 * x128 + x1222))
    result[8, 8] = numpy.sum(0.5 * x1171 * (2.0 * x1223 * x128 + x1224))
    result[8, 9] = numpy.sum(0.5 * x843 * (2.0 * x1120 * x423 + x1225))
    result[9, 0] = numpy.sum(
        x1128 * (x1232 + x3 * (x1227 + x1230) + x72 * (x1228 * x3 + 2.0 * x1233))
    )
    result[9, 1] = numpy.sum(
        x1134 * (x1239 + x3 * (x1234 * x3 + x1237) + x75 * (x1235 * x3 + x1240))
    )
    result[9, 2] = numpy.sum(
        x1134 * (x1247 + x3 * (x1242 + x1245) + x75 * (x1243 * x3 + x1248))
    )
    result[9, 3] = numpy.sum(x1134 * (x1250 * x3 + x1253 + x3 * (x1250 + x1251 * x3)))
    result[9, 4] = numpy.sum(x1148 * (x1255 * x3 + x1258 + x3 * (x1255 + x1256 * x3)))
    result[9, 5] = numpy.sum(x1134 * (x1260 * x3 + x1264 + x3 * (x1260 + x1262)))
    result[9, 6] = numpy.sum(x1128 * (x1265 * x128 + x1267))
    result[9, 7] = numpy.sum(x1134 * (x1268 * x128 + x1270))
    result[9, 8] = numpy.sum(x1134 * (x1271 * x128 + x1273))
    result[9, 9] = numpy.sum(x1128 * (x1274 * x128 + x1276))
    result[10, 0] = numpy.sum(x448 * (x1277 * x3 + 3.0 * x1278))
    result[10, 1] = numpy.sum(x490 * (x1279 * x3 + x1280))
    result[10, 2] = numpy.sum(x490 * (x1281 * x3 + 2.0 * x1282))
    result[10, 3] = numpy.sum(x490 * (x1283 + x1284 * x3))
    result[10, 4] = numpy.sum(x558 * (x1285 + x1286 * x3))
    result[10, 5] = numpy.sum(x490 * (x1287 + x1288 * x3))
    result[10, 6] = numpy.sum(x1289 * x1290)
    result[10, 7] = numpy.sum(x1291 * x1292)
    result[10, 8] = numpy.sum(x1292 * x1293)
    result[10, 9] = numpy.sum(x1290 * x1294)
    result[11, 0] = numpy.sum(x917 * (x1295 * x3 + 3.0 * x1296))
    result[11, 1] = numpy.sum(x880 * (x1297 * x3 + x1298))
    result[11, 2] = numpy.sum(x880 * (x1299 * x3 + 2.0 * x1300))
    result[11, 3] = numpy.sum(x880 * (x1301 + x1302 * x3))
    result[11, 4] = numpy.sum(x965 * (x1303 + x1304 * x3))
    result[11, 5] = numpy.sum(x880 * (x1305 + x1306 * x3))
    result[11, 6] = numpy.sum(x1307 * x1308)
    result[11, 7] = numpy.sum(x1309 * x1310)
    result[11, 8] = numpy.sum(x1310 * x1311)
    result[11, 9] = numpy.sum(x1308 * x1312)
    result[12, 0] = numpy.sum(x843 * (x1313 * x3 + 3.0 * x1314))
    result[12, 1] = numpy.sum(x1171 * (x1315 * x3 + x1316))
    result[12, 2] = numpy.sum(x1171 * (x1317 * x3 + 2.0 * x1318))
    result[12, 3] = numpy.sum(x1171 * (x1319 + x1320 * x3))
    result[12, 4] = numpy.sum(x1185 * (x1321 + x1322 * x3))
    result[12, 5] = numpy.sum(x1171 * (x1323 + x1324 * x3))
    result[12, 6] = numpy.sum(x1325 * x1326)
    result[12, 7] = numpy.sum(x1327 * x1328)
    result[12, 8] = numpy.sum(x1328 * x1329)
    result[12, 9] = numpy.sum(x1326 * x1330)
    result[13, 0] = numpy.sum(x423 * x917 * (x1227 + x1230))
    result[13, 1] = numpy.sum(x880 * (x1331 * x3 + x1332))
    result[13, 2] = numpy.sum(x423 * x880 * (x1242 + x1245))
    result[13, 3] = numpy.sum(x880 * (x1333 + x1334 * x3))
    result[13, 4] = numpy.sum(x965 * (x1335 + x1336 * x3))
    result[13, 5] = numpy.sum(x880 * (x1262 * x423 + x1337))
    result[13, 6] = numpy.sum(x1308 * x1338)
    result[13, 7] = numpy.sum(x1310 * x1340)
    result[13, 8] = numpy.sum(x1310 * x1341)
    result[13, 9] = numpy.sum(x1274 * x1308 * x423)
    result[14, 0] = numpy.sum(x448 * (x1342 * x3 + 3.0 * x1343))
    result[14, 1] = numpy.sum(x490 * (x1344 * x3 + x1346))
    result[14, 2] = numpy.sum(x490 * (x1347 * x3 + x1349))
    result[14, 3] = numpy.sum(x490 * (x1350 + x1351 * x3))
    result[14, 4] = numpy.sum(x558 * (x1352 + x1353 * x3))
    result[14, 5] = numpy.sum(x490 * (x1354 + x1355 * x3))
    result[14, 6] = numpy.sum(x1290 * x1356)
    result[14, 7] = numpy.sum(x1292 * x1357)
    result[14, 8] = numpy.sum(x1292 * x1358)
    result[14, 9] = numpy.sum(x1290 * x1359)
    result[15, 0] = numpy.sum(x131 * (x1126 * x21 + x1277 * x423))
    result[15, 1] = numpy.sum(x195 * (x1132 * x21 + x1278 + x1279 * x423))
    result[15, 2] = numpy.sum(x195 * (x1138 * x21 + x1281 * x423))
    result[15, 3] = numpy.sum(x195 * (x1143 * x21 + x1280 + x1284 * x423))
    result[15, 4] = numpy.sum(x327 * (x1147 * x21 + x1282 + x1286 * x423))
    result[15, 5] = numpy.sum(x195 * (x1152 * x21 + x1288 * x423))
    result[15, 6] = numpy.sum(x131 * (x1154 * x21 + 3.0 * x1283 + x1289 * x423))
    result[15, 7] = numpy.sum(x195 * (x1156 * x21 + 2.0 * x1285 + x1291 * x423))
    result[15, 8] = numpy.sum(x195 * (x1158 * x21 + x1287 + x1293 * x423))
    result[15, 9] = numpy.sum(x131 * (x1160 * x21 + x1294 * x423))
    result[16, 0] = numpy.sum(x448 * (x110 * x1164 + x1295 * x423))
    result[16, 1] = numpy.sum(x490 * (x110 * x1169 + x1296 + x1297 * x423))
    result[16, 2] = numpy.sum(x490 * (x110 * x1175 + x1299 * x423))
    result[16, 3] = numpy.sum(x490 * (x110 * x1180 + x1298 + x1302 * x423))
    result[16, 4] = numpy.sum(x558 * (x110 * x1184 + x1300 + x1304 * x423))
    result[16, 5] = numpy.sum(x490 * (x110 * x1189 + x1306 * x423))
    result[16, 6] = numpy.sum(x448 * (x110 * x1191 + 3.0 * x1301 + x1307 * x423))
    result[16, 7] = numpy.sum(x490 * (x110 * x1193 + 2.0 * x1303 + x1309 * x423))
    result[16, 8] = numpy.sum(x490 * (x110 * x1195 + x1305 + x1311 * x423))
    result[16, 9] = numpy.sum(x448 * (x110 * x1197 + x1312 * x423))
    result[17, 0] = numpy.sum(x824 * (x1198 + x1313 * x423))
    result[17, 1] = numpy.sum(x843 * (x1202 + x1314 + x1315 * x423))
    result[17, 2] = numpy.sum(x843 * (x1204 + x1317 * x423))
    result[17, 3] = numpy.sum(x843 * (x1209 + x1316 + x1320 * x423))
    result[17, 4] = numpy.sum(x880 * (x1213 + x1318 + x1322 * x423))
    result[17, 5] = numpy.sum(x843 * (x1215 + x1324 * x423))
    result[17, 6] = numpy.sum(x824 * (x1217 + 3.0 * x1319 + x1325 * x423))
    result[17, 7] = numpy.sum(x843 * (x1222 + 2.0 * x1321 + x1327 * x423))
    result[17, 8] = numpy.sum(x843 * (x1224 + x1323 + x1329 * x423))
    result[17, 9] = numpy.sum(x824 * (x1225 + x1330 * x423))
    result[18, 0] = numpy.sum(x1128 * (x1226 * x804 + x1232))
    result[18, 1] = numpy.sum(x1134 * (x1229 * x423 + x1239 + x1331 * x423))
    result[18, 2] = numpy.sum(x1134 * (x1241 * x804 + x1247))
    result[18, 3] = numpy.sum(x1134 * (x1253 + x1332 + x1334 * x423))
    result[18, 4] = numpy.sum(x1148 * (x1244 * x423 + x1258 + x1336 * x423))
    result[18, 5] = numpy.sum(x1134 * (x1261 * x804 + x1264))
    result[18, 6] = numpy.sum(x1128 * (x1267 + 3.0 * x1333 + x1338 * x423))
    result[18, 7] = numpy.sum(x1134 * (x1270 + 2.0 * x1335 + x1340 * x423))
    result[18, 8] = numpy.sum(x1134 * (x1273 + x1337 + x1341 * x423))
    result[18, 9] = numpy.sum(x1128 * (x1274 * x804 + x1276))
    result[19, 0] = numpy.sum(x1342 * x1360)
    result[19, 1] = numpy.sum(x490 * (x1343 + x1344 * x423))
    result[19, 2] = numpy.sum(x1347 * x1361)
    result[19, 3] = numpy.sum(x490 * (x1346 + x1351 * x423))
    result[19, 4] = numpy.sum(x558 * (x1348 + x1353 * x423))
    result[19, 5] = numpy.sum(x1355 * x1361)
    result[19, 6] = numpy.sum(x448 * (3.0 * x1350 + x1356 * x423))
    result[19, 7] = numpy.sum(x490 * (x1357 * x423 + x1362))
    result[19, 8] = numpy.sum(x490 * (x1354 + x1358 * x423))
    result[19, 9] = numpy.sum(x1359 * x1360)
    result[20, 0] = numpy.sum(x131 * (x1231 * x21 + x1342 * x608))
    result[20, 1] = numpy.sum(x195 * (x1238 * x21 + x1344 * x608))
    result[20, 2] = numpy.sum(x195 * (x1246 * x21 + x1343 + x1347 * x608))
    result[20, 3] = numpy.sum(x195 * (x1252 * x21 + x1351 * x608))
    result[20, 4] = numpy.sum(x327 * (x1257 * x21 + x1345 + x1353 * x608))
    result[20, 5] = numpy.sum(x195 * (x1263 * x21 + x1349 + x1355 * x608))
    result[20, 6] = numpy.sum(x131 * (x1266 * x21 + x1356 * x608))
    result[20, 7] = numpy.sum(x195 * (x1269 * x21 + x1350 + x1357 * x608))
    result[20, 8] = numpy.sum(x195 * (x1272 * x21 + x1358 * x608 + x1362))
    result[20, 9] = numpy.sum(x131 * (x1275 * x21 + 3.0 * x1354 + x1359 * x608))
    return result


def _2center2el3d_54(ax, da, A, bx, db, B):
    """Cartesian (h|g) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((21, 15), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = ax ** (-1.0)
    x5 = -x2 - B[0]
    x6 = bx ** (-1.0)
    x7 = ax * x1
    x8 = bx * x7 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x9 = boys(7, x8)
    x10 = 17.4934183276249
    x11 = x10 * x6
    x12 = x0 ** (-1.5)
    x13 = 2.0 * x12
    x14 = x11 * x13
    x15 = x14 * x9
    x16 = x15 * x5
    x17 = x0 ** (-0.5)
    x18 = boys(6, x8)
    x19 = 0.5 * x6
    x20 = x19 * (2.0 * x10 * x17 * x18 * x4 * x6 - x15)
    x21 = boys(8, x8)
    x22 = x5**2
    x23 = 2.0 * x4
    x24 = x17 * x23
    x25 = x11 * x24
    x26 = x22 * x25
    x27 = x21 * x26
    x28 = x20 + x27
    x29 = x28 * x5 + x6 * (2.0 * x10 * x17 * x18 * x4 * x5 * x6 - x16)
    x30 = x11 * x18
    x31 = x13 * x30
    x32 = boys(5, x8)
    x33 = x19 * (2.0 * x10 * x17 * x32 * x4 * x6 - x31)
    x34 = x25 * x9
    x35 = x22 * x34
    x36 = x33 + x35
    x37 = x14 * x32
    x38 = boys(4, x8)
    x39 = x19 * (2.0 * x10 * x17 * x38 * x4 * x6 - x37)
    x40 = x24 * x30
    x41 = x22 * x40
    x42 = x39 + x41
    x43 = 1.5 * x6
    x44 = x29 * x5 - x43 * (x36 * x7 - x42)
    x45 = x3 * x44
    x46 = 0.5 / (ax + bx)
    x47 = x31 * x5
    x48 = x36 * x5 + x6 * (2.0 * x10 * x17 * x32 * x4 * x5 * x6 - x47)
    x49 = x46 * x48
    x50 = 4.0 * x49
    x51 = x45 + x50
    x52 = bx * x1
    x53 = x14 * x38
    x54 = boys(3, x8)
    x55 = x19 * (2.0 * x10 * x17 * x4 * x54 * x6 - x53)
    x56 = x25 * x32
    x57 = x22 * x56
    x58 = x55 + x57
    x59 = -x43 * (x42 * x7 - x58) + x48 * x5
    x60 = x3 * x59
    x61 = x37 * x5
    x62 = -2.0 * x10 * x17 * x38 * x4 * x5 * x6
    x63 = x42 * x5 - x6 * (x61 + x62)
    x64 = x46 * x63
    x65 = 4.0 * x64
    x66 = x60 + x65
    x67 = x14 * x21
    x68 = x5 * x67
    x69 = x19 * (2.0 * x10 * x17 * x4 * x6 * x9 - x67)
    x70 = boys(9, x8)
    x71 = x26 * x70
    x72 = -x43 * (x28 * x7 - x36) + x5 * (
        x5 * (x69 + x71) + x6 * (2.0 * x10 * x17 * x4 * x5 * x6 * x9 - x68)
    )
    x73 = x3 * x72
    x74 = x29 * x46
    x75 = 4.0 * x74
    x76 = x44 * x52
    x77 = 0.5 * x4
    x78 = x77 * (x59 - x76)
    x79 = x29 * x3
    x80 = x36 * x46
    x81 = 3.0 * x80
    x82 = x79 + x81
    x83 = 4.0 * x46
    x84 = x48 * x52
    x85 = x77 * (x63 - x84)
    x86 = x3 * x36
    x87 = x5 * x83
    x88 = x17 * x4
    x89 = x30 * x88
    x90 = x87 * x89
    x91 = x86 + x90
    x92 = 3.0 * x46
    x93 = x3 * x82 + x85 + x91 * x92
    x94 = x3 * x48
    x95 = x42 * x46
    x96 = 3.0 * x95
    x97 = x94 + x96
    x98 = x3 * x63
    x99 = x46 * x58
    x100 = 3.0 * x99
    x101 = x100 + x98
    x102 = x42 * x52
    x103 = x77 * (-x102 + x58)
    x104 = 2.0 * x46
    x105 = x11 * x88
    x106 = x104 * x105
    x107 = x106 * x32
    x108 = x3 * x5
    x109 = x108 * x40
    x110 = x107 + x109
    x111 = x103 + x104 * x110 + x3 * x91
    x112 = x111 * x92 + x3 * x93 + x4 * (x101 - x52 * x97)
    x113 = x14 * x54
    x114 = boys(2, x8)
    x115 = x19 * (2.0 * x10 * x114 * x17 * x4 * x6 - x113)
    x116 = x25 * x38
    x117 = x116 * x22
    x118 = x115 + x117
    x119 = x43 * (x118 - x58 * x7) + x5 * x63
    x120 = x77 * (x119 - x52 * x59)
    x121 = x120 + x3 * x51 + x83 * x97
    x122 = x5 * x53
    x123 = x5 * x58 + x6 * (2.0 * x10 * x17 * x4 * x5 * x54 * x6 - x122)
    x124 = x19 * (2.0 * x10 * x17 * x4 * x6 * boys(1, x8) - x114 * x14)
    x125 = x25 * x54
    x126 = x123 * x5 + x43 * (-x118 * x7 + x124 + x125 * x22)
    x127 = x77 * (-x119 * x52 + x126)
    x128 = x101 * x83 + x127 + x3 * x66
    x129 = 1.5 * x4
    x130 = x5 * (x118 + x6 * (2.0 * x10 * x114 * x17 * x4 * x6 - x113))
    x131 = x130 * x46
    x132 = x123 * x46
    x133 = x119 * x3 + 4.0 * x132
    x134 = x77 * (-x123 * x52 + x130)
    x135 = x105 * x38
    x136 = x135 * x87 + x3 * x58
    x137 = x101 * x3 + x134 + x136 * x92
    x138 = x77 * (x123 - x52 * x63)
    x139 = x32 * x5
    x140 = x105 * x83
    x141 = x139 * x140
    x142 = x141 + x3 * x42
    x143 = x138 + x142 * x92 + x3 * x97
    x144 = x10 * x12 * x23
    x145 = x139 * x144
    x146 = -x77 * (x145 + x62)
    x147 = x107 * x3
    x148 = 0.179587122125167 * da * db * numpy.sqrt(ax**6.5) * numpy.sqrt(bx**5.5)
    x149 = 4.59731646942873 * x148
    x150 = -x1 * (ax * A[1] + bx * B[1])
    x151 = -x150 - B[1]
    x152 = x15 * x151
    x153 = x6 * (2.0 * x10 * x151 * x17 * x18 * x4 * x6 - x152)
    x154 = x151 * x27
    x155 = 0.5 * x153 + x154
    x156 = x151 * x6 * (2.0 * x10 * x17 * x18 * x4 * x5 * x6 - x16) + x155 * x5
    x157 = x156 * x3
    x158 = x151 * x31
    x159 = x6 * (2.0 * x10 * x151 * x17 * x32 * x4 * x6 - x158)
    x160 = x151 * x35
    x161 = 0.5 * x159 + x160
    x162 = x161 * x46
    x163 = 3.0 * x162
    x164 = x157 + x163
    x165 = -2.0 * x10 * x151 * x17 * x32 * x4 * x5 * x6
    x166 = x161 * x5 - x6 * (x151 * x47 + x165)
    x167 = x166 * x3
    x168 = x151 * x37
    x169 = -2.0 * x10 * x151 * x17 * x38 * x4 * x6
    x170 = -x6 * (x168 + x169)
    x171 = x151 * x41
    x172 = 0.5 * x170 + x171
    x173 = x172 * x46
    x174 = x167 + 3.0 * x173
    x175 = x151 * x67
    x176 = x6 * (2.0 * x10 * x151 * x17 * x4 * x6 * x9 - x175)
    x177 = x151 * x71
    x178 = x151 * x6 * (2.0 * x10 * x17 * x4 * x5 * x6 * x9 - x68) + 0.5 * x5 * (
        x176 + 2.0 * x177
    )
    x179 = x178 * x3
    x180 = x155 * x46
    x181 = 3.0 * x180
    x182 = x156 * x52
    x183 = x77 * (x166 - x182)
    x184 = x155 * x3
    x185 = x151 * x5
    x186 = x185 * x9
    x187 = x140 * x186
    x188 = x184 + x187
    x189 = x161 * x52
    x190 = x77 * (x172 - x189)
    x191 = x104 * x89
    x192 = x151 * x191
    x193 = x151 * x34
    x194 = x108 * x193
    x195 = x192 + x194
    x196 = x104 * x195 + x188 * x3 + x190
    x197 = x161 * x3
    x198 = x83 * x89
    x199 = x185 * x198
    x200 = x197 + x199
    x201 = x172 * x3
    x202 = x141 * x151
    x203 = x201 + x202
    x204 = x144 * x18
    x205 = x185 * x204
    x206 = -x77 * (x165 + x205)
    x207 = x192 * x3
    x208 = x195 * x3 + x206 + x207
    x209 = x104 * x208 + x196 * x3 - x4 * (x200 * x52 - x203)
    x210 = x151 * x6 * (2.0 * x10 * x17 * x38 * x4 * x5 * x6 - x61) + x172 * x5
    x211 = x77 * (-x166 * x52 + x210)
    x212 = x164 * x3 + x200 * x92 + x211
    x213 = x151 * x53
    x214 = x6 * (2.0 * x10 * x151 * x17 * x4 * x54 * x6 - x213)
    x215 = x151 * x57
    x216 = 0.5 * x214 + x215
    x217 = x151 * x6 * (2.0 * x10 * x17 * x4 * x5 * x54 * x6 - x122) + x216 * x5
    x218 = x77 * (-x210 * x52 + x217)
    x219 = x174 * x3 + x203 * x92 + x218
    x220 = x151 * x6 * (2.0 * x10 * x114 * x17 * x4 * x6 - x113)
    x221 = x117 * x151 + 0.5 * x220
    x222 = x221 * x46
    x223 = x216 * x46
    x224 = x210 * x3 + 3.0 * x223
    x225 = x77 * (-x216 * x52 + x221)
    x226 = x106 * x38
    x227 = x151 * x226
    x228 = x151 * x56
    x229 = x108 * x228 + x227
    x230 = x104 * x229 + x203 * x3 + x225
    x231 = x77 * (-x172 * x52 + x216)
    x232 = x107 * x151
    x233 = x151 * x40
    x234 = x108 * x233 + x232
    x235 = x104 * x234 + x200 * x3 + x231
    x236 = x144 * x32
    x237 = -x77 * (x151 * x236 + x169)
    x238 = x3**2
    x239 = 12.1633560763699 * x148
    x240 = -x1 * (ax * A[2] + bx * B[2])
    x241 = -x240 - B[2]
    x242 = x241 * x6 * (2.0 * x10 * x17 * x18 * x4 * x6 - x15)
    x243 = 0.5 * x242
    x244 = x241 * x27 + x243
    x245 = x241 * x6 * (2.0 * x10 * x17 * x18 * x4 * x5 * x6 - x16) + x244 * x5
    x246 = x245 * x3
    x247 = x241 * x6 * (2.0 * x10 * x17 * x32 * x4 * x6 - x31)
    x248 = 0.5 * x247
    x249 = x241 * x35 + x248
    x250 = x249 * x46
    x251 = 3.0 * x250
    x252 = x246 + x251
    x253 = -2.0 * x10 * x17 * x241 * x32 * x4 * x5 * x6
    x254 = x249 * x5 - x6 * (x241 * x47 + x253)
    x255 = x254 * x3
    x256 = -2.0 * x10 * x17 * x241 * x38 * x4 * x6
    x257 = -x6 * (x241 * x37 + x256)
    x258 = 0.5 * x257
    x259 = x241 * x41 + x258
    x260 = x259 * x46
    x261 = x255 + 3.0 * x260
    x262 = x241 * x6 * (2.0 * x10 * x17 * x4 * x6 * x9 - x67)
    x263 = 0.5 * x262
    x264 = x241 * x6 * (2.0 * x10 * x17 * x4 * x5 * x6 * x9 - x68) + x5 * (
        x241 * x71 + x263
    )
    x265 = x264 * x3
    x266 = x244 * x46
    x267 = 3.0 * x266
    x268 = x245 * x52
    x269 = x77 * (x254 - x268)
    x270 = x244 * x3
    x271 = x241 * x9
    x272 = x105 * x271 * x87
    x273 = x270 + x272
    x274 = x249 * x52
    x275 = x77 * (x259 - x274)
    x276 = x191 * x241
    x277 = x241 * x34
    x278 = x108 * x277
    x279 = x276 + x278
    x280 = x104 * x279 + x273 * x3 + x275
    x281 = x249 * x3
    x282 = x241 * x90
    x283 = x281 + x282
    x284 = x259 * x3
    x285 = x141 * x241
    x286 = x284 + x285
    x287 = x204 * x241
    x288 = -x77 * (x253 + x287 * x5)
    x289 = x276 * x3
    x290 = x279 * x3 + x288 + x289
    x291 = x104 * x290 + x280 * x3 - x4 * (x283 * x52 - x286)
    x292 = x241 * x6 * (2.0 * x10 * x17 * x38 * x4 * x5 * x6 - x61) + x259 * x5
    x293 = x77 * (-x254 * x52 + x292)
    x294 = x252 * x3 + x283 * x92 + x293
    x295 = x241 * x6 * (2.0 * x10 * x17 * x4 * x54 * x6 - x53)
    x296 = 0.5 * x295
    x297 = x241 * x57 + x296
    x298 = x241 * x6 * (2.0 * x10 * x17 * x4 * x5 * x54 * x6 - x122) + x297 * x5
    x299 = x77 * (-x292 * x52 + x298)
    x300 = x261 * x3 + x286 * x92 + x299
    x301 = x241 * x6 * (2.0 * x10 * x114 * x17 * x4 * x6 - x113)
    x302 = 0.5 * x301
    x303 = x117 * x241 + x302
    x304 = x303 * x46
    x305 = x297 * x46
    x306 = x292 * x3 + 3.0 * x305
    x307 = x77 * (-x297 * x52 + x303)
    x308 = x226 * x241
    x309 = x241 * x56
    x310 = x108 * x309 + x308
    x311 = x104 * x310 + x286 * x3 + x307
    x312 = x77 * (-x259 * x52 + x297)
    x313 = x107 * x241
    x314 = x241 * x40
    x315 = x108 * x314 + x313
    x316 = x104 * x315 + x283 * x3 + x312
    x317 = -x77 * (x236 * x241 + x256)
    x318 = x151**2
    x319 = x25 * x318
    x320 = x21 * x319
    x321 = x20 + x320
    x322 = x318 * x34 + x33
    x323 = x318 * x40 + x39
    x324 = -x322 * x7 + x323
    x325 = x19 * x324 + x22 * x321
    x326 = x3 * x325
    x327 = x322 * x46
    x328 = x327 * x5
    x329 = x326 + 2.0 * x328
    x330 = x318 * x56 + x55
    x331 = -x323 * x7 + x330
    x332 = x19 * x331 + x22 * x322
    x333 = x3 * x332
    x334 = x323 * x5
    x335 = x104 * x334 + x333
    x336 = x319 * x70
    x337 = x336 + x69
    x338 = -x321 * x7 + x322
    x339 = x19 * x338 + x22 * x337
    x340 = x3 * x339
    x341 = x321 * x5
    x342 = x104 * x341
    x343 = x325 * x52
    x344 = x77 * (x332 - x343)
    x345 = x108 * x321
    x346 = x327 + x345
    x347 = x5 * x52
    x348 = x77 * (-x322 * x347 + x323 * x5)
    x349 = x3 * x327
    x350 = x3 * x346 + x348 + x349
    x351 = x323 * x46
    x352 = x108 * x322
    x353 = x351 + x352
    x354 = x330 * x46
    x355 = x3 * x334
    x356 = x354 + x355
    x357 = x238 * x322
    x358 = x323 * x52
    x359 = x77 * (x330 - x358)
    x360 = x357 + x359
    x361 = x3 * x350 + x360 * x46 - x4 * (x353 * x52 - x356)
    x362 = x115 + x116 * x318
    x363 = -x330 * x7 + x362
    x364 = x19 * x363 + x22 * x323
    x365 = x77 * (-x332 * x52 + x364)
    x366 = x104 * x353 + x3 * x329 + x365
    x367 = x124 + x125 * x318 - x362 * x7
    x368 = x19 * x367 + x22 * x330
    x369 = x77 * (-x364 * x52 + x368)
    x370 = x104 * x356 + x3 * x335 + x369
    x371 = x354 * x5
    x372 = x3 * x364 + 2.0 * x371
    x373 = x362 * x5
    x374 = x5 * x77 * (-x330 * x52 + x362)
    x375 = x3 * x354 + x3 * x356 + x374
    x376 = x77 * (x330 * x5 - x334 * x52)
    x377 = x3 * x351 + x3 * x353 + x376
    x378 = 15.7028251725905 * x148
    x379 = x241 * x6 * (2.0 * x10 * x151 * x17 * x18 * x4 * x6 - x152)
    x380 = x154 * x241 + 0.5 * x379
    x381 = x186 * x241
    x382 = x140 * x381
    x383 = x3 * x380 + x382
    x384 = -2.0 * x10 * x151 * x17 * x241 * x32 * x4 * x6
    x385 = -x6 * (x158 * x241 + x384)
    x386 = x160 * x241 + 0.5 * x385
    x387 = x241 * x83
    x388 = x387 * x89
    x389 = x185 * x388
    x390 = x3 * x386 + x389
    x391 = x241 * x6 * (2.0 * x10 * x151 * x17 * x4 * x6 * x9 - x175)
    x392 = x177 * x241 + 0.5 * x391
    x393 = x140 * x185 * x21 * x241
    x394 = x77 * (-x380 * x52 + x386)
    x395 = x151 * x241
    x396 = x395 * x9
    x397 = x106 * x396
    x398 = x21 * x25 * x395
    x399 = x108 * x398 + x397
    x400 = x77 * (2.0 * x10 * x151 * x17 * x18 * x241 * x4 * x5 * x6 - x144 * x381)
    x401 = x3 * x397 + x3 * x399 + x400
    x402 = x193 * x241
    x403 = x108 * x402 + x192 * x241
    x404 = x233 * x241
    x405 = x108 * x404 + x232 * x241
    x406 = x151 * x287
    x407 = -x77 * (x384 + x406)
    x408 = x238 * x402 + x407
    x409 = x3 * x401 - x4 * (x403 * x52 - x405) + x408 * x46
    x410 = x241 * x6 * (2.0 * x10 * x151 * x17 * x38 * x4 * x6 - x168)
    x411 = x171 * x241 + 0.5 * x410
    x412 = x77 * (-x386 * x52 + x411)
    x413 = x104 * x403 + x3 * x383 + x412
    x414 = x241 * x6 * (2.0 * x10 * x151 * x17 * x4 * x54 * x6 - x213)
    x415 = x215 * x241 + 0.5 * x414
    x416 = x77 * (-x411 * x52 + x415)
    x417 = x104 * x405 + x3 * x390 + x416
    x418 = x141 * x395
    x419 = x3 * x411 + x418
    x420 = x135 * x387
    x421 = x185 * x420
    x422 = x77 * (2.0 * x10 * x151 * x17 * x241 * x38 * x4 * x5 * x6 - x145 * x395)
    x423 = x147 * x395 + x3 * x405 + x422
    x424 = x241 * x77 * (2.0 * x10 * x151 * x17 * x32 * x4 * x5 * x6 - x205)
    x425 = x207 * x241 + x3 * x403 + x424
    x426 = 27.1980910212982 * x148
    x427 = x241**2
    x428 = x25 * x427
    x429 = x20 + x21 * x428
    x430 = x33 + x34 * x427
    x431 = x39 + x40 * x427
    x432 = -x430 * x7 + x431
    x433 = x19 * x432
    x434 = x22 * x429 + x433
    x435 = x3 * x434
    x436 = x430 * x46
    x437 = x436 * x5
    x438 = x435 + 2.0 * x437
    x439 = x427 * x56 + x55
    x440 = -x431 * x7 + x439
    x441 = x19 * x440
    x442 = x22 * x430 + x441
    x443 = x3 * x442
    x444 = x431 * x5
    x445 = x104 * x444 + x443
    x446 = x428 * x70 + x69
    x447 = -x429 * x7 + x430
    x448 = x19 * x447
    x449 = x22 * x446 + x448
    x450 = x3 * x449
    x451 = x429 * x5
    x452 = x434 * x52
    x453 = x77 * (x442 - x452)
    x454 = x108 * x429
    x455 = x436 + x454
    x456 = x77 * (-x347 * x430 + x431 * x5)
    x457 = x3 * x436
    x458 = x3 * x455 + x456 + x457
    x459 = x431 * x46
    x460 = x108 * x430
    x461 = x459 + x460
    x462 = x439 * x46
    x463 = x3 * x444
    x464 = x462 + x463
    x465 = x238 * x430
    x466 = x431 * x52
    x467 = x77 * (x439 - x466)
    x468 = x465 + x467
    x469 = x3 * x458 - x4 * (x461 * x52 - x464) + x46 * x468
    x470 = x115 + x116 * x427
    x471 = -x439 * x7 + x470
    x472 = x19 * x471
    x473 = x22 * x431 + x472
    x474 = x77 * (-x442 * x52 + x473)
    x475 = x104 * x461 + x3 * x438 + x474
    x476 = x124 + x125 * x427 - x470 * x7
    x477 = x19 * x476
    x478 = x22 * x439 + x477
    x479 = x77 * (-x473 * x52 + x478)
    x480 = x104 * x464 + x3 * x445 + x479
    x481 = x462 * x5
    x482 = x3 * x473 + 2.0 * x481
    x483 = x470 * x5
    x484 = x439 * x5
    x485 = x77 * (x470 * x5 - x484 * x52)
    x486 = x3 * x462 + x3 * x464 + x485
    x487 = x77 * (x439 * x5 - x444 * x52)
    x488 = x3 * x459
    x489 = x3 * x461 + x487 + x488
    x490 = x151 * x323 + x170
    x491 = x151 * x322 + x159
    x492 = x491 * x52
    x493 = x3 * x492
    x494 = x151 * x321 + x153
    x495 = x238 * x494
    x496 = x77 * (x490 - x492)
    x497 = x495 + x496
    x498 = x3 * x497 + x4 * (x3 * x490 - x493)
    x499 = x151 * x330 + x214
    x500 = x77 * (-x490 * x52 + x499)
    x501 = x238 * x491 + x500
    x502 = x151 * x362 + x220
    x503 = x77 * (-x499 * x52 + x502)
    x504 = x238 * x490 + x503
    x505 = x46 * x491
    x506 = x108 * x494
    x507 = x505 + x506
    x508 = x46 * x490
    x509 = x491 * x5
    x510 = x3 * x509
    x511 = x508 + x510
    x512 = x46 * x494
    x513 = x151 * x337 + x176
    x514 = x108 * x513
    x515 = x77 * (-x347 * x494 + x491 * x5)
    x516 = x3 * x512
    x517 = x77 * (x490 * x5 - x509 * x52)
    x518 = x3 * x505
    x519 = x3 * x507 + x517 + x518
    x520 = x490 * x5
    x521 = x77 * (x499 * x5 - x52 * x520)
    x522 = x3 * x490
    x523 = x46 * x522
    x524 = x3 * x511 + x521 + x523
    x525 = x46 * x502
    x526 = x499 * x5
    x527 = x46 * x499
    x528 = x5 * x522 + x527
    x529 = x258 + x314 * x318
    x530 = x248 + x277 * x318
    x531 = x52 * x530
    x532 = x241 * x320 + x243
    x533 = x77 * (x529 - x531)
    x534 = x238 * x532 + x533
    x535 = x3 * (x4 * (x529 - x531) + x534)
    x536 = x296 + x309 * x318
    x537 = x77 * (-x52 * x529 + x536)
    x538 = x238 * x530 + x537
    x539 = x116 * x241 * x318 + x302
    x540 = x77 * (-x52 * x536 + x539)
    x541 = x238 * x529 + x540
    x542 = x46 * x530
    x543 = x108 * x532 + x542
    x544 = x46 * x529
    x545 = x5 * x530
    x546 = x3 * x545 + x544
    x547 = x46 * x532
    x548 = x241 * x336 + x263
    x549 = x77 * (-x347 * x532 + x5 * x530)
    x550 = x77 * (x5 * x529 - x52 * x545)
    x551 = x3 * x542 + x3 * x543 + x550
    x552 = x5 * x529
    x553 = x77 * (x5 * x536 - x52 * x552)
    x554 = x3 * x529
    x555 = x3 * x546 + x46 * x554 + x553
    x556 = x46 * x539
    x557 = x5 * x536
    x558 = x46 * x536
    x559 = x5 * x554 + x558
    x560 = x151 * x430
    x561 = x52 * x560
    x562 = x77 * (x151 * x431 - x561)
    x563 = x151 * x429
    x564 = x238 * x563 + x562
    x565 = x3 * (x4 * (x151 * x431 - x561) + x564)
    x566 = x151 * x431
    x567 = x77 * (x151 * x439 - x52 * x566)
    x568 = x151 * x465 + x567
    x569 = x151 * x77 * (-x439 * x52 + x470)
    x570 = x238 * x566 + x569
    x571 = x151 * x436
    x572 = x151 * x454 + x571
    x573 = x151 * x459
    x574 = x151 * x460 + x573
    x575 = x46 * x563
    x576 = x151 * x446
    x577 = x77 * (x151 * x430 * x5 - x347 * x563)
    x578 = x5 * x560
    x579 = x77 * (x151 * x431 * x5 - x52 * x578)
    x580 = x151 * x457 + x3 * x572 + x579
    x581 = x151 * x444
    x582 = x77 * (x151 * x439 * x5 - x52 * x581)
    x583 = x151 * x488 + x3 * x574 + x582
    x584 = x151 * x470
    x585 = x151 * x462
    x586 = x151 * x463 + x585
    x587 = x151 * x484
    x588 = x241 * x431 + x257
    x589 = x241 * x430 + x247
    x590 = x52 * x589
    x591 = x3 * x590
    x592 = x241 * x429 + x242
    x593 = x238 * x592
    x594 = x77 * (x588 - x590)
    x595 = x593 + x594
    x596 = x3 * x595 + x4 * (x3 * x588 - x591)
    x597 = x238 * x589
    x598 = x241 * x439 + x295
    x599 = x77 * (-x52 * x588 + x598)
    x600 = x597 + x599
    x601 = x241 * x470 + x301
    x602 = x77 * (-x52 * x598 + x601)
    x603 = x238 * x588 + x602
    x604 = x46 * x589
    x605 = x108 * x592
    x606 = x604 + x605
    x607 = x46 * x588
    x608 = x5 * x589
    x609 = x3 * x608
    x610 = x607 + x609
    x611 = x46 * x592
    x612 = x241 * x446 + x262
    x613 = x108 * x612
    x614 = x77 * (-x347 * x592 + x5 * x589)
    x615 = x3 * x611
    x616 = x77 * (x5 * x588 - x52 * x608)
    x617 = x3 * x604
    x618 = x3 * x606 + x616 + x617
    x619 = x5 * x588
    x620 = x77 * (x5 * x598 - x52 * x619)
    x621 = x3 * x588
    x622 = x46 * x621
    x623 = x3 * x610 + x620 + x622
    x624 = x46 * x601
    x625 = x5 * x598
    x626 = x46 * x598
    x627 = x5 * x621 + x626
    x628 = x151 * x491 + x331 * x43
    x629 = x151 * x494 + x324 * x43
    x630 = x52 * x629
    x631 = x3 * x630
    x632 = x151 * x513 + x338 * x43
    x633 = x238 * x632
    x634 = x77 * (x628 - x630)
    x635 = x238 * x629
    x636 = x151 * x490 + x363 * x43
    x637 = x77 * (-x52 * x628 + x636)
    x638 = x635 + x637
    x639 = x238 * x628
    x640 = x151 * x499 + x367 * x43
    x641 = x77 * (-x52 * x636 + x640)
    x642 = x639 + x641
    x643 = x151 * x530 + x385
    x644 = x151 * x532 + x379
    x645 = x52 * x644
    x646 = x151 * x548 + x391
    x647 = x77 * (x643 - x645)
    x648 = x151 * x529 + x410
    x649 = x77 * (-x52 * x643 + x648)
    x650 = x238 * x644 + x649
    x651 = x151 * x536 + x414
    x652 = x77 * (-x52 * x648 + x651)
    x653 = x238 * x643 + x652
    x654 = x318 * x430 + x441
    x655 = x318 * x429 + x433
    x656 = x52 * x655
    x657 = x318 * x446 + x448
    x658 = x77 * (x654 - x656)
    x659 = x318 * x431 + x472
    x660 = x77 * (-x52 * x654 + x659)
    x661 = x238 * x655 + x660
    x662 = x318 * x439 + x477
    x663 = x77 * (-x52 * x659 + x662)
    x664 = x238 * x654 + x663
    x665 = x151 * x592
    x666 = x52 * x665
    x667 = x77 * (x151 * x589 - x666)
    x668 = x151 * x612
    x669 = x151 * x589
    x670 = x77 * (x151 * x588 - x52 * x669)
    x671 = x151 * x593 + x670
    x672 = x151 * x588
    x673 = x77 * (x151 * x598 - x52 * x672)
    x674 = x151 * x597 + x673
    x675 = x241 * x589 + x43 * x440
    x676 = x241 * x592 + x43 * x432
    x677 = x52 * x676
    x678 = x3 * x677
    x679 = x241 * x612 + x43 * x447
    x680 = x238 * x679
    x681 = x77 * (x675 - x677)
    x682 = x238 * x676
    x683 = x241 * x588 + x43 * x471
    x684 = x77 * (-x52 * x675 + x683)
    x685 = x682 + x684
    x686 = x238 * x675
    x687 = x241 * x598 + x43 * x476
    x688 = x77 * (-x52 * x683 + x687)
    x689 = x686 + x688
    x690 = -x150 - A[1]
    x691 = x45 * x690
    x692 = x50 * x690
    x693 = x691 + x692
    x694 = x60 * x690
    x695 = x63 * x690
    x696 = x694 + x695 * x83
    x697 = x690 * x73
    x698 = x690 * x75
    x699 = x690 * x76
    x700 = x4 * (x59 * x690 - x699)
    x701 = x690 * x79
    x702 = x690 * x81
    x703 = x701 + x702
    x704 = x690 * x84
    x705 = x4 * (x63 * x690 - x704)
    x706 = x690 * x86
    x707 = x5 * x690
    x708 = x198 * x707
    x709 = x706 + x708
    x710 = x3 * x703 + 0.5 * x705 + x709 * x92
    x711 = x690 * (x94 + x96)
    x712 = x58 * x690
    x713 = x690 * x98 + x712 * x92
    x714 = x4 * x690 * (-x102 + x58)
    x715 = x107 * x690
    x716 = x4 * x690 * (-x119 * x52 + x126)
    x717 = x4 * x690 * (x119 - x52 * x59)
    x718 = 13.7919494082862 * x148
    x719 = x156 * x690
    x720 = x49 + x719
    x721 = x161 * x690
    x722 = x721 + x95
    x723 = x3 * x720 + x722 * x92
    x724 = x166 * x690
    x725 = x64 + x724
    x726 = x172 * x690
    x727 = x726 + x99
    x728 = x3 * x725 + x727 * x92
    x729 = x178 * x690
    x730 = x729 + x74
    x731 = x155 * x690
    x732 = x731 + x80
    x733 = x4 * (-x52 * x720 + x725)
    x734 = x191 * x5
    x735 = x193 * x707
    x736 = x734 + x735
    x737 = x104 * x736
    x738 = x3 * x732 + x737
    x739 = x4 * (-x52 * x722 + x727)
    x740 = x233 * x690
    x741 = x107 + x740
    x742 = x46 * x741
    x743 = x3 * x736 + x742
    x744 = x104 * x743 + x3 * x738 + 0.5 * x739
    x745 = x106 * x139
    x746 = x233 * x707 + x745
    x747 = x104 * x746
    x748 = x3 * x722 + x747
    x749 = x226 * x5
    x750 = x228 * x707 + x749
    x751 = x104 * x750
    x752 = x3 * x727 + x751
    x753 = x4 * (-x52 * x746 + x750)
    x754 = x132 + x210 * x690
    x755 = x4 * (x131 + x217 * x690 - x52 * x754)
    x756 = x4 * (-x52 * x725 + x754)
    x757 = 36.4900682291097 * x148
    x758 = x690 * (x246 + x251)
    x759 = x259 * x690
    x760 = x255 * x690 + x759 * x92
    x761 = x4 * x690 * (x254 - x268)
    x762 = x140 * x271 * x707 + x270 * x690
    x763 = x4 * x690 * (x259 - x274)
    x764 = x276 * x690
    x765 = x278 * x690 + x764
    x766 = x104 * x765 + x3 * x762 + 0.5 * x763
    x767 = x281 * x690 + x388 * x707
    x768 = x140 * x32
    x769 = x241 * x768
    x770 = x284 * x690 + x707 * x769
    x771 = x4 * (2.0 * x10 * x17 * x241 * x32 * x4 * x5 * x6 * x690 - x287 * x707)
    x772 = x4 * x690 * (-x292 * x52 + x298)
    x773 = x4 * x690 * (-x254 * x52 + x292)
    x774 = x325 * x690
    x775 = 2.0 * x162
    x776 = x774 + x775
    x777 = x322 * x690
    x778 = x5 * x777
    x779 = x199 + x778
    x780 = x104 * x779 + x3 * x776
    x781 = x332 * x690
    x782 = 2.0 * x173 + x781
    x783 = x334 * x690
    x784 = x202 + x783
    x785 = x104 * x784 + x3 * x782
    x786 = x339 * x690
    x787 = 2.0 * x180
    x788 = x786 + x787
    x789 = x321 * x707
    x790 = x187 + x789
    x791 = x4 * (-x52 * x776 + x782)
    x792 = x151 * x83
    x793 = x792 * x89
    x794 = x777 + x793
    x795 = x46 * x794
    x796 = x3 * x790 + x795
    x797 = x4 * (-x52 * x779 + x784)
    x798 = x3 * x795 + x3 * x796 + 0.5 * x797
    x799 = x151 * x768
    x800 = x323 * x690 + x799
    x801 = x46 * x800
    x802 = x3 * x779 + x801
    x803 = x135 * x792 + x330 * x690
    x804 = x46 * x803
    x805 = x3 * x784 + x804
    x806 = x4 * (-x52 * x800 + x803)
    x807 = 2.0 * x223 + x364 * x690
    x808 = x4 * (2.0 * x222 + x368 * x690 - x52 * x807)
    x809 = x4 * (-x52 * x782 + x807)
    x810 = 47.1084755177714 * x148
    x811 = x250 + x380 * x690
    x812 = x276 * x5 + x402 * x707
    x813 = x104 * x812
    x814 = x3 * x811 + x813
    x815 = x260 + x386 * x690
    x816 = x241 * x745 + x404 * x707
    x817 = x104 * x816
    x818 = x3 * x815 + x817
    x819 = x266 + x392 * x690
    x820 = x106 * x271
    x821 = x398 * x707 + x5 * x820
    x822 = x104 * x821
    x823 = x4 * (-x52 * x811 + x815)
    x824 = x276 + x402 * x690
    x825 = x46 * x824
    x826 = x3 * x821 + x825
    x827 = x4 * (-x52 * x812 + x816)
    x828 = x3 * x825 + x3 * x826 + 0.5 * x827
    x829 = x313 + x404 * x690
    x830 = x46 * x829
    x831 = x3 * x812 + x830
    x832 = x228 * x241
    x833 = x308 + x690 * x832
    x834 = x46 * x833
    x835 = x3 * x816 + x834
    x836 = x4 * (-x52 * x829 + x833)
    x837 = x305 + x411 * x690
    x838 = x4 * (x304 + x415 * x690 - x52 * x837)
    x839 = x4 * (-x52 * x815 + x837)
    x840 = 81.5942730638946 * x148
    x841 = x436 * x707
    x842 = x435 * x690 + 2.0 * x841
    x843 = x444 * x690
    x844 = x104 * x843 + x443 * x690
    x845 = x429 * x707
    x846 = x4 * x690 * (x442 - x452)
    x847 = x436 * x690
    x848 = x454 * x690 + x847
    x849 = x52 * x707
    x850 = x4 * (-x430 * x849 + x431 * x5 * x690)
    x851 = x3 * x848 + x457 * x690 + 0.5 * x850
    x852 = x459 * x690
    x853 = x460 * x690 + x852
    x854 = x462 * x690
    x855 = x463 * x690 + x854
    x856 = x4 * x690 * (x439 - x466)
    x857 = x4 * x690 * (-x473 * x52 + x478)
    x858 = x4 * x690 * (-x442 * x52 + x473)
    x859 = x490 * x690
    x860 = 3.0 * x354
    x861 = x859 + x860
    x862 = x491 * x690
    x863 = 3.0 * x351
    x864 = x862 + x863
    x865 = x52 * x864
    x866 = x494 * x690
    x867 = 3.0 * x327
    x868 = x866 + x867
    x869 = x4 * (x861 - x865)
    x870 = x238 * x868 + 0.5 * x869
    x871 = x46 * x864
    x872 = x5 * x866
    x873 = x5 * x867 + x872
    x874 = x3 * x873 + x871
    x875 = x46 * x861
    x876 = x509 * x690
    x877 = x334 * x92 + x876
    x878 = x3 * x877 + x875
    x879 = x46 * x868
    x880 = x513 * x707
    x881 = x341 * x92
    x882 = x880 + x881
    x883 = x4 * (-x52 * x873 + x877)
    x884 = x5 * x860 + x520 * x690
    x885 = x4 * (x373 * x92 - x52 * x884 + x526 * x690)
    x886 = x4 * (-x52 * x877 + x884)
    x887 = x395 * x768
    x888 = x529 * x690 + x887
    x889 = x198 * x395
    x890 = x530 * x690 + x889
    x891 = x52 * x890
    x892 = x532 * x690
    x893 = x140 * x396
    x894 = x892 + x893
    x895 = x4 * (x888 - x891)
    x896 = x238 * x894 + 0.5 * x895
    x897 = x46 * x890
    x898 = x382 + x5 * x892
    x899 = x3 * x898 + x897
    x900 = x46 * x888
    x901 = x389 + x545 * x690
    x902 = x3 * x901 + x900
    x903 = x46 * x894
    x904 = x393 + x548 * x707
    x905 = x4 * (-x52 * x898 + x901)
    x906 = x418 + x552 * x690
    x907 = x4 * (x421 - x52 * x906 + x557 * x690)
    x908 = x4 * (-x52 * x901 + x906)
    x909 = x462 + x566 * x690
    x910 = x459 + x560 * x690
    x911 = x52 * x910
    x912 = x436 + x563 * x690
    x913 = x4 * (x909 - x911)
    x914 = x238 * x912 + 0.5 * x913
    x915 = x46 * x910
    x916 = x437 + x563 * x707
    x917 = x3 * x916 + x915
    x918 = x46 * x909
    x919 = x444 * x46 + x578 * x690
    x920 = x3 * x919 + x918
    x921 = x46 * x912
    x922 = x451 * x46 + x576 * x707
    x923 = x4 * (-x52 * x916 + x919)
    x924 = x481 + x581 * x690
    x925 = x4 * (x46 * x483 - x52 * x924 + x587 * x690)
    x926 = x4 * (-x52 * x919 + x924)
    x927 = x4 * x690 * (x588 - x590)
    x928 = x593 * x690 + 0.5 * x927
    x929 = x604 * x690
    x930 = x605 * x690 + x929
    x931 = x607 * x690
    x932 = x609 * x690 + x931
    x933 = x611 * x690
    x934 = x4 * (x5 * x589 * x690 - x592 * x849)
    x935 = x4 * x690 * (x5 * x598 - x52 * x619)
    x936 = x4 * x690 * (x5 * x588 - x52 * x608)
    x937 = x628 * x690
    x938 = 4.0 * x508
    x939 = x937 + x938
    x940 = x629 * x690
    x941 = 4.0 * x505
    x942 = x940 + x941
    x943 = x52 * x942
    x944 = x632 * x690
    x945 = 4.0 * x512
    x946 = x944 + x945
    x947 = x4 * (x939 - x943)
    x948 = 4.0 * x527 + x636 * x690
    x949 = x4 * (-x52 * x948 + 4.0 * x525 + x640 * x690)
    x950 = x4 * (-x52 * x939 + x948)
    x951 = 3.0 * x544 + x643 * x690
    x952 = 3.0 * x542 + x644 * x690
    x953 = x52 * x952
    x954 = 3.0 * x547 + x646 * x690
    x955 = x4 * (x951 - x953)
    x956 = 3.0 * x558 + x648 * x690
    x957 = x4 * (-x52 * x956 + 3.0 * x556 + x651 * x690)
    x958 = x4 * (-x52 * x951 + x956)
    x959 = 2.0 * x573 + x654 * x690
    x960 = 2.0 * x571 + x655 * x690
    x961 = x52 * x960
    x962 = x104 * x563 + x657 * x690
    x963 = x4 * (x959 - x961)
    x964 = 2.0 * x585 + x659 * x690
    x965 = x4 * (x104 * x584 - x52 * x964 + x662 * x690)
    x966 = x4 * (-x52 * x959 + x964)
    x967 = x607 + x669 * x690
    x968 = x604 + x665 * x690
    x969 = x52 * x968
    x970 = x611 + x668 * x690
    x971 = x4 * (x967 - x969)
    x972 = x151 * x598
    x973 = x626 + x672 * x690
    x974 = x4 * (-x52 * x973 + x624 + x690 * x972)
    x975 = x4 * (-x52 * x967 + x973)
    x976 = x4 * x690 * (x675 - x677)
    x977 = x4 * x690 * (-x52 * x683 + x687)
    x978 = x4 * x690 * (-x52 * x675 + x683)
    x979 = -x240 - A[2]
    x980 = x979 * (x45 + x50)
    x981 = x979 * (x60 + x65)
    x982 = x4 * x979 * (x59 - x76)
    x983 = 0.5 * x982
    x984 = x979 * (x79 + x81)
    x985 = x4 * x979 * (x63 - x84)
    x986 = 0.5 * x985
    x987 = x5 * x979
    x988 = x198 * x987 + x86 * x979
    x989 = x3 * x984 + x92 * x988 + x986
    x990 = x979 * (x94 + x96)
    x991 = x979 * (x100 + x98)
    x992 = x4 * x979 * (-x102 + x58)
    x993 = 0.5 * x992
    x994 = x107 * x979
    x995 = x4 * x979 * (-x119 * x52 + x126)
    x996 = 0.5 * x995
    x997 = x59 * x979
    x998 = x4 * (x119 * x979 - x52 * x997)
    x999 = 0.5 * x998
    x1000 = x979 * (x157 + x163)
    x1001 = x172 * x979
    x1002 = x1001 * x92 + x167 * x979
    x1003 = x4 * x979 * (x166 - x182)
    x1004 = 0.5 * x1003
    x1005 = x151 * x987
    x1006 = x1005 * x9
    x1007 = x1006 * x140
    x1008 = x1007 + x184 * x979
    x1009 = x4 * x979 * (x172 - x189)
    x1010 = 0.5 * x1009
    x1011 = x192 * x979
    x1012 = x1011 + x194 * x979
    x1013 = x1008 * x3 + x1010 + x1012 * x104
    x1014 = x793 * x987
    x1015 = x1014 + x197 * x979
    x1016 = x799 * x987
    x1017 = x1016 + x201 * x979
    x1018 = x4 * (2.0 * x10 * x151 * x17 * x32 * x4 * x5 * x6 * x979 - x1005 * x204)
    x1019 = 0.5 * x1018
    x1020 = x4 * x979 * (-x210 * x52 + x217)
    x1021 = 0.5 * x1020
    x1022 = x4 * x979 * (-x166 * x52 + x210)
    x1023 = 0.5 * x1022
    x1024 = x245 * x979 + x49
    x1025 = x1024 * x3
    x1026 = x249 * x979 + x95
    x1027 = x1026 * x46
    x1028 = 3.0 * x1027
    x1029 = x1025 + x1028
    x1030 = x254 * x979 + x64
    x1031 = x1030 * x3
    x1032 = x259 * x979 + x99
    x1033 = x1032 * x46
    x1034 = x1031 + 3.0 * x1033
    x1035 = x264 * x979 + x74
    x1036 = x1035 * x3
    x1037 = x244 * x979 + x80
    x1038 = x1037 * x46
    x1039 = 3.0 * x1038
    x1040 = x1024 * x52
    x1041 = x4 * (x1030 - x1040)
    x1042 = 0.5 * x1041
    x1043 = x1037 * x3
    x1044 = x277 * x987 + x734
    x1045 = x1044 * x46
    x1046 = 2.0 * x1045
    x1047 = x1043 + x1046
    x1048 = x1026 * x52
    x1049 = x4 * (x1032 - x1048)
    x1050 = 0.5 * x1049
    x1051 = x107 + x314 * x979
    x1052 = x1051 * x46
    x1053 = x1044 * x3
    x1054 = x1052 + x1053
    x1055 = x104 * x1054 + x1047 * x3 + x1050
    x1056 = x314 * x987 + x745
    x1057 = x1056 * x46
    x1058 = 2.0 * x1057
    x1059 = x1026 * x3 + x1058
    x1060 = x309 * x987 + x749
    x1061 = x1060 * x46
    x1062 = 2.0 * x1061
    x1063 = x1032 * x3 + x1062
    x1064 = x4 * (-x1056 * x52 + x1060)
    x1065 = 0.5 * x1064
    x1066 = x132 + x292 * x979
    x1067 = x4 * (-x1066 * x52 + x131 + x298 * x979)
    x1068 = 0.5 * x1067
    x1069 = x4 * (-x1030 * x52 + x1066)
    x1070 = 0.5 * x1069
    x1071 = x327 * x979
    x1072 = x1071 * x5
    x1073 = 2.0 * x1072 + x326 * x979
    x1074 = x334 * x979
    x1075 = x104 * x1074 + x333 * x979
    x1076 = x4 * x979 * (x332 - x343)
    x1077 = 0.5 * x1076
    x1078 = x1071 + x345 * x979
    x1079 = x52 * x987
    x1080 = x4 * (-x1079 * x322 + x323 * x5 * x979)
    x1081 = 0.5 * x1080
    x1082 = x1078 * x3 + x1081 + x349 * x979
    x1083 = x351 * x979
    x1084 = x1083 + x352 * x979
    x1085 = x354 * x979
    x1086 = x1085 + x355 * x979
    x1087 = x4 * x979 * (x330 - x358)
    x1088 = 0.5 * x1087
    x1089 = x4 * x979 * (-x364 * x52 + x368)
    x1090 = 0.5 * x1089
    x1091 = x4 * x979 * (-x332 * x52 + x364)
    x1092 = 0.5 * x1091
    x1093 = x162 + x380 * x979
    x1094 = x185 * x191 + x402 * x987
    x1095 = x104 * x1094
    x1096 = x1093 * x3 + x1095
    x1097 = x173 + x386 * x979
    x1098 = x151 * x745 + x404 * x987
    x1099 = x104 * x1098
    x1100 = x1097 * x3 + x1099
    x1101 = x180 + x392 * x979
    x1102 = x106 * x186 + x398 * x987
    x1103 = x104 * x1102
    x1104 = x4 * (-x1093 * x52 + x1097)
    x1105 = 0.5 * x1104
    x1106 = x192 + x402 * x979
    x1107 = x1106 * x46
    x1108 = x1102 * x3 + x1107
    x1109 = x4 * (-x1094 * x52 + x1098)
    x1110 = 0.5 * x1109
    x1111 = x1107 * x3 + x1108 * x3 + x1110
    x1112 = x232 + x404 * x979
    x1113 = x1112 * x46
    x1114 = x1094 * x3 + x1113
    x1115 = x227 + x832 * x979
    x1116 = x1115 * x46
    x1117 = x1098 * x3 + x1116
    x1118 = x4 * (-x1112 * x52 + x1115)
    x1119 = 0.5 * x1118
    x1120 = x223 + x411 * x979
    x1121 = x4 * (-x1120 * x52 + x222 + x415 * x979)
    x1122 = 0.5 * x1121
    x1123 = x4 * (-x1097 * x52 + x1120)
    x1124 = 0.5 * x1123
    x1125 = 2.0 * x250 + x434 * x979
    x1126 = x1125 * x3
    x1127 = x430 * x979
    x1128 = x1127 * x5
    x1129 = x1128 + x282
    x1130 = x1129 * x46
    x1131 = 2.0 * x1130
    x1132 = x1126 + x1131
    x1133 = 2.0 * x260 + x442 * x979
    x1134 = x1133 * x3
    x1135 = x285 + x444 * x979
    x1136 = x1135 * x46
    x1137 = x1134 + 2.0 * x1136
    x1138 = 2.0 * x266 + x449 * x979
    x1139 = x1138 * x3
    x1140 = x272 + x451 * x979
    x1141 = x1140 * x46
    x1142 = 2.0 * x1141
    x1143 = x1125 * x52
    x1144 = x4 * (x1133 - x1143)
    x1145 = 0.5 * x1144
    x1146 = x1127 + x388
    x1147 = x1146 * x46
    x1148 = x1140 * x3
    x1149 = x1147 + x1148
    x1150 = x1129 * x52
    x1151 = x4 * (x1135 - x1150)
    x1152 = 0.5 * x1151
    x1153 = x1147 * x3
    x1154 = x1149 * x3 + x1152 + x1153
    x1155 = x431 * x979 + x769
    x1156 = x1155 * x46
    x1157 = x1129 * x3 + x1156
    x1158 = x420 + x439 * x979
    x1159 = x1158 * x46
    x1160 = x1135 * x3 + x1159
    x1161 = x4 * (-x1155 * x52 + x1158)
    x1162 = 0.5 * x1161
    x1163 = 2.0 * x305 + x473 * x979
    x1164 = x4 * (-x1163 * x52 + 2.0 * x304 + x478 * x979)
    x1165 = 0.5 * x1164
    x1166 = x4 * (-x1133 * x52 + x1163)
    x1167 = 0.5 * x1166
    x1168 = x4 * x979 * (x490 - x492)
    x1169 = 0.5 * x1168
    x1170 = x1169 + x495 * x979
    x1171 = x505 * x979
    x1172 = x1171 + x506 * x979
    x1173 = x508 * x979
    x1174 = x1173 + x510 * x979
    x1175 = x512 * x979
    x1176 = x4 * (-x1079 * x494 + x491 * x5 * x979)
    x1177 = 0.5 * x1176
    x1178 = x4 * x979 * (x499 * x5 - x52 * x520)
    x1179 = 0.5 * x1178
    x1180 = x4 * x979 * (x490 * x5 - x509 * x52)
    x1181 = 0.5 * x1180
    x1182 = x354 + x529 * x979
    x1183 = x351 + x530 * x979
    x1184 = x1183 * x52
    x1185 = x532 * x979
    x1186 = x1185 + x327
    x1187 = x4 * (x1182 - x1184)
    x1188 = 0.5 * x1187
    x1189 = x1186 * x238 + x1188
    x1190 = x1183 * x46
    x1191 = x1185 * x5 + x328
    x1192 = x1190 + x1191 * x3
    x1193 = x1182 * x46
    x1194 = x334 * x46 + x545 * x979
    x1195 = x1193 + x1194 * x3
    x1196 = x1186 * x46
    x1197 = x341 * x46
    x1198 = x1197 + x548 * x987
    x1199 = x4 * (-x1191 * x52 + x1194)
    x1200 = 0.5 * x1199
    x1201 = x371 + x552 * x979
    x1202 = x4 * (-x1201 * x52 + x373 * x46 + x557 * x979)
    x1203 = 0.5 * x1202
    x1204 = x4 * (-x1194 * x52 + x1201)
    x1205 = 0.5 * x1204
    x1206 = x566 * x979 + x887
    x1207 = x1127 * x151 + x889
    x1208 = x1207 * x52
    x1209 = x563 * x979 + x893
    x1210 = x4 * (x1206 - x1208)
    x1211 = 0.5 * x1210
    x1212 = x1209 * x238 + x1211
    x1213 = x1207 * x46
    x1214 = x382 + x563 * x987
    x1215 = x1213 + x1214 * x3
    x1216 = x1206 * x46
    x1217 = x1128 * x151 + x389
    x1218 = x1216 + x1217 * x3
    x1219 = x1209 * x46
    x1220 = x393 + x576 * x987
    x1221 = x4 * (-x1214 * x52 + x1217)
    x1222 = 0.5 * x1221
    x1223 = x418 + x581 * x979
    x1224 = x4 * (-x1223 * x52 + x421 + x587 * x979)
    x1225 = 0.5 * x1224
    x1226 = x4 * (-x1217 * x52 + x1223)
    x1227 = 0.5 * x1226
    x1228 = 3.0 * x462 + x588 * x979
    x1229 = 3.0 * x459 + x589 * x979
    x1230 = x1229 * x52
    x1231 = x592 * x979
    x1232 = x1231 + 3.0 * x436
    x1233 = x1232 * x238
    x1234 = x4 * (x1228 - x1230)
    x1235 = 0.5 * x1234
    x1236 = x1233 + x1235
    x1237 = x1229 * x46
    x1238 = x1231 * x5 + 3.0 * x437
    x1239 = x1238 * x3
    x1240 = x1237 + x1239
    x1241 = x1228 * x46
    x1242 = x444 * x92 + x608 * x979
    x1243 = x1242 * x3
    x1244 = x1241 + x1243
    x1245 = x1232 * x46
    x1246 = x451 * x92 + x612 * x987
    x1247 = x1246 * x3
    x1248 = x1238 * x52
    x1249 = x4 * (x1242 - x1248)
    x1250 = 0.5 * x1249
    x1251 = x1245 * x3
    x1252 = 3.0 * x481 + x619 * x979
    x1253 = x4 * (-x1252 * x52 + x483 * x92 + x625 * x979)
    x1254 = 0.5 * x1253
    x1255 = x4 * (-x1242 * x52 + x1252)
    x1256 = 0.5 * x1255
    x1257 = x4 * x979 * (x628 - x630)
    x1258 = 0.5 * x1257
    x1259 = x4 * x979 * (-x52 * x636 + x640)
    x1260 = 0.5 * x1259
    x1261 = x4 * x979 * (-x52 * x628 + x636)
    x1262 = 0.5 * x1261
    x1263 = x508 + x643 * x979
    x1264 = x505 + x644 * x979
    x1265 = x1264 * x52
    x1266 = x512 + x646 * x979
    x1267 = x4 * (x1263 - x1265)
    x1268 = 0.5 * x1267
    x1269 = x527 + x648 * x979
    x1270 = x4 * (-x1269 * x52 + x525 + x651 * x979)
    x1271 = 0.5 * x1270
    x1272 = x4 * (-x1263 * x52 + x1269)
    x1273 = 0.5 * x1272
    x1274 = 2.0 * x544 + x654 * x979
    x1275 = 2.0 * x542 + x655 * x979
    x1276 = x1275 * x52
    x1277 = 2.0 * x547 + x657 * x979
    x1278 = x4 * (x1274 - x1276)
    x1279 = 0.5 * x1278
    x1280 = 2.0 * x558 + x659 * x979
    x1281 = x4 * (-x1280 * x52 + 2.0 * x556 + x662 * x979)
    x1282 = 0.5 * x1281
    x1283 = x4 * (-x1274 * x52 + x1280)
    x1284 = 0.5 * x1283
    x1285 = 3.0 * x573 + x669 * x979
    x1286 = x1231 * x151 + 3.0 * x571
    x1287 = x1286 * x52
    x1288 = x563 * x92 + x668 * x979
    x1289 = x4 * (x1285 - x1287)
    x1290 = 0.5 * x1289
    x1291 = 3.0 * x585 + x672 * x979
    x1292 = x4 * (-x1291 * x52 + x584 * x92 + x972 * x979)
    x1293 = 0.5 * x1292
    x1294 = x4 * (-x1285 * x52 + x1291)
    x1295 = 0.5 * x1294
    x1296 = 4.0 * x607 + x675 * x979
    x1297 = 4.0 * x604 + x676 * x979
    x1298 = x1297 * x52
    x1299 = x1298 * x3
    x1300 = 4.0 * x611 + x679 * x979
    x1301 = x1300 * x238
    x1302 = x4 * (x1296 - x1298)
    x1303 = 0.5 * x1302
    x1304 = 4.0 * x626 + x683 * x979
    x1305 = x4 * (-x1304 * x52 + 4.0 * x624 + x687 * x979)
    x1306 = 0.5 * x1305
    x1307 = x4 * (-x1296 * x52 + x1304)
    x1308 = 0.5 * x1307
    x1309 = x690**2
    x1310 = x127 + x1309 * x59
    x1311 = x1309 * x63 + x134
    x1312 = x1311 * x46
    x1313 = x1309 * x44
    x1314 = x120 + x1313
    x1315 = x1309 * x48 + x138
    x1316 = x1315 * x46
    x1317 = x1309 * x72
    x1318 = x1317 + x78
    x1319 = x1309 * x29
    x1320 = x1319 + x85
    x1321 = x1320 * x46
    x1322 = x1310 - x1314 * x52
    x1323 = x1309 * x36
    x1324 = x103 + x1323
    x1325 = x1324 * x46
    x1326 = x1320 * x3 + 3.0 * x1325
    x1327 = x1311 - x1315 * x52
    x1328 = x1309 * x5
    x1329 = x46 * (x1328 * x40 + x146)
    x1330 = 21.0675507148243 * x148
    x1331 = x46 * x695
    x1332 = x1331 + x218 + x690 * x725
    x1333 = x225 + x46 * x712 + x690 * x727
    x1334 = x49 * x690
    x1335 = x1334 + x211 + x690 * x720
    x1336 = x231 + x690 * x722 + x690 * x95
    x1337 = x690 * x74
    x1338 = x1337 + x183 + x690 * x730
    x1339 = x690 * x80
    x1340 = x1339 + x190 + x690 * x732
    x1341 = x1332 - x1335 * x52
    x1342 = x191 * x707
    x1343 = x1342 + x206 + x690 * x736
    x1344 = x104 * x1343
    x1345 = x1340 * x3 + x1344
    x1346 = x1333 - x1336 * x52
    x1347 = x46 * (x237 + x690 * x741 + x715)
    x1348 = 55.7394999246661 * x148
    x1349 = x1309 * x254 + x299
    x1350 = x1309 * x259 + x307
    x1351 = x1350 * x46
    x1352 = x1309 * x245 + x293
    x1353 = x1309 * x249 + x312
    x1354 = x1353 * x46
    x1355 = x1309 * x264 + x269
    x1356 = x1309 * x244 + x275
    x1357 = x1356 * x46
    x1358 = x1349 - x1352 * x52
    x1359 = x1328 * x277 + x288
    x1360 = x1359 * x46
    x1361 = x1356 * x3 + 2.0 * x1360
    x1362 = x1350 - x1353 * x52
    x1363 = x46 * (x1309 * x314 + x317)
    x1364 = x104 * x727 + x369 + x690 * x782
    x1365 = x374 + x690 * x784 + x751
    x1366 = x104 * x722 + x365 + x690 * x776
    x1367 = x376 + x690 * x779 + x747
    x1368 = x104 * x732 + x344 + x690 * x788
    x1369 = x348 + x690 * x790 + x737
    x1370 = x1364 - x1366 * x52
    x1371 = x359 + x690 * x794 + 2.0 * x742
    x1372 = x1371 * x46
    x1373 = x1369 * x3 + x1372
    x1374 = x1365 - x1367 * x52
    x1375 = 71.9593849780538 * x148
    x1376 = x416 + x46 * x759 + x690 * x815
    x1377 = x313 * x707 + x422 + x690 * x816
    x1378 = x104 * x1377
    x1379 = x250 * x690 + x412 + x690 * x811
    x1380 = x276 * x707 + x424 + x690 * x812
    x1381 = x104 * x1380
    x1382 = x266 * x690 + x394 + x690 * x819
    x1383 = x400 + x690 * x821 + x707 * x820
    x1384 = x104 * x1383
    x1385 = x1376 - x1379 * x52
    x1386 = x407 + x690 * x824 + x764
    x1387 = x1386 * x46
    x1388 = x1383 * x3 + x1387
    x1389 = x1377 - x1380 * x52
    x1390 = 124.637310863398 * x148
    x1391 = x1309 * x442 + x479
    x1392 = x1309 * x444 + x485
    x1393 = x1392 * x46
    x1394 = x1309 * x434 + x474
    x1395 = x1309 * x430
    x1396 = x1395 * x5 + x487
    x1397 = x1396 * x46
    x1398 = x1309 * x449 + x453
    x1399 = x1309 * x451 + x456
    x1400 = x1399 * x46
    x1401 = x1391 - x1394 * x52
    x1402 = x1395 + x467
    x1403 = x1402 * x46
    x1404 = x1399 * x3 + x1403
    x1405 = x1392 - x1396 * x52
    x1406 = x503 + x690 * x861 + 3.0 * x804
    x1407 = x1406 * x46
    x1408 = x521 + x690 * x877 + x784 * x92
    x1409 = x500 + x690 * x864 + 3.0 * x801
    x1410 = x1409 * x46
    x1411 = x517 + x690 * x873 + x779 * x92
    x1412 = x496 + x690 * x868 + 3.0 * x795
    x1413 = x1406 - x1409 * x52
    x1414 = x1412 * x46
    x1415 = x515 + x690 * x882 + x790 * x92
    x1416 = x1408 - x1411 * x52
    x1417 = x540 + x690 * x888 + 2.0 * x834
    x1418 = x1417 * x46
    x1419 = x553 + x690 * x901 + x817
    x1420 = x537 + x690 * x890 + 2.0 * x830
    x1421 = x1420 * x46
    x1422 = x550 + x690 * x898 + x813
    x1423 = x533 + x690 * x894 + 2.0 * x825
    x1424 = x1417 - x1420 * x52
    x1425 = x1423 * x46
    x1426 = x549 + x690 * x904 + x822
    x1427 = x1419 - x1422 * x52
    x1428 = x569 + x690 * x909 + x854
    x1429 = x1428 * x46
    x1430 = x46 * x843 + x582 + x690 * x919
    x1431 = x567 + x690 * x910 + x852
    x1432 = x1431 * x46
    x1433 = x579 + x690 * x916 + x841
    x1434 = x562 + x690 * x912 + x847
    x1435 = x1428 - x1431 * x52
    x1436 = x1434 * x46
    x1437 = x46 * x845 + x577 + x690 * x922
    x1438 = x1430 - x1433 * x52
    x1439 = x1309 * x588 + x602
    x1440 = x1439 * x46
    x1441 = x1309 * x608 + x620
    x1442 = x1309 * x589 + x599
    x1443 = x1442 * x46
    x1444 = x1309 * x592
    x1445 = x1444 * x5 + x616
    x1446 = x1444 + x594
    x1447 = x1439 - x1442 * x52
    x1448 = x1446 * x46
    x1449 = x1328 * x612 + x614
    x1450 = x1441 - x1445 * x52
    x1451 = x641 + x690 * x939 + 4.0 * x875
    x1452 = x637 + x690 * x942 + 4.0 * x871
    x1453 = x1452 * x52
    x1454 = x634 + x690 * x946 + 4.0 * x879
    x1455 = x1451 - x1453
    x1456 = x652 + x690 * x951 + 3.0 * x900
    x1457 = x649 + x690 * x952 + 3.0 * x897
    x1458 = x1457 * x52
    x1459 = x647 + x690 * x954 + 3.0 * x903
    x1460 = x1456 - x1458
    x1461 = x663 + x690 * x959 + 2.0 * x918
    x1462 = x660 + x690 * x960 + 2.0 * x915
    x1463 = x1462 * x52
    x1464 = x658 + x690 * x962 + 2.0 * x921
    x1465 = x1461 - x1463
    x1466 = x673 + x690 * x967 + x931
    x1467 = x670 + x690 * x968 + x929
    x1468 = x1467 * x52
    x1469 = x667 + x690 * x970 + x933
    x1470 = x1466 - x1468
    x1471 = x1309 * x675 + x688
    x1472 = x1309 * x676 + x684
    x1473 = x1472 * x52
    x1474 = x1309 * x679 + x681
    x1475 = x1471 - x1473
    x1476 = x4 * x979 * (x59 * x690 - x699)
    x1477 = x979 * (x701 + x702)
    x1478 = x4 * x979 * (x63 * x690 - x704)
    x1479 = x64 * x979
    x1480 = x1479 + x724 * x979
    x1481 = x979 * x99
    x1482 = x1481 + x726 * x979
    x1483 = x49 * x979
    x1484 = x1483 + x719 * x979
    x1485 = x95 * x979
    x1486 = x1485 + x721 * x979
    x1487 = x74 * x979
    x1488 = x1487 + x729 * x979
    x1489 = x80 * x979
    x1490 = x1489 + x731 * x979
    x1491 = x4 * (x1480 - x1484 * x52)
    x1492 = x191 * x987
    x1493 = x1492 + x735 * x979
    x1494 = x104 * x1493
    x1495 = x1490 * x3 + x1494
    x1496 = x4 * (x1482 - x1486 * x52)
    x1497 = x46 * (x740 * x979 + x994)
    x1498 = 96.5436458580033 * x148
    x1499 = x1032 * x690
    x1500 = x4 * x690 * (x1030 - x1040)
    x1501 = x690 * (x1043 + x1046)
    x1502 = x4 * x690 * (x1032 - x1048)
    x1503 = x1052 * x690
    x1504 = x1001 * x104 + x781 * x979
    x1505 = x1016 + x783 * x979
    x1506 = x979 * (x774 + x775)
    x1507 = x1014 + x778 * x979
    x1508 = x979 * (x786 + x787)
    x1509 = x1007 + x789 * x979
    x1510 = x4 * (x1504 - x1506 * x52)
    x1511 = x979 * (x777 + x793)
    x1512 = x1511 * x46
    x1513 = x1509 * x3 + x1512
    x1514 = x4 * (x1505 - x1507 * x52)
    x1515 = x1033 + x1097 * x690
    x1516 = x1061 + x1098 * x690
    x1517 = x104 * x1516
    x1518 = x1027 + x1093 * x690
    x1519 = x1057 + x1094 * x690
    x1520 = x104 * x1519
    x1521 = x1038 + x1101 * x690
    x1522 = x1045 + x1102 * x690
    x1523 = x104 * x1522
    x1524 = x4 * (x1515 - x1518 * x52)
    x1525 = x1052 + x1106 * x690
    x1526 = x1525 * x46
    x1527 = x1522 * x3 + x1526
    x1528 = x4 * (x1516 - x1519 * x52)
    x1529 = 215.878154934161 * x148
    x1530 = x1135 * x690
    x1531 = x4 * x690 * (x1133 - x1143)
    x1532 = x1147 * x690
    x1533 = x1148 * x690 + x1532
    x1534 = x4 * x690 * (x1135 - x1150)
    x1535 = x979 * (x859 + x860)
    x1536 = x1535 * x46
    x1537 = x1074 * x92 + x876 * x979
    x1538 = x979 * (x862 + x863)
    x1539 = x1538 * x46
    x1540 = x867 * x987 + x872 * x979
    x1541 = x979 * (x866 + x867)
    x1542 = x4 * (x1535 - x1538 * x52)
    x1543 = x1541 * x46
    x1544 = x979 * (x880 + x881)
    x1545 = x4 * (x1537 - x1540 * x52)
    x1546 = 2.0 * x1116
    x1547 = x1182 * x690 + x1546
    x1548 = x1547 * x46
    x1549 = x1099 + x1194 * x690
    x1550 = 2.0 * x1113
    x1551 = x1183 * x690 + x1550
    x1552 = x1551 * x46
    x1553 = x1095 + x1191 * x690
    x1554 = 2.0 * x1107
    x1555 = x1186 * x690 + x1554
    x1556 = x4 * (x1547 - x1551 * x52)
    x1557 = x1555 * x46
    x1558 = x1103 + x1198 * x690
    x1559 = x4 * (x1549 - x1553 * x52)
    x1560 = x1159 + x1206 * x690
    x1561 = x1560 * x46
    x1562 = x1136 + x1217 * x690
    x1563 = x1156 + x1207 * x690
    x1564 = x1563 * x46
    x1565 = x1130 + x1214 * x690
    x1566 = x1147 + x1209 * x690
    x1567 = x4 * (x1560 - x1563 * x52)
    x1568 = x1566 * x46
    x1569 = x1141 + x1220 * x690
    x1570 = x4 * (x1562 - x1565 * x52)
    x1571 = x1241 * x690
    x1572 = x1237 * x690
    x1573 = x4 * x690 * (x1228 - x1230)
    x1574 = x1245 * x690
    x1575 = x4 * x690 * (x1242 - x1248)
    x1576 = x979 * (x937 + x938)
    x1577 = x979 * (x940 + x941)
    x1578 = x1577 * x52
    x1579 = x979 * (x944 + x945)
    x1580 = x4 * (x1576 - x1578)
    x1581 = 3.0 * x1193 + x1263 * x690
    x1582 = 3.0 * x1190 + x1264 * x690
    x1583 = x1582 * x52
    x1584 = 3.0 * x1196 + x1266 * x690
    x1585 = x4 * (x1581 - x1583)
    x1586 = 2.0 * x1216 + x1274 * x690
    x1587 = 2.0 * x1213 + x1275 * x690
    x1588 = x1587 * x52
    x1589 = 2.0 * x1219 + x1277 * x690
    x1590 = x4 * (x1586 - x1588)
    x1591 = x1241 + x1285 * x690
    x1592 = x1237 + x1286 * x690
    x1593 = x1592 * x52
    x1594 = x1245 + x1288 * x690
    x1595 = x4 * (x1591 - x1593)
    x1596 = x4 * x690 * (x1296 - x1298)
    x1597 = x979**2
    x1598 = x127 + x1597 * x59
    x1599 = x134 + x1597 * x63
    x1600 = x1599 * x46
    x1601 = x120 + x1597 * x44
    x1602 = x138 + x1597 * x48
    x1603 = x1602 * x46
    x1604 = x1597 * x72 + x78
    x1605 = x1604 * x3
    x1606 = x1597 * x29 + x85
    x1607 = x1606 * x46
    x1608 = 4.0 * x1607
    x1609 = x1601 * x52
    x1610 = x1598 - x1609
    x1611 = x1610 * x77
    x1612 = x1606 * x3
    x1613 = x103 + x1597 * x36
    x1614 = x1613 * x46
    x1615 = 3.0 * x1614
    x1616 = x1612 + x1615
    x1617 = x1599 - x1602 * x52
    x1618 = x1617 * x77
    x1619 = x1597 * x5
    x1620 = x46 * (x146 + x1619 * x40)
    x1621 = x1597 * x166 + x218
    x1622 = x1597 * x172 + x225
    x1623 = x1622 * x46
    x1624 = x156 * x1597 + x211
    x1625 = x1597 * x161 + x231
    x1626 = x1625 * x46
    x1627 = x1597 * x178 + x183
    x1628 = x155 * x1597 + x190
    x1629 = x1628 * x46
    x1630 = x1621 - x1624 * x52
    x1631 = x1630 * x77
    x1632 = x1619 * x193 + x206
    x1633 = x1632 * x46
    x1634 = 2.0 * x1633
    x1635 = x1628 * x3 + x1634
    x1636 = x1622 - x1625 * x52
    x1637 = x1636 * x77
    x1638 = x46 * (x1597 * x233 + x237)
    x1639 = x1030 * x979 + x1479 + x299
    x1640 = x1032 * x979 + x1481 + x307
    x1641 = x1640 * x46
    x1642 = x1024 * x979 + x1483 + x293
    x1643 = x1026 * x979 + x1485 + x312
    x1644 = x1643 * x46
    x1645 = x1035 * x979 + x1487 + x269
    x1646 = x1645 * x3
    x1647 = x1037 * x979 + x1489 + x275
    x1648 = x1647 * x46
    x1649 = 3.0 * x1648
    x1650 = x1642 * x52
    x1651 = x1639 - x1650
    x1652 = x1651 * x77
    x1653 = x1647 * x3
    x1654 = x1044 * x979 + x1492 + x288
    x1655 = x1654 * x46
    x1656 = 2.0 * x1655
    x1657 = x1653 + x1656
    x1658 = x1640 - x1643 * x52
    x1659 = x1658 * x77
    x1660 = x46 * (x1051 * x979 + x317 + x994)
    x1661 = x1597 * x332 + x369
    x1662 = x1597 * x334 + x374
    x1663 = x1662 * x46
    x1664 = x1597 * x325 + x365
    x1665 = x1597 * x322
    x1666 = x1665 * x5 + x376
    x1667 = x1666 * x46
    x1668 = x1597 * x339 + x344
    x1669 = x1597 * x341 + x348
    x1670 = x1669 * x46
    x1671 = x1661 - x1664 * x52
    x1672 = x1671 * x77
    x1673 = x1665 + x359
    x1674 = x1673 * x46
    x1675 = x1669 * x3 + x1674
    x1676 = x1662 - x1666 * x52
    x1677 = x1676 * x77
    x1678 = x1001 * x46 + x1097 * x979 + x416
    x1679 = x1098 * x979 + x232 * x987 + x422
    x1680 = x104 * x1679
    x1681 = x1093 * x979 + x162 * x979 + x412
    x1682 = x1094 * x979 + x192 * x987 + x424
    x1683 = x104 * x1682
    x1684 = x1101 * x979 + x180 * x979 + x394
    x1685 = x1006 * x106 + x1102 * x979 + x400
    x1686 = x104 * x1685
    x1687 = x1678 - x1681 * x52
    x1688 = x1687 * x77
    x1689 = x1011 + x1106 * x979 + x407
    x1690 = x1689 * x46
    x1691 = x1685 * x3 + x1690
    x1692 = x1679 - x1682 * x52
    x1693 = x1692 * x77
    x1694 = 2.0 * x1033 + x1133 * x979 + x479
    x1695 = x1062 + x1135 * x979 + x485
    x1696 = x1695 * x46
    x1697 = 2.0 * x1027 + x1125 * x979 + x474
    x1698 = x1058 + x1129 * x979 + x487
    x1699 = x1698 * x46
    x1700 = 2.0 * x1038 + x1138 * x979 + x453
    x1701 = x1700 * x3
    x1702 = x1046 + x1140 * x979 + x456
    x1703 = x1702 * x46
    x1704 = 2.0 * x1703
    x1705 = x1697 * x52
    x1706 = x1694 - x1705
    x1707 = x1706 * x77
    x1708 = 2.0 * x1052 + x1146 * x979 + x467
    x1709 = x1708 * x46
    x1710 = x1702 * x3
    x1711 = x1709 + x1710
    x1712 = x1695 - x1698 * x52
    x1713 = x1712 * x77
    x1714 = x1597 * x490 + x503
    x1715 = x1714 * x46
    x1716 = x1597 * x509 + x521
    x1717 = x1597 * x491 + x500
    x1718 = x1717 * x46
    x1719 = x1597 * x494
    x1720 = x1719 * x5 + x517
    x1721 = x1719 + x496
    x1722 = x1714 - x1717 * x52
    x1723 = x1722 * x77
    x1724 = x1721 * x46
    x1725 = x1619 * x513 + x515
    x1726 = x1716 - x1720 * x52
    x1727 = x1726 * x77
    x1728 = x1085 + x1182 * x979 + x540
    x1729 = x1728 * x46
    x1730 = x1074 * x46 + x1194 * x979 + x553
    x1731 = x1083 + x1183 * x979 + x537
    x1732 = x1731 * x46
    x1733 = x1072 + x1191 * x979 + x550
    x1734 = x1071 + x1186 * x979 + x533
    x1735 = x1728 - x1731 * x52
    x1736 = x1735 * x77
    x1737 = x1734 * x46
    x1738 = x1197 * x979 + x1198 * x979 + x549
    x1739 = x1730 - x1733 * x52
    x1740 = x1739 * x77
    x1741 = x1206 * x979 + x1546 + x569
    x1742 = x1741 * x46
    x1743 = x1099 + x1217 * x979 + x582
    x1744 = x1207 * x979 + x1550 + x567
    x1745 = x1744 * x46
    x1746 = x1095 + x1214 * x979 + x579
    x1747 = x1209 * x979 + x1554 + x562
    x1748 = x1741 - x1744 * x52
    x1749 = x1748 * x77
    x1750 = x1747 * x46
    x1751 = x1103 + x1220 * x979 + x577
    x1752 = x1743 - x1746 * x52
    x1753 = x1752 * x77
    x1754 = 3.0 * x1159 + x1228 * x979 + x602
    x1755 = x1754 * x46
    x1756 = 3.0 * x1136 + x1242 * x979 + x620
    x1757 = 3.0 * x1156 + x1229 * x979 + x599
    x1758 = x1757 * x46
    x1759 = 3.0 * x1130 + x1238 * x979 + x616
    x1760 = 3.0 * x1147 + x1232 * x979 + x594
    x1761 = x1754 - x1757 * x52
    x1762 = x1761 * x77
    x1763 = x1760 * x46
    x1764 = 3.0 * x1141 + x1246 * x979 + x614
    x1765 = x1764 * x3
    x1766 = x1759 * x52
    x1767 = x1756 - x1766
    x1768 = x1767 * x77
    x1769 = x1763 * x3
    x1770 = x1597 * x628 + x641
    x1771 = x1597 * x629 + x637
    x1772 = x1771 * x52
    x1773 = x1597 * x632 + x634
    x1774 = x1770 - x1772
    x1775 = x1774 * x77
    x1776 = x1173 + x1263 * x979 + x652
    x1777 = x1171 + x1264 * x979 + x649
    x1778 = x1777 * x52
    x1779 = x1175 + x1266 * x979 + x647
    x1780 = x1776 - x1778
    x1781 = x1780 * x77
    x1782 = 2.0 * x1193 + x1274 * x979 + x663
    x1783 = 2.0 * x1190 + x1275 * x979 + x660
    x1784 = x1783 * x52
    x1785 = 2.0 * x1196 + x1277 * x979 + x658
    x1786 = x1782 - x1784
    x1787 = x1786 * x77
    x1788 = 3.0 * x1216 + x1285 * x979 + x673
    x1789 = 3.0 * x1213 + x1286 * x979 + x670
    x1790 = x1789 * x52
    x1791 = 3.0 * x1219 + x1288 * x979 + x667
    x1792 = x1788 - x1790
    x1793 = x1792 * x77
    x1794 = 4.0 * x1241 + x1296 * x979 + x688
    x1795 = 4.0 * x1237 + x1297 * x979 + x684
    x1796 = x1795 * x52
    x1797 = 4.0 * x1245 + x1300 * x979 + x681
    x1798 = x1797 * x238
    x1799 = x1794 - x1796
    x1800 = x1799 * x77
    x1801 = x1318 * x690 + x700
    x1802 = x1320 * x690 + x705
    x1803 = x1802 * x46
    x1804 = x1310 * x690 - x52 * (x1314 * x690 + x717) + x716
    x1805 = x46 * (x1324 * x690 + x714)
    x1806 = 21.0675507148243 * x148
    x1807 = x1321 + x1338 * x690 + x733
    x1808 = x1325 + x1340 * x690 + x739
    x1809 = x1312 + x1332 * x690 - x52 * (x1316 + x1335 * x690 + x756) + x755
    x1810 = x104 * (x1329 + x1343 * x690 + x753)
    x1811 = 55.7394999246661 * x148
    x1812 = x1355 * x690 + x761
    x1813 = x1356 * x690 + x763
    x1814 = x1813 * x46
    x1815 = x1349 * x690 - x52 * (x1352 * x690 + x773) + x772
    x1816 = x46 * (x1359 * x690 + x771)
    x1817 = x104 * x1340 + x1368 * x690 + x791
    x1818 = x1344 + x1369 * x690 + x797
    x1819 = (
        x104 * x1333 + x1364 * x690 - x52 * (x104 * x1336 + x1366 * x690 + x809) + x808
    )
    x1820 = x46 * (2.0 * x1347 + x1371 * x690 + x806)
    x1821 = 71.9593849780538 * x148
    x1822 = x1357 + x1382 * x690 + x823
    x1823 = x1360 + x1383 * x690 + x827
    x1824 = x104 * x1823
    x1825 = x1351 + x1376 * x690 - x52 * (x1354 + x1379 * x690 + x839) + x838
    x1826 = x46 * (x1363 + x1386 * x690 + x836)
    x1827 = 124.637310863398 * x148
    x1828 = x1398 * x690 + x846
    x1829 = x1399 * x690 + x850
    x1830 = x1829 * x46
    x1831 = x1391 * x690 - x52 * (x1394 * x690 + x858) + x857
    x1832 = x46 * (x1402 * x690 + x856)
    x1833 = 3.0 * x1372 + x1412 * x690 + x869
    x1834 = x1833 * x46
    x1835 = x1369 * x92 + x1415 * x690 + x883
    x1836 = x1365 * x92 + x1408 * x690 - x52 * (x1367 * x92 + x1411 * x690 + x886) + x885
    x1837 = 2.0 * x1387 + x1423 * x690 + x895
    x1838 = x1837 * x46
    x1839 = x1384 + x1426 * x690 + x905
    x1840 = x1378 + x1419 * x690 - x52 * (x1381 + x1422 * x690 + x908) + x907
    x1841 = x1403 + x1434 * x690 + x913
    x1842 = x1841 * x46
    x1843 = x1400 + x1437 * x690 + x923
    x1844 = x1393 + x1430 * x690 - x52 * (x1397 + x1433 * x690 + x926) + x925
    x1845 = x1446 * x690 + x927
    x1846 = x1845 * x46
    x1847 = x1449 * x690 + x934
    x1848 = x1441 * x690 - x52 * (x1445 * x690 + x936) + x935
    x1849 = 4.0 * x1414 + x1454 * x690 + x947
    x1850 = 4.0 * x1407 + x1451 * x690 - x52 * (4.0 * x1410 + x1452 * x690 + x950) + x949
    x1851 = 3.0 * x1425 + x1459 * x690 + x955
    x1852 = 3.0 * x1418 + x1456 * x690 - x52 * (3.0 * x1421 + x1457 * x690 + x958) + x957
    x1853 = 2.0 * x1436 + x1464 * x690 + x963
    x1854 = 2.0 * x1429 + x1461 * x690 - x52 * (2.0 * x1432 + x1462 * x690 + x966) + x965
    x1855 = x1448 + x1469 * x690 + x971
    x1856 = x1440 + x1466 * x690 - x52 * (x1443 + x1467 * x690 + x975) + x974
    x1857 = x1474 * x690 + x976
    x1858 = x1471 * x690 - x52 * (x1472 * x690 + x978) + x977
    x1859 = x1317 * x979 + x983
    x1860 = x1319 * x979 + x986
    x1861 = x1860 * x46
    x1862 = x1309 * x997 - x52 * (x1313 * x979 + x999) + x996
    x1863 = x46 * (x1323 * x979 + x993)
    x1864 = x1004 + x1337 * x979 + x1488 * x690
    x1865 = x1010 + x1339 * x979 + x1490 * x690
    x1866 = (
        x1021 + x1331 * x979 + x1480 * x690 - x52 * (x1023 + x1334 * x979 + x1484 * x690)
    )
    x1867 = x104 * (x1019 + x1342 * x979 + x1493 * x690)
    x1868 = x1035 * x1309 + x1042
    x1869 = x1037 * x1309 + x1050
    x1870 = x1869 * x46
    x1871 = x1030 * x1309 + x1068 - x52 * (x1024 * x1309 + x1070)
    x1872 = x46 * (x1044 * x1309 + x1065)
    x1873 = x104 * x1490 + x1077 + x1508 * x690
    x1874 = x1081 + x1494 + x1509 * x690
    x1875 = (
        x104 * x1482 + x1090 + x1504 * x690 - x52 * (x104 * x1486 + x1092 + x1506 * x690)
    )
    x1876 = x46 * (x1088 + 2.0 * x1497 + x1511 * x690)
    x1877 = 160.906076430005 * x148
    x1878 = x1038 * x690 + x1105 + x1521 * x690
    x1879 = x1045 * x690 + x1110 + x1522 * x690
    x1880 = x104 * x1879
    x1881 = (
        x1122 + x1499 * x46 + x1515 * x690 - x52 * (x1027 * x690 + x1124 + x1518 * x690)
    )
    x1882 = x46 * (x1119 + x1503 + x1525 * x690)
    x1883 = 278.69749962333 * x148
    x1884 = x1138 * x1309 + x1145
    x1885 = x1140 * x1309 + x1152
    x1886 = x1885 * x46
    x1887 = x1133 * x1309 + x1165 - x52 * (x1125 * x1309 + x1167)
    x1888 = x46 * (x1146 * x1309 + x1162)
    x1889 = x1169 + 3.0 * x1512 + x1541 * x690
    x1890 = x1889 * x46
    x1891 = x1177 + x1509 * x92 + x1544 * x690
    x1892 = (
        x1179 + x1505 * x92 + x1537 * x690 - x52 * (x1181 + x1507 * x92 + x1540 * x690)
    )
    x1893 = x1188 + 2.0 * x1526 + x1555 * x690
    x1894 = x1893 * x46
    x1895 = x1200 + x1523 + x1558 * x690
    x1896 = x1203 + x1517 + x1549 * x690 - x52 * (x1205 + x1520 + x1553 * x690)
    x1897 = x1211 + x1532 + x1566 * x690
    x1898 = x1897 * x46
    x1899 = x1141 * x690 + x1222 + x1569 * x690
    x1900 = (
        x1225 + x1530 * x46 + x1562 * x690 - x52 * (x1130 * x690 + x1227 + x1565 * x690)
    )
    x1901 = x1232 * x1309 + x1235
    x1902 = x1901 * x46
    x1903 = x1246 * x1309 + x1250
    x1904 = x1242 * x1309 + x1254 - x52 * (x1238 * x1309 + x1256)
    x1905 = x1258 + 4.0 * x1543 + x1579 * x690
    x1906 = (
        x1260 + 4.0 * x1536 + x1576 * x690 - x52 * (x1262 + 4.0 * x1539 + x1577 * x690)
    )
    x1907 = x1268 + 3.0 * x1557 + x1584 * x690
    x1908 = (
        x1271 + 3.0 * x1548 + x1581 * x690 - x52 * (x1273 + 3.0 * x1552 + x1582 * x690)
    )
    x1909 = x1279 + 2.0 * x1568 + x1589 * x690
    x1910 = (
        x1282 + 2.0 * x1561 + x1586 * x690 - x52 * (x1284 + 2.0 * x1564 + x1587 * x690)
    )
    x1911 = x1290 + x1574 + x1594 * x690
    x1912 = x1293 + x1571 + x1591 * x690 - x52 * (x1295 + x1572 + x1592 * x690)
    x1913 = x1300 * x1309 + x1303
    x1914 = x1296 * x1309 + x1306 - x52 * (x1297 * x1309 + x1308)
    x1915 = x4 * x690 * (x1598 - x1609)
    x1916 = x1607 + x1627 * x690
    x1917 = x1614 + x1628 * x690
    x1918 = x4 * (x1600 + x1621 * x690 - x52 * (x1603 + x1624 * x690))
    x1919 = x104 * (x1620 + x1632 * x690)
    x1920 = x4 * x690 * (x1639 - x1650)
    x1921 = 2.0 * x1629 + x1668 * x690
    x1922 = x1634 + x1669 * x690
    x1923 = x4 * (2.0 * x1623 + x1661 * x690 - x52 * (2.0 * x1626 + x1664 * x690))
    x1924 = x46 * (2.0 * x1638 + x1673 * x690)
    x1925 = x1648 + x1684 * x690
    x1926 = x1655 + x1685 * x690
    x1927 = x104 * x1926
    x1928 = x4 * (x1641 + x1678 * x690 - x52 * (x1644 + x1681 * x690))
    x1929 = x46 * (x1660 + x1689 * x690)
    x1930 = x4 * x690 * (x1694 - x1705)
    x1931 = x1709 * x690
    x1932 = 3.0 * x1674 + x1721 * x690
    x1933 = x1932 * x46
    x1934 = 3.0 * x1670 + x1725 * x690
    x1935 = x4 * (3.0 * x1663 + x1716 * x690 - x52 * (3.0 * x1667 + x1720 * x690))
    x1936 = 2.0 * x1690
    x1937 = x1734 * x690 + x1936
    x1938 = x1937 * x46
    x1939 = x1686 + x1738 * x690
    x1940 = x4 * (x1680 + x1730 * x690 - x52 * (x1683 + x1733 * x690))
    x1941 = x1709 + x1747 * x690
    x1942 = x1941 * x46
    x1943 = x1703 + x1751 * x690
    x1944 = x4 * (x1696 + x1743 * x690 - x52 * (x1699 + x1746 * x690))
    x1945 = x1763 * x690
    x1946 = x4 * x690 * (x1756 - x1766)
    x1947 = 4.0 * x1724 + x1773 * x690
    x1948 = x4 * (4.0 * x1715 + x1770 * x690 - x52 * (4.0 * x1718 + x1771 * x690))
    x1949 = 3.0 * x1737 + x1779 * x690
    x1950 = x4 * (3.0 * x1729 + x1776 * x690 - x52 * (3.0 * x1732 + x1777 * x690))
    x1951 = 2.0 * x1750 + x1785 * x690
    x1952 = x4 * (2.0 * x1742 + x1782 * x690 - x52 * (2.0 * x1745 + x1783 * x690))
    x1953 = x1763 + x1791 * x690
    x1954 = x4 * (x1755 + x1788 * x690 - x52 * (x1758 + x1789 * x690))
    x1955 = x4 * x690 * (x1794 - x1796)
    x1956 = x1604 * x979 + x982
    x1957 = x1956 * x3
    x1958 = x1606 * x979 + x985
    x1959 = x1958 * x46
    x1960 = 4.0 * x1959
    x1961 = x1598 * x979 - x52 * (x1601 * x979 + x998) + x995
    x1962 = x1961 * x77
    x1963 = x46 * (x1613 * x979 + x992)
    x1964 = x1003 + x1627 * x979
    x1965 = x1009 + x1628 * x979
    x1966 = x1965 * x46
    x1967 = x1020 + x1621 * x979 - x52 * (x1022 + x1624 * x979)
    x1968 = x1967 * x77
    x1969 = x46 * (x1018 + x1632 * x979)
    x1970 = 2.0 * x1969
    x1971 = x1041 + x1607 + x1645 * x979
    x1972 = x1971 * x3
    x1973 = x1049 + x1614 + x1647 * x979
    x1974 = x1973 * x46
    x1975 = 3.0 * x1974
    x1976 = x1067 + x1600 + x1639 * x979 - x52 * (x1069 + x1603 + x1642 * x979)
    x1977 = x1976 * x77
    x1978 = x46 * (x1064 + x1620 + x1654 * x979)
    x1979 = 2.0 * x1978
    x1980 = x1076 + x1668 * x979
    x1981 = x1080 + x1669 * x979
    x1982 = x1981 * x46
    x1983 = x1089 + x1661 * x979 - x52 * (x1091 + x1664 * x979)
    x1984 = x1983 * x77
    x1985 = x46 * (x1087 + x1673 * x979)
    x1986 = x1104 + x1629 + x1684 * x979
    x1987 = x1109 + x1633 + x1685 * x979
    x1988 = x104 * x1987
    x1989 = x1121 + x1623 + x1678 * x979 - x52 * (x1123 + x1626 + x1681 * x979)
    x1990 = x1989 * x77
    x1991 = x46 * (x1118 + x1638 + x1689 * x979)
    x1992 = x1144 + 2.0 * x1648 + x1700 * x979
    x1993 = x1992 * x3
    x1994 = x1151 + x1656 + x1702 * x979
    x1995 = x1994 * x46
    x1996 = 2.0 * x1995
    x1997 = (
        x1164 + 2.0 * x1641 + x1694 * x979 - x52 * (x1166 + 2.0 * x1644 + x1697 * x979)
    )
    x1998 = x1997 * x77
    x1999 = x46 * (x1161 + 2.0 * x1660 + x1708 * x979)
    x2000 = x1168 + x1721 * x979
    x2001 = x2000 * x46
    x2002 = x1176 + x1725 * x979
    x2003 = x1178 + x1716 * x979 - x52 * (x1180 + x1720 * x979)
    x2004 = x2003 * x77
    x2005 = x1187 + x1674 + x1734 * x979
    x2006 = x2005 * x46
    x2007 = x1199 + x1670 + x1738 * x979
    x2008 = x1202 + x1663 + x1730 * x979 - x52 * (x1204 + x1667 + x1733 * x979)
    x2009 = x2008 * x77
    x2010 = x1210 + x1747 * x979 + x1936
    x2011 = x2010 * x46
    x2012 = x1221 + x1686 + x1751 * x979
    x2013 = x1224 + x1680 + x1743 * x979 - x52 * (x1226 + x1683 + x1746 * x979)
    x2014 = x2013 * x77
    x2015 = x1234 + 3.0 * x1709 + x1760 * x979
    x2016 = x2015 * x46
    x2017 = x1249 + 3.0 * x1703 + x1764 * x979
    x2018 = x2017 * x3
    x2019 = (
        x1253 + 3.0 * x1696 + x1756 * x979 - x52 * (x1255 + 3.0 * x1699 + x1759 * x979)
    )
    x2020 = x2019 * x77
    x2021 = x1257 + x1773 * x979
    x2022 = x1259 + x1770 * x979 - x52 * (x1261 + x1771 * x979)
    x2023 = x2022 * x77
    x2024 = x1267 + x1724 + x1779 * x979
    x2025 = x1270 + x1715 + x1776 * x979 - x52 * (x1272 + x1718 + x1777 * x979)
    x2026 = x2025 * x77
    x2027 = x1278 + 2.0 * x1737 + x1785 * x979
    x2028 = (
        x1281 + 2.0 * x1729 + x1782 * x979 - x52 * (x1283 + 2.0 * x1732 + x1783 * x979)
    )
    x2029 = x2028 * x77
    x2030 = x1289 + 3.0 * x1750 + x1791 * x979
    x2031 = (
        x1292 + 3.0 * x1742 + x1788 * x979 - x52 * (x1294 + 3.0 * x1745 + x1789 * x979)
    )
    x2032 = x2031 * x77
    x2033 = x1302 + 4.0 * x1763 + x1797 * x979
    x2034 = (
        x1305 + 4.0 * x1755 + x1794 * x979 - x52 * (x1307 + 4.0 * x1758 + x1795 * x979)
    )
    x2035 = x2034 * x77
    x2036 = x129 * x1322 + x1801 * x690
    x2037 = x46 * (x129 * x1327 + x1802 * x690)
    x2038 = x129 * x1341 + x1803 + x1807 * x690
    x2039 = x129 * x1346 + x1805 + x1808 * x690
    x2040 = x129 * x1358 + x1812 * x690
    x2041 = x46 * (x129 * x1362 + x1813 * x690)
    x2042 = x104 * x1808 + x129 * x1370 + x1817 * x690
    x2043 = x129 * x1374 + x1810 + x1818 * x690
    x2044 = x129 * x1385 + x1814 + x1822 * x690
    x2045 = x104 * (x129 * x1389 + x1816 + x1823 * x690)
    x2046 = x129 * x1401 + x1828 * x690
    x2047 = x46 * (x129 * x1405 + x1829 * x690)
    x2048 = x46 * (x129 * x1413 + 3.0 * x1820 + x1833 * x690)
    x2049 = x129 * x1416 + x1818 * x92 + x1835 * x690
    x2050 = x46 * (x129 * x1424 + 2.0 * x1826 + x1837 * x690)
    x2051 = x129 * x1427 + x1824 + x1839 * x690
    x2052 = x46 * (x129 * x1435 + x1832 + x1841 * x690)
    x2053 = x129 * x1438 + x1830 + x1843 * x690
    x2054 = x46 * (x129 * x1447 + x1845 * x690)
    x2055 = x129 * x1450 + x1847 * x690
    x2056 = x129 * x1455 + 4.0 * x1834 + x1849 * x690
    x2057 = x3 * x718
    x2058 = x129 * x1460 + 3.0 * x1838 + x1851 * x690
    x2059 = x3 * x757
    x2060 = x129 * x1465 + 2.0 * x1842 + x1853 * x690
    x2061 = x3 * x810
    x2062 = x129 * x1470 + x1846 + x1855 * x690
    x2063 = x129 * x1475 + x1857 * x690
    x2064 = x1476 + x1859 * x690
    x2065 = x46 * (x1478 + x1860 * x690)
    x2066 = x1491 + x1861 + x1864 * x690
    x2067 = x1496 + x1863 + x1865 * x690
    x2068 = x1500 + x1868 * x690
    x2069 = x46 * (x1502 + x1869 * x690)
    x2070 = x104 * x1865 + x1510 + x1873 * x690
    x2071 = x1514 + x1867 + x1874 * x690
    x2072 = x1524 + x1870 + x1878 * x690
    x2073 = x104 * (x1528 + x1872 + x1879 * x690)
    x2074 = x1531 + x1884 * x690
    x2075 = x46 * (x1534 + x1885 * x690)
    x2076 = x46 * (x1542 + 3.0 * x1876 + x1889 * x690)
    x2077 = x1545 + x1874 * x92 + x1891 * x690
    x2078 = x46 * (x1556 + 2.0 * x1882 + x1893 * x690)
    x2079 = x1559 + x1880 + x1895 * x690
    x2080 = x46 * (x1567 + x1888 + x1897 * x690)
    x2081 = x1570 + x1886 + x1899 * x690
    x2082 = x46 * (x1573 + x1901 * x690)
    x2083 = x1575 + x1903 * x690
    x2084 = x1580 + 4.0 * x1890 + x1905 * x690
    x2085 = x1585 + 3.0 * x1894 + x1907 * x690
    x2086 = x1498 * x3
    x2087 = x1590 + 2.0 * x1898 + x1909 * x690
    x2088 = x1390 * x3
    x2089 = x1595 + x1902 + x1911 * x690
    x2090 = x1596 + x1913 * x690
    x2091 = x1309 * x1604 + x1611
    x2092 = x46 * (x1309 * x1606 + x1618)
    x2093 = x1607 * x690 + x1631 + x1916 * x690
    x2094 = x1614 * x690 + x1637 + x1917 * x690
    x2095 = x1309 * x1645 + x1652
    x2096 = x46 * (x1309 * x1647 + x1659)
    x2097 = x104 * x1917 + x1672 + x1921 * x690
    x2098 = x1677 + x1919 + x1922 * x690
    x2099 = x1648 * x690 + x1688 + x1925 * x690
    x2100 = x104 * (x1655 * x690 + x1693 + x1926 * x690)
    x2101 = x1309 * x1700 + x1707
    x2102 = x46 * (x1309 * x1702 + x1713)
    x2103 = x46 * (x1723 + 3.0 * x1924 + x1932 * x690)
    x2104 = x1727 + x1922 * x92 + x1934 * x690
    x2105 = x46 * (x1736 + 2.0 * x1929 + x1937 * x690)
    x2106 = x1740 + x1927 + x1939 * x690
    x2107 = x46 * (x1749 + x1931 + x1941 * x690)
    x2108 = x1703 * x690 + x1753 + x1943 * x690
    x2109 = x46 * (x1309 * x1760 + x1762)
    x2110 = x1309 * x1764 + x1768
    x2111 = x1775 + 4.0 * x1933 + x1947 * x690
    x2112 = x1781 + 3.0 * x1938 + x1949 * x690
    x2113 = x1787 + 2.0 * x1942 + x1951 * x690
    x2114 = x1793 + x1945 + x1953 * x690
    x2115 = x1309 * x1797 + x1800
    x2116 = x1959 + x1964 * x690
    x2117 = x1963 + x1965 * x690
    x2118 = 2.0 * x1966 + x1980 * x690
    x2119 = x1970 + x1981 * x690
    x2120 = x1974 + x1986 * x690
    x2121 = x104 * (x1978 + x1987 * x690)
    x2122 = x46 * (3.0 * x1985 + x2000 * x690)
    x2123 = 3.0 * x1982 + x2002 * x690
    x2124 = 2.0 * x1991
    x2125 = x46 * (x2005 * x690 + x2124)
    x2126 = x1988 + x2007 * x690
    x2127 = x46 * (x1999 + x2010 * x690)
    x2128 = x1995 + x2012 * x690
    x2129 = x2016 * x690
    x2130 = 4.0 * x2001 + x2021 * x690
    x2131 = 3.0 * x2006 + x2024 * x690
    x2132 = 2.0 * x2011 + x2027 * x690
    x2133 = x2016 + x2030 * x690
    x2134 = x129 * x1610 + x1956 * x979
    x2135 = x46 * (x129 * x1617 + x1958 * x979)
    x2136 = x129 * x1630 + x1964 * x979
    x2137 = x46 * (x129 * x1636 + x1965 * x979)
    x2138 = x129 * x1651 + x1959 + x1971 * x979
    x2139 = x46 * (x129 * x1658 + x1963 + x1973 * x979)
    x2140 = x129 * x1671 + x1980 * x979
    x2141 = x46 * (x129 * x1676 + x1981 * x979)
    x2142 = x129 * x1687 + x1966 + x1986 * x979
    x2143 = x104 * (x129 * x1692 + x1969 + x1987 * x979)
    x2144 = x129 * x1706 + 2.0 * x1974 + x1992 * x979
    x2145 = x46 * (x129 * x1712 + x1979 + x1994 * x979)
    x2146 = x46 * (x129 * x1722 + x2000 * x979)
    x2147 = x129 * x1726 + x2002 * x979
    x2148 = x46 * (x129 * x1735 + x1985 + x2005 * x979)
    x2149 = x129 * x1739 + x1982 + x2007 * x979
    x2150 = x46 * (x129 * x1748 + x2010 * x979 + x2124)
    x2151 = x129 * x1752 + x1988 + x2012 * x979
    x2152 = x46 * (x129 * x1761 + 3.0 * x1999 + x2015 * x979)
    x2153 = x129 * x1767 + 3.0 * x1995 + x2017 * x979
    x2154 = x129 * x1774 + x2021 * x979
    x2155 = x129 * x1780 + x2001 + x2024 * x979
    x2156 = x129 * x1786 + 2.0 * x2006 + x2027 * x979
    x2157 = x129 * x1792 + 3.0 * x2011 + x2030 * x979
    x2158 = x129 * x1799 + 4.0 * x2016 + x2033 * x979
    x2159 = x690 * x718
    x2160 = x690 * x757

    # 315 item(s)
    result[0, 0] = numpy.sum(
        x149
        * (
            x23
            * (
                x128 * x3
                + x137 * x83
                + x4 * (x126 * x3 + 4.0 * x131 - x133 * x52)
                - x52 * (x121 * x3 + x143 * x83 + x4 * (x133 - x52 * x66))
            )
            + x3
            * (
                x112 * x83
                - x129 * (x121 * x52 - x128)
                + x3
                * (
                    x3 * (x3 * (x73 + x75) + x78 + x82 * x83)
                    - x4 * (x51 * x52 - x66)
                    + x83 * x93
                )
            )
            + x83
            * (
                x112 * x3
                + x129 * (x137 - x143 * x52)
                + x92
                * (
                    x104 * (x110 * x3 + x146 + x147)
                    + x111 * x3
                    + x4 * (x136 - x142 * x52)
                )
            )
        )
    )
    result[0, 1] = numpy.sum(
        x239
        * (
            x23
            * (
                x219 * x3
                + x230 * x92
                + x4 * (x217 * x3 + 3.0 * x222 - x224 * x52)
                - x52 * (x212 * x3 + x235 * x92 - x4 * (x174 * x52 - x224))
            )
            + x3
            * (
                -x129 * (x212 * x52 - x219)
                + x209 * x92
                + x3
                * (
                    x196 * x92
                    + x3 * (x183 + x188 * x92 + x3 * (x179 + x181))
                    - x4 * (x164 * x52 - x174)
                )
            )
            + x92
            * (
                x104 * (x208 * x3 + x4 * (x229 - x234 * x52) + x46 * (x233 * x238 + x237))
                + x129 * (x230 - x235 * x52)
                + x209 * x3
            )
        )
    )
    result[0, 2] = numpy.sum(
        x239
        * (
            x23
            * (
                x3 * x300
                + x311 * x92
                + x4 * (x298 * x3 + 3.0 * x304 - x306 * x52)
                - x52 * (x294 * x3 + x316 * x92 - x4 * (x261 * x52 - x306))
            )
            + x3
            * (
                -x129 * (x294 * x52 - x300)
                + x291 * x92
                + x3
                * (
                    x280 * x92
                    + x3 * (x269 + x273 * x92 + x3 * (x265 + x267))
                    - x4 * (x252 * x52 - x261)
                )
            )
            + x92
            * (
                x104 * (x290 * x3 + x4 * (x310 - x315 * x52) + x46 * (x238 * x314 + x317))
                + x129 * (x311 - x316 * x52)
                + x291 * x3
            )
        )
    )
    result[0, 3] = numpy.sum(
        x378
        * (
            x104
            * (
                x129 * (x375 - x377 * x52)
                + x3 * x361
                + x3 * x46 * (x360 + x4 * (x330 - x358))
            )
            + x23
            * (
                x104 * x375
                + x3 * x370
                + x4 * (x104 * x373 + x3 * x368 - x372 * x52)
                - x52 * (x104 * x377 + x3 * x366 - x4 * (x335 * x52 - x372))
            )
            + x3
            * (
                x104 * x361
                - x129 * (x366 * x52 - x370)
                + x3
                * (
                    x104 * x350
                    + x3 * (x104 * x346 + x3 * (x340 + x342) + x344)
                    - x4 * (x329 * x52 - x335)
                )
            )
        )
    )
    result[0, 4] = numpy.sum(
        x426
        * (
            x104
            * (
                x129 * (x423 - x425 * x52)
                + x3 * x409
                + x3
                * x46
                * (x4 * (2.0 * x10 * x151 * x17 * x241 * x32 * x4 * x6 - x406) + x408)
            )
            + x23
            * (
                x104 * x423
                + x3 * x417
                + x4 * (x3 * x415 - x419 * x52 + x421)
                - x52 * (x104 * x425 + x3 * x413 - x4 * (x390 * x52 - x419))
            )
            + x3
            * (
                x104 * x409
                - x129 * (x413 * x52 - x417)
                + x3
                * (
                    x104 * x401
                    + x3 * (x104 * x399 + x3 * (x3 * x392 + x393) + x394)
                    - x4 * (x383 * x52 - x390)
                )
            )
        )
    )
    result[0, 5] = numpy.sum(
        x378
        * (
            x104
            * (
                x129 * (x486 - x489 * x52)
                + x3 * x46 * (x4 * (x439 - x466) + x468)
                + x3 * x469
            )
            + x23
            * (
                x104 * x486
                + x3 * x480
                + x4 * (x104 * x483 + x3 * x478 - x482 * x52)
                - x52 * (x104 * x489 + x3 * x475 - x4 * (x445 * x52 - x482))
            )
            + x3
            * (
                x104 * x469
                - x129 * (x475 * x52 - x480)
                + x3
                * (
                    x104 * x458
                    + x3 * (x104 * x455 + x3 * (x104 * x451 + x450) + x453)
                    - x4 * (x438 * x52 - x445)
                )
            )
        )
    )
    result[0, 6] = numpy.sum(
        x239
        * (
            x23
            * (
                x3 * x524
                + x4 * (x3 * x526 - x52 * x528 + x525)
                + x46 * x504
                - x52 * (x3 * x519 - x4 * (x511 * x52 - x528) + x46 * x501)
            )
            + x3
            * (
                -x129 * (x519 * x52 - x524)
                + x3
                * (
                    x3 * (x3 * (x512 + x514) + x515 + x516)
                    - x4 * (x507 * x52 - x511)
                    + x46 * x497
                )
                + x46 * x498
            )
            - x46 * (x129 * (x501 * x52 - x504) - x3 * x498)
        )
    )
    result[0, 7] = numpy.sum(
        x426
        * (
            x23
            * (
                x3 * x555
                + x4 * (x3 * x557 - x52 * x559 + x556)
                + x46 * x541
                - x52 * (x3 * x551 - x4 * (x52 * x546 - x559) + x46 * x538)
            )
            + x3
            * (
                -x129 * (x52 * x551 - x555)
                + x3
                * (
                    x3 * (x3 * x547 + x3 * (x108 * x548 + x547) + x549)
                    - x4 * (x52 * x543 - x546)
                    + x46 * x534
                )
                + x46 * x535
            )
            - x46 * (x129 * (x52 * x538 - x541) - x3 * x535)
        )
    )
    result[0, 8] = numpy.sum(
        x426
        * (
            x23
            * (
                x3 * x583
                + x4 * (x3 * x587 + x46 * x584 - x52 * x586)
                + x46 * x570
                - x52 * (x3 * x580 - x4 * (x52 * x574 - x586) + x46 * x568)
            )
            + x3
            * (
                -x129 * (x52 * x580 - x583)
                + x3
                * (
                    x3 * (x3 * x575 + x3 * (x108 * x576 + x575) + x577)
                    - x4 * (x52 * x572 - x574)
                    + x46 * x564
                )
                + x46 * x565
            )
            - x46 * (x129 * (x52 * x568 - x570) - x3 * x565)
        )
    )
    result[0, 9] = numpy.sum(
        x239
        * (
            x23
            * (
                x3 * x623
                + x4 * (x3 * x625 - x52 * x627 + x624)
                + x46 * x603
                - x52 * (x3 * x618 - x4 * (x52 * x610 - x627) + x46 * x600)
            )
            + x3
            * (
                -x129 * (x52 * x618 - x623)
                + x3
                * (
                    x3 * (x3 * (x611 + x613) + x614 + x615)
                    - x4 * (x52 * x606 - x610)
                    + x46 * x595
                )
                + x46 * x596
            )
            - x46 * (x129 * (x52 * x600 - x603) - x3 * x596)
        )
    )
    result[0, 10] = numpy.sum(
        x149
        * x3
        * (
            -x129 * (x52 * x638 - x642)
            - x23
            * (x4 * (x52 * x636 - x640) + x52 * (-x4 * (x52 * x628 - x636) + x638) - x642)
            + x3 * (x3 * (x633 + x634) + x4 * (x3 * x628 - x631))
        )
    )
    result[0, 11] = numpy.sum(
        x239
        * x3
        * (
            -x129 * (x52 * x650 - x653)
            - x23
            * (x4 * (x52 * x648 - x651) + x52 * (-x4 * (x52 * x643 - x648) + x650) - x653)
            + x3**2 * (x238 * x646 + x4 * (x643 - x645) + x647)
        )
    )
    result[0, 12] = numpy.sum(
        x3
        * x378
        * (
            -x129 * (x52 * x661 - x664)
            - x23
            * (x4 * (x52 * x659 - x662) + x52 * (-x4 * (x52 * x654 - x659) + x661) - x664)
            + x3**2 * (x238 * x657 + x4 * (x654 - x656) + x658)
        )
    )
    result[0, 13] = numpy.sum(
        x239
        * (
            x23
            * (
                x151 * x4 * (x3 * x598 - x52 * x621)
                - x3 * x52 * (x4 * (x151 * x588 - x52 * x669) + x671)
                + x3 * x674
            )
            - x3
            * (
                x129 * (x52 * x671 - x674)
                - x3**2 * (x238 * x668 + x4 * (x151 * x589 - x666) + x667)
            )
        )
    )
    result[0, 14] = numpy.sum(
        x149
        * x3
        * (
            -x129 * (x52 * x685 - x689)
            - x23
            * (x4 * (x52 * x683 - x687) + x52 * (-x4 * (x52 * x675 - x683) + x685) - x689)
            + x3 * (x3 * (x680 + x681) + x4 * (x3 * x675 - x678))
        )
    )
    result[1, 0] = numpy.sum(
        0.5
        * x718
        * (
            x129
            * (
                2.0 * x3 * x696
                - x52 * (2.0 * x3 * x693 + 2.0 * x711 * x83 + x717)
                + 2.0 * x713 * x83
                + x716
            )
            + x3
            * (
                x3 * (2.0 * x3 * (x697 + x698) + x700 + 2.0 * x703 * x83)
                - 2.0 * x4 * (x52 * x693 - x696)
                + 2.0 * x710 * x83
            )
            + x83
            * (
                2.0 * x3 * x710
                - 2.0 * x4 * (x52 * x711 - x713)
                + x92 * (2.0 * x104 * (x109 * x690 + x715) + 2.0 * x3 * x709 + x714)
            )
        )
    )
    result[1, 1] = numpy.sum(
        0.5
        * x757
        * (
            x129
            * (
                2.0 * x3 * x728
                - x52 * (2.0 * x3 * x723 + 2.0 * x748 * x92 + x756)
                + 2.0 * x752 * x92
                + x755
            )
            + x3
            * (
                x3 * (2.0 * x3 * (x3 * x730 + x732 * x92) + x733 + 2.0 * x738 * x92)
                - 2.0 * x4 * (x52 * x723 - x728)
                + 2.0 * x744 * x92
            )
            + x92
            * (
                x104 * (2.0 * x3 * x742 + 2.0 * x3 * x743 + x753)
                + 2.0 * x3 * x744
                - 2.0 * x4 * (x52 * x748 - x752)
            )
        )
    )
    result[1, 2] = numpy.sum(
        0.5
        * x757
        * (
            x129
            * (
                2.0 * x3 * x760
                - x52 * (2.0 * x3 * x758 + 2.0 * x767 * x92 + x773)
                + 2.0 * x770 * x92
                + x772
            )
            + x3
            * (
                x3 * (2.0 * x3 * x690 * (x265 + x267) + x761 + 2.0 * x762 * x92)
                - 2.0 * x4 * (x52 * x758 - x760)
                + 2.0 * x766 * x92
            )
            + x92
            * (
                x104 * (2.0 * x289 * x690 + 2.0 * x3 * x765 + x771)
                + 2.0 * x3 * x766
                - 2.0 * x4 * (x52 * x767 - x770)
            )
        )
    )
    result[1, 3] = numpy.sum(
        0.5
        * x810
        * (
            x104
            * (
                2.0 * x3 * x798
                - 2.0 * x4 * (x52 * x802 - x805)
                + x46 * (2.0 * x238 * x794 + x806)
            )
            + x129
            * (
                2.0 * x104 * x805
                + 2.0 * x3 * x785
                - x52 * (2.0 * x104 * x802 + 2.0 * x3 * x780 + x809)
                + x808
            )
            + x3
            * (
                2.0 * x104 * x798
                + x3 * (2.0 * x104 * x796 + 2.0 * x3 * (x104 * x790 + x3 * x788) + x791)
                - 2.0 * x4 * (x52 * x780 - x785)
            )
        )
    )
    result[1, 4] = numpy.sum(
        0.5
        * x840
        * (
            x104
            * (
                2.0 * x3 * x828
                - 2.0 * x4 * (x52 * x831 - x835)
                + x46 * (2.0 * x238 * x824 + x836)
            )
            + x129
            * (
                2.0 * x104 * x835
                + 2.0 * x3 * x818
                - x52 * (2.0 * x104 * x831 + 2.0 * x3 * x814 + x839)
                + x838
            )
            + x3
            * (
                2.0 * x104 * x828
                + x3 * (2.0 * x104 * x826 + 2.0 * x3 * (x3 * x819 + x822) + x823)
                - 2.0 * x4 * (x52 * x814 - x818)
            )
        )
    )
    result[1, 5] = numpy.sum(
        0.5
        * x810
        * (
            x104
            * (
                2.0 * x3 * x851
                - 2.0 * x4 * (x52 * x853 - x855)
                + x46 * (2.0 * x465 * x690 + x856)
            )
            + x129
            * (
                2.0 * x104 * x855
                + 2.0 * x3 * x844
                - x52 * (2.0 * x104 * x853 + 2.0 * x3 * x842 + x858)
                + x857
            )
            + x3
            * (
                2.0 * x104 * x851
                + x3 * (2.0 * x104 * x848 + 2.0 * x3 * (x104 * x845 + x450 * x690) + x846)
                - 2.0 * x4 * (x52 * x842 - x844)
            )
        )
    )
    result[1, 6] = numpy.sum(
        0.5
        * x757
        * (
            x129
            * (
                2.0 * x3 * x46 * x861
                + 2.0 * x3 * x878
                - x52 * (2.0 * x3 * x871 + 2.0 * x3 * x874 + x886)
                + x885
            )
            + 2.0 * x3 * x46 * (x4 * (x861 - x865) + x870)
            + x3
            * (
                x3 * (2.0 * x3 * x879 + 2.0 * x3 * (x3 * x882 + x879) + x883)
                - 2.0 * x4 * (x52 * x874 - x878)
                + 2.0 * x46 * x870
            )
        )
    )
    result[1, 7] = numpy.sum(
        0.5
        * x840
        * (
            x129
            * (
                2.0 * x3 * x46 * x888
                + 2.0 * x3 * x902
                - x52 * (2.0 * x3 * x897 + 2.0 * x3 * x899 + x908)
                + x907
            )
            + 2.0 * x3 * x46 * (x4 * (x888 - x891) + x896)
            + x3
            * (
                x3 * (2.0 * x3 * x903 + 2.0 * x3 * (x3 * x904 + x903) + x905)
                - 2.0 * x4 * (x52 * x899 - x902)
                + 2.0 * x46 * x896
            )
        )
    )
    result[1, 8] = numpy.sum(
        0.5
        * x840
        * (
            x129
            * (
                2.0 * x3 * x46 * x909
                + 2.0 * x3 * x920
                - x52 * (2.0 * x3 * x915 + 2.0 * x3 * x917 + x926)
                + x925
            )
            + 2.0 * x3 * x46 * (x4 * (x909 - x911) + x914)
            + x3
            * (
                x3 * (2.0 * x3 * x921 + 2.0 * x3 * (x3 * x922 + x921) + x923)
                - 2.0 * x4 * (x52 * x917 - x920)
                + 2.0 * x46 * x914
            )
        )
    )
    result[1, 9] = numpy.sum(
        0.5
        * x757
        * (
            x129
            * (
                2.0 * x3 * x932
                - x52 * (2.0 * x3 * x930 + 2.0 * x617 * x690 + x936)
                + 2.0 * x622 * x690
                + x935
            )
            + x3
            * (
                x3 * (2.0 * x3 * (x613 * x690 + x933) + 2.0 * x615 * x690 + x934)
                - 2.0 * x4 * (x52 * x930 - x932)
                + 2.0 * x46 * x928
            )
            + 2.0 * x46 * (x3 * x928 + x4 * x690 * (x3 * x588 - x591))
        )
    )
    result[1, 10] = numpy.sum(
        0.5
        * x718
        * (
            x129 * (2.0 * x238 * x939 - x52 * (2.0 * x238 * x942 + x950) + x949)
            + x3**2 * (2.0 * x238 * x946 + 2.0 * x4 * (x939 - x943) + x947)
        )
    )
    result[1, 11] = numpy.sum(
        0.5
        * x757
        * (
            x129 * (2.0 * x238 * x951 - x52 * (2.0 * x238 * x952 + x958) + x957)
            + x3**2 * (2.0 * x238 * x954 + 2.0 * x4 * (x951 - x953) + x955)
        )
    )
    result[1, 12] = numpy.sum(
        0.5
        * x810
        * (
            x129 * (2.0 * x238 * x959 - x52 * (2.0 * x238 * x960 + x966) + x965)
            + x3**2 * (2.0 * x238 * x962 + 2.0 * x4 * (x959 - x961) + x963)
        )
    )
    result[1, 13] = numpy.sum(
        0.5
        * x757
        * (
            x129 * (2.0 * x238 * x967 - x52 * (2.0 * x238 * x968 + x975) + x974)
            + x3**2 * (2.0 * x238 * x970 + 2.0 * x4 * (x967 - x969) + x971)
        )
    )
    result[1, 14] = numpy.sum(
        0.5
        * x718
        * (
            x129 * (-x52 * (2.0 * x682 * x690 + x978) + 2.0 * x686 * x690 + x977)
            + x3
            * (x3 * (2.0 * x680 * x690 + x976) + 2.0 * x4 * x690 * (x3 * x675 - x678))
        )
    )
    result[2, 0] = numpy.sum(
        x718
        * (
            x129 * (x3 * x981 - x52 * (x3 * x980 + x83 * x990 + x999) + x83 * x991 + x996)
            + x3
            * (
                x3 * (x3 * x979 * (x73 + x75) + x83 * x984 + x983)
                - x4 * (x52 * x980 - x981)
                + x83 * x989
            )
            + x83
            * (
                x3 * x989
                - x4 * (x52 * x990 - x991)
                + x92 * (x104 * (x109 * x979 + x994) + x3 * x988 + x993)
            )
        )
    )
    result[2, 1] = numpy.sum(
        x757
        * (
            x129
            * (
                x1002 * x3
                + x1017 * x92
                + x1021
                - x52 * (x1000 * x3 + x1015 * x92 + x1023)
            )
            + x3
            * (
                x1013 * x92
                + x3 * (x1004 + x1008 * x92 + x3 * x979 * (x179 + x181))
                - x4 * (x1000 * x52 - x1002)
            )
            + x92
            * (
                x1013 * x3
                + x104 * (x1012 * x3 + x1019 + x207 * x979)
                - x4 * (x1015 * x52 - x1017)
            )
        )
    )
    result[2, 2] = numpy.sum(
        x757
        * (
            x129
            * (
                x1034 * x3
                + x1063 * x92
                + x1068
                - x52 * (x1029 * x3 + x1059 * x92 + x1070)
            )
            + x3
            * (
                x1055 * x92
                + x3 * (x1042 + x1047 * x92 + x3 * (x1036 + x1039))
                - x4 * (x1029 * x52 - x1034)
            )
            + x92
            * (
                x104 * (x1052 * x3 + x1054 * x3 + x1065)
                + x1055 * x3
                - x4 * (x1059 * x52 - x1063)
            )
        )
    )
    result[2, 3] = numpy.sum(
        x810
        * (
            x104 * (x1082 * x3 - x4 * (x1084 * x52 - x1086) + x46 * (x1088 + x357 * x979))
            + x129
            * (
                x104 * x1086
                + x1075 * x3
                + x1090
                - x52 * (x104 * x1084 + x1073 * x3 + x1092)
            )
            + x3
            * (
                x104 * x1082
                + x3 * (x104 * x1078 + x1077 + x3 * x979 * (x340 + x342))
                - x4 * (x1073 * x52 - x1075)
            )
        )
    )
    result[2, 4] = numpy.sum(
        x840
        * (
            x104
            * (x1111 * x3 - x4 * (x1114 * x52 - x1117) + x46 * (x1106 * x238 + x1119))
            + x129
            * (
                x104 * x1117
                + x1100 * x3
                + x1122
                - x52 * (x104 * x1114 + x1096 * x3 + x1124)
            )
            + x3
            * (
                x104 * x1111
                + x3 * (x104 * x1108 + x1105 + x3 * (x1101 * x3 + x1103))
                - x4 * (x1096 * x52 - x1100)
            )
        )
    )
    result[2, 5] = numpy.sum(
        x810
        * (
            x104
            * (x1154 * x3 - x4 * (x1157 * x52 - x1160) + x46 * (x1146 * x238 + x1162))
            + x129
            * (
                x104 * x1160
                + x1137 * x3
                + x1165
                - x52 * (x104 * x1157 + x1132 * x3 + x1167)
            )
            + x3
            * (
                x104 * x1154
                + x3 * (x104 * x1149 + x1145 + x3 * (x1139 + x1142))
                - x4 * (x1132 * x52 - x1137)
            )
        )
    )
    result[2, 6] = numpy.sum(
        x757
        * (
            x129
            * (
                x1174 * x3
                + x1179
                - x52 * (x1172 * x3 + x1181 + x518 * x979)
                + x523 * x979
            )
            + x3
            * (
                x1170 * x46
                + x3 * (x1177 + x3 * (x1175 + x514 * x979) + x516 * x979)
                - x4 * (x1172 * x52 - x1174)
            )
            + x46 * (x1170 * x3 + x4 * x979 * (x3 * x490 - x493))
        )
    )
    result[2, 7] = numpy.sum(
        x840
        * (
            x129
            * (
                x1182 * x3 * x46
                + x1195 * x3
                + x1203
                - x52 * (x1190 * x3 + x1192 * x3 + x1205)
            )
            + x3 * x46 * (x1189 + x4 * (x1182 - x1184))
            + x3
            * (
                x1189 * x46
                + x3 * (x1196 * x3 + x1200 + x3 * (x1196 + x1198 * x3))
                - x4 * (x1192 * x52 - x1195)
            )
        )
    )
    result[2, 8] = numpy.sum(
        x840
        * (
            x129
            * (
                x1206 * x3 * x46
                + x1218 * x3
                + x1225
                - x52 * (x1213 * x3 + x1215 * x3 + x1227)
            )
            + x3 * x46 * (x1212 + x4 * (x1206 - x1208))
            + x3
            * (
                x1212 * x46
                + x3 * (x1219 * x3 + x1222 + x3 * (x1219 + x1220 * x3))
                - x4 * (x1215 * x52 - x1218)
            )
        )
    )
    result[2, 9] = numpy.sum(
        x757
        * (
            x129
            * (
                x1228 * x3 * x46
                + x1244 * x3
                + x1254
                - x52 * (x1237 * x3 + x1240 * x3 + x1256)
            )
            + x3 * x46 * (x1236 + x4 * (x1228 - x1230))
            + x3
            * (
                x1236 * x46
                + x3 * (x1250 + x1251 + x3 * (x1245 + x1247))
                - x4 * (x1240 * x52 - x1244)
            )
        )
    )
    result[2, 10] = numpy.sum(
        x718
        * (
            x129 * (x1260 - x52 * (x1262 + x635 * x979) + x639 * x979)
            + x3 * (x3 * (x1258 + x633 * x979) + x4 * x979 * (x3 * x628 - x631))
        )
    )
    result[2, 11] = numpy.sum(
        x757
        * (
            x129 * (x1263 * x238 + x1271 - x52 * (x1264 * x238 + x1273))
            + x3**2 * (x1266 * x238 + x1268 + x4 * (x1263 - x1265))
        )
    )
    result[2, 12] = numpy.sum(
        x810
        * (
            x129 * (x1274 * x238 + x1282 - x52 * (x1275 * x238 + x1284))
            + x3**2 * (x1277 * x238 + x1279 + x4 * (x1274 - x1276))
        )
    )
    result[2, 13] = numpy.sum(
        x757
        * (
            x129 * (x1285 * x238 + x1293 - x52 * (x1286 * x238 + x1295))
            + x3**2 * (x1288 * x238 + x1290 + x4 * (x1285 - x1287))
        )
    )
    result[2, 14] = numpy.sum(
        x718
        * (
            x129 * (x1296 * x238 + x1306 - x52 * (x1297 * x238 + x1308))
            + x3 * (x3 * (x1301 + x1303) + x4 * (x1296 * x3 - x1299))
        )
    )
    result[3, 0] = numpy.sum(
        x1330
        * (
            x3 * (x1322 * x77 + x1326 * x83 + x3 * (x1318 * x3 + 4.0 * x1321))
            + x4 * (x1310 * x3 + 4.0 * x1312 - x52 * (x1314 * x3 + 4.0 * x1316))
            + x83 * (x1326 * x3 + x1327 * x77 + x92 * (x1324 * x3 + 2.0 * x1329))
        )
    )
    result[3, 1] = numpy.sum(
        x1348
        * (
            x3 * (x1341 * x77 + x1345 * x92 + x3 * (x1338 * x3 + x1340 * x92))
            + x4 * (x1332 * x3 + x1333 * x92 - x52 * (x1335 * x3 + x1336 * x92))
            + x92 * (x104 * (x1343 * x3 + x1347) + x1345 * x3 + x1346 * x77)
        )
    )
    result[3, 2] = numpy.sum(
        x1348
        * (
            x3 * (x1358 * x77 + x1361 * x92 + x3 * (x1355 * x3 + 3.0 * x1357))
            + x4 * (x1349 * x3 + 3.0 * x1351 - x52 * (x1352 * x3 + 3.0 * x1354))
            + x92 * (x104 * (x1359 * x3 + x1363) + x1361 * x3 + x1362 * x77)
        )
    )
    result[3, 3] = numpy.sum(
        x1375
        * (
            x104 * (x1372 * x3 + x1373 * x3 + x1374 * x77)
            + x3 * (x104 * x1373 + x1370 * x77 + x3 * (x104 * x1369 + x1368 * x3))
            + x4 * (x104 * x1365 + x1364 * x3 - x52 * (x104 * x1367 + x1366 * x3))
        )
    )
    result[3, 4] = numpy.sum(
        x1390
        * (
            x104 * (x1387 * x3 + x1388 * x3 + x1389 * x77)
            + x3 * (x104 * x1388 + x1385 * x77 + x3 * (x1382 * x3 + x1384))
            + x4 * (x1376 * x3 + x1378 - x52 * (x1379 * x3 + x1381))
        )
    )
    result[3, 5] = numpy.sum(
        x1375
        * (
            x104 * (x1403 * x3 + x1404 * x3 + x1405 * x77)
            + x3 * (x104 * x1404 + x1401 * x77 + x3 * (x1398 * x3 + 2.0 * x1400))
            + x4 * (x1391 * x3 + 2.0 * x1393 - x52 * (x1394 * x3 + 2.0 * x1397))
        )
    )
    result[3, 6] = numpy.sum(
        x1348
        * (
            x3 * (x1414 * x3 + x1416 * x77 + x3 * (x1414 + x1415 * x3))
            + x4 * (x1407 + x1408 * x3 - x52 * (x1410 + x1411 * x3))
            + x46 * (x1412 * x238 + x1413 * x77)
        )
    )
    result[3, 7] = numpy.sum(
        x1390
        * (
            x3 * (x1425 * x3 + x1427 * x77 + x3 * (x1425 + x1426 * x3))
            + x4 * (x1418 + x1419 * x3 - x52 * (x1421 + x1422 * x3))
            + x46 * (x1423 * x238 + x1424 * x77)
        )
    )
    result[3, 8] = numpy.sum(
        x1390
        * (
            x3 * (x1436 * x3 + x1438 * x77 + x3 * (x1436 + x1437 * x3))
            + x4 * (x1429 + x1430 * x3 - x52 * (x1432 + x1433 * x3))
            + x46 * (x1434 * x238 + x1435 * x77)
        )
    )
    result[3, 9] = numpy.sum(
        x1348
        * (
            x3 * (x1448 * x3 + x1450 * x77 + x3 * (x1448 + x1449 * x3))
            + x4 * (x1440 + x1441 * x3 - x52 * (x1443 + x1445 * x3))
            + x46 * (x1446 * x238 + x1447 * x77)
        )
    )
    result[3, 10] = numpy.sum(
        x1330 * x3 * (x1454 * x238 + x1455 * x77 + x4 * (x1451 - x1453))
    )
    result[3, 11] = numpy.sum(
        x1348 * x3 * (x1459 * x238 + x1460 * x77 + x4 * (x1456 - x1458))
    )
    result[3, 12] = numpy.sum(
        x1375 * x3 * (x1464 * x238 + x1465 * x77 + x4 * (x1461 - x1463))
    )
    result[3, 13] = numpy.sum(
        x1348 * x3 * (x1469 * x238 + x1470 * x77 + x4 * (x1466 - x1468))
    )
    result[3, 14] = numpy.sum(
        x1330 * x3 * (x1474 * x238 + x1475 * x77 + x4 * (x1471 - x1473))
    )
    result[4, 0] = numpy.sum(
        0.5
        * x757
        * (
            x3 * (x1476 + 2.0 * x1477 * x83 + 2.0 * x3 * x979 * (x697 + x698))
            + 2.0 * x4 * x979 * (-x52 * (x691 + x692) + x694 + x695 * x83)
            + x83 * (2.0 * x1477 * x3 + x1478 + 2.0 * x92 * x979 * (x706 + x708))
        )
    )
    result[4, 1] = numpy.sum(
        0.5
        * x1498
        * (
            x3 * (x1491 + 2.0 * x1495 * x92 + 2.0 * x3 * (x1488 * x3 + x1490 * x92))
            + 2.0 * x4 * (x1480 * x3 + x1482 * x92 - x52 * (x1484 * x3 + x1486 * x92))
            + x92 * (2.0 * x104 * (x1493 * x3 + x1497) + 2.0 * x1495 * x3 + x1496)
        )
    )
    result[4, 2] = numpy.sum(
        0.5
        * x1498
        * (
            x3 * (x1500 + 2.0 * x1501 * x92 + 2.0 * x3 * x690 * (x1036 + x1039))
            + 2.0 * x4 * (x1031 * x690 + x1499 * x92 - x52 * x690 * (x1025 + x1028))
            + x92 * (2.0 * x104 * (x1053 * x690 + x1503) + 2.0 * x1501 * x3 + x1502)
        )
    )
    result[4, 3] = numpy.sum(
        0.5
        * x1390
        * (
            x104 * (2.0 * x1512 * x3 + 2.0 * x1513 * x3 + x1514)
            + x3 * (2.0 * x104 * x1513 + x1510 + 2.0 * x3 * (x104 * x1509 + x1508 * x3))
            + 2.0 * x4 * (x104 * x1505 + x1504 * x3 - x52 * (x104 * x1507 + x1506 * x3))
        )
    )
    result[4, 4] = numpy.sum(
        0.5
        * x1529
        * (
            x104 * (2.0 * x1526 * x3 + 2.0 * x1527 * x3 + x1528)
            + x3 * (2.0 * x104 * x1527 + x1524 + 2.0 * x3 * (x1521 * x3 + x1523))
            + 2.0 * x4 * (x1515 * x3 + x1517 - x52 * (x1518 * x3 + x1520))
        )
    )
    result[4, 5] = numpy.sum(
        0.5
        * x1390
        * (
            x104 * (2.0 * x1153 * x690 + 2.0 * x1533 * x3 + x1534)
            + x3 * (2.0 * x104 * x1533 + x1531 + 2.0 * x3 * x690 * (x1139 + x1142))
            + 2.0 * x4 * (x104 * x1530 + x1134 * x690 - x52 * x690 * (x1126 + x1131))
        )
    )
    result[4, 6] = numpy.sum(
        0.5
        * x1498
        * (
            x3 * (2.0 * x1543 * x3 + x1545 + 2.0 * x3 * (x1543 + x1544 * x3))
            + 2.0 * x4 * (x1536 + x1537 * x3 - x52 * (x1539 + x1540 * x3))
            + x46 * (2.0 * x1541 * x238 + x1542)
        )
    )
    result[4, 7] = numpy.sum(
        0.5
        * x1529
        * (
            x3 * (2.0 * x1557 * x3 + x1559 + 2.0 * x3 * (x1557 + x1558 * x3))
            + 2.0 * x4 * (x1548 + x1549 * x3 - x52 * (x1552 + x1553 * x3))
            + x46 * (2.0 * x1555 * x238 + x1556)
        )
    )
    result[4, 8] = numpy.sum(
        0.5
        * x1529
        * (
            x3 * (2.0 * x1568 * x3 + x1570 + 2.0 * x3 * (x1568 + x1569 * x3))
            + 2.0 * x4 * (x1561 + x1562 * x3 - x52 * (x1564 + x1565 * x3))
            + x46 * (2.0 * x1566 * x238 + x1567)
        )
    )
    result[4, 9] = numpy.sum(
        0.5
        * x1498
        * (
            x3 * (2.0 * x1251 * x690 + x1575 + 2.0 * x3 * (x1247 * x690 + x1574))
            + 2.0 * x4 * (x1243 * x690 + x1571 - x52 * (x1239 * x690 + x1572))
            + x46 * (2.0 * x1233 * x690 + x1573)
        )
    )
    result[4, 10] = numpy.sum(
        0.5 * x3 * x757 * (2.0 * x1579 * x238 + x1580 + 2.0 * x4 * (x1576 - x1578))
    )
    result[4, 11] = numpy.sum(
        0.5 * x1498 * x3 * (2.0 * x1584 * x238 + x1585 + 2.0 * x4 * (x1581 - x1583))
    )
    result[4, 12] = numpy.sum(
        0.5 * x1390 * x3 * (2.0 * x1589 * x238 + x1590 + 2.0 * x4 * (x1586 - x1588))
    )
    result[4, 13] = numpy.sum(
        0.5 * x1498 * x3 * (2.0 * x1594 * x238 + x1595 + 2.0 * x4 * (x1591 - x1593))
    )
    result[4, 14] = numpy.sum(
        0.5
        * x757
        * (x3 * (2.0 * x1301 * x690 + x1596) + 2.0 * x4 * x690 * (x1296 * x3 - x1299))
    )
    result[5, 0] = numpy.sum(
        x1330
        * (
            x3 * (x1611 + x1616 * x83 + x3 * (x1605 + x1608))
            + x4 * (x1598 * x3 + 4.0 * x1600 - x52 * (x1601 * x3 + 4.0 * x1603))
            + x83 * (x1616 * x3 + x1618 + x92 * (x1613 * x3 + 2.0 * x1620))
        )
    )
    result[5, 1] = numpy.sum(
        x1348
        * (
            x3 * (x1631 + x1635 * x92 + x3 * (x1627 * x3 + 3.0 * x1629))
            + x4 * (x1621 * x3 + 3.0 * x1623 - x52 * (x1624 * x3 + 3.0 * x1626))
            + x92 * (x104 * (x1632 * x3 + x1638) + x1635 * x3 + x1637)
        )
    )
    result[5, 2] = numpy.sum(
        x1348
        * (
            x3 * (x1652 + x1657 * x92 + x3 * (x1646 + x1649))
            + x4 * (x1639 * x3 + 3.0 * x1641 - x52 * (x1642 * x3 + 3.0 * x1644))
            + x92 * (x104 * (x1654 * x3 + x1660) + x1657 * x3 + x1659)
        )
    )
    result[5, 3] = numpy.sum(
        x1375
        * (
            x104 * (x1674 * x3 + x1675 * x3 + x1677)
            + x3 * (x104 * x1675 + x1672 + x3 * (x1668 * x3 + 2.0 * x1670))
            + x4 * (x1661 * x3 + 2.0 * x1663 - x52 * (x1664 * x3 + 2.0 * x1667))
        )
    )
    result[5, 4] = numpy.sum(
        x1390
        * (
            x104 * (x1690 * x3 + x1691 * x3 + x1693)
            + x3 * (x104 * x1691 + x1688 + x3 * (x1684 * x3 + x1686))
            + x4 * (x1678 * x3 + x1680 - x52 * (x1681 * x3 + x1683))
        )
    )
    result[5, 5] = numpy.sum(
        x1375
        * (
            x104 * (x1709 * x3 + x1711 * x3 + x1713)
            + x3 * (x104 * x1711 + x1707 + x3 * (x1701 + x1704))
            + x4 * (x1694 * x3 + 2.0 * x1696 - x52 * (x1697 * x3 + 2.0 * x1699))
        )
    )
    result[5, 6] = numpy.sum(
        x1348
        * (
            x3 * (x1724 * x3 + x1727 + x3 * (x1724 + x1725 * x3))
            + x4 * (x1715 + x1716 * x3 - x52 * (x1718 + x1720 * x3))
            + x46 * (x1721 * x238 + x1723)
        )
    )
    result[5, 7] = numpy.sum(
        x1390
        * (
            x3 * (x1737 * x3 + x1740 + x3 * (x1737 + x1738 * x3))
            + x4 * (x1729 + x1730 * x3 - x52 * (x1732 + x1733 * x3))
            + x46 * (x1734 * x238 + x1736)
        )
    )
    result[5, 8] = numpy.sum(
        x1390
        * (
            x3 * (x1750 * x3 + x1753 + x3 * (x1750 + x1751 * x3))
            + x4 * (x1742 + x1743 * x3 - x52 * (x1745 + x1746 * x3))
            + x46 * (x1747 * x238 + x1749)
        )
    )
    result[5, 9] = numpy.sum(
        x1348
        * (
            x3 * (x1768 + x1769 + x3 * (x1763 + x1765))
            + x4 * (x1755 + x1756 * x3 - x52 * (x1758 + x1759 * x3))
            + x46 * (x1760 * x238 + x1762)
        )
    )
    result[5, 10] = numpy.sum(x1330 * x3 * (x1773 * x238 + x1775 + x4 * (x1770 - x1772)))
    result[5, 11] = numpy.sum(x1348 * x3 * (x1779 * x238 + x1781 + x4 * (x1776 - x1778)))
    result[5, 12] = numpy.sum(x1375 * x3 * (x1785 * x238 + x1787 + x4 * (x1782 - x1784)))
    result[5, 13] = numpy.sum(x1348 * x3 * (x1791 * x238 + x1793 + x4 * (x1788 - x1790)))
    result[5, 14] = numpy.sum(x1330 * x3 * (x1798 + x1800 + x4 * (x1794 - x1796)))
    result[6, 0] = numpy.sum(
        x1806
        * (
            x1804 * x77
            + x3 * (x1801 * x3 + 4.0 * x1803)
            + x83 * (x1802 * x3 + 3.0 * x1805)
        )
    )
    result[6, 1] = numpy.sum(
        x1811
        * (x1809 * x77 + x3 * (x1807 * x3 + x1808 * x92) + x92 * (x1808 * x3 + x1810))
    )
    result[6, 2] = numpy.sum(
        x1811
        * (
            x1815 * x77
            + x3 * (x1812 * x3 + 3.0 * x1814)
            + x92 * (x1813 * x3 + 2.0 * x1816)
        )
    )
    result[6, 3] = numpy.sum(
        x1821
        * (x104 * (x1818 * x3 + x1820) + x1819 * x77 + x3 * (x104 * x1818 + x1817 * x3))
    )
    result[6, 4] = numpy.sum(
        x1827 * (x104 * (x1823 * x3 + x1826) + x1825 * x77 + x3 * (x1822 * x3 + x1824))
    )
    result[6, 5] = numpy.sum(
        x1821
        * (x104 * (x1829 * x3 + x1832) + x1831 * x77 + x3 * (x1828 * x3 + 2.0 * x1830))
    )
    result[6, 6] = numpy.sum(
        x1811 * (x1834 * x3 + x1836 * x77 + x3 * (x1834 + x1835 * x3))
    )
    result[6, 7] = numpy.sum(
        x1827 * (x1838 * x3 + x1840 * x77 + x3 * (x1838 + x1839 * x3))
    )
    result[6, 8] = numpy.sum(
        x1827 * (x1842 * x3 + x1844 * x77 + x3 * (x1842 + x1843 * x3))
    )
    result[6, 9] = numpy.sum(
        x1811 * (x1846 * x3 + x1848 * x77 + x3 * (x1846 + x1847 * x3))
    )
    result[6, 10] = numpy.sum(x1806 * (x1849 * x238 + x1850 * x77))
    result[6, 11] = numpy.sum(x1811 * (x1851 * x238 + x1852 * x77))
    result[6, 12] = numpy.sum(x1821 * (x1853 * x238 + x1854 * x77))
    result[6, 13] = numpy.sum(x1811 * (x1855 * x238 + x1856 * x77))
    result[6, 14] = numpy.sum(x1806 * (x1857 * x238 + x1858 * x77))
    result[7, 0] = numpy.sum(
        x810
        * (
            x1862 * x77
            + x3 * (x1859 * x3 + 4.0 * x1861)
            + x83 * (x1860 * x3 + 3.0 * x1863)
        )
    )
    result[7, 1] = numpy.sum(
        x1390
        * (x1866 * x77 + x3 * (x1864 * x3 + x1865 * x92) + x92 * (x1865 * x3 + x1867))
    )
    result[7, 2] = numpy.sum(
        x1390
        * (
            x1871 * x77
            + x3 * (x1868 * x3 + 3.0 * x1870)
            + x92 * (x1869 * x3 + 2.0 * x1872)
        )
    )
    result[7, 3] = numpy.sum(
        x1877
        * (x104 * (x1874 * x3 + x1876) + x1875 * x77 + x3 * (x104 * x1874 + x1873 * x3))
    )
    result[7, 4] = numpy.sum(
        x1883 * (x104 * (x1879 * x3 + x1882) + x1881 * x77 + x3 * (x1878 * x3 + x1880))
    )
    result[7, 5] = numpy.sum(
        x1877
        * (x104 * (x1885 * x3 + x1888) + x1887 * x77 + x3 * (x1884 * x3 + 2.0 * x1886))
    )
    result[7, 6] = numpy.sum(
        x1390 * (x1890 * x3 + x1892 * x77 + x3 * (x1890 + x1891 * x3))
    )
    result[7, 7] = numpy.sum(
        x1883 * (x1894 * x3 + x1896 * x77 + x3 * (x1894 + x1895 * x3))
    )
    result[7, 8] = numpy.sum(
        x1883 * (x1898 * x3 + x1900 * x77 + x3 * (x1898 + x1899 * x3))
    )
    result[7, 9] = numpy.sum(
        x1390 * (x1902 * x3 + x1904 * x77 + x3 * (x1902 + x1903 * x3))
    )
    result[7, 10] = numpy.sum(x810 * (x1905 * x238 + x1906 * x77))
    result[7, 11] = numpy.sum(x1390 * (x1907 * x238 + x1908 * x77))
    result[7, 12] = numpy.sum(x1877 * (x1909 * x238 + x1910 * x77))
    result[7, 13] = numpy.sum(x1390 * (x1911 * x238 + x1912 * x77))
    result[7, 14] = numpy.sum(x810 * (x1913 * x238 + x1914 * x77))
    result[8, 0] = numpy.sum(
        0.5
        * x810
        * (x1915 + 2.0 * x3 * x690 * (x1605 + x1608) + 2.0 * x690 * x83 * (x1612 + x1615))
    )
    result[8, 1] = numpy.sum(
        0.5
        * x1390
        * (
            x1918
            + 2.0 * x3 * (x1916 * x3 + x1917 * x92)
            + 2.0 * x92 * (x1917 * x3 + x1919)
        )
    )
    result[8, 2] = numpy.sum(
        0.5
        * x1390
        * (x1920 + 2.0 * x3 * x690 * (x1646 + x1649) + 2.0 * x690 * x92 * (x1653 + x1656))
    )
    result[8, 3] = numpy.sum(
        0.5
        * x1877
        * (
            2.0 * x104 * (x1922 * x3 + x1924)
            + x1923
            + 2.0 * x3 * (x104 * x1922 + x1921 * x3)
        )
    )
    result[8, 4] = numpy.sum(
        0.5
        * x1883
        * (2.0 * x104 * (x1926 * x3 + x1929) + x1928 + 2.0 * x3 * (x1925 * x3 + x1927))
    )
    result[8, 5] = numpy.sum(
        0.5
        * x1877
        * (
            2.0 * x104 * (x1710 * x690 + x1931)
            + x1930
            + 2.0 * x3 * x690 * (x1701 + x1704)
        )
    )
    result[8, 6] = numpy.sum(
        0.5 * x1390 * (4.0 * x1933 * x3 + 2.0 * x1934 * x3**2 + x1935)
    )
    result[8, 7] = numpy.sum(
        0.5 * x1883 * (4.0 * x1938 * x3 + 2.0 * x1939 * x3**2 + x1940)
    )
    result[8, 8] = numpy.sum(
        0.5 * x1883 * (4.0 * x1942 * x3 + 2.0 * x1943 * x3**2 + x1944)
    )
    result[8, 9] = numpy.sum(
        0.5 * x1390 * (2.0 * x1769 * x690 + x1946 + 2.0 * x3 * (x1765 * x690 + x1945))
    )
    result[8, 10] = numpy.sum(0.5 * x810 * (2.0 * x1947 * x238 + x1948))
    result[8, 11] = numpy.sum(0.5 * x1390 * (2.0 * x1949 * x238 + x1950))
    result[8, 12] = numpy.sum(0.5 * x1877 * (2.0 * x1951 * x238 + x1952))
    result[8, 13] = numpy.sum(0.5 * x1390 * (2.0 * x1953 * x238 + x1954))
    result[8, 14] = numpy.sum(0.5 * x810 * (2.0 * x1798 * x690 + x1955))
    result[9, 0] = numpy.sum(
        x1806 * (x1962 + x3 * (x1957 + x1960) + x83 * (x1958 * x3 + 3.0 * x1963))
    )
    result[9, 1] = numpy.sum(
        x1811 * (x1968 + x3 * (x1964 * x3 + 3.0 * x1966) + x92 * (x1965 * x3 + x1970))
    )
    result[9, 2] = numpy.sum(
        x1811 * (x1977 + x3 * (x1972 + x1975) + x92 * (x1973 * x3 + x1979))
    )
    result[9, 3] = numpy.sum(
        x1821 * (x104 * (x1981 * x3 + x1985) + x1984 + x3 * (x1980 * x3 + 2.0 * x1982))
    )
    result[9, 4] = numpy.sum(
        x1827 * (x104 * (x1987 * x3 + x1991) + x1990 + x3 * (x1986 * x3 + x1988))
    )
    result[9, 5] = numpy.sum(
        x1821 * (x104 * (x1994 * x3 + x1999) + x1998 + x3 * (x1993 + x1996))
    )
    result[9, 6] = numpy.sum(x1811 * (x2001 * x3 + x2004 + x3 * (x2001 + x2002 * x3)))
    result[9, 7] = numpy.sum(x1827 * (x2006 * x3 + x2009 + x3 * (x2006 + x2007 * x3)))
    result[9, 8] = numpy.sum(x1827 * (x2011 * x3 + x2014 + x3 * (x2011 + x2012 * x3)))
    result[9, 9] = numpy.sum(x1811 * (x2016 * x3 + x2020 + x3 * (x2016 + x2018)))
    result[9, 10] = numpy.sum(x1806 * (x2021 * x238 + x2023))
    result[9, 11] = numpy.sum(x1811 * (x2024 * x238 + x2026))
    result[9, 12] = numpy.sum(x1821 * (x2027 * x238 + x2029))
    result[9, 13] = numpy.sum(x1811 * (x2030 * x238 + x2032))
    result[9, 14] = numpy.sum(x1806 * (x2033 * x238 + x2035))
    result[10, 0] = numpy.sum(x718 * (x2036 * x3 + 4.0 * x2037))
    result[10, 1] = numpy.sum(x757 * (x2038 * x3 + x2039 * x92))
    result[10, 2] = numpy.sum(x757 * (x2040 * x3 + 3.0 * x2041))
    result[10, 3] = numpy.sum(x810 * (x104 * x2043 + x2042 * x3))
    result[10, 4] = numpy.sum(x840 * (x2044 * x3 + x2045))
    result[10, 5] = numpy.sum(x810 * (x2046 * x3 + 2.0 * x2047))
    result[10, 6] = numpy.sum(x757 * (x2048 + x2049 * x3))
    result[10, 7] = numpy.sum(x840 * (x2050 + x2051 * x3))
    result[10, 8] = numpy.sum(x840 * (x2052 + x2053 * x3))
    result[10, 9] = numpy.sum(x757 * (x2054 + x2055 * x3))
    result[10, 10] = numpy.sum(x2056 * x2057)
    result[10, 11] = numpy.sum(x2058 * x2059)
    result[10, 12] = numpy.sum(x2060 * x2061)
    result[10, 13] = numpy.sum(x2059 * x2062)
    result[10, 14] = numpy.sum(x2057 * x2063)
    result[11, 0] = numpy.sum(x757 * (x2064 * x3 + 4.0 * x2065))
    result[11, 1] = numpy.sum(x1498 * (x2066 * x3 + x2067 * x92))
    result[11, 2] = numpy.sum(x1498 * (x2068 * x3 + 3.0 * x2069))
    result[11, 3] = numpy.sum(x1390 * (x104 * x2071 + x2070 * x3))
    result[11, 4] = numpy.sum(x1529 * (x2072 * x3 + x2073))
    result[11, 5] = numpy.sum(x1390 * (x2074 * x3 + 2.0 * x2075))
    result[11, 6] = numpy.sum(x1498 * (x2076 + x2077 * x3))
    result[11, 7] = numpy.sum(x1529 * (x2078 + x2079 * x3))
    result[11, 8] = numpy.sum(x1529 * (x2080 + x2081 * x3))
    result[11, 9] = numpy.sum(x1498 * (x2082 + x2083 * x3))
    result[11, 10] = numpy.sum(x2059 * x2084)
    result[11, 11] = numpy.sum(x2085 * x2086)
    result[11, 12] = numpy.sum(x2087 * x2088)
    result[11, 13] = numpy.sum(x2086 * x2089)
    result[11, 14] = numpy.sum(x2059 * x2090)
    result[12, 0] = numpy.sum(x810 * (x2091 * x3 + 4.0 * x2092))
    result[12, 1] = numpy.sum(x1390 * (x2093 * x3 + x2094 * x92))
    result[12, 2] = numpy.sum(x1390 * (x2095 * x3 + 3.0 * x2096))
    result[12, 3] = numpy.sum(x1877 * (x104 * x2098 + x2097 * x3))
    result[12, 4] = numpy.sum(x1883 * (x2099 * x3 + x2100))
    result[12, 5] = numpy.sum(x1877 * (x2101 * x3 + 2.0 * x2102))
    result[12, 6] = numpy.sum(x1390 * (x2103 + x2104 * x3))
    result[12, 7] = numpy.sum(x1883 * (x2105 + x2106 * x3))
    result[12, 8] = numpy.sum(x1883 * (x2107 + x2108 * x3))
    result[12, 9] = numpy.sum(x1390 * (x2109 + x2110 * x3))
    result[12, 10] = numpy.sum(x2061 * x2111)
    result[12, 11] = numpy.sum(x2088 * x2112)
    result[12, 12] = numpy.sum(x1877 * x2113 * x3)
    result[12, 13] = numpy.sum(x2088 * x2114)
    result[12, 14] = numpy.sum(x2061 * x2115)
    result[13, 0] = numpy.sum(x690 * x757 * (x1957 + x1960))
    result[13, 1] = numpy.sum(x1498 * (x2116 * x3 + x2117 * x92))
    result[13, 2] = numpy.sum(x1498 * x690 * (x1972 + x1975))
    result[13, 3] = numpy.sum(x1390 * (x104 * x2119 + x2118 * x3))
    result[13, 4] = numpy.sum(x1529 * (x2120 * x3 + x2121))
    result[13, 5] = numpy.sum(x1390 * x690 * (x1993 + x1996))
    result[13, 6] = numpy.sum(x1498 * (x2122 + x2123 * x3))
    result[13, 7] = numpy.sum(x1529 * (x2125 + x2126 * x3))
    result[13, 8] = numpy.sum(x1529 * (x2127 + x2128 * x3))
    result[13, 9] = numpy.sum(x1498 * (x2018 * x690 + x2129))
    result[13, 10] = numpy.sum(x2059 * x2130)
    result[13, 11] = numpy.sum(x2086 * x2131)
    result[13, 12] = numpy.sum(x2088 * x2132)
    result[13, 13] = numpy.sum(x2086 * x2133)
    result[13, 14] = numpy.sum(x2033 * x2059 * x690)
    result[14, 0] = numpy.sum(x718 * (x2134 * x3 + 4.0 * x2135))
    result[14, 1] = numpy.sum(x757 * (x2136 * x3 + 3.0 * x2137))
    result[14, 2] = numpy.sum(x757 * (x2138 * x3 + 3.0 * x2139))
    result[14, 3] = numpy.sum(x810 * (x2140 * x3 + 2.0 * x2141))
    result[14, 4] = numpy.sum(x840 * (x2142 * x3 + x2143))
    result[14, 5] = numpy.sum(x810 * (x2144 * x3 + 2.0 * x2145))
    result[14, 6] = numpy.sum(x757 * (x2146 + x2147 * x3))
    result[14, 7] = numpy.sum(x840 * (x2148 + x2149 * x3))
    result[14, 8] = numpy.sum(x840 * (x2150 + x2151 * x3))
    result[14, 9] = numpy.sum(x757 * (x2152 + x2153 * x3))
    result[14, 10] = numpy.sum(x2057 * x2154)
    result[14, 11] = numpy.sum(x2059 * x2155)
    result[14, 12] = numpy.sum(x2061 * x2156)
    result[14, 13] = numpy.sum(x2059 * x2157)
    result[14, 14] = numpy.sum(x2057 * x2158)
    result[15, 0] = numpy.sum(x149 * (x1804 * x23 + x2036 * x690))
    result[15, 1] = numpy.sum(x239 * (x1809 * x23 + x2037 + x2038 * x690))
    result[15, 2] = numpy.sum(x239 * (x1815 * x23 + x2040 * x690))
    result[15, 3] = numpy.sum(x378 * (x104 * x2039 + x1819 * x23 + x2042 * x690))
    result[15, 4] = numpy.sum(x426 * (x1825 * x23 + x2041 + x2044 * x690))
    result[15, 5] = numpy.sum(x378 * (x1831 * x23 + x2046 * x690))
    result[15, 6] = numpy.sum(x239 * (x1836 * x23 + x2043 * x92 + x2049 * x690))
    result[15, 7] = numpy.sum(x426 * (x1840 * x23 + x2045 + x2051 * x690))
    result[15, 8] = numpy.sum(x426 * (x1844 * x23 + x2047 + x2053 * x690))
    result[15, 9] = numpy.sum(x239 * (x1848 * x23 + x2055 * x690))
    result[15, 10] = numpy.sum(x149 * (x1850 * x23 + 4.0 * x2048 + x2056 * x690))
    result[15, 11] = numpy.sum(x239 * (x1852 * x23 + 3.0 * x2050 + x2058 * x690))
    result[15, 12] = numpy.sum(x378 * (x1854 * x23 + 2.0 * x2052 + x2060 * x690))
    result[15, 13] = numpy.sum(x239 * (x1856 * x23 + x2054 + x2062 * x690))
    result[15, 14] = numpy.sum(x149 * (x1858 * x23 + x2063 * x690))
    result[16, 0] = numpy.sum(x718 * (x129 * x1862 + x2064 * x690))
    result[16, 1] = numpy.sum(x757 * (x129 * x1866 + x2065 + x2066 * x690))
    result[16, 2] = numpy.sum(x757 * (x129 * x1871 + x2068 * x690))
    result[16, 3] = numpy.sum(x810 * (x104 * x2067 + x129 * x1875 + x2070 * x690))
    result[16, 4] = numpy.sum(x840 * (x129 * x1881 + x2069 + x2072 * x690))
    result[16, 5] = numpy.sum(x810 * (x129 * x1887 + x2074 * x690))
    result[16, 6] = numpy.sum(x757 * (x129 * x1892 + x2071 * x92 + x2077 * x690))
    result[16, 7] = numpy.sum(x840 * (x129 * x1896 + x2073 + x2079 * x690))
    result[16, 8] = numpy.sum(x840 * (x129 * x1900 + x2075 + x2081 * x690))
    result[16, 9] = numpy.sum(x757 * (x129 * x1904 + x2083 * x690))
    result[16, 10] = numpy.sum(x718 * (x129 * x1906 + 4.0 * x2076 + x2084 * x690))
    result[16, 11] = numpy.sum(x757 * (x129 * x1908 + 3.0 * x2078 + x2085 * x690))
    result[16, 12] = numpy.sum(x810 * (x129 * x1910 + 2.0 * x2080 + x2087 * x690))
    result[16, 13] = numpy.sum(x757 * (x129 * x1912 + x2082 + x2089 * x690))
    result[16, 14] = numpy.sum(x718 * (x129 * x1914 + x2090 * x690))
    result[17, 0] = numpy.sum(x1330 * (x1915 + x2091 * x690))
    result[17, 1] = numpy.sum(x1348 * (x1918 + x2092 + x2093 * x690))
    result[17, 2] = numpy.sum(x1348 * (x1920 + x2095 * x690))
    result[17, 3] = numpy.sum(x1375 * (x104 * x2094 + x1923 + x2097 * x690))
    result[17, 4] = numpy.sum(x1390 * (x1928 + x2096 + x2099 * x690))
    result[17, 5] = numpy.sum(x1375 * (x1930 + x2101 * x690))
    result[17, 6] = numpy.sum(x1348 * (x1935 + x2098 * x92 + x2104 * x690))
    result[17, 7] = numpy.sum(x1390 * (x1940 + x2100 + x2106 * x690))
    result[17, 8] = numpy.sum(x1390 * (x1944 + x2102 + x2108 * x690))
    result[17, 9] = numpy.sum(x1348 * (x1946 + x2110 * x690))
    result[17, 10] = numpy.sum(x1330 * (x1948 + 4.0 * x2103 + x2111 * x690))
    result[17, 11] = numpy.sum(x1348 * (x1950 + 3.0 * x2105 + x2112 * x690))
    result[17, 12] = numpy.sum(x1375 * (x1952 + 2.0 * x2107 + x2113 * x690))
    result[17, 13] = numpy.sum(x1348 * (x1954 + x2109 + x2114 * x690))
    result[17, 14] = numpy.sum(x1330 * (x1955 + x2115 * x690))
    result[18, 0] = numpy.sum(x1806 * (x1309 * x1956 + x1962))
    result[18, 1] = numpy.sum(x1811 * (x1959 * x690 + x1968 + x2116 * x690))
    result[18, 2] = numpy.sum(x1811 * (x1309 * x1971 + x1977))
    result[18, 3] = numpy.sum(x1821 * (x104 * x2117 + x1984 + x2118 * x690))
    result[18, 4] = numpy.sum(x1827 * (x1974 * x690 + x1990 + x2120 * x690))
    result[18, 5] = numpy.sum(x1821 * (x1309 * x1992 + x1998))
    result[18, 6] = numpy.sum(x1811 * (x2004 + x2119 * x92 + x2123 * x690))
    result[18, 7] = numpy.sum(x1827 * (x2009 + x2121 + x2126 * x690))
    result[18, 8] = numpy.sum(x1827 * (x1995 * x690 + x2014 + x2128 * x690))
    result[18, 9] = numpy.sum(x1811 * (x1309 * x2017 + x2020))
    result[18, 10] = numpy.sum(x1806 * (x2023 + 4.0 * x2122 + x2130 * x690))
    result[18, 11] = numpy.sum(x1811 * (x2026 + 3.0 * x2125 + x2131 * x690))
    result[18, 12] = numpy.sum(x1821 * (x2029 + 2.0 * x2127 + x2132 * x690))
    result[18, 13] = numpy.sum(x1811 * (x2032 + x2129 + x2133 * x690))
    result[18, 14] = numpy.sum(x1806 * (x1309 * x2033 + x2035))
    result[19, 0] = numpy.sum(x2134 * x2159)
    result[19, 1] = numpy.sum(x757 * (x2135 + x2136 * x690))
    result[19, 2] = numpy.sum(x2138 * x2160)
    result[19, 3] = numpy.sum(x810 * (2.0 * x2137 + x2140 * x690))
    result[19, 4] = numpy.sum(x840 * (x2139 + x2142 * x690))
    result[19, 5] = numpy.sum(x2144 * x690 * x810)
    result[19, 6] = numpy.sum(x757 * (3.0 * x2141 + x2147 * x690))
    result[19, 7] = numpy.sum(x840 * (x2143 + x2149 * x690))
    result[19, 8] = numpy.sum(x840 * (x2145 + x2151 * x690))
    result[19, 9] = numpy.sum(x2153 * x2160)
    result[19, 10] = numpy.sum(x718 * (4.0 * x2146 + x2154 * x690))
    result[19, 11] = numpy.sum(x757 * (3.0 * x2148 + x2155 * x690))
    result[19, 12] = numpy.sum(x810 * (2.0 * x2150 + x2156 * x690))
    result[19, 13] = numpy.sum(x757 * (x2152 + x2157 * x690))
    result[19, 14] = numpy.sum(x2158 * x2159)
    result[20, 0] = numpy.sum(x149 * (x1961 * x23 + x2134 * x979))
    result[20, 1] = numpy.sum(x239 * (x1967 * x23 + x2136 * x979))
    result[20, 2] = numpy.sum(x239 * (x1976 * x23 + x2135 + x2138 * x979))
    result[20, 3] = numpy.sum(x378 * (x1983 * x23 + x2140 * x979))
    result[20, 4] = numpy.sum(x426 * (x1989 * x23 + x2137 + x2142 * x979))
    result[20, 5] = numpy.sum(x378 * (x1997 * x23 + 2.0 * x2139 + x2144 * x979))
    result[20, 6] = numpy.sum(x239 * (x2003 * x23 + x2147 * x979))
    result[20, 7] = numpy.sum(x426 * (x2008 * x23 + x2141 + x2149 * x979))
    result[20, 8] = numpy.sum(x426 * (x2013 * x23 + x2143 + x2151 * x979))
    result[20, 9] = numpy.sum(x239 * (x2019 * x23 + 3.0 * x2145 + x2153 * x979))
    result[20, 10] = numpy.sum(x149 * (x2022 * x23 + x2154 * x979))
    result[20, 11] = numpy.sum(x239 * (x2025 * x23 + x2146 + x2155 * x979))
    result[20, 12] = numpy.sum(x378 * (x2028 * x23 + 2.0 * x2148 + x2156 * x979))
    result[20, 13] = numpy.sum(x239 * (x2031 * x23 + 3.0 * x2150 + x2157 * x979))
    result[20, 14] = numpy.sum(x149 * (x2034 * x23 + 4.0 * x2152 + x2158 * x979))
    return result


def _2center2el3d_55(ax, da, A, bx, db, B):
    """Cartesian (h|h) two-center two-electron repulsion integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((21, 21), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = ax ** (-1.0)
    x5 = -x2 - B[0]
    x6 = bx ** (-1.0)
    x7 = ax * x1
    x8 = bx * x7 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
    x9 = boys(8, x8)
    x10 = x0 ** (-1.5)
    x11 = 17.4934183276249
    x12 = 2.0 * x6
    x13 = x11 * x12
    x14 = x10 * x13
    x15 = x14 * x9
    x16 = x15 * x5
    x17 = x0 ** (-0.5)
    x18 = boys(7, x8)
    x19 = 0.5 * x6
    x20 = x19 * (2.0 * x11 * x17 * x18 * x4 * x6 - x15)
    x21 = boys(9, x8)
    x22 = x5**2
    x23 = x17 * x4
    x24 = x13 * x23
    x25 = x22 * x24
    x26 = x21 * x25
    x27 = x20 + x26
    x28 = x27 * x5 + x6 * (2.0 * x11 * x17 * x18 * x4 * x5 * x6 - x16)
    x29 = x14 * x18
    x30 = boys(6, x8)
    x31 = x19 * (2.0 * x11 * x17 * x30 * x4 * x6 - x29)
    x32 = x24 * x9
    x33 = x22 * x32
    x34 = x31 + x33
    x35 = x14 * x30
    x36 = boys(5, x8)
    x37 = x19 * (2.0 * x11 * x17 * x36 * x4 * x6 - x35)
    x38 = x18 * x24
    x39 = x22 * x38
    x40 = x37 + x39
    x41 = 1.5 * x6
    x42 = x28 * x5 - x41 * (x34 * x7 - x40)
    x43 = x29 * x5
    x44 = x34 * x5 + x6 * (2.0 * x11 * x17 * x30 * x4 * x5 * x6 - x43)
    x45 = x35 * x5
    x46 = x40 * x5 + x6 * (2.0 * x11 * x17 * x36 * x4 * x5 * x6 - x45)
    x47 = -x12 * (x44 * x7 - x46) + x42 * x5
    x48 = x3 * x47
    x49 = 0.5 / (ax + bx)
    x50 = x14 * x36
    x51 = boys(4, x8)
    x52 = x19 * (2.0 * x11 * x17 * x4 * x51 * x6 - x50)
    x53 = x24 * x30
    x54 = x22 * x53
    x55 = x52 + x54
    x56 = -x41 * (x40 * x7 - x55) + x44 * x5
    x57 = x49 * x56
    x58 = 5.0 * x57
    x59 = x48 + x58
    x60 = bx * x1
    x61 = x5 * x50
    x62 = x5 * x55 + x6 * (2.0 * x11 * x17 * x4 * x5 * x51 * x6 - x61)
    x63 = -x12 * (x46 * x7 - x62) + x5 * x56
    x64 = x3 * x63
    x65 = x14 * x51
    x66 = boys(3, x8)
    x67 = x19 * (2.0 * x11 * x17 * x4 * x6 * x66 - x65)
    x68 = x24 * x36
    x69 = x22 * x68
    x70 = x67 + x69
    x71 = -x41 * (x55 * x7 - x70) + x46 * x5
    x72 = x49 * x71
    x73 = 5.0 * x72
    x74 = x64 + x73
    x75 = x14 * x21
    x76 = x5 * x75
    x77 = x19 * (2.0 * x11 * x17 * x4 * x6 * x9 - x75)
    x78 = boys(10, x8)
    x79 = x25 * x78
    x80 = -x12 * (x28 * x7 - x44) - x5 * (
        x41 * (x27 * x7 - x34)
        - x5 * (x5 * (x77 + x79) + x6 * (2.0 * x11 * x17 * x4 * x5 * x6 * x9 - x76))
    )
    x81 = x3 * x80
    x82 = x42 * x49
    x83 = 5.0 * x82
    x84 = x47 * x60
    x85 = 0.5 * x4
    x86 = x85 * (x63 - x84)
    x87 = x3 * x42
    x88 = x44 * x49
    x89 = 4.0 * x88
    x90 = x87 + x89
    x91 = 5.0 * x49
    x92 = x56 * x60
    x93 = x85 * (x71 - x92)
    x94 = x3 * x44
    x95 = x40 * x49
    x96 = 3.0 * x95
    x97 = x94 + x96
    x98 = 4.0 * x49
    x99 = x3 * x90 + x93 + x97 * x98
    x100 = x3 * x56
    x101 = x46 * x49
    x102 = 4.0 * x101
    x103 = x100 + x102
    x104 = x3 * x71
    x105 = x49 * x62
    x106 = 4.0 * x105
    x107 = x104 + x106
    x108 = x46 * x60
    x109 = x85 * (-x108 + x62)
    x110 = x3 * x40
    x111 = x11 * x30
    x112 = x111 * x5
    x113 = x23 * x6 * x98
    x114 = x112 * x113
    x115 = x110 + x114
    x116 = 3.0 * x49
    x117 = x109 + x115 * x116 + x3 * x97
    x118 = x117 * x98 + x3 * x99 - x4 * (x103 * x60 - x107)
    x119 = x5 * x65
    x120 = x5 * x70 + x6 * (2.0 * x11 * x17 * x4 * x5 * x6 * x66 - x119)
    x121 = x12 * (x120 - x62 * x7) + x5 * x71
    x122 = x85 * (x121 - x60 * x63)
    x123 = x103 * x91 + x122 + x3 * x59
    x124 = x14 * x66
    x125 = boys(2, x8)
    x126 = x19 * (2.0 * x11 * x125 * x17 * x4 * x6 - x124)
    x127 = x24 * x51
    x128 = x127 * x22
    x129 = x126 + x128
    x130 = x41 * (x129 - x7 * x70) + x5 * x62
    x131 = (
        x12
        * (-x120 * x7 + x129 * x5 + x5 * x6 * (2.0 * x11 * x125 * x17 * x4 * x6 - x124))
        + x130 * x5
    )
    x132 = x85 * (-x121 * x60 + x131)
    x133 = x107 * x91 + x132 + x3 * x74
    x134 = 1.5 * x4
    x135 = x19 * (2.0 * x11 * x17 * x4 * x6 * boys(1, x8) - x125 * x14)
    x136 = x24 * x66
    x137 = x120 * x5 + x41 * (-x129 * x7 + x135 + x136 * x22)
    x138 = x137 * x49
    x139 = x130 * x49
    x140 = x121 * x3 + 5.0 * x139
    x141 = x85 * (-x130 * x60 + x137)
    x142 = x49 * x70
    x143 = 3.0 * x142 + x3 * x62
    x144 = x107 * x3 + x141 + x143 * x98
    x145 = x85 * (x130 - x60 * x71)
    x146 = x49 * x55
    x147 = 3.0 * x146 + x3 * x46
    x148 = x103 * x3 + x145 + x147 * x98
    x149 = 2.0 * x4
    x150 = x85 * (-x55 * x60 + x70)
    x151 = x49 * x68
    x152 = x3 * x5
    x153 = 2.0 * x49
    x154 = 0.179587122125167 * da * db * numpy.sqrt(ax**6.5) * numpy.sqrt(bx**6.5)
    x155 = 3.06487764628582 * x154
    x156 = -x1 * (ax * A[1] + bx * B[1])
    x157 = -x156 - B[1]
    x158 = x157 * x5
    x159 = x15 * x157
    x160 = x6 * (2.0 * x11 * x157 * x17 * x18 * x4 * x6 - x159)
    x161 = x157 * x26
    x162 = 0.5 * x160 + x161
    x163 = x162 * x5 + x6 * (2.0 * x11 * x157 * x17 * x18 * x4 * x5 * x6 - x15 * x158)
    x164 = x157 * x29
    x165 = x6 * (2.0 * x11 * x157 * x17 * x30 * x4 * x6 - x164)
    x166 = x157 * x33
    x167 = 0.5 * x165 + x166
    x168 = x157 * x35
    x169 = x6 * (2.0 * x11 * x157 * x17 * x36 * x4 * x6 - x168)
    x170 = x157 * x39
    x171 = 0.5 * x169 + x170
    x172 = x163 * x5 - x41 * (x167 * x7 - x171)
    x173 = x172 * x3
    x174 = x167 * x5 + x6 * (2.0 * x11 * x157 * x17 * x30 * x4 * x5 * x6 - x158 * x29)
    x175 = x174 * x49
    x176 = 4.0 * x175
    x177 = x173 + x176
    x178 = x157 * x50
    x179 = x6 * (2.0 * x11 * x157 * x17 * x4 * x51 * x6 - x178)
    x180 = x157 * x54
    x181 = 0.5 * x179 + x180
    x182 = x174 * x5 - x41 * (x171 * x7 - x181)
    x183 = x182 * x3
    x184 = -2.0 * x11 * x157 * x17 * x36 * x4 * x5 * x6
    x185 = x171 * x5 - x6 * (x158 * x35 + x184)
    x186 = x185 * x49
    x187 = x183 + 4.0 * x186
    x188 = x157 * x75
    x189 = x6 * (2.0 * x11 * x157 * x17 * x4 * x6 * x9 - x188)
    x190 = x157 * x79
    x191 = -x41 * (x162 * x7 - x167) + 0.5 * x5 * (
        x5 * (x189 + 2.0 * x190)
        + 2.0 * x6 * (2.0 * x11 * x157 * x17 * x4 * x5 * x6 * x9 - x158 * x75)
    )
    x192 = x191 * x3
    x193 = x163 * x49
    x194 = 4.0 * x193
    x195 = x172 * x60
    x196 = x85 * (x182 - x195)
    x197 = x163 * x3
    x198 = x167 * x49
    x199 = 3.0 * x198
    x200 = x197 + x199
    x201 = x174 * x60
    x202 = x85 * (x185 - x201)
    x203 = x167 * x3
    x204 = x11 * x158
    x205 = x18 * x204
    x206 = x113 * x205
    x207 = x203 + x206
    x208 = x116 * x207 + x200 * x3 + x202
    x209 = x174 * x3
    x210 = x171 * x49
    x211 = 3.0 * x210
    x212 = x209 + x211
    x213 = x185 * x3
    x214 = x181 * x49
    x215 = x213 + 3.0 * x214
    x216 = x171 * x60
    x217 = x85 * (x181 - x216)
    x218 = x157 * x53
    x219 = x218 * x49
    x220 = x152 * x157
    x221 = x220 * x38
    x222 = x219 + x221
    x223 = x153 * x222 + x207 * x3 + x217
    x224 = x116 * x223 + x208 * x3 - x4 * (x212 * x60 - x215)
    x225 = x157 * x65
    x226 = x6 * (2.0 * x11 * x157 * x17 * x4 * x6 * x66 - x225)
    x227 = x157 * x69
    x228 = 0.5 * x226 + x227
    x229 = x185 * x5 - x41 * (x181 * x7 - x228)
    x230 = x85 * (-x182 * x60 + x229)
    x231 = x177 * x3 + x212 * x98 + x230
    x232 = x181 * x5 + x6 * (2.0 * x11 * x157 * x17 * x4 * x5 * x51 * x6 - x158 * x50)
    x233 = x157 * x6 * (2.0 * x11 * x125 * x17 * x4 * x6 - x124)
    x234 = x128 * x157 * x41 - x228 * x41 * x7 + x232 * x5 + 0.5 * x233 * x41
    x235 = x85 * (-x229 * x60 + x234)
    x236 = x187 * x3 + x215 * x98 + x235
    x237 = x228 * x5 + x6 * (2.0 * x11 * x157 * x17 * x4 * x5 * x6 * x66 - x158 * x65)
    x238 = x237 * x49
    x239 = x232 * x49
    x240 = x229 * x3 + 4.0 * x239
    x241 = x85 * (-x232 * x60 + x237)
    x242 = x113 * x204 * x36
    x243 = x181 * x3 + x242
    x244 = x116 * x243 + x215 * x3 + x241
    x245 = x85 * (-x185 * x60 + x232)
    x246 = x111 * x158
    x247 = x113 * x246
    x248 = x171 * x3 + x247
    x249 = x116 * x248 + x212 * x3 + x245
    x250 = x10 * x149
    x251 = -x85 * (x184 + x246 * x250)
    x252 = 9.19463293885746 * x154
    x253 = -x1 * (ax * A[2] + bx * B[2])
    x254 = -x253 - B[2]
    x255 = x15 * x254
    x256 = x6 * (2.0 * x11 * x17 * x18 * x254 * x4 * x6 - x255)
    x257 = 0.5 * x256
    x258 = x254 * x26 + x257
    x259 = x254 * x6 * (2.0 * x11 * x17 * x18 * x4 * x5 * x6 - x16) + x258 * x5
    x260 = x254 * x29
    x261 = x6 * (2.0 * x11 * x17 * x254 * x30 * x4 * x6 - x260)
    x262 = 0.5 * x261
    x263 = x254 * x33 + x262
    x264 = x254 * x35
    x265 = x6 * (2.0 * x11 * x17 * x254 * x36 * x4 * x6 - x264)
    x266 = 0.5 * x265
    x267 = x254 * x39 + x266
    x268 = x259 * x5 - x41 * (x263 * x7 - x267)
    x269 = x268 * x3
    x270 = x254 * x6 * (2.0 * x11 * x17 * x30 * x4 * x5 * x6 - x43) + x263 * x5
    x271 = x270 * x49
    x272 = 4.0 * x271
    x273 = x269 + x272
    x274 = x254 * x50
    x275 = x6 * (2.0 * x11 * x17 * x254 * x4 * x51 * x6 - x274)
    x276 = 0.5 * x275
    x277 = x254 * x54 + x276
    x278 = x270 * x5 - x41 * (x267 * x7 - x277)
    x279 = x278 * x3
    x280 = -2.0 * x11 * x17 * x254 * x36 * x4 * x5 * x6
    x281 = x267 * x5 - x6 * (x254 * x45 + x280)
    x282 = x281 * x49
    x283 = x279 + 4.0 * x282
    x284 = x254 * x75
    x285 = x6 * (2.0 * x11 * x17 * x254 * x4 * x6 * x9 - x284)
    x286 = 0.5 * x285
    x287 = -x41 * (x258 * x7 - x263) + x5 * (
        x254 * x6 * (2.0 * x11 * x17 * x4 * x5 * x6 * x9 - x76) + x5 * (x254 * x79 + x286)
    )
    x288 = x287 * x3
    x289 = x259 * x49
    x290 = 4.0 * x289
    x291 = x268 * x60
    x292 = x85 * (x278 - x291)
    x293 = x259 * x3
    x294 = x263 * x49
    x295 = 3.0 * x294
    x296 = x293 + x295
    x297 = x270 * x60
    x298 = x85 * (x281 - x297)
    x299 = x263 * x3
    x300 = x113 * x254
    x301 = x11 * x300
    x302 = x301 * x5
    x303 = x18 * x302
    x304 = x299 + x303
    x305 = x116 * x304 + x296 * x3 + x298
    x306 = x270 * x3
    x307 = x267 * x49
    x308 = 3.0 * x307
    x309 = x306 + x308
    x310 = x281 * x3
    x311 = x277 * x49
    x312 = x310 + 3.0 * x311
    x313 = x267 * x60
    x314 = x85 * (x277 - x313)
    x315 = x254 * x53
    x316 = x315 * x49
    x317 = x254 * x38
    x318 = x152 * x317
    x319 = x316 + x318
    x320 = x153 * x319 + x3 * x304 + x314
    x321 = x116 * x320 + x3 * x305 - x4 * (x309 * x60 - x312)
    x322 = x254 * x6 * (2.0 * x11 * x17 * x4 * x6 * x66 - x65)
    x323 = 0.5 * x322
    x324 = x254 * x69 + x323
    x325 = x281 * x5 - x41 * (x277 * x7 - x324)
    x326 = x85 * (-x278 * x60 + x325)
    x327 = x273 * x3 + x309 * x98 + x326
    x328 = x254 * x6 * (2.0 * x11 * x17 * x4 * x5 * x51 * x6 - x61) + x277 * x5
    x329 = x254 * x6 * (2.0 * x11 * x125 * x17 * x4 * x6 - x124)
    x330 = 0.5 * x329
    x331 = x328 * x5 + x41 * (x128 * x254 - x324 * x7 + x330)
    x332 = x85 * (-x325 * x60 + x331)
    x333 = x283 * x3 + x312 * x98 + x332
    x334 = x254 * x6 * (2.0 * x11 * x17 * x4 * x5 * x6 * x66 - x119) + x324 * x5
    x335 = x334 * x49
    x336 = x328 * x49
    x337 = x3 * x325 + 4.0 * x336
    x338 = x85 * (-x328 * x60 + x334)
    x339 = x302 * x36
    x340 = x277 * x3 + x339
    x341 = x116 * x340 + x3 * x312 + x338
    x342 = x85 * (-x281 * x60 + x328)
    x343 = x114 * x254
    x344 = x267 * x3 + x343
    x345 = x116 * x344 + x3 * x309 + x342
    x346 = x250 * x254
    x347 = -x85 * (x112 * x346 + x280)
    x348 = x157**2
    x349 = x348 * x38 + x37
    x350 = x31 + x32 * x348
    x351 = x350 * x5
    x352 = x24 * x348
    x353 = x21 * x352
    x354 = x20 + x353
    x355 = x349 - x350 * x7
    x356 = x19 * x355 + x22 * x354
    x357 = x356 * x5 + x6 * (x349 * x5 - x351 * x7)
    x358 = x3 * x357
    x359 = x348 * x53 + x52
    x360 = -x349 * x7 + x359
    x361 = x19 * x360 + x22 * x350
    x362 = x361 * x49
    x363 = 3.0 * x362
    x364 = x358 + x363
    x365 = -x359 * x5
    x366 = x349 * x5
    x367 = x361 * x5 - x6 * (x365 + x366 * x7)
    x368 = x3 * x367
    x369 = x348 * x68 + x67
    x370 = -x359 * x7 + x369
    x371 = x19 * x370 + x22 * x349
    x372 = x371 * x49
    x373 = 3.0 * x372
    x374 = x368 + x373
    x375 = x354 * x7
    x376 = x352 * x78
    x377 = x376 + x77
    x378 = x350 - x375
    x379 = x5 * (x19 * x378 + x22 * x377 + x6 * (x350 - x375))
    x380 = x3 * x379
    x381 = x356 * x49
    x382 = 3.0 * x381
    x383 = x357 * x60
    x384 = x85 * (x367 - x383)
    x385 = x3 * x356
    x386 = x153 * x351
    x387 = x385 + x386
    x388 = x361 * x60
    x389 = x85 * (x371 - x388)
    x390 = x349 * x49
    x391 = x3 * x351
    x392 = x390 + x391
    x393 = x153 * x392 + x3 * x387 + x389
    x394 = x3 * x361
    x395 = x153 * x366
    x396 = x394 + x395
    x397 = x3 * x371
    x398 = x359 * x5
    x399 = x153 * x398 + x397
    x400 = x366 * x60
    x401 = -x85 * (x365 + x400)
    x402 = x3 * x390
    x403 = x3 * x392 + x401 + x402
    x404 = x153 * x403 + x3 * x393 - x4 * (x396 * x60 - x399)
    x405 = x371 * x5 + x6 * (x369 * x5 - x398 * x7)
    x406 = x85 * (-x367 * x60 + x405)
    x407 = x116 * x396 + x3 * x364 + x406
    x408 = x126 + x127 * x348
    x409 = -x369 * x7 + x408
    x410 = x19 * x409 + x22 * x359
    x411 = x5 * (x410 - x6 * (x369 * x7 - x408))
    x412 = x85 * (-x405 * x60 + x411)
    x413 = x116 * x399 + x3 * x374 + x412
    x414 = x135 + x136 * x348 - x408 * x7
    x415 = x19 * x414 + x22 * x369
    x416 = x415 * x49
    x417 = 3.0 * x416
    x418 = x410 * x49
    x419 = 3.0 * x418
    x420 = x3 * x405 + x419
    x421 = x85 * (-x410 * x60 + x415)
    x422 = x369 * x49
    x423 = x3 * x398 + x422
    x424 = x153 * x423 + x3 * x399 + x421
    x425 = x85 * (-x371 * x60 + x410)
    x426 = x359 * x49
    x427 = x3 * x366 + x426
    x428 = x153 * x427 + x3 * x396 + x425
    x429 = x3**2
    x430 = x85 * (-x359 * x60 + x369)
    x431 = 14.0450338098829 * x154
    x432 = x254 * x6 * (2.0 * x11 * x157 * x17 * x18 * x4 * x6 - x159)
    x433 = x161 * x254 + 0.5 * x432
    x434 = x433 * x5 + x6 * (
        2.0 * x11 * x157 * x17 * x18 * x254 * x4 * x5 * x6 - x158 * x255
    )
    x435 = x254 * x6 * (2.0 * x11 * x157 * x17 * x30 * x4 * x6 - x164)
    x436 = x166 * x254 + 0.5 * x435
    x437 = x116 * x436 + x3 * x434
    x438 = -2.0 * x11 * x157 * x17 * x254 * x30 * x4 * x5 * x6
    x439 = x436 * x5 - x6 * (x158 * x260 + x438)
    x440 = -2.0 * x11 * x157 * x17 * x254 * x36 * x4 * x6
    x441 = -x6 * (x168 * x254 + x440)
    x442 = x170 * x254 + 0.5 * x441
    x443 = x116 * x442 + x3 * x439
    x444 = x254 * x6 * (2.0 * x11 * x157 * x17 * x4 * x6 * x9 - x188)
    x445 = 0.5 * x5 * (2.0 * x190 * x254 + x444) + x6 * (
        2.0 * x11 * x157 * x17 * x254 * x4 * x5 * x6 * x9 - x158 * x284
    )
    x446 = x85 * (-x434 * x60 + x439)
    x447 = x204 * x300 * x9
    x448 = x3 * x433 + x447
    x449 = x85 * (-x436 * x60 + x442)
    x450 = x157 * x38
    x451 = x254 * x450
    x452 = x451 * x49
    x453 = x254 * x32
    x454 = x220 * x453 + x452
    x455 = x153 * x454 + x3 * x448 + x449
    x456 = x205 * x300
    x457 = x3 * x436 + x456
    x458 = x246 * x300
    x459 = x3 * x442 + x458
    x460 = -x85 * (x205 * x346 + x438)
    x461 = x3 * x452 + x3 * x454 + x460
    x462 = x153 * x461 + x3 * x455 - x4 * (x457 * x60 - x459)
    x463 = x442 * x5 + x6 * (
        2.0 * x11 * x157 * x17 * x254 * x36 * x4 * x5 * x6 - x158 * x264
    )
    x464 = x85 * (-x439 * x60 + x463)
    x465 = x116 * x457 + x3 * x437 + x464
    x466 = x254 * x6 * (2.0 * x11 * x157 * x17 * x4 * x51 * x6 - x178)
    x467 = x180 * x254 + 0.5 * x466
    x468 = x467 * x5 + x6 * (
        2.0 * x11 * x157 * x17 * x254 * x4 * x5 * x51 * x6 - x158 * x274
    )
    x469 = x85 * (-x463 * x60 + x468)
    x470 = x116 * x459 + x3 * x443 + x469
    x471 = x254 * x6 * (2.0 * x11 * x157 * x17 * x4 * x6 * x66 - x225)
    x472 = x227 * x254 + 0.5 * x471
    x473 = x116 * x467 + x3 * x463
    x474 = x85 * (-x467 * x60 + x472)
    x475 = x151 * x157 * x254 + x220 * x315
    x476 = x153 * x475 + x3 * x459 + x474
    x477 = x85 * (-x442 * x60 + x467)
    x478 = x219 * x254 + x220 * x317
    x479 = x153 * x478 + x3 * x457 + x477
    x480 = x111 * x157
    x481 = -x85 * (x346 * x480 + x440)
    x482 = 24.3267121527398 * x154
    x483 = x254**2
    x484 = x37 + x38 * x483
    x485 = x31 + x32 * x483
    x486 = x485 * x5
    x487 = x24 * x483
    x488 = x20 + x21 * x487
    x489 = x22 * x488
    x490 = x484 - x485 * x7
    x491 = x19 * x490
    x492 = x489 + x491
    x493 = x492 * x5 + x6 * (x484 * x5 - x486 * x7)
    x494 = x3 * x493
    x495 = x22 * x485
    x496 = x483 * x53 + x52
    x497 = -x484 * x7 + x496
    x498 = x19 * x497
    x499 = x495 + x498
    x500 = x49 * x499
    x501 = 3.0 * x500
    x502 = x494 + x501
    x503 = -x496 * x5
    x504 = x484 * x5
    x505 = x499 * x5 - x6 * (x503 + x504 * x7)
    x506 = x3 * x505
    x507 = x22 * x484
    x508 = x483 * x68 + x67
    x509 = -x496 * x7 + x508
    x510 = x19 * x509
    x511 = x507 + x510
    x512 = x49 * x511
    x513 = 3.0 * x512
    x514 = x506 + x513
    x515 = x488 * x7
    x516 = x487 * x78 + x77
    x517 = x22 * x516
    x518 = x485 - x515
    x519 = x19 * x518
    x520 = x5 * (x517 + x519 + x6 * (x485 - x515))
    x521 = x3 * x520
    x522 = x49 * x492
    x523 = 3.0 * x522
    x524 = x493 * x60
    x525 = x85 * (x505 - x524)
    x526 = x3 * x492
    x527 = x153 * x486
    x528 = x526 + x527
    x529 = x499 * x60
    x530 = x85 * (x511 - x529)
    x531 = x484 * x49
    x532 = x3 * x486
    x533 = x531 + x532
    x534 = x153 * x533 + x3 * x528 + x530
    x535 = x3 * x499
    x536 = x153 * x504
    x537 = x535 + x536
    x538 = x3 * x511
    x539 = x496 * x5
    x540 = x153 * x539 + x538
    x541 = x504 * x60
    x542 = -x85 * (x503 + x541)
    x543 = x3 * x531
    x544 = x3 * x533 + x542 + x543
    x545 = x153 * x544 + x3 * x534 - x4 * (x537 * x60 - x540)
    x546 = x5 * x511 + x6 * (x5 * x508 - x539 * x7)
    x547 = x85 * (-x505 * x60 + x546)
    x548 = x116 * x537 + x3 * x502 + x547
    x549 = x126 + x127 * x483
    x550 = x5 * x508
    x551 = x22 * x496
    x552 = -x508 * x7 + x549
    x553 = x19 * x552
    x554 = x551 + x553
    x555 = x5 * x554 + x6 * (x5 * x549 - x550 * x7)
    x556 = x85 * (-x546 * x60 + x555)
    x557 = x116 * x540 + x3 * x514 + x556
    x558 = x135 + x136 * x483 - x549 * x7
    x559 = x19 * x558
    x560 = x22 * x508 + x559
    x561 = x49 * x560
    x562 = 3.0 * x561
    x563 = x49 * x554
    x564 = 3.0 * x563
    x565 = x3 * x546 + x564
    x566 = x85 * (-x554 * x60 + x560)
    x567 = x49 * x508
    x568 = x3 * x539 + x567
    x569 = x153 * x568 + x3 * x540 + x566
    x570 = x85 * (-x511 * x60 + x554)
    x571 = x49 * x496
    x572 = x3 * x504
    x573 = x571 + x572
    x574 = x153 * x573 + x3 * x537 + x570
    x575 = x85 * (-x496 * x60 + x508)
    x576 = x157 * x354 + x160
    x577 = x157 * x350 + x165
    x578 = x157 * x349 + x169
    x579 = -x577 * x7 + x578
    x580 = x19 * x579 + x22 * x576
    x581 = x3 * x580
    x582 = x49 * x577
    x583 = x5 * x582
    x584 = x581 + 2.0 * x583
    x585 = x157 * x359 + x179
    x586 = -x578 * x7 + x585
    x587 = x19 * x586 + x22 * x577
    x588 = x3 * x587
    x589 = x5 * x578
    x590 = x153 * x589 + x588
    x591 = x157 * x377 + x189
    x592 = -x576 * x7 + x577
    x593 = x19 * x592 + x22 * x591
    x594 = x3 * x593
    x595 = x5 * x576
    x596 = x580 * x60
    x597 = x85 * (x587 - x596)
    x598 = x152 * x576
    x599 = x582 + x598
    x600 = x5 * x60
    x601 = x85 * (x5 * x578 - x577 * x600)
    x602 = x3 * x582
    x603 = x3 * x599 + x601 + x602
    x604 = x49 * x578
    x605 = x152 * x577
    x606 = x604 + x605
    x607 = x49 * x585
    x608 = x3 * x589
    x609 = x607 + x608
    x610 = x429 * x577
    x611 = x578 * x60
    x612 = x85 * (x585 - x611)
    x613 = x610 + x612
    x614 = x3 * x603 - x4 * (x60 * x606 - x609) + x49 * x613
    x615 = x157 * x369 + x226
    x616 = -x585 * x7 + x615
    x617 = x19 * x616 + x22 * x578
    x618 = x85 * (-x587 * x60 + x617)
    x619 = x153 * x606 + x3 * x584 + x618
    x620 = x157 * x408 + x233 - x615 * x7
    x621 = x19 * x620 + x22 * x585
    x622 = x85 * (-x60 * x617 + x621)
    x623 = x153 * x609 + x3 * x590 + x622
    x624 = x5 * x607
    x625 = x3 * x617 + 2.0 * x624
    x626 = x5 * x615
    x627 = x5 * x85 * (-x585 * x60 + x615)
    x628 = x3 * x607 + x3 * x609 + x627
    x629 = x85 * (x5 * x585 - x589 * x60)
    x630 = x3 * x604 + x3 * x606 + x629
    x631 = 14.0450338098829 * x154
    x632 = x254 * x353 + x257
    x633 = x262 + x348 * x453
    x634 = x266 + x317 * x348
    x635 = -x633 * x7 + x634
    x636 = x19 * x635 + x22 * x632
    x637 = x49 * x633
    x638 = 2.0 * x637
    x639 = x5 * x638
    x640 = x3 * x636 + x639
    x641 = x276 + x315 * x348
    x642 = -x634 * x7 + x641
    x643 = x19 * x642 + x22 * x633
    x644 = x5 * x634
    x645 = x153 * x644
    x646 = x3 * x643 + x645
    x647 = x254 * x376 + x286
    x648 = -x632 * x7 + x633
    x649 = x19 * x648 + x22 * x647
    x650 = x5 * x632
    x651 = x153 * x650
    x652 = x85 * (-x60 * x636 + x643)
    x653 = x152 * x632 + x637
    x654 = x85 * (x5 * x634 - x600 * x633)
    x655 = x3 * x637 + x3 * x653 + x654
    x656 = x49 * x634
    x657 = x152 * x633 + x656
    x658 = x49 * x641
    x659 = x3 * x644 + x658
    x660 = x60 * x634
    x661 = x85 * (x641 - x660)
    x662 = x429 * x633 + x661
    x663 = x3 * x655 - x4 * (x60 * x657 - x659) + x49 * x662
    x664 = x254 * x348 * x68 + x323
    x665 = -x641 * x7 + x664
    x666 = x19 * x665 + x22 * x634
    x667 = x85 * (-x60 * x643 + x666)
    x668 = x153 * x657 + x3 * x640 + x667
    x669 = x127 * x254 * x348 + x330 - x664 * x7
    x670 = x19 * x669 + x22 * x641
    x671 = x85 * (-x60 * x666 + x670)
    x672 = x153 * x659 + x3 * x646 + x671
    x673 = 2.0 * x658
    x674 = x5 * x673
    x675 = x3 * x666 + x674
    x676 = x5 * x664
    x677 = x153 * x676
    x678 = x5 * x85 * (-x60 * x641 + x664)
    x679 = x3 * x658 + x3 * x659 + x678
    x680 = x85 * (x5 * x641 - x60 * x644)
    x681 = x3 * x656 + x3 * x657 + x680
    x682 = 31.4056503451809 * x154
    x683 = x157 * x485
    x684 = x6 * (x157 * x484 - x683 * x7)
    x685 = x157 * x489 + 0.5 * x684
    x686 = x157 * x527
    x687 = x3 * x685 + x686
    x688 = -x157 * x496
    x689 = x157 * x484
    x690 = -x6 * (x688 + x689 * x7)
    x691 = x157 * x495 + 0.5 * x690
    x692 = x157 * x504
    x693 = x153 * x692
    x694 = x3 * x691 + x693
    x695 = x157 * x6 * (x485 - x515)
    x696 = x157 * x517 + 0.5 * x695
    x697 = x158 * x488
    x698 = x153 * x697
    x699 = x85 * (-x60 * x685 + x691)
    x700 = x49 * x683
    x701 = x220 * x488 + x700
    x702 = x157 * x486
    x703 = x85 * (x157 * x484 * x5 - x60 * x702)
    x704 = x3 * x700 + x3 * x701 + x703
    x705 = x157 * x531
    x706 = x157 * x532 + x705
    x707 = x157 * x571
    x708 = x157 * x572 + x707
    x709 = x60 * x689
    x710 = -x85 * (x688 + x709)
    x711 = x429 * x683 + x710
    x712 = x3 * x704 - x4 * (x60 * x706 - x708) + x49 * x711
    x713 = x157 * x496
    x714 = x6 * (x157 * x508 - x7 * x713)
    x715 = x157 * x507 + 0.5 * x714
    x716 = x85 * (-x60 * x691 + x715)
    x717 = x153 * x706 + x3 * x687 + x716
    x718 = x157 * x6 * (-x508 * x7 + x549)
    x719 = x157 * x551 + 0.5 * x718
    x720 = x85 * (-x60 * x715 + x719)
    x721 = x153 * x708 + x3 * x694 + x720
    x722 = x157 * x539
    x723 = x153 * x722
    x724 = x3 * x715 + x723
    x725 = x157 * x550
    x726 = x153 * x725
    x727 = x85 * (x157 * x5 * x508 - x60 * x722)
    x728 = x3 * x707 + x3 * x708 + x727
    x729 = x85 * (x157 * x496 * x5 - x60 * x692)
    x730 = x157 * x543 + x3 * x706 + x729
    x731 = x254 * x488 + x256
    x732 = x254 * x485 + x261
    x733 = x254 * x484 + x265
    x734 = -x7 * x732 + x733
    x735 = x19 * x734
    x736 = x22 * x731 + x735
    x737 = x3 * x736
    x738 = x49 * x732
    x739 = x5 * x738
    x740 = x737 + 2.0 * x739
    x741 = x254 * x496 + x275
    x742 = -x7 * x733 + x741
    x743 = x19 * x742
    x744 = x22 * x732 + x743
    x745 = x3 * x744
    x746 = x5 * x733
    x747 = x153 * x746 + x745
    x748 = x254 * x516 + x285
    x749 = -x7 * x731 + x732
    x750 = x19 * x749
    x751 = x22 * x748 + x750
    x752 = x3 * x751
    x753 = x5 * x731
    x754 = x60 * x736
    x755 = x85 * (x744 - x754)
    x756 = x152 * x731
    x757 = x738 + x756
    x758 = x85 * (x5 * x733 - x600 * x732)
    x759 = x3 * x738
    x760 = x3 * x757 + x758 + x759
    x761 = x49 * x733
    x762 = x152 * x732
    x763 = x761 + x762
    x764 = x49 * x741
    x765 = x3 * x746
    x766 = x764 + x765
    x767 = x429 * x732
    x768 = x60 * x733
    x769 = x85 * (x741 - x768)
    x770 = x767 + x769
    x771 = x3 * x760 - x4 * (x60 * x763 - x766) + x49 * x770
    x772 = x254 * x508 + x322
    x773 = -x7 * x741 + x772
    x774 = x19 * x773
    x775 = x22 * x733 + x774
    x776 = x85 * (-x60 * x744 + x775)
    x777 = x153 * x763 + x3 * x740 + x776
    x778 = x254 * x549 + x329 - x7 * x772
    x779 = x19 * x778
    x780 = x22 * x741 + x779
    x781 = x85 * (-x60 * x775 + x780)
    x782 = x153 * x766 + x3 * x747 + x781
    x783 = x5 * x764
    x784 = x3 * x775 + 2.0 * x783
    x785 = x5 * x772
    x786 = x5 * x741
    x787 = x85 * (x5 * x772 - x60 * x786)
    x788 = x3 * x764 + x3 * x766 + x787
    x789 = x85 * (x5 * x741 - x60 * x746)
    x790 = x3 * x761
    x791 = x3 * x763 + x789 + x790
    x792 = x157 * x578 + x370 * x41
    x793 = x157 * x577 + x360 * x41
    x794 = x60 * x793
    x795 = x3 * x794
    x796 = x157 * x576 + x355 * x41
    x797 = x429 * x796
    x798 = x85 * (x792 - x794)
    x799 = x797 + x798
    x800 = x3 * x799 + x4 * (x3 * x792 - x795)
    x801 = x157 * x585 + x409 * x41
    x802 = x85 * (-x60 * x792 + x801)
    x803 = x429 * x793 + x802
    x804 = x157 * x615 + x41 * x414
    x805 = x85 * (-x60 * x801 + x804)
    x806 = x429 * x792 + x805
    x807 = x49 * x793
    x808 = x152 * x796
    x809 = x807 + x808
    x810 = x49 * x792
    x811 = x5 * x793
    x812 = x3 * x811
    x813 = x810 + x812
    x814 = x49 * x796
    x815 = x157 * x591 + x378 * x41
    x816 = x152 * x815
    x817 = x85 * (x5 * x793 - x600 * x796)
    x818 = x3 * x814
    x819 = x85 * (x5 * x792 - x60 * x811)
    x820 = x3 * x807
    x821 = x3 * x809 + x819 + x820
    x822 = x5 * x792
    x823 = x85 * (x5 * x801 - x60 * x822)
    x824 = x3 * x792
    x825 = x49 * x824
    x826 = x3 * x813 + x823 + x825
    x827 = x49 * x804
    x828 = x5 * x801
    x829 = x49 * x801
    x830 = x5 * x824 + x829
    x831 = x157 * x634 + x441
    x832 = x157 * x633 + x435
    x833 = x60 * x832
    x834 = x157 * x632 + x432
    x835 = x85 * (x831 - x833)
    x836 = x429 * x834 + x835
    x837 = x3 * (x4 * (x831 - x833) + x836)
    x838 = x157 * x641 + x466
    x839 = x85 * (-x60 * x831 + x838)
    x840 = x429 * x832 + x839
    x841 = x157 * x664 + x471
    x842 = x85 * (-x60 * x838 + x841)
    x843 = x429 * x831 + x842
    x844 = x49 * x832
    x845 = x152 * x834 + x844
    x846 = x49 * x831
    x847 = x5 * x832
    x848 = x3 * x847 + x846
    x849 = x49 * x834
    x850 = x157 * x647 + x444
    x851 = x85 * (x5 * x832 - x600 * x834)
    x852 = x85 * (x5 * x831 - x60 * x847)
    x853 = x3 * x844 + x3 * x845 + x852
    x854 = x5 * x831
    x855 = x85 * (x5 * x838 - x60 * x854)
    x856 = x3 * x831
    x857 = x3 * x848 + x49 * x856 + x855
    x858 = x49 * x841
    x859 = x5 * x838
    x860 = x49 * x838
    x861 = x5 * x856 + x860
    x862 = x348 * x484 + x510
    x863 = x348 * x485 + x498
    x864 = x60 * x863
    x865 = x348 * x488 + x491
    x866 = x85 * (x862 - x864)
    x867 = x429 * x865 + x866
    x868 = x3 * (x4 * (x862 - x864) + x867)
    x869 = x348 * x496 + x553
    x870 = x85 * (-x60 * x862 + x869)
    x871 = x429 * x863 + x870
    x872 = x348 * x508 + x559
    x873 = x85 * (-x60 * x869 + x872)
    x874 = x429 * x862 + x873
    x875 = x49 * x863
    x876 = x152 * x865 + x875
    x877 = x49 * x862
    x878 = x5 * x863
    x879 = x3 * x878 + x877
    x880 = x49 * x865
    x881 = x348 * x516 + x519
    x882 = x85 * (x5 * x863 - x600 * x865)
    x883 = x85 * (x5 * x862 - x60 * x878)
    x884 = x3 * x875 + x3 * x876 + x883
    x885 = x5 * x862
    x886 = x85 * (x5 * x869 - x60 * x885)
    x887 = x3 * x862
    x888 = x3 * x879 + x49 * x887 + x886
    x889 = x49 * x872
    x890 = x5 * x869
    x891 = x49 * x869
    x892 = x5 * x887 + x891
    x893 = x157 * x732
    x894 = x60 * x893
    x895 = x85 * (x157 * x733 - x894)
    x896 = x157 * x731
    x897 = x429 * x896 + x895
    x898 = x3 * (x4 * (x157 * x733 - x894) + x897)
    x899 = x157 * x733
    x900 = x85 * (x157 * x741 - x60 * x899)
    x901 = x157 * x767 + x900
    x902 = x157 * x85 * (-x60 * x741 + x772)
    x903 = x429 * x899 + x902
    x904 = x157 * x738
    x905 = x157 * x756 + x904
    x906 = x157 * x761
    x907 = x157 * x762 + x906
    x908 = x49 * x896
    x909 = x85 * (x157 * x5 * x732 - x600 * x896)
    x910 = x5 * x893
    x911 = x85 * (x157 * x5 * x733 - x60 * x910)
    x912 = x157 * x759 + x3 * x905 + x911
    x913 = x157 * x746
    x914 = x85 * (x157 * x5 * x741 - x60 * x913)
    x915 = x157 * x790 + x3 * x907 + x914
    x916 = x157 * x772
    x917 = x157 * x764
    x918 = x157 * x765 + x917
    x919 = x157 * x786
    x920 = x254 * x733 + x41 * x509
    x921 = x254 * x732 + x41 * x497
    x922 = x60 * x921
    x923 = x3 * x922
    x924 = x254 * x731 + x41 * x490
    x925 = x429 * x924
    x926 = x85 * (x920 - x922)
    x927 = x925 + x926
    x928 = x3 * x927 + x4 * (x3 * x920 - x923)
    x929 = x429 * x921
    x930 = x254 * x741 + x41 * x552
    x931 = x85 * (-x60 * x920 + x930)
    x932 = x929 + x931
    x933 = x254 * x772 + x41 * x558
    x934 = x85 * (-x60 * x930 + x933)
    x935 = x429 * x920 + x934
    x936 = x49 * x921
    x937 = x152 * x924
    x938 = x936 + x937
    x939 = x49 * x920
    x940 = x5 * x921
    x941 = x3 * x940
    x942 = x939 + x941
    x943 = x49 * x924
    x944 = x254 * x748 + x41 * x518
    x945 = x152 * x944
    x946 = x85 * (x5 * x921 - x600 * x924)
    x947 = x3 * x943
    x948 = x85 * (x5 * x920 - x60 * x940)
    x949 = x3 * x936
    x950 = x3 * x938 + x948 + x949
    x951 = x5 * x920
    x952 = x85 * (x5 * x930 - x60 * x951)
    x953 = x3 * x920
    x954 = x49 * x953
    x955 = x3 * x942 + x952 + x954
    x956 = x49 * x933
    x957 = x5 * x930
    x958 = x49 * x930
    x959 = x5 * x953 + x958
    x960 = x12 * x586 + x157 * x793
    x961 = x12 * x579 + x157 * x796
    x962 = x60 * x961
    x963 = x3 * x962
    x964 = x12 * x592 + x157 * x815
    x965 = x429 * x964
    x966 = x85 * (x960 - x962)
    x967 = x429 * x961
    x968 = x12 * x616 + x157 * x792
    x969 = x85 * (-x60 * x960 + x968)
    x970 = x967 + x969
    x971 = x429 * x960
    x972 = x12 * x620 + x157 * x801
    x973 = x85 * (-x60 * x968 + x972)
    x974 = x971 + x973
    x975 = x157 * x832 + x41 * x642
    x976 = x157 * x834 + x41 * x635
    x977 = x60 * x976
    x978 = x157 * x850 + x41 * x648
    x979 = x85 * (x975 - x977)
    x980 = x157 * x831 + x41 * x665
    x981 = x85 * (-x60 * x975 + x980)
    x982 = x429 * x976 + x981
    x983 = x157 * x838 + x41 * x669
    x984 = x85 * (-x60 * x980 + x983)
    x985 = x429 * x975 + x984
    x986 = x157 * x863 + x690
    x987 = x157 * x865 + x684
    x988 = x60 * x987
    x989 = x157 * x881 + x695
    x990 = x85 * (x986 - x988)
    x991 = x157 * x862 + x714
    x992 = x85 * (-x60 * x986 + x991)
    x993 = x429 * x987 + x992
    x994 = x157 * x869 + x718
    x995 = x85 * (-x60 * x991 + x994)
    x996 = x429 * x986 + x995
    x997 = x348 * x732 + x743
    x998 = x348 * x731 + x735
    x999 = x60 * x998
    x1000 = x348 * x748 + x750
    x1001 = x85 * (x997 - x999)
    x1002 = x348 * x733 + x774
    x1003 = x85 * (x1002 - x60 * x997)
    x1004 = x1003 + x429 * x998
    x1005 = x348 * x741 + x779
    x1006 = x85 * (-x1002 * x60 + x1005)
    x1007 = x1006 + x429 * x997
    x1008 = x157 * x924
    x1009 = x1008 * x60
    x1010 = x85 * (-x1009 + x157 * x921)
    x1011 = x157 * x944
    x1012 = x157 * x921
    x1013 = x85 * (-x1012 * x60 + x157 * x920)
    x1014 = x1013 + x157 * x925
    x1015 = x157 * x920
    x1016 = x85 * (-x1015 * x60 + x157 * x930)
    x1017 = x1016 + x157 * x929
    x1018 = x12 * x742 + x254 * x921
    x1019 = x12 * x734 + x254 * x924
    x1020 = x1019 * x60
    x1021 = x1020 * x3
    x1022 = x12 * x749 + x254 * x944
    x1023 = x1022 * x429
    x1024 = x85 * (x1018 - x1020)
    x1025 = x1019 * x429
    x1026 = x12 * x773 + x254 * x920
    x1027 = x85 * (-x1018 * x60 + x1026)
    x1028 = x1025 + x1027
    x1029 = x1018 * x429
    x1030 = x12 * x778 + x254 * x930
    x1031 = x85 * (-x1026 * x60 + x1030)
    x1032 = x1029 + x1031
    x1033 = -x156 - A[1]
    x1034 = x1033 * x48
    x1035 = x1033 * x58
    x1036 = x1034 + x1035
    x1037 = x1033 * x64
    x1038 = x1033 * x71
    x1039 = x1037 + x1038 * x91
    x1040 = x1033 * x81
    x1041 = x1033 * x83
    x1042 = x1033 * x84
    x1043 = x4 * (x1033 * x63 - x1042)
    x1044 = x1033 * x87
    x1045 = x1033 * x89
    x1046 = x1044 + x1045
    x1047 = x1033 * x92
    x1048 = x4 * (x1033 * x71 - x1047)
    x1049 = x1033 * x94
    x1050 = x1033 * x96
    x1051 = x1049 + x1050
    x1052 = x1046 * x3 + 0.5 * x1048 + x1051 * x98
    x1053 = x1033 * (x100 + x102)
    x1054 = x1033 * x62
    x1055 = x1033 * x104 + x1054 * x98
    x1056 = x1033 * x4 * (-x108 + x62)
    x1057 = x1033 * x5
    x1058 = x111 * x113
    x1059 = x1033 * x4 * (-x121 * x60 + x131)
    x1060 = x1033 * x4 * (x121 - x60 * x63)
    x1061 = x1033 * x172
    x1062 = x1061 + x57
    x1063 = x1033 * x174
    x1064 = x101 + x1063
    x1065 = x1062 * x3 + x1064 * x98
    x1066 = x1033 * x182
    x1067 = x1066 + x72
    x1068 = x1033 * x185
    x1069 = x105 + x1068
    x1070 = x1067 * x3 + x1069 * x98
    x1071 = x1033 * x191
    x1072 = x1071 + x82
    x1073 = x1033 * x163
    x1074 = x1073 + x88
    x1075 = x4 * (-x1062 * x60 + x1067)
    x1076 = x1033 * x167
    x1077 = x1076 + x95
    x1078 = x1074 * x3 + x1077 * x116
    x1079 = x4 * (-x1064 * x60 + x1069)
    x1080 = x5 * x53
    x1081 = x1080 * x49
    x1082 = x1057 * x157
    x1083 = x1082 * x38
    x1084 = x1081 + x1083
    x1085 = x1084 * x153
    x1086 = x1077 * x3 + x1085
    x1087 = x1078 * x3 + 0.5 * x1079 + x1086 * x116
    x1088 = x1033 * x171 + x146
    x1089 = x1064 * x3 + x1088 * x116
    x1090 = x1033 * x181 + x142
    x1091 = x1069 * x3 + x1090 * x116
    x1092 = x4 * (-x1088 * x60 + x1090)
    x1093 = x49 * (x1033 * x218 + x151)
    x1094 = x1033 * x229 + x139
    x1095 = x4 * (x1033 * x234 - x1094 * x60 + x138)
    x1096 = x4 * (-x1067 * x60 + x1094)
    x1097 = 27.5838988165724 * x154
    x1098 = x1033 * (x269 + x272)
    x1099 = x1033 * x281
    x1100 = x1033 * x279 + x1099 * x98
    x1101 = x1033 * x4 * (x278 - x291)
    x1102 = x1033 * (x293 + x295)
    x1103 = x1033 * x4 * (x281 - x297)
    x1104 = x1033 * x299 + x1057 * x18 * x301
    x1105 = x1102 * x3 + 0.5 * x1103 + x1104 * x116
    x1106 = x1033 * (x306 + x308)
    x1107 = x1033 * x277
    x1108 = x1033 * x310 + x1107 * x116
    x1109 = x1033 * x4 * (x277 - x313)
    x1110 = x1033 * x316
    x1111 = x1033 * x4 * (-x325 * x60 + x331)
    x1112 = x1033 * x4 * (-x278 * x60 + x325)
    x1113 = x1033 * x357
    x1114 = 2.0 * x175
    x1115 = x1113 + x1114
    x1116 = x1033 * x361
    x1117 = 2.0 * x210
    x1118 = x1116 + x1117
    x1119 = x1118 * x116
    x1120 = x1115 * x3 + x1119
    x1121 = x1033 * x367
    x1122 = x1121 + 2.0 * x186
    x1123 = x1033 * x371
    x1124 = x1123 + 2.0 * x214
    x1125 = x1124 * x116
    x1126 = x1122 * x3 + x1125
    x1127 = x1033 * x379
    x1128 = 2.0 * x193
    x1129 = x1127 + x1128
    x1130 = x1033 * x356
    x1131 = 2.0 * x198
    x1132 = x1130 + x1131
    x1133 = x1132 * x116
    x1134 = x4 * (-x1115 * x60 + x1122)
    x1135 = x1033 * x351
    x1136 = x1135 + x206
    x1137 = x1132 * x3 + x1136 * x153
    x1138 = x4 * (-x1118 * x60 + x1124)
    x1139 = x1033 * x349
    x1140 = x113 * x480
    x1141 = x1139 + x1140
    x1142 = x1141 * x49
    x1143 = x1136 * x3 + x1142
    x1144 = x1137 * x3 + 0.5 * x1138 + x1143 * x153
    x1145 = x1033 * x366 + x247
    x1146 = x1118 * x3 + x1145 * x153
    x1147 = x1033 * x398 + x242
    x1148 = x1124 * x3 + x1147 * x153
    x1149 = x4 * (-x1145 * x60 + x1147)
    x1150 = x1033 * x405 + 2.0 * x239
    x1151 = x4 * (x1033 * x411 - x1150 * x60 + 2.0 * x238)
    x1152 = x4 * (-x1122 * x60 + x1150)
    x1153 = 42.1351014296486 * x154
    x1154 = x1033 * x434 + x271
    x1155 = x1033 * x436 + x307
    x1156 = x1154 * x3 + x1155 * x116
    x1157 = x1033 * x439 + x282
    x1158 = x1033 * x442 + x311
    x1159 = x1157 * x3 + x1158 * x116
    x1160 = x1033 * x445 + x289
    x1161 = x1033 * x433 + x294
    x1162 = x4 * (-x1154 * x60 + x1157)
    x1163 = x1082 * x453 + x254 * x38 * x49 * x5
    x1164 = x1163 * x153
    x1165 = x1161 * x3 + x1164
    x1166 = x4 * (-x1155 * x60 + x1158)
    x1167 = x1033 * x451 + x316
    x1168 = x1167 * x49
    x1169 = x1163 * x3 + x1168
    x1170 = x1165 * x3 + 0.5 * x1166 + x1169 * x153
    x1171 = x1080 * x254 * x49 + x1082 * x317
    x1172 = x1171 * x153
    x1173 = x1155 * x3 + x1172
    x1174 = x1082 * x315 + x254 * x49 * x5 * x68
    x1175 = x1174 * x153
    x1176 = x1158 * x3 + x1175
    x1177 = x4 * (-x1171 * x60 + x1174)
    x1178 = x1033 * x463 + x336
    x1179 = x4 * (x1033 * x468 - x1178 * x60 + x335)
    x1180 = x4 * (-x1157 * x60 + x1178)
    x1181 = 72.9801364582193 * x154
    x1182 = x1033 * (x494 + x501)
    x1183 = x1033 * x511
    x1184 = x1033 * x506 + x116 * x1183
    x1185 = x1033 * x4 * (x505 - x524)
    x1186 = x1033 * (x526 + x527)
    x1187 = x1033 * x4 * (x511 - x529)
    x1188 = x1033 * x531
    x1189 = x1033 * x532 + x1188
    x1190 = x1186 * x3 + 0.5 * x1187 + x1189 * x153
    x1191 = x1033 * (x535 + x536)
    x1192 = x1033 * x539
    x1193 = x1033 * x538 + x1192 * x153
    x1194 = x1033 * x4 * (x496 * x5 - x541)
    x1195 = x1033 * x4 * (-x546 * x60 + x555)
    x1196 = x1033 * x4 * (-x505 * x60 + x546)
    x1197 = x1033 * x580
    x1198 = x1197 + x363
    x1199 = x1033 * x577
    x1200 = x1199 * x5
    x1201 = x116 * x366
    x1202 = x1200 + x1201
    x1203 = x1198 * x3 + x1202 * x153
    x1204 = x1033 * x587
    x1205 = x1204 + x373
    x1206 = x1033 * x589
    x1207 = x116 * x398 + x1206
    x1208 = x1205 * x3 + x1207 * x153
    x1209 = x1033 * x593
    x1210 = x1209 + x382
    x1211 = x1057 * x576
    x1212 = x116 * x351
    x1213 = x1211 + x1212
    x1214 = x4 * (-x1198 * x60 + x1205)
    x1215 = 3.0 * x390
    x1216 = x1199 + x1215
    x1217 = x1216 * x49
    x1218 = x1213 * x3 + x1217
    x1219 = x4 * (-x1202 * x60 + x1207)
    x1220 = x1217 * x3 + x1218 * x3 + 0.5 * x1219
    x1221 = x1033 * x578 + 3.0 * x426
    x1222 = x1221 * x49
    x1223 = x1202 * x3 + x1222
    x1224 = x1033 * x585 + 3.0 * x422
    x1225 = x1224 * x49
    x1226 = x1207 * x3 + x1225
    x1227 = x4 * (-x1221 * x60 + x1224)
    x1228 = x1033 * x617 + x419
    x1229 = x4 * (x1033 * x621 - x1228 * x60 + x417)
    x1230 = x4 * (-x1205 * x60 + x1228)
    x1231 = 42.1351014296486 * x154
    x1232 = x153 * x436
    x1233 = x1033 * x636 + x1232
    x1234 = x1033 * x633
    x1235 = x1234 * x5 + x456
    x1236 = x1233 * x3 + x1235 * x153
    x1237 = x153 * x442
    x1238 = x1033 * x643 + x1237
    x1239 = x1033 * x644 + x458
    x1240 = x1238 * x3 + x1239 * x153
    x1241 = x153 * x433
    x1242 = x1033 * x649 + x1241
    x1243 = x1057 * x632 + x447
    x1244 = x4 * (-x1233 * x60 + x1238)
    x1245 = x157 * x301
    x1246 = x1245 * x18
    x1247 = x1234 + x1246
    x1248 = x1247 * x49
    x1249 = x1243 * x3 + x1248
    x1250 = x4 * (-x1235 * x60 + x1239)
    x1251 = x1248 * x3 + x1249 * x3 + 0.5 * x1250
    x1252 = x1140 * x254
    x1253 = x1033 * x634 + x1252
    x1254 = x1253 * x49
    x1255 = x1235 * x3 + x1254
    x1256 = x1245 * x36
    x1257 = x1033 * x641 + x1256
    x1258 = x1257 * x49
    x1259 = x1239 * x3 + x1258
    x1260 = x4 * (-x1253 * x60 + x1257)
    x1261 = x153 * x472
    x1262 = x153 * x467
    x1263 = x1033 * x666 + x1262
    x1264 = x4 * (x1033 * x670 + x1261 - x1263 * x60)
    x1265 = x4 * (-x1238 * x60 + x1263)
    x1266 = 94.2169510355428 * x154
    x1267 = x1033 * x685 + x500
    x1268 = x49 * x504
    x1269 = x1033 * x702 + x1268
    x1270 = x1269 * x153
    x1271 = x1267 * x3 + x1270
    x1272 = x1033 * x691 + x512
    x1273 = x1033 * x692 + x49 * x539
    x1274 = x1273 * x153
    x1275 = x1272 * x3 + x1274
    x1276 = x1033 * x696 + x522
    x1277 = x486 * x49
    x1278 = x1082 * x488 + x1277
    x1279 = x1278 * x153
    x1280 = x4 * (-x1267 * x60 + x1272)
    x1281 = x1033 * x683 + x531
    x1282 = x1281 * x49
    x1283 = x1278 * x3 + x1282
    x1284 = x4 * (-x1269 * x60 + x1273)
    x1285 = x1282 * x3 + x1283 * x3 + 0.5 * x1284
    x1286 = x1033 * x689 + x571
    x1287 = x1286 * x49
    x1288 = x1269 * x3 + x1287
    x1289 = x1033 * x713 + x567
    x1290 = x1289 * x49
    x1291 = x1273 * x3 + x1290
    x1292 = x4 * (-x1286 * x60 + x1289)
    x1293 = x1033 * x715 + x563
    x1294 = x4 * (x1033 * x719 - x1293 * x60 + x561)
    x1295 = x4 * (-x1272 * x60 + x1293)
    x1296 = x1057 * x738
    x1297 = x1033 * x737 + 2.0 * x1296
    x1298 = x1033 * x746
    x1299 = x1033 * x745 + x1298 * x153
    x1300 = x1057 * x731
    x1301 = x1033 * x4 * (x744 - x754)
    x1302 = x1033 * x738
    x1303 = x1033 * x756 + x1302
    x1304 = x1057 * x60
    x1305 = x4 * (x1033 * x5 * x733 - x1304 * x732)
    x1306 = x1033 * x759 + x1303 * x3 + 0.5 * x1305
    x1307 = x1033 * x761
    x1308 = x1033 * x762 + x1307
    x1309 = x1033 * x764
    x1310 = x1033 * x765 + x1309
    x1311 = x1033 * x4 * (x741 - x768)
    x1312 = x1033 * x4 * (-x60 * x775 + x780)
    x1313 = x1033 * x4 * (-x60 * x744 + x775)
    x1314 = x1033 * x792
    x1315 = 4.0 * x607
    x1316 = x1314 + x1315
    x1317 = x1033 * x793
    x1318 = 4.0 * x604
    x1319 = x1317 + x1318
    x1320 = x1319 * x60
    x1321 = x1033 * x796
    x1322 = 4.0 * x582
    x1323 = x1321 + x1322
    x1324 = x4 * (x1316 - x1320)
    x1325 = x1323 * x429 + 0.5 * x1324
    x1326 = x1319 * x49
    x1327 = x1321 * x5
    x1328 = x1322 * x5 + x1327
    x1329 = x1326 + x1328 * x3
    x1330 = x1316 * x49
    x1331 = x1033 * x811
    x1332 = x1331 + x589 * x98
    x1333 = x1330 + x1332 * x3
    x1334 = x1323 * x49
    x1335 = x1057 * x815
    x1336 = x1335 + x595 * x98
    x1337 = x4 * (-x1328 * x60 + x1332)
    x1338 = x1033 * x822 + x1315 * x5
    x1339 = x4 * (x1033 * x828 - x1338 * x60 + x626 * x98)
    x1340 = x4 * (-x1332 * x60 + x1338)
    x1341 = 3.0 * x658
    x1342 = x1033 * x831 + x1341
    x1343 = x1033 * x832 + 3.0 * x656
    x1344 = x1343 * x60
    x1345 = x1033 * x834
    x1346 = 3.0 * x637
    x1347 = x1345 + x1346
    x1348 = x4 * (x1342 - x1344)
    x1349 = x1347 * x429 + 0.5 * x1348
    x1350 = x1343 * x49
    x1351 = x5 * (x1345 + x1346)
    x1352 = x1350 + x1351 * x3
    x1353 = x1342 * x49
    x1354 = x1033 * x847 + x116 * x644
    x1355 = x1353 + x1354 * x3
    x1356 = x1347 * x49
    x1357 = x1057 * x850 + x116 * x650
    x1358 = x4 * (-x1351 * x60 + x1354)
    x1359 = x1033 * x854 + x1341 * x5
    x1360 = x4 * (x1033 * x859 + x116 * x676 - x1359 * x60)
    x1361 = x4 * (-x1354 * x60 + x1359)
    x1362 = x1033 * x862 + 2.0 * x707
    x1363 = x1033 * x863 + 2.0 * x705
    x1364 = x1363 * x60
    x1365 = x1033 * x865
    x1366 = x1365 + x153 * x683
    x1367 = x4 * (x1362 - x1364)
    x1368 = x1366 * x429 + 0.5 * x1367
    x1369 = x1363 * x49
    x1370 = x1365 * x5 + x686
    x1371 = x1369 + x1370 * x3
    x1372 = x1362 * x49
    x1373 = x1033 * x878 + x693
    x1374 = x1372 + x1373 * x3
    x1375 = x1366 * x49
    x1376 = x1057 * x881 + x698
    x1377 = x4 * (-x1370 * x60 + x1373)
    x1378 = x1033 * x885 + x723
    x1379 = x4 * (x1033 * x890 - x1378 * x60 + x726)
    x1380 = x4 * (-x1373 * x60 + x1378)
    x1381 = x1033 * x899 + x764
    x1382 = x1033 * x893 + x761
    x1383 = x1382 * x60
    x1384 = x1033 * x896 + x738
    x1385 = x4 * (x1381 - x1383)
    x1386 = x1384 * x429 + 0.5 * x1385
    x1387 = x1382 * x49
    x1388 = x1057 * x896 + x739
    x1389 = x1387 + x1388 * x3
    x1390 = x1381 * x49
    x1391 = x1033 * x910 + x49 * x746
    x1392 = x1390 + x1391 * x3
    x1393 = x1384 * x49
    x1394 = x1082 * x748 + x49 * x753
    x1395 = x4 * (-x1388 * x60 + x1391)
    x1396 = x1033 * x913 + x783
    x1397 = x4 * (x1033 * x919 - x1396 * x60 + x49 * x785)
    x1398 = x4 * (-x1391 * x60 + x1396)
    x1399 = x1033 * x4 * (x920 - x922)
    x1400 = x1033 * x925 + 0.5 * x1399
    x1401 = x1033 * x936
    x1402 = x1033 * x937 + x1401
    x1403 = x1033 * x939
    x1404 = x1033 * x941 + x1403
    x1405 = x1033 * x943
    x1406 = x4 * (x1033 * x5 * x921 - x1304 * x924)
    x1407 = x1033 * x4 * (x5 * x930 - x60 * x951)
    x1408 = x1033 * x4 * (x5 * x920 - x60 * x940)
    x1409 = x1033 * x960
    x1410 = 5.0 * x810
    x1411 = x1409 + x1410
    x1412 = x1033 * x961
    x1413 = 5.0 * x807
    x1414 = x1412 + x1413
    x1415 = x1414 * x60
    x1416 = x1033 * x964
    x1417 = 5.0 * x814
    x1418 = x1416 + x1417
    x1419 = x4 * (x1411 - x1415)
    x1420 = x1033 * x968 + 5.0 * x829
    x1421 = x4 * (x1033 * x972 - x1420 * x60 + 5.0 * x827)
    x1422 = x4 * (-x1411 * x60 + x1420)
    x1423 = x1033 * x975 + 4.0 * x846
    x1424 = x1033 * x976 + 4.0 * x844
    x1425 = x1424 * x60
    x1426 = x1033 * x978 + 4.0 * x849
    x1427 = x4 * (x1423 - x1425)
    x1428 = x1033 * x980 + 4.0 * x860
    x1429 = x4 * (x1033 * x983 - x1428 * x60 + 4.0 * x858)
    x1430 = x4 * (-x1423 * x60 + x1428)
    x1431 = 3.0 * x877
    x1432 = x1033 * x986 + x1431
    x1433 = 3.0 * x875
    x1434 = x1033 * x987 + x1433
    x1435 = x1434 * x60
    x1436 = 3.0 * x880
    x1437 = x1033 * x989 + x1436
    x1438 = x4 * (x1432 - x1435)
    x1439 = 3.0 * x889
    x1440 = 3.0 * x891
    x1441 = x1033 * x991 + x1440
    x1442 = x4 * (x1033 * x994 + x1439 - x1441 * x60)
    x1443 = x4 * (-x1432 * x60 + x1441)
    x1444 = x1033 * x997 + 2.0 * x906
    x1445 = x1033 * x998 + 2.0 * x904
    x1446 = x1445 * x60
    x1447 = x1000 * x1033 + x153 * x896
    x1448 = x4 * (x1444 - x1446)
    x1449 = x1002 * x1033 + 2.0 * x917
    x1450 = x4 * (x1005 * x1033 - x1449 * x60 + x153 * x916)
    x1451 = x4 * (-x1444 * x60 + x1449)
    x1452 = x1012 * x1033 + x939
    x1453 = x1008 * x1033 + x936
    x1454 = x1453 * x60
    x1455 = x1011 * x1033 + x943
    x1456 = x4 * (x1452 - x1454)
    x1457 = x157 * x930
    x1458 = x1015 * x1033 + x958
    x1459 = x4 * (x1033 * x1457 - x1458 * x60 + x956)
    x1460 = x4 * (-x1452 * x60 + x1458)
    x1461 = x1033 * x4 * (x1018 - x1020)
    x1462 = x1033 * x4 * (-x1026 * x60 + x1030)
    x1463 = x1033 * x4 * (-x1018 * x60 + x1026)
    x1464 = -x253 - A[2]
    x1465 = x1464 * (x48 + x58)
    x1466 = x1464 * (x64 + x73)
    x1467 = x1464 * x4 * (x63 - x84)
    x1468 = 0.5 * x1467
    x1469 = x1464 * (x87 + x89)
    x1470 = x1464 * x4 * (x71 - x92)
    x1471 = 0.5 * x1470
    x1472 = x1464 * (x94 + x96)
    x1473 = x1469 * x3 + x1471 + x1472 * x98
    x1474 = x1464 * (x100 + x102)
    x1475 = x1464 * (x104 + x106)
    x1476 = x1464 * x4 * (-x108 + x62)
    x1477 = 0.5 * x1476
    x1478 = x1464 * x5
    x1479 = x1464 * x4 * (-x121 * x60 + x131)
    x1480 = 0.5 * x1479
    x1481 = x1464 * x63
    x1482 = x4 * (x121 * x1464 - x1481 * x60)
    x1483 = 0.5 * x1482
    x1484 = x1464 * (x173 + x176)
    x1485 = x1464 * x185
    x1486 = x1464 * x183 + x1485 * x98
    x1487 = x1464 * x4 * (x182 - x195)
    x1488 = 0.5 * x1487
    x1489 = x1464 * (x197 + x199)
    x1490 = x1464 * x4 * (x185 - x201)
    x1491 = 0.5 * x1490
    x1492 = x1478 * x157
    x1493 = x11 * x113 * x1492 * x18
    x1494 = x1464 * x203 + x1493
    x1495 = x116 * x1494 + x1489 * x3 + x1491
    x1496 = x1464 * (x209 + x211)
    x1497 = x1464 * x181
    x1498 = x116 * x1497 + x1464 * x213
    x1499 = x1464 * x4 * (x181 - x216)
    x1500 = 0.5 * x1499
    x1501 = x1464 * x219
    x1502 = x1464 * x4 * (-x229 * x60 + x234)
    x1503 = 0.5 * x1502
    x1504 = x1464 * x4 * (-x182 * x60 + x229)
    x1505 = 0.5 * x1504
    x1506 = x1464 * x268 + x57
    x1507 = x1506 * x3
    x1508 = x101 + x1464 * x270
    x1509 = x1508 * x49
    x1510 = 4.0 * x1509
    x1511 = x1507 + x1510
    x1512 = x1464 * x278 + x72
    x1513 = x1512 * x3
    x1514 = x105 + x1464 * x281
    x1515 = x1514 * x49
    x1516 = x1513 + 4.0 * x1515
    x1517 = x1464 * x287 + x82
    x1518 = x1517 * x3
    x1519 = x1464 * x259 + x88
    x1520 = x1519 * x49
    x1521 = 4.0 * x1520
    x1522 = x1506 * x60
    x1523 = x4 * (x1512 - x1522)
    x1524 = 0.5 * x1523
    x1525 = x1519 * x3
    x1526 = x1464 * x263 + x95
    x1527 = x1526 * x49
    x1528 = 3.0 * x1527
    x1529 = x1525 + x1528
    x1530 = x1508 * x60
    x1531 = x4 * (x1514 - x1530)
    x1532 = 0.5 * x1531
    x1533 = x1526 * x3
    x1534 = x1081 + x1478 * x317
    x1535 = x1534 * x49
    x1536 = 2.0 * x1535
    x1537 = x1533 + x1536
    x1538 = x116 * x1537 + x1529 * x3 + x1532
    x1539 = x146 + x1464 * x267
    x1540 = x1539 * x49
    x1541 = x1508 * x3 + 3.0 * x1540
    x1542 = x142 + x1464 * x277
    x1543 = x1542 * x49
    x1544 = x1514 * x3 + 3.0 * x1543
    x1545 = x4 * (-x1539 * x60 + x1542)
    x1546 = 0.5 * x1545
    x1547 = x49 * (x1464 * x315 + x151)
    x1548 = x139 + x1464 * x325
    x1549 = x4 * (x138 + x1464 * x331 - x1548 * x60)
    x1550 = 0.5 * x1549
    x1551 = x4 * (-x1512 * x60 + x1548)
    x1552 = 0.5 * x1551
    x1553 = x1464 * x363
    x1554 = x1464 * x358 + x1553
    x1555 = x1464 * x371
    x1556 = x116 * x1555
    x1557 = x1464 * x368 + x1556
    x1558 = x1464 * x382
    x1559 = x1464 * x4 * (x367 - x383)
    x1560 = 0.5 * x1559
    x1561 = x1464 * (x385 + x386)
    x1562 = x1464 * x4 * (x371 - x388)
    x1563 = 0.5 * x1562
    x1564 = x1464 * x390
    x1565 = x1464 * x391 + x1564
    x1566 = x153 * x1565 + x1561 * x3 + x1563
    x1567 = x1464 * (x394 + x395)
    x1568 = x1464 * x398
    x1569 = x1464 * x397 + x153 * x1568
    x1570 = x1464 * x4 * (x359 * x5 - x400)
    x1571 = 0.5 * x1570
    x1572 = x1464 * x4 * (-x405 * x60 + x411)
    x1573 = 0.5 * x1572
    x1574 = x1464 * x4 * (-x367 * x60 + x405)
    x1575 = 0.5 * x1574
    x1576 = x1464 * x434 + x175
    x1577 = x1464 * x436 + x210
    x1578 = x116 * x1577 + x1576 * x3
    x1579 = x1464 * x439 + x186
    x1580 = x1464 * x442 + x214
    x1581 = x116 * x1580 + x1579 * x3
    x1582 = x1464 * x445 + x193
    x1583 = x1464 * x433 + x198
    x1584 = x4 * (-x1576 * x60 + x1579)
    x1585 = 0.5 * x1584
    x1586 = x1492 * x453 + x158 * x38 * x49
    x1587 = x153 * x1586
    x1588 = x1583 * x3 + x1587
    x1589 = x4 * (-x1577 * x60 + x1580)
    x1590 = 0.5 * x1589
    x1591 = x1464 * x451 + x219
    x1592 = x1591 * x49
    x1593 = x1586 * x3 + x1592
    x1594 = x153 * x1593 + x1588 * x3 + x1590
    x1595 = x1492 * x317 + x158 * x49 * x53
    x1596 = x153 * x1595
    x1597 = x1577 * x3 + x1596
    x1598 = x1492 * x315 + x151 * x158
    x1599 = x153 * x1598
    x1600 = x1580 * x3 + x1599
    x1601 = x4 * (-x1595 * x60 + x1598)
    x1602 = 0.5 * x1601
    x1603 = x1464 * x463 + x239
    x1604 = x4 * (x1464 * x468 - x1603 * x60 + x238)
    x1605 = 0.5 * x1604
    x1606 = x4 * (-x1579 * x60 + x1603)
    x1607 = 0.5 * x1606
    x1608 = x1464 * x493 + 2.0 * x271
    x1609 = x1608 * x3
    x1610 = x1464 * x499 + 2.0 * x307
    x1611 = x1610 * x49
    x1612 = 3.0 * x1611
    x1613 = x1609 + x1612
    x1614 = x1464 * x505 + 2.0 * x282
    x1615 = x1614 * x3
    x1616 = x1464 * x511 + 2.0 * x311
    x1617 = x1616 * x49
    x1618 = 3.0 * x1617
    x1619 = x1615 + x1618
    x1620 = x1464 * x520 + 2.0 * x289
    x1621 = x1620 * x3
    x1622 = x1464 * x492 + 2.0 * x294
    x1623 = x1622 * x49
    x1624 = 3.0 * x1623
    x1625 = x1608 * x60
    x1626 = x4 * (x1614 - x1625)
    x1627 = 0.5 * x1626
    x1628 = x1622 * x3
    x1629 = x1464 * x486
    x1630 = x1629 + x303
    x1631 = x1630 * x49
    x1632 = 2.0 * x1631
    x1633 = x1628 + x1632
    x1634 = x1610 * x60
    x1635 = x4 * (x1616 - x1634)
    x1636 = 0.5 * x1635
    x1637 = x111 * x300 + x1464 * x484
    x1638 = x1637 * x49
    x1639 = x1630 * x3
    x1640 = x1638 + x1639
    x1641 = x153 * x1640 + x1633 * x3 + x1636
    x1642 = x1464 * x504 + x343
    x1643 = x1642 * x49
    x1644 = x1610 * x3 + 2.0 * x1643
    x1645 = x1464 * x539 + x339
    x1646 = x1645 * x49
    x1647 = x1616 * x3 + 2.0 * x1646
    x1648 = x4 * (-x1642 * x60 + x1645)
    x1649 = 0.5 * x1648
    x1650 = x1464 * x546 + 2.0 * x336
    x1651 = x4 * (x1464 * x555 - x1650 * x60 + 2.0 * x335)
    x1652 = 0.5 * x1651
    x1653 = x4 * (-x1614 * x60 + x1650)
    x1654 = 0.5 * x1653
    x1655 = x1464 * x582
    x1656 = x1655 * x5
    x1657 = x1464 * x581 + 2.0 * x1656
    x1658 = x1464 * x589
    x1659 = x1464 * x588 + x153 * x1658
    x1660 = x1478 * x576
    x1661 = x1464 * x4 * (x587 - x596)
    x1662 = 0.5 * x1661
    x1663 = x1464 * x598 + x1655
    x1664 = x1478 * x60
    x1665 = x4 * (x1464 * x5 * x578 - x1664 * x577)
    x1666 = 0.5 * x1665
    x1667 = x1464 * x602 + x1663 * x3 + x1666
    x1668 = x1464 * x604
    x1669 = x1464 * x605 + x1668
    x1670 = x1464 * x607
    x1671 = x1464 * x608 + x1670
    x1672 = x1464 * x4 * (x585 - x611)
    x1673 = 0.5 * x1672
    x1674 = x1464 * x4 * (-x60 * x617 + x621)
    x1675 = 0.5 * x1674
    x1676 = x1464 * x4 * (-x587 * x60 + x617)
    x1677 = 0.5 * x1676
    x1678 = x1464 * x636 + x362
    x1679 = x366 * x49
    x1680 = x1464 * x633
    x1681 = x1679 + x1680 * x5
    x1682 = x153 * x1681
    x1683 = x1678 * x3 + x1682
    x1684 = x1464 * x643 + x372
    x1685 = x1464 * x644 + x398 * x49
    x1686 = x153 * x1685
    x1687 = x1684 * x3 + x1686
    x1688 = x1464 * x649 + x381
    x1689 = x351 * x49
    x1690 = x1478 * x632 + x1689
    x1691 = x153 * x1690
    x1692 = x4 * (-x1678 * x60 + x1684)
    x1693 = 0.5 * x1692
    x1694 = x1680 + x390
    x1695 = x1694 * x49
    x1696 = x1690 * x3 + x1695
    x1697 = x4 * (-x1681 * x60 + x1685)
    x1698 = 0.5 * x1697
    x1699 = x1695 * x3 + x1696 * x3 + x1698
    x1700 = x1464 * x634 + x426
    x1701 = x1700 * x49
    x1702 = x1681 * x3 + x1701
    x1703 = x1464 * x641 + x422
    x1704 = x1703 * x49
    x1705 = x1685 * x3 + x1704
    x1706 = x4 * (-x1700 * x60 + x1703)
    x1707 = 0.5 * x1706
    x1708 = x1464 * x666 + x418
    x1709 = x4 * (x1464 * x670 - x1708 * x60 + x416)
    x1710 = 0.5 * x1709
    x1711 = x4 * (-x1684 * x60 + x1708)
    x1712 = 0.5 * x1711
    x1713 = x1232 + x1464 * x685
    x1714 = x157 * x1629 + x456
    x1715 = x153 * x1714
    x1716 = x1713 * x3 + x1715
    x1717 = x1237 + x1464 * x691
    x1718 = x1464 * x692 + x458
    x1719 = x153 * x1718
    x1720 = x1717 * x3 + x1719
    x1721 = x1241 + x1464 * x696
    x1722 = x1492 * x488 + x447
    x1723 = x153 * x1722
    x1724 = x4 * (-x1713 * x60 + x1717)
    x1725 = 0.5 * x1724
    x1726 = x1246 + x1464 * x683
    x1727 = x1726 * x49
    x1728 = x1722 * x3 + x1727
    x1729 = x4 * (-x1714 * x60 + x1718)
    x1730 = 0.5 * x1729
    x1731 = x1727 * x3 + x1728 * x3 + x1730
    x1732 = x1252 + x1464 * x689
    x1733 = x1732 * x49
    x1734 = x1714 * x3 + x1733
    x1735 = x1256 + x1464 * x713
    x1736 = x1735 * x49
    x1737 = x1718 * x3 + x1736
    x1738 = x4 * (-x1732 * x60 + x1735)
    x1739 = 0.5 * x1738
    x1740 = x1262 + x1464 * x715
    x1741 = x4 * (x1261 + x1464 * x719 - x1740 * x60)
    x1742 = 0.5 * x1741
    x1743 = x4 * (-x1717 * x60 + x1740)
    x1744 = 0.5 * x1743
    x1745 = x1464 * x736 + x501
    x1746 = x1745 * x3
    x1747 = x1464 * x732
    x1748 = x1747 * x5
    x1749 = x116 * x504 + x1748
    x1750 = x1749 * x49
    x1751 = 2.0 * x1750
    x1752 = x1746 + x1751
    x1753 = x1464 * x744 + x513
    x1754 = x1753 * x3
    x1755 = x116 * x539 + x1464 * x746
    x1756 = x1755 * x49
    x1757 = x1754 + 2.0 * x1756
    x1758 = x1464 * x751 + x523
    x1759 = x1758 * x3
    x1760 = x116 * x486 + x1464 * x753
    x1761 = x1760 * x49
    x1762 = 2.0 * x1761
    x1763 = x1745 * x60
    x1764 = x4 * (x1753 - x1763)
    x1765 = 0.5 * x1764
    x1766 = x1747 + 3.0 * x531
    x1767 = x1766 * x49
    x1768 = x1760 * x3
    x1769 = x1767 + x1768
    x1770 = x1749 * x60
    x1771 = x4 * (x1755 - x1770)
    x1772 = 0.5 * x1771
    x1773 = x1767 * x3
    x1774 = x1769 * x3 + x1772 + x1773
    x1775 = x1464 * x733 + 3.0 * x571
    x1776 = x1775 * x49
    x1777 = x1749 * x3 + x1776
    x1778 = x1464 * x741 + 3.0 * x567
    x1779 = x1778 * x49
    x1780 = x1755 * x3 + x1779
    x1781 = x4 * (-x1775 * x60 + x1778)
    x1782 = 0.5 * x1781
    x1783 = x1464 * x775 + x564
    x1784 = x4 * (x1464 * x780 - x1783 * x60 + x562)
    x1785 = 0.5 * x1784
    x1786 = x4 * (-x1753 * x60 + x1783)
    x1787 = 0.5 * x1786
    x1788 = x1464 * x4 * (x792 - x794)
    x1789 = 0.5 * x1788
    x1790 = x1464 * x797 + x1789
    x1791 = x1464 * x807
    x1792 = x1464 * x808 + x1791
    x1793 = x1464 * x810
    x1794 = x1464 * x812 + x1793
    x1795 = x1464 * x814
    x1796 = x4 * (x1464 * x5 * x793 - x1664 * x796)
    x1797 = 0.5 * x1796
    x1798 = x1464 * x4 * (x5 * x801 - x60 * x822)
    x1799 = 0.5 * x1798
    x1800 = x1464 * x4 * (x5 * x792 - x60 * x811)
    x1801 = 0.5 * x1800
    x1802 = x1464 * x831 + x607
    x1803 = x1464 * x832 + x604
    x1804 = x1803 * x60
    x1805 = x1464 * x834
    x1806 = x1805 + x582
    x1807 = x4 * (x1802 - x1804)
    x1808 = 0.5 * x1807
    x1809 = x1806 * x429 + x1808
    x1810 = x1803 * x49
    x1811 = x1805 * x5 + x583
    x1812 = x1810 + x1811 * x3
    x1813 = x1802 * x49
    x1814 = x1464 * x847 + x49 * x589
    x1815 = x1813 + x1814 * x3
    x1816 = x1806 * x49
    x1817 = x1478 * x850 + x49 * x595
    x1818 = x4 * (-x1811 * x60 + x1814)
    x1819 = 0.5 * x1818
    x1820 = x1464 * x854 + x624
    x1821 = x4 * (x1464 * x859 - x1820 * x60 + x49 * x626)
    x1822 = 0.5 * x1821
    x1823 = x4 * (-x1814 * x60 + x1820)
    x1824 = 0.5 * x1823
    x1825 = x1464 * x862 + x673
    x1826 = x1464 * x863 + 2.0 * x656
    x1827 = x1826 * x60
    x1828 = x1464 * x865
    x1829 = x1828 + x638
    x1830 = x4 * (x1825 - x1827)
    x1831 = 0.5 * x1830
    x1832 = x1829 * x429 + x1831
    x1833 = x1826 * x49
    x1834 = x1828 * x5 + x639
    x1835 = x1833 + x1834 * x3
    x1836 = x1825 * x49
    x1837 = x1464 * x878 + x645
    x1838 = x1836 + x1837 * x3
    x1839 = x1829 * x49
    x1840 = x1478 * x881 + x651
    x1841 = x4 * (-x1834 * x60 + x1837)
    x1842 = 0.5 * x1841
    x1843 = x1464 * x885 + x674
    x1844 = x4 * (x1464 * x890 - x1843 * x60 + x677)
    x1845 = 0.5 * x1844
    x1846 = x4 * (-x1837 * x60 + x1843)
    x1847 = 0.5 * x1846
    x1848 = x1464 * x899 + 3.0 * x707
    x1849 = x157 * x1747 + 3.0 * x705
    x1850 = x1849 * x60
    x1851 = x116 * x683 + x1464 * x896
    x1852 = x4 * (x1848 - x1850)
    x1853 = 0.5 * x1852
    x1854 = x1851 * x429 + x1853
    x1855 = x1849 * x49
    x1856 = x116 * x702 + x1478 * x896
    x1857 = x1855 + x1856 * x3
    x1858 = x1848 * x49
    x1859 = x116 * x692 + x157 * x1748
    x1860 = x1858 + x1859 * x3
    x1861 = x1851 * x49
    x1862 = x116 * x697 + x1492 * x748
    x1863 = x4 * (-x1856 * x60 + x1859)
    x1864 = 0.5 * x1863
    x1865 = x116 * x722 + x1464 * x913
    x1866 = x4 * (x116 * x725 + x1464 * x919 - x1865 * x60)
    x1867 = 0.5 * x1866
    x1868 = x4 * (-x1859 * x60 + x1865)
    x1869 = 0.5 * x1868
    x1870 = x1464 * x920 + 4.0 * x764
    x1871 = x1464 * x921 + 4.0 * x761
    x1872 = x1871 * x60
    x1873 = x1464 * x924
    x1874 = x1873 + 4.0 * x738
    x1875 = x1874 * x429
    x1876 = x4 * (x1870 - x1872)
    x1877 = 0.5 * x1876
    x1878 = x1875 + x1877
    x1879 = x1871 * x49
    x1880 = x1873 * x5 + 4.0 * x739
    x1881 = x1880 * x3
    x1882 = x1879 + x1881
    x1883 = x1870 * x49
    x1884 = x1464 * x940 + x746 * x98
    x1885 = x1884 * x3
    x1886 = x1883 + x1885
    x1887 = x1874 * x49
    x1888 = x1478 * x944 + x753 * x98
    x1889 = x1888 * x3
    x1890 = x1880 * x60
    x1891 = x4 * (x1884 - x1890)
    x1892 = 0.5 * x1891
    x1893 = x1887 * x3
    x1894 = x1464 * x951 + 4.0 * x783
    x1895 = x4 * (x1464 * x957 - x1894 * x60 + x785 * x98)
    x1896 = 0.5 * x1895
    x1897 = x4 * (-x1884 * x60 + x1894)
    x1898 = 0.5 * x1897
    x1899 = x1464 * x4 * (x960 - x962)
    x1900 = 0.5 * x1899
    x1901 = x1464 * x4 * (-x60 * x968 + x972)
    x1902 = 0.5 * x1901
    x1903 = x1464 * x4 * (-x60 * x960 + x968)
    x1904 = 0.5 * x1903
    x1905 = x1464 * x975 + x810
    x1906 = x1464 * x976 + x807
    x1907 = x1906 * x60
    x1908 = x1464 * x978 + x814
    x1909 = x4 * (x1905 - x1907)
    x1910 = 0.5 * x1909
    x1911 = x1464 * x980 + x829
    x1912 = x4 * (x1464 * x983 - x1911 * x60 + x827)
    x1913 = 0.5 * x1912
    x1914 = x4 * (-x1905 * x60 + x1911)
    x1915 = 0.5 * x1914
    x1916 = x1464 * x986 + 2.0 * x846
    x1917 = x1464 * x987 + 2.0 * x844
    x1918 = x1917 * x60
    x1919 = x1464 * x989 + 2.0 * x849
    x1920 = x4 * (x1916 - x1918)
    x1921 = 0.5 * x1920
    x1922 = x1464 * x991 + 2.0 * x860
    x1923 = x4 * (x1464 * x994 - x1922 * x60 + 2.0 * x858)
    x1924 = 0.5 * x1923
    x1925 = x4 * (-x1916 * x60 + x1922)
    x1926 = 0.5 * x1925
    x1927 = x1431 + x1464 * x997
    x1928 = x1433 + x1464 * x998
    x1929 = x1928 * x60
    x1930 = x1000 * x1464 + x1436
    x1931 = x4 * (x1927 - x1929)
    x1932 = 0.5 * x1931
    x1933 = x1002 * x1464 + x1440
    x1934 = x4 * (x1005 * x1464 + x1439 - x1933 * x60)
    x1935 = 0.5 * x1934
    x1936 = x4 * (-x1927 * x60 + x1933)
    x1937 = 0.5 * x1936
    x1938 = x1012 * x1464 + 4.0 * x906
    x1939 = x157 * x1873 + 4.0 * x904
    x1940 = x1939 * x60
    x1941 = x1011 * x1464 + x896 * x98
    x1942 = x4 * (x1938 - x1940)
    x1943 = 0.5 * x1942
    x1944 = x1015 * x1464 + 4.0 * x917
    x1945 = x4 * (x1457 * x1464 - x1944 * x60 + x916 * x98)
    x1946 = 0.5 * x1945
    x1947 = x4 * (-x1938 * x60 + x1944)
    x1948 = 0.5 * x1947
    x1949 = x1018 * x1464 + 5.0 * x939
    x1950 = x1019 * x1464 + 5.0 * x936
    x1951 = x1950 * x60
    x1952 = x1951 * x3
    x1953 = x1022 * x1464 + 5.0 * x943
    x1954 = x1953 * x429
    x1955 = x4 * (x1949 - x1951)
    x1956 = 0.5 * x1955
    x1957 = x1026 * x1464 + 5.0 * x958
    x1958 = x4 * (x1030 * x1464 - x1957 * x60 + 5.0 * x956)
    x1959 = 0.5 * x1958
    x1960 = x4 * (-x1949 * x60 + x1957)
    x1961 = 0.5 * x1960
    x1962 = x1033**2
    x1963 = x132 + x1962 * x63
    x1964 = x141 + x1962 * x71
    x1965 = x1964 * x49
    x1966 = x1962 * x47
    x1967 = x122 + x1966
    x1968 = x145 + x1962 * x56
    x1969 = x1968 * x49
    x1970 = x1962 * x80
    x1971 = x1970 + x86
    x1972 = x1962 * x42
    x1973 = x1972 + x93
    x1974 = x1973 * x49
    x1975 = x1963 - x1967 * x60
    x1976 = x1962 * x44
    x1977 = x109 + x1976
    x1978 = x1977 * x49
    x1979 = x1973 * x3 + 4.0 * x1978
    x1980 = x1964 - x1968 * x60
    x1981 = x49 * (x150 + x1962 * x40)
    x1982 = x1038 * x49
    x1983 = x1033 * x1067 + x1982 + x235
    x1984 = x1033 * x1069 + x1054 * x49 + x241
    x1985 = x1033 * x57
    x1986 = x1033 * x1062 + x1985 + x230
    x1987 = x101 * x1033 + x1033 * x1064 + x245
    x1988 = x1033 * x82
    x1989 = x1033 * x1072 + x196 + x1988
    x1990 = x1033 * x88
    x1991 = x1033 * x1074 + x1990 + x202
    x1992 = x1983 - x1986 * x60
    x1993 = x1033 * x95
    x1994 = x1033 * x1077 + x1993 + x217
    x1995 = x116 * x1994 + x1991 * x3
    x1996 = x1984 - x1987 * x60
    x1997 = x49 * x53
    x1998 = x153 * (x1033 * x1084 + x1057 * x1997 + x251)
    x1999 = x1962 * x278 + x332
    x2000 = x1962 * x281 + x338
    x2001 = x2000 * x49
    x2002 = x1962 * x268 + x326
    x2003 = x1962 * x270 + x342
    x2004 = x2003 * x49
    x2005 = x1962 * x287 + x292
    x2006 = x1962 * x259 + x298
    x2007 = x2006 * x49
    x2008 = x1999 - x2002 * x60
    x2009 = x1962 * x263 + x314
    x2010 = x2009 * x49
    x2011 = x2006 * x3 + 3.0 * x2010
    x2012 = x2000 - x2003 * x60
    x2013 = x1962 * x5
    x2014 = x49 * (x2013 * x317 + x347)
    x2015 = x1033 * x1122 + x1069 * x153 + x412
    x2016 = x1033 * x1124 + x1090 * x153 + x421
    x2017 = x116 * x2016
    x2018 = x1033 * x1115 + x1064 * x153 + x406
    x2019 = x1033 * x1118 + x1088 * x153 + x425
    x2020 = x116 * x2019
    x2021 = x1033 * x1129 + x1074 * x153 + x384
    x2022 = x1033 * x1132 + x1077 * x153 + x389
    x2023 = x116 * x2022
    x2024 = x2015 - x2018 * x60
    x2025 = x1033 * x1136 + x1085 + x401
    x2026 = x153 * x2025 + x2022 * x3
    x2027 = x2016 - x2019 * x60
    x2028 = x49 * (x1033 * x1141 + 2.0 * x1093 + x430)
    x2029 = 64.3624305720022 * x154
    x2030 = x1033 * x1157 + x1099 * x49 + x469
    x2031 = x1033 * x1158 + x1107 * x49 + x474
    x2032 = x1033 * x1154 + x1033 * x271 + x464
    x2033 = x1033 * x1155 + x1033 * x307 + x477
    x2034 = x1033 * x1160 + x1033 * x289 + x446
    x2035 = x1033 * x1161 + x1033 * x294 + x449
    x2036 = x2030 - x2032 * x60
    x2037 = x1033 * x1163 + x1057 * x317 * x49 + x460
    x2038 = x153 * x2037
    x2039 = x2035 * x3 + x2038
    x2040 = x2031 - x2033 * x60
    x2041 = x49 * (x1033 * x1167 + x1110 + x481)
    x2042 = 111.478999849332 * x154
    x2043 = x1962 * x505 + x556
    x2044 = x1962 * x511 + x566
    x2045 = x2044 * x49
    x2046 = x1962 * x493 + x547
    x2047 = x1962 * x499 + x570
    x2048 = x2047 * x49
    x2049 = x1962 * x520 + x525
    x2050 = x1962 * x492 + x530
    x2051 = x2050 * x49
    x2052 = x2043 - x2046 * x60
    x2053 = x1962 * x486 + x542
    x2054 = x2053 * x49
    x2055 = x2050 * x3 + 2.0 * x2054
    x2056 = x2044 - x2047 * x60
    x2057 = x49 * (x1962 * x484 + x575)
    x2058 = x1033 * x1205 + x1125 + x622
    x2059 = x1033 * x1207 + x1147 * x116 + x627
    x2060 = x1033 * x1198 + x1119 + x618
    x2061 = x1033 * x1202 + x1145 * x116 + x629
    x2062 = x1033 * x1210 + x1133 + x597
    x2063 = x1033 * x1213 + x1136 * x116 + x601
    x2064 = x2058 - x2060 * x60
    x2065 = x1033 * x1216 + 3.0 * x1142 + x612
    x2066 = x2065 * x49
    x2067 = x2063 * x3 + x2066
    x2068 = x2059 - x2061 * x60
    x2069 = 64.3624305720022 * x154
    x2070 = x1033 * x1238 + x1158 * x153 + x671
    x2071 = x1033 * x1239 + x1175 + x678
    x2072 = x1033 * x1233 + x1155 * x153 + x667
    x2073 = x1033 * x1235 + x1172 + x680
    x2074 = x1033 * x1242 + x1161 * x153 + x652
    x2075 = x1033 * x1243 + x1164 + x654
    x2076 = x2070 - x2072 * x60
    x2077 = x1033 * x1247 + 2.0 * x1168 + x661
    x2078 = x2077 * x49
    x2079 = x2075 * x3 + x2078
    x2080 = x2071 - x2073 * x60
    x2081 = 143.918769956108 * x154
    x2082 = x1033 * x1272 + x1183 * x49 + x720
    x2083 = x1033 * x1273 + x1192 * x49 + x727
    x2084 = x153 * x2083
    x2085 = x1033 * x1267 + x1033 * x500 + x716
    x2086 = x1033 * x1268 + x1033 * x1269 + x729
    x2087 = x153 * x2086
    x2088 = x1033 * x1276 + x1033 * x522 + x699
    x2089 = x1033 * x1277 + x1033 * x1278 + x703
    x2090 = x153 * x2089
    x2091 = x2082 - x2085 * x60
    x2092 = x1033 * x1281 + x1188 + x710
    x2093 = x2092 * x49
    x2094 = x2089 * x3 + x2093
    x2095 = x2083 - x2086 * x60
    x2096 = x1962 * x744 + x781
    x2097 = x1962 * x746 + x787
    x2098 = x2097 * x49
    x2099 = x1962 * x736 + x776
    x2100 = x1962 * x732
    x2101 = x2100 * x5 + x789
    x2102 = x2101 * x49
    x2103 = x1962 * x751 + x755
    x2104 = x1962 * x753 + x758
    x2105 = x2104 * x49
    x2106 = x2096 - x2099 * x60
    x2107 = x2100 + x769
    x2108 = x2107 * x49
    x2109 = x2104 * x3 + x2108
    x2110 = x2097 - x2101 * x60
    x2111 = x1033 * x1316 + 4.0 * x1225 + x805
    x2112 = x2111 * x49
    x2113 = x1033 * x1332 + x1207 * x98 + x823
    x2114 = x1033 * x1319 + 4.0 * x1222 + x802
    x2115 = x2114 * x49
    x2116 = x1033 * x1328 + x1202 * x98 + x819
    x2117 = x1033 * x1323 + 4.0 * x1217 + x798
    x2118 = x2111 - x2114 * x60
    x2119 = x2117 * x49
    x2120 = x1033 * x1336 + x1213 * x98 + x817
    x2121 = x2113 - x2116 * x60
    x2122 = x1033 * x1342 + 3.0 * x1258 + x842
    x2123 = x2122 * x49
    x2124 = x1033 * x1354 + x116 * x1239 + x855
    x2125 = x1033 * x1343 + 3.0 * x1254 + x839
    x2126 = x2125 * x49
    x2127 = x1033 * x1351 + x116 * x1235 + x852
    x2128 = x1033 * x1347 + 3.0 * x1248 + x835
    x2129 = x2122 - x2125 * x60
    x2130 = x2128 * x49
    x2131 = x1033 * x1357 + x116 * x1243 + x851
    x2132 = x2124 - x2127 * x60
    x2133 = x1033 * x1362 + 2.0 * x1290 + x873
    x2134 = x2133 * x49
    x2135 = x1033 * x1373 + x1274 + x886
    x2136 = x1033 * x1363 + 2.0 * x1287 + x870
    x2137 = x2136 * x49
    x2138 = x1033 * x1370 + x1270 + x883
    x2139 = x1033 * x1366 + 2.0 * x1282 + x866
    x2140 = x2133 - x2136 * x60
    x2141 = x2139 * x49
    x2142 = x1033 * x1376 + x1279 + x882
    x2143 = x2135 - x2138 * x60
    x2144 = x1033 * x1381 + x1309 + x902
    x2145 = x2144 * x49
    x2146 = x1033 * x1391 + x1298 * x49 + x914
    x2147 = x1033 * x1382 + x1307 + x900
    x2148 = x2147 * x49
    x2149 = x1033 * x1388 + x1296 + x911
    x2150 = x1033 * x1384 + x1302 + x895
    x2151 = x2144 - x2147 * x60
    x2152 = x2150 * x49
    x2153 = x1033 * x1394 + x1300 * x49 + x909
    x2154 = x2146 - x2149 * x60
    x2155 = x1962 * x920 + x934
    x2156 = x2155 * x49
    x2157 = x1962 * x940 + x952
    x2158 = x1962 * x921 + x931
    x2159 = x2158 * x49
    x2160 = x1962 * x924
    x2161 = x2160 * x5 + x948
    x2162 = x2160 + x926
    x2163 = x2155 - x2158 * x60
    x2164 = x2162 * x49
    x2165 = x2013 * x944 + x946
    x2166 = x2157 - x2161 * x60
    x2167 = x1033 * x1411 + 5.0 * x1330 + x973
    x2168 = x1033 * x1414 + 5.0 * x1326 + x969
    x2169 = x2168 * x60
    x2170 = x1033 * x1418 + 5.0 * x1334 + x966
    x2171 = x2167 - x2169
    x2172 = x1033 * x1423 + 4.0 * x1353 + x984
    x2173 = x1033 * x1424 + 4.0 * x1350 + x981
    x2174 = x2173 * x60
    x2175 = x1033 * x1426 + 4.0 * x1356 + x979
    x2176 = x2172 - x2174
    x2177 = x1033 * x1432 + 3.0 * x1372 + x995
    x2178 = x1033 * x1434 + 3.0 * x1369 + x992
    x2179 = x2178 * x60
    x2180 = x1033 * x1437 + 3.0 * x1375 + x990
    x2181 = x2177 - x2179
    x2182 = x1006 + x1033 * x1444 + 2.0 * x1390
    x2183 = x1003 + x1033 * x1445 + 2.0 * x1387
    x2184 = x2183 * x60
    x2185 = x1001 + x1033 * x1447 + 2.0 * x1393
    x2186 = x2182 - x2184
    x2187 = x1016 + x1033 * x1452 + x1403
    x2188 = x1013 + x1033 * x1453 + x1401
    x2189 = x2188 * x60
    x2190 = x1010 + x1033 * x1455 + x1405
    x2191 = x2187 - x2189
    x2192 = x1018 * x1962 + x1031
    x2193 = x1019 * x1962 + x1027
    x2194 = x2193 * x60
    x2195 = x1022 * x1962 + x1024
    x2196 = x2192 - x2194
    x2197 = x1464 * x4 * (x1033 * x63 - x1042)
    x2198 = x1464 * (x1044 + x1045)
    x2199 = x1464 * x4 * (x1033 * x71 - x1047)
    x2200 = x1464 * x72
    x2201 = x1066 * x1464 + x2200
    x2202 = x105 * x1464
    x2203 = x1068 * x1464 + x2202
    x2204 = x1464 * x57
    x2205 = x1061 * x1464 + x2204
    x2206 = x101 * x1464
    x2207 = x1063 * x1464 + x2206
    x2208 = x1464 * x82
    x2209 = x1071 * x1464 + x2208
    x2210 = x1464 * x88
    x2211 = x1073 * x1464 + x2210
    x2212 = x4 * (x2201 - x2205 * x60)
    x2213 = x1464 * x95
    x2214 = x1076 * x1464 + x2213
    x2215 = x116 * x2214 + x2211 * x3
    x2216 = x4 * (x2203 - x2207 * x60)
    x2217 = x1478 * x1997
    x2218 = x153 * (x1083 * x1464 + x2217)
    x2219 = x1033 * x1514
    x2220 = x1033 * x4 * (x1512 - x1522)
    x2221 = x1033 * (x1525 + x1528)
    x2222 = x1033 * x4 * (x1514 - x1530)
    x2223 = x1121 * x1464 + x1485 * x153
    x2224 = x1123 * x1464 + x1497 * x153
    x2225 = x116 * x2224
    x2226 = x1464 * (x1113 + x1114)
    x2227 = x1464 * (x1116 + x1117)
    x2228 = x116 * x2227
    x2229 = x1464 * (x1127 + x1128)
    x2230 = x1464 * (x1130 + x1131)
    x2231 = x116 * x2230
    x2232 = x4 * (x2223 - x2226 * x60)
    x2233 = x1135 * x1464 + x1493
    x2234 = x153 * x2233 + x2230 * x3
    x2235 = x4 * (x2224 - x2227 * x60)
    x2236 = x1464 * x49 * (x1139 + x1140)
    x2237 = x1033 * x1579 + x1515
    x2238 = x1033 * x1580 + x1543
    x2239 = x1033 * x1576 + x1509
    x2240 = x1033 * x1577 + x1540
    x2241 = x1033 * x1582 + x1520
    x2242 = x1033 * x1583 + x1527
    x2243 = x4 * (x2237 - x2239 * x60)
    x2244 = x1033 * x1586 + x1535
    x2245 = x153 * x2244
    x2246 = x2242 * x3 + x2245
    x2247 = x4 * (x2238 - x2240 * x60)
    x2248 = x49 * (x1033 * x1591 + x1547)
    x2249 = 193.087291716007 * x154
    x2250 = x1033 * x1616
    x2251 = x1033 * x4 * (x1614 - x1625)
    x2252 = x1033 * (x1628 + x1632)
    x2253 = x1033 * x4 * (x1616 - x1634)
    x2254 = x1033 * x1638
    x2255 = x1204 * x1464 + x1556
    x2256 = x116 * x1568 + x1206 * x1464
    x2257 = x1197 * x1464 + x1553
    x2258 = x1464 * (x1200 + x1201)
    x2259 = x1209 * x1464 + x1558
    x2260 = x1464 * (x1211 + x1212)
    x2261 = x4 * (x2255 - x2257 * x60)
    x2262 = x1464 * (x1199 + x1215)
    x2263 = x2262 * x49
    x2264 = x2260 * x3 + x2263
    x2265 = x4 * (x2256 - x2258 * x60)
    x2266 = 111.478999849332 * x154
    x2267 = x153 * x1580
    x2268 = x1033 * x1684 + x2267
    x2269 = x1033 * x1685 + x1599
    x2270 = x153 * x1577
    x2271 = x1033 * x1678 + x2270
    x2272 = x1033 * x1681 + x1596
    x2273 = x153 * x1583
    x2274 = x1033 * x1688 + x2273
    x2275 = x1033 * x1690 + x1587
    x2276 = x4 * (x2268 - x2271 * x60)
    x2277 = 2.0 * x1592
    x2278 = x1033 * x1694 + x2277
    x2279 = x2278 * x49
    x2280 = x2275 * x3 + x2279
    x2281 = x4 * (x2269 - x2272 * x60)
    x2282 = 249.274621726796 * x154
    x2283 = x1033 * x1717 + x1617
    x2284 = x1033 * x1718 + x1646
    x2285 = x153 * x2284
    x2286 = x1033 * x1713 + x1611
    x2287 = x1033 * x1714 + x1643
    x2288 = x153 * x2287
    x2289 = x1033 * x1721 + x1623
    x2290 = x1033 * x1722 + x1631
    x2291 = x153 * x2290
    x2292 = x4 * (x2283 - x2286 * x60)
    x2293 = x1033 * x1726 + x1638
    x2294 = x2293 * x49
    x2295 = x2290 * x3 + x2294
    x2296 = x4 * (x2284 - x2287 * x60)
    x2297 = x1033 * x1755
    x2298 = x1033 * x4 * (x1753 - x1763)
    x2299 = x1033 * x1767
    x2300 = x1033 * x1768 + x2299
    x2301 = x1033 * x4 * (x1755 - x1770)
    x2302 = x1464 * (x1314 + x1315)
    x2303 = x2302 * x49
    x2304 = x1331 * x1464 + x1658 * x98
    x2305 = x1464 * (x1317 + x1318)
    x2306 = x2305 * x49
    x2307 = x1322 * x1478 + x1327 * x1464
    x2308 = x1464 * (x1321 + x1322)
    x2309 = x4 * (x2302 - x2305 * x60)
    x2310 = x2308 * x49
    x2311 = x1335 * x1464 + x1660 * x98
    x2312 = x4 * (x2304 - x2307 * x60)
    x2313 = x1033 * x1802 + 3.0 * x1704
    x2314 = x2313 * x49
    x2315 = x1033 * x1814 + x116 * x1685
    x2316 = x1033 * x1803 + 3.0 * x1701
    x2317 = x2316 * x49
    x2318 = x1033 * x1811 + x116 * x1681
    x2319 = x1033 * x1806 + 3.0 * x1695
    x2320 = x4 * (x2313 - x2316 * x60)
    x2321 = x2319 * x49
    x2322 = x1033 * x1817 + x116 * x1690
    x2323 = x4 * (x2315 - x2318 * x60)
    x2324 = x1033 * x1825 + 2.0 * x1736
    x2325 = x2324 * x49
    x2326 = x1033 * x1837 + x1719
    x2327 = x1033 * x1826 + 2.0 * x1733
    x2328 = x2327 * x49
    x2329 = x1033 * x1834 + x1715
    x2330 = x1033 * x1829 + 2.0 * x1727
    x2331 = x4 * (x2324 - x2327 * x60)
    x2332 = x2330 * x49
    x2333 = x1033 * x1840 + x1723
    x2334 = x4 * (x2326 - x2329 * x60)
    x2335 = x1033 * x1848 + x1779
    x2336 = x2335 * x49
    x2337 = x1033 * x1859 + x1756
    x2338 = x1033 * x1849 + x1776
    x2339 = x2338 * x49
    x2340 = x1033 * x1856 + x1750
    x2341 = x1033 * x1851 + x1767
    x2342 = x4 * (x2335 - x2338 * x60)
    x2343 = x2341 * x49
    x2344 = x1033 * x1862 + x1761
    x2345 = x4 * (x2337 - x2340 * x60)
    x2346 = x1033 * x1883
    x2347 = x1033 * x1879
    x2348 = x1033 * x4 * (x1870 - x1872)
    x2349 = x1033 * x1887
    x2350 = x1033 * x4 * (x1884 - x1890)
    x2351 = x1464 * (x1409 + x1410)
    x2352 = x1464 * (x1412 + x1413)
    x2353 = x2352 * x60
    x2354 = x1464 * (x1416 + x1417)
    x2355 = x4 * (x2351 - x2353)
    x2356 = x1033 * x1905 + 4.0 * x1813
    x2357 = x1033 * x1906 + 4.0 * x1810
    x2358 = x2357 * x60
    x2359 = x1033 * x1908 + 4.0 * x1816
    x2360 = x4 * (x2356 - x2358)
    x2361 = 3.0 * x1836
    x2362 = x1033 * x1916 + x2361
    x2363 = 3.0 * x1833
    x2364 = x1033 * x1917 + x2363
    x2365 = x2364 * x60
    x2366 = 3.0 * x1839
    x2367 = x1033 * x1919 + x2366
    x2368 = x4 * (x2362 - x2365)
    x2369 = x1033 * x1927 + 2.0 * x1858
    x2370 = x1033 * x1928 + 2.0 * x1855
    x2371 = x2370 * x60
    x2372 = x1033 * x1930 + 2.0 * x1861
    x2373 = x4 * (x2369 - x2371)
    x2374 = x1033 * x1938 + x1883
    x2375 = x1033 * x1939 + x1879
    x2376 = x2375 * x60
    x2377 = x1033 * x1941 + x1887
    x2378 = x4 * (x2374 - x2376)
    x2379 = x1033 * x4 * (x1949 - x1951)
    x2380 = x1464**2
    x2381 = x132 + x2380 * x63
    x2382 = x141 + x2380 * x71
    x2383 = x2382 * x49
    x2384 = x122 + x2380 * x47
    x2385 = x145 + x2380 * x56
    x2386 = x2385 * x49
    x2387 = x2380 * x80 + x86
    x2388 = x2387 * x3
    x2389 = x2380 * x42 + x93
    x2390 = x2389 * x49
    x2391 = 5.0 * x2390
    x2392 = x2384 * x60
    x2393 = x2381 - x2392
    x2394 = x2393 * x85
    x2395 = x2389 * x3
    x2396 = x109 + x2380 * x44
    x2397 = x2396 * x49
    x2398 = 4.0 * x2397
    x2399 = x2395 + x2398
    x2400 = x2382 - x2385 * x60
    x2401 = x2400 * x85
    x2402 = x49 * (x150 + x2380 * x40)
    x2403 = x182 * x2380 + x235
    x2404 = x185 * x2380 + x241
    x2405 = x2404 * x49
    x2406 = x172 * x2380 + x230
    x2407 = x174 * x2380 + x245
    x2408 = x2407 * x49
    x2409 = x191 * x2380 + x196
    x2410 = x163 * x2380 + x202
    x2411 = x2410 * x49
    x2412 = x2403 - x2406 * x60
    x2413 = x2412 * x85
    x2414 = x167 * x2380 + x217
    x2415 = x2414 * x49
    x2416 = x2410 * x3 + 3.0 * x2415
    x2417 = x2404 - x2407 * x60
    x2418 = x2417 * x85
    x2419 = x2380 * x5
    x2420 = x49 * (x2419 * x450 + x251)
    x2421 = 2.0 * x2420
    x2422 = x1464 * x1512 + x2200 + x332
    x2423 = x1464 * x1514 + x2202 + x338
    x2424 = x2423 * x49
    x2425 = x1464 * x1506 + x2204 + x326
    x2426 = x1464 * x1508 + x2206 + x342
    x2427 = x2426 * x49
    x2428 = x1464 * x1517 + x2208 + x292
    x2429 = x2428 * x3
    x2430 = x1464 * x1519 + x2210 + x298
    x2431 = x2430 * x49
    x2432 = 4.0 * x2431
    x2433 = x2425 * x60
    x2434 = x2422 - x2433
    x2435 = x2434 * x85
    x2436 = x2430 * x3
    x2437 = x1464 * x1526 + x2213 + x314
    x2438 = x2437 * x49
    x2439 = 3.0 * x2438
    x2440 = x2436 + x2439
    x2441 = x2423 - x2426 * x60
    x2442 = x2441 * x85
    x2443 = x49 * (x1464 * x1534 + x2217 + x347)
    x2444 = 2.0 * x2443
    x2445 = x2380 * x367 + x412
    x2446 = x2380 * x371 + x421
    x2447 = x2446 * x49
    x2448 = 3.0 * x2447
    x2449 = x2380 * x357 + x406
    x2450 = x2380 * x361 + x425
    x2451 = x2450 * x49
    x2452 = 3.0 * x2451
    x2453 = x2380 * x379 + x384
    x2454 = x2380 * x356 + x389
    x2455 = x2454 * x49
    x2456 = 3.0 * x2455
    x2457 = x2445 - x2449 * x60
    x2458 = x2457 * x85
    x2459 = x2380 * x351 + x401
    x2460 = x2459 * x49
    x2461 = x2454 * x3 + 2.0 * x2460
    x2462 = x2446 - x2450 * x60
    x2463 = x2462 * x85
    x2464 = x49 * (x2380 * x349 + x430)
    x2465 = x1464 * x1579 + x1485 * x49 + x469
    x2466 = x1464 * x1580 + x1497 * x49 + x474
    x2467 = x1464 * x1576 + x1464 * x175 + x464
    x2468 = x1464 * x1577 + x1464 * x210 + x477
    x2469 = x1464 * x1582 + x1464 * x193 + x446
    x2470 = x1464 * x1583 + x1464 * x198 + x449
    x2471 = x2465 - x2467 * x60
    x2472 = x2471 * x85
    x2473 = x1464 * x1586 + x1492 * x38 * x49 + x460
    x2474 = x153 * x2473
    x2475 = x2470 * x3 + x2474
    x2476 = x2466 - x2468 * x60
    x2477 = x2476 * x85
    x2478 = x49 * (x1464 * x1591 + x1501 + x481)
    x2479 = x1464 * x1614 + 2.0 * x1515 + x556
    x2480 = x1464 * x1616 + 2.0 * x1543 + x566
    x2481 = x2480 * x49
    x2482 = 3.0 * x2481
    x2483 = x1464 * x1608 + 2.0 * x1509 + x547
    x2484 = x1464 * x1610 + 2.0 * x1540 + x570
    x2485 = x2484 * x49
    x2486 = 3.0 * x2485
    x2487 = x1464 * x1620 + 2.0 * x1520 + x525
    x2488 = x2487 * x3
    x2489 = x1464 * x1622 + 2.0 * x1527 + x530
    x2490 = x2489 * x49
    x2491 = 3.0 * x2490
    x2492 = x2483 * x60
    x2493 = x2479 - x2492
    x2494 = x2493 * x85
    x2495 = x2489 * x3
    x2496 = x1464 * x1630 + x1536 + x542
    x2497 = x2496 * x49
    x2498 = 2.0 * x2497
    x2499 = x2495 + x2498
    x2500 = x2480 - x2484 * x60
    x2501 = x2500 * x85
    x2502 = x49 * (x1464 * x1637 + 2.0 * x1547 + x575)
    x2503 = x2380 * x587 + x622
    x2504 = x2380 * x589 + x627
    x2505 = x2504 * x49
    x2506 = x2380 * x580 + x618
    x2507 = x2380 * x577
    x2508 = x2507 * x5 + x629
    x2509 = x2508 * x49
    x2510 = x2380 * x593 + x597
    x2511 = x2380 * x595 + x601
    x2512 = x2511 * x49
    x2513 = x2503 - x2506 * x60
    x2514 = x2513 * x85
    x2515 = x2507 + x612
    x2516 = x2515 * x49
    x2517 = x2511 * x3 + x2516
    x2518 = x2504 - x2508 * x60
    x2519 = x2518 * x85
    x2520 = x1464 * x1684 + x1555 * x49 + x671
    x2521 = x1464 * x1685 + x1568 * x49 + x678
    x2522 = x153 * x2521
    x2523 = x1464 * x1678 + x1464 * x362 + x667
    x2524 = x1464 * x1679 + x1464 * x1681 + x680
    x2525 = x153 * x2524
    x2526 = x1464 * x1688 + x1464 * x381 + x652
    x2527 = x1464 * x1689 + x1464 * x1690 + x654
    x2528 = x153 * x2527
    x2529 = x2520 - x2523 * x60
    x2530 = x2529 * x85
    x2531 = x1464 * x1694 + x1564 + x661
    x2532 = x2531 * x49
    x2533 = x2527 * x3 + x2532
    x2534 = x2521 - x2524 * x60
    x2535 = x2534 * x85
    x2536 = x1464 * x1717 + x2267 + x720
    x2537 = x1464 * x1718 + x1599 + x727
    x2538 = x153 * x2537
    x2539 = x1464 * x1713 + x2270 + x716
    x2540 = x1464 * x1714 + x1596 + x729
    x2541 = x153 * x2540
    x2542 = x1464 * x1721 + x2273 + x699
    x2543 = x1464 * x1722 + x1587 + x703
    x2544 = x153 * x2543
    x2545 = x2536 - x2539 * x60
    x2546 = x2545 * x85
    x2547 = x1464 * x1726 + x2277 + x710
    x2548 = x2547 * x49
    x2549 = x2543 * x3 + x2548
    x2550 = x2537 - x2540 * x60
    x2551 = x2550 * x85
    x2552 = x1464 * x1753 + x1618 + x781
    x2553 = x1464 * x1755 + 3.0 * x1646 + x787
    x2554 = x2553 * x49
    x2555 = x1464 * x1745 + x1612 + x776
    x2556 = x1464 * x1749 + 3.0 * x1643 + x789
    x2557 = x2556 * x49
    x2558 = x1464 * x1758 + x1624 + x755
    x2559 = x2558 * x3
    x2560 = x1464 * x1760 + 3.0 * x1631 + x758
    x2561 = x2560 * x49
    x2562 = 2.0 * x2561
    x2563 = x2555 * x60
    x2564 = x2552 - x2563
    x2565 = x2564 * x85
    x2566 = x1464 * x1766 + 3.0 * x1638 + x769
    x2567 = x2566 * x49
    x2568 = x2560 * x3
    x2569 = x2567 + x2568
    x2570 = x2553 - x2556 * x60
    x2571 = x2570 * x85
    x2572 = x2380 * x792 + x805
    x2573 = x2572 * x49
    x2574 = x2380 * x811 + x823
    x2575 = x2380 * x793 + x802
    x2576 = x2575 * x49
    x2577 = x2380 * x796
    x2578 = x2577 * x5 + x819
    x2579 = x2577 + x798
    x2580 = x2572 - x2575 * x60
    x2581 = x2580 * x85
    x2582 = x2579 * x49
    x2583 = x2419 * x815 + x817
    x2584 = x2574 - x2578 * x60
    x2585 = x2584 * x85
    x2586 = x1464 * x1802 + x1670 + x842
    x2587 = x2586 * x49
    x2588 = x1464 * x1814 + x1658 * x49 + x855
    x2589 = x1464 * x1803 + x1668 + x839
    x2590 = x2589 * x49
    x2591 = x1464 * x1811 + x1656 + x852
    x2592 = x1464 * x1806 + x1655 + x835
    x2593 = x2586 - x2589 * x60
    x2594 = x2593 * x85
    x2595 = x2592 * x49
    x2596 = x1464 * x1817 + x1660 * x49 + x851
    x2597 = x2588 - x2591 * x60
    x2598 = x2597 * x85
    x2599 = x1464 * x1825 + 2.0 * x1704 + x873
    x2600 = x2599 * x49
    x2601 = x1464 * x1837 + x1686 + x886
    x2602 = x1464 * x1826 + 2.0 * x1701 + x870
    x2603 = x2602 * x49
    x2604 = x1464 * x1834 + x1682 + x883
    x2605 = x1464 * x1829 + 2.0 * x1695 + x866
    x2606 = x2599 - x2602 * x60
    x2607 = x2606 * x85
    x2608 = x2605 * x49
    x2609 = x1464 * x1840 + x1691 + x882
    x2610 = x2601 - x2604 * x60
    x2611 = x2610 * x85
    x2612 = x1464 * x1848 + 3.0 * x1736 + x902
    x2613 = x2612 * x49
    x2614 = x116 * x1718 + x1464 * x1859 + x914
    x2615 = x1464 * x1849 + 3.0 * x1733 + x900
    x2616 = x2615 * x49
    x2617 = x116 * x1714 + x1464 * x1856 + x911
    x2618 = x1464 * x1851 + 3.0 * x1727 + x895
    x2619 = x2612 - x2615 * x60
    x2620 = x2619 * x85
    x2621 = x2618 * x49
    x2622 = x116 * x1722 + x1464 * x1862 + x909
    x2623 = x2614 - x2617 * x60
    x2624 = x2623 * x85
    x2625 = x1464 * x1870 + 4.0 * x1779 + x934
    x2626 = x2625 * x49
    x2627 = x1464 * x1884 + 4.0 * x1756 + x952
    x2628 = x1464 * x1871 + 4.0 * x1776 + x931
    x2629 = x2628 * x49
    x2630 = x1464 * x1880 + 4.0 * x1750 + x948
    x2631 = x1464 * x1874 + 4.0 * x1767 + x926
    x2632 = x2625 - x2628 * x60
    x2633 = x2632 * x85
    x2634 = x2631 * x49
    x2635 = x1464 * x1888 + 4.0 * x1761 + x946
    x2636 = x2635 * x3
    x2637 = x2630 * x60
    x2638 = x2627 - x2637
    x2639 = x2638 * x85
    x2640 = x2634 * x3
    x2641 = x2380 * x960 + x973
    x2642 = x2380 * x961 + x969
    x2643 = x2642 * x60
    x2644 = x2380 * x964 + x966
    x2645 = x2641 - x2643
    x2646 = x2645 * x85
    x2647 = x1464 * x1905 + x1793 + x984
    x2648 = x1464 * x1906 + x1791 + x981
    x2649 = x2648 * x60
    x2650 = x1464 * x1908 + x1795 + x979
    x2651 = x2647 - x2649
    x2652 = x2651 * x85
    x2653 = x1464 * x1916 + 2.0 * x1813 + x995
    x2654 = x1464 * x1917 + 2.0 * x1810 + x992
    x2655 = x2654 * x60
    x2656 = x1464 * x1919 + 2.0 * x1816 + x990
    x2657 = x2653 - x2655
    x2658 = x2657 * x85
    x2659 = x1006 + x1464 * x1927 + x2361
    x2660 = x1003 + x1464 * x1928 + x2363
    x2661 = x2660 * x60
    x2662 = x1001 + x1464 * x1930 + x2366
    x2663 = x2659 - x2661
    x2664 = x2663 * x85
    x2665 = x1016 + x1464 * x1938 + 4.0 * x1858
    x2666 = x1013 + x1464 * x1939 + 4.0 * x1855
    x2667 = x2666 * x60
    x2668 = x1010 + x1464 * x1941 + 4.0 * x1861
    x2669 = x2665 - x2667
    x2670 = x2669 * x85
    x2671 = x1031 + x1464 * x1949 + 5.0 * x1883
    x2672 = x1027 + x1464 * x1950 + 5.0 * x1879
    x2673 = x2672 * x60
    x2674 = x1024 + x1464 * x1953 + 5.0 * x1887
    x2675 = x2674 * x429
    x2676 = x2671 - x2673
    x2677 = x2676 * x85
    x2678 = x1033 * x1971 + x1043
    x2679 = x1033 * x1973 + x1048
    x2680 = x2679 * x49
    x2681 = x1033 * x1963 + x1059 - x60 * (x1033 * x1967 + x1060)
    x2682 = x49 * (x1033 * x1977 + x1056)
    x2683 = x1033 * x1989 + x1075 + x1974
    x2684 = x1033 * x1991 + x1079 + x1978
    x2685 = x1033 * x1983 + x1095 + x1965 - x60 * (x1033 * x1986 + x1096 + x1969)
    x2686 = x1033 * x1994 + x1092 + x1981
    x2687 = x1033 * x2005 + x1101
    x2688 = x1033 * x2006 + x1103
    x2689 = x2688 * x49
    x2690 = x1033 * x1999 + x1111 - x60 * (x1033 * x2002 + x1112)
    x2691 = x49 * (x1033 * x2009 + x1109)
    x2692 = x1033 * x2021 + x1134 + x153 * x1991
    x2693 = x1033 * x2022 + x1138 + x153 * x1994
    x2694 = x116 * x2693
    x2695 = (
        x1033 * x2015
        + x1151
        + x153 * x1984
        - x60 * (x1033 * x2018 + x1152 + x153 * x1987)
    )
    x2696 = x1033 * x2025 + x1149 + x1998
    x2697 = x1033 * x2034 + x1162 + x2007
    x2698 = x1033 * x2035 + x1166 + x2010
    x2699 = x1033 * x2030 + x1179 + x2001 - x60 * (x1033 * x2032 + x1180 + x2004)
    x2700 = x153 * (x1033 * x2037 + x1177 + x2014)
    x2701 = x1033 * x2049 + x1185
    x2702 = x1033 * x2050 + x1187
    x2703 = x2702 * x49
    x2704 = x1033 * x2043 + x1195 - x60 * (x1033 * x2046 + x1196)
    x2705 = x49 * (x1033 * x2053 + x1194)
    x2706 = x1033 * x2062 + x1214 + x2023
    x2707 = x1033 * x2063 + x116 * x2025 + x1219
    x2708 = x1033 * x2058 + x1229 + x2017 - x60 * (x1033 * x2060 + x1230 + x2020)
    x2709 = x49 * (x1033 * x2065 + x1227 + 3.0 * x2028)
    x2710 = 64.3624305720022 * x154
    x2711 = x1033 * x2074 + x1244 + x153 * x2035
    x2712 = x1033 * x2075 + x1250 + x2038
    x2713 = (
        x1033 * x2070
        + x1264
        + x153 * x2031
        - x60 * (x1033 * x2072 + x1265 + x153 * x2033)
    )
    x2714 = x49 * (x1033 * x2077 + x1260 + 2.0 * x2041)
    x2715 = 143.918769956108 * x154
    x2716 = x1033 * x2088 + x1280 + x2051
    x2717 = x1033 * x2089 + x1284 + x2054
    x2718 = x153 * x2717
    x2719 = x1033 * x2082 + x1294 + x2045 - x60 * (x1033 * x2085 + x1295 + x2048)
    x2720 = x49 * (x1033 * x2092 + x1292 + x2057)
    x2721 = x1033 * x2103 + x1301
    x2722 = x1033 * x2104 + x1305
    x2723 = x2722 * x49
    x2724 = x1033 * x2096 + x1312 - x60 * (x1033 * x2099 + x1313)
    x2725 = x49 * (x1033 * x2107 + x1311)
    x2726 = x1033 * x2117 + x1324 + 4.0 * x2066
    x2727 = x2726 * x49
    x2728 = x1033 * x2120 + x1337 + x2063 * x98
    x2729 = (
        x1033 * x2113 + x1339 + x2059 * x98 - x60 * (x1033 * x2116 + x1340 + x2061 * x98)
    )
    x2730 = x1033 * x2128 + x1348 + 3.0 * x2078
    x2731 = x2730 * x49
    x2732 = x1033 * x2131 + x116 * x2075 + x1358
    x2733 = (
        x1033 * x2124
        + x116 * x2071
        + x1360
        - x60 * (x1033 * x2127 + x116 * x2073 + x1361)
    )
    x2734 = x1033 * x2139 + x1367 + 2.0 * x2093
    x2735 = x2734 * x49
    x2736 = x1033 * x2142 + x1377 + x2090
    x2737 = x1033 * x2135 + x1379 + x2084 - x60 * (x1033 * x2138 + x1380 + x2087)
    x2738 = x1033 * x2150 + x1385 + x2108
    x2739 = x2738 * x49
    x2740 = x1033 * x2153 + x1395 + x2105
    x2741 = x1033 * x2146 + x1397 + x2098 - x60 * (x1033 * x2149 + x1398 + x2102)
    x2742 = x1033 * x2162 + x1399
    x2743 = x2742 * x49
    x2744 = x1033 * x2165 + x1406
    x2745 = x1033 * x2157 + x1407 - x60 * (x1033 * x2161 + x1408)
    x2746 = x1033 * x2170 + x1419 + 5.0 * x2119
    x2747 = (
        x1033 * x2167 + x1421 + 5.0 * x2112 - x60 * (x1033 * x2168 + x1422 + 5.0 * x2115)
    )
    x2748 = x1033 * x2175 + x1427 + 4.0 * x2130
    x2749 = (
        x1033 * x2172 + x1429 + 4.0 * x2123 - x60 * (x1033 * x2173 + x1430 + 4.0 * x2126)
    )
    x2750 = x1033 * x2180 + x1438 + 3.0 * x2141
    x2751 = (
        x1033 * x2177 + x1442 + 3.0 * x2134 - x60 * (x1033 * x2178 + x1443 + 3.0 * x2137)
    )
    x2752 = x1033 * x2185 + x1448 + 2.0 * x2152
    x2753 = (
        x1033 * x2182 + x1450 + 2.0 * x2145 - x60 * (x1033 * x2183 + x1451 + 2.0 * x2148)
    )
    x2754 = x1033 * x2190 + x1456 + x2164
    x2755 = x1033 * x2187 + x1459 + x2156 - x60 * (x1033 * x2188 + x1460 + x2159)
    x2756 = x1033 * x2195 + x1461
    x2757 = x1033 * x2192 + x1462 - x60 * (x1033 * x2193 + x1463)
    x2758 = x1464 * x1970 + x1468
    x2759 = x1464 * x1972 + x1471
    x2760 = x2759 * x49
    x2761 = x1480 + x1481 * x1962 - x60 * (x1464 * x1966 + x1483)
    x2762 = x49 * (x1464 * x1976 + x1477)
    x2763 = x1033 * x2209 + x1464 * x1988 + x1488
    x2764 = x1033 * x2211 + x1464 * x1990 + x1491
    x2765 = (
        x1033 * x2201
        + x1464 * x1982
        + x1503
        - x60 * (x1033 * x2205 + x1464 * x1985 + x1505)
    )
    x2766 = x1033 * x2214 + x1464 * x1993 + x1500
    x2767 = x1517 * x1962 + x1524
    x2768 = x1519 * x1962 + x1532
    x2769 = x2768 * x49
    x2770 = x1512 * x1962 + x1550 - x60 * (x1506 * x1962 + x1552)
    x2771 = x49 * (x1526 * x1962 + x1546)
    x2772 = x1033 * x2229 + x153 * x2211 + x1560
    x2773 = x1033 * x2230 + x153 * x2214 + x1563
    x2774 = x116 * x2773
    x2775 = (
        x1033 * x2223
        + x153 * x2203
        + x1573
        - x60 * (x1033 * x2226 + x153 * x2207 + x1575)
    )
    x2776 = x1033 * x2233 + x1571 + x2218
    x2777 = x1033 * x1520 + x1033 * x2241 + x1585
    x2778 = x1033 * x1527 + x1033 * x2242 + x1590
    x2779 = (
        x1033 * x2237
        + x1605
        + x2219 * x49
        - x60 * (x1033 * x1509 + x1033 * x2239 + x1607)
    )
    x2780 = x153 * (x1033 * x1535 + x1033 * x2244 + x1602)
    x2781 = x1620 * x1962 + x1627
    x2782 = x1622 * x1962 + x1636
    x2783 = x2782 * x49
    x2784 = x1614 * x1962 + x1652 - x60 * (x1608 * x1962 + x1654)
    x2785 = x49 * (x1630 * x1962 + x1649)
    x2786 = x1033 * x2259 + x1662 + x2231
    x2787 = x1033 * x2260 + x116 * x2233 + x1666
    x2788 = x1033 * x2255 + x1675 + x2225 - x60 * (x1033 * x2257 + x1677 + x2228)
    x2789 = x49 * (x1033 * x2262 + x1673 + 3.0 * x2236)
    x2790 = x1033 * x2274 + x153 * x2242 + x1693
    x2791 = x1033 * x2275 + x1698 + x2245
    x2792 = (
        x1033 * x2268
        + x153 * x2238
        + x1710
        - x60 * (x1033 * x2271 + x153 * x2240 + x1712)
    )
    x2793 = x49 * (x1033 * x2278 + x1707 + 2.0 * x2248)
    x2794 = 321.812152860011 * x154
    x2795 = x1033 * x1623 + x1033 * x2289 + x1725
    x2796 = x1033 * x1631 + x1033 * x2290 + x1730
    x2797 = x153 * x2796
    x2798 = (
        x1033 * x2283
        + x1742
        + x2250 * x49
        - x60 * (x1033 * x1611 + x1033 * x2286 + x1744)
    )
    x2799 = x49 * (x1033 * x2293 + x1739 + x2254)
    x2800 = x1758 * x1962 + x1765
    x2801 = x1760 * x1962 + x1772
    x2802 = x2801 * x49
    x2803 = x1753 * x1962 + x1785 - x60 * (x1745 * x1962 + x1787)
    x2804 = x49 * (x1766 * x1962 + x1782)
    x2805 = x1033 * x2308 + x1789 + 4.0 * x2263
    x2806 = x2805 * x49
    x2807 = x1033 * x2311 + x1797 + x2260 * x98
    x2808 = (
        x1033 * x2304 + x1799 + x2256 * x98 - x60 * (x1033 * x2307 + x1801 + x2258 * x98)
    )
    x2809 = x1033 * x2319 + x1808 + 3.0 * x2279
    x2810 = x2809 * x49
    x2811 = x1033 * x2322 + x116 * x2275 + x1819
    x2812 = (
        x1033 * x2315
        + x116 * x2269
        + x1822
        - x60 * (x1033 * x2318 + x116 * x2272 + x1824)
    )
    x2813 = x1033 * x2330 + x1831 + 2.0 * x2294
    x2814 = x2813 * x49
    x2815 = x1033 * x2333 + x1842 + x2291
    x2816 = x1033 * x2326 + x1845 + x2285 - x60 * (x1033 * x2329 + x1847 + x2288)
    x2817 = x1033 * x2341 + x1853 + x2299
    x2818 = x2817 * x49
    x2819 = x1033 * x1761 + x1033 * x2344 + x1864
    x2820 = (
        x1033 * x2337
        + x1867
        + x2297 * x49
        - x60 * (x1033 * x1750 + x1033 * x2340 + x1869)
    )
    x2821 = x1874 * x1962 + x1877
    x2822 = x2821 * x49
    x2823 = x1888 * x1962 + x1892
    x2824 = x1884 * x1962 + x1896 - x60 * (x1880 * x1962 + x1898)
    x2825 = x1033 * x2354 + x1900 + 5.0 * x2310
    x2826 = (
        x1033 * x2351 + x1902 + 5.0 * x2303 - x60 * (x1033 * x2352 + x1904 + 5.0 * x2306)
    )
    x2827 = x1033 * x2359 + x1910 + 4.0 * x2321
    x2828 = (
        x1033 * x2356 + x1913 + 4.0 * x2314 - x60 * (x1033 * x2357 + x1915 + 4.0 * x2317)
    )
    x2829 = x1033 * x2367 + x1921 + 3.0 * x2332
    x2830 = (
        x1033 * x2362 + x1924 + 3.0 * x2325 - x60 * (x1033 * x2364 + x1926 + 3.0 * x2328)
    )
    x2831 = x1033 * x2372 + x1932 + 2.0 * x2343
    x2832 = (
        x1033 * x2369 + x1935 + 2.0 * x2336 - x60 * (x1033 * x2370 + x1937 + 2.0 * x2339)
    )
    x2833 = x1033 * x2377 + x1943 + x2349
    x2834 = x1033 * x2374 + x1946 + x2346 - x60 * (x1033 * x2375 + x1948 + x2347)
    x2835 = x1953 * x1962 + x1956
    x2836 = x1949 * x1962 + x1959 - x60 * (x1950 * x1962 + x1961)
    x2837 = x1033 * x4 * (x2381 - x2392)
    x2838 = x1033 * x2409 + x2390
    x2839 = x1033 * x2410 + x2397
    x2840 = x4 * (x1033 * x2403 + x2383 - x60 * (x1033 * x2406 + x2386))
    x2841 = x1033 * x2414 + x2402
    x2842 = x1033 * x4 * (x2422 - x2433)
    x2843 = x1033 * x2453 + 2.0 * x2411
    x2844 = x1033 * x2454 + 2.0 * x2415
    x2845 = x116 * x2844
    x2846 = x4 * (x1033 * x2445 + 2.0 * x2405 - x60 * (x1033 * x2449 + 2.0 * x2408))
    x2847 = x1033 * x2459 + x2421
    x2848 = x1033 * x2469 + x2431
    x2849 = x1033 * x2470 + x2438
    x2850 = x4 * (x1033 * x2465 + x2424 - x60 * (x1033 * x2467 + x2427))
    x2851 = x153 * (x1033 * x2473 + x2443)
    x2852 = x1033 * x4 * (x2479 - x2492)
    x2853 = x1033 * x2510 + x2456
    x2854 = x1033 * x2511 + 3.0 * x2460
    x2855 = x4 * (x1033 * x2503 + x2448 - x60 * (x1033 * x2506 + x2452))
    x2856 = x49 * (x1033 * x2515 + 3.0 * x2464)
    x2857 = x153 * x2470
    x2858 = x1033 * x2526 + x2857
    x2859 = x1033 * x2527 + x2474
    x2860 = x153 * x2466
    x2861 = x153 * x2468
    x2862 = x4 * (x1033 * x2520 + x2860 - x60 * (x1033 * x2523 + x2861))
    x2863 = 2.0 * x2478
    x2864 = x49 * (x1033 * x2531 + x2863)
    x2865 = x1033 * x2542 + x2490
    x2866 = x1033 * x2543 + x2497
    x2867 = x153 * x2866
    x2868 = x4 * (x1033 * x2536 + x2481 - x60 * (x1033 * x2539 + x2485))
    x2869 = x49 * (x1033 * x2547 + x2502)
    x2870 = x1033 * x4 * (x2552 - x2563)
    x2871 = x1033 * x2567
    x2872 = x1033 * x2579 + 4.0 * x2516
    x2873 = x2872 * x49
    x2874 = x1033 * x2583 + 4.0 * x2512
    x2875 = x4 * (x1033 * x2574 + 4.0 * x2505 - x60 * (x1033 * x2578 + 4.0 * x2509))
    x2876 = x1033 * x2592 + 3.0 * x2532
    x2877 = x2876 * x49
    x2878 = x1033 * x2596 + x116 * x2527
    x2879 = x4 * (x1033 * x2588 + x116 * x2521 - x60 * (x1033 * x2591 + x116 * x2524))
    x2880 = x1033 * x2605 + 2.0 * x2548
    x2881 = x2880 * x49
    x2882 = x1033 * x2609 + x2544
    x2883 = x4 * (x1033 * x2601 + x2538 - x60 * (x1033 * x2604 + x2541))
    x2884 = x1033 * x2618 + x2567
    x2885 = x2884 * x49
    x2886 = x1033 * x2622 + x2561
    x2887 = x4 * (x1033 * x2614 + x2554 - x60 * (x1033 * x2617 + x2557))
    x2888 = x1033 * x2634
    x2889 = x1033 * x4 * (x2627 - x2637)
    x2890 = x1033 * x2644 + 5.0 * x2582
    x2891 = x4 * (x1033 * x2641 + 5.0 * x2573 - x60 * (x1033 * x2642 + 5.0 * x2576))
    x2892 = x1033 * x2650 + 4.0 * x2595
    x2893 = x4 * (x1033 * x2647 + 4.0 * x2587 - x60 * (x1033 * x2648 + 4.0 * x2590))
    x2894 = 3.0 * x2608
    x2895 = x1033 * x2656 + x2894
    x2896 = 3.0 * x2600
    x2897 = 3.0 * x2603
    x2898 = x4 * (x1033 * x2653 + x2896 - x60 * (x1033 * x2654 + x2897))
    x2899 = x1033 * x2662 + 2.0 * x2621
    x2900 = x4 * (x1033 * x2659 + 2.0 * x2613 - x60 * (x1033 * x2660 + 2.0 * x2616))
    x2901 = x1033 * x2668 + x2634
    x2902 = x4 * (x1033 * x2665 + x2626 - x60 * (x1033 * x2666 + x2629))
    x2903 = x1033 * x4 * (x2671 - x2673)
    x2904 = x1464 * x2387 + x1467
    x2905 = x2904 * x3
    x2906 = x1464 * x2389 + x1470
    x2907 = x2906 * x49
    x2908 = 5.0 * x2907
    x2909 = x1464 * x2381 + x1479 - x60 * (x1464 * x2384 + x1482)
    x2910 = x2909 * x85
    x2911 = x49 * (x1464 * x2396 + x1476)
    x2912 = x1464 * x2409 + x1487
    x2913 = x1464 * x2410 + x1490
    x2914 = x2913 * x49
    x2915 = x1464 * x2403 + x1502 - x60 * (x1464 * x2406 + x1504)
    x2916 = x2915 * x85
    x2917 = x49 * (x1464 * x2414 + x1499)
    x2918 = x1464 * x2428 + x1523 + x2390
    x2919 = x2918 * x3
    x2920 = x1464 * x2430 + x1531 + x2397
    x2921 = x2920 * x49
    x2922 = 4.0 * x2921
    x2923 = x1464 * x2422 + x1549 + x2383 - x60 * (x1464 * x2425 + x1551 + x2386)
    x2924 = x2923 * x85
    x2925 = x49 * (x1464 * x2437 + x1545 + x2402)
    x2926 = x1464 * x2453 + x1559
    x2927 = x1464 * x2454 + x1562
    x2928 = x2927 * x49
    x2929 = 3.0 * x2928
    x2930 = x1464 * x2445 + x1572 - x60 * (x1464 * x2449 + x1574)
    x2931 = x2930 * x85
    x2932 = x49 * (x1464 * x2459 + x1570)
    x2933 = x1464 * x2469 + x1584 + x2411
    x2934 = x1464 * x2470 + x1589 + x2415
    x2935 = x1464 * x2465 + x1604 + x2405 - x60 * (x1464 * x2467 + x1606 + x2408)
    x2936 = x2935 * x85
    x2937 = x153 * (x1464 * x2473 + x1601 + x2420)
    x2938 = x1464 * x2487 + x1626 + 2.0 * x2431
    x2939 = x2938 * x3
    x2940 = x1464 * x2489 + x1635 + 2.0 * x2438
    x2941 = x2940 * x49
    x2942 = 3.0 * x2941
    x2943 = (
        x1464 * x2479 + x1651 + 2.0 * x2424 - x60 * (x1464 * x2483 + x1653 + 2.0 * x2427)
    )
    x2944 = x2943 * x85
    x2945 = x49 * (x1464 * x2496 + x1648 + x2444)
    x2946 = x1464 * x2510 + x1661
    x2947 = x1464 * x2511 + x1665
    x2948 = x2947 * x49
    x2949 = x1464 * x2503 + x1674 - x60 * (x1464 * x2506 + x1676)
    x2950 = x2949 * x85
    x2951 = x49 * (x1464 * x2515 + x1672)
    x2952 = x1464 * x2526 + x1692 + x2455
    x2953 = x1464 * x2527 + x1697 + x2460
    x2954 = x153 * x2953
    x2955 = x1464 * x2520 + x1709 + x2447 - x60 * (x1464 * x2523 + x1711 + x2451)
    x2956 = x2955 * x85
    x2957 = x49 * (x1464 * x2531 + x1706 + x2464)
    x2958 = x1464 * x2542 + x1724 + x2857
    x2959 = x1464 * x2543 + x1729 + x2474
    x2960 = x153 * x2959
    x2961 = x1464 * x2536 + x1741 + x2860 - x60 * (x1464 * x2539 + x1743 + x2861)
    x2962 = x2961 * x85
    x2963 = x49 * (x1464 * x2547 + x1738 + x2863)
    x2964 = x1464 * x2558 + x1764 + x2491
    x2965 = x2964 * x3
    x2966 = x1464 * x2560 + x1771 + 3.0 * x2497
    x2967 = x2966 * x49
    x2968 = 2.0 * x2967
    x2969 = x1464 * x2552 + x1784 + x2482 - x60 * (x1464 * x2555 + x1786 + x2486)
    x2970 = x2969 * x85
    x2971 = x49 * (x1464 * x2566 + x1781 + 3.0 * x2502)
    x2972 = x1464 * x2579 + x1788
    x2973 = x2972 * x49
    x2974 = x1464 * x2583 + x1796
    x2975 = x1464 * x2574 + x1798 - x60 * (x1464 * x2578 + x1800)
    x2976 = x2975 * x85
    x2977 = x1464 * x2592 + x1807 + x2516
    x2978 = x2977 * x49
    x2979 = x1464 * x2596 + x1818 + x2512
    x2980 = x1464 * x2588 + x1821 + x2505 - x60 * (x1464 * x2591 + x1823 + x2509)
    x2981 = x2980 * x85
    x2982 = x1464 * x2605 + x1830 + 2.0 * x2532
    x2983 = x2982 * x49
    x2984 = x1464 * x2609 + x1841 + x2528
    x2985 = x1464 * x2601 + x1844 + x2522 - x60 * (x1464 * x2604 + x1846 + x2525)
    x2986 = x2985 * x85
    x2987 = x1464 * x2618 + x1852 + 3.0 * x2548
    x2988 = x2987 * x49
    x2989 = x116 * x2543 + x1464 * x2622 + x1863
    x2990 = (
        x116 * x2537
        + x1464 * x2614
        + x1866
        - x60 * (x116 * x2540 + x1464 * x2617 + x1868)
    )
    x2991 = x2990 * x85
    x2992 = x1464 * x2631 + x1876 + 4.0 * x2567
    x2993 = x2992 * x49
    x2994 = x1464 * x2635 + x1891 + 4.0 * x2561
    x2995 = x2994 * x3
    x2996 = (
        x1464 * x2627 + x1895 + 4.0 * x2554 - x60 * (x1464 * x2630 + x1897 + 4.0 * x2557)
    )
    x2997 = x2996 * x85
    x2998 = x1464 * x2644 + x1899
    x2999 = x1464 * x2641 + x1901 - x60 * (x1464 * x2642 + x1903)
    x3000 = x2999 * x85
    x3001 = x1464 * x2650 + x1909 + x2582
    x3002 = x1464 * x2647 + x1912 + x2573 - x60 * (x1464 * x2648 + x1914 + x2576)
    x3003 = x3002 * x85
    x3004 = x1464 * x2656 + x1920 + 2.0 * x2595
    x3005 = (
        x1464 * x2653 + x1923 + 2.0 * x2587 - x60 * (x1464 * x2654 + x1925 + 2.0 * x2590)
    )
    x3006 = x3005 * x85
    x3007 = x1464 * x2662 + x1931 + x2894
    x3008 = x1464 * x2659 + x1934 + x2896 - x60 * (x1464 * x2660 + x1936 + x2897)
    x3009 = x3008 * x85
    x3010 = x1464 * x2668 + x1942 + 4.0 * x2621
    x3011 = (
        x1464 * x2665 + x1945 + 4.0 * x2613 - x60 * (x1464 * x2666 + x1947 + 4.0 * x2616)
    )
    x3012 = x3011 * x85
    x3013 = x1464 * x2674 + x1955 + 5.0 * x2634
    x3014 = (
        x1464 * x2671 + x1958 + 5.0 * x2626 - x60 * (x1464 * x2672 + x1960 + 5.0 * x2629)
    )
    x3015 = x3014 * x85
    x3016 = x1033 * x2678 + x134 * x1975
    x3017 = x49 * (x1033 * x2679 + x134 * x1980)
    x3018 = x1033 * x2683 + x134 * x1992 + x2680
    x3019 = x1033 * x2684 + x134 * x1996 + x2682
    x3020 = x1033 * x2687 + x134 * x2008
    x3021 = x49 * (x1033 * x2688 + x134 * x2012)
    x3022 = x1033 * x2692 + x134 * x2024 + x153 * x2684
    x3023 = x116 * (x1033 * x2693 + x134 * x2027 + x153 * x2686)
    x3024 = x1033 * x2697 + x134 * x2036 + x2689
    x3025 = x1033 * x2698 + x134 * x2040 + x2691
    x3026 = x1033 * x2701 + x134 * x2052
    x3027 = x49 * (x1033 * x2702 + x134 * x2056)
    x3028 = x1033 * x2706 + x134 * x2064 + x2694
    x3029 = x1033 * x2707 + x116 * x2696 + x134 * x2068
    x3030 = x1033 * x2711 + x134 * x2076 + x153 * x2698
    x3031 = x1033 * x2712 + x134 * x2080 + x2700
    x3032 = x1033 * x2716 + x134 * x2091 + x2703
    x3033 = x153 * (x1033 * x2717 + x134 * x2095 + x2705)
    x3034 = x1033 * x2721 + x134 * x2106
    x3035 = x49 * (x1033 * x2722 + x134 * x2110)
    x3036 = x49 * (x1033 * x2726 + x134 * x2118 + 4.0 * x2709)
    x3037 = x1033 * x2728 + x134 * x2121 + x2707 * x98
    x3038 = x49 * (x1033 * x2730 + x134 * x2129 + 3.0 * x2714)
    x3039 = x1033 * x2732 + x116 * x2712 + x134 * x2132
    x3040 = x49 * (x1033 * x2734 + x134 * x2140 + 2.0 * x2720)
    x3041 = x1033 * x2736 + x134 * x2143 + x2718
    x3042 = x49 * (x1033 * x2738 + x134 * x2151 + x2725)
    x3043 = x1033 * x2740 + x134 * x2154 + x2723
    x3044 = x49 * (x1033 * x2742 + x134 * x2163)
    x3045 = x1033 * x2744 + x134 * x2166
    x3046 = x1033 * x2746 + x134 * x2171 + 5.0 * x2727
    x3047 = x252 * x3
    x3048 = x1033 * x2748 + x134 * x2176 + 4.0 * x2731
    x3049 = x1097 * x3
    x3050 = x1033 * x2750 + x134 * x2181 + 3.0 * x2735
    x3051 = x1153 * x3
    x3052 = x1033 * x2752 + x134 * x2186 + 2.0 * x2739
    x3053 = x1231 * x3
    x3054 = x1033 * x2754 + x134 * x2191 + x2743
    x3055 = x1033 * x2756 + x134 * x2196
    x3056 = x1033 * x2758 + x2197
    x3057 = x49 * (x1033 * x2759 + x2199)
    x3058 = x1033 * x2763 + x2212 + x2760
    x3059 = x1033 * x2764 + x2216 + x2762
    x3060 = x1033 * x2767 + x2220
    x3061 = x49 * (x1033 * x2768 + x2222)
    x3062 = x1033 * x2772 + x153 * x2764 + x2232
    x3063 = x116 * (x1033 * x2773 + x153 * x2766 + x2235)
    x3064 = x1033 * x2777 + x2243 + x2769
    x3065 = x1033 * x2778 + x2247 + x2771
    x3066 = x1033 * x2781 + x2251
    x3067 = x49 * (x1033 * x2782 + x2253)
    x3068 = x1033 * x2786 + x2261 + x2774
    x3069 = x1033 * x2787 + x116 * x2776 + x2265
    x3070 = x1033 * x2790 + x153 * x2778 + x2276
    x3071 = x1033 * x2791 + x2281 + x2780
    x3072 = x1033 * x2795 + x2292 + x2783
    x3073 = x153 * (x1033 * x2796 + x2296 + x2785)
    x3074 = x1033 * x2800 + x2298
    x3075 = x49 * (x1033 * x2801 + x2301)
    x3076 = x49 * (x1033 * x2805 + x2309 + 4.0 * x2789)
    x3077 = x1033 * x2807 + x2312 + x2787 * x98
    x3078 = x49 * (x1033 * x2809 + x2320 + 3.0 * x2793)
    x3079 = x1033 * x2811 + x116 * x2791 + x2323
    x3080 = x49 * (x1033 * x2813 + x2331 + 2.0 * x2799)
    x3081 = x1033 * x2815 + x2334 + x2797
    x3082 = x49 * (x1033 * x2817 + x2342 + x2804)
    x3083 = x1033 * x2819 + x2345 + x2802
    x3084 = x49 * (x1033 * x2821 + x2348)
    x3085 = x1033 * x2823 + x2350
    x3086 = x1033 * x2825 + x2355 + 5.0 * x2806
    x3087 = x3 * x482
    x3088 = x1033 * x2827 + x2360 + 4.0 * x2810
    x3089 = x1181 * x3
    x3090 = x1033 * x2829 + x2368 + 3.0 * x2814
    x3091 = x2042 * x3
    x3092 = x1033 * x2831 + x2373 + 2.0 * x2818
    x3093 = x2266 * x3
    x3094 = x1033 * x2833 + x2378 + x2822
    x3095 = x1033 * x2835 + x2379
    x3096 = x1962 * x2387 + x2394
    x3097 = x49 * (x1962 * x2389 + x2401)
    x3098 = x1033 * x2390 + x1033 * x2838 + x2413
    x3099 = x1033 * x2397 + x1033 * x2839 + x2418
    x3100 = x1962 * x2428 + x2435
    x3101 = x49 * (x1962 * x2430 + x2442)
    x3102 = x1033 * x2843 + x153 * x2839 + x2458
    x3103 = x116 * (x1033 * x2844 + x153 * x2841 + x2463)
    x3104 = x1033 * x2431 + x1033 * x2848 + x2472
    x3105 = x1033 * x2438 + x1033 * x2849 + x2477
    x3106 = x1962 * x2487 + x2494
    x3107 = x49 * (x1962 * x2489 + x2501)
    x3108 = x1033 * x2853 + x2514 + x2845
    x3109 = x1033 * x2854 + x116 * x2847 + x2519
    x3110 = x1033 * x2858 + x153 * x2849 + x2530
    x3111 = x1033 * x2859 + x2535 + x2851
    x3112 = x1033 * x2490 + x1033 * x2865 + x2546
    x3113 = x153 * (x1033 * x2497 + x1033 * x2866 + x2551)
    x3114 = x1962 * x2558 + x2565
    x3115 = x49 * (x1962 * x2560 + x2571)
    x3116 = x49 * (x1033 * x2872 + x2581 + 4.0 * x2856)
    x3117 = x1033 * x2874 + x2585 + x2854 * x98
    x3118 = x49 * (x1033 * x2876 + x2594 + 3.0 * x2864)
    x3119 = x1033 * x2878 + x116 * x2859 + x2598
    x3120 = x49 * (x1033 * x2880 + x2607 + 2.0 * x2869)
    x3121 = x1033 * x2882 + x2611 + x2867
    x3122 = x49 * (x1033 * x2884 + x2620 + x2871)
    x3123 = x1033 * x2561 + x1033 * x2886 + x2624
    x3124 = x49 * (x1962 * x2631 + x2633)
    x3125 = x1962 * x2635 + x2639
    x3126 = x1033 * x2890 + x2646 + 5.0 * x2873
    x3127 = x3 * x682
    x3128 = x1033 * x2892 + x2652 + 4.0 * x2877
    x3129 = x1266 * x3
    x3130 = x1033 * x2895 + x2658 + 3.0 * x2881
    x3131 = x1033 * x2899 + x2664 + 2.0 * x2885
    x3132 = x1033 * x2901 + x2670 + x2888
    x3133 = x1962 * x2674 + x2677
    x3134 = x1033 * x2912 + x2907
    x3135 = x1033 * x2913 + x2911
    x3136 = x1033 * x2926 + 2.0 * x2914
    x3137 = x116 * (x1033 * x2927 + 2.0 * x2917)
    x3138 = x1033 * x2933 + x2921
    x3139 = x1033 * x2934 + x2925
    x3140 = x1033 * x2946 + x2929
    x3141 = x1033 * x2947 + 3.0 * x2932
    x3142 = x153 * x2934
    x3143 = x1033 * x2952 + x3142
    x3144 = x1033 * x2953 + x2937
    x3145 = x1033 * x2958 + x2941
    x3146 = x153 * (x1033 * x2959 + x2945)
    x3147 = x49 * (x1033 * x2972 + 4.0 * x2951)
    x3148 = x1033 * x2974 + 4.0 * x2948
    x3149 = x49 * (x1033 * x2977 + 3.0 * x2957)
    x3150 = x1033 * x2979 + x116 * x2953
    x3151 = x49 * (x1033 * x2982 + 2.0 * x2963)
    x3152 = x1033 * x2984 + x2960
    x3153 = x49 * (x1033 * x2987 + x2971)
    x3154 = x1033 * x2989 + x2967
    x3155 = x1033 * x2993
    x3156 = x1033 * x2998 + 5.0 * x2973
    x3157 = x1033 * x3001 + 4.0 * x2978
    x3158 = 3.0 * x2983
    x3159 = x1033 * x3004 + x3158
    x3160 = x1033 * x3007 + 2.0 * x2988
    x3161 = x1033 * x3010 + x2993
    x3162 = x134 * x2393 + x1464 * x2904
    x3163 = x49 * (x134 * x2400 + x1464 * x2906)
    x3164 = x134 * x2412 + x1464 * x2912
    x3165 = x49 * (x134 * x2417 + x1464 * x2913)
    x3166 = x134 * x2434 + x1464 * x2918 + x2907
    x3167 = x49 * (x134 * x2441 + x1464 * x2920 + x2911)
    x3168 = x134 * x2457 + x1464 * x2926
    x3169 = x49 * (x134 * x2462 + x1464 * x2927)
    x3170 = 3.0 * x3169
    x3171 = x134 * x2471 + x1464 * x2933 + x2914
    x3172 = x134 * x2476 + x1464 * x2934 + x2917
    x3173 = x134 * x2493 + x1464 * x2938 + 2.0 * x2921
    x3174 = x49 * (x134 * x2500 + x1464 * x2940 + 2.0 * x2925)
    x3175 = 3.0 * x3174
    x3176 = x134 * x2513 + x1464 * x2946
    x3177 = x49 * (x134 * x2518 + x1464 * x2947)
    x3178 = x134 * x2529 + x1464 * x2952 + x2928
    x3179 = x134 * x2534 + x1464 * x2953 + x2932
    x3180 = x153 * x3179
    x3181 = x134 * x2545 + x1464 * x2958 + x3142
    x3182 = x134 * x2550 + x1464 * x2959 + x2937
    x3183 = x153 * x3182
    x3184 = x134 * x2564 + x1464 * x2964 + x2942
    x3185 = x49 * (x134 * x2570 + x1464 * x2966 + 3.0 * x2945)
    x3186 = x49 * (x134 * x2580 + x1464 * x2972)
    x3187 = x134 * x2584 + x1464 * x2974
    x3188 = x49 * (x134 * x2593 + x1464 * x2977 + x2951)
    x3189 = x134 * x2597 + x1464 * x2979 + x2948
    x3190 = x49 * (x134 * x2606 + x1464 * x2982 + 2.0 * x2957)
    x3191 = x134 * x2610 + x1464 * x2984 + x2954
    x3192 = x49 * (x134 * x2619 + x1464 * x2987 + 3.0 * x2963)
    x3193 = x116 * x2959 + x134 * x2623 + x1464 * x2989
    x3194 = x49 * (x134 * x2632 + x1464 * x2992 + 4.0 * x2971)
    x3195 = x134 * x2638 + x1464 * x2994 + 4.0 * x2967
    x3196 = x134 * x2645 + x1464 * x2998
    x3197 = x134 * x2651 + x1464 * x3001 + x2973
    x3198 = x134 * x2657 + x1464 * x3004 + 2.0 * x2978
    x3199 = x134 * x2663 + x1464 * x3007 + x3158
    x3200 = x134 * x2669 + x1464 * x3010 + 4.0 * x2988
    x3201 = x134 * x2676 + x1464 * x3013 + 5.0 * x2993
    x3202 = x1033 * x252
    x3203 = x1033 * x1097
    x3204 = x153 * x3172
    x3205 = 3.0 * x3190

    # 441 item(s)
    result[0, 0] = numpy.sum(
        x155
        * (
            x149
            * (
                x133 * x3
                + x144 * x91
                + x4 * (x131 * x3 + 5.0 * x138 - x140 * x60)
                - x60 * (x123 * x3 + x148 * x91 + x4 * (x140 - x60 * x74))
            )
            + x3
            * (
                x118 * x91
                - x134 * (x123 * x60 - x133)
                + x3
                * (
                    x3 * (x3 * (x81 + x83) + x86 + x90 * x91)
                    - x4 * (x59 * x60 - x74)
                    + x91 * x99
                )
            )
            + x91
            * (
                x118 * x3
                + x134 * (x144 - x148 * x60)
                + x98
                * (
                    x116 * (x115 * x3 + x150 + x153 * (x151 + x152 * x53))
                    + x117 * x3
                    + x4 * (x143 - x147 * x60)
                )
            )
        )
    )
    result[0, 1] = numpy.sum(
        x252
        * (
            x149
            * (
                x236 * x3
                + x244 * x98
                + x4 * (x234 * x3 + 4.0 * x238 - x240 * x60)
                - x60 * (x231 * x3 + x249 * x98 - x4 * (x187 * x60 - x240))
            )
            + x3
            * (
                -x134 * (x231 * x60 - x236)
                + x224 * x98
                + x3
                * (
                    x208 * x98
                    + x3 * (x196 + x200 * x98 + x3 * (x192 + x194))
                    - x4 * (x177 * x60 - x187)
                )
            )
            + x98
            * (
                x116
                * (
                    x153 * (x219 * x3 + x222 * x3 + x251)
                    + x223 * x3
                    + x4 * (x243 - x248 * x60)
                )
                + x134 * (x244 - x249 * x60)
                + x224 * x3
            )
        )
    )
    result[0, 2] = numpy.sum(
        x252
        * (
            x149
            * (
                x3 * x333
                + x341 * x98
                + x4 * (x3 * x331 + 4.0 * x335 - x337 * x60)
                - x60 * (x3 * x327 + x345 * x98 - x4 * (x283 * x60 - x337))
            )
            + x3
            * (
                -x134 * (x327 * x60 - x333)
                + x3
                * (
                    x3 * (x292 + x296 * x98 + x3 * (x288 + x290))
                    + x305 * x98
                    - x4 * (x273 * x60 - x283)
                )
                + x321 * x98
            )
            + x98
            * (
                x116
                * (
                    x153 * (x3 * x316 + x3 * x319 + x347)
                    + x3 * x320
                    + x4 * (x340 - x344 * x60)
                )
                + x134 * (x341 - x345 * x60)
                + x3 * x321
            )
        )
    )
    result[0, 3] = numpy.sum(
        x431
        * (
            x116
            * (
                x134 * (x424 - x428 * x60)
                + x153
                * (x3 * x403 + x4 * (x423 - x427 * x60) + x49 * (x349 * x429 + x430))
                + x3 * x404
            )
            + x149
            * (
                x116 * x424
                + x3 * x413
                + x4 * (x3 * x411 + x417 - x420 * x60)
                - x60 * (x116 * x428 + x3 * x407 - x4 * (x374 * x60 - x420))
            )
            + x3
            * (
                x116 * x404
                - x134 * (x407 * x60 - x413)
                + x3
                * (
                    x116 * x393
                    + x3 * (x116 * x387 + x3 * (x380 + x382) + x384)
                    - x4 * (x364 * x60 - x374)
                )
            )
        )
    )
    result[0, 4] = numpy.sum(
        x482
        * (
            x116
            * (
                x134 * (x476 - x479 * x60)
                + x153
                * (x3 * x461 + x4 * (x475 - x478 * x60) + x49 * (x429 * x451 + x481))
                + x3 * x462
            )
            + x149
            * (
                x116 * x476
                + x3 * x470
                + x4 * (x116 * x472 + x3 * x468 - x473 * x60)
                - x60 * (x116 * x479 + x3 * x465 - x4 * (x443 * x60 - x473))
            )
            + x3
            * (
                x116 * x462
                - x134 * (x465 * x60 - x470)
                + x3
                * (
                    x116 * x455
                    + x3 * (x116 * x448 + x3 * (x116 * x433 + x3 * x445) + x446)
                    - x4 * (x437 * x60 - x443)
                )
            )
        )
    )
    result[0, 5] = numpy.sum(
        x431
        * (
            x116
            * (
                x134 * (x569 - x574 * x60)
                + x153
                * (x3 * x544 + x4 * (x568 - x573 * x60) + x49 * (x429 * x484 + x575))
                + x3 * x545
            )
            + x149
            * (
                x116 * x569
                + x3 * x557
                + x4 * (x3 * x555 + x562 - x565 * x60)
                - x60 * (x116 * x574 + x3 * x548 - x4 * (x514 * x60 - x565))
            )
            + x3
            * (
                x116 * x545
                - x134 * (x548 * x60 - x557)
                + x3
                * (
                    x116 * x534
                    + x3 * (x116 * x528 + x3 * (x521 + x523) + x525)
                    - x4 * (x502 * x60 - x514)
                )
            )
        )
    )
    result[0, 6] = numpy.sum(
        x631
        * (
            x149
            * (
                x153 * x628
                + x3 * x623
                + x4 * (x153 * x626 + x3 * x621 - x60 * x625)
                - x60 * (x153 * x630 + x3 * x619 - x4 * (x590 * x60 - x625))
            )
            + x153
            * (
                -x134 * (x60 * x630 - x628)
                + x3 * x49 * (x4 * (x585 - x611) + x613)
                + x3 * x614
            )
            + x3
            * (
                -x134 * (x60 * x619 - x623)
                + x153 * x614
                + x3
                * (
                    x153 * x603
                    + x3 * (x153 * x599 + x3 * (x153 * x595 + x594) + x597)
                    - x4 * (x584 * x60 - x590)
                )
            )
        )
    )
    result[0, 7] = numpy.sum(
        x682
        * (
            x149
            * (
                x153 * x679
                + x3 * x672
                + x4 * (x3 * x670 - x60 * x675 + x677)
                - x60 * (x153 * x681 + x3 * x668 - x4 * (x60 * x646 - x675))
            )
            + x153
            * (
                -x134 * (x60 * x681 - x679)
                + x3 * x49 * (x4 * (x641 - x660) + x662)
                + x3 * x663
            )
            + x3
            * (
                -x134 * (x60 * x668 - x672)
                + x153 * x663
                + x3
                * (
                    x153 * x655
                    + x3 * (x153 * x653 + x3 * (x3 * x649 + x651) + x652)
                    - x4 * (x60 * x640 - x646)
                )
            )
        )
    )
    result[0, 8] = numpy.sum(
        x682
        * (
            x149
            * (
                x153 * x728
                + x3 * x721
                + x4 * (x3 * x719 - x60 * x724 + x726)
                - x60 * (x153 * x730 + x3 * x717 - x4 * (x60 * x694 - x724))
            )
            + x153
            * (
                -x134 * (x60 * x730 - x728)
                + x3 * x49 * (x4 * (x157 * x496 - x709) + x711)
                + x3 * x712
            )
            + x3
            * (
                -x134 * (x60 * x717 - x721)
                + x153 * x712
                + x3
                * (
                    x153 * x704
                    + x3 * (x153 * x701 + x3 * (x3 * x696 + x698) + x699)
                    - x4 * (x60 * x687 - x694)
                )
            )
        )
    )
    result[0, 9] = numpy.sum(
        x631
        * (
            x149
            * (
                x153 * x788
                + x3 * x782
                + x4 * (x153 * x785 + x3 * x780 - x60 * x784)
                - x60 * (x153 * x791 + x3 * x777 - x4 * (x60 * x747 - x784))
            )
            + x153
            * (
                -x134 * (x60 * x791 - x788)
                + x3 * x49 * (x4 * (x741 - x768) + x770)
                + x3 * x771
            )
            + x3
            * (
                -x134 * (x60 * x777 - x782)
                + x153 * x771
                + x3
                * (
                    x153 * x760
                    + x3 * (x153 * x757 + x3 * (x153 * x753 + x752) + x755)
                    - x4 * (x60 * x740 - x747)
                )
            )
        )
    )
    result[0, 10] = numpy.sum(
        x252
        * (
            x149
            * (
                x3 * x826
                + x4 * (x3 * x828 - x60 * x830 + x827)
                + x49 * x806
                - x60 * (x3 * x821 - x4 * (x60 * x813 - x830) + x49 * x803)
            )
            + x3
            * (
                -x134 * (x60 * x821 - x826)
                + x3
                * (
                    x3 * (x3 * (x814 + x816) + x817 + x818)
                    - x4 * (x60 * x809 - x813)
                    + x49 * x799
                )
                + x49 * x800
            )
            - x49 * (x134 * (x60 * x803 - x806) - x3 * x800)
        )
    )
    result[0, 11] = numpy.sum(
        x482
        * (
            x149
            * (
                x3 * x857
                + x4 * (x3 * x859 - x60 * x861 + x858)
                + x49 * x843
                - x60 * (x3 * x853 - x4 * (x60 * x848 - x861) + x49 * x840)
            )
            + x3
            * (
                -x134 * (x60 * x853 - x857)
                + x3
                * (
                    x3 * (x3 * x849 + x3 * (x152 * x850 + x849) + x851)
                    - x4 * (x60 * x845 - x848)
                    + x49 * x836
                )
                + x49 * x837
            )
            - x49 * (x134 * (x60 * x840 - x843) - x3 * x837)
        )
    )
    result[0, 12] = numpy.sum(
        x682
        * (
            x149
            * (
                x3 * x888
                + x4 * (x3 * x890 - x60 * x892 + x889)
                + x49 * x874
                - x60 * (x3 * x884 - x4 * (x60 * x879 - x892) + x49 * x871)
            )
            + x3
            * (
                -x134 * (x60 * x884 - x888)
                + x3
                * (
                    x3 * (x3 * x880 + x3 * (x152 * x881 + x880) + x882)
                    - x4 * (x60 * x876 - x879)
                    + x49 * x867
                )
                + x49 * x868
            )
            - x49 * (x134 * (x60 * x871 - x874) - x3 * x868)
        )
    )
    result[0, 13] = numpy.sum(
        x482
        * (
            x149
            * (
                x3 * x915
                + x4 * (x3 * x919 + x49 * x916 - x60 * x918)
                + x49 * x903
                - x60 * (x3 * x912 - x4 * (x60 * x907 - x918) + x49 * x901)
            )
            + x3
            * (
                -x134 * (x60 * x912 - x915)
                + x3
                * (
                    x3 * (x3 * x908 + x3 * (x220 * x748 + x908) + x909)
                    - x4 * (x60 * x905 - x907)
                    + x49 * x897
                )
                + x49 * x898
            )
            - x49 * (x134 * (x60 * x901 - x903) - x3 * x898)
        )
    )
    result[0, 14] = numpy.sum(
        x252
        * (
            x149
            * (
                x3 * x955
                + x4 * (x3 * x957 - x60 * x959 + x956)
                + x49 * x935
                - x60 * (x3 * x950 - x4 * (x60 * x942 - x959) + x49 * x932)
            )
            + x3
            * (
                -x134 * (x60 * x950 - x955)
                + x3
                * (
                    x3 * (x3 * (x943 + x945) + x946 + x947)
                    - x4 * (x60 * x938 - x942)
                    + x49 * x927
                )
                + x49 * x928
            )
            - x49 * (x134 * (x60 * x932 - x935) - x3 * x928)
        )
    )
    result[0, 15] = numpy.sum(
        x155
        * x3
        * (
            -x134 * (x60 * x970 - x974)
            - x149
            * (x4 * (x60 * x968 - x972) + x60 * (-x4 * (x60 * x960 - x968) + x970) - x974)
            + x3 * (x3 * (x965 + x966) + x4 * (x3 * x960 - x963))
        )
    )
    result[0, 16] = numpy.sum(
        x252
        * x3
        * (
            -x134 * (x60 * x982 - x985)
            - x149
            * (x4 * (x60 * x980 - x983) + x60 * (-x4 * (x60 * x975 - x980) + x982) - x985)
            + x3**2 * (x4 * (x975 - x977) + x429 * x978 + x979)
        )
    )
    result[0, 17] = numpy.sum(
        x3
        * x431
        * (
            -x134 * (x60 * x993 - x996)
            - x149
            * (x4 * (x60 * x991 - x994) + x60 * (-x4 * (x60 * x986 - x991) + x993) - x996)
            + x3**2 * (x4 * (x986 - x988) + x429 * x989 + x990)
        )
    )
    result[0, 18] = numpy.sum(
        x3
        * x631
        * (
            -x134 * (x1004 * x60 - x1007)
            - x149
            * (
                -x1007
                + x4 * (x1002 * x60 - x1005)
                + x60 * (x1004 + x4 * (x1002 - x60 * x997))
            )
            + x3**2 * (x1000 * x429 + x1001 + x4 * (x997 - x999))
        )
    )
    result[0, 19] = numpy.sum(
        x252
        * (
            x149
            * (
                x1017 * x3
                + x157 * x4 * (x3 * x930 - x60 * x953)
                - x3 * x60 * (x1014 - x4 * (x1012 * x60 - x157 * x920))
            )
            - x3
            * (
                x134 * (x1014 * x60 - x1017)
                - x3**2 * (x1010 + x1011 * x429 - x4 * (x1009 - x157 * x921))
            )
        )
    )
    result[0, 20] = numpy.sum(
        x155
        * x3
        * (
            -x134 * (x1028 * x60 - x1032)
            - x149
            * (
                -x1032
                + x4 * (x1026 * x60 - x1030)
                + x60 * (x1028 - x4 * (x1018 * x60 - x1026))
            )
            + x3 * (x3 * (x1023 + x1024) + x4 * (x1018 * x3 - x1021))
        )
    )
    result[1, 0] = numpy.sum(
        0.5
        * x252
        * (
            x134
            * (
                2.0 * x1039 * x3
                + 2.0 * x1055 * x91
                + x1059
                - x60 * (2.0 * x1036 * x3 + 2.0 * x1053 * x91 + x1060)
            )
            + x3
            * (
                2.0 * x1052 * x91
                + x3 * (x1043 + 2.0 * x1046 * x91 + 2.0 * x3 * (x1040 + x1041))
                - 2.0 * x4 * (x1036 * x60 - x1039)
            )
            + x91
            * (
                2.0 * x1052 * x3
                - 2.0 * x4 * (x1053 * x60 - x1055)
                + x98
                * (2.0 * x1051 * x3 + x1056 + 2.0 * x116 * (x1033 * x110 + x1057 * x1058))
            )
        )
    )
    result[1, 1] = numpy.sum(
        0.5
        * x1097
        * (
            x134
            * (
                2.0 * x1070 * x3
                + 2.0 * x1091 * x98
                + x1095
                - x60 * (2.0 * x1065 * x3 + 2.0 * x1089 * x98 + x1096)
            )
            + x3
            * (
                2.0 * x1087 * x98
                + x3 * (x1075 + 2.0 * x1078 * x98 + 2.0 * x3 * (x1072 * x3 + x1074 * x98))
                - 2.0 * x4 * (x1065 * x60 - x1070)
            )
            + x98
            * (
                2.0 * x1087 * x3
                + x116 * (2.0 * x1086 * x3 + x1092 + 2.0 * x153 * (x1084 * x3 + x1093))
                - 2.0 * x4 * (x1089 * x60 - x1091)
            )
        )
    )
    result[1, 2] = numpy.sum(
        0.5
        * x1097
        * (
            x134
            * (
                2.0 * x1100 * x3
                + 2.0 * x1108 * x98
                + x1111
                - x60 * (2.0 * x1098 * x3 + 2.0 * x1106 * x98 + x1112)
            )
            + x3
            * (
                2.0 * x1105 * x98
                + x3 * (2.0 * x1033 * x3 * (x288 + x290) + x1101 + 2.0 * x1102 * x98)
                - 2.0 * x4 * (x1098 * x60 - x1100)
            )
            + x98
            * (
                2.0 * x1105 * x3
                + x116 * (2.0 * x1104 * x3 + x1109 + 2.0 * x153 * (x1033 * x318 + x1110))
                - 2.0 * x4 * (x1106 * x60 - x1108)
            )
        )
    )
    result[1, 3] = numpy.sum(
        0.5
        * x1153
        * (
            x116
            * (
                2.0 * x1144 * x3
                + x153 * (2.0 * x1142 * x3 + 2.0 * x1143 * x3 + x1149)
                - 2.0 * x4 * (x1146 * x60 - x1148)
            )
            + x134
            * (
                2.0 * x1126 * x3
                + 2.0 * x1148 * x116
                + x1151
                - x60 * (2.0 * x1120 * x3 + 2.0 * x1146 * x116 + x1152)
            )
            + x3
            * (
                2.0 * x1144 * x116
                + x3 * (x1134 + 2.0 * x1137 * x116 + 2.0 * x3 * (x1129 * x3 + x1133))
                - 2.0 * x4 * (x1120 * x60 - x1126)
            )
        )
    )
    result[1, 4] = numpy.sum(
        0.5
        * x1181
        * (
            x116
            * (
                2.0 * x1170 * x3
                + x153 * (2.0 * x1168 * x3 + 2.0 * x1169 * x3 + x1177)
                - 2.0 * x4 * (x1173 * x60 - x1176)
            )
            + x134
            * (
                2.0 * x1159 * x3
                + 2.0 * x116 * x1176
                + x1179
                - x60 * (2.0 * x1156 * x3 + 2.0 * x116 * x1173 + x1180)
            )
            + x3
            * (
                2.0 * x116 * x1170
                + x3
                * (2.0 * x116 * x1165 + x1162 + 2.0 * x3 * (x116 * x1161 + x1160 * x3))
                - 2.0 * x4 * (x1156 * x60 - x1159)
            )
        )
    )
    result[1, 5] = numpy.sum(
        0.5
        * x1153
        * (
            x116
            * (
                2.0 * x1190 * x3
                + x153 * (2.0 * x1033 * x543 + 2.0 * x1189 * x3 + x1194)
                - 2.0 * x4 * (x1191 * x60 - x1193)
            )
            + x134
            * (
                2.0 * x116 * x1193
                + 2.0 * x1184 * x3
                + x1195
                - x60 * (2.0 * x116 * x1191 + 2.0 * x1182 * x3 + x1196)
            )
            + x3
            * (
                2.0 * x116 * x1190
                + x3 * (2.0 * x1033 * x3 * (x521 + x523) + 2.0 * x116 * x1186 + x1185)
                - 2.0 * x4 * (x1182 * x60 - x1184)
            )
        )
    )
    result[1, 6] = numpy.sum(
        0.5
        * x1231
        * (
            x134
            * (
                2.0 * x1208 * x3
                + 2.0 * x1226 * x153
                + x1229
                - x60 * (2.0 * x1203 * x3 + 2.0 * x1223 * x153 + x1230)
            )
            + x153
            * (
                2.0 * x1220 * x3
                - 2.0 * x4 * (x1223 * x60 - x1226)
                + x49 * (2.0 * x1216 * x429 + x1227)
            )
            + x3
            * (
                2.0 * x1220 * x153
                + x3
                * (x1214 + 2.0 * x1218 * x153 + 2.0 * x3 * (x1210 * x3 + x1213 * x153))
                - 2.0 * x4 * (x1203 * x60 - x1208)
            )
        )
    )
    result[1, 7] = numpy.sum(
        0.5
        * x1266
        * (
            x134
            * (
                2.0 * x1240 * x3
                + 2.0 * x1259 * x153
                + x1264
                - x60 * (2.0 * x1236 * x3 + 2.0 * x1255 * x153 + x1265)
            )
            + x153
            * (
                2.0 * x1251 * x3
                - 2.0 * x4 * (x1255 * x60 - x1259)
                + x49 * (2.0 * x1247 * x429 + x1260)
            )
            + x3
            * (
                2.0 * x1251 * x153
                + x3
                * (x1244 + 2.0 * x1249 * x153 + 2.0 * x3 * (x1242 * x3 + x1243 * x153))
                - 2.0 * x4 * (x1236 * x60 - x1240)
            )
        )
    )
    result[1, 8] = numpy.sum(
        0.5
        * x1266
        * (
            x134
            * (
                2.0 * x1275 * x3
                + 2.0 * x1291 * x153
                + x1294
                - x60 * (2.0 * x1271 * x3 + 2.0 * x1288 * x153 + x1295)
            )
            + x153
            * (
                2.0 * x1285 * x3
                - 2.0 * x4 * (x1288 * x60 - x1291)
                + x49 * (2.0 * x1281 * x429 + x1292)
            )
            + x3
            * (
                2.0 * x1285 * x153
                + x3 * (x1280 + 2.0 * x1283 * x153 + 2.0 * x3 * (x1276 * x3 + x1279))
                - 2.0 * x4 * (x1271 * x60 - x1275)
            )
        )
    )
    result[1, 9] = numpy.sum(
        0.5
        * x1231
        * (
            x134
            * (
                2.0 * x1299 * x3
                + 2.0 * x1310 * x153
                + x1312
                - x60 * (2.0 * x1297 * x3 + 2.0 * x1308 * x153 + x1313)
            )
            + x153
            * (
                2.0 * x1306 * x3
                - 2.0 * x4 * (x1308 * x60 - x1310)
                + x49 * (2.0 * x1033 * x767 + x1311)
            )
            + x3
            * (
                2.0 * x1306 * x153
                + x3
                * (x1301 + 2.0 * x1303 * x153 + 2.0 * x3 * (x1033 * x752 + x1300 * x153))
                - 2.0 * x4 * (x1297 * x60 - x1299)
            )
        )
    )
    result[1, 10] = numpy.sum(
        0.5
        * x1097
        * (
            x134
            * (
                2.0 * x1316 * x3 * x49
                + 2.0 * x1333 * x3
                + x1339
                - x60 * (2.0 * x1326 * x3 + 2.0 * x1329 * x3 + x1340)
            )
            + 2.0 * x3 * x49 * (x1325 + x4 * (x1316 - x1320))
            + x3
            * (
                2.0 * x1325 * x49
                + x3 * (2.0 * x1334 * x3 + x1337 + 2.0 * x3 * (x1334 + x1336 * x3))
                - 2.0 * x4 * (x1329 * x60 - x1333)
            )
        )
    )
    result[1, 11] = numpy.sum(
        0.5
        * x1181
        * (
            x134
            * (
                2.0 * x1342 * x3 * x49
                + 2.0 * x1355 * x3
                + x1360
                - x60 * (2.0 * x1350 * x3 + 2.0 * x1352 * x3 + x1361)
            )
            + 2.0 * x3 * x49 * (x1349 + x4 * (x1342 - x1344))
            + x3
            * (
                2.0 * x1349 * x49
                + x3 * (2.0 * x1356 * x3 + x1358 + 2.0 * x3 * (x1356 + x1357 * x3))
                - 2.0 * x4 * (x1352 * x60 - x1355)
            )
        )
    )
    result[1, 12] = numpy.sum(
        0.5
        * x1266
        * (
            x134
            * (
                2.0 * x1362 * x3 * x49
                + 2.0 * x1374 * x3
                + x1379
                - x60 * (2.0 * x1369 * x3 + 2.0 * x1371 * x3 + x1380)
            )
            + 2.0 * x3 * x49 * (x1368 + x4 * (x1362 - x1364))
            + x3
            * (
                2.0 * x1368 * x49
                + x3 * (2.0 * x1375 * x3 + x1377 + 2.0 * x3 * (x1375 + x1376 * x3))
                - 2.0 * x4 * (x1371 * x60 - x1374)
            )
        )
    )
    result[1, 13] = numpy.sum(
        0.5
        * x1181
        * (
            x134
            * (
                2.0 * x1381 * x3 * x49
                + 2.0 * x1392 * x3
                + x1397
                - x60 * (2.0 * x1387 * x3 + 2.0 * x1389 * x3 + x1398)
            )
            + 2.0 * x3 * x49 * (x1386 + x4 * (x1381 - x1383))
            + x3
            * (
                2.0 * x1386 * x49
                + x3 * (2.0 * x1393 * x3 + x1395 + 2.0 * x3 * (x1393 + x1394 * x3))
                - 2.0 * x4 * (x1389 * x60 - x1392)
            )
        )
    )
    result[1, 14] = numpy.sum(
        0.5
        * x1097
        * (
            x134
            * (
                2.0 * x1033 * x954
                + 2.0 * x1404 * x3
                + x1407
                - x60 * (2.0 * x1033 * x949 + 2.0 * x1402 * x3 + x1408)
            )
            + x3
            * (
                2.0 * x1400 * x49
                + x3 * (2.0 * x1033 * x947 + x1406 + 2.0 * x3 * (x1033 * x945 + x1405))
                - 2.0 * x4 * (x1402 * x60 - x1404)
            )
            + 2.0 * x49 * (x1033 * x4 * (x3 * x920 - x923) + x1400 * x3)
        )
    )
    result[1, 15] = numpy.sum(
        0.5
        * x252
        * (
            x134 * (2.0 * x1411 * x429 + x1421 - x60 * (2.0 * x1414 * x429 + x1422))
            + x3**2 * (2.0 * x1418 * x429 + x1419 + 2.0 * x4 * (x1411 - x1415))
        )
    )
    result[1, 16] = numpy.sum(
        0.5
        * x1097
        * (
            x134 * (2.0 * x1423 * x429 + x1429 - x60 * (2.0 * x1424 * x429 + x1430))
            + x3**2 * (2.0 * x1426 * x429 + x1427 + 2.0 * x4 * (x1423 - x1425))
        )
    )
    result[1, 17] = numpy.sum(
        0.5
        * x1153
        * (
            x134 * (2.0 * x1432 * x429 + x1442 - x60 * (2.0 * x1434 * x429 + x1443))
            + x3**2 * (2.0 * x1437 * x429 + x1438 + 2.0 * x4 * (x1432 - x1435))
        )
    )
    result[1, 18] = numpy.sum(
        0.5
        * x1231
        * (
            x134 * (2.0 * x1444 * x429 + x1450 - x60 * (2.0 * x1445 * x429 + x1451))
            + x3**2 * (2.0 * x1447 * x429 + x1448 + 2.0 * x4 * (x1444 - x1446))
        )
    )
    result[1, 19] = numpy.sum(
        0.5
        * x1097
        * (
            x134 * (2.0 * x1452 * x429 + x1459 - x60 * (2.0 * x1453 * x429 + x1460))
            + x3**2 * (2.0 * x1455 * x429 + x1456 + 2.0 * x4 * (x1452 - x1454))
        )
    )
    result[1, 20] = numpy.sum(
        0.5
        * x252
        * (
            x134 * (2.0 * x1029 * x1033 + x1462 - x60 * (2.0 * x1025 * x1033 + x1463))
            + x3
            * (
                2.0 * x1033 * x4 * (x1018 * x3 - x1021)
                + x3 * (2.0 * x1023 * x1033 + x1461)
            )
        )
    )
    result[2, 0] = numpy.sum(
        x252
        * (
            x134
            * (
                x1466 * x3
                + x1475 * x91
                + x1480
                - x60 * (x1465 * x3 + x1474 * x91 + x1483)
            )
            + x3
            * (
                x1473 * x91
                + x3 * (x1464 * x3 * (x81 + x83) + x1468 + x1469 * x91)
                - x4 * (x1465 * x60 - x1466)
            )
            + x91
            * (
                x1473 * x3
                - x4 * (x1474 * x60 - x1475)
                + x98 * (x116 * (x1058 * x1478 + x110 * x1464) + x1472 * x3 + x1477)
            )
        )
    )
    result[2, 1] = numpy.sum(
        x1097
        * (
            x134
            * (
                x1486 * x3
                + x1498 * x98
                + x1503
                - x60 * (x1484 * x3 + x1496 * x98 + x1505)
            )
            + x3
            * (
                x1495 * x98
                + x3 * (x1464 * x3 * (x192 + x194) + x1488 + x1489 * x98)
                - x4 * (x1484 * x60 - x1486)
            )
            + x98
            * (
                x116 * (x1494 * x3 + x1500 + x153 * (x1464 * x221 + x1501))
                + x1495 * x3
                - x4 * (x1496 * x60 - x1498)
            )
        )
    )
    result[2, 2] = numpy.sum(
        x1097
        * (
            x134
            * (
                x1516 * x3
                + x1544 * x98
                + x1550
                - x60 * (x1511 * x3 + x1541 * x98 + x1552)
            )
            + x3
            * (
                x1538 * x98
                + x3 * (x1524 + x1529 * x98 + x3 * (x1518 + x1521))
                - x4 * (x1511 * x60 - x1516)
            )
            + x98
            * (
                x116 * (x153 * (x1534 * x3 + x1547) + x1537 * x3 + x1546)
                + x1538 * x3
                - x4 * (x1541 * x60 - x1544)
            )
        )
    )
    result[2, 3] = numpy.sum(
        x1153
        * (
            x116
            * (
                x153 * (x1464 * x402 + x1565 * x3 + x1571)
                + x1566 * x3
                - x4 * (x1567 * x60 - x1569)
            )
            + x134
            * (
                x116 * x1569
                + x1557 * x3
                + x1573
                - x60 * (x116 * x1567 + x1554 * x3 + x1575)
            )
            + x3
            * (
                x116 * x1566
                + x3 * (x116 * x1561 + x1560 + x3 * (x1464 * x380 + x1558))
                - x4 * (x1554 * x60 - x1557)
            )
        )
    )
    result[2, 4] = numpy.sum(
        x1181
        * (
            x116
            * (
                x153 * (x1592 * x3 + x1593 * x3 + x1602)
                + x1594 * x3
                - x4 * (x1597 * x60 - x1600)
            )
            + x134
            * (
                x116 * x1600
                + x1581 * x3
                + x1605
                - x60 * (x116 * x1597 + x1578 * x3 + x1607)
            )
            + x3
            * (
                x116 * x1594
                + x3 * (x116 * x1588 + x1585 + x3 * (x116 * x1583 + x1582 * x3))
                - x4 * (x1578 * x60 - x1581)
            )
        )
    )
    result[2, 5] = numpy.sum(
        x1153
        * (
            x116
            * (
                x153 * (x1638 * x3 + x1640 * x3 + x1649)
                + x1641 * x3
                - x4 * (x1644 * x60 - x1647)
            )
            + x134
            * (
                x116 * x1647
                + x1619 * x3
                + x1652
                - x60 * (x116 * x1644 + x1613 * x3 + x1654)
            )
            + x3
            * (
                x116 * x1641
                + x3 * (x116 * x1633 + x1627 + x3 * (x1621 + x1624))
                - x4 * (x1613 * x60 - x1619)
            )
        )
    )
    result[2, 6] = numpy.sum(
        x1231
        * (
            x134
            * (
                x153 * x1671
                + x1659 * x3
                + x1675
                - x60 * (x153 * x1669 + x1657 * x3 + x1677)
            )
            + x153
            * (x1667 * x3 - x4 * (x1669 * x60 - x1671) + x49 * (x1464 * x610 + x1673))
            + x3
            * (
                x153 * x1667
                + x3 * (x153 * x1663 + x1662 + x3 * (x1464 * x594 + x153 * x1660))
                - x4 * (x1657 * x60 - x1659)
            )
        )
    )
    result[2, 7] = numpy.sum(
        x1266
        * (
            x134
            * (
                x153 * x1705
                + x1687 * x3
                + x1710
                - x60 * (x153 * x1702 + x1683 * x3 + x1712)
            )
            + x153
            * (x1699 * x3 - x4 * (x1702 * x60 - x1705) + x49 * (x1694 * x429 + x1707))
            + x3
            * (
                x153 * x1699
                + x3 * (x153 * x1696 + x1693 + x3 * (x1688 * x3 + x1691))
                - x4 * (x1683 * x60 - x1687)
            )
        )
    )
    result[2, 8] = numpy.sum(
        x1266
        * (
            x134
            * (
                x153 * x1737
                + x1720 * x3
                + x1742
                - x60 * (x153 * x1734 + x1716 * x3 + x1744)
            )
            + x153
            * (x1731 * x3 - x4 * (x1734 * x60 - x1737) + x49 * (x1726 * x429 + x1739))
            + x3
            * (
                x153 * x1731
                + x3 * (x153 * x1728 + x1725 + x3 * (x1721 * x3 + x1723))
                - x4 * (x1716 * x60 - x1720)
            )
        )
    )
    result[2, 9] = numpy.sum(
        x1231
        * (
            x134
            * (
                x153 * x1780
                + x1757 * x3
                + x1785
                - x60 * (x153 * x1777 + x1752 * x3 + x1787)
            )
            + x153
            * (x1774 * x3 - x4 * (x1777 * x60 - x1780) + x49 * (x1766 * x429 + x1782))
            + x3
            * (
                x153 * x1774
                + x3 * (x153 * x1769 + x1765 + x3 * (x1759 + x1762))
                - x4 * (x1752 * x60 - x1757)
            )
        )
    )
    result[2, 10] = numpy.sum(
        x1097
        * (
            x134
            * (
                x1464 * x825
                + x1794 * x3
                + x1799
                - x60 * (x1464 * x820 + x1792 * x3 + x1801)
            )
            + x3
            * (
                x1790 * x49
                + x3 * (x1464 * x818 + x1797 + x3 * (x1464 * x816 + x1795))
                - x4 * (x1792 * x60 - x1794)
            )
            + x49 * (x1464 * x4 * (x3 * x792 - x795) + x1790 * x3)
        )
    )
    result[2, 11] = numpy.sum(
        x1181
        * (
            x134
            * (
                x1802 * x3 * x49
                + x1815 * x3
                + x1822
                - x60 * (x1810 * x3 + x1812 * x3 + x1824)
            )
            + x3 * x49 * (x1809 + x4 * (x1802 - x1804))
            + x3
            * (
                x1809 * x49
                + x3 * (x1816 * x3 + x1819 + x3 * (x1816 + x1817 * x3))
                - x4 * (x1812 * x60 - x1815)
            )
        )
    )
    result[2, 12] = numpy.sum(
        x1266
        * (
            x134
            * (
                x1825 * x3 * x49
                + x1838 * x3
                + x1845
                - x60 * (x1833 * x3 + x1835 * x3 + x1847)
            )
            + x3 * x49 * (x1832 + x4 * (x1825 - x1827))
            + x3
            * (
                x1832 * x49
                + x3 * (x1839 * x3 + x1842 + x3 * (x1839 + x1840 * x3))
                - x4 * (x1835 * x60 - x1838)
            )
        )
    )
    result[2, 13] = numpy.sum(
        x1181
        * (
            x134
            * (
                x1848 * x3 * x49
                + x1860 * x3
                + x1867
                - x60 * (x1855 * x3 + x1857 * x3 + x1869)
            )
            + x3 * x49 * (x1854 + x4 * (x1848 - x1850))
            + x3
            * (
                x1854 * x49
                + x3 * (x1861 * x3 + x1864 + x3 * (x1861 + x1862 * x3))
                - x4 * (x1857 * x60 - x1860)
            )
        )
    )
    result[2, 14] = numpy.sum(
        x1097
        * (
            x134
            * (
                x1870 * x3 * x49
                + x1886 * x3
                + x1896
                - x60 * (x1879 * x3 + x1882 * x3 + x1898)
            )
            + x3 * x49 * (x1878 + x4 * (x1870 - x1872))
            + x3
            * (
                x1878 * x49
                + x3 * (x1892 + x1893 + x3 * (x1887 + x1889))
                - x4 * (x1882 * x60 - x1886)
            )
        )
    )
    result[2, 15] = numpy.sum(
        x252
        * (
            x134 * (x1464 * x971 + x1902 - x60 * (x1464 * x967 + x1904))
            + x3 * (x1464 * x4 * (x3 * x960 - x963) + x3 * (x1464 * x965 + x1900))
        )
    )
    result[2, 16] = numpy.sum(
        x1097
        * (
            x134 * (x1905 * x429 + x1913 - x60 * (x1906 * x429 + x1915))
            + x3**2 * (x1908 * x429 + x1910 + x4 * (x1905 - x1907))
        )
    )
    result[2, 17] = numpy.sum(
        x1153
        * (
            x134 * (x1916 * x429 + x1924 - x60 * (x1917 * x429 + x1926))
            + x3**2 * (x1919 * x429 + x1921 + x4 * (x1916 - x1918))
        )
    )
    result[2, 18] = numpy.sum(
        x1231
        * (
            x134 * (x1927 * x429 + x1935 - x60 * (x1928 * x429 + x1937))
            + x3**2 * (x1930 * x429 + x1932 + x4 * (x1927 - x1929))
        )
    )
    result[2, 19] = numpy.sum(
        x1097
        * (
            x134 * (x1938 * x429 + x1946 - x60 * (x1939 * x429 + x1948))
            + x3**2 * (x1941 * x429 + x1943 + x4 * (x1938 - x1940))
        )
    )
    result[2, 20] = numpy.sum(
        x252
        * (
            x134 * (x1949 * x429 + x1959 - x60 * (x1950 * x429 + x1961))
            + x3 * (x3 * (x1954 + x1956) + x4 * (x1949 * x3 - x1952))
        )
    )
    result[3, 0] = numpy.sum(
        x431
        * (
            x3 * (x1975 * x85 + x1979 * x91 + x3 * (x1971 * x3 + 5.0 * x1974))
            + x4 * (x1963 * x3 + 5.0 * x1965 - x60 * (x1967 * x3 + 5.0 * x1969))
            + x91 * (x1979 * x3 + x1980 * x85 + x98 * (x1977 * x3 + 3.0 * x1981))
        )
    )
    result[3, 1] = numpy.sum(
        x1153
        * (
            x3 * (x1992 * x85 + x1995 * x98 + x3 * (x1989 * x3 + x1991 * x98))
            + x4 * (x1983 * x3 + x1984 * x98 - x60 * (x1986 * x3 + x1987 * x98))
            + x98 * (x116 * (x1994 * x3 + x1998) + x1995 * x3 + x1996 * x85)
        )
    )
    result[3, 2] = numpy.sum(
        x1153
        * (
            x3 * (x2008 * x85 + x2011 * x98 + x3 * (x2005 * x3 + 4.0 * x2007))
            + x4 * (x1999 * x3 + 4.0 * x2001 - x60 * (x2002 * x3 + 4.0 * x2004))
            + x98 * (x116 * (x2009 * x3 + 2.0 * x2014) + x2011 * x3 + x2012 * x85)
        )
    )
    result[3, 3] = numpy.sum(
        x2029
        * (
            x116 * (x153 * (x2025 * x3 + x2028) + x2026 * x3 + x2027 * x85)
            + x3 * (x116 * x2026 + x2024 * x85 + x3 * (x2021 * x3 + x2023))
            + x4 * (x2015 * x3 + x2017 - x60 * (x2018 * x3 + x2020))
        )
    )
    result[3, 4] = numpy.sum(
        x2042
        * (
            x116 * (x153 * (x2037 * x3 + x2041) + x2039 * x3 + x2040 * x85)
            + x3 * (x116 * x2039 + x2036 * x85 + x3 * (x116 * x2035 + x2034 * x3))
            + x4 * (x116 * x2031 + x2030 * x3 - x60 * (x116 * x2033 + x2032 * x3))
        )
    )
    result[3, 5] = numpy.sum(
        x2029
        * (
            x116 * (x153 * (x2053 * x3 + x2057) + x2055 * x3 + x2056 * x85)
            + x3 * (x116 * x2055 + x2052 * x85 + x3 * (x2049 * x3 + 3.0 * x2051))
            + x4 * (x2043 * x3 + 3.0 * x2045 - x60 * (x2046 * x3 + 3.0 * x2048))
        )
    )
    result[3, 6] = numpy.sum(
        x2069
        * (
            x153 * (x2066 * x3 + x2067 * x3 + x2068 * x85)
            + x3 * (x153 * x2067 + x2064 * x85 + x3 * (x153 * x2063 + x2062 * x3))
            + x4 * (x153 * x2059 + x2058 * x3 - x60 * (x153 * x2061 + x2060 * x3))
        )
    )
    result[3, 7] = numpy.sum(
        x2081
        * (
            x153 * (x2078 * x3 + x2079 * x3 + x2080 * x85)
            + x3 * (x153 * x2079 + x2076 * x85 + x3 * (x153 * x2075 + x2074 * x3))
            + x4 * (x153 * x2071 + x2070 * x3 - x60 * (x153 * x2073 + x2072 * x3))
        )
    )
    result[3, 8] = numpy.sum(
        x2081
        * (
            x153 * (x2093 * x3 + x2094 * x3 + x2095 * x85)
            + x3 * (x153 * x2094 + x2091 * x85 + x3 * (x2088 * x3 + x2090))
            + x4 * (x2082 * x3 + x2084 - x60 * (x2085 * x3 + x2087))
        )
    )
    result[3, 9] = numpy.sum(
        x2069
        * (
            x153 * (x2108 * x3 + x2109 * x3 + x2110 * x85)
            + x3 * (x153 * x2109 + x2106 * x85 + x3 * (x2103 * x3 + 2.0 * x2105))
            + x4 * (x2096 * x3 + 2.0 * x2098 - x60 * (x2099 * x3 + 2.0 * x2102))
        )
    )
    result[3, 10] = numpy.sum(
        x1153
        * (
            x3 * (x2119 * x3 + x2121 * x85 + x3 * (x2119 + x2120 * x3))
            + x4 * (x2112 + x2113 * x3 - x60 * (x2115 + x2116 * x3))
            + x49 * (x2117 * x429 + x2118 * x85)
        )
    )
    result[3, 11] = numpy.sum(
        x2042
        * (
            x3 * (x2130 * x3 + x2132 * x85 + x3 * (x2130 + x2131 * x3))
            + x4 * (x2123 + x2124 * x3 - x60 * (x2126 + x2127 * x3))
            + x49 * (x2128 * x429 + x2129 * x85)
        )
    )
    result[3, 12] = numpy.sum(
        x2081
        * (
            x3 * (x2141 * x3 + x2143 * x85 + x3 * (x2141 + x2142 * x3))
            + x4 * (x2134 + x2135 * x3 - x60 * (x2137 + x2138 * x3))
            + x49 * (x2139 * x429 + x2140 * x85)
        )
    )
    result[3, 13] = numpy.sum(
        x2042
        * (
            x3 * (x2152 * x3 + x2154 * x85 + x3 * (x2152 + x2153 * x3))
            + x4 * (x2145 + x2146 * x3 - x60 * (x2148 + x2149 * x3))
            + x49 * (x2150 * x429 + x2151 * x85)
        )
    )
    result[3, 14] = numpy.sum(
        x1153
        * (
            x3 * (x2164 * x3 + x2166 * x85 + x3 * (x2164 + x2165 * x3))
            + x4 * (x2156 + x2157 * x3 - x60 * (x2159 + x2161 * x3))
            + x49 * (x2162 * x429 + x2163 * x85)
        )
    )
    result[3, 15] = numpy.sum(
        x3 * x431 * (x2170 * x429 + x2171 * x85 + x4 * (x2167 - x2169))
    )
    result[3, 16] = numpy.sum(
        x1153 * x3 * (x2175 * x429 + x2176 * x85 + x4 * (x2172 - x2174))
    )
    result[3, 17] = numpy.sum(
        x2029 * x3 * (x2180 * x429 + x2181 * x85 + x4 * (x2177 - x2179))
    )
    result[3, 18] = numpy.sum(
        x2069 * x3 * (x2185 * x429 + x2186 * x85 + x4 * (x2182 - x2184))
    )
    result[3, 19] = numpy.sum(
        x1153 * x3 * (x2190 * x429 + x2191 * x85 + x4 * (x2187 - x2189))
    )
    result[3, 20] = numpy.sum(
        x3 * x431 * (x2195 * x429 + x2196 * x85 + x4 * (x2192 - x2194))
    )
    result[4, 0] = numpy.sum(
        0.5
        * x482
        * (
            2.0 * x1464 * x4 * (x1037 + x1038 * x91 - x60 * (x1034 + x1035))
            + x3 * (2.0 * x1464 * x3 * (x1040 + x1041) + x2197 + 2.0 * x2198 * x91)
            + x91 * (2.0 * x1464 * x98 * (x1049 + x1050) + 2.0 * x2198 * x3 + x2199)
        )
    )
    result[4, 1] = numpy.sum(
        0.5
        * x1181
        * (
            x3 * (x2212 + 2.0 * x2215 * x98 + 2.0 * x3 * (x2209 * x3 + x2211 * x98))
            + 2.0 * x4 * (x2201 * x3 + x2203 * x98 - x60 * (x2205 * x3 + x2207 * x98))
            + x98 * (2.0 * x116 * (x2214 * x3 + x2218) + 2.0 * x2215 * x3 + x2216)
        )
    )
    result[4, 2] = numpy.sum(
        0.5
        * x1181
        * (
            x3 * (2.0 * x1033 * x3 * (x1518 + x1521) + x2220 + 2.0 * x2221 * x98)
            + 2.0 * x4 * (x1033 * x1513 - x1033 * x60 * (x1507 + x1510) + x2219 * x98)
            + x98 * (2.0 * x1033 * x116 * (x1533 + x1536) + 2.0 * x2221 * x3 + x2222)
        )
    )
    result[4, 3] = numpy.sum(
        0.5
        * x2042
        * (
            x116 * (2.0 * x153 * (x2233 * x3 + x2236) + 2.0 * x2234 * x3 + x2235)
            + x3 * (2.0 * x116 * x2234 + x2232 + 2.0 * x3 * (x2229 * x3 + x2231))
            + 2.0 * x4 * (x2223 * x3 + x2225 - x60 * (x2226 * x3 + x2228))
        )
    )
    result[4, 4] = numpy.sum(
        0.5
        * x2249
        * (
            x116 * (2.0 * x153 * (x2244 * x3 + x2248) + 2.0 * x2246 * x3 + x2247)
            + x3 * (2.0 * x116 * x2246 + x2243 + 2.0 * x3 * (x116 * x2242 + x2241 * x3))
            + 2.0 * x4 * (x116 * x2238 + x2237 * x3 - x60 * (x116 * x2240 + x2239 * x3))
        )
    )
    result[4, 5] = numpy.sum(
        0.5
        * x2042
        * (
            x116 * (2.0 * x153 * (x1033 * x1639 + x2254) + 2.0 * x2252 * x3 + x2253)
            + x3 * (2.0 * x1033 * x3 * (x1621 + x1624) + 2.0 * x116 * x2252 + x2251)
            + 2.0 * x4 * (x1033 * x1615 - x1033 * x60 * (x1609 + x1612) + x116 * x2250)
        )
    )
    result[4, 6] = numpy.sum(
        0.5
        * x2266
        * (
            x153 * (2.0 * x2263 * x3 + 2.0 * x2264 * x3 + x2265)
            + x3 * (2.0 * x153 * x2264 + x2261 + 2.0 * x3 * (x153 * x2260 + x2259 * x3))
            + 2.0 * x4 * (x153 * x2256 + x2255 * x3 - x60 * (x153 * x2258 + x2257 * x3))
        )
    )
    result[4, 7] = numpy.sum(
        0.5
        * x2282
        * (
            x153 * (2.0 * x2279 * x3 + 2.0 * x2280 * x3 + x2281)
            + x3 * (2.0 * x153 * x2280 + x2276 + 2.0 * x3 * (x153 * x2275 + x2274 * x3))
            + 2.0 * x4 * (x153 * x2269 + x2268 * x3 - x60 * (x153 * x2272 + x2271 * x3))
        )
    )
    result[4, 8] = numpy.sum(
        0.5
        * x2282
        * (
            x153 * (2.0 * x2294 * x3 + 2.0 * x2295 * x3 + x2296)
            + x3 * (2.0 * x153 * x2295 + x2292 + 2.0 * x3 * (x2289 * x3 + x2291))
            + 2.0 * x4 * (x2283 * x3 + x2285 - x60 * (x2286 * x3 + x2288))
        )
    )
    result[4, 9] = numpy.sum(
        0.5
        * x2266
        * (
            x153 * (2.0 * x1033 * x1773 + 2.0 * x2300 * x3 + x2301)
            + x3 * (2.0 * x1033 * x3 * (x1759 + x1762) + 2.0 * x153 * x2300 + x2298)
            + 2.0 * x4 * (x1033 * x1754 - x1033 * x60 * (x1746 + x1751) + x153 * x2297)
        )
    )
    result[4, 10] = numpy.sum(
        0.5
        * x1181
        * (
            x3 * (2.0 * x2310 * x3 + x2312 + 2.0 * x3 * (x2310 + x2311 * x3))
            + 2.0 * x4 * (x2303 + x2304 * x3 - x60 * (x2306 + x2307 * x3))
            + x49 * (2.0 * x2308 * x429 + x2309)
        )
    )
    result[4, 11] = numpy.sum(
        0.5
        * x2249
        * (
            x3 * (2.0 * x2321 * x3 + x2323 + 2.0 * x3 * (x2321 + x2322 * x3))
            + 2.0 * x4 * (x2314 + x2315 * x3 - x60 * (x2317 + x2318 * x3))
            + x49 * (2.0 * x2319 * x429 + x2320)
        )
    )
    result[4, 12] = numpy.sum(
        0.5
        * x2282
        * (
            x3 * (2.0 * x2332 * x3 + x2334 + 2.0 * x3 * (x2332 + x2333 * x3))
            + 2.0 * x4 * (x2325 + x2326 * x3 - x60 * (x2328 + x2329 * x3))
            + x49 * (2.0 * x2330 * x429 + x2331)
        )
    )
    result[4, 13] = numpy.sum(
        0.5
        * x2249
        * (
            x3 * (2.0 * x2343 * x3 + x2345 + 2.0 * x3 * (x2343 + x2344 * x3))
            + 2.0 * x4 * (x2336 + x2337 * x3 - x60 * (x2339 + x2340 * x3))
            + x49 * (2.0 * x2341 * x429 + x2342)
        )
    )
    result[4, 14] = numpy.sum(
        0.5
        * x1181
        * (
            x3 * (2.0 * x1033 * x1893 + x2350 + 2.0 * x3 * (x1033 * x1889 + x2349))
            + 2.0 * x4 * (x1033 * x1885 + x2346 - x60 * (x1033 * x1881 + x2347))
            + x49 * (2.0 * x1033 * x1875 + x2348)
        )
    )
    result[4, 15] = numpy.sum(
        0.5 * x3 * x482 * (2.0 * x2354 * x429 + x2355 + 2.0 * x4 * (x2351 - x2353))
    )
    result[4, 16] = numpy.sum(
        0.5 * x1181 * x3 * (2.0 * x2359 * x429 + x2360 + 2.0 * x4 * (x2356 - x2358))
    )
    result[4, 17] = numpy.sum(
        0.5 * x2042 * x3 * (2.0 * x2367 * x429 + x2368 + 2.0 * x4 * (x2362 - x2365))
    )
    result[4, 18] = numpy.sum(
        0.5 * x2266 * x3 * (2.0 * x2372 * x429 + x2373 + 2.0 * x4 * (x2369 - x2371))
    )
    result[4, 19] = numpy.sum(
        0.5 * x1181 * x3 * (2.0 * x2377 * x429 + x2378 + 2.0 * x4 * (x2374 - x2376))
    )
    result[4, 20] = numpy.sum(
        0.5
        * x482
        * (2.0 * x1033 * x4 * (x1949 * x3 - x1952) + x3 * (2.0 * x1033 * x1954 + x2379))
    )
    result[5, 0] = numpy.sum(
        x431
        * (
            x3 * (x2394 + x2399 * x91 + x3 * (x2388 + x2391))
            + x4 * (x2381 * x3 + 5.0 * x2383 - x60 * (x2384 * x3 + 5.0 * x2386))
            + x91 * (x2399 * x3 + x2401 + x98 * (x2396 * x3 + 3.0 * x2402))
        )
    )
    result[5, 1] = numpy.sum(
        x1153
        * (
            x3 * (x2413 + x2416 * x98 + x3 * (x2409 * x3 + 4.0 * x2411))
            + x4 * (x2403 * x3 + 4.0 * x2405 - x60 * (x2406 * x3 + 4.0 * x2408))
            + x98 * (x116 * (x2414 * x3 + x2421) + x2416 * x3 + x2418)
        )
    )
    result[5, 2] = numpy.sum(
        x1153
        * (
            x3 * (x2435 + x2440 * x98 + x3 * (x2429 + x2432))
            + x4 * (x2422 * x3 + 4.0 * x2424 - x60 * (x2425 * x3 + 4.0 * x2427))
            + x98 * (x116 * (x2437 * x3 + x2444) + x2440 * x3 + x2442)
        )
    )
    result[5, 3] = numpy.sum(
        x2029
        * (
            x116 * (x153 * (x2459 * x3 + x2464) + x2461 * x3 + x2463)
            + x3 * (x116 * x2461 + x2458 + x3 * (x2453 * x3 + x2456))
            + x4 * (x2445 * x3 + x2448 - x60 * (x2449 * x3 + x2452))
        )
    )
    result[5, 4] = numpy.sum(
        x2042
        * (
            x116 * (x153 * (x2473 * x3 + x2478) + x2475 * x3 + x2477)
            + x3 * (x116 * x2475 + x2472 + x3 * (x116 * x2470 + x2469 * x3))
            + x4 * (x116 * x2466 + x2465 * x3 - x60 * (x116 * x2468 + x2467 * x3))
        )
    )
    result[5, 5] = numpy.sum(
        x2029
        * (
            x116 * (x153 * (x2496 * x3 + x2502) + x2499 * x3 + x2501)
            + x3 * (x116 * x2499 + x2494 + x3 * (x2488 + x2491))
            + x4 * (x2479 * x3 + x2482 - x60 * (x2483 * x3 + x2486))
        )
    )
    result[5, 6] = numpy.sum(
        x2069
        * (
            x153 * (x2516 * x3 + x2517 * x3 + x2519)
            + x3 * (x153 * x2517 + x2514 + x3 * (x2510 * x3 + 2.0 * x2512))
            + x4 * (x2503 * x3 + 2.0 * x2505 - x60 * (x2506 * x3 + 2.0 * x2509))
        )
    )
    result[5, 7] = numpy.sum(
        x2081
        * (
            x153 * (x2532 * x3 + x2533 * x3 + x2535)
            + x3 * (x153 * x2533 + x2530 + x3 * (x2526 * x3 + x2528))
            + x4 * (x2520 * x3 + x2522 - x60 * (x2523 * x3 + x2525))
        )
    )
    result[5, 8] = numpy.sum(
        x2081
        * (
            x153 * (x2548 * x3 + x2549 * x3 + x2551)
            + x3 * (x153 * x2549 + x2546 + x3 * (x2542 * x3 + x2544))
            + x4 * (x2536 * x3 + x2538 - x60 * (x2539 * x3 + x2541))
        )
    )
    result[5, 9] = numpy.sum(
        x2069
        * (
            x153 * (x2567 * x3 + x2569 * x3 + x2571)
            + x3 * (x153 * x2569 + x2565 + x3 * (x2559 + x2562))
            + x4 * (x2552 * x3 + 2.0 * x2554 - x60 * (x2555 * x3 + 2.0 * x2557))
        )
    )
    result[5, 10] = numpy.sum(
        x1153
        * (
            x3 * (x2582 * x3 + x2585 + x3 * (x2582 + x2583 * x3))
            + x4 * (x2573 + x2574 * x3 - x60 * (x2576 + x2578 * x3))
            + x49 * (x2579 * x429 + x2581)
        )
    )
    result[5, 11] = numpy.sum(
        x2042
        * (
            x3 * (x2595 * x3 + x2598 + x3 * (x2595 + x2596 * x3))
            + x4 * (x2587 + x2588 * x3 - x60 * (x2590 + x2591 * x3))
            + x49 * (x2592 * x429 + x2594)
        )
    )
    result[5, 12] = numpy.sum(
        x2081
        * (
            x3 * (x2608 * x3 + x2611 + x3 * (x2608 + x2609 * x3))
            + x4 * (x2600 + x2601 * x3 - x60 * (x2603 + x2604 * x3))
            + x49 * (x2605 * x429 + x2607)
        )
    )
    result[5, 13] = numpy.sum(
        x2042
        * (
            x3 * (x2621 * x3 + x2624 + x3 * (x2621 + x2622 * x3))
            + x4 * (x2613 + x2614 * x3 - x60 * (x2616 + x2617 * x3))
            + x49 * (x2618 * x429 + x2620)
        )
    )
    result[5, 14] = numpy.sum(
        x1153
        * (
            x3 * (x2639 + x2640 + x3 * (x2634 + x2636))
            + x4 * (x2626 + x2627 * x3 - x60 * (x2629 + x2630 * x3))
            + x49 * (x2631 * x429 + x2633)
        )
    )
    result[5, 15] = numpy.sum(x3 * x431 * (x2644 * x429 + x2646 + x4 * (x2641 - x2643)))
    result[5, 16] = numpy.sum(x1153 * x3 * (x2650 * x429 + x2652 + x4 * (x2647 - x2649)))
    result[5, 17] = numpy.sum(x2029 * x3 * (x2656 * x429 + x2658 + x4 * (x2653 - x2655)))
    result[5, 18] = numpy.sum(x2069 * x3 * (x2662 * x429 + x2664 + x4 * (x2659 - x2661)))
    result[5, 19] = numpy.sum(x1153 * x3 * (x2668 * x429 + x2670 + x4 * (x2665 - x2667)))
    result[5, 20] = numpy.sum(x3 * x431 * (x2675 + x2677 + x4 * (x2671 - x2673)))
    result[6, 0] = numpy.sum(
        x631
        * (
            x2681 * x85
            + x3 * (x2678 * x3 + 5.0 * x2680)
            + x91 * (x2679 * x3 + 4.0 * x2682)
        )
    )
    result[6, 1] = numpy.sum(
        x1231
        * (
            x2685 * x85
            + x3 * (x2683 * x3 + x2684 * x98)
            + x98 * (x116 * x2686 + x2684 * x3)
        )
    )
    result[6, 2] = numpy.sum(
        x1231
        * (
            x2690 * x85
            + x3 * (x2687 * x3 + 4.0 * x2689)
            + x98 * (x2688 * x3 + 3.0 * x2691)
        )
    )
    result[6, 3] = numpy.sum(
        x2069
        * (x116 * (x153 * x2696 + x2693 * x3) + x2695 * x85 + x3 * (x2692 * x3 + x2694))
    )
    result[6, 4] = numpy.sum(
        x2266
        * (x116 * (x2698 * x3 + x2700) + x2699 * x85 + x3 * (x116 * x2698 + x2697 * x3))
    )
    result[6, 5] = numpy.sum(
        x2069
        * (
            x116 * (x2702 * x3 + 2.0 * x2705)
            + x2704 * x85
            + x3 * (x2701 * x3 + 3.0 * x2703)
        )
    )
    result[6, 6] = numpy.sum(
        x2710
        * (x153 * (x2707 * x3 + x2709) + x2708 * x85 + x3 * (x153 * x2707 + x2706 * x3))
    )
    result[6, 7] = numpy.sum(
        x2715
        * (x153 * (x2712 * x3 + x2714) + x2713 * x85 + x3 * (x153 * x2712 + x2711 * x3))
    )
    result[6, 8] = numpy.sum(
        x2715 * (x153 * (x2717 * x3 + x2720) + x2719 * x85 + x3 * (x2716 * x3 + x2718))
    )
    result[6, 9] = numpy.sum(
        x2710
        * (x153 * (x2722 * x3 + x2725) + x2724 * x85 + x3 * (x2721 * x3 + 2.0 * x2723))
    )
    result[6, 10] = numpy.sum(
        x1231 * (x2727 * x3 + x2729 * x85 + x3 * (x2727 + x2728 * x3))
    )
    result[6, 11] = numpy.sum(
        x2266 * (x2731 * x3 + x2733 * x85 + x3 * (x2731 + x2732 * x3))
    )
    result[6, 12] = numpy.sum(
        x2715 * (x2735 * x3 + x2737 * x85 + x3 * (x2735 + x2736 * x3))
    )
    result[6, 13] = numpy.sum(
        x2266 * (x2739 * x3 + x2741 * x85 + x3 * (x2739 + x2740 * x3))
    )
    result[6, 14] = numpy.sum(
        x1231 * (x2743 * x3 + x2745 * x85 + x3 * (x2743 + x2744 * x3))
    )
    result[6, 15] = numpy.sum(x631 * (x2746 * x429 + x2747 * x85))
    result[6, 16] = numpy.sum(x1231 * (x2748 * x429 + x2749 * x85))
    result[6, 17] = numpy.sum(x2069 * (x2750 * x429 + x2751 * x85))
    result[6, 18] = numpy.sum(x2710 * (x2752 * x429 + x2753 * x85))
    result[6, 19] = numpy.sum(x1231 * (x2754 * x429 + x2755 * x85))
    result[6, 20] = numpy.sum(x631 * (x2756 * x429 + x2757 * x85))
    result[7, 0] = numpy.sum(
        x682
        * (
            x2761 * x85
            + x3 * (x2758 * x3 + 5.0 * x2760)
            + x91 * (x2759 * x3 + 4.0 * x2762)
        )
    )
    result[7, 1] = numpy.sum(
        x1266
        * (
            x2765 * x85
            + x3 * (x2763 * x3 + x2764 * x98)
            + x98 * (x116 * x2766 + x2764 * x3)
        )
    )
    result[7, 2] = numpy.sum(
        x1266
        * (
            x2770 * x85
            + x3 * (x2767 * x3 + 4.0 * x2769)
            + x98 * (x2768 * x3 + 3.0 * x2771)
        )
    )
    result[7, 3] = numpy.sum(
        x2081
        * (x116 * (x153 * x2776 + x2773 * x3) + x2775 * x85 + x3 * (x2772 * x3 + x2774))
    )
    result[7, 4] = numpy.sum(
        x2282
        * (x116 * (x2778 * x3 + x2780) + x2779 * x85 + x3 * (x116 * x2778 + x2777 * x3))
    )
    result[7, 5] = numpy.sum(
        x2081
        * (
            x116 * (x2782 * x3 + 2.0 * x2785)
            + x2784 * x85
            + x3 * (x2781 * x3 + 3.0 * x2783)
        )
    )
    result[7, 6] = numpy.sum(
        x2715
        * (x153 * (x2787 * x3 + x2789) + x2788 * x85 + x3 * (x153 * x2787 + x2786 * x3))
    )
    result[7, 7] = numpy.sum(
        x2794
        * (x153 * (x2791 * x3 + x2793) + x2792 * x85 + x3 * (x153 * x2791 + x2790 * x3))
    )
    result[7, 8] = numpy.sum(
        x2794 * (x153 * (x2796 * x3 + x2799) + x2798 * x85 + x3 * (x2795 * x3 + x2797))
    )
    result[7, 9] = numpy.sum(
        x2715
        * (x153 * (x2801 * x3 + x2804) + x2803 * x85 + x3 * (x2800 * x3 + 2.0 * x2802))
    )
    result[7, 10] = numpy.sum(
        x1266 * (x2806 * x3 + x2808 * x85 + x3 * (x2806 + x2807 * x3))
    )
    result[7, 11] = numpy.sum(
        x2282 * (x2810 * x3 + x2812 * x85 + x3 * (x2810 + x2811 * x3))
    )
    result[7, 12] = numpy.sum(
        x2794 * (x2814 * x3 + x2816 * x85 + x3 * (x2814 + x2815 * x3))
    )
    result[7, 13] = numpy.sum(
        x2282 * (x2818 * x3 + x2820 * x85 + x3 * (x2818 + x2819 * x3))
    )
    result[7, 14] = numpy.sum(
        x1266 * (x2822 * x3 + x2824 * x85 + x3 * (x2822 + x2823 * x3))
    )
    result[7, 15] = numpy.sum(x682 * (x2825 * x429 + x2826 * x85))
    result[7, 16] = numpy.sum(x1266 * (x2827 * x429 + x2828 * x85))
    result[7, 17] = numpy.sum(x2081 * (x2829 * x429 + x2830 * x85))
    result[7, 18] = numpy.sum(x2715 * (x2831 * x429 + x2832 * x85))
    result[7, 19] = numpy.sum(x1266 * (x2833 * x429 + x2834 * x85))
    result[7, 20] = numpy.sum(x682 * (x2835 * x429 + x2836 * x85))
    result[8, 0] = numpy.sum(
        0.5
        * x682
        * (
            2.0 * x1033 * x3 * (x2388 + x2391)
            + 2.0 * x1033 * x91 * (x2395 + x2398)
            + x2837
        )
    )
    result[8, 1] = numpy.sum(
        0.5
        * x1266
        * (
            x2840
            + 2.0 * x3 * (x2838 * x3 + x2839 * x98)
            + 2.0 * x98 * (x116 * x2841 + x2839 * x3)
        )
    )
    result[8, 2] = numpy.sum(
        0.5
        * x1266
        * (
            2.0 * x1033 * x3 * (x2429 + x2432)
            + 2.0 * x1033 * x98 * (x2436 + x2439)
            + x2842
        )
    )
    result[8, 3] = numpy.sum(
        0.5
        * x2081
        * (
            2.0 * x116 * (x153 * x2847 + x2844 * x3)
            + x2846
            + 2.0 * x3 * (x2843 * x3 + x2845)
        )
    )
    result[8, 4] = numpy.sum(
        0.5
        * x2282
        * (
            2.0 * x116 * (x2849 * x3 + x2851)
            + x2850
            + 2.0 * x3 * (x116 * x2849 + x2848 * x3)
        )
    )
    result[8, 5] = numpy.sum(
        0.5
        * x2081
        * (
            2.0 * x1033 * x116 * (x2495 + x2498)
            + 2.0 * x1033 * x3 * (x2488 + x2491)
            + x2852
        )
    )
    result[8, 6] = numpy.sum(
        0.5
        * x2715
        * (
            2.0 * x153 * (x2854 * x3 + x2856)
            + x2855
            + 2.0 * x3 * (x153 * x2854 + x2853 * x3)
        )
    )
    result[8, 7] = numpy.sum(
        0.5
        * x2794
        * (
            2.0 * x153 * (x2859 * x3 + x2864)
            + x2862
            + 2.0 * x3 * (x153 * x2859 + x2858 * x3)
        )
    )
    result[8, 8] = numpy.sum(
        0.5
        * x2794
        * (2.0 * x153 * (x2866 * x3 + x2869) + x2868 + 2.0 * x3 * (x2865 * x3 + x2867))
    )
    result[8, 9] = numpy.sum(
        0.5
        * x2715
        * (
            2.0 * x1033 * x3 * (x2559 + x2562)
            + 2.0 * x153 * (x1033 * x2568 + x2871)
            + x2870
        )
    )
    result[8, 10] = numpy.sum(
        0.5 * x1266 * (4.0 * x2873 * x3 + 2.0 * x2874 * x3**2 + x2875)
    )
    result[8, 11] = numpy.sum(
        0.5 * x2282 * (4.0 * x2877 * x3 + 2.0 * x2878 * x3**2 + x2879)
    )
    result[8, 12] = numpy.sum(
        0.5 * x2794 * (4.0 * x2881 * x3 + 2.0 * x2882 * x3**2 + x2883)
    )
    result[8, 13] = numpy.sum(
        0.5 * x2282 * (4.0 * x2885 * x3 + 2.0 * x2886 * x3**2 + x2887)
    )
    result[8, 14] = numpy.sum(
        0.5 * x1266 * (2.0 * x1033 * x2640 + x2889 + 2.0 * x3 * (x1033 * x2636 + x2888))
    )
    result[8, 15] = numpy.sum(0.5 * x682 * (2.0 * x2890 * x429 + x2891))
    result[8, 16] = numpy.sum(0.5 * x1266 * (2.0 * x2892 * x429 + x2893))
    result[8, 17] = numpy.sum(0.5 * x2081 * (2.0 * x2895 * x429 + x2898))
    result[8, 18] = numpy.sum(0.5 * x2715 * (2.0 * x2899 * x429 + x2900))
    result[8, 19] = numpy.sum(0.5 * x1266 * (2.0 * x2901 * x429 + x2902))
    result[8, 20] = numpy.sum(0.5 * x682 * (2.0 * x1033 * x2675 + x2903))
    result[9, 0] = numpy.sum(
        x631 * (x2910 + x3 * (x2905 + x2908) + x91 * (x2906 * x3 + 4.0 * x2911))
    )
    result[9, 1] = numpy.sum(
        x1231
        * (x2916 + x3 * (x2912 * x3 + 4.0 * x2914) + x98 * (x2913 * x3 + 3.0 * x2917))
    )
    result[9, 2] = numpy.sum(
        x1231 * (x2924 + x3 * (x2919 + x2922) + x98 * (x2920 * x3 + 3.0 * x2925))
    )
    result[9, 3] = numpy.sum(
        x2069 * (x116 * (x2927 * x3 + 2.0 * x2932) + x2931 + x3 * (x2926 * x3 + x2929))
    )
    result[9, 4] = numpy.sum(
        x2266 * (x116 * (x2934 * x3 + x2937) + x2936 + x3 * (x116 * x2934 + x2933 * x3))
    )
    result[9, 5] = numpy.sum(
        x2069 * (x116 * (x2940 * x3 + 2.0 * x2945) + x2944 + x3 * (x2939 + x2942))
    )
    result[9, 6] = numpy.sum(
        x2710 * (x153 * (x2947 * x3 + x2951) + x2950 + x3 * (x2946 * x3 + 2.0 * x2948))
    )
    result[9, 7] = numpy.sum(
        x2715 * (x153 * (x2953 * x3 + x2957) + x2956 + x3 * (x2952 * x3 + x2954))
    )
    result[9, 8] = numpy.sum(
        x2715 * (x153 * (x2959 * x3 + x2963) + x2962 + x3 * (x2958 * x3 + x2960))
    )
    result[9, 9] = numpy.sum(
        x2710 * (x153 * (x2966 * x3 + x2971) + x2970 + x3 * (x2965 + x2968))
    )
    result[9, 10] = numpy.sum(x1231 * (x2973 * x3 + x2976 + x3 * (x2973 + x2974 * x3)))
    result[9, 11] = numpy.sum(x2266 * (x2978 * x3 + x2981 + x3 * (x2978 + x2979 * x3)))
    result[9, 12] = numpy.sum(x2715 * (x2983 * x3 + x2986 + x3 * (x2983 + x2984 * x3)))
    result[9, 13] = numpy.sum(x2266 * (x2988 * x3 + x2991 + x3 * (x2988 + x2989 * x3)))
    result[9, 14] = numpy.sum(x1231 * (x2993 * x3 + x2997 + x3 * (x2993 + x2995)))
    result[9, 15] = numpy.sum(x631 * (x2998 * x429 + x3000))
    result[9, 16] = numpy.sum(x1231 * (x3001 * x429 + x3003))
    result[9, 17] = numpy.sum(x2069 * (x3004 * x429 + x3006))
    result[9, 18] = numpy.sum(x2710 * (x3007 * x429 + x3009))
    result[9, 19] = numpy.sum(x1231 * (x3010 * x429 + x3012))
    result[9, 20] = numpy.sum(x631 * (x3013 * x429 + x3015))
    result[10, 0] = numpy.sum(x252 * (x3 * x3016 + 5.0 * x3017))
    result[10, 1] = numpy.sum(x1097 * (x3 * x3018 + x3019 * x98))
    result[10, 2] = numpy.sum(x1097 * (x3 * x3020 + 4.0 * x3021))
    result[10, 3] = numpy.sum(x1153 * (x3 * x3022 + x3023))
    result[10, 4] = numpy.sum(x1181 * (x116 * x3025 + x3 * x3024))
    result[10, 5] = numpy.sum(x1153 * (x3 * x3026 + 3.0 * x3027))
    result[10, 6] = numpy.sum(x1231 * (x153 * x3029 + x3 * x3028))
    result[10, 7] = numpy.sum(x1266 * (x153 * x3031 + x3 * x3030))
    result[10, 8] = numpy.sum(x1266 * (x3 * x3032 + x3033))
    result[10, 9] = numpy.sum(x1231 * (x3 * x3034 + 2.0 * x3035))
    result[10, 10] = numpy.sum(x1097 * (x3 * x3037 + x3036))
    result[10, 11] = numpy.sum(x1181 * (x3 * x3039 + x3038))
    result[10, 12] = numpy.sum(x1266 * (x3 * x3041 + x3040))
    result[10, 13] = numpy.sum(x1181 * (x3 * x3043 + x3042))
    result[10, 14] = numpy.sum(x1097 * (x3 * x3045 + x3044))
    result[10, 15] = numpy.sum(x3046 * x3047)
    result[10, 16] = numpy.sum(x3048 * x3049)
    result[10, 17] = numpy.sum(x3050 * x3051)
    result[10, 18] = numpy.sum(x3052 * x3053)
    result[10, 19] = numpy.sum(x3049 * x3054)
    result[10, 20] = numpy.sum(x3047 * x3055)
    result[11, 0] = numpy.sum(x482 * (x3 * x3056 + 5.0 * x3057))
    result[11, 1] = numpy.sum(x1181 * (x3 * x3058 + x3059 * x98))
    result[11, 2] = numpy.sum(x1181 * (x3 * x3060 + 4.0 * x3061))
    result[11, 3] = numpy.sum(x2042 * (x3 * x3062 + x3063))
    result[11, 4] = numpy.sum(x2249 * (x116 * x3065 + x3 * x3064))
    result[11, 5] = numpy.sum(x2042 * (x3 * x3066 + 3.0 * x3067))
    result[11, 6] = numpy.sum(x2266 * (x153 * x3069 + x3 * x3068))
    result[11, 7] = numpy.sum(x2282 * (x153 * x3071 + x3 * x3070))
    result[11, 8] = numpy.sum(x2282 * (x3 * x3072 + x3073))
    result[11, 9] = numpy.sum(x2266 * (x3 * x3074 + 2.0 * x3075))
    result[11, 10] = numpy.sum(x1181 * (x3 * x3077 + x3076))
    result[11, 11] = numpy.sum(x2249 * (x3 * x3079 + x3078))
    result[11, 12] = numpy.sum(x2282 * (x3 * x3081 + x3080))
    result[11, 13] = numpy.sum(x2249 * (x3 * x3083 + x3082))
    result[11, 14] = numpy.sum(x1181 * (x3 * x3085 + x3084))
    result[11, 15] = numpy.sum(x3086 * x3087)
    result[11, 16] = numpy.sum(x3088 * x3089)
    result[11, 17] = numpy.sum(x3090 * x3091)
    result[11, 18] = numpy.sum(x3092 * x3093)
    result[11, 19] = numpy.sum(x3089 * x3094)
    result[11, 20] = numpy.sum(x3087 * x3095)
    result[12, 0] = numpy.sum(x682 * (x3 * x3096 + 5.0 * x3097))
    result[12, 1] = numpy.sum(x1266 * (x3 * x3098 + x3099 * x98))
    result[12, 2] = numpy.sum(x1266 * (x3 * x3100 + 4.0 * x3101))
    result[12, 3] = numpy.sum(x2081 * (x3 * x3102 + x3103))
    result[12, 4] = numpy.sum(x2282 * (x116 * x3105 + x3 * x3104))
    result[12, 5] = numpy.sum(x2081 * (x3 * x3106 + 3.0 * x3107))
    result[12, 6] = numpy.sum(x2715 * (x153 * x3109 + x3 * x3108))
    result[12, 7] = numpy.sum(x2794 * (x153 * x3111 + x3 * x3110))
    result[12, 8] = numpy.sum(x2794 * (x3 * x3112 + x3113))
    result[12, 9] = numpy.sum(x2715 * (x3 * x3114 + 2.0 * x3115))
    result[12, 10] = numpy.sum(x1266 * (x3 * x3117 + x3116))
    result[12, 11] = numpy.sum(x2282 * (x3 * x3119 + x3118))
    result[12, 12] = numpy.sum(x2794 * (x3 * x3121 + x3120))
    result[12, 13] = numpy.sum(x2282 * (x3 * x3123 + x3122))
    result[12, 14] = numpy.sum(x1266 * (x3 * x3125 + x3124))
    result[12, 15] = numpy.sum(x3126 * x3127)
    result[12, 16] = numpy.sum(x3128 * x3129)
    result[12, 17] = numpy.sum(x2081 * x3 * x3130)
    result[12, 18] = numpy.sum(x2715 * x3 * x3131)
    result[12, 19] = numpy.sum(x3129 * x3132)
    result[12, 20] = numpy.sum(x3127 * x3133)
    result[13, 0] = numpy.sum(x1033 * x482 * (x2905 + x2908))
    result[13, 1] = numpy.sum(x1181 * (x3 * x3134 + x3135 * x98))
    result[13, 2] = numpy.sum(x1033 * x1181 * (x2919 + x2922))
    result[13, 3] = numpy.sum(x2042 * (x3 * x3136 + x3137))
    result[13, 4] = numpy.sum(x2249 * (x116 * x3139 + x3 * x3138))
    result[13, 5] = numpy.sum(x1033 * x2042 * (x2939 + x2942))
    result[13, 6] = numpy.sum(x2266 * (x153 * x3141 + x3 * x3140))
    result[13, 7] = numpy.sum(x2282 * (x153 * x3144 + x3 * x3143))
    result[13, 8] = numpy.sum(x2282 * (x3 * x3145 + x3146))
    result[13, 9] = numpy.sum(x1033 * x2266 * (x2965 + x2968))
    result[13, 10] = numpy.sum(x1181 * (x3 * x3148 + x3147))
    result[13, 11] = numpy.sum(x2249 * (x3 * x3150 + x3149))
    result[13, 12] = numpy.sum(x2282 * (x3 * x3152 + x3151))
    result[13, 13] = numpy.sum(x2249 * (x3 * x3154 + x3153))
    result[13, 14] = numpy.sum(x1181 * (x1033 * x2995 + x3155))
    result[13, 15] = numpy.sum(x3087 * x3156)
    result[13, 16] = numpy.sum(x3089 * x3157)
    result[13, 17] = numpy.sum(x3091 * x3159)
    result[13, 18] = numpy.sum(x3093 * x3160)
    result[13, 19] = numpy.sum(x3089 * x3161)
    result[13, 20] = numpy.sum(x1033 * x3013 * x3087)
    result[14, 0] = numpy.sum(x252 * (x3 * x3162 + 5.0 * x3163))
    result[14, 1] = numpy.sum(x1097 * (x3 * x3164 + 4.0 * x3165))
    result[14, 2] = numpy.sum(x1097 * (x3 * x3166 + 4.0 * x3167))
    result[14, 3] = numpy.sum(x1153 * (x3 * x3168 + x3170))
    result[14, 4] = numpy.sum(x1181 * (x116 * x3172 + x3 * x3171))
    result[14, 5] = numpy.sum(x1153 * (x3 * x3173 + x3175))
    result[14, 6] = numpy.sum(x1231 * (x3 * x3176 + 2.0 * x3177))
    result[14, 7] = numpy.sum(x1266 * (x3 * x3178 + x3180))
    result[14, 8] = numpy.sum(x1266 * (x3 * x3181 + x3183))
    result[14, 9] = numpy.sum(x1231 * (x3 * x3184 + 2.0 * x3185))
    result[14, 10] = numpy.sum(x1097 * (x3 * x3187 + x3186))
    result[14, 11] = numpy.sum(x1181 * (x3 * x3189 + x3188))
    result[14, 12] = numpy.sum(x1266 * (x3 * x3191 + x3190))
    result[14, 13] = numpy.sum(x1181 * (x3 * x3193 + x3192))
    result[14, 14] = numpy.sum(x1097 * (x3 * x3195 + x3194))
    result[14, 15] = numpy.sum(x3047 * x3196)
    result[14, 16] = numpy.sum(x3049 * x3197)
    result[14, 17] = numpy.sum(x3051 * x3198)
    result[14, 18] = numpy.sum(x3053 * x3199)
    result[14, 19] = numpy.sum(x3049 * x3200)
    result[14, 20] = numpy.sum(x3047 * x3201)
    result[15, 0] = numpy.sum(x155 * (x1033 * x3016 + x149 * x2681))
    result[15, 1] = numpy.sum(x252 * (x1033 * x3018 + x149 * x2685 + x3017))
    result[15, 2] = numpy.sum(x252 * (x1033 * x3020 + x149 * x2690))
    result[15, 3] = numpy.sum(x431 * (x1033 * x3022 + x149 * x2695 + x153 * x3019))
    result[15, 4] = numpy.sum(x482 * (x1033 * x3024 + x149 * x2699 + x3021))
    result[15, 5] = numpy.sum(x431 * (x1033 * x3026 + x149 * x2704))
    result[15, 6] = numpy.sum(x631 * (x1033 * x3028 + x149 * x2708 + x3023))
    result[15, 7] = numpy.sum(x682 * (x1033 * x3030 + x149 * x2713 + x153 * x3025))
    result[15, 8] = numpy.sum(x682 * (x1033 * x3032 + x149 * x2719 + x3027))
    result[15, 9] = numpy.sum(x631 * (x1033 * x3034 + x149 * x2724))
    result[15, 10] = numpy.sum(x252 * (x1033 * x3037 + x149 * x2729 + x3029 * x98))
    result[15, 11] = numpy.sum(x482 * (x1033 * x3039 + x116 * x3031 + x149 * x2733))
    result[15, 12] = numpy.sum(x682 * (x1033 * x3041 + x149 * x2737 + x3033))
    result[15, 13] = numpy.sum(x482 * (x1033 * x3043 + x149 * x2741 + x3035))
    result[15, 14] = numpy.sum(x252 * (x1033 * x3045 + x149 * x2745))
    result[15, 15] = numpy.sum(x155 * (x1033 * x3046 + x149 * x2747 + 5.0 * x3036))
    result[15, 16] = numpy.sum(x252 * (x1033 * x3048 + x149 * x2749 + 4.0 * x3038))
    result[15, 17] = numpy.sum(x431 * (x1033 * x3050 + x149 * x2751 + 3.0 * x3040))
    result[15, 18] = numpy.sum(x631 * (x1033 * x3052 + x149 * x2753 + 2.0 * x3042))
    result[15, 19] = numpy.sum(x252 * (x1033 * x3054 + x149 * x2755 + x3044))
    result[15, 20] = numpy.sum(x155 * (x1033 * x3055 + x149 * x2757))
    result[16, 0] = numpy.sum(x252 * (x1033 * x3056 + x134 * x2761))
    result[16, 1] = numpy.sum(x1097 * (x1033 * x3058 + x134 * x2765 + x3057))
    result[16, 2] = numpy.sum(x1097 * (x1033 * x3060 + x134 * x2770))
    result[16, 3] = numpy.sum(x1153 * (x1033 * x3062 + x134 * x2775 + x153 * x3059))
    result[16, 4] = numpy.sum(x1181 * (x1033 * x3064 + x134 * x2779 + x3061))
    result[16, 5] = numpy.sum(x1153 * (x1033 * x3066 + x134 * x2784))
    result[16, 6] = numpy.sum(x1231 * (x1033 * x3068 + x134 * x2788 + x3063))
    result[16, 7] = numpy.sum(x1266 * (x1033 * x3070 + x134 * x2792 + x153 * x3065))
    result[16, 8] = numpy.sum(x1266 * (x1033 * x3072 + x134 * x2798 + x3067))
    result[16, 9] = numpy.sum(x1231 * (x1033 * x3074 + x134 * x2803))
    result[16, 10] = numpy.sum(x1097 * (x1033 * x3077 + x134 * x2808 + x3069 * x98))
    result[16, 11] = numpy.sum(x1181 * (x1033 * x3079 + x116 * x3071 + x134 * x2812))
    result[16, 12] = numpy.sum(x1266 * (x1033 * x3081 + x134 * x2816 + x3073))
    result[16, 13] = numpy.sum(x1181 * (x1033 * x3083 + x134 * x2820 + x3075))
    result[16, 14] = numpy.sum(x1097 * (x1033 * x3085 + x134 * x2824))
    result[16, 15] = numpy.sum(x252 * (x1033 * x3086 + x134 * x2826 + 5.0 * x3076))
    result[16, 16] = numpy.sum(x1097 * (x1033 * x3088 + x134 * x2828 + 4.0 * x3078))
    result[16, 17] = numpy.sum(x1153 * (x1033 * x3090 + x134 * x2830 + 3.0 * x3080))
    result[16, 18] = numpy.sum(x1231 * (x1033 * x3092 + x134 * x2832 + 2.0 * x3082))
    result[16, 19] = numpy.sum(x1097 * (x1033 * x3094 + x134 * x2834 + x3084))
    result[16, 20] = numpy.sum(x252 * (x1033 * x3095 + x134 * x2836))
    result[17, 0] = numpy.sum(x431 * (x1033 * x3096 + x2837))
    result[17, 1] = numpy.sum(x1153 * (x1033 * x3098 + x2840 + x3097))
    result[17, 2] = numpy.sum(x1153 * (x1033 * x3100 + x2842))
    result[17, 3] = numpy.sum(x2029 * (x1033 * x3102 + x153 * x3099 + x2846))
    result[17, 4] = numpy.sum(x2042 * (x1033 * x3104 + x2850 + x3101))
    result[17, 5] = numpy.sum(x2029 * (x1033 * x3106 + x2852))
    result[17, 6] = numpy.sum(x2069 * (x1033 * x3108 + x2855 + x3103))
    result[17, 7] = numpy.sum(x2081 * (x1033 * x3110 + x153 * x3105 + x2862))
    result[17, 8] = numpy.sum(x2081 * (x1033 * x3112 + x2868 + x3107))
    result[17, 9] = numpy.sum(x2069 * (x1033 * x3114 + x2870))
    result[17, 10] = numpy.sum(x1153 * (x1033 * x3117 + x2875 + x3109 * x98))
    result[17, 11] = numpy.sum(x2042 * (x1033 * x3119 + x116 * x3111 + x2879))
    result[17, 12] = numpy.sum(x2081 * (x1033 * x3121 + x2883 + x3113))
    result[17, 13] = numpy.sum(x2042 * (x1033 * x3123 + x2887 + x3115))
    result[17, 14] = numpy.sum(x1153 * (x1033 * x3125 + x2889))
    result[17, 15] = numpy.sum(x431 * (x1033 * x3126 + x2891 + 5.0 * x3116))
    result[17, 16] = numpy.sum(x1153 * (x1033 * x3128 + x2893 + 4.0 * x3118))
    result[17, 17] = numpy.sum(x2029 * (x1033 * x3130 + x2898 + 3.0 * x3120))
    result[17, 18] = numpy.sum(x2069 * (x1033 * x3131 + x2900 + 2.0 * x3122))
    result[17, 19] = numpy.sum(x1153 * (x1033 * x3132 + x2902 + x3124))
    result[17, 20] = numpy.sum(x431 * (x1033 * x3133 + x2903))
    result[18, 0] = numpy.sum(x631 * (x1962 * x2904 + x2910))
    result[18, 1] = numpy.sum(x1231 * (x1033 * x2907 + x1033 * x3134 + x2916))
    result[18, 2] = numpy.sum(x1231 * (x1962 * x2918 + x2924))
    result[18, 3] = numpy.sum(x2069 * (x1033 * x3136 + x153 * x3135 + x2931))
    result[18, 4] = numpy.sum(x2266 * (x1033 * x2921 + x1033 * x3138 + x2936))
    result[18, 5] = numpy.sum(x2069 * (x1962 * x2938 + x2944))
    result[18, 6] = numpy.sum(x2710 * (x1033 * x3140 + x2950 + x3137))
    result[18, 7] = numpy.sum(x2715 * (x1033 * x3143 + x153 * x3139 + x2956))
    result[18, 8] = numpy.sum(x2715 * (x1033 * x2941 + x1033 * x3145 + x2962))
    result[18, 9] = numpy.sum(x2710 * (x1962 * x2964 + x2970))
    result[18, 10] = numpy.sum(x1231 * (x1033 * x3148 + x2976 + x3141 * x98))
    result[18, 11] = numpy.sum(x2266 * (x1033 * x3150 + x116 * x3144 + x2981))
    result[18, 12] = numpy.sum(x2715 * (x1033 * x3152 + x2986 + x3146))
    result[18, 13] = numpy.sum(x2266 * (x1033 * x2967 + x1033 * x3154 + x2991))
    result[18, 14] = numpy.sum(x1231 * (x1962 * x2994 + x2997))
    result[18, 15] = numpy.sum(x631 * (x1033 * x3156 + x3000 + 5.0 * x3147))
    result[18, 16] = numpy.sum(x1231 * (x1033 * x3157 + x3003 + 4.0 * x3149))
    result[18, 17] = numpy.sum(x2069 * (x1033 * x3159 + x3006 + 3.0 * x3151))
    result[18, 18] = numpy.sum(x2710 * (x1033 * x3160 + x3009 + 2.0 * x3153))
    result[18, 19] = numpy.sum(x1231 * (x1033 * x3161 + x3012 + x3155))
    result[18, 20] = numpy.sum(x631 * (x1962 * x3013 + x3015))
    result[19, 0] = numpy.sum(x3162 * x3202)
    result[19, 1] = numpy.sum(x1097 * (x1033 * x3164 + x3163))
    result[19, 2] = numpy.sum(x3166 * x3203)
    result[19, 3] = numpy.sum(x1153 * (x1033 * x3168 + 2.0 * x3165))
    result[19, 4] = numpy.sum(x1181 * (x1033 * x3171 + x3167))
    result[19, 5] = numpy.sum(x1033 * x1153 * x3173)
    result[19, 6] = numpy.sum(x1231 * (x1033 * x3176 + x3170))
    result[19, 7] = numpy.sum(x1266 * (x1033 * x3178 + x3204))
    result[19, 8] = numpy.sum(x1266 * (x1033 * x3181 + x3174))
    result[19, 9] = numpy.sum(x1033 * x1231 * x3184)
    result[19, 10] = numpy.sum(x1097 * (x1033 * x3187 + 4.0 * x3177))
    result[19, 11] = numpy.sum(x1181 * (x1033 * x3189 + x116 * x3179))
    result[19, 12] = numpy.sum(x1266 * (x1033 * x3191 + x3183))
    result[19, 13] = numpy.sum(x1181 * (x1033 * x3193 + x3185))
    result[19, 14] = numpy.sum(x3195 * x3203)
    result[19, 15] = numpy.sum(x252 * (x1033 * x3196 + 5.0 * x3186))
    result[19, 16] = numpy.sum(x1097 * (x1033 * x3197 + 4.0 * x3188))
    result[19, 17] = numpy.sum(x1153 * (x1033 * x3198 + x3205))
    result[19, 18] = numpy.sum(x1231 * (x1033 * x3199 + 2.0 * x3192))
    result[19, 19] = numpy.sum(x1097 * (x1033 * x3200 + x3194))
    result[19, 20] = numpy.sum(x3201 * x3202)
    result[20, 0] = numpy.sum(x155 * (x1464 * x3162 + x149 * x2909))
    result[20, 1] = numpy.sum(x252 * (x1464 * x3164 + x149 * x2915))
    result[20, 2] = numpy.sum(x252 * (x1464 * x3166 + x149 * x2923 + x3163))
    result[20, 3] = numpy.sum(x431 * (x1464 * x3168 + x149 * x2930))
    result[20, 4] = numpy.sum(x482 * (x1464 * x3171 + x149 * x2935 + x3165))
    result[20, 5] = numpy.sum(x431 * (x1464 * x3173 + x149 * x2943 + 2.0 * x3167))
    result[20, 6] = numpy.sum(x631 * (x1464 * x3176 + x149 * x2949))
    result[20, 7] = numpy.sum(x682 * (x1464 * x3178 + x149 * x2955 + x3169))
    result[20, 8] = numpy.sum(x682 * (x1464 * x3181 + x149 * x2961 + x3204))
    result[20, 9] = numpy.sum(x631 * (x1464 * x3184 + x149 * x2969 + x3175))
    result[20, 10] = numpy.sum(x252 * (x1464 * x3187 + x149 * x2975))
    result[20, 11] = numpy.sum(x482 * (x1464 * x3189 + x149 * x2980 + x3177))
    result[20, 12] = numpy.sum(x682 * (x1464 * x3191 + x149 * x2985 + x3180))
    result[20, 13] = numpy.sum(x482 * (x116 * x3182 + x1464 * x3193 + x149 * x2990))
    result[20, 14] = numpy.sum(x252 * (x1464 * x3195 + x149 * x2996 + 4.0 * x3185))
    result[20, 15] = numpy.sum(x155 * (x1464 * x3196 + x149 * x2999))
    result[20, 16] = numpy.sum(x252 * (x1464 * x3197 + x149 * x3002 + x3186))
    result[20, 17] = numpy.sum(x431 * (x1464 * x3198 + x149 * x3005 + 2.0 * x3188))
    result[20, 18] = numpy.sum(x631 * (x1464 * x3199 + x149 * x3008 + x3205))
    result[20, 19] = numpy.sum(x252 * (x1464 * x3200 + x149 * x3011 + 4.0 * x3192))
    result[20, 20] = numpy.sum(x155 * (x1464 * x3201 + x149 * x3014 + 5.0 * x3194))
    return result


_2center2el3d = {
    (0, 0): _2center2el3d_00,
    (0, 1): _2center2el3d_01,
    (0, 2): _2center2el3d_02,
    (0, 3): _2center2el3d_03,
    (0, 4): _2center2el3d_04,
    (0, 5): _2center2el3d_05,
    (1, 0): _2center2el3d_10,
    (1, 1): _2center2el3d_11,
    (1, 2): _2center2el3d_12,
    (1, 3): _2center2el3d_13,
    (1, 4): _2center2el3d_14,
    (1, 5): _2center2el3d_15,
    (2, 0): _2center2el3d_20,
    (2, 1): _2center2el3d_21,
    (2, 2): _2center2el3d_22,
    (2, 3): _2center2el3d_23,
    (2, 4): _2center2el3d_24,
    (2, 5): _2center2el3d_25,
    (3, 0): _2center2el3d_30,
    (3, 1): _2center2el3d_31,
    (3, 2): _2center2el3d_32,
    (3, 3): _2center2el3d_33,
    (3, 4): _2center2el3d_34,
    (3, 5): _2center2el3d_35,
    (4, 0): _2center2el3d_40,
    (4, 1): _2center2el3d_41,
    (4, 2): _2center2el3d_42,
    (4, 3): _2center2el3d_43,
    (4, 4): _2center2el3d_44,
    (4, 5): _2center2el3d_45,
    (5, 0): _2center2el3d_50,
    (5, 1): _2center2el3d_51,
    (5, 2): _2center2el3d_52,
    (5, 3): _2center2el3d_53,
    (5, 4): _2center2el3d_54,
    (5, 5): _2center2el3d_55,
}
