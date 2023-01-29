import numpy

from pysisyphus.wavefunction.ints.boys import boys


_L_MAX = 4


def coulomb3d_00(ax, da, A, bx, db, B, C):
    """Cartesian (ss) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 1), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)

    # 1 item(s)
    result[0, 0] = numpy.sum(
        3.19153824321146
        * da
        * db
        * x1
        * numpy.sqrt(ax**1.5)
        * numpy.sqrt(bx**1.5)
        * boys(
            0,
            x0
            * (
                (-x1 * (ax * A[0] + bx * B[0]) + C[0]) ** 2
                + (-x1 * (ax * A[1] + bx * B[1]) + C[1]) ** 2
                + (-x1 * (ax * A[2] + bx * B[2]) + C[2]) ** 2
            ),
        )
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    return result


def coulomb3d_01(ax, da, A, bx, db, B, C):
    """Cartesian (sp) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 3), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - C[0]
    x4 = -x1 * (ax * A[1] + bx * B[1])
    x5 = -x4 - C[1]
    x6 = -x1 * (ax * A[2] + bx * B[2])
    x7 = -x6 - C[2]
    x8 = x0 * (x3**2 + x5**2 + x7**2)
    x9 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x10 = x9 * boys(1, x8)
    x11 = x9 * boys(0, x8)
    x12 = 1.01589817494786 * da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**2.5)

    # 3 item(s)
    result[0, 0] = numpy.sum(-x12 * (x10 * x3 + x11 * (x2 + B[0])))
    result[0, 1] = numpy.sum(-x12 * (x10 * x5 + x11 * (x4 + B[1])))
    result[0, 2] = numpy.sum(-x12 * (x10 * x7 + x11 * (x6 + B[2])))
    return result


def coulomb3d_02(ax, da, A, bx, db, B, C):
    """Cartesian (sd) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 6), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - C[0]
    x4 = -x1 * (ax * A[1] + bx * B[1])
    x5 = -x4 - C[1]
    x6 = -x1 * (ax * A[2] + bx * B[2])
    x7 = -x6 - C[2]
    x8 = x0 * (x3**2 + x5**2 + x7**2)
    x9 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x10 = x9 * boys(1, x8)
    x11 = x9 * boys(0, x8)
    x12 = 0.5 * (-x10 + x11) / (ax + bx)
    x13 = -x2 - B[0]
    x14 = x9 * boys(2, x8)
    x15 = 0.179587122125167 * da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**3.5)
    x16 = 6.53197264742181 * x15
    x17 = -x4 - B[1]
    x18 = -x10 * x5 + x11 * x17
    x19 = x10 * x17 - x14 * x5
    x20 = 11.3137084989848 * x15
    x21 = -x6 - B[2]
    x22 = -x10 * x7 + x11 * x21
    x23 = x10 * x21 - x14 * x7

    # 6 item(s)
    result[0, 0] = numpy.sum(
        x16 * (-2.0 * x10 * x13 * x3 + x11 * x13**2 + x12 + x14 * x3**2)
    )
    result[0, 1] = numpy.sum(x20 * (x13 * x18 - x19 * x3))
    result[0, 2] = numpy.sum(x20 * (x13 * x22 - x23 * x3))
    result[0, 3] = numpy.sum(x16 * (x12 + x17 * x18 - x19 * x5))
    result[0, 4] = numpy.sum(x20 * (x17 * x22 - x23 * x5))
    result[0, 5] = numpy.sum(x16 * (x12 + x21 * x22 - x23 * x7))
    return result


def coulomb3d_03(ax, da, A, bx, db, B, C):
    """Cartesian (sf) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 10), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = 0.5 / (ax + bx)
    x5 = -x2 - C[0]
    x6 = -x1 * (ax * A[1] + bx * B[1])
    x7 = -x6 - C[1]
    x8 = -x1 * (ax * A[2] + bx * B[2])
    x9 = -x8 - C[2]
    x10 = x0 * (x5**2 + x7**2 + x9**2)
    x11 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x12 = x11 * boys(1, x10)
    x13 = x11 * boys(0, x10)
    x14 = x4 * (-x12 + x13)
    x15 = -x12 * x5 + x13 * x3
    x16 = x11 * boys(2, x10)
    x17 = x16 * x5
    x18 = x12 * x3
    x19 = -x17 + x18
    x20 = x4 * (x12 - x16)
    x21 = x11 * boys(3, x10)
    x22 = 0.179587122125167 * da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**4.5)
    x23 = 5.84237394672177 * x22
    x24 = -x6 - B[1]
    x25 = x12 * x24
    x26 = x16 * x7
    x27 = -x12 * x7 + x13 * x24
    x28 = x4 * (-x25 + x26 + x27)
    x29 = x25 - x26
    x30 = x16 * x24 - x21 * x7
    x31 = 13.0639452948436 * x22
    x32 = -x8 - B[2]
    x33 = x12 * x32
    x34 = x16 * x9
    x35 = -x12 * x9 + x13 * x32
    x36 = x4 * (-x33 + x34 + x35)
    x37 = x33 - x34
    x38 = x16 * x32 - x21 * x9
    x39 = x14 + x24 * x27 - x29 * x7
    x40 = x20 + x24 * x29 - x30 * x7
    x41 = x24 * x35 - x37 * x7
    x42 = x24 * x37 - x38 * x7
    x43 = x14 + x32 * x35 - x37 * x9
    x44 = x20 + x32 * x37 - x38 * x9

    # 10 item(s)
    result[0, 0] = numpy.sum(
        x23
        * (
            x3 * (x14 + x15 * x3 - x19 * x5)
            + 2.0 * x4 * (x15 + x17 - x18)
            - x5 * (x19 * x3 + x20 - x5 * (x16 * x3 - x21 * x5))
        )
    )
    result[0, 1] = numpy.sum(
        x31 * (x28 + x3 * (x27 * x3 - x29 * x5) - x5 * (x29 * x3 - x30 * x5))
    )
    result[0, 2] = numpy.sum(
        x31 * (x3 * (x3 * x35 - x37 * x5) + x36 - x5 * (x3 * x37 - x38 * x5))
    )
    result[0, 3] = numpy.sum(x31 * (x3 * x39 - x40 * x5))
    result[0, 4] = numpy.sum(22.6274169979695 * x22 * (x3 * x41 - x42 * x5))
    result[0, 5] = numpy.sum(x31 * (x3 * x43 - x44 * x5))
    result[0, 6] = numpy.sum(x23 * (x24 * x39 + 2.0 * x28 - x40 * x7))
    result[0, 7] = numpy.sum(x31 * (x24 * x41 + x36 - x42 * x7))
    result[0, 8] = numpy.sum(x31 * (x24 * x43 - x44 * x7))
    result[0, 9] = numpy.sum(x23 * (x32 * x43 + 2.0 * x36 - x44 * x9))
    return result


def coulomb3d_04(ax, da, A, bx, db, B, C):
    """Cartesian (sg) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 15), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = 0.5 / (ax + bx)
    x5 = -x2 - C[0]
    x6 = -x1 * (ax * A[1] + bx * B[1])
    x7 = -x6 - C[1]
    x8 = -x1 * (ax * A[2] + bx * B[2])
    x9 = -x8 - C[2]
    x10 = x0 * (x5**2 + x7**2 + x9**2)
    x11 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x12 = x11 * boys(1, x10)
    x13 = x11 * boys(0, x10)
    x14 = x4 * (-x12 + x13)
    x15 = -x12 * x5 + x13 * x3
    x16 = x11 * boys(2, x10)
    x17 = x16 * x5
    x18 = x12 * x3
    x19 = -x17 + x18
    x20 = x14 + x15 * x3 - x19 * x5
    x21 = x4 * (x12 - x16)
    x22 = x19 * x3
    x23 = x11 * boys(3, x10)
    x24 = x23 * x5
    x25 = x16 * x3
    x26 = -x24 + x25
    x27 = x26 * x5
    x28 = x21 + x22 - x27
    x29 = 2.0 * x4
    x30 = x4 * (x16 - x23)
    x31 = x11 * boys(4, x10)
    x32 = -x21
    x33 = 0.179587122125167 * da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**5.5)
    x34 = 4.41641957979107 * x33
    x35 = -x6 - B[1]
    x36 = x12 * x35
    x37 = x16 * x7
    x38 = -x12 * x7 + x13 * x35
    x39 = x4 * (-x36 + x37 + x38)
    x40 = x36 - x37
    x41 = x3 * x38 - x40 * x5
    x42 = x3 * x40
    x43 = x23 * x7
    x44 = x16 * x35
    x45 = -x43 + x44
    x46 = x45 * x5
    x47 = x42 - x46
    x48 = x4 * (x40 + x43 - x44)
    x49 = x23 * x35 - x31 * x7
    x50 = 11.6847478934435 * x33
    x51 = -x8 - B[2]
    x52 = x12 * x51
    x53 = x16 * x9
    x54 = -x12 * x9 + x13 * x51
    x55 = x4 * (-x52 + x53 + x54)
    x56 = x52 - x53
    x57 = x3 * x54 - x5 * x56
    x58 = x3 * x56
    x59 = x23 * x9
    x60 = x16 * x51
    x61 = -x59 + x60
    x62 = x5 * x61
    x63 = x58 - x62
    x64 = x4 * (x56 + x59 - x60)
    x65 = x23 * x51 - x31 * x9
    x66 = x45 * x7
    x67 = x35 * x40
    x68 = x14 + x35 * x38 - x40 * x7
    x69 = x4 * (x32 + x66 - x67 + x68)
    x70 = x21 - x66 + x67
    x71 = x30 + x35 * x45 - x49 * x7
    x72 = 15.084944665313 * x33
    x73 = x61 * x7
    x74 = x35 * x56
    x75 = x35 * x54 - x56 * x7
    x76 = x4 * (x73 - x74 + x75)
    x77 = -x73 + x74
    x78 = x35 * x61 - x65 * x7
    x79 = 26.1278905896872 * x33
    x80 = x61 * x9
    x81 = x51 * x56
    x82 = x14 + x51 * x54 - x56 * x9
    x83 = x4 * (x32 + x80 - x81 + x82)
    x84 = x21 - x80 + x81
    x85 = x30 + x51 * x61 - x65 * x9
    x86 = x35 * x68 + 2.0 * x39 - x7 * x70
    x87 = x35 * x70 + 2.0 * x48 - x7 * x71
    x88 = x35 * x75 + x55 - x7 * x77
    x89 = x35 * x77 + x64 - x7 * x78
    x90 = x35 * x82 - x7 * x84
    x91 = x35 * x84 - x7 * x85
    x92 = x51 * x82 + 2.0 * x55 - x84 * x9
    x93 = x51 * x84 + 2.0 * x64 - x85 * x9

    # 15 item(s)
    result[0, 0] = numpy.sum(
        x34
        * (
            x3 * (x20 * x3 - x28 * x5 + x29 * (x15 + x17 - x18))
            + 3.0 * x4 * (x20 - x22 + x27 + x32)
            - x5
            * (
                x28 * x3
                + x29 * (x19 + x24 - x25)
                - x5 * (x26 * x3 + x30 - x5 * (x23 * x3 - x31 * x5))
            )
        )
    )
    result[0, 1] = numpy.sum(
        x50
        * (
            x29 * (x41 - x42 + x46)
            + x3 * (x3 * x41 + x39 - x47 * x5)
            - x5 * (x3 * x47 + x48 - x5 * (x3 * x45 - x49 * x5))
        )
    )
    result[0, 2] = numpy.sum(
        x50
        * (
            x29 * (x57 - x58 + x62)
            + x3 * (x3 * x57 - x5 * x63 + x55)
            - x5 * (x3 * x63 - x5 * (x3 * x61 - x5 * x65) + x64)
        )
    )
    result[0, 3] = numpy.sum(
        x72 * (x3 * (x3 * x68 - x5 * x70) - x5 * (x3 * x70 - x5 * x71) + x69)
    )
    result[0, 4] = numpy.sum(
        x79 * (x3 * (x3 * x75 - x5 * x77) - x5 * (x3 * x77 - x5 * x78) + x76)
    )
    result[0, 5] = numpy.sum(
        x72 * (x3 * (x3 * x82 - x5 * x84) - x5 * (x3 * x84 - x5 * x85) + x83)
    )
    result[0, 6] = numpy.sum(x50 * (x3 * x86 - x5 * x87))
    result[0, 7] = numpy.sum(x79 * (x3 * x88 - x5 * x89))
    result[0, 8] = numpy.sum(x79 * (x3 * x90 - x5 * x91))
    result[0, 9] = numpy.sum(x50 * (x3 * x92 - x5 * x93))
    result[0, 10] = numpy.sum(x34 * (x35 * x86 + 3.0 * x69 - x7 * x87))
    result[0, 11] = numpy.sum(x50 * (x35 * x88 - x7 * x89 + 2.0 * x76))
    result[0, 12] = numpy.sum(x72 * (x35 * x90 - x7 * x91 + x83))
    result[0, 13] = numpy.sum(x50 * (x35 * x92 - x7 * x93))
    result[0, 14] = numpy.sum(x34 * (x51 * x92 + 3.0 * x83 - x9 * x93))
    return result


def coulomb3d_10(ax, da, A, bx, db, B, C):
    """Cartesian (ps) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 1), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - C[0]
    x4 = -x1 * (ax * A[1] + bx * B[1])
    x5 = -x4 - C[1]
    x6 = -x1 * (ax * A[2] + bx * B[2])
    x7 = -x6 - C[2]
    x8 = x0 * (x3**2 + x5**2 + x7**2)
    x9 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x10 = x9 * boys(1, x8)
    x11 = x9 * boys(0, x8)
    x12 = 1.01589817494786 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**1.5)

    # 3 item(s)
    result[0, 0] = numpy.sum(-x12 * (x10 * x3 + x11 * (x2 + A[0])))
    result[1, 0] = numpy.sum(-x12 * (x10 * x5 + x11 * (x4 + A[1])))
    result[2, 0] = numpy.sum(-x12 * (x10 * x7 + x11 * (x6 + A[2])))
    return result


def coulomb3d_11(ax, da, A, bx, db, B, C):
    """Cartesian (pp) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 3), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - C[0]
    x4 = -x1 * (ax * A[1] + bx * B[1])
    x5 = -x4 - C[1]
    x6 = -x1 * (ax * A[2] + bx * B[2])
    x7 = -x6 - C[2]
    x8 = x0 * (x3**2 + x5**2 + x7**2)
    x9 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x10 = x9 * boys(1, x8)
    x11 = x9 * boys(0, x8)
    x12 = 0.5 * (-x10 + x11) / (ax + bx)
    x13 = -x2 - A[0]
    x14 = -x2 - B[0]
    x15 = x9 * boys(2, x8)
    x16 = 2.03179634989571 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**2.5)
    x17 = -x10 * x5
    x18 = -x4 - B[1]
    x19 = x11 * x18 + x17
    x20 = -x15 * x5
    x21 = x10 * x18 + x20
    x22 = -x10 * x7
    x23 = -x6 - B[2]
    x24 = x11 * x23 + x22
    x25 = -x15 * x7
    x26 = x10 * x23 + x25
    x27 = -x4 - A[1]
    x28 = -x6 - A[2]
    x29 = x11 * x28 + x22
    x30 = x10 * x28 + x25

    # 9 item(s)
    result[0, 0] = numpy.sum(
        -x16 * (-x12 + x13 * (x10 * x3 - x11 * x14) + x3 * (x10 * x14 - x15 * x3))
    )
    result[0, 1] = numpy.sum(x16 * (x13 * x19 - x21 * x3))
    result[0, 2] = numpy.sum(x16 * (x13 * x24 - x26 * x3))
    result[1, 0] = numpy.sum(x16 * (x14 * (x11 * x27 + x17) - x3 * (x10 * x27 + x20)))
    result[1, 1] = numpy.sum(x16 * (x12 + x19 * x27 - x21 * x5))
    result[1, 2] = numpy.sum(x16 * (x24 * x27 - x26 * x5))
    result[2, 0] = numpy.sum(x16 * (x14 * x29 - x3 * x30))
    result[2, 1] = numpy.sum(x16 * (x18 * x29 - x30 * x5))
    result[2, 2] = numpy.sum(x16 * (x12 + x24 * x28 - x26 * x7))
    return result


def coulomb3d_12(ax, da, A, bx, db, B, C):
    """Cartesian (pd) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 6), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = 0.5 / (ax + bx)
    x5 = -x2 - C[0]
    x6 = -x1 * (ax * A[1] + bx * B[1])
    x7 = -x6 - C[1]
    x8 = -x1 * (ax * A[2] + bx * B[2])
    x9 = -x8 - C[2]
    x10 = x0 * (x5**2 + x7**2 + x9**2)
    x11 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x12 = x11 * boys(1, x10)
    x13 = x11 * boys(0, x10)
    x14 = x4 * (-x12 + x13)
    x15 = -x2 - B[0]
    x16 = -x12 * x5 + x13 * x15
    x17 = x11 * boys(2, x10)
    x18 = x17 * x5
    x19 = x12 * x15
    x20 = -x18 + x19
    x21 = x4 * (x12 - x17)
    x22 = x11 * boys(3, x10)
    x23 = 0.179587122125167 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**3.5)
    x24 = 13.0639452948436 * x23
    x25 = -x6 - B[1]
    x26 = x12 * x25
    x27 = x17 * x7
    x28 = -x12 * x7
    x29 = x13 * x25 + x28
    x30 = x4 * (-x26 + x27 + x29)
    x31 = -x27
    x32 = x26 + x31
    x33 = -x22 * x7
    x34 = x17 * x25 + x33
    x35 = 22.6274169979695 * x23
    x36 = -x8 - B[2]
    x37 = x12 * x36
    x38 = x17 * x9
    x39 = -x12 * x9
    x40 = x13 * x36 + x39
    x41 = x4 * (-x37 + x38 + x40)
    x42 = -x38
    x43 = x37 + x42
    x44 = -x22 * x9
    x45 = x17 * x36 + x44
    x46 = x14 - x32 * x7
    x47 = x25 * x29 + x46
    x48 = x21 - x34 * x7
    x49 = x25 * x32 + x48
    x50 = -x43 * x7
    x51 = x25 * x40 + x50
    x52 = -x45 * x7
    x53 = x25 * x43 + x52
    x54 = x14 - x43 * x9
    x55 = x36 * x40 + x54
    x56 = x21 - x45 * x9
    x57 = x36 * x43 + x56
    x58 = -x6 - A[1]
    x59 = x12 * x58
    x60 = x13 * x58 + x28
    x61 = x31 + x59
    x62 = -x8 - A[2]
    x63 = x12 * x62
    x64 = x13 * x62 + x39
    x65 = x4 * (x38 - x63 + x64)
    x66 = x42 + x63
    x67 = x17 * x62 + x44
    x68 = x25 * x64 - x66 * x7
    x69 = x25 * x66 - x67 * x7
    x70 = x40 * x62 + x54
    x71 = x43 * x62 + x56

    # 18 item(s)
    result[0, 0] = numpy.sum(
        x24
        * (
            x3 * (x14 + x15 * x16 - x20 * x5)
            + 2.0 * x4 * (x16 + x18 - x19)
            - x5 * (x15 * x20 + x21 - x5 * (x15 * x17 - x22 * x5))
        )
    )
    result[0, 1] = numpy.sum(
        x35 * (x3 * (x15 * x29 - x32 * x5) + x30 - x5 * (x15 * x32 - x34 * x5))
    )
    result[0, 2] = numpy.sum(
        x35 * (x3 * (x15 * x40 - x43 * x5) + x41 - x5 * (x15 * x43 - x45 * x5))
    )
    result[0, 3] = numpy.sum(x24 * (x3 * x47 - x49 * x5))
    result[0, 4] = numpy.sum(x35 * (x3 * x51 - x5 * x53))
    result[0, 5] = numpy.sum(x24 * (x3 * x55 - x5 * x57))
    result[1, 0] = numpy.sum(
        x24
        * (
            x15 * (x15 * x60 - x5 * x61)
            + x4 * (x27 - x59 + x60)
            - x5 * (x15 * x61 - x5 * (x17 * x58 + x33))
        )
    )
    result[1, 1] = numpy.sum(x35 * (x15 * (x29 * x58 + x46) - x5 * (x32 * x58 + x48)))
    result[1, 2] = numpy.sum(x35 * (x15 * (x40 * x58 + x50) - x5 * (x43 * x58 + x52)))
    result[1, 3] = numpy.sum(x24 * (2.0 * x30 + x47 * x58 - x49 * x7))
    result[1, 4] = numpy.sum(x35 * (x41 + x51 * x58 - x53 * x7))
    result[1, 5] = numpy.sum(x24 * (x55 * x58 - x57 * x7))
    result[2, 0] = numpy.sum(
        x24 * (x15 * (x15 * x64 - x5 * x66) - x5 * (x15 * x66 - x5 * x67) + x65)
    )
    result[2, 1] = numpy.sum(x35 * (x15 * x68 - x5 * x69))
    result[2, 2] = numpy.sum(x35 * (x15 * x70 - x5 * x71))
    result[2, 3] = numpy.sum(x24 * (x25 * x68 + x65 - x69 * x7))
    result[2, 4] = numpy.sum(x35 * (x25 * x70 - x7 * x71))
    result[2, 5] = numpy.sum(x24 * (2.0 * x41 + x55 * x62 - x57 * x9))
    return result


def coulomb3d_13(ax, da, A, bx, db, B, C):
    """Cartesian (pf) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 10), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = -x2 - B[0]
    x5 = 0.5 / (ax + bx)
    x6 = -x2 - C[0]
    x7 = -x1 * (ax * A[1] + bx * B[1])
    x8 = -x7 - C[1]
    x9 = -x1 * (ax * A[2] + bx * B[2])
    x10 = -x9 - C[2]
    x11 = x0 * (x10**2 + x6**2 + x8**2)
    x12 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x13 = x12 * boys(1, x11)
    x14 = x12 * boys(0, x11)
    x15 = x5 * (-x13 + x14)
    x16 = -x13 * x6 + x14 * x4
    x17 = x12 * boys(2, x11)
    x18 = x17 * x6
    x19 = x13 * x4
    x20 = -x18 + x19
    x21 = x15 + x16 * x4 - x20 * x6
    x22 = x5 * (x13 - x17)
    x23 = x20 * x4
    x24 = x12 * boys(3, x11)
    x25 = x24 * x6
    x26 = x17 * x4
    x27 = -x25 + x26
    x28 = x27 * x6
    x29 = x22 + x23 - x28
    x30 = 2.0 * x5
    x31 = x5 * (x17 - x24)
    x32 = x12 * boys(4, x11)
    x33 = -x22
    x34 = 0.179587122125167 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**4.5)
    x35 = 11.6847478934435 * x34
    x36 = -x7 - B[1]
    x37 = x13 * x36
    x38 = x17 * x8
    x39 = -x13 * x8
    x40 = x14 * x36 + x39
    x41 = x5 * (-x37 + x38 + x40)
    x42 = -x38
    x43 = x37 + x42
    x44 = x4 * x40 - x43 * x6
    x45 = x4 * x43
    x46 = x24 * x8
    x47 = -x46
    x48 = x17 * x36
    x49 = x47 + x48
    x50 = x49 * x6
    x51 = x45 - x50
    x52 = x5 * (x43 + x46 - x48)
    x53 = -x32 * x8
    x54 = x24 * x36 + x53
    x55 = 26.1278905896872 * x34
    x56 = -x9 - B[2]
    x57 = x13 * x56
    x58 = x10 * x17
    x59 = -x10 * x13
    x60 = x14 * x56 + x59
    x61 = x5 * (-x57 + x58 + x60)
    x62 = -x58
    x63 = x57 + x62
    x64 = x4 * x60 - x6 * x63
    x65 = x4 * x63
    x66 = x10 * x24
    x67 = -x66
    x68 = x17 * x56
    x69 = x67 + x68
    x70 = x6 * x69
    x71 = x65 - x70
    x72 = x5 * (x63 + x66 - x68)
    x73 = -x10 * x32
    x74 = x24 * x56 + x73
    x75 = x36 * x43
    x76 = x15 - x43 * x8
    x77 = x36 * x40 + x76
    x78 = x49 * x8
    x79 = x33 + x78
    x80 = x5 * (-x75 + x77 + x79)
    x81 = x22 - x78
    x82 = x75 + x81
    x83 = x31 - x54 * x8
    x84 = x36 * x49 + x83
    x85 = x69 * x8
    x86 = x36 * x63
    x87 = -x63 * x8
    x88 = x36 * x60 + x87
    x89 = x5 * (x85 - x86 + x88)
    x90 = -x85
    x91 = x86 + x90
    x92 = -x74 * x8
    x93 = x36 * x69 + x92
    x94 = 45.2548339959391 * x34
    x95 = x56 * x63
    x96 = -x10 * x63 + x15
    x97 = x56 * x60 + x96
    x98 = x10 * x69
    x99 = x33 + x98
    x100 = x5 * (-x95 + x97 + x99)
    x101 = x22 - x98
    x102 = x101 + x95
    x103 = -x10 * x74 + x31
    x104 = x103 + x56 * x69
    x105 = 2.0 * x41 - x8 * x82
    x106 = x105 + x36 * x77
    x107 = 2.0 * x52 - x8 * x84
    x108 = x107 + x36 * x82
    x109 = x61 - x8 * x91
    x110 = x109 + x36 * x88
    x111 = x72 - x8 * x93
    x112 = x111 + x36 * x91
    x113 = -x102 * x8
    x114 = x113 + x36 * x97
    x115 = -x104 * x8
    x116 = x102 * x36 + x115
    x117 = -x10 * x102 + 2.0 * x61
    x118 = x117 + x56 * x97
    x119 = -x10 * x104 + 2.0 * x72
    x120 = x102 * x56 + x119
    x121 = -x7 - A[1]
    x122 = x121 * x13
    x123 = x121 * x14 + x39
    x124 = x122 + x42
    x125 = x123 * x4 - x124 * x6
    x126 = x124 * x4
    x127 = x121 * x17
    x128 = x127 + x47
    x129 = x128 * x6
    x130 = x126 - x129
    x131 = x121 * x43
    x132 = x121 * x40 + x76
    x133 = x131 + x81
    x134 = x121 * x63
    x135 = x121 * x60 + x87
    x136 = x134 + x90
    x137 = -x9 - A[2]
    x138 = x13 * x137
    x139 = x137 * x14 + x59
    x140 = x5 * (-x138 + x139 + x58)
    x141 = x138 + x62
    x142 = x139 * x4 - x141 * x6
    x143 = x141 * x4
    x144 = x137 * x17
    x145 = x144 + x67
    x146 = x145 * x6
    x147 = x143 - x146
    x148 = x5 * (x141 - x144 + x66)
    x149 = x137 * x24 + x73
    x150 = x145 * x8
    x151 = x141 * x36
    x152 = x139 * x36 - x141 * x8
    x153 = x5 * (x150 - x151 + x152)
    x154 = -x150 + x151
    x155 = x145 * x36 - x149 * x8
    x156 = x137 * x63
    x157 = x137 * x60 + x96
    x158 = x5 * (-x156 + x157 + x99)
    x159 = x101 + x156
    x160 = x103 + x137 * x69
    x161 = x140 + x152 * x36 - x154 * x8
    x162 = x148 + x154 * x36 - x155 * x8
    x163 = x157 * x36 - x159 * x8
    x164 = x159 * x36 - x160 * x8
    x165 = x117 + x137 * x97
    x166 = x102 * x137 + x119

    # 30 item(s)
    result[0, 0] = numpy.sum(
        x35
        * (
            x3 * (x21 * x4 - x29 * x6 + x30 * (x16 + x18 - x19))
            + 3.0 * x5 * (x21 - x23 + x28 + x33)
            - x6
            * (
                x29 * x4
                + x30 * (x20 + x25 - x26)
                - x6 * (x27 * x4 + x31 - x6 * (x24 * x4 - x32 * x6))
            )
        )
    )
    result[0, 1] = numpy.sum(
        x55
        * (
            x3 * (x4 * x44 + x41 - x51 * x6)
            + x30 * (x44 - x45 + x50)
            - x6 * (x4 * x51 + x52 - x6 * (x4 * x49 - x54 * x6))
        )
    )
    result[0, 2] = numpy.sum(
        x55
        * (
            x3 * (x4 * x64 - x6 * x71 + x61)
            + x30 * (x64 - x65 + x70)
            - x6 * (x4 * x71 - x6 * (x4 * x69 - x6 * x74) + x72)
        )
    )
    result[0, 3] = numpy.sum(
        x55 * (x3 * (x4 * x77 - x6 * x82) - x6 * (x4 * x82 - x6 * x84) + x80)
    )
    result[0, 4] = numpy.sum(
        x94 * (x3 * (x4 * x88 - x6 * x91) - x6 * (x4 * x91 - x6 * x93) + x89)
    )
    result[0, 5] = numpy.sum(
        -x55 * (-x100 + x3 * (x102 * x6 - x4 * x97) + x6 * (x102 * x4 - x104 * x6))
    )
    result[0, 6] = numpy.sum(x35 * (x106 * x3 - x108 * x6))
    result[0, 7] = numpy.sum(x55 * (x110 * x3 - x112 * x6))
    result[0, 8] = numpy.sum(x55 * (x114 * x3 - x116 * x6))
    result[0, 9] = numpy.sum(x35 * (x118 * x3 - x120 * x6))
    result[1, 0] = numpy.sum(
        x35
        * (
            x30 * (x125 - x126 + x129)
            + x4 * (x125 * x4 - x130 * x6 + x5 * (-x122 + x123 + x38))
            - x6
            * (
                x130 * x4
                + x5 * (x124 - x127 + x46)
                - x6 * (x128 * x4 - x6 * (x121 * x24 + x53))
            )
        )
    )
    result[1, 1] = numpy.sum(
        x55
        * (
            x4 * (x132 * x4 - x133 * x6)
            + x5 * (-x131 + x132 + x79)
            - x6 * (x133 * x4 - x6 * (x121 * x49 + x83))
        )
    )
    result[1, 2] = numpy.sum(
        x55
        * (
            x4 * (x135 * x4 - x136 * x6)
            + x5 * (-x134 + x135 + x85)
            - x6 * (x136 * x4 - x6 * (x121 * x69 + x92))
        )
    )
    result[1, 3] = numpy.sum(x55 * (x4 * (x105 + x121 * x77) - x6 * (x107 + x121 * x82)))
    result[1, 4] = numpy.sum(x94 * (x4 * (x109 + x121 * x88) - x6 * (x111 + x121 * x91)))
    result[1, 5] = numpy.sum(x55 * (x4 * (x113 + x121 * x97) - x6 * (x102 * x121 + x115)))
    result[1, 6] = numpy.sum(x35 * (x106 * x121 - x108 * x8 + 3.0 * x80))
    result[1, 7] = numpy.sum(x55 * (x110 * x121 - x112 * x8 + 2.0 * x89))
    result[1, 8] = numpy.sum(x55 * (x100 + x114 * x121 - x116 * x8))
    result[1, 9] = numpy.sum(x35 * (x118 * x121 - x120 * x8))
    result[2, 0] = numpy.sum(
        x35
        * (
            x30 * (x142 - x143 + x146)
            + x4 * (x140 + x142 * x4 - x147 * x6)
            - x6 * (x147 * x4 + x148 - x6 * (x145 * x4 - x149 * x6))
        )
    )
    result[2, 1] = numpy.sum(
        x55 * (x153 + x4 * (x152 * x4 - x154 * x6) - x6 * (x154 * x4 - x155 * x6))
    )
    result[2, 2] = numpy.sum(
        x55 * (x158 + x4 * (x157 * x4 - x159 * x6) - x6 * (x159 * x4 - x160 * x6))
    )
    result[2, 3] = numpy.sum(x55 * (x161 * x4 - x162 * x6))
    result[2, 4] = numpy.sum(x94 * (x163 * x4 - x164 * x6))
    result[2, 5] = numpy.sum(x55 * (x165 * x4 - x166 * x6))
    result[2, 6] = numpy.sum(x35 * (2.0 * x153 + x161 * x36 - x162 * x8))
    result[2, 7] = numpy.sum(x55 * (x158 + x163 * x36 - x164 * x8))
    result[2, 8] = numpy.sum(x55 * (x165 * x36 - x166 * x8))
    result[2, 9] = numpy.sum(x35 * (-x10 * x120 + 3.0 * x100 + x118 * x137))
    return result


def coulomb3d_14(ax, da, A, bx, db, B, C):
    """Cartesian (pg) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 15), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = -x2 - B[0]
    x5 = 0.5 / (ax + bx)
    x6 = -x2 - C[0]
    x7 = -x1 * (ax * A[1] + bx * B[1])
    x8 = -x7 - C[1]
    x9 = -x1 * (ax * A[2] + bx * B[2])
    x10 = -x9 - C[2]
    x11 = x0 * (x10**2 + x6**2 + x8**2)
    x12 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x13 = x12 * boys(1, x11)
    x14 = x12 * boys(0, x11)
    x15 = x5 * (-x13 + x14)
    x16 = -x13 * x6 + x14 * x4
    x17 = x12 * boys(2, x11)
    x18 = x17 * x6
    x19 = x13 * x4
    x20 = -x18 + x19
    x21 = x15 + x16 * x4 - x20 * x6
    x22 = x5 * (x13 - x17)
    x23 = x20 * x4
    x24 = x12 * boys(3, x11)
    x25 = x24 * x6
    x26 = x17 * x4
    x27 = -x25 + x26
    x28 = x27 * x6
    x29 = x22 + x23 - x28
    x30 = 2.0 * x5
    x31 = x21 * x4 - x29 * x6 + x30 * (x16 + x18 - x19)
    x32 = x29 * x4
    x33 = x5 * (x17 - x24)
    x34 = x27 * x4
    x35 = x12 * boys(4, x11)
    x36 = x35 * x6
    x37 = x24 * x4
    x38 = -x36 + x37
    x39 = x38 * x6
    x40 = x33 + x34 - x39
    x41 = x40 * x6
    x42 = x30 * (x20 + x25 - x26)
    x43 = x32 - x41 + x42
    x44 = -x22
    x45 = 3.0 * x5
    x46 = x5 * (x24 - x35)
    x47 = x12 * boys(5, x11)
    x48 = -x33
    x49 = 0.179587122125167 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**5.5)
    x50 = 8.83283915958214 * x49
    x51 = -x7 - B[1]
    x52 = x13 * x51
    x53 = x17 * x8
    x54 = -x13 * x8
    x55 = x14 * x51 + x54
    x56 = x5 * (-x52 + x53 + x55)
    x57 = -x53
    x58 = x52 + x57
    x59 = x4 * x55 - x58 * x6
    x60 = x4 * x58
    x61 = x24 * x8
    x62 = -x61
    x63 = x17 * x51
    x64 = x62 + x63
    x65 = x6 * x64
    x66 = x60 - x65
    x67 = x4 * x59 + x56 - x6 * x66
    x68 = x5 * (x58 + x61 - x63)
    x69 = x4 * x66
    x70 = x4 * x64
    x71 = x35 * x8
    x72 = -x71
    x73 = x24 * x51
    x74 = x72 + x73
    x75 = x6 * x74
    x76 = x70 - x75
    x77 = x6 * x76
    x78 = x68 + x69 - x77
    x79 = x5 * (x64 + x71 - x73)
    x80 = -x47 * x8
    x81 = x35 * x51 + x80
    x82 = 23.3694957868871 * x49
    x83 = -x9 - B[2]
    x84 = x13 * x83
    x85 = x10 * x17
    x86 = -x10 * x13
    x87 = x14 * x83 + x86
    x88 = x5 * (-x84 + x85 + x87)
    x89 = -x85
    x90 = x84 + x89
    x91 = x4 * x87 - x6 * x90
    x92 = x4 * x90
    x93 = x10 * x24
    x94 = -x93
    x95 = x17 * x83
    x96 = x94 + x95
    x97 = x6 * x96
    x98 = x92 - x97
    x99 = x4 * x91 - x6 * x98 + x88
    x100 = x5 * (x90 + x93 - x95)
    x101 = x4 * x98
    x102 = x4 * x96
    x103 = x10 * x35
    x104 = -x103
    x105 = x24 * x83
    x106 = x104 + x105
    x107 = x106 * x6
    x108 = x102 - x107
    x109 = x108 * x6
    x110 = x100 + x101 - x109
    x111 = x5 * (x103 - x105 + x96)
    x112 = -x10 * x47
    x113 = x112 + x35 * x83
    x114 = -x100
    x115 = x51 * x58
    x116 = x15 - x58 * x8
    x117 = x116 + x51 * x55
    x118 = x64 * x8
    x119 = x118 + x44
    x120 = x5 * (-x115 + x117 + x119)
    x121 = -x118 + x22
    x122 = x115 + x121
    x123 = x117 * x4 - x122 * x6
    x124 = x122 * x4
    x125 = x51 * x64
    x126 = x74 * x8
    x127 = -x126 + x33
    x128 = x125 + x127
    x129 = x128 * x6
    x130 = x124 - x129
    x131 = x126 + x48
    x132 = x5 * (x122 - x125 + x131)
    x133 = x46 - x8 * x81
    x134 = x133 + x51 * x74
    x135 = 30.169889330626 * x49
    x136 = x8 * x96
    x137 = x51 * x90
    x138 = -x8 * x90
    x139 = x138 + x51 * x87
    x140 = x5 * (x136 - x137 + x139)
    x141 = -x136
    x142 = x137 + x141
    x143 = x139 * x4 - x142 * x6
    x144 = x142 * x4
    x145 = x51 * x96
    x146 = x106 * x8
    x147 = -x146
    x148 = x145 + x147
    x149 = x148 * x6
    x150 = x144 - x149
    x151 = x5 * (x142 - x145 + x146)
    x152 = -x113 * x8
    x153 = x106 * x51 + x152
    x154 = 52.2557811793745 * x49
    x155 = x83 * x90
    x156 = -x10 * x90 + x15
    x157 = x156 + x83 * x87
    x158 = x10 * x96
    x159 = x158 + x44
    x160 = x5 * (-x155 + x157 + x159)
    x161 = -x158 + x22
    x162 = x155 + x161
    x163 = x157 * x4 - x162 * x6
    x164 = x162 * x4
    x165 = x83 * x96
    x166 = x10 * x106
    x167 = -x166 + x33
    x168 = x165 + x167
    x169 = x168 * x6
    x170 = x164 - x169
    x171 = x166 + x48
    x172 = x5 * (x162 - x165 + x171)
    x173 = -x10 * x113 + x46
    x174 = x106 * x83 + x173
    x175 = x122 * x51
    x176 = -x122 * x8 + 2.0 * x56
    x177 = x117 * x51 + x176
    x178 = x128 * x8
    x179 = 2.0 * x68
    x180 = x178 - x179
    x181 = x5 * (-x175 + x177 + x180)
    x182 = -x178 + x179
    x183 = x175 + x182
    x184 = -x134 * x8 + 2.0 * x79
    x185 = x128 * x51 + x184
    x186 = x142 * x51
    x187 = -x142 * x8 + x88
    x188 = x139 * x51 + x187
    x189 = x148 * x8
    x190 = x114 + x189
    x191 = x5 * (-x186 + x188 + x190)
    x192 = x100 - x189
    x193 = x186 + x192
    x194 = x111 - x153 * x8
    x195 = x148 * x51 + x194
    x196 = x168 * x8
    x197 = x162 * x51
    x198 = -x162 * x8
    x199 = x157 * x51 + x198
    x200 = x5 * (x196 - x197 + x199)
    x201 = -x196
    x202 = x197 + x201
    x203 = -x174 * x8
    x204 = x168 * x51 + x203
    x205 = x162 * x83
    x206 = -x10 * x162 + 2.0 * x88
    x207 = x157 * x83 + x206
    x208 = x10 * x168
    x209 = 2.0 * x100
    x210 = x208 - x209
    x211 = x5 * (-x205 + x207 + x210)
    x212 = -x208 + x209
    x213 = x205 + x212
    x214 = -x10 * x174 + 2.0 * x111
    x215 = x168 * x83 + x214
    x216 = 3.0 * x120 - x183 * x8
    x217 = x177 * x51 + x216
    x218 = 3.0 * x132 - x185 * x8
    x219 = x183 * x51 + x218
    x220 = 2.0 * x140 - x193 * x8
    x221 = x188 * x51 + x220
    x222 = 2.0 * x151 - x195 * x8
    x223 = x193 * x51 + x222
    x224 = x160 - x202 * x8
    x225 = x199 * x51 + x224
    x226 = x172 - x204 * x8
    x227 = x202 * x51 + x226
    x228 = -x213 * x8
    x229 = x207 * x51 + x228
    x230 = -x215 * x8
    x231 = x213 * x51 + x230
    x232 = -x10 * x213 + 3.0 * x160
    x233 = x207 * x83 + x232
    x234 = -x10 * x215 + 3.0 * x172
    x235 = x213 * x83 + x234
    x236 = -x7 - A[1]
    x237 = x13 * x236
    x238 = x14 * x236 + x54
    x239 = x237 + x57
    x240 = x238 * x4 - x239 * x6
    x241 = x239 * x4
    x242 = x17 * x236
    x243 = x242 + x62
    x244 = x243 * x6
    x245 = x241 - x244
    x246 = x240 * x4 - x245 * x6 + x5 * (-x237 + x238 + x53)
    x247 = x5 * (x239 - x242 + x61)
    x248 = x245 * x4
    x249 = x243 * x4
    x250 = x236 * x24
    x251 = x250 + x72
    x252 = x251 * x6
    x253 = x249 - x252
    x254 = x253 * x6
    x255 = x247 + x248 - x254
    x256 = x236 * x58
    x257 = x116 + x236 * x55
    x258 = x121 + x256
    x259 = x257 * x4 - x258 * x6
    x260 = x258 * x4
    x261 = x236 * x64
    x262 = x127 + x261
    x263 = x262 * x6
    x264 = x260 - x263
    x265 = x236 * x90
    x266 = x138 + x236 * x87
    x267 = x141 + x265
    x268 = x266 * x4 - x267 * x6
    x269 = x267 * x4
    x270 = x236 * x96
    x271 = x147 + x270
    x272 = x271 * x6
    x273 = x269 - x272
    x274 = x122 * x236
    x275 = x117 * x236 + x176
    x276 = x182 + x274
    x277 = x142 * x236
    x278 = x139 * x236 + x187
    x279 = x192 + x277
    x280 = x162 * x236
    x281 = x157 * x236 + x198
    x282 = x201 + x280
    x283 = -x9 - A[2]
    x284 = x13 * x283
    x285 = x14 * x283 + x86
    x286 = x5 * (-x284 + x285 + x85)
    x287 = x284 + x89
    x288 = x285 * x4 - x287 * x6
    x289 = x287 * x4
    x290 = x17 * x283
    x291 = x290 + x94
    x292 = x291 * x6
    x293 = x289 - x292
    x294 = x286 + x288 * x4 - x293 * x6
    x295 = x5 * (x287 - x290 + x93)
    x296 = x293 * x4
    x297 = x291 * x4
    x298 = x24 * x283
    x299 = x104 + x298
    x300 = x299 * x6
    x301 = x297 - x300
    x302 = x301 * x6
    x303 = x295 + x296 - x302
    x304 = x5 * (x103 + x291 - x298)
    x305 = x112 + x283 * x35
    x306 = -x295
    x307 = x291 * x8
    x308 = x287 * x51
    x309 = x285 * x51 - x287 * x8
    x310 = x5 * (x307 - x308 + x309)
    x311 = -x307 + x308
    x312 = x309 * x4 - x311 * x6
    x313 = x311 * x4
    x314 = x291 * x51
    x315 = x299 * x8
    x316 = x314 - x315
    x317 = x316 * x6
    x318 = x313 - x317
    x319 = x5 * (x311 - x314 + x315)
    x320 = x299 * x51 - x305 * x8
    x321 = x283 * x90
    x322 = x156 + x283 * x87
    x323 = x5 * (x159 - x321 + x322)
    x324 = x161 + x321
    x325 = x322 * x4 - x324 * x6
    x326 = x324 * x4
    x327 = x283 * x96
    x328 = x167 + x327
    x329 = x328 * x6
    x330 = x326 - x329
    x331 = x5 * (x171 + x324 - x327)
    x332 = x106 * x283 + x173
    x333 = x316 * x8
    x334 = x311 * x51
    x335 = x286 + x309 * x51 - x311 * x8
    x336 = x5 * (x306 + x333 - x334 + x335)
    x337 = x295 - x333 + x334
    x338 = x304 + x316 * x51 - x320 * x8
    x339 = x328 * x8
    x340 = x324 * x51
    x341 = x322 * x51 - x324 * x8
    x342 = x5 * (x339 - x340 + x341)
    x343 = -x339 + x340
    x344 = x328 * x51 - x332 * x8
    x345 = x162 * x283
    x346 = x157 * x283 + x206
    x347 = x5 * (x210 - x345 + x346)
    x348 = x212 + x345
    x349 = x168 * x283 + x214
    x350 = 2.0 * x310 + x335 * x51 - x337 * x8
    x351 = 2.0 * x319 + x337 * x51 - x338 * x8
    x352 = x323 + x341 * x51 - x343 * x8
    x353 = x331 + x343 * x51 - x344 * x8
    x354 = x346 * x51 - x348 * x8
    x355 = x348 * x51 - x349 * x8
    x356 = x207 * x283 + x232
    x357 = x213 * x283 + x234

    # 45 item(s)
    result[0, 0] = numpy.sum(
        x50
        * (
            x3 * (x31 * x4 - x43 * x6 + x45 * (x21 - x23 + x28 + x44))
            + 4.0 * x5 * (x31 - x32 + x41 - x42)
            - x6
            * (
                x4 * x43
                + x45 * (x29 - x34 + x39 + x48)
                - x6
                * (
                    x30 * (x27 + x36 - x37)
                    + x4 * x40
                    - x6 * (x38 * x4 + x46 - x6 * (x35 * x4 - x47 * x6))
                )
            )
        )
    )
    result[0, 1] = numpy.sum(
        x82
        * (
            x3 * (x30 * (x59 - x60 + x65) + x4 * x67 - x6 * x78)
            + x45 * (x67 - x68 - x69 + x77)
            - x6
            * (
                x30 * (x66 - x70 + x75)
                + x4 * x78
                - x6 * (x4 * x76 - x6 * (x4 * x74 - x6 * x81) + x79)
            )
        )
    )
    result[0, 2] = numpy.sum(
        x82
        * (
            x3 * (-x110 * x6 + x30 * (x91 - x92 + x97) + x4 * x99)
            + x45 * (-x101 + x109 + x114 + x99)
            - x6
            * (
                x110 * x4
                + x30 * (-x102 + x107 + x98)
                - x6 * (x108 * x4 + x111 - x6 * (x106 * x4 - x113 * x6))
            )
        )
    )
    result[0, 3] = numpy.sum(
        x135
        * (
            x3 * (x120 + x123 * x4 - x130 * x6)
            + x30 * (x123 - x124 + x129)
            - x6 * (x130 * x4 + x132 - x6 * (x128 * x4 - x134 * x6))
        )
    )
    result[0, 4] = numpy.sum(
        x154
        * (
            x3 * (x140 + x143 * x4 - x150 * x6)
            + x30 * (x143 - x144 + x149)
            - x6 * (x150 * x4 + x151 - x6 * (x148 * x4 - x153 * x6))
        )
    )
    result[0, 5] = numpy.sum(
        x135
        * (
            x3 * (x160 + x163 * x4 - x170 * x6)
            + x30 * (x163 - x164 + x169)
            - x6 * (x170 * x4 + x172 - x6 * (x168 * x4 - x174 * x6))
        )
    )
    result[0, 6] = numpy.sum(
        x82 * (x181 + x3 * (x177 * x4 - x183 * x6) - x6 * (x183 * x4 - x185 * x6))
    )
    result[0, 7] = numpy.sum(
        x154 * (x191 + x3 * (x188 * x4 - x193 * x6) - x6 * (x193 * x4 - x195 * x6))
    )
    result[0, 8] = numpy.sum(
        x154 * (x200 + x3 * (x199 * x4 - x202 * x6) - x6 * (x202 * x4 - x204 * x6))
    )
    result[0, 9] = numpy.sum(
        x82 * (x211 + x3 * (x207 * x4 - x213 * x6) - x6 * (x213 * x4 - x215 * x6))
    )
    result[0, 10] = numpy.sum(x50 * (x217 * x3 - x219 * x6))
    result[0, 11] = numpy.sum(x82 * (x221 * x3 - x223 * x6))
    result[0, 12] = numpy.sum(x135 * (x225 * x3 - x227 * x6))
    result[0, 13] = numpy.sum(x82 * (x229 * x3 - x231 * x6))
    result[0, 14] = numpy.sum(x50 * (x233 * x3 - x235 * x6))
    result[1, 0] = numpy.sum(
        x50
        * (
            x4 * (x246 * x4 - x255 * x6 + x30 * (x240 - x241 + x244))
            + x45 * (x246 - x247 - x248 + x254)
            - x6
            * (
                x255 * x4
                + x30 * (x245 - x249 + x252)
                - x6
                * (
                    x253 * x4
                    + x5 * (x243 - x250 + x71)
                    - x6 * (x251 * x4 - x6 * (x236 * x35 + x80))
                )
            )
        )
    )
    result[1, 1] = numpy.sum(
        x82
        * (
            x30 * (x259 - x260 + x263)
            + x4 * (x259 * x4 - x264 * x6 + x5 * (x119 - x256 + x257))
            - x6
            * (
                x264 * x4
                + x5 * (x131 + x258 - x261)
                - x6 * (x262 * x4 - x6 * (x133 + x236 * x74))
            )
        )
    )
    result[1, 2] = numpy.sum(
        x82
        * (
            x30 * (x268 - x269 + x272)
            + x4 * (x268 * x4 - x273 * x6 + x5 * (x136 - x265 + x266))
            - x6
            * (
                x273 * x4
                + x5 * (x146 + x267 - x270)
                - x6 * (x271 * x4 - x6 * (x106 * x236 + x152))
            )
        )
    )
    result[1, 3] = numpy.sum(
        x135
        * (
            x4 * (x275 * x4 - x276 * x6)
            + x5 * (x180 - x274 + x275)
            - x6 * (x276 * x4 - x6 * (x128 * x236 + x184))
        )
    )
    result[1, 4] = numpy.sum(
        x154
        * (
            x4 * (x278 * x4 - x279 * x6)
            + x5 * (x190 - x277 + x278)
            - x6 * (x279 * x4 - x6 * (x148 * x236 + x194))
        )
    )
    result[1, 5] = numpy.sum(
        x135
        * (
            x4 * (x281 * x4 - x282 * x6)
            + x5 * (x196 - x280 + x281)
            - x6 * (x282 * x4 - x6 * (x168 * x236 + x203))
        )
    )
    result[1, 6] = numpy.sum(
        x82 * (x4 * (x177 * x236 + x216) - x6 * (x183 * x236 + x218))
    )
    result[1, 7] = numpy.sum(
        x154 * (x4 * (x188 * x236 + x220) - x6 * (x193 * x236 + x222))
    )
    result[1, 8] = numpy.sum(
        x154 * (x4 * (x199 * x236 + x224) - x6 * (x202 * x236 + x226))
    )
    result[1, 9] = numpy.sum(
        x82 * (x4 * (x207 * x236 + x228) - x6 * (x213 * x236 + x230))
    )
    result[1, 10] = numpy.sum(x50 * (4.0 * x181 + x217 * x236 - x219 * x8))
    result[1, 11] = numpy.sum(x82 * (3.0 * x191 + x221 * x236 - x223 * x8))
    result[1, 12] = numpy.sum(x135 * (2.0 * x200 + x225 * x236 - x227 * x8))
    result[1, 13] = numpy.sum(x82 * (x211 + x229 * x236 - x231 * x8))
    result[1, 14] = numpy.sum(x50 * (x233 * x236 - x235 * x8))
    result[2, 0] = numpy.sum(
        x50
        * (
            x4 * (x294 * x4 + x30 * (x288 - x289 + x292) - x303 * x6)
            + x45 * (x294 - x296 + x302 + x306)
            - x6
            * (
                x30 * (x293 - x297 + x300)
                + x303 * x4
                - x6 * (x301 * x4 + x304 - x6 * (x299 * x4 - x305 * x6))
            )
        )
    )
    result[2, 1] = numpy.sum(
        x82
        * (
            x30 * (x312 - x313 + x317)
            + x4 * (x310 + x312 * x4 - x318 * x6)
            - x6 * (x318 * x4 + x319 - x6 * (x316 * x4 - x320 * x6))
        )
    )
    result[2, 2] = numpy.sum(
        x82
        * (
            x30 * (x325 - x326 + x329)
            + x4 * (x323 + x325 * x4 - x330 * x6)
            - x6 * (x330 * x4 + x331 - x6 * (x328 * x4 - x332 * x6))
        )
    )
    result[2, 3] = numpy.sum(
        x135 * (x336 + x4 * (x335 * x4 - x337 * x6) - x6 * (x337 * x4 - x338 * x6))
    )
    result[2, 4] = numpy.sum(
        x154 * (x342 + x4 * (x341 * x4 - x343 * x6) - x6 * (x343 * x4 - x344 * x6))
    )
    result[2, 5] = numpy.sum(
        x135 * (x347 + x4 * (x346 * x4 - x348 * x6) - x6 * (x348 * x4 - x349 * x6))
    )
    result[2, 6] = numpy.sum(x82 * (x350 * x4 - x351 * x6))
    result[2, 7] = numpy.sum(x154 * (x352 * x4 - x353 * x6))
    result[2, 8] = numpy.sum(x154 * (x354 * x4 - x355 * x6))
    result[2, 9] = numpy.sum(x82 * (x356 * x4 - x357 * x6))
    result[2, 10] = numpy.sum(x50 * (3.0 * x336 + x350 * x51 - x351 * x8))
    result[2, 11] = numpy.sum(x82 * (2.0 * x342 + x352 * x51 - x353 * x8))
    result[2, 12] = numpy.sum(x135 * (x347 + x354 * x51 - x355 * x8))
    result[2, 13] = numpy.sum(x82 * (x356 * x51 - x357 * x8))
    result[2, 14] = numpy.sum(x50 * (-x10 * x235 + 4.0 * x211 + x233 * x283))
    return result


def coulomb3d_20(ax, da, A, bx, db, B, C):
    """Cartesian (ds) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 1), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - C[0]
    x4 = -x1 * (ax * A[1] + bx * B[1])
    x5 = -x4 - C[1]
    x6 = -x1 * (ax * A[2] + bx * B[2])
    x7 = -x6 - C[2]
    x8 = x0 * (x3**2 + x5**2 + x7**2)
    x9 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x10 = x9 * boys(1, x8)
    x11 = x9 * boys(0, x8)
    x12 = 0.5 * (-x10 + x11) / (ax + bx)
    x13 = -x2 - A[0]
    x14 = x9 * boys(2, x8)
    x15 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**1.5)
    x16 = 6.53197264742181 * x15
    x17 = -x4 - A[1]
    x18 = -x10 * x5 + x11 * x17
    x19 = x10 * x17 - x14 * x5
    x20 = 11.3137084989848 * x15
    x21 = -x6 - A[2]
    x22 = -x10 * x7 + x11 * x21
    x23 = x10 * x21 - x14 * x7

    # 6 item(s)
    result[0, 0] = numpy.sum(
        x16 * (-2.0 * x10 * x13 * x3 + x11 * x13**2 + x12 + x14 * x3**2)
    )
    result[1, 0] = numpy.sum(x20 * (x13 * x18 - x19 * x3))
    result[2, 0] = numpy.sum(x20 * (x13 * x22 - x23 * x3))
    result[3, 0] = numpy.sum(x16 * (x12 + x17 * x18 - x19 * x5))
    result[4, 0] = numpy.sum(x20 * (x17 * x22 - x23 * x5))
    result[5, 0] = numpy.sum(x16 * (x12 + x21 * x22 - x23 * x7))
    return result


def coulomb3d_21(ax, da, A, bx, db, B, C):
    """Cartesian (dp) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 3), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = -x2 * (ax * A[0] + bx * B[0])
    x4 = -x3 - A[0]
    x5 = -x3 - C[0]
    x6 = -x2 * (ax * A[1] + bx * B[1])
    x7 = -x6 - C[1]
    x8 = -x2 * (ax * A[2] + bx * B[2])
    x9 = -x8 - C[2]
    x10 = x1 * (x5**2 + x7**2 + x9**2)
    x11 = (
        6.28318530717959
        * x2
        * numpy.exp(
            -ax * bx * x2 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x12 = x11 * boys(1, x10)
    x13 = -x12 * x5
    x14 = x11 * boys(0, x10)
    x15 = x11 * boys(2, x10)
    x16 = x15 * x5
    x17 = -x3 - B[0]
    x18 = x12 * x17
    x19 = x13 + x14 * x17
    x20 = x0 * (-x12 + x14)
    x21 = -x16 + x18
    x22 = x0 * (x12 - x15)
    x23 = x11 * boys(3, x10)
    x24 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**2.5)
    x25 = 13.0639452948436 * x24
    x26 = -x6 - B[1]
    x27 = x12 * x26
    x28 = x15 * x7
    x29 = -x12 * x7
    x30 = x14 * x26 + x29
    x31 = x0 * (-x27 + x28 + x30)
    x32 = -x28
    x33 = x27 + x32
    x34 = -x23 * x7
    x35 = x15 * x26 + x34
    x36 = -x8 - B[2]
    x37 = x12 * x36
    x38 = x15 * x9
    x39 = -x12 * x9
    x40 = x14 * x36 + x39
    x41 = x0 * (-x37 + x38 + x40)
    x42 = -x38
    x43 = x37 + x42
    x44 = -x23 * x9
    x45 = x15 * x36 + x44
    x46 = -x6 - A[1]
    x47 = x12 * x46
    x48 = x14 * x46 + x29
    x49 = x0 * (x28 - x47 + x48)
    x50 = x32 + x47
    x51 = x15 * x46 + x34
    x52 = 22.6274169979695 * x24
    x53 = x20 + x30 * x46 - x33 * x7
    x54 = x22 + x33 * x46 - x35 * x7
    x55 = x40 * x46 - x43 * x7
    x56 = x43 * x46 - x45 * x7
    x57 = -x8 - A[2]
    x58 = x12 * x57
    x59 = x14 * x57 + x39
    x60 = x0 * (x38 - x58 + x59)
    x61 = x42 + x58
    x62 = x15 * x57 + x44
    x63 = -x61 * x7
    x64 = x26 * x59 + x63
    x65 = -x62 * x7
    x66 = x26 * x61 + x65
    x67 = x20 + x40 * x57 - x43 * x9
    x68 = x22 + x43 * x57 - x45 * x9
    x69 = x20 + x57 * x59 - x61 * x9
    x70 = x22 + x57 * x61 - x62 * x9

    # 18 item(s)
    result[0, 0] = numpy.sum(
        x25
        * (
            x0 * (x16 - x18 + x19)
            + x0 * (-x12 * x4 + x13 + x14 * x4 + x16)
            + x4 * (x19 * x4 + x20 - x21 * x5)
            - x5 * (x21 * x4 + x22 - x5 * (x15 * x17 - x23 * x5))
        )
    )
    result[0, 1] = numpy.sum(
        x25 * (x31 + x4 * (x30 * x4 - x33 * x5) - x5 * (x33 * x4 - x35 * x5))
    )
    result[0, 2] = numpy.sum(
        x25 * (x4 * (x4 * x40 - x43 * x5) + x41 - x5 * (x4 * x43 - x45 * x5))
    )
    result[1, 0] = numpy.sum(
        x52 * (x4 * (x17 * x48 - x5 * x50) + x49 - x5 * (x17 * x50 - x5 * x51))
    )
    result[1, 1] = numpy.sum(x52 * (x4 * x53 - x5 * x54))
    result[1, 2] = numpy.sum(x52 * (x4 * x55 - x5 * x56))
    result[2, 0] = numpy.sum(
        x52 * (x4 * (x17 * x59 - x5 * x61) - x5 * (x17 * x61 - x5 * x62) + x60)
    )
    result[2, 1] = numpy.sum(x52 * (x4 * x64 - x5 * x66))
    result[2, 2] = numpy.sum(x52 * (x4 * x67 - x5 * x68))
    result[3, 0] = numpy.sum(
        x25 * (x17 * (x20 + x46 * x48 - x50 * x7) - x5 * (x22 + x46 * x50 - x51 * x7))
    )
    result[3, 1] = numpy.sum(x25 * (x31 + x46 * x53 + x49 - x54 * x7))
    result[3, 2] = numpy.sum(x25 * (x41 + x46 * x55 - x56 * x7))
    result[4, 0] = numpy.sum(x52 * (x17 * (x46 * x59 + x63) - x5 * (x46 * x61 + x65)))
    result[4, 1] = numpy.sum(x52 * (x46 * x64 + x60 - x66 * x7))
    result[4, 2] = numpy.sum(x52 * (x46 * x67 - x68 * x7))
    result[5, 0] = numpy.sum(x25 * (x17 * x69 - x5 * x70))
    result[5, 1] = numpy.sum(x25 * (x26 * x69 - x7 * x70))
    result[5, 2] = numpy.sum(x25 * (x41 + x57 * x67 + x60 - x68 * x9))
    return result


def coulomb3d_22(ax, da, A, bx, db, B, C):
    """Cartesian (dd) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 6), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = -x2 * (ax * A[0] + bx * B[0])
    x4 = -x3 - B[0]
    x5 = -x3 - C[0]
    x6 = -x2 * (ax * A[1] + bx * B[1])
    x7 = -x6 - C[1]
    x8 = -x2 * (ax * A[2] + bx * B[2])
    x9 = -x8 - C[2]
    x10 = x1 * (x5**2 + x7**2 + x9**2)
    x11 = (
        6.28318530717959
        * x2
        * numpy.exp(
            -ax * bx * x2 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x12 = x11 * boys(2, x10)
    x13 = x12 * x5
    x14 = x11 * boys(1, x10)
    x15 = x14 * x4
    x16 = -x13 + x15
    x17 = x16 * x4
    x18 = x11 * boys(0, x10)
    x19 = -x14 * x5 + x18 * x4
    x20 = x0 * (-x14 + x18)
    x21 = -x16 * x5 + x20
    x22 = x19 * x4 + x21
    x23 = x11 * boys(3, x10)
    x24 = x23 * x5
    x25 = x12 * x4
    x26 = -x24 + x25
    x27 = x26 * x5
    x28 = x0 * (-x12 + x14)
    x29 = -x28
    x30 = x27 + x29
    x31 = -x3 - A[0]
    x32 = x17 - x27 + x28
    x33 = 2.0 * x0
    x34 = x0 * (x12 - x23)
    x35 = x11 * boys(4, x10)
    x36 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**3.5)
    x37 = 15.084944665313 * x36
    x38 = -x14 * x7
    x39 = -x6 - B[1]
    x40 = x18 * x39 + x38
    x41 = x23 * x7
    x42 = -x41
    x43 = x12 * x39
    x44 = x42 + x43
    x45 = x44 * x5
    x46 = x12 * x7
    x47 = -x46
    x48 = x14 * x39
    x49 = x47 + x48
    x50 = -x49 * x5
    x51 = x4 * x49
    x52 = x4 * x40 + x50
    x53 = x0 * (x40 + x46 - x48)
    x54 = -x45 + x51
    x55 = x0 * (x41 - x43 + x49)
    x56 = -x35 * x7
    x57 = x23 * x39 + x56
    x58 = 26.1278905896872 * x36
    x59 = -x14 * x9
    x60 = -x8 - B[2]
    x61 = x18 * x60 + x59
    x62 = x23 * x9
    x63 = -x62
    x64 = x12 * x60
    x65 = x63 + x64
    x66 = x5 * x65
    x67 = x12 * x9
    x68 = -x67
    x69 = x14 * x60
    x70 = x68 + x69
    x71 = -x5 * x70
    x72 = x4 * x70
    x73 = x4 * x61 + x71
    x74 = x0 * (x61 + x67 - x69)
    x75 = -x66 + x72
    x76 = x0 * (x62 - x64 + x70)
    x77 = -x35 * x9
    x78 = x23 * x60 + x77
    x79 = x39 * x49
    x80 = x20 - x49 * x7
    x81 = x39 * x40 + x80
    x82 = x44 * x7
    x83 = x29 + x82
    x84 = x0 * (-x79 + x81 + x83)
    x85 = x28 - x82
    x86 = x79 + x85
    x87 = x34 - x57 * x7
    x88 = x39 * x44 + x87
    x89 = x65 * x7
    x90 = x39 * x70
    x91 = -x7 * x70
    x92 = x39 * x61 + x91
    x93 = x0 * (x89 - x90 + x92)
    x94 = -x89
    x95 = x90 + x94
    x96 = -x7 * x78
    x97 = x39 * x65 + x96
    x98 = x60 * x70
    x99 = x20 - x70 * x9
    x100 = x60 * x61 + x99
    x101 = x65 * x9
    x102 = x101 + x29
    x103 = x0 * (x100 + x102 - x98)
    x104 = -x101 + x28
    x105 = x104 + x98
    x106 = x34 - x78 * x9
    x107 = x106 + x60 * x65
    x108 = -x6 - A[1]
    x109 = x108 * x14
    x110 = x108 * x18 + x38
    x111 = x0 * (-x109 + x110 + x46)
    x112 = x109 + x47
    x113 = x110 * x4 - x112 * x5
    x114 = x112 * x4
    x115 = x108 * x12
    x116 = x115 + x42
    x117 = x116 * x5
    x118 = x114 - x117
    x119 = x0 * (x112 - x115 + x41)
    x120 = x108 * x23 + x56
    x121 = x108 * x49
    x122 = x108 * x40 + x80
    x123 = x0 * (-x121 + x122 + x83)
    x124 = x121 + x85
    x125 = x108 * x44 + x87
    x126 = 45.2548339959391 * x36
    x127 = x108 * x70
    x128 = x108 * x61 + x91
    x129 = x0 * (-x127 + x128 + x89)
    x130 = x127 + x94
    x131 = x108 * x65 + x96
    x132 = x108 * x81 + 2.0 * x53 - x7 * x86
    x133 = x108 * x86 + 2.0 * x55 - x7 * x88
    x134 = x108 * x92 - x7 * x95 + x74
    x135 = x108 * x95 - x7 * x97 + x76
    x136 = x100 * x108 - x105 * x7
    x137 = x105 * x108 - x107 * x7
    x138 = -x8 - A[2]
    x139 = x138 * x14
    x140 = x138 * x18 + x59
    x141 = x0 * (-x139 + x140 + x67)
    x142 = x139 + x68
    x143 = x140 * x4 - x142 * x5
    x144 = x142 * x4
    x145 = x12 * x138
    x146 = x145 + x63
    x147 = x146 * x5
    x148 = x144 - x147
    x149 = x0 * (x142 - x145 + x62)
    x150 = x138 * x23 + x77
    x151 = x146 * x7
    x152 = x142 * x39
    x153 = -x142 * x7
    x154 = x140 * x39 + x153
    x155 = x0 * (x151 - x152 + x154)
    x156 = -x151
    x157 = x152 + x156
    x158 = -x150 * x7
    x159 = x146 * x39 + x158
    x160 = x138 * x70
    x161 = x138 * x61 + x99
    x162 = x0 * (x102 - x160 + x161)
    x163 = x104 + x160
    x164 = x106 + x138 * x65
    x165 = x141 - x157 * x7
    x166 = x154 * x39 + x165
    x167 = x149 - x159 * x7
    x168 = x157 * x39 + x167
    x169 = -x163 * x7
    x170 = x161 * x39 + x169
    x171 = -x164 * x7
    x172 = x163 * x39 + x171
    x173 = x100 * x138 - x105 * x9 + 2.0 * x74
    x174 = x105 * x138 - x107 * x9 + 2.0 * x76
    x175 = x116 * x7
    x176 = x108 * x112
    x177 = x108 * x110 - x112 * x7 + x20
    x178 = -x175 + x176 + x28
    x179 = x108 * x142
    x180 = x108 * x140 + x153
    x181 = x156 + x179
    x182 = x146 * x9
    x183 = x138 * x142
    x184 = x138 * x140 - x142 * x9 + x20
    x185 = x0 * (x182 - x183 + x184 + x29)
    x186 = -x182 + x183 + x28
    x187 = x138 * x146 - x150 * x9 + x34
    x188 = x184 * x39 - x186 * x7
    x189 = x186 * x39 - x187 * x7
    x190 = x138 * x161 + x141 - x163 * x9 + x74
    x191 = x138 * x163 + x149 - x164 * x9 + x76

    # 36 item(s)
    result[0, 0] = numpy.sum(
        x37
        * (
            x0 * (-x17 + x22 + x30)
            + x31 * (x22 * x31 - x32 * x5 + x33 * (x13 - x15 + x19))
            + x33 * (-x16 * x31 + x19 * x31 + x21 + x30)
            - x5
            * (
                x31 * x32
                + x33 * (x16 + x24 - x25)
                - x5 * (x26 * x4 + x34 - x5 * (x23 * x4 - x35 * x5))
            )
        )
    )
    result[0, 1] = numpy.sum(
        x58
        * (
            x0 * (x45 - x51 + x52)
            + x0 * (x31 * x40 - x31 * x49 + x45 + x50)
            + x31 * (x31 * x52 - x5 * x54 + x53)
            - x5 * (x31 * x54 - x5 * (x4 * x44 - x5 * x57) + x55)
        )
    )
    result[0, 2] = numpy.sum(
        x58
        * (
            x0 * (x66 - x72 + x73)
            + x0 * (x31 * x61 - x31 * x70 + x66 + x71)
            + x31 * (x31 * x73 - x5 * x75 + x74)
            - x5 * (x31 * x75 - x5 * (x4 * x65 - x5 * x78) + x76)
        )
    )
    result[0, 3] = numpy.sum(
        x37 * (x31 * (x31 * x81 - x5 * x86) - x5 * (x31 * x86 - x5 * x88) + x84)
    )
    result[0, 4] = numpy.sum(
        x58 * (x31 * (x31 * x92 - x5 * x95) - x5 * (x31 * x95 - x5 * x97) + x93)
    )
    result[0, 5] = numpy.sum(
        x37 * (x103 + x31 * (x100 * x31 - x105 * x5) - x5 * (x105 * x31 - x107 * x5))
    )
    result[1, 0] = numpy.sum(
        x58
        * (
            x31 * (x111 + x113 * x4 - x118 * x5)
            + x33 * (x113 - x114 + x117)
            - x5 * (x118 * x4 + x119 - x5 * (x116 * x4 - x120 * x5))
        )
    )
    result[1, 1] = numpy.sum(
        x126 * (x123 + x31 * (x122 * x4 - x124 * x5) - x5 * (x124 * x4 - x125 * x5))
    )
    result[1, 2] = numpy.sum(
        x126 * (x129 + x31 * (x128 * x4 - x130 * x5) - x5 * (x130 * x4 - x131 * x5))
    )
    result[1, 3] = numpy.sum(x58 * (x132 * x31 - x133 * x5))
    result[1, 4] = numpy.sum(x126 * (x134 * x31 - x135 * x5))
    result[1, 5] = numpy.sum(x58 * (x136 * x31 - x137 * x5))
    result[2, 0] = numpy.sum(
        x58
        * (
            x31 * (x141 + x143 * x4 - x148 * x5)
            + x33 * (x143 - x144 + x147)
            - x5 * (x148 * x4 + x149 - x5 * (x146 * x4 - x150 * x5))
        )
    )
    result[2, 1] = numpy.sum(
        x126 * (x155 + x31 * (x154 * x4 - x157 * x5) - x5 * (x157 * x4 - x159 * x5))
    )
    result[2, 2] = numpy.sum(
        x126 * (x162 + x31 * (x161 * x4 - x163 * x5) - x5 * (x163 * x4 - x164 * x5))
    )
    result[2, 3] = numpy.sum(x58 * (x166 * x31 - x168 * x5))
    result[2, 4] = numpy.sum(x126 * (x170 * x31 - x172 * x5))
    result[2, 5] = numpy.sum(x58 * (x173 * x31 - x174 * x5))
    result[3, 0] = numpy.sum(
        x37
        * (
            x0 * (x175 - x176 + x177 + x29)
            + x4 * (x177 * x4 - x178 * x5)
            - x5 * (x178 * x4 - x5 * (x108 * x116 - x120 * x7 + x34))
        )
    )
    result[3, 1] = numpy.sum(
        x58
        * (
            x4 * (x108 * x122 + x111 - x124 * x7 + x53)
            - x5 * (x108 * x124 + x119 - x125 * x7 + x55)
        )
    )
    result[3, 2] = numpy.sum(
        x58
        * (x4 * (x108 * x128 - x130 * x7 + x74) - x5 * (x108 * x130 - x131 * x7 + x76))
    )
    result[3, 3] = numpy.sum(x37 * (x108 * x132 + 2.0 * x123 - x133 * x7 + x84))
    result[3, 4] = numpy.sum(x58 * (x108 * x134 + x129 - x135 * x7 + x93))
    result[3, 5] = numpy.sum(x37 * (x103 + x108 * x136 - x137 * x7))
    result[4, 0] = numpy.sum(
        x58
        * (
            x0 * (x151 - x179 + x180)
            + x4 * (x180 * x4 - x181 * x5)
            - x5 * (x181 * x4 - x5 * (x108 * x146 + x158))
        )
    )
    result[4, 1] = numpy.sum(
        x126 * (x4 * (x108 * x154 + x165) - x5 * (x108 * x157 + x167))
    )
    result[4, 2] = numpy.sum(
        x126 * (x4 * (x108 * x161 + x169) - x5 * (x108 * x163 + x171))
    )
    result[4, 3] = numpy.sum(x58 * (x108 * x166 + 2.0 * x155 - x168 * x7))
    result[4, 4] = numpy.sum(x126 * (x108 * x170 + x162 - x172 * x7))
    result[4, 5] = numpy.sum(x58 * (x108 * x173 - x174 * x7))
    result[5, 0] = numpy.sum(
        x37 * (x185 + x4 * (x184 * x4 - x186 * x5) - x5 * (x186 * x4 - x187 * x5))
    )
    result[5, 1] = numpy.sum(x58 * (x188 * x4 - x189 * x5))
    result[5, 2] = numpy.sum(x58 * (x190 * x4 - x191 * x5))
    result[5, 3] = numpy.sum(x37 * (x185 + x188 * x39 - x189 * x7))
    result[5, 4] = numpy.sum(x58 * (x190 * x39 - x191 * x7))
    result[5, 5] = numpy.sum(x37 * (x103 + x138 * x173 + 2.0 * x162 - x174 * x9))
    return result


def coulomb3d_23(ax, da, A, bx, db, B, C):
    """Cartesian (df) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 10), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = -x2 * (ax * A[0] + bx * B[0])
    x4 = -x3 - B[0]
    x5 = -x3 - C[0]
    x6 = -x2 * (ax * A[1] + bx * B[1])
    x7 = -x6 - C[1]
    x8 = -x2 * (ax * A[2] + bx * B[2])
    x9 = -x8 - C[2]
    x10 = x1 * (x5**2 + x7**2 + x9**2)
    x11 = (
        6.28318530717959
        * x2
        * numpy.exp(
            -ax * bx * x2 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x12 = x11 * boys(2, x10)
    x13 = x11 * boys(1, x10)
    x14 = x0 * (-x12 + x13)
    x15 = x12 * x5
    x16 = x13 * x4
    x17 = -x15 + x16
    x18 = x17 * x4
    x19 = x11 * boys(3, x10)
    x20 = x19 * x5
    x21 = x12 * x4
    x22 = -x20 + x21
    x23 = x22 * x5
    x24 = x14 + x18 - x23
    x25 = x24 * x4
    x26 = x11 * boys(0, x10)
    x27 = x0 * (-x13 + x26)
    x28 = -x13 * x5 + x26 * x4
    x29 = -x17 * x5 + x27 + x28 * x4
    x30 = 2.0 * x0
    x31 = -x24 * x5 + x30 * (x15 - x16 + x28)
    x32 = x29 * x4 + x31
    x33 = x0 * (x12 - x19)
    x34 = x22 * x4
    x35 = x11 * boys(4, x10)
    x36 = x35 * x5
    x37 = x19 * x4
    x38 = -x36 + x37
    x39 = x38 * x5
    x40 = x33 + x34 - x39
    x41 = x40 * x5
    x42 = x30 * (x17 + x20 - x21)
    x43 = x41 - x42
    x44 = -x3 - A[0]
    x45 = x25 - x41 + x42
    x46 = -x14
    x47 = 3.0 * x0
    x48 = x0 * (x19 - x35)
    x49 = x11 * boys(5, x10)
    x50 = -x33
    x51 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**4.5)
    x52 = 13.4923846833851 * x51
    x53 = x12 * x7
    x54 = -x53
    x55 = -x6 - B[1]
    x56 = x13 * x55
    x57 = x54 + x56
    x58 = x4 * x57
    x59 = x19 * x7
    x60 = -x59
    x61 = x12 * x55
    x62 = x60 + x61
    x63 = x5 * x62
    x64 = x58 - x63
    x65 = x4 * x64
    x66 = -x13 * x7
    x67 = x26 * x55 + x66
    x68 = x4 * x67 - x5 * x57
    x69 = x0 * (x53 - x56 + x67)
    x70 = -x5 * x64 + x69
    x71 = x4 * x68 + x70
    x72 = x4 * x62
    x73 = x35 * x7
    x74 = -x73
    x75 = x19 * x55
    x76 = x74 + x75
    x77 = x5 * x76
    x78 = x72 - x77
    x79 = x5 * x78
    x80 = x0 * (x57 + x59 - x61)
    x81 = -x80
    x82 = x79 + x81
    x83 = x65 - x79 + x80
    x84 = x0 * (x62 + x73 - x75)
    x85 = -x49 * x7
    x86 = x35 * x55 + x85
    x87 = 30.169889330626 * x51
    x88 = x12 * x9
    x89 = -x88
    x90 = -x8 - B[2]
    x91 = x13 * x90
    x92 = x89 + x91
    x93 = x4 * x92
    x94 = x19 * x9
    x95 = -x94
    x96 = x12 * x90
    x97 = x95 + x96
    x98 = x5 * x97
    x99 = x93 - x98
    x100 = x4 * x99
    x101 = -x13 * x9
    x102 = x101 + x26 * x90
    x103 = x102 * x4 - x5 * x92
    x104 = x0 * (x102 + x88 - x91)
    x105 = x104 - x5 * x99
    x106 = x103 * x4 + x105
    x107 = x4 * x97
    x108 = x35 * x9
    x109 = -x108
    x110 = x19 * x90
    x111 = x109 + x110
    x112 = x111 * x5
    x113 = x107 - x112
    x114 = x113 * x5
    x115 = x0 * (x92 + x94 - x96)
    x116 = -x115
    x117 = x114 + x116
    x118 = x100 - x114 + x115
    x119 = x0 * (x108 - x110 + x97)
    x120 = -x49 * x9
    x121 = x120 + x35 * x90
    x122 = x27 - x57 * x7
    x123 = x122 + x55 * x67
    x124 = x55 * x62
    x125 = x7 * x76
    x126 = -x125 + x33
    x127 = x124 + x126
    x128 = x127 * x5
    x129 = x55 * x57
    x130 = x62 * x7
    x131 = -x130 + x14
    x132 = x129 + x131
    x133 = -x132 * x5
    x134 = x132 * x4
    x135 = x123 * x4 + x133
    x136 = x130 + x46
    x137 = x0 * (x123 - x129 + x136)
    x138 = -x128 + x134
    x139 = x125 + x50
    x140 = x0 * (-x124 + x132 + x139)
    x141 = x48 - x7 * x86
    x142 = x141 + x55 * x76
    x143 = -x7 * x92
    x144 = x102 * x55 + x143
    x145 = x55 * x97
    x146 = x111 * x7
    x147 = -x146
    x148 = x145 + x147
    x149 = x148 * x5
    x150 = x55 * x92
    x151 = x7 * x97
    x152 = -x151
    x153 = x150 + x152
    x154 = -x153 * x5
    x155 = x153 * x4
    x156 = x144 * x4 + x154
    x157 = x0 * (x144 - x150 + x151)
    x158 = -x149 + x155
    x159 = x0 * (-x145 + x146 + x153)
    x160 = -x121 * x7
    x161 = x111 * x55 + x160
    x162 = 52.2557811793745 * x51
    x163 = x27 - x9 * x92
    x164 = x102 * x90 + x163
    x165 = x90 * x97
    x166 = x111 * x9
    x167 = -x166 + x33
    x168 = x165 + x167
    x169 = x168 * x5
    x170 = x90 * x92
    x171 = x9 * x97
    x172 = x14 - x171
    x173 = x170 + x172
    x174 = -x173 * x5
    x175 = x173 * x4
    x176 = x164 * x4 + x174
    x177 = x171 + x46
    x178 = x0 * (x164 - x170 + x177)
    x179 = -x169 + x175
    x180 = x166 + x50
    x181 = x0 * (-x165 + x173 + x180)
    x182 = -x121 * x9 + x48
    x183 = x111 * x90 + x182
    x184 = x132 * x55
    x185 = -x132 * x7 + 2.0 * x69
    x186 = x123 * x55 + x185
    x187 = x127 * x7
    x188 = 2.0 * x80
    x189 = x187 - x188
    x190 = x0 * (-x184 + x186 + x189)
    x191 = -x187 + x188
    x192 = x184 + x191
    x193 = -x142 * x7 + 2.0 * x84
    x194 = x127 * x55 + x193
    x195 = x153 * x55
    x196 = x104 - x153 * x7
    x197 = x144 * x55 + x196
    x198 = x148 * x7
    x199 = x116 + x198
    x200 = x0 * (-x195 + x197 + x199)
    x201 = x115 - x198
    x202 = x195 + x201
    x203 = x119 - x161 * x7
    x204 = x148 * x55 + x203
    x205 = x168 * x7
    x206 = x173 * x55
    x207 = -x173 * x7
    x208 = x164 * x55 + x207
    x209 = x0 * (x205 - x206 + x208)
    x210 = -x205
    x211 = x206 + x210
    x212 = -x183 * x7
    x213 = x168 * x55 + x212
    x214 = x173 * x90
    x215 = 2.0 * x104 - x173 * x9
    x216 = x164 * x90 + x215
    x217 = x168 * x9
    x218 = 2.0 * x115
    x219 = x217 - x218
    x220 = x0 * (-x214 + x216 + x219)
    x221 = -x217 + x218
    x222 = x214 + x221
    x223 = 2.0 * x119 - x183 * x9
    x224 = x168 * x90 + x223
    x225 = -x6 - A[1]
    x226 = x13 * x225
    x227 = x225 * x26 + x66
    x228 = x0 * (-x226 + x227 + x53)
    x229 = x226 + x54
    x230 = x227 * x4 - x229 * x5
    x231 = x229 * x4
    x232 = x12 * x225
    x233 = x232 + x60
    x234 = x233 * x5
    x235 = x231 - x234
    x236 = x228 + x230 * x4 - x235 * x5
    x237 = x0 * (x229 - x232 + x59)
    x238 = x235 * x4
    x239 = x233 * x4
    x240 = x19 * x225
    x241 = x240 + x74
    x242 = x241 * x5
    x243 = x239 - x242
    x244 = x243 * x5
    x245 = x237 + x238 - x244
    x246 = x0 * (x233 - x240 + x73)
    x247 = x225 * x35 + x85
    x248 = -x237
    x249 = 23.3694957868871 * x51
    x250 = x225 * x57
    x251 = x122 + x225 * x67
    x252 = x0 * (x136 - x250 + x251)
    x253 = x131 + x250
    x254 = x251 * x4 - x253 * x5
    x255 = x253 * x4
    x256 = x225 * x62
    x257 = x126 + x256
    x258 = x257 * x5
    x259 = x255 - x258
    x260 = x0 * (x139 + x253 - x256)
    x261 = x141 + x225 * x76
    x262 = x225 * x92
    x263 = x102 * x225 + x143
    x264 = x0 * (x151 - x262 + x263)
    x265 = x152 + x262
    x266 = x263 * x4 - x265 * x5
    x267 = x265 * x4
    x268 = x225 * x97
    x269 = x147 + x268
    x270 = x269 * x5
    x271 = x267 - x270
    x272 = x0 * (x146 + x265 - x268)
    x273 = x111 * x225 + x160
    x274 = x132 * x225
    x275 = x123 * x225 + x185
    x276 = x0 * (x189 - x274 + x275)
    x277 = x191 + x274
    x278 = x127 * x225 + x193
    x279 = x153 * x225
    x280 = x144 * x225 + x196
    x281 = x0 * (x199 - x279 + x280)
    x282 = x201 + x279
    x283 = x148 * x225 + x203
    x284 = 90.5096679918781 * x51
    x285 = x173 * x225
    x286 = x164 * x225 + x207
    x287 = x0 * (x205 - x285 + x286)
    x288 = x210 + x285
    x289 = x168 * x225 + x212
    x290 = 3.0 * x137 + x186 * x225 - x192 * x7
    x291 = 3.0 * x140 + x192 * x225 - x194 * x7
    x292 = 2.0 * x157 + x197 * x225 - x202 * x7
    x293 = 2.0 * x159 + x202 * x225 - x204 * x7
    x294 = x178 + x208 * x225 - x211 * x7
    x295 = x181 + x211 * x225 - x213 * x7
    x296 = x216 * x225 - x222 * x7
    x297 = x222 * x225 - x224 * x7
    x298 = -x8 - A[2]
    x299 = x13 * x298
    x300 = x101 + x26 * x298
    x301 = x0 * (-x299 + x300 + x88)
    x302 = x299 + x89
    x303 = x300 * x4 - x302 * x5
    x304 = x302 * x4
    x305 = x12 * x298
    x306 = x305 + x95
    x307 = x306 * x5
    x308 = x304 - x307
    x309 = x301 + x303 * x4 - x308 * x5
    x310 = x0 * (x302 - x305 + x94)
    x311 = x308 * x4
    x312 = x306 * x4
    x313 = x19 * x298
    x314 = x109 + x313
    x315 = x314 * x5
    x316 = x312 - x315
    x317 = x316 * x5
    x318 = x310 + x311 - x317
    x319 = x0 * (x108 + x306 - x313)
    x320 = x120 + x298 * x35
    x321 = -x310
    x322 = x306 * x7
    x323 = x302 * x55
    x324 = -x302 * x7
    x325 = x300 * x55 + x324
    x326 = x0 * (x322 - x323 + x325)
    x327 = -x322
    x328 = x323 + x327
    x329 = x325 * x4 - x328 * x5
    x330 = x328 * x4
    x331 = x306 * x55
    x332 = x314 * x7
    x333 = -x332
    x334 = x331 + x333
    x335 = x334 * x5
    x336 = x330 - x335
    x337 = x0 * (x328 - x331 + x332)
    x338 = -x320 * x7
    x339 = x314 * x55 + x338
    x340 = x298 * x92
    x341 = x102 * x298 + x163
    x342 = x0 * (x177 - x340 + x341)
    x343 = x172 + x340
    x344 = x341 * x4 - x343 * x5
    x345 = x343 * x4
    x346 = x298 * x97
    x347 = x167 + x346
    x348 = x347 * x5
    x349 = x345 - x348
    x350 = x0 * (x180 + x343 - x346)
    x351 = x111 * x298 + x182
    x352 = x328 * x55
    x353 = x301 - x328 * x7
    x354 = x325 * x55 + x353
    x355 = x334 * x7
    x356 = x321 + x355
    x357 = x0 * (-x352 + x354 + x356)
    x358 = x310 - x355
    x359 = x352 + x358
    x360 = x319 - x339 * x7
    x361 = x334 * x55 + x360
    x362 = x347 * x7
    x363 = x343 * x55
    x364 = -x343 * x7
    x365 = x341 * x55 + x364
    x366 = x0 * (x362 - x363 + x365)
    x367 = -x362
    x368 = x363 + x367
    x369 = -x351 * x7
    x370 = x347 * x55 + x369
    x371 = x173 * x298
    x372 = x164 * x298 + x215
    x373 = x0 * (x219 - x371 + x372)
    x374 = x221 + x371
    x375 = x168 * x298 + x223
    x376 = 2.0 * x326 - x359 * x7
    x377 = x354 * x55 + x376
    x378 = 2.0 * x337 - x361 * x7
    x379 = x359 * x55 + x378
    x380 = x342 - x368 * x7
    x381 = x365 * x55 + x380
    x382 = x350 - x370 * x7
    x383 = x368 * x55 + x382
    x384 = -x374 * x7
    x385 = x372 * x55 + x384
    x386 = -x375 * x7
    x387 = x374 * x55 + x386
    x388 = 3.0 * x178 + x216 * x298 - x222 * x9
    x389 = 3.0 * x181 + x222 * x298 - x224 * x9
    x390 = x233 * x7
    x391 = x225 * x229
    x392 = x225 * x227 - x229 * x7 + x27
    x393 = x14 - x390 + x391
    x394 = x392 * x4 - x393 * x5
    x395 = x393 * x4
    x396 = x225 * x233
    x397 = x241 * x7
    x398 = x33 + x396 - x397
    x399 = x398 * x5
    x400 = x395 - x399
    x401 = x257 * x7
    x402 = x225 * x253
    x403 = x225 * x251 + x228 - x253 * x7 + x69
    x404 = x237 - x401 + x402 + x80
    x405 = x269 * x7
    x406 = x225 * x265
    x407 = x104 + x225 * x263 - x265 * x7
    x408 = x115 - x405 + x406
    x409 = x225 * x302
    x410 = x225 * x300 + x324
    x411 = x327 + x409
    x412 = x4 * x410 - x411 * x5
    x413 = x4 * x411
    x414 = x225 * x306
    x415 = x333 + x414
    x416 = x415 * x5
    x417 = x413 - x416
    x418 = x225 * x328
    x419 = x225 * x325 + x353
    x420 = x358 + x418
    x421 = x225 * x343
    x422 = x225 * x341 + x364
    x423 = x367 + x421
    x424 = x306 * x9
    x425 = x298 * x302
    x426 = x27 + x298 * x300 - x302 * x9
    x427 = x0 * (x424 - x425 + x426 + x46)
    x428 = x14 - x424 + x425
    x429 = x4 * x426 - x428 * x5
    x430 = x4 * x428
    x431 = x298 * x306
    x432 = x314 * x9
    x433 = x33 + x431 - x432
    x434 = x433 * x5
    x435 = x430 - x434
    x436 = x0 * (x428 - x431 + x432 + x50)
    x437 = x298 * x314 - x320 * x9 + x48
    x438 = x433 * x7
    x439 = x428 * x55
    x440 = x426 * x55 - x428 * x7
    x441 = x0 * (x438 - x439 + x440)
    x442 = -x438 + x439
    x443 = x433 * x55 - x437 * x7
    x444 = x347 * x9
    x445 = x298 * x343
    x446 = x104 + x298 * x341 + x301 - x343 * x9
    x447 = x0 * (x116 + x321 + x444 - x445 + x446)
    x448 = x115 + x310 - x444 + x445
    x449 = x119 + x298 * x347 + x319 - x351 * x9
    x450 = x427 + x440 * x55 - x442 * x7
    x451 = x436 + x442 * x55 - x443 * x7
    x452 = x446 * x55 - x448 * x7
    x453 = x448 * x55 - x449 * x7
    x454 = x178 + x298 * x372 + 2.0 * x342 - x374 * x9
    x455 = x181 + x298 * x374 + 2.0 * x350 - x375 * x9

    # 60 item(s)
    result[0, 0] = numpy.sum(
        x52
        * (
            x0 * (-x25 + x32 + x43)
            + x44 * (x32 * x44 - x45 * x5 + x47 * (-x18 + x23 + x29 + x46))
            + x47 * (-x24 * x44 + x29 * x44 + x31 + x43)
            - x5
            * (
                x44 * x45
                + x47 * (x24 - x34 + x39 + x50)
                - x5
                * (
                    x30 * (x22 + x36 - x37)
                    + x4 * x40
                    - x5 * (x38 * x4 + x48 - x5 * (x35 * x4 - x49 * x5))
                )
            )
        )
    )
    result[0, 1] = numpy.sum(
        x87
        * (
            x0 * (-x65 + x71 + x82)
            + x30 * (-x44 * x64 + x44 * x68 + x70 + x82)
            + x44 * (x30 * (-x58 + x63 + x68) + x44 * x71 - x5 * x83)
            - x5
            * (
                x30 * (x64 - x72 + x77)
                + x44 * x83
                - x5 * (x4 * x78 - x5 * (x4 * x76 - x5 * x86) + x84)
            )
        )
    )
    result[0, 2] = numpy.sum(
        x87
        * (
            x0 * (-x100 + x106 + x117)
            + x30 * (x103 * x44 + x105 + x117 - x44 * x99)
            + x44 * (x106 * x44 - x118 * x5 + x30 * (x103 - x93 + x98))
            - x5
            * (
                x118 * x44
                + x30 * (-x107 + x112 + x99)
                - x5 * (x113 * x4 + x119 - x5 * (x111 * x4 - x121 * x5))
            )
        )
    )
    result[0, 3] = numpy.sum(
        x87
        * (
            x0 * (x128 - x134 + x135)
            + x0 * (x123 * x44 + x128 - x132 * x44 + x133)
            + x44 * (x135 * x44 + x137 - x138 * x5)
            - x5 * (x138 * x44 + x140 - x5 * (x127 * x4 - x142 * x5))
        )
    )
    result[0, 4] = numpy.sum(
        x162
        * (
            x0 * (x149 - x155 + x156)
            + x0 * (x144 * x44 + x149 - x153 * x44 + x154)
            + x44 * (x156 * x44 + x157 - x158 * x5)
            - x5 * (x158 * x44 + x159 - x5 * (x148 * x4 - x161 * x5))
        )
    )
    result[0, 5] = numpy.sum(
        x87
        * (
            x0 * (x169 - x175 + x176)
            + x0 * (x164 * x44 + x169 - x173 * x44 + x174)
            + x44 * (x176 * x44 + x178 - x179 * x5)
            - x5 * (x179 * x44 + x181 - x5 * (x168 * x4 - x183 * x5))
        )
    )
    result[0, 6] = numpy.sum(
        x52 * (x190 + x44 * (x186 * x44 - x192 * x5) - x5 * (x192 * x44 - x194 * x5))
    )
    result[0, 7] = numpy.sum(
        x87 * (x200 + x44 * (x197 * x44 - x202 * x5) - x5 * (x202 * x44 - x204 * x5))
    )
    result[0, 8] = numpy.sum(
        x87 * (x209 + x44 * (x208 * x44 - x211 * x5) - x5 * (x211 * x44 - x213 * x5))
    )
    result[0, 9] = numpy.sum(
        x52 * (x220 + x44 * (x216 * x44 - x222 * x5) - x5 * (x222 * x44 - x224 * x5))
    )
    result[1, 0] = numpy.sum(
        x249
        * (
            x44 * (x236 * x4 - x245 * x5 + x30 * (x230 - x231 + x234))
            + x47 * (x236 - x238 + x244 + x248)
            - x5
            * (
                x245 * x4
                + x30 * (x235 - x239 + x242)
                - x5 * (x243 * x4 + x246 - x5 * (x241 * x4 - x247 * x5))
            )
        )
    )
    result[1, 1] = numpy.sum(
        x162
        * (
            x30 * (x254 - x255 + x258)
            + x44 * (x252 + x254 * x4 - x259 * x5)
            - x5 * (x259 * x4 + x260 - x5 * (x257 * x4 - x261 * x5))
        )
    )
    result[1, 2] = numpy.sum(
        x162
        * (
            x30 * (x266 - x267 + x270)
            + x44 * (x264 + x266 * x4 - x271 * x5)
            - x5 * (x271 * x4 + x272 - x5 * (x269 * x4 - x273 * x5))
        )
    )
    result[1, 3] = numpy.sum(
        x162 * (x276 + x44 * (x275 * x4 - x277 * x5) - x5 * (x277 * x4 - x278 * x5))
    )
    result[1, 4] = numpy.sum(
        x284 * (x281 + x44 * (x280 * x4 - x282 * x5) - x5 * (x282 * x4 - x283 * x5))
    )
    result[1, 5] = numpy.sum(
        x162 * (x287 + x44 * (x286 * x4 - x288 * x5) - x5 * (x288 * x4 - x289 * x5))
    )
    result[1, 6] = numpy.sum(x249 * (x290 * x44 - x291 * x5))
    result[1, 7] = numpy.sum(x162 * (x292 * x44 - x293 * x5))
    result[1, 8] = numpy.sum(x162 * (x294 * x44 - x295 * x5))
    result[1, 9] = numpy.sum(x249 * (x296 * x44 - x297 * x5))
    result[2, 0] = numpy.sum(
        x249
        * (
            x44 * (x30 * (x303 - x304 + x307) + x309 * x4 - x318 * x5)
            + x47 * (x309 - x311 + x317 + x321)
            - x5
            * (
                x30 * (x308 - x312 + x315)
                + x318 * x4
                - x5 * (x316 * x4 + x319 - x5 * (x314 * x4 - x320 * x5))
            )
        )
    )
    result[2, 1] = numpy.sum(
        x162
        * (
            x30 * (x329 - x330 + x335)
            + x44 * (x326 + x329 * x4 - x336 * x5)
            - x5 * (x336 * x4 + x337 - x5 * (x334 * x4 - x339 * x5))
        )
    )
    result[2, 2] = numpy.sum(
        x162
        * (
            x30 * (x344 - x345 + x348)
            + x44 * (x342 + x344 * x4 - x349 * x5)
            - x5 * (x349 * x4 + x350 - x5 * (x347 * x4 - x351 * x5))
        )
    )
    result[2, 3] = numpy.sum(
        x162 * (x357 + x44 * (x354 * x4 - x359 * x5) - x5 * (x359 * x4 - x361 * x5))
    )
    result[2, 4] = numpy.sum(
        x284 * (x366 + x44 * (x365 * x4 - x368 * x5) - x5 * (x368 * x4 - x370 * x5))
    )
    result[2, 5] = numpy.sum(
        x162 * (x373 + x44 * (x372 * x4 - x374 * x5) - x5 * (x374 * x4 - x375 * x5))
    )
    result[2, 6] = numpy.sum(x249 * (x377 * x44 - x379 * x5))
    result[2, 7] = numpy.sum(x162 * (x381 * x44 - x383 * x5))
    result[2, 8] = numpy.sum(x162 * (x385 * x44 - x387 * x5))
    result[2, 9] = numpy.sum(x249 * (x388 * x44 - x389 * x5))
    result[3, 0] = numpy.sum(
        x52
        * (
            x30 * (x394 - x395 + x399)
            + x4 * (x0 * (x390 - x391 + x392 + x46) + x394 * x4 - x400 * x5)
            - x5
            * (
                x0 * (x393 - x396 + x397 + x50)
                + x4 * x400
                - x5 * (x398 * x4 - x5 * (x225 * x241 - x247 * x7 + x48))
            )
        )
    )
    result[3, 1] = numpy.sum(
        x87
        * (
            x0 * (x248 + x401 - x402 + x403 + x81)
            + x4 * (x4 * x403 - x404 * x5)
            - x5 * (x4 * x404 - x5 * (x225 * x257 + x246 - x261 * x7 + x84))
        )
    )
    result[3, 2] = numpy.sum(
        x87
        * (
            x0 * (x116 + x405 - x406 + x407)
            + x4 * (x4 * x407 - x408 * x5)
            - x5 * (x4 * x408 - x5 * (x119 + x225 * x269 - x273 * x7))
        )
    )
    result[3, 3] = numpy.sum(
        x87
        * (
            x4 * (x137 + x225 * x275 + 2.0 * x252 - x277 * x7)
            - x5 * (x140 + x225 * x277 + 2.0 * x260 - x278 * x7)
        )
    )
    result[3, 4] = numpy.sum(
        x162
        * (
            x4 * (x157 + x225 * x280 + x264 - x282 * x7)
            - x5 * (x159 + x225 * x282 + x272 - x283 * x7)
        )
    )
    result[3, 5] = numpy.sum(
        x87
        * (x4 * (x178 + x225 * x286 - x288 * x7) - x5 * (x181 + x225 * x288 - x289 * x7))
    )
    result[3, 6] = numpy.sum(x52 * (x190 + x225 * x290 + 3.0 * x276 - x291 * x7))
    result[3, 7] = numpy.sum(x87 * (x200 + x225 * x292 + 2.0 * x281 - x293 * x7))
    result[3, 8] = numpy.sum(x87 * (x209 + x225 * x294 + x287 - x295 * x7))
    result[3, 9] = numpy.sum(x52 * (x220 + x225 * x296 - x297 * x7))
    result[4, 0] = numpy.sum(
        x249
        * (
            x30 * (x412 - x413 + x416)
            + x4 * (x0 * (x322 - x409 + x410) + x4 * x412 - x417 * x5)
            - x5
            * (
                x0 * (x332 + x411 - x414)
                + x4 * x417
                - x5 * (x4 * x415 - x5 * (x225 * x314 + x338))
            )
        )
    )
    result[4, 1] = numpy.sum(
        x162
        * (
            x0 * (x356 - x418 + x419)
            + x4 * (x4 * x419 - x420 * x5)
            - x5 * (x4 * x420 - x5 * (x225 * x334 + x360))
        )
    )
    result[4, 2] = numpy.sum(
        x162
        * (
            x0 * (x362 - x421 + x422)
            + x4 * (x4 * x422 - x423 * x5)
            - x5 * (x4 * x423 - x5 * (x225 * x347 + x369))
        )
    )
    result[4, 3] = numpy.sum(
        x162 * (x4 * (x225 * x354 + x376) - x5 * (x225 * x359 + x378))
    )
    result[4, 4] = numpy.sum(
        x284 * (x4 * (x225 * x365 + x380) - x5 * (x225 * x368 + x382))
    )
    result[4, 5] = numpy.sum(
        x162 * (x4 * (x225 * x372 + x384) - x5 * (x225 * x374 + x386))
    )
    result[4, 6] = numpy.sum(x249 * (x225 * x377 + 3.0 * x357 - x379 * x7))
    result[4, 7] = numpy.sum(x162 * (x225 * x381 + 2.0 * x366 - x383 * x7))
    result[4, 8] = numpy.sum(x162 * (x225 * x385 + x373 - x387 * x7))
    result[4, 9] = numpy.sum(x249 * (x225 * x388 - x389 * x7))
    result[5, 0] = numpy.sum(
        x52
        * (
            x30 * (x429 - x430 + x434)
            + x4 * (x4 * x429 + x427 - x435 * x5)
            - x5 * (x4 * x435 + x436 - x5 * (x4 * x433 - x437 * x5))
        )
    )
    result[5, 1] = numpy.sum(
        x87 * (x4 * (x4 * x440 - x442 * x5) + x441 - x5 * (x4 * x442 - x443 * x5))
    )
    result[5, 2] = numpy.sum(
        x87 * (x4 * (x4 * x446 - x448 * x5) + x447 - x5 * (x4 * x448 - x449 * x5))
    )
    result[5, 3] = numpy.sum(x87 * (x4 * x450 - x451 * x5))
    result[5, 4] = numpy.sum(x162 * (x4 * x452 - x453 * x5))
    result[5, 5] = numpy.sum(x87 * (x4 * x454 - x455 * x5))
    result[5, 6] = numpy.sum(x52 * (2.0 * x441 + x450 * x55 - x451 * x7))
    result[5, 7] = numpy.sum(x87 * (x447 + x452 * x55 - x453 * x7))
    result[5, 8] = numpy.sum(x87 * (x454 * x55 - x455 * x7))
    result[5, 9] = numpy.sum(x52 * (x220 + x298 * x388 + 3.0 * x373 - x389 * x9))
    return result


def coulomb3d_24(ax, da, A, bx, db, B, C):
    """Cartesian (dg) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 15), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = -x2 * (ax * A[0] + bx * B[0])
    x4 = -x3 - B[0]
    x5 = -x3 - C[0]
    x6 = -x2 * (ax * A[1] + bx * B[1])
    x7 = -x6 - C[1]
    x8 = -x2 * (ax * A[2] + bx * B[2])
    x9 = -x8 - C[2]
    x10 = x1 * (x5**2 + x7**2 + x9**2)
    x11 = (
        6.28318530717959
        * x2
        * numpy.exp(
            -ax * bx * x2 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x12 = x11 * boys(2, x10)
    x13 = x11 * boys(1, x10)
    x14 = x0 * (-x12 + x13)
    x15 = x12 * x5
    x16 = x13 * x4
    x17 = -x15 + x16
    x18 = x17 * x4
    x19 = x11 * boys(3, x10)
    x20 = x19 * x5
    x21 = x12 * x4
    x22 = -x20 + x21
    x23 = x22 * x5
    x24 = x14 + x18 - x23
    x25 = x24 * x4
    x26 = x0 * (x12 - x19)
    x27 = x22 * x4
    x28 = x11 * boys(4, x10)
    x29 = x28 * x5
    x30 = x19 * x4
    x31 = -x29 + x30
    x32 = x31 * x5
    x33 = x26 + x27 - x32
    x34 = x33 * x5
    x35 = 2.0 * x0
    x36 = x35 * (x17 + x20 - x21)
    x37 = x25 - x34 + x36
    x38 = x37 * x4
    x39 = x11 * boys(0, x10)
    x40 = x0 * (-x13 + x39)
    x41 = -x13 * x5 + x39 * x4
    x42 = -x17 * x5 + x4 * x41 + x40
    x43 = -x24 * x5 + x35 * (x15 - x16 + x41) + x4 * x42
    x44 = -x14
    x45 = 3.0 * x0
    x46 = -x37 * x5 + x45 * (-x18 + x23 + x42 + x44)
    x47 = x4 * x43 + x46
    x48 = x33 * x4
    x49 = x0 * (x19 - x28)
    x50 = x31 * x4
    x51 = x11 * boys(5, x10)
    x52 = x5 * x51
    x53 = x28 * x4
    x54 = -x52 + x53
    x55 = x5 * x54
    x56 = x49 + x50 - x55
    x57 = x5 * x56
    x58 = x35 * (x22 + x29 - x30)
    x59 = x48 - x57 + x58
    x60 = x5 * x59
    x61 = -x26
    x62 = x45 * (x24 - x27 + x32 + x61)
    x63 = x60 - x62
    x64 = -x3 - A[0]
    x65 = x38 - x60 + x62
    x66 = 4.0 * x0
    x67 = x0 * (x28 - x51)
    x68 = x11 * boys(6, x10)
    x69 = -x49
    x70 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**5.5)
    x71 = 10.1992841329868 * x70
    x72 = -x6 - B[1]
    x73 = x12 * x72
    x74 = x19 * x7
    x75 = x12 * x7
    x76 = -x75
    x77 = x13 * x72
    x78 = x76 + x77
    x79 = x0 * (-x73 + x74 + x78)
    x80 = x4 * x78
    x81 = -x74
    x82 = x73 + x81
    x83 = x5 * x82
    x84 = x80 - x83
    x85 = x4 * x84
    x86 = x4 * x82
    x87 = x28 * x7
    x88 = -x87
    x89 = x19 * x72
    x90 = x88 + x89
    x91 = x5 * x90
    x92 = x86 - x91
    x93 = x5 * x92
    x94 = x79 + x85 - x93
    x95 = x4 * x94
    x96 = -x13 * x7
    x97 = x39 * x72 + x96
    x98 = x0 * (x75 - x77 + x97)
    x99 = x4 * x97 - x5 * x78
    x100 = x4 * x99 - x5 * x84 + x98
    x101 = x35 * (-x80 + x83 + x99) - x5 * x94
    x102 = x100 * x4 + x101
    x103 = x0 * (x82 + x87 - x89)
    x104 = x4 * x92
    x105 = x4 * x90
    x106 = x51 * x7
    x107 = -x106
    x108 = x28 * x72
    x109 = x107 + x108
    x110 = x109 * x5
    x111 = x105 - x110
    x112 = x111 * x5
    x113 = x103 + x104 - x112
    x114 = x113 * x5
    x115 = x35 * (x84 - x86 + x91)
    x116 = x114 - x115
    x117 = -x114 + x115 + x95
    x118 = -x79
    x119 = x0 * (x106 - x108 + x90)
    x120 = -x68 * x7
    x121 = x120 + x51 * x72
    x122 = -x103
    x123 = 26.9847693667702 * x70
    x124 = -x8 - B[2]
    x125 = x12 * x124
    x126 = x19 * x9
    x127 = x12 * x9
    x128 = -x127
    x129 = x124 * x13
    x130 = x128 + x129
    x131 = x0 * (-x125 + x126 + x130)
    x132 = x130 * x4
    x133 = -x126
    x134 = x125 + x133
    x135 = x134 * x5
    x136 = x132 - x135
    x137 = x136 * x4
    x138 = x134 * x4
    x139 = x28 * x9
    x140 = -x139
    x141 = x124 * x19
    x142 = x140 + x141
    x143 = x142 * x5
    x144 = x138 - x143
    x145 = x144 * x5
    x146 = x131 + x137 - x145
    x147 = x146 * x4
    x148 = -x13 * x9
    x149 = x124 * x39 + x148
    x150 = x0 * (x127 - x129 + x149)
    x151 = -x130 * x5 + x149 * x4
    x152 = -x136 * x5 + x150 + x151 * x4
    x153 = -x146 * x5 + x35 * (-x132 + x135 + x151)
    x154 = x152 * x4 + x153
    x155 = x0 * (x134 + x139 - x141)
    x156 = x144 * x4
    x157 = x142 * x4
    x158 = x51 * x9
    x159 = -x158
    x160 = x124 * x28
    x161 = x159 + x160
    x162 = x161 * x5
    x163 = x157 - x162
    x164 = x163 * x5
    x165 = x155 + x156 - x164
    x166 = x165 * x5
    x167 = x35 * (x136 - x138 + x143)
    x168 = x166 - x167
    x169 = x147 - x166 + x167
    x170 = -x131
    x171 = x0 * (x142 + x158 - x160)
    x172 = -x68 * x9
    x173 = x124 * x51 + x172
    x174 = -x155
    x175 = x72 * x78
    x176 = x7 * x82
    x177 = x14 - x176
    x178 = x175 + x177
    x179 = x178 * x4
    x180 = x72 * x82
    x181 = x7 * x90
    x182 = -x181 + x26
    x183 = x180 + x182
    x184 = x183 * x5
    x185 = x179 - x184
    x186 = x185 * x4
    x187 = x40 - x7 * x78
    x188 = x187 + x72 * x97
    x189 = -x178 * x5 + x188 * x4
    x190 = x176 + x44
    x191 = x0 * (-x175 + x188 + x190)
    x192 = -x185 * x5 + x191
    x193 = x189 * x4 + x192
    x194 = x183 * x4
    x195 = x72 * x90
    x196 = x109 * x7
    x197 = -x196 + x49
    x198 = x195 + x197
    x199 = x198 * x5
    x200 = x194 - x199
    x201 = x200 * x5
    x202 = x181 + x61
    x203 = x0 * (x178 - x180 + x202)
    x204 = -x203
    x205 = x201 + x204
    x206 = x186 - x201 + x203
    x207 = x196 + x69
    x208 = x0 * (x183 - x195 + x207)
    x209 = -x121 * x7 + x67
    x210 = x109 * x72 + x209
    x211 = 34.8371874529163 * x70
    x212 = x130 * x72
    x213 = x134 * x7
    x214 = -x213
    x215 = x212 + x214
    x216 = x215 * x4
    x217 = x134 * x72
    x218 = x142 * x7
    x219 = -x218
    x220 = x217 + x219
    x221 = x220 * x5
    x222 = x216 - x221
    x223 = x222 * x4
    x224 = -x130 * x7
    x225 = x149 * x72 + x224
    x226 = -x215 * x5 + x225 * x4
    x227 = x0 * (-x212 + x213 + x225)
    x228 = -x222 * x5 + x227
    x229 = x226 * x4 + x228
    x230 = x220 * x4
    x231 = x142 * x72
    x232 = x161 * x7
    x233 = -x232
    x234 = x231 + x233
    x235 = x234 * x5
    x236 = x230 - x235
    x237 = x236 * x5
    x238 = x0 * (x215 - x217 + x218)
    x239 = -x238
    x240 = x237 + x239
    x241 = x223 - x237 + x238
    x242 = x0 * (x220 - x231 + x232)
    x243 = -x173 * x7
    x244 = x161 * x72 + x243
    x245 = 60.3397786612521 * x70
    x246 = x124 * x130
    x247 = x134 * x9
    x248 = x14 - x247
    x249 = x246 + x248
    x250 = x249 * x4
    x251 = x124 * x134
    x252 = x142 * x9
    x253 = -x252 + x26
    x254 = x251 + x253
    x255 = x254 * x5
    x256 = x250 - x255
    x257 = x256 * x4
    x258 = -x130 * x9 + x40
    x259 = x124 * x149 + x258
    x260 = -x249 * x5 + x259 * x4
    x261 = x247 + x44
    x262 = x0 * (-x246 + x259 + x261)
    x263 = -x256 * x5 + x262
    x264 = x260 * x4 + x263
    x265 = x254 * x4
    x266 = x124 * x142
    x267 = x161 * x9
    x268 = -x267 + x49
    x269 = x266 + x268
    x270 = x269 * x5
    x271 = x265 - x270
    x272 = x271 * x5
    x273 = x252 + x61
    x274 = x0 * (x249 - x251 + x273)
    x275 = -x274
    x276 = x272 + x275
    x277 = x257 - x272 + x274
    x278 = x267 + x69
    x279 = x0 * (x254 - x266 + x278)
    x280 = -x173 * x9 + x67
    x281 = x124 * x161 + x280
    x282 = -x178 * x7 + 2.0 * x98
    x283 = x188 * x72 + x282
    x284 = x183 * x72
    x285 = x198 * x7
    x286 = 2.0 * x103
    x287 = -x285 + x286
    x288 = x284 + x287
    x289 = x288 * x5
    x290 = x178 * x72
    x291 = x183 * x7
    x292 = 2.0 * x79
    x293 = -x291 + x292
    x294 = x290 + x293
    x295 = -x294 * x5
    x296 = x294 * x4
    x297 = x283 * x4 + x295
    x298 = x291 - x292
    x299 = x0 * (x283 - x290 + x298)
    x300 = -x289 + x296
    x301 = x285 - x286
    x302 = x0 * (-x284 + x294 + x301)
    x303 = 2.0 * x119 - x210 * x7
    x304 = x198 * x72 + x303
    x305 = x150 - x215 * x7
    x306 = x225 * x72 + x305
    x307 = x220 * x72
    x308 = x234 * x7
    x309 = x155 - x308
    x310 = x307 + x309
    x311 = x310 * x5
    x312 = x215 * x72
    x313 = x220 * x7
    x314 = x131 - x313
    x315 = x312 + x314
    x316 = -x315 * x5
    x317 = x315 * x4
    x318 = x306 * x4 + x316
    x319 = x170 + x313
    x320 = x0 * (x306 - x312 + x319)
    x321 = -x311 + x317
    x322 = x174 + x308
    x323 = x0 * (-x307 + x315 + x322)
    x324 = x171 - x244 * x7
    x325 = x234 * x72 + x324
    x326 = -x249 * x7
    x327 = x259 * x72 + x326
    x328 = x254 * x72
    x329 = x269 * x7
    x330 = -x329
    x331 = x328 + x330
    x332 = x331 * x5
    x333 = x249 * x72
    x334 = x254 * x7
    x335 = -x334
    x336 = x333 + x335
    x337 = -x336 * x5
    x338 = x336 * x4
    x339 = x327 * x4 + x337
    x340 = x0 * (x327 - x333 + x334)
    x341 = -x332 + x338
    x342 = x0 * (-x328 + x329 + x336)
    x343 = -x281 * x7
    x344 = x269 * x72 + x343
    x345 = 2.0 * x150 - x249 * x9
    x346 = x124 * x259 + x345
    x347 = x124 * x254
    x348 = x269 * x9
    x349 = 2.0 * x155
    x350 = -x348 + x349
    x351 = x347 + x350
    x352 = x351 * x5
    x353 = x124 * x249
    x354 = x254 * x9
    x355 = 2.0 * x131
    x356 = -x354 + x355
    x357 = x353 + x356
    x358 = -x357 * x5
    x359 = x357 * x4
    x360 = x346 * x4 + x358
    x361 = x354 - x355
    x362 = x0 * (x346 - x353 + x361)
    x363 = -x352 + x359
    x364 = x348 - x349
    x365 = x0 * (-x347 + x357 + x364)
    x366 = 2.0 * x171 - x281 * x9
    x367 = x124 * x269 + x366
    x368 = x294 * x72
    x369 = 3.0 * x191 - x294 * x7
    x370 = x283 * x72 + x369
    x371 = x288 * x7
    x372 = 3.0 * x203
    x373 = x371 - x372
    x374 = x0 * (-x368 + x370 + x373)
    x375 = -x371 + x372
    x376 = x368 + x375
    x377 = 3.0 * x208 - x304 * x7
    x378 = x288 * x72 + x377
    x379 = x315 * x72
    x380 = 2.0 * x227 - x315 * x7
    x381 = x306 * x72 + x380
    x382 = x310 * x7
    x383 = 2.0 * x238
    x384 = x382 - x383
    x385 = x0 * (-x379 + x381 + x384)
    x386 = -x382 + x383
    x387 = x379 + x386
    x388 = 2.0 * x242 - x325 * x7
    x389 = x310 * x72 + x388
    x390 = x336 * x72
    x391 = x262 - x336 * x7
    x392 = x327 * x72 + x391
    x393 = x331 * x7
    x394 = x275 + x393
    x395 = x0 * (-x390 + x392 + x394)
    x396 = x274 - x393
    x397 = x390 + x396
    x398 = x279 - x344 * x7
    x399 = x331 * x72 + x398
    x400 = x351 * x7
    x401 = x357 * x72
    x402 = -x357 * x7
    x403 = x346 * x72 + x402
    x404 = x0 * (x400 - x401 + x403)
    x405 = -x400
    x406 = x401 + x405
    x407 = -x367 * x7
    x408 = x351 * x72 + x407
    x409 = x124 * x357
    x410 = 3.0 * x262 - x357 * x9
    x411 = x124 * x346 + x410
    x412 = x351 * x9
    x413 = 3.0 * x274
    x414 = x412 - x413
    x415 = x0 * (-x409 + x411 + x414)
    x416 = -x412 + x413
    x417 = x409 + x416
    x418 = 3.0 * x279 - x367 * x9
    x419 = x124 * x351 + x418
    x420 = -x6 - A[1]
    x421 = x13 * x420
    x422 = x39 * x420 + x96
    x423 = x0 * (-x421 + x422 + x75)
    x424 = x421 + x76
    x425 = x4 * x422 - x424 * x5
    x426 = x4 * x424
    x427 = x12 * x420
    x428 = x427 + x81
    x429 = x428 * x5
    x430 = x426 - x429
    x431 = x4 * x425 + x423 - x430 * x5
    x432 = x0 * (x424 - x427 + x74)
    x433 = x4 * x430
    x434 = x4 * x428
    x435 = x19 * x420
    x436 = x435 + x88
    x437 = x436 * x5
    x438 = x434 - x437
    x439 = x438 * x5
    x440 = x432 + x433 - x439
    x441 = x35 * (x425 - x426 + x429) + x4 * x431 - x440 * x5
    x442 = x4 * x440
    x443 = x0 * (x428 - x435 + x87)
    x444 = x4 * x438
    x445 = x4 * x436
    x446 = x28 * x420
    x447 = x107 + x446
    x448 = x447 * x5
    x449 = x445 - x448
    x450 = x449 * x5
    x451 = x443 + x444 - x450
    x452 = x451 * x5
    x453 = x35 * (x430 - x434 + x437)
    x454 = x442 - x452 + x453
    x455 = -x432
    x456 = x0 * (x106 + x436 - x446)
    x457 = x120 + x420 * x51
    x458 = -x443
    x459 = 17.6656783191643 * x70
    x460 = x420 * x78
    x461 = x187 + x420 * x97
    x462 = x0 * (x190 - x460 + x461)
    x463 = x177 + x460
    x464 = x4 * x461 - x463 * x5
    x465 = x4 * x463
    x466 = x420 * x82
    x467 = x182 + x466
    x468 = x467 * x5
    x469 = x465 - x468
    x470 = x4 * x464 + x462 - x469 * x5
    x471 = x0 * (x202 + x463 - x466)
    x472 = x4 * x469
    x473 = x4 * x467
    x474 = x420 * x90
    x475 = x197 + x474
    x476 = x475 * x5
    x477 = x473 - x476
    x478 = x477 * x5
    x479 = x471 + x472 - x478
    x480 = x0 * (x207 + x467 - x474)
    x481 = x109 * x420 + x209
    x482 = 46.7389915737742 * x70
    x483 = x130 * x420
    x484 = x149 * x420 + x224
    x485 = x0 * (x213 - x483 + x484)
    x486 = x214 + x483
    x487 = x4 * x484 - x486 * x5
    x488 = x4 * x486
    x489 = x134 * x420
    x490 = x219 + x489
    x491 = x490 * x5
    x492 = x488 - x491
    x493 = x4 * x487 + x485 - x492 * x5
    x494 = x0 * (x218 + x486 - x489)
    x495 = x4 * x492
    x496 = x4 * x490
    x497 = x142 * x420
    x498 = x233 + x497
    x499 = x498 * x5
    x500 = x496 - x499
    x501 = x5 * x500
    x502 = x494 + x495 - x501
    x503 = x0 * (x232 + x490 - x497)
    x504 = x161 * x420 + x243
    x505 = -x494
    x506 = x178 * x420
    x507 = x188 * x420 + x282
    x508 = x0 * (x298 - x506 + x507)
    x509 = x293 + x506
    x510 = x4 * x507 - x5 * x509
    x511 = x4 * x509
    x512 = x183 * x420
    x513 = x287 + x512
    x514 = x5 * x513
    x515 = x511 - x514
    x516 = x0 * (x301 + x509 - x512)
    x517 = x198 * x420 + x303
    x518 = 60.3397786612521 * x70
    x519 = x215 * x420
    x520 = x225 * x420 + x305
    x521 = x0 * (x319 - x519 + x520)
    x522 = x314 + x519
    x523 = x4 * x520 - x5 * x522
    x524 = x4 * x522
    x525 = x220 * x420
    x526 = x309 + x525
    x527 = x5 * x526
    x528 = x524 - x527
    x529 = x0 * (x322 + x522 - x525)
    x530 = x234 * x420 + x324
    x531 = 104.511562358749 * x70
    x532 = x249 * x420
    x533 = x259 * x420 + x326
    x534 = x0 * (x334 - x532 + x533)
    x535 = x335 + x532
    x536 = x4 * x533 - x5 * x535
    x537 = x4 * x535
    x538 = x254 * x420
    x539 = x330 + x538
    x540 = x5 * x539
    x541 = x537 - x540
    x542 = x0 * (x329 + x535 - x538)
    x543 = x269 * x420 + x343
    x544 = x294 * x420
    x545 = x283 * x420 + x369
    x546 = x0 * (x373 - x544 + x545)
    x547 = x375 + x544
    x548 = x288 * x420 + x377
    x549 = x315 * x420
    x550 = x306 * x420 + x380
    x551 = x0 * (x384 - x549 + x550)
    x552 = x386 + x549
    x553 = x310 * x420 + x388
    x554 = x336 * x420
    x555 = x327 * x420 + x391
    x556 = x0 * (x394 - x554 + x555)
    x557 = x396 + x554
    x558 = x331 * x420 + x398
    x559 = x357 * x420
    x560 = x346 * x420 + x402
    x561 = x0 * (x400 - x559 + x560)
    x562 = x405 + x559
    x563 = x351 * x420 + x407
    x564 = 4.0 * x299 + x370 * x420 - x376 * x7
    x565 = 4.0 * x302 + x376 * x420 - x378 * x7
    x566 = 3.0 * x320 + x381 * x420 - x387 * x7
    x567 = 3.0 * x323 + x387 * x420 - x389 * x7
    x568 = 2.0 * x340 + x392 * x420 - x397 * x7
    x569 = 2.0 * x342 + x397 * x420 - x399 * x7
    x570 = x362 + x403 * x420 - x406 * x7
    x571 = x365 + x406 * x420 - x408 * x7
    x572 = x411 * x420 - x417 * x7
    x573 = x417 * x420 - x419 * x7
    x574 = -x8 - A[2]
    x575 = x13 * x574
    x576 = x148 + x39 * x574
    x577 = x0 * (x127 - x575 + x576)
    x578 = x128 + x575
    x579 = x4 * x576 - x5 * x578
    x580 = x4 * x578
    x581 = x12 * x574
    x582 = x133 + x581
    x583 = x5 * x582
    x584 = x580 - x583
    x585 = x4 * x579 - x5 * x584 + x577
    x586 = x0 * (x126 + x578 - x581)
    x587 = x4 * x584
    x588 = x4 * x582
    x589 = x19 * x574
    x590 = x140 + x589
    x591 = x5 * x590
    x592 = x588 - x591
    x593 = x5 * x592
    x594 = x586 + x587 - x593
    x595 = x35 * (x579 - x580 + x583) + x4 * x585 - x5 * x594
    x596 = x4 * x594
    x597 = x0 * (x139 + x582 - x589)
    x598 = x4 * x592
    x599 = x4 * x590
    x600 = x28 * x574
    x601 = x159 + x600
    x602 = x5 * x601
    x603 = x599 - x602
    x604 = x5 * x603
    x605 = x597 + x598 - x604
    x606 = x5 * x605
    x607 = x35 * (x584 - x588 + x591)
    x608 = x596 - x606 + x607
    x609 = -x586
    x610 = x0 * (x158 + x590 - x600)
    x611 = x172 + x51 * x574
    x612 = -x597
    x613 = x582 * x7
    x614 = x578 * x72
    x615 = -x578 * x7
    x616 = x576 * x72 + x615
    x617 = x0 * (x613 - x614 + x616)
    x618 = -x613
    x619 = x614 + x618
    x620 = x4 * x616 - x5 * x619
    x621 = x4 * x619
    x622 = x582 * x72
    x623 = x590 * x7
    x624 = -x623
    x625 = x622 + x624
    x626 = x5 * x625
    x627 = x621 - x626
    x628 = x4 * x620 - x5 * x627 + x617
    x629 = x0 * (x619 - x622 + x623)
    x630 = x4 * x627
    x631 = x4 * x625
    x632 = x590 * x72
    x633 = x601 * x7
    x634 = -x633
    x635 = x632 + x634
    x636 = x5 * x635
    x637 = x631 - x636
    x638 = x5 * x637
    x639 = x629 + x630 - x638
    x640 = x0 * (x625 - x632 + x633)
    x641 = -x611 * x7
    x642 = x601 * x72 + x641
    x643 = x130 * x574
    x644 = x149 * x574 + x258
    x645 = x0 * (x261 - x643 + x644)
    x646 = x248 + x643
    x647 = x4 * x644 - x5 * x646
    x648 = x4 * x646
    x649 = x134 * x574
    x650 = x253 + x649
    x651 = x5 * x650
    x652 = x648 - x651
    x653 = x4 * x647 - x5 * x652 + x645
    x654 = x0 * (x273 + x646 - x649)
    x655 = x4 * x652
    x656 = x4 * x650
    x657 = x142 * x574
    x658 = x268 + x657
    x659 = x5 * x658
    x660 = x656 - x659
    x661 = x5 * x660
    x662 = x654 + x655 - x661
    x663 = x0 * (x278 + x650 - x657)
    x664 = x161 * x574 + x280
    x665 = -x654
    x666 = x619 * x72
    x667 = x577 - x619 * x7
    x668 = x616 * x72 + x667
    x669 = x625 * x7
    x670 = x609 + x669
    x671 = x0 * (-x666 + x668 + x670)
    x672 = x586 - x669
    x673 = x666 + x672
    x674 = x4 * x668 - x5 * x673
    x675 = x4 * x673
    x676 = x625 * x72
    x677 = x635 * x7
    x678 = x597 - x677
    x679 = x676 + x678
    x680 = x5 * x679
    x681 = x675 - x680
    x682 = x612 + x677
    x683 = x0 * (x673 - x676 + x682)
    x684 = x610 - x642 * x7
    x685 = x635 * x72 + x684
    x686 = x650 * x7
    x687 = x646 * x72
    x688 = -x646 * x7
    x689 = x644 * x72 + x688
    x690 = x0 * (x686 - x687 + x689)
    x691 = -x686
    x692 = x687 + x691
    x693 = x4 * x689 - x5 * x692
    x694 = x4 * x692
    x695 = x650 * x72
    x696 = x658 * x7
    x697 = -x696
    x698 = x695 + x697
    x699 = x5 * x698
    x700 = x694 - x699
    x701 = x0 * (x692 - x695 + x696)
    x702 = -x664 * x7
    x703 = x658 * x72 + x702
    x704 = x249 * x574
    x705 = x259 * x574 + x345
    x706 = x0 * (x361 - x704 + x705)
    x707 = x356 + x704
    x708 = x4 * x705 - x5 * x707
    x709 = x4 * x707
    x710 = x254 * x574
    x711 = x350 + x710
    x712 = x5 * x711
    x713 = x709 - x712
    x714 = x0 * (x364 + x707 - x710)
    x715 = x269 * x574 + x366
    x716 = x673 * x72
    x717 = 2.0 * x617 - x673 * x7
    x718 = x668 * x72 + x717
    x719 = x679 * x7
    x720 = 2.0 * x629
    x721 = x719 - x720
    x722 = x0 * (-x716 + x718 + x721)
    x723 = -x719 + x720
    x724 = x716 + x723
    x725 = 2.0 * x640 - x685 * x7
    x726 = x679 * x72 + x725
    x727 = x692 * x72
    x728 = x645 - x692 * x7
    x729 = x689 * x72 + x728
    x730 = x698 * x7
    x731 = x665 + x730
    x732 = x0 * (-x727 + x729 + x731)
    x733 = x654 - x730
    x734 = x727 + x733
    x735 = x663 - x7 * x703
    x736 = x698 * x72 + x735
    x737 = x7 * x711
    x738 = x707 * x72
    x739 = -x7 * x707
    x740 = x705 * x72 + x739
    x741 = x0 * (x737 - x738 + x740)
    x742 = -x737
    x743 = x738 + x742
    x744 = -x7 * x715
    x745 = x711 * x72 + x744
    x746 = x357 * x574
    x747 = x346 * x574 + x410
    x748 = x0 * (x414 - x746 + x747)
    x749 = x416 + x746
    x750 = x351 * x574 + x418
    x751 = 3.0 * x671 - x7 * x724
    x752 = x718 * x72 + x751
    x753 = 3.0 * x683 - x7 * x726
    x754 = x72 * x724 + x753
    x755 = 2.0 * x690 - x7 * x734
    x756 = x72 * x729 + x755
    x757 = -x7 * x736 + 2.0 * x701
    x758 = x72 * x734 + x757
    x759 = -x7 * x743 + x706
    x760 = x72 * x740 + x759
    x761 = -x7 * x745 + x714
    x762 = x72 * x743 + x761
    x763 = -x7 * x749
    x764 = x72 * x747 + x763
    x765 = -x7 * x750
    x766 = x72 * x749 + x765
    x767 = 4.0 * x362 + x411 * x574 - x417 * x9
    x768 = 4.0 * x365 + x417 * x574 - x419 * x9
    x769 = x428 * x7
    x770 = x420 * x424
    x771 = x40 + x420 * x422 - x424 * x7
    x772 = x14 - x769 + x770
    x773 = x4 * x771 - x5 * x772
    x774 = x4 * x772
    x775 = x420 * x428
    x776 = x436 * x7
    x777 = x26 + x775 - x776
    x778 = x5 * x777
    x779 = x774 - x778
    x780 = x0 * (x44 + x769 - x770 + x771) + x4 * x773 - x5 * x779
    x781 = x0 * (x61 + x772 - x775 + x776)
    x782 = x4 * x779
    x783 = x4 * x777
    x784 = x420 * x436
    x785 = x447 * x7
    x786 = x49 + x784 - x785
    x787 = x5 * x786
    x788 = x783 - x787
    x789 = x5 * x788
    x790 = x781 + x782 - x789
    x791 = x467 * x7
    x792 = x420 * x463
    x793 = x420 * x461 + x423 - x463 * x7 + x98
    x794 = x432 + x79 - x791 + x792
    x795 = x4 * x793 - x5 * x794
    x796 = x4 * x794
    x797 = x420 * x467
    x798 = x475 * x7
    x799 = x103 + x443 + x797 - x798
    x800 = x5 * x799
    x801 = x796 - x800
    x802 = x490 * x7
    x803 = x420 * x486
    x804 = x150 + x420 * x484 - x486 * x7
    x805 = x131 - x802 + x803
    x806 = x4 * x804 - x5 * x805
    x807 = x4 * x805
    x808 = x420 * x490
    x809 = x498 * x7
    x810 = x155 + x808 - x809
    x811 = x5 * x810
    x812 = x807 - x811
    x813 = x513 * x7
    x814 = x420 * x509
    x815 = 2.0 * x471
    x816 = x191 + x420 * x507 + 2.0 * x462 - x509 * x7
    x817 = x203 - x813 + x814 + x815
    x818 = x526 * x7
    x819 = x420 * x522
    x820 = x227 + x420 * x520 + x485 - x522 * x7
    x821 = x238 + x494 - x818 + x819
    x822 = x539 * x7
    x823 = x420 * x535
    x824 = x262 + x420 * x533 - x535 * x7
    x825 = x274 - x822 + x823
    x826 = x420 * x578
    x827 = x420 * x576 + x615
    x828 = x618 + x826
    x829 = x4 * x827 - x5 * x828
    x830 = x4 * x828
    x831 = x420 * x582
    x832 = x624 + x831
    x833 = x5 * x832
    x834 = x830 - x833
    x835 = x0 * (x613 - x826 + x827) + x4 * x829 - x5 * x834
    x836 = x0 * (x623 + x828 - x831)
    x837 = x4 * x834
    x838 = x4 * x832
    x839 = x420 * x590
    x840 = x634 + x839
    x841 = x5 * x840
    x842 = x838 - x841
    x843 = x5 * x842
    x844 = x836 + x837 - x843
    x845 = x420 * x619
    x846 = x420 * x616 + x667
    x847 = x672 + x845
    x848 = x4 * x846 - x5 * x847
    x849 = x4 * x847
    x850 = x420 * x625
    x851 = x678 + x850
    x852 = x5 * x851
    x853 = x849 - x852
    x854 = x420 * x646
    x855 = x420 * x644 + x688
    x856 = x691 + x854
    x857 = x4 * x855 - x5 * x856
    x858 = x4 * x856
    x859 = x420 * x650
    x860 = x697 + x859
    x861 = x5 * x860
    x862 = x858 - x861
    x863 = x420 * x673
    x864 = x420 * x668 + x717
    x865 = x723 + x863
    x866 = x420 * x692
    x867 = x420 * x689 + x728
    x868 = x733 + x866
    x869 = x420 * x707
    x870 = x420 * x705 + x739
    x871 = x742 + x869
    x872 = x582 * x9
    x873 = x574 * x578
    x874 = x40 + x574 * x576 - x578 * x9
    x875 = x0 * (x44 + x872 - x873 + x874)
    x876 = x14 - x872 + x873
    x877 = x4 * x874 - x5 * x876
    x878 = x4 * x876
    x879 = x574 * x582
    x880 = x590 * x9
    x881 = x26 + x879 - x880
    x882 = x5 * x881
    x883 = x878 - x882
    x884 = x4 * x877 - x5 * x883 + x875
    x885 = x0 * (x61 + x876 - x879 + x880)
    x886 = x4 * x883
    x887 = x4 * x881
    x888 = x574 * x590
    x889 = x601 * x9
    x890 = x49 + x888 - x889
    x891 = x5 * x890
    x892 = x887 - x891
    x893 = x5 * x892
    x894 = x885 + x886 - x893
    x895 = x0 * (x69 + x881 - x888 + x889)
    x896 = x574 * x601 - x611 * x9 + x67
    x897 = -x885
    x898 = x7 * x881
    x899 = x72 * x876
    x900 = -x7 * x876 + x72 * x874
    x901 = x0 * (x898 - x899 + x900)
    x902 = -x898 + x899
    x903 = x4 * x900 - x5 * x902
    x904 = x4 * x902
    x905 = x72 * x881
    x906 = x7 * x890
    x907 = x905 - x906
    x908 = x5 * x907
    x909 = x904 - x908
    x910 = x0 * (x902 - x905 + x906)
    x911 = -x7 * x896 + x72 * x890
    x912 = x650 * x9
    x913 = x574 * x646
    x914 = x150 + x574 * x644 + x577 - x646 * x9
    x915 = x0 * (x170 + x609 + x912 - x913 + x914)
    x916 = x131 + x586 - x912 + x913
    x917 = x4 * x914 - x5 * x916
    x918 = x4 * x916
    x919 = x574 * x650
    x920 = x658 * x9
    x921 = x155 + x597 + x919 - x920
    x922 = x5 * x921
    x923 = x918 - x922
    x924 = x0 * (x174 + x612 + x916 - x919 + x920)
    x925 = x171 + x574 * x658 + x610 - x664 * x9
    x926 = x7 * x907
    x927 = x72 * x902
    x928 = -x7 * x902 + x72 * x900 + x875
    x929 = x0 * (x897 + x926 - x927 + x928)
    x930 = x885 - x926 + x927
    x931 = -x7 * x911 + x72 * x907 + x895
    x932 = x7 * x921
    x933 = x72 * x916
    x934 = -x7 * x916 + x72 * x914
    x935 = x0 * (x932 - x933 + x934)
    x936 = -x932 + x933
    x937 = -x7 * x925 + x72 * x921
    x938 = x711 * x9
    x939 = x574 * x707
    x940 = 2.0 * x654
    x941 = x262 + x574 * x705 + 2.0 * x645 - x707 * x9
    x942 = x0 * (x275 + x938 - x939 - x940 + x941)
    x943 = x274 - x938 + x939 + x940
    x944 = x279 + x574 * x711 + 2.0 * x663 - x715 * x9
    x945 = -x7 * x930 + x72 * x928 + 2.0 * x901
    x946 = -x7 * x931 + x72 * x930 + 2.0 * x910
    x947 = -x7 * x936 + x72 * x934 + x915
    x948 = -x7 * x937 + x72 * x936 + x924
    x949 = -x7 * x943 + x72 * x941
    x950 = -x7 * x944 + x72 * x943
    x951 = x362 + x574 * x747 + 3.0 * x706 - x749 * x9
    x952 = x365 + x574 * x749 + 3.0 * x714 - x750 * x9

    # 90 item(s)
    result[0, 0] = numpy.sum(
        x71
        * (
            x0 * (-x38 + x47 + x63)
            - x5
            * (
                -x5
                * (
                    x4 * x59
                    + x45 * (x33 - x50 + x55 + x69)
                    - x5
                    * (
                        x35 * (x31 + x52 - x53)
                        + x4 * x56
                        - x5 * (x4 * x54 - x5 * (x4 * x51 - x5 * x68) + x67)
                    )
                )
                + x64 * x65
                + x66 * (x37 - x48 + x57 - x58)
            )
            - x64 * (-x47 * x64 + x5 * x65 + x66 * (x25 - x34 + x36 - x43))
            + x66 * (-x37 * x64 + x43 * x64 + x46 + x63)
        )
    )
    result[0, 1] = numpy.sum(
        x123
        * (
            x0 * (x102 + x116 - x95)
            + x45 * (x100 * x64 + x101 + x116 - x64 * x94)
            - x5
            * (
                x117 * x64
                + x45 * (-x104 + x112 + x122 + x94)
                - x5
                * (
                    x113 * x4
                    + x35 * (-x105 + x110 + x92)
                    - x5 * (x111 * x4 + x119 - x5 * (x109 * x4 - x121 * x5))
                )
            )
            + x64 * (x102 * x64 - x117 * x5 + x45 * (x100 + x118 - x85 + x93))
        )
    )
    result[0, 2] = numpy.sum(
        x123
        * (
            x0 * (-x147 + x154 + x168)
            + x45 * (-x146 * x64 + x152 * x64 + x153 + x168)
            - x5
            * (
                x169 * x64
                + x45 * (x146 - x156 + x164 + x174)
                - x5
                * (
                    x165 * x4
                    + x35 * (x144 - x157 + x162)
                    - x5 * (x163 * x4 + x171 - x5 * (x161 * x4 - x173 * x5))
                )
            )
            + x64 * (x154 * x64 - x169 * x5 + x45 * (-x137 + x145 + x152 + x170))
        )
    )
    result[0, 3] = numpy.sum(
        x211
        * (
            x0 * (-x186 + x193 + x205)
            + x35 * (-x185 * x64 + x189 * x64 + x192 + x205)
            - x5
            * (
                x206 * x64
                + x35 * (x185 - x194 + x199)
                - x5 * (x200 * x4 + x208 - x5 * (x198 * x4 - x210 * x5))
            )
            + x64 * (x193 * x64 - x206 * x5 + x35 * (-x179 + x184 + x189))
        )
    )
    result[0, 4] = numpy.sum(
        x245
        * (
            x0 * (-x223 + x229 + x240)
            + x35 * (-x222 * x64 + x226 * x64 + x228 + x240)
            - x5
            * (
                x241 * x64
                + x35 * (x222 - x230 + x235)
                - x5 * (x236 * x4 + x242 - x5 * (x234 * x4 - x244 * x5))
            )
            + x64 * (x229 * x64 - x241 * x5 + x35 * (-x216 + x221 + x226))
        )
    )
    result[0, 5] = numpy.sum(
        x211
        * (
            x0 * (-x257 + x264 + x276)
            + x35 * (-x256 * x64 + x260 * x64 + x263 + x276)
            - x5
            * (
                x277 * x64
                + x35 * (x256 - x265 + x270)
                - x5 * (x271 * x4 + x279 - x5 * (x269 * x4 - x281 * x5))
            )
            + x64 * (x264 * x64 - x277 * x5 + x35 * (-x250 + x255 + x260))
        )
    )
    result[0, 6] = numpy.sum(
        x123
        * (
            x0 * (x289 - x296 + x297)
            + x0 * (x283 * x64 + x289 - x294 * x64 + x295)
            - x5 * (x300 * x64 + x302 - x5 * (x288 * x4 - x304 * x5))
            + x64 * (x297 * x64 + x299 - x300 * x5)
        )
    )
    result[0, 7] = numpy.sum(
        x245
        * (
            x0 * (x311 - x317 + x318)
            + x0 * (x306 * x64 + x311 - x315 * x64 + x316)
            - x5 * (x321 * x64 + x323 - x5 * (x310 * x4 - x325 * x5))
            + x64 * (x318 * x64 + x320 - x321 * x5)
        )
    )
    result[0, 8] = numpy.sum(
        x245
        * (
            x0 * (x332 - x338 + x339)
            + x0 * (x327 * x64 + x332 - x336 * x64 + x337)
            - x5 * (x341 * x64 + x342 - x5 * (x331 * x4 - x344 * x5))
            + x64 * (x339 * x64 + x340 - x341 * x5)
        )
    )
    result[0, 9] = numpy.sum(
        x123
        * (
            x0 * (x352 - x359 + x360)
            + x0 * (x346 * x64 + x352 - x357 * x64 + x358)
            - x5 * (x363 * x64 + x365 - x5 * (x351 * x4 - x367 * x5))
            + x64 * (x360 * x64 + x362 - x363 * x5)
        )
    )
    result[0, 10] = numpy.sum(
        x71 * (x374 - x5 * (x376 * x64 - x378 * x5) + x64 * (x370 * x64 - x376 * x5))
    )
    result[0, 11] = numpy.sum(
        x123 * (x385 - x5 * (x387 * x64 - x389 * x5) + x64 * (x381 * x64 - x387 * x5))
    )
    result[0, 12] = numpy.sum(
        x211 * (x395 - x5 * (x397 * x64 - x399 * x5) + x64 * (x392 * x64 - x397 * x5))
    )
    result[0, 13] = numpy.sum(
        x123 * (x404 - x5 * (x406 * x64 - x408 * x5) + x64 * (x403 * x64 - x406 * x5))
    )
    result[0, 14] = numpy.sum(
        x71 * (x415 - x5 * (x417 * x64 - x419 * x5) + x64 * (x411 * x64 - x417 * x5))
    )
    result[1, 0] = numpy.sum(
        x459
        * (
            -x5
            * (
                x4 * x454
                + x45 * (x440 - x444 + x450 + x458)
                - x5
                * (
                    x35 * (x438 - x445 + x448)
                    + x4 * x451
                    - x5 * (x4 * x449 + x456 - x5 * (x4 * x447 - x457 * x5))
                )
            )
            + x64 * (x4 * x441 + x45 * (x431 - x433 + x439 + x455) - x454 * x5)
            + x66 * (x441 - x442 + x452 - x453)
        )
    )
    result[1, 1] = numpy.sum(
        x482
        * (
            x45 * (x470 - x471 - x472 + x478)
            - x5
            * (
                x35 * (x469 - x473 + x476)
                + x4 * x479
                - x5 * (x4 * x477 + x480 - x5 * (x4 * x475 - x481 * x5))
            )
            + x64 * (x35 * (x464 - x465 + x468) + x4 * x470 - x479 * x5)
        )
    )
    result[1, 2] = numpy.sum(
        x482
        * (
            x45 * (x493 - x495 + x501 + x505)
            - x5
            * (
                x35 * (x492 - x496 + x499)
                + x4 * x502
                - x5 * (x4 * x500 - x5 * (x4 * x498 - x5 * x504) + x503)
            )
            + x64 * (x35 * (x487 - x488 + x491) + x4 * x493 - x5 * x502)
        )
    )
    result[1, 3] = numpy.sum(
        x518
        * (
            x35 * (x510 - x511 + x514)
            - x5 * (x4 * x515 - x5 * (x4 * x513 - x5 * x517) + x516)
            + x64 * (x4 * x510 - x5 * x515 + x508)
        )
    )
    result[1, 4] = numpy.sum(
        x531
        * (
            x35 * (x523 - x524 + x527)
            - x5 * (x4 * x528 - x5 * (x4 * x526 - x5 * x530) + x529)
            + x64 * (x4 * x523 - x5 * x528 + x521)
        )
    )
    result[1, 5] = numpy.sum(
        x518
        * (
            x35 * (x536 - x537 + x540)
            - x5 * (x4 * x541 - x5 * (x4 * x539 - x5 * x543) + x542)
            + x64 * (x4 * x536 - x5 * x541 + x534)
        )
    )
    result[1, 6] = numpy.sum(
        x482 * (-x5 * (x4 * x547 - x5 * x548) + x546 + x64 * (x4 * x545 - x5 * x547))
    )
    result[1, 7] = numpy.sum(
        x531 * (-x5 * (x4 * x552 - x5 * x553) + x551 + x64 * (x4 * x550 - x5 * x552))
    )
    result[1, 8] = numpy.sum(
        x531 * (-x5 * (x4 * x557 - x5 * x558) + x556 + x64 * (x4 * x555 - x5 * x557))
    )
    result[1, 9] = numpy.sum(
        x482 * (-x5 * (x4 * x562 - x5 * x563) + x561 + x64 * (x4 * x560 - x5 * x562))
    )
    result[1, 10] = numpy.sum(x459 * (-x5 * x565 + x564 * x64))
    result[1, 11] = numpy.sum(x482 * (-x5 * x567 + x566 * x64))
    result[1, 12] = numpy.sum(x518 * (-x5 * x569 + x568 * x64))
    result[1, 13] = numpy.sum(x482 * (-x5 * x571 + x570 * x64))
    result[1, 14] = numpy.sum(x459 * (-x5 * x573 + x572 * x64))
    result[2, 0] = numpy.sum(
        x459
        * (
            -x5
            * (
                x4 * x608
                + x45 * (x594 - x598 + x604 + x612)
                - x5
                * (
                    x35 * (x592 - x599 + x602)
                    + x4 * x605
                    - x5 * (x4 * x603 - x5 * (x4 * x601 - x5 * x611) + x610)
                )
            )
            + x64 * (x4 * x595 + x45 * (x585 - x587 + x593 + x609) - x5 * x608)
            + x66 * (x595 - x596 + x606 - x607)
        )
    )
    result[2, 1] = numpy.sum(
        x482
        * (
            x45 * (x628 - x629 - x630 + x638)
            - x5
            * (
                x35 * (x627 - x631 + x636)
                + x4 * x639
                - x5 * (x4 * x637 - x5 * (x4 * x635 - x5 * x642) + x640)
            )
            + x64 * (x35 * (x620 - x621 + x626) + x4 * x628 - x5 * x639)
        )
    )
    result[2, 2] = numpy.sum(
        x482
        * (
            x45 * (x653 - x655 + x661 + x665)
            - x5
            * (
                x35 * (x652 - x656 + x659)
                + x4 * x662
                - x5 * (x4 * x660 - x5 * (x4 * x658 - x5 * x664) + x663)
            )
            + x64 * (x35 * (x647 - x648 + x651) + x4 * x653 - x5 * x662)
        )
    )
    result[2, 3] = numpy.sum(
        x518
        * (
            x35 * (x674 - x675 + x680)
            - x5 * (x4 * x681 - x5 * (x4 * x679 - x5 * x685) + x683)
            + x64 * (x4 * x674 - x5 * x681 + x671)
        )
    )
    result[2, 4] = numpy.sum(
        x531
        * (
            x35 * (x693 - x694 + x699)
            - x5 * (x4 * x700 - x5 * (x4 * x698 - x5 * x703) + x701)
            + x64 * (x4 * x693 - x5 * x700 + x690)
        )
    )
    result[2, 5] = numpy.sum(
        x518
        * (
            x35 * (x708 - x709 + x712)
            - x5 * (x4 * x713 - x5 * (x4 * x711 - x5 * x715) + x714)
            + x64 * (x4 * x708 - x5 * x713 + x706)
        )
    )
    result[2, 6] = numpy.sum(
        x482 * (-x5 * (x4 * x724 - x5 * x726) + x64 * (x4 * x718 - x5 * x724) + x722)
    )
    result[2, 7] = numpy.sum(
        x531 * (-x5 * (x4 * x734 - x5 * x736) + x64 * (x4 * x729 - x5 * x734) + x732)
    )
    result[2, 8] = numpy.sum(
        x531 * (-x5 * (x4 * x743 - x5 * x745) + x64 * (x4 * x740 - x5 * x743) + x741)
    )
    result[2, 9] = numpy.sum(
        x482 * (-x5 * (x4 * x749 - x5 * x750) + x64 * (x4 * x747 - x5 * x749) + x748)
    )
    result[2, 10] = numpy.sum(x459 * (-x5 * x754 + x64 * x752))
    result[2, 11] = numpy.sum(x482 * (-x5 * x758 + x64 * x756))
    result[2, 12] = numpy.sum(x518 * (-x5 * x762 + x64 * x760))
    result[2, 13] = numpy.sum(x482 * (-x5 * x766 + x64 * x764))
    result[2, 14] = numpy.sum(x459 * (-x5 * x768 + x64 * x767))
    result[3, 0] = numpy.sum(
        x71
        * (
            x4 * (x35 * (x773 - x774 + x778) + x4 * x780 - x5 * x790)
            + x45 * (x780 - x781 - x782 + x789)
            - x5
            * (
                x35 * (x779 - x783 + x787)
                + x4 * x790
                - x5
                * (
                    x0 * (x69 + x777 - x784 + x785)
                    + x4 * x788
                    - x5 * (x4 * x786 - x5 * (x420 * x447 - x457 * x7 + x67))
                )
            )
        )
    )
    result[3, 1] = numpy.sum(
        x123
        * (
            x35 * (x795 - x796 + x800)
            + x4 * (x0 * (x118 + x455 + x791 - x792 + x793) + x4 * x795 - x5 * x801)
            - x5
            * (
                x0 * (x122 + x458 + x794 - x797 + x798)
                + x4 * x801
                - x5 * (x4 * x799 - x5 * (x119 + x420 * x475 + x456 - x481 * x7))
            )
        )
    )
    result[3, 2] = numpy.sum(
        x123
        * (
            x35 * (x806 - x807 + x811)
            + x4 * (x0 * (x170 + x802 - x803 + x804) + x4 * x806 - x5 * x812)
            - x5
            * (
                x0 * (x174 + x805 - x808 + x809)
                + x4 * x812
                - x5 * (x4 * x810 - x5 * (x171 + x420 * x498 - x504 * x7))
            )
        )
    )
    result[3, 3] = numpy.sum(
        x211
        * (
            x0 * (x204 + x813 - x814 - x815 + x816)
            + x4 * (x4 * x816 - x5 * x817)
            - x5 * (x4 * x817 - x5 * (x208 + x420 * x513 + 2.0 * x480 - x517 * x7))
        )
    )
    result[3, 4] = numpy.sum(
        x245
        * (
            x0 * (x239 + x505 + x818 - x819 + x820)
            + x4 * (x4 * x820 - x5 * x821)
            - x5 * (x4 * x821 - x5 * (x242 + x420 * x526 + x503 - x530 * x7))
        )
    )
    result[3, 5] = numpy.sum(
        x211
        * (
            x0 * (x275 + x822 - x823 + x824)
            + x4 * (x4 * x824 - x5 * x825)
            - x5 * (x4 * x825 - x5 * (x279 + x420 * x539 - x543 * x7))
        )
    )
    result[3, 6] = numpy.sum(
        x123
        * (
            x4 * (x299 + x420 * x545 + 3.0 * x508 - x547 * x7)
            - x5 * (x302 + x420 * x547 + 3.0 * x516 - x548 * x7)
        )
    )
    result[3, 7] = numpy.sum(
        x245
        * (
            x4 * (x320 + x420 * x550 + 2.0 * x521 - x552 * x7)
            - x5 * (x323 + x420 * x552 + 2.0 * x529 - x553 * x7)
        )
    )
    result[3, 8] = numpy.sum(
        x245
        * (
            x4 * (x340 + x420 * x555 + x534 - x557 * x7)
            - x5 * (x342 + x420 * x557 + x542 - x558 * x7)
        )
    )
    result[3, 9] = numpy.sum(
        x123
        * (x4 * (x362 + x420 * x560 - x562 * x7) - x5 * (x365 + x420 * x562 - x563 * x7))
    )
    result[3, 10] = numpy.sum(x71 * (x374 + x420 * x564 + 4.0 * x546 - x565 * x7))
    result[3, 11] = numpy.sum(x123 * (x385 + x420 * x566 + 3.0 * x551 - x567 * x7))
    result[3, 12] = numpy.sum(x211 * (x395 + x420 * x568 + 2.0 * x556 - x569 * x7))
    result[3, 13] = numpy.sum(x123 * (x404 + x420 * x570 + x561 - x571 * x7))
    result[3, 14] = numpy.sum(x71 * (x415 + x420 * x572 - x573 * x7))
    result[4, 0] = numpy.sum(
        x459
        * (
            x4 * (x35 * (x829 - x830 + x833) + x4 * x835 - x5 * x844)
            + x45 * (x835 - x836 - x837 + x843)
            - x5
            * (
                x35 * (x834 - x838 + x841)
                + x4 * x844
                - x5
                * (
                    x0 * (x633 + x832 - x839)
                    + x4 * x842
                    - x5 * (x4 * x840 - x5 * (x420 * x601 + x641))
                )
            )
        )
    )
    result[4, 1] = numpy.sum(
        x482
        * (
            x35 * (x848 - x849 + x852)
            + x4 * (x0 * (x670 - x845 + x846) + x4 * x848 - x5 * x853)
            - x5
            * (
                x0 * (x682 + x847 - x850)
                + x4 * x853
                - x5 * (x4 * x851 - x5 * (x420 * x635 + x684))
            )
        )
    )
    result[4, 2] = numpy.sum(
        x482
        * (
            x35 * (x857 - x858 + x861)
            + x4 * (x0 * (x686 - x854 + x855) + x4 * x857 - x5 * x862)
            - x5
            * (
                x0 * (x696 + x856 - x859)
                + x4 * x862
                - x5 * (x4 * x860 - x5 * (x420 * x658 + x702))
            )
        )
    )
    result[4, 3] = numpy.sum(
        x518
        * (
            x0 * (x721 - x863 + x864)
            + x4 * (x4 * x864 - x5 * x865)
            - x5 * (x4 * x865 - x5 * (x420 * x679 + x725))
        )
    )
    result[4, 4] = numpy.sum(
        x531
        * (
            x0 * (x731 - x866 + x867)
            + x4 * (x4 * x867 - x5 * x868)
            - x5 * (x4 * x868 - x5 * (x420 * x698 + x735))
        )
    )
    result[4, 5] = numpy.sum(
        x518
        * (
            x0 * (x737 - x869 + x870)
            + x4 * (x4 * x870 - x5 * x871)
            - x5 * (x4 * x871 - x5 * (x420 * x711 + x744))
        )
    )
    result[4, 6] = numpy.sum(
        x482 * (x4 * (x420 * x718 + x751) - x5 * (x420 * x724 + x753))
    )
    result[4, 7] = numpy.sum(
        x531 * (x4 * (x420 * x729 + x755) - x5 * (x420 * x734 + x757))
    )
    result[4, 8] = numpy.sum(
        x531 * (x4 * (x420 * x740 + x759) - x5 * (x420 * x743 + x761))
    )
    result[4, 9] = numpy.sum(
        x482 * (x4 * (x420 * x747 + x763) - x5 * (x420 * x749 + x765))
    )
    result[4, 10] = numpy.sum(x459 * (x420 * x752 - x7 * x754 + 4.0 * x722))
    result[4, 11] = numpy.sum(x482 * (x420 * x756 - x7 * x758 + 3.0 * x732))
    result[4, 12] = numpy.sum(x518 * (x420 * x760 - x7 * x762 + 2.0 * x741))
    result[4, 13] = numpy.sum(x482 * (x420 * x764 - x7 * x766 + x748))
    result[4, 14] = numpy.sum(x459 * (x420 * x767 - x7 * x768))
    result[5, 0] = numpy.sum(
        x71
        * (
            x4 * (x35 * (x877 - x878 + x882) + x4 * x884 - x5 * x894)
            + x45 * (x884 - x886 + x893 + x897)
            - x5
            * (
                x35 * (x883 - x887 + x891)
                + x4 * x894
                - x5 * (x4 * x892 - x5 * (x4 * x890 - x5 * x896) + x895)
            )
        )
    )
    result[5, 1] = numpy.sum(
        x123
        * (
            x35 * (x903 - x904 + x908)
            + x4 * (x4 * x903 - x5 * x909 + x901)
            - x5 * (x4 * x909 - x5 * (x4 * x907 - x5 * x911) + x910)
        )
    )
    result[5, 2] = numpy.sum(
        x123
        * (
            x35 * (x917 - x918 + x922)
            + x4 * (x4 * x917 - x5 * x923 + x915)
            - x5 * (x4 * x923 - x5 * (x4 * x921 - x5 * x925) + x924)
        )
    )
    result[5, 3] = numpy.sum(
        x211 * (x4 * (x4 * x928 - x5 * x930) - x5 * (x4 * x930 - x5 * x931) + x929)
    )
    result[5, 4] = numpy.sum(
        x245 * (x4 * (x4 * x934 - x5 * x936) - x5 * (x4 * x936 - x5 * x937) + x935)
    )
    result[5, 5] = numpy.sum(
        x211 * (x4 * (x4 * x941 - x5 * x943) - x5 * (x4 * x943 - x5 * x944) + x942)
    )
    result[5, 6] = numpy.sum(x123 * (x4 * x945 - x5 * x946))
    result[5, 7] = numpy.sum(x245 * (x4 * x947 - x5 * x948))
    result[5, 8] = numpy.sum(x245 * (x4 * x949 - x5 * x950))
    result[5, 9] = numpy.sum(x123 * (x4 * x951 - x5 * x952))
    result[5, 10] = numpy.sum(x71 * (-x7 * x946 + x72 * x945 + 3.0 * x929))
    result[5, 11] = numpy.sum(x123 * (-x7 * x948 + x72 * x947 + 2.0 * x935))
    result[5, 12] = numpy.sum(x211 * (-x7 * x950 + x72 * x949 + x942))
    result[5, 13] = numpy.sum(x123 * (-x7 * x952 + x72 * x951))
    result[5, 14] = numpy.sum(x71 * (x415 + x574 * x767 + 4.0 * x748 - x768 * x9))
    return result


def coulomb3d_30(ax, da, A, bx, db, B, C):
    """Cartesian (fs) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 1), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = 0.5 / (ax + bx)
    x5 = -x2 - C[0]
    x6 = -x1 * (ax * A[1] + bx * B[1])
    x7 = -x6 - C[1]
    x8 = -x1 * (ax * A[2] + bx * B[2])
    x9 = -x8 - C[2]
    x10 = x0 * (x5**2 + x7**2 + x9**2)
    x11 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x12 = x11 * boys(1, x10)
    x13 = x11 * boys(0, x10)
    x14 = x4 * (-x12 + x13)
    x15 = -x12 * x5 + x13 * x3
    x16 = x11 * boys(2, x10)
    x17 = x16 * x5
    x18 = x12 * x3
    x19 = -x17 + x18
    x20 = x4 * (x12 - x16)
    x21 = x11 * boys(3, x10)
    x22 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**1.5)
    x23 = 5.84237394672177 * x22
    x24 = -x6 - A[1]
    x25 = x12 * x24
    x26 = x16 * x7
    x27 = -x12 * x7 + x13 * x24
    x28 = x4 * (-x25 + x26 + x27)
    x29 = x25 - x26
    x30 = x16 * x24 - x21 * x7
    x31 = 13.0639452948436 * x22
    x32 = -x8 - A[2]
    x33 = x12 * x32
    x34 = x16 * x9
    x35 = -x12 * x9 + x13 * x32
    x36 = x4 * (-x33 + x34 + x35)
    x37 = x33 - x34
    x38 = x16 * x32 - x21 * x9
    x39 = x14 + x24 * x27 - x29 * x7
    x40 = x20 + x24 * x29 - x30 * x7
    x41 = x24 * x35 - x37 * x7
    x42 = x24 * x37 - x38 * x7
    x43 = x14 + x32 * x35 - x37 * x9
    x44 = x20 + x32 * x37 - x38 * x9

    # 10 item(s)
    result[0, 0] = numpy.sum(
        x23
        * (
            x3 * (x14 + x15 * x3 - x19 * x5)
            + 2.0 * x4 * (x15 + x17 - x18)
            - x5 * (x19 * x3 + x20 - x5 * (x16 * x3 - x21 * x5))
        )
    )
    result[1, 0] = numpy.sum(
        x31 * (x28 + x3 * (x27 * x3 - x29 * x5) - x5 * (x29 * x3 - x30 * x5))
    )
    result[2, 0] = numpy.sum(
        x31 * (x3 * (x3 * x35 - x37 * x5) + x36 - x5 * (x3 * x37 - x38 * x5))
    )
    result[3, 0] = numpy.sum(x31 * (x3 * x39 - x40 * x5))
    result[4, 0] = numpy.sum(22.6274169979695 * x22 * (x3 * x41 - x42 * x5))
    result[5, 0] = numpy.sum(x31 * (x3 * x43 - x44 * x5))
    result[6, 0] = numpy.sum(x23 * (x24 * x39 + 2.0 * x28 - x40 * x7))
    result[7, 0] = numpy.sum(x31 * (x24 * x41 + x36 - x42 * x7))
    result[8, 0] = numpy.sum(x31 * (x24 * x43 - x44 * x7))
    result[9, 0] = numpy.sum(x23 * (x32 * x43 + 2.0 * x36 - x44 * x9))
    return result


def coulomb3d_31(ax, da, A, bx, db, B, C):
    """Cartesian (fp) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 3), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = -x2 * (ax * A[0] + bx * B[0])
    x4 = -x3 - C[0]
    x5 = -x2 * (ax * A[1] + bx * B[1])
    x6 = -x5 - C[1]
    x7 = -x2 * (ax * A[2] + bx * B[2])
    x8 = -x7 - C[2]
    x9 = x1 * (x4**2 + x6**2 + x8**2)
    x10 = (
        6.28318530717959
        * x2
        * numpy.exp(
            -ax * bx * x2 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x11 = x10 * boys(1, x9)
    x12 = x10 * boys(0, x9)
    x13 = x0 * (-x11 + x12)
    x14 = -x3 - A[0]
    x15 = -x11 * x4
    x16 = x12 * x14 + x15
    x17 = x10 * boys(3, x9)
    x18 = x17 * x4
    x19 = -x18
    x20 = x10 * boys(2, x9)
    x21 = x14 * x20
    x22 = x0 * (x11 - x20)
    x23 = -x22
    x24 = x20 * x4
    x25 = -x24
    x26 = x11 * x14
    x27 = x25 + x26
    x28 = -x3 - B[0]
    x29 = x11 * x28
    x30 = x12 * x28 + x15
    x31 = x25 + x29
    x32 = x13 + x14 * x30 - x31 * x4
    x33 = x14 * x31
    x34 = x20 * x28
    x35 = x19 + x34
    x36 = x35 * x4
    x37 = x22 + x33 - x36
    x38 = x0 * (-x17 + x20)
    x39 = x10 * boys(4, x9)
    x40 = 2.0 * x0
    x41 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**2.5)
    x42 = 11.6847478934435 * x41
    x43 = -x5 - B[1]
    x44 = x11 * x43
    x45 = x20 * x6
    x46 = -x11 * x6
    x47 = x12 * x43 + x46
    x48 = x0 * (-x44 + x45 + x47)
    x49 = -x45
    x50 = x44 + x49
    x51 = x14 * x47 - x4 * x50
    x52 = x14 * x50
    x53 = x17 * x6
    x54 = -x53
    x55 = x20 * x43
    x56 = x54 + x55
    x57 = x4 * x56
    x58 = x52 - x57
    x59 = x0 * (x50 + x53 - x55)
    x60 = -x39 * x6
    x61 = x17 * x43 + x60
    x62 = -x7 - B[2]
    x63 = x11 * x62
    x64 = x20 * x8
    x65 = -x11 * x8
    x66 = x12 * x62 + x65
    x67 = x0 * (-x63 + x64 + x66)
    x68 = -x64
    x69 = x63 + x68
    x70 = x14 * x66 - x4 * x69
    x71 = x14 * x69
    x72 = x17 * x8
    x73 = -x72
    x74 = x20 * x62
    x75 = x73 + x74
    x76 = x4 * x75
    x77 = x71 - x76
    x78 = x0 * (x69 + x72 - x74)
    x79 = -x39 * x8
    x80 = x17 * x62 + x79
    x81 = -x5 - A[1]
    x82 = x12 * x81 + x46
    x83 = x20 * x81
    x84 = x54 + x83
    x85 = x4 * x84
    x86 = x11 * x81
    x87 = x49 + x86
    x88 = -x4 * x87
    x89 = x28 * x87
    x90 = x28 * x82 + x88
    x91 = x0 * (x45 + x82 - x86)
    x92 = -x85 + x89
    x93 = x0 * (x53 - x83 + x87)
    x94 = x17 * x81 + x60
    x95 = 26.1278905896872 * x41
    x96 = x56 * x6
    x97 = x50 * x81
    x98 = x13 + x47 * x81 - x50 * x6
    x99 = x0 * (x23 + x96 - x97 + x98)
    x100 = x22 - x96 + x97
    x101 = x38 + x56 * x81 - x6 * x61
    x102 = x6 * x75
    x103 = x69 * x81
    x104 = -x6 * x69 + x66 * x81
    x105 = x0 * (x102 - x103 + x104)
    x106 = -x102 + x103
    x107 = -x6 * x80 + x75 * x81
    x108 = -x7 - A[2]
    x109 = x108 * x12 + x65
    x110 = x108 * x20
    x111 = x110 + x73
    x112 = x111 * x4
    x113 = x108 * x11
    x114 = x113 + x68
    x115 = -x114 * x4
    x116 = x114 * x28
    x117 = x109 * x28 + x115
    x118 = x0 * (x109 - x113 + x64)
    x119 = -x112 + x116
    x120 = x0 * (-x110 + x114 + x72)
    x121 = x108 * x17 + x79
    x122 = x111 * x6
    x123 = x114 * x43
    x124 = -x114 * x6
    x125 = x109 * x43 + x124
    x126 = x0 * (x122 - x123 + x125)
    x127 = -x122
    x128 = x123 + x127
    x129 = -x121 * x6
    x130 = x111 * x43 + x129
    x131 = x75 * x8
    x132 = x108 * x69
    x133 = x108 * x66 + x13 - x69 * x8
    x134 = x0 * (x131 - x132 + x133 + x23)
    x135 = -x131 + x132 + x22
    x136 = x108 * x75 + x38 - x8 * x80
    x137 = x6 * x84
    x138 = x81 * x87
    x139 = x13 - x6 * x87 + x81 * x82
    x140 = x0 * (x137 - x138 + x139 + x23)
    x141 = -x137 + x138 + x22
    x142 = x38 - x6 * x94 + x81 * x84
    x143 = -x100 * x6 + x48 + x81 * x98 + x91
    x144 = x100 * x81 - x101 * x6 + x59 + x93
    x145 = x104 * x81 - x106 * x6 + x67
    x146 = x106 * x81 - x107 * x6 + x78
    x147 = x114 * x81
    x148 = x109 * x81 + x124
    x149 = x0 * (x122 - x147 + x148)
    x150 = x127 + x147
    x151 = x111 * x81 + x129
    x152 = 45.2548339959391 * x41
    x153 = x118 + x125 * x81 - x128 * x6
    x154 = x120 + x128 * x81 - x130 * x6
    x155 = x133 * x81 - x135 * x6
    x156 = x135 * x81 - x136 * x6
    x157 = x111 * x8
    x158 = x108 * x114
    x159 = x108 * x109 - x114 * x8 + x13
    x160 = x0 * (x157 - x158 + x159 + x23)
    x161 = -x157 + x158 + x22
    x162 = x108 * x111 - x121 * x8 + x38
    x163 = -x161 * x6
    x164 = x159 * x43 + x163
    x165 = -x162 * x6
    x166 = x161 * x43 + x165
    x167 = x108 * x133 + x118 - x135 * x8 + x67
    x168 = x108 * x135 + x120 - x136 * x8 + x78
    x169 = x108 * x159 + 2.0 * x118 - x161 * x8
    x170 = x108 * x161 + 2.0 * x120 - x162 * x8

    # 30 item(s)
    result[0, 0] = numpy.sum(
        x42
        * (
            x0 * (x13 + x14 * x16 - x14 * x27 + x23 - x27 * x4 + x4 * (x19 + x21))
            + x14
            * (x0 * (x16 + x24 - x26) + x0 * (x24 - x29 + x30) + x14 * x32 - x37 * x4)
            - x4
            * (
                x0 * (x18 - x21 + x27)
                + x0 * (x18 + x31 - x34)
                + x14 * x37
                - x4 * (x14 * x35 + x38 - x4 * (x17 * x28 - x39 * x4))
            )
            + x40 * (x23 + x32 - x33 + x36)
        )
    )
    result[0, 1] = numpy.sum(
        x42
        * (
            x14 * (x14 * x51 - x4 * x58 + x48)
            - x4 * (x14 * x58 - x4 * (x14 * x56 - x4 * x61) + x59)
            + x40 * (x51 - x52 + x57)
        )
    )
    result[0, 2] = numpy.sum(
        x42
        * (
            x14 * (x14 * x70 - x4 * x77 + x67)
            - x4 * (x14 * x77 - x4 * (x14 * x75 - x4 * x80) + x78)
            + x40 * (x70 - x71 + x76)
        )
    )
    result[1, 0] = numpy.sum(
        x95
        * (
            x0 * (x85 - x89 + x90)
            + x0 * (x14 * x82 - x14 * x87 + x85 + x88)
            + x14 * (x14 * x90 - x4 * x92 + x91)
            - x4 * (x14 * x92 - x4 * (x28 * x84 - x4 * x94) + x93)
        )
    )
    result[1, 1] = numpy.sum(
        x95 * (-2.0 * x100 * x14 * x4 + x101 * x4**2 + x14**2 * x98 + x99)
    )
    result[1, 2] = numpy.sum(
        x95 * (x105 + x14 * (x104 * x14 - x106 * x4) - x4 * (x106 * x14 - x107 * x4))
    )
    result[2, 0] = numpy.sum(
        x95
        * (
            x0 * (x112 - x116 + x117)
            + x0 * (x109 * x14 + x112 - x114 * x14 + x115)
            + x14 * (x117 * x14 + x118 - x119 * x4)
            - x4 * (x119 * x14 + x120 - x4 * (x111 * x28 - x121 * x4))
        )
    )
    result[2, 1] = numpy.sum(
        x95 * (x126 + x14 * (x125 * x14 - x128 * x4) - x4 * (x128 * x14 - x130 * x4))
    )
    result[2, 2] = numpy.sum(
        x95 * (x134 + x14 * (x133 * x14 - x135 * x4) - x4 * (x135 * x14 - x136 * x4))
    )
    result[3, 0] = numpy.sum(
        x95 * (x14 * (x139 * x28 - x141 * x4) + x140 - x4 * (x141 * x28 - x142 * x4))
    )
    result[3, 1] = numpy.sum(x95 * (x14 * x143 - x144 * x4))
    result[3, 2] = numpy.sum(x95 * (x14 * x145 - x146 * x4))
    result[4, 0] = numpy.sum(
        x152 * (x14 * (x148 * x28 - x150 * x4) + x149 - x4 * (x150 * x28 - x151 * x4))
    )
    result[4, 1] = numpy.sum(x152 * (x14 * x153 - x154 * x4))
    result[4, 2] = numpy.sum(x152 * (x14 * x155 - x156 * x4))
    result[5, 0] = numpy.sum(
        x95 * (x14 * (x159 * x28 - x161 * x4) + x160 - x4 * (x161 * x28 - x162 * x4))
    )
    result[5, 1] = numpy.sum(x95 * (x14 * x164 - x166 * x4))
    result[5, 2] = numpy.sum(x95 * (x14 * x167 - x168 * x4))
    result[6, 0] = numpy.sum(
        x42
        * (
            x28 * (x139 * x81 - x141 * x6 + 2.0 * x91)
            - x4 * (x141 * x81 - x142 * x6 + 2.0 * x93)
        )
    )
    result[6, 1] = numpy.sum(x42 * (x140 + x143 * x81 - x144 * x6 + 2.0 * x99))
    result[6, 2] = numpy.sum(x42 * (2.0 * x105 + x145 * x81 - x146 * x6))
    result[7, 0] = numpy.sum(
        x95
        * (x28 * (x118 + x148 * x81 - x150 * x6) - x4 * (x120 + x150 * x81 - x151 * x6))
    )
    result[7, 1] = numpy.sum(x95 * (x126 + x149 + x153 * x81 - x154 * x6))
    result[7, 2] = numpy.sum(x95 * (x134 + x155 * x81 - x156 * x6))
    result[8, 0] = numpy.sum(x95 * (x28 * (x159 * x81 + x163) - x4 * (x161 * x81 + x165)))
    result[8, 1] = numpy.sum(x95 * (x160 + x164 * x81 - x166 * x6))
    result[8, 2] = numpy.sum(x95 * (x167 * x81 - x168 * x6))
    result[9, 0] = numpy.sum(x42 * (x169 * x28 - x170 * x4))
    result[9, 1] = numpy.sum(x42 * (x169 * x43 - x170 * x6))
    result[9, 2] = numpy.sum(x42 * (x108 * x167 + 2.0 * x134 + x160 - x168 * x8))
    return result


def coulomb3d_32(ax, da, A, bx, db, B, C):
    """Cartesian (fd) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 6), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = 0.5 / (ax + bx)
    x5 = -x2 - B[0]
    x6 = -x2 - C[0]
    x7 = -x1 * (ax * A[1] + bx * B[1])
    x8 = -x7 - C[1]
    x9 = -x1 * (ax * A[2] + bx * B[2])
    x10 = -x9 - C[2]
    x11 = x0 * (x10**2 + x6**2 + x8**2)
    x12 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x13 = x12 * boys(2, x11)
    x14 = x13 * x6
    x15 = -x14
    x16 = x12 * boys(1, x11)
    x17 = x16 * x5
    x18 = x15 + x17
    x19 = x18 * x5
    x20 = -x16 * x6
    x21 = x12 * boys(0, x11)
    x22 = x20 + x21 * x5
    x23 = x4 * (-x16 + x21)
    x24 = -x18 * x6 + x23
    x25 = x22 * x5 + x24
    x26 = x12 * boys(3, x11)
    x27 = x26 * x6
    x28 = x13 * x5
    x29 = -x27 + x28
    x30 = x29 * x6
    x31 = x4 * (-x13 + x16)
    x32 = -x31
    x33 = x30 + x32
    x34 = -x30 + x31
    x35 = x19 + x34
    x36 = x4 * (x14 - x17 + x22)
    x37 = x25 * x3 - x35 * x6 + 2.0 * x36
    x38 = x3 * x35
    x39 = x29 * x5
    x40 = x4 * (x13 - x26)
    x41 = x12 * boys(4, x11)
    x42 = x41 * x6
    x43 = x26 * x5
    x44 = -x42 + x43
    x45 = x44 * x6
    x46 = x40 - x45
    x47 = x39 + x46
    x48 = x47 * x6
    x49 = x4 * (x18 + x27 - x28)
    x50 = 2.0 * x49
    x51 = x38 - x48 + x50
    x52 = x18 * x3
    x53 = x22 * x3 + x24
    x54 = 2.0 * x4
    x55 = -x40
    x56 = x45 + x55
    x57 = x4 * (x26 - x41)
    x58 = x12 * boys(5, x11)
    x59 = x29 * x3
    x60 = x34 + x52
    x61 = x16 * x3
    x62 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**3.5)
    x63 = 13.4923846833851 * x62
    x64 = -x7 - B[1]
    x65 = x16 * x64
    x66 = x13 * x8
    x67 = -x16 * x8
    x68 = x21 * x64 + x67
    x69 = x4 * (-x65 + x66 + x68)
    x70 = -x66
    x71 = x65 + x70
    x72 = -x6 * x71
    x73 = x3 * x68 + x72
    x74 = x26 * x8
    x75 = -x74
    x76 = x13 * x64
    x77 = x75 + x76
    x78 = x3 * x77
    x79 = x41 * x8
    x80 = -x79
    x81 = x26 * x64
    x82 = x80 + x81
    x83 = x6 * x82
    x84 = -x83
    x85 = x4 * (x71 + x74 - x76)
    x86 = -x85
    x87 = x3 * x71
    x88 = x6 * x77
    x89 = -x88
    x90 = x87 + x89
    x91 = x5 * x71
    x92 = x5 * x68 + x72
    x93 = x89 + x91
    x94 = x3 * x92 - x6 * x93 + x69
    x95 = x3 * x93
    x96 = x5 * x77
    x97 = x84 + x96
    x98 = x6 * x97
    x99 = x85 + x95 - x98
    x100 = x4 * (x77 + x79 - x81)
    x101 = -x58 * x8
    x102 = x101 + x41 * x64
    x103 = 23.3694957868871 * x62
    x104 = -x9 - B[2]
    x105 = x104 * x16
    x106 = x10 * x13
    x107 = -x10 * x16
    x108 = x104 * x21 + x107
    x109 = x4 * (-x105 + x106 + x108)
    x110 = -x106
    x111 = x105 + x110
    x112 = -x111 * x6
    x113 = x108 * x3 + x112
    x114 = x10 * x26
    x115 = -x114
    x116 = x104 * x13
    x117 = x115 + x116
    x118 = x117 * x3
    x119 = x10 * x41
    x120 = -x119
    x121 = x104 * x26
    x122 = x120 + x121
    x123 = x122 * x6
    x124 = -x123
    x125 = x4 * (x111 + x114 - x116)
    x126 = -x125
    x127 = x111 * x3
    x128 = x117 * x6
    x129 = -x128
    x130 = x127 + x129
    x131 = x111 * x5
    x132 = x108 * x5 + x112
    x133 = x129 + x131
    x134 = x109 + x132 * x3 - x133 * x6
    x135 = x133 * x3
    x136 = x117 * x5
    x137 = x124 + x136
    x138 = x137 * x6
    x139 = x125 + x135 - x138
    x140 = x4 * (x117 + x119 - x121)
    x141 = -x10 * x58
    x142 = x104 * x41 + x141
    x143 = x64 * x71
    x144 = x23 - x71 * x8
    x145 = x144 + x64 * x68
    x146 = x77 * x8
    x147 = x146 + x32
    x148 = x4 * (-x143 + x145 + x147)
    x149 = -x146 + x31
    x150 = x143 + x149
    x151 = x145 * x3 - x150 * x6
    x152 = x150 * x3
    x153 = x64 * x77
    x154 = x8 * x82
    x155 = -x154 + x40
    x156 = x153 + x155
    x157 = x156 * x6
    x158 = x152 - x157
    x159 = x154 + x55
    x160 = x4 * (x150 - x153 + x159)
    x161 = -x102 * x8 + x57
    x162 = x161 + x64 * x82
    x163 = x117 * x8
    x164 = x111 * x64
    x165 = -x111 * x8
    x166 = x108 * x64 + x165
    x167 = x4 * (x163 - x164 + x166)
    x168 = -x163
    x169 = x164 + x168
    x170 = x166 * x3 - x169 * x6
    x171 = x169 * x3
    x172 = x117 * x64
    x173 = x122 * x8
    x174 = -x173
    x175 = x172 + x174
    x176 = x175 * x6
    x177 = x171 - x176
    x178 = x4 * (x169 - x172 + x173)
    x179 = -x142 * x8
    x180 = x122 * x64 + x179
    x181 = x104 * x111
    x182 = -x10 * x111 + x23
    x183 = x104 * x108 + x182
    x184 = x10 * x117
    x185 = x184 + x32
    x186 = x4 * (-x181 + x183 + x185)
    x187 = -x184 + x31
    x188 = x181 + x187
    x189 = x183 * x3 - x188 * x6
    x190 = x188 * x3
    x191 = x104 * x117
    x192 = x10 * x122
    x193 = -x192 + x40
    x194 = x191 + x193
    x195 = x194 * x6
    x196 = x190 - x195
    x197 = x192 + x55
    x198 = x4 * (x188 - x191 + x197)
    x199 = -x10 * x142 + x57
    x200 = x104 * x122 + x199
    x201 = -x7 - A[1]
    x202 = x16 * x201
    x203 = x202 + x70
    x204 = x203 * x5
    x205 = x13 * x201
    x206 = x205 + x75
    x207 = x206 * x6
    x208 = x204 - x207
    x209 = x208 * x5
    x210 = x201 * x21 + x67
    x211 = -x203 * x6 + x210 * x5
    x212 = x4 * (-x202 + x210 + x66)
    x213 = -x208 * x6 + x212
    x214 = x211 * x5 + x213
    x215 = x206 * x5
    x216 = x201 * x26
    x217 = x216 + x80
    x218 = x217 * x6
    x219 = x215 - x218
    x220 = x219 * x6
    x221 = x4 * (x203 - x205 + x74)
    x222 = -x221
    x223 = x220 + x222
    x224 = x209 - x220 + x221
    x225 = x4 * (x206 - x216 + x79)
    x226 = x101 + x201 * x41
    x227 = 30.169889330626 * x62
    x228 = x144 + x201 * x68
    x229 = x201 * x77
    x230 = x155 + x229
    x231 = x230 * x6
    x232 = x201 * x71
    x233 = x149 + x232
    x234 = -x233 * x6
    x235 = x233 * x5
    x236 = x228 * x5 + x234
    x237 = x4 * (x147 + x228 - x232)
    x238 = -x231 + x235
    x239 = x4 * (x159 - x229 + x233)
    x240 = x161 + x201 * x82
    x241 = 52.2557811793745 * x62
    x242 = x108 * x201 + x165
    x243 = x117 * x201
    x244 = x174 + x243
    x245 = x244 * x6
    x246 = x111 * x201
    x247 = x168 + x246
    x248 = -x247 * x6
    x249 = x247 * x5
    x250 = x242 * x5 + x248
    x251 = x4 * (x163 + x242 - x246)
    x252 = -x245 + x249
    x253 = x4 * (x173 - x243 + x247)
    x254 = x122 * x201 + x179
    x255 = x156 * x8
    x256 = x150 * x201
    x257 = 2.0 * x85
    x258 = x145 * x201 - x150 * x8 + 2.0 * x69
    x259 = x4 * (x255 - x256 - x257 + x258)
    x260 = -x255 + x256 + x257
    x261 = 2.0 * x100 + x156 * x201 - x162 * x8
    x262 = x175 * x8
    x263 = x169 * x201
    x264 = x109 + x166 * x201 - x169 * x8
    x265 = x4 * (x126 + x262 - x263 + x264)
    x266 = x125 - x262 + x263
    x267 = x140 + x175 * x201 - x180 * x8
    x268 = x194 * x8
    x269 = x188 * x201
    x270 = x183 * x201 - x188 * x8
    x271 = x4 * (x268 - x269 + x270)
    x272 = -x268 + x269
    x273 = x194 * x201 - x200 * x8
    x274 = -x9 - A[2]
    x275 = x16 * x274
    x276 = x110 + x275
    x277 = x276 * x5
    x278 = x13 * x274
    x279 = x115 + x278
    x280 = x279 * x6
    x281 = x277 - x280
    x282 = x281 * x5
    x283 = x107 + x21 * x274
    x284 = -x276 * x6 + x283 * x5
    x285 = x4 * (x106 - x275 + x283)
    x286 = -x281 * x6 + x285
    x287 = x284 * x5 + x286
    x288 = x279 * x5
    x289 = x26 * x274
    x290 = x120 + x289
    x291 = x290 * x6
    x292 = x288 - x291
    x293 = x292 * x6
    x294 = x4 * (x114 + x276 - x278)
    x295 = -x294
    x296 = x293 + x295
    x297 = x282 - x293 + x294
    x298 = x4 * (x119 + x279 - x289)
    x299 = x141 + x274 * x41
    x300 = -x276 * x8
    x301 = x283 * x64 + x300
    x302 = x279 * x64
    x303 = x290 * x8
    x304 = -x303
    x305 = x302 + x304
    x306 = x305 * x6
    x307 = x276 * x64
    x308 = x279 * x8
    x309 = -x308
    x310 = x307 + x309
    x311 = -x310 * x6
    x312 = x310 * x5
    x313 = x301 * x5 + x311
    x314 = x4 * (x301 - x307 + x308)
    x315 = -x306 + x312
    x316 = x4 * (-x302 + x303 + x310)
    x317 = -x299 * x8
    x318 = x290 * x64 + x317
    x319 = x108 * x274 + x182
    x320 = x117 * x274
    x321 = x193 + x320
    x322 = x321 * x6
    x323 = x111 * x274
    x324 = x187 + x323
    x325 = -x324 * x6
    x326 = x324 * x5
    x327 = x319 * x5 + x325
    x328 = x4 * (x185 + x319 - x323)
    x329 = -x322 + x326
    x330 = x4 * (x197 - x320 + x324)
    x331 = x122 * x274 + x199
    x332 = x310 * x64
    x333 = x285 - x310 * x8
    x334 = x301 * x64 + x333
    x335 = x305 * x8
    x336 = x295 + x335
    x337 = x4 * (-x332 + x334 + x336)
    x338 = x294 - x335
    x339 = x332 + x338
    x340 = x298 - x318 * x8
    x341 = x305 * x64 + x340
    x342 = x321 * x8
    x343 = x324 * x64
    x344 = -x324 * x8
    x345 = x319 * x64 + x344
    x346 = x4 * (x342 - x343 + x345)
    x347 = -x342
    x348 = x343 + x347
    x349 = -x331 * x8
    x350 = x321 * x64 + x349
    x351 = x10 * x194
    x352 = x188 * x274
    x353 = 2.0 * x125
    x354 = -x10 * x188 + 2.0 * x109 + x183 * x274
    x355 = x4 * (x351 - x352 - x353 + x354)
    x356 = -x351 + x352 + x353
    x357 = -x10 * x200 + 2.0 * x140 + x194 * x274
    x358 = x206 * x8
    x359 = x201 * x203
    x360 = x201 * x210 - x203 * x8 + x23
    x361 = x4 * (x32 + x358 - x359 + x360)
    x362 = x31 - x358 + x359
    x363 = x360 * x5 - x362 * x6
    x364 = x362 * x5
    x365 = x201 * x206
    x366 = x217 * x8
    x367 = x365 - x366 + x40
    x368 = x367 * x6
    x369 = x364 - x368
    x370 = x4 * (x362 - x365 + x366 + x55)
    x371 = x201 * x217 - x226 * x8 + x57
    x372 = x230 * x8
    x373 = x201 * x233
    x374 = x201 * x228 + x212 - x233 * x8 + x69
    x375 = x4 * (x222 + x372 - x373 + x374 + x86)
    x376 = x221 - x372 + x373 + x85
    x377 = x100 + x201 * x230 + x225 - x240 * x8
    x378 = x244 * x8
    x379 = x201 * x247
    x380 = x109 + x201 * x242 - x247 * x8
    x381 = x4 * (x126 + x378 - x379 + x380)
    x382 = x125 - x378 + x379
    x383 = x140 + x201 * x244 - x254 * x8
    x384 = 2.0 * x237
    x385 = x148 + x201 * x258 - x260 * x8 + x384
    x386 = 2.0 * x239
    x387 = x160 + x201 * x260 - x261 * x8 + x386
    x388 = x167 + x201 * x264 + x251 - x266 * x8
    x389 = x178 + x201 * x266 + x253 - x267 * x8
    x390 = x186 + x201 * x270 - x272 * x8
    x391 = x198 + x201 * x272 - x273 * x8
    x392 = x201 * x276
    x393 = x201 * x283 + x300
    x394 = x4 * (x308 - x392 + x393)
    x395 = x309 + x392
    x396 = x393 * x5 - x395 * x6
    x397 = x395 * x5
    x398 = x201 * x279
    x399 = x304 + x398
    x400 = x399 * x6
    x401 = x397 - x400
    x402 = x4 * (x303 + x395 - x398)
    x403 = x201 * x290 + x317
    x404 = x201 * x310
    x405 = x201 * x301 + x333
    x406 = x4 * (x336 - x404 + x405)
    x407 = x338 + x404
    x408 = x201 * x305 + x340
    x409 = 90.5096679918781 * x62
    x410 = x201 * x324
    x411 = x201 * x319 + x344
    x412 = x4 * (x342 - x410 + x411)
    x413 = x347 + x410
    x414 = x201 * x321 + x349
    x415 = x201 * x334 + 2.0 * x314 - x339 * x8
    x416 = x201 * x339 + 2.0 * x316 - x341 * x8
    x417 = x201 * x345 + x328 - x348 * x8
    x418 = x201 * x348 + x330 - x350 * x8
    x419 = x201 * x354 - x356 * x8
    x420 = x201 * x356 - x357 * x8
    x421 = x10 * x279
    x422 = x274 * x276
    x423 = -x10 * x276 + x23 + x274 * x283
    x424 = x4 * (x32 + x421 - x422 + x423)
    x425 = x31 - x421 + x422
    x426 = x423 * x5 - x425 * x6
    x427 = x425 * x5
    x428 = x274 * x279
    x429 = x10 * x290
    x430 = x40 + x428 - x429
    x431 = x430 * x6
    x432 = x427 - x431
    x433 = x4 * (x425 - x428 + x429 + x55)
    x434 = -x10 * x299 + x274 * x290 + x57
    x435 = x430 * x8
    x436 = x425 * x64
    x437 = -x425 * x8
    x438 = x423 * x64 + x437
    x439 = x4 * (x435 - x436 + x438)
    x440 = -x435
    x441 = x436 + x440
    x442 = -x434 * x8
    x443 = x430 * x64 + x442
    x444 = x10 * x321
    x445 = x274 * x324
    x446 = -x10 * x324 + x109 + x274 * x319 + x285
    x447 = x4 * (x126 + x295 + x444 - x445 + x446)
    x448 = x125 + x294 - x444 + x445
    x449 = -x10 * x331 + x140 + x274 * x321 + x298
    x450 = x424 - x441 * x8
    x451 = x438 * x64 + x450
    x452 = x433 - x443 * x8
    x453 = x441 * x64 + x452
    x454 = -x448 * x8
    x455 = x446 * x64 + x454
    x456 = -x449 * x8
    x457 = x448 * x64 + x456
    x458 = 2.0 * x328
    x459 = -x10 * x356 + x186 + x274 * x354 + x458
    x460 = 2.0 * x330
    x461 = -x10 * x357 + x198 + x274 * x356 + x460
    x462 = x367 * x8
    x463 = x201 * x362
    x464 = 2.0 * x221
    x465 = x201 * x360 + 2.0 * x212 - x362 * x8
    x466 = -x462 + x463 + x464
    x467 = x399 * x8
    x468 = x201 * x395
    x469 = x201 * x393 + x285 - x395 * x8
    x470 = x294 - x467 + x468
    x471 = x201 * x425
    x472 = x201 * x423 + x437
    x473 = x440 + x471
    x474 = x10 * x430
    x475 = x274 * x425
    x476 = 2.0 * x294
    x477 = -x10 * x425 + x274 * x423 + 2.0 * x285
    x478 = x4 * (x474 - x475 - x476 + x477)
    x479 = -x474 + x475 + x476
    x480 = -x10 * x434 + x274 * x430 + 2.0 * x298
    x481 = x477 * x64 - x479 * x8
    x482 = x479 * x64 - x480 * x8
    x483 = -x10 * x448 + x274 * x446 + x424 + x458
    x484 = -x10 * x449 + x274 * x448 + x433 + x460

    # 60 item(s)
    result[0, 0] = numpy.sum(
        x63
        * (
            x3 * (x3 * x37 + x4 * (-x19 + x25 + x33) - x51 * x6 + x54 * (x33 - x52 + x53))
            + x54 * (x37 - x38 + x48 - x50)
            + x54
            * (
                x3 * x53
                - x3 * x60
                + x36
                + x4 * (x14 + x20 + x21 * x3 - x61)
                - x4 * (-x13 * x3 + x15 + x27 + x61)
                - x49
                - x6 * x60
                + x6 * (x46 + x59)
            )
            - x6
            * (
                x3 * x51
                + x4 * (x35 - x39 + x56)
                + x54 * (x56 - x59 + x60)
                - x6
                * (
                    x3 * x47
                    + x54 * (x29 + x42 - x43)
                    - x6 * (x44 * x5 + x57 - x6 * (x41 * x5 - x58 * x6))
                )
            )
        )
    )
    result[0, 1] = numpy.sum(
        x103
        * (
            x3 * (x3 * x94 + x4 * (x73 - x87 + x88) + x4 * (x88 - x91 + x92) - x6 * x99)
            + x4 * (x3 * x73 - x3 * x90 - x6 * x90 + x6 * (x78 + x84) + x69 + x86)
            + x54 * (x86 + x94 - x95 + x98)
            - x6
            * (
                x3 * x99
                + x4 * (-x78 + x83 + x90)
                + x4 * (x83 + x93 - x96)
                - x6 * (x100 + x3 * x97 + x6 * (x102 * x6 - x5 * x82))
            )
        )
    )
    result[0, 2] = numpy.sum(
        x103
        * (
            x3
            * (
                x134 * x3
                - x139 * x6
                + x4 * (x113 - x127 + x128)
                + x4 * (x128 - x131 + x132)
            )
            + x4 * (x109 + x113 * x3 + x126 - x130 * x3 - x130 * x6 + x6 * (x118 + x124))
            + x54 * (x126 + x134 - x135 + x138)
            - x6
            * (
                x139 * x3
                + x4 * (-x118 + x123 + x130)
                + x4 * (x123 + x133 - x136)
                - x6 * (x137 * x3 + x140 - x6 * (x122 * x5 - x142 * x6))
            )
        )
    )
    result[0, 3] = numpy.sum(
        x63
        * (
            x3 * (x148 + x151 * x3 - x158 * x6)
            + x54 * (x151 - x152 + x157)
            - x6 * (x158 * x3 + x160 - x6 * (x156 * x3 - x162 * x6))
        )
    )
    result[0, 4] = numpy.sum(
        x103
        * (
            x3 * (x167 + x170 * x3 - x177 * x6)
            + x54 * (x170 - x171 + x176)
            - x6 * (x177 * x3 + x178 - x6 * (x175 * x3 - x180 * x6))
        )
    )
    result[0, 5] = numpy.sum(
        x63
        * (
            x3 * (x186 + x189 * x3 - x196 * x6)
            + x54 * (x189 - x190 + x195)
            - x6 * (x196 * x3 + x198 - x6 * (x194 * x3 - x200 * x6))
        )
    )
    result[1, 0] = numpy.sum(
        x227
        * (
            x3 * (x214 * x3 - x224 * x6 + x54 * (-x204 + x207 + x211))
            + x4 * (-x209 + x214 + x223)
            + x54 * (-x208 * x3 + x211 * x3 + x213 + x223)
            - x6
            * (
                x224 * x3
                + x54 * (x208 - x215 + x218)
                - x6 * (x219 * x5 + x225 - x6 * (x217 * x5 - x226 * x6))
            )
        )
    )
    result[1, 1] = numpy.sum(
        x241
        * (
            x3 * (x236 * x3 + x237 - x238 * x6)
            + x4 * (x231 - x235 + x236)
            + x4 * (x228 * x3 + x231 - x233 * x3 + x234)
            - x6 * (x238 * x3 + x239 - x6 * (x230 * x5 - x240 * x6))
        )
    )
    result[1, 2] = numpy.sum(
        x241
        * (
            x3 * (x250 * x3 + x251 - x252 * x6)
            + x4 * (x245 - x249 + x250)
            + x4 * (x242 * x3 + x245 - x247 * x3 + x248)
            - x6 * (x252 * x3 + x253 - x6 * (x244 * x5 - x254 * x6))
        )
    )
    result[1, 3] = numpy.sum(
        x227 * (x259 + x3 * (x258 * x3 - x260 * x6) - x6 * (x260 * x3 - x261 * x6))
    )
    result[1, 4] = numpy.sum(
        x241 * (x265 + x3 * (x264 * x3 - x266 * x6) - x6 * (x266 * x3 - x267 * x6))
    )
    result[1, 5] = numpy.sum(
        x227 * (x271 + x3 * (x270 * x3 - x272 * x6) - x6 * (x272 * x3 - x273 * x6))
    )
    result[2, 0] = numpy.sum(
        x227
        * (
            x3 * (x287 * x3 - x297 * x6 + x54 * (-x277 + x280 + x284))
            + x4 * (-x282 + x287 + x296)
            + x54 * (-x281 * x3 + x284 * x3 + x286 + x296)
            - x6
            * (
                x297 * x3
                + x54 * (x281 - x288 + x291)
                - x6 * (x292 * x5 + x298 - x6 * (x290 * x5 - x299 * x6))
            )
        )
    )
    result[2, 1] = numpy.sum(
        x241
        * (
            x3 * (x3 * x313 + x314 - x315 * x6)
            + x4 * (x306 - x312 + x313)
            + x4 * (x3 * x301 - x3 * x310 + x306 + x311)
            - x6 * (x3 * x315 + x316 - x6 * (x305 * x5 - x318 * x6))
        )
    )
    result[2, 2] = numpy.sum(
        x241
        * (
            x3 * (x3 * x327 + x328 - x329 * x6)
            + x4 * (x322 - x326 + x327)
            + x4 * (x3 * x319 - x3 * x324 + x322 + x325)
            - x6 * (x3 * x329 + x330 - x6 * (x321 * x5 - x331 * x6))
        )
    )
    result[2, 3] = numpy.sum(
        x227 * (x3 * (x3 * x334 - x339 * x6) + x337 - x6 * (x3 * x339 - x341 * x6))
    )
    result[2, 4] = numpy.sum(
        x241 * (x3 * (x3 * x345 - x348 * x6) + x346 - x6 * (x3 * x348 - x350 * x6))
    )
    result[2, 5] = numpy.sum(
        x227 * (x3 * (x3 * x354 - x356 * x6) + x355 - x6 * (x3 * x356 - x357 * x6))
    )
    result[3, 0] = numpy.sum(
        x227
        * (
            x3 * (x361 + x363 * x5 - x369 * x6)
            + x54 * (x363 - x364 + x368)
            - x6 * (x369 * x5 + x370 - x6 * (x367 * x5 - x371 * x6))
        )
    )
    result[3, 1] = numpy.sum(
        x241 * (x3 * (x374 * x5 - x376 * x6) + x375 - x6 * (x376 * x5 - x377 * x6))
    )
    result[3, 2] = numpy.sum(
        x241 * (x3 * (x380 * x5 - x382 * x6) + x381 - x6 * (x382 * x5 - x383 * x6))
    )
    result[3, 3] = numpy.sum(x227 * (x3 * x385 - x387 * x6))
    result[3, 4] = numpy.sum(x241 * (x3 * x388 - x389 * x6))
    result[3, 5] = numpy.sum(x227 * (x3 * x390 - x391 * x6))
    result[4, 0] = numpy.sum(
        x241
        * (
            x3 * (x394 + x396 * x5 - x401 * x6)
            + x54 * (x396 - x397 + x400)
            - x6 * (x401 * x5 + x402 - x6 * (x399 * x5 - x403 * x6))
        )
    )
    result[4, 1] = numpy.sum(
        x409 * (x3 * (x405 * x5 - x407 * x6) + x406 - x6 * (x407 * x5 - x408 * x6))
    )
    result[4, 2] = numpy.sum(
        x409 * (x3 * (x411 * x5 - x413 * x6) + x412 - x6 * (x413 * x5 - x414 * x6))
    )
    result[4, 3] = numpy.sum(x241 * (x3 * x415 - x416 * x6))
    result[4, 4] = numpy.sum(x409 * (x3 * x417 - x418 * x6))
    result[4, 5] = numpy.sum(x241 * (x3 * x419 - x420 * x6))
    result[5, 0] = numpy.sum(
        x227
        * (
            x3 * (x424 + x426 * x5 - x432 * x6)
            + x54 * (x426 - x427 + x431)
            - x6 * (x432 * x5 + x433 - x6 * (x430 * x5 - x434 * x6))
        )
    )
    result[5, 1] = numpy.sum(
        x241 * (x3 * (x438 * x5 - x441 * x6) + x439 - x6 * (x441 * x5 - x443 * x6))
    )
    result[5, 2] = numpy.sum(
        x241 * (x3 * (x446 * x5 - x448 * x6) + x447 - x6 * (x448 * x5 - x449 * x6))
    )
    result[5, 3] = numpy.sum(x227 * (x3 * x451 - x453 * x6))
    result[5, 4] = numpy.sum(x241 * (x3 * x455 - x457 * x6))
    result[5, 5] = numpy.sum(x227 * (x3 * x459 - x461 * x6))
    result[6, 0] = numpy.sum(
        x63
        * (
            x4 * (x462 - x463 - x464 + x465)
            + x5 * (x465 * x5 - x466 * x6)
            - x6 * (x466 * x5 - x6 * (x201 * x367 + 2.0 * x225 - x371 * x8))
        )
    )
    result[6, 1] = numpy.sum(
        x103
        * (
            x5 * (x201 * x374 + x361 - x376 * x8 + x384)
            - x6 * (x201 * x376 + x370 - x377 * x8 + x386)
        )
    )
    result[6, 2] = numpy.sum(
        x103
        * (
            x5 * (x201 * x380 + 2.0 * x251 - x382 * x8)
            - x6 * (x201 * x382 + 2.0 * x253 - x383 * x8)
        )
    )
    result[6, 3] = numpy.sum(x63 * (x201 * x385 + 2.0 * x259 + 2.0 * x375 - x387 * x8))
    result[6, 4] = numpy.sum(x103 * (x201 * x388 + 2.0 * x265 + x381 - x389 * x8))
    result[6, 5] = numpy.sum(x63 * (x201 * x390 + 2.0 * x271 - x391 * x8))
    result[7, 0] = numpy.sum(
        x227
        * (
            x4 * (x295 + x467 - x468 + x469)
            + x5 * (x469 * x5 - x470 * x6)
            - x6 * (x470 * x5 - x6 * (x201 * x399 + x298 - x403 * x8))
        )
    )
    result[7, 1] = numpy.sum(
        x241
        * (
            x5 * (x201 * x405 + x314 + x394 - x407 * x8)
            - x6 * (x201 * x407 + x316 + x402 - x408 * x8)
        )
    )
    result[7, 2] = numpy.sum(
        x241
        * (x5 * (x201 * x411 + x328 - x413 * x8) - x6 * (x201 * x413 + x330 - x414 * x8))
    )
    result[7, 3] = numpy.sum(x227 * (x201 * x415 + x337 + 2.0 * x406 - x416 * x8))
    result[7, 4] = numpy.sum(x241 * (x201 * x417 + x346 + x412 - x418 * x8))
    result[7, 5] = numpy.sum(x227 * (x201 * x419 + x355 - x420 * x8))
    result[8, 0] = numpy.sum(
        x227
        * (
            x4 * (x435 - x471 + x472)
            + x5 * (x472 * x5 - x473 * x6)
            - x6 * (x473 * x5 - x6 * (x201 * x430 + x442))
        )
    )
    result[8, 1] = numpy.sum(
        x241 * (x5 * (x201 * x438 + x450) - x6 * (x201 * x441 + x452))
    )
    result[8, 2] = numpy.sum(
        x241 * (x5 * (x201 * x446 + x454) - x6 * (x201 * x448 + x456))
    )
    result[8, 3] = numpy.sum(x227 * (x201 * x451 + 2.0 * x439 - x453 * x8))
    result[8, 4] = numpy.sum(x241 * (x201 * x455 + x447 - x457 * x8))
    result[8, 5] = numpy.sum(x227 * (x201 * x459 - x461 * x8))
    result[9, 0] = numpy.sum(
        x63 * (x478 + x5 * (x477 * x5 - x479 * x6) - x6 * (x479 * x5 - x480 * x6))
    )
    result[9, 1] = numpy.sum(x103 * (x481 * x5 - x482 * x6))
    result[9, 2] = numpy.sum(x103 * (x483 * x5 - x484 * x6))
    result[9, 3] = numpy.sum(x63 * (x478 + x481 * x64 - x482 * x8))
    result[9, 4] = numpy.sum(x103 * (x483 * x64 - x484 * x8))
    result[9, 5] = numpy.sum(x63 * (-x10 * x461 + x274 * x459 + 2.0 * x355 + 2.0 * x447))
    return result


def coulomb3d_33(ax, da, A, bx, db, B, C):
    """Cartesian (ff) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 10), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = 0.5 / (ax + bx)
    x5 = -x2 - B[0]
    x6 = -x2 - C[0]
    x7 = -x1 * (ax * A[1] + bx * B[1])
    x8 = -x7 - C[1]
    x9 = -x1 * (ax * A[2] + bx * B[2])
    x10 = -x9 - C[2]
    x11 = x0 * (x10**2 + x6**2 + x8**2)
    x12 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x13 = x12 * boys(2, x11)
    x14 = x13 * x6
    x15 = x12 * boys(1, x11)
    x16 = x15 * x5
    x17 = -x14 + x16
    x18 = x17 * x5
    x19 = x4 * (-x13 + x15)
    x20 = x12 * boys(3, x11)
    x21 = x20 * x6
    x22 = x13 * x5
    x23 = -x21 + x22
    x24 = x23 * x6
    x25 = x19 - x24
    x26 = x18 + x25
    x27 = x26 * x5
    x28 = x12 * boys(0, x11)
    x29 = -x15 * x6 + x28 * x5
    x30 = x4 * (-x15 + x28)
    x31 = -x17 * x6 + x30
    x32 = x29 * x5 + x31
    x33 = 2.0 * x4
    x34 = -x26 * x6 + x33 * (x14 - x16 + x29)
    x35 = x32 * x5 + x34
    x36 = x4 * (x13 - x20)
    x37 = x23 * x5
    x38 = x12 * boys(4, x11)
    x39 = x38 * x6
    x40 = x20 * x5
    x41 = -x39 + x40
    x42 = x41 * x6
    x43 = x36 + x37 - x42
    x44 = x43 * x6
    x45 = x33 * (x17 + x21 - x22)
    x46 = x44 - x45
    x47 = -x44 + x45
    x48 = x27 + x47
    x49 = -x19
    x50 = x24 + x49
    x51 = x4 * (-x18 + x32 + x50)
    x52 = x3 * x35 - x48 * x6 + 3.0 * x51
    x53 = x3 * x48
    x54 = x43 * x5
    x55 = x4 * (x20 - x38)
    x56 = x41 * x5
    x57 = x12 * boys(5, x11)
    x58 = x57 * x6
    x59 = x38 * x5
    x60 = -x58 + x59
    x61 = x6 * x60
    x62 = x55 + x56 - x61
    x63 = x6 * x62
    x64 = x33 * (x23 + x39 - x40)
    x65 = -x63 + x64
    x66 = x54 + x65
    x67 = x6 * x66
    x68 = -x36
    x69 = x42 + x68
    x70 = x4 * (x26 - x37 + x69)
    x71 = 3.0 * x70
    x72 = x53 - x67 + x71
    x73 = x26 * x3
    x74 = x3 * x32 + x34
    x75 = 3.0 * x4
    x76 = x63 - x64
    x77 = x4 * (x38 - x57)
    x78 = x12 * boys(6, x11)
    x79 = -x55
    x80 = x3 * x43
    x81 = x47 + x73
    x82 = x17 * x3
    x83 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**4.5)
    x84 = 12.0679557322504 * x83
    x85 = x13 * x8
    x86 = -x85
    x87 = -x7 - B[1]
    x88 = x15 * x87
    x89 = x86 + x88
    x90 = x5 * x89
    x91 = x20 * x8
    x92 = -x91
    x93 = x13 * x87
    x94 = x92 + x93
    x95 = x6 * x94
    x96 = -x95
    x97 = x90 + x96
    x98 = x5 * x97
    x99 = -x15 * x8
    x100 = x28 * x87 + x99
    x101 = -x6 * x89
    x102 = x100 * x5 + x101
    x103 = x4 * (x100 + x85 - x88)
    x104 = x103 - x6 * x97
    x105 = x102 * x5 + x104
    x106 = x5 * x94
    x107 = x38 * x8
    x108 = -x107
    x109 = x20 * x87
    x110 = x108 + x109
    x111 = x110 * x6
    x112 = x106 - x111
    x113 = x112 * x6
    x114 = x4 * (x89 + x91 - x93)
    x115 = -x114
    x116 = x113 + x115
    x117 = -x113 + x114
    x118 = x117 + x98
    x119 = x4 * (x102 - x90 + x95)
    x120 = x105 * x3 - x118 * x6 + 2.0 * x119
    x121 = x118 * x3
    x122 = x112 * x5
    x123 = x4 * (x107 - x109 + x94)
    x124 = x110 * x5
    x125 = x57 * x8
    x126 = -x125
    x127 = x38 * x87
    x128 = x126 + x127
    x129 = x128 * x6
    x130 = x124 - x129
    x131 = x130 * x6
    x132 = x123 - x131
    x133 = x122 + x132
    x134 = x133 * x6
    x135 = x4 * (-x106 + x111 + x97)
    x136 = 2.0 * x135
    x137 = x121 - x134 + x136
    x138 = x3 * x97
    x139 = x102 * x3 + x104
    x140 = -x123
    x141 = x131 + x140
    x142 = x4 * (x110 + x125 - x127)
    x143 = -x78 * x8
    x144 = x143 + x57 * x87
    x145 = x112 * x3
    x146 = x117 + x138
    x147 = x3 * x89
    x148 = 26.9847693667702 * x83
    x149 = x10 * x13
    x150 = -x149
    x151 = -x9 - B[2]
    x152 = x15 * x151
    x153 = x150 + x152
    x154 = x153 * x5
    x155 = x10 * x20
    x156 = -x155
    x157 = x13 * x151
    x158 = x156 + x157
    x159 = x158 * x6
    x160 = -x159
    x161 = x154 + x160
    x162 = x161 * x5
    x163 = -x10 * x15
    x164 = x151 * x28 + x163
    x165 = -x153 * x6
    x166 = x164 * x5 + x165
    x167 = x4 * (x149 - x152 + x164)
    x168 = -x161 * x6 + x167
    x169 = x166 * x5 + x168
    x170 = x158 * x5
    x171 = x10 * x38
    x172 = -x171
    x173 = x151 * x20
    x174 = x172 + x173
    x175 = x174 * x6
    x176 = x170 - x175
    x177 = x176 * x6
    x178 = x4 * (x153 + x155 - x157)
    x179 = -x178
    x180 = x177 + x179
    x181 = -x177 + x178
    x182 = x162 + x181
    x183 = x4 * (-x154 + x159 + x166)
    x184 = x169 * x3 - x182 * x6 + 2.0 * x183
    x185 = x182 * x3
    x186 = x176 * x5
    x187 = x4 * (x158 + x171 - x173)
    x188 = x174 * x5
    x189 = x10 * x57
    x190 = -x189
    x191 = x151 * x38
    x192 = x190 + x191
    x193 = x192 * x6
    x194 = x188 - x193
    x195 = x194 * x6
    x196 = x187 - x195
    x197 = x186 + x196
    x198 = x197 * x6
    x199 = x4 * (x161 - x170 + x175)
    x200 = 2.0 * x199
    x201 = x185 - x198 + x200
    x202 = x161 * x3
    x203 = x166 * x3 + x168
    x204 = -x187
    x205 = x195 + x204
    x206 = x4 * (x174 + x189 - x191)
    x207 = -x10 * x78
    x208 = x151 * x57 + x207
    x209 = x176 * x3
    x210 = x181 + x202
    x211 = x153 * x3
    x212 = x87 * x89
    x213 = x30 - x8 * x89
    x214 = x100 * x87 + x213
    x215 = x8 * x94
    x216 = x215 + x49
    x217 = x4 * (-x212 + x214 + x216)
    x218 = x19 - x215
    x219 = x212 + x218
    x220 = -x219 * x6
    x221 = x214 * x3 + x220
    x222 = x87 * x94
    x223 = x110 * x8
    x224 = -x223 + x36
    x225 = x222 + x224
    x226 = x225 * x3
    x227 = x110 * x87
    x228 = x128 * x8
    x229 = -x228 + x55
    x230 = x227 + x229
    x231 = x230 * x6
    x232 = -x231
    x233 = x223 + x68
    x234 = x4 * (x219 - x222 + x233)
    x235 = -x234
    x236 = x219 * x3
    x237 = x225 * x6
    x238 = -x237
    x239 = x236 + x238
    x240 = x219 * x5
    x241 = x214 * x5 + x220
    x242 = x238 + x240
    x243 = x217 + x241 * x3 - x242 * x6
    x244 = x242 * x3
    x245 = x225 * x5
    x246 = x232 + x245
    x247 = x246 * x6
    x248 = x234 + x244 - x247
    x249 = x228 + x79
    x250 = x4 * (x225 - x227 + x249)
    x251 = -x144 * x8 + x77
    x252 = x128 * x87 + x251
    x253 = x158 * x8
    x254 = x153 * x87
    x255 = -x153 * x8
    x256 = x164 * x87 + x255
    x257 = x4 * (x253 - x254 + x256)
    x258 = -x253
    x259 = x254 + x258
    x260 = -x259 * x6
    x261 = x256 * x3 + x260
    x262 = x158 * x87
    x263 = x174 * x8
    x264 = -x263
    x265 = x262 + x264
    x266 = x265 * x3
    x267 = x174 * x87
    x268 = x192 * x8
    x269 = -x268
    x270 = x267 + x269
    x271 = x270 * x6
    x272 = -x271
    x273 = x4 * (x259 - x262 + x263)
    x274 = -x273
    x275 = x259 * x3
    x276 = x265 * x6
    x277 = -x276
    x278 = x275 + x277
    x279 = x259 * x5
    x280 = x256 * x5 + x260
    x281 = x277 + x279
    x282 = x257 + x280 * x3 - x281 * x6
    x283 = x281 * x3
    x284 = x265 * x5
    x285 = x272 + x284
    x286 = x285 * x6
    x287 = x273 + x283 - x286
    x288 = x4 * (x265 - x267 + x268)
    x289 = -x208 * x8
    x290 = x192 * x87 + x289
    x291 = 46.7389915737742 * x83
    x292 = x151 * x153
    x293 = -x10 * x153 + x30
    x294 = x151 * x164 + x293
    x295 = x10 * x158
    x296 = x295 + x49
    x297 = x4 * (-x292 + x294 + x296)
    x298 = x19 - x295
    x299 = x292 + x298
    x300 = -x299 * x6
    x301 = x294 * x3 + x300
    x302 = x151 * x158
    x303 = x10 * x174
    x304 = -x303 + x36
    x305 = x302 + x304
    x306 = x3 * x305
    x307 = x151 * x174
    x308 = x10 * x192
    x309 = -x308 + x55
    x310 = x307 + x309
    x311 = x310 * x6
    x312 = -x311
    x313 = x303 + x68
    x314 = x4 * (x299 - x302 + x313)
    x315 = -x314
    x316 = x299 * x3
    x317 = x305 * x6
    x318 = -x317
    x319 = x316 + x318
    x320 = x299 * x5
    x321 = x294 * x5 + x300
    x322 = x318 + x320
    x323 = x297 + x3 * x321 - x322 * x6
    x324 = x3 * x322
    x325 = x305 * x5
    x326 = x312 + x325
    x327 = x326 * x6
    x328 = x314 + x324 - x327
    x329 = x308 + x79
    x330 = x4 * (x305 - x307 + x329)
    x331 = -x10 * x208 + x77
    x332 = x151 * x192 + x331
    x333 = x219 * x87
    x334 = 2.0 * x103 - x219 * x8
    x335 = x214 * x87 + x334
    x336 = x225 * x8
    x337 = 2.0 * x114
    x338 = x336 - x337
    x339 = x4 * (-x333 + x335 + x338)
    x340 = -x336 + x337
    x341 = x333 + x340
    x342 = x3 * x335 - x341 * x6
    x343 = x3 * x341
    x344 = x225 * x87
    x345 = x230 * x8
    x346 = 2.0 * x123
    x347 = -x345 + x346
    x348 = x344 + x347
    x349 = x348 * x6
    x350 = x343 - x349
    x351 = x345 - x346
    x352 = x4 * (x341 - x344 + x351)
    x353 = 2.0 * x142 - x252 * x8
    x354 = x230 * x87 + x353
    x355 = x259 * x87
    x356 = x167 - x259 * x8
    x357 = x256 * x87 + x356
    x358 = x265 * x8
    x359 = x179 + x358
    x360 = x4 * (-x355 + x357 + x359)
    x361 = x178 - x358
    x362 = x355 + x361
    x363 = x3 * x357 - x362 * x6
    x364 = x3 * x362
    x365 = x265 * x87
    x366 = x270 * x8
    x367 = x187 - x366
    x368 = x365 + x367
    x369 = x368 * x6
    x370 = x364 - x369
    x371 = x204 + x366
    x372 = x4 * (x362 - x365 + x371)
    x373 = x206 - x290 * x8
    x374 = x270 * x87 + x373
    x375 = x305 * x8
    x376 = x299 * x87
    x377 = -x299 * x8
    x378 = x294 * x87 + x377
    x379 = x4 * (x375 - x376 + x378)
    x380 = -x375
    x381 = x376 + x380
    x382 = x3 * x378 - x381 * x6
    x383 = x3 * x381
    x384 = x305 * x87
    x385 = x310 * x8
    x386 = -x385
    x387 = x384 + x386
    x388 = x387 * x6
    x389 = x383 - x388
    x390 = x4 * (x381 - x384 + x385)
    x391 = -x332 * x8
    x392 = x310 * x87 + x391
    x393 = x151 * x299
    x394 = -x10 * x299 + 2.0 * x167
    x395 = x151 * x294 + x394
    x396 = x10 * x305
    x397 = 2.0 * x178
    x398 = x396 - x397
    x399 = x4 * (-x393 + x395 + x398)
    x400 = -x396 + x397
    x401 = x393 + x400
    x402 = x3 * x395 - x401 * x6
    x403 = x3 * x401
    x404 = x151 * x305
    x405 = x10 * x310
    x406 = 2.0 * x187
    x407 = -x405 + x406
    x408 = x404 + x407
    x409 = x408 * x6
    x410 = x403 - x409
    x411 = x405 - x406
    x412 = x4 * (x401 - x404 + x411)
    x413 = -x10 * x332 + 2.0 * x206
    x414 = x151 * x310 + x413
    x415 = -x7 - A[1]
    x416 = x13 * x415
    x417 = x15 * x415
    x418 = x417 + x86
    x419 = x4 * (-x416 + x418 + x91)
    x420 = x418 * x5
    x421 = x416 + x92
    x422 = x421 * x6
    x423 = x420 - x422
    x424 = x423 * x5
    x425 = x421 * x5
    x426 = x20 * x415
    x427 = x108 + x426
    x428 = x427 * x6
    x429 = x425 - x428
    x430 = x429 * x6
    x431 = x419 + x424 - x430
    x432 = x431 * x5
    x433 = x28 * x415 + x99
    x434 = x4 * (-x417 + x433 + x85)
    x435 = -x418 * x6 + x433 * x5
    x436 = -x423 * x6 + x434 + x435 * x5
    x437 = x33 * (-x420 + x422 + x435) - x431 * x6
    x438 = x436 * x5 + x437
    x439 = x4 * (x107 + x421 - x426)
    x440 = x429 * x5
    x441 = x427 * x5
    x442 = x38 * x415
    x443 = x126 + x442
    x444 = x443 * x6
    x445 = x441 - x444
    x446 = x445 * x6
    x447 = x439 + x440 - x446
    x448 = x447 * x6
    x449 = x33 * (x423 - x425 + x428)
    x450 = x448 - x449
    x451 = x432 - x448 + x449
    x452 = -x419
    x453 = x4 * (x125 + x427 - x442)
    x454 = x143 + x415 * x57
    x455 = -x439
    x456 = x415 * x89
    x457 = x218 + x456
    x458 = x457 * x5
    x459 = x415 * x94
    x460 = x224 + x459
    x461 = x460 * x6
    x462 = x458 - x461
    x463 = x462 * x5
    x464 = x100 * x415 + x213
    x465 = -x457 * x6 + x464 * x5
    x466 = x4 * (x216 - x456 + x464)
    x467 = -x462 * x6 + x466
    x468 = x465 * x5 + x467
    x469 = x460 * x5
    x470 = x110 * x415
    x471 = x229 + x470
    x472 = x471 * x6
    x473 = x469 - x472
    x474 = x473 * x6
    x475 = x4 * (x233 + x457 - x459)
    x476 = x474 - x475
    x477 = x463 - x474 + x475
    x478 = x4 * (x249 + x460 - x470)
    x479 = x128 * x415 + x251
    x480 = 60.3397786612521 * x83
    x481 = x153 * x415
    x482 = x258 + x481
    x483 = x482 * x5
    x484 = x158 * x415
    x485 = x264 + x484
    x486 = x485 * x6
    x487 = x483 - x486
    x488 = x487 * x5
    x489 = x164 * x415 + x255
    x490 = -x482 * x6 + x489 * x5
    x491 = x4 * (x253 - x481 + x489)
    x492 = -x487 * x6 + x491
    x493 = x490 * x5 + x492
    x494 = x485 * x5
    x495 = x174 * x415
    x496 = x269 + x495
    x497 = x496 * x6
    x498 = x494 - x497
    x499 = x498 * x6
    x500 = x4 * (x263 + x482 - x484)
    x501 = -x500
    x502 = x499 + x501
    x503 = x488 - x499 + x500
    x504 = x4 * (x268 + x485 - x495)
    x505 = x192 * x415 + x289
    x506 = x214 * x415 + x334
    x507 = x225 * x415
    x508 = x347 + x507
    x509 = x508 * x6
    x510 = x219 * x415
    x511 = x340 + x510
    x512 = -x511 * x6
    x513 = x5 * x511
    x514 = x5 * x506 + x512
    x515 = x4 * (x338 + x506 - x510)
    x516 = -x509 + x513
    x517 = x4 * (x351 - x507 + x511)
    x518 = x230 * x415 + x353
    x519 = x256 * x415 + x356
    x520 = x265 * x415
    x521 = x367 + x520
    x522 = x521 * x6
    x523 = x259 * x415
    x524 = x361 + x523
    x525 = -x524 * x6
    x526 = x5 * x524
    x527 = x5 * x519 + x525
    x528 = x4 * (x359 + x519 - x523)
    x529 = -x522 + x526
    x530 = x4 * (x371 - x520 + x524)
    x531 = x270 * x415 + x373
    x532 = 104.511562358749 * x83
    x533 = x294 * x415 + x377
    x534 = x305 * x415
    x535 = x386 + x534
    x536 = x535 * x6
    x537 = x299 * x415
    x538 = x380 + x537
    x539 = -x538 * x6
    x540 = x5 * x538
    x541 = x5 * x533 + x539
    x542 = x4 * (x375 + x533 - x537)
    x543 = -x536 + x540
    x544 = x4 * (x385 - x534 + x538)
    x545 = x310 * x415 + x391
    x546 = x348 * x8
    x547 = x341 * x415
    x548 = 3.0 * x234
    x549 = 3.0 * x217 + x335 * x415 - x341 * x8
    x550 = x4 * (x546 - x547 - x548 + x549)
    x551 = -x546 + x547 + x548
    x552 = 3.0 * x250 + x348 * x415 - x354 * x8
    x553 = x368 * x8
    x554 = x362 * x415
    x555 = 2.0 * x273
    x556 = 2.0 * x257 + x357 * x415 - x362 * x8
    x557 = x4 * (x553 - x554 - x555 + x556)
    x558 = -x553 + x554 + x555
    x559 = 2.0 * x288 + x368 * x415 - x374 * x8
    x560 = x387 * x8
    x561 = x381 * x415
    x562 = x297 + x378 * x415 - x381 * x8
    x563 = x4 * (x315 + x560 - x561 + x562)
    x564 = x314 - x560 + x561
    x565 = x330 + x387 * x415 - x392 * x8
    x566 = x408 * x8
    x567 = x401 * x415
    x568 = x395 * x415 - x401 * x8
    x569 = x4 * (x566 - x567 + x568)
    x570 = -x566 + x567
    x571 = x408 * x415 - x414 * x8
    x572 = -x9 - A[2]
    x573 = x13 * x572
    x574 = x15 * x572
    x575 = x150 + x574
    x576 = x4 * (x155 - x573 + x575)
    x577 = x5 * x575
    x578 = x156 + x573
    x579 = x578 * x6
    x580 = x577 - x579
    x581 = x5 * x580
    x582 = x5 * x578
    x583 = x20 * x572
    x584 = x172 + x583
    x585 = x584 * x6
    x586 = x582 - x585
    x587 = x586 * x6
    x588 = x576 + x581 - x587
    x589 = x5 * x588
    x590 = x163 + x28 * x572
    x591 = x4 * (x149 - x574 + x590)
    x592 = x5 * x590 - x575 * x6
    x593 = x5 * x592 - x580 * x6 + x591
    x594 = x33 * (-x577 + x579 + x592) - x588 * x6
    x595 = x5 * x593 + x594
    x596 = x4 * (x171 + x578 - x583)
    x597 = x5 * x586
    x598 = x5 * x584
    x599 = x38 * x572
    x600 = x190 + x599
    x601 = x6 * x600
    x602 = x598 - x601
    x603 = x6 * x602
    x604 = x596 + x597 - x603
    x605 = x6 * x604
    x606 = x33 * (x580 - x582 + x585)
    x607 = x605 - x606
    x608 = x589 - x605 + x606
    x609 = -x576
    x610 = x4 * (x189 + x584 - x599)
    x611 = x207 + x57 * x572
    x612 = -x596
    x613 = x575 * x87
    x614 = x578 * x8
    x615 = -x614
    x616 = x613 + x615
    x617 = x5 * x616
    x618 = x578 * x87
    x619 = x584 * x8
    x620 = -x619
    x621 = x618 + x620
    x622 = x6 * x621
    x623 = x617 - x622
    x624 = x5 * x623
    x625 = -x575 * x8
    x626 = x590 * x87 + x625
    x627 = x5 * x626 - x6 * x616
    x628 = x4 * (-x613 + x614 + x626)
    x629 = -x6 * x623 + x628
    x630 = x5 * x627 + x629
    x631 = x5 * x621
    x632 = x584 * x87
    x633 = x600 * x8
    x634 = -x633
    x635 = x632 + x634
    x636 = x6 * x635
    x637 = x631 - x636
    x638 = x6 * x637
    x639 = x4 * (x616 - x618 + x619)
    x640 = -x639
    x641 = x638 + x640
    x642 = x624 - x638 + x639
    x643 = x4 * (x621 - x632 + x633)
    x644 = -x611 * x8
    x645 = x600 * x87 + x644
    x646 = x153 * x572
    x647 = x298 + x646
    x648 = x5 * x647
    x649 = x158 * x572
    x650 = x304 + x649
    x651 = x6 * x650
    x652 = x648 - x651
    x653 = x5 * x652
    x654 = x164 * x572 + x293
    x655 = x5 * x654 - x6 * x647
    x656 = x4 * (x296 - x646 + x654)
    x657 = -x6 * x652 + x656
    x658 = x5 * x655 + x657
    x659 = x5 * x650
    x660 = x174 * x572
    x661 = x309 + x660
    x662 = x6 * x661
    x663 = x659 - x662
    x664 = x6 * x663
    x665 = x4 * (x313 + x647 - x649)
    x666 = -x665
    x667 = x664 + x666
    x668 = x653 - x664 + x665
    x669 = x4 * (x329 + x650 - x660)
    x670 = x192 * x572 + x331
    x671 = x591 - x616 * x8
    x672 = x626 * x87 + x671
    x673 = x621 * x87
    x674 = x635 * x8
    x675 = x596 - x674
    x676 = x673 + x675
    x677 = x6 * x676
    x678 = x616 * x87
    x679 = x621 * x8
    x680 = x576 - x679
    x681 = x678 + x680
    x682 = -x6 * x681
    x683 = x5 * x681
    x684 = x5 * x672 + x682
    x685 = x609 + x679
    x686 = x4 * (x672 - x678 + x685)
    x687 = -x677 + x683
    x688 = x612 + x674
    x689 = x4 * (-x673 + x681 + x688)
    x690 = x610 - x645 * x8
    x691 = x635 * x87 + x690
    x692 = -x647 * x8
    x693 = x654 * x87 + x692
    x694 = x650 * x87
    x695 = x661 * x8
    x696 = -x695
    x697 = x694 + x696
    x698 = x6 * x697
    x699 = x647 * x87
    x700 = x650 * x8
    x701 = -x700
    x702 = x699 + x701
    x703 = -x6 * x702
    x704 = x5 * x702
    x705 = x5 * x693 + x703
    x706 = x4 * (x693 - x699 + x700)
    x707 = -x698 + x704
    x708 = x4 * (-x694 + x695 + x702)
    x709 = -x670 * x8
    x710 = x661 * x87 + x709
    x711 = x294 * x572 + x394
    x712 = x305 * x572
    x713 = x407 + x712
    x714 = x6 * x713
    x715 = x299 * x572
    x716 = x400 + x715
    x717 = -x6 * x716
    x718 = x5 * x716
    x719 = x5 * x711 + x717
    x720 = x4 * (x398 + x711 - x715)
    x721 = -x714 + x718
    x722 = x4 * (x411 - x712 + x716)
    x723 = x310 * x572 + x413
    x724 = x681 * x87
    x725 = 2.0 * x628 - x681 * x8
    x726 = x672 * x87 + x725
    x727 = x676 * x8
    x728 = 2.0 * x639
    x729 = x727 - x728
    x730 = x4 * (-x724 + x726 + x729)
    x731 = -x727 + x728
    x732 = x724 + x731
    x733 = 2.0 * x643 - x691 * x8
    x734 = x676 * x87 + x733
    x735 = x702 * x87
    x736 = x656 - x702 * x8
    x737 = x693 * x87 + x736
    x738 = x697 * x8
    x739 = x666 + x738
    x740 = x4 * (-x735 + x737 + x739)
    x741 = x665 - x738
    x742 = x735 + x741
    x743 = x669 - x710 * x8
    x744 = x697 * x87 + x743
    x745 = x713 * x8
    x746 = x716 * x87
    x747 = -x716 * x8
    x748 = x711 * x87 + x747
    x749 = x4 * (x745 - x746 + x748)
    x750 = -x745
    x751 = x746 + x750
    x752 = -x723 * x8
    x753 = x713 * x87 + x752
    x754 = x10 * x408
    x755 = x401 * x572
    x756 = 3.0 * x314
    x757 = -x10 * x401 + 3.0 * x297 + x395 * x572
    x758 = x4 * (x754 - x755 - x756 + x757)
    x759 = -x754 + x755 + x756
    x760 = -x10 * x414 + 3.0 * x330 + x408 * x572
    x761 = x421 * x8
    x762 = x415 * x418
    x763 = x30 + x415 * x433 - x418 * x8
    x764 = x4 * (x49 + x761 - x762 + x763)
    x765 = x19 - x761 + x762
    x766 = x5 * x763 - x6 * x765
    x767 = x5 * x765
    x768 = x415 * x421
    x769 = x427 * x8
    x770 = x36 + x768 - x769
    x771 = x6 * x770
    x772 = x767 - x771
    x773 = x5 * x766 - x6 * x772 + x764
    x774 = x4 * (x68 + x765 - x768 + x769)
    x775 = x5 * x772
    x776 = x5 * x770
    x777 = x415 * x427
    x778 = x443 * x8
    x779 = x55 + x777 - x778
    x780 = x6 * x779
    x781 = x776 - x780
    x782 = x6 * x781
    x783 = x774 + x775 - x782
    x784 = x4 * (x770 - x777 + x778 + x79)
    x785 = x415 * x443 - x454 * x8 + x77
    x786 = -x774
    x787 = x460 * x8
    x788 = x415 * x457
    x789 = x103 + x415 * x464 + x434 - x457 * x8
    x790 = x4 * (x115 + x452 + x787 - x788 + x789)
    x791 = x114 + x419 - x787 + x788
    x792 = x5 * x789 - x6 * x791
    x793 = x5 * x791
    x794 = x415 * x460
    x795 = x471 * x8
    x796 = x123 + x439 + x794 - x795
    x797 = x6 * x796
    x798 = x793 - x797
    x799 = x4 * (x140 + x455 + x791 - x794 + x795)
    x800 = x142 + x415 * x471 + x453 - x479 * x8
    x801 = x485 * x8
    x802 = x415 * x482
    x803 = x167 + x415 * x489 - x482 * x8
    x804 = x4 * (x179 + x801 - x802 + x803)
    x805 = x178 - x801 + x802
    x806 = x5 * x803 - x6 * x805
    x807 = x5 * x805
    x808 = x415 * x485
    x809 = x496 * x8
    x810 = x187 + x808 - x809
    x811 = x6 * x810
    x812 = x807 - x811
    x813 = x4 * (x204 + x805 - x808 + x809)
    x814 = x206 + x415 * x496 - x505 * x8
    x815 = x508 * x8
    x816 = x415 * x511
    x817 = 2.0 * x475
    x818 = -x817
    x819 = 2.0 * x466
    x820 = x217 + x415 * x506 - x511 * x8 + x819
    x821 = x4 * (x235 + x815 - x816 + x818 + x820)
    x822 = x234 - x815 + x816 + x817
    x823 = 2.0 * x478
    x824 = x250 + x415 * x508 - x518 * x8 + x823
    x825 = x521 * x8
    x826 = x415 * x524
    x827 = x257 + x415 * x519 + x491 - x524 * x8
    x828 = x4 * (x274 + x501 + x825 - x826 + x827)
    x829 = x273 + x500 - x825 + x826
    x830 = x288 + x415 * x521 + x504 - x531 * x8
    x831 = x535 * x8
    x832 = x415 * x538
    x833 = x297 + x415 * x533 - x538 * x8
    x834 = x4 * (x315 + x831 - x832 + x833)
    x835 = x314 - x831 + x832
    x836 = x330 + x415 * x535 - x545 * x8
    x837 = x339 + x415 * x549 + 3.0 * x515 - x551 * x8
    x838 = x352 + x415 * x551 + 3.0 * x517 - x552 * x8
    x839 = 2.0 * x528
    x840 = x360 + x415 * x556 - x558 * x8 + x839
    x841 = 2.0 * x530
    x842 = x372 + x415 * x558 - x559 * x8 + x841
    x843 = x379 + x415 * x562 + x542 - x564 * x8
    x844 = x390 + x415 * x564 + x544 - x565 * x8
    x845 = x399 + x415 * x568 - x570 * x8
    x846 = x412 + x415 * x570 - x571 * x8
    x847 = x415 * x575
    x848 = x415 * x590 + x625
    x849 = x4 * (x614 - x847 + x848)
    x850 = x615 + x847
    x851 = x5 * x848 - x6 * x850
    x852 = x5 * x850
    x853 = x415 * x578
    x854 = x620 + x853
    x855 = x6 * x854
    x856 = x852 - x855
    x857 = x5 * x851 - x6 * x856 + x849
    x858 = x4 * (x619 + x850 - x853)
    x859 = x5 * x856
    x860 = x5 * x854
    x861 = x415 * x584
    x862 = x634 + x861
    x863 = x6 * x862
    x864 = x860 - x863
    x865 = x6 * x864
    x866 = x858 + x859 - x865
    x867 = x4 * (x633 + x854 - x861)
    x868 = x415 * x600 + x644
    x869 = -x858
    x870 = x415 * x616
    x871 = x415 * x626 + x671
    x872 = x4 * (x685 - x870 + x871)
    x873 = x680 + x870
    x874 = x5 * x871 - x6 * x873
    x875 = x5 * x873
    x876 = x415 * x621
    x877 = x675 + x876
    x878 = x6 * x877
    x879 = x875 - x878
    x880 = x4 * (x688 + x873 - x876)
    x881 = x415 * x635 + x690
    x882 = x415 * x647
    x883 = x415 * x654 + x692
    x884 = x4 * (x700 - x882 + x883)
    x885 = x701 + x882
    x886 = x5 * x883 - x6 * x885
    x887 = x5 * x885
    x888 = x415 * x650
    x889 = x696 + x888
    x890 = x6 * x889
    x891 = x887 - x890
    x892 = x4 * (x695 + x885 - x888)
    x893 = x415 * x661 + x709
    x894 = x415 * x681
    x895 = x415 * x672 + x725
    x896 = x4 * (x729 - x894 + x895)
    x897 = x731 + x894
    x898 = x415 * x676 + x733
    x899 = x415 * x702
    x900 = x415 * x693 + x736
    x901 = x4 * (x739 - x899 + x900)
    x902 = x741 + x899
    x903 = x415 * x697 + x743
    x904 = x415 * x716
    x905 = x415 * x711 + x747
    x906 = x4 * (x745 - x904 + x905)
    x907 = x750 + x904
    x908 = x415 * x713 + x752
    x909 = x415 * x726 + 3.0 * x686 - x732 * x8
    x910 = x415 * x732 + 3.0 * x689 - x734 * x8
    x911 = x415 * x737 + 2.0 * x706 - x742 * x8
    x912 = x415 * x742 + 2.0 * x708 - x744 * x8
    x913 = x415 * x748 + x720 - x751 * x8
    x914 = x415 * x751 + x722 - x753 * x8
    x915 = x415 * x757 - x759 * x8
    x916 = x415 * x759 - x760 * x8
    x917 = x10 * x578
    x918 = x572 * x575
    x919 = -x10 * x575 + x30 + x572 * x590
    x920 = x4 * (x49 + x917 - x918 + x919)
    x921 = x19 - x917 + x918
    x922 = x5 * x919 - x6 * x921
    x923 = x5 * x921
    x924 = x572 * x578
    x925 = x10 * x584
    x926 = x36 + x924 - x925
    x927 = x6 * x926
    x928 = x923 - x927
    x929 = x5 * x922 - x6 * x928 + x920
    x930 = x4 * (x68 + x921 - x924 + x925)
    x931 = x5 * x928
    x932 = x5 * x926
    x933 = x572 * x584
    x934 = x10 * x600
    x935 = x55 + x933 - x934
    x936 = x6 * x935
    x937 = x932 - x936
    x938 = x6 * x937
    x939 = x930 + x931 - x938
    x940 = x4 * (x79 + x926 - x933 + x934)
    x941 = -x10 * x611 + x572 * x600 + x77
    x942 = -x930
    x943 = x8 * x926
    x944 = x87 * x921
    x945 = -x8 * x921
    x946 = x87 * x919 + x945
    x947 = x4 * (x943 - x944 + x946)
    x948 = -x943
    x949 = x944 + x948
    x950 = x5 * x946 - x6 * x949
    x951 = x5 * x949
    x952 = x87 * x926
    x953 = x8 * x935
    x954 = -x953
    x955 = x952 + x954
    x956 = x6 * x955
    x957 = x951 - x956
    x958 = x4 * (x949 - x952 + x953)
    x959 = -x8 * x941
    x960 = x87 * x935 + x959
    x961 = x10 * x650
    x962 = x572 * x647
    x963 = -x10 * x647 + x167 + x572 * x654 + x591
    x964 = x4 * (x179 + x609 + x961 - x962 + x963)
    x965 = x178 + x576 - x961 + x962
    x966 = x5 * x963 - x6 * x965
    x967 = x5 * x965
    x968 = x572 * x650
    x969 = x10 * x661
    x970 = x187 + x596 + x968 - x969
    x971 = x6 * x970
    x972 = x967 - x971
    x973 = x4 * (x204 + x612 + x965 - x968 + x969)
    x974 = -x10 * x670 + x206 + x572 * x661 + x610
    x975 = x87 * x949
    x976 = -x8 * x949 + x920
    x977 = x87 * x946 + x976
    x978 = x8 * x955
    x979 = x942 + x978
    x980 = x4 * (-x975 + x977 + x979)
    x981 = x930 - x978
    x982 = x975 + x981
    x983 = -x8 * x960 + x940
    x984 = x87 * x955 + x983
    x985 = x8 * x970
    x986 = x87 * x965
    x987 = -x8 * x965
    x988 = x87 * x963 + x987
    x989 = x4 * (x985 - x986 + x988)
    x990 = -x985
    x991 = x986 + x990
    x992 = -x8 * x974
    x993 = x87 * x970 + x992
    x994 = x10 * x713
    x995 = x572 * x716
    x996 = 2.0 * x665
    x997 = -x996
    x998 = 2.0 * x656
    x999 = -x10 * x716 + x297 + x572 * x711 + x998
    x1000 = x4 * (x315 + x994 - x995 + x997 + x999)
    x1001 = x314 - x994 + x995 + x996
    x1002 = 2.0 * x669
    x1003 = -x10 * x723 + x1002 + x330 + x572 * x713
    x1004 = -x8 * x982 + 2.0 * x947
    x1005 = x1004 + x87 * x977
    x1006 = -x8 * x984 + 2.0 * x958
    x1007 = x1006 + x87 * x982
    x1008 = -x8 * x991 + x964
    x1009 = x1008 + x87 * x988
    x1010 = -x8 * x993 + x973
    x1011 = x1010 + x87 * x991
    x1012 = -x1001 * x8
    x1013 = x1012 + x87 * x999
    x1014 = -x1003 * x8
    x1015 = x1001 * x87 + x1014
    x1016 = -x10 * x759 + x399 + x572 * x757 + 3.0 * x720
    x1017 = -x10 * x760 + x412 + x572 * x759 + 3.0 * x722
    x1018 = x770 * x8
    x1019 = x415 * x765
    x1020 = 2.0 * x419
    x1021 = x415 * x763 + 2.0 * x434 - x765 * x8
    x1022 = -x1018 + x1019 + x1020
    x1023 = x1021 * x5 - x1022 * x6
    x1024 = x1022 * x5
    x1025 = x415 * x770
    x1026 = x779 * x8
    x1027 = 2.0 * x439
    x1028 = x1025 - x1026 + x1027
    x1029 = x1028 * x6
    x1030 = x1024 - x1029
    x1031 = x796 * x8
    x1032 = x415 * x791
    x1033 = x415 * x789 + x764 - x791 * x8 + x819
    x1034 = -x1031 + x1032 + x774 + x817
    x1035 = x8 * x810
    x1036 = x415 * x805
    x1037 = 2.0 * x500
    x1038 = x415 * x803 + 2.0 * x491 - x8 * x805
    x1039 = -x1035 + x1036 + x1037
    x1040 = x8 * x854
    x1041 = x415 * x850
    x1042 = x415 * x848 + x591 - x8 * x850
    x1043 = -x1040 + x1041 + x576
    x1044 = x1042 * x5 - x1043 * x6
    x1045 = x1043 * x5
    x1046 = x415 * x854
    x1047 = x8 * x862
    x1048 = x1046 - x1047 + x596
    x1049 = x1048 * x6
    x1050 = x1045 - x1049
    x1051 = x8 * x877
    x1052 = x415 * x873
    x1053 = x415 * x871 + x628 - x8 * x873 + x849
    x1054 = -x1051 + x1052 + x639 + x858
    x1055 = x8 * x889
    x1056 = x415 * x885
    x1057 = x415 * x883 + x656 - x8 * x885
    x1058 = -x1055 + x1056 + x665
    x1059 = x415 * x921
    x1060 = x415 * x919 + x945
    x1061 = x1059 + x948
    x1062 = x1060 * x5 - x1061 * x6
    x1063 = x1061 * x5
    x1064 = x415 * x926
    x1065 = x1064 + x954
    x1066 = x1065 * x6
    x1067 = x1063 - x1066
    x1068 = x415 * x949
    x1069 = x415 * x946 + x976
    x1070 = x1068 + x981
    x1071 = x415 * x965
    x1072 = x415 * x963 + x987
    x1073 = x1071 + x990
    x1074 = x10 * x926
    x1075 = x572 * x921
    x1076 = 2.0 * x576
    x1077 = -x10 * x921 + x572 * x919 + 2.0 * x591
    x1078 = x4 * (x1074 - x1075 - x1076 + x1077)
    x1079 = -x1074 + x1075 + x1076
    x1080 = x1077 * x5 - x1079 * x6
    x1081 = x1079 * x5
    x1082 = x572 * x926
    x1083 = x10 * x935
    x1084 = 2.0 * x596
    x1085 = x1082 - x1083 + x1084
    x1086 = x1085 * x6
    x1087 = x1081 - x1086
    x1088 = x4 * (x1079 - x1082 + x1083 - x1084)
    x1089 = -x10 * x941 + x572 * x935 + 2.0 * x610
    x1090 = x1085 * x8
    x1091 = x1079 * x87
    x1092 = x1077 * x87 - x1079 * x8
    x1093 = x4 * (x1090 - x1091 + x1092)
    x1094 = -x1090 + x1091
    x1095 = x1085 * x87 - x1089 * x8
    x1096 = x10 * x970
    x1097 = x572 * x965
    x1098 = -x10 * x965 + x572 * x963 + x920 + x998
    x1099 = x4 * (x1096 - x1097 + x1098 + x942 + x997)
    x1100 = -x1096 + x1097 + x930 + x996
    x1101 = -x10 * x974 + x1002 + x572 * x970 + x940
    x1102 = x1078 + x1092 * x87 - x1094 * x8
    x1103 = x1088 + x1094 * x87 - x1095 * x8
    x1104 = x1098 * x87 - x1100 * x8
    x1105 = x1100 * x87 - x1101 * x8
    x1106 = -x10 * x1001 + x572 * x999 + 2.0 * x720 + 2.0 * x964
    x1107 = -x10 * x1003 + x1001 * x572 + 2.0 * x722 + 2.0 * x973

    # 100 item(s)
    result[0, 0] = numpy.sum(
        x84
        * (
            x3 * (x3 * x52 + x4 * (-x27 + x35 + x46) - x6 * x72 + x75 * (x46 - x73 + x74))
            + x33 * (x52 - x53 + x67 - x71)
            - x6
            * (
                x3 * x72
                + x4 * (x48 - x54 + x76)
                - x6
                * (
                    x3 * x66
                    - x6
                    * (
                        x33 * (x41 + x58 - x59)
                        + x5 * x62
                        - x6 * (x5 * x60 - x6 * (x5 * x57 - x6 * x78) + x77)
                    )
                    + x75 * (x43 - x56 + x61 + x79)
                )
                + x75 * (x76 - x80 + x81)
            )
            + x75
            * (
                x3 * x74
                - x3 * x81
                - x33 * (-x23 * x3 + x25 + x69 + x82)
                + x33 * (x29 * x3 + x31 + x50 - x82)
                + x51
                - x6 * x81
                + x6 * (x65 + x80)
                - x70
            )
        )
    )
    result[0, 1] = numpy.sum(
        x148
        * (
            x3
            * (
                x120 * x3
                - x137 * x6
                + x33 * (x116 - x138 + x139)
                + x4 * (x105 + x116 - x98)
            )
            + x33 * (x120 - x121 + x134 - x136)
            + x33
            * (
                x119
                - x135
                + x139 * x3
                - x146 * x3
                - x146 * x6
                - x4 * (x111 + x147 - x3 * x94 + x96)
                + x4 * (x100 * x3 + x101 - x147 + x95)
                + x6 * (x132 + x145)
            )
            - x6
            * (
                x137 * x3
                + x33 * (x141 - x145 + x146)
                + x4 * (x118 - x122 + x141)
                - x6
                * (
                    x133 * x3
                    + x33 * (x112 - x124 + x129)
                    - x6 * (x130 * x5 + x142 - x6 * (x128 * x5 - x144 * x6))
                )
            )
        )
    )
    result[0, 2] = numpy.sum(
        x148
        * (
            x3
            * (
                x184 * x3
                - x201 * x6
                + x33 * (x180 - x202 + x203)
                + x4 * (-x162 + x169 + x180)
            )
            + x33 * (x184 - x185 + x198 - x200)
            + x33
            * (
                x183
                - x199
                + x203 * x3
                - x210 * x3
                - x210 * x6
                + x4 * (x159 + x164 * x3 + x165 - x211)
                - x4 * (-x158 * x3 + x160 + x175 + x211)
                + x6 * (x196 + x209)
            )
            - x6
            * (
                x201 * x3
                + x33 * (x205 - x209 + x210)
                + x4 * (x182 - x186 + x205)
                - x6
                * (
                    x197 * x3
                    + x33 * (x176 - x188 + x193)
                    - x6 * (x194 * x5 + x206 - x6 * (x192 * x5 - x208 * x6))
                )
            )
        )
    )
    result[0, 3] = numpy.sum(
        x148
        * (
            x3
            * (
                x243 * x3
                - x248 * x6
                + x4 * (x221 - x236 + x237)
                + x4 * (x237 - x240 + x241)
            )
            + x33 * (x235 + x243 - x244 + x247)
            + x4 * (x217 + x221 * x3 + x235 - x239 * x3 - x239 * x6 + x6 * (x226 + x232))
            - x6
            * (
                x248 * x3
                + x4 * (-x226 + x231 + x239)
                + x4 * (x231 + x242 - x245)
                - x6 * (x246 * x3 + x250 - x6 * (x230 * x5 - x252 * x6))
            )
        )
    )
    result[0, 4] = numpy.sum(
        x291
        * (
            x3
            * (
                x282 * x3
                - x287 * x6
                + x4 * (x261 - x275 + x276)
                + x4 * (x276 - x279 + x280)
            )
            + x33 * (x274 + x282 - x283 + x286)
            + x4 * (x257 + x261 * x3 + x274 - x278 * x3 - x278 * x6 + x6 * (x266 + x272))
            - x6
            * (
                x287 * x3
                + x4 * (-x266 + x271 + x278)
                + x4 * (x271 + x281 - x284)
                - x6 * (x285 * x3 + x288 - x6 * (x270 * x5 - x290 * x6))
            )
        )
    )
    result[0, 5] = numpy.sum(
        x148
        * (
            x3
            * (
                x3 * x323
                - x328 * x6
                + x4 * (x301 - x316 + x317)
                + x4 * (x317 - x320 + x321)
            )
            + x33 * (x315 + x323 - x324 + x327)
            + x4 * (x297 + x3 * x301 - x3 * x319 + x315 - x319 * x6 + x6 * (x306 + x312))
            - x6
            * (
                x3 * x328
                + x4 * (-x306 + x311 + x319)
                + x4 * (x311 + x322 - x325)
                - x6 * (x3 * x326 + x330 - x6 * (x310 * x5 - x332 * x6))
            )
        )
    )
    result[0, 6] = numpy.sum(
        x84
        * (
            x3 * (x3 * x342 + x339 - x350 * x6)
            + x33 * (x342 - x343 + x349)
            - x6 * (x3 * x350 + x352 - x6 * (x3 * x348 - x354 * x6))
        )
    )
    result[0, 7] = numpy.sum(
        x148
        * (
            x3 * (x3 * x363 + x360 - x370 * x6)
            + x33 * (x363 - x364 + x369)
            - x6 * (x3 * x370 + x372 - x6 * (x3 * x368 - x374 * x6))
        )
    )
    result[0, 8] = numpy.sum(
        x148
        * (
            x3 * (x3 * x382 + x379 - x389 * x6)
            + x33 * (x382 - x383 + x388)
            - x6 * (x3 * x389 + x390 - x6 * (x3 * x387 - x392 * x6))
        )
    )
    result[0, 9] = numpy.sum(
        x84
        * (
            x3 * (x3 * x402 + x399 - x410 * x6)
            + x33 * (x402 - x403 + x409)
            - x6 * (x3 * x410 + x412 - x6 * (x3 * x408 - x414 * x6))
        )
    )
    result[1, 0] = numpy.sum(
        x148
        * (
            x3 * (x3 * x438 - x451 * x6 + x75 * (-x424 + x430 + x436 + x452))
            + x4 * (-x432 + x438 + x450)
            - x6
            * (
                x3 * x451
                - x6
                * (
                    x33 * (x429 - x441 + x444)
                    + x447 * x5
                    - x6 * (x445 * x5 + x453 - x6 * (x443 * x5 - x454 * x6))
                )
                + x75 * (x431 - x440 + x446 + x455)
            )
            + x75 * (-x3 * x431 + x3 * x436 + x437 + x450)
        )
    )
    result[1, 1] = numpy.sum(
        x480
        * (
            x3 * (x3 * x468 + x33 * (-x458 + x461 + x465) - x477 * x6)
            + x33 * (-x3 * x462 + x3 * x465 + x467 + x476)
            + x4 * (-x463 + x468 + x476)
            - x6
            * (
                x3 * x477
                + x33 * (x462 - x469 + x472)
                - x6 * (x473 * x5 + x478 - x6 * (x471 * x5 - x479 * x6))
            )
        )
    )
    result[1, 2] = numpy.sum(
        x480
        * (
            x3 * (x3 * x493 + x33 * (-x483 + x486 + x490) - x503 * x6)
            + x33 * (-x3 * x487 + x3 * x490 + x492 + x502)
            + x4 * (-x488 + x493 + x502)
            - x6
            * (
                x3 * x503
                + x33 * (x487 - x494 + x497)
                - x6 * (x498 * x5 + x504 - x6 * (x496 * x5 - x505 * x6))
            )
        )
    )
    result[1, 3] = numpy.sum(
        x480
        * (
            x3 * (x3 * x514 + x515 - x516 * x6)
            + x4 * (x509 - x513 + x514)
            + x4 * (x3 * x506 - x3 * x511 + x509 + x512)
            - x6 * (x3 * x516 + x517 - x6 * (x5 * x508 - x518 * x6))
        )
    )
    result[1, 4] = numpy.sum(
        x532
        * (
            x3 * (x3 * x527 + x528 - x529 * x6)
            + x4 * (x522 - x526 + x527)
            + x4 * (x3 * x519 - x3 * x524 + x522 + x525)
            - x6 * (x3 * x529 + x530 - x6 * (x5 * x521 - x531 * x6))
        )
    )
    result[1, 5] = numpy.sum(
        x480
        * (
            x3 * (x3 * x541 + x542 - x543 * x6)
            + x4 * (x536 - x540 + x541)
            + x4 * (x3 * x533 - x3 * x538 + x536 + x539)
            - x6 * (x3 * x543 + x544 - x6 * (x5 * x535 - x545 * x6))
        )
    )
    result[1, 6] = numpy.sum(
        x148 * (x3 * (x3 * x549 - x551 * x6) + x550 - x6 * (x3 * x551 - x552 * x6))
    )
    result[1, 7] = numpy.sum(
        x480 * (x3 * (x3 * x556 - x558 * x6) + x557 - x6 * (x3 * x558 - x559 * x6))
    )
    result[1, 8] = numpy.sum(
        x480 * (x3 * (x3 * x562 - x564 * x6) + x563 - x6 * (x3 * x564 - x565 * x6))
    )
    result[1, 9] = numpy.sum(
        x148 * (x3 * (x3 * x568 - x570 * x6) + x569 - x6 * (x3 * x570 - x571 * x6))
    )
    result[2, 0] = numpy.sum(
        x148
        * (
            x3 * (x3 * x595 - x6 * x608 + x75 * (-x581 + x587 + x593 + x609))
            + x4 * (-x589 + x595 + x607)
            - x6
            * (
                x3 * x608
                - x6
                * (
                    x33 * (x586 - x598 + x601)
                    + x5 * x604
                    - x6 * (x5 * x602 - x6 * (x5 * x600 - x6 * x611) + x610)
                )
                + x75 * (x588 - x597 + x603 + x612)
            )
            + x75 * (-x3 * x588 + x3 * x593 + x594 + x607)
        )
    )
    result[2, 1] = numpy.sum(
        x480
        * (
            x3 * (x3 * x630 + x33 * (-x617 + x622 + x627) - x6 * x642)
            + x33 * (-x3 * x623 + x3 * x627 + x629 + x641)
            + x4 * (-x624 + x630 + x641)
            - x6
            * (
                x3 * x642
                + x33 * (x623 - x631 + x636)
                - x6 * (x5 * x637 - x6 * (x5 * x635 - x6 * x645) + x643)
            )
        )
    )
    result[2, 2] = numpy.sum(
        x480
        * (
            x3 * (x3 * x658 + x33 * (-x648 + x651 + x655) - x6 * x668)
            + x33 * (-x3 * x652 + x3 * x655 + x657 + x667)
            + x4 * (-x653 + x658 + x667)
            - x6
            * (
                x3 * x668
                + x33 * (x652 - x659 + x662)
                - x6 * (x5 * x663 - x6 * (x5 * x661 - x6 * x670) + x669)
            )
        )
    )
    result[2, 3] = numpy.sum(
        x480
        * (
            x3 * (x3 * x684 - x6 * x687 + x686)
            + x4 * (x677 - x683 + x684)
            + x4 * (x3 * x672 - x3 * x681 + x677 + x682)
            - x6 * (x3 * x687 - x6 * (x5 * x676 - x6 * x691) + x689)
        )
    )
    result[2, 4] = numpy.sum(
        x532
        * (
            x3 * (x3 * x705 - x6 * x707 + x706)
            + x4 * (x698 - x704 + x705)
            + x4 * (x3 * x693 - x3 * x702 + x698 + x703)
            - x6 * (x3 * x707 - x6 * (x5 * x697 - x6 * x710) + x708)
        )
    )
    result[2, 5] = numpy.sum(
        x480
        * (
            x3 * (x3 * x719 - x6 * x721 + x720)
            + x4 * (x714 - x718 + x719)
            + x4 * (x3 * x711 - x3 * x716 + x714 + x717)
            - x6 * (x3 * x721 - x6 * (x5 * x713 - x6 * x723) + x722)
        )
    )
    result[2, 6] = numpy.sum(
        x148 * (x3 * (x3 * x726 - x6 * x732) - x6 * (x3 * x732 - x6 * x734) + x730)
    )
    result[2, 7] = numpy.sum(
        x480 * (x3 * (x3 * x737 - x6 * x742) - x6 * (x3 * x742 - x6 * x744) + x740)
    )
    result[2, 8] = numpy.sum(
        x480 * (x3 * (x3 * x748 - x6 * x751) - x6 * (x3 * x751 - x6 * x753) + x749)
    )
    result[2, 9] = numpy.sum(
        x148 * (x3 * (x3 * x757 - x6 * x759) - x6 * (x3 * x759 - x6 * x760) + x758)
    )
    result[3, 0] = numpy.sum(
        x148
        * (
            x3 * (x33 * (x766 - x767 + x771) + x5 * x773 - x6 * x783)
            - x6
            * (
                x33 * (x772 - x776 + x780)
                + x5 * x783
                - x6 * (x5 * x781 - x6 * (x5 * x779 - x6 * x785) + x784)
            )
            + x75 * (x773 - x775 + x782 + x786)
        )
    )
    result[3, 1] = numpy.sum(
        x480
        * (
            x3 * (x5 * x792 - x6 * x798 + x790)
            + x33 * (x792 - x793 + x797)
            - x6 * (x5 * x798 - x6 * (x5 * x796 - x6 * x800) + x799)
        )
    )
    result[3, 2] = numpy.sum(
        x480
        * (
            x3 * (x5 * x806 - x6 * x812 + x804)
            + x33 * (x806 - x807 + x811)
            - x6 * (x5 * x812 - x6 * (x5 * x810 - x6 * x814) + x813)
        )
    )
    result[3, 3] = numpy.sum(
        x480 * (x3 * (x5 * x820 - x6 * x822) - x6 * (x5 * x822 - x6 * x824) + x821)
    )
    result[3, 4] = numpy.sum(
        x532 * (x3 * (x5 * x827 - x6 * x829) - x6 * (x5 * x829 - x6 * x830) + x828)
    )
    result[3, 5] = numpy.sum(
        x480 * (x3 * (x5 * x833 - x6 * x835) - x6 * (x5 * x835 - x6 * x836) + x834)
    )
    result[3, 6] = numpy.sum(x148 * (x3 * x837 - x6 * x838))
    result[3, 7] = numpy.sum(x480 * (x3 * x840 - x6 * x842))
    result[3, 8] = numpy.sum(x480 * (x3 * x843 - x6 * x844))
    result[3, 9] = numpy.sum(x148 * (x3 * x845 - x6 * x846))
    result[4, 0] = numpy.sum(
        x291
        * (
            x3 * (x33 * (x851 - x852 + x855) + x5 * x857 - x6 * x866)
            - x6
            * (
                x33 * (x856 - x860 + x863)
                + x5 * x866
                - x6 * (x5 * x864 - x6 * (x5 * x862 - x6 * x868) + x867)
            )
            + x75 * (x857 - x859 + x865 + x869)
        )
    )
    result[4, 1] = numpy.sum(
        x532
        * (
            x3 * (x5 * x874 - x6 * x879 + x872)
            + x33 * (x874 - x875 + x878)
            - x6 * (x5 * x879 - x6 * (x5 * x877 - x6 * x881) + x880)
        )
    )
    result[4, 2] = numpy.sum(
        x532
        * (
            x3 * (x5 * x886 - x6 * x891 + x884)
            + x33 * (x886 - x887 + x890)
            - x6 * (x5 * x891 - x6 * (x5 * x889 - x6 * x893) + x892)
        )
    )
    result[4, 3] = numpy.sum(
        x532 * (x3 * (x5 * x895 - x6 * x897) - x6 * (x5 * x897 - x6 * x898) + x896)
    )
    result[4, 4] = numpy.sum(
        181.019335983756
        * x83
        * (x3 * (x5 * x900 - x6 * x902) - x6 * (x5 * x902 - x6 * x903) + x901)
    )
    result[4, 5] = numpy.sum(
        x532 * (x3 * (x5 * x905 - x6 * x907) - x6 * (x5 * x907 - x6 * x908) + x906)
    )
    result[4, 6] = numpy.sum(x291 * (x3 * x909 - x6 * x910))
    result[4, 7] = numpy.sum(x532 * (x3 * x911 - x6 * x912))
    result[4, 8] = numpy.sum(x532 * (x3 * x913 - x6 * x914))
    result[4, 9] = numpy.sum(x291 * (x3 * x915 - x6 * x916))
    result[5, 0] = numpy.sum(
        x148
        * (
            x3 * (x33 * (x922 - x923 + x927) + x5 * x929 - x6 * x939)
            - x6
            * (
                x33 * (x928 - x932 + x936)
                + x5 * x939
                - x6 * (x5 * x937 - x6 * (x5 * x935 - x6 * x941) + x940)
            )
            + x75 * (x929 - x931 + x938 + x942)
        )
    )
    result[5, 1] = numpy.sum(
        x480
        * (
            x3 * (x5 * x950 - x6 * x957 + x947)
            + x33 * (x950 - x951 + x956)
            - x6 * (x5 * x957 - x6 * (x5 * x955 - x6 * x960) + x958)
        )
    )
    result[5, 2] = numpy.sum(
        x480
        * (
            x3 * (x5 * x966 - x6 * x972 + x964)
            + x33 * (x966 - x967 + x971)
            - x6 * (x5 * x972 - x6 * (x5 * x970 - x6 * x974) + x973)
        )
    )
    result[5, 3] = numpy.sum(
        x480 * (x3 * (x5 * x977 - x6 * x982) - x6 * (x5 * x982 - x6 * x984) + x980)
    )
    result[5, 4] = numpy.sum(
        x532 * (x3 * (x5 * x988 - x6 * x991) - x6 * (x5 * x991 - x6 * x993) + x989)
    )
    result[5, 5] = numpy.sum(
        -x480 * (-x1000 + x3 * (x1001 * x6 - x5 * x999) + x6 * (x1001 * x5 - x1003 * x6))
    )
    result[5, 6] = numpy.sum(x148 * (x1005 * x3 - x1007 * x6))
    result[5, 7] = numpy.sum(x480 * (x1009 * x3 - x1011 * x6))
    result[5, 8] = numpy.sum(x480 * (x1013 * x3 - x1015 * x6))
    result[5, 9] = numpy.sum(x148 * (x1016 * x3 - x1017 * x6))
    result[6, 0] = numpy.sum(
        x84
        * (
            x33 * (x1023 - x1024 + x1029)
            + x5 * (x1023 * x5 - x1030 * x6 + x4 * (x1018 - x1019 - x1020 + x1021))
            - x6
            * (
                x1030 * x5
                + x4 * (x1022 - x1025 + x1026 - x1027)
                - x6 * (x1028 * x5 - x6 * (x415 * x779 + 2.0 * x453 - x785 * x8))
            )
        )
    )
    result[6, 1] = numpy.sum(
        x148
        * (
            x4 * (x1031 - x1032 + x1033 + x786 + x818)
            + x5 * (x1033 * x5 - x1034 * x6)
            - x6 * (x1034 * x5 - x6 * (x415 * x796 + x784 - x8 * x800 + x823))
        )
    )
    result[6, 2] = numpy.sum(
        x148
        * (
            x4 * (x1035 - x1036 - x1037 + x1038)
            + x5 * (x1038 * x5 - x1039 * x6)
            - x6 * (x1039 * x5 - x6 * (x415 * x810 + 2.0 * x504 - x8 * x814))
        )
    )
    result[6, 3] = numpy.sum(
        x148
        * (
            x5 * (x415 * x820 + 2.0 * x515 + 2.0 * x790 - x8 * x822)
            - x6 * (x415 * x822 + 2.0 * x517 + 2.0 * x799 - x8 * x824)
        )
    )
    result[6, 4] = numpy.sum(
        x291
        * (
            x5 * (x415 * x827 - x8 * x829 + x804 + x839)
            - x6 * (x415 * x829 - x8 * x830 + x813 + x841)
        )
    )
    result[6, 5] = numpy.sum(
        x148
        * (
            x5 * (x415 * x833 + 2.0 * x542 - x8 * x835)
            - x6 * (x415 * x835 + 2.0 * x544 - x8 * x836)
        )
    )
    result[6, 6] = numpy.sum(x84 * (x415 * x837 + 2.0 * x550 - x8 * x838 + 3.0 * x821))
    result[6, 7] = numpy.sum(x148 * (x415 * x840 + 2.0 * x557 - x8 * x842 + 2.0 * x828))
    result[6, 8] = numpy.sum(x148 * (x415 * x843 + 2.0 * x563 - x8 * x844 + x834))
    result[6, 9] = numpy.sum(x84 * (x415 * x845 + 2.0 * x569 - x8 * x846))
    result[7, 0] = numpy.sum(
        x148
        * (
            x33 * (x1044 - x1045 + x1049)
            + x5 * (x1044 * x5 - x1050 * x6 + x4 * (x1040 - x1041 + x1042 + x609))
            - x6
            * (
                x1050 * x5
                + x4 * (x1043 - x1046 + x1047 + x612)
                - x6 * (x1048 * x5 - x6 * (x415 * x862 + x610 - x8 * x868))
            )
        )
    )
    result[7, 1] = numpy.sum(
        x480
        * (
            x4 * (x1051 - x1052 + x1053 + x640 + x869)
            + x5 * (x1053 * x5 - x1054 * x6)
            - x6 * (x1054 * x5 - x6 * (x415 * x877 + x643 - x8 * x881 + x867))
        )
    )
    result[7, 2] = numpy.sum(
        x480
        * (
            x4 * (x1055 - x1056 + x1057 + x666)
            + x5 * (x1057 * x5 - x1058 * x6)
            - x6 * (x1058 * x5 - x6 * (x415 * x889 + x669 - x8 * x893))
        )
    )
    result[7, 3] = numpy.sum(
        x480
        * (
            x5 * (x415 * x895 + x686 - x8 * x897 + 2.0 * x872)
            - x6 * (x415 * x897 + x689 - x8 * x898 + 2.0 * x880)
        )
    )
    result[7, 4] = numpy.sum(
        x532
        * (
            x5 * (x415 * x900 + x706 - x8 * x902 + x884)
            - x6 * (x415 * x902 + x708 - x8 * x903 + x892)
        )
    )
    result[7, 5] = numpy.sum(
        x480
        * (x5 * (x415 * x905 + x720 - x8 * x907) - x6 * (x415 * x907 + x722 - x8 * x908))
    )
    result[7, 6] = numpy.sum(x148 * (x415 * x909 + x730 - x8 * x910 + 3.0 * x896))
    result[7, 7] = numpy.sum(x480 * (x415 * x911 + x740 - x8 * x912 + 2.0 * x901))
    result[7, 8] = numpy.sum(x480 * (x415 * x913 + x749 - x8 * x914 + x906))
    result[7, 9] = numpy.sum(x148 * (x415 * x915 + x758 - x8 * x916))
    result[8, 0] = numpy.sum(
        x148
        * (
            x33 * (x1062 - x1063 + x1066)
            + x5 * (x1062 * x5 - x1067 * x6 + x4 * (-x1059 + x1060 + x943))
            - x6
            * (
                x1067 * x5
                + x4 * (x1061 - x1064 + x953)
                - x6 * (x1065 * x5 - x6 * (x415 * x935 + x959))
            )
        )
    )
    result[8, 1] = numpy.sum(
        x480
        * (
            x4 * (-x1068 + x1069 + x979)
            + x5 * (x1069 * x5 - x1070 * x6)
            - x6 * (x1070 * x5 - x6 * (x415 * x955 + x983))
        )
    )
    result[8, 2] = numpy.sum(
        x480
        * (
            x4 * (-x1071 + x1072 + x985)
            + x5 * (x1072 * x5 - x1073 * x6)
            - x6 * (x1073 * x5 - x6 * (x415 * x970 + x992))
        )
    )
    result[8, 3] = numpy.sum(
        x480 * (x5 * (x1004 + x415 * x977) - x6 * (x1006 + x415 * x982))
    )
    result[8, 4] = numpy.sum(
        x532 * (x5 * (x1008 + x415 * x988) - x6 * (x1010 + x415 * x991))
    )
    result[8, 5] = numpy.sum(
        x480 * (x5 * (x1012 + x415 * x999) - x6 * (x1001 * x415 + x1014))
    )
    result[8, 6] = numpy.sum(x148 * (x1005 * x415 - x1007 * x8 + 3.0 * x980))
    result[8, 7] = numpy.sum(x480 * (x1009 * x415 - x1011 * x8 + 2.0 * x989))
    result[8, 8] = numpy.sum(x480 * (x1000 + x1013 * x415 - x1015 * x8))
    result[8, 9] = numpy.sum(x148 * (x1016 * x415 - x1017 * x8))
    result[9, 0] = numpy.sum(
        x84
        * (
            x33 * (x1080 - x1081 + x1086)
            + x5 * (x1078 + x1080 * x5 - x1087 * x6)
            - x6 * (x1087 * x5 + x1088 - x6 * (x1085 * x5 - x1089 * x6))
        )
    )
    result[9, 1] = numpy.sum(
        x148 * (x1093 + x5 * (x1092 * x5 - x1094 * x6) - x6 * (x1094 * x5 - x1095 * x6))
    )
    result[9, 2] = numpy.sum(
        x148 * (x1099 + x5 * (x1098 * x5 - x1100 * x6) - x6 * (x1100 * x5 - x1101 * x6))
    )
    result[9, 3] = numpy.sum(x148 * (x1102 * x5 - x1103 * x6))
    result[9, 4] = numpy.sum(x291 * (x1104 * x5 - x1105 * x6))
    result[9, 5] = numpy.sum(x148 * (x1106 * x5 - x1107 * x6))
    result[9, 6] = numpy.sum(x84 * (2.0 * x1093 + x1102 * x87 - x1103 * x8))
    result[9, 7] = numpy.sum(x148 * (x1099 + x1104 * x87 - x1105 * x8))
    result[9, 8] = numpy.sum(x148 * (x1106 * x87 - x1107 * x8))
    result[9, 9] = numpy.sum(
        x84 * (-x10 * x1017 + 3.0 * x1000 + x1016 * x572 + 2.0 * x758)
    )
    return result


def coulomb3d_34(ax, da, A, bx, db, B, C):
    """Cartesian (fg) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 15), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = 0.5 / (ax + bx)
    x5 = -x2 - B[0]
    x6 = -x2 - C[0]
    x7 = -x1 * (ax * A[1] + bx * B[1])
    x8 = -x7 - C[1]
    x9 = -x1 * (ax * A[2] + bx * B[2])
    x10 = -x9 - C[2]
    x11 = x0 * (x10**2 + x6**2 + x8**2)
    x12 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x13 = x12 * boys(2, x11)
    x14 = x12 * boys(1, x11)
    x15 = x4 * (-x13 + x14)
    x16 = x13 * x6
    x17 = x14 * x5
    x18 = -x16 + x17
    x19 = x18 * x5
    x20 = x12 * boys(3, x11)
    x21 = x20 * x6
    x22 = x13 * x5
    x23 = -x21 + x22
    x24 = x23 * x6
    x25 = x15 + x19 - x24
    x26 = x25 * x5
    x27 = x4 * (x13 - x20)
    x28 = x23 * x5
    x29 = x12 * boys(4, x11)
    x30 = x29 * x6
    x31 = x20 * x5
    x32 = -x30 + x31
    x33 = x32 * x6
    x34 = x27 + x28 - x33
    x35 = x34 * x6
    x36 = 2.0 * x4
    x37 = x36 * (x18 + x21 - x22)
    x38 = -x35 + x37
    x39 = x26 + x38
    x40 = x39 * x5
    x41 = x12 * boys(0, x11)
    x42 = x4 * (-x14 + x41)
    x43 = -x14 * x6 + x41 * x5
    x44 = -x18 * x6 + x42 + x43 * x5
    x45 = -x25 * x6 + x36 * (x16 - x17 + x43)
    x46 = x44 * x5 + x45
    x47 = -x15
    x48 = 3.0 * x4
    x49 = -x39 * x6 + x48 * (-x19 + x24 + x44 + x47)
    x50 = x46 * x5 + x49
    x51 = x34 * x5
    x52 = x4 * (x20 - x29)
    x53 = x32 * x5
    x54 = x12 * boys(5, x11)
    x55 = x54 * x6
    x56 = x29 * x5
    x57 = -x55 + x56
    x58 = x57 * x6
    x59 = x52 + x53 - x58
    x60 = x59 * x6
    x61 = x36 * (x23 + x30 - x31)
    x62 = x51 - x60 + x61
    x63 = x6 * x62
    x64 = -x27
    x65 = x48 * (x25 - x28 + x33 + x64)
    x66 = x63 - x65
    x67 = -x63 + x65
    x68 = x40 + x67
    x69 = x35 - x37
    x70 = x4 * (-x26 + x46 + x69)
    x71 = x3 * x50 - x6 * x68 + 4.0 * x70
    x72 = x3 * x68
    x73 = x5 * x62
    x74 = x5 * x59
    x75 = x4 * (x29 - x54)
    x76 = x5 * x57
    x77 = x12 * boys(6, x11)
    x78 = x6 * x77
    x79 = x5 * x54
    x80 = -x78 + x79
    x81 = x6 * x80
    x82 = x75 + x76 - x81
    x83 = x6 * x82
    x84 = x36 * (x32 + x55 - x56)
    x85 = x74 - x83 + x84
    x86 = x6 * x85
    x87 = -x52
    x88 = x48 * (x34 - x53 + x58 + x87)
    x89 = -x86 + x88
    x90 = x73 + x89
    x91 = x6 * x90
    x92 = x60 - x61
    x93 = x4 * (x39 - x51 + x92)
    x94 = 4.0 * x93
    x95 = x72 - x91 + x94
    x96 = x3 * x39
    x97 = x3 * x46 + x49
    x98 = 4.0 * x4
    x99 = x86 - x88
    x100 = x4 * (x54 - x77)
    x101 = x12 * boys(7, x11)
    x102 = -x75
    x103 = x3 * x62
    x104 = x67 + x96
    x105 = x25 * x3
    x106 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**5.5)
    x107 = 9.12251705727742 * x106
    x108 = x13 * x8
    x109 = -x108
    x110 = -x7 - B[1]
    x111 = x110 * x14
    x112 = x109 + x111
    x113 = x112 * x5
    x114 = x20 * x8
    x115 = -x114
    x116 = x110 * x13
    x117 = x115 + x116
    x118 = x117 * x6
    x119 = x113 - x118
    x120 = x119 * x5
    x121 = x4 * (x112 + x114 - x116)
    x122 = x117 * x5
    x123 = x29 * x8
    x124 = -x123
    x125 = x110 * x20
    x126 = x124 + x125
    x127 = x126 * x6
    x128 = x122 - x127
    x129 = x128 * x6
    x130 = x121 - x129
    x131 = x120 + x130
    x132 = x131 * x5
    x133 = -x14 * x8
    x134 = x110 * x41 + x133
    x135 = -x112 * x6 + x134 * x5
    x136 = x4 * (x108 - x111 + x134)
    x137 = -x119 * x6 + x136
    x138 = x135 * x5 + x137
    x139 = -x131 * x6 + x36 * (-x113 + x118 + x135)
    x140 = x138 * x5 + x139
    x141 = x4 * (x117 + x123 - x125)
    x142 = x128 * x5
    x143 = x126 * x5
    x144 = x54 * x8
    x145 = -x144
    x146 = x110 * x29
    x147 = x145 + x146
    x148 = x147 * x6
    x149 = x143 - x148
    x150 = x149 * x6
    x151 = x141 + x142 - x150
    x152 = x151 * x6
    x153 = x36 * (x119 - x122 + x127)
    x154 = x152 - x153
    x155 = -x152 + x153
    x156 = x132 + x155
    x157 = -x121
    x158 = x129 + x157
    x159 = x4 * (-x120 + x138 + x158)
    x160 = x140 * x3 - x156 * x6 + 3.0 * x159
    x161 = x156 * x3
    x162 = x151 * x5
    x163 = x4 * (x126 + x144 - x146)
    x164 = x149 * x5
    x165 = x147 * x5
    x166 = x77 * x8
    x167 = -x166
    x168 = x110 * x54
    x169 = x167 + x168
    x170 = x169 * x6
    x171 = x165 - x170
    x172 = x171 * x6
    x173 = x163 + x164 - x172
    x174 = x173 * x6
    x175 = x36 * (x128 - x143 + x148)
    x176 = -x174 + x175
    x177 = x162 + x176
    x178 = x177 * x6
    x179 = -x141
    x180 = x150 + x179
    x181 = x4 * (x131 - x142 + x180)
    x182 = 3.0 * x181
    x183 = x161 - x178 + x182
    x184 = x131 * x3
    x185 = x138 * x3 + x139
    x186 = x174 - x175
    x187 = x4 * (x147 + x166 - x168)
    x188 = -x101 * x8
    x189 = x110 * x77 + x188
    x190 = -x163
    x191 = x151 * x3
    x192 = x155 + x184
    x193 = x119 * x3
    x194 = 24.1359114645008 * x106
    x195 = x10 * x13
    x196 = -x195
    x197 = -x9 - B[2]
    x198 = x14 * x197
    x199 = x196 + x198
    x200 = x199 * x5
    x201 = x10 * x20
    x202 = -x201
    x203 = x13 * x197
    x204 = x202 + x203
    x205 = x204 * x6
    x206 = x200 - x205
    x207 = x206 * x5
    x208 = x4 * (x199 + x201 - x203)
    x209 = x204 * x5
    x210 = x10 * x29
    x211 = -x210
    x212 = x197 * x20
    x213 = x211 + x212
    x214 = x213 * x6
    x215 = x209 - x214
    x216 = x215 * x6
    x217 = x208 - x216
    x218 = x207 + x217
    x219 = x218 * x5
    x220 = -x10 * x14
    x221 = x197 * x41 + x220
    x222 = -x199 * x6 + x221 * x5
    x223 = x4 * (x195 - x198 + x221)
    x224 = -x206 * x6 + x223
    x225 = x222 * x5 + x224
    x226 = -x218 * x6 + x36 * (-x200 + x205 + x222)
    x227 = x225 * x5 + x226
    x228 = x4 * (x204 + x210 - x212)
    x229 = x215 * x5
    x230 = x213 * x5
    x231 = x10 * x54
    x232 = -x231
    x233 = x197 * x29
    x234 = x232 + x233
    x235 = x234 * x6
    x236 = x230 - x235
    x237 = x236 * x6
    x238 = x228 + x229 - x237
    x239 = x238 * x6
    x240 = x36 * (x206 - x209 + x214)
    x241 = x239 - x240
    x242 = -x239 + x240
    x243 = x219 + x242
    x244 = -x208
    x245 = x216 + x244
    x246 = x4 * (-x207 + x225 + x245)
    x247 = x227 * x3 - x243 * x6 + 3.0 * x246
    x248 = x243 * x3
    x249 = x238 * x5
    x250 = x4 * (x213 + x231 - x233)
    x251 = x236 * x5
    x252 = x234 * x5
    x253 = x10 * x77
    x254 = -x253
    x255 = x197 * x54
    x256 = x254 + x255
    x257 = x256 * x6
    x258 = x252 - x257
    x259 = x258 * x6
    x260 = x250 + x251 - x259
    x261 = x260 * x6
    x262 = x36 * (x215 - x230 + x235)
    x263 = -x261 + x262
    x264 = x249 + x263
    x265 = x264 * x6
    x266 = -x228
    x267 = x237 + x266
    x268 = x4 * (x218 - x229 + x267)
    x269 = 3.0 * x268
    x270 = x248 - x265 + x269
    x271 = x218 * x3
    x272 = x225 * x3 + x226
    x273 = x261 - x262
    x274 = x4 * (x234 + x253 - x255)
    x275 = -x10 * x101
    x276 = x197 * x77 + x275
    x277 = -x250
    x278 = x238 * x3
    x279 = x242 + x271
    x280 = x206 * x3
    x281 = x110 * x112
    x282 = x117 * x8
    x283 = x15 - x282
    x284 = x281 + x283
    x285 = x284 * x5
    x286 = x110 * x117
    x287 = x126 * x8
    x288 = x27 - x287
    x289 = x286 + x288
    x290 = x289 * x6
    x291 = -x290
    x292 = x285 + x291
    x293 = x292 * x5
    x294 = -x112 * x8 + x42
    x295 = x110 * x134 + x294
    x296 = -x284 * x6
    x297 = x295 * x5 + x296
    x298 = x282 + x47
    x299 = x4 * (-x281 + x295 + x298)
    x300 = -x292 * x6 + x299
    x301 = x297 * x5 + x300
    x302 = x289 * x5
    x303 = x110 * x126
    x304 = x147 * x8
    x305 = -x304 + x52
    x306 = x303 + x305
    x307 = x306 * x6
    x308 = x302 - x307
    x309 = x308 * x6
    x310 = x287 + x64
    x311 = x4 * (x284 - x286 + x310)
    x312 = -x311
    x313 = x309 + x312
    x314 = -x309 + x311
    x315 = x293 + x314
    x316 = x4 * (-x285 + x290 + x297)
    x317 = x3 * x301 - x315 * x6 + 2.0 * x316
    x318 = x3 * x315
    x319 = x308 * x5
    x320 = x304 + x87
    x321 = x4 * (x289 - x303 + x320)
    x322 = x306 * x5
    x323 = x110 * x147
    x324 = x169 * x8
    x325 = -x324 + x75
    x326 = x323 + x325
    x327 = x326 * x6
    x328 = x322 - x327
    x329 = x328 * x6
    x330 = x321 - x329
    x331 = x319 + x330
    x332 = x331 * x6
    x333 = x4 * (x292 - x302 + x307)
    x334 = 2.0 * x333
    x335 = x318 - x332 + x334
    x336 = x292 * x3
    x337 = x297 * x3 + x300
    x338 = -x321
    x339 = x329 + x338
    x340 = x102 + x324
    x341 = x4 * (x306 - x323 + x340)
    x342 = x100 - x189 * x8
    x343 = x110 * x169 + x342
    x344 = x3 * x308
    x345 = x314 + x336
    x346 = x284 * x3
    x347 = 31.1593277158494 * x106
    x348 = x110 * x199
    x349 = x204 * x8
    x350 = -x349
    x351 = x348 + x350
    x352 = x351 * x5
    x353 = x110 * x204
    x354 = x213 * x8
    x355 = -x354
    x356 = x353 + x355
    x357 = x356 * x6
    x358 = -x357
    x359 = x352 + x358
    x360 = x359 * x5
    x361 = -x199 * x8
    x362 = x110 * x221 + x361
    x363 = -x351 * x6
    x364 = x362 * x5 + x363
    x365 = x4 * (-x348 + x349 + x362)
    x366 = -x359 * x6 + x365
    x367 = x364 * x5 + x366
    x368 = x356 * x5
    x369 = x110 * x213
    x370 = x234 * x8
    x371 = -x370
    x372 = x369 + x371
    x373 = x372 * x6
    x374 = x368 - x373
    x375 = x374 * x6
    x376 = x4 * (x351 - x353 + x354)
    x377 = -x376
    x378 = x375 + x377
    x379 = -x375 + x376
    x380 = x360 + x379
    x381 = x4 * (-x352 + x357 + x364)
    x382 = x3 * x367 - x380 * x6 + 2.0 * x381
    x383 = x3 * x380
    x384 = x374 * x5
    x385 = x4 * (x356 - x369 + x370)
    x386 = x372 * x5
    x387 = x110 * x234
    x388 = x256 * x8
    x389 = -x388
    x390 = x387 + x389
    x391 = x390 * x6
    x392 = x386 - x391
    x393 = x392 * x6
    x394 = x385 - x393
    x395 = x384 + x394
    x396 = x395 * x6
    x397 = x4 * (x359 - x368 + x373)
    x398 = 2.0 * x397
    x399 = x383 - x396 + x398
    x400 = x3 * x359
    x401 = x3 * x364 + x366
    x402 = -x385
    x403 = x393 + x402
    x404 = x4 * (x372 - x387 + x388)
    x405 = -x276 * x8
    x406 = x110 * x256 + x405
    x407 = x3 * x374
    x408 = x379 + x400
    x409 = x3 * x351
    x410 = 53.9695387335403 * x106
    x411 = x197 * x199
    x412 = x10 * x204
    x413 = x15 - x412
    x414 = x411 + x413
    x415 = x414 * x5
    x416 = x197 * x204
    x417 = x10 * x213
    x418 = x27 - x417
    x419 = x416 + x418
    x420 = x419 * x6
    x421 = -x420
    x422 = x415 + x421
    x423 = x422 * x5
    x424 = -x10 * x199 + x42
    x425 = x197 * x221 + x424
    x426 = -x414 * x6
    x427 = x425 * x5 + x426
    x428 = x412 + x47
    x429 = x4 * (-x411 + x425 + x428)
    x430 = -x422 * x6 + x429
    x431 = x427 * x5 + x430
    x432 = x419 * x5
    x433 = x197 * x213
    x434 = x10 * x234
    x435 = -x434 + x52
    x436 = x433 + x435
    x437 = x436 * x6
    x438 = x432 - x437
    x439 = x438 * x6
    x440 = x417 + x64
    x441 = x4 * (x414 - x416 + x440)
    x442 = -x441
    x443 = x439 + x442
    x444 = -x439 + x441
    x445 = x423 + x444
    x446 = x4 * (-x415 + x420 + x427)
    x447 = x3 * x431 - x445 * x6 + 2.0 * x446
    x448 = x3 * x445
    x449 = x438 * x5
    x450 = x434 + x87
    x451 = x4 * (x419 - x433 + x450)
    x452 = x436 * x5
    x453 = x197 * x234
    x454 = x10 * x256
    x455 = -x454 + x75
    x456 = x453 + x455
    x457 = x456 * x6
    x458 = x452 - x457
    x459 = x458 * x6
    x460 = x451 - x459
    x461 = x449 + x460
    x462 = x461 * x6
    x463 = x4 * (x422 - x432 + x437)
    x464 = 2.0 * x463
    x465 = x448 - x462 + x464
    x466 = x3 * x422
    x467 = x3 * x427 + x430
    x468 = -x451
    x469 = x459 + x468
    x470 = x102 + x454
    x471 = x4 * (x436 - x453 + x470)
    x472 = -x10 * x276 + x100
    x473 = x197 * x256 + x472
    x474 = x3 * x438
    x475 = x444 + x466
    x476 = x3 * x414
    x477 = x110 * x284
    x478 = 2.0 * x136 - x284 * x8
    x479 = x110 * x295 + x478
    x480 = x289 * x8
    x481 = 2.0 * x121
    x482 = x480 - x481
    x483 = x4 * (-x477 + x479 + x482)
    x484 = -x480 + x481
    x485 = x477 + x484
    x486 = -x485 * x6
    x487 = x3 * x479 + x486
    x488 = x110 * x289
    x489 = x306 * x8
    x490 = 2.0 * x141
    x491 = -x489 + x490
    x492 = x488 + x491
    x493 = x3 * x492
    x494 = x110 * x306
    x495 = x326 * x8
    x496 = 2.0 * x163
    x497 = -x495 + x496
    x498 = x494 + x497
    x499 = x498 * x6
    x500 = -x499
    x501 = x489 - x490
    x502 = x4 * (x485 - x488 + x501)
    x503 = -x502
    x504 = x3 * x485
    x505 = x492 * x6
    x506 = -x505
    x507 = x504 + x506
    x508 = x485 * x5
    x509 = x479 * x5 + x486
    x510 = x506 + x508
    x511 = x3 * x509 + x483 - x510 * x6
    x512 = x3 * x510
    x513 = x492 * x5
    x514 = x500 + x513
    x515 = x514 * x6
    x516 = x502 + x512 - x515
    x517 = x495 - x496
    x518 = x4 * (x492 - x494 + x517)
    x519 = 2.0 * x187 - x343 * x8
    x520 = x110 * x326 + x519
    x521 = x110 * x351
    x522 = x223 - x351 * x8
    x523 = x110 * x362 + x522
    x524 = x356 * x8
    x525 = x244 + x524
    x526 = x4 * (-x521 + x523 + x525)
    x527 = x208 - x524
    x528 = x521 + x527
    x529 = -x528 * x6
    x530 = x3 * x523 + x529
    x531 = x110 * x356
    x532 = x372 * x8
    x533 = x228 - x532
    x534 = x531 + x533
    x535 = x3 * x534
    x536 = x110 * x372
    x537 = x390 * x8
    x538 = x250 - x537
    x539 = x536 + x538
    x540 = x539 * x6
    x541 = -x540
    x542 = x266 + x532
    x543 = x4 * (x528 - x531 + x542)
    x544 = -x543
    x545 = x3 * x528
    x546 = x534 * x6
    x547 = -x546
    x548 = x545 + x547
    x549 = x5 * x528
    x550 = x5 * x523 + x529
    x551 = x547 + x549
    x552 = x3 * x550 + x526 - x551 * x6
    x553 = x3 * x551
    x554 = x5 * x534
    x555 = x541 + x554
    x556 = x555 * x6
    x557 = x543 + x553 - x556
    x558 = x277 + x537
    x559 = x4 * (x534 - x536 + x558)
    x560 = x274 - x406 * x8
    x561 = x110 * x390 + x560
    x562 = x419 * x8
    x563 = x110 * x414
    x564 = -x414 * x8
    x565 = x110 * x425 + x564
    x566 = x4 * (x562 - x563 + x565)
    x567 = -x562
    x568 = x563 + x567
    x569 = -x568 * x6
    x570 = x3 * x565 + x569
    x571 = x110 * x419
    x572 = x436 * x8
    x573 = -x572
    x574 = x571 + x573
    x575 = x3 * x574
    x576 = x110 * x436
    x577 = x456 * x8
    x578 = -x577
    x579 = x576 + x578
    x580 = x579 * x6
    x581 = -x580
    x582 = x4 * (x568 - x571 + x572)
    x583 = -x582
    x584 = x3 * x568
    x585 = x574 * x6
    x586 = -x585
    x587 = x584 + x586
    x588 = x5 * x568
    x589 = x5 * x565 + x569
    x590 = x586 + x588
    x591 = x3 * x589 + x566 - x590 * x6
    x592 = x3 * x590
    x593 = x5 * x574
    x594 = x581 + x593
    x595 = x594 * x6
    x596 = x582 + x592 - x595
    x597 = x4 * (x574 - x576 + x577)
    x598 = -x473 * x8
    x599 = x110 * x456 + x598
    x600 = x197 * x414
    x601 = -x10 * x414 + 2.0 * x223
    x602 = x197 * x425 + x601
    x603 = x10 * x419
    x604 = 2.0 * x208
    x605 = x603 - x604
    x606 = x4 * (-x600 + x602 + x605)
    x607 = -x603 + x604
    x608 = x600 + x607
    x609 = -x6 * x608
    x610 = x3 * x602 + x609
    x611 = x197 * x419
    x612 = x10 * x436
    x613 = 2.0 * x228
    x614 = -x612 + x613
    x615 = x611 + x614
    x616 = x3 * x615
    x617 = x197 * x436
    x618 = x10 * x456
    x619 = 2.0 * x250
    x620 = -x618 + x619
    x621 = x617 + x620
    x622 = x6 * x621
    x623 = -x622
    x624 = x612 - x613
    x625 = x4 * (x608 - x611 + x624)
    x626 = -x625
    x627 = x3 * x608
    x628 = x6 * x615
    x629 = -x628
    x630 = x627 + x629
    x631 = x5 * x608
    x632 = x5 * x602 + x609
    x633 = x629 + x631
    x634 = x3 * x632 - x6 * x633 + x606
    x635 = x3 * x633
    x636 = x5 * x615
    x637 = x623 + x636
    x638 = x6 * x637
    x639 = x625 + x635 - x638
    x640 = x618 - x619
    x641 = x4 * (x615 - x617 + x640)
    x642 = -x10 * x473 + 2.0 * x274
    x643 = x197 * x456 + x642
    x644 = x110 * x485
    x645 = 3.0 * x299 - x485 * x8
    x646 = x110 * x479 + x645
    x647 = x492 * x8
    x648 = 3.0 * x311
    x649 = x647 - x648
    x650 = x4 * (-x644 + x646 + x649)
    x651 = -x647 + x648
    x652 = x644 + x651
    x653 = x3 * x646 - x6 * x652
    x654 = x3 * x652
    x655 = x110 * x492
    x656 = x498 * x8
    x657 = 3.0 * x321
    x658 = -x656 + x657
    x659 = x655 + x658
    x660 = x6 * x659
    x661 = x654 - x660
    x662 = x656 - x657
    x663 = x4 * (x652 - x655 + x662)
    x664 = 3.0 * x341 - x520 * x8
    x665 = x110 * x498 + x664
    x666 = x110 * x528
    x667 = 2.0 * x365 - x528 * x8
    x668 = x110 * x523 + x667
    x669 = x534 * x8
    x670 = 2.0 * x376
    x671 = x669 - x670
    x672 = x4 * (-x666 + x668 + x671)
    x673 = -x669 + x670
    x674 = x666 + x673
    x675 = x3 * x668 - x6 * x674
    x676 = x3 * x674
    x677 = x110 * x534
    x678 = x539 * x8
    x679 = 2.0 * x385
    x680 = -x678 + x679
    x681 = x677 + x680
    x682 = x6 * x681
    x683 = x676 - x682
    x684 = x678 - x679
    x685 = x4 * (x674 - x677 + x684)
    x686 = 2.0 * x404 - x561 * x8
    x687 = x110 * x539 + x686
    x688 = x110 * x568
    x689 = x429 - x568 * x8
    x690 = x110 * x565 + x689
    x691 = x574 * x8
    x692 = x442 + x691
    x693 = x4 * (-x688 + x690 + x692)
    x694 = x441 - x691
    x695 = x688 + x694
    x696 = x3 * x690 - x6 * x695
    x697 = x3 * x695
    x698 = x110 * x574
    x699 = x579 * x8
    x700 = x451 - x699
    x701 = x698 + x700
    x702 = x6 * x701
    x703 = x697 - x702
    x704 = x468 + x699
    x705 = x4 * (x695 - x698 + x704)
    x706 = x471 - x599 * x8
    x707 = x110 * x579 + x706
    x708 = x615 * x8
    x709 = x110 * x608
    x710 = -x608 * x8
    x711 = x110 * x602 + x710
    x712 = x4 * (x708 - x709 + x711)
    x713 = -x708
    x714 = x709 + x713
    x715 = x3 * x711 - x6 * x714
    x716 = x3 * x714
    x717 = x110 * x615
    x718 = x621 * x8
    x719 = -x718
    x720 = x717 + x719
    x721 = x6 * x720
    x722 = x716 - x721
    x723 = x4 * (x714 - x717 + x718)
    x724 = -x643 * x8
    x725 = x110 * x621 + x724
    x726 = x197 * x608
    x727 = -x10 * x608 + 3.0 * x429
    x728 = x197 * x602 + x727
    x729 = x10 * x615
    x730 = 3.0 * x441
    x731 = x729 - x730
    x732 = x4 * (-x726 + x728 + x731)
    x733 = -x729 + x730
    x734 = x726 + x733
    x735 = x3 * x728 - x6 * x734
    x736 = x3 * x734
    x737 = x197 * x615
    x738 = x10 * x621
    x739 = 3.0 * x451
    x740 = -x738 + x739
    x741 = x737 + x740
    x742 = x6 * x741
    x743 = x736 - x742
    x744 = x738 - x739
    x745 = x4 * (x734 - x737 + x744)
    x746 = -x10 * x643 + 3.0 * x471
    x747 = x197 * x621 + x746
    x748 = -x7 - A[1]
    x749 = x13 * x748
    x750 = x14 * x748
    x751 = x109 + x750
    x752 = x4 * (x114 - x749 + x751)
    x753 = x5 * x751
    x754 = x115 + x749
    x755 = x6 * x754
    x756 = x753 - x755
    x757 = x5 * x756
    x758 = x5 * x754
    x759 = x20 * x748
    x760 = x124 + x759
    x761 = x6 * x760
    x762 = x758 - x761
    x763 = x6 * x762
    x764 = x752 + x757 - x763
    x765 = x5 * x764
    x766 = x4 * (x123 + x754 - x759)
    x767 = x5 * x762
    x768 = x5 * x760
    x769 = x29 * x748
    x770 = x145 + x769
    x771 = x6 * x770
    x772 = x768 - x771
    x773 = x6 * x772
    x774 = x766 + x767 - x773
    x775 = x6 * x774
    x776 = x36 * (x756 - x758 + x761)
    x777 = x765 - x775 + x776
    x778 = x5 * x777
    x779 = x133 + x41 * x748
    x780 = x4 * (x108 - x750 + x779)
    x781 = x5 * x779 - x6 * x751
    x782 = x5 * x781 - x6 * x756 + x780
    x783 = x36 * (-x753 + x755 + x781) + x5 * x782 - x6 * x764
    x784 = -x752
    x785 = x48 * (-x757 + x763 + x782 + x784) - x6 * x777
    x786 = x5 * x783 + x785
    x787 = x5 * x774
    x788 = x4 * (x144 + x760 - x769)
    x789 = x5 * x772
    x790 = x5 * x770
    x791 = x54 * x748
    x792 = x167 + x791
    x793 = x6 * x792
    x794 = x790 - x793
    x795 = x6 * x794
    x796 = x788 + x789 - x795
    x797 = x6 * x796
    x798 = x36 * (x762 - x768 + x771)
    x799 = x787 - x797 + x798
    x800 = x6 * x799
    x801 = -x766
    x802 = x48 * (x764 - x767 + x773 + x801)
    x803 = x800 - x802
    x804 = x778 - x800 + x802
    x805 = x4 * (x166 + x770 - x791)
    x806 = x188 + x748 * x77
    x807 = -x788
    x808 = 20.3985682659737 * x106
    x809 = x117 * x748
    x810 = x112 * x748
    x811 = x283 + x810
    x812 = x4 * (x310 - x809 + x811)
    x813 = x5 * x811
    x814 = x288 + x809
    x815 = x6 * x814
    x816 = x813 - x815
    x817 = x5 * x816
    x818 = x5 * x814
    x819 = x126 * x748
    x820 = x305 + x819
    x821 = x6 * x820
    x822 = x818 - x821
    x823 = x6 * x822
    x824 = x812 + x817 - x823
    x825 = x5 * x824
    x826 = x134 * x748 + x294
    x827 = x4 * (x298 - x810 + x826)
    x828 = x5 * x826 - x6 * x811
    x829 = x5 * x828 - x6 * x816 + x827
    x830 = x36 * (-x813 + x815 + x828) - x6 * x824
    x831 = x5 * x829 + x830
    x832 = x4 * (x320 + x814 - x819)
    x833 = x5 * x822
    x834 = x5 * x820
    x835 = x147 * x748
    x836 = x325 + x835
    x837 = x6 * x836
    x838 = x834 - x837
    x839 = x6 * x838
    x840 = x832 + x833 - x839
    x841 = x6 * x840
    x842 = x36 * (x816 - x818 + x821)
    x843 = x841 - x842
    x844 = x825 - x841 + x842
    x845 = x4 * (x340 + x820 - x835)
    x846 = x169 * x748 + x342
    x847 = x204 * x748
    x848 = x199 * x748
    x849 = x350 + x848
    x850 = x4 * (x354 - x847 + x849)
    x851 = x5 * x849
    x852 = x355 + x847
    x853 = x6 * x852
    x854 = x851 - x853
    x855 = x5 * x854
    x856 = x5 * x852
    x857 = x213 * x748
    x858 = x371 + x857
    x859 = x6 * x858
    x860 = x856 - x859
    x861 = x6 * x860
    x862 = x850 + x855 - x861
    x863 = x5 * x862
    x864 = x221 * x748 + x361
    x865 = x4 * (x349 - x848 + x864)
    x866 = x5 * x864 - x6 * x849
    x867 = x5 * x866 - x6 * x854 + x865
    x868 = x36 * (-x851 + x853 + x866) - x6 * x862
    x869 = x5 * x867 + x868
    x870 = x4 * (x370 + x852 - x857)
    x871 = x5 * x860
    x872 = x5 * x858
    x873 = x234 * x748
    x874 = x389 + x873
    x875 = x6 * x874
    x876 = x872 - x875
    x877 = x6 * x876
    x878 = x870 + x871 - x877
    x879 = x6 * x878
    x880 = x36 * (x854 - x856 + x859)
    x881 = x879 - x880
    x882 = x863 - x879 + x880
    x883 = -x850
    x884 = x4 * (x388 + x858 - x873)
    x885 = x256 * x748 + x405
    x886 = -x870
    x887 = x284 * x748
    x888 = x484 + x887
    x889 = x5 * x888
    x890 = x289 * x748
    x891 = x491 + x890
    x892 = x6 * x891
    x893 = x889 - x892
    x894 = x5 * x893
    x895 = x295 * x748 + x478
    x896 = x5 * x895 - x6 * x888
    x897 = x4 * (x482 - x887 + x895)
    x898 = -x6 * x893 + x897
    x899 = x5 * x896 + x898
    x900 = x5 * x891
    x901 = x306 * x748
    x902 = x497 + x901
    x903 = x6 * x902
    x904 = x900 - x903
    x905 = x6 * x904
    x906 = x4 * (x501 + x888 - x890)
    x907 = x905 - x906
    x908 = x894 - x905 + x906
    x909 = x4 * (x517 + x891 - x901)
    x910 = x326 * x748 + x519
    x911 = 69.6743749058326 * x106
    x912 = x351 * x748
    x913 = x527 + x912
    x914 = x5 * x913
    x915 = x356 * x748
    x916 = x533 + x915
    x917 = x6 * x916
    x918 = x914 - x917
    x919 = x5 * x918
    x920 = x362 * x748 + x522
    x921 = x5 * x920 - x6 * x913
    x922 = x4 * (x525 - x912 + x920)
    x923 = -x6 * x918 + x922
    x924 = x5 * x921 + x923
    x925 = x5 * x916
    x926 = x372 * x748
    x927 = x538 + x926
    x928 = x6 * x927
    x929 = x925 - x928
    x930 = x6 * x929
    x931 = x4 * (x542 + x913 - x915)
    x932 = x930 - x931
    x933 = x919 - x930 + x931
    x934 = x4 * (x558 + x916 - x926)
    x935 = x390 * x748 + x560
    x936 = 120.679557322504 * x106
    x937 = x414 * x748
    x938 = x567 + x937
    x939 = x5 * x938
    x940 = x419 * x748
    x941 = x573 + x940
    x942 = x6 * x941
    x943 = x939 - x942
    x944 = x5 * x943
    x945 = x425 * x748 + x564
    x946 = x5 * x945 - x6 * x938
    x947 = x4 * (x562 - x937 + x945)
    x948 = -x6 * x943 + x947
    x949 = x5 * x946 + x948
    x950 = x5 * x941
    x951 = x436 * x748
    x952 = x578 + x951
    x953 = x6 * x952
    x954 = x950 - x953
    x955 = x6 * x954
    x956 = x4 * (x572 + x938 - x940)
    x957 = -x956
    x958 = x955 + x957
    x959 = x944 - x955 + x956
    x960 = x4 * (x577 + x941 - x951)
    x961 = x456 * x748 + x598
    x962 = x479 * x748 + x645
    x963 = x492 * x748
    x964 = x658 + x963
    x965 = x6 * x964
    x966 = x485 * x748
    x967 = x651 + x966
    x968 = -x6 * x967
    x969 = x5 * x967
    x970 = x5 * x962 + x968
    x971 = x4 * (x649 + x962 - x966)
    x972 = -x965 + x969
    x973 = x4 * (x662 - x963 + x967)
    x974 = x498 * x748 + x664
    x975 = x523 * x748 + x667
    x976 = x534 * x748
    x977 = x680 + x976
    x978 = x6 * x977
    x979 = x528 * x748
    x980 = x673 + x979
    x981 = -x6 * x980
    x982 = x5 * x980
    x983 = x5 * x975 + x981
    x984 = x4 * (x671 + x975 - x979)
    x985 = -x978 + x982
    x986 = x4 * (x684 - x976 + x980)
    x987 = x539 * x748 + x686
    x988 = x565 * x748 + x689
    x989 = x574 * x748
    x990 = x700 + x989
    x991 = x6 * x990
    x992 = x568 * x748
    x993 = x694 + x992
    x994 = -x6 * x993
    x995 = x5 * x993
    x996 = x5 * x988 + x994
    x997 = x4 * (x692 + x988 - x992)
    x998 = -x991 + x995
    x999 = x4 * (x704 - x989 + x993)
    x1000 = x579 * x748 + x706
    x1001 = x602 * x748 + x710
    x1002 = x615 * x748
    x1003 = x1002 + x719
    x1004 = x1003 * x6
    x1005 = x608 * x748
    x1006 = x1005 + x713
    x1007 = -x1006 * x6
    x1008 = x1006 * x5
    x1009 = x1001 * x5 + x1007
    x1010 = x4 * (x1001 - x1005 + x708)
    x1011 = -x1004 + x1008
    x1012 = x4 * (-x1002 + x1006 + x718)
    x1013 = x621 * x748 + x724
    x1014 = x659 * x8
    x1015 = x652 * x748
    x1016 = 4.0 * x502
    x1017 = 4.0 * x483 + x646 * x748 - x652 * x8
    x1018 = x4 * (x1014 - x1015 - x1016 + x1017)
    x1019 = -x1014 + x1015 + x1016
    x1020 = 4.0 * x518 + x659 * x748 - x665 * x8
    x1021 = x681 * x8
    x1022 = x674 * x748
    x1023 = 3.0 * x543
    x1024 = 3.0 * x526 + x668 * x748 - x674 * x8
    x1025 = x4 * (x1021 - x1022 - x1023 + x1024)
    x1026 = -x1021 + x1022 + x1023
    x1027 = 3.0 * x559 + x681 * x748 - x687 * x8
    x1028 = x701 * x8
    x1029 = x695 * x748
    x1030 = 2.0 * x582
    x1031 = 2.0 * x566 + x690 * x748 - x695 * x8
    x1032 = x4 * (x1028 - x1029 - x1030 + x1031)
    x1033 = -x1028 + x1029 + x1030
    x1034 = 2.0 * x597 + x701 * x748 - x707 * x8
    x1035 = x720 * x8
    x1036 = x714 * x748
    x1037 = x606 + x711 * x748 - x714 * x8
    x1038 = x4 * (x1035 - x1036 + x1037 + x626)
    x1039 = -x1035 + x1036 + x625
    x1040 = x641 + x720 * x748 - x725 * x8
    x1041 = x741 * x8
    x1042 = x734 * x748
    x1043 = x728 * x748 - x734 * x8
    x1044 = x4 * (x1041 - x1042 + x1043)
    x1045 = -x1041 + x1042
    x1046 = x741 * x748 - x747 * x8
    x1047 = -x9 - A[2]
    x1048 = x1047 * x13
    x1049 = x1047 * x14
    x1050 = x1049 + x196
    x1051 = x4 * (-x1048 + x1050 + x201)
    x1052 = x1050 * x5
    x1053 = x1048 + x202
    x1054 = x1053 * x6
    x1055 = x1052 - x1054
    x1056 = x1055 * x5
    x1057 = x1053 * x5
    x1058 = x1047 * x20
    x1059 = x1058 + x211
    x1060 = x1059 * x6
    x1061 = x1057 - x1060
    x1062 = x1061 * x6
    x1063 = x1051 + x1056 - x1062
    x1064 = x1063 * x5
    x1065 = x4 * (x1053 - x1058 + x210)
    x1066 = x1061 * x5
    x1067 = x1059 * x5
    x1068 = x1047 * x29
    x1069 = x1068 + x232
    x1070 = x1069 * x6
    x1071 = x1067 - x1070
    x1072 = x1071 * x6
    x1073 = x1065 + x1066 - x1072
    x1074 = x1073 * x6
    x1075 = x36 * (x1055 - x1057 + x1060)
    x1076 = x1064 - x1074 + x1075
    x1077 = x1076 * x5
    x1078 = x1047 * x41 + x220
    x1079 = x4 * (-x1049 + x1078 + x195)
    x1080 = -x1050 * x6 + x1078 * x5
    x1081 = -x1055 * x6 + x1079 + x1080 * x5
    x1082 = -x1063 * x6 + x1081 * x5 + x36 * (-x1052 + x1054 + x1080)
    x1083 = -x1051
    x1084 = -x1076 * x6 + x48 * (-x1056 + x1062 + x1081 + x1083)
    x1085 = x1082 * x5 + x1084
    x1086 = x1073 * x5
    x1087 = x4 * (x1059 - x1068 + x231)
    x1088 = x1071 * x5
    x1089 = x1069 * x5
    x1090 = x1047 * x54
    x1091 = x1090 + x254
    x1092 = x1091 * x6
    x1093 = x1089 - x1092
    x1094 = x1093 * x6
    x1095 = x1087 + x1088 - x1094
    x1096 = x1095 * x6
    x1097 = x36 * (x1061 - x1067 + x1070)
    x1098 = x1086 - x1096 + x1097
    x1099 = x1098 * x6
    x1100 = -x1065
    x1101 = x48 * (x1063 - x1066 + x1072 + x1100)
    x1102 = x1099 - x1101
    x1103 = x1077 - x1099 + x1101
    x1104 = x4 * (x1069 - x1090 + x253)
    x1105 = x1047 * x77 + x275
    x1106 = -x1087
    x1107 = x1059 * x8
    x1108 = x1053 * x110
    x1109 = x1050 * x110
    x1110 = x1053 * x8
    x1111 = -x1110
    x1112 = x1109 + x1111
    x1113 = x4 * (x1107 - x1108 + x1112)
    x1114 = x1112 * x5
    x1115 = -x1107
    x1116 = x1108 + x1115
    x1117 = x1116 * x6
    x1118 = x1114 - x1117
    x1119 = x1118 * x5
    x1120 = x1116 * x5
    x1121 = x1059 * x110
    x1122 = x1069 * x8
    x1123 = -x1122
    x1124 = x1121 + x1123
    x1125 = x1124 * x6
    x1126 = x1120 - x1125
    x1127 = x1126 * x6
    x1128 = x1113 + x1119 - x1127
    x1129 = x1128 * x5
    x1130 = -x1050 * x8
    x1131 = x1078 * x110 + x1130
    x1132 = x4 * (-x1109 + x1110 + x1131)
    x1133 = -x1112 * x6 + x1131 * x5
    x1134 = -x1118 * x6 + x1132 + x1133 * x5
    x1135 = -x1128 * x6 + x36 * (-x1114 + x1117 + x1133)
    x1136 = x1134 * x5 + x1135
    x1137 = x4 * (x1116 - x1121 + x1122)
    x1138 = x1126 * x5
    x1139 = x1124 * x5
    x1140 = x1069 * x110
    x1141 = x1091 * x8
    x1142 = -x1141
    x1143 = x1140 + x1142
    x1144 = x1143 * x6
    x1145 = x1139 - x1144
    x1146 = x1145 * x6
    x1147 = x1137 + x1138 - x1146
    x1148 = x1147 * x6
    x1149 = x36 * (x1118 - x1120 + x1125)
    x1150 = x1148 - x1149
    x1151 = x1129 - x1148 + x1149
    x1152 = -x1113
    x1153 = x4 * (x1124 - x1140 + x1141)
    x1154 = -x1105 * x8
    x1155 = x1091 * x110 + x1154
    x1156 = -x1137
    x1157 = x1047 * x204
    x1158 = x1047 * x199
    x1159 = x1158 + x413
    x1160 = x4 * (-x1157 + x1159 + x440)
    x1161 = x1159 * x5
    x1162 = x1157 + x418
    x1163 = x1162 * x6
    x1164 = x1161 - x1163
    x1165 = x1164 * x5
    x1166 = x1162 * x5
    x1167 = x1047 * x213
    x1168 = x1167 + x435
    x1169 = x1168 * x6
    x1170 = x1166 - x1169
    x1171 = x1170 * x6
    x1172 = x1160 + x1165 - x1171
    x1173 = x1172 * x5
    x1174 = x1047 * x221 + x424
    x1175 = x4 * (-x1158 + x1174 + x428)
    x1176 = -x1159 * x6 + x1174 * x5
    x1177 = -x1164 * x6 + x1175 + x1176 * x5
    x1178 = -x1172 * x6 + x36 * (-x1161 + x1163 + x1176)
    x1179 = x1177 * x5 + x1178
    x1180 = x4 * (x1162 - x1167 + x450)
    x1181 = x1170 * x5
    x1182 = x1168 * x5
    x1183 = x1047 * x234
    x1184 = x1183 + x455
    x1185 = x1184 * x6
    x1186 = x1182 - x1185
    x1187 = x1186 * x6
    x1188 = x1180 + x1181 - x1187
    x1189 = x1188 * x6
    x1190 = x36 * (x1164 - x1166 + x1169)
    x1191 = x1189 - x1190
    x1192 = x1173 - x1189 + x1190
    x1193 = -x1160
    x1194 = x4 * (x1168 - x1183 + x470)
    x1195 = x1047 * x256 + x472
    x1196 = -x1180
    x1197 = x110 * x1112
    x1198 = x1116 * x8
    x1199 = x1051 - x1198
    x1200 = x1197 + x1199
    x1201 = x1200 * x5
    x1202 = x110 * x1116
    x1203 = x1124 * x8
    x1204 = x1065 - x1203
    x1205 = x1202 + x1204
    x1206 = x1205 * x6
    x1207 = x1201 - x1206
    x1208 = x1207 * x5
    x1209 = x1079 - x1112 * x8
    x1210 = x110 * x1131 + x1209
    x1211 = -x1200 * x6 + x1210 * x5
    x1212 = x1083 + x1198
    x1213 = x4 * (-x1197 + x1210 + x1212)
    x1214 = -x1207 * x6 + x1213
    x1215 = x1211 * x5 + x1214
    x1216 = x1205 * x5
    x1217 = x110 * x1124
    x1218 = x1143 * x8
    x1219 = x1087 - x1218
    x1220 = x1217 + x1219
    x1221 = x1220 * x6
    x1222 = x1216 - x1221
    x1223 = x1222 * x6
    x1224 = x1100 + x1203
    x1225 = x4 * (x1200 - x1202 + x1224)
    x1226 = -x1225
    x1227 = x1223 + x1226
    x1228 = x1208 - x1223 + x1225
    x1229 = x1106 + x1218
    x1230 = x4 * (x1205 - x1217 + x1229)
    x1231 = x1104 - x1155 * x8
    x1232 = x110 * x1143 + x1231
    x1233 = x110 * x1159
    x1234 = x1162 * x8
    x1235 = -x1234
    x1236 = x1233 + x1235
    x1237 = x1236 * x5
    x1238 = x110 * x1162
    x1239 = x1168 * x8
    x1240 = -x1239
    x1241 = x1238 + x1240
    x1242 = x1241 * x6
    x1243 = x1237 - x1242
    x1244 = x1243 * x5
    x1245 = -x1159 * x8
    x1246 = x110 * x1174 + x1245
    x1247 = -x1236 * x6 + x1246 * x5
    x1248 = x4 * (-x1233 + x1234 + x1246)
    x1249 = -x1243 * x6 + x1248
    x1250 = x1247 * x5 + x1249
    x1251 = x1241 * x5
    x1252 = x110 * x1168
    x1253 = x1184 * x8
    x1254 = -x1253
    x1255 = x1252 + x1254
    x1256 = x1255 * x6
    x1257 = x1251 - x1256
    x1258 = x1257 * x6
    x1259 = x4 * (x1236 - x1238 + x1239)
    x1260 = -x1259
    x1261 = x1258 + x1260
    x1262 = x1244 - x1258 + x1259
    x1263 = x4 * (x1241 - x1252 + x1253)
    x1264 = -x1195 * x8
    x1265 = x110 * x1184 + x1264
    x1266 = x1047 * x414
    x1267 = x1266 + x607
    x1268 = x1267 * x5
    x1269 = x1047 * x419
    x1270 = x1269 + x614
    x1271 = x1270 * x6
    x1272 = x1268 - x1271
    x1273 = x1272 * x5
    x1274 = x1047 * x425 + x601
    x1275 = -x1267 * x6 + x1274 * x5
    x1276 = x4 * (-x1266 + x1274 + x605)
    x1277 = -x1272 * x6 + x1276
    x1278 = x1275 * x5 + x1277
    x1279 = x1270 * x5
    x1280 = x1047 * x436
    x1281 = x1280 + x620
    x1282 = x1281 * x6
    x1283 = x1279 - x1282
    x1284 = x1283 * x6
    x1285 = x4 * (x1267 - x1269 + x624)
    x1286 = -x1285
    x1287 = x1284 + x1286
    x1288 = x1273 - x1284 + x1285
    x1289 = x4 * (x1270 - x1280 + x640)
    x1290 = x1047 * x456 + x642
    x1291 = 2.0 * x1132 - x1200 * x8
    x1292 = x110 * x1210 + x1291
    x1293 = x110 * x1205
    x1294 = x1220 * x8
    x1295 = 2.0 * x1137
    x1296 = -x1294 + x1295
    x1297 = x1293 + x1296
    x1298 = x1297 * x6
    x1299 = x110 * x1200
    x1300 = x1205 * x8
    x1301 = 2.0 * x1113
    x1302 = -x1300 + x1301
    x1303 = x1299 + x1302
    x1304 = -x1303 * x6
    x1305 = x1303 * x5
    x1306 = x1292 * x5 + x1304
    x1307 = x1300 - x1301
    x1308 = x4 * (x1292 - x1299 + x1307)
    x1309 = -x1298 + x1305
    x1310 = x1294 - x1295
    x1311 = x4 * (-x1293 + x1303 + x1310)
    x1312 = 2.0 * x1153 - x1232 * x8
    x1313 = x110 * x1220 + x1312
    x1314 = x1175 - x1236 * x8
    x1315 = x110 * x1246 + x1314
    x1316 = x110 * x1241
    x1317 = x1255 * x8
    x1318 = x1180 - x1317
    x1319 = x1316 + x1318
    x1320 = x1319 * x6
    x1321 = x110 * x1236
    x1322 = x1241 * x8
    x1323 = x1160 - x1322
    x1324 = x1321 + x1323
    x1325 = -x1324 * x6
    x1326 = x1324 * x5
    x1327 = x1315 * x5 + x1325
    x1328 = x1193 + x1322
    x1329 = x4 * (x1315 - x1321 + x1328)
    x1330 = -x1320 + x1326
    x1331 = x1196 + x1317
    x1332 = x4 * (-x1316 + x1324 + x1331)
    x1333 = x1194 - x1265 * x8
    x1334 = x110 * x1255 + x1333
    x1335 = -x1267 * x8
    x1336 = x110 * x1274 + x1335
    x1337 = x110 * x1270
    x1338 = x1281 * x8
    x1339 = -x1338
    x1340 = x1337 + x1339
    x1341 = x1340 * x6
    x1342 = x110 * x1267
    x1343 = x1270 * x8
    x1344 = -x1343
    x1345 = x1342 + x1344
    x1346 = -x1345 * x6
    x1347 = x1345 * x5
    x1348 = x1336 * x5 + x1346
    x1349 = x4 * (x1336 - x1342 + x1343)
    x1350 = -x1341 + x1347
    x1351 = x4 * (-x1337 + x1338 + x1345)
    x1352 = -x1290 * x8
    x1353 = x110 * x1281 + x1352
    x1354 = x1047 * x602 + x727
    x1355 = x1047 * x615
    x1356 = x1355 + x740
    x1357 = x1356 * x6
    x1358 = x1047 * x608
    x1359 = x1358 + x733
    x1360 = -x1359 * x6
    x1361 = x1359 * x5
    x1362 = x1354 * x5 + x1360
    x1363 = x4 * (x1354 - x1358 + x731)
    x1364 = -x1357 + x1361
    x1365 = x4 * (-x1355 + x1359 + x744)
    x1366 = x1047 * x621 + x746
    x1367 = x110 * x1303
    x1368 = 3.0 * x1213 - x1303 * x8
    x1369 = x110 * x1292 + x1368
    x1370 = x1297 * x8
    x1371 = 3.0 * x1225
    x1372 = x1370 - x1371
    x1373 = x4 * (-x1367 + x1369 + x1372)
    x1374 = -x1370 + x1371
    x1375 = x1367 + x1374
    x1376 = 3.0 * x1230 - x1313 * x8
    x1377 = x110 * x1297 + x1376
    x1378 = x110 * x1324
    x1379 = 2.0 * x1248 - x1324 * x8
    x1380 = x110 * x1315 + x1379
    x1381 = x1319 * x8
    x1382 = 2.0 * x1259
    x1383 = x1381 - x1382
    x1384 = x4 * (-x1378 + x1380 + x1383)
    x1385 = -x1381 + x1382
    x1386 = x1378 + x1385
    x1387 = 2.0 * x1263 - x1334 * x8
    x1388 = x110 * x1319 + x1387
    x1389 = x110 * x1345
    x1390 = x1276 - x1345 * x8
    x1391 = x110 * x1336 + x1390
    x1392 = x1340 * x8
    x1393 = x1286 + x1392
    x1394 = x4 * (-x1389 + x1391 + x1393)
    x1395 = x1285 - x1392
    x1396 = x1389 + x1395
    x1397 = x1289 - x1353 * x8
    x1398 = x110 * x1340 + x1397
    x1399 = x1356 * x8
    x1400 = x110 * x1359
    x1401 = -x1359 * x8
    x1402 = x110 * x1354 + x1401
    x1403 = x4 * (x1399 - x1400 + x1402)
    x1404 = -x1399
    x1405 = x1400 + x1404
    x1406 = -x1366 * x8
    x1407 = x110 * x1356 + x1406
    x1408 = x10 * x741
    x1409 = x1047 * x734
    x1410 = 4.0 * x625
    x1411 = -x10 * x734 + x1047 * x728 + 4.0 * x606
    x1412 = x4 * (x1408 - x1409 - x1410 + x1411)
    x1413 = -x1408 + x1409 + x1410
    x1414 = -x10 * x747 + x1047 * x741 + 4.0 * x641
    x1415 = x754 * x8
    x1416 = x748 * x751
    x1417 = x42 + x748 * x779 - x751 * x8
    x1418 = x4 * (x1415 - x1416 + x1417 + x47)
    x1419 = -x1415 + x1416 + x15
    x1420 = x1417 * x5 - x1419 * x6
    x1421 = x1419 * x5
    x1422 = x748 * x754
    x1423 = x760 * x8
    x1424 = x1422 - x1423 + x27
    x1425 = x1424 * x6
    x1426 = x1421 - x1425
    x1427 = x1418 + x1420 * x5 - x1426 * x6
    x1428 = x4 * (x1419 - x1422 + x1423 + x64)
    x1429 = x1426 * x5
    x1430 = x1424 * x5
    x1431 = x748 * x760
    x1432 = x770 * x8
    x1433 = x1431 - x1432 + x52
    x1434 = x1433 * x6
    x1435 = x1430 - x1434
    x1436 = x1435 * x6
    x1437 = x1428 + x1429 - x1436
    x1438 = x1427 * x5 - x1437 * x6 + x36 * (x1420 - x1421 + x1425)
    x1439 = x1437 * x5
    x1440 = x4 * (x1424 - x1431 + x1432 + x87)
    x1441 = x1435 * x5
    x1442 = x1433 * x5
    x1443 = x748 * x770
    x1444 = x792 * x8
    x1445 = x1443 - x1444 + x75
    x1446 = x1445 * x6
    x1447 = x1442 - x1446
    x1448 = x1447 * x6
    x1449 = x1440 + x1441 - x1448
    x1450 = x1449 * x6
    x1451 = x36 * (x1426 - x1430 + x1434)
    x1452 = x1439 - x1450 + x1451
    x1453 = -x1428
    x1454 = x4 * (x102 + x1433 - x1443 + x1444)
    x1455 = x100 + x748 * x792 - x8 * x806
    x1456 = -x1440
    x1457 = x8 * x814
    x1458 = x748 * x811
    x1459 = x136 + x748 * x826 + x780 - x8 * x811
    x1460 = x4 * (x1457 - x1458 + x1459 + x157 + x784)
    x1461 = x121 - x1457 + x1458 + x752
    x1462 = x1459 * x5 - x1461 * x6
    x1463 = x1461 * x5
    x1464 = x748 * x814
    x1465 = x8 * x820
    x1466 = x141 + x1464 - x1465 + x766
    x1467 = x1466 * x6
    x1468 = x1463 - x1467
    x1469 = x1460 + x1462 * x5 - x1468 * x6
    x1470 = x4 * (x1461 - x1464 + x1465 + x179 + x801)
    x1471 = x1468 * x5
    x1472 = x1466 * x5
    x1473 = x748 * x820
    x1474 = x8 * x836
    x1475 = x1473 - x1474 + x163 + x788
    x1476 = x1475 * x6
    x1477 = x1472 - x1476
    x1478 = x1477 * x6
    x1479 = x1470 + x1471 - x1478
    x1480 = x4 * (x1466 - x1473 + x1474 + x190 + x807)
    x1481 = x187 + x748 * x836 - x8 * x846 + x805
    x1482 = x8 * x852
    x1483 = x748 * x849
    x1484 = x223 + x748 * x864 - x8 * x849
    x1485 = x4 * (x1482 - x1483 + x1484 + x244)
    x1486 = -x1482 + x1483 + x208
    x1487 = x1484 * x5 - x1486 * x6
    x1488 = x1486 * x5
    x1489 = x748 * x852
    x1490 = x8 * x858
    x1491 = x1489 - x1490 + x228
    x1492 = x1491 * x6
    x1493 = x1488 - x1492
    x1494 = x1485 + x1487 * x5 - x1493 * x6
    x1495 = x4 * (x1486 - x1489 + x1490 + x266)
    x1496 = x1493 * x5
    x1497 = x1491 * x5
    x1498 = x748 * x858
    x1499 = x8 * x874
    x1500 = x1498 - x1499 + x250
    x1501 = x1500 * x6
    x1502 = x1497 - x1501
    x1503 = x1502 * x6
    x1504 = x1495 + x1496 - x1503
    x1505 = x4 * (x1491 - x1498 + x1499 + x277)
    x1506 = x274 + x748 * x874 - x8 * x885
    x1507 = -x1495
    x1508 = x8 * x891
    x1509 = x748 * x888
    x1510 = 2.0 * x812
    x1511 = -x1510
    x1512 = 2.0 * x827
    x1513 = x1512 + x299 + x748 * x895 - x8 * x888
    x1514 = x4 * (x1508 - x1509 + x1511 + x1513 + x312)
    x1515 = -x1508 + x1509 + x1510 + x311
    x1516 = x1513 * x5 - x1515 * x6
    x1517 = x1515 * x5
    x1518 = x748 * x891
    x1519 = x8 * x902
    x1520 = 2.0 * x832
    x1521 = x1518 - x1519 + x1520 + x321
    x1522 = x1521 * x6
    x1523 = x1517 - x1522
    x1524 = -x1520
    x1525 = x4 * (x1515 - x1518 + x1519 + x1524 + x338)
    x1526 = 2.0 * x845
    x1527 = x1526 + x341 + x748 * x902 - x8 * x910
    x1528 = x8 * x916
    x1529 = x748 * x913
    x1530 = x365 + x748 * x920 - x8 * x913 + x865
    x1531 = x4 * (x1528 - x1529 + x1530 + x377 + x883)
    x1532 = -x1528 + x1529 + x376 + x850
    x1533 = x1530 * x5 - x1532 * x6
    x1534 = x1532 * x5
    x1535 = x748 * x916
    x1536 = x8 * x927
    x1537 = x1535 - x1536 + x385 + x870
    x1538 = x1537 * x6
    x1539 = x1534 - x1538
    x1540 = x4 * (x1532 - x1535 + x1536 + x402 + x886)
    x1541 = x404 + x748 * x927 - x8 * x935 + x884
    x1542 = x8 * x941
    x1543 = x748 * x938
    x1544 = x429 + x748 * x945 - x8 * x938
    x1545 = x4 * (x1542 - x1543 + x1544 + x442)
    x1546 = -x1542 + x1543 + x441
    x1547 = x1544 * x5 - x1546 * x6
    x1548 = x1546 * x5
    x1549 = x748 * x941
    x1550 = x8 * x952
    x1551 = x1549 - x1550 + x451
    x1552 = x1551 * x6
    x1553 = x1548 - x1552
    x1554 = x4 * (x1546 - x1549 + x1550 + x468)
    x1555 = x471 + x748 * x952 - x8 * x961
    x1556 = x8 * x964
    x1557 = x748 * x967
    x1558 = 3.0 * x906
    x1559 = x483 + x748 * x962 - x8 * x967 + 3.0 * x897
    x1560 = x4 * (x1556 - x1557 - x1558 + x1559 + x503)
    x1561 = -x1556 + x1557 + x1558 + x502
    x1562 = x518 + x748 * x964 - x8 * x974 + 3.0 * x909
    x1563 = x8 * x977
    x1564 = x748 * x980
    x1565 = 2.0 * x931
    x1566 = -x1565
    x1567 = 2.0 * x922
    x1568 = x1567 + x526 + x748 * x975 - x8 * x980
    x1569 = x4 * (x1563 - x1564 + x1566 + x1568 + x544)
    x1570 = -x1563 + x1564 + x1565 + x543
    x1571 = 2.0 * x934
    x1572 = x1571 + x559 + x748 * x977 - x8 * x987
    x1573 = x8 * x990
    x1574 = x748 * x993
    x1575 = x566 + x748 * x988 - x8 * x993 + x947
    x1576 = x4 * (x1573 - x1574 + x1575 + x583 + x957)
    x1577 = -x1573 + x1574 + x582 + x956
    x1578 = -x1000 * x8 + x597 + x748 * x990 + x960
    x1579 = x1003 * x8
    x1580 = x1006 * x748
    x1581 = x1001 * x748 - x1006 * x8 + x606
    x1582 = x4 * (x1579 - x1580 + x1581 + x626)
    x1583 = -x1579 + x1580 + x625
    x1584 = x1003 * x748 - x1013 * x8 + x641
    x1585 = x1017 * x748 - x1019 * x8 + x650 + 4.0 * x971
    x1586 = x1019 * x748 - x1020 * x8 + x663 + 4.0 * x973
    x1587 = x1024 * x748 - x1026 * x8 + x672 + 3.0 * x984
    x1588 = x1026 * x748 - x1027 * x8 + x685 + 3.0 * x986
    x1589 = 2.0 * x997
    x1590 = x1031 * x748 - x1033 * x8 + x1589 + x693
    x1591 = 2.0 * x999
    x1592 = x1033 * x748 - x1034 * x8 + x1591 + x705
    x1593 = x1010 + x1037 * x748 - x1039 * x8 + x712
    x1594 = x1012 + x1039 * x748 - x1040 * x8 + x723
    x1595 = x1043 * x748 - x1045 * x8 + x732
    x1596 = x1045 * x748 - x1046 * x8 + x745
    x1597 = x1050 * x748
    x1598 = x1078 * x748 + x1130
    x1599 = x4 * (x1110 - x1597 + x1598)
    x1600 = x1111 + x1597
    x1601 = x1598 * x5 - x1600 * x6
    x1602 = x1600 * x5
    x1603 = x1053 * x748
    x1604 = x1115 + x1603
    x1605 = x1604 * x6
    x1606 = x1602 - x1605
    x1607 = x1599 + x1601 * x5 - x1606 * x6
    x1608 = x4 * (x1107 + x1600 - x1603)
    x1609 = x1606 * x5
    x1610 = x1604 * x5
    x1611 = x1059 * x748
    x1612 = x1123 + x1611
    x1613 = x1612 * x6
    x1614 = x1610 - x1613
    x1615 = x1614 * x6
    x1616 = x1608 + x1609 - x1615
    x1617 = x1607 * x5 - x1616 * x6 + x36 * (x1601 - x1602 + x1605)
    x1618 = x1616 * x5
    x1619 = x4 * (x1122 + x1604 - x1611)
    x1620 = x1614 * x5
    x1621 = x1612 * x5
    x1622 = x1069 * x748
    x1623 = x1142 + x1622
    x1624 = x1623 * x6
    x1625 = x1621 - x1624
    x1626 = x1625 * x6
    x1627 = x1619 + x1620 - x1626
    x1628 = x1627 * x6
    x1629 = x36 * (x1606 - x1610 + x1613)
    x1630 = x1618 - x1628 + x1629
    x1631 = -x1608
    x1632 = x4 * (x1141 + x1612 - x1622)
    x1633 = x1091 * x748 + x1154
    x1634 = -x1619
    x1635 = 35.3313566383285 * x106
    x1636 = x1112 * x748
    x1637 = x1131 * x748 + x1209
    x1638 = x4 * (x1212 - x1636 + x1637)
    x1639 = x1199 + x1636
    x1640 = x1637 * x5 - x1639 * x6
    x1641 = x1639 * x5
    x1642 = x1116 * x748
    x1643 = x1204 + x1642
    x1644 = x1643 * x6
    x1645 = x1641 - x1644
    x1646 = x1638 + x1640 * x5 - x1645 * x6
    x1647 = x4 * (x1224 + x1639 - x1642)
    x1648 = x1645 * x5
    x1649 = x1643 * x5
    x1650 = x1124 * x748
    x1651 = x1219 + x1650
    x1652 = x1651 * x6
    x1653 = x1649 - x1652
    x1654 = x1653 * x6
    x1655 = x1647 + x1648 - x1654
    x1656 = x4 * (x1229 + x1643 - x1650)
    x1657 = x1143 * x748 + x1231
    x1658 = 93.4779831475484 * x106
    x1659 = x1159 * x748
    x1660 = x1174 * x748 + x1245
    x1661 = x4 * (x1234 - x1659 + x1660)
    x1662 = x1235 + x1659
    x1663 = x1660 * x5 - x1662 * x6
    x1664 = x1662 * x5
    x1665 = x1162 * x748
    x1666 = x1240 + x1665
    x1667 = x1666 * x6
    x1668 = x1664 - x1667
    x1669 = x1661 + x1663 * x5 - x1668 * x6
    x1670 = x4 * (x1239 + x1662 - x1665)
    x1671 = x1668 * x5
    x1672 = x1666 * x5
    x1673 = x1168 * x748
    x1674 = x1254 + x1673
    x1675 = x1674 * x6
    x1676 = x1672 - x1675
    x1677 = x1676 * x6
    x1678 = x1670 + x1671 - x1677
    x1679 = x4 * (x1253 + x1666 - x1673)
    x1680 = x1184 * x748 + x1264
    x1681 = -x1670
    x1682 = x1200 * x748
    x1683 = x1210 * x748 + x1291
    x1684 = x4 * (x1307 - x1682 + x1683)
    x1685 = x1302 + x1682
    x1686 = x1683 * x5 - x1685 * x6
    x1687 = x1685 * x5
    x1688 = x1205 * x748
    x1689 = x1296 + x1688
    x1690 = x1689 * x6
    x1691 = x1687 - x1690
    x1692 = x4 * (x1310 + x1685 - x1688)
    x1693 = x1220 * x748 + x1312
    x1694 = 120.679557322504 * x106
    x1695 = x1236 * x748
    x1696 = x1246 * x748 + x1314
    x1697 = x4 * (x1328 - x1695 + x1696)
    x1698 = x1323 + x1695
    x1699 = x1696 * x5 - x1698 * x6
    x1700 = x1698 * x5
    x1701 = x1241 * x748
    x1702 = x1318 + x1701
    x1703 = x1702 * x6
    x1704 = x1700 - x1703
    x1705 = x4 * (x1331 + x1698 - x1701)
    x1706 = x1255 * x748 + x1333
    x1707 = 209.023124717498 * x106
    x1708 = x1267 * x748
    x1709 = x1274 * x748 + x1335
    x1710 = x4 * (x1343 - x1708 + x1709)
    x1711 = x1344 + x1708
    x1712 = x1709 * x5 - x1711 * x6
    x1713 = x1711 * x5
    x1714 = x1270 * x748
    x1715 = x1339 + x1714
    x1716 = x1715 * x6
    x1717 = x1713 - x1716
    x1718 = x4 * (x1338 + x1711 - x1714)
    x1719 = x1281 * x748 + x1352
    x1720 = x1303 * x748
    x1721 = x1292 * x748 + x1368
    x1722 = x4 * (x1372 - x1720 + x1721)
    x1723 = x1374 + x1720
    x1724 = x1297 * x748 + x1376
    x1725 = x1324 * x748
    x1726 = x1315 * x748 + x1379
    x1727 = x4 * (x1383 - x1725 + x1726)
    x1728 = x1385 + x1725
    x1729 = x1319 * x748 + x1387
    x1730 = x1345 * x748
    x1731 = x1336 * x748 + x1390
    x1732 = x4 * (x1393 - x1730 + x1731)
    x1733 = x1395 + x1730
    x1734 = x1340 * x748 + x1397
    x1735 = x1359 * x748
    x1736 = x1354 * x748 + x1401
    x1737 = x4 * (x1399 - x1735 + x1736)
    x1738 = x1404 + x1735
    x1739 = x1356 * x748 + x1406
    x1740 = 4.0 * x1308 + x1369 * x748 - x1375 * x8
    x1741 = 4.0 * x1311 + x1375 * x748 - x1377 * x8
    x1742 = 3.0 * x1329 + x1380 * x748 - x1386 * x8
    x1743 = 3.0 * x1332 + x1386 * x748 - x1388 * x8
    x1744 = 2.0 * x1349 + x1391 * x748 - x1396 * x8
    x1745 = 2.0 * x1351 + x1396 * x748 - x1398 * x8
    x1746 = x1363 + x1402 * x748 - x1405 * x8
    x1747 = x1365 + x1405 * x748 - x1407 * x8
    x1748 = x1411 * x748 - x1413 * x8
    x1749 = x1413 * x748 - x1414 * x8
    x1750 = x10 * x1053
    x1751 = x1047 * x1050
    x1752 = -x10 * x1050 + x1047 * x1078 + x42
    x1753 = x4 * (x1750 - x1751 + x1752 + x47)
    x1754 = x15 - x1750 + x1751
    x1755 = x1752 * x5 - x1754 * x6
    x1756 = x1754 * x5
    x1757 = x1047 * x1053
    x1758 = x10 * x1059
    x1759 = x1757 - x1758 + x27
    x1760 = x1759 * x6
    x1761 = x1756 - x1760
    x1762 = x1753 + x1755 * x5 - x1761 * x6
    x1763 = x4 * (x1754 - x1757 + x1758 + x64)
    x1764 = x1761 * x5
    x1765 = x1759 * x5
    x1766 = x1047 * x1059
    x1767 = x10 * x1069
    x1768 = x1766 - x1767 + x52
    x1769 = x1768 * x6
    x1770 = x1765 - x1769
    x1771 = x1770 * x6
    x1772 = x1763 + x1764 - x1771
    x1773 = x1762 * x5 - x1772 * x6 + x36 * (x1755 - x1756 + x1760)
    x1774 = x1772 * x5
    x1775 = x4 * (x1759 - x1766 + x1767 + x87)
    x1776 = x1770 * x5
    x1777 = x1768 * x5
    x1778 = x1047 * x1069
    x1779 = x10 * x1091
    x1780 = x1778 - x1779 + x75
    x1781 = x1780 * x6
    x1782 = x1777 - x1781
    x1783 = x1782 * x6
    x1784 = x1775 + x1776 - x1783
    x1785 = x1784 * x6
    x1786 = x36 * (x1761 - x1765 + x1769)
    x1787 = x1774 - x1785 + x1786
    x1788 = -x1763
    x1789 = x4 * (x102 + x1768 - x1778 + x1779)
    x1790 = -x10 * x1105 + x100 + x1047 * x1091
    x1791 = -x1775
    x1792 = x1759 * x8
    x1793 = x110 * x1754
    x1794 = -x1754 * x8
    x1795 = x110 * x1752 + x1794
    x1796 = x4 * (x1792 - x1793 + x1795)
    x1797 = -x1792
    x1798 = x1793 + x1797
    x1799 = x1795 * x5 - x1798 * x6
    x1800 = x1798 * x5
    x1801 = x110 * x1759
    x1802 = x1768 * x8
    x1803 = -x1802
    x1804 = x1801 + x1803
    x1805 = x1804 * x6
    x1806 = x1800 - x1805
    x1807 = x1796 + x1799 * x5 - x1806 * x6
    x1808 = x4 * (x1798 - x1801 + x1802)
    x1809 = x1806 * x5
    x1810 = x1804 * x5
    x1811 = x110 * x1768
    x1812 = x1780 * x8
    x1813 = -x1812
    x1814 = x1811 + x1813
    x1815 = x1814 * x6
    x1816 = x1810 - x1815
    x1817 = x1816 * x6
    x1818 = x1808 + x1809 - x1817
    x1819 = x4 * (x1804 - x1811 + x1812)
    x1820 = -x1790 * x8
    x1821 = x110 * x1780 + x1820
    x1822 = x10 * x1162
    x1823 = x1047 * x1159
    x1824 = -x10 * x1159 + x1047 * x1174 + x1079 + x223
    x1825 = x4 * (x1083 + x1822 - x1823 + x1824 + x244)
    x1826 = x1051 - x1822 + x1823 + x208
    x1827 = x1824 * x5 - x1826 * x6
    x1828 = x1826 * x5
    x1829 = x1047 * x1162
    x1830 = x10 * x1168
    x1831 = x1065 + x1829 - x1830 + x228
    x1832 = x1831 * x6
    x1833 = x1828 - x1832
    x1834 = x1825 + x1827 * x5 - x1833 * x6
    x1835 = x4 * (x1100 + x1826 - x1829 + x1830 + x266)
    x1836 = x1833 * x5
    x1837 = x1831 * x5
    x1838 = x1047 * x1168
    x1839 = x10 * x1184
    x1840 = x1087 + x1838 - x1839 + x250
    x1841 = x1840 * x6
    x1842 = x1837 - x1841
    x1843 = x1842 * x6
    x1844 = x1835 + x1836 - x1843
    x1845 = x4 * (x1106 + x1831 - x1838 + x1839 + x277)
    x1846 = -x10 * x1195 + x1047 * x1184 + x1104 + x274
    x1847 = -x1835
    x1848 = x110 * x1798
    x1849 = x1753 - x1798 * x8
    x1850 = x110 * x1795 + x1849
    x1851 = x1804 * x8
    x1852 = x1788 + x1851
    x1853 = x4 * (-x1848 + x1850 + x1852)
    x1854 = x1763 - x1851
    x1855 = x1848 + x1854
    x1856 = x1850 * x5 - x1855 * x6
    x1857 = x1855 * x5
    x1858 = x110 * x1804
    x1859 = x1814 * x8
    x1860 = x1775 - x1859
    x1861 = x1858 + x1860
    x1862 = x1861 * x6
    x1863 = x1857 - x1862
    x1864 = x1791 + x1859
    x1865 = x4 * (x1855 - x1858 + x1864)
    x1866 = x1789 - x1821 * x8
    x1867 = x110 * x1814 + x1866
    x1868 = x1831 * x8
    x1869 = x110 * x1826
    x1870 = -x1826 * x8
    x1871 = x110 * x1824 + x1870
    x1872 = x4 * (x1868 - x1869 + x1871)
    x1873 = -x1868
    x1874 = x1869 + x1873
    x1875 = x1871 * x5 - x1874 * x6
    x1876 = x1874 * x5
    x1877 = x110 * x1831
    x1878 = x1840 * x8
    x1879 = -x1878
    x1880 = x1877 + x1879
    x1881 = x1880 * x6
    x1882 = x1876 - x1881
    x1883 = x4 * (x1874 - x1877 + x1878)
    x1884 = -x1846 * x8
    x1885 = x110 * x1840 + x1884
    x1886 = x10 * x1270
    x1887 = x1047 * x1267
    x1888 = 2.0 * x1160
    x1889 = -x1888
    x1890 = 2.0 * x1175
    x1891 = -x10 * x1267 + x1047 * x1274 + x1890 + x429
    x1892 = x4 * (x1886 - x1887 + x1889 + x1891 + x442)
    x1893 = -x1886 + x1887 + x1888 + x441
    x1894 = x1891 * x5 - x1893 * x6
    x1895 = x1893 * x5
    x1896 = x1047 * x1270
    x1897 = x10 * x1281
    x1898 = 2.0 * x1180
    x1899 = x1896 - x1897 + x1898 + x451
    x1900 = x1899 * x6
    x1901 = x1895 - x1900
    x1902 = -x1898
    x1903 = x4 * (x1893 - x1896 + x1897 + x1902 + x468)
    x1904 = 2.0 * x1194
    x1905 = -x10 * x1290 + x1047 * x1281 + x1904 + x471
    x1906 = x110 * x1855
    x1907 = 2.0 * x1796 - x1855 * x8
    x1908 = x110 * x1850 + x1907
    x1909 = x1861 * x8
    x1910 = 2.0 * x1808
    x1911 = x1909 - x1910
    x1912 = x4 * (-x1906 + x1908 + x1911)
    x1913 = -x1909 + x1910
    x1914 = x1906 + x1913
    x1915 = 2.0 * x1819 - x1867 * x8
    x1916 = x110 * x1861 + x1915
    x1917 = x110 * x1874
    x1918 = x1825 - x1874 * x8
    x1919 = x110 * x1871 + x1918
    x1920 = x1880 * x8
    x1921 = x1847 + x1920
    x1922 = x4 * (-x1917 + x1919 + x1921)
    x1923 = x1835 - x1920
    x1924 = x1917 + x1923
    x1925 = x1845 - x1885 * x8
    x1926 = x110 * x1880 + x1925
    x1927 = x1899 * x8
    x1928 = x110 * x1893
    x1929 = -x1893 * x8
    x1930 = x110 * x1891 + x1929
    x1931 = x4 * (x1927 - x1928 + x1930)
    x1932 = -x1927
    x1933 = x1928 + x1932
    x1934 = -x1905 * x8
    x1935 = x110 * x1899 + x1934
    x1936 = x10 * x1356
    x1937 = x1047 * x1359
    x1938 = 3.0 * x1285
    x1939 = -x10 * x1359 + x1047 * x1354 + 3.0 * x1276 + x606
    x1940 = x4 * (x1936 - x1937 - x1938 + x1939 + x626)
    x1941 = -x1936 + x1937 + x1938 + x625
    x1942 = -x10 * x1366 + x1047 * x1356 + 3.0 * x1289 + x641
    x1943 = 3.0 * x1853 - x1914 * x8
    x1944 = x110 * x1908 + x1943
    x1945 = 3.0 * x1865 - x1916 * x8
    x1946 = x110 * x1914 + x1945
    x1947 = 2.0 * x1872 - x1924 * x8
    x1948 = x110 * x1919 + x1947
    x1949 = 2.0 * x1883 - x1926 * x8
    x1950 = x110 * x1924 + x1949
    x1951 = x1892 - x1933 * x8
    x1952 = x110 * x1930 + x1951
    x1953 = x1903 - x1935 * x8
    x1954 = x110 * x1933 + x1953
    x1955 = -x1941 * x8
    x1956 = x110 * x1939 + x1955
    x1957 = -x1942 * x8
    x1958 = x110 * x1941 + x1957
    x1959 = -x10 * x1413 + x1047 * x1411 + 4.0 * x1363 + x732
    x1960 = -x10 * x1414 + x1047 * x1413 + 4.0 * x1365 + x745
    x1961 = x1424 * x8
    x1962 = x1419 * x748
    x1963 = 2.0 * x752
    x1964 = x1417 * x748 - x1419 * x8 + 2.0 * x780
    x1965 = -x1961 + x1962 + x1963
    x1966 = x1964 * x5 - x1965 * x6
    x1967 = x1965 * x5
    x1968 = x1424 * x748
    x1969 = x1433 * x8
    x1970 = 2.0 * x766
    x1971 = x1968 - x1969 + x1970
    x1972 = x1971 * x6
    x1973 = x1967 - x1972
    x1974 = x1966 * x5 - x1973 * x6 + x4 * (x1961 - x1962 - x1963 + x1964)
    x1975 = x4 * (x1965 - x1968 + x1969 - x1970)
    x1976 = x1973 * x5
    x1977 = x1971 * x5
    x1978 = x1433 * x748
    x1979 = x1445 * x8
    x1980 = 2.0 * x788
    x1981 = x1978 - x1979 + x1980
    x1982 = x1981 * x6
    x1983 = x1977 - x1982
    x1984 = x1983 * x6
    x1985 = x1975 + x1976 - x1984
    x1986 = x1466 * x8
    x1987 = x1461 * x748
    x1988 = x1418 + x1459 * x748 - x1461 * x8 + x1512
    x1989 = x1428 + x1510 - x1986 + x1987
    x1990 = x1988 * x5 - x1989 * x6
    x1991 = x1989 * x5
    x1992 = x1466 * x748
    x1993 = x1475 * x8
    x1994 = x1440 + x1520 + x1992 - x1993
    x1995 = x1994 * x6
    x1996 = x1991 - x1995
    x1997 = x1491 * x8
    x1998 = x1486 * x748
    x1999 = 2.0 * x850
    x2000 = x1484 * x748 - x1486 * x8 + 2.0 * x865
    x2001 = -x1997 + x1998 + x1999
    x2002 = x2000 * x5 - x2001 * x6
    x2003 = x2001 * x5
    x2004 = x1491 * x748
    x2005 = x1500 * x8
    x2006 = 2.0 * x870
    x2007 = x2004 - x2005 + x2006
    x2008 = x2007 * x6
    x2009 = x2003 - x2008
    x2010 = x1521 * x8
    x2011 = x1515 * x748
    x2012 = 2.0 * x906
    x2013 = 2.0 * x1470
    x2014 = 2.0 * x1460 + x1513 * x748 - x1515 * x8 + 2.0 * x897
    x2015 = -x2010 + x2011 + x2012 + x2013
    x2016 = x1537 * x8
    x2017 = x1532 * x748
    x2018 = x1485 + x1530 * x748 - x1532 * x8 + x1567
    x2019 = x1495 + x1565 - x2016 + x2017
    x2020 = x1551 * x8
    x2021 = x1546 * x748
    x2022 = 2.0 * x956
    x2023 = x1544 * x748 - x1546 * x8 + 2.0 * x947
    x2024 = -x2020 + x2021 + x2022
    x2025 = x1604 * x8
    x2026 = x1600 * x748
    x2027 = x1079 + x1598 * x748 - x1600 * x8
    x2028 = x1051 - x2025 + x2026
    x2029 = x2027 * x5 - x2028 * x6
    x2030 = x2028 * x5
    x2031 = x1604 * x748
    x2032 = x1612 * x8
    x2033 = x1065 + x2031 - x2032
    x2034 = x2033 * x6
    x2035 = x2030 - x2034
    x2036 = x2029 * x5 - x2035 * x6 + x4 * (x1083 + x2025 - x2026 + x2027)
    x2037 = x4 * (x1100 + x2028 - x2031 + x2032)
    x2038 = x2035 * x5
    x2039 = x2033 * x5
    x2040 = x1612 * x748
    x2041 = x1623 * x8
    x2042 = x1087 + x2040 - x2041
    x2043 = x2042 * x6
    x2044 = x2039 - x2043
    x2045 = x2044 * x6
    x2046 = x2037 + x2038 - x2045
    x2047 = x1643 * x8
    x2048 = x1639 * x748
    x2049 = x1132 + x1599 + x1637 * x748 - x1639 * x8
    x2050 = x1113 + x1608 - x2047 + x2048
    x2051 = x2049 * x5 - x2050 * x6
    x2052 = x2050 * x5
    x2053 = x1643 * x748
    x2054 = x1651 * x8
    x2055 = x1137 + x1619 + x2053 - x2054
    x2056 = x2055 * x6
    x2057 = x2052 - x2056
    x2058 = x1666 * x8
    x2059 = x1662 * x748
    x2060 = x1175 + x1660 * x748 - x1662 * x8
    x2061 = x1160 - x2058 + x2059
    x2062 = x2060 * x5 - x2061 * x6
    x2063 = x2061 * x5
    x2064 = x1666 * x748
    x2065 = x1674 * x8
    x2066 = x1180 + x2064 - x2065
    x2067 = x2066 * x6
    x2068 = x2063 - x2067
    x2069 = x1689 * x8
    x2070 = x1685 * x748
    x2071 = 2.0 * x1647
    x2072 = x1213 + 2.0 * x1638 + x1683 * x748 - x1685 * x8
    x2073 = x1225 - x2069 + x2070 + x2071
    x2074 = x1702 * x8
    x2075 = x1698 * x748
    x2076 = x1248 + x1661 + x1696 * x748 - x1698 * x8
    x2077 = x1259 + x1670 - x2074 + x2075
    x2078 = x1715 * x8
    x2079 = x1711 * x748
    x2080 = x1276 + x1709 * x748 - x1711 * x8
    x2081 = x1285 - x2078 + x2079
    x2082 = x1754 * x748
    x2083 = x1752 * x748 + x1794
    x2084 = x1797 + x2082
    x2085 = x2083 * x5 - x2084 * x6
    x2086 = x2084 * x5
    x2087 = x1759 * x748
    x2088 = x1803 + x2087
    x2089 = x2088 * x6
    x2090 = x2086 - x2089
    x2091 = x2085 * x5 - x2090 * x6 + x4 * (x1792 - x2082 + x2083)
    x2092 = x4 * (x1802 + x2084 - x2087)
    x2093 = x2090 * x5
    x2094 = x2088 * x5
    x2095 = x1768 * x748
    x2096 = x1813 + x2095
    x2097 = x2096 * x6
    x2098 = x2094 - x2097
    x2099 = x2098 * x6
    x2100 = x2092 + x2093 - x2099
    x2101 = x1798 * x748
    x2102 = x1795 * x748 + x1849
    x2103 = x1854 + x2101
    x2104 = x2102 * x5 - x2103 * x6
    x2105 = x2103 * x5
    x2106 = x1804 * x748
    x2107 = x1860 + x2106
    x2108 = x2107 * x6
    x2109 = x2105 - x2108
    x2110 = x1826 * x748
    x2111 = x1824 * x748 + x1870
    x2112 = x1873 + x2110
    x2113 = x2111 * x5 - x2112 * x6
    x2114 = x2112 * x5
    x2115 = x1831 * x748
    x2116 = x1879 + x2115
    x2117 = x2116 * x6
    x2118 = x2114 - x2117
    x2119 = x1855 * x748
    x2120 = x1850 * x748 + x1907
    x2121 = x1913 + x2119
    x2122 = x1874 * x748
    x2123 = x1871 * x748 + x1918
    x2124 = x1923 + x2122
    x2125 = x1893 * x748
    x2126 = x1891 * x748 + x1929
    x2127 = x1932 + x2125
    x2128 = x10 * x1759
    x2129 = x1047 * x1754
    x2130 = 2.0 * x1051
    x2131 = -x10 * x1754 + x1047 * x1752 + 2.0 * x1079
    x2132 = x4 * (x2128 - x2129 - x2130 + x2131)
    x2133 = -x2128 + x2129 + x2130
    x2134 = x2131 * x5 - x2133 * x6
    x2135 = x2133 * x5
    x2136 = x1047 * x1759
    x2137 = x10 * x1768
    x2138 = 2.0 * x1065
    x2139 = x2136 - x2137 + x2138
    x2140 = x2139 * x6
    x2141 = x2135 - x2140
    x2142 = x2132 + x2134 * x5 - x2141 * x6
    x2143 = x4 * (x2133 - x2136 + x2137 - x2138)
    x2144 = x2141 * x5
    x2145 = x2139 * x5
    x2146 = x1047 * x1768
    x2147 = x10 * x1780
    x2148 = 2.0 * x1087
    x2149 = x2146 - x2147 + x2148
    x2150 = x2149 * x6
    x2151 = x2145 - x2150
    x2152 = x2151 * x6
    x2153 = x2143 + x2144 - x2152
    x2154 = x4 * (x2139 - x2146 + x2147 - x2148)
    x2155 = -x10 * x1790 + x1047 * x1780 + 2.0 * x1104
    x2156 = -x2143
    x2157 = x2139 * x8
    x2158 = x110 * x2133
    x2159 = x110 * x2131 - x2133 * x8
    x2160 = x4 * (x2157 - x2158 + x2159)
    x2161 = -x2157 + x2158
    x2162 = x2159 * x5 - x2161 * x6
    x2163 = x2161 * x5
    x2164 = x110 * x2139
    x2165 = x2149 * x8
    x2166 = x2164 - x2165
    x2167 = x2166 * x6
    x2168 = x2163 - x2167
    x2169 = x4 * (x2161 - x2164 + x2165)
    x2170 = x110 * x2149 - x2155 * x8
    x2171 = x10 * x1831
    x2172 = x1047 * x1826
    x2173 = -x10 * x1826 + x1047 * x1824 + x1753 + x1890
    x2174 = x4 * (x1788 + x1889 + x2171 - x2172 + x2173)
    x2175 = x1763 + x1888 - x2171 + x2172
    x2176 = x2173 * x5 - x2175 * x6
    x2177 = x2175 * x5
    x2178 = x1047 * x1831
    x2179 = x10 * x1840
    x2180 = x1775 + x1898 + x2178 - x2179
    x2181 = x2180 * x6
    x2182 = x2177 - x2181
    x2183 = x4 * (x1791 + x1902 + x2175 - x2178 + x2179)
    x2184 = -x10 * x1846 + x1047 * x1840 + x1789 + x1904
    x2185 = x2166 * x8
    x2186 = x110 * x2161
    x2187 = x110 * x2159 + x2132 - x2161 * x8
    x2188 = x4 * (x2156 + x2185 - x2186 + x2187)
    x2189 = x2143 - x2185 + x2186
    x2190 = x110 * x2166 + x2154 - x2170 * x8
    x2191 = x2180 * x8
    x2192 = x110 * x2175
    x2193 = x110 * x2173 - x2175 * x8
    x2194 = x4 * (x2191 - x2192 + x2193)
    x2195 = -x2191 + x2192
    x2196 = x110 * x2180 - x2184 * x8
    x2197 = x10 * x1899
    x2198 = x1047 * x1893
    x2199 = 2.0 * x1285
    x2200 = 2.0 * x1835
    x2201 = -x10 * x1893 + x1047 * x1891 + 2.0 * x1276 + 2.0 * x1825
    x2202 = x4 * (x2197 - x2198 - x2199 - x2200 + x2201)
    x2203 = -x2197 + x2198 + x2199 + x2200
    x2204 = -x10 * x1905 + x1047 * x1899 + 2.0 * x1289 + 2.0 * x1845
    x2205 = x110 * x2187 + 2.0 * x2160 - x2189 * x8
    x2206 = x110 * x2189 + 2.0 * x2169 - x2190 * x8
    x2207 = x110 * x2193 + x2174 - x2195 * x8
    x2208 = x110 * x2195 + x2183 - x2196 * x8
    x2209 = x110 * x2201 - x2203 * x8
    x2210 = x110 * x2203 - x2204 * x8
    x2211 = -x10 * x1941 + x1047 * x1939 + 2.0 * x1363 + 3.0 * x1892
    x2212 = -x10 * x1942 + x1047 * x1941 + 2.0 * x1365 + 3.0 * x1903

    # 150 item(s)
    result[0, 0] = numpy.sum(
        x107
        * (
            x3 * (x3 * x71 + x4 * (-x40 + x50 + x66) - x6 * x95 + x98 * (x66 - x96 + x97))
            + x36 * (x71 - x72 + x91 - x94)
            - x6
            * (
                x3 * x95
                + x4 * (x68 - x73 + x99)
                - x6
                * (
                    x3 * x90
                    - x6
                    * (
                        x48 * (x102 + x59 - x76 + x81)
                        + x5 * x85
                        - x6
                        * (
                            x36 * (x57 + x78 - x79)
                            + x5 * x82
                            - x6 * (x100 + x5 * x80 + x6 * (x101 * x6 - x5 * x77))
                        )
                    )
                    + x98 * (x62 - x74 + x83 - x84)
                )
                + x98 * (-x103 + x104 + x99)
            )
            - x98
            * (
                x104 * x3
                + x104 * x6
                - x3 * x97
                - x48 * (-x105 + x3 * x44 + x45 + x69)
                + x48 * (x105 - x3 * x34 + x38 + x92)
                - x6 * (x103 + x89)
                - x70
                + x93
            )
        )
    )
    result[0, 1] = numpy.sum(
        x194
        * (
            x3
            * (
                x160 * x3
                - x183 * x6
                + x4 * (-x132 + x140 + x154)
                + x48 * (x154 - x184 + x185)
            )
            + x36 * (x160 - x161 + x178 - x182)
            + x48
            * (
                x159
                - x181
                + x185 * x3
                - x192 * x3
                - x192 * x6
                - x36 * (-x128 * x3 + x130 + x180 + x193)
                + x36 * (x135 * x3 + x137 + x158 - x193)
                + x6 * (x176 + x191)
            )
            - x6
            * (
                x183 * x3
                + x4 * (x156 - x162 + x186)
                + x48 * (x186 - x191 + x192)
                - x6
                * (
                    x177 * x3
                    + x48 * (x151 - x164 + x172 + x190)
                    - x6
                    * (
                        x173 * x5
                        + x36 * (x149 - x165 + x170)
                        - x6 * (x171 * x5 + x187 - x6 * (x169 * x5 - x189 * x6))
                    )
                )
            )
        )
    )
    result[0, 2] = numpy.sum(
        x194
        * (
            x3
            * (
                x247 * x3
                - x270 * x6
                + x4 * (-x219 + x227 + x241)
                + x48 * (x241 - x271 + x272)
            )
            + x36 * (x247 - x248 + x265 - x269)
            + x48
            * (
                x246
                - x268
                + x272 * x3
                - x279 * x3
                - x279 * x6
                - x36 * (-x215 * x3 + x217 + x267 + x280)
                + x36 * (x222 * x3 + x224 + x245 - x280)
                + x6 * (x263 + x278)
            )
            - x6
            * (
                x270 * x3
                + x4 * (x243 - x249 + x273)
                + x48 * (x273 - x278 + x279)
                - x6
                * (
                    x264 * x3
                    + x48 * (x238 - x251 + x259 + x277)
                    - x6
                    * (
                        x260 * x5
                        + x36 * (x236 - x252 + x257)
                        - x6 * (x258 * x5 + x274 - x6 * (x256 * x5 - x276 * x6))
                    )
                )
            )
        )
    )
    result[0, 3] = numpy.sum(
        x347
        * (
            x3
            * (
                x3 * x317
                - x335 * x6
                + x36 * (x313 - x336 + x337)
                + x4 * (-x293 + x301 + x313)
            )
            + x36 * (x317 - x318 + x332 - x334)
            + x36
            * (
                x3 * x337
                - x3 * x345
                + x316
                - x333
                - x345 * x6
                + x4 * (x290 + x295 * x3 + x296 - x346)
                - x4 * (-x289 * x3 + x291 + x307 + x346)
                + x6 * (x330 + x344)
            )
            - x6
            * (
                x3 * x335
                + x36 * (x339 - x344 + x345)
                + x4 * (x315 - x319 + x339)
                - x6
                * (
                    x3 * x331
                    + x36 * (x308 - x322 + x327)
                    - x6 * (x328 * x5 + x341 - x6 * (x326 * x5 - x343 * x6))
                )
            )
        )
    )
    result[0, 4] = numpy.sum(
        x410
        * (
            x3
            * (
                x3 * x382
                + x36 * (x378 - x400 + x401)
                - x399 * x6
                + x4 * (-x360 + x367 + x378)
            )
            + x36 * (x382 - x383 + x396 - x398)
            + x36
            * (
                x3 * x401
                - x3 * x408
                + x381
                - x397
                - x4 * (-x3 * x356 + x358 + x373 + x409)
                + x4 * (x3 * x362 + x357 + x363 - x409)
                - x408 * x6
                + x6 * (x394 + x407)
            )
            - x6
            * (
                x3 * x399
                + x36 * (x403 - x407 + x408)
                + x4 * (x380 - x384 + x403)
                - x6
                * (
                    x3 * x395
                    + x36 * (x374 - x386 + x391)
                    - x6 * (x392 * x5 + x404 - x6 * (x390 * x5 - x406 * x6))
                )
            )
        )
    )
    result[0, 5] = numpy.sum(
        x347
        * (
            x3
            * (
                x3 * x447
                + x36 * (x443 - x466 + x467)
                + x4 * (-x423 + x431 + x443)
                - x465 * x6
            )
            + x36 * (x447 - x448 + x462 - x464)
            + x36
            * (
                x3 * x467
                - x3 * x475
                - x4 * (-x3 * x419 + x421 + x437 + x476)
                + x4 * (x3 * x425 + x420 + x426 - x476)
                + x446
                - x463
                - x475 * x6
                + x6 * (x460 + x474)
            )
            - x6
            * (
                x3 * x465
                + x36 * (x469 - x474 + x475)
                + x4 * (x445 - x449 + x469)
                - x6
                * (
                    x3 * x461
                    + x36 * (x438 - x452 + x457)
                    - x6 * (x458 * x5 + x471 - x6 * (x456 * x5 - x473 * x6))
                )
            )
        )
    )
    result[0, 6] = numpy.sum(
        x194
        * (
            x3
            * (
                x3 * x511
                + x4 * (x487 - x504 + x505)
                + x4 * (x505 - x508 + x509)
                - x516 * x6
            )
            + x36 * (x503 + x511 - x512 + x515)
            + x4 * (x3 * x487 - x3 * x507 + x483 + x503 - x507 * x6 + x6 * (x493 + x500))
            - x6
            * (
                x3 * x516
                + x4 * (-x493 + x499 + x507)
                + x4 * (x499 + x510 - x513)
                - x6 * (x3 * x514 + x518 - x6 * (x498 * x5 - x520 * x6))
            )
        )
    )
    result[0, 7] = numpy.sum(
        x410
        * (
            x3
            * (
                x3 * x552
                + x4 * (x530 - x545 + x546)
                + x4 * (x546 - x549 + x550)
                - x557 * x6
            )
            + x36 * (x544 + x552 - x553 + x556)
            + x4 * (x3 * x530 - x3 * x548 + x526 + x544 - x548 * x6 + x6 * (x535 + x541))
            - x6
            * (
                x3 * x557
                + x4 * (-x535 + x540 + x548)
                + x4 * (x540 + x551 - x554)
                - x6 * (x3 * x555 + x559 - x6 * (x5 * x539 - x561 * x6))
            )
        )
    )
    result[0, 8] = numpy.sum(
        x410
        * (
            x3
            * (
                x3 * x591
                + x4 * (x570 - x584 + x585)
                + x4 * (x585 - x588 + x589)
                - x596 * x6
            )
            + x36 * (x583 + x591 - x592 + x595)
            + x4 * (x3 * x570 - x3 * x587 + x566 + x583 - x587 * x6 + x6 * (x575 + x581))
            - x6
            * (
                x3 * x596
                + x4 * (-x575 + x580 + x587)
                + x4 * (x580 + x590 - x593)
                - x6 * (x3 * x594 + x597 - x6 * (x5 * x579 - x599 * x6))
            )
        )
    )
    result[0, 9] = numpy.sum(
        x194
        * (
            x3
            * (
                x3 * x634
                + x4 * (x610 - x627 + x628)
                + x4 * (x628 - x631 + x632)
                - x6 * x639
            )
            + x36 * (x626 + x634 - x635 + x638)
            + x4 * (x3 * x610 - x3 * x630 - x6 * x630 + x6 * (x616 + x623) + x606 + x626)
            - x6
            * (
                x3 * x639
                + x4 * (-x616 + x622 + x630)
                + x4 * (x622 + x633 - x636)
                - x6 * (x3 * x637 - x6 * (x5 * x621 - x6 * x643) + x641)
            )
        )
    )
    result[0, 10] = numpy.sum(
        x107
        * (
            x3 * (x3 * x653 - x6 * x661 + x650)
            + x36 * (x653 - x654 + x660)
            - x6 * (x3 * x661 - x6 * (x3 * x659 - x6 * x665) + x663)
        )
    )
    result[0, 11] = numpy.sum(
        x194
        * (
            x3 * (x3 * x675 - x6 * x683 + x672)
            + x36 * (x675 - x676 + x682)
            - x6 * (x3 * x683 - x6 * (x3 * x681 - x6 * x687) + x685)
        )
    )
    result[0, 12] = numpy.sum(
        x347
        * (
            x3 * (x3 * x696 - x6 * x703 + x693)
            + x36 * (x696 - x697 + x702)
            - x6 * (x3 * x703 - x6 * (x3 * x701 - x6 * x707) + x705)
        )
    )
    result[0, 13] = numpy.sum(
        x194
        * (
            x3 * (x3 * x715 - x6 * x722 + x712)
            + x36 * (x715 - x716 + x721)
            - x6 * (x3 * x722 - x6 * (x3 * x720 - x6 * x725) + x723)
        )
    )
    result[0, 14] = numpy.sum(
        x107
        * (
            x3 * (x3 * x735 - x6 * x743 + x732)
            + x36 * (x735 - x736 + x742)
            - x6 * (x3 * x743 - x6 * (x3 * x741 - x6 * x747) + x745)
        )
    )
    result[1, 0] = numpy.sum(
        x808
        * (
            -x3 * (-x3 * x786 + x6 * x804 + x98 * (x765 - x775 + x776 - x783))
            + x4 * (-x778 + x786 + x803)
            - x6
            * (
                x3 * x804
                - x6
                * (
                    x48 * (x774 - x789 + x795 + x807)
                    + x5 * x799
                    - x6
                    * (
                        x36 * (x772 - x790 + x793)
                        + x5 * x796
                        - x6 * (x5 * x794 - x6 * (x5 * x792 - x6 * x806) + x805)
                    )
                )
                + x98 * (x777 - x787 + x797 - x798)
            )
            + x98 * (-x3 * x777 + x3 * x783 + x785 + x803)
        )
    )
    result[1, 1] = numpy.sum(
        x410
        * (
            -x3 * (-x3 * x831 + x48 * (x812 + x817 - x823 - x829) + x6 * x844)
            + x4 * (-x825 + x831 + x843)
            + x48 * (-x3 * x824 + x3 * x829 + x830 + x843)
            - x6
            * (
                x3 * x844
                + x48 * (x824 - x832 - x833 + x839)
                - x6
                * (
                    x36 * (x822 - x834 + x837)
                    + x5 * x840
                    - x6 * (x5 * x838 - x6 * (x5 * x836 - x6 * x846) + x845)
                )
            )
        )
    )
    result[1, 2] = numpy.sum(
        x410
        * (
            x3 * (x3 * x869 + x48 * (-x855 + x861 + x867 + x883) - x6 * x882)
            + x4 * (-x863 + x869 + x881)
            + x48 * (-x3 * x862 + x3 * x867 + x868 + x881)
            - x6
            * (
                x3 * x882
                + x48 * (x862 - x871 + x877 + x886)
                - x6
                * (
                    x36 * (x860 - x872 + x875)
                    + x5 * x878
                    - x6 * (x5 * x876 - x6 * (x5 * x874 - x6 * x885) + x884)
                )
            )
        )
    )
    result[1, 3] = numpy.sum(
        x911
        * (
            x3 * (x3 * x899 + x36 * (-x889 + x892 + x896) - x6 * x908)
            + x36 * (-x3 * x893 + x3 * x896 + x898 + x907)
            + x4 * (-x894 + x899 + x907)
            - x6
            * (
                x3 * x908
                + x36 * (x893 - x900 + x903)
                - x6 * (x5 * x904 - x6 * (x5 * x902 - x6 * x910) + x909)
            )
        )
    )
    result[1, 4] = numpy.sum(
        x936
        * (
            x3 * (x3 * x924 + x36 * (-x914 + x917 + x921) - x6 * x933)
            + x36 * (-x3 * x918 + x3 * x921 + x923 + x932)
            + x4 * (-x919 + x924 + x932)
            - x6
            * (
                x3 * x933
                + x36 * (x918 - x925 + x928)
                - x6 * (x5 * x929 - x6 * (x5 * x927 - x6 * x935) + x934)
            )
        )
    )
    result[1, 5] = numpy.sum(
        x911
        * (
            x3 * (x3 * x949 + x36 * (-x939 + x942 + x946) - x6 * x959)
            + x36 * (-x3 * x943 + x3 * x946 + x948 + x958)
            + x4 * (-x944 + x949 + x958)
            - x6
            * (
                x3 * x959
                + x36 * (x943 - x950 + x953)
                - x6 * (x5 * x954 - x6 * (x5 * x952 - x6 * x961) + x960)
            )
        )
    )
    result[1, 6] = numpy.sum(
        x410
        * (
            x3 * (x3 * x970 - x6 * x972 + x971)
            + x4 * (x965 - x969 + x970)
            + x4 * (x3 * x962 - x3 * x967 + x965 + x968)
            - x6 * (x3 * x972 - x6 * (x5 * x964 - x6 * x974) + x973)
        )
    )
    result[1, 7] = numpy.sum(
        x936
        * (
            x3 * (x3 * x983 - x6 * x985 + x984)
            + x4 * (x978 - x982 + x983)
            + x4 * (x3 * x975 - x3 * x980 + x978 + x981)
            - x6 * (x3 * x985 - x6 * (x5 * x977 - x6 * x987) + x986)
        )
    )
    result[1, 8] = numpy.sum(
        x936
        * (
            x3 * (x3 * x996 - x6 * x998 + x997)
            + x4 * (x991 - x995 + x996)
            + x4 * (x3 * x988 - x3 * x993 + x991 + x994)
            - x6 * (x3 * x998 + x6 * (x1000 * x6 - x5 * x990) + x999)
        )
    )
    result[1, 9] = numpy.sum(
        x410
        * (
            x3 * (x1009 * x3 + x1010 - x1011 * x6)
            + x4 * (x1004 - x1008 + x1009)
            + x4 * (x1001 * x3 + x1004 - x1006 * x3 + x1007)
            - x6 * (x1011 * x3 + x1012 - x6 * (x1003 * x5 - x1013 * x6))
        )
    )
    result[1, 10] = numpy.sum(
        x808 * (x1018 + x3 * (x1017 * x3 - x1019 * x6) - x6 * (x1019 * x3 - x1020 * x6))
    )
    result[1, 11] = numpy.sum(
        x410 * (x1025 + x3 * (x1024 * x3 - x1026 * x6) - x6 * (x1026 * x3 - x1027 * x6))
    )
    result[1, 12] = numpy.sum(
        x911 * (x1032 + x3 * (x1031 * x3 - x1033 * x6) - x6 * (x1033 * x3 - x1034 * x6))
    )
    result[1, 13] = numpy.sum(
        x410 * (x1038 + x3 * (x1037 * x3 - x1039 * x6) - x6 * (x1039 * x3 - x1040 * x6))
    )
    result[1, 14] = numpy.sum(
        x808 * (x1044 + x3 * (x1043 * x3 - x1045 * x6) - x6 * (x1045 * x3 - x1046 * x6))
    )
    result[2, 0] = numpy.sum(
        x808
        * (
            -x3 * (-x1085 * x3 + x1103 * x6 + x98 * (x1064 - x1074 + x1075 - x1082))
            + x4 * (-x1077 + x1085 + x1102)
            - x6
            * (
                x1103 * x3
                - x6
                * (
                    x1098 * x5
                    + x48 * (x1073 - x1088 + x1094 + x1106)
                    - x6
                    * (
                        x1095 * x5
                        + x36 * (x1071 - x1089 + x1092)
                        - x6 * (x1093 * x5 + x1104 - x6 * (x1091 * x5 - x1105 * x6))
                    )
                )
                + x98 * (x1076 - x1086 + x1096 - x1097)
            )
            + x98 * (-x1076 * x3 + x1082 * x3 + x1084 + x1102)
        )
    )
    result[2, 1] = numpy.sum(
        x410
        * (
            x3 * (x1136 * x3 - x1151 * x6 + x48 * (-x1119 + x1127 + x1134 + x1152))
            + x4 * (-x1129 + x1136 + x1150)
            + x48 * (-x1128 * x3 + x1134 * x3 + x1135 + x1150)
            - x6
            * (
                x1151 * x3
                + x48 * (x1128 - x1138 + x1146 + x1156)
                - x6
                * (
                    x1147 * x5
                    + x36 * (x1126 - x1139 + x1144)
                    - x6 * (x1145 * x5 + x1153 - x6 * (x1143 * x5 - x1155 * x6))
                )
            )
        )
    )
    result[2, 2] = numpy.sum(
        x410
        * (
            x3 * (x1179 * x3 - x1192 * x6 + x48 * (-x1165 + x1171 + x1177 + x1193))
            + x4 * (-x1173 + x1179 + x1191)
            + x48 * (-x1172 * x3 + x1177 * x3 + x1178 + x1191)
            - x6
            * (
                x1192 * x3
                + x48 * (x1172 - x1181 + x1187 + x1196)
                - x6
                * (
                    x1188 * x5
                    + x36 * (x1170 - x1182 + x1185)
                    - x6 * (x1186 * x5 + x1194 - x6 * (x1184 * x5 - x1195 * x6))
                )
            )
        )
    )
    result[2, 3] = numpy.sum(
        x911
        * (
            x3 * (x1215 * x3 - x1228 * x6 + x36 * (-x1201 + x1206 + x1211))
            + x36 * (-x1207 * x3 + x1211 * x3 + x1214 + x1227)
            + x4 * (-x1208 + x1215 + x1227)
            - x6
            * (
                x1228 * x3
                + x36 * (x1207 - x1216 + x1221)
                - x6 * (x1222 * x5 + x1230 - x6 * (x1220 * x5 - x1232 * x6))
            )
        )
    )
    result[2, 4] = numpy.sum(
        x936
        * (
            x3 * (x1250 * x3 - x1262 * x6 + x36 * (-x1237 + x1242 + x1247))
            + x36 * (-x1243 * x3 + x1247 * x3 + x1249 + x1261)
            + x4 * (-x1244 + x1250 + x1261)
            - x6
            * (
                x1262 * x3
                + x36 * (x1243 - x1251 + x1256)
                - x6 * (x1257 * x5 + x1263 - x6 * (x1255 * x5 - x1265 * x6))
            )
        )
    )
    result[2, 5] = numpy.sum(
        x911
        * (
            x3 * (x1278 * x3 - x1288 * x6 + x36 * (-x1268 + x1271 + x1275))
            + x36 * (-x1272 * x3 + x1275 * x3 + x1277 + x1287)
            + x4 * (-x1273 + x1278 + x1287)
            - x6
            * (
                x1288 * x3
                + x36 * (x1272 - x1279 + x1282)
                - x6 * (x1283 * x5 + x1289 - x6 * (x1281 * x5 - x1290 * x6))
            )
        )
    )
    result[2, 6] = numpy.sum(
        x410
        * (
            x3 * (x1306 * x3 + x1308 - x1309 * x6)
            + x4 * (x1298 - x1305 + x1306)
            + x4 * (x1292 * x3 + x1298 - x1303 * x3 + x1304)
            - x6 * (x1309 * x3 + x1311 - x6 * (x1297 * x5 - x1313 * x6))
        )
    )
    result[2, 7] = numpy.sum(
        x936
        * (
            x3 * (x1327 * x3 + x1329 - x1330 * x6)
            + x4 * (x1320 - x1326 + x1327)
            + x4 * (x1315 * x3 + x1320 - x1324 * x3 + x1325)
            - x6 * (x1330 * x3 + x1332 - x6 * (x1319 * x5 - x1334 * x6))
        )
    )
    result[2, 8] = numpy.sum(
        x936
        * (
            x3 * (x1348 * x3 + x1349 - x1350 * x6)
            + x4 * (x1341 - x1347 + x1348)
            + x4 * (x1336 * x3 + x1341 - x1345 * x3 + x1346)
            - x6 * (x1350 * x3 + x1351 - x6 * (x1340 * x5 - x1353 * x6))
        )
    )
    result[2, 9] = numpy.sum(
        x410
        * (
            x3 * (x1362 * x3 + x1363 - x1364 * x6)
            + x4 * (x1357 - x1361 + x1362)
            + x4 * (x1354 * x3 + x1357 - x1359 * x3 + x1360)
            - x6 * (x1364 * x3 + x1365 - x6 * (x1356 * x5 - x1366 * x6))
        )
    )
    result[2, 10] = numpy.sum(
        x808 * (x1373 + x3 * (x1369 * x3 - x1375 * x6) - x6 * (x1375 * x3 - x1377 * x6))
    )
    result[2, 11] = numpy.sum(
        x410 * (x1384 + x3 * (x1380 * x3 - x1386 * x6) - x6 * (x1386 * x3 - x1388 * x6))
    )
    result[2, 12] = numpy.sum(
        x911 * (x1394 + x3 * (x1391 * x3 - x1396 * x6) - x6 * (x1396 * x3 - x1398 * x6))
    )
    result[2, 13] = numpy.sum(
        x410 * (x1403 + x3 * (x1402 * x3 - x1405 * x6) - x6 * (x1405 * x3 - x1407 * x6))
    )
    result[2, 14] = numpy.sum(
        x808 * (x1412 + x3 * (x1411 * x3 - x1413 * x6) - x6 * (x1413 * x3 - x1414 * x6))
    )
    result[3, 0] = numpy.sum(
        x808
        * (
            x3 * (x1438 * x5 - x1452 * x6 + x48 * (x1427 - x1429 + x1436 + x1453))
            - x6
            * (
                x1452 * x5
                + x48 * (x1437 - x1441 + x1448 + x1456)
                - x6
                * (
                    x1449 * x5
                    + x36 * (x1435 - x1442 + x1446)
                    - x6 * (x1447 * x5 + x1454 - x6 * (x1445 * x5 - x1455 * x6))
                )
            )
            + x98 * (x1438 - x1439 + x1450 - x1451)
        )
    )
    result[3, 1] = numpy.sum(
        x410
        * (
            x3 * (x1469 * x5 - x1479 * x6 + x36 * (x1462 - x1463 + x1467))
            + x48 * (x1469 - x1470 - x1471 + x1478)
            - x6
            * (
                x1479 * x5
                + x36 * (x1468 - x1472 + x1476)
                - x6 * (x1477 * x5 + x1480 - x6 * (x1475 * x5 - x1481 * x6))
            )
        )
    )
    result[3, 2] = numpy.sum(
        x410
        * (
            x3 * (x1494 * x5 - x1504 * x6 + x36 * (x1487 - x1488 + x1492))
            + x48 * (x1494 - x1496 + x1503 + x1507)
            - x6
            * (
                x1504 * x5
                + x36 * (x1493 - x1497 + x1501)
                - x6 * (x1502 * x5 + x1505 - x6 * (x1500 * x5 - x1506 * x6))
            )
        )
    )
    result[3, 3] = numpy.sum(
        x911
        * (
            x3 * (x1514 + x1516 * x5 - x1523 * x6)
            + x36 * (x1516 - x1517 + x1522)
            - x6 * (x1523 * x5 + x1525 - x6 * (x1521 * x5 - x1527 * x6))
        )
    )
    result[3, 4] = numpy.sum(
        x936
        * (
            x3 * (x1531 + x1533 * x5 - x1539 * x6)
            + x36 * (x1533 - x1534 + x1538)
            - x6 * (x1539 * x5 + x1540 - x6 * (x1537 * x5 - x1541 * x6))
        )
    )
    result[3, 5] = numpy.sum(
        x911
        * (
            x3 * (x1545 + x1547 * x5 - x1553 * x6)
            + x36 * (x1547 - x1548 + x1552)
            - x6 * (x1553 * x5 + x1554 - x6 * (x1551 * x5 - x1555 * x6))
        )
    )
    result[3, 6] = numpy.sum(
        x410 * (x1560 + x3 * (x1559 * x5 - x1561 * x6) - x6 * (x1561 * x5 - x1562 * x6))
    )
    result[3, 7] = numpy.sum(
        x936 * (x1569 + x3 * (x1568 * x5 - x1570 * x6) - x6 * (x1570 * x5 - x1572 * x6))
    )
    result[3, 8] = numpy.sum(
        x936 * (x1576 + x3 * (x1575 * x5 - x1577 * x6) - x6 * (x1577 * x5 - x1578 * x6))
    )
    result[3, 9] = numpy.sum(
        x410 * (x1582 + x3 * (x1581 * x5 - x1583 * x6) - x6 * (x1583 * x5 - x1584 * x6))
    )
    result[3, 10] = numpy.sum(x808 * (x1585 * x3 - x1586 * x6))
    result[3, 11] = numpy.sum(x410 * (x1587 * x3 - x1588 * x6))
    result[3, 12] = numpy.sum(x911 * (x1590 * x3 - x1592 * x6))
    result[3, 13] = numpy.sum(x410 * (x1593 * x3 - x1594 * x6))
    result[3, 14] = numpy.sum(x808 * (x1595 * x3 - x1596 * x6))
    result[4, 0] = numpy.sum(
        x1635
        * (
            x3 * (x1617 * x5 - x1630 * x6 + x48 * (x1607 - x1609 + x1615 + x1631))
            - x6
            * (
                x1630 * x5
                + x48 * (x1616 - x1620 + x1626 + x1634)
                - x6
                * (
                    x1627 * x5
                    + x36 * (x1614 - x1621 + x1624)
                    - x6 * (x1625 * x5 + x1632 - x6 * (x1623 * x5 - x1633 * x6))
                )
            )
            + x98 * (x1617 - x1618 + x1628 - x1629)
        )
    )
    result[4, 1] = numpy.sum(
        x1658
        * (
            x3 * (x1646 * x5 - x1655 * x6 + x36 * (x1640 - x1641 + x1644))
            + x48 * (x1646 - x1647 - x1648 + x1654)
            - x6
            * (
                x1655 * x5
                + x36 * (x1645 - x1649 + x1652)
                - x6 * (x1653 * x5 + x1656 - x6 * (x1651 * x5 - x1657 * x6))
            )
        )
    )
    result[4, 2] = numpy.sum(
        x1658
        * (
            x3 * (x1669 * x5 - x1678 * x6 + x36 * (x1663 - x1664 + x1667))
            + x48 * (x1669 - x1671 + x1677 + x1681)
            - x6
            * (
                x1678 * x5
                + x36 * (x1668 - x1672 + x1675)
                - x6 * (x1676 * x5 + x1679 - x6 * (x1674 * x5 - x1680 * x6))
            )
        )
    )
    result[4, 3] = numpy.sum(
        x1694
        * (
            x3 * (x1684 + x1686 * x5 - x1691 * x6)
            + x36 * (x1686 - x1687 + x1690)
            - x6 * (x1691 * x5 + x1692 - x6 * (x1689 * x5 - x1693 * x6))
        )
    )
    result[4, 4] = numpy.sum(
        x1707
        * (
            x3 * (x1697 + x1699 * x5 - x1704 * x6)
            + x36 * (x1699 - x1700 + x1703)
            - x6 * (x1704 * x5 + x1705 - x6 * (x1702 * x5 - x1706 * x6))
        )
    )
    result[4, 5] = numpy.sum(
        x1694
        * (
            x3 * (x1710 + x1712 * x5 - x1717 * x6)
            + x36 * (x1712 - x1713 + x1716)
            - x6 * (x1717 * x5 + x1718 - x6 * (x1715 * x5 - x1719 * x6))
        )
    )
    result[4, 6] = numpy.sum(
        x1658 * (x1722 + x3 * (x1721 * x5 - x1723 * x6) - x6 * (x1723 * x5 - x1724 * x6))
    )
    result[4, 7] = numpy.sum(
        x1707 * (x1727 + x3 * (x1726 * x5 - x1728 * x6) - x6 * (x1728 * x5 - x1729 * x6))
    )
    result[4, 8] = numpy.sum(
        x1707 * (x1732 + x3 * (x1731 * x5 - x1733 * x6) - x6 * (x1733 * x5 - x1734 * x6))
    )
    result[4, 9] = numpy.sum(
        x1658 * (x1737 + x3 * (x1736 * x5 - x1738 * x6) - x6 * (x1738 * x5 - x1739 * x6))
    )
    result[4, 10] = numpy.sum(x1635 * (x1740 * x3 - x1741 * x6))
    result[4, 11] = numpy.sum(x1658 * (x1742 * x3 - x1743 * x6))
    result[4, 12] = numpy.sum(x1694 * (x1744 * x3 - x1745 * x6))
    result[4, 13] = numpy.sum(x1658 * (x1746 * x3 - x1747 * x6))
    result[4, 14] = numpy.sum(x1635 * (x1748 * x3 - x1749 * x6))
    result[5, 0] = numpy.sum(
        x808
        * (
            x3 * (x1773 * x5 - x1787 * x6 + x48 * (x1762 - x1764 + x1771 + x1788))
            - x6
            * (
                x1787 * x5
                + x48 * (x1772 - x1776 + x1783 + x1791)
                - x6
                * (
                    x1784 * x5
                    + x36 * (x1770 - x1777 + x1781)
                    - x6 * (x1782 * x5 + x1789 - x6 * (x1780 * x5 - x1790 * x6))
                )
            )
            + x98 * (x1773 - x1774 + x1785 - x1786)
        )
    )
    result[5, 1] = numpy.sum(
        x410
        * (
            x3 * (x1807 * x5 - x1818 * x6 + x36 * (x1799 - x1800 + x1805))
            + x48 * (x1807 - x1808 - x1809 + x1817)
            - x6
            * (
                x1818 * x5
                + x36 * (x1806 - x1810 + x1815)
                - x6 * (x1816 * x5 + x1819 - x6 * (x1814 * x5 - x1821 * x6))
            )
        )
    )
    result[5, 2] = numpy.sum(
        x410
        * (
            x3 * (x1834 * x5 - x1844 * x6 + x36 * (x1827 - x1828 + x1832))
            + x48 * (x1834 - x1836 + x1843 + x1847)
            - x6
            * (
                x1844 * x5
                + x36 * (x1833 - x1837 + x1841)
                - x6 * (x1842 * x5 + x1845 - x6 * (x1840 * x5 - x1846 * x6))
            )
        )
    )
    result[5, 3] = numpy.sum(
        x911
        * (
            x3 * (x1853 + x1856 * x5 - x1863 * x6)
            + x36 * (x1856 - x1857 + x1862)
            - x6 * (x1863 * x5 + x1865 - x6 * (x1861 * x5 - x1867 * x6))
        )
    )
    result[5, 4] = numpy.sum(
        x936
        * (
            x3 * (x1872 + x1875 * x5 - x1882 * x6)
            + x36 * (x1875 - x1876 + x1881)
            - x6 * (x1882 * x5 + x1883 - x6 * (x1880 * x5 - x1885 * x6))
        )
    )
    result[5, 5] = numpy.sum(
        x911
        * (
            x3 * (x1892 + x1894 * x5 - x1901 * x6)
            + x36 * (x1894 - x1895 + x1900)
            - x6 * (x1901 * x5 + x1903 - x6 * (x1899 * x5 - x1905 * x6))
        )
    )
    result[5, 6] = numpy.sum(
        x410 * (x1912 + x3 * (x1908 * x5 - x1914 * x6) - x6 * (x1914 * x5 - x1916 * x6))
    )
    result[5, 7] = numpy.sum(
        x936 * (x1922 + x3 * (x1919 * x5 - x1924 * x6) - x6 * (x1924 * x5 - x1926 * x6))
    )
    result[5, 8] = numpy.sum(
        x936 * (x1931 + x3 * (x1930 * x5 - x1933 * x6) - x6 * (x1933 * x5 - x1935 * x6))
    )
    result[5, 9] = numpy.sum(
        x410 * (x1940 + x3 * (x1939 * x5 - x1941 * x6) - x6 * (x1941 * x5 - x1942 * x6))
    )
    result[5, 10] = numpy.sum(x808 * (x1944 * x3 - x1946 * x6))
    result[5, 11] = numpy.sum(x410 * (x1948 * x3 - x1950 * x6))
    result[5, 12] = numpy.sum(x911 * (x1952 * x3 - x1954 * x6))
    result[5, 13] = numpy.sum(x410 * (x1956 * x3 - x1958 * x6))
    result[5, 14] = numpy.sum(x808 * (x1959 * x3 - x1960 * x6))
    result[6, 0] = numpy.sum(
        x107
        * (
            x48 * (x1974 - x1975 - x1976 + x1984)
            + x5 * (x1974 * x5 - x1985 * x6 + x36 * (x1966 - x1967 + x1972))
            - x6
            * (
                x1985 * x5
                + x36 * (x1973 - x1977 + x1982)
                - x6
                * (
                    x1983 * x5
                    + x4 * (x1971 - x1978 + x1979 - x1980)
                    - x6 * (x1981 * x5 - x6 * (x1445 * x748 - x1455 * x8 + 2.0 * x805))
                )
            )
        )
    )
    result[6, 1] = numpy.sum(
        x194
        * (
            x36 * (x1990 - x1991 + x1995)
            + x5
            * (x1990 * x5 - x1996 * x6 + x4 * (x1453 + x1511 + x1986 - x1987 + x1988))
            - x6
            * (
                x1996 * x5
                + x4 * (x1456 + x1524 + x1989 - x1992 + x1993)
                - x6 * (x1994 * x5 - x6 * (x1454 + x1475 * x748 - x1481 * x8 + x1526))
            )
        )
    )
    result[6, 2] = numpy.sum(
        x194
        * (
            x36 * (x2002 - x2003 + x2008)
            + x5 * (x2002 * x5 - x2009 * x6 + x4 * (x1997 - x1998 - x1999 + x2000))
            - x6
            * (
                x2009 * x5
                + x4 * (x2001 - x2004 + x2005 - x2006)
                - x6 * (x2007 * x5 - x6 * (x1500 * x748 - x1506 * x8 + 2.0 * x884))
            )
        )
    )
    result[6, 3] = numpy.sum(
        -x347
        * (
            x4 * (-x2010 + x2011 + x2012 + x2013 - x2014)
            - x5 * (x2014 * x5 - x2015 * x6)
            + x6
            * (x2015 * x5 - x6 * (2.0 * x1480 + x1521 * x748 - x1527 * x8 + 2.0 * x909))
        )
    )
    result[6, 4] = numpy.sum(
        x410
        * (
            x4 * (x1507 + x1566 + x2016 - x2017 + x2018)
            + x5 * (x2018 * x5 - x2019 * x6)
            - x6 * (x2019 * x5 - x6 * (x1505 + x1537 * x748 - x1541 * x8 + x1571))
        )
    )
    result[6, 5] = numpy.sum(
        x347
        * (
            x4 * (x2020 - x2021 - x2022 + x2023)
            + x5 * (x2023 * x5 - x2024 * x6)
            - x6 * (x2024 * x5 - x6 * (x1551 * x748 - x1555 * x8 + 2.0 * x960))
        )
    )
    result[6, 6] = numpy.sum(
        x194
        * (
            x5 * (3.0 * x1514 + x1559 * x748 - x1561 * x8 + 2.0 * x971)
            - x6 * (3.0 * x1525 + x1561 * x748 - x1562 * x8 + 2.0 * x973)
        )
    )
    result[6, 7] = numpy.sum(
        x410
        * (
            x5 * (2.0 * x1531 + x1568 * x748 - x1570 * x8 + 2.0 * x984)
            - x6 * (2.0 * x1540 + x1570 * x748 - x1572 * x8 + 2.0 * x986)
        )
    )
    result[6, 8] = numpy.sum(
        x410
        * (
            x5 * (x1545 + x1575 * x748 - x1577 * x8 + x1589)
            - x6 * (x1554 + x1577 * x748 - x1578 * x8 + x1591)
        )
    )
    result[6, 9] = numpy.sum(
        x194
        * (
            x5 * (2.0 * x1010 + x1581 * x748 - x1583 * x8)
            - x6 * (2.0 * x1012 + x1583 * x748 - x1584 * x8)
        )
    )
    result[6, 10] = numpy.sum(
        x107 * (2.0 * x1018 + 4.0 * x1560 + x1585 * x748 - x1586 * x8)
    )
    result[6, 11] = numpy.sum(
        x194 * (2.0 * x1025 + 3.0 * x1569 + x1587 * x748 - x1588 * x8)
    )
    result[6, 12] = numpy.sum(
        x347 * (2.0 * x1032 + 2.0 * x1576 + x1590 * x748 - x1592 * x8)
    )
    result[6, 13] = numpy.sum(x194 * (2.0 * x1038 + x1582 + x1593 * x748 - x1594 * x8))
    result[6, 14] = numpy.sum(x107 * (2.0 * x1044 + x1595 * x748 - x1596 * x8))
    result[7, 0] = numpy.sum(
        x808
        * (
            x48 * (x2036 - x2037 - x2038 + x2045)
            + x5 * (x2036 * x5 - x2046 * x6 + x36 * (x2029 - x2030 + x2034))
            - x6
            * (
                x2046 * x5
                + x36 * (x2035 - x2039 + x2043)
                - x6
                * (
                    x2044 * x5
                    + x4 * (x1106 + x2033 - x2040 + x2041)
                    - x6 * (x2042 * x5 - x6 * (x1104 + x1623 * x748 - x1633 * x8))
                )
            )
        )
    )
    result[7, 1] = numpy.sum(
        x410
        * (
            x36 * (x2051 - x2052 + x2056)
            + x5
            * (x2051 * x5 - x2057 * x6 + x4 * (x1152 + x1631 + x2047 - x2048 + x2049))
            - x6
            * (
                x2057 * x5
                + x4 * (x1156 + x1634 + x2050 - x2053 + x2054)
                - x6 * (x2055 * x5 - x6 * (x1153 + x1632 + x1651 * x748 - x1657 * x8))
            )
        )
    )
    result[7, 2] = numpy.sum(
        x410
        * (
            x36 * (x2062 - x2063 + x2067)
            + x5 * (x2062 * x5 - x2068 * x6 + x4 * (x1193 + x2058 - x2059 + x2060))
            - x6
            * (
                x2068 * x5
                + x4 * (x1196 + x2061 - x2064 + x2065)
                - x6 * (x2066 * x5 - x6 * (x1194 + x1674 * x748 - x1680 * x8))
            )
        )
    )
    result[7, 3] = numpy.sum(
        x911
        * (
            x4 * (x1226 + x2069 - x2070 - x2071 + x2072)
            + x5 * (x2072 * x5 - x2073 * x6)
            - x6 * (x2073 * x5 - x6 * (x1230 + 2.0 * x1656 + x1689 * x748 - x1693 * x8))
        )
    )
    result[7, 4] = numpy.sum(
        x936
        * (
            x4 * (x1260 + x1681 + x2074 - x2075 + x2076)
            + x5 * (x2076 * x5 - x2077 * x6)
            - x6 * (x2077 * x5 - x6 * (x1263 + x1679 + x1702 * x748 - x1706 * x8))
        )
    )
    result[7, 5] = numpy.sum(
        x911
        * (
            x4 * (x1286 + x2078 - x2079 + x2080)
            + x5 * (x2080 * x5 - x2081 * x6)
            - x6 * (x2081 * x5 - x6 * (x1289 + x1715 * x748 - x1719 * x8))
        )
    )
    result[7, 6] = numpy.sum(
        x410
        * (
            x5 * (x1308 + 3.0 * x1684 + x1721 * x748 - x1723 * x8)
            - x6 * (x1311 + 3.0 * x1692 + x1723 * x748 - x1724 * x8)
        )
    )
    result[7, 7] = numpy.sum(
        x936
        * (
            x5 * (x1329 + 2.0 * x1697 + x1726 * x748 - x1728 * x8)
            - x6 * (x1332 + 2.0 * x1705 + x1728 * x748 - x1729 * x8)
        )
    )
    result[7, 8] = numpy.sum(
        x936
        * (
            x5 * (x1349 + x1710 + x1731 * x748 - x1733 * x8)
            - x6 * (x1351 + x1718 + x1733 * x748 - x1734 * x8)
        )
    )
    result[7, 9] = numpy.sum(
        x410
        * (
            x5 * (x1363 + x1736 * x748 - x1738 * x8)
            - x6 * (x1365 + x1738 * x748 - x1739 * x8)
        )
    )
    result[7, 10] = numpy.sum(x808 * (x1373 + 4.0 * x1722 + x1740 * x748 - x1741 * x8))
    result[7, 11] = numpy.sum(x410 * (x1384 + 3.0 * x1727 + x1742 * x748 - x1743 * x8))
    result[7, 12] = numpy.sum(x911 * (x1394 + 2.0 * x1732 + x1744 * x748 - x1745 * x8))
    result[7, 13] = numpy.sum(x410 * (x1403 + x1737 + x1746 * x748 - x1747 * x8))
    result[7, 14] = numpy.sum(x808 * (x1412 + x1748 * x748 - x1749 * x8))
    result[8, 0] = numpy.sum(
        x808
        * (
            x48 * (x2091 - x2092 - x2093 + x2099)
            + x5 * (x2091 * x5 - x2100 * x6 + x36 * (x2085 - x2086 + x2089))
            - x6
            * (
                x2100 * x5
                + x36 * (x2090 - x2094 + x2097)
                - x6
                * (
                    x2098 * x5
                    + x4 * (x1812 + x2088 - x2095)
                    - x6 * (x2096 * x5 - x6 * (x1780 * x748 + x1820))
                )
            )
        )
    )
    result[8, 1] = numpy.sum(
        x410
        * (
            x36 * (x2104 - x2105 + x2108)
            + x5 * (x2104 * x5 - x2109 * x6 + x4 * (x1852 - x2101 + x2102))
            - x6
            * (
                x2109 * x5
                + x4 * (x1864 + x2103 - x2106)
                - x6 * (x2107 * x5 - x6 * (x1814 * x748 + x1866))
            )
        )
    )
    result[8, 2] = numpy.sum(
        x410
        * (
            x36 * (x2113 - x2114 + x2117)
            + x5 * (x2113 * x5 - x2118 * x6 + x4 * (x1868 - x2110 + x2111))
            - x6
            * (
                x2118 * x5
                + x4 * (x1878 + x2112 - x2115)
                - x6 * (x2116 * x5 - x6 * (x1840 * x748 + x1884))
            )
        )
    )
    result[8, 3] = numpy.sum(
        x911
        * (
            x4 * (x1911 - x2119 + x2120)
            + x5 * (x2120 * x5 - x2121 * x6)
            - x6 * (x2121 * x5 - x6 * (x1861 * x748 + x1915))
        )
    )
    result[8, 4] = numpy.sum(
        x936
        * (
            x4 * (x1921 - x2122 + x2123)
            + x5 * (x2123 * x5 - x2124 * x6)
            - x6 * (x2124 * x5 - x6 * (x1880 * x748 + x1925))
        )
    )
    result[8, 5] = numpy.sum(
        x911
        * (
            x4 * (x1927 - x2125 + x2126)
            + x5 * (x2126 * x5 - x2127 * x6)
            - x6 * (x2127 * x5 - x6 * (x1899 * x748 + x1934))
        )
    )
    result[8, 6] = numpy.sum(
        x410 * (x5 * (x1908 * x748 + x1943) - x6 * (x1914 * x748 + x1945))
    )
    result[8, 7] = numpy.sum(
        x936 * (x5 * (x1919 * x748 + x1947) - x6 * (x1924 * x748 + x1949))
    )
    result[8, 8] = numpy.sum(
        x936 * (x5 * (x1930 * x748 + x1951) - x6 * (x1933 * x748 + x1953))
    )
    result[8, 9] = numpy.sum(
        x410 * (x5 * (x1939 * x748 + x1955) - x6 * (x1941 * x748 + x1957))
    )
    result[8, 10] = numpy.sum(x808 * (4.0 * x1912 + x1944 * x748 - x1946 * x8))
    result[8, 11] = numpy.sum(x410 * (3.0 * x1922 + x1948 * x748 - x1950 * x8))
    result[8, 12] = numpy.sum(x911 * (2.0 * x1931 + x1952 * x748 - x1954 * x8))
    result[8, 13] = numpy.sum(x410 * (x1940 + x1956 * x748 - x1958 * x8))
    result[8, 14] = numpy.sum(x808 * (x1959 * x748 - x1960 * x8))
    result[9, 0] = numpy.sum(
        x107
        * (
            x48 * (x2142 - x2144 + x2152 + x2156)
            + x5 * (x2142 * x5 - x2153 * x6 + x36 * (x2134 - x2135 + x2140))
            - x6
            * (
                x2153 * x5
                + x36 * (x2141 - x2145 + x2150)
                - x6 * (x2151 * x5 + x2154 - x6 * (x2149 * x5 - x2155 * x6))
            )
        )
    )
    result[9, 1] = numpy.sum(
        x194
        * (
            x36 * (x2162 - x2163 + x2167)
            + x5 * (x2160 + x2162 * x5 - x2168 * x6)
            - x6 * (x2168 * x5 + x2169 - x6 * (x2166 * x5 - x2170 * x6))
        )
    )
    result[9, 2] = numpy.sum(
        x194
        * (
            x36 * (x2176 - x2177 + x2181)
            + x5 * (x2174 + x2176 * x5 - x2182 * x6)
            - x6 * (x2182 * x5 + x2183 - x6 * (x2180 * x5 - x2184 * x6))
        )
    )
    result[9, 3] = numpy.sum(
        x347 * (x2188 + x5 * (x2187 * x5 - x2189 * x6) - x6 * (x2189 * x5 - x2190 * x6))
    )
    result[9, 4] = numpy.sum(
        x410 * (x2194 + x5 * (x2193 * x5 - x2195 * x6) - x6 * (x2195 * x5 - x2196 * x6))
    )
    result[9, 5] = numpy.sum(
        x347 * (x2202 + x5 * (x2201 * x5 - x2203 * x6) - x6 * (x2203 * x5 - x2204 * x6))
    )
    result[9, 6] = numpy.sum(x194 * (x2205 * x5 - x2206 * x6))
    result[9, 7] = numpy.sum(x410 * (x2207 * x5 - x2208 * x6))
    result[9, 8] = numpy.sum(x410 * (x2209 * x5 - x2210 * x6))
    result[9, 9] = numpy.sum(x194 * (x2211 * x5 - x2212 * x6))
    result[9, 10] = numpy.sum(x107 * (x110 * x2205 + 3.0 * x2188 - x2206 * x8))
    result[9, 11] = numpy.sum(x194 * (x110 * x2207 + 2.0 * x2194 - x2208 * x8))
    result[9, 12] = numpy.sum(x347 * (x110 * x2209 + x2202 - x2210 * x8))
    result[9, 13] = numpy.sum(x194 * (x110 * x2211 - x2212 * x8))
    result[9, 14] = numpy.sum(
        x107 * (-x10 * x1960 + x1047 * x1959 + 2.0 * x1412 + 4.0 * x1940)
    )
    return result


def coulomb3d_40(ax, da, A, bx, db, B, C):
    """Cartesian (gs) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 1), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = 0.5 / (ax + bx)
    x5 = -x2 - C[0]
    x6 = -x1 * (ax * A[1] + bx * B[1])
    x7 = -x6 - C[1]
    x8 = -x1 * (ax * A[2] + bx * B[2])
    x9 = -x8 - C[2]
    x10 = x0 * (x5**2 + x7**2 + x9**2)
    x11 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x12 = x11 * boys(1, x10)
    x13 = x11 * boys(0, x10)
    x14 = x4 * (-x12 + x13)
    x15 = -x12 * x5 + x13 * x3
    x16 = x11 * boys(2, x10)
    x17 = x16 * x5
    x18 = x12 * x3
    x19 = -x17 + x18
    x20 = x14 + x15 * x3 - x19 * x5
    x21 = x4 * (x12 - x16)
    x22 = x19 * x3
    x23 = x11 * boys(3, x10)
    x24 = x23 * x5
    x25 = x16 * x3
    x26 = -x24 + x25
    x27 = x26 * x5
    x28 = x21 + x22 - x27
    x29 = 2.0 * x4
    x30 = x4 * (x16 - x23)
    x31 = x11 * boys(4, x10)
    x32 = -x21
    x33 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**1.5)
    x34 = 4.41641957979107 * x33
    x35 = -x6 - A[1]
    x36 = x12 * x35
    x37 = x16 * x7
    x38 = -x12 * x7 + x13 * x35
    x39 = x4 * (-x36 + x37 + x38)
    x40 = x36 - x37
    x41 = x3 * x38 - x40 * x5
    x42 = x3 * x40
    x43 = x23 * x7
    x44 = x16 * x35
    x45 = -x43 + x44
    x46 = x45 * x5
    x47 = x42 - x46
    x48 = x4 * (x40 + x43 - x44)
    x49 = x23 * x35 - x31 * x7
    x50 = 11.6847478934435 * x33
    x51 = -x8 - A[2]
    x52 = x12 * x51
    x53 = x16 * x9
    x54 = -x12 * x9 + x13 * x51
    x55 = x4 * (-x52 + x53 + x54)
    x56 = x52 - x53
    x57 = x3 * x54 - x5 * x56
    x58 = x3 * x56
    x59 = x23 * x9
    x60 = x16 * x51
    x61 = -x59 + x60
    x62 = x5 * x61
    x63 = x58 - x62
    x64 = x4 * (x56 + x59 - x60)
    x65 = x23 * x51 - x31 * x9
    x66 = x45 * x7
    x67 = x35 * x40
    x68 = x14 + x35 * x38 - x40 * x7
    x69 = x4 * (x32 + x66 - x67 + x68)
    x70 = x21 - x66 + x67
    x71 = x30 + x35 * x45 - x49 * x7
    x72 = 15.084944665313 * x33
    x73 = x61 * x7
    x74 = x35 * x56
    x75 = x35 * x54 - x56 * x7
    x76 = x4 * (x73 - x74 + x75)
    x77 = -x73 + x74
    x78 = x35 * x61 - x65 * x7
    x79 = 26.1278905896872 * x33
    x80 = x61 * x9
    x81 = x51 * x56
    x82 = x14 + x51 * x54 - x56 * x9
    x83 = x4 * (x32 + x80 - x81 + x82)
    x84 = x21 - x80 + x81
    x85 = x30 + x51 * x61 - x65 * x9
    x86 = x35 * x68 + 2.0 * x39 - x7 * x70
    x87 = x35 * x70 + 2.0 * x48 - x7 * x71
    x88 = x35 * x75 + x55 - x7 * x77
    x89 = x35 * x77 + x64 - x7 * x78
    x90 = x35 * x82 - x7 * x84
    x91 = x35 * x84 - x7 * x85
    x92 = x51 * x82 + 2.0 * x55 - x84 * x9
    x93 = x51 * x84 + 2.0 * x64 - x85 * x9

    # 15 item(s)
    result[0, 0] = numpy.sum(
        x34
        * (
            x3 * (x20 * x3 - x28 * x5 + x29 * (x15 + x17 - x18))
            + 3.0 * x4 * (x20 - x22 + x27 + x32)
            - x5
            * (
                x28 * x3
                + x29 * (x19 + x24 - x25)
                - x5 * (x26 * x3 + x30 - x5 * (x23 * x3 - x31 * x5))
            )
        )
    )
    result[1, 0] = numpy.sum(
        x50
        * (
            x29 * (x41 - x42 + x46)
            + x3 * (x3 * x41 + x39 - x47 * x5)
            - x5 * (x3 * x47 + x48 - x5 * (x3 * x45 - x49 * x5))
        )
    )
    result[2, 0] = numpy.sum(
        x50
        * (
            x29 * (x57 - x58 + x62)
            + x3 * (x3 * x57 - x5 * x63 + x55)
            - x5 * (x3 * x63 - x5 * (x3 * x61 - x5 * x65) + x64)
        )
    )
    result[3, 0] = numpy.sum(
        x72 * (x3 * (x3 * x68 - x5 * x70) - x5 * (x3 * x70 - x5 * x71) + x69)
    )
    result[4, 0] = numpy.sum(
        x79 * (x3 * (x3 * x75 - x5 * x77) - x5 * (x3 * x77 - x5 * x78) + x76)
    )
    result[5, 0] = numpy.sum(
        x72 * (x3 * (x3 * x82 - x5 * x84) - x5 * (x3 * x84 - x5 * x85) + x83)
    )
    result[6, 0] = numpy.sum(x50 * (x3 * x86 - x5 * x87))
    result[7, 0] = numpy.sum(x79 * (x3 * x88 - x5 * x89))
    result[8, 0] = numpy.sum(x79 * (x3 * x90 - x5 * x91))
    result[9, 0] = numpy.sum(x50 * (x3 * x92 - x5 * x93))
    result[10, 0] = numpy.sum(x34 * (x35 * x86 + 3.0 * x69 - x7 * x87))
    result[11, 0] = numpy.sum(x50 * (x35 * x88 - x7 * x89 + 2.0 * x76))
    result[12, 0] = numpy.sum(x72 * (x35 * x90 - x7 * x91 + x83))
    result[13, 0] = numpy.sum(x50 * (x35 * x92 - x7 * x93))
    result[14, 0] = numpy.sum(x34 * (x51 * x92 + 3.0 * x83 - x9 * x93))
    return result


def coulomb3d_41(ax, da, A, bx, db, B, C):
    """Cartesian (gp) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 3), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = ax + bx
    x2 = x1 ** (-1.0)
    x3 = -x2 * (ax * A[0] + bx * B[0])
    x4 = -x3 - A[0]
    x5 = -x3 - C[0]
    x6 = -x2 * (ax * A[1] + bx * B[1])
    x7 = -x6 - C[1]
    x8 = -x2 * (ax * A[2] + bx * B[2])
    x9 = -x8 - C[2]
    x10 = x1 * (x5**2 + x7**2 + x9**2)
    x11 = (
        6.28318530717959
        * x2
        * numpy.exp(
            -ax * bx * x2 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x12 = x11 * boys(1, x10)
    x13 = x11 * boys(0, x10)
    x14 = x0 * (-x12 + x13)
    x15 = -x12 * x5
    x16 = x13 * x4 + x15
    x17 = x11 * boys(2, x10)
    x18 = x17 * x5
    x19 = -x18
    x20 = x12 * x4
    x21 = x19 + x20
    x22 = x14 + x16 * x4 - x21 * x5
    x23 = x11 * boys(3, x10)
    x24 = x0 * (x17 - x23)
    x25 = x23 * x5
    x26 = -x25
    x27 = x17 * x4
    x28 = x26 + x27
    x29 = x28 * x4
    x30 = x11 * boys(4, x10)
    x31 = x30 * x5
    x32 = -x31
    x33 = x23 * x4
    x34 = x5 * (x32 + x33)
    x35 = x0 * (x12 - x17)
    x36 = x21 * x4
    x37 = x28 * x5
    x38 = x35 + x36 - x37
    x39 = x0 * (x21 + x25 - x27)
    x40 = x0 * (x16 + x18 - x20)
    x41 = -x35
    x42 = -x3 - B[0]
    x43 = x12 * x42
    x44 = x13 * x42 + x15
    x45 = x19 + x43
    x46 = x14 + x4 * x44 - x45 * x5
    x47 = x4 * x45
    x48 = x17 * x42
    x49 = x26 + x48
    x50 = x49 * x5
    x51 = x35 + x47 - x50
    x52 = x0 * (x18 - x43 + x44) + x4 * x46 + x40 - x5 * x51
    x53 = x0 * (x25 + x45 - x48)
    x54 = x4 * x51
    x55 = x4 * x49
    x56 = x23 * x42
    x57 = x32 + x56
    x58 = x5 * x57
    x59 = x24 + x55 - x58
    x60 = x5 * x59
    x61 = x39 + x53 + x54 - x60
    x62 = 2.0 * x0
    x63 = -x24
    x64 = x0 * (x23 - x30)
    x65 = x11 * boys(5, x10)
    x66 = 3.0 * x0
    x67 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**2.5)
    x68 = 8.83283915958214 * x67
    x69 = -x6 - B[1]
    x70 = x12 * x69
    x71 = x17 * x7
    x72 = -x12 * x7
    x73 = x13 * x69 + x72
    x74 = x0 * (-x70 + x71 + x73)
    x75 = -x71
    x76 = x70 + x75
    x77 = x4 * x73 - x5 * x76
    x78 = x4 * x76
    x79 = x23 * x7
    x80 = -x79
    x81 = x17 * x69
    x82 = x80 + x81
    x83 = x5 * x82
    x84 = x78 - x83
    x85 = x4 * x77 - x5 * x84 + x74
    x86 = x0 * (x76 + x79 - x81)
    x87 = x4 * x84
    x88 = x4 * x82
    x89 = x30 * x7
    x90 = -x89
    x91 = x23 * x69
    x92 = x90 + x91
    x93 = x5 * x92
    x94 = x88 - x93
    x95 = x5 * x94
    x96 = x86 + x87 - x95
    x97 = x0 * (x82 + x89 - x91)
    x98 = -x65 * x7
    x99 = x30 * x69 + x98
    x100 = -x86
    x101 = -x8 - B[2]
    x102 = x101 * x12
    x103 = x17 * x9
    x104 = -x12 * x9
    x105 = x101 * x13 + x104
    x106 = x0 * (-x102 + x103 + x105)
    x107 = -x103
    x108 = x102 + x107
    x109 = x105 * x4 - x108 * x5
    x110 = x108 * x4
    x111 = x23 * x9
    x112 = -x111
    x113 = x101 * x17
    x114 = x112 + x113
    x115 = x114 * x5
    x116 = x110 - x115
    x117 = x106 + x109 * x4 - x116 * x5
    x118 = x0 * (x108 + x111 - x113)
    x119 = x116 * x4
    x120 = x114 * x4
    x121 = x30 * x9
    x122 = -x121
    x123 = x101 * x23
    x124 = x122 + x123
    x125 = x124 * x5
    x126 = x120 - x125
    x127 = x126 * x5
    x128 = x118 + x119 - x127
    x129 = x0 * (x114 + x121 - x123)
    x130 = -x65 * x9
    x131 = x101 * x30 + x130
    x132 = -x118
    x133 = -x6 - A[1]
    x134 = x12 * x133
    x135 = x13 * x133 + x72
    x136 = x0 * (-x134 + x135 + x71)
    x137 = x134 + x75
    x138 = -x137 * x5
    x139 = x135 * x4 + x138
    x140 = x133 * x17
    x141 = x140 + x80
    x142 = x141 * x4
    x143 = x133 * x23
    x144 = x143 + x90
    x145 = x144 * x5
    x146 = -x145
    x147 = x0 * (x137 - x140 + x79)
    x148 = -x147
    x149 = x137 * x4
    x150 = x141 * x5
    x151 = -x150
    x152 = x149 + x151
    x153 = x137 * x42
    x154 = x135 * x42 + x138
    x155 = x151 + x153
    x156 = x136 + x154 * x4 - x155 * x5
    x157 = x155 * x4
    x158 = x141 * x42
    x159 = x146 + x158
    x160 = x159 * x5
    x161 = x147 + x157 - x160
    x162 = x0 * (x141 - x143 + x89)
    x163 = x133 * x30 + x98
    x164 = 23.3694957868871 * x67
    x165 = x7 * x82
    x166 = x133 * x76
    x167 = x133 * x73 + x14 - x7 * x76
    x168 = x0 * (x165 - x166 + x167 + x41)
    x169 = -x165 + x166 + x35
    x170 = x167 * x4 - x169 * x5
    x171 = x169 * x4
    x172 = x133 * x82
    x173 = x7 * x92
    x174 = x172 - x173 + x24
    x175 = x174 * x5
    x176 = x171 - x175
    x177 = x0 * (x169 - x172 + x173 + x63)
    x178 = x133 * x92 + x64 - x7 * x99
    x179 = x114 * x7
    x180 = x108 * x133
    x181 = x105 * x133 - x108 * x7
    x182 = x0 * (x179 - x180 + x181)
    x183 = -x179 + x180
    x184 = x181 * x4 - x183 * x5
    x185 = x183 * x4
    x186 = x114 * x133
    x187 = x124 * x7
    x188 = x186 - x187
    x189 = x188 * x5
    x190 = x185 - x189
    x191 = x0 * (x183 - x186 + x187)
    x192 = x124 * x133 - x131 * x7
    x193 = -x8 - A[2]
    x194 = x12 * x193
    x195 = x104 + x13 * x193
    x196 = x0 * (x103 - x194 + x195)
    x197 = x107 + x194
    x198 = -x197 * x5
    x199 = x195 * x4 + x198
    x200 = x17 * x193
    x201 = x112 + x200
    x202 = x201 * x4
    x203 = x193 * x23
    x204 = x122 + x203
    x205 = x204 * x5
    x206 = -x205
    x207 = x0 * (x111 + x197 - x200)
    x208 = -x207
    x209 = x197 * x4
    x210 = x201 * x5
    x211 = -x210
    x212 = x209 + x211
    x213 = x197 * x42
    x214 = x195 * x42 + x198
    x215 = x211 + x213
    x216 = x196 + x214 * x4 - x215 * x5
    x217 = x215 * x4
    x218 = x201 * x42
    x219 = x206 + x218
    x220 = x219 * x5
    x221 = x207 + x217 - x220
    x222 = x0 * (x121 + x201 - x203)
    x223 = x130 + x193 * x30
    x224 = x201 * x7
    x225 = x197 * x69
    x226 = -x197 * x7
    x227 = x195 * x69 + x226
    x228 = x0 * (x224 - x225 + x227)
    x229 = -x224
    x230 = x225 + x229
    x231 = x227 * x4 - x230 * x5
    x232 = x230 * x4
    x233 = x201 * x69
    x234 = x204 * x7
    x235 = -x234
    x236 = x233 + x235
    x237 = x236 * x5
    x238 = x232 - x237
    x239 = x0 * (x230 - x233 + x234)
    x240 = -x223 * x7
    x241 = x204 * x69 + x240
    x242 = x114 * x9
    x243 = x108 * x193
    x244 = x105 * x193 - x108 * x9 + x14
    x245 = x0 * (x242 - x243 + x244 + x41)
    x246 = -x242 + x243 + x35
    x247 = x244 * x4 - x246 * x5
    x248 = x246 * x4
    x249 = x114 * x193
    x250 = x124 * x9
    x251 = x24 + x249 - x250
    x252 = x251 * x5
    x253 = x248 - x252
    x254 = x0 * (x246 - x249 + x250 + x63)
    x255 = x124 * x193 - x131 * x9 + x64
    x256 = x133 * x135 - x137 * x7 + x14
    x257 = x133 * x141
    x258 = x144 * x7
    x259 = x24 + x257 - x258
    x260 = x259 * x5
    x261 = x133 * x137
    x262 = x141 * x7
    x263 = x261 - x262 + x35
    x264 = -x263 * x5
    x265 = x263 * x42
    x266 = x256 * x42 + x264
    x267 = x0 * (x256 - x261 + x262 + x41)
    x268 = -x260 + x265
    x269 = x0 * (-x257 + x258 + x263 + x63)
    x270 = x133 * x144 - x163 * x7 + x64
    x271 = 30.169889330626 * x67
    x272 = x174 * x7
    x273 = x133 * x169
    x274 = x133 * x167 + x136 - x169 * x7 + x74
    x275 = x0 * (x100 + x148 + x272 - x273 + x274)
    x276 = x147 - x272 + x273 + x86
    x277 = x133 * x174 + x162 - x178 * x7 + x97
    x278 = x188 * x7
    x279 = x133 * x183
    x280 = x106 + x133 * x181 - x183 * x7
    x281 = x0 * (x132 + x278 - x279 + x280)
    x282 = x118 - x278 + x279
    x283 = x129 + x133 * x188 - x192 * x7
    x284 = x133 * x195 + x226
    x285 = x133 * x201
    x286 = x235 + x285
    x287 = x286 * x5
    x288 = x133 * x197
    x289 = x229 + x288
    x290 = -x289 * x5
    x291 = x289 * x42
    x292 = x284 * x42 + x290
    x293 = x0 * (x224 + x284 - x288)
    x294 = -x287 + x291
    x295 = x0 * (x234 - x285 + x289)
    x296 = x133 * x204 + x240
    x297 = 52.2557811793745 * x67
    x298 = x236 * x7
    x299 = x133 * x230
    x300 = x133 * x227 + x196 - x230 * x7
    x301 = x0 * (x208 + x298 - x299 + x300)
    x302 = x207 - x298 + x299
    x303 = x133 * x236 + x222 - x241 * x7
    x304 = x251 * x7
    x305 = x133 * x246
    x306 = x133 * x244 - x246 * x7
    x307 = x0 * (x304 - x305 + x306)
    x308 = -x304 + x305
    x309 = x133 * x251 - x255 * x7
    x310 = x14 + x193 * x195 - x197 * x9
    x311 = x193 * x201
    x312 = x204 * x9
    x313 = x24 + x311 - x312
    x314 = x313 * x5
    x315 = x193 * x197
    x316 = x201 * x9
    x317 = x315 - x316 + x35
    x318 = -x317 * x5
    x319 = x317 * x42
    x320 = x310 * x42 + x318
    x321 = x0 * (x310 - x315 + x316 + x41)
    x322 = -x314 + x319
    x323 = x0 * (-x311 + x312 + x317 + x63)
    x324 = x193 * x204 - x223 * x9 + x64
    x325 = x313 * x7
    x326 = x317 * x69
    x327 = -x317 * x7
    x328 = x310 * x69 + x327
    x329 = x0 * (x325 - x326 + x328)
    x330 = -x325
    x331 = x326 + x330
    x332 = -x324 * x7
    x333 = x313 * x69 + x332
    x334 = x251 * x9
    x335 = x193 * x246
    x336 = x106 + x193 * x244 + x196 - x246 * x9
    x337 = x0 * (x132 + x208 + x334 - x335 + x336)
    x338 = x118 + x207 - x334 + x335
    x339 = x129 + x193 * x251 + x222 - x255 * x9
    x340 = x259 * x7
    x341 = x133 * x263
    x342 = 2.0 * x147
    x343 = x133 * x256 + 2.0 * x136 - x263 * x7
    x344 = x0 * (x340 - x341 - x342 + x343)
    x345 = -x340 + x341 + x342
    x346 = x133 * x259 + 2.0 * x162 - x270 * x7
    x347 = x133 * x274 + 2.0 * x168 + x267 - x276 * x7
    x348 = x133 * x276 + 2.0 * x177 + x269 - x277 * x7
    x349 = x133 * x280 + 2.0 * x182 - x282 * x7
    x350 = x133 * x282 + 2.0 * x191 - x283 * x7
    x351 = x286 * x7
    x352 = x133 * x289
    x353 = x133 * x284 + x196 - x289 * x7
    x354 = x0 * (x208 + x351 - x352 + x353)
    x355 = x207 - x351 + x352
    x356 = x133 * x286 + x222 - x296 * x7
    x357 = x133 * x300 + x228 + x293 - x302 * x7
    x358 = x133 * x302 + x239 + x295 - x303 * x7
    x359 = x133 * x306 + x245 - x308 * x7
    x360 = x133 * x308 + x254 - x309 * x7
    x361 = x133 * x317
    x362 = x133 * x310 + x327
    x363 = x0 * (x325 - x361 + x362)
    x364 = x330 + x361
    x365 = x133 * x313 + x332
    x366 = x133 * x328 + x321 - x331 * x7
    x367 = x133 * x331 + x323 - x333 * x7
    x368 = x133 * x336 - x338 * x7
    x369 = x133 * x338 - x339 * x7
    x370 = x313 * x9
    x371 = x193 * x317
    x372 = 2.0 * x207
    x373 = x193 * x310 + 2.0 * x196 - x317 * x9
    x374 = x0 * (x370 - x371 - x372 + x373)
    x375 = -x370 + x371 + x372
    x376 = x193 * x313 + 2.0 * x222 - x324 * x9
    x377 = -x375 * x7
    x378 = x373 * x69 + x377
    x379 = -x376 * x7
    x380 = x375 * x69 + x379
    x381 = x193 * x336 + 2.0 * x245 + x321 - x338 * x9
    x382 = x193 * x338 + 2.0 * x254 + x323 - x339 * x9
    x383 = x193 * x373 + 3.0 * x321 - x375 * x9
    x384 = x193 * x375 + 3.0 * x323 - x376 * x9

    # 45 item(s)
    result[0, 0] = numpy.sum(
        x68
        * (
            x0
            * (
                x22 * x4
                - x38 * x4
                - x38 * x5
                - 2.0 * x39
                + 2.0 * x40
                + x5 * (x24 + x29 - x34)
            )
            + x4
            * (
                x0 * (x22 - x36 + x37 + x41)
                + x4 * x52
                - x5 * x61
                + x62 * (x41 + x46 - x47 + x50)
            )
            - x5
            * (
                x0 * (-x29 + x34 + x38 + x63)
                + x4 * x61
                - x5
                * (
                    x0 * (x28 + x31 - x33)
                    + x0 * (x31 + x49 - x56)
                    + x4 * x59
                    - x5 * (x4 * x57 - x5 * (x30 * x42 - x5 * x65) + x64)
                )
                + x62 * (x51 - x55 + x58 + x63)
            )
            - x66 * (x39 - x52 + x53 + x54 - x60)
        )
    )
    result[0, 1] = numpy.sum(
        x68
        * (
            x4 * (x4 * x85 - x5 * x96 + x62 * (x77 - x78 + x83))
            - x5
            * (
                x4 * x96
                - x5 * (x4 * x94 - x5 * (x4 * x92 - x5 * x99) + x97)
                + x62 * (x84 - x88 + x93)
            )
            + x66 * (x100 + x85 - x87 + x95)
        )
    )
    result[0, 2] = numpy.sum(
        x68
        * (
            x4 * (x117 * x4 - x128 * x5 + x62 * (x109 - x110 + x115))
            - x5
            * (
                x128 * x4
                - x5 * (x126 * x4 + x129 - x5 * (x124 * x4 - x131 * x5))
                + x62 * (x116 - x120 + x125)
            )
            + x66 * (x117 - x119 + x127 + x132)
        )
    )
    result[1, 0] = numpy.sum(
        x164
        * (
            x0 * (x136 + x139 * x4 + x148 - x152 * x4 - x152 * x5 + x5 * (x142 + x146))
            + x4
            * (
                x0 * (x139 - x149 + x150)
                + x0 * (x150 - x153 + x154)
                + x156 * x4
                - x161 * x5
            )
            - x5
            * (
                x0 * (-x142 + x145 + x152)
                + x0 * (x145 + x155 - x158)
                + x161 * x4
                - x5 * (x159 * x4 + x162 - x5 * (x144 * x42 - x163 * x5))
            )
            + x62 * (x148 + x156 - x157 + x160)
        )
    )
    result[1, 1] = numpy.sum(
        x164
        * (
            x4 * (x168 + x170 * x4 - x176 * x5)
            - x5 * (x176 * x4 + x177 - x5 * (x174 * x4 - x178 * x5))
            + x62 * (x170 - x171 + x175)
        )
    )
    result[1, 2] = numpy.sum(
        x164
        * (
            x4 * (x182 + x184 * x4 - x190 * x5)
            - x5 * (x190 * x4 + x191 - x5 * (x188 * x4 - x192 * x5))
            + x62 * (x184 - x185 + x189)
        )
    )
    result[2, 0] = numpy.sum(
        x164
        * (
            x0 * (x196 + x199 * x4 + x208 - x212 * x4 - x212 * x5 + x5 * (x202 + x206))
            + x4
            * (
                x0 * (x199 - x209 + x210)
                + x0 * (x210 - x213 + x214)
                + x216 * x4
                - x221 * x5
            )
            - x5
            * (
                x0 * (-x202 + x205 + x212)
                + x0 * (x205 + x215 - x218)
                + x221 * x4
                - x5 * (x219 * x4 + x222 - x5 * (x204 * x42 - x223 * x5))
            )
            + x62 * (x208 + x216 - x217 + x220)
        )
    )
    result[2, 1] = numpy.sum(
        x164
        * (
            x4 * (x228 + x231 * x4 - x238 * x5)
            - x5 * (x238 * x4 + x239 - x5 * (x236 * x4 - x241 * x5))
            + x62 * (x231 - x232 + x237)
        )
    )
    result[2, 2] = numpy.sum(
        x164
        * (
            x4 * (x245 + x247 * x4 - x253 * x5)
            - x5 * (x253 * x4 + x254 - x5 * (x251 * x4 - x255 * x5))
            + x62 * (x247 - x248 + x252)
        )
    )
    result[3, 0] = numpy.sum(
        x271
        * (
            x0 * (x260 - x265 + x266)
            + x0 * (x256 * x4 + x260 - x263 * x4 + x264)
            + x4 * (x266 * x4 + x267 - x268 * x5)
            - x5 * (x268 * x4 + x269 - x5 * (x259 * x42 - x270 * x5))
        )
    )
    result[3, 1] = numpy.sum(
        x271 * (x275 + x4 * (x274 * x4 - x276 * x5) - x5 * (x276 * x4 - x277 * x5))
    )
    result[3, 2] = numpy.sum(
        x271 * (x281 + x4 * (x280 * x4 - x282 * x5) - x5 * (x282 * x4 - x283 * x5))
    )
    result[4, 0] = numpy.sum(
        x297
        * (
            x0 * (x287 - x291 + x292)
            + x0 * (x284 * x4 + x287 - x289 * x4 + x290)
            + x4 * (x292 * x4 + x293 - x294 * x5)
            - x5 * (x294 * x4 + x295 - x5 * (x286 * x42 - x296 * x5))
        )
    )
    result[4, 1] = numpy.sum(
        x297 * (x301 + x4 * (x300 * x4 - x302 * x5) - x5 * (x302 * x4 - x303 * x5))
    )
    result[4, 2] = numpy.sum(
        x297 * (x307 + x4 * (x306 * x4 - x308 * x5) - x5 * (x308 * x4 - x309 * x5))
    )
    result[5, 0] = numpy.sum(
        x271
        * (
            x0 * (x314 - x319 + x320)
            + x0 * (x310 * x4 + x314 - x317 * x4 + x318)
            + x4 * (x320 * x4 + x321 - x322 * x5)
            - x5 * (x322 * x4 + x323 - x5 * (x313 * x42 - x324 * x5))
        )
    )
    result[5, 1] = numpy.sum(
        x271 * (x329 + x4 * (x328 * x4 - x331 * x5) - x5 * (x331 * x4 - x333 * x5))
    )
    result[5, 2] = numpy.sum(
        x271 * (x337 + x4 * (x336 * x4 - x338 * x5) - x5 * (x338 * x4 - x339 * x5))
    )
    result[6, 0] = numpy.sum(
        x164 * (x344 + x4 * (x343 * x42 - x345 * x5) - x5 * (x345 * x42 - x346 * x5))
    )
    result[6, 1] = numpy.sum(x164 * (x347 * x4 - x348 * x5))
    result[6, 2] = numpy.sum(x164 * (x349 * x4 - x350 * x5))
    result[7, 0] = numpy.sum(
        x297 * (x354 + x4 * (x353 * x42 - x355 * x5) - x5 * (x355 * x42 - x356 * x5))
    )
    result[7, 1] = numpy.sum(x297 * (x357 * x4 - x358 * x5))
    result[7, 2] = numpy.sum(x297 * (x359 * x4 - x360 * x5))
    result[8, 0] = numpy.sum(
        x297 * (x363 + x4 * (x362 * x42 - x364 * x5) - x5 * (x364 * x42 - x365 * x5))
    )
    result[8, 1] = numpy.sum(x297 * (x366 * x4 - x367 * x5))
    result[8, 2] = numpy.sum(x297 * (x368 * x4 - x369 * x5))
    result[9, 0] = numpy.sum(
        x164 * (x374 + x4 * (x373 * x42 - x375 * x5) - x5 * (x375 * x42 - x376 * x5))
    )
    result[9, 1] = numpy.sum(x164 * (x378 * x4 - x380 * x5))
    result[9, 2] = numpy.sum(x164 * (x381 * x4 - x382 * x5))
    result[10, 0] = numpy.sum(
        x68
        * (
            x42 * (x133 * x343 + 3.0 * x267 - x345 * x7)
            - x5 * (x133 * x345 + 3.0 * x269 - x346 * x7)
        )
    )
    result[10, 1] = numpy.sum(x68 * (x133 * x347 + 3.0 * x275 + x344 - x348 * x7))
    result[10, 2] = numpy.sum(x68 * (x133 * x349 + 3.0 * x281 - x350 * x7))
    result[11, 0] = numpy.sum(
        x164
        * (
            x42 * (x133 * x353 + 2.0 * x293 - x355 * x7)
            - x5 * (x133 * x355 + 2.0 * x295 - x356 * x7)
        )
    )
    result[11, 1] = numpy.sum(x164 * (x133 * x357 + 2.0 * x301 + x354 - x358 * x7))
    result[11, 2] = numpy.sum(x164 * (x133 * x359 + 2.0 * x307 - x360 * x7))
    result[12, 0] = numpy.sum(
        x271
        * (x42 * (x133 * x362 + x321 - x364 * x7) - x5 * (x133 * x364 + x323 - x365 * x7))
    )
    result[12, 1] = numpy.sum(x271 * (x133 * x366 + x329 + x363 - x367 * x7))
    result[12, 2] = numpy.sum(x271 * (x133 * x368 + x337 - x369 * x7))
    result[13, 0] = numpy.sum(
        x164 * (x42 * (x133 * x373 + x377) - x5 * (x133 * x375 + x379))
    )
    result[13, 1] = numpy.sum(x164 * (x133 * x378 + x374 - x380 * x7))
    result[13, 2] = numpy.sum(x164 * (x133 * x381 - x382 * x7))
    result[14, 0] = numpy.sum(x68 * (x383 * x42 - x384 * x5))
    result[14, 1] = numpy.sum(x68 * (x383 * x69 - x384 * x7))
    result[14, 2] = numpy.sum(x68 * (x193 * x381 + 3.0 * x337 + x374 - x382 * x9))
    return result


def coulomb3d_42(ax, da, A, bx, db, B, C):
    """Cartesian (gd) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 6), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = 0.5 / (ax + bx)
    x5 = -x2 - B[0]
    x6 = -x2 - C[0]
    x7 = -x1 * (ax * A[1] + bx * B[1])
    x8 = -x7 - C[1]
    x9 = -x1 * (ax * A[2] + bx * B[2])
    x10 = -x9 - C[2]
    x11 = x0 * (x10**2 + x6**2 + x8**2)
    x12 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x13 = x12 * boys(2, x11)
    x14 = x13 * x6
    x15 = -x14
    x16 = x12 * boys(1, x11)
    x17 = x16 * x5
    x18 = x15 + x17
    x19 = x18 * x5
    x20 = -x16 * x6
    x21 = x12 * boys(0, x11)
    x22 = x20 + x21 * x5
    x23 = x4 * (-x16 + x21)
    x24 = -x18 * x6 + x23
    x25 = x22 * x5 + x24
    x26 = x12 * boys(3, x11)
    x27 = x26 * x6
    x28 = -x27
    x29 = x13 * x5
    x30 = x28 + x29
    x31 = x30 * x6
    x32 = x4 * (-x13 + x16)
    x33 = -x32
    x34 = x31 + x33
    x35 = -x31 + x32
    x36 = x19 + x35
    x37 = x4 * (x14 - x17 + x22)
    x38 = x25 * x3 - x36 * x6 + 2.0 * x37
    x39 = x3 * x36
    x40 = x30 * x5
    x41 = x4 * (x13 - x26)
    x42 = x12 * boys(4, x11)
    x43 = x42 * x6
    x44 = -x43
    x45 = x26 * x5
    x46 = x44 + x45
    x47 = x46 * x6
    x48 = x41 - x47
    x49 = x40 + x48
    x50 = x49 * x6
    x51 = x4 * (x18 + x27 - x29)
    x52 = 2.0 * x51
    x53 = x39 - x50 + x52
    x54 = x18 * x3
    x55 = x22 * x3 + x24
    x56 = 2.0 * x4
    x57 = x56 * (x34 - x54 + x55)
    x58 = x3 * x38 + x4 * (-x19 + x25 + x34) - x53 * x6 + x57
    x59 = -x41
    x60 = x47 + x59
    x61 = x4 * (x36 - x40 + x60)
    x62 = x3 * x53
    x63 = x3 * x49
    x64 = x46 * x5
    x65 = x4 * (x26 - x42)
    x66 = x12 * boys(5, x11)
    x67 = x6 * x66
    x68 = x42 * x5
    x69 = -x67 + x68
    x70 = x6 * x69
    x71 = x65 - x70
    x72 = x64 + x71
    x73 = x6 * x72
    x74 = x4 * (x30 + x43 - x45)
    x75 = 2.0 * x74
    x76 = x63 - x73 + x75
    x77 = x6 * x76
    x78 = x3 * x30
    x79 = x35 + x54
    x80 = x56 * (x60 - x78 + x79)
    x81 = x61 + x62 - x77 + x80
    x82 = x48 + x78
    x83 = x6 * x82
    x84 = x13 * x3
    x85 = x16 * x3
    x86 = x15 + x85
    x87 = x4 * (x27 - x84 + x86)
    x88 = x3 * x79
    x89 = x20 + x21 * x3
    x90 = x3 * x55 + x37 + x4 * (x14 - x85 + x89) - x6 * x79
    x91 = -x65
    x92 = x70 + x91
    x93 = x4 * (x42 - x66)
    x94 = x12 * boys(6, x11)
    x95 = x3 * x46
    x96 = x6 * (x71 + x95)
    x97 = x26 * x3
    x98 = x28 + x84
    x99 = x4 * (x43 - x97 + x98)
    x100 = x3 * x82
    x101 = x51 - x83 + x87 + x88
    x102 = x6 * x98
    x103 = x3 * x86
    x104 = -x80
    x105 = 3.0 * x4
    x106 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**3.5)
    x107 = 10.1992841329868 * x106
    x108 = -x7 - B[1]
    x109 = x108 * x16
    x110 = x13 * x8
    x111 = -x16 * x8
    x112 = x108 * x21 + x111
    x113 = x4 * (-x109 + x110 + x112)
    x114 = -x110
    x115 = x109 + x114
    x116 = -x115 * x6
    x117 = x112 * x3 + x116
    x118 = x115 * x3
    x119 = x26 * x8
    x120 = -x119
    x121 = x108 * x13
    x122 = x120 + x121
    x123 = x122 * x6
    x124 = -x123
    x125 = x118 + x124
    x126 = x113 + x117 * x3 - x125 * x6
    x127 = x108 * x26
    x128 = x42 * x8
    x129 = x4 * (x122 - x127 + x128)
    x130 = x122 * x3
    x131 = -x128
    x132 = x127 + x131
    x133 = x132 * x6
    x134 = -x133
    x135 = x130 + x134
    x136 = x135 * x3
    x137 = x132 * x3
    x138 = x66 * x8
    x139 = -x138
    x140 = x108 * x42
    x141 = x139 + x140
    x142 = x141 * x6
    x143 = -x142
    x144 = x6 * (x137 + x143)
    x145 = x4 * (x115 + x119 - x121)
    x146 = x125 * x3
    x147 = x135 * x6
    x148 = x145 + x146 - x147
    x149 = x4 * (x125 - x130 + x133)
    x150 = x4 * (x117 - x118 + x123)
    x151 = -x145
    x152 = x115 * x5
    x153 = x112 * x5 + x116
    x154 = x124 + x152
    x155 = x113 + x153 * x3 - x154 * x6
    x156 = x154 * x3
    x157 = x122 * x5
    x158 = x134 + x157
    x159 = x158 * x6
    x160 = x145 + x156 - x159
    x161 = x150 + x155 * x3 - x160 * x6 + x4 * (x123 - x152 + x153)
    x162 = x4 * (x133 + x154 - x157)
    x163 = x160 * x3
    x164 = x158 * x3
    x165 = x132 * x5
    x166 = x143 + x165
    x167 = x166 * x6
    x168 = x129 + x164 - x167
    x169 = x168 * x6
    x170 = x149 + x162 + x163 - x169
    x171 = -x129
    x172 = x4 * (x132 + x138 - x140)
    x173 = -x8 * x94
    x174 = x108 * x66 + x173
    x175 = 17.6656783191643 * x106
    x176 = -x9 - B[2]
    x177 = x16 * x176
    x178 = x10 * x13
    x179 = -x10 * x16
    x180 = x176 * x21 + x179
    x181 = x4 * (-x177 + x178 + x180)
    x182 = -x178
    x183 = x177 + x182
    x184 = -x183 * x6
    x185 = x180 * x3 + x184
    x186 = x183 * x3
    x187 = x10 * x26
    x188 = -x187
    x189 = x13 * x176
    x190 = x188 + x189
    x191 = x190 * x6
    x192 = -x191
    x193 = x186 + x192
    x194 = x181 + x185 * x3 - x193 * x6
    x195 = x176 * x26
    x196 = x10 * x42
    x197 = x4 * (x190 - x195 + x196)
    x198 = x190 * x3
    x199 = -x196
    x200 = x195 + x199
    x201 = x200 * x6
    x202 = -x201
    x203 = x198 + x202
    x204 = x203 * x3
    x205 = x200 * x3
    x206 = x10 * x66
    x207 = -x206
    x208 = x176 * x42
    x209 = x207 + x208
    x210 = x209 * x6
    x211 = -x210
    x212 = x6 * (x205 + x211)
    x213 = x4 * (x183 + x187 - x189)
    x214 = x193 * x3
    x215 = x203 * x6
    x216 = x213 + x214 - x215
    x217 = x4 * (x193 - x198 + x201)
    x218 = x4 * (x185 - x186 + x191)
    x219 = -x213
    x220 = x183 * x5
    x221 = x180 * x5 + x184
    x222 = x192 + x220
    x223 = x181 + x221 * x3 - x222 * x6
    x224 = x222 * x3
    x225 = x190 * x5
    x226 = x202 + x225
    x227 = x226 * x6
    x228 = x213 + x224 - x227
    x229 = x218 + x223 * x3 - x228 * x6 + x4 * (x191 - x220 + x221)
    x230 = x4 * (x201 + x222 - x225)
    x231 = x228 * x3
    x232 = x226 * x3
    x233 = x200 * x5
    x234 = x211 + x233
    x235 = x234 * x6
    x236 = x197 + x232 - x235
    x237 = x236 * x6
    x238 = x217 + x230 + x231 - x237
    x239 = -x197
    x240 = x4 * (x200 + x206 - x208)
    x241 = -x10 * x94
    x242 = x176 * x66 + x241
    x243 = x108 * x115
    x244 = -x115 * x8 + x23
    x245 = x108 * x112 + x244
    x246 = x122 * x8
    x247 = x246 + x33
    x248 = x4 * (-x243 + x245 + x247)
    x249 = -x246 + x32
    x250 = x243 + x249
    x251 = x245 * x3 - x250 * x6
    x252 = x250 * x3
    x253 = x108 * x122
    x254 = x132 * x8
    x255 = -x254 + x41
    x256 = x253 + x255
    x257 = x256 * x6
    x258 = x252 - x257
    x259 = x248 + x251 * x3 - x258 * x6
    x260 = x254 + x59
    x261 = x4 * (x250 - x253 + x260)
    x262 = x258 * x3
    x263 = x256 * x3
    x264 = x108 * x132
    x265 = x141 * x8
    x266 = -x265 + x65
    x267 = x264 + x266
    x268 = x267 * x6
    x269 = x263 - x268
    x270 = x269 * x6
    x271 = x261 + x262 - x270
    x272 = x265 + x91
    x273 = x4 * (x256 - x264 + x272)
    x274 = -x174 * x8 + x93
    x275 = x108 * x141 + x274
    x276 = -x261
    x277 = x190 * x8
    x278 = x108 * x183
    x279 = -x183 * x8
    x280 = x108 * x180 + x279
    x281 = x4 * (x277 - x278 + x280)
    x282 = -x277
    x283 = x278 + x282
    x284 = x280 * x3 - x283 * x6
    x285 = x283 * x3
    x286 = x108 * x190
    x287 = x200 * x8
    x288 = -x287
    x289 = x286 + x288
    x290 = x289 * x6
    x291 = x285 - x290
    x292 = x281 + x284 * x3 - x291 * x6
    x293 = x4 * (x283 - x286 + x287)
    x294 = x291 * x3
    x295 = x289 * x3
    x296 = x108 * x200
    x297 = x209 * x8
    x298 = -x297
    x299 = x296 + x298
    x300 = x299 * x6
    x301 = x295 - x300
    x302 = x301 * x6
    x303 = x293 + x294 - x302
    x304 = x4 * (x289 - x296 + x297)
    x305 = -x242 * x8
    x306 = x108 * x209 + x305
    x307 = -x293
    x308 = x176 * x183
    x309 = -x10 * x183 + x23
    x310 = x176 * x180 + x309
    x311 = x10 * x190
    x312 = x311 + x33
    x313 = x4 * (-x308 + x310 + x312)
    x314 = -x311 + x32
    x315 = x308 + x314
    x316 = x3 * x310 - x315 * x6
    x317 = x3 * x315
    x318 = x176 * x190
    x319 = x10 * x200
    x320 = -x319 + x41
    x321 = x318 + x320
    x322 = x321 * x6
    x323 = x317 - x322
    x324 = x3 * x316 + x313 - x323 * x6
    x325 = x319 + x59
    x326 = x4 * (x315 - x318 + x325)
    x327 = x3 * x323
    x328 = x3 * x321
    x329 = x176 * x200
    x330 = x10 * x209
    x331 = -x330 + x65
    x332 = x329 + x331
    x333 = x332 * x6
    x334 = x328 - x333
    x335 = x334 * x6
    x336 = x326 + x327 - x335
    x337 = x330 + x91
    x338 = x4 * (x321 - x329 + x337)
    x339 = -x10 * x242 + x93
    x340 = x176 * x209 + x339
    x341 = -x326
    x342 = -x7 - A[1]
    x343 = x16 * x342
    x344 = x114 + x343
    x345 = x344 * x5
    x346 = x13 * x342
    x347 = x120 + x346
    x348 = x347 * x6
    x349 = -x348
    x350 = x345 + x349
    x351 = x350 * x5
    x352 = x111 + x21 * x342
    x353 = -x344 * x6
    x354 = x352 * x5 + x353
    x355 = x4 * (x110 - x343 + x352)
    x356 = -x350 * x6 + x355
    x357 = x354 * x5 + x356
    x358 = x347 * x5
    x359 = x26 * x342
    x360 = x131 + x359
    x361 = x360 * x6
    x362 = x358 - x361
    x363 = x362 * x6
    x364 = x4 * (x119 + x344 - x346)
    x365 = -x364
    x366 = x363 + x365
    x367 = -x363 + x364
    x368 = x351 + x367
    x369 = x4 * (-x345 + x348 + x354)
    x370 = x3 * x357 - x368 * x6 + 2.0 * x369
    x371 = x3 * x368
    x372 = x362 * x5
    x373 = x4 * (x128 + x347 - x359)
    x374 = x360 * x5
    x375 = x342 * x42
    x376 = x139 + x375
    x377 = x376 * x6
    x378 = x374 - x377
    x379 = x378 * x6
    x380 = x373 - x379
    x381 = x372 + x380
    x382 = x381 * x6
    x383 = x4 * (x350 - x358 + x361)
    x384 = 2.0 * x383
    x385 = x371 - x382 + x384
    x386 = x3 * x350
    x387 = x3 * x354 + x356
    x388 = -x373
    x389 = x379 + x388
    x390 = x4 * (x138 + x360 - x375)
    x391 = x173 + x342 * x66
    x392 = x3 * x362
    x393 = x367 + x386
    x394 = x3 * x344
    x395 = 26.9847693667702 * x106
    x396 = x115 * x342
    x397 = x112 * x342 + x244
    x398 = x4 * (x247 - x396 + x397)
    x399 = x249 + x396
    x400 = -x399 * x6
    x401 = x3 * x397 + x400
    x402 = x122 * x342
    x403 = x255 + x402
    x404 = x3 * x403
    x405 = x132 * x342
    x406 = x266 + x405
    x407 = x406 * x6
    x408 = -x407
    x409 = x4 * (x260 + x399 - x402)
    x410 = -x409
    x411 = x3 * x399
    x412 = x403 * x6
    x413 = -x412
    x414 = x411 + x413
    x415 = x399 * x5
    x416 = x397 * x5 + x400
    x417 = x413 + x415
    x418 = x3 * x416 + x398 - x417 * x6
    x419 = x3 * x417
    x420 = x403 * x5
    x421 = x408 + x420
    x422 = x421 * x6
    x423 = x409 + x419 - x422
    x424 = x4 * (x272 + x403 - x405)
    x425 = x141 * x342 + x274
    x426 = 46.7389915737742 * x106
    x427 = x183 * x342
    x428 = x180 * x342 + x279
    x429 = x4 * (x277 - x427 + x428)
    x430 = x282 + x427
    x431 = -x430 * x6
    x432 = x3 * x428 + x431
    x433 = x190 * x342
    x434 = x288 + x433
    x435 = x3 * x434
    x436 = x200 * x342
    x437 = x298 + x436
    x438 = x437 * x6
    x439 = -x438
    x440 = x4 * (x287 + x430 - x433)
    x441 = -x440
    x442 = x3 * x430
    x443 = x434 * x6
    x444 = -x443
    x445 = x442 + x444
    x446 = x430 * x5
    x447 = x428 * x5 + x431
    x448 = x444 + x446
    x449 = x3 * x447 + x429 - x448 * x6
    x450 = x3 * x448
    x451 = x434 * x5
    x452 = x439 + x451
    x453 = x452 * x6
    x454 = x440 + x450 - x453
    x455 = x4 * (x297 + x434 - x436)
    x456 = x209 * x342 + x305
    x457 = x256 * x8
    x458 = x250 * x342
    x459 = 2.0 * x145
    x460 = 2.0 * x113 + x245 * x342 - x250 * x8
    x461 = x4 * (x457 - x458 - x459 + x460)
    x462 = -x457 + x458 + x459
    x463 = x3 * x460 - x462 * x6
    x464 = x3 * x462
    x465 = x256 * x342
    x466 = x267 * x8
    x467 = 2.0 * x129
    x468 = x465 - x466 + x467
    x469 = x468 * x6
    x470 = x464 - x469
    x471 = x4 * (x462 - x465 + x466 - x467)
    x472 = 2.0 * x172 + x267 * x342 - x275 * x8
    x473 = x289 * x8
    x474 = x283 * x342
    x475 = x181 + x280 * x342 - x283 * x8
    x476 = x4 * (x219 + x473 - x474 + x475)
    x477 = x213 - x473 + x474
    x478 = x3 * x475 - x477 * x6
    x479 = x3 * x477
    x480 = x289 * x342
    x481 = x299 * x8
    x482 = x197 + x480 - x481
    x483 = x482 * x6
    x484 = x479 - x483
    x485 = x4 * (x239 + x477 - x480 + x481)
    x486 = x240 + x299 * x342 - x306 * x8
    x487 = x321 * x8
    x488 = x315 * x342
    x489 = x310 * x342 - x315 * x8
    x490 = x4 * (x487 - x488 + x489)
    x491 = -x487 + x488
    x492 = x3 * x489 - x491 * x6
    x493 = x3 * x491
    x494 = x321 * x342
    x495 = x332 * x8
    x496 = x494 - x495
    x497 = x496 * x6
    x498 = x493 - x497
    x499 = x4 * (x491 - x494 + x495)
    x500 = x332 * x342 - x340 * x8
    x501 = -x9 - A[2]
    x502 = x16 * x501
    x503 = x182 + x502
    x504 = x5 * x503
    x505 = x13 * x501
    x506 = x188 + x505
    x507 = x506 * x6
    x508 = -x507
    x509 = x504 + x508
    x510 = x5 * x509
    x511 = x179 + x21 * x501
    x512 = -x503 * x6
    x513 = x5 * x511 + x512
    x514 = x4 * (x178 - x502 + x511)
    x515 = -x509 * x6 + x514
    x516 = x5 * x513 + x515
    x517 = x5 * x506
    x518 = x26 * x501
    x519 = x199 + x518
    x520 = x519 * x6
    x521 = x517 - x520
    x522 = x521 * x6
    x523 = x4 * (x187 + x503 - x505)
    x524 = -x523
    x525 = x522 + x524
    x526 = -x522 + x523
    x527 = x510 + x526
    x528 = x4 * (-x504 + x507 + x513)
    x529 = x3 * x516 - x527 * x6 + 2.0 * x528
    x530 = x3 * x527
    x531 = x5 * x521
    x532 = x4 * (x196 + x506 - x518)
    x533 = x5 * x519
    x534 = x42 * x501
    x535 = x207 + x534
    x536 = x535 * x6
    x537 = x533 - x536
    x538 = x537 * x6
    x539 = x532 - x538
    x540 = x531 + x539
    x541 = x540 * x6
    x542 = x4 * (x509 - x517 + x520)
    x543 = 2.0 * x542
    x544 = x530 - x541 + x543
    x545 = x3 * x509
    x546 = x3 * x513 + x515
    x547 = -x532
    x548 = x538 + x547
    x549 = x4 * (x206 + x519 - x534)
    x550 = x241 + x501 * x66
    x551 = x3 * x521
    x552 = x526 + x545
    x553 = x3 * x503
    x554 = x506 * x8
    x555 = x108 * x503
    x556 = -x503 * x8
    x557 = x108 * x511 + x556
    x558 = x4 * (x554 - x555 + x557)
    x559 = -x554
    x560 = x555 + x559
    x561 = -x560 * x6
    x562 = x3 * x557 + x561
    x563 = x108 * x506
    x564 = x519 * x8
    x565 = -x564
    x566 = x563 + x565
    x567 = x3 * x566
    x568 = x108 * x519
    x569 = x535 * x8
    x570 = -x569
    x571 = x568 + x570
    x572 = x571 * x6
    x573 = -x572
    x574 = x4 * (x560 - x563 + x564)
    x575 = -x574
    x576 = x3 * x560
    x577 = x566 * x6
    x578 = -x577
    x579 = x576 + x578
    x580 = x5 * x560
    x581 = x5 * x557 + x561
    x582 = x578 + x580
    x583 = x3 * x581 + x558 - x582 * x6
    x584 = x3 * x582
    x585 = x5 * x566
    x586 = x573 + x585
    x587 = x586 * x6
    x588 = x574 + x584 - x587
    x589 = x4 * (x566 - x568 + x569)
    x590 = -x550 * x8
    x591 = x108 * x535 + x590
    x592 = x183 * x501
    x593 = x180 * x501 + x309
    x594 = x4 * (x312 - x592 + x593)
    x595 = x314 + x592
    x596 = -x595 * x6
    x597 = x3 * x593 + x596
    x598 = x190 * x501
    x599 = x320 + x598
    x600 = x3 * x599
    x601 = x200 * x501
    x602 = x331 + x601
    x603 = x6 * x602
    x604 = -x603
    x605 = x4 * (x325 + x595 - x598)
    x606 = -x605
    x607 = x3 * x595
    x608 = x599 * x6
    x609 = -x608
    x610 = x607 + x609
    x611 = x5 * x595
    x612 = x5 * x593 + x596
    x613 = x609 + x611
    x614 = x3 * x612 + x594 - x6 * x613
    x615 = x3 * x613
    x616 = x5 * x599
    x617 = x604 + x616
    x618 = x6 * x617
    x619 = x605 + x615 - x618
    x620 = x4 * (x337 + x599 - x601)
    x621 = x209 * x501 + x339
    x622 = x108 * x560
    x623 = x514 - x560 * x8
    x624 = x108 * x557 + x623
    x625 = x566 * x8
    x626 = x524 + x625
    x627 = x4 * (-x622 + x624 + x626)
    x628 = x523 - x625
    x629 = x622 + x628
    x630 = x3 * x624 - x6 * x629
    x631 = x3 * x629
    x632 = x108 * x566
    x633 = x571 * x8
    x634 = x532 - x633
    x635 = x632 + x634
    x636 = x6 * x635
    x637 = x631 - x636
    x638 = x547 + x633
    x639 = x4 * (x629 - x632 + x638)
    x640 = x549 - x591 * x8
    x641 = x108 * x571 + x640
    x642 = x599 * x8
    x643 = x108 * x595
    x644 = -x595 * x8
    x645 = x108 * x593 + x644
    x646 = x4 * (x642 - x643 + x645)
    x647 = -x642
    x648 = x643 + x647
    x649 = x3 * x645 - x6 * x648
    x650 = x3 * x648
    x651 = x108 * x599
    x652 = x602 * x8
    x653 = -x652
    x654 = x651 + x653
    x655 = x6 * x654
    x656 = x650 - x655
    x657 = x4 * (x648 - x651 + x652)
    x658 = -x621 * x8
    x659 = x108 * x602 + x658
    x660 = x10 * x321
    x661 = x315 * x501
    x662 = 2.0 * x213
    x663 = -x10 * x315 + 2.0 * x181 + x310 * x501
    x664 = x4 * (x660 - x661 - x662 + x663)
    x665 = -x660 + x661 + x662
    x666 = x3 * x663 - x6 * x665
    x667 = x3 * x665
    x668 = x321 * x501
    x669 = x10 * x332
    x670 = 2.0 * x197
    x671 = x668 - x669 + x670
    x672 = x6 * x671
    x673 = x667 - x672
    x674 = x4 * (x665 - x668 + x669 - x670)
    x675 = -x10 * x340 + 2.0 * x240 + x332 * x501
    x676 = x342 * x344
    x677 = x347 * x8
    x678 = x32 + x676 - x677
    x679 = x5 * x678
    x680 = x342 * x347
    x681 = x360 * x8
    x682 = x41 + x680 - x681
    x683 = x6 * x682
    x684 = x679 - x683
    x685 = x5 * x684
    x686 = x23 + x342 * x352 - x344 * x8
    x687 = x5 * x686 - x6 * x678
    x688 = x4 * (x33 - x676 + x677 + x686)
    x689 = -x6 * x684 + x688
    x690 = x5 * x687 + x689
    x691 = x5 * x682
    x692 = x342 * x360
    x693 = x376 * x8
    x694 = x65 + x692 - x693
    x695 = x6 * x694
    x696 = x691 - x695
    x697 = x6 * x696
    x698 = x4 * (x59 + x678 - x680 + x681)
    x699 = -x698
    x700 = x697 + x699
    x701 = x685 - x697 + x698
    x702 = x4 * (x682 - x692 + x693 + x91)
    x703 = x342 * x376 - x391 * x8 + x93
    x704 = 34.8371874529163 * x106
    x705 = x113 + x342 * x397 + x355 - x399 * x8
    x706 = x342 * x403
    x707 = x406 * x8
    x708 = x129 + x373 + x706 - x707
    x709 = x6 * x708
    x710 = x342 * x399
    x711 = x403 * x8
    x712 = x145 + x364 + x710 - x711
    x713 = -x6 * x712
    x714 = x5 * x712
    x715 = x5 * x705 + x713
    x716 = x4 * (x151 + x365 + x705 - x710 + x711)
    x717 = -x709 + x714
    x718 = x4 * (x171 + x388 - x706 + x707 + x712)
    x719 = x172 + x342 * x406 + x390 - x425 * x8
    x720 = 60.3397786612521 * x106
    x721 = x181 + x342 * x428 - x430 * x8
    x722 = x342 * x434
    x723 = x437 * x8
    x724 = x197 + x722 - x723
    x725 = x6 * x724
    x726 = x342 * x430
    x727 = x434 * x8
    x728 = x213 + x726 - x727
    x729 = -x6 * x728
    x730 = x5 * x728
    x731 = x5 * x721 + x729
    x732 = x4 * (x219 + x721 - x726 + x727)
    x733 = -x725 + x730
    x734 = x4 * (x239 - x722 + x723 + x728)
    x735 = x240 + x342 * x437 - x456 * x8
    x736 = x468 * x8
    x737 = x342 * x462
    x738 = 2.0 * x409
    x739 = -x738
    x740 = 2.0 * x398
    x741 = x248 + x342 * x460 - x462 * x8 + x740
    x742 = x4 * (x276 + x736 - x737 + x739 + x741)
    x743 = x261 - x736 + x737 + x738
    x744 = 2.0 * x424
    x745 = x273 + x342 * x468 - x472 * x8 + x744
    x746 = x482 * x8
    x747 = x342 * x477
    x748 = x281 + x342 * x475 + x429 - x477 * x8
    x749 = x4 * (x307 + x441 + x746 - x747 + x748)
    x750 = x293 + x440 - x746 + x747
    x751 = x304 + x342 * x482 + x455 - x486 * x8
    x752 = x496 * x8
    x753 = x342 * x491
    x754 = x313 + x342 * x489 - x491 * x8
    x755 = x4 * (x341 + x752 - x753 + x754)
    x756 = x326 - x752 + x753
    x757 = x338 + x342 * x496 - x500 * x8
    x758 = x342 * x503
    x759 = x559 + x758
    x760 = x5 * x759
    x761 = x342 * x506
    x762 = x565 + x761
    x763 = x6 * x762
    x764 = x760 - x763
    x765 = x5 * x764
    x766 = x342 * x511 + x556
    x767 = x5 * x766 - x6 * x759
    x768 = x4 * (x554 - x758 + x766)
    x769 = -x6 * x764 + x768
    x770 = x5 * x767 + x769
    x771 = x5 * x762
    x772 = x342 * x519
    x773 = x570 + x772
    x774 = x6 * x773
    x775 = x771 - x774
    x776 = x6 * x775
    x777 = x4 * (x564 + x759 - x761)
    x778 = -x777
    x779 = x776 + x778
    x780 = x765 - x776 + x777
    x781 = x4 * (x569 + x762 - x772)
    x782 = x342 * x535 + x590
    x783 = 60.3397786612521 * x106
    x784 = x342 * x557 + x623
    x785 = x342 * x566
    x786 = x634 + x785
    x787 = x6 * x786
    x788 = x342 * x560
    x789 = x628 + x788
    x790 = -x6 * x789
    x791 = x5 * x789
    x792 = x5 * x784 + x790
    x793 = x4 * (x626 + x784 - x788)
    x794 = -x787 + x791
    x795 = x4 * (x638 - x785 + x789)
    x796 = x342 * x571 + x640
    x797 = 104.511562358749 * x106
    x798 = x342 * x593 + x644
    x799 = x342 * x599
    x800 = x653 + x799
    x801 = x6 * x800
    x802 = x342 * x595
    x803 = x647 + x802
    x804 = -x6 * x803
    x805 = x5 * x803
    x806 = x5 * x798 + x804
    x807 = x4 * (x642 + x798 - x802)
    x808 = -x801 + x805
    x809 = x4 * (x652 - x799 + x803)
    x810 = x342 * x602 + x658
    x811 = x635 * x8
    x812 = x342 * x629
    x813 = 2.0 * x574
    x814 = x342 * x624 + 2.0 * x558 - x629 * x8
    x815 = x4 * (x811 - x812 - x813 + x814)
    x816 = -x811 + x812 + x813
    x817 = x342 * x635 + 2.0 * x589 - x641 * x8
    x818 = x654 * x8
    x819 = x342 * x648
    x820 = x342 * x645 + x594 - x648 * x8
    x821 = x4 * (x606 + x818 - x819 + x820)
    x822 = x605 - x818 + x819
    x823 = x342 * x654 + x620 - x659 * x8
    x824 = x671 * x8
    x825 = x342 * x665
    x826 = x342 * x663 - x665 * x8
    x827 = x4 * (x824 - x825 + x826)
    x828 = -x824 + x825
    x829 = x342 * x671 - x675 * x8
    x830 = x501 * x503
    x831 = x10 * x506
    x832 = x32 + x830 - x831
    x833 = x5 * x832
    x834 = x501 * x506
    x835 = x10 * x519
    x836 = x41 + x834 - x835
    x837 = x6 * x836
    x838 = x833 - x837
    x839 = x5 * x838
    x840 = -x10 * x503 + x23 + x501 * x511
    x841 = x5 * x840 - x6 * x832
    x842 = x4 * (x33 - x830 + x831 + x840)
    x843 = -x6 * x838 + x842
    x844 = x5 * x841 + x843
    x845 = x5 * x836
    x846 = x501 * x519
    x847 = x10 * x535
    x848 = x65 + x846 - x847
    x849 = x6 * x848
    x850 = x845 - x849
    x851 = x6 * x850
    x852 = x4 * (x59 + x832 - x834 + x835)
    x853 = -x852
    x854 = x851 + x853
    x855 = x839 - x851 + x852
    x856 = x4 * (x836 - x846 + x847 + x91)
    x857 = -x10 * x550 + x501 * x535 + x93
    x858 = -x8 * x832
    x859 = x108 * x840 + x858
    x860 = x108 * x836
    x861 = x8 * x848
    x862 = -x861
    x863 = x860 + x862
    x864 = x6 * x863
    x865 = x108 * x832
    x866 = x8 * x836
    x867 = -x866
    x868 = x865 + x867
    x869 = -x6 * x868
    x870 = x5 * x868
    x871 = x5 * x859 + x869
    x872 = x4 * (x859 - x865 + x866)
    x873 = -x864 + x870
    x874 = x4 * (-x860 + x861 + x868)
    x875 = -x8 * x857
    x876 = x108 * x848 + x875
    x877 = -x10 * x595 + x181 + x501 * x593 + x514
    x878 = x501 * x599
    x879 = x10 * x602
    x880 = x197 + x532 + x878 - x879
    x881 = x6 * x880
    x882 = x501 * x595
    x883 = x10 * x599
    x884 = x213 + x523 + x882 - x883
    x885 = -x6 * x884
    x886 = x5 * x884
    x887 = x5 * x877 + x885
    x888 = x4 * (x219 + x524 + x877 - x882 + x883)
    x889 = -x881 + x886
    x890 = x4 * (x239 + x547 - x878 + x879 + x884)
    x891 = -x10 * x621 + x240 + x501 * x602 + x549
    x892 = x108 * x868
    x893 = -x8 * x868 + x842
    x894 = x108 * x859 + x893
    x895 = x8 * x863
    x896 = x853 + x895
    x897 = x4 * (-x892 + x894 + x896)
    x898 = x852 - x895
    x899 = x892 + x898
    x900 = -x8 * x876 + x856
    x901 = x108 * x863 + x900
    x902 = x8 * x880
    x903 = x108 * x884
    x904 = -x8 * x884
    x905 = x108 * x877 + x904
    x906 = x4 * (x902 - x903 + x905)
    x907 = -x902
    x908 = x903 + x907
    x909 = -x8 * x891
    x910 = x108 * x880 + x909
    x911 = x10 * x671
    x912 = x501 * x665
    x913 = 2.0 * x605
    x914 = -x913
    x915 = 2.0 * x594
    x916 = -x10 * x665 + x313 + x501 * x663 + x915
    x917 = x4 * (x341 + x911 - x912 + x914 + x916)
    x918 = x326 - x911 + x912 + x913
    x919 = 2.0 * x620
    x920 = -x10 * x675 + x338 + x501 * x671 + x919
    x921 = x682 * x8
    x922 = x342 * x678
    x923 = 2.0 * x364
    x924 = x342 * x686 + 2.0 * x355 - x678 * x8
    x925 = x4 * (x921 - x922 - x923 + x924)
    x926 = -x921 + x922 + x923
    x927 = x5 * x924 - x6 * x926
    x928 = x5 * x926
    x929 = x342 * x682
    x930 = x694 * x8
    x931 = 2.0 * x373
    x932 = x929 - x930 + x931
    x933 = x6 * x932
    x934 = x928 - x933
    x935 = x4 * (x926 - x929 + x930 - x931)
    x936 = x342 * x694 + 2.0 * x390 - x703 * x8
    x937 = x708 * x8
    x938 = x342 * x712
    x939 = x342 * x705 + x688 - x712 * x8 + x740
    x940 = x4 * (x699 + x739 + x937 - x938 + x939)
    x941 = x698 + x738 - x937 + x938
    x942 = x342 * x708 + x702 - x719 * x8 + x744
    x943 = x724 * x8
    x944 = x342 * x728
    x945 = 2.0 * x440
    x946 = x342 * x721 + 2.0 * x429 - x728 * x8
    x947 = x4 * (x943 - x944 - x945 + x946)
    x948 = -x943 + x944 + x945
    x949 = x342 * x724 + 2.0 * x455 - x735 * x8
    x950 = x342 * x741 + 2.0 * x461 + 2.0 * x716 - x743 * x8
    x951 = x342 * x743 + 2.0 * x471 + 2.0 * x718 - x745 * x8
    x952 = x342 * x748 + 2.0 * x476 + x732 - x750 * x8
    x953 = x342 * x750 + 2.0 * x485 + x734 - x751 * x8
    x954 = x342 * x754 + 2.0 * x490 - x756 * x8
    x955 = x342 * x756 + 2.0 * x499 - x757 * x8
    x956 = x762 * x8
    x957 = x342 * x759
    x958 = x342 * x766 + x514 - x759 * x8
    x959 = x4 * (x524 + x956 - x957 + x958)
    x960 = x523 - x956 + x957
    x961 = x5 * x958 - x6 * x960
    x962 = x5 * x960
    x963 = x342 * x762
    x964 = x773 * x8
    x965 = x532 + x963 - x964
    x966 = x6 * x965
    x967 = x962 - x966
    x968 = x4 * (x547 + x960 - x963 + x964)
    x969 = x342 * x773 + x549 - x782 * x8
    x970 = x786 * x8
    x971 = x342 * x789
    x972 = x342 * x784 + x558 + x768 - x789 * x8
    x973 = x4 * (x575 + x778 + x970 - x971 + x972)
    x974 = x574 + x777 - x970 + x971
    x975 = x342 * x786 + x589 + x781 - x796 * x8
    x976 = x8 * x800
    x977 = x342 * x803
    x978 = x342 * x798 + x594 - x8 * x803
    x979 = x4 * (x606 + x976 - x977 + x978)
    x980 = x605 - x976 + x977
    x981 = x342 * x800 + x620 - x8 * x810
    x982 = 2.0 * x793
    x983 = x342 * x814 + x627 - x8 * x816 + x982
    x984 = 2.0 * x795
    x985 = x342 * x816 + x639 - x8 * x817 + x984
    x986 = x342 * x820 + x646 - x8 * x822 + x807
    x987 = x342 * x822 + x657 - x8 * x823 + x809
    x988 = x342 * x826 + x664 - x8 * x828
    x989 = x342 * x828 + x674 - x8 * x829
    x990 = x342 * x832
    x991 = x342 * x840 + x858
    x992 = x4 * (x866 - x990 + x991)
    x993 = x867 + x990
    x994 = x5 * x991 - x6 * x993
    x995 = x5 * x993
    x996 = x342 * x836
    x997 = x862 + x996
    x998 = x6 * x997
    x999 = x995 - x998
    x1000 = x4 * (x861 + x993 - x996)
    x1001 = x342 * x848 + x875
    x1002 = x342 * x868
    x1003 = x342 * x859 + x893
    x1004 = x4 * (-x1002 + x1003 + x896)
    x1005 = x1002 + x898
    x1006 = x342 * x863 + x900
    x1007 = x342 * x884
    x1008 = x342 * x877 + x904
    x1009 = x4 * (-x1007 + x1008 + x902)
    x1010 = x1007 + x907
    x1011 = x342 * x880 + x909
    x1012 = x342 * x894 - x8 * x899 + 2.0 * x872
    x1013 = x342 * x899 - x8 * x901 + 2.0 * x874
    x1014 = x342 * x905 - x8 * x908 + x888
    x1015 = x342 * x908 - x8 * x910 + x890
    x1016 = x342 * x916 - x8 * x918
    x1017 = x342 * x918 - x8 * x920
    x1018 = x10 * x836
    x1019 = x501 * x832
    x1020 = 2.0 * x523
    x1021 = -x10 * x832 + x501 * x840 + 2.0 * x514
    x1022 = x4 * (x1018 - x1019 - x1020 + x1021)
    x1023 = -x1018 + x1019 + x1020
    x1024 = x1021 * x5 - x1023 * x6
    x1025 = x1023 * x5
    x1026 = x501 * x836
    x1027 = x10 * x848
    x1028 = 2.0 * x532
    x1029 = x1026 - x1027 + x1028
    x1030 = x1029 * x6
    x1031 = x1025 - x1030
    x1032 = x4 * (x1023 - x1026 + x1027 - x1028)
    x1033 = -x10 * x857 + x501 * x848 + 2.0 * x549
    x1034 = x1029 * x8
    x1035 = x1023 * x108
    x1036 = -x1023 * x8
    x1037 = x1021 * x108 + x1036
    x1038 = x4 * (x1034 - x1035 + x1037)
    x1039 = -x1034
    x1040 = x1035 + x1039
    x1041 = -x1033 * x8
    x1042 = x1029 * x108 + x1041
    x1043 = x10 * x880
    x1044 = x501 * x884
    x1045 = -x10 * x884 + x501 * x877 + x842 + x915
    x1046 = x4 * (x1043 - x1044 + x1045 + x853 + x914)
    x1047 = -x1043 + x1044 + x852 + x913
    x1048 = -x10 * x891 + x501 * x880 + x856 + x919
    x1049 = x1022 - x1040 * x8
    x1050 = x1037 * x108 + x1049
    x1051 = x1032 - x1042 * x8
    x1052 = x1040 * x108 + x1051
    x1053 = -x1047 * x8
    x1054 = x1045 * x108 + x1053
    x1055 = -x1048 * x8
    x1056 = x1047 * x108 + x1055
    x1057 = -x10 * x918 + x501 * x916 + 2.0 * x664 + 2.0 * x888
    x1058 = -x10 * x920 + x501 * x918 + 2.0 * x674 + 2.0 * x890
    x1059 = x8 * x932
    x1060 = x342 * x926
    x1061 = 3.0 * x698
    x1062 = x342 * x924 + 3.0 * x688 - x8 * x926
    x1063 = -x1059 + x1060 + x1061
    x1064 = x8 * x965
    x1065 = x342 * x960
    x1066 = 2.0 * x777
    x1067 = x342 * x958 + 2.0 * x768 - x8 * x960
    x1068 = -x1064 + x1065 + x1066
    x1069 = x8 * x997
    x1070 = x342 * x993
    x1071 = x342 * x991 - x8 * x993 + x842
    x1072 = -x1069 + x1070 + x852
    x1073 = x1023 * x342
    x1074 = x1021 * x342 + x1036
    x1075 = x1039 + x1073
    x1076 = x10 * x1029
    x1077 = x1023 * x501
    x1078 = 3.0 * x852
    x1079 = -x10 * x1023 + x1021 * x501 + 3.0 * x842
    x1080 = x4 * (x1076 - x1077 - x1078 + x1079)
    x1081 = -x1076 + x1077 + x1078
    x1082 = -x10 * x1033 + x1029 * x501 + 3.0 * x856
    x1083 = x1079 * x108 - x1081 * x8
    x1084 = x108 * x1081 - x1082 * x8
    x1085 = -x10 * x1047 + x1022 + x1045 * x501 + 3.0 * x888
    x1086 = -x10 * x1048 + x1032 + x1047 * x501 + 3.0 * x890

    # 90 item(s)
    result[0, 0] = numpy.sum(
        x107
        * (
            x105 * (x104 + x58 - x61 - x62 + x77)
            + x3
            * (
                x3 * x58
                + x56 * (x38 - x39 + x50 - x52)
                - x56 * (x51 - x83 + x87 + x88 - x90)
                - x6 * x81
            )
            + x56
            * (
                -x101 * x3
                - x101 * x6
                + x104
                + x3 * x90
                - x4 * (-x102 + x103 - x3 * x98 + x32 + x59 + x6 * (x44 + x97))
                + x4 * (x102 - x103 + x23 + x3 * x89 + x33 - x6 * x86)
                + x57
                + x6 * (x100 + x74 - x96 + x99)
            )
            - x6
            * (
                x3 * x81
                + x56 * (x53 - x63 + x73 - x75)
                - x56 * (x100 - x101 + x74 - x96 + x99)
                - x6
                * (
                    x3 * x76
                    + x4 * (x49 - x64 + x92)
                    + x56 * (x82 + x92 - x95)
                    - x6
                    * (
                        x3 * x72
                        + x56 * (x46 + x67 - x68)
                        - x6 * (x5 * x69 - x6 * (x5 * x66 - x6 * x94) + x93)
                    )
                )
            )
        )
    )
    result[0, 1] = numpy.sum(
        -x175
        * (
            x105 * (x149 - x161 + x162 + x163 - x169)
            - x3
            * (
                x161 * x3
                - x170 * x6
                + x4 * (x126 - x146 + x147 + x151)
                + x56 * (x151 + x155 - x156 + x159)
            )
            - x4
            * (
                x126 * x3
                - x148 * x3
                - x148 * x6
                - 2.0 * x149
                + 2.0 * x150
                + x6 * (x129 + x136 - x144)
            )
            + x6
            * (
                x170 * x3
                + x4 * (-x136 + x144 + x148 + x171)
                + x56 * (x160 - x164 + x167 + x171)
                - x6
                * (
                    x168 * x3
                    + x4 * (x135 - x137 + x142)
                    + x4 * (x142 + x158 - x165)
                    - x6 * (x166 * x3 + x172 - x6 * (x141 * x5 - x174 * x6))
                )
            )
        )
    )
    result[0, 2] = numpy.sum(
        -x175
        * (
            x105 * (x217 - x229 + x230 + x231 - x237)
            - x3
            * (
                x229 * x3
                - x238 * x6
                + x4 * (x194 - x214 + x215 + x219)
                + x56 * (x219 + x223 - x224 + x227)
            )
            - x4
            * (
                x194 * x3
                - x216 * x3
                - x216 * x6
                - 2.0 * x217
                + 2.0 * x218
                + x6 * (x197 + x204 - x212)
            )
            + x6
            * (
                x238 * x3
                + x4 * (-x204 + x212 + x216 + x239)
                + x56 * (x228 - x232 + x235 + x239)
                - x6
                * (
                    x236 * x3
                    + x4 * (x203 - x205 + x210)
                    + x4 * (x210 + x226 - x233)
                    - x6 * (x234 * x3 + x240 - x6 * (x209 * x5 - x242 * x6))
                )
            )
        )
    )
    result[0, 3] = numpy.sum(
        x107
        * (
            x105 * (x259 - x262 + x270 + x276)
            + x3 * (x259 * x3 - x271 * x6 + x56 * (x251 - x252 + x257))
            - x6
            * (
                x271 * x3
                + x56 * (x258 - x263 + x268)
                - x6 * (x269 * x3 + x273 - x6 * (x267 * x3 - x275 * x6))
            )
        )
    )
    result[0, 4] = numpy.sum(
        x175
        * (
            x105 * (x292 - x294 + x302 + x307)
            + x3 * (x292 * x3 - x303 * x6 + x56 * (x284 - x285 + x290))
            - x6
            * (
                x3 * x303
                + x56 * (x291 - x295 + x300)
                - x6 * (x3 * x301 + x304 - x6 * (x299 * x3 - x306 * x6))
            )
        )
    )
    result[0, 5] = numpy.sum(
        x107
        * (
            x105 * (x324 - x327 + x335 + x341)
            + x3 * (x3 * x324 - x336 * x6 + x56 * (x316 - x317 + x322))
            - x6
            * (
                x3 * x336
                + x56 * (x323 - x328 + x333)
                - x6 * (x3 * x334 + x338 - x6 * (x3 * x332 - x340 * x6))
            )
        )
    )
    result[1, 0] = numpy.sum(
        x395
        * (
            x3
            * (
                x3 * x370
                - x385 * x6
                + x4 * (-x351 + x357 + x366)
                + x56 * (x366 - x386 + x387)
            )
            + x56 * (x370 - x371 + x382 - x384)
            + x56
            * (
                x3 * x387
                - x3 * x393
                + x369
                - x383
                - x393 * x6
                - x4 * (-x3 * x347 + x349 + x361 + x394)
                + x4 * (x3 * x352 + x348 + x353 - x394)
                + x6 * (x380 + x392)
            )
            - x6
            * (
                x3 * x385
                + x4 * (x368 - x372 + x389)
                + x56 * (x389 - x392 + x393)
                - x6
                * (
                    x3 * x381
                    + x56 * (x362 - x374 + x377)
                    - x6 * (x378 * x5 + x390 - x6 * (x376 * x5 - x391 * x6))
                )
            )
        )
    )
    result[1, 1] = numpy.sum(
        x426
        * (
            x3
            * (
                x3 * x418
                + x4 * (x401 - x411 + x412)
                + x4 * (x412 - x415 + x416)
                - x423 * x6
            )
            + x4 * (x3 * x401 - x3 * x414 + x398 + x410 - x414 * x6 + x6 * (x404 + x408))
            + x56 * (x410 + x418 - x419 + x422)
            - x6
            * (
                x3 * x423
                + x4 * (-x404 + x407 + x414)
                + x4 * (x407 + x417 - x420)
                - x6 * (x3 * x421 + x424 - x6 * (x406 * x5 - x425 * x6))
            )
        )
    )
    result[1, 2] = numpy.sum(
        x426
        * (
            x3
            * (
                x3 * x449
                + x4 * (x432 - x442 + x443)
                + x4 * (x443 - x446 + x447)
                - x454 * x6
            )
            + x4 * (x3 * x432 - x3 * x445 + x429 + x441 - x445 * x6 + x6 * (x435 + x439))
            + x56 * (x441 + x449 - x450 + x453)
            - x6
            * (
                x3 * x454
                + x4 * (-x435 + x438 + x445)
                + x4 * (x438 + x448 - x451)
                - x6 * (x3 * x452 + x455 - x6 * (x437 * x5 - x456 * x6))
            )
        )
    )
    result[1, 3] = numpy.sum(
        x395
        * (
            x3 * (x3 * x463 + x461 - x470 * x6)
            + x56 * (x463 - x464 + x469)
            - x6 * (x3 * x470 + x471 - x6 * (x3 * x468 - x472 * x6))
        )
    )
    result[1, 4] = numpy.sum(
        x426
        * (
            x3 * (x3 * x478 + x476 - x484 * x6)
            + x56 * (x478 - x479 + x483)
            - x6 * (x3 * x484 + x485 - x6 * (x3 * x482 - x486 * x6))
        )
    )
    result[1, 5] = numpy.sum(
        x395
        * (
            x3 * (x3 * x492 + x490 - x498 * x6)
            + x56 * (x492 - x493 + x497)
            - x6 * (x3 * x498 + x499 - x6 * (x3 * x496 - x500 * x6))
        )
    )
    result[2, 0] = numpy.sum(
        x395
        * (
            x3
            * (
                x3 * x529
                + x4 * (-x510 + x516 + x525)
                - x544 * x6
                + x56 * (x525 - x545 + x546)
            )
            + x56 * (x529 - x530 + x541 - x543)
            + x56
            * (
                x3 * x546
                - x3 * x552
                - x4 * (-x3 * x506 + x508 + x520 + x553)
                + x4 * (x3 * x511 + x507 + x512 - x553)
                + x528
                - x542
                - x552 * x6
                + x6 * (x539 + x551)
            )
            - x6
            * (
                x3 * x544
                + x4 * (x527 - x531 + x548)
                + x56 * (x548 - x551 + x552)
                - x6
                * (
                    x3 * x540
                    + x56 * (x521 - x533 + x536)
                    - x6 * (x5 * x537 + x549 - x6 * (x5 * x535 - x550 * x6))
                )
            )
        )
    )
    result[2, 1] = numpy.sum(
        x426
        * (
            x3
            * (
                x3 * x583
                + x4 * (x562 - x576 + x577)
                + x4 * (x577 - x580 + x581)
                - x588 * x6
            )
            + x4 * (x3 * x562 - x3 * x579 + x558 + x575 - x579 * x6 + x6 * (x567 + x573))
            + x56 * (x575 + x583 - x584 + x587)
            - x6
            * (
                x3 * x588
                + x4 * (-x567 + x572 + x579)
                + x4 * (x572 + x582 - x585)
                - x6 * (x3 * x586 + x589 - x6 * (x5 * x571 - x591 * x6))
            )
        )
    )
    result[2, 2] = numpy.sum(
        x426
        * (
            x3
            * (
                x3 * x614
                + x4 * (x597 - x607 + x608)
                + x4 * (x608 - x611 + x612)
                - x6 * x619
            )
            + x4 * (x3 * x597 - x3 * x610 + x594 - x6 * x610 + x6 * (x600 + x604) + x606)
            + x56 * (x606 + x614 - x615 + x618)
            - x6
            * (
                x3 * x619
                + x4 * (-x600 + x603 + x610)
                + x4 * (x603 + x613 - x616)
                - x6 * (x3 * x617 - x6 * (x5 * x602 - x6 * x621) + x620)
            )
        )
    )
    result[2, 3] = numpy.sum(
        x395
        * (
            x3 * (x3 * x630 - x6 * x637 + x627)
            + x56 * (x630 - x631 + x636)
            - x6 * (x3 * x637 - x6 * (x3 * x635 - x6 * x641) + x639)
        )
    )
    result[2, 4] = numpy.sum(
        x426
        * (
            x3 * (x3 * x649 - x6 * x656 + x646)
            + x56 * (x649 - x650 + x655)
            - x6 * (x3 * x656 - x6 * (x3 * x654 - x6 * x659) + x657)
        )
    )
    result[2, 5] = numpy.sum(
        x395
        * (
            x3 * (x3 * x666 - x6 * x673 + x664)
            + x56 * (x666 - x667 + x672)
            - x6 * (x3 * x673 - x6 * (x3 * x671 - x6 * x675) + x674)
        )
    )
    result[3, 0] = numpy.sum(
        x704
        * (
            x3 * (x3 * x690 + x56 * (-x679 + x683 + x687) - x6 * x701)
            + x4 * (-x685 + x690 + x700)
            + x56 * (-x3 * x684 + x3 * x687 + x689 + x700)
            - x6
            * (
                x3 * x701
                + x56 * (x684 - x691 + x695)
                - x6 * (x5 * x696 - x6 * (x5 * x694 - x6 * x703) + x702)
            )
        )
    )
    result[3, 1] = numpy.sum(
        x720
        * (
            x3 * (x3 * x715 - x6 * x717 + x716)
            + x4 * (x709 - x714 + x715)
            + x4 * (x3 * x705 - x3 * x712 + x709 + x713)
            - x6 * (x3 * x717 - x6 * (x5 * x708 - x6 * x719) + x718)
        )
    )
    result[3, 2] = numpy.sum(
        x720
        * (
            x3 * (x3 * x731 - x6 * x733 + x732)
            + x4 * (x725 - x730 + x731)
            + x4 * (x3 * x721 - x3 * x728 + x725 + x729)
            - x6 * (x3 * x733 - x6 * (x5 * x724 - x6 * x735) + x734)
        )
    )
    result[3, 3] = numpy.sum(
        x704 * (x3 * (x3 * x741 - x6 * x743) - x6 * (x3 * x743 - x6 * x745) + x742)
    )
    result[3, 4] = numpy.sum(
        x720 * (x3 * (x3 * x748 - x6 * x750) - x6 * (x3 * x750 - x6 * x751) + x749)
    )
    result[3, 5] = numpy.sum(
        x704 * (x3 * (x3 * x754 - x6 * x756) - x6 * (x3 * x756 - x6 * x757) + x755)
    )
    result[4, 0] = numpy.sum(
        x783
        * (
            x3 * (x3 * x770 + x56 * (-x760 + x763 + x767) - x6 * x780)
            + x4 * (-x765 + x770 + x779)
            + x56 * (-x3 * x764 + x3 * x767 + x769 + x779)
            - x6
            * (
                x3 * x780
                + x56 * (x764 - x771 + x774)
                - x6 * (x5 * x775 - x6 * (x5 * x773 - x6 * x782) + x781)
            )
        )
    )
    result[4, 1] = numpy.sum(
        x797
        * (
            x3 * (x3 * x792 - x6 * x794 + x793)
            + x4 * (x787 - x791 + x792)
            + x4 * (x3 * x784 - x3 * x789 + x787 + x790)
            - x6 * (x3 * x794 - x6 * (x5 * x786 - x6 * x796) + x795)
        )
    )
    result[4, 2] = numpy.sum(
        x797
        * (
            x3 * (x3 * x806 - x6 * x808 + x807)
            + x4 * (x801 - x805 + x806)
            + x4 * (x3 * x798 - x3 * x803 + x801 + x804)
            - x6 * (x3 * x808 - x6 * (x5 * x800 - x6 * x810) + x809)
        )
    )
    result[4, 3] = numpy.sum(
        x783 * (x3 * (x3 * x814 - x6 * x816) - x6 * (x3 * x816 - x6 * x817) + x815)
    )
    result[4, 4] = numpy.sum(
        x797 * (x3 * (x3 * x820 - x6 * x822) - x6 * (x3 * x822 - x6 * x823) + x821)
    )
    result[4, 5] = numpy.sum(
        x783 * (x3 * (x3 * x826 - x6 * x828) - x6 * (x3 * x828 - x6 * x829) + x827)
    )
    result[5, 0] = numpy.sum(
        x704
        * (
            x3 * (x3 * x844 + x56 * (-x833 + x837 + x841) - x6 * x855)
            + x4 * (-x839 + x844 + x854)
            + x56 * (-x3 * x838 + x3 * x841 + x843 + x854)
            - x6
            * (
                x3 * x855
                + x56 * (x838 - x845 + x849)
                - x6 * (x5 * x850 - x6 * (x5 * x848 - x6 * x857) + x856)
            )
        )
    )
    result[5, 1] = numpy.sum(
        x720
        * (
            x3 * (x3 * x871 - x6 * x873 + x872)
            + x4 * (x864 - x870 + x871)
            + x4 * (x3 * x859 - x3 * x868 + x864 + x869)
            - x6 * (x3 * x873 - x6 * (x5 * x863 - x6 * x876) + x874)
        )
    )
    result[5, 2] = numpy.sum(
        x720
        * (
            x3 * (x3 * x887 - x6 * x889 + x888)
            + x4 * (x881 - x886 + x887)
            + x4 * (x3 * x877 - x3 * x884 + x881 + x885)
            - x6 * (x3 * x889 - x6 * (x5 * x880 - x6 * x891) + x890)
        )
    )
    result[5, 3] = numpy.sum(
        x704 * (x3 * (x3 * x894 - x6 * x899) - x6 * (x3 * x899 - x6 * x901) + x897)
    )
    result[5, 4] = numpy.sum(
        x720 * (x3 * (x3 * x905 - x6 * x908) - x6 * (x3 * x908 - x6 * x910) + x906)
    )
    result[5, 5] = numpy.sum(
        x704 * (x3 * (x3 * x916 - x6 * x918) - x6 * (x3 * x918 - x6 * x920) + x917)
    )
    result[6, 0] = numpy.sum(
        x395
        * (
            x3 * (x5 * x927 - x6 * x934 + x925)
            + x56 * (x927 - x928 + x933)
            - x6 * (x5 * x934 - x6 * (x5 * x932 - x6 * x936) + x935)
        )
    )
    result[6, 1] = numpy.sum(
        x426 * (x3 * (x5 * x939 - x6 * x941) - x6 * (x5 * x941 - x6 * x942) + x940)
    )
    result[6, 2] = numpy.sum(
        x426 * (x3 * (x5 * x946 - x6 * x948) - x6 * (x5 * x948 - x6 * x949) + x947)
    )
    result[6, 3] = numpy.sum(x395 * (x3 * x950 - x6 * x951))
    result[6, 4] = numpy.sum(x426 * (x3 * x952 - x6 * x953))
    result[6, 5] = numpy.sum(x395 * (x3 * x954 - x6 * x955))
    result[7, 0] = numpy.sum(
        x783
        * (
            x3 * (x5 * x961 - x6 * x967 + x959)
            + x56 * (x961 - x962 + x966)
            - x6 * (x5 * x967 - x6 * (x5 * x965 - x6 * x969) + x968)
        )
    )
    result[7, 1] = numpy.sum(
        x797 * (x3 * (x5 * x972 - x6 * x974) - x6 * (x5 * x974 - x6 * x975) + x973)
    )
    result[7, 2] = numpy.sum(
        x797 * (x3 * (x5 * x978 - x6 * x980) - x6 * (x5 * x980 - x6 * x981) + x979)
    )
    result[7, 3] = numpy.sum(x783 * (x3 * x983 - x6 * x985))
    result[7, 4] = numpy.sum(x797 * (x3 * x986 - x6 * x987))
    result[7, 5] = numpy.sum(x783 * (x3 * x988 - x6 * x989))
    result[8, 0] = numpy.sum(
        x783
        * (
            x3 * (x5 * x994 - x6 * x999 + x992)
            + x56 * (x994 - x995 + x998)
            - x6 * (x1000 + x5 * x999 + x6 * (x1001 * x6 - x5 * x997))
        )
    )
    result[8, 1] = numpy.sum(
        x797 * (x1004 + x3 * (x1003 * x5 - x1005 * x6) - x6 * (x1005 * x5 - x1006 * x6))
    )
    result[8, 2] = numpy.sum(
        x797 * (x1009 + x3 * (x1008 * x5 - x1010 * x6) - x6 * (x1010 * x5 - x1011 * x6))
    )
    result[8, 3] = numpy.sum(x783 * (x1012 * x3 - x1013 * x6))
    result[8, 4] = numpy.sum(x797 * (x1014 * x3 - x1015 * x6))
    result[8, 5] = numpy.sum(x783 * (x1016 * x3 - x1017 * x6))
    result[9, 0] = numpy.sum(
        x395
        * (
            x3 * (x1022 + x1024 * x5 - x1031 * x6)
            + x56 * (x1024 - x1025 + x1030)
            - x6 * (x1031 * x5 + x1032 - x6 * (x1029 * x5 - x1033 * x6))
        )
    )
    result[9, 1] = numpy.sum(
        x426 * (x1038 + x3 * (x1037 * x5 - x1040 * x6) - x6 * (x1040 * x5 - x1042 * x6))
    )
    result[9, 2] = numpy.sum(
        x426 * (x1046 + x3 * (x1045 * x5 - x1047 * x6) - x6 * (x1047 * x5 - x1048 * x6))
    )
    result[9, 3] = numpy.sum(x395 * (x1050 * x3 - x1052 * x6))
    result[9, 4] = numpy.sum(x426 * (x1054 * x3 - x1056 * x6))
    result[9, 5] = numpy.sum(x395 * (x1057 * x3 - x1058 * x6))
    result[10, 0] = numpy.sum(
        x107
        * (
            x4 * (x1059 - x1060 - x1061 + x1062)
            + x5 * (x1062 * x5 - x1063 * x6)
            - x6 * (x1063 * x5 - x6 * (x342 * x932 + 3.0 * x702 - x8 * x936))
        )
    )
    result[10, 1] = numpy.sum(
        x175
        * (
            x5 * (x342 * x939 + 3.0 * x716 - x8 * x941 + x925)
            - x6 * (x342 * x941 + 3.0 * x718 - x8 * x942 + x935)
        )
    )
    result[10, 2] = numpy.sum(
        x175
        * (
            x5 * (x342 * x946 + 3.0 * x732 - x8 * x948)
            - x6 * (x342 * x948 + 3.0 * x734 - x8 * x949)
        )
    )
    result[10, 3] = numpy.sum(x107 * (x342 * x950 + 3.0 * x742 - x8 * x951 + 2.0 * x940))
    result[10, 4] = numpy.sum(x175 * (x342 * x952 + 3.0 * x749 - x8 * x953 + x947))
    result[10, 5] = numpy.sum(x107 * (x342 * x954 + 3.0 * x755 - x8 * x955))
    result[11, 0] = numpy.sum(
        x395
        * (
            x4 * (x1064 - x1065 - x1066 + x1067)
            + x5 * (x1067 * x5 - x1068 * x6)
            - x6 * (x1068 * x5 - x6 * (x342 * x965 + 2.0 * x781 - x8 * x969))
        )
    )
    result[11, 1] = numpy.sum(
        x426
        * (
            x5 * (x342 * x972 - x8 * x974 + x959 + x982)
            - x6 * (x342 * x974 - x8 * x975 + x968 + x984)
        )
    )
    result[11, 2] = numpy.sum(
        x426
        * (
            x5 * (x342 * x978 - x8 * x980 + 2.0 * x807)
            - x6 * (x342 * x980 - x8 * x981 + 2.0 * x809)
        )
    )
    result[11, 3] = numpy.sum(x395 * (x342 * x983 - x8 * x985 + 2.0 * x815 + 2.0 * x973))
    result[11, 4] = numpy.sum(x426 * (x342 * x986 - x8 * x987 + 2.0 * x821 + x979))
    result[11, 5] = numpy.sum(x395 * (x342 * x988 - x8 * x989 + 2.0 * x827))
    result[12, 0] = numpy.sum(
        x704
        * (
            x4 * (x1069 - x1070 + x1071 + x853)
            + x5 * (x1071 * x5 - x1072 * x6)
            - x6 * (x1072 * x5 - x6 * (-x1001 * x8 + x342 * x997 + x856))
        )
    )
    result[12, 1] = numpy.sum(
        x720
        * (
            x5 * (x1003 * x342 - x1005 * x8 + x872 + x992)
            - x6 * (x1000 + x1005 * x342 - x1006 * x8 + x874)
        )
    )
    result[12, 2] = numpy.sum(
        x720
        * (
            x5 * (x1008 * x342 - x1010 * x8 + x888)
            - x6 * (x1010 * x342 - x1011 * x8 + x890)
        )
    )
    result[12, 3] = numpy.sum(x704 * (2.0 * x1004 + x1012 * x342 - x1013 * x8 + x897))
    result[12, 4] = numpy.sum(x720 * (x1009 + x1014 * x342 - x1015 * x8 + x906))
    result[12, 5] = numpy.sum(x704 * (x1016 * x342 - x1017 * x8 + x917))
    result[13, 0] = numpy.sum(
        x395
        * (
            x4 * (x1034 - x1073 + x1074)
            + x5 * (x1074 * x5 - x1075 * x6)
            - x6 * (x1075 * x5 - x6 * (x1029 * x342 + x1041))
        )
    )
    result[13, 1] = numpy.sum(
        x426 * (x5 * (x1037 * x342 + x1049) - x6 * (x1040 * x342 + x1051))
    )
    result[13, 2] = numpy.sum(
        x426 * (x5 * (x1045 * x342 + x1053) - x6 * (x1047 * x342 + x1055))
    )
    result[13, 3] = numpy.sum(x395 * (2.0 * x1038 + x1050 * x342 - x1052 * x8))
    result[13, 4] = numpy.sum(x426 * (x1046 + x1054 * x342 - x1056 * x8))
    result[13, 5] = numpy.sum(x395 * (x1057 * x342 - x1058 * x8))
    result[14, 0] = numpy.sum(
        x107 * (x1080 + x5 * (x1079 * x5 - x1081 * x6) - x6 * (x1081 * x5 - x1082 * x6))
    )
    result[14, 1] = numpy.sum(x175 * (x1083 * x5 - x1084 * x6))
    result[14, 2] = numpy.sum(x175 * (x1085 * x5 - x1086 * x6))
    result[14, 3] = numpy.sum(x107 * (x108 * x1083 + x1080 - x1084 * x8))
    result[14, 4] = numpy.sum(x175 * (x108 * x1085 - x1086 * x8))
    result[14, 5] = numpy.sum(
        x107 * (-x10 * x1058 + 2.0 * x1046 + x1057 * x501 + 3.0 * x917)
    )
    return result


def coulomb3d_43(ax, da, A, bx, db, B, C):
    """Cartesian (gf) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 10), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = 0.5 / (ax + bx)
    x5 = -x2 - B[0]
    x6 = -x2 - C[0]
    x7 = -x1 * (ax * A[1] + bx * B[1])
    x8 = -x7 - C[1]
    x9 = -x1 * (ax * A[2] + bx * B[2])
    x10 = -x9 - C[2]
    x11 = x0 * (x10**2 + x6**2 + x8**2)
    x12 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x13 = x12 * boys(2, x11)
    x14 = x13 * x6
    x15 = -x14
    x16 = x12 * boys(1, x11)
    x17 = x16 * x5
    x18 = x15 + x17
    x19 = x18 * x5
    x20 = x4 * (-x13 + x16)
    x21 = x12 * boys(3, x11)
    x22 = x21 * x6
    x23 = -x22
    x24 = x13 * x5
    x25 = x23 + x24
    x26 = x25 * x6
    x27 = x20 - x26
    x28 = x19 + x27
    x29 = x28 * x5
    x30 = -x16 * x6
    x31 = x12 * boys(0, x11)
    x32 = x30 + x31 * x5
    x33 = x4 * (-x16 + x31)
    x34 = -x18 * x6 + x33
    x35 = x32 * x5 + x34
    x36 = x4 * (x14 - x17 + x32)
    x37 = -x28 * x6 + 2.0 * x36
    x38 = x35 * x5 + x37
    x39 = x25 * x5
    x40 = x4 * (x13 - x21)
    x41 = x12 * boys(4, x11)
    x42 = x41 * x6
    x43 = x21 * x5
    x44 = -x42 + x43
    x45 = x44 * x6
    x46 = x40 - x45
    x47 = x39 + x46
    x48 = x47 * x6
    x49 = x4 * (x18 + x22 - x24)
    x50 = 2.0 * x49
    x51 = x48 - x50
    x52 = -x48 + x50
    x53 = x29 + x52
    x54 = -x20
    x55 = x26 + x54
    x56 = x4 * (-x19 + x35 + x55)
    x57 = x3 * x38 - x53 * x6 + 3.0 * x56
    x58 = x3 * x53
    x59 = x47 * x5
    x60 = x44 * x5
    x61 = x4 * (x21 - x41)
    x62 = x12 * boys(5, x11)
    x63 = x6 * x62
    x64 = x41 * x5
    x65 = -x63 + x64
    x66 = x6 * x65
    x67 = x61 - x66
    x68 = x60 + x67
    x69 = x6 * x68
    x70 = x4 * (x25 + x42 - x43)
    x71 = 2.0 * x70
    x72 = -x69 + x71
    x73 = x59 + x72
    x74 = x6 * x73
    x75 = -x40
    x76 = x45 + x75
    x77 = x4 * (x28 - x39 + x76)
    x78 = 3.0 * x77
    x79 = x58 - x74 + x78
    x80 = x28 * x3
    x81 = x3 * x35 + x37
    x82 = x51 - x80 + x81
    x83 = 3.0 * x4
    x84 = x3 * x57 + x4 * (-x29 + x38 + x51) - x6 * x79 + x82 * x83
    x85 = x69 - x71
    x86 = x4 * (x53 - x59 + x85)
    x87 = x3 * x79
    x88 = x3 * x73
    x89 = x5 * x68
    x90 = x4 * (x41 - x62)
    x91 = x5 * x65
    x92 = x12 * boys(6, x11)
    x93 = x6 * x92
    x94 = x5 * x62
    x95 = -x93 + x94
    x96 = x6 * x95
    x97 = x90 + x91 - x96
    x98 = x6 * x97
    x99 = 2.0 * x4
    x100 = x99 * (x44 + x63 - x64)
    x101 = x100 - x98
    x102 = x101 + x89
    x103 = x102 * x6
    x104 = -x61
    x105 = x104 + x66
    x106 = x4 * (x105 + x47 - x60)
    x107 = 3.0 * x106
    x108 = -x103 + x107 + x88
    x109 = x108 * x6
    x110 = x3 * x47
    x111 = x52 + x80
    x112 = -x110 + x111 + x85
    x113 = x112 * x83
    x114 = -x109 + x113 + x86 + x87
    x115 = x110 + x72
    x116 = x115 * x6
    x117 = x111 * x3
    x118 = x25 * x3
    x119 = x18 * x3
    x120 = x119 + x27
    x121 = x99 * (-x118 + x120 + x76)
    x122 = x3 * x32 + x34
    x123 = -x111 * x6 + x3 * x81 + x56 + x99 * (-x119 + x122 + x55)
    x124 = -x100 + x98
    x125 = x4 * (x62 - x92)
    x126 = x12 * boys(7, x11)
    x127 = -x90
    x128 = x3 * x68
    x129 = x6 * (x101 + x128)
    x130 = x115 * x3
    x131 = x3 * x44
    x132 = x118 + x46
    x133 = x99 * (x105 - x131 + x132)
    x134 = -x116 + x117 + x121 + x77
    x135 = x13 * x3
    x136 = x16 * x3
    x137 = x4 * (-x135 + x136 + x15 + x22)
    x138 = x120 * x3
    x139 = x132 * x6
    x140 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**4.5)
    x141 = 9.12251705727742 * x140
    x142 = x13 * x8
    x143 = -x142
    x144 = -x7 - B[1]
    x145 = x144 * x16
    x146 = x143 + x145
    x147 = x146 * x5
    x148 = x21 * x8
    x149 = -x148
    x150 = x13 * x144
    x151 = x149 + x150
    x152 = x151 * x6
    x153 = -x152
    x154 = x147 + x153
    x155 = x154 * x5
    x156 = -x16 * x8
    x157 = x144 * x31 + x156
    x158 = -x146 * x6
    x159 = x157 * x5 + x158
    x160 = x4 * (x142 - x145 + x157)
    x161 = -x154 * x6 + x160
    x162 = x159 * x5 + x161
    x163 = x151 * x5
    x164 = x41 * x8
    x165 = -x164
    x166 = x144 * x21
    x167 = x165 + x166
    x168 = x167 * x6
    x169 = -x168
    x170 = x163 + x169
    x171 = x170 * x6
    x172 = x4 * (x146 + x148 - x150)
    x173 = -x172
    x174 = x171 + x173
    x175 = -x171 + x172
    x176 = x155 + x175
    x177 = x4 * (-x147 + x152 + x159)
    x178 = x162 * x3 - x176 * x6 + 2.0 * x177
    x179 = x176 * x3
    x180 = x170 * x5
    x181 = x4 * (x151 + x164 - x166)
    x182 = x167 * x5
    x183 = x62 * x8
    x184 = -x183
    x185 = x144 * x41
    x186 = x184 + x185
    x187 = x186 * x6
    x188 = -x187
    x189 = x182 + x188
    x190 = x189 * x6
    x191 = x181 - x190
    x192 = x180 + x191
    x193 = x192 * x6
    x194 = x4 * (x154 - x163 + x168)
    x195 = 2.0 * x194
    x196 = x179 - x193 + x195
    x197 = x154 * x3
    x198 = x159 * x3 + x161
    x199 = x99 * (x174 - x197 + x198)
    x200 = x178 * x3 - x196 * x6 + x199 + x4 * (-x155 + x162 + x174)
    x201 = -x181
    x202 = x190 + x201
    x203 = x4 * (x176 - x180 + x202)
    x204 = x196 * x3
    x205 = x192 * x3
    x206 = x189 * x5
    x207 = x4 * (x167 + x183 - x185)
    x208 = x186 * x5
    x209 = x8 * x92
    x210 = -x209
    x211 = x144 * x62
    x212 = x210 + x211
    x213 = x212 * x6
    x214 = x208 - x213
    x215 = x214 * x6
    x216 = x207 - x215
    x217 = x206 + x216
    x218 = x217 * x6
    x219 = x4 * (x170 - x182 + x187)
    x220 = 2.0 * x219
    x221 = x205 - x218 + x220
    x222 = x221 * x6
    x223 = x170 * x3
    x224 = x175 + x197
    x225 = x99 * (x202 - x223 + x224)
    x226 = x203 + x204 - x222 + x225
    x227 = x191 + x223
    x228 = x227 * x6
    x229 = x151 * x3
    x230 = x146 * x3
    x231 = x153 + x230
    x232 = x4 * (x168 - x229 + x231)
    x233 = x224 * x3
    x234 = x157 * x3 + x158
    x235 = x177 + x198 * x3 - x224 * x6 + x4 * (x152 - x230 + x234)
    x236 = -x207
    x237 = x215 + x236
    x238 = x4 * (x186 + x209 - x211)
    x239 = -x126 * x8
    x240 = x144 * x92 + x239
    x241 = x189 * x3
    x242 = x6 * (x216 + x241)
    x243 = x167 * x3
    x244 = x169 + x229
    x245 = x4 * (x187 - x243 + x244)
    x246 = x227 * x3
    x247 = x194 - x228 + x232 + x233
    x248 = x244 * x6
    x249 = x231 * x3
    x250 = -x225
    x251 = 20.3985682659737 * x140
    x252 = x10 * x13
    x253 = -x252
    x254 = -x9 - B[2]
    x255 = x16 * x254
    x256 = x253 + x255
    x257 = x256 * x5
    x258 = x10 * x21
    x259 = -x258
    x260 = x13 * x254
    x261 = x259 + x260
    x262 = x261 * x6
    x263 = -x262
    x264 = x257 + x263
    x265 = x264 * x5
    x266 = -x10 * x16
    x267 = x254 * x31 + x266
    x268 = -x256 * x6
    x269 = x267 * x5 + x268
    x270 = x4 * (x252 - x255 + x267)
    x271 = -x264 * x6 + x270
    x272 = x269 * x5 + x271
    x273 = x261 * x5
    x274 = x10 * x41
    x275 = -x274
    x276 = x21 * x254
    x277 = x275 + x276
    x278 = x277 * x6
    x279 = -x278
    x280 = x273 + x279
    x281 = x280 * x6
    x282 = x4 * (x256 + x258 - x260)
    x283 = -x282
    x284 = x281 + x283
    x285 = -x281 + x282
    x286 = x265 + x285
    x287 = x4 * (-x257 + x262 + x269)
    x288 = x272 * x3 - x286 * x6 + 2.0 * x287
    x289 = x286 * x3
    x290 = x280 * x5
    x291 = x4 * (x261 + x274 - x276)
    x292 = x277 * x5
    x293 = x10 * x62
    x294 = -x293
    x295 = x254 * x41
    x296 = x294 + x295
    x297 = x296 * x6
    x298 = -x297
    x299 = x292 + x298
    x300 = x299 * x6
    x301 = x291 - x300
    x302 = x290 + x301
    x303 = x302 * x6
    x304 = x4 * (x264 - x273 + x278)
    x305 = 2.0 * x304
    x306 = x289 - x303 + x305
    x307 = x264 * x3
    x308 = x269 * x3 + x271
    x309 = x99 * (x284 - x307 + x308)
    x310 = x288 * x3 - x306 * x6 + x309 + x4 * (-x265 + x272 + x284)
    x311 = -x291
    x312 = x300 + x311
    x313 = x4 * (x286 - x290 + x312)
    x314 = x3 * x306
    x315 = x3 * x302
    x316 = x299 * x5
    x317 = x4 * (x277 + x293 - x295)
    x318 = x296 * x5
    x319 = x10 * x92
    x320 = -x319
    x321 = x254 * x62
    x322 = x320 + x321
    x323 = x322 * x6
    x324 = x318 - x323
    x325 = x324 * x6
    x326 = x317 - x325
    x327 = x316 + x326
    x328 = x327 * x6
    x329 = x4 * (x280 - x292 + x297)
    x330 = 2.0 * x329
    x331 = x315 - x328 + x330
    x332 = x331 * x6
    x333 = x280 * x3
    x334 = x285 + x307
    x335 = x99 * (x312 - x333 + x334)
    x336 = x313 + x314 - x332 + x335
    x337 = x301 + x333
    x338 = x337 * x6
    x339 = x261 * x3
    x340 = x256 * x3
    x341 = x263 + x340
    x342 = x4 * (x278 - x339 + x341)
    x343 = x3 * x334
    x344 = x267 * x3 + x268
    x345 = x287 + x3 * x308 - x334 * x6 + x4 * (x262 - x340 + x344)
    x346 = -x317
    x347 = x325 + x346
    x348 = x4 * (x296 + x319 - x321)
    x349 = -x10 * x126
    x350 = x254 * x92 + x349
    x351 = x299 * x3
    x352 = x6 * (x326 + x351)
    x353 = x277 * x3
    x354 = x279 + x339
    x355 = x4 * (x297 - x353 + x354)
    x356 = x3 * x337
    x357 = x304 - x338 + x342 + x343
    x358 = x354 * x6
    x359 = x3 * x341
    x360 = -x335
    x361 = x144 * x146
    x362 = -x146 * x8 + x33
    x363 = x144 * x157 + x362
    x364 = x151 * x8
    x365 = x364 + x54
    x366 = x4 * (-x361 + x363 + x365)
    x367 = x20 - x364
    x368 = x361 + x367
    x369 = -x368 * x6
    x370 = x3 * x363 + x369
    x371 = x3 * x368
    x372 = x144 * x151
    x373 = x167 * x8
    x374 = -x373 + x40
    x375 = x372 + x374
    x376 = x375 * x6
    x377 = -x376
    x378 = x371 + x377
    x379 = x3 * x370 + x366 - x378 * x6
    x380 = x144 * x167
    x381 = x186 * x8
    x382 = x104 + x381
    x383 = x4 * (x375 - x380 + x382)
    x384 = x3 * x375
    x385 = -x381 + x61
    x386 = x380 + x385
    x387 = x386 * x6
    x388 = -x387
    x389 = x384 + x388
    x390 = x3 * x389
    x391 = x3 * x386
    x392 = x144 * x186
    x393 = x212 * x8
    x394 = -x393 + x90
    x395 = x392 + x394
    x396 = x395 * x6
    x397 = -x396
    x398 = x6 * (x391 + x397)
    x399 = x373 + x75
    x400 = x4 * (x368 - x372 + x399)
    x401 = x3 * x378
    x402 = x389 * x6
    x403 = x400 + x401 - x402
    x404 = x4 * (x378 - x384 + x387)
    x405 = x4 * (x370 - x371 + x376)
    x406 = -x400
    x407 = x368 * x5
    x408 = x363 * x5 + x369
    x409 = x377 + x407
    x410 = x3 * x408 + x366 - x409 * x6
    x411 = x3 * x409
    x412 = x375 * x5
    x413 = x388 + x412
    x414 = x413 * x6
    x415 = x400 + x411 - x414
    x416 = x3 * x410 + x4 * (x376 - x407 + x408) + x405 - x415 * x6
    x417 = x4 * (x387 + x409 - x412)
    x418 = x3 * x415
    x419 = x3 * x413
    x420 = x386 * x5
    x421 = x397 + x420
    x422 = x421 * x6
    x423 = x383 + x419 - x422
    x424 = x423 * x6
    x425 = x404 + x417 + x418 - x424
    x426 = -x383
    x427 = x127 + x393
    x428 = x4 * (x386 - x392 + x427)
    x429 = x125 - x240 * x8
    x430 = x144 * x212 + x429
    x431 = x261 * x8
    x432 = x144 * x256
    x433 = -x256 * x8
    x434 = x144 * x267 + x433
    x435 = x4 * (x431 - x432 + x434)
    x436 = -x431
    x437 = x432 + x436
    x438 = -x437 * x6
    x439 = x3 * x434 + x438
    x440 = x3 * x437
    x441 = x144 * x261
    x442 = x277 * x8
    x443 = -x442
    x444 = x441 + x443
    x445 = x444 * x6
    x446 = -x445
    x447 = x440 + x446
    x448 = x3 * x439 + x435 - x447 * x6
    x449 = x296 * x8
    x450 = x144 * x277
    x451 = x4 * (x444 + x449 - x450)
    x452 = x3 * x444
    x453 = -x449
    x454 = x450 + x453
    x455 = x454 * x6
    x456 = -x455
    x457 = x452 + x456
    x458 = x3 * x457
    x459 = x3 * x454
    x460 = x144 * x296
    x461 = x322 * x8
    x462 = -x461
    x463 = x460 + x462
    x464 = x463 * x6
    x465 = -x464
    x466 = x6 * (x459 + x465)
    x467 = x4 * (x437 - x441 + x442)
    x468 = x3 * x447
    x469 = x457 * x6
    x470 = x467 + x468 - x469
    x471 = x4 * (x447 - x452 + x455)
    x472 = x4 * (x439 - x440 + x445)
    x473 = -x467
    x474 = x437 * x5
    x475 = x434 * x5 + x438
    x476 = x446 + x474
    x477 = x3 * x475 + x435 - x476 * x6
    x478 = x3 * x476
    x479 = x444 * x5
    x480 = x456 + x479
    x481 = x480 * x6
    x482 = x467 + x478 - x481
    x483 = x3 * x477 + x4 * (x445 - x474 + x475) + x472 - x482 * x6
    x484 = x4 * (x455 + x476 - x479)
    x485 = x3 * x482
    x486 = x3 * x480
    x487 = x454 * x5
    x488 = x465 + x487
    x489 = x488 * x6
    x490 = x451 + x486 - x489
    x491 = x490 * x6
    x492 = x471 + x484 + x485 - x491
    x493 = -x451
    x494 = x4 * (x454 - x460 + x461)
    x495 = -x350 * x8
    x496 = x144 * x322 + x495
    x497 = 35.3313566383285 * x140
    x498 = x254 * x256
    x499 = -x10 * x256 + x33
    x500 = x254 * x267 + x499
    x501 = x10 * x261
    x502 = x501 + x54
    x503 = x4 * (-x498 + x500 + x502)
    x504 = x20 - x501
    x505 = x498 + x504
    x506 = -x505 * x6
    x507 = x3 * x500 + x506
    x508 = x3 * x505
    x509 = x254 * x261
    x510 = x10 * x277
    x511 = x40 - x510
    x512 = x509 + x511
    x513 = x512 * x6
    x514 = -x513
    x515 = x508 + x514
    x516 = x3 * x507 + x503 - x515 * x6
    x517 = x254 * x277
    x518 = x10 * x296
    x519 = x104 + x518
    x520 = x4 * (x512 - x517 + x519)
    x521 = x3 * x512
    x522 = -x518 + x61
    x523 = x517 + x522
    x524 = x523 * x6
    x525 = -x524
    x526 = x521 + x525
    x527 = x3 * x526
    x528 = x3 * x523
    x529 = x254 * x296
    x530 = x10 * x322
    x531 = -x530 + x90
    x532 = x529 + x531
    x533 = x532 * x6
    x534 = -x533
    x535 = x6 * (x528 + x534)
    x536 = x510 + x75
    x537 = x4 * (x505 - x509 + x536)
    x538 = x3 * x515
    x539 = x526 * x6
    x540 = x537 + x538 - x539
    x541 = x4 * (x515 - x521 + x524)
    x542 = x4 * (x507 - x508 + x513)
    x543 = -x537
    x544 = x5 * x505
    x545 = x5 * x500 + x506
    x546 = x514 + x544
    x547 = x3 * x545 + x503 - x546 * x6
    x548 = x3 * x546
    x549 = x5 * x512
    x550 = x525 + x549
    x551 = x550 * x6
    x552 = x537 + x548 - x551
    x553 = x3 * x547 + x4 * (x513 - x544 + x545) + x542 - x552 * x6
    x554 = x4 * (x524 + x546 - x549)
    x555 = x3 * x552
    x556 = x3 * x550
    x557 = x5 * x523
    x558 = x534 + x557
    x559 = x558 * x6
    x560 = x520 + x556 - x559
    x561 = x560 * x6
    x562 = x541 + x554 + x555 - x561
    x563 = -x520
    x564 = x127 + x530
    x565 = x4 * (x523 - x529 + x564)
    x566 = -x10 * x350 + x125
    x567 = x254 * x322 + x566
    x568 = x144 * x368
    x569 = 2.0 * x160 - x368 * x8
    x570 = x144 * x363 + x569
    x571 = x375 * x8
    x572 = 2.0 * x172
    x573 = x571 - x572
    x574 = x4 * (-x568 + x570 + x573)
    x575 = -x571 + x572
    x576 = x568 + x575
    x577 = x3 * x570 - x576 * x6
    x578 = x3 * x576
    x579 = x144 * x375
    x580 = x386 * x8
    x581 = 2.0 * x181
    x582 = -x580 + x581
    x583 = x579 + x582
    x584 = x583 * x6
    x585 = x578 - x584
    x586 = x3 * x577 + x574 - x585 * x6
    x587 = x580 - x581
    x588 = x4 * (x576 - x579 + x587)
    x589 = x3 * x585
    x590 = x3 * x583
    x591 = x144 * x386
    x592 = x395 * x8
    x593 = 2.0 * x207
    x594 = -x592 + x593
    x595 = x591 + x594
    x596 = x595 * x6
    x597 = x590 - x596
    x598 = x597 * x6
    x599 = x588 + x589 - x598
    x600 = x592 - x593
    x601 = x4 * (x583 - x591 + x600)
    x602 = 2.0 * x238 - x430 * x8
    x603 = x144 * x395 + x602
    x604 = -x588
    x605 = x144 * x437
    x606 = x270 - x437 * x8
    x607 = x144 * x434 + x606
    x608 = x444 * x8
    x609 = x283 + x608
    x610 = x4 * (-x605 + x607 + x609)
    x611 = x282 - x608
    x612 = x605 + x611
    x613 = x3 * x607 - x6 * x612
    x614 = x3 * x612
    x615 = x144 * x444
    x616 = x454 * x8
    x617 = x291 - x616
    x618 = x615 + x617
    x619 = x6 * x618
    x620 = x614 - x619
    x621 = x3 * x613 - x6 * x620 + x610
    x622 = x311 + x616
    x623 = x4 * (x612 - x615 + x622)
    x624 = x3 * x620
    x625 = x3 * x618
    x626 = x144 * x454
    x627 = x463 * x8
    x628 = x317 - x627
    x629 = x626 + x628
    x630 = x6 * x629
    x631 = x625 - x630
    x632 = x6 * x631
    x633 = x623 + x624 - x632
    x634 = x346 + x627
    x635 = x4 * (x618 - x626 + x634)
    x636 = x348 - x496 * x8
    x637 = x144 * x463 + x636
    x638 = -x623
    x639 = x512 * x8
    x640 = x144 * x505
    x641 = -x505 * x8
    x642 = x144 * x500 + x641
    x643 = x4 * (x639 - x640 + x642)
    x644 = -x639
    x645 = x640 + x644
    x646 = x3 * x642 - x6 * x645
    x647 = x3 * x645
    x648 = x144 * x512
    x649 = x523 * x8
    x650 = -x649
    x651 = x648 + x650
    x652 = x6 * x651
    x653 = x647 - x652
    x654 = x3 * x646 - x6 * x653 + x643
    x655 = x4 * (x645 - x648 + x649)
    x656 = x3 * x653
    x657 = x3 * x651
    x658 = x144 * x523
    x659 = x532 * x8
    x660 = -x659
    x661 = x658 + x660
    x662 = x6 * x661
    x663 = x657 - x662
    x664 = x6 * x663
    x665 = x655 + x656 - x664
    x666 = x4 * (x651 - x658 + x659)
    x667 = -x567 * x8
    x668 = x144 * x532 + x667
    x669 = -x655
    x670 = x254 * x505
    x671 = -x10 * x505 + 2.0 * x270
    x672 = x254 * x500 + x671
    x673 = x10 * x512
    x674 = 2.0 * x282
    x675 = x673 - x674
    x676 = x4 * (-x670 + x672 + x675)
    x677 = -x673 + x674
    x678 = x670 + x677
    x679 = x3 * x672 - x6 * x678
    x680 = x3 * x678
    x681 = x254 * x512
    x682 = x10 * x523
    x683 = 2.0 * x291
    x684 = -x682 + x683
    x685 = x681 + x684
    x686 = x6 * x685
    x687 = x680 - x686
    x688 = x3 * x679 - x6 * x687 + x676
    x689 = x682 - x683
    x690 = x4 * (x678 - x681 + x689)
    x691 = x3 * x687
    x692 = x3 * x685
    x693 = x254 * x523
    x694 = x10 * x532
    x695 = 2.0 * x317
    x696 = -x694 + x695
    x697 = x693 + x696
    x698 = x6 * x697
    x699 = x692 - x698
    x700 = x6 * x699
    x701 = x690 + x691 - x700
    x702 = x694 - x695
    x703 = x4 * (x685 - x693 + x702)
    x704 = -x10 * x567 + 2.0 * x348
    x705 = x254 * x532 + x704
    x706 = -x690
    x707 = -x7 - A[1]
    x708 = x16 * x707
    x709 = x143 + x708
    x710 = x5 * x709
    x711 = x13 * x707
    x712 = x149 + x711
    x713 = x6 * x712
    x714 = x710 - x713
    x715 = x5 * x714
    x716 = x4 * (x148 + x709 - x711)
    x717 = x5 * x712
    x718 = x21 * x707
    x719 = x165 + x718
    x720 = x6 * x719
    x721 = x717 - x720
    x722 = x6 * x721
    x723 = x716 - x722
    x724 = x715 + x723
    x725 = x5 * x724
    x726 = x156 + x31 * x707
    x727 = x5 * x726 - x6 * x709
    x728 = x4 * (x142 - x708 + x726)
    x729 = -x6 * x714 + x728
    x730 = x5 * x727 + x729
    x731 = -x6 * x724 + x99 * (-x710 + x713 + x727)
    x732 = x5 * x730 + x731
    x733 = x4 * (x164 + x712 - x718)
    x734 = x5 * x721
    x735 = x5 * x719
    x736 = x41 * x707
    x737 = x184 + x736
    x738 = x6 * x737
    x739 = x735 - x738
    x740 = x6 * x739
    x741 = x733 + x734 - x740
    x742 = x6 * x741
    x743 = x99 * (x714 - x717 + x720)
    x744 = x742 - x743
    x745 = -x742 + x743
    x746 = x725 + x745
    x747 = -x716
    x748 = x722 + x747
    x749 = x4 * (-x715 + x730 + x748)
    x750 = x3 * x732 - x6 * x746 + 3.0 * x749
    x751 = x3 * x746
    x752 = x5 * x741
    x753 = x4 * (x183 + x719 - x736)
    x754 = x5 * x739
    x755 = x5 * x737
    x756 = x62 * x707
    x757 = x210 + x756
    x758 = x6 * x757
    x759 = x755 - x758
    x760 = x6 * x759
    x761 = x753 + x754 - x760
    x762 = x6 * x761
    x763 = x99 * (x721 - x735 + x738)
    x764 = -x762 + x763
    x765 = x752 + x764
    x766 = x6 * x765
    x767 = -x733
    x768 = x740 + x767
    x769 = x4 * (x724 - x734 + x768)
    x770 = 3.0 * x769
    x771 = x751 - x766 + x770
    x772 = x3 * x724
    x773 = x3 * x730 + x731
    x774 = x762 - x763
    x775 = x4 * (x209 + x737 - x756)
    x776 = x239 + x707 * x92
    x777 = -x753
    x778 = x3 * x741
    x779 = x745 + x772
    x780 = x3 * x714
    x781 = 24.1359114645008 * x140
    x782 = x146 * x707
    x783 = x367 + x782
    x784 = x5 * x783
    x785 = x151 * x707
    x786 = x374 + x785
    x787 = x6 * x786
    x788 = -x787
    x789 = x784 + x788
    x790 = x5 * x789
    x791 = x157 * x707 + x362
    x792 = -x6 * x783
    x793 = x5 * x791 + x792
    x794 = x4 * (x365 - x782 + x791)
    x795 = -x6 * x789 + x794
    x796 = x5 * x793 + x795
    x797 = x5 * x786
    x798 = x167 * x707
    x799 = x385 + x798
    x800 = x6 * x799
    x801 = x797 - x800
    x802 = x6 * x801
    x803 = x4 * (x399 + x783 - x785)
    x804 = x802 - x803
    x805 = -x802 + x803
    x806 = x790 + x805
    x807 = x4 * (-x784 + x787 + x793)
    x808 = x3 * x796 - x6 * x806 + 2.0 * x807
    x809 = x3 * x806
    x810 = x5 * x801
    x811 = x4 * (x382 + x786 - x798)
    x812 = x5 * x799
    x813 = x186 * x707
    x814 = x394 + x813
    x815 = x6 * x814
    x816 = x812 - x815
    x817 = x6 * x816
    x818 = x811 - x817
    x819 = x810 + x818
    x820 = x6 * x819
    x821 = x4 * (x789 - x797 + x800)
    x822 = 2.0 * x821
    x823 = x809 - x820 + x822
    x824 = x3 * x789
    x825 = x3 * x793 + x795
    x826 = -x811 + x817
    x827 = x4 * (x427 + x799 - x813)
    x828 = x212 * x707 + x429
    x829 = x3 * x801
    x830 = x805 + x824
    x831 = x3 * x783
    x832 = 53.9695387335403 * x140
    x833 = x256 * x707
    x834 = x436 + x833
    x835 = x5 * x834
    x836 = x261 * x707
    x837 = x443 + x836
    x838 = x6 * x837
    x839 = -x838
    x840 = x835 + x839
    x841 = x5 * x840
    x842 = x267 * x707 + x433
    x843 = -x6 * x834
    x844 = x5 * x842 + x843
    x845 = x4 * (x431 - x833 + x842)
    x846 = -x6 * x840 + x845
    x847 = x5 * x844 + x846
    x848 = x5 * x837
    x849 = x277 * x707
    x850 = x453 + x849
    x851 = x6 * x850
    x852 = x848 - x851
    x853 = x6 * x852
    x854 = x4 * (x442 + x834 - x836)
    x855 = -x854
    x856 = x853 + x855
    x857 = -x853 + x854
    x858 = x841 + x857
    x859 = x4 * (-x835 + x838 + x844)
    x860 = x3 * x847 - x6 * x858 + 2.0 * x859
    x861 = x3 * x858
    x862 = x5 * x852
    x863 = x4 * (x449 + x837 - x849)
    x864 = x5 * x850
    x865 = x296 * x707
    x866 = x462 + x865
    x867 = x6 * x866
    x868 = x864 - x867
    x869 = x6 * x868
    x870 = x863 - x869
    x871 = x862 + x870
    x872 = x6 * x871
    x873 = x4 * (x840 - x848 + x851)
    x874 = 2.0 * x873
    x875 = x861 - x872 + x874
    x876 = x3 * x840
    x877 = x3 * x844 + x846
    x878 = -x863
    x879 = x869 + x878
    x880 = x4 * (x461 + x850 - x865)
    x881 = x322 * x707 + x495
    x882 = x3 * x852
    x883 = x857 + x876
    x884 = x3 * x834
    x885 = x368 * x707
    x886 = x363 * x707 + x569
    x887 = x4 * (x573 - x885 + x886)
    x888 = x575 + x885
    x889 = -x6 * x888
    x890 = x3 * x886 + x889
    x891 = x375 * x707
    x892 = x582 + x891
    x893 = x3 * x892
    x894 = x386 * x707
    x895 = x594 + x894
    x896 = x6 * x895
    x897 = -x896
    x898 = x4 * (x587 + x888 - x891)
    x899 = -x898
    x900 = x3 * x888
    x901 = x6 * x892
    x902 = -x901
    x903 = x900 + x902
    x904 = x5 * x888
    x905 = x5 * x886 + x889
    x906 = x902 + x904
    x907 = x3 * x905 - x6 * x906 + x887
    x908 = x3 * x906
    x909 = x5 * x892
    x910 = x897 + x909
    x911 = x6 * x910
    x912 = x898 + x908 - x911
    x913 = x4 * (x600 + x892 - x894)
    x914 = x395 * x707 + x602
    x915 = x437 * x707
    x916 = x434 * x707 + x606
    x917 = x4 * (x609 - x915 + x916)
    x918 = x611 + x915
    x919 = -x6 * x918
    x920 = x3 * x916 + x919
    x921 = x444 * x707
    x922 = x617 + x921
    x923 = x3 * x922
    x924 = x454 * x707
    x925 = x628 + x924
    x926 = x6 * x925
    x927 = -x926
    x928 = x4 * (x622 + x918 - x921)
    x929 = -x928
    x930 = x3 * x918
    x931 = x6 * x922
    x932 = -x931
    x933 = x930 + x932
    x934 = x5 * x918
    x935 = x5 * x916 + x919
    x936 = x932 + x934
    x937 = x3 * x935 - x6 * x936 + x917
    x938 = x3 * x936
    x939 = x5 * x922
    x940 = x927 + x939
    x941 = x6 * x940
    x942 = x928 + x938 - x941
    x943 = x4 * (x634 + x922 - x924)
    x944 = x463 * x707 + x636
    x945 = 93.4779831475484 * x140
    x946 = x505 * x707
    x947 = x500 * x707 + x641
    x948 = x4 * (x639 - x946 + x947)
    x949 = x644 + x946
    x950 = -x6 * x949
    x951 = x3 * x947 + x950
    x952 = x512 * x707
    x953 = x650 + x952
    x954 = x3 * x953
    x955 = x523 * x707
    x956 = x660 + x955
    x957 = x6 * x956
    x958 = -x957
    x959 = x4 * (x649 + x949 - x952)
    x960 = -x959
    x961 = x3 * x949
    x962 = x6 * x953
    x963 = -x962
    x964 = x961 + x963
    x965 = x5 * x949
    x966 = x5 * x947 + x950
    x967 = x963 + x965
    x968 = x3 * x966 - x6 * x967 + x948
    x969 = x3 * x967
    x970 = x5 * x953
    x971 = x958 + x970
    x972 = x6 * x971
    x973 = x959 + x969 - x972
    x974 = x4 * (x659 + x953 - x955)
    x975 = x532 * x707 + x667
    x976 = x583 * x8
    x977 = x576 * x707
    x978 = 3.0 * x400
    x979 = 3.0 * x366 + x570 * x707 - x576 * x8
    x980 = x4 * (x976 - x977 - x978 + x979)
    x981 = -x976 + x977 + x978
    x982 = x3 * x979 - x6 * x981
    x983 = x3 * x981
    x984 = x583 * x707
    x985 = x595 * x8
    x986 = 3.0 * x383
    x987 = x984 - x985 + x986
    x988 = x6 * x987
    x989 = x983 - x988
    x990 = x4 * (x981 - x984 + x985 - x986)
    x991 = 3.0 * x428 + x595 * x707 - x603 * x8
    x992 = x618 * x8
    x993 = x612 * x707
    x994 = 2.0 * x467
    x995 = 2.0 * x435 + x607 * x707 - x612 * x8
    x996 = x4 * (x992 - x993 - x994 + x995)
    x997 = -x992 + x993 + x994
    x998 = x3 * x995 - x6 * x997
    x999 = x3 * x997
    x1000 = x618 * x707
    x1001 = x629 * x8
    x1002 = 2.0 * x451
    x1003 = x1000 - x1001 + x1002
    x1004 = x1003 * x6
    x1005 = -x1004 + x999
    x1006 = x4 * (-x1000 + x1001 - x1002 + x997)
    x1007 = 2.0 * x494 + x629 * x707 - x637 * x8
    x1008 = x651 * x8
    x1009 = x645 * x707
    x1010 = x503 + x642 * x707 - x645 * x8
    x1011 = x4 * (x1008 - x1009 + x1010 + x543)
    x1012 = -x1008 + x1009 + x537
    x1013 = x1010 * x3 - x1012 * x6
    x1014 = x1012 * x3
    x1015 = x651 * x707
    x1016 = x661 * x8
    x1017 = x1015 - x1016 + x520
    x1018 = x1017 * x6
    x1019 = x1014 - x1018
    x1020 = x4 * (x1012 - x1015 + x1016 + x563)
    x1021 = x565 + x661 * x707 - x668 * x8
    x1022 = x685 * x8
    x1023 = x678 * x707
    x1024 = x672 * x707 - x678 * x8
    x1025 = x4 * (x1022 - x1023 + x1024)
    x1026 = -x1022 + x1023
    x1027 = x1024 * x3 - x1026 * x6
    x1028 = x1026 * x3
    x1029 = x685 * x707
    x1030 = x697 * x8
    x1031 = x1029 - x1030
    x1032 = x1031 * x6
    x1033 = x1028 - x1032
    x1034 = x4 * (x1026 - x1029 + x1030)
    x1035 = x697 * x707 - x705 * x8
    x1036 = -x9 - A[2]
    x1037 = x1036 * x16
    x1038 = x1037 + x253
    x1039 = x1038 * x5
    x1040 = x1036 * x13
    x1041 = x1040 + x259
    x1042 = x1041 * x6
    x1043 = x1039 - x1042
    x1044 = x1043 * x5
    x1045 = x4 * (x1038 - x1040 + x258)
    x1046 = x1041 * x5
    x1047 = x1036 * x21
    x1048 = x1047 + x275
    x1049 = x1048 * x6
    x1050 = x1046 - x1049
    x1051 = x1050 * x6
    x1052 = x1045 - x1051
    x1053 = x1044 + x1052
    x1054 = x1053 * x5
    x1055 = x1036 * x31 + x266
    x1056 = -x1038 * x6 + x1055 * x5
    x1057 = x4 * (-x1037 + x1055 + x252)
    x1058 = -x1043 * x6 + x1057
    x1059 = x1056 * x5 + x1058
    x1060 = -x1053 * x6 + x99 * (-x1039 + x1042 + x1056)
    x1061 = x1059 * x5 + x1060
    x1062 = x4 * (x1041 - x1047 + x274)
    x1063 = x1050 * x5
    x1064 = x1048 * x5
    x1065 = x1036 * x41
    x1066 = x1065 + x294
    x1067 = x1066 * x6
    x1068 = x1064 - x1067
    x1069 = x1068 * x6
    x1070 = x1062 + x1063 - x1069
    x1071 = x1070 * x6
    x1072 = x99 * (x1043 - x1046 + x1049)
    x1073 = x1071 - x1072
    x1074 = -x1071 + x1072
    x1075 = x1054 + x1074
    x1076 = -x1045
    x1077 = x1051 + x1076
    x1078 = x4 * (-x1044 + x1059 + x1077)
    x1079 = x1061 * x3 - x1075 * x6 + 3.0 * x1078
    x1080 = x1075 * x3
    x1081 = x1070 * x5
    x1082 = x4 * (x1048 - x1065 + x293)
    x1083 = x1068 * x5
    x1084 = x1066 * x5
    x1085 = x1036 * x62
    x1086 = x1085 + x320
    x1087 = x1086 * x6
    x1088 = x1084 - x1087
    x1089 = x1088 * x6
    x1090 = x1082 + x1083 - x1089
    x1091 = x1090 * x6
    x1092 = x99 * (x1050 - x1064 + x1067)
    x1093 = -x1091 + x1092
    x1094 = x1081 + x1093
    x1095 = x1094 * x6
    x1096 = -x1062
    x1097 = x1069 + x1096
    x1098 = x4 * (x1053 - x1063 + x1097)
    x1099 = 3.0 * x1098
    x1100 = x1080 - x1095 + x1099
    x1101 = x1053 * x3
    x1102 = x1059 * x3 + x1060
    x1103 = x1091 - x1092
    x1104 = x4 * (x1066 - x1085 + x319)
    x1105 = x1036 * x92 + x349
    x1106 = -x1082
    x1107 = x1070 * x3
    x1108 = x1074 + x1101
    x1109 = x1043 * x3
    x1110 = x1038 * x144
    x1111 = x1041 * x8
    x1112 = -x1111
    x1113 = x1110 + x1112
    x1114 = x1113 * x5
    x1115 = x1041 * x144
    x1116 = x1048 * x8
    x1117 = -x1116
    x1118 = x1115 + x1117
    x1119 = x1118 * x6
    x1120 = -x1119
    x1121 = x1114 + x1120
    x1122 = x1121 * x5
    x1123 = -x1038 * x8
    x1124 = x1055 * x144 + x1123
    x1125 = -x1113 * x6
    x1126 = x1124 * x5 + x1125
    x1127 = x4 * (-x1110 + x1111 + x1124)
    x1128 = -x1121 * x6 + x1127
    x1129 = x1126 * x5 + x1128
    x1130 = x1118 * x5
    x1131 = x1048 * x144
    x1132 = x1066 * x8
    x1133 = -x1132
    x1134 = x1131 + x1133
    x1135 = x1134 * x6
    x1136 = x1130 - x1135
    x1137 = x1136 * x6
    x1138 = x4 * (x1113 - x1115 + x1116)
    x1139 = -x1138
    x1140 = x1137 + x1139
    x1141 = -x1137 + x1138
    x1142 = x1122 + x1141
    x1143 = x4 * (-x1114 + x1119 + x1126)
    x1144 = x1129 * x3 - x1142 * x6 + 2.0 * x1143
    x1145 = x1142 * x3
    x1146 = x1136 * x5
    x1147 = x4 * (x1118 - x1131 + x1132)
    x1148 = x1134 * x5
    x1149 = x1066 * x144
    x1150 = x1086 * x8
    x1151 = -x1150
    x1152 = x1149 + x1151
    x1153 = x1152 * x6
    x1154 = x1148 - x1153
    x1155 = x1154 * x6
    x1156 = x1147 - x1155
    x1157 = x1146 + x1156
    x1158 = x1157 * x6
    x1159 = x4 * (x1121 - x1130 + x1135)
    x1160 = 2.0 * x1159
    x1161 = x1145 - x1158 + x1160
    x1162 = x1121 * x3
    x1163 = x1126 * x3 + x1128
    x1164 = -x1147
    x1165 = x1155 + x1164
    x1166 = x4 * (x1134 - x1149 + x1150)
    x1167 = -x1105 * x8
    x1168 = x1086 * x144 + x1167
    x1169 = x1136 * x3
    x1170 = x1141 + x1162
    x1171 = x1113 * x3
    x1172 = x1036 * x256
    x1173 = x1172 + x504
    x1174 = x1173 * x5
    x1175 = x1036 * x261
    x1176 = x1175 + x511
    x1177 = x1176 * x6
    x1178 = -x1177
    x1179 = x1174 + x1178
    x1180 = x1179 * x5
    x1181 = x1036 * x267 + x499
    x1182 = -x1173 * x6
    x1183 = x1181 * x5 + x1182
    x1184 = x4 * (-x1172 + x1181 + x502)
    x1185 = -x1179 * x6 + x1184
    x1186 = x1183 * x5 + x1185
    x1187 = x1176 * x5
    x1188 = x1036 * x277
    x1189 = x1188 + x522
    x1190 = x1189 * x6
    x1191 = x1187 - x1190
    x1192 = x1191 * x6
    x1193 = x4 * (x1173 - x1175 + x536)
    x1194 = -x1193
    x1195 = x1192 + x1194
    x1196 = -x1192 + x1193
    x1197 = x1180 + x1196
    x1198 = x4 * (-x1174 + x1177 + x1183)
    x1199 = x1186 * x3 - x1197 * x6 + 2.0 * x1198
    x1200 = x1197 * x3
    x1201 = x1191 * x5
    x1202 = x4 * (x1176 - x1188 + x519)
    x1203 = x1189 * x5
    x1204 = x1036 * x296
    x1205 = x1204 + x531
    x1206 = x1205 * x6
    x1207 = x1203 - x1206
    x1208 = x1207 * x6
    x1209 = x1202 - x1208
    x1210 = x1201 + x1209
    x1211 = x1210 * x6
    x1212 = x4 * (x1179 - x1187 + x1190)
    x1213 = 2.0 * x1212
    x1214 = x1200 - x1211 + x1213
    x1215 = x1179 * x3
    x1216 = x1183 * x3 + x1185
    x1217 = -x1202
    x1218 = x1208 + x1217
    x1219 = x4 * (x1189 - x1204 + x564)
    x1220 = x1036 * x322 + x566
    x1221 = x1191 * x3
    x1222 = x1196 + x1215
    x1223 = x1173 * x3
    x1224 = x1113 * x144
    x1225 = x1057 - x1113 * x8
    x1226 = x1124 * x144 + x1225
    x1227 = x1118 * x8
    x1228 = x1076 + x1227
    x1229 = x4 * (-x1224 + x1226 + x1228)
    x1230 = x1045 - x1227
    x1231 = x1224 + x1230
    x1232 = -x1231 * x6
    x1233 = x1226 * x3 + x1232
    x1234 = x1118 * x144
    x1235 = x1134 * x8
    x1236 = x1062 - x1235
    x1237 = x1234 + x1236
    x1238 = x1237 * x3
    x1239 = x1134 * x144
    x1240 = x1152 * x8
    x1241 = x1082 - x1240
    x1242 = x1239 + x1241
    x1243 = x1242 * x6
    x1244 = -x1243
    x1245 = x1096 + x1235
    x1246 = x4 * (x1231 - x1234 + x1245)
    x1247 = -x1246
    x1248 = x1231 * x3
    x1249 = x1237 * x6
    x1250 = -x1249
    x1251 = x1248 + x1250
    x1252 = x1231 * x5
    x1253 = x1226 * x5 + x1232
    x1254 = x1250 + x1252
    x1255 = x1229 + x1253 * x3 - x1254 * x6
    x1256 = x1254 * x3
    x1257 = x1237 * x5
    x1258 = x1244 + x1257
    x1259 = x1258 * x6
    x1260 = x1246 + x1256 - x1259
    x1261 = x1106 + x1240
    x1262 = x4 * (x1237 - x1239 + x1261)
    x1263 = x1104 - x1168 * x8
    x1264 = x1152 * x144 + x1263
    x1265 = x1176 * x8
    x1266 = x1173 * x144
    x1267 = -x1173 * x8
    x1268 = x1181 * x144 + x1267
    x1269 = x4 * (x1265 - x1266 + x1268)
    x1270 = -x1265
    x1271 = x1266 + x1270
    x1272 = -x1271 * x6
    x1273 = x1268 * x3 + x1272
    x1274 = x1176 * x144
    x1275 = x1189 * x8
    x1276 = -x1275
    x1277 = x1274 + x1276
    x1278 = x1277 * x3
    x1279 = x1189 * x144
    x1280 = x1205 * x8
    x1281 = -x1280
    x1282 = x1279 + x1281
    x1283 = x1282 * x6
    x1284 = -x1283
    x1285 = x4 * (x1271 - x1274 + x1275)
    x1286 = -x1285
    x1287 = x1271 * x3
    x1288 = x1277 * x6
    x1289 = -x1288
    x1290 = x1287 + x1289
    x1291 = x1271 * x5
    x1292 = x1268 * x5 + x1272
    x1293 = x1289 + x1291
    x1294 = x1269 + x1292 * x3 - x1293 * x6
    x1295 = x1293 * x3
    x1296 = x1277 * x5
    x1297 = x1284 + x1296
    x1298 = x1297 * x6
    x1299 = x1285 + x1295 - x1298
    x1300 = x4 * (x1277 - x1279 + x1280)
    x1301 = -x1220 * x8
    x1302 = x1205 * x144 + x1301
    x1303 = x1036 * x505
    x1304 = x1036 * x500 + x671
    x1305 = x4 * (-x1303 + x1304 + x675)
    x1306 = x1303 + x677
    x1307 = -x1306 * x6
    x1308 = x1304 * x3 + x1307
    x1309 = x1036 * x512
    x1310 = x1309 + x684
    x1311 = x1310 * x3
    x1312 = x1036 * x523
    x1313 = x1312 + x696
    x1314 = x1313 * x6
    x1315 = -x1314
    x1316 = x4 * (x1306 - x1309 + x689)
    x1317 = -x1316
    x1318 = x1306 * x3
    x1319 = x1310 * x6
    x1320 = -x1319
    x1321 = x1318 + x1320
    x1322 = x1306 * x5
    x1323 = x1304 * x5 + x1307
    x1324 = x1320 + x1322
    x1325 = x1305 + x1323 * x3 - x1324 * x6
    x1326 = x1324 * x3
    x1327 = x1310 * x5
    x1328 = x1315 + x1327
    x1329 = x1328 * x6
    x1330 = x1316 + x1326 - x1329
    x1331 = x4 * (x1310 - x1312 + x702)
    x1332 = x1036 * x532 + x704
    x1333 = x1231 * x144
    x1334 = 2.0 * x1127 - x1231 * x8
    x1335 = x1226 * x144 + x1334
    x1336 = x1237 * x8
    x1337 = 2.0 * x1138
    x1338 = x1336 - x1337
    x1339 = x4 * (-x1333 + x1335 + x1338)
    x1340 = -x1336 + x1337
    x1341 = x1333 + x1340
    x1342 = x1335 * x3 - x1341 * x6
    x1343 = x1341 * x3
    x1344 = x1237 * x144
    x1345 = x1242 * x8
    x1346 = 2.0 * x1147
    x1347 = -x1345 + x1346
    x1348 = x1344 + x1347
    x1349 = x1348 * x6
    x1350 = x1343 - x1349
    x1351 = x1345 - x1346
    x1352 = x4 * (x1341 - x1344 + x1351)
    x1353 = 2.0 * x1166 - x1264 * x8
    x1354 = x1242 * x144 + x1353
    x1355 = x1271 * x144
    x1356 = x1184 - x1271 * x8
    x1357 = x1268 * x144 + x1356
    x1358 = x1277 * x8
    x1359 = x1194 + x1358
    x1360 = x4 * (-x1355 + x1357 + x1359)
    x1361 = x1193 - x1358
    x1362 = x1355 + x1361
    x1363 = x1357 * x3 - x1362 * x6
    x1364 = x1362 * x3
    x1365 = x1277 * x144
    x1366 = x1282 * x8
    x1367 = x1202 - x1366
    x1368 = x1365 + x1367
    x1369 = x1368 * x6
    x1370 = x1364 - x1369
    x1371 = x1217 + x1366
    x1372 = x4 * (x1362 - x1365 + x1371)
    x1373 = x1219 - x1302 * x8
    x1374 = x1282 * x144 + x1373
    x1375 = x1310 * x8
    x1376 = x1306 * x144
    x1377 = -x1306 * x8
    x1378 = x1304 * x144 + x1377
    x1379 = x4 * (x1375 - x1376 + x1378)
    x1380 = -x1375
    x1381 = x1376 + x1380
    x1382 = x1378 * x3 - x1381 * x6
    x1383 = x1381 * x3
    x1384 = x1310 * x144
    x1385 = x1313 * x8
    x1386 = -x1385
    x1387 = x1384 + x1386
    x1388 = x1387 * x6
    x1389 = x1383 - x1388
    x1390 = x4 * (x1381 - x1384 + x1385)
    x1391 = -x1332 * x8
    x1392 = x1313 * x144 + x1391
    x1393 = x10 * x685
    x1394 = x1036 * x678
    x1395 = 3.0 * x537
    x1396 = -x10 * x678 + x1036 * x672 + 3.0 * x503
    x1397 = x4 * (x1393 - x1394 - x1395 + x1396)
    x1398 = -x1393 + x1394 + x1395
    x1399 = x1396 * x3 - x1398 * x6
    x1400 = x1398 * x3
    x1401 = x1036 * x685
    x1402 = x10 * x697
    x1403 = 3.0 * x520
    x1404 = x1401 - x1402 + x1403
    x1405 = x1404 * x6
    x1406 = x1400 - x1405
    x1407 = x4 * (x1398 - x1401 + x1402 - x1403)
    x1408 = -x10 * x705 + x1036 * x697 + 3.0 * x565
    x1409 = x719 * x8
    x1410 = x707 * x712
    x1411 = x707 * x709
    x1412 = x712 * x8
    x1413 = x1411 - x1412 + x20
    x1414 = x4 * (x1409 - x1410 + x1413 + x75)
    x1415 = x1413 * x5
    x1416 = -x1409 + x1410 + x40
    x1417 = x1416 * x6
    x1418 = x1415 - x1417
    x1419 = x1418 * x5
    x1420 = x1416 * x5
    x1421 = x707 * x719
    x1422 = x737 * x8
    x1423 = x1421 - x1422 + x61
    x1424 = x1423 * x6
    x1425 = x1420 - x1424
    x1426 = x1425 * x6
    x1427 = x1414 + x1419 - x1426
    x1428 = x1427 * x5
    x1429 = x33 + x707 * x726 - x709 * x8
    x1430 = x4 * (-x1411 + x1412 + x1429 + x54)
    x1431 = -x1413 * x6 + x1429 * x5
    x1432 = -x1418 * x6 + x1430 + x1431 * x5
    x1433 = -x1427 * x6 + x99 * (-x1415 + x1417 + x1431)
    x1434 = x1432 * x5 + x1433
    x1435 = x4 * (x104 + x1416 - x1421 + x1422)
    x1436 = x1425 * x5
    x1437 = x1423 * x5
    x1438 = x707 * x737
    x1439 = x757 * x8
    x1440 = x1438 - x1439 + x90
    x1441 = x1440 * x6
    x1442 = x1437 - x1441
    x1443 = x1442 * x6
    x1444 = x1435 + x1436 - x1443
    x1445 = x1444 * x6
    x1446 = x99 * (x1418 - x1420 + x1424)
    x1447 = x1445 - x1446
    x1448 = x1428 - x1445 + x1446
    x1449 = -x1414
    x1450 = x4 * (x127 + x1423 - x1438 + x1439)
    x1451 = x125 + x707 * x757 - x776 * x8
    x1452 = -x1435
    x1453 = 31.1593277158494 * x140
    x1454 = x707 * x783
    x1455 = x786 * x8
    x1456 = x1454 - x1455 + x172 + x716
    x1457 = x1456 * x5
    x1458 = x707 * x786
    x1459 = x799 * x8
    x1460 = x1458 - x1459 + x181 + x733
    x1461 = x1460 * x6
    x1462 = x1457 - x1461
    x1463 = x1462 * x5
    x1464 = x160 + x707 * x791 + x728 - x783 * x8
    x1465 = -x1456 * x6 + x1464 * x5
    x1466 = x4 * (-x1454 + x1455 + x1464 + x173 + x747)
    x1467 = -x1462 * x6 + x1466
    x1468 = x1465 * x5 + x1467
    x1469 = x1460 * x5
    x1470 = x707 * x799
    x1471 = x8 * x814
    x1472 = x1470 - x1471 + x207 + x753
    x1473 = x1472 * x6
    x1474 = x1469 - x1473
    x1475 = x1474 * x6
    x1476 = x4 * (x1456 - x1458 + x1459 + x201 + x767)
    x1477 = x1475 - x1476
    x1478 = x1463 - x1475 + x1476
    x1479 = x4 * (x1460 - x1470 + x1471 + x236 + x777)
    x1480 = x238 + x707 * x814 + x775 - x8 * x828
    x1481 = 69.6743749058326 * x140
    x1482 = x707 * x834
    x1483 = x8 * x837
    x1484 = x1482 - x1483 + x282
    x1485 = x1484 * x5
    x1486 = x707 * x837
    x1487 = x8 * x850
    x1488 = x1486 - x1487 + x291
    x1489 = x1488 * x6
    x1490 = x1485 - x1489
    x1491 = x1490 * x5
    x1492 = x270 + x707 * x842 - x8 * x834
    x1493 = -x1484 * x6 + x1492 * x5
    x1494 = x4 * (-x1482 + x1483 + x1492 + x283)
    x1495 = -x1490 * x6 + x1494
    x1496 = x1493 * x5 + x1495
    x1497 = x1488 * x5
    x1498 = x707 * x850
    x1499 = x8 * x866
    x1500 = x1498 - x1499 + x317
    x1501 = x1500 * x6
    x1502 = x1497 - x1501
    x1503 = x1502 * x6
    x1504 = x4 * (x1484 - x1486 + x1487 + x311)
    x1505 = -x1504
    x1506 = x1503 + x1505
    x1507 = x1491 - x1503 + x1504
    x1508 = x4 * (x1488 - x1498 + x1499 + x346)
    x1509 = x348 + x707 * x866 - x8 * x881
    x1510 = 2.0 * x794
    x1511 = x1510 + x366 + x707 * x886 - x8 * x888
    x1512 = x707 * x892
    x1513 = x8 * x895
    x1514 = 2.0 * x811
    x1515 = x1512 - x1513 + x1514 + x383
    x1516 = x1515 * x6
    x1517 = x707 * x888
    x1518 = x8 * x892
    x1519 = 2.0 * x803
    x1520 = x1517 - x1518 + x1519 + x400
    x1521 = -x1520 * x6
    x1522 = x1520 * x5
    x1523 = x1511 * x5 + x1521
    x1524 = -x1519
    x1525 = x4 * (x1511 - x1517 + x1518 + x1524 + x406)
    x1526 = -x1516 + x1522
    x1527 = -x1514
    x1528 = x4 * (-x1512 + x1513 + x1520 + x1527 + x426)
    x1529 = 2.0 * x827
    x1530 = x1529 + x428 + x707 * x895 - x8 * x914
    x1531 = x435 + x707 * x916 - x8 * x918 + x845
    x1532 = x707 * x922
    x1533 = x8 * x925
    x1534 = x1532 - x1533 + x451 + x863
    x1535 = x1534 * x6
    x1536 = x707 * x918
    x1537 = x8 * x922
    x1538 = x1536 - x1537 + x467 + x854
    x1539 = -x1538 * x6
    x1540 = x1538 * x5
    x1541 = x1531 * x5 + x1539
    x1542 = x4 * (x1531 - x1536 + x1537 + x473 + x855)
    x1543 = -x1535 + x1540
    x1544 = x4 * (-x1532 + x1533 + x1538 + x493 + x878)
    x1545 = x494 + x707 * x925 - x8 * x944 + x880
    x1546 = 120.679557322504 * x140
    x1547 = x503 + x707 * x947 - x8 * x949
    x1548 = x707 * x953
    x1549 = x8 * x956
    x1550 = x1548 - x1549 + x520
    x1551 = x1550 * x6
    x1552 = x707 * x949
    x1553 = x8 * x953
    x1554 = x1552 - x1553 + x537
    x1555 = -x1554 * x6
    x1556 = x1554 * x5
    x1557 = x1547 * x5 + x1555
    x1558 = x4 * (x1547 - x1552 + x1553 + x543)
    x1559 = -x1551 + x1556
    x1560 = x4 * (-x1548 + x1549 + x1554 + x563)
    x1561 = x565 + x707 * x956 - x8 * x975
    x1562 = x8 * x987
    x1563 = x707 * x981
    x1564 = 3.0 * x898
    x1565 = x574 + x707 * x979 - x8 * x981 + 3.0 * x887
    x1566 = x4 * (x1562 - x1563 - x1564 + x1565 + x604)
    x1567 = -x1562 + x1563 + x1564 + x588
    x1568 = x601 + x707 * x987 - x8 * x991 + 3.0 * x913
    x1569 = x1003 * x8
    x1570 = x707 * x997
    x1571 = 2.0 * x928
    x1572 = -x1571
    x1573 = 2.0 * x917
    x1574 = x1573 + x610 + x707 * x995 - x8 * x997
    x1575 = x4 * (x1569 - x1570 + x1572 + x1574 + x638)
    x1576 = -x1569 + x1570 + x1571 + x623
    x1577 = 2.0 * x943
    x1578 = x1003 * x707 - x1007 * x8 + x1577 + x635
    x1579 = x1017 * x8
    x1580 = x1012 * x707
    x1581 = x1010 * x707 - x1012 * x8 + x643 + x948
    x1582 = x4 * (x1579 - x1580 + x1581 + x669 + x960)
    x1583 = -x1579 + x1580 + x655 + x959
    x1584 = x1017 * x707 - x1021 * x8 + x666 + x974
    x1585 = x1031 * x8
    x1586 = x1026 * x707
    x1587 = x1024 * x707 - x1026 * x8 + x676
    x1588 = x4 * (x1585 - x1586 + x1587 + x706)
    x1589 = -x1585 + x1586 + x690
    x1590 = x1031 * x707 - x1035 * x8 + x703
    x1591 = x1041 * x707
    x1592 = x1038 * x707
    x1593 = x1112 + x1592
    x1594 = x4 * (x1116 - x1591 + x1593)
    x1595 = x1593 * x5
    x1596 = x1117 + x1591
    x1597 = x1596 * x6
    x1598 = x1595 - x1597
    x1599 = x1598 * x5
    x1600 = x1596 * x5
    x1601 = x1048 * x707
    x1602 = x1133 + x1601
    x1603 = x1602 * x6
    x1604 = x1600 - x1603
    x1605 = x1604 * x6
    x1606 = x1594 + x1599 - x1605
    x1607 = x1606 * x5
    x1608 = x1055 * x707 + x1123
    x1609 = x4 * (x1111 - x1592 + x1608)
    x1610 = -x1593 * x6 + x1608 * x5
    x1611 = -x1598 * x6 + x1609 + x1610 * x5
    x1612 = -x1606 * x6 + x99 * (-x1595 + x1597 + x1610)
    x1613 = x1611 * x5 + x1612
    x1614 = x4 * (x1132 + x1596 - x1601)
    x1615 = x1604 * x5
    x1616 = x1602 * x5
    x1617 = x1066 * x707
    x1618 = x1151 + x1617
    x1619 = x1618 * x6
    x1620 = x1616 - x1619
    x1621 = x1620 * x6
    x1622 = x1614 + x1615 - x1621
    x1623 = x1622 * x6
    x1624 = x99 * (x1598 - x1600 + x1603)
    x1625 = x1623 - x1624
    x1626 = x1607 - x1623 + x1624
    x1627 = -x1594
    x1628 = x4 * (x1150 + x1602 - x1617)
    x1629 = x1086 * x707 + x1167
    x1630 = -x1614
    x1631 = x1113 * x707
    x1632 = x1230 + x1631
    x1633 = x1632 * x5
    x1634 = x1118 * x707
    x1635 = x1236 + x1634
    x1636 = x1635 * x6
    x1637 = x1633 - x1636
    x1638 = x1637 * x5
    x1639 = x1124 * x707 + x1225
    x1640 = -x1632 * x6 + x1639 * x5
    x1641 = x4 * (x1228 - x1631 + x1639)
    x1642 = -x1637 * x6 + x1641
    x1643 = x1640 * x5 + x1642
    x1644 = x1635 * x5
    x1645 = x1134 * x707
    x1646 = x1241 + x1645
    x1647 = x1646 * x6
    x1648 = x1644 - x1647
    x1649 = x1648 * x6
    x1650 = x4 * (x1245 + x1632 - x1634)
    x1651 = x1649 - x1650
    x1652 = x1638 - x1649 + x1650
    x1653 = x4 * (x1261 + x1635 - x1645)
    x1654 = x1152 * x707 + x1263
    x1655 = 120.679557322504 * x140
    x1656 = x1173 * x707
    x1657 = x1270 + x1656
    x1658 = x1657 * x5
    x1659 = x1176 * x707
    x1660 = x1276 + x1659
    x1661 = x1660 * x6
    x1662 = x1658 - x1661
    x1663 = x1662 * x5
    x1664 = x1181 * x707 + x1267
    x1665 = -x1657 * x6 + x1664 * x5
    x1666 = x4 * (x1265 - x1656 + x1664)
    x1667 = -x1662 * x6 + x1666
    x1668 = x1665 * x5 + x1667
    x1669 = x1660 * x5
    x1670 = x1189 * x707
    x1671 = x1281 + x1670
    x1672 = x1671 * x6
    x1673 = x1669 - x1672
    x1674 = x1673 * x6
    x1675 = x4 * (x1275 + x1657 - x1659)
    x1676 = -x1675
    x1677 = x1674 + x1676
    x1678 = x1663 - x1674 + x1675
    x1679 = x4 * (x1280 + x1660 - x1670)
    x1680 = x1205 * x707 + x1301
    x1681 = x1226 * x707 + x1334
    x1682 = x1237 * x707
    x1683 = x1347 + x1682
    x1684 = x1683 * x6
    x1685 = x1231 * x707
    x1686 = x1340 + x1685
    x1687 = -x1686 * x6
    x1688 = x1686 * x5
    x1689 = x1681 * x5 + x1687
    x1690 = x4 * (x1338 + x1681 - x1685)
    x1691 = -x1684 + x1688
    x1692 = x4 * (x1351 - x1682 + x1686)
    x1693 = x1242 * x707 + x1353
    x1694 = x1268 * x707 + x1356
    x1695 = x1277 * x707
    x1696 = x1367 + x1695
    x1697 = x1696 * x6
    x1698 = x1271 * x707
    x1699 = x1361 + x1698
    x1700 = -x1699 * x6
    x1701 = x1699 * x5
    x1702 = x1694 * x5 + x1700
    x1703 = x4 * (x1359 + x1694 - x1698)
    x1704 = -x1697 + x1701
    x1705 = x4 * (x1371 - x1695 + x1699)
    x1706 = x1282 * x707 + x1373
    x1707 = 209.023124717498 * x140
    x1708 = x1304 * x707 + x1377
    x1709 = x1310 * x707
    x1710 = x1386 + x1709
    x1711 = x1710 * x6
    x1712 = x1306 * x707
    x1713 = x1380 + x1712
    x1714 = -x1713 * x6
    x1715 = x1713 * x5
    x1716 = x1708 * x5 + x1714
    x1717 = x4 * (x1375 + x1708 - x1712)
    x1718 = -x1711 + x1715
    x1719 = x4 * (x1385 - x1709 + x1713)
    x1720 = x1313 * x707 + x1391
    x1721 = x1348 * x8
    x1722 = x1341 * x707
    x1723 = 3.0 * x1246
    x1724 = 3.0 * x1229 + x1335 * x707 - x1341 * x8
    x1725 = x4 * (x1721 - x1722 - x1723 + x1724)
    x1726 = -x1721 + x1722 + x1723
    x1727 = 3.0 * x1262 + x1348 * x707 - x1354 * x8
    x1728 = x1368 * x8
    x1729 = x1362 * x707
    x1730 = 2.0 * x1285
    x1731 = 2.0 * x1269 + x1357 * x707 - x1362 * x8
    x1732 = x4 * (x1728 - x1729 - x1730 + x1731)
    x1733 = -x1728 + x1729 + x1730
    x1734 = 2.0 * x1300 + x1368 * x707 - x1374 * x8
    x1735 = x1387 * x8
    x1736 = x1381 * x707
    x1737 = x1305 + x1378 * x707 - x1381 * x8
    x1738 = x4 * (x1317 + x1735 - x1736 + x1737)
    x1739 = x1316 - x1735 + x1736
    x1740 = x1331 + x1387 * x707 - x1392 * x8
    x1741 = x1404 * x8
    x1742 = x1398 * x707
    x1743 = x1396 * x707 - x1398 * x8
    x1744 = x4 * (x1741 - x1742 + x1743)
    x1745 = -x1741 + x1742
    x1746 = x1404 * x707 - x1408 * x8
    x1747 = x10 * x1048
    x1748 = x1036 * x1041
    x1749 = x1036 * x1038
    x1750 = x10 * x1041
    x1751 = x1749 - x1750 + x20
    x1752 = x4 * (x1747 - x1748 + x1751 + x75)
    x1753 = x1751 * x5
    x1754 = -x1747 + x1748 + x40
    x1755 = x1754 * x6
    x1756 = x1753 - x1755
    x1757 = x1756 * x5
    x1758 = x1754 * x5
    x1759 = x1036 * x1048
    x1760 = x10 * x1066
    x1761 = x1759 - x1760 + x61
    x1762 = x1761 * x6
    x1763 = x1758 - x1762
    x1764 = x1763 * x6
    x1765 = x1752 + x1757 - x1764
    x1766 = x1765 * x5
    x1767 = -x10 * x1038 + x1036 * x1055 + x33
    x1768 = x4 * (-x1749 + x1750 + x1767 + x54)
    x1769 = -x1751 * x6 + x1767 * x5
    x1770 = -x1756 * x6 + x1768 + x1769 * x5
    x1771 = -x1765 * x6 + x99 * (-x1753 + x1755 + x1769)
    x1772 = x1770 * x5 + x1771
    x1773 = x4 * (x104 + x1754 - x1759 + x1760)
    x1774 = x1763 * x5
    x1775 = x1761 * x5
    x1776 = x1036 * x1066
    x1777 = x10 * x1086
    x1778 = x1776 - x1777 + x90
    x1779 = x1778 * x6
    x1780 = x1775 - x1779
    x1781 = x1780 * x6
    x1782 = x1773 + x1774 - x1781
    x1783 = x1782 * x6
    x1784 = x99 * (x1756 - x1758 + x1762)
    x1785 = x1783 - x1784
    x1786 = x1766 - x1783 + x1784
    x1787 = -x1752
    x1788 = x4 * (x127 + x1761 - x1776 + x1777)
    x1789 = -x10 * x1105 + x1036 * x1086 + x125
    x1790 = -x1773
    x1791 = x144 * x1751
    x1792 = x1754 * x8
    x1793 = -x1792
    x1794 = x1791 + x1793
    x1795 = x1794 * x5
    x1796 = x144 * x1754
    x1797 = x1761 * x8
    x1798 = -x1797
    x1799 = x1796 + x1798
    x1800 = x1799 * x6
    x1801 = x1795 - x1800
    x1802 = x1801 * x5
    x1803 = -x1751 * x8
    x1804 = x144 * x1767 + x1803
    x1805 = -x1794 * x6 + x1804 * x5
    x1806 = x4 * (-x1791 + x1792 + x1804)
    x1807 = -x1801 * x6 + x1806
    x1808 = x1805 * x5 + x1807
    x1809 = x1799 * x5
    x1810 = x144 * x1761
    x1811 = x1778 * x8
    x1812 = -x1811
    x1813 = x1810 + x1812
    x1814 = x1813 * x6
    x1815 = x1809 - x1814
    x1816 = x1815 * x6
    x1817 = x4 * (x1794 - x1796 + x1797)
    x1818 = -x1817
    x1819 = x1816 + x1818
    x1820 = x1802 - x1816 + x1817
    x1821 = x4 * (x1799 - x1810 + x1811)
    x1822 = -x1789 * x8
    x1823 = x144 * x1778 + x1822
    x1824 = x1036 * x1173
    x1825 = x10 * x1176
    x1826 = x1045 + x1824 - x1825 + x282
    x1827 = x1826 * x5
    x1828 = x1036 * x1176
    x1829 = x10 * x1189
    x1830 = x1062 + x1828 - x1829 + x291
    x1831 = x1830 * x6
    x1832 = x1827 - x1831
    x1833 = x1832 * x5
    x1834 = -x10 * x1173 + x1036 * x1181 + x1057 + x270
    x1835 = -x1826 * x6 + x1834 * x5
    x1836 = x4 * (x1076 - x1824 + x1825 + x1834 + x283)
    x1837 = -x1832 * x6 + x1836
    x1838 = x1835 * x5 + x1837
    x1839 = x1830 * x5
    x1840 = x1036 * x1189
    x1841 = x10 * x1205
    x1842 = x1082 + x1840 - x1841 + x317
    x1843 = x1842 * x6
    x1844 = x1839 - x1843
    x1845 = x1844 * x6
    x1846 = x4 * (x1096 + x1826 - x1828 + x1829 + x311)
    x1847 = -x1846
    x1848 = x1845 + x1847
    x1849 = x1833 - x1845 + x1846
    x1850 = x4 * (x1106 + x1830 - x1840 + x1841 + x346)
    x1851 = -x10 * x1220 + x1036 * x1205 + x1104 + x348
    x1852 = x1768 - x1794 * x8
    x1853 = x144 * x1804 + x1852
    x1854 = x144 * x1799
    x1855 = x1813 * x8
    x1856 = x1773 - x1855
    x1857 = x1854 + x1856
    x1858 = x1857 * x6
    x1859 = x144 * x1794
    x1860 = x1799 * x8
    x1861 = x1752 - x1860
    x1862 = x1859 + x1861
    x1863 = -x1862 * x6
    x1864 = x1862 * x5
    x1865 = x1853 * x5 + x1863
    x1866 = x1787 + x1860
    x1867 = x4 * (x1853 - x1859 + x1866)
    x1868 = -x1858 + x1864
    x1869 = x1790 + x1855
    x1870 = x4 * (-x1854 + x1862 + x1869)
    x1871 = x1788 - x1823 * x8
    x1872 = x144 * x1813 + x1871
    x1873 = -x1826 * x8
    x1874 = x144 * x1834 + x1873
    x1875 = x144 * x1830
    x1876 = x1842 * x8
    x1877 = -x1876
    x1878 = x1875 + x1877
    x1879 = x1878 * x6
    x1880 = x144 * x1826
    x1881 = x1830 * x8
    x1882 = -x1881
    x1883 = x1880 + x1882
    x1884 = -x1883 * x6
    x1885 = x1883 * x5
    x1886 = x1874 * x5 + x1884
    x1887 = x4 * (x1874 - x1880 + x1881)
    x1888 = -x1879 + x1885
    x1889 = x4 * (-x1875 + x1876 + x1883)
    x1890 = -x1851 * x8
    x1891 = x144 * x1842 + x1890
    x1892 = 2.0 * x1184
    x1893 = -x10 * x1306 + x1036 * x1304 + x1892 + x503
    x1894 = x1036 * x1310
    x1895 = x10 * x1313
    x1896 = 2.0 * x1202
    x1897 = x1894 - x1895 + x1896 + x520
    x1898 = x1897 * x6
    x1899 = x1036 * x1306
    x1900 = x10 * x1310
    x1901 = 2.0 * x1193
    x1902 = x1899 - x1900 + x1901 + x537
    x1903 = -x1902 * x6
    x1904 = x1902 * x5
    x1905 = x1893 * x5 + x1903
    x1906 = -x1901
    x1907 = x4 * (x1893 - x1899 + x1900 + x1906 + x543)
    x1908 = -x1898 + x1904
    x1909 = -x1896
    x1910 = x4 * (-x1894 + x1895 + x1902 + x1909 + x563)
    x1911 = 2.0 * x1219
    x1912 = -x10 * x1332 + x1036 * x1313 + x1911 + x565
    x1913 = x144 * x1862
    x1914 = 2.0 * x1806 - x1862 * x8
    x1915 = x144 * x1853 + x1914
    x1916 = x1857 * x8
    x1917 = 2.0 * x1817
    x1918 = x1916 - x1917
    x1919 = x4 * (-x1913 + x1915 + x1918)
    x1920 = -x1916 + x1917
    x1921 = x1913 + x1920
    x1922 = 2.0 * x1821 - x1872 * x8
    x1923 = x144 * x1857 + x1922
    x1924 = x144 * x1883
    x1925 = x1836 - x1883 * x8
    x1926 = x144 * x1874 + x1925
    x1927 = x1878 * x8
    x1928 = x1847 + x1927
    x1929 = x4 * (-x1924 + x1926 + x1928)
    x1930 = x1846 - x1927
    x1931 = x1924 + x1930
    x1932 = x1850 - x1891 * x8
    x1933 = x144 * x1878 + x1932
    x1934 = x1897 * x8
    x1935 = x144 * x1902
    x1936 = -x1902 * x8
    x1937 = x144 * x1893 + x1936
    x1938 = x4 * (x1934 - x1935 + x1937)
    x1939 = -x1934
    x1940 = x1935 + x1939
    x1941 = -x1912 * x8
    x1942 = x144 * x1897 + x1941
    x1943 = x10 * x1404
    x1944 = x1036 * x1398
    x1945 = 3.0 * x1316
    x1946 = -x10 * x1398 + x1036 * x1396 + 3.0 * x1305 + x676
    x1947 = x4 * (x1943 - x1944 - x1945 + x1946 + x706)
    x1948 = -x1943 + x1944 + x1945 + x690
    x1949 = -x10 * x1408 + x1036 * x1404 + 3.0 * x1331 + x703
    x1950 = x1416 * x8
    x1951 = x1413 * x707
    x1952 = 2.0 * x716
    x1953 = -x1413 * x8 + x1429 * x707 + 2.0 * x728
    x1954 = x4 * (x1950 - x1951 - x1952 + x1953)
    x1955 = -x1950 + x1951 + x1952
    x1956 = x1953 * x5 - x1955 * x6
    x1957 = x1955 * x5
    x1958 = x1416 * x707
    x1959 = x1423 * x8
    x1960 = 2.0 * x733
    x1961 = x1958 - x1959 + x1960
    x1962 = x1961 * x6
    x1963 = x1957 - x1962
    x1964 = x1954 + x1956 * x5 - x1963 * x6
    x1965 = x4 * (x1955 - x1958 + x1959 - x1960)
    x1966 = x1963 * x5
    x1967 = x1961 * x5
    x1968 = x1423 * x707
    x1969 = x1440 * x8
    x1970 = 2.0 * x753
    x1971 = x1968 - x1969 + x1970
    x1972 = x1971 * x6
    x1973 = x1967 - x1972
    x1974 = x1973 * x6
    x1975 = x1965 + x1966 - x1974
    x1976 = x4 * (x1961 - x1968 + x1969 - x1970)
    x1977 = x1440 * x707 - x1451 * x8 + 2.0 * x775
    x1978 = -x1965
    x1979 = x1460 * x8
    x1980 = x1456 * x707
    x1981 = x1430 - x1456 * x8 + x1464 * x707 + x1510
    x1982 = x4 * (x1449 + x1524 + x1979 - x1980 + x1981)
    x1983 = x1414 + x1519 - x1979 + x1980
    x1984 = x1981 * x5 - x1983 * x6
    x1985 = x1983 * x5
    x1986 = x1460 * x707
    x1987 = x1472 * x8
    x1988 = x1435 + x1514 + x1986 - x1987
    x1989 = x1988 * x6
    x1990 = x1985 - x1989
    x1991 = x4 * (x1452 + x1527 + x1983 - x1986 + x1987)
    x1992 = x1450 + x1472 * x707 - x1480 * x8 + x1529
    x1993 = x1488 * x8
    x1994 = x1484 * x707
    x1995 = 2.0 * x854
    x1996 = -x1484 * x8 + x1492 * x707 + 2.0 * x845
    x1997 = x4 * (x1993 - x1994 - x1995 + x1996)
    x1998 = -x1993 + x1994 + x1995
    x1999 = x1996 * x5 - x1998 * x6
    x2000 = x1998 * x5
    x2001 = x1488 * x707
    x2002 = x1500 * x8
    x2003 = 2.0 * x863
    x2004 = x2001 - x2002 + x2003
    x2005 = x2004 * x6
    x2006 = x2000 - x2005
    x2007 = x4 * (x1998 - x2001 + x2002 - x2003)
    x2008 = x1500 * x707 - x1509 * x8 + 2.0 * x880
    x2009 = x1515 * x8
    x2010 = x1520 * x707
    x2011 = 2.0 * x898
    x2012 = 2.0 * x1476
    x2013 = 2.0 * x1466 + x1511 * x707 - x1520 * x8 + 2.0 * x887
    x2014 = x4 * (x2009 - x2010 - x2011 - x2012 + x2013)
    x2015 = -x2009 + x2010 + x2011 + x2012
    x2016 = 2.0 * x1479 + x1515 * x707 - x1530 * x8 + 2.0 * x913
    x2017 = x1534 * x8
    x2018 = x1538 * x707
    x2019 = x1494 + x1531 * x707 - x1538 * x8 + x1573
    x2020 = x4 * (x1505 + x1572 + x2017 - x2018 + x2019)
    x2021 = x1504 + x1571 - x2017 + x2018
    x2022 = x1508 + x1534 * x707 - x1545 * x8 + x1577
    x2023 = x1550 * x8
    x2024 = x1554 * x707
    x2025 = 2.0 * x959
    x2026 = x1547 * x707 - x1554 * x8 + 2.0 * x948
    x2027 = x4 * (x2023 - x2024 - x2025 + x2026)
    x2028 = -x2023 + x2024 + x2025
    x2029 = x1550 * x707 - x1561 * x8 + 2.0 * x974
    x2030 = 3.0 * x1525
    x2031 = x1565 * x707 - x1567 * x8 + x2030 + 2.0 * x980
    x2032 = 3.0 * x1528
    x2033 = x1567 * x707 - x1568 * x8 + x2032 + 2.0 * x990
    x2034 = 2.0 * x1542 + x1574 * x707 - x1576 * x8 + 2.0 * x996
    x2035 = 2.0 * x1006 + 2.0 * x1544 + x1576 * x707 - x1578 * x8
    x2036 = 2.0 * x1011 + x1558 + x1581 * x707 - x1583 * x8
    x2037 = 2.0 * x1020 + x1560 + x1583 * x707 - x1584 * x8
    x2038 = 2.0 * x1025 + x1587 * x707 - x1589 * x8
    x2039 = 2.0 * x1034 + x1589 * x707 - x1590 * x8
    x2040 = x1596 * x8
    x2041 = x1593 * x707
    x2042 = x1057 - x1593 * x8 + x1608 * x707
    x2043 = x4 * (x1076 + x2040 - x2041 + x2042)
    x2044 = x1045 - x2040 + x2041
    x2045 = x2042 * x5 - x2044 * x6
    x2046 = x2044 * x5
    x2047 = x1596 * x707
    x2048 = x1602 * x8
    x2049 = x1062 + x2047 - x2048
    x2050 = x2049 * x6
    x2051 = x2046 - x2050
    x2052 = x2043 + x2045 * x5 - x2051 * x6
    x2053 = x4 * (x1096 + x2044 - x2047 + x2048)
    x2054 = x2051 * x5
    x2055 = x2049 * x5
    x2056 = x1602 * x707
    x2057 = x1618 * x8
    x2058 = x1082 + x2056 - x2057
    x2059 = x2058 * x6
    x2060 = x2055 - x2059
    x2061 = x2060 * x6
    x2062 = x2053 + x2054 - x2061
    x2063 = x4 * (x1106 + x2049 - x2056 + x2057)
    x2064 = x1104 + x1618 * x707 - x1629 * x8
    x2065 = -x2053
    x2066 = x1635 * x8
    x2067 = x1632 * x707
    x2068 = x1127 + x1609 - x1632 * x8 + x1639 * x707
    x2069 = x4 * (x1139 + x1627 + x2066 - x2067 + x2068)
    x2070 = x1138 + x1594 - x2066 + x2067
    x2071 = x2068 * x5 - x2070 * x6
    x2072 = x2070 * x5
    x2073 = x1635 * x707
    x2074 = x1646 * x8
    x2075 = x1147 + x1614 + x2073 - x2074
    x2076 = x2075 * x6
    x2077 = x2072 - x2076
    x2078 = x4 * (x1164 + x1630 + x2070 - x2073 + x2074)
    x2079 = x1166 + x1628 + x1646 * x707 - x1654 * x8
    x2080 = x1660 * x8
    x2081 = x1657 * x707
    x2082 = x1184 - x1657 * x8 + x1664 * x707
    x2083 = x4 * (x1194 + x2080 - x2081 + x2082)
    x2084 = x1193 - x2080 + x2081
    x2085 = x2082 * x5 - x2084 * x6
    x2086 = x2084 * x5
    x2087 = x1660 * x707
    x2088 = x1671 * x8
    x2089 = x1202 + x2087 - x2088
    x2090 = x2089 * x6
    x2091 = x2086 - x2090
    x2092 = x4 * (x1217 + x2084 - x2087 + x2088)
    x2093 = x1219 + x1671 * x707 - x1680 * x8
    x2094 = x1683 * x8
    x2095 = x1686 * x707
    x2096 = 2.0 * x1650
    x2097 = -x2096
    x2098 = 2.0 * x1641
    x2099 = x1229 + x1681 * x707 - x1686 * x8 + x2098
    x2100 = x4 * (x1247 + x2094 - x2095 + x2097 + x2099)
    x2101 = x1246 - x2094 + x2095 + x2096
    x2102 = 2.0 * x1653
    x2103 = x1262 + x1683 * x707 - x1693 * x8 + x2102
    x2104 = x1696 * x8
    x2105 = x1699 * x707
    x2106 = x1269 + x1666 + x1694 * x707 - x1699 * x8
    x2107 = x4 * (x1286 + x1676 + x2104 - x2105 + x2106)
    x2108 = x1285 + x1675 - x2104 + x2105
    x2109 = x1300 + x1679 + x1696 * x707 - x1706 * x8
    x2110 = x1710 * x8
    x2111 = x1713 * x707
    x2112 = x1305 + x1708 * x707 - x1713 * x8
    x2113 = x4 * (x1317 + x2110 - x2111 + x2112)
    x2114 = x1316 - x2110 + x2111
    x2115 = x1331 + x1710 * x707 - x1720 * x8
    x2116 = x1339 + 3.0 * x1690 + x1724 * x707 - x1726 * x8
    x2117 = x1352 + 3.0 * x1692 + x1726 * x707 - x1727 * x8
    x2118 = 2.0 * x1703
    x2119 = x1360 + x1731 * x707 - x1733 * x8 + x2118
    x2120 = 2.0 * x1705
    x2121 = x1372 + x1733 * x707 - x1734 * x8 + x2120
    x2122 = x1379 + x1717 + x1737 * x707 - x1739 * x8
    x2123 = x1390 + x1719 + x1739 * x707 - x1740 * x8
    x2124 = x1397 + x1743 * x707 - x1745 * x8
    x2125 = x1407 + x1745 * x707 - x1746 * x8
    x2126 = x1751 * x707
    x2127 = x1767 * x707 + x1803
    x2128 = x4 * (x1792 - x2126 + x2127)
    x2129 = x1793 + x2126
    x2130 = x2127 * x5 - x2129 * x6
    x2131 = x2129 * x5
    x2132 = x1754 * x707
    x2133 = x1798 + x2132
    x2134 = x2133 * x6
    x2135 = x2131 - x2134
    x2136 = x2128 + x2130 * x5 - x2135 * x6
    x2137 = x4 * (x1797 + x2129 - x2132)
    x2138 = x2135 * x5
    x2139 = x2133 * x5
    x2140 = x1761 * x707
    x2141 = x1812 + x2140
    x2142 = x2141 * x6
    x2143 = x2139 - x2142
    x2144 = x2143 * x6
    x2145 = x2137 + x2138 - x2144
    x2146 = x4 * (x1811 + x2133 - x2140)
    x2147 = x1778 * x707 + x1822
    x2148 = -x2137
    x2149 = x1794 * x707
    x2150 = x1804 * x707 + x1852
    x2151 = x4 * (x1866 - x2149 + x2150)
    x2152 = x1861 + x2149
    x2153 = x2150 * x5 - x2152 * x6
    x2154 = x2152 * x5
    x2155 = x1799 * x707
    x2156 = x1856 + x2155
    x2157 = x2156 * x6
    x2158 = x2154 - x2157
    x2159 = x4 * (x1869 + x2152 - x2155)
    x2160 = x1813 * x707 + x1871
    x2161 = x1826 * x707
    x2162 = x1834 * x707 + x1873
    x2163 = x4 * (x1881 - x2161 + x2162)
    x2164 = x1882 + x2161
    x2165 = x2162 * x5 - x2164 * x6
    x2166 = x2164 * x5
    x2167 = x1830 * x707
    x2168 = x1877 + x2167
    x2169 = x2168 * x6
    x2170 = x2166 - x2169
    x2171 = x4 * (x1876 + x2164 - x2167)
    x2172 = x1842 * x707 + x1890
    x2173 = x1862 * x707
    x2174 = x1853 * x707 + x1914
    x2175 = x4 * (x1918 - x2173 + x2174)
    x2176 = x1920 + x2173
    x2177 = x1857 * x707 + x1922
    x2178 = x1883 * x707
    x2179 = x1874 * x707 + x1925
    x2180 = x4 * (x1928 - x2178 + x2179)
    x2181 = x1930 + x2178
    x2182 = x1878 * x707 + x1932
    x2183 = x1902 * x707
    x2184 = x1893 * x707 + x1936
    x2185 = x4 * (x1934 - x2183 + x2184)
    x2186 = x1939 + x2183
    x2187 = x1897 * x707 + x1941
    x2188 = 3.0 * x1867 + x1915 * x707 - x1921 * x8
    x2189 = 3.0 * x1870 + x1921 * x707 - x1923 * x8
    x2190 = 2.0 * x1887 + x1926 * x707 - x1931 * x8
    x2191 = 2.0 * x1889 + x1931 * x707 - x1933 * x8
    x2192 = x1907 + x1937 * x707 - x1940 * x8
    x2193 = x1910 + x1940 * x707 - x1942 * x8
    x2194 = x1946 * x707 - x1948 * x8
    x2195 = x1948 * x707 - x1949 * x8
    x2196 = x10 * x1754
    x2197 = x1036 * x1751
    x2198 = 2.0 * x1045
    x2199 = -x10 * x1751 + x1036 * x1767 + 2.0 * x1057
    x2200 = x4 * (x2196 - x2197 - x2198 + x2199)
    x2201 = -x2196 + x2197 + x2198
    x2202 = x2199 * x5 - x2201 * x6
    x2203 = x2201 * x5
    x2204 = x1036 * x1754
    x2205 = x10 * x1761
    x2206 = 2.0 * x1062
    x2207 = x2204 - x2205 + x2206
    x2208 = x2207 * x6
    x2209 = x2203 - x2208
    x2210 = x2200 + x2202 * x5 - x2209 * x6
    x2211 = x4 * (x2201 - x2204 + x2205 - x2206)
    x2212 = x2209 * x5
    x2213 = x2207 * x5
    x2214 = x1036 * x1761
    x2215 = x10 * x1778
    x2216 = 2.0 * x1082
    x2217 = x2214 - x2215 + x2216
    x2218 = x2217 * x6
    x2219 = x2213 - x2218
    x2220 = x2219 * x6
    x2221 = x2211 + x2212 - x2220
    x2222 = x4 * (x2207 - x2214 + x2215 - x2216)
    x2223 = -x10 * x1789 + x1036 * x1778 + 2.0 * x1104
    x2224 = -x2211
    x2225 = x2207 * x8
    x2226 = x144 * x2201
    x2227 = -x2201 * x8
    x2228 = x144 * x2199 + x2227
    x2229 = x4 * (x2225 - x2226 + x2228)
    x2230 = -x2225
    x2231 = x2226 + x2230
    x2232 = x2228 * x5 - x2231 * x6
    x2233 = x2231 * x5
    x2234 = x144 * x2207
    x2235 = x2217 * x8
    x2236 = -x2235
    x2237 = x2234 + x2236
    x2238 = x2237 * x6
    x2239 = x2233 - x2238
    x2240 = x4 * (x2231 - x2234 + x2235)
    x2241 = -x2223 * x8
    x2242 = x144 * x2217 + x2241
    x2243 = x10 * x1830
    x2244 = x1036 * x1826
    x2245 = -x10 * x1826 + x1036 * x1834 + x1768 + x1892
    x2246 = x4 * (x1787 + x1906 + x2243 - x2244 + x2245)
    x2247 = x1752 + x1901 - x2243 + x2244
    x2248 = x2245 * x5 - x2247 * x6
    x2249 = x2247 * x5
    x2250 = x1036 * x1830
    x2251 = x10 * x1842
    x2252 = x1773 + x1896 + x2250 - x2251
    x2253 = x2252 * x6
    x2254 = x2249 - x2253
    x2255 = x4 * (x1790 + x1909 + x2247 - x2250 + x2251)
    x2256 = -x10 * x1851 + x1036 * x1842 + x1788 + x1911
    x2257 = x144 * x2231
    x2258 = x2200 - x2231 * x8
    x2259 = x144 * x2228 + x2258
    x2260 = x2237 * x8
    x2261 = x2224 + x2260
    x2262 = x4 * (-x2257 + x2259 + x2261)
    x2263 = x2211 - x2260
    x2264 = x2257 + x2263
    x2265 = x2222 - x2242 * x8
    x2266 = x144 * x2237 + x2265
    x2267 = x2252 * x8
    x2268 = x144 * x2247
    x2269 = -x2247 * x8
    x2270 = x144 * x2245 + x2269
    x2271 = x4 * (x2267 - x2268 + x2270)
    x2272 = -x2267
    x2273 = x2268 + x2272
    x2274 = -x2256 * x8
    x2275 = x144 * x2252 + x2274
    x2276 = x10 * x1897
    x2277 = x1036 * x1902
    x2278 = 2.0 * x1316
    x2279 = 2.0 * x1846
    x2280 = -x10 * x1902 + x1036 * x1893 + 2.0 * x1305 + 2.0 * x1836
    x2281 = x4 * (x2276 - x2277 - x2278 - x2279 + x2280)
    x2282 = -x2276 + x2277 + x2278 + x2279
    x2283 = -x10 * x1912 + x1036 * x1897 + 2.0 * x1331 + 2.0 * x1850
    x2284 = 2.0 * x2229 - x2264 * x8
    x2285 = x144 * x2259 + x2284
    x2286 = 2.0 * x2240 - x2266 * x8
    x2287 = x144 * x2264 + x2286
    x2288 = x2246 - x2273 * x8
    x2289 = x144 * x2270 + x2288
    x2290 = x2255 - x2275 * x8
    x2291 = x144 * x2273 + x2290
    x2292 = -x2282 * x8
    x2293 = x144 * x2280 + x2292
    x2294 = -x2283 * x8
    x2295 = x144 * x2282 + x2294
    x2296 = 3.0 * x1907
    x2297 = -x10 * x1948 + x1036 * x1946 + 2.0 * x1397 + x2296
    x2298 = 3.0 * x1910
    x2299 = -x10 * x1949 + x1036 * x1948 + 2.0 * x1407 + x2298
    x2300 = x1961 * x8
    x2301 = x1955 * x707
    x2302 = 3.0 * x1414
    x2303 = 3.0 * x1430 + x1953 * x707 - x1955 * x8
    x2304 = -x2300 + x2301 + x2302
    x2305 = x2303 * x5 - x2304 * x6
    x2306 = x2304 * x5
    x2307 = x1961 * x707
    x2308 = x1971 * x8
    x2309 = 3.0 * x1435
    x2310 = x2307 - x2308 + x2309
    x2311 = x2310 * x6
    x2312 = x2306 - x2311
    x2313 = x1988 * x8
    x2314 = x1983 * x707
    x2315 = 3.0 * x1476
    x2316 = 3.0 * x1466 + x1954 + x1981 * x707 - x1983 * x8
    x2317 = x1965 - x2313 + x2314 + x2315
    x2318 = x2004 * x8
    x2319 = x1998 * x707
    x2320 = 3.0 * x1504
    x2321 = 3.0 * x1494 + x1996 * x707 - x1998 * x8
    x2322 = -x2318 + x2319 + x2320
    x2323 = x2049 * x8
    x2324 = x2044 * x707
    x2325 = 2.0 * x1594
    x2326 = 2.0 * x1609 + x2042 * x707 - x2044 * x8
    x2327 = -x2323 + x2324 + x2325
    x2328 = x2326 * x5 - x2327 * x6
    x2329 = x2327 * x5
    x2330 = x2049 * x707
    x2331 = x2058 * x8
    x2332 = 2.0 * x1614
    x2333 = x2330 - x2331 + x2332
    x2334 = x2333 * x6
    x2335 = x2329 - x2334
    x2336 = x2075 * x8
    x2337 = x2070 * x707
    x2338 = x2043 + x2068 * x707 - x2070 * x8 + x2098
    x2339 = x2053 + x2096 - x2336 + x2337
    x2340 = x2089 * x8
    x2341 = x2084 * x707
    x2342 = 2.0 * x1675
    x2343 = 2.0 * x1666 + x2082 * x707 - x2084 * x8
    x2344 = -x2340 + x2341 + x2342
    x2345 = x2133 * x8
    x2346 = x2129 * x707
    x2347 = x1768 + x2127 * x707 - x2129 * x8
    x2348 = x1752 - x2345 + x2346
    x2349 = x2347 * x5 - x2348 * x6
    x2350 = x2348 * x5
    x2351 = x2133 * x707
    x2352 = x2141 * x8
    x2353 = x1773 + x2351 - x2352
    x2354 = x2353 * x6
    x2355 = x2350 - x2354
    x2356 = x2156 * x8
    x2357 = x2152 * x707
    x2358 = x1806 + x2128 + x2150 * x707 - x2152 * x8
    x2359 = x1817 + x2137 - x2356 + x2357
    x2360 = x2168 * x8
    x2361 = x2164 * x707
    x2362 = x1836 + x2162 * x707 - x2164 * x8
    x2363 = x1846 - x2360 + x2361
    x2364 = x2201 * x707
    x2365 = x2199 * x707 + x2227
    x2366 = x2230 + x2364
    x2367 = x2365 * x5 - x2366 * x6
    x2368 = x2366 * x5
    x2369 = x2207 * x707
    x2370 = x2236 + x2369
    x2371 = x2370 * x6
    x2372 = x2368 - x2371
    x2373 = x2231 * x707
    x2374 = x2228 * x707 + x2258
    x2375 = x2263 + x2373
    x2376 = x2247 * x707
    x2377 = x2245 * x707 + x2269
    x2378 = x2272 + x2376
    x2379 = x10 * x2207
    x2380 = x1036 * x2201
    x2381 = 3.0 * x1752
    x2382 = -x10 * x2201 + x1036 * x2199 + 3.0 * x1768
    x2383 = x4 * (x2379 - x2380 - x2381 + x2382)
    x2384 = -x2379 + x2380 + x2381
    x2385 = x2382 * x5 - x2384 * x6
    x2386 = x2384 * x5
    x2387 = x1036 * x2207
    x2388 = x10 * x2217
    x2389 = 3.0 * x1773
    x2390 = x2387 - x2388 + x2389
    x2391 = x2390 * x6
    x2392 = x2386 - x2391
    x2393 = x4 * (x2384 - x2387 + x2388 - x2389)
    x2394 = -x10 * x2223 + x1036 * x2217 + 3.0 * x1788
    x2395 = x2390 * x8
    x2396 = x144 * x2384
    x2397 = x144 * x2382 - x2384 * x8
    x2398 = x4 * (x2395 - x2396 + x2397)
    x2399 = -x2395 + x2396
    x2400 = x144 * x2390 - x2394 * x8
    x2401 = x10 * x2252
    x2402 = x1036 * x2247
    x2403 = 3.0 * x1846
    x2404 = -x10 * x2247 + x1036 * x2245 + 3.0 * x1836 + x2200
    x2405 = x4 * (x2224 + x2401 - x2402 - x2403 + x2404)
    x2406 = x2211 - x2401 + x2402 + x2403
    x2407 = -x10 * x2256 + x1036 * x2252 + 3.0 * x1850 + x2222
    x2408 = x144 * x2397 + x2383 - x2399 * x8
    x2409 = x144 * x2399 + x2393 - x2400 * x8
    x2410 = x144 * x2404 - x2406 * x8
    x2411 = x144 * x2406 - x2407 * x8
    x2412 = -x10 * x2282 + x1036 * x2280 + 2.0 * x2246 + x2296
    x2413 = -x10 * x2283 + x1036 * x2282 + 2.0 * x2255 + x2298

    # 150 item(s)
    result[0, 0] = numpy.sum(
        -x141
        * (
            x3
            * (
                x114 * x6
                - x3 * x84
                + x83 * (-x116 + x117 + x121 - x123 + x77)
                - x99 * (x57 - x58 + x74 - x78)
            )
            + x6
            * (
                x114 * x3
                - x6
                * (
                    x108 * x3
                    + x4 * (x124 + x73 - x89)
                    - x6
                    * (
                        x102 * x3
                        - x6
                        * (
                            x5 * x97
                            - x6 * (x125 + x5 * x95 + x6 * (x126 * x6 - x5 * x92))
                            + x99 * (x65 + x93 - x94)
                        )
                        + x83 * (x127 + x68 - x91 + x96)
                    )
                    + x83 * (x115 + x124 - x128)
                )
                - x83 * (x106 - x129 + x130 + x133 - x134)
                + x99 * (x103 - x107 + x79 - x88)
            )
            + x83 * (-x109 + x113 - x84 + x86 + x87)
            + x83
            * (
                x112 * x99
                - x123 * x3
                + x134 * x3
                + x134 * x6
                - x6 * (x106 - x129 + x130 + x133)
                - x82 * x99
                + x99
                * (
                    x120 * x6
                    - x122 * x3
                    + x137
                    + x138
                    - x139
                    - x36
                    - x4 * (-x136 + x14 + x3 * x31 + x30)
                    + x49
                )
                - x99
                * (
                    x132 * x3
                    - x137
                    - x138
                    + x139
                    + x4 * (x135 - x21 * x3 + x23 + x42)
                    - x49
                    - x6 * (x131 + x67)
                    + x70
                )
            )
        )
    )
    result[0, 1] = numpy.sum(
        x251
        * (
            x3
            * (
                x200 * x3
                - x226 * x6
                + x99 * (x178 - x179 + x193 - x195)
                - x99 * (x194 - x228 + x232 + x233 - x235)
            )
            - x6
            * (
                x226 * x3
                - x6
                * (
                    x221 * x3
                    + x4 * (x192 - x206 + x237)
                    - x6
                    * (
                        x217 * x3
                        - x6 * (x214 * x5 + x238 - x6 * (x212 * x5 - x240 * x6))
                        + x99 * (x189 - x208 + x213)
                    )
                    + x99 * (x227 + x237 - x241)
                )
                + x99 * (x196 - x205 + x218 - x220)
                - x99 * (x219 - x242 + x245 + x246 - x247)
            )
            + x83 * (x200 - x203 - x204 + x222 + x250)
            + x99
            * (
                x199
                + x235 * x3
                - x247 * x3
                - x247 * x6
                + x250
                + x4 * (x160 + x173 - x231 * x6 + x234 * x3 + x248 - x249)
                - x4 * (x172 + x201 - x244 * x3 - x248 + x249 + x6 * (x188 + x243))
                + x6 * (x219 - x242 + x245 + x246)
            )
        )
    )
    result[0, 2] = numpy.sum(
        x251
        * (
            x3
            * (
                x3 * x310
                - x336 * x6
                + x99 * (x288 - x289 + x303 - x305)
                - x99 * (x304 - x338 + x342 + x343 - x345)
            )
            - x6
            * (
                x3 * x336
                - x6
                * (
                    x3 * x331
                    + x4 * (x302 - x316 + x347)
                    - x6
                    * (
                        x3 * x327
                        - x6 * (x324 * x5 + x348 - x6 * (x322 * x5 - x350 * x6))
                        + x99 * (x299 - x318 + x323)
                    )
                    + x99 * (x337 + x347 - x351)
                )
                + x99 * (x306 - x315 + x328 - x330)
                - x99 * (x329 - x352 + x355 + x356 - x357)
            )
            + x83 * (x310 - x313 - x314 + x332 + x360)
            + x99
            * (
                x3 * x345
                - x3 * x357
                + x309
                - x357 * x6
                + x360
                + x4 * (x270 + x283 + x3 * x344 - x341 * x6 + x358 - x359)
                - x4 * (x282 - x3 * x354 + x311 - x358 + x359 + x6 * (x298 + x353))
                + x6 * (x329 - x352 + x355 + x356)
            )
        )
    )
    result[0, 3] = numpy.sum(
        x251
        * (
            x3
            * (
                x3 * x416
                + x4 * (x379 - x401 + x402 + x406)
                - x425 * x6
                + x99 * (x406 + x410 - x411 + x414)
            )
            + x4
            * (
                x3 * x379
                - x3 * x403
                - x403 * x6
                - 2.0 * x404
                + 2.0 * x405
                + x6 * (x383 + x390 - x398)
            )
            - x6
            * (
                x3 * x425
                + x4 * (-x390 + x398 + x403 + x426)
                - x6
                * (
                    x3 * x423
                    + x4 * (x389 - x391 + x396)
                    + x4 * (x396 + x413 - x420)
                    - x6 * (x3 * x421 + x428 - x6 * (x395 * x5 - x430 * x6))
                )
                + x99 * (x415 - x419 + x422 + x426)
            )
            - x83 * (x404 - x416 + x417 + x418 - x424)
        )
    )
    result[0, 4] = numpy.sum(
        x497
        * (
            x3
            * (
                x3 * x483
                + x4 * (x448 - x468 + x469 + x473)
                - x492 * x6
                + x99 * (x473 + x477 - x478 + x481)
            )
            + x4
            * (
                x3 * x448
                - x3 * x470
                - x470 * x6
                - 2.0 * x471
                + 2.0 * x472
                + x6 * (x451 + x458 - x466)
            )
            - x6
            * (
                x3 * x492
                + x4 * (-x458 + x466 + x470 + x493)
                - x6
                * (
                    x3 * x490
                    + x4 * (x457 - x459 + x464)
                    + x4 * (x464 + x480 - x487)
                    - x6 * (x3 * x488 + x494 - x6 * (x463 * x5 - x496 * x6))
                )
                + x99 * (x482 - x486 + x489 + x493)
            )
            - x83 * (x471 - x483 + x484 + x485 - x491)
        )
    )
    result[0, 5] = numpy.sum(
        x251
        * (
            x3
            * (
                x3 * x553
                + x4 * (x516 - x538 + x539 + x543)
                - x562 * x6
                + x99 * (x543 + x547 - x548 + x551)
            )
            + x4
            * (
                x3 * x516
                - x3 * x540
                - x540 * x6
                - 2.0 * x541
                + 2.0 * x542
                + x6 * (x520 + x527 - x535)
            )
            - x6
            * (
                x3 * x562
                + x4 * (-x527 + x535 + x540 + x563)
                - x6
                * (
                    x3 * x560
                    + x4 * (x526 - x528 + x533)
                    + x4 * (x533 + x550 - x557)
                    - x6 * (x3 * x558 + x565 - x6 * (x5 * x532 - x567 * x6))
                )
                + x99 * (x552 - x556 + x559 + x563)
            )
            - x83 * (x541 - x553 + x554 + x555 - x561)
        )
    )
    result[0, 6] = numpy.sum(
        x141
        * (
            x3 * (x3 * x586 - x599 * x6 + x99 * (x577 - x578 + x584))
            - x6
            * (
                x3 * x599
                - x6 * (x3 * x597 - x6 * (x3 * x595 - x6 * x603) + x601)
                + x99 * (x585 - x590 + x596)
            )
            + x83 * (x586 - x589 + x598 + x604)
        )
    )
    result[0, 7] = numpy.sum(
        x251
        * (
            x3 * (x3 * x621 - x6 * x633 + x99 * (x613 - x614 + x619))
            - x6
            * (
                x3 * x633
                - x6 * (x3 * x631 - x6 * (x3 * x629 - x6 * x637) + x635)
                + x99 * (x620 - x625 + x630)
            )
            + x83 * (x621 - x624 + x632 + x638)
        )
    )
    result[0, 8] = numpy.sum(
        x251
        * (
            x3 * (x3 * x654 - x6 * x665 + x99 * (x646 - x647 + x652))
            - x6
            * (
                x3 * x665
                - x6 * (x3 * x663 - x6 * (x3 * x661 - x6 * x668) + x666)
                + x99 * (x653 - x657 + x662)
            )
            + x83 * (x654 - x656 + x664 + x669)
        )
    )
    result[0, 9] = numpy.sum(
        x141
        * (
            x3 * (x3 * x688 - x6 * x701 + x99 * (x679 - x680 + x686))
            - x6
            * (
                x3 * x701
                - x6 * (x3 * x699 - x6 * (x3 * x697 - x6 * x705) + x703)
                + x99 * (x687 - x692 + x698)
            )
            + x83 * (x688 - x691 + x700 + x706)
        )
    )
    result[1, 0] = numpy.sum(
        x781
        * (
            x3
            * (
                x3 * x750
                + x4 * (-x725 + x732 + x744)
                - x6 * x771
                + x83 * (x744 - x772 + x773)
            )
            - x6
            * (
                x3 * x771
                + x4 * (x746 - x752 + x774)
                - x6
                * (
                    x3 * x765
                    - x6
                    * (
                        x5 * x761
                        - x6 * (x5 * x759 - x6 * (x5 * x757 - x6 * x776) + x775)
                        + x99 * (x739 - x755 + x758)
                    )
                    + x83 * (x741 - x754 + x760 + x777)
                )
                + x83 * (x774 - x778 + x779)
            )
            + x83
            * (
                x3 * x773
                - x3 * x779
                - x6 * x779
                + x6 * (x764 + x778)
                + x749
                - x769
                - x99 * (-x3 * x721 + x723 + x768 + x780)
                + x99 * (x3 * x727 + x729 + x748 - x780)
            )
            + x99 * (x750 - x751 + x766 - x770)
        )
    )
    result[1, 1] = numpy.sum(
        x832
        * (
            x3
            * (
                x3 * x808
                + x4 * (-x790 + x796 + x804)
                - x6 * x823
                + x99 * (x804 - x824 + x825)
            )
            - x6
            * (
                x3 * x823
                + x4 * (x806 - x810 + x826)
                - x6
                * (
                    x3 * x819
                    - x6 * (x5 * x816 - x6 * (x5 * x814 - x6 * x828) + x827)
                    + x99 * (x801 - x812 + x815)
                )
                + x99 * (x826 - x829 + x830)
            )
            + x99 * (x808 - x809 + x820 - x822)
            + x99
            * (
                x3 * x825
                - x3 * x830
                - x4 * (-x3 * x786 + x788 + x800 + x831)
                + x4 * (x3 * x791 + x787 + x792 - x831)
                - x6 * x830
                + x6 * (x818 + x829)
                + x807
                - x821
            )
        )
    )
    result[1, 2] = numpy.sum(
        x832
        * (
            x3
            * (
                x3 * x860
                + x4 * (-x841 + x847 + x856)
                - x6 * x875
                + x99 * (x856 - x876 + x877)
            )
            - x6
            * (
                x3 * x875
                + x4 * (x858 - x862 + x879)
                - x6
                * (
                    x3 * x871
                    - x6 * (x5 * x868 - x6 * (x5 * x866 - x6 * x881) + x880)
                    + x99 * (x852 - x864 + x867)
                )
                + x99 * (x879 - x882 + x883)
            )
            + x99 * (x860 - x861 + x872 - x874)
            + x99
            * (
                x3 * x877
                - x3 * x883
                - x4 * (-x3 * x837 + x839 + x851 + x884)
                + x4 * (x3 * x842 + x838 + x843 - x884)
                - x6 * x883
                + x6 * (x870 + x882)
                + x859
                - x873
            )
        )
    )
    result[1, 3] = numpy.sum(
        x832
        * (
            x3
            * (
                x3 * x907
                + x4 * (x890 - x900 + x901)
                + x4 * (x901 - x904 + x905)
                - x6 * x912
            )
            + x4 * (x3 * x890 - x3 * x903 - x6 * x903 + x6 * (x893 + x897) + x887 + x899)
            - x6
            * (
                x3 * x912
                + x4 * (-x893 + x896 + x903)
                + x4 * (x896 + x906 - x909)
                - x6 * (x3 * x910 - x6 * (x5 * x895 - x6 * x914) + x913)
            )
            + x99 * (x899 + x907 - x908 + x911)
        )
    )
    result[1, 4] = numpy.sum(
        x945
        * (
            x3
            * (
                x3 * x937
                + x4 * (x920 - x930 + x931)
                + x4 * (x931 - x934 + x935)
                - x6 * x942
            )
            + x4 * (x3 * x920 - x3 * x933 - x6 * x933 + x6 * (x923 + x927) + x917 + x929)
            - x6
            * (
                x3 * x942
                + x4 * (-x923 + x926 + x933)
                + x4 * (x926 + x936 - x939)
                - x6 * (x3 * x940 - x6 * (x5 * x925 - x6 * x944) + x943)
            )
            + x99 * (x929 + x937 - x938 + x941)
        )
    )
    result[1, 5] = numpy.sum(
        x832
        * (
            x3
            * (
                x3 * x968
                + x4 * (x951 - x961 + x962)
                + x4 * (x962 - x965 + x966)
                - x6 * x973
            )
            + x4 * (x3 * x951 - x3 * x964 - x6 * x964 + x6 * (x954 + x958) + x948 + x960)
            - x6
            * (
                x3 * x973
                + x4 * (-x954 + x957 + x964)
                + x4 * (x957 + x967 - x970)
                - x6 * (x3 * x971 - x6 * (x5 * x956 - x6 * x975) + x974)
            )
            + x99 * (x960 + x968 - x969 + x972)
        )
    )
    result[1, 6] = numpy.sum(
        x781
        * (
            x3 * (x3 * x982 - x6 * x989 + x980)
            - x6 * (x3 * x989 - x6 * (x3 * x987 - x6 * x991) + x990)
            + x99 * (x982 - x983 + x988)
        )
    )
    result[1, 7] = numpy.sum(
        x832
        * (
            x3 * (-x1005 * x6 + x3 * x998 + x996)
            - x6 * (x1005 * x3 + x1006 - x6 * (x1003 * x3 - x1007 * x6))
            + x99 * (x1004 + x998 - x999)
        )
    )
    result[1, 8] = numpy.sum(
        x832
        * (
            x3 * (x1011 + x1013 * x3 - x1019 * x6)
            - x6 * (x1019 * x3 + x1020 - x6 * (x1017 * x3 - x1021 * x6))
            + x99 * (x1013 - x1014 + x1018)
        )
    )
    result[1, 9] = numpy.sum(
        x781
        * (
            x3 * (x1025 + x1027 * x3 - x1033 * x6)
            - x6 * (x1033 * x3 + x1034 - x6 * (x1031 * x3 - x1035 * x6))
            + x99 * (x1027 - x1028 + x1032)
        )
    )
    result[2, 0] = numpy.sum(
        x781
        * (
            x3
            * (
                x1079 * x3
                - x1100 * x6
                + x4 * (-x1054 + x1061 + x1073)
                + x83 * (x1073 - x1101 + x1102)
            )
            - x6
            * (
                x1100 * x3
                + x4 * (x1075 - x1081 + x1103)
                - x6
                * (
                    x1094 * x3
                    - x6
                    * (
                        x1090 * x5
                        - x6 * (x1088 * x5 + x1104 - x6 * (x1086 * x5 - x1105 * x6))
                        + x99 * (x1068 - x1084 + x1087)
                    )
                    + x83 * (x1070 - x1083 + x1089 + x1106)
                )
                + x83 * (x1103 - x1107 + x1108)
            )
            + x83
            * (
                x1078
                - x1098
                + x1102 * x3
                - x1108 * x3
                - x1108 * x6
                + x6 * (x1093 + x1107)
                - x99 * (-x1050 * x3 + x1052 + x1097 + x1109)
                + x99 * (x1056 * x3 + x1058 + x1077 - x1109)
            )
            + x99 * (x1079 - x1080 + x1095 - x1099)
        )
    )
    result[2, 1] = numpy.sum(
        x832
        * (
            x3
            * (
                x1144 * x3
                - x1161 * x6
                + x4 * (-x1122 + x1129 + x1140)
                + x99 * (x1140 - x1162 + x1163)
            )
            - x6
            * (
                x1161 * x3
                + x4 * (x1142 - x1146 + x1165)
                - x6
                * (
                    x1157 * x3
                    - x6 * (x1154 * x5 + x1166 - x6 * (x1152 * x5 - x1168 * x6))
                    + x99 * (x1136 - x1148 + x1153)
                )
                + x99 * (x1165 - x1169 + x1170)
            )
            + x99 * (x1144 - x1145 + x1158 - x1160)
            + x99
            * (
                x1143
                - x1159
                + x1163 * x3
                - x1170 * x3
                - x1170 * x6
                + x4 * (x1119 + x1124 * x3 + x1125 - x1171)
                - x4 * (-x1118 * x3 + x1120 + x1135 + x1171)
                + x6 * (x1156 + x1169)
            )
        )
    )
    result[2, 2] = numpy.sum(
        x832
        * (
            x3
            * (
                x1199 * x3
                - x1214 * x6
                + x4 * (-x1180 + x1186 + x1195)
                + x99 * (x1195 - x1215 + x1216)
            )
            - x6
            * (
                x1214 * x3
                + x4 * (x1197 - x1201 + x1218)
                - x6
                * (
                    x1210 * x3
                    - x6 * (x1207 * x5 + x1219 - x6 * (x1205 * x5 - x1220 * x6))
                    + x99 * (x1191 - x1203 + x1206)
                )
                + x99 * (x1218 - x1221 + x1222)
            )
            + x99 * (x1199 - x1200 + x1211 - x1213)
            + x99
            * (
                x1198
                - x1212
                + x1216 * x3
                - x1222 * x3
                - x1222 * x6
                + x4 * (x1177 + x1181 * x3 + x1182 - x1223)
                - x4 * (-x1176 * x3 + x1178 + x1190 + x1223)
                + x6 * (x1209 + x1221)
            )
        )
    )
    result[2, 3] = numpy.sum(
        x832
        * (
            x3
            * (
                x1255 * x3
                - x1260 * x6
                + x4 * (x1233 - x1248 + x1249)
                + x4 * (x1249 - x1252 + x1253)
            )
            + x4
            * (
                x1229
                + x1233 * x3
                + x1247
                - x1251 * x3
                - x1251 * x6
                + x6 * (x1238 + x1244)
            )
            - x6
            * (
                x1260 * x3
                + x4 * (-x1238 + x1243 + x1251)
                + x4 * (x1243 + x1254 - x1257)
                - x6 * (x1258 * x3 + x1262 - x6 * (x1242 * x5 - x1264 * x6))
            )
            + x99 * (x1247 + x1255 - x1256 + x1259)
        )
    )
    result[2, 4] = numpy.sum(
        x945
        * (
            x3
            * (
                x1294 * x3
                - x1299 * x6
                + x4 * (x1273 - x1287 + x1288)
                + x4 * (x1288 - x1291 + x1292)
            )
            + x4
            * (
                x1269
                + x1273 * x3
                + x1286
                - x1290 * x3
                - x1290 * x6
                + x6 * (x1278 + x1284)
            )
            - x6
            * (
                x1299 * x3
                + x4 * (-x1278 + x1283 + x1290)
                + x4 * (x1283 + x1293 - x1296)
                - x6 * (x1297 * x3 + x1300 - x6 * (x1282 * x5 - x1302 * x6))
            )
            + x99 * (x1286 + x1294 - x1295 + x1298)
        )
    )
    result[2, 5] = numpy.sum(
        x832
        * (
            x3
            * (
                x1325 * x3
                - x1330 * x6
                + x4 * (x1308 - x1318 + x1319)
                + x4 * (x1319 - x1322 + x1323)
            )
            + x4
            * (
                x1305
                + x1308 * x3
                + x1317
                - x1321 * x3
                - x1321 * x6
                + x6 * (x1311 + x1315)
            )
            - x6
            * (
                x1330 * x3
                + x4 * (-x1311 + x1314 + x1321)
                + x4 * (x1314 + x1324 - x1327)
                - x6 * (x1328 * x3 + x1331 - x6 * (x1313 * x5 - x1332 * x6))
            )
            + x99 * (x1317 + x1325 - x1326 + x1329)
        )
    )
    result[2, 6] = numpy.sum(
        x781
        * (
            x3 * (x1339 + x1342 * x3 - x1350 * x6)
            - x6 * (x1350 * x3 + x1352 - x6 * (x1348 * x3 - x1354 * x6))
            + x99 * (x1342 - x1343 + x1349)
        )
    )
    result[2, 7] = numpy.sum(
        x832
        * (
            x3 * (x1360 + x1363 * x3 - x1370 * x6)
            - x6 * (x1370 * x3 + x1372 - x6 * (x1368 * x3 - x1374 * x6))
            + x99 * (x1363 - x1364 + x1369)
        )
    )
    result[2, 8] = numpy.sum(
        x832
        * (
            x3 * (x1379 + x1382 * x3 - x1389 * x6)
            - x6 * (x1389 * x3 + x1390 - x6 * (x1387 * x3 - x1392 * x6))
            + x99 * (x1382 - x1383 + x1388)
        )
    )
    result[2, 9] = numpy.sum(
        x781
        * (
            x3 * (x1397 + x1399 * x3 - x1406 * x6)
            - x6 * (x1406 * x3 + x1407 - x6 * (x1404 * x3 - x1408 * x6))
            + x99 * (x1399 - x1400 + x1405)
        )
    )
    result[3, 0] = numpy.sum(
        x1453
        * (
            x3 * (x1434 * x3 - x1448 * x6 + x83 * (-x1419 + x1426 + x1432 + x1449))
            + x4 * (-x1428 + x1434 + x1447)
            - x6
            * (
                x1448 * x3
                - x6
                * (
                    x1444 * x5
                    - x6 * (x1442 * x5 + x1450 - x6 * (x1440 * x5 - x1451 * x6))
                    + x99 * (x1425 - x1437 + x1441)
                )
                + x83 * (x1427 - x1436 + x1443 + x1452)
            )
            + x83 * (-x1427 * x3 + x1432 * x3 + x1433 + x1447)
        )
    )
    result[3, 1] = numpy.sum(
        x1481
        * (
            x3 * (x1468 * x3 - x1478 * x6 + x99 * (-x1457 + x1461 + x1465))
            + x4 * (-x1463 + x1468 + x1477)
            - x6
            * (
                x1478 * x3
                - x6 * (x1474 * x5 + x1479 - x6 * (x1472 * x5 - x1480 * x6))
                + x99 * (x1462 - x1469 + x1473)
            )
            + x99 * (-x1462 * x3 + x1465 * x3 + x1467 + x1477)
        )
    )
    result[3, 2] = numpy.sum(
        x1481
        * (
            x3 * (x1496 * x3 - x1507 * x6 + x99 * (-x1485 + x1489 + x1493))
            + x4 * (-x1491 + x1496 + x1506)
            - x6
            * (
                x1507 * x3
                - x6 * (x1502 * x5 + x1508 - x6 * (x1500 * x5 - x1509 * x6))
                + x99 * (x1490 - x1497 + x1501)
            )
            + x99 * (-x1490 * x3 + x1493 * x3 + x1495 + x1506)
        )
    )
    result[3, 3] = numpy.sum(
        x1481
        * (
            x3 * (x1523 * x3 + x1525 - x1526 * x6)
            + x4 * (x1516 - x1522 + x1523)
            + x4 * (x1511 * x3 + x1516 - x1520 * x3 + x1521)
            - x6 * (x1526 * x3 + x1528 - x6 * (x1515 * x5 - x1530 * x6))
        )
    )
    result[3, 4] = numpy.sum(
        x1546
        * (
            x3 * (x1541 * x3 + x1542 - x1543 * x6)
            + x4 * (x1535 - x1540 + x1541)
            + x4 * (x1531 * x3 + x1535 - x1538 * x3 + x1539)
            - x6 * (x1543 * x3 + x1544 - x6 * (x1534 * x5 - x1545 * x6))
        )
    )
    result[3, 5] = numpy.sum(
        x1481
        * (
            x3 * (x1557 * x3 + x1558 - x1559 * x6)
            + x4 * (x1551 - x1556 + x1557)
            + x4 * (x1547 * x3 + x1551 - x1554 * x3 + x1555)
            - x6 * (x1559 * x3 + x1560 - x6 * (x1550 * x5 - x1561 * x6))
        )
    )
    result[3, 6] = numpy.sum(
        x1453 * (x1566 + x3 * (x1565 * x3 - x1567 * x6) - x6 * (x1567 * x3 - x1568 * x6))
    )
    result[3, 7] = numpy.sum(
        x1481 * (x1575 + x3 * (x1574 * x3 - x1576 * x6) - x6 * (x1576 * x3 - x1578 * x6))
    )
    result[3, 8] = numpy.sum(
        x1481 * (x1582 + x3 * (x1581 * x3 - x1583 * x6) - x6 * (x1583 * x3 - x1584 * x6))
    )
    result[3, 9] = numpy.sum(
        x1453 * (x1588 + x3 * (x1587 * x3 - x1589 * x6) - x6 * (x1589 * x3 - x1590 * x6))
    )
    result[4, 0] = numpy.sum(
        x832
        * (
            x3 * (x1613 * x3 - x1626 * x6 + x83 * (-x1599 + x1605 + x1611 + x1627))
            + x4 * (-x1607 + x1613 + x1625)
            - x6
            * (
                x1626 * x3
                - x6
                * (
                    x1622 * x5
                    - x6 * (x1620 * x5 + x1628 - x6 * (x1618 * x5 - x1629 * x6))
                    + x99 * (x1604 - x1616 + x1619)
                )
                + x83 * (x1606 - x1615 + x1621 + x1630)
            )
            + x83 * (-x1606 * x3 + x1611 * x3 + x1612 + x1625)
        )
    )
    result[4, 1] = numpy.sum(
        x1655
        * (
            x3 * (x1643 * x3 - x1652 * x6 + x99 * (-x1633 + x1636 + x1640))
            + x4 * (-x1638 + x1643 + x1651)
            - x6
            * (
                x1652 * x3
                - x6 * (x1648 * x5 + x1653 - x6 * (x1646 * x5 - x1654 * x6))
                + x99 * (x1637 - x1644 + x1647)
            )
            + x99 * (-x1637 * x3 + x1640 * x3 + x1642 + x1651)
        )
    )
    result[4, 2] = numpy.sum(
        x1655
        * (
            x3 * (x1668 * x3 - x1678 * x6 + x99 * (-x1658 + x1661 + x1665))
            + x4 * (-x1663 + x1668 + x1677)
            - x6
            * (
                x1678 * x3
                - x6 * (x1673 * x5 + x1679 - x6 * (x1671 * x5 - x1680 * x6))
                + x99 * (x1662 - x1669 + x1672)
            )
            + x99 * (-x1662 * x3 + x1665 * x3 + x1667 + x1677)
        )
    )
    result[4, 3] = numpy.sum(
        x1655
        * (
            x3 * (x1689 * x3 + x1690 - x1691 * x6)
            + x4 * (x1684 - x1688 + x1689)
            + x4 * (x1681 * x3 + x1684 - x1686 * x3 + x1687)
            - x6 * (x1691 * x3 + x1692 - x6 * (x1683 * x5 - x1693 * x6))
        )
    )
    result[4, 4] = numpy.sum(
        x1707
        * (
            x3 * (x1702 * x3 + x1703 - x1704 * x6)
            + x4 * (x1697 - x1701 + x1702)
            + x4 * (x1694 * x3 + x1697 - x1699 * x3 + x1700)
            - x6 * (x1704 * x3 + x1705 - x6 * (x1696 * x5 - x1706 * x6))
        )
    )
    result[4, 5] = numpy.sum(
        x1655
        * (
            x3 * (x1716 * x3 + x1717 - x1718 * x6)
            + x4 * (x1711 - x1715 + x1716)
            + x4 * (x1708 * x3 + x1711 - x1713 * x3 + x1714)
            - x6 * (x1718 * x3 + x1719 - x6 * (x1710 * x5 - x1720 * x6))
        )
    )
    result[4, 6] = numpy.sum(
        x832 * (x1725 + x3 * (x1724 * x3 - x1726 * x6) - x6 * (x1726 * x3 - x1727 * x6))
    )
    result[4, 7] = numpy.sum(
        x1655 * (x1732 + x3 * (x1731 * x3 - x1733 * x6) - x6 * (x1733 * x3 - x1734 * x6))
    )
    result[4, 8] = numpy.sum(
        x1655 * (x1738 + x3 * (x1737 * x3 - x1739 * x6) - x6 * (x1739 * x3 - x1740 * x6))
    )
    result[4, 9] = numpy.sum(
        x832 * (x1744 + x3 * (x1743 * x3 - x1745 * x6) - x6 * (x1745 * x3 - x1746 * x6))
    )
    result[5, 0] = numpy.sum(
        x1453
        * (
            x3 * (x1772 * x3 - x1786 * x6 + x83 * (-x1757 + x1764 + x1770 + x1787))
            + x4 * (-x1766 + x1772 + x1785)
            - x6
            * (
                x1786 * x3
                - x6
                * (
                    x1782 * x5
                    - x6 * (x1780 * x5 + x1788 - x6 * (x1778 * x5 - x1789 * x6))
                    + x99 * (x1763 - x1775 + x1779)
                )
                + x83 * (x1765 - x1774 + x1781 + x1790)
            )
            + x83 * (-x1765 * x3 + x1770 * x3 + x1771 + x1785)
        )
    )
    result[5, 1] = numpy.sum(
        x1481
        * (
            x3 * (x1808 * x3 - x1820 * x6 + x99 * (-x1795 + x1800 + x1805))
            + x4 * (-x1802 + x1808 + x1819)
            - x6
            * (
                x1820 * x3
                - x6 * (x1815 * x5 + x1821 - x6 * (x1813 * x5 - x1823 * x6))
                + x99 * (x1801 - x1809 + x1814)
            )
            + x99 * (-x1801 * x3 + x1805 * x3 + x1807 + x1819)
        )
    )
    result[5, 2] = numpy.sum(
        x1481
        * (
            x3 * (x1838 * x3 - x1849 * x6 + x99 * (-x1827 + x1831 + x1835))
            + x4 * (-x1833 + x1838 + x1848)
            - x6
            * (
                x1849 * x3
                - x6 * (x1844 * x5 + x1850 - x6 * (x1842 * x5 - x1851 * x6))
                + x99 * (x1832 - x1839 + x1843)
            )
            + x99 * (-x1832 * x3 + x1835 * x3 + x1837 + x1848)
        )
    )
    result[5, 3] = numpy.sum(
        x1481
        * (
            x3 * (x1865 * x3 + x1867 - x1868 * x6)
            + x4 * (x1858 - x1864 + x1865)
            + x4 * (x1853 * x3 + x1858 - x1862 * x3 + x1863)
            - x6 * (x1868 * x3 + x1870 - x6 * (x1857 * x5 - x1872 * x6))
        )
    )
    result[5, 4] = numpy.sum(
        x1546
        * (
            x3 * (x1886 * x3 + x1887 - x1888 * x6)
            + x4 * (x1879 - x1885 + x1886)
            + x4 * (x1874 * x3 + x1879 - x1883 * x3 + x1884)
            - x6 * (x1888 * x3 + x1889 - x6 * (x1878 * x5 - x1891 * x6))
        )
    )
    result[5, 5] = numpy.sum(
        x1481
        * (
            x3 * (x1905 * x3 + x1907 - x1908 * x6)
            + x4 * (x1898 - x1904 + x1905)
            + x4 * (x1893 * x3 + x1898 - x1902 * x3 + x1903)
            - x6 * (x1908 * x3 + x1910 - x6 * (x1897 * x5 - x1912 * x6))
        )
    )
    result[5, 6] = numpy.sum(
        x1453 * (x1919 + x3 * (x1915 * x3 - x1921 * x6) - x6 * (x1921 * x3 - x1923 * x6))
    )
    result[5, 7] = numpy.sum(
        x1481 * (x1929 + x3 * (x1926 * x3 - x1931 * x6) - x6 * (x1931 * x3 - x1933 * x6))
    )
    result[5, 8] = numpy.sum(
        x1481 * (x1938 + x3 * (x1937 * x3 - x1940 * x6) - x6 * (x1940 * x3 - x1942 * x6))
    )
    result[5, 9] = numpy.sum(
        x1453 * (x1947 + x3 * (x1946 * x3 - x1948 * x6) - x6 * (x1948 * x3 - x1949 * x6))
    )
    result[6, 0] = numpy.sum(
        x781
        * (
            x3 * (x1964 * x5 - x1975 * x6 + x99 * (x1956 - x1957 + x1962))
            - x6
            * (
                x1975 * x5
                - x6 * (x1973 * x5 + x1976 - x6 * (x1971 * x5 - x1977 * x6))
                + x99 * (x1963 - x1967 + x1972)
            )
            + x83 * (x1964 - x1966 + x1974 + x1978)
        )
    )
    result[6, 1] = numpy.sum(
        x832
        * (
            x3 * (x1982 + x1984 * x5 - x1990 * x6)
            - x6 * (x1990 * x5 + x1991 - x6 * (x1988 * x5 - x1992 * x6))
            + x99 * (x1984 - x1985 + x1989)
        )
    )
    result[6, 2] = numpy.sum(
        x832
        * (
            x3 * (x1997 + x1999 * x5 - x2006 * x6)
            - x6 * (x2006 * x5 + x2007 - x6 * (x2004 * x5 - x2008 * x6))
            + x99 * (x1999 - x2000 + x2005)
        )
    )
    result[6, 3] = numpy.sum(
        x832 * (x2014 + x3 * (x2013 * x5 - x2015 * x6) - x6 * (x2015 * x5 - x2016 * x6))
    )
    result[6, 4] = numpy.sum(
        x945 * (x2020 + x3 * (x2019 * x5 - x2021 * x6) - x6 * (x2021 * x5 - x2022 * x6))
    )
    result[6, 5] = numpy.sum(
        x832 * (x2027 + x3 * (x2026 * x5 - x2028 * x6) - x6 * (x2028 * x5 - x2029 * x6))
    )
    result[6, 6] = numpy.sum(x781 * (x2031 * x3 - x2033 * x6))
    result[6, 7] = numpy.sum(x832 * (x2034 * x3 - x2035 * x6))
    result[6, 8] = numpy.sum(x832 * (x2036 * x3 - x2037 * x6))
    result[6, 9] = numpy.sum(x781 * (x2038 * x3 - x2039 * x6))
    result[7, 0] = numpy.sum(
        x832
        * (
            x3 * (x2052 * x5 - x2062 * x6 + x99 * (x2045 - x2046 + x2050))
            - x6
            * (
                x2062 * x5
                - x6 * (x2060 * x5 + x2063 - x6 * (x2058 * x5 - x2064 * x6))
                + x99 * (x2051 - x2055 + x2059)
            )
            + x83 * (x2052 - x2054 + x2061 + x2065)
        )
    )
    result[7, 1] = numpy.sum(
        x1655
        * (
            x3 * (x2069 + x2071 * x5 - x2077 * x6)
            - x6 * (x2077 * x5 + x2078 - x6 * (x2075 * x5 - x2079 * x6))
            + x99 * (x2071 - x2072 + x2076)
        )
    )
    result[7, 2] = numpy.sum(
        x1655
        * (
            x3 * (x2083 + x2085 * x5 - x2091 * x6)
            - x6 * (x2091 * x5 + x2092 - x6 * (x2089 * x5 - x2093 * x6))
            + x99 * (x2085 - x2086 + x2090)
        )
    )
    result[7, 3] = numpy.sum(
        x1655 * (x2100 + x3 * (x2099 * x5 - x2101 * x6) - x6 * (x2101 * x5 - x2103 * x6))
    )
    result[7, 4] = numpy.sum(
        x1707 * (x2107 + x3 * (x2106 * x5 - x2108 * x6) - x6 * (x2108 * x5 - x2109 * x6))
    )
    result[7, 5] = numpy.sum(
        x1655 * (x2113 + x3 * (x2112 * x5 - x2114 * x6) - x6 * (x2114 * x5 - x2115 * x6))
    )
    result[7, 6] = numpy.sum(x832 * (x2116 * x3 - x2117 * x6))
    result[7, 7] = numpy.sum(x1655 * (x2119 * x3 - x2121 * x6))
    result[7, 8] = numpy.sum(x1655 * (x2122 * x3 - x2123 * x6))
    result[7, 9] = numpy.sum(x832 * (x2124 * x3 - x2125 * x6))
    result[8, 0] = numpy.sum(
        x832
        * (
            x3 * (x2136 * x5 - x2145 * x6 + x99 * (x2130 - x2131 + x2134))
            - x6
            * (
                x2145 * x5
                - x6 * (x2143 * x5 + x2146 - x6 * (x2141 * x5 - x2147 * x6))
                + x99 * (x2135 - x2139 + x2142)
            )
            + x83 * (x2136 - x2138 + x2144 + x2148)
        )
    )
    result[8, 1] = numpy.sum(
        x1655
        * (
            x3 * (x2151 + x2153 * x5 - x2158 * x6)
            - x6 * (x2158 * x5 + x2159 - x6 * (x2156 * x5 - x2160 * x6))
            + x99 * (x2153 - x2154 + x2157)
        )
    )
    result[8, 2] = numpy.sum(
        x1655
        * (
            x3 * (x2163 + x2165 * x5 - x2170 * x6)
            - x6 * (x2170 * x5 + x2171 - x6 * (x2168 * x5 - x2172 * x6))
            + x99 * (x2165 - x2166 + x2169)
        )
    )
    result[8, 3] = numpy.sum(
        x1655 * (x2175 + x3 * (x2174 * x5 - x2176 * x6) - x6 * (x2176 * x5 - x2177 * x6))
    )
    result[8, 4] = numpy.sum(
        x1707 * (x2180 + x3 * (x2179 * x5 - x2181 * x6) - x6 * (x2181 * x5 - x2182 * x6))
    )
    result[8, 5] = numpy.sum(
        x1655 * (x2185 + x3 * (x2184 * x5 - x2186 * x6) - x6 * (x2186 * x5 - x2187 * x6))
    )
    result[8, 6] = numpy.sum(x832 * (x2188 * x3 - x2189 * x6))
    result[8, 7] = numpy.sum(x1655 * (x2190 * x3 - x2191 * x6))
    result[8, 8] = numpy.sum(x1655 * (x2192 * x3 - x2193 * x6))
    result[8, 9] = numpy.sum(x832 * (x2194 * x3 - x2195 * x6))
    result[9, 0] = numpy.sum(
        x781
        * (
            x3 * (x2210 * x5 - x2221 * x6 + x99 * (x2202 - x2203 + x2208))
            - x6
            * (
                x2221 * x5
                - x6 * (x2219 * x5 + x2222 - x6 * (x2217 * x5 - x2223 * x6))
                + x99 * (x2209 - x2213 + x2218)
            )
            + x83 * (x2210 - x2212 + x2220 + x2224)
        )
    )
    result[9, 1] = numpy.sum(
        x832
        * (
            x3 * (x2229 + x2232 * x5 - x2239 * x6)
            - x6 * (x2239 * x5 + x2240 - x6 * (x2237 * x5 - x2242 * x6))
            + x99 * (x2232 - x2233 + x2238)
        )
    )
    result[9, 2] = numpy.sum(
        x832
        * (
            x3 * (x2246 + x2248 * x5 - x2254 * x6)
            - x6 * (x2254 * x5 + x2255 - x6 * (x2252 * x5 - x2256 * x6))
            + x99 * (x2248 - x2249 + x2253)
        )
    )
    result[9, 3] = numpy.sum(
        x832 * (x2262 + x3 * (x2259 * x5 - x2264 * x6) - x6 * (x2264 * x5 - x2266 * x6))
    )
    result[9, 4] = numpy.sum(
        x945 * (x2271 + x3 * (x2270 * x5 - x2273 * x6) - x6 * (x2273 * x5 - x2275 * x6))
    )
    result[9, 5] = numpy.sum(
        x832 * (x2281 + x3 * (x2280 * x5 - x2282 * x6) - x6 * (x2282 * x5 - x2283 * x6))
    )
    result[9, 6] = numpy.sum(x781 * (x2285 * x3 - x2287 * x6))
    result[9, 7] = numpy.sum(x832 * (x2289 * x3 - x2291 * x6))
    result[9, 8] = numpy.sum(x832 * (x2293 * x3 - x2295 * x6))
    result[9, 9] = numpy.sum(x781 * (x2297 * x3 - x2299 * x6))
    result[10, 0] = numpy.sum(
        x141
        * (
            x5 * (x2305 * x5 - x2312 * x6 + x4 * (x2300 - x2301 - x2302 + x2303))
            - x6
            * (
                x2312 * x5
                + x4 * (x2304 - x2307 + x2308 - x2309)
                - x6 * (x2310 * x5 - x6 * (3.0 * x1450 + x1971 * x707 - x1977 * x8))
            )
            + x99 * (x2305 - x2306 + x2311)
        )
    )
    result[10, 1] = numpy.sum(
        x251
        * (
            x4 * (x1978 + x2313 - x2314 - x2315 + x2316)
            + x5 * (x2316 * x5 - x2317 * x6)
            - x6 * (x2317 * x5 - x6 * (3.0 * x1479 + x1976 + x1988 * x707 - x1992 * x8))
        )
    )
    result[10, 2] = numpy.sum(
        x251
        * (
            x4 * (x2318 - x2319 - x2320 + x2321)
            + x5 * (x2321 * x5 - x2322 * x6)
            - x6 * (x2322 * x5 - x6 * (3.0 * x1508 + x2004 * x707 - x2008 * x8))
        )
    )
    result[10, 3] = numpy.sum(
        x251
        * (
            x5 * (2.0 * x1982 + x2013 * x707 - x2015 * x8 + x2030)
            - x6 * (2.0 * x1991 + x2015 * x707 - x2016 * x8 + x2032)
        )
    )
    result[10, 4] = numpy.sum(
        x497
        * (
            x5 * (3.0 * x1542 + x1997 + x2019 * x707 - x2021 * x8)
            - x6 * (3.0 * x1544 + x2007 + x2021 * x707 - x2022 * x8)
        )
    )
    result[10, 5] = numpy.sum(
        x251
        * (
            x5 * (3.0 * x1558 + x2026 * x707 - x2028 * x8)
            - x6 * (3.0 * x1560 + x2028 * x707 - x2029 * x8)
        )
    )
    result[10, 6] = numpy.sum(
        x141 * (3.0 * x1566 + 3.0 * x2014 + x2031 * x707 - x2033 * x8)
    )
    result[10, 7] = numpy.sum(
        x251 * (3.0 * x1575 + 2.0 * x2020 + x2034 * x707 - x2035 * x8)
    )
    result[10, 8] = numpy.sum(x251 * (3.0 * x1582 + x2027 + x2036 * x707 - x2037 * x8))
    result[10, 9] = numpy.sum(x141 * (3.0 * x1588 + x2038 * x707 - x2039 * x8))
    result[11, 0] = numpy.sum(
        x781
        * (
            x5 * (x2328 * x5 - x2335 * x6 + x4 * (x2323 - x2324 - x2325 + x2326))
            - x6
            * (
                x2335 * x5
                + x4 * (x2327 - x2330 + x2331 - x2332)
                - x6 * (x2333 * x5 - x6 * (2.0 * x1628 + x2058 * x707 - x2064 * x8))
            )
            + x99 * (x2328 - x2329 + x2334)
        )
    )
    result[11, 1] = numpy.sum(
        x832
        * (
            x4 * (x2065 + x2097 + x2336 - x2337 + x2338)
            + x5 * (x2338 * x5 - x2339 * x6)
            - x6 * (x2339 * x5 - x6 * (x2063 + x2075 * x707 - x2079 * x8 + x2102))
        )
    )
    result[11, 2] = numpy.sum(
        x832
        * (
            x4 * (x2340 - x2341 - x2342 + x2343)
            + x5 * (x2343 * x5 - x2344 * x6)
            - x6 * (x2344 * x5 - x6 * (2.0 * x1679 + x2089 * x707 - x2093 * x8))
        )
    )
    result[11, 3] = numpy.sum(
        x832
        * (
            x5 * (2.0 * x1690 + 2.0 * x2069 + x2099 * x707 - x2101 * x8)
            - x6 * (2.0 * x1692 + 2.0 * x2078 + x2101 * x707 - x2103 * x8)
        )
    )
    result[11, 4] = numpy.sum(
        x945
        * (
            x5 * (x2083 + x2106 * x707 - x2108 * x8 + x2118)
            - x6 * (x2092 + x2108 * x707 - x2109 * x8 + x2120)
        )
    )
    result[11, 5] = numpy.sum(
        x832
        * (
            x5 * (2.0 * x1717 + x2112 * x707 - x2114 * x8)
            - x6 * (2.0 * x1719 + x2114 * x707 - x2115 * x8)
        )
    )
    result[11, 6] = numpy.sum(
        x781 * (2.0 * x1725 + 3.0 * x2100 + x2116 * x707 - x2117 * x8)
    )
    result[11, 7] = numpy.sum(
        x832 * (2.0 * x1732 + 2.0 * x2107 + x2119 * x707 - x2121 * x8)
    )
    result[11, 8] = numpy.sum(x832 * (2.0 * x1738 + x2113 + x2122 * x707 - x2123 * x8))
    result[11, 9] = numpy.sum(x781 * (2.0 * x1744 + x2124 * x707 - x2125 * x8))
    result[12, 0] = numpy.sum(
        x1453
        * (
            x5 * (x2349 * x5 - x2355 * x6 + x4 * (x1787 + x2345 - x2346 + x2347))
            - x6
            * (
                x2355 * x5
                + x4 * (x1790 + x2348 - x2351 + x2352)
                - x6 * (x2353 * x5 - x6 * (x1788 + x2141 * x707 - x2147 * x8))
            )
            + x99 * (x2349 - x2350 + x2354)
        )
    )
    result[12, 1] = numpy.sum(
        x1481
        * (
            x4 * (x1818 + x2148 + x2356 - x2357 + x2358)
            + x5 * (x2358 * x5 - x2359 * x6)
            - x6 * (x2359 * x5 - x6 * (x1821 + x2146 + x2156 * x707 - x2160 * x8))
        )
    )
    result[12, 2] = numpy.sum(
        x1481
        * (
            x4 * (x1847 + x2360 - x2361 + x2362)
            + x5 * (x2362 * x5 - x2363 * x6)
            - x6 * (x2363 * x5 - x6 * (x1850 + x2168 * x707 - x2172 * x8))
        )
    )
    result[12, 3] = numpy.sum(
        x1481
        * (
            x5 * (x1867 + 2.0 * x2151 + x2174 * x707 - x2176 * x8)
            - x6 * (x1870 + 2.0 * x2159 + x2176 * x707 - x2177 * x8)
        )
    )
    result[12, 4] = numpy.sum(
        x1546
        * (
            x5 * (x1887 + x2163 + x2179 * x707 - x2181 * x8)
            - x6 * (x1889 + x2171 + x2181 * x707 - x2182 * x8)
        )
    )
    result[12, 5] = numpy.sum(
        x1481
        * (
            x5 * (x1907 + x2184 * x707 - x2186 * x8)
            - x6 * (x1910 + x2186 * x707 - x2187 * x8)
        )
    )
    result[12, 6] = numpy.sum(x1453 * (x1919 + 3.0 * x2175 + x2188 * x707 - x2189 * x8))
    result[12, 7] = numpy.sum(x1481 * (x1929 + 2.0 * x2180 + x2190 * x707 - x2191 * x8))
    result[12, 8] = numpy.sum(x1481 * (x1938 + x2185 + x2192 * x707 - x2193 * x8))
    result[12, 9] = numpy.sum(x1453 * (x1947 + x2194 * x707 - x2195 * x8))
    result[13, 0] = numpy.sum(
        x781
        * (
            x5 * (x2367 * x5 - x2372 * x6 + x4 * (x2225 - x2364 + x2365))
            - x6
            * (
                x2372 * x5
                + x4 * (x2235 + x2366 - x2369)
                - x6 * (x2370 * x5 - x6 * (x2217 * x707 + x2241))
            )
            + x99 * (x2367 - x2368 + x2371)
        )
    )
    result[13, 1] = numpy.sum(
        x832
        * (
            x4 * (x2261 - x2373 + x2374)
            + x5 * (x2374 * x5 - x2375 * x6)
            - x6 * (x2375 * x5 - x6 * (x2237 * x707 + x2265))
        )
    )
    result[13, 2] = numpy.sum(
        x832
        * (
            x4 * (x2267 - x2376 + x2377)
            + x5 * (x2377 * x5 - x2378 * x6)
            - x6 * (x2378 * x5 - x6 * (x2252 * x707 + x2274))
        )
    )
    result[13, 3] = numpy.sum(
        x832 * (x5 * (x2259 * x707 + x2284) - x6 * (x2264 * x707 + x2286))
    )
    result[13, 4] = numpy.sum(
        x945 * (x5 * (x2270 * x707 + x2288) - x6 * (x2273 * x707 + x2290))
    )
    result[13, 5] = numpy.sum(
        x832 * (x5 * (x2280 * x707 + x2292) - x6 * (x2282 * x707 + x2294))
    )
    result[13, 6] = numpy.sum(x781 * (3.0 * x2262 + x2285 * x707 - x2287 * x8))
    result[13, 7] = numpy.sum(x832 * (2.0 * x2271 + x2289 * x707 - x2291 * x8))
    result[13, 8] = numpy.sum(x832 * (x2281 + x2293 * x707 - x2295 * x8))
    result[13, 9] = numpy.sum(x781 * (x2297 * x707 - x2299 * x8))
    result[14, 0] = numpy.sum(
        x141
        * (
            x5 * (x2383 + x2385 * x5 - x2392 * x6)
            - x6 * (x2392 * x5 + x2393 - x6 * (x2390 * x5 - x2394 * x6))
            + x99 * (x2385 - x2386 + x2391)
        )
    )
    result[14, 1] = numpy.sum(
        x251 * (x2398 + x5 * (x2397 * x5 - x2399 * x6) - x6 * (x2399 * x5 - x2400 * x6))
    )
    result[14, 2] = numpy.sum(
        x251 * (x2405 + x5 * (x2404 * x5 - x2406 * x6) - x6 * (x2406 * x5 - x2407 * x6))
    )
    result[14, 3] = numpy.sum(x251 * (x2408 * x5 - x2409 * x6))
    result[14, 4] = numpy.sum(x497 * (x2410 * x5 - x2411 * x6))
    result[14, 5] = numpy.sum(x251 * (x2412 * x5 - x2413 * x6))
    result[14, 6] = numpy.sum(x141 * (x144 * x2408 + 2.0 * x2398 - x2409 * x8))
    result[14, 7] = numpy.sum(x251 * (x144 * x2410 + x2405 - x2411 * x8))
    result[14, 8] = numpy.sum(x251 * (x144 * x2412 - x2413 * x8))
    result[14, 9] = numpy.sum(
        x141 * (-x10 * x2299 + x1036 * x2297 + 3.0 * x1947 + 3.0 * x2281)
    )
    return result


def coulomb3d_44(ax, da, A, bx, db, B, C):
    """Cartesian (gg) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 15), dtype=float)

    x0 = ax + bx
    x1 = x0 ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = 0.5 / (ax + bx)
    x5 = -x2 - B[0]
    x6 = -x2 - C[0]
    x7 = -x1 * (ax * A[1] + bx * B[1])
    x8 = -x7 - C[1]
    x9 = -x1 * (ax * A[2] + bx * B[2])
    x10 = -x9 - C[2]
    x11 = x0 * (x10**2 + x6**2 + x8**2)
    x12 = (
        6.28318530717959
        * x1
        * numpy.exp(
            -ax * bx * x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2)
        )
    )
    x13 = x12 * boys(2, x11)
    x14 = x13 * x6
    x15 = x12 * boys(1, x11)
    x16 = x15 * x5
    x17 = -x14 + x16
    x18 = x17 * x5
    x19 = x4 * (-x13 + x15)
    x20 = x12 * boys(3, x11)
    x21 = x20 * x6
    x22 = x13 * x5
    x23 = -x21 + x22
    x24 = x23 * x6
    x25 = x19 - x24
    x26 = x18 + x25
    x27 = x26 * x5
    x28 = x23 * x5
    x29 = x4 * (x13 - x20)
    x30 = x12 * boys(4, x11)
    x31 = x30 * x6
    x32 = x20 * x5
    x33 = -x31 + x32
    x34 = x33 * x6
    x35 = x29 - x34
    x36 = x28 + x35
    x37 = x36 * x6
    x38 = 2.0 * x4
    x39 = x38 * (x17 + x21 - x22)
    x40 = -x37 + x39
    x41 = x27 + x40
    x42 = x41 * x5
    x43 = x12 * boys(0, x11)
    x44 = -x15 * x6 + x43 * x5
    x45 = x4 * (-x15 + x43)
    x46 = -x17 * x6 + x45
    x47 = x44 * x5 + x46
    x48 = -x26 * x6 + x38 * (x14 - x16 + x44)
    x49 = x47 * x5 + x48
    x50 = -x19
    x51 = x24 + x50
    x52 = x4 * (-x18 + x47 + x51)
    x53 = -x41 * x6 + 3.0 * x52
    x54 = x49 * x5 + x53
    x55 = x36 * x5
    x56 = x4 * (x20 - x30)
    x57 = x33 * x5
    x58 = x12 * boys(5, x11)
    x59 = x58 * x6
    x60 = x30 * x5
    x61 = -x59 + x60
    x62 = x6 * x61
    x63 = x56 + x57 - x62
    x64 = x6 * x63
    x65 = x38 * (x23 + x31 - x32)
    x66 = -x64 + x65
    x67 = x55 + x66
    x68 = x6 * x67
    x69 = -x29
    x70 = x34 + x69
    x71 = x4 * (x26 - x28 + x70)
    x72 = 3.0 * x71
    x73 = x68 - x72
    x74 = -x68 + x72
    x75 = x42 + x74
    x76 = x37 - x39
    x77 = x4 * (-x27 + x49 + x76)
    x78 = x3 * x54 - x6 * x75 + 4.0 * x77
    x79 = x3 * x75
    x80 = x5 * x67
    x81 = x5 * x63
    x82 = x4 * (x30 - x58)
    x83 = x5 * x61
    x84 = x12 * boys(6, x11)
    x85 = x6 * x84
    x86 = x5 * x58
    x87 = -x85 + x86
    x88 = x6 * x87
    x89 = x82 + x83 - x88
    x90 = x6 * x89
    x91 = x38 * (x33 + x59 - x60)
    x92 = -x90 + x91
    x93 = x81 + x92
    x94 = x6 * x93
    x95 = -x56
    x96 = x62 + x95
    x97 = x4 * (x36 - x57 + x96)
    x98 = 3.0 * x97
    x99 = -x94 + x98
    x100 = x80 + x99
    x101 = x100 * x6
    x102 = x64 - x65
    x103 = x4 * (x102 + x41 - x55)
    x104 = 4.0 * x103
    x105 = -x101 + x104 + x79
    x106 = x3 * x41
    x107 = x3 * x49 + x53
    x108 = -x106 + x107 + x73
    x109 = 4.0 * x4
    x110 = -x105 * x6 + x108 * x109 + x3 * x78 + x4 * (-x42 + x54 + x73)
    x111 = x94 - x98
    x112 = x4 * (x111 + x75 - x80)
    x113 = x105 * x3
    x114 = x100 * x3
    x115 = x5 * x93
    x116 = x5 * x89
    x117 = x4 * (x58 - x84)
    x118 = x5 * x87
    x119 = x12 * boys(7, x11)
    x120 = x119 * x6
    x121 = x5 * x84
    x122 = -x120 + x121
    x123 = x122 * x6
    x124 = x117 + x118 - x123
    x125 = x124 * x6
    x126 = x38 * (x61 + x85 - x86)
    x127 = x116 - x125 + x126
    x128 = x127 * x6
    x129 = -x82
    x130 = 3.0 * x4
    x131 = x130 * (x129 + x63 - x83 + x88)
    x132 = -x128 + x131
    x133 = x115 + x132
    x134 = x133 * x6
    x135 = x90 - x91
    x136 = x4 * (x135 + x67 - x81)
    x137 = 4.0 * x136
    x138 = x114 - x134 + x137
    x139 = x138 * x6
    x140 = x3 * x67
    x141 = x106 + x74
    x142 = x111 - x140 + x141
    x143 = x109 * x142
    x144 = x112 + x113 - x139 + x143
    x145 = x140 + x99
    x146 = x145 * x6
    x147 = x141 * x3
    x148 = x3 * x36
    x149 = x26 * x3
    x150 = x149 + x40
    x151 = x130 * (x102 - x148 + x150)
    x152 = x3 * x47 + x48
    x153 = x107 * x3 + x130 * (-x149 + x152 + x76) - x141 * x6 + x77
    x154 = x128 - x131
    x155 = x4 * (-x119 + x84)
    x156 = x12 * boys(8, x11)
    x157 = -x117
    x158 = x3 * x93
    x159 = x6 * (x132 + x158)
    x160 = x145 * x3
    x161 = x3 * x63
    x162 = x148 + x66
    x163 = x130 * (x135 - x161 + x162)
    x164 = x103 - x146 + x147 + x151
    x165 = x150 * x3
    x166 = x162 * x6
    x167 = x23 * x3
    x168 = x17 * x3
    x169 = x38 * (-x167 + x168 + x25 + x70)
    x170 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**5.5)
    x171 = 6.89597470414309 * x170
    x172 = x13 * x8
    x173 = -x172
    x174 = -x7 - B[1]
    x175 = x15 * x174
    x176 = x173 + x175
    x177 = x176 * x5
    x178 = x20 * x8
    x179 = -x178
    x180 = x13 * x174
    x181 = x179 + x180
    x182 = x181 * x6
    x183 = -x182
    x184 = x177 + x183
    x185 = x184 * x5
    x186 = x4 * (x176 + x178 - x180)
    x187 = x181 * x5
    x188 = x30 * x8
    x189 = -x188
    x190 = x174 * x20
    x191 = x189 + x190
    x192 = x191 * x6
    x193 = -x192
    x194 = x187 + x193
    x195 = x194 * x6
    x196 = x186 - x195
    x197 = x185 + x196
    x198 = x197 * x5
    x199 = -x15 * x8
    x200 = x174 * x43 + x199
    x201 = -x176 * x6
    x202 = x200 * x5 + x201
    x203 = x4 * (x172 - x175 + x200)
    x204 = -x184 * x6 + x203
    x205 = x202 * x5 + x204
    x206 = x4 * (-x177 + x182 + x202)
    x207 = -x197 * x6 + 2.0 * x206
    x208 = x205 * x5 + x207
    x209 = x194 * x5
    x210 = x4 * (x181 + x188 - x190)
    x211 = x191 * x5
    x212 = x58 * x8
    x213 = -x212
    x214 = x174 * x30
    x215 = x213 + x214
    x216 = x215 * x6
    x217 = x211 - x216
    x218 = x217 * x6
    x219 = x210 - x218
    x220 = x209 + x219
    x221 = x220 * x6
    x222 = x4 * (x184 - x187 + x192)
    x223 = 2.0 * x222
    x224 = x221 - x223
    x225 = -x221 + x223
    x226 = x198 + x225
    x227 = -x186
    x228 = x195 + x227
    x229 = x4 * (-x185 + x205 + x228)
    x230 = x208 * x3 - x226 * x6 + 3.0 * x229
    x231 = x226 * x3
    x232 = x220 * x5
    x233 = x217 * x5
    x234 = x4 * (x191 + x212 - x214)
    x235 = x215 * x5
    x236 = x8 * x84
    x237 = -x236
    x238 = x174 * x58
    x239 = x237 + x238
    x240 = x239 * x6
    x241 = x235 - x240
    x242 = x241 * x6
    x243 = x234 - x242
    x244 = x233 + x243
    x245 = x244 * x6
    x246 = x4 * (x194 - x211 + x216)
    x247 = 2.0 * x246
    x248 = -x245 + x247
    x249 = x232 + x248
    x250 = x249 * x6
    x251 = -x210
    x252 = x218 + x251
    x253 = x4 * (x197 - x209 + x252)
    x254 = 3.0 * x253
    x255 = x231 - x250 + x254
    x256 = x197 * x3
    x257 = x205 * x3 + x207
    x258 = x224 - x256 + x257
    x259 = x130 * x258 + x230 * x3 - x255 * x6 + x4 * (-x198 + x208 + x224)
    x260 = x245 - x247
    x261 = x4 * (x226 - x232 + x260)
    x262 = x255 * x3
    x263 = x249 * x3
    x264 = x244 * x5
    x265 = x4 * (x215 + x236 - x238)
    x266 = x241 * x5
    x267 = x239 * x5
    x268 = x119 * x8
    x269 = -x268
    x270 = x174 * x84
    x271 = x269 + x270
    x272 = x271 * x6
    x273 = x267 - x272
    x274 = x273 * x6
    x275 = x265 + x266 - x274
    x276 = x275 * x6
    x277 = x38 * (x217 - x235 + x240)
    x278 = -x276 + x277
    x279 = x264 + x278
    x280 = x279 * x6
    x281 = -x234
    x282 = x242 + x281
    x283 = x4 * (x220 - x233 + x282)
    x284 = 3.0 * x283
    x285 = x263 - x280 + x284
    x286 = x285 * x6
    x287 = x220 * x3
    x288 = x225 + x256
    x289 = x260 - x287 + x288
    x290 = x130 * x289
    x291 = x261 + x262 - x286 + x290
    x292 = x248 + x287
    x293 = x292 * x6
    x294 = x288 * x3
    x295 = x194 * x3
    x296 = x184 * x3
    x297 = x196 + x296
    x298 = x38 * (x252 - x295 + x297)
    x299 = x202 * x3 + x204
    x300 = x229 + x257 * x3 - x288 * x6 + x38 * (x228 - x296 + x299)
    x301 = x276 - x277
    x302 = x4 * (x239 + x268 - x270)
    x303 = -x156 * x8
    x304 = x119 * x174 + x303
    x305 = -x265
    x306 = x244 * x3
    x307 = x6 * (x278 + x306)
    x308 = x292 * x3
    x309 = x217 * x3
    x310 = x219 + x295
    x311 = x38 * (x282 - x309 + x310)
    x312 = x253 - x293 + x294 + x298
    x313 = x176 * x3
    x314 = x181 * x3
    x315 = x4 * (x183 + x192 + x313 - x314)
    x316 = x297 * x3
    x317 = x310 * x6
    x318 = 18.2450341145548 * x170
    x319 = x10 * x13
    x320 = -x319
    x321 = -x9 - B[2]
    x322 = x15 * x321
    x323 = x320 + x322
    x324 = x323 * x5
    x325 = x10 * x20
    x326 = -x325
    x327 = x13 * x321
    x328 = x326 + x327
    x329 = x328 * x6
    x330 = -x329
    x331 = x324 + x330
    x332 = x331 * x5
    x333 = x4 * (x323 + x325 - x327)
    x334 = x328 * x5
    x335 = x10 * x30
    x336 = -x335
    x337 = x20 * x321
    x338 = x336 + x337
    x339 = x338 * x6
    x340 = -x339
    x341 = x334 + x340
    x342 = x341 * x6
    x343 = x333 - x342
    x344 = x332 + x343
    x345 = x344 * x5
    x346 = -x10 * x15
    x347 = x321 * x43 + x346
    x348 = -x323 * x6
    x349 = x347 * x5 + x348
    x350 = x4 * (x319 - x322 + x347)
    x351 = -x331 * x6 + x350
    x352 = x349 * x5 + x351
    x353 = x4 * (-x324 + x329 + x349)
    x354 = -x344 * x6 + 2.0 * x353
    x355 = x352 * x5 + x354
    x356 = x341 * x5
    x357 = x4 * (x328 + x335 - x337)
    x358 = x338 * x5
    x359 = x10 * x58
    x360 = -x359
    x361 = x30 * x321
    x362 = x360 + x361
    x363 = x362 * x6
    x364 = x358 - x363
    x365 = x364 * x6
    x366 = x357 - x365
    x367 = x356 + x366
    x368 = x367 * x6
    x369 = x4 * (x331 - x334 + x339)
    x370 = 2.0 * x369
    x371 = x368 - x370
    x372 = -x368 + x370
    x373 = x345 + x372
    x374 = -x333
    x375 = x342 + x374
    x376 = x4 * (-x332 + x352 + x375)
    x377 = x3 * x355 - x373 * x6 + 3.0 * x376
    x378 = x3 * x373
    x379 = x367 * x5
    x380 = x364 * x5
    x381 = x4 * (x338 + x359 - x361)
    x382 = x362 * x5
    x383 = x10 * x84
    x384 = -x383
    x385 = x321 * x58
    x386 = x384 + x385
    x387 = x386 * x6
    x388 = x382 - x387
    x389 = x388 * x6
    x390 = x381 - x389
    x391 = x380 + x390
    x392 = x391 * x6
    x393 = x4 * (x341 - x358 + x363)
    x394 = 2.0 * x393
    x395 = -x392 + x394
    x396 = x379 + x395
    x397 = x396 * x6
    x398 = -x357
    x399 = x365 + x398
    x400 = x4 * (x344 - x356 + x399)
    x401 = 3.0 * x400
    x402 = x378 - x397 + x401
    x403 = x3 * x344
    x404 = x3 * x352 + x354
    x405 = x371 - x403 + x404
    x406 = x130 * x405 + x3 * x377 + x4 * (-x345 + x355 + x371) - x402 * x6
    x407 = x392 - x394
    x408 = x4 * (x373 - x379 + x407)
    x409 = x3 * x402
    x410 = x3 * x396
    x411 = x391 * x5
    x412 = x4 * (x362 + x383 - x385)
    x413 = x388 * x5
    x414 = x386 * x5
    x415 = x10 * x119
    x416 = -x415
    x417 = x321 * x84
    x418 = x416 + x417
    x419 = x418 * x6
    x420 = x414 - x419
    x421 = x420 * x6
    x422 = x412 + x413 - x421
    x423 = x422 * x6
    x424 = x38 * (x364 - x382 + x387)
    x425 = -x423 + x424
    x426 = x411 + x425
    x427 = x426 * x6
    x428 = -x381
    x429 = x389 + x428
    x430 = x4 * (x367 - x380 + x429)
    x431 = 3.0 * x430
    x432 = x410 - x427 + x431
    x433 = x432 * x6
    x434 = x3 * x367
    x435 = x372 + x403
    x436 = x407 - x434 + x435
    x437 = x130 * x436
    x438 = x408 + x409 - x433 + x437
    x439 = x395 + x434
    x440 = x439 * x6
    x441 = x3 * x435
    x442 = x3 * x341
    x443 = x3 * x331
    x444 = x343 + x443
    x445 = x38 * (x399 - x442 + x444)
    x446 = x3 * x349 + x351
    x447 = x3 * x404 + x376 + x38 * (x375 - x443 + x446) - x435 * x6
    x448 = x423 - x424
    x449 = x4 * (x386 + x415 - x417)
    x450 = -x10 * x156
    x451 = x119 * x321 + x450
    x452 = -x412
    x453 = x3 * x391
    x454 = x6 * (x425 + x453)
    x455 = x3 * x439
    x456 = x3 * x364
    x457 = x366 + x442
    x458 = x38 * (x429 - x456 + x457)
    x459 = x400 - x440 + x441 + x445
    x460 = x3 * x323
    x461 = x3 * x328
    x462 = x4 * (x330 + x339 + x460 - x461)
    x463 = x3 * x444
    x464 = x457 * x6
    x465 = x174 * x176
    x466 = x181 * x8
    x467 = x19 - x466
    x468 = x465 + x467
    x469 = x468 * x5
    x470 = x174 * x181
    x471 = x191 * x8
    x472 = x29 - x471
    x473 = x470 + x472
    x474 = x473 * x6
    x475 = -x474
    x476 = x469 + x475
    x477 = x476 * x5
    x478 = -x176 * x8 + x45
    x479 = x174 * x200 + x478
    x480 = -x468 * x6
    x481 = x479 * x5 + x480
    x482 = x466 + x50
    x483 = x4 * (-x465 + x479 + x482)
    x484 = -x476 * x6 + x483
    x485 = x481 * x5 + x484
    x486 = x473 * x5
    x487 = x174 * x191
    x488 = x215 * x8
    x489 = -x488 + x56
    x490 = x487 + x489
    x491 = x490 * x6
    x492 = -x491
    x493 = x486 + x492
    x494 = x493 * x6
    x495 = x471 + x69
    x496 = x4 * (x468 - x470 + x495)
    x497 = -x496
    x498 = x494 + x497
    x499 = -x494 + x496
    x500 = x477 + x499
    x501 = x4 * (-x469 + x474 + x481)
    x502 = x3 * x485 - x500 * x6 + 2.0 * x501
    x503 = x3 * x500
    x504 = x493 * x5
    x505 = x488 + x95
    x506 = x4 * (x473 - x487 + x505)
    x507 = x490 * x5
    x508 = x174 * x215
    x509 = x239 * x8
    x510 = -x509 + x82
    x511 = x508 + x510
    x512 = x511 * x6
    x513 = -x512
    x514 = x507 + x513
    x515 = x514 * x6
    x516 = x506 - x515
    x517 = x504 + x516
    x518 = x517 * x6
    x519 = x4 * (x476 - x486 + x491)
    x520 = 2.0 * x519
    x521 = x503 - x518 + x520
    x522 = x3 * x476
    x523 = x3 * x481 + x484
    x524 = x38 * (x498 - x522 + x523)
    x525 = x3 * x502 + x4 * (-x477 + x485 + x498) - x521 * x6 + x524
    x526 = -x506
    x527 = x515 + x526
    x528 = x4 * (x500 - x504 + x527)
    x529 = x3 * x521
    x530 = x3 * x517
    x531 = x5 * x514
    x532 = x129 + x509
    x533 = x4 * (x490 - x508 + x532)
    x534 = x5 * x511
    x535 = x174 * x239
    x536 = x271 * x8
    x537 = x117 - x536
    x538 = x535 + x537
    x539 = x538 * x6
    x540 = x534 - x539
    x541 = x540 * x6
    x542 = x533 - x541
    x543 = x531 + x542
    x544 = x543 * x6
    x545 = x4 * (x493 - x507 + x512)
    x546 = 2.0 * x545
    x547 = x530 - x544 + x546
    x548 = x547 * x6
    x549 = x3 * x493
    x550 = x499 + x522
    x551 = x38 * (x527 - x549 + x550)
    x552 = x528 + x529 - x548 + x551
    x553 = x516 + x549
    x554 = x553 * x6
    x555 = x3 * x473
    x556 = x3 * x468
    x557 = x475 + x556
    x558 = x4 * (x491 - x555 + x557)
    x559 = x3 * x550
    x560 = x3 * x479 + x480
    x561 = x3 * x523 + x4 * (x474 - x556 + x560) + x501 - x550 * x6
    x562 = -x533
    x563 = x541 + x562
    x564 = x157 + x536
    x565 = x4 * (x511 - x535 + x564)
    x566 = x155 - x304 * x8
    x567 = x174 * x271 + x566
    x568 = x3 * x514
    x569 = x6 * (x542 + x568)
    x570 = x3 * x490
    x571 = x492 + x555
    x572 = x4 * (x512 - x570 + x571)
    x573 = x3 * x553
    x574 = x519 - x554 + x558 + x559
    x575 = x571 * x6
    x576 = x3 * x557
    x577 = -x551
    x578 = 23.5542377588857 * x170
    x579 = x174 * x323
    x580 = x328 * x8
    x581 = -x580
    x582 = x579 + x581
    x583 = x5 * x582
    x584 = x174 * x328
    x585 = x338 * x8
    x586 = -x585
    x587 = x584 + x586
    x588 = x587 * x6
    x589 = -x588
    x590 = x583 + x589
    x591 = x5 * x590
    x592 = -x323 * x8
    x593 = x174 * x347 + x592
    x594 = -x582 * x6
    x595 = x5 * x593 + x594
    x596 = x4 * (-x579 + x580 + x593)
    x597 = -x590 * x6 + x596
    x598 = x5 * x595 + x597
    x599 = x5 * x587
    x600 = x174 * x338
    x601 = x362 * x8
    x602 = -x601
    x603 = x600 + x602
    x604 = x6 * x603
    x605 = -x604
    x606 = x599 + x605
    x607 = x6 * x606
    x608 = x4 * (x582 - x584 + x585)
    x609 = -x608
    x610 = x607 + x609
    x611 = -x607 + x608
    x612 = x591 + x611
    x613 = x4 * (-x583 + x588 + x595)
    x614 = x3 * x598 - x6 * x612 + 2.0 * x613
    x615 = x3 * x612
    x616 = x5 * x606
    x617 = x4 * (x587 - x600 + x601)
    x618 = x5 * x603
    x619 = x174 * x362
    x620 = x386 * x8
    x621 = -x620
    x622 = x619 + x621
    x623 = x6 * x622
    x624 = -x623
    x625 = x618 + x624
    x626 = x6 * x625
    x627 = x617 - x626
    x628 = x616 + x627
    x629 = x6 * x628
    x630 = x4 * (x590 - x599 + x604)
    x631 = 2.0 * x630
    x632 = x615 - x629 + x631
    x633 = x3 * x590
    x634 = x3 * x595 + x597
    x635 = x38 * (x610 - x633 + x634)
    x636 = x3 * x614 + x4 * (-x591 + x598 + x610) - x6 * x632 + x635
    x637 = -x617
    x638 = x626 + x637
    x639 = x4 * (x612 - x616 + x638)
    x640 = x3 * x632
    x641 = x3 * x628
    x642 = x5 * x625
    x643 = x4 * (x603 - x619 + x620)
    x644 = x5 * x622
    x645 = x174 * x386
    x646 = x418 * x8
    x647 = -x646
    x648 = x645 + x647
    x649 = x6 * x648
    x650 = x644 - x649
    x651 = x6 * x650
    x652 = x643 - x651
    x653 = x642 + x652
    x654 = x6 * x653
    x655 = x4 * (x606 - x618 + x623)
    x656 = 2.0 * x655
    x657 = x641 - x654 + x656
    x658 = x6 * x657
    x659 = x3 * x606
    x660 = x611 + x633
    x661 = x38 * (x638 - x659 + x660)
    x662 = x639 + x640 - x658 + x661
    x663 = x627 + x659
    x664 = x6 * x663
    x665 = x3 * x587
    x666 = x3 * x582
    x667 = x589 + x666
    x668 = x4 * (x604 - x665 + x667)
    x669 = x3 * x660
    x670 = x3 * x593 + x594
    x671 = x3 * x634 + x4 * (x588 - x666 + x670) - x6 * x660 + x613
    x672 = -x643
    x673 = x651 + x672
    x674 = x4 * (x622 - x645 + x646)
    x675 = -x451 * x8
    x676 = x174 * x418 + x675
    x677 = x3 * x625
    x678 = x6 * (x652 + x677)
    x679 = x3 * x603
    x680 = x605 + x665
    x681 = x4 * (x623 - x679 + x680)
    x682 = x3 * x663
    x683 = x630 - x664 + x668 + x669
    x684 = x6 * x680
    x685 = x3 * x667
    x686 = -x661
    x687 = 40.7971365319473 * x170
    x688 = x321 * x323
    x689 = x10 * x328
    x690 = x19 - x689
    x691 = x688 + x690
    x692 = x5 * x691
    x693 = x321 * x328
    x694 = x10 * x338
    x695 = x29 - x694
    x696 = x693 + x695
    x697 = x6 * x696
    x698 = -x697
    x699 = x692 + x698
    x700 = x5 * x699
    x701 = -x10 * x323 + x45
    x702 = x321 * x347 + x701
    x703 = -x6 * x691
    x704 = x5 * x702 + x703
    x705 = x50 + x689
    x706 = x4 * (-x688 + x702 + x705)
    x707 = -x6 * x699 + x706
    x708 = x5 * x704 + x707
    x709 = x5 * x696
    x710 = x321 * x338
    x711 = x10 * x362
    x712 = x56 - x711
    x713 = x710 + x712
    x714 = x6 * x713
    x715 = -x714
    x716 = x709 + x715
    x717 = x6 * x716
    x718 = x69 + x694
    x719 = x4 * (x691 - x693 + x718)
    x720 = -x719
    x721 = x717 + x720
    x722 = -x717 + x719
    x723 = x700 + x722
    x724 = x4 * (-x692 + x697 + x704)
    x725 = x3 * x708 - x6 * x723 + 2.0 * x724
    x726 = x3 * x723
    x727 = x5 * x716
    x728 = x711 + x95
    x729 = x4 * (x696 - x710 + x728)
    x730 = x5 * x713
    x731 = x321 * x362
    x732 = x10 * x386
    x733 = -x732 + x82
    x734 = x731 + x733
    x735 = x6 * x734
    x736 = -x735
    x737 = x730 + x736
    x738 = x6 * x737
    x739 = x729 - x738
    x740 = x727 + x739
    x741 = x6 * x740
    x742 = x4 * (x699 - x709 + x714)
    x743 = 2.0 * x742
    x744 = x726 - x741 + x743
    x745 = x3 * x699
    x746 = x3 * x704 + x707
    x747 = x38 * (x721 - x745 + x746)
    x748 = x3 * x725 + x4 * (-x700 + x708 + x721) - x6 * x744 + x747
    x749 = -x729
    x750 = x738 + x749
    x751 = x4 * (x723 - x727 + x750)
    x752 = x3 * x744
    x753 = x3 * x740
    x754 = x5 * x737
    x755 = x129 + x732
    x756 = x4 * (x713 - x731 + x755)
    x757 = x5 * x734
    x758 = x321 * x386
    x759 = x10 * x418
    x760 = x117 - x759
    x761 = x758 + x760
    x762 = x6 * x761
    x763 = x757 - x762
    x764 = x6 * x763
    x765 = x756 - x764
    x766 = x754 + x765
    x767 = x6 * x766
    x768 = x4 * (x716 - x730 + x735)
    x769 = 2.0 * x768
    x770 = x753 - x767 + x769
    x771 = x6 * x770
    x772 = x3 * x716
    x773 = x722 + x745
    x774 = x38 * (x750 - x772 + x773)
    x775 = x751 + x752 - x771 + x774
    x776 = x739 + x772
    x777 = x6 * x776
    x778 = x3 * x696
    x779 = x3 * x691
    x780 = x698 + x779
    x781 = x4 * (x714 - x778 + x780)
    x782 = x3 * x773
    x783 = x3 * x702 + x703
    x784 = x3 * x746 + x4 * (x697 - x779 + x783) - x6 * x773 + x724
    x785 = -x756
    x786 = x764 + x785
    x787 = x157 + x759
    x788 = x4 * (x734 - x758 + x787)
    x789 = -x10 * x451 + x155
    x790 = x321 * x418 + x789
    x791 = x3 * x737
    x792 = x6 * (x765 + x791)
    x793 = x3 * x713
    x794 = x715 + x778
    x795 = x4 * (x735 - x793 + x794)
    x796 = x3 * x776
    x797 = x742 - x777 + x781 + x782
    x798 = x6 * x794
    x799 = x3 * x780
    x800 = -x774
    x801 = x174 * x468
    x802 = 2.0 * x203 - x468 * x8
    x803 = x174 * x479 + x802
    x804 = x473 * x8
    x805 = 2.0 * x186
    x806 = x804 - x805
    x807 = x4 * (-x801 + x803 + x806)
    x808 = -x804 + x805
    x809 = x801 + x808
    x810 = -x6 * x809
    x811 = x3 * x803 + x810
    x812 = x3 * x809
    x813 = x174 * x473
    x814 = x490 * x8
    x815 = 2.0 * x210
    x816 = -x814 + x815
    x817 = x813 + x816
    x818 = x6 * x817
    x819 = -x818
    x820 = x812 + x819
    x821 = x3 * x811 - x6 * x820 + x807
    x822 = x174 * x490
    x823 = x511 * x8
    x824 = 2.0 * x234
    x825 = x823 - x824
    x826 = x4 * (x817 - x822 + x825)
    x827 = x3 * x817
    x828 = -x823 + x824
    x829 = x822 + x828
    x830 = x6 * x829
    x831 = -x830
    x832 = x827 + x831
    x833 = x3 * x832
    x834 = x3 * x829
    x835 = x174 * x511
    x836 = x538 * x8
    x837 = 2.0 * x265
    x838 = -x836 + x837
    x839 = x835 + x838
    x840 = x6 * x839
    x841 = -x840
    x842 = x6 * (x834 + x841)
    x843 = x814 - x815
    x844 = x4 * (x809 - x813 + x843)
    x845 = x3 * x820
    x846 = x6 * x832
    x847 = x844 + x845 - x846
    x848 = x4 * (x820 - x827 + x830)
    x849 = x4 * (x811 - x812 + x818)
    x850 = -x844
    x851 = x5 * x809
    x852 = x5 * x803 + x810
    x853 = x819 + x851
    x854 = x3 * x852 - x6 * x853 + x807
    x855 = x3 * x853
    x856 = x5 * x817
    x857 = x831 + x856
    x858 = x6 * x857
    x859 = x844 + x855 - x858
    x860 = x3 * x854 + x4 * (x818 - x851 + x852) - x6 * x859 + x849
    x861 = x4 * (x830 + x853 - x856)
    x862 = x3 * x859
    x863 = x3 * x857
    x864 = x5 * x829
    x865 = x841 + x864
    x866 = x6 * x865
    x867 = x826 + x863 - x866
    x868 = x6 * x867
    x869 = x848 + x861 + x862 - x868
    x870 = -x826
    x871 = x836 - x837
    x872 = x4 * (x829 - x835 + x871)
    x873 = 2.0 * x302 - x567 * x8
    x874 = x174 * x538 + x873
    x875 = x174 * x582
    x876 = x350 - x582 * x8
    x877 = x174 * x593 + x876
    x878 = x587 * x8
    x879 = x374 + x878
    x880 = x4 * (-x875 + x877 + x879)
    x881 = x333 - x878
    x882 = x875 + x881
    x883 = -x6 * x882
    x884 = x3 * x877 + x883
    x885 = x3 * x882
    x886 = x174 * x587
    x887 = x603 * x8
    x888 = x357 - x887
    x889 = x886 + x888
    x890 = x6 * x889
    x891 = -x890
    x892 = x885 + x891
    x893 = x3 * x884 - x6 * x892 + x880
    x894 = x174 * x603
    x895 = x622 * x8
    x896 = x428 + x895
    x897 = x4 * (x889 - x894 + x896)
    x898 = x3 * x889
    x899 = x381 - x895
    x900 = x894 + x899
    x901 = x6 * x900
    x902 = -x901
    x903 = x898 + x902
    x904 = x3 * x903
    x905 = x3 * x900
    x906 = x174 * x622
    x907 = x648 * x8
    x908 = x412 - x907
    x909 = x906 + x908
    x910 = x6 * x909
    x911 = -x910
    x912 = x6 * (x905 + x911)
    x913 = x398 + x887
    x914 = x4 * (x882 - x886 + x913)
    x915 = x3 * x892
    x916 = x6 * x903
    x917 = x914 + x915 - x916
    x918 = x4 * (x892 - x898 + x901)
    x919 = x4 * (x884 - x885 + x890)
    x920 = -x914
    x921 = x5 * x882
    x922 = x5 * x877 + x883
    x923 = x891 + x921
    x924 = x3 * x922 - x6 * x923 + x880
    x925 = x3 * x923
    x926 = x5 * x889
    x927 = x902 + x926
    x928 = x6 * x927
    x929 = x914 + x925 - x928
    x930 = x3 * x924 + x4 * (x890 - x921 + x922) - x6 * x929 + x919
    x931 = x4 * (x901 + x923 - x926)
    x932 = x3 * x929
    x933 = x3 * x927
    x934 = x5 * x900
    x935 = x911 + x934
    x936 = x6 * x935
    x937 = x897 + x933 - x936
    x938 = x6 * x937
    x939 = x918 + x931 + x932 - x938
    x940 = -x897
    x941 = x452 + x907
    x942 = x4 * (x900 - x906 + x941)
    x943 = x449 - x676 * x8
    x944 = x174 * x648 + x943
    x945 = x696 * x8
    x946 = x174 * x691
    x947 = -x691 * x8
    x948 = x174 * x702 + x947
    x949 = x4 * (x945 - x946 + x948)
    x950 = -x945
    x951 = x946 + x950
    x952 = -x6 * x951
    x953 = x3 * x948 + x952
    x954 = x3 * x951
    x955 = x174 * x696
    x956 = x713 * x8
    x957 = -x956
    x958 = x955 + x957
    x959 = x6 * x958
    x960 = -x959
    x961 = x954 + x960
    x962 = x3 * x953 - x6 * x961 + x949
    x963 = x734 * x8
    x964 = x174 * x713
    x965 = x4 * (x958 + x963 - x964)
    x966 = x3 * x958
    x967 = -x963
    x968 = x964 + x967
    x969 = x6 * x968
    x970 = -x969
    x971 = x966 + x970
    x972 = x3 * x971
    x973 = x3 * x968
    x974 = x174 * x734
    x975 = x761 * x8
    x976 = -x975
    x977 = x974 + x976
    x978 = x6 * x977
    x979 = -x978
    x980 = x6 * (x973 + x979)
    x981 = x4 * (x951 - x955 + x956)
    x982 = x3 * x961
    x983 = x6 * x971
    x984 = x981 + x982 - x983
    x985 = x4 * (x961 - x966 + x969)
    x986 = x4 * (x953 - x954 + x959)
    x987 = -x981
    x988 = x5 * x951
    x989 = x5 * x948 + x952
    x990 = x960 + x988
    x991 = x3 * x989 - x6 * x990 + x949
    x992 = x3 * x990
    x993 = x5 * x958
    x994 = x970 + x993
    x995 = x6 * x994
    x996 = x981 + x992 - x995
    x997 = x3 * x991 + x4 * (x959 - x988 + x989) - x6 * x996 + x986
    x998 = x4 * (x969 + x990 - x993)
    x999 = x3 * x996
    x1000 = x3 * x994
    x1001 = x5 * x968
    x1002 = x1001 + x979
    x1003 = x1002 * x6
    x1004 = x1000 - x1003 + x965
    x1005 = x1004 * x6
    x1006 = -x1005 + x985 + x998 + x999
    x1007 = -x965
    x1008 = x4 * (x968 - x974 + x975)
    x1009 = -x790 * x8
    x1010 = x1009 + x174 * x761
    x1011 = x321 * x691
    x1012 = -x10 * x691 + 2.0 * x350
    x1013 = x1012 + x321 * x702
    x1014 = x10 * x696
    x1015 = 2.0 * x333
    x1016 = x1014 - x1015
    x1017 = x4 * (-x1011 + x1013 + x1016)
    x1018 = -x1014 + x1015
    x1019 = x1011 + x1018
    x1020 = -x1019 * x6
    x1021 = x1013 * x3 + x1020
    x1022 = x1019 * x3
    x1023 = x321 * x696
    x1024 = x10 * x713
    x1025 = 2.0 * x357
    x1026 = -x1024 + x1025
    x1027 = x1023 + x1026
    x1028 = x1027 * x6
    x1029 = -x1028
    x1030 = x1022 + x1029
    x1031 = x1017 + x1021 * x3 - x1030 * x6
    x1032 = x321 * x713
    x1033 = x10 * x734
    x1034 = 2.0 * x381
    x1035 = x1033 - x1034
    x1036 = x4 * (x1027 - x1032 + x1035)
    x1037 = x1027 * x3
    x1038 = -x1033 + x1034
    x1039 = x1032 + x1038
    x1040 = x1039 * x6
    x1041 = -x1040
    x1042 = x1037 + x1041
    x1043 = x1042 * x3
    x1044 = x1039 * x3
    x1045 = x321 * x734
    x1046 = x10 * x761
    x1047 = 2.0 * x412
    x1048 = -x1046 + x1047
    x1049 = x1045 + x1048
    x1050 = x1049 * x6
    x1051 = -x1050
    x1052 = x6 * (x1044 + x1051)
    x1053 = x1024 - x1025
    x1054 = x4 * (x1019 - x1023 + x1053)
    x1055 = x1030 * x3
    x1056 = x1042 * x6
    x1057 = x1054 + x1055 - x1056
    x1058 = x4 * (x1030 - x1037 + x1040)
    x1059 = x4 * (x1021 - x1022 + x1028)
    x1060 = -x1054
    x1061 = x1019 * x5
    x1062 = x1013 * x5 + x1020
    x1063 = x1029 + x1061
    x1064 = x1017 + x1062 * x3 - x1063 * x6
    x1065 = x1063 * x3
    x1066 = x1027 * x5
    x1067 = x1041 + x1066
    x1068 = x1067 * x6
    x1069 = x1054 + x1065 - x1068
    x1070 = x1059 + x1064 * x3 - x1069 * x6 + x4 * (x1028 - x1061 + x1062)
    x1071 = x4 * (x1040 + x1063 - x1066)
    x1072 = x1069 * x3
    x1073 = x1067 * x3
    x1074 = x1039 * x5
    x1075 = x1051 + x1074
    x1076 = x1075 * x6
    x1077 = x1036 + x1073 - x1076
    x1078 = x1077 * x6
    x1079 = x1058 + x1071 + x1072 - x1078
    x1080 = -x1036
    x1081 = x1046 - x1047
    x1082 = x4 * (x1039 - x1045 + x1081)
    x1083 = -x10 * x790 + 2.0 * x449
    x1084 = x1083 + x321 * x761
    x1085 = x174 * x809
    x1086 = 3.0 * x483 - x8 * x809
    x1087 = x1086 + x174 * x803
    x1088 = x8 * x817
    x1089 = 3.0 * x496
    x1090 = x1088 - x1089
    x1091 = x4 * (-x1085 + x1087 + x1090)
    x1092 = -x1088 + x1089
    x1093 = x1085 + x1092
    x1094 = x1087 * x3 - x1093 * x6
    x1095 = x1093 * x3
    x1096 = x174 * x817
    x1097 = x8 * x829
    x1098 = 3.0 * x506
    x1099 = -x1097 + x1098
    x1100 = x1096 + x1099
    x1101 = x1100 * x6
    x1102 = x1095 - x1101
    x1103 = x1091 + x1094 * x3 - x1102 * x6
    x1104 = x1097 - x1098
    x1105 = x4 * (x1093 - x1096 + x1104)
    x1106 = x1102 * x3
    x1107 = x1100 * x3
    x1108 = x174 * x829
    x1109 = x8 * x839
    x1110 = 3.0 * x533
    x1111 = -x1109 + x1110
    x1112 = x1108 + x1111
    x1113 = x1112 * x6
    x1114 = x1107 - x1113
    x1115 = x1114 * x6
    x1116 = x1105 + x1106 - x1115
    x1117 = x1109 - x1110
    x1118 = x4 * (x1100 - x1108 + x1117)
    x1119 = 3.0 * x565 - x8 * x874
    x1120 = x1119 + x174 * x839
    x1121 = -x1105
    x1122 = x174 * x882
    x1123 = 2.0 * x596 - x8 * x882
    x1124 = x1123 + x174 * x877
    x1125 = x8 * x889
    x1126 = 2.0 * x608
    x1127 = x1125 - x1126
    x1128 = x4 * (-x1122 + x1124 + x1127)
    x1129 = -x1125 + x1126
    x1130 = x1122 + x1129
    x1131 = x1124 * x3 - x1130 * x6
    x1132 = x1130 * x3
    x1133 = x174 * x889
    x1134 = x8 * x900
    x1135 = 2.0 * x617
    x1136 = -x1134 + x1135
    x1137 = x1133 + x1136
    x1138 = x1137 * x6
    x1139 = x1132 - x1138
    x1140 = x1128 + x1131 * x3 - x1139 * x6
    x1141 = x1134 - x1135
    x1142 = x4 * (x1130 - x1133 + x1141)
    x1143 = x1139 * x3
    x1144 = x1137 * x3
    x1145 = x174 * x900
    x1146 = x8 * x909
    x1147 = 2.0 * x643
    x1148 = -x1146 + x1147
    x1149 = x1145 + x1148
    x1150 = x1149 * x6
    x1151 = x1144 - x1150
    x1152 = x1151 * x6
    x1153 = x1142 + x1143 - x1152
    x1154 = x1146 - x1147
    x1155 = x4 * (x1137 - x1145 + x1154)
    x1156 = 2.0 * x674 - x8 * x944
    x1157 = x1156 + x174 * x909
    x1158 = -x1142
    x1159 = x174 * x951
    x1160 = x706 - x8 * x951
    x1161 = x1160 + x174 * x948
    x1162 = x8 * x958
    x1163 = x1162 + x720
    x1164 = x4 * (-x1159 + x1161 + x1163)
    x1165 = -x1162 + x719
    x1166 = x1159 + x1165
    x1167 = x1161 * x3 - x1166 * x6
    x1168 = x1166 * x3
    x1169 = x174 * x958
    x1170 = x8 * x968
    x1171 = -x1170 + x729
    x1172 = x1169 + x1171
    x1173 = x1172 * x6
    x1174 = x1168 - x1173
    x1175 = x1164 + x1167 * x3 - x1174 * x6
    x1176 = x1170 + x749
    x1177 = x4 * (x1166 - x1169 + x1176)
    x1178 = x1174 * x3
    x1179 = x1172 * x3
    x1180 = x174 * x968
    x1181 = x8 * x977
    x1182 = -x1181 + x756
    x1183 = x1180 + x1182
    x1184 = x1183 * x6
    x1185 = x1179 - x1184
    x1186 = x1185 * x6
    x1187 = x1177 + x1178 - x1186
    x1188 = x1181 + x785
    x1189 = x4 * (x1172 - x1180 + x1188)
    x1190 = -x1010 * x8 + x788
    x1191 = x1190 + x174 * x977
    x1192 = -x1177
    x1193 = x1027 * x8
    x1194 = x1019 * x174
    x1195 = -x1019 * x8
    x1196 = x1013 * x174 + x1195
    x1197 = x4 * (x1193 - x1194 + x1196)
    x1198 = -x1193
    x1199 = x1194 + x1198
    x1200 = x1196 * x3 - x1199 * x6
    x1201 = x1199 * x3
    x1202 = x1027 * x174
    x1203 = x1039 * x8
    x1204 = -x1203
    x1205 = x1202 + x1204
    x1206 = x1205 * x6
    x1207 = x1201 - x1206
    x1208 = x1197 + x1200 * x3 - x1207 * x6
    x1209 = x4 * (x1199 - x1202 + x1203)
    x1210 = x1207 * x3
    x1211 = x1205 * x3
    x1212 = x1039 * x174
    x1213 = x1049 * x8
    x1214 = -x1213
    x1215 = x1212 + x1214
    x1216 = x1215 * x6
    x1217 = x1211 - x1216
    x1218 = x1217 * x6
    x1219 = x1209 + x1210 - x1218
    x1220 = x4 * (x1205 - x1212 + x1213)
    x1221 = -x1084 * x8
    x1222 = x1049 * x174 + x1221
    x1223 = -x1209
    x1224 = x1019 * x321
    x1225 = -x10 * x1019 + 3.0 * x706
    x1226 = x1013 * x321 + x1225
    x1227 = x10 * x1027
    x1228 = 3.0 * x719
    x1229 = x1227 - x1228
    x1230 = x4 * (-x1224 + x1226 + x1229)
    x1231 = -x1227 + x1228
    x1232 = x1224 + x1231
    x1233 = x1226 * x3 - x1232 * x6
    x1234 = x1232 * x3
    x1235 = x1027 * x321
    x1236 = x10 * x1039
    x1237 = 3.0 * x729
    x1238 = -x1236 + x1237
    x1239 = x1235 + x1238
    x1240 = x1239 * x6
    x1241 = x1234 - x1240
    x1242 = x1230 + x1233 * x3 - x1241 * x6
    x1243 = x1236 - x1237
    x1244 = x4 * (x1232 - x1235 + x1243)
    x1245 = x1241 * x3
    x1246 = x1239 * x3
    x1247 = x1039 * x321
    x1248 = x10 * x1049
    x1249 = 3.0 * x756
    x1250 = -x1248 + x1249
    x1251 = x1247 + x1250
    x1252 = x1251 * x6
    x1253 = x1246 - x1252
    x1254 = x1253 * x6
    x1255 = x1244 + x1245 - x1254
    x1256 = x1248 - x1249
    x1257 = x4 * (x1239 - x1247 + x1256)
    x1258 = -x10 * x1084 + 3.0 * x788
    x1259 = x1049 * x321 + x1258
    x1260 = -x1244
    x1261 = -x7 - A[1]
    x1262 = x1261 * x13
    x1263 = x1261 * x15
    x1264 = x1263 + x173
    x1265 = x4 * (-x1262 + x1264 + x178)
    x1266 = x1264 * x5
    x1267 = x1262 + x179
    x1268 = x1267 * x6
    x1269 = x1266 - x1268
    x1270 = x1269 * x5
    x1271 = x1267 * x5
    x1272 = x1261 * x20
    x1273 = x1272 + x189
    x1274 = x1273 * x6
    x1275 = x1271 - x1274
    x1276 = x1275 * x6
    x1277 = x1265 + x1270 - x1276
    x1278 = x1277 * x5
    x1279 = x4 * (x1267 - x1272 + x188)
    x1280 = x1275 * x5
    x1281 = x1273 * x5
    x1282 = x1261 * x30
    x1283 = x1282 + x213
    x1284 = x1283 * x6
    x1285 = x1281 - x1284
    x1286 = x1285 * x6
    x1287 = x1279 + x1280 - x1286
    x1288 = x1287 * x6
    x1289 = x38 * (x1269 - x1271 + x1274)
    x1290 = -x1288 + x1289
    x1291 = x1278 + x1290
    x1292 = x1291 * x5
    x1293 = x1261 * x43 + x199
    x1294 = x4 * (-x1263 + x1293 + x172)
    x1295 = -x1264 * x6 + x1293 * x5
    x1296 = -x1269 * x6 + x1294 + x1295 * x5
    x1297 = -x1277 * x6 + x38 * (-x1266 + x1268 + x1295)
    x1298 = x1296 * x5 + x1297
    x1299 = -x1265
    x1300 = -x1291 * x6 + x130 * (-x1270 + x1276 + x1296 + x1299)
    x1301 = x1298 * x5 + x1300
    x1302 = x1287 * x5
    x1303 = x4 * (x1273 - x1282 + x212)
    x1304 = x1285 * x5
    x1305 = x1283 * x5
    x1306 = x1261 * x58
    x1307 = x1306 + x237
    x1308 = x1307 * x6
    x1309 = x1305 - x1308
    x1310 = x1309 * x6
    x1311 = x1303 + x1304 - x1310
    x1312 = x1311 * x6
    x1313 = x38 * (x1275 - x1281 + x1284)
    x1314 = x1302 - x1312 + x1313
    x1315 = x1314 * x6
    x1316 = -x1279
    x1317 = x130 * (x1277 - x1280 + x1286 + x1316)
    x1318 = x1315 - x1317
    x1319 = -x1315 + x1317
    x1320 = x1292 + x1319
    x1321 = x1288 - x1289
    x1322 = x4 * (-x1278 + x1298 + x1321)
    x1323 = x1301 * x3 - x1320 * x6 + 4.0 * x1322
    x1324 = x1320 * x3
    x1325 = x1314 * x5
    x1326 = x1311 * x5
    x1327 = x4 * (x1283 - x1306 + x236)
    x1328 = x1309 * x5
    x1329 = x1307 * x5
    x1330 = x1261 * x84
    x1331 = x1330 + x269
    x1332 = x1331 * x6
    x1333 = x1329 - x1332
    x1334 = x1333 * x6
    x1335 = x1327 + x1328 - x1334
    x1336 = x1335 * x6
    x1337 = x38 * (x1285 - x1305 + x1308)
    x1338 = x1326 - x1336 + x1337
    x1339 = x1338 * x6
    x1340 = -x1303
    x1341 = x130 * (x1287 - x1304 + x1310 + x1340)
    x1342 = -x1339 + x1341
    x1343 = x1325 + x1342
    x1344 = x1343 * x6
    x1345 = x1312 - x1313
    x1346 = x4 * (x1291 - x1302 + x1345)
    x1347 = 4.0 * x1346
    x1348 = x1324 - x1344 + x1347
    x1349 = x1291 * x3
    x1350 = x1298 * x3 + x1300
    x1351 = x1339 - x1341
    x1352 = x4 * (x1307 - x1330 + x268)
    x1353 = x119 * x1261 + x303
    x1354 = -x1327
    x1355 = x1314 * x3
    x1356 = x1319 + x1349
    x1357 = x1277 * x3
    x1358 = x1261 * x176
    x1359 = x1358 + x467
    x1360 = x1359 * x5
    x1361 = x1261 * x181
    x1362 = x1361 + x472
    x1363 = x1362 * x6
    x1364 = x1360 - x1363
    x1365 = x1364 * x5
    x1366 = x4 * (x1359 - x1361 + x495)
    x1367 = x1362 * x5
    x1368 = x1261 * x191
    x1369 = x1368 + x489
    x1370 = x1369 * x6
    x1371 = x1367 - x1370
    x1372 = x1371 * x6
    x1373 = x1366 - x1372
    x1374 = x1365 + x1373
    x1375 = x1374 * x5
    x1376 = x1261 * x200 + x478
    x1377 = -x1359 * x6 + x1376 * x5
    x1378 = x4 * (-x1358 + x1376 + x482)
    x1379 = -x1364 * x6 + x1378
    x1380 = x1377 * x5 + x1379
    x1381 = -x1374 * x6 + x38 * (-x1360 + x1363 + x1377)
    x1382 = x1380 * x5 + x1381
    x1383 = x4 * (x1362 - x1368 + x505)
    x1384 = x1371 * x5
    x1385 = x1369 * x5
    x1386 = x1261 * x215
    x1387 = x1386 + x510
    x1388 = x1387 * x6
    x1389 = x1385 - x1388
    x1390 = x1389 * x6
    x1391 = x1383 + x1384 - x1390
    x1392 = x1391 * x6
    x1393 = x38 * (x1364 - x1367 + x1370)
    x1394 = x1392 - x1393
    x1395 = -x1392 + x1393
    x1396 = x1375 + x1395
    x1397 = -x1366 + x1372
    x1398 = x4 * (-x1365 + x1380 + x1397)
    x1399 = x1382 * x3 - x1396 * x6 + 3.0 * x1398
    x1400 = x1396 * x3
    x1401 = x1391 * x5
    x1402 = x4 * (x1369 - x1386 + x532)
    x1403 = x1389 * x5
    x1404 = x1387 * x5
    x1405 = x1261 * x239
    x1406 = x1405 + x537
    x1407 = x1406 * x6
    x1408 = x1404 - x1407
    x1409 = x1408 * x6
    x1410 = x1402 + x1403 - x1409
    x1411 = x1410 * x6
    x1412 = x38 * (x1371 - x1385 + x1388)
    x1413 = -x1411 + x1412
    x1414 = x1401 + x1413
    x1415 = x1414 * x6
    x1416 = -x1383 + x1390
    x1417 = x4 * (x1374 - x1384 + x1416)
    x1418 = 3.0 * x1417
    x1419 = x1400 - x1415 + x1418
    x1420 = x1374 * x3
    x1421 = x1380 * x3 + x1381
    x1422 = x1411 - x1412
    x1423 = x4 * (x1387 - x1405 + x564)
    x1424 = x1261 * x271 + x566
    x1425 = x1391 * x3
    x1426 = x1395 + x1420
    x1427 = x1364 * x3
    x1428 = 48.2718229290016 * x170
    x1429 = x1261 * x323
    x1430 = x1429 + x581
    x1431 = x1430 * x5
    x1432 = x1261 * x328
    x1433 = x1432 + x586
    x1434 = x1433 * x6
    x1435 = x1431 - x1434
    x1436 = x1435 * x5
    x1437 = x4 * (x1430 - x1432 + x585)
    x1438 = x1433 * x5
    x1439 = x1261 * x338
    x1440 = x1439 + x602
    x1441 = x1440 * x6
    x1442 = x1438 - x1441
    x1443 = x1442 * x6
    x1444 = x1437 - x1443
    x1445 = x1436 + x1444
    x1446 = x1445 * x5
    x1447 = x1261 * x347 + x592
    x1448 = -x1430 * x6 + x1447 * x5
    x1449 = x4 * (-x1429 + x1447 + x580)
    x1450 = -x1435 * x6 + x1449
    x1451 = x1448 * x5 + x1450
    x1452 = -x1445 * x6 + x38 * (-x1431 + x1434 + x1448)
    x1453 = x1451 * x5 + x1452
    x1454 = x4 * (x1433 - x1439 + x601)
    x1455 = x1442 * x5
    x1456 = x1440 * x5
    x1457 = x1261 * x362
    x1458 = x1457 + x621
    x1459 = x1458 * x6
    x1460 = x1456 - x1459
    x1461 = x1460 * x6
    x1462 = x1454 + x1455 - x1461
    x1463 = x1462 * x6
    x1464 = x38 * (x1435 - x1438 + x1441)
    x1465 = x1463 - x1464
    x1466 = -x1463 + x1464
    x1467 = x1446 + x1466
    x1468 = -x1437
    x1469 = x1443 + x1468
    x1470 = x4 * (-x1436 + x1451 + x1469)
    x1471 = x1453 * x3 - x1467 * x6 + 3.0 * x1470
    x1472 = x1467 * x3
    x1473 = x1462 * x5
    x1474 = x4 * (x1440 - x1457 + x620)
    x1475 = x1460 * x5
    x1476 = x1458 * x5
    x1477 = x1261 * x386
    x1478 = x1477 + x647
    x1479 = x1478 * x6
    x1480 = x1476 - x1479
    x1481 = x1480 * x6
    x1482 = x1474 + x1475 - x1481
    x1483 = x1482 * x6
    x1484 = x38 * (x1442 - x1456 + x1459)
    x1485 = -x1483 + x1484
    x1486 = x1473 + x1485
    x1487 = x1486 * x6
    x1488 = -x1454
    x1489 = x1461 + x1488
    x1490 = x4 * (x1445 - x1455 + x1489)
    x1491 = 3.0 * x1490
    x1492 = x1472 - x1487 + x1491
    x1493 = x1445 * x3
    x1494 = x1451 * x3 + x1452
    x1495 = x1483 - x1484
    x1496 = x4 * (x1458 - x1477 + x646)
    x1497 = x1261 * x418 + x675
    x1498 = -x1474
    x1499 = x1462 * x3
    x1500 = x1466 + x1493
    x1501 = x1435 * x3
    x1502 = x1261 * x468
    x1503 = x1502 + x808
    x1504 = x1503 * x5
    x1505 = x1261 * x473
    x1506 = x1505 + x816
    x1507 = x1506 * x6
    x1508 = -x1507
    x1509 = x1504 + x1508
    x1510 = x1509 * x5
    x1511 = x1261 * x479 + x802
    x1512 = -x1503 * x6
    x1513 = x1511 * x5 + x1512
    x1514 = x4 * (-x1502 + x1511 + x806)
    x1515 = -x1509 * x6 + x1514
    x1516 = x1513 * x5 + x1515
    x1517 = x1506 * x5
    x1518 = x1261 * x490
    x1519 = x1518 + x828
    x1520 = x1519 * x6
    x1521 = x1517 - x1520
    x1522 = x1521 * x6
    x1523 = x4 * (x1503 - x1505 + x843)
    x1524 = x1522 - x1523
    x1525 = -x1522 + x1523
    x1526 = x1510 + x1525
    x1527 = x4 * (-x1504 + x1507 + x1513)
    x1528 = x1516 * x3 - x1526 * x6 + 2.0 * x1527
    x1529 = x1526 * x3
    x1530 = x1521 * x5
    x1531 = x4 * (x1506 - x1518 + x825)
    x1532 = x1519 * x5
    x1533 = x1261 * x511
    x1534 = x1533 + x838
    x1535 = x1534 * x6
    x1536 = x1532 - x1535
    x1537 = x1536 * x6
    x1538 = x1531 - x1537
    x1539 = x1530 + x1538
    x1540 = x1539 * x6
    x1541 = x4 * (x1509 - x1517 + x1520)
    x1542 = 2.0 * x1541
    x1543 = x1529 - x1540 + x1542
    x1544 = x1509 * x3
    x1545 = x1513 * x3 + x1515
    x1546 = -x1531 + x1537
    x1547 = x4 * (x1519 - x1533 + x871)
    x1548 = x1261 * x538 + x873
    x1549 = x1521 * x3
    x1550 = x1525 + x1544
    x1551 = x1503 * x3
    x1552 = 62.3186554316989 * x170
    x1553 = x1261 * x582
    x1554 = x1553 + x881
    x1555 = x1554 * x5
    x1556 = x1261 * x587
    x1557 = x1556 + x888
    x1558 = x1557 * x6
    x1559 = -x1558
    x1560 = x1555 + x1559
    x1561 = x1560 * x5
    x1562 = x1261 * x593 + x876
    x1563 = -x1554 * x6
    x1564 = x1562 * x5 + x1563
    x1565 = x4 * (-x1553 + x1562 + x879)
    x1566 = -x1560 * x6 + x1565
    x1567 = x1564 * x5 + x1566
    x1568 = x1557 * x5
    x1569 = x1261 * x603
    x1570 = x1569 + x899
    x1571 = x1570 * x6
    x1572 = x1568 - x1571
    x1573 = x1572 * x6
    x1574 = x4 * (x1554 - x1556 + x913)
    x1575 = x1573 - x1574
    x1576 = -x1573 + x1574
    x1577 = x1561 + x1576
    x1578 = x4 * (-x1555 + x1558 + x1564)
    x1579 = x1567 * x3 - x1577 * x6 + 2.0 * x1578
    x1580 = x1577 * x3
    x1581 = x1572 * x5
    x1582 = x4 * (x1557 - x1569 + x896)
    x1583 = x1570 * x5
    x1584 = x1261 * x622
    x1585 = x1584 + x908
    x1586 = x1585 * x6
    x1587 = x1583 - x1586
    x1588 = x1587 * x6
    x1589 = x1582 - x1588
    x1590 = x1581 + x1589
    x1591 = x1590 * x6
    x1592 = x4 * (x1560 - x1568 + x1571)
    x1593 = 2.0 * x1592
    x1594 = x1580 - x1591 + x1593
    x1595 = x1560 * x3
    x1596 = x1564 * x3 + x1566
    x1597 = -x1582 + x1588
    x1598 = x4 * (x1570 - x1584 + x941)
    x1599 = x1261 * x648 + x943
    x1600 = x1572 * x3
    x1601 = x1576 + x1595
    x1602 = x1554 * x3
    x1603 = 107.939077467081 * x170
    x1604 = x1261 * x691
    x1605 = x1604 + x950
    x1606 = x1605 * x5
    x1607 = x1261 * x696
    x1608 = x1607 + x957
    x1609 = x1608 * x6
    x1610 = -x1609
    x1611 = x1606 + x1610
    x1612 = x1611 * x5
    x1613 = x1261 * x702 + x947
    x1614 = -x1605 * x6
    x1615 = x1613 * x5 + x1614
    x1616 = x4 * (-x1604 + x1613 + x945)
    x1617 = -x1611 * x6 + x1616
    x1618 = x1615 * x5 + x1617
    x1619 = x1608 * x5
    x1620 = x1261 * x713
    x1621 = x1620 + x967
    x1622 = x1621 * x6
    x1623 = x1619 - x1622
    x1624 = x1623 * x6
    x1625 = x4 * (x1605 - x1607 + x956)
    x1626 = -x1625
    x1627 = x1624 + x1626
    x1628 = -x1624 + x1625
    x1629 = x1612 + x1628
    x1630 = x4 * (-x1606 + x1609 + x1615)
    x1631 = x1618 * x3 - x1629 * x6 + 2.0 * x1630
    x1632 = x1629 * x3
    x1633 = x1623 * x5
    x1634 = x4 * (x1608 - x1620 + x963)
    x1635 = x1621 * x5
    x1636 = x1261 * x734
    x1637 = x1636 + x976
    x1638 = x1637 * x6
    x1639 = x1635 - x1638
    x1640 = x1639 * x6
    x1641 = x1634 - x1640
    x1642 = x1633 + x1641
    x1643 = x1642 * x6
    x1644 = x4 * (x1611 - x1619 + x1622)
    x1645 = 2.0 * x1644
    x1646 = x1632 - x1643 + x1645
    x1647 = x1611 * x3
    x1648 = x1615 * x3 + x1617
    x1649 = -x1634
    x1650 = x1640 + x1649
    x1651 = x4 * (x1621 - x1636 + x975)
    x1652 = x1009 + x1261 * x761
    x1653 = x1623 * x3
    x1654 = x1628 + x1647
    x1655 = x1605 * x3
    x1656 = x1261 * x809
    x1657 = x1086 + x1261 * x803
    x1658 = x4 * (x1090 - x1656 + x1657)
    x1659 = x1092 + x1656
    x1660 = -x1659 * x6
    x1661 = x1657 * x3 + x1660
    x1662 = x1261 * x817
    x1663 = x1099 + x1662
    x1664 = x1663 * x3
    x1665 = x1261 * x829
    x1666 = x1111 + x1665
    x1667 = x1666 * x6
    x1668 = -x1667
    x1669 = x4 * (x1104 + x1659 - x1662)
    x1670 = -x1669
    x1671 = x1659 * x3
    x1672 = x1663 * x6
    x1673 = -x1672
    x1674 = x1671 + x1673
    x1675 = x1659 * x5
    x1676 = x1657 * x5 + x1660
    x1677 = x1673 + x1675
    x1678 = x1658 + x1676 * x3 - x1677 * x6
    x1679 = x1677 * x3
    x1680 = x1663 * x5
    x1681 = x1668 + x1680
    x1682 = x1681 * x6
    x1683 = x1669 + x1679 - x1682
    x1684 = x4 * (x1117 + x1663 - x1665)
    x1685 = x1119 + x1261 * x839
    x1686 = x1261 * x882
    x1687 = x1123 + x1261 * x877
    x1688 = x4 * (x1127 - x1686 + x1687)
    x1689 = x1129 + x1686
    x1690 = -x1689 * x6
    x1691 = x1687 * x3 + x1690
    x1692 = x1261 * x889
    x1693 = x1136 + x1692
    x1694 = x1693 * x3
    x1695 = x1261 * x900
    x1696 = x1148 + x1695
    x1697 = x1696 * x6
    x1698 = -x1697
    x1699 = x4 * (x1141 + x1689 - x1692)
    x1700 = -x1699
    x1701 = x1689 * x3
    x1702 = x1693 * x6
    x1703 = -x1702
    x1704 = x1701 + x1703
    x1705 = x1689 * x5
    x1706 = x1687 * x5 + x1690
    x1707 = x1703 + x1705
    x1708 = x1688 + x1706 * x3 - x1707 * x6
    x1709 = x1707 * x3
    x1710 = x1693 * x5
    x1711 = x1698 + x1710
    x1712 = x1711 * x6
    x1713 = x1699 + x1709 - x1712
    x1714 = x4 * (x1154 + x1693 - x1695)
    x1715 = x1156 + x1261 * x909
    x1716 = x1261 * x951
    x1717 = x1160 + x1261 * x948
    x1718 = x4 * (x1163 - x1716 + x1717)
    x1719 = x1165 + x1716
    x1720 = -x1719 * x6
    x1721 = x1717 * x3 + x1720
    x1722 = x1261 * x958
    x1723 = x1171 + x1722
    x1724 = x1723 * x3
    x1725 = x1261 * x968
    x1726 = x1182 + x1725
    x1727 = x1726 * x6
    x1728 = -x1727
    x1729 = x4 * (x1176 + x1719 - x1722)
    x1730 = -x1729
    x1731 = x1719 * x3
    x1732 = x1723 * x6
    x1733 = -x1732
    x1734 = x1731 + x1733
    x1735 = x1719 * x5
    x1736 = x1717 * x5 + x1720
    x1737 = x1733 + x1735
    x1738 = x1718 + x1736 * x3 - x1737 * x6
    x1739 = x1737 * x3
    x1740 = x1723 * x5
    x1741 = x1728 + x1740
    x1742 = x1741 * x6
    x1743 = x1729 + x1739 - x1742
    x1744 = x4 * (x1188 + x1723 - x1725)
    x1745 = x1190 + x1261 * x977
    x1746 = x1019 * x1261
    x1747 = x1013 * x1261 + x1195
    x1748 = x4 * (x1193 - x1746 + x1747)
    x1749 = x1198 + x1746
    x1750 = -x1749 * x6
    x1751 = x1747 * x3 + x1750
    x1752 = x1027 * x1261
    x1753 = x1204 + x1752
    x1754 = x1753 * x3
    x1755 = x1039 * x1261
    x1756 = x1214 + x1755
    x1757 = x1756 * x6
    x1758 = -x1757
    x1759 = x4 * (x1203 + x1749 - x1752)
    x1760 = -x1759
    x1761 = x1749 * x3
    x1762 = x1753 * x6
    x1763 = -x1762
    x1764 = x1761 + x1763
    x1765 = x1749 * x5
    x1766 = x1747 * x5 + x1750
    x1767 = x1763 + x1765
    x1768 = x1748 + x1766 * x3 - x1767 * x6
    x1769 = x1767 * x3
    x1770 = x1753 * x5
    x1771 = x1758 + x1770
    x1772 = x1771 * x6
    x1773 = x1759 + x1769 - x1772
    x1774 = x4 * (x1213 + x1753 - x1755)
    x1775 = x1049 * x1261 + x1221
    x1776 = x1100 * x8
    x1777 = x1093 * x1261
    x1778 = 4.0 * x844
    x1779 = x1087 * x1261 - x1093 * x8 + 4.0 * x807
    x1780 = x4 * (x1776 - x1777 - x1778 + x1779)
    x1781 = -x1776 + x1777 + x1778
    x1782 = x1779 * x3 - x1781 * x6
    x1783 = x1781 * x3
    x1784 = x1100 * x1261
    x1785 = x1112 * x8
    x1786 = 4.0 * x826
    x1787 = x1784 - x1785 + x1786
    x1788 = x1787 * x6
    x1789 = x1783 - x1788
    x1790 = x4 * (x1781 - x1784 + x1785 - x1786)
    x1791 = x1112 * x1261 - x1120 * x8 + 4.0 * x872
    x1792 = x1137 * x8
    x1793 = x1130 * x1261
    x1794 = 3.0 * x914
    x1795 = x1124 * x1261 - x1130 * x8 + 3.0 * x880
    x1796 = x4 * (x1792 - x1793 - x1794 + x1795)
    x1797 = -x1792 + x1793 + x1794
    x1798 = x1795 * x3 - x1797 * x6
    x1799 = x1797 * x3
    x1800 = x1137 * x1261
    x1801 = x1149 * x8
    x1802 = 3.0 * x897
    x1803 = x1800 - x1801 + x1802
    x1804 = x1803 * x6
    x1805 = x1799 - x1804
    x1806 = x4 * (x1797 - x1800 + x1801 - x1802)
    x1807 = x1149 * x1261 - x1157 * x8 + 3.0 * x942
    x1808 = x1172 * x8
    x1809 = x1166 * x1261
    x1810 = 2.0 * x981
    x1811 = x1161 * x1261 - x1166 * x8 + 2.0 * x949
    x1812 = x4 * (x1808 - x1809 - x1810 + x1811)
    x1813 = -x1808 + x1809 + x1810
    x1814 = x1811 * x3 - x1813 * x6
    x1815 = x1813 * x3
    x1816 = x1172 * x1261
    x1817 = x1183 * x8
    x1818 = 2.0 * x965
    x1819 = x1816 - x1817 + x1818
    x1820 = x1819 * x6
    x1821 = x1815 - x1820
    x1822 = x4 * (x1813 - x1816 + x1817 - x1818)
    x1823 = 2.0 * x1008 + x1183 * x1261 - x1191 * x8
    x1824 = x1205 * x8
    x1825 = x1199 * x1261
    x1826 = x1017 + x1196 * x1261 - x1199 * x8
    x1827 = x4 * (x1060 + x1824 - x1825 + x1826)
    x1828 = x1054 - x1824 + x1825
    x1829 = x1826 * x3 - x1828 * x6
    x1830 = x1828 * x3
    x1831 = x1205 * x1261
    x1832 = x1215 * x8
    x1833 = x1036 + x1831 - x1832
    x1834 = x1833 * x6
    x1835 = x1830 - x1834
    x1836 = x4 * (x1080 + x1828 - x1831 + x1832)
    x1837 = x1082 + x1215 * x1261 - x1222 * x8
    x1838 = x1239 * x8
    x1839 = x1232 * x1261
    x1840 = x1226 * x1261 - x1232 * x8
    x1841 = x4 * (x1838 - x1839 + x1840)
    x1842 = -x1838 + x1839
    x1843 = x1840 * x3 - x1842 * x6
    x1844 = x1842 * x3
    x1845 = x1239 * x1261
    x1846 = x1251 * x8
    x1847 = x1845 - x1846
    x1848 = x1847 * x6
    x1849 = x1844 - x1848
    x1850 = x4 * (x1842 - x1845 + x1846)
    x1851 = x1251 * x1261 - x1259 * x8
    x1852 = -x9 - A[2]
    x1853 = x13 * x1852
    x1854 = x15 * x1852
    x1855 = x1854 + x320
    x1856 = x4 * (-x1853 + x1855 + x325)
    x1857 = x1855 * x5
    x1858 = x1853 + x326
    x1859 = x1858 * x6
    x1860 = x1857 - x1859
    x1861 = x1860 * x5
    x1862 = x1858 * x5
    x1863 = x1852 * x20
    x1864 = x1863 + x336
    x1865 = x1864 * x6
    x1866 = x1862 - x1865
    x1867 = x1866 * x6
    x1868 = x1856 + x1861 - x1867
    x1869 = x1868 * x5
    x1870 = x4 * (x1858 - x1863 + x335)
    x1871 = x1866 * x5
    x1872 = x1864 * x5
    x1873 = x1852 * x30
    x1874 = x1873 + x360
    x1875 = x1874 * x6
    x1876 = x1872 - x1875
    x1877 = x1876 * x6
    x1878 = x1870 + x1871 - x1877
    x1879 = x1878 * x6
    x1880 = x38 * (x1860 - x1862 + x1865)
    x1881 = -x1879 + x1880
    x1882 = x1869 + x1881
    x1883 = x1882 * x5
    x1884 = x1852 * x43 + x346
    x1885 = x4 * (-x1854 + x1884 + x319)
    x1886 = -x1855 * x6 + x1884 * x5
    x1887 = -x1860 * x6 + x1885 + x1886 * x5
    x1888 = -x1868 * x6 + x38 * (-x1857 + x1859 + x1886)
    x1889 = x1887 * x5 + x1888
    x1890 = -x1856
    x1891 = x130 * (-x1861 + x1867 + x1887 + x1890) - x1882 * x6
    x1892 = x1889 * x5 + x1891
    x1893 = x1878 * x5
    x1894 = x4 * (x1864 - x1873 + x359)
    x1895 = x1876 * x5
    x1896 = x1874 * x5
    x1897 = x1852 * x58
    x1898 = x1897 + x384
    x1899 = x1898 * x6
    x1900 = x1896 - x1899
    x1901 = x1900 * x6
    x1902 = x1894 + x1895 - x1901
    x1903 = x1902 * x6
    x1904 = x38 * (x1866 - x1872 + x1875)
    x1905 = x1893 - x1903 + x1904
    x1906 = x1905 * x6
    x1907 = -x1870
    x1908 = x130 * (x1868 - x1871 + x1877 + x1907)
    x1909 = x1906 - x1908
    x1910 = -x1906 + x1908
    x1911 = x1883 + x1910
    x1912 = x1879 - x1880
    x1913 = x4 * (-x1869 + x1889 + x1912)
    x1914 = x1892 * x3 - x1911 * x6 + 4.0 * x1913
    x1915 = x1911 * x3
    x1916 = x1905 * x5
    x1917 = x1902 * x5
    x1918 = x4 * (x1874 - x1897 + x383)
    x1919 = x1900 * x5
    x1920 = x1898 * x5
    x1921 = x1852 * x84
    x1922 = x1921 + x416
    x1923 = x1922 * x6
    x1924 = x1920 - x1923
    x1925 = x1924 * x6
    x1926 = x1918 + x1919 - x1925
    x1927 = x1926 * x6
    x1928 = x38 * (x1876 - x1896 + x1899)
    x1929 = x1917 - x1927 + x1928
    x1930 = x1929 * x6
    x1931 = -x1894
    x1932 = x130 * (x1878 - x1895 + x1901 + x1931)
    x1933 = -x1930 + x1932
    x1934 = x1916 + x1933
    x1935 = x1934 * x6
    x1936 = x1903 - x1904
    x1937 = x4 * (x1882 - x1893 + x1936)
    x1938 = 4.0 * x1937
    x1939 = x1915 - x1935 + x1938
    x1940 = x1882 * x3
    x1941 = x1889 * x3 + x1891
    x1942 = x1930 - x1932
    x1943 = x4 * (x1898 - x1921 + x415)
    x1944 = x119 * x1852 + x450
    x1945 = -x1918
    x1946 = x1905 * x3
    x1947 = x1910 + x1940
    x1948 = x1868 * x3
    x1949 = x174 * x1855
    x1950 = x1858 * x8
    x1951 = -x1950
    x1952 = x1949 + x1951
    x1953 = x1952 * x5
    x1954 = x174 * x1858
    x1955 = x1864 * x8
    x1956 = -x1955
    x1957 = x1954 + x1956
    x1958 = x1957 * x6
    x1959 = x1953 - x1958
    x1960 = x1959 * x5
    x1961 = x4 * (x1952 - x1954 + x1955)
    x1962 = x1957 * x5
    x1963 = x174 * x1864
    x1964 = x1874 * x8
    x1965 = -x1964
    x1966 = x1963 + x1965
    x1967 = x1966 * x6
    x1968 = x1962 - x1967
    x1969 = x1968 * x6
    x1970 = x1961 - x1969
    x1971 = x1960 + x1970
    x1972 = x1971 * x5
    x1973 = -x1855 * x8
    x1974 = x174 * x1884 + x1973
    x1975 = -x1952 * x6 + x1974 * x5
    x1976 = x4 * (-x1949 + x1950 + x1974)
    x1977 = -x1959 * x6 + x1976
    x1978 = x1975 * x5 + x1977
    x1979 = -x1971 * x6 + x38 * (-x1953 + x1958 + x1975)
    x1980 = x1978 * x5 + x1979
    x1981 = x4 * (x1957 - x1963 + x1964)
    x1982 = x1968 * x5
    x1983 = x1966 * x5
    x1984 = x174 * x1874
    x1985 = x1898 * x8
    x1986 = -x1985
    x1987 = x1984 + x1986
    x1988 = x1987 * x6
    x1989 = x1983 - x1988
    x1990 = x1989 * x6
    x1991 = x1981 + x1982 - x1990
    x1992 = x1991 * x6
    x1993 = x38 * (x1959 - x1962 + x1967)
    x1994 = x1992 - x1993
    x1995 = -x1992 + x1993
    x1996 = x1972 + x1995
    x1997 = -x1961
    x1998 = x1969 + x1997
    x1999 = x4 * (-x1960 + x1978 + x1998)
    x2000 = x1980 * x3 - x1996 * x6 + 3.0 * x1999
    x2001 = x1996 * x3
    x2002 = x1991 * x5
    x2003 = x4 * (x1966 - x1984 + x1985)
    x2004 = x1989 * x5
    x2005 = x1987 * x5
    x2006 = x174 * x1898
    x2007 = x1922 * x8
    x2008 = -x2007
    x2009 = x2006 + x2008
    x2010 = x2009 * x6
    x2011 = x2005 - x2010
    x2012 = x2011 * x6
    x2013 = x2003 + x2004 - x2012
    x2014 = x2013 * x6
    x2015 = x38 * (x1968 - x1983 + x1988)
    x2016 = -x2014 + x2015
    x2017 = x2002 + x2016
    x2018 = x2017 * x6
    x2019 = -x1981
    x2020 = x1990 + x2019
    x2021 = x4 * (x1971 - x1982 + x2020)
    x2022 = 3.0 * x2021
    x2023 = x2001 - x2018 + x2022
    x2024 = x1971 * x3
    x2025 = x1978 * x3 + x1979
    x2026 = x2014 - x2015
    x2027 = x4 * (x1987 - x2006 + x2007)
    x2028 = -x1944 * x8
    x2029 = x174 * x1922 + x2028
    x2030 = -x2003
    x2031 = x1991 * x3
    x2032 = x1995 + x2024
    x2033 = x1959 * x3
    x2034 = x1852 * x323
    x2035 = x2034 + x690
    x2036 = x2035 * x5
    x2037 = x1852 * x328
    x2038 = x2037 + x695
    x2039 = x2038 * x6
    x2040 = x2036 - x2039
    x2041 = x2040 * x5
    x2042 = x4 * (x2035 - x2037 + x718)
    x2043 = x2038 * x5
    x2044 = x1852 * x338
    x2045 = x2044 + x712
    x2046 = x2045 * x6
    x2047 = x2043 - x2046
    x2048 = x2047 * x6
    x2049 = x2042 - x2048
    x2050 = x2041 + x2049
    x2051 = x2050 * x5
    x2052 = x1852 * x347 + x701
    x2053 = -x2035 * x6 + x2052 * x5
    x2054 = x4 * (-x2034 + x2052 + x705)
    x2055 = -x2040 * x6 + x2054
    x2056 = x2053 * x5 + x2055
    x2057 = -x2050 * x6 + x38 * (-x2036 + x2039 + x2053)
    x2058 = x2056 * x5 + x2057
    x2059 = x4 * (x2038 - x2044 + x728)
    x2060 = x2047 * x5
    x2061 = x2045 * x5
    x2062 = x1852 * x362
    x2063 = x2062 + x733
    x2064 = x2063 * x6
    x2065 = x2061 - x2064
    x2066 = x2065 * x6
    x2067 = x2059 + x2060 - x2066
    x2068 = x2067 * x6
    x2069 = x38 * (x2040 - x2043 + x2046)
    x2070 = x2068 - x2069
    x2071 = -x2068 + x2069
    x2072 = x2051 + x2071
    x2073 = -x2042
    x2074 = x2048 + x2073
    x2075 = x4 * (-x2041 + x2056 + x2074)
    x2076 = x2058 * x3 - x2072 * x6 + 3.0 * x2075
    x2077 = x2072 * x3
    x2078 = x2067 * x5
    x2079 = x4 * (x2045 - x2062 + x755)
    x2080 = x2065 * x5
    x2081 = x2063 * x5
    x2082 = x1852 * x386
    x2083 = x2082 + x760
    x2084 = x2083 * x6
    x2085 = x2081 - x2084
    x2086 = x2085 * x6
    x2087 = x2079 + x2080 - x2086
    x2088 = x2087 * x6
    x2089 = x38 * (x2047 - x2061 + x2064)
    x2090 = -x2088 + x2089
    x2091 = x2078 + x2090
    x2092 = x2091 * x6
    x2093 = -x2059
    x2094 = x2066 + x2093
    x2095 = x4 * (x2050 - x2060 + x2094)
    x2096 = 3.0 * x2095
    x2097 = x2077 - x2092 + x2096
    x2098 = x2050 * x3
    x2099 = x2056 * x3 + x2057
    x2100 = x2088 - x2089
    x2101 = x4 * (x2063 - x2082 + x787)
    x2102 = x1852 * x418 + x789
    x2103 = -x2079
    x2104 = x2067 * x3
    x2105 = x2071 + x2098
    x2106 = x2040 * x3
    x2107 = x174 * x1952
    x2108 = x1957 * x8
    x2109 = x1856 - x2108
    x2110 = x2107 + x2109
    x2111 = x2110 * x5
    x2112 = x174 * x1957
    x2113 = x1966 * x8
    x2114 = x1870 - x2113
    x2115 = x2112 + x2114
    x2116 = x2115 * x6
    x2117 = -x2116
    x2118 = x2111 + x2117
    x2119 = x2118 * x5
    x2120 = x1885 - x1952 * x8
    x2121 = x174 * x1974 + x2120
    x2122 = -x2110 * x6
    x2123 = x2121 * x5 + x2122
    x2124 = x1890 + x2108
    x2125 = x4 * (-x2107 + x2121 + x2124)
    x2126 = -x2118 * x6 + x2125
    x2127 = x2123 * x5 + x2126
    x2128 = x2115 * x5
    x2129 = x174 * x1966
    x2130 = x1987 * x8
    x2131 = x1894 - x2130
    x2132 = x2129 + x2131
    x2133 = x2132 * x6
    x2134 = x2128 - x2133
    x2135 = x2134 * x6
    x2136 = x1907 + x2113
    x2137 = x4 * (x2110 - x2112 + x2136)
    x2138 = -x2137
    x2139 = x2135 + x2138
    x2140 = -x2135 + x2137
    x2141 = x2119 + x2140
    x2142 = x4 * (-x2111 + x2116 + x2123)
    x2143 = x2127 * x3 - x2141 * x6 + 2.0 * x2142
    x2144 = x2141 * x3
    x2145 = x2134 * x5
    x2146 = x1931 + x2130
    x2147 = x4 * (x2115 - x2129 + x2146)
    x2148 = x2132 * x5
    x2149 = x174 * x1987
    x2150 = x2009 * x8
    x2151 = x1918 - x2150
    x2152 = x2149 + x2151
    x2153 = x2152 * x6
    x2154 = x2148 - x2153
    x2155 = x2154 * x6
    x2156 = x2147 - x2155
    x2157 = x2145 + x2156
    x2158 = x2157 * x6
    x2159 = x4 * (x2118 - x2128 + x2133)
    x2160 = 2.0 * x2159
    x2161 = x2144 - x2158 + x2160
    x2162 = x2118 * x3
    x2163 = x2123 * x3 + x2126
    x2164 = -x2147
    x2165 = x2155 + x2164
    x2166 = x1945 + x2150
    x2167 = x4 * (x2132 - x2149 + x2166)
    x2168 = x1943 - x2029 * x8
    x2169 = x174 * x2009 + x2168
    x2170 = x2134 * x3
    x2171 = x2140 + x2162
    x2172 = x2110 * x3
    x2173 = x174 * x2035
    x2174 = x2038 * x8
    x2175 = -x2174
    x2176 = x2173 + x2175
    x2177 = x2176 * x5
    x2178 = x174 * x2038
    x2179 = x2045 * x8
    x2180 = -x2179
    x2181 = x2178 + x2180
    x2182 = x2181 * x6
    x2183 = -x2182
    x2184 = x2177 + x2183
    x2185 = x2184 * x5
    x2186 = -x2035 * x8
    x2187 = x174 * x2052 + x2186
    x2188 = -x2176 * x6
    x2189 = x2187 * x5 + x2188
    x2190 = x4 * (-x2173 + x2174 + x2187)
    x2191 = -x2184 * x6 + x2190
    x2192 = x2189 * x5 + x2191
    x2193 = x2181 * x5
    x2194 = x174 * x2045
    x2195 = x2063 * x8
    x2196 = -x2195
    x2197 = x2194 + x2196
    x2198 = x2197 * x6
    x2199 = x2193 - x2198
    x2200 = x2199 * x6
    x2201 = x4 * (x2176 - x2178 + x2179)
    x2202 = -x2201
    x2203 = x2200 + x2202
    x2204 = -x2200 + x2201
    x2205 = x2185 + x2204
    x2206 = x4 * (-x2177 + x2182 + x2189)
    x2207 = x2192 * x3 - x2205 * x6 + 2.0 * x2206
    x2208 = x2205 * x3
    x2209 = x2199 * x5
    x2210 = x4 * (x2181 - x2194 + x2195)
    x2211 = x2197 * x5
    x2212 = x174 * x2063
    x2213 = x2083 * x8
    x2214 = -x2213
    x2215 = x2212 + x2214
    x2216 = x2215 * x6
    x2217 = x2211 - x2216
    x2218 = x2217 * x6
    x2219 = x2210 - x2218
    x2220 = x2209 + x2219
    x2221 = x2220 * x6
    x2222 = x4 * (x2184 - x2193 + x2198)
    x2223 = 2.0 * x2222
    x2224 = x2208 - x2221 + x2223
    x2225 = x2184 * x3
    x2226 = x2189 * x3 + x2191
    x2227 = -x2210
    x2228 = x2218 + x2227
    x2229 = x4 * (x2197 - x2212 + x2213)
    x2230 = -x2102 * x8
    x2231 = x174 * x2083 + x2230
    x2232 = x2199 * x3
    x2233 = x2204 + x2225
    x2234 = x2176 * x3
    x2235 = x1852 * x691
    x2236 = x1018 + x2235
    x2237 = x2236 * x5
    x2238 = x1852 * x696
    x2239 = x1026 + x2238
    x2240 = x2239 * x6
    x2241 = -x2240
    x2242 = x2237 + x2241
    x2243 = x2242 * x5
    x2244 = x1012 + x1852 * x702
    x2245 = -x2236 * x6
    x2246 = x2244 * x5 + x2245
    x2247 = x4 * (x1016 - x2235 + x2244)
    x2248 = -x2242 * x6 + x2247
    x2249 = x2246 * x5 + x2248
    x2250 = x2239 * x5
    x2251 = x1852 * x713
    x2252 = x1038 + x2251
    x2253 = x2252 * x6
    x2254 = x2250 - x2253
    x2255 = x2254 * x6
    x2256 = x4 * (x1053 + x2236 - x2238)
    x2257 = -x2256
    x2258 = x2255 + x2257
    x2259 = -x2255 + x2256
    x2260 = x2243 + x2259
    x2261 = x4 * (-x2237 + x2240 + x2246)
    x2262 = x2249 * x3 - x2260 * x6 + 2.0 * x2261
    x2263 = x2260 * x3
    x2264 = x2254 * x5
    x2265 = x4 * (x1035 + x2239 - x2251)
    x2266 = x2252 * x5
    x2267 = x1852 * x734
    x2268 = x1048 + x2267
    x2269 = x2268 * x6
    x2270 = x2266 - x2269
    x2271 = x2270 * x6
    x2272 = x2265 - x2271
    x2273 = x2264 + x2272
    x2274 = x2273 * x6
    x2275 = x4 * (x2242 - x2250 + x2253)
    x2276 = 2.0 * x2275
    x2277 = x2263 - x2274 + x2276
    x2278 = x2242 * x3
    x2279 = x2246 * x3 + x2248
    x2280 = -x2265
    x2281 = x2271 + x2280
    x2282 = x4 * (x1081 + x2252 - x2267)
    x2283 = x1083 + x1852 * x761
    x2284 = x2254 * x3
    x2285 = x2259 + x2278
    x2286 = x2236 * x3
    x2287 = x174 * x2110
    x2288 = 2.0 * x1976 - x2110 * x8
    x2289 = x174 * x2121 + x2288
    x2290 = x2115 * x8
    x2291 = 2.0 * x1961
    x2292 = x2290 - x2291
    x2293 = x4 * (-x2287 + x2289 + x2292)
    x2294 = -x2290 + x2291
    x2295 = x2287 + x2294
    x2296 = -x2295 * x6
    x2297 = x2289 * x3 + x2296
    x2298 = x174 * x2115
    x2299 = x2132 * x8
    x2300 = 2.0 * x1981
    x2301 = -x2299 + x2300
    x2302 = x2298 + x2301
    x2303 = x2302 * x3
    x2304 = x174 * x2132
    x2305 = x2152 * x8
    x2306 = 2.0 * x2003
    x2307 = -x2305 + x2306
    x2308 = x2304 + x2307
    x2309 = x2308 * x6
    x2310 = -x2309
    x2311 = x2299 - x2300
    x2312 = x4 * (x2295 - x2298 + x2311)
    x2313 = -x2312
    x2314 = x2295 * x3
    x2315 = x2302 * x6
    x2316 = -x2315
    x2317 = x2314 + x2316
    x2318 = x2295 * x5
    x2319 = x2289 * x5 + x2296
    x2320 = x2316 + x2318
    x2321 = x2293 + x2319 * x3 - x2320 * x6
    x2322 = x2320 * x3
    x2323 = x2302 * x5
    x2324 = x2310 + x2323
    x2325 = x2324 * x6
    x2326 = x2312 + x2322 - x2325
    x2327 = x2305 - x2306
    x2328 = x4 * (x2302 - x2304 + x2327)
    x2329 = 2.0 * x2027 - x2169 * x8
    x2330 = x174 * x2152 + x2329
    x2331 = x174 * x2176
    x2332 = x2054 - x2176 * x8
    x2333 = x174 * x2187 + x2332
    x2334 = x2181 * x8
    x2335 = x2073 + x2334
    x2336 = x4 * (-x2331 + x2333 + x2335)
    x2337 = x2042 - x2334
    x2338 = x2331 + x2337
    x2339 = -x2338 * x6
    x2340 = x2333 * x3 + x2339
    x2341 = x174 * x2181
    x2342 = x2197 * x8
    x2343 = x2059 - x2342
    x2344 = x2341 + x2343
    x2345 = x2344 * x3
    x2346 = x174 * x2197
    x2347 = x2215 * x8
    x2348 = x2079 - x2347
    x2349 = x2346 + x2348
    x2350 = x2349 * x6
    x2351 = -x2350
    x2352 = x2093 + x2342
    x2353 = x4 * (x2338 - x2341 + x2352)
    x2354 = -x2353
    x2355 = x2338 * x3
    x2356 = x2344 * x6
    x2357 = -x2356
    x2358 = x2355 + x2357
    x2359 = x2338 * x5
    x2360 = x2333 * x5 + x2339
    x2361 = x2357 + x2359
    x2362 = x2336 + x2360 * x3 - x2361 * x6
    x2363 = x2361 * x3
    x2364 = x2344 * x5
    x2365 = x2351 + x2364
    x2366 = x2365 * x6
    x2367 = x2353 + x2363 - x2366
    x2368 = x2103 + x2347
    x2369 = x4 * (x2344 - x2346 + x2368)
    x2370 = x2101 - x2231 * x8
    x2371 = x174 * x2215 + x2370
    x2372 = x2239 * x8
    x2373 = x174 * x2236
    x2374 = -x2236 * x8
    x2375 = x174 * x2244 + x2374
    x2376 = x4 * (x2372 - x2373 + x2375)
    x2377 = -x2372
    x2378 = x2373 + x2377
    x2379 = -x2378 * x6
    x2380 = x2375 * x3 + x2379
    x2381 = x174 * x2239
    x2382 = x2252 * x8
    x2383 = -x2382
    x2384 = x2381 + x2383
    x2385 = x2384 * x3
    x2386 = x174 * x2252
    x2387 = x2268 * x8
    x2388 = -x2387
    x2389 = x2386 + x2388
    x2390 = x2389 * x6
    x2391 = -x2390
    x2392 = x4 * (x2378 - x2381 + x2382)
    x2393 = -x2392
    x2394 = x2378 * x3
    x2395 = x2384 * x6
    x2396 = -x2395
    x2397 = x2394 + x2396
    x2398 = x2378 * x5
    x2399 = x2375 * x5 + x2379
    x2400 = x2396 + x2398
    x2401 = x2376 + x2399 * x3 - x2400 * x6
    x2402 = x2400 * x3
    x2403 = x2384 * x5
    x2404 = x2391 + x2403
    x2405 = x2404 * x6
    x2406 = x2392 + x2402 - x2405
    x2407 = x4 * (x2384 - x2386 + x2387)
    x2408 = -x2283 * x8
    x2409 = x174 * x2268 + x2408
    x2410 = x1019 * x1852
    x2411 = x1013 * x1852 + x1225
    x2412 = x4 * (x1229 - x2410 + x2411)
    x2413 = x1231 + x2410
    x2414 = -x2413 * x6
    x2415 = x2411 * x3 + x2414
    x2416 = x1027 * x1852
    x2417 = x1238 + x2416
    x2418 = x2417 * x3
    x2419 = x1039 * x1852
    x2420 = x1250 + x2419
    x2421 = x2420 * x6
    x2422 = -x2421
    x2423 = x4 * (x1243 + x2413 - x2416)
    x2424 = -x2423
    x2425 = x2413 * x3
    x2426 = x2417 * x6
    x2427 = -x2426
    x2428 = x2425 + x2427
    x2429 = x2413 * x5
    x2430 = x2411 * x5 + x2414
    x2431 = x2427 + x2429
    x2432 = x2412 + x2430 * x3 - x2431 * x6
    x2433 = x2431 * x3
    x2434 = x2417 * x5
    x2435 = x2422 + x2434
    x2436 = x2435 * x6
    x2437 = x2423 + x2433 - x2436
    x2438 = x4 * (x1256 + x2417 - x2419)
    x2439 = x1049 * x1852 + x1258
    x2440 = x174 * x2295
    x2441 = 3.0 * x2125 - x2295 * x8
    x2442 = x174 * x2289 + x2441
    x2443 = x2302 * x8
    x2444 = 3.0 * x2137
    x2445 = x2443 - x2444
    x2446 = x4 * (-x2440 + x2442 + x2445)
    x2447 = -x2443 + x2444
    x2448 = x2440 + x2447
    x2449 = x2442 * x3 - x2448 * x6
    x2450 = x2448 * x3
    x2451 = x174 * x2302
    x2452 = x2308 * x8
    x2453 = 3.0 * x2147
    x2454 = -x2452 + x2453
    x2455 = x2451 + x2454
    x2456 = x2455 * x6
    x2457 = x2450 - x2456
    x2458 = x2452 - x2453
    x2459 = x4 * (x2448 - x2451 + x2458)
    x2460 = 3.0 * x2167 - x2330 * x8
    x2461 = x174 * x2308 + x2460
    x2462 = x174 * x2338
    x2463 = 2.0 * x2190 - x2338 * x8
    x2464 = x174 * x2333 + x2463
    x2465 = x2344 * x8
    x2466 = 2.0 * x2201
    x2467 = x2465 - x2466
    x2468 = x4 * (-x2462 + x2464 + x2467)
    x2469 = -x2465 + x2466
    x2470 = x2462 + x2469
    x2471 = x2464 * x3 - x2470 * x6
    x2472 = x2470 * x3
    x2473 = x174 * x2344
    x2474 = x2349 * x8
    x2475 = 2.0 * x2210
    x2476 = -x2474 + x2475
    x2477 = x2473 + x2476
    x2478 = x2477 * x6
    x2479 = x2472 - x2478
    x2480 = x2474 - x2475
    x2481 = x4 * (x2470 - x2473 + x2480)
    x2482 = 2.0 * x2229 - x2371 * x8
    x2483 = x174 * x2349 + x2482
    x2484 = x174 * x2378
    x2485 = x2247 - x2378 * x8
    x2486 = x174 * x2375 + x2485
    x2487 = x2384 * x8
    x2488 = x2257 + x2487
    x2489 = x4 * (-x2484 + x2486 + x2488)
    x2490 = x2256 - x2487
    x2491 = x2484 + x2490
    x2492 = x2486 * x3 - x2491 * x6
    x2493 = x2491 * x3
    x2494 = x174 * x2384
    x2495 = x2389 * x8
    x2496 = x2265 - x2495
    x2497 = x2494 + x2496
    x2498 = x2497 * x6
    x2499 = x2493 - x2498
    x2500 = x2280 + x2495
    x2501 = x4 * (x2491 - x2494 + x2500)
    x2502 = x2282 - x2409 * x8
    x2503 = x174 * x2389 + x2502
    x2504 = x2417 * x8
    x2505 = x174 * x2413
    x2506 = -x2413 * x8
    x2507 = x174 * x2411 + x2506
    x2508 = x4 * (x2504 - x2505 + x2507)
    x2509 = -x2504
    x2510 = x2505 + x2509
    x2511 = x2507 * x3 - x2510 * x6
    x2512 = x2510 * x3
    x2513 = x174 * x2417
    x2514 = x2420 * x8
    x2515 = -x2514
    x2516 = x2513 + x2515
    x2517 = x2516 * x6
    x2518 = x2512 - x2517
    x2519 = x4 * (x2510 - x2513 + x2514)
    x2520 = -x2439 * x8
    x2521 = x174 * x2420 + x2520
    x2522 = x10 * x1239
    x2523 = x1232 * x1852
    x2524 = 4.0 * x1054
    x2525 = -x10 * x1232 + 4.0 * x1017 + x1226 * x1852
    x2526 = x4 * (x2522 - x2523 - x2524 + x2525)
    x2527 = -x2522 + x2523 + x2524
    x2528 = x2525 * x3 - x2527 * x6
    x2529 = x2527 * x3
    x2530 = x1239 * x1852
    x2531 = x10 * x1251
    x2532 = 4.0 * x1036
    x2533 = x2530 - x2531 + x2532
    x2534 = x2533 * x6
    x2535 = x2529 - x2534
    x2536 = x4 * (x2527 - x2530 + x2531 - x2532)
    x2537 = -x10 * x1259 + 4.0 * x1082 + x1251 * x1852
    x2538 = x1273 * x8
    x2539 = x1261 * x1267
    x2540 = x1261 * x1264
    x2541 = x1267 * x8
    x2542 = x19 + x2540 - x2541
    x2543 = x4 * (x2538 - x2539 + x2542 + x69)
    x2544 = x2542 * x5
    x2545 = -x2538 + x2539 + x29
    x2546 = x2545 * x6
    x2547 = x2544 - x2546
    x2548 = x2547 * x5
    x2549 = x2545 * x5
    x2550 = x1261 * x1273
    x2551 = x1283 * x8
    x2552 = x2550 - x2551 + x56
    x2553 = x2552 * x6
    x2554 = x2549 - x2553
    x2555 = x2554 * x6
    x2556 = x2543 + x2548 - x2555
    x2557 = x2556 * x5
    x2558 = x4 * (x2545 - x2550 + x2551 + x95)
    x2559 = x2554 * x5
    x2560 = x2552 * x5
    x2561 = x1261 * x1283
    x2562 = x1307 * x8
    x2563 = x2561 - x2562 + x82
    x2564 = x2563 * x6
    x2565 = x2560 - x2564
    x2566 = x2565 * x6
    x2567 = x2558 + x2559 - x2566
    x2568 = x2567 * x6
    x2569 = x38 * (x2547 - x2549 + x2553)
    x2570 = x2557 - x2568 + x2569
    x2571 = x2570 * x5
    x2572 = x1261 * x1293 - x1264 * x8 + x45
    x2573 = x4 * (-x2540 + x2541 + x2572 + x50)
    x2574 = -x2542 * x6 + x2572 * x5
    x2575 = -x2547 * x6 + x2573 + x2574 * x5
    x2576 = -x2556 * x6 + x2575 * x5 + x38 * (-x2544 + x2546 + x2574)
    x2577 = -x2543
    x2578 = x130 * (-x2548 + x2555 + x2575 + x2577) - x2570 * x6
    x2579 = x2576 * x5 + x2578
    x2580 = x2567 * x5
    x2581 = x4 * (x129 + x2552 - x2561 + x2562)
    x2582 = x2565 * x5
    x2583 = x2563 * x5
    x2584 = x1261 * x1307
    x2585 = x1331 * x8
    x2586 = x117 + x2584 - x2585
    x2587 = x2586 * x6
    x2588 = x2583 - x2587
    x2589 = x2588 * x6
    x2590 = x2581 + x2582 - x2589
    x2591 = x2590 * x6
    x2592 = x38 * (x2554 - x2560 + x2564)
    x2593 = x2580 - x2591 + x2592
    x2594 = x2593 * x6
    x2595 = -x2558
    x2596 = x130 * (x2556 - x2559 + x2566 + x2595)
    x2597 = x2594 - x2596
    x2598 = x2571 - x2594 + x2596
    x2599 = x4 * (x157 + x2563 - x2584 + x2585)
    x2600 = x1261 * x1331 - x1353 * x8 + x155
    x2601 = -x2581
    x2602 = x1369 * x8
    x2603 = x1261 * x1362
    x2604 = x1261 * x1359
    x2605 = x1362 * x8
    x2606 = x1265 + x186 + x2604 - x2605
    x2607 = x4 * (x1316 + x251 + x2602 - x2603 + x2606)
    x2608 = x2606 * x5
    x2609 = x1279 + x210 - x2602 + x2603
    x2610 = x2609 * x6
    x2611 = x2608 - x2610
    x2612 = x2611 * x5
    x2613 = x2609 * x5
    x2614 = x1261 * x1369
    x2615 = x1387 * x8
    x2616 = x1303 + x234 + x2614 - x2615
    x2617 = x2616 * x6
    x2618 = x2613 - x2617
    x2619 = x2618 * x6
    x2620 = x2607 + x2612 - x2619
    x2621 = x2620 * x5
    x2622 = x1261 * x1376 + x1294 - x1359 * x8 + x203
    x2623 = x4 * (x1299 + x227 - x2604 + x2605 + x2622)
    x2624 = -x2606 * x6 + x2622 * x5
    x2625 = -x2611 * x6 + x2623 + x2624 * x5
    x2626 = -x2620 * x6 + x38 * (-x2608 + x2610 + x2624)
    x2627 = x2625 * x5 + x2626
    x2628 = x4 * (x1340 + x2609 - x2614 + x2615 + x281)
    x2629 = x2618 * x5
    x2630 = x2616 * x5
    x2631 = x1261 * x1387
    x2632 = x1406 * x8
    x2633 = x1327 + x2631 - x2632 + x265
    x2634 = x2633 * x6
    x2635 = x2630 - x2634
    x2636 = x2635 * x6
    x2637 = x2628 + x2629 - x2636
    x2638 = x2637 * x6
    x2639 = x38 * (x2611 - x2613 + x2617)
    x2640 = x2638 - x2639
    x2641 = x2621 - x2638 + x2639
    x2642 = x4 * (x1354 + x2616 - x2631 + x2632 + x305)
    x2643 = x1261 * x1406 + x1352 - x1424 * x8 + x302
    x2644 = x1440 * x8
    x2645 = x1261 * x1433
    x2646 = x1261 * x1430
    x2647 = x1433 * x8
    x2648 = x2646 - x2647 + x333
    x2649 = x4 * (x2644 - x2645 + x2648 + x398)
    x2650 = x2648 * x5
    x2651 = -x2644 + x2645 + x357
    x2652 = x2651 * x6
    x2653 = x2650 - x2652
    x2654 = x2653 * x5
    x2655 = x2651 * x5
    x2656 = x1261 * x1440
    x2657 = x1458 * x8
    x2658 = x2656 - x2657 + x381
    x2659 = x2658 * x6
    x2660 = x2655 - x2659
    x2661 = x2660 * x6
    x2662 = x2649 + x2654 - x2661
    x2663 = x2662 * x5
    x2664 = x1261 * x1447 - x1430 * x8 + x350
    x2665 = x4 * (-x2646 + x2647 + x2664 + x374)
    x2666 = -x2648 * x6 + x2664 * x5
    x2667 = -x2653 * x6 + x2665 + x2666 * x5
    x2668 = -x2662 * x6 + x38 * (-x2650 + x2652 + x2666)
    x2669 = x2667 * x5 + x2668
    x2670 = x4 * (x2651 - x2656 + x2657 + x428)
    x2671 = x2660 * x5
    x2672 = x2658 * x5
    x2673 = x1261 * x1458
    x2674 = x1478 * x8
    x2675 = x2673 - x2674 + x412
    x2676 = x2675 * x6
    x2677 = x2672 - x2676
    x2678 = x2677 * x6
    x2679 = x2670 + x2671 - x2678
    x2680 = x2679 * x6
    x2681 = x38 * (x2653 - x2655 + x2659)
    x2682 = x2680 - x2681
    x2683 = x2663 - x2680 + x2681
    x2684 = -x2649
    x2685 = x4 * (x2658 - x2673 + x2674 + x452)
    x2686 = x1261 * x1478 - x1497 * x8 + x449
    x2687 = -x2670
    x2688 = x1261 * x1503
    x2689 = x1506 * x8
    x2690 = 2.0 * x1366
    x2691 = x2688 - x2689 + x2690 + x496
    x2692 = x2691 * x5
    x2693 = x1261 * x1506
    x2694 = x1519 * x8
    x2695 = 2.0 * x1383
    x2696 = x2693 - x2694 + x2695 + x506
    x2697 = x2696 * x6
    x2698 = x2692 - x2697
    x2699 = x2698 * x5
    x2700 = 2.0 * x1378
    x2701 = x1261 * x1511 - x1503 * x8 + x2700 + x483
    x2702 = -x2691 * x6 + x2701 * x5
    x2703 = -x2690
    x2704 = x4 * (-x2688 + x2689 + x2701 + x2703 + x497)
    x2705 = -x2698 * x6 + x2704
    x2706 = x2702 * x5 + x2705
    x2707 = x2696 * x5
    x2708 = x1261 * x1519
    x2709 = x1534 * x8
    x2710 = 2.0 * x1402
    x2711 = x2708 - x2709 + x2710 + x533
    x2712 = x2711 * x6
    x2713 = x2707 - x2712
    x2714 = x2713 * x6
    x2715 = -x2695
    x2716 = x4 * (x2691 - x2693 + x2694 + x2715 + x526)
    x2717 = x2714 - x2716
    x2718 = x2699 - x2714 + x2716
    x2719 = -x2710
    x2720 = x4 * (x2696 - x2708 + x2709 + x2719 + x562)
    x2721 = 2.0 * x1423
    x2722 = x1261 * x1534 - x1548 * x8 + x2721 + x565
    x2723 = 80.4530382150027 * x170
    x2724 = x1261 * x1554
    x2725 = x1557 * x8
    x2726 = x1437 + x2724 - x2725 + x608
    x2727 = x2726 * x5
    x2728 = x1261 * x1557
    x2729 = x1570 * x8
    x2730 = x1454 + x2728 - x2729 + x617
    x2731 = x2730 * x6
    x2732 = x2727 - x2731
    x2733 = x2732 * x5
    x2734 = x1261 * x1562 + x1449 - x1554 * x8 + x596
    x2735 = -x2726 * x6 + x2734 * x5
    x2736 = x4 * (x1468 - x2724 + x2725 + x2734 + x609)
    x2737 = -x2732 * x6 + x2736
    x2738 = x2735 * x5 + x2737
    x2739 = x2730 * x5
    x2740 = x1261 * x1570
    x2741 = x1585 * x8
    x2742 = x1474 + x2740 - x2741 + x643
    x2743 = x2742 * x6
    x2744 = x2739 - x2743
    x2745 = x2744 * x6
    x2746 = x4 * (x1488 + x2726 - x2728 + x2729 + x637)
    x2747 = x2745 - x2746
    x2748 = x2733 - x2745 + x2746
    x2749 = x4 * (x1498 + x2730 - x2740 + x2741 + x672)
    x2750 = x1261 * x1585 + x1496 - x1599 * x8 + x674
    x2751 = 139.348749811665 * x170
    x2752 = x1261 * x1605
    x2753 = x1608 * x8
    x2754 = x2752 - x2753 + x719
    x2755 = x2754 * x5
    x2756 = x1261 * x1608
    x2757 = x1621 * x8
    x2758 = x2756 - x2757 + x729
    x2759 = x2758 * x6
    x2760 = x2755 - x2759
    x2761 = x2760 * x5
    x2762 = x1261 * x1613 - x1605 * x8 + x706
    x2763 = -x2754 * x6 + x2762 * x5
    x2764 = x4 * (-x2752 + x2753 + x2762 + x720)
    x2765 = -x2760 * x6 + x2764
    x2766 = x2763 * x5 + x2765
    x2767 = x2758 * x5
    x2768 = x1261 * x1621
    x2769 = x1637 * x8
    x2770 = x2768 - x2769 + x756
    x2771 = x2770 * x6
    x2772 = x2767 - x2771
    x2773 = x2772 * x6
    x2774 = x4 * (x2754 - x2756 + x2757 + x749)
    x2775 = -x2774
    x2776 = x2773 + x2775
    x2777 = x2761 - x2773 + x2774
    x2778 = x4 * (x2758 - x2768 + x2769 + x785)
    x2779 = x1261 * x1637 - x1652 * x8 + x788
    x2780 = x1261 * x1657 + 3.0 * x1514 - x1659 * x8 + x807
    x2781 = x1261 * x1663
    x2782 = x1666 * x8
    x2783 = 3.0 * x1531
    x2784 = x2781 - x2782 + x2783 + x826
    x2785 = x2784 * x6
    x2786 = x1261 * x1659
    x2787 = x1663 * x8
    x2788 = 3.0 * x1523
    x2789 = x2786 - x2787 + x2788 + x844
    x2790 = -x2789 * x6
    x2791 = x2789 * x5
    x2792 = x2780 * x5 + x2790
    x2793 = x4 * (x2780 - x2786 + x2787 - x2788 + x850)
    x2794 = -x2785 + x2791
    x2795 = x4 * (-x2781 + x2782 - x2783 + x2789 + x870)
    x2796 = x1261 * x1666 + 3.0 * x1547 - x1685 * x8 + x872
    x2797 = 2.0 * x1565
    x2798 = x1261 * x1687 - x1689 * x8 + x2797 + x880
    x2799 = x1261 * x1693
    x2800 = x1696 * x8
    x2801 = 2.0 * x1582
    x2802 = x2799 - x2800 + x2801 + x897
    x2803 = x2802 * x6
    x2804 = x1261 * x1689
    x2805 = x1693 * x8
    x2806 = 2.0 * x1574
    x2807 = x2804 - x2805 + x2806 + x914
    x2808 = -x2807 * x6
    x2809 = x2807 * x5
    x2810 = x2798 * x5 + x2808
    x2811 = -x2806
    x2812 = x4 * (x2798 - x2804 + x2805 + x2811 + x920)
    x2813 = -x2803 + x2809
    x2814 = -x2801
    x2815 = x4 * (-x2799 + x2800 + x2807 + x2814 + x940)
    x2816 = 2.0 * x1598
    x2817 = x1261 * x1696 - x1715 * x8 + x2816 + x942
    x2818 = x1261 * x1717 + x1616 - x1719 * x8 + x949
    x2819 = x1261 * x1723
    x2820 = x1726 * x8
    x2821 = x1634 + x2819 - x2820 + x965
    x2822 = x2821 * x6
    x2823 = x1261 * x1719
    x2824 = x1723 * x8
    x2825 = x1625 + x2823 - x2824 + x981
    x2826 = -x2825 * x6
    x2827 = x2825 * x5
    x2828 = x2818 * x5 + x2826
    x2829 = x4 * (x1626 + x2818 - x2823 + x2824 + x987)
    x2830 = -x2822 + x2827
    x2831 = x4 * (x1007 + x1649 - x2819 + x2820 + x2825)
    x2832 = x1008 + x1261 * x1726 + x1651 - x1745 * x8
    x2833 = x1017 + x1261 * x1747 - x1749 * x8
    x2834 = x1261 * x1753
    x2835 = x1756 * x8
    x2836 = x1036 + x2834 - x2835
    x2837 = x2836 * x6
    x2838 = x1261 * x1749
    x2839 = x1753 * x8
    x2840 = x1054 + x2838 - x2839
    x2841 = -x2840 * x6
    x2842 = x2840 * x5
    x2843 = x2833 * x5 + x2841
    x2844 = x4 * (x1060 + x2833 - x2838 + x2839)
    x2845 = -x2837 + x2842
    x2846 = x4 * (x1080 - x2834 + x2835 + x2840)
    x2847 = x1082 + x1261 * x1756 - x1775 * x8
    x2848 = x1787 * x8
    x2849 = x1261 * x1781
    x2850 = 4.0 * x1669
    x2851 = x1091 + x1261 * x1779 + 4.0 * x1658 - x1781 * x8
    x2852 = x4 * (x1121 + x2848 - x2849 - x2850 + x2851)
    x2853 = x1105 - x2848 + x2849 + x2850
    x2854 = x1118 + x1261 * x1787 + 4.0 * x1684 - x1791 * x8
    x2855 = x1803 * x8
    x2856 = x1261 * x1797
    x2857 = 3.0 * x1699
    x2858 = x1128 + x1261 * x1795 + 3.0 * x1688 - x1797 * x8
    x2859 = x4 * (x1158 + x2855 - x2856 - x2857 + x2858)
    x2860 = x1142 - x2855 + x2856 + x2857
    x2861 = x1155 + x1261 * x1803 + 3.0 * x1714 - x1807 * x8
    x2862 = x1819 * x8
    x2863 = x1261 * x1813
    x2864 = 2.0 * x1729
    x2865 = -x2864
    x2866 = 2.0 * x1718
    x2867 = x1164 + x1261 * x1811 - x1813 * x8 + x2866
    x2868 = x4 * (x1192 + x2862 - x2863 + x2865 + x2867)
    x2869 = x1177 - x2862 + x2863 + x2864
    x2870 = 2.0 * x1744
    x2871 = x1189 + x1261 * x1819 - x1823 * x8 + x2870
    x2872 = x1833 * x8
    x2873 = x1261 * x1828
    x2874 = x1197 + x1261 * x1826 + x1748 - x1828 * x8
    x2875 = x4 * (x1223 + x1760 + x2872 - x2873 + x2874)
    x2876 = x1209 + x1759 - x2872 + x2873
    x2877 = x1220 + x1261 * x1833 + x1774 - x1837 * x8
    x2878 = x1847 * x8
    x2879 = x1261 * x1842
    x2880 = x1230 + x1261 * x1840 - x1842 * x8
    x2881 = x4 * (x1260 + x2878 - x2879 + x2880)
    x2882 = x1244 - x2878 + x2879
    x2883 = x1257 + x1261 * x1847 - x1851 * x8
    x2884 = x1261 * x1858
    x2885 = x1261 * x1855
    x2886 = x1951 + x2885
    x2887 = x4 * (x1955 - x2884 + x2886)
    x2888 = x2886 * x5
    x2889 = x1956 + x2884
    x2890 = x2889 * x6
    x2891 = x2888 - x2890
    x2892 = x2891 * x5
    x2893 = x2889 * x5
    x2894 = x1261 * x1864
    x2895 = x1965 + x2894
    x2896 = x2895 * x6
    x2897 = x2893 - x2896
    x2898 = x2897 * x6
    x2899 = x2887 + x2892 - x2898
    x2900 = x2899 * x5
    x2901 = x4 * (x1964 + x2889 - x2894)
    x2902 = x2897 * x5
    x2903 = x2895 * x5
    x2904 = x1261 * x1874
    x2905 = x1986 + x2904
    x2906 = x2905 * x6
    x2907 = x2903 - x2906
    x2908 = x2907 * x6
    x2909 = x2901 + x2902 - x2908
    x2910 = x2909 * x6
    x2911 = x38 * (x2891 - x2893 + x2896)
    x2912 = x2900 - x2910 + x2911
    x2913 = x2912 * x5
    x2914 = x1261 * x1884 + x1973
    x2915 = x4 * (x1950 - x2885 + x2914)
    x2916 = -x2886 * x6 + x2914 * x5
    x2917 = -x2891 * x6 + x2915 + x2916 * x5
    x2918 = -x2899 * x6 + x2917 * x5 + x38 * (-x2888 + x2890 + x2916)
    x2919 = -x2887
    x2920 = x130 * (-x2892 + x2898 + x2917 + x2919) - x2912 * x6
    x2921 = x2918 * x5 + x2920
    x2922 = x2909 * x5
    x2923 = x4 * (x1985 + x2895 - x2904)
    x2924 = x2907 * x5
    x2925 = x2905 * x5
    x2926 = x1261 * x1898
    x2927 = x2008 + x2926
    x2928 = x2927 * x6
    x2929 = x2925 - x2928
    x2930 = x2929 * x6
    x2931 = x2923 + x2924 - x2930
    x2932 = x2931 * x6
    x2933 = x38 * (x2897 - x2903 + x2906)
    x2934 = x2922 - x2932 + x2933
    x2935 = x2934 * x6
    x2936 = -x2901
    x2937 = x130 * (x2899 - x2902 + x2908 + x2936)
    x2938 = x2935 - x2937
    x2939 = x2913 - x2935 + x2937
    x2940 = x4 * (x2007 + x2905 - x2926)
    x2941 = x1261 * x1922 + x2028
    x2942 = -x2923
    x2943 = x1261 * x1957
    x2944 = x1261 * x1952
    x2945 = x2109 + x2944
    x2946 = x4 * (x2136 - x2943 + x2945)
    x2947 = x2945 * x5
    x2948 = x2114 + x2943
    x2949 = x2948 * x6
    x2950 = x2947 - x2949
    x2951 = x2950 * x5
    x2952 = x2948 * x5
    x2953 = x1261 * x1966
    x2954 = x2131 + x2953
    x2955 = x2954 * x6
    x2956 = x2952 - x2955
    x2957 = x2956 * x6
    x2958 = x2946 + x2951 - x2957
    x2959 = x2958 * x5
    x2960 = x1261 * x1974 + x2120
    x2961 = x4 * (x2124 - x2944 + x2960)
    x2962 = -x2945 * x6 + x2960 * x5
    x2963 = -x2950 * x6 + x2961 + x2962 * x5
    x2964 = -x2958 * x6 + x38 * (-x2947 + x2949 + x2962)
    x2965 = x2963 * x5 + x2964
    x2966 = x4 * (x2146 + x2948 - x2953)
    x2967 = x2956 * x5
    x2968 = x2954 * x5
    x2969 = x1261 * x1987
    x2970 = x2151 + x2969
    x2971 = x2970 * x6
    x2972 = x2968 - x2971
    x2973 = x2972 * x6
    x2974 = x2966 + x2967 - x2973
    x2975 = x2974 * x6
    x2976 = x38 * (x2950 - x2952 + x2955)
    x2977 = x2975 - x2976
    x2978 = x2959 - x2975 + x2976
    x2979 = x4 * (x2166 + x2954 - x2969)
    x2980 = x1261 * x2009 + x2168
    x2981 = x1261 * x2038
    x2982 = x1261 * x2035
    x2983 = x2175 + x2982
    x2984 = x4 * (x2179 - x2981 + x2983)
    x2985 = x2983 * x5
    x2986 = x2180 + x2981
    x2987 = x2986 * x6
    x2988 = x2985 - x2987
    x2989 = x2988 * x5
    x2990 = x2986 * x5
    x2991 = x1261 * x2045
    x2992 = x2196 + x2991
    x2993 = x2992 * x6
    x2994 = x2990 - x2993
    x2995 = x2994 * x6
    x2996 = x2984 + x2989 - x2995
    x2997 = x2996 * x5
    x2998 = x1261 * x2052 + x2186
    x2999 = x4 * (x2174 - x2982 + x2998)
    x3000 = -x2983 * x6 + x2998 * x5
    x3001 = -x2988 * x6 + x2999 + x3000 * x5
    x3002 = -x2996 * x6 + x38 * (-x2985 + x2987 + x3000)
    x3003 = x3001 * x5 + x3002
    x3004 = x4 * (x2195 + x2986 - x2991)
    x3005 = x2994 * x5
    x3006 = x2992 * x5
    x3007 = x1261 * x2063
    x3008 = x2214 + x3007
    x3009 = x3008 * x6
    x3010 = x3006 - x3009
    x3011 = x3010 * x6
    x3012 = x3004 + x3005 - x3011
    x3013 = x3012 * x6
    x3014 = x38 * (x2988 - x2990 + x2993)
    x3015 = x3013 - x3014
    x3016 = x2997 - x3013 + x3014
    x3017 = -x2984
    x3018 = x4 * (x2213 + x2992 - x3007)
    x3019 = x1261 * x2083 + x2230
    x3020 = -x3004
    x3021 = x1261 * x2110
    x3022 = x2294 + x3021
    x3023 = x3022 * x5
    x3024 = x1261 * x2115
    x3025 = x2301 + x3024
    x3026 = x3025 * x6
    x3027 = x3023 - x3026
    x3028 = x3027 * x5
    x3029 = x1261 * x2121 + x2288
    x3030 = -x3022 * x6 + x3029 * x5
    x3031 = x4 * (x2292 - x3021 + x3029)
    x3032 = -x3027 * x6 + x3031
    x3033 = x3030 * x5 + x3032
    x3034 = x3025 * x5
    x3035 = x1261 * x2132
    x3036 = x2307 + x3035
    x3037 = x3036 * x6
    x3038 = x3034 - x3037
    x3039 = x3038 * x6
    x3040 = x4 * (x2311 + x3022 - x3024)
    x3041 = x3039 - x3040
    x3042 = x3028 - x3039 + x3040
    x3043 = x4 * (x2327 + x3025 - x3035)
    x3044 = x1261 * x2152 + x2329
    x3045 = x1261 * x2176
    x3046 = x2337 + x3045
    x3047 = x3046 * x5
    x3048 = x1261 * x2181
    x3049 = x2343 + x3048
    x3050 = x3049 * x6
    x3051 = x3047 - x3050
    x3052 = x3051 * x5
    x3053 = x1261 * x2187 + x2332
    x3054 = -x3046 * x6 + x3053 * x5
    x3055 = x4 * (x2335 - x3045 + x3053)
    x3056 = -x3051 * x6 + x3055
    x3057 = x3054 * x5 + x3056
    x3058 = x3049 * x5
    x3059 = x1261 * x2197
    x3060 = x2348 + x3059
    x3061 = x3060 * x6
    x3062 = x3058 - x3061
    x3063 = x3062 * x6
    x3064 = x4 * (x2352 + x3046 - x3048)
    x3065 = x3063 - x3064
    x3066 = x3052 - x3063 + x3064
    x3067 = x4 * (x2368 + x3049 - x3059)
    x3068 = x1261 * x2215 + x2370
    x3069 = 241.359114645008 * x170
    x3070 = x1261 * x2236
    x3071 = x2377 + x3070
    x3072 = x3071 * x5
    x3073 = x1261 * x2239
    x3074 = x2383 + x3073
    x3075 = x3074 * x6
    x3076 = x3072 - x3075
    x3077 = x3076 * x5
    x3078 = x1261 * x2244 + x2374
    x3079 = -x3071 * x6 + x3078 * x5
    x3080 = x4 * (x2372 - x3070 + x3078)
    x3081 = -x3076 * x6 + x3080
    x3082 = x3079 * x5 + x3081
    x3083 = x3074 * x5
    x3084 = x1261 * x2252
    x3085 = x2388 + x3084
    x3086 = x3085 * x6
    x3087 = x3083 - x3086
    x3088 = x3087 * x6
    x3089 = x4 * (x2382 + x3071 - x3073)
    x3090 = -x3089
    x3091 = x3088 + x3090
    x3092 = x3077 - x3088 + x3089
    x3093 = x4 * (x2387 + x3074 - x3084)
    x3094 = x1261 * x2268 + x2408
    x3095 = x1261 * x2289 + x2441
    x3096 = x1261 * x2302
    x3097 = x2454 + x3096
    x3098 = x3097 * x6
    x3099 = x1261 * x2295
    x3100 = x2447 + x3099
    x3101 = -x3100 * x6
    x3102 = x3100 * x5
    x3103 = x3095 * x5 + x3101
    x3104 = x4 * (x2445 + x3095 - x3099)
    x3105 = -x3098 + x3102
    x3106 = x4 * (x2458 - x3096 + x3100)
    x3107 = x1261 * x2308 + x2460
    x3108 = x1261 * x2333 + x2463
    x3109 = x1261 * x2344
    x3110 = x2476 + x3109
    x3111 = x3110 * x6
    x3112 = x1261 * x2338
    x3113 = x2469 + x3112
    x3114 = -x3113 * x6
    x3115 = x3113 * x5
    x3116 = x3108 * x5 + x3114
    x3117 = x4 * (x2467 + x3108 - x3112)
    x3118 = -x3111 + x3115
    x3119 = x4 * (x2480 - x3109 + x3113)
    x3120 = x1261 * x2349 + x2482
    x3121 = x1261 * x2375 + x2485
    x3122 = x1261 * x2384
    x3123 = x2496 + x3122
    x3124 = x3123 * x6
    x3125 = x1261 * x2378
    x3126 = x2490 + x3125
    x3127 = -x3126 * x6
    x3128 = x3126 * x5
    x3129 = x3121 * x5 + x3127
    x3130 = x4 * (x2488 + x3121 - x3125)
    x3131 = -x3124 + x3128
    x3132 = x4 * (x2500 - x3122 + x3126)
    x3133 = x1261 * x2389 + x2502
    x3134 = x1261 * x2411 + x2506
    x3135 = x1261 * x2417
    x3136 = x2515 + x3135
    x3137 = x3136 * x6
    x3138 = x1261 * x2413
    x3139 = x2509 + x3138
    x3140 = -x3139 * x6
    x3141 = x3139 * x5
    x3142 = x3134 * x5 + x3140
    x3143 = x4 * (x2504 + x3134 - x3138)
    x3144 = -x3137 + x3141
    x3145 = x4 * (x2514 - x3135 + x3139)
    x3146 = x1261 * x2420 + x2520
    x3147 = x2455 * x8
    x3148 = x1261 * x2448
    x3149 = 4.0 * x2312
    x3150 = x1261 * x2442 + 4.0 * x2293 - x2448 * x8
    x3151 = x4 * (x3147 - x3148 - x3149 + x3150)
    x3152 = -x3147 + x3148 + x3149
    x3153 = x1261 * x2455 + 4.0 * x2328 - x2461 * x8
    x3154 = x2477 * x8
    x3155 = x1261 * x2470
    x3156 = 3.0 * x2353
    x3157 = x1261 * x2464 + 3.0 * x2336 - x2470 * x8
    x3158 = x4 * (x3154 - x3155 - x3156 + x3157)
    x3159 = -x3154 + x3155 + x3156
    x3160 = x1261 * x2477 + 3.0 * x2369 - x2483 * x8
    x3161 = x2497 * x8
    x3162 = x1261 * x2491
    x3163 = 2.0 * x2392
    x3164 = x1261 * x2486 + 2.0 * x2376 - x2491 * x8
    x3165 = x4 * (x3161 - x3162 - x3163 + x3164)
    x3166 = -x3161 + x3162 + x3163
    x3167 = x1261 * x2497 + 2.0 * x2407 - x2503 * x8
    x3168 = x2516 * x8
    x3169 = x1261 * x2510
    x3170 = x1261 * x2507 + x2412 - x2510 * x8
    x3171 = x4 * (x2424 + x3168 - x3169 + x3170)
    x3172 = x2423 - x3168 + x3169
    x3173 = x1261 * x2516 + x2438 - x2521 * x8
    x3174 = x2533 * x8
    x3175 = x1261 * x2527
    x3176 = x1261 * x2525 - x2527 * x8
    x3177 = x4 * (x3174 - x3175 + x3176)
    x3178 = -x3174 + x3175
    x3179 = x1261 * x2533 - x2537 * x8
    x3180 = x10 * x1864
    x3181 = x1852 * x1858
    x3182 = x1852 * x1855
    x3183 = x10 * x1858
    x3184 = x19 + x3182 - x3183
    x3185 = x4 * (x3180 - x3181 + x3184 + x69)
    x3186 = x3184 * x5
    x3187 = x29 - x3180 + x3181
    x3188 = x3187 * x6
    x3189 = x3186 - x3188
    x3190 = x3189 * x5
    x3191 = x3187 * x5
    x3192 = x1852 * x1864
    x3193 = x10 * x1874
    x3194 = x3192 - x3193 + x56
    x3195 = x3194 * x6
    x3196 = x3191 - x3195
    x3197 = x3196 * x6
    x3198 = x3185 + x3190 - x3197
    x3199 = x3198 * x5
    x3200 = x4 * (x3187 - x3192 + x3193 + x95)
    x3201 = x3196 * x5
    x3202 = x3194 * x5
    x3203 = x1852 * x1874
    x3204 = x10 * x1898
    x3205 = x3203 - x3204 + x82
    x3206 = x3205 * x6
    x3207 = x3202 - x3206
    x3208 = x3207 * x6
    x3209 = x3200 + x3201 - x3208
    x3210 = x3209 * x6
    x3211 = x38 * (x3189 - x3191 + x3195)
    x3212 = x3199 - x3210 + x3211
    x3213 = x3212 * x5
    x3214 = -x10 * x1855 + x1852 * x1884 + x45
    x3215 = x4 * (-x3182 + x3183 + x3214 + x50)
    x3216 = -x3184 * x6 + x3214 * x5
    x3217 = -x3189 * x6 + x3215 + x3216 * x5
    x3218 = -x3198 * x6 + x3217 * x5 + x38 * (-x3186 + x3188 + x3216)
    x3219 = -x3185
    x3220 = x130 * (-x3190 + x3197 + x3217 + x3219) - x3212 * x6
    x3221 = x3218 * x5 + x3220
    x3222 = x3209 * x5
    x3223 = x4 * (x129 + x3194 - x3203 + x3204)
    x3224 = x3207 * x5
    x3225 = x3205 * x5
    x3226 = x1852 * x1898
    x3227 = x10 * x1922
    x3228 = x117 + x3226 - x3227
    x3229 = x3228 * x6
    x3230 = x3225 - x3229
    x3231 = x3230 * x6
    x3232 = x3223 + x3224 - x3231
    x3233 = x3232 * x6
    x3234 = x38 * (x3196 - x3202 + x3206)
    x3235 = x3222 - x3233 + x3234
    x3236 = x3235 * x6
    x3237 = -x3200
    x3238 = x130 * (x3198 - x3201 + x3208 + x3237)
    x3239 = x3236 - x3238
    x3240 = x3213 - x3236 + x3238
    x3241 = x4 * (x157 + x3205 - x3226 + x3227)
    x3242 = -x10 * x1944 + x155 + x1852 * x1922
    x3243 = -x3223
    x3244 = x3194 * x8
    x3245 = x174 * x3187
    x3246 = x174 * x3184
    x3247 = x3187 * x8
    x3248 = -x3247
    x3249 = x3246 + x3248
    x3250 = x4 * (x3244 - x3245 + x3249)
    x3251 = x3249 * x5
    x3252 = -x3244
    x3253 = x3245 + x3252
    x3254 = x3253 * x6
    x3255 = x3251 - x3254
    x3256 = x3255 * x5
    x3257 = x3253 * x5
    x3258 = x174 * x3194
    x3259 = x3205 * x8
    x3260 = -x3259
    x3261 = x3258 + x3260
    x3262 = x3261 * x6
    x3263 = x3257 - x3262
    x3264 = x3263 * x6
    x3265 = x3250 + x3256 - x3264
    x3266 = x3265 * x5
    x3267 = -x3184 * x8
    x3268 = x174 * x3214 + x3267
    x3269 = x4 * (-x3246 + x3247 + x3268)
    x3270 = -x3249 * x6 + x3268 * x5
    x3271 = -x3255 * x6 + x3269 + x3270 * x5
    x3272 = -x3265 * x6 + x38 * (-x3251 + x3254 + x3270)
    x3273 = x3271 * x5 + x3272
    x3274 = x4 * (x3253 - x3258 + x3259)
    x3275 = x3263 * x5
    x3276 = x3261 * x5
    x3277 = x174 * x3205
    x3278 = x3228 * x8
    x3279 = -x3278
    x3280 = x3277 + x3279
    x3281 = x3280 * x6
    x3282 = x3276 - x3281
    x3283 = x3282 * x6
    x3284 = x3274 + x3275 - x3283
    x3285 = x3284 * x6
    x3286 = x38 * (x3255 - x3257 + x3262)
    x3287 = x3285 - x3286
    x3288 = x3266 - x3285 + x3286
    x3289 = -x3250
    x3290 = x4 * (x3261 - x3277 + x3278)
    x3291 = -x3242 * x8
    x3292 = x174 * x3228 + x3291
    x3293 = -x3274
    x3294 = x10 * x2045
    x3295 = x1852 * x2038
    x3296 = x1852 * x2035
    x3297 = x10 * x2038
    x3298 = x1856 + x3296 - x3297 + x333
    x3299 = x4 * (x1907 + x3294 - x3295 + x3298 + x398)
    x3300 = x3298 * x5
    x3301 = x1870 - x3294 + x3295 + x357
    x3302 = x3301 * x6
    x3303 = x3300 - x3302
    x3304 = x3303 * x5
    x3305 = x3301 * x5
    x3306 = x1852 * x2045
    x3307 = x10 * x2063
    x3308 = x1894 + x3306 - x3307 + x381
    x3309 = x3308 * x6
    x3310 = x3305 - x3309
    x3311 = x3310 * x6
    x3312 = x3299 + x3304 - x3311
    x3313 = x3312 * x5
    x3314 = -x10 * x2035 + x1852 * x2052 + x1885 + x350
    x3315 = x4 * (x1890 - x3296 + x3297 + x3314 + x374)
    x3316 = -x3298 * x6 + x3314 * x5
    x3317 = -x3303 * x6 + x3315 + x3316 * x5
    x3318 = -x3312 * x6 + x38 * (-x3300 + x3302 + x3316)
    x3319 = x3317 * x5 + x3318
    x3320 = x4 * (x1931 + x3301 - x3306 + x3307 + x428)
    x3321 = x3310 * x5
    x3322 = x3308 * x5
    x3323 = x1852 * x2063
    x3324 = x10 * x2083
    x3325 = x1918 + x3323 - x3324 + x412
    x3326 = x3325 * x6
    x3327 = x3322 - x3326
    x3328 = x3327 * x6
    x3329 = x3320 + x3321 - x3328
    x3330 = x3329 * x6
    x3331 = x38 * (x3303 - x3305 + x3309)
    x3332 = x3330 - x3331
    x3333 = x3313 - x3330 + x3331
    x3334 = -x3299
    x3335 = x4 * (x1945 + x3308 - x3323 + x3324 + x452)
    x3336 = -x10 * x2102 + x1852 * x2083 + x1943 + x449
    x3337 = -x3320
    x3338 = x174 * x3249
    x3339 = x3253 * x8
    x3340 = x3185 - x3339
    x3341 = x3338 + x3340
    x3342 = x3341 * x5
    x3343 = x174 * x3253
    x3344 = x3261 * x8
    x3345 = x3200 - x3344
    x3346 = x3343 + x3345
    x3347 = x3346 * x6
    x3348 = x3342 - x3347
    x3349 = x3348 * x5
    x3350 = x3215 - x3249 * x8
    x3351 = x174 * x3268 + x3350
    x3352 = -x3341 * x6 + x3351 * x5
    x3353 = x3219 + x3339
    x3354 = x4 * (-x3338 + x3351 + x3353)
    x3355 = -x3348 * x6 + x3354
    x3356 = x3352 * x5 + x3355
    x3357 = x3346 * x5
    x3358 = x174 * x3261
    x3359 = x3280 * x8
    x3360 = x3223 - x3359
    x3361 = x3358 + x3360
    x3362 = x3361 * x6
    x3363 = x3357 - x3362
    x3364 = x3363 * x6
    x3365 = x3237 + x3344
    x3366 = x4 * (x3341 - x3343 + x3365)
    x3367 = -x3366
    x3368 = x3364 + x3367
    x3369 = x3349 - x3364 + x3366
    x3370 = x3243 + x3359
    x3371 = x4 * (x3346 - x3358 + x3370)
    x3372 = x3241 - x3292 * x8
    x3373 = x174 * x3280 + x3372
    x3374 = x174 * x3298
    x3375 = x3301 * x8
    x3376 = -x3375
    x3377 = x3374 + x3376
    x3378 = x3377 * x5
    x3379 = x174 * x3301
    x3380 = x3308 * x8
    x3381 = -x3380
    x3382 = x3379 + x3381
    x3383 = x3382 * x6
    x3384 = x3378 - x3383
    x3385 = x3384 * x5
    x3386 = -x3298 * x8
    x3387 = x174 * x3314 + x3386
    x3388 = -x3377 * x6 + x3387 * x5
    x3389 = x4 * (-x3374 + x3375 + x3387)
    x3390 = -x3384 * x6 + x3389
    x3391 = x3388 * x5 + x3390
    x3392 = x3382 * x5
    x3393 = x174 * x3308
    x3394 = x3325 * x8
    x3395 = -x3394
    x3396 = x3393 + x3395
    x3397 = x3396 * x6
    x3398 = x3392 - x3397
    x3399 = x3398 * x6
    x3400 = x4 * (x3377 - x3379 + x3380)
    x3401 = -x3400
    x3402 = x3399 + x3401
    x3403 = x3385 - x3399 + x3400
    x3404 = x4 * (x3382 - x3393 + x3394)
    x3405 = -x3336 * x8
    x3406 = x174 * x3325 + x3405
    x3407 = x1852 * x2236
    x3408 = x10 * x2239
    x3409 = 2.0 * x2042
    x3410 = x3407 - x3408 + x3409 + x719
    x3411 = x3410 * x5
    x3412 = x1852 * x2239
    x3413 = x10 * x2252
    x3414 = 2.0 * x2059
    x3415 = x3412 - x3413 + x3414 + x729
    x3416 = x3415 * x6
    x3417 = x3411 - x3416
    x3418 = x3417 * x5
    x3419 = 2.0 * x2054
    x3420 = -x10 * x2236 + x1852 * x2244 + x3419 + x706
    x3421 = -x3410 * x6 + x3420 * x5
    x3422 = -x3409
    x3423 = x4 * (-x3407 + x3408 + x3420 + x3422 + x720)
    x3424 = -x3417 * x6 + x3423
    x3425 = x3421 * x5 + x3424
    x3426 = x3415 * x5
    x3427 = x1852 * x2252
    x3428 = x10 * x2268
    x3429 = 2.0 * x2079
    x3430 = x3427 - x3428 + x3429 + x756
    x3431 = x3430 * x6
    x3432 = x3426 - x3431
    x3433 = x3432 * x6
    x3434 = -x3414
    x3435 = x4 * (x3410 - x3412 + x3413 + x3434 + x749)
    x3436 = -x3435
    x3437 = x3433 + x3436
    x3438 = x3418 - x3433 + x3435
    x3439 = -x3429
    x3440 = x4 * (x3415 - x3427 + x3428 + x3439 + x785)
    x3441 = 2.0 * x2101
    x3442 = -x10 * x2283 + x1852 * x2268 + x3441 + x788
    x3443 = 2.0 * x3269 - x3341 * x8
    x3444 = x174 * x3351 + x3443
    x3445 = x174 * x3346
    x3446 = x3361 * x8
    x3447 = 2.0 * x3274
    x3448 = -x3446 + x3447
    x3449 = x3445 + x3448
    x3450 = x3449 * x6
    x3451 = x174 * x3341
    x3452 = x3346 * x8
    x3453 = 2.0 * x3250
    x3454 = -x3452 + x3453
    x3455 = x3451 + x3454
    x3456 = -x3455 * x6
    x3457 = x3455 * x5
    x3458 = x3444 * x5 + x3456
    x3459 = x3452 - x3453
    x3460 = x4 * (x3444 - x3451 + x3459)
    x3461 = -x3450 + x3457
    x3462 = x3446 - x3447
    x3463 = x4 * (-x3445 + x3455 + x3462)
    x3464 = 2.0 * x3290 - x3373 * x8
    x3465 = x174 * x3361 + x3464
    x3466 = x3315 - x3377 * x8
    x3467 = x174 * x3387 + x3466
    x3468 = x174 * x3382
    x3469 = x3396 * x8
    x3470 = x3320 - x3469
    x3471 = x3468 + x3470
    x3472 = x3471 * x6
    x3473 = x174 * x3377
    x3474 = x3382 * x8
    x3475 = x3299 - x3474
    x3476 = x3473 + x3475
    x3477 = -x3476 * x6
    x3478 = x3476 * x5
    x3479 = x3467 * x5 + x3477
    x3480 = x3334 + x3474
    x3481 = x4 * (x3467 - x3473 + x3480)
    x3482 = -x3472 + x3478
    x3483 = x3337 + x3469
    x3484 = x4 * (-x3468 + x3476 + x3483)
    x3485 = x3335 - x3406 * x8
    x3486 = x174 * x3396 + x3485
    x3487 = -x3410 * x8
    x3488 = x174 * x3420 + x3487
    x3489 = x174 * x3415
    x3490 = x3430 * x8
    x3491 = -x3490
    x3492 = x3489 + x3491
    x3493 = x3492 * x6
    x3494 = x174 * x3410
    x3495 = x3415 * x8
    x3496 = -x3495
    x3497 = x3494 + x3496
    x3498 = -x3497 * x6
    x3499 = x3497 * x5
    x3500 = x3488 * x5 + x3498
    x3501 = x4 * (x3488 - x3494 + x3495)
    x3502 = -x3493 + x3499
    x3503 = x4 * (-x3489 + x3490 + x3497)
    x3504 = -x3442 * x8
    x3505 = x174 * x3430 + x3504
    x3506 = -x10 * x2413 + x1017 + x1852 * x2411 + 3.0 * x2247
    x3507 = x1852 * x2417
    x3508 = x10 * x2420
    x3509 = 3.0 * x2265
    x3510 = x1036 + x3507 - x3508 + x3509
    x3511 = x3510 * x6
    x3512 = x1852 * x2413
    x3513 = x10 * x2417
    x3514 = 3.0 * x2256
    x3515 = x1054 + x3512 - x3513 + x3514
    x3516 = -x3515 * x6
    x3517 = x3515 * x5
    x3518 = x3506 * x5 + x3516
    x3519 = x4 * (x1060 + x3506 - x3512 + x3513 - x3514)
    x3520 = -x3511 + x3517
    x3521 = x4 * (x1080 - x3507 + x3508 - x3509 + x3515)
    x3522 = -x10 * x2439 + x1082 + x1852 * x2420 + 3.0 * x2282
    x3523 = x174 * x3455
    x3524 = 3.0 * x3354 - x3455 * x8
    x3525 = x174 * x3444 + x3524
    x3526 = x3449 * x8
    x3527 = 3.0 * x3366
    x3528 = x3526 - x3527
    x3529 = x4 * (-x3523 + x3525 + x3528)
    x3530 = -x3526 + x3527
    x3531 = x3523 + x3530
    x3532 = 3.0 * x3371 - x3465 * x8
    x3533 = x174 * x3449 + x3532
    x3534 = x174 * x3476
    x3535 = 2.0 * x3389 - x3476 * x8
    x3536 = x174 * x3467 + x3535
    x3537 = x3471 * x8
    x3538 = 2.0 * x3400
    x3539 = x3537 - x3538
    x3540 = x4 * (-x3534 + x3536 + x3539)
    x3541 = -x3537 + x3538
    x3542 = x3534 + x3541
    x3543 = 2.0 * x3404 - x3486 * x8
    x3544 = x174 * x3471 + x3543
    x3545 = x174 * x3497
    x3546 = x3423 - x3497 * x8
    x3547 = x174 * x3488 + x3546
    x3548 = x3492 * x8
    x3549 = x3436 + x3548
    x3550 = x4 * (-x3545 + x3547 + x3549)
    x3551 = x3435 - x3548
    x3552 = x3545 + x3551
    x3553 = x3440 - x3505 * x8
    x3554 = x174 * x3492 + x3553
    x3555 = x3510 * x8
    x3556 = x174 * x3515
    x3557 = -x3515 * x8
    x3558 = x174 * x3506 + x3557
    x3559 = x4 * (x3555 - x3556 + x3558)
    x3560 = -x3555
    x3561 = x3556 + x3560
    x3562 = -x3522 * x8
    x3563 = x174 * x3510 + x3562
    x3564 = x10 * x2533
    x3565 = x1852 * x2527
    x3566 = 4.0 * x2423
    x3567 = -x10 * x2527 + x1230 + x1852 * x2525 + 4.0 * x2412
    x3568 = x4 * (x1260 + x3564 - x3565 - x3566 + x3567)
    x3569 = x1244 - x3564 + x3565 + x3566
    x3570 = -x10 * x2537 + x1257 + x1852 * x2533 + 4.0 * x2438
    x3571 = x2545 * x8
    x3572 = x1261 * x2542
    x3573 = 2.0 * x1265
    x3574 = x1261 * x2572 + 2.0 * x1294 - x2542 * x8
    x3575 = x4 * (x3571 - x3572 - x3573 + x3574)
    x3576 = -x3571 + x3572 + x3573
    x3577 = x3574 * x5 - x3576 * x6
    x3578 = x3576 * x5
    x3579 = x1261 * x2545
    x3580 = x2552 * x8
    x3581 = 2.0 * x1279
    x3582 = x3579 - x3580 + x3581
    x3583 = x3582 * x6
    x3584 = x3578 - x3583
    x3585 = x3575 + x3577 * x5 - x3584 * x6
    x3586 = x4 * (x3576 - x3579 + x3580 - x3581)
    x3587 = x3584 * x5
    x3588 = x3582 * x5
    x3589 = x1261 * x2552
    x3590 = x2563 * x8
    x3591 = 2.0 * x1303
    x3592 = x3589 - x3590 + x3591
    x3593 = x3592 * x6
    x3594 = x3588 - x3593
    x3595 = x3594 * x6
    x3596 = x3586 + x3587 - x3595
    x3597 = x3585 * x5 - x3596 * x6 + x38 * (x3577 - x3578 + x3583)
    x3598 = x3596 * x5
    x3599 = x4 * (x3582 - x3589 + x3590 - x3591)
    x3600 = x3594 * x5
    x3601 = x3592 * x5
    x3602 = x1261 * x2563
    x3603 = x2586 * x8
    x3604 = 2.0 * x1327
    x3605 = x3602 - x3603 + x3604
    x3606 = x3605 * x6
    x3607 = x3601 - x3606
    x3608 = x3607 * x6
    x3609 = x3599 + x3600 - x3608
    x3610 = x3609 * x6
    x3611 = x38 * (x3584 - x3588 + x3593)
    x3612 = x3598 - x3610 + x3611
    x3613 = -x3586
    x3614 = x4 * (x3592 - x3602 + x3603 - x3604)
    x3615 = x1261 * x2586 + 2.0 * x1352 - x2600 * x8
    x3616 = -x3599
    x3617 = x2609 * x8
    x3618 = x1261 * x2606
    x3619 = x1261 * x2622 + x2573 - x2606 * x8 + x2700
    x3620 = x4 * (x2577 + x2703 + x3617 - x3618 + x3619)
    x3621 = x2543 + x2690 - x3617 + x3618
    x3622 = x3619 * x5 - x3621 * x6
    x3623 = x3621 * x5
    x3624 = x1261 * x2609
    x3625 = x2616 * x8
    x3626 = x2558 + x2695 + x3624 - x3625
    x3627 = x3626 * x6
    x3628 = x3623 - x3627
    x3629 = x3620 + x3622 * x5 - x3628 * x6
    x3630 = x4 * (x2595 + x2715 + x3621 - x3624 + x3625)
    x3631 = x3628 * x5
    x3632 = x3626 * x5
    x3633 = x1261 * x2616
    x3634 = x2633 * x8
    x3635 = x2581 + x2710 + x3633 - x3634
    x3636 = x3635 * x6
    x3637 = x3632 - x3636
    x3638 = x3637 * x6
    x3639 = x3630 + x3631 - x3638
    x3640 = x4 * (x2601 + x2719 + x3626 - x3633 + x3634)
    x3641 = x1261 * x2633 + x2599 - x2643 * x8 + x2721
    x3642 = x2651 * x8
    x3643 = x1261 * x2648
    x3644 = 2.0 * x1437
    x3645 = x1261 * x2664 + 2.0 * x1449 - x2648 * x8
    x3646 = x4 * (x3642 - x3643 - x3644 + x3645)
    x3647 = -x3642 + x3643 + x3644
    x3648 = x3645 * x5 - x3647 * x6
    x3649 = x3647 * x5
    x3650 = x1261 * x2651
    x3651 = x2658 * x8
    x3652 = 2.0 * x1454
    x3653 = x3650 - x3651 + x3652
    x3654 = x3653 * x6
    x3655 = x3649 - x3654
    x3656 = x3646 + x3648 * x5 - x3655 * x6
    x3657 = x4 * (x3647 - x3650 + x3651 - x3652)
    x3658 = x3655 * x5
    x3659 = x3653 * x5
    x3660 = x1261 * x2658
    x3661 = x2675 * x8
    x3662 = 2.0 * x1474
    x3663 = x3660 - x3661 + x3662
    x3664 = x3663 * x6
    x3665 = x3659 - x3664
    x3666 = x3665 * x6
    x3667 = x3657 + x3658 - x3666
    x3668 = x4 * (x3653 - x3660 + x3661 - x3662)
    x3669 = x1261 * x2675 + 2.0 * x1496 - x2686 * x8
    x3670 = -x3657
    x3671 = x2696 * x8
    x3672 = x1261 * x2691
    x3673 = 2.0 * x1523
    x3674 = 2.0 * x2607
    x3675 = x1261 * x2701 + 2.0 * x1514 + 2.0 * x2623 - x2691 * x8
    x3676 = x4 * (x3671 - x3672 - x3673 - x3674 + x3675)
    x3677 = -x3671 + x3672 + x3673 + x3674
    x3678 = x3675 * x5 - x3677 * x6
    x3679 = x3677 * x5
    x3680 = x1261 * x2696
    x3681 = x2711 * x8
    x3682 = 2.0 * x1531
    x3683 = 2.0 * x2628
    x3684 = x3680 - x3681 + x3682 + x3683
    x3685 = x3684 * x6
    x3686 = x3679 - x3685
    x3687 = x4 * (x3677 - x3680 + x3681 - x3682 - x3683)
    x3688 = x1261 * x2711 + 2.0 * x1547 + 2.0 * x2642 - x2722 * x8
    x3689 = x2730 * x8
    x3690 = x1261 * x2726
    x3691 = x1261 * x2734 + x2665 - x2726 * x8 + x2797
    x3692 = x4 * (x2684 + x2811 + x3689 - x3690 + x3691)
    x3693 = x2649 + x2806 - x3689 + x3690
    x3694 = x3691 * x5 - x3693 * x6
    x3695 = x3693 * x5
    x3696 = x1261 * x2730
    x3697 = x2742 * x8
    x3698 = x2670 + x2801 + x3696 - x3697
    x3699 = x3698 * x6
    x3700 = x3695 - x3699
    x3701 = x4 * (x2687 + x2814 + x3693 - x3696 + x3697)
    x3702 = x1261 * x2742 + x2685 - x2750 * x8 + x2816
    x3703 = x2758 * x8
    x3704 = x1261 * x2754
    x3705 = 2.0 * x1625
    x3706 = x1261 * x2762 + 2.0 * x1616 - x2754 * x8
    x3707 = x4 * (x3703 - x3704 - x3705 + x3706)
    x3708 = -x3703 + x3704 + x3705
    x3709 = x3706 * x5 - x3708 * x6
    x3710 = x3708 * x5
    x3711 = x1261 * x2758
    x3712 = x2770 * x8
    x3713 = 2.0 * x1634
    x3714 = x3711 - x3712 + x3713
    x3715 = x3714 * x6
    x3716 = x3710 - x3715
    x3717 = x4 * (x3708 - x3711 + x3712 - x3713)
    x3718 = x1261 * x2770 + 2.0 * x1651 - x2779 * x8
    x3719 = x2784 * x8
    x3720 = x1261 * x2789
    x3721 = 3.0 * x2716
    x3722 = -x3721
    x3723 = 2.0 * x1669
    x3724 = 3.0 * x2704
    x3725 = x1261 * x2780 + 2.0 * x1658 - x2789 * x8 + x3724
    x3726 = x4 * (x3719 - x3720 + x3722 - x3723 + x3725)
    x3727 = -x3719 + x3720 + x3721 + x3723
    x3728 = 3.0 * x2720
    x3729 = x1261 * x2784 + 2.0 * x1684 - x2796 * x8 + x3728
    x3730 = x2802 * x8
    x3731 = x1261 * x2807
    x3732 = 2.0 * x1699
    x3733 = 2.0 * x2746
    x3734 = x1261 * x2798 + 2.0 * x1688 + 2.0 * x2736 - x2807 * x8
    x3735 = x4 * (x3730 - x3731 - x3732 - x3733 + x3734)
    x3736 = -x3730 + x3731 + x3732 + x3733
    x3737 = x1261 * x2802 + 2.0 * x1714 + 2.0 * x2749 - x2817 * x8
    x3738 = x2821 * x8
    x3739 = x1261 * x2825
    x3740 = x1261 * x2818 + x2764 - x2825 * x8 + x2866
    x3741 = x4 * (x2775 + x2865 + x3738 - x3739 + x3740)
    x3742 = x2774 + x2864 - x3738 + x3739
    x3743 = x1261 * x2821 + x2778 - x2832 * x8 + x2870
    x3744 = x2836 * x8
    x3745 = x1261 * x2840
    x3746 = 2.0 * x1759
    x3747 = x1261 * x2833 + 2.0 * x1748 - x2840 * x8
    x3748 = x4 * (x3744 - x3745 - x3746 + x3747)
    x3749 = -x3744 + x3745 + x3746
    x3750 = x1261 * x2836 + 2.0 * x1774 - x2847 * x8
    x3751 = x1261 * x2851 + 2.0 * x1780 + 4.0 * x2793 - x2853 * x8
    x3752 = x1261 * x2853 + 2.0 * x1790 + 4.0 * x2795 - x2854 * x8
    x3753 = 3.0 * x2812
    x3754 = x1261 * x2858 + 2.0 * x1796 - x2860 * x8 + x3753
    x3755 = 3.0 * x2815
    x3756 = x1261 * x2860 + 2.0 * x1806 - x2861 * x8 + x3755
    x3757 = x1261 * x2867 + 2.0 * x1812 + 2.0 * x2829 - x2869 * x8
    x3758 = x1261 * x2869 + 2.0 * x1822 + 2.0 * x2831 - x2871 * x8
    x3759 = x1261 * x2874 + 2.0 * x1827 + x2844 - x2876 * x8
    x3760 = x1261 * x2876 + 2.0 * x1836 + x2846 - x2877 * x8
    x3761 = x1261 * x2880 + 2.0 * x1841 - x2882 * x8
    x3762 = x1261 * x2882 + 2.0 * x1850 - x2883 * x8
    x3763 = x2889 * x8
    x3764 = x1261 * x2886
    x3765 = x1261 * x2914 + x1885 - x2886 * x8
    x3766 = x4 * (x1890 + x3763 - x3764 + x3765)
    x3767 = x1856 - x3763 + x3764
    x3768 = x3765 * x5 - x3767 * x6
    x3769 = x3767 * x5
    x3770 = x1261 * x2889
    x3771 = x2895 * x8
    x3772 = x1870 + x3770 - x3771
    x3773 = x3772 * x6
    x3774 = x3769 - x3773
    x3775 = x3766 + x3768 * x5 - x3774 * x6
    x3776 = x4 * (x1907 + x3767 - x3770 + x3771)
    x3777 = x3774 * x5
    x3778 = x3772 * x5
    x3779 = x1261 * x2895
    x3780 = x2905 * x8
    x3781 = x1894 + x3779 - x3780
    x3782 = x3781 * x6
    x3783 = x3778 - x3782
    x3784 = x3783 * x6
    x3785 = x3776 + x3777 - x3784
    x3786 = x3775 * x5 - x3785 * x6 + x38 * (x3768 - x3769 + x3773)
    x3787 = x3785 * x5
    x3788 = x4 * (x1931 + x3772 - x3779 + x3780)
    x3789 = x3783 * x5
    x3790 = x3781 * x5
    x3791 = x1261 * x2905
    x3792 = x2927 * x8
    x3793 = x1918 + x3791 - x3792
    x3794 = x3793 * x6
    x3795 = x3790 - x3794
    x3796 = x3795 * x6
    x3797 = x3788 + x3789 - x3796
    x3798 = x3797 * x6
    x3799 = x38 * (x3774 - x3778 + x3782)
    x3800 = x3787 - x3798 + x3799
    x3801 = -x3776
    x3802 = x4 * (x1945 + x3781 - x3791 + x3792)
    x3803 = x1261 * x2927 + x1943 - x2941 * x8
    x3804 = -x3788
    x3805 = x2948 * x8
    x3806 = x1261 * x2945
    x3807 = x1261 * x2960 + x1976 + x2915 - x2945 * x8
    x3808 = x4 * (x1997 + x2919 + x3805 - x3806 + x3807)
    x3809 = x1961 + x2887 - x3805 + x3806
    x3810 = x3807 * x5 - x3809 * x6
    x3811 = x3809 * x5
    x3812 = x1261 * x2948
    x3813 = x2954 * x8
    x3814 = x1981 + x2901 + x3812 - x3813
    x3815 = x3814 * x6
    x3816 = x3811 - x3815
    x3817 = x3808 + x3810 * x5 - x3816 * x6
    x3818 = x4 * (x2019 + x2936 + x3809 - x3812 + x3813)
    x3819 = x3816 * x5
    x3820 = x3814 * x5
    x3821 = x1261 * x2954
    x3822 = x2970 * x8
    x3823 = x2003 + x2923 + x3821 - x3822
    x3824 = x3823 * x6
    x3825 = x3820 - x3824
    x3826 = x3825 * x6
    x3827 = x3818 + x3819 - x3826
    x3828 = x4 * (x2030 + x2942 + x3814 - x3821 + x3822)
    x3829 = x1261 * x2970 + x2027 + x2940 - x2980 * x8
    x3830 = x2986 * x8
    x3831 = x1261 * x2983
    x3832 = x1261 * x2998 + x2054 - x2983 * x8
    x3833 = x4 * (x2073 + x3830 - x3831 + x3832)
    x3834 = x2042 - x3830 + x3831
    x3835 = x3832 * x5 - x3834 * x6
    x3836 = x3834 * x5
    x3837 = x1261 * x2986
    x3838 = x2992 * x8
    x3839 = x2059 + x3837 - x3838
    x3840 = x3839 * x6
    x3841 = x3836 - x3840
    x3842 = x3833 + x3835 * x5 - x3841 * x6
    x3843 = x4 * (x2093 + x3834 - x3837 + x3838)
    x3844 = x3841 * x5
    x3845 = x3839 * x5
    x3846 = x1261 * x2992
    x3847 = x3008 * x8
    x3848 = x2079 + x3846 - x3847
    x3849 = x3848 * x6
    x3850 = x3845 - x3849
    x3851 = x3850 * x6
    x3852 = x3843 + x3844 - x3851
    x3853 = x4 * (x2103 + x3839 - x3846 + x3847)
    x3854 = x1261 * x3008 + x2101 - x3019 * x8
    x3855 = -x3843
    x3856 = x3025 * x8
    x3857 = x1261 * x3022
    x3858 = 2.0 * x2946
    x3859 = -x3858
    x3860 = 2.0 * x2961
    x3861 = x1261 * x3029 + x2125 - x3022 * x8 + x3860
    x3862 = x4 * (x2138 + x3856 - x3857 + x3859 + x3861)
    x3863 = x2137 - x3856 + x3857 + x3858
    x3864 = x3861 * x5 - x3863 * x6
    x3865 = x3863 * x5
    x3866 = x1261 * x3025
    x3867 = x3036 * x8
    x3868 = 2.0 * x2966
    x3869 = x2147 + x3866 - x3867 + x3868
    x3870 = x3869 * x6
    x3871 = x3865 - x3870
    x3872 = -x3868
    x3873 = x4 * (x2164 + x3863 - x3866 + x3867 + x3872)
    x3874 = 2.0 * x2979
    x3875 = x1261 * x3036 + x2167 - x3044 * x8 + x3874
    x3876 = x3049 * x8
    x3877 = x1261 * x3046
    x3878 = x1261 * x3053 + x2190 + x2999 - x3046 * x8
    x3879 = x4 * (x2202 + x3017 + x3876 - x3877 + x3878)
    x3880 = x2201 + x2984 - x3876 + x3877
    x3881 = x3878 * x5 - x3880 * x6
    x3882 = x3880 * x5
    x3883 = x1261 * x3049
    x3884 = x3060 * x8
    x3885 = x2210 + x3004 + x3883 - x3884
    x3886 = x3885 * x6
    x3887 = x3882 - x3886
    x3888 = x4 * (x2227 + x3020 + x3880 - x3883 + x3884)
    x3889 = x1261 * x3060 + x2229 + x3018 - x3068 * x8
    x3890 = x3074 * x8
    x3891 = x1261 * x3071
    x3892 = x1261 * x3078 + x2247 - x3071 * x8
    x3893 = x4 * (x2257 + x3890 - x3891 + x3892)
    x3894 = x2256 - x3890 + x3891
    x3895 = x3892 * x5 - x3894 * x6
    x3896 = x3894 * x5
    x3897 = x1261 * x3074
    x3898 = x3085 * x8
    x3899 = x2265 + x3897 - x3898
    x3900 = x3899 * x6
    x3901 = x3896 - x3900
    x3902 = x4 * (x2280 + x3894 - x3897 + x3898)
    x3903 = x1261 * x3085 + x2282 - x3094 * x8
    x3904 = x3097 * x8
    x3905 = x1261 * x3100
    x3906 = 3.0 * x3040
    x3907 = x1261 * x3095 + x2293 + 3.0 * x3031 - x3100 * x8
    x3908 = x4 * (x2313 + x3904 - x3905 - x3906 + x3907)
    x3909 = x2312 - x3904 + x3905 + x3906
    x3910 = x1261 * x3097 + x2328 + 3.0 * x3043 - x3107 * x8
    x3911 = x3110 * x8
    x3912 = x1261 * x3113
    x3913 = 2.0 * x3064
    x3914 = -x3913
    x3915 = 2.0 * x3055
    x3916 = x1261 * x3108 + x2336 - x3113 * x8 + x3915
    x3917 = x4 * (x2354 + x3911 - x3912 + x3914 + x3916)
    x3918 = x2353 - x3911 + x3912 + x3913
    x3919 = 2.0 * x3067
    x3920 = x1261 * x3110 + x2369 - x3120 * x8 + x3919
    x3921 = x3123 * x8
    x3922 = x1261 * x3126
    x3923 = x1261 * x3121 + x2376 + x3080 - x3126 * x8
    x3924 = x4 * (x2393 + x3090 + x3921 - x3922 + x3923)
    x3925 = x2392 + x3089 - x3921 + x3922
    x3926 = x1261 * x3123 + x2407 + x3093 - x3133 * x8
    x3927 = x3136 * x8
    x3928 = x1261 * x3139
    x3929 = x1261 * x3134 + x2412 - x3139 * x8
    x3930 = x4 * (x2424 + x3927 - x3928 + x3929)
    x3931 = x2423 - x3927 + x3928
    x3932 = x1261 * x3136 + x2438 - x3146 * x8
    x3933 = x1261 * x3150 + x2446 + 4.0 * x3104 - x3152 * x8
    x3934 = x1261 * x3152 + x2459 + 4.0 * x3106 - x3153 * x8
    x3935 = x1261 * x3157 + x2468 + 3.0 * x3117 - x3159 * x8
    x3936 = x1261 * x3159 + x2481 + 3.0 * x3119 - x3160 * x8
    x3937 = 2.0 * x3130
    x3938 = x1261 * x3164 + x2489 - x3166 * x8 + x3937
    x3939 = 2.0 * x3132
    x3940 = x1261 * x3166 + x2501 - x3167 * x8 + x3939
    x3941 = x1261 * x3170 + x2508 + x3143 - x3172 * x8
    x3942 = x1261 * x3172 + x2519 + x3145 - x3173 * x8
    x3943 = x1261 * x3176 + x2526 - x3178 * x8
    x3944 = x1261 * x3178 + x2536 - x3179 * x8
    x3945 = x1261 * x3184
    x3946 = x1261 * x3214 + x3267
    x3947 = x4 * (x3247 - x3945 + x3946)
    x3948 = x3248 + x3945
    x3949 = x3946 * x5 - x3948 * x6
    x3950 = x3948 * x5
    x3951 = x1261 * x3187
    x3952 = x3252 + x3951
    x3953 = x3952 * x6
    x3954 = x3950 - x3953
    x3955 = x3947 + x3949 * x5 - x3954 * x6
    x3956 = x4 * (x3244 + x3948 - x3951)
    x3957 = x3954 * x5
    x3958 = x3952 * x5
    x3959 = x1261 * x3194
    x3960 = x3260 + x3959
    x3961 = x3960 * x6
    x3962 = x3958 - x3961
    x3963 = x3962 * x6
    x3964 = x3956 + x3957 - x3963
    x3965 = x38 * (x3949 - x3950 + x3953) + x3955 * x5 - x3964 * x6
    x3966 = x3964 * x5
    x3967 = x4 * (x3259 + x3952 - x3959)
    x3968 = x3962 * x5
    x3969 = x3960 * x5
    x3970 = x1261 * x3205
    x3971 = x3279 + x3970
    x3972 = x3971 * x6
    x3973 = x3969 - x3972
    x3974 = x3973 * x6
    x3975 = x3967 + x3968 - x3974
    x3976 = x3975 * x6
    x3977 = x38 * (x3954 - x3958 + x3961)
    x3978 = x3966 - x3976 + x3977
    x3979 = -x3956
    x3980 = x4 * (x3278 + x3960 - x3970)
    x3981 = x1261 * x3228 + x3291
    x3982 = -x3967
    x3983 = x1261 * x3249
    x3984 = x1261 * x3268 + x3350
    x3985 = x4 * (x3353 - x3983 + x3984)
    x3986 = x3340 + x3983
    x3987 = x3984 * x5 - x3986 * x6
    x3988 = x3986 * x5
    x3989 = x1261 * x3253
    x3990 = x3345 + x3989
    x3991 = x3990 * x6
    x3992 = x3988 - x3991
    x3993 = x3985 + x3987 * x5 - x3992 * x6
    x3994 = x4 * (x3365 + x3986 - x3989)
    x3995 = x3992 * x5
    x3996 = x3990 * x5
    x3997 = x1261 * x3261
    x3998 = x3360 + x3997
    x3999 = x3998 * x6
    x4000 = x3996 - x3999
    x4001 = x4000 * x6
    x4002 = x3994 + x3995 - x4001
    x4003 = x4 * (x3370 + x3990 - x3997)
    x4004 = x1261 * x3280 + x3372
    x4005 = x1261 * x3298
    x4006 = x1261 * x3314 + x3386
    x4007 = x4 * (x3375 - x4005 + x4006)
    x4008 = x3376 + x4005
    x4009 = x4006 * x5 - x4008 * x6
    x4010 = x4008 * x5
    x4011 = x1261 * x3301
    x4012 = x3381 + x4011
    x4013 = x4012 * x6
    x4014 = x4010 - x4013
    x4015 = x4007 + x4009 * x5 - x4014 * x6
    x4016 = x4 * (x3380 + x4008 - x4011)
    x4017 = x4014 * x5
    x4018 = x4012 * x5
    x4019 = x1261 * x3308
    x4020 = x3395 + x4019
    x4021 = x4020 * x6
    x4022 = x4018 - x4021
    x4023 = x4022 * x6
    x4024 = x4016 + x4017 - x4023
    x4025 = x4 * (x3394 + x4012 - x4019)
    x4026 = x1261 * x3325 + x3405
    x4027 = -x4016
    x4028 = x1261 * x3341
    x4029 = x1261 * x3351 + x3443
    x4030 = x4 * (x3459 - x4028 + x4029)
    x4031 = x3454 + x4028
    x4032 = x4029 * x5 - x4031 * x6
    x4033 = x4031 * x5
    x4034 = x1261 * x3346
    x4035 = x3448 + x4034
    x4036 = x4035 * x6
    x4037 = x4033 - x4036
    x4038 = x4 * (x3462 + x4031 - x4034)
    x4039 = x1261 * x3361 + x3464
    x4040 = x1261 * x3377
    x4041 = x1261 * x3387 + x3466
    x4042 = x4 * (x3480 - x4040 + x4041)
    x4043 = x3475 + x4040
    x4044 = x4041 * x5 - x4043 * x6
    x4045 = x4043 * x5
    x4046 = x1261 * x3382
    x4047 = x3470 + x4046
    x4048 = x4047 * x6
    x4049 = x4045 - x4048
    x4050 = x4 * (x3483 + x4043 - x4046)
    x4051 = x1261 * x3396 + x3485
    x4052 = x1261 * x3410
    x4053 = x1261 * x3420 + x3487
    x4054 = x4 * (x3495 - x4052 + x4053)
    x4055 = x3496 + x4052
    x4056 = x4053 * x5 - x4055 * x6
    x4057 = x4055 * x5
    x4058 = x1261 * x3415
    x4059 = x3491 + x4058
    x4060 = x4059 * x6
    x4061 = x4057 - x4060
    x4062 = x4 * (x3490 + x4055 - x4058)
    x4063 = x1261 * x3430 + x3504
    x4064 = x1261 * x3455
    x4065 = x1261 * x3444 + x3524
    x4066 = x4 * (x3528 - x4064 + x4065)
    x4067 = x3530 + x4064
    x4068 = x1261 * x3449 + x3532
    x4069 = x1261 * x3476
    x4070 = x1261 * x3467 + x3535
    x4071 = x4 * (x3539 - x4069 + x4070)
    x4072 = x3541 + x4069
    x4073 = x1261 * x3471 + x3543
    x4074 = x1261 * x3497
    x4075 = x1261 * x3488 + x3546
    x4076 = x4 * (x3549 - x4074 + x4075)
    x4077 = x3551 + x4074
    x4078 = x1261 * x3492 + x3553
    x4079 = x1261 * x3515
    x4080 = x1261 * x3506 + x3557
    x4081 = x4 * (x3555 - x4079 + x4080)
    x4082 = x3560 + x4079
    x4083 = x1261 * x3510 + x3562
    x4084 = x1261 * x3525 + 4.0 * x3460 - x3531 * x8
    x4085 = x1261 * x3531 + 4.0 * x3463 - x3533 * x8
    x4086 = x1261 * x3536 + 3.0 * x3481 - x3542 * x8
    x4087 = x1261 * x3542 + 3.0 * x3484 - x3544 * x8
    x4088 = x1261 * x3547 + 2.0 * x3501 - x3552 * x8
    x4089 = x1261 * x3552 + 2.0 * x3503 - x3554 * x8
    x4090 = x1261 * x3558 + x3519 - x3561 * x8
    x4091 = x1261 * x3561 + x3521 - x3563 * x8
    x4092 = x1261 * x3567 - x3569 * x8
    x4093 = x1261 * x3569 - x3570 * x8
    x4094 = x10 * x3187
    x4095 = x1852 * x3184
    x4096 = 2.0 * x1856
    x4097 = -x10 * x3184 + x1852 * x3214 + 2.0 * x1885
    x4098 = x4 * (x4094 - x4095 - x4096 + x4097)
    x4099 = -x4094 + x4095 + x4096
    x4100 = x4097 * x5 - x4099 * x6
    x4101 = x4099 * x5
    x4102 = x1852 * x3187
    x4103 = x10 * x3194
    x4104 = 2.0 * x1870
    x4105 = x4102 - x4103 + x4104
    x4106 = x4105 * x6
    x4107 = x4101 - x4106
    x4108 = x4098 + x4100 * x5 - x4107 * x6
    x4109 = x4 * (x4099 - x4102 + x4103 - x4104)
    x4110 = x4107 * x5
    x4111 = x4105 * x5
    x4112 = x1852 * x3194
    x4113 = x10 * x3205
    x4114 = 2.0 * x1894
    x4115 = x4112 - x4113 + x4114
    x4116 = x4115 * x6
    x4117 = x4111 - x4116
    x4118 = x4117 * x6
    x4119 = x4109 + x4110 - x4118
    x4120 = x38 * (x4100 - x4101 + x4106) + x4108 * x5 - x4119 * x6
    x4121 = x4119 * x5
    x4122 = x4 * (x4105 - x4112 + x4113 - x4114)
    x4123 = x4117 * x5
    x4124 = x4115 * x5
    x4125 = x1852 * x3205
    x4126 = x10 * x3228
    x4127 = 2.0 * x1918
    x4128 = x4125 - x4126 + x4127
    x4129 = x4128 * x6
    x4130 = x4124 - x4129
    x4131 = x4130 * x6
    x4132 = x4122 + x4123 - x4131
    x4133 = x4132 * x6
    x4134 = x38 * (x4107 - x4111 + x4116)
    x4135 = x4121 - x4133 + x4134
    x4136 = -x4109
    x4137 = x4 * (x4115 - x4125 + x4126 - x4127)
    x4138 = -x10 * x3242 + x1852 * x3228 + 2.0 * x1943
    x4139 = -x4122
    x4140 = x4105 * x8
    x4141 = x174 * x4099
    x4142 = -x4099 * x8
    x4143 = x174 * x4097 + x4142
    x4144 = x4 * (x4140 - x4141 + x4143)
    x4145 = -x4140
    x4146 = x4141 + x4145
    x4147 = x4143 * x5 - x4146 * x6
    x4148 = x4146 * x5
    x4149 = x174 * x4105
    x4150 = x4115 * x8
    x4151 = -x4150
    x4152 = x4149 + x4151
    x4153 = x4152 * x6
    x4154 = x4148 - x4153
    x4155 = x4144 + x4147 * x5 - x4154 * x6
    x4156 = x4 * (x4146 - x4149 + x4150)
    x4157 = x4154 * x5
    x4158 = x4152 * x5
    x4159 = x174 * x4115
    x4160 = x4128 * x8
    x4161 = -x4160
    x4162 = x4159 + x4161
    x4163 = x4162 * x6
    x4164 = x4158 - x4163
    x4165 = x4164 * x6
    x4166 = x4156 + x4157 - x4165
    x4167 = x4 * (x4152 - x4159 + x4160)
    x4168 = -x4138 * x8
    x4169 = x174 * x4128 + x4168
    x4170 = x10 * x3301
    x4171 = x1852 * x3298
    x4172 = -x10 * x3298 + x1852 * x3314 + x3215 + x3419
    x4173 = x4 * (x3219 + x3422 + x4170 - x4171 + x4172)
    x4174 = x3185 + x3409 - x4170 + x4171
    x4175 = x4172 * x5 - x4174 * x6
    x4176 = x4174 * x5
    x4177 = x1852 * x3301
    x4178 = x10 * x3308
    x4179 = x3200 + x3414 + x4177 - x4178
    x4180 = x4179 * x6
    x4181 = x4176 - x4180
    x4182 = x4173 + x4175 * x5 - x4181 * x6
    x4183 = x4 * (x3237 + x3434 + x4174 - x4177 + x4178)
    x4184 = x4181 * x5
    x4185 = x4179 * x5
    x4186 = x1852 * x3308
    x4187 = x10 * x3325
    x4188 = x3223 + x3429 + x4186 - x4187
    x4189 = x4188 * x6
    x4190 = x4185 - x4189
    x4191 = x4190 * x6
    x4192 = x4183 + x4184 - x4191
    x4193 = x4 * (x3243 + x3439 + x4179 - x4186 + x4187)
    x4194 = -x10 * x3336 + x1852 * x3325 + x3241 + x3441
    x4195 = -x4183
    x4196 = x174 * x4146
    x4197 = x4098 - x4146 * x8
    x4198 = x174 * x4143 + x4197
    x4199 = x4152 * x8
    x4200 = x4136 + x4199
    x4201 = x4 * (-x4196 + x4198 + x4200)
    x4202 = x4109 - x4199
    x4203 = x4196 + x4202
    x4204 = x4198 * x5 - x4203 * x6
    x4205 = x4203 * x5
    x4206 = x174 * x4152
    x4207 = x4162 * x8
    x4208 = x4122 - x4207
    x4209 = x4206 + x4208
    x4210 = x4209 * x6
    x4211 = x4205 - x4210
    x4212 = x4139 + x4207
    x4213 = x4 * (x4203 - x4206 + x4212)
    x4214 = x4137 - x4169 * x8
    x4215 = x174 * x4162 + x4214
    x4216 = x4179 * x8
    x4217 = x174 * x4174
    x4218 = -x4174 * x8
    x4219 = x174 * x4172 + x4218
    x4220 = x4 * (x4216 - x4217 + x4219)
    x4221 = -x4216
    x4222 = x4217 + x4221
    x4223 = x4219 * x5 - x4222 * x6
    x4224 = x4222 * x5
    x4225 = x174 * x4179
    x4226 = x4188 * x8
    x4227 = -x4226
    x4228 = x4225 + x4227
    x4229 = x4228 * x6
    x4230 = x4224 - x4229
    x4231 = x4 * (x4222 - x4225 + x4226)
    x4232 = -x4194 * x8
    x4233 = x174 * x4188 + x4232
    x4234 = x10 * x3415
    x4235 = x1852 * x3410
    x4236 = 2.0 * x2256
    x4237 = 2.0 * x3299
    x4238 = -x10 * x3410 + x1852 * x3420 + 2.0 * x2247 + 2.0 * x3315
    x4239 = x4 * (x4234 - x4235 - x4236 - x4237 + x4238)
    x4240 = -x4234 + x4235 + x4236 + x4237
    x4241 = x4238 * x5 - x4240 * x6
    x4242 = x4240 * x5
    x4243 = x1852 * x3415
    x4244 = x10 * x3430
    x4245 = 2.0 * x2265
    x4246 = 2.0 * x3320
    x4247 = x4243 - x4244 + x4245 + x4246
    x4248 = x4247 * x6
    x4249 = x4242 - x4248
    x4250 = x4 * (x4240 - x4243 + x4244 - x4245 - x4246)
    x4251 = -x10 * x3442 + x1852 * x3430 + 2.0 * x2282 + 2.0 * x3335
    x4252 = x174 * x4203
    x4253 = 2.0 * x4144 - x4203 * x8
    x4254 = x174 * x4198 + x4253
    x4255 = x4209 * x8
    x4256 = 2.0 * x4156
    x4257 = x4255 - x4256
    x4258 = x4 * (-x4252 + x4254 + x4257)
    x4259 = -x4255 + x4256
    x4260 = x4252 + x4259
    x4261 = 2.0 * x4167 - x4215 * x8
    x4262 = x174 * x4209 + x4261
    x4263 = x174 * x4222
    x4264 = x4173 - x4222 * x8
    x4265 = x174 * x4219 + x4264
    x4266 = x4228 * x8
    x4267 = x4195 + x4266
    x4268 = x4 * (-x4263 + x4265 + x4267)
    x4269 = x4183 - x4266
    x4270 = x4263 + x4269
    x4271 = x4193 - x4233 * x8
    x4272 = x174 * x4228 + x4271
    x4273 = x4247 * x8
    x4274 = x174 * x4240
    x4275 = -x4240 * x8
    x4276 = x174 * x4238 + x4275
    x4277 = x4 * (x4273 - x4274 + x4276)
    x4278 = -x4273
    x4279 = x4274 + x4278
    x4280 = -x4251 * x8
    x4281 = x174 * x4247 + x4280
    x4282 = x10 * x3510
    x4283 = x1852 * x3515
    x4284 = 3.0 * x3435
    x4285 = -x4284
    x4286 = 2.0 * x2423
    x4287 = 3.0 * x3423
    x4288 = -x10 * x3515 + x1852 * x3506 + 2.0 * x2412 + x4287
    x4289 = x4 * (x4282 - x4283 + x4285 - x4286 + x4288)
    x4290 = -x4282 + x4283 + x4284 + x4286
    x4291 = 3.0 * x3440
    x4292 = -x10 * x3522 + x1852 * x3510 + 2.0 * x2438 + x4291
    x4293 = 3.0 * x4201 - x4260 * x8
    x4294 = x174 * x4254 + x4293
    x4295 = 3.0 * x4213 - x4262 * x8
    x4296 = x174 * x4260 + x4295
    x4297 = 2.0 * x4220 - x4270 * x8
    x4298 = x174 * x4265 + x4297
    x4299 = 2.0 * x4231 - x4272 * x8
    x4300 = x174 * x4270 + x4299
    x4301 = x4239 - x4279 * x8
    x4302 = x174 * x4276 + x4301
    x4303 = x4250 - x4281 * x8
    x4304 = x174 * x4279 + x4303
    x4305 = -x4290 * x8
    x4306 = x174 * x4288 + x4305
    x4307 = -x4292 * x8
    x4308 = x174 * x4290 + x4307
    x4309 = -x10 * x3569 + x1852 * x3567 + 2.0 * x2526 + 4.0 * x3519
    x4310 = -x10 * x3570 + x1852 * x3569 + 2.0 * x2536 + 4.0 * x3521
    x4311 = x3582 * x8
    x4312 = x1261 * x3576
    x4313 = 3.0 * x2543
    x4314 = x1261 * x3574 + 3.0 * x2573 - x3576 * x8
    x4315 = -x4311 + x4312 + x4313
    x4316 = x4314 * x5 - x4315 * x6
    x4317 = x4315 * x5
    x4318 = x1261 * x3582
    x4319 = x3592 * x8
    x4320 = 3.0 * x2558
    x4321 = x4318 - x4319 + x4320
    x4322 = x4321 * x6
    x4323 = x4317 - x4322
    x4324 = x4 * (x4311 - x4312 - x4313 + x4314) + x4316 * x5 - x4323 * x6
    x4325 = x4 * (x4315 - x4318 + x4319 - x4320)
    x4326 = x4323 * x5
    x4327 = x4321 * x5
    x4328 = x1261 * x3592
    x4329 = x3605 * x8
    x4330 = 3.0 * x2581
    x4331 = x4328 - x4329 + x4330
    x4332 = x4331 * x6
    x4333 = x4327 - x4332
    x4334 = x4333 * x6
    x4335 = x4325 + x4326 - x4334
    x4336 = x3626 * x8
    x4337 = x1261 * x3621
    x4338 = 3.0 * x2607
    x4339 = x1261 * x3619 + 3.0 * x2623 + x3575 - x3621 * x8
    x4340 = x3586 - x4336 + x4337 + x4338
    x4341 = x4339 * x5 - x4340 * x6
    x4342 = x4340 * x5
    x4343 = x1261 * x3626
    x4344 = x3635 * x8
    x4345 = 3.0 * x2628
    x4346 = x3599 + x4343 - x4344 + x4345
    x4347 = x4346 * x6
    x4348 = x4342 - x4347
    x4349 = x3653 * x8
    x4350 = x1261 * x3647
    x4351 = 3.0 * x2649
    x4352 = x1261 * x3645 + 3.0 * x2665 - x3647 * x8
    x4353 = -x4349 + x4350 + x4351
    x4354 = x4352 * x5 - x4353 * x6
    x4355 = x4353 * x5
    x4356 = x1261 * x3653
    x4357 = x3663 * x8
    x4358 = 3.0 * x2670
    x4359 = x4356 - x4357 + x4358
    x4360 = x4359 * x6
    x4361 = x4355 - x4360
    x4362 = x3684 * x8
    x4363 = x1261 * x3677
    x4364 = 2.0 * x3630
    x4365 = x1261 * x3675 + 2.0 * x3620 - x3677 * x8 + x3724
    x4366 = x3721 - x4362 + x4363 + x4364
    x4367 = x3698 * x8
    x4368 = x1261 * x3693
    x4369 = 3.0 * x2746
    x4370 = x1261 * x3691 + 3.0 * x2736 + x3646 - x3693 * x8
    x4371 = x3657 - x4367 + x4368 + x4369
    x4372 = x3714 * x8
    x4373 = x1261 * x3708
    x4374 = 3.0 * x2774
    x4375 = x1261 * x3706 + 3.0 * x2764 - x3708 * x8
    x4376 = -x4372 + x4373 + x4374
    x4377 = x3772 * x8
    x4378 = x1261 * x3767
    x4379 = 2.0 * x2887
    x4380 = x1261 * x3765 + 2.0 * x2915 - x3767 * x8
    x4381 = -x4377 + x4378 + x4379
    x4382 = x4380 * x5 - x4381 * x6
    x4383 = x4381 * x5
    x4384 = x1261 * x3772
    x4385 = x3781 * x8
    x4386 = 2.0 * x2901
    x4387 = x4384 - x4385 + x4386
    x4388 = x4387 * x6
    x4389 = x4383 - x4388
    x4390 = x4 * (x4377 - x4378 - x4379 + x4380) + x4382 * x5 - x4389 * x6
    x4391 = x4 * (x4381 - x4384 + x4385 - x4386)
    x4392 = x4389 * x5
    x4393 = x4387 * x5
    x4394 = x1261 * x3781
    x4395 = x3793 * x8
    x4396 = 2.0 * x2923
    x4397 = x4394 - x4395 + x4396
    x4398 = x4397 * x6
    x4399 = x4393 - x4398
    x4400 = x4399 * x6
    x4401 = x4391 + x4392 - x4400
    x4402 = x3814 * x8
    x4403 = x1261 * x3809
    x4404 = x1261 * x3807 + x3766 - x3809 * x8 + x3860
    x4405 = x3776 + x3858 - x4402 + x4403
    x4406 = x4404 * x5 - x4405 * x6
    x4407 = x4405 * x5
    x4408 = x1261 * x3814
    x4409 = x3823 * x8
    x4410 = x3788 + x3868 + x4408 - x4409
    x4411 = x4410 * x6
    x4412 = x4407 - x4411
    x4413 = x3839 * x8
    x4414 = x1261 * x3834
    x4415 = 2.0 * x2984
    x4416 = x1261 * x3832 + 2.0 * x2999 - x3834 * x8
    x4417 = -x4413 + x4414 + x4415
    x4418 = x4416 * x5 - x4417 * x6
    x4419 = x4417 * x5
    x4420 = x1261 * x3839
    x4421 = x3848 * x8
    x4422 = 2.0 * x3004
    x4423 = x4420 - x4421 + x4422
    x4424 = x4423 * x6
    x4425 = x4419 - x4424
    x4426 = x3869 * x8
    x4427 = x1261 * x3863
    x4428 = 2.0 * x3040
    x4429 = 2.0 * x3818
    x4430 = x1261 * x3861 + 2.0 * x3031 + 2.0 * x3808 - x3863 * x8
    x4431 = -x4426 + x4427 + x4428 + x4429
    x4432 = x3885 * x8
    x4433 = x1261 * x3880
    x4434 = x1261 * x3878 + x3833 - x3880 * x8 + x3915
    x4435 = x3843 + x3913 - x4432 + x4433
    x4436 = x3899 * x8
    x4437 = x1261 * x3894
    x4438 = 2.0 * x3089
    x4439 = x1261 * x3892 + 2.0 * x3080 - x3894 * x8
    x4440 = -x4436 + x4437 + x4438
    x4441 = x3952 * x8
    x4442 = x1261 * x3948
    x4443 = x1261 * x3946 + x3215 - x3948 * x8
    x4444 = x3185 - x4441 + x4442
    x4445 = x4443 * x5 - x4444 * x6
    x4446 = x4444 * x5
    x4447 = x1261 * x3952
    x4448 = x3960 * x8
    x4449 = x3200 + x4447 - x4448
    x4450 = x4449 * x6
    x4451 = x4446 - x4450
    x4452 = x4 * (x3219 + x4441 - x4442 + x4443) + x4445 * x5 - x4451 * x6
    x4453 = x4 * (x3237 + x4444 - x4447 + x4448)
    x4454 = x4451 * x5
    x4455 = x4449 * x5
    x4456 = x1261 * x3960
    x4457 = x3971 * x8
    x4458 = x3223 + x4456 - x4457
    x4459 = x4458 * x6
    x4460 = x4455 - x4459
    x4461 = x4460 * x6
    x4462 = x4453 + x4454 - x4461
    x4463 = x3990 * x8
    x4464 = x1261 * x3986
    x4465 = x1261 * x3984 + x3269 + x3947 - x3986 * x8
    x4466 = x3250 + x3956 - x4463 + x4464
    x4467 = x4465 * x5 - x4466 * x6
    x4468 = x4466 * x5
    x4469 = x1261 * x3990
    x4470 = x3998 * x8
    x4471 = x3274 + x3967 + x4469 - x4470
    x4472 = x4471 * x6
    x4473 = x4468 - x4472
    x4474 = x4012 * x8
    x4475 = x1261 * x4008
    x4476 = x1261 * x4006 + x3315 - x4008 * x8
    x4477 = x3299 - x4474 + x4475
    x4478 = x4476 * x5 - x4477 * x6
    x4479 = x4477 * x5
    x4480 = x1261 * x4012
    x4481 = x4020 * x8
    x4482 = x3320 + x4480 - x4481
    x4483 = x4482 * x6
    x4484 = x4479 - x4483
    x4485 = x4035 * x8
    x4486 = x1261 * x4031
    x4487 = 2.0 * x3994
    x4488 = x1261 * x4029 + x3354 + 2.0 * x3985 - x4031 * x8
    x4489 = x3366 - x4485 + x4486 + x4487
    x4490 = x4047 * x8
    x4491 = x1261 * x4043
    x4492 = x1261 * x4041 + x3389 + x4007 - x4043 * x8
    x4493 = x3400 + x4016 - x4490 + x4491
    x4494 = x4059 * x8
    x4495 = x1261 * x4055
    x4496 = x1261 * x4053 + x3423 - x4055 * x8
    x4497 = x3435 - x4494 + x4495
    x4498 = x1261 * x4099
    x4499 = x1261 * x4097 + x4142
    x4500 = x4145 + x4498
    x4501 = x4499 * x5 - x4500 * x6
    x4502 = x4500 * x5
    x4503 = x1261 * x4105
    x4504 = x4151 + x4503
    x4505 = x4504 * x6
    x4506 = x4502 - x4505
    x4507 = x4 * (x4140 - x4498 + x4499) + x4501 * x5 - x4506 * x6
    x4508 = x4 * (x4150 + x4500 - x4503)
    x4509 = x4506 * x5
    x4510 = x4504 * x5
    x4511 = x1261 * x4115
    x4512 = x4161 + x4511
    x4513 = x4512 * x6
    x4514 = x4510 - x4513
    x4515 = x4514 * x6
    x4516 = x4508 + x4509 - x4515
    x4517 = x1261 * x4146
    x4518 = x1261 * x4143 + x4197
    x4519 = x4202 + x4517
    x4520 = x4518 * x5 - x4519 * x6
    x4521 = x4519 * x5
    x4522 = x1261 * x4152
    x4523 = x4208 + x4522
    x4524 = x4523 * x6
    x4525 = x4521 - x4524
    x4526 = x1261 * x4174
    x4527 = x1261 * x4172 + x4218
    x4528 = x4221 + x4526
    x4529 = x4527 * x5 - x4528 * x6
    x4530 = x4528 * x5
    x4531 = x1261 * x4179
    x4532 = x4227 + x4531
    x4533 = x4532 * x6
    x4534 = x4530 - x4533
    x4535 = x1261 * x4203
    x4536 = x1261 * x4198 + x4253
    x4537 = x4259 + x4535
    x4538 = x1261 * x4222
    x4539 = x1261 * x4219 + x4264
    x4540 = x4269 + x4538
    x4541 = x1261 * x4240
    x4542 = x1261 * x4238 + x4275
    x4543 = x4278 + x4541
    x4544 = x10 * x4105
    x4545 = x1852 * x4099
    x4546 = 3.0 * x3185
    x4547 = -x10 * x4099 + x1852 * x4097 + 3.0 * x3215
    x4548 = x4 * (x4544 - x4545 - x4546 + x4547)
    x4549 = -x4544 + x4545 + x4546
    x4550 = x4547 * x5 - x4549 * x6
    x4551 = x4549 * x5
    x4552 = x1852 * x4105
    x4553 = x10 * x4115
    x4554 = 3.0 * x3200
    x4555 = x4552 - x4553 + x4554
    x4556 = x4555 * x6
    x4557 = x4551 - x4556
    x4558 = x4548 + x4550 * x5 - x4557 * x6
    x4559 = x4 * (x4549 - x4552 + x4553 - x4554)
    x4560 = x4557 * x5
    x4561 = x4555 * x5
    x4562 = x1852 * x4115
    x4563 = x10 * x4128
    x4564 = 3.0 * x3223
    x4565 = x4562 - x4563 + x4564
    x4566 = x4565 * x6
    x4567 = x4561 - x4566
    x4568 = x4567 * x6
    x4569 = x4559 + x4560 - x4568
    x4570 = x4 * (x4555 - x4562 + x4563 - x4564)
    x4571 = -x10 * x4138 + x1852 * x4128 + 3.0 * x3241
    x4572 = -x4559
    x4573 = x4555 * x8
    x4574 = x174 * x4549
    x4575 = x174 * x4547 - x4549 * x8
    x4576 = x4 * (x4573 - x4574 + x4575)
    x4577 = -x4573 + x4574
    x4578 = x4575 * x5 - x4577 * x6
    x4579 = x4577 * x5
    x4580 = x174 * x4555
    x4581 = x4565 * x8
    x4582 = x4580 - x4581
    x4583 = x4582 * x6
    x4584 = x4579 - x4583
    x4585 = x4 * (x4577 - x4580 + x4581)
    x4586 = x174 * x4565 - x4571 * x8
    x4587 = x10 * x4179
    x4588 = x1852 * x4174
    x4589 = 3.0 * x3299
    x4590 = -x10 * x4174 + x1852 * x4172 + 3.0 * x3315 + x4098
    x4591 = x4 * (x4136 + x4587 - x4588 - x4589 + x4590)
    x4592 = x4109 - x4587 + x4588 + x4589
    x4593 = x4590 * x5 - x4592 * x6
    x4594 = x4592 * x5
    x4595 = x1852 * x4179
    x4596 = x10 * x4188
    x4597 = 3.0 * x3320
    x4598 = x4122 + x4595 - x4596 + x4597
    x4599 = x4598 * x6
    x4600 = x4594 - x4599
    x4601 = x4 * (x4139 + x4592 - x4595 + x4596 - x4597)
    x4602 = -x10 * x4194 + x1852 * x4188 + 3.0 * x3335 + x4137
    x4603 = x4582 * x8
    x4604 = x174 * x4577
    x4605 = x174 * x4575 + x4548 - x4577 * x8
    x4606 = x4 * (x4572 + x4603 - x4604 + x4605)
    x4607 = x4559 - x4603 + x4604
    x4608 = x174 * x4582 + x4570 - x4586 * x8
    x4609 = x4598 * x8
    x4610 = x174 * x4592
    x4611 = x174 * x4590 - x4592 * x8
    x4612 = x4 * (x4609 - x4610 + x4611)
    x4613 = -x4609 + x4610
    x4614 = x174 * x4598 - x4602 * x8
    x4615 = x10 * x4247
    x4616 = x1852 * x4240
    x4617 = 2.0 * x4183
    x4618 = -x10 * x4240 + x1852 * x4238 + 2.0 * x4173 + x4287
    x4619 = x4 * (x4285 + x4615 - x4616 - x4617 + x4618)
    x4620 = x4284 - x4615 + x4616 + x4617
    x4621 = -x10 * x4251 + x1852 * x4247 + 2.0 * x4193 + x4291
    x4622 = x174 * x4605 + 2.0 * x4576 - x4607 * x8
    x4623 = x174 * x4607 + 2.0 * x4585 - x4608 * x8
    x4624 = x174 * x4611 + x4591 - x4613 * x8
    x4625 = x174 * x4613 + x4601 - x4614 * x8
    x4626 = x174 * x4618 - x4620 * x8
    x4627 = x174 * x4620 - x4621 * x8
    x4628 = -x10 * x4290 + x1852 * x4288 + 3.0 * x3519 + 3.0 * x4239
    x4629 = -x10 * x4292 + x1852 * x4290 + 3.0 * x3521 + 3.0 * x4250

    # 225 item(s)
    result[0, 0] = numpy.sum(
        x171
        * (
            x109
            * (
                x108 * x38
                - x130
                * (
                    x150 * x6
                    - x152 * x3
                    + x165
                    - x166
                    + x169
                    - x38 * (-x168 + x3 * x44 + x46 + x51)
                    - x52
                    + x71
                )
                + x130
                * (
                    x162 * x3
                    - x165
                    + x166
                    - x169
                    + x38 * (x167 - x3 * x33 + x35 + x96)
                    - x6 * (x161 + x92)
                    - x71
                    + x97
                )
                - x142 * x38
                + x153 * x3
                - x164 * x3
                - x164 * x6
                + x6 * (x136 - x159 + x160 + x163)
            )
            - x130 * (-x110 + x112 + x113 - x139 + x143)
            - x3
            * (
                x109 * (x103 - x146 + x147 + x151 - x153)
                - x110 * x3
                + x144 * x6
                - x38 * (x101 - x104 + x78 - x79)
            )
            + x6
            * (
                x109 * (x136 - x159 + x160 + x163 - x164)
                - x144 * x3
                - x38 * (x105 - x114 + x134 - x137)
                + x6
                * (
                    x109 * (x145 + x154 - x158)
                    + x138 * x3
                    + x4 * (x100 - x115 + x154)
                    - x6
                    * (
                        -x109 * (x116 - x125 + x126 - x93)
                        + x133 * x3
                        - x6
                        * (
                            x127 * x5
                            + x130 * (-x118 + x123 + x157 + x89)
                            - x6
                            * (
                                x124 * x5
                                + x38 * (x120 - x121 + x87)
                                - x6 * (x122 * x5 + x155 - x6 * (x119 * x5 - x156 * x6))
                            )
                        )
                    )
                )
            )
        )
    )
    result[0, 1] = numpy.sum(
        -x318
        * (
            x130 * (-x259 + x261 + x262 - x286 + x290)
            - x130
            * (
                x258 * x38
                - x289 * x38
                + x3 * x300
                - x3 * x312
                - x312 * x6
                + x38
                * (
                    x206
                    - x222
                    - x297 * x6
                    + x299 * x3
                    - x315
                    - x316
                    + x317
                    + x4 * (x182 + x200 * x3 + x201 - x313)
                )
                - x38
                * (
                    x222
                    - x246
                    - x3 * x310
                    + x315
                    + x316
                    - x317
                    - x4 * (-x191 * x3 + x193 + x216 + x314)
                    + x6 * (x243 + x309)
                )
                + x6 * (x283 - x307 + x308 + x311)
            )
            + x3
            * (
                x130 * (x253 - x293 + x294 + x298 - x300)
                - x259 * x3
                + x291 * x6
                - x38 * (x230 - x231 + x250 - x254)
            )
            - x6
            * (
                x130 * (x283 - x307 + x308 + x311 - x312)
                - x291 * x3
                - x38 * (x255 - x263 + x280 - x284)
                + x6
                * (
                    x130 * (x292 + x301 - x306)
                    + x285 * x3
                    + x4 * (x249 - x264 + x301)
                    - x6
                    * (
                        x130 * (x244 - x266 + x274 + x305)
                        + x279 * x3
                        - x6
                        * (
                            x275 * x5
                            + x38 * (x241 - x267 + x272)
                            - x6 * (x273 * x5 + x302 - x6 * (x271 * x5 - x304 * x6))
                        )
                    )
                )
            )
        )
    )
    result[0, 2] = numpy.sum(
        -x318
        * (
            x130 * (-x406 + x408 + x409 - x433 + x437)
            - x130
            * (
                x3 * x447
                - x3 * x459
                + x38 * x405
                - x38 * x436
                + x38
                * (
                    x3 * x446
                    + x353
                    - x369
                    + x4 * (x3 * x347 + x329 + x348 - x460)
                    - x444 * x6
                    - x462
                    - x463
                    + x464
                )
                + x38
                * (
                    x3 * x457
                    - x369
                    + x393
                    + x4 * (-x3 * x338 + x340 + x363 + x461)
                    - x462
                    - x463
                    + x464
                    - x6 * (x390 + x456)
                )
                - x459 * x6
                + x6 * (x430 - x454 + x455 + x458)
            )
            + x3
            * (
                x130 * (x400 - x440 + x441 + x445 - x447)
                - x3 * x406
                - x38 * (x377 - x378 + x397 - x401)
                + x438 * x6
            )
            - x6
            * (
                x130 * (x430 - x454 + x455 + x458 - x459)
                - x3 * x438
                - x38 * (x402 - x410 + x427 - x431)
                + x6
                * (
                    x130 * (x439 + x448 - x453)
                    + x3 * x432
                    + x4 * (x396 - x411 + x448)
                    - x6
                    * (
                        x130 * (x391 - x413 + x421 + x452)
                        + x3 * x426
                        - x6
                        * (
                            x38 * (x388 - x414 + x419)
                            + x422 * x5
                            - x6 * (x420 * x5 + x449 - x6 * (x418 * x5 - x451 * x6))
                        )
                    )
                )
            )
        )
    )
    result[0, 3] = numpy.sum(
        x578
        * (
            x130 * (x525 - x528 - x529 + x548 + x577)
            + x3
            * (
                x3 * x525
                + x38 * (x502 - x503 + x518 - x520)
                - x38 * (x519 - x554 + x558 + x559 - x561)
                - x552 * x6
            )
            + x38
            * (
                x3 * x561
                - x3 * x574
                + x4 * (x3 * x560 + x483 + x497 - x557 * x6 + x575 - x576)
                - x4 * (-x3 * x571 + x496 + x526 - x575 + x576 + x6 * (x513 + x570))
                + x524
                - x574 * x6
                + x577
                + x6 * (x545 - x569 + x572 + x573)
            )
            - x6
            * (
                x3 * x552
                + x38 * (x521 - x530 + x544 - x546)
                - x38 * (x545 - x569 + x572 + x573 - x574)
                - x6
                * (
                    x3 * x547
                    + x38 * (x553 + x563 - x568)
                    + x4 * (x517 - x531 + x563)
                    - x6
                    * (
                        x3 * x543
                        + x38 * (x514 - x534 + x539)
                        - x6 * (x5 * x540 + x565 - x6 * (x5 * x538 - x567 * x6))
                    )
                )
            )
        )
    )
    result[0, 4] = numpy.sum(
        x687
        * (
            x130 * (x636 - x639 - x640 + x658 + x686)
            + x3
            * (
                x3 * x636
                + x38 * (x614 - x615 + x629 - x631)
                - x38 * (x630 - x664 + x668 + x669 - x671)
                - x6 * x662
            )
            + x38
            * (
                x3 * x671
                - x3 * x683
                + x4 * (x3 * x670 + x596 - x6 * x667 + x609 + x684 - x685)
                - x4 * (-x3 * x680 + x6 * (x624 + x679) + x608 + x637 - x684 + x685)
                - x6 * x683
                + x6 * (x655 - x678 + x681 + x682)
                + x635
                + x686
            )
            - x6
            * (
                x3 * x662
                + x38 * (x632 - x641 + x654 - x656)
                - x38 * (x655 - x678 + x681 + x682 - x683)
                - x6
                * (
                    x3 * x657
                    + x38 * (x663 + x673 - x677)
                    + x4 * (x628 - x642 + x673)
                    - x6
                    * (
                        x3 * x653
                        + x38 * (x625 - x644 + x649)
                        - x6 * (x5 * x650 - x6 * (x5 * x648 - x6 * x676) + x674)
                    )
                )
            )
        )
    )
    result[0, 5] = numpy.sum(
        x578
        * (
            x130 * (x748 - x751 - x752 + x771 + x800)
            + x3
            * (
                x3 * x748
                + x38 * (x725 - x726 + x741 - x743)
                - x38 * (x742 - x777 + x781 + x782 - x784)
                - x6 * x775
            )
            + x38
            * (
                x3 * x784
                - x3 * x797
                + x4 * (x3 * x783 - x6 * x780 + x706 + x720 + x798 - x799)
                - x4 * (-x3 * x794 + x6 * (x736 + x793) + x719 + x749 - x798 + x799)
                - x6 * x797
                + x6 * (x768 - x792 + x795 + x796)
                + x747
                + x800
            )
            - x6
            * (
                x3 * x775
                + x38 * (x744 - x753 + x767 - x769)
                - x38 * (x768 - x792 + x795 + x796 - x797)
                - x6
                * (
                    x3 * x770
                    + x38 * (x776 + x786 - x791)
                    + x4 * (x740 - x754 + x786)
                    - x6
                    * (
                        x3 * x766
                        + x38 * (x737 - x757 + x762)
                        - x6 * (x5 * x763 - x6 * (x5 * x761 - x6 * x790) + x788)
                    )
                )
            )
        )
    )
    result[0, 6] = numpy.sum(
        -x318
        * (
            x130 * (x848 - x860 + x861 + x862 - x868)
            - x3
            * (
                x3 * x860
                + x38 * (x850 + x854 - x855 + x858)
                + x4 * (x821 - x845 + x846 + x850)
                - x6 * x869
            )
            - x4
            * (
                x3 * x821
                - x3 * x847
                - x6 * x847
                + x6 * (x826 + x833 - x842)
                - 2.0 * x848
                + 2.0 * x849
            )
            + x6
            * (
                x3 * x869
                + x38 * (x859 - x863 + x866 + x870)
                + x4 * (-x833 + x842 + x847 + x870)
                - x6
                * (
                    x3 * x867
                    + x4 * (x832 - x834 + x840)
                    + x4 * (x840 + x857 - x864)
                    - x6 * (x3 * x865 - x6 * (x5 * x839 - x6 * x874) + x872)
                )
            )
        )
    )
    result[0, 7] = numpy.sum(
        -x687
        * (
            x130 * (x918 - x930 + x931 + x932 - x938)
            - x3
            * (
                x3 * x930
                + x38 * (x920 + x924 - x925 + x928)
                + x4 * (x893 - x915 + x916 + x920)
                - x6 * x939
            )
            - x4
            * (
                x3 * x893
                - x3 * x917
                - x6 * x917
                + x6 * (x897 + x904 - x912)
                - 2.0 * x918
                + 2.0 * x919
            )
            + x6
            * (
                x3 * x939
                + x38 * (x929 - x933 + x936 + x940)
                + x4 * (-x904 + x912 + x917 + x940)
                - x6
                * (
                    x3 * x937
                    + x4 * (x903 - x905 + x910)
                    + x4 * (x910 + x927 - x934)
                    - x6 * (x3 * x935 - x6 * (x5 * x909 - x6 * x944) + x942)
                )
            )
        )
    )
    result[0, 8] = numpy.sum(
        -x687
        * (
            x130 * (-x1005 + x985 - x997 + x998 + x999)
            - x3
            * (
                -x1006 * x6
                + x3 * x997
                + x38 * (x987 + x991 - x992 + x995)
                + x4 * (x962 - x982 + x983 + x987)
            )
            - x4
            * (
                x3 * x962
                - x3 * x984
                - x6 * x984
                + x6 * (x965 + x972 - x980)
                - 2.0 * x985
                + 2.0 * x986
            )
            + x6
            * (
                x1006 * x3
                + x38 * (-x1000 + x1003 + x1007 + x996)
                + x4 * (x1007 - x972 + x980 + x984)
                - x6
                * (
                    x1004 * x3
                    + x4 * (-x1001 + x978 + x994)
                    + x4 * (x971 - x973 + x978)
                    - x6 * (x1002 * x3 + x1008 + x6 * (x1010 * x6 - x5 * x977))
                )
            )
        )
    )
    result[0, 9] = numpy.sum(
        -x318
        * (
            x130 * (x1058 - x1070 + x1071 + x1072 - x1078)
            - x3
            * (
                x1070 * x3
                - x1079 * x6
                + x38 * (x1060 + x1064 - x1065 + x1068)
                + x4 * (x1031 - x1055 + x1056 + x1060)
            )
            - x4
            * (
                x1031 * x3
                - x1057 * x3
                - x1057 * x6
                - 2.0 * x1058
                + 2.0 * x1059
                + x6 * (x1036 + x1043 - x1052)
            )
            + x6
            * (
                x1079 * x3
                + x38 * (x1069 - x1073 + x1076 + x1080)
                + x4 * (-x1043 + x1052 + x1057 + x1080)
                - x6
                * (
                    x1077 * x3
                    + x4 * (x1042 - x1044 + x1050)
                    + x4 * (x1050 + x1067 - x1074)
                    - x6 * (x1075 * x3 + x1082 - x6 * (x1049 * x5 - x1084 * x6))
                )
            )
        )
    )
    result[0, 10] = numpy.sum(
        x171
        * (
            x130 * (x1103 - x1106 + x1115 + x1121)
            + x3 * (x1103 * x3 - x1116 * x6 + x38 * (x1094 - x1095 + x1101))
            - x6
            * (
                x1116 * x3
                + x38 * (x1102 - x1107 + x1113)
                - x6 * (x1114 * x3 + x1118 - x6 * (x1112 * x3 - x1120 * x6))
            )
        )
    )
    result[0, 11] = numpy.sum(
        x318
        * (
            x130 * (x1140 - x1143 + x1152 + x1158)
            + x3 * (x1140 * x3 - x1153 * x6 + x38 * (x1131 - x1132 + x1138))
            - x6
            * (
                x1153 * x3
                + x38 * (x1139 - x1144 + x1150)
                - x6 * (x1151 * x3 + x1155 - x6 * (x1149 * x3 - x1157 * x6))
            )
        )
    )
    result[0, 12] = numpy.sum(
        x578
        * (
            x130 * (x1175 - x1178 + x1186 + x1192)
            + x3 * (x1175 * x3 - x1187 * x6 + x38 * (x1167 - x1168 + x1173))
            - x6
            * (
                x1187 * x3
                + x38 * (x1174 - x1179 + x1184)
                - x6 * (x1185 * x3 + x1189 - x6 * (x1183 * x3 - x1191 * x6))
            )
        )
    )
    result[0, 13] = numpy.sum(
        x318
        * (
            x130 * (x1208 - x1210 + x1218 + x1223)
            + x3 * (x1208 * x3 - x1219 * x6 + x38 * (x1200 - x1201 + x1206))
            - x6
            * (
                x1219 * x3
                + x38 * (x1207 - x1211 + x1216)
                - x6 * (x1217 * x3 + x1220 - x6 * (x1215 * x3 - x1222 * x6))
            )
        )
    )
    result[0, 14] = numpy.sum(
        x171
        * (
            x130 * (x1242 - x1245 + x1254 + x1260)
            + x3 * (x1242 * x3 - x1255 * x6 + x38 * (x1233 - x1234 + x1240))
            - x6
            * (
                x1255 * x3
                + x38 * (x1241 - x1246 + x1252)
                - x6 * (x1253 * x3 + x1257 - x6 * (x1251 * x3 - x1259 * x6))
            )
        )
    )
    result[1, 0] = numpy.sum(
        -x318
        * (
            x109
            * (
                x130 * (-x1287 * x3 + x1290 + x1345 + x1357)
                - x130 * (x1296 * x3 + x1297 + x1321 - x1357)
                - x1322
                + x1346
                - x1350 * x3
                + x1356 * x3
                + x1356 * x6
                - x6 * (x1342 + x1355)
            )
            - x3
            * (
                x109 * (x1318 - x1349 + x1350)
                + x1323 * x3
                - x1348 * x6
                + x4 * (-x1292 + x1301 + x1318)
            )
            - x38 * (x1323 - x1324 + x1344 - x1347)
            + x6
            * (
                x109 * (x1351 - x1355 + x1356)
                + x1348 * x3
                + x4 * (x1320 - x1325 + x1351)
                - x6
                * (
                    x109 * (x1314 - x1326 + x1336 - x1337)
                    + x1343 * x3
                    - x6
                    * (
                        x130 * (x1311 - x1328 + x1334 + x1354)
                        + x1338 * x5
                        - x6
                        * (
                            x1335 * x5
                            + x38 * (x1309 - x1329 + x1332)
                            - x6 * (x1333 * x5 + x1352 - x6 * (x1331 * x5 - x1353 * x6))
                        )
                    )
                )
            )
        )
    )
    result[1, 1] = numpy.sum(
        x1428
        * (
            x130
            * (
                x1398
                - x1417
                + x1421 * x3
                - x1426 * x3
                - x1426 * x6
                - x38 * (-x1371 * x3 + x1373 + x1416 + x1427)
                + x38 * (x1377 * x3 + x1379 + x1397 - x1427)
                + x6 * (x1413 + x1425)
            )
            + x3
            * (
                x130 * (x1394 - x1420 + x1421)
                + x1399 * x3
                - x1419 * x6
                + x4 * (-x1375 + x1382 + x1394)
            )
            + x38 * (x1399 - x1400 + x1415 - x1418)
            - x6
            * (
                x130 * (x1422 - x1425 + x1426)
                + x1419 * x3
                + x4 * (x1396 - x1401 + x1422)
                - x6
                * (
                    x130 * (x1391 - x1402 - x1403 + x1409)
                    + x1414 * x3
                    - x6
                    * (
                        x1410 * x5
                        + x38 * (x1389 - x1404 + x1407)
                        - x6 * (x1408 * x5 + x1423 - x6 * (x1406 * x5 - x1424 * x6))
                    )
                )
            )
        )
    )
    result[1, 2] = numpy.sum(
        x1428
        * (
            x130
            * (
                x1470
                - x1490
                + x1494 * x3
                - x1500 * x3
                - x1500 * x6
                - x38 * (-x1442 * x3 + x1444 + x1489 + x1501)
                + x38 * (x1448 * x3 + x1450 + x1469 - x1501)
                + x6 * (x1485 + x1499)
            )
            + x3
            * (
                x130 * (x1465 - x1493 + x1494)
                + x1471 * x3
                - x1492 * x6
                + x4 * (-x1446 + x1453 + x1465)
            )
            + x38 * (x1471 - x1472 + x1487 - x1491)
            - x6
            * (
                x130 * (x1495 - x1499 + x1500)
                + x1492 * x3
                + x4 * (x1467 - x1473 + x1495)
                - x6
                * (
                    x130 * (x1462 - x1475 + x1481 + x1498)
                    + x1486 * x3
                    - x6
                    * (
                        x1482 * x5
                        + x38 * (x1460 - x1476 + x1479)
                        - x6 * (x1480 * x5 + x1496 - x6 * (x1478 * x5 - x1497 * x6))
                    )
                )
            )
        )
    )
    result[1, 3] = numpy.sum(
        x1552
        * (
            x3
            * (
                x1528 * x3
                - x1543 * x6
                + x38 * (x1524 - x1544 + x1545)
                + x4 * (-x1510 + x1516 + x1524)
            )
            + x38 * (x1528 - x1529 + x1540 - x1542)
            + x38
            * (
                x1527
                - x1541
                + x1545 * x3
                - x1550 * x3
                - x1550 * x6
                + x4 * (x1507 + x1511 * x3 + x1512 - x1551)
                - x4 * (-x1506 * x3 + x1508 + x1520 + x1551)
                + x6 * (x1538 + x1549)
            )
            - x6
            * (
                x1543 * x3
                + x38 * (x1546 - x1549 + x1550)
                + x4 * (x1526 - x1530 + x1546)
                - x6
                * (
                    x1539 * x3
                    + x38 * (x1521 - x1532 + x1535)
                    - x6 * (x1536 * x5 + x1547 - x6 * (x1534 * x5 - x1548 * x6))
                )
            )
        )
    )
    result[1, 4] = numpy.sum(
        x1603
        * (
            x3
            * (
                x1579 * x3
                - x1594 * x6
                + x38 * (x1575 - x1595 + x1596)
                + x4 * (-x1561 + x1567 + x1575)
            )
            + x38 * (x1579 - x1580 + x1591 - x1593)
            + x38
            * (
                x1578
                - x1592
                + x1596 * x3
                - x1601 * x3
                - x1601 * x6
                + x4 * (x1558 + x1562 * x3 + x1563 - x1602)
                - x4 * (-x1557 * x3 + x1559 + x1571 + x1602)
                + x6 * (x1589 + x1600)
            )
            - x6
            * (
                x1594 * x3
                + x38 * (x1597 - x1600 + x1601)
                + x4 * (x1577 - x1581 + x1597)
                - x6
                * (
                    x1590 * x3
                    + x38 * (x1572 - x1583 + x1586)
                    - x6 * (x1587 * x5 + x1598 - x6 * (x1585 * x5 - x1599 * x6))
                )
            )
        )
    )
    result[1, 5] = numpy.sum(
        x1552
        * (
            x3
            * (
                x1631 * x3
                - x1646 * x6
                + x38 * (x1627 - x1647 + x1648)
                + x4 * (-x1612 + x1618 + x1627)
            )
            + x38 * (x1631 - x1632 + x1643 - x1645)
            + x38
            * (
                x1630
                - x1644
                + x1648 * x3
                - x1654 * x3
                - x1654 * x6
                + x4 * (x1609 + x1613 * x3 + x1614 - x1655)
                - x4 * (-x1608 * x3 + x1610 + x1622 + x1655)
                + x6 * (x1641 + x1653)
            )
            - x6
            * (
                x1646 * x3
                + x38 * (x1650 - x1653 + x1654)
                + x4 * (x1629 - x1633 + x1650)
                - x6
                * (
                    x1642 * x3
                    + x38 * (x1623 - x1635 + x1638)
                    - x6 * (x1639 * x5 + x1651 - x6 * (x1637 * x5 - x1652 * x6))
                )
            )
        )
    )
    result[1, 6] = numpy.sum(
        x1428
        * (
            x3
            * (
                x1678 * x3
                - x1683 * x6
                + x4 * (x1661 - x1671 + x1672)
                + x4 * (x1672 - x1675 + x1676)
            )
            + x38 * (x1670 + x1678 - x1679 + x1682)
            + x4
            * (
                x1658
                + x1661 * x3
                + x1670
                - x1674 * x3
                - x1674 * x6
                + x6 * (x1664 + x1668)
            )
            - x6
            * (
                x1683 * x3
                + x4 * (-x1664 + x1667 + x1674)
                + x4 * (x1667 + x1677 - x1680)
                - x6 * (x1681 * x3 + x1684 - x6 * (x1666 * x5 - x1685 * x6))
            )
        )
    )
    result[1, 7] = numpy.sum(
        x1603
        * (
            x3
            * (
                x1708 * x3
                - x1713 * x6
                + x4 * (x1691 - x1701 + x1702)
                + x4 * (x1702 - x1705 + x1706)
            )
            + x38 * (x1700 + x1708 - x1709 + x1712)
            + x4
            * (
                x1688
                + x1691 * x3
                + x1700
                - x1704 * x3
                - x1704 * x6
                + x6 * (x1694 + x1698)
            )
            - x6
            * (
                x1713 * x3
                + x4 * (-x1694 + x1697 + x1704)
                + x4 * (x1697 + x1707 - x1710)
                - x6 * (x1711 * x3 + x1714 - x6 * (x1696 * x5 - x1715 * x6))
            )
        )
    )
    result[1, 8] = numpy.sum(
        x1603
        * (
            x3
            * (
                x1738 * x3
                - x1743 * x6
                + x4 * (x1721 - x1731 + x1732)
                + x4 * (x1732 - x1735 + x1736)
            )
            + x38 * (x1730 + x1738 - x1739 + x1742)
            + x4
            * (
                x1718
                + x1721 * x3
                + x1730
                - x1734 * x3
                - x1734 * x6
                + x6 * (x1724 + x1728)
            )
            - x6
            * (
                x1743 * x3
                + x4 * (-x1724 + x1727 + x1734)
                + x4 * (x1727 + x1737 - x1740)
                - x6 * (x1741 * x3 + x1744 - x6 * (x1726 * x5 - x1745 * x6))
            )
        )
    )
    result[1, 9] = numpy.sum(
        x1428
        * (
            x3
            * (
                x1768 * x3
                - x1773 * x6
                + x4 * (x1751 - x1761 + x1762)
                + x4 * (x1762 - x1765 + x1766)
            )
            + x38 * (x1760 + x1768 - x1769 + x1772)
            + x4
            * (
                x1748
                + x1751 * x3
                + x1760
                - x1764 * x3
                - x1764 * x6
                + x6 * (x1754 + x1758)
            )
            - x6
            * (
                x1773 * x3
                + x4 * (-x1754 + x1757 + x1764)
                + x4 * (x1757 + x1767 - x1770)
                - x6 * (x1771 * x3 + x1774 - x6 * (x1756 * x5 - x1775 * x6))
            )
        )
    )
    result[1, 10] = numpy.sum(
        x318
        * (
            x3 * (x1780 + x1782 * x3 - x1789 * x6)
            + x38 * (x1782 - x1783 + x1788)
            - x6 * (x1789 * x3 + x1790 - x6 * (x1787 * x3 - x1791 * x6))
        )
    )
    result[1, 11] = numpy.sum(
        x1428
        * (
            x3 * (x1796 + x1798 * x3 - x1805 * x6)
            + x38 * (x1798 - x1799 + x1804)
            - x6 * (x1805 * x3 + x1806 - x6 * (x1803 * x3 - x1807 * x6))
        )
    )
    result[1, 12] = numpy.sum(
        x1552
        * (
            x3 * (x1812 + x1814 * x3 - x1821 * x6)
            + x38 * (x1814 - x1815 + x1820)
            - x6 * (x1821 * x3 + x1822 - x6 * (x1819 * x3 - x1823 * x6))
        )
    )
    result[1, 13] = numpy.sum(
        x1428
        * (
            x3 * (x1827 + x1829 * x3 - x1835 * x6)
            + x38 * (x1829 - x1830 + x1834)
            - x6 * (x1835 * x3 + x1836 - x6 * (x1833 * x3 - x1837 * x6))
        )
    )
    result[1, 14] = numpy.sum(
        x318
        * (
            x3 * (x1841 + x1843 * x3 - x1849 * x6)
            + x38 * (x1843 - x1844 + x1848)
            - x6 * (x1849 * x3 + x1850 - x6 * (x1847 * x3 - x1851 * x6))
        )
    )
    result[2, 0] = numpy.sum(
        -x318
        * (
            x109
            * (
                x130 * (-x1878 * x3 + x1881 + x1936 + x1948)
                - x130 * (x1887 * x3 + x1888 + x1912 - x1948)
                - x1913
                + x1937
                - x1941 * x3
                + x1947 * x3
                + x1947 * x6
                - x6 * (x1933 + x1946)
            )
            - x3
            * (
                x109 * (x1909 - x1940 + x1941)
                + x1914 * x3
                - x1939 * x6
                + x4 * (-x1883 + x1892 + x1909)
            )
            - x38 * (x1914 - x1915 + x1935 - x1938)
            + x6
            * (
                x109 * (x1942 - x1946 + x1947)
                + x1939 * x3
                + x4 * (x1911 - x1916 + x1942)
                - x6
                * (
                    x109 * (x1905 - x1917 + x1927 - x1928)
                    + x1934 * x3
                    - x6
                    * (
                        x130 * (x1902 - x1919 + x1925 + x1945)
                        + x1929 * x5
                        - x6
                        * (
                            x1926 * x5
                            + x38 * (x1900 - x1920 + x1923)
                            - x6 * (x1924 * x5 + x1943 - x6 * (x1922 * x5 - x1944 * x6))
                        )
                    )
                )
            )
        )
    )
    result[2, 1] = numpy.sum(
        x1428
        * (
            x130
            * (
                x1999
                - x2021
                + x2025 * x3
                - x2032 * x3
                - x2032 * x6
                - x38 * (-x1968 * x3 + x1970 + x2020 + x2033)
                + x38 * (x1975 * x3 + x1977 + x1998 - x2033)
                + x6 * (x2016 + x2031)
            )
            + x3
            * (
                x130 * (x1994 - x2024 + x2025)
                + x2000 * x3
                - x2023 * x6
                + x4 * (-x1972 + x1980 + x1994)
            )
            + x38 * (x2000 - x2001 + x2018 - x2022)
            - x6
            * (
                x130 * (x2026 - x2031 + x2032)
                + x2023 * x3
                + x4 * (x1996 - x2002 + x2026)
                - x6
                * (
                    x130 * (x1991 - x2004 + x2012 + x2030)
                    + x2017 * x3
                    - x6
                    * (
                        x2013 * x5
                        + x38 * (x1989 - x2005 + x2010)
                        - x6 * (x2011 * x5 + x2027 - x6 * (x2009 * x5 - x2029 * x6))
                    )
                )
            )
        )
    )
    result[2, 2] = numpy.sum(
        x1428
        * (
            x130
            * (
                x2075
                - x2095
                + x2099 * x3
                - x2105 * x3
                - x2105 * x6
                - x38 * (-x2047 * x3 + x2049 + x2094 + x2106)
                + x38 * (x2053 * x3 + x2055 + x2074 - x2106)
                + x6 * (x2090 + x2104)
            )
            + x3
            * (
                x130 * (x2070 - x2098 + x2099)
                + x2076 * x3
                - x2097 * x6
                + x4 * (-x2051 + x2058 + x2070)
            )
            + x38 * (x2076 - x2077 + x2092 - x2096)
            - x6
            * (
                x130 * (x2100 - x2104 + x2105)
                + x2097 * x3
                + x4 * (x2072 - x2078 + x2100)
                - x6
                * (
                    x130 * (x2067 - x2080 + x2086 + x2103)
                    + x2091 * x3
                    - x6
                    * (
                        x2087 * x5
                        + x38 * (x2065 - x2081 + x2084)
                        - x6 * (x2085 * x5 + x2101 - x6 * (x2083 * x5 - x2102 * x6))
                    )
                )
            )
        )
    )
    result[2, 3] = numpy.sum(
        x1552
        * (
            x3
            * (
                x2143 * x3
                - x2161 * x6
                + x38 * (x2139 - x2162 + x2163)
                + x4 * (-x2119 + x2127 + x2139)
            )
            + x38 * (x2143 - x2144 + x2158 - x2160)
            + x38
            * (
                x2142
                - x2159
                + x2163 * x3
                - x2171 * x3
                - x2171 * x6
                + x4 * (x2116 + x2121 * x3 + x2122 - x2172)
                - x4 * (-x2115 * x3 + x2117 + x2133 + x2172)
                + x6 * (x2156 + x2170)
            )
            - x6
            * (
                x2161 * x3
                + x38 * (x2165 - x2170 + x2171)
                + x4 * (x2141 - x2145 + x2165)
                - x6
                * (
                    x2157 * x3
                    + x38 * (x2134 - x2148 + x2153)
                    - x6 * (x2154 * x5 + x2167 - x6 * (x2152 * x5 - x2169 * x6))
                )
            )
        )
    )
    result[2, 4] = numpy.sum(
        x1603
        * (
            x3
            * (
                x2207 * x3
                - x2224 * x6
                + x38 * (x2203 - x2225 + x2226)
                + x4 * (-x2185 + x2192 + x2203)
            )
            + x38 * (x2207 - x2208 + x2221 - x2223)
            + x38
            * (
                x2206
                - x2222
                + x2226 * x3
                - x2233 * x3
                - x2233 * x6
                + x4 * (x2182 + x2187 * x3 + x2188 - x2234)
                - x4 * (-x2181 * x3 + x2183 + x2198 + x2234)
                + x6 * (x2219 + x2232)
            )
            - x6
            * (
                x2224 * x3
                + x38 * (x2228 - x2232 + x2233)
                + x4 * (x2205 - x2209 + x2228)
                - x6
                * (
                    x2220 * x3
                    + x38 * (x2199 - x2211 + x2216)
                    - x6 * (x2217 * x5 + x2229 - x6 * (x2215 * x5 - x2231 * x6))
                )
            )
        )
    )
    result[2, 5] = numpy.sum(
        x1552
        * (
            x3
            * (
                x2262 * x3
                - x2277 * x6
                + x38 * (x2258 - x2278 + x2279)
                + x4 * (-x2243 + x2249 + x2258)
            )
            + x38 * (x2262 - x2263 + x2274 - x2276)
            + x38
            * (
                x2261
                - x2275
                + x2279 * x3
                - x2285 * x3
                - x2285 * x6
                + x4 * (x2240 + x2244 * x3 + x2245 - x2286)
                - x4 * (-x2239 * x3 + x2241 + x2253 + x2286)
                + x6 * (x2272 + x2284)
            )
            - x6
            * (
                x2277 * x3
                + x38 * (x2281 - x2284 + x2285)
                + x4 * (x2260 - x2264 + x2281)
                - x6
                * (
                    x2273 * x3
                    + x38 * (x2254 - x2266 + x2269)
                    - x6 * (x2270 * x5 + x2282 - x6 * (x2268 * x5 - x2283 * x6))
                )
            )
        )
    )
    result[2, 6] = numpy.sum(
        x1428
        * (
            x3
            * (
                x2321 * x3
                - x2326 * x6
                + x4 * (x2297 - x2314 + x2315)
                + x4 * (x2315 - x2318 + x2319)
            )
            + x38 * (x2313 + x2321 - x2322 + x2325)
            + x4
            * (
                x2293
                + x2297 * x3
                + x2313
                - x2317 * x3
                - x2317 * x6
                + x6 * (x2303 + x2310)
            )
            - x6
            * (
                x2326 * x3
                + x4 * (-x2303 + x2309 + x2317)
                + x4 * (x2309 + x2320 - x2323)
                - x6 * (x2324 * x3 + x2328 - x6 * (x2308 * x5 - x2330 * x6))
            )
        )
    )
    result[2, 7] = numpy.sum(
        x1603
        * (
            x3
            * (
                x2362 * x3
                - x2367 * x6
                + x4 * (x2340 - x2355 + x2356)
                + x4 * (x2356 - x2359 + x2360)
            )
            + x38 * (x2354 + x2362 - x2363 + x2366)
            + x4
            * (
                x2336
                + x2340 * x3
                + x2354
                - x2358 * x3
                - x2358 * x6
                + x6 * (x2345 + x2351)
            )
            - x6
            * (
                x2367 * x3
                + x4 * (-x2345 + x2350 + x2358)
                + x4 * (x2350 + x2361 - x2364)
                - x6 * (x2365 * x3 + x2369 - x6 * (x2349 * x5 - x2371 * x6))
            )
        )
    )
    result[2, 8] = numpy.sum(
        x1603
        * (
            x3
            * (
                x2401 * x3
                - x2406 * x6
                + x4 * (x2380 - x2394 + x2395)
                + x4 * (x2395 - x2398 + x2399)
            )
            + x38 * (x2393 + x2401 - x2402 + x2405)
            + x4
            * (
                x2376
                + x2380 * x3
                + x2393
                - x2397 * x3
                - x2397 * x6
                + x6 * (x2385 + x2391)
            )
            - x6
            * (
                x2406 * x3
                + x4 * (-x2385 + x2390 + x2397)
                + x4 * (x2390 + x2400 - x2403)
                - x6 * (x2404 * x3 + x2407 - x6 * (x2389 * x5 - x2409 * x6))
            )
        )
    )
    result[2, 9] = numpy.sum(
        x1428
        * (
            x3
            * (
                x2432 * x3
                - x2437 * x6
                + x4 * (x2415 - x2425 + x2426)
                + x4 * (x2426 - x2429 + x2430)
            )
            + x38 * (x2424 + x2432 - x2433 + x2436)
            + x4
            * (
                x2412
                + x2415 * x3
                + x2424
                - x2428 * x3
                - x2428 * x6
                + x6 * (x2418 + x2422)
            )
            - x6
            * (
                x2437 * x3
                + x4 * (-x2418 + x2421 + x2428)
                + x4 * (x2421 + x2431 - x2434)
                - x6 * (x2435 * x3 + x2438 - x6 * (x2420 * x5 - x2439 * x6))
            )
        )
    )
    result[2, 10] = numpy.sum(
        x318
        * (
            x3 * (x2446 + x2449 * x3 - x2457 * x6)
            + x38 * (x2449 - x2450 + x2456)
            - x6 * (x2457 * x3 + x2459 - x6 * (x2455 * x3 - x2461 * x6))
        )
    )
    result[2, 11] = numpy.sum(
        x1428
        * (
            x3 * (x2468 + x2471 * x3 - x2479 * x6)
            + x38 * (x2471 - x2472 + x2478)
            - x6 * (x2479 * x3 + x2481 - x6 * (x2477 * x3 - x2483 * x6))
        )
    )
    result[2, 12] = numpy.sum(
        x1552
        * (
            x3 * (x2489 + x2492 * x3 - x2499 * x6)
            + x38 * (x2492 - x2493 + x2498)
            - x6 * (x2499 * x3 + x2501 - x6 * (x2497 * x3 - x2503 * x6))
        )
    )
    result[2, 13] = numpy.sum(
        x1428
        * (
            x3 * (x2508 + x2511 * x3 - x2518 * x6)
            + x38 * (x2511 - x2512 + x2517)
            - x6 * (x2518 * x3 + x2519 - x6 * (x2516 * x3 - x2521 * x6))
        )
    )
    result[2, 14] = numpy.sum(
        x318
        * (
            x3 * (x2526 + x2528 * x3 - x2535 * x6)
            + x38 * (x2528 - x2529 + x2534)
            - x6 * (x2535 * x3 + x2536 - x6 * (x2533 * x3 - x2537 * x6))
        )
    )
    result[3, 0] = numpy.sum(
        x578
        * (
            x109 * (-x2570 * x3 + x2576 * x3 + x2578 + x2597)
            - x3 * (x109 * (x2557 - x2568 + x2569 - x2576) - x2579 * x3 + x2598 * x6)
            + x4 * (-x2571 + x2579 + x2597)
            - x6
            * (
                x109 * (x2570 - x2580 + x2591 - x2592)
                + x2598 * x3
                - x6
                * (
                    x130 * (x2567 - x2582 + x2589 + x2601)
                    + x2593 * x5
                    - x6
                    * (
                        x2590 * x5
                        + x38 * (x2565 - x2583 + x2587)
                        - x6 * (x2588 * x5 + x2599 - x6 * (x2586 * x5 - x2600 * x6))
                    )
                )
            )
        )
    )
    result[3, 1] = numpy.sum(
        x1552
        * (
            x130 * (-x2620 * x3 + x2625 * x3 + x2626 + x2640)
            - x3 * (x130 * (x2607 + x2612 - x2619 - x2625) - x2627 * x3 + x2641 * x6)
            + x4 * (-x2621 + x2627 + x2640)
            - x6
            * (
                x130 * (x2620 - x2628 - x2629 + x2636)
                + x2641 * x3
                - x6
                * (
                    x2637 * x5
                    + x38 * (x2618 - x2630 + x2634)
                    - x6 * (x2635 * x5 + x2642 - x6 * (x2633 * x5 - x2643 * x6))
                )
            )
        )
    )
    result[3, 2] = numpy.sum(
        x1552
        * (
            x130 * (-x2662 * x3 + x2667 * x3 + x2668 + x2682)
            + x3 * (x130 * (-x2654 + x2661 + x2667 + x2684) + x2669 * x3 - x2683 * x6)
            + x4 * (-x2663 + x2669 + x2682)
            - x6
            * (
                x130 * (x2662 - x2671 + x2678 + x2687)
                + x2683 * x3
                - x6
                * (
                    x2679 * x5
                    + x38 * (x2660 - x2672 + x2676)
                    - x6 * (x2677 * x5 + x2685 - x6 * (x2675 * x5 - x2686 * x6))
                )
            )
        )
    )
    result[3, 3] = numpy.sum(
        x2723
        * (
            x3 * (x2706 * x3 - x2718 * x6 + x38 * (-x2692 + x2697 + x2702))
            + x38 * (-x2698 * x3 + x2702 * x3 + x2705 + x2717)
            + x4 * (-x2699 + x2706 + x2717)
            - x6
            * (
                x2718 * x3
                + x38 * (x2698 - x2707 + x2712)
                - x6 * (x2713 * x5 + x2720 - x6 * (x2711 * x5 - x2722 * x6))
            )
        )
    )
    result[3, 4] = numpy.sum(
        x2751
        * (
            x3 * (x2738 * x3 - x2748 * x6 + x38 * (-x2727 + x2731 + x2735))
            + x38 * (-x2732 * x3 + x2735 * x3 + x2737 + x2747)
            + x4 * (-x2733 + x2738 + x2747)
            - x6
            * (
                x2748 * x3
                + x38 * (x2732 - x2739 + x2743)
                - x6 * (x2744 * x5 + x2749 - x6 * (x2742 * x5 - x2750 * x6))
            )
        )
    )
    result[3, 5] = numpy.sum(
        x2723
        * (
            x3 * (x2766 * x3 - x2777 * x6 + x38 * (-x2755 + x2759 + x2763))
            + x38 * (-x2760 * x3 + x2763 * x3 + x2765 + x2776)
            + x4 * (-x2761 + x2766 + x2776)
            - x6
            * (
                x2777 * x3
                + x38 * (x2760 - x2767 + x2771)
                - x6 * (x2772 * x5 + x2778 - x6 * (x2770 * x5 - x2779 * x6))
            )
        )
    )
    result[3, 6] = numpy.sum(
        x1552
        * (
            x3 * (x2792 * x3 + x2793 - x2794 * x6)
            + x4 * (x2785 - x2791 + x2792)
            + x4 * (x2780 * x3 + x2785 - x2789 * x3 + x2790)
            - x6 * (x2794 * x3 + x2795 - x6 * (x2784 * x5 - x2796 * x6))
        )
    )
    result[3, 7] = numpy.sum(
        x2751
        * (
            x3 * (x2810 * x3 + x2812 - x2813 * x6)
            + x4 * (x2803 - x2809 + x2810)
            + x4 * (x2798 * x3 + x2803 - x2807 * x3 + x2808)
            - x6 * (x2813 * x3 + x2815 - x6 * (x2802 * x5 - x2817 * x6))
        )
    )
    result[3, 8] = numpy.sum(
        x2751
        * (
            x3 * (x2828 * x3 + x2829 - x2830 * x6)
            + x4 * (x2822 - x2827 + x2828)
            + x4 * (x2818 * x3 + x2822 - x2825 * x3 + x2826)
            - x6 * (x2830 * x3 + x2831 - x6 * (x2821 * x5 - x2832 * x6))
        )
    )
    result[3, 9] = numpy.sum(
        x1552
        * (
            x3 * (x2843 * x3 + x2844 - x2845 * x6)
            + x4 * (x2837 - x2842 + x2843)
            + x4 * (x2833 * x3 + x2837 - x2840 * x3 + x2841)
            - x6 * (x2845 * x3 + x2846 - x6 * (x2836 * x5 - x2847 * x6))
        )
    )
    result[3, 10] = numpy.sum(
        x578 * (x2852 + x3 * (x2851 * x3 - x2853 * x6) - x6 * (x2853 * x3 - x2854 * x6))
    )
    result[3, 11] = numpy.sum(
        x1552 * (x2859 + x3 * (x2858 * x3 - x2860 * x6) - x6 * (x2860 * x3 - x2861 * x6))
    )
    result[3, 12] = numpy.sum(
        x2723 * (x2868 + x3 * (x2867 * x3 - x2869 * x6) - x6 * (x2869 * x3 - x2871 * x6))
    )
    result[3, 13] = numpy.sum(
        x1552 * (x2875 + x3 * (x2874 * x3 - x2876 * x6) - x6 * (x2876 * x3 - x2877 * x6))
    )
    result[3, 14] = numpy.sum(
        x578 * (x2881 + x3 * (x2880 * x3 - x2882 * x6) - x6 * (x2882 * x3 - x2883 * x6))
    )
    result[4, 0] = numpy.sum(
        x687
        * (
            x109 * (-x2912 * x3 + x2918 * x3 + x2920 + x2938)
            - x3 * (x109 * (x2900 - x2910 + x2911 - x2918) - x2921 * x3 + x2939 * x6)
            + x4 * (-x2913 + x2921 + x2938)
            - x6
            * (
                x109 * (x2912 - x2922 + x2932 - x2933)
                + x2939 * x3
                - x6
                * (
                    x130 * (x2909 - x2924 + x2930 + x2942)
                    + x2934 * x5
                    - x6
                    * (
                        x2931 * x5
                        + x38 * (x2907 - x2925 + x2928)
                        - x6 * (x2929 * x5 + x2940 - x6 * (x2927 * x5 - x2941 * x6))
                    )
                )
            )
        )
    )
    result[4, 1] = numpy.sum(
        x1603
        * (
            x130 * (-x2958 * x3 + x2963 * x3 + x2964 + x2977)
            - x3 * (x130 * (x2946 + x2951 - x2957 - x2963) - x2965 * x3 + x2978 * x6)
            + x4 * (-x2959 + x2965 + x2977)
            - x6
            * (
                x130 * (x2958 - x2966 - x2967 + x2973)
                + x2978 * x3
                - x6
                * (
                    x2974 * x5
                    + x38 * (x2956 - x2968 + x2971)
                    - x6 * (x2972 * x5 + x2979 - x6 * (x2970 * x5 - x2980 * x6))
                )
            )
        )
    )
    result[4, 2] = numpy.sum(
        x1603
        * (
            x130 * (-x2996 * x3 + x3 * x3001 + x3002 + x3015)
            + x3 * (x130 * (-x2989 + x2995 + x3001 + x3017) + x3 * x3003 - x3016 * x6)
            + x4 * (-x2997 + x3003 + x3015)
            - x6
            * (
                x130 * (x2996 - x3005 + x3011 + x3020)
                + x3 * x3016
                - x6
                * (
                    x3012 * x5
                    + x38 * (x2994 - x3006 + x3009)
                    - x6 * (x3010 * x5 + x3018 - x6 * (x3008 * x5 - x3019 * x6))
                )
            )
        )
    )
    result[4, 3] = numpy.sum(
        x2751
        * (
            x3 * (x3 * x3033 - x3042 * x6 + x38 * (-x3023 + x3026 + x3030))
            + x38 * (-x3 * x3027 + x3 * x3030 + x3032 + x3041)
            + x4 * (-x3028 + x3033 + x3041)
            - x6
            * (
                x3 * x3042
                + x38 * (x3027 - x3034 + x3037)
                - x6 * (x3038 * x5 + x3043 - x6 * (x3036 * x5 - x3044 * x6))
            )
        )
    )
    result[4, 4] = numpy.sum(
        x3069
        * (
            x3 * (x3 * x3057 - x3066 * x6 + x38 * (-x3047 + x3050 + x3054))
            + x38 * (-x3 * x3051 + x3 * x3054 + x3056 + x3065)
            + x4 * (-x3052 + x3057 + x3065)
            - x6
            * (
                x3 * x3066
                + x38 * (x3051 - x3058 + x3061)
                - x6 * (x3062 * x5 + x3067 - x6 * (x3060 * x5 - x3068 * x6))
            )
        )
    )
    result[4, 5] = numpy.sum(
        x2751
        * (
            x3 * (x3 * x3082 - x3092 * x6 + x38 * (-x3072 + x3075 + x3079))
            + x38 * (-x3 * x3076 + x3 * x3079 + x3081 + x3091)
            + x4 * (-x3077 + x3082 + x3091)
            - x6
            * (
                x3 * x3092
                + x38 * (x3076 - x3083 + x3086)
                - x6 * (x3087 * x5 + x3093 - x6 * (x3085 * x5 - x3094 * x6))
            )
        )
    )
    result[4, 6] = numpy.sum(
        x1603
        * (
            x3 * (x3 * x3103 + x3104 - x3105 * x6)
            + x4 * (x3098 - x3102 + x3103)
            + x4 * (x3 * x3095 - x3 * x3100 + x3098 + x3101)
            - x6 * (x3 * x3105 + x3106 - x6 * (x3097 * x5 - x3107 * x6))
        )
    )
    result[4, 7] = numpy.sum(
        x3069
        * (
            x3 * (x3 * x3116 + x3117 - x3118 * x6)
            + x4 * (x3111 - x3115 + x3116)
            + x4 * (x3 * x3108 - x3 * x3113 + x3111 + x3114)
            - x6 * (x3 * x3118 + x3119 - x6 * (x3110 * x5 - x3120 * x6))
        )
    )
    result[4, 8] = numpy.sum(
        x3069
        * (
            x3 * (x3 * x3129 + x3130 - x3131 * x6)
            + x4 * (x3124 - x3128 + x3129)
            + x4 * (x3 * x3121 - x3 * x3126 + x3124 + x3127)
            - x6 * (x3 * x3131 + x3132 - x6 * (x3123 * x5 - x3133 * x6))
        )
    )
    result[4, 9] = numpy.sum(
        x1603
        * (
            x3 * (x3 * x3142 + x3143 - x3144 * x6)
            + x4 * (x3137 - x3141 + x3142)
            + x4 * (x3 * x3134 - x3 * x3139 + x3137 + x3140)
            - x6 * (x3 * x3144 + x3145 - x6 * (x3136 * x5 - x3146 * x6))
        )
    )
    result[4, 10] = numpy.sum(
        x687 * (x3 * (x3 * x3150 - x3152 * x6) + x3151 - x6 * (x3 * x3152 - x3153 * x6))
    )
    result[4, 11] = numpy.sum(
        x1603 * (x3 * (x3 * x3157 - x3159 * x6) + x3158 - x6 * (x3 * x3159 - x3160 * x6))
    )
    result[4, 12] = numpy.sum(
        x2751 * (x3 * (x3 * x3164 - x3166 * x6) + x3165 - x6 * (x3 * x3166 - x3167 * x6))
    )
    result[4, 13] = numpy.sum(
        x1603 * (x3 * (x3 * x3170 - x3172 * x6) + x3171 - x6 * (x3 * x3172 - x3173 * x6))
    )
    result[4, 14] = numpy.sum(
        x687 * (x3 * (x3 * x3176 - x3178 * x6) + x3177 - x6 * (x3 * x3178 - x3179 * x6))
    )
    result[5, 0] = numpy.sum(
        x578
        * (
            x109 * (-x3 * x3212 + x3 * x3218 + x3220 + x3239)
            - x3 * (x109 * (x3199 - x3210 + x3211 - x3218) - x3 * x3221 + x3240 * x6)
            + x4 * (-x3213 + x3221 + x3239)
            - x6
            * (
                x109 * (x3212 - x3222 + x3233 - x3234)
                + x3 * x3240
                - x6
                * (
                    x130 * (x3209 - x3224 + x3231 + x3243)
                    + x3235 * x5
                    - x6
                    * (
                        x3232 * x5
                        + x38 * (x3207 - x3225 + x3229)
                        - x6 * (x3230 * x5 + x3241 - x6 * (x3228 * x5 - x3242 * x6))
                    )
                )
            )
        )
    )
    result[5, 1] = numpy.sum(
        x1552
        * (
            x130 * (-x3 * x3265 + x3 * x3271 + x3272 + x3287)
            + x3 * (x130 * (-x3256 + x3264 + x3271 + x3289) + x3 * x3273 - x3288 * x6)
            + x4 * (-x3266 + x3273 + x3287)
            - x6
            * (
                x130 * (x3265 - x3275 + x3283 + x3293)
                + x3 * x3288
                - x6
                * (
                    x3284 * x5
                    + x38 * (x3263 - x3276 + x3281)
                    - x6 * (x3282 * x5 + x3290 - x6 * (x3280 * x5 - x3292 * x6))
                )
            )
        )
    )
    result[5, 2] = numpy.sum(
        x1552
        * (
            x130 * (-x3 * x3312 + x3 * x3317 + x3318 + x3332)
            + x3 * (x130 * (-x3304 + x3311 + x3317 + x3334) + x3 * x3319 - x3333 * x6)
            + x4 * (-x3313 + x3319 + x3332)
            - x6
            * (
                x130 * (x3312 - x3321 + x3328 + x3337)
                + x3 * x3333
                - x6
                * (
                    x3329 * x5
                    + x38 * (x3310 - x3322 + x3326)
                    - x6 * (x3327 * x5 + x3335 - x6 * (x3325 * x5 - x3336 * x6))
                )
            )
        )
    )
    result[5, 3] = numpy.sum(
        x2723
        * (
            x3 * (x3 * x3356 - x3369 * x6 + x38 * (-x3342 + x3347 + x3352))
            + x38 * (-x3 * x3348 + x3 * x3352 + x3355 + x3368)
            + x4 * (-x3349 + x3356 + x3368)
            - x6
            * (
                x3 * x3369
                + x38 * (x3348 - x3357 + x3362)
                - x6 * (x3363 * x5 + x3371 - x6 * (x3361 * x5 - x3373 * x6))
            )
        )
    )
    result[5, 4] = numpy.sum(
        x2751
        * (
            x3 * (x3 * x3391 - x3403 * x6 + x38 * (-x3378 + x3383 + x3388))
            + x38 * (-x3 * x3384 + x3 * x3388 + x3390 + x3402)
            + x4 * (-x3385 + x3391 + x3402)
            - x6
            * (
                x3 * x3403
                + x38 * (x3384 - x3392 + x3397)
                - x6 * (x3398 * x5 + x3404 - x6 * (x3396 * x5 - x3406 * x6))
            )
        )
    )
    result[5, 5] = numpy.sum(
        x2723
        * (
            x3 * (x3 * x3425 - x3438 * x6 + x38 * (-x3411 + x3416 + x3421))
            + x38 * (-x3 * x3417 + x3 * x3421 + x3424 + x3437)
            + x4 * (-x3418 + x3425 + x3437)
            - x6
            * (
                x3 * x3438
                + x38 * (x3417 - x3426 + x3431)
                - x6 * (x3432 * x5 + x3440 - x6 * (x3430 * x5 - x3442 * x6))
            )
        )
    )
    result[5, 6] = numpy.sum(
        x1552
        * (
            x3 * (x3 * x3458 + x3460 - x3461 * x6)
            + x4 * (x3450 - x3457 + x3458)
            + x4 * (x3 * x3444 - x3 * x3455 + x3450 + x3456)
            - x6 * (x3 * x3461 + x3463 - x6 * (x3449 * x5 - x3465 * x6))
        )
    )
    result[5, 7] = numpy.sum(
        x2751
        * (
            x3 * (x3 * x3479 + x3481 - x3482 * x6)
            + x4 * (x3472 - x3478 + x3479)
            + x4 * (x3 * x3467 - x3 * x3476 + x3472 + x3477)
            - x6 * (x3 * x3482 + x3484 - x6 * (x3471 * x5 - x3486 * x6))
        )
    )
    result[5, 8] = numpy.sum(
        x2751
        * (
            x3 * (x3 * x3500 + x3501 - x3502 * x6)
            + x4 * (x3493 - x3499 + x3500)
            + x4 * (x3 * x3488 - x3 * x3497 + x3493 + x3498)
            - x6 * (x3 * x3502 + x3503 - x6 * (x3492 * x5 - x3505 * x6))
        )
    )
    result[5, 9] = numpy.sum(
        x1552
        * (
            x3 * (x3 * x3518 + x3519 - x3520 * x6)
            + x4 * (x3511 - x3517 + x3518)
            + x4 * (x3 * x3506 - x3 * x3515 + x3511 + x3516)
            - x6 * (x3 * x3520 + x3521 - x6 * (x3510 * x5 - x3522 * x6))
        )
    )
    result[5, 10] = numpy.sum(
        x578 * (x3 * (x3 * x3525 - x3531 * x6) + x3529 - x6 * (x3 * x3531 - x3533 * x6))
    )
    result[5, 11] = numpy.sum(
        x1552 * (x3 * (x3 * x3536 - x3542 * x6) + x3540 - x6 * (x3 * x3542 - x3544 * x6))
    )
    result[5, 12] = numpy.sum(
        x2723 * (x3 * (x3 * x3547 - x3552 * x6) + x3550 - x6 * (x3 * x3552 - x3554 * x6))
    )
    result[5, 13] = numpy.sum(
        x1552 * (x3 * (x3 * x3558 - x3561 * x6) + x3559 - x6 * (x3 * x3561 - x3563 * x6))
    )
    result[5, 14] = numpy.sum(
        x578 * (x3 * (x3 * x3567 - x3569 * x6) + x3568 - x6 * (x3 * x3569 - x3570 * x6))
    )
    result[6, 0] = numpy.sum(
        x318
        * (
            x109 * (x3597 - x3598 + x3610 - x3611)
            + x3 * (x130 * (x3585 - x3587 + x3595 + x3613) + x3597 * x5 - x3612 * x6)
            - x6
            * (
                x130 * (x3596 - x3600 + x3608 + x3616)
                + x3612 * x5
                - x6
                * (
                    x3609 * x5
                    + x38 * (x3594 - x3601 + x3606)
                    - x6 * (x3607 * x5 + x3614 - x6 * (x3605 * x5 - x3615 * x6))
                )
            )
        )
    )
    result[6, 1] = numpy.sum(
        x1428
        * (
            x130 * (x3629 - x3630 - x3631 + x3638)
            + x3 * (x3629 * x5 - x3639 * x6 + x38 * (x3622 - x3623 + x3627))
            - x6
            * (
                x3639 * x5
                + x38 * (x3628 - x3632 + x3636)
                - x6 * (x3637 * x5 + x3640 - x6 * (x3635 * x5 - x3641 * x6))
            )
        )
    )
    result[6, 2] = numpy.sum(
        x1428
        * (
            x130 * (x3656 - x3658 + x3666 + x3670)
            + x3 * (x3656 * x5 - x3667 * x6 + x38 * (x3648 - x3649 + x3654))
            - x6
            * (
                x3667 * x5
                + x38 * (x3655 - x3659 + x3664)
                - x6 * (x3665 * x5 + x3668 - x6 * (x3663 * x5 - x3669 * x6))
            )
        )
    )
    result[6, 3] = numpy.sum(
        x1552
        * (
            x3 * (x3676 + x3678 * x5 - x3686 * x6)
            + x38 * (x3678 - x3679 + x3685)
            - x6 * (x3686 * x5 + x3687 - x6 * (x3684 * x5 - x3688 * x6))
        )
    )
    result[6, 4] = numpy.sum(
        x1603
        * (
            x3 * (x3692 + x3694 * x5 - x3700 * x6)
            + x38 * (x3694 - x3695 + x3699)
            - x6 * (x3700 * x5 + x3701 - x6 * (x3698 * x5 - x3702 * x6))
        )
    )
    result[6, 5] = numpy.sum(
        x1552
        * (
            x3 * (x3707 + x3709 * x5 - x3716 * x6)
            + x38 * (x3709 - x3710 + x3715)
            - x6 * (x3716 * x5 + x3717 - x6 * (x3714 * x5 - x3718 * x6))
        )
    )
    result[6, 6] = numpy.sum(
        x1428 * (x3 * (x3725 * x5 - x3727 * x6) + x3726 - x6 * (x3727 * x5 - x3729 * x6))
    )
    result[6, 7] = numpy.sum(
        x1603 * (x3 * (x3734 * x5 - x3736 * x6) + x3735 - x6 * (x3736 * x5 - x3737 * x6))
    )
    result[6, 8] = numpy.sum(
        x1603 * (x3 * (x3740 * x5 - x3742 * x6) + x3741 - x6 * (x3742 * x5 - x3743 * x6))
    )
    result[6, 9] = numpy.sum(
        x1428 * (x3 * (x3747 * x5 - x3749 * x6) + x3748 - x6 * (x3749 * x5 - x3750 * x6))
    )
    result[6, 10] = numpy.sum(x318 * (x3 * x3751 - x3752 * x6))
    result[6, 11] = numpy.sum(x1428 * (x3 * x3754 - x3756 * x6))
    result[6, 12] = numpy.sum(x1552 * (x3 * x3757 - x3758 * x6))
    result[6, 13] = numpy.sum(x1428 * (x3 * x3759 - x3760 * x6))
    result[6, 14] = numpy.sum(x318 * (x3 * x3761 - x3762 * x6))
    result[7, 0] = numpy.sum(
        x687
        * (
            x109 * (x3786 - x3787 + x3798 - x3799)
            + x3 * (x130 * (x3775 - x3777 + x3784 + x3801) + x3786 * x5 - x3800 * x6)
            - x6
            * (
                x130 * (x3785 - x3789 + x3796 + x3804)
                + x3800 * x5
                - x6
                * (
                    x3797 * x5
                    + x38 * (x3783 - x3790 + x3794)
                    - x6 * (x3795 * x5 + x3802 - x6 * (x3793 * x5 - x3803 * x6))
                )
            )
        )
    )
    result[7, 1] = numpy.sum(
        x1603
        * (
            x130 * (x3817 - x3818 - x3819 + x3826)
            + x3 * (x38 * (x3810 - x3811 + x3815) + x3817 * x5 - x3827 * x6)
            - x6
            * (
                x38 * (x3816 - x3820 + x3824)
                + x3827 * x5
                - x6 * (x3825 * x5 + x3828 - x6 * (x3823 * x5 - x3829 * x6))
            )
        )
    )
    result[7, 2] = numpy.sum(
        x1603
        * (
            x130 * (x3842 - x3844 + x3851 + x3855)
            + x3 * (x38 * (x3835 - x3836 + x3840) + x3842 * x5 - x3852 * x6)
            - x6
            * (
                x38 * (x3841 - x3845 + x3849)
                + x3852 * x5
                - x6 * (x3850 * x5 + x3853 - x6 * (x3848 * x5 - x3854 * x6))
            )
        )
    )
    result[7, 3] = numpy.sum(
        x2751
        * (
            x3 * (x3862 + x3864 * x5 - x3871 * x6)
            + x38 * (x3864 - x3865 + x3870)
            - x6 * (x3871 * x5 + x3873 - x6 * (x3869 * x5 - x3875 * x6))
        )
    )
    result[7, 4] = numpy.sum(
        x3069
        * (
            x3 * (x3879 + x3881 * x5 - x3887 * x6)
            + x38 * (x3881 - x3882 + x3886)
            - x6 * (x3887 * x5 + x3888 - x6 * (x3885 * x5 - x3889 * x6))
        )
    )
    result[7, 5] = numpy.sum(
        x2751
        * (
            x3 * (x3893 + x3895 * x5 - x3901 * x6)
            + x38 * (x3895 - x3896 + x3900)
            - x6 * (x3901 * x5 + x3902 - x6 * (x3899 * x5 - x3903 * x6))
        )
    )
    result[7, 6] = numpy.sum(
        x1603 * (x3 * (x3907 * x5 - x3909 * x6) + x3908 - x6 * (x3909 * x5 - x3910 * x6))
    )
    result[7, 7] = numpy.sum(
        x3069 * (x3 * (x3916 * x5 - x3918 * x6) + x3917 - x6 * (x3918 * x5 - x3920 * x6))
    )
    result[7, 8] = numpy.sum(
        x3069 * (x3 * (x3923 * x5 - x3925 * x6) + x3924 - x6 * (x3925 * x5 - x3926 * x6))
    )
    result[7, 9] = numpy.sum(
        x1603 * (x3 * (x3929 * x5 - x3931 * x6) + x3930 - x6 * (x3931 * x5 - x3932 * x6))
    )
    result[7, 10] = numpy.sum(x687 * (x3 * x3933 - x3934 * x6))
    result[7, 11] = numpy.sum(x1603 * (x3 * x3935 - x3936 * x6))
    result[7, 12] = numpy.sum(x2751 * (x3 * x3938 - x3940 * x6))
    result[7, 13] = numpy.sum(x1603 * (x3 * x3941 - x3942 * x6))
    result[7, 14] = numpy.sum(x687 * (x3 * x3943 - x3944 * x6))
    result[8, 0] = numpy.sum(
        x687
        * (
            x109 * (x3965 - x3966 + x3976 - x3977)
            + x3 * (x130 * (x3955 - x3957 + x3963 + x3979) + x3965 * x5 - x3978 * x6)
            - x6
            * (
                x130 * (x3964 - x3968 + x3974 + x3982)
                + x3978 * x5
                - x6
                * (
                    x38 * (x3962 - x3969 + x3972)
                    + x3975 * x5
                    - x6 * (x3973 * x5 + x3980 - x6 * (x3971 * x5 - x3981 * x6))
                )
            )
        )
    )
    result[8, 1] = numpy.sum(
        x1603
        * (
            x130 * (x3993 - x3994 - x3995 + x4001)
            + x3 * (x38 * (x3987 - x3988 + x3991) + x3993 * x5 - x4002 * x6)
            - x6
            * (
                x38 * (x3992 - x3996 + x3999)
                + x4002 * x5
                - x6 * (x4000 * x5 + x4003 - x6 * (x3998 * x5 - x4004 * x6))
            )
        )
    )
    result[8, 2] = numpy.sum(
        x1603
        * (
            x130 * (x4015 - x4017 + x4023 + x4027)
            + x3 * (x38 * (x4009 - x4010 + x4013) + x4015 * x5 - x4024 * x6)
            - x6
            * (
                x38 * (x4014 - x4018 + x4021)
                + x4024 * x5
                - x6 * (x4022 * x5 + x4025 - x6 * (x4020 * x5 - x4026 * x6))
            )
        )
    )
    result[8, 3] = numpy.sum(
        x2751
        * (
            x3 * (x4030 + x4032 * x5 - x4037 * x6)
            + x38 * (x4032 - x4033 + x4036)
            - x6 * (x4037 * x5 + x4038 - x6 * (x4035 * x5 - x4039 * x6))
        )
    )
    result[8, 4] = numpy.sum(
        x3069
        * (
            x3 * (x4042 + x4044 * x5 - x4049 * x6)
            + x38 * (x4044 - x4045 + x4048)
            - x6 * (x4049 * x5 + x4050 - x6 * (x4047 * x5 - x4051 * x6))
        )
    )
    result[8, 5] = numpy.sum(
        x2751
        * (
            x3 * (x4054 + x4056 * x5 - x4061 * x6)
            + x38 * (x4056 - x4057 + x4060)
            - x6 * (x4061 * x5 + x4062 - x6 * (x4059 * x5 - x4063 * x6))
        )
    )
    result[8, 6] = numpy.sum(
        x1603 * (x3 * (x4065 * x5 - x4067 * x6) + x4066 - x6 * (x4067 * x5 - x4068 * x6))
    )
    result[8, 7] = numpy.sum(
        x3069 * (x3 * (x4070 * x5 - x4072 * x6) + x4071 - x6 * (x4072 * x5 - x4073 * x6))
    )
    result[8, 8] = numpy.sum(
        x3069 * (x3 * (x4075 * x5 - x4077 * x6) + x4076 - x6 * (x4077 * x5 - x4078 * x6))
    )
    result[8, 9] = numpy.sum(
        x1603 * (x3 * (x4080 * x5 - x4082 * x6) + x4081 - x6 * (x4082 * x5 - x4083 * x6))
    )
    result[8, 10] = numpy.sum(x687 * (x3 * x4084 - x4085 * x6))
    result[8, 11] = numpy.sum(x1603 * (x3 * x4086 - x4087 * x6))
    result[8, 12] = numpy.sum(x2751 * (x3 * x4088 - x4089 * x6))
    result[8, 13] = numpy.sum(x1603 * (x3 * x4090 - x4091 * x6))
    result[8, 14] = numpy.sum(x687 * (x3 * x4092 - x4093 * x6))
    result[9, 0] = numpy.sum(
        x318
        * (
            x109 * (x4120 - x4121 + x4133 - x4134)
            + x3 * (x130 * (x4108 - x4110 + x4118 + x4136) + x4120 * x5 - x4135 * x6)
            - x6
            * (
                x130 * (x4119 - x4123 + x4131 + x4139)
                + x4135 * x5
                - x6
                * (
                    x38 * (x4117 - x4124 + x4129)
                    + x4132 * x5
                    - x6 * (x4130 * x5 + x4137 - x6 * (x4128 * x5 - x4138 * x6))
                )
            )
        )
    )
    result[9, 1] = numpy.sum(
        x1428
        * (
            x130 * (x4155 - x4156 - x4157 + x4165)
            + x3 * (x38 * (x4147 - x4148 + x4153) + x4155 * x5 - x4166 * x6)
            - x6
            * (
                x38 * (x4154 - x4158 + x4163)
                + x4166 * x5
                - x6 * (x4164 * x5 + x4167 - x6 * (x4162 * x5 - x4169 * x6))
            )
        )
    )
    result[9, 2] = numpy.sum(
        x1428
        * (
            x130 * (x4182 - x4184 + x4191 + x4195)
            + x3 * (x38 * (x4175 - x4176 + x4180) + x4182 * x5 - x4192 * x6)
            - x6
            * (
                x38 * (x4181 - x4185 + x4189)
                + x4192 * x5
                - x6 * (x4190 * x5 + x4193 - x6 * (x4188 * x5 - x4194 * x6))
            )
        )
    )
    result[9, 3] = numpy.sum(
        x1552
        * (
            x3 * (x4201 + x4204 * x5 - x4211 * x6)
            + x38 * (x4204 - x4205 + x4210)
            - x6 * (x4211 * x5 + x4213 - x6 * (x4209 * x5 - x4215 * x6))
        )
    )
    result[9, 4] = numpy.sum(
        x1603
        * (
            x3 * (x4220 + x4223 * x5 - x4230 * x6)
            + x38 * (x4223 - x4224 + x4229)
            - x6 * (x4230 * x5 + x4231 - x6 * (x4228 * x5 - x4233 * x6))
        )
    )
    result[9, 5] = numpy.sum(
        x1552
        * (
            x3 * (x4239 + x4241 * x5 - x4249 * x6)
            + x38 * (x4241 - x4242 + x4248)
            - x6 * (x4249 * x5 + x4250 - x6 * (x4247 * x5 - x4251 * x6))
        )
    )
    result[9, 6] = numpy.sum(
        x1428 * (x3 * (x4254 * x5 - x4260 * x6) + x4258 - x6 * (x4260 * x5 - x4262 * x6))
    )
    result[9, 7] = numpy.sum(
        x1603 * (x3 * (x4265 * x5 - x4270 * x6) + x4268 - x6 * (x4270 * x5 - x4272 * x6))
    )
    result[9, 8] = numpy.sum(
        x1603 * (x3 * (x4276 * x5 - x4279 * x6) + x4277 - x6 * (x4279 * x5 - x4281 * x6))
    )
    result[9, 9] = numpy.sum(
        x1428 * (x3 * (x4288 * x5 - x4290 * x6) + x4289 - x6 * (x4290 * x5 - x4292 * x6))
    )
    result[9, 10] = numpy.sum(x318 * (x3 * x4294 - x4296 * x6))
    result[9, 11] = numpy.sum(x1428 * (x3 * x4298 - x4300 * x6))
    result[9, 12] = numpy.sum(x1552 * (x3 * x4302 - x4304 * x6))
    result[9, 13] = numpy.sum(x1428 * (x3 * x4306 - x4308 * x6))
    result[9, 14] = numpy.sum(x318 * (x3 * x4309 - x4310 * x6))
    result[10, 0] = numpy.sum(
        x171
        * (
            x130 * (x4324 - x4325 - x4326 + x4334)
            + x5 * (x38 * (x4316 - x4317 + x4322) + x4324 * x5 - x4335 * x6)
            - x6
            * (
                x38 * (x4323 - x4327 + x4332)
                + x4335 * x5
                - x6
                * (
                    x4 * (x4321 - x4328 + x4329 - x4330)
                    + x4333 * x5
                    - x6 * (x4331 * x5 - x6 * (x1261 * x3605 + 3.0 * x2599 - x3615 * x8))
                )
            )
        )
    )
    result[10, 1] = numpy.sum(
        x318
        * (
            x38 * (x4341 - x4342 + x4347)
            + x5
            * (x4 * (x3613 + x4336 - x4337 - x4338 + x4339) + x4341 * x5 - x4348 * x6)
            - x6
            * (
                x4 * (x3616 + x4340 - x4343 + x4344 - x4345)
                + x4348 * x5
                - x6
                * (x4346 * x5 - x6 * (x1261 * x3635 + 3.0 * x2642 + x3614 - x3641 * x8))
            )
        )
    )
    result[10, 2] = numpy.sum(
        x318
        * (
            x38 * (x4354 - x4355 + x4360)
            + x5 * (x4 * (x4349 - x4350 - x4351 + x4352) + x4354 * x5 - x4361 * x6)
            - x6
            * (
                x4 * (x4353 - x4356 + x4357 - x4358)
                + x4361 * x5
                - x6 * (x4359 * x5 - x6 * (x1261 * x3663 + 3.0 * x2685 - x3669 * x8))
            )
        )
    )
    result[10, 3] = numpy.sum(
        x578
        * (
            x4 * (x3722 + x4362 - x4363 - x4364 + x4365)
            + x5 * (x4365 * x5 - x4366 * x6)
            - x6 * (x4366 * x5 - x6 * (x1261 * x3684 + 2.0 * x3640 - x3688 * x8 + x3728))
        )
    )
    result[10, 4] = numpy.sum(
        x687
        * (
            x4 * (x3670 + x4367 - x4368 - x4369 + x4370)
            + x5 * (x4370 * x5 - x4371 * x6)
            - x6 * (x4371 * x5 - x6 * (x1261 * x3698 + 3.0 * x2749 + x3668 - x3702 * x8))
        )
    )
    result[10, 5] = numpy.sum(
        x578
        * (
            x4 * (x4372 - x4373 - x4374 + x4375)
            + x5 * (x4375 * x5 - x4376 * x6)
            - x6 * (x4376 * x5 - x6 * (x1261 * x3714 + 3.0 * x2778 - x3718 * x8))
        )
    )
    result[10, 6] = numpy.sum(
        x318
        * (
            x5 * (x1261 * x3725 + 3.0 * x2793 + 3.0 * x3676 - x3727 * x8)
            - x6 * (x1261 * x3727 + 3.0 * x2795 + 3.0 * x3687 - x3729 * x8)
        )
    )
    result[10, 7] = numpy.sum(
        x687
        * (
            x5 * (x1261 * x3734 + 2.0 * x3692 - x3736 * x8 + x3753)
            - x6 * (x1261 * x3736 + 2.0 * x3701 - x3737 * x8 + x3755)
        )
    )
    result[10, 8] = numpy.sum(
        x687
        * (
            x5 * (x1261 * x3740 + 3.0 * x2829 + x3707 - x3742 * x8)
            - x6 * (x1261 * x3742 + 3.0 * x2831 + x3717 - x3743 * x8)
        )
    )
    result[10, 9] = numpy.sum(
        x318
        * (
            x5 * (x1261 * x3747 + 3.0 * x2844 - x3749 * x8)
            - x6 * (x1261 * x3749 + 3.0 * x2846 - x3750 * x8)
        )
    )
    result[10, 10] = numpy.sum(
        x171 * (x1261 * x3751 + 3.0 * x2852 + 4.0 * x3726 - x3752 * x8)
    )
    result[10, 11] = numpy.sum(
        x318 * (x1261 * x3754 + 3.0 * x2859 + 3.0 * x3735 - x3756 * x8)
    )
    result[10, 12] = numpy.sum(
        x578 * (x1261 * x3757 + 3.0 * x2868 + 2.0 * x3741 - x3758 * x8)
    )
    result[10, 13] = numpy.sum(x318 * (x1261 * x3759 + 3.0 * x2875 + x3748 - x3760 * x8))
    result[10, 14] = numpy.sum(x171 * (x1261 * x3761 + 3.0 * x2881 - x3762 * x8))
    result[11, 0] = numpy.sum(
        x318
        * (
            x130 * (x4390 - x4391 - x4392 + x4400)
            + x5 * (x38 * (x4382 - x4383 + x4388) + x4390 * x5 - x4401 * x6)
            - x6
            * (
                x38 * (x4389 - x4393 + x4398)
                + x4401 * x5
                - x6
                * (
                    x4 * (x4387 - x4394 + x4395 - x4396)
                    + x4399 * x5
                    - x6 * (x4397 * x5 - x6 * (x1261 * x3793 + 2.0 * x2940 - x3803 * x8))
                )
            )
        )
    )
    result[11, 1] = numpy.sum(
        x1428
        * (
            x38 * (x4406 - x4407 + x4411)
            + x5
            * (x4 * (x3801 + x3859 + x4402 - x4403 + x4404) + x4406 * x5 - x4412 * x6)
            - x6
            * (
                x4 * (x3804 + x3872 + x4405 - x4408 + x4409)
                + x4412 * x5
                - x6 * (x4410 * x5 - x6 * (x1261 * x3823 + x3802 - x3829 * x8 + x3874))
            )
        )
    )
    result[11, 2] = numpy.sum(
        x1428
        * (
            x38 * (x4418 - x4419 + x4424)
            + x5 * (x4 * (x4413 - x4414 - x4415 + x4416) + x4418 * x5 - x4425 * x6)
            - x6
            * (
                x4 * (x4417 - x4420 + x4421 - x4422)
                + x4425 * x5
                - x6 * (x4423 * x5 - x6 * (x1261 * x3848 + 2.0 * x3018 - x3854 * x8))
            )
        )
    )
    result[11, 3] = numpy.sum(
        -x1552
        * (
            x4 * (-x4426 + x4427 + x4428 + x4429 - x4430)
            - x5 * (x4430 * x5 - x4431 * x6)
            + x6
            * (x4431 * x5 - x6 * (x1261 * x3869 + 2.0 * x3043 + 2.0 * x3828 - x3875 * x8))
        )
    )
    result[11, 4] = numpy.sum(
        x1603
        * (
            x4 * (x3855 + x3914 + x4432 - x4433 + x4434)
            + x5 * (x4434 * x5 - x4435 * x6)
            - x6 * (x4435 * x5 - x6 * (x1261 * x3885 + x3853 - x3889 * x8 + x3919))
        )
    )
    result[11, 5] = numpy.sum(
        x1552
        * (
            x4 * (x4436 - x4437 - x4438 + x4439)
            + x5 * (x4439 * x5 - x4440 * x6)
            - x6 * (x4440 * x5 - x6 * (x1261 * x3899 + 2.0 * x3093 - x3903 * x8))
        )
    )
    result[11, 6] = numpy.sum(
        x1428
        * (
            x5 * (x1261 * x3907 + 2.0 * x3104 + 3.0 * x3862 - x3909 * x8)
            - x6 * (x1261 * x3909 + 2.0 * x3106 + 3.0 * x3873 - x3910 * x8)
        )
    )
    result[11, 7] = numpy.sum(
        x1603
        * (
            x5 * (x1261 * x3916 + 2.0 * x3117 + 2.0 * x3879 - x3918 * x8)
            - x6 * (x1261 * x3918 + 2.0 * x3119 + 2.0 * x3888 - x3920 * x8)
        )
    )
    result[11, 8] = numpy.sum(
        x1603
        * (
            x5 * (x1261 * x3923 + x3893 - x3925 * x8 + x3937)
            - x6 * (x1261 * x3925 + x3902 - x3926 * x8 + x3939)
        )
    )
    result[11, 9] = numpy.sum(
        x1428
        * (
            x5 * (x1261 * x3929 + 2.0 * x3143 - x3931 * x8)
            - x6 * (x1261 * x3931 + 2.0 * x3145 - x3932 * x8)
        )
    )
    result[11, 10] = numpy.sum(
        x318 * (x1261 * x3933 + 2.0 * x3151 + 4.0 * x3908 - x3934 * x8)
    )
    result[11, 11] = numpy.sum(
        x1428 * (x1261 * x3935 + 2.0 * x3158 + 3.0 * x3917 - x3936 * x8)
    )
    result[11, 12] = numpy.sum(
        x1552 * (x1261 * x3938 + 2.0 * x3165 + 2.0 * x3924 - x3940 * x8)
    )
    result[11, 13] = numpy.sum(x1428 * (x1261 * x3941 + 2.0 * x3171 + x3930 - x3942 * x8))
    result[11, 14] = numpy.sum(x318 * (x1261 * x3943 + 2.0 * x3177 - x3944 * x8))
    result[12, 0] = numpy.sum(
        x578
        * (
            x130 * (x4452 - x4453 - x4454 + x4461)
            + x5 * (x38 * (x4445 - x4446 + x4450) + x4452 * x5 - x4462 * x6)
            - x6
            * (
                x38 * (x4451 - x4455 + x4459)
                + x4462 * x5
                - x6
                * (
                    x4 * (x3243 + x4449 - x4456 + x4457)
                    + x4460 * x5
                    - x6 * (x4458 * x5 - x6 * (x1261 * x3971 + x3241 - x3981 * x8))
                )
            )
        )
    )
    result[12, 1] = numpy.sum(
        x1552
        * (
            x38 * (x4467 - x4468 + x4472)
            + x5
            * (x4 * (x3289 + x3979 + x4463 - x4464 + x4465) + x4467 * x5 - x4473 * x6)
            - x6
            * (
                x4 * (x3293 + x3982 + x4466 - x4469 + x4470)
                + x4473 * x5
                - x6 * (x4471 * x5 - x6 * (x1261 * x3998 + x3290 + x3980 - x4004 * x8))
            )
        )
    )
    result[12, 2] = numpy.sum(
        x1552
        * (
            x38 * (x4478 - x4479 + x4483)
            + x5 * (x4 * (x3334 + x4474 - x4475 + x4476) + x4478 * x5 - x4484 * x6)
            - x6
            * (
                x4 * (x3337 + x4477 - x4480 + x4481)
                + x4484 * x5
                - x6 * (x4482 * x5 - x6 * (x1261 * x4020 + x3335 - x4026 * x8))
            )
        )
    )
    result[12, 3] = numpy.sum(
        x2723
        * (
            x4 * (x3367 + x4485 - x4486 - x4487 + x4488)
            + x5 * (x4488 * x5 - x4489 * x6)
            - x6 * (x4489 * x5 - x6 * (x1261 * x4035 + x3371 + 2.0 * x4003 - x4039 * x8))
        )
    )
    result[12, 4] = numpy.sum(
        x2751
        * (
            x4 * (x3401 + x4027 + x4490 - x4491 + x4492)
            + x5 * (x4492 * x5 - x4493 * x6)
            - x6 * (x4493 * x5 - x6 * (x1261 * x4047 + x3404 + x4025 - x4051 * x8))
        )
    )
    result[12, 5] = numpy.sum(
        x2723
        * (
            x4 * (x3436 + x4494 - x4495 + x4496)
            + x5 * (x4496 * x5 - x4497 * x6)
            - x6 * (x4497 * x5 - x6 * (x1261 * x4059 + x3440 - x4063 * x8))
        )
    )
    result[12, 6] = numpy.sum(
        x1552
        * (
            x5 * (x1261 * x4065 + x3460 + 3.0 * x4030 - x4067 * x8)
            - x6 * (x1261 * x4067 + x3463 + 3.0 * x4038 - x4068 * x8)
        )
    )
    result[12, 7] = numpy.sum(
        x2751
        * (
            x5 * (x1261 * x4070 + x3481 + 2.0 * x4042 - x4072 * x8)
            - x6 * (x1261 * x4072 + x3484 + 2.0 * x4050 - x4073 * x8)
        )
    )
    result[12, 8] = numpy.sum(
        x2751
        * (
            x5 * (x1261 * x4075 + x3501 + x4054 - x4077 * x8)
            - x6 * (x1261 * x4077 + x3503 + x4062 - x4078 * x8)
        )
    )
    result[12, 9] = numpy.sum(
        x1552
        * (
            x5 * (x1261 * x4080 + x3519 - x4082 * x8)
            - x6 * (x1261 * x4082 + x3521 - x4083 * x8)
        )
    )
    result[12, 10] = numpy.sum(x578 * (x1261 * x4084 + x3529 + 4.0 * x4066 - x4085 * x8))
    result[12, 11] = numpy.sum(x1552 * (x1261 * x4086 + x3540 + 3.0 * x4071 - x4087 * x8))
    result[12, 12] = numpy.sum(x2723 * (x1261 * x4088 + x3550 + 2.0 * x4076 - x4089 * x8))
    result[12, 13] = numpy.sum(x1552 * (x1261 * x4090 + x3559 + x4081 - x4091 * x8))
    result[12, 14] = numpy.sum(x578 * (x1261 * x4092 + x3568 - x4093 * x8))
    result[13, 0] = numpy.sum(
        x318
        * (
            x130 * (x4507 - x4508 - x4509 + x4515)
            + x5 * (x38 * (x4501 - x4502 + x4505) + x4507 * x5 - x4516 * x6)
            - x6
            * (
                x38 * (x4506 - x4510 + x4513)
                + x4516 * x5
                - x6
                * (
                    x4 * (x4160 + x4504 - x4511)
                    + x4514 * x5
                    - x6 * (x4512 * x5 - x6 * (x1261 * x4128 + x4168))
                )
            )
        )
    )
    result[13, 1] = numpy.sum(
        x1428
        * (
            x38 * (x4520 - x4521 + x4524)
            + x5 * (x4 * (x4200 - x4517 + x4518) + x4520 * x5 - x4525 * x6)
            - x6
            * (
                x4 * (x4212 + x4519 - x4522)
                + x4525 * x5
                - x6 * (x4523 * x5 - x6 * (x1261 * x4162 + x4214))
            )
        )
    )
    result[13, 2] = numpy.sum(
        x1428
        * (
            x38 * (x4529 - x4530 + x4533)
            + x5 * (x4 * (x4216 - x4526 + x4527) + x4529 * x5 - x4534 * x6)
            - x6
            * (
                x4 * (x4226 + x4528 - x4531)
                + x4534 * x5
                - x6 * (x4532 * x5 - x6 * (x1261 * x4188 + x4232))
            )
        )
    )
    result[13, 3] = numpy.sum(
        x1552
        * (
            x4 * (x4257 - x4535 + x4536)
            + x5 * (x4536 * x5 - x4537 * x6)
            - x6 * (x4537 * x5 - x6 * (x1261 * x4209 + x4261))
        )
    )
    result[13, 4] = numpy.sum(
        x1603
        * (
            x4 * (x4267 - x4538 + x4539)
            + x5 * (x4539 * x5 - x4540 * x6)
            - x6 * (x4540 * x5 - x6 * (x1261 * x4228 + x4271))
        )
    )
    result[13, 5] = numpy.sum(
        x1552
        * (
            x4 * (x4273 - x4541 + x4542)
            + x5 * (x4542 * x5 - x4543 * x6)
            - x6 * (x4543 * x5 - x6 * (x1261 * x4247 + x4280))
        )
    )
    result[13, 6] = numpy.sum(
        x1428 * (x5 * (x1261 * x4254 + x4293) - x6 * (x1261 * x4260 + x4295))
    )
    result[13, 7] = numpy.sum(
        x1603 * (x5 * (x1261 * x4265 + x4297) - x6 * (x1261 * x4270 + x4299))
    )
    result[13, 8] = numpy.sum(
        x1603 * (x5 * (x1261 * x4276 + x4301) - x6 * (x1261 * x4279 + x4303))
    )
    result[13, 9] = numpy.sum(
        x1428 * (x5 * (x1261 * x4288 + x4305) - x6 * (x1261 * x4290 + x4307))
    )
    result[13, 10] = numpy.sum(x318 * (x1261 * x4294 + 4.0 * x4258 - x4296 * x8))
    result[13, 11] = numpy.sum(x1428 * (x1261 * x4298 + 3.0 * x4268 - x4300 * x8))
    result[13, 12] = numpy.sum(x1552 * (x1261 * x4302 + 2.0 * x4277 - x4304 * x8))
    result[13, 13] = numpy.sum(x1428 * (x1261 * x4306 + x4289 - x4308 * x8))
    result[13, 14] = numpy.sum(x318 * (x1261 * x4309 - x4310 * x8))
    result[14, 0] = numpy.sum(
        x171
        * (
            x130 * (x4558 - x4560 + x4568 + x4572)
            + x5 * (x38 * (x4550 - x4551 + x4556) + x4558 * x5 - x4569 * x6)
            - x6
            * (
                x38 * (x4557 - x4561 + x4566)
                + x4569 * x5
                - x6 * (x4567 * x5 + x4570 - x6 * (x4565 * x5 - x4571 * x6))
            )
        )
    )
    result[14, 1] = numpy.sum(
        x318
        * (
            x38 * (x4578 - x4579 + x4583)
            + x5 * (x4576 + x4578 * x5 - x4584 * x6)
            - x6 * (x4584 * x5 + x4585 - x6 * (x4582 * x5 - x4586 * x6))
        )
    )
    result[14, 2] = numpy.sum(
        x318
        * (
            x38 * (x4593 - x4594 + x4599)
            + x5 * (x4591 + x4593 * x5 - x4600 * x6)
            - x6 * (x4600 * x5 + x4601 - x6 * (x4598 * x5 - x4602 * x6))
        )
    )
    result[14, 3] = numpy.sum(
        x578 * (x4606 + x5 * (x4605 * x5 - x4607 * x6) - x6 * (x4607 * x5 - x4608 * x6))
    )
    result[14, 4] = numpy.sum(
        x687 * (x4612 + x5 * (x4611 * x5 - x4613 * x6) - x6 * (x4613 * x5 - x4614 * x6))
    )
    result[14, 5] = numpy.sum(
        x578 * (x4619 + x5 * (x4618 * x5 - x4620 * x6) - x6 * (x4620 * x5 - x4621 * x6))
    )
    result[14, 6] = numpy.sum(x318 * (x4622 * x5 - x4623 * x6))
    result[14, 7] = numpy.sum(x687 * (x4624 * x5 - x4625 * x6))
    result[14, 8] = numpy.sum(x687 * (x4626 * x5 - x4627 * x6))
    result[14, 9] = numpy.sum(x318 * (x4628 * x5 - x4629 * x6))
    result[14, 10] = numpy.sum(x171 * (x174 * x4622 + 3.0 * x4606 - x4623 * x8))
    result[14, 11] = numpy.sum(x318 * (x174 * x4624 + 2.0 * x4612 - x4625 * x8))
    result[14, 12] = numpy.sum(x578 * (x174 * x4626 + x4619 - x4627 * x8))
    result[14, 13] = numpy.sum(x318 * (x174 * x4628 - x4629 * x8))
    result[14, 14] = numpy.sum(
        x171 * (-x10 * x4310 + x1852 * x4309 + 3.0 * x3568 + 4.0 * x4289)
    )
    return result


coulomb3d = {
    (0, 0): coulomb3d_00,
    (0, 1): coulomb3d_01,
    (0, 2): coulomb3d_02,
    (0, 3): coulomb3d_03,
    (0, 4): coulomb3d_04,
    (1, 0): coulomb3d_10,
    (1, 1): coulomb3d_11,
    (1, 2): coulomb3d_12,
    (1, 3): coulomb3d_13,
    (1, 4): coulomb3d_14,
    (2, 0): coulomb3d_20,
    (2, 1): coulomb3d_21,
    (2, 2): coulomb3d_22,
    (2, 3): coulomb3d_23,
    (2, 4): coulomb3d_24,
    (3, 0): coulomb3d_30,
    (3, 1): coulomb3d_31,
    (3, 2): coulomb3d_32,
    (3, 3): coulomb3d_33,
    (3, 4): coulomb3d_34,
    (4, 0): coulomb3d_40,
    (4, 1): coulomb3d_41,
    (4, 2): coulomb3d_42,
    (4, 3): coulomb3d_43,
    (4, 4): coulomb3d_44,
}
