import numpy


_L_MAX = 4


def ovlp3d_00(ax, da, A, bx, db, B):
    """Cartesian 3D (ss) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 1), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = ax * bx * x0

    # 1 item(s)
    result[0, 0] = numpy.sum(
        2.82842712474619
        * da
        * db
        * x0**1.5
        * numpy.sqrt(ax**1.5)
        * numpy.sqrt(bx**1.5)
        * numpy.exp(-x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2))
    )
    return result


def ovlp3d_01(ax, da, A, bx, db, B):
    """Cartesian 3D (sp) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 3), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = ax * bx * x0
    x2 = (
        5.65685424949238
        * da
        * db
        * x0**1.5
        * numpy.sqrt(ax**1.5)
        * numpy.sqrt(bx**2.5)
        * numpy.exp(-x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2))
    )

    # 3 item(s)
    result[0, 0] = numpy.sum(x2 * (x0 * (ax * A[0] + bx * B[0]) - B[0]))
    result[0, 1] = numpy.sum(x2 * (x0 * (ax * A[1] + bx * B[1]) - B[1]))
    result[0, 2] = numpy.sum(x2 * (x0 * (ax * A[2] + bx * B[2]) - B[2]))
    return result


def ovlp3d_02(ax, da, A, bx, db, B):
    """Cartesian 3D (sd) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 6), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = x0 * (ax * A[0] + bx * B[0]) - B[0]
    x2 = ax * bx * x0
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(x0)
    x5 = 1.77245385090552 * x4
    x6 = x3 * x5
    x7 = 0.5 / (ax + bx)
    x8 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x9 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x10 = numpy.sqrt(ax**1.5)
    x11 = numpy.sqrt(bx**3.5)
    x12 = da * db * x0 * x10 * x11 * x8 * x9
    x13 = 3.68527092769425
    x14 = x0 * (ax * A[1] + bx * B[1]) - B[1]
    x15 = 11.3137084989848 * x1 * x12 * x3 * x4
    x16 = x0 * (ax * A[2] + bx * B[2]) - B[2]
    x17 = x5 * x8
    x18 = da * db * x0 * x10 * x11 * x13 * x3
    x19 = x5 * x9

    # 6 item(s)
    result[0, 0] = numpy.sum(x12 * x13 * x6 * (x1**2 + x7))
    result[0, 1] = numpy.sum(x14 * x15)
    result[0, 2] = numpy.sum(x15 * x16)
    result[0, 3] = numpy.sum(x17 * x18 * x9 * (x14**2 + x7))
    result[0, 4] = numpy.sum(11.3137084989848 * x12 * x14 * x16 * x3 * x4)
    result[0, 5] = numpy.sum(x18 * x19 * x8 * (x16**2 + x7))
    return result


def ovlp3d_03(ax, da, A, bx, db, B):
    """Cartesian 3D (sf) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 10), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = ax * bx * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = x0 * (ax * A[0] + bx * B[0]) - B[0]
    x4 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x5 = numpy.sqrt(x0)
    x6 = 1.77245385090552 * x5
    x7 = x4 * x6
    x8 = 0.5 / (ax + bx)
    x9 = x7 * x8
    x10 = x3**2 * x7 + x9
    x11 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x12 = 0.564189583547756
    x13 = numpy.sqrt(ax**1.5)
    x14 = numpy.sqrt(bx**4.5)
    x15 = da * db * x0 * x11 * x12 * x13 * x14
    x16 = 5.84237394672177 * x15
    x17 = x0 * (ax * A[1] + bx * B[1]) - B[1]
    x18 = 13.0639452948436 * x15
    x19 = x10 * x18 * x2
    x20 = x0 * (ax * A[2] + bx * B[2]) - B[2]
    x21 = x2 * x6
    x22 = x21 * x8
    x23 = x17**2 * x21 + x22
    x24 = x18 * x23 * x4
    x25 = da * db * x0 * x13 * x14 * x2 * x4
    x26 = x11 * x6
    x27 = x26 * x8
    x28 = x20**2 * x26 + x27
    x29 = x12 * x25
    x30 = 13.0639452948436 * x28 * x29

    # 10 item(s)
    result[0, 0] = numpy.sum(x16 * x2 * x3 * (x10 + 2.0 * x9))
    result[0, 1] = numpy.sum(x17 * x19)
    result[0, 2] = numpy.sum(x19 * x20)
    result[0, 3] = numpy.sum(x24 * x3)
    result[0, 4] = numpy.sum(22.6274169979695 * x11 * x17 * x20 * x25 * x3 * x5)
    result[0, 5] = numpy.sum(x3 * x30)
    result[0, 6] = numpy.sum(x16 * x17 * x4 * (2.0 * x22 + x23))
    result[0, 7] = numpy.sum(x20 * x24)
    result[0, 8] = numpy.sum(x17 * x30)
    result[0, 9] = numpy.sum(5.84237394672177 * x20 * x29 * (2.0 * x27 + x28))
    return result


def ovlp3d_04(ax, da, A, bx, db, B):
    """Cartesian 3D (sg) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 15), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = ax * bx * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = 0.5 / (ax + bx)
    x4 = x0 * (ax * A[0] + bx * B[0]) - B[0]
    x5 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x0)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x4**2 * x8
    x10 = x3 * x8
    x11 = x10 + x9
    x12 = x4 * (2.0 * x10 + x11)
    x13 = 0.564189583547756
    x14 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x15 = da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**5.5)
    x16 = x14 * x15
    x17 = x0 * x13 * x16
    x18 = 4.41641957979107 * x17
    x19 = x0 * (ax * A[1] + bx * B[1]) - B[1]
    x20 = x19 * x2
    x21 = 11.6847478934435 * x17
    x22 = x12 * x21
    x23 = x0 * (ax * A[2] + bx * B[2]) - B[2]
    x24 = x2 * x7
    x25 = x19**2 * x24
    x26 = x24 * x3
    x27 = x25 + x26
    x28 = 0.318309886183791
    x29 = 15.084944665313 * x11 * x28 * x6
    x30 = 26.1278905896872 * x17 * x23
    x31 = x14 * x7
    x32 = x23**2 * x31
    x33 = x3 * x31
    x34 = x32 + x33
    x35 = x15 * x34
    x36 = x2 * x35
    x37 = x19 * (2.0 * x26 + x27)
    x38 = x21 * x37 * x5
    x39 = x27 * x5
    x40 = x0 * x13 * x5
    x41 = x23 * (2.0 * x33 + x34)
    x42 = x15 * x2 * x40
    x43 = 11.6847478934435 * x41 * x42

    # 15 item(s)
    result[0, 0] = numpy.sum(x18 * x2 * (x12 * x4 + 3.0 * x3 * (x10 + x9)))
    result[0, 1] = numpy.sum(x20 * x22)
    result[0, 2] = numpy.sum(x2 * x22 * x23)
    result[0, 3] = numpy.sum(x16 * x27 * x29)
    result[0, 4] = numpy.sum(x11 * x20 * x30)
    result[0, 5] = numpy.sum(x29 * x36)
    result[0, 6] = numpy.sum(x38 * x4)
    result[0, 7] = numpy.sum(x30 * x39 * x4)
    result[0, 8] = numpy.sum(26.1278905896872 * x19 * x36 * x4 * x40)
    result[0, 9] = numpy.sum(x4 * x43)
    result[0, 10] = numpy.sum(x18 * x5 * (x19 * x37 + 3.0 * x3 * (x25 + x26)))
    result[0, 11] = numpy.sum(x23 * x38)
    result[0, 12] = numpy.sum(15.084944665313 * x28 * x35 * x39 * x6)
    result[0, 13] = numpy.sum(x19 * x43)
    result[0, 14] = numpy.sum(
        4.41641957979107 * x42 * (x23 * x41 + 3.0 * x3 * (x32 + x33))
    )
    return result


def ovlp3d_10(ax, da, A, bx, db, B):
    """Cartesian 3D (ps) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 1), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = ax * bx * x0
    x2 = (
        5.65685424949238
        * da
        * db
        * x0**1.5
        * numpy.sqrt(ax**2.5)
        * numpy.sqrt(bx**1.5)
        * numpy.exp(-x1 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2))
    )

    # 3 item(s)
    result[0, 0] = numpy.sum(x2 * (x0 * (ax * A[0] + bx * B[0]) - A[0]))
    result[1, 0] = numpy.sum(x2 * (x0 * (ax * A[1] + bx * B[1]) - A[1]))
    result[2, 0] = numpy.sum(x2 * (x0 * (ax * A[2] + bx * B[2]) - A[2]))
    return result


def ovlp3d_11(ax, da, A, bx, db, B):
    """Cartesian 3D (pp) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 3), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = ax * bx * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(x0)
    x5 = 1.77245385090552 * x4
    x6 = 0.5 * x5 / (ax + bx)
    x7 = -x0 * (ax * A[0] + bx * B[0])
    x8 = -x7 - B[0]
    x9 = -x3 * (x7 + A[0])
    x10 = 0.564189583547756
    x11 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x12 = numpy.sqrt(ax**2.5)
    x13 = numpy.sqrt(bx**2.5)
    x14 = 11.3137084989848 * da * db * x0 * x11 * x12 * x13
    x15 = x10 * x14
    x16 = x4 * x9
    x17 = -x0 * (ax * A[1] + bx * B[1])
    x18 = -x2 * (x17 + B[1])
    x19 = x14 * x18
    x20 = -x0 * (ax * A[2] + bx * B[2])
    x21 = -x11 * (x20 + B[2])
    x22 = 11.3137084989848 * da * db * x0 * x12 * x13
    x23 = -x17 - A[1]
    x24 = x23 * x4
    x25 = x2 * x3
    x26 = x14 * x25 * x8
    x27 = x22 * x25
    x28 = -x20 - A[2]
    x29 = x28 * x4

    # 9 item(s)
    result[0, 0] = numpy.sum(x15 * x2 * (x3 * x6 + x5 * x8 * x9))
    result[0, 1] = numpy.sum(x16 * x19)
    result[0, 2] = numpy.sum(x16 * x2 * x21 * x22)
    result[1, 0] = numpy.sum(x24 * x26)
    result[1, 1] = numpy.sum(x15 * x3 * (x18 * x23 * x5 + x2 * x6))
    result[1, 2] = numpy.sum(x21 * x24 * x27)
    result[2, 0] = numpy.sum(x26 * x29)
    result[2, 1] = numpy.sum(x19 * x29 * x3)
    result[2, 2] = numpy.sum(x10 * x27 * (x11 * x6 + x21 * x28 * x5))
    return result


def ovlp3d_12(ax, da, A, bx, db, B):
    """Cartesian 3D (pd) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 6), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = -x0 * (ax * A[0] + bx * B[0])
    x2 = -x1 - A[0]
    x3 = -x1 - B[0]
    x4 = ax * bx * x0
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x0)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = 0.5 / (ax + bx)
    x10 = x8 * x9
    x11 = x10 + x3**2 * x8
    x12 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x13 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x14 = 0.564189583547756
    x15 = numpy.sqrt(ax**2.5)
    x16 = numpy.sqrt(bx**3.5)
    x17 = da * db * x0 * x13 * x14 * x15 * x16
    x18 = 13.0639452948436 * x17
    x19 = x12 * x18
    x20 = -x0 * (ax * A[1] + bx * B[1])
    x21 = -x20 - B[1]
    x22 = 22.6274169979695 * x17
    x23 = x12 * x22 * (x10 + x2 * x3 * x8)
    x24 = -x0 * (ax * A[2] + bx * B[2])
    x25 = -x24 - B[2]
    x26 = x12 * x7
    x27 = x26 * x9
    x28 = x21**2 * x26 + x27
    x29 = x18 * x5
    x30 = x28 * x29
    x31 = 22.6274169979695 * x21
    x32 = da * db * x0 * x12 * x15 * x16 * x5
    x33 = x13 * x32 * x6
    x34 = x25 * x33
    x35 = x13 * x7
    x36 = x35 * x9
    x37 = x25**2 * x35 + x36
    x38 = x14 * x32
    x39 = 13.0639452948436 * x38
    x40 = x37 * x39
    x41 = -x20 - A[1]
    x42 = x11 * x19
    x43 = x22 * x5 * (x21 * x26 * x41 + x27)
    x44 = 22.6274169979695 * x3
    x45 = -x24 - A[2]
    x46 = x38 * (x25 * x35 * x45 + x36)

    # 18 item(s)
    result[0, 0] = numpy.sum(x19 * (2.0 * x10 * x3 + x11 * x2))
    result[0, 1] = numpy.sum(x21 * x23)
    result[0, 2] = numpy.sum(x23 * x25)
    result[0, 3] = numpy.sum(x2 * x30)
    result[0, 4] = numpy.sum(x2 * x31 * x34)
    result[0, 5] = numpy.sum(x2 * x40)
    result[1, 0] = numpy.sum(x41 * x42)
    result[1, 1] = numpy.sum(x3 * x43)
    result[1, 2] = numpy.sum(x34 * x41 * x44)
    result[1, 3] = numpy.sum(x29 * (2.0 * x21 * x27 + x28 * x41))
    result[1, 4] = numpy.sum(x25 * x43)
    result[1, 5] = numpy.sum(x40 * x41)
    result[2, 0] = numpy.sum(x42 * x45)
    result[2, 1] = numpy.sum(x21 * x33 * x44 * x45)
    result[2, 2] = numpy.sum(x44 * x46)
    result[2, 3] = numpy.sum(x30 * x45)
    result[2, 4] = numpy.sum(x31 * x46)
    result[2, 5] = numpy.sum(x39 * (2.0 * x25 * x36 + x37 * x45))
    return result


def ovlp3d_13(ax, da, A, bx, db, B):
    """Cartesian 3D (pf) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 10), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = -x2 - A[0]
    x12 = x10 + x9
    x13 = 2.0 * x10 * x3
    x14 = x12 * x3 + x13
    x15 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x16 = 0.564189583547756
    x17 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x18 = da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**4.5)
    x19 = x17 * x18
    x20 = x1 * x16 * x19
    x21 = 11.6847478934435 * x20
    x22 = x15 * x21
    x23 = -x1 * (ax * A[1] + bx * B[1])
    x24 = -x23 - B[1]
    x25 = x15 * x24
    x26 = 26.1278905896872 * x20
    x27 = x26 * (x11 * x12 + x13)
    x28 = -x1 * (ax * A[2] + bx * B[2])
    x29 = -x28 - B[2]
    x30 = x15 * x29
    x31 = x15 * x7
    x32 = x24**2 * x31
    x33 = x0 * x31
    x34 = x32 + x33
    x35 = 26.1278905896872 * x34
    x36 = x10 + x11 * x3 * x8
    x37 = 0.318309886183791 * x6
    x38 = x36 * x37
    x39 = x20 * x25
    x40 = 45.2548339959391 * x29
    x41 = x15 * x18
    x42 = x17 * x7
    x43 = x29**2 * x42
    x44 = x0 * x42
    x45 = x43 + x44
    x46 = 26.1278905896872 * x45
    x47 = x41 * x46
    x48 = 2.0 * x24 * x33
    x49 = x24 * x34 + x48
    x50 = x21 * x5
    x51 = x49 * x50
    x52 = x29 * x5
    x53 = x20 * x35
    x54 = x1 * x16 * x5
    x55 = x47 * x54
    x56 = 2.0 * x29 * x44
    x57 = x29 * x45 + x56
    x58 = x41 * x54
    x59 = 11.6847478934435 * x58
    x60 = x57 * x59
    x61 = -x23 - A[1]
    x62 = x14 * x22
    x63 = x24 * x31 * x61 + x33
    x64 = x37 * x63
    x65 = 26.1278905896872 * x12
    x66 = x3 * x5
    x67 = x26 * (x34 * x61 + x48)
    x68 = x18 * x5
    x69 = -x28 - A[2]
    x70 = x29 * x42 * x69 + x44
    x71 = x37 * x70
    x72 = x3 * x58
    x73 = 26.1278905896872 * x45 * x69 + 26.1278905896872 * x56

    # 30 item(s)
    result[0, 0] = numpy.sum(x22 * (3.0 * x0 * (x10 + x9) + x11 * x14))
    result[0, 1] = numpy.sum(x25 * x27)
    result[0, 2] = numpy.sum(x27 * x30)
    result[0, 3] = numpy.sum(x19 * x35 * x38)
    result[0, 4] = numpy.sum(x36 * x39 * x40)
    result[0, 5] = numpy.sum(x38 * x47)
    result[0, 6] = numpy.sum(x11 * x51)
    result[0, 7] = numpy.sum(x11 * x52 * x53)
    result[0, 8] = numpy.sum(x11 * x24 * x55)
    result[0, 9] = numpy.sum(x11 * x60)
    result[1, 0] = numpy.sum(x61 * x62)
    result[1, 1] = numpy.sum(x19 * x64 * x65)
    result[1, 2] = numpy.sum(x20 * x30 * x61 * x65)
    result[1, 3] = numpy.sum(x66 * x67)
    result[1, 4] = numpy.sum(x20 * x40 * x63 * x66)
    result[1, 5] = numpy.sum(x3 * x55 * x61)
    result[1, 6] = numpy.sum(x50 * (3.0 * x0 * (x32 + x33) + x49 * x61))
    result[1, 7] = numpy.sum(x52 * x67)
    result[1, 8] = numpy.sum(x46 * x64 * x68)
    result[1, 9] = numpy.sum(x60 * x61)
    result[2, 0] = numpy.sum(x62 * x69)
    result[2, 1] = numpy.sum(x39 * x65 * x69)
    result[2, 2] = numpy.sum(x41 * x65 * x71)
    result[2, 3] = numpy.sum(x53 * x66 * x69)
    result[2, 4] = numpy.sum(45.2548339959391 * x24 * x70 * x72)
    result[2, 5] = numpy.sum(x72 * x73)
    result[2, 6] = numpy.sum(x51 * x69)
    result[2, 7] = numpy.sum(x35 * x68 * x71)
    result[2, 8] = numpy.sum(x24 * x58 * x73)
    result[2, 9] = numpy.sum(x59 * (3.0 * x0 * (x43 + x44) + x57 * x69))
    return result


def ovlp3d_14(ax, da, A, bx, db, B):
    """Cartesian 3D (pg) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 15), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = x10 + x9
    x12 = x11 * x3
    x13 = x10 * x3
    x14 = -x2 - A[0]
    x15 = 3.0 * x0 * (x10 + x9)
    x16 = 2.0 * x13
    x17 = x12 + x16
    x18 = x15 + x17 * x3
    x19 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x20 = 0.564189583547756
    x21 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x22 = da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**5.5)
    x23 = x21 * x22
    x24 = x1 * x20 * x23
    x25 = 8.83283915958214 * x24
    x26 = x19 * x25
    x27 = -x1 * (ax * A[1] + bx * B[1])
    x28 = -x27 - B[1]
    x29 = x19 * x28
    x30 = 23.3694957868871 * x24
    x31 = x30 * (x14 * x17 + x15)
    x32 = -x1 * (ax * A[2] + bx * B[2])
    x33 = -x32 - B[2]
    x34 = x19 * x33
    x35 = x11 * x14 + x16
    x36 = 30.169889330626 * x35
    x37 = x19 * x7
    x38 = x28**2 * x37
    x39 = x0 * x37
    x40 = x38 + x39
    x41 = 0.318309886183791 * x6
    x42 = x23 * x41
    x43 = x40 * x42
    x44 = 52.2557811793745 * x24 * x33
    x45 = x21 * x7
    x46 = x33**2 * x45
    x47 = x0 * x45
    x48 = x46 + x47
    x49 = x22 * x41
    x50 = x19 * x49
    x51 = x48 * x50
    x52 = x28 * x40
    x53 = x28 * x39
    x54 = 2.0 * x53
    x55 = x52 + x54
    x56 = x10 + x14 * x3 * x8
    x57 = 23.3694957868871 * x56
    x58 = 52.2557811793745 * x56
    x59 = x33 * x48
    x60 = x33 * x47
    x61 = 2.0 * x60
    x62 = x59 + x61
    x63 = 3.0 * x0 * (x38 + x39)
    x64 = x28 * x55 + x63
    x65 = x25 * x5
    x66 = x64 * x65
    x67 = x30 * x5
    x68 = x55 * x67
    x69 = x49 * x5
    x70 = x48 * x69
    x71 = x1 * x19 * x20 * x22 * x5
    x72 = 23.3694957868871 * x71
    x73 = x62 * x72
    x74 = 3.0 * x0 * (x46 + x47)
    x75 = x33 * x62 + x74
    x76 = 8.83283915958214 * x71
    x77 = x75 * x76
    x78 = -x27 - A[1]
    x79 = x18 * x26
    x80 = x28 * x37 * x78 + x39
    x81 = 23.3694957868871 * x80
    x82 = x17 * x30
    x83 = x40 * x78 + x54
    x84 = 30.169889330626 * x83
    x85 = x11 * x42
    x86 = 52.2557811793745 * x80
    x87 = 30.169889330626 * x11
    x88 = x67 * (x55 * x78 + x63)
    x89 = -x32 - A[2]
    x90 = x33 * x45 * x89 + x47
    x91 = 23.3694957868871 * x90
    x92 = x11 * x50
    x93 = 52.2557811793745 * x90
    x94 = x48 * x89 + x61
    x95 = 30.169889330626 * x94
    x96 = x40 * x69
    x97 = x72 * (x62 * x89 + x74)

    # 45 item(s)
    result[0, 0] = numpy.sum(x26 * (4.0 * x0 * (x12 + 2.0 * x13) + x14 * x18))
    result[0, 1] = numpy.sum(x29 * x31)
    result[0, 2] = numpy.sum(x31 * x34)
    result[0, 3] = numpy.sum(x36 * x43)
    result[0, 4] = numpy.sum(x29 * x35 * x44)
    result[0, 5] = numpy.sum(x36 * x51)
    result[0, 6] = numpy.sum(x42 * x55 * x57)
    result[0, 7] = numpy.sum(x33 * x43 * x58)
    result[0, 8] = numpy.sum(x28 * x51 * x58)
    result[0, 9] = numpy.sum(x50 * x57 * x62)
    result[0, 10] = numpy.sum(x14 * x66)
    result[0, 11] = numpy.sum(x14 * x33 * x68)
    result[0, 12] = numpy.sum(30.169889330626 * x14 * x40 * x70)
    result[0, 13] = numpy.sum(x14 * x28 * x73)
    result[0, 14] = numpy.sum(x14 * x77)
    result[1, 0] = numpy.sum(x78 * x79)
    result[1, 1] = numpy.sum(x17 * x42 * x81)
    result[1, 2] = numpy.sum(x34 * x78 * x82)
    result[1, 3] = numpy.sum(x84 * x85)
    result[1, 4] = numpy.sum(x33 * x85 * x86)
    result[1, 5] = numpy.sum(x51 * x78 * x87)
    result[1, 6] = numpy.sum(x3 * x88)
    result[1, 7] = numpy.sum(x3 * x44 * x5 * x83)
    result[1, 8] = numpy.sum(x3 * x70 * x86)
    result[1, 9] = numpy.sum(x3 * x73 * x78)
    result[1, 10] = numpy.sum(x65 * (4.0 * x0 * (x52 + 2.0 * x53) + x64 * x78))
    result[1, 11] = numpy.sum(x33 * x88)
    result[1, 12] = numpy.sum(x70 * x84)
    result[1, 13] = numpy.sum(x62 * x69 * x81)
    result[1, 14] = numpy.sum(x77 * x78)
    result[2, 0] = numpy.sum(x79 * x89)
    result[2, 1] = numpy.sum(x29 * x82 * x89)
    result[2, 2] = numpy.sum(x17 * x50 * x91)
    result[2, 3] = numpy.sum(x43 * x87 * x89)
    result[2, 4] = numpy.sum(x28 * x92 * x93)
    result[2, 5] = numpy.sum(x92 * x95)
    result[2, 6] = numpy.sum(x3 * x68 * x89)
    result[2, 7] = numpy.sum(x3 * x93 * x96)
    result[2, 8] = numpy.sum(52.2557811793745 * x28 * x3 * x71 * x94)
    result[2, 9] = numpy.sum(x3 * x97)
    result[2, 10] = numpy.sum(x66 * x89)
    result[2, 11] = numpy.sum(x55 * x69 * x91)
    result[2, 12] = numpy.sum(x95 * x96)
    result[2, 13] = numpy.sum(x28 * x97)
    result[2, 14] = numpy.sum(x76 * (4.0 * x0 * (x59 + 2.0 * x60) + x75 * x89))
    return result


def ovlp3d_20(ax, da, A, bx, db, B):
    """Cartesian 3D (ds) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 1), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = x0 * (ax * A[0] + bx * B[0]) - A[0]
    x2 = ax * bx * x0
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(x0)
    x5 = 1.77245385090552 * x4
    x6 = x3 * x5
    x7 = 0.5 / (ax + bx)
    x8 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x9 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x10 = numpy.sqrt(ax**3.5)
    x11 = numpy.sqrt(bx**1.5)
    x12 = da * db * x0 * x10 * x11 * x8 * x9
    x13 = 3.68527092769425
    x14 = x0 * (ax * A[1] + bx * B[1]) - A[1]
    x15 = 11.3137084989848 * x1 * x12 * x3 * x4
    x16 = x0 * (ax * A[2] + bx * B[2]) - A[2]
    x17 = x5 * x8
    x18 = da * db * x0 * x10 * x11 * x13 * x3
    x19 = x5 * x9

    # 6 item(s)
    result[0, 0] = numpy.sum(x12 * x13 * x6 * (x1**2 + x7))
    result[1, 0] = numpy.sum(x14 * x15)
    result[2, 0] = numpy.sum(x15 * x16)
    result[3, 0] = numpy.sum(x17 * x18 * x9 * (x14**2 + x7))
    result[4, 0] = numpy.sum(11.3137084989848 * x12 * x14 * x16 * x3 * x4)
    result[5, 0] = numpy.sum(x18 * x19 * x8 * (x16**2 + x7))
    return result


def ovlp3d_21(ax, da, A, bx, db, B):
    """Cartesian 3D (dp) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 3), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = numpy.sqrt(x1)
    x3 = 1.77245385090552 * x2
    x4 = -x1 * (ax * A[0] + bx * B[0])
    x5 = -x4 - A[0]
    x6 = ax * bx * x1
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x4 - B[0]
    x11 = x3 * x7
    x12 = x0 * x11
    x13 = x10 * x9 + x12
    x14 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x15 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x16 = 0.564189583547756
    x17 = numpy.sqrt(ax**3.5)
    x18 = numpy.sqrt(bx**2.5)
    x19 = da * db * x1 * x15 * x16 * x17 * x18
    x20 = 13.0639452948436 * x19
    x21 = x14 * x20
    x22 = -x1 * (ax * A[1] + bx * B[1])
    x23 = -x22 - B[1]
    x24 = x21 * (x11 * x5**2 + x12)
    x25 = -x1 * (ax * A[2] + bx * B[2])
    x26 = -x25 - B[2]
    x27 = -x22 - A[1]
    x28 = 22.6274169979695 * x19
    x29 = x13 * x14 * x28
    x30 = x0 * x3
    x31 = x14 * x30
    x32 = x14 * x3
    x33 = x27 * x32
    x34 = x23 * x33 + x31
    x35 = x28 * x34
    x36 = da * db * x1 * x14 * x17 * x18
    x37 = 22.6274169979695 * x15 * x2 * x36 * x8
    x38 = -x25 - A[2]
    x39 = x16 * x36
    x40 = x15 * x30
    x41 = x15 * x3
    x42 = x38 * x41
    x43 = x26 * x42 + x40
    x44 = 22.6274169979695 * x43
    x45 = x20 * x7
    x46 = x45 * (x27**2 * x32 + x31)
    x47 = x38 * x7
    x48 = x39 * x7
    x49 = 13.0639452948436 * x48
    x50 = x49 * (x38**2 * x41 + x40)

    # 18 item(s)
    result[0, 0] = numpy.sum(x21 * (x0 * (x10 * x11 + x9) + x13 * x5))
    result[0, 1] = numpy.sum(x23 * x24)
    result[0, 2] = numpy.sum(x24 * x26)
    result[1, 0] = numpy.sum(x27 * x29)
    result[1, 1] = numpy.sum(x35 * x8)
    result[1, 2] = numpy.sum(x26 * x27 * x37)
    result[2, 0] = numpy.sum(x29 * x38)
    result[2, 1] = numpy.sum(x23 * x37 * x38)
    result[2, 2] = numpy.sum(x39 * x44 * x8)
    result[3, 0] = numpy.sum(x10 * x46)
    result[3, 1] = numpy.sum(x45 * (x0 * (x23 * x32 + x33) + x27 * x34))
    result[3, 2] = numpy.sum(x26 * x46)
    result[4, 0] = numpy.sum(22.6274169979695 * x10 * x15 * x2 * x27 * x36 * x47)
    result[4, 1] = numpy.sum(x35 * x47)
    result[4, 2] = numpy.sum(x27 * x44 * x48)
    result[5, 0] = numpy.sum(x10 * x50)
    result[5, 1] = numpy.sum(x23 * x50)
    result[5, 2] = numpy.sum(x49 * (x0 * (x26 * x41 + x42) + x38 * x43))
    return result


def ovlp3d_22(ax, da, A, bx, db, B):
    """Cartesian 3D (dd) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 6), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = -x2 - A[0]
    x12 = x11 * x8
    x13 = x12 * x3
    x14 = x10 + x9
    x15 = 2.0 * x10 * x3 + x11 * x14
    x16 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x17 = da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**3.5)
    x18 = x16 * x17
    x19 = 15.084944665313 * x18
    x20 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x21 = 0.564189583547756 * x1
    x22 = x20 * x21
    x23 = -x1 * (ax * A[1] + bx * B[1])
    x24 = -x23 - B[1]
    x25 = x10 + x13
    x26 = 26.1278905896872 * x18
    x27 = x22 * x26
    x28 = x27 * (x0 * (x12 + x3 * x8) + x11 * x25)
    x29 = -x1 * (ax * A[2] + bx * B[2])
    x30 = -x29 - B[2]
    x31 = x20 * x7
    x32 = x24**2 * x31
    x33 = x0 * x31
    x34 = x32 + x33
    x35 = x10 + x11**2 * x8
    x36 = 0.318309886183791 * x6
    x37 = x35 * x36
    x38 = x17 * x20
    x39 = x16 * x7
    x40 = x30**2 * x39
    x41 = x0 * x39
    x42 = x40 + x41
    x43 = 15.084944665313 * x42
    x44 = -x23 - A[1]
    x45 = x15 * x27
    x46 = 45.2548339959391 * x25 * x36
    x47 = x31 * x44
    x48 = x24 * x47
    x49 = x33 + x48
    x50 = x18 * x49
    x51 = 45.2548339959391 * x30
    x52 = x18 * x22 * x25
    x53 = 2.0 * x24 * x33 + x34 * x44
    x54 = x21 * x5
    x55 = x26 * x54
    x56 = x53 * x55
    x57 = x50 * x54
    x58 = x17 * x5
    x59 = x22 * x58
    x60 = 26.1278905896872 * x59
    x61 = -x29 - A[2]
    x62 = 45.2548339959391 * x61
    x63 = x39 * x61
    x64 = x30 * x63
    x65 = x41 + x64
    x66 = 45.2548339959391 * x65
    x67 = x59 * x66
    x68 = 2.0 * x30 * x41 + x42 * x61
    x69 = x60 * x68
    x70 = x31 * x44**2 + x33
    x71 = x36 * x70
    x72 = x55 * (x0 * (x24 * x31 + x47) + x44 * x49)
    x73 = x36 * x58
    x74 = x39 * x61**2 + x41
    x75 = 15.084944665313 * x74
    x76 = x60 * (x0 * (x30 * x39 + x63) + x61 * x65)

    # 36 item(s)
    result[0, 0] = numpy.sum(x19 * x22 * (x0 * (3.0 * x10 + 2.0 * x13 + x9) + x11 * x15))
    result[0, 1] = numpy.sum(x24 * x28)
    result[0, 2] = numpy.sum(x28 * x30)
    result[0, 3] = numpy.sum(x19 * x34 * x37)
    result[0, 4] = numpy.sum(x24 * x27 * x30 * x35)
    result[0, 5] = numpy.sum(x37 * x38 * x43)
    result[1, 0] = numpy.sum(x44 * x45)
    result[1, 1] = numpy.sum(x46 * x50)
    result[1, 2] = numpy.sum(x44 * x51 * x52)
    result[1, 3] = numpy.sum(x11 * x56)
    result[1, 4] = numpy.sum(x11 * x51 * x57)
    result[1, 5] = numpy.sum(x11 * x42 * x44 * x60)
    result[2, 0] = numpy.sum(x45 * x61)
    result[2, 1] = numpy.sum(x24 * x52 * x62)
    result[2, 2] = numpy.sum(x38 * x46 * x65)
    result[2, 3] = numpy.sum(x11 * x34 * x55 * x61)
    result[2, 4] = numpy.sum(x11 * x24 * x67)
    result[2, 5] = numpy.sum(x11 * x69)
    result[3, 0] = numpy.sum(x14 * x19 * x71)
    result[3, 1] = numpy.sum(x3 * x72)
    result[3, 2] = numpy.sum(x3 * x30 * x55 * x70)
    result[3, 3] = numpy.sum(x19 * x54 * (x0 * (x32 + 3.0 * x33 + 2.0 * x48) + x44 * x53))
    result[3, 4] = numpy.sum(x30 * x72)
    result[3, 5] = numpy.sum(x43 * x58 * x71)
    result[4, 0] = numpy.sum(x14 * x27 * x44 * x61)
    result[4, 1] = numpy.sum(x3 * x57 * x62)
    result[4, 2] = numpy.sum(x3 * x44 * x67)
    result[4, 3] = numpy.sum(x56 * x61)
    result[4, 4] = numpy.sum(x49 * x66 * x73)
    result[4, 5] = numpy.sum(x44 * x69)
    result[5, 0] = numpy.sum(x14 * x36 * x38 * x75)
    result[5, 1] = numpy.sum(x24 * x3 * x60 * x74)
    result[5, 2] = numpy.sum(x3 * x76)
    result[5, 3] = numpy.sum(x34 * x73 * x75)
    result[5, 4] = numpy.sum(x24 * x76)
    result[5, 5] = numpy.sum(
        15.084944665313 * x59 * (x0 * (x40 + 3.0 * x41 + 2.0 * x64) + x61 * x68)
    )
    return result


def ovlp3d_23(ax, da, A, bx, db, B):
    """Cartesian 3D (df) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 10), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = x10 + x9
    x12 = x11 * x3
    x13 = -x2 - A[0]
    x14 = x11 * x13
    x15 = x10 * x3
    x16 = 3.0 * x10
    x17 = 2.0 * x15
    x18 = x12 + x17
    x19 = x0 * (x16 + 3.0 * x9) + x13 * x18
    x20 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x21 = da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**4.5)
    x22 = x20 * x21
    x23 = 13.4923846833851 * x22
    x24 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x25 = 0.564189583547756 * x1
    x26 = x24 * x25
    x27 = -x1 * (ax * A[1] + bx * B[1])
    x28 = -x27 - B[1]
    x29 = x26 * x28
    x30 = x13 * x8
    x31 = x3 * x30
    x32 = x14 + x17
    x33 = 30.169889330626 * x22
    x34 = x33 * (x0 * (x16 + 2.0 * x31 + x9) + x13 * x32)
    x35 = -x1 * (ax * A[2] + bx * B[2])
    x36 = -x35 - B[2]
    x37 = x26 * x36
    x38 = x10 + x31
    x39 = x0 * (x3 * x8 + x30) + x13 * x38
    x40 = x24 * x7
    x41 = x28**2 * x40
    x42 = x0 * x40
    x43 = x41 + x42
    x44 = 0.318309886183791 * x6
    x45 = x43 * x44
    x46 = x33 * x45
    x47 = 52.2557811793745 * x22 * x36
    x48 = x21 * x44
    x49 = x24 * x48
    x50 = x20 * x7
    x51 = x36**2 * x50
    x52 = x0 * x50
    x53 = x51 + x52
    x54 = 30.169889330626 * x53
    x55 = x49 * x54
    x56 = x28 * x43
    x57 = x28 * x42
    x58 = 2.0 * x57
    x59 = x56 + x58
    x60 = x10 + x13**2 * x8
    x61 = x23 * x44
    x62 = x36 * x53
    x63 = x36 * x52
    x64 = 2.0 * x63
    x65 = x62 + x64
    x66 = 13.4923846833851 * x65
    x67 = -x27 - A[1]
    x68 = 23.3694957868871 * x22
    x69 = x19 * x26 * x68
    x70 = 52.2557811793745 * x32
    x71 = x40 * x67
    x72 = x28 * x71
    x73 = x42 + x72
    x74 = x22 * x44
    x75 = x73 * x74
    x76 = x43 * x67
    x77 = x58 + x76
    x78 = 52.2557811793745 * x38
    x79 = 90.5096679918781 * x38
    x80 = x49 * x78
    x81 = x25 * x5
    x82 = x13 * x81
    x83 = 3.0 * x42
    x84 = x0 * (3.0 * x41 + x83) + x59 * x67
    x85 = x68 * x84
    x86 = x48 * x5
    x87 = 52.2557811793745 * x86
    x88 = x73 * x87
    x89 = x21 * x5
    x90 = x26 * x89
    x91 = -x35 - A[2]
    x92 = x22 * x91
    x93 = x50 * x91
    x94 = x36 * x93
    x95 = x52 + x94
    x96 = x49 * x95
    x97 = x53 * x91
    x98 = x64 + x97
    x99 = x68 * x91
    x100 = x45 * x89
    x101 = 52.2557811793745 * x13
    x102 = x28 * x90
    x103 = 3.0 * x52
    x104 = x0 * (x103 + 3.0 * x51) + x65 * x91
    x105 = 23.3694957868871 * x104 * x90
    x106 = x40 * x67**2 + x42
    x107 = x0 * (x28 * x40 + x71) + x67 * x73
    x108 = x11 * x33 * x44
    x109 = x3 * x81
    x110 = x33 * (x0 * (x41 + 2.0 * x72 + x83) + x67 * x77)
    x111 = x54 * x86
    x112 = 52.2557811793745 * x11
    x113 = x3 * x90
    x114 = 52.2557811793745 * x113
    x115 = x50 * x91**2 + x52
    x116 = 13.4923846833851 * x115
    x117 = x11 * x49
    x118 = 30.169889330626 * x115
    x119 = x0 * (x36 * x50 + x93) + x91 * x95
    x120 = 30.169889330626 * x119
    x121 = 30.169889330626 * x0 * (x103 + x51 + 2.0 * x94) + 30.169889330626 * x91 * x98

    # 60 item(s)
    result[0, 0] = numpy.sum(x23 * x26 * (x0 * (x12 + 3.0 * x14 + 8.0 * x15) + x13 * x19))
    result[0, 1] = numpy.sum(x29 * x34)
    result[0, 2] = numpy.sum(x34 * x37)
    result[0, 3] = numpy.sum(x39 * x46)
    result[0, 4] = numpy.sum(x29 * x39 * x47)
    result[0, 5] = numpy.sum(x39 * x55)
    result[0, 6] = numpy.sum(x59 * x60 * x61)
    result[0, 7] = numpy.sum(x36 * x46 * x60)
    result[0, 8] = numpy.sum(x28 * x55 * x60)
    result[0, 9] = numpy.sum(x49 * x60 * x66)
    result[1, 0] = numpy.sum(x67 * x69)
    result[1, 1] = numpy.sum(x70 * x75)
    result[1, 2] = numpy.sum(x22 * x37 * x67 * x70)
    result[1, 3] = numpy.sum(x74 * x77 * x78)
    result[1, 4] = numpy.sum(x36 * x75 * x79)
    result[1, 5] = numpy.sum(x53 * x67 * x80)
    result[1, 6] = numpy.sum(x82 * x85)
    result[1, 7] = numpy.sum(x47 * x77 * x82)
    result[1, 8] = numpy.sum(x13 * x53 * x88)
    result[1, 9] = numpy.sum(23.3694957868871 * x13 * x65 * x67 * x90)
    result[2, 0] = numpy.sum(x69 * x91)
    result[2, 1] = numpy.sum(x29 * x70 * x92)
    result[2, 2] = numpy.sum(x70 * x96)
    result[2, 3] = numpy.sum(x45 * x78 * x92)
    result[2, 4] = numpy.sum(x28 * x79 * x96)
    result[2, 5] = numpy.sum(x80 * x98)
    result[2, 6] = numpy.sum(x59 * x82 * x99)
    result[2, 7] = numpy.sum(x100 * x101 * x95)
    result[2, 8] = numpy.sum(x101 * x102 * x98)
    result[2, 9] = numpy.sum(x105 * x13)
    result[3, 0] = numpy.sum(x106 * x18 * x61)
    result[3, 1] = numpy.sum(x107 * x108)
    result[3, 2] = numpy.sum(x106 * x108 * x36)
    result[3, 3] = numpy.sum(x109 * x110)
    result[3, 4] = numpy.sum(x107 * x109 * x47)
    result[3, 5] = numpy.sum(x106 * x111 * x3)
    result[3, 6] = numpy.sum(x23 * x81 * (x0 * (x56 + 8.0 * x57 + 3.0 * x76) + x67 * x84))
    result[3, 7] = numpy.sum(x110 * x36 * x81)
    result[3, 8] = numpy.sum(x107 * x111)
    result[3, 9] = numpy.sum(x106 * x66 * x86)
    result[4, 0] = numpy.sum(x18 * x26 * x67 * x99)
    result[4, 1] = numpy.sum(x112 * x75 * x91)
    result[4, 2] = numpy.sum(x112 * x67 * x96)
    result[4, 3] = numpy.sum(52.2557811793745 * x109 * x77 * x92)
    result[4, 4] = numpy.sum(90.5096679918781 * x3 * x73 * x86 * x95)
    result[4, 5] = numpy.sum(x114 * x67 * x98)
    result[4, 6] = numpy.sum(x81 * x85 * x91)
    result[4, 7] = numpy.sum(x77 * x87 * x95)
    result[4, 8] = numpy.sum(x88 * x98)
    result[4, 9] = numpy.sum(x105 * x67)
    result[5, 0] = numpy.sum(x116 * x18 * x49)
    result[5, 1] = numpy.sum(x117 * x118 * x28)
    result[5, 2] = numpy.sum(x117 * x120)
    result[5, 3] = numpy.sum(x100 * x118 * x3)
    result[5, 4] = numpy.sum(x114 * x119 * x28)
    result[5, 5] = numpy.sum(x113 * x121)
    result[5, 6] = numpy.sum(x116 * x59 * x86)
    result[5, 7] = numpy.sum(x100 * x120)
    result[5, 8] = numpy.sum(x102 * x121)
    result[5, 9] = numpy.sum(
        13.4923846833851 * x90 * (x0 * (x62 + 8.0 * x63 + 3.0 * x97) + x104 * x91)
    )
    return result


def ovlp3d_24(ax, da, A, bx, db, B):
    """Cartesian 3D (dg) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 15), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = x10 + x9
    x12 = x11 * x3
    x13 = x10 * x3
    x14 = 2.0 * x13
    x15 = x12 + x14
    x16 = x15 * x3
    x17 = -x2 - A[0]
    x18 = x15 * x17
    x19 = 3.0 * x10
    x20 = x0 * (x19 + 3.0 * x9)
    x21 = 8.0 * x13
    x22 = x16 + x20
    x23 = x0 * (4.0 * x12 + x21) + x17 * x22
    x24 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x25 = da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**5.5)
    x26 = x24 * x25
    x27 = 10.1992841329868 * x26
    x28 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x29 = 0.564189583547756 * x1
    x30 = x28 * x29
    x31 = -x1 * (ax * A[1] + bx * B[1])
    x32 = -x31 - B[1]
    x33 = 26.9847693667702 * x32
    x34 = x11 * x17
    x35 = x18 + x20
    x36 = x26 * x30
    x37 = x36 * (x0 * (x12 + x21 + 3.0 * x34) + x17 * x35)
    x38 = -x1 * (ax * A[2] + bx * B[2])
    x39 = -x38 - B[2]
    x40 = 26.9847693667702 * x39
    x41 = x17 * x8
    x42 = x3 * x41
    x43 = x14 + x34
    x44 = x0 * (x19 + 2.0 * x42 + x9) + x17 * x43
    x45 = x26 * x44
    x46 = x28 * x7
    x47 = x32**2 * x46
    x48 = x0 * x46
    x49 = x47 + x48
    x50 = 34.8371874529163 * x49
    x51 = 0.318309886183791 * x6
    x52 = x50 * x51
    x53 = 60.3397786612521 * x39
    x54 = x30 * x32
    x55 = x25 * x51
    x56 = x28 * x55
    x57 = x24 * x7
    x58 = x39**2 * x57
    x59 = x0 * x57
    x60 = x58 + x59
    x61 = 34.8371874529163 * x60
    x62 = x10 + x42
    x63 = x0 * (x3 * x8 + x41) + x17 * x62
    x64 = 26.9847693667702 * x63
    x65 = x32 * x49
    x66 = x32 * x48
    x67 = 2.0 * x66
    x68 = x65 + x67
    x69 = x26 * x51
    x70 = x68 * x69
    x71 = 60.3397786612521 * x63
    x72 = x39 * x69
    x73 = x32 * x56
    x74 = x39 * x60
    x75 = x39 * x59
    x76 = 2.0 * x75
    x77 = x74 + x76
    x78 = x56 * x77
    x79 = 3.0 * x48
    x80 = x0 * (3.0 * x47 + x79)
    x81 = x32 * x68
    x82 = x80 + x81
    x83 = x10 + x17**2 * x8
    x84 = x27 * x51
    x85 = 0.179587122125167 * x25
    x86 = x60 * x85
    x87 = 3.0 * x59
    x88 = x0 * (3.0 * x58 + x87)
    x89 = x39 * x77
    x90 = x88 + x89
    x91 = 10.1992841329868 * x90
    x92 = -x31 - A[1]
    x93 = 17.6656783191643 * x92
    x94 = x23 * x36
    x95 = x46 * x92
    x96 = x32 * x95
    x97 = x48 + x96
    x98 = 46.7389915737742 * x69
    x99 = x97 * x98
    x100 = 46.7389915737742 * x92
    x101 = x35 * x36
    x102 = 60.3397786612521 * x43
    x103 = x49 * x92
    x104 = x103 + x67
    x105 = x104 * x69
    x106 = 104.511562358749 * x97
    x107 = x102 * x56
    x108 = x68 * x92
    x109 = x108 + x80
    x110 = 104.511562358749 * x62
    x111 = 46.7389915737742 * x62
    x112 = 17.6656783191643 * x17
    x113 = 8.0 * x66
    x114 = x0 * (x113 + 4.0 * x65) + x82 * x92
    x115 = x29 * x5
    x116 = x115 * x26
    x117 = x114 * x116
    x118 = 46.7389915737742 * x17
    x119 = x109 * x116
    x120 = x5 * x55
    x121 = x104 * x120
    x122 = x120 * x77
    x123 = x25 * x5
    x124 = x123 * x30
    x125 = -x38 - A[2]
    x126 = 17.6656783191643 * x125
    x127 = 46.7389915737742 * x125
    x128 = x125 * x57
    x129 = x128 * x39
    x130 = x129 + x59
    x131 = 46.7389915737742 * x56
    x132 = x130 * x131
    x133 = 104.511562358749 * x130
    x134 = x125 * x60
    x135 = x134 + x76
    x136 = x130 * x85
    x137 = x125 * x77
    x138 = x137 + x88
    x139 = 46.7389915737742 * x120
    x140 = x130 * x139
    x141 = 60.3397786612521 * x135
    x142 = x120 * x49
    x143 = x124 * x138
    x144 = 8.0 * x75
    x145 = x0 * (x144 + 4.0 * x74) + x125 * x90
    x146 = x124 * x145
    x147 = x46 * x92**2 + x48
    x148 = x0 * (x32 * x46 + x95) + x92 * x97
    x149 = 26.9847693667702 * x148
    x150 = x15 * x69
    x151 = x0 * (x47 + x79 + 2.0 * x96) + x104 * x92
    x152 = 34.8371874529163 * x11
    x153 = 60.3397786612521 * x148
    x154 = 26.9847693667702 * x3
    x155 = x116 * (x0 * (3.0 * x103 + x113 + x65) + x109 * x92)
    x156 = x120 * x3
    x157 = x125**2 * x57 + x59
    x158 = 10.1992841329868 * x157
    x159 = x15 * x56
    x160 = x0 * (x128 + x39 * x57) + x125 * x130
    x161 = 26.9847693667702 * x160
    x162 = 60.3397786612521 * x160
    x163 = x0 * (2.0 * x129 + x58 + x87) + x125 * x135
    x164 = x120 * x68
    x165 = x123 * x163
    x166 = x124 * (x0 * (3.0 * x134 + x144 + x74) + x125 * x138)

    # 90 item(s)
    result[0, 0] = numpy.sum(x27 * x30 * (x0 * (x16 + 4.0 * x18 + 5.0 * x20) + x17 * x23))
    result[0, 1] = numpy.sum(x33 * x37)
    result[0, 2] = numpy.sum(x37 * x40)
    result[0, 3] = numpy.sum(x45 * x52)
    result[0, 4] = numpy.sum(x45 * x53 * x54)
    result[0, 5] = numpy.sum(x44 * x56 * x61)
    result[0, 6] = numpy.sum(x64 * x70)
    result[0, 7] = numpy.sum(x49 * x71 * x72)
    result[0, 8] = numpy.sum(x60 * x71 * x73)
    result[0, 9] = numpy.sum(x64 * x78)
    result[0, 10] = numpy.sum(x82 * x83 * x84)
    result[0, 11] = numpy.sum(x40 * x70 * x83)
    result[0, 12] = numpy.sum(x50 * x83 * x86)
    result[0, 13] = numpy.sum(x33 * x78 * x83)
    result[0, 14] = numpy.sum(x56 * x83 * x91)
    result[1, 0] = numpy.sum(x93 * x94)
    result[1, 1] = numpy.sum(x35 * x99)
    result[1, 2] = numpy.sum(x100 * x101 * x39)
    result[1, 3] = numpy.sum(x102 * x105)
    result[1, 4] = numpy.sum(x106 * x43 * x72)
    result[1, 5] = numpy.sum(x107 * x60 * x92)
    result[1, 6] = numpy.sum(x109 * x62 * x98)
    result[1, 7] = numpy.sum(x105 * x110 * x39)
    result[1, 8] = numpy.sum(x110 * x86 * x97)
    result[1, 9] = numpy.sum(x111 * x78 * x92)
    result[1, 10] = numpy.sum(x112 * x117)
    result[1, 11] = numpy.sum(x118 * x119 * x39)
    result[1, 12] = numpy.sum(60.3397786612521 * x121 * x17 * x60)
    result[1, 13] = numpy.sum(x118 * x122 * x97)
    result[1, 14] = numpy.sum(x124 * x17 * x90 * x93)
    result[2, 0] = numpy.sum(x126 * x94)
    result[2, 1] = numpy.sum(x101 * x127 * x32)
    result[2, 2] = numpy.sum(x132 * x35)
    result[2, 3] = numpy.sum(x102 * x125 * x49 * x69)
    result[2, 4] = numpy.sum(x133 * x43 * x73)
    result[2, 5] = numpy.sum(x107 * x135)
    result[2, 6] = numpy.sum(x111 * x125 * x70)
    result[2, 7] = numpy.sum(x110 * x136 * x49)
    result[2, 8] = numpy.sum(x110 * x135 * x73)
    result[2, 9] = numpy.sum(x131 * x138 * x62)
    result[2, 10] = numpy.sum(x112 * x116 * x125 * x82)
    result[2, 11] = numpy.sum(x140 * x17 * x68)
    result[2, 12] = numpy.sum(x141 * x142 * x17)
    result[2, 13] = numpy.sum(x118 * x143 * x32)
    result[2, 14] = numpy.sum(x112 * x146)
    result[3, 0] = numpy.sum(x147 * x22 * x84)
    result[3, 1] = numpy.sum(x149 * x150)
    result[3, 2] = numpy.sum(x147 * x150 * x40)
    result[3, 3] = numpy.sum(x151 * x152 * x69)
    result[3, 4] = numpy.sum(x11 * x153 * x72)
    result[3, 5] = numpy.sum(x147 * x152 * x86)
    result[3, 6] = numpy.sum(x154 * x155)
    result[3, 7] = numpy.sum(x116 * x151 * x3 * x53)
    result[3, 8] = numpy.sum(x153 * x156 * x60)
    result[3, 9] = numpy.sum(x122 * x147 * x154)
    result[3, 10] = numpy.sum(
        x115 * x27 * (x0 * (4.0 * x108 + 5.0 * x80 + x81) + x114 * x92)
    )
    result[3, 11] = numpy.sum(x155 * x40)
    result[3, 12] = numpy.sum(x120 * x151 * x61)
    result[3, 13] = numpy.sum(x122 * x149)
    result[3, 14] = numpy.sum(x120 * x147 * x91)
    result[4, 0] = numpy.sum(x125 * x22 * x36 * x93)
    result[4, 1] = numpy.sum(x125 * x15 * x99)
    result[4, 2] = numpy.sum(x132 * x15 * x92)
    result[4, 3] = numpy.sum(60.3397786612521 * x105 * x11 * x125)
    result[4, 4] = numpy.sum(x106 * x11 * x136)
    result[4, 5] = numpy.sum(x11 * x141 * x56 * x92)
    result[4, 6] = numpy.sum(x119 * x127 * x3)
    result[4, 7] = numpy.sum(x121 * x133 * x3)
    result[4, 8] = numpy.sum(x106 * x135 * x156)
    result[4, 9] = numpy.sum(x100 * x143 * x3)
    result[4, 10] = numpy.sum(x117 * x126)
    result[4, 11] = numpy.sum(x109 * x140)
    result[4, 12] = numpy.sum(x121 * x141)
    result[4, 13] = numpy.sum(x138 * x139 * x97)
    result[4, 14] = numpy.sum(x146 * x93)
    result[5, 0] = numpy.sum(x158 * x22 * x56)
    result[5, 1] = numpy.sum(x157 * x159 * x33)
    result[5, 2] = numpy.sum(x159 * x161)
    result[5, 3] = numpy.sum(x11 * x157 * x50 * x85)
    result[5, 4] = numpy.sum(x11 * x162 * x73)
    result[5, 5] = numpy.sum(x152 * x163 * x56)
    result[5, 6] = numpy.sum(x154 * x157 * x164)
    result[5, 7] = numpy.sum(x142 * x162 * x3)
    result[5, 8] = numpy.sum(60.3397786612521 * x165 * x3 * x54)
    result[5, 9] = numpy.sum(x154 * x166)
    result[5, 10] = numpy.sum(x120 * x158 * x82)
    result[5, 11] = numpy.sum(x161 * x164)
    result[5, 12] = numpy.sum(x165 * x52)
    result[5, 13] = numpy.sum(x166 * x33)
    result[5, 14] = numpy.sum(
        10.1992841329868 * x124 * (x0 * (4.0 * x137 + 5.0 * x88 + x89) + x125 * x145)
    )
    return result


def ovlp3d_30(ax, da, A, bx, db, B):
    """Cartesian 3D (fs) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 1), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = ax * bx * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = x0 * (ax * A[0] + bx * B[0]) - A[0]
    x4 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x5 = numpy.sqrt(x0)
    x6 = 1.77245385090552 * x5
    x7 = x4 * x6
    x8 = 0.5 / (ax + bx)
    x9 = x7 * x8
    x10 = x3**2 * x7 + x9
    x11 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x12 = 0.564189583547756
    x13 = numpy.sqrt(ax**4.5)
    x14 = numpy.sqrt(bx**1.5)
    x15 = da * db * x0 * x11 * x12 * x13 * x14
    x16 = 5.84237394672177 * x15
    x17 = x0 * (ax * A[1] + bx * B[1]) - A[1]
    x18 = 13.0639452948436 * x15
    x19 = x10 * x18 * x2
    x20 = x0 * (ax * A[2] + bx * B[2]) - A[2]
    x21 = x2 * x6
    x22 = x21 * x8
    x23 = x17**2 * x21 + x22
    x24 = x18 * x23 * x4
    x25 = da * db * x0 * x13 * x14 * x2 * x4
    x26 = x11 * x6
    x27 = x26 * x8
    x28 = x20**2 * x26 + x27
    x29 = x12 * x25
    x30 = 13.0639452948436 * x28 * x29

    # 10 item(s)
    result[0, 0] = numpy.sum(x16 * x2 * x3 * (x10 + 2.0 * x9))
    result[1, 0] = numpy.sum(x17 * x19)
    result[2, 0] = numpy.sum(x19 * x20)
    result[3, 0] = numpy.sum(x24 * x3)
    result[4, 0] = numpy.sum(22.6274169979695 * x11 * x17 * x20 * x25 * x3 * x5)
    result[5, 0] = numpy.sum(x3 * x30)
    result[6, 0] = numpy.sum(x16 * x17 * x4 * (2.0 * x22 + x23))
    result[7, 0] = numpy.sum(x20 * x24)
    result[8, 0] = numpy.sum(x17 * x30)
    result[9, 0] = numpy.sum(5.84237394672177 * x20 * x29 * (2.0 * x27 + x28))
    return result


def ovlp3d_31(ax, da, A, bx, db, B):
    """Cartesian 3D (fp) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 3), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = -x2 - B[0]
    x12 = x3 * x8
    x13 = x11 * x12
    x14 = x10 + x13
    x15 = x0 * (x11 * x8 + x12) + x14 * x3
    x16 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x17 = 0.564189583547756
    x18 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x19 = da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**2.5)
    x20 = x18 * x19
    x21 = x1 * x17 * x20
    x22 = 11.6847478934435 * x21
    x23 = x16 * x22
    x24 = -x1 * (ax * A[1] + bx * B[1])
    x25 = -x24 - B[1]
    x26 = x10 + x9
    x27 = x23 * (2.0 * x0 * x12 + x26 * x3)
    x28 = -x1 * (ax * A[2] + bx * B[2])
    x29 = -x28 - B[2]
    x30 = -x24 - A[1]
    x31 = x16 * x30
    x32 = 26.1278905896872 * x21
    x33 = x15 * x32
    x34 = x0 * x7
    x35 = x16 * x34
    x36 = x16 * x7
    x37 = x30 * x36
    x38 = x25 * x37
    x39 = x35 + x38
    x40 = 26.1278905896872 * x26
    x41 = 0.318309886183791 * x6
    x42 = x40 * x41
    x43 = x21 * x40
    x44 = -x28 - A[2]
    x45 = x16 * x44
    x46 = x18 * x34
    x47 = x18 * x7
    x48 = x44 * x47
    x49 = x29 * x48
    x50 = x46 + x49
    x51 = x16 * x19
    x52 = x50 * x51
    x53 = x30**2 * x36
    x54 = x35 + x53
    x55 = 26.1278905896872 * x41
    x56 = x14 * x55
    x57 = x3 * x32
    x58 = x0 * (x25 * x36 + x37) + x30 * x39
    x59 = x5 * x58
    x60 = x5 * x54
    x61 = 45.2548339959391 * x21 * x44
    x62 = x39 * x5
    x63 = x1 * x17 * x5
    x64 = x44**2 * x47
    x65 = x46 + x64
    x66 = x51 * x65
    x67 = 26.1278905896872 * x3
    x68 = x63 * x66
    x69 = x0 * (x29 * x47 + x48) + x44 * x50
    x70 = x51 * x63
    x71 = x69 * x70
    x72 = x22 * x5
    x73 = x30 * x72 * (2.0 * x35 + x54)
    x74 = x32 * x44
    x75 = x19 * x55
    x76 = 26.1278905896872 * x30
    x77 = 11.6847478934435 * x70
    x78 = x44 * x77 * (2.0 * x46 + x65)

    # 30 item(s)
    result[0, 0] = numpy.sum(x23 * (x0 * (3.0 * x10 + 2.0 * x13 + x9) + x15 * x3))
    result[0, 1] = numpy.sum(x25 * x27)
    result[0, 2] = numpy.sum(x27 * x29)
    result[1, 0] = numpy.sum(x31 * x33)
    result[1, 1] = numpy.sum(x20 * x39 * x42)
    result[1, 2] = numpy.sum(x29 * x31 * x43)
    result[2, 0] = numpy.sum(x33 * x45)
    result[2, 1] = numpy.sum(x25 * x43 * x45)
    result[2, 2] = numpy.sum(x42 * x52)
    result[3, 0] = numpy.sum(x20 * x54 * x56)
    result[3, 1] = numpy.sum(x57 * x59)
    result[3, 2] = numpy.sum(x29 * x57 * x60)
    result[4, 0] = numpy.sum(x14 * x31 * x61)
    result[4, 1] = numpy.sum(x3 * x61 * x62)
    result[4, 2] = numpy.sum(45.2548339959391 * x3 * x30 * x52 * x63)
    result[5, 0] = numpy.sum(x56 * x66)
    result[5, 1] = numpy.sum(x25 * x67 * x68)
    result[5, 2] = numpy.sum(x67 * x71)
    result[6, 0] = numpy.sum(x11 * x73)
    result[6, 1] = numpy.sum(x72 * (x0 * (3.0 * x35 + 2.0 * x38 + x53) + x30 * x58))
    result[6, 2] = numpy.sum(x29 * x73)
    result[7, 0] = numpy.sum(x11 * x60 * x74)
    result[7, 1] = numpy.sum(x59 * x74)
    result[7, 2] = numpy.sum(x50 * x60 * x75)
    result[8, 0] = numpy.sum(x11 * x68 * x76)
    result[8, 1] = numpy.sum(x62 * x65 * x75)
    result[8, 2] = numpy.sum(x71 * x76)
    result[9, 0] = numpy.sum(x11 * x78)
    result[9, 1] = numpy.sum(x25 * x78)
    result[9, 2] = numpy.sum(x77 * (x0 * (3.0 * x46 + 2.0 * x49 + x64) + x44 * x69))
    return result


def ovlp3d_32(ax, da, A, bx, db, B):
    """Cartesian 3D (fd) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 6), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = numpy.sqrt(x1)
    x3 = 1.77245385090552 * x2
    x4 = -x1 * (ax * A[0] + bx * B[0])
    x5 = -x4 - A[0]
    x6 = ax * bx * x1
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x4 - B[0]
    x11 = x3 * x7
    x12 = x10 * x11
    x13 = x0 * (x12 + x9)
    x14 = x10**2 * x11
    x15 = x0 * x11
    x16 = x14 + x15
    x17 = x16 * x5
    x18 = x10 * x9
    x19 = x15 + x18
    x20 = x19 * x5
    x21 = x0 * x12
    x22 = 3.0 * x15 + 2.0 * x18
    x23 = x17 + 2.0 * x21
    x24 = x0 * (x14 + x22) + x23 * x5
    x25 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x26 = da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**3.5)
    x27 = x25 * x26
    x28 = 13.4923846833851 * x27
    x29 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x30 = 0.564189583547756 * x1
    x31 = x29 * x30
    x32 = -x1 * (ax * A[1] + bx * B[1])
    x33 = -x32 - B[1]
    x34 = 23.3694957868871 * x33
    x35 = x11 * x5**2
    x36 = x13 + x20
    x37 = x27 * x31
    x38 = x37 * (x0 * (x22 + x35) + x36 * x5)
    x39 = -x1 * (ax * A[2] + bx * B[2])
    x40 = -x39 - B[2]
    x41 = 23.3694957868871 * x40
    x42 = x29 * x3
    x43 = x33**2 * x42
    x44 = x0 * x42
    x45 = x43 + x44
    x46 = x15 + x35
    x47 = 2.0 * x0 * x9 + x46 * x5
    x48 = 0.318309886183791 * x2
    x49 = x47 * x48
    x50 = x37 * x40
    x51 = x26 * x29
    x52 = x25 * x3
    x53 = x40**2 * x52
    x54 = x0 * x52
    x55 = x53 + x54
    x56 = 13.4923846833851 * x55
    x57 = -x32 - A[1]
    x58 = 30.169889330626 * x57
    x59 = x24 * x37
    x60 = 52.2557811793745 * x36
    x61 = x42 * x57
    x62 = x33 * x61
    x63 = x44 + x62
    x64 = x27 * x48
    x65 = x63 * x64
    x66 = x45 * x57
    x67 = 2.0 * x44
    x68 = x33 * x67 + x66
    x69 = 30.169889330626 * x46
    x70 = x64 * x69
    x71 = 52.2557811793745 * x46
    x72 = x48 * x51
    x73 = x69 * x72
    x74 = -x39 - A[2]
    x75 = 30.169889330626 * x74
    x76 = x37 * x74
    x77 = x52 * x74
    x78 = x40 * x77
    x79 = x54 + x78
    x80 = x72 * x79
    x81 = x55 * x74
    x82 = 2.0 * x54
    x83 = x40 * x82 + x81
    x84 = 30.169889330626 * x23
    x85 = x42 * x57**2
    x86 = x44 + x85
    x87 = x64 * x86
    x88 = x0 * (x33 * x42 + x61)
    x89 = x57 * x63
    x90 = x88 + x89
    x91 = 52.2557811793745 * x19
    x92 = 30.169889330626 * x8
    x93 = 3.0 * x44 + 2.0 * x62
    x94 = x0 * (x43 + x93) + x57 * x68
    x95 = x27 * x94
    x96 = 52.2557811793745 * x90
    x97 = x27 * x96
    x98 = x30 * x8
    x99 = x26 * x48
    x100 = x92 * x99
    x101 = 52.2557811793745 * x57
    x102 = 90.5096679918781 * x19
    x103 = x26 * x31
    x104 = x103 * x8
    x105 = x52 * x74**2
    x106 = x105 + x54
    x107 = x106 * x72
    x108 = x0 * (x40 * x52 + x77)
    x109 = x74 * x79
    x110 = x108 + x109
    x111 = 52.2557811793745 * x110
    x112 = 3.0 * x54 + 2.0 * x78
    x113 = x0 * (x112 + x53) + x74 * x83
    x114 = x57 * (x67 + x86)
    x115 = x114 * x48
    x116 = 23.3694957868871 * x10
    x117 = x30 * x7
    x118 = x117 * x27 * (x0 * (x85 + x93) + x57 * x90)
    x119 = x10 * x117
    x120 = x26 * x7
    x121 = x120 * x48
    x122 = x121 * x79
    x123 = 52.2557811793745 * x10
    x124 = 30.169889330626 * x121
    x125 = x121 * x63
    x126 = x120 * x31
    x127 = x10 * x126
    x128 = x74 * (x106 + x82)
    x129 = 13.4923846833851 * x128
    x130 = x126 * (x0 * (x105 + x112) + x110 * x74)

    # 60 item(s)
    result[0, 0] = numpy.sum(
        x28 * x31 * (2.0 * x0 * (x13 + x17 + x20 + 2.0 * x21) + x24 * x5)
    )
    result[0, 1] = numpy.sum(x34 * x38)
    result[0, 2] = numpy.sum(x38 * x41)
    result[0, 3] = numpy.sum(x28 * x45 * x49)
    result[0, 4] = numpy.sum(x34 * x47 * x50)
    result[0, 5] = numpy.sum(x49 * x51 * x56)
    result[1, 0] = numpy.sum(x58 * x59)
    result[1, 1] = numpy.sum(x60 * x65)
    result[1, 2] = numpy.sum(x50 * x57 * x60)
    result[1, 3] = numpy.sum(x68 * x70)
    result[1, 4] = numpy.sum(x40 * x65 * x71)
    result[1, 5] = numpy.sum(x55 * x57 * x73)
    result[2, 0] = numpy.sum(x59 * x75)
    result[2, 1] = numpy.sum(x33 * x60 * x76)
    result[2, 2] = numpy.sum(x60 * x80)
    result[2, 3] = numpy.sum(x45 * x70 * x74)
    result[2, 4] = numpy.sum(x33 * x71 * x80)
    result[2, 5] = numpy.sum(x73 * x83)
    result[3, 0] = numpy.sum(x84 * x87)
    result[3, 1] = numpy.sum(x64 * x90 * x91)
    result[3, 2] = numpy.sum(x40 * x87 * x91)
    result[3, 3] = numpy.sum(x30 * x92 * x95)
    result[3, 4] = numpy.sum(x40 * x97 * x98)
    result[3, 5] = numpy.sum(x100 * x55 * x86)
    result[4, 0] = numpy.sum(x101 * x23 * x76)
    result[4, 1] = numpy.sum(x102 * x65 * x74)
    result[4, 2] = numpy.sum(x102 * x57 * x80)
    result[4, 3] = numpy.sum(52.2557811793745 * x27 * x68 * x74 * x98)
    result[4, 4] = numpy.sum(90.5096679918781 * x63 * x79 * x8 * x99)
    result[4, 5] = numpy.sum(x101 * x104 * x83)
    result[5, 0] = numpy.sum(x107 * x84)
    result[5, 1] = numpy.sum(x107 * x33 * x91)
    result[5, 2] = numpy.sum(x110 * x72 * x91)
    result[5, 3] = numpy.sum(x100 * x106 * x45)
    result[5, 4] = numpy.sum(x104 * x111 * x33)
    result[5, 5] = numpy.sum(x103 * x113 * x92)
    result[6, 0] = numpy.sum(x115 * x16 * x28)
    result[6, 1] = numpy.sum(x116 * x118)
    result[6, 2] = numpy.sum(x114 * x119 * x27 * x41)
    result[6, 3] = numpy.sum(
        x117 * x28 * (2.0 * x0 * (2.0 * x33 * x44 + x66 + x88 + x89) + x57 * x94)
    )
    result[6, 4] = numpy.sum(x118 * x41)
    result[6, 5] = numpy.sum(x115 * x120 * x56)
    result[7, 0] = numpy.sum(x16 * x75 * x87)
    result[7, 1] = numpy.sum(x119 * x74 * x97)
    result[7, 2] = numpy.sum(x122 * x123 * x86)
    result[7, 3] = numpy.sum(x117 * x75 * x95)
    result[7, 4] = numpy.sum(x122 * x96)
    result[7, 5] = numpy.sum(x124 * x83 * x86)
    result[8, 0] = numpy.sum(x107 * x16 * x58)
    result[8, 1] = numpy.sum(x106 * x123 * x125)
    result[8, 2] = numpy.sum(x111 * x127 * x57)
    result[8, 3] = numpy.sum(x106 * x124 * x68)
    result[8, 4] = numpy.sum(x111 * x125)
    result[8, 5] = numpy.sum(x113 * x126 * x58)
    result[9, 0] = numpy.sum(x129 * x16 * x72)
    result[9, 1] = numpy.sum(x127 * x128 * x34)
    result[9, 2] = numpy.sum(x116 * x130)
    result[9, 3] = numpy.sum(x121 * x129 * x45)
    result[9, 4] = numpy.sum(x130 * x34)
    result[9, 5] = numpy.sum(
        13.4923846833851
        * x126
        * (2.0 * x0 * (x108 + x109 + 2.0 * x40 * x54 + x81) + x113 * x74)
    )
    return result


def ovlp3d_33(ax, da, A, bx, db, B):
    """Cartesian 3D (ff) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 10), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = 3.0 * x10
    x12 = x0 * (x11 + 3.0 * x9)
    x13 = -x2 - A[0]
    x14 = x10 + x9
    x15 = x14 * x3
    x16 = x10 * x3
    x17 = 2.0 * x16
    x18 = x15 + x17
    x19 = x13 * x18
    x20 = x13 * x8
    x21 = x20 * x3
    x22 = x11 + 2.0 * x21
    x23 = x0 * (x22 + x9)
    x24 = x13 * x14
    x25 = x17 + x24
    x26 = x13 * x25
    x27 = x12 + x19
    x28 = x0 * (x15 + 8.0 * x16 + 3.0 * x24) + x13 * x27
    x29 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x30 = da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**4.5)
    x31 = x29 * x30
    x32 = 12.0679557322504 * x31
    x33 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x34 = 0.564189583547756 * x1
    x35 = x33 * x34
    x36 = -x1 * (ax * A[1] + bx * B[1])
    x37 = -x36 - B[1]
    x38 = x0 * (x20 + x3 * x8)
    x39 = x10 + x21
    x40 = x13 * x39
    x41 = x23 + x26
    x42 = 26.9847693667702 * x31
    x43 = x35 * x42
    x44 = x43 * (2.0 * x0 * (2.0 * x16 + x24 + x38 + x40) + x13 * x41)
    x45 = -x1 * (ax * A[2] + bx * B[2])
    x46 = -x45 - B[2]
    x47 = x13**2 * x8
    x48 = x38 + x40
    x49 = x0 * (x22 + x47) + x13 * x48
    x50 = x33 * x7
    x51 = x37**2 * x50
    x52 = x0 * x50
    x53 = x51 + x52
    x54 = 0.318309886183791 * x6
    x55 = x53 * x54
    x56 = x42 * x55
    x57 = x35 * x37
    x58 = x31 * x46
    x59 = 46.7389915737742 * x58
    x60 = x29 * x7
    x61 = x46**2 * x60
    x62 = x0 * x60
    x63 = x61 + x62
    x64 = x30 * x63
    x65 = x33 * x54
    x66 = 26.9847693667702 * x65
    x67 = x64 * x66
    x68 = x37 * x53
    x69 = x37 * x52
    x70 = 2.0 * x69
    x71 = x68 + x70
    x72 = x10 + x47
    x73 = x13 * (2.0 * x10 + x72)
    x74 = x32 * x54
    x75 = x46 * x63
    x76 = x46 * x62
    x77 = 2.0 * x76
    x78 = x75 + x77
    x79 = 12.0679557322504 * x30
    x80 = x78 * x79
    x81 = -x36 - A[1]
    x82 = x28 * x43
    x83 = x50 * x81
    x84 = x37 * x83
    x85 = x52 + x84
    x86 = 60.3397786612521 * x30
    x87 = x54 * x86
    x88 = x29 * x87
    x89 = x29 * x86
    x90 = x46 * x89
    x91 = x35 * x81
    x92 = x53 * x81
    x93 = x70 + x92
    x94 = x88 * x93
    x95 = 104.511562358749 * x48
    x96 = x54 * x58
    x97 = 60.3397786612521 * x64
    x98 = x65 * x81
    x99 = 3.0 * x52
    x100 = x0 * (3.0 * x51 + x99)
    x101 = x71 * x81
    x102 = x100 + x101
    x103 = x42 * x54
    x104 = x103 * x72
    x105 = 0.179587122125167
    x106 = x105 * x85
    x107 = x30 * x66
    x108 = x107 * x72
    x109 = -x45 - A[2]
    x110 = x109 * x89
    x111 = x109 * x60
    x112 = x111 * x46
    x113 = x112 + x62
    x114 = x113 * x30
    x115 = 60.3397786612521 * x114
    x116 = x37 * x65
    x117 = x109 * x63
    x118 = x117 + x77
    x119 = x65 * x86
    x120 = x118 * x119
    x121 = x105 * x53
    x122 = 3.0 * x62
    x123 = x0 * (x122 + 3.0 * x61)
    x124 = x109 * x78
    x125 = x123 + x124
    x126 = x50 * x81**2
    x127 = x126 + x52
    x128 = x103 * x127
    x129 = x0 * (x37 * x50 + x83)
    x130 = x81 * x85
    x131 = x129 + x130
    x132 = x25 * x88
    x133 = 2.0 * x84 + x99
    x134 = x0 * (x133 + x51)
    x135 = x81 * x93
    x136 = x134 + x135
    x137 = 104.511562358749 * x39
    x138 = x105 * x127
    x139 = x0 * (x68 + 8.0 * x69 + 3.0 * x92) + x102 * x81
    x140 = x34 * x5
    x141 = x140 * x42
    x142 = x139 * x141
    x143 = x13 * x140
    x144 = x5 * x54
    x145 = x13 * x144
    x146 = 26.9847693667702 * x144
    x147 = x146 * x30
    x148 = x127 * x147
    x149 = x109 * x31
    x150 = 46.7389915737742 * x149
    x151 = 104.511562358749 * x25
    x152 = x149 * x54
    x153 = x118 * x30
    x154 = 104.511562358749 * x145
    x155 = x35 * x5
    x156 = x155 * x81
    x157 = 46.7389915737742 * x30
    x158 = x109**2 * x60
    x159 = x158 + x62
    x160 = x107 * x159
    x161 = x159 * x86
    x162 = x0 * (x111 + x46 * x60)
    x163 = x109 * x113
    x164 = x162 + x163
    x165 = x119 * x164
    x166 = x164 * x30
    x167 = 2.0 * x112 + x122
    x168 = x0 * (x167 + x61)
    x169 = x109 * x118
    x170 = x168 + x169
    x171 = x147 * x159
    x172 = x5 * x55
    x173 = x13 * x86
    x174 = x155 * x37
    x175 = x0 * (3.0 * x117 + x75 + 8.0 * x76) + x109 * x125
    x176 = 26.9847693667702 * x30
    x177 = x155 * x176
    x178 = x175 * x177
    x179 = x81 * (x127 + 2.0 * x52)
    x180 = x0 * (x126 + x133) + x131 * x81
    x181 = x103 * x14
    x182 = x141 * (2.0 * x0 * (x129 + x130 + 2.0 * x69 + x92) + x136 * x81)
    x183 = x140 * x3
    x184 = x146 * x64
    x185 = x144 * x3
    x186 = 104.511562358749 * x185
    x187 = x5 * x87
    x188 = x118 * x187
    x189 = x109 * (x159 + 2.0 * x62)
    x190 = x189 * x79
    x191 = x107 * x14
    x192 = x0 * (x158 + x167) + x109 * x164
    x193 = x172 * x176
    x194 = x177 * (2.0 * x0 * (x117 + x162 + x163 + 2.0 * x76) + x109 * x170)

    # 100 item(s)
    result[0, 0] = numpy.sum(
        x32 * x35 * (x0 * (2.0 * x12 + 2.0 * x19 + 3.0 * x23 + 3.0 * x26) + x13 * x28)
    )
    result[0, 1] = numpy.sum(x37 * x44)
    result[0, 2] = numpy.sum(x44 * x46)
    result[0, 3] = numpy.sum(x49 * x56)
    result[0, 4] = numpy.sum(x49 * x57 * x59)
    result[0, 5] = numpy.sum(x49 * x67)
    result[0, 6] = numpy.sum(x71 * x73 * x74)
    result[0, 7] = numpy.sum(x46 * x56 * x73)
    result[0, 8] = numpy.sum(x37 * x67 * x73)
    result[0, 9] = numpy.sum(x65 * x73 * x80)
    result[1, 0] = numpy.sum(x81 * x82)
    result[1, 1] = numpy.sum(x41 * x85 * x88)
    result[1, 2] = numpy.sum(x41 * x90 * x91)
    result[1, 3] = numpy.sum(x48 * x94)
    result[1, 4] = numpy.sum(x85 * x95 * x96)
    result[1, 5] = numpy.sum(x48 * x97 * x98)
    result[1, 6] = numpy.sum(x102 * x104)
    result[1, 7] = numpy.sum(x46 * x72 * x94)
    result[1, 8] = numpy.sum(x106 * x72 * x97)
    result[1, 9] = numpy.sum(x108 * x78 * x81)
    result[2, 0] = numpy.sum(x109 * x82)
    result[2, 1] = numpy.sum(x110 * x41 * x57)
    result[2, 2] = numpy.sum(x115 * x41 * x65)
    result[2, 3] = numpy.sum(x110 * x48 * x55)
    result[2, 4] = numpy.sum(x114 * x116 * x95)
    result[2, 5] = numpy.sum(x120 * x48)
    result[2, 6] = numpy.sum(x104 * x109 * x71)
    result[2, 7] = numpy.sum(x115 * x121 * x72)
    result[2, 8] = numpy.sum(x120 * x37 * x72)
    result[2, 9] = numpy.sum(x108 * x125)
    result[3, 0] = numpy.sum(x128 * x27)
    result[3, 1] = numpy.sum(x131 * x132)
    result[3, 2] = numpy.sum(x127 * x132 * x46)
    result[3, 3] = numpy.sum(x136 * x39 * x88)
    result[3, 4] = numpy.sum(x131 * x137 * x96)
    result[3, 5] = numpy.sum(x138 * x39 * x97)
    result[3, 6] = numpy.sum(x13 * x142)
    result[3, 7] = numpy.sum(x136 * x143 * x90)
    result[3, 8] = numpy.sum(x131 * x145 * x97)
    result[3, 9] = numpy.sum(x13 * x148 * x78)
    result[4, 0] = numpy.sum(x150 * x27 * x91)
    result[4, 1] = numpy.sum(x151 * x152 * x85)
    result[4, 2] = numpy.sum(x114 * x151 * x98)
    result[4, 3] = numpy.sum(x137 * x152 * x93)
    result[4, 4] = numpy.sum(181.019335983756 * x106 * x114 * x39)
    result[4, 5] = numpy.sum(x137 * x153 * x98)
    result[4, 6] = numpy.sum(x102 * x143 * x150)
    result[4, 7] = numpy.sum(x114 * x154 * x93)
    result[4, 8] = numpy.sum(x153 * x154 * x85)
    result[4, 9] = numpy.sum(x125 * x13 * x156 * x157)
    result[5, 0] = numpy.sum(x160 * x27)
    result[5, 1] = numpy.sum(x116 * x161 * x25)
    result[5, 2] = numpy.sum(x165 * x25)
    result[5, 3] = numpy.sum(x121 * x161 * x39)
    result[5, 4] = numpy.sum(x116 * x137 * x166)
    result[5, 5] = numpy.sum(x119 * x170 * x39)
    result[5, 6] = numpy.sum(x13 * x171 * x71)
    result[5, 7] = numpy.sum(x164 * x172 * x173)
    result[5, 8] = numpy.sum(x170 * x173 * x174)
    result[5, 9] = numpy.sum(x13 * x178)
    result[6, 0] = numpy.sum(x179 * x18 * x74)
    result[6, 1] = numpy.sum(x180 * x181)
    result[6, 2] = numpy.sum(x179 * x181 * x46)
    result[6, 3] = numpy.sum(x182 * x3)
    result[6, 4] = numpy.sum(x180 * x183 * x59)
    result[6, 5] = numpy.sum(x179 * x184 * x3)
    result[6, 6] = numpy.sum(
        x140
        * x32
        * (x0 * (2.0 * x100 + 2.0 * x101 + 3.0 * x134 + 3.0 * x135) + x139 * x81)
    )
    result[6, 7] = numpy.sum(x182 * x46)
    result[6, 8] = numpy.sum(x180 * x184)
    result[6, 9] = numpy.sum(x144 * x179 * x80)
    result[7, 0] = numpy.sum(x109 * x128 * x18)
    result[7, 1] = numpy.sum(x109 * x131 * x14 * x88)
    result[7, 2] = numpy.sum(x115 * x138 * x14)
    result[7, 3] = numpy.sum(x110 * x136 * x183)
    result[7, 4] = numpy.sum(x114 * x131 * x186)
    result[7, 5] = numpy.sum(x127 * x188 * x3)
    result[7, 6] = numpy.sum(x109 * x142)
    result[7, 7] = numpy.sum(x115 * x136 * x144)
    result[7, 8] = numpy.sum(x131 * x188)
    result[7, 9] = numpy.sum(x125 * x148)
    result[8, 0] = numpy.sum(x160 * x18 * x81)
    result[8, 1] = numpy.sum(x106 * x14 * x161)
    result[8, 2] = numpy.sum(x14 * x165 * x81)
    result[8, 3] = numpy.sum(x161 * x185 * x93)
    result[8, 4] = numpy.sum(x166 * x186 * x85)
    result[8, 5] = numpy.sum(x156 * x170 * x3 * x86)
    result[8, 6] = numpy.sum(x102 * x171)
    result[8, 7] = numpy.sum(x164 * x187 * x93)
    result[8, 8] = numpy.sum(x170 * x187 * x85)
    result[8, 9] = numpy.sum(x178 * x81)
    result[9, 0] = numpy.sum(x18 * x190 * x65)
    result[9, 1] = numpy.sum(x189 * x191 * x37)
    result[9, 2] = numpy.sum(x191 * x192)
    result[9, 3] = numpy.sum(x189 * x193 * x3)
    result[9, 4] = numpy.sum(x157 * x174 * x192 * x3)
    result[9, 5] = numpy.sum(x194 * x3)
    result[9, 6] = numpy.sum(x144 * x190 * x71)
    result[9, 7] = numpy.sum(x192 * x193)
    result[9, 8] = numpy.sum(x194 * x37)
    result[9, 9] = numpy.sum(
        x155
        * x79
        * (x0 * (2.0 * x123 + 2.0 * x124 + 3.0 * x168 + 3.0 * x169) + x109 * x175)
    )
    return result


def ovlp3d_34(ax, da, A, bx, db, B):
    """Cartesian 3D (fg) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 15), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = x10 + x9
    x12 = x11 * x3
    x13 = x10 * x3
    x14 = 8.0 * x13
    x15 = x0 * (4.0 * x12 + x14)
    x16 = -x2 - A[0]
    x17 = 3.0 * x10
    x18 = x0 * (x17 + 3.0 * x9)
    x19 = 2.0 * x13
    x20 = x12 + x19
    x21 = x20 * x3
    x22 = x18 + x21
    x23 = x16 * x22
    x24 = x11 * x16
    x25 = x0 * (x12 + x14 + 3.0 * x24)
    x26 = x16 * x20
    x27 = x18 + x26
    x28 = x16 * x27
    x29 = x15 + x23
    x30 = x0 * (5.0 * x18 + x21 + 4.0 * x26) + x16 * x29
    x31 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x32 = da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**5.5)
    x33 = x31 * x32
    x34 = 9.12251705727742 * x33
    x35 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x36 = 0.564189583547756 * x1
    x37 = x35 * x36
    x38 = -x1 * (ax * A[1] + bx * B[1])
    x39 = -x38 - B[1]
    x40 = 24.1359114645008 * x39
    x41 = x16 * x8
    x42 = x3 * x41
    x43 = x17 + 2.0 * x42
    x44 = x0 * (x43 + x9)
    x45 = x19 + x24
    x46 = x16 * x45
    x47 = x25 + x28
    x48 = x33 * x37
    x49 = x48 * (x0 * (2.0 * x18 + 2.0 * x26 + 3.0 * x44 + 3.0 * x46) + x16 * x47)
    x50 = -x1 * (ax * A[2] + bx * B[2])
    x51 = -x50 - B[2]
    x52 = 24.1359114645008 * x51
    x53 = x0 * (x3 * x8 + x41)
    x54 = x10 + x42
    x55 = x16 * x54
    x56 = x44 + x46
    x57 = 2.0 * x0 * (2.0 * x13 + x24 + x53 + x55) + x16 * x56
    x58 = x33 * x57
    x59 = x35 * x7
    x60 = x39**2 * x59
    x61 = x0 * x59
    x62 = x60 + x61
    x63 = 31.1593277158494 * x62
    x64 = 0.318309886183791 * x6
    x65 = x63 * x64
    x66 = 53.9695387335403 * x51
    x67 = x37 * x39
    x68 = x32 * x64
    x69 = x35 * x68
    x70 = x31 * x7
    x71 = x51**2 * x70
    x72 = x0 * x70
    x73 = x71 + x72
    x74 = 31.1593277158494 * x73
    x75 = x16**2 * x8
    x76 = x53 + x55
    x77 = x0 * (x43 + x75) + x16 * x76
    x78 = 24.1359114645008 * x77
    x79 = x39 * x62
    x80 = x39 * x61
    x81 = 2.0 * x80
    x82 = x79 + x81
    x83 = x33 * x64
    x84 = x82 * x83
    x85 = 53.9695387335403 * x83
    x86 = x51 * x85
    x87 = 53.9695387335403 * x69
    x88 = x39 * x87
    x89 = x51 * x73
    x90 = x51 * x72
    x91 = 2.0 * x90
    x92 = x89 + x91
    x93 = x69 * x92
    x94 = 3.0 * x61
    x95 = x0 * (3.0 * x60 + x94)
    x96 = x39 * x82
    x97 = x95 + x96
    x98 = x10 + x75
    x99 = x16 * (2.0 * x10 + x98)
    x100 = x34 * x64
    x101 = 0.179587122125167 * x32
    x102 = x101 * x73
    x103 = 3.0 * x72
    x104 = x0 * (x103 + 3.0 * x71)
    x105 = x51 * x92
    x106 = x104 + x105
    x107 = 9.12251705727742 * x106
    x108 = -x38 - A[1]
    x109 = 20.3985682659737 * x108
    x110 = x30 * x48
    x111 = x108 * x59
    x112 = x111 * x39
    x113 = x112 + x61
    x114 = x47 * x48
    x115 = x108 * x62
    x116 = x115 + x81
    x117 = 69.6743749058326 * x116
    x118 = x56 * x83
    x119 = 120.679557322504 * x51
    x120 = x56 * x69
    x121 = 69.6743749058326 * x73
    x122 = x108 * x82
    x123 = x122 + x95
    x124 = x123 * x85
    x125 = 120.679557322504 * x76
    x126 = x116 * x83
    x127 = 53.9695387335403 * x92
    x128 = x108 * x69
    x129 = 8.0 * x80
    x130 = x0 * (x129 + 4.0 * x79)
    x131 = x108 * x97
    x132 = x130 + x131
    x133 = 20.3985682659737 * x98
    x134 = x133 * x83
    x135 = 69.6743749058326 * x98
    x136 = x101 * x113
    x137 = x133 * x69
    x138 = -x50 - A[2]
    x139 = 20.3985682659737 * x138
    x140 = 53.9695387335403 * x39
    x141 = x138 * x70
    x142 = x141 * x51
    x143 = x142 + x72
    x144 = 69.6743749058326 * x62
    x145 = 120.679557322504 * x39
    x146 = x138 * x73
    x147 = x146 + x91
    x148 = 69.6743749058326 * x147
    x149 = 53.9695387335403 * x82
    x150 = x138 * x83
    x151 = x101 * x143
    x152 = x138 * x92
    x153 = x104 + x152
    x154 = x153 * x87
    x155 = x101 * x62
    x156 = 8.0 * x90
    x157 = x0 * (x156 + 4.0 * x89)
    x158 = x106 * x138
    x159 = x157 + x158
    x160 = 20.3985682659737 * x29
    x161 = x108**2 * x59
    x162 = x161 + x61
    x163 = x162 * x83
    x164 = x0 * (x111 + x39 * x59)
    x165 = x108 * x113
    x166 = x164 + x165
    x167 = x27 * x85
    x168 = 69.6743749058326 * x45
    x169 = 2.0 * x112 + x94
    x170 = x0 * (x169 + x60)
    x171 = x108 * x116
    x172 = x170 + x171
    x173 = x172 * x83
    x174 = x0 * (3.0 * x115 + x129 + x79)
    x175 = x108 * x123
    x176 = x174 + x175
    x177 = 120.679557322504 * x54
    x178 = x101 * x54
    x179 = 20.3985682659737 * x16
    x180 = x0 * (4.0 * x122 + 5.0 * x95 + x96) + x108 * x132
    x181 = x36 * x5
    x182 = x181 * x33
    x183 = x180 * x182
    x184 = x16 * x182
    x185 = x5 * x68
    x186 = x172 * x185
    x187 = x16 * x185
    x188 = 20.3985682659737 * x185
    x189 = x162 * x188
    x190 = 35.3313566383285 * x138
    x191 = 93.4779831475484 * x27
    x192 = 120.679557322504 * x45
    x193 = 209.023124717498 * x136
    x194 = 93.4779831475484 * x54
    x195 = 93.4779831475484 * x187
    x196 = x32 * x5
    x197 = x196 * x37
    x198 = x16 * x197
    x199 = x138**2 * x70
    x200 = x199 + x72
    x201 = x200 * x69
    x202 = x27 * x87
    x203 = x0 * (x141 + x51 * x70)
    x204 = x138 * x143
    x205 = x203 + x204
    x206 = x103 + 2.0 * x142
    x207 = x0 * (x206 + x71)
    x208 = x138 * x147
    x209 = x207 + x208
    x210 = x209 * x69
    x211 = x0 * (3.0 * x146 + x156 + x89)
    x212 = x138 * x153
    x213 = x211 + x212
    x214 = x188 * x200
    x215 = x185 * x209
    x216 = x0 * (5.0 * x104 + x105 + 4.0 * x152) + x138 * x159
    x217 = x197 * x216
    x218 = x108 * (x162 + 2.0 * x61)
    x219 = x0 * (x161 + x169) + x108 * x166
    x220 = 24.1359114645008 * x219
    x221 = x20 * x83
    x222 = 2.0 * x0 * (x115 + x164 + x165 + 2.0 * x80) + x108 * x172
    x223 = 31.1593277158494 * x11
    x224 = 24.1359114645008 * x3
    x225 = x182 * (x0 * (2.0 * x122 + 3.0 * x170 + 3.0 * x171 + 2.0 * x95) + x108 * x176)
    x226 = x182 * x3
    x227 = 53.9695387335403 * x185
    x228 = x227 * x3
    x229 = x185 * x92
    x230 = 53.9695387335403 * x20
    x231 = 69.6743749058326 * x11
    x232 = 120.679557322504 * x11
    x233 = x101 * x11
    x234 = 120.679557322504 * x3
    x235 = x185 * x234
    x236 = x153 * x227
    x237 = x123 * x227
    x238 = 53.9695387335403 * x3
    x239 = x138 * (x200 + 2.0 * x72)
    x240 = 9.12251705727742 * x239
    x241 = x20 * x69
    x242 = x0 * (x199 + x206) + x138 * x205
    x243 = 24.1359114645008 * x242
    x244 = 2.0 * x0 * (x146 + x203 + x204 + 2.0 * x90) + x138 * x209
    x245 = x185 * x82
    x246 = x196 * x244
    x247 = x197 * (x0 * (2.0 * x104 + 2.0 * x152 + 3.0 * x207 + 3.0 * x208) + x138 * x213)

    # 150 item(s)
    result[0, 0] = numpy.sum(
        x34 * x37 * (2.0 * x0 * (x15 + x23 + 2.0 * x25 + 2.0 * x28) + x16 * x30)
    )
    result[0, 1] = numpy.sum(x40 * x49)
    result[0, 2] = numpy.sum(x49 * x52)
    result[0, 3] = numpy.sum(x58 * x65)
    result[0, 4] = numpy.sum(x58 * x66 * x67)
    result[0, 5] = numpy.sum(x57 * x69 * x74)
    result[0, 6] = numpy.sum(x78 * x84)
    result[0, 7] = numpy.sum(x62 * x77 * x86)
    result[0, 8] = numpy.sum(x73 * x77 * x88)
    result[0, 9] = numpy.sum(x78 * x93)
    result[0, 10] = numpy.sum(x100 * x97 * x99)
    result[0, 11] = numpy.sum(x52 * x84 * x99)
    result[0, 12] = numpy.sum(x102 * x63 * x99)
    result[0, 13] = numpy.sum(x40 * x93 * x99)
    result[0, 14] = numpy.sum(x107 * x69 * x99)
    result[1, 0] = numpy.sum(x109 * x110)
    result[1, 1] = numpy.sum(x113 * x47 * x85)
    result[1, 2] = numpy.sum(x108 * x114 * x66)
    result[1, 3] = numpy.sum(x117 * x118)
    result[1, 4] = numpy.sum(x113 * x118 * x119)
    result[1, 5] = numpy.sum(x108 * x120 * x121)
    result[1, 6] = numpy.sum(x124 * x76)
    result[1, 7] = numpy.sum(x125 * x126 * x51)
    result[1, 8] = numpy.sum(x102 * x113 * x125)
    result[1, 9] = numpy.sum(x127 * x128 * x76)
    result[1, 10] = numpy.sum(x132 * x134)
    result[1, 11] = numpy.sum(x124 * x51 * x98)
    result[1, 12] = numpy.sum(x102 * x116 * x135)
    result[1, 13] = numpy.sum(x127 * x136 * x98)
    result[1, 14] = numpy.sum(x106 * x108 * x137)
    result[2, 0] = numpy.sum(x110 * x139)
    result[2, 1] = numpy.sum(x114 * x138 * x140)
    result[2, 2] = numpy.sum(x143 * x47 * x87)
    result[2, 3] = numpy.sum(x118 * x138 * x144)
    result[2, 4] = numpy.sum(x120 * x143 * x145)
    result[2, 5] = numpy.sum(x120 * x148)
    result[2, 6] = numpy.sum(x149 * x150 * x76)
    result[2, 7] = numpy.sum(x125 * x151 * x62)
    result[2, 8] = numpy.sum(x125 * x147 * x39 * x69)
    result[2, 9] = numpy.sum(x154 * x76)
    result[2, 10] = numpy.sum(x134 * x138 * x97)
    result[2, 11] = numpy.sum(x149 * x151 * x98)
    result[2, 12] = numpy.sum(x135 * x147 * x155)
    result[2, 13] = numpy.sum(x154 * x39 * x98)
    result[2, 14] = numpy.sum(x137 * x159)
    result[3, 0] = numpy.sum(x160 * x163)
    result[3, 1] = numpy.sum(x166 * x167)
    result[3, 2] = numpy.sum(x162 * x167 * x51)
    result[3, 3] = numpy.sum(x168 * x173)
    result[3, 4] = numpy.sum(x119 * x166 * x45 * x83)
    result[3, 5] = numpy.sum(x102 * x162 * x168)
    result[3, 6] = numpy.sum(x176 * x54 * x85)
    result[3, 7] = numpy.sum(x173 * x177 * x51)
    result[3, 8] = numpy.sum(x102 * x166 * x177)
    result[3, 9] = numpy.sum(x127 * x162 * x178)
    result[3, 10] = numpy.sum(x179 * x183)
    result[3, 11] = numpy.sum(x176 * x184 * x66)
    result[3, 12] = numpy.sum(x121 * x16 * x186)
    result[3, 13] = numpy.sum(x127 * x166 * x187)
    result[3, 14] = numpy.sum(x106 * x16 * x189)
    result[4, 0] = numpy.sum(x108 * x190 * x29 * x48)
    result[4, 1] = numpy.sum(x113 * x150 * x191)
    result[4, 2] = numpy.sum(x128 * x143 * x191)
    result[4, 3] = numpy.sum(x126 * x138 * x192)
    result[4, 4] = numpy.sum(x143 * x193 * x45)
    result[4, 5] = numpy.sum(x128 * x147 * x192)
    result[4, 6] = numpy.sum(x123 * x150 * x194)
    result[4, 7] = numpy.sum(209.023124717498 * x116 * x151 * x54)
    result[4, 8] = numpy.sum(x147 * x193 * x54)
    result[4, 9] = numpy.sum(x128 * x153 * x194)
    result[4, 10] = numpy.sum(x132 * x184 * x190)
    result[4, 11] = numpy.sum(x123 * x143 * x195)
    result[4, 12] = numpy.sum(120.679557322504 * x116 * x147 * x187)
    result[4, 13] = numpy.sum(x113 * x153 * x195)
    result[4, 14] = numpy.sum(35.3313566383285 * x108 * x159 * x198)
    result[5, 0] = numpy.sum(x160 * x201)
    result[5, 1] = numpy.sum(x200 * x202 * x39)
    result[5, 2] = numpy.sum(x202 * x205)
    result[5, 3] = numpy.sum(x155 * x168 * x200)
    result[5, 4] = numpy.sum(x145 * x205 * x45 * x69)
    result[5, 5] = numpy.sum(x168 * x210)
    result[5, 6] = numpy.sum(x149 * x178 * x200)
    result[5, 7] = numpy.sum(x155 * x177 * x205)
    result[5, 8] = numpy.sum(x177 * x210 * x39)
    result[5, 9] = numpy.sum(x213 * x54 * x87)
    result[5, 10] = numpy.sum(x16 * x214 * x97)
    result[5, 11] = numpy.sum(x149 * x187 * x205)
    result[5, 12] = numpy.sum(x144 * x16 * x215)
    result[5, 13] = numpy.sum(x140 * x198 * x213)
    result[5, 14] = numpy.sum(x179 * x217)
    result[6, 0] = numpy.sum(x100 * x218 * x22)
    result[6, 1] = numpy.sum(x220 * x221)
    result[6, 2] = numpy.sum(x218 * x221 * x52)
    result[6, 3] = numpy.sum(x222 * x223 * x83)
    result[6, 4] = numpy.sum(x11 * x219 * x86)
    result[6, 5] = numpy.sum(x102 * x218 * x223)
    result[6, 6] = numpy.sum(x224 * x225)
    result[6, 7] = numpy.sum(x222 * x226 * x66)
    result[6, 8] = numpy.sum(x219 * x228 * x73)
    result[6, 9] = numpy.sum(x218 * x224 * x229)
    result[6, 10] = numpy.sum(
        x181 * x34 * (2.0 * x0 * (x130 + x131 + 2.0 * x174 + 2.0 * x175) + x108 * x180)
    )
    result[6, 11] = numpy.sum(x225 * x52)
    result[6, 12] = numpy.sum(x185 * x222 * x74)
    result[6, 13] = numpy.sum(x220 * x229)
    result[6, 14] = numpy.sum(x107 * x185 * x218)
    result[7, 0] = numpy.sum(x139 * x163 * x22)
    result[7, 1] = numpy.sum(x150 * x166 * x230)
    result[7, 2] = numpy.sum(x151 * x162 * x230)
    result[7, 3] = numpy.sum(x138 * x173 * x231)
    result[7, 4] = numpy.sum(x151 * x166 * x232)
    result[7, 5] = numpy.sum(x148 * x162 * x233)
    result[7, 6] = numpy.sum(53.9695387335403 * x138 * x176 * x226)
    result[7, 7] = numpy.sum(x143 * x186 * x234)
    result[7, 8] = numpy.sum(x147 * x166 * x235)
    result[7, 9] = numpy.sum(x162 * x236 * x3)
    result[7, 10] = numpy.sum(x139 * x183)
    result[7, 11] = numpy.sum(x143 * x176 * x227)
    result[7, 12] = numpy.sum(x148 * x186)
    result[7, 13] = numpy.sum(x166 * x236)
    result[7, 14] = numpy.sum(x159 * x189)
    result[8, 0] = numpy.sum(x109 * x201 * x22)
    result[8, 1] = numpy.sum(x136 * x200 * x230)
    result[8, 2] = numpy.sum(x128 * x205 * x230)
    result[8, 3] = numpy.sum(x117 * x200 * x233)
    result[8, 4] = numpy.sum(x136 * x205 * x232)
    result[8, 5] = numpy.sum(x108 * x210 * x231)
    result[8, 6] = numpy.sum(x200 * x237 * x3)
    result[8, 7] = numpy.sum(x116 * x205 * x235)
    result[8, 8] = numpy.sum(x113 * x215 * x234)
    result[8, 9] = numpy.sum(x108 * x197 * x213 * x238)
    result[8, 10] = numpy.sum(x132 * x214)
    result[8, 11] = numpy.sum(x205 * x237)
    result[8, 12] = numpy.sum(x117 * x215)
    result[8, 13] = numpy.sum(x113 * x213 * x227)
    result[8, 14] = numpy.sum(x109 * x217)
    result[9, 0] = numpy.sum(x22 * x240 * x69)
    result[9, 1] = numpy.sum(x239 * x241 * x40)
    result[9, 2] = numpy.sum(x241 * x243)
    result[9, 3] = numpy.sum(x233 * x239 * x63)
    result[9, 4] = numpy.sum(x11 * x242 * x88)
    result[9, 5] = numpy.sum(x223 * x244 * x69)
    result[9, 6] = numpy.sum(x224 * x239 * x245)
    result[9, 7] = numpy.sum(x228 * x242 * x62)
    result[9, 8] = numpy.sum(x238 * x246 * x67)
    result[9, 9] = numpy.sum(x224 * x247)
    result[9, 10] = numpy.sum(x185 * x240 * x97)
    result[9, 11] = numpy.sum(x243 * x245)
    result[9, 12] = numpy.sum(x246 * x65)
    result[9, 13] = numpy.sum(x247 * x40)
    result[9, 14] = numpy.sum(
        9.12251705727742
        * x197
        * (2.0 * x0 * (x157 + x158 + 2.0 * x211 + 2.0 * x212) + x138 * x216)
    )
    return result


def ovlp3d_40(ax, da, A, bx, db, B):
    """Cartesian 3D (gs) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 1), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = ax * bx * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = 0.5 / (ax + bx)
    x4 = x0 * (ax * A[0] + bx * B[0]) - A[0]
    x5 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x0)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x4**2 * x8
    x10 = x3 * x8
    x11 = x10 + x9
    x12 = x4 * (2.0 * x10 + x11)
    x13 = 0.564189583547756
    x14 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x15 = da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**1.5)
    x16 = x14 * x15
    x17 = x0 * x13 * x16
    x18 = 4.41641957979107 * x17
    x19 = x0 * (ax * A[1] + bx * B[1]) - A[1]
    x20 = x19 * x2
    x21 = 11.6847478934435 * x17
    x22 = x12 * x21
    x23 = x0 * (ax * A[2] + bx * B[2]) - A[2]
    x24 = x2 * x7
    x25 = x19**2 * x24
    x26 = x24 * x3
    x27 = x25 + x26
    x28 = 0.318309886183791
    x29 = 15.084944665313 * x11 * x28 * x6
    x30 = 26.1278905896872 * x17 * x23
    x31 = x14 * x7
    x32 = x23**2 * x31
    x33 = x3 * x31
    x34 = x32 + x33
    x35 = x15 * x34
    x36 = x2 * x35
    x37 = x19 * (2.0 * x26 + x27)
    x38 = x21 * x37 * x5
    x39 = x27 * x5
    x40 = x0 * x13 * x5
    x41 = x23 * (2.0 * x33 + x34)
    x42 = x15 * x2 * x40
    x43 = 11.6847478934435 * x41 * x42

    # 15 item(s)
    result[0, 0] = numpy.sum(x18 * x2 * (x12 * x4 + 3.0 * x3 * (x10 + x9)))
    result[1, 0] = numpy.sum(x20 * x22)
    result[2, 0] = numpy.sum(x2 * x22 * x23)
    result[3, 0] = numpy.sum(x16 * x27 * x29)
    result[4, 0] = numpy.sum(x11 * x20 * x30)
    result[5, 0] = numpy.sum(x29 * x36)
    result[6, 0] = numpy.sum(x38 * x4)
    result[7, 0] = numpy.sum(x30 * x39 * x4)
    result[8, 0] = numpy.sum(26.1278905896872 * x19 * x36 * x4 * x40)
    result[9, 0] = numpy.sum(x4 * x43)
    result[10, 0] = numpy.sum(x18 * x5 * (x19 * x37 + 3.0 * x3 * (x25 + x26)))
    result[11, 0] = numpy.sum(x23 * x38)
    result[12, 0] = numpy.sum(15.084944665313 * x28 * x35 * x39 * x6)
    result[13, 0] = numpy.sum(x19 * x43)
    result[14, 0] = numpy.sum(
        4.41641957979107 * x42 * (x23 * x41 + 3.0 * x3 * (x32 + x33))
    )
    return result


def ovlp3d_41(ax, da, A, bx, db, B):
    """Cartesian 3D (gp) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 3), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x2 - B[0]
    x11 = x0 * (x10 * x8 + x9)
    x12 = x0 * x8
    x13 = x10 * x9
    x14 = x12 + x13
    x15 = x14 * x3
    x16 = x3**2 * x8
    x17 = x12 + x16
    x18 = x3 * (2.0 * x12 + x17)
    x19 = 3.0 * x12
    x20 = x11 + x15
    x21 = x0 * (2.0 * x13 + x16 + x19) + x20 * x3
    x22 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x23 = 0.564189583547756
    x24 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x25 = da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**2.5)
    x26 = x24 * x25
    x27 = x1 * x23 * x26
    x28 = 8.83283915958214 * x27
    x29 = x22 * x28
    x30 = -x1 * (ax * A[1] + bx * B[1])
    x31 = -x30 - B[1]
    x32 = x29 * (x0 * (3.0 * x16 + x19) + x18 * x3)
    x33 = -x1 * (ax * A[2] + bx * B[2])
    x34 = -x33 - B[2]
    x35 = -x30 - A[1]
    x36 = x22 * x35
    x37 = x27 * x36
    x38 = 23.3694957868871 * x21
    x39 = 23.3694957868871 * x18
    x40 = x0 * x7
    x41 = x22 * x40
    x42 = x22 * x7
    x43 = x35 * x42
    x44 = x31 * x43
    x45 = x41 + x44
    x46 = 0.318309886183791 * x6
    x47 = x26 * x46
    x48 = x45 * x47
    x49 = x22 * x39
    x50 = x27 * x34
    x51 = -x33 - A[2]
    x52 = x27 * x51
    x53 = x24 * x40
    x54 = x24 * x7
    x55 = x51 * x54
    x56 = x34 * x55
    x57 = x53 + x56
    x58 = x25 * x46
    x59 = x57 * x58
    x60 = x35**2 * x42
    x61 = x41 + x60
    x62 = 30.169889330626 * x47
    x63 = x61 * x62
    x64 = x0 * (x31 * x42 + x43)
    x65 = x35 * x45
    x66 = x64 + x65
    x67 = 52.2557811793745 * x17
    x68 = x51**2 * x54
    x69 = x53 + x68
    x70 = x22 * x58
    x71 = 30.169889330626 * x70
    x72 = x69 * x71
    x73 = x0 * (x34 * x54 + x55)
    x74 = x51 * x57
    x75 = x73 + x74
    x76 = x35 * (2.0 * x41 + x61)
    x77 = 23.3694957868871 * x14
    x78 = 3.0 * x41
    x79 = x0 * (2.0 * x44 + x60 + x78) + x35 * x66
    x80 = 23.3694957868871 * x5
    x81 = x3 * x80
    x82 = x76 * x80
    x83 = 52.2557811793745 * x14
    x84 = x5 * x66
    x85 = 52.2557811793745 * x3
    x86 = x5 * x61
    x87 = x45 * x58
    x88 = x1 * x22 * x23 * x25
    x89 = x5 * x88
    x90 = x51 * (2.0 * x53 + x69)
    x91 = x80 * x90
    x92 = x88 * x91
    x93 = 3.0 * x53
    x94 = x0 * (2.0 * x56 + x68 + x93) + x51 * x75
    x95 = x88 * x94
    x96 = x28 * x5
    x97 = x96 * (x0 * (3.0 * x60 + x78) + x35 * x76)
    x98 = 30.169889330626 * x58
    x99 = x69 * x98
    x100 = 8.83283915958214 * x89
    x101 = x100 * (x0 * (3.0 * x68 + x93) + x51 * x90)

    # 45 item(s)
    result[0, 0] = numpy.sum(x29 * (x0 * (3.0 * x11 + 3.0 * x15 + x18) + x21 * x3))
    result[0, 1] = numpy.sum(x31 * x32)
    result[0, 2] = numpy.sum(x32 * x34)
    result[1, 0] = numpy.sum(x37 * x38)
    result[1, 1] = numpy.sum(x39 * x48)
    result[1, 2] = numpy.sum(x35 * x49 * x50)
    result[2, 0] = numpy.sum(x22 * x38 * x52)
    result[2, 1] = numpy.sum(x31 * x49 * x52)
    result[2, 2] = numpy.sum(x49 * x59)
    result[3, 0] = numpy.sum(x20 * x63)
    result[3, 1] = numpy.sum(x17 * x62 * x66)
    result[3, 2] = numpy.sum(x17 * x34 * x63)
    result[4, 0] = numpy.sum(52.2557811793745 * x20 * x37 * x51)
    result[4, 1] = numpy.sum(x48 * x51 * x67)
    result[4, 2] = numpy.sum(x36 * x59 * x67)
    result[5, 0] = numpy.sum(x20 * x72)
    result[5, 1] = numpy.sum(x17 * x31 * x72)
    result[5, 2] = numpy.sum(x17 * x71 * x75)
    result[6, 0] = numpy.sum(x47 * x76 * x77)
    result[6, 1] = numpy.sum(x27 * x79 * x81)
    result[6, 2] = numpy.sum(x3 * x50 * x82)
    result[7, 0] = numpy.sum(x47 * x51 * x61 * x83)
    result[7, 1] = numpy.sum(x52 * x84 * x85)
    result[7, 2] = numpy.sum(x59 * x85 * x86)
    result[8, 0] = numpy.sum(x35 * x69 * x70 * x83)
    result[8, 1] = numpy.sum(x5 * x69 * x85 * x87)
    result[8, 2] = numpy.sum(x35 * x75 * x85 * x89)
    result[9, 0] = numpy.sum(x70 * x77 * x90)
    result[9, 1] = numpy.sum(x3 * x31 * x92)
    result[9, 2] = numpy.sum(x81 * x95)
    result[10, 0] = numpy.sum(x10 * x97)
    result[10, 1] = numpy.sum(x96 * (x0 * (3.0 * x64 + 3.0 * x65 + x76) + x35 * x79))
    result[10, 2] = numpy.sum(x34 * x97)
    result[11, 0] = numpy.sum(x10 * x52 * x82)
    result[11, 1] = numpy.sum(x52 * x79 * x80)
    result[11, 2] = numpy.sum(x59 * x82)
    result[12, 0] = numpy.sum(x10 * x86 * x99)
    result[12, 1] = numpy.sum(x84 * x99)
    result[12, 2] = numpy.sum(x75 * x86 * x98)
    result[13, 0] = numpy.sum(x10 * x35 * x92)
    result[13, 1] = numpy.sum(x87 * x91)
    result[13, 2] = numpy.sum(x35 * x80 * x95)
    result[14, 0] = numpy.sum(x10 * x101)
    result[14, 1] = numpy.sum(x101 * x31)
    result[14, 2] = numpy.sum(x100 * (x0 * (3.0 * x73 + 3.0 * x74 + x90) + x51 * x94))
    return result


def ovlp3d_42(ax, da, A, bx, db, B):
    """Cartesian 3D (gd) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 6), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = 3.0 * x10
    x12 = -x2 - B[0]
    x13 = x3 * x8
    x14 = x12 * x13
    x15 = x11 + 2.0 * x14
    x16 = x0 * (x15 + x9)
    x17 = x12 * x8
    x18 = x0 * (x13 + x17)
    x19 = x10 + x14
    x20 = x19 * x3
    x21 = x18 + x20
    x22 = x21 * x3
    x23 = x12**2 * x8
    x24 = x0 * (x15 + x23)
    x25 = x10 + x23
    x26 = x25 * x3
    x27 = x0 * x17
    x28 = x26 + 2.0 * x27
    x29 = x28 * x3
    x30 = x24 + x29
    x31 = 2.0 * x0 * (x18 + x20 + x26 + 2.0 * x27) + x3 * x30
    x32 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x33 = da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**3.5)
    x34 = x32 * x33
    x35 = 10.1992841329868 * x34
    x36 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x37 = 0.564189583547756 * x1
    x38 = x36 * x37
    x39 = -x1 * (ax * A[1] + bx * B[1])
    x40 = -x39 - B[1]
    x41 = 17.6656783191643 * x40
    x42 = x10 + x9
    x43 = 2.0 * x0 * x13 + x3 * x42
    x44 = x16 + x22
    x45 = x34 * x38
    x46 = x45 * (x0 * (3.0 * x18 + 3.0 * x20 + x43) + x3 * x44)
    x47 = -x1 * (ax * A[2] + bx * B[2])
    x48 = -x47 - B[2]
    x49 = 17.6656783191643 * x48
    x50 = x36 * x7
    x51 = x40**2 * x50
    x52 = x0 * x50
    x53 = x51 + x52
    x54 = x0 * (x11 + 3.0 * x9) + x3 * x43
    x55 = 0.318309886183791 * x6
    x56 = x54 * x55
    x57 = x45 * x48
    x58 = x33 * x36
    x59 = x32 * x7
    x60 = x48**2 * x59
    x61 = x0 * x59
    x62 = x60 + x61
    x63 = 10.1992841329868 * x62
    x64 = -x39 - A[1]
    x65 = 26.9847693667702 * x64
    x66 = x31 * x45
    x67 = 46.7389915737742 * x44
    x68 = x50 * x64
    x69 = x40 * x68
    x70 = x52 + x69
    x71 = x34 * x55
    x72 = x70 * x71
    x73 = 26.9847693667702 * x43
    x74 = x53 * x64
    x75 = 2.0 * x52
    x76 = x40 * x75 + x74
    x77 = x71 * x76
    x78 = 46.7389915737742 * x43
    x79 = x55 * x58
    x80 = x73 * x79
    x81 = -x47 - A[2]
    x82 = 26.9847693667702 * x81
    x83 = x45 * x81
    x84 = x59 * x81
    x85 = x48 * x84
    x86 = x61 + x85
    x87 = x79 * x86
    x88 = x71 * x81
    x89 = x62 * x81
    x90 = 2.0 * x61
    x91 = x48 * x90 + x89
    x92 = x50 * x64**2
    x93 = x52 + x92
    x94 = 34.8371874529163 * x93
    x95 = 60.3397786612521 * x21
    x96 = x0 * (x40 * x50 + x68)
    x97 = x64 * x70
    x98 = x96 + x97
    x99 = x71 * x98
    x100 = 3.0 * x52
    x101 = x100 + 2.0 * x69
    x102 = x0 * (x101 + x51)
    x103 = x64 * x76
    x104 = x102 + x103
    x105 = 34.8371874529163 * x42
    x106 = 0.179587122125167 * x33
    x107 = x106 * x42
    x108 = 60.3397786612521 * x64
    x109 = 104.511562358749 * x70
    x110 = 104.511562358749 * x86
    x111 = x64 * x79
    x112 = 60.3397786612521 * x42
    x113 = x59 * x81**2
    x114 = x113 + x61
    x115 = 34.8371874529163 * x114
    x116 = x79 * x95
    x117 = x0 * (x48 * x59 + x84)
    x118 = x81 * x86
    x119 = x117 + x118
    x120 = 60.3397786612521 * x119
    x121 = 3.0 * x61
    x122 = x121 + 2.0 * x85
    x123 = x0 * (x122 + x60)
    x124 = x81 * x91
    x125 = x123 + x124
    x126 = 26.9847693667702 * x28
    x127 = x64 * (x75 + x93)
    x128 = x127 * x71
    x129 = x0 * (x101 + x92)
    x130 = x64 * x98
    x131 = x129 + x130
    x132 = 46.7389915737742 * x19
    x133 = 26.9847693667702 * x3
    x134 = 2.0 * x0 * (2.0 * x40 * x52 + x74 + x96 + x97) + x104 * x64
    x135 = x37 * x5
    x136 = x135 * x34
    x137 = x134 * x136
    x138 = 46.7389915737742 * x131
    x139 = x136 * x3
    x140 = x33 * x5
    x141 = x140 * x55
    x142 = 26.9847693667702 * x141
    x143 = x127 * x142
    x144 = 60.3397786612521 * x93
    x145 = 104.511562358749 * x19
    x146 = x106 * x19
    x147 = x141 * x98
    x148 = x141 * x3
    x149 = 60.3397786612521 * x114
    x150 = x140 * x38
    x151 = x150 * x3
    x152 = x81 * (x114 + x90)
    x153 = x152 * x79
    x154 = x0 * (x113 + x122)
    x155 = x119 * x81
    x156 = x154 + x155
    x157 = x142 * x152
    x158 = 46.7389915737742 * x156
    x159 = 2.0 * x0 * (x117 + x118 + 2.0 * x48 * x61 + x89) + x125 * x81
    x160 = x150 * x159
    x161 = x0 * (x100 + 3.0 * x92) + x127 * x64
    x162 = x161 * x55
    x163 = x136 * (x0 * (x127 + 3.0 * x96 + 3.0 * x97) + x131 * x64)
    x164 = 17.6656783191643 * x12
    x165 = x12 * x136
    x166 = x141 * x86
    x167 = 46.7389915737742 * x12
    x168 = x141 * x70
    x169 = x12 * x150
    x170 = x0 * (3.0 * x113 + x121) + x152 * x81
    x171 = 10.1992841329868 * x170
    x172 = x150 * (x0 * (3.0 * x117 + 3.0 * x118 + x152) + x156 * x81)

    # 90 item(s)
    result[0, 0] = numpy.sum(
        x35 * x38 * (x0 * (2.0 * x16 + 2.0 * x22 + 3.0 * x24 + 3.0 * x29) + x3 * x31)
    )
    result[0, 1] = numpy.sum(x41 * x46)
    result[0, 2] = numpy.sum(x46 * x49)
    result[0, 3] = numpy.sum(x35 * x53 * x56)
    result[0, 4] = numpy.sum(x41 * x54 * x57)
    result[0, 5] = numpy.sum(x56 * x58 * x63)
    result[1, 0] = numpy.sum(x65 * x66)
    result[1, 1] = numpy.sum(x67 * x72)
    result[1, 2] = numpy.sum(x57 * x64 * x67)
    result[1, 3] = numpy.sum(x73 * x77)
    result[1, 4] = numpy.sum(x48 * x72 * x78)
    result[1, 5] = numpy.sum(x62 * x64 * x80)
    result[2, 0] = numpy.sum(x66 * x82)
    result[2, 1] = numpy.sum(x40 * x67 * x83)
    result[2, 2] = numpy.sum(x67 * x87)
    result[2, 3] = numpy.sum(x53 * x73 * x88)
    result[2, 4] = numpy.sum(x40 * x78 * x87)
    result[2, 5] = numpy.sum(x80 * x91)
    result[3, 0] = numpy.sum(x30 * x71 * x94)
    result[3, 1] = numpy.sum(x95 * x99)
    result[3, 2] = numpy.sum(x48 * x71 * x93 * x95)
    result[3, 3] = numpy.sum(x104 * x105 * x71)
    result[3, 4] = numpy.sum(60.3397786612521 * x42 * x48 * x99)
    result[3, 5] = numpy.sum(x107 * x62 * x94)
    result[4, 0] = numpy.sum(x108 * x30 * x83)
    result[4, 1] = numpy.sum(x109 * x21 * x88)
    result[4, 2] = numpy.sum(x110 * x111 * x21)
    result[4, 3] = numpy.sum(x112 * x77 * x81)
    result[4, 4] = numpy.sum(x107 * x110 * x70)
    result[4, 5] = numpy.sum(x111 * x112 * x91)
    result[5, 0] = numpy.sum(x115 * x30 * x79)
    result[5, 1] = numpy.sum(x114 * x116 * x40)
    result[5, 2] = numpy.sum(x116 * x119)
    result[5, 3] = numpy.sum(x107 * x115 * x53)
    result[5, 4] = numpy.sum(x120 * x40 * x42 * x79)
    result[5, 5] = numpy.sum(x105 * x125 * x79)
    result[6, 0] = numpy.sum(x126 * x128)
    result[6, 1] = numpy.sum(x131 * x132 * x71)
    result[6, 2] = numpy.sum(x128 * x132 * x48)
    result[6, 3] = numpy.sum(x133 * x137)
    result[6, 4] = numpy.sum(x138 * x139 * x48)
    result[6, 5] = numpy.sum(x143 * x3 * x62)
    result[7, 0] = numpy.sum(x144 * x28 * x88)
    result[7, 1] = numpy.sum(x145 * x81 * x99)
    result[7, 2] = numpy.sum(x110 * x146 * x93)
    result[7, 3] = numpy.sum(60.3397786612521 * x104 * x139 * x81)
    result[7, 4] = numpy.sum(x110 * x147 * x3)
    result[7, 5] = numpy.sum(x144 * x148 * x91)
    result[8, 0] = numpy.sum(x111 * x149 * x28)
    result[8, 1] = numpy.sum(x109 * x114 * x146)
    result[8, 2] = numpy.sum(x111 * x119 * x145)
    result[8, 3] = numpy.sum(x148 * x149 * x76)
    result[8, 4] = numpy.sum(x109 * x119 * x148)
    result[8, 5] = numpy.sum(x108 * x125 * x151)
    result[9, 0] = numpy.sum(x126 * x153)
    result[9, 1] = numpy.sum(x132 * x153 * x40)
    result[9, 2] = numpy.sum(x132 * x156 * x79)
    result[9, 3] = numpy.sum(x157 * x3 * x53)
    result[9, 4] = numpy.sum(x151 * x158 * x40)
    result[9, 5] = numpy.sum(x133 * x160)
    result[10, 0] = numpy.sum(x162 * x25 * x35)
    result[10, 1] = numpy.sum(x163 * x164)
    result[10, 2] = numpy.sum(x161 * x165 * x49)
    result[10, 3] = numpy.sum(
        x135
        * x35
        * (x0 * (3.0 * x102 + 3.0 * x103 + 2.0 * x129 + 2.0 * x130) + x134 * x64)
    )
    result[10, 4] = numpy.sum(x163 * x49)
    result[10, 5] = numpy.sum(x140 * x162 * x63)
    result[11, 0] = numpy.sum(x128 * x25 * x82)
    result[11, 1] = numpy.sum(x138 * x165 * x81)
    result[11, 2] = numpy.sum(x127 * x166 * x167)
    result[11, 3] = numpy.sum(x137 * x82)
    result[11, 4] = numpy.sum(x138 * x166)
    result[11, 5] = numpy.sum(x143 * x91)
    result[12, 0] = numpy.sum(x106 * x114 * x25 * x94)
    result[12, 1] = numpy.sum(60.3397786612521 * x114 * x12 * x147)
    result[12, 2] = numpy.sum(x12 * x120 * x141 * x93)
    result[12, 3] = numpy.sum(x104 * x115 * x141)
    result[12, 4] = numpy.sum(x120 * x147)
    result[12, 5] = numpy.sum(x125 * x141 * x94)
    result[13, 0] = numpy.sum(x153 * x25 * x65)
    result[13, 1] = numpy.sum(x152 * x167 * x168)
    result[13, 2] = numpy.sum(x158 * x169 * x64)
    result[13, 3] = numpy.sum(x157 * x76)
    result[13, 4] = numpy.sum(x158 * x168)
    result[13, 5] = numpy.sum(x160 * x65)
    result[14, 0] = numpy.sum(x171 * x25 * x79)
    result[14, 1] = numpy.sum(x169 * x170 * x41)
    result[14, 2] = numpy.sum(x164 * x172)
    result[14, 3] = numpy.sum(x141 * x171 * x53)
    result[14, 4] = numpy.sum(x172 * x41)
    result[14, 5] = numpy.sum(
        10.1992841329868
        * x150
        * (x0 * (3.0 * x123 + 3.0 * x124 + 2.0 * x154 + 2.0 * x155) + x159 * x81)
    )
    return result


def ovlp3d_43(ax, da, A, bx, db, B):
    """Cartesian 3D (gf) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 10), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = x10 + x9
    x12 = x11 * x3
    x13 = -x2 - A[0]
    x14 = x11 * x13
    x15 = x10 * x3
    x16 = x0 * (x12 + 3.0 * x14 + 8.0 * x15)
    x17 = x13 * x8
    x18 = x0 * (x17 + x3 * x8)
    x19 = x17 * x3
    x20 = x10 + x19
    x21 = x13 * x20
    x22 = 2.0 * x0 * (x14 + 2.0 * x15 + x18 + x21)
    x23 = 3.0 * x10
    x24 = x0 * (x23 + 3.0 * x9)
    x25 = 2.0 * x15
    x26 = x12 + x25
    x27 = x13 * x26
    x28 = x24 + x27
    x29 = x13 * x28
    x30 = 2.0 * x19 + x23
    x31 = x0 * (x30 + x9)
    x32 = x14 + x25
    x33 = x13 * x32
    x34 = x31 + x33
    x35 = x13 * x34
    x36 = 3.0 * x31 + 3.0 * x33
    x37 = x16 + x29
    x38 = x0 * (2.0 * x24 + 2.0 * x27 + x36) + x13 * x37
    x39 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x40 = da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**4.5)
    x41 = x39 * x40
    x42 = 9.12251705727742 * x41
    x43 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x44 = 0.564189583547756 * x1
    x45 = x43 * x44
    x46 = -x1 * (ax * A[1] + bx * B[1])
    x47 = -x46 - B[1]
    x48 = x45 * x47
    x49 = x13**2 * x8
    x50 = x0 * (x30 + x49)
    x51 = x18 + x21
    x52 = x13 * x51
    x53 = x22 + x35
    x54 = 20.3985682659737 * x41
    x55 = x54 * (x0 * (x36 + 2.0 * x50 + 2.0 * x52) + x13 * x53)
    x56 = -x1 * (ax * A[2] + bx * B[2])
    x57 = -x56 - B[2]
    x58 = x45 * x57
    x59 = x10 + x49
    x60 = x13 * (2.0 * x10 + x59)
    x61 = x50 + x52
    x62 = x0 * (3.0 * x18 + 3.0 * x21 + x60) + x13 * x61
    x63 = x43 * x7
    x64 = x47**2 * x63
    x65 = x0 * x63
    x66 = x64 + x65
    x67 = 0.318309886183791 * x6
    x68 = x66 * x67
    x69 = x54 * x68
    x70 = x41 * x57
    x71 = 35.3313566383285 * x70
    x72 = x40 * x67
    x73 = x43 * x72
    x74 = x39 * x7
    x75 = x57**2 * x74
    x76 = x0 * x74
    x77 = x75 + x76
    x78 = 20.3985682659737 * x77
    x79 = x73 * x78
    x80 = x47 * x66
    x81 = x47 * x65
    x82 = 2.0 * x81
    x83 = x80 + x82
    x84 = x0 * (x23 + 3.0 * x49) + x13 * x60
    x85 = x42 * x67
    x86 = x57 * x77
    x87 = x57 * x76
    x88 = 2.0 * x87
    x89 = x86 + x88
    x90 = 9.12251705727742 * x89
    x91 = -x46 - A[1]
    x92 = x45 * x91
    x93 = 24.1359114645008 * x41
    x94 = x38 * x93
    x95 = x63 * x91
    x96 = x47 * x95
    x97 = x65 + x96
    x98 = 53.9695387335403 * x97
    x99 = x41 * x67
    x100 = 53.9695387335403 * x53
    x101 = 53.9695387335403 * x61
    x102 = x66 * x91
    x103 = x102 + x82
    x104 = x103 * x99
    x105 = 93.4779831475484 * x61
    x106 = x57 * x99
    x107 = x101 * x73
    x108 = 24.1359114645008 * x60
    x109 = 3.0 * x65
    x110 = x0 * (x109 + 3.0 * x64)
    x111 = x83 * x91
    x112 = x110 + x111
    x113 = x112 * x99
    x114 = 53.9695387335403 * x60
    x115 = 0.179587122125167 * x40
    x116 = x115 * x77
    x117 = x108 * x73
    x118 = -x56 - A[2]
    x119 = x118 * x41
    x120 = x118 * x74
    x121 = x120 * x57
    x122 = x121 + x76
    x123 = 53.9695387335403 * x73
    x124 = x47 * x73
    x125 = x118 * x77
    x126 = x125 + x88
    x127 = x118 * x99
    x128 = x115 * x122
    x129 = 3.0 * x76
    x130 = x0 * (x129 + 3.0 * x75)
    x131 = x118 * x89
    x132 = x130 + x131
    x133 = x63 * x91**2
    x134 = x133 + x65
    x135 = 31.1593277158494 * x134
    x136 = x0 * (x47 * x63 + x95)
    x137 = x91 * x97
    x138 = x136 + x137
    x139 = 69.6743749058326 * x99
    x140 = x139 * x34
    x141 = x109 + 2.0 * x96
    x142 = x0 * (x141 + x64)
    x143 = x103 * x91
    x144 = x142 + x143
    x145 = x139 * x144
    x146 = 120.679557322504 * x51
    x147 = 69.6743749058326 * x116
    x148 = x0 * (3.0 * x102 + x80 + 8.0 * x81)
    x149 = x112 * x91
    x150 = x148 + x149
    x151 = 31.1593277158494 * x59
    x152 = x115 * x59
    x153 = 53.9695387335403 * x119
    x154 = 120.679557322504 * x34
    x155 = x73 * x91
    x156 = 120.679557322504 * x118
    x157 = 209.023124717498 * x128
    x158 = 120.679557322504 * x126
    x159 = 120.679557322504 * x128
    x160 = x123 * x91
    x161 = x118**2 * x74
    x162 = x161 + x76
    x163 = 31.1593277158494 * x162
    x164 = 69.6743749058326 * x73
    x165 = x164 * x34
    x166 = x0 * (x120 + x57 * x74)
    x167 = x118 * x122
    x168 = x166 + x167
    x169 = 69.6743749058326 * x66
    x170 = x115 * x162
    x171 = 2.0 * x121 + x129
    x172 = x0 * (x171 + x75)
    x173 = x118 * x126
    x174 = x172 + x173
    x175 = x164 * x174
    x176 = x0 * (3.0 * x125 + x86 + 8.0 * x87)
    x177 = x118 * x132
    x178 = x176 + x177
    x179 = 24.1359114645008 * x28
    x180 = x91 * (x134 + 2.0 * x65)
    x181 = x180 * x99
    x182 = x0 * (x133 + x141)
    x183 = x138 * x91
    x184 = x182 + x183
    x185 = 53.9695387335403 * x184
    x186 = x32 * x99
    x187 = 53.9695387335403 * x180
    x188 = 2.0 * x0 * (x102 + x136 + x137 + 2.0 * x81)
    x189 = x144 * x91
    x190 = x188 + x189
    x191 = 53.9695387335403 * x20
    x192 = 93.4779831475484 * x20
    x193 = x44 * x5
    x194 = x13 * x193
    x195 = 3.0 * x142 + 3.0 * x143
    x196 = x0 * (2.0 * x110 + 2.0 * x111 + x195) + x150 * x91
    x197 = x196 * x93
    x198 = x5 * x72
    x199 = x185 * x198
    x200 = 24.1359114645008 * x198
    x201 = x180 * x200
    x202 = 120.679557322504 * x144
    x203 = x115 * x20
    x204 = x13 * x198
    x205 = 53.9695387335403 * x198
    x206 = x13 * x205
    x207 = 120.679557322504 * x170
    x208 = 120.679557322504 * x155
    x209 = 120.679557322504 * x204
    x210 = x40 * x5
    x211 = x210 * x45
    x212 = 53.9695387335403 * x13
    x213 = x118 * (x162 + 2.0 * x76)
    x214 = x213 * x73
    x215 = x123 * x32
    x216 = x0 * (x161 + x171)
    x217 = x118 * x168
    x218 = x216 + x217
    x219 = x115 * x213
    x220 = 2.0 * x0 * (x125 + x166 + x167 + 2.0 * x87)
    x221 = x118 * x174
    x222 = x220 + x221
    x223 = x200 * x213
    x224 = x210 * x68
    x225 = x211 * x47
    x226 = 3.0 * x172 + 3.0 * x173
    x227 = x0 * (2.0 * x130 + 2.0 * x131 + x226) + x118 * x178
    x228 = 24.1359114645008 * x211 * x227
    x229 = x0 * (x109 + 3.0 * x133) + x180 * x91
    x230 = x0 * (3.0 * x136 + 3.0 * x137 + x180) + x184 * x91
    x231 = x11 * x54 * x67
    x232 = x193 * x3
    x233 = x54 * (x0 * (2.0 * x182 + 2.0 * x183 + x195) + x190 * x91)
    x234 = x198 * x78
    x235 = 24.1359114645008 * x26
    x236 = x198 * x3
    x237 = 93.4779831475484 * x236
    x238 = 69.6743749058326 * x11
    x239 = 69.6743749058326 * x198
    x240 = x144 * x239
    x241 = x174 * x239
    x242 = x103 * x205
    x243 = x211 * x3
    x244 = x0 * (x129 + 3.0 * x161) + x118 * x213
    x245 = 9.12251705727742 * x244
    x246 = x11 * x73
    x247 = 20.3985682659737 * x244
    x248 = x0 * (3.0 * x166 + 3.0 * x167 + x213) + x118 * x218
    x249 = 20.3985682659737 * x248
    x250 = (
        20.3985682659737 * x0 * (2.0 * x216 + 2.0 * x217 + x226)
        + 20.3985682659737 * x118 * x222
    )

    # 150 item(s)
    result[0, 0] = numpy.sum(x42 * x45 * (3.0 * x0 * (x16 + x22 + x29 + x35) + x13 * x38))
    result[0, 1] = numpy.sum(x48 * x55)
    result[0, 2] = numpy.sum(x55 * x58)
    result[0, 3] = numpy.sum(x62 * x69)
    result[0, 4] = numpy.sum(x48 * x62 * x71)
    result[0, 5] = numpy.sum(x62 * x79)
    result[0, 6] = numpy.sum(x83 * x84 * x85)
    result[0, 7] = numpy.sum(x57 * x69 * x84)
    result[0, 8] = numpy.sum(x47 * x79 * x84)
    result[0, 9] = numpy.sum(x73 * x84 * x90)
    result[1, 0] = numpy.sum(x92 * x94)
    result[1, 1] = numpy.sum(x53 * x98 * x99)
    result[1, 2] = numpy.sum(x100 * x41 * x58 * x91)
    result[1, 3] = numpy.sum(x101 * x104)
    result[1, 4] = numpy.sum(x105 * x106 * x97)
    result[1, 5] = numpy.sum(x107 * x77 * x91)
    result[1, 6] = numpy.sum(x108 * x113)
    result[1, 7] = numpy.sum(x104 * x114 * x57)
    result[1, 8] = numpy.sum(x116 * x60 * x98)
    result[1, 9] = numpy.sum(x117 * x89 * x91)
    result[2, 0] = numpy.sum(x118 * x45 * x94)
    result[2, 1] = numpy.sum(x100 * x119 * x48)
    result[2, 2] = numpy.sum(x122 * x123 * x53)
    result[2, 3] = numpy.sum(x101 * x119 * x68)
    result[2, 4] = numpy.sum(x105 * x122 * x124)
    result[2, 5] = numpy.sum(x107 * x126)
    result[2, 6] = numpy.sum(x108 * x127 * x83)
    result[2, 7] = numpy.sum(x114 * x128 * x66)
    result[2, 8] = numpy.sum(x114 * x124 * x126)
    result[2, 9] = numpy.sum(x117 * x132)
    result[3, 0] = numpy.sum(x135 * x37 * x99)
    result[3, 1] = numpy.sum(x138 * x140)
    result[3, 2] = numpy.sum(x134 * x140 * x57)
    result[3, 3] = numpy.sum(x145 * x51)
    result[3, 4] = numpy.sum(x106 * x138 * x146)
    result[3, 5] = numpy.sum(x134 * x147 * x51)
    result[3, 6] = numpy.sum(x150 * x151 * x99)
    result[3, 7] = numpy.sum(x145 * x57 * x59)
    result[3, 8] = numpy.sum(x138 * x147 * x59)
    result[3, 9] = numpy.sum(x135 * x152 * x89)
    result[4, 0] = numpy.sum(x153 * x37 * x92)
    result[4, 1] = numpy.sum(x127 * x154 * x97)
    result[4, 2] = numpy.sum(x122 * x154 * x155)
    result[4, 3] = numpy.sum(x104 * x156 * x51)
    result[4, 4] = numpy.sum(x157 * x51 * x97)
    result[4, 5] = numpy.sum(x155 * x158 * x51)
    result[4, 6] = numpy.sum(53.9695387335403 * x113 * x118 * x59)
    result[4, 7] = numpy.sum(x103 * x159 * x59)
    result[4, 8] = numpy.sum(x152 * x158 * x97)
    result[4, 9] = numpy.sum(x132 * x160 * x59)
    result[5, 0] = numpy.sum(x163 * x37 * x73)
    result[5, 1] = numpy.sum(x162 * x165 * x47)
    result[5, 2] = numpy.sum(x165 * x168)
    result[5, 3] = numpy.sum(x169 * x170 * x51)
    result[5, 4] = numpy.sum(x124 * x146 * x168)
    result[5, 5] = numpy.sum(x175 * x51)
    result[5, 6] = numpy.sum(x152 * x163 * x83)
    result[5, 7] = numpy.sum(x152 * x168 * x169)
    result[5, 8] = numpy.sum(x175 * x47 * x59)
    result[5, 9] = numpy.sum(x151 * x178 * x73)
    result[6, 0] = numpy.sum(x179 * x181)
    result[6, 1] = numpy.sum(x185 * x186)
    result[6, 2] = numpy.sum(x186 * x187 * x57)
    result[6, 3] = numpy.sum(x190 * x191 * x99)
    result[6, 4] = numpy.sum(x106 * x184 * x192)
    result[6, 5] = numpy.sum(x116 * x180 * x191)
    result[6, 6] = numpy.sum(x194 * x197)
    result[6, 7] = numpy.sum(53.9695387335403 * x190 * x194 * x70)
    result[6, 8] = numpy.sum(x13 * x199 * x77)
    result[6, 9] = numpy.sum(x13 * x201 * x89)
    result[7, 0] = numpy.sum(53.9695387335403 * x127 * x134 * x28)
    result[7, 1] = numpy.sum(x138 * x156 * x186)
    result[7, 2] = numpy.sum(x134 * x159 * x32)
    result[7, 3] = numpy.sum(x127 * x20 * x202)
    result[7, 4] = numpy.sum(x138 * x157 * x20)
    result[7, 5] = numpy.sum(x134 * x158 * x203)
    result[7, 6] = numpy.sum(x150 * x153 * x194)
    result[7, 7] = numpy.sum(x122 * x202 * x204)
    result[7, 8] = numpy.sum(x138 * x158 * x204)
    result[7, 9] = numpy.sum(x132 * x134 * x206)
    result[8, 0] = numpy.sum(x160 * x162 * x28)
    result[8, 1] = numpy.sum(x207 * x32 * x97)
    result[8, 2] = numpy.sum(x168 * x208 * x32)
    result[8, 3] = numpy.sum(x103 * x20 * x207)
    result[8, 4] = numpy.sum(209.023124717498 * x168 * x203 * x97)
    result[8, 5] = numpy.sum(x174 * x20 * x208)
    result[8, 6] = numpy.sum(x112 * x162 * x206)
    result[8, 7] = numpy.sum(x103 * x168 * x209)
    result[8, 8] = numpy.sum(x174 * x209 * x97)
    result[8, 9] = numpy.sum(x178 * x211 * x212 * x91)
    result[9, 0] = numpy.sum(x179 * x214)
    result[9, 1] = numpy.sum(x213 * x215 * x47)
    result[9, 2] = numpy.sum(x215 * x218)
    result[9, 3] = numpy.sum(x191 * x219 * x66)
    result[9, 4] = numpy.sum(x124 * x192 * x218)
    result[9, 5] = numpy.sum(x191 * x222 * x73)
    result[9, 6] = numpy.sum(x13 * x223 * x83)
    result[9, 7] = numpy.sum(x212 * x218 * x224)
    result[9, 8] = numpy.sum(x212 * x222 * x225)
    result[9, 9] = numpy.sum(x13 * x228)
    result[10, 0] = numpy.sum(x229 * x26 * x85)
    result[10, 1] = numpy.sum(x230 * x231)
    result[10, 2] = numpy.sum(x229 * x231 * x57)
    result[10, 3] = numpy.sum(x232 * x233)
    result[10, 4] = numpy.sum(x230 * x232 * x71)
    result[10, 5] = numpy.sum(x229 * x234 * x3)
    result[10, 6] = numpy.sum(
        x193 * x42 * (3.0 * x0 * (x148 + x149 + x188 + x189) + x196 * x91)
    )
    result[10, 7] = numpy.sum(x193 * x233 * x57)
    result[10, 8] = numpy.sum(x230 * x234)
    result[10, 9] = numpy.sum(x198 * x229 * x90)
    result[11, 0] = numpy.sum(x118 * x181 * x235)
    result[11, 1] = numpy.sum(x11 * x127 * x185)
    result[11, 2] = numpy.sum(x11 * x128 * x187)
    result[11, 3] = numpy.sum(x153 * x190 * x232)
    result[11, 4] = numpy.sum(x122 * x184 * x237)
    result[11, 5] = numpy.sum(x126 * x187 * x236)
    result[11, 6] = numpy.sum(x118 * x193 * x197)
    result[11, 7] = numpy.sum(x122 * x190 * x205)
    result[11, 8] = numpy.sum(x126 * x199)
    result[11, 9] = numpy.sum(x132 * x201)
    result[12, 0] = numpy.sum(x135 * x170 * x26)
    result[12, 1] = numpy.sum(x138 * x170 * x238)
    result[12, 2] = numpy.sum(x115 * x134 * x168 * x238)
    result[12, 3] = numpy.sum(x162 * x240 * x3)
    result[12, 4] = numpy.sum(120.679557322504 * x138 * x168 * x236)
    result[12, 5] = numpy.sum(x134 * x241 * x3)
    result[12, 6] = numpy.sum(x150 * x163 * x198)
    result[12, 7] = numpy.sum(x168 * x240)
    result[12, 8] = numpy.sum(x138 * x241)
    result[12, 9] = numpy.sum(x135 * x178 * x198)
    result[13, 0] = numpy.sum(x214 * x235 * x91)
    result[13, 1] = numpy.sum(x11 * x219 * x98)
    result[13, 2] = numpy.sum(x11 * x160 * x218)
    result[13, 3] = numpy.sum(x213 * x242 * x3)
    result[13, 4] = numpy.sum(x218 * x237 * x97)
    result[13, 5] = numpy.sum(53.9695387335403 * x222 * x243 * x91)
    result[13, 6] = numpy.sum(x112 * x223)
    result[13, 7] = numpy.sum(x218 * x242)
    result[13, 8] = numpy.sum(x198 * x222 * x98)
    result[13, 9] = numpy.sum(x228 * x91)
    result[14, 0] = numpy.sum(x245 * x26 * x73)
    result[14, 1] = numpy.sum(x246 * x247 * x47)
    result[14, 2] = numpy.sum(x246 * x249)
    result[14, 3] = numpy.sum(x224 * x247 * x3)
    result[14, 4] = numpy.sum(35.3313566383285 * x243 * x248 * x47)
    result[14, 5] = numpy.sum(x243 * x250)
    result[14, 6] = numpy.sum(x198 * x245 * x83)
    result[14, 7] = numpy.sum(x224 * x249)
    result[14, 8] = numpy.sum(x225 * x250)
    result[14, 9] = numpy.sum(
        9.12251705727742 * x211 * (3.0 * x0 * (x176 + x177 + x220 + x221) + x118 * x227)
    )
    return result


def ovlp3d_44(ax, da, A, bx, db, B):
    """Cartesian 3D (gg) overlap integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 15), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = x10 + x9
    x12 = x11 * x3
    x13 = x10 * x3
    x14 = 2.0 * x13
    x15 = x12 + x14
    x16 = x15 * x3
    x17 = -x2 - A[0]
    x18 = x15 * x17
    x19 = 3.0 * x10
    x20 = x0 * (x19 + 3.0 * x9)
    x21 = x0 * (x16 + 4.0 * x18 + 5.0 * x20)
    x22 = 8.0 * x13
    x23 = x0 * (4.0 * x12 + x22)
    x24 = x16 + x20
    x25 = x17 * x24
    x26 = x23 + x25
    x27 = x17 * x26
    x28 = x17 * x8
    x29 = x28 * x3
    x30 = x19 + 2.0 * x29
    x31 = x0 * (x30 + x9)
    x32 = x11 * x17
    x33 = x14 + x32
    x34 = x17 * x33
    x35 = 3.0 * x31 + 3.0 * x34
    x36 = x0 * (2.0 * x18 + 2.0 * x20 + x35)
    x37 = x0 * (x12 + x22 + 3.0 * x32)
    x38 = x18 + x20
    x39 = x17 * x38
    x40 = x37 + x39
    x41 = x17 * x40
    x42 = x21 + x27
    x43 = 2.0 * x0 * (x23 + x25 + 2.0 * x37 + 2.0 * x39) + x17 * x42
    x44 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x45 = da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**5.5)
    x46 = x44 * x45
    x47 = 6.89597470414309 * x46
    x48 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x49 = 0.564189583547756 * x1
    x50 = x48 * x49
    x51 = -x1 * (ax * A[1] + bx * B[1])
    x52 = -x51 - B[1]
    x53 = 18.2450341145548 * x52
    x54 = x0 * (x28 + x3 * x8)
    x55 = x10 + x29
    x56 = x17 * x55
    x57 = 2.0 * x0 * (2.0 * x13 + x32 + x54 + x56)
    x58 = x31 + x34
    x59 = x17 * x58
    x60 = x36 + x41
    x61 = x46 * x50
    x62 = x61 * (3.0 * x0 * (x37 + x39 + x57 + x59) + x17 * x60)
    x63 = -x1 * (ax * A[2] + bx * B[2])
    x64 = -x63 - B[2]
    x65 = 18.2450341145548 * x64
    x66 = x17**2 * x8
    x67 = x0 * (x30 + x66)
    x68 = x54 + x56
    x69 = x17 * x68
    x70 = x57 + x59
    x71 = x0 * (x35 + 2.0 * x67 + 2.0 * x69) + x17 * x70
    x72 = x46 * x71
    x73 = x48 * x7
    x74 = x52**2 * x73
    x75 = x0 * x73
    x76 = x74 + x75
    x77 = 23.5542377588857 * x76
    x78 = 0.318309886183791 * x6
    x79 = x77 * x78
    x80 = 40.7971365319473 * x52
    x81 = x50 * x80
    x82 = x44 * x7
    x83 = x64**2 * x82
    x84 = x0 * x82
    x85 = x83 + x84
    x86 = 23.5542377588857 * x78
    x87 = x45 * x48
    x88 = x86 * x87
    x89 = x10 + x66
    x90 = x17 * (2.0 * x10 + x89)
    x91 = x67 + x69
    x92 = x0 * (3.0 * x54 + 3.0 * x56 + x90) + x17 * x91
    x93 = 18.2450341145548 * x92
    x94 = x52 * x76
    x95 = x52 * x75
    x96 = 2.0 * x95
    x97 = x94 + x96
    x98 = x46 * x78
    x99 = x97 * x98
    x100 = 40.7971365319473 * x92
    x101 = x64 * x98
    x102 = x78 * x87
    x103 = x102 * x52
    x104 = x64 * x85
    x105 = x64 * x84
    x106 = 2.0 * x105
    x107 = x104 + x106
    x108 = x102 * x107
    x109 = 3.0 * x75
    x110 = x0 * (x109 + 3.0 * x74)
    x111 = x52 * x97
    x112 = x110 + x111
    x113 = x0 * (x19 + 3.0 * x66) + x17 * x90
    x114 = x47 * x78
    x115 = 0.179587122125167 * x45
    x116 = x115 * x85
    x117 = 3.0 * x84
    x118 = x0 * (x117 + 3.0 * x83)
    x119 = x107 * x64
    x120 = x118 + x119
    x121 = 6.89597470414309 * x120
    x122 = -x51 - A[1]
    x123 = 18.2450341145548 * x122
    x124 = x43 * x61
    x125 = x122 * x73
    x126 = x125 * x52
    x127 = x126 + x75
    x128 = 48.2718229290016 * x127
    x129 = 48.2718229290016 * x64
    x130 = x60 * x61
    x131 = 62.3186554316989 * x70
    x132 = x122 * x76
    x133 = x132 + x96
    x134 = x133 * x98
    x135 = 107.939077467081 * x127
    x136 = x102 * x131
    x137 = 48.2718229290016 * x91
    x138 = x122 * x97
    x139 = x110 + x138
    x140 = x139 * x98
    x141 = 107.939077467081 * x91
    x142 = 18.2450341145548 * x90
    x143 = 8.0 * x95
    x144 = x0 * (x143 + 4.0 * x94)
    x145 = x112 * x122
    x146 = x144 + x145
    x147 = x146 * x98
    x148 = 62.3186554316989 * x90
    x149 = x107 * x115
    x150 = x102 * x142
    x151 = -x63 - A[2]
    x152 = 18.2450341145548 * x151
    x153 = 48.2718229290016 * x151
    x154 = x151 * x82
    x155 = x154 * x64
    x156 = x155 + x84
    x157 = 48.2718229290016 * x102
    x158 = x151 * x98
    x159 = 107.939077467081 * x156
    x160 = x151 * x85
    x161 = x106 + x160
    x162 = 48.2718229290016 * x97
    x163 = x115 * x156
    x164 = x107 * x151
    x165 = x118 + x164
    x166 = x115 * x161
    x167 = 8.0 * x105
    x168 = x0 * (4.0 * x104 + x167)
    x169 = x120 * x151
    x170 = x168 + x169
    x171 = x122**2 * x73
    x172 = x171 + x75
    x173 = 23.5542377588857 * x172
    x174 = x0 * (x125 + x52 * x73)
    x175 = x122 * x127
    x176 = x174 + x175
    x177 = 62.3186554316989 * x98
    x178 = x177 * x40
    x179 = 80.4530382150027 * x58
    x180 = x109 + 2.0 * x126
    x181 = x0 * (x180 + x74)
    x182 = x122 * x133
    x183 = x181 + x182
    x184 = x183 * x98
    x185 = 139.348749811665 * x176
    x186 = 80.4530382150027 * x116
    x187 = x0 * (3.0 * x132 + x143 + x94)
    x188 = x122 * x139
    x189 = x187 + x188
    x190 = x177 * x189
    x191 = 139.348749811665 * x68
    x192 = 62.3186554316989 * x149
    x193 = x0 * (5.0 * x110 + x111 + 4.0 * x138)
    x194 = x122 * x146
    x195 = x193 + x194
    x196 = x195 * x46
    x197 = x115 * x89
    x198 = 40.7971365319473 * x151
    x199 = x102 * x122
    x200 = 139.348749811665 * x133
    x201 = 241.359114645008 * x163
    x202 = 139.348749811665 * x58
    x203 = 241.359114645008 * x166
    x204 = 107.939077467081 * x165
    x205 = 40.7971365319473 * x89
    x206 = 107.939077467081 * x163
    x207 = 139.348749811665 * x166
    x208 = x151**2 * x82
    x209 = x208 + x84
    x210 = 23.5542377588857 * x209
    x211 = 62.3186554316989 * x102
    x212 = x211 * x40
    x213 = x0 * (x154 + x64 * x82)
    x214 = x151 * x156
    x215 = x213 + x214
    x216 = x115 * x209
    x217 = 80.4530382150027 * x76
    x218 = x117 + 2.0 * x155
    x219 = x0 * (x218 + x83)
    x220 = x151 * x161
    x221 = x219 + x220
    x222 = x102 * x221
    x223 = 62.3186554316989 * x97
    x224 = x115 * x215
    x225 = x0 * (x104 + 3.0 * x160 + x167)
    x226 = x151 * x165
    x227 = x225 + x226
    x228 = x211 * x227
    x229 = x0 * (5.0 * x118 + x119 + 4.0 * x164)
    x230 = x151 * x170
    x231 = x229 + x230
    x232 = 18.2450341145548 * x26
    x233 = x122 * (x172 + 2.0 * x75)
    x234 = x233 * x98
    x235 = x0 * (x171 + x180)
    x236 = x122 * x176
    x237 = x235 + x236
    x238 = 48.2718229290016 * x98
    x239 = x237 * x238
    x240 = 48.2718229290016 * x233
    x241 = 62.3186554316989 * x33
    x242 = 2.0 * x0 * (x132 + x174 + x175 + 2.0 * x95)
    x243 = x122 * x183
    x244 = x242 + x243
    x245 = x244 * x98
    x246 = 107.939077467081 * x33
    x247 = 3.0 * x181 + 3.0 * x182
    x248 = x0 * (2.0 * x110 + 2.0 * x138 + x247)
    x249 = x122 * x189
    x250 = x248 + x249
    x251 = 107.939077467081 * x55
    x252 = 18.2450341145548 * x17
    x253 = 2.0 * x0 * (x144 + x145 + 2.0 * x187 + 2.0 * x188) + x122 * x195
    x254 = x49 * x5
    x255 = x254 * x46
    x256 = x253 * x255
    x257 = x250 * x255
    x258 = x45 * x5
    x259 = x258 * x85
    x260 = x259 * x78
    x261 = x258 * x78
    x262 = x107 * x261
    x263 = 48.2718229290016 * x17
    x264 = 18.2450341145548 * x261
    x265 = x233 * x264
    x266 = 40.7971365319473 * x172
    x267 = 107.939077467081 * x38
    x268 = 139.348749811665 * x33
    x269 = x115 * x55
    x270 = x17 * x261
    x271 = x183 * x261
    x272 = 139.348749811665 * x271
    x273 = 40.7971365319473 * x209
    x274 = 241.359114645008 * x224
    x275 = x258 * x50
    x276 = x122 * x275
    x277 = x151 * (x209 + 2.0 * x84)
    x278 = x102 * x277
    x279 = x157 * x38
    x280 = x0 * (x208 + x218)
    x281 = x151 * x215
    x282 = x280 + x281
    x283 = x115 * x76
    x284 = 2.0 * x0 * (2.0 * x105 + x160 + x213 + x214)
    x285 = x151 * x221
    x286 = x284 + x285
    x287 = x102 * x286
    x288 = 3.0 * x219 + 3.0 * x220
    x289 = x0 * (2.0 * x118 + 2.0 * x164 + x288)
    x290 = x151 * x227
    x291 = x289 + x290
    x292 = x264 * x277
    x293 = 62.3186554316989 * x261
    x294 = x286 * x293
    x295 = 2.0 * x0 * (x168 + x169 + 2.0 * x225 + 2.0 * x226) + x151 * x231
    x296 = x275 * x295
    x297 = x0 * (x109 + 3.0 * x171) + x122 * x233
    x298 = x0 * (3.0 * x174 + 3.0 * x175 + x233) + x122 * x237
    x299 = 18.2450341145548 * x298
    x300 = x15 * x98
    x301 = x0 * (2.0 * x235 + 2.0 * x236 + x247) + x122 * x244
    x302 = 23.5542377588857 * x11
    x303 = 40.7971365319473 * x298
    x304 = 18.2450341145548 * x3
    x305 = x255 * (3.0 * x0 * (x187 + x188 + x242 + x243) + x122 * x250)
    x306 = 62.3186554316989 * x11
    x307 = x261 * x3
    x308 = 107.939077467081 * x307
    x309 = 48.2718229290016 * x261
    x310 = 62.3186554316989 * x15
    x311 = 80.4530382150027 * x11
    x312 = x189 * x293
    x313 = x227 * x293
    x314 = x115 * x277
    x315 = x11 * x115
    x316 = x139 * x309
    x317 = x0 * (x117 + 3.0 * x208) + x151 * x277
    x318 = 6.89597470414309 * x317
    x319 = 18.2450341145548 * x15
    x320 = x0 * (3.0 * x213 + 3.0 * x214 + x277) + x151 * x282
    x321 = x102 * x320
    x322 = x0 * (2.0 * x280 + 2.0 * x281 + x288) + x151 * x286
    x323 = x264 * x97
    x324 = x258 * x322
    x325 = x275 * (3.0 * x0 * (x225 + x226 + x284 + x285) + x151 * x291)

    # 225 item(s)
    result[0, 0] = numpy.sum(
        x47 * x50 * (x0 * (3.0 * x21 + 3.0 * x27 + 4.0 * x36 + 4.0 * x41) + x17 * x43)
    )
    result[0, 1] = numpy.sum(x53 * x62)
    result[0, 2] = numpy.sum(x62 * x65)
    result[0, 3] = numpy.sum(x72 * x79)
    result[0, 4] = numpy.sum(x64 * x72 * x81)
    result[0, 5] = numpy.sum(x71 * x85 * x88)
    result[0, 6] = numpy.sum(x93 * x99)
    result[0, 7] = numpy.sum(x100 * x101 * x76)
    result[0, 8] = numpy.sum(x100 * x103 * x85)
    result[0, 9] = numpy.sum(x108 * x93)
    result[0, 10] = numpy.sum(x112 * x113 * x114)
    result[0, 11] = numpy.sum(x113 * x65 * x99)
    result[0, 12] = numpy.sum(x113 * x116 * x77)
    result[0, 13] = numpy.sum(x108 * x113 * x53)
    result[0, 14] = numpy.sum(x102 * x113 * x121)
    result[1, 0] = numpy.sum(x123 * x124)
    result[1, 1] = numpy.sum(x128 * x60 * x98)
    result[1, 2] = numpy.sum(x122 * x129 * x130)
    result[1, 3] = numpy.sum(x131 * x134)
    result[1, 4] = numpy.sum(x101 * x135 * x70)
    result[1, 5] = numpy.sum(x122 * x136 * x85)
    result[1, 6] = numpy.sum(x137 * x140)
    result[1, 7] = numpy.sum(x134 * x141 * x64)
    result[1, 8] = numpy.sum(x116 * x127 * x141)
    result[1, 9] = numpy.sum(x108 * x122 * x137)
    result[1, 10] = numpy.sum(x142 * x147)
    result[1, 11] = numpy.sum(x129 * x140 * x90)
    result[1, 12] = numpy.sum(x116 * x133 * x148)
    result[1, 13] = numpy.sum(x128 * x149 * x90)
    result[1, 14] = numpy.sum(x120 * x122 * x150)
    result[2, 0] = numpy.sum(x124 * x152)
    result[2, 1] = numpy.sum(x130 * x153 * x52)
    result[2, 2] = numpy.sum(x156 * x157 * x60)
    result[2, 3] = numpy.sum(x131 * x158 * x76)
    result[2, 4] = numpy.sum(x103 * x159 * x70)
    result[2, 5] = numpy.sum(x136 * x161)
    result[2, 6] = numpy.sum(x158 * x162 * x91)
    result[2, 7] = numpy.sum(x141 * x163 * x76)
    result[2, 8] = numpy.sum(x103 * x141 * x161)
    result[2, 9] = numpy.sum(x102 * x137 * x165)
    result[2, 10] = numpy.sum(x112 * x142 * x158)
    result[2, 11] = numpy.sum(x162 * x163 * x90)
    result[2, 12] = numpy.sum(x148 * x166 * x76)
    result[2, 13] = numpy.sum(x157 * x165 * x52 * x90)
    result[2, 14] = numpy.sum(x150 * x170)
    result[3, 0] = numpy.sum(x173 * x42 * x98)
    result[3, 1] = numpy.sum(x176 * x178)
    result[3, 2] = numpy.sum(x172 * x178 * x64)
    result[3, 3] = numpy.sum(x179 * x184)
    result[3, 4] = numpy.sum(x101 * x185 * x58)
    result[3, 5] = numpy.sum(x172 * x186 * x58)
    result[3, 6] = numpy.sum(x190 * x68)
    result[3, 7] = numpy.sum(x184 * x191 * x64)
    result[3, 8] = numpy.sum(x116 * x176 * x191)
    result[3, 9] = numpy.sum(x172 * x192 * x68)
    result[3, 10] = numpy.sum(x196 * x86 * x89)
    result[3, 11] = numpy.sum(x190 * x64 * x89)
    result[3, 12] = numpy.sum(x183 * x186 * x89)
    result[3, 13] = numpy.sum(x176 * x192 * x89)
    result[3, 14] = numpy.sum(x120 * x173 * x197)
    result[4, 0] = numpy.sum(x122 * x198 * x42 * x61)
    result[4, 1] = numpy.sum(x135 * x158 * x40)
    result[4, 2] = numpy.sum(x159 * x199 * x40)
    result[4, 3] = numpy.sum(x158 * x200 * x58)
    result[4, 4] = numpy.sum(x127 * x201 * x58)
    result[4, 5] = numpy.sum(x161 * x199 * x202)
    result[4, 6] = numpy.sum(107.939077467081 * x140 * x151 * x68)
    result[4, 7] = numpy.sum(x133 * x201 * x68)
    result[4, 8] = numpy.sum(x127 * x203 * x68)
    result[4, 9] = numpy.sum(x199 * x204 * x68)
    result[4, 10] = numpy.sum(x147 * x151 * x205)
    result[4, 11] = numpy.sum(x139 * x206 * x89)
    result[4, 12] = numpy.sum(x133 * x207 * x89)
    result[4, 13] = numpy.sum(x127 * x197 * x204)
    result[4, 14] = numpy.sum(x170 * x199 * x205)
    result[5, 0] = numpy.sum(x102 * x210 * x42)
    result[5, 1] = numpy.sum(x209 * x212 * x52)
    result[5, 2] = numpy.sum(x212 * x215)
    result[5, 3] = numpy.sum(x216 * x217 * x58)
    result[5, 4] = numpy.sum(x103 * x202 * x215)
    result[5, 5] = numpy.sum(x179 * x222)
    result[5, 6] = numpy.sum(x216 * x223 * x68)
    result[5, 7] = numpy.sum(x191 * x224 * x76)
    result[5, 8] = numpy.sum(x191 * x222 * x52)
    result[5, 9] = numpy.sum(x228 * x68)
    result[5, 10] = numpy.sum(x112 * x197 * x210)
    result[5, 11] = numpy.sum(x197 * x215 * x223)
    result[5, 12] = numpy.sum(x197 * x217 * x221)
    result[5, 13] = numpy.sum(x228 * x52 * x89)
    result[5, 14] = numpy.sum(x231 * x88 * x89)
    result[6, 0] = numpy.sum(x232 * x234)
    result[6, 1] = numpy.sum(x239 * x38)
    result[6, 2] = numpy.sum(x101 * x240 * x38)
    result[6, 3] = numpy.sum(x241 * x245)
    result[6, 4] = numpy.sum(x101 * x237 * x246)
    result[6, 5] = numpy.sum(x116 * x233 * x241)
    result[6, 6] = numpy.sum(x238 * x250 * x55)
    result[6, 7] = numpy.sum(x245 * x251 * x64)
    result[6, 8] = numpy.sum(x116 * x237 * x251)
    result[6, 9] = numpy.sum(x149 * x240 * x55)
    result[6, 10] = numpy.sum(x252 * x256)
    result[6, 11] = numpy.sum(x129 * x17 * x257)
    result[6, 12] = numpy.sum(62.3186554316989 * x17 * x244 * x260)
    result[6, 13] = numpy.sum(x237 * x262 * x263)
    result[6, 14] = numpy.sum(x120 * x17 * x265)
    result[7, 0] = numpy.sum(x158 * x26 * x266)
    result[7, 1] = numpy.sum(x158 * x176 * x267)
    result[7, 2] = numpy.sum(x172 * x206 * x38)
    result[7, 3] = numpy.sum(x151 * x184 * x268)
    result[7, 4] = numpy.sum(x176 * x201 * x33)
    result[7, 5] = numpy.sum(x172 * x207 * x33)
    result[7, 6] = numpy.sum(x158 * x189 * x251)
    result[7, 7] = numpy.sum(x183 * x201 * x55)
    result[7, 8] = numpy.sum(x176 * x203 * x55)
    result[7, 9] = numpy.sum(x172 * x204 * x269)
    result[7, 10] = numpy.sum(x17 * x196 * x198 * x254)
    result[7, 11] = numpy.sum(x159 * x189 * x270)
    result[7, 12] = numpy.sum(x161 * x17 * x272)
    result[7, 13] = numpy.sum(x176 * x204 * x270)
    result[7, 14] = numpy.sum(x170 * x266 * x270)
    result[8, 0] = numpy.sum(x199 * x26 * x273)
    result[8, 1] = numpy.sum(x135 * x216 * x38)
    result[8, 2] = numpy.sum(x199 * x215 * x267)
    result[8, 3] = numpy.sum(x200 * x216 * x33)
    result[8, 4] = numpy.sum(x127 * x274 * x33)
    result[8, 5] = numpy.sum(x122 * x222 * x268)
    result[8, 6] = numpy.sum(x139 * x216 * x251)
    result[8, 7] = numpy.sum(x133 * x274 * x55)
    result[8, 8] = numpy.sum(241.359114645008 * x127 * x221 * x269)
    result[8, 9] = numpy.sum(x199 * x227 * x251)
    result[8, 10] = numpy.sum(x146 * x270 * x273)
    result[8, 11] = numpy.sum(107.939077467081 * x139 * x215 * x270)
    result[8, 12] = numpy.sum(x200 * x221 * x270)
    result[8, 13] = numpy.sum(x135 * x227 * x270)
    result[8, 14] = numpy.sum(40.7971365319473 * x17 * x231 * x276)
    result[9, 0] = numpy.sum(x232 * x278)
    result[9, 1] = numpy.sum(x277 * x279 * x52)
    result[9, 2] = numpy.sum(x279 * x282)
    result[9, 3] = numpy.sum(x241 * x277 * x283)
    result[9, 4] = numpy.sum(x103 * x246 * x282)
    result[9, 5] = numpy.sum(x241 * x287)
    result[9, 6] = numpy.sum(x162 * x269 * x277)
    result[9, 7] = numpy.sum(x251 * x282 * x283)
    result[9, 8] = numpy.sum(x251 * x287 * x52)
    result[9, 9] = numpy.sum(x157 * x291 * x55)
    result[9, 10] = numpy.sum(x112 * x17 * x292)
    result[9, 11] = numpy.sum(x162 * x270 * x282)
    result[9, 12] = numpy.sum(x17 * x294 * x76)
    result[9, 13] = numpy.sum(x263 * x275 * x291 * x52)
    result[9, 14] = numpy.sum(x252 * x296)
    result[10, 0] = numpy.sum(x114 * x24 * x297)
    result[10, 1] = numpy.sum(x299 * x300)
    result[10, 2] = numpy.sum(x297 * x300 * x65)
    result[10, 3] = numpy.sum(x301 * x302 * x98)
    result[10, 4] = numpy.sum(x101 * x11 * x303)
    result[10, 5] = numpy.sum(x116 * x297 * x302)
    result[10, 6] = numpy.sum(x304 * x305)
    result[10, 7] = numpy.sum(40.7971365319473 * x255 * x3 * x301 * x64)
    result[10, 8] = numpy.sum(x260 * x3 * x303)
    result[10, 9] = numpy.sum(x262 * x297 * x304)
    result[10, 10] = numpy.sum(
        x254
        * x47
        * (x0 * (3.0 * x193 + 3.0 * x194 + 4.0 * x248 + 4.0 * x249) + x122 * x253)
    )
    result[10, 11] = numpy.sum(x305 * x65)
    result[10, 12] = numpy.sum(x259 * x301 * x86)
    result[10, 13] = numpy.sum(x262 * x299)
    result[10, 14] = numpy.sum(x121 * x261 * x297)
    result[11, 0] = numpy.sum(x152 * x234 * x24)
    result[11, 1] = numpy.sum(x15 * x151 * x239)
    result[11, 2] = numpy.sum(x15 * x163 * x240)
    result[11, 3] = numpy.sum(x151 * x245 * x306)
    result[11, 4] = numpy.sum(x11 * x206 * x237)
    result[11, 5] = numpy.sum(x166 * x233 * x306)
    result[11, 6] = numpy.sum(x153 * x257 * x3)
    result[11, 7] = numpy.sum(x159 * x244 * x307)
    result[11, 8] = numpy.sum(x161 * x237 * x308)
    result[11, 9] = numpy.sum(x165 * x240 * x307)
    result[11, 10] = numpy.sum(x152 * x256)
    result[11, 11] = numpy.sum(x156 * x250 * x309)
    result[11, 12] = numpy.sum(x161 * x244 * x293)
    result[11, 13] = numpy.sum(x165 * x237 * x309)
    result[11, 14] = numpy.sum(x170 * x265)
    result[12, 0] = numpy.sum(x173 * x216 * x24)
    result[12, 1] = numpy.sum(x176 * x216 * x310)
    result[12, 2] = numpy.sum(x172 * x224 * x310)
    result[12, 3] = numpy.sum(x183 * x216 * x311)
    result[12, 4] = numpy.sum(x11 * x185 * x224)
    result[12, 5] = numpy.sum(x115 * x172 * x221 * x311)
    result[12, 6] = numpy.sum(x209 * x3 * x312)
    result[12, 7] = numpy.sum(x215 * x272 * x3)
    result[12, 8] = numpy.sum(x185 * x221 * x307)
    result[12, 9] = numpy.sum(x172 * x3 * x313)
    result[12, 10] = numpy.sum(x195 * x210 * x261)
    result[12, 11] = numpy.sum(x215 * x312)
    result[12, 12] = numpy.sum(80.4530382150027 * x221 * x271)
    result[12, 13] = numpy.sum(x176 * x313)
    result[12, 14] = numpy.sum(x173 * x231 * x261)
    result[13, 0] = numpy.sum(x123 * x24 * x278)
    result[13, 1] = numpy.sum(x128 * x15 * x314)
    result[13, 2] = numpy.sum(x122 * x15 * x157 * x282)
    result[13, 3] = numpy.sum(x133 * x306 * x314)
    result[13, 4] = numpy.sum(x135 * x282 * x315)
    result[13, 5] = numpy.sum(x122 * x287 * x306)
    result[13, 6] = numpy.sum(x277 * x3 * x316)
    result[13, 7] = numpy.sum(x133 * x282 * x308)
    result[13, 8] = numpy.sum(x135 * x286 * x307)
    result[13, 9] = numpy.sum(48.2718229290016 * x276 * x291 * x3)
    result[13, 10] = numpy.sum(x146 * x292)
    result[13, 11] = numpy.sum(x282 * x316)
    result[13, 12] = numpy.sum(x133 * x294)
    result[13, 13] = numpy.sum(x128 * x261 * x291)
    result[13, 14] = numpy.sum(x123 * x296)
    result[14, 0] = numpy.sum(x102 * x24 * x318)
    result[14, 1] = numpy.sum(x103 * x317 * x319)
    result[14, 2] = numpy.sum(x319 * x321)
    result[14, 3] = numpy.sum(x315 * x317 * x77)
    result[14, 4] = numpy.sum(x11 * x321 * x80)
    result[14, 5] = numpy.sum(x102 * x302 * x322)
    result[14, 6] = numpy.sum(x3 * x317 * x323)
    result[14, 7] = numpy.sum(40.7971365319473 * x307 * x320 * x76)
    result[14, 8] = numpy.sum(x3 * x324 * x81)
    result[14, 9] = numpy.sum(x304 * x325)
    result[14, 10] = numpy.sum(x112 * x261 * x318)
    result[14, 11] = numpy.sum(x320 * x323)
    result[14, 12] = numpy.sum(x324 * x79)
    result[14, 13] = numpy.sum(x325 * x53)
    result[14, 14] = numpy.sum(
        6.89597470414309
        * x275
        * (x0 * (3.0 * x229 + 3.0 * x230 + 4.0 * x289 + 4.0 * x290) + x151 * x295)
    )
    return result


ovlp3d = {
    (0, 0): ovlp3d_00,
    (0, 1): ovlp3d_01,
    (0, 2): ovlp3d_02,
    (0, 3): ovlp3d_03,
    (0, 4): ovlp3d_04,
    (1, 0): ovlp3d_10,
    (1, 1): ovlp3d_11,
    (1, 2): ovlp3d_12,
    (1, 3): ovlp3d_13,
    (1, 4): ovlp3d_14,
    (2, 0): ovlp3d_20,
    (2, 1): ovlp3d_21,
    (2, 2): ovlp3d_22,
    (2, 3): ovlp3d_23,
    (2, 4): ovlp3d_24,
    (3, 0): ovlp3d_30,
    (3, 1): ovlp3d_31,
    (3, 2): ovlp3d_32,
    (3, 3): ovlp3d_33,
    (3, 4): ovlp3d_34,
    (4, 0): ovlp3d_40,
    (4, 1): ovlp3d_41,
    (4, 2): ovlp3d_42,
    (4, 3): ovlp3d_43,
    (4, 4): ovlp3d_44,
}
