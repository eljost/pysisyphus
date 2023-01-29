import numpy


_L_MAX = 4


def kinetic3d_00(ax, da, A, bx, db, B):
    """Cartesian 3D (ss) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 1), dtype=float)

    x0 = -ax
    x1 = (ax + bx) ** (-1.0)
    x2 = 0.5 / (ax + bx)
    x3 = 2.0 * ax**2
    x4 = ax * bx * x1
    x5 = (
        5.56832799683171
        * x1**1.5
        * numpy.exp(-x4 * ((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 + (A[2] - B[2]) ** 2))
    )

    # 1 item(s)
    result[0, 0] = numpy.sum(
        -0.507949087473928
        * da
        * db
        * x5
        * (
            3.0 * x0
            + x3 * (x2 + (-x1 * (ax * A[0] + bx * B[0]) + A[0]) ** 2)
            + x3 * (x2 + (-x1 * (ax * A[1] + bx * B[1]) + A[1]) ** 2)
            + x3 * (x2 + (-x1 * (ax * A[2] + bx * B[2]) + A[2]) ** 2)
        )
        * numpy.sqrt(ax**1.5)
        * numpy.sqrt(bx**1.5)
    )
    return result


def kinetic3d_01(ax, da, A, bx, db, B):
    """Cartesian 3D (sp) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 3), dtype=float)

    x0 = -ax
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = 2.0 * ax
    x4 = (2.0 * bx + x3) ** (-1.0)
    x5 = 2.0 * ax**2
    x6 = -x0 - x5 * (x4 + (x2 + A[0]) ** 2)
    x7 = numpy.sqrt(x1)
    x8 = 1.77245385090552 * x7
    x9 = bx * x1
    x10 = ax * x9
    x11 = numpy.exp(-x10 * (A[0] - B[0]) ** 2)
    x12 = -x11 * (x2 + B[0])
    x13 = x12 * x8
    x14 = x3 * x9
    x15 = numpy.exp(-x10 * (A[1] - B[1]) ** 2)
    x16 = numpy.exp(-x10 * (A[2] - B[2]) ** 2)
    x17 = 3.14159265358979 * x1 * x16
    x18 = 5.56832799683171 * x7
    x19 = x12 * x18
    x20 = -x1 * (ax * A[1] + bx * B[1])
    x21 = -x15 * (x0 + x5 * (x4 + (x20 + A[1]) ** 2))
    x22 = x1 * x16 * x21
    x23 = -x1 * (ax * A[2] + bx * B[2])
    x24 = -x16 * (x0 + x5 * (x4 + (x23 + A[2]) ** 2))
    x25 = x1 * x15
    x26 = 1.01589817494786 * da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**2.5)
    x27 = -x20 - B[1]
    x28 = x27 * x8
    x29 = x11 * x25
    x30 = x18 * x27 * x29
    x31 = x16 * x6
    x32 = -x23 - B[2]
    x33 = x32 * x8
    x34 = x18 * x32

    # 3 item(s)
    result[0, 0] = numpy.sum(
        x26 * (x13 * x15 * x17 * (x14 + x6) + x19 * x22 + x19 * x24 * x25)
    )
    result[0, 1] = numpy.sum(
        x26 * (x11 * x17 * x28 * (x14 * x15 + x21) + x24 * x30 + x30 * x31)
    )
    result[0, 2] = numpy.sum(
        x26
        * (
            x11 * x22 * x34
            + x29 * x31 * x34
            + 3.14159265358979 * x29 * x33 * (x14 * x16 + x24)
        )
    )
    return result


def kinetic3d_02(ax, da, A, bx, db, B):
    """Cartesian 3D (sd) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 6), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = -x0 * (ax * A[0] + bx * B[0])
    x2 = -x1 - B[0]
    x3 = -ax
    x4 = 2.0 * ax
    x5 = 2.0 * bx
    x6 = (x4 + x5) ** (-1.0)
    x7 = 2.0 * ax**2
    x8 = -x3 - x7 * (x6 + (x1 + A[0]) ** 2)
    x9 = ax * x0
    x10 = bx * x9
    x11 = numpy.exp(-x10 * (A[0] - B[0]) ** 2)
    x12 = numpy.sqrt(x0)
    x13 = 1.77245385090552 * x12
    x14 = x11 * x13
    x15 = x14 * x2
    x16 = bx * x0 * x4
    x17 = x15 * (x16 + x8)
    x18 = x14 * x6
    x19 = x14 * x2**2 + x18
    x20 = numpy.exp(-x10 * (A[1] - B[1]) ** 2)
    x21 = numpy.exp(-x10 * (A[2] - B[2]) ** 2)
    x22 = 3.14159265358979 * x0 * x21
    x23 = x20 * x22
    x24 = -x0 * (ax * A[1] + bx * B[1])
    x25 = -x3 - x7 * (x6 + (x24 + A[1]) ** 2)
    x26 = x19 * x23
    x27 = -x0 * (ax * A[2] + bx * B[2])
    x28 = -x3 - x7 * (x6 + (x27 + A[2]) ** 2)
    x29 = 0.179587122125167 * da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**3.5)
    x30 = 6.53197264742181 * x29
    x31 = -x24 - B[1]
    x32 = x13 * x20
    x33 = x31 * x32
    x34 = x33 * (x16 + x25)
    x35 = x11 * x22
    x36 = x34 * x35
    x37 = x17 * x23
    x38 = 5.56832799683171
    x39 = x0 * x11 * x20
    x40 = x12 * x2 * x21 * x38 * x39
    x41 = 11.3137084989848 * x29
    x42 = -x27 - B[2]
    x43 = x13 * x21
    x44 = x42 * x43
    x45 = x44 * (x16 + x28)
    x46 = 3.14159265358979 * x39
    x47 = x45 * x46
    x48 = x32 * x6
    x49 = x31**2 * x32 + x48
    x50 = x35 * x49
    x51 = x43 * x6
    x52 = x42**2 * x43 + x51
    x53 = x46 * x52

    # 6 item(s)
    result[0, 0] = numpy.sum(
        x30
        * (
            x23 * (x17 * x2 + x18 * x8 - x9 * (2.0 * x14 - x19 * x5))
            + x25 * x26
            + x26 * x28
        )
    )
    result[0, 1] = numpy.sum(x41 * (x2 * x36 + x28 * x31 * x40 + x31 * x37))
    result[0, 2] = numpy.sum(x41 * (x2 * x47 + x25 * x40 * x42 + x37 * x42))
    result[0, 3] = numpy.sum(
        x30
        * (
            x28 * x50
            + x35 * (x25 * x48 + x31 * x34 - x9 * (2.0 * x32 - x49 * x5))
            + x50 * x8
        )
    )
    result[0, 4] = numpy.sum(
        x41 * (x12 * x21 * x31 * x38 * x39 * x42 * x8 + x31 * x47 + x36 * x42)
    )
    result[0, 5] = numpy.sum(
        x30
        * (
            x25 * x53
            + x46 * (x28 * x51 + x42 * x45 - x9 * (2.0 * x43 - x5 * x52))
            + x53 * x8
        )
    )
    return result


def kinetic3d_03(ax, da, A, bx, db, B):
    """Cartesian 3D (sf) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 10), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - B[0]
    x6 = -ax
    x7 = 2.0 * ax**2
    x8 = -x6 - x7 * (x2 + (x4 + A[0]) ** 2)
    x9 = ax * x3
    x10 = bx * x9
    x11 = numpy.exp(-x10 * (A[0] - B[0]) ** 2)
    x12 = 1.77245385090552 * numpy.sqrt(x3)
    x13 = x11 * x12
    x14 = 2.0 * x13
    x15 = x13 * x5
    x16 = 4.0 * x10
    x17 = bx * x0 * x3
    x18 = x15 * (x17 + x8)
    x19 = x13 * x2
    x20 = x13 * x5**2 + x19
    x21 = x18 * x5 + x19 * x8 + x9 * (x1 * x20 - x14)
    x22 = x5 * (2.0 * x19 + x20)
    x23 = numpy.exp(-x10 * (A[1] - B[1]) ** 2)
    x24 = numpy.exp(-x10 * (A[2] - B[2]) ** 2)
    x25 = 3.14159265358979 * x24 * x3
    x26 = x23 * x25
    x27 = -x3 * (ax * A[1] + bx * B[1])
    x28 = -x6 - x7 * (x2 + (x27 + A[1]) ** 2)
    x29 = x22 * x26
    x30 = -x3 * (ax * A[2] + bx * B[2])
    x31 = -x6 - x7 * (x2 + (x30 + A[2]) ** 2)
    x32 = 0.179587122125167 * da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**4.5)
    x33 = 5.84237394672177 * x32
    x34 = -x27 - B[1]
    x35 = x12 * x23
    x36 = x34 * x35
    x37 = x36 * (x17 + x28)
    x38 = x12 * x24
    x39 = x21 * x26
    x40 = x20 * x26
    x41 = 13.0639452948436 * x32
    x42 = -x30 - B[2]
    x43 = x38 * x42
    x44 = x43 * (x17 + x31)
    x45 = x2 * x35
    x46 = x34**2 * x35 + x45
    x47 = 2.0 * x35
    x48 = x28 * x45 + x34 * x37 + x9 * (x1 * x46 - x47)
    x49 = x11 * x25
    x50 = x48 * x49
    x51 = x49 * x5
    x52 = 3.14159265358979 * x11 * x23 * x3
    x53 = x5 * x52
    x54 = x2 * x38
    x55 = x38 * x42**2 + x54
    x56 = 2.0 * x38
    x57 = x31 * x54 + x42 * x44 + x9 * (x1 * x55 - x56)
    x58 = x52 * x57
    x59 = x34 * (2.0 * x45 + x46)
    x60 = x49 * x59
    x61 = x42 * (2.0 * x54 + x55)
    x62 = x52 * x61

    # 10 item(s)
    result[0, 0] = numpy.sum(
        x33
        * (
            x26
            * (x2 * (x14 * x5 * x8 + x15 * x16) + x21 * x5 + x9 * (x1 * x22 - 3.0 * x15))
            + x28 * x29
            + x29 * x31
        )
    )
    result[0, 1] = numpy.sum(x41 * (x20 * x37 * x38 + x31 * x34 * x40 + x34 * x39))
    result[0, 2] = numpy.sum(x41 * (x20 * x35 * x44 + x28 * x40 * x42 + x39 * x42))
    result[0, 3] = numpy.sum(x41 * (x18 * x38 * x46 + x31 * x46 * x51 + x5 * x50))
    result[0, 4] = numpy.sum(
        22.6274169979695
        * x32
        * (x18 * x26 * x34 * x42 + x34 * x44 * x53 + x37 * x42 * x51)
    )
    result[0, 5] = numpy.sum(x41 * (x18 * x35 * x55 + x28 * x53 * x55 + x5 * x58))
    result[0, 6] = numpy.sum(
        x33
        * (
            x31 * x60
            + x49
            * (
                x2 * (x16 * x36 + x28 * x34 * x47)
                + x34 * x48
                + x9 * (x1 * x59 - 3.0 * x36)
            )
            + x60 * x8
        )
    )
    result[0, 7] = numpy.sum(x41 * (x13 * x44 * x46 + x42 * x46 * x49 * x8 + x42 * x50))
    result[0, 8] = numpy.sum(x41 * (x13 * x37 * x55 + x34 * x52 * x55 * x8 + x34 * x58))
    result[0, 9] = numpy.sum(
        x33
        * (
            x28 * x62
            + x52
            * (
                x2 * (x16 * x43 + x31 * x42 * x56)
                + x42 * x57
                + x9 * (x1 * x61 - 3.0 * x43)
            )
            + x62 * x8
        )
    )
    return result


def kinetic3d_04(ax, da, A, bx, db, B):
    """Cartesian 3D (sg) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((1, 15), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - B[0]
    x6 = -ax
    x7 = 2.0 * ax**2
    x8 = -x6 - x7 * (x2 + (x4 + A[0]) ** 2)
    x9 = ax * x3
    x10 = bx * x9
    x11 = numpy.exp(-x10 * (A[0] - B[0]) ** 2)
    x12 = 1.77245385090552 * numpy.sqrt(x3)
    x13 = x11 * x12
    x14 = x13 * x5
    x15 = bx * x0 * x3
    x16 = x14 * (x15 + x8)
    x17 = x16 * x5
    x18 = x13 * x5**2
    x19 = x13 * x2
    x20 = x18 + x19
    x21 = 2.0 * x13
    x22 = x9 * (x1 * x20 - x21)
    x23 = x19 * x8
    x24 = 4.0 * x10
    x25 = x17 + x22 + x23
    x26 = x5 * (2.0 * x19 + x20)
    x27 = x2 * (x14 * x24 + x21 * x5 * x8) + x25 * x5 + x9 * (x1 * x26 - 3.0 * x14)
    x28 = 3.0 * x2 * (x18 + x19) + x26 * x5
    x29 = numpy.exp(-x10 * (A[1] - B[1]) ** 2)
    x30 = numpy.exp(-x10 * (A[2] - B[2]) ** 2)
    x31 = 3.14159265358979 * x3 * x30
    x32 = x29 * x31
    x33 = -x3 * (ax * A[1] + bx * B[1])
    x34 = -x6 - x7 * (x2 + (x33 + A[1]) ** 2)
    x35 = x28 * x32
    x36 = -x3 * (ax * A[2] + bx * B[2])
    x37 = -x6 - x7 * (x2 + (x36 + A[2]) ** 2)
    x38 = 0.179587122125167 * da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**5.5)
    x39 = 4.41641957979107 * x38
    x40 = -x33 - B[1]
    x41 = x12 * x29
    x42 = x40 * x41
    x43 = x42 * (x15 + x34)
    x44 = x12 * x30
    x45 = x43 * x44
    x46 = x27 * x32
    x47 = x26 * x32
    x48 = 11.6847478934435 * x38
    x49 = -x36 - B[2]
    x50 = x44 * x49
    x51 = x50 * (x15 + x37)
    x52 = x40 * x43
    x53 = x40**2 * x41
    x54 = x2 * x41
    x55 = x53 + x54
    x56 = 2.0 * x41
    x57 = x9 * (x1 * x55 - x56)
    x58 = x34 * x54
    x59 = x52 + x57 + x58
    x60 = x20 * x44
    x61 = 15.084944665313 * x38
    x62 = 26.1278905896872 * x38
    x63 = x49 * x51
    x64 = x44 * x49**2
    x65 = x2 * x44
    x66 = x64 + x65
    x67 = 2.0 * x44
    x68 = x9 * (x1 * x66 - x67)
    x69 = x37 * x65
    x70 = x63 + x68 + x69
    x71 = x20 * x41
    x72 = x40 * (2.0 * x54 + x55)
    x73 = x2 * (x24 * x42 + x34 * x40 * x56) + x40 * x59 + x9 * (x1 * x72 - 3.0 * x42)
    x74 = x11 * x31
    x75 = x73 * x74
    x76 = x5 * x74
    x77 = 3.14159265358979 * x11 * x29 * x3
    x78 = x5 * x77
    x79 = x49 * (2.0 * x65 + x66)
    x80 = x2 * (x24 * x50 + x37 * x49 * x67) + x49 * x70 + x9 * (x1 * x79 - 3.0 * x50)
    x81 = x77 * x80
    x82 = 3.0 * x2 * (x53 + x54) + x40 * x72
    x83 = x74 * x82
    x84 = x13 * x55
    x85 = 3.0 * x2 * (x64 + x65) + x49 * x79
    x86 = x77 * x85

    # 15 item(s)
    result[0, 0] = numpy.sum(
        x39
        * (
            x32
            * (
                3.0 * x2 * (x17 + x22 + x23)
                + x27 * x5
                - 2.0 * x9 * (-bx * x28 + 2.0 * x18 + 2.0 * x19)
            )
            + x34 * x35
            + x35 * x37
        )
    )
    result[0, 1] = numpy.sum(x48 * (x26 * x45 + x37 * x40 * x47 + x40 * x46))
    result[0, 2] = numpy.sum(x48 * (x26 * x41 * x51 + x34 * x47 * x49 + x46 * x49))
    result[0, 3] = numpy.sum(x61 * (x25 * x44 * x55 + x37 * x55 * x60 + x59 * x60))
    result[0, 4] = numpy.sum(
        x62 * (x20 * x42 * x51 + x20 * x45 * x49 + x25 * x32 * x40 * x49)
    )
    result[0, 5] = numpy.sum(x61 * (x25 * x41 * x66 + x34 * x66 * x71 + x70 * x71))
    result[0, 6] = numpy.sum(x48 * (x16 * x44 * x72 + x37 * x72 * x76 + x5 * x75))
    result[0, 7] = numpy.sum(x62 * (x14 * x51 * x55 + x16 * x50 * x55 + x49 * x59 * x76))
    result[0, 8] = numpy.sum(x62 * (x14 * x43 * x66 + x16 * x42 * x66 + x40 * x70 * x78))
    result[0, 9] = numpy.sum(x48 * (x16 * x41 * x79 + x34 * x78 * x79 + x5 * x81))
    result[0, 10] = numpy.sum(
        x39
        * (
            x37 * x83
            + x74
            * (
                3.0 * x2 * (x52 + x57 + x58)
                + x40 * x73
                - 2.0 * x9 * (-bx * x82 + 2.0 * x53 + 2.0 * x54)
            )
            + x8 * x83
        )
    )
    result[0, 11] = numpy.sum(x48 * (x13 * x51 * x72 + x49 * x72 * x74 * x8 + x49 * x75))
    result[0, 12] = numpy.sum(x61 * (x13 * x59 * x66 + x66 * x8 * x84 + x70 * x84))
    result[0, 13] = numpy.sum(x48 * (x13 * x43 * x79 + x40 * x77 * x79 * x8 + x40 * x81))
    result[0, 14] = numpy.sum(
        x39
        * (
            x34 * x86
            + x77
            * (
                3.0 * x2 * (x63 + x68 + x69)
                + x49 * x80
                - 2.0 * x9 * (-bx * x85 + 2.0 * x64 + 2.0 * x65)
            )
            + x8 * x86
        )
    )
    return result


def kinetic3d_10(ax, da, A, bx, db, B):
    """Cartesian 3D (ps) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 1), dtype=float)

    x0 = -ax
    x1 = (ax + bx) ** (-1.0)
    x2 = x1 * (ax * A[0] + bx * B[0]) - A[0]
    x3 = 2.0 * ax
    x4 = (2.0 * bx + x3) ** (-1.0)
    x5 = 2.0 * ax**2
    x6 = -x0 - x5 * (x2**2 + x4)
    x7 = numpy.sqrt(x1)
    x8 = 1.77245385090552 * x7
    x9 = bx * x1
    x10 = ax * x9
    x11 = numpy.exp(-x10 * (A[0] - B[0]) ** 2)
    x12 = x11 * x2
    x13 = x12 * x8
    x14 = x3 * x9
    x15 = numpy.exp(-x10 * (A[1] - B[1]) ** 2)
    x16 = numpy.exp(-x10 * (A[2] - B[2]) ** 2)
    x17 = 3.14159265358979 * x1 * x16
    x18 = 5.56832799683171 * x7
    x19 = x12 * x18
    x20 = x1 * (ax * A[1] + bx * B[1]) - A[1]
    x21 = -x15 * (x0 + x5 * (x20**2 + x4))
    x22 = x1 * x16 * x21
    x23 = x1 * (ax * A[2] + bx * B[2]) - A[2]
    x24 = -x16 * (x0 + x5 * (x23**2 + x4))
    x25 = x1 * x15
    x26 = 1.01589817494786 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**1.5)
    x27 = x20 * x8
    x28 = x11 * x25
    x29 = x18 * x20 * x28
    x30 = x16 * x6
    x31 = x23 * x8
    x32 = x18 * x23

    # 3 item(s)
    result[0, 0] = numpy.sum(
        x26 * (x13 * x15 * x17 * (x14 + x6) + x19 * x22 + x19 * x24 * x25)
    )
    result[1, 0] = numpy.sum(
        x26 * (x11 * x17 * x27 * (x14 * x15 + x21) + x24 * x29 + x29 * x30)
    )
    result[2, 0] = numpy.sum(
        x26
        * (
            x11 * x22 * x32
            + x28 * x30 * x32
            + 3.14159265358979 * x28 * x31 * (x14 * x16 + x24)
        )
    )
    return result


def kinetic3d_11(ax, da, A, bx, db, B):
    """Cartesian 3D (pp) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 3), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = -x0 * (ax * A[0] + bx * B[0])
    x2 = -x1 - A[0]
    x3 = -ax
    x4 = 2.0 * ax
    x5 = (2.0 * bx + x4) ** (-1.0)
    x6 = 2.0 * ax**2
    x7 = -x3 - x6 * (x2**2 + x5)
    x8 = numpy.sqrt(x0)
    x9 = 1.77245385090552 * x8
    x10 = bx * x0
    x11 = ax * x10
    x12 = numpy.exp(-x11 * (A[0] - B[0]) ** 2)
    x13 = -x12 * (x1 + B[0])
    x14 = x13 * x9
    x15 = x10 * x4
    x16 = x14 * (x15 + x7)
    x17 = x5 * x9
    x18 = x12 * x17
    x19 = x14 * x2 + x18
    x20 = numpy.exp(-x11 * (A[1] - B[1]) ** 2)
    x21 = numpy.exp(-x11 * (A[2] - B[2]) ** 2)
    x22 = 3.14159265358979 * x0 * x21
    x23 = x20 * x22
    x24 = -x0 * (ax * A[1] + bx * B[1])
    x25 = -x24 - A[1]
    x26 = -x3 - x6 * (x25**2 + x5)
    x27 = x19 * x23
    x28 = -x0 * (ax * A[2] + bx * B[2])
    x29 = -x28 - A[2]
    x30 = -x3 - x6 * (x29**2 + x5)
    x31 = 2.03179634989571 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**2.5)
    x32 = -x24 - B[1]
    x33 = x20 * x32
    x34 = x26 * x9
    x35 = x33 * x9
    x36 = x15 * x35 + x33 * x34
    x37 = x12 * x2
    x38 = x37 * x9
    x39 = x23 * x38 * (x15 + x7)
    x40 = x0 * x33
    x41 = 5.56832799683171 * x21 * x8
    x42 = x37 * x41
    x43 = -x28 - B[2]
    x44 = x21 * x9
    x45 = x30 * x44
    x46 = x43 * x44
    x47 = x15 * x46 + x43 * x45
    x48 = x0 * x20
    x49 = 3.14159265358979 * x48
    x50 = x26 * x48
    x51 = x16 * x23
    x52 = x20 * x25
    x53 = x52 * (x15 * x9 + x34)
    x54 = x0 * x52
    x55 = x41 * x54
    x56 = x17 * x20
    x57 = x25 * x35 + x56
    x58 = x12 * x22
    x59 = x57 * x58
    x60 = 3.14159265358979 * x12
    x61 = x12 * x7
    x62 = x29 * (x15 * x44 + x45)
    x63 = x29 * x41
    x64 = x17 * x21
    x65 = x29 * x46 + x64
    x66 = x12 * x49
    x67 = x65 * x66

    # 9 item(s)
    result[0, 0] = numpy.sum(
        x31 * (x23 * (x15 * x19 + x16 * x2 + x18 * x7) + x26 * x27 + x27 * x30)
    )
    result[0, 1] = numpy.sum(x31 * (x22 * x36 * x37 + x30 * x40 * x42 + x32 * x39))
    result[0, 2] = numpy.sum(x31 * (x37 * x47 * x49 + x39 * x43 + x42 * x43 * x50))
    result[1, 0] = numpy.sum(x31 * (x13 * x22 * x53 + x13 * x30 * x55 + x25 * x51))
    result[1, 1] = numpy.sum(
        x31 * (x30 * x59 + x58 * (x15 * x57 + x25 * x36 + x26 * x56) + x59 * x7)
    )
    result[1, 2] = numpy.sum(x31 * (x43 * x53 * x58 + x43 * x55 * x61 + x47 * x54 * x60))
    result[2, 0] = numpy.sum(x31 * (x13 * x49 * x62 + x13 * x50 * x63 + x29 * x51))
    result[2, 1] = numpy.sum(x31 * (x29 * x36 * x58 + x40 * x60 * x62 + x40 * x61 * x63))
    result[2, 2] = numpy.sum(
        x31 * (x26 * x67 + x66 * (x15 * x65 + x29 * x47 + x30 * x64) + x67 * x7)
    )
    return result


def kinetic3d_12(ax, da, A, bx, db, B):
    """Cartesian 3D (pd) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 6), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - B[0]
    x6 = -ax
    x7 = -x4 - A[0]
    x8 = 2.0 * ax**2
    x9 = -x6 - x8 * (x2 + x7**2)
    x10 = ax * x3
    x11 = bx * x10
    x12 = numpy.exp(-x11 * (A[0] - B[0]) ** 2)
    x13 = 1.77245385090552 * numpy.sqrt(x3)
    x14 = x12 * x13
    x15 = 2.0 * x14
    x16 = x14 * x5
    x17 = 4.0 * x11
    x18 = bx * x0
    x19 = x18 * x3
    x20 = x16 * (x19 + x9)
    x21 = x14 * x2
    x22 = x14 * x5**2 + x21
    x23 = x21 * x9
    x24 = x10 * (x1 * x22 - x15) + x20 * x5 + x23
    x25 = x3 * (2.0 * x21 * x5 + x22 * x7)
    x26 = numpy.exp(-x11 * (A[1] - B[1]) ** 2)
    x27 = numpy.exp(-x11 * (A[2] - B[2]) ** 2)
    x28 = 3.14159265358979 * x27 * x3
    x29 = x26 * x28
    x30 = -x3 * (ax * A[1] + bx * B[1])
    x31 = -x30 - A[1]
    x32 = -x6 - x8 * (x2 + x31**2)
    x33 = 3.14159265358979 * x26
    x34 = x25 * x27 * x33
    x35 = -x3 * (ax * A[2] + bx * B[2])
    x36 = -x35 - A[2]
    x37 = -x6 - x8 * (x2 + x36**2)
    x38 = 0.179587122125167 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**3.5)
    x39 = 13.0639452948436 * x38
    x40 = -x30 - B[1]
    x41 = x13 * x26
    x42 = x32 * x41
    x43 = x40 * x41
    x44 = x19 * x43 + x40 * x42
    x45 = x16 * x7 + x21
    x46 = x13 * x27
    x47 = x29 * (x19 * x45 + x20 * x7 + x23)
    x48 = x29 * x45
    x49 = 22.6274169979695 * x38
    x50 = -x35 - B[2]
    x51 = x37 * x46
    x52 = x46 * x50
    x53 = x19 * x52 + x50 * x51
    x54 = x2 * x41
    x55 = x40**2 * x41 + x54
    x56 = x14 * x7
    x57 = x56 * (x19 + x9)
    x58 = 2.0 * x41
    x59 = x32 * x54
    x60 = x10 * (x1 * x55 - x58) + x40 * x44 + x59
    x61 = x12 * x28
    x62 = x60 * x61
    x63 = x61 * x7
    x64 = x12 * x3 * x33
    x65 = x64 * x7
    x66 = x29 * x50
    x67 = x2 * x46
    x68 = x46 * x50**2 + x67
    x69 = 2.0 * x46
    x70 = x37 * x67
    x71 = x10 * (x1 * x68 - x69) + x50 * x53 + x70
    x72 = x64 * x71
    x73 = x31 * (x19 * x41 + x42)
    x74 = x24 * x29
    x75 = x22 * x29
    x76 = x31 * x43 + x54
    x77 = x61 * (x19 * x76 + x31 * x44 + x59)
    x78 = x5 * x61
    x79 = x31 * x64
    x80 = x31 * x55 + 2.0 * x40 * x54
    x81 = x61 * x80
    x82 = x61 * x9
    x83 = x36 * (x19 * x46 + x51)
    x84 = x5 * x64
    x85 = x36 * x52 + x67
    x86 = x64 * (x19 * x85 + x36 * x53 + x70)
    x87 = x36 * x68 + 2.0 * x50 * x67
    x88 = x64 * x87

    # 18 item(s)
    result[0, 0] = numpy.sum(
        x39
        * (
            x29 * (x18 * x25 + x2 * (x15 * x5 * x9 + x16 * x17) + x24 * x7)
            + x32 * x34
            + x34 * x37
        )
    )
    result[0, 1] = numpy.sum(x49 * (x37 * x40 * x48 + x40 * x47 + x44 * x45 * x46))
    result[0, 2] = numpy.sum(x49 * (x32 * x48 * x50 + x41 * x45 * x53 + x47 * x50))
    result[0, 3] = numpy.sum(x39 * (x37 * x55 * x63 + x46 * x55 * x57 + x62 * x7))
    result[0, 4] = numpy.sum(x49 * (x40 * x53 * x65 + x40 * x57 * x66 + x44 * x50 * x63))
    result[0, 5] = numpy.sum(x39 * (x32 * x65 * x68 + x41 * x57 * x68 + x7 * x72))
    result[1, 0] = numpy.sum(x39 * (x22 * x46 * x73 + x31 * x37 * x75 + x31 * x74))
    result[1, 1] = numpy.sum(x49 * (x20 * x46 * x76 + x37 * x76 * x78 + x5 * x77))
    result[1, 2] = numpy.sum(x49 * (x20 * x31 * x66 + x5 * x53 * x79 + x50 * x73 * x78))
    result[1, 3] = numpy.sum(
        x39
        * (
            x37 * x81
            + x61 * (x19 * x80 + x2 * (x17 * x43 + x32 * x40 * x58) + x31 * x60)
            + x81 * x9
        )
    )
    result[1, 4] = numpy.sum(x49 * (x14 * x53 * x76 + x50 * x76 * x82 + x50 * x77))
    result[1, 5] = numpy.sum(x39 * (x14 * x68 * x73 + x31 * x72 + x68 * x79 * x9))
    result[2, 0] = numpy.sum(x39 * (x22 * x41 * x83 + x32 * x36 * x75 + x36 * x74))
    result[2, 1] = numpy.sum(
        x49 * (x20 * x29 * x36 * x40 + x36 * x44 * x78 + x40 * x83 * x84)
    )
    result[2, 2] = numpy.sum(x49 * (x20 * x41 * x85 + x32 * x84 * x85 + x5 * x86))
    result[2, 3] = numpy.sum(x39 * (x14 * x55 * x83 + x36 * x55 * x82 + x36 * x62))
    result[2, 4] = numpy.sum(x49 * (x14 * x44 * x85 + x40 * x64 * x85 * x9 + x40 * x86))
    result[2, 5] = numpy.sum(
        x39
        * (
            x32 * x88
            + x64 * (x19 * x87 + x2 * (x17 * x52 + x37 * x50 * x69) + x36 * x71)
            + x88 * x9
        )
    )
    return result


def kinetic3d_13(ax, da, A, bx, db, B):
    """Cartesian 3D (pf) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 10), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - B[0]
    x6 = -ax
    x7 = -x4 - A[0]
    x8 = 2.0 * ax**2
    x9 = -x6 - x8 * (x2 + x7**2)
    x10 = ax * x3
    x11 = bx * x10
    x12 = numpy.exp(-x11 * (A[0] - B[0]) ** 2)
    x13 = 1.77245385090552 * numpy.sqrt(x3)
    x14 = x12 * x13
    x15 = x14 * x5
    x16 = bx * x0
    x17 = x16 * x3
    x18 = x15 * (x17 + x9)
    x19 = x18 * x5
    x20 = x14 * x5**2
    x21 = x14 * x2
    x22 = x20 + x21
    x23 = 2.0 * x14
    x24 = x10 * (x1 * x22 - x23)
    x25 = x21 * x9
    x26 = 4.0 * x11
    x27 = x2 * (x15 * x26 + x23 * x5 * x9)
    x28 = x19 + x24 + x25
    x29 = 2.0 * x21 * x5
    x30 = x22 * x5 + x29
    x31 = x10 * (x1 * x30 - 3.0 * x15) + x27 + x28 * x5
    x32 = x3 * (3.0 * x2 * (x20 + x21) + x30 * x7)
    x33 = numpy.exp(-x11 * (A[1] - B[1]) ** 2)
    x34 = numpy.exp(-x11 * (A[2] - B[2]) ** 2)
    x35 = 3.14159265358979 * x3 * x34
    x36 = x33 * x35
    x37 = -x3 * (ax * A[1] + bx * B[1])
    x38 = -x37 - A[1]
    x39 = -x6 - x8 * (x2 + x38**2)
    x40 = 3.14159265358979 * x33
    x41 = x32 * x34 * x40
    x42 = -x3 * (ax * A[2] + bx * B[2])
    x43 = -x42 - A[2]
    x44 = -x6 - x8 * (x2 + x43**2)
    x45 = 0.179587122125167 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**4.5)
    x46 = 11.6847478934435 * x45
    x47 = x22 * x7 + x29
    x48 = -x37 - B[1]
    x49 = x13 * x33
    x50 = x48 * x49
    x51 = x50 * (x17 + x39)
    x52 = x13 * x34
    x53 = x51 * x52
    x54 = x36 * (x17 * x47 + x27 + x28 * x7)
    x55 = x36 * x47
    x56 = 26.1278905896872 * x45
    x57 = -x42 - B[2]
    x58 = x52 * x57
    x59 = x58 * (x17 + x44)
    x60 = x49 * x59
    x61 = x15 * x7 + x21
    x62 = x17 * x61 + x18 * x7 + x25
    x63 = x48**2 * x49
    x64 = x2 * x49
    x65 = x63 + x64
    x66 = x52 * x65
    x67 = x48 * x51
    x68 = 2.0 * x49
    x69 = x10 * (x1 * x65 - x68)
    x70 = x39 * x64
    x71 = x67 + x69 + x70
    x72 = x36 * x57
    x73 = 45.2548339959391 * x45
    x74 = x52 * x57**2
    x75 = x2 * x52
    x76 = x74 + x75
    x77 = x49 * x76
    x78 = x57 * x59
    x79 = 2.0 * x52
    x80 = x10 * (x1 * x76 - x79)
    x81 = x44 * x75
    x82 = x78 + x80 + x81
    x83 = 2.0 * x48 * x64
    x84 = x48 * x65 + x83
    x85 = x14 * x7
    x86 = x85 * (x17 + x9)
    x87 = x2 * (x26 * x50 + x39 * x48 * x68)
    x88 = x10 * (x1 * x84 - 3.0 * x50) + x48 * x71 + x87
    x89 = x12 * x35
    x90 = x88 * x89
    x91 = x7 * x89
    x92 = x12 * x3 * x40
    x93 = x7 * x92
    x94 = 2.0 * x57 * x75
    x95 = x57 * x76 + x94
    x96 = x2 * (x26 * x58 + x44 * x57 * x79)
    x97 = x10 * (x1 * x95 - 3.0 * x58) + x57 * x82 + x96
    x98 = x92 * x97
    x99 = x39 * x49
    x100 = x38 * (x17 * x49 + x99)
    x101 = x31 * x36
    x102 = x30 * x36
    x103 = x38 * x50 + x64
    x104 = x103 * x17 + x38 * x51 + x70
    x105 = x22 * x52
    x106 = x38 * x65 + x83
    x107 = x89 * (x106 * x17 + x38 * x71 + x87)
    x108 = x5 * x89
    x109 = x38 * x92
    x110 = 3.0 * x2 * (x63 + x64) + x38 * x84
    x111 = x110 * x89
    x112 = x89 * x9
    x113 = x14 * x76
    x114 = x43 * x52
    x115 = x114 * (x17 + x44)
    x116 = x43 * x58 + x75
    x117 = x116 * x17 + x43 * x59 + x81
    x118 = x5 * x92
    x119 = x43 * x76 + x94
    x120 = x92 * (x119 * x17 + x43 * x82 + x96)
    x121 = x14 * x65
    x122 = 3.0 * x2 * (x74 + x75) + x43 * x95
    x123 = x122 * x92

    # 30 item(s)
    result[0, 0] = numpy.sum(
        x46
        * (
            x36 * (x16 * x32 + 3.0 * x2 * (x19 + x24 + x25) + x31 * x7)
            + x39 * x41
            + x41 * x44
        )
    )
    result[0, 1] = numpy.sum(x56 * (x44 * x48 * x55 + x47 * x53 + x48 * x54))
    result[0, 2] = numpy.sum(x56 * (x39 * x55 * x57 + x47 * x60 + x54 * x57))
    result[0, 3] = numpy.sum(x56 * (x44 * x61 * x66 + x52 * x61 * x71 + x62 * x66))
    result[0, 4] = numpy.sum(x73 * (x48 * x62 * x72 + x50 * x59 * x61 + x53 * x57 * x61))
    result[0, 5] = numpy.sum(x56 * (x39 * x61 * x77 + x49 * x61 * x82 + x62 * x77))
    result[0, 6] = numpy.sum(x46 * (x44 * x84 * x91 + x52 * x84 * x86 + x7 * x90))
    result[0, 7] = numpy.sum(x56 * (x57 * x71 * x91 + x58 * x65 * x86 + x59 * x65 * x85))
    result[0, 8] = numpy.sum(x56 * (x48 * x82 * x93 + x50 * x76 * x86 + x51 * x76 * x85))
    result[0, 9] = numpy.sum(x46 * (x39 * x93 * x95 + x49 * x86 * x95 + x7 * x98))
    result[1, 0] = numpy.sum(x46 * (x100 * x30 * x52 + x101 * x38 + x102 * x38 * x44))
    result[1, 1] = numpy.sum(x56 * (x103 * x105 * x44 + x103 * x28 * x52 + x104 * x105))
    result[1, 2] = numpy.sum(x56 * (x100 * x22 * x58 + x22 * x38 * x60 + x28 * x38 * x72))
    result[1, 3] = numpy.sum(x56 * (x106 * x108 * x44 + x106 * x18 * x52 + x107 * x5))
    result[1, 4] = numpy.sum(
        x73 * (x103 * x15 * x59 + x103 * x18 * x58 + x104 * x108 * x57)
    )
    result[1, 5] = numpy.sum(x56 * (x100 * x15 * x76 + x109 * x5 * x82 + x18 * x38 * x77))
    result[1, 6] = numpy.sum(
        x46
        * (
            x111 * x44
            + x111 * x9
            + x89 * (x110 * x17 + 3.0 * x2 * (x67 + x69 + x70) + x38 * x88)
        )
    )
    result[1, 7] = numpy.sum(x56 * (x106 * x112 * x57 + x106 * x14 * x59 + x107 * x57))
    result[1, 8] = numpy.sum(x56 * (x103 * x113 * x9 + x103 * x14 * x82 + x104 * x113))
    result[1, 9] = numpy.sum(x46 * (x100 * x14 * x95 + x109 * x9 * x95 + x38 * x98))
    result[2, 0] = numpy.sum(x46 * (x101 * x43 + x102 * x39 * x43 + x115 * x30 * x49))
    result[2, 1] = numpy.sum(
        x56 * (x115 * x22 * x50 + x22 * x43 * x53 + x28 * x36 * x43 * x48)
    )
    result[2, 2] = numpy.sum(
        x56 * (x116 * x22 * x99 + x116 * x28 * x49 + x117 * x22 * x49)
    )
    result[2, 3] = numpy.sum(
        x56 * (x108 * x43 * x71 + x115 * x15 * x65 + x18 * x43 * x66)
    )
    result[2, 4] = numpy.sum(
        x73 * (x116 * x15 * x51 + x116 * x18 * x50 + x117 * x118 * x48)
    )
    result[2, 5] = numpy.sum(x56 * (x118 * x119 * x39 + x119 * x18 * x49 + x120 * x5))
    result[2, 6] = numpy.sum(x46 * (x112 * x43 * x84 + x115 * x14 * x84 + x43 * x90))
    result[2, 7] = numpy.sum(x56 * (x116 * x121 * x9 + x116 * x14 * x71 + x117 * x121))
    result[2, 8] = numpy.sum(
        x56 * (x119 * x14 * x51 + x119 * x48 * x9 * x92 + x120 * x48)
    )
    result[2, 9] = numpy.sum(
        x46
        * (
            x123 * x39
            + x123 * x9
            + x92 * (x122 * x17 + 3.0 * x2 * (x78 + x80 + x81) + x43 * x97)
        )
    )
    return result


def kinetic3d_14(ax, da, A, bx, db, B):
    """Cartesian 3D (pg) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 15), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - B[0]
    x6 = -ax
    x7 = -x4 - A[0]
    x8 = 2.0 * ax**2
    x9 = -x6 - x8 * (x2 + x7**2)
    x10 = ax * x3
    x11 = bx * x10
    x12 = numpy.exp(-x11 * (A[0] - B[0]) ** 2)
    x13 = 1.77245385090552 * numpy.sqrt(x3)
    x14 = x12 * x13
    x15 = 2.0 * x14
    x16 = x14 * x5
    x17 = 4.0 * x11
    x18 = x2 * (x15 * x5 * x9 + x16 * x17)
    x19 = bx * x0
    x20 = x19 * x3
    x21 = x16 * (x20 + x9)
    x22 = x21 * x5
    x23 = x14 * x5**2
    x24 = x14 * x2
    x25 = x23 + x24
    x26 = x10 * (x1 * x25 - x15)
    x27 = x24 * x9
    x28 = x22 + x26 + x27
    x29 = x28 * x5
    x30 = x25 * x5
    x31 = x24 * x5
    x32 = 2.0 * x31
    x33 = x30 + x32
    x34 = x10 * (x1 * x33 - 3.0 * x16)
    x35 = 3.0 * x2 * (x22 + x26 + x27)
    x36 = x18 + x29 + x34
    x37 = 3.0 * x2 * (x23 + x24)
    x38 = x33 * x5 + x37
    x39 = -2.0 * x10 * (-bx * x38 + 2.0 * x23 + 2.0 * x24) + x35 + x36 * x5
    x40 = x3 * (4.0 * x2 * (x30 + 2.0 * x31) + x38 * x7)
    x41 = numpy.exp(-x11 * (A[1] - B[1]) ** 2)
    x42 = numpy.exp(-x11 * (A[2] - B[2]) ** 2)
    x43 = 3.14159265358979 * x3 * x42
    x44 = x41 * x43
    x45 = -x3 * (ax * A[1] + bx * B[1])
    x46 = -x45 - A[1]
    x47 = -x6 - x8 * (x2 + x46**2)
    x48 = 3.14159265358979 * x41
    x49 = x40 * x42 * x48
    x50 = -x3 * (ax * A[2] + bx * B[2])
    x51 = -x50 - A[2]
    x52 = -x6 - x8 * (x2 + x51**2)
    x53 = 0.179587122125167 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**5.5)
    x54 = 8.83283915958214 * x53
    x55 = x33 * x7 + x37
    x56 = -x45 - B[1]
    x57 = x13 * x41
    x58 = x56 * x57
    x59 = x58 * (x20 + x47)
    x60 = x13 * x42
    x61 = x59 * x60
    x62 = x44 * (x20 * x55 + x35 + x36 * x7)
    x63 = x44 * x55
    x64 = 23.3694957868871 * x53
    x65 = -x50 - B[2]
    x66 = x60 * x65
    x67 = x66 * (x20 + x52)
    x68 = x57 * x67
    x69 = x56 * x59
    x70 = x56**2 * x57
    x71 = x2 * x57
    x72 = x70 + x71
    x73 = 2.0 * x57
    x74 = x10 * (x1 * x72 - x73)
    x75 = x47 * x71
    x76 = x69 + x74 + x75
    x77 = x25 * x7 + x32
    x78 = x60 * x77
    x79 = x18 + x20 * x77 + x28 * x7
    x80 = x60 * x72
    x81 = 30.169889330626 * x53
    x82 = x44 * x65
    x83 = 52.2557811793745 * x53
    x84 = x65 * x67
    x85 = x60 * x65**2
    x86 = x2 * x60
    x87 = x85 + x86
    x88 = 2.0 * x60
    x89 = x10 * (x1 * x87 - x88)
    x90 = x52 * x86
    x91 = x84 + x89 + x90
    x92 = x57 * x77
    x93 = x57 * x87
    x94 = x16 * x7 + x24
    x95 = x20 * x94 + x21 * x7 + x27
    x96 = x56 * x72
    x97 = x56 * x71
    x98 = 2.0 * x97
    x99 = x96 + x98
    x100 = x60 * x99
    x101 = x2 * (x17 * x58 + x47 * x56 * x73)
    x102 = x56 * x76
    x103 = x10 * (x1 * x99 - 3.0 * x58)
    x104 = x101 + x102 + x103
    x105 = x65 * x87
    x106 = x65 * x86
    x107 = 2.0 * x106
    x108 = x105 + x107
    x109 = x108 * x57
    x110 = x2 * (x17 * x66 + x52 * x65 * x88)
    x111 = x65 * x91
    x112 = x10 * (x1 * x108 - 3.0 * x66)
    x113 = x110 + x111 + x112
    x114 = 3.0 * x2 * (x70 + x71)
    x115 = x114 + x56 * x99
    x116 = x14 * x7
    x117 = x116 * (x20 + x9)
    x118 = 3.0 * x2 * (x69 + x74 + x75)
    x119 = -2.0 * x10 * (-bx * x115 + 2.0 * x70 + 2.0 * x71) + x104 * x56 + x118
    x120 = x12 * x43
    x121 = x119 * x120
    x122 = x120 * x7
    x123 = x12 * x3 * x48
    x124 = x123 * x7
    x125 = 3.0 * x2 * (x85 + x86)
    x126 = x108 * x65 + x125
    x127 = 3.0 * x2 * (x84 + x89 + x90)
    x128 = -2.0 * x10 * (-bx * x126 + 2.0 * x85 + 2.0 * x86) + x113 * x65 + x127
    x129 = x123 * x128
    x130 = x46 * x57
    x131 = x130 * (x20 + x47)
    x132 = x39 * x44
    x133 = x38 * x44
    x134 = x46 * x58 + x71
    x135 = x134 * x20 + x46 * x59 + x75
    x136 = x33 * x60
    x137 = x46 * x72 + x98
    x138 = x137 * x60
    x139 = x101 + x137 * x20 + x46 * x76
    x140 = x25 * x60
    x141 = x114 + x46 * x99
    x142 = x120 * (x104 * x46 + x118 + x141 * x20)
    x143 = x120 * x5
    x144 = x123 * x46
    x145 = x115 * x46 + 4.0 * x2 * (x96 + 2.0 * x97)
    x146 = x120 * x145
    x147 = x120 * x9
    x148 = x137 * x14
    x149 = x108 * x14
    x150 = x51 * x60
    x151 = x150 * (x20 + x52)
    x152 = x51 * x66 + x86
    x153 = x152 * x20 + x51 * x67 + x90
    x154 = x33 * x57
    x155 = x107 + x51 * x87
    x156 = x155 * x57
    x157 = x110 + x155 * x20 + x51 * x91
    x158 = x123 * x5
    x159 = x108 * x51 + x125
    x160 = x123 * (x113 * x51 + x127 + x159 * x20)
    x161 = x14 * x99
    x162 = x14 * x155
    x163 = x126 * x51 + 4.0 * x2 * (x105 + 2.0 * x106)
    x164 = x123 * x163

    # 45 item(s)
    result[0, 0] = numpy.sum(
        x54
        * (
            x44 * (x19 * x40 + 4.0 * x2 * (x18 + x29 + x34) + x39 * x7)
            + x47 * x49
            + x49 * x52
        )
    )
    result[0, 1] = numpy.sum(x64 * (x52 * x56 * x63 + x55 * x61 + x56 * x62))
    result[0, 2] = numpy.sum(x64 * (x47 * x63 * x65 + x55 * x68 + x62 * x65))
    result[0, 3] = numpy.sum(x81 * (x52 * x72 * x78 + x76 * x78 + x79 * x80))
    result[0, 4] = numpy.sum(x83 * (x56 * x79 * x82 + x58 * x67 * x77 + x61 * x65 * x77))
    result[0, 5] = numpy.sum(x81 * (x47 * x87 * x92 + x79 * x93 + x91 * x92))
    result[0, 6] = numpy.sum(x64 * (x100 * x52 * x94 + x100 * x95 + x104 * x60 * x94))
    result[0, 7] = numpy.sum(x83 * (x66 * x72 * x95 + x66 * x76 * x94 + x67 * x72 * x94))
    result[0, 8] = numpy.sum(x83 * (x58 * x87 * x95 + x58 * x91 * x94 + x59 * x87 * x94))
    result[0, 9] = numpy.sum(x64 * (x109 * x47 * x94 + x109 * x95 + x113 * x57 * x94))
    result[0, 10] = numpy.sum(x54 * (x115 * x117 * x60 + x115 * x122 * x52 + x121 * x7))
    result[0, 11] = numpy.sum(
        x64 * (x104 * x122 * x65 + x116 * x67 * x99 + x117 * x66 * x99)
    )
    result[0, 12] = numpy.sum(
        x81 * (x116 * x72 * x91 + x116 * x76 * x87 + x117 * x72 * x87)
    )
    result[0, 13] = numpy.sum(
        x64 * (x108 * x116 * x59 + x108 * x117 * x58 + x113 * x124 * x56)
    )
    result[0, 14] = numpy.sum(x54 * (x117 * x126 * x57 + x124 * x126 * x47 + x129 * x7))
    result[1, 0] = numpy.sum(x54 * (x131 * x38 * x60 + x132 * x46 + x133 * x46 * x52))
    result[1, 1] = numpy.sum(x64 * (x134 * x136 * x52 + x134 * x36 * x60 + x135 * x136))
    result[1, 2] = numpy.sum(x64 * (x131 * x33 * x66 + x33 * x46 * x68 + x36 * x46 * x82))
    result[1, 3] = numpy.sum(x81 * (x138 * x25 * x52 + x138 * x28 + x139 * x140))
    result[1, 4] = numpy.sum(
        x83 * (x134 * x25 * x67 + x134 * x28 * x66 + x135 * x25 * x66)
    )
    result[1, 5] = numpy.sum(
        x81 * (x130 * x25 * x91 + x131 * x25 * x87 + x28 * x46 * x93)
    )
    result[1, 6] = numpy.sum(x64 * (x141 * x143 * x52 + x141 * x21 * x60 + x142 * x5))
    result[1, 7] = numpy.sum(
        x83 * (x137 * x16 * x67 + x137 * x21 * x66 + x139 * x143 * x65)
    )
    result[1, 8] = numpy.sum(
        x83 * (x134 * x16 * x91 + x134 * x21 * x87 + x135 * x16 * x87)
    )
    result[1, 9] = numpy.sum(
        x64 * (x108 * x131 * x16 + x109 * x21 * x46 + x113 * x144 * x5)
    )
    result[1, 10] = numpy.sum(
        x54
        * (
            x120 * (x119 * x46 + x145 * x20 + 4.0 * x2 * (x101 + x102 + x103))
            + x146 * x52
            + x146 * x9
        )
    )
    result[1, 11] = numpy.sum(x64 * (x14 * x141 * x67 + x141 * x147 * x65 + x142 * x65))
    result[1, 12] = numpy.sum(x81 * (x139 * x14 * x87 + x148 * x87 * x9 + x148 * x91))
    result[1, 13] = numpy.sum(x64 * (x113 * x134 * x14 + x134 * x149 * x9 + x135 * x149))
    result[1, 14] = numpy.sum(x54 * (x126 * x131 * x14 + x126 * x144 * x9 + x129 * x46))
    result[2, 0] = numpy.sum(x54 * (x132 * x51 + x133 * x47 * x51 + x151 * x38 * x57))
    result[2, 1] = numpy.sum(
        x64 * (x151 * x33 * x58 + x33 * x51 * x61 + x36 * x44 * x51 * x56)
    )
    result[2, 2] = numpy.sum(x64 * (x152 * x154 * x47 + x152 * x36 * x57 + x153 * x154))
    result[2, 3] = numpy.sum(
        x81 * (x140 * x51 * x76 + x151 * x25 * x72 + x28 * x51 * x80)
    )
    result[2, 4] = numpy.sum(
        x83 * (x152 * x25 * x59 + x152 * x28 * x58 + x153 * x25 * x58)
    )
    result[2, 5] = numpy.sum(x81 * (x156 * x25 * x47 + x156 * x28 + x157 * x25 * x57))
    result[2, 6] = numpy.sum(
        x64 * (x100 * x21 * x51 + x104 * x143 * x51 + x151 * x16 * x99)
    )
    result[2, 7] = numpy.sum(
        x83 * (x152 * x16 * x76 + x152 * x21 * x72 + x153 * x16 * x72)
    )
    result[2, 8] = numpy.sum(
        x83 * (x155 * x16 * x59 + x155 * x21 * x58 + x157 * x158 * x56)
    )
    result[2, 9] = numpy.sum(x64 * (x158 * x159 * x47 + x159 * x21 * x57 + x160 * x5))
    result[2, 10] = numpy.sum(x54 * (x115 * x14 * x151 + x115 * x147 * x51 + x121 * x51))
    result[2, 11] = numpy.sum(x64 * (x104 * x14 * x152 + x152 * x161 * x9 + x153 * x161))
    result[2, 12] = numpy.sum(x81 * (x14 * x157 * x72 + x162 * x72 * x9 + x162 * x76))
    result[2, 13] = numpy.sum(
        x64 * (x123 * x159 * x56 * x9 + x14 * x159 * x59 + x160 * x56)
    )
    result[2, 14] = numpy.sum(
        x54
        * (
            x123 * (x128 * x51 + x163 * x20 + 4.0 * x2 * (x110 + x111 + x112))
            + x164 * x47
            + x164 * x9
        )
    )
    return result


def kinetic3d_20(ax, da, A, bx, db, B):
    """Cartesian 3D (ds) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 1), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = x0 * (ax * A[0] + bx * B[0]) - A[0]
    x2 = -ax
    x3 = x1**2
    x4 = 2.0 * ax
    x5 = (2.0 * bx + x4) ** (-1.0)
    x6 = 2.0 * ax**2
    x7 = -x2 - x6 * (x3 + x5)
    x8 = bx * x0
    x9 = ax * x8
    x10 = numpy.exp(-x9 * (A[0] - B[0]) ** 2)
    x11 = numpy.sqrt(x0)
    x12 = 1.77245385090552 * x11
    x13 = x10 * x12
    x14 = x1 * x13
    x15 = x4 * x8
    x16 = x14 * (x15 + x7)
    x17 = x13 * x5
    x18 = x13 * x3 + x17
    x19 = numpy.exp(-x9 * (A[1] - B[1]) ** 2)
    x20 = numpy.exp(-x9 * (A[2] - B[2]) ** 2)
    x21 = 3.14159265358979 * x0 * x20
    x22 = x19 * x21
    x23 = x0 * (ax * A[1] + bx * B[1]) - A[1]
    x24 = x23**2
    x25 = -x2 - x6 * (x24 + x5)
    x26 = x18 * x22
    x27 = x0 * (ax * A[2] + bx * B[2]) - A[2]
    x28 = x27**2
    x29 = -x2 - x6 * (x28 + x5)
    x30 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**1.5)
    x31 = 6.53197264742181 * x30
    x32 = x12 * x19
    x33 = x23 * x32
    x34 = x33 * (x15 + x25)
    x35 = x10 * x21
    x36 = x34 * x35
    x37 = x16 * x22
    x38 = 5.56832799683171
    x39 = x0 * x10 * x19
    x40 = x1 * x11 * x20 * x38 * x39
    x41 = 11.3137084989848 * x30
    x42 = x12 * x20
    x43 = x27 * x42
    x44 = x43 * (x15 + x29)
    x45 = 3.14159265358979 * x39
    x46 = x44 * x45
    x47 = x32 * x5
    x48 = x24 * x32 + x47
    x49 = x35 * x48
    x50 = x42 * x5
    x51 = x28 * x42 + x50
    x52 = x45 * x51

    # 6 item(s)
    result[0, 0] = numpy.sum(
        x31
        * (
            x22 * (x1 * x16 + x17 * x7 - x8 * (2.0 * x13 - x18 * x4))
            + x25 * x26
            + x26 * x29
        )
    )
    result[1, 0] = numpy.sum(x41 * (x1 * x36 + x23 * x29 * x40 + x23 * x37))
    result[2, 0] = numpy.sum(x41 * (x1 * x46 + x25 * x27 * x40 + x27 * x37))
    result[3, 0] = numpy.sum(
        x31
        * (
            x29 * x49
            + x35 * (x23 * x34 + x25 * x47 - x8 * (2.0 * x32 - x4 * x48))
            + x49 * x7
        )
    )
    result[4, 0] = numpy.sum(
        x41 * (x11 * x20 * x23 * x27 * x38 * x39 * x7 + x23 * x46 + x27 * x36)
    )
    result[5, 0] = numpy.sum(
        x31
        * (
            x25 * x52
            + x45 * (x27 * x44 + x29 * x50 + x8 * (x4 * x51 - 2.0 * x42))
            + x52 * x7
        )
    )
    return result


def kinetic3d_21(ax, da, A, bx, db, B):
    """Cartesian 3D (dp) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 3), dtype=float)

    x0 = 2.0 * ax
    x1 = (2.0 * bx + x0) ** (-1.0)
    x2 = -ax
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - A[0]
    x6 = x5**2
    x7 = 2.0 * ax**2
    x8 = -x2 - x7 * (x1 + x6)
    x9 = -x4 - B[0]
    x10 = bx * x3
    x11 = ax * x10
    x12 = numpy.exp(-x11 * (A[0] - B[0]) ** 2)
    x13 = 1.77245385090552 * numpy.sqrt(x3)
    x14 = x12 * x13
    x15 = x14 * x9
    x16 = x0 * x10
    x17 = x15 * (x16 + x8)
    x18 = x14 * x5
    x19 = x18 * (x16 + x8)
    x20 = x1 * x14
    x21 = x18 * x9 + x20
    x22 = x20 * x8
    x23 = x16 * x21 + x17 * x5 + x22
    x24 = x1 * (x15 + x18) + x21 * x5
    x25 = numpy.exp(-x11 * (A[1] - B[1]) ** 2)
    x26 = numpy.exp(-x11 * (A[2] - B[2]) ** 2)
    x27 = 3.14159265358979 * x26 * x3
    x28 = x25 * x27
    x29 = -x3 * (ax * A[1] + bx * B[1])
    x30 = -x29 - A[1]
    x31 = x30**2
    x32 = -x2 - x7 * (x1 + x31)
    x33 = x24 * x28
    x34 = -x3 * (ax * A[2] + bx * B[2])
    x35 = -x34 - A[2]
    x36 = x35**2
    x37 = -x2 - x7 * (x1 + x36)
    x38 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**2.5)
    x39 = 13.0639452948436 * x38
    x40 = -x29 - B[1]
    x41 = x13 * x25
    x42 = x40 * x41
    x43 = x42 * (x16 + x32)
    x44 = x14 * x6 + x20
    x45 = x13 * x26
    x46 = x28 * (x10 * (x0 * x44 - 2.0 * x14) + x19 * x5 + x22)
    x47 = x28 * x44
    x48 = -x34 - B[2]
    x49 = x45 * x48
    x50 = x49 * (x16 + x37)
    x51 = x30 * x41
    x52 = x51 * (x16 + x32)
    x53 = x23 * x28
    x54 = x28 * x30
    x55 = 22.6274169979695 * x38
    x56 = x1 * x13
    x57 = x25 * x56
    x58 = x40 * x51 + x57
    x59 = x32 * x57
    x60 = x16 * x58 + x30 * x43 + x59
    x61 = x12 * x27
    x62 = x60 * x61
    x63 = x5 * x61
    x64 = 3.14159265358979 * x12 * x25 * x3
    x65 = x5 * x64
    x66 = x35 * x45
    x67 = x66 * (x16 + x37)
    x68 = x28 * x35
    x69 = x26 * x56
    x70 = x48 * x66 + x69
    x71 = x37 * x69
    x72 = x16 * x70 + x35 * x50 + x71
    x73 = x64 * x72
    x74 = x31 * x41 + x57
    x75 = 2.0 * x41
    x76 = x61 * (x10 * (x0 * x74 - x75) + x30 * x52 + x59)
    x77 = x61 * x74
    x78 = x1 * (x42 + x51) + x30 * x58
    x79 = x61 * x78
    x80 = x30 * x64
    x81 = x35 * x61
    x82 = x36 * x45 + x69
    x83 = 2.0 * x45
    x84 = x64 * (x10 * (x0 * x82 - x83) + x35 * x67 + x71)
    x85 = x64 * x82
    x86 = x1 * (x49 + x66) + x35 * x70
    x87 = x64 * x86

    # 18 item(s)
    result[0, 0] = numpy.sum(
        x39
        * (
            x28 * (x1 * (x17 + x19) + x10 * (x0 * x24 - 2.0 * x15) + x23 * x5)
            + x32 * x33
            + x33 * x37
        )
    )
    result[0, 1] = numpy.sum(x39 * (x37 * x40 * x47 + x40 * x46 + x43 * x44 * x45))
    result[0, 2] = numpy.sum(x39 * (x32 * x47 * x48 + x41 * x44 * x50 + x46 * x48))
    result[1, 0] = numpy.sum(x55 * (x21 * x37 * x54 + x21 * x45 * x52 + x30 * x53))
    result[1, 1] = numpy.sum(x55 * (x19 * x45 * x58 + x37 * x58 * x63 + x5 * x62))
    result[1, 2] = numpy.sum(x55 * (x19 * x48 * x54 + x30 * x50 * x65 + x48 * x52 * x63))
    result[2, 0] = numpy.sum(x55 * (x21 * x32 * x68 + x21 * x41 * x67 + x35 * x53))
    result[2, 1] = numpy.sum(x55 * (x19 * x40 * x68 + x35 * x43 * x63 + x40 * x65 * x67))
    result[2, 2] = numpy.sum(x55 * (x19 * x41 * x70 + x32 * x65 * x70 + x5 * x73))
    result[3, 0] = numpy.sum(x39 * (x17 * x45 * x74 + x37 * x77 * x9 + x76 * x9))
    result[3, 1] = numpy.sum(
        x39
        * (
            x37 * x79
            + x61 * (x1 * (x43 + x52) + x10 * (x0 * x78 - x40 * x75) + x30 * x60)
            + x79 * x8
        )
    )
    result[3, 2] = numpy.sum(x39 * (x14 * x50 * x74 + x48 * x76 + x48 * x77 * x8))
    result[4, 0] = numpy.sum(x55 * (x17 * x35 * x54 + x52 * x81 * x9 + x67 * x80 * x9))
    result[4, 1] = numpy.sum(x55 * (x14 * x58 * x67 + x35 * x62 + x58 * x8 * x81))
    result[4, 2] = numpy.sum(x55 * (x14 * x52 * x70 + x30 * x73 + x70 * x8 * x80))
    result[5, 0] = numpy.sum(x39 * (x17 * x41 * x82 + x32 * x85 * x9 + x84 * x9))
    result[5, 1] = numpy.sum(x39 * (x14 * x43 * x82 + x40 * x8 * x85 + x40 * x84))
    result[5, 2] = numpy.sum(
        x39
        * (
            x32 * x87
            + x64 * (x1 * (x50 + x67) + x10 * (x0 * x86 - x48 * x83) + x35 * x72)
            + x8 * x87
        )
    )
    return result


def kinetic3d_22(ax, da, A, bx, db, B):
    """Cartesian 3D (dd) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 6), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - A[0]
    x6 = -ax
    x7 = x5**2
    x8 = 2.0 * ax**2
    x9 = -x6 - x8 * (x2 + x7)
    x10 = -x4 - B[0]
    x11 = ax * x3
    x12 = bx * x11
    x13 = numpy.exp(-x12 * (A[0] - B[0]) ** 2)
    x14 = 1.77245385090552 * numpy.sqrt(x3)
    x15 = x13 * x14
    x16 = x10 * x15
    x17 = bx * x3
    x18 = x0 * x17
    x19 = x16 * (x18 + x9)
    x20 = x19 * x5
    x21 = x15 * x2
    x22 = x15 * x5
    x23 = x10 * x22 + x21
    x24 = 4.0 * x12
    x25 = x21 * x9
    x26 = x10**2 * x15
    x27 = x21 + x26
    x28 = 2.0 * x15
    x29 = -x28
    x30 = x10 * x19 + x11 * (x1 * x27 + x29)
    x31 = x10 * x28
    x32 = x25 + x30
    x33 = 2.0 * x21
    x34 = x10 * x33 + x27 * x5
    x35 = x18 * x34 + x2 * (x16 * x24 + x31 * x9) + x32 * x5
    x36 = x2 * (3.0 * x21 + x26 + x31 * x5) + x34 * x5
    x37 = numpy.exp(-x12 * (A[1] - B[1]) ** 2)
    x38 = numpy.exp(-x12 * (A[2] - B[2]) ** 2)
    x39 = 3.14159265358979 * x3 * x38
    x40 = x37 * x39
    x41 = -x3 * (ax * A[1] + bx * B[1])
    x42 = -x41 - A[1]
    x43 = x42**2
    x44 = -x6 - x8 * (x2 + x43)
    x45 = x36 * x40
    x46 = -x3 * (ax * A[2] + bx * B[2])
    x47 = -x46 - A[2]
    x48 = x47**2
    x49 = -x6 - x8 * (x2 + x48)
    x50 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**3.5)
    x51 = 15.084944665313 * x50
    x52 = -x41 - B[1]
    x53 = x14 * x37
    x54 = x52 * x53
    x55 = x54 * (x18 + x44)
    x56 = x2 * (x16 + x22) + x23 * x5
    x57 = x14 * x38
    x58 = x22 * (x18 + x9)
    x59 = x18 * x23 + x20 + x25
    x60 = x40 * (x17 * (x0 * x56 - x31) + x2 * (x19 + x58) + x5 * x59)
    x61 = x40 * x56
    x62 = 26.1278905896872 * x50
    x63 = -x46 - B[2]
    x64 = x57 * x63
    x65 = x64 * (x18 + x49)
    x66 = x2 * x53
    x67 = x44 * x66
    x68 = x52**2 * x53
    x69 = x66 + x68
    x70 = 2.0 * x53
    x71 = -x70
    x72 = x11 * (x1 * x69 + x71) + x52 * x55
    x73 = x67 + x72
    x74 = x15 * x7 + x21
    x75 = x57 * x74
    x76 = x17 * (x0 * x74 + x29) + x25 + x5 * x58
    x77 = x40 * x63
    x78 = x2 * x57
    x79 = x49 * x78
    x80 = x57 * x63**2
    x81 = x78 + x80
    x82 = 2.0 * x57
    x83 = -x82
    x84 = x11 * (x1 * x81 + x83) + x63 * x65
    x85 = x79 + x84
    x86 = x53 * x74
    x87 = x42 * x53
    x88 = x87 * (x18 + x44)
    x89 = x35 * x40
    x90 = x34 * x40
    x91 = x42 * x55
    x92 = x52 * x87 + x66
    x93 = x18 * x92 + x67 + x91
    x94 = x23 * x57
    x95 = 45.2548339959391 * x50
    x96 = 2.0 * x66
    x97 = x42 * x69 + x52 * x96
    x98 = x52 * x70
    x99 = x18 * x97 + x2 * (x24 * x54 + x44 * x98) + x42 * x73
    x100 = x13 * x39
    x101 = x100 * x99
    x102 = x100 * x5
    x103 = 3.14159265358979 * x13 * x3 * x37
    x104 = x103 * x5
    x105 = x47 * x57
    x106 = x105 * (x18 + x49)
    x107 = x40 * x47
    x108 = x47 * x65
    x109 = x105 * x63 + x78
    x110 = x108 + x109 * x18 + x79
    x111 = x23 * x53
    x112 = 2.0 * x78
    x113 = x112 * x63 + x47 * x81
    x114 = x63 * x82
    x115 = x113 * x18 + x2 * (x114 * x49 + x24 * x64) + x47 * x85
    x116 = x103 * x115
    x117 = x43 * x53 + x66
    x118 = x117 * x57
    x119 = x17 * (x0 * x117 + x71) + x42 * x88 + x67
    x120 = x2 * (x54 + x87) + x42 * x92
    x121 = x100 * (x17 * (x0 * x120 - x98) + x2 * (x55 + x88) + x42 * x93)
    x122 = x10 * x100
    x123 = x2 * (x42 * x98 + 3.0 * x66 + x68) + x42 * x97
    x124 = x100 * x123
    x125 = x100 * x9
    x126 = x117 * x15
    x127 = x103 * x42
    x128 = x15 * x92
    x129 = x48 * x57 + x78
    x130 = x129 * x53
    x131 = x106 * x47 + x17 * (x0 * x129 + x83) + x79
    x132 = x10 * x103
    x133 = x109 * x47 + x2 * (x105 + x64)
    x134 = x103 * (x110 * x47 + x17 * (x0 * x133 - x114) + x2 * (x106 + x65))
    x135 = x129 * x15
    x136 = x113 * x47 + x2 * (x114 * x47 + 3.0 * x78 + x80)
    x137 = x103 * x136

    # 36 item(s)
    result[0, 0] = numpy.sum(
        x51
        * (
            x40
            * (
                -x17 * (-2.0 * ax * x36 + 2.0 * x26 + x33)
                + x2 * (2.0 * x20 + x23 * x24 + 3.0 * x25 + x30)
                + x35 * x5
            )
            + x44 * x45
            + x45 * x49
        )
    )
    result[0, 1] = numpy.sum(x62 * (x49 * x52 * x61 + x52 * x60 + x55 * x56 * x57))
    result[0, 2] = numpy.sum(x62 * (x44 * x61 * x63 + x53 * x56 * x65 + x60 * x63))
    result[0, 3] = numpy.sum(x51 * (x49 * x69 * x75 + x57 * x69 * x76 + x73 * x75))
    result[0, 4] = numpy.sum(x62 * (x52 * x76 * x77 + x54 * x65 * x74 + x55 * x64 * x74))
    result[0, 5] = numpy.sum(x51 * (x44 * x81 * x86 + x53 * x76 * x81 + x85 * x86))
    result[1, 0] = numpy.sum(x62 * (x34 * x57 * x88 + x42 * x49 * x90 + x42 * x89))
    result[1, 1] = numpy.sum(x95 * (x49 * x92 * x94 + x57 * x59 * x92 + x93 * x94))
    result[1, 2] = numpy.sum(x95 * (x23 * x64 * x88 + x23 * x65 * x87 + x42 * x59 * x77))
    result[1, 3] = numpy.sum(x62 * (x101 * x5 + x102 * x49 * x97 + x57 * x58 * x97))
    result[1, 4] = numpy.sum(x95 * (x102 * x63 * x93 + x22 * x65 * x92 + x58 * x64 * x92))
    result[1, 5] = numpy.sum(x62 * (x104 * x42 * x85 + x22 * x81 * x88 + x58 * x81 * x87))
    result[2, 0] = numpy.sum(x62 * (x106 * x34 * x53 + x44 * x47 * x90 + x47 * x89))
    result[2, 1] = numpy.sum(
        x95 * (x105 * x23 * x55 + x106 * x23 * x54 + x107 * x52 * x59)
    )
    result[2, 2] = numpy.sum(x95 * (x109 * x111 * x44 + x109 * x53 * x59 + x110 * x111))
    result[2, 3] = numpy.sum(
        x62 * (x102 * x47 * x73 + x105 * x58 * x69 + x106 * x22 * x69)
    )
    result[2, 4] = numpy.sum(
        x95 * (x104 * x110 * x52 + x109 * x22 * x55 + x109 * x54 * x58)
    )
    result[2, 5] = numpy.sum(x62 * (x104 * x113 * x44 + x113 * x53 * x58 + x116 * x5))
    result[3, 0] = numpy.sum(x51 * (x118 * x27 * x49 + x118 * x32 + x119 * x27 * x57))
    result[3, 1] = numpy.sum(x62 * (x10 * x121 + x120 * x122 * x49 + x120 * x19 * x57))
    result[3, 2] = numpy.sum(
        x62 * (x117 * x16 * x65 + x117 * x19 * x64 + x119 * x122 * x63)
    )
    result[3, 3] = numpy.sum(
        x51
        * (
            x100
            * (
                -x17 * (-2.0 * ax * x123 + 2.0 * x68 + x96)
                + x2 * (x24 * x92 + 3.0 * x67 + x72 + 2.0 * x91)
                + x42 * x99
            )
            + x124 * x49
            + x124 * x9
        )
    )
    result[3, 4] = numpy.sum(x62 * (x120 * x125 * x63 + x120 * x15 * x65 + x121 * x63))
    result[3, 5] = numpy.sum(x51 * (x119 * x15 * x81 + x126 * x81 * x9 + x126 * x85))
    result[4, 0] = numpy.sum(
        x62 * (x105 * x27 * x88 + x106 * x27 * x87 + x107 * x32 * x42)
    )
    result[4, 1] = numpy.sum(
        x95 * (x105 * x19 * x92 + x106 * x16 * x92 + x122 * x47 * x93)
    )
    result[4, 2] = numpy.sum(
        x95 * (x10 * x110 * x127 + x109 * x16 * x88 + x109 * x19 * x87)
    )
    result[4, 3] = numpy.sum(x62 * (x101 * x47 + x106 * x15 * x97 + x125 * x47 * x97))
    result[4, 4] = numpy.sum(x95 * (x109 * x128 * x9 + x109 * x15 * x93 + x110 * x128))
    result[4, 5] = numpy.sum(x62 * (x113 * x127 * x9 + x113 * x15 * x88 + x116 * x42))
    result[5, 0] = numpy.sum(x51 * (x130 * x27 * x44 + x130 * x32 + x131 * x27 * x53))
    result[5, 1] = numpy.sum(
        x62 * (x129 * x16 * x55 + x129 * x19 * x54 + x131 * x132 * x52)
    )
    result[5, 2] = numpy.sum(x62 * (x10 * x134 + x132 * x133 * x44 + x133 * x19 * x53))
    result[5, 3] = numpy.sum(x51 * (x131 * x15 * x69 + x135 * x69 * x9 + x135 * x73))
    result[5, 4] = numpy.sum(
        x62 * (x103 * x133 * x52 * x9 + x133 * x15 * x55 + x134 * x52)
    )
    result[5, 5] = numpy.sum(
        x51
        * (
            x103
            * (
                x115 * x47
                - x17 * (-2.0 * ax * x136 + x112 + 2.0 * x80)
                + x2 * (2.0 * x108 + x109 * x24 + 3.0 * x79 + x84)
            )
            + x137 * x44
            + x137 * x9
        )
    )
    return result


def kinetic3d_23(ax, da, A, bx, db, B):
    """Cartesian 3D (df) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 10), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - A[0]
    x6 = -ax
    x7 = x5**2
    x8 = 2.0 * ax**2
    x9 = -x6 - x8 * (x2 + x7)
    x10 = ax * x3
    x11 = bx * x10
    x12 = numpy.exp(-x11 * (A[0] - B[0]) ** 2)
    x13 = 1.77245385090552 * numpy.sqrt(x3)
    x14 = x12 * x13
    x15 = x14 * x2
    x16 = x15 * x9
    x17 = -x4 - B[0]
    x18 = x14 * x17
    x19 = bx * x3
    x20 = x0 * x19
    x21 = x18 * (x20 + x9)
    x22 = x17 * x21
    x23 = x14 * x17**2
    x24 = x15 + x23
    x25 = 2.0 * x14
    x26 = -x25
    x27 = x10 * (x1 * x24 + x26)
    x28 = x22 + x27
    x29 = x16 + x28
    x30 = x29 * x5
    x31 = x17 * x25
    x32 = 4.0 * x11
    x33 = x2 * (x18 * x32 + x31 * x9)
    x34 = x24 * x5
    x35 = 2.0 * x15
    x36 = x17 * x35
    x37 = x34 + x36
    x38 = 6.0 * x11
    x39 = x17 * x24
    x40 = x36 + x39
    x41 = x10 * (x1 * x40 - 3.0 * x18) + x17 * x29
    x42 = 3.0 * x16
    x43 = x33 + x41
    x44 = 3.0 * x15
    x45 = x2 * (3.0 * x23 + x44) + x40 * x5
    x46 = x2 * (3.0 * x22 + 3.0 * x27 + x42) + x20 * x45 + x43 * x5
    x47 = x15 * x17
    x48 = x2 * (3.0 * x34 + x39 + 8.0 * x47) + x45 * x5
    x49 = numpy.exp(-x11 * (A[1] - B[1]) ** 2)
    x50 = numpy.exp(-x11 * (A[2] - B[2]) ** 2)
    x51 = 3.14159265358979 * x3 * x50
    x52 = x49 * x51
    x53 = -x3 * (ax * A[1] + bx * B[1])
    x54 = -x53 - A[1]
    x55 = x54**2
    x56 = -x6 - x8 * (x2 + x55)
    x57 = x48 * x52
    x58 = -x3 * (ax * A[2] + bx * B[2])
    x59 = -x58 - A[2]
    x60 = x59**2
    x61 = -x6 - x8 * (x2 + x60)
    x62 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**4.5)
    x63 = 13.4923846833851 * x62
    x64 = -x53 - B[1]
    x65 = x13 * x49
    x66 = x64 * x65
    x67 = x66 * (x20 + x56)
    x68 = x2 * (x23 + x31 * x5 + x44) + x37 * x5
    x69 = x13 * x50
    x70 = x21 * x5
    x71 = x14 * x5
    x72 = x15 + x17 * x71
    x73 = x20 * x37 + x30 + x33
    x74 = x52 * (
        -x19 * (-2.0 * ax * x68 + 2.0 * x23 + x35)
        + x2 * (x28 + x32 * x72 + x42 + 2.0 * x70)
        + x5 * x73
    )
    x75 = x52 * x68
    x76 = 30.169889330626 * x62
    x77 = -x58 - B[2]
    x78 = x69 * x77
    x79 = x78 * (x20 + x61)
    x80 = x2 * x65
    x81 = x56 * x80
    x82 = x64 * x67
    x83 = x64**2 * x65
    x84 = x80 + x83
    x85 = 2.0 * x65
    x86 = -x85
    x87 = x10 * (x1 * x84 + x86)
    x88 = x82 + x87
    x89 = x81 + x88
    x90 = x2 * (x18 + x71) + x5 * x72
    x91 = x69 * x90
    x92 = x71 * (x20 + x9)
    x93 = x16 + x20 * x72 + x70
    x94 = x19 * (x0 * x90 - x31) + x2 * (x21 + x92) + x5 * x93
    x95 = x52 * x77
    x96 = 52.2557811793745 * x62
    x97 = x2 * x69
    x98 = x61 * x97
    x99 = x77 * x79
    x100 = x69 * x77**2
    x101 = x100 + x97
    x102 = 2.0 * x69
    x103 = -x102
    x104 = x10 * (x1 * x101 + x103)
    x105 = x104 + x99
    x106 = x105 + x98
    x107 = x65 * x90
    x108 = x14 * x7 + x15
    x109 = x16 + x19 * (x0 * x108 + x26) + x5 * x92
    x110 = x64 * x84
    x111 = 2.0 * x80
    x112 = x111 * x64
    x113 = x110 + x112
    x114 = x113 * x69
    x115 = x64 * x85
    x116 = x2 * (x115 * x56 + x32 * x66)
    x117 = x10 * (x1 * x113 - 3.0 * x66) + x64 * x89
    x118 = x116 + x117
    x119 = x101 * x77
    x120 = 2.0 * x97
    x121 = x120 * x77
    x122 = x119 + x121
    x123 = x122 * x65
    x124 = x102 * x77
    x125 = x2 * (x124 * x61 + x32 * x78)
    x126 = x10 * (x1 * x122 - 3.0 * x78) + x106 * x77
    x127 = x125 + x126
    x128 = x54 * x65
    x129 = x128 * (x20 + x56)
    x130 = x46 * x52
    x131 = x45 * x52
    x132 = 23.3694957868871 * x62
    x133 = x54 * x67
    x134 = x128 * x64 + x80
    x135 = x133 + x134 * x20 + x81
    x136 = x37 * x69
    x137 = x54 * x84
    x138 = x112 + x137
    x139 = x138 * x69
    x140 = x54 * x89
    x141 = x116 + x138 * x20 + x140
    x142 = 90.5096679918781 * x62
    x143 = 3.0 * x80
    x144 = x113 * x54 + x2 * (x143 + 3.0 * x83)
    x145 = 3.0 * x81
    x146 = x118 * x54 + x144 * x20 + x2 * (x145 + 3.0 * x82 + 3.0 * x87)
    x147 = x12 * x51
    x148 = x146 * x147
    x149 = x147 * x5
    x150 = 3.14159265358979 * x12 * x3 * x49
    x151 = x150 * x5
    x152 = x59 * x69
    x153 = x152 * (x20 + x61)
    x154 = x52 * x59
    x155 = x59 * x79
    x156 = x152 * x77 + x97
    x157 = x155 + x156 * x20 + x98
    x158 = x37 * x65
    x159 = x101 * x59
    x160 = x121 + x159
    x161 = x160 * x65
    x162 = x106 * x59
    x163 = x125 + x160 * x20 + x162
    x164 = 3.0 * x97
    x165 = x122 * x59 + x2 * (3.0 * x100 + x164)
    x166 = 3.0 * x98
    x167 = x127 * x59 + x165 * x20 + x2 * (3.0 * x104 + x166 + 3.0 * x99)
    x168 = x150 * x167
    x169 = x55 * x65 + x80
    x170 = x129 * x54 + x19 * (x0 * x169 + x86) + x81
    x171 = x40 * x69
    x172 = x134 * x54 + x2 * (x128 + x66)
    x173 = x172 * x69
    x174 = x135 * x54 + x19 * (x0 * x172 - x115) + x2 * (x129 + x67)
    x175 = x138 * x54 + x2 * (x115 * x54 + x143 + x83)
    x176 = x147 * (
        x141 * x54
        - x19 * (-2.0 * ax * x175 + x111 + 2.0 * x83)
        + x2 * (2.0 * x133 + x134 * x32 + x145 + x88)
    )
    x177 = x147 * x17
    x178 = x64 * x80
    x179 = x144 * x54 + x2 * (x110 + 3.0 * x137 + 8.0 * x178)
    x180 = x147 * x179
    x181 = x147 * x9
    x182 = x14 * x172
    x183 = x122 * x14
    x184 = x150 * x54
    x185 = x138 * x14
    x186 = x14 * x160
    x187 = x60 * x69 + x97
    x188 = x153 * x59 + x19 * (x0 * x187 + x103) + x98
    x189 = x40 * x65
    x190 = x156 * x59 + x2 * (x152 + x78)
    x191 = x190 * x65
    x192 = x157 * x59 + x19 * (x0 * x190 - x124) + x2 * (x153 + x79)
    x193 = x150 * x17
    x194 = x160 * x59 + x2 * (x100 + x124 * x59 + x164)
    x195 = x150 * (
        x163 * x59
        - x19 * (-2.0 * ax * x194 + 2.0 * x100 + x120)
        + x2 * (x105 + 2.0 * x155 + x156 * x32 + x166)
    )
    x196 = x113 * x14
    x197 = x14 * x190
    x198 = x77 * x97
    x199 = x165 * x59 + x2 * (x119 + 3.0 * x159 + 8.0 * x198)
    x200 = x150 * x199

    # 60 item(s)
    result[0, 0] = numpy.sum(
        x63
        * (
            x52
            * (
                -2.0 * x19 * (-ax * x48 + x39 + 2.0 * x47)
                + x2 * (3.0 * x30 + 4.0 * x33 + x37 * x38 + x41)
                + x46 * x5
            )
            + x56 * x57
            + x57 * x61
        )
    )
    result[0, 1] = numpy.sum(x76 * (x61 * x64 * x75 + x64 * x74 + x67 * x68 * x69))
    result[0, 2] = numpy.sum(x76 * (x56 * x75 * x77 + x65 * x68 * x79 + x74 * x77))
    result[0, 3] = numpy.sum(x76 * (x61 * x84 * x91 + x69 * x84 * x94 + x89 * x91))
    result[0, 4] = numpy.sum(x96 * (x64 * x94 * x95 + x66 * x79 * x90 + x67 * x78 * x90))
    result[0, 5] = numpy.sum(x76 * (x101 * x107 * x56 + x101 * x65 * x94 + x106 * x107))
    result[0, 6] = numpy.sum(x63 * (x108 * x114 * x61 + x108 * x118 * x69 + x109 * x114))
    result[0, 7] = numpy.sum(
        x76 * (x108 * x78 * x89 + x108 * x79 * x84 + x109 * x78 * x84)
    )
    result[0, 8] = numpy.sum(
        x76 * (x101 * x108 * x67 + x101 * x109 * x66 + x106 * x108 * x66)
    )
    result[0, 9] = numpy.sum(x63 * (x108 * x123 * x56 + x108 * x127 * x65 + x109 * x123))
    result[1, 0] = numpy.sum(x132 * (x129 * x45 * x69 + x130 * x54 + x131 * x54 * x61))
    result[1, 1] = numpy.sum(x96 * (x134 * x136 * x61 + x134 * x69 * x73 + x135 * x136))
    result[1, 2] = numpy.sum(
        x96 * (x128 * x37 * x79 + x129 * x37 * x78 + x54 * x73 * x95)
    )
    result[1, 3] = numpy.sum(x96 * (x139 * x61 * x72 + x139 * x93 + x141 * x69 * x72))
    result[1, 4] = numpy.sum(
        x142 * (x134 * x72 * x79 + x134 * x78 * x93 + x135 * x72 * x78)
    )
    result[1, 5] = numpy.sum(
        x96 * (x101 * x128 * x93 + x101 * x129 * x72 + x106 * x128 * x72)
    )
    result[1, 6] = numpy.sum(x132 * (x144 * x149 * x61 + x144 * x69 * x92 + x148 * x5))
    result[1, 7] = numpy.sum(
        x96 * (x138 * x71 * x79 + x138 * x78 * x92 + x141 * x149 * x77)
    )
    result[1, 8] = numpy.sum(
        x96 * (x101 * x134 * x92 + x101 * x135 * x71 + x106 * x134 * x71)
    )
    result[1, 9] = numpy.sum(
        x132 * (x122 * x128 * x92 + x122 * x129 * x71 + x127 * x151 * x54)
    )
    result[2, 0] = numpy.sum(x132 * (x130 * x59 + x131 * x56 * x59 + x153 * x45 * x65))
    result[2, 1] = numpy.sum(
        x96 * (x152 * x37 * x67 + x153 * x37 * x66 + x154 * x64 * x73)
    )
    result[2, 2] = numpy.sum(x96 * (x156 * x158 * x56 + x156 * x65 * x73 + x157 * x158))
    result[2, 3] = numpy.sum(
        x96 * (x152 * x72 * x89 + x152 * x84 * x93 + x153 * x72 * x84)
    )
    result[2, 4] = numpy.sum(
        x142 * (x156 * x66 * x93 + x156 * x67 * x72 + x157 * x66 * x72)
    )
    result[2, 5] = numpy.sum(x96 * (x161 * x56 * x72 + x161 * x93 + x163 * x65 * x72))
    result[2, 6] = numpy.sum(
        x132 * (x113 * x152 * x92 + x113 * x153 * x71 + x118 * x149 * x59)
    )
    result[2, 7] = numpy.sum(
        x96 * (x156 * x71 * x89 + x156 * x84 * x92 + x157 * x71 * x84)
    )
    result[2, 8] = numpy.sum(
        x96 * (x151 * x163 * x64 + x160 * x66 * x92 + x160 * x67 * x71)
    )
    result[2, 9] = numpy.sum(x132 * (x151 * x165 * x56 + x165 * x65 * x92 + x168 * x5))
    result[3, 0] = numpy.sum(x63 * (x169 * x171 * x61 + x169 * x43 * x69 + x170 * x171))
    result[3, 1] = numpy.sum(x76 * (x173 * x24 * x61 + x173 * x29 + x174 * x24 * x69))
    result[3, 2] = numpy.sum(
        x76 * (x169 * x24 * x79 + x169 * x29 * x78 + x170 * x24 * x78)
    )
    result[3, 3] = numpy.sum(x76 * (x17 * x176 + x175 * x177 * x61 + x175 * x21 * x69))
    result[3, 4] = numpy.sum(
        x96 * (x172 * x18 * x79 + x172 * x21 * x78 + x174 * x177 * x77)
    )
    result[3, 5] = numpy.sum(
        x76 * (x101 * x169 * x21 + x101 * x170 * x18 + x106 * x169 * x18)
    )
    result[3, 6] = numpy.sum(
        x63
        * (
            x147
            * (
                x146 * x54
                - 2.0 * x19 * (-ax * x179 + x110 + 2.0 * x178)
                + x2 * (4.0 * x116 + x117 + x138 * x38 + 3.0 * x140)
            )
            + x180 * x61
            + x180 * x9
        )
    )
    result[3, 7] = numpy.sum(x76 * (x14 * x175 * x79 + x175 * x181 * x77 + x176 * x77))
    result[3, 8] = numpy.sum(x76 * (x101 * x14 * x174 + x101 * x182 * x9 + x106 * x182))
    result[3, 9] = numpy.sum(x63 * (x127 * x14 * x169 + x169 * x183 * x9 + x170 * x183))
    result[4, 0] = numpy.sum(
        x132 * (x128 * x153 * x40 + x129 * x152 * x40 + x154 * x43 * x54)
    )
    result[4, 1] = numpy.sum(
        x96 * (x134 * x152 * x29 + x134 * x153 * x24 + x135 * x152 * x24)
    )
    result[4, 2] = numpy.sum(
        x96 * (x128 * x156 * x29 + x128 * x157 * x24 + x129 * x156 * x24)
    )
    result[4, 3] = numpy.sum(
        x96 * (x138 * x152 * x21 + x138 * x153 * x18 + x141 * x177 * x59)
    )
    result[4, 4] = numpy.sum(
        x142 * (x134 * x156 * x21 + x134 * x157 * x18 + x135 * x156 * x18)
    )
    result[4, 5] = numpy.sum(
        x96 * (x128 * x160 * x21 + x129 * x160 * x18 + x163 * x17 * x184)
    )
    result[4, 6] = numpy.sum(x132 * (x14 * x144 * x153 + x144 * x181 * x59 + x148 * x59))
    result[4, 7] = numpy.sum(x96 * (x14 * x141 * x156 + x156 * x185 * x9 + x157 * x185))
    result[4, 8] = numpy.sum(x96 * (x134 * x14 * x163 + x134 * x186 * x9 + x135 * x186))
    result[4, 9] = numpy.sum(x132 * (x129 * x14 * x165 + x165 * x184 * x9 + x168 * x54))
    result[5, 0] = numpy.sum(x63 * (x187 * x189 * x56 + x187 * x43 * x65 + x188 * x189))
    result[5, 1] = numpy.sum(
        x76 * (x187 * x24 * x67 + x187 * x29 * x66 + x188 * x24 * x66)
    )
    result[5, 2] = numpy.sum(x76 * (x191 * x24 * x56 + x191 * x29 + x192 * x24 * x65))
    result[5, 3] = numpy.sum(
        x76 * (x18 * x187 * x89 + x18 * x188 * x84 + x187 * x21 * x84)
    )
    result[5, 4] = numpy.sum(
        x96 * (x18 * x190 * x67 + x190 * x21 * x66 + x192 * x193 * x64)
    )
    result[5, 5] = numpy.sum(x76 * (x17 * x195 + x193 * x194 * x56 + x194 * x21 * x65))
    result[5, 6] = numpy.sum(x63 * (x118 * x14 * x187 + x187 * x196 * x9 + x188 * x196))
    result[5, 7] = numpy.sum(x76 * (x14 * x192 * x84 + x197 * x84 * x9 + x197 * x89))
    result[5, 8] = numpy.sum(
        x76 * (x14 * x194 * x67 + x150 * x194 * x64 * x9 + x195 * x64)
    )
    result[5, 9] = numpy.sum(
        x63
        * (
            x150
            * (
                x167 * x59
                - 2.0 * x19 * (-ax * x199 + x119 + 2.0 * x198)
                + x2 * (4.0 * x125 + x126 + x160 * x38 + 3.0 * x162)
            )
            + x200 * x56
            + x200 * x9
        )
    )
    return result


def kinetic3d_24(ax, da, A, bx, db, B):
    """Cartesian 3D (dg) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((6, 15), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - A[0]
    x6 = -ax
    x7 = x5**2
    x8 = 2.0 * ax**2
    x9 = -x6 - x8 * (x2 + x7)
    x10 = -x4 - B[0]
    x11 = ax * x3
    x12 = bx * x11
    x13 = numpy.exp(-x12 * (A[0] - B[0]) ** 2)
    x14 = 1.77245385090552 * numpy.sqrt(x3)
    x15 = x13 * x14
    x16 = 2.0 * x15
    x17 = x10 * x16
    x18 = x10 * x15
    x19 = 4.0 * x12
    x20 = x2 * (x17 * x9 + x18 * x19)
    x21 = x15 * x2
    x22 = x21 * x9
    x23 = bx * x3
    x24 = x0 * x23
    x25 = x18 * (x24 + x9)
    x26 = x10 * x25
    x27 = x10**2 * x15
    x28 = x21 + x27
    x29 = -x16
    x30 = x11 * (x1 * x28 + x29)
    x31 = x26 + x30
    x32 = x22 + x31
    x33 = x10 * x32
    x34 = x10 * x28
    x35 = 2.0 * x21
    x36 = x10 * x35
    x37 = x34 + x36
    x38 = x11 * (x1 * x37 - 3.0 * x18)
    x39 = x33 + x38
    x40 = x20 + x39
    x41 = x40 * x5
    x42 = 3.0 * x22
    x43 = x2 * (3.0 * x26 + 3.0 * x30 + x42)
    x44 = 3.0 * x21
    x45 = x2 * (3.0 * x27 + x44)
    x46 = x37 * x5
    x47 = x45 + x46
    x48 = 8.0 * x12
    x49 = x10 * x37
    x50 = x45 + x49
    x51 = 4.0 * x21
    x52 = x10 * x40 - x11 * (-2.0 * bx * x50 + 4.0 * x27 + x51)
    x53 = 4.0 * x20
    x54 = x43 + x52
    x55 = 8.0 * x10 * x21
    x56 = x2 * (4.0 * x34 + x55) + x5 * x50
    x57 = x2 * (4.0 * x33 + 4.0 * x38 + x53) + x24 * x56 + x5 * x54
    x58 = x2 * (5.0 * x45 + 4.0 * x46 + x49) + x5 * x56
    x59 = numpy.exp(-x12 * (A[1] - B[1]) ** 2)
    x60 = numpy.exp(-x12 * (A[2] - B[2]) ** 2)
    x61 = 3.14159265358979 * x3 * x60
    x62 = x59 * x61
    x63 = -x3 * (ax * A[1] + bx * B[1])
    x64 = -x63 - A[1]
    x65 = x64**2
    x66 = -x6 - x8 * (x2 + x65)
    x67 = x58 * x62
    x68 = -x3 * (ax * A[2] + bx * B[2])
    x69 = -x68 - A[2]
    x70 = x69**2
    x71 = -x6 - x8 * (x2 + x70)
    x72 = 0.179587122125167 * da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**5.5)
    x73 = 10.1992841329868 * x72
    x74 = -x63 - B[1]
    x75 = x14 * x59
    x76 = x74 * x75
    x77 = x76 * (x24 + x66)
    x78 = x28 * x5
    x79 = x2 * (x34 + x55 + 3.0 * x78) + x47 * x5
    x80 = x14 * x60
    x81 = x32 * x5
    x82 = x36 + x78
    x83 = 6.0 * x12
    x84 = x24 * x47 + x41 + x43
    x85 = x62 * (
        x2 * (x39 + x53 + 3.0 * x81 + x82 * x83)
        - x23 * (-2.0 * ax * x79 + x10 * x51 + 2.0 * x34)
        + x5 * x84
    )
    x86 = x62 * x79
    x87 = 26.9847693667702 * x72
    x88 = -x68 - B[2]
    x89 = x80 * x88
    x90 = x89 * (x24 + x71)
    x91 = x2 * x75
    x92 = x66 * x91
    x93 = x74 * x77
    x94 = x74**2 * x75
    x95 = x91 + x94
    x96 = 2.0 * x75
    x97 = -x96
    x98 = x11 * (x1 * x95 + x97)
    x99 = x93 + x98
    x100 = x92 + x99
    x101 = x2 * (x17 * x5 + x27 + x44) + x5 * x82
    x102 = x101 * x80
    x103 = x25 * x5
    x104 = x15 * x5
    x105 = x10 * x104 + x21
    x106 = x20 + x24 * x82 + x81
    x107 = (
        x106 * x5
        + x2 * (2.0 * x103 + x105 * x19 + x31 + x42)
        - x23 * (-2.0 * ax * x101 + 2.0 * x27 + x35)
    )
    x108 = 34.8371874529163 * x72
    x109 = x62 * x88
    x110 = 60.3397786612521 * x72
    x111 = x2 * x80
    x112 = x111 * x71
    x113 = x88 * x90
    x114 = x80 * x88**2
    x115 = x111 + x114
    x116 = 2.0 * x80
    x117 = -x116
    x118 = x11 * (x1 * x115 + x117)
    x119 = x113 + x118
    x120 = x112 + x119
    x121 = x101 * x75
    x122 = x74 * x96
    x123 = x2 * (x122 * x66 + x19 * x76)
    x124 = x100 * x74
    x125 = x74 * x95
    x126 = 2.0 * x91
    x127 = x126 * x74
    x128 = x125 + x127
    x129 = x11 * (x1 * x128 - 3.0 * x76)
    x130 = x124 + x129
    x131 = x123 + x130
    x132 = x105 * x5 + x2 * (x104 + x18)
    x133 = x132 * x80
    x134 = x104 * (x24 + x9)
    x135 = x103 + x105 * x24 + x22
    x136 = x135 * x5 + x2 * (x134 + x25) + x23 * (x0 * x132 - x17)
    x137 = x116 * x88
    x138 = x2 * (x137 * x71 + x19 * x89)
    x139 = x120 * x88
    x140 = x115 * x88
    x141 = 2.0 * x111
    x142 = x141 * x88
    x143 = x140 + x142
    x144 = x11 * (x1 * x143 - 3.0 * x89)
    x145 = x139 + x144
    x146 = x138 + x145
    x147 = x132 * x75
    x148 = x15 * x7 + x21
    x149 = x134 * x5 + x22 + x23 * (x0 * x148 + x29)
    x150 = 3.0 * x91
    x151 = x2 * (x150 + 3.0 * x94)
    x152 = x128 * x74
    x153 = x151 + x152
    x154 = x153 * x80
    x155 = 3.0 * x92
    x156 = x2 * (x155 + 3.0 * x93 + 3.0 * x98)
    x157 = 4.0 * x91
    x158 = -x11 * (-2.0 * bx * x153 + x157 + 4.0 * x94) + x131 * x74
    x159 = x156 + x158
    x160 = 3.0 * x111
    x161 = x2 * (3.0 * x114 + x160)
    x162 = x143 * x88
    x163 = x161 + x162
    x164 = x163 * x75
    x165 = 3.0 * x112
    x166 = x2 * (3.0 * x113 + 3.0 * x118 + x165)
    x167 = 4.0 * x111
    x168 = -x11 * (-2.0 * bx * x163 + 4.0 * x114 + x167) + x146 * x88
    x169 = x166 + x168
    x170 = x64 * x75
    x171 = x170 * (x24 + x66)
    x172 = x57 * x62
    x173 = x56 * x62
    x174 = 17.6656783191643 * x72
    x175 = x64 * x77
    x176 = x170 * x74 + x91
    x177 = x175 + x176 * x24 + x92
    x178 = x47 * x80
    x179 = 46.7389915737742 * x72
    x180 = x100 * x64
    x181 = x64 * x95
    x182 = x127 + x181
    x183 = x123 + x180 + x182 * x24
    x184 = x80 * x82
    x185 = 60.3397786612521 * x72
    x186 = 104.511562358749 * x72
    x187 = x128 * x64
    x188 = x151 + x187
    x189 = x188 * x80
    x190 = x131 * x64
    x191 = x156 + x188 * x24 + x190
    x192 = 8.0 * x74 * x91
    x193 = x153 * x64 + x2 * (4.0 * x125 + x192)
    x194 = 4.0 * x123
    x195 = x159 * x64 + x193 * x24 + x2 * (4.0 * x124 + 4.0 * x129 + x194)
    x196 = x13 * x61
    x197 = x195 * x196
    x198 = x196 * x5
    x199 = 3.14159265358979 * x13 * x3 * x59
    x200 = x199 * x5
    x201 = x69 * x80
    x202 = x201 * (x24 + x71)
    x203 = x62 * x69
    x204 = x69 * x90
    x205 = x111 + x201 * x88
    x206 = x112 + x204 + x205 * x24
    x207 = x47 * x75
    x208 = x120 * x69
    x209 = x115 * x69
    x210 = x142 + x209
    x211 = x138 + x208 + x210 * x24
    x212 = x75 * x82
    x213 = x143 * x69
    x214 = x161 + x213
    x215 = x214 * x75
    x216 = x146 * x69
    x217 = x166 + x214 * x24 + x216
    x218 = 8.0 * x111 * x88
    x219 = x163 * x69 + x2 * (4.0 * x140 + x218)
    x220 = 4.0 * x138
    x221 = x169 * x69 + x2 * (4.0 * x139 + 4.0 * x144 + x220) + x219 * x24
    x222 = x199 * x221
    x223 = x65 * x75 + x91
    x224 = x171 * x64 + x23 * (x0 * x223 + x97) + x92
    x225 = x50 * x80
    x226 = x176 * x64 + x2 * (x170 + x76)
    x227 = x226 * x80
    x228 = x177 * x64 + x2 * (x171 + x77) + x23 * (x0 * x226 - x122)
    x229 = x182 * x64 + x2 * (x122 * x64 + x150 + x94)
    x230 = x229 * x80
    x231 = (
        x183 * x64
        + x2 * (x155 + 2.0 * x175 + x176 * x19 + x99)
        - x23 * (-2.0 * ax * x229 + x126 + 2.0 * x94)
    )
    x232 = x188 * x64 + x2 * (x125 + 3.0 * x181 + x192)
    x233 = x196 * (
        x191 * x64
        + x2 * (x130 + 3.0 * x180 + x182 * x83 + x194)
        - x23 * (-2.0 * ax * x232 + 2.0 * x125 + x157 * x74)
    )
    x234 = x10 * x196
    x235 = x193 * x64 + x2 * (5.0 * x151 + x152 + 4.0 * x187)
    x236 = x196 * x235
    x237 = x196 * x9
    x238 = x15 * x229
    x239 = x15 * x226
    x240 = x15 * x163
    x241 = x199 * x64
    x242 = x15 * x188
    x243 = x15 * x182
    x244 = x15 * x214
    x245 = x111 + x70 * x80
    x246 = x112 + x202 * x69 + x23 * (x0 * x245 + x117)
    x247 = x50 * x75
    x248 = x2 * (x201 + x89) + x205 * x69
    x249 = x248 * x75
    x250 = x2 * (x202 + x90) + x206 * x69 + x23 * (x0 * x248 - x137)
    x251 = x2 * (x114 + x137 * x69 + x160) + x210 * x69
    x252 = x251 * x75
    x253 = (
        x2 * (x119 + x165 + x19 * x205 + 2.0 * x204)
        + x211 * x69
        - x23 * (-2.0 * ax * x251 + 2.0 * x114 + x141)
    )
    x254 = x10 * x199
    x255 = x2 * (x140 + 3.0 * x209 + x218) + x214 * x69
    x256 = x199 * (
        x2 * (x145 + 3.0 * x208 + x210 * x83 + x220)
        + x217 * x69
        - x23 * (-2.0 * ax * x255 + 2.0 * x140 + x167 * x88)
    )
    x257 = x15 * x153
    x258 = x15 * x248
    x259 = x15 * x251
    x260 = x2 * (5.0 * x161 + x162 + 4.0 * x213) + x219 * x69
    x261 = x199 * x260

    # 90 item(s)
    result[0, 0] = numpy.sum(
        x73
        * (
            x62
            * (
                x2 * (4.0 * x41 + 5.0 * x43 + x47 * x48 + x52)
                - 2.0 * x23 * (-ax * x58 + x45 + x49)
                + x5 * x57
            )
            + x66 * x67
            + x67 * x71
        )
    )
    result[0, 1] = numpy.sum(x87 * (x71 * x74 * x86 + x74 * x85 + x77 * x79 * x80))
    result[0, 2] = numpy.sum(x87 * (x66 * x86 * x88 + x75 * x79 * x90 + x85 * x88))
    result[0, 3] = numpy.sum(x108 * (x100 * x102 + x102 * x71 * x95 + x107 * x80 * x95))
    result[0, 4] = numpy.sum(
        x110 * (x101 * x76 * x90 + x101 * x77 * x89 + x107 * x109 * x74)
    )
    result[0, 5] = numpy.sum(x108 * (x107 * x115 * x75 + x115 * x121 * x66 + x120 * x121))
    result[0, 6] = numpy.sum(x87 * (x128 * x133 * x71 + x128 * x136 * x80 + x131 * x133))
    result[0, 7] = numpy.sum(
        x110 * (x100 * x132 * x89 + x132 * x90 * x95 + x136 * x89 * x95)
    )
    result[0, 8] = numpy.sum(
        x110 * (x115 * x132 * x77 + x115 * x136 * x76 + x120 * x132 * x76)
    )
    result[0, 9] = numpy.sum(x87 * (x136 * x143 * x75 + x143 * x147 * x66 + x146 * x147))
    result[0, 10] = numpy.sum(x73 * (x148 * x154 * x71 + x148 * x159 * x80 + x149 * x154))
    result[0, 11] = numpy.sum(
        x87 * (x128 * x148 * x90 + x128 * x149 * x89 + x131 * x148 * x89)
    )
    result[0, 12] = numpy.sum(
        x108 * (x100 * x115 * x148 + x115 * x149 * x95 + x120 * x148 * x95)
    )
    result[0, 13] = numpy.sum(
        x87 * (x143 * x148 * x77 + x143 * x149 * x76 + x146 * x148 * x76)
    )
    result[0, 14] = numpy.sum(x73 * (x148 * x164 * x66 + x148 * x169 * x75 + x149 * x164))
    result[1, 0] = numpy.sum(x174 * (x171 * x56 * x80 + x172 * x64 + x173 * x64 * x71))
    result[1, 1] = numpy.sum(x179 * (x176 * x178 * x71 + x176 * x80 * x84 + x177 * x178))
    result[1, 2] = numpy.sum(
        x179 * (x109 * x64 * x84 + x170 * x47 * x90 + x171 * x47 * x89)
    )
    result[1, 3] = numpy.sum(x185 * (x106 * x182 * x80 + x182 * x184 * x71 + x183 * x184))
    result[1, 4] = numpy.sum(
        x186 * (x106 * x176 * x89 + x176 * x82 * x90 + x177 * x82 * x89)
    )
    result[1, 5] = numpy.sum(
        x185 * (x106 * x115 * x170 + x115 * x171 * x82 + x120 * x170 * x82)
    )
    result[1, 6] = numpy.sum(x179 * (x105 * x189 * x71 + x105 * x191 * x80 + x135 * x189))
    result[1, 7] = numpy.sum(
        x186 * (x105 * x182 * x90 + x105 * x183 * x89 + x135 * x182 * x89)
    )
    result[1, 8] = numpy.sum(
        x186 * (x105 * x115 * x177 + x105 * x120 * x176 + x115 * x135 * x176)
    )
    result[1, 9] = numpy.sum(
        x179 * (x105 * x143 * x171 + x105 * x146 * x170 + x135 * x143 * x170)
    )
    result[1, 10] = numpy.sum(x174 * (x134 * x193 * x80 + x193 * x198 * x71 + x197 * x5))
    result[1, 11] = numpy.sum(
        x179 * (x104 * x188 * x90 + x134 * x188 * x89 + x191 * x198 * x88)
    )
    result[1, 12] = numpy.sum(
        x185 * (x104 * x115 * x183 + x104 * x120 * x182 + x115 * x134 * x182)
    )
    result[1, 13] = numpy.sum(
        x179 * (x104 * x143 * x177 + x104 * x146 * x176 + x134 * x143 * x176)
    )
    result[1, 14] = numpy.sum(
        x174 * (x104 * x163 * x171 + x134 * x163 * x170 + x169 * x200 * x64)
    )
    result[2, 0] = numpy.sum(x174 * (x172 * x69 + x173 * x66 * x69 + x202 * x56 * x75))
    result[2, 1] = numpy.sum(
        x179 * (x201 * x47 * x77 + x202 * x47 * x76 + x203 * x74 * x84)
    )
    result[2, 2] = numpy.sum(x179 * (x205 * x207 * x66 + x205 * x75 * x84 + x206 * x207))
    result[2, 3] = numpy.sum(
        x185 * (x100 * x201 * x82 + x106 * x201 * x95 + x202 * x82 * x95)
    )
    result[2, 4] = numpy.sum(
        x186 * (x106 * x205 * x76 + x205 * x77 * x82 + x206 * x76 * x82)
    )
    result[2, 5] = numpy.sum(x185 * (x106 * x210 * x75 + x210 * x212 * x66 + x211 * x212))
    result[2, 6] = numpy.sum(
        x179 * (x105 * x128 * x202 + x105 * x131 * x201 + x128 * x135 * x201)
    )
    result[2, 7] = numpy.sum(
        x186 * (x100 * x105 * x205 + x105 * x206 * x95 + x135 * x205 * x95)
    )
    result[2, 8] = numpy.sum(
        x186 * (x105 * x210 * x77 + x105 * x211 * x76 + x135 * x210 * x76)
    )
    result[2, 9] = numpy.sum(x179 * (x105 * x215 * x66 + x105 * x217 * x75 + x135 * x215))
    result[2, 10] = numpy.sum(
        x174 * (x104 * x153 * x202 + x134 * x153 * x201 + x159 * x198 * x69)
    )
    result[2, 11] = numpy.sum(
        x179 * (x104 * x128 * x206 + x104 * x131 * x205 + x128 * x134 * x205)
    )
    result[2, 12] = numpy.sum(
        x185 * (x100 * x104 * x210 + x104 * x211 * x95 + x134 * x210 * x95)
    )
    result[2, 13] = numpy.sum(
        x179 * (x104 * x214 * x77 + x134 * x214 * x76 + x200 * x217 * x74)
    )
    result[2, 14] = numpy.sum(x174 * (x134 * x219 * x75 + x200 * x219 * x66 + x222 * x5))
    result[3, 0] = numpy.sum(x73 * (x223 * x225 * x71 + x223 * x54 * x80 + x224 * x225))
    result[3, 1] = numpy.sum(x87 * (x227 * x37 * x71 + x227 * x40 + x228 * x37 * x80))
    result[3, 2] = numpy.sum(
        x87 * (x223 * x37 * x90 + x223 * x40 * x89 + x224 * x37 * x89)
    )
    result[3, 3] = numpy.sum(x108 * (x230 * x28 * x71 + x230 * x32 + x231 * x28 * x80))
    result[3, 4] = numpy.sum(
        x110 * (x226 * x28 * x90 + x226 * x32 * x89 + x228 * x28 * x89)
    )
    result[3, 5] = numpy.sum(
        x108 * (x115 * x223 * x32 + x115 * x224 * x28 + x120 * x223 * x28)
    )
    result[3, 6] = numpy.sum(x87 * (x10 * x233 + x232 * x234 * x71 + x232 * x25 * x80))
    result[3, 7] = numpy.sum(
        x110 * (x18 * x229 * x90 + x229 * x25 * x89 + x231 * x234 * x88)
    )
    result[3, 8] = numpy.sum(
        x110 * (x115 * x18 * x228 + x115 * x226 * x25 + x120 * x18 * x226)
    )
    result[3, 9] = numpy.sum(
        x87 * (x143 * x18 * x224 + x143 * x223 * x25 + x146 * x18 * x223)
    )
    result[3, 10] = numpy.sum(
        x73
        * (
            x196
            * (
                x195 * x64
                + x2 * (5.0 * x156 + x158 + x188 * x48 + 4.0 * x190)
                - 2.0 * x23 * (-ax * x235 + x151 + x152)
            )
            + x236 * x71
            + x236 * x9
        )
    )
    result[3, 11] = numpy.sum(x87 * (x15 * x232 * x90 + x232 * x237 * x88 + x233 * x88))
    result[3, 12] = numpy.sum(x108 * (x115 * x15 * x231 + x115 * x238 * x9 + x120 * x238))
    result[3, 13] = numpy.sum(x87 * (x143 * x15 * x228 + x143 * x239 * x9 + x146 * x239))
    result[3, 14] = numpy.sum(x73 * (x15 * x169 * x223 + x223 * x240 * x9 + x224 * x240))
    result[4, 0] = numpy.sum(
        x174 * (x170 * x202 * x50 + x171 * x201 * x50 + x203 * x54 * x64)
    )
    result[4, 1] = numpy.sum(
        x179 * (x176 * x201 * x40 + x176 * x202 * x37 + x177 * x201 * x37)
    )
    result[4, 2] = numpy.sum(
        x179 * (x170 * x205 * x40 + x170 * x206 * x37 + x171 * x205 * x37)
    )
    result[4, 3] = numpy.sum(
        x185 * (x182 * x201 * x32 + x182 * x202 * x28 + x183 * x201 * x28)
    )
    result[4, 4] = numpy.sum(
        x186 * (x176 * x205 * x32 + x176 * x206 * x28 + x177 * x205 * x28)
    )
    result[4, 5] = numpy.sum(
        x185 * (x170 * x210 * x32 + x170 * x211 * x28 + x171 * x210 * x28)
    )
    result[4, 6] = numpy.sum(
        x179 * (x18 * x188 * x202 + x188 * x201 * x25 + x191 * x234 * x69)
    )
    result[4, 7] = numpy.sum(
        x186 * (x18 * x182 * x206 + x18 * x183 * x205 + x182 * x205 * x25)
    )
    result[4, 8] = numpy.sum(
        x186 * (x176 * x18 * x211 + x176 * x210 * x25 + x177 * x18 * x210)
    )
    result[4, 9] = numpy.sum(
        x179 * (x10 * x217 * x241 + x170 * x214 * x25 + x171 * x18 * x214)
    )
    result[4, 10] = numpy.sum(x174 * (x15 * x193 * x202 + x193 * x237 * x69 + x197 * x69))
    result[4, 11] = numpy.sum(x179 * (x15 * x191 * x205 + x205 * x242 * x9 + x206 * x242))
    result[4, 12] = numpy.sum(x185 * (x15 * x183 * x210 + x210 * x243 * x9 + x211 * x243))
    result[4, 13] = numpy.sum(x179 * (x15 * x176 * x217 + x176 * x244 * x9 + x177 * x244))
    result[4, 14] = numpy.sum(x174 * (x15 * x171 * x219 + x219 * x241 * x9 + x222 * x64))
    result[5, 0] = numpy.sum(x73 * (x245 * x247 * x66 + x245 * x54 * x75 + x246 * x247))
    result[5, 1] = numpy.sum(
        x87 * (x245 * x37 * x77 + x245 * x40 * x76 + x246 * x37 * x76)
    )
    result[5, 2] = numpy.sum(x87 * (x249 * x37 * x66 + x249 * x40 + x250 * x37 * x75))
    result[5, 3] = numpy.sum(
        x108 * (x100 * x245 * x28 + x245 * x32 * x95 + x246 * x28 * x95)
    )
    result[5, 4] = numpy.sum(
        x110 * (x248 * x28 * x77 + x248 * x32 * x76 + x250 * x28 * x76)
    )
    result[5, 5] = numpy.sum(x108 * (x252 * x28 * x66 + x252 * x32 + x253 * x28 * x75))
    result[5, 6] = numpy.sum(
        x87 * (x128 * x18 * x246 + x128 * x245 * x25 + x131 * x18 * x245)
    )
    result[5, 7] = numpy.sum(
        x110 * (x100 * x18 * x248 + x18 * x250 * x95 + x248 * x25 * x95)
    )
    result[5, 8] = numpy.sum(
        x110 * (x18 * x251 * x77 + x25 * x251 * x76 + x253 * x254 * x74)
    )
    result[5, 9] = numpy.sum(x87 * (x10 * x256 + x25 * x255 * x75 + x254 * x255 * x66))
    result[5, 10] = numpy.sum(x73 * (x15 * x159 * x245 + x245 * x257 * x9 + x246 * x257))
    result[5, 11] = numpy.sum(x87 * (x128 * x15 * x250 + x128 * x258 * x9 + x131 * x258))
    result[5, 12] = numpy.sum(x108 * (x100 * x259 + x15 * x253 * x95 + x259 * x9 * x95))
    result[5, 13] = numpy.sum(
        x87 * (x15 * x255 * x77 + x199 * x255 * x74 * x9 + x256 * x74)
    )
    result[5, 14] = numpy.sum(
        x73
        * (
            x199
            * (
                x2 * (5.0 * x166 + x168 + x214 * x48 + 4.0 * x216)
                + x221 * x69
                - 2.0 * x23 * (-ax * x260 + x161 + x162)
            )
            + x261 * x66
            + x261 * x9
        )
    )
    return result


def kinetic3d_30(ax, da, A, bx, db, B):
    """Cartesian 3D (fs) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 1), dtype=float)

    x0 = 2.0 * ax
    x1 = (2.0 * bx + x0) ** (-1.0)
    x2 = (ax + bx) ** (-1.0)
    x3 = x2 * (ax * A[0] + bx * B[0]) - A[0]
    x4 = -ax
    x5 = x3**2
    x6 = 2.0 * ax**2
    x7 = -x4 - x6 * (x1 + x5)
    x8 = bx * x2
    x9 = ax * x8
    x10 = numpy.exp(-x9 * (A[0] - B[0]) ** 2)
    x11 = 1.77245385090552 * numpy.sqrt(x2)
    x12 = x10 * x11
    x13 = 2.0 * x12
    x14 = x12 * x3
    x15 = 4.0 * x9
    x16 = x0 * x8
    x17 = x14 * (x16 + x7)
    x18 = x1 * x12
    x19 = x12 * x5 + x18
    x20 = x17 * x3 + x18 * x7 + x8 * (x0 * x19 - x13)
    x21 = x3 * (2.0 * x18 + x19)
    x22 = numpy.exp(-x9 * (A[1] - B[1]) ** 2)
    x23 = numpy.exp(-x9 * (A[2] - B[2]) ** 2)
    x24 = 3.14159265358979 * x2 * x23
    x25 = x22 * x24
    x26 = x2 * (ax * A[1] + bx * B[1]) - A[1]
    x27 = x26**2
    x28 = -x4 - x6 * (x1 + x27)
    x29 = x21 * x25
    x30 = x2 * (ax * A[2] + bx * B[2]) - A[2]
    x31 = x30**2
    x32 = -x4 - x6 * (x1 + x31)
    x33 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**1.5)
    x34 = 5.84237394672177 * x33
    x35 = x11 * x22
    x36 = x26 * x35
    x37 = x36 * (x16 + x28)
    x38 = x11 * x23
    x39 = x20 * x25
    x40 = x19 * x25
    x41 = 13.0639452948436 * x33
    x42 = x30 * x38
    x43 = x42 * (x16 + x32)
    x44 = x1 * x35
    x45 = x27 * x35 + x44
    x46 = 2.0 * x35
    x47 = x26 * x37 + x28 * x44 + x8 * (x0 * x45 - x46)
    x48 = x10 * x24
    x49 = x47 * x48
    x50 = x3 * x48
    x51 = 3.14159265358979 * x10 * x2 * x22
    x52 = x3 * x51
    x53 = x1 * x38
    x54 = x31 * x38 + x53
    x55 = 2.0 * x38
    x56 = x30 * x43 + x32 * x53 + x8 * (x0 * x54 - x55)
    x57 = x51 * x56
    x58 = x26 * (2.0 * x44 + x45)
    x59 = x48 * x58
    x60 = x30 * (2.0 * x53 + x54)
    x61 = x51 * x60

    # 10 item(s)
    result[0, 0] = numpy.sum(
        x34
        * (
            x25
            * (x1 * (x13 * x3 * x7 + x14 * x15) + x20 * x3 + x8 * (x0 * x21 - 3.0 * x14))
            + x28 * x29
            + x29 * x32
        )
    )
    result[1, 0] = numpy.sum(x41 * (x19 * x37 * x38 + x26 * x32 * x40 + x26 * x39))
    result[2, 0] = numpy.sum(x41 * (x19 * x35 * x43 + x28 * x30 * x40 + x30 * x39))
    result[3, 0] = numpy.sum(x41 * (x17 * x38 * x45 + x3 * x49 + x32 * x45 * x50))
    result[4, 0] = numpy.sum(
        22.6274169979695
        * x33
        * (x17 * x25 * x26 * x30 + x26 * x43 * x52 + x30 * x37 * x50)
    )
    result[5, 0] = numpy.sum(x41 * (x17 * x35 * x54 + x28 * x52 * x54 + x3 * x57))
    result[6, 0] = numpy.sum(
        x34
        * (
            x32 * x59
            + x48
            * (
                x1 * (x15 * x36 + x26 * x28 * x46)
                + x26 * x47
                + x8 * (x0 * x58 - 3.0 * x36)
            )
            + x59 * x7
        )
    )
    result[7, 0] = numpy.sum(x41 * (x12 * x43 * x45 + x30 * x45 * x48 * x7 + x30 * x49))
    result[8, 0] = numpy.sum(x41 * (x12 * x37 * x54 + x26 * x51 * x54 * x7 + x26 * x57))
    result[9, 0] = numpy.sum(
        x34
        * (
            x28 * x61
            + x51
            * (
                x1 * (x15 * x42 + x30 * x32 * x55)
                + x30 * x56
                + x8 * (x0 * x60 - 3.0 * x42)
            )
            + x61 * x7
        )
    )
    return result


def kinetic3d_31(ax, da, A, bx, db, B):
    """Cartesian 3D (fp) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 3), dtype=float)

    x0 = 2.0 * ax
    x1 = (2.0 * bx + x0) ** (-1.0)
    x2 = (ax + bx) ** (-1.0)
    x3 = -x2 * (ax * A[0] + bx * B[0])
    x4 = -x3 - A[0]
    x5 = -ax
    x6 = x4**2
    x7 = 2.0 * ax**2
    x8 = -x5 - x7 * (x1 + x6)
    x9 = -x3 - B[0]
    x10 = bx * x2
    x11 = ax * x10
    x12 = numpy.exp(-x11 * (A[0] - B[0]) ** 2)
    x13 = 1.77245385090552 * numpy.sqrt(x2)
    x14 = x12 * x13
    x15 = x14 * x9
    x16 = x0 * x10
    x17 = x15 * (x16 + x8)
    x18 = x17 * x4
    x19 = x1 * x14
    x20 = x14 * x4
    x21 = x20 * x9
    x22 = x19 + x21
    x23 = 4.0 * x11
    x24 = x19 * x8
    x25 = x20 * (x16 + x8)
    x26 = x14 * x6
    x27 = x19 + x26
    x28 = 2.0 * x14
    x29 = x10 * (x0 * x27 - x28) + x25 * x4
    x30 = x16 * x22 + x18 + x24
    x31 = x1 * (x15 + x20) + x22 * x4
    x32 = x28 * x9
    x33 = x1 * (x17 + x25) + x10 * (x0 * x31 - x32) + x30 * x4
    x34 = 3.0 * x19
    x35 = x1 * (x26 + x32 * x4 + x34) + x31 * x4
    x36 = numpy.exp(-x11 * (A[1] - B[1]) ** 2)
    x37 = numpy.exp(-x11 * (A[2] - B[2]) ** 2)
    x38 = 3.14159265358979 * x2 * x37
    x39 = x36 * x38
    x40 = -x2 * (ax * A[1] + bx * B[1])
    x41 = -x40 - A[1]
    x42 = x41**2
    x43 = -x5 - x7 * (x1 + x42)
    x44 = x35 * x39
    x45 = -x2 * (ax * A[2] + bx * B[2])
    x46 = -x45 - A[2]
    x47 = x46**2
    x48 = -x5 - x7 * (x1 + x47)
    x49 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**2.5)
    x50 = 11.6847478934435 * x49
    x51 = -x40 - B[1]
    x52 = x13 * x36
    x53 = x51 * x52
    x54 = x53 * (x16 + x43)
    x55 = x4 * (2.0 * x19 + x27)
    x56 = x13 * x37
    x57 = x24 + x29
    x58 = x39 * (
        x1 * (x20 * x23 + x28 * x4 * x8) + x10 * (x0 * x55 - 3.0 * x20) + x4 * x57
    )
    x59 = x39 * x55
    x60 = -x45 - B[2]
    x61 = x56 * x60
    x62 = x61 * (x16 + x48)
    x63 = x41 * x52
    x64 = x63 * (x16 + x43)
    x65 = x33 * x39
    x66 = x39 * x41
    x67 = 26.1278905896872 * x49
    x68 = x41 * x54
    x69 = x1 * x13
    x70 = x36 * x69
    x71 = x51 * x63
    x72 = x70 + x71
    x73 = x43 * x70
    x74 = x16 * x72 + x68 + x73
    x75 = x27 * x56
    x76 = x46 * x56
    x77 = x76 * (x16 + x48)
    x78 = x39 * x46
    x79 = x46 * x62
    x80 = x37 * x69
    x81 = x60 * x76
    x82 = x80 + x81
    x83 = x48 * x80
    x84 = x16 * x82 + x79 + x83
    x85 = x27 * x52
    x86 = x42 * x52
    x87 = x70 + x86
    x88 = x56 * x87
    x89 = 2.0 * x52
    x90 = x10 * (x0 * x87 - x89) + x41 * x64
    x91 = x73 + x90
    x92 = x1 * (x53 + x63) + x41 * x72
    x93 = x51 * x89
    x94 = x1 * (x54 + x64) + x10 * (x0 * x92 - x93) + x41 * x74
    x95 = x12 * x38
    x96 = x94 * x95
    x97 = x4 * x95
    x98 = 45.2548339959391 * x49
    x99 = 3.14159265358979 * x12 * x2 * x36
    x100 = x4 * x99
    x101 = x47 * x56
    x102 = x101 + x80
    x103 = x102 * x52
    x104 = 2.0 * x56
    x105 = x10 * (x0 * x102 - x104) + x46 * x77
    x106 = x105 + x83
    x107 = x1 * (x61 + x76) + x46 * x82
    x108 = x104 * x60
    x109 = x1 * (x62 + x77) + x10 * (x0 * x107 - x108) + x46 * x84
    x110 = x109 * x99
    x111 = x41 * (2.0 * x70 + x87)
    x112 = x95 * (
        x1 * (x23 * x63 + x41 * x43 * x89) + x10 * (x0 * x111 - 3.0 * x63) + x41 * x91
    )
    x113 = x111 * x95
    x114 = 3.0 * x70
    x115 = x1 * (x114 + x41 * x93 + x86) + x41 * x92
    x116 = x115 * x95
    x117 = x46 * x95
    x118 = x14 * x87
    x119 = x41 * x99
    x120 = x102 * x14
    x121 = x46 * (x102 + 2.0 * x80)
    x122 = x99 * (
        x1 * (x104 * x46 * x48 + x23 * x76) + x10 * (x0 * x121 - 3.0 * x76) + x106 * x46
    )
    x123 = x121 * x99
    x124 = 3.0 * x80
    x125 = x1 * (x101 + x108 * x46 + x124) + x107 * x46
    x126 = x125 * x99

    # 30 item(s)
    result[0, 0] = numpy.sum(
        x50
        * (
            x39
            * (
                x1 * (2.0 * x18 + x22 * x23 + 3.0 * x24 + x29)
                - x10 * (-2.0 * ax * x35 + 3.0 * x21 + x34)
                + x33 * x4
            )
            + x43 * x44
            + x44 * x48
        )
    )
    result[0, 1] = numpy.sum(x50 * (x48 * x51 * x59 + x51 * x58 + x54 * x55 * x56))
    result[0, 2] = numpy.sum(x50 * (x43 * x59 * x60 + x52 * x55 * x62 + x58 * x60))
    result[1, 0] = numpy.sum(x67 * (x31 * x48 * x66 + x31 * x56 * x64 + x41 * x65))
    result[1, 1] = numpy.sum(x67 * (x48 * x72 * x75 + x56 * x57 * x72 + x74 * x75))
    result[1, 2] = numpy.sum(x67 * (x27 * x61 * x64 + x27 * x62 * x63 + x57 * x60 * x66))
    result[2, 0] = numpy.sum(x67 * (x31 * x43 * x78 + x31 * x52 * x77 + x46 * x65))
    result[2, 1] = numpy.sum(x67 * (x27 * x53 * x77 + x27 * x54 * x76 + x51 * x57 * x78))
    result[2, 2] = numpy.sum(x67 * (x43 * x82 * x85 + x52 * x57 * x82 + x84 * x85))
    result[3, 0] = numpy.sum(x67 * (x22 * x48 * x88 + x22 * x56 * x91 + x30 * x88))
    result[3, 1] = numpy.sum(x67 * (x25 * x56 * x92 + x4 * x96 + x48 * x92 * x97))
    result[3, 2] = numpy.sum(x67 * (x20 * x62 * x87 + x25 * x61 * x87 + x60 * x91 * x97))
    result[4, 0] = numpy.sum(x98 * (x22 * x63 * x77 + x22 * x64 * x76 + x30 * x46 * x66))
    result[4, 1] = numpy.sum(x98 * (x20 * x72 * x77 + x25 * x72 * x76 + x46 * x74 * x97))
    result[4, 2] = numpy.sum(x98 * (x100 * x41 * x84 + x20 * x64 * x82 + x25 * x63 * x82))
    result[5, 0] = numpy.sum(x67 * (x103 * x22 * x43 + x103 * x30 + x106 * x22 * x52))
    result[5, 1] = numpy.sum(
        x67 * (x100 * x106 * x51 + x102 * x20 * x54 + x102 * x25 * x53)
    )
    result[5, 2] = numpy.sum(x67 * (x100 * x107 * x43 + x107 * x25 * x52 + x110 * x4))
    result[6, 0] = numpy.sum(x50 * (x111 * x17 * x56 + x112 * x9 + x113 * x48 * x9))
    result[6, 1] = numpy.sum(
        x50
        * (
            x116 * x48
            + x116 * x8
            + x95
            * (
                x1 * (x23 * x72 + 2.0 * x68 + 3.0 * x73 + x90)
                - x10 * (-2.0 * ax * x115 + x114 + 3.0 * x71)
                + x41 * x94
            )
        )
    )
    result[6, 2] = numpy.sum(x50 * (x111 * x14 * x62 + x112 * x60 + x113 * x60 * x8))
    result[7, 0] = numpy.sum(x67 * (x117 * x9 * x91 + x15 * x77 * x87 + x17 * x76 * x87))
    result[7, 1] = numpy.sum(x67 * (x117 * x8 * x92 + x14 * x77 * x92 + x46 * x96))
    result[7, 2] = numpy.sum(x67 * (x118 * x8 * x82 + x118 * x84 + x14 * x82 * x91))
    result[8, 0] = numpy.sum(
        x67 * (x102 * x15 * x64 + x102 * x17 * x63 + x106 * x119 * x9)
    )
    result[8, 1] = numpy.sum(x67 * (x106 * x14 * x72 + x120 * x72 * x8 + x120 * x74))
    result[8, 2] = numpy.sum(x67 * (x107 * x119 * x8 + x107 * x14 * x64 + x110 * x41))
    result[9, 0] = numpy.sum(x50 * (x121 * x17 * x52 + x122 * x9 + x123 * x43 * x9))
    result[9, 1] = numpy.sum(x50 * (x121 * x14 * x54 + x122 * x51 + x123 * x51 * x8))
    result[9, 2] = numpy.sum(
        x50
        * (
            x126 * x43
            + x126 * x8
            + x99
            * (
                x1 * (x105 + x23 * x82 + 2.0 * x79 + 3.0 * x83)
                - x10 * (-2.0 * ax * x125 + x124 + 3.0 * x81)
                + x109 * x46
            )
        )
    )
    return result


def kinetic3d_32(ax, da, A, bx, db, B):
    """Cartesian 3D (fd) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 6), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = -ax
    x4 = (ax + bx) ** (-1.0)
    x5 = -x4 * (ax * A[0] + bx * B[0])
    x6 = -x5 - A[0]
    x7 = x6**2
    x8 = 2.0 * ax**2
    x9 = -x3 - x8 * (x2 + x7)
    x10 = -x5 - B[0]
    x11 = ax * x4
    x12 = bx * x11
    x13 = numpy.exp(-x12 * (A[0] - B[0]) ** 2)
    x14 = 1.77245385090552 * numpy.sqrt(x4)
    x15 = x13 * x14
    x16 = 2.0 * x15
    x17 = x10 * x16
    x18 = x10 * x15
    x19 = 4.0 * x12
    x20 = x2 * (x17 * x9 + x18 * x19)
    x21 = bx * x4
    x22 = x0 * x21
    x23 = x18 * (x22 + x9)
    x24 = x15 * x6
    x25 = x24 * (x22 + x9)
    x26 = x2 * (x23 + x25)
    x27 = x23 * x6
    x28 = x15 * x2
    x29 = x10 * x24
    x30 = x28 + x29
    x31 = x28 * x9
    x32 = x22 * x30 + x27 + x31
    x33 = x32 * x6
    x34 = x10**2 * x15
    x35 = x28 + x34
    x36 = -x16
    x37 = x10 * x23 + x11 * (x1 * x35 + x36)
    x38 = x31 + x37
    x39 = x38 * x6
    x40 = x2 * (x18 + x24)
    x41 = x30 * x6
    x42 = x40 + x41
    x43 = x0 * x42 - x17
    x44 = x1 * x4
    x45 = x35 * x6
    x46 = 2.0 * x28
    x47 = x10 * x46 + x45
    x48 = x19 * x30 + 2.0 * x27 + 3.0 * x31
    x49 = x20 + x22 * x47 + x39
    x50 = 3.0 * x28
    x51 = x17 * x6 + x50
    x52 = x2 * (x34 + x51) + x47 * x6
    x53 = x2 * (x37 + x48) - x21 * (-2.0 * ax * x52 + 2.0 * x34 + x46) + x49 * x6
    x54 = x10 * x28
    x55 = 2.0 * x2 * (x40 + x41 + x45 + 2.0 * x54) + x52 * x6
    x56 = numpy.exp(-x12 * (A[1] - B[1]) ** 2)
    x57 = numpy.exp(-x12 * (A[2] - B[2]) ** 2)
    x58 = 3.14159265358979 * x4 * x57
    x59 = x56 * x58
    x60 = -x4 * (ax * A[1] + bx * B[1])
    x61 = -x60 - A[1]
    x62 = x61**2
    x63 = -x3 - x8 * (x2 + x62)
    x64 = x55 * x59
    x65 = -x4 * (ax * A[2] + bx * B[2])
    x66 = -x65 - A[2]
    x67 = x66**2
    x68 = -x3 - x8 * (x2 + x67)
    x69 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**3.5)
    x70 = 13.4923846833851 * x69
    x71 = -x60 - B[1]
    x72 = x14 * x56
    x73 = x71 * x72
    x74 = x73 * (x22 + x63)
    x75 = x15 * x7
    x76 = x2 * (x51 + x75) + x42 * x6
    x77 = x14 * x57
    x78 = x28 + x75
    x79 = x21 * (x0 * x78 + x36) + x25 * x6
    x80 = x21 * x43 + x26 + x33
    x81 = x59 * (x2 * (x48 + x79) - x21 * (-2.0 * ax * x76 + 3.0 * x29 + x50) + x6 * x80)
    x82 = x59 * x76
    x83 = 23.3694957868871 * x69
    x84 = -x65 - B[2]
    x85 = x77 * x84
    x86 = x85 * (x22 + x68)
    x87 = x2 * x72
    x88 = x63 * x87
    x89 = x71**2 * x72
    x90 = x87 + x89
    x91 = 2.0 * x72
    x92 = -x91
    x93 = x11 * (x1 * x90 + x92) + x71 * x74
    x94 = x88 + x93
    x95 = x6 * (x46 + x78)
    x96 = x77 * x95
    x97 = x31 + x79
    x98 = x2 * (x16 * x6 * x9 + x19 * x24) + x21 * (x0 * x95 - 3.0 * x24) + x6 * x97
    x99 = x59 * x84
    x100 = x2 * x77
    x101 = x100 * x68
    x102 = x77 * x84**2
    x103 = x100 + x102
    x104 = 2.0 * x77
    x105 = -x104
    x106 = x11 * (x1 * x103 + x105) + x84 * x86
    x107 = x101 + x106
    x108 = x72 * x95
    x109 = x61 * x72
    x110 = x109 * (x22 + x63)
    x111 = x53 * x59
    x112 = x52 * x59
    x113 = 30.169889330626 * x69
    x114 = x61 * x74
    x115 = x109 * x71
    x116 = x115 + x87
    x117 = x114 + x116 * x22 + x88
    x118 = x42 * x77
    x119 = 52.2557811793745 * x69
    x120 = x61 * x90
    x121 = 2.0 * x87
    x122 = x120 + x121 * x71
    x123 = x122 * x77
    x124 = x71 * x91
    x125 = x2 * (x124 * x63 + x19 * x73)
    x126 = x61 * x94
    x127 = x122 * x22 + x125 + x126
    x128 = x66 * x77
    x129 = x128 * (x22 + x68)
    x130 = x59 * x66
    x131 = x66 * x86
    x132 = x128 * x84
    x133 = x100 + x132
    x134 = x101 + x131 + x133 * x22
    x135 = x42 * x72
    x136 = x103 * x66
    x137 = 2.0 * x100
    x138 = x136 + x137 * x84
    x139 = x138 * x72
    x140 = x104 * x84
    x141 = x2 * (x140 * x68 + x19 * x85)
    x142 = x107 * x66
    x143 = x138 * x22 + x141 + x142
    x144 = x62 * x72
    x145 = x144 + x87
    x146 = x110 * x61 + x21 * (x0 * x145 + x92)
    x147 = x146 + x88
    x148 = x47 * x77
    x149 = x2 * (x109 + x73)
    x150 = x116 * x61
    x151 = x149 + x150
    x152 = x151 * x77
    x153 = x2 * (x110 + x74)
    x154 = x117 * x61
    x155 = x0 * x151 - x124
    x156 = x153 + x154 + x155 * x21
    x157 = 3.0 * x87
    x158 = x124 * x61 + x157
    x159 = x122 * x61 + x2 * (x158 + x89)
    x160 = 2.0 * x114 + x116 * x19 + 3.0 * x88
    x161 = x127 * x61 + x2 * (x160 + x93) - x21 * (-2.0 * ax * x159 + x121 + 2.0 * x89)
    x162 = x13 * x58
    x163 = x161 * x162
    x164 = x162 * x6
    x165 = 90.5096679918781 * x69
    x166 = 3.14159265358979 * x13 * x4 * x56
    x167 = x166 * x6
    x168 = x67 * x77
    x169 = x100 + x168
    x170 = x129 * x66 + x21 * (x0 * x169 + x105)
    x171 = x101 + x170
    x172 = x47 * x72
    x173 = x2 * (x128 + x85)
    x174 = x133 * x66
    x175 = x173 + x174
    x176 = x175 * x72
    x177 = x2 * (x129 + x86)
    x178 = x134 * x66
    x179 = x0 * x175 - x140
    x180 = x177 + x178 + x179 * x21
    x181 = 3.0 * x100
    x182 = x140 * x66 + x181
    x183 = x138 * x66 + x2 * (x102 + x182)
    x184 = 3.0 * x101 + 2.0 * x131 + x133 * x19
    x185 = x143 * x66 + x2 * (x106 + x184) - x21 * (-2.0 * ax * x183 + 2.0 * x102 + x137)
    x186 = x166 * x185
    x187 = x61 * (x121 + x145)
    x188 = x187 * x77
    x189 = (
        x147 * x61 + x2 * (x109 * x19 + x61 * x63 * x91) + x21 * (x0 * x187 - 3.0 * x109)
    )
    x190 = x151 * x61 + x2 * (x144 + x158)
    x191 = x162 * (
        x156 * x61 + x2 * (x146 + x160) - x21 * (-2.0 * ax * x190 + 3.0 * x115 + x157)
    )
    x192 = x10 * x162
    x193 = x71 * x87
    x194 = x159 * x61 + 2.0 * x2 * (x120 + x149 + x150 + 2.0 * x193)
    x195 = x162 * x194
    x196 = x162 * x9
    x197 = x15 * x187
    x198 = x15 * x151
    x199 = x138 * x15
    x200 = x166 * x61
    x201 = x122 * x15
    x202 = x15 * x175
    x203 = x66 * (x137 + x169)
    x204 = x203 * x72
    x205 = (
        x171 * x66 + x2 * (x104 * x66 * x68 + x128 * x19) + x21 * (x0 * x203 - 3.0 * x128)
    )
    x206 = x10 * x166
    x207 = x175 * x66 + x2 * (x168 + x182)
    x208 = x166 * (
        x180 * x66 + x2 * (x170 + x184) - x21 * (-2.0 * ax * x207 + 3.0 * x132 + x181)
    )
    x209 = x15 * x203
    x210 = x100 * x84
    x211 = x183 * x66 + 2.0 * x2 * (x136 + x173 + x174 + 2.0 * x210)
    x212 = x166 * x211

    # 60 item(s)
    result[0, 0] = numpy.sum(
        x70
        * (
            x59
            * (
                x2
                * (x19 * x47 + 2.0 * x20 + 2.0 * x26 + 2.0 * x33 + 2.0 * x39 + x43 * x44)
                - x21 * (-2.0 * ax * x55 + 3.0 * x45 + 6.0 * x54)
                + x53 * x6
            )
            + x63 * x64
            + x64 * x68
        )
    )
    result[0, 1] = numpy.sum(x83 * (x68 * x71 * x82 + x71 * x81 + x74 * x76 * x77))
    result[0, 2] = numpy.sum(x83 * (x63 * x82 * x84 + x72 * x76 * x86 + x81 * x84))
    result[0, 3] = numpy.sum(x70 * (x68 * x90 * x96 + x77 * x90 * x98 + x94 * x96))
    result[0, 4] = numpy.sum(x83 * (x71 * x98 * x99 + x73 * x86 * x95 + x74 * x85 * x95))
    result[0, 5] = numpy.sum(x70 * (x103 * x108 * x63 + x103 * x72 * x98 + x107 * x108))
    result[1, 0] = numpy.sum(x113 * (x110 * x52 * x77 + x111 * x61 + x112 * x61 * x68))
    result[1, 1] = numpy.sum(x119 * (x116 * x118 * x68 + x116 * x77 * x80 + x117 * x118))
    result[1, 2] = numpy.sum(
        x119 * (x109 * x42 * x86 + x110 * x42 * x85 + x61 * x80 * x99)
    )
    result[1, 3] = numpy.sum(x113 * (x123 * x68 * x78 + x123 * x97 + x127 * x77 * x78))
    result[1, 4] = numpy.sum(
        x119 * (x116 * x78 * x86 + x116 * x85 * x97 + x117 * x78 * x85)
    )
    result[1, 5] = numpy.sum(
        x113 * (x103 * x109 * x97 + x103 * x110 * x78 + x107 * x109 * x78)
    )
    result[2, 0] = numpy.sum(x113 * (x111 * x66 + x112 * x63 * x66 + x129 * x52 * x72))
    result[2, 1] = numpy.sum(
        x119 * (x128 * x42 * x74 + x129 * x42 * x73 + x130 * x71 * x80)
    )
    result[2, 2] = numpy.sum(x119 * (x133 * x135 * x63 + x133 * x72 * x80 + x134 * x135))
    result[2, 3] = numpy.sum(
        x113 * (x128 * x78 * x94 + x128 * x90 * x97 + x129 * x78 * x90)
    )
    result[2, 4] = numpy.sum(
        x119 * (x133 * x73 * x97 + x133 * x74 * x78 + x134 * x73 * x78)
    )
    result[2, 5] = numpy.sum(x113 * (x139 * x63 * x78 + x139 * x97 + x143 * x72 * x78))
    result[3, 0] = numpy.sum(x113 * (x145 * x148 * x68 + x145 * x49 * x77 + x147 * x148))
    result[3, 1] = numpy.sum(x119 * (x152 * x30 * x68 + x152 * x32 + x156 * x30 * x77))
    result[3, 2] = numpy.sum(
        x119 * (x145 * x30 * x86 + x145 * x32 * x85 + x147 * x30 * x85)
    )
    result[3, 3] = numpy.sum(x113 * (x159 * x164 * x68 + x159 * x25 * x77 + x163 * x6))
    result[3, 4] = numpy.sum(
        x119 * (x151 * x24 * x86 + x151 * x25 * x85 + x156 * x164 * x84)
    )
    result[3, 5] = numpy.sum(
        x113 * (x103 * x145 * x25 + x103 * x147 * x24 + x107 * x145 * x24)
    )
    result[4, 0] = numpy.sum(
        x119 * (x109 * x129 * x47 + x110 * x128 * x47 + x130 * x49 * x61)
    )
    result[4, 1] = numpy.sum(
        x165 * (x116 * x128 * x32 + x116 * x129 * x30 + x117 * x128 * x30)
    )
    result[4, 2] = numpy.sum(
        x165 * (x109 * x133 * x32 + x109 * x134 * x30 + x110 * x133 * x30)
    )
    result[4, 3] = numpy.sum(
        x119 * (x122 * x128 * x25 + x122 * x129 * x24 + x127 * x164 * x66)
    )
    result[4, 4] = numpy.sum(
        x165 * (x116 * x133 * x25 + x116 * x134 * x24 + x117 * x133 * x24)
    )
    result[4, 5] = numpy.sum(
        x119 * (x109 * x138 * x25 + x110 * x138 * x24 + x143 * x167 * x61)
    )
    result[5, 0] = numpy.sum(x113 * (x169 * x172 * x63 + x169 * x49 * x72 + x171 * x172))
    result[5, 1] = numpy.sum(
        x119 * (x169 * x30 * x74 + x169 * x32 * x73 + x171 * x30 * x73)
    )
    result[5, 2] = numpy.sum(x119 * (x176 * x30 * x63 + x176 * x32 + x180 * x30 * x72))
    result[5, 3] = numpy.sum(
        x113 * (x169 * x24 * x94 + x169 * x25 * x90 + x171 * x24 * x90)
    )
    result[5, 4] = numpy.sum(
        x119 * (x167 * x180 * x71 + x175 * x24 * x74 + x175 * x25 * x73)
    )
    result[5, 5] = numpy.sum(x113 * (x167 * x183 * x63 + x183 * x25 * x72 + x186 * x6))
    result[6, 0] = numpy.sum(x70 * (x188 * x35 * x68 + x188 * x38 + x189 * x35 * x77))
    result[6, 1] = numpy.sum(x83 * (x10 * x191 + x190 * x192 * x68 + x190 * x23 * x77))
    result[6, 2] = numpy.sum(
        x83 * (x18 * x187 * x86 + x187 * x23 * x85 + x189 * x192 * x84)
    )
    result[6, 3] = numpy.sum(
        x70
        * (
            x162
            * (
                x161 * x61
                + x2
                * (
                    x122 * x19
                    + 2.0 * x125
                    + 2.0 * x126
                    + 2.0 * x153
                    + 2.0 * x154
                    + x155 * x44
                )
                - x21 * (-2.0 * ax * x194 + 3.0 * x120 + 6.0 * x193)
            )
            + x195 * x68
            + x195 * x9
        )
    )
    result[6, 4] = numpy.sum(x83 * (x15 * x190 * x86 + x190 * x196 * x84 + x191 * x84))
    result[6, 5] = numpy.sum(x70 * (x103 * x15 * x189 + x103 * x197 * x9 + x107 * x197))
    result[7, 0] = numpy.sum(
        x113 * (x128 * x145 * x38 + x128 * x147 * x35 + x129 * x145 * x35)
    )
    result[7, 1] = numpy.sum(
        x119 * (x128 * x151 * x23 + x129 * x151 * x18 + x156 * x192 * x66)
    )
    result[7, 2] = numpy.sum(
        x119 * (x133 * x145 * x23 + x133 * x147 * x18 + x134 * x145 * x18)
    )
    result[7, 3] = numpy.sum(x113 * (x129 * x15 * x159 + x159 * x196 * x66 + x163 * x66))
    result[7, 4] = numpy.sum(x119 * (x133 * x15 * x156 + x133 * x198 * x9 + x134 * x198))
    result[7, 5] = numpy.sum(x113 * (x143 * x145 * x15 + x145 * x199 * x9 + x147 * x199))
    result[8, 0] = numpy.sum(
        x113 * (x109 * x169 * x38 + x109 * x171 * x35 + x110 * x169 * x35)
    )
    result[8, 1] = numpy.sum(
        x119 * (x116 * x169 * x23 + x116 * x171 * x18 + x117 * x169 * x18)
    )
    result[8, 2] = numpy.sum(
        x119 * (x10 * x180 * x200 + x109 * x175 * x23 + x110 * x175 * x18)
    )
    result[8, 3] = numpy.sum(x113 * (x127 * x15 * x169 + x169 * x201 * x9 + x171 * x201))
    result[8, 4] = numpy.sum(x119 * (x116 * x15 * x180 + x116 * x202 * x9 + x117 * x202))
    result[8, 5] = numpy.sum(x113 * (x110 * x15 * x183 + x183 * x200 * x9 + x186 * x61))
    result[9, 0] = numpy.sum(x70 * (x204 * x35 * x63 + x204 * x38 + x205 * x35 * x72))
    result[9, 1] = numpy.sum(
        x83 * (x18 * x203 * x74 + x203 * x23 * x73 + x205 * x206 * x71)
    )
    result[9, 2] = numpy.sum(x83 * (x10 * x208 + x206 * x207 * x63 + x207 * x23 * x72))
    result[9, 3] = numpy.sum(x70 * (x15 * x205 * x90 + x209 * x9 * x90 + x209 * x94))
    result[9, 4] = numpy.sum(
        x83 * (x15 * x207 * x74 + x166 * x207 * x71 * x9 + x208 * x71)
    )
    result[9, 5] = numpy.sum(
        x70
        * (
            x166
            * (
                x185 * x66
                + x2
                * (
                    x138 * x19
                    + 2.0 * x141
                    + 2.0 * x142
                    + 2.0 * x177
                    + 2.0 * x178
                    + x179 * x44
                )
                - x21 * (-2.0 * ax * x211 + 3.0 * x136 + 6.0 * x210)
            )
            + x212 * x63
            + x212 * x9
        )
    )
    return result


def kinetic3d_33(ax, da, A, bx, db, B):
    """Cartesian 3D (ff) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 10), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - B[0]
    x6 = -ax
    x7 = -x4 - A[0]
    x8 = x7**2
    x9 = 2.0 * ax**2
    x10 = -x6 - x9 * (x2 + x8)
    x11 = ax * x3
    x12 = bx * x11
    x13 = numpy.exp(-x12 * (A[0] - B[0]) ** 2)
    x14 = 1.77245385090552 * numpy.sqrt(x3)
    x15 = x13 * x14
    x16 = x15 * x5
    x17 = bx * x3
    x18 = x0 * x17
    x19 = x16 * (x10 + x18)
    x20 = x19 * x5
    x21 = x15 * x5**2
    x22 = x15 * x2
    x23 = x21 + x22
    x24 = 2.0 * x15
    x25 = -x24
    x26 = x11 * (x1 * x23 + x25)
    x27 = x10 * x22
    x28 = 3.0 * x27
    x29 = x2 * (3.0 * x20 + 3.0 * x26 + x28)
    x30 = x24 * x5
    x31 = 4.0 * x12
    x32 = x2 * (x10 * x30 + x16 * x31)
    x33 = x20 + x26
    x34 = x27 + x33
    x35 = x23 * x5
    x36 = 2.0 * x22
    x37 = x36 * x5
    x38 = x35 + x37
    x39 = x11 * (x1 * x38 - 3.0 * x16) + x34 * x5
    x40 = x32 + x39
    x41 = x40 * x7
    x42 = x19 * x7
    x43 = x15 * x7
    x44 = x43 * x5
    x45 = x22 + x44
    x46 = x28 + x31 * x45 + 2.0 * x42
    x47 = x2 * (x33 + x46)
    x48 = x34 * x7
    x49 = x23 * x7
    x50 = x37 + x49
    x51 = x18 * x50 + x32 + x48
    x52 = x51 * x7
    x53 = 3.0 * x22
    x54 = x30 * x7 + x53
    x55 = x2 * (x21 + x54)
    x56 = x50 * x7
    x57 = x55 + x56
    x58 = x17 * (2.0 * ax * x57 - 2.0 * x21 - x36)
    x59 = x2 * (3.0 * x21 + x53)
    x60 = x38 * x7
    x61 = x59 + x60
    x62 = 6.0 * x12
    x63 = x18 * x61 + x29 + x41
    x64 = 3.0 * x49
    x65 = x22 * x5
    x66 = x2 * (x35 + x64 + 8.0 * x65) + x61 * x7
    x67 = 4.0 * x65
    x68 = (
        -x17 * (-2.0 * ax * x66 + 2.0 * x35 + x67)
        + x2 * (4.0 * x32 + x39 + 3.0 * x48 + x50 * x62)
        + x63 * x7
    )
    x69 = x2 * (3.0 * x55 + 3.0 * x56 + 2.0 * x59 + 2.0 * x60) + x66 * x7
    x70 = numpy.exp(-x12 * (A[1] - B[1]) ** 2)
    x71 = numpy.exp(-x12 * (A[2] - B[2]) ** 2)
    x72 = 3.14159265358979 * x3 * x71
    x73 = x70 * x72
    x74 = -x3 * (ax * A[1] + bx * B[1])
    x75 = -x74 - A[1]
    x76 = x75**2
    x77 = -x6 - x9 * (x2 + x76)
    x78 = x69 * x73
    x79 = -x3 * (ax * A[2] + bx * B[2])
    x80 = -x79 - A[2]
    x81 = x80**2
    x82 = -x6 - x9 * (x2 + x81)
    x83 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**4.5)
    x84 = 12.0679557322504 * x83
    x85 = -x74 - B[1]
    x86 = x14 * x70
    x87 = x85 * x86
    x88 = x87 * (x18 + x77)
    x89 = x2 * (x16 + x43)
    x90 = x45 * x7
    x91 = x2 * (2.0 * x49 + x67 + 2.0 * x89 + 2.0 * x90) + x57 * x7
    x92 = x14 * x71
    x93 = x43 * (x10 + x18)
    x94 = x2 * (x19 + x93)
    x95 = x18 * x45 + x27 + x42
    x96 = x7 * x95
    x97 = x89 + x90
    x98 = x0 * x97 - x30
    x99 = x1 * x3
    x100 = x47 + x52 + x58
    x101 = x73 * (
        x100 * x7
        - x17 * (-2.0 * ax * x91 + x64 + 6.0 * x65)
        + x2 * (x31 * x50 + 2.0 * x32 + 2.0 * x48 + 2.0 * x94 + 2.0 * x96 + x98 * x99)
    )
    x102 = x73 * x91
    x103 = 26.9847693667702 * x83
    x104 = -x79 - B[2]
    x105 = x104 * x92
    x106 = x105 * (x18 + x82)
    x107 = x2 * x86
    x108 = x107 * x77
    x109 = x85 * x88
    x110 = x85**2 * x86
    x111 = x107 + x110
    x112 = 2.0 * x86
    x113 = -x112
    x114 = x11 * (x1 * x111 + x113)
    x115 = x109 + x114
    x116 = x108 + x115
    x117 = x15 * x8
    x118 = x2 * (x117 + x54) + x7 * x97
    x119 = x118 * x92
    x120 = x117 + x22
    x121 = x17 * (x0 * x120 + x25) + x7 * x93
    x122 = x17 * x98 + x94 + x96
    x123 = x122 * x7 - x17 * (-2.0 * ax * x118 + 3.0 * x44 + x53) + x2 * (x121 + x46)
    x124 = x104 * x73
    x125 = 46.7389915737742 * x83
    x126 = x2 * x92
    x127 = x126 * x82
    x128 = x104 * x106
    x129 = x104**2 * x92
    x130 = x126 + x129
    x131 = 2.0 * x92
    x132 = -x131
    x133 = x11 * (x1 * x130 + x132)
    x134 = x128 + x133
    x135 = x127 + x134
    x136 = x118 * x86
    x137 = x112 * x85
    x138 = x2 * (x137 * x77 + x31 * x87)
    x139 = x111 * x85
    x140 = 2.0 * x107
    x141 = x140 * x85
    x142 = x139 + x141
    x143 = x11 * (x1 * x142 - 3.0 * x87) + x116 * x85
    x144 = x138 + x143
    x145 = x7 * (x120 + x36)
    x146 = x145 * x92
    x147 = x121 + x27
    x148 = x147 * x7 + x17 * (x0 * x145 - 3.0 * x43) + x2 * (x10 * x24 * x7 + x31 * x43)
    x149 = x104 * x131
    x150 = x2 * (x105 * x31 + x149 * x82)
    x151 = x104 * x130
    x152 = 2.0 * x126
    x153 = x104 * x152
    x154 = x151 + x153
    x155 = x104 * x135 + x11 * (x1 * x154 - 3.0 * x105)
    x156 = x150 + x155
    x157 = x145 * x86
    x158 = x75 * x86
    x159 = x158 * (x18 + x77)
    x160 = x68 * x73
    x161 = x66 * x73
    x162 = x75 * x88
    x163 = x158 * x85
    x164 = x107 + x163
    x165 = x108 + x162 + x164 * x18
    x166 = x57 * x92
    x167 = 60.3397786612521 * x83
    x168 = x116 * x75
    x169 = x111 * x75
    x170 = x141 + x169
    x171 = x138 + x168 + x170 * x18
    x172 = x92 * x97
    x173 = 104.511562358749 * x83
    x174 = 3.0 * x107
    x175 = x2 * (3.0 * x110 + x174)
    x176 = x142 * x75
    x177 = x175 + x176
    x178 = x177 * x92
    x179 = 3.0 * x108
    x180 = x2 * (3.0 * x109 + 3.0 * x114 + x179)
    x181 = x144 * x75
    x182 = x177 * x18 + x180 + x181
    x183 = x80 * x92
    x184 = x183 * (x18 + x82)
    x185 = x73 * x80
    x186 = x106 * x80
    x187 = x104 * x183
    x188 = x126 + x187
    x189 = x127 + x18 * x188 + x186
    x190 = x57 * x86
    x191 = x135 * x80
    x192 = x130 * x80
    x193 = x153 + x192
    x194 = x150 + x18 * x193 + x191
    x195 = x86 * x97
    x196 = 3.0 * x126
    x197 = x2 * (3.0 * x129 + x196)
    x198 = x154 * x80
    x199 = x197 + x198
    x200 = x199 * x86
    x201 = 3.0 * x127
    x202 = x2 * (3.0 * x128 + 3.0 * x133 + x201)
    x203 = x156 * x80
    x204 = x18 * x199 + x202 + x203
    x205 = x76 * x86
    x206 = x107 + x205
    x207 = x159 * x75 + x17 * (x0 * x206 + x113)
    x208 = x108 + x207
    x209 = x61 * x92
    x210 = x2 * (x158 + x87)
    x211 = x164 * x75
    x212 = x210 + x211
    x213 = x212 * x92
    x214 = x2 * (x159 + x88)
    x215 = x165 * x75
    x216 = x0 * x212 - x137
    x217 = x17 * x216 + x214 + x215
    x218 = x137 * x75 + x174
    x219 = x2 * (x110 + x218)
    x220 = x170 * x75
    x221 = x219 + x220
    x222 = x221 * x92
    x223 = 2.0 * x162 + x164 * x31 + x179
    x224 = x2 * (x115 + x223)
    x225 = x171 * x75
    x226 = x17 * (2.0 * ax * x221 - 2.0 * x110 - x140)
    x227 = x224 + x225 + x226
    x228 = 3.0 * x169
    x229 = x107 * x85
    x230 = x177 * x75 + x2 * (x139 + x228 + 8.0 * x229)
    x231 = 4.0 * x229
    x232 = (
        -x17 * (-2.0 * ax * x230 + 2.0 * x139 + x231)
        + x182 * x75
        + x2 * (4.0 * x138 + x143 + 3.0 * x168 + x170 * x62)
    )
    x233 = x13 * x72
    x234 = x232 * x233
    x235 = x233 * x7
    x236 = 3.14159265358979 * x13 * x3 * x70
    x237 = x236 * x7
    x238 = x81 * x92
    x239 = x126 + x238
    x240 = x17 * (x0 * x239 + x132) + x184 * x80
    x241 = x127 + x240
    x242 = x61 * x86
    x243 = x2 * (x105 + x183)
    x244 = x188 * x80
    x245 = x243 + x244
    x246 = x245 * x86
    x247 = x2 * (x106 + x184)
    x248 = x189 * x80
    x249 = x0 * x245 - x149
    x250 = x17 * x249 + x247 + x248
    x251 = x149 * x80 + x196
    x252 = x2 * (x129 + x251)
    x253 = x193 * x80
    x254 = x252 + x253
    x255 = x254 * x86
    x256 = 2.0 * x186 + x188 * x31 + x201
    x257 = x2 * (x134 + x256)
    x258 = x194 * x80
    x259 = x17 * (2.0 * ax * x254 - 2.0 * x129 - x152)
    x260 = x257 + x258 + x259
    x261 = 3.0 * x192
    x262 = x104 * x126
    x263 = x199 * x80 + x2 * (x151 + x261 + 8.0 * x262)
    x264 = 4.0 * x262
    x265 = (
        -x17 * (-2.0 * ax * x263 + 2.0 * x151 + x264)
        + x2 * (4.0 * x150 + x155 + 3.0 * x191 + x193 * x62)
        + x204 * x80
    )
    x266 = x236 * x265
    x267 = x75 * (x140 + x206)
    x268 = x267 * x92
    x269 = (
        x17 * (x0 * x267 - 3.0 * x158) + x2 * (x112 * x75 * x77 + x158 * x31) + x208 * x75
    )
    x270 = x2 * (x205 + x218) + x212 * x75
    x271 = x270 * x92
    x272 = -x17 * (-2.0 * ax * x270 + 3.0 * x163 + x174) + x2 * (x207 + x223) + x217 * x75
    x273 = x2 * (2.0 * x169 + 2.0 * x210 + 2.0 * x211 + x231) + x221 * x75
    x274 = x233 * (
        -x17 * (-2.0 * ax * x273 + x228 + 6.0 * x229)
        + x2
        * (2.0 * x138 + 2.0 * x168 + x170 * x31 + 2.0 * x214 + 2.0 * x215 + x216 * x99)
        + x227 * x75
    )
    x275 = x233 * x5
    x276 = x2 * (2.0 * x175 + 2.0 * x176 + 3.0 * x219 + 3.0 * x220) + x230 * x75
    x277 = x233 * x276
    x278 = x10 * x233
    x279 = x15 * x270
    x280 = x15 * x267
    x281 = x15 * x221
    x282 = x15 * x212
    x283 = x15 * x199
    x284 = x236 * x75
    x285 = x15 * x177
    x286 = x15 * x245
    x287 = x15 * x254
    x288 = x80 * (x152 + x239)
    x289 = x288 * x86
    x290 = (
        x17 * (x0 * x288 - 3.0 * x183) + x2 * (x131 * x80 * x82 + x183 * x31) + x241 * x80
    )
    x291 = x2 * (x238 + x251) + x245 * x80
    x292 = x291 * x86
    x293 = -x17 * (-2.0 * ax * x291 + 3.0 * x187 + x196) + x2 * (x240 + x256) + x250 * x80
    x294 = x236 * x5
    x295 = x2 * (2.0 * x192 + 2.0 * x243 + 2.0 * x244 + x264) + x254 * x80
    x296 = x236 * (
        -x17 * (-2.0 * ax * x295 + x261 + 6.0 * x262)
        + x2
        * (2.0 * x150 + 2.0 * x191 + x193 * x31 + 2.0 * x247 + 2.0 * x248 + x249 * x99)
        + x260 * x80
    )
    x297 = x15 * x288
    x298 = x15 * x291
    x299 = x2 * (2.0 * x197 + 2.0 * x198 + 3.0 * x252 + 3.0 * x253) + x263 * x80
    x300 = x236 * x299

    # 100 item(s)
    result[0, 0] = numpy.sum(
        x84
        * (
            x73
            * (
                -x17 * (-2.0 * ax * x69 + 3.0 * x59 + 3.0 * x60)
                + x2
                * (2.0 * x29 + x31 * x61 + 2.0 * x41 + 3.0 * x47 + 3.0 * x52 + 3.0 * x58)
                + x68 * x7
            )
            + x77 * x78
            + x78 * x82
        )
    )
    result[0, 1] = numpy.sum(x103 * (x101 * x85 + x102 * x82 * x85 + x88 * x91 * x92))
    result[0, 2] = numpy.sum(x103 * (x101 * x104 + x102 * x104 * x77 + x106 * x86 * x91))
    result[0, 3] = numpy.sum(x103 * (x111 * x119 * x82 + x111 * x123 * x92 + x116 * x119))
    result[0, 4] = numpy.sum(
        x125 * (x105 * x118 * x88 + x106 * x118 * x87 + x123 * x124 * x85)
    )
    result[0, 5] = numpy.sum(x103 * (x123 * x130 * x86 + x130 * x136 * x77 + x135 * x136))
    result[0, 6] = numpy.sum(x84 * (x142 * x146 * x82 + x142 * x148 * x92 + x144 * x146))
    result[0, 7] = numpy.sum(
        x103 * (x105 * x111 * x148 + x105 * x116 * x145 + x106 * x111 * x145)
    )
    result[0, 8] = numpy.sum(
        x103 * (x130 * x145 * x88 + x130 * x148 * x87 + x135 * x145 * x87)
    )
    result[0, 9] = numpy.sum(x84 * (x148 * x154 * x86 + x154 * x157 * x77 + x156 * x157))
    result[1, 0] = numpy.sum(x103 * (x159 * x66 * x92 + x160 * x75 + x161 * x75 * x82))
    result[1, 1] = numpy.sum(x167 * (x100 * x164 * x92 + x164 * x166 * x82 + x165 * x166))
    result[1, 2] = numpy.sum(
        x167 * (x100 * x124 * x75 + x105 * x159 * x57 + x106 * x158 * x57)
    )
    result[1, 3] = numpy.sum(x167 * (x122 * x170 * x92 + x170 * x172 * x82 + x171 * x172))
    result[1, 4] = numpy.sum(
        x173 * (x105 * x122 * x164 + x105 * x165 * x97 + x106 * x164 * x97)
    )
    result[1, 5] = numpy.sum(
        x167 * (x122 * x130 * x158 + x130 * x159 * x97 + x135 * x158 * x97)
    )
    result[1, 6] = numpy.sum(x103 * (x120 * x178 * x82 + x120 * x182 * x92 + x147 * x178))
    result[1, 7] = numpy.sum(
        x167 * (x105 * x120 * x171 + x105 * x147 * x170 + x106 * x120 * x170)
    )
    result[1, 8] = numpy.sum(
        x167 * (x120 * x130 * x165 + x120 * x135 * x164 + x130 * x147 * x164)
    )
    result[1, 9] = numpy.sum(
        x103 * (x120 * x154 * x159 + x120 * x156 * x158 + x147 * x154 * x158)
    )
    result[2, 0] = numpy.sum(x103 * (x160 * x80 + x161 * x77 * x80 + x184 * x66 * x86))
    result[2, 1] = numpy.sum(
        x167 * (x100 * x185 * x85 + x183 * x57 * x88 + x184 * x57 * x87)
    )
    result[2, 2] = numpy.sum(x167 * (x100 * x188 * x86 + x188 * x190 * x77 + x189 * x190))
    result[2, 3] = numpy.sum(
        x167 * (x111 * x122 * x183 + x111 * x184 * x97 + x116 * x183 * x97)
    )
    result[2, 4] = numpy.sum(
        x173 * (x122 * x188 * x87 + x188 * x88 * x97 + x189 * x87 * x97)
    )
    result[2, 5] = numpy.sum(x167 * (x122 * x193 * x86 + x193 * x195 * x77 + x194 * x195))
    result[2, 6] = numpy.sum(
        x103 * (x120 * x142 * x184 + x120 * x144 * x183 + x142 * x147 * x183)
    )
    result[2, 7] = numpy.sum(
        x167 * (x111 * x120 * x189 + x111 * x147 * x188 + x116 * x120 * x188)
    )
    result[2, 8] = numpy.sum(
        x167 * (x120 * x193 * x88 + x120 * x194 * x87 + x147 * x193 * x87)
    )
    result[2, 9] = numpy.sum(x103 * (x120 * x200 * x77 + x120 * x204 * x86 + x147 * x200))
    result[3, 0] = numpy.sum(x103 * (x206 * x209 * x82 + x206 * x63 * x92 + x208 * x209))
    result[3, 1] = numpy.sum(x167 * (x213 * x50 * x82 + x213 * x51 + x217 * x50 * x92))
    result[3, 2] = numpy.sum(
        x167 * (x105 * x206 * x51 + x105 * x208 * x50 + x106 * x206 * x50)
    )
    result[3, 3] = numpy.sum(x167 * (x222 * x45 * x82 + x222 * x95 + x227 * x45 * x92))
    result[3, 4] = numpy.sum(
        x173 * (x105 * x212 * x95 + x105 * x217 * x45 + x106 * x212 * x45)
    )
    result[3, 5] = numpy.sum(
        x167 * (x130 * x206 * x95 + x130 * x208 * x45 + x135 * x206 * x45)
    )
    result[3, 6] = numpy.sum(x103 * (x230 * x235 * x82 + x230 * x92 * x93 + x234 * x7))
    result[3, 7] = numpy.sum(
        x167 * (x104 * x227 * x235 + x105 * x221 * x93 + x106 * x221 * x43)
    )
    result[3, 8] = numpy.sum(
        x167 * (x130 * x212 * x93 + x130 * x217 * x43 + x135 * x212 * x43)
    )
    result[3, 9] = numpy.sum(
        x103 * (x154 * x206 * x93 + x154 * x208 * x43 + x156 * x206 * x43)
    )
    result[4, 0] = numpy.sum(
        x125 * (x158 * x184 * x61 + x159 * x183 * x61 + x185 * x63 * x75)
    )
    result[4, 1] = numpy.sum(
        x173 * (x164 * x183 * x51 + x164 * x184 * x50 + x165 * x183 * x50)
    )
    result[4, 2] = numpy.sum(
        x173 * (x158 * x188 * x51 + x158 * x189 * x50 + x159 * x188 * x50)
    )
    result[4, 3] = numpy.sum(
        x173 * (x170 * x183 * x95 + x170 * x184 * x45 + x171 * x183 * x45)
    )
    result[4, 4] = numpy.sum(
        181.019335983756
        * x83
        * (x164 * x188 * x95 + x164 * x189 * x45 + x165 * x188 * x45)
    )
    result[4, 5] = numpy.sum(
        x173 * (x158 * x193 * x95 + x158 * x194 * x45 + x159 * x193 * x45)
    )
    result[4, 6] = numpy.sum(
        x125 * (x177 * x183 * x93 + x177 * x184 * x43 + x182 * x235 * x80)
    )
    result[4, 7] = numpy.sum(
        x173 * (x170 * x188 * x93 + x170 * x189 * x43 + x171 * x188 * x43)
    )
    result[4, 8] = numpy.sum(
        x173 * (x164 * x193 * x93 + x164 * x194 * x43 + x165 * x193 * x43)
    )
    result[4, 9] = numpy.sum(
        x125 * (x158 * x199 * x93 + x159 * x199 * x43 + x204 * x237 * x75)
    )
    result[5, 0] = numpy.sum(x103 * (x239 * x242 * x77 + x239 * x63 * x86 + x241 * x242))
    result[5, 1] = numpy.sum(
        x167 * (x239 * x50 * x88 + x239 * x51 * x87 + x241 * x50 * x87)
    )
    result[5, 2] = numpy.sum(x167 * (x246 * x50 * x77 + x246 * x51 + x250 * x50 * x86))
    result[5, 3] = numpy.sum(
        x167 * (x111 * x239 * x95 + x111 * x241 * x45 + x116 * x239 * x45)
    )
    result[5, 4] = numpy.sum(
        x173 * (x245 * x45 * x88 + x245 * x87 * x95 + x250 * x45 * x87)
    )
    result[5, 5] = numpy.sum(x167 * (x255 * x45 * x77 + x255 * x95 + x260 * x45 * x86))
    result[5, 6] = numpy.sum(
        x103 * (x142 * x239 * x93 + x142 * x241 * x43 + x144 * x239 * x43)
    )
    result[5, 7] = numpy.sum(
        x167 * (x111 * x245 * x93 + x111 * x250 * x43 + x116 * x245 * x43)
    )
    result[5, 8] = numpy.sum(
        x167 * (x237 * x260 * x85 + x254 * x43 * x88 + x254 * x87 * x93)
    )
    result[5, 9] = numpy.sum(x103 * (x237 * x263 * x77 + x263 * x86 * x93 + x266 * x7))
    result[6, 0] = numpy.sum(x84 * (x268 * x38 * x82 + x268 * x40 + x269 * x38 * x92))
    result[6, 1] = numpy.sum(x103 * (x23 * x271 * x82 + x23 * x272 * x92 + x271 * x34))
    result[6, 2] = numpy.sum(
        x103 * (x105 * x23 * x269 + x105 * x267 * x34 + x106 * x23 * x267)
    )
    result[6, 3] = numpy.sum(x103 * (x19 * x273 * x92 + x273 * x275 * x82 + x274 * x5))
    result[6, 4] = numpy.sum(
        x125 * (x104 * x272 * x275 + x105 * x19 * x270 + x106 * x16 * x270)
    )
    result[6, 5] = numpy.sum(
        x103 * (x130 * x16 * x269 + x130 * x19 * x267 + x135 * x16 * x267)
    )
    result[6, 6] = numpy.sum(
        x84
        * (
            x10 * x277
            + x233
            * (
                -x17 * (-2.0 * ax * x276 + 3.0 * x175 + 3.0 * x176)
                + x2
                * (
                    x177 * x31
                    + 2.0 * x180
                    + 2.0 * x181
                    + 3.0 * x224
                    + 3.0 * x225
                    + 3.0 * x226
                )
                + x232 * x75
            )
            + x277 * x82
        )
    )
    result[6, 7] = numpy.sum(
        x103 * (x104 * x273 * x278 + x104 * x274 + x106 * x15 * x273)
    )
    result[6, 8] = numpy.sum(x103 * (x10 * x130 * x279 + x130 * x15 * x272 + x135 * x279))
    result[6, 9] = numpy.sum(x84 * (x10 * x154 * x280 + x15 * x154 * x269 + x156 * x280))
    result[7, 0] = numpy.sum(
        x103 * (x183 * x206 * x40 + x183 * x208 * x38 + x184 * x206 * x38)
    )
    result[7, 1] = numpy.sum(
        x167 * (x183 * x212 * x34 + x183 * x217 * x23 + x184 * x212 * x23)
    )
    result[7, 2] = numpy.sum(
        x167 * (x188 * x206 * x34 + x188 * x208 * x23 + x189 * x206 * x23)
    )
    result[7, 3] = numpy.sum(
        x167 * (x16 * x184 * x221 + x183 * x19 * x221 + x227 * x275 * x80)
    )
    result[7, 4] = numpy.sum(
        x173 * (x16 * x188 * x217 + x16 * x189 * x212 + x188 * x19 * x212)
    )
    result[7, 5] = numpy.sum(
        x167 * (x16 * x193 * x208 + x16 * x194 * x206 + x19 * x193 * x206)
    )
    result[7, 6] = numpy.sum(x103 * (x15 * x184 * x230 + x230 * x278 * x80 + x234 * x80))
    result[7, 7] = numpy.sum(x167 * (x10 * x188 * x281 + x15 * x188 * x227 + x189 * x281))
    result[7, 8] = numpy.sum(x167 * (x10 * x193 * x282 + x15 * x193 * x217 + x194 * x282))
    result[7, 9] = numpy.sum(x103 * (x10 * x206 * x283 + x15 * x204 * x206 + x208 * x283))
    result[8, 0] = numpy.sum(
        x103 * (x158 * x239 * x40 + x158 * x241 * x38 + x159 * x239 * x38)
    )
    result[8, 1] = numpy.sum(
        x167 * (x164 * x23 * x241 + x164 * x239 * x34 + x165 * x23 * x239)
    )
    result[8, 2] = numpy.sum(
        x167 * (x158 * x23 * x250 + x158 * x245 * x34 + x159 * x23 * x245)
    )
    result[8, 3] = numpy.sum(
        x167 * (x16 * x170 * x241 + x16 * x171 * x239 + x170 * x19 * x239)
    )
    result[8, 4] = numpy.sum(
        x173 * (x16 * x164 * x250 + x16 * x165 * x245 + x164 * x19 * x245)
    )
    result[8, 5] = numpy.sum(
        x167 * (x158 * x19 * x254 + x159 * x16 * x254 + x260 * x284 * x5)
    )
    result[8, 6] = numpy.sum(x103 * (x10 * x239 * x285 + x15 * x182 * x239 + x241 * x285))
    result[8, 7] = numpy.sum(x167 * (x10 * x170 * x286 + x15 * x170 * x250 + x171 * x286))
    result[8, 8] = numpy.sum(x167 * (x10 * x164 * x287 + x15 * x164 * x260 + x165 * x287))
    result[8, 9] = numpy.sum(x103 * (x10 * x263 * x284 + x15 * x159 * x263 + x266 * x75))
    result[9, 0] = numpy.sum(x84 * (x289 * x38 * x77 + x289 * x40 + x290 * x38 * x86))
    result[9, 1] = numpy.sum(
        x103 * (x23 * x288 * x88 + x23 * x290 * x87 + x288 * x34 * x87)
    )
    result[9, 2] = numpy.sum(x103 * (x23 * x292 * x77 + x23 * x293 * x86 + x292 * x34))
    result[9, 3] = numpy.sum(
        x103 * (x111 * x16 * x290 + x111 * x19 * x288 + x116 * x16 * x288)
    )
    result[9, 4] = numpy.sum(
        x125 * (x16 * x291 * x88 + x19 * x291 * x87 + x293 * x294 * x85)
    )
    result[9, 5] = numpy.sum(x103 * (x19 * x295 * x86 + x294 * x295 * x77 + x296 * x5))
    result[9, 6] = numpy.sum(x84 * (x10 * x142 * x297 + x142 * x15 * x290 + x144 * x297))
    result[9, 7] = numpy.sum(x103 * (x10 * x111 * x298 + x111 * x15 * x293 + x116 * x298))
    result[9, 8] = numpy.sum(
        x103 * (x10 * x236 * x295 * x85 + x15 * x295 * x88 + x296 * x85)
    )
    result[9, 9] = numpy.sum(
        x84
        * (
            x10 * x300
            + x236
            * (
                -x17 * (-2.0 * ax * x299 + 3.0 * x197 + 3.0 * x198)
                + x2
                * (
                    x199 * x31
                    + 2.0 * x202
                    + 2.0 * x203
                    + 3.0 * x257
                    + 3.0 * x258
                    + 3.0 * x259
                )
                + x265 * x80
            )
            + x300 * x77
        )
    )
    return result


def kinetic3d_34(ax, da, A, bx, db, B):
    """Cartesian 3D (fg) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((10, 15), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = -ax
    x4 = (ax + bx) ** (-1.0)
    x5 = -x4 * (ax * A[0] + bx * B[0])
    x6 = -x5 - A[0]
    x7 = x6**2
    x8 = 2.0 * ax**2
    x9 = -x3 - x8 * (x2 + x7)
    x10 = -x5 - B[0]
    x11 = ax * x4
    x12 = bx * x11
    x13 = numpy.exp(-x12 * (A[0] - B[0]) ** 2)
    x14 = 1.77245385090552 * numpy.sqrt(x4)
    x15 = x13 * x14
    x16 = 2.0 * x15
    x17 = x10 * x16
    x18 = x10 * x15
    x19 = 4.0 * x12
    x20 = x2 * (x17 * x9 + x18 * x19)
    x21 = 4.0 * x20
    x22 = x15 * x2
    x23 = x22 * x9
    x24 = bx * x4
    x25 = x0 * x24
    x26 = x18 * (x25 + x9)
    x27 = x10 * x26
    x28 = x10**2 * x15
    x29 = x22 + x28
    x30 = -x16
    x31 = x11 * (x1 * x29 + x30)
    x32 = x27 + x31
    x33 = x23 + x32
    x34 = x10 * x33
    x35 = x10 * x29
    x36 = 2.0 * x22
    x37 = x10 * x36
    x38 = x35 + x37
    x39 = x11 * (x1 * x38 - 3.0 * x18)
    x40 = x2 * (x21 + 4.0 * x34 + 4.0 * x39)
    x41 = 3.0 * x23
    x42 = x2 * (3.0 * x27 + 3.0 * x31 + x41)
    x43 = x34 + x39
    x44 = x20 + x43
    x45 = 3.0 * x22
    x46 = x2 * (3.0 * x28 + x45)
    x47 = x10 * x38
    x48 = x46 + x47
    x49 = 4.0 * x22
    x50 = x10 * x44 - x11 * (-2.0 * bx * x48 + 4.0 * x28 + x49)
    x51 = x42 + x50
    x52 = x51 * x6
    x53 = x33 * x6
    x54 = x29 * x6
    x55 = x37 + x54
    x56 = x12 * x55
    x57 = x2 * (x21 + x43 + 3.0 * x53 + 6.0 * x56)
    x58 = x44 * x6
    x59 = x38 * x6
    x60 = x46 + x59
    x61 = x25 * x60 + x42 + x58
    x62 = x6 * x61
    x63 = 3.0 * x54
    x64 = x10 * x22
    x65 = 8.0 * x64
    x66 = x2 * (x35 + x63 + x65)
    x67 = x6 * x60
    x68 = x66 + x67
    x69 = x10 * x49
    x70 = x24 * (2.0 * ax * x68 - 2.0 * x35 - x69)
    x71 = x2 * (4.0 * x35 + x65)
    x72 = x48 * x6
    x73 = x71 + x72
    x74 = 8.0 * x12
    x75 = x25 * x73 + x40 + x52
    x76 = x2 * (5.0 * x46 + x47 + 4.0 * x59) + x6 * x73
    x77 = 2.0 * x46
    x78 = (
        x2 * (5.0 * x42 + x50 + 4.0 * x58 + x60 * x74)
        - x24 * (-2.0 * ax * x76 + 2.0 * x47 + x77)
        + x6 * x75
    )
    x79 = 2.0 * x2 * (2.0 * x66 + 2.0 * x67 + x71 + x72) + x6 * x76
    x80 = numpy.exp(-x12 * (A[1] - B[1]) ** 2)
    x81 = numpy.exp(-x12 * (A[2] - B[2]) ** 2)
    x82 = 3.14159265358979 * x4 * x81
    x83 = x80 * x82
    x84 = -x4 * (ax * A[1] + bx * B[1])
    x85 = -x84 - A[1]
    x86 = x85**2
    x87 = -x3 - x8 * (x2 + x86)
    x88 = x79 * x83
    x89 = -x4 * (ax * A[2] + bx * B[2])
    x90 = -x89 - A[2]
    x91 = x90**2
    x92 = -x3 - x8 * (x2 + x91)
    x93 = 0.179587122125167 * da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**5.5)
    x94 = 9.12251705727742 * x93
    x95 = -x84 - B[1]
    x96 = x14 * x80
    x97 = x95 * x96
    x98 = x97 * (x25 + x87)
    x99 = x17 * x6 + x45
    x100 = x2 * (x28 + x99)
    x101 = x55 * x6
    x102 = x2 * (3.0 * x100 + 3.0 * x101 + 2.0 * x59 + x77) + x6 * x68
    x103 = x14 * x81
    x104 = x26 * x6
    x105 = x15 * x6
    x106 = x10 * x105
    x107 = x106 + x22
    x108 = 2.0 * x104 + x107 * x19 + x41
    x109 = x2 * (x108 + x32)
    x110 = x20 + x25 * x55 + x53
    x111 = x110 * x6
    x112 = x100 + x101
    x113 = x24 * (2.0 * ax * x112 - 2.0 * x28 - x36)
    x114 = x57 + x62 + x70
    x115 = x83 * (
        x114 * x6
        + x2 * (3.0 * x109 + 3.0 * x111 + 3.0 * x113 + x19 * x60 + 2.0 * x42 + 2.0 * x58)
        - x24 * (-2.0 * ax * x102 + 3.0 * x46 + 3.0 * x59)
    )
    x116 = x102 * x83
    x117 = 24.1359114645008 * x93
    x118 = -x89 - B[2]
    x119 = x103 * x118
    x120 = x119 * (x25 + x92)
    x121 = x2 * x96
    x122 = x121 * x87
    x123 = x95 * x98
    x124 = x95**2 * x96
    x125 = x121 + x124
    x126 = 2.0 * x96
    x127 = -x126
    x128 = x11 * (x1 * x125 + x127)
    x129 = x123 + x128
    x130 = x122 + x129
    x131 = x2 * (x105 + x18)
    x132 = x107 * x6
    x133 = x112 * x6 + x2 * (2.0 * x131 + 2.0 * x132 + 2.0 * x54 + x69)
    x134 = x103 * x133
    x135 = x105 * (x25 + x9)
    x136 = x2 * (x135 + x26)
    x137 = x104 + x107 * x25 + x23
    x138 = x137 * x6
    x139 = x131 + x132
    x140 = x0 * x139 - x17
    x141 = x1 * x4
    x142 = x109 + x111 + x113
    x143 = (
        x142 * x6
        + x2 * (2.0 * x136 + 2.0 * x138 + x140 * x141 + 2.0 * x20 + 2.0 * x53 + 4.0 * x56)
        - x24 * (-2.0 * ax * x133 + x63 + 6.0 * x64)
    )
    x144 = 31.1593277158494 * x93
    x145 = x118 * x83
    x146 = 53.9695387335403 * x93
    x147 = x103 * x2
    x148 = x147 * x92
    x149 = x118 * x120
    x150 = x103 * x118**2
    x151 = x147 + x150
    x152 = 2.0 * x103
    x153 = -x152
    x154 = x11 * (x1 * x151 + x153)
    x155 = x149 + x154
    x156 = x148 + x155
    x157 = x133 * x96
    x158 = x126 * x95
    x159 = x2 * (x158 * x87 + x19 * x97)
    x160 = x130 * x95
    x161 = x125 * x95
    x162 = 2.0 * x121
    x163 = x162 * x95
    x164 = x161 + x163
    x165 = x11 * (x1 * x164 - 3.0 * x97)
    x166 = x160 + x165
    x167 = x159 + x166
    x168 = x15 * x7
    x169 = x139 * x6 + x2 * (x168 + x99)
    x170 = x103 * x169
    x171 = x168 + x22
    x172 = x135 * x6 + x24 * (x0 * x171 + x30)
    x173 = x136 + x138 + x140 * x24
    x174 = x173 * x6 + x2 * (x108 + x172) - x24 * (-2.0 * ax * x169 + 3.0 * x106 + x45)
    x175 = x118 * x152
    x176 = x2 * (x119 * x19 + x175 * x92)
    x177 = x118 * x156
    x178 = x118 * x151
    x179 = 2.0 * x147
    x180 = x118 * x179
    x181 = x178 + x180
    x182 = x11 * (x1 * x181 - 3.0 * x119)
    x183 = x177 + x182
    x184 = x176 + x183
    x185 = x169 * x96
    x186 = x172 + x23
    x187 = x6 * (x171 + x36)
    x188 = x186 * x6 + x2 * (x105 * x19 + x16 * x6 * x9) + x24 * (x0 * x187 - 3.0 * x105)
    x189 = 3.0 * x121
    x190 = x2 * (3.0 * x124 + x189)
    x191 = x164 * x95
    x192 = x190 + x191
    x193 = x103 * x192
    x194 = 3.0 * x122
    x195 = x2 * (3.0 * x123 + 3.0 * x128 + x194)
    x196 = 4.0 * x121
    x197 = -x11 * (-2.0 * bx * x192 + 4.0 * x124 + x196) + x167 * x95
    x198 = x195 + x197
    x199 = 3.0 * x147
    x200 = x2 * (3.0 * x150 + x199)
    x201 = x118 * x181
    x202 = x200 + x201
    x203 = x202 * x96
    x204 = 3.0 * x148
    x205 = x2 * (3.0 * x149 + 3.0 * x154 + x204)
    x206 = 4.0 * x147
    x207 = -x11 * (-2.0 * bx * x202 + 4.0 * x150 + x206) + x118 * x184
    x208 = x205 + x207
    x209 = x85 * x96
    x210 = x209 * (x25 + x87)
    x211 = x78 * x83
    x212 = x76 * x83
    x213 = 20.3985682659737 * x93
    x214 = x85 * x98
    x215 = x209 * x95
    x216 = x121 + x215
    x217 = x122 + x214 + x216 * x25
    x218 = x103 * x68
    x219 = x130 * x85
    x220 = x125 * x85
    x221 = x163 + x220
    x222 = x159 + x219 + x221 * x25
    x223 = x103 * x112
    x224 = 69.6743749058326 * x93
    x225 = 120.679557322504 * x93
    x226 = x167 * x85
    x227 = x164 * x85
    x228 = x190 + x227
    x229 = x195 + x226 + x228 * x25
    x230 = x103 * x139
    x231 = x121 * x95
    x232 = 8.0 * x231
    x233 = x2 * (4.0 * x161 + x232)
    x234 = x192 * x85
    x235 = x233 + x234
    x236 = x103 * x235
    x237 = 4.0 * x159
    x238 = x2 * (4.0 * x160 + 4.0 * x165 + x237)
    x239 = x198 * x85
    x240 = x235 * x25 + x238 + x239
    x241 = x103 * x90
    x242 = x241 * (x25 + x92)
    x243 = x83 * x90
    x244 = x120 * x90
    x245 = x118 * x241
    x246 = x147 + x245
    x247 = x148 + x244 + x246 * x25
    x248 = x68 * x96
    x249 = x156 * x90
    x250 = x151 * x90
    x251 = x180 + x250
    x252 = x176 + x249 + x25 * x251
    x253 = x112 * x96
    x254 = x184 * x90
    x255 = x181 * x90
    x256 = x200 + x255
    x257 = x205 + x25 * x256 + x254
    x258 = x139 * x96
    x259 = x118 * x147
    x260 = 8.0 * x259
    x261 = x2 * (4.0 * x178 + x260)
    x262 = x202 * x90
    x263 = x261 + x262
    x264 = x263 * x96
    x265 = 4.0 * x176
    x266 = x2 * (4.0 * x177 + 4.0 * x182 + x265)
    x267 = x208 * x90
    x268 = x25 * x263 + x266 + x267
    x269 = x86 * x96
    x270 = x121 + x269
    x271 = x210 * x85 + x24 * (x0 * x270 + x127)
    x272 = x122 + x271
    x273 = x103 * x73
    x274 = x2 * (x209 + x97)
    x275 = x216 * x85
    x276 = x274 + x275
    x277 = x103 * x276
    x278 = x2 * (x210 + x98)
    x279 = x217 * x85
    x280 = x0 * x276 - x158
    x281 = x24 * x280 + x278 + x279
    x282 = x158 * x85 + x189
    x283 = x2 * (x124 + x282)
    x284 = x221 * x85
    x285 = x283 + x284
    x286 = x103 * x285
    x287 = x19 * x216 + x194 + 2.0 * x214
    x288 = x2 * (x129 + x287)
    x289 = x222 * x85
    x290 = x24 * (2.0 * ax * x285 - 2.0 * x124 - x162)
    x291 = x288 + x289 + x290
    x292 = 3.0 * x220
    x293 = x2 * (x161 + x232 + x292)
    x294 = x228 * x85
    x295 = x293 + x294
    x296 = x103 * x295
    x297 = 6.0 * x12
    x298 = x2 * (x166 + 3.0 * x219 + x221 * x297 + x237)
    x299 = x229 * x85
    x300 = x196 * x95
    x301 = x24 * (2.0 * ax * x295 - 2.0 * x161 - x300)
    x302 = x298 + x299 + x301
    x303 = x2 * (5.0 * x190 + x191 + 4.0 * x227) + x235 * x85
    x304 = 2.0 * x190
    x305 = (
        x2 * (5.0 * x195 + x197 + 4.0 * x226 + x228 * x74)
        - x24 * (-2.0 * ax * x303 + 2.0 * x191 + x304)
        + x240 * x85
    )
    x306 = x13 * x82
    x307 = x305 * x306
    x308 = x306 * x6
    x309 = 35.3313566383285 * x93
    x310 = 93.4779831475484 * x93
    x311 = 120.679557322504 * x93
    x312 = 209.023124717498 * x93
    x313 = 3.14159265358979 * x13 * x4 * x80
    x314 = x313 * x6
    x315 = x103 * x91
    x316 = x147 + x315
    x317 = x24 * (x0 * x316 + x153) + x242 * x90
    x318 = x148 + x317
    x319 = x73 * x96
    x320 = x2 * (x119 + x241)
    x321 = x246 * x90
    x322 = x320 + x321
    x323 = x322 * x96
    x324 = x2 * (x120 + x242)
    x325 = x247 * x90
    x326 = x0 * x322 - x175
    x327 = x24 * x326 + x324 + x325
    x328 = x175 * x90 + x199
    x329 = x2 * (x150 + x328)
    x330 = x251 * x90
    x331 = x329 + x330
    x332 = x331 * x96
    x333 = x19 * x246 + x204 + 2.0 * x244
    x334 = x2 * (x155 + x333)
    x335 = x252 * x90
    x336 = x24 * (2.0 * ax * x331 - 2.0 * x150 - x179)
    x337 = x334 + x335 + x336
    x338 = 3.0 * x250
    x339 = x2 * (x178 + x260 + x338)
    x340 = x256 * x90
    x341 = x339 + x340
    x342 = x341 * x96
    x343 = x2 * (x183 + 3.0 * x249 + x251 * x297 + x265)
    x344 = x257 * x90
    x345 = x118 * x206
    x346 = x24 * (2.0 * ax * x341 - 2.0 * x178 - x345)
    x347 = x343 + x344 + x346
    x348 = x2 * (5.0 * x200 + x201 + 4.0 * x255) + x263 * x90
    x349 = 2.0 * x200
    x350 = (
        x2 * (5.0 * x205 + x207 + 4.0 * x254 + x256 * x74)
        - x24 * (-2.0 * ax * x348 + 2.0 * x201 + x349)
        + x268 * x90
    )
    x351 = x313 * x350
    x352 = x85 * (x162 + x270)
    x353 = (
        x2 * (x126 * x85 * x87 + x19 * x209) + x24 * (x0 * x352 - 3.0 * x209) + x272 * x85
    )
    x354 = x103 * x48
    x355 = x2 * (x269 + x282) + x276 * x85
    x356 = x103 * x355
    x357 = x2 * (x271 + x287) - x24 * (-2.0 * ax * x355 + x189 + 3.0 * x215) + x281 * x85
    x358 = x2 * (2.0 * x220 + 2.0 * x274 + 2.0 * x275 + x300) + x285 * x85
    x359 = x103 * x358
    x360 = (
        x2
        * (x141 * x280 + 2.0 * x159 + x19 * x221 + 2.0 * x219 + 2.0 * x278 + 2.0 * x279)
        - x24 * (-2.0 * ax * x358 + 6.0 * x231 + x292)
        + x291 * x85
    )
    x361 = x2 * (2.0 * x227 + 3.0 * x283 + 3.0 * x284 + x304) + x295 * x85
    x362 = x306 * (
        x2 * (x19 * x228 + 2.0 * x195 + 2.0 * x226 + 3.0 * x288 + 3.0 * x289 + 3.0 * x290)
        - x24 * (-2.0 * ax * x361 + 3.0 * x190 + 3.0 * x227)
        + x302 * x85
    )
    x363 = x10 * x306
    x364 = 2.0 * x2 * (x233 + x234 + 2.0 * x293 + 2.0 * x294) + x303 * x85
    x365 = x306 * x364
    x366 = x306 * x9
    x367 = x15 * x358
    x368 = x15 * x355
    x369 = x15 * x202
    x370 = x15 * x295
    x371 = x15 * x285
    x372 = x15 * x276
    x373 = x15 * x263
    x374 = x313 * x85
    x375 = x15 * x235
    x376 = x15 * x322
    x377 = x15 * x331
    x378 = x15 * x341
    x379 = x90 * (x179 + x316)
    x380 = (
        x2 * (x152 * x90 * x92 + x19 * x241) + x24 * (x0 * x379 - 3.0 * x241) + x318 * x90
    )
    x381 = x48 * x96
    x382 = x2 * (x315 + x328) + x322 * x90
    x383 = x382 * x96
    x384 = x2 * (x317 + x333) - x24 * (-2.0 * ax * x382 + x199 + 3.0 * x245) + x327 * x90
    x385 = x2 * (2.0 * x250 + 2.0 * x320 + 2.0 * x321 + x345) + x331 * x90
    x386 = x385 * x96
    x387 = (
        x2
        * (x141 * x326 + 2.0 * x176 + x19 * x251 + 2.0 * x249 + 2.0 * x324 + 2.0 * x325)
        - x24 * (-2.0 * ax * x385 + 6.0 * x259 + x338)
        + x337 * x90
    )
    x388 = x10 * x313
    x389 = x2 * (2.0 * x255 + 3.0 * x329 + 3.0 * x330 + x349) + x341 * x90
    x390 = x313 * (
        x2 * (x19 * x256 + 2.0 * x205 + 2.0 * x254 + 3.0 * x334 + 3.0 * x335 + 3.0 * x336)
        - x24 * (-2.0 * ax * x389 + 3.0 * x200 + 3.0 * x255)
        + x347 * x90
    )
    x391 = x15 * x192
    x392 = x15 * x382
    x393 = x15 * x385
    x394 = 2.0 * x2 * (x261 + x262 + 2.0 * x339 + 2.0 * x340) + x348 * x90
    x395 = x313 * x394

    # 150 item(s)
    result[0, 0] = numpy.sum(
        x94
        * (
            x83
            * (
                x2
                * (x19 * x73 + 2.0 * x40 + 2.0 * x52 + 4.0 * x57 + 4.0 * x62 + 4.0 * x70)
                - x24 * (-2.0 * ax * x79 + 3.0 * x71 + 3.0 * x72)
                + x6 * x78
            )
            + x87 * x88
            + x88 * x92
        )
    )
    result[0, 1] = numpy.sum(x117 * (x102 * x103 * x98 + x115 * x95 + x116 * x92 * x95))
    result[0, 2] = numpy.sum(x117 * (x102 * x120 * x96 + x115 * x118 + x116 * x118 * x87))
    result[0, 3] = numpy.sum(
        x144 * (x103 * x125 * x143 + x125 * x134 * x92 + x130 * x134)
    )
    result[0, 4] = numpy.sum(
        x146 * (x119 * x133 * x98 + x120 * x133 * x97 + x143 * x145 * x95)
    )
    result[0, 5] = numpy.sum(x144 * (x143 * x151 * x96 + x151 * x157 * x87 + x156 * x157))
    result[0, 6] = numpy.sum(
        x117 * (x103 * x164 * x174 + x164 * x170 * x92 + x167 * x170)
    )
    result[0, 7] = numpy.sum(
        x146 * (x119 * x125 * x174 + x119 * x130 * x169 + x120 * x125 * x169)
    )
    result[0, 8] = numpy.sum(
        x146 * (x151 * x169 * x98 + x151 * x174 * x97 + x156 * x169 * x97)
    )
    result[0, 9] = numpy.sum(x117 * (x174 * x181 * x96 + x181 * x185 * x87 + x184 * x185))
    result[0, 10] = numpy.sum(
        x94 * (x103 * x187 * x198 + x187 * x193 * x92 + x188 * x193)
    )
    result[0, 11] = numpy.sum(
        x117 * (x119 * x164 * x188 + x119 * x167 * x187 + x120 * x164 * x187)
    )
    result[0, 12] = numpy.sum(
        x144 * (x125 * x151 * x188 + x125 * x156 * x187 + x130 * x151 * x187)
    )
    result[0, 13] = numpy.sum(
        x117 * (x181 * x187 * x98 + x181 * x188 * x97 + x184 * x187 * x97)
    )
    result[0, 14] = numpy.sum(x94 * (x187 * x203 * x87 + x187 * x208 * x96 + x188 * x203))
    result[1, 0] = numpy.sum(x213 * (x103 * x210 * x76 + x211 * x85 + x212 * x85 * x92))
    result[1, 1] = numpy.sum(
        x146 * (x103 * x114 * x216 + x216 * x218 * x92 + x217 * x218)
    )
    result[1, 2] = numpy.sum(
        x146 * (x114 * x145 * x85 + x119 * x210 * x68 + x120 * x209 * x68)
    )
    result[1, 3] = numpy.sum(
        x224 * (x103 * x142 * x221 + x221 * x223 * x92 + x222 * x223)
    )
    result[1, 4] = numpy.sum(
        x225 * (x112 * x119 * x217 + x112 * x120 * x216 + x119 * x142 * x216)
    )
    result[1, 5] = numpy.sum(
        x224 * (x112 * x151 * x210 + x112 * x156 * x209 + x142 * x151 * x209)
    )
    result[1, 6] = numpy.sum(
        x146 * (x103 * x173 * x228 + x228 * x230 * x92 + x229 * x230)
    )
    result[1, 7] = numpy.sum(
        x225 * (x119 * x139 * x222 + x119 * x173 * x221 + x120 * x139 * x221)
    )
    result[1, 8] = numpy.sum(
        x225 * (x139 * x151 * x217 + x139 * x156 * x216 + x151 * x173 * x216)
    )
    result[1, 9] = numpy.sum(
        x146 * (x139 * x181 * x210 + x139 * x184 * x209 + x173 * x181 * x209)
    )
    result[1, 10] = numpy.sum(
        x213 * (x103 * x171 * x240 + x171 * x236 * x92 + x186 * x236)
    )
    result[1, 11] = numpy.sum(
        x146 * (x119 * x171 * x229 + x119 * x186 * x228 + x120 * x171 * x228)
    )
    result[1, 12] = numpy.sum(
        x224 * (x151 * x171 * x222 + x151 * x186 * x221 + x156 * x171 * x221)
    )
    result[1, 13] = numpy.sum(
        x146 * (x171 * x181 * x217 + x171 * x184 * x216 + x181 * x186 * x216)
    )
    result[1, 14] = numpy.sum(
        x213 * (x171 * x202 * x210 + x171 * x208 * x209 + x186 * x202 * x209)
    )
    result[2, 0] = numpy.sum(x213 * (x211 * x90 + x212 * x87 * x90 + x242 * x76 * x96))
    result[2, 1] = numpy.sum(
        x146 * (x114 * x243 * x95 + x241 * x68 * x98 + x242 * x68 * x97)
    )
    result[2, 2] = numpy.sum(x146 * (x114 * x246 * x96 + x246 * x248 * x87 + x247 * x248))
    result[2, 3] = numpy.sum(
        x224 * (x112 * x125 * x242 + x112 * x130 * x241 + x125 * x142 * x241)
    )
    result[2, 4] = numpy.sum(
        x225 * (x112 * x246 * x98 + x112 * x247 * x97 + x142 * x246 * x97)
    )
    result[2, 5] = numpy.sum(x224 * (x142 * x251 * x96 + x251 * x253 * x87 + x252 * x253))
    result[2, 6] = numpy.sum(
        x146 * (x139 * x164 * x242 + x139 * x167 * x241 + x164 * x173 * x241)
    )
    result[2, 7] = numpy.sum(
        x225 * (x125 * x139 * x247 + x125 * x173 * x246 + x130 * x139 * x246)
    )
    result[2, 8] = numpy.sum(
        x225 * (x139 * x251 * x98 + x139 * x252 * x97 + x173 * x251 * x97)
    )
    result[2, 9] = numpy.sum(x146 * (x173 * x256 * x96 + x256 * x258 * x87 + x257 * x258))
    result[2, 10] = numpy.sum(
        x213 * (x171 * x192 * x242 + x171 * x198 * x241 + x186 * x192 * x241)
    )
    result[2, 11] = numpy.sum(
        x146 * (x164 * x171 * x247 + x164 * x186 * x246 + x167 * x171 * x246)
    )
    result[2, 12] = numpy.sum(
        x224 * (x125 * x171 * x252 + x125 * x186 * x251 + x130 * x171 * x251)
    )
    result[2, 13] = numpy.sum(
        x146 * (x171 * x256 * x98 + x171 * x257 * x97 + x186 * x256 * x97)
    )
    result[2, 14] = numpy.sum(
        x213 * (x171 * x264 * x87 + x171 * x268 * x96 + x186 * x264)
    )
    result[3, 0] = numpy.sum(x213 * (x103 * x270 * x75 + x270 * x273 * x92 + x272 * x273))
    result[3, 1] = numpy.sum(x146 * (x103 * x281 * x60 + x277 * x60 * x92 + x277 * x61))
    result[3, 2] = numpy.sum(
        x146 * (x119 * x270 * x61 + x119 * x272 * x60 + x120 * x270 * x60)
    )
    result[3, 3] = numpy.sum(x224 * (x103 * x291 * x55 + x110 * x286 + x286 * x55 * x92))
    result[3, 4] = numpy.sum(
        x225 * (x110 * x119 * x276 + x119 * x281 * x55 + x120 * x276 * x55)
    )
    result[3, 5] = numpy.sum(
        x224 * (x110 * x151 * x270 + x151 * x272 * x55 + x156 * x270 * x55)
    )
    result[3, 6] = numpy.sum(
        x146 * (x103 * x107 * x302 + x107 * x296 * x92 + x137 * x296)
    )
    result[3, 7] = numpy.sum(
        x225 * (x107 * x119 * x291 + x107 * x120 * x285 + x119 * x137 * x285)
    )
    result[3, 8] = numpy.sum(
        x225 * (x107 * x151 * x281 + x107 * x156 * x276 + x137 * x151 * x276)
    )
    result[3, 9] = numpy.sum(
        x146 * (x107 * x181 * x272 + x107 * x184 * x270 + x137 * x181 * x270)
    )
    result[3, 10] = numpy.sum(x213 * (x103 * x135 * x303 + x303 * x308 * x92 + x307 * x6))
    result[3, 11] = numpy.sum(
        x146 * (x105 * x120 * x295 + x118 * x302 * x308 + x119 * x135 * x295)
    )
    result[3, 12] = numpy.sum(
        x224 * (x105 * x151 * x291 + x105 * x156 * x285 + x135 * x151 * x285)
    )
    result[3, 13] = numpy.sum(
        x146 * (x105 * x181 * x281 + x105 * x184 * x276 + x135 * x181 * x276)
    )
    result[3, 14] = numpy.sum(
        x213 * (x105 * x202 * x272 + x105 * x208 * x270 + x135 * x202 * x270)
    )
    result[4, 0] = numpy.sum(
        x309 * (x209 * x242 * x73 + x210 * x241 * x73 + x243 * x75 * x85)
    )
    result[4, 1] = numpy.sum(
        x310 * (x216 * x241 * x61 + x216 * x242 * x60 + x217 * x241 * x60)
    )
    result[4, 2] = numpy.sum(
        x310 * (x209 * x246 * x61 + x209 * x247 * x60 + x210 * x246 * x60)
    )
    result[4, 3] = numpy.sum(
        x311 * (x110 * x221 * x241 + x221 * x242 * x55 + x222 * x241 * x55)
    )
    result[4, 4] = numpy.sum(
        x312 * (x110 * x216 * x246 + x216 * x247 * x55 + x217 * x246 * x55)
    )
    result[4, 5] = numpy.sum(
        x311 * (x110 * x209 * x251 + x209 * x252 * x55 + x210 * x251 * x55)
    )
    result[4, 6] = numpy.sum(
        x310 * (x107 * x228 * x242 + x107 * x229 * x241 + x137 * x228 * x241)
    )
    result[4, 7] = numpy.sum(
        x312 * (x107 * x221 * x247 + x107 * x222 * x246 + x137 * x221 * x246)
    )
    result[4, 8] = numpy.sum(
        x312 * (x107 * x216 * x252 + x107 * x217 * x251 + x137 * x216 * x251)
    )
    result[4, 9] = numpy.sum(
        x310 * (x107 * x209 * x257 + x107 * x210 * x256 + x137 * x209 * x256)
    )
    result[4, 10] = numpy.sum(
        x309 * (x105 * x235 * x242 + x135 * x235 * x241 + x240 * x308 * x90)
    )
    result[4, 11] = numpy.sum(
        x310 * (x105 * x228 * x247 + x105 * x229 * x246 + x135 * x228 * x246)
    )
    result[4, 12] = numpy.sum(
        x311 * (x105 * x221 * x252 + x105 * x222 * x251 + x135 * x221 * x251)
    )
    result[4, 13] = numpy.sum(
        x310 * (x105 * x216 * x257 + x105 * x217 * x256 + x135 * x216 * x256)
    )
    result[4, 14] = numpy.sum(
        x309 * (x105 * x210 * x263 + x135 * x209 * x263 + x268 * x314 * x85)
    )
    result[5, 0] = numpy.sum(x213 * (x316 * x319 * x87 + x316 * x75 * x96 + x318 * x319))
    result[5, 1] = numpy.sum(
        x146 * (x316 * x60 * x98 + x316 * x61 * x97 + x318 * x60 * x97)
    )
    result[5, 2] = numpy.sum(x146 * (x323 * x60 * x87 + x323 * x61 + x327 * x60 * x96))
    result[5, 3] = numpy.sum(
        x224 * (x110 * x125 * x316 + x125 * x318 * x55 + x130 * x316 * x55)
    )
    result[5, 4] = numpy.sum(
        x225 * (x110 * x322 * x97 + x322 * x55 * x98 + x327 * x55 * x97)
    )
    result[5, 5] = numpy.sum(x224 * (x110 * x332 + x332 * x55 * x87 + x337 * x55 * x96))
    result[5, 6] = numpy.sum(
        x146 * (x107 * x164 * x318 + x107 * x167 * x316 + x137 * x164 * x316)
    )
    result[5, 7] = numpy.sum(
        x225 * (x107 * x125 * x327 + x107 * x130 * x322 + x125 * x137 * x322)
    )
    result[5, 8] = numpy.sum(
        x225 * (x107 * x331 * x98 + x107 * x337 * x97 + x137 * x331 * x97)
    )
    result[5, 9] = numpy.sum(x146 * (x107 * x342 * x87 + x107 * x347 * x96 + x137 * x342))
    result[5, 10] = numpy.sum(
        x213 * (x105 * x192 * x318 + x105 * x198 * x316 + x135 * x192 * x316)
    )
    result[5, 11] = numpy.sum(
        x146 * (x105 * x164 * x327 + x105 * x167 * x322 + x135 * x164 * x322)
    )
    result[5, 12] = numpy.sum(
        x224 * (x105 * x125 * x337 + x105 * x130 * x331 + x125 * x135 * x331)
    )
    result[5, 13] = numpy.sum(
        x146 * (x105 * x341 * x98 + x135 * x341 * x97 + x314 * x347 * x95)
    )
    result[5, 14] = numpy.sum(x213 * (x135 * x348 * x96 + x314 * x348 * x87 + x351 * x6))
    result[6, 0] = numpy.sum(x94 * (x103 * x352 * x51 + x352 * x354 * x92 + x353 * x354))
    result[6, 1] = numpy.sum(x117 * (x103 * x357 * x38 + x356 * x38 * x92 + x356 * x44))
    result[6, 2] = numpy.sum(
        x117 * (x119 * x352 * x44 + x119 * x353 * x38 + x120 * x352 * x38)
    )
    result[6, 3] = numpy.sum(x144 * (x103 * x29 * x360 + x29 * x359 * x92 + x33 * x359))
    result[6, 4] = numpy.sum(
        x146 * (x119 * x29 * x357 + x119 * x33 * x355 + x120 * x29 * x355)
    )
    result[6, 5] = numpy.sum(
        x144 * (x151 * x29 * x353 + x151 * x33 * x352 + x156 * x29 * x352)
    )
    result[6, 6] = numpy.sum(x117 * (x10 * x362 + x103 * x26 * x361 + x361 * x363 * x92))
    result[6, 7] = numpy.sum(
        x146 * (x118 * x360 * x363 + x119 * x26 * x358 + x120 * x18 * x358)
    )
    result[6, 8] = numpy.sum(
        x146 * (x151 * x18 * x357 + x151 * x26 * x355 + x156 * x18 * x355)
    )
    result[6, 9] = numpy.sum(
        x117 * (x18 * x181 * x353 + x18 * x184 * x352 + x181 * x26 * x352)
    )
    result[6, 10] = numpy.sum(
        x94
        * (
            x306
            * (
                x2
                * (
                    x19 * x235
                    + 2.0 * x238
                    + 2.0 * x239
                    + 4.0 * x298
                    + 4.0 * x299
                    + 4.0 * x301
                )
                - x24 * (-2.0 * ax * x364 + 3.0 * x233 + 3.0 * x234)
                + x305 * x85
            )
            + x365 * x9
            + x365 * x92
        )
    )
    result[6, 11] = numpy.sum(
        x117 * (x118 * x361 * x366 + x118 * x362 + x120 * x15 * x361)
    )
    result[6, 12] = numpy.sum(x144 * (x15 * x151 * x360 + x151 * x367 * x9 + x156 * x367))
    result[6, 13] = numpy.sum(x117 * (x15 * x181 * x357 + x181 * x368 * x9 + x184 * x368))
    result[6, 14] = numpy.sum(x94 * (x15 * x208 * x352 + x352 * x369 * x9 + x353 * x369))
    result[7, 0] = numpy.sum(
        x213 * (x241 * x270 * x51 + x241 * x272 * x48 + x242 * x270 * x48)
    )
    result[7, 1] = numpy.sum(
        x146 * (x241 * x276 * x44 + x241 * x281 * x38 + x242 * x276 * x38)
    )
    result[7, 2] = numpy.sum(
        x146 * (x246 * x270 * x44 + x246 * x272 * x38 + x247 * x270 * x38)
    )
    result[7, 3] = numpy.sum(
        x224 * (x241 * x285 * x33 + x241 * x29 * x291 + x242 * x285 * x29)
    )
    result[7, 4] = numpy.sum(
        x225 * (x246 * x276 * x33 + x246 * x281 * x29 + x247 * x276 * x29)
    )
    result[7, 5] = numpy.sum(
        x224 * (x251 * x270 * x33 + x251 * x272 * x29 + x252 * x270 * x29)
    )
    result[7, 6] = numpy.sum(
        x146 * (x18 * x242 * x295 + x241 * x26 * x295 + x302 * x363 * x90)
    )
    result[7, 7] = numpy.sum(
        x225 * (x18 * x246 * x291 + x18 * x247 * x285 + x246 * x26 * x285)
    )
    result[7, 8] = numpy.sum(
        x225 * (x18 * x251 * x281 + x18 * x252 * x276 + x251 * x26 * x276)
    )
    result[7, 9] = numpy.sum(
        x146 * (x18 * x256 * x272 + x18 * x257 * x270 + x256 * x26 * x270)
    )
    result[7, 10] = numpy.sum(x213 * (x15 * x242 * x303 + x303 * x366 * x90 + x307 * x90))
    result[7, 11] = numpy.sum(x146 * (x15 * x246 * x302 + x246 * x370 * x9 + x247 * x370))
    result[7, 12] = numpy.sum(x224 * (x15 * x251 * x291 + x251 * x371 * x9 + x252 * x371))
    result[7, 13] = numpy.sum(x146 * (x15 * x256 * x281 + x256 * x372 * x9 + x257 * x372))
    result[7, 14] = numpy.sum(x213 * (x15 * x268 * x270 + x270 * x373 * x9 + x272 * x373))
    result[8, 0] = numpy.sum(
        x213 * (x209 * x316 * x51 + x209 * x318 * x48 + x210 * x316 * x48)
    )
    result[8, 1] = numpy.sum(
        x146 * (x216 * x316 * x44 + x216 * x318 * x38 + x217 * x316 * x38)
    )
    result[8, 2] = numpy.sum(
        x146 * (x209 * x322 * x44 + x209 * x327 * x38 + x210 * x322 * x38)
    )
    result[8, 3] = numpy.sum(
        x224 * (x221 * x29 * x318 + x221 * x316 * x33 + x222 * x29 * x316)
    )
    result[8, 4] = numpy.sum(
        x225 * (x216 * x29 * x327 + x216 * x322 * x33 + x217 * x29 * x322)
    )
    result[8, 5] = numpy.sum(
        x224 * (x209 * x29 * x337 + x209 * x33 * x331 + x210 * x29 * x331)
    )
    result[8, 6] = numpy.sum(
        x146 * (x18 * x228 * x318 + x18 * x229 * x316 + x228 * x26 * x316)
    )
    result[8, 7] = numpy.sum(
        x225 * (x18 * x221 * x327 + x18 * x222 * x322 + x221 * x26 * x322)
    )
    result[8, 8] = numpy.sum(
        x225 * (x18 * x216 * x337 + x18 * x217 * x331 + x216 * x26 * x331)
    )
    result[8, 9] = numpy.sum(
        x146 * (x10 * x347 * x374 + x18 * x210 * x341 + x209 * x26 * x341)
    )
    result[8, 10] = numpy.sum(x213 * (x15 * x240 * x316 + x316 * x375 * x9 + x318 * x375))
    result[8, 11] = numpy.sum(x146 * (x15 * x228 * x327 + x228 * x376 * x9 + x229 * x376))
    result[8, 12] = numpy.sum(x224 * (x15 * x221 * x337 + x221 * x377 * x9 + x222 * x377))
    result[8, 13] = numpy.sum(x146 * (x15 * x216 * x347 + x216 * x378 * x9 + x217 * x378))
    result[8, 14] = numpy.sum(x213 * (x15 * x210 * x348 + x348 * x374 * x9 + x351 * x85))
    result[9, 0] = numpy.sum(x94 * (x379 * x381 * x87 + x379 * x51 * x96 + x380 * x381))
    result[9, 1] = numpy.sum(
        x117 * (x379 * x38 * x98 + x379 * x44 * x97 + x38 * x380 * x97)
    )
    result[9, 2] = numpy.sum(x117 * (x38 * x383 * x87 + x38 * x384 * x96 + x383 * x44))
    result[9, 3] = numpy.sum(
        x144 * (x125 * x29 * x380 + x125 * x33 * x379 + x130 * x29 * x379)
    )
    result[9, 4] = numpy.sum(
        x146 * (x29 * x382 * x98 + x29 * x384 * x97 + x33 * x382 * x97)
    )
    result[9, 5] = numpy.sum(x144 * (x29 * x386 * x87 + x29 * x387 * x96 + x33 * x386))
    result[9, 6] = numpy.sum(
        x117 * (x164 * x18 * x380 + x164 * x26 * x379 + x167 * x18 * x379)
    )
    result[9, 7] = numpy.sum(
        x146 * (x125 * x18 * x384 + x125 * x26 * x382 + x130 * x18 * x382)
    )
    result[9, 8] = numpy.sum(
        x146 * (x18 * x385 * x98 + x26 * x385 * x97 + x387 * x388 * x95)
    )
    result[9, 9] = numpy.sum(x117 * (x10 * x390 + x26 * x389 * x96 + x388 * x389 * x87))
    result[9, 10] = numpy.sum(x94 * (x15 * x198 * x379 + x379 * x391 * x9 + x380 * x391))
    result[9, 11] = numpy.sum(x117 * (x15 * x164 * x384 + x164 * x392 * x9 + x167 * x392))
    result[9, 12] = numpy.sum(x144 * (x125 * x15 * x387 + x125 * x393 * x9 + x130 * x393))
    result[9, 13] = numpy.sum(
        x117 * (x15 * x389 * x98 + x313 * x389 * x9 * x95 + x390 * x95)
    )
    result[9, 14] = numpy.sum(
        x94
        * (
            x313
            * (
                x2
                * (
                    x19 * x263
                    + 2.0 * x266
                    + 2.0 * x267
                    + 4.0 * x343
                    + 4.0 * x344
                    + 4.0 * x346
                )
                - x24 * (-2.0 * ax * x394 + 3.0 * x261 + 3.0 * x262)
                + x350 * x90
            )
            + x395 * x87
            + x395 * x9
        )
    )
    return result


def kinetic3d_40(ax, da, A, bx, db, B):
    """Cartesian 3D (gs) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 1), dtype=float)

    x0 = 2.0 * ax
    x1 = (2.0 * bx + x0) ** (-1.0)
    x2 = (ax + bx) ** (-1.0)
    x3 = x2 * (ax * A[0] + bx * B[0]) - A[0]
    x4 = -ax
    x5 = x3**2
    x6 = 2.0 * ax**2
    x7 = -x4 - x6 * (x1 + x5)
    x8 = bx * x2
    x9 = ax * x8
    x10 = numpy.exp(-x9 * (A[0] - B[0]) ** 2)
    x11 = 1.77245385090552 * numpy.sqrt(x2)
    x12 = x10 * x11
    x13 = x12 * x3
    x14 = x0 * x8
    x15 = x13 * (x14 + x7)
    x16 = x15 * x3
    x17 = x12 * x5
    x18 = x1 * x12
    x19 = x17 + x18
    x20 = 2.0 * x12
    x21 = x8 * (x0 * x19 - x20)
    x22 = x18 * x7
    x23 = 4.0 * x9
    x24 = x16 + x21 + x22
    x25 = x3 * (2.0 * x18 + x19)
    x26 = x1 * (x13 * x23 + x20 * x3 * x7) + x24 * x3 + x8 * (x0 * x25 - 3.0 * x13)
    x27 = 3.0 * x1 * (x17 + x18) + x25 * x3
    x28 = numpy.exp(-x9 * (A[1] - B[1]) ** 2)
    x29 = numpy.exp(-x9 * (A[2] - B[2]) ** 2)
    x30 = 3.14159265358979 * x2 * x29
    x31 = x28 * x30
    x32 = x2 * (ax * A[1] + bx * B[1]) - A[1]
    x33 = x32**2
    x34 = -x4 - x6 * (x1 + x33)
    x35 = x27 * x31
    x36 = x2 * (ax * A[2] + bx * B[2]) - A[2]
    x37 = x36**2
    x38 = -x4 - x6 * (x1 + x37)
    x39 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**1.5)
    x40 = 4.41641957979107 * x39
    x41 = x11 * x28
    x42 = x32 * x41
    x43 = x42 * (x14 + x34)
    x44 = x11 * x29
    x45 = x43 * x44
    x46 = x26 * x31
    x47 = x25 * x31
    x48 = 11.6847478934435 * x39
    x49 = x36 * x44
    x50 = x49 * (x14 + x38)
    x51 = x32 * x43
    x52 = x33 * x41
    x53 = x1 * x41
    x54 = x52 + x53
    x55 = 2.0 * x41
    x56 = x8 * (x0 * x54 - x55)
    x57 = x34 * x53
    x58 = x51 + x56 + x57
    x59 = x19 * x44
    x60 = 15.084944665313 * x39
    x61 = 26.1278905896872 * x39
    x62 = x36 * x50
    x63 = x37 * x44
    x64 = x1 * x44
    x65 = x63 + x64
    x66 = 2.0 * x44
    x67 = x8 * (x0 * x65 - x66)
    x68 = x38 * x64
    x69 = x62 + x67 + x68
    x70 = x19 * x41
    x71 = x32 * (2.0 * x53 + x54)
    x72 = x1 * (x23 * x42 + x32 * x34 * x55) + x32 * x58 + x8 * (x0 * x71 - 3.0 * x42)
    x73 = x10 * x30
    x74 = x72 * x73
    x75 = x3 * x73
    x76 = 3.14159265358979 * x10 * x2 * x28
    x77 = x3 * x76
    x78 = x36 * (2.0 * x64 + x65)
    x79 = x1 * (x23 * x49 + x36 * x38 * x66) + x36 * x69 + x8 * (x0 * x78 - 3.0 * x49)
    x80 = x76 * x79
    x81 = 3.0 * x1 * (x52 + x53) + x32 * x71
    x82 = x73 * x81
    x83 = x12 * x54
    x84 = 3.0 * x1 * (x63 + x64) + x36 * x78
    x85 = x76 * x84

    # 15 item(s)
    result[0, 0] = numpy.sum(
        x40
        * (
            x31
            * (
                3.0 * x1 * (x16 + x21 + x22)
                + x26 * x3
                - 2.0 * x8 * (-ax * x27 + 2.0 * x17 + 2.0 * x18)
            )
            + x34 * x35
            + x35 * x38
        )
    )
    result[1, 0] = numpy.sum(x48 * (x25 * x45 + x32 * x38 * x47 + x32 * x46))
    result[2, 0] = numpy.sum(x48 * (x25 * x41 * x50 + x34 * x36 * x47 + x36 * x46))
    result[3, 0] = numpy.sum(x60 * (x24 * x44 * x54 + x38 * x54 * x59 + x58 * x59))
    result[4, 0] = numpy.sum(
        x61 * (x19 * x36 * x45 + x19 * x42 * x50 + x24 * x31 * x32 * x36)
    )
    result[5, 0] = numpy.sum(x60 * (x24 * x41 * x65 + x34 * x65 * x70 + x69 * x70))
    result[6, 0] = numpy.sum(x48 * (x15 * x44 * x71 + x3 * x74 + x38 * x71 * x75))
    result[7, 0] = numpy.sum(x61 * (x13 * x50 * x54 + x15 * x49 * x54 + x36 * x58 * x75))
    result[8, 0] = numpy.sum(x61 * (x13 * x43 * x65 + x15 * x42 * x65 + x32 * x69 * x77))
    result[9, 0] = numpy.sum(x48 * (x15 * x41 * x78 + x3 * x80 + x34 * x77 * x78))
    result[10, 0] = numpy.sum(
        x40
        * (
            x38 * x82
            + x7 * x82
            + x73
            * (
                3.0 * x1 * (x51 + x56 + x57)
                + x32 * x72
                - 2.0 * x8 * (-ax * x81 + 2.0 * x52 + 2.0 * x53)
            )
        )
    )
    result[11, 0] = numpy.sum(x48 * (x12 * x50 * x71 + x36 * x7 * x71 * x73 + x36 * x74))
    result[12, 0] = numpy.sum(x60 * (x12 * x58 * x65 + x65 * x7 * x83 + x69 * x83))
    result[13, 0] = numpy.sum(x48 * (x12 * x43 * x78 + x32 * x7 * x76 * x78 + x32 * x80))
    result[14, 0] = numpy.sum(
        x40
        * (
            x34 * x85
            + x7 * x85
            + x76
            * (
                3.0 * x1 * (x62 + x67 + x68)
                + x36 * x79
                - 2.0 * x8 * (-ax * x84 + 2.0 * x63 + 2.0 * x64)
            )
        )
    )
    return result


def kinetic3d_41(ax, da, A, bx, db, B):
    """Cartesian 3D (gp) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 3), dtype=float)

    x0 = 2.0 * ax
    x1 = (2.0 * bx + x0) ** (-1.0)
    x2 = -ax
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - A[0]
    x6 = x5**2
    x7 = 2.0 * ax**2
    x8 = -x2 - x7 * (x1 + x6)
    x9 = bx * x3
    x10 = ax * x9
    x11 = numpy.exp(-x10 * (A[0] - B[0]) ** 2)
    x12 = 1.77245385090552 * numpy.sqrt(x3)
    x13 = x11 * x12
    x14 = x13 * x5
    x15 = x0 * x9
    x16 = x14 * (x15 + x8)
    x17 = -x4 - B[0]
    x18 = x13 * x17
    x19 = x18 * (x15 + x8)
    x20 = x1 * (x16 + x19)
    x21 = x19 * x5
    x22 = x1 * x13
    x23 = x14 * x17
    x24 = x22 + x23
    x25 = x22 * x8
    x26 = x15 * x24 + x21 + x25
    x27 = x26 * x5
    x28 = x1 * (x14 + x18)
    x29 = x24 * x5
    x30 = x28 + x29
    x31 = 2.0 * x13
    x32 = x17 * x31
    x33 = x9 * (x0 * x30 - x32)
    x34 = 4.0 * x10
    x35 = x16 * x5
    x36 = x13 * x6
    x37 = x22 + x36
    x38 = x9 * (x0 * x37 - x31)
    x39 = x35 + x38
    x40 = x25 + x39
    x41 = x5 * (2.0 * x22 + x37)
    x42 = x1 * (x14 * x34 + x31 * x5 * x8) + x40 * x5 + x9 * (x0 * x41 - 3.0 * x14)
    x43 = 3.0 * x25
    x44 = x20 + x27 + x33
    x45 = 3.0 * x22
    x46 = x1 * (x32 * x5 + x36 + x45) + x30 * x5
    x47 = (
        x1 * (2.0 * x21 + x24 * x34 + x39 + x43)
        + x44 * x5
        - x9 * (-2.0 * ax * x46 + 3.0 * x23 + x45)
    )
    x48 = x1 * (3.0 * x28 + 3.0 * x29 + x41) + x46 * x5
    x49 = numpy.exp(-x10 * (A[1] - B[1]) ** 2)
    x50 = numpy.exp(-x10 * (A[2] - B[2]) ** 2)
    x51 = 3.14159265358979 * x3 * x50
    x52 = x49 * x51
    x53 = -x3 * (ax * A[1] + bx * B[1])
    x54 = -x53 - A[1]
    x55 = x54**2
    x56 = -x2 - x7 * (x1 + x55)
    x57 = x48 * x52
    x58 = -x3 * (ax * A[2] + bx * B[2])
    x59 = -x58 - A[2]
    x60 = x59**2
    x61 = -x2 - x7 * (x1 + x60)
    x62 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**2.5)
    x63 = 8.83283915958214 * x62
    x64 = -x53 - B[1]
    x65 = x12 * x49
    x66 = x64 * x65
    x67 = x66 * (x15 + x56)
    x68 = x1 * (3.0 * x36 + x45) + x41 * x5
    x69 = x12 * x50
    x70 = x52 * (
        x1 * (3.0 * x35 + 3.0 * x38 + x43)
        + x42 * x5
        - 2.0 * x9 * (-ax * x68 + 2.0 * x22 + 2.0 * x36)
    )
    x71 = x52 * x68
    x72 = -x58 - B[2]
    x73 = x69 * x72
    x74 = x73 * (x15 + x61)
    x75 = x54 * x65
    x76 = x75 * (x15 + x56)
    x77 = x47 * x52
    x78 = x52 * x54
    x79 = 23.3694957868871 * x62
    x80 = x54 * x67
    x81 = x1 * x12
    x82 = x49 * x81
    x83 = x64 * x75
    x84 = x82 + x83
    x85 = x56 * x82
    x86 = x15 * x84 + x80 + x85
    x87 = x41 * x69
    x88 = x59 * x69
    x89 = x88 * (x15 + x61)
    x90 = x52 * x59
    x91 = x59 * x74
    x92 = x50 * x81
    x93 = x72 * x88
    x94 = x92 + x93
    x95 = x61 * x92
    x96 = x15 * x94 + x91 + x95
    x97 = x41 * x65
    x98 = x54 * x76
    x99 = x55 * x65
    x100 = x82 + x99
    x101 = 2.0 * x65
    x102 = x9 * (x0 * x100 - x101)
    x103 = x102 + x98
    x104 = x103 + x85
    x105 = x30 * x69
    x106 = 30.169889330626 * x62
    x107 = x1 * (x66 + x75)
    x108 = x54 * x84
    x109 = x107 + x108
    x110 = x109 * x69
    x111 = x1 * (x67 + x76)
    x112 = x54 * x86
    x113 = x101 * x64
    x114 = x9 * (x0 * x109 - x113)
    x115 = x111 + x112 + x114
    x116 = 52.2557811793745 * x62
    x117 = x59 * x89
    x118 = x60 * x69
    x119 = x118 + x92
    x120 = 2.0 * x69
    x121 = x9 * (x0 * x119 - x120)
    x122 = x117 + x121
    x123 = x122 + x95
    x124 = x30 * x65
    x125 = x1 * (x73 + x88)
    x126 = x59 * x94
    x127 = x125 + x126
    x128 = x127 * x65
    x129 = x1 * (x74 + x89)
    x130 = x59 * x96
    x131 = x120 * x72
    x132 = x9 * (x0 * x127 - x131)
    x133 = x129 + x130 + x132
    x134 = x54 * (x100 + 2.0 * x82)
    x135 = x134 * x69
    x136 = x1 * (x101 * x54 * x56 + x34 * x75) + x104 * x54 + x9 * (x0 * x134 - 3.0 * x75)
    x137 = 3.0 * x82
    x138 = x1 * (x113 * x54 + x137 + x99) + x109 * x54
    x139 = 3.0 * x85
    x140 = (
        x1 * (x103 + x139 + x34 * x84 + 2.0 * x80)
        + x115 * x54
        - x9 * (-2.0 * ax * x138 + x137 + 3.0 * x83)
    )
    x141 = x11 * x51
    x142 = x140 * x141
    x143 = x141 * x5
    x144 = 3.14159265358979 * x11 * x3 * x49
    x145 = x144 * x5
    x146 = x59 * (x119 + 2.0 * x92)
    x147 = x146 * x65
    x148 = x1 * (x120 * x59 * x61 + x34 * x88) + x123 * x59 + x9 * (x0 * x146 - 3.0 * x88)
    x149 = 3.0 * x92
    x150 = x1 * (x118 + x131 * x59 + x149) + x127 * x59
    x151 = 3.0 * x95
    x152 = (
        x1 * (x122 + x151 + x34 * x94 + 2.0 * x91)
        + x133 * x59
        - x9 * (-2.0 * ax * x150 + x149 + 3.0 * x93)
    )
    x153 = x144 * x152
    x154 = x1 * (x137 + 3.0 * x99) + x134 * x54
    x155 = x141 * (
        x1 * (3.0 * x102 + x139 + 3.0 * x98)
        + x136 * x54
        - 2.0 * x9 * (-ax * x154 + 2.0 * x82 + 2.0 * x99)
    )
    x156 = x141 * x154
    x157 = x1 * (3.0 * x107 + 3.0 * x108 + x134) + x138 * x54
    x158 = x141 * x157
    x159 = x141 * x59
    x160 = x13 * x134
    x161 = x109 * x13
    x162 = x127 * x13
    x163 = x144 * x54
    x164 = x13 * x146
    x165 = x1 * (3.0 * x118 + x149) + x146 * x59
    x166 = x144 * (
        x1 * (3.0 * x117 + 3.0 * x121 + x151)
        + x148 * x59
        - 2.0 * x9 * (-ax * x165 + 2.0 * x118 + 2.0 * x92)
    )
    x167 = x144 * x165
    x168 = x1 * (3.0 * x125 + 3.0 * x126 + x146) + x150 * x59
    x169 = x144 * x168

    # 45 item(s)
    result[0, 0] = numpy.sum(
        x63
        * (
            x52
            * (
                x1 * (3.0 * x20 + 3.0 * x27 + 3.0 * x33 + x42)
                + x47 * x5
                - 2.0 * x9 * (-ax * x48 + 2.0 * x28 + 2.0 * x29)
            )
            + x56 * x57
            + x57 * x61
        )
    )
    result[0, 1] = numpy.sum(x63 * (x61 * x64 * x71 + x64 * x70 + x67 * x68 * x69))
    result[0, 2] = numpy.sum(x63 * (x56 * x71 * x72 + x65 * x68 * x74 + x70 * x72))
    result[1, 0] = numpy.sum(x79 * (x46 * x61 * x78 + x46 * x69 * x76 + x54 * x77))
    result[1, 1] = numpy.sum(x79 * (x42 * x69 * x84 + x61 * x84 * x87 + x86 * x87))
    result[1, 2] = numpy.sum(x79 * (x41 * x73 * x76 + x41 * x74 * x75 + x42 * x72 * x78))
    result[2, 0] = numpy.sum(x79 * (x46 * x56 * x90 + x46 * x65 * x89 + x59 * x77))
    result[2, 1] = numpy.sum(x79 * (x41 * x66 * x89 + x41 * x67 * x88 + x42 * x64 * x90))
    result[2, 2] = numpy.sum(x79 * (x42 * x65 * x94 + x56 * x94 * x97 + x96 * x97))
    result[3, 0] = numpy.sum(x106 * (x100 * x105 * x61 + x100 * x44 * x69 + x104 * x105))
    result[3, 1] = numpy.sum(x106 * (x110 * x37 * x61 + x110 * x40 + x115 * x37 * x69))
    result[3, 2] = numpy.sum(
        x106 * (x100 * x37 * x74 + x100 * x40 * x73 + x104 * x37 * x73)
    )
    result[4, 0] = numpy.sum(x116 * (x30 * x75 * x89 + x30 * x76 * x88 + x44 * x59 * x78))
    result[4, 1] = numpy.sum(x116 * (x37 * x84 * x89 + x37 * x86 * x88 + x40 * x84 * x88))
    result[4, 2] = numpy.sum(x116 * (x37 * x75 * x96 + x37 * x76 * x94 + x40 * x75 * x94))
    result[5, 0] = numpy.sum(x106 * (x119 * x124 * x56 + x119 * x44 * x65 + x123 * x124))
    result[5, 1] = numpy.sum(
        x106 * (x119 * x37 * x67 + x119 * x40 * x66 + x123 * x37 * x66)
    )
    result[5, 2] = numpy.sum(x106 * (x128 * x37 * x56 + x128 * x40 + x133 * x37 * x65))
    result[6, 0] = numpy.sum(x79 * (x135 * x24 * x61 + x135 * x26 + x136 * x24 * x69))
    result[6, 1] = numpy.sum(x79 * (x138 * x143 * x61 + x138 * x16 * x69 + x142 * x5))
    result[6, 2] = numpy.sum(
        x79 * (x134 * x14 * x74 + x134 * x16 * x73 + x136 * x143 * x72)
    )
    result[7, 0] = numpy.sum(
        x116 * (x100 * x24 * x89 + x100 * x26 * x88 + x104 * x24 * x88)
    )
    result[7, 1] = numpy.sum(
        x116 * (x109 * x14 * x89 + x109 * x16 * x88 + x115 * x143 * x59)
    )
    result[7, 2] = numpy.sum(
        x116 * (x100 * x14 * x96 + x100 * x16 * x94 + x104 * x14 * x94)
    )
    result[8, 0] = numpy.sum(
        x116 * (x119 * x24 * x76 + x119 * x26 * x75 + x123 * x24 * x75)
    )
    result[8, 1] = numpy.sum(
        x116 * (x119 * x14 * x86 + x119 * x16 * x84 + x123 * x14 * x84)
    )
    result[8, 2] = numpy.sum(
        x116 * (x127 * x14 * x76 + x127 * x16 * x75 + x133 * x145 * x54)
    )
    result[9, 0] = numpy.sum(x79 * (x147 * x24 * x56 + x147 * x26 + x148 * x24 * x65))
    result[9, 1] = numpy.sum(
        x79 * (x14 * x146 * x67 + x145 * x148 * x64 + x146 * x16 * x66)
    )
    result[9, 2] = numpy.sum(x79 * (x145 * x150 * x56 + x150 * x16 * x65 + x153 * x5))
    result[10, 0] = numpy.sum(x63 * (x154 * x19 * x69 + x155 * x17 + x156 * x17 * x61))
    result[10, 1] = numpy.sum(
        x63
        * (
            x141
            * (
                x1 * (3.0 * x111 + 3.0 * x112 + 3.0 * x114 + x136)
                + x140 * x54
                - 2.0 * x9 * (-ax * x157 + 2.0 * x107 + 2.0 * x108)
            )
            + x158 * x61
            + x158 * x8
        )
    )
    result[10, 2] = numpy.sum(x63 * (x13 * x154 * x74 + x155 * x72 + x156 * x72 * x8))
    result[11, 0] = numpy.sum(
        x79 * (x134 * x18 * x89 + x134 * x19 * x88 + x136 * x159 * x17)
    )
    result[11, 1] = numpy.sum(x79 * (x13 * x138 * x89 + x138 * x159 * x8 + x142 * x59))
    result[11, 2] = numpy.sum(x79 * (x13 * x136 * x94 + x160 * x8 * x94 + x160 * x96))
    result[12, 0] = numpy.sum(
        x106 * (x100 * x119 * x19 + x100 * x123 * x18 + x104 * x119 * x18)
    )
    result[12, 1] = numpy.sum(x106 * (x115 * x119 * x13 + x119 * x161 * x8 + x123 * x161))
    result[12, 2] = numpy.sum(x106 * (x100 * x13 * x133 + x100 * x162 * x8 + x104 * x162))
    result[13, 0] = numpy.sum(
        x79 * (x146 * x18 * x76 + x146 * x19 * x75 + x148 * x163 * x17)
    )
    result[13, 1] = numpy.sum(x79 * (x13 * x148 * x84 + x164 * x8 * x84 + x164 * x86))
    result[13, 2] = numpy.sum(x79 * (x13 * x150 * x76 + x150 * x163 * x8 + x153 * x54))
    result[14, 0] = numpy.sum(x63 * (x165 * x19 * x65 + x166 * x17 + x167 * x17 * x56))
    result[14, 1] = numpy.sum(x63 * (x13 * x165 * x67 + x166 * x64 + x167 * x64 * x8))
    result[14, 2] = numpy.sum(
        x63
        * (
            x144
            * (
                x1 * (3.0 * x129 + 3.0 * x130 + 3.0 * x132 + x148)
                + x152 * x59
                - 2.0 * x9 * (-ax * x168 + 2.0 * x125 + 2.0 * x126)
            )
            + x169 * x56
            + x169 * x8
        )
    )
    return result


def kinetic3d_42(ax, da, A, bx, db, B):
    """Cartesian 3D (gd) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 6), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - A[0]
    x6 = -ax
    x7 = x5**2
    x8 = 2.0 * ax**2
    x9 = -x6 - x8 * (x2 + x7)
    x10 = bx * x3
    x11 = ax * x10
    x12 = numpy.exp(-x11 * (A[0] - B[0]) ** 2)
    x13 = 1.77245385090552 * numpy.sqrt(x3)
    x14 = x12 * x13
    x15 = x14 * x5
    x16 = x0 * x10
    x17 = x15 * (x16 + x9)
    x18 = x17 * x5
    x19 = x14 * x7
    x20 = x14 * x2
    x21 = x19 + x20
    x22 = 2.0 * x14
    x23 = -x22
    x24 = x10 * (x0 * x21 + x23)
    x25 = x18 + x24
    x26 = x20 * x9
    x27 = 3.0 * x26
    x28 = -x4 - B[0]
    x29 = x14 * x28
    x30 = x29 * (x16 + x9)
    x31 = x30 * x5
    x32 = x15 * x28
    x33 = x20 + x32
    x34 = 4.0 * x11
    x35 = x27 + 2.0 * x31 + x33 * x34
    x36 = x2 * (x25 + x35)
    x37 = x2 * (x17 + x30)
    x38 = x16 * x33 + x26 + x31
    x39 = x38 * x5
    x40 = x2 * (x15 + x29)
    x41 = x33 * x5
    x42 = x40 + x41
    x43 = x22 * x28
    x44 = x0 * x42 - x43
    x45 = x10 * x44
    x46 = x37 + x39 + x45
    x47 = x46 * x5
    x48 = x28**2
    x49 = x14 * x48
    x50 = x20 + x49
    x51 = ax * x3
    x52 = x28 * x30 + x51 * (x1 * x50 + x23)
    x53 = x2 * (x35 + x52)
    x54 = x2 * (x29 * x34 + x43 * x9)
    x55 = x26 + x52
    x56 = x5 * x55
    x57 = x5 * x50
    x58 = 2.0 * x20
    x59 = x28 * x58 + x57
    x60 = x16 * x59 + x54 + x56
    x61 = x5 * x60
    x62 = 3.0 * x20
    x63 = x43 * x5 + x62
    x64 = x2 * (x19 + x63)
    x65 = x42 * x5
    x66 = x64 + x65
    x67 = 2.0 * ax * x66 - 3.0 * x32 - x62
    x68 = x1 * x3
    x69 = x2 * (x49 + x63)
    x70 = x5 * x59
    x71 = x69 + x70
    x72 = x10 * (2.0 * ax * x71 - x22 * x48 - x58)
    x73 = x53 + x61 + x72
    x74 = 4.0 * x20
    x75 = x2 * (x28 * x74 + 2.0 * x40 + 2.0 * x41 + 2.0 * x57) + x5 * x71
    x76 = (
        -x10 * (-2.0 * ax * x75 + 6.0 * x20 * x28 + 3.0 * x57)
        + x2 * (x34 * x59 + 2.0 * x37 + 2.0 * x39 + x44 * x68 + 2.0 * x54 + 2.0 * x56)
        + x5 * x73
    )
    x77 = x2 * (2.0 * x64 + 2.0 * x65 + 3.0 * x69 + 3.0 * x70) + x5 * x75
    x78 = numpy.exp(-x11 * (A[1] - B[1]) ** 2)
    x79 = numpy.exp(-x11 * (A[2] - B[2]) ** 2)
    x80 = 3.14159265358979 * x3 * x79
    x81 = x78 * x80
    x82 = -x3 * (ax * A[1] + bx * B[1])
    x83 = -x82 - A[1]
    x84 = x83**2
    x85 = -x6 - x8 * (x2 + x84)
    x86 = x77 * x81
    x87 = -x3 * (ax * A[2] + bx * B[2])
    x88 = -x87 - A[2]
    x89 = x88**2
    x90 = -x6 - x8 * (x2 + x89)
    x91 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**3.5)
    x92 = 10.1992841329868 * x91
    x93 = -x82 - B[1]
    x94 = x13 * x78
    x95 = x93 * x94
    x96 = x95 * (x16 + x85)
    x97 = x5 * (x21 + x58)
    x98 = x2 * (3.0 * x40 + 3.0 * x41 + x97) + x5 * x66
    x99 = x13 * x79
    x100 = x25 + x26
    x101 = x10 * (x0 * x97 - 3.0 * x15) + x100 * x5 + x2 * (x15 * x34 + x22 * x5 * x9)
    x102 = x10 * x67 + x36 + x47
    x103 = x81 * (
        -2.0 * x10 * (-ax * x98 + 2.0 * x40 + 2.0 * x41)
        + x102 * x5
        + x2 * (x101 + 3.0 * x37 + 3.0 * x39 + 3.0 * x45)
    )
    x104 = x81 * x98
    x105 = 17.6656783191643 * x91
    x106 = -x87 - B[2]
    x107 = x106 * x99
    x108 = x107 * (x16 + x90)
    x109 = x2 * x94
    x110 = x109 * x85
    x111 = x93**2 * x94
    x112 = x109 + x111
    x113 = 2.0 * x94
    x114 = -x113
    x115 = x51 * (x1 * x112 + x114) + x93 * x96
    x116 = x110 + x115
    x117 = x2 * (3.0 * x19 + x62) + x5 * x97
    x118 = x117 * x99
    x119 = (
        -x10 * (-2.0 * ax * x117 + 4.0 * x19 + x74)
        + x101 * x5
        + x2 * (3.0 * x18 + 3.0 * x24 + x27)
    )
    x120 = x106 * x81
    x121 = x2 * x99
    x122 = x121 * x90
    x123 = x106**2 * x99
    x124 = x121 + x123
    x125 = 2.0 * x99
    x126 = -x125
    x127 = x106 * x108 + x51 * (x1 * x124 + x126)
    x128 = x122 + x127
    x129 = x117 * x94
    x130 = x83 * x94
    x131 = x130 * (x16 + x85)
    x132 = x76 * x81
    x133 = x75 * x81
    x134 = 26.9847693667702 * x91
    x135 = x83 * x96
    x136 = x130 * x93
    x137 = x109 + x136
    x138 = x110 + x135 + x137 * x16
    x139 = x66 * x99
    x140 = 46.7389915737742 * x91
    x141 = x113 * x93
    x142 = x2 * (x141 * x85 + x34 * x95)
    x143 = x116 * x83
    x144 = x112 * x83
    x145 = 2.0 * x109
    x146 = x144 + x145 * x93
    x147 = x142 + x143 + x146 * x16
    x148 = x97 * x99
    x149 = x88 * x99
    x150 = x149 * (x16 + x90)
    x151 = x81 * x88
    x152 = x108 * x88
    x153 = x106 * x149
    x154 = x121 + x153
    x155 = x122 + x152 + x154 * x16
    x156 = x66 * x94
    x157 = x106 * x125
    x158 = x2 * (x107 * x34 + x157 * x90)
    x159 = x128 * x88
    x160 = x124 * x88
    x161 = 2.0 * x121
    x162 = x106 * x161 + x160
    x163 = x158 + x159 + x16 * x162
    x164 = x94 * x97
    x165 = x131 * x83
    x166 = x84 * x94
    x167 = x109 + x166
    x168 = x10 * (x0 * x167 + x114)
    x169 = x165 + x168
    x170 = x110 + x169
    x171 = x71 * x99
    x172 = 34.8371874529163 * x91
    x173 = x2 * (x131 + x96)
    x174 = x138 * x83
    x175 = x2 * (x130 + x95)
    x176 = x137 * x83
    x177 = x175 + x176
    x178 = x0 * x177 - x141
    x179 = x10 * x178
    x180 = x173 + x174 + x179
    x181 = x42 * x99
    x182 = 60.3397786612521 * x91
    x183 = 3.0 * x109
    x184 = x141 * x83 + x183
    x185 = x2 * (x111 + x184)
    x186 = x146 * x83
    x187 = x185 + x186
    x188 = x187 * x99
    x189 = 3.0 * x110
    x190 = 2.0 * x135 + x137 * x34 + x189
    x191 = x2 * (x115 + x190)
    x192 = x147 * x83
    x193 = x10 * (2.0 * ax * x187 - 2.0 * x111 - x145)
    x194 = x191 + x192 + x193
    x195 = 60.3397786612521 * x91
    x196 = 104.511562358749 * x91
    x197 = x150 * x88
    x198 = x89 * x99
    x199 = x121 + x198
    x200 = x10 * (x0 * x199 + x126)
    x201 = x197 + x200
    x202 = x122 + x201
    x203 = x71 * x94
    x204 = x2 * (x108 + x150)
    x205 = x155 * x88
    x206 = x2 * (x107 + x149)
    x207 = x154 * x88
    x208 = x206 + x207
    x209 = x0 * x208 - x157
    x210 = x10 * x209
    x211 = x204 + x205 + x210
    x212 = x42 * x94
    x213 = 3.0 * x121
    x214 = x157 * x88 + x213
    x215 = x2 * (x123 + x214)
    x216 = x162 * x88
    x217 = x215 + x216
    x218 = x217 * x94
    x219 = 3.0 * x122
    x220 = 2.0 * x152 + x154 * x34 + x219
    x221 = x2 * (x127 + x220)
    x222 = x163 * x88
    x223 = x10 * (2.0 * ax * x217 - 2.0 * x123 - x161)
    x224 = x221 + x222 + x223
    x225 = x83 * (x145 + x167)
    x226 = (
        x10 * (x0 * x225 - 3.0 * x130) + x170 * x83 + x2 * (x113 * x83 * x85 + x130 * x34)
    )
    x227 = x59 * x99
    x228 = x2 * (x166 + x184)
    x229 = x177 * x83
    x230 = x228 + x229
    x231 = x230 * x99
    x232 = x2 * (x169 + x190)
    x233 = x180 * x83
    x234 = 2.0 * ax * x230 - 3.0 * x136 - x183
    x235 = x10 * x234 + x232 + x233
    x236 = 4.0 * x109
    x237 = x187 * x83 + x2 * (2.0 * x144 + 2.0 * x175 + 2.0 * x176 + x236 * x93)
    x238 = (
        -x10 * (-2.0 * ax * x237 + 6.0 * x109 * x93 + 3.0 * x144)
        + x194 * x83
        + x2
        * (2.0 * x142 + 2.0 * x143 + x146 * x34 + 2.0 * x173 + 2.0 * x174 + x178 * x68)
    )
    x239 = x12 * x80
    x240 = x238 * x239
    x241 = x239 * x5
    x242 = 3.14159265358979 * x12 * x3 * x78
    x243 = x242 * x5
    x244 = x88 * (x161 + x199)
    x245 = (
        x10 * (x0 * x244 - 3.0 * x149) + x2 * (x125 * x88 * x90 + x149 * x34) + x202 * x88
    )
    x246 = x59 * x94
    x247 = x2 * (x198 + x214)
    x248 = x208 * x88
    x249 = x247 + x248
    x250 = x249 * x94
    x251 = x2 * (x201 + x220)
    x252 = x211 * x88
    x253 = 2.0 * ax * x249 - 3.0 * x153 - x213
    x254 = x10 * x253 + x251 + x252
    x255 = 4.0 * x121
    x256 = x2 * (x106 * x255 + 2.0 * x160 + 2.0 * x206 + 2.0 * x207) + x217 * x88
    x257 = (
        -x10 * (-2.0 * ax * x256 + 6.0 * x106 * x121 + 3.0 * x160)
        + x2
        * (2.0 * x158 + 2.0 * x159 + x162 * x34 + 2.0 * x204 + 2.0 * x205 + x209 * x68)
        + x224 * x88
    )
    x258 = x242 * x257
    x259 = x2 * (3.0 * x166 + x183) + x225 * x83
    x260 = x259 * x99
    x261 = (
        -x10 * (-2.0 * ax * x259 + 4.0 * x166 + x236)
        + x2 * (3.0 * x165 + 3.0 * x168 + x189)
        + x226 * x83
    )
    x262 = x2 * (3.0 * x175 + 3.0 * x176 + x225) + x230 * x83
    x263 = x239 * (
        -2.0 * x10 * (-ax * x262 + 2.0 * x175 + 2.0 * x176)
        + x2 * (3.0 * x173 + 3.0 * x174 + 3.0 * x179 + x226)
        + x235 * x83
    )
    x264 = x239 * x28
    x265 = x2 * (3.0 * x185 + 3.0 * x186 + 2.0 * x228 + 2.0 * x229) + x237 * x83
    x266 = x239 * x265
    x267 = x239 * x9
    x268 = x14 * x259
    x269 = x14 * x230
    x270 = x14 * x225
    x271 = x14 * x187
    x272 = x14 * x177
    x273 = x14 * x217
    x274 = x242 * x83
    x275 = x14 * x146
    x276 = x14 * x249
    x277 = x2 * (3.0 * x198 + x213) + x244 * x88
    x278 = x277 * x94
    x279 = (
        -x10 * (-2.0 * ax * x277 + 4.0 * x198 + x255)
        + x2 * (3.0 * x197 + 3.0 * x200 + x219)
        + x245 * x88
    )
    x280 = x242 * x28
    x281 = x2 * (3.0 * x206 + 3.0 * x207 + x244) + x249 * x88
    x282 = x242 * (
        -2.0 * x10 * (-ax * x281 + 2.0 * x206 + 2.0 * x207)
        + x2 * (3.0 * x204 + 3.0 * x205 + 3.0 * x210 + x245)
        + x254 * x88
    )
    x283 = x14 * x277
    x284 = x2 * (3.0 * x215 + 3.0 * x216 + 2.0 * x247 + 2.0 * x248) + x256 * x88
    x285 = x242 * x284

    # 90 item(s)
    result[0, 0] = numpy.sum(
        x92
        * (
            x81
            * (
                -2.0 * x10 * (-ax * x77 + 2.0 * x69 + 2.0 * x70)
                + x2
                * (2.0 * x36 + 2.0 * x47 + 3.0 * x53 + 3.0 * x61 + x67 * x68 + 3.0 * x72)
                + x5 * x76
            )
            + x85 * x86
            + x86 * x90
        )
    )
    result[0, 1] = numpy.sum(x105 * (x103 * x93 + x104 * x90 * x93 + x96 * x98 * x99))
    result[0, 2] = numpy.sum(x105 * (x103 * x106 + x104 * x106 * x85 + x108 * x94 * x98))
    result[0, 3] = numpy.sum(x92 * (x112 * x118 * x90 + x112 * x119 * x99 + x116 * x118))
    result[0, 4] = numpy.sum(
        x105 * (x107 * x117 * x96 + x108 * x117 * x95 + x119 * x120 * x93)
    )
    result[0, 5] = numpy.sum(x92 * (x119 * x124 * x94 + x124 * x129 * x85 + x128 * x129))
    result[1, 0] = numpy.sum(x134 * (x131 * x75 * x99 + x132 * x83 + x133 * x83 * x90))
    result[1, 1] = numpy.sum(x140 * (x102 * x137 * x99 + x137 * x139 * x90 + x138 * x139))
    result[1, 2] = numpy.sum(
        x140 * (x102 * x120 * x83 + x107 * x131 * x66 + x108 * x130 * x66)
    )
    result[1, 3] = numpy.sum(x134 * (x101 * x146 * x99 + x146 * x148 * x90 + x147 * x148))
    result[1, 4] = numpy.sum(
        x140 * (x101 * x107 * x137 + x107 * x138 * x97 + x108 * x137 * x97)
    )
    result[1, 5] = numpy.sum(
        x134 * (x101 * x124 * x130 + x124 * x131 * x97 + x128 * x130 * x97)
    )
    result[2, 0] = numpy.sum(x134 * (x132 * x88 + x133 * x85 * x88 + x150 * x75 * x94))
    result[2, 1] = numpy.sum(
        x140 * (x102 * x151 * x93 + x149 * x66 * x96 + x150 * x66 * x95)
    )
    result[2, 2] = numpy.sum(x140 * (x102 * x154 * x94 + x154 * x156 * x85 + x155 * x156))
    result[2, 3] = numpy.sum(
        x134 * (x101 * x112 * x149 + x112 * x150 * x97 + x116 * x149 * x97)
    )
    result[2, 4] = numpy.sum(
        x140 * (x101 * x154 * x95 + x154 * x96 * x97 + x155 * x95 * x97)
    )
    result[2, 5] = numpy.sum(x134 * (x101 * x162 * x94 + x162 * x164 * x85 + x163 * x164))
    result[3, 0] = numpy.sum(x172 * (x167 * x171 * x90 + x167 * x73 * x99 + x170 * x171))
    result[3, 1] = numpy.sum(x182 * (x177 * x181 * x90 + x177 * x46 * x99 + x180 * x181))
    result[3, 2] = numpy.sum(
        x182 * (x107 * x167 * x46 + x107 * x170 * x42 + x108 * x167 * x42)
    )
    result[3, 3] = numpy.sum(x172 * (x100 * x188 + x188 * x21 * x90 + x194 * x21 * x99))
    result[3, 4] = numpy.sum(
        x182 * (x100 * x107 * x177 + x107 * x180 * x21 + x108 * x177 * x21)
    )
    result[3, 5] = numpy.sum(
        x172 * (x100 * x124 * x167 + x124 * x170 * x21 + x128 * x167 * x21)
    )
    result[4, 0] = numpy.sum(
        x195 * (x130 * x150 * x71 + x131 * x149 * x71 + x151 * x73 * x83)
    )
    result[4, 1] = numpy.sum(
        x196 * (x137 * x149 * x46 + x137 * x150 * x42 + x138 * x149 * x42)
    )
    result[4, 2] = numpy.sum(
        x196 * (x130 * x154 * x46 + x130 * x155 * x42 + x131 * x154 * x42)
    )
    result[4, 3] = numpy.sum(
        x195 * (x100 * x146 * x149 + x146 * x150 * x21 + x147 * x149 * x21)
    )
    result[4, 4] = numpy.sum(
        x196 * (x100 * x137 * x154 + x137 * x155 * x21 + x138 * x154 * x21)
    )
    result[4, 5] = numpy.sum(
        x195 * (x100 * x130 * x162 + x130 * x163 * x21 + x131 * x162 * x21)
    )
    result[5, 0] = numpy.sum(x172 * (x199 * x203 * x85 + x199 * x73 * x94 + x202 * x203))
    result[5, 1] = numpy.sum(
        x182 * (x199 * x42 * x96 + x199 * x46 * x95 + x202 * x42 * x95)
    )
    result[5, 2] = numpy.sum(x182 * (x208 * x212 * x85 + x208 * x46 * x94 + x211 * x212))
    result[5, 3] = numpy.sum(
        x172 * (x100 * x112 * x199 + x112 * x202 * x21 + x116 * x199 * x21)
    )
    result[5, 4] = numpy.sum(
        x182 * (x100 * x208 * x95 + x208 * x21 * x96 + x21 * x211 * x95)
    )
    result[5, 5] = numpy.sum(x172 * (x100 * x218 + x21 * x218 * x85 + x21 * x224 * x94))
    result[6, 0] = numpy.sum(x134 * (x225 * x227 * x90 + x225 * x60 * x99 + x226 * x227))
    result[6, 1] = numpy.sum(x140 * (x231 * x33 * x90 + x231 * x38 + x235 * x33 * x99))
    result[6, 2] = numpy.sum(
        x140 * (x107 * x225 * x38 + x107 * x226 * x33 + x108 * x225 * x33)
    )
    result[6, 3] = numpy.sum(x134 * (x17 * x237 * x99 + x237 * x241 * x90 + x240 * x5))
    result[6, 4] = numpy.sum(
        x140 * (x106 * x235 * x241 + x107 * x17 * x230 + x108 * x15 * x230)
    )
    result[6, 5] = numpy.sum(
        x134 * (x124 * x15 * x226 + x124 * x17 * x225 + x128 * x15 * x225)
    )
    result[7, 0] = numpy.sum(
        x195 * (x149 * x167 * x60 + x149 * x170 * x59 + x150 * x167 * x59)
    )
    result[7, 1] = numpy.sum(
        x196 * (x149 * x177 * x38 + x149 * x180 * x33 + x150 * x177 * x33)
    )
    result[7, 2] = numpy.sum(
        x196 * (x154 * x167 * x38 + x154 * x170 * x33 + x155 * x167 * x33)
    )
    result[7, 3] = numpy.sum(
        x195 * (x149 * x17 * x187 + x15 * x150 * x187 + x194 * x241 * x88)
    )
    result[7, 4] = numpy.sum(
        x196 * (x15 * x154 * x180 + x15 * x155 * x177 + x154 * x17 * x177)
    )
    result[7, 5] = numpy.sum(
        x195 * (x15 * x162 * x170 + x15 * x163 * x167 + x162 * x167 * x17)
    )
    result[8, 0] = numpy.sum(
        x195 * (x130 * x199 * x60 + x130 * x202 * x59 + x131 * x199 * x59)
    )
    result[8, 1] = numpy.sum(
        x196 * (x137 * x199 * x38 + x137 * x202 * x33 + x138 * x199 * x33)
    )
    result[8, 2] = numpy.sum(
        x196 * (x130 * x208 * x38 + x130 * x211 * x33 + x131 * x208 * x33)
    )
    result[8, 3] = numpy.sum(
        x195 * (x146 * x15 * x202 + x146 * x17 * x199 + x147 * x15 * x199)
    )
    result[8, 4] = numpy.sum(
        x196 * (x137 * x15 * x211 + x137 * x17 * x208 + x138 * x15 * x208)
    )
    result[8, 5] = numpy.sum(
        x195 * (x130 * x17 * x217 + x131 * x15 * x217 + x224 * x243 * x83)
    )
    result[9, 0] = numpy.sum(x134 * (x244 * x246 * x85 + x244 * x60 * x94 + x245 * x246))
    result[9, 1] = numpy.sum(
        x140 * (x244 * x33 * x96 + x244 * x38 * x95 + x245 * x33 * x95)
    )
    result[9, 2] = numpy.sum(x140 * (x250 * x33 * x85 + x250 * x38 + x254 * x33 * x94))
    result[9, 3] = numpy.sum(
        x134 * (x112 * x15 * x245 + x112 * x17 * x244 + x116 * x15 * x244)
    )
    result[9, 4] = numpy.sum(
        x140 * (x15 * x249 * x96 + x17 * x249 * x95 + x243 * x254 * x93)
    )
    result[9, 5] = numpy.sum(x134 * (x17 * x256 * x94 + x243 * x256 * x85 + x258 * x5))
    result[10, 0] = numpy.sum(x92 * (x260 * x50 * x90 + x260 * x55 + x261 * x50 * x99))
    result[10, 1] = numpy.sum(x105 * (x262 * x264 * x90 + x262 * x30 * x99 + x263 * x28))
    result[10, 2] = numpy.sum(
        x105 * (x106 * x261 * x264 + x107 * x259 * x30 + x108 * x259 * x29)
    )
    result[10, 3] = numpy.sum(
        x92
        * (
            x239
            * (
                -2.0 * x10 * (-ax * x265 + 2.0 * x185 + 2.0 * x186)
                + x2
                * (
                    3.0 * x191
                    + 3.0 * x192
                    + 3.0 * x193
                    + 2.0 * x232
                    + 2.0 * x233
                    + x234 * x68
                )
                + x238 * x83
            )
            + x266 * x9
            + x266 * x90
        )
    )
    result[10, 4] = numpy.sum(
        x105 * (x106 * x262 * x267 + x106 * x263 + x108 * x14 * x262)
    )
    result[10, 5] = numpy.sum(x92 * (x124 * x14 * x261 + x124 * x268 * x9 + x128 * x268))
    result[11, 0] = numpy.sum(
        x134 * (x149 * x225 * x55 + x149 * x226 * x50 + x150 * x225 * x50)
    )
    result[11, 1] = numpy.sum(
        x140 * (x149 * x230 * x30 + x150 * x230 * x29 + x235 * x264 * x88)
    )
    result[11, 2] = numpy.sum(
        x140 * (x154 * x225 * x30 + x154 * x226 * x29 + x155 * x225 * x29)
    )
    result[11, 3] = numpy.sum(x134 * (x14 * x150 * x237 + x237 * x267 * x88 + x240 * x88))
    result[11, 4] = numpy.sum(x140 * (x14 * x154 * x235 + x154 * x269 * x9 + x155 * x269))
    result[11, 5] = numpy.sum(x134 * (x14 * x162 * x226 + x162 * x270 * x9 + x163 * x270))
    result[12, 0] = numpy.sum(
        x172 * (x167 * x199 * x55 + x167 * x202 * x50 + x170 * x199 * x50)
    )
    result[12, 1] = numpy.sum(
        x182 * (x177 * x199 * x30 + x177 * x202 * x29 + x180 * x199 * x29)
    )
    result[12, 2] = numpy.sum(
        x182 * (x167 * x208 * x30 + x167 * x211 * x29 + x170 * x208 * x29)
    )
    result[12, 3] = numpy.sum(x172 * (x14 * x194 * x199 + x199 * x271 * x9 + x202 * x271))
    result[12, 4] = numpy.sum(x182 * (x14 * x180 * x208 + x208 * x272 * x9 + x211 * x272))
    result[12, 5] = numpy.sum(x172 * (x14 * x167 * x224 + x167 * x273 * x9 + x170 * x273))
    result[13, 0] = numpy.sum(
        x134 * (x130 * x244 * x55 + x130 * x245 * x50 + x131 * x244 * x50)
    )
    result[13, 1] = numpy.sum(
        x140 * (x137 * x244 * x30 + x137 * x245 * x29 + x138 * x244 * x29)
    )
    result[13, 2] = numpy.sum(
        x140 * (x130 * x249 * x30 + x131 * x249 * x29 + x254 * x274 * x28)
    )
    result[13, 3] = numpy.sum(x134 * (x14 * x147 * x244 + x244 * x275 * x9 + x245 * x275))
    result[13, 4] = numpy.sum(x140 * (x137 * x14 * x254 + x137 * x276 * x9 + x138 * x276))
    result[13, 5] = numpy.sum(x134 * (x131 * x14 * x256 + x256 * x274 * x9 + x258 * x83))
    result[14, 0] = numpy.sum(x92 * (x278 * x50 * x85 + x278 * x55 + x279 * x50 * x94))
    result[14, 1] = numpy.sum(
        x105 * (x277 * x29 * x96 + x277 * x30 * x95 + x279 * x280 * x93)
    )
    result[14, 2] = numpy.sum(x105 * (x28 * x282 + x280 * x281 * x85 + x281 * x30 * x94))
    result[14, 3] = numpy.sum(x92 * (x112 * x14 * x279 + x112 * x283 * x9 + x116 * x283))
    result[14, 4] = numpy.sum(
        x105 * (x14 * x281 * x96 + x242 * x281 * x9 * x93 + x282 * x93)
    )
    result[14, 5] = numpy.sum(
        x92
        * (
            x242
            * (
                -2.0 * x10 * (-ax * x284 + 2.0 * x215 + 2.0 * x216)
                + x2
                * (
                    3.0 * x221
                    + 3.0 * x222
                    + 3.0 * x223
                    + 2.0 * x251
                    + 2.0 * x252
                    + x253 * x68
                )
                + x257 * x88
            )
            + x285 * x85
            + x285 * x9
        )
    )
    return result


def kinetic3d_43(ax, da, A, bx, db, B):
    """Cartesian 3D (gf) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 10), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - A[0]
    x6 = -ax
    x7 = x5**2
    x8 = 2.0 * ax**2
    x9 = -x6 - x8 * (x2 + x7)
    x10 = ax * x3
    x11 = bx * x10
    x12 = numpy.exp(-x11 * (A[0] - B[0]) ** 2)
    x13 = 1.77245385090552 * numpy.sqrt(x3)
    x14 = x12 * x13
    x15 = x14 * x2
    x16 = x15 * x9
    x17 = -x4 - B[0]
    x18 = x14 * x17
    x19 = bx * x3
    x20 = x0 * x19
    x21 = x18 * (x20 + x9)
    x22 = x17 * x21
    x23 = x14 * x17**2
    x24 = x15 + x23
    x25 = 2.0 * x14
    x26 = -x25
    x27 = x10 * (x1 * x24 + x26)
    x28 = x22 + x27
    x29 = x16 + x28
    x30 = x29 * x5
    x31 = x17 * x25
    x32 = 4.0 * x11
    x33 = x2 * (x18 * x32 + x31 * x9)
    x34 = x24 * x5
    x35 = 2.0 * x15
    x36 = x17 * x35
    x37 = x34 + x36
    x38 = x11 * x37
    x39 = x17 * x24
    x40 = x36 + x39
    x41 = x10 * (x1 * x40 - 3.0 * x18) + x17 * x29
    x42 = x2 * (3.0 * x30 + 4.0 * x33 + 6.0 * x38 + x41)
    x43 = x14 * x5
    x44 = x43 * (x20 + x9)
    x45 = x2 * (x21 + x44)
    x46 = x21 * x5
    x47 = x17 * x43
    x48 = x15 + x47
    x49 = x16 + x20 * x48 + x46
    x50 = x49 * x5
    x51 = x2 * (x18 + x43)
    x52 = x48 * x5
    x53 = x51 + x52
    x54 = x0 * x53 - x31
    x55 = x1 * x3
    x56 = x2 * (2.0 * x30 + 2.0 * x33 + 4.0 * x38 + 2.0 * x45 + 2.0 * x50 + x54 * x55)
    x57 = 3.0 * x16
    x58 = x2 * (3.0 * x22 + 3.0 * x27 + x57)
    x59 = x33 + x41
    x60 = x5 * x59
    x61 = 3.0 * x15
    x62 = x2 * (3.0 * x23 + x61)
    x63 = x40 * x5
    x64 = x62 + x63
    x65 = x20 * x64 + x58 + x60
    x66 = x5 * x65
    x67 = x32 * x48 + 2.0 * x46 + x57
    x68 = x2 * (x28 + x67)
    x69 = x20 * x37 + x30 + x33
    x70 = x5 * x69
    x71 = x31 * x5 + x61
    x72 = x2 * (x23 + x71)
    x73 = x37 * x5
    x74 = x72 + x73
    x75 = x19 * (2.0 * ax * x74 - 2.0 * x23 - x35)
    x76 = x68 + x70 + x75
    x77 = x5 * x76
    x78 = 4.0 * x15
    x79 = x17 * x78
    x80 = x2 * (2.0 * x34 + 2.0 * x51 + 2.0 * x52 + x79)
    x81 = x5 * x74
    x82 = x80 + x81
    x83 = 3.0 * x34
    x84 = x15 * x17
    x85 = x19 * (2.0 * ax * x82 - x83 - 6.0 * x84)
    x86 = x2 * (x39 + x83 + 8.0 * x84)
    x87 = x5 * x64
    x88 = x86 + x87
    x89 = x19 * (2.0 * ax * x88 - 2.0 * x39 - x79)
    x90 = 3.0 * x68 + 3.0 * x70 + 3.0 * x75
    x91 = x42 + x66 + x89
    x92 = 3.0 * x72 + 3.0 * x73
    x93 = x2 * (2.0 * x62 + 2.0 * x63 + x92) + x5 * x88
    x94 = (
        -x19 * (-2.0 * ax * x93 + 3.0 * x62 + 3.0 * x63)
        + x2 * (x32 * x64 + 2.0 * x58 + 2.0 * x60 + x90)
        + x5 * x91
    )
    x95 = 3.0 * x2 * (x80 + x81 + x86 + x87) + x5 * x93
    x96 = numpy.exp(-x11 * (A[1] - B[1]) ** 2)
    x97 = numpy.exp(-x11 * (A[2] - B[2]) ** 2)
    x98 = 3.14159265358979 * x3 * x97
    x99 = x96 * x98
    x100 = -x3 * (ax * A[1] + bx * B[1])
    x101 = -x100 - A[1]
    x102 = x101**2
    x103 = -x6 - x8 * (x102 + x2)
    x104 = x95 * x99
    x105 = -x3 * (ax * A[2] + bx * B[2])
    x106 = -x105 - A[2]
    x107 = x106**2
    x108 = -x6 - x8 * (x107 + x2)
    x109 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**4.5)
    x110 = 9.12251705727742 * x109
    x111 = -x100 - B[1]
    x112 = x13 * x96
    x113 = x111 * x112
    x114 = x113 * (x103 + x20)
    x115 = x14 * x7
    x116 = x2 * (x115 + x71)
    x117 = x5 * x53
    x118 = x2 * (2.0 * x116 + 2.0 * x117 + x92) + x5 * x82
    x119 = x13 * x97
    x120 = x44 * x5
    x121 = x115 + x15
    x122 = x19 * (x0 * x121 + x26)
    x123 = x120 + x122
    x124 = x2 * (x123 + x67)
    x125 = x19 * x54
    x126 = x125 + x45 + x50
    x127 = x126 * x5
    x128 = x116 + x117
    x129 = 2.0 * ax * x128 - 3.0 * x47 - x61
    x130 = x56 + x77 + x85
    x131 = x99 * (
        x130 * x5
        - 2.0 * x19 * (-ax * x118 + 2.0 * x72 + 2.0 * x73)
        + x2 * (2.0 * x124 + 2.0 * x127 + x129 * x55 + x90)
    )
    x132 = x118 * x99
    x133 = 20.3985682659737 * x109
    x134 = -x105 - B[2]
    x135 = x119 * x134
    x136 = x135 * (x108 + x20)
    x137 = x112 * x2
    x138 = x103 * x137
    x139 = x111 * x114
    x140 = x111**2 * x112
    x141 = x137 + x140
    x142 = 2.0 * x112
    x143 = -x142
    x144 = x10 * (x1 * x141 + x143)
    x145 = x139 + x144
    x146 = x138 + x145
    x147 = x5 * (x121 + x35)
    x148 = x128 * x5 + x2 * (x147 + 3.0 * x51 + 3.0 * x52)
    x149 = x119 * x148
    x150 = x123 + x16
    x151 = x150 * x5 + x19 * (x0 * x147 - 3.0 * x43) + x2 * (x25 * x5 * x9 + x32 * x43)
    x152 = x124 + x127 + x129 * x19
    x153 = (
        x152 * x5
        - 2.0 * x19 * (-ax * x148 + 2.0 * x51 + 2.0 * x52)
        + x2 * (3.0 * x125 + x151 + 3.0 * x45 + 3.0 * x50)
    )
    x154 = x134 * x99
    x155 = 35.3313566383285 * x109
    x156 = x119 * x2
    x157 = x108 * x156
    x158 = x134 * x136
    x159 = x119 * x134**2
    x160 = x156 + x159
    x161 = 2.0 * x119
    x162 = -x161
    x163 = x10 * (x1 * x160 + x162)
    x164 = x158 + x163
    x165 = x157 + x164
    x166 = x112 * x148
    x167 = x111 * x142
    x168 = x2 * (x103 * x167 + x113 * x32)
    x169 = x111 * x141
    x170 = 2.0 * x137
    x171 = x111 * x170
    x172 = x169 + x171
    x173 = x10 * (x1 * x172 - 3.0 * x113) + x111 * x146
    x174 = x168 + x173
    x175 = x147 * x5 + x2 * (3.0 * x115 + x61)
    x176 = x119 * x175
    x177 = (
        x151 * x5
        - x19 * (-2.0 * ax * x175 + 4.0 * x115 + x78)
        + x2 * (3.0 * x120 + 3.0 * x122 + x57)
    )
    x178 = x134 * x161
    x179 = x2 * (x108 * x178 + x135 * x32)
    x180 = x134 * x160
    x181 = 2.0 * x156
    x182 = x134 * x181
    x183 = x180 + x182
    x184 = x10 * (x1 * x183 - 3.0 * x135) + x134 * x165
    x185 = x179 + x184
    x186 = x112 * x175
    x187 = x101 * x112
    x188 = x187 * (x103 + x20)
    x189 = x94 * x99
    x190 = x93 * x99
    x191 = 24.1359114645008 * x109
    x192 = x101 * x114
    x193 = x111 * x187
    x194 = x137 + x193
    x195 = x138 + x192 + x194 * x20
    x196 = x119 * x82
    x197 = 53.9695387335403 * x109
    x198 = x101 * x146
    x199 = x101 * x141
    x200 = x171 + x199
    x201 = x168 + x198 + x20 * x200
    x202 = x119 * x128
    x203 = 93.4779831475484 * x109
    x204 = 3.0 * x137
    x205 = x2 * (3.0 * x140 + x204)
    x206 = x101 * x172
    x207 = x205 + x206
    x208 = x119 * x207
    x209 = 3.0 * x138
    x210 = x2 * (3.0 * x139 + 3.0 * x144 + x209)
    x211 = x101 * x174
    x212 = x20 * x207 + x210 + x211
    x213 = x106 * x119
    x214 = x213 * (x108 + x20)
    x215 = x106 * x99
    x216 = x106 * x136
    x217 = x134 * x213
    x218 = x156 + x217
    x219 = x157 + x20 * x218 + x216
    x220 = x112 * x82
    x221 = x106 * x165
    x222 = x106 * x160
    x223 = x182 + x222
    x224 = x179 + x20 * x223 + x221
    x225 = x112 * x128
    x226 = 3.0 * x156
    x227 = x2 * (3.0 * x159 + x226)
    x228 = x106 * x183
    x229 = x227 + x228
    x230 = x112 * x229
    x231 = 3.0 * x157
    x232 = x2 * (3.0 * x158 + 3.0 * x163 + x231)
    x233 = x106 * x185
    x234 = x20 * x229 + x232 + x233
    x235 = x101 * x188
    x236 = x102 * x112
    x237 = x137 + x236
    x238 = x19 * (x0 * x237 + x143)
    x239 = x235 + x238
    x240 = x138 + x239
    x241 = x119 * x88
    x242 = 31.1593277158494 * x109
    x243 = x2 * (x113 + x187)
    x244 = x101 * x194
    x245 = x243 + x244
    x246 = x119 * x245
    x247 = x2 * (x114 + x188)
    x248 = x101 * x195
    x249 = x0 * x245 - x167
    x250 = x19 * x249
    x251 = x247 + x248 + x250
    x252 = 69.6743749058326 * x109
    x253 = 2.0 * x192 + x194 * x32 + x209
    x254 = x2 * (x145 + x253)
    x255 = x101 * x201
    x256 = x101 * x167 + x204
    x257 = x2 * (x140 + x256)
    x258 = x101 * x200
    x259 = x257 + x258
    x260 = x19 * (2.0 * ax * x259 - 2.0 * x140 - x170)
    x261 = x254 + x255 + x260
    x262 = x119 * x53
    x263 = 120.679557322504 * x109
    x264 = 3.0 * x199
    x265 = x111 * x137
    x266 = x2 * (x169 + x264 + 8.0 * x265)
    x267 = x101 * x207
    x268 = x266 + x267
    x269 = x119 * x268
    x270 = 6.0 * x11
    x271 = x2 * (4.0 * x168 + x173 + 3.0 * x198 + x200 * x270)
    x272 = x101 * x212
    x273 = 4.0 * x137
    x274 = x111 * x273
    x275 = x19 * (2.0 * ax * x268 - 2.0 * x169 - x274)
    x276 = x271 + x272 + x275
    x277 = 120.679557322504 * x109
    x278 = 209.023124717498 * x109
    x279 = x106 * x214
    x280 = x107 * x119
    x281 = x156 + x280
    x282 = x19 * (x0 * x281 + x162)
    x283 = x279 + x282
    x284 = x157 + x283
    x285 = x112 * x88
    x286 = x2 * (x135 + x213)
    x287 = x106 * x218
    x288 = x286 + x287
    x289 = x112 * x288
    x290 = x2 * (x136 + x214)
    x291 = x106 * x219
    x292 = x0 * x288 - x178
    x293 = x19 * x292
    x294 = x290 + x291 + x293
    x295 = 2.0 * x216 + x218 * x32 + x231
    x296 = x2 * (x164 + x295)
    x297 = x106 * x224
    x298 = x106 * x178 + x226
    x299 = x2 * (x159 + x298)
    x300 = x106 * x223
    x301 = x299 + x300
    x302 = x19 * (2.0 * ax * x301 - 2.0 * x159 - x181)
    x303 = x296 + x297 + x302
    x304 = x112 * x53
    x305 = 3.0 * x222
    x306 = x134 * x156
    x307 = x2 * (x180 + x305 + 8.0 * x306)
    x308 = x106 * x229
    x309 = x307 + x308
    x310 = x112 * x309
    x311 = x2 * (4.0 * x179 + x184 + 3.0 * x221 + x223 * x270)
    x312 = x106 * x234
    x313 = 4.0 * x156
    x314 = x134 * x313
    x315 = x19 * (2.0 * ax * x309 - 2.0 * x180 - x314)
    x316 = x311 + x312 + x315
    x317 = x101 * (x170 + x237)
    x318 = (
        x101 * x240
        + x19 * (x0 * x317 - 3.0 * x187)
        + x2 * (x101 * x103 * x142 + x187 * x32)
    )
    x319 = x119 * x64
    x320 = x2 * (x236 + x256)
    x321 = x101 * x245
    x322 = x320 + x321
    x323 = x119 * x322
    x324 = x2 * (x239 + x253)
    x325 = x101 * x251
    x326 = 2.0 * ax * x322 - 3.0 * x193 - x204
    x327 = x19 * x326 + x324 + x325
    x328 = x2 * (2.0 * x199 + 2.0 * x243 + 2.0 * x244 + x274)
    x329 = x101 * x259
    x330 = x328 + x329
    x331 = x119 * x330
    x332 = x2 * (
        2.0 * x168 + 2.0 * x198 + x200 * x32 + 2.0 * x247 + 2.0 * x248 + x249 * x55
    )
    x333 = x101 * x261
    x334 = x19 * (2.0 * ax * x330 - x264 - 6.0 * x265)
    x335 = x332 + x333 + x334
    x336 = 3.0 * x257 + 3.0 * x258
    x337 = x101 * x268 + x2 * (2.0 * x205 + 2.0 * x206 + x336)
    x338 = 3.0 * x254 + 3.0 * x255 + 3.0 * x260
    x339 = (
        x101 * x276
        - x19 * (-2.0 * ax * x337 + 3.0 * x205 + 3.0 * x206)
        + x2 * (x207 * x32 + 2.0 * x210 + 2.0 * x211 + x338)
    )
    x340 = x12 * x98
    x341 = x339 * x340
    x342 = x340 * x5
    x343 = 3.14159265358979 * x12 * x3 * x96
    x344 = x343 * x5
    x345 = x106 * (x181 + x281)
    x346 = (
        x106 * x284
        + x19 * (x0 * x345 - 3.0 * x213)
        + x2 * (x106 * x108 * x161 + x213 * x32)
    )
    x347 = x112 * x64
    x348 = x2 * (x280 + x298)
    x349 = x106 * x288
    x350 = x348 + x349
    x351 = x112 * x350
    x352 = x2 * (x283 + x295)
    x353 = x106 * x294
    x354 = 2.0 * ax * x350 - 3.0 * x217 - x226
    x355 = x19 * x354 + x352 + x353
    x356 = x2 * (2.0 * x222 + 2.0 * x286 + 2.0 * x287 + x314)
    x357 = x106 * x301
    x358 = x356 + x357
    x359 = x112 * x358
    x360 = x2 * (
        2.0 * x179 + 2.0 * x221 + x223 * x32 + 2.0 * x290 + 2.0 * x291 + x292 * x55
    )
    x361 = x106 * x303
    x362 = x19 * (2.0 * ax * x358 - x305 - 6.0 * x306)
    x363 = x360 + x361 + x362
    x364 = 3.0 * x299 + 3.0 * x300
    x365 = x106 * x309 + x2 * (2.0 * x227 + 2.0 * x228 + x364)
    x366 = 3.0 * x296 + 3.0 * x297 + 3.0 * x302
    x367 = (
        x106 * x316
        - x19 * (-2.0 * ax * x365 + 3.0 * x227 + 3.0 * x228)
        + x2 * (x229 * x32 + 2.0 * x232 + 2.0 * x233 + x366)
    )
    x368 = x343 * x367
    x369 = x101 * x317 + x2 * (x204 + 3.0 * x236)
    x370 = x119 * x369
    x371 = (
        x101 * x318
        - x19 * (-2.0 * ax * x369 + 4.0 * x236 + x273)
        + x2 * (x209 + 3.0 * x235 + 3.0 * x238)
    )
    x372 = x101 * x322 + x2 * (3.0 * x243 + 3.0 * x244 + x317)
    x373 = x119 * x372
    x374 = (
        x101 * x327
        - 2.0 * x19 * (-ax * x372 + 2.0 * x243 + 2.0 * x244)
        + x2 * (3.0 * x247 + 3.0 * x248 + 3.0 * x250 + x318)
    )
    x375 = x101 * x330 + x2 * (2.0 * x320 + 2.0 * x321 + x336)
    x376 = x340 * (
        x101 * x335
        - 2.0 * x19 * (-ax * x375 + 2.0 * x257 + 2.0 * x258)
        + x2 * (2.0 * x324 + 2.0 * x325 + x326 * x55 + x338)
    )
    x377 = x17 * x340
    x378 = x101 * x337 + 3.0 * x2 * (x266 + x267 + x328 + x329)
    x379 = x340 * x378
    x380 = x340 * x9
    x381 = x14 * x372
    x382 = x14 * x369
    x383 = x14 * x330
    x384 = x14 * x322
    x385 = x14 * x229
    x386 = x14 * x268
    x387 = x14 * x288
    x388 = x14 * x245
    x389 = x14 * x309
    x390 = x101 * x343
    x391 = x14 * x207
    x392 = x14 * x350
    x393 = x14 * x358
    x394 = x106 * x345 + x2 * (x226 + 3.0 * x280)
    x395 = x112 * x394
    x396 = (
        x106 * x346
        - x19 * (-2.0 * ax * x394 + 4.0 * x280 + x313)
        + x2 * (x231 + 3.0 * x279 + 3.0 * x282)
    )
    x397 = x106 * x350 + x2 * (3.0 * x286 + 3.0 * x287 + x345)
    x398 = x112 * x397
    x399 = (
        x106 * x355
        - 2.0 * x19 * (-ax * x397 + 2.0 * x286 + 2.0 * x287)
        + x2 * (3.0 * x290 + 3.0 * x291 + 3.0 * x293 + x346)
    )
    x400 = x17 * x343
    x401 = x106 * x358 + x2 * (2.0 * x348 + 2.0 * x349 + x364)
    x402 = x343 * (
        x106 * x363
        - 2.0 * x19 * (-ax * x401 + 2.0 * x299 + 2.0 * x300)
        + x2 * (2.0 * x352 + 2.0 * x353 + x354 * x55 + x366)
    )
    x403 = x14 * x394
    x404 = x14 * x397
    x405 = x106 * x365 + 3.0 * x2 * (x307 + x308 + x356 + x357)
    x406 = x343 * x405

    # 150 item(s)
    result[0, 0] = numpy.sum(
        x110
        * (
            x103 * x104
            + x104 * x108
            + x99
            * (
                -2.0 * x19 * (-ax * x95 + 2.0 * x86 + 2.0 * x87)
                + 3.0 * x2 * (x42 + x56 + x66 + x77 + x85 + x89)
                + x5 * x94
            )
        )
    )
    result[0, 1] = numpy.sum(
        x133 * (x108 * x111 * x132 + x111 * x131 + x114 * x118 * x119)
    )
    result[0, 2] = numpy.sum(
        x133 * (x103 * x132 * x134 + x112 * x118 * x136 + x131 * x134)
    )
    result[0, 3] = numpy.sum(
        x133 * (x108 * x141 * x149 + x119 * x141 * x153 + x146 * x149)
    )
    result[0, 4] = numpy.sum(
        x155 * (x111 * x153 * x154 + x113 * x136 * x148 + x114 * x135 * x148)
    )
    result[0, 5] = numpy.sum(
        x133 * (x103 * x160 * x166 + x112 * x153 * x160 + x165 * x166)
    )
    result[0, 6] = numpy.sum(
        x110 * (x108 * x172 * x176 + x119 * x172 * x177 + x174 * x176)
    )
    result[0, 7] = numpy.sum(
        x133 * (x135 * x141 * x177 + x135 * x146 * x175 + x136 * x141 * x175)
    )
    result[0, 8] = numpy.sum(
        x133 * (x113 * x160 * x177 + x113 * x165 * x175 + x114 * x160 * x175)
    )
    result[0, 9] = numpy.sum(
        x110 * (x103 * x183 * x186 + x112 * x177 * x183 + x185 * x186)
    )
    result[1, 0] = numpy.sum(
        x191 * (x101 * x108 * x190 + x101 * x189 + x119 * x188 * x93)
    )
    result[1, 1] = numpy.sum(
        x197 * (x108 * x194 * x196 + x119 * x130 * x194 + x195 * x196)
    )
    result[1, 2] = numpy.sum(
        x197 * (x101 * x130 * x154 + x135 * x188 * x82 + x136 * x187 * x82)
    )
    result[1, 3] = numpy.sum(
        x197 * (x108 * x200 * x202 + x119 * x152 * x200 + x201 * x202)
    )
    result[1, 4] = numpy.sum(
        x203 * (x128 * x135 * x195 + x128 * x136 * x194 + x135 * x152 * x194)
    )
    result[1, 5] = numpy.sum(
        x197 * (x128 * x160 * x188 + x128 * x165 * x187 + x152 * x160 * x187)
    )
    result[1, 6] = numpy.sum(
        x191 * (x108 * x147 * x208 + x119 * x147 * x212 + x151 * x208)
    )
    result[1, 7] = numpy.sum(
        x197 * (x135 * x147 * x201 + x135 * x151 * x200 + x136 * x147 * x200)
    )
    result[1, 8] = numpy.sum(
        x197 * (x147 * x160 * x195 + x147 * x165 * x194 + x151 * x160 * x194)
    )
    result[1, 9] = numpy.sum(
        x191 * (x147 * x183 * x188 + x147 * x185 * x187 + x151 * x183 * x187)
    )
    result[2, 0] = numpy.sum(
        x191 * (x103 * x106 * x190 + x106 * x189 + x112 * x214 * x93)
    )
    result[2, 1] = numpy.sum(
        x197 * (x111 * x130 * x215 + x113 * x214 * x82 + x114 * x213 * x82)
    )
    result[2, 2] = numpy.sum(
        x197 * (x103 * x218 * x220 + x112 * x130 * x218 + x219 * x220)
    )
    result[2, 3] = numpy.sum(
        x197 * (x128 * x141 * x214 + x128 * x146 * x213 + x141 * x152 * x213)
    )
    result[2, 4] = numpy.sum(
        x203 * (x113 * x128 * x219 + x113 * x152 * x218 + x114 * x128 * x218)
    )
    result[2, 5] = numpy.sum(
        x197 * (x103 * x223 * x225 + x112 * x152 * x223 + x224 * x225)
    )
    result[2, 6] = numpy.sum(
        x191 * (x147 * x172 * x214 + x147 * x174 * x213 + x151 * x172 * x213)
    )
    result[2, 7] = numpy.sum(
        x197 * (x141 * x147 * x219 + x141 * x151 * x218 + x146 * x147 * x218)
    )
    result[2, 8] = numpy.sum(
        x197 * (x113 * x147 * x224 + x113 * x151 * x223 + x114 * x147 * x223)
    )
    result[2, 9] = numpy.sum(
        x191 * (x103 * x147 * x230 + x112 * x147 * x234 + x151 * x230)
    )
    result[3, 0] = numpy.sum(
        x242 * (x108 * x237 * x241 + x119 * x237 * x91 + x240 * x241)
    )
    result[3, 1] = numpy.sum(x252 * (x108 * x246 * x74 + x119 * x251 * x74 + x246 * x76))
    result[3, 2] = numpy.sum(
        x252 * (x135 * x237 * x76 + x135 * x240 * x74 + x136 * x237 * x74)
    )
    result[3, 3] = numpy.sum(
        x252 * (x108 * x259 * x262 + x119 * x126 * x259 + x261 * x262)
    )
    result[3, 4] = numpy.sum(
        x263 * (x126 * x135 * x245 + x135 * x251 * x53 + x136 * x245 * x53)
    )
    result[3, 5] = numpy.sum(
        x252 * (x126 * x160 * x237 + x160 * x240 * x53 + x165 * x237 * x53)
    )
    result[3, 6] = numpy.sum(
        x242 * (x108 * x121 * x269 + x119 * x121 * x276 + x150 * x269)
    )
    result[3, 7] = numpy.sum(
        x252 * (x121 * x135 * x261 + x121 * x136 * x259 + x135 * x150 * x259)
    )
    result[3, 8] = numpy.sum(
        x252 * (x121 * x160 * x251 + x121 * x165 * x245 + x150 * x160 * x245)
    )
    result[3, 9] = numpy.sum(
        x242 * (x121 * x183 * x240 + x121 * x185 * x237 + x150 * x183 * x237)
    )
    result[4, 0] = numpy.sum(
        x197 * (x101 * x215 * x91 + x187 * x214 * x88 + x188 * x213 * x88)
    )
    result[4, 1] = numpy.sum(
        x277 * (x194 * x213 * x76 + x194 * x214 * x74 + x195 * x213 * x74)
    )
    result[4, 2] = numpy.sum(
        x277 * (x187 * x218 * x76 + x187 * x219 * x74 + x188 * x218 * x74)
    )
    result[4, 3] = numpy.sum(
        x277 * (x126 * x200 * x213 + x200 * x214 * x53 + x201 * x213 * x53)
    )
    result[4, 4] = numpy.sum(
        x278 * (x126 * x194 * x218 + x194 * x219 * x53 + x195 * x218 * x53)
    )
    result[4, 5] = numpy.sum(
        x277 * (x126 * x187 * x223 + x187 * x224 * x53 + x188 * x223 * x53)
    )
    result[4, 6] = numpy.sum(
        x197 * (x121 * x207 * x214 + x121 * x212 * x213 + x150 * x207 * x213)
    )
    result[4, 7] = numpy.sum(
        x277 * (x121 * x200 * x219 + x121 * x201 * x218 + x150 * x200 * x218)
    )
    result[4, 8] = numpy.sum(
        x277 * (x121 * x194 * x224 + x121 * x195 * x223 + x150 * x194 * x223)
    )
    result[4, 9] = numpy.sum(
        x197 * (x121 * x187 * x234 + x121 * x188 * x229 + x150 * x187 * x229)
    )
    result[5, 0] = numpy.sum(
        x242 * (x103 * x281 * x285 + x112 * x281 * x91 + x284 * x285)
    )
    result[5, 1] = numpy.sum(
        x252 * (x113 * x281 * x76 + x113 * x284 * x74 + x114 * x281 * x74)
    )
    result[5, 2] = numpy.sum(x252 * (x103 * x289 * x74 + x112 * x294 * x74 + x289 * x76))
    result[5, 3] = numpy.sum(
        x252 * (x126 * x141 * x281 + x141 * x284 * x53 + x146 * x281 * x53)
    )
    result[5, 4] = numpy.sum(
        x263 * (x113 * x126 * x288 + x113 * x294 * x53 + x114 * x288 * x53)
    )
    result[5, 5] = numpy.sum(
        x252 * (x103 * x301 * x304 + x112 * x126 * x301 + x303 * x304)
    )
    result[5, 6] = numpy.sum(
        x242 * (x121 * x172 * x284 + x121 * x174 * x281 + x150 * x172 * x281)
    )
    result[5, 7] = numpy.sum(
        x252 * (x121 * x141 * x294 + x121 * x146 * x288 + x141 * x150 * x288)
    )
    result[5, 8] = numpy.sum(
        x252 * (x113 * x121 * x303 + x113 * x150 * x301 + x114 * x121 * x301)
    )
    result[5, 9] = numpy.sum(
        x242 * (x103 * x121 * x310 + x112 * x121 * x316 + x150 * x310)
    )
    result[6, 0] = numpy.sum(
        x191 * (x108 * x317 * x319 + x119 * x317 * x65 + x318 * x319)
    )
    result[6, 1] = numpy.sum(x197 * (x108 * x323 * x37 + x119 * x327 * x37 + x323 * x69))
    result[6, 2] = numpy.sum(
        x197 * (x135 * x317 * x69 + x135 * x318 * x37 + x136 * x317 * x37)
    )
    result[6, 3] = numpy.sum(x197 * (x108 * x331 * x48 + x119 * x335 * x48 + x331 * x49))
    result[6, 4] = numpy.sum(
        x203 * (x135 * x322 * x49 + x135 * x327 * x48 + x136 * x322 * x48)
    )
    result[6, 5] = numpy.sum(
        x197 * (x160 * x317 * x49 + x160 * x318 * x48 + x165 * x317 * x48)
    )
    result[6, 6] = numpy.sum(x191 * (x108 * x337 * x342 + x119 * x337 * x44 + x341 * x5))
    result[6, 7] = numpy.sum(
        x197 * (x134 * x335 * x342 + x135 * x330 * x44 + x136 * x330 * x43)
    )
    result[6, 8] = numpy.sum(
        x197 * (x160 * x322 * x44 + x160 * x327 * x43 + x165 * x322 * x43)
    )
    result[6, 9] = numpy.sum(
        x191 * (x183 * x317 * x44 + x183 * x318 * x43 + x185 * x317 * x43)
    )
    result[7, 0] = numpy.sum(
        x197 * (x213 * x237 * x65 + x213 * x240 * x64 + x214 * x237 * x64)
    )
    result[7, 1] = numpy.sum(
        x277 * (x213 * x245 * x69 + x213 * x251 * x37 + x214 * x245 * x37)
    )
    result[7, 2] = numpy.sum(
        x277 * (x218 * x237 * x69 + x218 * x240 * x37 + x219 * x237 * x37)
    )
    result[7, 3] = numpy.sum(
        x277 * (x213 * x259 * x49 + x213 * x261 * x48 + x214 * x259 * x48)
    )
    result[7, 4] = numpy.sum(
        x278 * (x218 * x245 * x49 + x218 * x251 * x48 + x219 * x245 * x48)
    )
    result[7, 5] = numpy.sum(
        x277 * (x223 * x237 * x49 + x223 * x240 * x48 + x224 * x237 * x48)
    )
    result[7, 6] = numpy.sum(
        x197 * (x106 * x276 * x342 + x213 * x268 * x44 + x214 * x268 * x43)
    )
    result[7, 7] = numpy.sum(
        x277 * (x218 * x259 * x44 + x218 * x261 * x43 + x219 * x259 * x43)
    )
    result[7, 8] = numpy.sum(
        x277 * (x223 * x245 * x44 + x223 * x251 * x43 + x224 * x245 * x43)
    )
    result[7, 9] = numpy.sum(
        x197 * (x229 * x237 * x44 + x229 * x240 * x43 + x234 * x237 * x43)
    )
    result[8, 0] = numpy.sum(
        x197 * (x187 * x281 * x65 + x187 * x284 * x64 + x188 * x281 * x64)
    )
    result[8, 1] = numpy.sum(
        x277 * (x194 * x281 * x69 + x194 * x284 * x37 + x195 * x281 * x37)
    )
    result[8, 2] = numpy.sum(
        x277 * (x187 * x288 * x69 + x187 * x294 * x37 + x188 * x288 * x37)
    )
    result[8, 3] = numpy.sum(
        x277 * (x200 * x281 * x49 + x200 * x284 * x48 + x201 * x281 * x48)
    )
    result[8, 4] = numpy.sum(
        x278 * (x194 * x288 * x49 + x194 * x294 * x48 + x195 * x288 * x48)
    )
    result[8, 5] = numpy.sum(
        x277 * (x187 * x301 * x49 + x187 * x303 * x48 + x188 * x301 * x48)
    )
    result[8, 6] = numpy.sum(
        x197 * (x207 * x281 * x44 + x207 * x284 * x43 + x212 * x281 * x43)
    )
    result[8, 7] = numpy.sum(
        x277 * (x200 * x288 * x44 + x200 * x294 * x43 + x201 * x288 * x43)
    )
    result[8, 8] = numpy.sum(
        x277 * (x194 * x301 * x44 + x194 * x303 * x43 + x195 * x301 * x43)
    )
    result[8, 9] = numpy.sum(
        x197 * (x101 * x316 * x344 + x187 * x309 * x44 + x188 * x309 * x43)
    )
    result[9, 0] = numpy.sum(
        x191 * (x103 * x345 * x347 + x112 * x345 * x65 + x346 * x347)
    )
    result[9, 1] = numpy.sum(
        x197 * (x113 * x345 * x69 + x113 * x346 * x37 + x114 * x345 * x37)
    )
    result[9, 2] = numpy.sum(x197 * (x103 * x351 * x37 + x112 * x355 * x37 + x351 * x69))
    result[9, 3] = numpy.sum(
        x197 * (x141 * x345 * x49 + x141 * x346 * x48 + x146 * x345 * x48)
    )
    result[9, 4] = numpy.sum(
        x203 * (x113 * x350 * x49 + x113 * x355 * x48 + x114 * x350 * x48)
    )
    result[9, 5] = numpy.sum(x197 * (x103 * x359 * x48 + x112 * x363 * x48 + x359 * x49))
    result[9, 6] = numpy.sum(
        x191 * (x172 * x345 * x44 + x172 * x346 * x43 + x174 * x345 * x43)
    )
    result[9, 7] = numpy.sum(
        x197 * (x141 * x350 * x44 + x141 * x355 * x43 + x146 * x350 * x43)
    )
    result[9, 8] = numpy.sum(
        x197 * (x111 * x344 * x363 + x113 * x358 * x44 + x114 * x358 * x43)
    )
    result[9, 9] = numpy.sum(x191 * (x103 * x344 * x365 + x112 * x365 * x44 + x368 * x5))
    result[10, 0] = numpy.sum(x110 * (x108 * x370 * x40 + x119 * x371 * x40 + x370 * x59))
    result[10, 1] = numpy.sum(x133 * (x108 * x24 * x373 + x119 * x24 * x374 + x29 * x373))
    result[10, 2] = numpy.sum(
        x133 * (x135 * x24 * x371 + x135 * x29 * x369 + x136 * x24 * x369)
    )
    result[10, 3] = numpy.sum(
        x133 * (x108 * x375 * x377 + x119 * x21 * x375 + x17 * x376)
    )
    result[10, 4] = numpy.sum(
        x155 * (x134 * x374 * x377 + x135 * x21 * x372 + x136 * x18 * x372)
    )
    result[10, 5] = numpy.sum(
        x133 * (x160 * x18 * x371 + x160 * x21 * x369 + x165 * x18 * x369)
    )
    result[10, 6] = numpy.sum(
        x110
        * (
            x108 * x379
            + x340
            * (
                x101 * x339
                - 2.0 * x19 * (-ax * x378 + 2.0 * x266 + 2.0 * x267)
                + 3.0 * x2 * (x271 + x272 + x275 + x332 + x333 + x334)
            )
            + x379 * x9
        )
    )
    result[10, 7] = numpy.sum(
        x133 * (x134 * x375 * x380 + x134 * x376 + x136 * x14 * x375)
    )
    result[10, 8] = numpy.sum(x133 * (x14 * x160 * x374 + x160 * x381 * x9 + x165 * x381))
    result[10, 9] = numpy.sum(x110 * (x14 * x183 * x371 + x183 * x382 * x9 + x185 * x382))
    result[11, 0] = numpy.sum(
        x191 * (x213 * x317 * x59 + x213 * x318 * x40 + x214 * x317 * x40)
    )
    result[11, 1] = numpy.sum(
        x197 * (x213 * x24 * x327 + x213 * x29 * x322 + x214 * x24 * x322)
    )
    result[11, 2] = numpy.sum(
        x197 * (x218 * x24 * x318 + x218 * x29 * x317 + x219 * x24 * x317)
    )
    result[11, 3] = numpy.sum(
        x197 * (x106 * x335 * x377 + x18 * x214 * x330 + x21 * x213 * x330)
    )
    result[11, 4] = numpy.sum(
        x203 * (x18 * x218 * x327 + x18 * x219 * x322 + x21 * x218 * x322)
    )
    result[11, 5] = numpy.sum(
        x197 * (x18 * x223 * x318 + x18 * x224 * x317 + x21 * x223 * x317)
    )
    result[11, 6] = numpy.sum(
        x191 * (x106 * x337 * x380 + x106 * x341 + x14 * x214 * x337)
    )
    result[11, 7] = numpy.sum(x197 * (x14 * x218 * x335 + x218 * x383 * x9 + x219 * x383))
    result[11, 8] = numpy.sum(x197 * (x14 * x223 * x327 + x223 * x384 * x9 + x224 * x384))
    result[11, 9] = numpy.sum(x191 * (x14 * x234 * x317 + x317 * x385 * x9 + x318 * x385))
    result[12, 0] = numpy.sum(
        x242 * (x237 * x281 * x59 + x237 * x284 * x40 + x240 * x281 * x40)
    )
    result[12, 1] = numpy.sum(
        x252 * (x24 * x245 * x284 + x24 * x251 * x281 + x245 * x281 * x29)
    )
    result[12, 2] = numpy.sum(
        x252 * (x237 * x24 * x294 + x237 * x288 * x29 + x24 * x240 * x288)
    )
    result[12, 3] = numpy.sum(
        x252 * (x18 * x259 * x284 + x18 * x261 * x281 + x21 * x259 * x281)
    )
    result[12, 4] = numpy.sum(
        x263 * (x18 * x245 * x294 + x18 * x251 * x288 + x21 * x245 * x288)
    )
    result[12, 5] = numpy.sum(
        x252 * (x18 * x237 * x303 + x18 * x240 * x301 + x21 * x237 * x301)
    )
    result[12, 6] = numpy.sum(x242 * (x14 * x276 * x281 + x281 * x386 * x9 + x284 * x386))
    result[12, 7] = numpy.sum(x252 * (x14 * x259 * x294 + x259 * x387 * x9 + x261 * x387))
    result[12, 8] = numpy.sum(x252 * (x14 * x251 * x301 + x301 * x388 * x9 + x303 * x388))
    result[12, 9] = numpy.sum(x242 * (x14 * x237 * x316 + x237 * x389 * x9 + x240 * x389))
    result[13, 0] = numpy.sum(
        x191 * (x187 * x345 * x59 + x187 * x346 * x40 + x188 * x345 * x40)
    )
    result[13, 1] = numpy.sum(
        x197 * (x194 * x24 * x346 + x194 * x29 * x345 + x195 * x24 * x345)
    )
    result[13, 2] = numpy.sum(
        x197 * (x187 * x24 * x355 + x187 * x29 * x350 + x188 * x24 * x350)
    )
    result[13, 3] = numpy.sum(
        x197 * (x18 * x200 * x346 + x18 * x201 * x345 + x200 * x21 * x345)
    )
    result[13, 4] = numpy.sum(
        x203 * (x18 * x194 * x355 + x18 * x195 * x350 + x194 * x21 * x350)
    )
    result[13, 5] = numpy.sum(
        x197 * (x17 * x363 * x390 + x18 * x188 * x358 + x187 * x21 * x358)
    )
    result[13, 6] = numpy.sum(x191 * (x14 * x212 * x345 + x345 * x391 * x9 + x346 * x391))
    result[13, 7] = numpy.sum(x197 * (x14 * x200 * x355 + x200 * x392 * x9 + x201 * x392))
    result[13, 8] = numpy.sum(x197 * (x14 * x194 * x363 + x194 * x393 * x9 + x195 * x393))
    result[13, 9] = numpy.sum(x191 * (x101 * x368 + x14 * x188 * x365 + x365 * x390 * x9))
    result[14, 0] = numpy.sum(x110 * (x103 * x395 * x40 + x112 * x396 * x40 + x395 * x59))
    result[14, 1] = numpy.sum(
        x133 * (x113 * x24 * x396 + x113 * x29 * x394 + x114 * x24 * x394)
    )
    result[14, 2] = numpy.sum(x133 * (x103 * x24 * x398 + x112 * x24 * x399 + x29 * x398))
    result[14, 3] = numpy.sum(
        x133 * (x141 * x18 * x396 + x141 * x21 * x394 + x146 * x18 * x394)
    )
    result[14, 4] = numpy.sum(
        x155 * (x111 * x399 * x400 + x113 * x21 * x397 + x114 * x18 * x397)
    )
    result[14, 5] = numpy.sum(
        x133 * (x103 * x400 * x401 + x112 * x21 * x401 + x17 * x402)
    )
    result[14, 6] = numpy.sum(x110 * (x14 * x172 * x396 + x172 * x403 * x9 + x174 * x403))
    result[14, 7] = numpy.sum(x133 * (x14 * x141 * x399 + x141 * x404 * x9 + x146 * x404))
    result[14, 8] = numpy.sum(
        x133 * (x111 * x343 * x401 * x9 + x111 * x402 + x114 * x14 * x401)
    )
    result[14, 9] = numpy.sum(
        x110
        * (
            x103 * x406
            + x343
            * (
                x106 * x367
                - 2.0 * x19 * (-ax * x405 + 2.0 * x307 + 2.0 * x308)
                + 3.0 * x2 * (x311 + x312 + x315 + x360 + x361 + x362)
            )
            + x406 * x9
        )
    )
    return result


def kinetic3d_44(ax, da, A, bx, db, B):
    """Cartesian 3D (gg) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((15, 15), dtype=float)

    x0 = 2.0 * ax
    x1 = 2.0 * bx
    x2 = (x0 + x1) ** (-1.0)
    x3 = (ax + bx) ** (-1.0)
    x4 = -x3 * (ax * A[0] + bx * B[0])
    x5 = -x4 - A[0]
    x6 = -ax
    x7 = x5**2
    x8 = 2.0 * ax**2
    x9 = -x6 - x8 * (x2 + x7)
    x10 = -x4 - B[0]
    x11 = ax * x3
    x12 = bx * x11
    x13 = numpy.exp(-x12 * (A[0] - B[0]) ** 2)
    x14 = 1.77245385090552 * numpy.sqrt(x3)
    x15 = x13 * x14
    x16 = 2.0 * x15
    x17 = x10 * x16
    x18 = x10 * x15
    x19 = 4.0 * x12
    x20 = x2 * (x17 * x9 + x18 * x19)
    x21 = x15 * x2
    x22 = x21 * x9
    x23 = bx * x3
    x24 = x0 * x23
    x25 = x18 * (x24 + x9)
    x26 = x10 * x25
    x27 = x10**2 * x15
    x28 = x21 + x27
    x29 = -x16
    x30 = x11 * (x1 * x28 + x29)
    x31 = x26 + x30
    x32 = x22 + x31
    x33 = x10 * x32
    x34 = x10 * x28
    x35 = 2.0 * x21
    x36 = x10 * x35
    x37 = x34 + x36
    x38 = x11 * (x1 * x37 - 3.0 * x18)
    x39 = x33 + x38
    x40 = x20 + x39
    x41 = x40 * x5
    x42 = 3.0 * x22
    x43 = x2 * (3.0 * x26 + 3.0 * x30 + x42)
    x44 = 3.0 * x21
    x45 = x2 * (3.0 * x27 + x44)
    x46 = x37 * x5
    x47 = x45 + x46
    x48 = x12 * x47
    x49 = x10 * x37
    x50 = x45 + x49
    x51 = 4.0 * x21
    x52 = x10 * x40 - x11 * (-2.0 * bx * x50 + 4.0 * x27 + x51)
    x53 = x2 * (4.0 * x41 + 5.0 * x43 + 8.0 * x48 + x52)
    x54 = 4.0 * x20
    x55 = x2 * (4.0 * x33 + 4.0 * x38 + x54)
    x56 = x43 + x52
    x57 = x5 * x56
    x58 = x10 * x21
    x59 = 8.0 * x58
    x60 = x2 * (4.0 * x34 + x59)
    x61 = x5 * x50
    x62 = x60 + x61
    x63 = x24 * x62 + x55 + x57
    x64 = x5 * x63
    x65 = x25 * x5
    x66 = x15 * x5
    x67 = x10 * x66
    x68 = x21 + x67
    x69 = x19 * x68 + x42 + 2.0 * x65
    x70 = x2 * (x31 + x69)
    x71 = x32 * x5
    x72 = x28 * x5
    x73 = x36 + x72
    x74 = x20 + x24 * x73 + x71
    x75 = x5 * x74
    x76 = x17 * x5 + x44
    x77 = x2 * (x27 + x76)
    x78 = x5 * x73
    x79 = x77 + x78
    x80 = x23 * (2.0 * ax * x79 - 2.0 * x27 - x35)
    x81 = 3.0 * x70 + 3.0 * x75 + 3.0 * x80
    x82 = x2 * (2.0 * x41 + 2.0 * x43 + 4.0 * x48 + x81)
    x83 = 6.0 * x12
    x84 = x2 * (x39 + x54 + 3.0 * x71 + x73 * x83)
    x85 = x24 * x47 + x41 + x43
    x86 = x5 * x85
    x87 = 3.0 * x72
    x88 = x2 * (x34 + x59 + x87)
    x89 = x47 * x5
    x90 = x88 + x89
    x91 = x10 * x51
    x92 = x23 * (2.0 * ax * x90 - 2.0 * x34 - x91)
    x93 = x84 + x86 + x92
    x94 = x5 * x93
    x95 = x2 * (5.0 * x45 + 4.0 * x46 + x49)
    x96 = x5 * x62
    x97 = x95 + x96
    x98 = 2.0 * x45
    x99 = x23 * (2.0 * ax * x97 - 2.0 * x49 - x98)
    x100 = 3.0 * x77 + 3.0 * x78
    x101 = x2 * (x100 + 2.0 * x46 + x98)
    x102 = x5 * x90
    x103 = x101 + x102
    x104 = x23 * (2.0 * ax * x103 - 3.0 * x45 - 3.0 * x46)
    x105 = x53 + x64 + x99
    x106 = 4.0 * x88 + 4.0 * x89
    x107 = x2 * (x106 + 2.0 * x60 + 2.0 * x61) + x5 * x97
    x108 = (
        x105 * x5
        + x2 * (x19 * x62 + 2.0 * x55 + 2.0 * x57 + 4.0 * x84 + 4.0 * x86 + 4.0 * x92)
        - x23 * (-2.0 * ax * x107 + 3.0 * x60 + 3.0 * x61)
    )
    x109 = x107 * x5 + x2 * (4.0 * x101 + 4.0 * x102 + 3.0 * x95 + 3.0 * x96)
    x110 = numpy.exp(-x12 * (A[1] - B[1]) ** 2)
    x111 = numpy.exp(-x12 * (A[2] - B[2]) ** 2)
    x112 = 3.14159265358979 * x111 * x3
    x113 = x110 * x112
    x114 = -x3 * (ax * A[1] + bx * B[1])
    x115 = -x114 - A[1]
    x116 = x115**2
    x117 = -x6 - x8 * (x116 + x2)
    x118 = x109 * x113
    x119 = -x3 * (ax * A[2] + bx * B[2])
    x120 = -x119 - A[2]
    x121 = x120**2
    x122 = -x6 - x8 * (x121 + x2)
    x123 = 0.179587122125167 * da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**5.5)
    x124 = 6.89597470414309 * x123
    x125 = -x114 - B[1]
    x126 = x110 * x14
    x127 = x125 * x126
    x128 = x127 * (x117 + x24)
    x129 = x2 * (x18 + x66)
    x130 = x5 * x68
    x131 = x2 * (2.0 * x129 + 2.0 * x130 + 2.0 * x72 + x91)
    x132 = x5 * x79
    x133 = x103 * x5 + 3.0 * x2 * (x131 + x132 + x88 + x89)
    x134 = x111 * x14
    x135 = x66 * (x24 + x9)
    x136 = x2 * (x135 + x25)
    x137 = x22 + x24 * x68 + x65
    x138 = x137 * x5
    x139 = x129 + x130
    x140 = x0 * x139 - x17
    x141 = x1 * x3
    x142 = x2 * (
        2.0 * x136 + 2.0 * x138 + x140 * x141 + x19 * x73 + 2.0 * x20 + 2.0 * x71
    )
    x143 = x70 + x75 + x80
    x144 = x143 * x5
    x145 = x131 + x132
    x146 = x23 * (2.0 * ax * x145 - 6.0 * x58 - x87)
    x147 = x104 + x82 + x94
    x148 = x113 * (
        x147 * x5
        + 3.0 * x2 * (x142 + x144 + x146 + x84 + x86 + x92)
        + x23 * (2.0 * ax * x133 - x106)
    )
    x149 = x113 * x133
    x150 = 18.2450341145548 * x123
    x151 = -x119 - B[2]
    x152 = x134 * x151
    x153 = x152 * (x122 + x24)
    x154 = x126 * x2
    x155 = x117 * x154
    x156 = x125 * x128
    x157 = x125**2 * x126
    x158 = x154 + x157
    x159 = 2.0 * x126
    x160 = -x159
    x161 = x11 * (x1 * x158 + x160)
    x162 = x156 + x161
    x163 = x155 + x162
    x164 = x15 * x7
    x165 = x2 * (x164 + x76)
    x166 = x139 * x5
    x167 = x145 * x5 + x2 * (x100 + 2.0 * x165 + 2.0 * x166)
    x168 = x134 * x167
    x169 = x135 * x5
    x170 = x164 + x21
    x171 = x23 * (x0 * x170 + x29)
    x172 = x169 + x171
    x173 = x2 * (x172 + x69)
    x174 = x140 * x23
    x175 = x136 + x138 + x174
    x176 = x175 * x5
    x177 = x165 + x166
    x178 = 2.0 * ax * x177 - x44 - 3.0 * x67
    x179 = x142 + x144 + x146
    x180 = (
        x179 * x5
        + x2 * (x141 * x178 + 2.0 * x173 + 2.0 * x176 + x81)
        - 2.0 * x23 * (-ax * x167 + 2.0 * x77 + 2.0 * x78)
    )
    x181 = 23.5542377588857 * x123
    x182 = x113 * x151
    x183 = 40.7971365319473 * x123
    x184 = x134 * x2
    x185 = x122 * x184
    x186 = x151 * x153
    x187 = x134 * x151**2
    x188 = x184 + x187
    x189 = 2.0 * x134
    x190 = -x189
    x191 = x11 * (x1 * x188 + x190)
    x192 = x186 + x191
    x193 = x185 + x192
    x194 = x126 * x167
    x195 = x125 * x159
    x196 = x2 * (x117 * x195 + x127 * x19)
    x197 = x125 * x163
    x198 = x125 * x158
    x199 = 2.0 * x154
    x200 = x125 * x199
    x201 = x198 + x200
    x202 = x11 * (x1 * x201 - 3.0 * x127)
    x203 = x197 + x202
    x204 = x196 + x203
    x205 = x5 * (x170 + x35)
    x206 = x177 * x5 + x2 * (3.0 * x129 + 3.0 * x130 + x205)
    x207 = x134 * x206
    x208 = x172 + x22
    x209 = x2 * (x16 * x5 * x9 + x19 * x66) + x208 * x5 + x23 * (x0 * x205 - 3.0 * x66)
    x210 = x173 + x176 + x178 * x23
    x211 = (
        x2 * (3.0 * x136 + 3.0 * x138 + 3.0 * x174 + x209)
        + x210 * x5
        - 2.0 * x23 * (-ax * x206 + 2.0 * x129 + 2.0 * x130)
    )
    x212 = x151 * x189
    x213 = x2 * (x122 * x212 + x152 * x19)
    x214 = x151 * x193
    x215 = x151 * x188
    x216 = 2.0 * x184
    x217 = x151 * x216
    x218 = x215 + x217
    x219 = x11 * (x1 * x218 - 3.0 * x152)
    x220 = x214 + x219
    x221 = x213 + x220
    x222 = x126 * x206
    x223 = 3.0 * x155
    x224 = x2 * (3.0 * x156 + 3.0 * x161 + x223)
    x225 = 3.0 * x154
    x226 = x2 * (3.0 * x157 + x225)
    x227 = x125 * x201
    x228 = x226 + x227
    x229 = 4.0 * x154
    x230 = -x11 * (-2.0 * bx * x228 + 4.0 * x157 + x229) + x125 * x204
    x231 = x224 + x230
    x232 = x2 * (3.0 * x164 + x44) + x205 * x5
    x233 = x134 * x232
    x234 = (
        x2 * (3.0 * x169 + 3.0 * x171 + x42)
        + x209 * x5
        - x23 * (-2.0 * ax * x232 + 4.0 * x164 + x51)
    )
    x235 = 3.0 * x185
    x236 = x2 * (3.0 * x186 + 3.0 * x191 + x235)
    x237 = 3.0 * x184
    x238 = x2 * (3.0 * x187 + x237)
    x239 = x151 * x218
    x240 = x238 + x239
    x241 = 4.0 * x184
    x242 = -x11 * (-2.0 * bx * x240 + 4.0 * x187 + x241) + x151 * x221
    x243 = x236 + x242
    x244 = x126 * x232
    x245 = x115 * x126
    x246 = x245 * (x117 + x24)
    x247 = x108 * x113
    x248 = x107 * x113
    x249 = x115 * x128
    x250 = x125 * x245
    x251 = x154 + x250
    x252 = x155 + x24 * x251 + x249
    x253 = x103 * x134
    x254 = 48.2718229290016 * x123
    x255 = x115 * x163
    x256 = x115 * x158
    x257 = x200 + x256
    x258 = x196 + x24 * x257 + x255
    x259 = x134 * x145
    x260 = 62.3186554316989 * x123
    x261 = 107.939077467081 * x123
    x262 = x115 * x201
    x263 = x226 + x262
    x264 = x134 * x263
    x265 = x115 * x204
    x266 = x224 + x24 * x263 + x265
    x267 = x125 * x154
    x268 = 8.0 * x267
    x269 = x2 * (4.0 * x198 + x268)
    x270 = x115 * x228
    x271 = x269 + x270
    x272 = x134 * x271
    x273 = 4.0 * x196
    x274 = x2 * (4.0 * x197 + 4.0 * x202 + x273)
    x275 = x115 * x231
    x276 = x24 * x271 + x274 + x275
    x277 = x120 * x134
    x278 = x277 * (x122 + x24)
    x279 = x113 * x120
    x280 = x120 * x153
    x281 = x151 * x277
    x282 = x184 + x281
    x283 = x185 + x24 * x282 + x280
    x284 = x103 * x126
    x285 = x120 * x193
    x286 = x120 * x188
    x287 = x217 + x286
    x288 = x213 + x24 * x287 + x285
    x289 = x126 * x145
    x290 = x120 * x218
    x291 = x238 + x290
    x292 = x126 * x291
    x293 = x120 * x221
    x294 = x236 + x24 * x291 + x293
    x295 = x151 * x184
    x296 = 8.0 * x295
    x297 = x2 * (4.0 * x215 + x296)
    x298 = x120 * x240
    x299 = x297 + x298
    x300 = x126 * x299
    x301 = 4.0 * x213
    x302 = x2 * (4.0 * x214 + 4.0 * x219 + x301)
    x303 = x120 * x243
    x304 = x24 * x299 + x302 + x303
    x305 = x115 * x246
    x306 = x116 * x126
    x307 = x154 + x306
    x308 = x23 * (x0 * x307 + x160)
    x309 = x305 + x308
    x310 = x155 + x309
    x311 = x134 * x97
    x312 = x2 * (x127 + x245)
    x313 = x115 * x251
    x314 = x312 + x313
    x315 = x134 * x314
    x316 = x2 * (x128 + x246)
    x317 = x115 * x252
    x318 = x0 * x314 - x195
    x319 = x23 * x318
    x320 = x316 + x317 + x319
    x321 = x19 * x251 + x223 + 2.0 * x249
    x322 = x2 * (x162 + x321)
    x323 = x115 * x258
    x324 = x115 * x195 + x225
    x325 = x2 * (x157 + x324)
    x326 = x115 * x257
    x327 = x325 + x326
    x328 = x23 * (2.0 * ax * x327 - 2.0 * x157 - x199)
    x329 = x322 + x323 + x328
    x330 = x134 * x79
    x331 = 80.4530382150027 * x123
    x332 = 139.348749811665 * x123
    x333 = x2 * (x203 + 3.0 * x255 + x257 * x83 + x273)
    x334 = x115 * x266
    x335 = 3.0 * x256
    x336 = x2 * (x198 + x268 + x335)
    x337 = x115 * x263
    x338 = x336 + x337
    x339 = x125 * x229
    x340 = x23 * (2.0 * ax * x338 - 2.0 * x198 - x339)
    x341 = x333 + x334 + x340
    x342 = x134 * x139
    x343 = x2 * (5.0 * x226 + x227 + 4.0 * x262)
    x344 = x115 * x271
    x345 = x343 + x344
    x346 = x134 * x345
    x347 = 8.0 * x12
    x348 = x2 * (5.0 * x224 + x230 + x263 * x347 + 4.0 * x265)
    x349 = x115 * x276
    x350 = 2.0 * x226
    x351 = x23 * (2.0 * ax * x345 - 2.0 * x227 - x350)
    x352 = x348 + x349 + x351
    x353 = 241.359114645008 * x123
    x354 = x120 * x278
    x355 = x121 * x134
    x356 = x184 + x355
    x357 = x23 * (x0 * x356 + x190)
    x358 = x354 + x357
    x359 = x185 + x358
    x360 = x126 * x97
    x361 = x2 * (x152 + x277)
    x362 = x120 * x282
    x363 = x361 + x362
    x364 = x126 * x363
    x365 = x2 * (x153 + x278)
    x366 = x120 * x283
    x367 = x0 * x363 - x212
    x368 = x23 * x367
    x369 = x365 + x366 + x368
    x370 = x19 * x282 + x235 + 2.0 * x280
    x371 = x2 * (x192 + x370)
    x372 = x120 * x288
    x373 = x120 * x212 + x237
    x374 = x2 * (x187 + x373)
    x375 = x120 * x287
    x376 = x374 + x375
    x377 = x23 * (2.0 * ax * x376 - 2.0 * x187 - x216)
    x378 = x371 + x372 + x377
    x379 = x126 * x79
    x380 = x2 * (x220 + 3.0 * x285 + x287 * x83 + x301)
    x381 = x120 * x294
    x382 = 3.0 * x286
    x383 = x2 * (x215 + x296 + x382)
    x384 = x120 * x291
    x385 = x383 + x384
    x386 = x151 * x241
    x387 = x23 * (2.0 * ax * x385 - 2.0 * x215 - x386)
    x388 = x380 + x381 + x387
    x389 = x126 * x139
    x390 = x2 * (5.0 * x238 + x239 + 4.0 * x290)
    x391 = x120 * x299
    x392 = x390 + x391
    x393 = x126 * x392
    x394 = x2 * (5.0 * x236 + x242 + x291 * x347 + 4.0 * x293)
    x395 = x120 * x304
    x396 = 2.0 * x238
    x397 = x23 * (2.0 * ax * x392 - 2.0 * x239 - x396)
    x398 = x394 + x395 + x397
    x399 = x115 * (x199 + x307)
    x400 = (
        x115 * x310
        + x2 * (x115 * x117 * x159 + x19 * x245)
        + x23 * (x0 * x399 - 3.0 * x245)
    )
    x401 = x134 * x62
    x402 = x2 * (x309 + x321)
    x403 = x115 * x320
    x404 = x2 * (x306 + x324)
    x405 = x115 * x314
    x406 = x404 + x405
    x407 = 2.0 * ax * x406 - x225 - 3.0 * x250
    x408 = x23 * x407 + x402 + x403
    x409 = x134 * x47
    x410 = x2 * (2.0 * x256 + 2.0 * x312 + 2.0 * x313 + x339)
    x411 = x115 * x327
    x412 = x410 + x411
    x413 = x134 * x412
    x414 = x2 * (
        x141 * x318 + x19 * x257 + 2.0 * x196 + 2.0 * x255 + 2.0 * x316 + 2.0 * x317
    )
    x415 = x115 * x329
    x416 = x23 * (2.0 * ax * x412 - 6.0 * x267 - x335)
    x417 = x414 + x415 + x416
    x418 = 3.0 * x325 + 3.0 * x326
    x419 = x2 * (2.0 * x262 + x350 + x418)
    x420 = x115 * x338
    x421 = x419 + x420
    x422 = x134 * x421
    x423 = 3.0 * x322 + 3.0 * x323 + 3.0 * x328
    x424 = x2 * (x19 * x263 + 2.0 * x224 + 2.0 * x265 + x423)
    x425 = x115 * x341
    x426 = x23 * (2.0 * ax * x421 - 3.0 * x226 - 3.0 * x262)
    x427 = x424 + x425 + x426
    x428 = 4.0 * x336 + 4.0 * x337
    x429 = x115 * x345 + x2 * (2.0 * x269 + 2.0 * x270 + x428)
    x430 = (
        x115 * x352
        + x2
        * (x19 * x271 + 2.0 * x274 + 2.0 * x275 + 4.0 * x333 + 4.0 * x334 + 4.0 * x340)
        - x23 * (-2.0 * ax * x429 + 3.0 * x269 + 3.0 * x270)
    )
    x431 = x112 * x13
    x432 = x430 * x431
    x433 = x431 * x5
    x434 = 3.14159265358979 * x110 * x13 * x3
    x435 = x434 * x5
    x436 = x120 * (x216 + x356)
    x437 = (
        x120 * x359
        + x2 * (x120 * x122 * x189 + x19 * x277)
        + x23 * (x0 * x436 - 3.0 * x277)
    )
    x438 = x126 * x62
    x439 = x2 * (x358 + x370)
    x440 = x120 * x369
    x441 = x2 * (x355 + x373)
    x442 = x120 * x363
    x443 = x441 + x442
    x444 = 2.0 * ax * x443 - x237 - 3.0 * x281
    x445 = x23 * x444 + x439 + x440
    x446 = x126 * x47
    x447 = x2 * (2.0 * x286 + 2.0 * x361 + 2.0 * x362 + x386)
    x448 = x120 * x376
    x449 = x447 + x448
    x450 = x126 * x449
    x451 = x2 * (
        x141 * x367 + x19 * x287 + 2.0 * x213 + 2.0 * x285 + 2.0 * x365 + 2.0 * x366
    )
    x452 = x120 * x378
    x453 = x23 * (2.0 * ax * x449 - 6.0 * x295 - x382)
    x454 = x451 + x452 + x453
    x455 = 3.0 * x374 + 3.0 * x375
    x456 = x2 * (2.0 * x290 + x396 + x455)
    x457 = x120 * x385
    x458 = x456 + x457
    x459 = x126 * x458
    x460 = 3.0 * x371 + 3.0 * x372 + 3.0 * x377
    x461 = x2 * (x19 * x291 + 2.0 * x236 + 2.0 * x293 + x460)
    x462 = x120 * x388
    x463 = x23 * (2.0 * ax * x458 - 3.0 * x238 - 3.0 * x290)
    x464 = x461 + x462 + x463
    x465 = 4.0 * x383 + 4.0 * x384
    x466 = x120 * x392 + x2 * (2.0 * x297 + 2.0 * x298 + x465)
    x467 = (
        x120 * x398
        + x2
        * (x19 * x299 + 2.0 * x302 + 2.0 * x303 + 4.0 * x380 + 4.0 * x381 + 4.0 * x387)
        - x23 * (-2.0 * ax * x466 + 3.0 * x297 + 3.0 * x298)
    )
    x468 = x434 * x467
    x469 = x115 * x399 + x2 * (x225 + 3.0 * x306)
    x470 = x134 * x469
    x471 = (
        x115 * x400
        + x2 * (x223 + 3.0 * x305 + 3.0 * x308)
        - x23 * (-2.0 * ax * x469 + x229 + 4.0 * x306)
    )
    x472 = x115 * x406 + x2 * (3.0 * x312 + 3.0 * x313 + x399)
    x473 = x134 * x472
    x474 = (
        x115 * x408
        + x2 * (3.0 * x316 + 3.0 * x317 + 3.0 * x319 + x400)
        - 2.0 * x23 * (-ax * x472 + 2.0 * x312 + 2.0 * x313)
    )
    x475 = x115 * x412 + x2 * (2.0 * x404 + 2.0 * x405 + x418)
    x476 = x134 * x475
    x477 = (
        x115 * x417
        + x2 * (x141 * x407 + 2.0 * x402 + 2.0 * x403 + x423)
        - 2.0 * x23 * (-ax * x475 + 2.0 * x325 + 2.0 * x326)
    )
    x478 = x115 * x421 + 3.0 * x2 * (x336 + x337 + x410 + x411)
    x479 = x431 * (
        x115 * x427
        + 3.0 * x2 * (x333 + x334 + x340 + x414 + x415 + x416)
        + x23 * (2.0 * ax * x478 - x428)
    )
    x480 = x10 * x431
    x481 = x115 * x429 + x2 * (3.0 * x343 + 3.0 * x344 + 4.0 * x419 + 4.0 * x420)
    x482 = x431 * x481
    x483 = x431 * x9
    x484 = x15 * x475
    x485 = x15 * x472
    x486 = x15 * x469
    x487 = x15 * x421
    x488 = x15 * x412
    x489 = x15 * x291
    x490 = x15 * x299
    x491 = x15 * x345
    x492 = x15 * x363
    x493 = x15 * x327
    x494 = x15 * x314
    x495 = x15 * x392
    x496 = x115 * x434
    x497 = x15 * x271
    x498 = x15 * x263
    x499 = x15 * x449
    x500 = x15 * x458
    x501 = x120 * x436 + x2 * (x237 + 3.0 * x355)
    x502 = x126 * x501
    x503 = (
        x120 * x437
        + x2 * (x235 + 3.0 * x354 + 3.0 * x357)
        - x23 * (-2.0 * ax * x501 + x241 + 4.0 * x355)
    )
    x504 = x120 * x443 + x2 * (3.0 * x361 + 3.0 * x362 + x436)
    x505 = x126 * x504
    x506 = (
        x120 * x445
        + x2 * (3.0 * x365 + 3.0 * x366 + 3.0 * x368 + x437)
        - 2.0 * x23 * (-ax * x504 + 2.0 * x361 + 2.0 * x362)
    )
    x507 = x120 * x449 + x2 * (2.0 * x441 + 2.0 * x442 + x455)
    x508 = x126 * x507
    x509 = (
        x120 * x454
        + x2 * (x141 * x444 + 2.0 * x439 + 2.0 * x440 + x460)
        - 2.0 * x23 * (-ax * x507 + 2.0 * x374 + 2.0 * x375)
    )
    x510 = x10 * x434
    x511 = x120 * x458 + 3.0 * x2 * (x383 + x384 + x447 + x448)
    x512 = x434 * (
        x120 * x464
        + 3.0 * x2 * (x380 + x381 + x387 + x451 + x452 + x453)
        + x23 * (2.0 * ax * x511 - x465)
    )
    x513 = x15 * x501
    x514 = x15 * x504
    x515 = x15 * x507
    x516 = x120 * x466 + x2 * (3.0 * x390 + 3.0 * x391 + 4.0 * x456 + 4.0 * x457)
    x517 = x434 * x516

    # 225 item(s)
    result[0, 0] = numpy.sum(
        x124
        * (
            x113
            * (
                x108 * x5
                + x2
                * (4.0 * x104 + 3.0 * x53 + 3.0 * x64 + 4.0 * x82 + 4.0 * x94 + 3.0 * x99)
                - 2.0 * x23 * (-ax * x109 + 2.0 * x95 + 2.0 * x96)
            )
            + x117 * x118
            + x118 * x122
        )
    )
    result[0, 1] = numpy.sum(
        x150 * (x122 * x125 * x149 + x125 * x148 + x128 * x133 * x134)
    )
    result[0, 2] = numpy.sum(
        x150 * (x117 * x149 * x151 + x126 * x133 * x153 + x148 * x151)
    )
    result[0, 3] = numpy.sum(
        x181 * (x122 * x158 * x168 + x134 * x158 * x180 + x163 * x168)
    )
    result[0, 4] = numpy.sum(
        x183 * (x125 * x180 * x182 + x127 * x153 * x167 + x128 * x152 * x167)
    )
    result[0, 5] = numpy.sum(
        x181 * (x117 * x188 * x194 + x126 * x180 * x188 + x193 * x194)
    )
    result[0, 6] = numpy.sum(
        x150 * (x122 * x201 * x207 + x134 * x201 * x211 + x204 * x207)
    )
    result[0, 7] = numpy.sum(
        x183 * (x152 * x158 * x211 + x152 * x163 * x206 + x153 * x158 * x206)
    )
    result[0, 8] = numpy.sum(
        x183 * (x127 * x188 * x211 + x127 * x193 * x206 + x128 * x188 * x206)
    )
    result[0, 9] = numpy.sum(
        x150 * (x117 * x218 * x222 + x126 * x211 * x218 + x221 * x222)
    )
    result[0, 10] = numpy.sum(
        x124 * (x122 * x228 * x233 + x134 * x228 * x234 + x231 * x233)
    )
    result[0, 11] = numpy.sum(
        x150 * (x152 * x201 * x234 + x152 * x204 * x232 + x153 * x201 * x232)
    )
    result[0, 12] = numpy.sum(
        x181 * (x158 * x188 * x234 + x158 * x193 * x232 + x163 * x188 * x232)
    )
    result[0, 13] = numpy.sum(
        x150 * (x127 * x218 * x234 + x127 * x221 * x232 + x128 * x218 * x232)
    )
    result[0, 14] = numpy.sum(
        x124 * (x117 * x240 * x244 + x126 * x234 * x240 + x243 * x244)
    )
    result[1, 0] = numpy.sum(
        x150 * (x107 * x134 * x246 + x115 * x122 * x248 + x115 * x247)
    )
    result[1, 1] = numpy.sum(
        x254 * (x122 * x251 * x253 + x134 * x147 * x251 + x252 * x253)
    )
    result[1, 2] = numpy.sum(
        x254 * (x103 * x152 * x246 + x103 * x153 * x245 + x115 * x147 * x182)
    )
    result[1, 3] = numpy.sum(
        x260 * (x122 * x257 * x259 + x134 * x179 * x257 + x258 * x259)
    )
    result[1, 4] = numpy.sum(
        x261 * (x145 * x152 * x252 + x145 * x153 * x251 + x152 * x179 * x251)
    )
    result[1, 5] = numpy.sum(
        x260 * (x145 * x188 * x246 + x145 * x193 * x245 + x179 * x188 * x245)
    )
    result[1, 6] = numpy.sum(
        x254 * (x122 * x177 * x264 + x134 * x177 * x266 + x210 * x264)
    )
    result[1, 7] = numpy.sum(
        x261 * (x152 * x177 * x258 + x152 * x210 * x257 + x153 * x177 * x257)
    )
    result[1, 8] = numpy.sum(
        x261 * (x177 * x188 * x252 + x177 * x193 * x251 + x188 * x210 * x251)
    )
    result[1, 9] = numpy.sum(
        x254 * (x177 * x218 * x246 + x177 * x221 * x245 + x210 * x218 * x245)
    )
    result[1, 10] = numpy.sum(
        x150 * (x122 * x205 * x272 + x134 * x205 * x276 + x209 * x272)
    )
    result[1, 11] = numpy.sum(
        x254 * (x152 * x205 * x266 + x152 * x209 * x263 + x153 * x205 * x263)
    )
    result[1, 12] = numpy.sum(
        x260 * (x188 * x205 * x258 + x188 * x209 * x257 + x193 * x205 * x257)
    )
    result[1, 13] = numpy.sum(
        x254 * (x205 * x218 * x252 + x205 * x221 * x251 + x209 * x218 * x251)
    )
    result[1, 14] = numpy.sum(
        x150 * (x205 * x240 * x246 + x205 * x243 * x245 + x209 * x240 * x245)
    )
    result[2, 0] = numpy.sum(
        x150 * (x107 * x126 * x278 + x117 * x120 * x248 + x120 * x247)
    )
    result[2, 1] = numpy.sum(
        x254 * (x103 * x127 * x278 + x103 * x128 * x277 + x125 * x147 * x279)
    )
    result[2, 2] = numpy.sum(
        x254 * (x117 * x282 * x284 + x126 * x147 * x282 + x283 * x284)
    )
    result[2, 3] = numpy.sum(
        x260 * (x145 * x158 * x278 + x145 * x163 * x277 + x158 * x179 * x277)
    )
    result[2, 4] = numpy.sum(
        x261 * (x127 * x145 * x283 + x127 * x179 * x282 + x128 * x145 * x282)
    )
    result[2, 5] = numpy.sum(
        x260 * (x117 * x287 * x289 + x126 * x179 * x287 + x288 * x289)
    )
    result[2, 6] = numpy.sum(
        x254 * (x177 * x201 * x278 + x177 * x204 * x277 + x201 * x210 * x277)
    )
    result[2, 7] = numpy.sum(
        x261 * (x158 * x177 * x283 + x158 * x210 * x282 + x163 * x177 * x282)
    )
    result[2, 8] = numpy.sum(
        x261 * (x127 * x177 * x288 + x127 * x210 * x287 + x128 * x177 * x287)
    )
    result[2, 9] = numpy.sum(
        x254 * (x117 * x177 * x292 + x126 * x177 * x294 + x210 * x292)
    )
    result[2, 10] = numpy.sum(
        x150 * (x205 * x228 * x278 + x205 * x231 * x277 + x209 * x228 * x277)
    )
    result[2, 11] = numpy.sum(
        x254 * (x201 * x205 * x283 + x201 * x209 * x282 + x204 * x205 * x282)
    )
    result[2, 12] = numpy.sum(
        x260 * (x158 * x205 * x288 + x158 * x209 * x287 + x163 * x205 * x287)
    )
    result[2, 13] = numpy.sum(
        x254 * (x127 * x205 * x294 + x127 * x209 * x291 + x128 * x205 * x291)
    )
    result[2, 14] = numpy.sum(
        x150 * (x117 * x205 * x300 + x126 * x205 * x304 + x209 * x300)
    )
    result[3, 0] = numpy.sum(
        x181 * (x105 * x134 * x307 + x122 * x307 * x311 + x310 * x311)
    )
    result[3, 1] = numpy.sum(x260 * (x122 * x315 * x90 + x134 * x320 * x90 + x315 * x93))
    result[3, 2] = numpy.sum(
        x260 * (x152 * x307 * x93 + x152 * x310 * x90 + x153 * x307 * x90)
    )
    result[3, 3] = numpy.sum(
        x331 * (x122 * x327 * x330 + x134 * x143 * x327 + x329 * x330)
    )
    result[3, 4] = numpy.sum(
        x332 * (x143 * x152 * x314 + x152 * x320 * x79 + x153 * x314 * x79)
    )
    result[3, 5] = numpy.sum(
        x331 * (x143 * x188 * x307 + x188 * x310 * x79 + x193 * x307 * x79)
    )
    result[3, 6] = numpy.sum(
        x260 * (x122 * x338 * x342 + x134 * x175 * x338 + x341 * x342)
    )
    result[3, 7] = numpy.sum(
        x332 * (x139 * x152 * x329 + x139 * x153 * x327 + x152 * x175 * x327)
    )
    result[3, 8] = numpy.sum(
        x332 * (x139 * x188 * x320 + x139 * x193 * x314 + x175 * x188 * x314)
    )
    result[3, 9] = numpy.sum(
        x260 * (x139 * x218 * x310 + x139 * x221 * x307 + x175 * x218 * x307)
    )
    result[3, 10] = numpy.sum(
        x181 * (x122 * x170 * x346 + x134 * x170 * x352 + x208 * x346)
    )
    result[3, 11] = numpy.sum(
        x260 * (x152 * x170 * x341 + x152 * x208 * x338 + x153 * x170 * x338)
    )
    result[3, 12] = numpy.sum(
        x331 * (x170 * x188 * x329 + x170 * x193 * x327 + x188 * x208 * x327)
    )
    result[3, 13] = numpy.sum(
        x260 * (x170 * x218 * x320 + x170 * x221 * x314 + x208 * x218 * x314)
    )
    result[3, 14] = numpy.sum(
        x181 * (x170 * x240 * x310 + x170 * x243 * x307 + x208 * x240 * x307)
    )
    result[4, 0] = numpy.sum(
        x183 * (x105 * x115 * x279 + x245 * x278 * x97 + x246 * x277 * x97)
    )
    result[4, 1] = numpy.sum(
        x261 * (x251 * x277 * x93 + x251 * x278 * x90 + x252 * x277 * x90)
    )
    result[4, 2] = numpy.sum(
        x261 * (x245 * x282 * x93 + x245 * x283 * x90 + x246 * x282 * x90)
    )
    result[4, 3] = numpy.sum(
        x332 * (x143 * x257 * x277 + x257 * x278 * x79 + x258 * x277 * x79)
    )
    result[4, 4] = numpy.sum(
        x353 * (x143 * x251 * x282 + x251 * x283 * x79 + x252 * x282 * x79)
    )
    result[4, 5] = numpy.sum(
        x332 * (x143 * x245 * x287 + x245 * x288 * x79 + x246 * x287 * x79)
    )
    result[4, 6] = numpy.sum(
        x261 * (x139 * x263 * x278 + x139 * x266 * x277 + x175 * x263 * x277)
    )
    result[4, 7] = numpy.sum(
        x353 * (x139 * x257 * x283 + x139 * x258 * x282 + x175 * x257 * x282)
    )
    result[4, 8] = numpy.sum(
        x353 * (x139 * x251 * x288 + x139 * x252 * x287 + x175 * x251 * x287)
    )
    result[4, 9] = numpy.sum(
        x261 * (x139 * x245 * x294 + x139 * x246 * x291 + x175 * x245 * x291)
    )
    result[4, 10] = numpy.sum(
        x183 * (x170 * x271 * x278 + x170 * x276 * x277 + x208 * x271 * x277)
    )
    result[4, 11] = numpy.sum(
        x261 * (x170 * x263 * x283 + x170 * x266 * x282 + x208 * x263 * x282)
    )
    result[4, 12] = numpy.sum(
        x332 * (x170 * x257 * x288 + x170 * x258 * x287 + x208 * x257 * x287)
    )
    result[4, 13] = numpy.sum(
        x261 * (x170 * x251 * x294 + x170 * x252 * x291 + x208 * x251 * x291)
    )
    result[4, 14] = numpy.sum(
        x183 * (x170 * x245 * x304 + x170 * x246 * x299 + x208 * x245 * x299)
    )
    result[5, 0] = numpy.sum(
        x181 * (x105 * x126 * x356 + x117 * x356 * x360 + x359 * x360)
    )
    result[5, 1] = numpy.sum(
        x260 * (x127 * x356 * x93 + x127 * x359 * x90 + x128 * x356 * x90)
    )
    result[5, 2] = numpy.sum(x260 * (x117 * x364 * x90 + x126 * x369 * x90 + x364 * x93))
    result[5, 3] = numpy.sum(
        x331 * (x143 * x158 * x356 + x158 * x359 * x79 + x163 * x356 * x79)
    )
    result[5, 4] = numpy.sum(
        x332 * (x127 * x143 * x363 + x127 * x369 * x79 + x128 * x363 * x79)
    )
    result[5, 5] = numpy.sum(
        x331 * (x117 * x376 * x379 + x126 * x143 * x376 + x378 * x379)
    )
    result[5, 6] = numpy.sum(
        x260 * (x139 * x201 * x359 + x139 * x204 * x356 + x175 * x201 * x356)
    )
    result[5, 7] = numpy.sum(
        x332 * (x139 * x158 * x369 + x139 * x163 * x363 + x158 * x175 * x363)
    )
    result[5, 8] = numpy.sum(
        x332 * (x127 * x139 * x378 + x127 * x175 * x376 + x128 * x139 * x376)
    )
    result[5, 9] = numpy.sum(
        x260 * (x117 * x385 * x389 + x126 * x175 * x385 + x388 * x389)
    )
    result[5, 10] = numpy.sum(
        x181 * (x170 * x228 * x359 + x170 * x231 * x356 + x208 * x228 * x356)
    )
    result[5, 11] = numpy.sum(
        x260 * (x170 * x201 * x369 + x170 * x204 * x363 + x201 * x208 * x363)
    )
    result[5, 12] = numpy.sum(
        x331 * (x158 * x170 * x378 + x158 * x208 * x376 + x163 * x170 * x376)
    )
    result[5, 13] = numpy.sum(
        x260 * (x127 * x170 * x388 + x127 * x208 * x385 + x128 * x170 * x385)
    )
    result[5, 14] = numpy.sum(
        x181 * (x117 * x170 * x393 + x126 * x170 * x398 + x208 * x393)
    )
    result[6, 0] = numpy.sum(
        x150 * (x122 * x399 * x401 + x134 * x399 * x63 + x400 * x401)
    )
    result[6, 1] = numpy.sum(
        x254 * (x122 * x406 * x409 + x134 * x406 * x85 + x408 * x409)
    )
    result[6, 2] = numpy.sum(
        x254 * (x152 * x399 * x85 + x152 * x400 * x47 + x153 * x399 * x47)
    )
    result[6, 3] = numpy.sum(x260 * (x122 * x413 * x73 + x134 * x417 * x73 + x413 * x74))
    result[6, 4] = numpy.sum(
        x261 * (x152 * x406 * x74 + x152 * x408 * x73 + x153 * x406 * x73)
    )
    result[6, 5] = numpy.sum(
        x260 * (x188 * x399 * x74 + x188 * x400 * x73 + x193 * x399 * x73)
    )
    result[6, 6] = numpy.sum(x254 * (x122 * x422 * x68 + x134 * x427 * x68 + x137 * x422))
    result[6, 7] = numpy.sum(
        x261 * (x137 * x152 * x412 + x152 * x417 * x68 + x153 * x412 * x68)
    )
    result[6, 8] = numpy.sum(
        x261 * (x137 * x188 * x406 + x188 * x408 * x68 + x193 * x406 * x68)
    )
    result[6, 9] = numpy.sum(
        x254 * (x137 * x218 * x399 + x218 * x400 * x68 + x221 * x399 * x68)
    )
    result[6, 10] = numpy.sum(
        x150 * (x122 * x429 * x433 + x134 * x135 * x429 + x432 * x5)
    )
    result[6, 11] = numpy.sum(
        x254 * (x135 * x152 * x421 + x151 * x427 * x433 + x153 * x421 * x66)
    )
    result[6, 12] = numpy.sum(
        x260 * (x135 * x188 * x412 + x188 * x417 * x66 + x193 * x412 * x66)
    )
    result[6, 13] = numpy.sum(
        x254 * (x135 * x218 * x406 + x218 * x408 * x66 + x221 * x406 * x66)
    )
    result[6, 14] = numpy.sum(
        x150 * (x135 * x240 * x399 + x240 * x400 * x66 + x243 * x399 * x66)
    )
    result[7, 0] = numpy.sum(
        x183 * (x277 * x307 * x63 + x277 * x310 * x62 + x278 * x307 * x62)
    )
    result[7, 1] = numpy.sum(
        x261 * (x277 * x314 * x85 + x277 * x320 * x47 + x278 * x314 * x47)
    )
    result[7, 2] = numpy.sum(
        x261 * (x282 * x307 * x85 + x282 * x310 * x47 + x283 * x307 * x47)
    )
    result[7, 3] = numpy.sum(
        x332 * (x277 * x327 * x74 + x277 * x329 * x73 + x278 * x327 * x73)
    )
    result[7, 4] = numpy.sum(
        x353 * (x282 * x314 * x74 + x282 * x320 * x73 + x283 * x314 * x73)
    )
    result[7, 5] = numpy.sum(
        x332 * (x287 * x307 * x74 + x287 * x310 * x73 + x288 * x307 * x73)
    )
    result[7, 6] = numpy.sum(
        x261 * (x137 * x277 * x338 + x277 * x341 * x68 + x278 * x338 * x68)
    )
    result[7, 7] = numpy.sum(
        x353 * (x137 * x282 * x327 + x282 * x329 * x68 + x283 * x327 * x68)
    )
    result[7, 8] = numpy.sum(
        x353 * (x137 * x287 * x314 + x287 * x320 * x68 + x288 * x314 * x68)
    )
    result[7, 9] = numpy.sum(
        x261 * (x137 * x291 * x307 + x291 * x310 * x68 + x294 * x307 * x68)
    )
    result[7, 10] = numpy.sum(
        x183 * (x120 * x352 * x433 + x135 * x277 * x345 + x278 * x345 * x66)
    )
    result[7, 11] = numpy.sum(
        x261 * (x135 * x282 * x338 + x282 * x341 * x66 + x283 * x338 * x66)
    )
    result[7, 12] = numpy.sum(
        x332 * (x135 * x287 * x327 + x287 * x329 * x66 + x288 * x327 * x66)
    )
    result[7, 13] = numpy.sum(
        x261 * (x135 * x291 * x314 + x291 * x320 * x66 + x294 * x314 * x66)
    )
    result[7, 14] = numpy.sum(
        x183 * (x135 * x299 * x307 + x299 * x310 * x66 + x304 * x307 * x66)
    )
    result[8, 0] = numpy.sum(
        x183 * (x245 * x356 * x63 + x245 * x359 * x62 + x246 * x356 * x62)
    )
    result[8, 1] = numpy.sum(
        x261 * (x251 * x356 * x85 + x251 * x359 * x47 + x252 * x356 * x47)
    )
    result[8, 2] = numpy.sum(
        x261 * (x245 * x363 * x85 + x245 * x369 * x47 + x246 * x363 * x47)
    )
    result[8, 3] = numpy.sum(
        x332 * (x257 * x356 * x74 + x257 * x359 * x73 + x258 * x356 * x73)
    )
    result[8, 4] = numpy.sum(
        x353 * (x251 * x363 * x74 + x251 * x369 * x73 + x252 * x363 * x73)
    )
    result[8, 5] = numpy.sum(
        x332 * (x245 * x376 * x74 + x245 * x378 * x73 + x246 * x376 * x73)
    )
    result[8, 6] = numpy.sum(
        x261 * (x137 * x263 * x356 + x263 * x359 * x68 + x266 * x356 * x68)
    )
    result[8, 7] = numpy.sum(
        x353 * (x137 * x257 * x363 + x257 * x369 * x68 + x258 * x363 * x68)
    )
    result[8, 8] = numpy.sum(
        x353 * (x137 * x251 * x376 + x251 * x378 * x68 + x252 * x376 * x68)
    )
    result[8, 9] = numpy.sum(
        x261 * (x137 * x245 * x385 + x245 * x388 * x68 + x246 * x385 * x68)
    )
    result[8, 10] = numpy.sum(
        x183 * (x135 * x271 * x356 + x271 * x359 * x66 + x276 * x356 * x66)
    )
    result[8, 11] = numpy.sum(
        x261 * (x135 * x263 * x363 + x263 * x369 * x66 + x266 * x363 * x66)
    )
    result[8, 12] = numpy.sum(
        x332 * (x135 * x257 * x376 + x257 * x378 * x66 + x258 * x376 * x66)
    )
    result[8, 13] = numpy.sum(
        x261 * (x135 * x251 * x385 + x251 * x388 * x66 + x252 * x385 * x66)
    )
    result[8, 14] = numpy.sum(
        x183 * (x115 * x398 * x435 + x135 * x245 * x392 + x246 * x392 * x66)
    )
    result[9, 0] = numpy.sum(
        x150 * (x117 * x436 * x438 + x126 * x436 * x63 + x437 * x438)
    )
    result[9, 1] = numpy.sum(
        x254 * (x127 * x436 * x85 + x127 * x437 * x47 + x128 * x436 * x47)
    )
    result[9, 2] = numpy.sum(
        x254 * (x117 * x443 * x446 + x126 * x443 * x85 + x445 * x446)
    )
    result[9, 3] = numpy.sum(
        x260 * (x158 * x436 * x74 + x158 * x437 * x73 + x163 * x436 * x73)
    )
    result[9, 4] = numpy.sum(
        x261 * (x127 * x443 * x74 + x127 * x445 * x73 + x128 * x443 * x73)
    )
    result[9, 5] = numpy.sum(x260 * (x117 * x450 * x73 + x126 * x454 * x73 + x450 * x74))
    result[9, 6] = numpy.sum(
        x254 * (x137 * x201 * x436 + x201 * x437 * x68 + x204 * x436 * x68)
    )
    result[9, 7] = numpy.sum(
        x261 * (x137 * x158 * x443 + x158 * x445 * x68 + x163 * x443 * x68)
    )
    result[9, 8] = numpy.sum(
        x261 * (x127 * x137 * x449 + x127 * x454 * x68 + x128 * x449 * x68)
    )
    result[9, 9] = numpy.sum(x254 * (x117 * x459 * x68 + x126 * x464 * x68 + x137 * x459))
    result[9, 10] = numpy.sum(
        x150 * (x135 * x228 * x436 + x228 * x437 * x66 + x231 * x436 * x66)
    )
    result[9, 11] = numpy.sum(
        x254 * (x135 * x201 * x443 + x201 * x445 * x66 + x204 * x443 * x66)
    )
    result[9, 12] = numpy.sum(
        x260 * (x135 * x158 * x449 + x158 * x454 * x66 + x163 * x449 * x66)
    )
    result[9, 13] = numpy.sum(
        x254 * (x125 * x435 * x464 + x127 * x135 * x458 + x128 * x458 * x66)
    )
    result[9, 14] = numpy.sum(
        x150 * (x117 * x435 * x466 + x126 * x135 * x466 + x468 * x5)
    )
    result[10, 0] = numpy.sum(x124 * (x122 * x470 * x50 + x134 * x471 * x50 + x470 * x56))
    result[10, 1] = numpy.sum(x150 * (x122 * x37 * x473 + x134 * x37 * x474 + x40 * x473))
    result[10, 2] = numpy.sum(
        x150 * (x152 * x37 * x471 + x152 * x40 * x469 + x153 * x37 * x469)
    )
    result[10, 3] = numpy.sum(x181 * (x122 * x28 * x476 + x134 * x28 * x477 + x32 * x476))
    result[10, 4] = numpy.sum(
        x183 * (x152 * x28 * x474 + x152 * x32 * x472 + x153 * x28 * x472)
    )
    result[10, 5] = numpy.sum(
        x181 * (x188 * x28 * x471 + x188 * x32 * x469 + x193 * x28 * x469)
    )
    result[10, 6] = numpy.sum(
        x150 * (x10 * x479 + x122 * x478 * x480 + x134 * x25 * x478)
    )
    result[10, 7] = numpy.sum(
        x183 * (x151 * x477 * x480 + x152 * x25 * x475 + x153 * x18 * x475)
    )
    result[10, 8] = numpy.sum(
        x183 * (x18 * x188 * x474 + x18 * x193 * x472 + x188 * x25 * x472)
    )
    result[10, 9] = numpy.sum(
        x150 * (x18 * x218 * x471 + x18 * x221 * x469 + x218 * x25 * x469)
    )
    result[10, 10] = numpy.sum(
        x124
        * (
            x122 * x482
            + x431
            * (
                x115 * x430
                + x2
                * (
                    3.0 * x348
                    + 3.0 * x349
                    + 3.0 * x351
                    + 4.0 * x424
                    + 4.0 * x425
                    + 4.0 * x426
                )
                - 2.0 * x23 * (-ax * x481 + 2.0 * x343 + 2.0 * x344)
            )
            + x482 * x9
        )
    )
    result[10, 11] = numpy.sum(
        x150 * (x15 * x153 * x478 + x151 * x478 * x483 + x151 * x479)
    )
    result[10, 12] = numpy.sum(
        x181 * (x15 * x188 * x477 + x188 * x484 * x9 + x193 * x484)
    )
    result[10, 13] = numpy.sum(
        x150 * (x15 * x218 * x474 + x218 * x485 * x9 + x221 * x485)
    )
    result[10, 14] = numpy.sum(
        x124 * (x15 * x240 * x471 + x240 * x486 * x9 + x243 * x486)
    )
    result[11, 0] = numpy.sum(
        x150 * (x277 * x399 * x56 + x277 * x400 * x50 + x278 * x399 * x50)
    )
    result[11, 1] = numpy.sum(
        x254 * (x277 * x37 * x408 + x277 * x40 * x406 + x278 * x37 * x406)
    )
    result[11, 2] = numpy.sum(
        x254 * (x282 * x37 * x400 + x282 * x399 * x40 + x283 * x37 * x399)
    )
    result[11, 3] = numpy.sum(
        x260 * (x277 * x28 * x417 + x277 * x32 * x412 + x278 * x28 * x412)
    )
    result[11, 4] = numpy.sum(
        x261 * (x28 * x282 * x408 + x28 * x283 * x406 + x282 * x32 * x406)
    )
    result[11, 5] = numpy.sum(
        x260 * (x28 * x287 * x400 + x28 * x288 * x399 + x287 * x32 * x399)
    )
    result[11, 6] = numpy.sum(
        x254 * (x120 * x427 * x480 + x18 * x278 * x421 + x25 * x277 * x421)
    )
    result[11, 7] = numpy.sum(
        x261 * (x18 * x282 * x417 + x18 * x283 * x412 + x25 * x282 * x412)
    )
    result[11, 8] = numpy.sum(
        x261 * (x18 * x287 * x408 + x18 * x288 * x406 + x25 * x287 * x406)
    )
    result[11, 9] = numpy.sum(
        x254 * (x18 * x291 * x400 + x18 * x294 * x399 + x25 * x291 * x399)
    )
    result[11, 10] = numpy.sum(
        x150 * (x120 * x429 * x483 + x120 * x432 + x15 * x278 * x429)
    )
    result[11, 11] = numpy.sum(
        x254 * (x15 * x282 * x427 + x282 * x487 * x9 + x283 * x487)
    )
    result[11, 12] = numpy.sum(
        x260 * (x15 * x287 * x417 + x287 * x488 * x9 + x288 * x488)
    )
    result[11, 13] = numpy.sum(
        x254 * (x15 * x294 * x406 + x406 * x489 * x9 + x408 * x489)
    )
    result[11, 14] = numpy.sum(
        x150 * (x15 * x304 * x399 + x399 * x490 * x9 + x400 * x490)
    )
    result[12, 0] = numpy.sum(
        x181 * (x307 * x356 * x56 + x307 * x359 * x50 + x310 * x356 * x50)
    )
    result[12, 1] = numpy.sum(
        x260 * (x314 * x356 * x40 + x314 * x359 * x37 + x320 * x356 * x37)
    )
    result[12, 2] = numpy.sum(
        x260 * (x307 * x363 * x40 + x307 * x369 * x37 + x310 * x363 * x37)
    )
    result[12, 3] = numpy.sum(
        x331 * (x28 * x327 * x359 + x28 * x329 * x356 + x32 * x327 * x356)
    )
    result[12, 4] = numpy.sum(
        x332 * (x28 * x314 * x369 + x28 * x320 * x363 + x314 * x32 * x363)
    )
    result[12, 5] = numpy.sum(
        x331 * (x28 * x307 * x378 + x28 * x310 * x376 + x307 * x32 * x376)
    )
    result[12, 6] = numpy.sum(
        x260 * (x18 * x338 * x359 + x18 * x341 * x356 + x25 * x338 * x356)
    )
    result[12, 7] = numpy.sum(
        x332 * (x18 * x327 * x369 + x18 * x329 * x363 + x25 * x327 * x363)
    )
    result[12, 8] = numpy.sum(
        x332 * (x18 * x314 * x378 + x18 * x320 * x376 + x25 * x314 * x376)
    )
    result[12, 9] = numpy.sum(
        x260 * (x18 * x307 * x388 + x18 * x310 * x385 + x25 * x307 * x385)
    )
    result[12, 10] = numpy.sum(
        x181 * (x15 * x352 * x356 + x356 * x491 * x9 + x359 * x491)
    )
    result[12, 11] = numpy.sum(
        x260 * (x15 * x338 * x369 + x338 * x492 * x9 + x341 * x492)
    )
    result[12, 12] = numpy.sum(
        x331 * (x15 * x329 * x376 + x376 * x493 * x9 + x378 * x493)
    )
    result[12, 13] = numpy.sum(
        x260 * (x15 * x320 * x385 + x385 * x494 * x9 + x388 * x494)
    )
    result[12, 14] = numpy.sum(
        x181 * (x15 * x307 * x398 + x307 * x495 * x9 + x310 * x495)
    )
    result[13, 0] = numpy.sum(
        x150 * (x245 * x436 * x56 + x245 * x437 * x50 + x246 * x436 * x50)
    )
    result[13, 1] = numpy.sum(
        x254 * (x251 * x37 * x437 + x251 * x40 * x436 + x252 * x37 * x436)
    )
    result[13, 2] = numpy.sum(
        x254 * (x245 * x37 * x445 + x245 * x40 * x443 + x246 * x37 * x443)
    )
    result[13, 3] = numpy.sum(
        x260 * (x257 * x28 * x437 + x257 * x32 * x436 + x258 * x28 * x436)
    )
    result[13, 4] = numpy.sum(
        x261 * (x251 * x28 * x445 + x251 * x32 * x443 + x252 * x28 * x443)
    )
    result[13, 5] = numpy.sum(
        x260 * (x245 * x28 * x454 + x245 * x32 * x449 + x246 * x28 * x449)
    )
    result[13, 6] = numpy.sum(
        x254 * (x18 * x263 * x437 + x18 * x266 * x436 + x25 * x263 * x436)
    )
    result[13, 7] = numpy.sum(
        x261 * (x18 * x257 * x445 + x18 * x258 * x443 + x25 * x257 * x443)
    )
    result[13, 8] = numpy.sum(
        x261 * (x18 * x251 * x454 + x18 * x252 * x449 + x25 * x251 * x449)
    )
    result[13, 9] = numpy.sum(
        x254 * (x10 * x464 * x496 + x18 * x246 * x458 + x245 * x25 * x458)
    )
    result[13, 10] = numpy.sum(
        x150 * (x15 * x276 * x436 + x436 * x497 * x9 + x437 * x497)
    )
    result[13, 11] = numpy.sum(
        x254 * (x15 * x266 * x443 + x443 * x498 * x9 + x445 * x498)
    )
    result[13, 12] = numpy.sum(
        x260 * (x15 * x257 * x454 + x257 * x499 * x9 + x258 * x499)
    )
    result[13, 13] = numpy.sum(
        x254 * (x15 * x251 * x464 + x251 * x500 * x9 + x252 * x500)
    )
    result[13, 14] = numpy.sum(
        x150 * (x115 * x468 + x15 * x246 * x466 + x466 * x496 * x9)
    )
    result[14, 0] = numpy.sum(x124 * (x117 * x50 * x502 + x126 * x50 * x503 + x502 * x56))
    result[14, 1] = numpy.sum(
        x150 * (x127 * x37 * x503 + x127 * x40 * x501 + x128 * x37 * x501)
    )
    result[14, 2] = numpy.sum(x150 * (x117 * x37 * x505 + x126 * x37 * x506 + x40 * x505))
    result[14, 3] = numpy.sum(
        x181 * (x158 * x28 * x503 + x158 * x32 * x501 + x163 * x28 * x501)
    )
    result[14, 4] = numpy.sum(
        x183 * (x127 * x28 * x506 + x127 * x32 * x504 + x128 * x28 * x504)
    )
    result[14, 5] = numpy.sum(x181 * (x117 * x28 * x508 + x126 * x28 * x509 + x32 * x508))
    result[14, 6] = numpy.sum(
        x150 * (x18 * x201 * x503 + x18 * x204 * x501 + x201 * x25 * x501)
    )
    result[14, 7] = numpy.sum(
        x183 * (x158 * x18 * x506 + x158 * x25 * x504 + x163 * x18 * x504)
    )
    result[14, 8] = numpy.sum(
        x183 * (x125 * x509 * x510 + x127 * x25 * x507 + x128 * x18 * x507)
    )
    result[14, 9] = numpy.sum(
        x150 * (x10 * x512 + x117 * x510 * x511 + x126 * x25 * x511)
    )
    result[14, 10] = numpy.sum(
        x124 * (x15 * x228 * x503 + x228 * x513 * x9 + x231 * x513)
    )
    result[14, 11] = numpy.sum(
        x150 * (x15 * x201 * x506 + x201 * x514 * x9 + x204 * x514)
    )
    result[14, 12] = numpy.sum(
        x181 * (x15 * x158 * x509 + x158 * x515 * x9 + x163 * x515)
    )
    result[14, 13] = numpy.sum(
        x150 * (x125 * x434 * x511 * x9 + x125 * x512 + x128 * x15 * x511)
    )
    result[14, 14] = numpy.sum(
        x124
        * (
            x117 * x517
            + x434
            * (
                x120 * x467
                + x2
                * (
                    3.0 * x394
                    + 3.0 * x395
                    + 3.0 * x397
                    + 4.0 * x461
                    + 4.0 * x462
                    + 4.0 * x463
                )
                - 2.0 * x23 * (-ax * x516 + 2.0 * x390 + 2.0 * x391)
            )
            + x517 * x9
        )
    )
    return result


kinetic3d = {
    (0, 0): kinetic3d_00,
    (0, 1): kinetic3d_01,
    (0, 2): kinetic3d_02,
    (0, 3): kinetic3d_03,
    (0, 4): kinetic3d_04,
    (1, 0): kinetic3d_10,
    (1, 1): kinetic3d_11,
    (1, 2): kinetic3d_12,
    (1, 3): kinetic3d_13,
    (1, 4): kinetic3d_14,
    (2, 0): kinetic3d_20,
    (2, 1): kinetic3d_21,
    (2, 2): kinetic3d_22,
    (2, 3): kinetic3d_23,
    (2, 4): kinetic3d_24,
    (3, 0): kinetic3d_30,
    (3, 1): kinetic3d_31,
    (3, 2): kinetic3d_32,
    (3, 3): kinetic3d_33,
    (3, 4): kinetic3d_34,
    (4, 0): kinetic3d_40,
    (4, 1): kinetic3d_41,
    (4, 2): kinetic3d_42,
    (4, 3): kinetic3d_43,
    (4, 4): kinetic3d_44,
}
