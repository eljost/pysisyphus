import numpy


def kinetic3d_00(a, A, b, B):
    """Cartesian 3D (ss) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = (2 * a + 2 * b) ** (-1.0)
    x2 = (a + b) ** (-1.0)
    x3 = 2 * a**2
    x4 = a * b * x2
    x5 = (
        numpy.pi ** (3 / 2)
        * x2 ** (3 / 2)
        * numpy.exp(-x4 * (A[0] - B[0]) ** 2)
        * numpy.exp(-x4 * (A[1] - B[1]) ** 2)
        * numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    )

    # 1 item(s)
    S = numpy.array(
        [
            x5 * (-x0 - x3 * (x1 + (x2 * (a * A[0] + b * B[0]) - A[0]) ** 2))
            + x5 * (-x0 - x3 * (x1 + (x2 * (a * A[1] + b * B[1]) - A[1]) ** 2))
            + x5 * (-x0 - x3 * (x1 + (x2 * (a * A[2] + b * B[2]) - A[2]) ** 2))
        ]
    )
    return S


def kinetic3d_01(a, A, b, B):
    """Cartesian 3D (sp) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = b * x0
    x3 = a * x2
    x4 = numpy.exp(-x3 * (A[0] - B[0]) ** 2)
    x5 = x4 * (-x1 - B[0])
    x6 = numpy.sqrt(x0)
    x7 = numpy.pi ** (3 / 2) * x6
    x8 = x5 * x7
    x9 = numpy.exp(-x3 * (A[2] - B[2]) ** 2)
    x10 = numpy.exp(-x3 * (A[1] - B[1]) ** 2)
    x11 = -a
    x12 = 2 * a
    x13 = (2 * b + x12) ** (-1.0)
    x14 = -x0 * (a * A[1] + b * B[1])
    x15 = 2 * a**2
    x16 = x10 * (-x11 - x15 * (x13 + (-x14 - A[1]) ** 2))
    x17 = x0 * x16 * x9
    x18 = -x0 * (a * A[2] + b * B[2])
    x19 = x9 * (-x11 - x15 * (x13 + (-x18 - A[2]) ** 2))
    x20 = x0 * x10
    x21 = numpy.sqrt(numpy.pi) * x6
    x22 = x21 * x5
    x23 = x12 * x2
    x24 = -x11 - x15 * (x13 + (-x1 - A[0]) ** 2)
    x25 = numpy.pi * x0 * x9
    x26 = -x14 - B[1]
    x27 = x20 * x4
    x28 = x26 * x27 * x7
    x29 = x24 * x9
    x30 = x21 * x26
    x31 = -x18 - B[2]
    x32 = x31 * x7
    x33 = x21 * x31

    # 3 item(s)
    S = numpy.array(
        [
            x10 * x25 * (x22 * x23 + x22 * x24) + x17 * x8 + x19 * x20 * x8,
            x19 * x28 + x25 * x4 * (x10 * x23 * x30 + x16 * x30) + x28 * x29,
            x17 * x32 * x4
            + x27 * x29 * x32
            + numpy.pi * x27 * (x19 * x33 + x23 * x33 * x9),
        ]
    )
    return S


def kinetic3d_02(a, A, b, B):
    """Cartesian 3D (sd) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = 2 * a**2
    x7 = -x0 - x6 * (x3 + (-x5 - A[1]) ** 2)
    x8 = a * x4
    x9 = b * x8
    x10 = numpy.exp(-x9 * (A[0] - B[0]) ** 2)
    x11 = numpy.sqrt(x4)
    x12 = numpy.sqrt(numpy.pi) * x11
    x13 = x10 * x12
    x14 = x13 * x3
    x15 = -x4 * (a * A[0] + b * B[0])
    x16 = -x15 - B[0]
    x17 = x13 * x16**2 + x14
    x18 = numpy.exp(-x9 * (A[1] - B[1]) ** 2)
    x19 = numpy.exp(-x9 * (A[2] - B[2]) ** 2)
    x20 = numpy.pi * x19 * x4
    x21 = x18 * x20
    x22 = x17 * x21
    x23 = -x4 * (a * A[2] + b * B[2])
    x24 = -x0 - x6 * (x3 + (-x23 - A[2]) ** 2)
    x25 = -x0 - x6 * (x3 + (-x15 - A[0]) ** 2)
    x26 = x13 * x16
    x27 = b * x1 * x4
    x28 = x25 * x26 + x26 * x27
    x29 = -x5 - B[1]
    x30 = numpy.pi ** (3 / 2)
    x31 = x10 * x18 * x4
    x32 = x11 * x16 * x19 * x30 * x31
    x33 = x12 * x18
    x34 = x29 * x33
    x35 = x27 * x34 + x34 * x7
    x36 = x10 * x20
    x37 = x35 * x36
    x38 = x21 * x28
    x39 = -x23 - B[2]
    x40 = x12 * x19
    x41 = x39 * x40
    x42 = x24 * x41 + x27 * x41
    x43 = numpy.pi * x31
    x44 = x42 * x43
    x45 = x3 * x33
    x46 = x29**2 * x33 + x45
    x47 = x36 * x46
    x48 = x3 * x40
    x49 = x39**2 * x40 + x48
    x50 = x43 * x49

    # 6 item(s)
    S = numpy.array(
        [
            x21 * (x14 * x25 + x16 * x28 + x8 * (-2 * x13 + x17 * x2))
            + x22 * x24
            + x22 * x7,
            x16 * x37 + x24 * x29 * x32 + x29 * x38,
            x16 * x44 + x32 * x39 * x7 + x38 * x39,
            x24 * x47
            + x25 * x47
            + x36 * (x29 * x35 + x45 * x7 + x8 * (x2 * x46 - 2 * x33)),
            x11 * x19 * x25 * x29 * x30 * x31 * x39 + x29 * x44 + x37 * x39,
            x25 * x50
            + x43 * (x24 * x48 + x39 * x42 + x8 * (x2 * x49 - 2 * x40))
            + x50 * x7,
        ]
    )
    return S


def kinetic3d_03(a, A, b, B):
    """Cartesian 3D (sf) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = 2 * a**2
    x7 = -x0 - x6 * (x3 + (-x5 - A[1]) ** 2)
    x8 = -x4 * (a * A[0] + b * B[0])
    x9 = -x8 - B[0]
    x10 = a * x4
    x11 = b * x10
    x12 = numpy.exp(-x11 * (A[0] - B[0]) ** 2)
    x13 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x14 = x12 * x13
    x15 = x14 * x3
    x16 = x14 * x9**2 + x15
    x17 = 2 * x15 * x9 + x16 * x9
    x18 = numpy.exp(-x11 * (A[1] - B[1]) ** 2)
    x19 = numpy.exp(-x11 * (A[2] - B[2]) ** 2)
    x20 = numpy.pi * x19 * x4
    x21 = x18 * x20
    x22 = x17 * x21
    x23 = -x4 * (a * A[2] + b * B[2])
    x24 = -x0 - x6 * (x3 + (-x23 - A[2]) ** 2)
    x25 = x14 * x9
    x26 = 4 * x11
    x27 = -x0 - x6 * (x3 + (-x8 - A[0]) ** 2)
    x28 = 2 * x14
    x29 = b * x1 * x4
    x30 = x25 * x27 + x25 * x29
    x31 = x10 * (x16 * x2 - x28) + x15 * x27 + x30 * x9
    x32 = -x5 - B[1]
    x33 = x16 * x21
    x34 = x13 * x18
    x35 = x32 * x34
    x36 = x29 * x35 + x35 * x7
    x37 = x13 * x19
    x38 = x21 * x31
    x39 = -x23 - B[2]
    x40 = x37 * x39
    x41 = x24 * x40 + x29 * x40
    x42 = x3 * x34
    x43 = x32**2 * x34 + x42
    x44 = x12 * x20
    x45 = x44 * x9
    x46 = 2 * x34
    x47 = x10 * (x2 * x43 - x46) + x32 * x36 + x42 * x7
    x48 = x44 * x47
    x49 = numpy.pi * x12 * x18 * x4
    x50 = x49 * x9
    x51 = x3 * x37
    x52 = x37 * x39**2 + x51
    x53 = 2 * x37
    x54 = x10 * (x2 * x52 - x53) + x24 * x51 + x39 * x41
    x55 = x49 * x54
    x56 = 2 * x32 * x42 + x32 * x43
    x57 = x44 * x56
    x58 = 2 * x39 * x51 + x39 * x52
    x59 = x49 * x58

    # 10 item(s)
    S = numpy.array(
        [
            x21
            * (x10 * (x17 * x2 - 3 * x25) + x3 * (x25 * x26 + x27 * x28 * x9) + x31 * x9)
            + x22 * x24
            + x22 * x7,
            x16 * x36 * x37 + x24 * x32 * x33 + x32 * x38,
            x16 * x34 * x41 + x33 * x39 * x7 + x38 * x39,
            x24 * x43 * x45 + x30 * x37 * x43 + x48 * x9,
            x21 * x30 * x32 * x39 + x32 * x41 * x50 + x36 * x39 * x45,
            x30 * x34 * x52 + x50 * x52 * x7 + x55 * x9,
            x24 * x57
            + x27 * x57
            + x44
            * (
                x10 * (x2 * x56 - 3 * x35) + x3 * (x26 * x35 + x32 * x46 * x7) + x32 * x47
            ),
            x14 * x41 * x43 + x27 * x39 * x43 * x44 + x39 * x48,
            x14 * x36 * x52 + x27 * x32 * x49 * x52 + x32 * x55,
            x27 * x59
            + x49
            * (
                x10 * (x2 * x58 - 3 * x40)
                + x3 * (x24 * x39 * x53 + x26 * x40)
                + x39 * x54
            )
            + x59 * x7,
        ]
    )
    return S


def kinetic3d_04(a, A, b, B):
    """Cartesian 3D (sg) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = 2 * a**2
    x7 = -x0 - x6 * (x3 + (-x5 - A[1]) ** 2)
    x8 = a * x4
    x9 = b * x8
    x10 = numpy.exp(-x9 * (A[0] - B[0]) ** 2)
    x11 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x12 = x10 * x11
    x13 = x12 * x3
    x14 = -x4 * (a * A[0] + b * B[0])
    x15 = -x14 - B[0]
    x16 = x12 * x15**2
    x17 = x13 + x16
    x18 = 2 * x13 * x15 + x15 * x17
    x19 = x15 * x18 + x3 * (3 * x13 + 3 * x16)
    x20 = numpy.exp(-x9 * (A[1] - B[1]) ** 2)
    x21 = numpy.exp(-x9 * (A[2] - B[2]) ** 2)
    x22 = numpy.pi * x21 * x4
    x23 = x20 * x22
    x24 = x19 * x23
    x25 = -x4 * (a * A[2] + b * B[2])
    x26 = -x0 - x6 * (x3 + (-x25 - A[2]) ** 2)
    x27 = -x0 - x6 * (x3 + (-x14 - A[0]) ** 2)
    x28 = x13 * x27
    x29 = 2 * x12
    x30 = x8 * (x17 * x2 - x29)
    x31 = x12 * x15
    x32 = b * x1 * x4
    x33 = x27 * x31 + x31 * x32
    x34 = x15 * x33
    x35 = 4 * x9
    x36 = x28 + x30 + x34
    x37 = x15 * x36 + x3 * (x15 * x27 * x29 + x31 * x35) + x8 * (x18 * x2 - 3 * x31)
    x38 = -x5 - B[1]
    x39 = x18 * x23
    x40 = x11 * x20
    x41 = x38 * x40
    x42 = x32 * x41 + x41 * x7
    x43 = x11 * x21
    x44 = x42 * x43
    x45 = x23 * x37
    x46 = -x25 - B[2]
    x47 = x43 * x46
    x48 = x26 * x47 + x32 * x47
    x49 = x3 * x40
    x50 = x38**2 * x40
    x51 = x49 + x50
    x52 = x17 * x43
    x53 = x49 * x7
    x54 = 2 * x40
    x55 = x8 * (x2 * x51 - x54)
    x56 = x38 * x42
    x57 = x53 + x55 + x56
    x58 = x3 * x43
    x59 = x43 * x46**2
    x60 = x58 + x59
    x61 = x17 * x40
    x62 = x26 * x58
    x63 = 2 * x43
    x64 = x8 * (x2 * x60 - x63)
    x65 = x46 * x48
    x66 = x62 + x64 + x65
    x67 = 2 * x38 * x49 + x38 * x51
    x68 = x10 * x22
    x69 = x15 * x68
    x70 = x3 * (x35 * x41 + x38 * x54 * x7) + x38 * x57 + x8 * (x2 * x67 - 3 * x41)
    x71 = x68 * x70
    x72 = numpy.pi * x10 * x20 * x4
    x73 = x15 * x72
    x74 = 2 * x46 * x58 + x46 * x60
    x75 = x3 * (x26 * x46 * x63 + x35 * x47) + x46 * x66 + x8 * (x2 * x74 - 3 * x47)
    x76 = x72 * x75
    x77 = x3 * (3 * x49 + 3 * x50) + x38 * x67
    x78 = x68 * x77
    x79 = x12 * x51
    x80 = x3 * (3 * x58 + 3 * x59) + x46 * x74
    x81 = x72 * x80

    # 15 item(s)
    S = numpy.array(
        [
            x23
            * (
                x15 * x37
                + x3 * (3 * x28 + 3 * x30 + 3 * x34)
                + x8 * (2 * b * x19 - 4 * x13 - 4 * x16)
            )
            + x24 * x26
            + x24 * x7,
            x18 * x44 + x26 * x38 * x39 + x38 * x45,
            x18 * x40 * x48 + x39 * x46 * x7 + x45 * x46,
            x26 * x51 * x52 + x36 * x43 * x51 + x52 * x57,
            x17 * x41 * x48 + x17 * x44 * x46 + x23 * x36 * x38 * x46,
            x36 * x40 * x60 + x60 * x61 * x7 + x61 * x66,
            x15 * x71 + x26 * x67 * x69 + x33 * x43 * x67,
            x31 * x48 * x51 + x33 * x47 * x51 + x46 * x57 * x69,
            x31 * x42 * x60 + x33 * x41 * x60 + x38 * x66 * x73,
            x15 * x76 + x33 * x40 * x74 + x7 * x73 * x74,
            x26 * x78
            + x27 * x78
            + x68
            * (
                x3 * (3 * x53 + 3 * x55 + 3 * x56)
                + x38 * x70
                + x8 * (2 * b * x77 - 4 * x49 - 4 * x50)
            ),
            x12 * x48 * x67 + x27 * x46 * x67 * x68 + x46 * x71,
            x12 * x57 * x60 + x27 * x60 * x79 + x66 * x79,
            x12 * x42 * x74 + x27 * x38 * x72 * x74 + x38 * x76,
            x27 * x81
            + x7 * x81
            + x72
            * (
                x3 * (3 * x62 + 3 * x64 + 3 * x65)
                + x46 * x75
                + x8 * (2 * b * x80 - 4 * x58 - 4 * x59)
            ),
        ]
    )
    return S


def kinetic3d_10(a, A, b, B):
    """Cartesian 3D (ps) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = x0 * (a * A[0] + b * B[0]) - A[0]
    x2 = b * x0
    x3 = a * x2
    x4 = numpy.exp(-x3 * (A[0] - B[0]) ** 2)
    x5 = x1 * x4
    x6 = numpy.sqrt(x0)
    x7 = numpy.pi ** (3 / 2) * x6
    x8 = x5 * x7
    x9 = numpy.exp(-x3 * (A[2] - B[2]) ** 2)
    x10 = numpy.exp(-x3 * (A[1] - B[1]) ** 2)
    x11 = -a
    x12 = 2 * a
    x13 = (2 * b + x12) ** (-1.0)
    x14 = x0 * (a * A[1] + b * B[1]) - A[1]
    x15 = 2 * a**2
    x16 = x10 * (-x11 - x15 * (x13 + x14**2))
    x17 = x0 * x16 * x9
    x18 = x0 * (a * A[2] + b * B[2]) - A[2]
    x19 = x9 * (-x11 - x15 * (x13 + x18**2))
    x20 = x0 * x10
    x21 = numpy.sqrt(numpy.pi) * x6
    x22 = x21 * x5
    x23 = x12 * x2
    x24 = -x11 - x15 * (x1**2 + x13)
    x25 = numpy.pi * x0 * x9
    x26 = x20 * x4
    x27 = x14 * x26 * x7
    x28 = x24 * x9
    x29 = x14 * x21
    x30 = x18 * x7
    x31 = x18 * x21

    # 3 item(s)
    S = numpy.array(
        [
            x10 * x25 * (x22 * x23 + x22 * x24) + x17 * x8 + x19 * x20 * x8,
            x19 * x27 + x25 * x4 * (x10 * x23 * x29 + x16 * x29) + x27 * x28,
            x17 * x30 * x4
            + x26 * x28 * x30
            + numpy.pi * x26 * (x19 * x31 + x23 * x31 * x9),
        ]
    )
    return S


def kinetic3d_11(a, A, b, B):
    """Cartesian 3D (pp) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = (2 * b + x1) ** (-1.0)
    x3 = (a + b) ** (-1.0)
    x4 = -x3 * (a * A[1] + b * B[1])
    x5 = -x4 - A[1]
    x6 = 2 * a**2
    x7 = -x0 - x6 * (x2 + x5**2)
    x8 = b * x3
    x9 = a * x8
    x10 = numpy.exp(-x9 * (A[0] - B[0]) ** 2)
    x11 = numpy.sqrt(x3)
    x12 = numpy.sqrt(numpy.pi) * x11
    x13 = x12 * x2
    x14 = x10 * x13
    x15 = -x3 * (a * A[0] + b * B[0])
    x16 = -x15 - A[0]
    x17 = x10 * (-x15 - B[0])
    x18 = x12 * x17
    x19 = x14 + x16 * x18
    x20 = numpy.exp(-x9 * (A[1] - B[1]) ** 2)
    x21 = numpy.exp(-x9 * (A[2] - B[2]) ** 2)
    x22 = numpy.pi * x21 * x3
    x23 = x20 * x22
    x24 = x19 * x23
    x25 = -x3 * (a * A[2] + b * B[2])
    x26 = -x25 - A[2]
    x27 = -x0 - x6 * (x2 + x26**2)
    x28 = -x0 - x6 * (x16**2 + x2)
    x29 = x1 * x8
    x30 = x18 * x28 + x18 * x29
    x31 = -x4 - B[1]
    x32 = x20 * x31
    x33 = x3 * x32
    x34 = x10 * x16
    x35 = numpy.pi ** (3 / 2) * x11 * x21
    x36 = x34 * x35
    x37 = x12 * x32
    x38 = x12 * x7
    x39 = x29 * x37 + x32 * x38
    x40 = x12 * x34
    x41 = x23 * (x28 * x40 + x29 * x40)
    x42 = -x25 - B[2]
    x43 = x20 * x3
    x44 = x43 * x7
    x45 = x12 * x21
    x46 = x42 * x45
    x47 = x27 * x45
    x48 = x29 * x46 + x42 * x47
    x49 = numpy.pi * x43
    x50 = x20 * x5
    x51 = x3 * x50
    x52 = x35 * x51
    x53 = x23 * x30
    x54 = x12 * x29 * x50 + x38 * x50
    x55 = x13 * x20
    x56 = x37 * x5 + x55
    x57 = x10 * x22
    x58 = x56 * x57
    x59 = x10 * x28
    x60 = numpy.pi * x10
    x61 = x26 * x35
    x62 = x26 * x29 * x45 + x26 * x47
    x63 = x13 * x21
    x64 = x26 * x46 + x63
    x65 = x10 * x49
    x66 = x64 * x65

    # 9 item(s)
    S = numpy.array(
        [
            x23 * (x14 * x28 + x16 * x30 + x19 * x29) + x24 * x27 + x24 * x7,
            x22 * x34 * x39 + x27 * x33 * x36 + x31 * x41,
            x34 * x48 * x49 + x36 * x42 * x44 + x41 * x42,
            x17 * x22 * x54 + x17 * x27 * x52 + x5 * x53,
            x27 * x58 + x28 * x58 + x57 * (x29 * x56 + x39 * x5 + x55 * x7),
            x42 * x52 * x59 + x42 * x54 * x57 + x48 * x51 * x60,
            x17 * x44 * x61 + x17 * x49 * x62 + x26 * x53,
            x26 * x39 * x57 + x33 * x59 * x61 + x33 * x60 * x62,
            x28 * x66 + x65 * (x26 * x48 + x27 * x63 + x29 * x64) + x66 * x7,
        ]
    )
    return S


def kinetic3d_12(a, A, b, B):
    """Cartesian 3D (pd) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = -x5 - A[1]
    x7 = 2 * a**2
    x8 = -x0 - x7 * (x3 + x6**2)
    x9 = a * x4
    x10 = b * x9
    x11 = numpy.exp(-x10 * (A[2] - B[2]) ** 2)
    x12 = -x4 * (a * A[0] + b * B[0])
    x13 = -x12 - B[0]
    x14 = numpy.exp(-x10 * (A[0] - B[0]) ** 2)
    x15 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x16 = x14 * x15
    x17 = x16 * x3
    x18 = -x12 - A[0]
    x19 = x13**2 * x16 + x17
    x20 = x4 * (2 * x13 * x17 + x18 * x19)
    x21 = numpy.exp(-x10 * (A[1] - B[1]) ** 2)
    x22 = numpy.pi * x21
    x23 = x11 * x20 * x22
    x24 = -x4 * (a * A[2] + b * B[2])
    x25 = -x24 - A[2]
    x26 = -x0 - x7 * (x25**2 + x3)
    x27 = x13 * x16
    x28 = 4 * x10
    x29 = -x0 - x7 * (x18**2 + x3)
    x30 = 2 * x16
    x31 = b * x1
    x32 = x17 * x29
    x33 = x31 * x4
    x34 = x27 * x29 + x27 * x33
    x35 = x13 * x34 + x32 + x9 * (x19 * x2 - x30)
    x36 = numpy.pi * x11 * x4
    x37 = x21 * x36
    x38 = -x5 - B[1]
    x39 = x17 + x18 * x27
    x40 = x37 * x39
    x41 = x15 * x21
    x42 = x38 * x41
    x43 = x41 * x8
    x44 = x33 * x42 + x38 * x43
    x45 = x11 * x15
    x46 = x37 * (x18 * x34 + x32 + x33 * x39)
    x47 = -x24 - B[2]
    x48 = x45 * x47
    x49 = x26 * x45
    x50 = x33 * x48 + x47 * x49
    x51 = x3 * x41
    x52 = x38**2 * x41 + x51
    x53 = x14 * x36
    x54 = x18 * x53
    x55 = x16 * x18
    x56 = x29 * x55 + x33 * x55
    x57 = x51 * x8
    x58 = 2 * x41
    x59 = x38 * x44 + x57 + x9 * (x2 * x52 - x58)
    x60 = x53 * x59
    x61 = x14 * x22 * x4
    x62 = x18 * x61
    x63 = x37 * x47
    x64 = x3 * x45
    x65 = x45 * x47**2 + x64
    x66 = x26 * x64
    x67 = 2 * x45
    x68 = x47 * x50 + x66 + x9 * (x2 * x65 - x67)
    x69 = x61 * x68
    x70 = x19 * x37
    x71 = x33 * x41 * x6 + x43 * x6
    x72 = x35 * x37
    x73 = x42 * x6 + x51
    x74 = x13 * x53
    x75 = x53 * (x33 * x73 + x44 * x6 + x57)
    x76 = x6 * x61
    x77 = 2 * x38 * x51 + x52 * x6
    x78 = x53 * x77
    x79 = x29 * x53
    x80 = x25 * x33 * x45 + x25 * x49
    x81 = x13 * x61
    x82 = x25 * x48 + x64
    x83 = x61 * (x25 * x50 + x33 * x82 + x66)
    x84 = x25 * x65 + 2 * x47 * x64
    x85 = x61 * x84

    # 18 item(s)
    S = numpy.array(
        [
            x23 * x26
            + x23 * x8
            + x37 * (x18 * x35 + x20 * x31 + x3 * (x13 * x29 * x30 + x27 * x28)),
            x26 * x38 * x40 + x38 * x46 + x39 * x44 * x45,
            x39 * x41 * x50 + x40 * x47 * x8 + x46 * x47,
            x18 * x60 + x26 * x52 * x54 + x45 * x52 * x56,
            x38 * x50 * x62 + x38 * x56 * x63 + x44 * x47 * x54,
            x18 * x69 + x41 * x56 * x65 + x62 * x65 * x8,
            x19 * x45 * x71 + x26 * x6 * x70 + x6 * x72,
            x13 * x75 + x26 * x73 * x74 + x34 * x45 * x73,
            x13 * x50 * x76 + x34 * x6 * x63 + x47 * x71 * x74,
            x26 * x78
            + x29 * x78
            + x53 * (x3 * (x28 * x42 + x38 * x58 * x8) + x33 * x77 + x59 * x6),
            x16 * x50 * x73 + x47 * x73 * x79 + x47 * x75,
            x16 * x65 * x71 + x29 * x65 * x76 + x6 * x69,
            x19 * x41 * x80 + x25 * x70 * x8 + x25 * x72,
            x25 * x34 * x37 * x38 + x25 * x44 * x74 + x38 * x80 * x81,
            x13 * x83 + x34 * x41 * x82 + x8 * x81 * x82,
            x16 * x52 * x80 + x25 * x52 * x79 + x25 * x60,
            x16 * x44 * x82 + x29 * x38 * x61 * x82 + x38 * x83,
            x29 * x85
            + x61 * (x25 * x68 + x3 * (x26 * x47 * x67 + x28 * x48) + x33 * x84)
            + x8 * x85,
        ]
    )
    return S


def kinetic3d_13(a, A, b, B):
    """Cartesian 3D (pf) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = -x5 - A[1]
    x7 = 2 * a**2
    x8 = -x0 - x7 * (x3 + x6**2)
    x9 = a * x4
    x10 = b * x9
    x11 = numpy.exp(-x10 * (A[2] - B[2]) ** 2)
    x12 = numpy.exp(-x10 * (A[0] - B[0]) ** 2)
    x13 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x14 = x12 * x13
    x15 = x14 * x3
    x16 = -x4 * (a * A[0] + b * B[0])
    x17 = -x16 - B[0]
    x18 = x14 * x17**2
    x19 = -x16 - A[0]
    x20 = 2 * x15 * x17
    x21 = x15 + x18
    x22 = x17 * x21 + x20
    x23 = x4 * (x19 * x22 + x3 * (3 * x15 + 3 * x18))
    x24 = numpy.exp(-x10 * (A[1] - B[1]) ** 2)
    x25 = numpy.pi * x24
    x26 = x11 * x23 * x25
    x27 = -x4 * (a * A[2] + b * B[2])
    x28 = -x27 - A[2]
    x29 = -x0 - x7 * (x28**2 + x3)
    x30 = b * x1
    x31 = -x0 - x7 * (x19**2 + x3)
    x32 = x15 * x31
    x33 = 2 * x14
    x34 = x9 * (x2 * x21 - x33)
    x35 = x14 * x17
    x36 = x30 * x4
    x37 = x31 * x35 + x35 * x36
    x38 = x17 * x37
    x39 = 4 * x10
    x40 = x3 * (x17 * x31 * x33 + x35 * x39)
    x41 = x32 + x34 + x38
    x42 = x17 * x41 + x40 + x9 * (x2 * x22 - 3 * x35)
    x43 = numpy.pi * x11 * x4
    x44 = x24 * x43
    x45 = -x5 - B[1]
    x46 = x19 * x21 + x20
    x47 = x44 * x46
    x48 = x13 * x24
    x49 = x45 * x48
    x50 = x36 * x49 + x49 * x8
    x51 = x11 * x13
    x52 = x50 * x51
    x53 = x44 * (x19 * x41 + x36 * x46 + x40)
    x54 = -x27 - B[2]
    x55 = x51 * x54
    x56 = x29 * x55 + x36 * x55
    x57 = x48 * x56
    x58 = x15 + x19 * x35
    x59 = x3 * x48
    x60 = x45**2 * x48
    x61 = x59 + x60
    x62 = x51 * x61
    x63 = x19 * x37 + x32 + x36 * x58
    x64 = x59 * x8
    x65 = 2 * x48
    x66 = x9 * (x2 * x61 - x65)
    x67 = x45 * x50
    x68 = x64 + x66 + x67
    x69 = x44 * x54
    x70 = x3 * x51
    x71 = x51 * x54**2
    x72 = x70 + x71
    x73 = x48 * x72
    x74 = x29 * x70
    x75 = 2 * x51
    x76 = x9 * (x2 * x72 - x75)
    x77 = x54 * x56
    x78 = x74 + x76 + x77
    x79 = 2 * x45 * x59
    x80 = x45 * x61 + x79
    x81 = x12 * x43
    x82 = x19 * x81
    x83 = x14 * x19
    x84 = x31 * x83 + x36 * x83
    x85 = x3 * (x39 * x49 + x45 * x65 * x8)
    x86 = x45 * x68 + x85 + x9 * (x2 * x80 - 3 * x49)
    x87 = x81 * x86
    x88 = x12 * x25 * x4
    x89 = x19 * x88
    x90 = 2 * x54 * x70
    x91 = x54 * x72 + x90
    x92 = x3 * (x29 * x54 * x75 + x39 * x55)
    x93 = x54 * x78 + x9 * (x2 * x91 - 3 * x55) + x92
    x94 = x88 * x93
    x95 = x22 * x44
    x96 = x48 * x8
    x97 = x36 * x48 * x6 + x6 * x96
    x98 = x42 * x44
    x99 = x49 * x6 + x59
    x100 = x21 * x51
    x101 = x36 * x99 + x50 * x6 + x64
    x102 = x6 * x61 + x79
    x103 = x17 * x81
    x104 = x81 * (x102 * x36 + x6 * x68 + x85)
    x105 = x6 * x88
    x106 = x3 * (3 * x59 + 3 * x60) + x6 * x80
    x107 = x106 * x81
    x108 = x31 * x81
    x109 = x14 * x72
    x110 = x28 * x51
    x111 = x110 * x29 + x110 * x36
    x112 = x28 * x55 + x70
    x113 = x112 * x36 + x28 * x56 + x74
    x114 = x17 * x88
    x115 = x28 * x72 + x90
    x116 = x88 * (x115 * x36 + x28 * x78 + x92)
    x117 = x14 * x61
    x118 = x28 * x91 + x3 * (3 * x70 + 3 * x71)
    x119 = x118 * x88

    # 30 item(s)
    S = numpy.array(
        [
            x26 * x29
            + x26 * x8
            + x44 * (x19 * x42 + x23 * x30 + x3 * (3 * x32 + 3 * x34 + 3 * x38)),
            x29 * x45 * x47 + x45 * x53 + x46 * x52,
            x46 * x57 + x47 * x54 * x8 + x53 * x54,
            x29 * x58 * x62 + x51 * x58 * x68 + x62 * x63,
            x45 * x63 * x69 + x49 * x56 * x58 + x52 * x54 * x58,
            x48 * x58 * x78 + x58 * x73 * x8 + x63 * x73,
            x19 * x87 + x29 * x80 * x82 + x51 * x80 * x84,
            x54 * x68 * x82 + x55 * x61 * x84 + x56 * x61 * x83,
            x45 * x78 * x89 + x49 * x72 * x84 + x50 * x72 * x83,
            x19 * x94 + x48 * x84 * x91 + x8 * x89 * x91,
            x22 * x51 * x97 + x29 * x6 * x95 + x6 * x98,
            x100 * x101 + x100 * x29 * x99 + x41 * x51 * x99,
            x21 * x55 * x97 + x21 * x57 * x6 + x41 * x6 * x69,
            x102 * x103 * x29 + x102 * x37 * x51 + x104 * x17,
            x101 * x103 * x54 + x35 * x56 * x99 + x37 * x55 * x99,
            x105 * x17 * x78 + x35 * x72 * x97 + x37 * x6 * x73,
            x107 * x29
            + x107 * x31
            + x81 * (x106 * x36 + x3 * (3 * x64 + 3 * x66 + 3 * x67) + x6 * x86),
            x102 * x108 * x54 + x102 * x14 * x56 + x104 * x54,
            x101 * x109 + x109 * x31 * x99 + x14 * x78 * x99,
            x105 * x31 * x91 + x14 * x91 * x97 + x6 * x94,
            x111 * x22 * x48 + x28 * x8 * x95 + x28 * x98,
            x111 * x21 * x49 + x21 * x28 * x52 + x28 * x41 * x44 * x45,
            x112 * x21 * x96 + x112 * x41 * x48 + x113 * x21 * x48,
            x103 * x28 * x68 + x111 * x35 * x61 + x28 * x37 * x62,
            x112 * x35 * x50 + x112 * x37 * x49 + x113 * x114 * x45,
            x114 * x115 * x8 + x115 * x37 * x48 + x116 * x17,
            x108 * x28 * x80 + x111 * x14 * x80 + x28 * x87,
            x112 * x117 * x31 + x112 * x14 * x68 + x113 * x117,
            x115 * x14 * x50 + x115 * x31 * x45 * x88 + x116 * x45,
            x119 * x31
            + x119 * x8
            + x88 * (x118 * x36 + x28 * x93 + x3 * (3 * x74 + 3 * x76 + 3 * x77)),
        ]
    )
    return S


def kinetic3d_14(a, A, b, B):
    """Cartesian 3D (pg) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = -x5 - A[1]
    x7 = 2 * a**2
    x8 = -x0 - x7 * (x3 + x6**2)
    x9 = a * x4
    x10 = b * x9
    x11 = numpy.exp(-x10 * (A[2] - B[2]) ** 2)
    x12 = -x4 * (a * A[0] + b * B[0])
    x13 = -x12 - B[0]
    x14 = numpy.exp(-x10 * (A[0] - B[0]) ** 2)
    x15 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x16 = x14 * x15
    x17 = x16 * x3
    x18 = x13 * x17
    x19 = x13**2 * x16
    x20 = x17 + x19
    x21 = x13 * x20
    x22 = -x12 - A[0]
    x23 = x3 * (3 * x17 + 3 * x19)
    x24 = 2 * x18
    x25 = x21 + x24
    x26 = x13 * x25 + x23
    x27 = x4 * (x22 * x26 + x3 * (8 * x18 + 4 * x21))
    x28 = numpy.exp(-x10 * (A[1] - B[1]) ** 2)
    x29 = numpy.pi * x28
    x30 = x11 * x27 * x29
    x31 = -x4 * (a * A[2] + b * B[2])
    x32 = -x31 - A[2]
    x33 = -x0 - x7 * (x3 + x32**2)
    x34 = b * x1
    x35 = x13 * x16
    x36 = 4 * x10
    x37 = -x0 - x7 * (x22**2 + x3)
    x38 = 2 * x16
    x39 = x3 * (x13 * x37 * x38 + x35 * x36)
    x40 = x9 * (x2 * x25 - 3 * x35)
    x41 = x17 * x37
    x42 = x9 * (x2 * x20 - x38)
    x43 = x34 * x4
    x44 = x35 * x37 + x35 * x43
    x45 = x13 * x44
    x46 = x41 + x42 + x45
    x47 = x13 * x46
    x48 = x3 * (3 * x41 + 3 * x42 + 3 * x45)
    x49 = x39 + x40 + x47
    x50 = x13 * x49 + x48 + x9 * (2 * b * x26 - 4 * x17 - 4 * x19)
    x51 = numpy.pi * x11 * x4
    x52 = x28 * x51
    x53 = -x5 - B[1]
    x54 = x22 * x25 + x23
    x55 = x52 * x54
    x56 = x15 * x28
    x57 = x53 * x56
    x58 = x43 * x57 + x57 * x8
    x59 = x11 * x15
    x60 = x58 * x59
    x61 = x52 * (x22 * x49 + x43 * x54 + x48)
    x62 = -x31 - B[2]
    x63 = x59 * x62
    x64 = x33 * x63 + x43 * x63
    x65 = x56 * x64
    x66 = x3 * x56
    x67 = x53**2 * x56
    x68 = x66 + x67
    x69 = x20 * x22 + x24
    x70 = x59 * x69
    x71 = x66 * x8
    x72 = 2 * x56
    x73 = x9 * (x2 * x68 - x72)
    x74 = x53 * x58
    x75 = x71 + x73 + x74
    x76 = x22 * x46 + x39 + x43 * x69
    x77 = x59 * x68
    x78 = x52 * x62
    x79 = x3 * x59
    x80 = x59 * x62**2
    x81 = x79 + x80
    x82 = x56 * x69
    x83 = x33 * x79
    x84 = 2 * x59
    x85 = x9 * (x2 * x81 - x84)
    x86 = x62 * x64
    x87 = x83 + x85 + x86
    x88 = x56 * x81
    x89 = x17 + x22 * x35
    x90 = x53 * x66
    x91 = 2 * x90
    x92 = x53 * x68
    x93 = x91 + x92
    x94 = x59 * x93
    x95 = x22 * x44 + x41 + x43 * x89
    x96 = x3 * (x36 * x57 + x53 * x72 * x8)
    x97 = x9 * (x2 * x93 - 3 * x57)
    x98 = x53 * x75
    x99 = x96 + x97 + x98
    x100 = x62 * x79
    x101 = 2 * x100
    x102 = x62 * x81
    x103 = x101 + x102
    x104 = x103 * x56
    x105 = x3 * (x33 * x62 * x84 + x36 * x63)
    x106 = x9 * (x103 * x2 - 3 * x63)
    x107 = x62 * x87
    x108 = x105 + x106 + x107
    x109 = x3 * (3 * x66 + 3 * x67)
    x110 = x109 + x53 * x93
    x111 = x14 * x51
    x112 = x111 * x22
    x113 = x16 * x22
    x114 = x113 * x37 + x113 * x43
    x115 = x3 * (3 * x71 + 3 * x73 + 3 * x74)
    x116 = x115 + x53 * x99 + x9 * (2 * b * x110 - 4 * x66 - 4 * x67)
    x117 = x111 * x116
    x118 = x14 * x29 * x4
    x119 = x118 * x22
    x120 = x3 * (3 * x79 + 3 * x80)
    x121 = x103 * x62 + x120
    x122 = x3 * (3 * x83 + 3 * x85 + 3 * x86)
    x123 = x108 * x62 + x122 + x9 * (2 * b * x121 - 4 * x79 - 4 * x80)
    x124 = x118 * x123
    x125 = x26 * x52
    x126 = x56 * x6
    x127 = x126 * x43 + x126 * x8
    x128 = x50 * x52
    x129 = x57 * x6 + x66
    x130 = x25 * x59
    x131 = x129 * x43 + x58 * x6 + x71
    x132 = x6 * x68 + x91
    x133 = x132 * x59
    x134 = x132 * x43 + x6 * x75 + x96
    x135 = x20 * x59
    x136 = x109 + x6 * x93
    x137 = x111 * x13
    x138 = x111 * (x115 + x136 * x43 + x6 * x99)
    x139 = x118 * x6
    x140 = x110 * x6 + x3 * (8 * x90 + 4 * x92)
    x141 = x111 * x140
    x142 = x111 * x37
    x143 = x132 * x16
    x144 = x103 * x16
    x145 = x32 * x59
    x146 = x145 * x33 + x145 * x43
    x147 = x32 * x63 + x79
    x148 = x25 * x56
    x149 = x147 * x43 + x32 * x64 + x83
    x150 = x101 + x32 * x81
    x151 = x150 * x56
    x152 = x105 + x150 * x43 + x32 * x87
    x153 = x118 * x13
    x154 = x103 * x32 + x120
    x155 = x118 * (x108 * x32 + x122 + x154 * x43)
    x156 = x16 * x93
    x157 = x150 * x16
    x158 = x121 * x32 + x3 * (8 * x100 + 4 * x102)
    x159 = x118 * x158

    # 45 item(s)
    S = numpy.array(
        [
            x30 * x33
            + x30 * x8
            + x52 * (x22 * x50 + x27 * x34 + x3 * (4 * x39 + 4 * x40 + 4 * x47)),
            x33 * x53 * x55 + x53 * x61 + x54 * x60,
            x54 * x65 + x55 * x62 * x8 + x61 * x62,
            x33 * x68 * x70 + x70 * x75 + x76 * x77,
            x53 * x76 * x78 + x57 * x64 * x69 + x60 * x62 * x69,
            x76 * x88 + x8 * x81 * x82 + x82 * x87,
            x33 * x89 * x94 + x59 * x89 * x99 + x94 * x95,
            x63 * x68 * x95 + x63 * x75 * x89 + x64 * x68 * x89,
            x57 * x81 * x95 + x57 * x87 * x89 + x58 * x81 * x89,
            x104 * x8 * x89 + x104 * x95 + x108 * x56 * x89,
            x110 * x112 * x33 + x110 * x114 * x59 + x117 * x22,
            x112 * x62 * x99 + x113 * x64 * x93 + x114 * x63 * x93,
            x113 * x68 * x87 + x113 * x75 * x81 + x114 * x68 * x81,
            x103 * x113 * x58 + x103 * x114 * x57 + x108 * x119 * x53,
            x114 * x121 * x56 + x119 * x121 * x8 + x124 * x22,
            x125 * x33 * x6 + x127 * x26 * x59 + x128 * x6,
            x129 * x130 * x33 + x129 * x49 * x59 + x130 * x131,
            x127 * x25 * x63 + x25 * x6 * x65 + x49 * x6 * x78,
            x133 * x20 * x33 + x133 * x46 + x134 * x135,
            x129 * x20 * x64 + x129 * x46 * x63 + x131 * x20 * x63,
            x126 * x20 * x87 + x127 * x20 * x81 + x46 * x6 * x88,
            x13 * x138 + x136 * x137 * x33 + x136 * x44 * x59,
            x132 * x35 * x64 + x132 * x44 * x63 + x134 * x137 * x62,
            x129 * x35 * x87 + x129 * x44 * x81 + x131 * x35 * x81,
            x103 * x127 * x35 + x104 * x44 * x6 + x108 * x13 * x139,
            x111 * (x116 * x6 + x140 * x43 + x3 * (4 * x96 + 4 * x97 + 4 * x98))
            + x141 * x33
            + x141 * x37,
            x136 * x142 * x62 + x136 * x16 * x64 + x138 * x62,
            x134 * x16 * x81 + x143 * x37 * x81 + x143 * x87,
            x108 * x129 * x16 + x129 * x144 * x37 + x131 * x144,
            x121 * x127 * x16 + x121 * x139 * x37 + x124 * x6,
            x125 * x32 * x8 + x128 * x32 + x146 * x26 * x56,
            x146 * x25 * x57 + x25 * x32 * x60 + x32 * x49 * x52 * x53,
            x147 * x148 * x8 + x147 * x49 * x56 + x148 * x149,
            x135 * x32 * x75 + x146 * x20 * x68 + x32 * x46 * x77,
            x147 * x20 * x58 + x147 * x46 * x57 + x149 * x20 * x57,
            x151 * x20 * x8 + x151 * x46 + x152 * x20 * x56,
            x137 * x32 * x99 + x146 * x35 * x93 + x32 * x44 * x94,
            x147 * x35 * x75 + x147 * x44 * x68 + x149 * x35 * x68,
            x150 * x35 * x58 + x150 * x44 * x57 + x152 * x153 * x53,
            x13 * x155 + x153 * x154 * x8 + x154 * x44 * x56,
            x110 * x142 * x32 + x110 * x146 * x16 + x117 * x32,
            x147 * x156 * x37 + x147 * x16 * x99 + x149 * x156,
            x152 * x16 * x68 + x157 * x37 * x68 + x157 * x75,
            x118 * x154 * x37 * x53 + x154 * x16 * x58 + x155 * x53,
            x118 * (x123 * x32 + x158 * x43 + x3 * (4 * x105 + 4 * x106 + 4 * x107))
            + x159 * x37
            + x159 * x8,
        ]
    )
    return S


def kinetic3d_20(a, A, b, B):
    """Cartesian 3D (ds) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = (2 * b + x1) ** (-1.0)
    x3 = (a + b) ** (-1.0)
    x4 = x3 * (a * A[1] + b * B[1]) - A[1]
    x5 = x4**2
    x6 = 2 * a**2
    x7 = -x0 - x6 * (x2 + x5)
    x8 = b * x3
    x9 = a * x8
    x10 = numpy.exp(-x9 * (A[0] - B[0]) ** 2)
    x11 = numpy.sqrt(x3)
    x12 = numpy.sqrt(numpy.pi) * x11
    x13 = x10 * x12
    x14 = x13 * x2
    x15 = x3 * (a * A[0] + b * B[0]) - A[0]
    x16 = x15**2
    x17 = x13 * x16 + x14
    x18 = numpy.exp(-x9 * (A[1] - B[1]) ** 2)
    x19 = numpy.exp(-x9 * (A[2] - B[2]) ** 2)
    x20 = numpy.pi * x19 * x3
    x21 = x18 * x20
    x22 = x17 * x21
    x23 = x3 * (a * A[2] + b * B[2]) - A[2]
    x24 = x23**2
    x25 = -x0 - x6 * (x2 + x24)
    x26 = -x0 - x6 * (x16 + x2)
    x27 = x13 * x15
    x28 = x1 * x8
    x29 = x26 * x27 + x27 * x28
    x30 = numpy.pi ** (3 / 2)
    x31 = x10 * x18 * x3
    x32 = x11 * x15 * x19 * x30 * x31
    x33 = x12 * x18
    x34 = x33 * x4
    x35 = x28 * x34 + x34 * x7
    x36 = x10 * x20
    x37 = x35 * x36
    x38 = x21 * x29
    x39 = x12 * x19
    x40 = x23 * x39
    x41 = x25 * x40 + x28 * x40
    x42 = numpy.pi * x31
    x43 = x41 * x42
    x44 = x2 * x33
    x45 = x33 * x5 + x44
    x46 = x36 * x45
    x47 = x2 * x39
    x48 = x24 * x39 + x47
    x49 = x42 * x48

    # 6 item(s)
    S = numpy.array(
        [
            x21 * (x14 * x26 + x15 * x29 + x8 * (x1 * x17 - 2 * x13))
            + x22 * x25
            + x22 * x7,
            x15 * x37 + x25 * x32 * x4 + x38 * x4,
            x15 * x43 + x23 * x32 * x7 + x23 * x38,
            x25 * x46
            + x26 * x46
            + x36 * (x35 * x4 + x44 * x7 + x8 * (x1 * x45 - 2 * x33)),
            x11 * x19 * x23 * x26 * x30 * x31 * x4 + x23 * x37 + x4 * x43,
            x26 * x49
            + x42 * (x23 * x41 + x25 * x47 + x8 * (x1 * x48 - 2 * x39))
            + x49 * x7,
        ]
    )
    return S


def kinetic3d_21(a, A, b, B):
    """Cartesian 3D (dp) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = (2 * b + x1) ** (-1.0)
    x3 = (a + b) ** (-1.0)
    x4 = -x3 * (a * A[1] + b * B[1])
    x5 = -x4 - A[1]
    x6 = x5**2
    x7 = 2 * a**2
    x8 = -x0 - x7 * (x2 + x6)
    x9 = -x3 * (a * A[0] + b * B[0])
    x10 = -x9 - A[0]
    x11 = b * x3
    x12 = a * x11
    x13 = numpy.exp(-x12 * (A[0] - B[0]) ** 2)
    x14 = numpy.sqrt(numpy.pi) * numpy.sqrt(x3)
    x15 = x13 * x14
    x16 = x10 * x15
    x17 = -x9 - B[0]
    x18 = x15 * x17
    x19 = x15 * x2
    x20 = x16 * x17 + x19
    x21 = x10 * x20 + x2 * (x16 + x18)
    x22 = numpy.exp(-x12 * (A[1] - B[1]) ** 2)
    x23 = numpy.exp(-x12 * (A[2] - B[2]) ** 2)
    x24 = numpy.pi * x23 * x3
    x25 = x22 * x24
    x26 = x21 * x25
    x27 = -x3 * (a * A[2] + b * B[2])
    x28 = -x27 - A[2]
    x29 = x28**2
    x30 = -x0 - x7 * (x2 + x29)
    x31 = x1 * x11
    x32 = x10**2
    x33 = -x0 - x7 * (x2 + x32)
    x34 = x18 * x31 + x18 * x33
    x35 = x16 * x31 + x16 * x33
    x36 = x19 * x33
    x37 = x10 * x34 + x20 * x31 + x36
    x38 = -x4 - B[1]
    x39 = x15 * x32 + x19
    x40 = x25 * x39
    x41 = x14 * x22
    x42 = x38 * x41
    x43 = x31 * x42 + x42 * x8
    x44 = x14 * x23
    x45 = x25 * (x10 * x35 + x11 * (x1 * x39 - 2 * x15) + x36)
    x46 = -x27 - B[2]
    x47 = x44 * x46
    x48 = x30 * x47 + x31 * x47
    x49 = x25 * x5
    x50 = x41 * x5
    x51 = x31 * x50 + x50 * x8
    x52 = x25 * x37
    x53 = x14 * x2
    x54 = x22 * x53
    x55 = x38 * x50 + x54
    x56 = x13 * x24
    x57 = x10 * x56
    x58 = x54 * x8
    x59 = x31 * x55 + x43 * x5 + x58
    x60 = x56 * x59
    x61 = numpy.pi * x13 * x22 * x3
    x62 = x10 * x61
    x63 = x25 * x28
    x64 = x28 * x44
    x65 = x30 * x64 + x31 * x64
    x66 = x23 * x53
    x67 = x46 * x64 + x66
    x68 = x30 * x66
    x69 = x28 * x48 + x31 * x67 + x68
    x70 = x61 * x69
    x71 = x41 * x6 + x54
    x72 = x56 * x71
    x73 = 2 * x41
    x74 = x56 * (x11 * (x1 * x71 - x73) + x5 * x51 + x58)
    x75 = x2 * (x42 + x50) + x5 * x55
    x76 = x56 * x75
    x77 = x5 * x61
    x78 = x28 * x56
    x79 = x29 * x44 + x66
    x80 = x61 * x79
    x81 = 2 * x44
    x82 = x61 * (x11 * (x1 * x79 - x81) + x28 * x65 + x68)
    x83 = x2 * (x47 + x64) + x28 * x67
    x84 = x61 * x83

    # 18 item(s)
    S = numpy.array(
        [
            x25 * (x10 * x37 + x11 * (x1 * x21 - 2 * x18) + x2 * (x34 + x35))
            + x26 * x30
            + x26 * x8,
            x30 * x38 * x40 + x38 * x45 + x39 * x43 * x44,
            x39 * x41 * x48 + x40 * x46 * x8 + x45 * x46,
            x20 * x30 * x49 + x20 * x44 * x51 + x5 * x52,
            x10 * x60 + x30 * x55 * x57 + x35 * x44 * x55,
            x35 * x46 * x49 + x46 * x51 * x57 + x48 * x5 * x62,
            x20 * x41 * x65 + x20 * x63 * x8 + x28 * x52,
            x28 * x43 * x57 + x35 * x38 * x63 + x38 * x62 * x65,
            x10 * x70 + x35 * x41 * x67 + x62 * x67 * x8,
            x17 * x30 * x72 + x17 * x74 + x34 * x44 * x71,
            x30 * x76
            + x33 * x76
            + x56 * (x11 * (x1 * x75 - x38 * x73) + x2 * (x43 + x51) + x5 * x59),
            x15 * x48 * x71 + x33 * x46 * x72 + x46 * x74,
            x17 * x51 * x78 + x17 * x65 * x77 + x28 * x34 * x49,
            x15 * x55 * x65 + x28 * x60 + x33 * x55 * x78,
            x15 * x51 * x67 + x33 * x67 * x77 + x5 * x70,
            x17 * x8 * x80 + x17 * x82 + x34 * x41 * x79,
            x15 * x43 * x79 + x33 * x38 * x80 + x38 * x82,
            x33 * x84
            + x61 * (x11 * (x1 * x83 - x46 * x81) + x2 * (x48 + x65) + x28 * x69)
            + x8 * x84,
        ]
    )
    return S


def kinetic3d_22(a, A, b, B):
    """Cartesian 3D (dd) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = -x5 - A[1]
    x7 = x6**2
    x8 = 2 * a**2
    x9 = -x0 - x8 * (x3 + x7)
    x10 = a * x4
    x11 = b * x10
    x12 = numpy.exp(-x11 * (A[0] - B[0]) ** 2)
    x13 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x14 = x12 * x13
    x15 = x14 * x3
    x16 = -x4 * (a * A[0] + b * B[0])
    x17 = -x16 - B[0]
    x18 = x14 * x17**2
    x19 = -x16 - A[0]
    x20 = 2 * x14
    x21 = x17 * x20
    x22 = 2 * x15
    x23 = x15 + x18
    x24 = x17 * x22 + x19 * x23
    x25 = x19 * x24 + x3 * (3 * x15 + x18 + x19 * x21)
    x26 = numpy.exp(-x11 * (A[1] - B[1]) ** 2)
    x27 = numpy.exp(-x11 * (A[2] - B[2]) ** 2)
    x28 = numpy.pi * x27 * x4
    x29 = x26 * x28
    x30 = x25 * x29
    x31 = -x4 * (a * A[2] + b * B[2])
    x32 = -x31 - A[2]
    x33 = x32**2
    x34 = -x0 - x8 * (x3 + x33)
    x35 = b * x4
    x36 = x14 * x17
    x37 = x1 * x35
    x38 = x19**2
    x39 = -x0 - x8 * (x3 + x38)
    x40 = x36 * x37 + x36 * x39
    x41 = x19 * x40
    x42 = x14 * x19
    x43 = x15 + x17 * x42
    x44 = 4 * x11
    x45 = x15 * x39
    x46 = -x20
    x47 = x10 * (x2 * x23 + x46) + x17 * x40
    x48 = x45 + x47
    x49 = x19 * x48 + x24 * x37 + x3 * (x21 * x39 + x36 * x44)
    x50 = -x5 - B[1]
    x51 = x19 * x43 + x3 * (x36 + x42)
    x52 = x29 * x51
    x53 = x13 * x26
    x54 = x50 * x53
    x55 = x37 * x54 + x54 * x9
    x56 = x13 * x27
    x57 = x37 * x42 + x39 * x42
    x58 = x37 * x43 + x41 + x45
    x59 = x29 * (x19 * x58 + x3 * (x40 + x57) + x35 * (x1 * x51 - x21))
    x60 = -x31 - B[2]
    x61 = x56 * x60
    x62 = x34 * x61 + x37 * x61
    x63 = x3 * x53
    x64 = x50**2 * x53
    x65 = x63 + x64
    x66 = x14 * x38 + x15
    x67 = x56 * x66
    x68 = x63 * x9
    x69 = 2 * x53
    x70 = -x69
    x71 = x10 * (x2 * x65 + x70) + x50 * x55
    x72 = x68 + x71
    x73 = x19 * x57 + x35 * (x1 * x66 + x46) + x45
    x74 = x29 * x60
    x75 = x3 * x56
    x76 = x56 * x60**2
    x77 = x75 + x76
    x78 = x53 * x66
    x79 = x34 * x75
    x80 = 2 * x56
    x81 = -x80
    x82 = x10 * (x2 * x77 + x81) + x60 * x62
    x83 = x79 + x82
    x84 = x24 * x29
    x85 = x53 * x6
    x86 = x37 * x85 + x85 * x9
    x87 = x29 * x49
    x88 = x50 * x85 + x63
    x89 = x43 * x56
    x90 = x55 * x6
    x91 = x37 * x88 + x68 + x90
    x92 = 2 * x63
    x93 = x50 * x92 + x6 * x65
    x94 = x12 * x28
    x95 = x19 * x94
    x96 = x50 * x69
    x97 = x3 * (x44 * x54 + x9 * x96) + x37 * x93 + x6 * x72
    x98 = x94 * x97
    x99 = numpy.pi * x12 * x26 * x4
    x100 = x19 * x99
    x101 = x32 * x56
    x102 = x101 * x34 + x101 * x37
    x103 = x29 * x32
    x104 = x101 * x60 + x75
    x105 = x43 * x53
    x106 = x32 * x62
    x107 = x104 * x37 + x106 + x79
    x108 = 2 * x75
    x109 = x108 * x60 + x32 * x77
    x110 = x60 * x80
    x111 = x109 * x37 + x3 * (x110 * x34 + x44 * x61) + x32 * x83
    x112 = x111 * x99
    x113 = x53 * x7 + x63
    x114 = x113 * x56
    x115 = x35 * (x1 * x113 + x70) + x6 * x86 + x68
    x116 = x3 * (x54 + x85) + x6 * x88
    x117 = x17 * x94
    x118 = x94 * (x3 * (x55 + x86) + x35 * (x1 * x116 - x96) + x6 * x91)
    x119 = x3 * (x6 * x96 + 3 * x63 + x64) + x6 * x93
    x120 = x119 * x94
    x121 = x39 * x94
    x122 = x113 * x14
    x123 = x6 * x99
    x124 = x14 * x88
    x125 = x33 * x56 + x75
    x126 = x125 * x53
    x127 = x102 * x32 + x35 * (x1 * x125 + x81) + x79
    x128 = x17 * x99
    x129 = x104 * x32 + x3 * (x101 + x61)
    x130 = x99 * (x107 * x32 + x3 * (x102 + x62) + x35 * (x1 * x129 - x110))
    x131 = x125 * x14
    x132 = x109 * x32 + x3 * (x110 * x32 + 3 * x75 + x76)
    x133 = x132 * x99

    # 36 item(s)
    S = numpy.array(
        [
            x29
            * (
                x19 * x49
                + x3 * (2 * x41 + x43 * x44 + 3 * x45 + x47)
                + x35 * (2 * a * x25 - 2 * x18 - x22)
            )
            + x30 * x34
            + x30 * x9,
            x34 * x50 * x52 + x50 * x59 + x51 * x55 * x56,
            x51 * x53 * x62 + x52 * x60 * x9 + x59 * x60,
            x34 * x65 * x67 + x56 * x65 * x73 + x67 * x72,
            x50 * x73 * x74 + x54 * x62 * x66 + x55 * x61 * x66,
            x53 * x73 * x77 + x77 * x78 * x9 + x78 * x83,
            x24 * x56 * x86 + x34 * x6 * x84 + x6 * x87,
            x34 * x88 * x89 + x56 * x58 * x88 + x89 * x91,
            x43 * x61 * x86 + x43 * x62 * x85 + x58 * x6 * x74,
            x19 * x98 + x34 * x93 * x95 + x56 * x57 * x93,
            x42 * x62 * x88 + x57 * x61 * x88 + x60 * x91 * x95,
            x100 * x6 * x83 + x42 * x77 * x86 + x57 * x77 * x85,
            x102 * x24 * x53 + x32 * x84 * x9 + x32 * x87,
            x101 * x43 * x55 + x102 * x43 * x54 + x103 * x50 * x58,
            x104 * x105 * x9 + x104 * x53 * x58 + x105 * x107,
            x101 * x57 * x65 + x102 * x42 * x65 + x32 * x72 * x95,
            x100 * x107 * x50 + x104 * x42 * x55 + x104 * x54 * x57,
            x100 * x109 * x9 + x109 * x53 * x57 + x112 * x19,
            x114 * x23 * x34 + x114 * x48 + x115 * x23 * x56,
            x116 * x117 * x34 + x116 * x40 * x56 + x118 * x17,
            x113 * x36 * x62 + x113 * x40 * x61 + x115 * x117 * x60,
            x120 * x34
            + x120 * x39
            + x94
            * (
                x3 * (x44 * x88 + 3 * x68 + x71 + 2 * x90)
                + x35 * (2 * a * x119 - 2 * x64 - x92)
                + x6 * x97
            ),
            x116 * x121 * x60 + x116 * x14 * x62 + x118 * x60,
            x115 * x14 * x77 + x122 * x39 * x77 + x122 * x83,
            x101 * x23 * x86 + x102 * x23 * x85 + x103 * x48 * x6,
            x101 * x40 * x88 + x102 * x36 * x88 + x117 * x32 * x91,
            x104 * x36 * x86 + x104 * x40 * x85 + x107 * x123 * x17,
            x102 * x14 * x93 + x121 * x32 * x93 + x32 * x98,
            x104 * x124 * x39 + x104 * x14 * x91 + x107 * x124,
            x109 * x123 * x39 + x109 * x14 * x86 + x112 * x6,
            x126 * x23 * x9 + x126 * x48 + x127 * x23 * x53,
            x125 * x36 * x55 + x125 * x40 * x54 + x127 * x128 * x50,
            x128 * x129 * x9 + x129 * x40 * x53 + x130 * x17,
            x127 * x14 * x65 + x131 * x39 * x65 + x131 * x72,
            x129 * x14 * x55 + x129 * x39 * x50 * x99 + x130 * x50,
            x133 * x39
            + x133 * x9
            + x99
            * (
                x111 * x32
                + x3 * (x104 * x44 + 2 * x106 + 3 * x79 + x82)
                + x35 * (2 * a * x132 - x108 - 2 * x76)
            ),
        ]
    )
    return S


def kinetic3d_23(a, A, b, B):
    """Cartesian 3D (df) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = -x5 - A[1]
    x7 = x6**2
    x8 = 2 * a**2
    x9 = -x0 - x8 * (x3 + x7)
    x10 = -x4 * (a * A[0] + b * B[0])
    x11 = -x10 - B[0]
    x12 = a * x4
    x13 = b * x12
    x14 = numpy.exp(-x13 * (A[0] - B[0]) ** 2)
    x15 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x16 = x14 * x15
    x17 = x16 * x3
    x18 = x11 * x17
    x19 = x11**2 * x16
    x20 = x17 + x19
    x21 = x11 * x20
    x22 = -x10 - A[0]
    x23 = x20 * x22
    x24 = 3 * x17
    x25 = 2 * x17
    x26 = x11 * x25
    x27 = x21 + x26
    x28 = x22 * x27 + x3 * (3 * x19 + x24)
    x29 = x22 * x28 + x3 * (8 * x18 + x21 + 3 * x23)
    x30 = numpy.exp(-x13 * (A[1] - B[1]) ** 2)
    x31 = numpy.exp(-x13 * (A[2] - B[2]) ** 2)
    x32 = numpy.pi * x31 * x4
    x33 = x30 * x32
    x34 = x29 * x33
    x35 = -x4 * (a * A[2] + b * B[2])
    x36 = -x35 - A[2]
    x37 = x36**2
    x38 = -x0 - x8 * (x3 + x37)
    x39 = b * x4
    x40 = x22**2
    x41 = -x0 - x8 * (x3 + x40)
    x42 = x17 * x41
    x43 = x11 * x16
    x44 = x1 * x39
    x45 = x41 * x43 + x43 * x44
    x46 = x11 * x45
    x47 = 2 * x16
    x48 = -x47
    x49 = x12 * (x2 * x20 + x48)
    x50 = x46 + x49
    x51 = x42 + x50
    x52 = x22 * x51
    x53 = 4 * x13
    x54 = x11 * x47
    x55 = x3 * (x41 * x54 + x43 * x53)
    x56 = x23 + x26
    x57 = 6 * x13
    x58 = x11 * x51 + x12 * (x2 * x27 - 3 * x43)
    x59 = 3 * x42
    x60 = x55 + x58
    x61 = x22 * x60 + x28 * x44 + x3 * (3 * x46 + 3 * x49 + x59)
    x62 = -x5 - B[1]
    x63 = x22 * x56 + x3 * (x19 + x22 * x54 + x24)
    x64 = x33 * x63
    x65 = x15 * x30
    x66 = x62 * x65
    x67 = x44 * x66 + x66 * x9
    x68 = x15 * x31
    x69 = x22 * x45
    x70 = x16 * x22
    x71 = x11 * x70 + x17
    x72 = x44 * x56 + x52 + x55
    x73 = x33 * (
        x22 * x72
        + x3 * (x50 + x53 * x71 + x59 + 2 * x69)
        + x39 * (2 * a * x63 - 2 * x19 - x25)
    )
    x74 = -x35 - B[2]
    x75 = x68 * x74
    x76 = x38 * x75 + x44 * x75
    x77 = x3 * x65
    x78 = x62**2 * x65
    x79 = x77 + x78
    x80 = x22 * x71 + x3 * (x43 + x70)
    x81 = x68 * x80
    x82 = x77 * x9
    x83 = x62 * x67
    x84 = 2 * x65
    x85 = -x84
    x86 = x12 * (x2 * x79 + x85)
    x87 = x83 + x86
    x88 = x82 + x87
    x89 = x41 * x70 + x44 * x70
    x90 = x42 + x44 * x71 + x69
    x91 = x22 * x90 + x3 * (x45 + x89) + x39 * (x1 * x80 - x54)
    x92 = x33 * x74
    x93 = x3 * x68
    x94 = x68 * x74**2
    x95 = x93 + x94
    x96 = x65 * x80
    x97 = x38 * x93
    x98 = x74 * x76
    x99 = 2 * x68
    x100 = -x99
    x101 = x12 * (x100 + x2 * x95)
    x102 = x101 + x98
    x103 = x102 + x97
    x104 = x16 * x40 + x17
    x105 = 2 * x77
    x106 = x105 * x62
    x107 = x62 * x79
    x108 = x106 + x107
    x109 = x108 * x68
    x110 = x22 * x89 + x39 * (x1 * x104 + x48) + x42
    x111 = x62 * x84
    x112 = x3 * (x111 * x9 + x53 * x66)
    x113 = x12 * (x108 * x2 - 3 * x66) + x62 * x88
    x114 = x112 + x113
    x115 = 2 * x93
    x116 = x115 * x74
    x117 = x74 * x95
    x118 = x116 + x117
    x119 = x118 * x65
    x120 = x74 * x99
    x121 = x3 * (x120 * x38 + x53 * x75)
    x122 = x103 * x74 + x12 * (x118 * x2 - 3 * x75)
    x123 = x121 + x122
    x124 = x28 * x33
    x125 = x6 * x65
    x126 = x125 * x44 + x125 * x9
    x127 = x33 * x61
    x128 = x125 * x62 + x77
    x129 = x56 * x68
    x130 = x6 * x67
    x131 = x128 * x44 + x130 + x82
    x132 = x6 * x79
    x133 = x106 + x132
    x134 = x133 * x68
    x135 = x6 * x88
    x136 = x112 + x133 * x44 + x135
    x137 = 3 * x77
    x138 = x108 * x6 + x3 * (x137 + 3 * x78)
    x139 = x14 * x32
    x140 = x139 * x22
    x141 = 3 * x82
    x142 = x114 * x6 + x138 * x44 + x3 * (x141 + 3 * x83 + 3 * x86)
    x143 = x139 * x142
    x144 = numpy.pi * x14 * x30 * x4
    x145 = x144 * x22
    x146 = x36 * x68
    x147 = x146 * x38 + x146 * x44
    x148 = x33 * x36
    x149 = x146 * x74 + x93
    x150 = x56 * x65
    x151 = x36 * x76
    x152 = x149 * x44 + x151 + x97
    x153 = x36 * x95
    x154 = x116 + x153
    x155 = x154 * x65
    x156 = x103 * x36
    x157 = x121 + x154 * x44 + x156
    x158 = 3 * x93
    x159 = x118 * x36 + x3 * (x158 + 3 * x94)
    x160 = 3 * x97
    x161 = x123 * x36 + x159 * x44 + x3 * (3 * x101 + x160 + 3 * x98)
    x162 = x144 * x161
    x163 = x65 * x7 + x77
    x164 = x27 * x68
    x165 = x126 * x6 + x39 * (x1 * x163 + x85) + x82
    x166 = x128 * x6 + x3 * (x125 + x66)
    x167 = x166 * x68
    x168 = x131 * x6 + x3 * (x126 + x67) + x39 * (x1 * x166 - x111)
    x169 = x133 * x6 + x3 * (x111 * x6 + x137 + x78)
    x170 = x11 * x139
    x171 = x139 * (
        x136 * x6
        + x3 * (x128 * x53 + 2 * x130 + x141 + x87)
        + x39 * (2 * a * x169 - x105 - 2 * x78)
    )
    x172 = x62 * x77
    x173 = x138 * x6 + x3 * (x107 + 3 * x132 + 8 * x172)
    x174 = x139 * x173
    x175 = x139 * x41
    x176 = x16 * x166
    x177 = x118 * x16
    x178 = x144 * x6
    x179 = x133 * x16
    x180 = x154 * x16
    x181 = x37 * x68 + x93
    x182 = x27 * x65
    x183 = x147 * x36 + x39 * (x1 * x181 + x100) + x97
    x184 = x149 * x36 + x3 * (x146 + x75)
    x185 = x184 * x65
    x186 = x152 * x36 + x3 * (x147 + x76) + x39 * (x1 * x184 - x120)
    x187 = x11 * x144
    x188 = x154 * x36 + x3 * (x120 * x36 + x158 + x94)
    x189 = x144 * (
        x157 * x36
        + x3 * (x102 + x149 * x53 + 2 * x151 + x160)
        + x39 * (2 * a * x188 - x115 - 2 * x94)
    )
    x190 = x108 * x16
    x191 = x16 * x184
    x192 = x74 * x93
    x193 = x159 * x36 + x3 * (x117 + 3 * x153 + 8 * x192)
    x194 = x144 * x193

    # 60 item(s)
    S = numpy.array(
        [
            x33
            * (
                x22 * x61
                + x3 * (3 * x52 + 4 * x55 + x56 * x57 + x58)
                + x39 * (2 * a * x29 - 4 * x18 - 2 * x21)
            )
            + x34 * x38
            + x34 * x9,
            x38 * x62 * x64 + x62 * x73 + x63 * x67 * x68,
            x63 * x65 * x76 + x64 * x74 * x9 + x73 * x74,
            x38 * x79 * x81 + x68 * x79 * x91 + x81 * x88,
            x62 * x91 * x92 + x66 * x76 * x80 + x67 * x75 * x80,
            x103 * x96 + x65 * x91 * x95 + x9 * x95 * x96,
            x104 * x109 * x38 + x104 * x114 * x68 + x109 * x110,
            x104 * x75 * x88 + x104 * x76 * x79 + x110 * x75 * x79,
            x103 * x104 * x66 + x104 * x67 * x95 + x110 * x66 * x95,
            x104 * x119 * x9 + x104 * x123 * x65 + x110 * x119,
            x124 * x38 * x6 + x126 * x28 * x68 + x127 * x6,
            x128 * x129 * x38 + x128 * x68 * x72 + x129 * x131,
            x125 * x56 * x76 + x126 * x56 * x75 + x6 * x72 * x92,
            x134 * x38 * x71 + x134 * x90 + x136 * x68 * x71,
            x128 * x71 * x76 + x128 * x75 * x90 + x131 * x71 * x75,
            x103 * x125 * x71 + x125 * x90 * x95 + x126 * x71 * x95,
            x138 * x140 * x38 + x138 * x68 * x89 + x143 * x22,
            x133 * x70 * x76 + x133 * x75 * x89 + x136 * x140 * x74,
            x103 * x128 * x70 + x128 * x89 * x95 + x131 * x70 * x95,
            x118 * x125 * x89 + x118 * x126 * x70 + x123 * x145 * x6,
            x124 * x36 * x9 + x127 * x36 + x147 * x28 * x65,
            x146 * x56 * x67 + x147 * x56 * x66 + x148 * x62 * x72,
            x149 * x150 * x9 + x149 * x65 * x72 + x150 * x152,
            x146 * x71 * x88 + x146 * x79 * x90 + x147 * x71 * x79,
            x149 * x66 * x90 + x149 * x67 * x71 + x152 * x66 * x71,
            x155 * x71 * x9 + x155 * x90 + x157 * x65 * x71,
            x108 * x146 * x89 + x108 * x147 * x70 + x114 * x140 * x36,
            x149 * x70 * x88 + x149 * x79 * x89 + x152 * x70 * x79,
            x145 * x157 * x62 + x154 * x66 * x89 + x154 * x67 * x70,
            x145 * x159 * x9 + x159 * x65 * x89 + x162 * x22,
            x163 * x164 * x38 + x163 * x60 * x68 + x164 * x165,
            x167 * x20 * x38 + x167 * x51 + x168 * x20 * x68,
            x163 * x20 * x76 + x163 * x51 * x75 + x165 * x20 * x75,
            x11 * x171 + x169 * x170 * x38 + x169 * x45 * x68,
            x166 * x43 * x76 + x166 * x45 * x75 + x168 * x170 * x74,
            x103 * x163 * x43 + x163 * x45 * x95 + x165 * x43 * x95,
            x139
            * (
                x142 * x6
                + x3 * (4 * x112 + x113 + x133 * x57 + 3 * x135)
                + x39 * (2 * a * x173 - 2 * x107 - 4 * x172)
            )
            + x174 * x38
            + x174 * x41,
            x16 * x169 * x76 + x169 * x175 * x74 + x171 * x74,
            x103 * x176 + x16 * x168 * x95 + x176 * x41 * x95,
            x123 * x16 * x163 + x163 * x177 * x41 + x165 * x177,
            x125 * x147 * x27 + x126 * x146 * x27 + x148 * x6 * x60,
            x128 * x146 * x51 + x128 * x147 * x20 + x131 * x146 * x20,
            x125 * x149 * x51 + x125 * x152 * x20 + x126 * x149 * x20,
            x133 * x146 * x45 + x133 * x147 * x43 + x136 * x170 * x36,
            x128 * x149 * x45 + x128 * x152 * x43 + x131 * x149 * x43,
            x11 * x157 * x178 + x125 * x154 * x45 + x126 * x154 * x43,
            x138 * x147 * x16 + x138 * x175 * x36 + x143 * x36,
            x136 * x149 * x16 + x149 * x179 * x41 + x152 * x179,
            x128 * x157 * x16 + x128 * x180 * x41 + x131 * x180,
            x126 * x159 * x16 + x159 * x178 * x41 + x162 * x6,
            x181 * x182 * x9 + x181 * x60 * x65 + x182 * x183,
            x181 * x20 * x67 + x181 * x51 * x66 + x183 * x20 * x66,
            x185 * x20 * x9 + x185 * x51 + x186 * x20 * x65,
            x181 * x43 * x88 + x181 * x45 * x79 + x183 * x43 * x79,
            x184 * x43 * x67 + x184 * x45 * x66 + x186 * x187 * x62,
            x11 * x189 + x187 * x188 * x9 + x188 * x45 * x65,
            x114 * x16 * x181 + x181 * x190 * x41 + x183 * x190,
            x16 * x186 * x79 + x191 * x41 * x79 + x191 * x88,
            x144 * x188 * x41 * x62 + x16 * x188 * x67 + x189 * x62,
            x144
            * (
                x161 * x36
                + x3 * (4 * x121 + x122 + x154 * x57 + 3 * x156)
                + x39 * (2 * a * x193 - 2 * x117 - 4 * x192)
            )
            + x194 * x41
            + x194 * x9,
        ]
    )
    return S


def kinetic3d_24(a, A, b, B):
    """Cartesian 3D (dg) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = -x5 - A[1]
    x7 = x6**2
    x8 = 2 * a**2
    x9 = -x0 - x8 * (x3 + x7)
    x10 = a * x4
    x11 = b * x10
    x12 = numpy.exp(-x11 * (A[0] - B[0]) ** 2)
    x13 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x14 = x12 * x13
    x15 = x14 * x3
    x16 = 3 * x15
    x17 = -x4 * (a * A[0] + b * B[0])
    x18 = -x17 - B[0]
    x19 = x14 * x18**2
    x20 = x3 * (x16 + 3 * x19)
    x21 = 2 * x15
    x22 = x18 * x21
    x23 = x15 + x19
    x24 = x18 * x23
    x25 = x22 + x24
    x26 = x18 * x25
    x27 = -x17 - A[0]
    x28 = x25 * x27
    x29 = 8 * x15 * x18
    x30 = x20 + x26
    x31 = x27 * x30 + x3 * (4 * x24 + x29)
    x32 = x27 * x31 + x3 * (5 * x20 + x26 + 4 * x28)
    x33 = numpy.exp(-x11 * (A[1] - B[1]) ** 2)
    x34 = numpy.exp(-x11 * (A[2] - B[2]) ** 2)
    x35 = numpy.pi * x34 * x4
    x36 = x33 * x35
    x37 = x32 * x36
    x38 = -x4 * (a * A[2] + b * B[2])
    x39 = -x38 - A[2]
    x40 = x39**2
    x41 = -x0 - x8 * (x3 + x40)
    x42 = b * x4
    x43 = x14 * x18
    x44 = 4 * x11
    x45 = x27**2
    x46 = -x0 - x8 * (x3 + x45)
    x47 = 2 * x14
    x48 = x18 * x47
    x49 = x3 * (x43 * x44 + x46 * x48)
    x50 = x15 * x46
    x51 = x1 * x42
    x52 = x43 * x46 + x43 * x51
    x53 = x18 * x52
    x54 = -x47
    x55 = x10 * (x2 * x23 + x54)
    x56 = x53 + x55
    x57 = x50 + x56
    x58 = x18 * x57
    x59 = x10 * (x2 * x25 - 3 * x43)
    x60 = x58 + x59
    x61 = x49 + x60
    x62 = x27 * x61
    x63 = 3 * x50
    x64 = x3 * (3 * x53 + 3 * x55 + x63)
    x65 = x20 + x28
    x66 = 8 * x11
    x67 = 4 * x15
    x68 = x10 * (2 * b * x30 - 4 * x19 - x67) + x18 * x61
    x69 = 4 * x49
    x70 = x64 + x68
    x71 = x27 * x70 + x3 * (4 * x58 + 4 * x59 + x69) + x31 * x51
    x72 = -x5 - B[1]
    x73 = x23 * x27
    x74 = x27 * x65 + x3 * (x24 + x29 + 3 * x73)
    x75 = x36 * x74
    x76 = x13 * x33
    x77 = x72 * x76
    x78 = x51 * x77 + x77 * x9
    x79 = x13 * x34
    x80 = x27 * x57
    x81 = x22 + x73
    x82 = 6 * x11
    x83 = x51 * x65 + x62 + x64
    x84 = x36 * (
        x27 * x83
        + x3 * (x60 + x69 + 3 * x80 + x81 * x82)
        + x42 * (2 * a * x74 - x18 * x67 - 2 * x24)
    )
    x85 = -x38 - B[2]
    x86 = x79 * x85
    x87 = x41 * x86 + x51 * x86
    x88 = x3 * x76
    x89 = x72**2 * x76
    x90 = x88 + x89
    x91 = x27 * x81 + x3 * (x16 + x19 + x27 * x48)
    x92 = x79 * x91
    x93 = x88 * x9
    x94 = x72 * x78
    x95 = 2 * x76
    x96 = -x95
    x97 = x10 * (x2 * x90 + x96)
    x98 = x94 + x97
    x99 = x93 + x98
    x100 = x27 * x52
    x101 = x14 * x27
    x102 = x101 * x18 + x15
    x103 = x49 + x51 * x81 + x80
    x104 = (
        x103 * x27
        + x3 * (2 * x100 + x102 * x44 + x56 + x63)
        + x42 * (2 * a * x91 - 2 * x19 - x21)
    )
    x105 = x36 * x85
    x106 = x3 * x79
    x107 = x79 * x85**2
    x108 = x106 + x107
    x109 = x76 * x91
    x110 = x106 * x41
    x111 = x85 * x87
    x112 = 2 * x79
    x113 = -x112
    x114 = x10 * (x108 * x2 + x113)
    x115 = x111 + x114
    x116 = x110 + x115
    x117 = 2 * x88
    x118 = x117 * x72
    x119 = x72 * x90
    x120 = x118 + x119
    x121 = x102 * x27 + x3 * (x101 + x43)
    x122 = x121 * x79
    x123 = x72 * x95
    x124 = x3 * (x123 * x9 + x44 * x77)
    x125 = x72 * x99
    x126 = x10 * (x120 * x2 - 3 * x77)
    x127 = x125 + x126
    x128 = x124 + x127
    x129 = x101 * x46 + x101 * x51
    x130 = x100 + x102 * x51 + x50
    x131 = x130 * x27 + x3 * (x129 + x52) + x42 * (x1 * x121 - x48)
    x132 = 2 * x106
    x133 = x132 * x85
    x134 = x108 * x85
    x135 = x133 + x134
    x136 = x121 * x76
    x137 = x112 * x85
    x138 = x3 * (x137 * x41 + x44 * x86)
    x139 = x116 * x85
    x140 = x10 * (x135 * x2 - 3 * x86)
    x141 = x139 + x140
    x142 = x138 + x141
    x143 = x14 * x45 + x15
    x144 = 3 * x88
    x145 = x3 * (x144 + 3 * x89)
    x146 = x120 * x72
    x147 = x145 + x146
    x148 = x147 * x79
    x149 = x129 * x27 + x42 * (x1 * x143 + x54) + x50
    x150 = 3 * x93
    x151 = x3 * (x150 + 3 * x94 + 3 * x97)
    x152 = 4 * x88
    x153 = x10 * (2 * b * x147 - x152 - 4 * x89) + x128 * x72
    x154 = x151 + x153
    x155 = 3 * x106
    x156 = x3 * (3 * x107 + x155)
    x157 = x135 * x85
    x158 = x156 + x157
    x159 = x158 * x76
    x160 = 3 * x110
    x161 = x3 * (3 * x111 + 3 * x114 + x160)
    x162 = 4 * x106
    x163 = x10 * (2 * b * x158 - 4 * x107 - x162) + x142 * x85
    x164 = x161 + x163
    x165 = x31 * x36
    x166 = x6 * x76
    x167 = x166 * x51 + x166 * x9
    x168 = x36 * x71
    x169 = x166 * x72 + x88
    x170 = x65 * x79
    x171 = x6 * x78
    x172 = x169 * x51 + x171 + x93
    x173 = x6 * x90
    x174 = x118 + x173
    x175 = x79 * x81
    x176 = x6 * x99
    x177 = x124 + x174 * x51 + x176
    x178 = x120 * x6
    x179 = x145 + x178
    x180 = x179 * x79
    x181 = x128 * x6
    x182 = x151 + x179 * x51 + x181
    x183 = 8 * x72 * x88
    x184 = x147 * x6 + x3 * (4 * x119 + x183)
    x185 = x12 * x35
    x186 = x185 * x27
    x187 = 4 * x124
    x188 = x154 * x6 + x184 * x51 + x3 * (4 * x125 + 4 * x126 + x187)
    x189 = x185 * x188
    x190 = numpy.pi * x12 * x33 * x4
    x191 = x190 * x27
    x192 = x39 * x79
    x193 = x192 * x41 + x192 * x51
    x194 = x36 * x39
    x195 = x106 + x192 * x85
    x196 = x65 * x76
    x197 = x39 * x87
    x198 = x110 + x195 * x51 + x197
    x199 = x108 * x39
    x200 = x133 + x199
    x201 = x76 * x81
    x202 = x116 * x39
    x203 = x138 + x200 * x51 + x202
    x204 = x135 * x39
    x205 = x156 + x204
    x206 = x205 * x76
    x207 = x142 * x39
    x208 = x161 + x205 * x51 + x207
    x209 = 8 * x106 * x85
    x210 = x158 * x39 + x3 * (4 * x134 + x209)
    x211 = 4 * x138
    x212 = x164 * x39 + x210 * x51 + x3 * (4 * x139 + 4 * x140 + x211)
    x213 = x190 * x212
    x214 = x7 * x76 + x88
    x215 = x30 * x79
    x216 = x167 * x6 + x42 * (x1 * x214 + x96) + x93
    x217 = x169 * x6 + x3 * (x166 + x77)
    x218 = x217 * x79
    x219 = x172 * x6 + x3 * (x167 + x78) + x42 * (x1 * x217 - x123)
    x220 = x174 * x6 + x3 * (x123 * x6 + x144 + x89)
    x221 = x220 * x79
    x222 = (
        x177 * x6
        + x3 * (x150 + x169 * x44 + 2 * x171 + x98)
        + x42 * (2 * a * x220 - x117 - 2 * x89)
    )
    x223 = x179 * x6 + x3 * (x119 + 3 * x173 + x183)
    x224 = x18 * x185
    x225 = x185 * (
        x182 * x6
        + x3 * (x127 + x174 * x82 + 3 * x176 + x187)
        + x42 * (2 * a * x223 - 2 * x119 - x152 * x72)
    )
    x226 = x184 * x6 + x3 * (5 * x145 + x146 + 4 * x178)
    x227 = x185 * x226
    x228 = x185 * x46
    x229 = x14 * x220
    x230 = x14 * x217
    x231 = x14 * x158
    x232 = x190 * x6
    x233 = x14 * x179
    x234 = x14 * x174
    x235 = x14 * x205
    x236 = x106 + x40 * x79
    x237 = x30 * x76
    x238 = x110 + x193 * x39 + x42 * (x1 * x236 + x113)
    x239 = x195 * x39 + x3 * (x192 + x86)
    x240 = x239 * x76
    x241 = x198 * x39 + x3 * (x193 + x87) + x42 * (x1 * x239 - x137)
    x242 = x200 * x39 + x3 * (x107 + x137 * x39 + x155)
    x243 = x242 * x76
    x244 = (
        x203 * x39
        + x3 * (x115 + x160 + x195 * x44 + 2 * x197)
        + x42 * (2 * a * x242 - 2 * x107 - x132)
    )
    x245 = x18 * x190
    x246 = x205 * x39 + x3 * (x134 + 3 * x199 + x209)
    x247 = x190 * (
        x208 * x39
        + x3 * (x141 + x200 * x82 + 3 * x202 + x211)
        + x42 * (2 * a * x246 - 2 * x134 - x162 * x85)
    )
    x248 = x14 * x147
    x249 = x14 * x239
    x250 = x14 * x242
    x251 = x210 * x39 + x3 * (5 * x156 + x157 + 4 * x204)
    x252 = x190 * x251

    # 90 item(s)
    S = numpy.array(
        [
            x36
            * (
                x27 * x71
                + x3 * (4 * x62 + 5 * x64 + x65 * x66 + x68)
                + x42 * (2 * a * x32 - 2 * x20 - 2 * x26)
            )
            + x37 * x41
            + x37 * x9,
            x41 * x72 * x75 + x72 * x84 + x74 * x78 * x79,
            x74 * x76 * x87 + x75 * x85 * x9 + x84 * x85,
            x104 * x79 * x90 + x41 * x90 * x92 + x92 * x99,
            x104 * x105 * x72 + x77 * x87 * x91 + x78 * x86 * x91,
            x104 * x108 * x76 + x108 * x109 * x9 + x109 * x116,
            x120 * x122 * x41 + x120 * x131 * x79 + x122 * x128,
            x121 * x86 * x99 + x121 * x87 * x90 + x131 * x86 * x90,
            x108 * x121 * x78 + x108 * x131 * x77 + x116 * x121 * x77,
            x131 * x135 * x76 + x135 * x136 * x9 + x136 * x142,
            x143 * x148 * x41 + x143 * x154 * x79 + x148 * x149,
            x120 * x143 * x87 + x120 * x149 * x86 + x128 * x143 * x86,
            x108 * x143 * x99 + x108 * x149 * x90 + x116 * x143 * x90,
            x135 * x143 * x78 + x135 * x149 * x77 + x142 * x143 * x77,
            x143 * x159 * x9 + x143 * x164 * x76 + x149 * x159,
            x165 * x41 * x6 + x167 * x31 * x79 + x168 * x6,
            x169 * x170 * x41 + x169 * x79 * x83 + x170 * x172,
            x105 * x6 * x83 + x166 * x65 * x87 + x167 * x65 * x86,
            x103 * x174 * x79 + x174 * x175 * x41 + x175 * x177,
            x103 * x169 * x86 + x169 * x81 * x87 + x172 * x81 * x86,
            x103 * x108 * x166 + x108 * x167 * x81 + x116 * x166 * x81,
            x102 * x180 * x41 + x102 * x182 * x79 + x130 * x180,
            x102 * x174 * x87 + x102 * x177 * x86 + x130 * x174 * x86,
            x102 * x108 * x172 + x102 * x116 * x169 + x108 * x130 * x169,
            x102 * x135 * x167 + x102 * x142 * x166 + x130 * x135 * x166,
            x129 * x184 * x79 + x184 * x186 * x41 + x189 * x27,
            x101 * x179 * x87 + x129 * x179 * x86 + x182 * x186 * x85,
            x101 * x108 * x177 + x101 * x116 * x174 + x108 * x129 * x174,
            x101 * x135 * x172 + x101 * x142 * x169 + x129 * x135 * x169,
            x101 * x158 * x167 + x129 * x158 * x166 + x164 * x191 * x6,
            x165 * x39 * x9 + x168 * x39 + x193 * x31 * x76,
            x192 * x65 * x78 + x193 * x65 * x77 + x194 * x72 * x83,
            x195 * x196 * x9 + x195 * x76 * x83 + x196 * x198,
            x103 * x192 * x90 + x192 * x81 * x99 + x193 * x81 * x90,
            x103 * x195 * x77 + x195 * x78 * x81 + x198 * x77 * x81,
            x103 * x200 * x76 + x200 * x201 * x9 + x201 * x203,
            x102 * x120 * x193 + x102 * x128 * x192 + x120 * x130 * x192,
            x102 * x195 * x99 + x102 * x198 * x90 + x130 * x195 * x90,
            x102 * x200 * x78 + x102 * x203 * x77 + x130 * x200 * x77,
            x102 * x206 * x9 + x102 * x208 * x76 + x130 * x206,
            x101 * x147 * x193 + x129 * x147 * x192 + x154 * x186 * x39,
            x101 * x120 * x198 + x101 * x128 * x195 + x120 * x129 * x195,
            x101 * x200 * x99 + x101 * x203 * x90 + x129 * x200 * x90,
            x101 * x205 * x78 + x129 * x205 * x77 + x191 * x208 * x72,
            x129 * x210 * x76 + x191 * x210 * x9 + x213 * x27,
            x214 * x215 * x41 + x214 * x70 * x79 + x215 * x216,
            x218 * x25 * x41 + x218 * x61 + x219 * x25 * x79,
            x214 * x25 * x87 + x214 * x61 * x86 + x216 * x25 * x86,
            x221 * x23 * x41 + x221 * x57 + x222 * x23 * x79,
            x217 * x23 * x87 + x217 * x57 * x86 + x219 * x23 * x86,
            x108 * x214 * x57 + x108 * x216 * x23 + x116 * x214 * x23,
            x18 * x225 + x223 * x224 * x41 + x223 * x52 * x79,
            x220 * x43 * x87 + x220 * x52 * x86 + x222 * x224 * x85,
            x108 * x217 * x52 + x108 * x219 * x43 + x116 * x217 * x43,
            x135 * x214 * x52 + x135 * x216 * x43 + x142 * x214 * x43,
            x185
            * (
                x188 * x6
                + x3 * (5 * x151 + x153 + x179 * x66 + 4 * x181)
                + x42 * (2 * a * x226 - 2 * x145 - 2 * x146)
            )
            + x227 * x41
            + x227 * x46,
            x14 * x223 * x87 + x223 * x228 * x85 + x225 * x85,
            x108 * x14 * x222 + x108 * x229 * x46 + x116 * x229,
            x135 * x14 * x219 + x135 * x230 * x46 + x142 * x230,
            x14 * x164 * x214 + x214 * x231 * x46 + x216 * x231,
            x166 * x193 * x30 + x167 * x192 * x30 + x194 * x6 * x70,
            x169 * x192 * x61 + x169 * x193 * x25 + x172 * x192 * x25,
            x166 * x195 * x61 + x166 * x198 * x25 + x167 * x195 * x25,
            x174 * x192 * x57 + x174 * x193 * x23 + x177 * x192 * x23,
            x169 * x195 * x57 + x169 * x198 * x23 + x172 * x195 * x23,
            x166 * x200 * x57 + x166 * x203 * x23 + x167 * x200 * x23,
            x179 * x192 * x52 + x179 * x193 * x43 + x182 * x224 * x39,
            x174 * x195 * x52 + x174 * x198 * x43 + x177 * x195 * x43,
            x169 * x200 * x52 + x169 * x203 * x43 + x172 * x200 * x43,
            x166 * x205 * x52 + x167 * x205 * x43 + x18 * x208 * x232,
            x14 * x184 * x193 + x184 * x228 * x39 + x189 * x39,
            x14 * x182 * x195 + x195 * x233 * x46 + x198 * x233,
            x14 * x177 * x200 + x200 * x234 * x46 + x203 * x234,
            x14 * x169 * x208 + x169 * x235 * x46 + x172 * x235,
            x14 * x167 * x210 + x210 * x232 * x46 + x213 * x6,
            x236 * x237 * x9 + x236 * x70 * x76 + x237 * x238,
            x236 * x25 * x78 + x236 * x61 * x77 + x238 * x25 * x77,
            x240 * x25 * x9 + x240 * x61 + x241 * x25 * x76,
            x23 * x236 * x99 + x23 * x238 * x90 + x236 * x57 * x90,
            x23 * x239 * x78 + x23 * x241 * x77 + x239 * x57 * x77,
            x23 * x243 * x9 + x23 * x244 * x76 + x243 * x57,
            x120 * x236 * x52 + x120 * x238 * x43 + x128 * x236 * x43,
            x239 * x43 * x99 + x239 * x52 * x90 + x241 * x43 * x90,
            x242 * x43 * x78 + x242 * x52 * x77 + x244 * x245 * x72,
            x18 * x247 + x245 * x246 * x9 + x246 * x52 * x76,
            x14 * x154 * x236 + x236 * x248 * x46 + x238 * x248,
            x120 * x14 * x241 + x120 * x249 * x46 + x128 * x249,
            x14 * x244 * x90 + x250 * x46 * x90 + x250 * x99,
            x14 * x246 * x78 + x190 * x246 * x46 * x72 + x247 * x72,
            x190
            * (
                x212 * x39
                + x3 * (5 * x161 + x163 + x205 * x66 + 4 * x207)
                + x42 * (2 * a * x251 - 2 * x156 - 2 * x157)
            )
            + x252 * x46
            + x252 * x9,
        ]
    )
    return S


def kinetic3d_30(a, A, b, B):
    """Cartesian 3D (fs) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = (2 * b + x1) ** (-1.0)
    x3 = (a + b) ** (-1.0)
    x4 = x3 * (a * A[1] + b * B[1]) - A[1]
    x5 = x4**2
    x6 = 2 * a**2
    x7 = -x0 - x6 * (x2 + x5)
    x8 = x3 * (a * A[0] + b * B[0]) - A[0]
    x9 = b * x3
    x10 = a * x9
    x11 = numpy.exp(-x10 * (A[0] - B[0]) ** 2)
    x12 = numpy.sqrt(numpy.pi) * numpy.sqrt(x3)
    x13 = x11 * x12
    x14 = x13 * x2
    x15 = x8**2
    x16 = x13 * x15 + x14
    x17 = 2 * x14 * x8 + x16 * x8
    x18 = numpy.exp(-x10 * (A[1] - B[1]) ** 2)
    x19 = numpy.exp(-x10 * (A[2] - B[2]) ** 2)
    x20 = numpy.pi * x19 * x3
    x21 = x18 * x20
    x22 = x17 * x21
    x23 = x3 * (a * A[2] + b * B[2]) - A[2]
    x24 = x23**2
    x25 = -x0 - x6 * (x2 + x24)
    x26 = x13 * x8
    x27 = 4 * x10
    x28 = -x0 - x6 * (x15 + x2)
    x29 = 2 * x13
    x30 = x1 * x9
    x31 = x26 * x28 + x26 * x30
    x32 = x14 * x28 + x31 * x8 + x9 * (x1 * x16 - x29)
    x33 = x16 * x21
    x34 = x12 * x18
    x35 = x34 * x4
    x36 = x30 * x35 + x35 * x7
    x37 = x12 * x19
    x38 = x21 * x32
    x39 = x23 * x37
    x40 = x25 * x39 + x30 * x39
    x41 = x2 * x34
    x42 = x34 * x5 + x41
    x43 = x11 * x20
    x44 = x43 * x8
    x45 = 2 * x34
    x46 = x36 * x4 + x41 * x7 + x9 * (x1 * x42 - x45)
    x47 = x43 * x46
    x48 = numpy.pi * x11 * x18 * x3
    x49 = x48 * x8
    x50 = x2 * x37
    x51 = x24 * x37 + x50
    x52 = 2 * x37
    x53 = x23 * x40 + x25 * x50 + x9 * (x1 * x51 - x52)
    x54 = x48 * x53
    x55 = 2 * x4 * x41 + x4 * x42
    x56 = x43 * x55
    x57 = 2 * x23 * x50 + x23 * x51
    x58 = x48 * x57

    # 10 item(s)
    S = numpy.array(
        [
            x21
            * (x2 * (x26 * x27 + x28 * x29 * x8) + x32 * x8 + x9 * (x1 * x17 - 3 * x26))
            + x22 * x25
            + x22 * x7,
            x16 * x36 * x37 + x25 * x33 * x4 + x38 * x4,
            x16 * x34 * x40 + x23 * x33 * x7 + x23 * x38,
            x25 * x42 * x44 + x31 * x37 * x42 + x47 * x8,
            x21 * x23 * x31 * x4 + x23 * x36 * x44 + x4 * x40 * x49,
            x31 * x34 * x51 + x49 * x51 * x7 + x54 * x8,
            x25 * x56
            + x28 * x56
            + x43
            * (x2 * (x27 * x35 + x4 * x45 * x7) + x4 * x46 + x9 * (x1 * x55 - 3 * x35)),
            x13 * x40 * x42 + x23 * x28 * x42 * x43 + x23 * x47,
            x13 * x36 * x51 + x28 * x4 * x48 * x51 + x4 * x54,
            x28 * x58
            + x48
            * (x2 * (x23 * x25 * x52 + x27 * x39) + x23 * x53 + x9 * (x1 * x57 - 3 * x39))
            + x58 * x7,
        ]
    )
    return S


def kinetic3d_31(a, A, b, B):
    """Cartesian 3D (fp) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = (2 * b + x1) ** (-1.0)
    x3 = (a + b) ** (-1.0)
    x4 = -x3 * (a * A[1] + b * B[1])
    x5 = -x4 - A[1]
    x6 = x5**2
    x7 = 2 * a**2
    x8 = -x0 - x7 * (x2 + x6)
    x9 = b * x3
    x10 = a * x9
    x11 = numpy.exp(-x10 * (A[0] - B[0]) ** 2)
    x12 = numpy.sqrt(numpy.pi) * numpy.sqrt(x3)
    x13 = x11 * x12
    x14 = x13 * x2
    x15 = 3 * x14
    x16 = -x3 * (a * A[0] + b * B[0])
    x17 = -x16 - A[0]
    x18 = x17**2
    x19 = x13 * x18
    x20 = -x16 - B[0]
    x21 = 2 * x13
    x22 = x20 * x21
    x23 = x13 * x17
    x24 = x13 * x20
    x25 = x20 * x23
    x26 = x14 + x25
    x27 = x17 * x26 + x2 * (x23 + x24)
    x28 = x17 * x27 + x2 * (x15 + x17 * x22 + x19)
    x29 = numpy.exp(-x10 * (A[1] - B[1]) ** 2)
    x30 = numpy.exp(-x10 * (A[2] - B[2]) ** 2)
    x31 = numpy.pi * x3 * x30
    x32 = x29 * x31
    x33 = x28 * x32
    x34 = -x3 * (a * A[2] + b * B[2])
    x35 = -x34 - A[2]
    x36 = x35**2
    x37 = -x0 - x7 * (x2 + x36)
    x38 = x1 * x9
    x39 = -x0 - x7 * (x18 + x2)
    x40 = x24 * x38 + x24 * x39
    x41 = x17 * x40
    x42 = 4 * x10
    x43 = x14 * x39
    x44 = x23 * x38 + x23 * x39
    x45 = x14 + x19
    x46 = x17 * x44 + x9 * (x1 * x45 - x21)
    x47 = x26 * x38 + x41 + x43
    x48 = x17 * x47 + x2 * (x40 + x44) + x9 * (x1 * x27 - x22)
    x49 = -x4 - B[1]
    x50 = 2 * x14 * x17 + x17 * x45
    x51 = x32 * x50
    x52 = x12 * x29
    x53 = x49 * x52
    x54 = x38 * x53 + x53 * x8
    x55 = x12 * x30
    x56 = x43 + x46
    x57 = x32 * (
        x17 * x56 + x2 * (x17 * x21 * x39 + x23 * x42) + x9 * (x1 * x50 - 3 * x23)
    )
    x58 = -x34 - B[2]
    x59 = x55 * x58
    x60 = x37 * x59 + x38 * x59
    x61 = x32 * x5
    x62 = x5 * x52
    x63 = x38 * x62 + x62 * x8
    x64 = x32 * x48
    x65 = x12 * x2
    x66 = x29 * x65
    x67 = x49 * x62
    x68 = x66 + x67
    x69 = x45 * x55
    x70 = x66 * x8
    x71 = x5 * x54
    x72 = x38 * x68 + x70 + x71
    x73 = x32 * x35
    x74 = x35 * x55
    x75 = x37 * x74 + x38 * x74
    x76 = x30 * x65
    x77 = x58 * x74
    x78 = x76 + x77
    x79 = x45 * x52
    x80 = x37 * x76
    x81 = x35 * x60
    x82 = x38 * x78 + x80 + x81
    x83 = x52 * x6
    x84 = x66 + x83
    x85 = x55 * x84
    x86 = 2 * x52
    x87 = x5 * x63 + x9 * (x1 * x84 - x86)
    x88 = x70 + x87
    x89 = x2 * (x53 + x62) + x5 * x68
    x90 = x11 * x31
    x91 = x17 * x90
    x92 = x49 * x86
    x93 = x2 * (x54 + x63) + x5 * x72 + x9 * (x1 * x89 - x92)
    x94 = x90 * x93
    x95 = numpy.pi * x11 * x29 * x3
    x96 = x17 * x95
    x97 = x36 * x55
    x98 = x76 + x97
    x99 = x52 * x98
    x100 = 2 * x55
    x101 = x35 * x75 + x9 * (x1 * x98 - x100)
    x102 = x101 + x80
    x103 = x2 * (x59 + x74) + x35 * x78
    x104 = x100 * x58
    x105 = x2 * (x60 + x75) + x35 * x82 + x9 * (x1 * x103 - x104)
    x106 = x105 * x95
    x107 = 2 * x5 * x66 + x5 * x84
    x108 = x107 * x90
    x109 = x90 * (
        x2 * (x42 * x62 + x5 * x8 * x86) + x5 * x88 + x9 * (x1 * x107 - 3 * x62)
    )
    x110 = 3 * x66
    x111 = x2 * (x110 + x5 * x92 + x83) + x5 * x89
    x112 = x111 * x90
    x113 = x35 * x90
    x114 = x13 * x84
    x115 = x5 * x95
    x116 = x13 * x98
    x117 = 2 * x35 * x76 + x35 * x98
    x118 = x117 * x95
    x119 = x95 * (
        x102 * x35 + x2 * (x100 * x35 * x37 + x42 * x74) + x9 * (x1 * x117 - 3 * x74)
    )
    x120 = 3 * x76
    x121 = x103 * x35 + x2 * (x104 * x35 + x120 + x97)
    x122 = x121 * x95

    # 30 item(s)
    S = numpy.array(
        [
            x32
            * (
                x17 * x48
                + x2 * (x26 * x42 + 2 * x41 + 3 * x43 + x46)
                + x9 * (2 * a * x28 - x15 - 3 * x25)
            )
            + x33 * x37
            + x33 * x8,
            x37 * x49 * x51 + x49 * x57 + x50 * x54 * x55,
            x50 * x52 * x60 + x51 * x58 * x8 + x57 * x58,
            x27 * x37 * x61 + x27 * x55 * x63 + x5 * x64,
            x37 * x68 * x69 + x55 * x56 * x68 + x69 * x72,
            x45 * x59 * x63 + x45 * x60 * x62 + x56 * x58 * x61,
            x27 * x52 * x75 + x27 * x73 * x8 + x35 * x64,
            x45 * x53 * x75 + x45 * x54 * x74 + x49 * x56 * x73,
            x52 * x56 * x78 + x78 * x79 * x8 + x79 * x82,
            x26 * x37 * x85 + x26 * x55 * x88 + x47 * x85,
            x17 * x94 + x37 * x89 * x91 + x44 * x55 * x89,
            x23 * x60 * x84 + x44 * x59 * x84 + x58 * x88 * x91,
            x26 * x62 * x75 + x26 * x63 * x74 + x35 * x47 * x61,
            x23 * x68 * x75 + x35 * x72 * x91 + x44 * x68 * x74,
            x23 * x63 * x78 + x44 * x62 * x78 + x5 * x82 * x96,
            x102 * x26 * x52 + x26 * x8 * x99 + x47 * x99,
            x102 * x49 * x96 + x23 * x54 * x98 + x44 * x53 * x98,
            x103 * x44 * x52 + x103 * x8 * x96 + x106 * x17,
            x107 * x40 * x55 + x108 * x20 * x37 + x109 * x20,
            x112 * x37
            + x112 * x39
            + x90
            * (
                x2 * (x42 * x68 + 3 * x70 + 2 * x71 + x87)
                + x5 * x93
                + x9 * (2 * a * x111 - x110 - 3 * x67)
            ),
            x107 * x13 * x60 + x108 * x39 * x58 + x109 * x58,
            x113 * x20 * x88 + x24 * x75 * x84 + x40 * x74 * x84,
            x113 * x39 * x89 + x13 * x75 * x89 + x35 * x94,
            x114 * x39 * x78 + x114 * x82 + x13 * x78 * x88,
            x102 * x115 * x20 + x24 * x63 * x98 + x40 * x62 * x98,
            x102 * x13 * x68 + x116 * x39 * x68 + x116 * x72,
            x103 * x115 * x39 + x103 * x13 * x63 + x106 * x5,
            x117 * x40 * x52 + x118 * x20 * x8 + x119 * x20,
            x117 * x13 * x54 + x118 * x39 * x49 + x119 * x49,
            x122 * x39
            + x122 * x8
            + x95
            * (
                x105 * x35
                + x2 * (x101 + x42 * x78 + 3 * x80 + 2 * x81)
                + x9 * (2 * a * x121 - x120 - 3 * x77)
            ),
        ]
    )
    return S


def kinetic3d_32(a, A, b, B):
    """Cartesian 3D (fd) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = -x5 - A[1]
    x7 = x6**2
    x8 = 2 * a**2
    x9 = -x0 - x8 * (x3 + x7)
    x10 = -x4 * (a * A[0] + b * B[0])
    x11 = -x10 - A[0]
    x12 = -x10 - B[0]
    x13 = a * x4
    x14 = b * x13
    x15 = numpy.exp(-x14 * (A[0] - B[0]) ** 2)
    x16 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x17 = x15 * x16
    x18 = x12**2 * x17
    x19 = x17 * x3
    x20 = 3 * x19
    x21 = 2 * x17
    x22 = x12 * x21
    x23 = x11 * x22 + x20
    x24 = 2 * x19
    x25 = x18 + x19
    x26 = x11 * x25
    x27 = x12 * x24 + x26
    x28 = x11 * x27 + x3 * (x18 + x23)
    x29 = x12 * x19
    x30 = x11 * x17
    x31 = x12 * x17
    x32 = x3 * (x30 + x31)
    x33 = x12 * x30
    x34 = x19 + x33
    x35 = x11 * x34
    x36 = x11 * x28 + x3 * (2 * x26 + 4 * x29 + 2 * x32 + 2 * x35)
    x37 = numpy.exp(-x14 * (A[1] - B[1]) ** 2)
    x38 = numpy.exp(-x14 * (A[2] - B[2]) ** 2)
    x39 = numpy.pi * x38 * x4
    x40 = x37 * x39
    x41 = x36 * x40
    x42 = -x4 * (a * A[2] + b * B[2])
    x43 = -x42 - A[2]
    x44 = x43**2
    x45 = -x0 - x8 * (x3 + x44)
    x46 = b * x4
    x47 = 4 * x14
    x48 = x11**2
    x49 = -x0 - x8 * (x3 + x48)
    x50 = x3 * (x22 * x49 + x31 * x47)
    x51 = x32 + x35
    x52 = x1 * x51 - x22
    x53 = x2 * x4
    x54 = x1 * x46
    x55 = x31 * x49 + x31 * x54
    x56 = x30 * x49 + x30 * x54
    x57 = x3 * (x55 + x56)
    x58 = x19 * x49
    x59 = x11 * x55
    x60 = x34 * x54 + x58 + x59
    x61 = x11 * x60
    x62 = -x21
    x63 = x12 * x55 + x13 * (x2 * x25 + x62)
    x64 = x58 + x63
    x65 = x11 * x64
    x66 = x34 * x47 + 3 * x58 + 2 * x59
    x67 = x27 * x54 + x50 + x65
    x68 = x11 * x67 + x3 * (x63 + x66) + x46 * (2 * a * x28 - 2 * x18 - x24)
    x69 = -x5 - B[1]
    x70 = x17 * x48
    x71 = x11 * x51 + x3 * (x23 + x70)
    x72 = x40 * x71
    x73 = x16 * x37
    x74 = x69 * x73
    x75 = x54 * x74 + x74 * x9
    x76 = x16 * x38
    x77 = x19 + x70
    x78 = x11 * x56 + x46 * (x1 * x77 + x62)
    x79 = x46 * x52 + x57 + x61
    x80 = x40 * (x11 * x79 + x3 * (x66 + x78) + x46 * (2 * a * x71 - x20 - 3 * x33))
    x81 = -x42 - B[2]
    x82 = x76 * x81
    x83 = x45 * x82 + x54 * x82
    x84 = x3 * x73
    x85 = x69**2 * x73
    x86 = x84 + x85
    x87 = x11 * x24 + x11 * x77
    x88 = x76 * x87
    x89 = x84 * x9
    x90 = 2 * x73
    x91 = -x90
    x92 = x13 * (x2 * x86 + x91) + x69 * x75
    x93 = x89 + x92
    x94 = x58 + x78
    x95 = x11 * x94 + x3 * (x11 * x21 * x49 + x30 * x47) + x46 * (x1 * x87 - 3 * x30)
    x96 = x40 * x81
    x97 = x3 * x76
    x98 = x76 * x81**2
    x99 = x97 + x98
    x100 = x73 * x87
    x101 = x45 * x97
    x102 = 2 * x76
    x103 = -x102
    x104 = x13 * (x103 + x2 * x99) + x81 * x83
    x105 = x101 + x104
    x106 = x28 * x40
    x107 = x6 * x73
    x108 = x107 * x54 + x107 * x9
    x109 = x40 * x68
    x110 = x107 * x69
    x111 = x110 + x84
    x112 = x51 * x76
    x113 = x6 * x75
    x114 = x111 * x54 + x113 + x89
    x115 = 2 * x84
    x116 = x6 * x86
    x117 = x115 * x69 + x116
    x118 = x117 * x76
    x119 = x69 * x90
    x120 = x3 * (x119 * x9 + x47 * x74)
    x121 = x6 * x93
    x122 = x117 * x54 + x120 + x121
    x123 = x43 * x76
    x124 = x123 * x45 + x123 * x54
    x125 = x40 * x43
    x126 = x123 * x81
    x127 = x126 + x97
    x128 = x51 * x73
    x129 = x43 * x83
    x130 = x101 + x127 * x54 + x129
    x131 = 2 * x97
    x132 = x43 * x99
    x133 = x131 * x81 + x132
    x134 = x133 * x73
    x135 = x102 * x81
    x136 = x3 * (x135 * x45 + x47 * x82)
    x137 = x105 * x43
    x138 = x133 * x54 + x136 + x137
    x139 = x7 * x73
    x140 = x139 + x84
    x141 = x27 * x76
    x142 = x108 * x6 + x46 * (x1 * x140 + x91)
    x143 = x142 + x89
    x144 = x3 * (x107 + x74)
    x145 = x111 * x6
    x146 = x144 + x145
    x147 = x146 * x76
    x148 = x1 * x146 - x119
    x149 = x3 * (x108 + x75)
    x150 = x114 * x6
    x151 = x148 * x46 + x149 + x150
    x152 = 3 * x84
    x153 = x119 * x6 + x152
    x154 = x117 * x6 + x3 * (x153 + x85)
    x155 = x15 * x39
    x156 = x11 * x155
    x157 = x111 * x47 + 2 * x113 + 3 * x89
    x158 = x122 * x6 + x3 * (x157 + x92) + x46 * (2 * a * x154 - x115 - 2 * x85)
    x159 = x155 * x158
    x160 = numpy.pi * x15 * x37 * x4
    x161 = x11 * x160
    x162 = x44 * x76
    x163 = x162 + x97
    x164 = x27 * x73
    x165 = x124 * x43 + x46 * (x1 * x163 + x103)
    x166 = x101 + x165
    x167 = x3 * (x123 + x82)
    x168 = x127 * x43
    x169 = x167 + x168
    x170 = x169 * x73
    x171 = x1 * x169 - x135
    x172 = x3 * (x124 + x83)
    x173 = x130 * x43
    x174 = x171 * x46 + x172 + x173
    x175 = 3 * x97
    x176 = x135 * x43 + x175
    x177 = x133 * x43 + x3 * (x176 + x98)
    x178 = 3 * x101 + x127 * x47 + 2 * x129
    x179 = x138 * x43 + x3 * (x104 + x178) + x46 * (2 * a * x177 - x131 - 2 * x98)
    x180 = x160 * x179
    x181 = x115 * x6 + x140 * x6
    x182 = x181 * x76
    x183 = x143 * x6 + x3 * (x107 * x47 + x6 * x9 * x90) + x46 * (x1 * x181 - 3 * x107)
    x184 = x146 * x6 + x3 * (x139 + x153)
    x185 = x12 * x155
    x186 = x155 * (
        x151 * x6 + x3 * (x142 + x157) + x46 * (2 * a * x184 - 3 * x110 - x152)
    )
    x187 = x69 * x84
    x188 = x154 * x6 + x3 * (2 * x116 + 2 * x144 + 2 * x145 + 4 * x187)
    x189 = x155 * x188
    x190 = x155 * x49
    x191 = x17 * x181
    x192 = x146 * x17
    x193 = x133 * x17
    x194 = x160 * x6
    x195 = x117 * x17
    x196 = x169 * x17
    x197 = x131 * x43 + x163 * x43
    x198 = x197 * x73
    x199 = (
        x166 * x43 + x3 * (x102 * x43 * x45 + x123 * x47) + x46 * (x1 * x197 - 3 * x123)
    )
    x200 = x12 * x160
    x201 = x169 * x43 + x3 * (x162 + x176)
    x202 = x160 * (
        x174 * x43 + x3 * (x165 + x178) + x46 * (2 * a * x201 - 3 * x126 - x175)
    )
    x203 = x17 * x197
    x204 = x81 * x97
    x205 = x177 * x43 + x3 * (2 * x132 + 2 * x167 + 2 * x168 + 4 * x204)
    x206 = x160 * x205

    # 60 item(s)
    S = numpy.array(
        [
            x40
            * (
                x11 * x68
                + x3 * (x27 * x47 + 2 * x50 + x52 * x53 + 2 * x57 + 2 * x61 + 2 * x65)
                + x46 * (2 * a * x36 - 3 * x26 - 6 * x29)
            )
            + x41 * x45
            + x41 * x9,
            x45 * x69 * x72 + x69 * x80 + x71 * x75 * x76,
            x71 * x73 * x83 + x72 * x81 * x9 + x80 * x81,
            x45 * x86 * x88 + x76 * x86 * x95 + x88 * x93,
            x69 * x95 * x96 + x74 * x83 * x87 + x75 * x82 * x87,
            x100 * x105 + x100 * x9 * x99 + x73 * x95 * x99,
            x106 * x45 * x6 + x108 * x28 * x76 + x109 * x6,
            x111 * x112 * x45 + x111 * x76 * x79 + x112 * x114,
            x107 * x51 * x83 + x108 * x51 * x82 + x6 * x79 * x96,
            x118 * x45 * x77 + x118 * x94 + x122 * x76 * x77,
            x111 * x77 * x83 + x111 * x82 * x94 + x114 * x77 * x82,
            x105 * x107 * x77 + x107 * x94 * x99 + x108 * x77 * x99,
            x106 * x43 * x9 + x109 * x43 + x124 * x28 * x73,
            x123 * x51 * x75 + x124 * x51 * x74 + x125 * x69 * x79,
            x127 * x128 * x9 + x127 * x73 * x79 + x128 * x130,
            x123 * x77 * x93 + x123 * x86 * x94 + x124 * x77 * x86,
            x127 * x74 * x94 + x127 * x75 * x77 + x130 * x74 * x77,
            x134 * x77 * x9 + x134 * x94 + x138 * x73 * x77,
            x140 * x141 * x45 + x140 * x67 * x76 + x141 * x143,
            x147 * x34 * x45 + x147 * x60 + x151 * x34 * x76,
            x140 * x34 * x83 + x140 * x60 * x82 + x143 * x34 * x82,
            x11 * x159 + x154 * x156 * x45 + x154 * x56 * x76,
            x146 * x30 * x83 + x146 * x56 * x82 + x151 * x156 * x81,
            x105 * x140 * x30 + x140 * x56 * x99 + x143 * x30 * x99,
            x107 * x124 * x27 + x108 * x123 * x27 + x125 * x6 * x67,
            x111 * x123 * x60 + x111 * x124 * x34 + x114 * x123 * x34,
            x107 * x127 * x60 + x107 * x130 * x34 + x108 * x127 * x34,
            x117 * x123 * x56 + x117 * x124 * x30 + x122 * x156 * x43,
            x111 * x127 * x56 + x111 * x130 * x30 + x114 * x127 * x30,
            x107 * x133 * x56 + x108 * x133 * x30 + x138 * x161 * x6,
            x163 * x164 * x9 + x163 * x67 * x73 + x164 * x166,
            x163 * x34 * x75 + x163 * x60 * x74 + x166 * x34 * x74,
            x170 * x34 * x9 + x170 * x60 + x174 * x34 * x73,
            x163 * x30 * x93 + x163 * x56 * x86 + x166 * x30 * x86,
            x161 * x174 * x69 + x169 * x30 * x75 + x169 * x56 * x74,
            x11 * x180 + x161 * x177 * x9 + x177 * x56 * x73,
            x182 * x25 * x45 + x182 * x64 + x183 * x25 * x76,
            x12 * x186 + x184 * x185 * x45 + x184 * x55 * x76,
            x181 * x31 * x83 + x181 * x55 * x82 + x183 * x185 * x81,
            x155
            * (
                x158 * x6
                + x3
                * (x117 * x47 + 2 * x120 + 2 * x121 + x148 * x53 + 2 * x149 + 2 * x150)
                + x46 * (2 * a * x188 - 3 * x116 - 6 * x187)
            )
            + x189 * x45
            + x189 * x49,
            x17 * x184 * x83 + x184 * x190 * x81 + x186 * x81,
            x105 * x191 + x17 * x183 * x99 + x191 * x49 * x99,
            x123 * x140 * x64 + x123 * x143 * x25 + x124 * x140 * x25,
            x123 * x146 * x55 + x124 * x146 * x31 + x151 * x185 * x43,
            x127 * x140 * x55 + x127 * x143 * x31 + x130 * x140 * x31,
            x124 * x154 * x17 + x154 * x190 * x43 + x159 * x43,
            x127 * x151 * x17 + x127 * x192 * x49 + x130 * x192,
            x138 * x140 * x17 + x140 * x193 * x49 + x143 * x193,
            x107 * x163 * x64 + x107 * x166 * x25 + x108 * x163 * x25,
            x111 * x163 * x55 + x111 * x166 * x31 + x114 * x163 * x31,
            x107 * x169 * x55 + x108 * x169 * x31 + x12 * x174 * x194,
            x122 * x163 * x17 + x163 * x195 * x49 + x166 * x195,
            x111 * x17 * x174 + x111 * x196 * x49 + x114 * x196,
            x108 * x17 * x177 + x177 * x194 * x49 + x180 * x6,
            x198 * x25 * x9 + x198 * x64 + x199 * x25 * x73,
            x197 * x31 * x75 + x197 * x55 * x74 + x199 * x200 * x69,
            x12 * x202 + x200 * x201 * x9 + x201 * x55 * x73,
            x17 * x199 * x86 + x203 * x49 * x86 + x203 * x93,
            x160 * x201 * x49 * x69 + x17 * x201 * x75 + x202 * x69,
            x160
            * (
                x179 * x43
                + x3
                * (x133 * x47 + 2 * x136 + 2 * x137 + x171 * x53 + 2 * x172 + 2 * x173)
                + x46 * (2 * a * x205 - 3 * x132 - 6 * x204)
            )
            + x206 * x49
            + x206 * x9,
        ]
    )
    return S


def kinetic3d_33(a, A, b, B):
    """Cartesian 3D (ff) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = -x5 - A[1]
    x7 = x6**2
    x8 = 2 * a**2
    x9 = -x0 - x8 * (x3 + x7)
    x10 = -x4 * (a * A[0] + b * B[0])
    x11 = -x10 - A[0]
    x12 = -x10 - B[0]
    x13 = a * x4
    x14 = b * x13
    x15 = numpy.exp(-x14 * (A[0] - B[0]) ** 2)
    x16 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x17 = x15 * x16
    x18 = x17 * x3
    x19 = x12 * x18
    x20 = x12**2 * x17
    x21 = x18 + x20
    x22 = x12 * x21
    x23 = x11 * x21
    x24 = 3 * x23
    x25 = 3 * x18
    x26 = x3 * (3 * x20 + x25)
    x27 = 2 * x18
    x28 = x12 * x27
    x29 = x22 + x28
    x30 = x11 * x29
    x31 = x26 + x30
    x32 = x11 * x31 + x3 * (8 * x19 + x22 + x24)
    x33 = 2 * x17
    x34 = x12 * x33
    x35 = x11 * x34 + x25
    x36 = x3 * (x20 + x35)
    x37 = x23 + x28
    x38 = x11 * x37
    x39 = x11 * x32 + x3 * (2 * x26 + 2 * x30 + 3 * x36 + 3 * x38)
    x40 = numpy.exp(-x14 * (A[1] - B[1]) ** 2)
    x41 = numpy.exp(-x14 * (A[2] - B[2]) ** 2)
    x42 = numpy.pi * x4 * x41
    x43 = x40 * x42
    x44 = x39 * x43
    x45 = -x4 * (a * A[2] + b * B[2])
    x46 = -x45 - A[2]
    x47 = x46**2
    x48 = -x0 - x8 * (x3 + x47)
    x49 = b * x4
    x50 = 4 * x14
    x51 = x11**2
    x52 = -x0 - x8 * (x3 + x51)
    x53 = x18 * x52
    x54 = 3 * x53
    x55 = -x33
    x56 = x13 * (x2 * x21 + x55)
    x57 = x12 * x17
    x58 = x1 * x49
    x59 = x52 * x57 + x57 * x58
    x60 = x12 * x59
    x61 = x3 * (x54 + 3 * x56 + 3 * x60)
    x62 = x36 + x38
    x63 = x49 * (2 * a * x62 - 2 * x20 - x27)
    x64 = x56 + x60
    x65 = x11 * x59
    x66 = x11 * x17
    x67 = x12 * x66
    x68 = x18 + x67
    x69 = x50 * x68 + x54 + 2 * x65
    x70 = x3 * (x64 + x69)
    x71 = x3 * (x34 * x52 + x50 * x57)
    x72 = x53 + x64
    x73 = x11 * x72
    x74 = x37 * x58 + x71 + x73
    x75 = x11 * x74
    x76 = x12 * x72 + x13 * (x2 * x29 - 3 * x57)
    x77 = x71 + x76
    x78 = x11 * x77
    x79 = 4 * x19
    x80 = 6 * x14
    x81 = x31 * x58 + x61 + x78
    x82 = (
        x11 * x81
        + x3 * (x37 * x80 + 4 * x71 + 3 * x73 + x76)
        + x49 * (2 * a * x32 - 2 * x22 - x79)
    )
    x83 = -x5 - B[1]
    x84 = x3 * (x57 + x66)
    x85 = x11 * x68
    x86 = x11 * x62 + x3 * (2 * x23 + x79 + 2 * x84 + 2 * x85)
    x87 = x43 * x86
    x88 = x16 * x40
    x89 = x83 * x88
    x90 = x58 * x89 + x89 * x9
    x91 = x16 * x41
    x92 = x84 + x85
    x93 = x1 * x92 - x34
    x94 = x2 * x4
    x95 = x52 * x66 + x58 * x66
    x96 = x3 * (x59 + x95)
    x97 = x53 + x58 * x68 + x65
    x98 = x11 * x97
    x99 = x63 + x70 + x75
    x100 = x43 * (
        x11 * x99
        + x3 * (x37 * x50 + 2 * x71 + 2 * x73 + x93 * x94 + 2 * x96 + 2 * x98)
        + x49 * (2 * a * x86 - 6 * x19 - x24)
    )
    x101 = -x45 - B[2]
    x102 = x101 * x91
    x103 = x102 * x48 + x102 * x58
    x104 = x3 * x88
    x105 = x83**2 * x88
    x106 = x104 + x105
    x107 = x17 * x51
    x108 = x11 * x92 + x3 * (x107 + x35)
    x109 = x108 * x91
    x110 = x104 * x9
    x111 = x83 * x90
    x112 = 2 * x88
    x113 = -x112
    x114 = x13 * (x106 * x2 + x113)
    x115 = x111 + x114
    x116 = x110 + x115
    x117 = x107 + x18
    x118 = x11 * x95 + x49 * (x1 * x117 + x55)
    x119 = x49 * x93 + x96 + x98
    x120 = x11 * x119 + x3 * (x118 + x69) + x49 * (2 * a * x108 - x25 - 3 * x67)
    x121 = x101 * x43
    x122 = x3 * x91
    x123 = x101**2 * x91
    x124 = x122 + x123
    x125 = x108 * x88
    x126 = x122 * x48
    x127 = x101 * x103
    x128 = 2 * x91
    x129 = -x128
    x130 = x13 * (x124 * x2 + x129)
    x131 = x127 + x130
    x132 = x126 + x131
    x133 = 2 * x104
    x134 = x133 * x83
    x135 = x106 * x83
    x136 = x134 + x135
    x137 = x11 * x117 + x11 * x27
    x138 = x137 * x91
    x139 = x112 * x83
    x140 = x3 * (x139 * x9 + x50 * x89)
    x141 = x116 * x83 + x13 * (x136 * x2 - 3 * x89)
    x142 = x140 + x141
    x143 = x118 + x53
    x144 = x11 * x143 + x3 * (x11 * x33 * x52 + x50 * x66) + x49 * (x1 * x137 - 3 * x66)
    x145 = 2 * x122
    x146 = x101 * x145
    x147 = x101 * x124
    x148 = x146 + x147
    x149 = x137 * x88
    x150 = x101 * x128
    x151 = x3 * (x102 * x50 + x150 * x48)
    x152 = x101 * x132 + x13 * (-3 * x102 + x148 * x2)
    x153 = x151 + x152
    x154 = x32 * x43
    x155 = x6 * x88
    x156 = x155 * x58 + x155 * x9
    x157 = x43 * x82
    x158 = x155 * x83
    x159 = x104 + x158
    x160 = x62 * x91
    x161 = x6 * x90
    x162 = x110 + x159 * x58 + x161
    x163 = x106 * x6
    x164 = x134 + x163
    x165 = x91 * x92
    x166 = x116 * x6
    x167 = x140 + x164 * x58 + x166
    x168 = 3 * x104
    x169 = x3 * (3 * x105 + x168)
    x170 = x136 * x6
    x171 = x169 + x170
    x172 = x171 * x91
    x173 = 3 * x110
    x174 = x3 * (3 * x111 + 3 * x114 + x173)
    x175 = x142 * x6
    x176 = x171 * x58 + x174 + x175
    x177 = x46 * x91
    x178 = x177 * x48 + x177 * x58
    x179 = x43 * x46
    x180 = x101 * x177
    x181 = x122 + x180
    x182 = x62 * x88
    x183 = x103 * x46
    x184 = x126 + x181 * x58 + x183
    x185 = x124 * x46
    x186 = x146 + x185
    x187 = x88 * x92
    x188 = x132 * x46
    x189 = x151 + x186 * x58 + x188
    x190 = 3 * x122
    x191 = x3 * (3 * x123 + x190)
    x192 = x148 * x46
    x193 = x191 + x192
    x194 = x193 * x88
    x195 = 3 * x126
    x196 = x3 * (3 * x127 + 3 * x130 + x195)
    x197 = x153 * x46
    x198 = x193 * x58 + x196 + x197
    x199 = x7 * x88
    x200 = x104 + x199
    x201 = x31 * x91
    x202 = x156 * x6 + x49 * (x1 * x200 + x113)
    x203 = x110 + x202
    x204 = x3 * (x155 + x89)
    x205 = x159 * x6
    x206 = x204 + x205
    x207 = x206 * x91
    x208 = x1 * x206 - x139
    x209 = x3 * (x156 + x90)
    x210 = x162 * x6
    x211 = x208 * x49 + x209 + x210
    x212 = x139 * x6 + x168
    x213 = x3 * (x105 + x212)
    x214 = x164 * x6
    x215 = x213 + x214
    x216 = x215 * x91
    x217 = x49 * (2 * a * x215 - 2 * x105 - x133)
    x218 = x159 * x50 + 2 * x161 + x173
    x219 = x3 * (x115 + x218)
    x220 = x167 * x6
    x221 = x217 + x219 + x220
    x222 = x104 * x83
    x223 = 3 * x163
    x224 = x171 * x6 + x3 * (x135 + 8 * x222 + x223)
    x225 = x15 * x42
    x226 = x11 * x225
    x227 = 4 * x222
    x228 = (
        x176 * x6
        + x3 * (4 * x140 + x141 + x164 * x80 + 3 * x166)
        + x49 * (2 * a * x224 - 2 * x135 - x227)
    )
    x229 = x225 * x228
    x230 = numpy.pi * x15 * x4 * x40
    x231 = x11 * x230
    x232 = x47 * x91
    x233 = x122 + x232
    x234 = x31 * x88
    x235 = x178 * x46 + x49 * (x1 * x233 + x129)
    x236 = x126 + x235
    x237 = x3 * (x102 + x177)
    x238 = x181 * x46
    x239 = x237 + x238
    x240 = x239 * x88
    x241 = x1 * x239 - x150
    x242 = x3 * (x103 + x178)
    x243 = x184 * x46
    x244 = x241 * x49 + x242 + x243
    x245 = x150 * x46 + x190
    x246 = x3 * (x123 + x245)
    x247 = x186 * x46
    x248 = x246 + x247
    x249 = x248 * x88
    x250 = x49 * (2 * a * x248 - 2 * x123 - x145)
    x251 = x181 * x50 + 2 * x183 + x195
    x252 = x3 * (x131 + x251)
    x253 = x189 * x46
    x254 = x250 + x252 + x253
    x255 = x101 * x122
    x256 = 3 * x185
    x257 = x193 * x46 + x3 * (x147 + 8 * x255 + x256)
    x258 = 4 * x255
    x259 = (
        x198 * x46
        + x3 * (4 * x151 + x152 + x186 * x80 + 3 * x188)
        + x49 * (2 * a * x257 - 2 * x147 - x258)
    )
    x260 = x230 * x259
    x261 = x133 * x6 + x200 * x6
    x262 = x261 * x91
    x263 = x203 * x6 + x3 * (x112 * x6 * x9 + x155 * x50) + x49 * (x1 * x261 - 3 * x155)
    x264 = x206 * x6 + x3 * (x199 + x212)
    x265 = x264 * x91
    x266 = x211 * x6 + x3 * (x202 + x218) + x49 * (2 * a * x264 - 3 * x158 - x168)
    x267 = x215 * x6 + x3 * (2 * x163 + 2 * x204 + 2 * x205 + x227)
    x268 = x12 * x225
    x269 = x225 * (
        x221 * x6
        + x3 * (2 * x140 + x164 * x50 + 2 * x166 + x208 * x94 + 2 * x209 + 2 * x210)
        + x49 * (2 * a * x267 - 6 * x222 - x223)
    )
    x270 = x224 * x6 + x3 * (2 * x169 + 2 * x170 + 3 * x213 + 3 * x214)
    x271 = x225 * x270
    x272 = x225 * x52
    x273 = x17 * x264
    x274 = x17 * x261
    x275 = x17 * x215
    x276 = x17 * x206
    x277 = x17 * x193
    x278 = x230 * x6
    x279 = x17 * x171
    x280 = x17 * x239
    x281 = x17 * x248
    x282 = x145 * x46 + x233 * x46
    x283 = x282 * x88
    x284 = (
        x236 * x46 + x3 * (x128 * x46 * x48 + x177 * x50) + x49 * (x1 * x282 - 3 * x177)
    )
    x285 = x239 * x46 + x3 * (x232 + x245)
    x286 = x285 * x88
    x287 = x244 * x46 + x3 * (x235 + x251) + x49 * (2 * a * x285 - 3 * x180 - x190)
    x288 = x12 * x230
    x289 = x248 * x46 + x3 * (2 * x185 + 2 * x237 + 2 * x238 + x258)
    x290 = x230 * (
        x254 * x46
        + x3 * (2 * x151 + x186 * x50 + 2 * x188 + x241 * x94 + 2 * x242 + 2 * x243)
        + x49 * (2 * a * x289 - 6 * x255 - x256)
    )
    x291 = x17 * x282
    x292 = x17 * x285
    x293 = x257 * x46 + x3 * (2 * x191 + 2 * x192 + 3 * x246 + 3 * x247)
    x294 = x230 * x293

    # 100 item(s)
    S = numpy.array(
        [
            x43
            * (
                x11 * x82
                + x3 * (x31 * x50 + 2 * x61 + 3 * x63 + 3 * x70 + 3 * x75 + 2 * x78)
                + x49 * (2 * a * x39 - 3 * x26 - 3 * x30)
            )
            + x44 * x48
            + x44 * x9,
            x100 * x83 + x48 * x83 * x87 + x86 * x90 * x91,
            x100 * x101 + x101 * x87 * x9 + x103 * x86 * x88,
            x106 * x109 * x48 + x106 * x120 * x91 + x109 * x116,
            x102 * x108 * x90 + x103 * x108 * x89 + x120 * x121 * x83,
            x120 * x124 * x88 + x124 * x125 * x9 + x125 * x132,
            x136 * x138 * x48 + x136 * x144 * x91 + x138 * x142,
            x102 * x106 * x144 + x102 * x116 * x137 + x103 * x106 * x137,
            x124 * x137 * x90 + x124 * x144 * x89 + x132 * x137 * x89,
            x144 * x148 * x88 + x148 * x149 * x9 + x149 * x153,
            x154 * x48 * x6 + x156 * x32 * x91 + x157 * x6,
            x159 * x160 * x48 + x159 * x91 * x99 + x160 * x162,
            x102 * x156 * x62 + x103 * x155 * x62 + x121 * x6 * x99,
            x119 * x164 * x91 + x164 * x165 * x48 + x165 * x167,
            x102 * x119 * x159 + x102 * x162 * x92 + x103 * x159 * x92,
            x119 * x124 * x155 + x124 * x156 * x92 + x132 * x155 * x92,
            x117 * x172 * x48 + x117 * x176 * x91 + x143 * x172,
            x102 * x117 * x167 + x102 * x143 * x164 + x103 * x117 * x164,
            x117 * x124 * x162 + x117 * x132 * x159 + x124 * x143 * x159,
            x117 * x148 * x156 + x117 * x153 * x155 + x143 * x148 * x155,
            x154 * x46 * x9 + x157 * x46 + x178 * x32 * x88,
            x177 * x62 * x90 + x178 * x62 * x89 + x179 * x83 * x99,
            x181 * x182 * x9 + x181 * x88 * x99 + x182 * x184,
            x106 * x119 * x177 + x106 * x178 * x92 + x116 * x177 * x92,
            x119 * x181 * x89 + x181 * x90 * x92 + x184 * x89 * x92,
            x119 * x186 * x88 + x186 * x187 * x9 + x187 * x189,
            x117 * x136 * x178 + x117 * x142 * x177 + x136 * x143 * x177,
            x106 * x117 * x184 + x106 * x143 * x181 + x116 * x117 * x181,
            x117 * x186 * x90 + x117 * x189 * x89 + x143 * x186 * x89,
            x117 * x194 * x9 + x117 * x198 * x88 + x143 * x194,
            x200 * x201 * x48 + x200 * x81 * x91 + x201 * x203,
            x207 * x37 * x48 + x207 * x74 + x211 * x37 * x91,
            x102 * x200 * x74 + x102 * x203 * x37 + x103 * x200 * x37,
            x216 * x48 * x68 + x216 * x97 + x221 * x68 * x91,
            x102 * x206 * x97 + x102 * x211 * x68 + x103 * x206 * x68,
            x124 * x200 * x97 + x124 * x203 * x68 + x132 * x200 * x68,
            x11 * x229 + x224 * x226 * x48 + x224 * x91 * x95,
            x101 * x221 * x226 + x102 * x215 * x95 + x103 * x215 * x66,
            x124 * x206 * x95 + x124 * x211 * x66 + x132 * x206 * x66,
            x148 * x200 * x95 + x148 * x203 * x66 + x153 * x200 * x66,
            x155 * x178 * x31 + x156 * x177 * x31 + x179 * x6 * x81,
            x159 * x177 * x74 + x159 * x178 * x37 + x162 * x177 * x37,
            x155 * x181 * x74 + x155 * x184 * x37 + x156 * x181 * x37,
            x164 * x177 * x97 + x164 * x178 * x68 + x167 * x177 * x68,
            x159 * x181 * x97 + x159 * x184 * x68 + x162 * x181 * x68,
            x155 * x186 * x97 + x155 * x189 * x68 + x156 * x186 * x68,
            x171 * x177 * x95 + x171 * x178 * x66 + x176 * x226 * x46,
            x164 * x181 * x95 + x164 * x184 * x66 + x167 * x181 * x66,
            x159 * x186 * x95 + x159 * x189 * x66 + x162 * x186 * x66,
            x155 * x193 * x95 + x156 * x193 * x66 + x198 * x231 * x6,
            x233 * x234 * x9 + x233 * x81 * x88 + x234 * x236,
            x233 * x37 * x90 + x233 * x74 * x89 + x236 * x37 * x89,
            x240 * x37 * x9 + x240 * x74 + x244 * x37 * x88,
            x106 * x233 * x97 + x106 * x236 * x68 + x116 * x233 * x68,
            x239 * x68 * x90 + x239 * x89 * x97 + x244 * x68 * x89,
            x249 * x68 * x9 + x249 * x97 + x254 * x68 * x88,
            x136 * x233 * x95 + x136 * x236 * x66 + x142 * x233 * x66,
            x106 * x239 * x95 + x106 * x244 * x66 + x116 * x239 * x66,
            x231 * x254 * x83 + x248 * x66 * x90 + x248 * x89 * x95,
            x11 * x260 + x231 * x257 * x9 + x257 * x88 * x95,
            x262 * x29 * x48 + x262 * x77 + x263 * x29 * x91,
            x21 * x265 * x48 + x21 * x266 * x91 + x265 * x72,
            x102 * x21 * x263 + x102 * x261 * x72 + x103 * x21 * x261,
            x12 * x269 + x267 * x268 * x48 + x267 * x59 * x91,
            x101 * x266 * x268 + x102 * x264 * x59 + x103 * x264 * x57,
            x124 * x261 * x59 + x124 * x263 * x57 + x132 * x261 * x57,
            x225
            * (
                x228 * x6
                + x3 * (x171 * x50 + 2 * x174 + 2 * x175 + 3 * x217 + 3 * x219 + 3 * x220)
                + x49 * (2 * a * x270 - 3 * x169 - 3 * x170)
            )
            + x271 * x48
            + x271 * x52,
            x101 * x267 * x272 + x101 * x269 + x103 * x17 * x267,
            x124 * x17 * x266 + x124 * x273 * x52 + x132 * x273,
            x148 * x17 * x263 + x148 * x274 * x52 + x153 * x274,
            x177 * x200 * x77 + x177 * x203 * x29 + x178 * x200 * x29,
            x177 * x206 * x72 + x177 * x21 * x211 + x178 * x206 * x21,
            x181 * x200 * x72 + x181 * x203 * x21 + x184 * x200 * x21,
            x177 * x215 * x59 + x178 * x215 * x57 + x221 * x268 * x46,
            x181 * x206 * x59 + x181 * x211 * x57 + x184 * x206 * x57,
            x186 * x200 * x59 + x186 * x203 * x57 + x189 * x200 * x57,
            x17 * x178 * x224 + x224 * x272 * x46 + x229 * x46,
            x17 * x181 * x221 + x181 * x275 * x52 + x184 * x275,
            x17 * x186 * x211 + x186 * x276 * x52 + x189 * x276,
            x17 * x198 * x200 + x200 * x277 * x52 + x203 * x277,
            x155 * x233 * x77 + x155 * x236 * x29 + x156 * x233 * x29,
            x159 * x21 * x236 + x159 * x233 * x72 + x162 * x21 * x233,
            x155 * x21 * x244 + x155 * x239 * x72 + x156 * x21 * x239,
            x164 * x233 * x59 + x164 * x236 * x57 + x167 * x233 * x57,
            x159 * x239 * x59 + x159 * x244 * x57 + x162 * x239 * x57,
            x12 * x254 * x278 + x155 * x248 * x59 + x156 * x248 * x57,
            x17 * x176 * x233 + x233 * x279 * x52 + x236 * x279,
            x164 * x17 * x244 + x164 * x280 * x52 + x167 * x280,
            x159 * x17 * x254 + x159 * x281 * x52 + x162 * x281,
            x156 * x17 * x257 + x257 * x278 * x52 + x260 * x6,
            x283 * x29 * x9 + x283 * x77 + x284 * x29 * x88,
            x21 * x282 * x90 + x21 * x284 * x89 + x282 * x72 * x89,
            x21 * x286 * x9 + x21 * x287 * x88 + x286 * x72,
            x106 * x282 * x59 + x106 * x284 * x57 + x116 * x282 * x57,
            x285 * x57 * x90 + x285 * x59 * x89 + x287 * x288 * x83,
            x12 * x290 + x288 * x289 * x9 + x289 * x59 * x88,
            x136 * x17 * x284 + x136 * x291 * x52 + x142 * x291,
            x106 * x17 * x287 + x106 * x292 * x52 + x116 * x292,
            x17 * x289 * x90 + x230 * x289 * x52 * x83 + x290 * x83,
            x230
            * (
                x259 * x46
                + x3 * (x193 * x50 + 2 * x196 + 2 * x197 + 3 * x250 + 3 * x252 + 3 * x253)
                + x49 * (2 * a * x293 - 3 * x191 - 3 * x192)
            )
            + x294 * x52
            + x294 * x9,
        ]
    )
    return S


def kinetic3d_34(a, A, b, B):
    """Cartesian 3D (fg) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = -x5 - A[1]
    x7 = x6**2
    x8 = 2 * a**2
    x9 = -x0 - x8 * (x3 + x7)
    x10 = -x4 * (a * A[0] + b * B[0])
    x11 = -x10 - A[0]
    x12 = a * x4
    x13 = b * x12
    x14 = numpy.exp(-x13 * (A[0] - B[0]) ** 2)
    x15 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x16 = x14 * x15
    x17 = x16 * x3
    x18 = 3 * x17
    x19 = -x10 - B[0]
    x20 = x16 * x19**2
    x21 = x3 * (x18 + 3 * x20)
    x22 = 2 * x17
    x23 = x19 * x22
    x24 = x17 + x20
    x25 = x19 * x24
    x26 = x23 + x25
    x27 = x19 * x26
    x28 = x11 * x26
    x29 = x17 * x19
    x30 = 8 * x29
    x31 = x3 * (4 * x25 + x30)
    x32 = x21 + x27
    x33 = x11 * x32
    x34 = x31 + x33
    x35 = x11 * x34 + x3 * (5 * x21 + x27 + 4 * x28)
    x36 = x11 * x24
    x37 = 3 * x36
    x38 = x3 * (x25 + x30 + x37)
    x39 = x21 + x28
    x40 = x11 * x39
    x41 = x11 * x35 + x3 * (2 * x31 + 2 * x33 + 4 * x38 + 4 * x40)
    x42 = numpy.exp(-x13 * (A[1] - B[1]) ** 2)
    x43 = numpy.exp(-x13 * (A[2] - B[2]) ** 2)
    x44 = numpy.pi * x4 * x43
    x45 = x42 * x44
    x46 = x41 * x45
    x47 = -x4 * (a * A[2] + b * B[2])
    x48 = -x47 - A[2]
    x49 = x48**2
    x50 = -x0 - x8 * (x3 + x49)
    x51 = b * x4
    x52 = 4 * x13
    x53 = 4 * x17
    x54 = x19 * x53
    x55 = x38 + x40
    x56 = x51 * (2 * a * x55 - 2 * x25 - x54)
    x57 = x16 * x19
    x58 = x11**2
    x59 = -x0 - x8 * (x3 + x58)
    x60 = 2 * x16
    x61 = x19 * x60
    x62 = x3 * (x52 * x57 + x59 * x61)
    x63 = 4 * x62
    x64 = x12 * (x2 * x26 - 3 * x57)
    x65 = x17 * x59
    x66 = x1 * x51
    x67 = x57 * x59 + x57 * x66
    x68 = x19 * x67
    x69 = -x60
    x70 = x12 * (x2 * x24 + x69)
    x71 = x68 + x70
    x72 = x65 + x71
    x73 = x19 * x72
    x74 = x3 * (x63 + 4 * x64 + 4 * x73)
    x75 = x11 * x72
    x76 = x23 + x36
    x77 = x13 * x76
    x78 = x64 + x73
    x79 = x3 * (x63 + 3 * x75 + 6 * x77 + x78)
    x80 = 3 * x65
    x81 = x3 * (3 * x68 + 3 * x70 + x80)
    x82 = x62 + x78
    x83 = x11 * x82
    x84 = x39 * x66 + x81 + x83
    x85 = x11 * x84
    x86 = x12 * (2 * b * x32 - 4 * x20 - x53) + x19 * x82
    x87 = x81 + x86
    x88 = x11 * x87
    x89 = 2 * x21
    x90 = 8 * x13
    x91 = x34 * x66 + x74 + x88
    x92 = (
        x11 * x91
        + x3 * (x39 * x90 + 5 * x81 + 4 * x83 + x86)
        + x51 * (2 * a * x35 - 2 * x27 - x89)
    )
    x93 = -x5 - B[1]
    x94 = x11 * x61 + x18
    x95 = x3 * (x20 + x94)
    x96 = x11 * x76
    x97 = x11 * x55 + x3 * (2 * x28 + x89 + 3 * x95 + 3 * x96)
    x98 = x45 * x97
    x99 = x15 * x42
    x100 = x93 * x99
    x101 = x100 * x66 + x100 * x9
    x102 = x15 * x43
    x103 = x95 + x96
    x104 = x51 * (2 * a * x103 - 2 * x20 - x22)
    x105 = x11 * x67
    x106 = x11 * x16
    x107 = x106 * x19
    x108 = x107 + x17
    x109 = 2 * x105 + x108 * x52 + x80
    x110 = x3 * (x109 + x71)
    x111 = x62 + x66 * x76 + x75
    x112 = x11 * x111
    x113 = x56 + x79 + x85
    x114 = x45 * (
        x11 * x113
        + x3 * (3 * x104 + 3 * x110 + 3 * x112 + x39 * x52 + 2 * x81 + 2 * x83)
        + x51 * (2 * a * x97 - 3 * x21 - 3 * x28)
    )
    x115 = -x47 - B[2]
    x116 = x102 * x115
    x117 = x116 * x50 + x116 * x66
    x118 = x3 * x99
    x119 = x93**2 * x99
    x120 = x118 + x119
    x121 = x3 * (x106 + x57)
    x122 = x108 * x11
    x123 = x103 * x11 + x3 * (2 * x121 + 2 * x122 + 2 * x36 + x54)
    x124 = x102 * x123
    x125 = x118 * x9
    x126 = x101 * x93
    x127 = 2 * x99
    x128 = -x127
    x129 = x12 * (x120 * x2 + x128)
    x130 = x126 + x129
    x131 = x125 + x130
    x132 = x121 + x122
    x133 = x1 * x132 - x61
    x134 = x2 * x4
    x135 = x106 * x59 + x106 * x66
    x136 = x3 * (x135 + x67)
    x137 = x105 + x108 * x66 + x65
    x138 = x11 * x137
    x139 = x104 + x110 + x112
    x140 = (
        x11 * x139
        + x3 * (x133 * x134 + 2 * x136 + 2 * x138 + 2 * x62 + 2 * x75 + 4 * x77)
        + x51 * (2 * a * x123 - 6 * x29 - x37)
    )
    x141 = x115 * x45
    x142 = x102 * x3
    x143 = x102 * x115**2
    x144 = x142 + x143
    x145 = x123 * x99
    x146 = x142 * x50
    x147 = x115 * x117
    x148 = 2 * x102
    x149 = -x148
    x150 = x12 * (x144 * x2 + x149)
    x151 = x147 + x150
    x152 = x146 + x151
    x153 = 2 * x118
    x154 = x153 * x93
    x155 = x120 * x93
    x156 = x154 + x155
    x157 = x16 * x58
    x158 = x11 * x132 + x3 * (x157 + x94)
    x159 = x102 * x158
    x160 = x127 * x93
    x161 = x3 * (x100 * x52 + x160 * x9)
    x162 = x131 * x93
    x163 = x12 * (-3 * x100 + x156 * x2)
    x164 = x162 + x163
    x165 = x161 + x164
    x166 = x157 + x17
    x167 = x11 * x135 + x51 * (x1 * x166 + x69)
    x168 = x133 * x51 + x136 + x138
    x169 = x11 * x168 + x3 * (x109 + x167) + x51 * (2 * a * x158 - 3 * x107 - x18)
    x170 = 2 * x142
    x171 = x115 * x170
    x172 = x115 * x144
    x173 = x171 + x172
    x174 = x158 * x99
    x175 = x115 * x148
    x176 = x3 * (x116 * x52 + x175 * x50)
    x177 = x115 * x152
    x178 = x12 * (-3 * x116 + x173 * x2)
    x179 = x177 + x178
    x180 = x176 + x179
    x181 = x11 * x166 + x11 * x22
    x182 = 3 * x118
    x183 = x3 * (3 * x119 + x182)
    x184 = x156 * x93
    x185 = x183 + x184
    x186 = x102 * x185
    x187 = x167 + x65
    x188 = x11 * x187 + x3 * (x106 * x52 + x11 * x59 * x60) + x51 * (x1 * x181 - 3 * x106)
    x189 = 3 * x125
    x190 = x3 * (3 * x126 + 3 * x129 + x189)
    x191 = 4 * x118
    x192 = x12 * (2 * b * x185 - 4 * x119 - x191) + x165 * x93
    x193 = x190 + x192
    x194 = 3 * x142
    x195 = x3 * (3 * x143 + x194)
    x196 = x115 * x173
    x197 = x195 + x196
    x198 = x197 * x99
    x199 = 3 * x146
    x200 = x3 * (3 * x147 + 3 * x150 + x199)
    x201 = 4 * x142
    x202 = x115 * x180 + x12 * (2 * b * x197 - 4 * x143 - x201)
    x203 = x200 + x202
    x204 = x35 * x45
    x205 = x6 * x99
    x206 = x205 * x66 + x205 * x9
    x207 = x45 * x92
    x208 = x205 * x93
    x209 = x118 + x208
    x210 = x102 * x55
    x211 = x101 * x6
    x212 = x125 + x209 * x66 + x211
    x213 = x120 * x6
    x214 = x154 + x213
    x215 = x102 * x103
    x216 = x131 * x6
    x217 = x161 + x214 * x66 + x216
    x218 = x156 * x6
    x219 = x183 + x218
    x220 = x102 * x132
    x221 = x165 * x6
    x222 = x190 + x219 * x66 + x221
    x223 = x118 * x93
    x224 = 8 * x223
    x225 = x3 * (4 * x155 + x224)
    x226 = x185 * x6
    x227 = x225 + x226
    x228 = x102 * x227
    x229 = 4 * x161
    x230 = x3 * (4 * x162 + 4 * x163 + x229)
    x231 = x193 * x6
    x232 = x227 * x66 + x230 + x231
    x233 = x102 * x48
    x234 = x233 * x50 + x233 * x66
    x235 = x45 * x48
    x236 = x115 * x233
    x237 = x142 + x236
    x238 = x55 * x99
    x239 = x117 * x48
    x240 = x146 + x237 * x66 + x239
    x241 = x144 * x48
    x242 = x171 + x241
    x243 = x103 * x99
    x244 = x152 * x48
    x245 = x176 + x242 * x66 + x244
    x246 = x173 * x48
    x247 = x195 + x246
    x248 = x132 * x99
    x249 = x180 * x48
    x250 = x200 + x247 * x66 + x249
    x251 = x115 * x142
    x252 = 8 * x251
    x253 = x3 * (4 * x172 + x252)
    x254 = x197 * x48
    x255 = x253 + x254
    x256 = x255 * x99
    x257 = 4 * x176
    x258 = x3 * (4 * x177 + 4 * x178 + x257)
    x259 = x203 * x48
    x260 = x255 * x66 + x258 + x259
    x261 = x7 * x99
    x262 = x118 + x261
    x263 = x102 * x34
    x264 = x206 * x6 + x51 * (x1 * x262 + x128)
    x265 = x125 + x264
    x266 = x3 * (x100 + x205)
    x267 = x209 * x6
    x268 = x266 + x267
    x269 = x102 * x268
    x270 = x1 * x268 - x160
    x271 = x3 * (x101 + x206)
    x272 = x212 * x6
    x273 = x270 * x51 + x271 + x272
    x274 = x160 * x6 + x182
    x275 = x3 * (x119 + x274)
    x276 = x214 * x6
    x277 = x275 + x276
    x278 = x102 * x277
    x279 = x51 * (2 * a * x277 - 2 * x119 - x153)
    x280 = x189 + x209 * x52 + 2 * x211
    x281 = x3 * (x130 + x280)
    x282 = x217 * x6
    x283 = x279 + x281 + x282
    x284 = 3 * x213
    x285 = x3 * (x155 + x224 + x284)
    x286 = x219 * x6
    x287 = x285 + x286
    x288 = x102 * x287
    x289 = x191 * x93
    x290 = x51 * (2 * a * x287 - 2 * x155 - x289)
    x291 = 6 * x13
    x292 = x3 * (x164 + x214 * x291 + 3 * x216 + x229)
    x293 = x222 * x6
    x294 = x290 + x292 + x293
    x295 = x227 * x6 + x3 * (5 * x183 + x184 + 4 * x218)
    x296 = x14 * x44
    x297 = x11 * x296
    x298 = 2 * x183
    x299 = (
        x232 * x6
        + x3 * (5 * x190 + x192 + x219 * x90 + 4 * x221)
        + x51 * (2 * a * x295 - 2 * x184 - x298)
    )
    x300 = x296 * x299
    x301 = numpy.pi * x14 * x4 * x42
    x302 = x11 * x301
    x303 = x102 * x49
    x304 = x142 + x303
    x305 = x34 * x99
    x306 = x234 * x48 + x51 * (x1 * x304 + x149)
    x307 = x146 + x306
    x308 = x3 * (x116 + x233)
    x309 = x237 * x48
    x310 = x308 + x309
    x311 = x310 * x99
    x312 = x1 * x310 - x175
    x313 = x3 * (x117 + x234)
    x314 = x240 * x48
    x315 = x312 * x51 + x313 + x314
    x316 = x175 * x48 + x194
    x317 = x3 * (x143 + x316)
    x318 = x242 * x48
    x319 = x317 + x318
    x320 = x319 * x99
    x321 = x51 * (2 * a * x319 - 2 * x143 - x170)
    x322 = x199 + x237 * x52 + 2 * x239
    x323 = x3 * (x151 + x322)
    x324 = x245 * x48
    x325 = x321 + x323 + x324
    x326 = 3 * x241
    x327 = x3 * (x172 + x252 + x326)
    x328 = x247 * x48
    x329 = x327 + x328
    x330 = x329 * x99
    x331 = x115 * x201
    x332 = x51 * (2 * a * x329 - 2 * x172 - x331)
    x333 = x3 * (x179 + x242 * x291 + 3 * x244 + x257)
    x334 = x250 * x48
    x335 = x332 + x333 + x334
    x336 = x255 * x48 + x3 * (5 * x195 + x196 + 4 * x246)
    x337 = 2 * x195
    x338 = (
        x260 * x48
        + x3 * (5 * x200 + x202 + x247 * x90 + 4 * x249)
        + x51 * (2 * a * x336 - 2 * x196 - x337)
    )
    x339 = x301 * x338
    x340 = x153 * x6 + x262 * x6
    x341 = x102 * x32
    x342 = x265 * x6 + x3 * (x127 * x6 * x9 + x205 * x52) + x51 * (x1 * x340 - 3 * x205)
    x343 = x268 * x6 + x3 * (x261 + x274)
    x344 = x102 * x343
    x345 = x273 * x6 + x3 * (x264 + x280) + x51 * (2 * a * x343 - x182 - 3 * x208)
    x346 = x277 * x6 + x3 * (2 * x213 + 2 * x266 + 2 * x267 + x289)
    x347 = x102 * x346
    x348 = (
        x283 * x6
        + x3 * (x134 * x270 + 2 * x161 + x214 * x52 + 2 * x216 + 2 * x271 + 2 * x272)
        + x51 * (2 * a * x346 - 6 * x223 - x284)
    )
    x349 = x287 * x6 + x3 * (2 * x218 + 3 * x275 + 3 * x276 + x298)
    x350 = x19 * x296
    x351 = x296 * (
        x294 * x6
        + x3 * (2 * x190 + x219 * x52 + 2 * x221 + 3 * x279 + 3 * x281 + 3 * x282)
        + x51 * (2 * a * x349 - 3 * x183 - 3 * x218)
    )
    x352 = x295 * x6 + x3 * (2 * x225 + 2 * x226 + 4 * x285 + 4 * x286)
    x353 = x296 * x352
    x354 = x296 * x59
    x355 = x16 * x346
    x356 = x16 * x343
    x357 = x16 * x197
    x358 = x16 * x287
    x359 = x16 * x277
    x360 = x16 * x268
    x361 = x16 * x255
    x362 = x301 * x6
    x363 = x16 * x227
    x364 = x16 * x310
    x365 = x16 * x319
    x366 = x16 * x329
    x367 = x170 * x48 + x304 * x48
    x368 = x32 * x99
    x369 = (
        x3 * (x148 * x48 * x50 + x233 * x52) + x307 * x48 + x51 * (x1 * x367 - 3 * x233)
    )
    x370 = x3 * (x303 + x316) + x310 * x48
    x371 = x370 * x99
    x372 = x3 * (x306 + x322) + x315 * x48 + x51 * (2 * a * x370 - x194 - 3 * x236)
    x373 = x3 * (2 * x241 + 2 * x308 + 2 * x309 + x331) + x319 * x48
    x374 = x373 * x99
    x375 = (
        x3 * (x134 * x312 + 2 * x176 + x242 * x52 + 2 * x244 + 2 * x313 + 2 * x314)
        + x325 * x48
        + x51 * (2 * a * x373 - 6 * x251 - x326)
    )
    x376 = x19 * x301
    x377 = x3 * (2 * x246 + 3 * x317 + 3 * x318 + x337) + x329 * x48
    x378 = x301 * (
        x3 * (2 * x200 + x247 * x52 + 2 * x249 + 3 * x321 + 3 * x323 + 3 * x324)
        + x335 * x48
        + x51 * (2 * a * x377 - 3 * x195 - 3 * x246)
    )
    x379 = x16 * x185
    x380 = x16 * x370
    x381 = x16 * x373
    x382 = x3 * (2 * x253 + 2 * x254 + 4 * x327 + 4 * x328) + x336 * x48
    x383 = x301 * x382

    # 150 item(s)
    S = numpy.array(
        [
            x45
            * (
                x11 * x92
                + x3 * (x34 * x52 + 4 * x56 + 2 * x74 + 4 * x79 + 4 * x85 + 2 * x88)
                + x51 * (2 * a * x41 - 3 * x31 - 3 * x33)
            )
            + x46 * x50
            + x46 * x9,
            x101 * x102 * x97 + x114 * x93 + x50 * x93 * x98,
            x114 * x115 + x115 * x9 * x98 + x117 * x97 * x99,
            x102 * x120 * x140 + x120 * x124 * x50 + x124 * x131,
            x100 * x117 * x123 + x101 * x116 * x123 + x140 * x141 * x93,
            x140 * x144 * x99 + x144 * x145 * x9 + x145 * x152,
            x102 * x156 * x169 + x156 * x159 * x50 + x159 * x165,
            x116 * x120 * x169 + x116 * x131 * x158 + x117 * x120 * x158,
            x100 * x144 * x169 + x100 * x152 * x158 + x101 * x144 * x158,
            x169 * x173 * x99 + x173 * x174 * x9 + x174 * x180,
            x102 * x181 * x193 + x181 * x186 * x50 + x186 * x188,
            x116 * x156 * x188 + x116 * x165 * x181 + x117 * x156 * x181,
            x120 * x144 * x188 + x120 * x152 * x181 + x131 * x144 * x181,
            x100 * x173 * x188 + x100 * x180 * x181 + x101 * x173 * x181,
            x181 * x198 * x9 + x181 * x203 * x99 + x188 * x198,
            x102 * x206 * x35 + x204 * x50 * x6 + x207 * x6,
            x102 * x113 * x209 + x209 * x210 * x50 + x210 * x212,
            x113 * x141 * x6 + x116 * x206 * x55 + x117 * x205 * x55,
            x102 * x139 * x214 + x214 * x215 * x50 + x215 * x217,
            x103 * x116 * x212 + x103 * x117 * x209 + x116 * x139 * x209,
            x103 * x144 * x206 + x103 * x152 * x205 + x139 * x144 * x205,
            x102 * x168 * x219 + x219 * x220 * x50 + x220 * x222,
            x116 * x132 * x217 + x116 * x168 * x214 + x117 * x132 * x214,
            x132 * x144 * x212 + x132 * x152 * x209 + x144 * x168 * x209,
            x132 * x173 * x206 + x132 * x180 * x205 + x168 * x173 * x205,
            x102 * x166 * x232 + x166 * x228 * x50 + x187 * x228,
            x116 * x166 * x222 + x116 * x187 * x219 + x117 * x166 * x219,
            x144 * x166 * x217 + x144 * x187 * x214 + x152 * x166 * x214,
            x166 * x173 * x212 + x166 * x180 * x209 + x173 * x187 * x209,
            x166 * x197 * x206 + x166 * x203 * x205 + x187 * x197 * x205,
            x204 * x48 * x9 + x207 * x48 + x234 * x35 * x99,
            x100 * x234 * x55 + x101 * x233 * x55 + x113 * x235 * x93,
            x113 * x237 * x99 + x237 * x238 * x9 + x238 * x240,
            x103 * x120 * x234 + x103 * x131 * x233 + x120 * x139 * x233,
            x100 * x103 * x240 + x100 * x139 * x237 + x101 * x103 * x237,
            x139 * x242 * x99 + x242 * x243 * x9 + x243 * x245,
            x132 * x156 * x234 + x132 * x165 * x233 + x156 * x168 * x233,
            x120 * x132 * x240 + x120 * x168 * x237 + x131 * x132 * x237,
            x100 * x132 * x245 + x100 * x168 * x242 + x101 * x132 * x242,
            x168 * x247 * x99 + x247 * x248 * x9 + x248 * x250,
            x166 * x185 * x234 + x166 * x193 * x233 + x185 * x187 * x233,
            x156 * x166 * x240 + x156 * x187 * x237 + x165 * x166 * x237,
            x120 * x166 * x245 + x120 * x187 * x242 + x131 * x166 * x242,
            x100 * x166 * x250 + x100 * x187 * x247 + x101 * x166 * x247,
            x166 * x256 * x9 + x166 * x260 * x99 + x187 * x256,
            x102 * x262 * x91 + x262 * x263 * x50 + x263 * x265,
            x102 * x273 * x39 + x269 * x39 * x50 + x269 * x84,
            x116 * x262 * x84 + x116 * x265 * x39 + x117 * x262 * x39,
            x102 * x283 * x76 + x111 * x278 + x278 * x50 * x76,
            x111 * x116 * x268 + x116 * x273 * x76 + x117 * x268 * x76,
            x111 * x144 * x262 + x144 * x265 * x76 + x152 * x262 * x76,
            x102 * x108 * x294 + x108 * x288 * x50 + x137 * x288,
            x108 * x116 * x283 + x108 * x117 * x277 + x116 * x137 * x277,
            x108 * x144 * x273 + x108 * x152 * x268 + x137 * x144 * x268,
            x108 * x173 * x265 + x108 * x180 * x262 + x137 * x173 * x262,
            x102 * x135 * x295 + x11 * x300 + x295 * x297 * x50,
            x106 * x117 * x287 + x115 * x294 * x297 + x116 * x135 * x287,
            x106 * x144 * x283 + x106 * x152 * x277 + x135 * x144 * x277,
            x106 * x173 * x273 + x106 * x180 * x268 + x135 * x173 * x268,
            x106 * x197 * x265 + x106 * x203 * x262 + x135 * x197 * x262,
            x205 * x234 * x34 + x206 * x233 * x34 + x235 * x6 * x91,
            x209 * x233 * x84 + x209 * x234 * x39 + x212 * x233 * x39,
            x205 * x237 * x84 + x205 * x240 * x39 + x206 * x237 * x39,
            x111 * x214 * x233 + x214 * x234 * x76 + x217 * x233 * x76,
            x111 * x209 * x237 + x209 * x240 * x76 + x212 * x237 * x76,
            x111 * x205 * x242 + x205 * x245 * x76 + x206 * x242 * x76,
            x108 * x219 * x234 + x108 * x222 * x233 + x137 * x219 * x233,
            x108 * x214 * x240 + x108 * x217 * x237 + x137 * x214 * x237,
            x108 * x209 * x245 + x108 * x212 * x242 + x137 * x209 * x242,
            x108 * x205 * x250 + x108 * x206 * x247 + x137 * x205 * x247,
            x106 * x227 * x234 + x135 * x227 * x233 + x232 * x297 * x48,
            x106 * x219 * x240 + x106 * x222 * x237 + x135 * x219 * x237,
            x106 * x214 * x245 + x106 * x217 * x242 + x135 * x214 * x242,
            x106 * x209 * x250 + x106 * x212 * x247 + x135 * x209 * x247,
            x106 * x206 * x255 + x135 * x205 * x255 + x260 * x302 * x6,
            x304 * x305 * x9 + x304 * x91 * x99 + x305 * x307,
            x100 * x304 * x84 + x100 * x307 * x39 + x101 * x304 * x39,
            x311 * x39 * x9 + x311 * x84 + x315 * x39 * x99,
            x111 * x120 * x304 + x120 * x307 * x76 + x131 * x304 * x76,
            x100 * x111 * x310 + x100 * x315 * x76 + x101 * x310 * x76,
            x111 * x320 + x320 * x76 * x9 + x325 * x76 * x99,
            x108 * x156 * x307 + x108 * x165 * x304 + x137 * x156 * x304,
            x108 * x120 * x315 + x108 * x131 * x310 + x120 * x137 * x310,
            x100 * x108 * x325 + x100 * x137 * x319 + x101 * x108 * x319,
            x108 * x330 * x9 + x108 * x335 * x99 + x137 * x330,
            x106 * x185 * x307 + x106 * x193 * x304 + x135 * x185 * x304,
            x106 * x156 * x315 + x106 * x165 * x310 + x135 * x156 * x310,
            x106 * x120 * x325 + x106 * x131 * x319 + x120 * x135 * x319,
            x100 * x135 * x329 + x101 * x106 * x329 + x302 * x335 * x93,
            x11 * x339 + x135 * x336 * x99 + x302 * x336 * x9,
            x102 * x340 * x87 + x340 * x341 * x50 + x341 * x342,
            x102 * x26 * x345 + x26 * x344 * x50 + x344 * x82,
            x116 * x26 * x342 + x116 * x340 * x82 + x117 * x26 * x340,
            x102 * x24 * x348 + x24 * x347 * x50 + x347 * x72,
            x116 * x24 * x345 + x116 * x343 * x72 + x117 * x24 * x343,
            x144 * x24 * x342 + x144 * x340 * x72 + x152 * x24 * x340,
            x102 * x349 * x67 + x19 * x351 + x349 * x350 * x50,
            x115 * x348 * x350 + x116 * x346 * x67 + x117 * x346 * x57,
            x144 * x343 * x67 + x144 * x345 * x57 + x152 * x343 * x57,
            x173 * x340 * x67 + x173 * x342 * x57 + x180 * x340 * x57,
            x296
            * (
                x299 * x6
                + x3 * (x227 * x52 + 2 * x230 + 2 * x231 + 4 * x290 + 4 * x292 + 4 * x293)
                + x51 * (2 * a * x352 - 3 * x225 - 3 * x226)
            )
            + x353 * x50
            + x353 * x59,
            x115 * x349 * x354 + x115 * x351 + x117 * x16 * x349,
            x144 * x16 * x348 + x144 * x355 * x59 + x152 * x355,
            x16 * x173 * x345 + x173 * x356 * x59 + x180 * x356,
            x16 * x203 * x340 + x340 * x357 * x59 + x342 * x357,
            x233 * x262 * x87 + x233 * x265 * x32 + x234 * x262 * x32,
            x233 * x26 * x273 + x233 * x268 * x82 + x234 * x26 * x268,
            x237 * x26 * x265 + x237 * x262 * x82 + x240 * x26 * x262,
            x233 * x24 * x283 + x233 * x277 * x72 + x234 * x24 * x277,
            x237 * x24 * x273 + x237 * x268 * x72 + x24 * x240 * x268,
            x24 * x242 * x265 + x24 * x245 * x262 + x242 * x262 * x72,
            x233 * x287 * x67 + x234 * x287 * x57 + x294 * x350 * x48,
            x237 * x277 * x67 + x237 * x283 * x57 + x240 * x277 * x57,
            x242 * x268 * x67 + x242 * x273 * x57 + x245 * x268 * x57,
            x247 * x262 * x67 + x247 * x265 * x57 + x250 * x262 * x57,
            x16 * x234 * x295 + x295 * x354 * x48 + x300 * x48,
            x16 * x237 * x294 + x237 * x358 * x59 + x240 * x358,
            x16 * x242 * x283 + x242 * x359 * x59 + x245 * x359,
            x16 * x247 * x273 + x247 * x360 * x59 + x250 * x360,
            x16 * x260 * x262 + x262 * x361 * x59 + x265 * x361,
            x205 * x304 * x87 + x205 * x307 * x32 + x206 * x304 * x32,
            x209 * x26 * x307 + x209 * x304 * x82 + x212 * x26 * x304,
            x205 * x26 * x315 + x205 * x310 * x82 + x206 * x26 * x310,
            x214 * x24 * x307 + x214 * x304 * x72 + x217 * x24 * x304,
            x209 * x24 * x315 + x209 * x310 * x72 + x212 * x24 * x310,
            x205 * x24 * x325 + x205 * x319 * x72 + x206 * x24 * x319,
            x219 * x304 * x67 + x219 * x307 * x57 + x222 * x304 * x57,
            x214 * x310 * x67 + x214 * x315 * x57 + x217 * x310 * x57,
            x209 * x319 * x67 + x209 * x325 * x57 + x212 * x319 * x57,
            x19 * x335 * x362 + x205 * x329 * x67 + x206 * x329 * x57,
            x16 * x232 * x304 + x304 * x363 * x59 + x307 * x363,
            x16 * x219 * x315 + x219 * x364 * x59 + x222 * x364,
            x16 * x214 * x325 + x214 * x365 * x59 + x217 * x365,
            x16 * x209 * x335 + x209 * x366 * x59 + x212 * x366,
            x16 * x206 * x336 + x336 * x362 * x59 + x339 * x6,
            x367 * x368 * x9 + x367 * x87 * x99 + x368 * x369,
            x100 * x26 * x369 + x100 * x367 * x82 + x101 * x26 * x367,
            x26 * x371 * x9 + x26 * x372 * x99 + x371 * x82,
            x120 * x24 * x369 + x120 * x367 * x72 + x131 * x24 * x367,
            x100 * x24 * x372 + x100 * x370 * x72 + x101 * x24 * x370,
            x24 * x374 * x9 + x24 * x375 * x99 + x374 * x72,
            x156 * x367 * x67 + x156 * x369 * x57 + x165 * x367 * x57,
            x120 * x370 * x67 + x120 * x372 * x57 + x131 * x370 * x57,
            x100 * x373 * x67 + x101 * x373 * x57 + x375 * x376 * x93,
            x19 * x378 + x376 * x377 * x9 + x377 * x67 * x99,
            x16 * x193 * x367 + x367 * x379 * x59 + x369 * x379,
            x156 * x16 * x372 + x156 * x380 * x59 + x165 * x380,
            x120 * x16 * x375 + x120 * x381 * x59 + x131 * x381,
            x101 * x16 * x377 + x301 * x377 * x59 * x93 + x378 * x93,
            x301
            * (
                x3 * (x255 * x52 + 2 * x258 + 2 * x259 + 4 * x332 + 4 * x333 + 4 * x334)
                + x338 * x48
                + x51 * (2 * a * x382 - 3 * x253 - 3 * x254)
            )
            + x383 * x59
            + x383 * x9,
        ]
    )
    return S


def kinetic3d_40(a, A, b, B):
    """Cartesian 3D (gs) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = (2 * b + x1) ** (-1.0)
    x3 = (a + b) ** (-1.0)
    x4 = x3 * (a * A[1] + b * B[1]) - A[1]
    x5 = x4**2
    x6 = 2 * a**2
    x7 = -x0 - x6 * (x2 + x5)
    x8 = b * x3
    x9 = a * x8
    x10 = numpy.exp(-x9 * (A[0] - B[0]) ** 2)
    x11 = numpy.sqrt(numpy.pi) * numpy.sqrt(x3)
    x12 = x10 * x11
    x13 = x12 * x2
    x14 = x3 * (a * A[0] + b * B[0]) - A[0]
    x15 = x14**2
    x16 = x12 * x15
    x17 = x13 + x16
    x18 = 2 * x13 * x14 + x14 * x17
    x19 = x14 * x18 + x2 * (3 * x13 + 3 * x16)
    x20 = numpy.exp(-x9 * (A[1] - B[1]) ** 2)
    x21 = numpy.exp(-x9 * (A[2] - B[2]) ** 2)
    x22 = numpy.pi * x21 * x3
    x23 = x20 * x22
    x24 = x19 * x23
    x25 = x3 * (a * A[2] + b * B[2]) - A[2]
    x26 = x25**2
    x27 = -x0 - x6 * (x2 + x26)
    x28 = -x0 - x6 * (x15 + x2)
    x29 = x13 * x28
    x30 = 2 * x12
    x31 = x8 * (x1 * x17 - x30)
    x32 = x12 * x14
    x33 = x1 * x8
    x34 = x28 * x32 + x32 * x33
    x35 = x14 * x34
    x36 = 4 * x9
    x37 = x29 + x31 + x35
    x38 = x14 * x37 + x2 * (x14 * x28 * x30 + x32 * x36) + x8 * (x1 * x18 - 3 * x32)
    x39 = x18 * x23
    x40 = x11 * x20
    x41 = x4 * x40
    x42 = x33 * x41 + x41 * x7
    x43 = x11 * x21
    x44 = x42 * x43
    x45 = x23 * x38
    x46 = x25 * x43
    x47 = x27 * x46 + x33 * x46
    x48 = x2 * x40
    x49 = x40 * x5
    x50 = x48 + x49
    x51 = x17 * x43
    x52 = x48 * x7
    x53 = 2 * x40
    x54 = x8 * (x1 * x50 - x53)
    x55 = x4 * x42
    x56 = x52 + x54 + x55
    x57 = x2 * x43
    x58 = x26 * x43
    x59 = x57 + x58
    x60 = x17 * x40
    x61 = x27 * x57
    x62 = 2 * x43
    x63 = x8 * (x1 * x59 - x62)
    x64 = x25 * x47
    x65 = x61 + x63 + x64
    x66 = 2 * x4 * x48 + x4 * x50
    x67 = x10 * x22
    x68 = x14 * x67
    x69 = x2 * (x36 * x41 + x4 * x53 * x7) + x4 * x56 + x8 * (x1 * x66 - 3 * x41)
    x70 = x67 * x69
    x71 = numpy.pi * x10 * x20 * x3
    x72 = x14 * x71
    x73 = 2 * x25 * x57 + x25 * x59
    x74 = x2 * (x25 * x27 * x62 + x36 * x46) + x25 * x65 + x8 * (x1 * x73 - 3 * x46)
    x75 = x71 * x74
    x76 = x2 * (3 * x48 + 3 * x49) + x4 * x66
    x77 = x67 * x76
    x78 = x12 * x50
    x79 = x2 * (3 * x57 + 3 * x58) + x25 * x73
    x80 = x71 * x79

    # 15 item(s)
    S = numpy.array(
        [
            x23
            * (
                x14 * x38
                + x2 * (3 * x29 + 3 * x31 + 3 * x35)
                + x8 * (2 * a * x19 - 4 * x13 - 4 * x16)
            )
            + x24 * x27
            + x24 * x7,
            x18 * x44 + x27 * x39 * x4 + x4 * x45,
            x18 * x40 * x47 + x25 * x39 * x7 + x25 * x45,
            x27 * x50 * x51 + x37 * x43 * x50 + x51 * x56,
            x17 * x25 * x44 + x17 * x41 * x47 + x23 * x25 * x37 * x4,
            x37 * x40 * x59 + x59 * x60 * x7 + x60 * x65,
            x14 * x70 + x27 * x66 * x68 + x34 * x43 * x66,
            x25 * x56 * x68 + x32 * x47 * x50 + x34 * x46 * x50,
            x32 * x42 * x59 + x34 * x41 * x59 + x4 * x65 * x72,
            x14 * x75 + x34 * x40 * x73 + x7 * x72 * x73,
            x27 * x77
            + x28 * x77
            + x67
            * (
                x2 * (3 * x52 + 3 * x54 + 3 * x55)
                + x4 * x69
                + x8 * (2 * a * x76 - 4 * x48 - 4 * x49)
            ),
            x12 * x47 * x66 + x25 * x28 * x66 * x67 + x25 * x70,
            x12 * x56 * x59 + x28 * x59 * x78 + x65 * x78,
            x12 * x42 * x73 + x28 * x4 * x71 * x73 + x4 * x75,
            x28 * x80
            + x7 * x80
            + x71
            * (
                x2 * (3 * x61 + 3 * x63 + 3 * x64)
                + x25 * x74
                + x8 * (2 * a * x79 - 4 * x57 - 4 * x58)
            ),
        ]
    )
    return S


def kinetic3d_41(a, A, b, B):
    """Cartesian 3D (gp) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = (2 * b + x1) ** (-1.0)
    x3 = (a + b) ** (-1.0)
    x4 = -x3 * (a * A[1] + b * B[1])
    x5 = -x4 - A[1]
    x6 = x5**2
    x7 = 2 * a**2
    x8 = -x0 - x7 * (x2 + x6)
    x9 = -x3 * (a * A[0] + b * B[0])
    x10 = -x9 - A[0]
    x11 = b * x3
    x12 = a * x11
    x13 = numpy.exp(-x12 * (A[0] - B[0]) ** 2)
    x14 = numpy.sqrt(numpy.pi) * numpy.sqrt(x3)
    x15 = x13 * x14
    x16 = x10 * x15
    x17 = -x9 - B[0]
    x18 = x15 * x17
    x19 = x2 * (x16 + x18)
    x20 = x15 * x2
    x21 = x16 * x17
    x22 = x20 + x21
    x23 = x10 * x22
    x24 = x10**2
    x25 = x15 * x24
    x26 = x20 + x25
    x27 = 2 * x10 * x20 + x10 * x26
    x28 = 3 * x20
    x29 = 2 * x15
    x30 = x17 * x29
    x31 = x19 + x23
    x32 = x10 * x31 + x2 * (x10 * x30 + x25 + x28)
    x33 = x10 * x32 + x2 * (3 * x19 + 3 * x23 + x27)
    x34 = numpy.exp(-x12 * (A[1] - B[1]) ** 2)
    x35 = numpy.exp(-x12 * (A[2] - B[2]) ** 2)
    x36 = numpy.pi * x3 * x35
    x37 = x34 * x36
    x38 = x33 * x37
    x39 = -x3 * (a * A[2] + b * B[2])
    x40 = -x39 - A[2]
    x41 = x40**2
    x42 = -x0 - x7 * (x2 + x41)
    x43 = x1 * x11
    x44 = -x0 - x7 * (x2 + x24)
    x45 = x16 * x43 + x16 * x44
    x46 = x18 * x43 + x18 * x44
    x47 = x2 * (x45 + x46)
    x48 = x20 * x44
    x49 = x10 * x46
    x50 = x22 * x43 + x48 + x49
    x51 = x10 * x50
    x52 = x11 * (x1 * x31 - x30)
    x53 = 4 * x12
    x54 = x10 * x45
    x55 = x11 * (x1 * x26 - x29)
    x56 = x54 + x55
    x57 = x48 + x56
    x58 = x10 * x57 + x11 * (x1 * x27 - 3 * x16) + x2 * (x10 * x29 * x44 + x16 * x53)
    x59 = 3 * x48
    x60 = x47 + x51 + x52
    x61 = (
        x10 * x60
        + x11 * (2 * a * x32 - 3 * x21 - x28)
        + x2 * (x22 * x53 + 2 * x49 + x56 + x59)
    )
    x62 = -x4 - B[1]
    x63 = x10 * x27 + x2 * (3 * x25 + x28)
    x64 = x37 * x63
    x65 = x14 * x34
    x66 = x62 * x65
    x67 = x43 * x66 + x66 * x8
    x68 = x14 * x35
    x69 = x37 * (
        x10 * x58
        + x11 * (2 * a * x63 - 4 * x20 - 4 * x25)
        + x2 * (3 * x54 + 3 * x55 + x59)
    )
    x70 = -x39 - B[2]
    x71 = x68 * x70
    x72 = x42 * x71 + x43 * x71
    x73 = x37 * x5
    x74 = x5 * x65
    x75 = x43 * x74 + x74 * x8
    x76 = x37 * x61
    x77 = x14 * x2
    x78 = x34 * x77
    x79 = x62 * x74
    x80 = x78 + x79
    x81 = x27 * x68
    x82 = x78 * x8
    x83 = x5 * x67
    x84 = x43 * x80 + x82 + x83
    x85 = x37 * x40
    x86 = x40 * x68
    x87 = x42 * x86 + x43 * x86
    x88 = x35 * x77
    x89 = x70 * x86
    x90 = x88 + x89
    x91 = x27 * x65
    x92 = x42 * x88
    x93 = x40 * x72
    x94 = x43 * x90 + x92 + x93
    x95 = x6 * x65
    x96 = x78 + x95
    x97 = x31 * x68
    x98 = x5 * x75
    x99 = 2 * x65
    x100 = x11 * (x1 * x96 - x99)
    x101 = x100 + x98
    x102 = x101 + x82
    x103 = x2 * (x66 + x74)
    x104 = x5 * x80
    x105 = x103 + x104
    x106 = x105 * x68
    x107 = x62 * x99
    x108 = x11 * (x1 * x105 - x107)
    x109 = x2 * (x67 + x75)
    x110 = x5 * x84
    x111 = x108 + x109 + x110
    x112 = x41 * x68
    x113 = x112 + x88
    x114 = x31 * x65
    x115 = x40 * x87
    x116 = 2 * x68
    x117 = x11 * (x1 * x113 - x116)
    x118 = x115 + x117
    x119 = x118 + x92
    x120 = x2 * (x71 + x86)
    x121 = x40 * x90
    x122 = x120 + x121
    x123 = x122 * x65
    x124 = x116 * x70
    x125 = x11 * (x1 * x122 - x124)
    x126 = x2 * (x72 + x87)
    x127 = x40 * x94
    x128 = x125 + x126 + x127
    x129 = 2 * x5 * x78 + x5 * x96
    x130 = x129 * x68
    x131 = x102 * x5 + x11 * (x1 * x129 - 3 * x74) + x2 * (x5 * x8 * x99 + x53 * x74)
    x132 = 3 * x78
    x133 = x105 * x5 + x2 * (x107 * x5 + x132 + x95)
    x134 = x13 * x36
    x135 = x10 * x134
    x136 = 3 * x82
    x137 = (
        x11 * (2 * a * x133 - x132 - 3 * x79)
        + x111 * x5
        + x2 * (x101 + x136 + x53 * x80 + 2 * x83)
    )
    x138 = x134 * x137
    x139 = numpy.pi * x13 * x3 * x34
    x140 = x10 * x139
    x141 = x113 * x40 + 2 * x40 * x88
    x142 = x141 * x65
    x143 = x11 * (x1 * x141 - 3 * x86) + x119 * x40 + x2 * (x116 * x40 * x42 + x53 * x86)
    x144 = 3 * x88
    x145 = x122 * x40 + x2 * (x112 + x124 * x40 + x144)
    x146 = 3 * x92
    x147 = (
        x11 * (2 * a * x145 - x144 - 3 * x89)
        + x128 * x40
        + x2 * (x118 + x146 + x53 * x90 + 2 * x93)
    )
    x148 = x139 * x147
    x149 = x129 * x5 + x2 * (x132 + 3 * x95)
    x150 = x134 * x149
    x151 = x134 * (
        x11 * (2 * a * x149 - 4 * x78 - 4 * x95)
        + x131 * x5
        + x2 * (3 * x100 + x136 + 3 * x98)
    )
    x152 = x133 * x5 + x2 * (3 * x103 + 3 * x104 + x129)
    x153 = x134 * x152
    x154 = x134 * x40
    x155 = x129 * x15
    x156 = x105 * x15
    x157 = x122 * x15
    x158 = x139 * x5
    x159 = x141 * x15
    x160 = x141 * x40 + x2 * (3 * x112 + x144)
    x161 = x139 * x160
    x162 = x139 * (
        x11 * (2 * a * x160 - 4 * x112 - 4 * x88)
        + x143 * x40
        + x2 * (3 * x115 + 3 * x117 + x146)
    )
    x163 = x145 * x40 + x2 * (3 * x120 + 3 * x121 + x141)
    x164 = x139 * x163

    # 45 item(s)
    S = numpy.array(
        [
            x37
            * (
                x10 * x61
                + x11 * (2 * a * x33 - 4 * x19 - 4 * x23)
                + x2 * (3 * x47 + 3 * x51 + 3 * x52 + x58)
            )
            + x38 * x42
            + x38 * x8,
            x42 * x62 * x64 + x62 * x69 + x63 * x67 * x68,
            x63 * x65 * x72 + x64 * x70 * x8 + x69 * x70,
            x32 * x42 * x73 + x32 * x68 * x75 + x5 * x76,
            x42 * x80 * x81 + x58 * x68 * x80 + x81 * x84,
            x27 * x71 * x75 + x27 * x72 * x74 + x58 * x70 * x73,
            x32 * x65 * x87 + x32 * x8 * x85 + x40 * x76,
            x27 * x66 * x87 + x27 * x67 * x86 + x58 * x62 * x85,
            x58 * x65 * x90 + x8 * x90 * x91 + x91 * x94,
            x102 * x97 + x42 * x96 * x97 + x60 * x68 * x96,
            x106 * x26 * x42 + x106 * x57 + x111 * x26 * x68,
            x102 * x26 * x71 + x26 * x72 * x96 + x57 * x71 * x96,
            x31 * x74 * x87 + x31 * x75 * x86 + x40 * x60 * x73,
            x26 * x80 * x87 + x26 * x84 * x86 + x57 * x80 * x86,
            x26 * x74 * x94 + x26 * x75 * x90 + x57 * x74 * x90,
            x113 * x114 * x8 + x113 * x60 * x65 + x114 * x119,
            x113 * x26 * x67 + x113 * x57 * x66 + x119 * x26 * x66,
            x123 * x26 * x8 + x123 * x57 + x128 * x26 * x65,
            x130 * x22 * x42 + x130 * x50 + x131 * x22 * x68,
            x10 * x138 + x133 * x135 * x42 + x133 * x45 * x68,
            x129 * x16 * x72 + x129 * x45 * x71 + x131 * x135 * x70,
            x102 * x22 * x86 + x22 * x87 * x96 + x50 * x86 * x96,
            x105 * x16 * x87 + x105 * x45 * x86 + x111 * x135 * x40,
            x102 * x16 * x90 + x16 * x94 * x96 + x45 * x90 * x96,
            x113 * x22 * x75 + x113 * x50 * x74 + x119 * x22 * x74,
            x113 * x16 * x84 + x113 * x45 * x80 + x119 * x16 * x80,
            x122 * x16 * x75 + x122 * x45 * x74 + x128 * x140 * x5,
            x142 * x22 * x8 + x142 * x50 + x143 * x22 * x65,
            x140 * x143 * x62 + x141 * x16 * x67 + x141 * x45 * x66,
            x10 * x148 + x140 * x145 * x8 + x145 * x45 * x65,
            x149 * x46 * x68 + x150 * x17 * x42 + x151 * x17,
            x134
            * (
                x11 * (2 * a * x152 - 4 * x103 - 4 * x104)
                + x137 * x5
                + x2 * (3 * x108 + 3 * x109 + 3 * x110 + x131)
            )
            + x153 * x42
            + x153 * x44,
            x149 * x15 * x72 + x150 * x44 * x70 + x151 * x70,
            x129 * x18 * x87 + x129 * x46 * x86 + x131 * x154 * x17,
            x133 * x15 * x87 + x133 * x154 * x44 + x138 * x40,
            x131 * x15 * x90 + x155 * x44 * x90 + x155 * x94,
            x102 * x113 * x18 + x113 * x46 * x96 + x119 * x18 * x96,
            x111 * x113 * x15 + x113 * x156 * x44 + x119 * x156,
            x102 * x157 + x128 * x15 * x96 + x157 * x44 * x96,
            x141 * x18 * x75 + x141 * x46 * x74 + x143 * x158 * x17,
            x143 * x15 * x80 + x159 * x44 * x80 + x159 * x84,
            x145 * x15 * x75 + x145 * x158 * x44 + x148 * x5,
            x160 * x46 * x65 + x161 * x17 * x8 + x162 * x17,
            x15 * x160 * x67 + x161 * x44 * x62 + x162 * x62,
            x139
            * (
                x11 * (2 * a * x163 - 4 * x120 - 4 * x121)
                + x147 * x40
                + x2 * (3 * x125 + 3 * x126 + 3 * x127 + x143)
            )
            + x164 * x44
            + x164 * x8,
        ]
    )
    return S


def kinetic3d_42(a, A, b, B):
    """Cartesian 3D (gd) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = -x5 - A[1]
    x7 = x6**2
    x8 = 2 * a**2
    x9 = -x0 - x8 * (x3 + x7)
    x10 = -x4 * (a * A[0] + b * B[0])
    x11 = -x10 - A[0]
    x12 = x11**2
    x13 = b * x4
    x14 = a * x13
    x15 = numpy.exp(-x14 * (A[0] - B[0]) ** 2)
    x16 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x17 = x15 * x16
    x18 = x12 * x17
    x19 = x17 * x3
    x20 = 3 * x19
    x21 = -x10 - B[0]
    x22 = 2 * x17
    x23 = x21 * x22
    x24 = x11 * x23 + x20
    x25 = x3 * (x18 + x24)
    x26 = x21**2
    x27 = x17 * x26
    x28 = x3 * (x24 + x27)
    x29 = 2 * x19
    x30 = x19 + x27
    x31 = x11 * x30
    x32 = x21 * x29 + x31
    x33 = x11 * x32
    x34 = x11 * x17
    x35 = x17 * x21
    x36 = x3 * (x34 + x35)
    x37 = x21 * x34
    x38 = x19 + x37
    x39 = x11 * x38
    x40 = x36 + x39
    x41 = x11 * x40
    x42 = x28 + x33
    x43 = 4 * x19
    x44 = x11 * x42 + x3 * (x21 * x43 + 2 * x31 + 2 * x36 + 2 * x39)
    x45 = x11 * x44 + x3 * (2 * x25 + 3 * x28 + 3 * x33 + 2 * x41)
    x46 = numpy.exp(-x14 * (A[1] - B[1]) ** 2)
    x47 = numpy.exp(-x14 * (A[2] - B[2]) ** 2)
    x48 = numpy.pi * x4 * x47
    x49 = x46 * x48
    x50 = x45 * x49
    x51 = -x4 * (a * A[2] + b * B[2])
    x52 = -x51 - A[2]
    x53 = x52**2
    x54 = -x0 - x8 * (x3 + x53)
    x55 = x13 * (2 * a * x42 - x22 * x26 - x29)
    x56 = x25 + x41
    x57 = 2 * a * x56 - x20 - 3 * x37
    x58 = x2 * x4
    x59 = x1 * x13
    x60 = -x0 - x8 * (x12 + x3)
    x61 = x35 * x59 + x35 * x60
    x62 = -x22
    x63 = a * x4
    x64 = x21 * x61 + x63 * (x2 * x30 + x62)
    x65 = x19 * x60
    x66 = 3 * x65
    x67 = x11 * x61
    x68 = 4 * x14
    x69 = x38 * x68 + x66 + 2 * x67
    x70 = x3 * (x64 + x69)
    x71 = x34 * x59 + x34 * x60
    x72 = x11 * x71
    x73 = x18 + x19
    x74 = x13 * (x1 * x73 + x62)
    x75 = x72 + x74
    x76 = x3 * (x69 + x75)
    x77 = x3 * (x23 * x60 + x35 * x68)
    x78 = x64 + x65
    x79 = x11 * x78
    x80 = x32 * x59 + x77 + x79
    x81 = x11 * x80
    x82 = x1 * x40 - x23
    x83 = x13 * x82
    x84 = x3 * (x61 + x71)
    x85 = x38 * x59 + x65 + x67
    x86 = x11 * x85
    x87 = x83 + x84 + x86
    x88 = x11 * x87
    x89 = x55 + x70 + x81
    x90 = (
        x11 * x89
        + x13 * (2 * a * x44 - 6 * x19 * x21 - 3 * x31)
        + x3 * (x32 * x68 + x58 * x82 + 2 * x77 + 2 * x79 + 2 * x84 + 2 * x86)
    )
    x91 = -x5 - B[1]
    x92 = x11 * x29 + x11 * x73
    x93 = x11 * x56 + x3 * (3 * x36 + 3 * x39 + x92)
    x94 = x49 * x93
    x95 = x16 * x46
    x96 = x91 * x95
    x97 = x59 * x96 + x9 * x96
    x98 = x16 * x47
    x99 = x65 + x75
    x100 = x11 * x99 + x13 * (x1 * x92 - 3 * x34) + x3 * (x11 * x22 * x60 + x34 * x68)
    x101 = x13 * x57 + x76 + x88
    x102 = x49 * (
        x101 * x11
        + x13 * (2 * a * x93 - 4 * x36 - 4 * x39)
        + x3 * (x100 + 3 * x83 + 3 * x84 + 3 * x86)
    )
    x103 = -x51 - B[2]
    x104 = x103 * x98
    x105 = x104 * x54 + x104 * x59
    x106 = x3 * x95
    x107 = x91**2 * x95
    x108 = x106 + x107
    x109 = x11 * x92 + x3 * (3 * x18 + x20)
    x110 = x109 * x98
    x111 = x106 * x9
    x112 = 2 * x95
    x113 = -x112
    x114 = x63 * (x108 * x2 + x113) + x91 * x97
    x115 = x111 + x114
    x116 = (
        x100 * x11 + x13 * (2 * a * x109 - 4 * x18 - x43) + x3 * (x66 + 3 * x72 + 3 * x74)
    )
    x117 = x103 * x49
    x118 = x3 * x98
    x119 = x103**2 * x98
    x120 = x118 + x119
    x121 = x109 * x95
    x122 = x118 * x54
    x123 = 2 * x98
    x124 = -x123
    x125 = x103 * x105 + x63 * (x120 * x2 + x124)
    x126 = x122 + x125
    x127 = x44 * x49
    x128 = x6 * x95
    x129 = x128 * x59 + x128 * x9
    x130 = x49 * x90
    x131 = x128 * x91
    x132 = x106 + x131
    x133 = x56 * x98
    x134 = x6 * x97
    x135 = x111 + x132 * x59 + x134
    x136 = 2 * x106
    x137 = x108 * x6
    x138 = x136 * x91 + x137
    x139 = x92 * x98
    x140 = x112 * x91
    x141 = x3 * (x140 * x9 + x68 * x96)
    x142 = x115 * x6
    x143 = x138 * x59 + x141 + x142
    x144 = x52 * x98
    x145 = x144 * x54 + x144 * x59
    x146 = x49 * x52
    x147 = x103 * x144
    x148 = x118 + x147
    x149 = x56 * x95
    x150 = x105 * x52
    x151 = x122 + x148 * x59 + x150
    x152 = 2 * x118
    x153 = x120 * x52
    x154 = x103 * x152 + x153
    x155 = x92 * x95
    x156 = x103 * x123
    x157 = x3 * (x104 * x68 + x156 * x54)
    x158 = x126 * x52
    x159 = x154 * x59 + x157 + x158
    x160 = x7 * x95
    x161 = x106 + x160
    x162 = x42 * x98
    x163 = x129 * x6
    x164 = x13 * (x1 * x161 + x113)
    x165 = x163 + x164
    x166 = x111 + x165
    x167 = x3 * (x128 + x96)
    x168 = x132 * x6
    x169 = x167 + x168
    x170 = x40 * x98
    x171 = x1 * x169 - x140
    x172 = x13 * x171
    x173 = x3 * (x129 + x97)
    x174 = x135 * x6
    x175 = x172 + x173 + x174
    x176 = 3 * x106
    x177 = x140 * x6 + x176
    x178 = x3 * (x107 + x177)
    x179 = x138 * x6
    x180 = x178 + x179
    x181 = x180 * x98
    x182 = x13 * (2 * a * x180 - 2 * x107 - x136)
    x183 = 3 * x111
    x184 = x132 * x68 + 2 * x134 + x183
    x185 = x3 * (x114 + x184)
    x186 = x143 * x6
    x187 = x182 + x185 + x186
    x188 = x53 * x98
    x189 = x118 + x188
    x190 = x42 * x95
    x191 = x145 * x52
    x192 = x13 * (x1 * x189 + x124)
    x193 = x191 + x192
    x194 = x122 + x193
    x195 = x3 * (x104 + x144)
    x196 = x148 * x52
    x197 = x195 + x196
    x198 = x40 * x95
    x199 = x1 * x197 - x156
    x200 = x13 * x199
    x201 = x3 * (x105 + x145)
    x202 = x151 * x52
    x203 = x200 + x201 + x202
    x204 = 3 * x118
    x205 = x156 * x52 + x204
    x206 = x3 * (x119 + x205)
    x207 = x154 * x52
    x208 = x206 + x207
    x209 = x208 * x95
    x210 = x13 * (2 * a * x208 - 2 * x119 - x152)
    x211 = 3 * x122
    x212 = x148 * x68 + 2 * x150 + x211
    x213 = x3 * (x125 + x212)
    x214 = x159 * x52
    x215 = x210 + x213 + x214
    x216 = x136 * x6 + x161 * x6
    x217 = x32 * x98
    x218 = x13 * (x1 * x216 - 3 * x128) + x166 * x6 + x3 * (x112 * x6 * x9 + x128 * x68)
    x219 = x3 * (x160 + x177)
    x220 = x169 * x6
    x221 = x219 + x220
    x222 = x221 * x98
    x223 = 2 * a * x221 - 3 * x131 - x176
    x224 = x3 * (x165 + x184)
    x225 = x175 * x6
    x226 = x13 * x223 + x224 + x225
    x227 = 4 * x106
    x228 = x180 * x6 + x3 * (2 * x137 + 2 * x167 + 2 * x168 + x227 * x91)
    x229 = x15 * x48
    x230 = x11 * x229
    x231 = (
        x13 * (2 * a * x228 - 6 * x106 * x91 - 3 * x137)
        + x187 * x6
        + x3 * (x138 * x68 + 2 * x141 + 2 * x142 + x171 * x58 + 2 * x173 + 2 * x174)
    )
    x232 = x229 * x231
    x233 = numpy.pi * x15 * x4 * x46
    x234 = x11 * x233
    x235 = x152 * x52 + x189 * x52
    x236 = x32 * x95
    x237 = (
        x13 * (x1 * x235 - 3 * x144) + x194 * x52 + x3 * (x123 * x52 * x54 + x144 * x68)
    )
    x238 = x3 * (x188 + x205)
    x239 = x197 * x52
    x240 = x238 + x239
    x241 = x240 * x95
    x242 = 2 * a * x240 - 3 * x147 - x204
    x243 = x3 * (x193 + x212)
    x244 = x203 * x52
    x245 = x13 * x242 + x243 + x244
    x246 = 4 * x118
    x247 = x208 * x52 + x3 * (x103 * x246 + 2 * x153 + 2 * x195 + 2 * x196)
    x248 = (
        x13 * (2 * a * x247 - 6 * x103 * x118 - 3 * x153)
        + x215 * x52
        + x3 * (x154 * x68 + 2 * x157 + 2 * x158 + x199 * x58 + 2 * x201 + 2 * x202)
    )
    x249 = x233 * x248
    x250 = x216 * x6 + x3 * (3 * x160 + x176)
    x251 = x250 * x98
    x252 = (
        x13 * (2 * a * x250 - 4 * x160 - x227)
        + x218 * x6
        + x3 * (3 * x163 + 3 * x164 + x183)
    )
    x253 = x221 * x6 + x3 * (3 * x167 + 3 * x168 + x216)
    x254 = x21 * x229
    x255 = x229 * (
        x13 * (2 * a * x253 - 4 * x167 - 4 * x168)
        + x226 * x6
        + x3 * (3 * x172 + 3 * x173 + 3 * x174 + x218)
    )
    x256 = x228 * x6 + x3 * (3 * x178 + 3 * x179 + 2 * x219 + 2 * x220)
    x257 = x229 * x256
    x258 = x229 * x60
    x259 = x17 * x250
    x260 = x17 * x221
    x261 = x17 * x216
    x262 = x17 * x180
    x263 = x169 * x17
    x264 = x17 * x208
    x265 = x233 * x6
    x266 = x138 * x17
    x267 = x17 * x240
    x268 = x235 * x52 + x3 * (3 * x188 + x204)
    x269 = x268 * x95
    x270 = (
        x13 * (2 * a * x268 - 4 * x188 - x246)
        + x237 * x52
        + x3 * (3 * x191 + 3 * x192 + x211)
    )
    x271 = x21 * x233
    x272 = x240 * x52 + x3 * (3 * x195 + 3 * x196 + x235)
    x273 = x233 * (
        x13 * (2 * a * x272 - 4 * x195 - 4 * x196)
        + x245 * x52
        + x3 * (3 * x200 + 3 * x201 + 3 * x202 + x237)
    )
    x274 = x17 * x268
    x275 = x247 * x52 + x3 * (3 * x206 + 3 * x207 + 2 * x238 + 2 * x239)
    x276 = x233 * x275

    # 90 item(s)
    S = numpy.array(
        [
            x49
            * (
                x11 * x90
                + x13 * (2 * a * x45 - 4 * x28 - 4 * x33)
                + x3 * (3 * x55 + x57 * x58 + 3 * x70 + 2 * x76 + 3 * x81 + 2 * x88)
            )
            + x50 * x54
            + x50 * x9,
            x102 * x91 + x54 * x91 * x94 + x93 * x97 * x98,
            x102 * x103 + x103 * x9 * x94 + x105 * x93 * x95,
            x108 * x110 * x54 + x108 * x116 * x98 + x110 * x115,
            x104 * x109 * x97 + x105 * x109 * x96 + x116 * x117 * x91,
            x116 * x120 * x95 + x120 * x121 * x9 + x121 * x126,
            x127 * x54 * x6 + x129 * x44 * x98 + x130 * x6,
            x101 * x132 * x98 + x132 * x133 * x54 + x133 * x135,
            x101 * x117 * x6 + x104 * x129 * x56 + x105 * x128 * x56,
            x100 * x138 * x98 + x138 * x139 * x54 + x139 * x143,
            x100 * x104 * x132 + x104 * x135 * x92 + x105 * x132 * x92,
            x100 * x120 * x128 + x120 * x129 * x92 + x126 * x128 * x92,
            x127 * x52 * x9 + x130 * x52 + x145 * x44 * x95,
            x101 * x146 * x91 + x144 * x56 * x97 + x145 * x56 * x96,
            x101 * x148 * x95 + x148 * x149 * x9 + x149 * x151,
            x100 * x108 * x144 + x108 * x145 * x92 + x115 * x144 * x92,
            x100 * x148 * x96 + x148 * x92 * x97 + x151 * x92 * x96,
            x100 * x154 * x95 + x154 * x155 * x9 + x155 * x159,
            x161 * x162 * x54 + x161 * x89 * x98 + x162 * x166,
            x169 * x170 * x54 + x169 * x87 * x98 + x170 * x175,
            x104 * x161 * x87 + x104 * x166 * x40 + x105 * x161 * x40,
            x181 * x54 * x73 + x181 * x99 + x187 * x73 * x98,
            x104 * x169 * x99 + x104 * x175 * x73 + x105 * x169 * x73,
            x120 * x161 * x99 + x120 * x166 * x73 + x126 * x161 * x73,
            x128 * x145 * x42 + x129 * x144 * x42 + x146 * x6 * x89,
            x132 * x144 * x87 + x132 * x145 * x40 + x135 * x144 * x40,
            x128 * x148 * x87 + x128 * x151 * x40 + x129 * x148 * x40,
            x138 * x144 * x99 + x138 * x145 * x73 + x143 * x144 * x73,
            x132 * x148 * x99 + x132 * x151 * x73 + x135 * x148 * x73,
            x128 * x154 * x99 + x128 * x159 * x73 + x129 * x154 * x73,
            x189 * x190 * x9 + x189 * x89 * x95 + x190 * x194,
            x189 * x40 * x97 + x189 * x87 * x96 + x194 * x40 * x96,
            x197 * x198 * x9 + x197 * x87 * x95 + x198 * x203,
            x108 * x189 * x99 + x108 * x194 * x73 + x115 * x189 * x73,
            x197 * x73 * x97 + x197 * x96 * x99 + x203 * x73 * x96,
            x209 * x73 * x9 + x209 * x99 + x215 * x73 * x95,
            x216 * x217 * x54 + x216 * x80 * x98 + x217 * x218,
            x222 * x38 * x54 + x222 * x85 + x226 * x38 * x98,
            x104 * x216 * x85 + x104 * x218 * x38 + x105 * x216 * x38,
            x11 * x232 + x228 * x230 * x54 + x228 * x71 * x98,
            x103 * x226 * x230 + x104 * x221 * x71 + x105 * x221 * x34,
            x120 * x216 * x71 + x120 * x218 * x34 + x126 * x216 * x34,
            x144 * x161 * x80 + x144 * x166 * x32 + x145 * x161 * x32,
            x144 * x169 * x85 + x144 * x175 * x38 + x145 * x169 * x38,
            x148 * x161 * x85 + x148 * x166 * x38 + x151 * x161 * x38,
            x144 * x180 * x71 + x145 * x180 * x34 + x187 * x230 * x52,
            x148 * x169 * x71 + x148 * x175 * x34 + x151 * x169 * x34,
            x154 * x161 * x71 + x154 * x166 * x34 + x159 * x161 * x34,
            x128 * x189 * x80 + x128 * x194 * x32 + x129 * x189 * x32,
            x132 * x189 * x85 + x132 * x194 * x38 + x135 * x189 * x38,
            x128 * x197 * x85 + x128 * x203 * x38 + x129 * x197 * x38,
            x138 * x189 * x71 + x138 * x194 * x34 + x143 * x189 * x34,
            x132 * x197 * x71 + x132 * x203 * x34 + x135 * x197 * x34,
            x128 * x208 * x71 + x129 * x208 * x34 + x215 * x234 * x6,
            x235 * x236 * x9 + x235 * x80 * x95 + x236 * x237,
            x235 * x38 * x97 + x235 * x85 * x96 + x237 * x38 * x96,
            x241 * x38 * x9 + x241 * x85 + x245 * x38 * x95,
            x108 * x235 * x71 + x108 * x237 * x34 + x115 * x235 * x34,
            x234 * x245 * x91 + x240 * x34 * x97 + x240 * x71 * x96,
            x11 * x249 + x234 * x247 * x9 + x247 * x71 * x95,
            x251 * x30 * x54 + x251 * x78 + x252 * x30 * x98,
            x21 * x255 + x253 * x254 * x54 + x253 * x61 * x98,
            x103 * x252 * x254 + x104 * x250 * x61 + x105 * x250 * x35,
            x229
            * (
                x13 * (2 * a * x256 - 4 * x178 - 4 * x179)
                + x231 * x6
                + x3 * (3 * x182 + 3 * x185 + 3 * x186 + x223 * x58 + 2 * x224 + 2 * x225)
            )
            + x257 * x54
            + x257 * x60,
            x103 * x253 * x258 + x103 * x255 + x105 * x17 * x253,
            x120 * x17 * x252 + x120 * x259 * x60 + x126 * x259,
            x144 * x216 * x78 + x144 * x218 * x30 + x145 * x216 * x30,
            x144 * x221 * x61 + x145 * x221 * x35 + x226 * x254 * x52,
            x148 * x216 * x61 + x148 * x218 * x35 + x151 * x216 * x35,
            x145 * x17 * x228 + x228 * x258 * x52 + x232 * x52,
            x148 * x17 * x226 + x148 * x260 * x60 + x151 * x260,
            x154 * x17 * x218 + x154 * x261 * x60 + x159 * x261,
            x161 * x189 * x78 + x161 * x194 * x30 + x166 * x189 * x30,
            x169 * x189 * x61 + x169 * x194 * x35 + x175 * x189 * x35,
            x161 * x197 * x61 + x161 * x203 * x35 + x166 * x197 * x35,
            x17 * x187 * x189 + x189 * x262 * x60 + x194 * x262,
            x17 * x175 * x197 + x197 * x263 * x60 + x203 * x263,
            x161 * x17 * x215 + x161 * x264 * x60 + x166 * x264,
            x128 * x235 * x78 + x128 * x237 * x30 + x129 * x235 * x30,
            x132 * x235 * x61 + x132 * x237 * x35 + x135 * x235 * x35,
            x128 * x240 * x61 + x129 * x240 * x35 + x21 * x245 * x265,
            x143 * x17 * x235 + x235 * x266 * x60 + x237 * x266,
            x132 * x17 * x245 + x132 * x267 * x60 + x135 * x267,
            x129 * x17 * x247 + x247 * x265 * x60 + x249 * x6,
            x269 * x30 * x9 + x269 * x78 + x270 * x30 * x95,
            x268 * x35 * x97 + x268 * x61 * x96 + x270 * x271 * x91,
            x21 * x273 + x271 * x272 * x9 + x272 * x61 * x95,
            x108 * x17 * x270 + x108 * x274 * x60 + x115 * x274,
            x17 * x272 * x97 + x233 * x272 * x60 * x91 + x273 * x91,
            x233
            * (
                x13 * (2 * a * x275 - 4 * x206 - 4 * x207)
                + x248 * x52
                + x3 * (3 * x210 + 3 * x213 + 3 * x214 + x242 * x58 + 2 * x243 + 2 * x244)
            )
            + x276 * x60
            + x276 * x9,
        ]
    )
    return S


def kinetic3d_43(a, A, b, B):
    """Cartesian 3D (gf) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = -x5 - A[1]
    x7 = x6**2
    x8 = 2 * a**2
    x9 = -x0 - x8 * (x3 + x7)
    x10 = -x4 * (a * A[0] + b * B[0])
    x11 = -x10 - A[0]
    x12 = -x10 - B[0]
    x13 = a * x4
    x14 = b * x13
    x15 = numpy.exp(-x14 * (A[0] - B[0]) ** 2)
    x16 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x17 = x15 * x16
    x18 = x17 * x3
    x19 = x12 * x18
    x20 = x12**2 * x17
    x21 = x18 + x20
    x22 = x12 * x21
    x23 = x11 * x21
    x24 = 3 * x23
    x25 = x3 * (8 * x19 + x22 + x24)
    x26 = 3 * x18
    x27 = x3 * (3 * x20 + x26)
    x28 = 2 * x18
    x29 = x12 * x28
    x30 = x22 + x29
    x31 = x11 * x30
    x32 = x27 + x31
    x33 = x11 * x32
    x34 = x25 + x33
    x35 = 2 * x17
    x36 = x12 * x35
    x37 = x11 * x36 + x26
    x38 = x3 * (x20 + x37)
    x39 = x23 + x29
    x40 = x11 * x39
    x41 = 3 * x38 + 3 * x40
    x42 = x11 * x34 + x3 * (2 * x27 + 2 * x31 + x41)
    x43 = x38 + x40
    x44 = x11 * x43
    x45 = 4 * x18
    x46 = x12 * x45
    x47 = x11 * x17
    x48 = x12 * x17
    x49 = x3 * (x47 + x48)
    x50 = x12 * x47
    x51 = x18 + x50
    x52 = x11 * x51
    x53 = x3 * (2 * x23 + x46 + 2 * x49 + 2 * x52)
    x54 = x11 * x42 + x3 * (3 * x25 + 3 * x33 + 3 * x44 + 3 * x53)
    x55 = numpy.exp(-x14 * (A[1] - B[1]) ** 2)
    x56 = numpy.exp(-x14 * (A[2] - B[2]) ** 2)
    x57 = numpy.pi * x4 * x56
    x58 = x55 * x57
    x59 = x54 * x58
    x60 = -x4 * (a * A[2] + b * B[2])
    x61 = -x60 - A[2]
    x62 = x61**2
    x63 = -x0 - x8 * (x3 + x62)
    x64 = b * x4
    x65 = x64 * (2 * a * x34 - 2 * x22 - x46)
    x66 = x44 + x53
    x67 = x64 * (2 * a * x66 - 6 * x19 - x24)
    x68 = x11**2
    x69 = -x0 - x8 * (x3 + x68)
    x70 = x18 * x69
    x71 = x1 * x64
    x72 = x48 * x69 + x48 * x71
    x73 = x12 * x72
    x74 = -x35
    x75 = x13 * (x2 * x21 + x74)
    x76 = x73 + x75
    x77 = x70 + x76
    x78 = x11 * x77
    x79 = 4 * x14
    x80 = x3 * (x36 * x69 + x48 * x79)
    x81 = x14 * x39
    x82 = x12 * x77 + x13 * (x2 * x30 - 3 * x48)
    x83 = x3 * (3 * x78 + 4 * x80 + 6 * x81 + x82)
    x84 = 3 * x70
    x85 = x3 * (3 * x73 + 3 * x75 + x84)
    x86 = x80 + x82
    x87 = x11 * x86
    x88 = x32 * x71 + x85 + x87
    x89 = x11 * x88
    x90 = x49 + x52
    x91 = x1 * x90 - x36
    x92 = x2 * x4
    x93 = x47 * x69 + x47 * x71
    x94 = x3 * (x72 + x93)
    x95 = x11 * x72
    x96 = x51 * x71 + x70 + x95
    x97 = x11 * x96
    x98 = x3 * (2 * x78 + 2 * x80 + 4 * x81 + x91 * x92 + 2 * x94 + 2 * x97)
    x99 = x64 * (2 * a * x43 - 2 * x20 - x28)
    x100 = x51 * x79 + x84 + 2 * x95
    x101 = x3 * (x100 + x76)
    x102 = x39 * x71 + x78 + x80
    x103 = x102 * x11
    x104 = x101 + x103 + x99
    x105 = x104 * x11
    x106 = 3 * x101 + 3 * x103 + 3 * x99
    x107 = x65 + x83 + x89
    x108 = (
        x107 * x11
        + x3 * (x106 + x32 * x79 + 2 * x85 + 2 * x87)
        + x64 * (2 * a * x42 - 3 * x27 - 3 * x31)
    )
    x109 = -x5 - B[1]
    x110 = x17 * x68
    x111 = x3 * (x110 + x37)
    x112 = x11 * x90
    x113 = x11 * x66 + x3 * (2 * x111 + 2 * x112 + x41)
    x114 = x113 * x58
    x115 = x16 * x55
    x116 = x109 * x115
    x117 = x116 * x71 + x116 * x9
    x118 = x16 * x56
    x119 = x11 * x93
    x120 = x110 + x18
    x121 = x64 * (x1 * x120 + x74)
    x122 = x119 + x121
    x123 = x3 * (x100 + x122)
    x124 = x64 * x91
    x125 = x124 + x94 + x97
    x126 = x11 * x125
    x127 = x111 + x112
    x128 = 2 * a * x127 - x26 - 3 * x50
    x129 = x105 + x67 + x98
    x130 = x58 * (
        x11 * x129
        + x3 * (x106 + 2 * x123 + 2 * x126 + x128 * x92)
        + x64 * (2 * a * x113 - 4 * x38 - 4 * x40)
    )
    x131 = -x60 - B[2]
    x132 = x118 * x131
    x133 = x132 * x63 + x132 * x71
    x134 = x115 * x3
    x135 = x109**2 * x115
    x136 = x134 + x135
    x137 = x11 * x120 + x11 * x28
    x138 = x11 * x127 + x3 * (x137 + 3 * x49 + 3 * x52)
    x139 = x118 * x138
    x140 = x134 * x9
    x141 = x109 * x117
    x142 = 2 * x115
    x143 = -x142
    x144 = x13 * (x136 * x2 + x143)
    x145 = x141 + x144
    x146 = x140 + x145
    x147 = x122 + x70
    x148 = x11 * x147 + x3 * (x11 * x35 * x69 + x47 * x79) + x64 * (x1 * x137 - 3 * x47)
    x149 = x123 + x126 + x128 * x64
    x150 = (
        x11 * x149
        + x3 * (3 * x124 + x148 + 3 * x94 + 3 * x97)
        + x64 * (2 * a * x138 - 4 * x49 - 4 * x52)
    )
    x151 = x131 * x58
    x152 = x118 * x3
    x153 = x118 * x131**2
    x154 = x152 + x153
    x155 = x115 * x138
    x156 = x152 * x63
    x157 = x131 * x133
    x158 = 2 * x118
    x159 = -x158
    x160 = x13 * (x154 * x2 + x159)
    x161 = x157 + x160
    x162 = x156 + x161
    x163 = 2 * x134
    x164 = x109 * x163
    x165 = x109 * x136
    x166 = x164 + x165
    x167 = x11 * x137 + x3 * (3 * x110 + x26)
    x168 = x118 * x167
    x169 = x109 * x142
    x170 = x3 * (x116 * x79 + x169 * x9)
    x171 = x109 * x146 + x13 * (-3 * x116 + x166 * x2)
    x172 = x170 + x171
    x173 = (
        x11 * x148
        + x3 * (3 * x119 + 3 * x121 + x84)
        + x64 * (2 * a * x167 - 4 * x110 - x45)
    )
    x174 = 2 * x152
    x175 = x131 * x174
    x176 = x131 * x154
    x177 = x175 + x176
    x178 = x115 * x167
    x179 = x131 * x158
    x180 = x3 * (x132 * x79 + x179 * x63)
    x181 = x13 * (-3 * x132 + x177 * x2) + x131 * x162
    x182 = x180 + x181
    x183 = x42 * x58
    x184 = x115 * x6
    x185 = x184 * x71 + x184 * x9
    x186 = x108 * x58
    x187 = x109 * x184
    x188 = x134 + x187
    x189 = x118 * x66
    x190 = x117 * x6
    x191 = x140 + x188 * x71 + x190
    x192 = x136 * x6
    x193 = x164 + x192
    x194 = x118 * x127
    x195 = x146 * x6
    x196 = x170 + x193 * x71 + x195
    x197 = 3 * x134
    x198 = x3 * (3 * x135 + x197)
    x199 = x166 * x6
    x200 = x198 + x199
    x201 = x118 * x200
    x202 = 3 * x140
    x203 = x3 * (3 * x141 + 3 * x144 + x202)
    x204 = x172 * x6
    x205 = x200 * x71 + x203 + x204
    x206 = x118 * x61
    x207 = x206 * x63 + x206 * x71
    x208 = x58 * x61
    x209 = x131 * x206
    x210 = x152 + x209
    x211 = x115 * x66
    x212 = x133 * x61
    x213 = x156 + x210 * x71 + x212
    x214 = x154 * x61
    x215 = x175 + x214
    x216 = x115 * x127
    x217 = x162 * x61
    x218 = x180 + x215 * x71 + x217
    x219 = 3 * x152
    x220 = x3 * (3 * x153 + x219)
    x221 = x177 * x61
    x222 = x220 + x221
    x223 = x115 * x222
    x224 = 3 * x156
    x225 = x3 * (3 * x157 + 3 * x160 + x224)
    x226 = x182 * x61
    x227 = x222 * x71 + x225 + x226
    x228 = x115 * x7
    x229 = x134 + x228
    x230 = x118 * x34
    x231 = x185 * x6
    x232 = x64 * (x1 * x229 + x143)
    x233 = x231 + x232
    x234 = x140 + x233
    x235 = x3 * (x116 + x184)
    x236 = x188 * x6
    x237 = x235 + x236
    x238 = x118 * x237
    x239 = x1 * x237 - x169
    x240 = x239 * x64
    x241 = x3 * (x117 + x185)
    x242 = x191 * x6
    x243 = x240 + x241 + x242
    x244 = x169 * x6 + x197
    x245 = x3 * (x135 + x244)
    x246 = x193 * x6
    x247 = x245 + x246
    x248 = x118 * x90
    x249 = x64 * (2 * a * x247 - 2 * x135 - x163)
    x250 = x188 * x79 + 2 * x190 + x202
    x251 = x3 * (x145 + x250)
    x252 = x196 * x6
    x253 = x249 + x251 + x252
    x254 = x109 * x134
    x255 = 3 * x192
    x256 = x3 * (x165 + 8 * x254 + x255)
    x257 = x200 * x6
    x258 = x256 + x257
    x259 = x118 * x258
    x260 = 4 * x134
    x261 = x109 * x260
    x262 = x64 * (2 * a * x258 - 2 * x165 - x261)
    x263 = 6 * x14
    x264 = x3 * (4 * x170 + x171 + x193 * x263 + 3 * x195)
    x265 = x205 * x6
    x266 = x262 + x264 + x265
    x267 = x118 * x62
    x268 = x152 + x267
    x269 = x115 * x34
    x270 = x207 * x61
    x271 = x64 * (x1 * x268 + x159)
    x272 = x270 + x271
    x273 = x156 + x272
    x274 = x3 * (x132 + x206)
    x275 = x210 * x61
    x276 = x274 + x275
    x277 = x115 * x276
    x278 = x1 * x276 - x179
    x279 = x278 * x64
    x280 = x3 * (x133 + x207)
    x281 = x213 * x61
    x282 = x279 + x280 + x281
    x283 = x179 * x61 + x219
    x284 = x3 * (x153 + x283)
    x285 = x215 * x61
    x286 = x284 + x285
    x287 = x115 * x90
    x288 = x64 * (2 * a * x286 - 2 * x153 - x174)
    x289 = x210 * x79 + 2 * x212 + x224
    x290 = x3 * (x161 + x289)
    x291 = x218 * x61
    x292 = x288 + x290 + x291
    x293 = x131 * x152
    x294 = 3 * x214
    x295 = x3 * (x176 + 8 * x293 + x294)
    x296 = x222 * x61
    x297 = x295 + x296
    x298 = x115 * x297
    x299 = 4 * x152
    x300 = x131 * x299
    x301 = x64 * (2 * a * x297 - 2 * x176 - x300)
    x302 = x3 * (4 * x180 + x181 + x215 * x263 + 3 * x217)
    x303 = x227 * x61
    x304 = x301 + x302 + x303
    x305 = x163 * x6 + x229 * x6
    x306 = x118 * x32
    x307 = x234 * x6 + x3 * (x142 * x6 * x9 + x184 * x79) + x64 * (x1 * x305 - 3 * x184)
    x308 = x3 * (x228 + x244)
    x309 = x237 * x6
    x310 = x308 + x309
    x311 = x118 * x310
    x312 = 2 * a * x310 - 3 * x187 - x197
    x313 = x3 * (x233 + x250)
    x314 = x243 * x6
    x315 = x312 * x64 + x313 + x314
    x316 = x247 * x6
    x317 = x3 * (2 * x192 + 2 * x235 + 2 * x236 + x261)
    x318 = x316 + x317
    x319 = x118 * x318
    x320 = x64 * (2 * a * x318 - 6 * x254 - x255)
    x321 = x3 * (2 * x170 + x193 * x79 + 2 * x195 + x239 * x92 + 2 * x241 + 2 * x242)
    x322 = x253 * x6
    x323 = x320 + x321 + x322
    x324 = 3 * x245 + 3 * x246
    x325 = x258 * x6 + x3 * (2 * x198 + 2 * x199 + x324)
    x326 = x15 * x57
    x327 = x11 * x326
    x328 = 3 * x249 + 3 * x251 + 3 * x252
    x329 = (
        x266 * x6
        + x3 * (x200 * x79 + 2 * x203 + 2 * x204 + x328)
        + x64 * (2 * a * x325 - 3 * x198 - 3 * x199)
    )
    x330 = x326 * x329
    x331 = numpy.pi * x15 * x4 * x55
    x332 = x11 * x331
    x333 = x174 * x61 + x268 * x61
    x334 = x115 * x32
    x335 = (
        x273 * x61 + x3 * (x158 * x61 * x63 + x206 * x79) + x64 * (x1 * x333 - 3 * x206)
    )
    x336 = x3 * (x267 + x283)
    x337 = x276 * x61
    x338 = x336 + x337
    x339 = x115 * x338
    x340 = 2 * a * x338 - 3 * x209 - x219
    x341 = x3 * (x272 + x289)
    x342 = x282 * x61
    x343 = x340 * x64 + x341 + x342
    x344 = x286 * x61
    x345 = x3 * (2 * x214 + 2 * x274 + 2 * x275 + x300)
    x346 = x344 + x345
    x347 = x115 * x346
    x348 = x64 * (2 * a * x346 - 6 * x293 - x294)
    x349 = x3 * (2 * x180 + x215 * x79 + 2 * x217 + x278 * x92 + 2 * x280 + 2 * x281)
    x350 = x292 * x61
    x351 = x348 + x349 + x350
    x352 = 3 * x284 + 3 * x285
    x353 = x297 * x61 + x3 * (2 * x220 + 2 * x221 + x352)
    x354 = 3 * x288 + 3 * x290 + 3 * x291
    x355 = (
        x3 * (x222 * x79 + 2 * x225 + 2 * x226 + x354)
        + x304 * x61
        + x64 * (2 * a * x353 - 3 * x220 - 3 * x221)
    )
    x356 = x331 * x355
    x357 = x3 * (x197 + 3 * x228) + x305 * x6
    x358 = x118 * x357
    x359 = (
        x3 * (x202 + 3 * x231 + 3 * x232)
        + x307 * x6
        + x64 * (2 * a * x357 - 4 * x228 - x260)
    )
    x360 = x3 * (3 * x235 + 3 * x236 + x305) + x310 * x6
    x361 = x118 * x360
    x362 = (
        x3 * (3 * x240 + 3 * x241 + 3 * x242 + x307)
        + x315 * x6
        + x64 * (2 * a * x360 - 4 * x235 - 4 * x236)
    )
    x363 = x3 * (2 * x308 + 2 * x309 + x324) + x318 * x6
    x364 = x12 * x326
    x365 = x326 * (
        x3 * (x312 * x92 + 2 * x313 + 2 * x314 + x328)
        + x323 * x6
        + x64 * (2 * a * x363 - 4 * x245 - 4 * x246)
    )
    x366 = x3 * (3 * x256 + 3 * x257 + 3 * x316 + 3 * x317) + x325 * x6
    x367 = x326 * x366
    x368 = x326 * x69
    x369 = x17 * x360
    x370 = x17 * x357
    x371 = x17 * x318
    x372 = x17 * x310
    x373 = x17 * x222
    x374 = x17 * x258
    x375 = x17 * x276
    x376 = x17 * x237
    x377 = x17 * x297
    x378 = x331 * x6
    x379 = x17 * x200
    x380 = x17 * x338
    x381 = x17 * x346
    x382 = x3 * (x219 + 3 * x267) + x333 * x61
    x383 = x115 * x382
    x384 = (
        x3 * (x224 + 3 * x270 + 3 * x271)
        + x335 * x61
        + x64 * (2 * a * x382 - 4 * x267 - x299)
    )
    x385 = x3 * (3 * x274 + 3 * x275 + x333) + x338 * x61
    x386 = x115 * x385
    x387 = (
        x3 * (3 * x279 + 3 * x280 + 3 * x281 + x335)
        + x343 * x61
        + x64 * (2 * a * x385 - 4 * x274 - 4 * x275)
    )
    x388 = x12 * x331
    x389 = x3 * (2 * x336 + 2 * x337 + x352) + x346 * x61
    x390 = x331 * (
        x3 * (x340 * x92 + 2 * x341 + 2 * x342 + x354)
        + x351 * x61
        + x64 * (2 * a * x389 - 4 * x284 - 4 * x285)
    )
    x391 = x17 * x382
    x392 = x17 * x385
    x393 = x3 * (3 * x295 + 3 * x296 + 3 * x344 + 3 * x345) + x353 * x61
    x394 = x331 * x393

    # 150 item(s)
    S = numpy.array(
        [
            x58
            * (
                x108 * x11
                + x3 * (3 * x105 + 3 * x65 + 3 * x67 + 3 * x83 + 3 * x89 + 3 * x98)
                + x64 * (2 * a * x54 - 4 * x25 - 4 * x33)
            )
            + x59 * x63
            + x59 * x9,
            x109 * x114 * x63 + x109 * x130 + x113 * x117 * x118,
            x113 * x115 * x133 + x114 * x131 * x9 + x130 * x131,
            x118 * x136 * x150 + x136 * x139 * x63 + x139 * x146,
            x109 * x150 * x151 + x116 * x133 * x138 + x117 * x132 * x138,
            x115 * x150 * x154 + x154 * x155 * x9 + x155 * x162,
            x118 * x166 * x173 + x166 * x168 * x63 + x168 * x172,
            x132 * x136 * x173 + x132 * x146 * x167 + x133 * x136 * x167,
            x116 * x154 * x173 + x116 * x162 * x167 + x117 * x154 * x167,
            x115 * x173 * x177 + x177 * x178 * x9 + x178 * x182,
            x118 * x185 * x42 + x183 * x6 * x63 + x186 * x6,
            x118 * x129 * x188 + x188 * x189 * x63 + x189 * x191,
            x129 * x151 * x6 + x132 * x185 * x66 + x133 * x184 * x66,
            x118 * x149 * x193 + x193 * x194 * x63 + x194 * x196,
            x127 * x132 * x191 + x127 * x133 * x188 + x132 * x149 * x188,
            x127 * x154 * x185 + x127 * x162 * x184 + x149 * x154 * x184,
            x118 * x137 * x205 + x137 * x201 * x63 + x148 * x201,
            x132 * x137 * x196 + x132 * x148 * x193 + x133 * x137 * x193,
            x137 * x154 * x191 + x137 * x162 * x188 + x148 * x154 * x188,
            x137 * x177 * x185 + x137 * x182 * x184 + x148 * x177 * x184,
            x115 * x207 * x42 + x183 * x61 * x9 + x186 * x61,
            x109 * x129 * x208 + x116 * x207 * x66 + x117 * x206 * x66,
            x115 * x129 * x210 + x210 * x211 * x9 + x211 * x213,
            x127 * x136 * x207 + x127 * x146 * x206 + x136 * x149 * x206,
            x116 * x127 * x213 + x116 * x149 * x210 + x117 * x127 * x210,
            x115 * x149 * x215 + x215 * x216 * x9 + x216 * x218,
            x137 * x166 * x207 + x137 * x172 * x206 + x148 * x166 * x206,
            x136 * x137 * x213 + x136 * x148 * x210 + x137 * x146 * x210,
            x116 * x137 * x218 + x116 * x148 * x215 + x117 * x137 * x215,
            x115 * x137 * x227 + x137 * x223 * x9 + x148 * x223,
            x107 * x118 * x229 + x229 * x230 * x63 + x230 * x234,
            x104 * x238 + x118 * x243 * x43 + x238 * x43 * x63,
            x104 * x132 * x229 + x132 * x234 * x43 + x133 * x229 * x43,
            x118 * x125 * x247 + x247 * x248 * x63 + x248 * x253,
            x125 * x132 * x237 + x132 * x243 * x90 + x133 * x237 * x90,
            x125 * x154 * x229 + x154 * x234 * x90 + x162 * x229 * x90,
            x118 * x120 * x266 + x120 * x259 * x63 + x147 * x259,
            x120 * x132 * x253 + x120 * x133 * x247 + x132 * x147 * x247,
            x120 * x154 * x243 + x120 * x162 * x237 + x147 * x154 * x237,
            x120 * x177 * x234 + x120 * x182 * x229 + x147 * x177 * x229,
            x107 * x208 * x6 + x184 * x207 * x34 + x185 * x206 * x34,
            x104 * x188 * x206 + x188 * x207 * x43 + x191 * x206 * x43,
            x104 * x184 * x210 + x184 * x213 * x43 + x185 * x210 * x43,
            x125 * x193 * x206 + x193 * x207 * x90 + x196 * x206 * x90,
            x125 * x188 * x210 + x188 * x213 * x90 + x191 * x210 * x90,
            x125 * x184 * x215 + x184 * x218 * x90 + x185 * x215 * x90,
            x120 * x200 * x207 + x120 * x205 * x206 + x147 * x200 * x206,
            x120 * x193 * x213 + x120 * x196 * x210 + x147 * x193 * x210,
            x120 * x188 * x218 + x120 * x191 * x215 + x147 * x188 * x215,
            x120 * x184 * x227 + x120 * x185 * x222 + x147 * x184 * x222,
            x107 * x115 * x268 + x268 * x269 * x9 + x269 * x273,
            x104 * x116 * x268 + x116 * x273 * x43 + x117 * x268 * x43,
            x104 * x277 + x115 * x282 * x43 + x277 * x43 * x9,
            x125 * x136 * x268 + x136 * x273 * x90 + x146 * x268 * x90,
            x116 * x125 * x276 + x116 * x282 * x90 + x117 * x276 * x90,
            x115 * x125 * x286 + x286 * x287 * x9 + x287 * x292,
            x120 * x166 * x273 + x120 * x172 * x268 + x147 * x166 * x268,
            x120 * x136 * x282 + x120 * x146 * x276 + x136 * x147 * x276,
            x116 * x120 * x292 + x116 * x147 * x286 + x117 * x120 * x286,
            x115 * x120 * x304 + x120 * x298 * x9 + x147 * x298,
            x118 * x305 * x88 + x305 * x306 * x63 + x306 * x307,
            x102 * x311 + x118 * x315 * x39 + x311 * x39 * x63,
            x102 * x132 * x305 + x132 * x307 * x39 + x133 * x305 * x39,
            x118 * x323 * x51 + x319 * x51 * x63 + x319 * x96,
            x132 * x310 * x96 + x132 * x315 * x51 + x133 * x310 * x51,
            x154 * x305 * x96 + x154 * x307 * x51 + x162 * x305 * x51,
            x11 * x330 + x118 * x325 * x93 + x325 * x327 * x63,
            x131 * x323 * x327 + x132 * x318 * x93 + x133 * x318 * x47,
            x154 * x310 * x93 + x154 * x315 * x47 + x162 * x310 * x47,
            x177 * x305 * x93 + x177 * x307 * x47 + x182 * x305 * x47,
            x206 * x229 * x88 + x206 * x234 * x32 + x207 * x229 * x32,
            x102 * x206 * x237 + x206 * x243 * x39 + x207 * x237 * x39,
            x102 * x210 * x229 + x210 * x234 * x39 + x213 * x229 * x39,
            x206 * x247 * x96 + x206 * x253 * x51 + x207 * x247 * x51,
            x210 * x237 * x96 + x210 * x243 * x51 + x213 * x237 * x51,
            x215 * x229 * x96 + x215 * x234 * x51 + x218 * x229 * x51,
            x206 * x258 * x93 + x207 * x258 * x47 + x266 * x327 * x61,
            x210 * x247 * x93 + x210 * x253 * x47 + x213 * x247 * x47,
            x215 * x237 * x93 + x215 * x243 * x47 + x218 * x237 * x47,
            x222 * x229 * x93 + x222 * x234 * x47 + x227 * x229 * x47,
            x184 * x268 * x88 + x184 * x273 * x32 + x185 * x268 * x32,
            x102 * x188 * x268 + x188 * x273 * x39 + x191 * x268 * x39,
            x102 * x184 * x276 + x184 * x282 * x39 + x185 * x276 * x39,
            x193 * x268 * x96 + x193 * x273 * x51 + x196 * x268 * x51,
            x188 * x276 * x96 + x188 * x282 * x51 + x191 * x276 * x51,
            x184 * x286 * x96 + x184 * x292 * x51 + x185 * x286 * x51,
            x200 * x268 * x93 + x200 * x273 * x47 + x205 * x268 * x47,
            x193 * x276 * x93 + x193 * x282 * x47 + x196 * x276 * x47,
            x188 * x286 * x93 + x188 * x292 * x47 + x191 * x286 * x47,
            x184 * x297 * x93 + x185 * x297 * x47 + x304 * x332 * x6,
            x115 * x333 * x88 + x333 * x334 * x9 + x334 * x335,
            x102 * x116 * x333 + x116 * x335 * x39 + x117 * x333 * x39,
            x102 * x339 + x115 * x343 * x39 + x339 * x39 * x9,
            x136 * x333 * x96 + x136 * x335 * x51 + x146 * x333 * x51,
            x116 * x338 * x96 + x116 * x343 * x51 + x117 * x338 * x51,
            x115 * x351 * x51 + x347 * x51 * x9 + x347 * x96,
            x166 * x333 * x93 + x166 * x335 * x47 + x172 * x333 * x47,
            x136 * x338 * x93 + x136 * x343 * x47 + x146 * x338 * x47,
            x109 * x332 * x351 + x116 * x346 * x93 + x117 * x346 * x47,
            x11 * x356 + x115 * x353 * x93 + x332 * x353 * x9,
            x118 * x30 * x359 + x30 * x358 * x63 + x358 * x86,
            x118 * x21 * x362 + x21 * x361 * x63 + x361 * x77,
            x132 * x21 * x359 + x132 * x357 * x77 + x133 * x21 * x357,
            x118 * x363 * x72 + x12 * x365 + x363 * x364 * x63,
            x131 * x362 * x364 + x132 * x360 * x72 + x133 * x360 * x48,
            x154 * x357 * x72 + x154 * x359 * x48 + x162 * x357 * x48,
            x326
            * (
                x3 * (3 * x262 + 3 * x264 + 3 * x265 + 3 * x320 + 3 * x321 + 3 * x322)
                + x329 * x6
                + x64 * (2 * a * x366 - 4 * x256 - 4 * x257)
            )
            + x367 * x63
            + x367 * x69,
            x131 * x363 * x368 + x131 * x365 + x133 * x17 * x363,
            x154 * x17 * x362 + x154 * x369 * x69 + x162 * x369,
            x17 * x177 * x359 + x177 * x370 * x69 + x182 * x370,
            x206 * x30 * x307 + x206 * x305 * x86 + x207 * x30 * x305,
            x206 * x21 * x315 + x206 * x310 * x77 + x207 * x21 * x310,
            x21 * x210 * x307 + x21 * x213 * x305 + x210 * x305 * x77,
            x206 * x318 * x72 + x207 * x318 * x48 + x323 * x364 * x61,
            x210 * x310 * x72 + x210 * x315 * x48 + x213 * x310 * x48,
            x215 * x305 * x72 + x215 * x307 * x48 + x218 * x305 * x48,
            x17 * x207 * x325 + x325 * x368 * x61 + x330 * x61,
            x17 * x210 * x323 + x210 * x371 * x69 + x213 * x371,
            x17 * x215 * x315 + x215 * x372 * x69 + x218 * x372,
            x17 * x227 * x305 + x305 * x373 * x69 + x307 * x373,
            x229 * x268 * x86 + x229 * x273 * x30 + x234 * x268 * x30,
            x21 * x237 * x273 + x21 * x243 * x268 + x237 * x268 * x77,
            x21 * x229 * x282 + x21 * x234 * x276 + x229 * x276 * x77,
            x247 * x268 * x72 + x247 * x273 * x48 + x253 * x268 * x48,
            x237 * x276 * x72 + x237 * x282 * x48 + x243 * x276 * x48,
            x229 * x286 * x72 + x229 * x292 * x48 + x234 * x286 * x48,
            x17 * x266 * x268 + x268 * x374 * x69 + x273 * x374,
            x17 * x247 * x282 + x247 * x375 * x69 + x253 * x375,
            x17 * x243 * x286 + x286 * x376 * x69 + x292 * x376,
            x17 * x229 * x304 + x229 * x377 * x69 + x234 * x377,
            x184 * x30 * x335 + x184 * x333 * x86 + x185 * x30 * x333,
            x188 * x21 * x335 + x188 * x333 * x77 + x191 * x21 * x333,
            x184 * x21 * x343 + x184 * x338 * x77 + x185 * x21 * x338,
            x193 * x333 * x72 + x193 * x335 * x48 + x196 * x333 * x48,
            x188 * x338 * x72 + x188 * x343 * x48 + x191 * x338 * x48,
            x12 * x351 * x378 + x184 * x346 * x72 + x185 * x346 * x48,
            x17 * x205 * x333 + x333 * x379 * x69 + x335 * x379,
            x17 * x193 * x343 + x193 * x380 * x69 + x196 * x380,
            x17 * x188 * x351 + x188 * x381 * x69 + x191 * x381,
            x17 * x185 * x353 + x353 * x378 * x69 + x356 * x6,
            x115 * x30 * x384 + x30 * x383 * x9 + x383 * x86,
            x116 * x21 * x384 + x116 * x382 * x77 + x117 * x21 * x382,
            x115 * x21 * x387 + x21 * x386 * x9 + x386 * x77,
            x136 * x382 * x72 + x136 * x384 * x48 + x146 * x382 * x48,
            x109 * x387 * x388 + x116 * x385 * x72 + x117 * x385 * x48,
            x115 * x389 * x72 + x12 * x390 + x388 * x389 * x9,
            x166 * x17 * x384 + x166 * x391 * x69 + x172 * x391,
            x136 * x17 * x387 + x136 * x392 * x69 + x146 * x392,
            x109 * x331 * x389 * x69 + x109 * x390 + x117 * x17 * x389,
            x331
            * (
                x3 * (3 * x301 + 3 * x302 + 3 * x303 + 3 * x348 + 3 * x349 + 3 * x350)
                + x355 * x61
                + x64 * (2 * a * x393 - 4 * x295 - 4 * x296)
            )
            + x394 * x69
            + x394 * x9,
        ]
    )
    return S


def kinetic3d_44(a, A, b, B):
    """Cartesian 3D (gg) kinetic energy integral.

    Generated code; DO NOT modify by hand!"""

    x0 = -a
    x1 = 2 * a
    x2 = 2 * b
    x3 = (x1 + x2) ** (-1.0)
    x4 = (a + b) ** (-1.0)
    x5 = -x4 * (a * A[1] + b * B[1])
    x6 = -x5 - A[1]
    x7 = x6**2
    x8 = 2 * a**2
    x9 = -x0 - x8 * (x3 + x7)
    x10 = -x4 * (a * A[0] + b * B[0])
    x11 = -x10 - A[0]
    x12 = a * x4
    x13 = b * x12
    x14 = numpy.exp(-x13 * (A[0] - B[0]) ** 2)
    x15 = numpy.sqrt(numpy.pi) * numpy.sqrt(x4)
    x16 = x14 * x15
    x17 = x16 * x3
    x18 = 3 * x17
    x19 = -x10 - B[0]
    x20 = x16 * x19**2
    x21 = x3 * (x18 + 3 * x20)
    x22 = 2 * x17
    x23 = x19 * x22
    x24 = x17 + x20
    x25 = x19 * x24
    x26 = x23 + x25
    x27 = x19 * x26
    x28 = x11 * x26
    x29 = x3 * (5 * x21 + x27 + 4 * x28)
    x30 = x17 * x19
    x31 = 8 * x30
    x32 = x3 * (4 * x25 + x31)
    x33 = x21 + x27
    x34 = x11 * x33
    x35 = x32 + x34
    x36 = x11 * x35
    x37 = x29 + x36
    x38 = x11 * x24
    x39 = 3 * x38
    x40 = x3 * (x25 + x31 + x39)
    x41 = x21 + x28
    x42 = x11 * x41
    x43 = 4 * x40 + 4 * x42
    x44 = x11 * x37 + x3 * (2 * x32 + 2 * x34 + x43)
    x45 = x40 + x42
    x46 = x11 * x45
    x47 = 2 * x21
    x48 = 2 * x16
    x49 = x19 * x48
    x50 = x11 * x49 + x18
    x51 = x3 * (x20 + x50)
    x52 = x23 + x38
    x53 = x11 * x52
    x54 = 3 * x51 + 3 * x53
    x55 = x3 * (2 * x28 + x47 + x54)
    x56 = x11 * x44 + x3 * (3 * x29 + 3 * x36 + 4 * x46 + 4 * x55)
    x57 = numpy.exp(-x13 * (A[1] - B[1]) ** 2)
    x58 = numpy.exp(-x13 * (A[2] - B[2]) ** 2)
    x59 = numpy.pi * x4 * x58
    x60 = x57 * x59
    x61 = x56 * x60
    x62 = -x4 * (a * A[2] + b * B[2])
    x63 = -x62 - A[2]
    x64 = x63**2
    x65 = -x0 - x8 * (x3 + x64)
    x66 = b * x4
    x67 = x66 * (2 * a * x37 - 2 * x27 - x47)
    x68 = x46 + x55
    x69 = x66 * (2 * a * x68 - 3 * x21 - 3 * x28)
    x70 = x16 * x19
    x71 = 4 * x13
    x72 = x11**2
    x73 = -x0 - x8 * (x3 + x72)
    x74 = x3 * (x49 * x73 + x70 * x71)
    x75 = x17 * x73
    x76 = x1 * x66
    x77 = x70 * x73 + x70 * x76
    x78 = x19 * x77
    x79 = -x48
    x80 = x12 * (x2 * x24 + x79)
    x81 = x78 + x80
    x82 = x75 + x81
    x83 = x19 * x82
    x84 = x12 * (x2 * x26 - 3 * x70)
    x85 = x83 + x84
    x86 = x74 + x85
    x87 = x11 * x86
    x88 = 3 * x75
    x89 = x3 * (3 * x78 + 3 * x80 + x88)
    x90 = x13 * x41
    x91 = 4 * x17
    x92 = x12 * (2 * b * x33 - 4 * x20 - x91) + x19 * x86
    x93 = x3 * (4 * x87 + 5 * x89 + 8 * x90 + x92)
    x94 = 4 * x74
    x95 = x3 * (4 * x83 + 4 * x84 + x94)
    x96 = x89 + x92
    x97 = x11 * x96
    x98 = x35 * x76 + x95 + x97
    x99 = x11 * x98
    x100 = x11 * x77
    x101 = x11 * x16
    x102 = x101 * x19
    x103 = x102 + x17
    x104 = 2 * x100 + x103 * x71 + x88
    x105 = x3 * (x104 + x81)
    x106 = x11 * x82
    x107 = x106 + x52 * x76 + x74
    x108 = x107 * x11
    x109 = x51 + x53
    x110 = x66 * (2 * a * x109 - 2 * x20 - x22)
    x111 = 3 * x105 + 3 * x108 + 3 * x110
    x112 = x3 * (x111 + 2 * x87 + 2 * x89 + 4 * x90)
    x113 = x19 * x91
    x114 = x66 * (2 * a * x45 - x113 - 2 * x25)
    x115 = 6 * x13
    x116 = x3 * (3 * x106 + x115 * x52 + x85 + x94)
    x117 = x41 * x76 + x87 + x89
    x118 = x11 * x117
    x119 = x114 + x116 + x118
    x120 = x11 * x119
    x121 = x67 + x93 + x99
    x122 = (
        x11 * x121
        + x3 * (4 * x114 + 4 * x116 + 4 * x118 + x35 * x71 + 2 * x95 + 2 * x97)
        + x66 * (2 * a * x44 - 3 * x32 - 3 * x34)
    )
    x123 = -x5 - B[1]
    x124 = x109 * x11
    x125 = x3 * (x101 + x70)
    x126 = x103 * x11
    x127 = x3 * (x113 + 2 * x125 + 2 * x126 + 2 * x38)
    x128 = x11 * x68 + x3 * (3 * x124 + 3 * x127 + 3 * x40 + 3 * x42)
    x129 = x128 * x60
    x130 = x15 * x57
    x131 = x123 * x130
    x132 = x131 * x76 + x131 * x9
    x133 = x15 * x58
    x134 = x124 + x127
    x135 = x66 * (2 * a * x134 - 6 * x30 - x39)
    x136 = x125 + x126
    x137 = x1 * x136 - x49
    x138 = x2 * x4
    x139 = x101 * x73 + x101 * x76
    x140 = x3 * (x139 + x77)
    x141 = x100 + x103 * x76 + x75
    x142 = x11 * x141
    x143 = x3 * (2 * x106 + x137 * x138 + 2 * x140 + 2 * x142 + x52 * x71 + 2 * x74)
    x144 = x105 + x108 + x110
    x145 = x11 * x144
    x146 = x112 + x120 + x69
    x147 = x60 * (
        x11 * x146
        + x3 * (3 * x114 + 3 * x116 + 3 * x118 + 3 * x135 + 3 * x143 + 3 * x145)
        + x66 * (2 * a * x128 - x43)
    )
    x148 = -x62 - B[2]
    x149 = x133 * x148
    x150 = x149 * x65 + x149 * x76
    x151 = x130 * x3
    x152 = x123**2 * x130
    x153 = x151 + x152
    x154 = x16 * x72
    x155 = x3 * (x154 + x50)
    x156 = x11 * x136
    x157 = x11 * x134 + x3 * (2 * x155 + 2 * x156 + x54)
    x158 = x133 * x157
    x159 = x151 * x9
    x160 = x123 * x132
    x161 = 2 * x130
    x162 = -x161
    x163 = x12 * (x153 * x2 + x162)
    x164 = x160 + x163
    x165 = x159 + x164
    x166 = x11 * x139
    x167 = x154 + x17
    x168 = x66 * (x1 * x167 + x79)
    x169 = x166 + x168
    x170 = x3 * (x104 + x169)
    x171 = x137 * x66
    x172 = x140 + x142 + x171
    x173 = x11 * x172
    x174 = x155 + x156
    x175 = 2 * a * x174 - 3 * x102 - x18
    x176 = x135 + x143 + x145
    x177 = (
        x11 * x176
        + x3 * (x111 + x138 * x175 + 2 * x170 + 2 * x173)
        + x66 * (2 * a * x157 - 4 * x51 - 4 * x53)
    )
    x178 = x148 * x60
    x179 = x133 * x3
    x180 = x133 * x148**2
    x181 = x179 + x180
    x182 = x130 * x157
    x183 = x179 * x65
    x184 = x148 * x150
    x185 = 2 * x133
    x186 = -x185
    x187 = x12 * (x181 * x2 + x186)
    x188 = x184 + x187
    x189 = x183 + x188
    x190 = 2 * x151
    x191 = x123 * x190
    x192 = x123 * x153
    x193 = x191 + x192
    x194 = x11 * x167 + x11 * x22
    x195 = x11 * x174 + x3 * (3 * x125 + 3 * x126 + x194)
    x196 = x133 * x195
    x197 = x123 * x161
    x198 = x3 * (x131 * x71 + x197 * x9)
    x199 = x123 * x165
    x200 = x12 * (-3 * x131 + x193 * x2)
    x201 = x199 + x200
    x202 = x198 + x201
    x203 = x169 + x75
    x204 = x11 * x203 + x3 * (x101 * x71 + x11 * x48 * x73) + x66 * (x1 * x194 - 3 * x101)
    x205 = x170 + x173 + x175 * x66
    x206 = (
        x11 * x205
        + x3 * (3 * x140 + 3 * x142 + 3 * x171 + x204)
        + x66 * (2 * a * x195 - 4 * x125 - 4 * x126)
    )
    x207 = 2 * x179
    x208 = x148 * x207
    x209 = x148 * x181
    x210 = x208 + x209
    x211 = x130 * x195
    x212 = x148 * x185
    x213 = x3 * (x149 * x71 + x212 * x65)
    x214 = x148 * x189
    x215 = x12 * (-3 * x149 + x2 * x210)
    x216 = x214 + x215
    x217 = x213 + x216
    x218 = 3 * x151
    x219 = x3 * (3 * x152 + x218)
    x220 = x123 * x193
    x221 = x219 + x220
    x222 = x11 * x194 + x3 * (3 * x154 + x18)
    x223 = x133 * x222
    x224 = 3 * x159
    x225 = x3 * (3 * x160 + 3 * x163 + x224)
    x226 = 4 * x151
    x227 = x12 * (2 * b * x221 - 4 * x152 - x226) + x123 * x202
    x228 = x225 + x227
    x229 = (
        x11 * x204
        + x3 * (3 * x166 + 3 * x168 + x88)
        + x66 * (2 * a * x222 - 4 * x154 - x91)
    )
    x230 = 3 * x179
    x231 = x3 * (3 * x180 + x230)
    x232 = x148 * x210
    x233 = x231 + x232
    x234 = x130 * x222
    x235 = 3 * x183
    x236 = x3 * (3 * x184 + 3 * x187 + x235)
    x237 = 4 * x179
    x238 = x12 * (2 * b * x233 - 4 * x180 - x237) + x148 * x217
    x239 = x236 + x238
    x240 = x44 * x60
    x241 = x130 * x6
    x242 = x241 * x76 + x241 * x9
    x243 = x122 * x60
    x244 = x123 * x241
    x245 = x151 + x244
    x246 = x133 * x68
    x247 = x132 * x6
    x248 = x159 + x245 * x76 + x247
    x249 = x153 * x6
    x250 = x191 + x249
    x251 = x133 * x134
    x252 = x165 * x6
    x253 = x198 + x250 * x76 + x252
    x254 = x193 * x6
    x255 = x219 + x254
    x256 = x133 * x255
    x257 = x202 * x6
    x258 = x225 + x255 * x76 + x257
    x259 = x123 * x151
    x260 = 8 * x259
    x261 = x3 * (4 * x192 + x260)
    x262 = x221 * x6
    x263 = x261 + x262
    x264 = x133 * x263
    x265 = 4 * x198
    x266 = x3 * (4 * x199 + 4 * x200 + x265)
    x267 = x228 * x6
    x268 = x263 * x76 + x266 + x267
    x269 = x133 * x63
    x270 = x269 * x65 + x269 * x76
    x271 = x60 * x63
    x272 = x148 * x269
    x273 = x179 + x272
    x274 = x130 * x68
    x275 = x150 * x63
    x276 = x183 + x273 * x76 + x275
    x277 = x181 * x63
    x278 = x208 + x277
    x279 = x130 * x134
    x280 = x189 * x63
    x281 = x213 + x278 * x76 + x280
    x282 = x210 * x63
    x283 = x231 + x282
    x284 = x130 * x283
    x285 = x217 * x63
    x286 = x236 + x283 * x76 + x285
    x287 = x148 * x179
    x288 = 8 * x287
    x289 = x3 * (4 * x209 + x288)
    x290 = x233 * x63
    x291 = x289 + x290
    x292 = x130 * x291
    x293 = 4 * x213
    x294 = x3 * (4 * x214 + 4 * x215 + x293)
    x295 = x239 * x63
    x296 = x291 * x76 + x294 + x295
    x297 = x130 * x7
    x298 = x151 + x297
    x299 = x133 * x37
    x300 = x242 * x6
    x301 = x66 * (x1 * x298 + x162)
    x302 = x300 + x301
    x303 = x159 + x302
    x304 = x3 * (x131 + x241)
    x305 = x245 * x6
    x306 = x304 + x305
    x307 = x133 * x306
    x308 = x1 * x306 - x197
    x309 = x308 * x66
    x310 = x3 * (x132 + x242)
    x311 = x248 * x6
    x312 = x309 + x310 + x311
    x313 = x197 * x6 + x218
    x314 = x3 * (x152 + x313)
    x315 = x250 * x6
    x316 = x314 + x315
    x317 = x109 * x133
    x318 = x66 * (2 * a * x316 - 2 * x152 - x190)
    x319 = x224 + x245 * x71 + 2 * x247
    x320 = x3 * (x164 + x319)
    x321 = x253 * x6
    x322 = x318 + x320 + x321
    x323 = 3 * x249
    x324 = x3 * (x192 + x260 + x323)
    x325 = x255 * x6
    x326 = x324 + x325
    x327 = x133 * x136
    x328 = x123 * x226
    x329 = x66 * (2 * a * x326 - 2 * x192 - x328)
    x330 = x3 * (x115 * x250 + x201 + 3 * x252 + x265)
    x331 = x258 * x6
    x332 = x329 + x330 + x331
    x333 = x3 * (5 * x219 + x220 + 4 * x254)
    x334 = x263 * x6
    x335 = x333 + x334
    x336 = x133 * x335
    x337 = 2 * x219
    x338 = x66 * (2 * a * x335 - 2 * x220 - x337)
    x339 = 8 * x13
    x340 = x3 * (5 * x225 + x227 + x255 * x339 + 4 * x257)
    x341 = x268 * x6
    x342 = x338 + x340 + x341
    x343 = x133 * x64
    x344 = x179 + x343
    x345 = x130 * x37
    x346 = x270 * x63
    x347 = x66 * (x1 * x344 + x186)
    x348 = x346 + x347
    x349 = x183 + x348
    x350 = x3 * (x149 + x269)
    x351 = x273 * x63
    x352 = x350 + x351
    x353 = x130 * x352
    x354 = x1 * x352 - x212
    x355 = x354 * x66
    x356 = x3 * (x150 + x270)
    x357 = x276 * x63
    x358 = x355 + x356 + x357
    x359 = x212 * x63 + x230
    x360 = x3 * (x180 + x359)
    x361 = x278 * x63
    x362 = x360 + x361
    x363 = x109 * x130
    x364 = x66 * (2 * a * x362 - 2 * x180 - x207)
    x365 = x235 + x273 * x71 + 2 * x275
    x366 = x3 * (x188 + x365)
    x367 = x281 * x63
    x368 = x364 + x366 + x367
    x369 = 3 * x277
    x370 = x3 * (x209 + x288 + x369)
    x371 = x283 * x63
    x372 = x370 + x371
    x373 = x130 * x136
    x374 = x148 * x237
    x375 = x66 * (2 * a * x372 - 2 * x209 - x374)
    x376 = x3 * (x115 * x278 + x216 + 3 * x280 + x293)
    x377 = x286 * x63
    x378 = x375 + x376 + x377
    x379 = x3 * (5 * x231 + x232 + 4 * x282)
    x380 = x291 * x63
    x381 = x379 + x380
    x382 = x130 * x381
    x383 = 2 * x231
    x384 = x66 * (2 * a * x381 - 2 * x232 - x383)
    x385 = x3 * (5 * x236 + x238 + x283 * x339 + 4 * x285)
    x386 = x296 * x63
    x387 = x384 + x385 + x386
    x388 = x190 * x6 + x298 * x6
    x389 = x133 * x35
    x390 = x3 * (x161 * x6 * x9 + x241 * x71) + x303 * x6 + x66 * (x1 * x388 - 3 * x241)
    x391 = x3 * (x297 + x313)
    x392 = x306 * x6
    x393 = x391 + x392
    x394 = x133 * x41
    x395 = 2 * a * x393 - x218 - 3 * x244
    x396 = x3 * (x302 + x319)
    x397 = x312 * x6
    x398 = x395 * x66 + x396 + x397
    x399 = x316 * x6
    x400 = x3 * (2 * x249 + 2 * x304 + 2 * x305 + x328)
    x401 = x399 + x400
    x402 = x133 * x401
    x403 = x66 * (2 * a * x401 - 6 * x259 - x323)
    x404 = x3 * (x138 * x308 + 2 * x198 + x250 * x71 + 2 * x252 + 2 * x310 + 2 * x311)
    x405 = x322 * x6
    x406 = x403 + x404 + x405
    x407 = x326 * x6
    x408 = 3 * x314 + 3 * x315
    x409 = x3 * (2 * x254 + x337 + x408)
    x410 = x407 + x409
    x411 = x133 * x410
    x412 = x66 * (2 * a * x410 - 3 * x219 - 3 * x254)
    x413 = 3 * x318 + 3 * x320 + 3 * x321
    x414 = x3 * (2 * x225 + x255 * x71 + 2 * x257 + x413)
    x415 = x332 * x6
    x416 = x412 + x414 + x415
    x417 = 4 * x324 + 4 * x325
    x418 = x3 * (2 * x261 + 2 * x262 + x417) + x335 * x6
    x419 = x14 * x59
    x420 = x11 * x419
    x421 = (
        x3 * (x263 * x71 + 2 * x266 + 2 * x267 + 4 * x329 + 4 * x330 + 4 * x331)
        + x342 * x6
        + x66 * (2 * a * x418 - 3 * x261 - 3 * x262)
    )
    x422 = x419 * x421
    x423 = numpy.pi * x14 * x4 * x57
    x424 = x11 * x423
    x425 = x207 * x63 + x344 * x63
    x426 = x130 * x35
    x427 = (
        x3 * (x185 * x63 * x65 + x269 * x71) + x349 * x63 + x66 * (x1 * x425 - 3 * x269)
    )
    x428 = x3 * (x343 + x359)
    x429 = x352 * x63
    x430 = x428 + x429
    x431 = x130 * x41
    x432 = 2 * a * x430 - x230 - 3 * x272
    x433 = x3 * (x348 + x365)
    x434 = x358 * x63
    x435 = x432 * x66 + x433 + x434
    x436 = x362 * x63
    x437 = x3 * (2 * x277 + 2 * x350 + 2 * x351 + x374)
    x438 = x436 + x437
    x439 = x130 * x438
    x440 = x66 * (2 * a * x438 - 6 * x287 - x369)
    x441 = x3 * (x138 * x354 + 2 * x213 + x278 * x71 + 2 * x280 + 2 * x356 + 2 * x357)
    x442 = x368 * x63
    x443 = x440 + x441 + x442
    x444 = x372 * x63
    x445 = 3 * x360 + 3 * x361
    x446 = x3 * (2 * x282 + x383 + x445)
    x447 = x444 + x446
    x448 = x130 * x447
    x449 = x66 * (2 * a * x447 - 3 * x231 - 3 * x282)
    x450 = 3 * x364 + 3 * x366 + 3 * x367
    x451 = x3 * (2 * x236 + x283 * x71 + 2 * x285 + x450)
    x452 = x378 * x63
    x453 = x449 + x451 + x452
    x454 = 4 * x370 + 4 * x371
    x455 = x3 * (2 * x289 + 2 * x290 + x454) + x381 * x63
    x456 = (
        x3 * (x291 * x71 + 2 * x294 + 2 * x295 + 4 * x375 + 4 * x376 + 4 * x377)
        + x387 * x63
        + x66 * (2 * a * x455 - 3 * x289 - 3 * x290)
    )
    x457 = x423 * x456
    x458 = x3 * (x218 + 3 * x297) + x388 * x6
    x459 = x133 * x458
    x460 = (
        x3 * (x224 + 3 * x300 + 3 * x301)
        + x390 * x6
        + x66 * (2 * a * x458 - x226 - 4 * x297)
    )
    x461 = x3 * (3 * x304 + 3 * x305 + x388) + x393 * x6
    x462 = x133 * x461
    x463 = (
        x3 * (3 * x309 + 3 * x310 + 3 * x311 + x390)
        + x398 * x6
        + x66 * (2 * a * x461 - 4 * x304 - 4 * x305)
    )
    x464 = x3 * (2 * x391 + 2 * x392 + x408) + x401 * x6
    x465 = x133 * x464
    x466 = (
        x3 * (x138 * x395 + 2 * x396 + 2 * x397 + x413)
        + x406 * x6
        + x66 * (2 * a * x464 - 4 * x314 - 4 * x315)
    )
    x467 = x3 * (3 * x324 + 3 * x325 + 3 * x399 + 3 * x400) + x410 * x6
    x468 = x19 * x419
    x469 = x419 * (
        x3 * (3 * x329 + 3 * x330 + 3 * x331 + 3 * x403 + 3 * x404 + 3 * x405)
        + x416 * x6
        + x66 * (2 * a * x467 - x417)
    )
    x470 = x3 * (3 * x333 + 3 * x334 + 4 * x407 + 4 * x409) + x418 * x6
    x471 = x419 * x470
    x472 = x419 * x73
    x473 = x16 * x464
    x474 = x16 * x461
    x475 = x16 * x458
    x476 = x16 * x410
    x477 = x16 * x401
    x478 = x16 * x283
    x479 = x16 * x291
    x480 = x16 * x335
    x481 = x16 * x352
    x482 = x16 * x316
    x483 = x16 * x306
    x484 = x16 * x381
    x485 = x423 * x6
    x486 = x16 * x263
    x487 = x16 * x255
    x488 = x16 * x438
    x489 = x16 * x447
    x490 = x3 * (x230 + 3 * x343) + x425 * x63
    x491 = x130 * x490
    x492 = (
        x3 * (x235 + 3 * x346 + 3 * x347)
        + x427 * x63
        + x66 * (2 * a * x490 - x237 - 4 * x343)
    )
    x493 = x3 * (3 * x350 + 3 * x351 + x425) + x430 * x63
    x494 = x130 * x493
    x495 = (
        x3 * (3 * x355 + 3 * x356 + 3 * x357 + x427)
        + x435 * x63
        + x66 * (2 * a * x493 - 4 * x350 - 4 * x351)
    )
    x496 = x3 * (2 * x428 + 2 * x429 + x445) + x438 * x63
    x497 = x130 * x496
    x498 = (
        x3 * (x138 * x432 + 2 * x433 + 2 * x434 + x450)
        + x443 * x63
        + x66 * (2 * a * x496 - 4 * x360 - 4 * x361)
    )
    x499 = x19 * x423
    x500 = x3 * (3 * x370 + 3 * x371 + 3 * x436 + 3 * x437) + x447 * x63
    x501 = x423 * (
        x3 * (3 * x375 + 3 * x376 + 3 * x377 + 3 * x440 + 3 * x441 + 3 * x442)
        + x453 * x63
        + x66 * (2 * a * x500 - x454)
    )
    x502 = x16 * x490
    x503 = x16 * x493
    x504 = x16 * x496
    x505 = x3 * (3 * x379 + 3 * x380 + 4 * x444 + 4 * x446) + x455 * x63
    x506 = x423 * x505

    # 225 item(s)
    S = numpy.array(
        [
            x60
            * (
                x11 * x122
                + x3 * (4 * x112 + 4 * x120 + 3 * x67 + 4 * x69 + 3 * x93 + 3 * x99)
                + x66 * (2 * a * x56 - 4 * x29 - 4 * x36)
            )
            + x61 * x65
            + x61 * x9,
            x123 * x129 * x65 + x123 * x147 + x128 * x132 * x133,
            x128 * x130 * x150 + x129 * x148 * x9 + x147 * x148,
            x133 * x153 * x177 + x153 * x158 * x65 + x158 * x165,
            x123 * x177 * x178 + x131 * x150 * x157 + x132 * x149 * x157,
            x130 * x177 * x181 + x181 * x182 * x9 + x182 * x189,
            x133 * x193 * x206 + x193 * x196 * x65 + x196 * x202,
            x149 * x153 * x206 + x149 * x165 * x195 + x150 * x153 * x195,
            x131 * x181 * x206 + x131 * x189 * x195 + x132 * x181 * x195,
            x130 * x206 * x210 + x210 * x211 * x9 + x211 * x217,
            x133 * x221 * x229 + x221 * x223 * x65 + x223 * x228,
            x149 * x193 * x229 + x149 * x202 * x222 + x150 * x193 * x222,
            x153 * x181 * x229 + x153 * x189 * x222 + x165 * x181 * x222,
            x131 * x210 * x229 + x131 * x217 * x222 + x132 * x210 * x222,
            x130 * x229 * x233 + x233 * x234 * x9 + x234 * x239,
            x133 * x242 * x44 + x240 * x6 * x65 + x243 * x6,
            x133 * x146 * x245 + x245 * x246 * x65 + x246 * x248,
            x146 * x178 * x6 + x149 * x242 * x68 + x150 * x241 * x68,
            x133 * x176 * x250 + x250 * x251 * x65 + x251 * x253,
            x134 * x149 * x248 + x134 * x150 * x245 + x149 * x176 * x245,
            x134 * x181 * x242 + x134 * x189 * x241 + x176 * x181 * x241,
            x133 * x174 * x258 + x174 * x256 * x65 + x205 * x256,
            x149 * x174 * x253 + x149 * x205 * x250 + x150 * x174 * x250,
            x174 * x181 * x248 + x174 * x189 * x245 + x181 * x205 * x245,
            x174 * x210 * x242 + x174 * x217 * x241 + x205 * x210 * x241,
            x133 * x194 * x268 + x194 * x264 * x65 + x204 * x264,
            x149 * x194 * x258 + x149 * x204 * x255 + x150 * x194 * x255,
            x181 * x194 * x253 + x181 * x204 * x250 + x189 * x194 * x250,
            x194 * x210 * x248 + x194 * x217 * x245 + x204 * x210 * x245,
            x194 * x233 * x242 + x194 * x239 * x241 + x204 * x233 * x241,
            x130 * x270 * x44 + x240 * x63 * x9 + x243 * x63,
            x123 * x146 * x271 + x131 * x270 * x68 + x132 * x269 * x68,
            x130 * x146 * x273 + x273 * x274 * x9 + x274 * x276,
            x134 * x153 * x270 + x134 * x165 * x269 + x153 * x176 * x269,
            x131 * x134 * x276 + x131 * x176 * x273 + x132 * x134 * x273,
            x130 * x176 * x278 + x278 * x279 * x9 + x279 * x281,
            x174 * x193 * x270 + x174 * x202 * x269 + x193 * x205 * x269,
            x153 * x174 * x276 + x153 * x205 * x273 + x165 * x174 * x273,
            x131 * x174 * x281 + x131 * x205 * x278 + x132 * x174 * x278,
            x130 * x174 * x286 + x174 * x284 * x9 + x205 * x284,
            x194 * x221 * x270 + x194 * x228 * x269 + x204 * x221 * x269,
            x193 * x194 * x276 + x193 * x204 * x273 + x194 * x202 * x273,
            x153 * x194 * x281 + x153 * x204 * x278 + x165 * x194 * x278,
            x131 * x194 * x286 + x131 * x204 * x283 + x132 * x194 * x283,
            x130 * x194 * x296 + x194 * x292 * x9 + x204 * x292,
            x121 * x133 * x298 + x298 * x299 * x65 + x299 * x303,
            x119 * x307 + x133 * x312 * x45 + x307 * x45 * x65,
            x119 * x149 * x298 + x149 * x303 * x45 + x150 * x298 * x45,
            x133 * x144 * x316 + x316 * x317 * x65 + x317 * x322,
            x109 * x149 * x312 + x109 * x150 * x306 + x144 * x149 * x306,
            x109 * x181 * x303 + x109 * x189 * x298 + x144 * x181 * x298,
            x133 * x172 * x326 + x326 * x327 * x65 + x327 * x332,
            x136 * x149 * x322 + x136 * x150 * x316 + x149 * x172 * x316,
            x136 * x181 * x312 + x136 * x189 * x306 + x172 * x181 * x306,
            x136 * x210 * x303 + x136 * x217 * x298 + x172 * x210 * x298,
            x133 * x167 * x342 + x167 * x336 * x65 + x203 * x336,
            x149 * x167 * x332 + x149 * x203 * x326 + x150 * x167 * x326,
            x167 * x181 * x322 + x167 * x189 * x316 + x181 * x203 * x316,
            x167 * x210 * x312 + x167 * x217 * x306 + x203 * x210 * x306,
            x167 * x233 * x303 + x167 * x239 * x298 + x203 * x233 * x298,
            x121 * x271 * x6 + x241 * x270 * x37 + x242 * x269 * x37,
            x119 * x245 * x269 + x245 * x270 * x45 + x248 * x269 * x45,
            x119 * x241 * x273 + x241 * x276 * x45 + x242 * x273 * x45,
            x109 * x250 * x270 + x109 * x253 * x269 + x144 * x250 * x269,
            x109 * x245 * x276 + x109 * x248 * x273 + x144 * x245 * x273,
            x109 * x241 * x281 + x109 * x242 * x278 + x144 * x241 * x278,
            x136 * x255 * x270 + x136 * x258 * x269 + x172 * x255 * x269,
            x136 * x250 * x276 + x136 * x253 * x273 + x172 * x250 * x273,
            x136 * x245 * x281 + x136 * x248 * x278 + x172 * x245 * x278,
            x136 * x241 * x286 + x136 * x242 * x283 + x172 * x241 * x283,
            x167 * x263 * x270 + x167 * x268 * x269 + x203 * x263 * x269,
            x167 * x255 * x276 + x167 * x258 * x273 + x203 * x255 * x273,
            x167 * x250 * x281 + x167 * x253 * x278 + x203 * x250 * x278,
            x167 * x245 * x286 + x167 * x248 * x283 + x203 * x245 * x283,
            x167 * x241 * x296 + x167 * x242 * x291 + x203 * x241 * x291,
            x121 * x130 * x344 + x344 * x345 * x9 + x345 * x349,
            x119 * x131 * x344 + x131 * x349 * x45 + x132 * x344 * x45,
            x119 * x353 + x130 * x358 * x45 + x353 * x45 * x9,
            x109 * x153 * x349 + x109 * x165 * x344 + x144 * x153 * x344,
            x109 * x131 * x358 + x109 * x132 * x352 + x131 * x144 * x352,
            x130 * x144 * x362 + x362 * x363 * x9 + x363 * x368,
            x136 * x193 * x349 + x136 * x202 * x344 + x172 * x193 * x344,
            x136 * x153 * x358 + x136 * x165 * x352 + x153 * x172 * x352,
            x131 * x136 * x368 + x131 * x172 * x362 + x132 * x136 * x362,
            x130 * x172 * x372 + x372 * x373 * x9 + x373 * x378,
            x167 * x221 * x349 + x167 * x228 * x344 + x203 * x221 * x344,
            x167 * x193 * x358 + x167 * x202 * x352 + x193 * x203 * x352,
            x153 * x167 * x368 + x153 * x203 * x362 + x165 * x167 * x362,
            x131 * x167 * x378 + x131 * x203 * x372 + x132 * x167 * x372,
            x130 * x167 * x387 + x167 * x382 * x9 + x203 * x382,
            x133 * x388 * x98 + x388 * x389 * x65 + x389 * x390,
            x117 * x133 * x393 + x393 * x394 * x65 + x394 * x398,
            x117 * x149 * x388 + x149 * x390 * x41 + x150 * x388 * x41,
            x107 * x402 + x133 * x406 * x52 + x402 * x52 * x65,
            x107 * x149 * x393 + x149 * x398 * x52 + x150 * x393 * x52,
            x107 * x181 * x388 + x181 * x390 * x52 + x189 * x388 * x52,
            x103 * x133 * x416 + x103 * x411 * x65 + x141 * x411,
            x103 * x149 * x406 + x103 * x150 * x401 + x141 * x149 * x401,
            x103 * x181 * x398 + x103 * x189 * x393 + x141 * x181 * x393,
            x103 * x210 * x390 + x103 * x217 * x388 + x141 * x210 * x388,
            x11 * x422 + x133 * x139 * x418 + x418 * x420 * x65,
            x101 * x150 * x410 + x139 * x149 * x410 + x148 * x416 * x420,
            x101 * x181 * x406 + x101 * x189 * x401 + x139 * x181 * x401,
            x101 * x210 * x398 + x101 * x217 * x393 + x139 * x210 * x393,
            x101 * x233 * x390 + x101 * x239 * x388 + x139 * x233 * x388,
            x269 * x298 * x98 + x269 * x303 * x35 + x270 * x298 * x35,
            x117 * x269 * x306 + x269 * x312 * x41 + x270 * x306 * x41,
            x117 * x273 * x298 + x273 * x303 * x41 + x276 * x298 * x41,
            x107 * x269 * x316 + x269 * x322 * x52 + x270 * x316 * x52,
            x107 * x273 * x306 + x273 * x312 * x52 + x276 * x306 * x52,
            x107 * x278 * x298 + x278 * x303 * x52 + x281 * x298 * x52,
            x103 * x269 * x332 + x103 * x270 * x326 + x141 * x269 * x326,
            x103 * x273 * x322 + x103 * x276 * x316 + x141 * x273 * x316,
            x103 * x278 * x312 + x103 * x281 * x306 + x141 * x278 * x306,
            x103 * x283 * x303 + x103 * x286 * x298 + x141 * x283 * x298,
            x101 * x270 * x335 + x139 * x269 * x335 + x342 * x420 * x63,
            x101 * x273 * x332 + x101 * x276 * x326 + x139 * x273 * x326,
            x101 * x278 * x322 + x101 * x281 * x316 + x139 * x278 * x316,
            x101 * x283 * x312 + x101 * x286 * x306 + x139 * x283 * x306,
            x101 * x291 * x303 + x101 * x296 * x298 + x139 * x291 * x298,
            x241 * x344 * x98 + x241 * x349 * x35 + x242 * x344 * x35,
            x117 * x245 * x344 + x245 * x349 * x41 + x248 * x344 * x41,
            x117 * x241 * x352 + x241 * x358 * x41 + x242 * x352 * x41,
            x107 * x250 * x344 + x250 * x349 * x52 + x253 * x344 * x52,
            x107 * x245 * x352 + x245 * x358 * x52 + x248 * x352 * x52,
            x107 * x241 * x362 + x241 * x368 * x52 + x242 * x362 * x52,
            x103 * x255 * x349 + x103 * x258 * x344 + x141 * x255 * x344,
            x103 * x250 * x358 + x103 * x253 * x352 + x141 * x250 * x352,
            x103 * x245 * x368 + x103 * x248 * x362 + x141 * x245 * x362,
            x103 * x241 * x378 + x103 * x242 * x372 + x141 * x241 * x372,
            x101 * x263 * x349 + x101 * x268 * x344 + x139 * x263 * x344,
            x101 * x255 * x358 + x101 * x258 * x352 + x139 * x255 * x352,
            x101 * x250 * x368 + x101 * x253 * x362 + x139 * x250 * x362,
            x101 * x245 * x378 + x101 * x248 * x372 + x139 * x245 * x372,
            x101 * x242 * x381 + x139 * x241 * x381 + x387 * x424 * x6,
            x130 * x425 * x98 + x425 * x426 * x9 + x426 * x427,
            x117 * x131 * x425 + x131 * x41 * x427 + x132 * x41 * x425,
            x117 * x130 * x430 + x430 * x431 * x9 + x431 * x435,
            x107 * x153 * x425 + x153 * x427 * x52 + x165 * x425 * x52,
            x107 * x131 * x430 + x131 * x435 * x52 + x132 * x430 * x52,
            x107 * x439 + x130 * x443 * x52 + x439 * x52 * x9,
            x103 * x193 * x427 + x103 * x202 * x425 + x141 * x193 * x425,
            x103 * x153 * x435 + x103 * x165 * x430 + x141 * x153 * x430,
            x103 * x131 * x443 + x103 * x132 * x438 + x131 * x141 * x438,
            x103 * x130 * x453 + x103 * x448 * x9 + x141 * x448,
            x101 * x221 * x427 + x101 * x228 * x425 + x139 * x221 * x425,
            x101 * x193 * x435 + x101 * x202 * x430 + x139 * x193 * x430,
            x101 * x153 * x443 + x101 * x165 * x438 + x139 * x153 * x438,
            x101 * x132 * x447 + x123 * x424 * x453 + x131 * x139 * x447,
            x11 * x457 + x130 * x139 * x455 + x424 * x455 * x9,
            x133 * x33 * x460 + x33 * x459 * x65 + x459 * x96,
            x133 * x26 * x463 + x26 * x462 * x65 + x462 * x86,
            x149 * x26 * x460 + x149 * x458 * x86 + x150 * x26 * x458,
            x133 * x24 * x466 + x24 * x465 * x65 + x465 * x82,
            x149 * x24 * x463 + x149 * x461 * x82 + x150 * x24 * x461,
            x181 * x24 * x460 + x181 * x458 * x82 + x189 * x24 * x458,
            x133 * x467 * x77 + x19 * x469 + x467 * x468 * x65,
            x148 * x466 * x468 + x149 * x464 * x77 + x150 * x464 * x70,
            x181 * x461 * x77 + x181 * x463 * x70 + x189 * x461 * x70,
            x210 * x458 * x77 + x210 * x460 * x70 + x217 * x458 * x70,
            x419
            * (
                x3 * (3 * x338 + 3 * x340 + 3 * x341 + 4 * x412 + 4 * x414 + 4 * x415)
                + x421 * x6
                + x66 * (2 * a * x470 - 4 * x333 - 4 * x334)
            )
            + x471 * x65
            + x471 * x73,
            x148 * x467 * x472 + x148 * x469 + x150 * x16 * x467,
            x16 * x181 * x466 + x181 * x473 * x73 + x189 * x473,
            x16 * x210 * x463 + x210 * x474 * x73 + x217 * x474,
            x16 * x233 * x460 + x233 * x475 * x73 + x239 * x475,
            x269 * x33 * x390 + x269 * x388 * x96 + x270 * x33 * x388,
            x26 * x269 * x398 + x26 * x270 * x393 + x269 * x393 * x86,
            x26 * x273 * x390 + x26 * x276 * x388 + x273 * x388 * x86,
            x24 * x269 * x406 + x24 * x270 * x401 + x269 * x401 * x82,
            x24 * x273 * x398 + x24 * x276 * x393 + x273 * x393 * x82,
            x24 * x278 * x390 + x24 * x281 * x388 + x278 * x388 * x82,
            x269 * x410 * x77 + x270 * x410 * x70 + x416 * x468 * x63,
            x273 * x401 * x77 + x273 * x406 * x70 + x276 * x401 * x70,
            x278 * x393 * x77 + x278 * x398 * x70 + x281 * x393 * x70,
            x283 * x388 * x77 + x283 * x390 * x70 + x286 * x388 * x70,
            x16 * x270 * x418 + x418 * x472 * x63 + x422 * x63,
            x16 * x273 * x416 + x273 * x476 * x73 + x276 * x476,
            x16 * x278 * x406 + x278 * x477 * x73 + x281 * x477,
            x16 * x286 * x393 + x393 * x478 * x73 + x398 * x478,
            x16 * x296 * x388 + x388 * x479 * x73 + x390 * x479,
            x298 * x33 * x349 + x298 * x344 * x96 + x303 * x33 * x344,
            x26 * x306 * x349 + x26 * x312 * x344 + x306 * x344 * x86,
            x26 * x298 * x358 + x26 * x303 * x352 + x298 * x352 * x86,
            x24 * x316 * x349 + x24 * x322 * x344 + x316 * x344 * x82,
            x24 * x306 * x358 + x24 * x312 * x352 + x306 * x352 * x82,
            x24 * x298 * x368 + x24 * x303 * x362 + x298 * x362 * x82,
            x326 * x344 * x77 + x326 * x349 * x70 + x332 * x344 * x70,
            x316 * x352 * x77 + x316 * x358 * x70 + x322 * x352 * x70,
            x306 * x362 * x77 + x306 * x368 * x70 + x312 * x362 * x70,
            x298 * x372 * x77 + x298 * x378 * x70 + x303 * x372 * x70,
            x16 * x342 * x344 + x344 * x480 * x73 + x349 * x480,
            x16 * x326 * x358 + x326 * x481 * x73 + x332 * x481,
            x16 * x322 * x362 + x362 * x482 * x73 + x368 * x482,
            x16 * x312 * x372 + x372 * x483 * x73 + x378 * x483,
            x16 * x298 * x387 + x298 * x484 * x73 + x303 * x484,
            x241 * x33 * x427 + x241 * x425 * x96 + x242 * x33 * x425,
            x245 * x26 * x427 + x245 * x425 * x86 + x248 * x26 * x425,
            x241 * x26 * x435 + x241 * x430 * x86 + x242 * x26 * x430,
            x24 * x250 * x427 + x24 * x253 * x425 + x250 * x425 * x82,
            x24 * x245 * x435 + x24 * x248 * x430 + x245 * x430 * x82,
            x24 * x241 * x443 + x24 * x242 * x438 + x241 * x438 * x82,
            x255 * x425 * x77 + x255 * x427 * x70 + x258 * x425 * x70,
            x250 * x430 * x77 + x250 * x435 * x70 + x253 * x430 * x70,
            x245 * x438 * x77 + x245 * x443 * x70 + x248 * x438 * x70,
            x19 * x453 * x485 + x241 * x447 * x77 + x242 * x447 * x70,
            x16 * x268 * x425 + x425 * x486 * x73 + x427 * x486,
            x16 * x258 * x430 + x430 * x487 * x73 + x435 * x487,
            x16 * x250 * x443 + x250 * x488 * x73 + x253 * x488,
            x16 * x245 * x453 + x245 * x489 * x73 + x248 * x489,
            x16 * x242 * x455 + x455 * x485 * x73 + x457 * x6,
            x130 * x33 * x492 + x33 * x491 * x9 + x491 * x96,
            x131 * x26 * x492 + x131 * x490 * x86 + x132 * x26 * x490,
            x130 * x26 * x495 + x26 * x494 * x9 + x494 * x86,
            x153 * x24 * x492 + x153 * x490 * x82 + x165 * x24 * x490,
            x131 * x24 * x495 + x131 * x493 * x82 + x132 * x24 * x493,
            x130 * x24 * x498 + x24 * x497 * x9 + x497 * x82,
            x193 * x490 * x77 + x193 * x492 * x70 + x202 * x490 * x70,
            x153 * x493 * x77 + x153 * x495 * x70 + x165 * x493 * x70,
            x123 * x498 * x499 + x131 * x496 * x77 + x132 * x496 * x70,
            x130 * x500 * x77 + x19 * x501 + x499 * x500 * x9,
            x16 * x221 * x492 + x221 * x502 * x73 + x228 * x502,
            x16 * x193 * x495 + x193 * x503 * x73 + x202 * x503,
            x153 * x16 * x498 + x153 * x504 * x73 + x165 * x504,
            x123 * x423 * x500 * x73 + x123 * x501 + x132 * x16 * x500,
            x423
            * (
                x3 * (3 * x384 + 3 * x385 + 3 * x386 + 4 * x449 + 4 * x451 + 4 * x452)
                + x456 * x63
                + x66 * (2 * a * x505 - 4 * x379 - 4 * x380)
            )
            + x506 * x73
            + x506 * x9,
        ]
    )
    return S
