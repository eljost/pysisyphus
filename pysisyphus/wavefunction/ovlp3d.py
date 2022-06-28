import numpy


def ovlp3d_00(a, A, b, B):
    """Cartesian 3D (ss) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = a * b * x0

    # 1 item(s)
    S = numpy.array(
        [
            numpy.pi ** (3 / 2)
            * x0 ** (3 / 2)
            * numpy.exp(-x1 * (A[0] - B[0]) ** 2)
            * numpy.exp(-x1 * (A[1] - B[1]) ** 2)
            * numpy.exp(-x1 * (A[2] - B[2]) ** 2)
        ]
    )
    return S


def ovlp3d_01(a, A, b, B):
    """Cartesian 3D (sp) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = a * b * x0
    x2 = (
        numpy.pi ** (3 / 2)
        * x0 ** (3 / 2)
        * numpy.exp(-x1 * (A[0] - B[0]) ** 2)
        * numpy.exp(-x1 * (A[1] - B[1]) ** 2)
        * numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    )

    # 3 item(s)
    S = numpy.array(
        [
            x2 * (x0 * (a * A[0] + b * B[0]) - B[0]),
            x2 * (x0 * (a * A[1] + b * B[1]) - B[1]),
            x2 * (x0 * (a * A[2] + b * B[2]) - B[2]),
        ]
    )
    return S


def ovlp3d_02(a, A, b, B):
    """Cartesian 3D (sd) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = a * b * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x5 = numpy.sqrt(x0)
    x6 = numpy.sqrt(numpy.pi) * x5
    x7 = x4 * x6
    x8 = x0 * (a * A[0] + b * B[0]) - B[0]
    x9 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x10 = numpy.pi * x0 * x9
    x11 = x0 * (a * A[1] + b * B[1]) - B[1]
    x12 = numpy.pi ** (3 / 2)
    x13 = x0 * x2 * x4
    x14 = x12 * x13 * x5 * x8 * x9
    x15 = x0 * (a * A[2] + b * B[2]) - B[2]
    x16 = x2 * x6
    x17 = x6 * x9

    # 6 item(s)
    S = numpy.array(
        [
            x10 * x2 * (x3 * x7 + x7 * x8**2),
            x11 * x14,
            x14 * x15,
            x10 * x4 * (x11**2 * x16 + x16 * x3),
            x11 * x12 * x13 * x15 * x5 * x9,
            numpy.pi * x13 * (x15**2 * x17 + x17 * x3),
        ]
    )
    return S


def ovlp3d_03(a, A, b, B):
    """Cartesian 3D (sf) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = x0 * (a * A[0] + b * B[0]) - B[0]
    x2 = (2 * a + 2 * b) ** (-1.0)
    x3 = a * b * x0
    x4 = numpy.exp(-x3 * (A[0] - B[0]) ** 2)
    x5 = numpy.sqrt(x0)
    x6 = numpy.sqrt(numpy.pi) * x5
    x7 = x4 * x6
    x8 = x2 * x7
    x9 = x1**2 * x7 + x8
    x10 = numpy.exp(-x3 * (A[1] - B[1]) ** 2)
    x11 = numpy.exp(-x3 * (A[2] - B[2]) ** 2)
    x12 = numpy.pi * x0 * x11
    x13 = x10 * x12
    x14 = x0 * (a * A[1] + b * B[1]) - B[1]
    x15 = x13 * x9
    x16 = x0 * (a * A[2] + b * B[2]) - B[2]
    x17 = x10 * x6
    x18 = x17 * x2
    x19 = x14**2 * x17 + x18
    x20 = x12 * x4
    x21 = x19 * x20
    x22 = x0 * x10 * x4
    x23 = x11 * x6
    x24 = x2 * x23
    x25 = x16**2 * x23 + x24
    x26 = numpy.pi * x22
    x27 = x25 * x26

    # 10 item(s)
    S = numpy.array(
        [
            x13 * (2 * x1 * x8 + x1 * x9),
            x14 * x15,
            x15 * x16,
            x1 * x21,
            numpy.pi ** (3 / 2) * x1 * x11 * x14 * x16 * x22 * x5,
            x1 * x27,
            x20 * (2 * x14 * x18 + x14 * x19),
            x16 * x21,
            x14 * x27,
            x26 * (2 * x16 * x24 + x16 * x25),
        ]
    )
    return S


def ovlp3d_04(a, A, b, B):
    """Cartesian 3D (sg) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = a * b * x1
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x5 = x3 * x4
    x6 = x0 * x5
    x7 = x1 * (a * A[0] + b * B[0]) - B[0]
    x8 = x5 * x7**2
    x9 = x6 + x8
    x10 = 2 * x6 * x7 + x7 * x9
    x11 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x12 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x13 = numpy.pi * x1 * x12
    x14 = x11 * x13
    x15 = x1 * (a * A[1] + b * B[1]) - B[1]
    x16 = x10 * x14
    x17 = x1 * (a * A[2] + b * B[2]) - B[2]
    x18 = x11 * x4
    x19 = x0 * x18
    x20 = x15**2 * x18
    x21 = x19 + x20
    x22 = x12 * x4
    x23 = x0 * x22
    x24 = x17**2 * x22
    x25 = x23 + x24
    x26 = 2 * x15 * x19 + x15 * x21
    x27 = x13 * x3
    x28 = x26 * x27
    x29 = numpy.pi * x1 * x11 * x3
    x30 = 2 * x17 * x23 + x17 * x25
    x31 = x29 * x30

    # 15 item(s)
    S = numpy.array(
        [
            x14 * (x0 * (3 * x6 + 3 * x8) + x10 * x7),
            x15 * x16,
            x16 * x17,
            x21 * x22 * x9,
            x14 * x15 * x17 * x9,
            x18 * x25 * x9,
            x28 * x7,
            x17 * x21 * x27 * x7,
            x15 * x25 * x29 * x7,
            x31 * x7,
            x27 * (x0 * (3 * x19 + 3 * x20) + x15 * x26),
            x17 * x28,
            x21 * x25 * x5,
            x15 * x31,
            x29 * (x0 * (3 * x23 + 3 * x24) + x17 * x30),
        ]
    )
    return S


def ovlp3d_10(a, A, b, B):
    """Cartesian 3D (ps) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = a * b * x0
    x2 = (
        numpy.pi ** (3 / 2)
        * x0 ** (3 / 2)
        * numpy.exp(-x1 * (A[0] - B[0]) ** 2)
        * numpy.exp(-x1 * (A[1] - B[1]) ** 2)
        * numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    )

    # 3 item(s)
    S = numpy.array(
        [
            x2 * (x0 * (a * A[0] + b * B[0]) - A[0]),
            x2 * (x0 * (a * A[1] + b * B[1]) - A[1]),
            x2 * (x0 * (a * A[2] + b * B[2]) - A[2]),
        ]
    )
    return S


def ovlp3d_11(a, A, b, B):
    """Cartesian 3D (pp) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = a * b * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(x0)
    x5 = numpy.sqrt(numpy.pi) * x4
    x6 = x5 / (2 * a + 2 * b)
    x7 = -x0 * (a * A[0] + b * B[0])
    x8 = -x7 - B[0]
    x9 = x3 * (-x7 - A[0])
    x10 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x11 = numpy.pi * x0 * x10
    x12 = numpy.pi ** (3 / 2) * x4
    x13 = x12 * x9
    x14 = -x0 * (a * A[1] + b * B[1])
    x15 = x2 * (-x14 - B[1])
    x16 = x0 * x10 * x15
    x17 = -x0 * (a * A[2] + b * B[2])
    x18 = x10 * (-x17 - B[2])
    x19 = x0 * x2
    x20 = -x14 - A[1]
    x21 = x19 * x3
    x22 = x12 * x20 * x21
    x23 = x10 * x8
    x24 = -x17 - A[2]
    x25 = x12 * x24

    # 9 item(s)
    S = numpy.array(
        [
            x11 * x2 * (x3 * x6 + x5 * x8 * x9),
            x13 * x16,
            x13 * x18 * x19,
            x22 * x23,
            x11 * x3 * (x15 * x20 * x5 + x2 * x6),
            x18 * x22,
            x21 * x23 * x25,
            x16 * x25 * x3,
            numpy.pi * x21 * (x10 * x6 + x18 * x24 * x5),
        ]
    )
    return S


def ovlp3d_12(a, A, b, B):
    """Cartesian 3D (pd) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - B[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = a * b * x0
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x0)
    x7 = numpy.sqrt(numpy.pi) * x6
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x1 - A[0]
    x11 = x2**2 * x8 + x9
    x12 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x13 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x14 = numpy.pi * x0 * x13
    x15 = x12 * x14
    x16 = -x0 * (a * A[1] + b * B[1])
    x17 = -x16 - B[1]
    x18 = x15 * (x10 * x2 * x8 + x9)
    x19 = -x0 * (a * A[2] + b * B[2])
    x20 = -x19 - B[2]
    x21 = x12 * x7
    x22 = x21 * x3
    x23 = x17**2 * x21 + x22
    x24 = x14 * x5
    x25 = x23 * x24
    x26 = numpy.pi ** (3 / 2)
    x27 = x0 * x12 * x5
    x28 = x13 * x20 * x26 * x27 * x6
    x29 = x13 * x7
    x30 = x29 * x3
    x31 = x20**2 * x29 + x30
    x32 = numpy.pi * x27
    x33 = x31 * x32
    x34 = -x16 - A[1]
    x35 = x11 * x15
    x36 = x24 * (x17 * x21 * x34 + x22)
    x37 = -x19 - A[2]
    x38 = x32 * (x20 * x29 * x37 + x30)

    # 18 item(s)
    S = numpy.array(
        [
            x15 * (x10 * x11 + 2 * x2 * x9),
            x17 * x18,
            x18 * x20,
            x10 * x25,
            x10 * x17 * x28,
            x10 * x33,
            x34 * x35,
            x2 * x36,
            x2 * x28 * x34,
            x24 * (2 * x17 * x22 + x23 * x34),
            x20 * x36,
            x33 * x34,
            x35 * x37,
            x13 * x17 * x2 * x26 * x27 * x37 * x6,
            x2 * x38,
            x25 * x37,
            x17 * x38,
            x32 * (2 * x20 * x30 + x31 * x37),
        ]
    )
    return S


def ovlp3d_13(a, A, b, B):
    """Cartesian 3D (pf) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = a * b * x1
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x5 = x3 * x4
    x6 = x0 * x5
    x7 = -x1 * (a * A[0] + b * B[0])
    x8 = -x7 - B[0]
    x9 = x5 * x8**2
    x10 = -x7 - A[0]
    x11 = 2 * x6 * x8
    x12 = x6 + x9
    x13 = x11 + x12 * x8
    x14 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x15 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x16 = numpy.pi * x1 * x15
    x17 = x14 * x16
    x18 = -x1 * (a * A[1] + b * B[1])
    x19 = -x18 - B[1]
    x20 = x17 * (x10 * x12 + x11)
    x21 = -x1 * (a * A[2] + b * B[2])
    x22 = -x21 - B[2]
    x23 = x14 * x4
    x24 = x0 * x23
    x25 = x19**2 * x23
    x26 = x24 + x25
    x27 = x10 * x5 * x8 + x6
    x28 = x15 * x4
    x29 = x17 * x22
    x30 = x0 * x28
    x31 = x22**2 * x28
    x32 = x30 + x31
    x33 = 2 * x19 * x24
    x34 = x19 * x26 + x33
    x35 = x16 * x3
    x36 = x34 * x35
    x37 = x22 * x35
    x38 = numpy.pi * x1 * x14 * x3
    x39 = x32 * x38
    x40 = 2 * x22 * x30
    x41 = x22 * x32 + x40
    x42 = x38 * x41
    x43 = -x18 - A[1]
    x44 = x13 * x17
    x45 = x19 * x23 * x43 + x24
    x46 = x35 * (x26 * x43 + x33)
    x47 = -x21 - A[2]
    x48 = x22 * x28 * x47 + x30
    x49 = x38 * (x32 * x47 + x40)

    # 30 item(s)
    S = numpy.array(
        [
            x17 * (x0 * (3 * x6 + 3 * x9) + x10 * x13),
            x19 * x20,
            x20 * x22,
            x26 * x27 * x28,
            x19 * x27 * x29,
            x23 * x27 * x32,
            x10 * x36,
            x10 * x26 * x37,
            x10 * x19 * x39,
            x10 * x42,
            x43 * x44,
            x12 * x28 * x45,
            x12 * x29 * x43,
            x46 * x8,
            x37 * x45 * x8,
            x39 * x43 * x8,
            x35 * (x0 * (3 * x24 + 3 * x25) + x34 * x43),
            x22 * x46,
            x32 * x45 * x5,
            x42 * x43,
            x44 * x47,
            x12 * x17 * x19 * x47,
            x12 * x23 * x48,
            x26 * x35 * x47 * x8,
            x19 * x38 * x48 * x8,
            x49 * x8,
            x36 * x47,
            x26 * x48 * x5,
            x19 * x49,
            x38 * (x0 * (3 * x30 + 3 * x31) + x41 * x47),
        ]
    )
    return S


def ovlp3d_14(a, A, b, B):
    """Cartesian 3D (pg) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x0 * x7
    x9 = x3 * x8
    x10 = x3**2 * x7
    x11 = x10 + x8
    x12 = x11 * x3
    x13 = -x2 - A[0]
    x14 = x0 * (3 * x10 + 3 * x8)
    x15 = 2 * x9
    x16 = x12 + x15
    x17 = x14 + x16 * x3
    x18 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x19 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x20 = numpy.pi * x1 * x19
    x21 = x18 * x20
    x22 = -x1 * (a * A[1] + b * B[1])
    x23 = -x22 - B[1]
    x24 = x21 * (x13 * x16 + x14)
    x25 = -x1 * (a * A[2] + b * B[2])
    x26 = -x25 - B[2]
    x27 = x11 * x13 + x15
    x28 = x18 * x6
    x29 = x0 * x28
    x30 = x23**2 * x28
    x31 = x29 + x30
    x32 = x19 * x6
    x33 = x31 * x32
    x34 = x21 * x26
    x35 = x0 * x32
    x36 = x26**2 * x32
    x37 = x35 + x36
    x38 = x28 * x37
    x39 = x23 * x29
    x40 = 2 * x39
    x41 = x23 * x31
    x42 = x40 + x41
    x43 = x13 * x7
    x44 = x3 * x43 + x8
    x45 = x26 * x35
    x46 = 2 * x45
    x47 = x26 * x37
    x48 = x46 + x47
    x49 = x0 * (3 * x29 + 3 * x30)
    x50 = x23 * x42 + x49
    x51 = x20 * x5
    x52 = x50 * x51
    x53 = x26 * x51
    x54 = numpy.pi * x1 * x18 * x5
    x55 = x48 * x54
    x56 = x0 * (3 * x35 + 3 * x36)
    x57 = x26 * x48 + x56
    x58 = x54 * x57
    x59 = -x22 - A[1]
    x60 = x17 * x21
    x61 = x23 * x28
    x62 = x29 + x59 * x61
    x63 = x32 * x62
    x64 = x31 * x59 + x40
    x65 = x51 * (x42 * x59 + x49)
    x66 = x37 * x7
    x67 = -x25 - A[2]
    x68 = x26 * x32 * x67 + x35
    x69 = x37 * x67 + x46
    x70 = x68 * x7
    x71 = x54 * (x48 * x67 + x56)

    # 45 item(s)
    S = numpy.array(
        [
            x21 * (x0 * (4 * x12 + 8 * x9) + x13 * x17),
            x23 * x24,
            x24 * x26,
            x27 * x33,
            x23 * x27 * x34,
            x27 * x38,
            x32 * x42 * x44,
            x26 * x33 * x44,
            x23 * x38 * x44,
            x28 * x44 * x48,
            x13 * x52,
            x13 * x42 * x53,
            x31 * x37 * x43,
            x13 * x23 * x55,
            x13 * x58,
            x59 * x60,
            x16 * x63,
            x16 * x34 * x59,
            x11 * x32 * x64,
            x11 * x26 * x63,
            x11 * x38 * x59,
            x3 * x65,
            x3 * x53 * x64,
            x3 * x62 * x66,
            x3 * x55 * x59,
            x51 * (x0 * (8 * x39 + 4 * x41) + x50 * x59),
            x26 * x65,
            x64 * x66,
            x48 * x62 * x7,
            x58 * x59,
            x60 * x67,
            x16 * x21 * x23 * x67,
            x16 * x28 * x68,
            x11 * x33 * x67,
            x11 * x61 * x68,
            x11 * x28 * x69,
            x3 * x42 * x51 * x67,
            x3 * x31 * x70,
            x23 * x3 * x54 * x69,
            x3 * x71,
            x52 * x67,
            x42 * x70,
            x31 * x69 * x7,
            x23 * x71,
            x54 * (x0 * (8 * x45 + 4 * x47) + x57 * x67),
        ]
    )
    return S


def ovlp3d_20(a, A, b, B):
    """Cartesian 3D (ds) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = a * b * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x5 = numpy.sqrt(x0)
    x6 = numpy.sqrt(numpy.pi) * x5
    x7 = x4 * x6
    x8 = x0 * (a * A[0] + b * B[0]) - A[0]
    x9 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x10 = numpy.pi * x0 * x9
    x11 = x0 * (a * A[1] + b * B[1]) - A[1]
    x12 = numpy.pi ** (3 / 2)
    x13 = x0 * x2 * x4
    x14 = x12 * x13 * x5 * x8 * x9
    x15 = x0 * (a * A[2] + b * B[2]) - A[2]
    x16 = x2 * x6
    x17 = x6 * x9

    # 6 item(s)
    S = numpy.array(
        [
            x10 * x2 * (x3 * x7 + x7 * x8**2),
            x11 * x14,
            x14 * x15,
            x10 * x4 * (x11**2 * x16 + x16 * x3),
            x11 * x12 * x13 * x15 * x5 * x9,
            numpy.pi * x13 * (x15**2 * x17 + x17 * x3),
        ]
    )
    return S


def ovlp3d_21(a, A, b, B):
    """Cartesian 3D (dp) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = numpy.sqrt(x1)
    x3 = numpy.sqrt(numpy.pi) * x2
    x4 = -x1 * (a * A[0] + b * B[0])
    x5 = -x4 - A[0]
    x6 = a * b * x1
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x4 - B[0]
    x11 = x3 * x7
    x12 = x0 * x11
    x13 = x10 * x9 + x12
    x14 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x15 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x16 = numpy.pi * x1 * x15
    x17 = x14 * x16
    x18 = -x1 * (a * A[1] + b * B[1])
    x19 = -x18 - B[1]
    x20 = x17 * (x11 * x5**2 + x12)
    x21 = -x1 * (a * A[2] + b * B[2])
    x22 = -x21 - B[2]
    x23 = -x18 - A[1]
    x24 = x13 * x17
    x25 = x0 * x3
    x26 = x14 * x25
    x27 = x14 * x3
    x28 = x23 * x27
    x29 = x19 * x28 + x26
    x30 = numpy.pi ** (3 / 2)
    x31 = x1 * x14
    x32 = x15 * x2 * x30 * x31 * x8
    x33 = -x21 - A[2]
    x34 = x15 * x25
    x35 = x15 * x3
    x36 = x33 * x35
    x37 = x22 * x36 + x34
    x38 = numpy.pi * x31
    x39 = x16 * x7
    x40 = x39 * (x23**2 * x27 + x26)
    x41 = x38 * x7
    x42 = x41 * (x33**2 * x35 + x34)

    # 18 item(s)
    S = numpy.array(
        [
            x17 * (x0 * (x10 * x11 + x9) + x13 * x5),
            x19 * x20,
            x20 * x22,
            x23 * x24,
            x16 * x29 * x8,
            x22 * x23 * x32,
            x24 * x33,
            x19 * x32 * x33,
            x37 * x38 * x8,
            x10 * x40,
            x39 * (x0 * (x19 * x27 + x28) + x23 * x29),
            x22 * x40,
            x10 * x15 * x2 * x23 * x30 * x31 * x33 * x7,
            x29 * x33 * x39,
            x23 * x37 * x41,
            x10 * x42,
            x19 * x42,
            x41 * (x0 * (x22 * x35 + x36) + x33 * x37),
        ]
    )
    return S


def ovlp3d_22(a, A, b, B):
    """Cartesian 3D (dd) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = a * b * x1
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x5 = x3 * x4
    x6 = x0 * x5
    x7 = -x1 * (a * A[0] + b * B[0])
    x8 = -x7 - B[0]
    x9 = x5 * x8**2
    x10 = -x7 - A[0]
    x11 = x10 * x5
    x12 = x11 * x8
    x13 = x6 + x9
    x14 = x10 * x13 + 2 * x6 * x8
    x15 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x16 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x17 = numpy.pi * x1 * x16
    x18 = x15 * x17
    x19 = -x1 * (a * A[1] + b * B[1])
    x20 = -x19 - B[1]
    x21 = x12 + x6
    x22 = x18 * (x0 * (x11 + x5 * x8) + x10 * x21)
    x23 = -x1 * (a * A[2] + b * B[2])
    x24 = -x23 - B[2]
    x25 = x15 * x4
    x26 = x0 * x25
    x27 = x20**2 * x25
    x28 = x26 + x27
    x29 = x10**2 * x5 + x6
    x30 = x16 * x4
    x31 = x18 * x24
    x32 = x0 * x30
    x33 = x24**2 * x30
    x34 = x32 + x33
    x35 = -x19 - A[1]
    x36 = x14 * x18
    x37 = x25 * x35
    x38 = x20 * x37
    x39 = x26 + x38
    x40 = 2 * x20 * x26 + x28 * x35
    x41 = x17 * x3
    x42 = x40 * x41
    x43 = x10 * x41
    x44 = numpy.pi * x1 * x15 * x3
    x45 = x10 * x44
    x46 = -x23 - A[2]
    x47 = x18 * x46
    x48 = x30 * x46
    x49 = x24 * x48
    x50 = x32 + x49
    x51 = 2 * x24 * x32 + x34 * x46
    x52 = x44 * x51
    x53 = x25 * x35**2 + x26
    x54 = x41 * (x0 * (x20 * x25 + x37) + x35 * x39)
    x55 = x41 * x8
    x56 = x44 * x8
    x57 = x30 * x46**2 + x32
    x58 = x44 * (x0 * (x24 * x30 + x48) + x46 * x50)

    # 36 item(s)
    S = numpy.array(
        [
            x18 * (x0 * (2 * x12 + 3 * x6 + x9) + x10 * x14),
            x20 * x22,
            x22 * x24,
            x28 * x29 * x30,
            x20 * x29 * x31,
            x25 * x29 * x34,
            x35 * x36,
            x21 * x30 * x39,
            x21 * x31 * x35,
            x10 * x42,
            x24 * x39 * x43,
            x34 * x35 * x45,
            x36 * x46,
            x20 * x21 * x47,
            x21 * x25 * x50,
            x28 * x43 * x46,
            x20 * x45 * x50,
            x10 * x52,
            x13 * x30 * x53,
            x54 * x8,
            x24 * x53 * x55,
            x41 * (x0 * (3 * x26 + x27 + 2 * x38) + x35 * x40),
            x24 * x54,
            x34 * x5 * x53,
            x13 * x35 * x47,
            x39 * x46 * x55,
            x35 * x50 * x56,
            x42 * x46,
            x39 * x5 * x50,
            x35 * x52,
            x13 * x25 * x57,
            x20 * x56 * x57,
            x58 * x8,
            x28 * x5 * x57,
            x20 * x58,
            x44 * (x0 * (3 * x32 + x33 + 2 * x49) + x46 * x51),
        ]
    )
    return S


def ovlp3d_23(a, A, b, B):
    """Cartesian 3D (df) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x0 * x7
    x9 = x3 * x8
    x10 = x3**2 * x7
    x11 = x10 + x8
    x12 = x11 * x3
    x13 = -x2 - A[0]
    x14 = x11 * x13
    x15 = 3 * x8
    x16 = 2 * x9
    x17 = x12 + x16
    x18 = x0 * (3 * x10 + x15) + x13 * x17
    x19 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x20 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x21 = numpy.pi * x1 * x20
    x22 = x19 * x21
    x23 = -x1 * (a * A[1] + b * B[1])
    x24 = -x23 - B[1]
    x25 = x13 * x7
    x26 = x25 * x3
    x27 = x14 + x16
    x28 = x22 * (x0 * (x10 + x15 + 2 * x26) + x13 * x27)
    x29 = -x1 * (a * A[2] + b * B[2])
    x30 = -x29 - B[2]
    x31 = x19 * x6
    x32 = x0 * x31
    x33 = x24**2 * x31
    x34 = x32 + x33
    x35 = x3 * x7
    x36 = x26 + x8
    x37 = x0 * (x25 + x35) + x13 * x36
    x38 = x20 * x6
    x39 = x22 * x30
    x40 = x0 * x38
    x41 = x30**2 * x38
    x42 = x40 + x41
    x43 = x24 * x32
    x44 = 2 * x43
    x45 = x24 * x34
    x46 = x44 + x45
    x47 = x13**2 * x7 + x8
    x48 = x30 * x38
    x49 = x24 * x31
    x50 = x30 * x40
    x51 = 2 * x50
    x52 = x30 * x42
    x53 = x51 + x52
    x54 = -x23 - A[1]
    x55 = x18 * x22
    x56 = x31 * x54
    x57 = x24 * x56
    x58 = x32 + x57
    x59 = x34 * x54
    x60 = x44 + x59
    x61 = 3 * x32
    x62 = x0 * (3 * x33 + x61) + x46 * x54
    x63 = x21 * x5
    x64 = x62 * x63
    x65 = x13 * x63
    x66 = numpy.pi * x1 * x19 * x5
    x67 = x13 * x66
    x68 = -x29 - A[2]
    x69 = x22 * x68
    x70 = x38 * x68
    x71 = x30 * x70
    x72 = x40 + x71
    x73 = x42 * x68
    x74 = x51 + x73
    x75 = 3 * x40
    x76 = x0 * (3 * x41 + x75) + x53 * x68
    x77 = x66 * x76
    x78 = x31 * x54**2 + x32
    x79 = x0 * (x49 + x56) + x54 * x58
    x80 = x63 * (x0 * (x33 + 2 * x57 + x61) + x54 * x60)
    x81 = x3 * x63
    x82 = x3 * x66
    x83 = x38 * x68**2 + x40
    x84 = x0 * (x48 + x70) + x68 * x72
    x85 = x66 * (x0 * (x41 + 2 * x71 + x75) + x68 * x74)

    # 60 item(s)
    S = numpy.array(
        [
            x22 * (x0 * (x12 + 3 * x14 + 8 * x9) + x13 * x18),
            x24 * x28,
            x28 * x30,
            x34 * x37 * x38,
            x24 * x37 * x39,
            x31 * x37 * x42,
            x38 * x46 * x47,
            x34 * x47 * x48,
            x42 * x47 * x49,
            x31 * x47 * x53,
            x54 * x55,
            x27 * x38 * x58,
            x27 * x39 * x54,
            x36 * x38 * x60,
            x36 * x48 * x58,
            x36 * x42 * x56,
            x13 * x64,
            x30 * x60 * x65,
            x25 * x42 * x58,
            x53 * x54 * x67,
            x55 * x68,
            x24 * x27 * x69,
            x27 * x31 * x72,
            x34 * x36 * x70,
            x36 * x49 * x72,
            x31 * x36 * x74,
            x46 * x65 * x68,
            x25 * x34 * x72,
            x24 * x67 * x74,
            x13 * x77,
            x17 * x38 * x78,
            x11 * x38 * x79,
            x11 * x48 * x78,
            x3 * x80,
            x30 * x79 * x81,
            x35 * x42 * x78,
            x63 * (x0 * (8 * x43 + x45 + 3 * x59) + x54 * x62),
            x30 * x80,
            x42 * x7 * x79,
            x53 * x7 * x78,
            x17 * x54 * x69,
            x11 * x58 * x70,
            x11 * x56 * x72,
            x60 * x68 * x81,
            x35 * x58 * x72,
            x54 * x74 * x82,
            x64 * x68,
            x60 * x7 * x72,
            x58 * x7 * x74,
            x54 * x77,
            x17 * x31 * x83,
            x11 * x49 * x83,
            x11 * x31 * x84,
            x34 * x35 * x83,
            x24 * x82 * x84,
            x3 * x85,
            x46 * x7 * x83,
            x34 * x7 * x84,
            x24 * x85,
            x66 * (x0 * (8 * x50 + x52 + 3 * x73) + x68 * x76),
        ]
    )
    return S


def ovlp3d_24(a, A, b, B):
    """Cartesian 3D (dg) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = a * b * x1
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x5 = x3 * x4
    x6 = x0 * x5
    x7 = 3 * x6
    x8 = -x1 * (a * A[0] + b * B[0])
    x9 = -x8 - B[0]
    x10 = x5 * x9**2
    x11 = x0 * (3 * x10 + x7)
    x12 = x6 * x9
    x13 = 2 * x12
    x14 = x10 + x6
    x15 = x14 * x9
    x16 = x13 + x15
    x17 = x16 * x9
    x18 = -x8 - A[0]
    x19 = x16 * x18
    x20 = 8 * x12
    x21 = x11 + x17
    x22 = x0 * (4 * x15 + x20) + x18 * x21
    x23 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x24 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x25 = numpy.pi * x1 * x24
    x26 = x23 * x25
    x27 = -x1 * (a * A[1] + b * B[1])
    x28 = -x27 - B[1]
    x29 = x14 * x18
    x30 = x11 + x19
    x31 = x26 * (x0 * (x15 + x20 + 3 * x29) + x18 * x30)
    x32 = -x1 * (a * A[2] + b * B[2])
    x33 = -x32 - B[2]
    x34 = x23 * x4
    x35 = x0 * x34
    x36 = x28**2 * x34
    x37 = x35 + x36
    x38 = x18 * x5
    x39 = x38 * x9
    x40 = x13 + x29
    x41 = x0 * (x10 + 2 * x39 + x7) + x18 * x40
    x42 = x24 * x4
    x43 = x26 * x33
    x44 = x0 * x42
    x45 = x33**2 * x42
    x46 = x44 + x45
    x47 = x28 * x35
    x48 = 2 * x47
    x49 = x28 * x37
    x50 = x48 + x49
    x51 = x5 * x9
    x52 = x39 + x6
    x53 = x0 * (x38 + x51) + x18 * x52
    x54 = x33 * x42
    x55 = x28 * x34
    x56 = x33 * x44
    x57 = 2 * x56
    x58 = x33 * x46
    x59 = x57 + x58
    x60 = 3 * x35
    x61 = x0 * (3 * x36 + x60)
    x62 = x28 * x50
    x63 = x61 + x62
    x64 = x18**2 * x5 + x6
    x65 = 3 * x44
    x66 = x0 * (3 * x45 + x65)
    x67 = x33 * x59
    x68 = x66 + x67
    x69 = -x27 - A[1]
    x70 = x22 * x26
    x71 = x34 * x69
    x72 = x28 * x71
    x73 = x35 + x72
    x74 = x37 * x69
    x75 = x48 + x74
    x76 = x50 * x69
    x77 = x61 + x76
    x78 = 8 * x47
    x79 = x0 * (4 * x49 + x78) + x63 * x69
    x80 = x25 * x3
    x81 = x79 * x80
    x82 = x18 * x80
    x83 = numpy.pi * x1 * x23 * x3
    x84 = x18 * x83
    x85 = -x32 - A[2]
    x86 = x26 * x85
    x87 = x42 * x85
    x88 = x33 * x87
    x89 = x44 + x88
    x90 = x46 * x85
    x91 = x57 + x90
    x92 = x59 * x85
    x93 = x66 + x92
    x94 = 8 * x56
    x95 = x0 * (4 * x58 + x94) + x68 * x85
    x96 = x83 * x95
    x97 = x34 * x69**2 + x35
    x98 = x0 * (x55 + x71) + x69 * x73
    x99 = x0 * (x36 + x60 + 2 * x72) + x69 * x75
    x100 = x80 * (x0 * (x49 + 3 * x74 + x78) + x69 * x77)
    x101 = x80 * x9
    x102 = x83 * x9
    x103 = x42 * x85**2 + x44
    x104 = x0 * (x54 + x87) + x85 * x89
    x105 = x0 * (x45 + x65 + 2 * x88) + x85 * x91
    x106 = x83 * (x0 * (x58 + 3 * x90 + x94) + x85 * x93)

    # 90 item(s)
    S = numpy.array(
        [
            x26 * (x0 * (5 * x11 + x17 + 4 * x19) + x18 * x22),
            x28 * x31,
            x31 * x33,
            x37 * x41 * x42,
            x28 * x41 * x43,
            x34 * x41 * x46,
            x42 * x50 * x53,
            x37 * x53 * x54,
            x46 * x53 * x55,
            x34 * x53 * x59,
            x42 * x63 * x64,
            x50 * x54 * x64,
            x37 * x46 * x64,
            x55 * x59 * x64,
            x34 * x64 * x68,
            x69 * x70,
            x30 * x42 * x73,
            x30 * x43 * x69,
            x40 * x42 * x75,
            x40 * x54 * x73,
            x40 * x46 * x71,
            x42 * x52 * x77,
            x52 * x54 * x75,
            x46 * x52 * x73,
            x52 * x59 * x71,
            x18 * x81,
            x33 * x77 * x82,
            x38 * x46 * x75,
            x38 * x59 * x73,
            x68 * x69 * x84,
            x70 * x85,
            x28 * x30 * x86,
            x30 * x34 * x89,
            x37 * x40 * x87,
            x40 * x55 * x89,
            x34 * x40 * x91,
            x50 * x52 * x87,
            x37 * x52 * x89,
            x52 * x55 * x91,
            x34 * x52 * x93,
            x63 * x82 * x85,
            x38 * x50 * x89,
            x37 * x38 * x91,
            x28 * x84 * x93,
            x18 * x96,
            x21 * x42 * x97,
            x16 * x42 * x98,
            x16 * x54 * x97,
            x14 * x42 * x99,
            x14 * x54 * x98,
            x14 * x46 * x97,
            x100 * x9,
            x101 * x33 * x99,
            x46 * x51 * x98,
            x51 * x59 * x97,
            x80 * (x0 * (5 * x61 + x62 + 4 * x76) + x69 * x79),
            x100 * x33,
            x46 * x5 * x99,
            x5 * x59 * x98,
            x5 * x68 * x97,
            x21 * x69 * x86,
            x16 * x73 * x87,
            x16 * x71 * x89,
            x14 * x75 * x87,
            x14 * x73 * x89,
            x14 * x71 * x91,
            x101 * x77 * x85,
            x51 * x75 * x89,
            x51 * x73 * x91,
            x102 * x69 * x93,
            x81 * x85,
            x5 * x77 * x89,
            x5 * x75 * x91,
            x5 * x73 * x93,
            x69 * x96,
            x103 * x21 * x34,
            x103 * x16 * x55,
            x104 * x16 * x34,
            x103 * x14 * x37,
            x104 * x14 * x55,
            x105 * x14 * x34,
            x103 * x50 * x51,
            x104 * x37 * x51,
            x102 * x105 * x28,
            x106 * x9,
            x103 * x5 * x63,
            x104 * x5 * x50,
            x105 * x37 * x5,
            x106 * x28,
            x83 * (x0 * (5 * x66 + x67 + 4 * x92) + x85 * x95),
        ]
    )
    return S


def ovlp3d_30(a, A, b, B):
    """Cartesian 3D (fs) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = x0 * (a * A[0] + b * B[0]) - A[0]
    x2 = (2 * a + 2 * b) ** (-1.0)
    x3 = a * b * x0
    x4 = numpy.exp(-x3 * (A[0] - B[0]) ** 2)
    x5 = numpy.sqrt(x0)
    x6 = numpy.sqrt(numpy.pi) * x5
    x7 = x4 * x6
    x8 = x2 * x7
    x9 = x1**2 * x7 + x8
    x10 = numpy.exp(-x3 * (A[1] - B[1]) ** 2)
    x11 = numpy.exp(-x3 * (A[2] - B[2]) ** 2)
    x12 = numpy.pi * x0 * x11
    x13 = x10 * x12
    x14 = x0 * (a * A[1] + b * B[1]) - A[1]
    x15 = x13 * x9
    x16 = x0 * (a * A[2] + b * B[2]) - A[2]
    x17 = x10 * x6
    x18 = x17 * x2
    x19 = x14**2 * x17 + x18
    x20 = x12 * x4
    x21 = x19 * x20
    x22 = x0 * x10 * x4
    x23 = x11 * x6
    x24 = x2 * x23
    x25 = x16**2 * x23 + x24
    x26 = numpy.pi * x22
    x27 = x25 * x26

    # 10 item(s)
    S = numpy.array(
        [
            x13 * (2 * x1 * x8 + x1 * x9),
            x14 * x15,
            x15 * x16,
            x1 * x21,
            numpy.pi ** (3 / 2) * x1 * x11 * x14 * x16 * x22 * x5,
            x1 * x27,
            x20 * (2 * x14 * x18 + x14 * x19),
            x16 * x21,
            x14 * x27,
            x26 * (2 * x16 * x24 + x16 * x25),
        ]
    )
    return S


def ovlp3d_31(a, A, b, B):
    """Cartesian 3D (fp) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = a * b * x1
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x5 = x3 * x4
    x6 = x0 * x5
    x7 = -x1 * (a * A[0] + b * B[0])
    x8 = -x7 - A[0]
    x9 = x5 * x8**2
    x10 = -x7 - B[0]
    x11 = x5 * x8
    x12 = x10 * x11
    x13 = x12 + x6
    x14 = x0 * (x10 * x5 + x11) + x13 * x8
    x15 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x16 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x17 = numpy.pi * x1 * x16
    x18 = x15 * x17
    x19 = -x1 * (a * A[1] + b * B[1])
    x20 = -x19 - B[1]
    x21 = x6 + x9
    x22 = x18 * (2 * x0 * x11 + x21 * x8)
    x23 = -x1 * (a * A[2] + b * B[2])
    x24 = -x23 - B[2]
    x25 = -x19 - A[1]
    x26 = x14 * x18
    x27 = x0 * x4
    x28 = x15 * x27
    x29 = x15 * x4
    x30 = x25 * x29
    x31 = x20 * x30
    x32 = x28 + x31
    x33 = x16 * x4
    x34 = x18 * x21
    x35 = -x23 - A[2]
    x36 = x16 * x27
    x37 = x33 * x35
    x38 = x24 * x37
    x39 = x36 + x38
    x40 = x25**2 * x29
    x41 = x28 + x40
    x42 = x0 * (x20 * x29 + x30) + x25 * x32
    x43 = x17 * x3
    x44 = x42 * x43
    x45 = x43 * x8
    x46 = numpy.pi * x1 * x15 * x3
    x47 = x46 * x8
    x48 = x33 * x35**2
    x49 = x36 + x48
    x50 = x0 * (x24 * x33 + x37) + x35 * x39
    x51 = x46 * x50
    x52 = x43 * (2 * x25 * x28 + x25 * x41)
    x53 = x46 * (2 * x35 * x36 + x35 * x49)

    # 30 item(s)
    S = numpy.array(
        [
            x18 * (x0 * (2 * x12 + 3 * x6 + x9) + x14 * x8),
            x20 * x22,
            x22 * x24,
            x25 * x26,
            x21 * x32 * x33,
            x24 * x25 * x34,
            x26 * x35,
            x20 * x34 * x35,
            x21 * x29 * x39,
            x13 * x33 * x41,
            x44 * x8,
            x24 * x41 * x45,
            x13 * x18 * x25 * x35,
            x32 * x35 * x45,
            x25 * x39 * x47,
            x13 * x29 * x49,
            x20 * x47 * x49,
            x51 * x8,
            x10 * x52,
            x43 * (x0 * (3 * x28 + 2 * x31 + x40) + x25 * x42),
            x24 * x52,
            x10 * x35 * x41 * x43,
            x35 * x44,
            x39 * x41 * x5,
            x10 * x25 * x46 * x49,
            x32 * x49 * x5,
            x25 * x51,
            x10 * x53,
            x20 * x53,
            x46 * (x0 * (3 * x36 + 2 * x38 + x48) + x35 * x50),
        ]
    )
    return S


def ovlp3d_32(a, A, b, B):
    """Cartesian 3D (fd) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - B[0]
    x5 = a * b * x0
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x8 = x6 * x7
    x9 = x4**2 * x8
    x10 = x3 * x8
    x11 = x2 * x6
    x12 = x11 * x7
    x13 = x12 * x4
    x14 = 3 * x10 + 2 * x13
    x15 = x4 * x8
    x16 = x15 * x3
    x17 = x10 + x9
    x18 = x17 * x2
    x19 = 2 * x16 + x18
    x20 = x19 * x2 + x3 * (x14 + x9)
    x21 = x3 * (x12 + x15)
    x22 = x10 + x13
    x23 = x2 * x22
    x24 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x25 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x26 = numpy.pi * x0 * x25
    x27 = x24 * x26
    x28 = -x0 * (a * A[1] + b * B[1])
    x29 = -x28 - B[1]
    x30 = x2**2 * x8
    x31 = x21 + x23
    x32 = x27 * (x2 * x31 + x3 * (x14 + x30))
    x33 = -x0 * (a * A[2] + b * B[2])
    x34 = -x33 - B[2]
    x35 = x24 * x7
    x36 = x3 * x35
    x37 = x29**2 * x35
    x38 = x36 + x37
    x39 = x10 + x30
    x40 = 2 * x12 * x3 + x2 * x39
    x41 = x25 * x7
    x42 = x27 * x34
    x43 = x3 * x41
    x44 = x34**2 * x41
    x45 = x43 + x44
    x46 = -x28 - A[1]
    x47 = x20 * x27
    x48 = x35 * x46
    x49 = x29 * x48
    x50 = x36 + x49
    x51 = 2 * x36
    x52 = x38 * x46
    x53 = x29 * x51 + x52
    x54 = x34 * x41
    x55 = -x33 - A[2]
    x56 = x27 * x55
    x57 = x41 * x55
    x58 = x34 * x57
    x59 = x43 + x58
    x60 = x29 * x35
    x61 = 2 * x43
    x62 = x45 * x55
    x63 = x34 * x61 + x62
    x64 = x35 * x46**2
    x65 = x36 + x64
    x66 = x3 * (x48 + x60)
    x67 = x46 * x50
    x68 = x66 + x67
    x69 = 3 * x36 + 2 * x49
    x70 = x3 * (x37 + x69) + x46 * x53
    x71 = x11 * x26
    x72 = numpy.pi * x0 * x24
    x73 = x11 * x72
    x74 = x41 * x55**2
    x75 = x43 + x74
    x76 = x3 * (x54 + x57)
    x77 = x55 * x59
    x78 = x76 + x77
    x79 = 3 * x43 + 2 * x58
    x80 = x3 * (x44 + x79) + x55 * x63
    x81 = x46 * x51 + x46 * x65
    x82 = x26 * x6
    x83 = x82 * (x3 * (x64 + x69) + x46 * x68)
    x84 = x55 * x82
    x85 = x6 * x72
    x86 = x46 * x85
    x87 = x55 * x61 + x55 * x75
    x88 = x85 * (x3 * (x74 + x79) + x55 * x78)

    # 60 item(s)
    S = numpy.array(
        [
            x27 * (x2 * x20 + x3 * (4 * x16 + 2 * x18 + 2 * x21 + 2 * x23)),
            x29 * x32,
            x32 * x34,
            x38 * x40 * x41,
            x29 * x40 * x42,
            x35 * x40 * x45,
            x46 * x47,
            x31 * x41 * x50,
            x31 * x42 * x46,
            x39 * x41 * x53,
            x39 * x50 * x54,
            x39 * x45 * x48,
            x47 * x55,
            x29 * x31 * x56,
            x31 * x35 * x59,
            x38 * x39 * x57,
            x39 * x59 * x60,
            x35 * x39 * x63,
            x19 * x41 * x65,
            x22 * x41 * x68,
            x22 * x54 * x65,
            x70 * x71,
            x34 * x68 * x71,
            x12 * x45 * x65,
            x19 * x46 * x56,
            x22 * x50 * x57,
            x22 * x48 * x59,
            x53 * x55 * x71,
            x12 * x50 * x59,
            x46 * x63 * x73,
            x19 * x35 * x75,
            x22 * x60 * x75,
            x22 * x35 * x78,
            x12 * x38 * x75,
            x29 * x73 * x78,
            x73 * x80,
            x17 * x41 * x81,
            x4 * x83,
            x34 * x4 * x81 * x82,
            x82 * (x3 * (4 * x29 * x36 + 2 * x52 + 2 * x66 + 2 * x67) + x46 * x70),
            x34 * x83,
            x45 * x8 * x81,
            x17 * x57 * x65,
            x4 * x68 * x84,
            x15 * x59 * x65,
            x70 * x84,
            x59 * x68 * x8,
            x63 * x65 * x8,
            x17 * x48 * x75,
            x15 * x50 * x75,
            x4 * x78 * x86,
            x53 * x75 * x8,
            x50 * x78 * x8,
            x80 * x86,
            x17 * x35 * x87,
            x29 * x4 * x85 * x87,
            x4 * x88,
            x38 * x8 * x87,
            x29 * x88,
            x85 * (x3 * (4 * x34 * x43 + 2 * x62 + 2 * x76 + 2 * x77) + x55 * x80),
        ]
    )
    return S


def ovlp3d_33(a, A, b, B):
    """Cartesian 3D (ff) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - B[0]
    x5 = a * b * x0
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x8 = x6 * x7
    x9 = x3 * x8
    x10 = x4 * x9
    x11 = x4**2 * x8
    x12 = x11 + x9
    x13 = x12 * x4
    x14 = x12 * x2
    x15 = 3 * x9
    x16 = x3 * (3 * x11 + x15)
    x17 = 2 * x10
    x18 = x13 + x17
    x19 = x18 * x2
    x20 = x16 + x19
    x21 = x2 * x20 + x3 * (8 * x10 + x13 + 3 * x14)
    x22 = x2 * x8
    x23 = x22 * x4
    x24 = x15 + 2 * x23
    x25 = x3 * (x11 + x24)
    x26 = x14 + x17
    x27 = x2 * x26
    x28 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x29 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x30 = numpy.pi * x0 * x29
    x31 = x28 * x30
    x32 = -x0 * (a * A[1] + b * B[1])
    x33 = -x32 - B[1]
    x34 = x25 + x27
    x35 = x4 * x8
    x36 = x3 * (x22 + x35)
    x37 = x23 + x9
    x38 = x2 * x37
    x39 = x31 * (x2 * x34 + x3 * (4 * x10 + 2 * x14 + 2 * x36 + 2 * x38))
    x40 = -x0 * (a * A[2] + b * B[2])
    x41 = -x40 - B[2]
    x42 = x28 * x7
    x43 = x3 * x42
    x44 = x33**2 * x42
    x45 = x43 + x44
    x46 = x2**2 * x8
    x47 = x36 + x38
    x48 = x2 * x47 + x3 * (x24 + x46)
    x49 = x29 * x7
    x50 = x31 * x41
    x51 = x3 * x49
    x52 = x41**2 * x49
    x53 = x51 + x52
    x54 = x33 * x43
    x55 = 2 * x54
    x56 = x33 * x45
    x57 = x55 + x56
    x58 = x46 + x9
    x59 = x2 * x58 + 2 * x2 * x9
    x60 = x41 * x49
    x61 = x33 * x42
    x62 = x41 * x51
    x63 = 2 * x62
    x64 = x41 * x53
    x65 = x63 + x64
    x66 = -x32 - A[1]
    x67 = x21 * x31
    x68 = x42 * x66
    x69 = x33 * x68
    x70 = x43 + x69
    x71 = x45 * x66
    x72 = x55 + x71
    x73 = 3 * x43
    x74 = x3 * (3 * x44 + x73)
    x75 = x57 * x66
    x76 = x74 + x75
    x77 = -x40 - A[2]
    x78 = x31 * x77
    x79 = x49 * x77
    x80 = x41 * x79
    x81 = x51 + x80
    x82 = x53 * x77
    x83 = x63 + x82
    x84 = 3 * x51
    x85 = x3 * (3 * x52 + x84)
    x86 = x65 * x77
    x87 = x85 + x86
    x88 = x42 * x66**2
    x89 = x43 + x88
    x90 = x3 * (x61 + x68)
    x91 = x66 * x70
    x92 = x90 + x91
    x93 = 2 * x69 + x73
    x94 = x3 * (x44 + x93)
    x95 = x66 * x72
    x96 = x94 + x95
    x97 = x3 * (8 * x54 + x56 + 3 * x71) + x66 * x76
    x98 = x30 * x6
    x99 = x97 * x98
    x100 = x2 * x98
    x101 = numpy.pi * x0 * x28 * x6
    x102 = x101 * x2
    x103 = x49 * x77**2
    x104 = x103 + x51
    x105 = x3 * (x60 + x79)
    x106 = x77 * x81
    x107 = x105 + x106
    x108 = 2 * x80 + x84
    x109 = x3 * (x108 + x52)
    x110 = x77 * x83
    x111 = x109 + x110
    x112 = x3 * (8 * x62 + x64 + 3 * x82) + x77 * x87
    x113 = x101 * x112
    x114 = 2 * x43 * x66 + x66 * x89
    x115 = x3 * (x88 + x93) + x66 * x92
    x116 = x98 * (x3 * (4 * x54 + 2 * x71 + 2 * x90 + 2 * x91) + x66 * x96)
    x117 = x4 * x98
    x118 = x101 * x4
    x119 = x104 * x77 + 2 * x51 * x77
    x120 = x107 * x77 + x3 * (x103 + x108)
    x121 = x101 * (x111 * x77 + x3 * (2 * x105 + 2 * x106 + 4 * x62 + 2 * x82))

    # 100 item(s)
    S = numpy.array(
        [
            x31 * (x2 * x21 + x3 * (2 * x16 + 2 * x19 + 3 * x25 + 3 * x27)),
            x33 * x39,
            x39 * x41,
            x45 * x48 * x49,
            x33 * x48 * x50,
            x42 * x48 * x53,
            x49 * x57 * x59,
            x45 * x59 * x60,
            x53 * x59 * x61,
            x42 * x59 * x65,
            x66 * x67,
            x34 * x49 * x70,
            x34 * x50 * x66,
            x47 * x49 * x72,
            x47 * x60 * x70,
            x47 * x53 * x68,
            x49 * x58 * x76,
            x58 * x60 * x72,
            x53 * x58 * x70,
            x58 * x65 * x68,
            x67 * x77,
            x33 * x34 * x78,
            x34 * x42 * x81,
            x45 * x47 * x79,
            x47 * x61 * x81,
            x42 * x47 * x83,
            x57 * x58 * x79,
            x45 * x58 * x81,
            x58 * x61 * x83,
            x42 * x58 * x87,
            x20 * x49 * x89,
            x26 * x49 * x92,
            x26 * x60 * x89,
            x37 * x49 * x96,
            x37 * x60 * x92,
            x37 * x53 * x89,
            x2 * x99,
            x100 * x41 * x96,
            x22 * x53 * x92,
            x22 * x65 * x89,
            x20 * x66 * x78,
            x26 * x70 * x79,
            x26 * x68 * x81,
            x37 * x72 * x79,
            x37 * x70 * x81,
            x37 * x68 * x83,
            x100 * x76 * x77,
            x22 * x72 * x81,
            x22 * x70 * x83,
            x102 * x66 * x87,
            x104 * x20 * x42,
            x104 * x26 * x61,
            x107 * x26 * x42,
            x104 * x37 * x45,
            x107 * x37 * x61,
            x111 * x37 * x42,
            x104 * x22 * x57,
            x107 * x22 * x45,
            x102 * x111 * x33,
            x113 * x2,
            x114 * x18 * x49,
            x115 * x12 * x49,
            x114 * x12 * x60,
            x116 * x4,
            x115 * x117 * x41,
            x114 * x35 * x53,
            x98 * (x3 * (2 * x74 + 2 * x75 + 3 * x94 + 3 * x95) + x66 * x97),
            x116 * x41,
            x115 * x53 * x8,
            x114 * x65 * x8,
            x18 * x79 * x89,
            x12 * x79 * x92,
            x12 * x81 * x89,
            x117 * x77 * x96,
            x35 * x81 * x92,
            x35 * x83 * x89,
            x77 * x99,
            x8 * x81 * x96,
            x8 * x83 * x92,
            x8 * x87 * x89,
            x104 * x18 * x68,
            x104 * x12 * x70,
            x107 * x12 * x68,
            x104 * x35 * x72,
            x107 * x35 * x70,
            x111 * x118 * x66,
            x104 * x76 * x8,
            x107 * x72 * x8,
            x111 * x70 * x8,
            x113 * x66,
            x119 * x18 * x42,
            x119 * x12 * x61,
            x12 * x120 * x42,
            x119 * x35 * x45,
            x118 * x120 * x33,
            x121 * x4,
            x119 * x57 * x8,
            x120 * x45 * x8,
            x121 * x33,
            x101 * (x112 * x77 + x3 * (3 * x109 + 3 * x110 + 2 * x85 + 2 * x86)),
        ]
    )
    return S


def ovlp3d_34(a, A, b, B):
    """Cartesian 3D (fg) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = a * b * x0
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x7 = x5 * x6
    x8 = x3 * x7
    x9 = 3 * x8
    x10 = -x1 - B[0]
    x11 = x10**2 * x7
    x12 = x3 * (3 * x11 + x9)
    x13 = x10 * x8
    x14 = 2 * x13
    x15 = x11 + x8
    x16 = x10 * x15
    x17 = x14 + x16
    x18 = x10 * x17
    x19 = x17 * x2
    x20 = 8 * x13
    x21 = x3 * (4 * x16 + x20)
    x22 = x12 + x18
    x23 = x2 * x22
    x24 = x21 + x23
    x25 = x2 * x24 + x3 * (5 * x12 + x18 + 4 * x19)
    x26 = x15 * x2
    x27 = x3 * (x16 + x20 + 3 * x26)
    x28 = x12 + x19
    x29 = x2 * x28
    x30 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x31 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x32 = numpy.pi * x0 * x31
    x33 = x30 * x32
    x34 = -x0 * (a * A[1] + b * B[1])
    x35 = -x34 - B[1]
    x36 = x27 + x29
    x37 = x2 * x7
    x38 = x10 * x37
    x39 = 2 * x38 + x9
    x40 = x3 * (x11 + x39)
    x41 = x14 + x26
    x42 = x2 * x41
    x43 = x33 * (x2 * x36 + x3 * (2 * x12 + 2 * x19 + 3 * x40 + 3 * x42))
    x44 = -x0 * (a * A[2] + b * B[2])
    x45 = -x44 - B[2]
    x46 = x30 * x6
    x47 = x3 * x46
    x48 = x35**2 * x46
    x49 = x47 + x48
    x50 = x40 + x42
    x51 = x10 * x7
    x52 = x3 * (x37 + x51)
    x53 = x38 + x8
    x54 = x2 * x53
    x55 = x2 * x50 + x3 * (4 * x13 + 2 * x26 + 2 * x52 + 2 * x54)
    x56 = x31 * x6
    x57 = x33 * x45
    x58 = x3 * x56
    x59 = x45**2 * x56
    x60 = x58 + x59
    x61 = x35 * x47
    x62 = 2 * x61
    x63 = x35 * x49
    x64 = x62 + x63
    x65 = x2**2 * x7
    x66 = x52 + x54
    x67 = x2 * x66 + x3 * (x39 + x65)
    x68 = x45 * x56
    x69 = x35 * x46
    x70 = x45 * x58
    x71 = 2 * x70
    x72 = x45 * x60
    x73 = x71 + x72
    x74 = 3 * x47
    x75 = x3 * (3 * x48 + x74)
    x76 = x35 * x64
    x77 = x75 + x76
    x78 = x65 + x8
    x79 = x2 * x78 + 2 * x2 * x8
    x80 = 3 * x58
    x81 = x3 * (3 * x59 + x80)
    x82 = x45 * x73
    x83 = x81 + x82
    x84 = -x34 - A[1]
    x85 = x25 * x33
    x86 = x46 * x84
    x87 = x35 * x86
    x88 = x47 + x87
    x89 = x49 * x84
    x90 = x62 + x89
    x91 = x64 * x84
    x92 = x75 + x91
    x93 = 8 * x61
    x94 = x3 * (4 * x63 + x93)
    x95 = x77 * x84
    x96 = x94 + x95
    x97 = -x44 - A[2]
    x98 = x33 * x97
    x99 = x56 * x97
    x100 = x45 * x99
    x101 = x100 + x58
    x102 = x60 * x97
    x103 = x102 + x71
    x104 = x73 * x97
    x105 = x104 + x81
    x106 = 8 * x70
    x107 = x3 * (x106 + 4 * x72)
    x108 = x83 * x97
    x109 = x107 + x108
    x110 = x46 * x84**2
    x111 = x110 + x47
    x112 = x3 * (x69 + x86)
    x113 = x84 * x88
    x114 = x112 + x113
    x115 = x74 + 2 * x87
    x116 = x3 * (x115 + x48)
    x117 = x84 * x90
    x118 = x116 + x117
    x119 = x3 * (x63 + 3 * x89 + x93)
    x120 = x84 * x92
    x121 = x119 + x120
    x122 = x3 * (5 * x75 + x76 + 4 * x91) + x84 * x96
    x123 = x32 * x5
    x124 = x122 * x123
    x125 = x123 * x2
    x126 = numpy.pi * x0 * x30 * x5
    x127 = x126 * x2
    x128 = x56 * x97**2
    x129 = x128 + x58
    x130 = x3 * (x68 + x99)
    x131 = x101 * x97
    x132 = x130 + x131
    x133 = 2 * x100 + x80
    x134 = x3 * (x133 + x59)
    x135 = x103 * x97
    x136 = x134 + x135
    x137 = x3 * (3 * x102 + x106 + x72)
    x138 = x105 * x97
    x139 = x137 + x138
    x140 = x109 * x97 + x3 * (4 * x104 + 5 * x81 + x82)
    x141 = x126 * x140
    x142 = x111 * x84 + 2 * x47 * x84
    x143 = x114 * x84 + x3 * (x110 + x115)
    x144 = x118 * x84 + x3 * (2 * x112 + 2 * x113 + 4 * x61 + 2 * x89)
    x145 = x123 * (x121 * x84 + x3 * (3 * x116 + 3 * x117 + 2 * x75 + 2 * x91))
    x146 = x10 * x123
    x147 = x10 * x126
    x148 = x129 * x97 + 2 * x58 * x97
    x149 = x132 * x97 + x3 * (x128 + x133)
    x150 = x136 * x97 + x3 * (2 * x102 + 2 * x130 + 2 * x131 + 4 * x70)
    x151 = x126 * (x139 * x97 + x3 * (2 * x104 + 3 * x134 + 3 * x135 + 2 * x81))

    # 150 item(s)
    S = numpy.array(
        [
            x33 * (x2 * x25 + x3 * (2 * x21 + 2 * x23 + 4 * x27 + 4 * x29)),
            x35 * x43,
            x43 * x45,
            x49 * x55 * x56,
            x35 * x55 * x57,
            x46 * x55 * x60,
            x56 * x64 * x67,
            x49 * x67 * x68,
            x60 * x67 * x69,
            x46 * x67 * x73,
            x56 * x77 * x79,
            x64 * x68 * x79,
            x49 * x60 * x79,
            x69 * x73 * x79,
            x46 * x79 * x83,
            x84 * x85,
            x36 * x56 * x88,
            x36 * x57 * x84,
            x50 * x56 * x90,
            x50 * x68 * x88,
            x50 * x60 * x86,
            x56 * x66 * x92,
            x66 * x68 * x90,
            x60 * x66 * x88,
            x66 * x73 * x86,
            x56 * x78 * x96,
            x68 * x78 * x92,
            x60 * x78 * x90,
            x73 * x78 * x88,
            x78 * x83 * x86,
            x85 * x97,
            x35 * x36 * x98,
            x101 * x36 * x46,
            x49 * x50 * x99,
            x101 * x50 * x69,
            x103 * x46 * x50,
            x64 * x66 * x99,
            x101 * x49 * x66,
            x103 * x66 * x69,
            x105 * x46 * x66,
            x77 * x78 * x99,
            x101 * x64 * x78,
            x103 * x49 * x78,
            x105 * x69 * x78,
            x109 * x46 * x78,
            x111 * x24 * x56,
            x114 * x28 * x56,
            x111 * x28 * x68,
            x118 * x41 * x56,
            x114 * x41 * x68,
            x111 * x41 * x60,
            x121 * x53 * x56,
            x118 * x53 * x68,
            x114 * x53 * x60,
            x111 * x53 * x73,
            x124 * x2,
            x121 * x125 * x45,
            x118 * x37 * x60,
            x114 * x37 * x73,
            x111 * x37 * x83,
            x24 * x84 * x98,
            x28 * x88 * x99,
            x101 * x28 * x86,
            x41 * x90 * x99,
            x101 * x41 * x88,
            x103 * x41 * x86,
            x53 * x92 * x99,
            x101 * x53 * x90,
            x103 * x53 * x88,
            x105 * x53 * x86,
            x125 * x96 * x97,
            x101 * x37 * x92,
            x103 * x37 * x90,
            x105 * x37 * x88,
            x109 * x127 * x84,
            x129 * x24 * x46,
            x129 * x28 * x69,
            x132 * x28 * x46,
            x129 * x41 * x49,
            x132 * x41 * x69,
            x136 * x41 * x46,
            x129 * x53 * x64,
            x132 * x49 * x53,
            x136 * x53 * x69,
            x139 * x46 * x53,
            x129 * x37 * x77,
            x132 * x37 * x64,
            x136 * x37 * x49,
            x127 * x139 * x35,
            x141 * x2,
            x142 * x22 * x56,
            x143 * x17 * x56,
            x142 * x17 * x68,
            x144 * x15 * x56,
            x143 * x15 * x68,
            x142 * x15 * x60,
            x10 * x145,
            x144 * x146 * x45,
            x143 * x51 * x60,
            x142 * x51 * x73,
            x123 * (x122 * x84 + x3 * (4 * x119 + 4 * x120 + 2 * x94 + 2 * x95)),
            x145 * x45,
            x144 * x60 * x7,
            x143 * x7 * x73,
            x142 * x7 * x83,
            x111 * x22 * x99,
            x114 * x17 * x99,
            x101 * x111 * x17,
            x118 * x15 * x99,
            x101 * x114 * x15,
            x103 * x111 * x15,
            x121 * x146 * x97,
            x101 * x118 * x51,
            x103 * x114 * x51,
            x105 * x111 * x51,
            x124 * x97,
            x101 * x121 * x7,
            x103 * x118 * x7,
            x105 * x114 * x7,
            x109 * x111 * x7,
            x129 * x22 * x86,
            x129 * x17 * x88,
            x132 * x17 * x86,
            x129 * x15 * x90,
            x132 * x15 * x88,
            x136 * x15 * x86,
            x129 * x51 * x92,
            x132 * x51 * x90,
            x136 * x51 * x88,
            x139 * x147 * x84,
            x129 * x7 * x96,
            x132 * x7 * x92,
            x136 * x7 * x90,
            x139 * x7 * x88,
            x141 * x84,
            x148 * x22 * x46,
            x148 * x17 * x69,
            x149 * x17 * x46,
            x148 * x15 * x49,
            x149 * x15 * x69,
            x15 * x150 * x46,
            x148 * x51 * x64,
            x149 * x49 * x51,
            x147 * x150 * x35,
            x10 * x151,
            x148 * x7 * x77,
            x149 * x64 * x7,
            x150 * x49 * x7,
            x151 * x35,
            x126 * (x140 * x97 + x3 * (2 * x107 + 2 * x108 + 4 * x137 + 4 * x138)),
        ]
    )
    return S


def ovlp3d_40(a, A, b, B):
    """Cartesian 3D (gs) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = a * b * x1
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x5 = x3 * x4
    x6 = x0 * x5
    x7 = x1 * (a * A[0] + b * B[0]) - A[0]
    x8 = x5 * x7**2
    x9 = x6 + x8
    x10 = 2 * x6 * x7 + x7 * x9
    x11 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x12 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x13 = numpy.pi * x1 * x12
    x14 = x11 * x13
    x15 = x1 * (a * A[1] + b * B[1]) - A[1]
    x16 = x10 * x14
    x17 = x1 * (a * A[2] + b * B[2]) - A[2]
    x18 = x11 * x4
    x19 = x0 * x18
    x20 = x15**2 * x18
    x21 = x19 + x20
    x22 = x12 * x4
    x23 = x0 * x22
    x24 = x17**2 * x22
    x25 = x23 + x24
    x26 = 2 * x15 * x19 + x15 * x21
    x27 = x13 * x3
    x28 = x26 * x27
    x29 = numpy.pi * x1 * x11 * x3
    x30 = 2 * x17 * x23 + x17 * x25
    x31 = x29 * x30

    # 15 item(s)
    S = numpy.array(
        [
            x14 * (x0 * (3 * x6 + 3 * x8) + x10 * x7),
            x15 * x16,
            x16 * x17,
            x21 * x22 * x9,
            x14 * x15 * x17 * x9,
            x18 * x25 * x9,
            x28 * x7,
            x17 * x21 * x27 * x7,
            x15 * x25 * x29 * x7,
            x31 * x7,
            x27 * (x0 * (3 * x19 + 3 * x20) + x15 * x26),
            x17 * x28,
            x21 * x25 * x5,
            x15 * x31,
            x29 * (x0 * (3 * x23 + 3 * x24) + x17 * x30),
        ]
    )
    return S


def ovlp3d_41(a, A, b, B):
    """Cartesian 3D (gp) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - A[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3 * x7
    x9 = -x2 - B[0]
    x10 = x7 * x9
    x11 = x0 * (x10 + x8)
    x12 = x0 * x7
    x13 = x8 * x9
    x14 = x12 + x13
    x15 = x14 * x3
    x16 = x3**2 * x7
    x17 = x12 + x16
    x18 = 2 * x12 * x3 + x17 * x3
    x19 = 3 * x12
    x20 = x11 + x15
    x21 = x0 * (2 * x13 + x16 + x19) + x20 * x3
    x22 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x23 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x24 = numpy.pi * x1 * x23
    x25 = x22 * x24
    x26 = -x1 * (a * A[1] + b * B[1])
    x27 = -x26 - B[1]
    x28 = x25 * (x0 * (3 * x16 + x19) + x18 * x3)
    x29 = -x1 * (a * A[2] + b * B[2])
    x30 = -x29 - B[2]
    x31 = -x26 - A[1]
    x32 = x21 * x25
    x33 = x0 * x6
    x34 = x22 * x33
    x35 = x22 * x6
    x36 = x31 * x35
    x37 = x27 * x36
    x38 = x34 + x37
    x39 = x23 * x6
    x40 = x18 * x25
    x41 = -x29 - A[2]
    x42 = x23 * x33
    x43 = x39 * x41
    x44 = x30 * x43
    x45 = x42 + x44
    x46 = x31**2 * x35
    x47 = x34 + x46
    x48 = x27 * x35
    x49 = x0 * (x36 + x48)
    x50 = x31 * x38
    x51 = x49 + x50
    x52 = x30 * x39
    x53 = x39 * x41**2
    x54 = x42 + x53
    x55 = x0 * (x43 + x52)
    x56 = x41 * x45
    x57 = x55 + x56
    x58 = 2 * x31 * x34 + x31 * x47
    x59 = 3 * x34
    x60 = x0 * (2 * x37 + x46 + x59) + x31 * x51
    x61 = x24 * x5
    x62 = x60 * x61
    x63 = x3 * x61
    x64 = numpy.pi * x1 * x22 * x5
    x65 = x3 * x64
    x66 = 2 * x41 * x42 + x41 * x54
    x67 = 3 * x42
    x68 = x0 * (2 * x44 + x53 + x67) + x41 * x57
    x69 = x64 * x68
    x70 = x61 * (x0 * (3 * x46 + x59) + x31 * x58)
    x71 = x64 * (x0 * (3 * x53 + x67) + x41 * x66)

    # 45 item(s)
    S = numpy.array(
        [
            x25 * (x0 * (3 * x11 + 3 * x15 + x18) + x21 * x3),
            x27 * x28,
            x28 * x30,
            x31 * x32,
            x18 * x38 * x39,
            x30 * x31 * x40,
            x32 * x41,
            x27 * x40 * x41,
            x18 * x35 * x45,
            x20 * x39 * x47,
            x17 * x39 * x51,
            x17 * x47 * x52,
            x20 * x25 * x31 * x41,
            x17 * x38 * x43,
            x17 * x36 * x45,
            x20 * x35 * x54,
            x17 * x48 * x54,
            x17 * x35 * x57,
            x14 * x39 * x58,
            x3 * x62,
            x30 * x58 * x63,
            x14 * x43 * x47,
            x41 * x51 * x63,
            x45 * x47 * x8,
            x14 * x36 * x54,
            x38 * x54 * x8,
            x31 * x57 * x65,
            x14 * x35 * x66,
            x27 * x65 * x66,
            x3 * x69,
            x70 * x9,
            x61 * (x0 * (3 * x49 + 3 * x50 + x58) + x31 * x60),
            x30 * x70,
            x41 * x58 * x61 * x9,
            x41 * x62,
            x45 * x58 * x7,
            x10 * x47 * x54,
            x51 * x54 * x7,
            x47 * x57 * x7,
            x31 * x64 * x66 * x9,
            x38 * x66 * x7,
            x31 * x69,
            x71 * x9,
            x27 * x71,
            x64 * (x0 * (3 * x55 + 3 * x56 + x66) + x41 * x68),
        ]
    )
    return S


def ovlp3d_42(a, A, b, B):
    """Cartesian 3D (gd) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - A[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3**2 * x7
    x9 = x0 * x7
    x10 = 3 * x9
    x11 = -x2 - B[0]
    x12 = x3 * x7
    x13 = x11 * x12
    x14 = x10 + 2 * x13
    x15 = x0 * (x14 + x8)
    x16 = x11**2 * x7
    x17 = x0 * (x14 + x16)
    x18 = x11 * x7
    x19 = x0 * x18
    x20 = x16 + x9
    x21 = x20 * x3
    x22 = 2 * x19 + x21
    x23 = x22 * x3
    x24 = x0 * (x12 + x18)
    x25 = x13 + x9
    x26 = x25 * x3
    x27 = x24 + x26
    x28 = x27 * x3
    x29 = x17 + x23
    x30 = x0 * (4 * x19 + 2 * x21 + 2 * x24 + 2 * x26) + x29 * x3
    x31 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x32 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x33 = numpy.pi * x1 * x32
    x34 = x31 * x33
    x35 = -x1 * (a * A[1] + b * B[1])
    x36 = -x35 - B[1]
    x37 = x8 + x9
    x38 = 2 * x0 * x12 + x3 * x37
    x39 = x15 + x28
    x40 = x34 * (x0 * (3 * x24 + 3 * x26 + x38) + x3 * x39)
    x41 = -x1 * (a * A[2] + b * B[2])
    x42 = -x41 - B[2]
    x43 = x31 * x6
    x44 = x0 * x43
    x45 = x36**2 * x43
    x46 = x44 + x45
    x47 = x0 * (x10 + 3 * x8) + x3 * x38
    x48 = x32 * x6
    x49 = x34 * x42
    x50 = x0 * x48
    x51 = x42**2 * x48
    x52 = x50 + x51
    x53 = -x35 - A[1]
    x54 = x30 * x34
    x55 = x43 * x53
    x56 = x36 * x55
    x57 = x44 + x56
    x58 = 2 * x44
    x59 = x46 * x53
    x60 = x36 * x58 + x59
    x61 = x42 * x48
    x62 = -x41 - A[2]
    x63 = x34 * x62
    x64 = x48 * x62
    x65 = x42 * x64
    x66 = x50 + x65
    x67 = x36 * x43
    x68 = 2 * x50
    x69 = x52 * x62
    x70 = x42 * x68 + x69
    x71 = x43 * x53**2
    x72 = x44 + x71
    x73 = x0 * (x55 + x67)
    x74 = x53 * x57
    x75 = x73 + x74
    x76 = 3 * x44
    x77 = 2 * x56 + x76
    x78 = x0 * (x45 + x77)
    x79 = x53 * x60
    x80 = x78 + x79
    x81 = x48 * x62**2
    x82 = x50 + x81
    x83 = x0 * (x61 + x64)
    x84 = x62 * x66
    x85 = x83 + x84
    x86 = 3 * x50
    x87 = 2 * x65 + x86
    x88 = x0 * (x51 + x87)
    x89 = x62 * x70
    x90 = x88 + x89
    x91 = x53 * x58 + x53 * x72
    x92 = x0 * (x71 + x77)
    x93 = x53 * x75
    x94 = x92 + x93
    x95 = x0 * (4 * x36 * x44 + 2 * x59 + 2 * x73 + 2 * x74) + x53 * x80
    x96 = x33 * x5
    x97 = x95 * x96
    x98 = x3 * x96
    x99 = numpy.pi * x1 * x31 * x5
    x100 = x3 * x99
    x101 = x62 * x68 + x62 * x82
    x102 = x0 * (x81 + x87)
    x103 = x62 * x85
    x104 = x102 + x103
    x105 = x0 * (4 * x42 * x50 + 2 * x69 + 2 * x83 + 2 * x84) + x62 * x90
    x106 = x105 * x99
    x107 = x0 * (3 * x71 + x76) + x53 * x91
    x108 = x96 * (x0 * (3 * x73 + 3 * x74 + x91) + x53 * x94)
    x109 = x11 * x96
    x110 = x11 * x99
    x111 = x0 * (3 * x81 + x86) + x101 * x62
    x112 = x99 * (x0 * (x101 + 3 * x83 + 3 * x84) + x104 * x62)

    # 90 item(s)
    S = numpy.array(
        [
            x34 * (x0 * (2 * x15 + 3 * x17 + 3 * x23 + 2 * x28) + x3 * x30),
            x36 * x40,
            x40 * x42,
            x46 * x47 * x48,
            x36 * x47 * x49,
            x43 * x47 * x52,
            x53 * x54,
            x39 * x48 * x57,
            x39 * x49 * x53,
            x38 * x48 * x60,
            x38 * x57 * x61,
            x38 * x52 * x55,
            x54 * x62,
            x36 * x39 * x63,
            x39 * x43 * x66,
            x38 * x46 * x64,
            x38 * x66 * x67,
            x38 * x43 * x70,
            x29 * x48 * x72,
            x27 * x48 * x75,
            x27 * x61 * x72,
            x37 * x48 * x80,
            x37 * x61 * x75,
            x37 * x52 * x72,
            x29 * x53 * x63,
            x27 * x57 * x64,
            x27 * x55 * x66,
            x37 * x60 * x64,
            x37 * x57 * x66,
            x37 * x55 * x70,
            x29 * x43 * x82,
            x27 * x67 * x82,
            x27 * x43 * x85,
            x37 * x46 * x82,
            x37 * x67 * x85,
            x37 * x43 * x90,
            x22 * x48 * x91,
            x25 * x48 * x94,
            x25 * x61 * x91,
            x3 * x97,
            x42 * x94 * x98,
            x12 * x52 * x91,
            x22 * x64 * x72,
            x25 * x64 * x75,
            x25 * x66 * x72,
            x62 * x80 * x98,
            x12 * x66 * x75,
            x12 * x70 * x72,
            x22 * x55 * x82,
            x25 * x57 * x82,
            x25 * x55 * x85,
            x12 * x60 * x82,
            x12 * x57 * x85,
            x100 * x53 * x90,
            x101 * x22 * x43,
            x101 * x25 * x67,
            x104 * x25 * x43,
            x101 * x12 * x46,
            x100 * x104 * x36,
            x106 * x3,
            x107 * x20 * x48,
            x108 * x11,
            x107 * x109 * x42,
            x96 * (x0 * (3 * x78 + 3 * x79 + 2 * x92 + 2 * x93) + x53 * x95),
            x108 * x42,
            x107 * x52 * x7,
            x20 * x64 * x91,
            x109 * x62 * x94,
            x18 * x66 * x91,
            x62 * x97,
            x66 * x7 * x94,
            x7 * x70 * x91,
            x20 * x72 * x82,
            x18 * x75 * x82,
            x18 * x72 * x85,
            x7 * x80 * x82,
            x7 * x75 * x85,
            x7 * x72 * x90,
            x101 * x20 * x55,
            x101 * x18 * x57,
            x104 * x110 * x53,
            x101 * x60 * x7,
            x104 * x57 * x7,
            x106 * x53,
            x111 * x20 * x43,
            x110 * x111 * x36,
            x11 * x112,
            x111 * x46 * x7,
            x112 * x36,
            x99 * (x0 * (2 * x102 + 2 * x103 + 3 * x88 + 3 * x89) + x105 * x62),
        ]
    )
    return S


def ovlp3d_43(a, A, b, B):
    """Cartesian 3D (gf) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - B[0]
    x5 = a * b * x0
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x8 = x6 * x7
    x9 = x3 * x8
    x10 = x4 * x9
    x11 = x4**2 * x8
    x12 = x11 + x9
    x13 = x12 * x4
    x14 = x12 * x2
    x15 = x3 * (8 * x10 + x13 + 3 * x14)
    x16 = 3 * x9
    x17 = x3 * (3 * x11 + x16)
    x18 = 2 * x10
    x19 = x13 + x18
    x20 = x19 * x2
    x21 = x17 + x20
    x22 = x2 * x21
    x23 = x15 + x22
    x24 = x2 * x8
    x25 = x24 * x4
    x26 = x16 + 2 * x25
    x27 = x3 * (x11 + x26)
    x28 = x14 + x18
    x29 = x2 * x28
    x30 = 3 * x27 + 3 * x29
    x31 = x2 * x23 + x3 * (2 * x17 + 2 * x20 + x30)
    x32 = x27 + x29
    x33 = x2 * x32
    x34 = x4 * x8
    x35 = x3 * (x24 + x34)
    x36 = x25 + x9
    x37 = x2 * x36
    x38 = x3 * (4 * x10 + 2 * x14 + 2 * x35 + 2 * x37)
    x39 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x40 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x41 = numpy.pi * x0 * x40
    x42 = x39 * x41
    x43 = -x0 * (a * A[1] + b * B[1])
    x44 = -x43 - B[1]
    x45 = x2**2 * x8
    x46 = x3 * (x26 + x45)
    x47 = x35 + x37
    x48 = x2 * x47
    x49 = x33 + x38
    x50 = x42 * (x2 * x49 + x3 * (x30 + 2 * x46 + 2 * x48))
    x51 = -x0 * (a * A[2] + b * B[2])
    x52 = -x51 - B[2]
    x53 = x39 * x7
    x54 = x3 * x53
    x55 = x44**2 * x53
    x56 = x54 + x55
    x57 = x45 + x9
    x58 = x2 * x57 + 2 * x2 * x9
    x59 = x46 + x48
    x60 = x2 * x59 + x3 * (3 * x35 + 3 * x37 + x58)
    x61 = x40 * x7
    x62 = x42 * x52
    x63 = x3 * x61
    x64 = x52**2 * x61
    x65 = x63 + x64
    x66 = x44 * x54
    x67 = 2 * x66
    x68 = x44 * x56
    x69 = x67 + x68
    x70 = x2 * x58 + x3 * (x16 + 3 * x45)
    x71 = x52 * x61
    x72 = x44 * x53
    x73 = x52 * x63
    x74 = 2 * x73
    x75 = x52 * x65
    x76 = x74 + x75
    x77 = -x43 - A[1]
    x78 = x31 * x42
    x79 = x53 * x77
    x80 = x44 * x79
    x81 = x54 + x80
    x82 = x56 * x77
    x83 = x67 + x82
    x84 = 3 * x54
    x85 = x3 * (3 * x55 + x84)
    x86 = x69 * x77
    x87 = x85 + x86
    x88 = -x51 - A[2]
    x89 = x42 * x88
    x90 = x61 * x88
    x91 = x52 * x90
    x92 = x63 + x91
    x93 = x65 * x88
    x94 = x74 + x93
    x95 = 3 * x63
    x96 = x3 * (3 * x64 + x95)
    x97 = x76 * x88
    x98 = x96 + x97
    x99 = x53 * x77**2
    x100 = x54 + x99
    x101 = x3 * (x72 + x79)
    x102 = x77 * x81
    x103 = x101 + x102
    x104 = 2 * x80 + x84
    x105 = x3 * (x104 + x55)
    x106 = x77 * x83
    x107 = x105 + x106
    x108 = x3 * (8 * x66 + x68 + 3 * x82)
    x109 = x77 * x87
    x110 = x108 + x109
    x111 = x61 * x88**2
    x112 = x111 + x63
    x113 = x3 * (x71 + x90)
    x114 = x88 * x92
    x115 = x113 + x114
    x116 = 2 * x91 + x95
    x117 = x3 * (x116 + x64)
    x118 = x88 * x94
    x119 = x117 + x118
    x120 = x3 * (8 * x73 + x75 + 3 * x93)
    x121 = x88 * x98
    x122 = x120 + x121
    x123 = x100 * x77 + 2 * x54 * x77
    x124 = x3 * (x104 + x99)
    x125 = x103 * x77
    x126 = x124 + x125
    x127 = x107 * x77
    x128 = x3 * (2 * x101 + 2 * x102 + 4 * x66 + 2 * x82)
    x129 = x127 + x128
    x130 = 3 * x105 + 3 * x106
    x131 = x110 * x77 + x3 * (x130 + 2 * x85 + 2 * x86)
    x132 = x41 * x6
    x133 = x131 * x132
    x134 = x132 * x2
    x135 = numpy.pi * x0 * x39 * x6
    x136 = x135 * x2
    x137 = x112 * x88 + 2 * x63 * x88
    x138 = x3 * (x111 + x116)
    x139 = x115 * x88
    x140 = x138 + x139
    x141 = x119 * x88
    x142 = x3 * (2 * x113 + 2 * x114 + 4 * x73 + 2 * x93)
    x143 = x141 + x142
    x144 = 3 * x117 + 3 * x118
    x145 = x122 * x88 + x3 * (x144 + 2 * x96 + 2 * x97)
    x146 = x135 * x145
    x147 = x123 * x77 + x3 * (x84 + 3 * x99)
    x148 = x126 * x77 + x3 * (3 * x101 + 3 * x102 + x123)
    x149 = x132 * (x129 * x77 + x3 * (2 * x124 + 2 * x125 + x130))
    x150 = x132 * x4
    x151 = x135 * x4
    x152 = x137 * x88 + x3 * (3 * x111 + x95)
    x153 = x140 * x88 + x3 * (3 * x113 + 3 * x114 + x137)
    x154 = x135 * (x143 * x88 + x3 * (2 * x138 + 2 * x139 + x144))

    # 150 item(s)
    S = numpy.array(
        [
            x42 * (x2 * x31 + x3 * (3 * x15 + 3 * x22 + 3 * x33 + 3 * x38)),
            x44 * x50,
            x50 * x52,
            x56 * x60 * x61,
            x44 * x60 * x62,
            x53 * x60 * x65,
            x61 * x69 * x70,
            x56 * x70 * x71,
            x65 * x70 * x72,
            x53 * x70 * x76,
            x77 * x78,
            x49 * x61 * x81,
            x49 * x62 * x77,
            x59 * x61 * x83,
            x59 * x71 * x81,
            x59 * x65 * x79,
            x58 * x61 * x87,
            x58 * x71 * x83,
            x58 * x65 * x81,
            x58 * x76 * x79,
            x78 * x88,
            x44 * x49 * x89,
            x49 * x53 * x92,
            x56 * x59 * x90,
            x59 * x72 * x92,
            x53 * x59 * x94,
            x58 * x69 * x90,
            x56 * x58 * x92,
            x58 * x72 * x94,
            x53 * x58 * x98,
            x100 * x23 * x61,
            x103 * x32 * x61,
            x100 * x32 * x71,
            x107 * x47 * x61,
            x103 * x47 * x71,
            x100 * x47 * x65,
            x110 * x57 * x61,
            x107 * x57 * x71,
            x103 * x57 * x65,
            x100 * x57 * x76,
            x23 * x77 * x89,
            x32 * x81 * x90,
            x32 * x79 * x92,
            x47 * x83 * x90,
            x47 * x81 * x92,
            x47 * x79 * x94,
            x57 * x87 * x90,
            x57 * x83 * x92,
            x57 * x81 * x94,
            x57 * x79 * x98,
            x112 * x23 * x53,
            x112 * x32 * x72,
            x115 * x32 * x53,
            x112 * x47 * x56,
            x115 * x47 * x72,
            x119 * x47 * x53,
            x112 * x57 * x69,
            x115 * x56 * x57,
            x119 * x57 * x72,
            x122 * x53 * x57,
            x123 * x21 * x61,
            x126 * x28 * x61,
            x123 * x28 * x71,
            x129 * x36 * x61,
            x126 * x36 * x71,
            x123 * x36 * x65,
            x133 * x2,
            x129 * x134 * x52,
            x126 * x24 * x65,
            x123 * x24 * x76,
            x100 * x21 * x90,
            x103 * x28 * x90,
            x100 * x28 * x92,
            x107 * x36 * x90,
            x103 * x36 * x92,
            x100 * x36 * x94,
            x110 * x134 * x88,
            x107 * x24 * x92,
            x103 * x24 * x94,
            x100 * x24 * x98,
            x112 * x21 * x79,
            x112 * x28 * x81,
            x115 * x28 * x79,
            x112 * x36 * x83,
            x115 * x36 * x81,
            x119 * x36 * x79,
            x112 * x24 * x87,
            x115 * x24 * x83,
            x119 * x24 * x81,
            x122 * x136 * x77,
            x137 * x21 * x53,
            x137 * x28 * x72,
            x140 * x28 * x53,
            x137 * x36 * x56,
            x140 * x36 * x72,
            x143 * x36 * x53,
            x137 * x24 * x69,
            x140 * x24 * x56,
            x136 * x143 * x44,
            x146 * x2,
            x147 * x19 * x61,
            x12 * x148 * x61,
            x12 * x147 * x71,
            x149 * x4,
            x148 * x150 * x52,
            x147 * x34 * x65,
            x132 * (x131 * x77 + x3 * (3 * x108 + 3 * x109 + 3 * x127 + 3 * x128)),
            x149 * x52,
            x148 * x65 * x8,
            x147 * x76 * x8,
            x123 * x19 * x90,
            x12 * x126 * x90,
            x12 * x123 * x92,
            x129 * x150 * x88,
            x126 * x34 * x92,
            x123 * x34 * x94,
            x133 * x88,
            x129 * x8 * x92,
            x126 * x8 * x94,
            x123 * x8 * x98,
            x100 * x112 * x19,
            x103 * x112 * x12,
            x100 * x115 * x12,
            x107 * x112 * x34,
            x103 * x115 * x34,
            x100 * x119 * x34,
            x110 * x112 * x8,
            x107 * x115 * x8,
            x103 * x119 * x8,
            x100 * x122 * x8,
            x137 * x19 * x79,
            x12 * x137 * x81,
            x12 * x140 * x79,
            x137 * x34 * x83,
            x140 * x34 * x81,
            x143 * x151 * x77,
            x137 * x8 * x87,
            x140 * x8 * x83,
            x143 * x8 * x81,
            x146 * x77,
            x152 * x19 * x53,
            x12 * x152 * x72,
            x12 * x153 * x53,
            x152 * x34 * x56,
            x151 * x153 * x44,
            x154 * x4,
            x152 * x69 * x8,
            x153 * x56 * x8,
            x154 * x44,
            x135 * (x145 * x88 + x3 * (3 * x120 + 3 * x121 + 3 * x141 + 3 * x142)),
        ]
    )
    return S


def ovlp3d_44(a, A, b, B):
    """Cartesian 3D (gg) overlap integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = a * b * x0
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x7 = x5 * x6
    x8 = x3 * x7
    x9 = 3 * x8
    x10 = -x1 - B[0]
    x11 = x10**2 * x7
    x12 = x3 * (3 * x11 + x9)
    x13 = x10 * x8
    x14 = 2 * x13
    x15 = x11 + x8
    x16 = x10 * x15
    x17 = x14 + x16
    x18 = x10 * x17
    x19 = x17 * x2
    x20 = x3 * (5 * x12 + x18 + 4 * x19)
    x21 = 8 * x13
    x22 = x3 * (4 * x16 + x21)
    x23 = x12 + x18
    x24 = x2 * x23
    x25 = x22 + x24
    x26 = x2 * x25
    x27 = x20 + x26
    x28 = x15 * x2
    x29 = x3 * (x16 + x21 + 3 * x28)
    x30 = x12 + x19
    x31 = x2 * x30
    x32 = x2 * x27 + x3 * (2 * x22 + 2 * x24 + 4 * x29 + 4 * x31)
    x33 = x29 + x31
    x34 = x2 * x33
    x35 = x2 * x7
    x36 = x10 * x35
    x37 = 2 * x36 + x9
    x38 = x3 * (x11 + x37)
    x39 = x14 + x28
    x40 = x2 * x39
    x41 = 3 * x38 + 3 * x40
    x42 = x3 * (2 * x12 + 2 * x19 + x41)
    x43 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x44 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x45 = numpy.pi * x0 * x44
    x46 = x43 * x45
    x47 = -x0 * (a * A[1] + b * B[1])
    x48 = -x47 - B[1]
    x49 = x34 + x42
    x50 = x38 + x40
    x51 = x2 * x50
    x52 = x10 * x7
    x53 = x3 * (x35 + x52)
    x54 = x36 + x8
    x55 = x2 * x54
    x56 = x3 * (4 * x13 + 2 * x28 + 2 * x53 + 2 * x55)
    x57 = x46 * (x2 * x49 + x3 * (3 * x29 + 3 * x31 + 3 * x51 + 3 * x56))
    x58 = -x0 * (a * A[2] + b * B[2])
    x59 = -x58 - B[2]
    x60 = x43 * x6
    x61 = x3 * x60
    x62 = x48**2 * x60
    x63 = x61 + x62
    x64 = x2**2 * x7
    x65 = x3 * (x37 + x64)
    x66 = x53 + x55
    x67 = x2 * x66
    x68 = x51 + x56
    x69 = x2 * x68 + x3 * (x41 + 2 * x65 + 2 * x67)
    x70 = x44 * x6
    x71 = x46 * x59
    x72 = x3 * x70
    x73 = x59**2 * x70
    x74 = x72 + x73
    x75 = x48 * x61
    x76 = 2 * x75
    x77 = x48 * x63
    x78 = x76 + x77
    x79 = x64 + x8
    x80 = x2 * x79 + 2 * x2 * x8
    x81 = x65 + x67
    x82 = x2 * x81 + x3 * (3 * x53 + 3 * x55 + x80)
    x83 = x59 * x70
    x84 = x48 * x60
    x85 = x59 * x72
    x86 = 2 * x85
    x87 = x59 * x74
    x88 = x86 + x87
    x89 = 3 * x61
    x90 = x3 * (3 * x62 + x89)
    x91 = x48 * x78
    x92 = x90 + x91
    x93 = x2 * x80 + x3 * (3 * x64 + x9)
    x94 = 3 * x72
    x95 = x3 * (3 * x73 + x94)
    x96 = x59 * x88
    x97 = x95 + x96
    x98 = -x47 - A[1]
    x99 = x32 * x46
    x100 = x60 * x98
    x101 = x100 * x48
    x102 = x101 + x61
    x103 = x63 * x98
    x104 = x103 + x76
    x105 = x78 * x98
    x106 = x105 + x90
    x107 = 8 * x75
    x108 = x3 * (x107 + 4 * x77)
    x109 = x92 * x98
    x110 = x108 + x109
    x111 = -x58 - A[2]
    x112 = x111 * x46
    x113 = x111 * x70
    x114 = x113 * x59
    x115 = x114 + x72
    x116 = x111 * x74
    x117 = x116 + x86
    x118 = x111 * x88
    x119 = x118 + x95
    x120 = 8 * x85
    x121 = x3 * (x120 + 4 * x87)
    x122 = x111 * x97
    x123 = x121 + x122
    x124 = x60 * x98**2
    x125 = x124 + x61
    x126 = x3 * (x100 + x84)
    x127 = x102 * x98
    x128 = x126 + x127
    x129 = 2 * x101 + x89
    x130 = x3 * (x129 + x62)
    x131 = x104 * x98
    x132 = x130 + x131
    x133 = x3 * (3 * x103 + x107 + x77)
    x134 = x106 * x98
    x135 = x133 + x134
    x136 = x3 * (4 * x105 + 5 * x90 + x91)
    x137 = x110 * x98
    x138 = x136 + x137
    x139 = x111**2 * x70
    x140 = x139 + x72
    x141 = x3 * (x113 + x83)
    x142 = x111 * x115
    x143 = x141 + x142
    x144 = 2 * x114 + x94
    x145 = x3 * (x144 + x73)
    x146 = x111 * x117
    x147 = x145 + x146
    x148 = x3 * (3 * x116 + x120 + x87)
    x149 = x111 * x119
    x150 = x148 + x149
    x151 = x3 * (4 * x118 + 5 * x95 + x96)
    x152 = x111 * x123
    x153 = x151 + x152
    x154 = x125 * x98 + 2 * x61 * x98
    x155 = x3 * (x124 + x129)
    x156 = x128 * x98
    x157 = x155 + x156
    x158 = x132 * x98
    x159 = x3 * (2 * x103 + 2 * x126 + 2 * x127 + 4 * x75)
    x160 = x158 + x159
    x161 = x135 * x98
    x162 = 3 * x130 + 3 * x131
    x163 = x3 * (2 * x105 + x162 + 2 * x90)
    x164 = x161 + x163
    x165 = x138 * x98 + x3 * (2 * x108 + 2 * x109 + 4 * x133 + 4 * x134)
    x166 = x45 * x5
    x167 = x165 * x166
    x168 = x166 * x2
    x169 = numpy.pi * x0 * x43 * x5
    x170 = x169 * x2
    x171 = x111 * x140 + 2 * x111 * x72
    x172 = x3 * (x139 + x144)
    x173 = x111 * x143
    x174 = x172 + x173
    x175 = x111 * x147
    x176 = x3 * (2 * x116 + 2 * x141 + 2 * x142 + 4 * x85)
    x177 = x175 + x176
    x178 = x111 * x150
    x179 = 3 * x145 + 3 * x146
    x180 = x3 * (2 * x118 + x179 + 2 * x95)
    x181 = x178 + x180
    x182 = x111 * x153 + x3 * (2 * x121 + 2 * x122 + 4 * x148 + 4 * x149)
    x183 = x169 * x182
    x184 = x154 * x98 + x3 * (3 * x124 + x89)
    x185 = x157 * x98 + x3 * (3 * x126 + 3 * x127 + x154)
    x186 = x160 * x98 + x3 * (2 * x155 + 2 * x156 + x162)
    x187 = x166 * (x164 * x98 + x3 * (3 * x133 + 3 * x134 + 3 * x158 + 3 * x159))
    x188 = x10 * x166
    x189 = x10 * x169
    x190 = x111 * x171 + x3 * (3 * x139 + x94)
    x191 = x111 * x174 + x3 * (3 * x141 + 3 * x142 + x171)
    x192 = x111 * x177 + x3 * (2 * x172 + 2 * x173 + x179)
    x193 = x169 * (x111 * x181 + x3 * (3 * x148 + 3 * x149 + 3 * x175 + 3 * x176))

    # 225 item(s)
    S = numpy.array(
        [
            x46 * (x2 * x32 + x3 * (3 * x20 + 3 * x26 + 4 * x34 + 4 * x42)),
            x48 * x57,
            x57 * x59,
            x63 * x69 * x70,
            x48 * x69 * x71,
            x60 * x69 * x74,
            x70 * x78 * x82,
            x63 * x82 * x83,
            x74 * x82 * x84,
            x60 * x82 * x88,
            x70 * x92 * x93,
            x78 * x83 * x93,
            x63 * x74 * x93,
            x84 * x88 * x93,
            x60 * x93 * x97,
            x98 * x99,
            x102 * x49 * x70,
            x49 * x71 * x98,
            x104 * x68 * x70,
            x102 * x68 * x83,
            x100 * x68 * x74,
            x106 * x70 * x81,
            x104 * x81 * x83,
            x102 * x74 * x81,
            x100 * x81 * x88,
            x110 * x70 * x80,
            x106 * x80 * x83,
            x104 * x74 * x80,
            x102 * x80 * x88,
            x100 * x80 * x97,
            x111 * x99,
            x112 * x48 * x49,
            x115 * x49 * x60,
            x113 * x63 * x68,
            x115 * x68 * x84,
            x117 * x60 * x68,
            x113 * x78 * x81,
            x115 * x63 * x81,
            x117 * x81 * x84,
            x119 * x60 * x81,
            x113 * x80 * x92,
            x115 * x78 * x80,
            x117 * x63 * x80,
            x119 * x80 * x84,
            x123 * x60 * x80,
            x125 * x27 * x70,
            x128 * x33 * x70,
            x125 * x33 * x83,
            x132 * x50 * x70,
            x128 * x50 * x83,
            x125 * x50 * x74,
            x135 * x66 * x70,
            x132 * x66 * x83,
            x128 * x66 * x74,
            x125 * x66 * x88,
            x138 * x70 * x79,
            x135 * x79 * x83,
            x132 * x74 * x79,
            x128 * x79 * x88,
            x125 * x79 * x97,
            x112 * x27 * x98,
            x102 * x113 * x33,
            x100 * x115 * x33,
            x104 * x113 * x50,
            x102 * x115 * x50,
            x100 * x117 * x50,
            x106 * x113 * x66,
            x104 * x115 * x66,
            x102 * x117 * x66,
            x100 * x119 * x66,
            x110 * x113 * x79,
            x106 * x115 * x79,
            x104 * x117 * x79,
            x102 * x119 * x79,
            x100 * x123 * x79,
            x140 * x27 * x60,
            x140 * x33 * x84,
            x143 * x33 * x60,
            x140 * x50 * x63,
            x143 * x50 * x84,
            x147 * x50 * x60,
            x140 * x66 * x78,
            x143 * x63 * x66,
            x147 * x66 * x84,
            x150 * x60 * x66,
            x140 * x79 * x92,
            x143 * x78 * x79,
            x147 * x63 * x79,
            x150 * x79 * x84,
            x153 * x60 * x79,
            x154 * x25 * x70,
            x157 * x30 * x70,
            x154 * x30 * x83,
            x160 * x39 * x70,
            x157 * x39 * x83,
            x154 * x39 * x74,
            x164 * x54 * x70,
            x160 * x54 * x83,
            x157 * x54 * x74,
            x154 * x54 * x88,
            x167 * x2,
            x164 * x168 * x59,
            x160 * x35 * x74,
            x157 * x35 * x88,
            x154 * x35 * x97,
            x113 * x125 * x25,
            x113 * x128 * x30,
            x115 * x125 * x30,
            x113 * x132 * x39,
            x115 * x128 * x39,
            x117 * x125 * x39,
            x113 * x135 * x54,
            x115 * x132 * x54,
            x117 * x128 * x54,
            x119 * x125 * x54,
            x111 * x138 * x168,
            x115 * x135 * x35,
            x117 * x132 * x35,
            x119 * x128 * x35,
            x123 * x125 * x35,
            x100 * x140 * x25,
            x102 * x140 * x30,
            x100 * x143 * x30,
            x104 * x140 * x39,
            x102 * x143 * x39,
            x100 * x147 * x39,
            x106 * x140 * x54,
            x104 * x143 * x54,
            x102 * x147 * x54,
            x100 * x150 * x54,
            x110 * x140 * x35,
            x106 * x143 * x35,
            x104 * x147 * x35,
            x102 * x150 * x35,
            x153 * x170 * x98,
            x171 * x25 * x60,
            x171 * x30 * x84,
            x174 * x30 * x60,
            x171 * x39 * x63,
            x174 * x39 * x84,
            x177 * x39 * x60,
            x171 * x54 * x78,
            x174 * x54 * x63,
            x177 * x54 * x84,
            x181 * x54 * x60,
            x171 * x35 * x92,
            x174 * x35 * x78,
            x177 * x35 * x63,
            x170 * x181 * x48,
            x183 * x2,
            x184 * x23 * x70,
            x17 * x185 * x70,
            x17 * x184 * x83,
            x15 * x186 * x70,
            x15 * x185 * x83,
            x15 * x184 * x74,
            x10 * x187,
            x186 * x188 * x59,
            x185 * x52 * x74,
            x184 * x52 * x88,
            x166 * (x165 * x98 + x3 * (3 * x136 + 3 * x137 + 4 * x161 + 4 * x163)),
            x187 * x59,
            x186 * x7 * x74,
            x185 * x7 * x88,
            x184 * x7 * x97,
            x113 * x154 * x23,
            x113 * x157 * x17,
            x115 * x154 * x17,
            x113 * x15 * x160,
            x115 * x15 * x157,
            x117 * x15 * x154,
            x111 * x164 * x188,
            x115 * x160 * x52,
            x117 * x157 * x52,
            x119 * x154 * x52,
            x111 * x167,
            x115 * x164 * x7,
            x117 * x160 * x7,
            x119 * x157 * x7,
            x123 * x154 * x7,
            x125 * x140 * x23,
            x128 * x140 * x17,
            x125 * x143 * x17,
            x132 * x140 * x15,
            x128 * x143 * x15,
            x125 * x147 * x15,
            x135 * x140 * x52,
            x132 * x143 * x52,
            x128 * x147 * x52,
            x125 * x150 * x52,
            x138 * x140 * x7,
            x135 * x143 * x7,
            x132 * x147 * x7,
            x128 * x150 * x7,
            x125 * x153 * x7,
            x100 * x171 * x23,
            x102 * x17 * x171,
            x100 * x17 * x174,
            x104 * x15 * x171,
            x102 * x15 * x174,
            x100 * x15 * x177,
            x106 * x171 * x52,
            x104 * x174 * x52,
            x102 * x177 * x52,
            x181 * x189 * x98,
            x110 * x171 * x7,
            x106 * x174 * x7,
            x104 * x177 * x7,
            x102 * x181 * x7,
            x183 * x98,
            x190 * x23 * x60,
            x17 * x190 * x84,
            x17 * x191 * x60,
            x15 * x190 * x63,
            x15 * x191 * x84,
            x15 * x192 * x60,
            x190 * x52 * x78,
            x191 * x52 * x63,
            x189 * x192 * x48,
            x10 * x193,
            x190 * x7 * x92,
            x191 * x7 * x78,
            x192 * x63 * x7,
            x193 * x48,
            x169 * (x111 * x182 + x3 * (3 * x151 + 3 * x152 + 4 * x178 + 4 * x180)),
        ]
    )
    return S
