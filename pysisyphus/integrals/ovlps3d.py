import numpy


def ovlp3d_00(a, A, b, B):
    """Cartesian 3d (ss) overlap.

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
    """Cartesian 3d (sp) overlap.

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
    """Cartesian 3d (sd) overlap.

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
            x10 * x2 * (x3 * x7 + x7 * x8 ** 2),
            x11 * x14,
            x14 * x15,
            x10 * x4 * (x11 ** 2 * x16 + x16 * x3),
            x11 * x12 * x13 * x15 * x5 * x9,
            numpy.pi * x13 * (x15 ** 2 * x17 + x17 * x3),
        ]
    )
    return S


def ovlp3d_10(a, A, b, B):
    """Cartesian 3d (ps) overlap.

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
    """Cartesian 3d (pp) overlap.

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
    """Cartesian 3d (pd) overlap.

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
    x11 = x2 ** 2 * x8 + x9
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
    x23 = x17 ** 2 * x21 + x22
    x24 = x14 * x5
    x25 = x23 * x24
    x26 = numpy.pi ** (3 / 2)
    x27 = x0 * x12 * x5
    x28 = x13 * x20 * x26 * x27 * x6
    x29 = x13 * x7
    x30 = x29 * x3
    x31 = x20 ** 2 * x29 + x30
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


def ovlp3d_20(a, A, b, B):
    """Cartesian 3d (ds) overlap.

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
            x10 * x2 * (x3 * x7 + x7 * x8 ** 2),
            x11 * x14,
            x14 * x15,
            x10 * x4 * (x11 ** 2 * x16 + x16 * x3),
            x11 * x12 * x13 * x15 * x5 * x9,
            numpy.pi * x13 * (x15 ** 2 * x17 + x17 * x3),
        ]
    )
    return S


def ovlp3d_21(a, A, b, B):
    """Cartesian 3d (dp) overlap.

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
    x20 = x17 * (x11 * x5 ** 2 + x12)
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
    x40 = x39 * (x23 ** 2 * x27 + x26)
    x41 = x38 * x7
    x42 = x41 * (x33 ** 2 * x35 + x34)

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
    """Cartesian 3d (dd) overlap.

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
    x9 = x5 * x8 ** 2
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
    x27 = x20 ** 2 * x25
    x28 = x26 + x27
    x29 = x10 ** 2 * x5 + x6
    x30 = x16 * x4
    x31 = x18 * x24
    x32 = x0 * x30
    x33 = x24 ** 2 * x30
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
    x53 = x25 * x35 ** 2 + x26
    x54 = x41 * (x0 * (x20 * x25 + x37) + x35 * x39)
    x55 = x41 * x8
    x56 = x44 * x8
    x57 = x30 * x46 ** 2 + x32
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
