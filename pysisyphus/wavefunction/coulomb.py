import numpy

from pysisyphus.wavefunction.boys import neville_boys as boys


def coulomb_00(a, A, b, B, C):
    """Cartesian (ss) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = a + b
    x1 = x0 ** (-1.0)
    x2 = (
        (x1 * (a * A[0] + b * B[0]) - C[0]) ** 2
        + (x1 * (a * A[1] + b * B[1]) - C[1]) ** 2
        + (x1 * (a * A[2] + b * B[2]) - C[2]) ** 2
    )

    # 1 item(s)
    S = numpy.array([2 * numpy.pi * x1 * boys(0, x0 * x2) * numpy.exp(-a * b * x1 * x2)])
    return S


def coulomb_01(a, A, b, B, C):
    """Cartesian (sp) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = a + b
    x1 = x0 ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x1 * (a * A[1] + b * B[1])
    x4 = -x1 * (a * A[2] + b * B[2])
    x5 = x2 + C[0]
    x6 = x3 + C[1]
    x7 = x4 + C[2]
    x8 = x5 ** 2 + x6 ** 2 + x7 ** 2
    x9 = x0 * x8
    x10 = 2 * numpy.pi * x1 * numpy.exp(-a * b * x1 * x8)
    x11 = x10 * boys(0, x9)
    x12 = -x10 * numpy.sqrt(abs(x5) ** 2 + abs(x6) ** 2 + abs(x7) ** 2) * boys(1, x9)
    x13 = (
        x11 * numpy.sqrt(abs(x2 + B[0]) ** 2 + abs(x3 + B[1]) ** 2 + abs(x4 + B[2]) ** 2)
        + x12
    )

    # 3 item(s)
    S = numpy.array(
        [
            x13,
            x11
            * numpy.sqrt(abs(x2 + A[0]) ** 2 + abs(x3 + A[1]) ** 2 + abs(x4 + A[2]) ** 2)
            + x12,
            x13,
        ]
    )
    return S


def coulomb_02(a, A, b, B, C):
    """Cartesian (sd) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = a + b
    x1 = x0 ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x1 * (a * A[1] + b * B[1])
    x4 = -x1 * (a * A[2] + b * B[2])
    x5 = numpy.sqrt(abs(x2 + B[0]) ** 2 + abs(x3 + B[1]) ** 2 + abs(x4 + B[2]) ** 2)
    x6 = x2 + C[0]
    x7 = x3 + C[1]
    x8 = x4 + C[2]
    x9 = x6 ** 2 + x7 ** 2 + x8 ** 2
    x10 = x0 * x9
    x11 = boys(0, x10)
    x12 = 2 * numpy.pi * x1 * numpy.exp(-a * b * x1 * x9)
    x13 = x12 * x5
    x14 = boys(1, x10)
    x15 = numpy.sqrt(abs(x6) ** 2 + abs(x7) ** 2 + abs(x8) ** 2)
    x16 = x12 * x15
    x17 = -x14 * x16
    x18 = x11 * x13 + x17
    x19 = -x16 * boys(2, x10)
    x20 = -x15 * (x13 * x14 + x19)
    x21 = x18 * x5 + x20
    x22 = numpy.sqrt(abs(x2 + A[0]) ** 2 + abs(x3 + A[1]) ** 2 + abs(x4 + A[2]) ** 2)
    x23 = x12 * x22
    x24 = x11 * x23 + x17
    x25 = x14 * x23 + x19
    x26 = -x15 * x25
    x27 = (2 * a + 2 * b) ** (-1.0)
    x28 = x12 * x27

    # 6 item(s)
    S = numpy.array(
        [
            x21,
            x24 * x5 + x26,
            x21,
            x11 * x28 - x14 * x28 + x22 * x24 + x26,
            x18 * x22 + x20 + x24 * x27 - x25 * x27,
            x21,
        ]
    )
    return S


def coulomb_03(a, A, b, B, C):
    """Cartesian (sf) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = a + b
    x1 = x0 ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x1 * (a * A[1] + b * B[1])
    x4 = -x1 * (a * A[2] + b * B[2])
    x5 = numpy.sqrt(abs(x2 + B[0]) ** 2 + abs(x3 + B[1]) ** 2 + abs(x4 + B[2]) ** 2)
    x6 = x2 + C[0]
    x7 = x3 + C[1]
    x8 = x4 + C[2]
    x9 = x6 ** 2 + x7 ** 2 + x8 ** 2
    x10 = x0 * x9
    x11 = boys(0, x10)
    x12 = numpy.pi * x1 * numpy.exp(-a * b * x1 * x9)
    x13 = 2 * x12
    x14 = x13 * x5
    x15 = boys(1, x10)
    x16 = numpy.sqrt(abs(x6) ** 2 + abs(x7) ** 2 + abs(x8) ** 2)
    x17 = x13 * x16
    x18 = -x15 * x17
    x19 = x11 * x14 + x18
    x20 = boys(2, x10)
    x21 = -x17 * x20
    x22 = x14 * x15 + x21
    x23 = x16 * x22
    x24 = -x23
    x25 = x19 * x5 + x24
    x26 = -x17 * boys(3, x10)
    x27 = x16 * (x14 * x20 + x26)
    x28 = -x27
    x29 = -x16 * (x22 * x5 + x28)
    x30 = x25 * x5 + x29
    x31 = numpy.sqrt(abs(x2 + A[0]) ** 2 + abs(x3 + A[1]) ** 2 + abs(x4 + A[2]) ** 2)
    x32 = x13 * x31
    x33 = x11 * x32 + x18
    x34 = x15 * x32 + x21
    x35 = -x16 * x34
    x36 = x20 * x32 + x26
    x37 = -x16 * x36
    x38 = (2 * a + 2 * b) ** (-1.0)
    x39 = x13 * x38
    x40 = x15 * x39
    x41 = x11 * x39 + x31 * x33 + x35 - x40
    x42 = -x20 * x39 + x31 * x34 + x37 + x40
    x43 = -x16 * x42
    x44 = x33 * x38
    x45 = x34 * x38
    x46 = x19 * x31
    x47 = x24 + x44 - x45 + x46
    x48 = x36 * x38
    x49 = x22 * x31
    x50 = -x16 * (x28 + x45 - x48 + x49)
    x51 = 4 * x12
    x52 = x31 * x51
    x53 = x16 * x51
    x54 = 2 * x45

    # 10 item(s)
    S = numpy.array(
        [
            x30,
            -x16 * (x34 * x5 + x37) + x5 * (x33 * x5 + x35),
            x30,
            x41 * x5 + x43,
            x47 * x5 + x50,
            x30,
            x31 * x41
            + x38 * (x11 * x52 - x15 * x53)
            - x38 * (x15 * x52 - x20 * x53)
            + x43,
            x31 * x47 + x38 * (x19 + x41) - x38 * (x22 + x42) + x50,
            x25 * x31
            + x29
            + x38 * (-2 * x23 + 2 * x44 + 2 * x46 - x54)
            - x38 * (-2 * x27 - 2 * x48 + 2 * x49 + x54),
            x30,
        ]
    )
    return S


def coulomb_04(a, A, b, B, C):
    """Cartesian (sg) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = a + b
    x1 = x0 ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x1 * (a * A[1] + b * B[1])
    x4 = -x1 * (a * A[2] + b * B[2])
    x5 = numpy.sqrt(abs(x2 + B[0]) ** 2 + abs(x3 + B[1]) ** 2 + abs(x4 + B[2]) ** 2)
    x6 = x2 + C[0]
    x7 = x3 + C[1]
    x8 = x4 + C[2]
    x9 = x6 ** 2 + x7 ** 2 + x8 ** 2
    x10 = x0 * x9
    x11 = boys(0, x10)
    x12 = numpy.pi * x1 * numpy.exp(-a * b * x1 * x9)
    x13 = 2 * x12
    x14 = x13 * x5
    x15 = boys(1, x10)
    x16 = numpy.sqrt(abs(x6) ** 2 + abs(x7) ** 2 + abs(x8) ** 2)
    x17 = x13 * x16
    x18 = -x15 * x17
    x19 = x11 * x14 + x18
    x20 = boys(2, x10)
    x21 = -x17 * x20
    x22 = x14 * x15 + x21
    x23 = x16 * x22
    x24 = -x23
    x25 = x19 * x5 + x24
    x26 = boys(3, x10)
    x27 = -x17 * x26
    x28 = x14 * x20 + x27
    x29 = x16 * x28
    x30 = -x29
    x31 = x22 * x5 + x30
    x32 = x16 * x31
    x33 = -x32
    x34 = x25 * x5 + x33
    x35 = -x17 * boys(4, x10)
    x36 = x16 * (x14 * x26 + x35)
    x37 = -x36
    x38 = x16 * (x28 * x5 + x37)
    x39 = -x38
    x40 = -x16 * (x31 * x5 + x39)
    x41 = x34 * x5 + x40
    x42 = numpy.sqrt(abs(x2 + A[0]) ** 2 + abs(x3 + A[1]) ** 2 + abs(x4 + A[2]) ** 2)
    x43 = x13 * x42
    x44 = x11 * x43 + x18
    x45 = x15 * x43 + x21
    x46 = x16 * x45
    x47 = -x46
    x48 = x20 * x43 + x27
    x49 = x16 * x48
    x50 = -x49
    x51 = x45 * x5 + x50
    x52 = x26 * x43 + x35
    x53 = -x16 * x52
    x54 = (2 * a + 2 * b) ** (-1.0)
    x55 = x13 * x54
    x56 = x15 * x55
    x57 = x42 * x44
    x58 = x11 * x55 + x47 - x56 + x57
    x59 = x20 * x55
    x60 = x42 * x45
    x61 = x50 + x56 - x59 + x60
    x62 = -x16 * x61
    x63 = -x26 * x55 + x42 * x48 + x53 + x59
    x64 = -x16 * x63
    x65 = x44 * x54
    x66 = x45 * x54
    x67 = x19 * x42
    x68 = x24 + x65 - x66 + x67
    x69 = x48 * x54
    x70 = x22 * x42
    x71 = x30 + x66 - x69 + x70
    x72 = x16 * x71
    x73 = -x72
    x74 = x52 * x54
    x75 = x28 * x42
    x76 = x16 * (x37 + x69 - x74 + x75)
    x77 = -x76
    x78 = 4 * x12
    x79 = x42 * x78
    x80 = x16 * x78
    x81 = x54 * (x15 * x79 - x20 * x80)
    x82 = x42 * x58 + x54 * (x11 * x79 - x15 * x80) + x62 - x81
    x83 = x42 * x61 - x54 * (x20 * x79 - x26 * x80) + x64 + x81
    x84 = -x16 * x83
    x85 = x54 * (x19 + x58)
    x86 = x54 * (x22 + x61)
    x87 = x42 * x68
    x88 = x73 + x85 - x86 + x87
    x89 = x54 * (x28 + x63)
    x90 = x42 * x71
    x91 = -x16 * (x77 + x86 - x89 + x90)
    x92 = x25 * x42
    x93 = 2 * x66
    x94 = -2 * x23 + 2 * x65 + 2 * x67 - x93
    x95 = x54 * x94
    x96 = 2 * x69
    x97 = -2 * x29 + 2 * x70 + x93 - x96
    x98 = x54 * x97
    x99 = x33 + x92 + x95 - x98
    x100 = x31 * x42
    x101 = x54 * (-2 * x36 - 2 * x74 + 2 * x75 + x96)
    x102 = -x16 * (x100 - x101 + x39 + x98)
    x103 = 6 * x12 * x54
    x104 = x103 * x15
    x105 = 2 * x86
    x106 = 3 * x98

    # 15 item(s)
    S = numpy.array(
        [
            x41,
            -x16 * (-x16 * (x48 * x5 + x53) + x5 * x51)
            + x5 * (-x16 * x51 + x5 * (x44 * x5 + x47)),
            x41,
            -x16 * (x5 * x61 + x64) + x5 * (x5 * x58 + x62),
            -x16 * (x5 * x71 + x77) + x5 * (x5 * x68 + x73),
            x41,
            x5 * x82 + x84,
            x5 * x88 + x91,
            x102 + x5 * x99,
            x41,
            x42 * x82
            + x54 * (x103 * x11 - x104 - 3 * x46 + 3 * x57)
            - x54 * (-x103 * x20 + x104 - 3 * x49 + 3 * x60)
            + x84,
            x42 * x88 + x54 * (x82 + x94) - x54 * (x83 + x97) + x91,
            x102
            + x42 * x99
            + x54 * (-x105 + x25 - 2 * x72 + 2 * x85 + 2 * x87)
            - x54 * (x105 + x31 - 2 * x76 - 2 * x89 + 2 * x90),
            x34 * x42
            + x40
            - x54 * (3 * x100 - 3 * x101 + x106 - 3 * x38)
            + x54 * (-x106 - 3 * x32 + 3 * x92 + 3 * x95),
            x41,
        ]
    )
    return S


def coulomb_10(a, A, b, B, C):
    """Cartesian (ps) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = a + b
    x1 = x0 ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x1 * (a * A[1] + b * B[1])
    x4 = -x1 * (a * A[2] + b * B[2])
    x5 = x2 + C[0]
    x6 = x3 + C[1]
    x7 = x4 + C[2]
    x8 = x5 ** 2 + x6 ** 2 + x7 ** 2
    x9 = x0 * x8
    x10 = 2 * numpy.pi * x1 * numpy.exp(-a * b * x1 * x8)
    x11 = x10 * boys(0, x9)
    x12 = -x10 * numpy.sqrt(abs(x5) ** 2 + abs(x6) ** 2 + abs(x7) ** 2) * boys(1, x9)
    x13 = (
        x11 * numpy.sqrt(abs(x2 + A[0]) ** 2 + abs(x3 + A[1]) ** 2 + abs(x4 + A[2]) ** 2)
        + x12
    )

    # 3 item(s)
    S = numpy.array(
        [
            x13,
            x11
            * numpy.sqrt(abs(x2 + B[0]) ** 2 + abs(x3 + B[1]) ** 2 + abs(x4 + B[2]) ** 2)
            + x12,
            x13,
        ]
    )
    return S


def coulomb_11(a, A, b, B, C):
    """Cartesian (pp) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = a + b
    x1 = x0 ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x1 * (a * A[1] + b * B[1])
    x4 = -x1 * (a * A[2] + b * B[2])
    x5 = numpy.sqrt(abs(x2 + A[0]) ** 2 + abs(x3 + A[1]) ** 2 + abs(x4 + A[2]) ** 2)
    x6 = x2 + C[0]
    x7 = x3 + C[1]
    x8 = x4 + C[2]
    x9 = x6 ** 2 + x7 ** 2 + x8 ** 2
    x10 = x0 * x9
    x11 = boys(0, x10)
    x12 = numpy.sqrt(abs(x2 + B[0]) ** 2 + abs(x3 + B[1]) ** 2 + abs(x4 + B[2]) ** 2)
    x13 = 2 * numpy.pi * x1 * numpy.exp(-a * b * x1 * x9)
    x14 = x12 * x13
    x15 = boys(1, x10)
    x16 = numpy.sqrt(abs(x6) ** 2 + abs(x7) ** 2 + abs(x8) ** 2)
    x17 = x13 * x16
    x18 = -x15 * x17
    x19 = x11 * x14 + x18
    x20 = -x17 * boys(2, x10)
    x21 = -x16 * (x14 * x15 + x20)
    x22 = x19 * x5 + x21
    x23 = x13 * x5
    x24 = x11 * x23 + x18
    x25 = x15 * x23 + x20
    x26 = -x16 * x25
    x27 = x24 * x5 + x26
    x28 = x12 * x19 + x21
    x29 = (2 * a + 2 * b) ** (-1.0)

    # 9 item(s)
    S = numpy.array(
        [x22, x27, x22, x28, x12 * x24 + x26, x28, x22 + x24 * x29 - x25 * x29, x27, x22]
    )
    return S


def coulomb_12(a, A, b, B, C):
    """Cartesian (pd) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = a + b
    x1 = x0 ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x1 * (a * A[1] + b * B[1])
    x4 = -x1 * (a * A[2] + b * B[2])
    x5 = numpy.sqrt(abs(x2 + A[0]) ** 2 + abs(x3 + A[1]) ** 2 + abs(x4 + A[2]) ** 2)
    x6 = numpy.sqrt(abs(x2 + B[0]) ** 2 + abs(x3 + B[1]) ** 2 + abs(x4 + B[2]) ** 2)
    x7 = x2 + C[0]
    x8 = x3 + C[1]
    x9 = x4 + C[2]
    x10 = x7 ** 2 + x8 ** 2 + x9 ** 2
    x11 = x0 * x10
    x12 = boys(0, x11)
    x13 = 2 * numpy.pi * x1 * numpy.exp(-a * b * x1 * x10)
    x14 = x13 * x6
    x15 = boys(1, x11)
    x16 = numpy.sqrt(abs(x7) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x17 = x13 * x16
    x18 = -x15 * x17
    x19 = x12 * x14 + x18
    x20 = boys(2, x11)
    x21 = -x17 * x20
    x22 = x14 * x15 + x21
    x23 = x16 * x22
    x24 = -x23
    x25 = x19 * x6 + x24
    x26 = -x17 * boys(3, x11)
    x27 = x16 * (x14 * x20 + x26)
    x28 = -x27
    x29 = -x16 * (x22 * x6 + x28)
    x30 = x25 * x5 + x29
    x31 = x13 * x5
    x32 = x12 * x31 + x18
    x33 = x15 * x31 + x21
    x34 = -x16 * x33
    x35 = x32 * x6 + x34
    x36 = x20 * x31 + x26
    x37 = -x16 * x36
    x38 = -x16 * (x33 * x6 + x37)
    x39 = x35 * x5 + x38
    x40 = (2 * a + 2 * b) ** (-1.0)
    x41 = x13 * x40
    x42 = x15 * x41
    x43 = x32 * x5 + x34
    x44 = x12 * x41 - x42 + x43
    x45 = x33 * x5 + x37
    x46 = -x16 * (-x20 * x41 + x42 + x45)
    x47 = x44 * x5 + x46
    x48 = x32 * x40
    x49 = x33 * x40
    x50 = x19 * x5
    x51 = x24 + x50
    x52 = x48 - x49 + x51
    x53 = x36 * x40
    x54 = x22 * x5
    x55 = x28 + x54
    x56 = -x16 * (x49 - x53 + x55)
    x57 = x5 * x52 + x56
    x58 = x25 * x6 + x29
    x59 = 2 * x49

    # 18 item(s)
    S = numpy.array(
        [
            x30,
            x39,
            x30,
            x47,
            x57,
            x30,
            x58,
            x35 * x6 + x38,
            x58,
            x44 * x6 + x46,
            x52 * x6 + x56,
            x58,
            x30
            + x40 * (-2 * x23 + 2 * x48 + 2 * x50 - x59)
            - x40 * (-2 * x27 - 2 * x53 + 2 * x54 + x59),
            x39 + x40 * x43 - x40 * x45,
            x30 + x40 * x51 - x40 * x55,
            x47,
            x57,
            x30,
        ]
    )
    return S


def coulomb_13(a, A, b, B, C):
    """Cartesian (pf) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = a + b
    x1 = x0 ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x1 * (a * A[1] + b * B[1])
    x4 = -x1 * (a * A[2] + b * B[2])
    x5 = numpy.sqrt(abs(x2 + A[0]) ** 2 + abs(x3 + A[1]) ** 2 + abs(x4 + A[2]) ** 2)
    x6 = numpy.sqrt(abs(x2 + B[0]) ** 2 + abs(x3 + B[1]) ** 2 + abs(x4 + B[2]) ** 2)
    x7 = x2 + C[0]
    x8 = x3 + C[1]
    x9 = x4 + C[2]
    x10 = x7 ** 2 + x8 ** 2 + x9 ** 2
    x11 = x0 * x10
    x12 = boys(0, x11)
    x13 = numpy.pi * x1 * numpy.exp(-a * b * x1 * x10)
    x14 = 2 * x13
    x15 = x14 * x6
    x16 = boys(1, x11)
    x17 = numpy.sqrt(abs(x7) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x18 = x14 * x17
    x19 = -x16 * x18
    x20 = x12 * x15 + x19
    x21 = boys(2, x11)
    x22 = -x18 * x21
    x23 = x15 * x16 + x22
    x24 = x17 * x23
    x25 = -x24
    x26 = x20 * x6 + x25
    x27 = boys(3, x11)
    x28 = -x18 * x27
    x29 = x15 * x21 + x28
    x30 = x17 * x29
    x31 = -x30
    x32 = x23 * x6 + x31
    x33 = x17 * x32
    x34 = -x33
    x35 = x26 * x6 + x34
    x36 = -x18 * boys(4, x11)
    x37 = x17 * (x15 * x27 + x36)
    x38 = -x37
    x39 = x17 * (x29 * x6 + x38)
    x40 = -x39
    x41 = -x17 * (x32 * x6 + x40)
    x42 = x35 * x5 + x41
    x43 = 2 * x5
    x44 = x13 * x43
    x45 = x12 * x44 + x19
    x46 = x16 * x44 + x22
    x47 = -x17 * x46
    x48 = x45 * x6 + x47
    x49 = x21 * x44 + x28
    x50 = -x17 * x49
    x51 = x46 * x6 + x50
    x52 = x17 * x51
    x53 = x48 * x6 - x52
    x54 = x27 * x44 + x36
    x55 = -x17 * x54
    x56 = x17 * (x49 * x6 + x55)
    x57 = -x17 * (x51 * x6 - x56)
    x58 = x5 * x53 + x57
    x59 = (2 * a + 2 * b) ** (-1.0)
    x60 = 2 * x59
    x61 = x13 * x60
    x62 = x16 * x61
    x63 = x45 * x5 + x47
    x64 = x12 * x61 - x62 + x63
    x65 = x21 * x61
    x66 = x46 * x5 + x50
    x67 = x62 - x65 + x66
    x68 = -x17 * x67
    x69 = x6 * x64 + x68
    x70 = x49 * x5 + x55
    x71 = -x27 * x61 + x65 + x70
    x72 = -x17 * x71
    x73 = -x17 * (x6 * x67 + x72)
    x74 = x5 * x69 + x73
    x75 = x45 * x59
    x76 = x46 * x59
    x77 = x20 * x5
    x78 = x25 + x77
    x79 = x75 - x76 + x78
    x80 = x49 * x59
    x81 = x23 * x5
    x82 = x31 + x81
    x83 = x76 - x80 + x82
    x84 = -x17 * x83
    x85 = x6 * x79 + x84
    x86 = x54 * x59
    x87 = x29 * x5
    x88 = x38 + x87
    x89 = -x17 * (x80 - x86 + x88)
    x90 = -x17 * (x6 * x83 + x89)
    x91 = x5 * x85 + x90
    x92 = 4 * x13
    x93 = x5 * x92
    x94 = x17 * x92
    x95 = x59 * (x16 * x93 - x21 * x94)
    x96 = x5 * x64 + x68
    x97 = x59 * (x12 * x93 - x16 * x94) - x95 + x96
    x98 = x5 * x67 + x72
    x99 = -x17 * (-x59 * (x21 * x93 - x27 * x94) + x95 + x98)
    x100 = x5 * x97 + x99
    x101 = x59 * (x23 + x67)
    x102 = x5 * x79 + x84
    x103 = -x101 + x102 + x59 * (x20 + x64)
    x104 = x5 * x83 + x89
    x105 = -x17 * (x101 + x104 - x59 * (x29 + x71))
    x106 = x103 * x5 + x105
    x107 = 2 * x76
    x108 = x59 * (-x107 - 2 * x24 + 2 * x75 + 2 * x77)
    x109 = 2 * x80
    x110 = x59 * (x107 - x109 - 2 * x30 + 2 * x81)
    x111 = x26 * x5
    x112 = x111 + x34
    x113 = x108 - x110 + x112
    x114 = x59 * (x109 - 2 * x37 - 2 * x86 + 2 * x87)
    x115 = x32 * x5
    x116 = x115 + x40
    x117 = -x17 * (x110 - x114 + x116)
    x118 = x113 * x5 + x117
    x119 = x35 * x6 + x41
    x120 = 3 * x110
    x121 = x60 * x66
    x122 = x60 * x82

    # 30 item(s)
    S = numpy.array(
        [
            x42,
            x58,
            x42,
            x74,
            x91,
            x42,
            x100,
            x106,
            x118,
            x42,
            x119,
            x53 * x6 + x57,
            x119,
            x6 * x69 + x73,
            x6 * x85 + x90,
            x119,
            x6 * x97 + x99,
            x103 * x6 + x105,
            x113 * x6 + x117,
            x119,
            x42
            + x59 * (3 * x108 + 3 * x111 - x120 - 3 * x33)
            - x59 * (-3 * x114 + 3 * x115 + x120 - 3 * x39),
            x58
            + x59 * (-x121 + x43 * x48 - 2 * x52 + x60 * x63)
            - x59 * (x121 + x43 * x51 - 2 * x56 - x60 * x70),
            x42
            + x59 * (2 * x111 - x122 - 2 * x33 + x60 * x78)
            - x59 * (2 * x115 + x122 - 2 * x39 - x60 * x88),
            x59 * x96 - x59 * x98 + x74,
            x102 * x59 - x104 * x59 + x91,
            x112 * x59 - x116 * x59 + x42,
            x100,
            x106,
            x118,
            x42,
        ]
    )
    return S


def coulomb_14(a, A, b, B, C):
    """Cartesian (pg) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = a + b
    x1 = x0 ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x1 * (a * A[1] + b * B[1])
    x4 = -x1 * (a * A[2] + b * B[2])
    x5 = numpy.sqrt(abs(x2 + A[0]) ** 2 + abs(x3 + A[1]) ** 2 + abs(x4 + A[2]) ** 2)
    x6 = numpy.sqrt(abs(x2 + B[0]) ** 2 + abs(x3 + B[1]) ** 2 + abs(x4 + B[2]) ** 2)
    x7 = x2 + C[0]
    x8 = x3 + C[1]
    x9 = x4 + C[2]
    x10 = x7 ** 2 + x8 ** 2 + x9 ** 2
    x11 = x0 * x10
    x12 = boys(0, x11)
    x13 = numpy.pi * x1 * numpy.exp(-a * b * x1 * x10)
    x14 = 2 * x13
    x15 = x14 * x6
    x16 = boys(1, x11)
    x17 = numpy.sqrt(abs(x7) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x18 = x14 * x17
    x19 = -x16 * x18
    x20 = x12 * x15 + x19
    x21 = boys(2, x11)
    x22 = -x18 * x21
    x23 = x15 * x16 + x22
    x24 = x17 * x23
    x25 = -x24
    x26 = x20 * x6 + x25
    x27 = boys(3, x11)
    x28 = -x18 * x27
    x29 = x15 * x21 + x28
    x30 = x17 * x29
    x31 = -x30
    x32 = x23 * x6 + x31
    x33 = x17 * x32
    x34 = -x33
    x35 = x26 * x6 + x34
    x36 = boys(4, x11)
    x37 = -x18 * x36
    x38 = x15 * x27 + x37
    x39 = x17 * x38
    x40 = -x39
    x41 = x29 * x6 + x40
    x42 = x17 * x41
    x43 = -x42
    x44 = x32 * x6 + x43
    x45 = x17 * x44
    x46 = -x45
    x47 = x35 * x6 + x46
    x48 = -x18 * boys(5, x11)
    x49 = x17 * (x15 * x36 + x48)
    x50 = -x49
    x51 = x17 * (x38 * x6 + x50)
    x52 = -x51
    x53 = x17 * (x41 * x6 + x52)
    x54 = -x53
    x55 = -x17 * (x44 * x6 + x54)
    x56 = x47 * x5 + x55
    x57 = 2 * x5
    x58 = x13 * x57
    x59 = x12 * x58 + x19
    x60 = x16 * x58 + x22
    x61 = x17 * x60
    x62 = -x61
    x63 = x59 * x6 + x62
    x64 = x21 * x58 + x28
    x65 = x17 * x64
    x66 = -x65
    x67 = x6 * x60 + x66
    x68 = x17 * x67
    x69 = x6 * x63 - x68
    x70 = x27 * x58 + x37
    x71 = x17 * x70
    x72 = -x71
    x73 = x6 * x64 + x72
    x74 = x17 * x73
    x75 = x6 * x67 - x74
    x76 = x17 * x75
    x77 = x6 * x69 - x76
    x78 = x36 * x58 + x48
    x79 = -x17 * x78
    x80 = x17 * (x6 * x70 + x79)
    x81 = x17 * (x6 * x73 - x80)
    x82 = -x17 * (x6 * x75 - x81)
    x83 = x5 * x77 + x82
    x84 = (2 * a + 2 * b) ** (-1.0)
    x85 = 2 * x84
    x86 = x13 * x85
    x87 = x16 * x86
    x88 = x5 * x59
    x89 = x62 + x88
    x90 = x12 * x86 - x87 + x89
    x91 = x21 * x86
    x92 = x5 * x60
    x93 = x66 + x92
    x94 = x87 - x91 + x93
    x95 = -x17 * x94
    x96 = x6 * x90 + x95
    x97 = x27 * x86
    x98 = x5 * x64
    x99 = x72 + x98
    x100 = x91 - x97 + x99
    x101 = -x100 * x17
    x102 = x101 + x6 * x94
    x103 = x102 * x17
    x104 = -x103 + x6 * x96
    x105 = x5 * x70 + x79
    x106 = x105 - x36 * x86 + x97
    x107 = -x106 * x17
    x108 = x17 * (x100 * x6 + x107)
    x109 = -x17 * (x102 * x6 - x108)
    x110 = x104 * x5 + x109
    x111 = x59 * x84
    x112 = x60 * x84
    x113 = x20 * x5
    x114 = x113 + x25
    x115 = x111 - x112 + x114
    x116 = x64 * x84
    x117 = x23 * x5
    x118 = x117 + x31
    x119 = x112 - x116 + x118
    x120 = x119 * x17
    x121 = -x120
    x122 = x115 * x6 + x121
    x123 = x70 * x84
    x124 = x29 * x5
    x125 = x124 + x40
    x126 = x116 - x123 + x125
    x127 = x126 * x17
    x128 = -x127
    x129 = x119 * x6 + x128
    x130 = x129 * x17
    x131 = x122 * x6 - x130
    x132 = x78 * x84
    x133 = x38 * x5
    x134 = x133 + x50
    x135 = x17 * (x123 - x132 + x134)
    x136 = -x135
    x137 = x17 * (x126 * x6 + x136)
    x138 = -x17 * (x129 * x6 - x137)
    x139 = x131 * x5 + x138
    x140 = 4 * x13
    x141 = x140 * x5
    x142 = x140 * x17
    x143 = x84 * (x141 * x16 - x142 * x21)
    x144 = x5 * x90 + x95
    x145 = -x143 + x144 + x84 * (x12 * x141 - x142 * x16)
    x146 = x84 * (x141 * x21 - x142 * x27)
    x147 = x101 + x5 * x94
    x148 = x143 - x146 + x147
    x149 = -x148 * x17
    x150 = x145 * x6 + x149
    x151 = x100 * x5 + x107
    x152 = x146 + x151 - x84 * (x141 * x27 - x142 * x36)
    x153 = -x152 * x17
    x154 = -x17 * (x148 * x6 + x153)
    x155 = x150 * x5 + x154
    x156 = x84 * (x20 + x90)
    x157 = x84 * (x23 + x94)
    x158 = x115 * x5
    x159 = x121 + x158
    x160 = x156 - x157 + x159
    x161 = x84 * (x100 + x29)
    x162 = x119 * x5
    x163 = x128 + x162
    x164 = x157 - x161 + x163
    x165 = -x164 * x17
    x166 = x160 * x6 + x165
    x167 = x84 * (x106 + x38)
    x168 = x126 * x5
    x169 = x136 + x168
    x170 = -x17 * (x161 - x167 + x169)
    x171 = -x17 * (x164 * x6 + x170)
    x172 = x166 * x5 + x171
    x173 = 2 * x112
    x174 = 2 * x111 + 2 * x113 - x173 - 2 * x24
    x175 = x174 * x84
    x176 = 2 * x116
    x177 = 2 * x117 + x173 - x176 - 2 * x30
    x178 = x177 * x84
    x179 = x26 * x5
    x180 = x179 + x34
    x181 = x175 - x178 + x180
    x182 = 2 * x123
    x183 = 2 * x124 + x176 - x182 - 2 * x39
    x184 = x183 * x84
    x185 = x32 * x5
    x186 = x185 + x43
    x187 = x178 - x184 + x186
    x188 = -x17 * x187
    x189 = x181 * x6 + x188
    x190 = x84 * (-2 * x132 + 2 * x133 + x182 - 2 * x49)
    x191 = x41 * x5
    x192 = x191 + x52
    x193 = -x17 * (x184 - x190 + x192)
    x194 = -x17 * (x187 * x6 + x193)
    x195 = x189 * x5 + x194
    x196 = 6 * x13 * x84
    x197 = x16 * x196
    x198 = x196 * x21
    x199 = x84 * (x197 - x198 - 3 * x65 + 3 * x92)
    x200 = x145 * x5 + x149
    x201 = -x199 + x200 + x84 * (x12 * x196 - x197 - 3 * x61 + 3 * x88)
    x202 = x148 * x5 + x153
    x203 = -x17 * (x199 + x202 - x84 * (-x196 * x27 + x198 - 3 * x71 + 3 * x98))
    x204 = x201 * x5 + x203
    x205 = x84 * (x148 + x177)
    x206 = x160 * x5 + x165
    x207 = -x205 + x206 + x84 * (x145 + x174)
    x208 = x164 * x5 + x170
    x209 = -x17 * (x205 + x208 - x84 * (x152 + x183))
    x210 = x207 * x5 + x209
    x211 = 2 * x157
    x212 = 2 * x161
    x213 = x84 * (-2 * x127 + 2 * x162 + x211 - x212 + x32)
    x214 = x181 * x5 + x188
    x215 = -x213 + x214 + x84 * (-2 * x120 + 2 * x156 + 2 * x158 - x211 + x26)
    x216 = x187 * x5 + x193
    x217 = -x17 * (x213 + x216 - x84 * (-2 * x135 - 2 * x167 + 2 * x168 + x212 + x41))
    x218 = x215 * x5 + x217
    x219 = 3 * x178
    x220 = x84 * (3 * x175 + 3 * x179 - x219 - 3 * x33)
    x221 = 3 * x184
    x222 = x84 * (3 * x185 + x219 - x221 - 3 * x42)
    x223 = x35 * x5
    x224 = x223 + x46
    x225 = x220 - x222 + x224
    x226 = x84 * (-3 * x190 + 3 * x191 + x221 - 3 * x51)
    x227 = x44 * x5
    x228 = x227 + x54
    x229 = -x17 * (x222 - x226 + x228)
    x230 = x225 * x5 + x229
    x231 = x47 * x6 + x55
    x232 = 4 * x222
    x233 = 3 * x5
    x234 = x85 * x93
    x235 = 3 * x84
    x236 = x85 * x99
    x237 = x235 * (x234 - x236 + x57 * x67 - 2 * x74)
    x238 = x118 * x85
    x239 = x125 * x85
    x240 = x235 * (2 * x185 + x238 - x239 - 2 * x42)
    x241 = x147 * x85
    x242 = x163 * x85
    x243 = x186 * x85

    # 45 item(s)
    S = numpy.array(
        [
            x56,
            x83,
            x56,
            x110,
            x139,
            x56,
            x155,
            x172,
            x195,
            x56,
            x204,
            x210,
            x218,
            x230,
            x56,
            x231,
            x6 * x77 + x82,
            x231,
            x104 * x6 + x109,
            x131 * x6 + x138,
            x231,
            x150 * x6 + x154,
            x166 * x6 + x171,
            x189 * x6 + x194,
            x231,
            x201 * x6 + x203,
            x207 * x6 + x209,
            x215 * x6 + x217,
            x225 * x6 + x229,
            x231,
            x56
            + x84 * (4 * x220 + 4 * x223 - x232 - 4 * x45)
            - x84 * (-4 * x226 + 4 * x227 + x232 - 4 * x53),
            x83
            + x84
            * (
                x233 * x69
                + x235 * (-x234 + x57 * x63 - 2 * x68 + x85 * x89)
                - x237
                - 3 * x76
            )
            - x84
            * (
                x233 * x75
                - x235 * (-x105 * x85 + x236 + x57 * x73 - 2 * x80)
                + x237
                - 3 * x81
            ),
            x56
            + x84
            * (
                3 * x223
                + x235 * (x114 * x85 + 2 * x179 - x238 - 2 * x33)
                - x240
                - 3 * x45
            )
            - x84
            * (
                3 * x227
                - x235 * (-x134 * x85 + 2 * x191 + x239 - 2 * x51)
                + x240
                - 3 * x53
            ),
            x110
            + x84 * (-2 * x103 + x144 * x85 - x241 + x57 * x96)
            - x84 * (x102 * x57 - 2 * x108 - x151 * x85 + x241),
            x139
            + x84 * (x122 * x57 - 2 * x130 + x159 * x85 - x242)
            - x84 * (x129 * x57 - 2 * x137 - x169 * x85 + x242),
            x56
            + x84 * (x180 * x85 + 2 * x223 - x243 - 2 * x45)
            - x84 * (-x192 * x85 + 2 * x227 + x243 - 2 * x53),
            x155 + x200 * x84 - x202 * x84,
            x172 + x206 * x84 - x208 * x84,
            x195 + x214 * x84 - x216 * x84,
            x224 * x84 - x228 * x84 + x56,
            x204,
            x210,
            x218,
            x230,
            x56,
        ]
    )
    return S


def coulomb_20(a, A, b, B, C):
    """Cartesian (ds) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = a + b
    x1 = x0 ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = x2 + C[0]
    x4 = -x1 * (a * A[1] + b * B[1])
    x5 = x4 + C[1]
    x6 = -x1 * (a * A[2] + b * B[2])
    x7 = x6 + C[2]
    x8 = x3 ** 2 + x5 ** 2 + x7 ** 2
    x9 = x0 * x8
    x10 = boys(1, x9)
    x11 = (2 * a + 2 * b) ** (-1.0)
    x12 = 2 * numpy.pi * x1 * numpy.exp(-a * b * x1 * x8)
    x13 = x11 * x12
    x14 = boys(0, x9)
    x15 = numpy.sqrt(abs(x2 + A[0]) ** 2 + abs(x4 + A[1]) ** 2 + abs(x6 + A[2]) ** 2)
    x16 = x12 * x15
    x17 = numpy.sqrt(abs(x3) ** 2 + abs(x5) ** 2 + abs(x7) ** 2)
    x18 = x12 * x17
    x19 = -x10 * x18
    x20 = x14 * x16 + x19
    x21 = -x18 * boys(2, x9)
    x22 = x10 * x16 + x21
    x23 = -x17 * x22
    x24 = x15 * x20 + x23
    x25 = -x10 * x13 + x13 * x14 + x24
    x26 = numpy.sqrt(abs(x2 + B[0]) ** 2 + abs(x4 + B[1]) ** 2 + abs(x6 + B[2]) ** 2)
    x27 = x12 * x26
    x28 = x14 * x27 + x19
    x29 = -x17 * (x10 * x27 + x21)

    # 6 item(s)
    S = numpy.array(
        [
            x25,
            x11 * x20 - x11 * x22 + x15 * x28 + x29,
            x24,
            x26 * x28 + x29,
            x20 * x26 + x23,
            x25,
        ]
    )
    return S


def coulomb_21(a, A, b, B, C):
    """Cartesian (dp) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = x3 + C[0]
    x5 = -x2 * (a * A[1] + b * B[1])
    x6 = x5 + C[1]
    x7 = -x2 * (a * A[2] + b * B[2])
    x8 = x7 + C[2]
    x9 = x4 ** 2 + x6 ** 2 + x8 ** 2
    x10 = x1 * x9
    x11 = boys(0, x10)
    x12 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x5 + B[1]) ** 2 + abs(x7 + B[2]) ** 2)
    x13 = 2 * numpy.pi * x2 * numpy.exp(-a * b * x2 * x9)
    x14 = x12 * x13
    x15 = boys(1, x10)
    x16 = numpy.sqrt(abs(x4) ** 2 + abs(x6) ** 2 + abs(x8) ** 2)
    x17 = x13 * x16
    x18 = -x15 * x17
    x19 = x11 * x14 + x18
    x20 = boys(2, x10)
    x21 = -x17 * x20
    x22 = x14 * x15 + x21
    x23 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x5 + A[1]) ** 2 + abs(x7 + A[2]) ** 2)
    x24 = -x16 * x22
    x25 = x19 * x23 + x24
    x26 = -x17 * boys(3, x10)
    x27 = -x16 * (x14 * x20 + x26)
    x28 = x22 * x23 + x27
    x29 = -x16 * x28
    x30 = x23 * x25 + x29
    x31 = x0 * x19 - x0 * x22 + x30
    x32 = x13 * x23
    x33 = x11 * x32 + x18
    x34 = x15 * x32 + x21
    x35 = -x16 * x34
    x36 = x23 * x33 + x35
    x37 = x20 * x32 + x26
    x38 = -x16 * x37
    x39 = x23 * x34 + x38
    x40 = -x16 * x39
    x41 = x23 * x36 + x40
    x42 = x0 * x34
    x43 = x0 * x33 - x42
    x44 = x41 + x43
    x45 = x12 * x19 + x24
    x46 = -x16 * (x12 * x22 + x27)
    x47 = x0 * x25 - x0 * x28 + x23 * x45 + x46
    x48 = x12 * x33 + x35
    x49 = -x16 * (x12 * x34 + x38)
    x50 = x25 + x43
    x51 = -x16 * (-x0 * x37 + x28 + x42)
    x52 = x23 * x50 + x51
    x53 = x12 * x45 + x46
    x54 = x0 * x13
    x55 = x15 * x54

    # 18 item(s)
    S = numpy.array(
        [
            x31,
            x44,
            x31,
            x47,
            x0 * x36 - x0 * x39 + x23 * x48 + x49,
            x47,
            x52,
            x41,
            x30,
            x53,
            x12 * x48 + x49,
            x53,
            x12 * x50 + x51,
            x12 * x36 + x40,
            x12 * x25 + x29,
            x0 * (x11 * x54 + x19 + x36 - x55)
            - x0 * (-x20 * x54 + x22 + x39 + x55)
            + x52,
            x44,
            x31,
        ]
    )
    return S


def coulomb_22(a, A, b, B, C):
    """Cartesian (dd) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = -x2 * (a * A[1] + b * B[1])
    x5 = -x2 * (a * A[2] + b * B[2])
    x6 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x4 + B[1]) ** 2 + abs(x5 + B[2]) ** 2)
    x7 = x3 + C[0]
    x8 = x4 + C[1]
    x9 = x5 + C[2]
    x10 = x7 ** 2 + x8 ** 2 + x9 ** 2
    x11 = x1 * x10
    x12 = boys(0, x11)
    x13 = numpy.pi * x2 * numpy.exp(-a * b * x10 * x2)
    x14 = 2 * x13
    x15 = x14 * x6
    x16 = boys(1, x11)
    x17 = numpy.sqrt(abs(x7) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x18 = x14 * x17
    x19 = -x16 * x18
    x20 = x12 * x15 + x19
    x21 = boys(2, x11)
    x22 = -x18 * x21
    x23 = x15 * x16 + x22
    x24 = x17 * x23
    x25 = -x24
    x26 = x20 * x6 + x25
    x27 = boys(3, x11)
    x28 = -x18 * x27
    x29 = x15 * x21 + x28
    x30 = x17 * x29
    x31 = -x30
    x32 = x23 * x6 + x31
    x33 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x4 + A[1]) ** 2 + abs(x5 + A[2]) ** 2)
    x34 = -x17 * x32
    x35 = x26 * x33 + x34
    x36 = -x18 * boys(4, x11)
    x37 = x17 * (x15 * x27 + x36)
    x38 = -x37
    x39 = -x17 * (x29 * x6 + x38)
    x40 = x32 * x33 + x39
    x41 = -x17 * x40
    x42 = x33 * x35 + x41
    x43 = x0 * x26 - x0 * x32 + x42
    x44 = x14 * x33
    x45 = x12 * x44 + x19
    x46 = x16 * x44 + x22
    x47 = -x17 * x46
    x48 = x45 * x6 + x47
    x49 = x21 * x44 + x28
    x50 = -x17 * x49
    x51 = x46 * x6 + x50
    x52 = -x17 * x51
    x53 = x33 * x48 + x52
    x54 = x27 * x44 + x36
    x55 = -x17 * x54
    x56 = -x17 * (x49 * x6 + x55)
    x57 = x33 * x51 + x56
    x58 = 2 * x0
    x59 = x13 * x58
    x60 = x16 * x59
    x61 = x33 * x45 + x47
    x62 = x12 * x59 - x60 + x61
    x63 = x21 * x59
    x64 = x33 * x46 + x50
    x65 = x60 - x63 + x64
    x66 = -x17 * x65
    x67 = x33 * x62 + x66
    x68 = x33 * x49 + x55
    x69 = -x27 * x59 + x63 + x68
    x70 = -x17 * x69
    x71 = x33 * x65 + x70
    x72 = -x17 * x71
    x73 = x33 * x67 + x72
    x74 = x0 * x62 - x0 * x65 + x73
    x75 = x20 * x33
    x76 = x25 + x75
    x77 = x0 * x45
    x78 = x0 * x46
    x79 = x77 - x78
    x80 = x76 + x79
    x81 = x23 * x33
    x82 = x31 + x81
    x83 = x0 * x49
    x84 = x78 - x83
    x85 = x82 + x84
    x86 = x33 * x80
    x87 = x17 * x85
    x88 = -x87
    x89 = x86 + x88
    x90 = x33 * x85
    x91 = x0 * x54
    x92 = x29 * x33
    x93 = x38 + x92
    x94 = x17 * (x83 - x91 + x93)
    x95 = -x94
    x96 = x90 + x95
    x97 = -x17 * x96
    x98 = x33 * x89 + x97
    x99 = x0 * x80 - x0 * x85 + x98
    x100 = x26 * x6 + x34
    x101 = -x17 * (x32 * x6 + x39)
    x102 = x0 * x35 - x0 * x40 + x100 * x33 + x101
    x103 = x48 * x6 + x52
    x104 = -x17 * (x51 * x6 + x56)
    x105 = x6 * x62 + x66
    x106 = -x17 * (x6 * x65 + x70)
    x107 = x6 * x80 + x88
    x108 = -x17 * (x6 * x85 + x95)
    x109 = 2 * x78
    x110 = 2 * x83
    x111 = x0 * (x109 - x110 - 2 * x30 + 2 * x81)
    x112 = x0 * (-x109 - 2 * x24 + 2 * x75 + 2 * x77) - x111 + x35
    x113 = -x17 * (-x0 * (x110 - 2 * x37 - 2 * x91 + 2 * x92) + x111 + x40)
    x114 = x112 * x33 + x113
    x115 = x0 * x64
    x116 = x0 * x61 - x115 + x53
    x117 = -x17 * (-x0 * x68 + x115 + x57)
    x118 = x116 * x33 + x117
    x119 = x0 * x82
    x120 = x0 * x76 - x119 + x35
    x121 = -x17 * (-x0 * x93 + x119 + x40)
    x122 = x120 * x33 + x121
    x123 = x100 * x6 + x101
    x124 = x58 * (x23 + x65)
    x125 = x0 * x23

    # 36 item(s)
    S = numpy.array(
        [
            x43,
            x0 * x48 - x0 * x51 - x17 * x57 + x33 * x53,
            x43,
            x74,
            x99,
            x43,
            x102,
            x0 * x53 - x0 * x57 + x103 * x33 + x104,
            x102,
            x0 * x67 - x0 * x71 + x105 * x33 + x106,
            x0 * x89 - x0 * x96 + x107 * x33 + x108,
            x102,
            x114,
            x118,
            x122,
            x73,
            x98,
            x42,
            x123,
            x103 * x6 + x104,
            x123,
            x105 * x6 + x106,
            x107 * x6 + x108,
            x123,
            x112 * x6 + x113,
            x116 * x6 + x117,
            x120 * x6 + x121,
            x6 * x67 + x72,
            x6 * x89 + x97,
            x35 * x6 + x41,
            x0 * (-x124 + x26 + x58 * (x20 + x62) + 2 * x86 - 2 * x87)
            - x0 * (x124 + x32 - x58 * (x29 + x69) + 2 * x90 - 2 * x94)
            + x114,
            x0 * (-x17 * x64 + x33 * x61 + x48 + x79)
            - x0 * (-x17 * x68 + x33 * x64 + x51 + x84)
            + x118,
            x0 * (x0 * x20 - x125 - x17 * x82 + x26 + x33 * x76)
            - x0 * (-x0 * x29 + x125 - x17 * x93 + x32 + x33 * x82)
            + x122,
            x74,
            x99,
            x43,
        ]
    )
    return S


def coulomb_23(a, A, b, B, C):
    """Cartesian (df) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = -x2 * (a * A[1] + b * B[1])
    x5 = -x2 * (a * A[2] + b * B[2])
    x6 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x4 + B[1]) ** 2 + abs(x5 + B[2]) ** 2)
    x7 = x3 + C[0]
    x8 = x4 + C[1]
    x9 = x5 + C[2]
    x10 = x7 ** 2 + x8 ** 2 + x9 ** 2
    x11 = x1 * x10
    x12 = numpy.pi * x2 * numpy.exp(-a * b * x10 * x2)
    x13 = x12 * boys(0, x11)
    x14 = 2 * x6
    x15 = boys(1, x11)
    x16 = numpy.sqrt(abs(x7) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x17 = 2 * x16
    x18 = x12 * x17
    x19 = -x15 * x18
    x20 = x13 * x14 + x19
    x21 = x12 * x14
    x22 = boys(2, x11)
    x23 = -x18 * x22
    x24 = x15 * x21 + x23
    x25 = x16 * x24
    x26 = -x25
    x27 = x20 * x6 + x26
    x28 = boys(3, x11)
    x29 = -x18 * x28
    x30 = x21 * x22 + x29
    x31 = x16 * x30
    x32 = -x31
    x33 = x24 * x6 + x32
    x34 = x16 * x33
    x35 = -x34
    x36 = x27 * x6 + x35
    x37 = boys(4, x11)
    x38 = -x18 * x37
    x39 = x21 * x28 + x38
    x40 = x16 * x39
    x41 = -x40
    x42 = x30 * x6 + x41
    x43 = x16 * x42
    x44 = -x43
    x45 = x33 * x6 + x44
    x46 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x4 + A[1]) ** 2 + abs(x5 + A[2]) ** 2)
    x47 = -x16 * x45
    x48 = x36 * x46 + x47
    x49 = -x18 * boys(5, x11)
    x50 = x16 * (x21 * x37 + x49)
    x51 = -x50
    x52 = x16 * (x39 * x6 + x51)
    x53 = -x52
    x54 = -x16 * (x42 * x6 + x53)
    x55 = x45 * x46 + x54
    x56 = -x16 * x55
    x57 = x46 * x48 + x56
    x58 = x0 * x36 - x0 * x45 + x57
    x59 = 2 * x46
    x60 = x13 * x59 + x19
    x61 = x12 * x59
    x62 = x15 * x61 + x23
    x63 = -x16 * x62
    x64 = x6 * x60 + x63
    x65 = x22 * x61 + x29
    x66 = -x16 * x65
    x67 = x6 * x62 + x66
    x68 = x16 * x67
    x69 = -x68
    x70 = x6 * x64 + x69
    x71 = x28 * x61 + x38
    x72 = -x16 * x71
    x73 = x6 * x65 + x72
    x74 = x16 * x73
    x75 = -x74
    x76 = x6 * x67 + x75
    x77 = -x16 * x76
    x78 = x46 * x70 + x77
    x79 = x37 * x61 + x49
    x80 = -x16 * x79
    x81 = x16 * (x6 * x71 + x80)
    x82 = -x81
    x83 = -x16 * (x6 * x73 + x82)
    x84 = x46 * x76 + x83
    x85 = 2 * x0
    x86 = x12 * x85
    x87 = x15 * x86
    x88 = x46 * x60 + x63
    x89 = x13 * x85 - x87 + x88
    x90 = x22 * x86
    x91 = x46 * x62 + x66
    x92 = x87 - x90 + x91
    x93 = -x16 * x92
    x94 = x6 * x89 + x93
    x95 = x28 * x86
    x96 = x46 * x65 + x72
    x97 = x90 - x95 + x96
    x98 = -x16 * x97
    x99 = x6 * x92 + x98
    x100 = -x16 * x99
    x101 = x100 + x46 * x94
    x102 = x46 * x71 + x80
    x103 = x102 - x37 * x86 + x95
    x104 = -x103 * x16
    x105 = -x16 * (x104 + x6 * x97)
    x106 = x105 + x46 * x99
    x107 = x20 * x46
    x108 = x107 + x26
    x109 = x0 * x60
    x110 = x0 * x62
    x111 = x109 - x110
    x112 = x108 + x111
    x113 = x24 * x46
    x114 = x113 + x32
    x115 = x0 * x65
    x116 = x110 - x115
    x117 = x114 + x116
    x118 = x117 * x16
    x119 = -x118
    x120 = x112 * x6 + x119
    x121 = x30 * x46
    x122 = x121 + x41
    x123 = x0 * x71
    x124 = x115 - x123
    x125 = x122 + x124
    x126 = x125 * x16
    x127 = -x126
    x128 = x117 * x6 + x127
    x129 = -x128 * x16
    x130 = x120 * x46 + x129
    x131 = x0 * x79
    x132 = x39 * x46
    x133 = x132 + x51
    x134 = x16 * (x123 - x131 + x133)
    x135 = -x134
    x136 = -x16 * (x125 * x6 + x135)
    x137 = x128 * x46 + x136
    x138 = 4 * x12
    x139 = x138 * x16
    x140 = x138 * x46
    x141 = x0 * (-x139 * x22 + x140 * x15)
    x142 = x46 * x89 + x93
    x143 = x0 * (4 * x13 * x46 - x139 * x15) - x141 + x142
    x144 = x0 * (-x139 * x28 + x140 * x22)
    x145 = x46 * x92 + x98
    x146 = x141 - x144 + x145
    x147 = -x146 * x16
    x148 = x143 * x46 + x147
    x149 = x104 + x46 * x97
    x150 = -x16 * (-x0 * (-x139 * x37 + x140 * x28) + x144 + x149)
    x151 = x146 * x46 + x150
    x152 = -x151 * x16
    x153 = x148 * x46 + x152
    x154 = x0 * x143 - x0 * x146 + x153
    x155 = x0 * (x20 + x89)
    x156 = x0 * (x24 + x92)
    x157 = x112 * x46
    x158 = x119 + x157
    x159 = x155 - x156 + x158
    x160 = x0 * (x30 + x97)
    x161 = x117 * x46
    x162 = x127 + x161
    x163 = x156 - x160 + x162
    x164 = -x16 * x163
    x165 = x159 * x46 + x164
    x166 = x0 * (x103 + x39)
    x167 = x125 * x46
    x168 = x135 + x167
    x169 = -x16 * (x160 - x166 + x168)
    x170 = x163 * x46 + x169
    x171 = -x16 * x170
    x172 = x165 * x46 + x171
    x173 = x0 * x159 - x0 * x163 + x172
    x174 = 2 * x110
    x175 = x0 * (2 * x107 + 2 * x109 - x174 - 2 * x25)
    x176 = 2 * x115
    x177 = x0 * (2 * x113 + x174 - x176 - 2 * x31)
    x178 = x27 * x46
    x179 = x178 + x35
    x180 = x175 - x177 + x179
    x181 = 2 * x123
    x182 = x0 * (2 * x121 + x176 - x181 - 2 * x40)
    x183 = x33 * x46
    x184 = x183 + x44
    x185 = x177 - x182 + x184
    x186 = x180 * x46
    x187 = x16 * x185
    x188 = -x187
    x189 = x186 + x188
    x190 = x185 * x46
    x191 = x0 * (-2 * x131 + 2 * x132 + x181 - 2 * x50)
    x192 = x42 * x46
    x193 = x192 + x53
    x194 = x16 * (x182 - x191 + x193)
    x195 = -x194
    x196 = x190 + x195
    x197 = -x16 * x196
    x198 = x189 * x46 + x197
    x199 = x0 * x180 - x0 * x185 + x198
    x200 = x36 * x6 + x47
    x201 = -x16 * (x45 * x6 + x54)
    x202 = x0 * x48 - x0 * x55 + x200 * x46 + x201
    x203 = x6 * x70 + x77
    x204 = -x16 * (x6 * x76 + x83)
    x205 = x100 + x6 * x94
    x206 = -x16 * (x105 + x6 * x99)
    x207 = x120 * x6 + x129
    x208 = -x16 * (x128 * x6 + x136)
    x209 = x143 * x6 + x147
    x210 = -x16 * (x146 * x6 + x150)
    x211 = x159 * x6 + x164
    x212 = -x16 * (x163 * x6 + x169)
    x213 = x180 * x6 + x188
    x214 = -x16 * (x185 * x6 + x195)
    x215 = 3 * x177
    x216 = 3 * x182
    x217 = x0 * (3 * x183 + x215 - x216 - 3 * x43)
    x218 = x0 * (3 * x175 + 3 * x178 - x215 - 3 * x34) - x217 + x48
    x219 = -x16 * (-x0 * (-3 * x191 + 3 * x192 + x216 - 3 * x52) + x217 + x55)
    x220 = x218 * x46 + x219
    x221 = x0 * x88
    x222 = x0 * x91
    x223 = 2 * x222
    x224 = x46 * x64
    x225 = x0 * x96
    x226 = 2 * x225
    x227 = x46 * x67
    x228 = x0 * (x223 - x226 + 2 * x227 - 2 * x74)
    x229 = x0 * (2 * x221 - x223 + 2 * x224 - 2 * x68) - x228 + x78
    x230 = x0 * x102
    x231 = x46 * x73
    x232 = -x16 * (-x0 * (x226 - 2 * x230 + 2 * x231 - 2 * x81) + x228 + x84)
    x233 = x229 * x46 + x232
    x234 = x0 * x108
    x235 = x0 * x114
    x236 = 2 * x235
    x237 = x0 * x122
    x238 = 2 * x237
    x239 = x0 * (2 * x183 + x236 - x238 - 2 * x43)
    x240 = x0 * (2 * x178 + 2 * x234 - x236 - 2 * x34) - x239 + x48
    x241 = x0 * x133
    x242 = -x16 * (-x0 * (2 * x192 + x238 - 2 * x241 - 2 * x52) + x239 + x55)
    x243 = x240 * x46 + x242
    x244 = x0 * x145
    x245 = x0 * x142 + x101 - x244
    x246 = -x16 * (-x0 * x149 + x106 + x244)
    x247 = x245 * x46 + x246
    x248 = x0 * x162
    x249 = x0 * x158 + x130 - x248
    x250 = -x16 * (-x0 * x168 + x137 + x248)
    x251 = x249 * x46 + x250
    x252 = x0 * x184
    x253 = x0 * x179 - x252 + x48
    x254 = -x16 * (-x0 * x193 + x252 + x55)
    x255 = x253 * x46 + x254
    x256 = x200 * x6 + x201
    x257 = 2 * x160
    x258 = 2 * x156
    x259 = 3 * x0
    x260 = x259 * (-2 * x126 + 2 * x161 - x257 + x258 + x33)
    x261 = x85 * (x116 - x16 * x96 + x46 * x91 + x67)
    x262 = x222 - x225 + x227 + x75
    x263 = x0 * x24
    x264 = x0 * x30
    x265 = x85 * (x114 * x46 - x122 * x16 + x263 - x264 + x33)
    x266 = x184 + x235 - x237
    x267 = x0 * x92
    x268 = x0 * x117
    x269 = x0 * x33

    # 60 item(s)
    S = numpy.array(
        [
            x58,
            x0 * x70 - x0 * x76 - x16 * x84 + x46 * x78,
            x58,
            x0 * x94 - x0 * x99 + x101 * x46 - x106 * x16,
            x0 * x120 - x0 * x128 + x130 * x46 - x137 * x16,
            x58,
            x154,
            x173,
            x199,
            x58,
            x202,
            x0 * x78 - x0 * x84 + x203 * x46 + x204,
            x202,
            x0 * x101 - x0 * x106 + x205 * x46 + x206,
            x0 * x130 - x0 * x137 + x207 * x46 + x208,
            x202,
            x0 * x148 - x0 * x151 + x209 * x46 + x210,
            x0 * x165 - x0 * x170 + x211 * x46 + x212,
            x0 * x189 - x0 * x196 + x213 * x46 + x214,
            x202,
            x220,
            x233,
            x243,
            x247,
            x251,
            x255,
            x153,
            x172,
            x198,
            x57,
            x256,
            x203 * x6 + x204,
            x256,
            x205 * x6 + x206,
            x207 * x6 + x208,
            x256,
            x209 * x6 + x210,
            x211 * x6 + x212,
            x213 * x6 + x214,
            x256,
            x218 * x6 + x219,
            x229 * x6 + x232,
            x240 * x6 + x242,
            x245 * x6 + x246,
            x249 * x6 + x250,
            x253 * x6 + x254,
            x148 * x6 + x152,
            x165 * x6 + x171,
            x189 * x6 + x197,
            x48 * x6 + x56,
            x0
            * (
                3 * x186
                - 3 * x187
                + x259 * (-2 * x118 + 2 * x155 + 2 * x157 - x258 + x27)
                - x260
                + x36
            )
            - x0
            * (
                3 * x190
                - 3 * x194
                - x259 * (-2 * x134 - 2 * x166 + 2 * x167 + x257 + x42)
                + x260
                + x45
            )
            + x220,
            x0
            * (
                -x17 * x262
                - x261
                + x59 * (x221 - x222 + x224 + x69)
                + x70
                + x85 * (x111 - x16 * x91 + x46 * x88 + x64)
            )
            - x0
            * (
                -x17 * (x225 - x230 + x231 + x82)
                + x261
                + x262 * x59
                + x76
                - x85 * (-x102 * x16 + x124 + x46 * x96 + x73)
            )
            + x233,
            x0
            * (
                -x17 * x266
                - x265
                + x36
                + x59 * (x179 + x234 - x235)
                + x85 * (x0 * x20 + x108 * x46 - x114 * x16 - x263 + x27)
            )
            - x0
            * (
                -x17 * (x193 + x237 - x241)
                + x265
                + x266 * x59
                + x45
                - x85 * (-x0 * x39 + x122 * x46 - x133 * x16 + x264 + x42)
            )
            + x243,
            x0 * (x0 * x89 + x142 * x46 - x145 * x16 - x267 + x94)
            - x0 * (-x0 * x97 + x145 * x46 - x149 * x16 + x267 + x99)
            + x247,
            x0 * (x0 * x112 + x120 + x158 * x46 - x16 * x162 - x268)
            - x0 * (-x0 * x125 + x128 - x16 * x168 + x162 * x46 + x268)
            + x251,
            x0 * (x0 * x27 - x16 * x184 + x179 * x46 - x269 + x36)
            - x0 * (-x0 * x42 - x16 * x193 + x184 * x46 + x269 + x45)
            + x255,
            x154,
            x173,
            x199,
            x58,
        ]
    )
    return S


def coulomb_24(a, A, b, B, C):
    """Cartesian (dg) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = -x2 * (a * A[1] + b * B[1])
    x5 = -x2 * (a * A[2] + b * B[2])
    x6 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x4 + B[1]) ** 2 + abs(x5 + B[2]) ** 2)
    x7 = x3 + C[0]
    x8 = x4 + C[1]
    x9 = x5 + C[2]
    x10 = x7 ** 2 + x8 ** 2 + x9 ** 2
    x11 = x1 * x10
    x12 = numpy.pi * x2 * numpy.exp(-a * b * x10 * x2)
    x13 = x12 * boys(0, x11)
    x14 = 2 * x6
    x15 = boys(1, x11)
    x16 = numpy.sqrt(abs(x7) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x17 = 2 * x16
    x18 = x12 * x17
    x19 = -x15 * x18
    x20 = x13 * x14 + x19
    x21 = x12 * x14
    x22 = boys(2, x11)
    x23 = -x18 * x22
    x24 = x15 * x21 + x23
    x25 = x16 * x24
    x26 = -x25
    x27 = x20 * x6 + x26
    x28 = boys(3, x11)
    x29 = -x18 * x28
    x30 = x21 * x22 + x29
    x31 = x16 * x30
    x32 = -x31
    x33 = x24 * x6 + x32
    x34 = x16 * x33
    x35 = -x34
    x36 = x27 * x6 + x35
    x37 = boys(4, x11)
    x38 = -x18 * x37
    x39 = x21 * x28 + x38
    x40 = x16 * x39
    x41 = -x40
    x42 = x30 * x6 + x41
    x43 = x16 * x42
    x44 = -x43
    x45 = x33 * x6 + x44
    x46 = x16 * x45
    x47 = -x46
    x48 = x36 * x6 + x47
    x49 = boys(5, x11)
    x50 = -x18 * x49
    x51 = x21 * x37 + x50
    x52 = x16 * x51
    x53 = -x52
    x54 = x39 * x6 + x53
    x55 = x16 * x54
    x56 = -x55
    x57 = x42 * x6 + x56
    x58 = x16 * x57
    x59 = -x58
    x60 = x45 * x6 + x59
    x61 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x4 + A[1]) ** 2 + abs(x5 + A[2]) ** 2)
    x62 = -x16 * x60
    x63 = x48 * x61 + x62
    x64 = -x18 * boys(6, x11)
    x65 = x16 * (x21 * x49 + x64)
    x66 = -x65
    x67 = x16 * (x51 * x6 + x66)
    x68 = -x67
    x69 = x16 * (x54 * x6 + x68)
    x70 = -x69
    x71 = -x16 * (x57 * x6 + x70)
    x72 = x60 * x61 + x71
    x73 = -x16 * x72
    x74 = x61 * x63 + x73
    x75 = x0 * x48 - x0 * x60 + x74
    x76 = 2 * x61
    x77 = x13 * x76 + x19
    x78 = x12 * x76
    x79 = x15 * x78 + x23
    x80 = x16 * x79
    x81 = -x80
    x82 = x6 * x77 + x81
    x83 = x22 * x78 + x29
    x84 = x16 * x83
    x85 = -x84
    x86 = x6 * x79 + x85
    x87 = x16 * x86
    x88 = -x87
    x89 = x6 * x82 + x88
    x90 = x28 * x78 + x38
    x91 = x16 * x90
    x92 = -x91
    x93 = x6 * x83 + x92
    x94 = x16 * x93
    x95 = -x94
    x96 = x6 * x86 + x95
    x97 = x16 * x96
    x98 = -x97
    x99 = x6 * x89 + x98
    x100 = x37 * x78 + x50
    x101 = x100 * x16
    x102 = -x101
    x103 = x102 + x6 * x90
    x104 = x103 * x16
    x105 = -x104
    x106 = x105 + x6 * x93
    x107 = x106 * x16
    x108 = -x107
    x109 = x108 + x6 * x96
    x110 = -x109 * x16
    x111 = x110 + x61 * x99
    x112 = x49 * x78 + x64
    x113 = -x112 * x16
    x114 = x16 * (x100 * x6 + x113)
    x115 = -x114
    x116 = x16 * (x103 * x6 + x115)
    x117 = -x116
    x118 = -x16 * (x106 * x6 + x117)
    x119 = x109 * x61 + x118
    x120 = 2 * x0
    x121 = x12 * x120
    x122 = x121 * x15
    x123 = x61 * x77
    x124 = x123 + x81
    x125 = x120 * x13 - x122 + x124
    x126 = x121 * x22
    x127 = x61 * x79
    x128 = x127 + x85
    x129 = x122 - x126 + x128
    x130 = -x129 * x16
    x131 = x125 * x6 + x130
    x132 = x121 * x28
    x133 = x61 * x83
    x134 = x133 + x92
    x135 = x126 - x132 + x134
    x136 = -x135 * x16
    x137 = x129 * x6 + x136
    x138 = x137 * x16
    x139 = -x138
    x140 = x131 * x6 + x139
    x141 = x121 * x37
    x142 = x61 * x90
    x143 = x102 + x142
    x144 = x132 - x141 + x143
    x145 = -x144 * x16
    x146 = x135 * x6 + x145
    x147 = x146 * x16
    x148 = -x147
    x149 = x137 * x6 + x148
    x150 = -x149 * x16
    x151 = x140 * x61 + x150
    x152 = x100 * x61 + x113
    x153 = -x121 * x49 + x141 + x152
    x154 = -x153 * x16
    x155 = x16 * (x144 * x6 + x154)
    x156 = -x155
    x157 = -x16 * (x146 * x6 + x156)
    x158 = x149 * x61 + x157
    x159 = x20 * x61
    x160 = x159 + x26
    x161 = x0 * x77
    x162 = x0 * x79
    x163 = x161 - x162
    x164 = x160 + x163
    x165 = x24 * x61
    x166 = x165 + x32
    x167 = x0 * x83
    x168 = x162 - x167
    x169 = x166 + x168
    x170 = x16 * x169
    x171 = -x170
    x172 = x164 * x6 + x171
    x173 = x30 * x61
    x174 = x173 + x41
    x175 = x0 * x90
    x176 = x167 - x175
    x177 = x174 + x176
    x178 = x16 * x177
    x179 = -x178
    x180 = x169 * x6 + x179
    x181 = x16 * x180
    x182 = -x181
    x183 = x172 * x6 + x182
    x184 = x39 * x61
    x185 = x184 + x53
    x186 = x0 * x100
    x187 = x175 - x186
    x188 = x185 + x187
    x189 = x16 * x188
    x190 = -x189
    x191 = x177 * x6 + x190
    x192 = x16 * x191
    x193 = -x192
    x194 = x180 * x6 + x193
    x195 = -x16 * x194
    x196 = x183 * x61 + x195
    x197 = x0 * x112
    x198 = x51 * x61
    x199 = x198 + x66
    x200 = x16 * (x186 - x197 + x199)
    x201 = -x200
    x202 = x16 * (x188 * x6 + x201)
    x203 = -x202
    x204 = -x16 * (x191 * x6 + x203)
    x205 = x194 * x61 + x204
    x206 = 4 * x12
    x207 = x16 * x206
    x208 = x206 * x61
    x209 = x0 * (x15 * x208 - x207 * x22)
    x210 = x125 * x61 + x130
    x211 = x0 * (4 * x13 * x61 - x15 * x207) - x209 + x210
    x212 = x0 * (-x207 * x28 + x208 * x22)
    x213 = x129 * x61 + x136
    x214 = x209 - x212 + x213
    x215 = -x16 * x214
    x216 = x211 * x6 + x215
    x217 = x0 * (-x207 * x37 + x208 * x28)
    x218 = x135 * x61 + x145
    x219 = x212 - x217 + x218
    x220 = -x16 * x219
    x221 = x214 * x6 + x220
    x222 = -x16 * x221
    x223 = x216 * x61 + x222
    x224 = x144 * x61 + x154
    x225 = -x0 * (-x207 * x49 + x208 * x37) + x217 + x224
    x226 = -x16 * x225
    x227 = -x16 * (x219 * x6 + x226)
    x228 = x221 * x61 + x227
    x229 = x0 * (x125 + x20)
    x230 = x0 * (x129 + x24)
    x231 = x164 * x61
    x232 = x171 + x231
    x233 = x229 - x230 + x232
    x234 = x0 * (x135 + x30)
    x235 = x169 * x61
    x236 = x179 + x235
    x237 = x230 - x234 + x236
    x238 = -x16 * x237
    x239 = x233 * x6 + x238
    x240 = x0 * (x144 + x39)
    x241 = x177 * x61
    x242 = x190 + x241
    x243 = x234 - x240 + x242
    x244 = -x16 * x243
    x245 = x237 * x6 + x244
    x246 = -x16 * x245
    x247 = x239 * x61 + x246
    x248 = x0 * (x153 + x51)
    x249 = x188 * x61
    x250 = x201 + x249
    x251 = -x16 * (x240 - x248 + x250)
    x252 = -x16 * (x243 * x6 + x251)
    x253 = x245 * x61 + x252
    x254 = 2 * x162
    x255 = 2 * x159 + 2 * x161 - 2 * x25 - x254
    x256 = x0 * x255
    x257 = 2 * x167
    x258 = 2 * x165 + x254 - x257 - 2 * x31
    x259 = x0 * x258
    x260 = x27 * x61
    x261 = x260 + x35
    x262 = x256 - x259 + x261
    x263 = 2 * x175
    x264 = 2 * x173 + x257 - x263 - 2 * x40
    x265 = x0 * x264
    x266 = x33 * x61
    x267 = x266 + x44
    x268 = x259 - x265 + x267
    x269 = x16 * x268
    x270 = -x269
    x271 = x262 * x6 + x270
    x272 = 2 * x186
    x273 = 2 * x184 + x263 - x272 - 2 * x52
    x274 = x0 * x273
    x275 = x42 * x61
    x276 = x275 + x56
    x277 = x265 - x274 + x276
    x278 = x16 * x277
    x279 = -x278
    x280 = x268 * x6 + x279
    x281 = -x16 * x280
    x282 = x271 * x61 + x281
    x283 = x0 * (-2 * x197 + 2 * x198 + x272 - 2 * x65)
    x284 = x54 * x61
    x285 = x284 + x68
    x286 = x16 * (x274 - x283 + x285)
    x287 = -x286
    x288 = -x16 * (x277 * x6 + x287)
    x289 = x280 * x61 + x288
    x290 = 6 * x0
    x291 = x12 * x290
    x292 = x15 * x291
    x293 = x22 * x291
    x294 = x0 * (3 * x127 + x292 - x293 - 3 * x84)
    x295 = x211 * x61 + x215
    x296 = x0 * (3 * x123 + x13 * x290 - x292 - 3 * x80) - x294 + x295
    x297 = x28 * x291
    x298 = x0 * (3 * x133 + x293 - x297 - 3 * x91)
    x299 = x214 * x61 + x220
    x300 = x294 - x298 + x299
    x301 = -x16 * x300
    x302 = x296 * x61 + x301
    x303 = x219 * x61 + x226
    x304 = -x16 * (-x0 * (-3 * x101 + 3 * x142 - x291 * x37 + x297) + x298 + x303)
    x305 = x300 * x61 + x304
    x306 = -x16 * x305
    x307 = x302 * x61 + x306
    x308 = x0 * x296 - x0 * x300 + x307
    x309 = x0 * (x214 + x258)
    x310 = x233 * x61 + x238
    x311 = x0 * (x211 + x255) - x309 + x310
    x312 = x0 * (x219 + x264)
    x313 = x237 * x61 + x244
    x314 = x309 - x312 + x313
    x315 = -x16 * x314
    x316 = x311 * x61 + x315
    x317 = x243 * x61 + x251
    x318 = -x16 * (-x0 * (x225 + x273) + x312 + x317)
    x319 = x314 * x61 + x318
    x320 = -x16 * x319
    x321 = x316 * x61 + x320
    x322 = x0 * x311 - x0 * x314 + x321
    x323 = 2 * x230
    x324 = x0 * (-2 * x170 + 2 * x229 + 2 * x231 + x27 - x323)
    x325 = 2 * x234
    x326 = x0 * (-2 * x178 + 2 * x235 + x323 - x325 + x33)
    x327 = x262 * x61
    x328 = x270 + x327
    x329 = x324 - x326 + x328
    x330 = 2 * x240
    x331 = x0 * (-2 * x189 + 2 * x241 + x325 - x330 + x42)
    x332 = x268 * x61
    x333 = x279 + x332
    x334 = x326 - x331 + x333
    x335 = -x16 * x334
    x336 = x329 * x61 + x335
    x337 = x0 * (-2 * x200 - 2 * x248 + 2 * x249 + x330 + x54)
    x338 = x277 * x61
    x339 = x287 + x338
    x340 = -x16 * (x331 - x337 + x339)
    x341 = x334 * x61 + x340
    x342 = -x16 * x341
    x343 = x336 * x61 + x342
    x344 = x0 * x329 - x0 * x334 + x343
    x345 = 3 * x259
    x346 = x0 * (3 * x256 + 3 * x260 - 3 * x34 - x345)
    x347 = 3 * x265
    x348 = x0 * (3 * x266 + x345 - x347 - 3 * x43)
    x349 = x36 * x61
    x350 = x349 + x47
    x351 = x346 - x348 + x350
    x352 = 3 * x274
    x353 = x0 * (3 * x275 + x347 - x352 - 3 * x55)
    x354 = x45 * x61
    x355 = x354 + x59
    x356 = x348 - x353 + x355
    x357 = x351 * x61
    x358 = x16 * x356
    x359 = -x358
    x360 = x357 + x359
    x361 = x356 * x61
    x362 = x0 * (-3 * x283 + 3 * x284 + x352 - 3 * x67)
    x363 = x57 * x61
    x364 = x363 + x70
    x365 = x16 * (x353 - x362 + x364)
    x366 = -x365
    x367 = x361 + x366
    x368 = -x16 * x367
    x369 = x360 * x61 + x368
    x370 = x0 * x351 - x0 * x356 + x369
    x371 = x48 * x6 + x62
    x372 = -x16 * (x6 * x60 + x71)
    x373 = x0 * x63 - x0 * x72 + x371 * x61 + x372
    x374 = x110 + x6 * x99
    x375 = -x16 * (x109 * x6 + x118)
    x376 = x140 * x6 + x150
    x377 = -x16 * (x149 * x6 + x157)
    x378 = x183 * x6 + x195
    x379 = -x16 * (x194 * x6 + x204)
    x380 = x216 * x6 + x222
    x381 = -x16 * (x221 * x6 + x227)
    x382 = x239 * x6 + x246
    x383 = -x16 * (x245 * x6 + x252)
    x384 = x271 * x6 + x281
    x385 = -x16 * (x280 * x6 + x288)
    x386 = x296 * x6 + x301
    x387 = -x16 * (x300 * x6 + x304)
    x388 = x311 * x6 + x315
    x389 = -x16 * (x314 * x6 + x318)
    x390 = x329 * x6 + x335
    x391 = -x16 * (x334 * x6 + x340)
    x392 = x351 * x6 + x359
    x393 = -x16 * (x356 * x6 + x366)
    x394 = 4 * x348
    x395 = 4 * x353
    x396 = x0 * (4 * x354 + x394 - x395 - 4 * x58)
    x397 = x0 * (4 * x346 + 4 * x349 - x394 - 4 * x46) - x396 + x63
    x398 = -x16 * (-x0 * (-4 * x362 + 4 * x363 + x395 - 4 * x69) + x396 + x72)
    x399 = x397 * x61 + x398
    x400 = x61 * x89
    x401 = x0 * x124
    x402 = x0 * x128
    x403 = 2 * x402
    x404 = x61 * x82
    x405 = x0 * (2 * x401 - x403 + 2 * x404 - 2 * x87)
    x406 = x0 * x134
    x407 = 2 * x406
    x408 = x61 * x86
    x409 = x0 * (x403 - x407 + 2 * x408 - 2 * x94)
    x410 = 3 * x409
    x411 = x61 * x96
    x412 = x0 * x143
    x413 = 2 * x412
    x414 = x61 * x93
    x415 = x0 * (-2 * x104 + x407 - x413 + 2 * x414)
    x416 = 3 * x415
    x417 = x0 * (-3 * x107 + x410 + 3 * x411 - x416)
    x418 = x0 * (3 * x400 + 3 * x405 - x410 - 3 * x97) + x111 - x417
    x419 = x106 * x61
    x420 = x0 * x152
    x421 = x103 * x61
    x422 = x0 * (-2 * x114 + x413 - 2 * x420 + 2 * x421)
    x423 = -x16 * (-x0 * (-3 * x116 + x416 + 3 * x419 - 3 * x422) + x119 + x417)
    x424 = x418 * x61 + x423
    x425 = x0 * x160
    x426 = x0 * x166
    x427 = 2 * x426
    x428 = x0 * (2 * x260 - 2 * x34 + 2 * x425 - x427)
    x429 = x0 * x174
    x430 = 2 * x429
    x431 = x0 * (2 * x266 + x427 - 2 * x43 - x430)
    x432 = 3 * x431
    x433 = x0 * x185
    x434 = 2 * x433
    x435 = x0 * (2 * x275 + x430 - x434 - 2 * x55)
    x436 = 3 * x435
    x437 = x0 * (3 * x354 + x432 - x436 - 3 * x58)
    x438 = x0 * (3 * x349 + 3 * x428 - x432 - 3 * x46) - x437 + x63
    x439 = x0 * x199
    x440 = x0 * (2 * x284 + x434 - 2 * x439 - 2 * x67)
    x441 = -x16 * (-x0 * (3 * x363 + x436 - 3 * x440 - 3 * x69) + x437 + x72)
    x442 = x438 * x61 + x441
    x443 = x0 * x210
    x444 = x0 * x213
    x445 = 2 * x444
    x446 = x131 * x61
    x447 = x0 * x218
    x448 = 2 * x447
    x449 = x137 * x61
    x450 = x0 * (-2 * x147 + x445 - x448 + 2 * x449)
    x451 = x0 * (-2 * x138 + 2 * x443 - x445 + 2 * x446) + x151 - x450
    x452 = x0 * x224
    x453 = x146 * x61
    x454 = -x16 * (-x0 * (-2 * x155 + x448 - 2 * x452 + 2 * x453) + x158 + x450)
    x455 = x451 * x61 + x454
    x456 = x0 * x232
    x457 = x0 * x236
    x458 = 2 * x457
    x459 = x172 * x61
    x460 = x0 * x242
    x461 = 2 * x460
    x462 = x180 * x61
    x463 = x0 * (-2 * x192 + x458 - x461 + 2 * x462)
    x464 = x0 * (-2 * x181 + 2 * x456 - x458 + 2 * x459) + x196 - x463
    x465 = x0 * x250
    x466 = x191 * x61
    x467 = -x16 * (-x0 * (-2 * x202 + x461 - 2 * x465 + 2 * x466) + x205 + x463)
    x468 = x464 * x61 + x467
    x469 = x0 * x261
    x470 = x0 * x267
    x471 = 2 * x470
    x472 = x0 * x276
    x473 = 2 * x472
    x474 = x0 * (2 * x354 + x471 - x473 - 2 * x58)
    x475 = x0 * (2 * x349 - 2 * x46 + 2 * x469 - x471) - x474 + x63
    x476 = x0 * x285
    x477 = -x16 * (-x0 * (2 * x363 + x473 - 2 * x476 - 2 * x69) + x474 + x72)
    x478 = x475 * x61 + x477
    x479 = x0 * x299
    x480 = x0 * x295 + x223 - x479
    x481 = -x16 * (-x0 * x303 + x228 + x479)
    x482 = x480 * x61 + x481
    x483 = x0 * x313
    x484 = x0 * x310 + x247 - x483
    x485 = -x16 * (-x0 * x317 + x253 + x483)
    x486 = x484 * x61 + x485
    x487 = x0 * x333
    x488 = x0 * x328 + x282 - x487
    x489 = -x16 * (-x0 * x339 + x289 + x487)
    x490 = x488 * x61 + x489
    x491 = x0 * x355
    x492 = x0 * x350 - x491 + x63
    x493 = -x16 * (-x0 * x364 + x491 + x72)
    x494 = x492 * x61 + x493
    x495 = x371 * x6 + x372
    x496 = 3 * x331
    x497 = 3 * x326
    x498 = 4 * x0
    x499 = x498 * (-3 * x278 + 3 * x332 + x45 - x496 + x497)
    x500 = x120 * (x134 * x61 - x143 * x16 + x176 + x93)
    x501 = x105 + x406 - x412 + x414
    x502 = x120 * (x128 * x61 - x134 * x16 + x168 + x86)
    x503 = x402 - x406 + x408 + x95
    x504 = 3 * x0
    x505 = x504 * (-x17 * x501 - x500 + x502 + x503 * x76 + x96)
    x506 = x108 + x409 + x411 - x415
    x507 = 3 * x16
    x508 = 3 * x61
    x509 = x0 * x30
    x510 = x0 * x39
    x511 = x120 * (-x16 * x185 + x174 * x61 + x42 + x509 - x510)
    x512 = x276 + x429 - x433
    x513 = x0 * x24
    x514 = x120 * (-x16 * x174 + x166 * x61 + x33 - x509 + x513)
    x515 = x267 + x426 - x429
    x516 = x504 * (-x17 * x512 + x45 - x511 + x514 + x515 * x76)
    x517 = x355 + x431 - x435
    x518 = x0 * x129
    x519 = x0 * x135
    x520 = x120 * (x137 - x16 * x218 + x213 * x61 + x518 - x519)
    x521 = x148 + x444 - x447 + x449
    x522 = x0 * x169
    x523 = x0 * x177
    x524 = x120 * (-x16 * x242 + x180 + x236 * x61 + x522 - x523)
    x525 = x193 + x457 - x460 + x462
    x526 = x0 * x33
    x527 = x0 * x42
    x528 = x120 * (-x16 * x276 + x267 * x61 + x45 + x526 - x527)
    x529 = x355 + x470 - x472
    x530 = x0 * x214
    x531 = x0 * x237
    x532 = x0 * x268
    x533 = x0 * x45

    # 90 item(s)
    S = numpy.array(
        [
            x75,
            -x0 * x109 + x0 * x99 + x111 * x61 - x119 * x16,
            x75,
            x0 * x140 - x0 * x149 + x151 * x61 - x158 * x16,
            x0 * x183 - x0 * x194 - x16 * x205 + x196 * x61,
            x75,
            x0 * x216 - x0 * x221 - x16 * x228 + x223 * x61,
            x0 * x239 - x0 * x245 - x16 * x253 + x247 * x61,
            x0 * x271 - x0 * x280 - x16 * x289 + x282 * x61,
            x75,
            x308,
            x322,
            x344,
            x370,
            x75,
            x373,
            x0 * x111 - x0 * x119 + x374 * x61 + x375,
            x373,
            x0 * x151 - x0 * x158 + x376 * x61 + x377,
            x0 * x196 - x0 * x205 + x378 * x61 + x379,
            x373,
            x0 * x223 - x0 * x228 + x380 * x61 + x381,
            x0 * x247 - x0 * x253 + x382 * x61 + x383,
            x0 * x282 - x0 * x289 + x384 * x61 + x385,
            x373,
            x0 * x302 - x0 * x305 + x386 * x61 + x387,
            x0 * x316 - x0 * x319 + x388 * x61 + x389,
            x0 * x336 - x0 * x341 + x390 * x61 + x391,
            x0 * x360 - x0 * x367 + x392 * x61 + x393,
            x373,
            x399,
            x424,
            x442,
            x455,
            x468,
            x478,
            x482,
            x486,
            x490,
            x494,
            x307,
            x321,
            x343,
            x369,
            x74,
            x495,
            x374 * x6 + x375,
            x495,
            x376 * x6 + x377,
            x378 * x6 + x379,
            x495,
            x380 * x6 + x381,
            x382 * x6 + x383,
            x384 * x6 + x385,
            x495,
            x386 * x6 + x387,
            x388 * x6 + x389,
            x390 * x6 + x391,
            x392 * x6 + x393,
            x495,
            x397 * x6 + x398,
            x418 * x6 + x423,
            x438 * x6 + x441,
            x451 * x6 + x454,
            x464 * x6 + x467,
            x475 * x6 + x477,
            x480 * x6 + x481,
            x484 * x6 + x485,
            x488 * x6 + x489,
            x492 * x6 + x493,
            x302 * x6 + x306,
            x316 * x6 + x320,
            x336 * x6 + x342,
            x360 * x6 + x368,
            x6 * x63 + x73,
            x0
            * (
                4 * x357
                - 4 * x358
                + x48
                + x498 * (-3 * x269 + 3 * x324 + 3 * x327 + x36 - x497)
                - x499
            )
            - x0
            * (
                4 * x361
                - 4 * x365
                - x498 * (-3 * x286 - 3 * x337 + 3 * x338 + x496 + x57)
                + x499
                + x60
            )
            + x399,
            -x0
            * (
                x109
                - x504
                * (
                    x106
                    - x120 * (x103 + x143 * x61 - x152 * x16 + x187)
                    - x17 * (x115 + x412 - x420 + x421)
                    + x500
                    + x501 * x76
                )
                + x505
                + x506 * x508
                - x507 * (x117 + x415 + x419 - x422)
            )
            + x0
            * (
                x504
                * (
                    x120 * (x124 * x61 - x128 * x16 + x163 + x82)
                    - x17 * x503
                    - x502
                    + x76 * (x401 - x402 + x404 + x88)
                    + x89
                )
                - x505
                - x506 * x507
                + x508 * (x400 + x405 - x409 + x98)
                + x99
            )
            + x424,
            x0
            * (
                x48
                + x504
                * (
                    x120 * (x0 * x20 - x16 * x166 + x160 * x61 + x27 - x513)
                    - x17 * x515
                    + x36
                    - x514
                    + x76 * (x261 + x425 - x426)
                )
                - x507 * x517
                + x508 * (x350 + x428 - x431)
                - x516
            )
            - x0
            * (
                -x504
                * (
                    -x120 * (-x0 * x51 - x16 * x199 + x185 * x61 + x510 + x54)
                    - x17 * (x285 + x433 - x439)
                    + x511
                    + x512 * x76
                    + x57
                )
                - x507 * (x364 + x435 - x440)
                + x508 * x517
                + x516
                + x60
            )
            + x442,
            x0
            * (
                x120 * (x0 * x125 + x131 - x16 * x213 + x210 * x61 - x518)
                + x140
                - x17 * x521
                - x520
                + x76 * (x139 + x443 - x444 + x446)
            )
            - x0
            * (
                -x120 * (-x0 * x144 + x146 - x16 * x224 + x218 * x61 + x519)
                + x149
                - x17 * (x156 + x447 - x452 + x453)
                + x520
                + x521 * x76
            )
            + x455,
            x0
            * (
                x120 * (x0 * x164 - x16 * x236 + x172 + x232 * x61 - x522)
                - x17 * x525
                + x183
                - x524
                + x76 * (x182 + x456 - x457 + x459)
            )
            - x0
            * (
                -x120 * (-x0 * x188 - x16 * x250 + x191 + x242 * x61 + x523)
                - x17 * (x203 + x460 - x465 + x466)
                + x194
                + x524
                + x525 * x76
            )
            + x468,
            x0
            * (
                x120 * (x0 * x27 - x16 * x267 + x261 * x61 + x36 - x526)
                - x17 * x529
                + x48
                - x528
                + x76 * (x350 + x469 - x470)
            )
            - x0
            * (
                -x120 * (-x0 * x54 - x16 * x285 + x276 * x61 + x527 + x57)
                - x17 * (x364 + x472 - x476)
                + x528
                + x529 * x76
                + x60
            )
            + x478,
            x0 * (x0 * x211 - x16 * x299 + x216 + x295 * x61 - x530)
            - x0 * (-x0 * x219 - x16 * x303 + x221 + x299 * x61 + x530)
            + x482,
            x0 * (x0 * x233 - x16 * x313 + x239 + x310 * x61 - x531)
            - x0 * (-x0 * x243 - x16 * x317 + x245 + x313 * x61 + x531)
            + x486,
            x0 * (x0 * x262 - x16 * x333 + x271 + x328 * x61 - x532)
            - x0 * (-x0 * x277 - x16 * x339 + x280 + x333 * x61 + x532)
            + x490,
            x0 * (x0 * x36 - x16 * x355 + x350 * x61 + x48 - x533)
            - x0 * (-x0 * x57 - x16 * x364 + x355 * x61 + x533 + x60)
            + x494,
            x308,
            x322,
            x344,
            x370,
            x75,
        ]
    )
    return S


def coulomb_30(a, A, b, B, C):
    """Cartesian (fs) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = x3 + C[0]
    x5 = -x2 * (a * A[1] + b * B[1])
    x6 = x5 + C[1]
    x7 = -x2 * (a * A[2] + b * B[2])
    x8 = x7 + C[2]
    x9 = x4 ** 2 + x6 ** 2 + x8 ** 2
    x10 = numpy.pi * x2 * numpy.exp(-a * b * x2 * x9)
    x11 = 4 * x10
    x12 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x5 + A[1]) ** 2 + abs(x7 + A[2]) ** 2)
    x13 = x1 * x9
    x14 = boys(0, x13)
    x15 = x12 * x14
    x16 = numpy.sqrt(abs(x4) ** 2 + abs(x6) ** 2 + abs(x8) ** 2)
    x17 = boys(1, x13)
    x18 = x16 * x17
    x19 = boys(2, x13)
    x20 = x16 * x19
    x21 = 2 * x10
    x22 = x17 * x21
    x23 = x0 * x22
    x24 = x0 * x21
    x25 = -x18 * x21
    x26 = x15 * x21 + x25
    x27 = -x20 * x21
    x28 = x12 * x22 + x27
    x29 = -x16 * x28
    x30 = x12 * x26 + x29
    x31 = x14 * x24 - x23 + x30
    x32 = x19 * x21
    x33 = -x16 * x21 * boys(3, x13)
    x34 = x12 * x32 + x33
    x35 = -x16 * x34
    x36 = x12 * x28 + x35
    x37 = -x19 * x24 + x23 + x36
    x38 = -x16 * x37
    x39 = x12 * x31 + x38
    x40 = x0 * (x11 * x15 - x11 * x18) - x0 * (x11 * x12 * x17 - x11 * x20) + x39
    x41 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x5 + B[1]) ** 2 + abs(x7 + B[2]) ** 2)
    x42 = x14 * x21 * x41 + x25
    x43 = x22 * x41 + x27
    x44 = x16 * x43
    x45 = -x44
    x46 = x12 * x42
    x47 = x0 * x26
    x48 = x0 * x28
    x49 = x47 - x48
    x50 = x0 * x34
    x51 = x12 * x43
    x52 = x16 * (x32 * x41 + x33)
    x53 = -x52
    x54 = x41 * x42 + x45
    x55 = -x16 * (x41 * x43 + x53)
    x56 = 2 * x48
    x57 = x26 * x41 + x29
    x58 = -x16 * (x28 * x41 + x35)

    # 10 item(s)
    S = numpy.array(
        [
            x40,
            x0 * (x31 + x42)
            - x0 * (x37 + x43)
            + x12 * (x45 + x46 + x49)
            - x16 * (x48 - x50 + x51 + x53),
            x12 * x30 - x16 * x36 + x49,
            x0 * (-2 * x44 + 2 * x46 + 2 * x47 - x56)
            - x0 * (-2 * x50 + 2 * x51 - 2 * x52 + x56)
            + x12 * x54
            + x55,
            x0 * x30 - x0 * x36 + x12 * x57 + x58,
            x39,
            x41 * x54 + x55,
            x41 * x57 + x58,
            x31 * x41 + x38,
            x40,
        ]
    )
    return S


def coulomb_31(a, A, b, B, C):
    """Cartesian (fp) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = -x2 * (a * A[1] + b * B[1])
    x5 = -x2 * (a * A[2] + b * B[2])
    x6 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x4 + A[1]) ** 2 + abs(x5 + A[2]) ** 2)
    x7 = x3 + C[0]
    x8 = x4 + C[1]
    x9 = x5 + C[2]
    x10 = x7 ** 2 + x8 ** 2 + x9 ** 2
    x11 = x1 * x10
    x12 = boys(0, x11)
    x13 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x4 + B[1]) ** 2 + abs(x5 + B[2]) ** 2)
    x14 = numpy.pi * x2 * numpy.exp(-a * b * x10 * x2)
    x15 = 2 * x14
    x16 = x13 * x15
    x17 = boys(1, x11)
    x18 = numpy.sqrt(abs(x7) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x19 = x15 * x18
    x20 = -x17 * x19
    x21 = x12 * x16 + x20
    x22 = x21 * x6
    x23 = boys(2, x11)
    x24 = -x19 * x23
    x25 = x16 * x17 + x24
    x26 = x18 * x25
    x27 = 2 * x22 - 2 * x26
    x28 = x25 * x6
    x29 = boys(3, x11)
    x30 = -x19 * x29
    x31 = x16 * x23 + x30
    x32 = x18 * x31
    x33 = 2 * x28 - 2 * x32
    x34 = x0 * x25
    x35 = -x26
    x36 = x22 + x35
    x37 = -x32
    x38 = x28 + x37
    x39 = -x18 * x38
    x40 = x36 * x6 + x39
    x41 = x0 * x21 - x34 + x40
    x42 = -x19 * boys(4, x11)
    x43 = -x18 * (x16 * x29 + x42)
    x44 = x31 * x6 + x43
    x45 = -x18 * x44
    x46 = x38 * x6 + x45
    x47 = -x0 * x31 + x34 + x46
    x48 = -x18 * x47
    x49 = x41 * x6 + x48
    x50 = x0 * x27 - x0 * x33 + x49
    x51 = x15 * x6
    x52 = x12 * x51 + x20
    x53 = x52 * x6
    x54 = x17 * x51 + x24
    x55 = x18 * x54
    x56 = x54 * x6
    x57 = x23 * x51 + x30
    x58 = x18 * x57
    x59 = -x55
    x60 = x53 + x59
    x61 = -x58
    x62 = x56 + x61
    x63 = -x18 * x62
    x64 = x6 * x60 + x63
    x65 = x0 * x52
    x66 = x0 * x54
    x67 = x65 - x66
    x68 = x64 + x67
    x69 = x29 * x51 + x42
    x70 = -x18 * x69
    x71 = x57 * x6 + x70
    x72 = -x18 * x71
    x73 = x6 * x62 + x72
    x74 = x0 * x57
    x75 = x66 - x74
    x76 = x73 + x75
    x77 = -x18 * x76
    x78 = x6 * x68 + x77
    x79 = x0 * (2 * x53 - 2 * x55) - x0 * (2 * x56 - 2 * x58) + x78
    x80 = x13 * x21 + x35
    x81 = x13 * x25 + x37
    x82 = x18 * x81
    x83 = -x82
    x84 = x6 * x80
    x85 = x0 * x36
    x86 = x0 * x38
    x87 = x85 - x86
    x88 = x0 * x44
    x89 = x6 * x81
    x90 = x18 * (x13 * x31 + x43)
    x91 = -x90
    x92 = (
        x0 * (x41 + x80)
        - x0 * (x47 + x81)
        - x18 * (x86 - x88 + x89 + x91)
        + x6 * (x83 + x84 + x87)
    )
    x93 = x13 * x52 + x59
    x94 = x13 * x54 + x61
    x95 = x18 * x94
    x96 = -x95
    x97 = x6 * x93
    x98 = x0 * x60
    x99 = x0 * x62
    x100 = x98 - x99
    x101 = x0 * x71
    x102 = x6 * x94
    x103 = x18 * (x13 * x57 + x70)
    x104 = -x103
    x105 = x36 + x67
    x106 = x38 + x75
    x107 = -x106 * x18
    x108 = x105 * x6 + x107
    x109 = -x18 * (-x0 * x69 + x44 + x74)
    x110 = x106 * x6 + x109
    x111 = x13 * x80 + x83
    x112 = -x18 * (x13 * x81 + x91)
    x113 = 2 * x86
    x114 = (
        x0 * (-x113 - 2 * x82 + 2 * x84 + 2 * x85)
        - x0 * (x113 - 2 * x88 + 2 * x89 - 2 * x90)
        + x111 * x6
        + x112
    )
    x115 = x13 * x93 + x96
    x116 = -x18 * (x104 + x13 * x94)
    x117 = 2 * x99
    x118 = x105 * x13 + x107
    x119 = -x18 * (x106 * x13 + x109)
    x120 = x13 * x60 + x63
    x121 = -x18 * (x13 * x62 + x72)
    x122 = x13 * x36 + x39
    x123 = -x18 * (x13 * x38 + x45)
    x124 = x0 * x15
    x125 = x124 * x17
    x126 = x12 * x124 - x125 + x60
    x127 = x124 * x23
    x128 = x125 - x127 + x62
    x129 = x0 * (x128 + x25)
    x130 = x0 * (x126 + x21) + x108 - x129
    x131 = -x124 * x29 + x127 + x71
    x132 = -x18 * (-x0 * (x131 + x31) + x110 + x129)
    x133 = x130 * x6 + x132
    x134 = x111 * x13 + x112
    x135 = 4 * x14
    x136 = x135 * x6
    x137 = x135 * x18
    x138 = x0 * (x136 * x17 - x137 * x23)
    x139 = 2 * x66

    # 30 item(s)
    S = numpy.array(
        [
            x50,
            x79,
            x50,
            x92,
            x0 * (x68 + x93)
            - x0 * (x76 + x94)
            - x18 * (-x101 + x102 + x104 + x99)
            + x6 * (x100 + x96 + x97),
            x92,
            x0 * x105 - x0 * x106 + x108 * x6 - x110 * x18,
            x100 - x18 * x73 + x6 * x64,
            -x18 * x46 + x40 * x6 + x87,
            x114,
            -x0 * (-2 * x101 + 2 * x102 - 2 * x103 + x117)
            + x0 * (-x117 - 2 * x95 + 2 * x97 + 2 * x98)
            + x115 * x6
            + x116,
            x114,
            x0 * x108 - x0 * x110 + x118 * x6 + x119,
            x0 * x64 - x0 * x73 + x120 * x6 + x121,
            x0 * x40 - x0 * x46 + x122 * x6 + x123,
            x133,
            x78,
            x49,
            x134,
            x115 * x13 + x116,
            x134,
            x118 * x13 + x119,
            x120 * x13 + x121,
            x122 * x13 + x123,
            x13 * x130 + x132,
            x13 * x68 + x77,
            x13 * x41 + x48,
            x0
            * (
                x0 * (x12 * x136 - x137 * x17)
                + x126 * x6
                - x128 * x18
                - x138
                - x139
                + x27
                + 2 * x65
            )
            - x0
            * (
                -x0 * (x136 * x23 - x137 * x29)
                + x128 * x6
                - x131 * x18
                + x138
                + x139
                + x33
                - 2 * x74
            )
            + x133,
            x79,
            x50,
        ]
    )
    return S


def coulomb_32(a, A, b, B, C):
    """Cartesian (fd) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = -x2 * (a * A[1] + b * B[1])
    x5 = -x2 * (a * A[2] + b * B[2])
    x6 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x4 + A[1]) ** 2 + abs(x5 + A[2]) ** 2)
    x7 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x4 + B[1]) ** 2 + abs(x5 + B[2]) ** 2)
    x8 = x3 + C[0]
    x9 = x4 + C[1]
    x10 = x5 + C[2]
    x11 = x10 ** 2 + x8 ** 2 + x9 ** 2
    x12 = x1 * x11
    x13 = boys(0, x12)
    x14 = numpy.pi * x2 * numpy.exp(-a * b * x11 * x2)
    x15 = 2 * x7
    x16 = x14 * x15
    x17 = boys(1, x12)
    x18 = numpy.sqrt(abs(x10) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x19 = 2 * x18
    x20 = x14 * x19
    x21 = -x17 * x20
    x22 = x13 * x16 + x21
    x23 = x14 * x17
    x24 = boys(2, x12)
    x25 = -x20 * x24
    x26 = x15 * x23 + x25
    x27 = x18 * x26
    x28 = -x27
    x29 = x22 * x7 + x28
    x30 = x29 * x6
    x31 = boys(3, x12)
    x32 = -x20 * x31
    x33 = x16 * x24 + x32
    x34 = x18 * x33
    x35 = -x34
    x36 = x26 * x7 + x35
    x37 = x18 * x36
    x38 = 2 * x30 - 2 * x37
    x39 = x36 * x6
    x40 = boys(4, x12)
    x41 = -x20 * x40
    x42 = x16 * x31 + x41
    x43 = x18 * x42
    x44 = -x43
    x45 = x33 * x7 + x44
    x46 = x18 * x45
    x47 = 2 * x39 - 2 * x46
    x48 = x0 * x36
    x49 = -x37
    x50 = x30 + x49
    x51 = -x46
    x52 = x39 + x51
    x53 = -x18 * x52
    x54 = x50 * x6 + x53
    x55 = x0 * x29 - x48 + x54
    x56 = -x20 * boys(5, x12)
    x57 = x18 * (x16 * x40 + x56)
    x58 = -x57
    x59 = -x18 * (x42 * x7 + x58)
    x60 = x45 * x6 + x59
    x61 = -x18 * x60
    x62 = x52 * x6 + x61
    x63 = -x0 * x45 + x48 + x62
    x64 = -x18 * x63
    x65 = x55 * x6 + x64
    x66 = x0 * x38 - x0 * x47 + x65
    x67 = 2 * x6
    x68 = x14 * x67
    x69 = x13 * x68 + x21
    x70 = x23 * x67 + x25
    x71 = x18 * x70
    x72 = -x71
    x73 = x69 * x7 + x72
    x74 = x6 * x73
    x75 = x24 * x68 + x32
    x76 = x18 * x75
    x77 = -x76
    x78 = x7 * x70 + x77
    x79 = x18 * x78
    x80 = 2 * x74 - 2 * x79
    x81 = x6 * x78
    x82 = x31 * x68 + x41
    x83 = x18 * x82
    x84 = -x83
    x85 = x7 * x75 + x84
    x86 = x18 * x85
    x87 = 2 * x81 - 2 * x86
    x88 = x0 * x78
    x89 = -x79
    x90 = x74 + x89
    x91 = -x86
    x92 = x81 + x91
    x93 = x0 * x73 - x18 * x92 + x6 * x90 - x88
    x94 = x40 * x68 + x56
    x95 = -x18 * x94
    x96 = -x18 * (x7 * x82 + x95)
    x97 = x6 * x85 + x96
    x98 = -x0 * x85 - x18 * x97 + x6 * x92 + x88
    x99 = 2 * x0
    x100 = x23 * x99
    x101 = x14 * x99
    x102 = x6 * x69
    x103 = x102 + x72
    x104 = -x100 + x101 * x13 + x103
    x105 = x104 * x6
    x106 = x101 * x24
    x107 = x6 * x70
    x108 = x107 + x77
    x109 = x100 - x106 + x108
    x110 = x109 * x18
    x111 = x109 * x6
    x112 = x101 * x31
    x113 = x6 * x75
    x114 = x113 + x84
    x115 = x106 - x112 + x114
    x116 = x115 * x18
    x117 = x0 * x109
    x118 = -x110
    x119 = x105 + x118
    x120 = -x116
    x121 = x111 + x120
    x122 = -x121 * x18
    x123 = x119 * x6 + x122
    x124 = x0 * x104 - x117 + x123
    x125 = x6 * x82 + x95
    x126 = -x101 * x40 + x112 + x125
    x127 = -x126 * x18
    x128 = x115 * x6 + x127
    x129 = -x128 * x18
    x130 = x121 * x6 + x129
    x131 = -x0 * x115 + x117 + x130
    x132 = -x131 * x18
    x133 = x124 * x6 + x132
    x134 = x0 * (2 * x105 - 2 * x110) - x0 * (2 * x111 - 2 * x116) + x133
    x135 = x22 * x6
    x136 = x135 + x28
    x137 = x0 * x69
    x138 = x0 * x70
    x139 = x137 - x138
    x140 = x136 + x139
    x141 = x140 * x6
    x142 = x26 * x6
    x143 = x142 + x35
    x144 = x0 * x75
    x145 = x138 - x144
    x146 = x143 + x145
    x147 = x146 * x18
    x148 = 2 * x141 - 2 * x147
    x149 = x146 * x6
    x150 = x33 * x6
    x151 = x150 + x44
    x152 = x0 * x82
    x153 = x144 - x152
    x154 = x151 + x153
    x155 = x154 * x18
    x156 = 2 * x149 - 2 * x155
    x157 = x0 * x146
    x158 = -x147
    x159 = x141 + x158
    x160 = -x155
    x161 = x149 + x160
    x162 = -x161 * x18
    x163 = x159 * x6 + x162
    x164 = x0 * x140 - x157 + x163
    x165 = x154 * x6
    x166 = x0 * x94
    x167 = x42 * x6
    x168 = x167 + x58
    x169 = x18 * (x152 - x166 + x168)
    x170 = -x169
    x171 = x165 + x170
    x172 = -x171 * x18
    x173 = x161 * x6 + x172
    x174 = -x0 * x154 + x157 + x173
    x175 = -x174 * x18
    x176 = x164 * x6 + x175
    x177 = x0 * x148 - x0 * x156 + x176
    x178 = x29 * x7 + x49
    x179 = x36 * x7 + x51
    x180 = x179 * x18
    x181 = -x180
    x182 = x178 * x6
    x183 = x0 * x50
    x184 = x0 * x52
    x185 = x183 - x184
    x186 = x0 * x60
    x187 = x179 * x6
    x188 = x18 * (x45 * x7 + x59)
    x189 = -x188
    x190 = (
        x0 * (x178 + x55)
        - x0 * (x179 + x63)
        - x18 * (x184 - x186 + x187 + x189)
        + x6 * (x181 + x182 + x185)
    )
    x191 = x7 * x73 + x89
    x192 = x7 * x78 + x91
    x193 = x0 * x90
    x194 = x0 * x92
    x195 = x191 * x6
    x196 = x18 * x192
    x197 = -x196
    x198 = x0 * x97
    x199 = x192 * x6
    x200 = x18 * (x7 * x85 + x96)
    x201 = -x200
    x202 = x104 * x7 + x118
    x203 = x109 * x7 + x120
    x204 = x18 * x203
    x205 = -x204
    x206 = x202 * x6
    x207 = x0 * x119
    x208 = x0 * x121
    x209 = x207 - x208
    x210 = x0 * x128
    x211 = x203 * x6
    x212 = x18 * (x115 * x7 + x127)
    x213 = -x212
    x214 = x140 * x7 + x158
    x215 = x146 * x7 + x160
    x216 = x18 * x215
    x217 = -x216
    x218 = x214 * x6
    x219 = x0 * x159
    x220 = x0 * x161
    x221 = x219 - x220
    x222 = x0 * x171
    x223 = x215 * x6
    x224 = x18 * (x154 * x7 + x170)
    x225 = -x224
    x226 = 2 * x138
    x227 = 2 * x135 - 2 * x27
    x228 = 2 * x137 - x226 + x227
    x229 = x0 * x228
    x230 = 2 * x144
    x231 = 2 * x142 - 2 * x34
    x232 = x226 - x230 + x231
    x233 = x0 * x232
    x234 = x229 - x233 + x50
    x235 = 2 * x152
    x236 = 2 * x150 - 2 * x43
    x237 = x230 - x235 + x236
    x238 = x0 * x237
    x239 = x233 - x238 + x52
    x240 = -x18 * x239
    x241 = x234 * x6 + x240
    x242 = -x18 * (-x0 * (-2 * x166 + 2 * x167 + x235 - 2 * x57) + x238 + x60)
    x243 = x239 * x6 + x242
    x244 = x0 * x103
    x245 = x0 * x108
    x246 = x244 - x245 + x90
    x247 = x0 * x114
    x248 = x245 - x247 + x92
    x249 = -x18 * x248
    x250 = x246 * x6 + x249
    x251 = -x18 * (-x0 * x125 + x247 + x97)
    x252 = x248 * x6 + x251
    x253 = x0 * x136
    x254 = x0 * x143
    x255 = x253 - x254 + x50
    x256 = x0 * x151
    x257 = x254 - x256 + x52
    x258 = -x18 * x257
    x259 = x255 * x6 + x258
    x260 = -x18 * (-x0 * x168 + x256 + x60)
    x261 = x257 * x6 + x260
    x262 = x178 * x7 + x181
    x263 = -x18 * (x179 * x7 + x189)
    x264 = 2 * x184
    x265 = (
        x0 * (-2 * x180 + 2 * x182 + 2 * x183 - x264)
        - x0 * (-2 * x186 + 2 * x187 - 2 * x188 + x264)
        + x262 * x6
        + x263
    )
    x266 = x191 * x7 + x197
    x267 = -x18 * (x192 * x7 + x201)
    x268 = 2 * x194
    x269 = x202 * x7 + x205
    x270 = -x18 * (x203 * x7 + x213)
    x271 = 2 * x208
    x272 = x214 * x7 + x217
    x273 = -x18 * (x215 * x7 + x225)
    x274 = 2 * x220
    x275 = x234 * x7 + x240
    x276 = -x18 * (x239 * x7 + x242)
    x277 = x246 * x7 + x249
    x278 = -x18 * (x248 * x7 + x251)
    x279 = x255 * x7 + x258
    x280 = -x18 * (x257 * x7 + x260)
    x281 = x119 * x7 + x122
    x282 = -x18 * (x121 * x7 + x129)
    x283 = x159 * x7 + x162
    x284 = -x18 * (x161 * x7 + x172)
    x285 = x50 * x7 + x53
    x286 = -x18 * (x52 * x7 + x61)
    x287 = x0 * (x109 + x26)
    x288 = 2 * x287
    x289 = x0 * (x104 + x22)
    x290 = x0 * (x115 + x33)
    x291 = 2 * x290
    x292 = x0 * (x156 + x288 - x291 + x36)
    x293 = x0 * (x148 - x288 + 2 * x289 + x29) + x241 - x292
    x294 = x0 * (x126 + x42)
    x295 = -x18 * (-x0 * (2 * x165 - 2 * x169 + x291 - 2 * x294 + x45) + x243 + x292)
    x296 = x293 * x6 + x295
    x297 = x103 * x6 - x108 * x18 + x139
    x298 = x108 * x6 - x114 * x18 + x145
    x299 = x0 * (x298 + x78)
    x300 = x0 * (x297 + x73) + x250 - x299
    x301 = x114 * x6 - x125 * x18 + x153
    x302 = -x18 * (-x0 * (x301 + x85) + x252 + x299)
    x303 = x300 * x6 + x302
    x304 = x0 * x26
    x305 = x0 * x22 + x136 * x6 - x143 * x18 - x304
    x306 = x0 * x33
    x307 = x143 * x6 - x151 * x18 + x304 - x306
    x308 = x0 * (x307 + x36)
    x309 = x0 * (x29 + x305) + x259 - x308
    x310 = -x0 * x42 + x151 * x6 - x168 * x18 + x306
    x311 = -x18 * (-x0 * (x310 + x45) + x261 + x308)
    x312 = x309 * x6 + x311
    x313 = x262 * x7 + x263
    x314 = 2 * x233
    x315 = 4 * x6
    x316 = 4 * x18
    x317 = x14 * x316
    x318 = x0 * (x23 * x315 - x24 * x317)
    x319 = x14 * x315
    x320 = x0 * (x24 * x319 - x31 * x317)
    x321 = x99 * (x121 + x232 + x318 - x320)
    x322 = x161 + x287 - x290
    x323 = x0 * (2 * x107 - 2 * x76)
    x324 = 2 * x245
    x325 = x0 * x231
    x326 = 2 * x254

    # 60 item(s)
    S = numpy.array(
        [
            x66,
            x0 * x80 - x0 * x87 - x18 * x98 + x6 * x93,
            x66,
            x134,
            x177,
            x66,
            x190,
            x0 * (x191 + x93)
            - x0 * (x192 + x98)
            - x18 * (x194 - x198 + x199 + x201)
            + x6 * (x193 - x194 + x195 + x197),
            x190,
            x0 * (x124 + x202)
            - x0 * (x131 + x203)
            - x18 * (x208 - x210 + x211 + x213)
            + x6 * (x205 + x206 + x209),
            x0 * (x164 + x214)
            - x0 * (x174 + x215)
            - x18 * (x220 - x222 + x223 + x225)
            + x6 * (x217 + x218 + x221),
            x190,
            x0 * x234 - x0 * x239 - x18 * x243 + x241 * x6,
            x0 * x246 - x0 * x248 - x18 * x252 + x250 * x6,
            x0 * x255 - x0 * x257 - x18 * x261 + x259 * x6,
            x123 * x6 - x130 * x18 + x209,
            x163 * x6 - x173 * x18 + x221,
            -x18 * x62 + x185 + x54 * x6,
            x265,
            x0 * (2 * x193 + 2 * x195 - 2 * x196 - x268)
            - x0 * (-2 * x198 + 2 * x199 - 2 * x200 + x268)
            + x266 * x6
            + x267,
            x265,
            x0 * (-2 * x204 + 2 * x206 + 2 * x207 - x271)
            - x0 * (-2 * x210 + 2 * x211 - 2 * x212 + x271)
            + x269 * x6
            + x270,
            x0 * (-2 * x216 + 2 * x218 + 2 * x219 - x274)
            - x0 * (-2 * x222 + 2 * x223 - 2 * x224 + x274)
            + x272 * x6
            + x273,
            x265,
            x0 * x241 - x0 * x243 + x275 * x6 + x276,
            x0 * x250 - x0 * x252 + x277 * x6 + x278,
            x0 * x259 - x0 * x261 + x279 * x6 + x280,
            x0 * x123 - x0 * x130 + x281 * x6 + x282,
            x0 * x163 - x0 * x173 + x283 * x6 + x284,
            x0 * x54 - x0 * x62 + x285 * x6 + x286,
            x296,
            x303,
            x312,
            x133,
            x176,
            x65,
            x313,
            x266 * x7 + x267,
            x313,
            x269 * x7 + x270,
            x272 * x7 + x273,
            x313,
            x275 * x7 + x276,
            x277 * x7 + x278,
            x279 * x7 + x280,
            x281 * x7 + x282,
            x283 * x7 + x284,
            x285 * x7 + x286,
            x293 * x7 + x295,
            x300 * x7 + x302,
            x309 * x7 + x311,
            x124 * x7 + x132,
            x164 * x7 + x175,
            x55 * x7 + x64,
            x0
            * (
                -x19 * x322
                + 2 * x229
                - x314
                - x321
                + x38
                + x67 * (x159 - x287 + x289)
                + x99 * (x0 * (x13 * x319 - x23 * x316) + x119 + x228 - x318)
            )
            - x0
            * (
                -x19 * (x171 + x290 - x294)
                - 2 * x238
                + x314
                + x321
                + x322 * x67
                + x47
                - x99 * (-x0 * (x31 * x319 - x317 * x40) + x128 + x237 + x320)
            )
            + x296,
            x0
            * (
                x0 * (2 * x102 - 2 * x71)
                - x18 * x298
                + 2 * x244
                + x297 * x6
                - x323
                - x324
                + x80
            )
            - x0
            * (
                -x0 * (2 * x113 - 2 * x83)
                - x18 * x301
                - 2 * x247
                + x298 * x6
                + x323
                + x324
                + x87
            )
            + x303,
            x0 * (x0 * x227 - x18 * x307 + 2 * x253 + x305 * x6 - x325 - x326 + x38)
            - x0 * (-x0 * x236 - x18 * x310 - 2 * x256 + x307 * x6 + x325 + x326 + x47)
            + x312,
            x134,
            x177,
            x66,
        ]
    )
    return S


def coulomb_33(a, A, b, B, C):
    """Cartesian (ff) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = -x2 * (a * A[1] + b * B[1])
    x5 = -x2 * (a * A[2] + b * B[2])
    x6 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x4 + A[1]) ** 2 + abs(x5 + A[2]) ** 2)
    x7 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x4 + B[1]) ** 2 + abs(x5 + B[2]) ** 2)
    x8 = x3 + C[0]
    x9 = x4 + C[1]
    x10 = x5 + C[2]
    x11 = x10 ** 2 + x8 ** 2 + x9 ** 2
    x12 = x1 * x11
    x13 = boys(0, x12)
    x14 = numpy.pi * x2 * numpy.exp(-a * b * x11 * x2)
    x15 = 2 * x7
    x16 = x14 * x15
    x17 = boys(1, x12)
    x18 = numpy.sqrt(abs(x10) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x19 = 2 * x18
    x20 = x14 * x19
    x21 = -x17 * x20
    x22 = x13 * x16 + x21
    x23 = x14 * x17
    x24 = boys(2, x12)
    x25 = -x20 * x24
    x26 = x15 * x23 + x25
    x27 = x18 * x26
    x28 = -x27
    x29 = x22 * x7 + x28
    x30 = boys(3, x12)
    x31 = -x20 * x30
    x32 = x16 * x24 + x31
    x33 = x18 * x32
    x34 = -x33
    x35 = x26 * x7 + x34
    x36 = x18 * x35
    x37 = -x36
    x38 = x29 * x7 + x37
    x39 = x38 * x6
    x40 = boys(4, x12)
    x41 = -x20 * x40
    x42 = x16 * x30 + x41
    x43 = x18 * x42
    x44 = -x43
    x45 = x32 * x7 + x44
    x46 = x18 * x45
    x47 = -x46
    x48 = x35 * x7 + x47
    x49 = x18 * x48
    x50 = 2 * x39 - 2 * x49
    x51 = x48 * x6
    x52 = boys(5, x12)
    x53 = -x20 * x52
    x54 = x16 * x40 + x53
    x55 = x18 * x54
    x56 = -x55
    x57 = x42 * x7 + x56
    x58 = x18 * x57
    x59 = -x58
    x60 = x45 * x7 + x59
    x61 = x18 * x60
    x62 = 2 * x51 - 2 * x61
    x63 = x0 * x48
    x64 = -x49
    x65 = x39 + x64
    x66 = -x61
    x67 = x51 + x66
    x68 = -x18 * x67
    x69 = x6 * x65 + x68
    x70 = x0 * x38 - x63 + x69
    x71 = -x20 * boys(6, x12)
    x72 = x18 * (x16 * x52 + x71)
    x73 = -x72
    x74 = x18 * (x54 * x7 + x73)
    x75 = -x74
    x76 = -x18 * (x57 * x7 + x75)
    x77 = x6 * x60 + x76
    x78 = -x18 * x77
    x79 = x6 * x67 + x78
    x80 = -x0 * x60 + x63 + x79
    x81 = -x18 * x80
    x82 = x6 * x70 + x81
    x83 = x0 * x50 - x0 * x62 + x82
    x84 = 2 * x6
    x85 = x14 * x84
    x86 = x13 * x85 + x21
    x87 = x23 * x84 + x25
    x88 = x18 * x87
    x89 = -x88
    x90 = x7 * x86 + x89
    x91 = x24 * x85 + x31
    x92 = x18 * x91
    x93 = -x92
    x94 = x7 * x87 + x93
    x95 = x18 * x94
    x96 = -x95
    x97 = x7 * x90 + x96
    x98 = x6 * x97
    x99 = x30 * x85 + x41
    x100 = x18 * x99
    x101 = -x100
    x102 = x101 + x7 * x91
    x103 = x102 * x18
    x104 = -x103
    x105 = x104 + x7 * x94
    x106 = x105 * x18
    x107 = -2 * x106 + 2 * x98
    x108 = x105 * x6
    x109 = x40 * x85 + x53
    x110 = x109 * x18
    x111 = -x110
    x112 = x111 + x7 * x99
    x113 = x112 * x18
    x114 = -x113
    x115 = x102 * x7 + x114
    x116 = x115 * x18
    x117 = 2 * x108 - 2 * x116
    x118 = x0 * x105
    x119 = -x106
    x120 = x119 + x98
    x121 = -x116
    x122 = x108 + x121
    x123 = x0 * x97 - x118 + x120 * x6 - x122 * x18
    x124 = x52 * x85 + x71
    x125 = -x124 * x18
    x126 = x18 * (x109 * x7 + x125)
    x127 = -x126
    x128 = -x18 * (x112 * x7 + x127)
    x129 = x115 * x6 + x128
    x130 = -x0 * x115 + x118 + x122 * x6 - x129 * x18
    x131 = 2 * x0
    x132 = x131 * x23
    x133 = x131 * x14
    x134 = x6 * x86
    x135 = x134 + x89
    x136 = x13 * x133 - x132 + x135
    x137 = x133 * x24
    x138 = x6 * x87
    x139 = x138 + x93
    x140 = x132 - x137 + x139
    x141 = x140 * x18
    x142 = -x141
    x143 = x136 * x7 + x142
    x144 = x143 * x6
    x145 = x133 * x30
    x146 = x6 * x91
    x147 = x101 + x146
    x148 = x137 - x145 + x147
    x149 = x148 * x18
    x150 = -x149
    x151 = x140 * x7 + x150
    x152 = x151 * x18
    x153 = 2 * x144 - 2 * x152
    x154 = x151 * x6
    x155 = x133 * x40
    x156 = x6 * x99
    x157 = x111 + x156
    x158 = x145 - x155 + x157
    x159 = x158 * x18
    x160 = -x159
    x161 = x148 * x7 + x160
    x162 = x161 * x18
    x163 = 2 * x154 - 2 * x162
    x164 = x0 * x151
    x165 = -x152
    x166 = x144 + x165
    x167 = -x162
    x168 = x154 + x167
    x169 = x0 * x143 - x164 + x166 * x6 - x168 * x18
    x170 = x109 * x6 + x125
    x171 = -x133 * x52 + x155 + x170
    x172 = -x171 * x18
    x173 = -x18 * (x158 * x7 + x172)
    x174 = x161 * x6 + x173
    x175 = -x0 * x161 + x164 + x168 * x6 - x174 * x18
    x176 = x22 * x6
    x177 = x176 + x28
    x178 = x0 * x86
    x179 = x0 * x87
    x180 = x178 - x179
    x181 = x177 + x180
    x182 = x26 * x6
    x183 = x182 + x34
    x184 = x0 * x91
    x185 = x179 - x184
    x186 = x183 + x185
    x187 = x18 * x186
    x188 = -x187
    x189 = x181 * x7 + x188
    x190 = x189 * x6
    x191 = x32 * x6
    x192 = x191 + x44
    x193 = x0 * x99
    x194 = x184 - x193
    x195 = x192 + x194
    x196 = x18 * x195
    x197 = -x196
    x198 = x186 * x7 + x197
    x199 = x18 * x198
    x200 = 2 * x190 - 2 * x199
    x201 = x198 * x6
    x202 = x42 * x6
    x203 = x202 + x56
    x204 = x0 * x109
    x205 = x193 - x204
    x206 = x203 + x205
    x207 = x18 * x206
    x208 = -x207
    x209 = x195 * x7 + x208
    x210 = x18 * x209
    x211 = 2 * x201 - 2 * x210
    x212 = x0 * x198
    x213 = -x199
    x214 = x190 + x213
    x215 = -x210
    x216 = x201 + x215
    x217 = x0 * x189 - x18 * x216 - x212 + x214 * x6
    x218 = x0 * x124
    x219 = x54 * x6
    x220 = x219 + x73
    x221 = x18 * (x204 - x218 + x220)
    x222 = -x221
    x223 = -x18 * (x206 * x7 + x222)
    x224 = x209 * x6 + x223
    x225 = -x0 * x209 - x18 * x224 + x212 + x216 * x6
    x226 = 4 * x6
    x227 = x14 * x226
    x228 = 4 * x18
    x229 = x14 * x228
    x230 = x0 * (x226 * x23 - x229 * x24)
    x231 = x136 * x6
    x232 = x142 + x231
    x233 = x0 * (x13 * x227 - x228 * x23) - x230 + x232
    x234 = x233 * x6
    x235 = x0 * (x227 * x24 - x229 * x30)
    x236 = x140 * x6
    x237 = x150 + x236
    x238 = x230 - x235 + x237
    x239 = x18 * x238
    x240 = x238 * x6
    x241 = x0 * (x227 * x30 - x229 * x40)
    x242 = x148 * x6
    x243 = x160 + x242
    x244 = x235 - x241 + x243
    x245 = x18 * x244
    x246 = x0 * x238
    x247 = -x239
    x248 = x234 + x247
    x249 = -x245
    x250 = x240 + x249
    x251 = -x18 * x250
    x252 = x248 * x6 + x251
    x253 = x0 * x233 - x246 + x252
    x254 = x158 * x6 + x172
    x255 = -x0 * (x227 * x40 - x229 * x52) + x241 + x254
    x256 = -x18 * x255
    x257 = x244 * x6 + x256
    x258 = -x18 * x257
    x259 = x250 * x6 + x258
    x260 = -x0 * x244 + x246 + x259
    x261 = -x18 * x260
    x262 = x253 * x6 + x261
    x263 = x0 * (2 * x234 - 2 * x239) - x0 * (2 * x240 - 2 * x245) + x262
    x264 = x0 * (x136 + x22)
    x265 = x0 * (x140 + x26)
    x266 = x181 * x6
    x267 = x188 + x266
    x268 = x264 - x265 + x267
    x269 = x268 * x6
    x270 = x0 * (x148 + x32)
    x271 = x186 * x6
    x272 = x197 + x271
    x273 = x265 - x270 + x272
    x274 = x18 * x273
    x275 = 2 * x269 - 2 * x274
    x276 = x273 * x6
    x277 = x0 * (x158 + x42)
    x278 = x195 * x6
    x279 = x208 + x278
    x280 = x270 - x277 + x279
    x281 = x18 * x280
    x282 = 2 * x276 - 2 * x281
    x283 = x0 * x273
    x284 = -x274
    x285 = x269 + x284
    x286 = -x281
    x287 = x276 + x286
    x288 = -x18 * x287
    x289 = x285 * x6 + x288
    x290 = x0 * x268 - x283 + x289
    x291 = x280 * x6
    x292 = x0 * (x171 + x54)
    x293 = x206 * x6
    x294 = x222 + x293
    x295 = x18 * (x277 - x292 + x294)
    x296 = -x295
    x297 = x291 + x296
    x298 = -x18 * x297
    x299 = x287 * x6 + x298
    x300 = -x0 * x280 + x283 + x299
    x301 = -x18 * x300
    x302 = x290 * x6 + x301
    x303 = x0 * x275 - x0 * x282 + x302
    x304 = 2 * x179
    x305 = 2 * x176 - 2 * x27
    x306 = 2 * x178 - x304 + x305
    x307 = x0 * x306
    x308 = 2 * x184
    x309 = 2 * x182 - 2 * x33
    x310 = x304 - x308 + x309
    x311 = x0 * x310
    x312 = x29 * x6
    x313 = x312 + x37
    x314 = x307 - x311 + x313
    x315 = x314 * x6
    x316 = 2 * x193
    x317 = 2 * x191 - 2 * x43
    x318 = x308 - x316 + x317
    x319 = x0 * x318
    x320 = x35 * x6
    x321 = x320 + x47
    x322 = x311 - x319 + x321
    x323 = x18 * x322
    x324 = x322 * x6
    x325 = 2 * x204
    x326 = 2 * x202 - 2 * x55
    x327 = x316 - x325 + x326
    x328 = x0 * x327
    x329 = x45 * x6
    x330 = x329 + x59
    x331 = x319 - x328 + x330
    x332 = x18 * x331
    x333 = x0 * x322
    x334 = -x323
    x335 = x315 + x334
    x336 = -x332
    x337 = x324 + x336
    x338 = -x18 * x337
    x339 = x335 * x6 + x338
    x340 = x0 * x314 - x333 + x339
    x341 = x331 * x6
    x342 = x0 * (-2 * x218 + 2 * x219 + x325 - 2 * x72)
    x343 = x57 * x6
    x344 = x343 + x75
    x345 = x18 * (x328 - x342 + x344)
    x346 = -x345
    x347 = x341 + x346
    x348 = -x18 * x347
    x349 = x337 * x6 + x348
    x350 = -x0 * x331 + x333 + x349
    x351 = -x18 * x350
    x352 = x340 * x6 + x351
    x353 = x0 * (2 * x315 - 2 * x323) - x0 * (2 * x324 - 2 * x332) + x352
    x354 = x38 * x7 + x64
    x355 = x48 * x7 + x66
    x356 = x18 * x355
    x357 = -x356
    x358 = x354 * x6
    x359 = x0 * x65
    x360 = x0 * x67
    x361 = x359 - x360
    x362 = x0 * x77
    x363 = x355 * x6
    x364 = x18 * (x60 * x7 + x76)
    x365 = -x364
    x366 = (
        x0 * (x354 + x70)
        - x0 * (x355 + x80)
        - x18 * (x360 - x362 + x363 + x365)
        + x6 * (x357 + x358 + x361)
    )
    x367 = x119 + x7 * x97
    x368 = x105 * x7 + x121
    x369 = x0 * x120
    x370 = x0 * x122
    x371 = x367 * x6
    x372 = x18 * x368
    x373 = -x372
    x374 = x0 * x129
    x375 = x368 * x6
    x376 = x18 * (x115 * x7 + x128)
    x377 = -x376
    x378 = x143 * x7 + x165
    x379 = x151 * x7 + x167
    x380 = x0 * x166
    x381 = x0 * x168
    x382 = x378 * x6
    x383 = x18 * x379
    x384 = -x383
    x385 = x0 * x174
    x386 = x379 * x6
    x387 = x18 * (x161 * x7 + x173)
    x388 = -x387
    x389 = x189 * x7 + x213
    x390 = x198 * x7 + x215
    x391 = x0 * x214
    x392 = x0 * x216
    x393 = x389 * x6
    x394 = x18 * x390
    x395 = -x394
    x396 = x0 * x224
    x397 = x390 * x6
    x398 = x18 * (x209 * x7 + x223)
    x399 = -x398
    x400 = x233 * x7 + x247
    x401 = x238 * x7 + x249
    x402 = x18 * x401
    x403 = -x402
    x404 = x400 * x6
    x405 = x0 * x248
    x406 = x0 * x250
    x407 = x405 - x406
    x408 = x0 * x257
    x409 = x401 * x6
    x410 = x18 * (x244 * x7 + x256)
    x411 = -x410
    x412 = x268 * x7 + x284
    x413 = x273 * x7 + x286
    x414 = x18 * x413
    x415 = -x414
    x416 = x412 * x6
    x417 = x0 * x285
    x418 = x0 * x287
    x419 = x417 - x418
    x420 = x0 * x297
    x421 = x413 * x6
    x422 = x18 * (x280 * x7 + x296)
    x423 = -x422
    x424 = x314 * x7 + x334
    x425 = x322 * x7 + x336
    x426 = x18 * x425
    x427 = -x426
    x428 = x424 * x6
    x429 = x0 * x335
    x430 = x0 * x337
    x431 = x429 - x430
    x432 = x0 * x347
    x433 = x425 * x6
    x434 = x18 * (x331 * x7 + x346)
    x435 = -x434
    x436 = 3 * x311
    x437 = x0 * (3 * x307 + 3 * x312 - 3 * x36 - x436)
    x438 = 3 * x319
    x439 = x0 * (3 * x320 + x436 - x438 - 3 * x46)
    x440 = x437 - x439 + x65
    x441 = 3 * x328
    x442 = x0 * (3 * x329 + x438 - x441 - 3 * x58)
    x443 = x439 - x442 + x67
    x444 = -x18 * x443
    x445 = x440 * x6 + x444
    x446 = -x18 * (-x0 * (-3 * x342 + 3 * x343 + x441 - 3 * x74) + x442 + x77)
    x447 = x443 * x6 + x446
    x448 = x0 * x135
    x449 = x0 * x139
    x450 = 2 * x449
    x451 = x6 * x90
    x452 = 2 * x448 - x450 + 2 * x451 - 2 * x95
    x453 = x0 * x452
    x454 = x0 * x147
    x455 = 2 * x454
    x456 = x6 * x94
    x457 = -2 * x103 + x450 - x455 + 2 * x456
    x458 = x0 * x457
    x459 = x120 + x453 - x458
    x460 = x0 * x157
    x461 = 2 * x460
    x462 = x102 * x6
    x463 = -2 * x113 + x455 - x461 + 2 * x462
    x464 = x0 * x463
    x465 = x122 + x458 - x464
    x466 = -x18 * x465
    x467 = x459 * x6 + x466
    x468 = x0 * x170
    x469 = x112 * x6
    x470 = -x18 * (-x0 * (-2 * x126 + x461 - 2 * x468 + 2 * x469) + x129 + x464)
    x471 = x465 * x6 + x470
    x472 = x0 * x183
    x473 = 2 * x472
    x474 = x0 * x177
    x475 = 2 * x312 - 2 * x36
    x476 = -x473 + 2 * x474 + x475
    x477 = x0 * x476
    x478 = x0 * x192
    x479 = 2 * x478
    x480 = 2 * x320 - 2 * x46
    x481 = x473 - x479 + x480
    x482 = x0 * x481
    x483 = x477 - x482 + x65
    x484 = x0 * x203
    x485 = 2 * x484
    x486 = 2 * x329 - 2 * x58
    x487 = x479 - x485 + x486
    x488 = x0 * x487
    x489 = x482 - x488 + x67
    x490 = -x18 * x489
    x491 = x483 * x6 + x490
    x492 = x0 * x220
    x493 = -x18 * (-x0 * (2 * x343 + x485 - 2 * x492 - 2 * x74) + x488 + x77)
    x494 = x489 * x6 + x493
    x495 = x0 * x232
    x496 = x0 * x237
    x497 = x166 + x495 - x496
    x498 = x0 * x243
    x499 = x168 + x496 - x498
    x500 = -x18 * x499
    x501 = x497 * x6 + x500
    x502 = -x18 * (-x0 * x254 + x174 + x498)
    x503 = x499 * x6 + x502
    x504 = x0 * x267
    x505 = x0 * x272
    x506 = x214 + x504 - x505
    x507 = x0 * x279
    x508 = x216 + x505 - x507
    x509 = -x18 * x508
    x510 = x506 * x6 + x509
    x511 = -x18 * (-x0 * x294 + x224 + x507)
    x512 = x508 * x6 + x511
    x513 = x0 * x313
    x514 = x0 * x321
    x515 = x513 - x514 + x65
    x516 = x0 * x330
    x517 = x514 - x516 + x67
    x518 = -x18 * x517
    x519 = x515 * x6 + x518
    x520 = -x18 * (-x0 * x344 + x516 + x77)
    x521 = x517 * x6 + x520
    x522 = x354 * x7 + x357
    x523 = -x18 * (x355 * x7 + x365)
    x524 = 2 * x360
    x525 = (
        x0 * (-2 * x356 + 2 * x358 + 2 * x359 - x524)
        - x0 * (-2 * x362 + 2 * x363 - 2 * x364 + x524)
        + x522 * x6
        + x523
    )
    x526 = x367 * x7 + x373
    x527 = -x18 * (x368 * x7 + x377)
    x528 = 2 * x370
    x529 = x378 * x7 + x384
    x530 = -x18 * (x379 * x7 + x388)
    x531 = 2 * x381
    x532 = x389 * x7 + x395
    x533 = -x18 * (x390 * x7 + x399)
    x534 = 2 * x392
    x535 = x400 * x7 + x403
    x536 = -x18 * (x401 * x7 + x411)
    x537 = 2 * x406
    x538 = x412 * x7 + x415
    x539 = -x18 * (x413 * x7 + x423)
    x540 = 2 * x418
    x541 = x424 * x7 + x427
    x542 = -x18 * (x425 * x7 + x435)
    x543 = 2 * x430
    x544 = x440 * x7 + x444
    x545 = -x18 * (x443 * x7 + x446)
    x546 = x459 * x7 + x466
    x547 = -x18 * (x465 * x7 + x470)
    x548 = x483 * x7 + x490
    x549 = -x18 * (x489 * x7 + x493)
    x550 = x497 * x7 + x500
    x551 = -x18 * (x499 * x7 + x502)
    x552 = x506 * x7 + x509
    x553 = -x18 * (x508 * x7 + x511)
    x554 = x515 * x7 + x518
    x555 = -x18 * (x517 * x7 + x520)
    x556 = x248 * x7 + x251
    x557 = -x18 * (x250 * x7 + x258)
    x558 = x285 * x7 + x288
    x559 = -x18 * (x287 * x7 + x298)
    x560 = x335 * x7 + x338
    x561 = -x18 * (x337 * x7 + x348)
    x562 = x65 * x7 + x68
    x563 = -x18 * (x67 * x7 + x78)
    x564 = 2 * x270
    x565 = 2 * x265
    x566 = -2 * x196 + 2 * x271
    x567 = x0 * (x35 - x564 + x565 + x566)
    x568 = 3 * x567
    x569 = -2 * x187 + 2 * x266
    x570 = x0 * (2 * x264 + x29 - x565 + x569)
    x571 = 2 * x277
    x572 = -2 * x207 + 2 * x278
    x573 = x0 * (x45 + x564 - x571 + x572)
    x574 = 3 * x573
    x575 = x0 * (3 * x324 - 3 * x332 + x48 + x568 - x574)
    x576 = x0 * (3 * x315 - 3 * x323 + x38 - x568 + 3 * x570) + x445 - x575
    x577 = x0 * (-2 * x221 - 2 * x292 + 2 * x293 + x57 + x571)
    x578 = -x18 * (-x0 * (3 * x341 - 3 * x345 + x574 - 3 * x577 + x60) + x447 + x575)
    x579 = x576 * x6 + x578
    x580 = x139 * x6 - x147 * x18 + x185
    x581 = x0 * (x580 + x94)
    x582 = 2 * x581
    x583 = x104 + x449 - x454 + x456
    x584 = x18 * x583
    x585 = x135 * x6 - x139 * x18 + x180
    x586 = x0 * (x585 + x90)
    x587 = x6 * (x448 - x449 + x451 + x96)
    x588 = x147 * x6 - x157 * x18 + x194
    x589 = x0 * (x102 + x588)
    x590 = 2 * x589
    x591 = x114 + x454 - x460 + x462
    x592 = x18 * x591
    x593 = x583 * x6
    x594 = x0 * (x105 + x582 - x590 - 2 * x592 + 2 * x593)
    x595 = x0 * (-x582 - 2 * x584 + 2 * x586 + 2 * x587 + x97) + x467 - x594
    x596 = x157 * x6 - x170 * x18 + x205
    x597 = x0 * (x112 + x596)
    x598 = x18 * (x127 + x460 - x468 + x469)
    x599 = x591 * x6
    x600 = -x18 * (-x0 * (x115 + x590 - 2 * x597 - 2 * x598 + 2 * x599) + x471 + x594)
    x601 = x595 * x6 + x600
    x602 = x0 * x26
    x603 = x0 * x32
    x604 = -x18 * x192 + x183 * x6 + x602 - x603
    x605 = x0 * (x35 + x604)
    x606 = 2 * x605
    x607 = x321 + x472 - x478
    x608 = x18 * x607
    x609 = x0 * x22 + x177 * x6 - x18 * x183 - x602
    x610 = x0 * (x29 + x609)
    x611 = x6 * (x313 - x472 + x474)
    x612 = x0 * x42
    x613 = -x18 * x203 + x192 * x6 + x603 - x612
    x614 = x0 * (x45 + x613)
    x615 = 2 * x614
    x616 = x330 + x478 - x484
    x617 = x18 * x616
    x618 = x6 * x607
    x619 = x0 * (x48 + x606 - x615 - 2 * x617 + 2 * x618)
    x620 = x0 * (x38 - x606 - 2 * x608 + 2 * x610 + 2 * x611) + x491 - x619
    x621 = -x0 * x54 - x18 * x220 + x203 * x6 + x612
    x622 = x0 * (x57 + x621)
    x623 = x18 * (x344 + x484 - x492)
    x624 = x6 * x616
    x625 = -x18 * (-x0 * (x60 + x615 - 2 * x622 - 2 * x623 + 2 * x624) + x494 + x619)
    x626 = x6 * x620 + x625
    x627 = x0 * x140
    x628 = x0 * x136 - x18 * x237 + x232 * x6 - x627
    x629 = x0 * x148
    x630 = -x18 * x243 + x237 * x6 + x627 - x629
    x631 = x0 * (x151 + x630)
    x632 = x0 * (x143 + x628) + x501 - x631
    x633 = -x0 * x158 - x18 * x254 + x243 * x6 + x629
    x634 = -x18 * (-x0 * (x161 + x633) + x503 + x631)
    x635 = x6 * x632 + x634
    x636 = x0 * x186
    x637 = x0 * x181 - x18 * x272 + x267 * x6 - x636
    x638 = x0 * x195
    x639 = -x18 * x279 + x272 * x6 + x636 - x638
    x640 = x0 * (x198 + x639)
    x641 = x0 * (x189 + x637) + x510 - x640
    x642 = -x0 * x206 - x18 * x294 + x279 * x6 + x638
    x643 = -x18 * (-x0 * (x209 + x642) + x512 + x640)
    x644 = x6 * x641 + x643
    x645 = x0 * x35
    x646 = x0 * x29 - x18 * x321 + x313 * x6 - x645
    x647 = x0 * x45
    x648 = -x18 * x330 + x321 * x6 + x645 - x647
    x649 = x0 * (x48 + x648)
    x650 = x0 * (x38 + x646) + x519 - x649
    x651 = -x0 * x57 - x18 * x344 + x330 * x6 + x647
    x652 = -x18 * (-x0 * (x60 + x651) + x521 + x649)
    x653 = x6 * x650 + x652
    x654 = x522 * x7 + x523
    x655 = 2 * x319
    x656 = x131 * (x244 + x318)
    x657 = 2 * x311
    x658 = x131 * (x238 + x310)
    x659 = 3 * x0
    x660 = x659 * (x282 + x480 - x655 - x656 + x657 + x658)
    x661 = x337 + x567 - x573
    x662 = 3 * x18
    x663 = 2 * x439
    x664 = 3 * x6
    x665 = 2 * x458
    x666 = x0 * (2 * x138 - 2 * x92)
    x667 = x0 * (-2 * x100 + 2 * x146)
    x668 = x131 * (-x18 * x588 + x457 + x580 * x6 + x666 - x667)
    x669 = x581 - x589 - x592 + x593
    x670 = 2 * x482
    x671 = x0 * x309
    x672 = x0 * x317
    x673 = x131 * (-x18 * x613 + x481 + x6 * x604 + x671 - x672)
    x674 = x605 - x614 - x617 + x618
    x675 = x0 * (-2 * x149 + 2 * x236)
    x676 = 2 * x496
    x677 = x0 * x566
    x678 = 2 * x505
    x679 = x0 * x480
    x680 = 2 * x514

    # 100 item(s)
    S = numpy.array(
        [
            x83,
            x0 * x107 - x0 * x117 + x123 * x6 - x130 * x18,
            x83,
            x0 * x153 - x0 * x163 + x169 * x6 - x175 * x18,
            x0 * x200 - x0 * x211 - x18 * x225 + x217 * x6,
            x83,
            x263,
            x303,
            x353,
            x83,
            x366,
            x0 * (x123 + x367)
            - x0 * (x130 + x368)
            - x18 * (x370 - x374 + x375 + x377)
            + x6 * (x369 - x370 + x371 + x373),
            x366,
            x0 * (x169 + x378)
            - x0 * (x175 + x379)
            - x18 * (x381 - x385 + x386 + x388)
            + x6 * (x380 - x381 + x382 + x384),
            x0 * (x217 + x389)
            - x0 * (x225 + x390)
            - x18 * (x392 - x396 + x397 + x399)
            + x6 * (x391 - x392 + x393 + x395),
            x366,
            x0 * (x253 + x400)
            - x0 * (x260 + x401)
            - x18 * (x406 - x408 + x409 + x411)
            + x6 * (x403 + x404 + x407),
            x0 * (x290 + x412)
            - x0 * (x300 + x413)
            - x18 * (x418 - x420 + x421 + x423)
            + x6 * (x415 + x416 + x419),
            x0 * (x340 + x424)
            - x0 * (x350 + x425)
            - x18 * (x430 - x432 + x433 + x435)
            + x6 * (x427 + x428 + x431),
            x366,
            x0 * x440 - x0 * x443 - x18 * x447 + x445 * x6,
            x0 * x459 - x0 * x465 - x18 * x471 + x467 * x6,
            x0 * x483 - x0 * x489 - x18 * x494 + x491 * x6,
            x0 * x497 - x0 * x499 - x18 * x503 + x501 * x6,
            x0 * x506 - x0 * x508 - x18 * x512 + x510 * x6,
            x0 * x515 - x0 * x517 - x18 * x521 + x519 * x6,
            -x18 * x259 + x252 * x6 + x407,
            -x18 * x299 + x289 * x6 + x419,
            -x18 * x349 + x339 * x6 + x431,
            -x18 * x79 + x361 + x6 * x69,
            x525,
            x0 * (2 * x369 + 2 * x371 - 2 * x372 - x528)
            - x0 * (-2 * x374 + 2 * x375 - 2 * x376 + x528)
            + x526 * x6
            + x527,
            x525,
            x0 * (2 * x380 + 2 * x382 - 2 * x383 - x531)
            - x0 * (-2 * x385 + 2 * x386 - 2 * x387 + x531)
            + x529 * x6
            + x530,
            x0 * (2 * x391 + 2 * x393 - 2 * x394 - x534)
            - x0 * (-2 * x396 + 2 * x397 - 2 * x398 + x534)
            + x532 * x6
            + x533,
            x525,
            x0 * (-2 * x402 + 2 * x404 + 2 * x405 - x537)
            - x0 * (-2 * x408 + 2 * x409 - 2 * x410 + x537)
            + x535 * x6
            + x536,
            x0 * (-2 * x414 + 2 * x416 + 2 * x417 - x540)
            - x0 * (-2 * x420 + 2 * x421 - 2 * x422 + x540)
            + x538 * x6
            + x539,
            x0 * (-2 * x426 + 2 * x428 + 2 * x429 - x543)
            - x0 * (-2 * x432 + 2 * x433 - 2 * x434 + x543)
            + x541 * x6
            + x542,
            x525,
            x0 * x445 - x0 * x447 + x544 * x6 + x545,
            x0 * x467 - x0 * x471 + x546 * x6 + x547,
            x0 * x491 - x0 * x494 + x548 * x6 + x549,
            x0 * x501 - x0 * x503 + x550 * x6 + x551,
            x0 * x510 - x0 * x512 + x552 * x6 + x553,
            x0 * x519 - x0 * x521 + x554 * x6 + x555,
            x0 * x252 - x0 * x259 + x556 * x6 + x557,
            x0 * x289 - x0 * x299 + x558 * x6 + x559,
            x0 * x339 - x0 * x349 + x560 * x6 + x561,
            x0 * x69 - x0 * x79 + x562 * x6 + x563,
            x579,
            x601,
            x626,
            x635,
            x644,
            x653,
            x262,
            x302,
            x352,
            x82,
            x654,
            x526 * x7 + x527,
            x654,
            x529 * x7 + x530,
            x532 * x7 + x533,
            x654,
            x535 * x7 + x536,
            x538 * x7 + x539,
            x541 * x7 + x542,
            x654,
            x544 * x7 + x545,
            x546 * x7 + x547,
            x548 * x7 + x549,
            x550 * x7 + x551,
            x552 * x7 + x553,
            x554 * x7 + x555,
            x556 * x7 + x557,
            x558 * x7 + x559,
            x560 * x7 + x561,
            x562 * x7 + x563,
            x576 * x7 + x578,
            x595 * x7 + x600,
            x620 * x7 + x625,
            x632 * x7 + x634,
            x641 * x7 + x643,
            x650 * x7 + x652,
            x253 * x7 + x261,
            x290 * x7 + x301,
            x340 * x7 + x351,
            x7 * x70 + x81,
            x0
            * (
                2 * x437
                + x50
                + x659 * (x131 * (x233 + x306) + x275 + 2 * x307 + x475 - x657 - x658)
                - x660
                - x661 * x662
                - x663
                + x664 * (x335 - x567 + x570)
            )
            - x0
            * (
                -2 * x442
                + x62
                - x659
                * (
                    -x131 * (x255 + x327)
                    + 2 * x291
                    - 2 * x295
                    - 2 * x328
                    + x486
                    + x655
                    + x656
                )
                + x660
                + x661 * x664
                - x662 * (x347 + x573 - x577)
                + x663
            )
            + x579,
            x0
            * (
                x107
                + x131
                * (x0 * (2 * x134 - 2 * x88) - x18 * x580 + x452 + x585 * x6 - x666)
                - x19 * x669
                + 2 * x453
                - x665
                - x668
                + x84 * (-x581 - x584 + x586 + x587)
            )
            - x0
            * (
                x117
                - x131
                * (-x0 * (-2 * x110 + 2 * x156) - x18 * x596 + x463 + x588 * x6 + x667)
                - x19 * (x589 - x597 - x598 + x599)
                - 2 * x464
                + x665
                + x668
                + x669 * x84
            )
            + x601,
            x0
            * (
                x131 * (x0 * x305 - x18 * x604 + x476 + x6 * x609 - x671)
                - x19 * x674
                + 2 * x477
                + x50
                - x670
                - x673
                + x84 * (-x605 - x608 + x610 + x611)
            )
            - x0
            * (
                -x131 * (-x0 * x326 - x18 * x621 + x487 + x6 * x613 + x672)
                - x19 * (x614 - x622 - x623 + x624)
                - 2 * x488
                + x62
                + x670
                + x673
                + x674 * x84
            )
            + x626,
            x0
            * (
                x0 * (-2 * x141 + 2 * x231)
                + x153
                - x18 * x630
                + 2 * x495
                + x6 * x628
                - x675
                - x676
            )
            - x0
            * (
                -x0 * (-2 * x159 + 2 * x242)
                + x163
                - x18 * x633
                - 2 * x498
                + x6 * x630
                + x675
                + x676
            )
            + x635,
            x0 * (x0 * x569 - x18 * x639 + x200 + 2 * x504 + x6 * x637 - x677 - x678)
            - x0 * (-x0 * x572 - x18 * x642 + x211 - 2 * x507 + x6 * x639 + x677 + x678)
            + x644,
            x0 * (x0 * x475 - x18 * x648 + x50 + 2 * x513 + x6 * x646 - x679 - x680)
            - x0 * (-x0 * x486 - x18 * x651 - 2 * x516 + x6 * x648 + x62 + x679 + x680)
            + x653,
            x263,
            x303,
            x353,
            x83,
        ]
    )
    return S


def coulomb_34(a, A, b, B, C):
    """Cartesian (fg) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = -x2 * (a * A[1] + b * B[1])
    x5 = -x2 * (a * A[2] + b * B[2])
    x6 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x4 + A[1]) ** 2 + abs(x5 + A[2]) ** 2)
    x7 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x4 + B[1]) ** 2 + abs(x5 + B[2]) ** 2)
    x8 = x3 + C[0]
    x9 = x4 + C[1]
    x10 = x5 + C[2]
    x11 = x10 ** 2 + x8 ** 2 + x9 ** 2
    x12 = x1 * x11
    x13 = boys(0, x12)
    x14 = numpy.pi * x2 * numpy.exp(-a * b * x11 * x2)
    x15 = 2 * x7
    x16 = x14 * x15
    x17 = boys(1, x12)
    x18 = numpy.sqrt(abs(x10) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x19 = 2 * x18
    x20 = x14 * x19
    x21 = -x17 * x20
    x22 = x13 * x16 + x21
    x23 = x14 * x17
    x24 = boys(2, x12)
    x25 = -x20 * x24
    x26 = x15 * x23 + x25
    x27 = x18 * x26
    x28 = -x27
    x29 = x22 * x7 + x28
    x30 = boys(3, x12)
    x31 = -x20 * x30
    x32 = x16 * x24 + x31
    x33 = x18 * x32
    x34 = -x33
    x35 = x26 * x7 + x34
    x36 = x18 * x35
    x37 = -x36
    x38 = x29 * x7 + x37
    x39 = boys(4, x12)
    x40 = -x20 * x39
    x41 = x16 * x30 + x40
    x42 = x18 * x41
    x43 = -x42
    x44 = x32 * x7 + x43
    x45 = x18 * x44
    x46 = -x45
    x47 = x35 * x7 + x46
    x48 = x18 * x47
    x49 = -x48
    x50 = x38 * x7 + x49
    x51 = x50 * x6
    x52 = boys(5, x12)
    x53 = -x20 * x52
    x54 = x16 * x39 + x53
    x55 = x18 * x54
    x56 = -x55
    x57 = x41 * x7 + x56
    x58 = x18 * x57
    x59 = -x58
    x60 = x44 * x7 + x59
    x61 = x18 * x60
    x62 = -x61
    x63 = x47 * x7 + x62
    x64 = x18 * x63
    x65 = 2 * x51 - 2 * x64
    x66 = x6 * x63
    x67 = boys(6, x12)
    x68 = -x20 * x67
    x69 = x16 * x52 + x68
    x70 = x18 * x69
    x71 = -x70
    x72 = x54 * x7 + x71
    x73 = x18 * x72
    x74 = -x73
    x75 = x57 * x7 + x74
    x76 = x18 * x75
    x77 = -x76
    x78 = x60 * x7 + x77
    x79 = x18 * x78
    x80 = 2 * x66 - 2 * x79
    x81 = x0 * x63
    x82 = -x64
    x83 = x51 + x82
    x84 = -x79
    x85 = x66 + x84
    x86 = -x18 * x85
    x87 = x6 * x83 + x86
    x88 = x0 * x50 - x81 + x87
    x89 = -x20 * boys(7, x12)
    x90 = x18 * (x16 * x67 + x89)
    x91 = -x90
    x92 = x18 * (x69 * x7 + x91)
    x93 = -x92
    x94 = x18 * (x7 * x72 + x93)
    x95 = -x94
    x96 = -x18 * (x7 * x75 + x95)
    x97 = x6 * x78 + x96
    x98 = -x18 * x97
    x99 = x6 * x85 + x98
    x100 = -x0 * x78 + x81 + x99
    x101 = -x100 * x18
    x102 = x101 + x6 * x88
    x103 = x0 * x65 - x0 * x80 + x102
    x104 = 2 * x6
    x105 = x104 * x14
    x106 = x105 * x13 + x21
    x107 = x104 * x23 + x25
    x108 = x107 * x18
    x109 = -x108
    x110 = x106 * x7 + x109
    x111 = x105 * x24 + x31
    x112 = x111 * x18
    x113 = -x112
    x114 = x107 * x7 + x113
    x115 = x114 * x18
    x116 = -x115
    x117 = x110 * x7 + x116
    x118 = x105 * x30 + x40
    x119 = x118 * x18
    x120 = -x119
    x121 = x111 * x7 + x120
    x122 = x121 * x18
    x123 = -x122
    x124 = x114 * x7 + x123
    x125 = x124 * x18
    x126 = -x125
    x127 = x117 * x7 + x126
    x128 = x127 * x6
    x129 = x105 * x39 + x53
    x130 = x129 * x18
    x131 = -x130
    x132 = x118 * x7 + x131
    x133 = x132 * x18
    x134 = -x133
    x135 = x121 * x7 + x134
    x136 = x135 * x18
    x137 = -x136
    x138 = x124 * x7 + x137
    x139 = x138 * x18
    x140 = 2 * x128 - 2 * x139
    x141 = x138 * x6
    x142 = x105 * x52 + x68
    x143 = x142 * x18
    x144 = -x143
    x145 = x129 * x7 + x144
    x146 = x145 * x18
    x147 = -x146
    x148 = x132 * x7 + x147
    x149 = x148 * x18
    x150 = -x149
    x151 = x135 * x7 + x150
    x152 = x151 * x18
    x153 = 2 * x141 - 2 * x152
    x154 = x0 * x138
    x155 = -x139
    x156 = x128 + x155
    x157 = -x152
    x158 = x141 + x157
    x159 = x0 * x127 - x154 + x156 * x6 - x158 * x18
    x160 = x105 * x67 + x89
    x161 = -x160 * x18
    x162 = x18 * (x142 * x7 + x161)
    x163 = -x162
    x164 = x18 * (x145 * x7 + x163)
    x165 = -x164
    x166 = -x18 * (x148 * x7 + x165)
    x167 = x151 * x6 + x166
    x168 = -x0 * x151 + x154 + x158 * x6 - x167 * x18
    x169 = 2 * x0
    x170 = x169 * x23
    x171 = x14 * x169
    x172 = x106 * x6
    x173 = x109 + x172
    x174 = x13 * x171 - x170 + x173
    x175 = x171 * x24
    x176 = x107 * x6
    x177 = x113 + x176
    x178 = x170 - x175 + x177
    x179 = x178 * x18
    x180 = -x179
    x181 = x174 * x7 + x180
    x182 = x171 * x30
    x183 = x111 * x6
    x184 = x120 + x183
    x185 = x175 - x182 + x184
    x186 = x18 * x185
    x187 = -x186
    x188 = x178 * x7 + x187
    x189 = x18 * x188
    x190 = -x189
    x191 = x181 * x7 + x190
    x192 = x191 * x6
    x193 = x171 * x39
    x194 = x118 * x6
    x195 = x131 + x194
    x196 = x182 - x193 + x195
    x197 = x18 * x196
    x198 = -x197
    x199 = x185 * x7 + x198
    x200 = x18 * x199
    x201 = -x200
    x202 = x188 * x7 + x201
    x203 = x18 * x202
    x204 = 2 * x192 - 2 * x203
    x205 = x202 * x6
    x206 = x171 * x52
    x207 = x129 * x6
    x208 = x144 + x207
    x209 = x193 - x206 + x208
    x210 = x18 * x209
    x211 = -x210
    x212 = x196 * x7 + x211
    x213 = x18 * x212
    x214 = -x213
    x215 = x199 * x7 + x214
    x216 = x18 * x215
    x217 = 2 * x205 - 2 * x216
    x218 = x0 * x202
    x219 = -x203
    x220 = x192 + x219
    x221 = -x216
    x222 = x205 + x221
    x223 = x0 * x191 - x18 * x222 - x218 + x220 * x6
    x224 = x142 * x6 + x161
    x225 = -x171 * x67 + x206 + x224
    x226 = -x18 * x225
    x227 = x18 * (x209 * x7 + x226)
    x228 = -x227
    x229 = -x18 * (x212 * x7 + x228)
    x230 = x215 * x6 + x229
    x231 = -x0 * x215 - x18 * x230 + x218 + x222 * x6
    x232 = x22 * x6
    x233 = x232 + x28
    x234 = x0 * x106
    x235 = x0 * x107
    x236 = x234 - x235
    x237 = x233 + x236
    x238 = x26 * x6
    x239 = x238 + x34
    x240 = x0 * x111
    x241 = x235 - x240
    x242 = x239 + x241
    x243 = x18 * x242
    x244 = -x243
    x245 = x237 * x7 + x244
    x246 = x32 * x6
    x247 = x246 + x43
    x248 = x0 * x118
    x249 = x240 - x248
    x250 = x247 + x249
    x251 = x18 * x250
    x252 = -x251
    x253 = x242 * x7 + x252
    x254 = x18 * x253
    x255 = -x254
    x256 = x245 * x7 + x255
    x257 = x256 * x6
    x258 = x41 * x6
    x259 = x258 + x56
    x260 = x0 * x129
    x261 = x248 - x260
    x262 = x259 + x261
    x263 = x18 * x262
    x264 = -x263
    x265 = x250 * x7 + x264
    x266 = x18 * x265
    x267 = -x266
    x268 = x253 * x7 + x267
    x269 = x18 * x268
    x270 = 2 * x257 - 2 * x269
    x271 = x268 * x6
    x272 = x54 * x6
    x273 = x272 + x71
    x274 = x0 * x142
    x275 = x260 - x274
    x276 = x273 + x275
    x277 = x18 * x276
    x278 = -x277
    x279 = x262 * x7 + x278
    x280 = x18 * x279
    x281 = -x280
    x282 = x265 * x7 + x281
    x283 = x18 * x282
    x284 = 2 * x271 - 2 * x283
    x285 = x0 * x268
    x286 = -x269
    x287 = x257 + x286
    x288 = -x283
    x289 = x271 + x288
    x290 = x0 * x256 - x18 * x289 - x285 + x287 * x6
    x291 = x0 * x160
    x292 = x6 * x69
    x293 = x292 + x91
    x294 = x18 * (x274 - x291 + x293)
    x295 = -x294
    x296 = x18 * (x276 * x7 + x295)
    x297 = -x296
    x298 = -x18 * (x279 * x7 + x297)
    x299 = x282 * x6 + x298
    x300 = -x0 * x282 - x18 * x299 + x285 + x289 * x6
    x301 = 4 * x6
    x302 = x14 * x301
    x303 = 4 * x18
    x304 = x14 * x303
    x305 = x0 * (x23 * x301 - x24 * x304)
    x306 = x174 * x6
    x307 = x180 + x306
    x308 = x0 * (x13 * x302 - x23 * x303) - x305 + x307
    x309 = x0 * (x24 * x302 - x30 * x304)
    x310 = x178 * x6
    x311 = x187 + x310
    x312 = x305 - x309 + x311
    x313 = x18 * x312
    x314 = -x313
    x315 = x308 * x7 + x314
    x316 = x315 * x6
    x317 = x0 * (x30 * x302 - x304 * x39)
    x318 = x185 * x6
    x319 = x198 + x318
    x320 = x309 - x317 + x319
    x321 = x18 * x320
    x322 = -x321
    x323 = x312 * x7 + x322
    x324 = x18 * x323
    x325 = 2 * x316 - 2 * x324
    x326 = x323 * x6
    x327 = x0 * (x302 * x39 - x304 * x52)
    x328 = x196 * x6
    x329 = x211 + x328
    x330 = x317 - x327 + x329
    x331 = x18 * x330
    x332 = -x331
    x333 = x320 * x7 + x332
    x334 = x18 * x333
    x335 = 2 * x326 - 2 * x334
    x336 = x0 * x323
    x337 = -x324
    x338 = x316 + x337
    x339 = -x334
    x340 = x326 + x339
    x341 = x0 * x315 - x18 * x340 - x336 + x338 * x6
    x342 = x209 * x6 + x226
    x343 = -x0 * (x302 * x52 - x304 * x67) + x327 + x342
    x344 = -x18 * x343
    x345 = -x18 * (x330 * x7 + x344)
    x346 = x333 * x6 + x345
    x347 = -x0 * x333 - x18 * x346 + x336 + x340 * x6
    x348 = x0 * (x174 + x22)
    x349 = x0 * (x178 + x26)
    x350 = x237 * x6
    x351 = x244 + x350
    x352 = x348 - x349 + x351
    x353 = x0 * (x185 + x32)
    x354 = x242 * x6
    x355 = x252 + x354
    x356 = x349 - x353 + x355
    x357 = x18 * x356
    x358 = -x357
    x359 = x352 * x7 + x358
    x360 = x359 * x6
    x361 = x0 * (x196 + x41)
    x362 = x250 * x6
    x363 = x264 + x362
    x364 = x353 - x361 + x363
    x365 = x18 * x364
    x366 = -x365
    x367 = x356 * x7 + x366
    x368 = x18 * x367
    x369 = 2 * x360 - 2 * x368
    x370 = x367 * x6
    x371 = x0 * (x209 + x54)
    x372 = x262 * x6
    x373 = x278 + x372
    x374 = x361 - x371 + x373
    x375 = x18 * x374
    x376 = -x375
    x377 = x364 * x7 + x376
    x378 = x18 * x377
    x379 = 2 * x370 - 2 * x378
    x380 = x0 * x367
    x381 = -x368
    x382 = x360 + x381
    x383 = -x378
    x384 = x370 + x383
    x385 = x0 * x359 - x18 * x384 - x380 + x382 * x6
    x386 = x0 * (x225 + x69)
    x387 = x276 * x6
    x388 = x295 + x387
    x389 = x18 * (x371 - x386 + x388)
    x390 = -x389
    x391 = -x18 * (x374 * x7 + x390)
    x392 = x377 * x6 + x391
    x393 = -x0 * x377 - x18 * x392 + x380 + x384 * x6
    x394 = 2 * x235
    x395 = 2 * x232 - 2 * x27
    x396 = 2 * x234 - x394 + x395
    x397 = x0 * x396
    x398 = 2 * x240
    x399 = 2 * x238 - 2 * x33
    x400 = x394 - x398 + x399
    x401 = x0 * x400
    x402 = x29 * x6
    x403 = x37 + x402
    x404 = x397 - x401 + x403
    x405 = 2 * x248
    x406 = 2 * x246 - 2 * x42
    x407 = x398 - x405 + x406
    x408 = x0 * x407
    x409 = x35 * x6
    x410 = x409 + x46
    x411 = x401 - x408 + x410
    x412 = x18 * x411
    x413 = -x412
    x414 = x404 * x7 + x413
    x415 = x414 * x6
    x416 = 2 * x260
    x417 = 2 * x258 - 2 * x55
    x418 = x405 - x416 + x417
    x419 = x0 * x418
    x420 = x44 * x6
    x421 = x420 + x59
    x422 = x408 - x419 + x421
    x423 = x18 * x422
    x424 = -x423
    x425 = x411 * x7 + x424
    x426 = x18 * x425
    x427 = 2 * x415 - 2 * x426
    x428 = x425 * x6
    x429 = 2 * x274
    x430 = 2 * x272 - 2 * x70
    x431 = x416 - x429 + x430
    x432 = x0 * x431
    x433 = x57 * x6
    x434 = x433 + x74
    x435 = x419 - x432 + x434
    x436 = x18 * x435
    x437 = -x436
    x438 = x422 * x7 + x437
    x439 = x18 * x438
    x440 = 2 * x428 - 2 * x439
    x441 = x0 * x425
    x442 = -x426
    x443 = x415 + x442
    x444 = -x439
    x445 = x428 + x444
    x446 = x0 * x414 - x18 * x445 - x441 + x443 * x6
    x447 = x0 * (-2 * x291 + 2 * x292 + x429 - 2 * x90)
    x448 = x6 * x72
    x449 = x448 + x93
    x450 = x18 * (x432 - x447 + x449)
    x451 = -x450
    x452 = -x18 * (x435 * x7 + x451)
    x453 = x438 * x6 + x452
    x454 = -x0 * x438 - x18 * x453 + x441 + x445 * x6
    x455 = 6 * x0
    x456 = x14 * x455
    x457 = x23 * x455
    x458 = x24 * x456
    x459 = x0 * (-3 * x112 + 3 * x176 + x457 - x458)
    x460 = x308 * x6
    x461 = x314 + x460
    x462 = x0 * (-3 * x108 + x13 * x456 + 3 * x172 - x457) - x459 + x461
    x463 = x462 * x6
    x464 = x30 * x456
    x465 = x0 * (-3 * x119 + 3 * x183 + x458 - x464)
    x466 = x312 * x6
    x467 = x322 + x466
    x468 = x459 - x465 + x467
    x469 = x18 * x468
    x470 = x468 * x6
    x471 = x39 * x456
    x472 = x0 * (-3 * x130 + 3 * x194 + x464 - x471)
    x473 = x320 * x6
    x474 = x332 + x473
    x475 = x465 - x472 + x474
    x476 = x18 * x475
    x477 = x0 * x468
    x478 = -x469
    x479 = x463 + x478
    x480 = -x476
    x481 = x470 + x480
    x482 = -x18 * x481
    x483 = x479 * x6 + x482
    x484 = x0 * x462 - x477 + x483
    x485 = x330 * x6 + x344
    x486 = -x18 * (-x0 * (-3 * x143 + 3 * x207 - x456 * x52 + x471) + x472 + x485)
    x487 = x475 * x6 + x486
    x488 = -x18 * x487
    x489 = x481 * x6 + x488
    x490 = -x0 * x475 + x477 + x489
    x491 = -x18 * x490
    x492 = x484 * x6 + x491
    x493 = x0 * (2 * x463 - 2 * x469) - x0 * (2 * x470 - 2 * x476) + x492
    x494 = x0 * (x308 + x396)
    x495 = x0 * (x312 + x400)
    x496 = x352 * x6
    x497 = x358 + x496
    x498 = x494 - x495 + x497
    x499 = x498 * x6
    x500 = x0 * (x320 + x407)
    x501 = x356 * x6
    x502 = x366 + x501
    x503 = x495 - x500 + x502
    x504 = x18 * x503
    x505 = x503 * x6
    x506 = x0 * (x330 + x418)
    x507 = x364 * x6
    x508 = x376 + x507
    x509 = x500 - x506 + x508
    x510 = x18 * x509
    x511 = x0 * x503
    x512 = -x504
    x513 = x499 + x512
    x514 = -x510
    x515 = x505 + x514
    x516 = -x18 * x515
    x517 = x513 * x6 + x516
    x518 = x0 * x498 - x511 + x517
    x519 = x0 * (x343 + x431)
    x520 = x374 * x6
    x521 = x390 + x520
    x522 = -x18 * (x506 - x519 + x521)
    x523 = x509 * x6 + x522
    x524 = -x18 * x523
    x525 = x515 * x6 + x524
    x526 = -x0 * x509 + x511 + x525
    x527 = -x18 * x526
    x528 = x518 * x6 + x527
    x529 = x0 * (2 * x499 - 2 * x504) - x0 * (2 * x505 - 2 * x510) + x528
    x530 = 2 * x349
    x531 = -2 * x243 + 2 * x350
    x532 = x0 * (x29 + 2 * x348 - x530 + x531)
    x533 = 2 * x353
    x534 = -2 * x251 + 2 * x354
    x535 = x0 * (x35 + x530 - x533 + x534)
    x536 = x404 * x6
    x537 = x413 + x536
    x538 = x532 - x535 + x537
    x539 = x538 * x6
    x540 = 2 * x361
    x541 = -2 * x263 + 2 * x362
    x542 = x0 * (x44 + x533 - x540 + x541)
    x543 = x411 * x6
    x544 = x424 + x543
    x545 = x535 - x542 + x544
    x546 = x18 * x545
    x547 = x545 * x6
    x548 = 2 * x371
    x549 = -2 * x277 + 2 * x372
    x550 = x0 * (x540 - x548 + x549 + x57)
    x551 = x422 * x6
    x552 = x437 + x551
    x553 = x542 - x550 + x552
    x554 = x18 * x553
    x555 = x0 * x545
    x556 = -x546
    x557 = x539 + x556
    x558 = -x554
    x559 = x547 + x558
    x560 = -x18 * x559
    x561 = x557 * x6 + x560
    x562 = x0 * x538 - x555 + x561
    x563 = x553 * x6
    x564 = x0 * (-2 * x294 - 2 * x386 + 2 * x387 + x548 + x72)
    x565 = x435 * x6
    x566 = x451 + x565
    x567 = x18 * (x550 - x564 + x566)
    x568 = -x567
    x569 = x563 + x568
    x570 = -x18 * x569
    x571 = x559 * x6 + x570
    x572 = -x0 * x553 + x555 + x571
    x573 = -x18 * x572
    x574 = x562 * x6 + x573
    x575 = x0 * (2 * x539 - 2 * x546) - x0 * (2 * x547 - 2 * x554) + x574
    x576 = 3 * x401
    x577 = x0 * (-3 * x36 + 3 * x397 + 3 * x402 - x576)
    x578 = 3 * x408
    x579 = x0 * (3 * x409 - 3 * x45 + x576 - x578)
    x580 = x38 * x6
    x581 = x49 + x580
    x582 = x577 - x579 + x581
    x583 = x582 * x6
    x584 = 3 * x419
    x585 = x0 * (3 * x420 + x578 - 3 * x58 - x584)
    x586 = x47 * x6
    x587 = x586 + x62
    x588 = x579 - x585 + x587
    x589 = x18 * x588
    x590 = x588 * x6
    x591 = 3 * x432
    x592 = x0 * (3 * x433 + x584 - x591 - 3 * x73)
    x593 = x6 * x60
    x594 = x593 + x77
    x595 = x585 - x592 + x594
    x596 = x18 * x595
    x597 = x0 * x588
    x598 = -x589
    x599 = x583 + x598
    x600 = -x596
    x601 = x590 + x600
    x602 = -x18 * x601
    x603 = x599 * x6 + x602
    x604 = x0 * x582 - x597 + x603
    x605 = x595 * x6
    x606 = x0 * (-3 * x447 + 3 * x448 + x591 - 3 * x92)
    x607 = x6 * x75
    x608 = x607 + x95
    x609 = x18 * (x592 - x606 + x608)
    x610 = -x609
    x611 = x605 + x610
    x612 = -x18 * x611
    x613 = x6 * x601 + x612
    x614 = -x0 * x595 + x597 + x613
    x615 = -x18 * x614
    x616 = x6 * x604 + x615
    x617 = x0 * (2 * x583 - 2 * x589) - x0 * (2 * x590 - 2 * x596) + x616
    x618 = x50 * x7 + x82
    x619 = x63 * x7 + x84
    x620 = x18 * x619
    x621 = -x620
    x622 = x6 * x618
    x623 = x0 * x83
    x624 = x0 * x85
    x625 = x623 - x624
    x626 = x0 * x97
    x627 = x6 * x619
    x628 = x18 * (x7 * x78 + x96)
    x629 = -x628
    x630 = (
        -x0 * (x100 + x619)
        + x0 * (x618 + x88)
        - x18 * (x624 - x626 + x627 + x629)
        + x6 * (x621 + x622 + x625)
    )
    x631 = x127 * x7 + x155
    x632 = x138 * x7 + x157
    x633 = x0 * x156
    x634 = x0 * x158
    x635 = x6 * x631
    x636 = x18 * x632
    x637 = -x636
    x638 = x0 * x167
    x639 = x6 * x632
    x640 = x18 * (x151 * x7 + x166)
    x641 = -x640
    x642 = x191 * x7 + x219
    x643 = x202 * x7 + x221
    x644 = x0 * x220
    x645 = x0 * x222
    x646 = x6 * x642
    x647 = x18 * x643
    x648 = -x647
    x649 = x0 * x230
    x650 = x6 * x643
    x651 = x18 * (x215 * x7 + x229)
    x652 = -x651
    x653 = x256 * x7 + x286
    x654 = x268 * x7 + x288
    x655 = x0 * x287
    x656 = x0 * x289
    x657 = x6 * x653
    x658 = x18 * x654
    x659 = -x658
    x660 = x0 * x299
    x661 = x6 * x654
    x662 = x18 * (x282 * x7 + x298)
    x663 = -x662
    x664 = x315 * x7 + x337
    x665 = x323 * x7 + x339
    x666 = x0 * x338
    x667 = x0 * x340
    x668 = x6 * x664
    x669 = x18 * x665
    x670 = -x669
    x671 = x0 * x346
    x672 = x6 * x665
    x673 = x18 * (x333 * x7 + x345)
    x674 = -x673
    x675 = x359 * x7 + x381
    x676 = x367 * x7 + x383
    x677 = x0 * x382
    x678 = x0 * x384
    x679 = x6 * x675
    x680 = x18 * x676
    x681 = -x680
    x682 = x0 * x392
    x683 = x6 * x676
    x684 = x18 * (x377 * x7 + x391)
    x685 = -x684
    x686 = x414 * x7 + x442
    x687 = x425 * x7 + x444
    x688 = x0 * x443
    x689 = x0 * x445
    x690 = x6 * x686
    x691 = x18 * x687
    x692 = -x691
    x693 = x0 * x453
    x694 = x6 * x687
    x695 = x18 * (x438 * x7 + x452)
    x696 = -x695
    x697 = x462 * x7 + x478
    x698 = x468 * x7 + x480
    x699 = x18 * x698
    x700 = -x699
    x701 = x6 * x697
    x702 = x0 * x479
    x703 = x0 * x481
    x704 = x702 - x703
    x705 = x0 * x487
    x706 = x6 * x698
    x707 = x18 * (x475 * x7 + x486)
    x708 = -x707
    x709 = x498 * x7 + x512
    x710 = x503 * x7 + x514
    x711 = x18 * x710
    x712 = -x711
    x713 = x6 * x709
    x714 = x0 * x513
    x715 = x0 * x515
    x716 = x714 - x715
    x717 = x0 * x523
    x718 = x6 * x710
    x719 = x18 * (x509 * x7 + x522)
    x720 = -x719
    x721 = x538 * x7 + x556
    x722 = x545 * x7 + x558
    x723 = x18 * x722
    x724 = -x723
    x725 = x6 * x721
    x726 = x0 * x557
    x727 = x0 * x559
    x728 = x726 - x727
    x729 = x0 * x569
    x730 = x6 * x722
    x731 = x18 * (x553 * x7 + x568)
    x732 = -x731
    x733 = x582 * x7 + x598
    x734 = x588 * x7 + x600
    x735 = x18 * x734
    x736 = -x735
    x737 = x6 * x733
    x738 = x0 * x599
    x739 = x0 * x601
    x740 = x738 - x739
    x741 = x0 * x611
    x742 = x6 * x734
    x743 = x18 * (x595 * x7 + x610)
    x744 = -x743
    x745 = 4 * x579
    x746 = x0 * (-4 * x48 + 4 * x577 + 4 * x580 - x745)
    x747 = 4 * x585
    x748 = x0 * (4 * x586 - 4 * x61 + x745 - x747)
    x749 = x746 - x748 + x83
    x750 = 4 * x592
    x751 = x0 * (4 * x593 + x747 - x750 - 4 * x76)
    x752 = x748 - x751 + x85
    x753 = -x18 * x752
    x754 = x6 * x749 + x753
    x755 = -x18 * (-x0 * (-4 * x606 + 4 * x607 + x750 - 4 * x94) + x751 + x97)
    x756 = x6 * x752 + x755
    x757 = x117 * x6
    x758 = x0 * x173
    x759 = x0 * x177
    x760 = 2 * x759
    x761 = x110 * x6
    x762 = -2 * x115 + 2 * x758 - x760 + 2 * x761
    x763 = x0 * x762
    x764 = x0 * x184
    x765 = 2 * x764
    x766 = x114 * x6
    x767 = -2 * x122 + x760 - x765 + 2 * x766
    x768 = x0 * x767
    x769 = 3 * x768
    x770 = x0 * (-3 * x125 + 3 * x757 + 3 * x763 - x769)
    x771 = x124 * x6
    x772 = x0 * x195
    x773 = 2 * x772
    x774 = x121 * x6
    x775 = -2 * x133 + x765 - x773 + 2 * x774
    x776 = x0 * x775
    x777 = 3 * x776
    x778 = x0 * (-3 * x136 + x769 + 3 * x771 - x777)
    x779 = x156 + x770 - x778
    x780 = x135 * x6
    x781 = x0 * x208
    x782 = 2 * x781
    x783 = x132 * x6
    x784 = -2 * x146 + x773 - x782 + 2 * x783
    x785 = x0 * x784
    x786 = 3 * x785
    x787 = x0 * (-3 * x149 + x777 + 3 * x780 - x786)
    x788 = x158 + x778 - x787
    x789 = -x18 * x788
    x790 = x6 * x779 + x789
    x791 = x148 * x6
    x792 = x0 * x224
    x793 = x145 * x6
    x794 = x0 * (-2 * x162 + x782 - 2 * x792 + 2 * x793)
    x795 = -x18 * (-x0 * (-3 * x164 + x786 + 3 * x791 - 3 * x794) + x167 + x787)
    x796 = x6 * x788 + x795
    x797 = x0 * x239
    x798 = 2 * x797
    x799 = x0 * x233
    x800 = -2 * x36 + 2 * x402
    x801 = -x798 + 2 * x799 + x800
    x802 = x0 * x801
    x803 = x0 * x247
    x804 = 2 * x803
    x805 = 2 * x409 - 2 * x45
    x806 = x798 - x804 + x805
    x807 = x0 * x806
    x808 = 3 * x807
    x809 = x0 * (-3 * x48 + 3 * x580 + 3 * x802 - x808)
    x810 = x0 * x259
    x811 = 2 * x810
    x812 = 2 * x420 - 2 * x58
    x813 = x804 - x811 + x812
    x814 = x0 * x813
    x815 = 3 * x814
    x816 = x0 * (3 * x586 - 3 * x61 + x808 - x815)
    x817 = x809 - x816 + x83
    x818 = x0 * x273
    x819 = 2 * x818
    x820 = 2 * x433 - 2 * x73
    x821 = x811 - x819 + x820
    x822 = x0 * x821
    x823 = 3 * x822
    x824 = x0 * (3 * x593 - 3 * x76 + x815 - x823)
    x825 = x816 - x824 + x85
    x826 = -x18 * x825
    x827 = x6 * x817 + x826
    x828 = x0 * x293
    x829 = x0 * (2 * x448 + x819 - 2 * x828 - 2 * x92)
    x830 = -x18 * (-x0 * (3 * x607 + x823 - 3 * x829 - 3 * x94) + x824 + x97)
    x831 = x6 * x825 + x830
    x832 = x0 * x307
    x833 = x0 * x311
    x834 = 2 * x833
    x835 = x181 * x6
    x836 = -2 * x189 + 2 * x832 - x834 + 2 * x835
    x837 = x0 * x836
    x838 = x0 * x319
    x839 = 2 * x838
    x840 = x188 * x6
    x841 = -2 * x200 + x834 - x839 + 2 * x840
    x842 = x0 * x841
    x843 = x220 + x837 - x842
    x844 = x0 * x329
    x845 = 2 * x844
    x846 = x199 * x6
    x847 = -2 * x213 + x839 - x845 + 2 * x846
    x848 = x0 * x847
    x849 = x222 + x842 - x848
    x850 = -x18 * x849
    x851 = x6 * x843 + x850
    x852 = x0 * x342
    x853 = x212 * x6
    x854 = -x18 * (-x0 * (-2 * x227 + x845 - 2 * x852 + 2 * x853) + x230 + x848)
    x855 = x6 * x849 + x854
    x856 = x0 * x351
    x857 = x0 * x355
    x858 = 2 * x857
    x859 = x245 * x6
    x860 = -2 * x254 + 2 * x856 - x858 + 2 * x859
    x861 = x0 * x860
    x862 = x0 * x363
    x863 = 2 * x862
    x864 = x253 * x6
    x865 = -2 * x266 + x858 - x863 + 2 * x864
    x866 = x0 * x865
    x867 = x287 + x861 - x866
    x868 = x0 * x373
    x869 = 2 * x868
    x870 = x265 * x6
    x871 = -2 * x280 + x863 - x869 + 2 * x870
    x872 = x0 * x871
    x873 = x289 + x866 - x872
    x874 = -x18 * x873
    x875 = x6 * x867 + x874
    x876 = x0 * x388
    x877 = x279 * x6
    x878 = -x18 * (-x0 * (-2 * x296 + x869 - 2 * x876 + 2 * x877) + x299 + x872)
    x879 = x6 * x873 + x878
    x880 = x0 * x410
    x881 = 2 * x880
    x882 = x0 * x403
    x883 = -2 * x48 + 2 * x580
    x884 = -x881 + 2 * x882 + x883
    x885 = x0 * x884
    x886 = x0 * x421
    x887 = 2 * x886
    x888 = 2 * x586 - 2 * x61
    x889 = x881 - x887 + x888
    x890 = x0 * x889
    x891 = x83 + x885 - x890
    x892 = x0 * x434
    x893 = 2 * x892
    x894 = 2 * x593 - 2 * x76
    x895 = x887 - x893 + x894
    x896 = x0 * x895
    x897 = x85 + x890 - x896
    x898 = -x18 * x897
    x899 = x6 * x891 + x898
    x900 = x0 * x449
    x901 = -x18 * (-x0 * (2 * x607 + x893 - 2 * x900 - 2 * x94) + x896 + x97)
    x902 = x6 * x897 + x901
    x903 = x0 * x461
    x904 = x0 * x467
    x905 = x338 + x903 - x904
    x906 = x0 * x474
    x907 = x340 + x904 - x906
    x908 = -x18 * x907
    x909 = x6 * x905 + x908
    x910 = -x18 * (-x0 * x485 + x346 + x906)
    x911 = x6 * x907 + x910
    x912 = x0 * x497
    x913 = x0 * x502
    x914 = x382 + x912 - x913
    x915 = x0 * x508
    x916 = x384 + x913 - x915
    x917 = -x18 * x916
    x918 = x6 * x914 + x917
    x919 = -x18 * (-x0 * x521 + x392 + x915)
    x920 = x6 * x916 + x919
    x921 = x0 * x537
    x922 = x0 * x544
    x923 = x443 + x921 - x922
    x924 = x0 * x552
    x925 = x445 + x922 - x924
    x926 = -x18 * x925
    x927 = x6 * x923 + x926
    x928 = -x18 * (-x0 * x566 + x453 + x924)
    x929 = x6 * x925 + x928
    x930 = x0 * x581
    x931 = x0 * x587
    x932 = x83 + x930 - x931
    x933 = x0 * x594
    x934 = x85 + x931 - x933
    x935 = -x18 * x934
    x936 = x6 * x932 + x935
    x937 = -x18 * (-x0 * x608 + x933 + x97)
    x938 = x6 * x934 + x937
    x939 = x618 * x7 + x621
    x940 = -x18 * (x619 * x7 + x629)
    x941 = 2 * x624
    x942 = (
        x0 * (-2 * x620 + 2 * x622 + 2 * x623 - x941)
        - x0 * (-2 * x626 + 2 * x627 - 2 * x628 + x941)
        + x6 * x939
        + x940
    )
    x943 = x631 * x7 + x637
    x944 = -x18 * (x632 * x7 + x641)
    x945 = 2 * x634
    x946 = x642 * x7 + x648
    x947 = -x18 * (x643 * x7 + x652)
    x948 = 2 * x645
    x949 = x653 * x7 + x659
    x950 = -x18 * (x654 * x7 + x663)
    x951 = 2 * x656
    x952 = x664 * x7 + x670
    x953 = -x18 * (x665 * x7 + x674)
    x954 = 2 * x667
    x955 = x675 * x7 + x681
    x956 = -x18 * (x676 * x7 + x685)
    x957 = 2 * x678
    x958 = x686 * x7 + x692
    x959 = -x18 * (x687 * x7 + x696)
    x960 = 2 * x689
    x961 = x697 * x7 + x700
    x962 = -x18 * (x698 * x7 + x708)
    x963 = 2 * x703
    x964 = x7 * x709 + x712
    x965 = -x18 * (x7 * x710 + x720)
    x966 = 2 * x715
    x967 = x7 * x721 + x724
    x968 = -x18 * (x7 * x722 + x732)
    x969 = 2 * x727
    x970 = x7 * x733 + x736
    x971 = -x18 * (x7 * x734 + x744)
    x972 = 2 * x739
    x973 = x7 * x749 + x753
    x974 = -x18 * (x7 * x752 + x755)
    x975 = x7 * x779 + x789
    x976 = -x18 * (x7 * x788 + x795)
    x977 = x7 * x817 + x826
    x978 = -x18 * (x7 * x825 + x830)
    x979 = x7 * x843 + x850
    x980 = -x18 * (x7 * x849 + x854)
    x981 = x7 * x867 + x874
    x982 = -x18 * (x7 * x873 + x878)
    x983 = x7 * x891 + x898
    x984 = -x18 * (x7 * x897 + x901)
    x985 = x7 * x905 + x908
    x986 = -x18 * (x7 * x907 + x910)
    x987 = x7 * x914 + x917
    x988 = -x18 * (x7 * x916 + x919)
    x989 = x7 * x923 + x926
    x990 = -x18 * (x7 * x925 + x928)
    x991 = x7 * x932 + x935
    x992 = -x18 * (x7 * x934 + x937)
    x993 = x479 * x7 + x482
    x994 = -x18 * (x481 * x7 + x488)
    x995 = x513 * x7 + x516
    x996 = -x18 * (x515 * x7 + x524)
    x997 = x557 * x7 + x560
    x998 = -x18 * (x559 * x7 + x570)
    x999 = x599 * x7 + x602
    x1000 = -x18 * (x601 * x7 + x612)
    x1001 = x7 * x83 + x86
    x1002 = -x18 * (x7 * x85 + x98)
    x1003 = 3 * x542
    x1004 = 3 * x535
    x1005 = x0 * (-x1003 + x1004 - 3 * x423 + x47 + 3 * x543)
    x1006 = 4 * x1005
    x1007 = x0 * (-x1004 + x38 - 3 * x412 + 3 * x532 + 3 * x536)
    x1008 = 3 * x550
    x1009 = x0 * (x1003 - x1008 - 3 * x436 + 3 * x551 + x60)
    x1010 = 4 * x1009
    x1011 = x0 * (x1006 - x1010 + 4 * x590 - 4 * x596 + x63)
    x1012 = x0 * (-x1006 + 4 * x1007 + x50 + 4 * x583 - 4 * x589) - x1011 + x754
    x1013 = x0 * (x1008 - 3 * x450 - 3 * x564 + 3 * x565 + x75)
    x1014 = -x18 * (-x0 * (x1010 - 4 * x1013 + 4 * x605 - 4 * x609 + x78) + x1011 + x756)
    x1015 = x1012 * x6 + x1014
    x1016 = -x18 * x195 + x184 * x6 + x249
    x1017 = x0 * (x1016 + x121)
    x1018 = 2 * x1017
    x1019 = x134 + x764 - x772 + x774
    x1020 = x1019 * x18
    x1021 = x177 * x6 - x18 * x184 + x241
    x1022 = x0 * (x1021 + x114)
    x1023 = 2 * x1022
    x1024 = x123 + x759 - x764 + x766
    x1025 = x1024 * x6
    x1026 = x0 * (-x1018 - 2 * x1020 + x1023 + 2 * x1025 + x124)
    x1027 = 3 * x1026
    x1028 = x137 + x768 + x771 - x776
    x1029 = x1028 * x18
    x1030 = x1024 * x18
    x1031 = x173 * x6 - x177 * x18 + x236
    x1032 = x0 * (x1031 + x110)
    x1033 = x6 * (x116 + x758 - x759 + x761)
    x1034 = x0 * (-x1023 - 2 * x1030 + 2 * x1032 + 2 * x1033 + x117)
    x1035 = x6 * (x126 + x757 + x763 - x768)
    x1036 = -x18 * x208 + x195 * x6 + x261
    x1037 = x0 * (x1036 + x132)
    x1038 = 2 * x1037
    x1039 = x147 + x772 - x781 + x783
    x1040 = x1039 * x18
    x1041 = x1019 * x6
    x1042 = x0 * (x1018 - x1038 - 2 * x1040 + 2 * x1041 + x135)
    x1043 = 3 * x1042
    x1044 = x150 + x776 + x780 - x785
    x1045 = x1044 * x18
    x1046 = x1028 * x6
    x1047 = x0 * (x1027 - x1043 - 3 * x1045 + 3 * x1046 + x138)
    x1048 = x0 * (-x1027 - 3 * x1029 + 3 * x1034 + 3 * x1035 + x127) - x1047 + x790
    x1049 = -x18 * x224 + x208 * x6 + x275
    x1050 = x0 * (x1049 + x145)
    x1051 = x18 * (x163 + x781 - x792 + x793)
    x1052 = x1039 * x6
    x1053 = x0 * (x1038 - 2 * x1050 - 2 * x1051 + 2 * x1052 + x148)
    x1054 = x18 * (x165 + x785 + x791 - x794)
    x1055 = x1044 * x6
    x1056 = -x18 * (
        -x0 * (x1043 - 3 * x1053 - 3 * x1054 + 3 * x1055 + x151) + x1047 + x796
    )
    x1057 = x1048 * x6 + x1056
    x1058 = x0 * x32
    x1059 = x0 * x41
    x1060 = x1058 - x1059 - x18 * x259 + x247 * x6
    x1061 = x0 * (x1060 + x44)
    x1062 = 2 * x1061
    x1063 = x421 + x803 - x810
    x1064 = x1063 * x18
    x1065 = x0 * x26
    x1066 = -x1058 + x1065 - x18 * x247 + x239 * x6
    x1067 = x0 * (x1066 + x35)
    x1068 = 2 * x1067
    x1069 = x410 + x797 - x803
    x1070 = x1069 * x6
    x1071 = x0 * (-x1062 - 2 * x1064 + x1068 + 2 * x1070 + x47)
    x1072 = 3 * x1071
    x1073 = x587 + x807 - x814
    x1074 = x1073 * x18
    x1075 = x1069 * x18
    x1076 = x0 * x22 - x1065 - x18 * x239 + x233 * x6
    x1077 = x0 * (x1076 + x29)
    x1078 = x6 * (x403 - x797 + x799)
    x1079 = x0 * (-x1068 - 2 * x1075 + 2 * x1077 + 2 * x1078 + x38)
    x1080 = x6 * (x581 + x802 - x807)
    x1081 = x0 * x54
    x1082 = x1059 - x1081 - x18 * x273 + x259 * x6
    x1083 = x0 * (x1082 + x57)
    x1084 = 2 * x1083
    x1085 = x434 + x810 - x818
    x1086 = x1085 * x18
    x1087 = x1063 * x6
    x1088 = x0 * (x1062 - x1084 - 2 * x1086 + 2 * x1087 + x60)
    x1089 = 3 * x1088
    x1090 = x594 + x814 - x822
    x1091 = x1090 * x18
    x1092 = x1073 * x6
    x1093 = x0 * (x1072 - x1089 - 3 * x1091 + 3 * x1092 + x63)
    x1094 = x0 * (-x1072 - 3 * x1074 + 3 * x1079 + 3 * x1080 + x50) - x1093 + x827
    x1095 = -x0 * x69 + x1081 - x18 * x293 + x273 * x6
    x1096 = x0 * (x1095 + x72)
    x1097 = x18 * (x449 + x818 - x828)
    x1098 = x1085 * x6
    x1099 = x0 * (x1084 - 2 * x1096 - 2 * x1097 + 2 * x1098 + x75)
    x1100 = x18 * (x608 + x822 - x829)
    x1101 = x1090 * x6
    x1102 = -x18 * (
        -x0 * (x1089 - 3 * x1099 - 3 * x1100 + 3 * x1101 + x78) + x1093 + x831
    )
    x1103 = x1094 * x6 + x1102
    x1104 = x0 * x178
    x1105 = x0 * x185
    x1106 = x1104 - x1105 - x18 * x319 + x311 * x6
    x1107 = x0 * (x1106 + x188)
    x1108 = 2 * x1107
    x1109 = x201 + x833 - x838 + x840
    x1110 = x1109 * x18
    x1111 = x0 * x174 - x1104 - x18 * x311 + x307 * x6
    x1112 = x0 * (x1111 + x181)
    x1113 = x6 * (x190 + x832 - x833 + x835)
    x1114 = x0 * x196
    x1115 = x1105 - x1114 - x18 * x329 + x319 * x6
    x1116 = x0 * (x1115 + x199)
    x1117 = 2 * x1116
    x1118 = x214 + x838 - x844 + x846
    x1119 = x1118 * x18
    x1120 = x1109 * x6
    x1121 = x0 * (x1108 - x1117 - 2 * x1119 + 2 * x1120 + x202)
    x1122 = x0 * (-x1108 - 2 * x1110 + 2 * x1112 + 2 * x1113 + x191) - x1121 + x851
    x1123 = -x0 * x209 + x1114 - x18 * x342 + x329 * x6
    x1124 = x0 * (x1123 + x212)
    x1125 = x18 * (x228 + x844 - x852 + x853)
    x1126 = x1118 * x6
    x1127 = -x18 * (
        -x0 * (x1117 - 2 * x1124 - 2 * x1125 + 2 * x1126 + x215) + x1121 + x855
    )
    x1128 = x1122 * x6 + x1127
    x1129 = x0 * x242
    x1130 = x0 * x250
    x1131 = x1129 - x1130 - x18 * x363 + x355 * x6
    x1132 = x0 * (x1131 + x253)
    x1133 = 2 * x1132
    x1134 = x267 + x857 - x862 + x864
    x1135 = x1134 * x18
    x1136 = x0 * x237 - x1129 - x18 * x355 + x351 * x6
    x1137 = x0 * (x1136 + x245)
    x1138 = x6 * (x255 + x856 - x857 + x859)
    x1139 = x0 * x262
    x1140 = x1130 - x1139 - x18 * x373 + x363 * x6
    x1141 = x0 * (x1140 + x265)
    x1142 = 2 * x1141
    x1143 = x281 + x862 - x868 + x870
    x1144 = x1143 * x18
    x1145 = x1134 * x6
    x1146 = x0 * (x1133 - x1142 - 2 * x1144 + 2 * x1145 + x268)
    x1147 = x0 * (-x1133 - 2 * x1135 + 2 * x1137 + 2 * x1138 + x256) - x1146 + x875
    x1148 = -x0 * x276 + x1139 - x18 * x388 + x373 * x6
    x1149 = x0 * (x1148 + x279)
    x1150 = x18 * (x297 + x868 - x876 + x877)
    x1151 = x1143 * x6
    x1152 = -x18 * (
        -x0 * (x1142 - 2 * x1149 - 2 * x1150 + 2 * x1151 + x282) + x1146 + x879
    )
    x1153 = x1147 * x6 + x1152
    x1154 = x0 * x35
    x1155 = x0 * x44
    x1156 = x1154 - x1155 - x18 * x421 + x410 * x6
    x1157 = x0 * (x1156 + x47)
    x1158 = 2 * x1157
    x1159 = x587 + x880 - x886
    x1160 = x1159 * x18
    x1161 = x0 * x29 - x1154 - x18 * x410 + x403 * x6
    x1162 = x0 * (x1161 + x38)
    x1163 = x6 * (x581 - x880 + x882)
    x1164 = x0 * x57
    x1165 = x1155 - x1164 - x18 * x434 + x421 * x6
    x1166 = x0 * (x1165 + x60)
    x1167 = 2 * x1166
    x1168 = x594 + x886 - x892
    x1169 = x1168 * x18
    x1170 = x1159 * x6
    x1171 = x0 * (x1158 - x1167 - 2 * x1169 + 2 * x1170 + x63)
    x1172 = x0 * (-x1158 - 2 * x1160 + 2 * x1162 + 2 * x1163 + x50) - x1171 + x899
    x1173 = -x0 * x72 + x1164 - x18 * x449 + x434 * x6
    x1174 = x0 * (x1173 + x75)
    x1175 = x18 * (x608 + x892 - x900)
    x1176 = x1168 * x6
    x1177 = -x18 * (
        -x0 * (x1167 - 2 * x1174 - 2 * x1175 + 2 * x1176 + x78) + x1171 + x902
    )
    x1178 = x1172 * x6 + x1177
    x1179 = x0 * x312
    x1180 = x0 * x308 - x1179 - x18 * x467 + x461 * x6
    x1181 = x0 * x320
    x1182 = x1179 - x1181 - x18 * x474 + x467 * x6
    x1183 = x0 * (x1182 + x323)
    x1184 = x0 * (x1180 + x315) - x1183 + x909
    x1185 = -x0 * x330 + x1181 - x18 * x485 + x474 * x6
    x1186 = -x18 * (-x0 * (x1185 + x333) + x1183 + x911)
    x1187 = x1184 * x6 + x1186
    x1188 = x0 * x356
    x1189 = x0 * x352 - x1188 - x18 * x502 + x497 * x6
    x1190 = x0 * x364
    x1191 = x1188 - x1190 - x18 * x508 + x502 * x6
    x1192 = x0 * (x1191 + x367)
    x1193 = x0 * (x1189 + x359) - x1192 + x918
    x1194 = -x0 * x374 + x1190 - x18 * x521 + x508 * x6
    x1195 = -x18 * (-x0 * (x1194 + x377) + x1192 + x920)
    x1196 = x1193 * x6 + x1195
    x1197 = x0 * x411
    x1198 = x0 * x404 - x1197 - x18 * x544 + x537 * x6
    x1199 = x0 * x422
    x1200 = x1197 - x1199 - x18 * x552 + x544 * x6
    x1201 = x0 * (x1200 + x425)
    x1202 = x0 * (x1198 + x414) - x1201 + x927
    x1203 = -x0 * x435 + x1199 - x18 * x566 + x552 * x6
    x1204 = -x18 * (-x0 * (x1203 + x438) + x1201 + x929)
    x1205 = x1202 * x6 + x1204
    x1206 = x0 * x47
    x1207 = x0 * x38 - x1206 - x18 * x587 + x581 * x6
    x1208 = x0 * x60
    x1209 = x1206 - x1208 - x18 * x594 + x587 * x6
    x1210 = x0 * (x1209 + x63)
    x1211 = x0 * (x1207 + x50) - x1210 + x936
    x1212 = -x0 * x75 + x1208 - x18 * x608 + x594 * x6
    x1213 = -x18 * (-x0 * (x1212 + x78) + x1210 + x938)
    x1214 = x1211 * x6 + x1213
    x1215 = x7 * x939 + x940
    x1216 = 2 * x419
    x1217 = 2 * x506
    x1218 = 2 * x408
    x1219 = 2 * x500
    x1220 = -2 * x375 + 2 * x507
    x1221 = 3 * x0
    x1222 = x1221 * (-x1216 - x1217 + x1218 + x1219 + x1220 + x812)
    x1223 = 2 * x585
    x1224 = 2 * x579
    x1225 = 2 * x401
    x1226 = 2 * x495
    x1227 = -2 * x365 + 2 * x501
    x1228 = x1221 * (-x1218 - x1219 + x1225 + x1226 + x1227 + x805)
    x1229 = 4 * x0
    x1230 = x1229 * (-x1222 - x1223 + x1224 + x1228 + 3 * x547 - 3 * x554 + x888)
    x1231 = x1005 - x1009 + x601
    x1232 = 2 * x748
    x1233 = -2 * x357 + 2 * x496
    x1234 = 2 * x768
    x1235 = 2 * x776
    x1236 = x0 * (-2 * x112 + 2 * x176)
    x1237 = x0 * (-2 * x119 + 2 * x183)
    x1238 = x169 * (-x1016 * x18 + x1021 * x6 + x1236 - x1237 + x767)
    x1239 = x0 * (-2 * x130 + 2 * x194)
    x1240 = x169 * (x1016 * x6 - x1036 * x18 + x1237 - x1239 + x775)
    x1241 = -x1017 - x1020 + x1022 + x1025
    x1242 = x1017 - x1037 - x1040 + x1041
    x1243 = x1221 * (
        x104 * x1241 + x1234 - x1235 + x1238 - x1240 - x1242 * x19 - 2 * x136 + 2 * x771
    )
    x1244 = x1026 - x1042 - x1045 + x1046
    x1245 = 3 * x18
    x1246 = 2 * x778
    x1247 = 3 * x6
    x1248 = 2 * x814
    x1249 = x0 * x406
    x1250 = x0 * x417
    x1251 = x169 * (x1060 * x6 - x1082 * x18 + x1249 - x1250 + x813)
    x1252 = x1061 - x1083 - x1086 + x1087
    x1253 = 2 * x807
    x1254 = x0 * x399
    x1255 = x169 * (-x1060 * x18 + x1066 * x6 - x1249 + x1254 + x806)
    x1256 = -x1061 - x1064 + x1067 + x1070
    x1257 = x1221 * (x104 * x1256 - x1248 - x1251 - x1252 * x19 + x1253 + x1255 + x888)
    x1258 = x1071 - x1088 - x1091 + x1092
    x1259 = 2 * x816
    x1260 = 2 * x842
    x1261 = x0 * (-2 * x186 + 2 * x310)
    x1262 = x0 * (-2 * x197 + 2 * x318)
    x1263 = x169 * (x1106 * x6 - x1115 * x18 + x1261 - x1262 + x841)
    x1264 = x1107 - x1116 - x1119 + x1120
    x1265 = 2 * x866
    x1266 = x0 * x534
    x1267 = x0 * x541
    x1268 = x169 * (x1131 * x6 - x1140 * x18 + x1266 - x1267 + x865)
    x1269 = x1132 - x1141 - x1144 + x1145
    x1270 = 2 * x890
    x1271 = x0 * x805
    x1272 = x0 * x812
    x1273 = x169 * (x1156 * x6 - x1165 * x18 + x1271 - x1272 + x889)
    x1274 = x1157 - x1166 - x1169 + x1170
    x1275 = x0 * (-2 * x321 + 2 * x466)
    x1276 = 2 * x904
    x1277 = x0 * x1227
    x1278 = 2 * x913
    x1279 = x0 * (-2 * x423 + 2 * x543)
    x1280 = 2 * x922
    x1281 = x0 * x888
    x1282 = 2 * x931

    # 150 item(s)
    S = numpy.array(
        [
            x103,
            x0 * x140 - x0 * x153 + x159 * x6 - x168 * x18,
            x103,
            x0 * x204 - x0 * x217 - x18 * x231 + x223 * x6,
            x0 * x270 - x0 * x284 - x18 * x300 + x290 * x6,
            x103,
            x0 * x325 - x0 * x335 - x18 * x347 + x341 * x6,
            x0 * x369 - x0 * x379 - x18 * x393 + x385 * x6,
            x0 * x427 - x0 * x440 - x18 * x454 + x446 * x6,
            x103,
            x493,
            x529,
            x575,
            x617,
            x103,
            x630,
            x0 * (x159 + x631)
            - x0 * (x168 + x632)
            - x18 * (x634 - x638 + x639 + x641)
            + x6 * (x633 - x634 + x635 + x637),
            x630,
            x0 * (x223 + x642)
            - x0 * (x231 + x643)
            - x18 * (x645 - x649 + x650 + x652)
            + x6 * (x644 - x645 + x646 + x648),
            x0 * (x290 + x653)
            - x0 * (x300 + x654)
            - x18 * (x656 - x660 + x661 + x663)
            + x6 * (x655 - x656 + x657 + x659),
            x630,
            x0 * (x341 + x664)
            - x0 * (x347 + x665)
            - x18 * (x667 - x671 + x672 + x674)
            + x6 * (x666 - x667 + x668 + x670),
            x0 * (x385 + x675)
            - x0 * (x393 + x676)
            - x18 * (x678 - x682 + x683 + x685)
            + x6 * (x677 - x678 + x679 + x681),
            x0 * (x446 + x686)
            - x0 * (x454 + x687)
            - x18 * (x689 - x693 + x694 + x696)
            + x6 * (x688 - x689 + x690 + x692),
            x630,
            x0 * (x484 + x697)
            - x0 * (x490 + x698)
            - x18 * (x703 - x705 + x706 + x708)
            + x6 * (x700 + x701 + x704),
            x0 * (x518 + x709)
            - x0 * (x526 + x710)
            - x18 * (x715 - x717 + x718 + x720)
            + x6 * (x712 + x713 + x716),
            x0 * (x562 + x721)
            - x0 * (x572 + x722)
            - x18 * (x727 - x729 + x730 + x732)
            + x6 * (x724 + x725 + x728),
            x0 * (x604 + x733)
            - x0 * (x614 + x734)
            - x18 * (x739 - x741 + x742 + x744)
            + x6 * (x736 + x737 + x740),
            x630,
            x0 * x749 - x0 * x752 - x18 * x756 + x6 * x754,
            x0 * x779 - x0 * x788 - x18 * x796 + x6 * x790,
            x0 * x817 - x0 * x825 - x18 * x831 + x6 * x827,
            x0 * x843 - x0 * x849 - x18 * x855 + x6 * x851,
            x0 * x867 - x0 * x873 - x18 * x879 + x6 * x875,
            x0 * x891 - x0 * x897 - x18 * x902 + x6 * x899,
            x0 * x905 - x0 * x907 - x18 * x911 + x6 * x909,
            x0 * x914 - x0 * x916 - x18 * x920 + x6 * x918,
            x0 * x923 - x0 * x925 - x18 * x929 + x6 * x927,
            x0 * x932 - x0 * x934 - x18 * x938 + x6 * x936,
            -x18 * x489 + x483 * x6 + x704,
            -x18 * x525 + x517 * x6 + x716,
            -x18 * x571 + x561 * x6 + x728,
            -x18 * x613 + x6 * x603 + x740,
            -x18 * x99 + x6 * x87 + x625,
            x942,
            x0 * (2 * x633 + 2 * x635 - 2 * x636 - x945)
            - x0 * (-2 * x638 + 2 * x639 - 2 * x640 + x945)
            + x6 * x943
            + x944,
            x942,
            x0 * (2 * x644 + 2 * x646 - 2 * x647 - x948)
            - x0 * (-2 * x649 + 2 * x650 - 2 * x651 + x948)
            + x6 * x946
            + x947,
            x0 * (2 * x655 + 2 * x657 - 2 * x658 - x951)
            - x0 * (-2 * x660 + 2 * x661 - 2 * x662 + x951)
            + x6 * x949
            + x950,
            x942,
            x0 * (2 * x666 + 2 * x668 - 2 * x669 - x954)
            - x0 * (-2 * x671 + 2 * x672 - 2 * x673 + x954)
            + x6 * x952
            + x953,
            x0 * (2 * x677 + 2 * x679 - 2 * x680 - x957)
            - x0 * (-2 * x682 + 2 * x683 - 2 * x684 + x957)
            + x6 * x955
            + x956,
            x0 * (2 * x688 + 2 * x690 - 2 * x691 - x960)
            - x0 * (-2 * x693 + 2 * x694 - 2 * x695 + x960)
            + x6 * x958
            + x959,
            x942,
            x0 * (-2 * x699 + 2 * x701 + 2 * x702 - x963)
            - x0 * (-2 * x705 + 2 * x706 - 2 * x707 + x963)
            + x6 * x961
            + x962,
            x0 * (-2 * x711 + 2 * x713 + 2 * x714 - x966)
            - x0 * (-2 * x717 + 2 * x718 - 2 * x719 + x966)
            + x6 * x964
            + x965,
            x0 * (-2 * x723 + 2 * x725 + 2 * x726 - x969)
            - x0 * (-2 * x729 + 2 * x730 - 2 * x731 + x969)
            + x6 * x967
            + x968,
            x0 * (-2 * x735 + 2 * x737 + 2 * x738 - x972)
            - x0 * (-2 * x741 + 2 * x742 - 2 * x743 + x972)
            + x6 * x970
            + x971,
            x942,
            x0 * x754 - x0 * x756 + x6 * x973 + x974,
            x0 * x790 - x0 * x796 + x6 * x975 + x976,
            x0 * x827 - x0 * x831 + x6 * x977 + x978,
            x0 * x851 - x0 * x855 + x6 * x979 + x980,
            x0 * x875 - x0 * x879 + x6 * x981 + x982,
            x0 * x899 - x0 * x902 + x6 * x983 + x984,
            x0 * x909 - x0 * x911 + x6 * x985 + x986,
            x0 * x918 - x0 * x920 + x6 * x987 + x988,
            x0 * x927 - x0 * x929 + x6 * x989 + x990,
            x0 * x936 - x0 * x938 + x6 * x991 + x992,
            x0 * x483 - x0 * x489 + x6 * x993 + x994,
            x0 * x517 - x0 * x525 + x6 * x995 + x996,
            x0 * x561 - x0 * x571 + x6 * x997 + x998,
            x0 * x603 - x0 * x613 + x1000 + x6 * x999,
            x0 * x87 - x0 * x99 + x1001 * x6 + x1002,
            x1015,
            x1057,
            x1103,
            x1128,
            x1153,
            x1178,
            x1187,
            x1196,
            x1205,
            x1214,
            x492,
            x528,
            x574,
            x616,
            x102,
            x1215,
            x7 * x943 + x944,
            x1215,
            x7 * x946 + x947,
            x7 * x949 + x950,
            x1215,
            x7 * x952 + x953,
            x7 * x955 + x956,
            x7 * x958 + x959,
            x1215,
            x7 * x961 + x962,
            x7 * x964 + x965,
            x7 * x967 + x968,
            x7 * x970 + x971,
            x1215,
            x7 * x973 + x974,
            x7 * x975 + x976,
            x7 * x977 + x978,
            x7 * x979 + x980,
            x7 * x981 + x982,
            x7 * x983 + x984,
            x7 * x985 + x986,
            x7 * x987 + x988,
            x7 * x989 + x990,
            x7 * x991 + x992,
            x7 * x993 + x994,
            x7 * x995 + x996,
            x7 * x997 + x998,
            x1000 + x7 * x999,
            x1001 * x7 + x1002,
            x1012 * x7 + x1014,
            x1048 * x7 + x1056,
            x1094 * x7 + x1102,
            x1122 * x7 + x1127,
            x1147 * x7 + x1152,
            x1172 * x7 + x1177,
            x1184 * x7 + x1186,
            x1193 * x7 + x1195,
            x1202 * x7 + x1204,
            x1211 * x7 + x1213,
            x484 * x7 + x491,
            x518 * x7 + x527,
            x562 * x7 + x573,
            x604 * x7 + x615,
            x101 + x7 * x88,
            x0
            * (
                x1229
                * (
                    x1221 * (-x1225 - x1226 + x1233 + 2 * x397 + 2 * x494 + x800)
                    - x1224
                    - x1228
                    + 3 * x539
                    - 3 * x546
                    + 2 * x577
                    + x883
                )
                - x1230
                - x1231 * x303
                - x1232
                + x301 * (-x1005 + x1007 + x599)
                + x65
                + 2 * x746
            )
            - x0
            * (
                -x1229
                * (
                    -x1221
                    * (x1216 + x1217 - 2 * x389 - 2 * x432 - 2 * x519 + 2 * x520 + x820)
                    + x1222
                    + x1223
                    + 3 * x563
                    - 3 * x567
                    - 2 * x592
                    + x894
                )
                + x1230
                + x1231 * x301
                + x1232
                - x303 * (x1009 - x1013 + x611)
                - 2 * x751
                + x80
            )
            + x1015,
            -x0
            * (
                -x1221
                * (
                    x104 * x1242
                    + x1235
                    + x1240
                    - 2 * x149
                    - x169
                    * (
                        -x0 * (-2 * x143 + 2 * x207)
                        + x1036 * x6
                        - x1049 * x18
                        + x1239
                        + x784
                    )
                    - x19 * (x1037 - x1050 - x1051 + x1052)
                    + 2 * x780
                    - 2 * x785
                )
                + x1243
                + x1244 * x1247
                - x1245 * (x1042 - x1053 - x1054 + x1055)
                + x1246
                + x153
                - 2 * x787
            )
            + x0
            * (
                x1221
                * (
                    x104 * (-x1022 - x1030 + x1032 + x1033)
                    - x1234
                    - x1238
                    - x1241 * x19
                    - 2 * x125
                    + x169
                    * (
                        x0 * (-2 * x108 + 2 * x172)
                        - x1021 * x18
                        + x1031 * x6
                        - x1236
                        + x762
                    )
                    + 2 * x757
                    + 2 * x763
                )
                - x1243
                - x1244 * x1245
                - x1246
                + x1247 * (-x1026 - x1029 + x1034 + x1035)
                + x140
                + 2 * x770
            )
            + x1057,
            -x0
            * (
                -x1221
                * (
                    x104 * x1252
                    + x1248
                    + x1251
                    - x169 * (-x0 * x430 + x1082 * x6 - x1095 * x18 + x1250 + x821)
                    - x19 * (x1083 - x1096 - x1097 + x1098)
                    - 2 * x822
                    + x894
                )
                - x1245 * (x1088 - x1099 - x1100 + x1101)
                + x1247 * x1258
                + x1257
                + x1259
                + x80
                - 2 * x824
            )
            + x0
            * (
                x1221
                * (
                    x104 * (-x1067 - x1075 + x1077 + x1078)
                    - x1253
                    - x1255
                    - x1256 * x19
                    + x169 * (x0 * x395 - x1066 * x18 + x1076 * x6 - x1254 + x801)
                    + 2 * x802
                    + x883
                )
                - x1245 * x1258
                + x1247 * (-x1071 - x1074 + x1079 + x1080)
                - x1257
                - x1259
                + x65
                + 2 * x809
            )
            + x1103,
            -x0
            * (
                x104 * x1264
                + x1260
                + x1263
                - x169
                * (-x0 * (-2 * x210 + 2 * x328) + x1115 * x6 - x1123 * x18 + x1262 + x847)
                - x19 * (x1116 - x1124 - x1125 + x1126)
                + x217
                - 2 * x848
            )
            + x0
            * (
                x104 * (-x1107 - x1110 + x1112 + x1113)
                - x1260
                - x1263
                - x1264 * x19
                + x169
                * (x0 * (-2 * x179 + 2 * x306) - x1106 * x18 + x1111 * x6 - x1261 + x836)
                + x204
                + 2 * x837
            )
            + x1128,
            -x0
            * (
                x104 * x1269
                + x1265
                + x1268
                - x169 * (-x0 * x549 + x1140 * x6 - x1148 * x18 + x1267 + x871)
                - x19 * (x1141 - x1149 - x1150 + x1151)
                + x284
                - 2 * x872
            )
            + x0
            * (
                x104 * (-x1132 - x1135 + x1137 + x1138)
                - x1265
                - x1268
                - x1269 * x19
                + x169 * (x0 * x531 - x1131 * x18 + x1136 * x6 - x1266 + x860)
                + x270
                + 2 * x861
            )
            + x1153,
            -x0
            * (
                x104 * x1274
                + x1270
                + x1273
                - x169 * (-x0 * x820 + x1165 * x6 - x1173 * x18 + x1272 + x895)
                - x19 * (x1166 - x1174 - x1175 + x1176)
                + x80
                - 2 * x896
            )
            + x0
            * (
                x104 * (-x1157 - x1160 + x1162 + x1163)
                - x1270
                - x1273
                - x1274 * x19
                + x169 * (x0 * x800 - x1156 * x18 + x1161 * x6 - x1271 + x884)
                + x65
                + 2 * x885
            )
            + x1178,
            x0
            * (
                x0 * (-2 * x313 + 2 * x460)
                + x1180 * x6
                - x1182 * x18
                - x1275
                - x1276
                + x325
                + 2 * x903
            )
            - x0
            * (
                -x0 * (-2 * x331 + 2 * x473)
                + x1182 * x6
                - x1185 * x18
                + x1275
                + x1276
                + x335
                - 2 * x906
            )
            + x1187,
            -x0
            * (-x0 * x1220 + x1191 * x6 - x1194 * x18 + x1277 + x1278 + x379 - 2 * x915)
            + x0
            * (x0 * x1233 + x1189 * x6 - x1191 * x18 - x1277 - x1278 + x369 + 2 * x912)
            + x1196,
            x0
            * (
                x0 * (-2 * x412 + 2 * x536)
                + x1198 * x6
                - x1200 * x18
                - x1279
                - x1280
                + x427
                + 2 * x921
            )
            - x0
            * (
                -x0 * (-2 * x436 + 2 * x551)
                + x1200 * x6
                - x1203 * x18
                + x1279
                + x1280
                + x440
                - 2 * x924
            )
            + x1205,
            x0 * (x0 * x883 + x1207 * x6 - x1209 * x18 - x1281 - x1282 + x65 + 2 * x930)
            - x0
            * (-x0 * x894 + x1209 * x6 - x1212 * x18 + x1281 + x1282 + x80 - 2 * x933)
            + x1214,
            x493,
            x529,
            x575,
            x617,
            x103,
        ]
    )
    return S


def coulomb_40(a, A, b, B, C):
    """Cartesian (gs) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = 6 * x0
    x2 = a + b
    x3 = x2 ** (-1.0)
    x4 = -x3 * (a * A[0] + b * B[0])
    x5 = x4 + C[0]
    x6 = -x3 * (a * A[1] + b * B[1])
    x7 = x6 + C[1]
    x8 = -x3 * (a * A[2] + b * B[2])
    x9 = x8 + C[2]
    x10 = x5 ** 2 + x7 ** 2 + x9 ** 2
    x11 = x10 * x2
    x12 = boys(0, x11)
    x13 = numpy.pi * x3 * numpy.exp(-a * b * x10 * x3)
    x14 = x12 * x13
    x15 = boys(1, x11)
    x16 = x13 * x15
    x17 = x1 * x16
    x18 = numpy.sqrt(abs(x4 + A[0]) ** 2 + abs(x6 + A[1]) ** 2 + abs(x8 + A[2]) ** 2)
    x19 = 2 * x13
    x20 = x18 * x19
    x21 = numpy.sqrt(abs(x5) ** 2 + abs(x7) ** 2 + abs(x9) ** 2)
    x22 = x19 * x21
    x23 = -x15 * x22
    x24 = x12 * x20 + x23
    x25 = x18 * x24
    x26 = boys(2, x11)
    x27 = -x22 * x26
    x28 = x15 * x20 + x27
    x29 = x21 * x28
    x30 = x13 * x26
    x31 = x18 * x28
    x32 = boys(3, x11)
    x33 = -x22 * x32
    x34 = x20 * x26 + x33
    x35 = x21 * x34
    x36 = 4 * x18
    x37 = 4 * x21
    x38 = x0 * (x16 * x36 - x30 * x37)
    x39 = x0 * x19
    x40 = x15 * x39
    x41 = -x29
    x42 = x25 + x41
    x43 = x12 * x39 - x40 + x42
    x44 = x26 * x39
    x45 = -x35
    x46 = x31 + x45
    x47 = x40 - x44 + x46
    x48 = -x21 * x47
    x49 = x18 * x43 + x48
    x50 = x0 * (x14 * x36 - x16 * x37) - x38 + x49
    x51 = -x22 * boys(4, x11)
    x52 = x20 * x32 + x51
    x53 = -x21 * x52
    x54 = x18 * x34 + x53
    x55 = -x32 * x39 + x44 + x54
    x56 = -x21 * x55
    x57 = x18 * x47 + x56
    x58 = -x0 * (-x13 * x32 * x37 + x30 * x36) + x38 + x57
    x59 = -x21 * x58
    x60 = x18 * x50 + x59
    x61 = (
        x0 * (x1 * x14 - x17 + 3 * x25 - 3 * x29)
        - x0 * (-x1 * x30 + x17 + 3 * x31 - 3 * x35)
        + x60
    )
    x62 = x0 * x24
    x63 = x0 * x28
    x64 = 2 * x63
    x65 = numpy.sqrt(abs(x4 + B[0]) ** 2 + abs(x6 + B[1]) ** 2 + abs(x8 + B[2]) ** 2)
    x66 = x19 * x65
    x67 = x12 * x66 + x23
    x68 = x18 * x67
    x69 = x15 * x66 + x27
    x70 = x21 * x69
    x71 = 2 * x62 - x64 + 2 * x68 - 2 * x70
    x72 = x0 * x34
    x73 = 2 * x72
    x74 = x18 * x69
    x75 = x26 * x66 + x33
    x76 = x21 * x75
    x77 = x64 - x73 + 2 * x74 - 2 * x76
    x78 = x0 * (x43 + x67)
    x79 = x0 * (x47 + x69)
    x80 = -x70
    x81 = x62 - x63
    x82 = x18 * (x68 + x80 + x81)
    x83 = -x76
    x84 = x63 - x72
    x85 = x74 + x83 + x84
    x86 = x21 * x85
    x87 = x0 * (x55 + x75)
    x88 = x18 * x85
    x89 = x0 * x52
    x90 = x18 * x75
    x91 = x21 * (x32 * x66 + x51)
    x92 = -x91
    x93 = x21 * (x72 - x89 + x90 + x92)
    x94 = x18 * x42 - x21 * x46 + x81
    x95 = x18 * x46 - x21 * x54 + x84
    x96 = x65 * x67 + x80
    x97 = x18 * x96
    x98 = x65 * x69 + x83
    x99 = x21 * x98
    x100 = -x99
    x101 = x0 * x71
    x102 = x0 * x77
    x103 = x18 * x98
    x104 = x21 * (x65 * x75 + x92)
    x105 = -x104
    x106 = x0 * (x73 - 2 * x89 + 2 * x90 - 2 * x91)
    x107 = 2 * x79
    x108 = x24 * x65 + x41
    x109 = x28 * x65 + x45
    x110 = x0 * x42
    x111 = x0 * x46
    x112 = x108 * x18
    x113 = x109 * x21
    x114 = -x113
    x115 = x0 * x54
    x116 = x109 * x18
    x117 = x21 * (x34 * x65 + x53)
    x118 = -x117
    x119 = x100 + x65 * x96
    x120 = -x21 * (x105 + x65 * x98)
    x121 = 3 * x102
    x122 = x108 * x65 + x114
    x123 = -x21 * (x109 * x65 + x118)
    x124 = 2 * x111
    x125 = x43 * x65 + x48
    x126 = -x21 * (x47 * x65 + x56)

    # 15 item(s)
    S = numpy.array(
        [
            x61,
            x0 * (x50 + x71)
            - x0 * (x58 + x77)
            + x18 * (x78 - x79 + x82 - x86)
            - x21 * (x79 - x87 + x88 - x93),
            x0 * (2 * x25 - 2 * x29) - x0 * (2 * x31 - 2 * x35) + x18 * x94 - x21 * x95,
            x0 * (-x107 + 2 * x78 + 2 * x82 - 2 * x86 + x96)
            - x0 * (x107 - 2 * x87 + 2 * x88 - 2 * x93 + x98)
            + x18 * (x100 + x101 - x102 + x97)
            - x21 * (x102 + x103 + x105 - x106),
            x0 * (x108 + x94)
            - x0 * (x109 + x95)
            + x18 * (x110 - x111 + x112 + x114)
            - x21 * (x111 - x115 + x116 + x118),
            x0 * x43 - x0 * x47 + x18 * x49 - x21 * x57,
            x0 * (3 * x101 - x121 + 3 * x97 - 3 * x99)
            - x0 * (3 * x103 - 3 * x104 - 3 * x106 + x121)
            + x119 * x18
            + x120,
            x0 * (2 * x110 + 2 * x112 - 2 * x113 - x124)
            - x0 * (-2 * x115 + 2 * x116 - 2 * x117 + x124)
            + x122 * x18
            + x123,
            x0 * x49 - x0 * x57 + x125 * x18 + x126,
            x60,
            x119 * x65 + x120,
            x122 * x65 + x123,
            x125 * x65 + x126,
            x50 * x65 + x59,
            x61,
        ]
    )
    return S


def coulomb_41(a, A, b, B, C):
    """Cartesian (gp) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = x3 + C[0]
    x5 = -x2 * (a * A[1] + b * B[1])
    x6 = x5 + C[1]
    x7 = -x2 * (a * A[2] + b * B[2])
    x8 = x7 + C[2]
    x9 = x4 ** 2 + x6 ** 2 + x8 ** 2
    x10 = x1 * x9
    x11 = boys(0, x10)
    x12 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x5 + B[1]) ** 2 + abs(x7 + B[2]) ** 2)
    x13 = numpy.pi * x2 * numpy.exp(-a * b * x2 * x9)
    x14 = 2 * x13
    x15 = x12 * x14
    x16 = boys(1, x10)
    x17 = numpy.sqrt(abs(x4) ** 2 + abs(x6) ** 2 + abs(x8) ** 2)
    x18 = x14 * x17
    x19 = -x16 * x18
    x20 = x11 * x15 + x19
    x21 = x0 * x20
    x22 = boys(2, x10)
    x23 = -x18 * x22
    x24 = x15 * x16 + x23
    x25 = x0 * x24
    x26 = 3 * x25
    x27 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x5 + A[1]) ** 2 + abs(x7 + A[2]) ** 2)
    x28 = x20 * x27
    x29 = x17 * x24
    x30 = -x29
    x31 = x28 + x30
    x32 = x27 * x31
    x33 = x24 * x27
    x34 = boys(3, x10)
    x35 = -x18 * x34
    x36 = x15 * x22 + x35
    x37 = x17 * x36
    x38 = -x37
    x39 = x33 + x38
    x40 = x17 * x39
    x41 = x0 * x36
    x42 = x27 * x39
    x43 = x27 * x36
    x44 = boys(4, x10)
    x45 = -x18 * x44
    x46 = x15 * x34 + x45
    x47 = x17 * x46
    x48 = -x47
    x49 = x43 + x48
    x50 = x17 * x49
    x51 = 2 * x28 - 2 * x29
    x52 = 2 * x33 - 2 * x37
    x53 = x0 * x52
    x54 = -x40
    x55 = x32 + x54
    x56 = x21 - x25 + x55
    x57 = -x50
    x58 = x42 + x57
    x59 = x25 - x41 + x58
    x60 = -x17 * x59
    x61 = x27 * x56 + x60
    x62 = x0 * x51 - x53 + x61
    x63 = 2 * x43 - 2 * x47
    x64 = -x18 * boys(5, x10)
    x65 = -x17 * (x15 * x44 + x64)
    x66 = x27 * x46 + x65
    x67 = -x17 * x66
    x68 = x27 * x49 + x67
    x69 = -x0 * x46 + x41 + x68
    x70 = -x17 * x69
    x71 = x27 * x59 + x70
    x72 = -x0 * x63 + x53 + x71
    x73 = -x17 * x72
    x74 = x27 * x62 + x73
    x75 = (
        x0 * (3 * x21 - x26 + 3 * x32 - 3 * x40)
        - x0 * (x26 - 3 * x41 + 3 * x42 - 3 * x50)
        + x74
    )
    x76 = x14 * x27
    x77 = x11 * x76 + x19
    x78 = x0 * x77
    x79 = x16 * x76 + x23
    x80 = x0 * x79
    x81 = 3 * x80
    x82 = x27 * x77
    x83 = x17 * x79
    x84 = -x83
    x85 = x82 + x84
    x86 = x27 * x85
    x87 = x27 * x79
    x88 = x22 * x76 + x35
    x89 = x17 * x88
    x90 = -x89
    x91 = x87 + x90
    x92 = x17 * x91
    x93 = x0 * x88
    x94 = x27 * x91
    x95 = x27 * x88
    x96 = x34 * x76 + x45
    x97 = x17 * x96
    x98 = -x97
    x99 = x95 + x98
    x100 = x17 * x99
    x101 = x0 * (2 * x87 - 2 * x89)
    x102 = -x92
    x103 = x102 + x86
    x104 = x78 - x80
    x105 = x103 + x104
    x106 = -x100
    x107 = x106 + x94
    x108 = x80 - x93
    x109 = x107 + x108
    x110 = -x109 * x17
    x111 = x105 * x27 + x110
    x112 = x0 * (2 * x82 - 2 * x83) - x101 + x111
    x113 = x44 * x76 + x64
    x114 = -x113 * x17
    x115 = x114 + x27 * x96
    x116 = -x115 * x17
    x117 = x116 + x27 * x99
    x118 = x0 * x96
    x119 = -x118 + x93
    x120 = x117 + x119
    x121 = -x120 * x17
    x122 = x109 * x27 + x121
    x123 = -x0 * (2 * x95 - 2 * x97) + x101 + x122
    x124 = -x123 * x17
    x125 = x112 * x27 + x124
    x126 = (
        -x0 * (-3 * x100 + x81 - 3 * x93 + 3 * x94)
        + x0 * (3 * x78 - x81 + 3 * x86 - 3 * x92)
        + x125
    )
    x127 = x0 * x31
    x128 = x0 * x39
    x129 = 2 * x128
    x130 = x12 * x20 + x30
    x131 = x130 * x27
    x132 = x12 * x24 + x38
    x133 = x132 * x17
    x134 = 2 * x127 - x129 + 2 * x131 - 2 * x133
    x135 = x0 * x49
    x136 = 2 * x135
    x137 = x132 * x27
    x138 = x12 * x36 + x48
    x139 = x138 * x17
    x140 = x129 - x136 + 2 * x137 - 2 * x139
    x141 = x0 * (x130 + x56)
    x142 = x0 * (x132 + x59)
    x143 = -x133
    x144 = x127 - x128
    x145 = x27 * (x131 + x143 + x144)
    x146 = -x139
    x147 = x128 - x135
    x148 = x137 + x146 + x147
    x149 = x148 * x17
    x150 = x0 * (x138 + x69)
    x151 = x148 * x27
    x152 = x0 * x66
    x153 = x138 * x27
    x154 = x17 * (x12 * x46 + x65)
    x155 = -x154
    x156 = x17 * (x135 - x152 + x153 + x155)
    x157 = (
        x0 * (x134 + x62)
        - x0 * (x140 + x72)
        - x17 * (x142 - x150 + x151 - x156)
        + x27 * (x141 - x142 + x145 - x149)
    )
    x158 = x0 * x85
    x159 = x0 * x91
    x160 = 2 * x159
    x161 = x12 * x77 + x84
    x162 = x161 * x27
    x163 = x12 * x79 + x90
    x164 = x163 * x17
    x165 = 2 * x158 - x160 + 2 * x162 - 2 * x164
    x166 = x0 * x99
    x167 = 2 * x166
    x168 = x163 * x27
    x169 = x12 * x88 + x98
    x170 = x169 * x17
    x171 = x160 - x167 + 2 * x168 - 2 * x170
    x172 = x0 * (x105 + x161)
    x173 = x0 * (x109 + x163)
    x174 = -x164
    x175 = x158 - x159
    x176 = x27 * (x162 + x174 + x175)
    x177 = -x170
    x178 = x159 - x166
    x179 = x168 + x177 + x178
    x180 = x17 * x179
    x181 = x0 * (x120 + x169)
    x182 = x179 * x27
    x183 = x0 * x115
    x184 = x169 * x27
    x185 = x17 * (x114 + x12 * x96)
    x186 = -x185
    x187 = x17 * (x166 - x183 + x184 + x186)
    x188 = x104 + x31
    x189 = x188 * x27
    x190 = x108 + x39
    x191 = x17 * x190
    x192 = x190 * x27
    x193 = x119 + x49
    x194 = x17 * x193
    x195 = x0 * x190
    x196 = -x191
    x197 = x189 + x196
    x198 = -x194
    x199 = x192 + x198
    x200 = x0 * x188 - x17 * x199 - x195 + x197 * x27
    x201 = -x17 * (-x0 * x113 + x118 + x66)
    x202 = x193 * x27 + x201
    x203 = -x0 * x193 - x17 * x202 + x195 + x199 * x27
    x204 = x103 * x27 - x107 * x17 + x175
    x205 = x107 * x27 - x117 * x17 + x178
    x206 = x144 - x17 * x58 + x27 * x55
    x207 = x147 - x17 * x68 + x27 * x58
    x208 = x12 * x130 + x143
    x209 = x208 * x27
    x210 = x12 * x132 + x146
    x211 = x17 * x210
    x212 = -x211
    x213 = x0 * x134
    x214 = x0 * x140
    x215 = x210 * x27
    x216 = x17 * (x12 * x138 + x155)
    x217 = -x216
    x218 = x0 * (x136 - 2 * x152 + 2 * x153 - 2 * x154)
    x219 = 2 * x142
    x220 = (
        x0 * (2 * x141 + 2 * x145 - 2 * x149 + x208 - x219)
        - x0 * (-2 * x150 + 2 * x151 - 2 * x156 + x210 + x219)
        - x17 * (x214 + x215 + x217 - x218)
        + x27 * (x209 + x212 + x213 - x214)
    )
    x221 = x12 * x161 + x174
    x222 = x221 * x27
    x223 = x12 * x163 + x177
    x224 = x17 * x223
    x225 = -x224
    x226 = x0 * x165
    x227 = x0 * x171
    x228 = x223 * x27
    x229 = x17 * (x12 * x169 + x186)
    x230 = -x229
    x231 = x0 * (x167 - 2 * x183 + 2 * x184 - 2 * x185)
    x232 = 2 * x173
    x233 = x12 * x188 + x196
    x234 = x12 * x190 + x198
    x235 = x0 * x197
    x236 = x0 * x199
    x237 = x233 * x27
    x238 = x17 * x234
    x239 = -x238
    x240 = x0 * x202
    x241 = x234 * x27
    x242 = x17 * (x12 * x193 + x201)
    x243 = -x242
    x244 = x102 + x12 * x85
    x245 = x106 + x12 * x91
    x246 = x0 * x103
    x247 = x0 * x107
    x248 = x244 * x27
    x249 = x17 * x245
    x250 = -x249
    x251 = x0 * x117
    x252 = x245 * x27
    x253 = x17 * (x116 + x12 * x99)
    x254 = -x253
    x255 = x12 * x31 + x54
    x256 = x12 * x39 + x57
    x257 = x0 * x55
    x258 = x0 * x58
    x259 = x255 * x27
    x260 = x17 * x256
    x261 = -x260
    x262 = x0 * x68
    x263 = x256 * x27
    x264 = x17 * (x12 * x49 + x67)
    x265 = -x264
    x266 = x0 * x14
    x267 = x16 * x266
    x268 = x11 * x266 - x267 + x85
    x269 = x0 * (x20 + x268)
    x270 = x22 * x266
    x271 = x267 - x270 + x91
    x272 = x0 * (x24 + x271)
    x273 = x197 + x269 - x272
    x274 = x266 * x34
    x275 = x270 - x274 + x99
    x276 = x0 * (x275 + x36)
    x277 = x199 + x272 - x276
    x278 = -x17 * x277
    x279 = x27 * x273 + x278
    x280 = x115 - x266 * x44 + x274
    x281 = -x17 * (-x0 * (x280 + x46) + x202 + x276)
    x282 = x27 * x277 + x281
    x283 = x12 * x208 + x212
    x284 = -x17 * (x12 * x210 + x217)
    x285 = 3 * x214
    x286 = (
        x0 * (3 * x209 - 3 * x211 + 3 * x213 - x285)
        - x0 * (3 * x215 - 3 * x216 - 3 * x218 + x285)
        + x27 * x283
        + x284
    )
    x287 = x12 * x221 + x225
    x288 = -x17 * (x12 * x223 + x230)
    x289 = 3 * x227
    x290 = x12 * x233 + x239
    x291 = -x17 * (x12 * x234 + x243)
    x292 = 2 * x236
    x293 = x12 * x244 + x250
    x294 = -x17 * (x12 * x245 + x254)
    x295 = 2 * x247
    x296 = x12 * x255 + x261
    x297 = -x17 * (x12 * x256 + x265)
    x298 = 2 * x258
    x299 = x12 * x273 + x278
    x300 = -x17 * (x12 * x277 + x281)
    x301 = x105 * x12 + x110
    x302 = -x17 * (x109 * x12 + x121)
    x303 = x12 * x56 + x60
    x304 = -x17 * (x12 * x59 + x70)
    x305 = 2 * x80
    x306 = 4 * x13
    x307 = x27 * x306
    x308 = x17 * x306
    x309 = x0 * (x16 * x307 - x22 * x308)
    x310 = x0 * (x11 * x307 - x16 * x308) - x17 * x271 + x268 * x27 - x309
    x311 = 2 * x93
    x312 = x0 * (x22 * x307 - x308 * x34)
    x313 = -x17 * x275 + x27 * x271 + x309 - x312
    x314 = x0 * (x305 - x311 + x313 + x52)
    x315 = x0 * (-x305 + x310 + x51 + 2 * x78) + x279 - x314
    x316 = -x0 * (x307 * x34 - x308 * x44) - x17 * x280 + x27 * x275 + x312
    x317 = -x17 * (-x0 * (-2 * x118 + x311 + x316 + x63) + x282 + x314)
    x318 = x27 * x315 + x317
    x319 = x12 * x283 + x284
    x320 = 6 * x0 * x13
    x321 = x16 * x320
    x322 = x22 * x320
    x323 = x0 * (x321 - x322 + 3 * x87 - 3 * x89)
    x324 = 3 * x272

    # 45 item(s)
    S = numpy.array(
        [
            x75,
            x126,
            x75,
            x157,
            x0 * (x112 + x165)
            - x0 * (x123 + x171)
            - x17 * (x173 - x181 + x182 - x187)
            + x27 * (x172 - x173 + x176 - x180),
            x157,
            x0 * (2 * x189 - 2 * x191)
            - x0 * (2 * x192 - 2 * x194)
            - x17 * x203
            + x200 * x27,
            -x0 * (-2 * x100 + 2 * x94)
            + x0 * (2 * x86 - 2 * x92)
            - x17 * x205
            + x204 * x27,
            x0 * (2 * x32 - 2 * x40) - x0 * (2 * x42 - 2 * x50) - x17 * x207 + x206 * x27,
            x220,
            x0 * (2 * x172 + 2 * x176 - 2 * x180 + x221 - x232)
            - x0 * (-2 * x181 + 2 * x182 - 2 * x187 + x223 + x232)
            - x17 * (x227 + x228 + x230 - x231)
            + x27 * (x222 + x225 + x226 - x227),
            x220,
            x0 * (x200 + x233)
            - x0 * (x203 + x234)
            - x17 * (x236 - x240 + x241 + x243)
            + x27 * (x235 - x236 + x237 + x239),
            x0 * (x204 + x244)
            - x0 * (x205 + x245)
            - x17 * (x247 - x251 + x252 + x254)
            + x27 * (x246 - x247 + x248 + x250),
            x0 * (x206 + x255)
            - x0 * (x207 + x256)
            - x17 * (x258 - x262 + x263 + x265)
            + x27 * (x257 - x258 + x259 + x261),
            x0 * x273 - x0 * x277 - x17 * x282 + x27 * x279,
            x0 * x105 - x0 * x109 + x111 * x27 - x122 * x17,
            x0 * x56 - x0 * x59 - x17 * x71 + x27 * x61,
            x286,
            x0 * (3 * x222 - 3 * x224 + 3 * x226 - x289)
            - x0 * (3 * x228 - 3 * x229 - 3 * x231 + x289)
            + x27 * x287
            + x288,
            x286,
            x0 * (2 * x235 + 2 * x237 - 2 * x238 - x292)
            - x0 * (-2 * x240 + 2 * x241 - 2 * x242 + x292)
            + x27 * x290
            + x291,
            x0 * (2 * x246 + 2 * x248 - 2 * x249 - x295)
            - x0 * (-2 * x251 + 2 * x252 - 2 * x253 + x295)
            + x27 * x293
            + x294,
            x0 * (2 * x257 + 2 * x259 - 2 * x260 - x298)
            - x0 * (-2 * x262 + 2 * x263 - 2 * x264 + x298)
            + x27 * x296
            + x297,
            x0 * x279 - x0 * x282 + x27 * x299 + x300,
            x0 * x111 - x0 * x122 + x27 * x301 + x302,
            x0 * x61 - x0 * x71 + x27 * x303 + x304,
            x318,
            x125,
            x74,
            x319,
            x12 * x287 + x288,
            x319,
            x12 * x290 + x291,
            x12 * x293 + x294,
            x12 * x296 + x297,
            x12 * x299 + x300,
            x12 * x301 + x302,
            x12 * x303 + x304,
            x12 * x315 + x317,
            x112 * x12 + x124,
            x12 * x62 + x73,
            x0
            * (
                x0 * (x11 * x320 - x321 + 3 * x82 - 3 * x83)
                - x17 * x313
                + 3 * x189
                - 3 * x191
                + 3 * x269
                + x27 * x310
                - x323
                - x324
            )
            - x0
            * (
                -x0 * (-x320 * x34 + x322 + 3 * x95 - 3 * x97)
                - x17 * x316
                + 3 * x192
                - 3 * x194
                + x27 * x313
                - 3 * x276
                + x323
                + x324
            )
            + x318,
            x126,
            x75,
        ]
    )
    return S


def coulomb_42(a, A, b, B, C):
    """Cartesian (gd) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = -x2 * (a * A[1] + b * B[1])
    x5 = -x2 * (a * A[2] + b * B[2])
    x6 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x4 + B[1]) ** 2 + abs(x5 + B[2]) ** 2)
    x7 = x3 + C[0]
    x8 = x4 + C[1]
    x9 = x5 + C[2]
    x10 = x7 ** 2 + x8 ** 2 + x9 ** 2
    x11 = x1 * x10
    x12 = boys(0, x11)
    x13 = numpy.pi * x2 * numpy.exp(-a * b * x10 * x2)
    x14 = 2 * x6
    x15 = x13 * x14
    x16 = boys(1, x11)
    x17 = numpy.sqrt(abs(x7) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x18 = 2 * x17
    x19 = x13 * x18
    x20 = -x16 * x19
    x21 = x12 * x15 + x20
    x22 = x13 * x16
    x23 = boys(2, x11)
    x24 = -x19 * x23
    x25 = x14 * x22 + x24
    x26 = x17 * x25
    x27 = -x26
    x28 = x21 * x6 + x27
    x29 = x0 * x28
    x30 = boys(3, x11)
    x31 = -x19 * x30
    x32 = x15 * x23 + x31
    x33 = x17 * x32
    x34 = -x33
    x35 = x25 * x6 + x34
    x36 = x0 * x35
    x37 = 3 * x36
    x38 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x4 + A[1]) ** 2 + abs(x5 + A[2]) ** 2)
    x39 = x28 * x38
    x40 = x17 * x35
    x41 = -x40
    x42 = x39 + x41
    x43 = x38 * x42
    x44 = x35 * x38
    x45 = boys(4, x11)
    x46 = -x19 * x45
    x47 = x15 * x30 + x46
    x48 = x17 * x47
    x49 = -x48
    x50 = x32 * x6 + x49
    x51 = x17 * x50
    x52 = -x51
    x53 = x44 + x52
    x54 = x17 * x53
    x55 = x0 * x50
    x56 = x38 * x53
    x57 = x38 * x50
    x58 = boys(5, x11)
    x59 = -x19 * x58
    x60 = x15 * x45 + x59
    x61 = x17 * x60
    x62 = -x61
    x63 = x47 * x6 + x62
    x64 = x17 * x63
    x65 = -x64
    x66 = x57 + x65
    x67 = x17 * x66
    x68 = 2 * x39 - 2 * x40
    x69 = 2 * x44 - 2 * x51
    x70 = x0 * x69
    x71 = -x54
    x72 = x43 + x71
    x73 = x29 - x36 + x72
    x74 = -x67
    x75 = x56 + x74
    x76 = x36 - x55 + x75
    x77 = -x17 * x76
    x78 = x38 * x73 + x77
    x79 = x0 * x68 - x70 + x78
    x80 = 2 * x57 - 2 * x64
    x81 = -x19 * boys(6, x11)
    x82 = x17 * (x15 * x58 + x81)
    x83 = -x82
    x84 = -x17 * (x6 * x60 + x83)
    x85 = x38 * x63 + x84
    x86 = -x17 * x85
    x87 = x38 * x66 + x86
    x88 = -x0 * x63 + x55 + x87
    x89 = -x17 * x88
    x90 = x38 * x76 + x89
    x91 = -x0 * x80 + x70 + x90
    x92 = -x17 * x91
    x93 = x38 * x79 + x92
    x94 = (
        x0 * (3 * x29 - x37 + 3 * x43 - 3 * x54)
        - x0 * (x37 - 3 * x55 + 3 * x56 - 3 * x67)
        + x93
    )
    x95 = 2 * x38
    x96 = x13 * x95
    x97 = x12 * x96 + x20
    x98 = x22 * x95 + x24
    x99 = x17 * x98
    x100 = -x99
    x101 = x100 + x6 * x97
    x102 = x0 * x101
    x103 = x23 * x96 + x31
    x104 = x103 * x17
    x105 = -x104
    x106 = x105 + x6 * x98
    x107 = x0 * x106
    x108 = 3 * x107
    x109 = x101 * x38
    x110 = x106 * x17
    x111 = -x110
    x112 = x109 + x111
    x113 = x112 * x38
    x114 = x106 * x38
    x115 = x30 * x96 + x46
    x116 = x115 * x17
    x117 = -x116
    x118 = x103 * x6 + x117
    x119 = x118 * x17
    x120 = -x119
    x121 = x114 + x120
    x122 = x121 * x17
    x123 = x0 * x118
    x124 = x121 * x38
    x125 = x118 * x38
    x126 = x45 * x96 + x59
    x127 = x126 * x17
    x128 = -x127
    x129 = x115 * x6 + x128
    x130 = x129 * x17
    x131 = -x130
    x132 = x125 + x131
    x133 = x132 * x17
    x134 = 2 * x109 - 2 * x110
    x135 = 2 * x114 - 2 * x119
    x136 = x0 * x135
    x137 = x102 - x107 + x113 - x122
    x138 = x107 - x123 + x124 - x133
    x139 = x0 * x134 - x136 + x137 * x38 - x138 * x17
    x140 = 2 * x125 - 2 * x130
    x141 = x58 * x96 + x81
    x142 = -x141 * x17
    x143 = -x17 * (x126 * x6 + x142)
    x144 = x129 * x38 + x143
    x145 = -x0 * x129 + x123 + x132 * x38 - x144 * x17
    x146 = -x0 * x140 + x136 + x138 * x38 - x145 * x17
    x147 = 2 * x0
    x148 = x147 * x22
    x149 = x13 * x147
    x150 = x38 * x97
    x151 = x100 + x150
    x152 = x12 * x149 - x148 + x151
    x153 = x0 * x152
    x154 = x149 * x23
    x155 = x38 * x98
    x156 = x105 + x155
    x157 = x148 - x154 + x156
    x158 = x0 * x157
    x159 = 3 * x158
    x160 = x152 * x38
    x161 = x157 * x17
    x162 = -x161
    x163 = x160 + x162
    x164 = x163 * x38
    x165 = x157 * x38
    x166 = x149 * x30
    x167 = x103 * x38
    x168 = x117 + x167
    x169 = x154 - x166 + x168
    x170 = x169 * x17
    x171 = -x170
    x172 = x165 + x171
    x173 = x17 * x172
    x174 = x0 * x169
    x175 = x172 * x38
    x176 = x169 * x38
    x177 = x149 * x45
    x178 = x115 * x38
    x179 = x128 + x178
    x180 = x166 - x177 + x179
    x181 = x17 * x180
    x182 = -x181
    x183 = x176 + x182
    x184 = x17 * x183
    x185 = x0 * (2 * x165 - 2 * x170)
    x186 = -x173
    x187 = x164 + x186
    x188 = x153 - x158 + x187
    x189 = -x184
    x190 = x175 + x189
    x191 = x158 - x174 + x190
    x192 = -x17 * x191
    x193 = x188 * x38 + x192
    x194 = x0 * (2 * x160 - 2 * x161) - x185 + x193
    x195 = x126 * x38 + x142
    x196 = -x149 * x58 + x177 + x195
    x197 = -x17 * x196
    x198 = x180 * x38 + x197
    x199 = -x17 * x198
    x200 = x183 * x38 + x199
    x201 = -x0 * x180 + x174 + x200
    x202 = -x17 * x201
    x203 = x191 * x38 + x202
    x204 = -x0 * (2 * x176 - 2 * x181) + x185 + x203
    x205 = -x17 * x204
    x206 = x194 * x38 + x205
    x207 = (
        x0 * (3 * x153 - x159 + 3 * x164 - 3 * x173)
        - x0 * (x159 - 3 * x174 + 3 * x175 - 3 * x184)
        + x206
    )
    x208 = x21 * x38
    x209 = x208 + x27
    x210 = x0 * x97
    x211 = x0 * x98
    x212 = x210 - x211
    x213 = x209 + x212
    x214 = x0 * x213
    x215 = x25 * x38
    x216 = x215 + x34
    x217 = x0 * x103
    x218 = x211 - x217
    x219 = x216 + x218
    x220 = x0 * x219
    x221 = 3 * x220
    x222 = x213 * x38
    x223 = x17 * x219
    x224 = -x223
    x225 = x222 + x224
    x226 = x225 * x38
    x227 = x219 * x38
    x228 = x32 * x38
    x229 = x228 + x49
    x230 = x0 * x115
    x231 = x217 - x230
    x232 = x229 + x231
    x233 = x17 * x232
    x234 = -x233
    x235 = x227 + x234
    x236 = x17 * x235
    x237 = x0 * x232
    x238 = x235 * x38
    x239 = x232 * x38
    x240 = x38 * x47
    x241 = x240 + x62
    x242 = x0 * x126
    x243 = x230 - x242
    x244 = x241 + x243
    x245 = x17 * x244
    x246 = -x245
    x247 = x239 + x246
    x248 = x17 * x247
    x249 = 2 * x222 - 2 * x223
    x250 = 2 * x227 - 2 * x233
    x251 = x0 * x250
    x252 = -x236
    x253 = x226 + x252
    x254 = x214 - x220 + x253
    x255 = -x248
    x256 = x238 + x255
    x257 = x220 - x237 + x256
    x258 = -x17 * x257
    x259 = x254 * x38 + x258
    x260 = x0 * x249 - x251 + x259
    x261 = 2 * x239 - 2 * x245
    x262 = x244 * x38
    x263 = x0 * x141
    x264 = x38 * x60
    x265 = x264 + x83
    x266 = x17 * (x242 - x263 + x265)
    x267 = -x266
    x268 = x262 + x267
    x269 = -x17 * x268
    x270 = x247 * x38 + x269
    x271 = -x0 * x244 + x237 + x270
    x272 = -x17 * x271
    x273 = x257 * x38 + x272
    x274 = -x0 * x261 + x251 + x273
    x275 = -x17 * x274
    x276 = x260 * x38 + x275
    x277 = (
        x0 * (3 * x214 - x221 + 3 * x226 - 3 * x236)
        - x0 * (x221 - 3 * x237 + 3 * x238 - 3 * x248)
        + x276
    )
    x278 = x0 * x42
    x279 = x0 * x53
    x280 = 2 * x279
    x281 = x28 * x6 + x41
    x282 = x281 * x38
    x283 = x35 * x6 + x52
    x284 = x17 * x283
    x285 = 2 * x278 - x280 + 2 * x282 - 2 * x284
    x286 = x0 * x66
    x287 = 2 * x286
    x288 = x283 * x38
    x289 = x50 * x6 + x65
    x290 = x17 * x289
    x291 = x280 - x287 + 2 * x288 - 2 * x290
    x292 = x0 * (x281 + x73)
    x293 = x0 * (x283 + x76)
    x294 = -x284
    x295 = x278 - x279
    x296 = x38 * (x282 + x294 + x295)
    x297 = -x290
    x298 = x279 - x286
    x299 = x288 + x297 + x298
    x300 = x17 * x299
    x301 = x0 * (x289 + x88)
    x302 = x299 * x38
    x303 = x0 * x85
    x304 = x289 * x38
    x305 = x17 * (x6 * x63 + x84)
    x306 = -x305
    x307 = x17 * (x286 - x303 + x304 + x306)
    x308 = (
        x0 * (x285 + x79)
        - x0 * (x291 + x91)
        - x17 * (x293 - x301 + x302 - x307)
        + x38 * (x292 - x293 + x296 - x300)
    )
    x309 = x0 * x112
    x310 = x0 * x121
    x311 = 2 * x310
    x312 = x101 * x6 + x111
    x313 = x312 * x38
    x314 = x106 * x6 + x120
    x315 = x17 * x314
    x316 = 2 * x309 - x311 + 2 * x313 - 2 * x315
    x317 = x0 * x132
    x318 = 2 * x317
    x319 = x314 * x38
    x320 = x118 * x6 + x131
    x321 = x17 * x320
    x322 = x311 - x318 + 2 * x319 - 2 * x321
    x323 = x0 * (x137 + x312)
    x324 = x0 * (x138 + x314)
    x325 = -x315
    x326 = x38 * (x309 - x310 + x313 + x325)
    x327 = -x321
    x328 = x310 - x317 + x319 + x327
    x329 = x17 * x328
    x330 = x0 * (x145 + x320)
    x331 = x328 * x38
    x332 = x0 * x144
    x333 = x320 * x38
    x334 = x17 * (x129 * x6 + x143)
    x335 = -x334
    x336 = x17 * (x317 - x332 + x333 + x335)
    x337 = x0 * x163
    x338 = x0 * x172
    x339 = 2 * x338
    x340 = x152 * x6 + x162
    x341 = x340 * x38
    x342 = x157 * x6 + x171
    x343 = x17 * x342
    x344 = 2 * x337 - x339 + 2 * x341 - 2 * x343
    x345 = x0 * x183
    x346 = 2 * x345
    x347 = x342 * x38
    x348 = x169 * x6 + x182
    x349 = x17 * x348
    x350 = x339 - x346 + 2 * x347 - 2 * x349
    x351 = x0 * (x188 + x340)
    x352 = x0 * (x191 + x342)
    x353 = -x343
    x354 = x337 - x338
    x355 = x38 * (x341 + x353 + x354)
    x356 = -x349
    x357 = x338 - x345
    x358 = x347 + x356 + x357
    x359 = x17 * x358
    x360 = x0 * (x201 + x348)
    x361 = x358 * x38
    x362 = x0 * x198
    x363 = x348 * x38
    x364 = x17 * (x180 * x6 + x197)
    x365 = -x364
    x366 = x17 * (x345 - x362 + x363 + x365)
    x367 = x0 * x225
    x368 = x0 * x235
    x369 = 2 * x368
    x370 = x213 * x6 + x224
    x371 = x370 * x38
    x372 = x219 * x6 + x234
    x373 = x17 * x372
    x374 = 2 * x367 - x369 + 2 * x371 - 2 * x373
    x375 = x0 * x247
    x376 = 2 * x375
    x377 = x372 * x38
    x378 = x232 * x6 + x246
    x379 = x17 * x378
    x380 = x369 - x376 + 2 * x377 - 2 * x379
    x381 = x0 * (x254 + x370)
    x382 = x0 * (x257 + x372)
    x383 = -x373
    x384 = x367 - x368
    x385 = x38 * (x371 + x383 + x384)
    x386 = -x379
    x387 = x368 - x375
    x388 = x377 + x386 + x387
    x389 = x17 * x388
    x390 = x0 * (x271 + x378)
    x391 = x38 * x388
    x392 = x0 * x268
    x393 = x378 * x38
    x394 = x17 * (x244 * x6 + x267)
    x395 = -x394
    x396 = x17 * (x375 - x392 + x393 + x395)
    x397 = 2 * x211
    x398 = 2 * x208 - 2 * x26
    x399 = 2 * x210 - x397 + x398
    x400 = x0 * x399
    x401 = 2 * x217
    x402 = 2 * x215 - 2 * x33
    x403 = x397 - x401 + x402
    x404 = x0 * x403
    x405 = x400 - x404 + x42
    x406 = x38 * x405
    x407 = 2 * x230
    x408 = 2 * x228 - 2 * x48
    x409 = x401 - x407 + x408
    x410 = x0 * x409
    x411 = x404 - x410 + x53
    x412 = x17 * x411
    x413 = x38 * x411
    x414 = 2 * x242
    x415 = 2 * x240 - 2 * x61
    x416 = x407 - x414 + x415
    x417 = x0 * x416
    x418 = x410 - x417 + x66
    x419 = x17 * x418
    x420 = x0 * x411
    x421 = -x412
    x422 = x406 + x421
    x423 = -x419
    x424 = x413 + x423
    x425 = x0 * x405 - x17 * x424 + x38 * x422 - x420
    x426 = -x17 * (-x0 * (-2 * x263 + 2 * x264 + x414 - 2 * x82) + x417 + x85)
    x427 = x38 * x418 + x426
    x428 = -x0 * x418 - x17 * x427 + x38 * x424 + x420
    x429 = x0 * x151
    x430 = x0 * x156
    x431 = x112 + x429 - x430
    x432 = x38 * x431
    x433 = x0 * x168
    x434 = x121 + x430 - x433
    x435 = x17 * x434
    x436 = x38 * x434
    x437 = x0 * x179
    x438 = x132 + x433 - x437
    x439 = x17 * x438
    x440 = x0 * x434
    x441 = -x435
    x442 = x432 + x441
    x443 = -x439
    x444 = x436 + x443
    x445 = x0 * x431 - x17 * x444 + x38 * x442 - x440
    x446 = -x17 * (-x0 * x195 + x144 + x437)
    x447 = x38 * x438 + x446
    x448 = -x0 * x438 - x17 * x447 + x38 * x444 + x440
    x449 = x0 * x209
    x450 = x0 * x216
    x451 = x42 + x449 - x450
    x452 = x38 * x451
    x453 = x0 * x229
    x454 = x450 - x453 + x53
    x455 = x17 * x454
    x456 = x38 * x454
    x457 = x0 * x241
    x458 = x453 - x457 + x66
    x459 = x17 * x458
    x460 = x0 * x454
    x461 = -x455
    x462 = x452 + x461
    x463 = -x459
    x464 = x456 + x463
    x465 = x0 * x451 - x17 * x464 + x38 * x462 - x460
    x466 = -x17 * (-x0 * x265 + x457 + x85)
    x467 = x38 * x458 + x466
    x468 = -x0 * x458 - x17 * x467 + x38 * x464 + x460
    x469 = -x17 * x190 + x187 * x38 + x354
    x470 = -x17 * x200 + x190 * x38 + x357
    x471 = -x17 * x256 + x253 * x38 + x384
    x472 = -x17 * x270 + x256 * x38 + x387
    x473 = -x17 * x75 + x295 + x38 * x72
    x474 = -x17 * x87 + x298 + x38 * x75
    x475 = x281 * x6 + x294
    x476 = x38 * x475
    x477 = x283 * x6 + x297
    x478 = x17 * x477
    x479 = -x478
    x480 = x0 * x285
    x481 = x0 * x291
    x482 = x38 * x477
    x483 = x17 * (x289 * x6 + x306)
    x484 = -x483
    x485 = x0 * (x287 - 2 * x303 + 2 * x304 - 2 * x305)
    x486 = 2 * x293
    x487 = (
        x0 * (2 * x292 + 2 * x296 - 2 * x300 + x475 - x486)
        - x0 * (-2 * x301 + 2 * x302 - 2 * x307 + x477 + x486)
        - x17 * (x481 + x482 + x484 - x485)
        + x38 * (x476 + x479 + x480 - x481)
    )
    x488 = x312 * x6 + x325
    x489 = x38 * x488
    x490 = x314 * x6 + x327
    x491 = x17 * x490
    x492 = -x491
    x493 = x0 * x316
    x494 = x0 * x322
    x495 = x38 * x490
    x496 = x17 * (x320 * x6 + x335)
    x497 = -x496
    x498 = x0 * (x318 - 2 * x332 + 2 * x333 - 2 * x334)
    x499 = 2 * x324
    x500 = x340 * x6 + x353
    x501 = x38 * x500
    x502 = x342 * x6 + x356
    x503 = x17 * x502
    x504 = -x503
    x505 = x0 * x344
    x506 = x0 * x350
    x507 = x38 * x502
    x508 = x17 * (x348 * x6 + x365)
    x509 = -x508
    x510 = x0 * (x346 - 2 * x362 + 2 * x363 - 2 * x364)
    x511 = 2 * x352
    x512 = x370 * x6 + x383
    x513 = x38 * x512
    x514 = x372 * x6 + x386
    x515 = x17 * x514
    x516 = -x515
    x517 = x0 * x374
    x518 = x0 * x380
    x519 = x38 * x514
    x520 = x17 * (x378 * x6 + x395)
    x521 = -x520
    x522 = x0 * (x376 - 2 * x392 + 2 * x393 - 2 * x394)
    x523 = 2 * x382
    x524 = x405 * x6 + x421
    x525 = x411 * x6 + x423
    x526 = x0 * x422
    x527 = x0 * x424
    x528 = x38 * x524
    x529 = x17 * x525
    x530 = -x529
    x531 = x0 * x427
    x532 = x38 * x525
    x533 = x17 * (x418 * x6 + x426)
    x534 = -x533
    x535 = x431 * x6 + x441
    x536 = x434 * x6 + x443
    x537 = x0 * x442
    x538 = x0 * x444
    x539 = x38 * x535
    x540 = x17 * x536
    x541 = -x540
    x542 = x0 * x447
    x543 = x38 * x536
    x544 = x17 * (x438 * x6 + x446)
    x545 = -x544
    x546 = x451 * x6 + x461
    x547 = x454 * x6 + x463
    x548 = x0 * x462
    x549 = x0 * x464
    x550 = x38 * x546
    x551 = x17 * x547
    x552 = -x551
    x553 = x0 * x467
    x554 = x38 * x547
    x555 = x17 * (x458 * x6 + x466)
    x556 = -x555
    x557 = x163 * x6 + x186
    x558 = x172 * x6 + x189
    x559 = x0 * x187
    x560 = x0 * x190
    x561 = x38 * x557
    x562 = x17 * x558
    x563 = -x562
    x564 = x0 * x200
    x565 = x38 * x558
    x566 = x17 * (x183 * x6 + x199)
    x567 = -x566
    x568 = x225 * x6 + x252
    x569 = x235 * x6 + x255
    x570 = x0 * x253
    x571 = x0 * x256
    x572 = x38 * x568
    x573 = x17 * x569
    x574 = -x573
    x575 = x0 * x270
    x576 = x38 * x569
    x577 = x17 * (x247 * x6 + x269)
    x578 = -x577
    x579 = x42 * x6 + x71
    x580 = x53 * x6 + x74
    x581 = x0 * x72
    x582 = x0 * x75
    x583 = x38 * x579
    x584 = x17 * x580
    x585 = -x584
    x586 = x0 * x87
    x587 = x38 * x580
    x588 = x17 * (x6 * x66 + x86)
    x589 = -x588
    x590 = x0 * (x157 + x25)
    x591 = 2 * x590
    x592 = x0 * (x152 + x21)
    x593 = x0 * (x249 + x28 - x591 + 2 * x592)
    x594 = x0 * (x169 + x32)
    x595 = 2 * x594
    x596 = x0 * (x250 + x35 + x591 - x595)
    x597 = x422 + x593 - x596
    x598 = x0 * (x180 + x47)
    x599 = 2 * x598
    x600 = x0 * (x261 + x50 + x595 - x599)
    x601 = x424 + x596 - x600
    x602 = -x17 * x601
    x603 = x38 * x597 + x602
    x604 = x0 * (x196 + x60)
    x605 = -x17 * (-x0 * (2 * x262 - 2 * x266 + x599 - 2 * x604 + x63) + x427 + x600)
    x606 = x38 * x601 + x605
    x607 = x151 * x38
    x608 = x156 * x17
    x609 = x212 + x607 - x608
    x610 = x0 * (x101 + x609)
    x611 = x156 * x38
    x612 = x168 * x17
    x613 = x218 + x611 - x612
    x614 = x0 * (x106 + x613)
    x615 = x442 + x610 - x614
    x616 = x168 * x38
    x617 = x17 * x179
    x618 = x231 + x616 - x617
    x619 = x0 * (x118 + x618)
    x620 = x444 + x614 - x619
    x621 = -x17 * x620
    x622 = x38 * x615 + x621
    x623 = -x17 * x195 + x179 * x38 + x243
    x624 = -x17 * (-x0 * (x129 + x623) + x447 + x619)
    x625 = x38 * x620 + x624
    x626 = x0 * x21
    x627 = x0 * x25
    x628 = x209 * x38
    x629 = x17 * x216
    x630 = x626 - x627 + x628 - x629
    x631 = x0 * (x28 + x630)
    x632 = x0 * x32
    x633 = x216 * x38
    x634 = x17 * x229
    x635 = x627 - x632 + x633 - x634
    x636 = x0 * (x35 + x635)
    x637 = x462 + x631 - x636
    x638 = x0 * x47
    x639 = x229 * x38
    x640 = x17 * x241
    x641 = x632 - x638 + x639 - x640
    x642 = x0 * (x50 + x641)
    x643 = x464 + x636 - x642
    x644 = -x17 * x643
    x645 = x38 * x637 + x644
    x646 = -x0 * x60 - x17 * x265 + x241 * x38 + x638
    x647 = -x17 * (-x0 * (x63 + x646) + x467 + x642)
    x648 = x38 * x643 + x647
    x649 = x475 * x6 + x479
    x650 = -x17 * (x477 * x6 + x484)
    x651 = 3 * x481
    x652 = (
        x0 * (3 * x476 - 3 * x478 + 3 * x480 - x651)
        - x0 * (3 * x482 - 3 * x483 - 3 * x485 + x651)
        + x38 * x649
        + x650
    )
    x653 = x488 * x6 + x492
    x654 = -x17 * (x490 * x6 + x497)
    x655 = 3 * x494
    x656 = x500 * x6 + x504
    x657 = -x17 * (x502 * x6 + x509)
    x658 = 3 * x506
    x659 = x512 * x6 + x516
    x660 = -x17 * (x514 * x6 + x521)
    x661 = 3 * x518
    x662 = x524 * x6 + x530
    x663 = -x17 * (x525 * x6 + x534)
    x664 = 2 * x527
    x665 = x535 * x6 + x541
    x666 = -x17 * (x536 * x6 + x545)
    x667 = 2 * x538
    x668 = x546 * x6 + x552
    x669 = -x17 * (x547 * x6 + x556)
    x670 = 2 * x549
    x671 = x557 * x6 + x563
    x672 = -x17 * (x558 * x6 + x567)
    x673 = 2 * x560
    x674 = x568 * x6 + x574
    x675 = -x17 * (x569 * x6 + x578)
    x676 = 2 * x571
    x677 = x579 * x6 + x585
    x678 = -x17 * (x580 * x6 + x589)
    x679 = 2 * x582
    x680 = x597 * x6 + x602
    x681 = -x17 * (x6 * x601 + x605)
    x682 = x6 * x615 + x621
    x683 = -x17 * (x6 * x620 + x624)
    x684 = x6 * x637 + x644
    x685 = -x17 * (x6 * x643 + x647)
    x686 = x188 * x6 + x192
    x687 = -x17 * (x191 * x6 + x202)
    x688 = x254 * x6 + x258
    x689 = -x17 * (x257 * x6 + x272)
    x690 = x6 * x73 + x77
    x691 = -x17 * (x6 * x76 + x89)
    x692 = 2 * x404
    x693 = 4 * x38
    x694 = 4 * x17
    x695 = x13 * x23
    x696 = x0 * (x22 * x693 - x694 * x695)
    x697 = x13 * x694
    x698 = x0 * (-x30 * x697 + x693 * x695)
    x699 = x172 + x696 - x698
    x700 = x0 * (x403 + x699)
    x701 = 2 * x700
    x702 = x235 + x590 - x594
    x703 = x17 * x702
    x704 = x13 * x693
    x705 = x0 * (x12 * x704 - x22 * x694) + x163 - x696
    x706 = x0 * (x399 + x705)
    x707 = x38 * (x225 - x590 + x592)
    x708 = 2 * x410
    x709 = x0 * (x30 * x704 - x45 * x697)
    x710 = x183 + x698 - x709
    x711 = x0 * (x409 + x710)
    x712 = 2 * x711
    x713 = x247 + x594 - x598
    x714 = x17 * x713
    x715 = x38 * x702
    x716 = x0 * (x69 + x692 + x701 - x708 - x712 - 2 * x714 + 2 * x715)
    x717 = (
        x0 * (2 * x400 + x68 - x692 - x701 - 2 * x703 + 2 * x706 + 2 * x707) + x603 - x716
    )
    x718 = -x0 * (x45 * x704 - x58 * x697) + x198 + x709
    x719 = x0 * (x416 + x718)
    x720 = x17 * (x268 + x598 - x604)
    x721 = x38 * x713
    x722 = -x17 * (
        -x0 * (-2 * x417 + x708 + x712 - 2 * x719 - 2 * x720 + 2 * x721 + x80)
        + x606
        + x716
    )
    x723 = x38 * x717 + x722
    x724 = 2 * x430
    x725 = x0 * (-2 * x104 + 2 * x155)
    x726 = x0 * (2 * x150 - 2 * x99) - x17 * x613 + x38 * x609 - x725
    x727 = 2 * x433
    x728 = x0 * (-2 * x116 + 2 * x167)
    x729 = -x17 * x618 + x38 * x613 + x725 - x728
    x730 = x0 * (x135 + x724 - x727 + x729)
    x731 = x0 * (x134 + 2 * x429 - x724 + x726) + x622 - x730
    x732 = -x0 * (-2 * x127 + 2 * x178) - x17 * x623 + x38 * x618 + x728
    x733 = -x17 * (-x0 * (x140 - 2 * x437 + x727 + x732) + x625 + x730)
    x734 = x38 * x731 + x733
    x735 = 2 * x450
    x736 = x0 * x402
    x737 = x0 * x398 - x17 * x635 + x38 * x630 - x736
    x738 = 2 * x453
    x739 = x0 * x408
    x740 = -x17 * x641 + x38 * x635 + x736 - x739
    x741 = x0 * (x69 + x735 - x738 + x740)
    x742 = x0 * (2 * x449 + x68 - x735 + x737) + x645 - x741
    x743 = -x0 * x415 - x17 * x646 + x38 * x641 + x739
    x744 = -x17 * (-x0 * (-2 * x457 + x738 + x743 + x80) + x648 + x741)
    x745 = x38 * x742 + x744
    x746 = x6 * x649 + x650
    x747 = 3 * x596
    x748 = 6 * x0
    x749 = x13 * x748
    x750 = x22 * x748
    x751 = x695 * x748
    x752 = x0 * (-3 * x104 + 3 * x155 + x750 - x751)
    x753 = 3 * x590
    x754 = x30 * x749
    x755 = x0 * (-3 * x116 + 3 * x167 + x751 - x754)
    x756 = 3 * x594
    x757 = x147 * (
        -x17 * x710 + 3 * x227 - 3 * x233 + x38 * x699 + x752 + x753 - x755 - x756
    )
    x758 = x700 - x711 - x714 + x715
    x759 = 3 * x211
    x760 = 3 * x217
    x761 = x0 * (3 * x611 - 3 * x612 + x759 - x760)
    x762 = 3 * x614
    x763 = 3 * x627
    x764 = 3 * x632
    x765 = x0 * (3 * x633 - 3 * x634 + x763 - x764)
    x766 = 3 * x636

    # 90 item(s)
    S = numpy.array(
        [
            x94,
            x0 * (3 * x102 - x108 + 3 * x113 - 3 * x122)
            - x0 * (x108 - 3 * x123 + 3 * x124 - 3 * x133)
            + x139 * x38
            - x146 * x17,
            x94,
            x207,
            x277,
            x94,
            x308,
            x0 * (x139 + x316)
            - x0 * (x146 + x322)
            - x17 * (x324 - x330 + x331 - x336)
            + x38 * (x323 - x324 + x326 - x329),
            x308,
            x0 * (x194 + x344)
            - x0 * (x204 + x350)
            - x17 * (x352 - x360 + x361 - x366)
            + x38 * (x351 - x352 + x355 - x359),
            x0 * (x260 + x374)
            - x0 * (x274 + x380)
            - x17 * (x382 - x390 + x391 - x396)
            + x38 * (x381 - x382 + x385 - x389),
            x308,
            x0 * (2 * x406 - 2 * x412)
            - x0 * (2 * x413 - 2 * x419)
            - x17 * x428
            + x38 * x425,
            x0 * (2 * x432 - 2 * x435)
            - x0 * (2 * x436 - 2 * x439)
            - x17 * x448
            + x38 * x445,
            x0 * (2 * x452 - 2 * x455)
            - x0 * (2 * x456 - 2 * x459)
            - x17 * x468
            + x38 * x465,
            x0 * (2 * x164 - 2 * x173)
            - x0 * (2 * x175 - 2 * x184)
            - x17 * x470
            + x38 * x469,
            x0 * (2 * x226 - 2 * x236)
            - x0 * (2 * x238 - 2 * x248)
            - x17 * x472
            + x38 * x471,
            x0 * (2 * x43 - 2 * x54) - x0 * (2 * x56 - 2 * x67) - x17 * x474 + x38 * x473,
            x487,
            x0 * (2 * x323 + 2 * x326 - 2 * x329 + x488 - x499)
            - x0 * (-2 * x330 + 2 * x331 - 2 * x336 + x490 + x499)
            - x17 * (x494 + x495 + x497 - x498)
            + x38 * (x489 + x492 + x493 - x494),
            x487,
            x0 * (2 * x351 + 2 * x355 - 2 * x359 + x500 - x511)
            - x0 * (-2 * x360 + 2 * x361 - 2 * x366 + x502 + x511)
            - x17 * (x506 + x507 + x509 - x510)
            + x38 * (x501 + x504 + x505 - x506),
            x0 * (2 * x381 + 2 * x385 - 2 * x389 + x512 - x523)
            - x0 * (-2 * x390 + 2 * x391 - 2 * x396 + x514 + x523)
            - x17 * (x518 + x519 + x521 - x522)
            + x38 * (x513 + x516 + x517 - x518),
            x487,
            x0 * (x425 + x524)
            - x0 * (x428 + x525)
            - x17 * (x527 - x531 + x532 + x534)
            + x38 * (x526 - x527 + x528 + x530),
            x0 * (x445 + x535)
            - x0 * (x448 + x536)
            - x17 * (x538 - x542 + x543 + x545)
            + x38 * (x537 - x538 + x539 + x541),
            x0 * (x465 + x546)
            - x0 * (x468 + x547)
            - x17 * (x549 - x553 + x554 + x556)
            + x38 * (x548 - x549 + x550 + x552),
            x0 * (x469 + x557)
            - x0 * (x470 + x558)
            - x17 * (x560 - x564 + x565 + x567)
            + x38 * (x559 - x560 + x561 + x563),
            x0 * (x471 + x568)
            - x0 * (x472 + x569)
            - x17 * (x571 - x575 + x576 + x578)
            + x38 * (x570 - x571 + x572 + x574),
            x0 * (x473 + x579)
            - x0 * (x474 + x580)
            - x17 * (x582 - x586 + x587 + x589)
            + x38 * (x581 - x582 + x583 + x585),
            x0 * x597 - x0 * x601 - x17 * x606 + x38 * x603,
            x0 * x615 - x0 * x620 - x17 * x625 + x38 * x622,
            x0 * x637 - x0 * x643 - x17 * x648 + x38 * x645,
            x0 * x188 - x0 * x191 - x17 * x203 + x193 * x38,
            x0 * x254 - x0 * x257 - x17 * x273 + x259 * x38,
            x0 * x73 - x0 * x76 - x17 * x90 + x38 * x78,
            x652,
            x0 * (3 * x489 - 3 * x491 + 3 * x493 - x655)
            - x0 * (3 * x495 - 3 * x496 - 3 * x498 + x655)
            + x38 * x653
            + x654,
            x652,
            x0 * (3 * x501 - 3 * x503 + 3 * x505 - x658)
            - x0 * (3 * x507 - 3 * x508 - 3 * x510 + x658)
            + x38 * x656
            + x657,
            x0 * (3 * x513 - 3 * x515 + 3 * x517 - x661)
            - x0 * (3 * x519 - 3 * x520 - 3 * x522 + x661)
            + x38 * x659
            + x660,
            x652,
            x0 * (2 * x526 + 2 * x528 - 2 * x529 - x664)
            - x0 * (-2 * x531 + 2 * x532 - 2 * x533 + x664)
            + x38 * x662
            + x663,
            x0 * (2 * x537 + 2 * x539 - 2 * x540 - x667)
            - x0 * (-2 * x542 + 2 * x543 - 2 * x544 + x667)
            + x38 * x665
            + x666,
            x0 * (2 * x548 + 2 * x550 - 2 * x551 - x670)
            - x0 * (-2 * x553 + 2 * x554 - 2 * x555 + x670)
            + x38 * x668
            + x669,
            x0 * (2 * x559 + 2 * x561 - 2 * x562 - x673)
            - x0 * (-2 * x564 + 2 * x565 - 2 * x566 + x673)
            + x38 * x671
            + x672,
            x0 * (2 * x570 + 2 * x572 - 2 * x573 - x676)
            - x0 * (-2 * x575 + 2 * x576 - 2 * x577 + x676)
            + x38 * x674
            + x675,
            x0 * (2 * x581 + 2 * x583 - 2 * x584 - x679)
            - x0 * (-2 * x586 + 2 * x587 - 2 * x588 + x679)
            + x38 * x677
            + x678,
            x0 * x603 - x0 * x606 + x38 * x680 + x681,
            x0 * x622 - x0 * x625 + x38 * x682 + x683,
            x0 * x645 - x0 * x648 + x38 * x684 + x685,
            x0 * x193 - x0 * x203 + x38 * x686 + x687,
            x0 * x259 - x0 * x273 + x38 * x688 + x689,
            x0 * x78 - x0 * x90 + x38 * x690 + x691,
            x723,
            x734,
            x745,
            x206,
            x276,
            x93,
            x746,
            x6 * x653 + x654,
            x746,
            x6 * x656 + x657,
            x6 * x659 + x660,
            x746,
            x6 * x662 + x663,
            x6 * x665 + x666,
            x6 * x668 + x669,
            x6 * x671 + x672,
            x6 * x674 + x675,
            x6 * x677 + x678,
            x6 * x680 + x681,
            x6 * x682 + x683,
            x6 * x684 + x685,
            x6 * x686 + x687,
            x6 * x688 + x689,
            x6 * x690 + x691,
            x6 * x717 + x722,
            x6 * x731 + x733,
            x6 * x742 + x744,
            x194 * x6 + x205,
            x260 * x6 + x275,
            x6 * x79 + x92,
            -x0
            * (
                -x147
                * (
                    -x0 * (-3 * x127 + 3 * x178 - x45 * x749 + x754)
                    - x17 * x718
                    + 3 * x239
                    - 3 * x245
                    + x38 * x710
                    - 3 * x598
                    + x755
                    + x756
                )
                - x18 * (x711 - x719 - x720 + x721)
                + 3 * x413
                - 3 * x419
                - 3 * x600
                + x747
                + x757
                + x758 * x95
            )
            + x0
            * (
                x147
                * (
                    x0 * (x12 * x749 + 3 * x150 - x750 - 3 * x99)
                    - x17 * x699
                    + 3 * x222
                    - 3 * x223
                    + x38 * x705
                    + 3 * x592
                    - x752
                    - x753
                )
                - x18 * x758
                + 3 * x406
                - 3 * x412
                + 3 * x593
                - x747
                - x757
                + x95 * (-x700 - x703 + x706 + x707)
            )
            + x723,
            x0
            * (
                x0 * (3 * x210 + 3 * x607 - 3 * x608 - x759)
                - x17 * x729
                + x38 * x726
                + 3 * x432
                - 3 * x435
                + 3 * x610
                - x761
                - x762
            )
            - x0
            * (
                -x0 * (-3 * x230 + 3 * x616 - 3 * x617 + x760)
                - x17 * x732
                + x38 * x729
                + 3 * x436
                - 3 * x439
                - 3 * x619
                + x761
                + x762
            )
            + x734,
            x0
            * (
                x0 * (3 * x626 + 3 * x628 - 3 * x629 - x763)
                - x17 * x740
                + x38 * x737
                + 3 * x452
                - 3 * x455
                + 3 * x631
                - x765
                - x766
            )
            - x0
            * (
                -x0 * (-3 * x638 + 3 * x639 - 3 * x640 + x764)
                - x17 * x743
                + x38 * x740
                + 3 * x456
                - 3 * x459
                - 3 * x642
                + x765
                + x766
            )
            + x745,
            x207,
            x277,
            x94,
        ]
    )
    return S


def coulomb_43(a, A, b, B, C):
    """Cartesian (gf) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = -x2 * (a * A[1] + b * B[1])
    x5 = -x2 * (a * A[2] + b * B[2])
    x6 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x4 + B[1]) ** 2 + abs(x5 + B[2]) ** 2)
    x7 = x3 + C[0]
    x8 = x4 + C[1]
    x9 = x5 + C[2]
    x10 = x7 ** 2 + x8 ** 2 + x9 ** 2
    x11 = x1 * x10
    x12 = boys(0, x11)
    x13 = numpy.pi * x2 * numpy.exp(-a * b * x10 * x2)
    x14 = 2 * x6
    x15 = x13 * x14
    x16 = boys(1, x11)
    x17 = numpy.sqrt(abs(x7) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x18 = 2 * x17
    x19 = x13 * x18
    x20 = -x16 * x19
    x21 = x12 * x15 + x20
    x22 = x13 * x16
    x23 = boys(2, x11)
    x24 = -x19 * x23
    x25 = x14 * x22 + x24
    x26 = x17 * x25
    x27 = -x26
    x28 = x21 * x6 + x27
    x29 = boys(3, x11)
    x30 = -x19 * x29
    x31 = x15 * x23 + x30
    x32 = x17 * x31
    x33 = -x32
    x34 = x25 * x6 + x33
    x35 = x17 * x34
    x36 = -x35
    x37 = x28 * x6 + x36
    x38 = x0 * x37
    x39 = boys(4, x11)
    x40 = -x19 * x39
    x41 = x15 * x29 + x40
    x42 = x17 * x41
    x43 = -x42
    x44 = x31 * x6 + x43
    x45 = x17 * x44
    x46 = -x45
    x47 = x34 * x6 + x46
    x48 = x0 * x47
    x49 = 3 * x48
    x50 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x4 + A[1]) ** 2 + abs(x5 + A[2]) ** 2)
    x51 = x37 * x50
    x52 = x17 * x47
    x53 = -x52
    x54 = x51 + x53
    x55 = x50 * x54
    x56 = x47 * x50
    x57 = boys(5, x11)
    x58 = -x19 * x57
    x59 = x15 * x39 + x58
    x60 = x17 * x59
    x61 = -x60
    x62 = x41 * x6 + x61
    x63 = x17 * x62
    x64 = -x63
    x65 = x44 * x6 + x64
    x66 = x17 * x65
    x67 = -x66
    x68 = x56 + x67
    x69 = x17 * x68
    x70 = x0 * x65
    x71 = x50 * x68
    x72 = x50 * x65
    x73 = boys(6, x11)
    x74 = -x19 * x73
    x75 = x15 * x57 + x74
    x76 = x17 * x75
    x77 = -x76
    x78 = x59 * x6 + x77
    x79 = x17 * x78
    x80 = -x79
    x81 = x6 * x62 + x80
    x82 = x17 * x81
    x83 = -x82
    x84 = x72 + x83
    x85 = x17 * x84
    x86 = 2 * x51 - 2 * x52
    x87 = 2 * x56 - 2 * x66
    x88 = x0 * x87
    x89 = -x69
    x90 = x55 + x89
    x91 = x38 - x48 + x90
    x92 = -x85
    x93 = x71 + x92
    x94 = x48 - x70 + x93
    x95 = -x17 * x94
    x96 = x50 * x91 + x95
    x97 = x0 * x86 - x88 + x96
    x98 = 2 * x72 - 2 * x82
    x99 = -x19 * boys(7, x11)
    x100 = x17 * (x15 * x73 + x99)
    x101 = -x100
    x102 = x17 * (x101 + x6 * x75)
    x103 = -x102
    x104 = -x17 * (x103 + x6 * x78)
    x105 = x104 + x50 * x81
    x106 = -x105 * x17
    x107 = x106 + x50 * x84
    x108 = -x0 * x81 + x107 + x70
    x109 = -x108 * x17
    x110 = x109 + x50 * x94
    x111 = -x0 * x98 + x110 + x88
    x112 = -x111 * x17
    x113 = x112 + x50 * x97
    x114 = (
        x0 * (3 * x38 - x49 + 3 * x55 - 3 * x69)
        - x0 * (x49 - 3 * x70 + 3 * x71 - 3 * x85)
        + x113
    )
    x115 = 2 * x50
    x116 = x115 * x13
    x117 = x116 * x12 + x20
    x118 = x115 * x22 + x24
    x119 = x118 * x17
    x120 = -x119
    x121 = x117 * x6 + x120
    x122 = x116 * x23 + x30
    x123 = x122 * x17
    x124 = -x123
    x125 = x118 * x6 + x124
    x126 = x125 * x17
    x127 = -x126
    x128 = x121 * x6 + x127
    x129 = x0 * x128
    x130 = x116 * x29 + x40
    x131 = x130 * x17
    x132 = -x131
    x133 = x122 * x6 + x132
    x134 = x133 * x17
    x135 = -x134
    x136 = x125 * x6 + x135
    x137 = x0 * x136
    x138 = 3 * x137
    x139 = x128 * x50
    x140 = x136 * x17
    x141 = -x140
    x142 = x139 + x141
    x143 = x142 * x50
    x144 = x136 * x50
    x145 = x116 * x39 + x58
    x146 = x145 * x17
    x147 = -x146
    x148 = x130 * x6 + x147
    x149 = x148 * x17
    x150 = -x149
    x151 = x133 * x6 + x150
    x152 = x151 * x17
    x153 = -x152
    x154 = x144 + x153
    x155 = x154 * x17
    x156 = x0 * x151
    x157 = x154 * x50
    x158 = x151 * x50
    x159 = x116 * x57 + x74
    x160 = x159 * x17
    x161 = -x160
    x162 = x145 * x6 + x161
    x163 = x162 * x17
    x164 = -x163
    x165 = x148 * x6 + x164
    x166 = x165 * x17
    x167 = -x166
    x168 = x158 + x167
    x169 = x168 * x17
    x170 = 2 * x139 - 2 * x140
    x171 = 2 * x144 - 2 * x152
    x172 = x0 * x171
    x173 = x129 - x137 + x143 - x155
    x174 = x137 - x156 + x157 - x169
    x175 = x0 * x170 - x17 * x174 - x172 + x173 * x50
    x176 = 2 * x158 - 2 * x166
    x177 = x116 * x73 + x99
    x178 = -x17 * x177
    x179 = x17 * (x159 * x6 + x178)
    x180 = -x179
    x181 = -x17 * (x162 * x6 + x180)
    x182 = x165 * x50 + x181
    x183 = -x0 * x165 + x156 + x168 * x50 - x17 * x182
    x184 = -x0 * x176 - x17 * x183 + x172 + x174 * x50
    x185 = 2 * x0
    x186 = x185 * x22
    x187 = x13 * x185
    x188 = x117 * x50
    x189 = x120 + x188
    x190 = x12 * x187 - x186 + x189
    x191 = x187 * x23
    x192 = x118 * x50
    x193 = x124 + x192
    x194 = x186 - x191 + x193
    x195 = x17 * x194
    x196 = -x195
    x197 = x190 * x6 + x196
    x198 = x0 * x197
    x199 = x187 * x29
    x200 = x122 * x50
    x201 = x132 + x200
    x202 = x191 - x199 + x201
    x203 = x17 * x202
    x204 = -x203
    x205 = x194 * x6 + x204
    x206 = x0 * x205
    x207 = 3 * x206
    x208 = x197 * x50
    x209 = x17 * x205
    x210 = -x209
    x211 = x208 + x210
    x212 = x211 * x50
    x213 = x205 * x50
    x214 = x187 * x39
    x215 = x130 * x50
    x216 = x147 + x215
    x217 = x199 - x214 + x216
    x218 = x17 * x217
    x219 = -x218
    x220 = x202 * x6 + x219
    x221 = x17 * x220
    x222 = -x221
    x223 = x213 + x222
    x224 = x17 * x223
    x225 = x0 * x220
    x226 = x223 * x50
    x227 = x220 * x50
    x228 = x187 * x57
    x229 = x145 * x50
    x230 = x161 + x229
    x231 = x214 - x228 + x230
    x232 = x17 * x231
    x233 = -x232
    x234 = x217 * x6 + x233
    x235 = x17 * x234
    x236 = -x235
    x237 = x227 + x236
    x238 = x17 * x237
    x239 = 2 * x208 - 2 * x209
    x240 = 2 * x213 - 2 * x221
    x241 = x0 * x240
    x242 = x198 - x206 + x212 - x224
    x243 = x206 - x225 + x226 - x238
    x244 = x0 * x239 - x17 * x243 - x241 + x242 * x50
    x245 = 2 * x227 - 2 * x235
    x246 = x159 * x50 + x178
    x247 = -x187 * x73 + x228 + x246
    x248 = -x17 * x247
    x249 = -x17 * (x231 * x6 + x248)
    x250 = x234 * x50 + x249
    x251 = -x0 * x234 - x17 * x250 + x225 + x237 * x50
    x252 = -x0 * x245 - x17 * x251 + x241 + x243 * x50
    x253 = x21 * x50
    x254 = x253 + x27
    x255 = x0 * x117
    x256 = x0 * x118
    x257 = x255 - x256
    x258 = x254 + x257
    x259 = x25 * x50
    x260 = x259 + x33
    x261 = x0 * x122
    x262 = x256 - x261
    x263 = x260 + x262
    x264 = x17 * x263
    x265 = -x264
    x266 = x258 * x6 + x265
    x267 = x0 * x266
    x268 = x31 * x50
    x269 = x268 + x43
    x270 = x0 * x130
    x271 = x261 - x270
    x272 = x269 + x271
    x273 = x17 * x272
    x274 = -x273
    x275 = x263 * x6 + x274
    x276 = x0 * x275
    x277 = 3 * x276
    x278 = x266 * x50
    x279 = x17 * x275
    x280 = -x279
    x281 = x278 + x280
    x282 = x281 * x50
    x283 = x275 * x50
    x284 = x41 * x50
    x285 = x284 + x61
    x286 = x0 * x145
    x287 = x270 - x286
    x288 = x285 + x287
    x289 = x17 * x288
    x290 = -x289
    x291 = x272 * x6 + x290
    x292 = x17 * x291
    x293 = -x292
    x294 = x283 + x293
    x295 = x17 * x294
    x296 = x0 * x291
    x297 = x294 * x50
    x298 = x291 * x50
    x299 = x50 * x59
    x300 = x299 + x77
    x301 = x0 * x159
    x302 = x286 - x301
    x303 = x300 + x302
    x304 = x17 * x303
    x305 = -x304
    x306 = x288 * x6 + x305
    x307 = x17 * x306
    x308 = -x307
    x309 = x298 + x308
    x310 = x17 * x309
    x311 = 2 * x278 - 2 * x279
    x312 = 2 * x283 - 2 * x292
    x313 = x0 * x312
    x314 = x267 - x276 + x282 - x295
    x315 = x276 - x296 + x297 - x310
    x316 = x0 * x311 - x17 * x315 - x313 + x314 * x50
    x317 = 2 * x298 - 2 * x307
    x318 = x0 * x177
    x319 = x50 * x75
    x320 = x101 + x319
    x321 = x17 * (x301 - x318 + x320)
    x322 = -x321
    x323 = -x17 * (x303 * x6 + x322)
    x324 = x306 * x50 + x323
    x325 = -x0 * x306 - x17 * x324 + x296 + x309 * x50
    x326 = -x0 * x317 - x17 * x325 + x313 + x315 * x50
    x327 = 4 * x13
    x328 = x327 * x50
    x329 = 4 * x22
    x330 = x17 * x327
    x331 = x0 * (-x23 * x330 + x329 * x50)
    x332 = x190 * x50
    x333 = x196 + x332
    x334 = x0 * (x12 * x328 - x17 * x329) - x331 + x333
    x335 = x0 * x334
    x336 = x0 * (x23 * x328 - x29 * x330)
    x337 = x194 * x50
    x338 = x204 + x337
    x339 = x331 - x336 + x338
    x340 = x0 * x339
    x341 = 3 * x340
    x342 = x334 * x50
    x343 = x17 * x339
    x344 = -x343
    x345 = x342 + x344
    x346 = x345 * x50
    x347 = x339 * x50
    x348 = x0 * (x29 * x328 - x330 * x39)
    x349 = x202 * x50
    x350 = x219 + x349
    x351 = x336 - x348 + x350
    x352 = x17 * x351
    x353 = -x352
    x354 = x347 + x353
    x355 = x17 * x354
    x356 = x0 * x351
    x357 = x354 * x50
    x358 = x351 * x50
    x359 = x0 * (x328 * x39 - x330 * x57)
    x360 = x217 * x50
    x361 = x233 + x360
    x362 = x348 - x359 + x361
    x363 = x17 * x362
    x364 = -x363
    x365 = x358 + x364
    x366 = x17 * x365
    x367 = x0 * (2 * x347 - 2 * x352)
    x368 = -x355
    x369 = x346 + x368
    x370 = x335 - x340 + x369
    x371 = -x366
    x372 = x357 + x371
    x373 = x340 - x356 + x372
    x374 = -x17 * x373
    x375 = x370 * x50 + x374
    x376 = x0 * (2 * x342 - 2 * x343) - x367 + x375
    x377 = x231 * x50 + x248
    x378 = -x0 * (x328 * x57 - x330 * x73) + x359 + x377
    x379 = -x17 * x378
    x380 = x362 * x50 + x379
    x381 = -x17 * x380
    x382 = x365 * x50 + x381
    x383 = -x0 * x362 + x356 + x382
    x384 = -x17 * x383
    x385 = x373 * x50 + x384
    x386 = -x0 * (2 * x358 - 2 * x363) + x367 + x385
    x387 = -x17 * x386
    x388 = x376 * x50 + x387
    x389 = (
        x0 * (3 * x335 - x341 + 3 * x346 - 3 * x355)
        - x0 * (x341 - 3 * x356 + 3 * x357 - 3 * x366)
        + x388
    )
    x390 = x0 * (x190 + x21)
    x391 = x0 * (x194 + x25)
    x392 = x258 * x50
    x393 = x265 + x392
    x394 = x390 - x391 + x393
    x395 = x0 * x394
    x396 = x0 * (x202 + x31)
    x397 = x263 * x50
    x398 = x274 + x397
    x399 = x391 - x396 + x398
    x400 = x0 * x399
    x401 = 3 * x400
    x402 = x394 * x50
    x403 = x17 * x399
    x404 = -x403
    x405 = x402 + x404
    x406 = x405 * x50
    x407 = x399 * x50
    x408 = x0 * (x217 + x41)
    x409 = x272 * x50
    x410 = x290 + x409
    x411 = x396 - x408 + x410
    x412 = x17 * x411
    x413 = -x412
    x414 = x407 + x413
    x415 = x17 * x414
    x416 = x0 * x411
    x417 = x414 * x50
    x418 = x411 * x50
    x419 = x0 * (x231 + x59)
    x420 = x288 * x50
    x421 = x305 + x420
    x422 = x408 - x419 + x421
    x423 = x17 * x422
    x424 = -x423
    x425 = x418 + x424
    x426 = x17 * x425
    x427 = 2 * x402 - 2 * x403
    x428 = 2 * x407 - 2 * x412
    x429 = x0 * x428
    x430 = -x415
    x431 = x406 + x430
    x432 = x395 - x400 + x431
    x433 = -x426
    x434 = x417 + x433
    x435 = x400 - x416 + x434
    x436 = -x17 * x435
    x437 = x432 * x50 + x436
    x438 = x0 * x427 - x429 + x437
    x439 = 2 * x418 - 2 * x423
    x440 = x422 * x50
    x441 = x0 * (x247 + x75)
    x442 = x303 * x50
    x443 = x322 + x442
    x444 = x17 * (x419 - x441 + x443)
    x445 = -x444
    x446 = x440 + x445
    x447 = -x17 * x446
    x448 = x425 * x50 + x447
    x449 = -x0 * x422 + x416 + x448
    x450 = -x17 * x449
    x451 = x435 * x50 + x450
    x452 = -x0 * x439 + x429 + x451
    x453 = -x17 * x452
    x454 = x438 * x50 + x453
    x455 = (
        x0 * (3 * x395 - x401 + 3 * x406 - 3 * x415)
        - x0 * (x401 - 3 * x416 + 3 * x417 - 3 * x426)
        + x454
    )
    x456 = 2 * x256
    x457 = 2 * x253 - 2 * x26
    x458 = 2 * x255 - x456 + x457
    x459 = x0 * x458
    x460 = 2 * x261
    x461 = 2 * x259 - 2 * x32
    x462 = x456 - x460 + x461
    x463 = x0 * x462
    x464 = x28 * x50
    x465 = x36 + x464
    x466 = x459 - x463 + x465
    x467 = x0 * x466
    x468 = 2 * x270
    x469 = 2 * x268 - 2 * x42
    x470 = x460 - x468 + x469
    x471 = x0 * x470
    x472 = x34 * x50
    x473 = x46 + x472
    x474 = x463 - x471 + x473
    x475 = x0 * x474
    x476 = 3 * x475
    x477 = x466 * x50
    x478 = x17 * x474
    x479 = -x478
    x480 = x477 + x479
    x481 = x480 * x50
    x482 = x474 * x50
    x483 = 2 * x286
    x484 = 2 * x284 - 2 * x60
    x485 = x468 - x483 + x484
    x486 = x0 * x485
    x487 = x44 * x50
    x488 = x487 + x64
    x489 = x471 - x486 + x488
    x490 = x17 * x489
    x491 = -x490
    x492 = x482 + x491
    x493 = x17 * x492
    x494 = x0 * x489
    x495 = x492 * x50
    x496 = x489 * x50
    x497 = 2 * x301
    x498 = 2 * x299 - 2 * x76
    x499 = x483 - x497 + x498
    x500 = x0 * x499
    x501 = x50 * x62
    x502 = x501 + x80
    x503 = x486 - x500 + x502
    x504 = x17 * x503
    x505 = -x504
    x506 = x496 + x505
    x507 = x17 * x506
    x508 = x0 * (2 * x482 - 2 * x490)
    x509 = -x493
    x510 = x481 + x509
    x511 = x467 - x475 + x510
    x512 = -x507
    x513 = x495 + x512
    x514 = x475 - x494 + x513
    x515 = -x17 * x514
    x516 = x50 * x511 + x515
    x517 = x0 * (2 * x477 - 2 * x478) - x508 + x516
    x518 = x50 * x503
    x519 = x0 * (-2 * x100 - 2 * x318 + 2 * x319 + x497)
    x520 = x50 * x78
    x521 = x103 + x520
    x522 = x17 * (x500 - x519 + x521)
    x523 = -x522
    x524 = x518 + x523
    x525 = -x17 * x524
    x526 = x50 * x506 + x525
    x527 = -x0 * x503 + x494 + x526
    x528 = -x17 * x527
    x529 = x50 * x514 + x528
    x530 = -x0 * (2 * x496 - 2 * x504) + x508 + x529
    x531 = -x17 * x530
    x532 = x50 * x517 + x531
    x533 = (
        x0 * (3 * x467 - x476 + 3 * x481 - 3 * x493)
        - x0 * (x476 - 3 * x494 + 3 * x495 - 3 * x507)
        + x532
    )
    x534 = x0 * x54
    x535 = x0 * x68
    x536 = 2 * x535
    x537 = x37 * x6 + x53
    x538 = x50 * x537
    x539 = x47 * x6 + x67
    x540 = x17 * x539
    x541 = 2 * x534 - x536 + 2 * x538 - 2 * x540
    x542 = x0 * x84
    x543 = 2 * x542
    x544 = x50 * x539
    x545 = x6 * x65 + x83
    x546 = x17 * x545
    x547 = x536 - x543 + 2 * x544 - 2 * x546
    x548 = x0 * (x537 + x91)
    x549 = x0 * (x539 + x94)
    x550 = -x540
    x551 = x534 - x535
    x552 = x50 * (x538 + x550 + x551)
    x553 = -x546
    x554 = x535 - x542
    x555 = x544 + x553 + x554
    x556 = x17 * x555
    x557 = x0 * (x108 + x545)
    x558 = x50 * x555
    x559 = x0 * x105
    x560 = x50 * x545
    x561 = x17 * (x104 + x6 * x81)
    x562 = -x561
    x563 = x17 * (x542 - x559 + x560 + x562)
    x564 = (
        -x0 * (x111 + x547)
        + x0 * (x541 + x97)
        - x17 * (x549 - x557 + x558 - x563)
        + x50 * (x548 - x549 + x552 - x556)
    )
    x565 = x0 * x142
    x566 = x0 * x154
    x567 = 2 * x566
    x568 = x128 * x6 + x141
    x569 = x50 * x568
    x570 = x136 * x6 + x153
    x571 = x17 * x570
    x572 = 2 * x565 - x567 + 2 * x569 - 2 * x571
    x573 = x0 * x168
    x574 = 2 * x573
    x575 = x50 * x570
    x576 = x151 * x6 + x167
    x577 = x17 * x576
    x578 = x567 - x574 + 2 * x575 - 2 * x577
    x579 = x0 * (x173 + x568)
    x580 = x0 * (x174 + x570)
    x581 = -x571
    x582 = x50 * (x565 - x566 + x569 + x581)
    x583 = -x577
    x584 = x566 - x573 + x575 + x583
    x585 = x17 * x584
    x586 = x0 * (x183 + x576)
    x587 = x50 * x584
    x588 = x0 * x182
    x589 = x50 * x576
    x590 = x17 * (x165 * x6 + x181)
    x591 = -x590
    x592 = x17 * (x573 - x588 + x589 + x591)
    x593 = x0 * x211
    x594 = x0 * x223
    x595 = 2 * x594
    x596 = x197 * x6 + x210
    x597 = x50 * x596
    x598 = x205 * x6 + x222
    x599 = x17 * x598
    x600 = 2 * x593 - x595 + 2 * x597 - 2 * x599
    x601 = x0 * x237
    x602 = 2 * x601
    x603 = x50 * x598
    x604 = x220 * x6 + x236
    x605 = x17 * x604
    x606 = x595 - x602 + 2 * x603 - 2 * x605
    x607 = x0 * (x242 + x596)
    x608 = x0 * (x243 + x598)
    x609 = -x599
    x610 = x50 * (x593 - x594 + x597 + x609)
    x611 = -x605
    x612 = x594 - x601 + x603 + x611
    x613 = x17 * x612
    x614 = x0 * (x251 + x604)
    x615 = x50 * x612
    x616 = x0 * x250
    x617 = x50 * x604
    x618 = x17 * (x234 * x6 + x249)
    x619 = -x618
    x620 = x17 * (x601 - x616 + x617 + x619)
    x621 = x0 * x281
    x622 = x0 * x294
    x623 = 2 * x622
    x624 = x266 * x6 + x280
    x625 = x50 * x624
    x626 = x275 * x6 + x293
    x627 = x17 * x626
    x628 = 2 * x621 - x623 + 2 * x625 - 2 * x627
    x629 = x0 * x309
    x630 = 2 * x629
    x631 = x50 * x626
    x632 = x291 * x6 + x308
    x633 = x17 * x632
    x634 = x623 - x630 + 2 * x631 - 2 * x633
    x635 = x0 * (x314 + x624)
    x636 = x0 * (x315 + x626)
    x637 = -x627
    x638 = x50 * (x621 - x622 + x625 + x637)
    x639 = -x633
    x640 = x622 - x629 + x631 + x639
    x641 = x17 * x640
    x642 = x0 * (x325 + x632)
    x643 = x50 * x640
    x644 = x0 * x324
    x645 = x50 * x632
    x646 = x17 * (x306 * x6 + x323)
    x647 = -x646
    x648 = x17 * (x629 - x644 + x645 + x647)
    x649 = x0 * x345
    x650 = x0 * x354
    x651 = 2 * x650
    x652 = x334 * x6 + x344
    x653 = x50 * x652
    x654 = x339 * x6 + x353
    x655 = x17 * x654
    x656 = 2 * x649 - x651 + 2 * x653 - 2 * x655
    x657 = x0 * x365
    x658 = 2 * x657
    x659 = x50 * x654
    x660 = x351 * x6 + x364
    x661 = x17 * x660
    x662 = x651 - x658 + 2 * x659 - 2 * x661
    x663 = x0 * (x370 + x652)
    x664 = x0 * (x373 + x654)
    x665 = -x655
    x666 = x649 - x650
    x667 = x50 * (x653 + x665 + x666)
    x668 = -x661
    x669 = x650 - x657
    x670 = x659 + x668 + x669
    x671 = x17 * x670
    x672 = x0 * (x383 + x660)
    x673 = x50 * x670
    x674 = x0 * x380
    x675 = x50 * x660
    x676 = x17 * (x362 * x6 + x379)
    x677 = -x676
    x678 = x17 * (x657 - x674 + x675 + x677)
    x679 = x0 * x405
    x680 = x0 * x414
    x681 = 2 * x680
    x682 = x394 * x6 + x404
    x683 = x50 * x682
    x684 = x399 * x6 + x413
    x685 = x17 * x684
    x686 = 2 * x679 - x681 + 2 * x683 - 2 * x685
    x687 = x0 * x425
    x688 = 2 * x687
    x689 = x50 * x684
    x690 = x411 * x6 + x424
    x691 = x17 * x690
    x692 = x681 - x688 + 2 * x689 - 2 * x691
    x693 = x0 * (x432 + x682)
    x694 = x0 * (x435 + x684)
    x695 = -x685
    x696 = x679 - x680
    x697 = x50 * (x683 + x695 + x696)
    x698 = -x691
    x699 = x680 - x687
    x700 = x689 + x698 + x699
    x701 = x17 * x700
    x702 = x0 * (x449 + x690)
    x703 = x50 * x700
    x704 = x0 * x446
    x705 = x50 * x690
    x706 = x17 * (x422 * x6 + x445)
    x707 = -x706
    x708 = x17 * (x687 - x704 + x705 + x707)
    x709 = x0 * x480
    x710 = x0 * x492
    x711 = 2 * x710
    x712 = x466 * x6 + x479
    x713 = x50 * x712
    x714 = x474 * x6 + x491
    x715 = x17 * x714
    x716 = 2 * x709 - x711 + 2 * x713 - 2 * x715
    x717 = x0 * x506
    x718 = 2 * x717
    x719 = x50 * x714
    x720 = x489 * x6 + x505
    x721 = x17 * x720
    x722 = x711 - x718 + 2 * x719 - 2 * x721
    x723 = x0 * (x511 + x712)
    x724 = x0 * (x514 + x714)
    x725 = -x715
    x726 = x709 - x710
    x727 = x50 * (x713 + x725 + x726)
    x728 = -x721
    x729 = x710 - x717
    x730 = x719 + x728 + x729
    x731 = x17 * x730
    x732 = x0 * (x527 + x720)
    x733 = x50 * x730
    x734 = x0 * x524
    x735 = x50 * x720
    x736 = x17 * (x503 * x6 + x523)
    x737 = -x736
    x738 = x17 * (x717 - x734 + x735 + x737)
    x739 = 3 * x463
    x740 = x0 * (-3 * x35 + 3 * x459 + 3 * x464 - x739)
    x741 = 3 * x471
    x742 = x0 * (-3 * x45 + 3 * x472 + x739 - x741)
    x743 = x54 + x740 - x742
    x744 = x50 * x743
    x745 = 3 * x486
    x746 = x0 * (3 * x487 - 3 * x63 + x741 - x745)
    x747 = x68 + x742 - x746
    x748 = x17 * x747
    x749 = x50 * x747
    x750 = 3 * x500
    x751 = x0 * (3 * x501 + x745 - x750 - 3 * x79)
    x752 = x746 - x751 + x84
    x753 = x17 * x752
    x754 = x0 * x747
    x755 = -x748
    x756 = x744 + x755
    x757 = -x753
    x758 = x749 + x757
    x759 = x0 * x743 - x17 * x758 + x50 * x756 - x754
    x760 = -x17 * (-x0 * (-3 * x102 - 3 * x519 + 3 * x520 + x750) + x105 + x751)
    x761 = x50 * x752 + x760
    x762 = -x0 * x752 - x17 * x761 + x50 * x758 + x754
    x763 = x0 * x189
    x764 = x0 * x193
    x765 = 2 * x764
    x766 = x121 * x50
    x767 = -2 * x126 + 2 * x763 - x765 + 2 * x766
    x768 = x0 * x767
    x769 = x0 * x201
    x770 = 2 * x769
    x771 = x125 * x50
    x772 = -2 * x134 + x765 - x770 + 2 * x771
    x773 = x0 * x772
    x774 = x142 + x768 - x773
    x775 = x50 * x774
    x776 = x0 * x216
    x777 = 2 * x776
    x778 = x133 * x50
    x779 = -2 * x149 + x770 - x777 + 2 * x778
    x780 = x0 * x779
    x781 = x154 + x773 - x780
    x782 = x17 * x781
    x783 = x50 * x781
    x784 = x0 * x230
    x785 = 2 * x784
    x786 = x148 * x50
    x787 = -2 * x163 + x777 - x785 + 2 * x786
    x788 = x0 * x787
    x789 = x168 + x780 - x788
    x790 = x17 * x789
    x791 = x0 * x781
    x792 = -x782
    x793 = x775 + x792
    x794 = -x790
    x795 = x783 + x794
    x796 = x0 * x774 - x17 * x795 + x50 * x793 - x791
    x797 = x0 * x246
    x798 = x162 * x50
    x799 = -x17 * (-x0 * (-2 * x179 + x785 - 2 * x797 + 2 * x798) + x182 + x788)
    x800 = x50 * x789 + x799
    x801 = -x0 * x789 - x17 * x800 + x50 * x795 + x791
    x802 = x0 * x260
    x803 = 2 * x802
    x804 = x0 * x254
    x805 = -2 * x35 + 2 * x464
    x806 = -x803 + 2 * x804 + x805
    x807 = x0 * x806
    x808 = x0 * x269
    x809 = 2 * x808
    x810 = -2 * x45 + 2 * x472
    x811 = x803 - x809 + x810
    x812 = x0 * x811
    x813 = x54 + x807 - x812
    x814 = x50 * x813
    x815 = x0 * x285
    x816 = 2 * x815
    x817 = 2 * x487 - 2 * x63
    x818 = x809 - x816 + x817
    x819 = x0 * x818
    x820 = x68 + x812 - x819
    x821 = x17 * x820
    x822 = x50 * x820
    x823 = x0 * x300
    x824 = 2 * x823
    x825 = 2 * x501 - 2 * x79
    x826 = x816 - x824 + x825
    x827 = x0 * x826
    x828 = x819 - x827 + x84
    x829 = x17 * x828
    x830 = x0 * x820
    x831 = -x821
    x832 = x814 + x831
    x833 = -x829
    x834 = x822 + x833
    x835 = x0 * x813 - x17 * x834 + x50 * x832 - x830
    x836 = x0 * x320
    x837 = -x17 * (-x0 * (-2 * x102 + 2 * x520 + x824 - 2 * x836) + x105 + x827)
    x838 = x50 * x828 + x837
    x839 = -x0 * x828 - x17 * x838 + x50 * x834 + x830
    x840 = x0 * x333
    x841 = x0 * x338
    x842 = x211 + x840 - x841
    x843 = x50 * x842
    x844 = x0 * x350
    x845 = x223 + x841 - x844
    x846 = x17 * x845
    x847 = x50 * x845
    x848 = x0 * x361
    x849 = x237 + x844 - x848
    x850 = x17 * x849
    x851 = x0 * x845
    x852 = -x846
    x853 = x843 + x852
    x854 = -x850
    x855 = x847 + x854
    x856 = x0 * x842 - x17 * x855 + x50 * x853 - x851
    x857 = -x17 * (-x0 * x377 + x250 + x848)
    x858 = x50 * x849 + x857
    x859 = -x0 * x849 - x17 * x858 + x50 * x855 + x851
    x860 = x0 * x393
    x861 = x0 * x398
    x862 = x281 + x860 - x861
    x863 = x50 * x862
    x864 = x0 * x410
    x865 = x294 + x861 - x864
    x866 = x17 * x865
    x867 = x50 * x865
    x868 = x0 * x421
    x869 = x309 + x864 - x868
    x870 = x17 * x869
    x871 = x0 * x865
    x872 = -x866
    x873 = x863 + x872
    x874 = -x870
    x875 = x867 + x874
    x876 = x0 * x862 - x17 * x875 + x50 * x873 - x871
    x877 = -x17 * (-x0 * x443 + x324 + x868)
    x878 = x50 * x869 + x877
    x879 = -x0 * x869 - x17 * x878 + x50 * x875 + x871
    x880 = x0 * x465
    x881 = x0 * x473
    x882 = x54 + x880 - x881
    x883 = x50 * x882
    x884 = x0 * x488
    x885 = x68 + x881 - x884
    x886 = x17 * x885
    x887 = x50 * x885
    x888 = x0 * x502
    x889 = x84 + x884 - x888
    x890 = x17 * x889
    x891 = x0 * x885
    x892 = -x886
    x893 = x883 + x892
    x894 = -x890
    x895 = x887 + x894
    x896 = x0 * x882 - x17 * x895 + x50 * x893 - x891
    x897 = -x17 * (-x0 * x521 + x105 + x888)
    x898 = x50 * x889 + x897
    x899 = -x0 * x889 - x17 * x898 + x50 * x895 + x891
    x900 = -x17 * x372 + x369 * x50 + x666
    x901 = -x17 * x382 + x372 * x50 + x669
    x902 = -x17 * x434 + x431 * x50 + x696
    x903 = -x17 * x448 + x434 * x50 + x699
    x904 = -x17 * x513 + x50 * x510 + x726
    x905 = -x17 * x526 + x50 * x513 + x729
    x906 = -x17 * x93 + x50 * x90 + x551
    x907 = -x107 * x17 + x50 * x93 + x554
    x908 = x537 * x6 + x550
    x909 = x50 * x908
    x910 = x539 * x6 + x553
    x911 = x17 * x910
    x912 = -x911
    x913 = x0 * x541
    x914 = x0 * x547
    x915 = x50 * x910
    x916 = x17 * (x545 * x6 + x562)
    x917 = -x916
    x918 = x0 * (x543 - 2 * x559 + 2 * x560 - 2 * x561)
    x919 = 2 * x549
    x920 = (
        x0 * (2 * x548 + 2 * x552 - 2 * x556 + x908 - x919)
        - x0 * (-2 * x557 + 2 * x558 - 2 * x563 + x910 + x919)
        - x17 * (x914 + x915 + x917 - x918)
        + x50 * (x909 + x912 + x913 - x914)
    )
    x921 = x568 * x6 + x581
    x922 = x50 * x921
    x923 = x570 * x6 + x583
    x924 = x17 * x923
    x925 = -x924
    x926 = x0 * x572
    x927 = x0 * x578
    x928 = x50 * x923
    x929 = x17 * (x576 * x6 + x591)
    x930 = -x929
    x931 = x0 * (x574 - 2 * x588 + 2 * x589 - 2 * x590)
    x932 = 2 * x580
    x933 = x596 * x6 + x609
    x934 = x50 * x933
    x935 = x598 * x6 + x611
    x936 = x17 * x935
    x937 = -x936
    x938 = x0 * x600
    x939 = x0 * x606
    x940 = x50 * x935
    x941 = x17 * (x6 * x604 + x619)
    x942 = -x941
    x943 = x0 * (x602 - 2 * x616 + 2 * x617 - 2 * x618)
    x944 = 2 * x608
    x945 = x6 * x624 + x637
    x946 = x50 * x945
    x947 = x6 * x626 + x639
    x948 = x17 * x947
    x949 = -x948
    x950 = x0 * x628
    x951 = x0 * x634
    x952 = x50 * x947
    x953 = x17 * (x6 * x632 + x647)
    x954 = -x953
    x955 = x0 * (x630 - 2 * x644 + 2 * x645 - 2 * x646)
    x956 = 2 * x636
    x957 = x6 * x652 + x665
    x958 = x50 * x957
    x959 = x6 * x654 + x668
    x960 = x17 * x959
    x961 = -x960
    x962 = x0 * x656
    x963 = x0 * x662
    x964 = x50 * x959
    x965 = x17 * (x6 * x660 + x677)
    x966 = -x965
    x967 = x0 * (x658 - 2 * x674 + 2 * x675 - 2 * x676)
    x968 = 2 * x664
    x969 = x6 * x682 + x695
    x970 = x50 * x969
    x971 = x6 * x684 + x698
    x972 = x17 * x971
    x973 = -x972
    x974 = x0 * x686
    x975 = x0 * x692
    x976 = x50 * x971
    x977 = x17 * (x6 * x690 + x707)
    x978 = -x977
    x979 = x0 * (x688 - 2 * x704 + 2 * x705 - 2 * x706)
    x980 = 2 * x694
    x981 = x6 * x712 + x725
    x982 = x50 * x981
    x983 = x6 * x714 + x728
    x984 = x17 * x983
    x985 = -x984
    x986 = x0 * x716
    x987 = x0 * x722
    x988 = x50 * x983
    x989 = x17 * (x6 * x720 + x737)
    x990 = -x989
    x991 = x0 * (x718 - 2 * x734 + 2 * x735 - 2 * x736)
    x992 = 2 * x724
    x993 = x6 * x743 + x755
    x994 = x6 * x747 + x757
    x995 = x0 * x756
    x996 = x0 * x758
    x997 = x50 * x993
    x998 = x17 * x994
    x999 = -x998
    x1000 = x0 * x761
    x1001 = x50 * x994
    x1002 = x17 * (x6 * x752 + x760)
    x1003 = -x1002
    x1004 = x6 * x774 + x792
    x1005 = x6 * x781 + x794
    x1006 = x0 * x793
    x1007 = x0 * x795
    x1008 = x1004 * x50
    x1009 = x1005 * x17
    x1010 = -x1009
    x1011 = x0 * x800
    x1012 = x1005 * x50
    x1013 = x17 * (x6 * x789 + x799)
    x1014 = -x1013
    x1015 = x6 * x813 + x831
    x1016 = x6 * x820 + x833
    x1017 = x0 * x832
    x1018 = x0 * x834
    x1019 = x1015 * x50
    x1020 = x1016 * x17
    x1021 = -x1020
    x1022 = x0 * x838
    x1023 = x1016 * x50
    x1024 = x17 * (x6 * x828 + x837)
    x1025 = -x1024
    x1026 = x6 * x842 + x852
    x1027 = x6 * x845 + x854
    x1028 = x0 * x853
    x1029 = x0 * x855
    x1030 = x1026 * x50
    x1031 = x1027 * x17
    x1032 = -x1031
    x1033 = x0 * x858
    x1034 = x1027 * x50
    x1035 = x17 * (x6 * x849 + x857)
    x1036 = -x1035
    x1037 = x6 * x862 + x872
    x1038 = x6 * x865 + x874
    x1039 = x0 * x873
    x1040 = x0 * x875
    x1041 = x1037 * x50
    x1042 = x1038 * x17
    x1043 = -x1042
    x1044 = x0 * x878
    x1045 = x1038 * x50
    x1046 = x17 * (x6 * x869 + x877)
    x1047 = -x1046
    x1048 = x6 * x882 + x892
    x1049 = x6 * x885 + x894
    x1050 = x0 * x893
    x1051 = x0 * x895
    x1052 = x1048 * x50
    x1053 = x1049 * x17
    x1054 = -x1053
    x1055 = x0 * x898
    x1056 = x1049 * x50
    x1057 = x17 * (x6 * x889 + x897)
    x1058 = -x1057
    x1059 = x345 * x6 + x368
    x1060 = x354 * x6 + x371
    x1061 = x0 * x369
    x1062 = x0 * x372
    x1063 = x1059 * x50
    x1064 = x1060 * x17
    x1065 = -x1064
    x1066 = x0 * x382
    x1067 = x1060 * x50
    x1068 = x17 * (x365 * x6 + x381)
    x1069 = -x1068
    x1070 = x405 * x6 + x430
    x1071 = x414 * x6 + x433
    x1072 = x0 * x431
    x1073 = x0 * x434
    x1074 = x1070 * x50
    x1075 = x1071 * x17
    x1076 = -x1075
    x1077 = x0 * x448
    x1078 = x1071 * x50
    x1079 = x17 * (x425 * x6 + x447)
    x1080 = -x1079
    x1081 = x480 * x6 + x509
    x1082 = x492 * x6 + x512
    x1083 = x0 * x510
    x1084 = x0 * x513
    x1085 = x1081 * x50
    x1086 = x1082 * x17
    x1087 = -x1086
    x1088 = x0 * x526
    x1089 = x1082 * x50
    x1090 = x17 * (x506 * x6 + x525)
    x1091 = -x1090
    x1092 = x54 * x6 + x89
    x1093 = x6 * x68 + x92
    x1094 = x0 * x90
    x1095 = x0 * x93
    x1096 = x1092 * x50
    x1097 = x1093 * x17
    x1098 = -x1097
    x1099 = x0 * x107
    x1100 = x1093 * x50
    x1101 = x17 * (x106 + x6 * x84)
    x1102 = -x1101
    x1103 = 2 * x396
    x1104 = 2 * x391
    x1105 = -2 * x273 + 2 * x397
    x1106 = x0 * (-x1103 + x1104 + x1105 + x34)
    x1107 = 3 * x1106
    x1108 = -2 * x264 + 2 * x392
    x1109 = x0 * (-x1104 + x1108 + x28 + 2 * x390)
    x1110 = -x1107 + 3 * x1109 + 3 * x477 - 3 * x478
    x1111 = x0 * (x1110 + x37)
    x1112 = 2 * x408
    x1113 = -2 * x289 + 2 * x409
    x1114 = x0 * (x1103 - x1112 + x1113 + x44)
    x1115 = 3 * x1114
    x1116 = x1107 - x1115 + 3 * x482 - 3 * x490
    x1117 = x0 * (x1116 + x47)
    x1118 = x1111 - x1117 + x756
    x1119 = 2 * x419
    x1120 = -2 * x304 + 2 * x420
    x1121 = x0 * (x1112 - x1119 + x1120 + x62)
    x1122 = 3 * x1121
    x1123 = x1115 - x1122 + 3 * x496 - 3 * x504
    x1124 = x0 * (x1123 + x65)
    x1125 = x1117 - x1124 + x758
    x1126 = -x1125 * x17
    x1127 = x1118 * x50 + x1126
    x1128 = x0 * (x1119 - 2 * x321 - 2 * x441 + 2 * x442 + x78)
    x1129 = -x17 * (-x0 * (x1122 - 3 * x1128 + 3 * x518 - 3 * x522 + x81) + x1124 + x761)
    x1130 = x1125 * x50 + x1129
    x1131 = x193 * x50
    x1132 = x17 * x201
    x1133 = x1131 - x1132 + x262
    x1134 = x0 * (x1133 + x125)
    x1135 = 2 * x1134
    x1136 = x135 + x764 - x769 + x771
    x1137 = x1136 * x17
    x1138 = x189 * x50
    x1139 = x17 * x193
    x1140 = x1138 - x1139 + x257
    x1141 = x0 * (x1140 + x121)
    x1142 = x50 * (x127 + x763 - x764 + x766)
    x1143 = x0 * (-x1135 - 2 * x1137 + 2 * x1141 + 2 * x1142 + x128)
    x1144 = x201 * x50
    x1145 = x17 * x216
    x1146 = x1144 - x1145 + x271
    x1147 = x0 * (x1146 + x133)
    x1148 = 2 * x1147
    x1149 = x150 + x769 - x776 + x778
    x1150 = x1149 * x17
    x1151 = x1136 * x50
    x1152 = x0 * (x1135 - x1148 - 2 * x1150 + 2 * x1151 + x136)
    x1153 = x1143 - x1152 + x793
    x1154 = x216 * x50
    x1155 = x17 * x230
    x1156 = x1154 - x1155 + x287
    x1157 = x0 * (x1156 + x148)
    x1158 = 2 * x1157
    x1159 = x164 + x776 - x784 + x786
    x1160 = x1159 * x17
    x1161 = x1149 * x50
    x1162 = x0 * (x1148 - x1158 - 2 * x1160 + 2 * x1161 + x151)
    x1163 = x1152 - x1162 + x795
    x1164 = -x1163 * x17
    x1165 = x1153 * x50 + x1164
    x1166 = -x17 * x246 + x230 * x50 + x302
    x1167 = x0 * (x1166 + x162)
    x1168 = x17 * (x180 + x784 - x797 + x798)
    x1169 = x1159 * x50
    x1170 = -x17 * (
        -x0 * (x1158 - 2 * x1167 - 2 * x1168 + 2 * x1169 + x165) + x1162 + x800
    )
    x1171 = x1163 * x50 + x1170
    x1172 = x0 * x25
    x1173 = x0 * x31
    x1174 = x260 * x50
    x1175 = x17 * x269
    x1176 = x1172 - x1173 + x1174 - x1175
    x1177 = x0 * (x1176 + x34)
    x1178 = 2 * x1177
    x1179 = x473 + x802 - x808
    x1180 = x1179 * x17
    x1181 = x0 * x21
    x1182 = x254 * x50
    x1183 = x17 * x260
    x1184 = -x1172 + x1181 + x1182 - x1183
    x1185 = x0 * (x1184 + x28)
    x1186 = x50 * (x465 - x802 + x804)
    x1187 = x0 * (-x1178 - 2 * x1180 + 2 * x1185 + 2 * x1186 + x37)
    x1188 = x0 * x41
    x1189 = x269 * x50
    x1190 = x17 * x285
    x1191 = x1173 - x1188 + x1189 - x1190
    x1192 = x0 * (x1191 + x44)
    x1193 = 2 * x1192
    x1194 = x488 + x808 - x815
    x1195 = x1194 * x17
    x1196 = x1179 * x50
    x1197 = x0 * (x1178 - x1193 - 2 * x1195 + 2 * x1196 + x47)
    x1198 = x1187 - x1197 + x832
    x1199 = x0 * x59
    x1200 = x285 * x50
    x1201 = x17 * x300
    x1202 = x1188 - x1199 + x1200 - x1201
    x1203 = x0 * (x1202 + x62)
    x1204 = 2 * x1203
    x1205 = x502 + x815 - x823
    x1206 = x1205 * x17
    x1207 = x1194 * x50
    x1208 = x0 * (x1193 - x1204 - 2 * x1206 + 2 * x1207 + x65)
    x1209 = x1197 - x1208 + x834
    x1210 = -x1209 * x17
    x1211 = x1198 * x50 + x1210
    x1212 = -x0 * x75 + x1199 - x17 * x320 + x300 * x50
    x1213 = x0 * (x1212 + x78)
    x1214 = x17 * (x521 + x823 - x836)
    x1215 = x1205 * x50
    x1216 = -x17 * (
        -x0 * (x1204 - 2 * x1213 - 2 * x1214 + 2 * x1215 + x81) + x1208 + x838
    )
    x1217 = x1209 * x50 + x1216
    x1218 = x0 * x190
    x1219 = x0 * x194
    x1220 = x333 * x50
    x1221 = x17 * x338
    x1222 = x1218 - x1219 + x1220 - x1221
    x1223 = x0 * (x1222 + x197)
    x1224 = x0 * x202
    x1225 = x338 * x50
    x1226 = x17 * x350
    x1227 = x1219 - x1224 + x1225 - x1226
    x1228 = x0 * (x1227 + x205)
    x1229 = x1223 - x1228 + x853
    x1230 = x0 * x217
    x1231 = x350 * x50
    x1232 = x17 * x361
    x1233 = x1224 - x1230 + x1231 - x1232
    x1234 = x0 * (x1233 + x220)
    x1235 = x1228 - x1234 + x855
    x1236 = -x1235 * x17
    x1237 = x1229 * x50 + x1236
    x1238 = -x0 * x231 + x1230 - x17 * x377 + x361 * x50
    x1239 = -x17 * (-x0 * (x1238 + x234) + x1234 + x858)
    x1240 = x1235 * x50 + x1239
    x1241 = x0 * x258
    x1242 = x0 * x263
    x1243 = x393 * x50
    x1244 = x17 * x398
    x1245 = x1241 - x1242 + x1243 - x1244
    x1246 = x0 * (x1245 + x266)
    x1247 = x0 * x272
    x1248 = x398 * x50
    x1249 = x17 * x410
    x1250 = x1242 - x1247 + x1248 - x1249
    x1251 = x0 * (x1250 + x275)
    x1252 = x1246 - x1251 + x873
    x1253 = x0 * x288
    x1254 = x410 * x50
    x1255 = x17 * x421
    x1256 = x1247 - x1253 + x1254 - x1255
    x1257 = x0 * (x1256 + x291)
    x1258 = x1251 - x1257 + x875
    x1259 = -x1258 * x17
    x1260 = x1252 * x50 + x1259
    x1261 = -x0 * x303 + x1253 - x17 * x443 + x421 * x50
    x1262 = -x17 * (-x0 * (x1261 + x306) + x1257 + x878)
    x1263 = x1258 * x50 + x1262
    x1264 = x0 * x28
    x1265 = x0 * x34
    x1266 = x465 * x50
    x1267 = x17 * x473
    x1268 = x1264 - x1265 + x1266 - x1267
    x1269 = x0 * (x1268 + x37)
    x1270 = x0 * x44
    x1271 = x473 * x50
    x1272 = x17 * x488
    x1273 = x1265 - x1270 + x1271 - x1272
    x1274 = x0 * (x1273 + x47)
    x1275 = x1269 - x1274 + x893
    x1276 = x0 * x62
    x1277 = x488 * x50
    x1278 = x17 * x502
    x1279 = x1270 - x1276 + x1277 - x1278
    x1280 = x0 * (x1279 + x65)
    x1281 = x1274 - x1280 + x895
    x1282 = -x1281 * x17
    x1283 = x1275 * x50 + x1282
    x1284 = -x0 * x78 + x1276 - x17 * x521 + x50 * x502
    x1285 = -x17 * (-x0 * (x1284 + x81) + x1280 + x898)
    x1286 = x1281 * x50 + x1285
    x1287 = x6 * x908 + x912
    x1288 = -x17 * (x6 * x910 + x917)
    x1289 = 3 * x914
    x1290 = (
        x0 * (-x1289 + 3 * x909 - 3 * x911 + 3 * x913)
        - x0 * (x1289 + 3 * x915 - 3 * x916 - 3 * x918)
        + x1287 * x50
        + x1288
    )
    x1291 = x6 * x921 + x925
    x1292 = -x17 * (x6 * x923 + x930)
    x1293 = 3 * x927
    x1294 = x6 * x933 + x937
    x1295 = -x17 * (x6 * x935 + x942)
    x1296 = 3 * x939
    x1297 = x6 * x945 + x949
    x1298 = -x17 * (x6 * x947 + x954)
    x1299 = 3 * x951
    x1300 = x6 * x957 + x961
    x1301 = -x17 * (x6 * x959 + x966)
    x1302 = 3 * x963
    x1303 = x6 * x969 + x973
    x1304 = -x17 * (x6 * x971 + x978)
    x1305 = 3 * x975
    x1306 = x6 * x981 + x985
    x1307 = -x17 * (x6 * x983 + x990)
    x1308 = 3 * x987
    x1309 = x6 * x993 + x999
    x1310 = -x17 * (x1003 + x6 * x994)
    x1311 = 2 * x996
    x1312 = x1004 * x6 + x1010
    x1313 = -x17 * (x1005 * x6 + x1014)
    x1314 = 2 * x1007
    x1315 = x1015 * x6 + x1021
    x1316 = -x17 * (x1016 * x6 + x1025)
    x1317 = 2 * x1018
    x1318 = x1026 * x6 + x1032
    x1319 = -x17 * (x1027 * x6 + x1036)
    x1320 = 2 * x1029
    x1321 = x1037 * x6 + x1043
    x1322 = -x17 * (x1038 * x6 + x1047)
    x1323 = 2 * x1040
    x1324 = x1048 * x6 + x1054
    x1325 = -x17 * (x1049 * x6 + x1058)
    x1326 = 2 * x1051
    x1327 = x1059 * x6 + x1065
    x1328 = -x17 * (x1060 * x6 + x1069)
    x1329 = 2 * x1062
    x1330 = x1070 * x6 + x1076
    x1331 = -x17 * (x1071 * x6 + x1080)
    x1332 = 2 * x1073
    x1333 = x1081 * x6 + x1087
    x1334 = -x17 * (x1082 * x6 + x1091)
    x1335 = 2 * x1084
    x1336 = x1092 * x6 + x1098
    x1337 = -x17 * (x1093 * x6 + x1102)
    x1338 = 2 * x1095
    x1339 = x1118 * x6 + x1126
    x1340 = -x17 * (x1125 * x6 + x1129)
    x1341 = x1153 * x6 + x1164
    x1342 = -x17 * (x1163 * x6 + x1170)
    x1343 = x1198 * x6 + x1210
    x1344 = -x17 * (x1209 * x6 + x1216)
    x1345 = x1229 * x6 + x1236
    x1346 = -x17 * (x1235 * x6 + x1239)
    x1347 = x1252 * x6 + x1259
    x1348 = -x17 * (x1258 * x6 + x1262)
    x1349 = x1275 * x6 + x1282
    x1350 = -x17 * (x1281 * x6 + x1285)
    x1351 = x370 * x6 + x374
    x1352 = -x17 * (x373 * x6 + x384)
    x1353 = x432 * x6 + x436
    x1354 = -x17 * (x435 * x6 + x450)
    x1355 = x511 * x6 + x515
    x1356 = -x17 * (x514 * x6 + x528)
    x1357 = x6 * x91 + x95
    x1358 = -x17 * (x109 + x6 * x94)
    x1359 = 2 * x471
    x1360 = x0 * (x351 + x470)
    x1361 = 2 * x1360
    x1362 = 2 * x463
    x1363 = x0 * (x339 + x462)
    x1364 = 2 * x1363
    x1365 = x0 * (-x1359 - x1361 + x1362 + x1364 + x428 + x810)
    x1366 = 3 * x1365
    x1367 = x1106 - x1114 + x492
    x1368 = x1367 * x17
    x1369 = 2 * x742
    x1370 = x0 * (x334 + x458)
    x1371 = x0 * (-x1362 - x1364 + 2 * x1370 + x427 + 2 * x459 + x805)
    x1372 = x50 * (-x1106 + x1109 + x480)
    x1373 = 2 * x486
    x1374 = x0 * (x362 + x485)
    x1375 = 2 * x1374
    x1376 = x0 * (x1359 + x1361 - x1373 - x1375 + x439 + x817)
    x1377 = 3 * x1376
    x1378 = x1114 - x1121 + x506
    x1379 = x1378 * x17
    x1380 = 2 * x746
    x1381 = x1367 * x50
    x1382 = x0 * (x1366 + x1369 - x1377 - 3 * x1379 - x1380 + 3 * x1381 + x87)
    x1383 = (
        x0 * (-x1366 - 3 * x1368 - x1369 + 3 * x1371 + 3 * x1372 + 2 * x740 + x86)
        + x1127
        - x1382
    )
    x1384 = x0 * (x378 + x499)
    x1385 = x0 * (x1373 + x1375 - 2 * x1384 + 2 * x440 - 2 * x444 - 2 * x500 + x825)
    x1386 = x17 * (x1121 - x1128 + x524)
    x1387 = x1378 * x50
    x1388 = -x17 * (
        -x0 * (x1377 + x1380 - 3 * x1385 - 3 * x1386 + 3 * x1387 - 2 * x751 + x98)
        + x1130
        + x1382
    )
    x1389 = x1383 * x50 + x1388
    x1390 = 2 * x773
    x1391 = x0 * (-2 * x123 + 2 * x192)
    x1392 = x0 * (-2 * x131 + 2 * x200)
    x1393 = x1133 * x50 - x1146 * x17 + x1391 - x1392
    x1394 = x0 * (x1393 + x772)
    x1395 = 2 * x1394
    x1396 = x1134 - x1147 - x1150 + x1151
    x1397 = x1396 * x17
    x1398 = x0 * (-2 * x119 + 2 * x188) - x1133 * x17 + x1140 * x50 - x1391
    x1399 = x0 * (x1398 + x767)
    x1400 = x50 * (-x1134 - x1137 + x1141 + x1142)
    x1401 = 2 * x780
    x1402 = x0 * (-2 * x146 + 2 * x215)
    x1403 = x1146 * x50 - x1156 * x17 + x1392 - x1402
    x1404 = x0 * (x1403 + x779)
    x1405 = 2 * x1404
    x1406 = x1147 - x1157 - x1160 + x1161
    x1407 = x1406 * x17
    x1408 = x1396 * x50
    x1409 = x0 * (x1390 + x1395 - x1401 - x1405 - 2 * x1407 + 2 * x1408 + x171)
    x1410 = (
        x0 * (-x1390 - x1395 - 2 * x1397 + 2 * x1399 + 2 * x1400 + x170 + 2 * x768)
        + x1165
        - x1409
    )
    x1411 = -x0 * (-2 * x160 + 2 * x229) + x1156 * x50 - x1166 * x17 + x1402
    x1412 = x0 * (x1411 + x787)
    x1413 = x17 * (x1157 - x1167 - x1168 + x1169)
    x1414 = x1406 * x50
    x1415 = -x17 * (
        -x0 * (x1401 + x1405 - 2 * x1412 - 2 * x1413 + 2 * x1414 + x176 - 2 * x788)
        + x1171
        + x1409
    )
    x1416 = x1410 * x50 + x1415
    x1417 = 2 * x812
    x1418 = x0 * x461
    x1419 = x0 * x469
    x1420 = x1176 * x50 - x1191 * x17 + x1418 - x1419
    x1421 = x0 * (x1420 + x811)
    x1422 = 2 * x1421
    x1423 = x1177 - x1192 - x1195 + x1196
    x1424 = x1423 * x17
    x1425 = x0 * x457 - x1176 * x17 + x1184 * x50 - x1418
    x1426 = x0 * (x1425 + x806)
    x1427 = x50 * (-x1177 - x1180 + x1185 + x1186)
    x1428 = 2 * x819
    x1429 = x0 * x484
    x1430 = x1191 * x50 - x1202 * x17 + x1419 - x1429
    x1431 = x0 * (x1430 + x818)
    x1432 = 2 * x1431
    x1433 = x1192 - x1203 - x1206 + x1207
    x1434 = x1433 * x17
    x1435 = x1423 * x50
    x1436 = x0 * (x1417 + x1422 - x1428 - x1432 - 2 * x1434 + 2 * x1435 + x87)
    x1437 = (
        x0 * (-x1417 - x1422 - 2 * x1424 + 2 * x1426 + 2 * x1427 + 2 * x807 + x86)
        + x1211
        - x1436
    )
    x1438 = -x0 * x498 + x1202 * x50 - x1212 * x17 + x1429
    x1439 = x0 * (x1438 + x826)
    x1440 = x17 * (x1203 - x1213 - x1214 + x1215)
    x1441 = x1433 * x50
    x1442 = -x17 * (
        -x0 * (x1428 + x1432 - 2 * x1439 - 2 * x1440 + 2 * x1441 - 2 * x827 + x98)
        + x1217
        + x1436
    )
    x1443 = x1437 * x50 + x1442
    x1444 = 2 * x841
    x1445 = x0 * (-2 * x203 + 2 * x337)
    x1446 = x0 * (-2 * x195 + 2 * x332) + x1222 * x50 - x1227 * x17 - x1445
    x1447 = 2 * x844
    x1448 = x0 * (-2 * x218 + 2 * x349)
    x1449 = x1227 * x50 - x1233 * x17 + x1445 - x1448
    x1450 = x0 * (x1444 - x1447 + x1449 + x240)
    x1451 = x0 * (-x1444 + x1446 + x239 + 2 * x840) + x1237 - x1450
    x1452 = -x0 * (-2 * x232 + 2 * x360) + x1233 * x50 - x1238 * x17 + x1448
    x1453 = -x17 * (-x0 * (x1447 + x1452 + x245 - 2 * x848) + x1240 + x1450)
    x1454 = x1451 * x50 + x1453
    x1455 = 2 * x861
    x1456 = x0 * x1105
    x1457 = x0 * x1108 + x1245 * x50 - x1250 * x17 - x1456
    x1458 = 2 * x864
    x1459 = x0 * x1113
    x1460 = x1250 * x50 - x1256 * x17 + x1456 - x1459
    x1461 = x0 * (x1455 - x1458 + x1460 + x312)
    x1462 = x0 * (-x1455 + x1457 + x311 + 2 * x860) + x1260 - x1461
    x1463 = -x0 * x1120 + x1256 * x50 - x1261 * x17 + x1459
    x1464 = -x17 * (-x0 * (x1458 + x1463 + x317 - 2 * x868) + x1263 + x1461)
    x1465 = x1462 * x50 + x1464
    x1466 = 2 * x881
    x1467 = x0 * x810
    x1468 = x0 * x805 + x1268 * x50 - x1273 * x17 - x1467
    x1469 = 2 * x884
    x1470 = x0 * x817
    x1471 = x1273 * x50 - x1279 * x17 + x1467 - x1470
    x1472 = x0 * (x1466 - x1469 + x1471 + x87)
    x1473 = x0 * (-x1466 + x1468 + x86 + 2 * x880) + x1283 - x1472
    x1474 = -x0 * x825 + x1279 * x50 - x1284 * x17 + x1470
    x1475 = -x17 * (-x0 * (x1469 + x1474 - 2 * x888 + x98) + x1286 + x1472)
    x1476 = x1473 * x50 + x1475
    x1477 = x1287 * x6 + x1288
    x1478 = 3 * x1117
    x1479 = 6 * x0
    x1480 = x1479 * x22
    x1481 = x13 * x1479
    x1482 = x1481 * x23
    x1483 = x0 * (-3 * x123 + x1480 - x1482 + 3 * x192)
    x1484 = x1481 * x29
    x1485 = x0 * (-3 * x131 + x1482 - x1484 + 3 * x200)
    x1486 = 3 * x396
    x1487 = 3 * x391
    x1488 = x185 * (x1483 - x1485 - x1486 + x1487 - 3 * x273 + x354 + 3 * x397)
    x1489 = -x1360 + x1363 + x414
    x1490 = 3 * x0
    x1491 = x1481 * x39
    x1492 = x0 * (-3 * x146 + x1484 - x1491 + 3 * x215)
    x1493 = 3 * x408
    x1494 = x185 * (x1485 + x1486 - x1492 - x1493 - 3 * x289 + x365 + 3 * x409)
    x1495 = x1360 - x1374 + x425
    x1496 = x1490 * (x1116 + x115 * x1489 + x1488 - x1494 - x1495 * x18)
    x1497 = 3 * x50
    x1498 = x1365 - x1376 - x1379 + x1381
    x1499 = 3 * x17
    x1500 = 3 * x1152
    x1501 = 3 * x256
    x1502 = 3 * x261
    x1503 = x0 * (3 * x1131 - 3 * x1132 + x1501 - x1502)
    x1504 = 3 * x1134
    x1505 = 3 * x270
    x1506 = x0 * (3 * x1144 - 3 * x1145 + x1502 - x1505)
    x1507 = 3 * x1147
    x1508 = x185 * (
        -3 * x1150 + 3 * x1151 + x1393 * x50 - x1403 * x17 + x1503 + x1504 - x1506 - x1507
    )
    x1509 = x1394 - x1404 - x1407 + x1408
    x1510 = 3 * x1197
    x1511 = 3 * x1172
    x1512 = 3 * x1173
    x1513 = x0 * (3 * x1174 - 3 * x1175 + x1511 - x1512)
    x1514 = 3 * x1177
    x1515 = 3 * x1188
    x1516 = x0 * (3 * x1189 - 3 * x1190 + x1512 - x1515)
    x1517 = 3 * x1192
    x1518 = x185 * (
        -3 * x1195 + 3 * x1196 + x1420 * x50 - x1430 * x17 + x1513 + x1514 - x1516 - x1517
    )
    x1519 = x1421 - x1431 - x1434 + x1435
    x1520 = 3 * x1219
    x1521 = 3 * x1224
    x1522 = x0 * (3 * x1225 - 3 * x1226 + x1520 - x1521)
    x1523 = 3 * x1228
    x1524 = 3 * x1242
    x1525 = 3 * x1247
    x1526 = x0 * (3 * x1248 - 3 * x1249 + x1524 - x1525)
    x1527 = 3 * x1251
    x1528 = 3 * x1265
    x1529 = 3 * x1270
    x1530 = x0 * (3 * x1271 - 3 * x1272 + x1528 - x1529)
    x1531 = 3 * x1274

    # 150 item(s)
    S = numpy.array(
        [
            x114,
            x0 * (3 * x129 - x138 + 3 * x143 - 3 * x155)
            - x0 * (x138 - 3 * x156 + 3 * x157 - 3 * x169)
            - x17 * x184
            + x175 * x50,
            x114,
            x0 * (3 * x198 - x207 + 3 * x212 - 3 * x224)
            - x0 * (x207 - 3 * x225 + 3 * x226 - 3 * x238)
            - x17 * x252
            + x244 * x50,
            x0 * (3 * x267 - x277 + 3 * x282 - 3 * x295)
            - x0 * (x277 - 3 * x296 + 3 * x297 - 3 * x310)
            - x17 * x326
            + x316 * x50,
            x114,
            x389,
            x455,
            x533,
            x114,
            x564,
            x0 * (x175 + x572)
            - x0 * (x184 + x578)
            - x17 * (x580 - x586 + x587 - x592)
            + x50 * (x579 - x580 + x582 - x585),
            x564,
            x0 * (x244 + x600)
            - x0 * (x252 + x606)
            - x17 * (x608 - x614 + x615 - x620)
            + x50 * (x607 - x608 + x610 - x613),
            x0 * (x316 + x628)
            - x0 * (x326 + x634)
            - x17 * (x636 - x642 + x643 - x648)
            + x50 * (x635 - x636 + x638 - x641),
            x564,
            x0 * (x376 + x656)
            - x0 * (x386 + x662)
            - x17 * (x664 - x672 + x673 - x678)
            + x50 * (x663 - x664 + x667 - x671),
            x0 * (x438 + x686)
            - x0 * (x452 + x692)
            - x17 * (x694 - x702 + x703 - x708)
            + x50 * (x693 - x694 + x697 - x701),
            x0 * (x517 + x716)
            - x0 * (x530 + x722)
            - x17 * (x724 - x732 + x733 - x738)
            + x50 * (x723 - x724 + x727 - x731),
            x564,
            x0 * (2 * x744 - 2 * x748)
            - x0 * (2 * x749 - 2 * x753)
            - x17 * x762
            + x50 * x759,
            x0 * (2 * x775 - 2 * x782)
            - x0 * (2 * x783 - 2 * x790)
            - x17 * x801
            + x50 * x796,
            x0 * (2 * x814 - 2 * x821)
            - x0 * (2 * x822 - 2 * x829)
            - x17 * x839
            + x50 * x835,
            x0 * (2 * x843 - 2 * x846)
            - x0 * (2 * x847 - 2 * x850)
            - x17 * x859
            + x50 * x856,
            x0 * (2 * x863 - 2 * x866)
            - x0 * (2 * x867 - 2 * x870)
            - x17 * x879
            + x50 * x876,
            x0 * (2 * x883 - 2 * x886)
            - x0 * (2 * x887 - 2 * x890)
            - x17 * x899
            + x50 * x896,
            x0 * (2 * x346 - 2 * x355)
            - x0 * (2 * x357 - 2 * x366)
            - x17 * x901
            + x50 * x900,
            x0 * (2 * x406 - 2 * x415)
            - x0 * (2 * x417 - 2 * x426)
            - x17 * x903
            + x50 * x902,
            x0 * (2 * x481 - 2 * x493)
            - x0 * (2 * x495 - 2 * x507)
            - x17 * x905
            + x50 * x904,
            x0 * (2 * x55 - 2 * x69) - x0 * (2 * x71 - 2 * x85) - x17 * x907 + x50 * x906,
            x920,
            x0 * (2 * x579 + 2 * x582 - 2 * x585 + x921 - x932)
            - x0 * (-2 * x586 + 2 * x587 - 2 * x592 + x923 + x932)
            - x17 * (x927 + x928 + x930 - x931)
            + x50 * (x922 + x925 + x926 - x927),
            x920,
            x0 * (2 * x607 + 2 * x610 - 2 * x613 + x933 - x944)
            - x0 * (-2 * x614 + 2 * x615 - 2 * x620 + x935 + x944)
            - x17 * (x939 + x940 + x942 - x943)
            + x50 * (x934 + x937 + x938 - x939),
            x0 * (2 * x635 + 2 * x638 - 2 * x641 + x945 - x956)
            - x0 * (-2 * x642 + 2 * x643 - 2 * x648 + x947 + x956)
            - x17 * (x951 + x952 + x954 - x955)
            + x50 * (x946 + x949 + x950 - x951),
            x920,
            x0 * (2 * x663 + 2 * x667 - 2 * x671 + x957 - x968)
            - x0 * (-2 * x672 + 2 * x673 - 2 * x678 + x959 + x968)
            - x17 * (x963 + x964 + x966 - x967)
            + x50 * (x958 + x961 + x962 - x963),
            x0 * (2 * x693 + 2 * x697 - 2 * x701 + x969 - x980)
            - x0 * (-2 * x702 + 2 * x703 - 2 * x708 + x971 + x980)
            - x17 * (x975 + x976 + x978 - x979)
            + x50 * (x970 + x973 + x974 - x975),
            x0 * (2 * x723 + 2 * x727 - 2 * x731 + x981 - x992)
            - x0 * (-2 * x732 + 2 * x733 - 2 * x738 + x983 + x992)
            - x17 * (x987 + x988 + x990 - x991)
            + x50 * (x982 + x985 + x986 - x987),
            x920,
            x0 * (x759 + x993)
            - x0 * (x762 + x994)
            - x17 * (-x1000 + x1001 + x1003 + x996)
            + x50 * (x995 - x996 + x997 + x999),
            x0 * (x1004 + x796)
            - x0 * (x1005 + x801)
            - x17 * (x1007 - x1011 + x1012 + x1014)
            + x50 * (x1006 - x1007 + x1008 + x1010),
            x0 * (x1015 + x835)
            - x0 * (x1016 + x839)
            - x17 * (x1018 - x1022 + x1023 + x1025)
            + x50 * (x1017 - x1018 + x1019 + x1021),
            x0 * (x1026 + x856)
            - x0 * (x1027 + x859)
            - x17 * (x1029 - x1033 + x1034 + x1036)
            + x50 * (x1028 - x1029 + x1030 + x1032),
            x0 * (x1037 + x876)
            - x0 * (x1038 + x879)
            - x17 * (x1040 - x1044 + x1045 + x1047)
            + x50 * (x1039 - x1040 + x1041 + x1043),
            x0 * (x1048 + x896)
            - x0 * (x1049 + x899)
            - x17 * (x1051 - x1055 + x1056 + x1058)
            + x50 * (x1050 - x1051 + x1052 + x1054),
            x0 * (x1059 + x900)
            - x0 * (x1060 + x901)
            - x17 * (x1062 - x1066 + x1067 + x1069)
            + x50 * (x1061 - x1062 + x1063 + x1065),
            x0 * (x1070 + x902)
            - x0 * (x1071 + x903)
            - x17 * (x1073 - x1077 + x1078 + x1080)
            + x50 * (x1072 - x1073 + x1074 + x1076),
            x0 * (x1081 + x904)
            - x0 * (x1082 + x905)
            - x17 * (x1084 - x1088 + x1089 + x1091)
            + x50 * (x1083 - x1084 + x1085 + x1087),
            x0 * (x1092 + x906)
            - x0 * (x1093 + x907)
            - x17 * (x1095 - x1099 + x1100 + x1102)
            + x50 * (x1094 - x1095 + x1096 + x1098),
            x0 * x1118 - x0 * x1125 + x1127 * x50 - x1130 * x17,
            x0 * x1153 - x0 * x1163 + x1165 * x50 - x1171 * x17,
            x0 * x1198 - x0 * x1209 + x1211 * x50 - x1217 * x17,
            x0 * x1229 - x0 * x1235 + x1237 * x50 - x1240 * x17,
            x0 * x1252 - x0 * x1258 + x1260 * x50 - x1263 * x17,
            x0 * x1275 - x0 * x1281 + x1283 * x50 - x1286 * x17,
            x0 * x370 - x0 * x373 - x17 * x385 + x375 * x50,
            x0 * x432 - x0 * x435 - x17 * x451 + x437 * x50,
            x0 * x511 - x0 * x514 - x17 * x529 + x50 * x516,
            x0 * x91 - x0 * x94 - x110 * x17 + x50 * x96,
            x1290,
            x0 * (-x1293 + 3 * x922 - 3 * x924 + 3 * x926)
            - x0 * (x1293 + 3 * x928 - 3 * x929 - 3 * x931)
            + x1291 * x50
            + x1292,
            x1290,
            x0 * (-x1296 + 3 * x934 - 3 * x936 + 3 * x938)
            - x0 * (x1296 + 3 * x940 - 3 * x941 - 3 * x943)
            + x1294 * x50
            + x1295,
            x0 * (-x1299 + 3 * x946 - 3 * x948 + 3 * x950)
            - x0 * (x1299 + 3 * x952 - 3 * x953 - 3 * x955)
            + x1297 * x50
            + x1298,
            x1290,
            x0 * (-x1302 + 3 * x958 - 3 * x960 + 3 * x962)
            - x0 * (x1302 + 3 * x964 - 3 * x965 - 3 * x967)
            + x1300 * x50
            + x1301,
            x0 * (-x1305 + 3 * x970 - 3 * x972 + 3 * x974)
            - x0 * (x1305 + 3 * x976 - 3 * x977 - 3 * x979)
            + x1303 * x50
            + x1304,
            x0 * (-x1308 + 3 * x982 - 3 * x984 + 3 * x986)
            - x0 * (x1308 + 3 * x988 - 3 * x989 - 3 * x991)
            + x1306 * x50
            + x1307,
            x1290,
            -x0 * (-2 * x1000 + 2 * x1001 - 2 * x1002 + x1311)
            + x0 * (-x1311 + 2 * x995 + 2 * x997 - 2 * x998)
            + x1309 * x50
            + x1310,
            x0 * (2 * x1006 + 2 * x1008 - 2 * x1009 - x1314)
            - x0 * (-2 * x1011 + 2 * x1012 - 2 * x1013 + x1314)
            + x1312 * x50
            + x1313,
            x0 * (2 * x1017 + 2 * x1019 - 2 * x1020 - x1317)
            - x0 * (-2 * x1022 + 2 * x1023 - 2 * x1024 + x1317)
            + x1315 * x50
            + x1316,
            x0 * (2 * x1028 + 2 * x1030 - 2 * x1031 - x1320)
            - x0 * (-2 * x1033 + 2 * x1034 - 2 * x1035 + x1320)
            + x1318 * x50
            + x1319,
            x0 * (2 * x1039 + 2 * x1041 - 2 * x1042 - x1323)
            - x0 * (-2 * x1044 + 2 * x1045 - 2 * x1046 + x1323)
            + x1321 * x50
            + x1322,
            x0 * (2 * x1050 + 2 * x1052 - 2 * x1053 - x1326)
            - x0 * (-2 * x1055 + 2 * x1056 - 2 * x1057 + x1326)
            + x1324 * x50
            + x1325,
            x0 * (2 * x1061 + 2 * x1063 - 2 * x1064 - x1329)
            - x0 * (-2 * x1066 + 2 * x1067 - 2 * x1068 + x1329)
            + x1327 * x50
            + x1328,
            x0 * (2 * x1072 + 2 * x1074 - 2 * x1075 - x1332)
            - x0 * (-2 * x1077 + 2 * x1078 - 2 * x1079 + x1332)
            + x1330 * x50
            + x1331,
            x0 * (2 * x1083 + 2 * x1085 - 2 * x1086 - x1335)
            - x0 * (-2 * x1088 + 2 * x1089 - 2 * x1090 + x1335)
            + x1333 * x50
            + x1334,
            x0 * (2 * x1094 + 2 * x1096 - 2 * x1097 - x1338)
            - x0 * (-2 * x1099 + 2 * x1100 - 2 * x1101 + x1338)
            + x1336 * x50
            + x1337,
            x0 * x1127 - x0 * x1130 + x1339 * x50 + x1340,
            x0 * x1165 - x0 * x1171 + x1341 * x50 + x1342,
            x0 * x1211 - x0 * x1217 + x1343 * x50 + x1344,
            x0 * x1237 - x0 * x1240 + x1345 * x50 + x1346,
            x0 * x1260 - x0 * x1263 + x1347 * x50 + x1348,
            x0 * x1283 - x0 * x1286 + x1349 * x50 + x1350,
            x0 * x375 - x0 * x385 + x1351 * x50 + x1352,
            x0 * x437 - x0 * x451 + x1353 * x50 + x1354,
            x0 * x516 - x0 * x529 + x1355 * x50 + x1356,
            -x0 * x110 + x0 * x96 + x1357 * x50 + x1358,
            x1389,
            x1416,
            x1443,
            x1454,
            x1465,
            x1476,
            x388,
            x454,
            x532,
            x113,
            x1477,
            x1291 * x6 + x1292,
            x1477,
            x1294 * x6 + x1295,
            x1297 * x6 + x1298,
            x1477,
            x1300 * x6 + x1301,
            x1303 * x6 + x1304,
            x1306 * x6 + x1307,
            x1477,
            x1309 * x6 + x1310,
            x1312 * x6 + x1313,
            x1315 * x6 + x1316,
            x1318 * x6 + x1319,
            x1321 * x6 + x1322,
            x1324 * x6 + x1325,
            x1327 * x6 + x1328,
            x1330 * x6 + x1331,
            x1333 * x6 + x1334,
            x1336 * x6 + x1337,
            x1339 * x6 + x1340,
            x1341 * x6 + x1342,
            x1343 * x6 + x1344,
            x1345 * x6 + x1346,
            x1347 * x6 + x1348,
            x1349 * x6 + x1350,
            x1351 * x6 + x1352,
            x1353 * x6 + x1354,
            x1355 * x6 + x1356,
            x1357 * x6 + x1358,
            x1383 * x6 + x1388,
            x1410 * x6 + x1415,
            x1437 * x6 + x1442,
            x1451 * x6 + x1453,
            x1462 * x6 + x1464,
            x1473 * x6 + x1475,
            x376 * x6 + x387,
            x438 * x6 + x453,
            x517 * x6 + x531,
            x112 + x6 * x97,
            x0
            * (
                3 * x1111
                - x1478
                + x1490
                * (
                    x1110
                    + x115 * (-x1363 + x1370 + x405)
                    - x1488
                    - x1489 * x18
                    + x185
                    * (
                        x0 * (-3 * x119 + x12 * x1481 - x1480 + 3 * x188)
                        - x1483
                        - x1487
                        - 3 * x264
                        + x345
                        + 3 * x390
                        + 3 * x392
                    )
                )
                - x1496
                + x1497 * (-x1365 - x1368 + x1371 + x1372)
                - x1498 * x1499
                + 3 * x744
                - 3 * x748
            )
            - x0
            * (
                -3 * x1124
                + x1478
                - x1490
                * (
                    x1123
                    + x115 * x1495
                    + x1494
                    - x18 * (x1374 - x1384 + x446)
                    - x185
                    * (
                        -x0 * (-x1481 * x57 + x1491 - 3 * x160 + 3 * x229)
                        + x1492
                        + x1493
                        - 3 * x304
                        + x380
                        - 3 * x419
                        + 3 * x420
                    )
                )
                + x1496
                + x1497 * x1498
                - x1499 * (x1376 - x1385 - x1386 + x1387)
                + 3 * x749
                - 3 * x753
            )
            + x1389,
            x0
            * (
                3 * x1143
                + x115 * (-x1394 - x1397 + x1399 + x1400)
                - x1500
                - x1508
                - x1509 * x18
                + x185
                * (
                    x0 * (3 * x1138 - 3 * x1139 - x1501 + 3 * x255)
                    - 3 * x1137
                    + 3 * x1141
                    + 3 * x1142
                    - x1393 * x17
                    + x1398 * x50
                    - x1503
                    - x1504
                )
                + 3 * x775
                - 3 * x782
            )
            - x0
            * (
                x115 * x1509
                - 3 * x1162
                + x1500
                + x1508
                - x18 * (x1404 - x1412 - x1413 + x1414)
                - x185
                * (
                    -x0 * (3 * x1154 - 3 * x1155 + x1505 - 3 * x286)
                    - 3 * x1157
                    - 3 * x1160
                    + 3 * x1161
                    + x1403 * x50
                    - x1411 * x17
                    + x1506
                    + x1507
                )
                + 3 * x783
                - 3 * x790
            )
            + x1416,
            -x0
            * (
                x115 * x1519
                - 3 * x1208
                + x1510
                + x1518
                - x18 * (x1431 - x1439 - x1440 + x1441)
                - x185
                * (
                    -x0 * (-3 * x1199 + 3 * x1200 - 3 * x1201 + x1515)
                    - 3 * x1203
                    - 3 * x1206
                    + 3 * x1207
                    + x1430 * x50
                    - x1438 * x17
                    + x1516
                    + x1517
                )
                + 3 * x822
                - 3 * x829
            )
            + x0
            * (
                x115 * (-x1421 - x1424 + x1426 + x1427)
                + 3 * x1187
                - x1510
                - x1518
                - x1519 * x18
                + x185
                * (
                    x0 * (3 * x1181 + 3 * x1182 - 3 * x1183 - x1511)
                    - 3 * x1180
                    + 3 * x1185
                    + 3 * x1186
                    - x1420 * x17
                    + x1425 * x50
                    - x1513
                    - x1514
                )
                + 3 * x814
                - 3 * x821
            )
            + x1443,
            x0
            * (
                x0 * (3 * x1218 + 3 * x1220 - 3 * x1221 - x1520)
                + 3 * x1223
                + x1446 * x50
                - x1449 * x17
                - x1522
                - x1523
                + 3 * x843
                - 3 * x846
            )
            - x0
            * (
                -x0 * (-3 * x1230 + 3 * x1231 - 3 * x1232 + x1521)
                - 3 * x1234
                + x1449 * x50
                - x1452 * x17
                + x1522
                + x1523
                + 3 * x847
                - 3 * x850
            )
            + x1454,
            x0
            * (
                x0 * (3 * x1241 + 3 * x1243 - 3 * x1244 - x1524)
                + 3 * x1246
                + x1457 * x50
                - x1460 * x17
                - x1526
                - x1527
                + 3 * x863
                - 3 * x866
            )
            - x0
            * (
                -x0 * (-3 * x1253 + 3 * x1254 - 3 * x1255 + x1525)
                - 3 * x1257
                + x1460 * x50
                - x1463 * x17
                + x1526
                + x1527
                + 3 * x867
                - 3 * x870
            )
            + x1465,
            x0
            * (
                x0 * (3 * x1264 + 3 * x1266 - 3 * x1267 - x1528)
                + 3 * x1269
                + x1468 * x50
                - x1471 * x17
                - x1530
                - x1531
                + 3 * x883
                - 3 * x886
            )
            - x0
            * (
                -x0 * (-3 * x1276 + 3 * x1277 - 3 * x1278 + x1529)
                - 3 * x1280
                + x1471 * x50
                - x1474 * x17
                + x1530
                + x1531
                + 3 * x887
                - 3 * x890
            )
            + x1476,
            x389,
            x455,
            x533,
            x114,
        ]
    )
    return S


def coulomb_44(a, A, b, B, C):
    """Cartesian (gg) 1-electron Coulomb integral.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = a + b
    x2 = x1 ** (-1.0)
    x3 = -x2 * (a * A[0] + b * B[0])
    x4 = -x2 * (a * A[1] + b * B[1])
    x5 = -x2 * (a * A[2] + b * B[2])
    x6 = numpy.sqrt(abs(x3 + B[0]) ** 2 + abs(x4 + B[1]) ** 2 + abs(x5 + B[2]) ** 2)
    x7 = x3 + C[0]
    x8 = x4 + C[1]
    x9 = x5 + C[2]
    x10 = x7 ** 2 + x8 ** 2 + x9 ** 2
    x11 = x1 * x10
    x12 = boys(0, x11)
    x13 = numpy.pi * x2 * numpy.exp(-a * b * x10 * x2)
    x14 = 2 * x6
    x15 = x13 * x14
    x16 = boys(1, x11)
    x17 = numpy.sqrt(abs(x7) ** 2 + abs(x8) ** 2 + abs(x9) ** 2)
    x18 = 2 * x17
    x19 = x13 * x18
    x20 = -x16 * x19
    x21 = x12 * x15 + x20
    x22 = x13 * x16
    x23 = boys(2, x11)
    x24 = -x19 * x23
    x25 = x14 * x22 + x24
    x26 = x17 * x25
    x27 = -x26
    x28 = x21 * x6 + x27
    x29 = boys(3, x11)
    x30 = -x19 * x29
    x31 = x15 * x23 + x30
    x32 = x17 * x31
    x33 = -x32
    x34 = x25 * x6 + x33
    x35 = x17 * x34
    x36 = -x35
    x37 = x28 * x6 + x36
    x38 = boys(4, x11)
    x39 = -x19 * x38
    x40 = x15 * x29 + x39
    x41 = x17 * x40
    x42 = -x41
    x43 = x31 * x6 + x42
    x44 = x17 * x43
    x45 = -x44
    x46 = x34 * x6 + x45
    x47 = x17 * x46
    x48 = -x47
    x49 = x37 * x6 + x48
    x50 = x0 * x49
    x51 = boys(5, x11)
    x52 = -x19 * x51
    x53 = x15 * x38 + x52
    x54 = x17 * x53
    x55 = -x54
    x56 = x40 * x6 + x55
    x57 = x17 * x56
    x58 = -x57
    x59 = x43 * x6 + x58
    x60 = x17 * x59
    x61 = -x60
    x62 = x46 * x6 + x61
    x63 = x0 * x62
    x64 = 3 * x63
    x65 = numpy.sqrt(abs(x3 + A[0]) ** 2 + abs(x4 + A[1]) ** 2 + abs(x5 + A[2]) ** 2)
    x66 = x49 * x65
    x67 = x17 * x62
    x68 = -x67
    x69 = x66 + x68
    x70 = x65 * x69
    x71 = x62 * x65
    x72 = boys(6, x11)
    x73 = -x19 * x72
    x74 = x15 * x51 + x73
    x75 = x17 * x74
    x76 = -x75
    x77 = x53 * x6 + x76
    x78 = x17 * x77
    x79 = -x78
    x80 = x56 * x6 + x79
    x81 = x17 * x80
    x82 = -x81
    x83 = x59 * x6 + x82
    x84 = x17 * x83
    x85 = -x84
    x86 = x71 + x85
    x87 = x17 * x86
    x88 = x0 * x83
    x89 = x65 * x86
    x90 = x65 * x83
    x91 = boys(7, x11)
    x92 = -x19 * x91
    x93 = x15 * x72 + x92
    x94 = x17 * x93
    x95 = -x94
    x96 = x6 * x74 + x95
    x97 = x17 * x96
    x98 = -x97
    x99 = x6 * x77 + x98
    x100 = x17 * x99
    x101 = -x100
    x102 = x101 + x6 * x80
    x103 = x102 * x17
    x104 = -x103
    x105 = x104 + x90
    x106 = x105 * x17
    x107 = 2 * x66 - 2 * x67
    x108 = 2 * x71 - 2 * x84
    x109 = x0 * x108
    x110 = -x87
    x111 = x110 + x70
    x112 = x111 + x50 - x63
    x113 = -x106
    x114 = x113 + x89
    x115 = x114 + x63 - x88
    x116 = -x115 * x17
    x117 = x112 * x65 + x116
    x118 = x0 * x107 - x109 + x117
    x119 = -2 * x103 + 2 * x90
    x120 = -x19 * boys(8, x11)
    x121 = x17 * (x120 + x15 * x91)
    x122 = -x121
    x123 = x17 * (x122 + x6 * x93)
    x124 = -x123
    x125 = x17 * (x124 + x6 * x96)
    x126 = -x125
    x127 = -x17 * (x126 + x6 * x99)
    x128 = x102 * x65 + x127
    x129 = -x128 * x17
    x130 = x105 * x65 + x129
    x131 = -x0 * x102 + x130 + x88
    x132 = -x131 * x17
    x133 = x115 * x65 + x132
    x134 = -x0 * x119 + x109 + x133
    x135 = -x134 * x17
    x136 = x118 * x65 + x135
    x137 = (
        -x0 * (-3 * x106 + x64 - 3 * x88 + 3 * x89)
        + x0 * (3 * x50 - x64 + 3 * x70 - 3 * x87)
        + x136
    )
    x138 = 2 * x65
    x139 = x13 * x138
    x140 = x12 * x139 + x20
    x141 = x138 * x22 + x24
    x142 = x141 * x17
    x143 = -x142
    x144 = x140 * x6 + x143
    x145 = x139 * x23 + x30
    x146 = x145 * x17
    x147 = -x146
    x148 = x141 * x6 + x147
    x149 = x148 * x17
    x150 = -x149
    x151 = x144 * x6 + x150
    x152 = x139 * x29 + x39
    x153 = x152 * x17
    x154 = -x153
    x155 = x145 * x6 + x154
    x156 = x155 * x17
    x157 = -x156
    x158 = x148 * x6 + x157
    x159 = x158 * x17
    x160 = -x159
    x161 = x151 * x6 + x160
    x162 = x0 * x161
    x163 = x139 * x38 + x52
    x164 = x163 * x17
    x165 = -x164
    x166 = x152 * x6 + x165
    x167 = x166 * x17
    x168 = -x167
    x169 = x155 * x6 + x168
    x170 = x169 * x17
    x171 = -x170
    x172 = x158 * x6 + x171
    x173 = x0 * x172
    x174 = 3 * x173
    x175 = x161 * x65
    x176 = x17 * x172
    x177 = -x176
    x178 = x175 + x177
    x179 = x178 * x65
    x180 = x172 * x65
    x181 = x139 * x51 + x73
    x182 = x17 * x181
    x183 = -x182
    x184 = x163 * x6 + x183
    x185 = x17 * x184
    x186 = -x185
    x187 = x166 * x6 + x186
    x188 = x17 * x187
    x189 = -x188
    x190 = x169 * x6 + x189
    x191 = x17 * x190
    x192 = -x191
    x193 = x180 + x192
    x194 = x17 * x193
    x195 = x0 * x190
    x196 = x193 * x65
    x197 = x190 * x65
    x198 = x139 * x72 + x92
    x199 = x17 * x198
    x200 = -x199
    x201 = x181 * x6 + x200
    x202 = x17 * x201
    x203 = -x202
    x204 = x184 * x6 + x203
    x205 = x17 * x204
    x206 = -x205
    x207 = x187 * x6 + x206
    x208 = x17 * x207
    x209 = -x208
    x210 = x197 + x209
    x211 = x17 * x210
    x212 = 2 * x175 - 2 * x176
    x213 = 2 * x180 - 2 * x191
    x214 = x0 * x213
    x215 = x162 - x173 + x179 - x194
    x216 = x173 - x195 + x196 - x211
    x217 = x0 * x212 - x17 * x216 - x214 + x215 * x65
    x218 = 2 * x197 - 2 * x208
    x219 = x120 + x139 * x91
    x220 = -x17 * x219
    x221 = x17 * (x198 * x6 + x220)
    x222 = -x221
    x223 = x17 * (x201 * x6 + x222)
    x224 = -x223
    x225 = -x17 * (x204 * x6 + x224)
    x226 = x207 * x65 + x225
    x227 = -x0 * x207 - x17 * x226 + x195 + x210 * x65
    x228 = -x0 * x218 - x17 * x227 + x214 + x216 * x65
    x229 = 2 * x0
    x230 = x22 * x229
    x231 = x13 * x229
    x232 = x140 * x65
    x233 = x143 + x232
    x234 = x12 * x231 - x230 + x233
    x235 = x23 * x231
    x236 = x141 * x65
    x237 = x147 + x236
    x238 = x230 - x235 + x237
    x239 = x17 * x238
    x240 = -x239
    x241 = x234 * x6 + x240
    x242 = x231 * x29
    x243 = x145 * x65
    x244 = x154 + x243
    x245 = x235 - x242 + x244
    x246 = x17 * x245
    x247 = -x246
    x248 = x238 * x6 + x247
    x249 = x17 * x248
    x250 = -x249
    x251 = x241 * x6 + x250
    x252 = x0 * x251
    x253 = x231 * x38
    x254 = x152 * x65
    x255 = x165 + x254
    x256 = x242 - x253 + x255
    x257 = x17 * x256
    x258 = -x257
    x259 = x245 * x6 + x258
    x260 = x17 * x259
    x261 = -x260
    x262 = x248 * x6 + x261
    x263 = x0 * x262
    x264 = 3 * x263
    x265 = x251 * x65
    x266 = x17 * x262
    x267 = -x266
    x268 = x265 + x267
    x269 = x268 * x65
    x270 = x262 * x65
    x271 = x231 * x51
    x272 = x163 * x65
    x273 = x183 + x272
    x274 = x253 - x271 + x273
    x275 = x17 * x274
    x276 = -x275
    x277 = x256 * x6 + x276
    x278 = x17 * x277
    x279 = -x278
    x280 = x259 * x6 + x279
    x281 = x17 * x280
    x282 = -x281
    x283 = x270 + x282
    x284 = x17 * x283
    x285 = x0 * x280
    x286 = x283 * x65
    x287 = x280 * x65
    x288 = x231 * x72
    x289 = x181 * x65
    x290 = x200 + x289
    x291 = x271 - x288 + x290
    x292 = x17 * x291
    x293 = -x292
    x294 = x274 * x6 + x293
    x295 = x17 * x294
    x296 = -x295
    x297 = x277 * x6 + x296
    x298 = x17 * x297
    x299 = -x298
    x300 = x287 + x299
    x301 = x17 * x300
    x302 = 2 * x265 - 2 * x266
    x303 = 2 * x270 - 2 * x281
    x304 = x0 * x303
    x305 = x252 - x263 + x269 - x284
    x306 = x263 - x285 + x286 - x301
    x307 = x0 * x302 - x17 * x306 - x304 + x305 * x65
    x308 = 2 * x287 - 2 * x298
    x309 = x198 * x65 + x220
    x310 = -x231 * x91 + x288 + x309
    x311 = -x17 * x310
    x312 = x17 * (x291 * x6 + x311)
    x313 = -x312
    x314 = -x17 * (x294 * x6 + x313)
    x315 = x297 * x65 + x314
    x316 = -x0 * x297 - x17 * x315 + x285 + x300 * x65
    x317 = -x0 * x308 - x17 * x316 + x304 + x306 * x65
    x318 = x21 * x65
    x319 = x27 + x318
    x320 = x0 * x140
    x321 = x0 * x141
    x322 = x320 - x321
    x323 = x319 + x322
    x324 = x25 * x65
    x325 = x324 + x33
    x326 = x0 * x145
    x327 = x321 - x326
    x328 = x325 + x327
    x329 = x17 * x328
    x330 = -x329
    x331 = x323 * x6 + x330
    x332 = x31 * x65
    x333 = x332 + x42
    x334 = x0 * x152
    x335 = x326 - x334
    x336 = x333 + x335
    x337 = x17 * x336
    x338 = -x337
    x339 = x328 * x6 + x338
    x340 = x17 * x339
    x341 = -x340
    x342 = x331 * x6 + x341
    x343 = x0 * x342
    x344 = x40 * x65
    x345 = x344 + x55
    x346 = x0 * x163
    x347 = x334 - x346
    x348 = x345 + x347
    x349 = x17 * x348
    x350 = -x349
    x351 = x336 * x6 + x350
    x352 = x17 * x351
    x353 = -x352
    x354 = x339 * x6 + x353
    x355 = x0 * x354
    x356 = 3 * x355
    x357 = x342 * x65
    x358 = x17 * x354
    x359 = -x358
    x360 = x357 + x359
    x361 = x360 * x65
    x362 = x354 * x65
    x363 = x53 * x65
    x364 = x363 + x76
    x365 = x0 * x181
    x366 = x346 - x365
    x367 = x364 + x366
    x368 = x17 * x367
    x369 = -x368
    x370 = x348 * x6 + x369
    x371 = x17 * x370
    x372 = -x371
    x373 = x351 * x6 + x372
    x374 = x17 * x373
    x375 = -x374
    x376 = x362 + x375
    x377 = x17 * x376
    x378 = x0 * x373
    x379 = x376 * x65
    x380 = x373 * x65
    x381 = x65 * x74
    x382 = x381 + x95
    x383 = x0 * x198
    x384 = x365 - x383
    x385 = x382 + x384
    x386 = x17 * x385
    x387 = -x386
    x388 = x367 * x6 + x387
    x389 = x17 * x388
    x390 = -x389
    x391 = x370 * x6 + x390
    x392 = x17 * x391
    x393 = -x392
    x394 = x380 + x393
    x395 = x17 * x394
    x396 = 2 * x357 - 2 * x358
    x397 = 2 * x362 - 2 * x374
    x398 = x0 * x397
    x399 = x343 - x355 + x361 - x377
    x400 = x355 - x378 + x379 - x395
    x401 = x0 * x396 - x17 * x400 - x398 + x399 * x65
    x402 = 2 * x380 - 2 * x392
    x403 = x0 * x219
    x404 = x65 * x93
    x405 = x122 + x404
    x406 = x17 * (x383 - x403 + x405)
    x407 = -x406
    x408 = x17 * (x385 * x6 + x407)
    x409 = -x408
    x410 = -x17 * (x388 * x6 + x409)
    x411 = x391 * x65 + x410
    x412 = -x0 * x391 - x17 * x411 + x378 + x394 * x65
    x413 = -x0 * x402 - x17 * x412 + x398 + x400 * x65
    x414 = 4 * x65
    x415 = x13 * x414
    x416 = 4 * x17
    x417 = x13 * x416
    x418 = x0 * (x22 * x414 - x23 * x417)
    x419 = x234 * x65
    x420 = x240 + x419
    x421 = x0 * (x12 * x415 - x22 * x416) - x418 + x420
    x422 = x0 * (x23 * x415 - x29 * x417)
    x423 = x238 * x65
    x424 = x247 + x423
    x425 = x418 - x422 + x424
    x426 = x17 * x425
    x427 = -x426
    x428 = x421 * x6 + x427
    x429 = x0 * x428
    x430 = x0 * (x29 * x415 - x38 * x417)
    x431 = x245 * x65
    x432 = x258 + x431
    x433 = x422 - x430 + x432
    x434 = x17 * x433
    x435 = -x434
    x436 = x425 * x6 + x435
    x437 = x0 * x436
    x438 = 3 * x437
    x439 = x428 * x65
    x440 = x17 * x436
    x441 = -x440
    x442 = x439 + x441
    x443 = x442 * x65
    x444 = x436 * x65
    x445 = x0 * (x38 * x415 - x417 * x51)
    x446 = x256 * x65
    x447 = x276 + x446
    x448 = x430 - x445 + x447
    x449 = x17 * x448
    x450 = -x449
    x451 = x433 * x6 + x450
    x452 = x17 * x451
    x453 = -x452
    x454 = x444 + x453
    x455 = x17 * x454
    x456 = x0 * x451
    x457 = x454 * x65
    x458 = x451 * x65
    x459 = x0 * (x415 * x51 - x417 * x72)
    x460 = x274 * x65
    x461 = x293 + x460
    x462 = x445 - x459 + x461
    x463 = x17 * x462
    x464 = -x463
    x465 = x448 * x6 + x464
    x466 = x17 * x465
    x467 = -x466
    x468 = x458 + x467
    x469 = x17 * x468
    x470 = 2 * x439 - 2 * x440
    x471 = 2 * x444 - 2 * x452
    x472 = x0 * x471
    x473 = x429 - x437 + x443 - x455
    x474 = x437 - x456 + x457 - x469
    x475 = x0 * x470 - x17 * x474 - x472 + x473 * x65
    x476 = 2 * x458 - 2 * x466
    x477 = x291 * x65 + x311
    x478 = -x0 * (x415 * x72 - x417 * x91) + x459 + x477
    x479 = -x17 * x478
    x480 = -x17 * (x462 * x6 + x479)
    x481 = x465 * x65 + x480
    x482 = -x0 * x465 - x17 * x481 + x456 + x468 * x65
    x483 = -x0 * x476 - x17 * x482 + x472 + x474 * x65
    x484 = x0 * (x21 + x234)
    x485 = x0 * (x238 + x25)
    x486 = x323 * x65
    x487 = x330 + x486
    x488 = x484 - x485 + x487
    x489 = x0 * (x245 + x31)
    x490 = x328 * x65
    x491 = x338 + x490
    x492 = x485 - x489 + x491
    x493 = x17 * x492
    x494 = -x493
    x495 = x488 * x6 + x494
    x496 = x0 * x495
    x497 = x0 * (x256 + x40)
    x498 = x336 * x65
    x499 = x350 + x498
    x500 = x489 - x497 + x499
    x501 = x17 * x500
    x502 = -x501
    x503 = x492 * x6 + x502
    x504 = x0 * x503
    x505 = 3 * x504
    x506 = x495 * x65
    x507 = x17 * x503
    x508 = -x507
    x509 = x506 + x508
    x510 = x509 * x65
    x511 = x503 * x65
    x512 = x0 * (x274 + x53)
    x513 = x348 * x65
    x514 = x369 + x513
    x515 = x497 - x512 + x514
    x516 = x17 * x515
    x517 = -x516
    x518 = x500 * x6 + x517
    x519 = x17 * x518
    x520 = -x519
    x521 = x511 + x520
    x522 = x17 * x521
    x523 = x0 * x518
    x524 = x521 * x65
    x525 = x518 * x65
    x526 = x0 * (x291 + x74)
    x527 = x367 * x65
    x528 = x387 + x527
    x529 = x512 - x526 + x528
    x530 = x17 * x529
    x531 = -x530
    x532 = x515 * x6 + x531
    x533 = x17 * x532
    x534 = -x533
    x535 = x525 + x534
    x536 = x17 * x535
    x537 = 2 * x506 - 2 * x507
    x538 = 2 * x511 - 2 * x519
    x539 = x0 * x538
    x540 = x496 - x504 + x510 - x522
    x541 = x504 - x523 + x524 - x536
    x542 = x0 * x537 - x17 * x541 - x539 + x540 * x65
    x543 = 2 * x525 - 2 * x533
    x544 = x0 * (x310 + x93)
    x545 = x385 * x65
    x546 = x407 + x545
    x547 = x17 * (x526 - x544 + x546)
    x548 = -x547
    x549 = -x17 * (x529 * x6 + x548)
    x550 = x532 * x65 + x549
    x551 = -x0 * x532 - x17 * x550 + x523 + x535 * x65
    x552 = -x0 * x543 - x17 * x551 + x539 + x541 * x65
    x553 = 2 * x321
    x554 = -2 * x26 + 2 * x318
    x555 = 2 * x320 - x553 + x554
    x556 = x0 * x555
    x557 = 2 * x326
    x558 = -2 * x32 + 2 * x324
    x559 = x553 - x557 + x558
    x560 = x0 * x559
    x561 = x28 * x65
    x562 = x36 + x561
    x563 = x556 - x560 + x562
    x564 = 2 * x334
    x565 = 2 * x332 - 2 * x41
    x566 = x557 - x564 + x565
    x567 = x0 * x566
    x568 = x34 * x65
    x569 = x45 + x568
    x570 = x560 - x567 + x569
    x571 = x17 * x570
    x572 = -x571
    x573 = x563 * x6 + x572
    x574 = x0 * x573
    x575 = 2 * x346
    x576 = 2 * x344 - 2 * x54
    x577 = x564 - x575 + x576
    x578 = x0 * x577
    x579 = x43 * x65
    x580 = x579 + x58
    x581 = x567 - x578 + x580
    x582 = x17 * x581
    x583 = -x582
    x584 = x570 * x6 + x583
    x585 = x0 * x584
    x586 = 3 * x585
    x587 = x573 * x65
    x588 = x17 * x584
    x589 = -x588
    x590 = x587 + x589
    x591 = x590 * x65
    x592 = x584 * x65
    x593 = 2 * x365
    x594 = 2 * x363 - 2 * x75
    x595 = x575 - x593 + x594
    x596 = x0 * x595
    x597 = x56 * x65
    x598 = x597 + x79
    x599 = x578 - x596 + x598
    x600 = x17 * x599
    x601 = -x600
    x602 = x581 * x6 + x601
    x603 = x17 * x602
    x604 = -x603
    x605 = x592 + x604
    x606 = x17 * x605
    x607 = x0 * x602
    x608 = x605 * x65
    x609 = x602 * x65
    x610 = 2 * x383
    x611 = 2 * x381 - 2 * x94
    x612 = x593 - x610 + x611
    x613 = x0 * x612
    x614 = x65 * x77
    x615 = x614 + x98
    x616 = x596 - x613 + x615
    x617 = x17 * x616
    x618 = -x617
    x619 = x599 * x6 + x618
    x620 = x17 * x619
    x621 = -x620
    x622 = x609 + x621
    x623 = x17 * x622
    x624 = 2 * x587 - 2 * x588
    x625 = 2 * x592 - 2 * x603
    x626 = x0 * x625
    x627 = x574 - x585 + x591 - x606
    x628 = x585 - x607 + x608 - x623
    x629 = x0 * x624 - x17 * x628 - x626 + x627 * x65
    x630 = 2 * x609 - 2 * x620
    x631 = x0 * (-2 * x121 - 2 * x403 + 2 * x404 + x610)
    x632 = x65 * x96
    x633 = x124 + x632
    x634 = x17 * (x613 - x631 + x633)
    x635 = -x634
    x636 = -x17 * (x6 * x616 + x635)
    x637 = x619 * x65 + x636
    x638 = -x0 * x619 - x17 * x637 + x607 + x622 * x65
    x639 = -x0 * x630 - x17 * x638 + x626 + x628 * x65
    x640 = 6 * x0
    x641 = x13 * x640
    x642 = x22 * x640
    x643 = x23 * x641
    x644 = x0 * (-3 * x146 + 3 * x236 + x642 - x643)
    x645 = x421 * x65
    x646 = x427 + x645
    x647 = x0 * (x12 * x641 - 3 * x142 + 3 * x232 - x642) - x644 + x646
    x648 = x0 * x647
    x649 = x29 * x641
    x650 = x0 * (-3 * x153 + 3 * x243 + x643 - x649)
    x651 = x425 * x65
    x652 = x435 + x651
    x653 = x644 - x650 + x652
    x654 = x0 * x653
    x655 = 3 * x654
    x656 = x647 * x65
    x657 = x17 * x653
    x658 = -x657
    x659 = x656 + x658
    x660 = x65 * x659
    x661 = x65 * x653
    x662 = x38 * x641
    x663 = x0 * (-3 * x164 + 3 * x254 + x649 - x662)
    x664 = x433 * x65
    x665 = x450 + x664
    x666 = x650 - x663 + x665
    x667 = x17 * x666
    x668 = -x667
    x669 = x661 + x668
    x670 = x17 * x669
    x671 = x0 * x666
    x672 = x65 * x669
    x673 = x65 * x666
    x674 = x51 * x641
    x675 = x0 * (-3 * x182 + 3 * x272 + x662 - x674)
    x676 = x448 * x65
    x677 = x464 + x676
    x678 = x663 - x675 + x677
    x679 = x17 * x678
    x680 = -x679
    x681 = x673 + x680
    x682 = x17 * x681
    x683 = x0 * (2 * x661 - 2 * x667)
    x684 = -x670
    x685 = x660 + x684
    x686 = x648 - x654 + x685
    x687 = -x682
    x688 = x672 + x687
    x689 = x654 - x671 + x688
    x690 = -x17 * x689
    x691 = x65 * x686 + x690
    x692 = x0 * (2 * x656 - 2 * x657) - x683 + x691
    x693 = x462 * x65 + x479
    x694 = -x0 * (-3 * x199 + 3 * x289 - x641 * x72 + x674) + x675 + x693
    x695 = -x17 * x694
    x696 = x65 * x678 + x695
    x697 = -x17 * x696
    x698 = x65 * x681 + x697
    x699 = -x0 * x678 + x671 + x698
    x700 = -x17 * x699
    x701 = x65 * x689 + x700
    x702 = -x0 * (2 * x673 - 2 * x679) + x683 + x701
    x703 = -x17 * x702
    x704 = x65 * x692 + x703
    x705 = (
        x0 * (3 * x648 - x655 + 3 * x660 - 3 * x670)
        - x0 * (x655 - 3 * x671 + 3 * x672 - 3 * x682)
        + x704
    )
    x706 = x0 * (x421 + x555)
    x707 = x0 * (x425 + x559)
    x708 = x488 * x65
    x709 = x494 + x708
    x710 = x706 - x707 + x709
    x711 = x0 * x710
    x712 = x0 * (x433 + x566)
    x713 = x492 * x65
    x714 = x502 + x713
    x715 = x707 - x712 + x714
    x716 = x0 * x715
    x717 = 3 * x716
    x718 = x65 * x710
    x719 = x17 * x715
    x720 = -x719
    x721 = x718 + x720
    x722 = x65 * x721
    x723 = x65 * x715
    x724 = x0 * (x448 + x577)
    x725 = x500 * x65
    x726 = x517 + x725
    x727 = x712 - x724 + x726
    x728 = x17 * x727
    x729 = -x728
    x730 = x723 + x729
    x731 = x17 * x730
    x732 = x0 * x727
    x733 = x65 * x730
    x734 = x65 * x727
    x735 = x0 * (x462 + x595)
    x736 = x515 * x65
    x737 = x531 + x736
    x738 = x724 - x735 + x737
    x739 = x17 * x738
    x740 = -x739
    x741 = x734 + x740
    x742 = x17 * x741
    x743 = 2 * x718 - 2 * x719
    x744 = 2 * x723 - 2 * x728
    x745 = x0 * x744
    x746 = -x731
    x747 = x722 + x746
    x748 = x711 - x716 + x747
    x749 = -x742
    x750 = x733 + x749
    x751 = x716 - x732 + x750
    x752 = -x17 * x751
    x753 = x65 * x748 + x752
    x754 = x0 * x743 - x745 + x753
    x755 = 2 * x734 - 2 * x739
    x756 = x65 * x738
    x757 = x0 * (x478 + x612)
    x758 = x529 * x65
    x759 = x548 + x758
    x760 = x17 * (x735 - x757 + x759)
    x761 = -x760
    x762 = x756 + x761
    x763 = -x17 * x762
    x764 = x65 * x741 + x763
    x765 = -x0 * x738 + x732 + x764
    x766 = -x17 * x765
    x767 = x65 * x751 + x766
    x768 = -x0 * x755 + x745 + x767
    x769 = -x17 * x768
    x770 = x65 * x754 + x769
    x771 = (
        x0 * (3 * x711 - x717 + 3 * x722 - 3 * x731)
        - x0 * (x717 - 3 * x732 + 3 * x733 - 3 * x742)
        + x770
    )
    x772 = 2 * x485
    x773 = -2 * x329 + 2 * x486
    x774 = x0 * (x28 + 2 * x484 - x772 + x773)
    x775 = 2 * x489
    x776 = -2 * x337 + 2 * x490
    x777 = x0 * (x34 + x772 - x775 + x776)
    x778 = x563 * x65
    x779 = x572 + x778
    x780 = x774 - x777 + x779
    x781 = x0 * x780
    x782 = 2 * x497
    x783 = -2 * x349 + 2 * x498
    x784 = x0 * (x43 + x775 - x782 + x783)
    x785 = x570 * x65
    x786 = x583 + x785
    x787 = x777 - x784 + x786
    x788 = x0 * x787
    x789 = 3 * x788
    x790 = x65 * x780
    x791 = x17 * x787
    x792 = -x791
    x793 = x790 + x792
    x794 = x65 * x793
    x795 = x65 * x787
    x796 = 2 * x512
    x797 = -2 * x368 + 2 * x513
    x798 = x0 * (x56 + x782 - x796 + x797)
    x799 = x581 * x65
    x800 = x601 + x799
    x801 = x784 - x798 + x800
    x802 = x17 * x801
    x803 = -x802
    x804 = x795 + x803
    x805 = x17 * x804
    x806 = x0 * x801
    x807 = x65 * x804
    x808 = x65 * x801
    x809 = 2 * x526
    x810 = -2 * x386 + 2 * x527
    x811 = x0 * (x77 + x796 - x809 + x810)
    x812 = x599 * x65
    x813 = x618 + x812
    x814 = x798 - x811 + x813
    x815 = x17 * x814
    x816 = -x815
    x817 = x808 + x816
    x818 = x17 * x817
    x819 = x0 * (2 * x795 - 2 * x802)
    x820 = -x805
    x821 = x794 + x820
    x822 = x781 - x788 + x821
    x823 = -x818
    x824 = x807 + x823
    x825 = x788 - x806 + x824
    x826 = -x17 * x825
    x827 = x65 * x822 + x826
    x828 = x0 * (2 * x790 - 2 * x791) - x819 + x827
    x829 = x65 * x814
    x830 = x0 * (-2 * x406 - 2 * x544 + 2 * x545 + x809 + x96)
    x831 = x616 * x65
    x832 = x635 + x831
    x833 = x17 * (x811 - x830 + x832)
    x834 = -x833
    x835 = x829 + x834
    x836 = -x17 * x835
    x837 = x65 * x817 + x836
    x838 = -x0 * x814 + x806 + x837
    x839 = -x17 * x838
    x840 = x65 * x825 + x839
    x841 = -x0 * (2 * x808 - 2 * x815) + x819 + x840
    x842 = -x17 * x841
    x843 = x65 * x828 + x842
    x844 = (
        x0 * (3 * x781 - x789 + 3 * x794 - 3 * x805)
        - x0 * (x789 - 3 * x806 + 3 * x807 - 3 * x818)
        + x843
    )
    x845 = 3 * x560
    x846 = x0 * (-3 * x35 + 3 * x556 + 3 * x561 - x845)
    x847 = 3 * x567
    x848 = x0 * (-3 * x44 + 3 * x568 + x845 - x847)
    x849 = x37 * x65
    x850 = x48 + x849
    x851 = x846 - x848 + x850
    x852 = x0 * x851
    x853 = 3 * x578
    x854 = x0 * (-3 * x57 + 3 * x579 + x847 - x853)
    x855 = x46 * x65
    x856 = x61 + x855
    x857 = x848 - x854 + x856
    x858 = x0 * x857
    x859 = 3 * x858
    x860 = x65 * x851
    x861 = x17 * x857
    x862 = -x861
    x863 = x860 + x862
    x864 = x65 * x863
    x865 = x65 * x857
    x866 = 3 * x596
    x867 = x0 * (3 * x597 - 3 * x78 + x853 - x866)
    x868 = x59 * x65
    x869 = x82 + x868
    x870 = x854 - x867 + x869
    x871 = x17 * x870
    x872 = -x871
    x873 = x865 + x872
    x874 = x17 * x873
    x875 = x0 * x870
    x876 = x65 * x873
    x877 = x65 * x870
    x878 = 3 * x613
    x879 = x0 * (3 * x614 + x866 - x878 - 3 * x97)
    x880 = x65 * x80
    x881 = x101 + x880
    x882 = x867 - x879 + x881
    x883 = x17 * x882
    x884 = -x883
    x885 = x877 + x884
    x886 = x17 * x885
    x887 = x0 * (2 * x865 - 2 * x871)
    x888 = -x874
    x889 = x864 + x888
    x890 = x852 - x858 + x889
    x891 = -x886
    x892 = x876 + x891
    x893 = x858 - x875 + x892
    x894 = -x17 * x893
    x895 = x65 * x890 + x894
    x896 = x0 * (2 * x860 - 2 * x861) - x887 + x895
    x897 = x65 * x882
    x898 = x0 * (-3 * x123 - 3 * x631 + 3 * x632 + x878)
    x899 = x65 * x99
    x900 = x126 + x899
    x901 = x17 * (x879 - x898 + x900)
    x902 = -x901
    x903 = x897 + x902
    x904 = -x17 * x903
    x905 = x65 * x885 + x904
    x906 = -x0 * x882 + x875 + x905
    x907 = -x17 * x906
    x908 = x65 * x893 + x907
    x909 = -x0 * (2 * x877 - 2 * x883) + x887 + x908
    x910 = -x17 * x909
    x911 = x65 * x896 + x910
    x912 = (
        x0 * (3 * x852 - x859 + 3 * x864 - 3 * x874)
        - x0 * (x859 - 3 * x875 + 3 * x876 - 3 * x886)
        + x911
    )
    x913 = x0 * x69
    x914 = x0 * x86
    x915 = 2 * x914
    x916 = x49 * x6 + x68
    x917 = x65 * x916
    x918 = x6 * x62 + x85
    x919 = x17 * x918
    x920 = 2 * x913 - x915 + 2 * x917 - 2 * x919
    x921 = x0 * x105
    x922 = 2 * x921
    x923 = x65 * x918
    x924 = x104 + x6 * x83
    x925 = x17 * x924
    x926 = x915 - x922 + 2 * x923 - 2 * x925
    x927 = x0 * (x112 + x916)
    x928 = x0 * (x115 + x918)
    x929 = -x919
    x930 = x913 - x914
    x931 = x65 * (x917 + x929 + x930)
    x932 = -x925
    x933 = x914 - x921
    x934 = x923 + x932 + x933
    x935 = x17 * x934
    x936 = x0 * (x131 + x924)
    x937 = x65 * x934
    x938 = x0 * x128
    x939 = x65 * x924
    x940 = x17 * (x102 * x6 + x127)
    x941 = -x940
    x942 = x17 * (x921 - x938 + x939 + x941)
    x943 = (
        x0 * (x118 + x920)
        - x0 * (x134 + x926)
        - x17 * (x928 - x936 + x937 - x942)
        + x65 * (x927 - x928 + x931 - x935)
    )
    x944 = x0 * x178
    x945 = x0 * x193
    x946 = 2 * x945
    x947 = x161 * x6 + x177
    x948 = x65 * x947
    x949 = x172 * x6 + x192
    x950 = x17 * x949
    x951 = 2 * x944 - x946 + 2 * x948 - 2 * x950
    x952 = x0 * x210
    x953 = 2 * x952
    x954 = x65 * x949
    x955 = x190 * x6 + x209
    x956 = x17 * x955
    x957 = x946 - x953 + 2 * x954 - 2 * x956
    x958 = x0 * (x215 + x947)
    x959 = x0 * (x216 + x949)
    x960 = -x950
    x961 = x65 * (x944 - x945 + x948 + x960)
    x962 = -x956
    x963 = x945 - x952 + x954 + x962
    x964 = x17 * x963
    x965 = x0 * (x227 + x955)
    x966 = x65 * x963
    x967 = x0 * x226
    x968 = x65 * x955
    x969 = x17 * (x207 * x6 + x225)
    x970 = -x969
    x971 = x17 * (x952 - x967 + x968 + x970)
    x972 = x0 * x268
    x973 = x0 * x283
    x974 = 2 * x973
    x975 = x251 * x6 + x267
    x976 = x65 * x975
    x977 = x262 * x6 + x282
    x978 = x17 * x977
    x979 = 2 * x972 - x974 + 2 * x976 - 2 * x978
    x980 = x0 * x300
    x981 = 2 * x980
    x982 = x65 * x977
    x983 = x280 * x6 + x299
    x984 = x17 * x983
    x985 = x974 - x981 + 2 * x982 - 2 * x984
    x986 = x0 * (x305 + x975)
    x987 = x0 * (x306 + x977)
    x988 = -x978
    x989 = x65 * (x972 - x973 + x976 + x988)
    x990 = -x984
    x991 = x973 - x980 + x982 + x990
    x992 = x17 * x991
    x993 = x0 * (x316 + x983)
    x994 = x65 * x991
    x995 = x0 * x315
    x996 = x65 * x983
    x997 = x17 * (x297 * x6 + x314)
    x998 = -x997
    x999 = x17 * (x980 - x995 + x996 + x998)
    x1000 = x0 * x360
    x1001 = x0 * x376
    x1002 = 2 * x1001
    x1003 = x342 * x6 + x359
    x1004 = x1003 * x65
    x1005 = x354 * x6 + x375
    x1006 = x1005 * x17
    x1007 = 2 * x1000 - x1002 + 2 * x1004 - 2 * x1006
    x1008 = x0 * x394
    x1009 = 2 * x1008
    x1010 = x1005 * x65
    x1011 = x373 * x6 + x393
    x1012 = x1011 * x17
    x1013 = x1002 - x1009 + 2 * x1010 - 2 * x1012
    x1014 = x0 * (x1003 + x399)
    x1015 = x0 * (x1005 + x400)
    x1016 = -x1006
    x1017 = x65 * (x1000 - x1001 + x1004 + x1016)
    x1018 = -x1012
    x1019 = x1001 - x1008 + x1010 + x1018
    x1020 = x1019 * x17
    x1021 = x0 * (x1011 + x412)
    x1022 = x1019 * x65
    x1023 = x0 * x411
    x1024 = x1011 * x65
    x1025 = x17 * (x391 * x6 + x410)
    x1026 = -x1025
    x1027 = x17 * (x1008 - x1023 + x1024 + x1026)
    x1028 = x0 * x442
    x1029 = x0 * x454
    x1030 = 2 * x1029
    x1031 = x428 * x6 + x441
    x1032 = x1031 * x65
    x1033 = x436 * x6 + x453
    x1034 = x1033 * x17
    x1035 = 2 * x1028 - x1030 + 2 * x1032 - 2 * x1034
    x1036 = x0 * x468
    x1037 = 2 * x1036
    x1038 = x1033 * x65
    x1039 = x451 * x6 + x467
    x1040 = x1039 * x17
    x1041 = x1030 - x1037 + 2 * x1038 - 2 * x1040
    x1042 = x0 * (x1031 + x473)
    x1043 = x0 * (x1033 + x474)
    x1044 = -x1034
    x1045 = x65 * (x1028 - x1029 + x1032 + x1044)
    x1046 = -x1040
    x1047 = x1029 - x1036 + x1038 + x1046
    x1048 = x1047 * x17
    x1049 = x0 * (x1039 + x482)
    x1050 = x1047 * x65
    x1051 = x0 * x481
    x1052 = x1039 * x65
    x1053 = x17 * (x465 * x6 + x480)
    x1054 = -x1053
    x1055 = x17 * (x1036 - x1051 + x1052 + x1054)
    x1056 = x0 * x509
    x1057 = x0 * x521
    x1058 = 2 * x1057
    x1059 = x495 * x6 + x508
    x1060 = x1059 * x65
    x1061 = x503 * x6 + x520
    x1062 = x1061 * x17
    x1063 = 2 * x1056 - x1058 + 2 * x1060 - 2 * x1062
    x1064 = x0 * x535
    x1065 = 2 * x1064
    x1066 = x1061 * x65
    x1067 = x518 * x6 + x534
    x1068 = x1067 * x17
    x1069 = x1058 - x1065 + 2 * x1066 - 2 * x1068
    x1070 = x0 * (x1059 + x540)
    x1071 = x0 * (x1061 + x541)
    x1072 = -x1062
    x1073 = x65 * (x1056 - x1057 + x1060 + x1072)
    x1074 = -x1068
    x1075 = x1057 - x1064 + x1066 + x1074
    x1076 = x1075 * x17
    x1077 = x0 * (x1067 + x551)
    x1078 = x1075 * x65
    x1079 = x0 * x550
    x1080 = x1067 * x65
    x1081 = x17 * (x532 * x6 + x549)
    x1082 = -x1081
    x1083 = x17 * (x1064 - x1079 + x1080 + x1082)
    x1084 = x0 * x590
    x1085 = x0 * x605
    x1086 = 2 * x1085
    x1087 = x573 * x6 + x589
    x1088 = x1087 * x65
    x1089 = x584 * x6 + x604
    x1090 = x1089 * x17
    x1091 = 2 * x1084 - x1086 + 2 * x1088 - 2 * x1090
    x1092 = x0 * x622
    x1093 = 2 * x1092
    x1094 = x1089 * x65
    x1095 = x6 * x602 + x621
    x1096 = x1095 * x17
    x1097 = x1086 - x1093 + 2 * x1094 - 2 * x1096
    x1098 = x0 * (x1087 + x627)
    x1099 = x0 * (x1089 + x628)
    x1100 = -x1090
    x1101 = x65 * (x1084 - x1085 + x1088 + x1100)
    x1102 = -x1096
    x1103 = x1085 - x1092 + x1094 + x1102
    x1104 = x1103 * x17
    x1105 = x0 * (x1095 + x638)
    x1106 = x1103 * x65
    x1107 = x0 * x637
    x1108 = x1095 * x65
    x1109 = x17 * (x6 * x619 + x636)
    x1110 = -x1109
    x1111 = x17 * (x1092 - x1107 + x1108 + x1110)
    x1112 = x0 * x659
    x1113 = x0 * x669
    x1114 = 2 * x1113
    x1115 = x6 * x647 + x658
    x1116 = x1115 * x65
    x1117 = x6 * x653 + x668
    x1118 = x1117 * x17
    x1119 = 2 * x1112 - x1114 + 2 * x1116 - 2 * x1118
    x1120 = x0 * x681
    x1121 = 2 * x1120
    x1122 = x1117 * x65
    x1123 = x6 * x666 + x680
    x1124 = x1123 * x17
    x1125 = x1114 - x1121 + 2 * x1122 - 2 * x1124
    x1126 = x0 * (x1115 + x686)
    x1127 = x0 * (x1117 + x689)
    x1128 = -x1118
    x1129 = x1112 - x1113
    x1130 = x65 * (x1116 + x1128 + x1129)
    x1131 = -x1124
    x1132 = x1113 - x1120
    x1133 = x1122 + x1131 + x1132
    x1134 = x1133 * x17
    x1135 = x0 * (x1123 + x699)
    x1136 = x1133 * x65
    x1137 = x0 * x696
    x1138 = x1123 * x65
    x1139 = x17 * (x6 * x678 + x695)
    x1140 = -x1139
    x1141 = x17 * (x1120 - x1137 + x1138 + x1140)
    x1142 = x0 * x721
    x1143 = x0 * x730
    x1144 = 2 * x1143
    x1145 = x6 * x710 + x720
    x1146 = x1145 * x65
    x1147 = x6 * x715 + x729
    x1148 = x1147 * x17
    x1149 = 2 * x1142 - x1144 + 2 * x1146 - 2 * x1148
    x1150 = x0 * x741
    x1151 = 2 * x1150
    x1152 = x1147 * x65
    x1153 = x6 * x727 + x740
    x1154 = x1153 * x17
    x1155 = x1144 - x1151 + 2 * x1152 - 2 * x1154
    x1156 = x0 * (x1145 + x748)
    x1157 = x0 * (x1147 + x751)
    x1158 = -x1148
    x1159 = x1142 - x1143
    x1160 = x65 * (x1146 + x1158 + x1159)
    x1161 = -x1154
    x1162 = x1143 - x1150
    x1163 = x1152 + x1161 + x1162
    x1164 = x1163 * x17
    x1165 = x0 * (x1153 + x765)
    x1166 = x1163 * x65
    x1167 = x0 * x762
    x1168 = x1153 * x65
    x1169 = x17 * (x6 * x738 + x761)
    x1170 = -x1169
    x1171 = x17 * (x1150 - x1167 + x1168 + x1170)
    x1172 = x0 * x793
    x1173 = x0 * x804
    x1174 = 2 * x1173
    x1175 = x6 * x780 + x792
    x1176 = x1175 * x65
    x1177 = x6 * x787 + x803
    x1178 = x1177 * x17
    x1179 = 2 * x1172 - x1174 + 2 * x1176 - 2 * x1178
    x1180 = x0 * x817
    x1181 = 2 * x1180
    x1182 = x1177 * x65
    x1183 = x6 * x801 + x816
    x1184 = x1183 * x17
    x1185 = x1174 - x1181 + 2 * x1182 - 2 * x1184
    x1186 = x0 * (x1175 + x822)
    x1187 = x0 * (x1177 + x825)
    x1188 = -x1178
    x1189 = x1172 - x1173
    x1190 = x65 * (x1176 + x1188 + x1189)
    x1191 = -x1184
    x1192 = x1173 - x1180
    x1193 = x1182 + x1191 + x1192
    x1194 = x1193 * x17
    x1195 = x0 * (x1183 + x838)
    x1196 = x1193 * x65
    x1197 = x0 * x835
    x1198 = x1183 * x65
    x1199 = x17 * (x6 * x814 + x834)
    x1200 = -x1199
    x1201 = x17 * (x1180 - x1197 + x1198 + x1200)
    x1202 = x0 * x863
    x1203 = x0 * x873
    x1204 = 2 * x1203
    x1205 = x6 * x851 + x862
    x1206 = x1205 * x65
    x1207 = x6 * x857 + x872
    x1208 = x1207 * x17
    x1209 = 2 * x1202 - x1204 + 2 * x1206 - 2 * x1208
    x1210 = x0 * x885
    x1211 = 2 * x1210
    x1212 = x1207 * x65
    x1213 = x6 * x870 + x884
    x1214 = x1213 * x17
    x1215 = x1204 - x1211 + 2 * x1212 - 2 * x1214
    x1216 = x0 * (x1205 + x890)
    x1217 = x0 * (x1207 + x893)
    x1218 = -x1208
    x1219 = x1202 - x1203
    x1220 = x65 * (x1206 + x1218 + x1219)
    x1221 = -x1214
    x1222 = x1203 - x1210
    x1223 = x1212 + x1221 + x1222
    x1224 = x1223 * x17
    x1225 = x0 * (x1213 + x906)
    x1226 = x1223 * x65
    x1227 = x0 * x903
    x1228 = x1213 * x65
    x1229 = x17 * (x6 * x882 + x902)
    x1230 = -x1229
    x1231 = x17 * (x1210 - x1227 + x1228 + x1230)
    x1232 = 4 * x848
    x1233 = x0 * (-x1232 - 4 * x47 + 4 * x846 + 4 * x849)
    x1234 = 4 * x854
    x1235 = x0 * (x1232 - x1234 - 4 * x60 + 4 * x855)
    x1236 = x1233 - x1235 + x69
    x1237 = x1236 * x65
    x1238 = 4 * x867
    x1239 = x0 * (x1234 - x1238 - 4 * x81 + 4 * x868)
    x1240 = x1235 - x1239 + x86
    x1241 = x1240 * x17
    x1242 = x1240 * x65
    x1243 = 4 * x879
    x1244 = x0 * (-4 * x100 + x1238 - x1243 + 4 * x880)
    x1245 = x105 + x1239 - x1244
    x1246 = x1245 * x17
    x1247 = x0 * x1240
    x1248 = -x1241
    x1249 = x1237 + x1248
    x1250 = -x1246
    x1251 = x1242 + x1250
    x1252 = x0 * x1236 - x1247 + x1249 * x65 - x1251 * x17
    x1253 = -x17 * (-x0 * (x1243 - 4 * x125 - 4 * x898 + 4 * x899) + x1244 + x128)
    x1254 = x1245 * x65 + x1253
    x1255 = -x0 * x1245 + x1247 + x1251 * x65 - x1254 * x17
    x1256 = x151 * x65
    x1257 = x0 * x233
    x1258 = x0 * x237
    x1259 = 2 * x1258
    x1260 = x144 * x65
    x1261 = 2 * x1257 - x1259 + 2 * x1260 - 2 * x149
    x1262 = x0 * x1261
    x1263 = x0 * x244
    x1264 = 2 * x1263
    x1265 = x148 * x65
    x1266 = x1259 - x1264 + 2 * x1265 - 2 * x156
    x1267 = x0 * x1266
    x1268 = 3 * x1267
    x1269 = x0 * (3 * x1256 + 3 * x1262 - x1268 - 3 * x159)
    x1270 = x158 * x65
    x1271 = x0 * x255
    x1272 = 2 * x1271
    x1273 = x155 * x65
    x1274 = x1264 - x1272 + 2 * x1273 - 2 * x167
    x1275 = x0 * x1274
    x1276 = 3 * x1275
    x1277 = x0 * (x1268 + 3 * x1270 - x1276 - 3 * x170)
    x1278 = x1269 - x1277 + x178
    x1279 = x1278 * x65
    x1280 = x169 * x65
    x1281 = x0 * x273
    x1282 = 2 * x1281
    x1283 = x166 * x65
    x1284 = x1272 - x1282 + 2 * x1283 - 2 * x185
    x1285 = x0 * x1284
    x1286 = 3 * x1285
    x1287 = x0 * (x1276 + 3 * x1280 - x1286 - 3 * x188)
    x1288 = x1277 - x1287 + x193
    x1289 = x1288 * x17
    x1290 = x1288 * x65
    x1291 = x187 * x65
    x1292 = x0 * x290
    x1293 = 2 * x1292
    x1294 = x184 * x65
    x1295 = x1282 - x1293 + 2 * x1294 - 2 * x202
    x1296 = x0 * x1295
    x1297 = 3 * x1296
    x1298 = x0 * (x1286 + 3 * x1291 - x1297 - 3 * x205)
    x1299 = x1287 - x1298 + x210
    x1300 = x1299 * x17
    x1301 = x0 * x1288
    x1302 = -x1289
    x1303 = x1279 + x1302
    x1304 = -x1300
    x1305 = x1290 + x1304
    x1306 = x0 * x1278 - x1301 + x1303 * x65 - x1305 * x17
    x1307 = x204 * x65
    x1308 = x0 * x309
    x1309 = x201 * x65
    x1310 = x0 * (x1293 - 2 * x1308 + 2 * x1309 - 2 * x221)
    x1311 = -x17 * (-x0 * (x1297 + 3 * x1307 - 3 * x1310 - 3 * x223) + x1298 + x226)
    x1312 = x1299 * x65 + x1311
    x1313 = -x0 * x1299 + x1301 + x1305 * x65 - x1312 * x17
    x1314 = x0 * x325
    x1315 = 2 * x1314
    x1316 = x0 * x319
    x1317 = -2 * x35 + 2 * x561
    x1318 = -x1315 + 2 * x1316 + x1317
    x1319 = x0 * x1318
    x1320 = x0 * x333
    x1321 = 2 * x1320
    x1322 = -2 * x44 + 2 * x568
    x1323 = x1315 - x1321 + x1322
    x1324 = x0 * x1323
    x1325 = 3 * x1324
    x1326 = x0 * (3 * x1319 - x1325 - 3 * x47 + 3 * x849)
    x1327 = x0 * x345
    x1328 = 2 * x1327
    x1329 = -2 * x57 + 2 * x579
    x1330 = x1321 - x1328 + x1329
    x1331 = x0 * x1330
    x1332 = 3 * x1331
    x1333 = x0 * (x1325 - x1332 - 3 * x60 + 3 * x855)
    x1334 = x1326 - x1333 + x69
    x1335 = x1334 * x65
    x1336 = x0 * x364
    x1337 = 2 * x1336
    x1338 = 2 * x597 - 2 * x78
    x1339 = x1328 - x1337 + x1338
    x1340 = x0 * x1339
    x1341 = 3 * x1340
    x1342 = x0 * (x1332 - x1341 - 3 * x81 + 3 * x868)
    x1343 = x1333 - x1342 + x86
    x1344 = x1343 * x17
    x1345 = x1343 * x65
    x1346 = x0 * x382
    x1347 = 2 * x1346
    x1348 = 2 * x614 - 2 * x97
    x1349 = x1337 - x1347 + x1348
    x1350 = x0 * x1349
    x1351 = 3 * x1350
    x1352 = x0 * (-3 * x100 + x1341 - x1351 + 3 * x880)
    x1353 = x105 + x1342 - x1352
    x1354 = x1353 * x17
    x1355 = x0 * x1343
    x1356 = -x1344
    x1357 = x1335 + x1356
    x1358 = -x1354
    x1359 = x1345 + x1358
    x1360 = x0 * x1334 - x1355 + x1357 * x65 - x1359 * x17
    x1361 = x0 * x405
    x1362 = x0 * (-2 * x123 + x1347 - 2 * x1361 + 2 * x632)
    x1363 = -x17 * (-x0 * (-3 * x125 + x1351 - 3 * x1362 + 3 * x899) + x128 + x1352)
    x1364 = x1353 * x65 + x1363
    x1365 = -x0 * x1353 + x1355 + x1359 * x65 - x1364 * x17
    x1366 = x0 * x420
    x1367 = x0 * x424
    x1368 = 2 * x1367
    x1369 = x241 * x65
    x1370 = 2 * x1366 - x1368 + 2 * x1369 - 2 * x249
    x1371 = x0 * x1370
    x1372 = x0 * x432
    x1373 = 2 * x1372
    x1374 = x248 * x65
    x1375 = x1368 - x1373 + 2 * x1374 - 2 * x260
    x1376 = x0 * x1375
    x1377 = x1371 - x1376 + x268
    x1378 = x1377 * x65
    x1379 = x0 * x447
    x1380 = 2 * x1379
    x1381 = x259 * x65
    x1382 = x1373 - x1380 + 2 * x1381 - 2 * x278
    x1383 = x0 * x1382
    x1384 = x1376 - x1383 + x283
    x1385 = x1384 * x17
    x1386 = x1384 * x65
    x1387 = x0 * x461
    x1388 = 2 * x1387
    x1389 = x277 * x65
    x1390 = x1380 - x1388 + 2 * x1389 - 2 * x295
    x1391 = x0 * x1390
    x1392 = x1383 - x1391 + x300
    x1393 = x1392 * x17
    x1394 = x0 * x1384
    x1395 = -x1385
    x1396 = x1378 + x1395
    x1397 = -x1393
    x1398 = x1386 + x1397
    x1399 = x0 * x1377 - x1394 + x1396 * x65 - x1398 * x17
    x1400 = x0 * x477
    x1401 = x294 * x65
    x1402 = -x17 * (-x0 * (x1388 - 2 * x1400 + 2 * x1401 - 2 * x312) + x1391 + x315)
    x1403 = x1392 * x65 + x1402
    x1404 = -x0 * x1392 + x1394 + x1398 * x65 - x1403 * x17
    x1405 = x0 * x487
    x1406 = x0 * x491
    x1407 = 2 * x1406
    x1408 = x331 * x65
    x1409 = 2 * x1405 - x1407 + 2 * x1408 - 2 * x340
    x1410 = x0 * x1409
    x1411 = x0 * x499
    x1412 = 2 * x1411
    x1413 = x339 * x65
    x1414 = x1407 - x1412 + 2 * x1413 - 2 * x352
    x1415 = x0 * x1414
    x1416 = x1410 - x1415 + x360
    x1417 = x1416 * x65
    x1418 = x0 * x514
    x1419 = 2 * x1418
    x1420 = x351 * x65
    x1421 = x1412 - x1419 + 2 * x1420 - 2 * x371
    x1422 = x0 * x1421
    x1423 = x1415 - x1422 + x376
    x1424 = x1423 * x17
    x1425 = x1423 * x65
    x1426 = x0 * x528
    x1427 = 2 * x1426
    x1428 = x370 * x65
    x1429 = x1419 - x1427 + 2 * x1428 - 2 * x389
    x1430 = x0 * x1429
    x1431 = x1422 - x1430 + x394
    x1432 = x1431 * x17
    x1433 = x0 * x1423
    x1434 = -x1424
    x1435 = x1417 + x1434
    x1436 = -x1432
    x1437 = x1425 + x1436
    x1438 = x0 * x1416 - x1433 + x1435 * x65 - x1437 * x17
    x1439 = x0 * x546
    x1440 = x388 * x65
    x1441 = -x17 * (-x0 * (x1427 - 2 * x1439 + 2 * x1440 - 2 * x408) + x1430 + x411)
    x1442 = x1431 * x65 + x1441
    x1443 = -x0 * x1431 + x1433 + x1437 * x65 - x1442 * x17
    x1444 = x0 * x569
    x1445 = 2 * x1444
    x1446 = x0 * x562
    x1447 = -2 * x47 + 2 * x849
    x1448 = -x1445 + 2 * x1446 + x1447
    x1449 = x0 * x1448
    x1450 = x0 * x580
    x1451 = 2 * x1450
    x1452 = -2 * x60 + 2 * x855
    x1453 = x1445 - x1451 + x1452
    x1454 = x0 * x1453
    x1455 = x1449 - x1454 + x69
    x1456 = x1455 * x65
    x1457 = x0 * x598
    x1458 = 2 * x1457
    x1459 = -2 * x81 + 2 * x868
    x1460 = x1451 - x1458 + x1459
    x1461 = x0 * x1460
    x1462 = x1454 - x1461 + x86
    x1463 = x1462 * x17
    x1464 = x1462 * x65
    x1465 = x0 * x615
    x1466 = 2 * x1465
    x1467 = -2 * x100 + 2 * x880
    x1468 = x1458 - x1466 + x1467
    x1469 = x0 * x1468
    x1470 = x105 + x1461 - x1469
    x1471 = x1470 * x17
    x1472 = x0 * x1462
    x1473 = -x1463
    x1474 = x1456 + x1473
    x1475 = -x1471
    x1476 = x1464 + x1475
    x1477 = x0 * x1455 - x1472 + x1474 * x65 - x1476 * x17
    x1478 = x0 * x633
    x1479 = -x17 * (-x0 * (-2 * x125 + x1466 - 2 * x1478 + 2 * x899) + x128 + x1469)
    x1480 = x1470 * x65 + x1479
    x1481 = -x0 * x1470 + x1472 + x1476 * x65 - x1480 * x17
    x1482 = x0 * x646
    x1483 = x0 * x652
    x1484 = x1482 - x1483 + x442
    x1485 = x1484 * x65
    x1486 = x0 * x665
    x1487 = x1483 - x1486 + x454
    x1488 = x1487 * x17
    x1489 = x1487 * x65
    x1490 = x0 * x677
    x1491 = x1486 - x1490 + x468
    x1492 = x1491 * x17
    x1493 = x0 * x1487
    x1494 = -x1488
    x1495 = x1485 + x1494
    x1496 = -x1492
    x1497 = x1489 + x1496
    x1498 = x0 * x1484 - x1493 + x1495 * x65 - x1497 * x17
    x1499 = -x17 * (-x0 * x693 + x1490 + x481)
    x1500 = x1491 * x65 + x1499
    x1501 = -x0 * x1491 + x1493 + x1497 * x65 - x1500 * x17
    x1502 = x0 * x709
    x1503 = x0 * x714
    x1504 = x1502 - x1503 + x509
    x1505 = x1504 * x65
    x1506 = x0 * x726
    x1507 = x1503 - x1506 + x521
    x1508 = x1507 * x17
    x1509 = x1507 * x65
    x1510 = x0 * x737
    x1511 = x1506 - x1510 + x535
    x1512 = x1511 * x17
    x1513 = x0 * x1507
    x1514 = -x1508
    x1515 = x1505 + x1514
    x1516 = -x1512
    x1517 = x1509 + x1516
    x1518 = x0 * x1504 - x1513 + x1515 * x65 - x1517 * x17
    x1519 = -x17 * (-x0 * x759 + x1510 + x550)
    x1520 = x1511 * x65 + x1519
    x1521 = -x0 * x1511 + x1513 + x1517 * x65 - x1520 * x17
    x1522 = x0 * x779
    x1523 = x0 * x786
    x1524 = x1522 - x1523 + x590
    x1525 = x1524 * x65
    x1526 = x0 * x800
    x1527 = x1523 - x1526 + x605
    x1528 = x1527 * x17
    x1529 = x1527 * x65
    x1530 = x0 * x813
    x1531 = x1526 - x1530 + x622
    x1532 = x1531 * x17
    x1533 = x0 * x1527
    x1534 = -x1528
    x1535 = x1525 + x1534
    x1536 = -x1532
    x1537 = x1529 + x1536
    x1538 = x0 * x1524 - x1533 + x1535 * x65 - x1537 * x17
    x1539 = -x17 * (-x0 * x832 + x1530 + x637)
    x1540 = x1531 * x65 + x1539
    x1541 = -x0 * x1531 + x1533 + x1537 * x65 - x1540 * x17
    x1542 = x0 * x850
    x1543 = x0 * x856
    x1544 = x1542 - x1543 + x69
    x1545 = x1544 * x65
    x1546 = x0 * x869
    x1547 = x1543 - x1546 + x86
    x1548 = x1547 * x17
    x1549 = x1547 * x65
    x1550 = x0 * x881
    x1551 = x105 + x1546 - x1550
    x1552 = x1551 * x17
    x1553 = x0 * x1547
    x1554 = -x1548
    x1555 = x1545 + x1554
    x1556 = -x1552
    x1557 = x1549 + x1556
    x1558 = x0 * x1544 - x1553 + x1555 * x65 - x1557 * x17
    x1559 = -x17 * (-x0 * x900 + x128 + x1550)
    x1560 = x1551 * x65 + x1559
    x1561 = -x0 * x1551 + x1553 + x1557 * x65 - x1560 * x17
    x1562 = x1129 - x17 * x688 + x65 * x685
    x1563 = x1132 - x17 * x698 + x65 * x688
    x1564 = x1159 - x17 * x750 + x65 * x747
    x1565 = x1162 - x17 * x764 + x65 * x750
    x1566 = x1189 - x17 * x824 + x65 * x821
    x1567 = x1192 - x17 * x837 + x65 * x824
    x1568 = x1219 - x17 * x892 + x65 * x889
    x1569 = x1222 - x17 * x905 + x65 * x892
    x1570 = x111 * x65 - x114 * x17 + x930
    x1571 = x114 * x65 - x130 * x17 + x933
    x1572 = x6 * x916 + x929
    x1573 = x1572 * x65
    x1574 = x6 * x918 + x932
    x1575 = x1574 * x17
    x1576 = -x1575
    x1577 = x0 * x920
    x1578 = x0 * x926
    x1579 = x1574 * x65
    x1580 = x17 * (x6 * x924 + x941)
    x1581 = -x1580
    x1582 = x0 * (x922 - 2 * x938 + 2 * x939 - 2 * x940)
    x1583 = 2 * x928
    x1584 = (
        x0 * (x1572 - x1583 + 2 * x927 + 2 * x931 - 2 * x935)
        - x0 * (x1574 + x1583 - 2 * x936 + 2 * x937 - 2 * x942)
        - x17 * (x1578 + x1579 + x1581 - x1582)
        + x65 * (x1573 + x1576 + x1577 - x1578)
    )
    x1585 = x6 * x947 + x960
    x1586 = x1585 * x65
    x1587 = x6 * x949 + x962
    x1588 = x1587 * x17
    x1589 = -x1588
    x1590 = x0 * x951
    x1591 = x0 * x957
    x1592 = x1587 * x65
    x1593 = x17 * (x6 * x955 + x970)
    x1594 = -x1593
    x1595 = x0 * (x953 - 2 * x967 + 2 * x968 - 2 * x969)
    x1596 = 2 * x959
    x1597 = x6 * x975 + x988
    x1598 = x1597 * x65
    x1599 = x6 * x977 + x990
    x1600 = x1599 * x17
    x1601 = -x1600
    x1602 = x0 * x979
    x1603 = x0 * x985
    x1604 = x1599 * x65
    x1605 = x17 * (x6 * x983 + x998)
    x1606 = -x1605
    x1607 = x0 * (x981 - 2 * x995 + 2 * x996 - 2 * x997)
    x1608 = 2 * x987
    x1609 = x1003 * x6 + x1016
    x1610 = x1609 * x65
    x1611 = x1005 * x6 + x1018
    x1612 = x1611 * x17
    x1613 = -x1612
    x1614 = x0 * x1007
    x1615 = x0 * x1013
    x1616 = x1611 * x65
    x1617 = x17 * (x1011 * x6 + x1026)
    x1618 = -x1617
    x1619 = x0 * (x1009 - 2 * x1023 + 2 * x1024 - 2 * x1025)
    x1620 = 2 * x1015
    x1621 = x1031 * x6 + x1044
    x1622 = x1621 * x65
    x1623 = x1033 * x6 + x1046
    x1624 = x1623 * x17
    x1625 = -x1624
    x1626 = x0 * x1035
    x1627 = x0 * x1041
    x1628 = x1623 * x65
    x1629 = x17 * (x1039 * x6 + x1054)
    x1630 = -x1629
    x1631 = x0 * (x1037 - 2 * x1051 + 2 * x1052 - 2 * x1053)
    x1632 = 2 * x1043
    x1633 = x1059 * x6 + x1072
    x1634 = x1633 * x65
    x1635 = x1061 * x6 + x1074
    x1636 = x1635 * x17
    x1637 = -x1636
    x1638 = x0 * x1063
    x1639 = x0 * x1069
    x1640 = x1635 * x65
    x1641 = x17 * (x1067 * x6 + x1082)
    x1642 = -x1641
    x1643 = x0 * (x1065 - 2 * x1079 + 2 * x1080 - 2 * x1081)
    x1644 = 2 * x1071
    x1645 = x1087 * x6 + x1100
    x1646 = x1645 * x65
    x1647 = x1089 * x6 + x1102
    x1648 = x1647 * x17
    x1649 = -x1648
    x1650 = x0 * x1091
    x1651 = x0 * x1097
    x1652 = x1647 * x65
    x1653 = x17 * (x1095 * x6 + x1110)
    x1654 = -x1653
    x1655 = x0 * (x1093 - 2 * x1107 + 2 * x1108 - 2 * x1109)
    x1656 = 2 * x1099
    x1657 = x1115 * x6 + x1128
    x1658 = x1657 * x65
    x1659 = x1117 * x6 + x1131
    x1660 = x1659 * x17
    x1661 = -x1660
    x1662 = x0 * x1119
    x1663 = x0 * x1125
    x1664 = x1659 * x65
    x1665 = x17 * (x1123 * x6 + x1140)
    x1666 = -x1665
    x1667 = x0 * (x1121 - 2 * x1137 + 2 * x1138 - 2 * x1139)
    x1668 = 2 * x1127
    x1669 = x1145 * x6 + x1158
    x1670 = x1669 * x65
    x1671 = x1147 * x6 + x1161
    x1672 = x1671 * x17
    x1673 = -x1672
    x1674 = x0 * x1149
    x1675 = x0 * x1155
    x1676 = x1671 * x65
    x1677 = x17 * (x1153 * x6 + x1170)
    x1678 = -x1677
    x1679 = x0 * (x1151 - 2 * x1167 + 2 * x1168 - 2 * x1169)
    x1680 = 2 * x1157
    x1681 = x1175 * x6 + x1188
    x1682 = x1681 * x65
    x1683 = x1177 * x6 + x1191
    x1684 = x1683 * x17
    x1685 = -x1684
    x1686 = x0 * x1179
    x1687 = x0 * x1185
    x1688 = x1683 * x65
    x1689 = x17 * (x1183 * x6 + x1200)
    x1690 = -x1689
    x1691 = x0 * (x1181 - 2 * x1197 + 2 * x1198 - 2 * x1199)
    x1692 = 2 * x1187
    x1693 = x1205 * x6 + x1218
    x1694 = x1693 * x65
    x1695 = x1207 * x6 + x1221
    x1696 = x1695 * x17
    x1697 = -x1696
    x1698 = x0 * x1209
    x1699 = x0 * x1215
    x1700 = x1695 * x65
    x1701 = x17 * (x1213 * x6 + x1230)
    x1702 = -x1701
    x1703 = x0 * (x1211 - 2 * x1227 + 2 * x1228 - 2 * x1229)
    x1704 = 2 * x1217
    x1705 = x1236 * x6 + x1248
    x1706 = x1240 * x6 + x1250
    x1707 = x0 * x1249
    x1708 = x0 * x1251
    x1709 = x1705 * x65
    x1710 = x17 * x1706
    x1711 = -x1710
    x1712 = x0 * x1254
    x1713 = x1706 * x65
    x1714 = x17 * (x1245 * x6 + x1253)
    x1715 = -x1714
    x1716 = x1278 * x6 + x1302
    x1717 = x1288 * x6 + x1304
    x1718 = x0 * x1303
    x1719 = x0 * x1305
    x1720 = x1716 * x65
    x1721 = x17 * x1717
    x1722 = -x1721
    x1723 = x0 * x1312
    x1724 = x1717 * x65
    x1725 = x17 * (x1299 * x6 + x1311)
    x1726 = -x1725
    x1727 = x1334 * x6 + x1356
    x1728 = x1343 * x6 + x1358
    x1729 = x0 * x1357
    x1730 = x0 * x1359
    x1731 = x1727 * x65
    x1732 = x17 * x1728
    x1733 = -x1732
    x1734 = x0 * x1364
    x1735 = x1728 * x65
    x1736 = x17 * (x1353 * x6 + x1363)
    x1737 = -x1736
    x1738 = x1377 * x6 + x1395
    x1739 = x1384 * x6 + x1397
    x1740 = x0 * x1396
    x1741 = x0 * x1398
    x1742 = x1738 * x65
    x1743 = x17 * x1739
    x1744 = -x1743
    x1745 = x0 * x1403
    x1746 = x1739 * x65
    x1747 = x17 * (x1392 * x6 + x1402)
    x1748 = -x1747
    x1749 = x1416 * x6 + x1434
    x1750 = x1423 * x6 + x1436
    x1751 = x0 * x1435
    x1752 = x0 * x1437
    x1753 = x1749 * x65
    x1754 = x17 * x1750
    x1755 = -x1754
    x1756 = x0 * x1442
    x1757 = x1750 * x65
    x1758 = x17 * (x1431 * x6 + x1441)
    x1759 = -x1758
    x1760 = x1455 * x6 + x1473
    x1761 = x1462 * x6 + x1475
    x1762 = x0 * x1474
    x1763 = x0 * x1476
    x1764 = x1760 * x65
    x1765 = x17 * x1761
    x1766 = -x1765
    x1767 = x0 * x1480
    x1768 = x1761 * x65
    x1769 = x17 * (x1470 * x6 + x1479)
    x1770 = -x1769
    x1771 = x1484 * x6 + x1494
    x1772 = x1487 * x6 + x1496
    x1773 = x0 * x1495
    x1774 = x0 * x1497
    x1775 = x1771 * x65
    x1776 = x17 * x1772
    x1777 = -x1776
    x1778 = x0 * x1500
    x1779 = x1772 * x65
    x1780 = x17 * (x1491 * x6 + x1499)
    x1781 = -x1780
    x1782 = x1504 * x6 + x1514
    x1783 = x1507 * x6 + x1516
    x1784 = x0 * x1515
    x1785 = x0 * x1517
    x1786 = x1782 * x65
    x1787 = x17 * x1783
    x1788 = -x1787
    x1789 = x0 * x1520
    x1790 = x1783 * x65
    x1791 = x17 * (x1511 * x6 + x1519)
    x1792 = -x1791
    x1793 = x1524 * x6 + x1534
    x1794 = x1527 * x6 + x1536
    x1795 = x0 * x1535
    x1796 = x0 * x1537
    x1797 = x1793 * x65
    x1798 = x17 * x1794
    x1799 = -x1798
    x1800 = x0 * x1540
    x1801 = x1794 * x65
    x1802 = x17 * (x1531 * x6 + x1539)
    x1803 = -x1802
    x1804 = x1544 * x6 + x1554
    x1805 = x1547 * x6 + x1556
    x1806 = x0 * x1555
    x1807 = x0 * x1557
    x1808 = x1804 * x65
    x1809 = x17 * x1805
    x1810 = -x1809
    x1811 = x0 * x1560
    x1812 = x1805 * x65
    x1813 = x17 * (x1551 * x6 + x1559)
    x1814 = -x1813
    x1815 = x6 * x659 + x684
    x1816 = x6 * x669 + x687
    x1817 = x0 * x685
    x1818 = x0 * x688
    x1819 = x1815 * x65
    x1820 = x17 * x1816
    x1821 = -x1820
    x1822 = x0 * x698
    x1823 = x1816 * x65
    x1824 = x17 * (x6 * x681 + x697)
    x1825 = -x1824
    x1826 = x6 * x721 + x746
    x1827 = x6 * x730 + x749
    x1828 = x0 * x747
    x1829 = x0 * x750
    x1830 = x1826 * x65
    x1831 = x17 * x1827
    x1832 = -x1831
    x1833 = x0 * x764
    x1834 = x1827 * x65
    x1835 = x17 * (x6 * x741 + x763)
    x1836 = -x1835
    x1837 = x6 * x793 + x820
    x1838 = x6 * x804 + x823
    x1839 = x0 * x821
    x1840 = x0 * x824
    x1841 = x1837 * x65
    x1842 = x17 * x1838
    x1843 = -x1842
    x1844 = x0 * x837
    x1845 = x1838 * x65
    x1846 = x17 * (x6 * x817 + x836)
    x1847 = -x1846
    x1848 = x6 * x863 + x888
    x1849 = x6 * x873 + x891
    x1850 = x0 * x889
    x1851 = x0 * x892
    x1852 = x1848 * x65
    x1853 = x17 * x1849
    x1854 = -x1853
    x1855 = x0 * x905
    x1856 = x1849 * x65
    x1857 = x17 * (x6 * x885 + x904)
    x1858 = -x1857
    x1859 = x110 + x6 * x69
    x1860 = x113 + x6 * x86
    x1861 = x0 * x111
    x1862 = x0 * x114
    x1863 = x1859 * x65
    x1864 = x17 * x1860
    x1865 = -x1864
    x1866 = x0 * x130
    x1867 = x1860 * x65
    x1868 = x17 * (x105 * x6 + x129)
    x1869 = -x1868
    x1870 = 3 * x784
    x1871 = 3 * x777
    x1872 = -x1870 + x1871 - 3 * x582 + 3 * x785
    x1873 = x0 * (x1872 + x46)
    x1874 = 4 * x1873
    x1875 = -x1871 - 3 * x571 + 3 * x774 + 3 * x778
    x1876 = x0 * (x1875 + x37)
    x1877 = x0 * (-x1874 + 4 * x1876 + x49 + 4 * x860 - 4 * x861)
    x1878 = 3 * x798
    x1879 = x1870 - x1878 - 3 * x600 + 3 * x799
    x1880 = x0 * (x1879 + x59)
    x1881 = 4 * x1880
    x1882 = x0 * (x1874 - x1881 + x62 + 4 * x865 - 4 * x871)
    x1883 = x1249 + x1877 - x1882
    x1884 = 3 * x811
    x1885 = x1878 - x1884 - 3 * x617 + 3 * x812
    x1886 = x0 * (x1885 + x80)
    x1887 = 4 * x1886
    x1888 = x0 * (x1881 - x1887 + x83 + 4 * x877 - 4 * x883)
    x1889 = x1251 + x1882 - x1888
    x1890 = -x17 * x1889
    x1891 = x1883 * x65 + x1890
    x1892 = x0 * (x1884 - 3 * x634 - 3 * x830 + 3 * x831 + x99)
    x1893 = -x17 * (
        -x0 * (x102 + x1887 - 4 * x1892 + 4 * x897 - 4 * x901) + x1254 + x1888
    )
    x1894 = x1889 * x65 + x1893
    x1895 = x244 * x65
    x1896 = x17 * x255
    x1897 = x1895 - x1896 + x335
    x1898 = x0 * (x155 + x1897)
    x1899 = 2 * x1898
    x1900 = x1263 - x1271 + x1273 + x168
    x1901 = x17 * x1900
    x1902 = x237 * x65
    x1903 = x17 * x244
    x1904 = x1902 - x1903 + x327
    x1905 = x0 * (x148 + x1904)
    x1906 = 2 * x1905
    x1907 = x1258 - x1263 + x1265 + x157
    x1908 = x1907 * x65
    x1909 = x0 * (x158 - x1899 - 2 * x1901 + x1906 + 2 * x1908)
    x1910 = 3 * x1909
    x1911 = x1267 + x1270 - x1275 + x171
    x1912 = x17 * x1911
    x1913 = x17 * x1907
    x1914 = x233 * x65
    x1915 = x17 * x237
    x1916 = x1914 - x1915 + x322
    x1917 = x0 * (x144 + x1916)
    x1918 = x65 * (x1257 - x1258 + x1260 + x150)
    x1919 = x0 * (x151 - x1906 - 2 * x1913 + 2 * x1917 + 2 * x1918)
    x1920 = x65 * (x1256 + x1262 - x1267 + x160)
    x1921 = -x1910 - 3 * x1912 + 3 * x1919 + 3 * x1920
    x1922 = x0 * (x161 + x1921)
    x1923 = x255 * x65
    x1924 = x17 * x273
    x1925 = x1923 - x1924 + x347
    x1926 = x0 * (x166 + x1925)
    x1927 = 2 * x1926
    x1928 = x1271 - x1281 + x1283 + x186
    x1929 = x17 * x1928
    x1930 = x1900 * x65
    x1931 = x0 * (x169 + x1899 - x1927 - 2 * x1929 + 2 * x1930)
    x1932 = 3 * x1931
    x1933 = x1275 + x1280 - x1285 + x189
    x1934 = x17 * x1933
    x1935 = x1911 * x65
    x1936 = x1910 - x1932 - 3 * x1934 + 3 * x1935
    x1937 = x0 * (x172 + x1936)
    x1938 = x1303 + x1922 - x1937
    x1939 = x273 * x65
    x1940 = x17 * x290
    x1941 = x1939 - x1940 + x366
    x1942 = x0 * (x184 + x1941)
    x1943 = 2 * x1942
    x1944 = x1281 - x1292 + x1294 + x203
    x1945 = x17 * x1944
    x1946 = x1928 * x65
    x1947 = x0 * (x187 + x1927 - x1943 - 2 * x1945 + 2 * x1946)
    x1948 = 3 * x1947
    x1949 = x1285 + x1291 - x1296 + x206
    x1950 = x17 * x1949
    x1951 = x1933 * x65
    x1952 = x1932 - x1948 - 3 * x1950 + 3 * x1951
    x1953 = x0 * (x190 + x1952)
    x1954 = x1305 + x1937 - x1953
    x1955 = -x17 * x1954
    x1956 = x1938 * x65 + x1955
    x1957 = -x17 * x309 + x290 * x65 + x384
    x1958 = x0 * (x1957 + x201)
    x1959 = x17 * (x1292 - x1308 + x1309 + x222)
    x1960 = x1944 * x65
    x1961 = x0 * (x1943 - 2 * x1958 - 2 * x1959 + 2 * x1960 + x204)
    x1962 = x17 * (x1296 + x1307 - x1310 + x224)
    x1963 = x1949 * x65
    x1964 = -x17 * (
        -x0 * (x1948 - 3 * x1961 - 3 * x1962 + 3 * x1963 + x207) + x1312 + x1953
    )
    x1965 = x1954 * x65 + x1964
    x1966 = x0 * x31
    x1967 = x0 * x40
    x1968 = x333 * x65
    x1969 = x17 * x345
    x1970 = x1966 - x1967 + x1968 - x1969
    x1971 = x0 * (x1970 + x43)
    x1972 = 2 * x1971
    x1973 = x1320 - x1327 + x580
    x1974 = x17 * x1973
    x1975 = x0 * x25
    x1976 = x325 * x65
    x1977 = x17 * x333
    x1978 = -x1966 + x1975 + x1976 - x1977
    x1979 = x0 * (x1978 + x34)
    x1980 = 2 * x1979
    x1981 = x1314 - x1320 + x569
    x1982 = x1981 * x65
    x1983 = x0 * (-x1972 - 2 * x1974 + x1980 + 2 * x1982 + x46)
    x1984 = 3 * x1983
    x1985 = x1324 - x1331 + x856
    x1986 = x17 * x1985
    x1987 = x17 * x1981
    x1988 = x0 * x21
    x1989 = x319 * x65
    x1990 = x17 * x325
    x1991 = -x1975 + x1988 + x1989 - x1990
    x1992 = x0 * (x1991 + x28)
    x1993 = x65 * (-x1314 + x1316 + x562)
    x1994 = x0 * (-x1980 - 2 * x1987 + 2 * x1992 + 2 * x1993 + x37)
    x1995 = x65 * (x1319 - x1324 + x850)
    x1996 = -x1984 - 3 * x1986 + 3 * x1994 + 3 * x1995
    x1997 = x0 * (x1996 + x49)
    x1998 = x0 * x53
    x1999 = x345 * x65
    x2000 = x17 * x364
    x2001 = x1967 - x1998 + x1999 - x2000
    x2002 = x0 * (x2001 + x56)
    x2003 = 2 * x2002
    x2004 = x1327 - x1336 + x598
    x2005 = x17 * x2004
    x2006 = x1973 * x65
    x2007 = x0 * (x1972 - x2003 - 2 * x2005 + 2 * x2006 + x59)
    x2008 = 3 * x2007
    x2009 = x1331 - x1340 + x869
    x2010 = x17 * x2009
    x2011 = x1985 * x65
    x2012 = x1984 - x2008 - 3 * x2010 + 3 * x2011
    x2013 = x0 * (x2012 + x62)
    x2014 = x1357 + x1997 - x2013
    x2015 = x0 * x74
    x2016 = x364 * x65
    x2017 = x17 * x382
    x2018 = x1998 - x2015 + x2016 - x2017
    x2019 = x0 * (x2018 + x77)
    x2020 = 2 * x2019
    x2021 = x1336 - x1346 + x615
    x2022 = x17 * x2021
    x2023 = x2004 * x65
    x2024 = x0 * (x2003 - x2020 - 2 * x2022 + 2 * x2023 + x80)
    x2025 = 3 * x2024
    x2026 = x1340 - x1350 + x881
    x2027 = x17 * x2026
    x2028 = x2009 * x65
    x2029 = x2008 - x2025 - 3 * x2027 + 3 * x2028
    x2030 = x0 * (x2029 + x83)
    x2031 = x1359 + x2013 - x2030
    x2032 = -x17 * x2031
    x2033 = x2014 * x65 + x2032
    x2034 = -x0 * x93 - x17 * x405 + x2015 + x382 * x65
    x2035 = x0 * (x2034 + x96)
    x2036 = x17 * (x1346 - x1361 + x633)
    x2037 = x2021 * x65
    x2038 = x0 * (x2020 - 2 * x2035 - 2 * x2036 + 2 * x2037 + x99)
    x2039 = x17 * (x1350 - x1362 + x900)
    x2040 = x2026 * x65
    x2041 = -x17 * (
        -x0 * (x102 + x2025 - 3 * x2038 - 3 * x2039 + 3 * x2040) + x1364 + x2030
    )
    x2042 = x2031 * x65 + x2041
    x2043 = x0 * x238
    x2044 = x0 * x245
    x2045 = x424 * x65
    x2046 = x17 * x432
    x2047 = x2043 - x2044 + x2045 - x2046
    x2048 = x0 * (x2047 + x248)
    x2049 = 2 * x2048
    x2050 = x1367 - x1372 + x1374 + x261
    x2051 = x17 * x2050
    x2052 = x0 * x234
    x2053 = x420 * x65
    x2054 = x17 * x424
    x2055 = -x2043 + x2052 + x2053 - x2054
    x2056 = x0 * (x2055 + x241)
    x2057 = x65 * (x1366 - x1367 + x1369 + x250)
    x2058 = x0 * (-x2049 - 2 * x2051 + 2 * x2056 + 2 * x2057 + x251)
    x2059 = x0 * x256
    x2060 = x432 * x65
    x2061 = x17 * x447
    x2062 = x2044 - x2059 + x2060 - x2061
    x2063 = x0 * (x2062 + x259)
    x2064 = 2 * x2063
    x2065 = x1372 - x1379 + x1381 + x279
    x2066 = x17 * x2065
    x2067 = x2050 * x65
    x2068 = x0 * (x2049 - x2064 - 2 * x2066 + 2 * x2067 + x262)
    x2069 = x1396 + x2058 - x2068
    x2070 = x0 * x274
    x2071 = x447 * x65
    x2072 = x17 * x461
    x2073 = x2059 - x2070 + x2071 - x2072
    x2074 = x0 * (x2073 + x277)
    x2075 = 2 * x2074
    x2076 = x1379 - x1387 + x1389 + x296
    x2077 = x17 * x2076
    x2078 = x2065 * x65
    x2079 = x0 * (x2064 - x2075 - 2 * x2077 + 2 * x2078 + x280)
    x2080 = x1398 + x2068 - x2079
    x2081 = -x17 * x2080
    x2082 = x2069 * x65 + x2081
    x2083 = -x0 * x291 - x17 * x477 + x2070 + x461 * x65
    x2084 = x0 * (x2083 + x294)
    x2085 = x17 * (x1387 - x1400 + x1401 + x313)
    x2086 = x2076 * x65
    x2087 = -x17 * (
        -x0 * (x2075 - 2 * x2084 - 2 * x2085 + 2 * x2086 + x297) + x1403 + x2079
    )
    x2088 = x2080 * x65 + x2087
    x2089 = x0 * x328
    x2090 = x0 * x336
    x2091 = x491 * x65
    x2092 = x17 * x499
    x2093 = x2089 - x2090 + x2091 - x2092
    x2094 = x0 * (x2093 + x339)
    x2095 = 2 * x2094
    x2096 = x1406 - x1411 + x1413 + x353
    x2097 = x17 * x2096
    x2098 = x0 * x323
    x2099 = x487 * x65
    x2100 = x17 * x491
    x2101 = -x2089 + x2098 + x2099 - x2100
    x2102 = x0 * (x2101 + x331)
    x2103 = x65 * (x1405 - x1406 + x1408 + x341)
    x2104 = x0 * (-x2095 - 2 * x2097 + 2 * x2102 + 2 * x2103 + x342)
    x2105 = x0 * x348
    x2106 = x499 * x65
    x2107 = x17 * x514
    x2108 = x2090 - x2105 + x2106 - x2107
    x2109 = x0 * (x2108 + x351)
    x2110 = 2 * x2109
    x2111 = x1411 - x1418 + x1420 + x372
    x2112 = x17 * x2111
    x2113 = x2096 * x65
    x2114 = x0 * (x2095 - x2110 - 2 * x2112 + 2 * x2113 + x354)
    x2115 = x1435 + x2104 - x2114
    x2116 = x0 * x367
    x2117 = x514 * x65
    x2118 = x17 * x528
    x2119 = x2105 - x2116 + x2117 - x2118
    x2120 = x0 * (x2119 + x370)
    x2121 = 2 * x2120
    x2122 = x1418 - x1426 + x1428 + x390
    x2123 = x17 * x2122
    x2124 = x2111 * x65
    x2125 = x0 * (x2110 - x2121 - 2 * x2123 + 2 * x2124 + x373)
    x2126 = x1437 + x2114 - x2125
    x2127 = -x17 * x2126
    x2128 = x2115 * x65 + x2127
    x2129 = -x0 * x385 - x17 * x546 + x2116 + x528 * x65
    x2130 = x0 * (x2129 + x388)
    x2131 = x17 * (x1426 - x1439 + x1440 + x409)
    x2132 = x2122 * x65
    x2133 = -x17 * (
        -x0 * (x2121 - 2 * x2130 - 2 * x2131 + 2 * x2132 + x391) + x1442 + x2125
    )
    x2134 = x2126 * x65 + x2133
    x2135 = x0 * x34
    x2136 = x0 * x43
    x2137 = x569 * x65
    x2138 = x17 * x580
    x2139 = x2135 - x2136 + x2137 - x2138
    x2140 = x0 * (x2139 + x46)
    x2141 = 2 * x2140
    x2142 = x1444 - x1450 + x856
    x2143 = x17 * x2142
    x2144 = x0 * x28
    x2145 = x562 * x65
    x2146 = x17 * x569
    x2147 = -x2135 + x2144 + x2145 - x2146
    x2148 = x0 * (x2147 + x37)
    x2149 = x65 * (-x1444 + x1446 + x850)
    x2150 = x0 * (-x2141 - 2 * x2143 + 2 * x2148 + 2 * x2149 + x49)
    x2151 = x0 * x56
    x2152 = x580 * x65
    x2153 = x17 * x598
    x2154 = x2136 - x2151 + x2152 - x2153
    x2155 = x0 * (x2154 + x59)
    x2156 = 2 * x2155
    x2157 = x1450 - x1457 + x869
    x2158 = x17 * x2157
    x2159 = x2142 * x65
    x2160 = x0 * (x2141 - x2156 - 2 * x2158 + 2 * x2159 + x62)
    x2161 = x1474 + x2150 - x2160
    x2162 = x0 * x77
    x2163 = x598 * x65
    x2164 = x17 * x615
    x2165 = x2151 - x2162 + x2163 - x2164
    x2166 = x0 * (x2165 + x80)
    x2167 = 2 * x2166
    x2168 = x1457 - x1465 + x881
    x2169 = x17 * x2168
    x2170 = x2157 * x65
    x2171 = x0 * (x2156 - x2167 - 2 * x2169 + 2 * x2170 + x83)
    x2172 = x1476 + x2160 - x2171
    x2173 = -x17 * x2172
    x2174 = x2161 * x65 + x2173
    x2175 = -x0 * x96 - x17 * x633 + x2162 + x615 * x65
    x2176 = x0 * (x2175 + x99)
    x2177 = x17 * (x1465 - x1478 + x900)
    x2178 = x2168 * x65
    x2179 = -x17 * (
        -x0 * (x102 + x2167 - 2 * x2176 - 2 * x2177 + 2 * x2178) + x1480 + x2171
    )
    x2180 = x2172 * x65 + x2179
    x2181 = x0 * x421
    x2182 = x0 * x425
    x2183 = x646 * x65
    x2184 = x17 * x652
    x2185 = x2181 - x2182 + x2183 - x2184
    x2186 = x0 * (x2185 + x428)
    x2187 = x0 * x433
    x2188 = x65 * x652
    x2189 = x17 * x665
    x2190 = x2182 - x2187 + x2188 - x2189
    x2191 = x0 * (x2190 + x436)
    x2192 = x1495 + x2186 - x2191
    x2193 = x0 * x448
    x2194 = x65 * x665
    x2195 = x17 * x677
    x2196 = x2187 - x2193 + x2194 - x2195
    x2197 = x0 * (x2196 + x451)
    x2198 = x1497 + x2191 - x2197
    x2199 = -x17 * x2198
    x2200 = x2192 * x65 + x2199
    x2201 = -x0 * x462 - x17 * x693 + x2193 + x65 * x677
    x2202 = -x17 * (-x0 * (x2201 + x465) + x1500 + x2197)
    x2203 = x2198 * x65 + x2202
    x2204 = x0 * x488
    x2205 = x0 * x492
    x2206 = x65 * x709
    x2207 = x17 * x714
    x2208 = x2204 - x2205 + x2206 - x2207
    x2209 = x0 * (x2208 + x495)
    x2210 = x0 * x500
    x2211 = x65 * x714
    x2212 = x17 * x726
    x2213 = x2205 - x2210 + x2211 - x2212
    x2214 = x0 * (x2213 + x503)
    x2215 = x1515 + x2209 - x2214
    x2216 = x0 * x515
    x2217 = x65 * x726
    x2218 = x17 * x737
    x2219 = x2210 - x2216 + x2217 - x2218
    x2220 = x0 * (x2219 + x518)
    x2221 = x1517 + x2214 - x2220
    x2222 = -x17 * x2221
    x2223 = x2215 * x65 + x2222
    x2224 = -x0 * x529 - x17 * x759 + x2216 + x65 * x737
    x2225 = -x17 * (-x0 * (x2224 + x532) + x1520 + x2220)
    x2226 = x2221 * x65 + x2225
    x2227 = x0 * x563
    x2228 = x0 * x570
    x2229 = x65 * x779
    x2230 = x17 * x786
    x2231 = x2227 - x2228 + x2229 - x2230
    x2232 = x0 * (x2231 + x573)
    x2233 = x0 * x581
    x2234 = x65 * x786
    x2235 = x17 * x800
    x2236 = x2228 - x2233 + x2234 - x2235
    x2237 = x0 * (x2236 + x584)
    x2238 = x1535 + x2232 - x2237
    x2239 = x0 * x599
    x2240 = x65 * x800
    x2241 = x17 * x813
    x2242 = x2233 - x2239 + x2240 - x2241
    x2243 = x0 * (x2242 + x602)
    x2244 = x1537 + x2237 - x2243
    x2245 = -x17 * x2244
    x2246 = x2238 * x65 + x2245
    x2247 = -x0 * x616 - x17 * x832 + x2239 + x65 * x813
    x2248 = -x17 * (-x0 * (x2247 + x619) + x1540 + x2243)
    x2249 = x2244 * x65 + x2248
    x2250 = x0 * x37
    x2251 = x0 * x46
    x2252 = x65 * x850
    x2253 = x17 * x856
    x2254 = x2250 - x2251 + x2252 - x2253
    x2255 = x0 * (x2254 + x49)
    x2256 = x0 * x59
    x2257 = x65 * x856
    x2258 = x17 * x869
    x2259 = x2251 - x2256 + x2257 - x2258
    x2260 = x0 * (x2259 + x62)
    x2261 = x1555 + x2255 - x2260
    x2262 = x0 * x80
    x2263 = x65 * x869
    x2264 = x17 * x881
    x2265 = x2256 - x2262 + x2263 - x2264
    x2266 = x0 * (x2265 + x83)
    x2267 = x1557 + x2260 - x2266
    x2268 = -x17 * x2267
    x2269 = x2261 * x65 + x2268
    x2270 = -x0 * x99 - x17 * x900 + x2262 + x65 * x881
    x2271 = -x17 * (-x0 * (x102 + x2270) + x1560 + x2266)
    x2272 = x2267 * x65 + x2271
    x2273 = x1572 * x6 + x1576
    x2274 = -x17 * (x1574 * x6 + x1581)
    x2275 = 3 * x1578
    x2276 = (
        x0 * (3 * x1573 - 3 * x1575 + 3 * x1577 - x2275)
        - x0 * (3 * x1579 - 3 * x1580 - 3 * x1582 + x2275)
        + x2273 * x65
        + x2274
    )
    x2277 = x1585 * x6 + x1589
    x2278 = -x17 * (x1587 * x6 + x1594)
    x2279 = 3 * x1591
    x2280 = x1597 * x6 + x1601
    x2281 = -x17 * (x1599 * x6 + x1606)
    x2282 = 3 * x1603
    x2283 = x1609 * x6 + x1613
    x2284 = -x17 * (x1611 * x6 + x1618)
    x2285 = 3 * x1615
    x2286 = x1621 * x6 + x1625
    x2287 = -x17 * (x1623 * x6 + x1630)
    x2288 = 3 * x1627
    x2289 = x1633 * x6 + x1637
    x2290 = -x17 * (x1635 * x6 + x1642)
    x2291 = 3 * x1639
    x2292 = x1645 * x6 + x1649
    x2293 = -x17 * (x1647 * x6 + x1654)
    x2294 = 3 * x1651
    x2295 = x1657 * x6 + x1661
    x2296 = -x17 * (x1659 * x6 + x1666)
    x2297 = 3 * x1663
    x2298 = x1669 * x6 + x1673
    x2299 = -x17 * (x1671 * x6 + x1678)
    x2300 = 3 * x1675
    x2301 = x1681 * x6 + x1685
    x2302 = -x17 * (x1683 * x6 + x1690)
    x2303 = 3 * x1687
    x2304 = x1693 * x6 + x1697
    x2305 = -x17 * (x1695 * x6 + x1702)
    x2306 = 3 * x1699
    x2307 = x1705 * x6 + x1711
    x2308 = -x17 * (x1706 * x6 + x1715)
    x2309 = 2 * x1708
    x2310 = x1716 * x6 + x1722
    x2311 = -x17 * (x1717 * x6 + x1726)
    x2312 = 2 * x1719
    x2313 = x1727 * x6 + x1733
    x2314 = -x17 * (x1728 * x6 + x1737)
    x2315 = 2 * x1730
    x2316 = x1738 * x6 + x1744
    x2317 = -x17 * (x1739 * x6 + x1748)
    x2318 = 2 * x1741
    x2319 = x1749 * x6 + x1755
    x2320 = -x17 * (x1750 * x6 + x1759)
    x2321 = 2 * x1752
    x2322 = x1760 * x6 + x1766
    x2323 = -x17 * (x1761 * x6 + x1770)
    x2324 = 2 * x1763
    x2325 = x1771 * x6 + x1777
    x2326 = -x17 * (x1772 * x6 + x1781)
    x2327 = 2 * x1774
    x2328 = x1782 * x6 + x1788
    x2329 = -x17 * (x1783 * x6 + x1792)
    x2330 = 2 * x1785
    x2331 = x1793 * x6 + x1799
    x2332 = -x17 * (x1794 * x6 + x1803)
    x2333 = 2 * x1796
    x2334 = x1804 * x6 + x1810
    x2335 = -x17 * (x1805 * x6 + x1814)
    x2336 = 2 * x1807
    x2337 = x1815 * x6 + x1821
    x2338 = -x17 * (x1816 * x6 + x1825)
    x2339 = 2 * x1818
    x2340 = x1826 * x6 + x1832
    x2341 = -x17 * (x1827 * x6 + x1836)
    x2342 = 2 * x1829
    x2343 = x1837 * x6 + x1843
    x2344 = -x17 * (x1838 * x6 + x1847)
    x2345 = 2 * x1840
    x2346 = x1848 * x6 + x1854
    x2347 = -x17 * (x1849 * x6 + x1858)
    x2348 = 2 * x1851
    x2349 = x1859 * x6 + x1865
    x2350 = -x17 * (x1860 * x6 + x1869)
    x2351 = 2 * x1862
    x2352 = x1883 * x6 + x1890
    x2353 = -x17 * (x1889 * x6 + x1893)
    x2354 = x1938 * x6 + x1955
    x2355 = -x17 * (x1954 * x6 + x1964)
    x2356 = x2014 * x6 + x2032
    x2357 = -x17 * (x2031 * x6 + x2041)
    x2358 = x2069 * x6 + x2081
    x2359 = -x17 * (x2080 * x6 + x2087)
    x2360 = x2115 * x6 + x2127
    x2361 = -x17 * (x2126 * x6 + x2133)
    x2362 = x2161 * x6 + x2173
    x2363 = -x17 * (x2172 * x6 + x2179)
    x2364 = x2192 * x6 + x2199
    x2365 = -x17 * (x2198 * x6 + x2202)
    x2366 = x2215 * x6 + x2222
    x2367 = -x17 * (x2221 * x6 + x2225)
    x2368 = x2238 * x6 + x2245
    x2369 = -x17 * (x2244 * x6 + x2248)
    x2370 = x2261 * x6 + x2268
    x2371 = -x17 * (x2267 * x6 + x2271)
    x2372 = x6 * x686 + x690
    x2373 = -x17 * (x6 * x689 + x700)
    x2374 = x6 * x748 + x752
    x2375 = -x17 * (x6 * x751 + x766)
    x2376 = x6 * x822 + x826
    x2377 = -x17 * (x6 * x825 + x839)
    x2378 = x6 * x890 + x894
    x2379 = -x17 * (x6 * x893 + x907)
    x2380 = x112 * x6 + x116
    x2381 = -x17 * (x115 * x6 + x132)
    x2382 = 2 * x578
    x2383 = 2 * x724
    x2384 = 2 * x567
    x2385 = 2 * x712
    x2386 = -2 * x516 + 2 * x725
    x2387 = x0 * (x1329 - x2382 - x2383 + x2384 + x2385 + x2386)
    x2388 = 3 * x2387
    x2389 = 2 * x854
    x2390 = 2 * x848
    x2391 = 2 * x560
    x2392 = 2 * x707
    x2393 = -2 * x501 + 2 * x713
    x2394 = x0 * (x1322 - x2384 - x2385 + x2391 + x2392 + x2393)
    x2395 = 3 * x2394
    x2396 = x0 * (x1452 - x2388 - x2389 + x2390 + x2395 + 3 * x795 - 3 * x802)
    x2397 = 4 * x2396
    x2398 = x1873 - x1880 + x873
    x2399 = x17 * x2398
    x2400 = 2 * x1235
    x2401 = -2 * x493 + 2 * x708
    x2402 = x0 * (x1317 - x2391 - x2392 + x2401 + 2 * x556 + 2 * x706)
    x2403 = x0 * (x1447 - x2390 - x2395 + 3 * x2402 + 3 * x790 - 3 * x791 + 2 * x846)
    x2404 = x65 * (-x1873 + x1876 + x863)
    x2405 = 2 * x596
    x2406 = 2 * x735
    x2407 = -2 * x530 + 2 * x736
    x2408 = x0 * (x1338 + x2382 + x2383 - x2405 - x2406 + x2407)
    x2409 = 3 * x2408
    x2410 = 2 * x867
    x2411 = x0 * (x1459 + x2388 + x2389 - x2409 - x2410 + 3 * x808 - 3 * x815)
    x2412 = 4 * x2411
    x2413 = x1880 - x1886 + x885
    x2414 = x17 * x2413
    x2415 = 2 * x1239
    x2416 = x2398 * x65
    x2417 = x0 * (x108 + x2397 + x2400 - x2412 - 4 * x2414 - x2415 + 4 * x2416)
    x2418 = (
        x0 * (x107 + 2 * x1233 - x2397 - 4 * x2399 - x2400 + 4 * x2403 + 4 * x2404)
        + x1891
        - x2417
    )
    x2419 = x0 * (x1348 + x2405 + x2406 - 2 * x547 - 2 * x613 - 2 * x757 + 2 * x758)
    x2420 = x0 * (x1467 + x2409 + x2410 - 3 * x2419 + 3 * x829 - 3 * x833 - 2 * x879)
    x2421 = x17 * (x1886 - x1892 + x903)
    x2422 = x2413 * x65
    x2423 = -x17 * (
        -x0 * (x119 - 2 * x1244 + x2412 + x2415 - 4 * x2420 - 4 * x2421 + 4 * x2422)
        + x1894
        + x2417
    )
    x2424 = x2418 * x65 + x2423
    x2425 = 2 * x1267
    x2426 = 2 * x1275
    x2427 = x0 * (-2 * x146 + 2 * x236)
    x2428 = x0 * (-2 * x153 + 2 * x243)
    x2429 = -x17 * x1897 + x1904 * x65 + x2427 - x2428
    x2430 = x0 * (x1266 + x2429)
    x2431 = 2 * x2430
    x2432 = x0 * (-2 * x164 + 2 * x254)
    x2433 = -x17 * x1925 + x1897 * x65 + x2428 - x2432
    x2434 = x0 * (x1274 + x2433)
    x2435 = 2 * x2434
    x2436 = -x1898 - x1901 + x1905 + x1908
    x2437 = x2436 * x65
    x2438 = x1898 - x1926 - x1929 + x1930
    x2439 = x17 * x2438
    x2440 = x0 * (
        2 * x1270 - 2 * x170 + x2425 - x2426 + x2431 - x2435 + 2 * x2437 - 2 * x2439
    )
    x2441 = 3 * x2440
    x2442 = x1909 - x1931 - x1934 + x1935
    x2443 = x17 * x2442
    x2444 = 2 * x1277
    x2445 = x0 * (-2 * x142 + 2 * x232) - x17 * x1904 + x1916 * x65 - x2427
    x2446 = x0 * (x1261 + x2445)
    x2447 = x65 * (-x1905 - x1913 + x1917 + x1918)
    x2448 = x17 * x2436
    x2449 = x0 * (
        2 * x1256
        + 2 * x1262
        - 2 * x159
        - x2425
        - x2431
        + 2 * x2446
        + 2 * x2447
        - 2 * x2448
    )
    x2450 = x65 * (-x1909 - x1912 + x1919 + x1920)
    x2451 = 2 * x1285
    x2452 = x0 * (-2 * x182 + 2 * x272)
    x2453 = -x17 * x1941 + x1925 * x65 + x2432 - x2452
    x2454 = x0 * (x1284 + x2453)
    x2455 = 2 * x2454
    x2456 = x2438 * x65
    x2457 = x1926 - x1942 - x1945 + x1946
    x2458 = x17 * x2457
    x2459 = x0 * (
        2 * x1280 - 2 * x188 + x2426 + x2435 - x2451 - x2455 + 2 * x2456 - 2 * x2458
    )
    x2460 = 3 * x2459
    x2461 = x1931 - x1947 - x1950 + x1951
    x2462 = x17 * x2461
    x2463 = 2 * x1287
    x2464 = x2442 * x65
    x2465 = x0 * (x213 + x2441 + x2444 - x2460 - 3 * x2462 - x2463 + 3 * x2464)
    x2466 = (
        x0 * (2 * x1269 + x212 - x2441 - 3 * x2443 - x2444 + 3 * x2449 + 3 * x2450)
        + x1956
        - x2465
    )
    x2467 = -x0 * (-2 * x199 + 2 * x289) - x17 * x1957 + x1941 * x65 + x2452
    x2468 = x0 * (x1295 + x2467)
    x2469 = x2457 * x65
    x2470 = x17 * (x1942 - x1958 - x1959 + x1960)
    x2471 = x0 * (
        2 * x1291
        - 2 * x1296
        - 2 * x205
        + x2451
        + x2455
        - 2 * x2468
        + 2 * x2469
        - 2 * x2470
    )
    x2472 = x17 * (x1947 - x1961 - x1962 + x1963)
    x2473 = x2461 * x65
    x2474 = -x17 * (
        -x0 * (-2 * x1298 + x218 + x2460 + x2463 - 3 * x2471 - 3 * x2472 + 3 * x2473)
        + x1965
        + x2465
    )
    x2475 = x2466 * x65 + x2474
    x2476 = 2 * x1331
    x2477 = x0 * x565
    x2478 = x0 * x576
    x2479 = -x17 * x2001 + x1970 * x65 + x2477 - x2478
    x2480 = x0 * (x1330 + x2479)
    x2481 = 2 * x2480
    x2482 = x1971 - x2002 - x2005 + x2006
    x2483 = x17 * x2482
    x2484 = 2 * x1324
    x2485 = x0 * x558
    x2486 = -x17 * x1970 + x1978 * x65 - x2477 + x2485
    x2487 = x0 * (x1323 + x2486)
    x2488 = 2 * x2487
    x2489 = -x1971 - x1974 + x1979 + x1982
    x2490 = x2489 * x65
    x2491 = x0 * (x1452 - x2476 - x2481 - 2 * x2483 + x2484 + x2488 + 2 * x2490)
    x2492 = 3 * x2491
    x2493 = x1983 - x2007 - x2010 + x2011
    x2494 = x17 * x2493
    x2495 = 2 * x1333
    x2496 = x17 * x2489
    x2497 = x0 * x554 - x17 * x1978 + x1991 * x65 - x2485
    x2498 = x0 * (x1318 + x2497)
    x2499 = x65 * (-x1979 - x1987 + x1992 + x1993)
    x2500 = x0 * (2 * x1319 + x1447 - x2484 - x2488 - 2 * x2496 + 2 * x2498 + 2 * x2499)
    x2501 = x65 * (-x1983 - x1986 + x1994 + x1995)
    x2502 = 2 * x1340
    x2503 = x0 * x594
    x2504 = -x17 * x2018 + x2001 * x65 + x2478 - x2503
    x2505 = x0 * (x1339 + x2504)
    x2506 = 2 * x2505
    x2507 = x2002 - x2019 - x2022 + x2023
    x2508 = x17 * x2507
    x2509 = x2482 * x65
    x2510 = x0 * (x1459 + x2476 + x2481 - x2502 - x2506 - 2 * x2508 + 2 * x2509)
    x2511 = 3 * x2510
    x2512 = x2007 - x2024 - x2027 + x2028
    x2513 = x17 * x2512
    x2514 = 2 * x1342
    x2515 = x2493 * x65
    x2516 = x0 * (x108 + x2492 + x2495 - x2511 - 3 * x2513 - x2514 + 3 * x2515)
    x2517 = (
        x0 * (x107 + 2 * x1326 - x2492 - 3 * x2494 - x2495 + 3 * x2500 + 3 * x2501)
        + x2033
        - x2516
    )
    x2518 = -x0 * x611 - x17 * x2034 + x2018 * x65 + x2503
    x2519 = x0 * (x1349 + x2518)
    x2520 = x17 * (x2019 - x2035 - x2036 + x2037)
    x2521 = x2507 * x65
    x2522 = x0 * (-2 * x1350 + x1467 + x2502 + x2506 - 2 * x2519 - 2 * x2520 + 2 * x2521)
    x2523 = x17 * (x2024 - x2038 - x2039 + x2040)
    x2524 = x2512 * x65
    x2525 = -x17 * (
        -x0 * (x119 - 2 * x1352 + x2511 + x2514 - 3 * x2522 - 3 * x2523 + 3 * x2524)
        + x2042
        + x2516
    )
    x2526 = x2517 * x65 + x2525
    x2527 = 2 * x1376
    x2528 = x0 * (-2 * x246 + 2 * x423)
    x2529 = x0 * (-2 * x257 + 2 * x431)
    x2530 = -x17 * x2062 + x2047 * x65 + x2528 - x2529
    x2531 = x0 * (x1375 + x2530)
    x2532 = 2 * x2531
    x2533 = x2048 - x2063 - x2066 + x2067
    x2534 = x17 * x2533
    x2535 = x0 * (-2 * x239 + 2 * x419) - x17 * x2047 + x2055 * x65 - x2528
    x2536 = x0 * (x1370 + x2535)
    x2537 = x65 * (-x2048 - x2051 + x2056 + x2057)
    x2538 = 2 * x1383
    x2539 = x0 * (-2 * x275 + 2 * x446)
    x2540 = -x17 * x2073 + x2062 * x65 + x2529 - x2539
    x2541 = x0 * (x1382 + x2540)
    x2542 = 2 * x2541
    x2543 = x2063 - x2074 - x2077 + x2078
    x2544 = x17 * x2543
    x2545 = x2533 * x65
    x2546 = x0 * (x2527 + x2532 - x2538 - x2542 - 2 * x2544 + 2 * x2545 + x303)
    x2547 = (
        x0 * (2 * x1371 - x2527 - x2532 - 2 * x2534 + 2 * x2536 + 2 * x2537 + x302)
        + x2082
        - x2546
    )
    x2548 = -x0 * (-2 * x292 + 2 * x460) - x17 * x2083 + x2073 * x65 + x2539
    x2549 = x0 * (x1390 + x2548)
    x2550 = x17 * (x2074 - x2084 - x2085 + x2086)
    x2551 = x2543 * x65
    x2552 = -x17 * (
        -x0 * (-2 * x1391 + x2538 + x2542 - 2 * x2549 - 2 * x2550 + 2 * x2551 + x308)
        + x2088
        + x2546
    )
    x2553 = x2547 * x65 + x2552
    x2554 = 2 * x1415
    x2555 = x0 * x776
    x2556 = x0 * x783
    x2557 = -x17 * x2108 + x2093 * x65 + x2555 - x2556
    x2558 = x0 * (x1414 + x2557)
    x2559 = 2 * x2558
    x2560 = x2094 - x2109 - x2112 + x2113
    x2561 = x17 * x2560
    x2562 = x0 * x773 - x17 * x2093 + x2101 * x65 - x2555
    x2563 = x0 * (x1409 + x2562)
    x2564 = x65 * (-x2094 - x2097 + x2102 + x2103)
    x2565 = 2 * x1422
    x2566 = x0 * x797
    x2567 = -x17 * x2119 + x2108 * x65 + x2556 - x2566
    x2568 = x0 * (x1421 + x2567)
    x2569 = 2 * x2568
    x2570 = x2109 - x2120 - x2123 + x2124
    x2571 = x17 * x2570
    x2572 = x2560 * x65
    x2573 = x0 * (x2554 + x2559 - x2565 - x2569 - 2 * x2571 + 2 * x2572 + x397)
    x2574 = (
        x0 * (2 * x1410 - x2554 - x2559 - 2 * x2561 + 2 * x2563 + 2 * x2564 + x396)
        + x2128
        - x2573
    )
    x2575 = -x0 * x810 - x17 * x2129 + x2119 * x65 + x2566
    x2576 = x0 * (x1429 + x2575)
    x2577 = x17 * (x2120 - x2130 - x2131 + x2132)
    x2578 = x2570 * x65
    x2579 = -x17 * (
        -x0 * (-2 * x1430 + x2565 + x2569 - 2 * x2576 - 2 * x2577 + 2 * x2578 + x402)
        + x2134
        + x2573
    )
    x2580 = x2574 * x65 + x2579
    x2581 = 2 * x1454
    x2582 = x0 * x1322
    x2583 = x0 * x1329
    x2584 = -x17 * x2154 + x2139 * x65 + x2582 - x2583
    x2585 = x0 * (x1453 + x2584)
    x2586 = 2 * x2585
    x2587 = x2140 - x2155 - x2158 + x2159
    x2588 = x17 * x2587
    x2589 = x0 * x1317 - x17 * x2139 + x2147 * x65 - x2582
    x2590 = x0 * (x1448 + x2589)
    x2591 = x65 * (-x2140 - x2143 + x2148 + x2149)
    x2592 = 2 * x1461
    x2593 = x0 * x1338
    x2594 = -x17 * x2165 + x2154 * x65 + x2583 - x2593
    x2595 = x0 * (x1460 + x2594)
    x2596 = 2 * x2595
    x2597 = x2155 - x2166 - x2169 + x2170
    x2598 = x17 * x2597
    x2599 = x2587 * x65
    x2600 = x0 * (x108 + x2581 + x2586 - x2592 - x2596 - 2 * x2598 + 2 * x2599)
    x2601 = (
        x0 * (x107 + 2 * x1449 - x2581 - x2586 - 2 * x2588 + 2 * x2590 + 2 * x2591)
        + x2174
        - x2600
    )
    x2602 = -x0 * x1348 - x17 * x2175 + x2165 * x65 + x2593
    x2603 = x0 * (x1468 + x2602)
    x2604 = x17 * (x2166 - x2176 - x2177 + x2178)
    x2605 = x2597 * x65
    x2606 = -x17 * (
        -x0 * (x119 - 2 * x1469 + x2592 + x2596 - 2 * x2603 - 2 * x2604 + 2 * x2605)
        + x2180
        + x2600
    )
    x2607 = x2601 * x65 + x2606
    x2608 = 2 * x1483
    x2609 = x0 * (-2 * x434 + 2 * x651)
    x2610 = x0 * (-2 * x426 + 2 * x645) - x17 * x2190 + x2185 * x65 - x2609
    x2611 = 2 * x1486
    x2612 = x0 * (-2 * x449 + 2 * x664)
    x2613 = -x17 * x2196 + x2190 * x65 + x2609 - x2612
    x2614 = x0 * (x2608 - x2611 + x2613 + x471)
    x2615 = x0 * (2 * x1482 - x2608 + x2610 + x470) + x2200 - x2614
    x2616 = -x0 * (-2 * x463 + 2 * x676) - x17 * x2201 + x2196 * x65 + x2612
    x2617 = -x17 * (-x0 * (-2 * x1490 + x2611 + x2616 + x476) + x2203 + x2614)
    x2618 = x2615 * x65 + x2617
    x2619 = 2 * x1503
    x2620 = x0 * x2393
    x2621 = x0 * x2401 - x17 * x2213 + x2208 * x65 - x2620
    x2622 = 2 * x1506
    x2623 = x0 * x2386
    x2624 = -x17 * x2219 + x2213 * x65 + x2620 - x2623
    x2625 = x0 * (x2619 - x2622 + x2624 + x538)
    x2626 = x0 * (2 * x1502 - x2619 + x2621 + x537) + x2223 - x2625
    x2627 = -x0 * x2407 - x17 * x2224 + x2219 * x65 + x2623
    x2628 = -x17 * (-x0 * (-2 * x1510 + x2622 + x2627 + x543) + x2226 + x2625)
    x2629 = x2626 * x65 + x2628
    x2630 = 2 * x1523
    x2631 = x0 * (-2 * x582 + 2 * x785)
    x2632 = x0 * (-2 * x571 + 2 * x778) - x17 * x2236 + x2231 * x65 - x2631
    x2633 = 2 * x1526
    x2634 = x0 * (-2 * x600 + 2 * x799)
    x2635 = -x17 * x2242 + x2236 * x65 + x2631 - x2634
    x2636 = x0 * (x2630 - x2633 + x2635 + x625)
    x2637 = x0 * (2 * x1522 - x2630 + x2632 + x624) + x2246 - x2636
    x2638 = -x0 * (-2 * x617 + 2 * x812) - x17 * x2247 + x2242 * x65 + x2634
    x2639 = -x17 * (-x0 * (-2 * x1530 + x2633 + x2638 + x630) + x2249 + x2636)
    x2640 = x2637 * x65 + x2639
    x2641 = 2 * x1543
    x2642 = x0 * x1452
    x2643 = x0 * x1447 - x17 * x2259 + x2254 * x65 - x2642
    x2644 = 2 * x1546
    x2645 = x0 * x1459
    x2646 = -x17 * x2265 + x2259 * x65 + x2642 - x2645
    x2647 = x0 * (x108 + x2641 - x2644 + x2646)
    x2648 = x0 * (x107 + 2 * x1542 - x2641 + x2643) + x2269 - x2647
    x2649 = -x0 * x1467 - x17 * x2270 + x2265 * x65 + x2645
    x2650 = -x17 * (-x0 * (x119 - 2 * x1550 + x2644 + x2649) + x2272 + x2647)
    x2651 = x2648 * x65 + x2650
    x2652 = x2273 * x6 + x2274
    x2653 = 3 * x1882
    x2654 = x2396 - x2411 - x2414 + x2416
    x2655 = 3 * x1873
    x2656 = 3 * x489
    x2657 = 3 * x485
    x2658 = x229 * (-x2656 + x2657 - 3 * x337 + 3 * x490 + x653)
    x2659 = 3 * x0
    x2660 = 3 * x497
    x2661 = x229 * (x2656 - x2660 - 3 * x349 + 3 * x498 + x666)
    x2662 = x2659 * (x1872 + x2658 - x2661 + x744)
    x2663 = 3 * x65
    x2664 = -x2387 + x2394 + x804
    x2665 = 3 * x17
    x2666 = 4 * x0
    x2667 = 3 * x1880
    x2668 = 3 * x512
    x2669 = x229 * (x2660 - x2668 - 3 * x368 + 3 * x513 + x678)
    x2670 = x2659 * (x1879 + x2661 - x2669 + x755)
    x2671 = x2387 - x2408 + x817
    x2672 = x2666 * (
        x2655
        + x2662
        + x2663 * x2664
        - x2665 * x2671
        - x2667
        - x2670
        + 3 * x865
        - 3 * x871
    )
    x2673 = 3 * x1937
    x2674 = 3 * x321
    x2675 = 3 * x326
    x2676 = x0 * (3 * x1902 - 3 * x1903 + x2674 - x2675)
    x2677 = 3 * x334
    x2678 = x0 * (3 * x1895 - 3 * x1896 + x2675 - x2677)
    x2679 = 3 * x1905
    x2680 = 3 * x1898
    x2681 = x229 * (
        -x17 * x2433 - 3 * x1901 + 3 * x1908 + x2429 * x65 + x2676 - x2678 + x2679 - x2680
    )
    x2682 = x2430 - x2434 + x2437 - x2439
    x2683 = 3 * x346
    x2684 = x0 * (3 * x1923 - 3 * x1924 + x2677 - x2683)
    x2685 = 3 * x1926
    x2686 = x229 * (
        -x17 * x2453 - 3 * x1929 + 3 * x1930 + x2433 * x65 + x2678 + x2680 - x2684 - x2685
    )
    x2687 = x2434 - x2454 + x2456 - x2458
    x2688 = x2659 * (x138 * x2682 - x18 * x2687 + x1936 + x2681 - x2686)
    x2689 = x2440 - x2459 - x2462 + x2464
    x2690 = 3 * x2013
    x2691 = 3 * x1975
    x2692 = 3 * x1966
    x2693 = x0 * (3 * x1976 - 3 * x1977 + x2691 - x2692)
    x2694 = 3 * x1967
    x2695 = x0 * (3 * x1968 - 3 * x1969 + x2692 - x2694)
    x2696 = 3 * x1979
    x2697 = 3 * x1971
    x2698 = x229 * (
        -x17 * x2479 - 3 * x1974 + 3 * x1982 + x2486 * x65 + x2693 - x2695 + x2696 - x2697
    )
    x2699 = -x2480 - x2483 + x2487 + x2490
    x2700 = 3 * x1998
    x2701 = x0 * (3 * x1999 - 3 * x2000 + x2694 - x2700)
    x2702 = 3 * x2002
    x2703 = x229 * (
        -x17 * x2504 - 3 * x2005 + 3 * x2006 + x2479 * x65 + x2695 + x2697 - x2701 - x2702
    )
    x2704 = x2480 - x2505 - x2508 + x2509
    x2705 = x2659 * (x138 * x2699 - x18 * x2704 + x2012 + x2698 - x2703)
    x2706 = x2491 - x2510 - x2513 + x2515
    x2707 = 3 * x2068
    x2708 = 3 * x2043
    x2709 = 3 * x2044
    x2710 = x0 * (3 * x2045 - 3 * x2046 + x2708 - x2709)
    x2711 = 3 * x2048
    x2712 = 3 * x2059
    x2713 = x0 * (3 * x2060 - 3 * x2061 + x2709 - x2712)
    x2714 = 3 * x2063
    x2715 = x229 * (
        -x17 * x2540 - 3 * x2066 + 3 * x2067 + x2530 * x65 + x2710 + x2711 - x2713 - x2714
    )
    x2716 = x2531 - x2541 - x2544 + x2545
    x2717 = 3 * x2114
    x2718 = 3 * x2089
    x2719 = 3 * x2090
    x2720 = x0 * (3 * x2091 - 3 * x2092 + x2718 - x2719)
    x2721 = 3 * x2094
    x2722 = 3 * x2105
    x2723 = x0 * (3 * x2106 - 3 * x2107 + x2719 - x2722)
    x2724 = 3 * x2109
    x2725 = x229 * (
        -x17 * x2567 - 3 * x2112 + 3 * x2113 + x2557 * x65 + x2720 + x2721 - x2723 - x2724
    )
    x2726 = x2558 - x2568 - x2571 + x2572
    x2727 = 3 * x2160
    x2728 = 3 * x2135
    x2729 = 3 * x2136
    x2730 = x0 * (3 * x2137 - 3 * x2138 + x2728 - x2729)
    x2731 = 3 * x2140
    x2732 = 3 * x2151
    x2733 = x0 * (3 * x2152 - 3 * x2153 + x2729 - x2732)
    x2734 = 3 * x2155
    x2735 = x229 * (
        -x17 * x2594 - 3 * x2158 + 3 * x2159 + x2584 * x65 + x2730 + x2731 - x2733 - x2734
    )
    x2736 = x2585 - x2595 - x2598 + x2599
    x2737 = 3 * x2182
    x2738 = 3 * x2187
    x2739 = x0 * (3 * x2188 - 3 * x2189 + x2737 - x2738)
    x2740 = 3 * x2191
    x2741 = 3 * x2205
    x2742 = 3 * x2210
    x2743 = x0 * (3 * x2211 - 3 * x2212 + x2741 - x2742)
    x2744 = 3 * x2214
    x2745 = 3 * x2228
    x2746 = 3 * x2233
    x2747 = x0 * (3 * x2234 - 3 * x2235 + x2745 - x2746)
    x2748 = 3 * x2237
    x2749 = 3 * x2251
    x2750 = 3 * x2256
    x2751 = x0 * (3 * x2257 - 3 * x2258 + x2749 - x2750)
    x2752 = 3 * x2260

    # 225 item(s)
    S = numpy.array(
        [
            x137,
            x0 * (3 * x162 - x174 + 3 * x179 - 3 * x194)
            - x0 * (x174 - 3 * x195 + 3 * x196 - 3 * x211)
            - x17 * x228
            + x217 * x65,
            x137,
            x0 * (3 * x252 - x264 + 3 * x269 - 3 * x284)
            - x0 * (x264 - 3 * x285 + 3 * x286 - 3 * x301)
            - x17 * x317
            + x307 * x65,
            x0 * (3 * x343 - x356 + 3 * x361 - 3 * x377)
            - x0 * (x356 - 3 * x378 + 3 * x379 - 3 * x395)
            - x17 * x413
            + x401 * x65,
            x137,
            x0 * (3 * x429 - x438 + 3 * x443 - 3 * x455)
            - x0 * (x438 - 3 * x456 + 3 * x457 - 3 * x469)
            - x17 * x483
            + x475 * x65,
            x0 * (3 * x496 - x505 + 3 * x510 - 3 * x522)
            - x0 * (x505 - 3 * x523 + 3 * x524 - 3 * x536)
            - x17 * x552
            + x542 * x65,
            x0 * (3 * x574 - x586 + 3 * x591 - 3 * x606)
            - x0 * (x586 - 3 * x607 + 3 * x608 - 3 * x623)
            - x17 * x639
            + x629 * x65,
            x137,
            x705,
            x771,
            x844,
            x912,
            x137,
            x943,
            x0 * (x217 + x951)
            - x0 * (x228 + x957)
            - x17 * (x959 - x965 + x966 - x971)
            + x65 * (x958 - x959 + x961 - x964),
            x943,
            x0 * (x307 + x979)
            - x0 * (x317 + x985)
            - x17 * (x987 - x993 + x994 - x999)
            + x65 * (x986 - x987 + x989 - x992),
            x0 * (x1007 + x401)
            - x0 * (x1013 + x413)
            - x17 * (x1015 - x1021 + x1022 - x1027)
            + x65 * (x1014 - x1015 + x1017 - x1020),
            x943,
            x0 * (x1035 + x475)
            - x0 * (x1041 + x483)
            - x17 * (x1043 - x1049 + x1050 - x1055)
            + x65 * (x1042 - x1043 + x1045 - x1048),
            x0 * (x1063 + x542)
            - x0 * (x1069 + x552)
            - x17 * (x1071 - x1077 + x1078 - x1083)
            + x65 * (x1070 - x1071 + x1073 - x1076),
            x0 * (x1091 + x629)
            - x0 * (x1097 + x639)
            - x17 * (x1099 - x1105 + x1106 - x1111)
            + x65 * (x1098 - x1099 + x1101 - x1104),
            x943,
            x0 * (x1119 + x692)
            - x0 * (x1125 + x702)
            - x17 * (x1127 - x1135 + x1136 - x1141)
            + x65 * (x1126 - x1127 + x1130 - x1134),
            x0 * (x1149 + x754)
            - x0 * (x1155 + x768)
            - x17 * (x1157 - x1165 + x1166 - x1171)
            + x65 * (x1156 - x1157 + x1160 - x1164),
            x0 * (x1179 + x828)
            - x0 * (x1185 + x841)
            - x17 * (x1187 - x1195 + x1196 - x1201)
            + x65 * (x1186 - x1187 + x1190 - x1194),
            x0 * (x1209 + x896)
            - x0 * (x1215 + x909)
            - x17 * (x1217 - x1225 + x1226 - x1231)
            + x65 * (x1216 - x1217 + x1220 - x1224),
            x943,
            x0 * (2 * x1237 - 2 * x1241)
            - x0 * (2 * x1242 - 2 * x1246)
            + x1252 * x65
            - x1255 * x17,
            x0 * (2 * x1279 - 2 * x1289)
            - x0 * (2 * x1290 - 2 * x1300)
            + x1306 * x65
            - x1313 * x17,
            x0 * (2 * x1335 - 2 * x1344)
            - x0 * (2 * x1345 - 2 * x1354)
            + x1360 * x65
            - x1365 * x17,
            x0 * (2 * x1378 - 2 * x1385)
            - x0 * (2 * x1386 - 2 * x1393)
            + x1399 * x65
            - x1404 * x17,
            x0 * (2 * x1417 - 2 * x1424)
            - x0 * (2 * x1425 - 2 * x1432)
            + x1438 * x65
            - x1443 * x17,
            x0 * (2 * x1456 - 2 * x1463)
            - x0 * (2 * x1464 - 2 * x1471)
            + x1477 * x65
            - x1481 * x17,
            x0 * (2 * x1485 - 2 * x1488)
            - x0 * (2 * x1489 - 2 * x1492)
            + x1498 * x65
            - x1501 * x17,
            x0 * (2 * x1505 - 2 * x1508)
            - x0 * (2 * x1509 - 2 * x1512)
            + x1518 * x65
            - x1521 * x17,
            x0 * (2 * x1525 - 2 * x1528)
            - x0 * (2 * x1529 - 2 * x1532)
            + x1538 * x65
            - x1541 * x17,
            x0 * (2 * x1545 - 2 * x1548)
            - x0 * (2 * x1549 - 2 * x1552)
            + x1558 * x65
            - x1561 * x17,
            x0 * (2 * x660 - 2 * x670)
            - x0 * (2 * x672 - 2 * x682)
            + x1562 * x65
            - x1563 * x17,
            x0 * (2 * x722 - 2 * x731)
            - x0 * (2 * x733 - 2 * x742)
            + x1564 * x65
            - x1565 * x17,
            x0 * (2 * x794 - 2 * x805)
            - x0 * (2 * x807 - 2 * x818)
            + x1566 * x65
            - x1567 * x17,
            x0 * (2 * x864 - 2 * x874)
            - x0 * (2 * x876 - 2 * x886)
            + x1568 * x65
            - x1569 * x17,
            -x0 * (-2 * x106 + 2 * x89)
            + x0 * (2 * x70 - 2 * x87)
            + x1570 * x65
            - x1571 * x17,
            x1584,
            x0 * (x1585 - x1596 + 2 * x958 + 2 * x961 - 2 * x964)
            - x0 * (x1587 + x1596 - 2 * x965 + 2 * x966 - 2 * x971)
            - x17 * (x1591 + x1592 + x1594 - x1595)
            + x65 * (x1586 + x1589 + x1590 - x1591),
            x1584,
            x0 * (x1597 - x1608 + 2 * x986 + 2 * x989 - 2 * x992)
            - x0 * (x1599 + x1608 - 2 * x993 + 2 * x994 - 2 * x999)
            - x17 * (x1603 + x1604 + x1606 - x1607)
            + x65 * (x1598 + x1601 + x1602 - x1603),
            x0 * (2 * x1014 + 2 * x1017 - 2 * x1020 + x1609 - x1620)
            - x0 * (-2 * x1021 + 2 * x1022 - 2 * x1027 + x1611 + x1620)
            - x17 * (x1615 + x1616 + x1618 - x1619)
            + x65 * (x1610 + x1613 + x1614 - x1615),
            x1584,
            x0 * (2 * x1042 + 2 * x1045 - 2 * x1048 + x1621 - x1632)
            - x0 * (-2 * x1049 + 2 * x1050 - 2 * x1055 + x1623 + x1632)
            - x17 * (x1627 + x1628 + x1630 - x1631)
            + x65 * (x1622 + x1625 + x1626 - x1627),
            x0 * (2 * x1070 + 2 * x1073 - 2 * x1076 + x1633 - x1644)
            - x0 * (-2 * x1077 + 2 * x1078 - 2 * x1083 + x1635 + x1644)
            - x17 * (x1639 + x1640 + x1642 - x1643)
            + x65 * (x1634 + x1637 + x1638 - x1639),
            x0 * (2 * x1098 + 2 * x1101 - 2 * x1104 + x1645 - x1656)
            - x0 * (-2 * x1105 + 2 * x1106 - 2 * x1111 + x1647 + x1656)
            - x17 * (x1651 + x1652 + x1654 - x1655)
            + x65 * (x1646 + x1649 + x1650 - x1651),
            x1584,
            x0 * (2 * x1126 + 2 * x1130 - 2 * x1134 + x1657 - x1668)
            - x0 * (-2 * x1135 + 2 * x1136 - 2 * x1141 + x1659 + x1668)
            - x17 * (x1663 + x1664 + x1666 - x1667)
            + x65 * (x1658 + x1661 + x1662 - x1663),
            x0 * (2 * x1156 + 2 * x1160 - 2 * x1164 + x1669 - x1680)
            - x0 * (-2 * x1165 + 2 * x1166 - 2 * x1171 + x1671 + x1680)
            - x17 * (x1675 + x1676 + x1678 - x1679)
            + x65 * (x1670 + x1673 + x1674 - x1675),
            x0 * (2 * x1186 + 2 * x1190 - 2 * x1194 + x1681 - x1692)
            - x0 * (-2 * x1195 + 2 * x1196 - 2 * x1201 + x1683 + x1692)
            - x17 * (x1687 + x1688 + x1690 - x1691)
            + x65 * (x1682 + x1685 + x1686 - x1687),
            x0 * (2 * x1216 + 2 * x1220 - 2 * x1224 + x1693 - x1704)
            - x0 * (-2 * x1225 + 2 * x1226 - 2 * x1231 + x1695 + x1704)
            - x17 * (x1699 + x1700 + x1702 - x1703)
            + x65 * (x1694 + x1697 + x1698 - x1699),
            x1584,
            x0 * (x1252 + x1705)
            - x0 * (x1255 + x1706)
            - x17 * (x1708 - x1712 + x1713 + x1715)
            + x65 * (x1707 - x1708 + x1709 + x1711),
            x0 * (x1306 + x1716)
            - x0 * (x1313 + x1717)
            - x17 * (x1719 - x1723 + x1724 + x1726)
            + x65 * (x1718 - x1719 + x1720 + x1722),
            x0 * (x1360 + x1727)
            - x0 * (x1365 + x1728)
            - x17 * (x1730 - x1734 + x1735 + x1737)
            + x65 * (x1729 - x1730 + x1731 + x1733),
            x0 * (x1399 + x1738)
            - x0 * (x1404 + x1739)
            - x17 * (x1741 - x1745 + x1746 + x1748)
            + x65 * (x1740 - x1741 + x1742 + x1744),
            x0 * (x1438 + x1749)
            - x0 * (x1443 + x1750)
            - x17 * (x1752 - x1756 + x1757 + x1759)
            + x65 * (x1751 - x1752 + x1753 + x1755),
            x0 * (x1477 + x1760)
            - x0 * (x1481 + x1761)
            - x17 * (x1763 - x1767 + x1768 + x1770)
            + x65 * (x1762 - x1763 + x1764 + x1766),
            x0 * (x1498 + x1771)
            - x0 * (x1501 + x1772)
            - x17 * (x1774 - x1778 + x1779 + x1781)
            + x65 * (x1773 - x1774 + x1775 + x1777),
            x0 * (x1518 + x1782)
            - x0 * (x1521 + x1783)
            - x17 * (x1785 - x1789 + x1790 + x1792)
            + x65 * (x1784 - x1785 + x1786 + x1788),
            x0 * (x1538 + x1793)
            - x0 * (x1541 + x1794)
            - x17 * (x1796 - x1800 + x1801 + x1803)
            + x65 * (x1795 - x1796 + x1797 + x1799),
            x0 * (x1558 + x1804)
            - x0 * (x1561 + x1805)
            - x17 * (x1807 - x1811 + x1812 + x1814)
            + x65 * (x1806 - x1807 + x1808 + x1810),
            x0 * (x1562 + x1815)
            - x0 * (x1563 + x1816)
            - x17 * (x1818 - x1822 + x1823 + x1825)
            + x65 * (x1817 - x1818 + x1819 + x1821),
            x0 * (x1564 + x1826)
            - x0 * (x1565 + x1827)
            - x17 * (x1829 - x1833 + x1834 + x1836)
            + x65 * (x1828 - x1829 + x1830 + x1832),
            x0 * (x1566 + x1837)
            - x0 * (x1567 + x1838)
            - x17 * (x1840 - x1844 + x1845 + x1847)
            + x65 * (x1839 - x1840 + x1841 + x1843),
            x0 * (x1568 + x1848)
            - x0 * (x1569 + x1849)
            - x17 * (x1851 - x1855 + x1856 + x1858)
            + x65 * (x1850 - x1851 + x1852 + x1854),
            x0 * (x1570 + x1859)
            - x0 * (x1571 + x1860)
            - x17 * (x1862 - x1866 + x1867 + x1869)
            + x65 * (x1861 - x1862 + x1863 + x1865),
            x0 * x1883 - x0 * x1889 - x17 * x1894 + x1891 * x65,
            x0 * x1938 - x0 * x1954 - x17 * x1965 + x1956 * x65,
            x0 * x2014 - x0 * x2031 - x17 * x2042 + x2033 * x65,
            x0 * x2069 - x0 * x2080 - x17 * x2088 + x2082 * x65,
            x0 * x2115 - x0 * x2126 - x17 * x2134 + x2128 * x65,
            x0 * x2161 - x0 * x2172 - x17 * x2180 + x2174 * x65,
            x0 * x2192 - x0 * x2198 - x17 * x2203 + x2200 * x65,
            x0 * x2215 - x0 * x2221 - x17 * x2226 + x2223 * x65,
            x0 * x2238 - x0 * x2244 - x17 * x2249 + x2246 * x65,
            x0 * x2261 - x0 * x2267 - x17 * x2272 + x2269 * x65,
            x0 * x686 - x0 * x689 - x17 * x701 + x65 * x691,
            x0 * x748 - x0 * x751 - x17 * x767 + x65 * x753,
            x0 * x822 - x0 * x825 - x17 * x840 + x65 * x827,
            x0 * x890 - x0 * x893 - x17 * x908 + x65 * x895,
            x0 * x112 - x0 * x115 + x117 * x65 - x133 * x17,
            x2276,
            x0 * (3 * x1586 - 3 * x1588 + 3 * x1590 - x2279)
            - x0 * (3 * x1592 - 3 * x1593 - 3 * x1595 + x2279)
            + x2277 * x65
            + x2278,
            x2276,
            x0 * (3 * x1598 - 3 * x1600 + 3 * x1602 - x2282)
            - x0 * (3 * x1604 - 3 * x1605 - 3 * x1607 + x2282)
            + x2280 * x65
            + x2281,
            x0 * (3 * x1610 - 3 * x1612 + 3 * x1614 - x2285)
            - x0 * (3 * x1616 - 3 * x1617 - 3 * x1619 + x2285)
            + x2283 * x65
            + x2284,
            x2276,
            x0 * (3 * x1622 - 3 * x1624 + 3 * x1626 - x2288)
            - x0 * (3 * x1628 - 3 * x1629 - 3 * x1631 + x2288)
            + x2286 * x65
            + x2287,
            x0 * (3 * x1634 - 3 * x1636 + 3 * x1638 - x2291)
            - x0 * (3 * x1640 - 3 * x1641 - 3 * x1643 + x2291)
            + x2289 * x65
            + x2290,
            x0 * (3 * x1646 - 3 * x1648 + 3 * x1650 - x2294)
            - x0 * (3 * x1652 - 3 * x1653 - 3 * x1655 + x2294)
            + x2292 * x65
            + x2293,
            x2276,
            x0 * (3 * x1658 - 3 * x1660 + 3 * x1662 - x2297)
            - x0 * (3 * x1664 - 3 * x1665 - 3 * x1667 + x2297)
            + x2295 * x65
            + x2296,
            x0 * (3 * x1670 - 3 * x1672 + 3 * x1674 - x2300)
            - x0 * (3 * x1676 - 3 * x1677 - 3 * x1679 + x2300)
            + x2298 * x65
            + x2299,
            x0 * (3 * x1682 - 3 * x1684 + 3 * x1686 - x2303)
            - x0 * (3 * x1688 - 3 * x1689 - 3 * x1691 + x2303)
            + x2301 * x65
            + x2302,
            x0 * (3 * x1694 - 3 * x1696 + 3 * x1698 - x2306)
            - x0 * (3 * x1700 - 3 * x1701 - 3 * x1703 + x2306)
            + x2304 * x65
            + x2305,
            x2276,
            x0 * (2 * x1707 + 2 * x1709 - 2 * x1710 - x2309)
            - x0 * (-2 * x1712 + 2 * x1713 - 2 * x1714 + x2309)
            + x2307 * x65
            + x2308,
            x0 * (2 * x1718 + 2 * x1720 - 2 * x1721 - x2312)
            - x0 * (-2 * x1723 + 2 * x1724 - 2 * x1725 + x2312)
            + x2310 * x65
            + x2311,
            x0 * (2 * x1729 + 2 * x1731 - 2 * x1732 - x2315)
            - x0 * (-2 * x1734 + 2 * x1735 - 2 * x1736 + x2315)
            + x2313 * x65
            + x2314,
            x0 * (2 * x1740 + 2 * x1742 - 2 * x1743 - x2318)
            - x0 * (-2 * x1745 + 2 * x1746 - 2 * x1747 + x2318)
            + x2316 * x65
            + x2317,
            x0 * (2 * x1751 + 2 * x1753 - 2 * x1754 - x2321)
            - x0 * (-2 * x1756 + 2 * x1757 - 2 * x1758 + x2321)
            + x2319 * x65
            + x2320,
            x0 * (2 * x1762 + 2 * x1764 - 2 * x1765 - x2324)
            - x0 * (-2 * x1767 + 2 * x1768 - 2 * x1769 + x2324)
            + x2322 * x65
            + x2323,
            x0 * (2 * x1773 + 2 * x1775 - 2 * x1776 - x2327)
            - x0 * (-2 * x1778 + 2 * x1779 - 2 * x1780 + x2327)
            + x2325 * x65
            + x2326,
            x0 * (2 * x1784 + 2 * x1786 - 2 * x1787 - x2330)
            - x0 * (-2 * x1789 + 2 * x1790 - 2 * x1791 + x2330)
            + x2328 * x65
            + x2329,
            x0 * (2 * x1795 + 2 * x1797 - 2 * x1798 - x2333)
            - x0 * (-2 * x1800 + 2 * x1801 - 2 * x1802 + x2333)
            + x2331 * x65
            + x2332,
            x0 * (2 * x1806 + 2 * x1808 - 2 * x1809 - x2336)
            - x0 * (-2 * x1811 + 2 * x1812 - 2 * x1813 + x2336)
            + x2334 * x65
            + x2335,
            x0 * (2 * x1817 + 2 * x1819 - 2 * x1820 - x2339)
            - x0 * (-2 * x1822 + 2 * x1823 - 2 * x1824 + x2339)
            + x2337 * x65
            + x2338,
            x0 * (2 * x1828 + 2 * x1830 - 2 * x1831 - x2342)
            - x0 * (-2 * x1833 + 2 * x1834 - 2 * x1835 + x2342)
            + x2340 * x65
            + x2341,
            x0 * (2 * x1839 + 2 * x1841 - 2 * x1842 - x2345)
            - x0 * (-2 * x1844 + 2 * x1845 - 2 * x1846 + x2345)
            + x2343 * x65
            + x2344,
            x0 * (2 * x1850 + 2 * x1852 - 2 * x1853 - x2348)
            - x0 * (-2 * x1855 + 2 * x1856 - 2 * x1857 + x2348)
            + x2346 * x65
            + x2347,
            x0 * (2 * x1861 + 2 * x1863 - 2 * x1864 - x2351)
            - x0 * (-2 * x1866 + 2 * x1867 - 2 * x1868 + x2351)
            + x2349 * x65
            + x2350,
            x0 * x1891 - x0 * x1894 + x2352 * x65 + x2353,
            x0 * x1956 - x0 * x1965 + x2354 * x65 + x2355,
            x0 * x2033 - x0 * x2042 + x2356 * x65 + x2357,
            x0 * x2082 - x0 * x2088 + x2358 * x65 + x2359,
            x0 * x2128 - x0 * x2134 + x2360 * x65 + x2361,
            x0 * x2174 - x0 * x2180 + x2362 * x65 + x2363,
            x0 * x2200 - x0 * x2203 + x2364 * x65 + x2365,
            x0 * x2223 - x0 * x2226 + x2366 * x65 + x2367,
            x0 * x2246 - x0 * x2249 + x2368 * x65 + x2369,
            x0 * x2269 - x0 * x2272 + x2370 * x65 + x2371,
            x0 * x691 - x0 * x701 + x2372 * x65 + x2373,
            x0 * x753 - x0 * x767 + x2374 * x65 + x2375,
            x0 * x827 - x0 * x840 + x2376 * x65 + x2377,
            x0 * x895 - x0 * x908 + x2378 * x65 + x2379,
            x0 * x117 - x0 * x133 + x2380 * x65 + x2381,
            x2424,
            x2475,
            x2526,
            x2553,
            x2580,
            x2607,
            x2618,
            x2629,
            x2640,
            x2651,
            x704,
            x770,
            x843,
            x911,
            x136,
            x2652,
            x2277 * x6 + x2278,
            x2652,
            x2280 * x6 + x2281,
            x2283 * x6 + x2284,
            x2652,
            x2286 * x6 + x2287,
            x2289 * x6 + x2290,
            x2292 * x6 + x2293,
            x2652,
            x2295 * x6 + x2296,
            x2298 * x6 + x2299,
            x2301 * x6 + x2302,
            x2304 * x6 + x2305,
            x2652,
            x2307 * x6 + x2308,
            x2310 * x6 + x2311,
            x2313 * x6 + x2314,
            x2316 * x6 + x2317,
            x2319 * x6 + x2320,
            x2322 * x6 + x2323,
            x2325 * x6 + x2326,
            x2328 * x6 + x2329,
            x2331 * x6 + x2332,
            x2334 * x6 + x2335,
            x2337 * x6 + x2338,
            x2340 * x6 + x2341,
            x2343 * x6 + x2344,
            x2346 * x6 + x2347,
            x2349 * x6 + x2350,
            x2352 * x6 + x2353,
            x2354 * x6 + x2355,
            x2356 * x6 + x2357,
            x2358 * x6 + x2359,
            x2360 * x6 + x2361,
            x2362 * x6 + x2363,
            x2364 * x6 + x2365,
            x2366 * x6 + x2367,
            x2368 * x6 + x2369,
            x2370 * x6 + x2371,
            x2372 * x6 + x2373,
            x2374 * x6 + x2375,
            x2376 * x6 + x2377,
            x2378 * x6 + x2379,
            x2380 * x6 + x2381,
            x2418 * x6 + x2423,
            x2466 * x6 + x2474,
            x2517 * x6 + x2525,
            x2547 * x6 + x2552,
            x2574 * x6 + x2579,
            x2601 * x6 + x2606,
            x2615 * x6 + x2617,
            x2626 * x6 + x2628,
            x2637 * x6 + x2639,
            x2648 * x6 + x2650,
            x6 * x692 + x703,
            x6 * x754 + x769,
            x6 * x828 + x842,
            x6 * x896 + x910,
            x118 * x6 + x135,
            x0
            * (
                3 * x1237
                - 3 * x1241
                + 3 * x1877
                - x2653
                - x2654 * x416
                + x2666
                * (
                    3 * x1876
                    - x2655
                    + x2659
                    * (
                        x1875
                        + x229 * (-x2657 - 3 * x329 + 3 * x484 + 3 * x486 + x647)
                        - x2658
                        + x743
                    )
                    - x2662
                    + x2663 * (-x2394 + x2402 + x793)
                    - x2664 * x2665
                    + 3 * x860
                    - 3 * x861
                )
                - x2672
                + x414 * (-x2396 - x2399 + x2403 + x2404)
            )
            - x0
            * (
                3 * x1242
                - 3 * x1246
                - 3 * x1888
                + x2653
                + x2654 * x414
                - x2666
                * (
                    -3 * x1886
                    - x2659
                    * (
                        x1885
                        - x229 * (x2668 - 3 * x386 - 3 * x526 + 3 * x527 + x694)
                        + x2669
                        + 2 * x756
                        - 2 * x760
                    )
                    + x2663 * x2671
                    - x2665 * (x2408 - x2419 + x835)
                    + x2667
                    + x2670
                    + 3 * x877
                    - 3 * x883
                )
                + x2672
                - x416 * (x2411 - x2420 - x2421 + x2422)
            )
            + x2424,
            x0
            * (
                3 * x1279
                - 3 * x1289
                + 3 * x1922
                + x2659
                * (
                    x138 * (-x2430 + x2446 + x2447 - x2448)
                    - x18 * x2682
                    + x1921
                    + x229
                    * (
                        x0 * (3 * x1914 - 3 * x1915 - x2674 + 3 * x320)
                        - x17 * x2429
                        - 3 * x1913
                        + 3 * x1917
                        + 3 * x1918
                        + x2445 * x65
                        - x2676
                        - x2679
                    )
                    - x2681
                )
                + x2663 * (-x2440 - x2443 + x2449 + x2450)
                - x2665 * x2689
                - x2673
                - x2688
            )
            - x0
            * (
                3 * x1290
                - 3 * x1300
                - 3 * x1953
                - x2659
                * (
                    x138 * x2687
                    - x18 * (x2454 - x2468 + x2469 - x2470)
                    + x1952
                    - x229
                    * (
                        -x0 * (3 * x1939 - 3 * x1940 + x2683 - 3 * x365)
                        - x17 * x2467
                        - 3 * x1942
                        - 3 * x1945
                        + 3 * x1946
                        + x2453 * x65
                        + x2684
                        + x2685
                    )
                    + x2686
                )
                + x2663 * x2689
                - x2665 * (x2459 - x2471 - x2472 + x2473)
                + x2673
                + x2688
            )
            + x2475,
            x0
            * (
                3 * x1335
                - 3 * x1344
                + 3 * x1997
                + x2659
                * (
                    x138 * (-x2487 - x2496 + x2498 + x2499)
                    - x18 * x2699
                    + x1996
                    + x229
                    * (
                        x0 * (3 * x1988 + 3 * x1989 - 3 * x1990 - x2691)
                        - x17 * x2486
                        - 3 * x1987
                        + 3 * x1992
                        + 3 * x1993
                        + x2497 * x65
                        - x2693
                        - x2696
                    )
                    - x2698
                )
                + x2663 * (-x2491 - x2494 + x2500 + x2501)
                - x2665 * x2706
                - x2690
                - x2705
            )
            - x0
            * (
                3 * x1345
                - 3 * x1354
                - 3 * x2030
                - x2659
                * (
                    x138 * x2704
                    - x18 * (x2505 - x2519 - x2520 + x2521)
                    + x2029
                    - x229
                    * (
                        -x0 * (-3 * x2015 + 3 * x2016 - 3 * x2017 + x2700)
                        - x17 * x2518
                        - 3 * x2019
                        - 3 * x2022
                        + 3 * x2023
                        + x2504 * x65
                        + x2701
                        + x2702
                    )
                    + x2703
                )
                + x2663 * x2706
                - x2665 * (x2510 - x2522 - x2523 + x2524)
                + x2690
                + x2705
            )
            + x2526,
            x0
            * (
                3 * x1378
                + x138 * (-x2531 - x2534 + x2536 + x2537)
                - 3 * x1385
                - x18 * x2716
                + 3 * x2058
                + x229
                * (
                    x0 * (3 * x2052 + 3 * x2053 - 3 * x2054 - x2708)
                    - x17 * x2530
                    - 3 * x2051
                    + 3 * x2056
                    + 3 * x2057
                    + x2535 * x65
                    - x2710
                    - x2711
                )
                - x2707
                - x2715
            )
            - x0
            * (
                x138 * x2716
                + 3 * x1386
                - 3 * x1393
                - x18 * (x2541 - x2549 - x2550 + x2551)
                - 3 * x2079
                - x229
                * (
                    -x0 * (-3 * x2070 + 3 * x2071 - 3 * x2072 + x2712)
                    - x17 * x2548
                    - 3 * x2074
                    - 3 * x2077
                    + 3 * x2078
                    + x2540 * x65
                    + x2713
                    + x2714
                )
                + x2707
                + x2715
            )
            + x2553,
            -x0
            * (
                x138 * x2726
                + 3 * x1425
                - 3 * x1432
                - x18 * (x2568 - x2576 - x2577 + x2578)
                - 3 * x2125
                - x229
                * (
                    -x0 * (-3 * x2116 + 3 * x2117 - 3 * x2118 + x2722)
                    - x17 * x2575
                    - 3 * x2120
                    - 3 * x2123
                    + 3 * x2124
                    + x2567 * x65
                    + x2723
                    + x2724
                )
                + x2717
                + x2725
            )
            + x0
            * (
                x138 * (-x2558 - x2561 + x2563 + x2564)
                + 3 * x1417
                - 3 * x1424
                - x18 * x2726
                + 3 * x2104
                + x229
                * (
                    x0 * (3 * x2098 + 3 * x2099 - 3 * x2100 - x2718)
                    - x17 * x2557
                    - 3 * x2097
                    + 3 * x2102
                    + 3 * x2103
                    + x2562 * x65
                    - x2720
                    - x2721
                )
                - x2717
                - x2725
            )
            + x2580,
            -x0
            * (
                x138 * x2736
                + 3 * x1464
                - 3 * x1471
                - x18 * (x2595 - x2603 - x2604 + x2605)
                - 3 * x2171
                - x229
                * (
                    -x0 * (-3 * x2162 + 3 * x2163 - 3 * x2164 + x2732)
                    - x17 * x2602
                    - 3 * x2166
                    - 3 * x2169
                    + 3 * x2170
                    + x2594 * x65
                    + x2733
                    + x2734
                )
                + x2727
                + x2735
            )
            + x0
            * (
                x138 * (-x2585 - x2588 + x2590 + x2591)
                + 3 * x1456
                - 3 * x1463
                - x18 * x2736
                + 3 * x2150
                + x229
                * (
                    x0 * (3 * x2144 + 3 * x2145 - 3 * x2146 - x2728)
                    - x17 * x2584
                    - 3 * x2143
                    + 3 * x2148
                    + 3 * x2149
                    + x2589 * x65
                    - x2730
                    - x2731
                )
                - x2727
                - x2735
            )
            + x2607,
            x0
            * (
                x0 * (3 * x2181 + 3 * x2183 - 3 * x2184 - x2737)
                + 3 * x1485
                - 3 * x1488
                - x17 * x2613
                + 3 * x2186
                + x2610 * x65
                - x2739
                - x2740
            )
            - x0
            * (
                -x0 * (-3 * x2193 + 3 * x2194 - 3 * x2195 + x2738)
                + 3 * x1489
                - 3 * x1492
                - x17 * x2616
                - 3 * x2197
                + x2613 * x65
                + x2739
                + x2740
            )
            + x2618,
            x0
            * (
                x0 * (3 * x2204 + 3 * x2206 - 3 * x2207 - x2741)
                + 3 * x1505
                - 3 * x1508
                - x17 * x2624
                + 3 * x2209
                + x2621 * x65
                - x2743
                - x2744
            )
            - x0
            * (
                -x0 * (-3 * x2216 + 3 * x2217 - 3 * x2218 + x2742)
                + 3 * x1509
                - 3 * x1512
                - x17 * x2627
                - 3 * x2220
                + x2624 * x65
                + x2743
                + x2744
            )
            + x2629,
            x0
            * (
                x0 * (3 * x2227 + 3 * x2229 - 3 * x2230 - x2745)
                + 3 * x1525
                - 3 * x1528
                - x17 * x2635
                + 3 * x2232
                + x2632 * x65
                - x2747
                - x2748
            )
            - x0
            * (
                -x0 * (-3 * x2239 + 3 * x2240 - 3 * x2241 + x2746)
                + 3 * x1529
                - 3 * x1532
                - x17 * x2638
                - 3 * x2243
                + x2635 * x65
                + x2747
                + x2748
            )
            + x2640,
            x0
            * (
                x0 * (3 * x2250 + 3 * x2252 - 3 * x2253 - x2749)
                + 3 * x1545
                - 3 * x1548
                - x17 * x2646
                + 3 * x2255
                + x2643 * x65
                - x2751
                - x2752
            )
            - x0
            * (
                -x0 * (-3 * x2262 + 3 * x2263 - 3 * x2264 + x2750)
                + 3 * x1549
                - 3 * x1552
                - x17 * x2649
                - 3 * x2266
                + x2646 * x65
                + x2751
                + x2752
            )
            + x2651,
            x705,
            x771,
            x844,
            x912,
            x137,
        ]
    )
    return S
