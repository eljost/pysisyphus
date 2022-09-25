import numpy

"""

    Quadrupole integrals contain the upper triangular part of the symmetric
    3x3 quadrupole matrix.

    / xx xy xz \
    |    yy yz |
    \       zz /

    for cart_dir1 in (x, y, z):
        for bf_a in basis_functions_a:
            for bf_b in basis_functions_b:
                    quadrupole_integrals(bf_a, bf_b, cart_dir1)

"""


def quadrupole3d_00(a, A, b, B, C):
    """Cartesian 3D (ss) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = a * b * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x5 = numpy.sqrt(x0)
    x6 = numpy.sqrt(numpy.pi) * x5
    x7 = x4 * x6
    x8 = x0 * (a * A[0] + b * B[0]) - C[0]
    x9 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x10 = numpy.pi * x0 * x9
    x11 = x0 * (a * A[1] + b * B[1]) - C[1]
    x12 = numpy.pi ** (3 / 2)
    x13 = x0 * x2 * x4
    x14 = x12 * x13 * x5 * x8 * x9
    x15 = x0 * (a * A[2] + b * B[2]) - C[2]
    x16 = x2 * x6
    x17 = x6 * x9

    # 6 item(s)
    return numpy.array(
        [
            x10 * x2 * (x3 * x7 + x7 * x8 ** 2),
            x11 * x14,
            x14 * x15,
            x10 * x4 * (x11 ** 2 * x16 + x16 * x3),
            x11 * x12 * x13 * x15 * x5 * x9,
            numpy.pi * x13 * (x15 ** 2 * x17 + x17 * x3),
        ]
    )


def quadrupole3d_01(a, A, b, B, C):
    """Cartesian 3D (sp) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - C[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = a * b * x0
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x0)
    x7 = numpy.sqrt(numpy.pi) * x6
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x1 - B[0]
    x11 = x2 ** 2 * x8 + x9
    x12 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x13 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x14 = numpy.pi * x0 * x13
    x15 = x12 * x14
    x16 = -x0 * (a * A[1] + b * B[1])
    x17 = -x16 - B[1]
    x18 = x11 * x15
    x19 = -x0 * (a * A[2] + b * B[2])
    x20 = -x19 - B[2]
    x21 = -x16 - C[1]
    x22 = x15 * (x10 * x2 * x8 + x9)
    x23 = x3 * x7
    x24 = x12 * x23
    x25 = x12 * x7
    x26 = x14 * x5
    x27 = x26 * (x17 * x21 * x25 + x24)
    x28 = numpy.pi ** (3 / 2)
    x29 = x0 * x12 * x5
    x30 = x13 * x2 * x28 * x29 * x6
    x31 = -x19 - C[2]
    x32 = x13 * x23
    x33 = x13 * x7
    x34 = numpy.pi * x29
    x35 = x34 * (x20 * x31 * x33 + x32)
    x36 = x21 ** 2 * x25 + x24
    x37 = x26 * x36
    x38 = x31 ** 2 * x33 + x32
    x39 = x34 * x38

    # 18 item(s)
    return numpy.array(
        [
            x15 * (x10 * x11 + 2 * x2 * x9),
            x17 * x18,
            x18 * x20,
            x21 * x22,
            x2 * x27,
            x20 * x21 * x30,
            x22 * x31,
            x17 * x30 * x31,
            x2 * x35,
            x10 * x37,
            x26 * (x17 * x36 + 2 * x21 * x24),
            x20 * x37,
            x10 * x13 * x21 * x28 * x29 * x31 * x6,
            x27 * x31,
            x21 * x35,
            x10 * x39,
            x17 * x39,
            x34 * (x20 * x38 + 2 * x31 * x32),
        ]
    )


def quadrupole3d_02(a, A, b, B, C):
    """Cartesian 3D (sd) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = a * b * x1
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x5 = x3 * x4
    x6 = x0 * x5
    x7 = -x1 * (a * A[0] + b * B[0])
    x8 = -x7 - C[0]
    x9 = x5 * x8 ** 2
    x10 = -x7 - B[0]
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
    x21 = x14 * x18
    x22 = -x1 * (a * A[2] + b * B[2])
    x23 = -x22 - B[2]
    x24 = x15 * x4
    x25 = x0 * x24
    x26 = x20 ** 2 * x24 + x25
    x27 = x16 * x4
    x28 = x18 * x23
    x29 = x0 * x27
    x30 = x23 ** 2 * x27 + x29
    x31 = -x19 - C[1]
    x32 = x12 + x6
    x33 = x18 * (x0 * (x11 + x5 * x8) + x10 * x32)
    x34 = x20 * x24
    x35 = x31 * x34
    x36 = x25 + x35
    x37 = x17 * x3
    x38 = x37 * (x0 * (x24 * x31 + x34) + x20 * x36)
    x39 = x37 * x8
    x40 = numpy.pi * x1 * x15 * x3
    x41 = x40 * x8
    x42 = -x22 - C[2]
    x43 = x18 * x42
    x44 = x23 * x27
    x45 = x42 * x44
    x46 = x29 + x45
    x47 = x40 * (x0 * (x27 * x42 + x44) + x23 * x46)
    x48 = x24 * x31 ** 2
    x49 = x25 + x48
    x50 = x10 ** 2 * x5 + x6
    x51 = x20 * x49 + 2 * x25 * x31
    x52 = x37 * x51
    x53 = x10 * x37
    x54 = x10 * x40
    x55 = x27 * x42 ** 2
    x56 = x29 + x55
    x57 = x23 * x56 + 2 * x29 * x42
    x58 = x40 * x57

    # 36 item(s)
    return numpy.array(
        [
            x18 * (x0 * (2 * x12 + 3 * x6 + x9) + x10 * x14),
            x20 * x21,
            x21 * x23,
            x13 * x26 * x27,
            x13 * x20 * x28,
            x13 * x24 * x30,
            x31 * x33,
            x27 * x32 * x36,
            x28 * x31 * x32,
            x38 * x8,
            x23 * x36 * x39,
            x30 * x31 * x41,
            x33 * x42,
            x20 * x32 * x43,
            x24 * x32 * x46,
            x26 * x39 * x42,
            x20 * x41 * x46,
            x47 * x8,
            x27 * x49 * x50,
            x10 * x52,
            x23 * x49 * x53,
            x37 * (x0 * (3 * x25 + 2 * x35 + x48) + x20 * x51),
            x23 * x52,
            x30 * x49 * x5,
            x31 * x43 * x50,
            x36 * x42 * x53,
            x31 * x46 * x54,
            x38 * x42,
            x36 * x46 * x5,
            x31 * x47,
            x24 * x50 * x56,
            x20 * x54 * x56,
            x10 * x58,
            x26 * x5 * x56,
            x20 * x58,
            x40 * (x0 * (3 * x29 + 2 * x45 + x55) + x23 * x57),
        ]
    )


def quadrupole3d_03(a, A, b, B, C):
    """Cartesian 3D (sf) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - B[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - C[0]
    x5 = a * b * x0
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x8 = x6 * x7
    x9 = x4 ** 2 * x8
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
    x30 = x20 * x27
    x31 = -x0 * (a * A[2] + b * B[2])
    x32 = -x31 - B[2]
    x33 = x24 * x7
    x34 = x3 * x33
    x35 = x29 ** 2 * x33
    x36 = x34 + x35
    x37 = x25 * x7
    x38 = x27 * x32
    x39 = x3 * x37
    x40 = x32 ** 2 * x37
    x41 = x39 + x40
    x42 = 2 * x34
    x43 = x29 * x36 + x29 * x42
    x44 = x32 * x37
    x45 = x29 * x33
    x46 = 2 * x39
    x47 = x32 * x41 + x32 * x46
    x48 = -x28 - C[1]
    x49 = x2 ** 2 * x8
    x50 = x21 + x23
    x51 = x27 * (x2 * x50 + x3 * (x14 + x49))
    x52 = x45 * x48
    x53 = x34 + x52
    x54 = x33 * x48
    x55 = x3 * (x45 + x54)
    x56 = x29 * x53
    x57 = x55 + x56
    x58 = 3 * x34 + 2 * x52
    x59 = x26 * x6
    x60 = x59 * (x29 * x57 + x3 * (x35 + x58))
    x61 = x32 * x59
    x62 = numpy.pi * x0 * x24
    x63 = x6 * x62
    x64 = -x31 - C[2]
    x65 = x27 * x64
    x66 = x44 * x64
    x67 = x39 + x66
    x68 = x37 * x64
    x69 = x3 * (x44 + x68)
    x70 = x32 * x67
    x71 = x69 + x70
    x72 = x29 * x63
    x73 = 3 * x39 + 2 * x66
    x74 = x63 * (x3 * (x40 + x73) + x32 * x71)
    x75 = x33 * x48 ** 2
    x76 = x34 + x75
    x77 = x10 + x49
    x78 = 2 * x12 * x3 + x2 * x77
    x79 = x29 * x76
    x80 = x42 * x48 + x79
    x81 = x29 * x80 + x3 * (x58 + x75)
    x82 = x11 * x26
    x83 = x11 * x62
    x84 = x37 * x64 ** 2
    x85 = x39 + x84
    x86 = x32 * x85
    x87 = x46 * x64 + x86
    x88 = x3 * (x73 + x84) + x32 * x87

    # 60 item(s)
    return numpy.array(
        [
            x27 * (x2 * x20 + x3 * (4 * x16 + 2 * x18 + 2 * x21 + 2 * x23)),
            x29 * x30,
            x30 * x32,
            x19 * x36 * x37,
            x19 * x29 * x38,
            x19 * x33 * x41,
            x17 * x37 * x43,
            x17 * x36 * x44,
            x17 * x41 * x45,
            x17 * x33 * x47,
            x48 * x51,
            x37 * x50 * x53,
            x38 * x48 * x50,
            x22 * x37 * x57,
            x22 * x44 * x53,
            x22 * x41 * x54,
            x4 * x60,
            x4 * x57 * x61,
            x15 * x41 * x53,
            x4 * x47 * x48 * x63,
            x51 * x64,
            x29 * x50 * x65,
            x33 * x50 * x67,
            x22 * x36 * x68,
            x22 * x45 * x67,
            x22 * x33 * x71,
            x4 * x43 * x59 * x64,
            x15 * x36 * x67,
            x4 * x71 * x72,
            x4 * x74,
            x37 * x76 * x78,
            x37 * x77 * x80,
            x44 * x76 * x77,
            x81 * x82,
            x32 * x80 * x82,
            x12 * x41 * x76,
            x59 * (x29 * x81 + x3 * (4 * x34 * x48 + 2 * x55 + 2 * x56 + 2 * x79)),
            x61 * x81,
            x41 * x8 * x80,
            x47 * x76 * x8,
            x48 * x65 * x78,
            x53 * x68 * x77,
            x54 * x67 * x77,
            x57 * x64 * x82,
            x12 * x53 * x67,
            x48 * x71 * x83,
            x60 * x64,
            x57 * x67 * x8,
            x53 * x71 * x8,
            x48 * x74,
            x33 * x78 * x85,
            x45 * x77 * x85,
            x33 * x77 * x87,
            x12 * x36 * x85,
            x29 * x83 * x87,
            x83 * x88,
            x43 * x8 * x85,
            x36 * x8 * x87,
            x72 * x88,
            x63 * (x3 * (4 * x39 * x64 + 2 * x69 + 2 * x70 + 2 * x86) + x32 * x88),
        ]
    )


def quadrupole3d_04(a, A, b, B, C):
    """Cartesian 3D (sg) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3 ** 2 * x7
    x9 = x0 * x7
    x10 = 3 * x9
    x11 = -x2 - C[0]
    x12 = x3 * x7
    x13 = x11 * x12
    x14 = x10 + 2 * x13
    x15 = x0 * (x14 + x8)
    x16 = x11 ** 2 * x7
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
    x37 = x30 * x34
    x38 = -x1 * (a * A[2] + b * B[2])
    x39 = -x38 - B[2]
    x40 = x31 * x6
    x41 = x0 * x40
    x42 = x36 ** 2 * x40
    x43 = x41 + x42
    x44 = x32 * x6
    x45 = x34 * x39
    x46 = x0 * x44
    x47 = x39 ** 2 * x44
    x48 = x46 + x47
    x49 = 2 * x41
    x50 = x36 * x43 + x36 * x49
    x51 = x39 * x44
    x52 = x36 * x40
    x53 = 2 * x46
    x54 = x39 * x48 + x39 * x53
    x55 = 3 * x41
    x56 = x0 * (3 * x42 + x55) + x36 * x50
    x57 = 3 * x46
    x58 = x0 * (3 * x47 + x57) + x39 * x54
    x59 = -x35 - C[1]
    x60 = x8 + x9
    x61 = 2 * x0 * x12 + x3 * x60
    x62 = x15 + x28
    x63 = x34 * (x0 * (3 * x24 + 3 * x26 + x61) + x3 * x62)
    x64 = x52 * x59
    x65 = x41 + x64
    x66 = x40 * x59
    x67 = x0 * (x52 + x66)
    x68 = x36 * x65
    x69 = x67 + x68
    x70 = x55 + 2 * x64
    x71 = x0 * (x42 + x70)
    x72 = x36 * x69
    x73 = x71 + x72
    x74 = x33 * x5
    x75 = x74 * (x0 * (x50 + 3 * x67 + 3 * x68) + x36 * x73)
    x76 = x11 * x74
    x77 = numpy.pi * x1 * x31 * x5
    x78 = x11 * x77
    x79 = -x38 - C[2]
    x80 = x34 * x79
    x81 = x51 * x79
    x82 = x46 + x81
    x83 = x44 * x79
    x84 = x0 * (x51 + x83)
    x85 = x39 * x82
    x86 = x84 + x85
    x87 = x57 + 2 * x81
    x88 = x0 * (x47 + x87)
    x89 = x39 * x86
    x90 = x88 + x89
    x91 = x77 * (x0 * (x54 + 3 * x84 + 3 * x85) + x39 * x90)
    x92 = x40 * x59 ** 2
    x93 = x41 + x92
    x94 = x0 * (x10 + 3 * x8) + x3 * x61
    x95 = x36 * x93
    x96 = x49 * x59 + x95
    x97 = x0 * (x70 + x92)
    x98 = x36 * x96
    x99 = x97 + x98
    x100 = x0 * (4 * x41 * x59 + 2 * x67 + 2 * x68 + 2 * x95) + x36 * x99
    x101 = x100 * x74
    x102 = x3 * x74
    x103 = x3 * x77
    x104 = x44 * x79 ** 2
    x105 = x104 + x46
    x106 = x105 * x39
    x107 = x106 + x53 * x79
    x108 = x0 * (x104 + x87)
    x109 = x107 * x39
    x110 = x108 + x109
    x111 = x0 * (2 * x106 + 4 * x46 * x79 + 2 * x84 + 2 * x85) + x110 * x39
    x112 = x111 * x77

    # 90 item(s)
    return numpy.array(
        [
            x34 * (x0 * (2 * x15 + 3 * x17 + 3 * x23 + 2 * x28) + x3 * x30),
            x36 * x37,
            x37 * x39,
            x29 * x43 * x44,
            x29 * x36 * x45,
            x29 * x40 * x48,
            x22 * x44 * x50,
            x22 * x43 * x51,
            x22 * x48 * x52,
            x22 * x40 * x54,
            x20 * x44 * x56,
            x20 * x50 * x51,
            x20 * x43 * x48,
            x20 * x52 * x54,
            x20 * x40 * x58,
            x59 * x63,
            x44 * x62 * x65,
            x45 * x59 * x62,
            x27 * x44 * x69,
            x27 * x51 * x65,
            x27 * x48 * x66,
            x25 * x44 * x73,
            x25 * x51 * x69,
            x25 * x48 * x65,
            x25 * x54 * x66,
            x11 * x75,
            x39 * x73 * x76,
            x18 * x48 * x69,
            x18 * x54 * x65,
            x58 * x59 * x78,
            x63 * x79,
            x36 * x62 * x80,
            x40 * x62 * x82,
            x27 * x43 * x83,
            x27 * x52 * x82,
            x27 * x40 * x86,
            x25 * x50 * x83,
            x25 * x43 * x82,
            x25 * x52 * x86,
            x25 * x40 * x90,
            x56 * x76 * x79,
            x18 * x50 * x82,
            x18 * x43 * x86,
            x36 * x78 * x90,
            x11 * x91,
            x44 * x93 * x94,
            x44 * x61 * x96,
            x51 * x61 * x93,
            x44 * x60 * x99,
            x51 * x60 * x96,
            x48 * x60 * x93,
            x101 * x3,
            x102 * x39 * x99,
            x12 * x48 * x96,
            x12 * x54 * x93,
            x74 * (x0 * (2 * x71 + 2 * x72 + 3 * x97 + 3 * x98) + x100 * x36),
            x101 * x39,
            x48 * x7 * x99,
            x54 * x7 * x96,
            x58 * x7 * x93,
            x59 * x80 * x94,
            x61 * x65 * x83,
            x61 * x66 * x82,
            x60 * x69 * x83,
            x60 * x65 * x82,
            x60 * x66 * x86,
            x102 * x73 * x79,
            x12 * x69 * x82,
            x12 * x65 * x86,
            x103 * x59 * x90,
            x75 * x79,
            x7 * x73 * x82,
            x69 * x7 * x86,
            x65 * x7 * x90,
            x59 * x91,
            x105 * x40 * x94,
            x105 * x52 * x61,
            x107 * x40 * x61,
            x105 * x43 * x60,
            x107 * x52 * x60,
            x110 * x40 * x60,
            x105 * x12 * x50,
            x107 * x12 * x43,
            x103 * x110 * x36,
            x112 * x3,
            x105 * x56 * x7,
            x107 * x50 * x7,
            x110 * x43 * x7,
            x112 * x36,
            x77 * (x0 * (3 * x108 + 3 * x109 + 2 * x88 + 2 * x89) + x111 * x39),
        ]
    )


def quadrupole3d_10(a, A, b, B, C):
    """Cartesian 3D (ps) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - C[0]
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
    x17 = -x16 - A[1]
    x18 = x11 * x15
    x19 = -x0 * (a * A[2] + b * B[2])
    x20 = -x19 - A[2]
    x21 = -x16 - C[1]
    x22 = x15 * (x10 * x2 * x8 + x9)
    x23 = x3 * x7
    x24 = x12 * x23
    x25 = x12 * x7
    x26 = x14 * x5
    x27 = x26 * (x17 * x21 * x25 + x24)
    x28 = numpy.pi ** (3 / 2)
    x29 = x0 * x12 * x5
    x30 = x13 * x2 * x28 * x29 * x6
    x31 = -x19 - C[2]
    x32 = x13 * x23
    x33 = x13 * x7
    x34 = numpy.pi * x29
    x35 = x34 * (x20 * x31 * x33 + x32)
    x36 = x21 ** 2 * x25 + x24
    x37 = x26 * x36
    x38 = x31 ** 2 * x33 + x32
    x39 = x34 * x38

    # 18 item(s)
    return numpy.array(
        [
            x15 * (x10 * x11 + 2 * x2 * x9),
            x17 * x18,
            x18 * x20,
            x21 * x22,
            x2 * x27,
            x20 * x21 * x30,
            x22 * x31,
            x17 * x30 * x31,
            x2 * x35,
            x10 * x37,
            x26 * (x17 * x36 + 2 * x21 * x24),
            x20 * x37,
            x10 * x13 * x21 * x28 * x29 * x31 * x6,
            x27 * x31,
            x21 * x35,
            x10 * x39,
            x17 * x39,
            x34 * (x20 * x38 + 2 * x31 * x32),
        ]
    )


def quadrupole3d_11(a, A, b, B, C):
    """Cartesian 3D (pp) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = a * b * x1
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x5 = x3 * x4
    x6 = x0 * x5
    x7 = -x1 * (a * A[0] + b * B[0])
    x8 = -x7 - C[0]
    x9 = x5 * x8 ** 2
    x10 = -x7 - B[0]
    x11 = x10 * x5
    x12 = x11 * x8
    x13 = -x7 - A[0]
    x14 = 2 * x6 * x8
    x15 = x6 + x9
    x16 = x10 * x15 + x14
    x17 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x18 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x19 = numpy.pi * x1 * x18
    x20 = x17 * x19
    x21 = -x1 * (a * A[1] + b * B[1])
    x22 = -x21 - B[1]
    x23 = x20 * (x13 * x15 + x14)
    x24 = -x1 * (a * A[2] + b * B[2])
    x25 = -x24 - B[2]
    x26 = -x21 - A[1]
    x27 = x16 * x20
    x28 = x0 * x4
    x29 = x17 * x28
    x30 = x17 * x4
    x31 = x22 * x30
    x32 = x26 * x31 + x29
    x33 = x18 * x4
    x34 = x15 * x20
    x35 = -x24 - A[2]
    x36 = x18 * x28
    x37 = x25 * x33
    x38 = x35 * x37 + x36
    x39 = -x21 - C[1]
    x40 = x5 * x8
    x41 = x12 + x6
    x42 = x20 * (x0 * (x11 + x40) + x13 * x41)
    x43 = x31 * x39
    x44 = x29 + x43
    x45 = x13 * x40 + x6
    x46 = x20 * x39
    x47 = x30 * x39
    x48 = x26 * x47 + x29
    x49 = x19 * x3
    x50 = x49 * (x0 * (x31 + x47) + x26 * x44)
    x51 = x49 * x8
    x52 = numpy.pi * x1 * x17 * x3
    x53 = x52 * x8
    x54 = -x24 - C[2]
    x55 = x20 * x54
    x56 = x37 * x54
    x57 = x36 + x56
    x58 = x33 * x54
    x59 = x35 * x58 + x36
    x60 = x52 * (x0 * (x37 + x58) + x35 * x57)
    x61 = x30 * x39 ** 2
    x62 = x29 + x61
    x63 = x11 * x13 + x6
    x64 = 2 * x29 * x39
    x65 = x22 * x62 + x64
    x66 = x49 * x65
    x67 = x49 * x62
    x68 = x49 * (x26 * x62 + x64)
    x69 = x49 * x54
    x70 = x39 * x52
    x71 = x33 * x54 ** 2
    x72 = x36 + x71
    x73 = x52 * x72
    x74 = 2 * x36 * x54
    x75 = x25 * x72 + x74
    x76 = x52 * x75
    x77 = x52 * (x35 * x72 + x74)

    # 54 item(s)
    return numpy.array(
        [
            x20 * (x0 * (2 * x12 + 3 * x6 + x9) + x13 * x16),
            x22 * x23,
            x23 * x25,
            x26 * x27,
            x15 * x32 * x33,
            x25 * x26 * x34,
            x27 * x35,
            x22 * x34 * x35,
            x15 * x30 * x38,
            x39 * x42,
            x33 * x44 * x45,
            x25 * x45 * x46,
            x33 * x41 * x48,
            x50 * x8,
            x25 * x48 * x51,
            x35 * x41 * x46,
            x35 * x44 * x51,
            x38 * x39 * x53,
            x42 * x54,
            x22 * x45 * x55,
            x30 * x45 * x57,
            x26 * x41 * x55,
            x32 * x51 * x54,
            x26 * x53 * x57,
            x30 * x41 * x59,
            x22 * x53 * x59,
            x60 * x8,
            x33 * x62 * x63,
            x13 * x66,
            x13 * x25 * x67,
            x10 * x68,
            x49 * (x0 * (3 * x29 + 2 * x43 + x61) + x26 * x65),
            x25 * x68,
            x10 * x35 * x67,
            x35 * x66,
            x38 * x5 * x62,
            x46 * x54 * x63,
            x13 * x44 * x69,
            x13 * x57 * x70,
            x10 * x48 * x69,
            x50 * x54,
            x48 * x5 * x57,
            x10 * x59 * x70,
            x44 * x5 * x59,
            x39 * x60,
            x30 * x63 * x72,
            x13 * x22 * x73,
            x13 * x76,
            x10 * x26 * x73,
            x32 * x5 * x72,
            x26 * x76,
            x10 * x77,
            x22 * x77,
            x52 * (x0 * (3 * x36 + 2 * x56 + x71) + x35 * x75),
        ]
    )


def quadrupole3d_12(a, A, b, B, C):
    """Cartesian 3D (pd) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - C[0]
    x5 = a * b * x0
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x8 = x6 * x7
    x9 = x4 ** 2 * x8
    x10 = x3 * x8
    x11 = -x1 - B[0]
    x12 = x11 * x6
    x13 = x12 * x7
    x14 = x13 * x4
    x15 = 3 * x10 + 2 * x14
    x16 = x3 * (x15 + x9)
    x17 = x4 * x8
    x18 = x17 * x3
    x19 = 2 * x18
    x20 = x10 + x9
    x21 = x11 * x20
    x22 = x19 + x21
    x23 = x11 * x22 + x16
    x24 = x3 * (x13 + x17)
    x25 = x10 + x14
    x26 = x11 * x25
    x27 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x28 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x29 = numpy.pi * x0 * x28
    x30 = x27 * x29
    x31 = -x0 * (a * A[1] + b * B[1])
    x32 = -x31 - B[1]
    x33 = x30 * (x16 + x2 * x22)
    x34 = -x0 * (a * A[2] + b * B[2])
    x35 = -x34 - B[2]
    x36 = x19 + x2 * x20
    x37 = x27 * x7
    x38 = x3 * x37
    x39 = x32 ** 2 * x37
    x40 = x38 + x39
    x41 = x28 * x7
    x42 = x40 * x41
    x43 = x30 * x35
    x44 = x3 * x41
    x45 = x35 ** 2 * x41
    x46 = x44 + x45
    x47 = x37 * x46
    x48 = -x31 - A[1]
    x49 = x23 * x30
    x50 = x32 * x37
    x51 = x38 + x48 * x50
    x52 = 2 * x38
    x53 = x32 * x52 + x40 * x48
    x54 = x35 * x41
    x55 = -x34 - A[2]
    x56 = x30 * x55
    x57 = x44 + x54 * x55
    x58 = 2 * x44
    x59 = x35 * x58 + x46 * x55
    x60 = -x31 - C[1]
    x61 = x11 ** 2 * x8
    x62 = x24 + x26
    x63 = x30 * (x2 * x62 + x3 * (x15 + x61))
    x64 = x2 * x25 + x24
    x65 = x50 * x60
    x66 = x38 + x65
    x67 = x41 * x66
    x68 = x37 * x60
    x69 = x3 * (x50 + x68)
    x70 = x32 * x66
    x71 = x69 + x70
    x72 = x10 + x17 * x2
    x73 = x38 + x48 * x68
    x74 = x48 * x66 + x69
    x75 = 3 * x38 + 2 * x65
    x76 = x29 * x6
    x77 = x76 * (x3 * (x39 + x75) + x48 * x71)
    x78 = x35 * x76
    x79 = x4 * x76
    x80 = numpy.pi * x0 * x27
    x81 = x6 * x80
    x82 = x4 * x81
    x83 = -x34 - C[2]
    x84 = x30 * x83
    x85 = x54 * x83
    x86 = x44 + x85
    x87 = x37 * x86
    x88 = x41 * x83
    x89 = x3 * (x54 + x88)
    x90 = x35 * x86
    x91 = x89 + x90
    x92 = x44 + x55 * x88
    x93 = x55 * x86 + x89
    x94 = x32 * x81
    x95 = 3 * x44 + 2 * x85
    x96 = x81 * (x3 * (x45 + x95) + x55 * x91)
    x97 = x10 + x61
    x98 = 2 * x13 * x3 + x2 * x97
    x99 = x37 * x60 ** 2
    x100 = x38 + x99
    x101 = x100 * x41
    x102 = x52 * x60
    x103 = x100 * x32
    x104 = x102 + x103
    x105 = x10 + x13 * x2
    x106 = x3 * (x75 + x99)
    x107 = x104 * x32 + x106
    x108 = x107 * x76
    x109 = x46 * x8
    x110 = x100 * x48 + x102
    x111 = x104 * x48 + x106
    x112 = x12 * x29
    x113 = x8 * x86
    x114 = x12 * x80
    x115 = x41 * x83 ** 2
    x116 = x115 + x44
    x117 = x116 * x37
    x118 = x58 * x83
    x119 = x116 * x35
    x120 = x118 + x119
    x121 = x116 * x8
    x122 = x3 * (x115 + x95)
    x123 = x120 * x35 + x122
    x124 = x123 * x81
    x125 = x116 * x55 + x118
    x126 = x120 * x55 + x122

    # 108 item(s)
    return numpy.array(
        [
            x30 * (x2 * x23 + x3 * (4 * x18 + 2 * x21 + 2 * x24 + 2 * x26)),
            x32 * x33,
            x33 * x35,
            x36 * x42,
            x32 * x36 * x43,
            x36 * x47,
            x48 * x49,
            x22 * x41 * x51,
            x22 * x43 * x48,
            x20 * x41 * x53,
            x20 * x51 * x54,
            x20 * x47 * x48,
            x49 * x55,
            x22 * x32 * x56,
            x22 * x37 * x57,
            x20 * x42 * x55,
            x20 * x50 * x57,
            x20 * x37 * x59,
            x60 * x63,
            x64 * x67,
            x43 * x60 * x64,
            x41 * x71 * x72,
            x54 * x66 * x72,
            x46 * x68 * x72,
            x41 * x62 * x73,
            x25 * x41 * x74,
            x25 * x54 * x73,
            x4 * x77,
            x4 * x74 * x78,
            x17 * x46 * x73,
            x56 * x60 * x62,
            x25 * x55 * x67,
            x25 * x57 * x68,
            x55 * x71 * x79,
            x17 * x57 * x66,
            x59 * x60 * x82,
            x63 * x83,
            x32 * x64 * x84,
            x64 * x87,
            x40 * x72 * x88,
            x50 * x72 * x86,
            x37 * x72 * x91,
            x48 * x62 * x84,
            x25 * x51 * x88,
            x25 * x48 * x87,
            x53 * x79 * x83,
            x17 * x51 * x86,
            x48 * x82 * x91,
            x37 * x62 * x92,
            x25 * x50 * x92,
            x25 * x37 * x93,
            x17 * x40 * x92,
            x4 * x93 * x94,
            x4 * x96,
            x101 * x98,
            x104 * x105 * x41,
            x100 * x105 * x54,
            x108 * x2,
            x104 * x2 * x78,
            x100 * x109 * x2,
            x110 * x41 * x97,
            x111 * x112,
            x110 * x112 * x35,
            x76 * (x107 * x48 + x3 * (2 * x103 + 4 * x38 * x60 + 2 * x69 + 2 * x70)),
            x111 * x78,
            x109 * x110,
            x101 * x55 * x97,
            x104 * x112 * x55,
            x100 * x13 * x57,
            x108 * x55,
            x104 * x57 * x8,
            x100 * x59 * x8,
            x60 * x84 * x98,
            x105 * x66 * x88,
            x105 * x68 * x86,
            x2 * x71 * x76 * x83,
            x113 * x2 * x66,
            x2 * x60 * x81 * x91,
            x73 * x88 * x97,
            x112 * x74 * x83,
            x13 * x73 * x86,
            x77 * x83,
            x113 * x74,
            x73 * x8 * x91,
            x68 * x92 * x97,
            x13 * x66 * x92,
            x114 * x60 * x93,
            x71 * x8 * x92,
            x66 * x8 * x93,
            x60 * x96,
            x117 * x98,
            x105 * x116 * x50,
            x105 * x120 * x37,
            x121 * x2 * x40,
            x120 * x2 * x94,
            x124 * x2,
            x117 * x48 * x97,
            x116 * x13 * x51,
            x114 * x120 * x48,
            x121 * x53,
            x120 * x51 * x8,
            x124 * x48,
            x125 * x37 * x97,
            x114 * x125 * x32,
            x114 * x126,
            x125 * x40 * x8,
            x126 * x94,
            x81 * (x123 * x55 + x3 * (2 * x119 + 4 * x44 * x83 + 2 * x89 + 2 * x90)),
        ]
    )


def quadrupole3d_13(a, A, b, B, C):
    """Cartesian 3D (pf) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3 ** 2 * x7
    x9 = x0 * x7
    x10 = 3 * x9
    x11 = -x2 - C[0]
    x12 = x3 * x7
    x13 = x11 * x12
    x14 = x10 + 2 * x13
    x15 = x0 * (x14 + x8)
    x16 = x11 ** 2 * x7
    x17 = x0 * (x14 + x16)
    x18 = x11 * x7
    x19 = x0 * x18
    x20 = 2 * x19
    x21 = x16 + x9
    x22 = x21 * x3
    x23 = x20 + x22
    x24 = x23 * x3
    x25 = x0 * (x12 + x18)
    x26 = x13 + x9
    x27 = x26 * x3
    x28 = x25 + x27
    x29 = x28 * x3
    x30 = -x2 - A[0]
    x31 = x17 + x24
    x32 = x0 * (4 * x19 + 2 * x22 + 2 * x25 + 2 * x27)
    x33 = x3 * x31 + x32
    x34 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x35 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x36 = numpy.pi * x1 * x35
    x37 = x34 * x36
    x38 = -x1 * (a * A[1] + b * B[1])
    x39 = -x38 - B[1]
    x40 = x37 * (x30 * x31 + x32)
    x41 = -x1 * (a * A[2] + b * B[2])
    x42 = -x41 - B[2]
    x43 = x17 + x23 * x30
    x44 = x34 * x6
    x45 = x0 * x44
    x46 = x39 ** 2 * x44
    x47 = x45 + x46
    x48 = x35 * x6
    x49 = x47 * x48
    x50 = x37 * x42
    x51 = x0 * x48
    x52 = x42 ** 2 * x48
    x53 = x51 + x52
    x54 = x44 * x53
    x55 = x20 + x21 * x30
    x56 = 2 * x45
    x57 = x39 * x56
    x58 = x39 * x47 + x57
    x59 = x48 * x58
    x60 = x42 * x48
    x61 = x39 * x44
    x62 = 2 * x51
    x63 = x42 * x62
    x64 = x42 * x53 + x63
    x65 = x44 * x64
    x66 = -x38 - A[1]
    x67 = x33 * x37
    x68 = x45 + x61 * x66
    x69 = x47 * x66 + x57
    x70 = 3 * x45
    x71 = x0 * (3 * x46 + x70) + x58 * x66
    x72 = -x41 - A[2]
    x73 = x37 * x72
    x74 = x51 + x60 * x72
    x75 = x53 * x72 + x63
    x76 = 3 * x51
    x77 = x0 * (3 * x52 + x76) + x64 * x72
    x78 = -x38 - C[1]
    x79 = 2 * x0 * x12
    x80 = x8 + x9
    x81 = x3 * x80 + x79
    x82 = x15 + x29
    x83 = x37 * (x0 * (3 * x25 + 3 * x27 + x81) + x30 * x82)
    x84 = x15 + x28 * x30
    x85 = x61 * x78
    x86 = x45 + x85
    x87 = x48 * x86
    x88 = x25 + x26 * x30
    x89 = x44 * x78
    x90 = x0 * (x61 + x89)
    x91 = x39 * x86
    x92 = x90 + x91
    x93 = x48 * x92
    x94 = x70 + 2 * x85
    x95 = x0 * (x46 + x94)
    x96 = x39 * x92
    x97 = x95 + x96
    x98 = x18 * x30 + x9
    x99 = x45 + x66 * x89
    x100 = x66 * x86 + x90
    x101 = x66 * x92 + x95
    x102 = x36 * x5
    x103 = x102 * (x0 * (x58 + 3 * x90 + 3 * x91) + x66 * x97)
    x104 = x102 * x11
    x105 = numpy.pi * x1 * x34 * x5
    x106 = x105 * x11
    x107 = -x41 - C[2]
    x108 = x107 * x37
    x109 = x107 * x60
    x110 = x109 + x51
    x111 = x110 * x44
    x112 = x107 * x48
    x113 = x0 * (x112 + x60)
    x114 = x110 * x42
    x115 = x113 + x114
    x116 = x115 * x44
    x117 = 2 * x109 + x76
    x118 = x0 * (x117 + x52)
    x119 = x115 * x42
    x120 = x118 + x119
    x121 = x112 * x72 + x51
    x122 = x110 * x72 + x113
    x123 = x115 * x72 + x118
    x124 = x105 * (x0 * (3 * x113 + 3 * x114 + x64) + x120 * x72)
    x125 = x0 * (x10 + 3 * x8) + x30 * x81
    x126 = x44 * x78 ** 2
    x127 = x126 + x45
    x128 = x127 * x48
    x129 = x30 * x80 + x79
    x130 = x56 * x78
    x131 = x127 * x39
    x132 = x130 + x131
    x133 = x132 * x48
    x134 = x0 * (x126 + x94)
    x135 = x132 * x39
    x136 = x134 + x135
    x137 = x12 * x30 + x9
    x138 = x0 * (2 * x131 + 4 * x45 * x78 + 2 * x90 + 2 * x91)
    x139 = x136 * x39 + x138
    x140 = x102 * x139
    x141 = x102 * x42
    x142 = x53 * x7
    x143 = x64 * x7
    x144 = x127 * x66 + x130
    x145 = x132 * x66 + x134
    x146 = x102 * (x136 * x66 + x138)
    x147 = x102 * x3
    x148 = x110 * x7
    x149 = x115 * x7
    x150 = x105 * x78
    x151 = x107 ** 2 * x48
    x152 = x151 + x51
    x153 = x152 * x44
    x154 = x107 * x62
    x155 = x152 * x42
    x156 = x154 + x155
    x157 = x156 * x44
    x158 = x0 * (x117 + x151)
    x159 = x156 * x42
    x160 = x158 + x159
    x161 = x152 * x7
    x162 = x156 * x7
    x163 = x105 * x160
    x164 = x0 * (4 * x107 * x51 + 2 * x113 + 2 * x114 + 2 * x155)
    x165 = x160 * x42 + x164
    x166 = x105 * x165
    x167 = x152 * x72 + x154
    x168 = x156 * x72 + x158
    x169 = x105 * (x160 * x72 + x164)

    # 180 item(s)
    return numpy.array(
        [
            x37 * (x0 * (2 * x15 + 3 * x17 + 3 * x24 + 2 * x29) + x30 * x33),
            x39 * x40,
            x40 * x42,
            x43 * x49,
            x39 * x43 * x50,
            x43 * x54,
            x55 * x59,
            x47 * x55 * x60,
            x53 * x55 * x61,
            x55 * x65,
            x66 * x67,
            x31 * x48 * x68,
            x31 * x50 * x66,
            x23 * x48 * x69,
            x23 * x60 * x68,
            x23 * x54 * x66,
            x21 * x48 * x71,
            x21 * x60 * x69,
            x21 * x53 * x68,
            x21 * x65 * x66,
            x67 * x72,
            x31 * x39 * x73,
            x31 * x44 * x74,
            x23 * x49 * x72,
            x23 * x61 * x74,
            x23 * x44 * x75,
            x21 * x59 * x72,
            x21 * x47 * x74,
            x21 * x61 * x75,
            x21 * x44 * x77,
            x78 * x83,
            x84 * x87,
            x50 * x78 * x84,
            x88 * x93,
            x60 * x86 * x88,
            x53 * x88 * x89,
            x48 * x97 * x98,
            x60 * x92 * x98,
            x53 * x86 * x98,
            x64 * x89 * x98,
            x48 * x82 * x99,
            x100 * x28 * x48,
            x28 * x60 * x99,
            x101 * x26 * x48,
            x100 * x26 * x60,
            x26 * x53 * x99,
            x103 * x11,
            x101 * x104 * x42,
            x100 * x18 * x53,
            x18 * x64 * x99,
            x73 * x78 * x82,
            x28 * x72 * x87,
            x28 * x74 * x89,
            x26 * x72 * x93,
            x26 * x74 * x86,
            x26 * x75 * x89,
            x104 * x72 * x97,
            x18 * x74 * x92,
            x18 * x75 * x86,
            x106 * x77 * x78,
            x107 * x83,
            x108 * x39 * x84,
            x111 * x84,
            x112 * x47 * x88,
            x110 * x61 * x88,
            x116 * x88,
            x112 * x58 * x98,
            x110 * x47 * x98,
            x115 * x61 * x98,
            x120 * x44 * x98,
            x108 * x66 * x82,
            x112 * x28 * x68,
            x111 * x28 * x66,
            x112 * x26 * x69,
            x110 * x26 * x68,
            x116 * x26 * x66,
            x104 * x107 * x71,
            x110 * x18 * x69,
            x115 * x18 * x68,
            x106 * x120 * x66,
            x121 * x44 * x82,
            x121 * x28 * x61,
            x122 * x28 * x44,
            x121 * x26 * x47,
            x122 * x26 * x61,
            x123 * x26 * x44,
            x121 * x18 * x58,
            x122 * x18 * x47,
            x106 * x123 * x39,
            x11 * x124,
            x125 * x128,
            x129 * x133,
            x127 * x129 * x60,
            x136 * x137 * x48,
            x132 * x137 * x60,
            x127 * x137 * x53,
            x140 * x30,
            x136 * x141 * x30,
            x132 * x142 * x30,
            x127 * x143 * x30,
            x144 * x48 * x81,
            x145 * x48 * x80,
            x144 * x60 * x80,
            x146 * x3,
            x141 * x145 * x3,
            x12 * x144 * x53,
            x102 * (x0 * (3 * x134 + 3 * x135 + 2 * x95 + 2 * x96) + x139 * x66),
            x146 * x42,
            x142 * x145,
            x143 * x144,
            x128 * x72 * x81,
            x133 * x72 * x80,
            x127 * x74 * x80,
            x136 * x147 * x72,
            x12 * x132 * x74,
            x12 * x127 * x75,
            x140 * x72,
            x136 * x7 * x74,
            x132 * x7 * x75,
            x127 * x7 * x77,
            x108 * x125 * x78,
            x112 * x129 * x86,
            x110 * x129 * x89,
            x112 * x137 * x92,
            x110 * x137 * x86,
            x115 * x137 * x89,
            x102 * x107 * x30 * x97,
            x148 * x30 * x92,
            x149 * x30 * x86,
            x120 * x150 * x30,
            x112 * x81 * x99,
            x100 * x112 * x80,
            x110 * x80 * x99,
            x101 * x107 * x147,
            x100 * x110 * x12,
            x115 * x12 * x99,
            x103 * x107,
            x101 * x148,
            x100 * x149,
            x120 * x7 * x99,
            x121 * x81 * x89,
            x121 * x80 * x86,
            x122 * x80 * x89,
            x12 * x121 * x92,
            x12 * x122 * x86,
            x123 * x150 * x3,
            x121 * x7 * x97,
            x122 * x7 * x92,
            x123 * x7 * x86,
            x124 * x78,
            x125 * x153,
            x129 * x152 * x61,
            x129 * x157,
            x137 * x152 * x47,
            x137 * x156 * x61,
            x137 * x160 * x44,
            x161 * x30 * x58,
            x162 * x30 * x47,
            x163 * x30 * x39,
            x166 * x30,
            x153 * x66 * x81,
            x152 * x68 * x80,
            x157 * x66 * x80,
            x12 * x152 * x69,
            x12 * x156 * x68,
            x163 * x3 * x66,
            x161 * x71,
            x162 * x69,
            x160 * x68 * x7,
            x166 * x66,
            x167 * x44 * x81,
            x167 * x61 * x80,
            x168 * x44 * x80,
            x12 * x167 * x47,
            x105 * x168 * x3 * x39,
            x169 * x3,
            x167 * x58 * x7,
            x168 * x47 * x7,
            x169 * x39,
            x105 * (x0 * (2 * x118 + 2 * x119 + 3 * x158 + 3 * x159) + x165 * x72),
        ]
    )


def quadrupole3d_14(a, A, b, B, C):
    """Cartesian 3D (pg) quadrupole moment integrals.
    The origin is at C.

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
    x9 = x4 ** 2 * x8
    x10 = x3 * x8
    x11 = 3 * x10
    x12 = -x1 - C[0]
    x13 = x4 * x8
    x14 = x12 * x13
    x15 = x11 + 2 * x14
    x16 = x3 * (x15 + x9)
    x17 = x12 ** 2 * x8
    x18 = x3 * (x15 + x17)
    x19 = 2 * x10
    x20 = x12 * x19
    x21 = x10 + x17
    x22 = x21 * x4
    x23 = x20 + x22
    x24 = x23 * x4
    x25 = x12 * x8
    x26 = x3 * (x13 + x25)
    x27 = x10 + x14
    x28 = x27 * x4
    x29 = x26 + x28
    x30 = x29 * x4
    x31 = x3 * (2 * x16 + 3 * x18 + 3 * x24 + 2 * x30)
    x32 = x18 + x24
    x33 = x32 * x4
    x34 = x3 * (4 * x10 * x12 + 2 * x22 + 2 * x26 + 2 * x28)
    x35 = x33 + x34
    x36 = x31 + x35 * x4
    x37 = x19 * x4
    x38 = x10 + x9
    x39 = x38 * x4
    x40 = x37 + x39
    x41 = x3 * (3 * x26 + 3 * x28 + x40)
    x42 = x16 + x30
    x43 = x4 * x42
    x44 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x45 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x46 = numpy.pi * x0 * x45
    x47 = x44 * x46
    x48 = -x0 * (a * A[1] + b * B[1])
    x49 = -x48 - B[1]
    x50 = x47 * (x2 * x35 + x31)
    x51 = -x0 * (a * A[2] + b * B[2])
    x52 = -x51 - B[2]
    x53 = x2 * x32 + x34
    x54 = x44 * x7
    x55 = x3 * x54
    x56 = x49 ** 2 * x54
    x57 = x55 + x56
    x58 = x45 * x7
    x59 = x57 * x58
    x60 = x47 * x52
    x61 = x3 * x58
    x62 = x52 ** 2 * x58
    x63 = x61 + x62
    x64 = x54 * x63
    x65 = x18 + x2 * x23
    x66 = x49 * x55
    x67 = 2 * x66
    x68 = x49 * x57
    x69 = x67 + x68
    x70 = x58 * x69
    x71 = x52 * x58
    x72 = x49 * x54
    x73 = x52 * x61
    x74 = 2 * x73
    x75 = x52 * x63
    x76 = x74 + x75
    x77 = x54 * x76
    x78 = x2 * x21 + x20
    x79 = 3 * x55
    x80 = x3 * (3 * x56 + x79)
    x81 = x49 * x69 + x80
    x82 = x58 * x81
    x83 = 3 * x61
    x84 = x3 * (3 * x62 + x83)
    x85 = x52 * x76 + x84
    x86 = x54 * x85
    x87 = -x48 - A[1]
    x88 = x36 * x47
    x89 = x55 + x72 * x87
    x90 = x57 * x87 + x67
    x91 = x69 * x87 + x80
    x92 = x3 * (8 * x66 + 4 * x68) + x81 * x87
    x93 = -x51 - A[2]
    x94 = x47 * x93
    x95 = x61 + x71 * x93
    x96 = x63 * x93 + x74
    x97 = x76 * x93 + x84
    x98 = x3 * (8 * x73 + 4 * x75) + x85 * x93
    x99 = -x48 - C[1]
    x100 = x3 * (x11 + 3 * x9)
    x101 = x100 + x4 * x40
    x102 = x41 + x43
    x103 = x47 * (x102 * x2 + x3 * (x101 + 4 * x16 + 4 * x30))
    x104 = x2 * x42 + x41
    x105 = x72 * x99
    x106 = x105 + x55
    x107 = x106 * x58
    x108 = x16 + x2 * x29
    x109 = x54 * x99
    x110 = x3 * (x109 + x72)
    x111 = x106 * x49
    x112 = x110 + x111
    x113 = x112 * x58
    x114 = x2 * x27 + x26
    x115 = 2 * x105 + x79
    x116 = x3 * (x115 + x56)
    x117 = x112 * x49
    x118 = x116 + x117
    x119 = x118 * x58
    x120 = x3 * (3 * x110 + 3 * x111 + x69)
    x121 = x118 * x49
    x122 = x120 + x121
    x123 = x10 + x2 * x25
    x124 = x109 * x87 + x55
    x125 = x106 * x87 + x110
    x126 = x112 * x87 + x116
    x127 = x118 * x87 + x120
    x128 = x46 * x6
    x129 = x128 * (x122 * x87 + x3 * (4 * x116 + 4 * x117 + x81))
    x130 = x12 * x128
    x131 = numpy.pi * x0 * x44 * x6
    x132 = x12 * x131
    x133 = -x51 - C[2]
    x134 = x133 * x47
    x135 = x133 * x71
    x136 = x135 + x61
    x137 = x136 * x54
    x138 = x133 * x58
    x139 = x3 * (x138 + x71)
    x140 = x136 * x52
    x141 = x139 + x140
    x142 = x141 * x54
    x143 = 2 * x135 + x83
    x144 = x3 * (x143 + x62)
    x145 = x141 * x52
    x146 = x144 + x145
    x147 = x146 * x54
    x148 = x3 * (3 * x139 + 3 * x140 + x76)
    x149 = x146 * x52
    x150 = x148 + x149
    x151 = x138 * x93 + x61
    x152 = x136 * x93 + x139
    x153 = x141 * x93 + x144
    x154 = x146 * x93 + x148
    x155 = x131 * (x150 * x93 + x3 * (4 * x144 + 4 * x145 + x85))
    x156 = x101 * x2 + x3 * (8 * x10 * x4 + 4 * x39)
    x157 = x54 * x99 ** 2
    x158 = x157 + x55
    x159 = x158 * x58
    x160 = x100 + x2 * x40
    x161 = x55 * x99
    x162 = 2 * x161
    x163 = x158 * x49
    x164 = x162 + x163
    x165 = x164 * x58
    x166 = x2 * x38 + x37
    x167 = x3 * (x115 + x157)
    x168 = x164 * x49
    x169 = x167 + x168
    x170 = x169 * x58
    x171 = x169 * x49
    x172 = x3 * (2 * x110 + 2 * x111 + 4 * x161 + 2 * x163)
    x173 = x171 + x172
    x174 = x10 + x13 * x2
    x175 = x3 * (2 * x116 + 2 * x117 + 3 * x167 + 3 * x168)
    x176 = x173 * x49 + x175
    x177 = x128 * x176
    x178 = x128 * x52
    x179 = x63 * x8
    x180 = x76 * x8
    x181 = x8 * x85
    x182 = x158 * x87 + x162
    x183 = x164 * x87 + x167
    x184 = x169 * x87 + x172
    x185 = x128 * (x173 * x87 + x175)
    x186 = x128 * x4
    x187 = x136 * x8
    x188 = x141 * x8
    x189 = x146 * x8
    x190 = x131 * x99
    x191 = x133 ** 2 * x58
    x192 = x191 + x61
    x193 = x192 * x54
    x194 = x133 * x61
    x195 = 2 * x194
    x196 = x192 * x52
    x197 = x195 + x196
    x198 = x197 * x54
    x199 = x3 * (x143 + x191)
    x200 = x197 * x52
    x201 = x199 + x200
    x202 = x201 * x54
    x203 = x201 * x52
    x204 = x3 * (2 * x139 + 2 * x140 + 4 * x194 + 2 * x196)
    x205 = x203 + x204
    x206 = x192 * x8
    x207 = x197 * x8
    x208 = x201 * x8
    x209 = x131 * x205
    x210 = x3 * (2 * x144 + 2 * x145 + 3 * x199 + 3 * x200)
    x211 = x205 * x52 + x210
    x212 = x131 * x211
    x213 = x192 * x93 + x195
    x214 = x197 * x93 + x199
    x215 = x201 * x93 + x204
    x216 = x131 * (x205 * x93 + x210)

    # 270 item(s)
    return numpy.array(
        [
            x47 * (x2 * x36 + x3 * (4 * x33 + 4 * x34 + 2 * x41 + 2 * x43)),
            x49 * x50,
            x50 * x52,
            x53 * x59,
            x49 * x53 * x60,
            x53 * x64,
            x65 * x70,
            x57 * x65 * x71,
            x63 * x65 * x72,
            x65 * x77,
            x78 * x82,
            x69 * x71 * x78,
            x57 * x63 * x78,
            x72 * x76 * x78,
            x78 * x86,
            x87 * x88,
            x35 * x58 * x89,
            x35 * x60 * x87,
            x32 * x58 * x90,
            x32 * x71 * x89,
            x32 * x64 * x87,
            x23 * x58 * x91,
            x23 * x71 * x90,
            x23 * x63 * x89,
            x23 * x77 * x87,
            x21 * x58 * x92,
            x21 * x71 * x91,
            x21 * x63 * x90,
            x21 * x76 * x89,
            x21 * x86 * x87,
            x88 * x93,
            x35 * x49 * x94,
            x35 * x54 * x95,
            x32 * x59 * x93,
            x32 * x72 * x95,
            x32 * x54 * x96,
            x23 * x70 * x93,
            x23 * x57 * x95,
            x23 * x72 * x96,
            x23 * x54 * x97,
            x21 * x82 * x93,
            x21 * x69 * x95,
            x21 * x57 * x96,
            x21 * x72 * x97,
            x21 * x54 * x98,
            x103 * x99,
            x104 * x107,
            x104 * x60 * x99,
            x108 * x113,
            x106 * x108 * x71,
            x108 * x109 * x63,
            x114 * x119,
            x112 * x114 * x71,
            x106 * x114 * x63,
            x109 * x114 * x76,
            x122 * x123 * x58,
            x118 * x123 * x71,
            x112 * x123 * x63,
            x106 * x123 * x76,
            x109 * x123 * x85,
            x102 * x124 * x58,
            x125 * x42 * x58,
            x124 * x42 * x71,
            x126 * x29 * x58,
            x125 * x29 * x71,
            x124 * x29 * x63,
            x127 * x27 * x58,
            x126 * x27 * x71,
            x125 * x27 * x63,
            x124 * x27 * x76,
            x12 * x129,
            x127 * x130 * x52,
            x126 * x25 * x63,
            x125 * x25 * x76,
            x124 * x25 * x85,
            x102 * x94 * x99,
            x107 * x42 * x93,
            x109 * x42 * x95,
            x113 * x29 * x93,
            x106 * x29 * x95,
            x109 * x29 * x96,
            x119 * x27 * x93,
            x112 * x27 * x95,
            x106 * x27 * x96,
            x109 * x27 * x97,
            x122 * x130 * x93,
            x118 * x25 * x95,
            x112 * x25 * x96,
            x106 * x25 * x97,
            x132 * x98 * x99,
            x103 * x133,
            x104 * x134 * x49,
            x104 * x137,
            x108 * x138 * x57,
            x108 * x136 * x72,
            x108 * x142,
            x114 * x138 * x69,
            x114 * x136 * x57,
            x114 * x141 * x72,
            x114 * x147,
            x123 * x138 * x81,
            x123 * x136 * x69,
            x123 * x141 * x57,
            x123 * x146 * x72,
            x123 * x150 * x54,
            x102 * x134 * x87,
            x138 * x42 * x89,
            x137 * x42 * x87,
            x138 * x29 * x90,
            x136 * x29 * x89,
            x142 * x29 * x87,
            x138 * x27 * x91,
            x136 * x27 * x90,
            x141 * x27 * x89,
            x147 * x27 * x87,
            x130 * x133 * x92,
            x136 * x25 * x91,
            x141 * x25 * x90,
            x146 * x25 * x89,
            x132 * x150 * x87,
            x102 * x151 * x54,
            x151 * x42 * x72,
            x152 * x42 * x54,
            x151 * x29 * x57,
            x152 * x29 * x72,
            x153 * x29 * x54,
            x151 * x27 * x69,
            x152 * x27 * x57,
            x153 * x27 * x72,
            x154 * x27 * x54,
            x151 * x25 * x81,
            x152 * x25 * x69,
            x153 * x25 * x57,
            x132 * x154 * x49,
            x12 * x155,
            x156 * x159,
            x160 * x165,
            x158 * x160 * x71,
            x166 * x170,
            x164 * x166 * x71,
            x158 * x166 * x63,
            x173 * x174 * x58,
            x169 * x174 * x71,
            x164 * x174 * x63,
            x158 * x174 * x76,
            x177 * x2,
            x173 * x178 * x2,
            x169 * x179 * x2,
            x164 * x180 * x2,
            x158 * x181 * x2,
            x101 * x182 * x58,
            x183 * x40 * x58,
            x182 * x40 * x71,
            x184 * x38 * x58,
            x183 * x38 * x71,
            x182 * x38 * x63,
            x185 * x4,
            x178 * x184 * x4,
            x13 * x183 * x63,
            x13 * x182 * x76,
            x128 * (x176 * x87 + x3 * (2 * x120 + 2 * x121 + 4 * x171 + 4 * x172)),
            x185 * x52,
            x179 * x184,
            x180 * x183,
            x181 * x182,
            x101 * x159 * x93,
            x165 * x40 * x93,
            x158 * x40 * x95,
            x170 * x38 * x93,
            x164 * x38 * x95,
            x158 * x38 * x96,
            x173 * x186 * x93,
            x13 * x169 * x95,
            x13 * x164 * x96,
            x13 * x158 * x97,
            x177 * x93,
            x173 * x8 * x95,
            x169 * x8 * x96,
            x164 * x8 * x97,
            x158 * x8 * x98,
            x134 * x156 * x99,
            x106 * x138 * x160,
            x109 * x136 * x160,
            x112 * x138 * x166,
            x106 * x136 * x166,
            x109 * x141 * x166,
            x118 * x138 * x174,
            x112 * x136 * x174,
            x106 * x141 * x174,
            x109 * x146 * x174,
            x122 * x128 * x133 * x2,
            x118 * x187 * x2,
            x112 * x188 * x2,
            x106 * x189 * x2,
            x150 * x190 * x2,
            x101 * x124 * x138,
            x125 * x138 * x40,
            x124 * x136 * x40,
            x126 * x138 * x38,
            x125 * x136 * x38,
            x124 * x141 * x38,
            x127 * x133 * x186,
            x126 * x13 * x136,
            x125 * x13 * x141,
            x124 * x13 * x146,
            x129 * x133,
            x127 * x187,
            x126 * x188,
            x125 * x189,
            x124 * x150 * x8,
            x101 * x109 * x151,
            x106 * x151 * x40,
            x109 * x152 * x40,
            x112 * x151 * x38,
            x106 * x152 * x38,
            x109 * x153 * x38,
            x118 * x13 * x151,
            x112 * x13 * x152,
            x106 * x13 * x153,
            x154 * x190 * x4,
            x122 * x151 * x8,
            x118 * x152 * x8,
            x112 * x153 * x8,
            x106 * x154 * x8,
            x155 * x99,
            x156 * x193,
            x160 * x192 * x72,
            x160 * x198,
            x166 * x192 * x57,
            x166 * x197 * x72,
            x166 * x202,
            x174 * x192 * x69,
            x174 * x197 * x57,
            x174 * x201 * x72,
            x174 * x205 * x54,
            x2 * x206 * x81,
            x2 * x207 * x69,
            x2 * x208 * x57,
            x2 * x209 * x49,
            x2 * x212,
            x101 * x193 * x87,
            x192 * x40 * x89,
            x198 * x40 * x87,
            x192 * x38 * x90,
            x197 * x38 * x89,
            x202 * x38 * x87,
            x13 * x192 * x91,
            x13 * x197 * x90,
            x13 * x201 * x89,
            x209 * x4 * x87,
            x206 * x92,
            x207 * x91,
            x208 * x90,
            x205 * x8 * x89,
            x212 * x87,
            x101 * x213 * x54,
            x213 * x40 * x72,
            x214 * x40 * x54,
            x213 * x38 * x57,
            x214 * x38 * x72,
            x215 * x38 * x54,
            x13 * x213 * x69,
            x13 * x214 * x57,
            x131 * x215 * x4 * x49,
            x216 * x4,
            x213 * x8 * x81,
            x214 * x69 * x8,
            x215 * x57 * x8,
            x216 * x49,
            x131 * (x211 * x93 + x3 * (2 * x148 + 2 * x149 + 4 * x203 + 4 * x204)),
        ]
    )


def quadrupole3d_20(a, A, b, B, C):
    """Cartesian 3D (ds) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = a * b * x1
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x5 = x3 * x4
    x6 = x0 * x5
    x7 = -x1 * (a * A[0] + b * B[0])
    x8 = -x7 - C[0]
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
    x20 = -x19 - A[1]
    x21 = x14 * x18
    x22 = -x1 * (a * A[2] + b * B[2])
    x23 = -x22 - A[2]
    x24 = x15 * x4
    x25 = x0 * x24
    x26 = x20 ** 2 * x24 + x25
    x27 = x16 * x4
    x28 = x18 * x23
    x29 = x0 * x27
    x30 = x23 ** 2 * x27 + x29
    x31 = -x19 - C[1]
    x32 = x12 + x6
    x33 = x18 * (x0 * (x11 + x5 * x8) + x10 * x32)
    x34 = x20 * x24
    x35 = x31 * x34
    x36 = x25 + x35
    x37 = x17 * x3
    x38 = x37 * (x0 * (x24 * x31 + x34) + x20 * x36)
    x39 = x37 * x8
    x40 = numpy.pi * x1 * x15 * x3
    x41 = x40 * x8
    x42 = -x22 - C[2]
    x43 = x18 * x42
    x44 = x23 * x27
    x45 = x42 * x44
    x46 = x29 + x45
    x47 = x40 * (x0 * (x27 * x42 + x44) + x23 * x46)
    x48 = x24 * x31 ** 2
    x49 = x25 + x48
    x50 = x10 ** 2 * x5 + x6
    x51 = x20 * x49 + 2 * x25 * x31
    x52 = x37 * x51
    x53 = x10 * x37
    x54 = x10 * x40
    x55 = x27 * x42 ** 2
    x56 = x29 + x55
    x57 = x23 * x56 + 2 * x29 * x42
    x58 = x40 * x57

    # 36 item(s)
    return numpy.array(
        [
            x18 * (x0 * (2 * x12 + 3 * x6 + x9) + x10 * x14),
            x20 * x21,
            x21 * x23,
            x13 * x26 * x27,
            x13 * x20 * x28,
            x13 * x24 * x30,
            x31 * x33,
            x27 * x32 * x36,
            x28 * x31 * x32,
            x38 * x8,
            x23 * x36 * x39,
            x30 * x31 * x41,
            x33 * x42,
            x20 * x32 * x43,
            x24 * x32 * x46,
            x26 * x39 * x42,
            x20 * x41 * x46,
            x47 * x8,
            x27 * x49 * x50,
            x10 * x52,
            x23 * x49 * x53,
            x37 * (x0 * (3 * x25 + 2 * x35 + x48) + x20 * x51),
            x23 * x52,
            x30 * x49 * x5,
            x31 * x43 * x50,
            x36 * x42 * x53,
            x31 * x46 * x54,
            x38 * x42,
            x36 * x46 * x5,
            x31 * x47,
            x24 * x50 * x56,
            x20 * x54 * x56,
            x10 * x58,
            x26 * x5 * x56,
            x20 * x58,
            x40 * (x0 * (3 * x29 + 2 * x45 + x55) + x23 * x57),
        ]
    )


def quadrupole3d_21(a, A, b, B, C):
    """Cartesian 3D (dp) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - C[0]
    x5 = -x1 - B[0]
    x6 = a * b * x0
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x9 = x7 * x8
    x10 = x5 * x9
    x11 = x10 * x4
    x12 = x4 ** 2 * x9
    x13 = x3 * x9
    x14 = 3 * x13
    x15 = x12 + x14
    x16 = x13 * x4
    x17 = 2 * x16
    x18 = x12 + x13
    x19 = x18 * x5
    x20 = x17 + x19
    x21 = x2 * x20 + x3 * (2 * x11 + x15)
    x22 = x4 * x9
    x23 = x3 * (x10 + x22)
    x24 = x18 * x2
    x25 = x11 + x13
    x26 = x2 * x25
    x27 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x28 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x29 = numpy.pi * x0 * x28
    x30 = x27 * x29
    x31 = -x0 * (a * A[1] + b * B[1])
    x32 = -x31 - B[1]
    x33 = x2 * x22
    x34 = x17 + x24
    x35 = x30 * (x2 * x34 + x3 * (x15 + 2 * x33))
    x36 = -x0 * (a * A[2] + b * B[2])
    x37 = -x36 - B[2]
    x38 = -x31 - A[1]
    x39 = x21 * x30
    x40 = x3 * x8
    x41 = x27 * x40
    x42 = x27 * x8
    x43 = x38 * x42
    x44 = x32 * x43
    x45 = x41 + x44
    x46 = x28 * x8
    x47 = x30 * x34
    x48 = -x36 - A[2]
    x49 = x28 * x40
    x50 = x46 * x48
    x51 = x37 * x50
    x52 = x49 + x51
    x53 = x38 ** 2 * x42 + x41
    x54 = x32 * x42
    x55 = x3 * (x43 + x54) + x38 * x45
    x56 = x37 * x46
    x57 = x30 * x48
    x58 = x46 * x48 ** 2 + x49
    x59 = x3 * (x50 + x56) + x48 * x52
    x60 = -x31 - C[1]
    x61 = x10 * x2
    x62 = x23 + x26
    x63 = x30 * (x2 * x62 + x3 * (x11 + x14 + x33 + x61))
    x64 = x54 * x60
    x65 = x41 + x64
    x66 = x2 * x9
    x67 = x13 + x33
    x68 = x2 * x67 + x3 * (x22 + x66)
    x69 = x30 * x68
    x70 = x43 * x60
    x71 = x41 + x70
    x72 = x42 * x60
    x73 = x3 * (x54 + x72)
    x74 = x38 * x65
    x75 = x73 + x74
    x76 = x3 * (x43 + x72) + x38 * x71
    x77 = 3 * x41
    x78 = x29 * x7
    x79 = x78 * (x3 * (x44 + x64 + x70 + x77) + x38 * x75)
    x80 = x4 * x78
    x81 = numpy.pi * x0 * x27 * x7
    x82 = x4 * x81
    x83 = -x36 - C[2]
    x84 = x56 * x83
    x85 = x49 + x84
    x86 = x30 * x83
    x87 = x46 * x83
    x88 = x50 * x83
    x89 = x49 + x88
    x90 = x3 * (x56 + x87)
    x91 = x48 * x85
    x92 = x90 + x91
    x93 = x3 * (x50 + x87) + x48 * x89
    x94 = 3 * x49
    x95 = x81 * (x3 * (x51 + x84 + x88 + x94) + x48 * x92)
    x96 = x42 * x60 ** 2
    x97 = x41 + x96
    x98 = x13 + x61
    x99 = x2 * x98 + x3 * (x10 + x66)
    x100 = x41 * x60
    x101 = 2 * x100
    x102 = x32 * x97
    x103 = x101 + x102
    x104 = x13 + x2 ** 2 * x9
    x105 = x38 * x97
    x106 = x101 + x105
    x107 = x77 + x96
    x108 = x103 * x38 + x3 * (x107 + 2 * x64)
    x109 = x108 * x78
    x110 = x2 * x78
    x111 = x78 * (x106 * x38 + x3 * (x107 + 2 * x70))
    x112 = x5 * x78
    x113 = x60 * x81
    x114 = x46 * x83 ** 2
    x115 = x114 + x49
    x116 = x49 * x83
    x117 = 2 * x116
    x118 = x115 * x37
    x119 = x117 + x118
    x120 = x2 * x81
    x121 = x115 * x48
    x122 = x117 + x121
    x123 = x114 + x94
    x124 = x119 * x48 + x3 * (x123 + 2 * x84)
    x125 = x124 * x81
    x126 = x81 * (x122 * x48 + x3 * (x123 + 2 * x88))

    # 108 item(s)
    return numpy.array(
        [
            x30 * (x2 * x21 + x3 * (4 * x16 + x19 + 2 * x23 + x24 + 2 * x26)),
            x32 * x35,
            x35 * x37,
            x38 * x39,
            x34 * x45 * x46,
            x37 * x38 * x47,
            x39 * x48,
            x32 * x47 * x48,
            x34 * x42 * x52,
            x20 * x46 * x53,
            x18 * x46 * x55,
            x18 * x53 * x56,
            x20 * x38 * x57,
            x18 * x45 * x50,
            x18 * x43 * x52,
            x20 * x42 * x58,
            x18 * x54 * x58,
            x18 * x42 * x59,
            x60 * x63,
            x46 * x65 * x68,
            x37 * x60 * x69,
            x46 * x62 * x71,
            x46 * x67 * x75,
            x56 * x67 * x71,
            x57 * x60 * x62,
            x50 * x65 * x67,
            x52 * x67 * x72,
            x25 * x46 * x76,
            x4 * x79,
            x37 * x76 * x80,
            x25 * x50 * x71,
            x48 * x75 * x80,
            x22 * x52 * x71,
            x25 * x58 * x72,
            x22 * x58 * x65,
            x59 * x60 * x82,
            x63 * x83,
            x32 * x69 * x83,
            x42 * x68 * x85,
            x38 * x62 * x86,
            x45 * x67 * x87,
            x43 * x67 * x85,
            x42 * x62 * x89,
            x54 * x67 * x89,
            x42 * x67 * x92,
            x25 * x53 * x87,
            x55 * x80 * x83,
            x22 * x53 * x85,
            x25 * x43 * x89,
            x22 * x45 * x89,
            x38 * x82 * x92,
            x25 * x42 * x93,
            x32 * x82 * x93,
            x4 * x95,
            x46 * x97 * x99,
            x103 * x104 * x46,
            x104 * x56 * x97,
            x106 * x46 * x98,
            x109 * x2,
            x106 * x110 * x37,
            x50 * x97 * x98,
            x103 * x110 * x48,
            x52 * x66 * x97,
            x111 * x5,
            x78 * (x108 * x38 + x3 * (4 * x100 + x102 + x105 + 2 * x73 + 2 * x74)),
            x111 * x37,
            x106 * x112 * x48,
            x109 * x48,
            x106 * x52 * x9,
            x10 * x58 * x97,
            x103 * x58 * x9,
            x59 * x9 * x97,
            x60 * x86 * x99,
            x104 * x65 * x87,
            x104 * x72 * x85,
            x71 * x87 * x98,
            x110 * x75 * x83,
            x66 * x71 * x85,
            x72 * x89 * x98,
            x65 * x66 * x89,
            x113 * x2 * x92,
            x112 * x76 * x83,
            x79 * x83,
            x76 * x85 * x9,
            x10 * x71 * x89,
            x75 * x89 * x9,
            x71 * x9 * x92,
            x113 * x5 * x93,
            x65 * x9 * x93,
            x60 * x95,
            x115 * x42 * x99,
            x104 * x115 * x54,
            x104 * x119 * x42,
            x115 * x43 * x98,
            x115 * x45 * x66,
            x119 * x120 * x38,
            x122 * x42 * x98,
            x120 * x122 * x32,
            x125 * x2,
            x10 * x115 * x53,
            x115 * x55 * x9,
            x119 * x53 * x9,
            x122 * x38 * x5 * x81,
            x122 * x45 * x9,
            x125 * x38,
            x126 * x5,
            x126 * x32,
            x81 * (x124 * x48 + x3 * (4 * x116 + x118 + x121 + 2 * x90 + 2 * x91)),
        ]
    )


def quadrupole3d_22(a, A, b, B, C):
    """Cartesian 3D (dd) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - C[0]
    x5 = a * b * x0
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x8 = x6 * x7
    x9 = x4 ** 2 * x8
    x10 = x3 * x8
    x11 = 3 * x10
    x12 = -x1 - B[0]
    x13 = x12 * x8
    x14 = x13 * x4
    x15 = x11 + 2 * x14
    x16 = x3 * (x15 + x9)
    x17 = x10 * x4
    x18 = 2 * x17
    x19 = x10 + x9
    x20 = x12 * x19
    x21 = x18 + x20
    x22 = x12 * x21
    x23 = x16 + x22
    x24 = x10 + x14
    x25 = x12 * x24
    x26 = x4 * x8
    x27 = x3 * (x13 + x26)
    x28 = 4 * x17 + 2 * x27
    x29 = x2 * x23 + x3 * (2 * x20 + 2 * x25 + x28)
    x30 = x12 ** 2 * x8
    x31 = x3 * (x15 + x30)
    x32 = x2 * x21
    x33 = x25 + x27
    x34 = x2 * x33
    x35 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x36 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x37 = numpy.pi * x0 * x36
    x38 = x35 * x37
    x39 = -x0 * (a * A[1] + b * B[1])
    x40 = -x39 - B[1]
    x41 = x16 + x32
    x42 = x19 * x2
    x43 = x2 * x24
    x44 = 2 * x43
    x45 = x38 * (x2 * x41 + x3 * (x20 + x28 + x42 + x44))
    x46 = -x0 * (a * A[2] + b * B[2])
    x47 = -x46 - B[2]
    x48 = x35 * x7
    x49 = x3 * x48
    x50 = x40 ** 2 * x48
    x51 = x49 + x50
    x52 = x2 * x26
    x53 = x18 + x42
    x54 = x2 * x53 + x3 * (x11 + 2 * x52 + x9)
    x55 = x36 * x7
    x56 = x38 * x47
    x57 = x3 * x55
    x58 = x47 ** 2 * x55
    x59 = x57 + x58
    x60 = -x39 - A[1]
    x61 = x29 * x38
    x62 = x48 * x60
    x63 = x40 * x62
    x64 = x49 + x63
    x65 = 2 * x49
    x66 = x40 * x65 + x51 * x60
    x67 = x47 * x55
    x68 = -x46 - A[2]
    x69 = x38 * x68
    x70 = x55 * x68
    x71 = x47 * x70
    x72 = x57 + x71
    x73 = x40 * x48
    x74 = 2 * x57
    x75 = x47 * x74 + x59 * x68
    x76 = x48 * x60 ** 2 + x49
    x77 = x3 * (x62 + x73) + x60 * x64
    x78 = 3 * x49
    x79 = x50 + x78
    x80 = x3 * (2 * x63 + x79) + x60 * x66
    x81 = x55 * x68 ** 2 + x57
    x82 = x3 * (x67 + x70) + x68 * x72
    x83 = 3 * x57
    x84 = x58 + x83
    x85 = x3 * (2 * x71 + x84) + x68 * x75
    x86 = -x39 - C[1]
    x87 = x31 + x34
    x88 = x10 + x30
    x89 = 2 * x10 * x12 + x2 * x88
    x90 = x38 * (x2 * x87 + x3 * (x25 + 3 * x27 + x44 + x89))
    x91 = x73 * x86
    x92 = x49 + x91
    x93 = x13 * x2
    x94 = x27 + x43
    x95 = x2 * x94 + x3 * (x11 + x14 + x52 + x93)
    x96 = x48 * x86
    x97 = x3 * (x73 + x96)
    x98 = x40 * x92
    x99 = x97 + x98
    x100 = x2 * x8
    x101 = x10 + x52
    x102 = x101 * x2 + x3 * (x100 + x26)
    x103 = x62 * x86
    x104 = x103 + x49
    x105 = x60 * x92
    x106 = x105 + x97
    x107 = 2 * x91
    x108 = x3 * (x107 + x79)
    x109 = x60 * x99
    x110 = x108 + x109
    x111 = x104 * x60 + x3 * (x62 + x96)
    x112 = x106 * x60 + x3 * (x103 + x63 + x78 + x91)
    x113 = 2 * x105
    x114 = x37 * x6
    x115 = x114 * (x110 * x60 + x3 * (x113 + x66 + 3 * x97 + x98))
    x116 = x114 * x4
    x117 = numpy.pi * x0 * x35 * x6
    x118 = x117 * x4
    x119 = -x46 - C[2]
    x120 = x119 * x38
    x121 = x119 * x67
    x122 = x121 + x57
    x123 = x119 * x55
    x124 = x3 * (x123 + x67)
    x125 = x122 * x47
    x126 = x124 + x125
    x127 = x119 * x70
    x128 = x127 + x57
    x129 = x122 * x68
    x130 = x124 + x129
    x131 = 2 * x121
    x132 = x3 * (x131 + x84)
    x133 = x126 * x68
    x134 = x132 + x133
    x135 = x128 * x68 + x3 * (x123 + x70)
    x136 = x130 * x68 + x3 * (x121 + x127 + x71 + x83)
    x137 = 2 * x129
    x138 = x117 * (x134 * x68 + x3 * (3 * x124 + x125 + x137 + x75))
    x139 = x48 * x86 ** 2
    x140 = x139 + x49
    x141 = x2 * x89 + x3 * (x11 + x30 + 2 * x93)
    x142 = x65 * x86
    x143 = x140 * x40
    x144 = x142 + x143
    x145 = x10 + x93
    x146 = x145 * x2 + x3 * (x100 + x13)
    x147 = x139 + x78
    x148 = x3 * (x107 + x147)
    x149 = x144 * x40
    x150 = x148 + x149
    x151 = x10 + x2 ** 2 * x8
    x152 = x140 * x60
    x153 = x142 + x152
    x154 = x144 * x60
    x155 = x148 + x154
    x156 = 4 * x49 * x86 + 2 * x97
    x157 = x150 * x60 + x3 * (2 * x143 + x156 + 2 * x98)
    x158 = x114 * x157
    x159 = x114 * x2
    x160 = x153 * x60 + x3 * (2 * x103 + x147)
    x161 = x114 * (x155 * x60 + x3 * (x113 + x143 + x152 + x156))
    x162 = x114 * x12
    x163 = x117 * x86
    x164 = x119 ** 2 * x55
    x165 = x164 + x57
    x166 = x119 * x74
    x167 = x165 * x47
    x168 = x166 + x167
    x169 = x164 + x83
    x170 = x3 * (x131 + x169)
    x171 = x168 * x47
    x172 = x170 + x171
    x173 = x117 * x2
    x174 = x165 * x68
    x175 = x166 + x174
    x176 = x168 * x68
    x177 = x170 + x176
    x178 = 4 * x119 * x57 + 2 * x124
    x179 = x172 * x68 + x3 * (2 * x125 + 2 * x167 + x178)
    x180 = x117 * x179
    x181 = x117 * x12
    x182 = x175 * x68 + x3 * (2 * x127 + x169)
    x183 = x117 * (x177 * x68 + x3 * (x137 + x167 + x174 + x178))

    # 216 item(s)
    return numpy.array(
        [
            x38 * (x2 * x29 + x3 * (3 * x16 + x22 + 2 * x31 + 2 * x32 + 2 * x34)),
            x40 * x45,
            x45 * x47,
            x51 * x54 * x55,
            x40 * x54 * x56,
            x48 * x54 * x59,
            x60 * x61,
            x41 * x55 * x64,
            x41 * x56 * x60,
            x53 * x55 * x66,
            x53 * x64 * x67,
            x53 * x59 * x62,
            x61 * x68,
            x40 * x41 * x69,
            x41 * x48 * x72,
            x51 * x53 * x70,
            x53 * x72 * x73,
            x48 * x53 * x75,
            x23 * x55 * x76,
            x21 * x55 * x77,
            x21 * x67 * x76,
            x19 * x55 * x80,
            x19 * x67 * x77,
            x19 * x59 * x76,
            x23 * x60 * x69,
            x21 * x64 * x70,
            x21 * x62 * x72,
            x19 * x66 * x70,
            x19 * x64 * x72,
            x19 * x62 * x75,
            x23 * x48 * x81,
            x21 * x73 * x81,
            x21 * x48 * x82,
            x19 * x51 * x81,
            x19 * x73 * x82,
            x19 * x48 * x85,
            x86 * x90,
            x55 * x92 * x95,
            x56 * x86 * x95,
            x102 * x55 * x99,
            x102 * x67 * x92,
            x102 * x59 * x96,
            x104 * x55 * x87,
            x106 * x55 * x94,
            x104 * x67 * x94,
            x101 * x110 * x55,
            x101 * x106 * x67,
            x101 * x104 * x59,
            x69 * x86 * x87,
            x70 * x92 * x94,
            x72 * x94 * x96,
            x101 * x70 * x99,
            x101 * x72 * x92,
            x101 * x75 * x96,
            x111 * x33 * x55,
            x112 * x24 * x55,
            x111 * x24 * x67,
            x115 * x4,
            x112 * x116 * x47,
            x111 * x26 * x59,
            x104 * x33 * x70,
            x106 * x24 * x70,
            x104 * x24 * x72,
            x110 * x116 * x68,
            x106 * x26 * x72,
            x104 * x26 * x75,
            x33 * x81 * x96,
            x24 * x81 * x92,
            x24 * x82 * x96,
            x26 * x81 * x99,
            x26 * x82 * x92,
            x118 * x85 * x86,
            x119 * x90,
            x120 * x40 * x95,
            x122 * x48 * x95,
            x102 * x123 * x51,
            x102 * x122 * x73,
            x102 * x126 * x48,
            x120 * x60 * x87,
            x123 * x64 * x94,
            x122 * x62 * x94,
            x101 * x123 * x66,
            x101 * x122 * x64,
            x101 * x126 * x62,
            x128 * x48 * x87,
            x128 * x73 * x94,
            x130 * x48 * x94,
            x101 * x128 * x51,
            x101 * x130 * x73,
            x101 * x134 * x48,
            x123 * x33 * x76,
            x123 * x24 * x77,
            x122 * x24 * x76,
            x116 * x119 * x80,
            x122 * x26 * x77,
            x126 * x26 * x76,
            x128 * x33 * x62,
            x128 * x24 * x64,
            x130 * x24 * x62,
            x128 * x26 * x66,
            x130 * x26 * x64,
            x118 * x134 * x60,
            x135 * x33 * x48,
            x135 * x24 * x73,
            x136 * x24 * x48,
            x135 * x26 * x51,
            x118 * x136 * x40,
            x138 * x4,
            x140 * x141 * x55,
            x144 * x146 * x55,
            x140 * x146 * x67,
            x150 * x151 * x55,
            x144 * x151 * x67,
            x140 * x151 * x59,
            x153 * x55 * x89,
            x145 * x155 * x55,
            x145 * x153 * x67,
            x158 * x2,
            x155 * x159 * x47,
            x100 * x153 * x59,
            x140 * x70 * x89,
            x144 * x145 * x70,
            x140 * x145 * x72,
            x150 * x159 * x68,
            x100 * x144 * x72,
            x100 * x140 * x75,
            x160 * x55 * x88,
            x12 * x161,
            x160 * x162 * x47,
            x114 * (x157 * x60 + x3 * (2 * x108 + 2 * x109 + 3 * x148 + x149 + 2 * x154)),
            x161 * x47,
            x160 * x59 * x8,
            x153 * x70 * x88,
            x155 * x162 * x68,
            x13 * x153 * x72,
            x158 * x68,
            x155 * x72 * x8,
            x153 * x75 * x8,
            x140 * x81 * x88,
            x13 * x144 * x81,
            x13 * x140 * x82,
            x150 * x8 * x81,
            x144 * x8 * x82,
            x140 * x8 * x85,
            x120 * x141 * x86,
            x123 * x146 * x92,
            x122 * x146 * x96,
            x123 * x151 * x99,
            x122 * x151 * x92,
            x126 * x151 * x96,
            x104 * x123 * x89,
            x106 * x123 * x145,
            x104 * x122 * x145,
            x110 * x119 * x159,
            x100 * x106 * x122,
            x100 * x104 * x126,
            x128 * x89 * x96,
            x128 * x145 * x92,
            x130 * x145 * x96,
            x100 * x128 * x99,
            x100 * x130 * x92,
            x134 * x163 * x2,
            x111 * x123 * x88,
            x112 * x119 * x162,
            x111 * x122 * x13,
            x115 * x119,
            x112 * x122 * x8,
            x111 * x126 * x8,
            x104 * x128 * x88,
            x106 * x128 * x13,
            x104 * x13 * x130,
            x110 * x128 * x8,
            x106 * x130 * x8,
            x104 * x134 * x8,
            x135 * x88 * x96,
            x13 * x135 * x92,
            x12 * x136 * x163,
            x135 * x8 * x99,
            x136 * x8 * x92,
            x138 * x86,
            x141 * x165 * x48,
            x146 * x165 * x73,
            x146 * x168 * x48,
            x151 * x165 * x51,
            x151 * x168 * x73,
            x151 * x172 * x48,
            x165 * x62 * x89,
            x145 * x165 * x64,
            x145 * x168 * x62,
            x100 * x165 * x66,
            x100 * x168 * x64,
            x172 * x173 * x60,
            x175 * x48 * x89,
            x145 * x175 * x73,
            x145 * x177 * x48,
            x100 * x175 * x51,
            x173 * x177 * x40,
            x180 * x2,
            x165 * x76 * x88,
            x13 * x165 * x77,
            x13 * x168 * x76,
            x165 * x8 * x80,
            x168 * x77 * x8,
            x172 * x76 * x8,
            x175 * x62 * x88,
            x13 * x175 * x64,
            x177 * x181 * x60,
            x175 * x66 * x8,
            x177 * x64 * x8,
            x180 * x60,
            x182 * x48 * x88,
            x181 * x182 * x40,
            x12 * x183,
            x182 * x51 * x8,
            x183 * x40,
            x117 * (x179 * x68 + x3 * (2 * x132 + 2 * x133 + 3 * x170 + x171 + 2 * x176)),
        ]
    )


def quadrupole3d_23(a, A, b, B, C):
    """Cartesian 3D (df) quadrupole moment integrals.
    The origin is at C.

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
    x9 = x4 * x8
    x10 = -x1 - C[0]
    x11 = x10 * x8
    x12 = x3 * (x11 + x9)
    x13 = x3 * x8
    x14 = x10 * x9
    x15 = x13 + x14
    x16 = x15 * x4
    x17 = x12 + x16
    x18 = x17 * x4
    x19 = 2 * x13
    x20 = x10 * x19
    x21 = x10 ** 2 * x8
    x22 = x13 + x21
    x23 = x22 * x4
    x24 = x20 + x23
    x25 = x24 * x4
    x26 = x4 ** 2 * x8
    x27 = 3 * x13
    x28 = 2 * x14 + x27
    x29 = x3 * (x26 + x28)
    x30 = x3 * (x21 + x28)
    x31 = 2 * x29 + 3 * x30
    x32 = x25 + x30
    x33 = x32 * x4
    x34 = 4 * x10 * x13 + 2 * x12
    x35 = x3 * (2 * x16 + 2 * x23 + x34)
    x36 = x33 + x35
    x37 = x2 * x36 + x3 * (2 * x18 + 3 * x25 + x31)
    x38 = x2 * x32
    x39 = 3 * x12
    x40 = x19 * x4
    x41 = x13 + x26
    x42 = x4 * x41
    x43 = x40 + x42
    x44 = x3 * (3 * x16 + x39 + x43)
    x45 = x18 + x29
    x46 = x2 * x45
    x47 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x48 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x49 = numpy.pi * x0 * x48
    x50 = x47 * x49
    x51 = -x0 * (a * A[1] + b * B[1])
    x52 = -x51 - B[1]
    x53 = x35 + x38
    x54 = x17 * x2
    x55 = x2 * x24
    x56 = x50 * (x2 * x53 + x3 * (x25 + x31 + 2 * x54 + 2 * x55))
    x57 = -x0 * (a * A[2] + b * B[2])
    x58 = -x57 - B[2]
    x59 = x47 * x7
    x60 = x3 * x59
    x61 = x52 ** 2 * x59
    x62 = x60 + x61
    x63 = x30 + x55
    x64 = x2 * x22
    x65 = x15 * x2
    x66 = 2 * x65
    x67 = x2 * x63 + x3 * (x23 + x34 + x64 + x66)
    x68 = x48 * x7
    x69 = x50 * x58
    x70 = x3 * x68
    x71 = x58 ** 2 * x68
    x72 = x70 + x71
    x73 = x52 * x60
    x74 = 2 * x73
    x75 = x52 * x62
    x76 = x74 + x75
    x77 = x11 * x2
    x78 = x20 + x64
    x79 = x2 * x78 + x3 * (x21 + x27 + 2 * x77)
    x80 = x58 * x68
    x81 = x52 * x59
    x82 = x58 * x70
    x83 = 2 * x82
    x84 = x58 * x72
    x85 = x83 + x84
    x86 = -x51 - A[1]
    x87 = x37 * x50
    x88 = x59 * x86
    x89 = x52 * x88
    x90 = x60 + x89
    x91 = x62 * x86
    x92 = x74 + x91
    x93 = 3 * x60
    x94 = x3 * (3 * x61 + x93) + x76 * x86
    x95 = -x57 - A[2]
    x96 = x50 * x95
    x97 = x68 * x95
    x98 = x58 * x97
    x99 = x70 + x98
    x100 = x72 * x95
    x101 = x100 + x83
    x102 = 3 * x70
    x103 = x3 * (x102 + 3 * x71) + x85 * x95
    x104 = x59 * x86 ** 2 + x60
    x105 = x3 * (x81 + x88) + x86 * x90
    x106 = x61 + x93
    x107 = x3 * (x106 + 2 * x89) + x86 * x92
    x108 = x3 * (8 * x73 + x75 + 3 * x91) + x86 * x94
    x109 = x68 * x95 ** 2 + x70
    x110 = x3 * (x80 + x97) + x95 * x99
    x111 = x102 + x71
    x112 = x101 * x95 + x3 * (x111 + 2 * x98)
    x113 = x103 * x95 + x3 * (3 * x100 + 8 * x82 + x84)
    x114 = -x51 - C[1]
    x115 = x44 + x46
    x116 = x2 * x43 + x3 * (3 * x26 + x27)
    x117 = x50 * (x115 * x2 + x3 * (x116 + x18 + 4 * x29 + 3 * x54))
    x118 = x114 * x81
    x119 = x118 + x60
    x120 = x29 + x54
    x121 = x2 * x41
    x122 = x121 + x40
    x123 = x120 * x2 + x3 * (x122 + x16 + x39 + x66)
    x124 = x114 * x59
    x125 = x3 * (x124 + x81)
    x126 = x119 * x52
    x127 = x125 + x126
    x128 = x2 * x9
    x129 = x12 + x65
    x130 = x129 * x2 + x3 * (x128 + x14 + x27 + x77)
    x131 = 2 * x118
    x132 = x3 * (x106 + x131)
    x133 = x127 * x52
    x134 = x132 + x133
    x135 = x2 * x8
    x136 = x13 + x77
    x137 = x136 * x2 + x3 * (x11 + x135)
    x138 = x114 * x88
    x139 = x138 + x60
    x140 = x119 * x86
    x141 = x125 + x140
    x142 = x127 * x86
    x143 = x132 + x142
    x144 = 3 * x125
    x145 = x3 * (3 * x126 + x144 + x76)
    x146 = x134 * x86
    x147 = x145 + x146
    x148 = x139 * x86 + x3 * (x124 + x88)
    x149 = x141 * x86 + x3 * (x118 + x138 + x89 + x93)
    x150 = 2 * x140
    x151 = x143 * x86 + x3 * (x126 + x144 + x150 + x92)
    x152 = x49 * x6
    x153 = x152 * (x147 * x86 + x3 * (4 * x132 + x133 + 3 * x142 + x94))
    x154 = x10 * x152
    x155 = numpy.pi * x0 * x47 * x6
    x156 = x10 * x155
    x157 = -x57 - C[2]
    x158 = x157 * x50
    x159 = x157 * x80
    x160 = x159 + x70
    x161 = x157 * x68
    x162 = x3 * (x161 + x80)
    x163 = x160 * x58
    x164 = x162 + x163
    x165 = 2 * x159
    x166 = x3 * (x111 + x165)
    x167 = x164 * x58
    x168 = x166 + x167
    x169 = x157 * x97
    x170 = x169 + x70
    x171 = x160 * x95
    x172 = x162 + x171
    x173 = x164 * x95
    x174 = x166 + x173
    x175 = 3 * x162
    x176 = x3 * (3 * x163 + x175 + x85)
    x177 = x168 * x95
    x178 = x176 + x177
    x179 = x170 * x95 + x3 * (x161 + x97)
    x180 = x172 * x95 + x3 * (x102 + x159 + x169 + x98)
    x181 = 2 * x171
    x182 = x174 * x95 + x3 * (x101 + x163 + x175 + x181)
    x183 = x155 * (x178 * x95 + x3 * (x103 + 4 * x166 + x167 + 3 * x173))
    x184 = x114 ** 2 * x59
    x185 = x184 + x60
    x186 = x116 * x2 + x3 * (3 * x121 + 8 * x13 * x4 + x42)
    x187 = x114 * x60
    x188 = 2 * x187
    x189 = x185 * x52
    x190 = x188 + x189
    x191 = x122 * x2 + x3 * (2 * x128 + x26 + x27)
    x192 = x184 + x93
    x193 = x3 * (x131 + x192)
    x194 = x190 * x52
    x195 = x193 + x194
    x196 = x128 + x13
    x197 = x196 * x2 + x3 * (x135 + x9)
    x198 = x195 * x52
    x199 = 2 * x125 + 4 * x187
    x200 = x3 * (2 * x126 + 2 * x189 + x199)
    x201 = x198 + x200
    x202 = x13 + x2 ** 2 * x8
    x203 = x185 * x86
    x204 = x188 + x203
    x205 = x190 * x86
    x206 = x193 + x205
    x207 = x195 * x86
    x208 = x200 + x207
    x209 = 2 * x132 + 3 * x193
    x210 = x201 * x86 + x3 * (2 * x133 + 3 * x194 + x209)
    x211 = x152 * x210
    x212 = x152 * x2
    x213 = x204 * x86 + x3 * (2 * x138 + x192)
    x214 = x206 * x86 + x3 * (x150 + x189 + x199 + x203)
    x215 = x152 * (x208 * x86 + x3 * (2 * x142 + x194 + 2 * x205 + x209))
    x216 = x152 * x4
    x217 = x114 * x155
    x218 = x157 ** 2 * x68
    x219 = x218 + x70
    x220 = x157 * x70
    x221 = 2 * x220
    x222 = x219 * x58
    x223 = x221 + x222
    x224 = x102 + x218
    x225 = x3 * (x165 + x224)
    x226 = x223 * x58
    x227 = x225 + x226
    x228 = x227 * x58
    x229 = 2 * x162 + 4 * x220
    x230 = x3 * (2 * x163 + 2 * x222 + x229)
    x231 = x228 + x230
    x232 = x155 * x2
    x233 = x219 * x95
    x234 = x221 + x233
    x235 = x223 * x95
    x236 = x225 + x235
    x237 = x227 * x95
    x238 = x230 + x237
    x239 = 2 * x166 + 3 * x225
    x240 = x231 * x95 + x3 * (2 * x167 + 3 * x226 + x239)
    x241 = x155 * x240
    x242 = x155 * x4
    x243 = x234 * x95 + x3 * (2 * x169 + x224)
    x244 = x236 * x95 + x3 * (x181 + x222 + x229 + x233)
    x245 = x155 * (x238 * x95 + x3 * (2 * x173 + x226 + 2 * x235 + x239))

    # 360 item(s)
    return numpy.array(
        [
            x50 * (x2 * x37 + x3 * (x33 + 4 * x35 + 3 * x38 + 2 * x44 + 2 * x46)),
            x52 * x56,
            x56 * x58,
            x62 * x67 * x68,
            x52 * x67 * x69,
            x59 * x67 * x72,
            x68 * x76 * x79,
            x62 * x79 * x80,
            x72 * x79 * x81,
            x59 * x79 * x85,
            x86 * x87,
            x53 * x68 * x90,
            x53 * x69 * x86,
            x63 * x68 * x92,
            x63 * x80 * x90,
            x63 * x72 * x88,
            x68 * x78 * x94,
            x78 * x80 * x92,
            x72 * x78 * x90,
            x78 * x85 * x88,
            x87 * x95,
            x52 * x53 * x96,
            x53 * x59 * x99,
            x62 * x63 * x97,
            x63 * x81 * x99,
            x101 * x59 * x63,
            x76 * x78 * x97,
            x62 * x78 * x99,
            x101 * x78 * x81,
            x103 * x59 * x78,
            x104 * x36 * x68,
            x105 * x32 * x68,
            x104 * x32 * x80,
            x107 * x24 * x68,
            x105 * x24 * x80,
            x104 * x24 * x72,
            x108 * x22 * x68,
            x107 * x22 * x80,
            x105 * x22 * x72,
            x104 * x22 * x85,
            x36 * x86 * x96,
            x32 * x90 * x97,
            x32 * x88 * x99,
            x24 * x92 * x97,
            x24 * x90 * x99,
            x101 * x24 * x88,
            x22 * x94 * x97,
            x22 * x92 * x99,
            x101 * x22 * x90,
            x103 * x22 * x88,
            x109 * x36 * x59,
            x109 * x32 * x81,
            x110 * x32 * x59,
            x109 * x24 * x62,
            x110 * x24 * x81,
            x112 * x24 * x59,
            x109 * x22 * x76,
            x110 * x22 * x62,
            x112 * x22 * x81,
            x113 * x22 * x59,
            x114 * x117,
            x119 * x123 * x68,
            x114 * x123 * x69,
            x127 * x130 * x68,
            x119 * x130 * x80,
            x124 * x130 * x72,
            x134 * x137 * x68,
            x127 * x137 * x80,
            x119 * x137 * x72,
            x124 * x137 * x85,
            x115 * x139 * x68,
            x120 * x141 * x68,
            x120 * x139 * x80,
            x129 * x143 * x68,
            x129 * x141 * x80,
            x129 * x139 * x72,
            x136 * x147 * x68,
            x136 * x143 * x80,
            x136 * x141 * x72,
            x136 * x139 * x85,
            x114 * x115 * x96,
            x119 * x120 * x97,
            x120 * x124 * x99,
            x127 * x129 * x97,
            x119 * x129 * x99,
            x101 * x124 * x129,
            x134 * x136 * x97,
            x127 * x136 * x99,
            x101 * x119 * x136,
            x103 * x124 * x136,
            x148 * x45 * x68,
            x149 * x17 * x68,
            x148 * x17 * x80,
            x15 * x151 * x68,
            x149 * x15 * x80,
            x148 * x15 * x72,
            x10 * x153,
            x151 * x154 * x58,
            x11 * x149 * x72,
            x11 * x148 * x85,
            x139 * x45 * x97,
            x141 * x17 * x97,
            x139 * x17 * x99,
            x143 * x15 * x97,
            x141 * x15 * x99,
            x101 * x139 * x15,
            x147 * x154 * x95,
            x11 * x143 * x99,
            x101 * x11 * x141,
            x103 * x11 * x139,
            x109 * x124 * x45,
            x109 * x119 * x17,
            x110 * x124 * x17,
            x109 * x127 * x15,
            x110 * x119 * x15,
            x112 * x124 * x15,
            x109 * x11 * x134,
            x11 * x110 * x127,
            x11 * x112 * x119,
            x113 * x114 * x156,
            x117 * x157,
            x123 * x158 * x52,
            x123 * x160 * x59,
            x130 * x161 * x62,
            x130 * x160 * x81,
            x130 * x164 * x59,
            x137 * x161 * x76,
            x137 * x160 * x62,
            x137 * x164 * x81,
            x137 * x168 * x59,
            x115 * x158 * x86,
            x120 * x161 * x90,
            x120 * x160 * x88,
            x129 * x161 * x92,
            x129 * x160 * x90,
            x129 * x164 * x88,
            x136 * x161 * x94,
            x136 * x160 * x92,
            x136 * x164 * x90,
            x136 * x168 * x88,
            x115 * x170 * x59,
            x120 * x170 * x81,
            x120 * x172 * x59,
            x129 * x170 * x62,
            x129 * x172 * x81,
            x129 * x174 * x59,
            x136 * x170 * x76,
            x136 * x172 * x62,
            x136 * x174 * x81,
            x136 * x178 * x59,
            x104 * x161 * x45,
            x105 * x161 * x17,
            x104 * x160 * x17,
            x107 * x15 * x161,
            x105 * x15 * x160,
            x104 * x15 * x164,
            x108 * x154 * x157,
            x107 * x11 * x160,
            x105 * x11 * x164,
            x104 * x11 * x168,
            x170 * x45 * x88,
            x17 * x170 * x90,
            x17 * x172 * x88,
            x15 * x170 * x92,
            x15 * x172 * x90,
            x15 * x174 * x88,
            x11 * x170 * x94,
            x11 * x172 * x92,
            x11 * x174 * x90,
            x156 * x178 * x86,
            x179 * x45 * x59,
            x17 * x179 * x81,
            x17 * x180 * x59,
            x15 * x179 * x62,
            x15 * x180 * x81,
            x15 * x182 * x59,
            x11 * x179 * x76,
            x11 * x180 * x62,
            x156 * x182 * x52,
            x10 * x183,
            x185 * x186 * x68,
            x190 * x191 * x68,
            x185 * x191 * x80,
            x195 * x197 * x68,
            x190 * x197 * x80,
            x185 * x197 * x72,
            x201 * x202 * x68,
            x195 * x202 * x80,
            x190 * x202 * x72,
            x185 * x202 * x85,
            x116 * x204 * x68,
            x122 * x206 * x68,
            x122 * x204 * x80,
            x196 * x208 * x68,
            x196 * x206 * x80,
            x196 * x204 * x72,
            x2 * x211,
            x208 * x212 * x58,
            x135 * x206 * x72,
            x135 * x204 * x85,
            x116 * x185 * x97,
            x122 * x190 * x97,
            x122 * x185 * x99,
            x195 * x196 * x97,
            x190 * x196 * x99,
            x101 * x185 * x196,
            x201 * x212 * x95,
            x135 * x195 * x99,
            x101 * x135 * x190,
            x103 * x135 * x185,
            x213 * x43 * x68,
            x214 * x41 * x68,
            x213 * x41 * x80,
            x215 * x4,
            x214 * x216 * x58,
            x213 * x72 * x9,
            x152 * (x210 * x86 + x3 * (2 * x145 + 2 * x146 + x198 + 4 * x200 + 3 * x207)),
            x215 * x58,
            x214 * x72 * x8,
            x213 * x8 * x85,
            x204 * x43 * x97,
            x206 * x41 * x97,
            x204 * x41 * x99,
            x208 * x216 * x95,
            x206 * x9 * x99,
            x101 * x204 * x9,
            x211 * x95,
            x208 * x8 * x99,
            x101 * x206 * x8,
            x103 * x204 * x8,
            x109 * x185 * x43,
            x109 * x190 * x41,
            x110 * x185 * x41,
            x109 * x195 * x9,
            x110 * x190 * x9,
            x112 * x185 * x9,
            x109 * x201 * x8,
            x110 * x195 * x8,
            x112 * x190 * x8,
            x113 * x185 * x8,
            x114 * x158 * x186,
            x119 * x161 * x191,
            x124 * x160 * x191,
            x127 * x161 * x197,
            x119 * x160 * x197,
            x124 * x164 * x197,
            x134 * x161 * x202,
            x127 * x160 * x202,
            x119 * x164 * x202,
            x124 * x168 * x202,
            x116 * x139 * x161,
            x122 * x141 * x161,
            x122 * x139 * x160,
            x143 * x161 * x196,
            x141 * x160 * x196,
            x139 * x164 * x196,
            x147 * x157 * x212,
            x135 * x143 * x160,
            x135 * x141 * x164,
            x135 * x139 * x168,
            x116 * x124 * x170,
            x119 * x122 * x170,
            x122 * x124 * x172,
            x127 * x170 * x196,
            x119 * x172 * x196,
            x124 * x174 * x196,
            x134 * x135 * x170,
            x127 * x135 * x172,
            x119 * x135 * x174,
            x178 * x2 * x217,
            x148 * x161 * x43,
            x149 * x161 * x41,
            x148 * x160 * x41,
            x151 * x157 * x216,
            x149 * x160 * x9,
            x148 * x164 * x9,
            x153 * x157,
            x151 * x160 * x8,
            x149 * x164 * x8,
            x148 * x168 * x8,
            x139 * x170 * x43,
            x141 * x170 * x41,
            x139 * x172 * x41,
            x143 * x170 * x9,
            x141 * x172 * x9,
            x139 * x174 * x9,
            x147 * x170 * x8,
            x143 * x172 * x8,
            x141 * x174 * x8,
            x139 * x178 * x8,
            x124 * x179 * x43,
            x119 * x179 * x41,
            x124 * x180 * x41,
            x127 * x179 * x9,
            x119 * x180 * x9,
            x182 * x217 * x4,
            x134 * x179 * x8,
            x127 * x180 * x8,
            x119 * x182 * x8,
            x114 * x183,
            x186 * x219 * x59,
            x191 * x219 * x81,
            x191 * x223 * x59,
            x197 * x219 * x62,
            x197 * x223 * x81,
            x197 * x227 * x59,
            x202 * x219 * x76,
            x202 * x223 * x62,
            x202 * x227 * x81,
            x202 * x231 * x59,
            x116 * x219 * x88,
            x122 * x219 * x90,
            x122 * x223 * x88,
            x196 * x219 * x92,
            x196 * x223 * x90,
            x196 * x227 * x88,
            x135 * x219 * x94,
            x135 * x223 * x92,
            x135 * x227 * x90,
            x231 * x232 * x86,
            x116 * x234 * x59,
            x122 * x234 * x81,
            x122 * x236 * x59,
            x196 * x234 * x62,
            x196 * x236 * x81,
            x196 * x238 * x59,
            x135 * x234 * x76,
            x135 * x236 * x62,
            x232 * x238 * x52,
            x2 * x241,
            x104 * x219 * x43,
            x105 * x219 * x41,
            x104 * x223 * x41,
            x107 * x219 * x9,
            x105 * x223 * x9,
            x104 * x227 * x9,
            x108 * x219 * x8,
            x107 * x223 * x8,
            x105 * x227 * x8,
            x104 * x231 * x8,
            x234 * x43 * x88,
            x234 * x41 * x90,
            x236 * x41 * x88,
            x234 * x9 * x92,
            x236 * x9 * x90,
            x238 * x242 * x86,
            x234 * x8 * x94,
            x236 * x8 * x92,
            x238 * x8 * x90,
            x241 * x86,
            x243 * x43 * x59,
            x243 * x41 * x81,
            x244 * x41 * x59,
            x243 * x62 * x9,
            x242 * x244 * x52,
            x245 * x4,
            x243 * x76 * x8,
            x244 * x62 * x8,
            x245 * x52,
            x155 * (x240 * x95 + x3 * (2 * x176 + 2 * x177 + x228 + 4 * x230 + 3 * x237)),
        ]
    )


def quadrupole3d_24(a, A, b, B, C):
    """Cartesian 3D (dg) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - B[0]
    x5 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x6 = a * b * x0
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = x4 * x7
    x9 = x5 * x8
    x10 = -x1 - C[0]
    x11 = x5 * x7
    x12 = x10 * x11
    x13 = x3 * (x12 + x9)
    x14 = x11 * x3
    x15 = x10 * x9
    x16 = x14 + x15
    x17 = x16 * x4
    x18 = x13 + x17
    x19 = x18 * x4
    x20 = x12 * x3
    x21 = 2 * x20
    x22 = x10 ** 2 * x11
    x23 = x14 + x22
    x24 = x23 * x4
    x25 = x21 + x24
    x26 = x25 * x4
    x27 = x11 * x4 ** 2
    x28 = 3 * x14
    x29 = 2 * x15 + x28
    x30 = x3 * (x27 + x29)
    x31 = x3 * (x22 + x29)
    x32 = 2 * x30 + 3 * x31
    x33 = x3 * (2 * x19 + 3 * x26 + x32)
    x34 = x26 + x31
    x35 = x34 * x4
    x36 = 2 * x13 + 4 * x20
    x37 = x3 * (2 * x17 + 2 * x24 + x36)
    x38 = x35 + x37
    x39 = x38 * x4
    x40 = x33 + x39
    x41 = x19 + x30
    x42 = x4 * x41
    x43 = 3 * x13
    x44 = x3 * x9
    x45 = 2 * x44
    x46 = x14 + x27
    x47 = x4 * x46
    x48 = x45 + x47
    x49 = x3 * (3 * x17 + x43 + x48)
    x50 = 4 * x37 + 2 * x49
    x51 = x2 * x40 + x3 * (4 * x35 + 2 * x42 + x50)
    x52 = 4 * x30
    x53 = x3 * (3 * x27 + x28)
    x54 = x4 * x48
    x55 = x53 + x54
    x56 = x3 * (4 * x19 + x52 + x55)
    x57 = x2 * x38
    x58 = x42 + x49
    x59 = x2 * x58
    x60 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x61 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x62 = numpy.pi * x0 * x61
    x63 = x60 * x62
    x64 = -x0 * (a * A[1] + b * B[1])
    x65 = -x64 - B[1]
    x66 = x33 + x57
    x67 = x2 * x41
    x68 = x2 * x34
    x69 = x63 * (x2 * x66 + x3 * (x35 + x50 + 2 * x67 + 3 * x68))
    x70 = -x0 * (a * A[2] + b * B[2])
    x71 = -x70 - B[2]
    x72 = x5 * x60
    x73 = x3 * x72
    x74 = x65 ** 2 * x72
    x75 = x73 + x74
    x76 = x37 + x68
    x77 = x18 * x2
    x78 = x2 * x25
    x79 = x2 * x76 + x3 * (x26 + x32 + 2 * x77 + 2 * x78)
    x80 = x5 * x61
    x81 = x63 * x71
    x82 = x3 * x80
    x83 = x71 ** 2 * x80
    x84 = x82 + x83
    x85 = x65 * x73
    x86 = 2 * x85
    x87 = x65 * x75
    x88 = x86 + x87
    x89 = x31 + x78
    x90 = x2 * x23
    x91 = x16 * x2
    x92 = 2 * x91
    x93 = x2 * x89 + x3 * (x24 + x36 + x90 + x92)
    x94 = x71 * x80
    x95 = x65 * x72
    x96 = x71 * x82
    x97 = 2 * x96
    x98 = x71 * x84
    x99 = x97 + x98
    x100 = 3 * x73
    x101 = x3 * (x100 + 3 * x74)
    x102 = x65 * x88
    x103 = x101 + x102
    x104 = x12 * x2
    x105 = x21 + x90
    x106 = x105 * x2 + x3 * (2 * x104 + x22 + x28)
    x107 = 3 * x82
    x108 = x3 * (x107 + 3 * x83)
    x109 = x71 * x99
    x110 = x108 + x109
    x111 = -x64 - A[1]
    x112 = x51 * x63
    x113 = x111 * x72
    x114 = x113 * x65
    x115 = x114 + x73
    x116 = x111 * x75
    x117 = x116 + x86
    x118 = x111 * x88
    x119 = x101 + x118
    x120 = 8 * x85
    x121 = x103 * x111 + x3 * (x120 + 4 * x87)
    x122 = -x70 - A[2]
    x123 = x122 * x63
    x124 = x122 * x80
    x125 = x124 * x71
    x126 = x125 + x82
    x127 = x122 * x84
    x128 = x127 + x97
    x129 = x122 * x99
    x130 = x108 + x129
    x131 = 8 * x96
    x132 = x110 * x122 + x3 * (x131 + 4 * x98)
    x133 = x111 ** 2 * x72 + x73
    x134 = x111 * x115 + x3 * (x113 + x95)
    x135 = x100 + x74
    x136 = x111 * x117 + x3 * (2 * x114 + x135)
    x137 = x111 * x119 + x3 * (3 * x116 + x120 + x87)
    x138 = x111 * x121 + x3 * (5 * x101 + x102 + 4 * x118)
    x139 = x122 ** 2 * x80 + x82
    x140 = x122 * x126 + x3 * (x124 + x94)
    x141 = x107 + x83
    x142 = x122 * x128 + x3 * (2 * x125 + x141)
    x143 = x122 * x130 + x3 * (3 * x127 + x131 + x98)
    x144 = x122 * x132 + x3 * (5 * x108 + x109 + 4 * x129)
    x145 = -x64 - C[1]
    x146 = x56 + x59
    x147 = 8 * x44
    x148 = x2 * x55 + x3 * (x147 + 4 * x47)
    x149 = x63 * (x146 * x2 + x3 * (x148 + x42 + 5 * x49 + 4 * x67))
    x150 = x145 * x95
    x151 = x150 + x73
    x152 = x49 + x67
    x153 = x2 * x48
    x154 = x153 + x53
    x155 = x152 * x2 + x3 * (x154 + x19 + x52 + 3 * x77)
    x156 = x145 * x72
    x157 = x3 * (x156 + x95)
    x158 = x151 * x65
    x159 = x157 + x158
    x160 = x30 + x77
    x161 = x2 * x46
    x162 = x161 + x45
    x163 = x160 * x2 + x3 * (x162 + x17 + x43 + x92)
    x164 = 2 * x150
    x165 = x3 * (x135 + x164)
    x166 = x159 * x65
    x167 = x165 + x166
    x168 = x2 * x9
    x169 = x13 + x91
    x170 = x169 * x2 + x3 * (x104 + x15 + x168 + x28)
    x171 = 3 * x157
    x172 = x3 * (3 * x158 + x171 + x88)
    x173 = x167 * x65
    x174 = x172 + x173
    x175 = x11 * x2
    x176 = x104 + x14
    x177 = x176 * x2 + x3 * (x12 + x175)
    x178 = x113 * x145
    x179 = x178 + x73
    x180 = x111 * x151
    x181 = x157 + x180
    x182 = x111 * x159
    x183 = x165 + x182
    x184 = x111 * x167
    x185 = x172 + x184
    x186 = 4 * x165
    x187 = x3 * (x103 + 4 * x166 + x186)
    x188 = x111 * x174
    x189 = x187 + x188
    x190 = x111 * x179 + x3 * (x113 + x156)
    x191 = x111 * x181 + x3 * (x100 + x114 + x150 + x178)
    x192 = 2 * x180
    x193 = x111 * x183 + x3 * (x117 + x158 + x171 + x192)
    x194 = x111 * x185 + x3 * (x119 + x166 + 3 * x182 + x186)
    x195 = x62 * x7
    x196 = x195 * (x111 * x189 + x3 * (x121 + 5 * x172 + x173 + 4 * x184))
    x197 = x195 * x71
    x198 = x10 * x195
    x199 = numpy.pi * x0 * x60
    x200 = x199 * x7
    x201 = x10 * x200
    x202 = -x70 - C[2]
    x203 = x202 * x63
    x204 = x202 * x94
    x205 = x204 + x82
    x206 = x202 * x80
    x207 = x3 * (x206 + x94)
    x208 = x205 * x71
    x209 = x207 + x208
    x210 = 2 * x204
    x211 = x3 * (x141 + x210)
    x212 = x209 * x71
    x213 = x211 + x212
    x214 = 3 * x207
    x215 = x3 * (3 * x208 + x214 + x99)
    x216 = x213 * x71
    x217 = x215 + x216
    x218 = x124 * x202
    x219 = x218 + x82
    x220 = x122 * x205
    x221 = x207 + x220
    x222 = x122 * x209
    x223 = x211 + x222
    x224 = x122 * x213
    x225 = x215 + x224
    x226 = 4 * x211
    x227 = x3 * (x110 + 4 * x212 + x226)
    x228 = x122 * x217
    x229 = x227 + x228
    x230 = x122 * x219 + x3 * (x124 + x206)
    x231 = x122 * x221 + x3 * (x107 + x125 + x204 + x218)
    x232 = 2 * x220
    x233 = x122 * x223 + x3 * (x128 + x208 + x214 + x232)
    x234 = x122 * x225 + x3 * (x130 + x212 + 3 * x222 + x226)
    x235 = x200 * x65
    x236 = x200 * (x122 * x229 + x3 * (x132 + 5 * x215 + x216 + 4 * x224))
    x237 = x145 ** 2 * x72
    x238 = x237 + x73
    x239 = x148 * x2 + x3 * (4 * x153 + 5 * x53 + x54)
    x240 = x145 * x73
    x241 = 2 * x240
    x242 = x238 * x65
    x243 = x241 + x242
    x244 = x154 * x2 + x3 * (x147 + 3 * x161 + x47)
    x245 = x100 + x237
    x246 = x3 * (x164 + x245)
    x247 = x243 * x65
    x248 = x246 + x247
    x249 = x162 * x2 + x3 * (2 * x168 + x27 + x28)
    x250 = x248 * x65
    x251 = 2 * x157 + 4 * x240
    x252 = x3 * (2 * x158 + 2 * x242 + x251)
    x253 = x250 + x252
    x254 = x14 + x168
    x255 = x2 * x254 + x3 * (x175 + x9)
    x256 = 2 * x165 + 3 * x246
    x257 = x3 * (2 * x166 + 3 * x247 + x256)
    x258 = x253 * x65
    x259 = x257 + x258
    x260 = x11 * x2 ** 2 + x14
    x261 = x111 * x238
    x262 = x241 + x261
    x263 = x111 * x243
    x264 = x246 + x263
    x265 = x111 * x248
    x266 = x252 + x265
    x267 = x111 * x253
    x268 = x257 + x267
    x269 = 2 * x172 + 4 * x252
    x270 = x111 * x259 + x3 * (2 * x173 + 4 * x250 + x269)
    x271 = x195 * x270
    x272 = x195 * x2
    x273 = x111 * x262 + x3 * (2 * x178 + x245)
    x274 = x111 * x264 + x3 * (x192 + x242 + x251 + x261)
    x275 = x111 * x266 + x3 * (2 * x182 + x247 + x256 + 2 * x263)
    x276 = x111 * x268 + x3 * (2 * x184 + x250 + 3 * x265 + x269)
    x277 = x62 * x8
    x278 = x2 * x200
    x279 = x199 * x8
    x280 = x202 ** 2 * x80
    x281 = x280 + x82
    x282 = x202 * x82
    x283 = 2 * x282
    x284 = x281 * x71
    x285 = x283 + x284
    x286 = x107 + x280
    x287 = x3 * (x210 + x286)
    x288 = x285 * x71
    x289 = x287 + x288
    x290 = x289 * x71
    x291 = 2 * x207 + 4 * x282
    x292 = x3 * (2 * x208 + 2 * x284 + x291)
    x293 = x290 + x292
    x294 = 2 * x211 + 3 * x287
    x295 = x3 * (2 * x212 + 3 * x288 + x294)
    x296 = x293 * x71
    x297 = x295 + x296
    x298 = x122 * x281
    x299 = x283 + x298
    x300 = x122 * x285
    x301 = x287 + x300
    x302 = x122 * x289
    x303 = x292 + x302
    x304 = x122 * x293
    x305 = x295 + x304
    x306 = 2 * x215 + 4 * x292
    x307 = x122 * x297 + x3 * (2 * x216 + 4 * x290 + x306)
    x308 = x200 * x307
    x309 = x122 * x299 + x3 * (2 * x218 + x286)
    x310 = x122 * x301 + x3 * (x232 + x284 + x291 + x298)
    x311 = x122 * x303 + x3 * (2 * x222 + x288 + x294 + 2 * x300)
    x312 = x122 * x305 + x3 * (2 * x224 + x290 + 3 * x302 + x306)

    # 540 item(s)
    return numpy.array(
        [
            x63 * (x2 * x51 + x3 * (5 * x33 + x39 + 2 * x56 + 4 * x57 + 2 * x59)),
            x65 * x69,
            x69 * x71,
            x75 * x79 * x80,
            x65 * x79 * x81,
            x72 * x79 * x84,
            x80 * x88 * x93,
            x75 * x93 * x94,
            x84 * x93 * x95,
            x72 * x93 * x99,
            x103 * x106 * x80,
            x106 * x88 * x94,
            x106 * x75 * x84,
            x106 * x95 * x99,
            x106 * x110 * x72,
            x111 * x112,
            x115 * x66 * x80,
            x111 * x66 * x81,
            x117 * x76 * x80,
            x115 * x76 * x94,
            x113 * x76 * x84,
            x119 * x80 * x89,
            x117 * x89 * x94,
            x115 * x84 * x89,
            x113 * x89 * x99,
            x105 * x121 * x80,
            x105 * x119 * x94,
            x105 * x117 * x84,
            x105 * x115 * x99,
            x105 * x110 * x113,
            x112 * x122,
            x123 * x65 * x66,
            x126 * x66 * x72,
            x124 * x75 * x76,
            x126 * x76 * x95,
            x128 * x72 * x76,
            x124 * x88 * x89,
            x126 * x75 * x89,
            x128 * x89 * x95,
            x130 * x72 * x89,
            x103 * x105 * x124,
            x105 * x126 * x88,
            x105 * x128 * x75,
            x105 * x130 * x95,
            x105 * x132 * x72,
            x133 * x40 * x80,
            x134 * x38 * x80,
            x133 * x38 * x94,
            x136 * x34 * x80,
            x134 * x34 * x94,
            x133 * x34 * x84,
            x137 * x25 * x80,
            x136 * x25 * x94,
            x134 * x25 * x84,
            x133 * x25 * x99,
            x138 * x23 * x80,
            x137 * x23 * x94,
            x136 * x23 * x84,
            x134 * x23 * x99,
            x110 * x133 * x23,
            x111 * x123 * x40,
            x115 * x124 * x38,
            x113 * x126 * x38,
            x117 * x124 * x34,
            x115 * x126 * x34,
            x113 * x128 * x34,
            x119 * x124 * x25,
            x117 * x126 * x25,
            x115 * x128 * x25,
            x113 * x130 * x25,
            x121 * x124 * x23,
            x119 * x126 * x23,
            x117 * x128 * x23,
            x115 * x130 * x23,
            x113 * x132 * x23,
            x139 * x40 * x72,
            x139 * x38 * x95,
            x140 * x38 * x72,
            x139 * x34 * x75,
            x140 * x34 * x95,
            x142 * x34 * x72,
            x139 * x25 * x88,
            x140 * x25 * x75,
            x142 * x25 * x95,
            x143 * x25 * x72,
            x103 * x139 * x23,
            x140 * x23 * x88,
            x142 * x23 * x75,
            x143 * x23 * x95,
            x144 * x23 * x72,
            x145 * x149,
            x151 * x155 * x80,
            x145 * x155 * x81,
            x159 * x163 * x80,
            x151 * x163 * x94,
            x156 * x163 * x84,
            x167 * x170 * x80,
            x159 * x170 * x94,
            x151 * x170 * x84,
            x156 * x170 * x99,
            x174 * x177 * x80,
            x167 * x177 * x94,
            x159 * x177 * x84,
            x151 * x177 * x99,
            x110 * x156 * x177,
            x146 * x179 * x80,
            x152 * x181 * x80,
            x152 * x179 * x94,
            x160 * x183 * x80,
            x160 * x181 * x94,
            x160 * x179 * x84,
            x169 * x185 * x80,
            x169 * x183 * x94,
            x169 * x181 * x84,
            x169 * x179 * x99,
            x176 * x189 * x80,
            x176 * x185 * x94,
            x176 * x183 * x84,
            x176 * x181 * x99,
            x110 * x176 * x179,
            x123 * x145 * x146,
            x124 * x151 * x152,
            x126 * x152 * x156,
            x124 * x159 * x160,
            x126 * x151 * x160,
            x128 * x156 * x160,
            x124 * x167 * x169,
            x126 * x159 * x169,
            x128 * x151 * x169,
            x130 * x156 * x169,
            x124 * x174 * x176,
            x126 * x167 * x176,
            x128 * x159 * x176,
            x130 * x151 * x176,
            x132 * x156 * x176,
            x190 * x58 * x80,
            x191 * x41 * x80,
            x190 * x41 * x94,
            x18 * x193 * x80,
            x18 * x191 * x94,
            x18 * x190 * x84,
            x16 * x194 * x80,
            x16 * x193 * x94,
            x16 * x191 * x84,
            x16 * x190 * x99,
            x10 * x196,
            x10 * x194 * x197,
            x12 * x193 * x84,
            x12 * x191 * x99,
            x110 * x12 * x190,
            x124 * x179 * x58,
            x124 * x181 * x41,
            x126 * x179 * x41,
            x124 * x18 * x183,
            x126 * x18 * x181,
            x128 * x179 * x18,
            x124 * x16 * x185,
            x126 * x16 * x183,
            x128 * x16 * x181,
            x130 * x16 * x179,
            x122 * x189 * x198,
            x12 * x126 * x185,
            x12 * x128 * x183,
            x12 * x130 * x181,
            x12 * x132 * x179,
            x139 * x156 * x58,
            x139 * x151 * x41,
            x140 * x156 * x41,
            x139 * x159 * x18,
            x140 * x151 * x18,
            x142 * x156 * x18,
            x139 * x16 * x167,
            x140 * x159 * x16,
            x142 * x151 * x16,
            x143 * x156 * x16,
            x12 * x139 * x174,
            x12 * x140 * x167,
            x12 * x142 * x159,
            x12 * x143 * x151,
            x144 * x145 * x201,
            x149 * x202,
            x155 * x203 * x65,
            x155 * x205 * x72,
            x163 * x206 * x75,
            x163 * x205 * x95,
            x163 * x209 * x72,
            x170 * x206 * x88,
            x170 * x205 * x75,
            x170 * x209 * x95,
            x170 * x213 * x72,
            x103 * x177 * x206,
            x177 * x205 * x88,
            x177 * x209 * x75,
            x177 * x213 * x95,
            x177 * x217 * x72,
            x111 * x146 * x203,
            x115 * x152 * x206,
            x113 * x152 * x205,
            x117 * x160 * x206,
            x115 * x160 * x205,
            x113 * x160 * x209,
            x119 * x169 * x206,
            x117 * x169 * x205,
            x115 * x169 * x209,
            x113 * x169 * x213,
            x121 * x176 * x206,
            x119 * x176 * x205,
            x117 * x176 * x209,
            x115 * x176 * x213,
            x113 * x176 * x217,
            x146 * x219 * x72,
            x152 * x219 * x95,
            x152 * x221 * x72,
            x160 * x219 * x75,
            x160 * x221 * x95,
            x160 * x223 * x72,
            x169 * x219 * x88,
            x169 * x221 * x75,
            x169 * x223 * x95,
            x169 * x225 * x72,
            x103 * x176 * x219,
            x176 * x221 * x88,
            x176 * x223 * x75,
            x176 * x225 * x95,
            x176 * x229 * x72,
            x133 * x206 * x58,
            x134 * x206 * x41,
            x133 * x205 * x41,
            x136 * x18 * x206,
            x134 * x18 * x205,
            x133 * x18 * x209,
            x137 * x16 * x206,
            x136 * x16 * x205,
            x134 * x16 * x209,
            x133 * x16 * x213,
            x138 * x198 * x202,
            x12 * x137 * x205,
            x12 * x136 * x209,
            x12 * x134 * x213,
            x12 * x133 * x217,
            x113 * x219 * x58,
            x115 * x219 * x41,
            x113 * x221 * x41,
            x117 * x18 * x219,
            x115 * x18 * x221,
            x113 * x18 * x223,
            x119 * x16 * x219,
            x117 * x16 * x221,
            x115 * x16 * x223,
            x113 * x16 * x225,
            x12 * x121 * x219,
            x119 * x12 * x221,
            x117 * x12 * x223,
            x115 * x12 * x225,
            x111 * x201 * x229,
            x230 * x58 * x72,
            x230 * x41 * x95,
            x231 * x41 * x72,
            x18 * x230 * x75,
            x18 * x231 * x95,
            x18 * x233 * x72,
            x16 * x230 * x88,
            x16 * x231 * x75,
            x16 * x233 * x95,
            x16 * x234 * x72,
            x103 * x12 * x230,
            x12 * x231 * x88,
            x12 * x233 * x75,
            x10 * x234 * x235,
            x10 * x236,
            x238 * x239 * x80,
            x243 * x244 * x80,
            x238 * x244 * x94,
            x248 * x249 * x80,
            x243 * x249 * x94,
            x238 * x249 * x84,
            x253 * x255 * x80,
            x248 * x255 * x94,
            x243 * x255 * x84,
            x238 * x255 * x99,
            x259 * x260 * x80,
            x253 * x260 * x94,
            x248 * x260 * x84,
            x243 * x260 * x99,
            x110 * x238 * x260,
            x148 * x262 * x80,
            x154 * x264 * x80,
            x154 * x262 * x94,
            x162 * x266 * x80,
            x162 * x264 * x94,
            x162 * x262 * x84,
            x254 * x268 * x80,
            x254 * x266 * x94,
            x254 * x264 * x84,
            x254 * x262 * x99,
            x2 * x271,
            x197 * x2 * x268,
            x175 * x266 * x84,
            x175 * x264 * x99,
            x110 * x175 * x262,
            x124 * x148 * x238,
            x124 * x154 * x243,
            x126 * x154 * x238,
            x124 * x162 * x248,
            x126 * x162 * x243,
            x128 * x162 * x238,
            x124 * x253 * x254,
            x126 * x248 * x254,
            x128 * x243 * x254,
            x130 * x238 * x254,
            x122 * x259 * x272,
            x126 * x175 * x253,
            x128 * x175 * x248,
            x130 * x175 * x243,
            x132 * x175 * x238,
            x273 * x55 * x80,
            x274 * x48 * x80,
            x273 * x48 * x94,
            x275 * x46 * x80,
            x274 * x46 * x94,
            x273 * x46 * x84,
            x276 * x277,
            x275 * x277 * x71,
            x274 * x84 * x9,
            x273 * x9 * x99,
            x195
            * (x111 * x270 + x3 * (2 * x187 + 2 * x188 + 5 * x257 + x258 + 4 * x267)),
            x197 * x276,
            x11 * x275 * x84,
            x11 * x274 * x99,
            x11 * x110 * x273,
            x124 * x262 * x55,
            x124 * x264 * x48,
            x126 * x262 * x48,
            x124 * x266 * x46,
            x126 * x264 * x46,
            x128 * x262 * x46,
            x122 * x268 * x277,
            x126 * x266 * x9,
            x128 * x264 * x9,
            x130 * x262 * x9,
            x122 * x271,
            x11 * x126 * x268,
            x11 * x128 * x266,
            x11 * x130 * x264,
            x11 * x132 * x262,
            x139 * x238 * x55,
            x139 * x243 * x48,
            x140 * x238 * x48,
            x139 * x248 * x46,
            x140 * x243 * x46,
            x142 * x238 * x46,
            x139 * x253 * x9,
            x140 * x248 * x9,
            x142 * x243 * x9,
            x143 * x238 * x9,
            x11 * x139 * x259,
            x11 * x140 * x253,
            x11 * x142 * x248,
            x11 * x143 * x243,
            x11 * x144 * x238,
            x145 * x203 * x239,
            x151 * x206 * x244,
            x156 * x205 * x244,
            x159 * x206 * x249,
            x151 * x205 * x249,
            x156 * x209 * x249,
            x167 * x206 * x255,
            x159 * x205 * x255,
            x151 * x209 * x255,
            x156 * x213 * x255,
            x174 * x206 * x260,
            x167 * x205 * x260,
            x159 * x209 * x260,
            x151 * x213 * x260,
            x156 * x217 * x260,
            x148 * x179 * x206,
            x154 * x181 * x206,
            x154 * x179 * x205,
            x162 * x183 * x206,
            x162 * x181 * x205,
            x162 * x179 * x209,
            x185 * x206 * x254,
            x183 * x205 * x254,
            x181 * x209 * x254,
            x179 * x213 * x254,
            x189 * x202 * x272,
            x175 * x185 * x205,
            x175 * x183 * x209,
            x175 * x181 * x213,
            x175 * x179 * x217,
            x148 * x156 * x219,
            x151 * x154 * x219,
            x154 * x156 * x221,
            x159 * x162 * x219,
            x151 * x162 * x221,
            x156 * x162 * x223,
            x167 * x219 * x254,
            x159 * x221 * x254,
            x151 * x223 * x254,
            x156 * x225 * x254,
            x174 * x175 * x219,
            x167 * x175 * x221,
            x159 * x175 * x223,
            x151 * x175 * x225,
            x145 * x229 * x278,
            x190 * x206 * x55,
            x191 * x206 * x48,
            x190 * x205 * x48,
            x193 * x206 * x46,
            x191 * x205 * x46,
            x190 * x209 * x46,
            x194 * x202 * x277,
            x193 * x205 * x9,
            x191 * x209 * x9,
            x190 * x213 * x9,
            x196 * x202,
            x11 * x194 * x205,
            x11 * x193 * x209,
            x11 * x191 * x213,
            x11 * x190 * x217,
            x179 * x219 * x55,
            x181 * x219 * x48,
            x179 * x221 * x48,
            x183 * x219 * x46,
            x181 * x221 * x46,
            x179 * x223 * x46,
            x185 * x219 * x9,
            x183 * x221 * x9,
            x181 * x223 * x9,
            x179 * x225 * x9,
            x11 * x189 * x219,
            x11 * x185 * x221,
            x11 * x183 * x223,
            x11 * x181 * x225,
            x11 * x179 * x229,
            x156 * x230 * x55,
            x151 * x230 * x48,
            x156 * x231 * x48,
            x159 * x230 * x46,
            x151 * x231 * x46,
            x156 * x233 * x46,
            x167 * x230 * x9,
            x159 * x231 * x9,
            x151 * x233 * x9,
            x145 * x234 * x279,
            x11 * x174 * x230,
            x11 * x167 * x231,
            x11 * x159 * x233,
            x11 * x151 * x234,
            x145 * x236,
            x239 * x281 * x72,
            x244 * x281 * x95,
            x244 * x285 * x72,
            x249 * x281 * x75,
            x249 * x285 * x95,
            x249 * x289 * x72,
            x255 * x281 * x88,
            x255 * x285 * x75,
            x255 * x289 * x95,
            x255 * x293 * x72,
            x103 * x260 * x281,
            x260 * x285 * x88,
            x260 * x289 * x75,
            x260 * x293 * x95,
            x260 * x297 * x72,
            x113 * x148 * x281,
            x115 * x154 * x281,
            x113 * x154 * x285,
            x117 * x162 * x281,
            x115 * x162 * x285,
            x113 * x162 * x289,
            x119 * x254 * x281,
            x117 * x254 * x285,
            x115 * x254 * x289,
            x113 * x254 * x293,
            x121 * x175 * x281,
            x119 * x175 * x285,
            x117 * x175 * x289,
            x115 * x175 * x293,
            x111 * x278 * x297,
            x148 * x299 * x72,
            x154 * x299 * x95,
            x154 * x301 * x72,
            x162 * x299 * x75,
            x162 * x301 * x95,
            x162 * x303 * x72,
            x254 * x299 * x88,
            x254 * x301 * x75,
            x254 * x303 * x95,
            x254 * x305 * x72,
            x103 * x175 * x299,
            x175 * x301 * x88,
            x175 * x303 * x75,
            x2 * x235 * x305,
            x2 * x308,
            x133 * x281 * x55,
            x134 * x281 * x48,
            x133 * x285 * x48,
            x136 * x281 * x46,
            x134 * x285 * x46,
            x133 * x289 * x46,
            x137 * x281 * x9,
            x136 * x285 * x9,
            x134 * x289 * x9,
            x133 * x293 * x9,
            x11 * x138 * x281,
            x11 * x137 * x285,
            x11 * x136 * x289,
            x11 * x134 * x293,
            x11 * x133 * x297,
            x113 * x299 * x55,
            x115 * x299 * x48,
            x113 * x301 * x48,
            x117 * x299 * x46,
            x115 * x301 * x46,
            x113 * x303 * x46,
            x119 * x299 * x9,
            x117 * x301 * x9,
            x115 * x303 * x9,
            x111 * x279 * x305,
            x11 * x121 * x299,
            x11 * x119 * x301,
            x11 * x117 * x303,
            x11 * x115 * x305,
            x111 * x308,
            x309 * x55 * x72,
            x309 * x48 * x95,
            x310 * x48 * x72,
            x309 * x46 * x75,
            x310 * x46 * x95,
            x311 * x46 * x72,
            x309 * x88 * x9,
            x310 * x75 * x9,
            x279 * x311 * x65,
            x279 * x312,
            x103 * x11 * x309,
            x11 * x310 * x88,
            x11 * x311 * x75,
            x235 * x312,
            x200
            * (x122 * x307 + x3 * (2 * x227 + 2 * x228 + 5 * x295 + x296 + 4 * x304)),
        ]
    )


def quadrupole3d_30(a, A, b, B, C):
    """Cartesian 3D (fs) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - C[0]
    x5 = a * b * x0
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x8 = x6 * x7
    x9 = x4 ** 2 * x8
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
    x29 = -x28 - A[1]
    x30 = x20 * x27
    x31 = -x0 * (a * A[2] + b * B[2])
    x32 = -x31 - A[2]
    x33 = x24 * x7
    x34 = x3 * x33
    x35 = x29 ** 2 * x33
    x36 = x34 + x35
    x37 = x25 * x7
    x38 = x27 * x32
    x39 = x3 * x37
    x40 = x32 ** 2 * x37
    x41 = x39 + x40
    x42 = 2 * x34
    x43 = x29 * x36 + x29 * x42
    x44 = x32 * x37
    x45 = x29 * x33
    x46 = 2 * x39
    x47 = x32 * x41 + x32 * x46
    x48 = -x28 - C[1]
    x49 = x2 ** 2 * x8
    x50 = x21 + x23
    x51 = x27 * (x2 * x50 + x3 * (x14 + x49))
    x52 = x45 * x48
    x53 = x34 + x52
    x54 = x33 * x48
    x55 = x3 * (x45 + x54)
    x56 = x29 * x53
    x57 = x55 + x56
    x58 = 3 * x34 + 2 * x52
    x59 = x26 * x6
    x60 = x59 * (x29 * x57 + x3 * (x35 + x58))
    x61 = x32 * x59
    x62 = numpy.pi * x0 * x24
    x63 = x6 * x62
    x64 = -x31 - C[2]
    x65 = x27 * x64
    x66 = x44 * x64
    x67 = x39 + x66
    x68 = x37 * x64
    x69 = x3 * (x44 + x68)
    x70 = x32 * x67
    x71 = x69 + x70
    x72 = x29 * x63
    x73 = 3 * x39 + 2 * x66
    x74 = x63 * (x3 * (x40 + x73) + x32 * x71)
    x75 = x33 * x48 ** 2
    x76 = x34 + x75
    x77 = x10 + x49
    x78 = 2 * x12 * x3 + x2 * x77
    x79 = x29 * x76
    x80 = x42 * x48 + x79
    x81 = x29 * x80 + x3 * (x58 + x75)
    x82 = x11 * x26
    x83 = x11 * x62
    x84 = x37 * x64 ** 2
    x85 = x39 + x84
    x86 = x32 * x85
    x87 = x46 * x64 + x86
    x88 = x3 * (x73 + x84) + x32 * x87

    # 60 item(s)
    return numpy.array(
        [
            x27 * (x2 * x20 + x3 * (4 * x16 + 2 * x18 + 2 * x21 + 2 * x23)),
            x29 * x30,
            x30 * x32,
            x19 * x36 * x37,
            x19 * x29 * x38,
            x19 * x33 * x41,
            x17 * x37 * x43,
            x17 * x36 * x44,
            x17 * x41 * x45,
            x17 * x33 * x47,
            x48 * x51,
            x37 * x50 * x53,
            x38 * x48 * x50,
            x22 * x37 * x57,
            x22 * x44 * x53,
            x22 * x41 * x54,
            x4 * x60,
            x4 * x57 * x61,
            x15 * x41 * x53,
            x4 * x47 * x48 * x63,
            x51 * x64,
            x29 * x50 * x65,
            x33 * x50 * x67,
            x22 * x36 * x68,
            x22 * x45 * x67,
            x22 * x33 * x71,
            x4 * x43 * x59 * x64,
            x15 * x36 * x67,
            x4 * x71 * x72,
            x4 * x74,
            x37 * x76 * x78,
            x37 * x77 * x80,
            x44 * x76 * x77,
            x81 * x82,
            x32 * x80 * x82,
            x12 * x41 * x76,
            x59 * (x29 * x81 + x3 * (4 * x34 * x48 + 2 * x55 + 2 * x56 + 2 * x79)),
            x61 * x81,
            x41 * x8 * x80,
            x47 * x76 * x8,
            x48 * x65 * x78,
            x53 * x68 * x77,
            x54 * x67 * x77,
            x57 * x64 * x82,
            x12 * x53 * x67,
            x48 * x71 * x83,
            x60 * x64,
            x57 * x67 * x8,
            x53 * x71 * x8,
            x48 * x74,
            x33 * x78 * x85,
            x45 * x77 * x85,
            x33 * x77 * x87,
            x12 * x36 * x85,
            x29 * x83 * x87,
            x83 * x88,
            x43 * x8 * x85,
            x36 * x8 * x87,
            x72 * x88,
            x63 * (x3 * (4 * x39 * x64 + 2 * x69 + 2 * x70 + 2 * x86) + x32 * x88),
        ]
    )


def quadrupole3d_31(a, A, b, B, C):
    """Cartesian 3D (fp) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - C[0]
    x5 = -x1 - B[0]
    x6 = a * b * x0
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x9 = x7 * x8
    x10 = x5 * x9
    x11 = x10 * x4
    x12 = x4 ** 2 * x9
    x13 = x3 * x9
    x14 = 3 * x13
    x15 = x12 + x14
    x16 = x3 * (2 * x11 + x15)
    x17 = x13 * x4
    x18 = 2 * x17
    x19 = x12 + x13
    x20 = x19 * x5
    x21 = x18 + x20
    x22 = x2 * x21
    x23 = x16 + x22
    x24 = x19 * x2
    x25 = 4 * x17
    x26 = x4 * x9
    x27 = x3 * (x10 + x26)
    x28 = x11 + x13
    x29 = x2 * x28
    x30 = 2 * x27 + 2 * x29
    x31 = x2 * x23 + x3 * (x20 + x24 + x25 + x30)
    x32 = x10 * x2
    x33 = x2 * x26
    x34 = x3 * (x11 + x14 + x32 + x33)
    x35 = x27 + x29
    x36 = x2 * x35
    x37 = 2 * x33
    x38 = x18 + x24
    x39 = x2 * x38 + x3 * (x15 + x37)
    x40 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x41 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x42 = numpy.pi * x0 * x41
    x43 = x40 * x42
    x44 = -x0 * (a * A[1] + b * B[1])
    x45 = -x44 - B[1]
    x46 = x2 * x9
    x47 = x3 * (x26 + x46)
    x48 = x13 + x33
    x49 = x2 * x48
    x50 = x43 * (x2 * x39 + x3 * (2 * x24 + x25 + 2 * x47 + 2 * x49))
    x51 = -x0 * (a * A[2] + b * B[2])
    x52 = -x51 - B[2]
    x53 = -x44 - A[1]
    x54 = x31 * x43
    x55 = x3 * x8
    x56 = x40 * x55
    x57 = x40 * x8
    x58 = x53 * x57
    x59 = x45 * x58
    x60 = x56 + x59
    x61 = x41 * x8
    x62 = x39 * x43
    x63 = -x51 - A[2]
    x64 = x41 * x55
    x65 = x61 * x63
    x66 = x52 * x65
    x67 = x64 + x66
    x68 = x53 ** 2 * x57
    x69 = x56 + x68
    x70 = x45 * x57
    x71 = x3 * (x58 + x70) + x53 * x60
    x72 = x52 * x61
    x73 = x43 * x63
    x74 = x61 * x63 ** 2
    x75 = x64 + x74
    x76 = x3 * (x65 + x72) + x63 * x67
    x77 = 2 * x56
    x78 = x53 * x69 + x53 * x77
    x79 = 3 * x56
    x80 = x68 + x79
    x81 = x3 * (2 * x59 + x80) + x53 * x71
    x82 = 2 * x64
    x83 = x63 * x75 + x63 * x82
    x84 = 3 * x64
    x85 = x74 + x84
    x86 = x3 * (2 * x66 + x85) + x63 * x76
    x87 = -x44 - C[1]
    x88 = x34 + x36
    x89 = x47 + x49
    x90 = x13 + x32
    x91 = x2 * x90 + x3 * (x10 + x46)
    x92 = x43 * (x2 * x88 + x3 * (x30 + x89 + x91))
    x93 = x70 * x87
    x94 = x56 + x93
    x95 = x2 ** 2 * x9
    x96 = x14 + x95
    x97 = x2 * x89 + x3 * (x37 + x96)
    x98 = x43 * x97
    x99 = x58 * x87
    x100 = x56 + x99
    x101 = x57 * x87
    x102 = x3 * (x101 + x70)
    x103 = x53 * x94
    x104 = x102 + x103
    x105 = x3 * (x101 + x58)
    x106 = x100 * x53
    x107 = x105 + x106
    x108 = x3 * (x59 + x79 + x93 + x99)
    x109 = x104 * x53
    x110 = x108 + x109
    x111 = 2 * x99
    x112 = x107 * x53 + x3 * (x111 + x80)
    x113 = 2 * x102 + 2 * x103
    x114 = x42 * x7
    x115 = x114 * (x110 * x53 + x3 * (x107 + x113 + x71))
    x116 = x114 * x4
    x117 = numpy.pi * x0 * x40 * x7
    x118 = x117 * x4
    x119 = -x51 - C[2]
    x120 = x119 * x72
    x121 = x120 + x64
    x122 = x119 * x43
    x123 = x119 * x61
    x124 = x119 * x65
    x125 = x124 + x64
    x126 = x3 * (x123 + x72)
    x127 = x121 * x63
    x128 = x126 + x127
    x129 = x3 * (x123 + x65)
    x130 = x125 * x63
    x131 = x129 + x130
    x132 = x3 * (x120 + x124 + x66 + x84)
    x133 = x128 * x63
    x134 = x132 + x133
    x135 = 2 * x124
    x136 = x131 * x63 + x3 * (x135 + x85)
    x137 = 2 * x126 + 2 * x127
    x138 = x117 * (x134 * x63 + x3 * (x131 + x137 + x76))
    x139 = x57 * x87 ** 2
    x140 = x139 + x56
    x141 = x2 * x91 + x3 * (2 * x32 + x96)
    x142 = x77 * x87
    x143 = x140 * x45
    x144 = x142 + x143
    x145 = x13 + x95
    x146 = 2 * x13 * x2 + x145 * x2
    x147 = x140 * x53
    x148 = x142 + x147
    x149 = x139 + x79
    x150 = x3 * (x149 + 2 * x93)
    x151 = x144 * x53
    x152 = x150 + x151
    x153 = x148 * x53 + x3 * (x111 + x149)
    x154 = 4 * x56 * x87
    x155 = x152 * x53 + x3 * (x113 + x143 + x147 + x154)
    x156 = x114 * x155
    x157 = x114 * x2
    x158 = x114 * (x153 * x53 + x3 * (2 * x105 + 2 * x106 + 2 * x147 + x154))
    x159 = x114 * x5
    x160 = x117 * x87
    x161 = x119 ** 2 * x61
    x162 = x161 + x64
    x163 = x119 * x82
    x164 = x162 * x52
    x165 = x163 + x164
    x166 = x162 * x63
    x167 = x163 + x166
    x168 = x161 + x84
    x169 = x3 * (2 * x120 + x168)
    x170 = x165 * x63
    x171 = x169 + x170
    x172 = x117 * x2
    x173 = x167 * x63 + x3 * (x135 + x168)
    x174 = 4 * x119 * x64
    x175 = x171 * x63 + x3 * (x137 + x164 + x166 + x174)
    x176 = x117 * x175
    x177 = x117 * (x173 * x63 + x3 * (2 * x129 + 2 * x130 + 2 * x166 + x174))

    # 180 item(s)
    return numpy.array(
        [
            x43 * (x2 * x31 + x3 * (2 * x16 + 2 * x22 + 2 * x34 + 2 * x36 + x39)),
            x45 * x50,
            x50 * x52,
            x53 * x54,
            x39 * x60 * x61,
            x52 * x53 * x62,
            x54 * x63,
            x45 * x62 * x63,
            x39 * x57 * x67,
            x23 * x61 * x69,
            x38 * x61 * x71,
            x38 * x69 * x72,
            x23 * x53 * x73,
            x38 * x60 * x65,
            x38 * x58 * x67,
            x23 * x57 * x75,
            x38 * x70 * x75,
            x38 * x57 * x76,
            x21 * x61 * x78,
            x19 * x61 * x81,
            x19 * x72 * x78,
            x21 * x65 * x69,
            x19 * x65 * x71,
            x19 * x67 * x69,
            x21 * x58 * x75,
            x19 * x60 * x75,
            x19 * x58 * x76,
            x21 * x57 * x83,
            x19 * x70 * x83,
            x19 * x57 * x86,
            x87 * x92,
            x61 * x94 * x97,
            x52 * x87 * x98,
            x100 * x61 * x88,
            x104 * x61 * x89,
            x100 * x72 * x89,
            x73 * x87 * x88,
            x65 * x89 * x94,
            x101 * x67 * x89,
            x107 * x35 * x61,
            x110 * x48 * x61,
            x107 * x48 * x72,
            x100 * x35 * x65,
            x104 * x48 * x65,
            x100 * x48 * x67,
            x101 * x35 * x75,
            x48 * x75 * x94,
            x101 * x48 * x76,
            x112 * x28 * x61,
            x115 * x4,
            x112 * x116 * x52,
            x107 * x28 * x65,
            x110 * x116 * x63,
            x107 * x26 * x67,
            x100 * x28 * x75,
            x104 * x26 * x75,
            x100 * x26 * x76,
            x101 * x28 * x83,
            x26 * x83 * x94,
            x118 * x86 * x87,
            x119 * x92,
            x119 * x45 * x98,
            x121 * x57 * x97,
            x122 * x53 * x88,
            x123 * x60 * x89,
            x121 * x58 * x89,
            x125 * x57 * x88,
            x125 * x70 * x89,
            x128 * x57 * x89,
            x123 * x35 * x69,
            x123 * x48 * x71,
            x121 * x48 * x69,
            x125 * x35 * x58,
            x125 * x48 * x60,
            x128 * x48 * x58,
            x131 * x35 * x57,
            x131 * x48 * x70,
            x134 * x48 * x57,
            x123 * x28 * x78,
            x116 * x119 * x81,
            x121 * x26 * x78,
            x125 * x28 * x69,
            x125 * x26 * x71,
            x128 * x26 * x69,
            x131 * x28 * x58,
            x131 * x26 * x60,
            x118 * x134 * x53,
            x136 * x28 * x57,
            x118 * x136 * x45,
            x138 * x4,
            x140 * x141 * x61,
            x144 * x146 * x61,
            x140 * x146 * x72,
            x148 * x61 * x91,
            x145 * x152 * x61,
            x145 * x148 * x72,
            x140 * x65 * x91,
            x144 * x145 * x65,
            x140 * x145 * x67,
            x153 * x61 * x90,
            x156 * x2,
            x153 * x157 * x52,
            x148 * x65 * x90,
            x152 * x157 * x63,
            x148 * x46 * x67,
            x140 * x75 * x90,
            x144 * x46 * x75,
            x140 * x46 * x76,
            x158 * x5,
            x114 * (x155 * x53 + x3 * (2 * x108 + 2 * x109 + 2 * x150 + 2 * x151 + x153)),
            x158 * x52,
            x153 * x159 * x63,
            x156 * x63,
            x153 * x67 * x9,
            x10 * x148 * x75,
            x152 * x75 * x9,
            x148 * x76 * x9,
            x10 * x140 * x83,
            x144 * x83 * x9,
            x140 * x86 * x9,
            x122 * x141 * x87,
            x123 * x146 * x94,
            x101 * x121 * x146,
            x100 * x123 * x91,
            x104 * x123 * x145,
            x100 * x121 * x145,
            x101 * x125 * x91,
            x125 * x145 * x94,
            x101 * x128 * x145,
            x107 * x123 * x90,
            x110 * x119 * x157,
            x107 * x121 * x46,
            x100 * x125 * x90,
            x104 * x125 * x46,
            x100 * x128 * x46,
            x101 * x131 * x90,
            x131 * x46 * x94,
            x134 * x160 * x2,
            x112 * x119 * x159,
            x115 * x119,
            x112 * x121 * x9,
            x10 * x107 * x125,
            x110 * x125 * x9,
            x107 * x128 * x9,
            x10 * x100 * x131,
            x104 * x131 * x9,
            x100 * x134 * x9,
            x136 * x160 * x5,
            x136 * x9 * x94,
            x138 * x87,
            x141 * x162 * x57,
            x146 * x162 * x70,
            x146 * x165 * x57,
            x162 * x58 * x91,
            x145 * x162 * x60,
            x145 * x165 * x58,
            x167 * x57 * x91,
            x145 * x167 * x70,
            x145 * x171 * x57,
            x162 * x69 * x90,
            x162 * x46 * x71,
            x165 * x46 * x69,
            x167 * x58 * x90,
            x167 * x46 * x60,
            x171 * x172 * x53,
            x173 * x57 * x90,
            x172 * x173 * x45,
            x176 * x2,
            x10 * x162 * x78,
            x162 * x81 * x9,
            x165 * x78 * x9,
            x10 * x167 * x69,
            x167 * x71 * x9,
            x171 * x69 * x9,
            x117 * x173 * x5 * x53,
            x173 * x60 * x9,
            x176 * x53,
            x177 * x5,
            x177 * x45,
            x117 * (x175 * x63 + x3 * (2 * x132 + 2 * x133 + 2 * x169 + 2 * x170 + x173)),
        ]
    )


def quadrupole3d_32(a, A, b, B, C):
    """Cartesian 3D (fd) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - C[0]
    x5 = a * b * x0
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x8 = x6 * x7
    x9 = x4 ** 2 * x8
    x10 = x3 * x8
    x11 = 3 * x10
    x12 = -x1 - B[0]
    x13 = x12 * x6
    x14 = x13 * x7
    x15 = x14 * x4
    x16 = x11 + 2 * x15
    x17 = x3 * (x16 + x9)
    x18 = x4 * x8
    x19 = x18 * x3
    x20 = 2 * x19
    x21 = x10 + x9
    x22 = x12 * x21
    x23 = x20 + x22
    x24 = x12 * x23
    x25 = x17 + x24
    x26 = x2 * x25
    x27 = x10 + x15
    x28 = x12 * x27
    x29 = x3 * (x14 + x18)
    x30 = 2 * x29
    x31 = 4 * x19
    x32 = x30 + x31
    x33 = x3 * (2 * x22 + 2 * x28 + x32)
    x34 = x26 + x33
    x35 = x2 * x23
    x36 = 2 * x35
    x37 = x12 ** 2 * x8
    x38 = x3 * (x16 + x37)
    x39 = x28 + x29
    x40 = x2 * x39
    x41 = 2 * x38 + 2 * x40
    x42 = x2 * x34 + x3 * (3 * x17 + x24 + x36 + x41)
    x43 = x17 + x35
    x44 = x2 * x43
    x45 = x38 + x40
    x46 = x2 * x45
    x47 = x2 * x21
    x48 = x2 * x27
    x49 = 2 * x48
    x50 = x3 * (x22 + x32 + x47 + x49)
    x51 = x14 * x3
    x52 = x10 + x37
    x53 = x2 * x52
    x54 = 2 * x51 + x53
    x55 = x3 * (x28 + 3 * x29 + x49 + x54)
    x56 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x57 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x58 = numpy.pi * x0 * x57
    x59 = x56 * x58
    x60 = -x0 * (a * A[1] + b * B[1])
    x61 = -x60 - B[1]
    x62 = x44 + x50
    x63 = x18 * x2
    x64 = x11 + 2 * x63
    x65 = x20 + x47
    x66 = x2 * x65 + x3 * (x64 + x9)
    x67 = x14 * x2
    x68 = x3 * (x11 + x15 + x63 + x67)
    x69 = x29 + x48
    x70 = x2 * x69
    x71 = 2 * x68 + 2 * x70
    x72 = x59 * (x2 * x62 + x3 * (2 * x17 + x36 + x66 + x71))
    x73 = -x0 * (a * A[2] + b * B[2])
    x74 = -x73 - B[2]
    x75 = x56 * x7
    x76 = x3 * x75
    x77 = x61 ** 2 * x75
    x78 = x76 + x77
    x79 = x2 * x8
    x80 = x3 * (x18 + x79)
    x81 = x10 + x63
    x82 = x2 * x81
    x83 = x2 * x66 + x3 * (x31 + 2 * x47 + 2 * x80 + 2 * x82)
    x84 = x57 * x7
    x85 = x59 * x74
    x86 = x3 * x84
    x87 = x74 ** 2 * x84
    x88 = x86 + x87
    x89 = -x60 - A[1]
    x90 = x42 * x59
    x91 = x75 * x89
    x92 = x61 * x91
    x93 = x76 + x92
    x94 = 2 * x76
    x95 = x78 * x89
    x96 = x61 * x94 + x95
    x97 = x74 * x84
    x98 = -x73 - A[2]
    x99 = x59 * x98
    x100 = x84 * x98
    x101 = x100 * x74
    x102 = x101 + x86
    x103 = x61 * x75
    x104 = 2 * x86
    x105 = x88 * x98
    x106 = x104 * x74 + x105
    x107 = x75 * x89 ** 2
    x108 = x107 + x76
    x109 = x3 * (x103 + x91)
    x110 = x89 * x93
    x111 = x109 + x110
    x112 = 3 * x76
    x113 = x112 + 2 * x92
    x114 = x3 * (x113 + x77) + x89 * x96
    x115 = x84 * x98 ** 2
    x116 = x115 + x86
    x117 = x3 * (x100 + x97)
    x118 = x102 * x98
    x119 = x117 + x118
    x120 = 3 * x86
    x121 = 2 * x101 + x120
    x122 = x106 * x98 + x3 * (x121 + x87)
    x123 = x108 * x89 + x89 * x94
    x124 = x111 * x89 + x3 * (x107 + x113)
    x125 = 4 * x76
    x126 = x114 * x89 + x3 * (2 * x109 + 2 * x110 + x125 * x61 + 2 * x95)
    x127 = x104 * x98 + x116 * x98
    x128 = x119 * x98 + x3 * (x115 + x121)
    x129 = 4 * x86
    x130 = x122 * x98 + x3 * (2 * x105 + 2 * x117 + 2 * x118 + x129 * x74)
    x131 = -x60 - C[1]
    x132 = x46 + x55
    x133 = x11 + 2 * x67
    x134 = x2 * x54 + x3 * (x133 + x37)
    x135 = x59 * (x132 * x2 + x3 * (x134 + x41 + x71))
    x136 = x103 * x131
    x137 = x136 + x76
    x138 = x68 + x70
    x139 = x80 + x82
    x140 = x3 * (x14 + x79)
    x141 = x10 + x67
    x142 = x141 * x2
    x143 = x140 + x142
    x144 = x138 * x2 + x3 * (x139 + x143 + x30 + x49)
    x145 = x131 * x75
    x146 = x3 * (x103 + x145)
    x147 = x137 * x61
    x148 = x146 + x147
    x149 = x2 ** 2 * x8
    x150 = x139 * x2 + x3 * (x149 + x64)
    x151 = x131 * x91
    x152 = x151 + x76
    x153 = x137 * x89
    x154 = x146 + x153
    x155 = x112 + 2 * x136
    x156 = x3 * (x155 + x77)
    x157 = x148 * x89
    x158 = x156 + x157
    x159 = x3 * (x145 + x91)
    x160 = x152 * x89
    x161 = x159 + x160
    x162 = x3 * (x112 + x136 + x151 + x92)
    x163 = x154 * x89
    x164 = x162 + x163
    x165 = x158 * x89
    x166 = 2 * x153
    x167 = x3 * (3 * x146 + x147 + x166 + x96)
    x168 = x165 + x167
    x169 = x112 + 2 * x151
    x170 = x161 * x89 + x3 * (x107 + x169)
    x171 = 2 * x146
    x172 = x164 * x89 + x3 * (x111 + x161 + x166 + x171)
    x173 = 2 * x156 + 2 * x157
    x174 = 2 * x162 + 2 * x163
    x175 = x58 * x6
    x176 = x175 * (x168 * x89 + x3 * (x114 + x173 + x174))
    x177 = x175 * x74
    x178 = x175 * x4
    x179 = numpy.pi * x0 * x56
    x180 = x179 * x6
    x181 = x180 * x4
    x182 = -x73 - C[2]
    x183 = x182 * x59
    x184 = x182 * x97
    x185 = x184 + x86
    x186 = x182 * x84
    x187 = x3 * (x186 + x97)
    x188 = x185 * x74
    x189 = x187 + x188
    x190 = x100 * x182
    x191 = x190 + x86
    x192 = x185 * x98
    x193 = x187 + x192
    x194 = x120 + 2 * x184
    x195 = x3 * (x194 + x87)
    x196 = x189 * x98
    x197 = x195 + x196
    x198 = x3 * (x100 + x186)
    x199 = x191 * x98
    x200 = x198 + x199
    x201 = x3 * (x101 + x120 + x184 + x190)
    x202 = x193 * x98
    x203 = x201 + x202
    x204 = x197 * x98
    x205 = 2 * x192
    x206 = x3 * (x106 + 3 * x187 + x188 + x205)
    x207 = x204 + x206
    x208 = x120 + 2 * x190
    x209 = x200 * x98 + x3 * (x115 + x208)
    x210 = 2 * x187
    x211 = x203 * x98 + x3 * (x119 + x200 + x205 + x210)
    x212 = x180 * x61
    x213 = 2 * x195 + 2 * x196
    x214 = 2 * x201 + 2 * x202
    x215 = x180 * (x207 * x98 + x3 * (x122 + x213 + x214))
    x216 = x131 ** 2 * x75
    x217 = x216 + x76
    x218 = x134 * x2 + x3 * (2 * x140 + 2 * x142 + 4 * x51 + 2 * x53)
    x219 = x131 * x94
    x220 = x217 * x61
    x221 = x219 + x220
    x222 = x143 * x2 + x3 * (x133 + x149)
    x223 = x3 * (x155 + x216)
    x224 = x221 * x61
    x225 = x223 + x224
    x226 = x10 + x149
    x227 = 2 * x10 * x2 + x2 * x226
    x228 = x217 * x89
    x229 = x219 + x228
    x230 = x221 * x89
    x231 = x223 + x230
    x232 = x225 * x89
    x233 = x125 * x131
    x234 = x171 + x233
    x235 = x3 * (2 * x147 + 2 * x220 + x234)
    x236 = x232 + x235
    x237 = x229 * x89 + x3 * (x169 + x216)
    x238 = x231 * x89
    x239 = x3 * (x166 + x220 + x228 + x234)
    x240 = x238 + x239
    x241 = 2 * x230
    x242 = x236 * x89 + x3 * (x173 + 3 * x223 + x224 + x241)
    x243 = x175 * x242
    x244 = x175 * x2
    x245 = x237 * x89 + x3 * (2 * x159 + 2 * x160 + 2 * x228 + x233)
    x246 = x240 * x89 + x3 * (x174 + 2 * x223 + x237 + x241)
    x247 = x13 * x58
    x248 = x180 * x2
    x249 = x13 * x179
    x250 = x182 ** 2 * x84
    x251 = x250 + x86
    x252 = x104 * x182
    x253 = x251 * x74
    x254 = x252 + x253
    x255 = x3 * (x194 + x250)
    x256 = x254 * x74
    x257 = x255 + x256
    x258 = x251 * x98
    x259 = x252 + x258
    x260 = x254 * x98
    x261 = x255 + x260
    x262 = x257 * x98
    x263 = x129 * x182
    x264 = x210 + x263
    x265 = x3 * (2 * x188 + 2 * x253 + x264)
    x266 = x262 + x265
    x267 = x259 * x98 + x3 * (x208 + x250)
    x268 = x261 * x98
    x269 = x3 * (x205 + x253 + x258 + x264)
    x270 = x268 + x269
    x271 = 2 * x260
    x272 = x266 * x98 + x3 * (x213 + 3 * x255 + x256 + x271)
    x273 = x180 * x272
    x274 = x267 * x98 + x3 * (2 * x198 + 2 * x199 + 2 * x258 + x263)
    x275 = x270 * x98 + x3 * (x214 + 2 * x255 + x267 + x271)

    # 360 item(s)
    return numpy.array(
        [
            x59
            * (
                x2 * x42
                + x3 * (2 * x26 + 2 * x33 + 2 * x44 + 2 * x46 + 2 * x50 + 2 * x55)
            ),
            x61 * x72,
            x72 * x74,
            x78 * x83 * x84,
            x61 * x83 * x85,
            x75 * x83 * x88,
            x89 * x90,
            x62 * x84 * x93,
            x62 * x85 * x89,
            x66 * x84 * x96,
            x66 * x93 * x97,
            x66 * x88 * x91,
            x90 * x98,
            x61 * x62 * x99,
            x102 * x62 * x75,
            x100 * x66 * x78,
            x102 * x103 * x66,
            x106 * x66 * x75,
            x108 * x34 * x84,
            x111 * x43 * x84,
            x108 * x43 * x97,
            x114 * x65 * x84,
            x111 * x65 * x97,
            x108 * x65 * x88,
            x34 * x89 * x99,
            x100 * x43 * x93,
            x102 * x43 * x91,
            x100 * x65 * x96,
            x102 * x65 * x93,
            x106 * x65 * x91,
            x116 * x34 * x75,
            x103 * x116 * x43,
            x119 * x43 * x75,
            x116 * x65 * x78,
            x103 * x119 * x65,
            x122 * x65 * x75,
            x123 * x25 * x84,
            x124 * x23 * x84,
            x123 * x23 * x97,
            x126 * x21 * x84,
            x124 * x21 * x97,
            x123 * x21 * x88,
            x100 * x108 * x25,
            x100 * x111 * x23,
            x102 * x108 * x23,
            x100 * x114 * x21,
            x102 * x111 * x21,
            x106 * x108 * x21,
            x116 * x25 * x91,
            x116 * x23 * x93,
            x119 * x23 * x91,
            x116 * x21 * x96,
            x119 * x21 * x93,
            x122 * x21 * x91,
            x127 * x25 * x75,
            x103 * x127 * x23,
            x128 * x23 * x75,
            x127 * x21 * x78,
            x103 * x128 * x21,
            x130 * x21 * x75,
            x131 * x135,
            x137 * x144 * x84,
            x131 * x144 * x85,
            x148 * x150 * x84,
            x137 * x150 * x97,
            x145 * x150 * x88,
            x132 * x152 * x84,
            x138 * x154 * x84,
            x138 * x152 * x97,
            x139 * x158 * x84,
            x139 * x154 * x97,
            x139 * x152 * x88,
            x131 * x132 * x99,
            x100 * x137 * x138,
            x102 * x138 * x145,
            x100 * x139 * x148,
            x102 * x137 * x139,
            x106 * x139 * x145,
            x161 * x45 * x84,
            x164 * x69 * x84,
            x161 * x69 * x97,
            x168 * x81 * x84,
            x164 * x81 * x97,
            x161 * x81 * x88,
            x100 * x152 * x45,
            x100 * x154 * x69,
            x102 * x152 * x69,
            x100 * x158 * x81,
            x102 * x154 * x81,
            x106 * x152 * x81,
            x116 * x145 * x45,
            x116 * x137 * x69,
            x119 * x145 * x69,
            x116 * x148 * x81,
            x119 * x137 * x81,
            x122 * x145 * x81,
            x170 * x39 * x84,
            x172 * x27 * x84,
            x170 * x27 * x97,
            x176 * x4,
            x172 * x177 * x4,
            x170 * x18 * x88,
            x100 * x161 * x39,
            x100 * x164 * x27,
            x102 * x161 * x27,
            x168 * x178 * x98,
            x102 * x164 * x18,
            x106 * x161 * x18,
            x116 * x152 * x39,
            x116 * x154 * x27,
            x119 * x152 * x27,
            x116 * x158 * x18,
            x119 * x154 * x18,
            x122 * x152 * x18,
            x127 * x145 * x39,
            x127 * x137 * x27,
            x128 * x145 * x27,
            x127 * x148 * x18,
            x128 * x137 * x18,
            x130 * x131 * x181,
            x135 * x182,
            x144 * x183 * x61,
            x144 * x185 * x75,
            x150 * x186 * x78,
            x103 * x150 * x185,
            x150 * x189 * x75,
            x132 * x183 * x89,
            x138 * x186 * x93,
            x138 * x185 * x91,
            x139 * x186 * x96,
            x139 * x185 * x93,
            x139 * x189 * x91,
            x132 * x191 * x75,
            x103 * x138 * x191,
            x138 * x193 * x75,
            x139 * x191 * x78,
            x103 * x139 * x193,
            x139 * x197 * x75,
            x108 * x186 * x45,
            x111 * x186 * x69,
            x108 * x185 * x69,
            x114 * x186 * x81,
            x111 * x185 * x81,
            x108 * x189 * x81,
            x191 * x45 * x91,
            x191 * x69 * x93,
            x193 * x69 * x91,
            x191 * x81 * x96,
            x193 * x81 * x93,
            x197 * x81 * x91,
            x200 * x45 * x75,
            x103 * x200 * x69,
            x203 * x69 * x75,
            x200 * x78 * x81,
            x103 * x203 * x81,
            x207 * x75 * x81,
            x123 * x186 * x39,
            x124 * x186 * x27,
            x123 * x185 * x27,
            x126 * x178 * x182,
            x124 * x18 * x185,
            x123 * x18 * x189,
            x108 * x191 * x39,
            x111 * x191 * x27,
            x108 * x193 * x27,
            x114 * x18 * x191,
            x111 * x18 * x193,
            x108 * x18 * x197,
            x200 * x39 * x91,
            x200 * x27 * x93,
            x203 * x27 * x91,
            x18 * x200 * x96,
            x18 * x203 * x93,
            x181 * x207 * x89,
            x209 * x39 * x75,
            x103 * x209 * x27,
            x211 * x27 * x75,
            x18 * x209 * x78,
            x211 * x212 * x4,
            x215 * x4,
            x217 * x218 * x84,
            x221 * x222 * x84,
            x217 * x222 * x97,
            x225 * x227 * x84,
            x221 * x227 * x97,
            x217 * x227 * x88,
            x134 * x229 * x84,
            x143 * x231 * x84,
            x143 * x229 * x97,
            x226 * x236 * x84,
            x226 * x231 * x97,
            x226 * x229 * x88,
            x100 * x134 * x217,
            x100 * x143 * x221,
            x102 * x143 * x217,
            x100 * x225 * x226,
            x102 * x221 * x226,
            x106 * x217 * x226,
            x237 * x54 * x84,
            x141 * x240 * x84,
            x141 * x237 * x97,
            x2 * x243,
            x177 * x2 * x240,
            x237 * x79 * x88,
            x100 * x229 * x54,
            x100 * x141 * x231,
            x102 * x141 * x229,
            x236 * x244 * x98,
            x102 * x231 * x79,
            x106 * x229 * x79,
            x116 * x217 * x54,
            x116 * x141 * x221,
            x119 * x141 * x217,
            x116 * x225 * x79,
            x119 * x221 * x79,
            x122 * x217 * x79,
            x245 * x52 * x84,
            x246 * x247,
            x245 * x247 * x74,
            x175
            * (
                x242 * x89
                + x3 * (2 * x165 + 2 * x167 + 2 * x232 + 2 * x235 + 2 * x238 + 2 * x239)
            ),
            x177 * x246,
            x245 * x8 * x88,
            x100 * x237 * x52,
            x240 * x247 * x98,
            x102 * x14 * x237,
            x243 * x98,
            x102 * x240 * x8,
            x106 * x237 * x8,
            x116 * x229 * x52,
            x116 * x14 * x231,
            x119 * x14 * x229,
            x116 * x236 * x8,
            x119 * x231 * x8,
            x122 * x229 * x8,
            x127 * x217 * x52,
            x127 * x14 * x221,
            x128 * x14 * x217,
            x127 * x225 * x8,
            x128 * x221 * x8,
            x130 * x217 * x8,
            x131 * x183 * x218,
            x137 * x186 * x222,
            x145 * x185 * x222,
            x148 * x186 * x227,
            x137 * x185 * x227,
            x145 * x189 * x227,
            x134 * x152 * x186,
            x143 * x154 * x186,
            x143 * x152 * x185,
            x158 * x186 * x226,
            x154 * x185 * x226,
            x152 * x189 * x226,
            x134 * x145 * x191,
            x137 * x143 * x191,
            x143 * x145 * x193,
            x148 * x191 * x226,
            x137 * x193 * x226,
            x145 * x197 * x226,
            x161 * x186 * x54,
            x141 * x164 * x186,
            x141 * x161 * x185,
            x168 * x182 * x244,
            x164 * x185 * x79,
            x161 * x189 * x79,
            x152 * x191 * x54,
            x141 * x154 * x191,
            x141 * x152 * x193,
            x158 * x191 * x79,
            x154 * x193 * x79,
            x152 * x197 * x79,
            x145 * x200 * x54,
            x137 * x141 * x200,
            x141 * x145 * x203,
            x148 * x200 * x79,
            x137 * x203 * x79,
            x131 * x207 * x248,
            x170 * x186 * x52,
            x172 * x182 * x247,
            x14 * x170 * x185,
            x176 * x182,
            x172 * x185 * x8,
            x170 * x189 * x8,
            x161 * x191 * x52,
            x14 * x164 * x191,
            x14 * x161 * x193,
            x168 * x191 * x8,
            x164 * x193 * x8,
            x161 * x197 * x8,
            x152 * x200 * x52,
            x14 * x154 * x200,
            x14 * x152 * x203,
            x158 * x200 * x8,
            x154 * x203 * x8,
            x152 * x207 * x8,
            x145 * x209 * x52,
            x137 * x14 * x209,
            x131 * x211 * x249,
            x148 * x209 * x8,
            x137 * x211 * x8,
            x131 * x215,
            x218 * x251 * x75,
            x103 * x222 * x251,
            x222 * x254 * x75,
            x227 * x251 * x78,
            x103 * x227 * x254,
            x227 * x257 * x75,
            x134 * x251 * x91,
            x143 * x251 * x93,
            x143 * x254 * x91,
            x226 * x251 * x96,
            x226 * x254 * x93,
            x226 * x257 * x91,
            x134 * x259 * x75,
            x103 * x143 * x259,
            x143 * x261 * x75,
            x226 * x259 * x78,
            x103 * x226 * x261,
            x226 * x266 * x75,
            x108 * x251 * x54,
            x111 * x141 * x251,
            x108 * x141 * x254,
            x114 * x251 * x79,
            x111 * x254 * x79,
            x108 * x257 * x79,
            x259 * x54 * x91,
            x141 * x259 * x93,
            x141 * x261 * x91,
            x259 * x79 * x96,
            x261 * x79 * x93,
            x248 * x266 * x89,
            x267 * x54 * x75,
            x103 * x141 * x267,
            x141 * x270 * x75,
            x267 * x78 * x79,
            x2 * x212 * x270,
            x2 * x273,
            x123 * x251 * x52,
            x124 * x14 * x251,
            x123 * x14 * x254,
            x126 * x251 * x8,
            x124 * x254 * x8,
            x123 * x257 * x8,
            x108 * x259 * x52,
            x111 * x14 * x259,
            x108 * x14 * x261,
            x114 * x259 * x8,
            x111 * x261 * x8,
            x108 * x266 * x8,
            x267 * x52 * x91,
            x14 * x267 * x93,
            x249 * x270 * x89,
            x267 * x8 * x96,
            x270 * x8 * x93,
            x273 * x89,
            x274 * x52 * x75,
            x249 * x274 * x61,
            x249 * x275,
            x274 * x78 * x8,
            x212 * x275,
            x180
            * (
                x272 * x98
                + x3 * (2 * x204 + 2 * x206 + 2 * x262 + 2 * x265 + 2 * x268 + 2 * x269)
            ),
        ]
    )


def quadrupole3d_33(a, A, b, B, C):
    """Cartesian 3D (ff) quadrupole moment integrals.
    The origin is at C.

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
    x9 = x4 * x8
    x10 = -x1 - C[0]
    x11 = x10 * x8
    x12 = x3 * (x11 + x9)
    x13 = x3 * x8
    x14 = x10 * x9
    x15 = x13 + x14
    x16 = x15 * x4
    x17 = x12 + x16
    x18 = x17 * x4
    x19 = 2 * x3
    x20 = x11 * x19
    x21 = x10 ** 2 * x8
    x22 = x13 + x21
    x23 = x22 * x4
    x24 = x20 + x23
    x25 = x24 * x4
    x26 = x4 ** 2 * x8
    x27 = 3 * x13
    x28 = 2 * x14 + x27
    x29 = x3 * (x26 + x28)
    x30 = 2 * x29
    x31 = x3 * (x21 + x28)
    x32 = x30 + 3 * x31
    x33 = x3 * (2 * x18 + 3 * x25 + x32)
    x34 = x25 + x31
    x35 = x34 * x4
    x36 = 2 * x12
    x37 = 4 * x3
    x38 = x11 * x37
    x39 = x36 + x38
    x40 = x3 * (2 * x16 + 2 * x23 + x39)
    x41 = x35 + x40
    x42 = x2 * x41
    x43 = x33 + x42
    x44 = x2 * x34
    x45 = 3 * x12
    x46 = x19 * x9
    x47 = x13 + x26
    x48 = x4 * x47
    x49 = x46 + x48
    x50 = x3 * (3 * x16 + x45 + x49)
    x51 = x18 + x29
    x52 = x2 * x51
    x53 = 2 * x50 + 2 * x52
    x54 = x2 * x43 + x3 * (x35 + 4 * x40 + 3 * x44 + x53)
    x55 = x40 + x44
    x56 = x2 * x55
    x57 = x50 + x52
    x58 = x2 * x57
    x59 = x17 * x2
    x60 = x3 * (3 * x26 + x27)
    x61 = x2 * x49
    x62 = x60 + x61
    x63 = x3 * (x18 + 4 * x29 + 3 * x59 + x62)
    x64 = 2 * x59
    x65 = x2 * x24
    x66 = 2 * x65
    x67 = x3 * (x25 + x32 + x64 + x66)
    x68 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x69 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x70 = numpy.pi * x0 * x69
    x71 = x68 * x70
    x72 = -x0 * (a * A[1] + b * B[1])
    x73 = -x72 - B[1]
    x74 = x56 + x67
    x75 = x31 + x65
    x76 = x2 * x75
    x77 = x29 + x59
    x78 = x2 * x77
    x79 = x2 * x22
    x80 = x15 * x2
    x81 = 2 * x80
    x82 = x3 * (x23 + x39 + x79 + x81)
    x83 = x2 * x47
    x84 = x46 + x83
    x85 = x3 * (x16 + x45 + x81 + x84)
    x86 = x71 * (
        x2 * x74 + x3 * (2 * x40 + 2 * x44 + 2 * x76 + 2 * x78 + 2 * x82 + 2 * x85)
    )
    x87 = -x0 * (a * A[2] + b * B[2])
    x88 = -x87 - B[2]
    x89 = x68 * x7
    x90 = x3 * x89
    x91 = x73 ** 2 * x89
    x92 = x90 + x91
    x93 = x76 + x82
    x94 = x11 * x2
    x95 = x27 + 2 * x94
    x96 = x20 + x79
    x97 = x2 * x96 + x3 * (x21 + x95)
    x98 = x2 * x9
    x99 = x3 * (x14 + x27 + x94 + x98)
    x100 = x12 + x80
    x101 = x100 * x2
    x102 = 2 * x101 + 2 * x99
    x103 = x2 * x93 + x3 * (x102 + 2 * x31 + x66 + x97)
    x104 = x69 * x7
    x105 = x71 * x88
    x106 = x104 * x3
    x107 = x104 * x88 ** 2
    x108 = x106 + x107
    x109 = x73 * x90
    x110 = 2 * x109
    x111 = x73 * x92
    x112 = x110 + x111
    x113 = x2 * x8
    x114 = x3 * (x11 + x113)
    x115 = x13 + x94
    x116 = x115 * x2
    x117 = x2 * x97 + x3 * (2 * x114 + 2 * x116 + x38 + 2 * x79)
    x118 = x104 * x88
    x119 = x73 * x89
    x120 = x106 * x88
    x121 = 2 * x120
    x122 = x108 * x88
    x123 = x121 + x122
    x124 = -x72 - A[1]
    x125 = x54 * x71
    x126 = x124 * x89
    x127 = x126 * x73
    x128 = x127 + x90
    x129 = x124 * x92
    x130 = x110 + x129
    x131 = 3 * x90
    x132 = x3 * (x131 + 3 * x91)
    x133 = x112 * x124
    x134 = x132 + x133
    x135 = -x87 - A[2]
    x136 = x135 * x71
    x137 = x104 * x135
    x138 = x137 * x88
    x139 = x106 + x138
    x140 = x108 * x135
    x141 = x121 + x140
    x142 = 3 * x106
    x143 = x3 * (3 * x107 + x142)
    x144 = x123 * x135
    x145 = x143 + x144
    x146 = x124 ** 2 * x89
    x147 = x146 + x90
    x148 = x3 * (x119 + x126)
    x149 = x124 * x128
    x150 = x148 + x149
    x151 = 2 * x127 + x131
    x152 = x3 * (x151 + x91)
    x153 = x124 * x130
    x154 = x152 + x153
    x155 = x124 * x134 + x3 * (8 * x109 + x111 + 3 * x129)
    x156 = x104 * x135 ** 2
    x157 = x106 + x156
    x158 = x3 * (x118 + x137)
    x159 = x135 * x139
    x160 = x158 + x159
    x161 = 2 * x138 + x142
    x162 = x3 * (x107 + x161)
    x163 = x135 * x141
    x164 = x162 + x163
    x165 = x135 * x145 + x3 * (8 * x120 + x122 + 3 * x140)
    x166 = 2 * x90
    x167 = x124 * x147 + x124 * x166
    x168 = x124 * x150 + x3 * (x146 + x151)
    x169 = x124 * x154 + x3 * (4 * x109 + 2 * x129 + 2 * x148 + 2 * x149)
    x170 = x124 * x155 + x3 * (2 * x132 + 2 * x133 + 3 * x152 + 3 * x153)
    x171 = 2 * x106
    x172 = x135 * x157 + x135 * x171
    x173 = x135 * x160 + x3 * (x156 + x161)
    x174 = x135 * x164 + x3 * (4 * x120 + 2 * x140 + 2 * x158 + 2 * x159)
    x175 = x135 * x165 + x3 * (2 * x143 + 2 * x144 + 3 * x162 + 3 * x163)
    x176 = -x72 - C[1]
    x177 = x58 + x63
    x178 = x2 * x62 + x3 * (8 * x3 * x9 + x48 + 3 * x83)
    x179 = x71 * (x177 * x2 + x3 * (x178 + x53 + 3 * x78 + 3 * x85))
    x180 = x119 * x176
    x181 = x180 + x90
    x182 = x78 + x85
    x183 = x27 + 2 * x98
    x184 = x3 * (x183 + x26)
    x185 = x2 * x84
    x186 = x184 + x185
    x187 = x182 * x2 + x3 * (x102 + x186 + x30 + x64)
    x188 = x176 * x89
    x189 = x3 * (x119 + x188)
    x190 = x181 * x73
    x191 = x189 + x190
    x192 = x101 + x99
    x193 = x114 + x116
    x194 = x3 * (x113 + x9)
    x195 = x13 + x98
    x196 = x195 * x2
    x197 = x194 + x196
    x198 = x192 * x2 + x3 * (x193 + x197 + x36 + x81)
    x199 = x131 + 2 * x180
    x200 = x3 * (x199 + x91)
    x201 = x191 * x73
    x202 = x200 + x201
    x203 = x2 ** 2 * x8
    x204 = x193 * x2 + x3 * (x203 + x95)
    x205 = x126 * x176
    x206 = x205 + x90
    x207 = x124 * x181
    x208 = x189 + x207
    x209 = x124 * x191
    x210 = x200 + x209
    x211 = 3 * x189
    x212 = x3 * (x112 + 3 * x190 + x211)
    x213 = x124 * x202
    x214 = x212 + x213
    x215 = x3 * (x126 + x188)
    x216 = x124 * x206
    x217 = x215 + x216
    x218 = x3 * (x127 + x131 + x180 + x205)
    x219 = x124 * x208
    x220 = x218 + x219
    x221 = x124 * x210
    x222 = 2 * x207
    x223 = x3 * (x130 + x190 + x211 + x222)
    x224 = x221 + x223
    x225 = x124 * x214
    x226 = x3 * (x134 + 4 * x200 + x201 + 3 * x209)
    x227 = x225 + x226
    x228 = x131 + 2 * x205
    x229 = x124 * x217 + x3 * (x146 + x228)
    x230 = 2 * x189
    x231 = x124 * x220 + x3 * (x150 + x217 + x222 + x230)
    x232 = 2 * x200
    x233 = 2 * x209
    x234 = 2 * x218 + 2 * x219
    x235 = x124 * x224 + x3 * (x154 + x232 + x233 + x234)
    x236 = 2 * x212 + 2 * x213
    x237 = x6 * x70
    x238 = x237 * (x124 * x227 + x3 * (x155 + 3 * x221 + 3 * x223 + x236))
    x239 = x10 * x237
    x240 = numpy.pi * x0 * x6 * x68
    x241 = x10 * x240
    x242 = -x87 - C[2]
    x243 = x242 * x71
    x244 = x118 * x242
    x245 = x106 + x244
    x246 = x104 * x242
    x247 = x3 * (x118 + x246)
    x248 = x245 * x88
    x249 = x247 + x248
    x250 = x142 + 2 * x244
    x251 = x3 * (x107 + x250)
    x252 = x249 * x88
    x253 = x251 + x252
    x254 = x137 * x242
    x255 = x106 + x254
    x256 = x135 * x245
    x257 = x247 + x256
    x258 = x135 * x249
    x259 = x251 + x258
    x260 = 3 * x247
    x261 = x3 * (x123 + 3 * x248 + x260)
    x262 = x135 * x253
    x263 = x261 + x262
    x264 = x3 * (x137 + x246)
    x265 = x135 * x255
    x266 = x264 + x265
    x267 = x3 * (x138 + x142 + x244 + x254)
    x268 = x135 * x257
    x269 = x267 + x268
    x270 = x135 * x259
    x271 = 2 * x256
    x272 = x3 * (x141 + x248 + x260 + x271)
    x273 = x270 + x272
    x274 = x135 * x263
    x275 = x3 * (x145 + 4 * x251 + x252 + 3 * x258)
    x276 = x274 + x275
    x277 = x142 + 2 * x254
    x278 = x135 * x266 + x3 * (x156 + x277)
    x279 = 2 * x247
    x280 = x135 * x269 + x3 * (x160 + x266 + x271 + x279)
    x281 = 2 * x251
    x282 = 2 * x258
    x283 = 2 * x267 + 2 * x268
    x284 = x135 * x273 + x3 * (x164 + x281 + x282 + x283)
    x285 = 2 * x261 + 2 * x262
    x286 = x240 * (x135 * x276 + x3 * (x165 + 3 * x270 + 3 * x272 + x285))
    x287 = x176 ** 2 * x89
    x288 = x287 + x90
    x289 = x178 * x2 + x3 * (3 * x184 + 3 * x185 + 2 * x60 + 2 * x61)
    x290 = x166 * x176
    x291 = x288 * x73
    x292 = x290 + x291
    x293 = x186 * x2 + x3 * (2 * x194 + 2 * x196 + x37 * x9 + 2 * x83)
    x294 = x3 * (x199 + x287)
    x295 = x292 * x73
    x296 = x294 + x295
    x297 = x197 * x2 + x3 * (x183 + x203)
    x298 = x296 * x73
    x299 = 4 * x176 * x90
    x300 = x230 + x299
    x301 = x3 * (2 * x190 + 2 * x291 + x300)
    x302 = x298 + x301
    x303 = x13 + x203
    x304 = 2 * x13 * x2 + x2 * x303
    x305 = x124 * x288
    x306 = x290 + x305
    x307 = x124 * x292
    x308 = x294 + x307
    x309 = x124 * x296
    x310 = x301 + x309
    x311 = x232 + 3 * x294
    x312 = x3 * (2 * x201 + 3 * x295 + x311)
    x313 = x124 * x302
    x314 = x312 + x313
    x315 = x124 * x306 + x3 * (x228 + x287)
    x316 = x124 * x308
    x317 = x3 * (x222 + x291 + x300 + x305)
    x318 = x316 + x317
    x319 = x124 * x310
    x320 = 2 * x307
    x321 = x3 * (x233 + x295 + x311 + x320)
    x322 = x319 + x321
    x323 = x124 * x314 + x3 * (x236 + x298 + 4 * x301 + 3 * x309)
    x324 = x237 * x323
    x325 = x2 * x237
    x326 = x124 * x315 + x3 * (2 * x215 + 2 * x216 + x299 + 2 * x305)
    x327 = x124 * x318 + x3 * (x234 + 2 * x294 + x315 + x320)
    x328 = x237 * (
        x124 * x322
        + x3 * (2 * x221 + 2 * x223 + 2 * x301 + 2 * x309 + 2 * x316 + 2 * x317)
    )
    x329 = x237 * x4
    x330 = x176 * x240
    x331 = x104 * x242 ** 2
    x332 = x106 + x331
    x333 = x171 * x242
    x334 = x332 * x88
    x335 = x333 + x334
    x336 = x3 * (x250 + x331)
    x337 = x335 * x88
    x338 = x336 + x337
    x339 = x338 * x88
    x340 = 4 * x106 * x242
    x341 = x279 + x340
    x342 = x3 * (2 * x248 + 2 * x334 + x341)
    x343 = x339 + x342
    x344 = x135 * x332
    x345 = x333 + x344
    x346 = x135 * x335
    x347 = x336 + x346
    x348 = x135 * x338
    x349 = x342 + x348
    x350 = x281 + 3 * x336
    x351 = x3 * (2 * x252 + 3 * x337 + x350)
    x352 = x135 * x343
    x353 = x351 + x352
    x354 = x2 * x240
    x355 = x135 * x345 + x3 * (x277 + x331)
    x356 = x135 * x347
    x357 = x3 * (x271 + x334 + x341 + x344)
    x358 = x356 + x357
    x359 = x135 * x349
    x360 = 2 * x346
    x361 = x3 * (x282 + x337 + x350 + x360)
    x362 = x359 + x361
    x363 = x135 * x353 + x3 * (x285 + x339 + 4 * x342 + 3 * x348)
    x364 = x240 * x363
    x365 = x240 * x4
    x366 = x135 * x355 + x3 * (2 * x264 + 2 * x265 + x340 + 2 * x344)
    x367 = x135 * x358 + x3 * (x283 + 2 * x336 + x355 + x360)
    x368 = x240 * (
        x135 * x362
        + x3 * (2 * x270 + 2 * x272 + 2 * x342 + 2 * x348 + 2 * x356 + 2 * x357)
    )

    # 600 item(s)
    return numpy.array(
        [
            x71
            * (
                x2 * x54
                + x3 * (2 * x33 + 2 * x42 + 3 * x56 + 2 * x58 + 2 * x63 + 3 * x67)
            ),
            x73 * x86,
            x86 * x88,
            x103 * x104 * x92,
            x103 * x105 * x73,
            x103 * x108 * x89,
            x104 * x112 * x117,
            x117 * x118 * x92,
            x108 * x117 * x119,
            x117 * x123 * x89,
            x124 * x125,
            x104 * x128 * x74,
            x105 * x124 * x74,
            x104 * x130 * x93,
            x118 * x128 * x93,
            x108 * x126 * x93,
            x104 * x134 * x97,
            x118 * x130 * x97,
            x108 * x128 * x97,
            x123 * x126 * x97,
            x125 * x135,
            x136 * x73 * x74,
            x139 * x74 * x89,
            x137 * x92 * x93,
            x119 * x139 * x93,
            x141 * x89 * x93,
            x112 * x137 * x97,
            x139 * x92 * x97,
            x119 * x141 * x97,
            x145 * x89 * x97,
            x104 * x147 * x43,
            x104 * x150 * x55,
            x118 * x147 * x55,
            x104 * x154 * x75,
            x118 * x150 * x75,
            x108 * x147 * x75,
            x104 * x155 * x96,
            x118 * x154 * x96,
            x108 * x150 * x96,
            x123 * x147 * x96,
            x124 * x136 * x43,
            x128 * x137 * x55,
            x126 * x139 * x55,
            x130 * x137 * x75,
            x128 * x139 * x75,
            x126 * x141 * x75,
            x134 * x137 * x96,
            x130 * x139 * x96,
            x128 * x141 * x96,
            x126 * x145 * x96,
            x157 * x43 * x89,
            x119 * x157 * x55,
            x160 * x55 * x89,
            x157 * x75 * x92,
            x119 * x160 * x75,
            x164 * x75 * x89,
            x112 * x157 * x96,
            x160 * x92 * x96,
            x119 * x164 * x96,
            x165 * x89 * x96,
            x104 * x167 * x41,
            x104 * x168 * x34,
            x118 * x167 * x34,
            x104 * x169 * x24,
            x118 * x168 * x24,
            x108 * x167 * x24,
            x104 * x170 * x22,
            x118 * x169 * x22,
            x108 * x168 * x22,
            x123 * x167 * x22,
            x137 * x147 * x41,
            x137 * x150 * x34,
            x139 * x147 * x34,
            x137 * x154 * x24,
            x139 * x150 * x24,
            x141 * x147 * x24,
            x137 * x155 * x22,
            x139 * x154 * x22,
            x141 * x150 * x22,
            x145 * x147 * x22,
            x126 * x157 * x41,
            x128 * x157 * x34,
            x126 * x160 * x34,
            x130 * x157 * x24,
            x128 * x160 * x24,
            x126 * x164 * x24,
            x134 * x157 * x22,
            x130 * x160 * x22,
            x128 * x164 * x22,
            x126 * x165 * x22,
            x172 * x41 * x89,
            x119 * x172 * x34,
            x173 * x34 * x89,
            x172 * x24 * x92,
            x119 * x173 * x24,
            x174 * x24 * x89,
            x112 * x172 * x22,
            x173 * x22 * x92,
            x119 * x174 * x22,
            x175 * x22 * x89,
            x176 * x179,
            x104 * x181 * x187,
            x105 * x176 * x187,
            x104 * x191 * x198,
            x118 * x181 * x198,
            x108 * x188 * x198,
            x104 * x202 * x204,
            x118 * x191 * x204,
            x108 * x181 * x204,
            x123 * x188 * x204,
            x104 * x177 * x206,
            x104 * x182 * x208,
            x118 * x182 * x206,
            x104 * x192 * x210,
            x118 * x192 * x208,
            x108 * x192 * x206,
            x104 * x193 * x214,
            x118 * x193 * x210,
            x108 * x193 * x208,
            x123 * x193 * x206,
            x136 * x176 * x177,
            x137 * x181 * x182,
            x139 * x182 * x188,
            x137 * x191 * x192,
            x139 * x181 * x192,
            x141 * x188 * x192,
            x137 * x193 * x202,
            x139 * x191 * x193,
            x141 * x181 * x193,
            x145 * x188 * x193,
            x104 * x217 * x57,
            x104 * x220 * x77,
            x118 * x217 * x77,
            x100 * x104 * x224,
            x100 * x118 * x220,
            x100 * x108 * x217,
            x104 * x115 * x227,
            x115 * x118 * x224,
            x108 * x115 * x220,
            x115 * x123 * x217,
            x137 * x206 * x57,
            x137 * x208 * x77,
            x139 * x206 * x77,
            x100 * x137 * x210,
            x100 * x139 * x208,
            x100 * x141 * x206,
            x115 * x137 * x214,
            x115 * x139 * x210,
            x115 * x141 * x208,
            x115 * x145 * x206,
            x157 * x188 * x57,
            x157 * x181 * x77,
            x160 * x188 * x77,
            x100 * x157 * x191,
            x100 * x160 * x181,
            x100 * x164 * x188,
            x115 * x157 * x202,
            x115 * x160 * x191,
            x115 * x164 * x181,
            x115 * x165 * x188,
            x104 * x229 * x51,
            x104 * x17 * x231,
            x118 * x17 * x229,
            x104 * x15 * x235,
            x118 * x15 * x231,
            x108 * x15 * x229,
            x10 * x238,
            x235 * x239 * x88,
            x108 * x11 * x231,
            x11 * x123 * x229,
            x137 * x217 * x51,
            x137 * x17 * x220,
            x139 * x17 * x217,
            x137 * x15 * x224,
            x139 * x15 * x220,
            x141 * x15 * x217,
            x135 * x227 * x239,
            x11 * x139 * x224,
            x11 * x141 * x220,
            x11 * x145 * x217,
            x157 * x206 * x51,
            x157 * x17 * x208,
            x160 * x17 * x206,
            x15 * x157 * x210,
            x15 * x160 * x208,
            x15 * x164 * x206,
            x11 * x157 * x214,
            x11 * x160 * x210,
            x11 * x164 * x208,
            x11 * x165 * x206,
            x172 * x188 * x51,
            x17 * x172 * x181,
            x17 * x173 * x188,
            x15 * x172 * x191,
            x15 * x173 * x181,
            x15 * x174 * x188,
            x11 * x172 * x202,
            x11 * x173 * x191,
            x11 * x174 * x181,
            x175 * x176 * x241,
            x179 * x242,
            x187 * x243 * x73,
            x187 * x245 * x89,
            x198 * x246 * x92,
            x119 * x198 * x245,
            x198 * x249 * x89,
            x112 * x204 * x246,
            x204 * x245 * x92,
            x119 * x204 * x249,
            x204 * x253 * x89,
            x124 * x177 * x243,
            x128 * x182 * x246,
            x126 * x182 * x245,
            x130 * x192 * x246,
            x128 * x192 * x245,
            x126 * x192 * x249,
            x134 * x193 * x246,
            x130 * x193 * x245,
            x128 * x193 * x249,
            x126 * x193 * x253,
            x177 * x255 * x89,
            x119 * x182 * x255,
            x182 * x257 * x89,
            x192 * x255 * x92,
            x119 * x192 * x257,
            x192 * x259 * x89,
            x112 * x193 * x255,
            x193 * x257 * x92,
            x119 * x193 * x259,
            x193 * x263 * x89,
            x147 * x246 * x57,
            x150 * x246 * x77,
            x147 * x245 * x77,
            x100 * x154 * x246,
            x100 * x150 * x245,
            x100 * x147 * x249,
            x115 * x155 * x246,
            x115 * x154 * x245,
            x115 * x150 * x249,
            x115 * x147 * x253,
            x126 * x255 * x57,
            x128 * x255 * x77,
            x126 * x257 * x77,
            x100 * x130 * x255,
            x100 * x128 * x257,
            x100 * x126 * x259,
            x115 * x134 * x255,
            x115 * x130 * x257,
            x115 * x128 * x259,
            x115 * x126 * x263,
            x266 * x57 * x89,
            x119 * x266 * x77,
            x269 * x77 * x89,
            x100 * x266 * x92,
            x100 * x119 * x269,
            x100 * x273 * x89,
            x112 * x115 * x266,
            x115 * x269 * x92,
            x115 * x119 * x273,
            x115 * x276 * x89,
            x167 * x246 * x51,
            x168 * x17 * x246,
            x167 * x17 * x245,
            x15 * x169 * x246,
            x15 * x168 * x245,
            x15 * x167 * x249,
            x170 * x239 * x242,
            x11 * x169 * x245,
            x11 * x168 * x249,
            x11 * x167 * x253,
            x147 * x255 * x51,
            x150 * x17 * x255,
            x147 * x17 * x257,
            x15 * x154 * x255,
            x15 * x150 * x257,
            x147 * x15 * x259,
            x11 * x155 * x255,
            x11 * x154 * x257,
            x11 * x150 * x259,
            x11 * x147 * x263,
            x126 * x266 * x51,
            x128 * x17 * x266,
            x126 * x17 * x269,
            x130 * x15 * x266,
            x128 * x15 * x269,
            x126 * x15 * x273,
            x11 * x134 * x266,
            x11 * x130 * x269,
            x11 * x128 * x273,
            x124 * x241 * x276,
            x278 * x51 * x89,
            x119 * x17 * x278,
            x17 * x280 * x89,
            x15 * x278 * x92,
            x119 * x15 * x280,
            x15 * x284 * x89,
            x11 * x112 * x278,
            x11 * x280 * x92,
            x241 * x284 * x73,
            x10 * x286,
            x104 * x288 * x289,
            x104 * x292 * x293,
            x118 * x288 * x293,
            x104 * x296 * x297,
            x118 * x292 * x297,
            x108 * x288 * x297,
            x104 * x302 * x304,
            x118 * x296 * x304,
            x108 * x292 * x304,
            x123 * x288 * x304,
            x104 * x178 * x306,
            x104 * x186 * x308,
            x118 * x186 * x306,
            x104 * x197 * x310,
            x118 * x197 * x308,
            x108 * x197 * x306,
            x104 * x303 * x314,
            x118 * x303 * x310,
            x108 * x303 * x308,
            x123 * x303 * x306,
            x137 * x178 * x288,
            x137 * x186 * x292,
            x139 * x186 * x288,
            x137 * x197 * x296,
            x139 * x197 * x292,
            x141 * x197 * x288,
            x137 * x302 * x303,
            x139 * x296 * x303,
            x141 * x292 * x303,
            x145 * x288 * x303,
            x104 * x315 * x62,
            x104 * x318 * x84,
            x118 * x315 * x84,
            x104 * x195 * x322,
            x118 * x195 * x318,
            x108 * x195 * x315,
            x2 * x324,
            x322 * x325 * x88,
            x108 * x113 * x318,
            x113 * x123 * x315,
            x137 * x306 * x62,
            x137 * x308 * x84,
            x139 * x306 * x84,
            x137 * x195 * x310,
            x139 * x195 * x308,
            x141 * x195 * x306,
            x135 * x314 * x325,
            x113 * x139 * x310,
            x113 * x141 * x308,
            x113 * x145 * x306,
            x157 * x288 * x62,
            x157 * x292 * x84,
            x160 * x288 * x84,
            x157 * x195 * x296,
            x160 * x195 * x292,
            x164 * x195 * x288,
            x113 * x157 * x302,
            x113 * x160 * x296,
            x113 * x164 * x292,
            x113 * x165 * x288,
            x104 * x326 * x49,
            x104 * x327 * x47,
            x118 * x326 * x47,
            x328 * x4,
            x327 * x329 * x88,
            x108 * x326 * x9,
            x237
            * (
                x124 * x323
                + x3 * (2 * x225 + 2 * x226 + 2 * x312 + 2 * x313 + 3 * x319 + 3 * x321)
            ),
            x328 * x88,
            x108 * x327 * x8,
            x123 * x326 * x8,
            x137 * x315 * x49,
            x137 * x318 * x47,
            x139 * x315 * x47,
            x135 * x322 * x329,
            x139 * x318 * x9,
            x141 * x315 * x9,
            x135 * x324,
            x139 * x322 * x8,
            x141 * x318 * x8,
            x145 * x315 * x8,
            x157 * x306 * x49,
            x157 * x308 * x47,
            x160 * x306 * x47,
            x157 * x310 * x9,
            x160 * x308 * x9,
            x164 * x306 * x9,
            x157 * x314 * x8,
            x160 * x310 * x8,
            x164 * x308 * x8,
            x165 * x306 * x8,
            x172 * x288 * x49,
            x172 * x292 * x47,
            x173 * x288 * x47,
            x172 * x296 * x9,
            x173 * x292 * x9,
            x174 * x288 * x9,
            x172 * x302 * x8,
            x173 * x296 * x8,
            x174 * x292 * x8,
            x175 * x288 * x8,
            x176 * x243 * x289,
            x181 * x246 * x293,
            x188 * x245 * x293,
            x191 * x246 * x297,
            x181 * x245 * x297,
            x188 * x249 * x297,
            x202 * x246 * x304,
            x191 * x245 * x304,
            x181 * x249 * x304,
            x188 * x253 * x304,
            x178 * x206 * x246,
            x186 * x208 * x246,
            x186 * x206 * x245,
            x197 * x210 * x246,
            x197 * x208 * x245,
            x197 * x206 * x249,
            x214 * x246 * x303,
            x210 * x245 * x303,
            x208 * x249 * x303,
            x206 * x253 * x303,
            x178 * x188 * x255,
            x181 * x186 * x255,
            x186 * x188 * x257,
            x191 * x197 * x255,
            x181 * x197 * x257,
            x188 * x197 * x259,
            x202 * x255 * x303,
            x191 * x257 * x303,
            x181 * x259 * x303,
            x188 * x263 * x303,
            x217 * x246 * x62,
            x220 * x246 * x84,
            x217 * x245 * x84,
            x195 * x224 * x246,
            x195 * x220 * x245,
            x195 * x217 * x249,
            x227 * x242 * x325,
            x113 * x224 * x245,
            x113 * x220 * x249,
            x113 * x217 * x253,
            x206 * x255 * x62,
            x208 * x255 * x84,
            x206 * x257 * x84,
            x195 * x210 * x255,
            x195 * x208 * x257,
            x195 * x206 * x259,
            x113 * x214 * x255,
            x113 * x210 * x257,
            x113 * x208 * x259,
            x113 * x206 * x263,
            x188 * x266 * x62,
            x181 * x266 * x84,
            x188 * x269 * x84,
            x191 * x195 * x266,
            x181 * x195 * x269,
            x188 * x195 * x273,
            x113 * x202 * x266,
            x113 * x191 * x269,
            x113 * x181 * x273,
            x2 * x276 * x330,
            x229 * x246 * x49,
            x231 * x246 * x47,
            x229 * x245 * x47,
            x235 * x242 * x329,
            x231 * x245 * x9,
            x229 * x249 * x9,
            x238 * x242,
            x235 * x245 * x8,
            x231 * x249 * x8,
            x229 * x253 * x8,
            x217 * x255 * x49,
            x220 * x255 * x47,
            x217 * x257 * x47,
            x224 * x255 * x9,
            x220 * x257 * x9,
            x217 * x259 * x9,
            x227 * x255 * x8,
            x224 * x257 * x8,
            x220 * x259 * x8,
            x217 * x263 * x8,
            x206 * x266 * x49,
            x208 * x266 * x47,
            x206 * x269 * x47,
            x210 * x266 * x9,
            x208 * x269 * x9,
            x206 * x273 * x9,
            x214 * x266 * x8,
            x210 * x269 * x8,
            x208 * x273 * x8,
            x206 * x276 * x8,
            x188 * x278 * x49,
            x181 * x278 * x47,
            x188 * x280 * x47,
            x191 * x278 * x9,
            x181 * x280 * x9,
            x284 * x330 * x4,
            x202 * x278 * x8,
            x191 * x280 * x8,
            x181 * x284 * x8,
            x176 * x286,
            x289 * x332 * x89,
            x119 * x293 * x332,
            x293 * x335 * x89,
            x297 * x332 * x92,
            x119 * x297 * x335,
            x297 * x338 * x89,
            x112 * x304 * x332,
            x304 * x335 * x92,
            x119 * x304 * x338,
            x304 * x343 * x89,
            x126 * x178 * x332,
            x128 * x186 * x332,
            x126 * x186 * x335,
            x130 * x197 * x332,
            x128 * x197 * x335,
            x126 * x197 * x338,
            x134 * x303 * x332,
            x130 * x303 * x335,
            x128 * x303 * x338,
            x126 * x303 * x343,
            x178 * x345 * x89,
            x119 * x186 * x345,
            x186 * x347 * x89,
            x197 * x345 * x92,
            x119 * x197 * x347,
            x197 * x349 * x89,
            x112 * x303 * x345,
            x303 * x347 * x92,
            x119 * x303 * x349,
            x303 * x353 * x89,
            x147 * x332 * x62,
            x150 * x332 * x84,
            x147 * x335 * x84,
            x154 * x195 * x332,
            x150 * x195 * x335,
            x147 * x195 * x338,
            x113 * x155 * x332,
            x113 * x154 * x335,
            x113 * x150 * x338,
            x113 * x147 * x343,
            x126 * x345 * x62,
            x128 * x345 * x84,
            x126 * x347 * x84,
            x130 * x195 * x345,
            x128 * x195 * x347,
            x126 * x195 * x349,
            x113 * x134 * x345,
            x113 * x130 * x347,
            x113 * x128 * x349,
            x124 * x353 * x354,
            x355 * x62 * x89,
            x119 * x355 * x84,
            x358 * x84 * x89,
            x195 * x355 * x92,
            x119 * x195 * x358,
            x195 * x362 * x89,
            x112 * x113 * x355,
            x113 * x358 * x92,
            x354 * x362 * x73,
            x2 * x364,
            x167 * x332 * x49,
            x168 * x332 * x47,
            x167 * x335 * x47,
            x169 * x332 * x9,
            x168 * x335 * x9,
            x167 * x338 * x9,
            x170 * x332 * x8,
            x169 * x335 * x8,
            x168 * x338 * x8,
            x167 * x343 * x8,
            x147 * x345 * x49,
            x150 * x345 * x47,
            x147 * x347 * x47,
            x154 * x345 * x9,
            x150 * x347 * x9,
            x147 * x349 * x9,
            x155 * x345 * x8,
            x154 * x347 * x8,
            x150 * x349 * x8,
            x147 * x353 * x8,
            x126 * x355 * x49,
            x128 * x355 * x47,
            x126 * x358 * x47,
            x130 * x355 * x9,
            x128 * x358 * x9,
            x124 * x362 * x365,
            x134 * x355 * x8,
            x130 * x358 * x8,
            x128 * x362 * x8,
            x124 * x364,
            x366 * x49 * x89,
            x119 * x366 * x47,
            x367 * x47 * x89,
            x366 * x9 * x92,
            x365 * x367 * x73,
            x368 * x4,
            x112 * x366 * x8,
            x367 * x8 * x92,
            x368 * x73,
            x240
            * (
                x135 * x363
                + x3 * (2 * x274 + 2 * x275 + 2 * x351 + 2 * x352 + 3 * x359 + 3 * x361)
            ),
        ]
    )


def quadrupole3d_34(a, A, b, B, C):
    """Cartesian 3D (fg) quadrupole moment integrals.
    The origin is at C.

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
    x9 = x4 * x8
    x10 = -x1 - C[0]
    x11 = x10 * x8
    x12 = x3 * (x11 + x9)
    x13 = x3 * x8
    x14 = x10 * x9
    x15 = x13 + x14
    x16 = x15 * x4
    x17 = x12 + x16
    x18 = x17 * x4
    x19 = 2 * x13
    x20 = x10 * x19
    x21 = x10 ** 2 * x8
    x22 = x13 + x21
    x23 = x22 * x4
    x24 = x20 + x23
    x25 = x24 * x4
    x26 = x4 ** 2 * x8
    x27 = 3 * x13
    x28 = 2 * x14 + x27
    x29 = x3 * (x26 + x28)
    x30 = 2 * x29
    x31 = x3 * (x21 + x28)
    x32 = x30 + 3 * x31
    x33 = x3 * (2 * x18 + 3 * x25 + x32)
    x34 = x25 + x31
    x35 = x34 * x4
    x36 = 2 * x12
    x37 = 4 * x13
    x38 = x10 * x37
    x39 = x36 + x38
    x40 = x3 * (2 * x16 + 2 * x23 + x39)
    x41 = x35 + x40
    x42 = x4 * x41
    x43 = x33 + x42
    x44 = x2 * x43
    x45 = x18 + x29
    x46 = x4 * x45
    x47 = 3 * x12
    x48 = x19 * x4
    x49 = x13 + x26
    x50 = x4 * x49
    x51 = x48 + x50
    x52 = x3 * (3 * x16 + x47 + x51)
    x53 = 2 * x52
    x54 = 4 * x40 + x53
    x55 = x3 * (4 * x35 + 2 * x46 + x54)
    x56 = x44 + x55
    x57 = x2 * x41
    x58 = 4 * x29
    x59 = x3 * (3 * x26 + x27)
    x60 = x4 * x51
    x61 = x59 + x60
    x62 = x3 * (4 * x18 + x58 + x61)
    x63 = x46 + x52
    x64 = x2 * x63
    x65 = 2 * x62 + 2 * x64
    x66 = x2 * x56 + x3 * (5 * x33 + x42 + 4 * x57 + x65)
    x67 = x62 + x64
    x68 = x2 * x67
    x69 = x33 + x57
    x70 = x2 * x69
    x71 = x2 * x45
    x72 = 8 * x13 * x4
    x73 = x3 * (4 * x50 + x72)
    x74 = x2 * x61
    x75 = x73 + x74
    x76 = x3 * (x46 + 5 * x52 + 4 * x71 + x75)
    x77 = 2 * x71
    x78 = x2 * x34
    x79 = x3 * (x35 + x54 + x77 + 3 * x78)
    x80 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x81 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x82 = numpy.pi * x0 * x81
    x83 = x80 * x82
    x84 = -x0 * (a * A[1] + b * B[1])
    x85 = -x84 - B[1]
    x86 = x70 + x79
    x87 = x40 + x78
    x88 = x2 * x87
    x89 = x52 + x71
    x90 = x2 * x89
    x91 = x17 * x2
    x92 = x2 * x51
    x93 = x59 + x92
    x94 = x3 * (x18 + x58 + 3 * x91 + x93)
    x95 = 2 * x91
    x96 = x2 * x24
    x97 = 2 * x96
    x98 = x3 * (x25 + x32 + x95 + x97)
    x99 = x83 * (
        x2 * x86 + x3 * (2 * x33 + 2 * x57 + 3 * x88 + 2 * x90 + 2 * x94 + 3 * x98)
    )
    x100 = -x0 * (a * A[2] + b * B[2])
    x101 = -x100 - B[2]
    x102 = x7 * x80
    x103 = x102 * x3
    x104 = x102 * x85 ** 2
    x105 = x103 + x104
    x106 = x88 + x98
    x107 = x31 + x96
    x108 = x107 * x2
    x109 = x29 + x91
    x110 = x109 * x2
    x111 = x2 * x22
    x112 = x15 * x2
    x113 = 2 * x112
    x114 = x3 * (x111 + x113 + x23 + x39)
    x115 = x2 * x49
    x116 = x115 + x48
    x117 = x3 * (x113 + x116 + x16 + x47)
    x118 = x106 * x2 + x3 * (
        2 * x108 + 2 * x110 + 2 * x114 + 2 * x117 + 2 * x40 + 2 * x78
    )
    x119 = x7 * x81
    x120 = x101 * x83
    x121 = x119 * x3
    x122 = x101 ** 2 * x119
    x123 = x121 + x122
    x124 = x103 * x85
    x125 = 2 * x124
    x126 = x105 * x85
    x127 = x125 + x126
    x128 = x108 + x114
    x129 = x11 * x2
    x130 = 2 * x129 + x27
    x131 = x111 + x20
    x132 = x131 * x2 + x3 * (x130 + x21)
    x133 = x2 * x9
    x134 = x3 * (x129 + x133 + x14 + x27)
    x135 = x112 + x12
    x136 = x135 * x2
    x137 = 2 * x134 + 2 * x136
    x138 = x128 * x2 + x3 * (x132 + x137 + 2 * x31 + x97)
    x139 = x101 * x119
    x140 = x102 * x85
    x141 = x101 * x121
    x142 = 2 * x141
    x143 = x101 * x123
    x144 = x142 + x143
    x145 = 3 * x103
    x146 = x3 * (3 * x104 + x145)
    x147 = x127 * x85
    x148 = x146 + x147
    x149 = x2 * x8
    x150 = x3 * (x11 + x149)
    x151 = x129 + x13
    x152 = x151 * x2
    x153 = x132 * x2 + x3 * (2 * x111 + 2 * x150 + 2 * x152 + x38)
    x154 = 3 * x121
    x155 = x3 * (3 * x122 + x154)
    x156 = x101 * x144
    x157 = x155 + x156
    x158 = -x84 - A[1]
    x159 = x66 * x83
    x160 = x102 * x158
    x161 = x160 * x85
    x162 = x103 + x161
    x163 = x105 * x158
    x164 = x125 + x163
    x165 = x127 * x158
    x166 = x146 + x165
    x167 = 8 * x124
    x168 = x3 * (4 * x126 + x167)
    x169 = x148 * x158
    x170 = x168 + x169
    x171 = -x100 - A[2]
    x172 = x171 * x83
    x173 = x119 * x171
    x174 = x101 * x173
    x175 = x121 + x174
    x176 = x123 * x171
    x177 = x142 + x176
    x178 = x144 * x171
    x179 = x155 + x178
    x180 = 8 * x141
    x181 = x3 * (4 * x143 + x180)
    x182 = x157 * x171
    x183 = x181 + x182
    x184 = x102 * x158 ** 2
    x185 = x103 + x184
    x186 = x3 * (x140 + x160)
    x187 = x158 * x162
    x188 = x186 + x187
    x189 = x145 + 2 * x161
    x190 = x3 * (x104 + x189)
    x191 = x158 * x164
    x192 = x190 + x191
    x193 = x3 * (x126 + 3 * x163 + x167)
    x194 = x158 * x166
    x195 = x193 + x194
    x196 = x158 * x170 + x3 * (5 * x146 + x147 + 4 * x165)
    x197 = x119 * x171 ** 2
    x198 = x121 + x197
    x199 = x3 * (x139 + x173)
    x200 = x171 * x175
    x201 = x199 + x200
    x202 = x154 + 2 * x174
    x203 = x3 * (x122 + x202)
    x204 = x171 * x177
    x205 = x203 + x204
    x206 = x3 * (x143 + 3 * x176 + x180)
    x207 = x171 * x179
    x208 = x206 + x207
    x209 = x171 * x183 + x3 * (5 * x155 + x156 + 4 * x178)
    x210 = 2 * x103
    x211 = x158 * x185 + x158 * x210
    x212 = x158 * x188 + x3 * (x184 + x189)
    x213 = x158 * x192 + x3 * (4 * x124 + 2 * x163 + 2 * x186 + 2 * x187)
    x214 = x158 * x195 + x3 * (2 * x146 + 2 * x165 + 3 * x190 + 3 * x191)
    x215 = x158 * x196 + x3 * (2 * x168 + 2 * x169 + 4 * x193 + 4 * x194)
    x216 = 2 * x121
    x217 = x171 * x198 + x171 * x216
    x218 = x171 * x201 + x3 * (x197 + x202)
    x219 = x171 * x205 + x3 * (4 * x141 + 2 * x176 + 2 * x199 + 2 * x200)
    x220 = x171 * x208 + x3 * (2 * x155 + 2 * x178 + 3 * x203 + 3 * x204)
    x221 = x171 * x209 + x3 * (2 * x181 + 2 * x182 + 4 * x206 + 4 * x207)
    x222 = -x84 - C[1]
    x223 = x68 + x76
    x224 = x2 * x75 + x3 * (5 * x59 + x60 + 4 * x92)
    x225 = x83 * (x2 * x223 + x3 * (x224 + x65 + 4 * x90 + 4 * x94))
    x226 = x140 * x222
    x227 = x103 + x226
    x228 = x90 + x94
    x229 = x3 * (3 * x115 + x50 + x72)
    x230 = x2 * x93
    x231 = x229 + x230
    x232 = x2 * x228 + x3 * (3 * x110 + 3 * x117 + x231 + x53 + x77)
    x233 = x102 * x222
    x234 = x3 * (x140 + x233)
    x235 = x227 * x85
    x236 = x234 + x235
    x237 = x110 + x117
    x238 = 2 * x133 + x27
    x239 = x3 * (x238 + x26)
    x240 = x116 * x2
    x241 = x239 + x240
    x242 = x2 * x237 + x3 * (x137 + x241 + x30 + x95)
    x243 = x145 + 2 * x226
    x244 = x3 * (x104 + x243)
    x245 = x236 * x85
    x246 = x244 + x245
    x247 = x134 + x136
    x248 = x150 + x152
    x249 = x3 * (x149 + x9)
    x250 = x13 + x133
    x251 = x2 * x250
    x252 = x249 + x251
    x253 = x2 * x247 + x3 * (x113 + x248 + x252 + x36)
    x254 = 3 * x234
    x255 = x3 * (x127 + 3 * x235 + x254)
    x256 = x246 * x85
    x257 = x255 + x256
    x258 = x2 ** 2 * x8
    x259 = x2 * x248 + x3 * (x130 + x258)
    x260 = x160 * x222
    x261 = x103 + x260
    x262 = x158 * x227
    x263 = x234 + x262
    x264 = x158 * x236
    x265 = x244 + x264
    x266 = x158 * x246
    x267 = x255 + x266
    x268 = 4 * x244
    x269 = x3 * (x148 + 4 * x245 + x268)
    x270 = x158 * x257
    x271 = x269 + x270
    x272 = x3 * (x160 + x233)
    x273 = x158 * x261
    x274 = x272 + x273
    x275 = x3 * (x145 + x161 + x226 + x260)
    x276 = x158 * x263
    x277 = x275 + x276
    x278 = x158 * x265
    x279 = 2 * x262
    x280 = x3 * (x164 + x235 + x254 + x279)
    x281 = x278 + x280
    x282 = x158 * x267
    x283 = x3 * (x166 + x245 + 3 * x264 + x268)
    x284 = x282 + x283
    x285 = x158 * x271
    x286 = x3 * (x170 + 5 * x255 + x256 + 4 * x266)
    x287 = x285 + x286
    x288 = x145 + 2 * x260
    x289 = x158 * x274 + x3 * (x184 + x288)
    x290 = 2 * x234
    x291 = x158 * x277 + x3 * (x188 + x274 + x279 + x290)
    x292 = 2 * x244
    x293 = 2 * x264
    x294 = 2 * x275 + 2 * x276
    x295 = x158 * x281 + x3 * (x192 + x292 + x293 + x294)
    x296 = 2 * x255
    x297 = 2 * x266
    x298 = x158 * x284 + x3 * (x195 + 3 * x278 + 3 * x280 + x296 + x297)
    x299 = 2 * x269 + 2 * x270
    x300 = x6 * x82
    x301 = x300 * (x158 * x287 + x3 * (x196 + 4 * x282 + 4 * x283 + x299))
    x302 = x10 * x300
    x303 = numpy.pi * x0 * x6 * x80
    x304 = x10 * x303
    x305 = -x100 - C[2]
    x306 = x305 * x83
    x307 = x139 * x305
    x308 = x121 + x307
    x309 = x119 * x305
    x310 = x3 * (x139 + x309)
    x311 = x101 * x308
    x312 = x310 + x311
    x313 = x154 + 2 * x307
    x314 = x3 * (x122 + x313)
    x315 = x101 * x312
    x316 = x314 + x315
    x317 = 3 * x310
    x318 = x3 * (x144 + 3 * x311 + x317)
    x319 = x101 * x316
    x320 = x318 + x319
    x321 = x173 * x305
    x322 = x121 + x321
    x323 = x171 * x308
    x324 = x310 + x323
    x325 = x171 * x312
    x326 = x314 + x325
    x327 = x171 * x316
    x328 = x318 + x327
    x329 = 4 * x314
    x330 = x3 * (x157 + 4 * x315 + x329)
    x331 = x171 * x320
    x332 = x330 + x331
    x333 = x3 * (x173 + x309)
    x334 = x171 * x322
    x335 = x333 + x334
    x336 = x3 * (x154 + x174 + x307 + x321)
    x337 = x171 * x324
    x338 = x336 + x337
    x339 = x171 * x326
    x340 = 2 * x323
    x341 = x3 * (x177 + x311 + x317 + x340)
    x342 = x339 + x341
    x343 = x171 * x328
    x344 = x3 * (x179 + x315 + 3 * x325 + x329)
    x345 = x343 + x344
    x346 = x171 * x332
    x347 = x3 * (x183 + 5 * x318 + x319 + 4 * x327)
    x348 = x346 + x347
    x349 = x154 + 2 * x321
    x350 = x171 * x335 + x3 * (x197 + x349)
    x351 = 2 * x310
    x352 = x171 * x338 + x3 * (x201 + x335 + x340 + x351)
    x353 = 2 * x314
    x354 = 2 * x325
    x355 = 2 * x336 + 2 * x337
    x356 = x171 * x342 + x3 * (x205 + x353 + x354 + x355)
    x357 = 2 * x318
    x358 = 2 * x327
    x359 = x171 * x345 + x3 * (x208 + 3 * x339 + 3 * x341 + x357 + x358)
    x360 = 2 * x330 + 2 * x331
    x361 = x303 * (x171 * x348 + x3 * (x209 + 4 * x343 + 4 * x344 + x360))
    x362 = x102 * x222 ** 2
    x363 = x103 + x362
    x364 = x2 * x224 + x3 * (4 * x229 + 4 * x230 + 2 * x73 + 2 * x74)
    x365 = x210 * x222
    x366 = x363 * x85
    x367 = x365 + x366
    x368 = x2 * x231 + x3 * (3 * x239 + 3 * x240 + 2 * x59 + 2 * x92)
    x369 = x3 * (x243 + x362)
    x370 = x367 * x85
    x371 = x369 + x370
    x372 = x2 * x241 + x3 * (2 * x115 + 2 * x249 + 2 * x251 + x37 * x4)
    x373 = x371 * x85
    x374 = 4 * x103 * x222
    x375 = x290 + x374
    x376 = x3 * (2 * x235 + 2 * x366 + x375)
    x377 = x373 + x376
    x378 = x2 * x252 + x3 * (x238 + x258)
    x379 = x292 + 3 * x369
    x380 = x3 * (2 * x245 + 3 * x370 + x379)
    x381 = x377 * x85
    x382 = x380 + x381
    x383 = x13 + x258
    x384 = x19 * x2 + x2 * x383
    x385 = x158 * x363
    x386 = x365 + x385
    x387 = x158 * x367
    x388 = x369 + x387
    x389 = x158 * x371
    x390 = x376 + x389
    x391 = x158 * x377
    x392 = x380 + x391
    x393 = x158 * x382
    x394 = x296 + 4 * x376
    x395 = x3 * (2 * x256 + 4 * x373 + x394)
    x396 = x393 + x395
    x397 = x158 * x386 + x3 * (x288 + x362)
    x398 = x158 * x388
    x399 = x3 * (x279 + x366 + x375 + x385)
    x400 = x398 + x399
    x401 = x158 * x390
    x402 = 2 * x387
    x403 = x3 * (x293 + x370 + x379 + x402)
    x404 = x401 + x403
    x405 = x158 * x392
    x406 = x3 * (x297 + x373 + 3 * x389 + x394)
    x407 = x405 + x406
    x408 = x158 * x396 + x3 * (x299 + 5 * x380 + x381 + 4 * x391)
    x409 = x300 * x408
    x410 = x2 * x300
    x411 = x158 * x397 + x3 * (2 * x272 + 2 * x273 + x374 + 2 * x385)
    x412 = x158 * x400 + x3 * (x294 + 2 * x369 + x397 + x402)
    x413 = x158 * x404 + x3 * (
        2 * x278 + 2 * x280 + 2 * x376 + 2 * x389 + 2 * x398 + 2 * x399
    )
    x414 = x300 * (
        x158 * x407
        + x3 * (2 * x282 + 2 * x283 + 2 * x380 + 2 * x391 + 3 * x401 + 3 * x403)
    )
    x415 = x300 * x4
    x416 = x222 * x303
    x417 = x119 * x305 ** 2
    x418 = x121 + x417
    x419 = x216 * x305
    x420 = x101 * x418
    x421 = x419 + x420
    x422 = x3 * (x313 + x417)
    x423 = x101 * x421
    x424 = x422 + x423
    x425 = x101 * x424
    x426 = 4 * x121 * x305
    x427 = x351 + x426
    x428 = x3 * (2 * x311 + 2 * x420 + x427)
    x429 = x425 + x428
    x430 = x353 + 3 * x422
    x431 = x3 * (2 * x315 + 3 * x423 + x430)
    x432 = x101 * x429
    x433 = x431 + x432
    x434 = x171 * x418
    x435 = x419 + x434
    x436 = x171 * x421
    x437 = x422 + x436
    x438 = x171 * x424
    x439 = x428 + x438
    x440 = x171 * x429
    x441 = x431 + x440
    x442 = x171 * x433
    x443 = x357 + 4 * x428
    x444 = x3 * (2 * x319 + 4 * x425 + x443)
    x445 = x442 + x444
    x446 = x2 * x303
    x447 = x171 * x435 + x3 * (x349 + x417)
    x448 = x171 * x437
    x449 = x3 * (x340 + x420 + x427 + x434)
    x450 = x448 + x449
    x451 = x171 * x439
    x452 = 2 * x436
    x453 = x3 * (x354 + x423 + x430 + x452)
    x454 = x451 + x453
    x455 = x171 * x441
    x456 = x3 * (x358 + x425 + 3 * x438 + x443)
    x457 = x455 + x456
    x458 = x171 * x445 + x3 * (x360 + 5 * x431 + x432 + 4 * x440)
    x459 = x303 * x458
    x460 = x303 * x4
    x461 = x171 * x447 + x3 * (2 * x333 + 2 * x334 + x426 + 2 * x434)
    x462 = x171 * x450 + x3 * (x355 + 2 * x422 + x447 + x452)
    x463 = x171 * x454 + x3 * (
        2 * x339 + 2 * x341 + 2 * x428 + 2 * x438 + 2 * x448 + 2 * x449
    )
    x464 = x303 * (
        x171 * x457
        + x3 * (2 * x343 + 2 * x344 + 2 * x431 + 2 * x440 + 3 * x451 + 3 * x453)
    )

    # 900 item(s)
    return numpy.array(
        [
            x83
            * (
                x2 * x66
                + x3 * (2 * x44 + 2 * x55 + 2 * x68 + 4 * x70 + 2 * x76 + 4 * x79)
            ),
            x85 * x99,
            x101 * x99,
            x105 * x118 * x119,
            x118 * x120 * x85,
            x102 * x118 * x123,
            x119 * x127 * x138,
            x105 * x138 * x139,
            x123 * x138 * x140,
            x102 * x138 * x144,
            x119 * x148 * x153,
            x127 * x139 * x153,
            x105 * x123 * x153,
            x140 * x144 * x153,
            x102 * x153 * x157,
            x158 * x159,
            x119 * x162 * x86,
            x120 * x158 * x86,
            x106 * x119 * x164,
            x106 * x139 * x162,
            x106 * x123 * x160,
            x119 * x128 * x166,
            x128 * x139 * x164,
            x123 * x128 * x162,
            x128 * x144 * x160,
            x119 * x132 * x170,
            x132 * x139 * x166,
            x123 * x132 * x164,
            x132 * x144 * x162,
            x132 * x157 * x160,
            x159 * x171,
            x172 * x85 * x86,
            x102 * x175 * x86,
            x105 * x106 * x173,
            x106 * x140 * x175,
            x102 * x106 * x177,
            x127 * x128 * x173,
            x105 * x128 * x175,
            x128 * x140 * x177,
            x102 * x128 * x179,
            x132 * x148 * x173,
            x127 * x132 * x175,
            x105 * x132 * x177,
            x132 * x140 * x179,
            x102 * x132 * x183,
            x119 * x185 * x56,
            x119 * x188 * x69,
            x139 * x185 * x69,
            x119 * x192 * x87,
            x139 * x188 * x87,
            x123 * x185 * x87,
            x107 * x119 * x195,
            x107 * x139 * x192,
            x107 * x123 * x188,
            x107 * x144 * x185,
            x119 * x131 * x196,
            x131 * x139 * x195,
            x123 * x131 * x192,
            x131 * x144 * x188,
            x131 * x157 * x185,
            x158 * x172 * x56,
            x162 * x173 * x69,
            x160 * x175 * x69,
            x164 * x173 * x87,
            x162 * x175 * x87,
            x160 * x177 * x87,
            x107 * x166 * x173,
            x107 * x164 * x175,
            x107 * x162 * x177,
            x107 * x160 * x179,
            x131 * x170 * x173,
            x131 * x166 * x175,
            x131 * x164 * x177,
            x131 * x162 * x179,
            x131 * x160 * x183,
            x102 * x198 * x56,
            x140 * x198 * x69,
            x102 * x201 * x69,
            x105 * x198 * x87,
            x140 * x201 * x87,
            x102 * x205 * x87,
            x107 * x127 * x198,
            x105 * x107 * x201,
            x107 * x140 * x205,
            x102 * x107 * x208,
            x131 * x148 * x198,
            x127 * x131 * x201,
            x105 * x131 * x205,
            x131 * x140 * x208,
            x102 * x131 * x209,
            x119 * x211 * x43,
            x119 * x212 * x41,
            x139 * x211 * x41,
            x119 * x213 * x34,
            x139 * x212 * x34,
            x123 * x211 * x34,
            x119 * x214 * x24,
            x139 * x213 * x24,
            x123 * x212 * x24,
            x144 * x211 * x24,
            x119 * x215 * x22,
            x139 * x214 * x22,
            x123 * x213 * x22,
            x144 * x212 * x22,
            x157 * x211 * x22,
            x173 * x185 * x43,
            x173 * x188 * x41,
            x175 * x185 * x41,
            x173 * x192 * x34,
            x175 * x188 * x34,
            x177 * x185 * x34,
            x173 * x195 * x24,
            x175 * x192 * x24,
            x177 * x188 * x24,
            x179 * x185 * x24,
            x173 * x196 * x22,
            x175 * x195 * x22,
            x177 * x192 * x22,
            x179 * x188 * x22,
            x183 * x185 * x22,
            x160 * x198 * x43,
            x162 * x198 * x41,
            x160 * x201 * x41,
            x164 * x198 * x34,
            x162 * x201 * x34,
            x160 * x205 * x34,
            x166 * x198 * x24,
            x164 * x201 * x24,
            x162 * x205 * x24,
            x160 * x208 * x24,
            x170 * x198 * x22,
            x166 * x201 * x22,
            x164 * x205 * x22,
            x162 * x208 * x22,
            x160 * x209 * x22,
            x102 * x217 * x43,
            x140 * x217 * x41,
            x102 * x218 * x41,
            x105 * x217 * x34,
            x140 * x218 * x34,
            x102 * x219 * x34,
            x127 * x217 * x24,
            x105 * x218 * x24,
            x140 * x219 * x24,
            x102 * x220 * x24,
            x148 * x217 * x22,
            x127 * x218 * x22,
            x105 * x219 * x22,
            x140 * x22 * x220,
            x102 * x22 * x221,
            x222 * x225,
            x119 * x227 * x232,
            x120 * x222 * x232,
            x119 * x236 * x242,
            x139 * x227 * x242,
            x123 * x233 * x242,
            x119 * x246 * x253,
            x139 * x236 * x253,
            x123 * x227 * x253,
            x144 * x233 * x253,
            x119 * x257 * x259,
            x139 * x246 * x259,
            x123 * x236 * x259,
            x144 * x227 * x259,
            x157 * x233 * x259,
            x119 * x223 * x261,
            x119 * x228 * x263,
            x139 * x228 * x261,
            x119 * x237 * x265,
            x139 * x237 * x263,
            x123 * x237 * x261,
            x119 * x247 * x267,
            x139 * x247 * x265,
            x123 * x247 * x263,
            x144 * x247 * x261,
            x119 * x248 * x271,
            x139 * x248 * x267,
            x123 * x248 * x265,
            x144 * x248 * x263,
            x157 * x248 * x261,
            x172 * x222 * x223,
            x173 * x227 * x228,
            x175 * x228 * x233,
            x173 * x236 * x237,
            x175 * x227 * x237,
            x177 * x233 * x237,
            x173 * x246 * x247,
            x175 * x236 * x247,
            x177 * x227 * x247,
            x179 * x233 * x247,
            x173 * x248 * x257,
            x175 * x246 * x248,
            x177 * x236 * x248,
            x179 * x227 * x248,
            x183 * x233 * x248,
            x119 * x274 * x67,
            x119 * x277 * x89,
            x139 * x274 * x89,
            x109 * x119 * x281,
            x109 * x139 * x277,
            x109 * x123 * x274,
            x119 * x135 * x284,
            x135 * x139 * x281,
            x123 * x135 * x277,
            x135 * x144 * x274,
            x119 * x151 * x287,
            x139 * x151 * x284,
            x123 * x151 * x281,
            x144 * x151 * x277,
            x151 * x157 * x274,
            x173 * x261 * x67,
            x173 * x263 * x89,
            x175 * x261 * x89,
            x109 * x173 * x265,
            x109 * x175 * x263,
            x109 * x177 * x261,
            x135 * x173 * x267,
            x135 * x175 * x265,
            x135 * x177 * x263,
            x135 * x179 * x261,
            x151 * x173 * x271,
            x151 * x175 * x267,
            x151 * x177 * x265,
            x151 * x179 * x263,
            x151 * x183 * x261,
            x198 * x233 * x67,
            x198 * x227 * x89,
            x201 * x233 * x89,
            x109 * x198 * x236,
            x109 * x201 * x227,
            x109 * x205 * x233,
            x135 * x198 * x246,
            x135 * x201 * x236,
            x135 * x205 * x227,
            x135 * x208 * x233,
            x151 * x198 * x257,
            x151 * x201 * x246,
            x151 * x205 * x236,
            x151 * x208 * x227,
            x151 * x209 * x233,
            x119 * x289 * x63,
            x119 * x291 * x45,
            x139 * x289 * x45,
            x119 * x17 * x295,
            x139 * x17 * x291,
            x123 * x17 * x289,
            x119 * x15 * x298,
            x139 * x15 * x295,
            x123 * x15 * x291,
            x144 * x15 * x289,
            x10 * x301,
            x101 * x298 * x302,
            x11 * x123 * x295,
            x11 * x144 * x291,
            x11 * x157 * x289,
            x173 * x274 * x63,
            x173 * x277 * x45,
            x175 * x274 * x45,
            x17 * x173 * x281,
            x17 * x175 * x277,
            x17 * x177 * x274,
            x15 * x173 * x284,
            x15 * x175 * x281,
            x15 * x177 * x277,
            x15 * x179 * x274,
            x171 * x287 * x302,
            x11 * x175 * x284,
            x11 * x177 * x281,
            x11 * x179 * x277,
            x11 * x183 * x274,
            x198 * x261 * x63,
            x198 * x263 * x45,
            x201 * x261 * x45,
            x17 * x198 * x265,
            x17 * x201 * x263,
            x17 * x205 * x261,
            x15 * x198 * x267,
            x15 * x201 * x265,
            x15 * x205 * x263,
            x15 * x208 * x261,
            x11 * x198 * x271,
            x11 * x201 * x267,
            x11 * x205 * x265,
            x11 * x208 * x263,
            x11 * x209 * x261,
            x217 * x233 * x63,
            x217 * x227 * x45,
            x218 * x233 * x45,
            x17 * x217 * x236,
            x17 * x218 * x227,
            x17 * x219 * x233,
            x15 * x217 * x246,
            x15 * x218 * x236,
            x15 * x219 * x227,
            x15 * x220 * x233,
            x11 * x217 * x257,
            x11 * x218 * x246,
            x11 * x219 * x236,
            x11 * x220 * x227,
            x221 * x222 * x304,
            x225 * x305,
            x232 * x306 * x85,
            x102 * x232 * x308,
            x105 * x242 * x309,
            x140 * x242 * x308,
            x102 * x242 * x312,
            x127 * x253 * x309,
            x105 * x253 * x308,
            x140 * x253 * x312,
            x102 * x253 * x316,
            x148 * x259 * x309,
            x127 * x259 * x308,
            x105 * x259 * x312,
            x140 * x259 * x316,
            x102 * x259 * x320,
            x158 * x223 * x306,
            x162 * x228 * x309,
            x160 * x228 * x308,
            x164 * x237 * x309,
            x162 * x237 * x308,
            x160 * x237 * x312,
            x166 * x247 * x309,
            x164 * x247 * x308,
            x162 * x247 * x312,
            x160 * x247 * x316,
            x170 * x248 * x309,
            x166 * x248 * x308,
            x164 * x248 * x312,
            x162 * x248 * x316,
            x160 * x248 * x320,
            x102 * x223 * x322,
            x140 * x228 * x322,
            x102 * x228 * x324,
            x105 * x237 * x322,
            x140 * x237 * x324,
            x102 * x237 * x326,
            x127 * x247 * x322,
            x105 * x247 * x324,
            x140 * x247 * x326,
            x102 * x247 * x328,
            x148 * x248 * x322,
            x127 * x248 * x324,
            x105 * x248 * x326,
            x140 * x248 * x328,
            x102 * x248 * x332,
            x185 * x309 * x67,
            x188 * x309 * x89,
            x185 * x308 * x89,
            x109 * x192 * x309,
            x109 * x188 * x308,
            x109 * x185 * x312,
            x135 * x195 * x309,
            x135 * x192 * x308,
            x135 * x188 * x312,
            x135 * x185 * x316,
            x151 * x196 * x309,
            x151 * x195 * x308,
            x151 * x192 * x312,
            x151 * x188 * x316,
            x151 * x185 * x320,
            x160 * x322 * x67,
            x162 * x322 * x89,
            x160 * x324 * x89,
            x109 * x164 * x322,
            x109 * x162 * x324,
            x109 * x160 * x326,
            x135 * x166 * x322,
            x135 * x164 * x324,
            x135 * x162 * x326,
            x135 * x160 * x328,
            x151 * x170 * x322,
            x151 * x166 * x324,
            x151 * x164 * x326,
            x151 * x162 * x328,
            x151 * x160 * x332,
            x102 * x335 * x67,
            x140 * x335 * x89,
            x102 * x338 * x89,
            x105 * x109 * x335,
            x109 * x140 * x338,
            x102 * x109 * x342,
            x127 * x135 * x335,
            x105 * x135 * x338,
            x135 * x140 * x342,
            x102 * x135 * x345,
            x148 * x151 * x335,
            x127 * x151 * x338,
            x105 * x151 * x342,
            x140 * x151 * x345,
            x102 * x151 * x348,
            x211 * x309 * x63,
            x212 * x309 * x45,
            x211 * x308 * x45,
            x17 * x213 * x309,
            x17 * x212 * x308,
            x17 * x211 * x312,
            x15 * x214 * x309,
            x15 * x213 * x308,
            x15 * x212 * x312,
            x15 * x211 * x316,
            x215 * x302 * x305,
            x11 * x214 * x308,
            x11 * x213 * x312,
            x11 * x212 * x316,
            x11 * x211 * x320,
            x185 * x322 * x63,
            x188 * x322 * x45,
            x185 * x324 * x45,
            x17 * x192 * x322,
            x17 * x188 * x324,
            x17 * x185 * x326,
            x15 * x195 * x322,
            x15 * x192 * x324,
            x15 * x188 * x326,
            x15 * x185 * x328,
            x11 * x196 * x322,
            x11 * x195 * x324,
            x11 * x192 * x326,
            x11 * x188 * x328,
            x11 * x185 * x332,
            x160 * x335 * x63,
            x162 * x335 * x45,
            x160 * x338 * x45,
            x164 * x17 * x335,
            x162 * x17 * x338,
            x160 * x17 * x342,
            x15 * x166 * x335,
            x15 * x164 * x338,
            x15 * x162 * x342,
            x15 * x160 * x345,
            x11 * x170 * x335,
            x11 * x166 * x338,
            x11 * x164 * x342,
            x11 * x162 * x345,
            x158 * x304 * x348,
            x102 * x350 * x63,
            x140 * x350 * x45,
            x102 * x352 * x45,
            x105 * x17 * x350,
            x140 * x17 * x352,
            x102 * x17 * x356,
            x127 * x15 * x350,
            x105 * x15 * x352,
            x140 * x15 * x356,
            x102 * x15 * x359,
            x11 * x148 * x350,
            x11 * x127 * x352,
            x105 * x11 * x356,
            x304 * x359 * x85,
            x10 * x361,
            x119 * x363 * x364,
            x119 * x367 * x368,
            x139 * x363 * x368,
            x119 * x371 * x372,
            x139 * x367 * x372,
            x123 * x363 * x372,
            x119 * x377 * x378,
            x139 * x371 * x378,
            x123 * x367 * x378,
            x144 * x363 * x378,
            x119 * x382 * x384,
            x139 * x377 * x384,
            x123 * x371 * x384,
            x144 * x367 * x384,
            x157 * x363 * x384,
            x119 * x224 * x386,
            x119 * x231 * x388,
            x139 * x231 * x386,
            x119 * x241 * x390,
            x139 * x241 * x388,
            x123 * x241 * x386,
            x119 * x252 * x392,
            x139 * x252 * x390,
            x123 * x252 * x388,
            x144 * x252 * x386,
            x119 * x383 * x396,
            x139 * x383 * x392,
            x123 * x383 * x390,
            x144 * x383 * x388,
            x157 * x383 * x386,
            x173 * x224 * x363,
            x173 * x231 * x367,
            x175 * x231 * x363,
            x173 * x241 * x371,
            x175 * x241 * x367,
            x177 * x241 * x363,
            x173 * x252 * x377,
            x175 * x252 * x371,
            x177 * x252 * x367,
            x179 * x252 * x363,
            x173 * x382 * x383,
            x175 * x377 * x383,
            x177 * x371 * x383,
            x179 * x367 * x383,
            x183 * x363 * x383,
            x119 * x397 * x75,
            x119 * x400 * x93,
            x139 * x397 * x93,
            x116 * x119 * x404,
            x116 * x139 * x400,
            x116 * x123 * x397,
            x119 * x250 * x407,
            x139 * x250 * x404,
            x123 * x250 * x400,
            x144 * x250 * x397,
            x2 * x409,
            x101 * x407 * x410,
            x123 * x149 * x404,
            x144 * x149 * x400,
            x149 * x157 * x397,
            x173 * x386 * x75,
            x173 * x388 * x93,
            x175 * x386 * x93,
            x116 * x173 * x390,
            x116 * x175 * x388,
            x116 * x177 * x386,
            x173 * x250 * x392,
            x175 * x250 * x390,
            x177 * x250 * x388,
            x179 * x250 * x386,
            x171 * x396 * x410,
            x149 * x175 * x392,
            x149 * x177 * x390,
            x149 * x179 * x388,
            x149 * x183 * x386,
            x198 * x363 * x75,
            x198 * x367 * x93,
            x201 * x363 * x93,
            x116 * x198 * x371,
            x116 * x201 * x367,
            x116 * x205 * x363,
            x198 * x250 * x377,
            x201 * x250 * x371,
            x205 * x250 * x367,
            x208 * x250 * x363,
            x149 * x198 * x382,
            x149 * x201 * x377,
            x149 * x205 * x371,
            x149 * x208 * x367,
            x149 * x209 * x363,
            x119 * x411 * x61,
            x119 * x412 * x51,
            x139 * x411 * x51,
            x119 * x413 * x49,
            x139 * x412 * x49,
            x123 * x411 * x49,
            x4 * x414,
            x101 * x413 * x415,
            x123 * x412 * x9,
            x144 * x411 * x9,
            x300
            * (
                x158 * x408
                + x3 * (2 * x285 + 2 * x286 + 2 * x393 + 2 * x395 + 4 * x405 + 4 * x406)
            ),
            x101 * x414,
            x123 * x413 * x8,
            x144 * x412 * x8,
            x157 * x411 * x8,
            x173 * x397 * x61,
            x173 * x400 * x51,
            x175 * x397 * x51,
            x173 * x404 * x49,
            x175 * x400 * x49,
            x177 * x397 * x49,
            x171 * x407 * x415,
            x175 * x404 * x9,
            x177 * x400 * x9,
            x179 * x397 * x9,
            x171 * x409,
            x175 * x407 * x8,
            x177 * x404 * x8,
            x179 * x400 * x8,
            x183 * x397 * x8,
            x198 * x386 * x61,
            x198 * x388 * x51,
            x201 * x386 * x51,
            x198 * x390 * x49,
            x201 * x388 * x49,
            x205 * x386 * x49,
            x198 * x392 * x9,
            x201 * x390 * x9,
            x205 * x388 * x9,
            x208 * x386 * x9,
            x198 * x396 * x8,
            x201 * x392 * x8,
            x205 * x390 * x8,
            x208 * x388 * x8,
            x209 * x386 * x8,
            x217 * x363 * x61,
            x217 * x367 * x51,
            x218 * x363 * x51,
            x217 * x371 * x49,
            x218 * x367 * x49,
            x219 * x363 * x49,
            x217 * x377 * x9,
            x218 * x371 * x9,
            x219 * x367 * x9,
            x220 * x363 * x9,
            x217 * x382 * x8,
            x218 * x377 * x8,
            x219 * x371 * x8,
            x220 * x367 * x8,
            x221 * x363 * x8,
            x222 * x306 * x364,
            x227 * x309 * x368,
            x233 * x308 * x368,
            x236 * x309 * x372,
            x227 * x308 * x372,
            x233 * x312 * x372,
            x246 * x309 * x378,
            x236 * x308 * x378,
            x227 * x312 * x378,
            x233 * x316 * x378,
            x257 * x309 * x384,
            x246 * x308 * x384,
            x236 * x312 * x384,
            x227 * x316 * x384,
            x233 * x320 * x384,
            x224 * x261 * x309,
            x231 * x263 * x309,
            x231 * x261 * x308,
            x241 * x265 * x309,
            x241 * x263 * x308,
            x241 * x261 * x312,
            x252 * x267 * x309,
            x252 * x265 * x308,
            x252 * x263 * x312,
            x252 * x261 * x316,
            x271 * x309 * x383,
            x267 * x308 * x383,
            x265 * x312 * x383,
            x263 * x316 * x383,
            x261 * x320 * x383,
            x224 * x233 * x322,
            x227 * x231 * x322,
            x231 * x233 * x324,
            x236 * x241 * x322,
            x227 * x241 * x324,
            x233 * x241 * x326,
            x246 * x252 * x322,
            x236 * x252 * x324,
            x227 * x252 * x326,
            x233 * x252 * x328,
            x257 * x322 * x383,
            x246 * x324 * x383,
            x236 * x326 * x383,
            x227 * x328 * x383,
            x233 * x332 * x383,
            x274 * x309 * x75,
            x277 * x309 * x93,
            x274 * x308 * x93,
            x116 * x281 * x309,
            x116 * x277 * x308,
            x116 * x274 * x312,
            x250 * x284 * x309,
            x250 * x281 * x308,
            x250 * x277 * x312,
            x250 * x274 * x316,
            x287 * x305 * x410,
            x149 * x284 * x308,
            x149 * x281 * x312,
            x149 * x277 * x316,
            x149 * x274 * x320,
            x261 * x322 * x75,
            x263 * x322 * x93,
            x261 * x324 * x93,
            x116 * x265 * x322,
            x116 * x263 * x324,
            x116 * x261 * x326,
            x250 * x267 * x322,
            x250 * x265 * x324,
            x250 * x263 * x326,
            x250 * x261 * x328,
            x149 * x271 * x322,
            x149 * x267 * x324,
            x149 * x265 * x326,
            x149 * x263 * x328,
            x149 * x261 * x332,
            x233 * x335 * x75,
            x227 * x335 * x93,
            x233 * x338 * x93,
            x116 * x236 * x335,
            x116 * x227 * x338,
            x116 * x233 * x342,
            x246 * x250 * x335,
            x236 * x250 * x338,
            x227 * x250 * x342,
            x233 * x250 * x345,
            x149 * x257 * x335,
            x149 * x246 * x338,
            x149 * x236 * x342,
            x149 * x227 * x345,
            x2 * x348 * x416,
            x289 * x309 * x61,
            x291 * x309 * x51,
            x289 * x308 * x51,
            x295 * x309 * x49,
            x291 * x308 * x49,
            x289 * x312 * x49,
            x298 * x305 * x415,
            x295 * x308 * x9,
            x291 * x312 * x9,
            x289 * x316 * x9,
            x301 * x305,
            x298 * x308 * x8,
            x295 * x312 * x8,
            x291 * x316 * x8,
            x289 * x320 * x8,
            x274 * x322 * x61,
            x277 * x322 * x51,
            x274 * x324 * x51,
            x281 * x322 * x49,
            x277 * x324 * x49,
            x274 * x326 * x49,
            x284 * x322 * x9,
            x281 * x324 * x9,
            x277 * x326 * x9,
            x274 * x328 * x9,
            x287 * x322 * x8,
            x284 * x324 * x8,
            x281 * x326 * x8,
            x277 * x328 * x8,
            x274 * x332 * x8,
            x261 * x335 * x61,
            x263 * x335 * x51,
            x261 * x338 * x51,
            x265 * x335 * x49,
            x263 * x338 * x49,
            x261 * x342 * x49,
            x267 * x335 * x9,
            x265 * x338 * x9,
            x263 * x342 * x9,
            x261 * x345 * x9,
            x271 * x335 * x8,
            x267 * x338 * x8,
            x265 * x342 * x8,
            x263 * x345 * x8,
            x261 * x348 * x8,
            x233 * x350 * x61,
            x227 * x350 * x51,
            x233 * x352 * x51,
            x236 * x350 * x49,
            x227 * x352 * x49,
            x233 * x356 * x49,
            x246 * x350 * x9,
            x236 * x352 * x9,
            x227 * x356 * x9,
            x359 * x4 * x416,
            x257 * x350 * x8,
            x246 * x352 * x8,
            x236 * x356 * x8,
            x227 * x359 * x8,
            x222 * x361,
            x102 * x364 * x418,
            x140 * x368 * x418,
            x102 * x368 * x421,
            x105 * x372 * x418,
            x140 * x372 * x421,
            x102 * x372 * x424,
            x127 * x378 * x418,
            x105 * x378 * x421,
            x140 * x378 * x424,
            x102 * x378 * x429,
            x148 * x384 * x418,
            x127 * x384 * x421,
            x105 * x384 * x424,
            x140 * x384 * x429,
            x102 * x384 * x433,
            x160 * x224 * x418,
            x162 * x231 * x418,
            x160 * x231 * x421,
            x164 * x241 * x418,
            x162 * x241 * x421,
            x160 * x241 * x424,
            x166 * x252 * x418,
            x164 * x252 * x421,
            x162 * x252 * x424,
            x160 * x252 * x429,
            x170 * x383 * x418,
            x166 * x383 * x421,
            x164 * x383 * x424,
            x162 * x383 * x429,
            x160 * x383 * x433,
            x102 * x224 * x435,
            x140 * x231 * x435,
            x102 * x231 * x437,
            x105 * x241 * x435,
            x140 * x241 * x437,
            x102 * x241 * x439,
            x127 * x252 * x435,
            x105 * x252 * x437,
            x140 * x252 * x439,
            x102 * x252 * x441,
            x148 * x383 * x435,
            x127 * x383 * x437,
            x105 * x383 * x439,
            x140 * x383 * x441,
            x102 * x383 * x445,
            x185 * x418 * x75,
            x188 * x418 * x93,
            x185 * x421 * x93,
            x116 * x192 * x418,
            x116 * x188 * x421,
            x116 * x185 * x424,
            x195 * x250 * x418,
            x192 * x250 * x421,
            x188 * x250 * x424,
            x185 * x250 * x429,
            x149 * x196 * x418,
            x149 * x195 * x421,
            x149 * x192 * x424,
            x149 * x188 * x429,
            x149 * x185 * x433,
            x160 * x435 * x75,
            x162 * x435 * x93,
            x160 * x437 * x93,
            x116 * x164 * x435,
            x116 * x162 * x437,
            x116 * x160 * x439,
            x166 * x250 * x435,
            x164 * x250 * x437,
            x162 * x250 * x439,
            x160 * x250 * x441,
            x149 * x170 * x435,
            x149 * x166 * x437,
            x149 * x164 * x439,
            x149 * x162 * x441,
            x158 * x445 * x446,
            x102 * x447 * x75,
            x140 * x447 * x93,
            x102 * x450 * x93,
            x105 * x116 * x447,
            x116 * x140 * x450,
            x102 * x116 * x454,
            x127 * x250 * x447,
            x105 * x250 * x450,
            x140 * x250 * x454,
            x102 * x250 * x457,
            x148 * x149 * x447,
            x127 * x149 * x450,
            x105 * x149 * x454,
            x446 * x457 * x85,
            x2 * x459,
            x211 * x418 * x61,
            x212 * x418 * x51,
            x211 * x421 * x51,
            x213 * x418 * x49,
            x212 * x421 * x49,
            x211 * x424 * x49,
            x214 * x418 * x9,
            x213 * x421 * x9,
            x212 * x424 * x9,
            x211 * x429 * x9,
            x215 * x418 * x8,
            x214 * x421 * x8,
            x213 * x424 * x8,
            x212 * x429 * x8,
            x211 * x433 * x8,
            x185 * x435 * x61,
            x188 * x435 * x51,
            x185 * x437 * x51,
            x192 * x435 * x49,
            x188 * x437 * x49,
            x185 * x439 * x49,
            x195 * x435 * x9,
            x192 * x437 * x9,
            x188 * x439 * x9,
            x185 * x441 * x9,
            x196 * x435 * x8,
            x195 * x437 * x8,
            x192 * x439 * x8,
            x188 * x441 * x8,
            x185 * x445 * x8,
            x160 * x447 * x61,
            x162 * x447 * x51,
            x160 * x450 * x51,
            x164 * x447 * x49,
            x162 * x450 * x49,
            x160 * x454 * x49,
            x166 * x447 * x9,
            x164 * x450 * x9,
            x162 * x454 * x9,
            x158 * x457 * x460,
            x170 * x447 * x8,
            x166 * x450 * x8,
            x164 * x454 * x8,
            x162 * x457 * x8,
            x158 * x459,
            x102 * x461 * x61,
            x140 * x461 * x51,
            x102 * x462 * x51,
            x105 * x461 * x49,
            x140 * x462 * x49,
            x102 * x463 * x49,
            x127 * x461 * x9,
            x105 * x462 * x9,
            x460 * x463 * x85,
            x4 * x464,
            x148 * x461 * x8,
            x127 * x462 * x8,
            x105 * x463 * x8,
            x464 * x85,
            x303
            * (
                x171 * x458
                + x3 * (2 * x346 + 2 * x347 + 2 * x442 + 2 * x444 + 4 * x455 + 4 * x456)
            ),
        ]
    )


def quadrupole3d_40(a, A, b, B, C):
    """Cartesian 3D (gs) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - A[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3 ** 2 * x7
    x9 = x0 * x7
    x10 = 3 * x9
    x11 = -x2 - C[0]
    x12 = x3 * x7
    x13 = x11 * x12
    x14 = x10 + 2 * x13
    x15 = x0 * (x14 + x8)
    x16 = x11 ** 2 * x7
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
    x36 = -x35 - A[1]
    x37 = x30 * x34
    x38 = -x1 * (a * A[2] + b * B[2])
    x39 = -x38 - A[2]
    x40 = x31 * x6
    x41 = x0 * x40
    x42 = x36 ** 2 * x40
    x43 = x41 + x42
    x44 = x32 * x6
    x45 = x34 * x39
    x46 = x0 * x44
    x47 = x39 ** 2 * x44
    x48 = x46 + x47
    x49 = 2 * x41
    x50 = x36 * x43 + x36 * x49
    x51 = x39 * x44
    x52 = x36 * x40
    x53 = 2 * x46
    x54 = x39 * x48 + x39 * x53
    x55 = 3 * x41
    x56 = x0 * (3 * x42 + x55) + x36 * x50
    x57 = 3 * x46
    x58 = x0 * (3 * x47 + x57) + x39 * x54
    x59 = -x35 - C[1]
    x60 = x8 + x9
    x61 = 2 * x0 * x12 + x3 * x60
    x62 = x15 + x28
    x63 = x34 * (x0 * (3 * x24 + 3 * x26 + x61) + x3 * x62)
    x64 = x52 * x59
    x65 = x41 + x64
    x66 = x40 * x59
    x67 = x0 * (x52 + x66)
    x68 = x36 * x65
    x69 = x67 + x68
    x70 = x55 + 2 * x64
    x71 = x0 * (x42 + x70)
    x72 = x36 * x69
    x73 = x71 + x72
    x74 = x33 * x5
    x75 = x74 * (x0 * (x50 + 3 * x67 + 3 * x68) + x36 * x73)
    x76 = x11 * x74
    x77 = numpy.pi * x1 * x31 * x5
    x78 = x11 * x77
    x79 = -x38 - C[2]
    x80 = x34 * x79
    x81 = x51 * x79
    x82 = x46 + x81
    x83 = x44 * x79
    x84 = x0 * (x51 + x83)
    x85 = x39 * x82
    x86 = x84 + x85
    x87 = x57 + 2 * x81
    x88 = x0 * (x47 + x87)
    x89 = x39 * x86
    x90 = x88 + x89
    x91 = x77 * (x0 * (x54 + 3 * x84 + 3 * x85) + x39 * x90)
    x92 = x40 * x59 ** 2
    x93 = x41 + x92
    x94 = x0 * (x10 + 3 * x8) + x3 * x61
    x95 = x36 * x93
    x96 = x49 * x59 + x95
    x97 = x0 * (x70 + x92)
    x98 = x36 * x96
    x99 = x97 + x98
    x100 = x0 * (4 * x41 * x59 + 2 * x67 + 2 * x68 + 2 * x95) + x36 * x99
    x101 = x100 * x74
    x102 = x3 * x74
    x103 = x3 * x77
    x104 = x44 * x79 ** 2
    x105 = x104 + x46
    x106 = x105 * x39
    x107 = x106 + x53 * x79
    x108 = x0 * (x104 + x87)
    x109 = x107 * x39
    x110 = x108 + x109
    x111 = x0 * (2 * x106 + 4 * x46 * x79 + 2 * x84 + 2 * x85) + x110 * x39
    x112 = x111 * x77

    # 90 item(s)
    return numpy.array(
        [
            x34 * (x0 * (2 * x15 + 3 * x17 + 3 * x23 + 2 * x28) + x3 * x30),
            x36 * x37,
            x37 * x39,
            x29 * x43 * x44,
            x29 * x36 * x45,
            x29 * x40 * x48,
            x22 * x44 * x50,
            x22 * x43 * x51,
            x22 * x48 * x52,
            x22 * x40 * x54,
            x20 * x44 * x56,
            x20 * x50 * x51,
            x20 * x43 * x48,
            x20 * x52 * x54,
            x20 * x40 * x58,
            x59 * x63,
            x44 * x62 * x65,
            x45 * x59 * x62,
            x27 * x44 * x69,
            x27 * x51 * x65,
            x27 * x48 * x66,
            x25 * x44 * x73,
            x25 * x51 * x69,
            x25 * x48 * x65,
            x25 * x54 * x66,
            x11 * x75,
            x39 * x73 * x76,
            x18 * x48 * x69,
            x18 * x54 * x65,
            x58 * x59 * x78,
            x63 * x79,
            x36 * x62 * x80,
            x40 * x62 * x82,
            x27 * x43 * x83,
            x27 * x52 * x82,
            x27 * x40 * x86,
            x25 * x50 * x83,
            x25 * x43 * x82,
            x25 * x52 * x86,
            x25 * x40 * x90,
            x56 * x76 * x79,
            x18 * x50 * x82,
            x18 * x43 * x86,
            x36 * x78 * x90,
            x11 * x91,
            x44 * x93 * x94,
            x44 * x61 * x96,
            x51 * x61 * x93,
            x44 * x60 * x99,
            x51 * x60 * x96,
            x48 * x60 * x93,
            x101 * x3,
            x102 * x39 * x99,
            x12 * x48 * x96,
            x12 * x54 * x93,
            x74 * (x0 * (2 * x71 + 2 * x72 + 3 * x97 + 3 * x98) + x100 * x36),
            x101 * x39,
            x48 * x7 * x99,
            x54 * x7 * x96,
            x58 * x7 * x93,
            x59 * x80 * x94,
            x61 * x65 * x83,
            x61 * x66 * x82,
            x60 * x69 * x83,
            x60 * x65 * x82,
            x60 * x66 * x86,
            x102 * x73 * x79,
            x12 * x69 * x82,
            x12 * x65 * x86,
            x103 * x59 * x90,
            x75 * x79,
            x7 * x73 * x82,
            x69 * x7 * x86,
            x65 * x7 * x90,
            x59 * x91,
            x105 * x40 * x94,
            x105 * x52 * x61,
            x107 * x40 * x61,
            x105 * x43 * x60,
            x107 * x52 * x60,
            x110 * x40 * x60,
            x105 * x12 * x50,
            x107 * x12 * x43,
            x103 * x110 * x36,
            x112 * x3,
            x105 * x56 * x7,
            x107 * x50 * x7,
            x110 * x43 * x7,
            x112 * x36,
            x77 * (x0 * (3 * x108 + 3 * x109 + 2 * x88 + 2 * x89) + x111 * x39),
        ]
    )


def quadrupole3d_41(a, A, b, B, C):
    """Cartesian 3D (gp) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - B[0]
    x5 = -x1 - C[0]
    x6 = a * b * x0
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x9 = x7 * x8
    x10 = x5 * x9
    x11 = x10 * x4
    x12 = x5 ** 2 * x9
    x13 = x3 * x9
    x14 = 3 * x13
    x15 = x12 + x14
    x16 = x3 * (2 * x11 + x15)
    x17 = x10 * x3
    x18 = 2 * x17
    x19 = x12 + x13
    x20 = x19 * x4
    x21 = x18 + x20
    x22 = x2 * x21
    x23 = x16 + x22
    x24 = x2 * x23
    x25 = x19 * x2
    x26 = 4 * x17
    x27 = x4 * x9
    x28 = x3 * (x10 + x27)
    x29 = x11 + x13
    x30 = x2 * x29
    x31 = 2 * x28 + 2 * x30
    x32 = x3 * (x20 + x25 + x26 + x31)
    x33 = x24 + x32
    x34 = x2 * x7
    x35 = x34 * x8
    x36 = x35 * x4
    x37 = x35 * x5
    x38 = x3 * (x11 + x14 + x36 + x37)
    x39 = x28 + x30
    x40 = x2 * x39
    x41 = 2 * x37
    x42 = x3 * (x15 + x41)
    x43 = x18 + x25
    x44 = x2 * x43
    x45 = x42 + x44
    x46 = x2 * x33 + x3 * (2 * x16 + 2 * x22 + 2 * x38 + 2 * x40 + x45)
    x47 = x3 * (x10 + x35)
    x48 = x13 + x37
    x49 = x2 * x48
    x50 = x47 + x49
    x51 = x3 * (x27 + x35)
    x52 = x13 + x36
    x53 = x2 * x52
    x54 = x51 + x53
    x55 = x3 * (x31 + x50 + x54)
    x56 = x38 + x40
    x57 = x2 * x56
    x58 = x2 * x45 + x3 * (2 * x25 + x26 + 2 * x47 + 2 * x49)
    x59 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x60 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x61 = numpy.pi * x0 * x60
    x62 = x59 * x61
    x63 = -x0 * (a * A[1] + b * B[1])
    x64 = -x63 - B[1]
    x65 = x2 ** 2 * x9
    x66 = x14 + x65
    x67 = x3 * (x41 + x66)
    x68 = x2 * x50
    x69 = x62 * (x2 * x58 + x3 * (3 * x42 + 3 * x44 + 2 * x67 + 2 * x68))
    x70 = -x0 * (a * A[2] + b * B[2])
    x71 = -x70 - B[2]
    x72 = -x63 - A[1]
    x73 = x46 * x62
    x74 = x3 * x8
    x75 = x59 * x74
    x76 = x59 * x8
    x77 = x72 * x76
    x78 = x64 * x77
    x79 = x75 + x78
    x80 = x60 * x8
    x81 = x58 * x62
    x82 = -x70 - A[2]
    x83 = x60 * x74
    x84 = x80 * x82
    x85 = x71 * x84
    x86 = x83 + x85
    x87 = x72 ** 2 * x76
    x88 = x75 + x87
    x89 = x64 * x76
    x90 = x3 * (x77 + x89)
    x91 = x72 * x79
    x92 = x90 + x91
    x93 = x71 * x80
    x94 = x62 * x82
    x95 = x80 * x82 ** 2
    x96 = x83 + x95
    x97 = x3 * (x84 + x93)
    x98 = x82 * x86
    x99 = x97 + x98
    x100 = 2 * x75
    x101 = x100 * x72 + x72 * x88
    x102 = 3 * x75
    x103 = x102 + x87
    x104 = x3 * (x103 + 2 * x78) + x72 * x92
    x105 = 2 * x83
    x106 = x105 * x82 + x82 * x96
    x107 = 3 * x83
    x108 = x107 + x95
    x109 = x3 * (x108 + 2 * x85) + x82 * x99
    x110 = x101 * x72 + x3 * (x102 + 3 * x87)
    x111 = x104 * x72 + x3 * (x101 + 3 * x90 + 3 * x91)
    x112 = x106 * x82 + x3 * (x107 + 3 * x95)
    x113 = x109 * x82 + x3 * (x106 + 3 * x97 + 3 * x98)
    x114 = -x63 - C[1]
    x115 = x55 + x57
    x116 = x67 + x68
    x117 = x2 * x54 + x3 * (2 * x36 + x66)
    x118 = x62 * (x115 * x2 + x3 * (x116 + x117 + 3 * x38 + 3 * x40))
    x119 = x114 * x89
    x120 = x119 + x75
    x121 = x13 + x65
    x122 = x121 * x2 + 2 * x3 * x35
    x123 = x116 * x2 + x3 * (x122 + 3 * x47 + 3 * x49)
    x124 = x123 * x62
    x125 = x114 * x77
    x126 = x125 + x75
    x127 = x114 * x76
    x128 = x3 * (x127 + x89)
    x129 = x120 * x72
    x130 = x128 + x129
    x131 = x3 * (x127 + x77)
    x132 = x126 * x72
    x133 = x131 + x132
    x134 = x3 * (x102 + x119 + x125 + x78)
    x135 = x130 * x72
    x136 = x134 + x135
    x137 = 2 * x125
    x138 = x3 * (x103 + x137)
    x139 = x133 * x72
    x140 = x138 + x139
    x141 = x136 * x72
    x142 = 2 * x128 + 2 * x129
    x143 = x3 * (x133 + x142 + x92)
    x144 = x141 + x143
    x145 = x140 * x72 + x3 * (x101 + 3 * x131 + 3 * x132)
    x146 = x61 * x7
    x147 = x146 * (x144 * x72 + x3 * (x104 + 3 * x134 + 3 * x135 + x140))
    x148 = x146 * x5
    x149 = x146 * x82
    x150 = numpy.pi * x0 * x59
    x151 = x150 * x7
    x152 = x151 * x5
    x153 = -x70 - C[2]
    x154 = x153 * x93
    x155 = x154 + x83
    x156 = x153 * x62
    x157 = x153 * x80
    x158 = x153 * x84
    x159 = x158 + x83
    x160 = x3 * (x157 + x93)
    x161 = x155 * x82
    x162 = x160 + x161
    x163 = x3 * (x157 + x84)
    x164 = x159 * x82
    x165 = x163 + x164
    x166 = x3 * (x107 + x154 + x158 + x85)
    x167 = x162 * x82
    x168 = x166 + x167
    x169 = 2 * x158
    x170 = x3 * (x108 + x169)
    x171 = x165 * x82
    x172 = x170 + x171
    x173 = x168 * x82
    x174 = 2 * x160 + 2 * x161
    x175 = x3 * (x165 + x174 + x99)
    x176 = x173 + x175
    x177 = x151 * x72
    x178 = x172 * x82 + x3 * (x106 + 3 * x163 + 3 * x164)
    x179 = x151 * (x176 * x82 + x3 * (x109 + 3 * x166 + 3 * x167 + x172))
    x180 = x114 ** 2 * x76
    x181 = x180 + x75
    x182 = x117 * x2 + x3 * (x122 + 3 * x51 + 3 * x53)
    x183 = x100 * x114
    x184 = x181 * x64
    x185 = x183 + x184
    x186 = x122 * x2 + x3 * (x14 + 3 * x65)
    x187 = x181 * x72
    x188 = x183 + x187
    x189 = x102 + x180
    x190 = x3 * (2 * x119 + x189)
    x191 = x185 * x72
    x192 = x190 + x191
    x193 = x3 * (x137 + x189)
    x194 = x188 * x72
    x195 = x193 + x194
    x196 = x192 * x72
    x197 = 4 * x114 * x75
    x198 = x3 * (x142 + x184 + x187 + x197)
    x199 = x196 + x198
    x200 = x195 * x72 + x3 * (2 * x131 + 2 * x132 + 2 * x187 + x197)
    x201 = x199 * x72 + x3 * (2 * x134 + 2 * x135 + 2 * x190 + 2 * x191 + x195)
    x202 = x34 * x61
    x203 = x146 * (x200 * x72 + x3 * (2 * x138 + 2 * x139 + 3 * x193 + 3 * x194))
    x204 = x150 * x34
    x205 = x153 ** 2 * x80
    x206 = x205 + x83
    x207 = x105 * x153
    x208 = x206 * x71
    x209 = x207 + x208
    x210 = x206 * x82
    x211 = x207 + x210
    x212 = x107 + x205
    x213 = x3 * (2 * x154 + x212)
    x214 = x209 * x82
    x215 = x213 + x214
    x216 = x3 * (x169 + x212)
    x217 = x211 * x82
    x218 = x216 + x217
    x219 = x215 * x82
    x220 = 4 * x153 * x83
    x221 = x3 * (x174 + x208 + x210 + x220)
    x222 = x219 + x221
    x223 = x218 * x82 + x3 * (2 * x163 + 2 * x164 + 2 * x210 + x220)
    x224 = x222 * x82 + x3 * (2 * x166 + 2 * x167 + 2 * x213 + 2 * x214 + x218)
    x225 = x151 * (x223 * x82 + x3 * (2 * x170 + 2 * x171 + 3 * x216 + 3 * x217))

    # 270 item(s)
    return numpy.array(
        [
            x62 * (x2 * x46 + x3 * (3 * x24 + 3 * x32 + 2 * x55 + 2 * x57 + x58)),
            x64 * x69,
            x69 * x71,
            x72 * x73,
            x58 * x79 * x80,
            x71 * x72 * x81,
            x73 * x82,
            x64 * x81 * x82,
            x58 * x76 * x86,
            x33 * x80 * x88,
            x45 * x80 * x92,
            x45 * x88 * x93,
            x33 * x72 * x94,
            x45 * x79 * x84,
            x45 * x77 * x86,
            x33 * x76 * x96,
            x45 * x89 * x96,
            x45 * x76 * x99,
            x101 * x23 * x80,
            x104 * x43 * x80,
            x101 * x43 * x93,
            x23 * x84 * x88,
            x43 * x84 * x92,
            x43 * x86 * x88,
            x23 * x77 * x96,
            x43 * x79 * x96,
            x43 * x77 * x99,
            x106 * x23 * x76,
            x106 * x43 * x89,
            x109 * x43 * x76,
            x110 * x21 * x80,
            x111 * x19 * x80,
            x110 * x19 * x93,
            x101 * x21 * x84,
            x104 * x19 * x84,
            x101 * x19 * x86,
            x21 * x88 * x96,
            x19 * x92 * x96,
            x19 * x88 * x99,
            x106 * x21 * x77,
            x106 * x19 * x79,
            x109 * x19 * x77,
            x112 * x21 * x76,
            x112 * x19 * x89,
            x113 * x19 * x76,
            x114 * x118,
            x120 * x123 * x80,
            x114 * x124 * x71,
            x115 * x126 * x80,
            x116 * x130 * x80,
            x116 * x126 * x93,
            x114 * x115 * x94,
            x116 * x120 * x84,
            x116 * x127 * x86,
            x133 * x56 * x80,
            x136 * x50 * x80,
            x133 * x50 * x93,
            x126 * x56 * x84,
            x130 * x50 * x84,
            x126 * x50 * x86,
            x127 * x56 * x96,
            x120 * x50 * x96,
            x127 * x50 * x99,
            x140 * x39 * x80,
            x144 * x48 * x80,
            x140 * x48 * x93,
            x133 * x39 * x84,
            x136 * x48 * x84,
            x133 * x48 * x86,
            x126 * x39 * x96,
            x130 * x48 * x96,
            x126 * x48 * x99,
            x106 * x127 * x39,
            x106 * x120 * x48,
            x109 * x127 * x48,
            x145 * x29 * x80,
            x147 * x5,
            x145 * x148 * x71,
            x140 * x29 * x84,
            x144 * x149 * x5,
            x10 * x140 * x86,
            x133 * x29 * x96,
            x10 * x136 * x96,
            x10 * x133 * x99,
            x106 * x126 * x29,
            x10 * x106 * x130,
            x10 * x109 * x126,
            x112 * x127 * x29,
            x10 * x112 * x120,
            x113 * x114 * x152,
            x118 * x153,
            x124 * x153 * x64,
            x123 * x155 * x76,
            x115 * x156 * x72,
            x116 * x157 * x79,
            x116 * x155 * x77,
            x115 * x159 * x76,
            x116 * x159 * x89,
            x116 * x162 * x76,
            x157 * x56 * x88,
            x157 * x50 * x92,
            x155 * x50 * x88,
            x159 * x56 * x77,
            x159 * x50 * x79,
            x162 * x50 * x77,
            x165 * x56 * x76,
            x165 * x50 * x89,
            x168 * x50 * x76,
            x101 * x157 * x39,
            x104 * x157 * x48,
            x101 * x155 * x48,
            x159 * x39 * x88,
            x159 * x48 * x92,
            x162 * x48 * x88,
            x165 * x39 * x77,
            x165 * x48 * x79,
            x168 * x48 * x77,
            x172 * x39 * x76,
            x172 * x48 * x89,
            x176 * x48 * x76,
            x110 * x157 * x29,
            x111 * x148 * x153,
            x10 * x110 * x155,
            x101 * x159 * x29,
            x10 * x104 * x159,
            x10 * x101 * x162,
            x165 * x29 * x88,
            x10 * x165 * x92,
            x10 * x168 * x88,
            x172 * x29 * x77,
            x10 * x172 * x79,
            x176 * x177 * x5,
            x178 * x29 * x76,
            x152 * x178 * x64,
            x179 * x5,
            x181 * x182 * x80,
            x185 * x186 * x80,
            x181 * x186 * x93,
            x117 * x188 * x80,
            x122 * x192 * x80,
            x122 * x188 * x93,
            x117 * x181 * x84,
            x122 * x185 * x84,
            x122 * x181 * x86,
            x195 * x54 * x80,
            x121 * x199 * x80,
            x121 * x195 * x93,
            x188 * x54 * x84,
            x121 * x192 * x84,
            x121 * x188 * x86,
            x181 * x54 * x96,
            x121 * x185 * x96,
            x121 * x181 * x99,
            x200 * x52 * x80,
            x201 * x202,
            x200 * x202 * x71,
            x195 * x52 * x84,
            x199 * x202 * x82,
            x195 * x35 * x86,
            x188 * x52 * x96,
            x192 * x35 * x96,
            x188 * x35 * x99,
            x106 * x181 * x52,
            x106 * x185 * x35,
            x109 * x181 * x35,
            x203 * x4,
            x146 * (x201 * x72 + x3 * (2 * x141 + 2 * x143 + 3 * x196 + 3 * x198 + x200)),
            x203 * x71,
            x149 * x200 * x4,
            x149 * x201,
            x200 * x86 * x9,
            x195 * x27 * x96,
            x199 * x9 * x96,
            x195 * x9 * x99,
            x106 * x188 * x27,
            x106 * x192 * x9,
            x109 * x188 * x9,
            x112 * x181 * x27,
            x112 * x185 * x9,
            x113 * x181 * x9,
            x114 * x156 * x182,
            x120 * x157 * x186,
            x127 * x155 * x186,
            x117 * x126 * x157,
            x122 * x130 * x157,
            x122 * x126 * x155,
            x117 * x127 * x159,
            x120 * x122 * x159,
            x122 * x127 * x162,
            x133 * x157 * x54,
            x121 * x136 * x157,
            x121 * x133 * x155,
            x126 * x159 * x54,
            x121 * x130 * x159,
            x121 * x126 * x162,
            x127 * x165 * x54,
            x120 * x121 * x165,
            x121 * x127 * x168,
            x140 * x157 * x52,
            x144 * x153 * x202,
            x140 * x155 * x35,
            x133 * x159 * x52,
            x136 * x159 * x35,
            x133 * x162 * x35,
            x126 * x165 * x52,
            x130 * x165 * x35,
            x126 * x168 * x35,
            x127 * x172 * x52,
            x120 * x172 * x35,
            x114 * x176 * x204,
            x145 * x146 * x153 * x4,
            x147 * x153,
            x145 * x155 * x9,
            x140 * x159 * x27,
            x144 * x159 * x9,
            x140 * x162 * x9,
            x133 * x165 * x27,
            x136 * x165 * x9,
            x133 * x168 * x9,
            x126 * x172 * x27,
            x130 * x172 * x9,
            x126 * x176 * x9,
            x114 * x151 * x178 * x4,
            x120 * x178 * x9,
            x114 * x179,
            x182 * x206 * x76,
            x186 * x206 * x89,
            x186 * x209 * x76,
            x117 * x206 * x77,
            x122 * x206 * x79,
            x122 * x209 * x77,
            x117 * x211 * x76,
            x122 * x211 * x89,
            x122 * x215 * x76,
            x206 * x54 * x88,
            x121 * x206 * x92,
            x121 * x209 * x88,
            x211 * x54 * x77,
            x121 * x211 * x79,
            x121 * x215 * x77,
            x218 * x54 * x76,
            x121 * x218 * x89,
            x121 * x222 * x76,
            x101 * x206 * x52,
            x104 * x206 * x35,
            x101 * x209 * x35,
            x211 * x52 * x88,
            x211 * x35 * x92,
            x215 * x35 * x88,
            x218 * x52 * x77,
            x218 * x35 * x79,
            x204 * x222 * x72,
            x223 * x52 * x76,
            x204 * x223 * x64,
            x204 * x224,
            x110 * x206 * x27,
            x111 * x206 * x9,
            x110 * x209 * x9,
            x101 * x211 * x27,
            x104 * x211 * x9,
            x101 * x215 * x9,
            x218 * x27 * x88,
            x218 * x9 * x92,
            x222 * x88 * x9,
            x177 * x223 * x4,
            x223 * x79 * x9,
            x177 * x224,
            x225 * x4,
            x225 * x64,
            x151 * (x224 * x82 + x3 * (2 * x173 + 2 * x175 + 3 * x219 + 3 * x221 + x223)),
        ]
    )


def quadrupole3d_42(a, A, b, B, C):
    """Cartesian 3D (gd) quadrupole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - A[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = -x1 - C[0]
    x5 = -x1 - B[0]
    x6 = a * b * x0
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x9 = x7 * x8
    x10 = x5 * x9
    x11 = x10 * x4
    x12 = 2 * x11
    x13 = x4 ** 2 * x9
    x14 = x3 * x9
    x15 = 3 * x14
    x16 = x13 + x15
    x17 = x3 * (x12 + x16)
    x18 = 2 * x14
    x19 = x18 * x4
    x20 = x13 + x14
    x21 = x20 * x5
    x22 = x19 + x21
    x23 = x22 * x5
    x24 = x17 + x23
    x25 = x2 * x24
    x26 = x11 + x14
    x27 = x26 * x5
    x28 = x4 * x9
    x29 = x3 * (x10 + x28)
    x30 = 2 * x29
    x31 = 4 * x14
    x32 = x31 * x4
    x33 = x30 + x32
    x34 = x3 * (2 * x21 + 2 * x27 + x33)
    x35 = x25 + x34
    x36 = x2 * x35
    x37 = x2 * x22
    x38 = 2 * x37
    x39 = x5 ** 2 * x9
    x40 = x15 + x39
    x41 = x3 * (x12 + x40)
    x42 = x27 + x29
    x43 = x2 * x42
    x44 = 2 * x41 + 2 * x43
    x45 = x3 * (3 * x17 + x23 + x38 + x44)
    x46 = x36 + x45
    x47 = x17 + x37
    x48 = x2 * x47
    x49 = x41 + x43
    x50 = x2 * x49
    x51 = x2 * x20
    x52 = x2 * x26
    x53 = 2 * x52
    x54 = x3 * (x21 + x33 + x51 + x53)
    x55 = x14 + x39
    x56 = x2 * x55
    x57 = x18 * x5 + x56
    x58 = x3 * (x27 + 3 * x29 + x53 + x57)
    x59 = x2 * x46 + x3 * (2 * x25 + 2 * x34 + 2 * x48 + 2 * x50 + 2 * x54 + 2 * x58)
    x60 = x48 + x54
    x61 = x2 * x60
    x62 = x50 + x58
    x63 = x2 * x62
    x64 = x2 * x28
    x65 = 2 * x64
    x66 = x3 * (x16 + x65)
    x67 = x19 + x51
    x68 = x2 * x67
    x69 = x66 + x68
    x70 = x10 * x2
    x71 = x3 * (x11 + x15 + x64 + x70)
    x72 = x29 + x52
    x73 = x2 * x72
    x74 = 2 * x71 + 2 * x73
    x75 = x3 * (2 * x17 + x38 + x69 + x74)
    x76 = 2 * x70
    x77 = x3 * (x40 + x76)
    x78 = x2 * x57
    x79 = x77 + x78
    x80 = x3 * (x44 + x74 + x79)
    x81 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x82 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x83 = numpy.pi * x0 * x82
    x84 = x81 * x83
    x85 = -x0 * (a * A[1] + b * B[1])
    x86 = -x85 - B[1]
    x87 = x61 + x75
    x88 = x2 * x9
    x89 = x3 * (x28 + x88)
    x90 = x14 + x64
    x91 = x2 * x90
    x92 = x2 * x69 + x3 * (x32 + 2 * x51 + 2 * x89 + 2 * x91)
    x93 = x89 + x91
    x94 = x3 * (x10 + x88)
    x95 = x14 + x70
    x96 = x2 * x95
    x97 = x94 + x96
    x98 = x3 * (x30 + x53 + x93 + x97)
    x99 = x71 + x73
    x100 = x2 * x99
    x101 = 2 * x100 + 2 * x98
    x102 = x84 * (x2 * x87 + x3 * (x101 + 3 * x48 + 3 * x54 + x92))
    x103 = -x0 * (a * A[2] + b * B[2])
    x104 = -x103 - B[2]
    x105 = x8 * x81
    x106 = x105 * x3
    x107 = x105 * x86 ** 2
    x108 = x106 + x107
    x109 = x2 ** 2 * x9
    x110 = x109 + x15
    x111 = x3 * (x110 + x65)
    x112 = x2 * x93
    x113 = x2 * x92 + x3 * (2 * x111 + 2 * x112 + 3 * x66 + 3 * x68)
    x114 = x8 * x82
    x115 = x104 * x84
    x116 = x114 * x3
    x117 = x104 ** 2 * x114
    x118 = x116 + x117
    x119 = -x85 - A[1]
    x120 = x59 * x84
    x121 = x105 * x119
    x122 = x121 * x86
    x123 = x106 + x122
    x124 = 2 * x106
    x125 = x108 * x119
    x126 = x124 * x86 + x125
    x127 = x104 * x114
    x128 = -x103 - A[2]
    x129 = x128 * x84
    x130 = x114 * x128
    x131 = x104 * x130
    x132 = x116 + x131
    x133 = x105 * x86
    x134 = 2 * x116
    x135 = x118 * x128
    x136 = x104 * x134 + x135
    x137 = x105 * x119 ** 2
    x138 = x106 + x137
    x139 = x3 * (x121 + x133)
    x140 = x119 * x123
    x141 = x139 + x140
    x142 = 3 * x106
    x143 = 2 * x122 + x142
    x144 = x3 * (x107 + x143)
    x145 = x119 * x126
    x146 = x144 + x145
    x147 = x114 * x128 ** 2
    x148 = x116 + x147
    x149 = x3 * (x127 + x130)
    x150 = x128 * x132
    x151 = x149 + x150
    x152 = 3 * x116
    x153 = 2 * x131 + x152
    x154 = x3 * (x117 + x153)
    x155 = x128 * x136
    x156 = x154 + x155
    x157 = x119 * x124 + x119 * x138
    x158 = x3 * (x137 + x143)
    x159 = x119 * x141
    x160 = x158 + x159
    x161 = 4 * x106
    x162 = x119 * x146 + x3 * (2 * x125 + 2 * x139 + 2 * x140 + x161 * x86)
    x163 = x128 * x134 + x128 * x148
    x164 = x3 * (x147 + x153)
    x165 = x128 * x151
    x166 = x164 + x165
    x167 = 4 * x116
    x168 = x128 * x156 + x3 * (x104 * x167 + 2 * x135 + 2 * x149 + 2 * x150)
    x169 = x119 * x157 + x3 * (3 * x137 + x142)
    x170 = x119 * x160 + x3 * (3 * x139 + 3 * x140 + x157)
    x171 = x119 * x162 + x3 * (3 * x144 + 3 * x145 + 2 * x158 + 2 * x159)
    x172 = x128 * x163 + x3 * (3 * x147 + x152)
    x173 = x128 * x166 + x3 * (3 * x149 + 3 * x150 + x163)
    x174 = x128 * x168 + x3 * (3 * x154 + 3 * x155 + 2 * x164 + 2 * x165)
    x175 = -x85 - C[1]
    x176 = x63 + x80
    x177 = x2 * x79 + x3 * (x31 * x5 + 2 * x56 + 2 * x94 + 2 * x96)
    x178 = x84 * (x176 * x2 + x3 * (x101 + x177 + 3 * x50 + 3 * x58))
    x179 = x133 * x175
    x180 = x106 + x179
    x181 = x100 + x98
    x182 = x111 + x112
    x183 = x3 * (x110 + x76)
    x184 = x2 * x97
    x185 = x183 + x184
    x186 = x181 * x2 + x3 * (x182 + x185 + 3 * x71 + 3 * x73)
    x187 = x105 * x175
    x188 = x3 * (x133 + x187)
    x189 = x180 * x86
    x190 = x188 + x189
    x191 = x109 + x14
    x192 = x18 * x2 + x191 * x2
    x193 = x182 * x2 + x3 * (x192 + 3 * x89 + 3 * x91)
    x194 = x121 * x175
    x195 = x106 + x194
    x196 = x119 * x180
    x197 = x188 + x196
    x198 = x142 + 2 * x179
    x199 = x3 * (x107 + x198)
    x200 = x119 * x190
    x201 = x199 + x200
    x202 = x3 * (x121 + x187)
    x203 = x119 * x195
    x204 = x202 + x203
    x205 = x3 * (x122 + x142 + x179 + x194)
    x206 = x119 * x197
    x207 = x205 + x206
    x208 = x119 * x201
    x209 = 2 * x196
    x210 = x3 * (x126 + 3 * x188 + x189 + x209)
    x211 = x208 + x210
    x212 = x142 + 2 * x194
    x213 = x3 * (x137 + x212)
    x214 = x119 * x204
    x215 = x213 + x214
    x216 = x119 * x207
    x217 = 2 * x188
    x218 = x3 * (x141 + x204 + x209 + x217)
    x219 = x216 + x218
    x220 = x119 * x211
    x221 = 2 * x199 + 2 * x200
    x222 = 2 * x205 + 2 * x206
    x223 = x3 * (x146 + x221 + x222)
    x224 = x220 + x223
    x225 = x119 * x215 + x3 * (x157 + 3 * x202 + 3 * x203)
    x226 = x119 * x219 + x3 * (x160 + 3 * x205 + 3 * x206 + x215)
    x227 = 2 * x216 + 2 * x218
    x228 = x7 * x83
    x229 = x228 * (x119 * x224 + x3 * (x162 + 3 * x208 + 3 * x210 + x227))
    x230 = x228 * x4
    x231 = numpy.pi * x0 * x7 * x81
    x232 = x231 * x4
    x233 = -x103 - C[2]
    x234 = x233 * x84
    x235 = x127 * x233
    x236 = x116 + x235
    x237 = x114 * x233
    x238 = x3 * (x127 + x237)
    x239 = x104 * x236
    x240 = x238 + x239
    x241 = x130 * x233
    x242 = x116 + x241
    x243 = x128 * x236
    x244 = x238 + x243
    x245 = x152 + 2 * x235
    x246 = x3 * (x117 + x245)
    x247 = x128 * x240
    x248 = x246 + x247
    x249 = x3 * (x130 + x237)
    x250 = x128 * x242
    x251 = x249 + x250
    x252 = x3 * (x131 + x152 + x235 + x241)
    x253 = x128 * x244
    x254 = x252 + x253
    x255 = x128 * x248
    x256 = 2 * x243
    x257 = x3 * (x136 + 3 * x238 + x239 + x256)
    x258 = x255 + x257
    x259 = x152 + 2 * x241
    x260 = x3 * (x147 + x259)
    x261 = x128 * x251
    x262 = x260 + x261
    x263 = x128 * x254
    x264 = 2 * x238
    x265 = x3 * (x151 + x251 + x256 + x264)
    x266 = x263 + x265
    x267 = x128 * x258
    x268 = 2 * x246 + 2 * x247
    x269 = 2 * x252 + 2 * x253
    x270 = x3 * (x156 + x268 + x269)
    x271 = x267 + x270
    x272 = x128 * x262 + x3 * (x163 + 3 * x249 + 3 * x250)
    x273 = x128 * x266 + x3 * (x166 + 3 * x252 + 3 * x253 + x262)
    x274 = 2 * x263 + 2 * x265
    x275 = x231 * (x128 * x271 + x3 * (x168 + 3 * x255 + 3 * x257 + x274))
    x276 = x105 * x175 ** 2
    x277 = x106 + x276
    x278 = x177 * x2 + x3 * (2 * x183 + 2 * x184 + 3 * x77 + 3 * x78)
    x279 = x124 * x175
    x280 = x277 * x86
    x281 = x279 + x280
    x282 = x185 * x2 + x3 * (x192 + 3 * x94 + 3 * x96)
    x283 = x3 * (x198 + x276)
    x284 = x281 * x86
    x285 = x283 + x284
    x286 = x192 * x2 + x3 * (3 * x109 + x15)
    x287 = x119 * x277
    x288 = x279 + x287
    x289 = x119 * x281
    x290 = x283 + x289
    x291 = x119 * x285
    x292 = x161 * x175
    x293 = x217 + x292
    x294 = x3 * (2 * x189 + 2 * x280 + x293)
    x295 = x291 + x294
    x296 = x3 * (x212 + x276)
    x297 = x119 * x288
    x298 = x296 + x297
    x299 = x119 * x290
    x300 = x3 * (x209 + x280 + x287 + x293)
    x301 = x299 + x300
    x302 = x119 * x295
    x303 = 2 * x289
    x304 = x3 * (x221 + 3 * x283 + x284 + x303)
    x305 = x302 + x304
    x306 = x119 * x298 + x3 * (2 * x202 + 2 * x203 + 2 * x287 + x292)
    x307 = x119 * x301
    x308 = x3 * (x222 + 2 * x283 + x298 + x303)
    x309 = x307 + x308
    x310 = x119 * x305 + x3 * (
        2 * x208 + 2 * x210 + 2 * x291 + 2 * x294 + 2 * x299 + 2 * x300
    )
    x311 = x228 * x310
    x312 = x2 * x228
    x313 = x119 * x306 + x3 * (2 * x213 + 2 * x214 + 3 * x296 + 3 * x297)
    x314 = x228 * (x119 * x309 + x3 * (x227 + 3 * x299 + 3 * x300 + x306))
    x315 = x228 * x5
    x316 = x175 * x231
    x317 = x114 * x233 ** 2
    x318 = x116 + x317
    x319 = x134 * x233
    x320 = x104 * x318
    x321 = x319 + x320
    x322 = x3 * (x245 + x317)
    x323 = x104 * x321
    x324 = x322 + x323
    x325 = x128 * x318
    x326 = x319 + x325
    x327 = x128 * x321
    x328 = x322 + x327
    x329 = x128 * x324
    x330 = x167 * x233
    x331 = x264 + x330
    x332 = x3 * (2 * x239 + 2 * x320 + x331)
    x333 = x329 + x332
    x334 = x3 * (x259 + x317)
    x335 = x128 * x326
    x336 = x334 + x335
    x337 = x128 * x328
    x338 = x3 * (x256 + x320 + x325 + x331)
    x339 = x337 + x338
    x340 = x128 * x333
    x341 = 2 * x327
    x342 = x3 * (x268 + 3 * x322 + x323 + x341)
    x343 = x340 + x342
    x344 = x2 * x231
    x345 = x128 * x336 + x3 * (2 * x249 + 2 * x250 + 2 * x325 + x330)
    x346 = x128 * x339
    x347 = x3 * (x269 + 2 * x322 + x336 + x341)
    x348 = x346 + x347
    x349 = x128 * x343 + x3 * (
        2 * x255 + 2 * x257 + 2 * x329 + 2 * x332 + 2 * x337 + 2 * x338
    )
    x350 = x231 * x349
    x351 = x231 * x5
    x352 = x128 * x345 + x3 * (2 * x260 + 2 * x261 + 3 * x334 + 3 * x335)
    x353 = x231 * (x128 * x348 + x3 * (x274 + 3 * x337 + 3 * x338 + x345))

    # 540 item(s)
    return numpy.array(
        [
            x84
            * (
                x2 * x59
                + x3 * (3 * x36 + 3 * x45 + 2 * x61 + 2 * x63 + 2 * x75 + 2 * x80)
            ),
            x102 * x86,
            x102 * x104,
            x108 * x113 * x114,
            x113 * x115 * x86,
            x105 * x113 * x118,
            x119 * x120,
            x114 * x123 * x87,
            x115 * x119 * x87,
            x114 * x126 * x92,
            x123 * x127 * x92,
            x118 * x121 * x92,
            x120 * x128,
            x129 * x86 * x87,
            x105 * x132 * x87,
            x108 * x130 * x92,
            x132 * x133 * x92,
            x105 * x136 * x92,
            x114 * x138 * x46,
            x114 * x141 * x60,
            x127 * x138 * x60,
            x114 * x146 * x69,
            x127 * x141 * x69,
            x118 * x138 * x69,
            x119 * x129 * x46,
            x123 * x130 * x60,
            x121 * x132 * x60,
            x126 * x130 * x69,
            x123 * x132 * x69,
            x121 * x136 * x69,
            x105 * x148 * x46,
            x133 * x148 * x60,
            x105 * x151 * x60,
            x108 * x148 * x69,
            x133 * x151 * x69,
            x105 * x156 * x69,
            x114 * x157 * x35,
            x114 * x160 * x47,
            x127 * x157 * x47,
            x114 * x162 * x67,
            x127 * x160 * x67,
            x118 * x157 * x67,
            x130 * x138 * x35,
            x130 * x141 * x47,
            x132 * x138 * x47,
            x130 * x146 * x67,
            x132 * x141 * x67,
            x136 * x138 * x67,
            x121 * x148 * x35,
            x123 * x148 * x47,
            x121 * x151 * x47,
            x126 * x148 * x67,
            x123 * x151 * x67,
            x121 * x156 * x67,
            x105 * x163 * x35,
            x133 * x163 * x47,
            x105 * x166 * x47,
            x108 * x163 * x67,
            x133 * x166 * x67,
            x105 * x168 * x67,
            x114 * x169 * x24,
            x114 * x170 * x22,
            x127 * x169 * x22,
            x114 * x171 * x20,
            x127 * x170 * x20,
            x118 * x169 * x20,
            x130 * x157 * x24,
            x130 * x160 * x22,
            x132 * x157 * x22,
            x130 * x162 * x20,
            x132 * x160 * x20,
            x136 * x157 * x20,
            x138 * x148 * x24,
            x141 * x148 * x22,
            x138 * x151 * x22,
            x146 * x148 * x20,
            x141 * x151 * x20,
            x138 * x156 * x20,
            x121 * x163 * x24,
            x123 * x163 * x22,
            x121 * x166 * x22,
            x126 * x163 * x20,
            x123 * x166 * x20,
            x121 * x168 * x20,
            x105 * x172 * x24,
            x133 * x172 * x22,
            x105 * x173 * x22,
            x108 * x172 * x20,
            x133 * x173 * x20,
            x105 * x174 * x20,
            x175 * x178,
            x114 * x180 * x186,
            x115 * x175 * x186,
            x114 * x190 * x193,
            x127 * x180 * x193,
            x118 * x187 * x193,
            x114 * x176 * x195,
            x114 * x181 * x197,
            x127 * x181 * x195,
            x114 * x182 * x201,
            x127 * x182 * x197,
            x118 * x182 * x195,
            x129 * x175 * x176,
            x130 * x180 * x181,
            x132 * x181 * x187,
            x130 * x182 * x190,
            x132 * x180 * x182,
            x136 * x182 * x187,
            x114 * x204 * x62,
            x114 * x207 * x99,
            x127 * x204 * x99,
            x114 * x211 * x93,
            x127 * x207 * x93,
            x118 * x204 * x93,
            x130 * x195 * x62,
            x130 * x197 * x99,
            x132 * x195 * x99,
            x130 * x201 * x93,
            x132 * x197 * x93,
            x136 * x195 * x93,
            x148 * x187 * x62,
            x148 * x180 * x99,
            x151 * x187 * x99,
            x148 * x190 * x93,
            x151 * x180 * x93,
            x156 * x187 * x93,
            x114 * x215 * x49,
            x114 * x219 * x72,
            x127 * x215 * x72,
            x114 * x224 * x90,
            x127 * x219 * x90,
            x118 * x215 * x90,
            x130 * x204 * x49,
            x130 * x207 * x72,
            x132 * x204 * x72,
            x130 * x211 * x90,
            x132 * x207 * x90,
            x136 * x204 * x90,
            x148 * x195 * x49,
            x148 * x197 * x72,
            x151 * x195 * x72,
            x148 * x201 * x90,
            x151 * x197 * x90,
            x156 * x195 * x90,
            x163 * x187 * x49,
            x163 * x180 * x72,
            x166 * x187 * x72,
            x163 * x190 * x90,
            x166 * x180 * x90,
            x168 * x187 * x90,
            x114 * x225 * x42,
            x114 * x226 * x26,
            x127 * x225 * x26,
            x229 * x4,
            x104 * x226 * x230,
            x118 * x225 * x28,
            x130 * x215 * x42,
            x130 * x219 * x26,
            x132 * x215 * x26,
            x128 * x224 * x230,
            x132 * x219 * x28,
            x136 * x215 * x28,
            x148 * x204 * x42,
            x148 * x207 * x26,
            x151 * x204 * x26,
            x148 * x211 * x28,
            x151 * x207 * x28,
            x156 * x204 * x28,
            x163 * x195 * x42,
            x163 * x197 * x26,
            x166 * x195 * x26,
            x163 * x201 * x28,
            x166 * x197 * x28,
            x168 * x195 * x28,
            x172 * x187 * x42,
            x172 * x180 * x26,
            x173 * x187 * x26,
            x172 * x190 * x28,
            x173 * x180 * x28,
            x174 * x175 * x232,
            x178 * x233,
            x186 * x234 * x86,
            x105 * x186 * x236,
            x108 * x193 * x237,
            x133 * x193 * x236,
            x105 * x193 * x240,
            x119 * x176 * x234,
            x123 * x181 * x237,
            x121 * x181 * x236,
            x126 * x182 * x237,
            x123 * x182 * x236,
            x121 * x182 * x240,
            x105 * x176 * x242,
            x133 * x181 * x242,
            x105 * x181 * x244,
            x108 * x182 * x242,
            x133 * x182 * x244,
            x105 * x182 * x248,
            x138 * x237 * x62,
            x141 * x237 * x99,
            x138 * x236 * x99,
            x146 * x237 * x93,
            x141 * x236 * x93,
            x138 * x240 * x93,
            x121 * x242 * x62,
            x123 * x242 * x99,
            x121 * x244 * x99,
            x126 * x242 * x93,
            x123 * x244 * x93,
            x121 * x248 * x93,
            x105 * x251 * x62,
            x133 * x251 * x99,
            x105 * x254 * x99,
            x108 * x251 * x93,
            x133 * x254 * x93,
            x105 * x258 * x93,
            x157 * x237 * x49,
            x160 * x237 * x72,
            x157 * x236 * x72,
            x162 * x237 * x90,
            x160 * x236 * x90,
            x157 * x240 * x90,
            x138 * x242 * x49,
            x141 * x242 * x72,
            x138 * x244 * x72,
            x146 * x242 * x90,
            x141 * x244 * x90,
            x138 * x248 * x90,
            x121 * x251 * x49,
            x123 * x251 * x72,
            x121 * x254 * x72,
            x126 * x251 * x90,
            x123 * x254 * x90,
            x121 * x258 * x90,
            x105 * x262 * x49,
            x133 * x262 * x72,
            x105 * x266 * x72,
            x108 * x262 * x90,
            x133 * x266 * x90,
            x105 * x271 * x90,
            x169 * x237 * x42,
            x170 * x237 * x26,
            x169 * x236 * x26,
            x171 * x230 * x233,
            x170 * x236 * x28,
            x169 * x240 * x28,
            x157 * x242 * x42,
            x160 * x242 * x26,
            x157 * x244 * x26,
            x162 * x242 * x28,
            x160 * x244 * x28,
            x157 * x248 * x28,
            x138 * x251 * x42,
            x141 * x251 * x26,
            x138 * x254 * x26,
            x146 * x251 * x28,
            x141 * x254 * x28,
            x138 * x258 * x28,
            x121 * x262 * x42,
            x123 * x26 * x262,
            x121 * x26 * x266,
            x126 * x262 * x28,
            x123 * x266 * x28,
            x119 * x232 * x271,
            x105 * x272 * x42,
            x133 * x26 * x272,
            x105 * x26 * x273,
            x108 * x272 * x28,
            x232 * x273 * x86,
            x275 * x4,
            x114 * x277 * x278,
            x114 * x281 * x282,
            x127 * x277 * x282,
            x114 * x285 * x286,
            x127 * x281 * x286,
            x118 * x277 * x286,
            x114 * x177 * x288,
            x114 * x185 * x290,
            x127 * x185 * x288,
            x114 * x192 * x295,
            x127 * x192 * x290,
            x118 * x192 * x288,
            x130 * x177 * x277,
            x130 * x185 * x281,
            x132 * x185 * x277,
            x130 * x192 * x285,
            x132 * x192 * x281,
            x136 * x192 * x277,
            x114 * x298 * x79,
            x114 * x301 * x97,
            x127 * x298 * x97,
            x114 * x191 * x305,
            x127 * x191 * x301,
            x118 * x191 * x298,
            x130 * x288 * x79,
            x130 * x290 * x97,
            x132 * x288 * x97,
            x130 * x191 * x295,
            x132 * x191 * x290,
            x136 * x191 * x288,
            x148 * x277 * x79,
            x148 * x281 * x97,
            x151 * x277 * x97,
            x148 * x191 * x285,
            x151 * x191 * x281,
            x156 * x191 * x277,
            x114 * x306 * x57,
            x114 * x309 * x95,
            x127 * x306 * x95,
            x2 * x311,
            x104 * x309 * x312,
            x118 * x306 * x88,
            x130 * x298 * x57,
            x130 * x301 * x95,
            x132 * x298 * x95,
            x128 * x305 * x312,
            x132 * x301 * x88,
            x136 * x298 * x88,
            x148 * x288 * x57,
            x148 * x290 * x95,
            x151 * x288 * x95,
            x148 * x295 * x88,
            x151 * x290 * x88,
            x156 * x288 * x88,
            x163 * x277 * x57,
            x163 * x281 * x95,
            x166 * x277 * x95,
            x163 * x285 * x88,
            x166 * x281 * x88,
            x168 * x277 * x88,
            x114 * x313 * x55,
            x314 * x5,
            x104 * x313 * x315,
            x228
            * (
                x119 * x310
                + x3 * (2 * x220 + 2 * x223 + 3 * x302 + 3 * x304 + 2 * x307 + 2 * x308)
            ),
            x104 * x314,
            x118 * x313 * x9,
            x130 * x306 * x55,
            x128 * x309 * x315,
            x10 * x132 * x306,
            x128 * x311,
            x132 * x309 * x9,
            x136 * x306 * x9,
            x148 * x298 * x55,
            x10 * x148 * x301,
            x10 * x151 * x298,
            x148 * x305 * x9,
            x151 * x301 * x9,
            x156 * x298 * x9,
            x163 * x288 * x55,
            x10 * x163 * x290,
            x10 * x166 * x288,
            x163 * x295 * x9,
            x166 * x290 * x9,
            x168 * x288 * x9,
            x172 * x277 * x55,
            x10 * x172 * x281,
            x10 * x173 * x277,
            x172 * x285 * x9,
            x173 * x281 * x9,
            x174 * x277 * x9,
            x175 * x234 * x278,
            x180 * x237 * x282,
            x187 * x236 * x282,
            x190 * x237 * x286,
            x180 * x236 * x286,
            x187 * x240 * x286,
            x177 * x195 * x237,
            x185 * x197 * x237,
            x185 * x195 * x236,
            x192 * x201 * x237,
            x192 * x197 * x236,
            x192 * x195 * x240,
            x177 * x187 * x242,
            x180 * x185 * x242,
            x185 * x187 * x244,
            x190 * x192 * x242,
            x180 * x192 * x244,
            x187 * x192 * x248,
            x204 * x237 * x79,
            x207 * x237 * x97,
            x204 * x236 * x97,
            x191 * x211 * x237,
            x191 * x207 * x236,
            x191 * x204 * x240,
            x195 * x242 * x79,
            x197 * x242 * x97,
            x195 * x244 * x97,
            x191 * x201 * x242,
            x191 * x197 * x244,
            x191 * x195 * x248,
            x187 * x251 * x79,
            x180 * x251 * x97,
            x187 * x254 * x97,
            x190 * x191 * x251,
            x180 * x191 * x254,
            x187 * x191 * x258,
            x215 * x237 * x57,
            x219 * x237 * x95,
            x215 * x236 * x95,
            x224 * x233 * x312,
            x219 * x236 * x88,
            x215 * x240 * x88,
            x204 * x242 * x57,
            x207 * x242 * x95,
            x204 * x244 * x95,
            x211 * x242 * x88,
            x207 * x244 * x88,
            x204 * x248 * x88,
            x195 * x251 * x57,
            x197 * x251 * x95,
            x195 * x254 * x95,
            x201 * x251 * x88,
            x197 * x254 * x88,
            x195 * x258 * x88,
            x187 * x262 * x57,
            x180 * x262 * x95,
            x187 * x266 * x95,
            x190 * x262 * x88,
            x180 * x266 * x88,
            x2 * x271 * x316,
            x225 * x237 * x55,
            x226 * x233 * x315,
            x10 * x225 * x236,
            x229 * x233,
            x226 * x236 * x9,
            x225 * x240 * x9,
            x215 * x242 * x55,
            x10 * x219 * x242,
            x10 * x215 * x244,
            x224 * x242 * x9,
            x219 * x244 * x9,
            x215 * x248 * x9,
            x204 * x251 * x55,
            x10 * x207 * x251,
            x10 * x204 * x254,
            x211 * x251 * x9,
            x207 * x254 * x9,
            x204 * x258 * x9,
            x195 * x262 * x55,
            x10 * x197 * x262,
            x10 * x195 * x266,
            x201 * x262 * x9,
            x197 * x266 * x9,
            x195 * x271 * x9,
            x187 * x272 * x55,
            x10 * x180 * x272,
            x273 * x316 * x5,
            x190 * x272 * x9,
            x180 * x273 * x9,
            x175 * x275,
            x105 * x278 * x318,
            x133 * x282 * x318,
            x105 * x282 * x321,
            x108 * x286 * x318,
            x133 * x286 * x321,
            x105 * x286 * x324,
            x121 * x177 * x318,
            x123 * x185 * x318,
            x121 * x185 * x321,
            x126 * x192 * x318,
            x123 * x192 * x321,
            x121 * x192 * x324,
            x105 * x177 * x326,
            x133 * x185 * x326,
            x105 * x185 * x328,
            x108 * x192 * x326,
            x133 * x192 * x328,
            x105 * x192 * x333,
            x138 * x318 * x79,
            x141 * x318 * x97,
            x138 * x321 * x97,
            x146 * x191 * x318,
            x141 * x191 * x321,
            x138 * x191 * x324,
            x121 * x326 * x79,
            x123 * x326 * x97,
            x121 * x328 * x97,
            x126 * x191 * x326,
            x123 * x191 * x328,
            x121 * x191 * x333,
            x105 * x336 * x79,
            x133 * x336 * x97,
            x105 * x339 * x97,
            x108 * x191 * x336,
            x133 * x191 * x339,
            x105 * x191 * x343,
            x157 * x318 * x57,
            x160 * x318 * x95,
            x157 * x321 * x95,
            x162 * x318 * x88,
            x160 * x321 * x88,
            x157 * x324 * x88,
            x138 * x326 * x57,
            x141 * x326 * x95,
            x138 * x328 * x95,
            x146 * x326 * x88,
            x141 * x328 * x88,
            x138 * x333 * x88,
            x121 * x336 * x57,
            x123 * x336 * x95,
            x121 * x339 * x95,
            x126 * x336 * x88,
            x123 * x339 * x88,
            x119 * x343 * x344,
            x105 * x345 * x57,
            x133 * x345 * x95,
            x105 * x348 * x95,
            x108 * x345 * x88,
            x344 * x348 * x86,
            x2 * x350,
            x169 * x318 * x55,
            x10 * x170 * x318,
            x10 * x169 * x321,
            x171 * x318 * x9,
            x170 * x321 * x9,
            x169 * x324 * x9,
            x157 * x326 * x55,
            x10 * x160 * x326,
            x10 * x157 * x328,
            x162 * x326 * x9,
            x160 * x328 * x9,
            x157 * x333 * x9,
            x138 * x336 * x55,
            x10 * x141 * x336,
            x10 * x138 * x339,
            x146 * x336 * x9,
            x141 * x339 * x9,
            x138 * x343 * x9,
            x121 * x345 * x55,
            x10 * x123 * x345,
            x119 * x348 * x351,
            x126 * x345 * x9,
            x123 * x348 * x9,
            x119 * x350,
            x105 * x352 * x55,
            x351 * x352 * x86,
            x353 * x5,
            x108 * x352 * x9,
            x353 * x86,
            x231
            * (
                x128 * x349
                + x3 * (2 * x267 + 2 * x270 + 3 * x340 + 3 * x342 + 2 * x346 + 2 * x347)
            ),
        ]
    )


def quadrupole3d_43(a, A, b, B, C):
    """Cartesian 3D (gf) quadrupole moment integrals.
    The origin is at C.

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
    x9 = x4 * x8
    x10 = -x1 - C[0]
    x11 = x10 * x8
    x12 = x3 * (x11 + x9)
    x13 = x3 * x8
    x14 = x10 * x9
    x15 = x13 + x14
    x16 = x15 * x4
    x17 = x12 + x16
    x18 = x17 * x4
    x19 = x10 * x13
    x20 = 2 * x19
    x21 = x10 ** 2 * x8
    x22 = x13 + x21
    x23 = x22 * x4
    x24 = x20 + x23
    x25 = x24 * x4
    x26 = x4 ** 2 * x8
    x27 = 3 * x13
    x28 = 2 * x14 + x27
    x29 = x3 * (x26 + x28)
    x30 = 2 * x29
    x31 = x3 * (x21 + x28)
    x32 = x30 + 3 * x31
    x33 = x3 * (2 * x18 + 3 * x25 + x32)
    x34 = x25 + x31
    x35 = x34 * x4
    x36 = 2 * x12
    x37 = 4 * x19
    x38 = x36 + x37
    x39 = x3 * (2 * x16 + 2 * x23 + x38)
    x40 = x35 + x39
    x41 = x2 * x40
    x42 = x33 + x41
    x43 = x2 * x42
    x44 = x2 * x34
    x45 = 3 * x12
    x46 = x13 * x4
    x47 = 2 * x46
    x48 = x13 + x26
    x49 = x4 * x48
    x50 = x47 + x49
    x51 = x3 * (3 * x16 + x45 + x50)
    x52 = x18 + x29
    x53 = x2 * x52
    x54 = 2 * x51 + 2 * x53
    x55 = x3 * (x35 + 4 * x39 + 3 * x44 + x54)
    x56 = x43 + x55
    x57 = x17 * x2
    x58 = x3 * (3 * x26 + x27)
    x59 = x2 * x50
    x60 = x58 + x59
    x61 = x3 * (x18 + 4 * x29 + 3 * x57 + x60)
    x62 = x51 + x53
    x63 = x2 * x62
    x64 = 2 * x57
    x65 = x2 * x24
    x66 = 2 * x65
    x67 = x3 * (x25 + x32 + x64 + x66)
    x68 = x39 + x44
    x69 = x2 * x68
    x70 = 3 * x67 + 3 * x69
    x71 = x2 * x56 + x3 * (2 * x33 + 2 * x41 + 2 * x61 + 2 * x63 + x70)
    x72 = x67 + x69
    x73 = x2 * x72
    x74 = x61 + x63
    x75 = x2 * x74
    x76 = x2 * x48
    x77 = x3 * (8 * x46 + x49 + 3 * x76)
    x78 = x2 * x60
    x79 = x77 + x78
    x80 = x15 * x2
    x81 = 2 * x80
    x82 = x47 + x76
    x83 = x3 * (x16 + x45 + x81 + x82)
    x84 = x29 + x57
    x85 = x2 * x84
    x86 = 3 * x83 + 3 * x85
    x87 = x3 * (x54 + x79 + x86)
    x88 = x31 + x65
    x89 = x2 * x88
    x90 = x2 * x22
    x91 = x3 * (x23 + x38 + x81 + x90)
    x92 = x3 * (2 * x39 + 2 * x44 + 2 * x83 + 2 * x85 + 2 * x89 + 2 * x91)
    x93 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x94 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x95 = numpy.pi * x0 * x94
    x96 = x93 * x95
    x97 = -x0 * (a * A[1] + b * B[1])
    x98 = -x97 - B[1]
    x99 = x73 + x92
    x100 = x2 * x9
    x101 = 2 * x100 + x27
    x102 = x3 * (x101 + x26)
    x103 = x2 * x82
    x104 = x102 + x103
    x105 = x11 * x2
    x106 = x3 * (x100 + x105 + x14 + x27)
    x107 = x12 + x80
    x108 = x107 * x2
    x109 = 2 * x106 + 2 * x108
    x110 = x3 * (x104 + x109 + x30 + x64)
    x111 = 2 * x105 + x27
    x112 = x3 * (x111 + x21)
    x113 = x20 + x90
    x114 = x113 * x2
    x115 = x112 + x114
    x116 = x3 * (x109 + x115 + 2 * x31 + x66)
    x117 = x83 + x85
    x118 = x117 * x2
    x119 = x89 + x91
    x120 = x119 * x2
    x121 = x96 * (x2 * x99 + x3 * (2 * x110 + 2 * x116 + 2 * x118 + 2 * x120 + x70))
    x122 = -x0 * (a * A[2] + b * B[2])
    x123 = -x122 - B[2]
    x124 = x7 * x93
    x125 = x124 * x3
    x126 = x124 * x98 ** 2
    x127 = x125 + x126
    x128 = x116 + x120
    x129 = x2 * x8
    x130 = x3 * (x11 + x129)
    x131 = x105 + x13
    x132 = x131 * x2
    x133 = x115 * x2 + x3 * (2 * x130 + 2 * x132 + x37 + 2 * x90)
    x134 = x130 + x132
    x135 = x3 * (x129 + x9)
    x136 = x100 + x13
    x137 = x136 * x2
    x138 = x135 + x137
    x139 = x3 * (x134 + x138 + x36 + x81)
    x140 = x106 + x108
    x141 = x140 * x2
    x142 = 2 * x139 + 2 * x141
    x143 = x128 * x2 + x3 * (x133 + x142 + 3 * x89 + 3 * x91)
    x144 = x7 * x94
    x145 = x123 * x96
    x146 = x144 * x3
    x147 = x123 ** 2 * x144
    x148 = x146 + x147
    x149 = x125 * x98
    x150 = 2 * x149
    x151 = x127 * x98
    x152 = x150 + x151
    x153 = x2 ** 2 * x8
    x154 = x3 * (x111 + x153)
    x155 = x134 * x2
    x156 = x133 * x2 + x3 * (3 * x112 + 3 * x114 + 2 * x154 + 2 * x155)
    x157 = x123 * x144
    x158 = x124 * x98
    x159 = x123 * x146
    x160 = 2 * x159
    x161 = x123 * x148
    x162 = x160 + x161
    x163 = -x97 - A[1]
    x164 = x71 * x96
    x165 = x124 * x163
    x166 = x165 * x98
    x167 = x125 + x166
    x168 = x127 * x163
    x169 = x150 + x168
    x170 = 3 * x125
    x171 = x3 * (3 * x126 + x170)
    x172 = x152 * x163
    x173 = x171 + x172
    x174 = -x122 - A[2]
    x175 = x174 * x96
    x176 = x144 * x174
    x177 = x123 * x176
    x178 = x146 + x177
    x179 = x148 * x174
    x180 = x160 + x179
    x181 = 3 * x146
    x182 = x3 * (3 * x147 + x181)
    x183 = x162 * x174
    x184 = x182 + x183
    x185 = x124 * x163 ** 2
    x186 = x125 + x185
    x187 = x3 * (x158 + x165)
    x188 = x163 * x167
    x189 = x187 + x188
    x190 = 2 * x166 + x170
    x191 = x3 * (x126 + x190)
    x192 = x163 * x169
    x193 = x191 + x192
    x194 = x3 * (8 * x149 + x151 + 3 * x168)
    x195 = x163 * x173
    x196 = x194 + x195
    x197 = x144 * x174 ** 2
    x198 = x146 + x197
    x199 = x3 * (x157 + x176)
    x200 = x174 * x178
    x201 = x199 + x200
    x202 = 2 * x177 + x181
    x203 = x3 * (x147 + x202)
    x204 = x174 * x180
    x205 = x203 + x204
    x206 = x3 * (8 * x159 + x161 + 3 * x179)
    x207 = x174 * x184
    x208 = x206 + x207
    x209 = 2 * x125
    x210 = x163 * x186 + x163 * x209
    x211 = x3 * (x185 + x190)
    x212 = x163 * x189
    x213 = x211 + x212
    x214 = x163 * x193
    x215 = x3 * (4 * x149 + 2 * x168 + 2 * x187 + 2 * x188)
    x216 = x214 + x215
    x217 = 3 * x191 + 3 * x192
    x218 = x163 * x196 + x3 * (2 * x171 + 2 * x172 + x217)
    x219 = 2 * x146
    x220 = x174 * x198 + x174 * x219
    x221 = x3 * (x197 + x202)
    x222 = x174 * x201
    x223 = x221 + x222
    x224 = x174 * x205
    x225 = x3 * (4 * x159 + 2 * x179 + 2 * x199 + 2 * x200)
    x226 = x224 + x225
    x227 = 3 * x203 + 3 * x204
    x228 = x174 * x208 + x3 * (2 * x182 + 2 * x183 + x227)
    x229 = x163 * x210 + x3 * (x170 + 3 * x185)
    x230 = x163 * x213 + x3 * (3 * x187 + 3 * x188 + x210)
    x231 = x163 * x216 + x3 * (2 * x211 + 2 * x212 + x217)
    x232 = x163 * x218 + x3 * (3 * x194 + 3 * x195 + 3 * x214 + 3 * x215)
    x233 = x174 * x220 + x3 * (x181 + 3 * x197)
    x234 = x174 * x223 + x3 * (3 * x199 + 3 * x200 + x220)
    x235 = x174 * x226 + x3 * (2 * x221 + 2 * x222 + x227)
    x236 = x174 * x228 + x3 * (3 * x206 + 3 * x207 + 3 * x224 + 3 * x225)
    x237 = -x97 - C[1]
    x238 = x75 + x87
    x239 = 3 * x102 + 3 * x103
    x240 = x2 * x79 + x3 * (x239 + 2 * x58 + 2 * x59)
    x241 = x96 * (x2 * x238 + x3 * (3 * x110 + 3 * x118 + x240 + 3 * x61 + 3 * x63))
    x242 = x158 * x237
    x243 = x125 + x242
    x244 = x110 + x118
    x245 = x104 * x2
    x246 = x3 * (2 * x135 + 2 * x137 + 4 * x46 + 2 * x76)
    x247 = x245 + x246
    x248 = x2 * x244 + x3 * (x142 + x247 + x86)
    x249 = x124 * x237
    x250 = x3 * (x158 + x249)
    x251 = x243 * x98
    x252 = x250 + x251
    x253 = x139 + x141
    x254 = x154 + x155
    x255 = x3 * (x101 + x153)
    x256 = x138 * x2
    x257 = x255 + x256
    x258 = x2 * x253 + x3 * (3 * x106 + 3 * x108 + x254 + x257)
    x259 = x170 + 2 * x242
    x260 = x3 * (x126 + x259)
    x261 = x252 * x98
    x262 = x260 + x261
    x263 = x13 + x153
    x264 = 2 * x13 * x2 + x2 * x263
    x265 = x2 * x254 + x3 * (3 * x130 + 3 * x132 + x264)
    x266 = x165 * x237
    x267 = x125 + x266
    x268 = x163 * x243
    x269 = x250 + x268
    x270 = x163 * x252
    x271 = x260 + x270
    x272 = 3 * x250
    x273 = x3 * (x152 + 3 * x251 + x272)
    x274 = x163 * x262
    x275 = x273 + x274
    x276 = x3 * (x165 + x249)
    x277 = x163 * x267
    x278 = x276 + x277
    x279 = x3 * (x166 + x170 + x242 + x266)
    x280 = x163 * x269
    x281 = x279 + x280
    x282 = x163 * x271
    x283 = 2 * x268
    x284 = x3 * (x169 + x251 + x272 + x283)
    x285 = x282 + x284
    x286 = x163 * x275
    x287 = x3 * (x173 + 4 * x260 + x261 + 3 * x270)
    x288 = x286 + x287
    x289 = x170 + 2 * x266
    x290 = x3 * (x185 + x289)
    x291 = x163 * x278
    x292 = x290 + x291
    x293 = x163 * x281
    x294 = 2 * x250
    x295 = x3 * (x189 + x278 + x283 + x294)
    x296 = x293 + x295
    x297 = x163 * x285
    x298 = 2 * x260
    x299 = 2 * x270
    x300 = 2 * x279 + 2 * x280
    x301 = x3 * (x193 + x298 + x299 + x300)
    x302 = x297 + x301
    x303 = x163 * x288
    x304 = 2 * x273 + 2 * x274
    x305 = 3 * x282 + 3 * x284
    x306 = x3 * (x196 + x304 + x305)
    x307 = x303 + x306
    x308 = x163 * x292 + x3 * (x210 + 3 * x276 + 3 * x277)
    x309 = x163 * x296 + x3 * (x213 + 3 * x279 + 3 * x280 + x292)
    x310 = 2 * x293 + 2 * x295
    x311 = x163 * x302 + x3 * (x216 + x305 + x310)
    x312 = x6 * x95
    x313 = x312 * (x163 * x307 + x3 * (x218 + 3 * x286 + 3 * x287 + 3 * x297 + 3 * x301))
    x314 = x10 * x312
    x315 = numpy.pi * x0 * x6 * x93
    x316 = x10 * x315
    x317 = -x122 - C[2]
    x318 = x317 * x96
    x319 = x157 * x317
    x320 = x146 + x319
    x321 = x144 * x317
    x322 = x3 * (x157 + x321)
    x323 = x123 * x320
    x324 = x322 + x323
    x325 = x181 + 2 * x319
    x326 = x3 * (x147 + x325)
    x327 = x123 * x324
    x328 = x326 + x327
    x329 = x176 * x317
    x330 = x146 + x329
    x331 = x174 * x320
    x332 = x322 + x331
    x333 = x174 * x324
    x334 = x326 + x333
    x335 = 3 * x322
    x336 = x3 * (x162 + 3 * x323 + x335)
    x337 = x174 * x328
    x338 = x336 + x337
    x339 = x3 * (x176 + x321)
    x340 = x174 * x330
    x341 = x339 + x340
    x342 = x3 * (x177 + x181 + x319 + x329)
    x343 = x174 * x332
    x344 = x342 + x343
    x345 = x174 * x334
    x346 = 2 * x331
    x347 = x3 * (x180 + x323 + x335 + x346)
    x348 = x345 + x347
    x349 = x174 * x338
    x350 = x3 * (x184 + 4 * x326 + x327 + 3 * x333)
    x351 = x349 + x350
    x352 = x181 + 2 * x329
    x353 = x3 * (x197 + x352)
    x354 = x174 * x341
    x355 = x353 + x354
    x356 = x174 * x344
    x357 = 2 * x322
    x358 = x3 * (x201 + x341 + x346 + x357)
    x359 = x356 + x358
    x360 = x174 * x348
    x361 = 2 * x326
    x362 = 2 * x333
    x363 = 2 * x342 + 2 * x343
    x364 = x3 * (x205 + x361 + x362 + x363)
    x365 = x360 + x364
    x366 = x174 * x351
    x367 = 2 * x336 + 2 * x337
    x368 = 3 * x345 + 3 * x347
    x369 = x3 * (x208 + x367 + x368)
    x370 = x366 + x369
    x371 = x174 * x355 + x3 * (x220 + 3 * x339 + 3 * x340)
    x372 = x174 * x359 + x3 * (x223 + 3 * x342 + 3 * x343 + x355)
    x373 = 2 * x356 + 2 * x358
    x374 = x174 * x365 + x3 * (x226 + x368 + x373)
    x375 = x315 * (x174 * x370 + x3 * (x228 + 3 * x349 + 3 * x350 + 3 * x360 + 3 * x364))
    x376 = x124 * x237 ** 2
    x377 = x125 + x376
    x378 = x2 * x240 + x3 * (3 * x245 + 3 * x246 + 3 * x77 + 3 * x78)
    x379 = x209 * x237
    x380 = x377 * x98
    x381 = x379 + x380
    x382 = x2 * x247 + x3 * (x239 + 2 * x255 + 2 * x256)
    x383 = x3 * (x259 + x376)
    x384 = x381 * x98
    x385 = x383 + x384
    x386 = x2 * x257 + x3 * (3 * x135 + 3 * x137 + x264)
    x387 = x385 * x98
    x388 = 4 * x125 * x237
    x389 = x294 + x388
    x390 = x3 * (2 * x251 + 2 * x380 + x389)
    x391 = x387 + x390
    x392 = x2 * x264 + x3 * (3 * x153 + x27)
    x393 = x163 * x377
    x394 = x379 + x393
    x395 = x163 * x381
    x396 = x383 + x395
    x397 = x163 * x385
    x398 = x390 + x397
    x399 = x298 + 3 * x383
    x400 = x3 * (2 * x261 + 3 * x384 + x399)
    x401 = x163 * x391
    x402 = x400 + x401
    x403 = x3 * (x289 + x376)
    x404 = x163 * x394
    x405 = x403 + x404
    x406 = x163 * x396
    x407 = x3 * (x283 + x380 + x389 + x393)
    x408 = x406 + x407
    x409 = x163 * x398
    x410 = 2 * x395
    x411 = x3 * (x299 + x384 + x399 + x410)
    x412 = x409 + x411
    x413 = x163 * x402
    x414 = x3 * (x304 + x387 + 4 * x390 + 3 * x397)
    x415 = x413 + x414
    x416 = x163 * x405 + x3 * (2 * x276 + 2 * x277 + x388 + 2 * x393)
    x417 = x163 * x408
    x418 = x3 * (x300 + 2 * x383 + x405 + x410)
    x419 = x417 + x418
    x420 = x163 * x412
    x421 = x3 * (2 * x282 + 2 * x284 + 2 * x390 + 2 * x397 + 2 * x406 + 2 * x407)
    x422 = x420 + x421
    x423 = 3 * x409 + 3 * x411
    x424 = x163 * x415 + x3 * (2 * x286 + 2 * x287 + 2 * x400 + 2 * x401 + x423)
    x425 = x312 * x424
    x426 = x2 * x312
    x427 = x163 * x416 + x3 * (2 * x290 + 2 * x291 + 3 * x403 + 3 * x404)
    x428 = x163 * x419 + x3 * (x310 + 3 * x406 + 3 * x407 + x416)
    x429 = x312 * (x163 * x422 + x3 * (2 * x297 + 2 * x301 + 2 * x417 + 2 * x418 + x423))
    x430 = x312 * x4
    x431 = x237 * x315
    x432 = x144 * x317 ** 2
    x433 = x146 + x432
    x434 = x219 * x317
    x435 = x123 * x433
    x436 = x434 + x435
    x437 = x3 * (x325 + x432)
    x438 = x123 * x436
    x439 = x437 + x438
    x440 = x123 * x439
    x441 = 4 * x146 * x317
    x442 = x357 + x441
    x443 = x3 * (2 * x323 + 2 * x435 + x442)
    x444 = x440 + x443
    x445 = x174 * x433
    x446 = x434 + x445
    x447 = x174 * x436
    x448 = x437 + x447
    x449 = x174 * x439
    x450 = x443 + x449
    x451 = x361 + 3 * x437
    x452 = x3 * (2 * x327 + 3 * x438 + x451)
    x453 = x174 * x444
    x454 = x452 + x453
    x455 = x3 * (x352 + x432)
    x456 = x174 * x446
    x457 = x455 + x456
    x458 = x174 * x448
    x459 = x3 * (x346 + x435 + x442 + x445)
    x460 = x458 + x459
    x461 = x174 * x450
    x462 = 2 * x447
    x463 = x3 * (x362 + x438 + x451 + x462)
    x464 = x461 + x463
    x465 = x174 * x454
    x466 = x3 * (x367 + x440 + 4 * x443 + 3 * x449)
    x467 = x465 + x466
    x468 = x2 * x315
    x469 = x174 * x457 + x3 * (2 * x339 + 2 * x340 + x441 + 2 * x445)
    x470 = x174 * x460
    x471 = x3 * (x363 + 2 * x437 + x457 + x462)
    x472 = x470 + x471
    x473 = x174 * x464
    x474 = x3 * (2 * x345 + 2 * x347 + 2 * x443 + 2 * x449 + 2 * x458 + 2 * x459)
    x475 = x473 + x474
    x476 = 3 * x461 + 3 * x463
    x477 = x174 * x467 + x3 * (2 * x349 + 2 * x350 + 2 * x452 + 2 * x453 + x476)
    x478 = x315 * x477
    x479 = x315 * x4
    x480 = x174 * x469 + x3 * (2 * x353 + 2 * x354 + 3 * x455 + 3 * x456)
    x481 = x174 * x472 + x3 * (x373 + 3 * x458 + 3 * x459 + x469)
    x482 = x315 * (x174 * x475 + x3 * (2 * x360 + 2 * x364 + 2 * x470 + 2 * x471 + x476))

    # 900 item(s)
    return numpy.array(
        [
            x96
            * (
                x2 * x71
                + x3 * (3 * x43 + 3 * x55 + 3 * x73 + 2 * x75 + 2 * x87 + 3 * x92)
            ),
            x121 * x98,
            x121 * x123,
            x127 * x143 * x144,
            x143 * x145 * x98,
            x124 * x143 * x148,
            x144 * x152 * x156,
            x127 * x156 * x157,
            x148 * x156 * x158,
            x124 * x156 * x162,
            x163 * x164,
            x144 * x167 * x99,
            x145 * x163 * x99,
            x128 * x144 * x169,
            x128 * x157 * x167,
            x128 * x148 * x165,
            x133 * x144 * x173,
            x133 * x157 * x169,
            x133 * x148 * x167,
            x133 * x162 * x165,
            x164 * x174,
            x175 * x98 * x99,
            x124 * x178 * x99,
            x127 * x128 * x176,
            x128 * x158 * x178,
            x124 * x128 * x180,
            x133 * x152 * x176,
            x127 * x133 * x178,
            x133 * x158 * x180,
            x124 * x133 * x184,
            x144 * x186 * x56,
            x144 * x189 * x72,
            x157 * x186 * x72,
            x119 * x144 * x193,
            x119 * x157 * x189,
            x119 * x148 * x186,
            x115 * x144 * x196,
            x115 * x157 * x193,
            x115 * x148 * x189,
            x115 * x162 * x186,
            x163 * x175 * x56,
            x167 * x176 * x72,
            x165 * x178 * x72,
            x119 * x169 * x176,
            x119 * x167 * x178,
            x119 * x165 * x180,
            x115 * x173 * x176,
            x115 * x169 * x178,
            x115 * x167 * x180,
            x115 * x165 * x184,
            x124 * x198 * x56,
            x158 * x198 * x72,
            x124 * x201 * x72,
            x119 * x127 * x198,
            x119 * x158 * x201,
            x119 * x124 * x205,
            x115 * x152 * x198,
            x115 * x127 * x201,
            x115 * x158 * x205,
            x115 * x124 * x208,
            x144 * x210 * x42,
            x144 * x213 * x68,
            x157 * x210 * x68,
            x144 * x216 * x88,
            x157 * x213 * x88,
            x148 * x210 * x88,
            x113 * x144 * x218,
            x113 * x157 * x216,
            x113 * x148 * x213,
            x113 * x162 * x210,
            x176 * x186 * x42,
            x176 * x189 * x68,
            x178 * x186 * x68,
            x176 * x193 * x88,
            x178 * x189 * x88,
            x180 * x186 * x88,
            x113 * x176 * x196,
            x113 * x178 * x193,
            x113 * x180 * x189,
            x113 * x184 * x186,
            x165 * x198 * x42,
            x167 * x198 * x68,
            x165 * x201 * x68,
            x169 * x198 * x88,
            x167 * x201 * x88,
            x165 * x205 * x88,
            x113 * x173 * x198,
            x113 * x169 * x201,
            x113 * x167 * x205,
            x113 * x165 * x208,
            x124 * x220 * x42,
            x158 * x220 * x68,
            x124 * x223 * x68,
            x127 * x220 * x88,
            x158 * x223 * x88,
            x124 * x226 * x88,
            x113 * x152 * x220,
            x113 * x127 * x223,
            x113 * x158 * x226,
            x113 * x124 * x228,
            x144 * x229 * x40,
            x144 * x230 * x34,
            x157 * x229 * x34,
            x144 * x231 * x24,
            x157 * x230 * x24,
            x148 * x229 * x24,
            x144 * x22 * x232,
            x157 * x22 * x231,
            x148 * x22 * x230,
            x162 * x22 * x229,
            x176 * x210 * x40,
            x176 * x213 * x34,
            x178 * x210 * x34,
            x176 * x216 * x24,
            x178 * x213 * x24,
            x180 * x210 * x24,
            x176 * x218 * x22,
            x178 * x216 * x22,
            x180 * x213 * x22,
            x184 * x210 * x22,
            x186 * x198 * x40,
            x189 * x198 * x34,
            x186 * x201 * x34,
            x193 * x198 * x24,
            x189 * x201 * x24,
            x186 * x205 * x24,
            x196 * x198 * x22,
            x193 * x201 * x22,
            x189 * x205 * x22,
            x186 * x208 * x22,
            x165 * x220 * x40,
            x167 * x220 * x34,
            x165 * x223 * x34,
            x169 * x220 * x24,
            x167 * x223 * x24,
            x165 * x226 * x24,
            x173 * x22 * x220,
            x169 * x22 * x223,
            x167 * x22 * x226,
            x165 * x22 * x228,
            x124 * x233 * x40,
            x158 * x233 * x34,
            x124 * x234 * x34,
            x127 * x233 * x24,
            x158 * x234 * x24,
            x124 * x235 * x24,
            x152 * x22 * x233,
            x127 * x22 * x234,
            x158 * x22 * x235,
            x124 * x22 * x236,
            x237 * x241,
            x144 * x243 * x248,
            x145 * x237 * x248,
            x144 * x252 * x258,
            x157 * x243 * x258,
            x148 * x249 * x258,
            x144 * x262 * x265,
            x157 * x252 * x265,
            x148 * x243 * x265,
            x162 * x249 * x265,
            x144 * x238 * x267,
            x144 * x244 * x269,
            x157 * x244 * x267,
            x144 * x253 * x271,
            x157 * x253 * x269,
            x148 * x253 * x267,
            x144 * x254 * x275,
            x157 * x254 * x271,
            x148 * x254 * x269,
            x162 * x254 * x267,
            x175 * x237 * x238,
            x176 * x243 * x244,
            x178 * x244 * x249,
            x176 * x252 * x253,
            x178 * x243 * x253,
            x180 * x249 * x253,
            x176 * x254 * x262,
            x178 * x252 * x254,
            x180 * x243 * x254,
            x184 * x249 * x254,
            x144 * x278 * x74,
            x117 * x144 * x281,
            x117 * x157 * x278,
            x140 * x144 * x285,
            x140 * x157 * x281,
            x140 * x148 * x278,
            x134 * x144 * x288,
            x134 * x157 * x285,
            x134 * x148 * x281,
            x134 * x162 * x278,
            x176 * x267 * x74,
            x117 * x176 * x269,
            x117 * x178 * x267,
            x140 * x176 * x271,
            x140 * x178 * x269,
            x140 * x180 * x267,
            x134 * x176 * x275,
            x134 * x178 * x271,
            x134 * x180 * x269,
            x134 * x184 * x267,
            x198 * x249 * x74,
            x117 * x198 * x243,
            x117 * x201 * x249,
            x140 * x198 * x252,
            x140 * x201 * x243,
            x140 * x205 * x249,
            x134 * x198 * x262,
            x134 * x201 * x252,
            x134 * x205 * x243,
            x134 * x208 * x249,
            x144 * x292 * x62,
            x144 * x296 * x84,
            x157 * x292 * x84,
            x107 * x144 * x302,
            x107 * x157 * x296,
            x107 * x148 * x292,
            x131 * x144 * x307,
            x131 * x157 * x302,
            x131 * x148 * x296,
            x131 * x162 * x292,
            x176 * x278 * x62,
            x176 * x281 * x84,
            x178 * x278 * x84,
            x107 * x176 * x285,
            x107 * x178 * x281,
            x107 * x180 * x278,
            x131 * x176 * x288,
            x131 * x178 * x285,
            x131 * x180 * x281,
            x131 * x184 * x278,
            x198 * x267 * x62,
            x198 * x269 * x84,
            x201 * x267 * x84,
            x107 * x198 * x271,
            x107 * x201 * x269,
            x107 * x205 * x267,
            x131 * x198 * x275,
            x131 * x201 * x271,
            x131 * x205 * x269,
            x131 * x208 * x267,
            x220 * x249 * x62,
            x220 * x243 * x84,
            x223 * x249 * x84,
            x107 * x220 * x252,
            x107 * x223 * x243,
            x107 * x226 * x249,
            x131 * x220 * x262,
            x131 * x223 * x252,
            x131 * x226 * x243,
            x131 * x228 * x249,
            x144 * x308 * x52,
            x144 * x17 * x309,
            x157 * x17 * x308,
            x144 * x15 * x311,
            x15 * x157 * x309,
            x148 * x15 * x308,
            x10 * x313,
            x123 * x311 * x314,
            x11 * x148 * x309,
            x11 * x162 * x308,
            x176 * x292 * x52,
            x17 * x176 * x296,
            x17 * x178 * x292,
            x15 * x176 * x302,
            x15 * x178 * x296,
            x15 * x180 * x292,
            x174 * x307 * x314,
            x11 * x178 * x302,
            x11 * x180 * x296,
            x11 * x184 * x292,
            x198 * x278 * x52,
            x17 * x198 * x281,
            x17 * x201 * x278,
            x15 * x198 * x285,
            x15 * x201 * x281,
            x15 * x205 * x278,
            x11 * x198 * x288,
            x11 * x201 * x285,
            x11 * x205 * x281,
            x11 * x208 * x278,
            x220 * x267 * x52,
            x17 * x220 * x269,
            x17 * x223 * x267,
            x15 * x220 * x271,
            x15 * x223 * x269,
            x15 * x226 * x267,
            x11 * x220 * x275,
            x11 * x223 * x271,
            x11 * x226 * x269,
            x11 * x228 * x267,
            x233 * x249 * x52,
            x17 * x233 * x243,
            x17 * x234 * x249,
            x15 * x233 * x252,
            x15 * x234 * x243,
            x15 * x235 * x249,
            x11 * x233 * x262,
            x11 * x234 * x252,
            x11 * x235 * x243,
            x236 * x237 * x316,
            x241 * x317,
            x248 * x318 * x98,
            x124 * x248 * x320,
            x127 * x258 * x321,
            x158 * x258 * x320,
            x124 * x258 * x324,
            x152 * x265 * x321,
            x127 * x265 * x320,
            x158 * x265 * x324,
            x124 * x265 * x328,
            x163 * x238 * x318,
            x167 * x244 * x321,
            x165 * x244 * x320,
            x169 * x253 * x321,
            x167 * x253 * x320,
            x165 * x253 * x324,
            x173 * x254 * x321,
            x169 * x254 * x320,
            x167 * x254 * x324,
            x165 * x254 * x328,
            x124 * x238 * x330,
            x158 * x244 * x330,
            x124 * x244 * x332,
            x127 * x253 * x330,
            x158 * x253 * x332,
            x124 * x253 * x334,
            x152 * x254 * x330,
            x127 * x254 * x332,
            x158 * x254 * x334,
            x124 * x254 * x338,
            x186 * x321 * x74,
            x117 * x189 * x321,
            x117 * x186 * x320,
            x140 * x193 * x321,
            x140 * x189 * x320,
            x140 * x186 * x324,
            x134 * x196 * x321,
            x134 * x193 * x320,
            x134 * x189 * x324,
            x134 * x186 * x328,
            x165 * x330 * x74,
            x117 * x167 * x330,
            x117 * x165 * x332,
            x140 * x169 * x330,
            x140 * x167 * x332,
            x140 * x165 * x334,
            x134 * x173 * x330,
            x134 * x169 * x332,
            x134 * x167 * x334,
            x134 * x165 * x338,
            x124 * x341 * x74,
            x117 * x158 * x341,
            x117 * x124 * x344,
            x127 * x140 * x341,
            x140 * x158 * x344,
            x124 * x140 * x348,
            x134 * x152 * x341,
            x127 * x134 * x344,
            x134 * x158 * x348,
            x124 * x134 * x351,
            x210 * x321 * x62,
            x213 * x321 * x84,
            x210 * x320 * x84,
            x107 * x216 * x321,
            x107 * x213 * x320,
            x107 * x210 * x324,
            x131 * x218 * x321,
            x131 * x216 * x320,
            x131 * x213 * x324,
            x131 * x210 * x328,
            x186 * x330 * x62,
            x189 * x330 * x84,
            x186 * x332 * x84,
            x107 * x193 * x330,
            x107 * x189 * x332,
            x107 * x186 * x334,
            x131 * x196 * x330,
            x131 * x193 * x332,
            x131 * x189 * x334,
            x131 * x186 * x338,
            x165 * x341 * x62,
            x167 * x341 * x84,
            x165 * x344 * x84,
            x107 * x169 * x341,
            x107 * x167 * x344,
            x107 * x165 * x348,
            x131 * x173 * x341,
            x131 * x169 * x344,
            x131 * x167 * x348,
            x131 * x165 * x351,
            x124 * x355 * x62,
            x158 * x355 * x84,
            x124 * x359 * x84,
            x107 * x127 * x355,
            x107 * x158 * x359,
            x107 * x124 * x365,
            x131 * x152 * x355,
            x127 * x131 * x359,
            x131 * x158 * x365,
            x124 * x131 * x370,
            x229 * x321 * x52,
            x17 * x230 * x321,
            x17 * x229 * x320,
            x15 * x231 * x321,
            x15 * x230 * x320,
            x15 * x229 * x324,
            x232 * x314 * x317,
            x11 * x231 * x320,
            x11 * x230 * x324,
            x11 * x229 * x328,
            x210 * x330 * x52,
            x17 * x213 * x330,
            x17 * x210 * x332,
            x15 * x216 * x330,
            x15 * x213 * x332,
            x15 * x210 * x334,
            x11 * x218 * x330,
            x11 * x216 * x332,
            x11 * x213 * x334,
            x11 * x210 * x338,
            x186 * x341 * x52,
            x17 * x189 * x341,
            x17 * x186 * x344,
            x15 * x193 * x341,
            x15 * x189 * x344,
            x15 * x186 * x348,
            x11 * x196 * x341,
            x11 * x193 * x344,
            x11 * x189 * x348,
            x11 * x186 * x351,
            x165 * x355 * x52,
            x167 * x17 * x355,
            x165 * x17 * x359,
            x15 * x169 * x355,
            x15 * x167 * x359,
            x15 * x165 * x365,
            x11 * x173 * x355,
            x11 * x169 * x359,
            x11 * x167 * x365,
            x163 * x316 * x370,
            x124 * x371 * x52,
            x158 * x17 * x371,
            x124 * x17 * x372,
            x127 * x15 * x371,
            x15 * x158 * x372,
            x124 * x15 * x374,
            x11 * x152 * x371,
            x11 * x127 * x372,
            x316 * x374 * x98,
            x10 * x375,
            x144 * x377 * x378,
            x144 * x381 * x382,
            x157 * x377 * x382,
            x144 * x385 * x386,
            x157 * x381 * x386,
            x148 * x377 * x386,
            x144 * x391 * x392,
            x157 * x385 * x392,
            x148 * x381 * x392,
            x162 * x377 * x392,
            x144 * x240 * x394,
            x144 * x247 * x396,
            x157 * x247 * x394,
            x144 * x257 * x398,
            x157 * x257 * x396,
            x148 * x257 * x394,
            x144 * x264 * x402,
            x157 * x264 * x398,
            x148 * x264 * x396,
            x162 * x264 * x394,
            x176 * x240 * x377,
            x176 * x247 * x381,
            x178 * x247 * x377,
            x176 * x257 * x385,
            x178 * x257 * x381,
            x180 * x257 * x377,
            x176 * x264 * x391,
            x178 * x264 * x385,
            x180 * x264 * x381,
            x184 * x264 * x377,
            x144 * x405 * x79,
            x104 * x144 * x408,
            x104 * x157 * x405,
            x138 * x144 * x412,
            x138 * x157 * x408,
            x138 * x148 * x405,
            x144 * x263 * x415,
            x157 * x263 * x412,
            x148 * x263 * x408,
            x162 * x263 * x405,
            x176 * x394 * x79,
            x104 * x176 * x396,
            x104 * x178 * x394,
            x138 * x176 * x398,
            x138 * x178 * x396,
            x138 * x180 * x394,
            x176 * x263 * x402,
            x178 * x263 * x398,
            x180 * x263 * x396,
            x184 * x263 * x394,
            x198 * x377 * x79,
            x104 * x198 * x381,
            x104 * x201 * x377,
            x138 * x198 * x385,
            x138 * x201 * x381,
            x138 * x205 * x377,
            x198 * x263 * x391,
            x201 * x263 * x385,
            x205 * x263 * x381,
            x208 * x263 * x377,
            x144 * x416 * x60,
            x144 * x419 * x82,
            x157 * x416 * x82,
            x136 * x144 * x422,
            x136 * x157 * x419,
            x136 * x148 * x416,
            x2 * x425,
            x123 * x422 * x426,
            x129 * x148 * x419,
            x129 * x162 * x416,
            x176 * x405 * x60,
            x176 * x408 * x82,
            x178 * x405 * x82,
            x136 * x176 * x412,
            x136 * x178 * x408,
            x136 * x180 * x405,
            x174 * x415 * x426,
            x129 * x178 * x412,
            x129 * x180 * x408,
            x129 * x184 * x405,
            x198 * x394 * x60,
            x198 * x396 * x82,
            x201 * x394 * x82,
            x136 * x198 * x398,
            x136 * x201 * x396,
            x136 * x205 * x394,
            x129 * x198 * x402,
            x129 * x201 * x398,
            x129 * x205 * x396,
            x129 * x208 * x394,
            x220 * x377 * x60,
            x220 * x381 * x82,
            x223 * x377 * x82,
            x136 * x220 * x385,
            x136 * x223 * x381,
            x136 * x226 * x377,
            x129 * x220 * x391,
            x129 * x223 * x385,
            x129 * x226 * x381,
            x129 * x228 * x377,
            x144 * x427 * x50,
            x144 * x428 * x48,
            x157 * x427 * x48,
            x4 * x429,
            x123 * x428 * x430,
            x148 * x427 * x9,
            x312
            * (
                x163 * x424
                + x3 * (2 * x303 + 2 * x306 + 3 * x413 + 3 * x414 + 3 * x420 + 3 * x421)
            ),
            x123 * x429,
            x148 * x428 * x8,
            x162 * x427 * x8,
            x176 * x416 * x50,
            x176 * x419 * x48,
            x178 * x416 * x48,
            x174 * x422 * x430,
            x178 * x419 * x9,
            x180 * x416 * x9,
            x174 * x425,
            x178 * x422 * x8,
            x180 * x419 * x8,
            x184 * x416 * x8,
            x198 * x405 * x50,
            x198 * x408 * x48,
            x201 * x405 * x48,
            x198 * x412 * x9,
            x201 * x408 * x9,
            x205 * x405 * x9,
            x198 * x415 * x8,
            x201 * x412 * x8,
            x205 * x408 * x8,
            x208 * x405 * x8,
            x220 * x394 * x50,
            x220 * x396 * x48,
            x223 * x394 * x48,
            x220 * x398 * x9,
            x223 * x396 * x9,
            x226 * x394 * x9,
            x220 * x402 * x8,
            x223 * x398 * x8,
            x226 * x396 * x8,
            x228 * x394 * x8,
            x233 * x377 * x50,
            x233 * x381 * x48,
            x234 * x377 * x48,
            x233 * x385 * x9,
            x234 * x381 * x9,
            x235 * x377 * x9,
            x233 * x391 * x8,
            x234 * x385 * x8,
            x235 * x381 * x8,
            x236 * x377 * x8,
            x237 * x318 * x378,
            x243 * x321 * x382,
            x249 * x320 * x382,
            x252 * x321 * x386,
            x243 * x320 * x386,
            x249 * x324 * x386,
            x262 * x321 * x392,
            x252 * x320 * x392,
            x243 * x324 * x392,
            x249 * x328 * x392,
            x240 * x267 * x321,
            x247 * x269 * x321,
            x247 * x267 * x320,
            x257 * x271 * x321,
            x257 * x269 * x320,
            x257 * x267 * x324,
            x264 * x275 * x321,
            x264 * x271 * x320,
            x264 * x269 * x324,
            x264 * x267 * x328,
            x240 * x249 * x330,
            x243 * x247 * x330,
            x247 * x249 * x332,
            x252 * x257 * x330,
            x243 * x257 * x332,
            x249 * x257 * x334,
            x262 * x264 * x330,
            x252 * x264 * x332,
            x243 * x264 * x334,
            x249 * x264 * x338,
            x278 * x321 * x79,
            x104 * x281 * x321,
            x104 * x278 * x320,
            x138 * x285 * x321,
            x138 * x281 * x320,
            x138 * x278 * x324,
            x263 * x288 * x321,
            x263 * x285 * x320,
            x263 * x281 * x324,
            x263 * x278 * x328,
            x267 * x330 * x79,
            x104 * x269 * x330,
            x104 * x267 * x332,
            x138 * x271 * x330,
            x138 * x269 * x332,
            x138 * x267 * x334,
            x263 * x275 * x330,
            x263 * x271 * x332,
            x263 * x269 * x334,
            x263 * x267 * x338,
            x249 * x341 * x79,
            x104 * x243 * x341,
            x104 * x249 * x344,
            x138 * x252 * x341,
            x138 * x243 * x344,
            x138 * x249 * x348,
            x262 * x263 * x341,
            x252 * x263 * x344,
            x243 * x263 * x348,
            x249 * x263 * x351,
            x292 * x321 * x60,
            x296 * x321 * x82,
            x292 * x320 * x82,
            x136 * x302 * x321,
            x136 * x296 * x320,
            x136 * x292 * x324,
            x307 * x317 * x426,
            x129 * x302 * x320,
            x129 * x296 * x324,
            x129 * x292 * x328,
            x278 * x330 * x60,
            x281 * x330 * x82,
            x278 * x332 * x82,
            x136 * x285 * x330,
            x136 * x281 * x332,
            x136 * x278 * x334,
            x129 * x288 * x330,
            x129 * x285 * x332,
            x129 * x281 * x334,
            x129 * x278 * x338,
            x267 * x341 * x60,
            x269 * x341 * x82,
            x267 * x344 * x82,
            x136 * x271 * x341,
            x136 * x269 * x344,
            x136 * x267 * x348,
            x129 * x275 * x341,
            x129 * x271 * x344,
            x129 * x269 * x348,
            x129 * x267 * x351,
            x249 * x355 * x60,
            x243 * x355 * x82,
            x249 * x359 * x82,
            x136 * x252 * x355,
            x136 * x243 * x359,
            x136 * x249 * x365,
            x129 * x262 * x355,
            x129 * x252 * x359,
            x129 * x243 * x365,
            x2 * x370 * x431,
            x308 * x321 * x50,
            x309 * x321 * x48,
            x308 * x320 * x48,
            x311 * x317 * x430,
            x309 * x320 * x9,
            x308 * x324 * x9,
            x313 * x317,
            x311 * x320 * x8,
            x309 * x324 * x8,
            x308 * x328 * x8,
            x292 * x330 * x50,
            x296 * x330 * x48,
            x292 * x332 * x48,
            x302 * x330 * x9,
            x296 * x332 * x9,
            x292 * x334 * x9,
            x307 * x330 * x8,
            x302 * x332 * x8,
            x296 * x334 * x8,
            x292 * x338 * x8,
            x278 * x341 * x50,
            x281 * x341 * x48,
            x278 * x344 * x48,
            x285 * x341 * x9,
            x281 * x344 * x9,
            x278 * x348 * x9,
            x288 * x341 * x8,
            x285 * x344 * x8,
            x281 * x348 * x8,
            x278 * x351 * x8,
            x267 * x355 * x50,
            x269 * x355 * x48,
            x267 * x359 * x48,
            x271 * x355 * x9,
            x269 * x359 * x9,
            x267 * x365 * x9,
            x275 * x355 * x8,
            x271 * x359 * x8,
            x269 * x365 * x8,
            x267 * x370 * x8,
            x249 * x371 * x50,
            x243 * x371 * x48,
            x249 * x372 * x48,
            x252 * x371 * x9,
            x243 * x372 * x9,
            x374 * x4 * x431,
            x262 * x371 * x8,
            x252 * x372 * x8,
            x243 * x374 * x8,
            x237 * x375,
            x124 * x378 * x433,
            x158 * x382 * x433,
            x124 * x382 * x436,
            x127 * x386 * x433,
            x158 * x386 * x436,
            x124 * x386 * x439,
            x152 * x392 * x433,
            x127 * x392 * x436,
            x158 * x392 * x439,
            x124 * x392 * x444,
            x165 * x240 * x433,
            x167 * x247 * x433,
            x165 * x247 * x436,
            x169 * x257 * x433,
            x167 * x257 * x436,
            x165 * x257 * x439,
            x173 * x264 * x433,
            x169 * x264 * x436,
            x167 * x264 * x439,
            x165 * x264 * x444,
            x124 * x240 * x446,
            x158 * x247 * x446,
            x124 * x247 * x448,
            x127 * x257 * x446,
            x158 * x257 * x448,
            x124 * x257 * x450,
            x152 * x264 * x446,
            x127 * x264 * x448,
            x158 * x264 * x450,
            x124 * x264 * x454,
            x186 * x433 * x79,
            x104 * x189 * x433,
            x104 * x186 * x436,
            x138 * x193 * x433,
            x138 * x189 * x436,
            x138 * x186 * x439,
            x196 * x263 * x433,
            x193 * x263 * x436,
            x189 * x263 * x439,
            x186 * x263 * x444,
            x165 * x446 * x79,
            x104 * x167 * x446,
            x104 * x165 * x448,
            x138 * x169 * x446,
            x138 * x167 * x448,
            x138 * x165 * x450,
            x173 * x263 * x446,
            x169 * x263 * x448,
            x167 * x263 * x450,
            x165 * x263 * x454,
            x124 * x457 * x79,
            x104 * x158 * x457,
            x104 * x124 * x460,
            x127 * x138 * x457,
            x138 * x158 * x460,
            x124 * x138 * x464,
            x152 * x263 * x457,
            x127 * x263 * x460,
            x158 * x263 * x464,
            x124 * x263 * x467,
            x210 * x433 * x60,
            x213 * x433 * x82,
            x210 * x436 * x82,
            x136 * x216 * x433,
            x136 * x213 * x436,
            x136 * x210 * x439,
            x129 * x218 * x433,
            x129 * x216 * x436,
            x129 * x213 * x439,
            x129 * x210 * x444,
            x186 * x446 * x60,
            x189 * x446 * x82,
            x186 * x448 * x82,
            x136 * x193 * x446,
            x136 * x189 * x448,
            x136 * x186 * x450,
            x129 * x196 * x446,
            x129 * x193 * x448,
            x129 * x189 * x450,
            x129 * x186 * x454,
            x165 * x457 * x60,
            x167 * x457 * x82,
            x165 * x460 * x82,
            x136 * x169 * x457,
            x136 * x167 * x460,
            x136 * x165 * x464,
            x129 * x173 * x457,
            x129 * x169 * x460,
            x129 * x167 * x464,
            x163 * x467 * x468,
            x124 * x469 * x60,
            x158 * x469 * x82,
            x124 * x472 * x82,
            x127 * x136 * x469,
            x136 * x158 * x472,
            x124 * x136 * x475,
            x129 * x152 * x469,
            x127 * x129 * x472,
            x468 * x475 * x98,
            x2 * x478,
            x229 * x433 * x50,
            x230 * x433 * x48,
            x229 * x436 * x48,
            x231 * x433 * x9,
            x230 * x436 * x9,
            x229 * x439 * x9,
            x232 * x433 * x8,
            x231 * x436 * x8,
            x230 * x439 * x8,
            x229 * x444 * x8,
            x210 * x446 * x50,
            x213 * x446 * x48,
            x210 * x448 * x48,
            x216 * x446 * x9,
            x213 * x448 * x9,
            x210 * x450 * x9,
            x218 * x446 * x8,
            x216 * x448 * x8,
            x213 * x450 * x8,
            x210 * x454 * x8,
            x186 * x457 * x50,
            x189 * x457 * x48,
            x186 * x460 * x48,
            x193 * x457 * x9,
            x189 * x460 * x9,
            x186 * x464 * x9,
            x196 * x457 * x8,
            x193 * x460 * x8,
            x189 * x464 * x8,
            x186 * x467 * x8,
            x165 * x469 * x50,
            x167 * x469 * x48,
            x165 * x472 * x48,
            x169 * x469 * x9,
            x167 * x472 * x9,
            x163 * x475 * x479,
            x173 * x469 * x8,
            x169 * x472 * x8,
            x167 * x475 * x8,
            x163 * x478,
            x124 * x480 * x50,
            x158 * x48 * x480,
            x124 * x48 * x481,
            x127 * x480 * x9,
            x479 * x481 * x98,
            x4 * x482,
            x152 * x480 * x8,
            x127 * x481 * x8,
            x482 * x98,
            x315
            * (
                x174 * x477
                + x3 * (2 * x366 + 2 * x369 + 3 * x465 + 3 * x466 + 3 * x473 + 3 * x474)
            ),
        ]
    )


def quadrupole3d_44(a, A, b, B, C):
    """Cartesian 3D (gg) quadrupole moment integrals.
    The origin is at C.

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
    x9 = x4 * x8
    x10 = -x1 - C[0]
    x11 = x10 * x8
    x12 = x3 * (x11 + x9)
    x13 = x3 * x8
    x14 = x10 * x9
    x15 = x13 + x14
    x16 = x15 * x4
    x17 = x12 + x16
    x18 = x17 * x4
    x19 = x10 * x13
    x20 = 2 * x19
    x21 = x10 ** 2 * x8
    x22 = x13 + x21
    x23 = x22 * x4
    x24 = x20 + x23
    x25 = x24 * x4
    x26 = x4 ** 2 * x8
    x27 = 3 * x13
    x28 = 2 * x14 + x27
    x29 = x3 * (x26 + x28)
    x30 = 2 * x29
    x31 = x3 * (x21 + x28)
    x32 = x30 + 3 * x31
    x33 = x3 * (2 * x18 + 3 * x25 + x32)
    x34 = x25 + x31
    x35 = x34 * x4
    x36 = 2 * x12
    x37 = 4 * x19
    x38 = x36 + x37
    x39 = x3 * (2 * x16 + 2 * x23 + x38)
    x40 = x35 + x39
    x41 = x4 * x40
    x42 = x33 + x41
    x43 = x2 * x42
    x44 = x18 + x29
    x45 = x4 * x44
    x46 = 3 * x12
    x47 = x13 * x4
    x48 = 2 * x47
    x49 = x13 + x26
    x50 = x4 * x49
    x51 = x48 + x50
    x52 = x3 * (3 * x16 + x46 + x51)
    x53 = 2 * x52
    x54 = 4 * x39 + x53
    x55 = x3 * (4 * x35 + 2 * x45 + x54)
    x56 = x43 + x55
    x57 = x2 * x56
    x58 = x2 * x40
    x59 = 4 * x29
    x60 = x3 * (3 * x26 + x27)
    x61 = x4 * x51
    x62 = x60 + x61
    x63 = x3 * (4 * x18 + x59 + x62)
    x64 = x45 + x52
    x65 = x2 * x64
    x66 = 2 * x63 + 2 * x65
    x67 = x3 * (5 * x33 + x41 + 4 * x58 + x66)
    x68 = x57 + x67
    x69 = x63 + x65
    x70 = x2 * x69
    x71 = x33 + x58
    x72 = x2 * x71
    x73 = x2 * x44
    x74 = 8 * x47
    x75 = x3 * (4 * x50 + x74)
    x76 = x2 * x62
    x77 = x75 + x76
    x78 = x3 * (x45 + 5 * x52 + 4 * x73 + x77)
    x79 = 2 * x73
    x80 = x2 * x34
    x81 = x3 * (x35 + x54 + x79 + 3 * x80)
    x82 = x2 * x68 + x3 * (2 * x43 + 2 * x55 + 2 * x70 + 4 * x72 + 2 * x78 + 4 * x81)
    x83 = x70 + x78
    x84 = x2 * x83
    x85 = x72 + x81
    x86 = x2 * x85
    x87 = x17 * x2
    x88 = x2 * x51
    x89 = x60 + x88
    x90 = x3 * (x18 + x59 + 3 * x87 + x89)
    x91 = x52 + x73
    x92 = x2 * x91
    x93 = x3 * (5 * x60 + x61 + 4 * x88)
    x94 = x2 * x77
    x95 = x93 + x94
    x96 = x3 * (x66 + 4 * x90 + 4 * x92 + x95)
    x97 = 2 * x87
    x98 = x2 * x24
    x99 = 2 * x98
    x100 = x3 * (x25 + x32 + x97 + x99)
    x101 = x39 + x80
    x102 = x101 * x2
    x103 = 3 * x100 + 3 * x102
    x104 = x3 * (x103 + 2 * x33 + 2 * x58 + 2 * x90 + 2 * x92)
    x105 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x106 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x107 = numpy.pi * x0 * x106
    x108 = x105 * x107
    x109 = -x0 * (a * A[1] + b * B[1])
    x110 = -x109 - B[1]
    x111 = x104 + x86
    x112 = x100 + x102
    x113 = x112 * x2
    x114 = x90 + x92
    x115 = x114 * x2
    x116 = x2 * x49
    x117 = x3 * (3 * x116 + x50 + x74)
    x118 = x2 * x89
    x119 = x117 + x118
    x120 = x15 * x2
    x121 = 2 * x120
    x122 = x116 + x48
    x123 = x3 * (x121 + x122 + x16 + x46)
    x124 = x29 + x87
    x125 = x124 * x2
    x126 = 3 * x123 + 3 * x125
    x127 = x3 * (x119 + x126 + x53 + x79)
    x128 = x31 + x98
    x129 = x128 * x2
    x130 = x2 * x22
    x131 = x3 * (x121 + x130 + x23 + x38)
    x132 = x3 * (2 * x123 + 2 * x125 + 2 * x129 + 2 * x131 + 2 * x39 + 2 * x80)
    x133 = x108 * (
        x111 * x2 + x3 * (3 * x113 + 2 * x115 + 2 * x127 + 3 * x132 + 3 * x72 + 3 * x81)
    )
    x134 = -x0 * (a * A[2] + b * B[2])
    x135 = -x134 - B[2]
    x136 = x105 * x7
    x137 = x136 * x3
    x138 = x110 ** 2 * x136
    x139 = x137 + x138
    x140 = x113 + x132
    x141 = x2 * x9
    x142 = 2 * x141 + x27
    x143 = x3 * (x142 + x26)
    x144 = x122 * x2
    x145 = x143 + x144
    x146 = x11 * x2
    x147 = x3 * (x14 + x141 + x146 + x27)
    x148 = x12 + x120
    x149 = x148 * x2
    x150 = 2 * x147 + 2 * x149
    x151 = x3 * (x145 + x150 + x30 + x97)
    x152 = 2 * x146 + x27
    x153 = x3 * (x152 + x21)
    x154 = x130 + x20
    x155 = x154 * x2
    x156 = x153 + x155
    x157 = x3 * (x150 + x156 + 2 * x31 + x99)
    x158 = x123 + x125
    x159 = x158 * x2
    x160 = x129 + x131
    x161 = x160 * x2
    x162 = x140 * x2 + x3 * (x103 + 2 * x151 + 2 * x157 + 2 * x159 + 2 * x161)
    x163 = x106 * x7
    x164 = x108 * x135
    x165 = x163 * x3
    x166 = x135 ** 2 * x163
    x167 = x165 + x166
    x168 = x110 * x137
    x169 = 2 * x168
    x170 = x110 * x139
    x171 = x169 + x170
    x172 = x157 + x161
    x173 = x2 * x8
    x174 = x3 * (x11 + x173)
    x175 = x13 + x146
    x176 = x175 * x2
    x177 = x156 * x2 + x3 * (2 * x130 + 2 * x174 + 2 * x176 + x37)
    x178 = x174 + x176
    x179 = x3 * (x173 + x9)
    x180 = x13 + x141
    x181 = x180 * x2
    x182 = x179 + x181
    x183 = x3 * (x121 + x178 + x182 + x36)
    x184 = x147 + x149
    x185 = x184 * x2
    x186 = 2 * x183 + 2 * x185
    x187 = x172 * x2 + x3 * (3 * x129 + 3 * x131 + x177 + x186)
    x188 = x135 * x163
    x189 = x110 * x136
    x190 = x135 * x165
    x191 = 2 * x190
    x192 = x135 * x167
    x193 = x191 + x192
    x194 = 3 * x137
    x195 = x3 * (3 * x138 + x194)
    x196 = x110 * x171
    x197 = x195 + x196
    x198 = x2 ** 2 * x8
    x199 = x3 * (x152 + x198)
    x200 = x178 * x2
    x201 = x177 * x2 + x3 * (3 * x153 + 3 * x155 + 2 * x199 + 2 * x200)
    x202 = 3 * x165
    x203 = x3 * (3 * x166 + x202)
    x204 = x135 * x193
    x205 = x203 + x204
    x206 = -x109 - A[1]
    x207 = x108 * x82
    x208 = x136 * x206
    x209 = x110 * x208
    x210 = x137 + x209
    x211 = x139 * x206
    x212 = x169 + x211
    x213 = x171 * x206
    x214 = x195 + x213
    x215 = 8 * x168
    x216 = x3 * (4 * x170 + x215)
    x217 = x197 * x206
    x218 = x216 + x217
    x219 = -x134 - A[2]
    x220 = x108 * x219
    x221 = x163 * x219
    x222 = x135 * x221
    x223 = x165 + x222
    x224 = x167 * x219
    x225 = x191 + x224
    x226 = x193 * x219
    x227 = x203 + x226
    x228 = 8 * x190
    x229 = x3 * (4 * x192 + x228)
    x230 = x205 * x219
    x231 = x229 + x230
    x232 = x136 * x206 ** 2
    x233 = x137 + x232
    x234 = x3 * (x189 + x208)
    x235 = x206 * x210
    x236 = x234 + x235
    x237 = x194 + 2 * x209
    x238 = x3 * (x138 + x237)
    x239 = x206 * x212
    x240 = x238 + x239
    x241 = x3 * (x170 + 3 * x211 + x215)
    x242 = x206 * x214
    x243 = x241 + x242
    x244 = x3 * (5 * x195 + x196 + 4 * x213)
    x245 = x206 * x218
    x246 = x244 + x245
    x247 = x163 * x219 ** 2
    x248 = x165 + x247
    x249 = x3 * (x188 + x221)
    x250 = x219 * x223
    x251 = x249 + x250
    x252 = x202 + 2 * x222
    x253 = x3 * (x166 + x252)
    x254 = x219 * x225
    x255 = x253 + x254
    x256 = x3 * (x192 + 3 * x224 + x228)
    x257 = x219 * x227
    x258 = x256 + x257
    x259 = x3 * (5 * x203 + x204 + 4 * x226)
    x260 = x219 * x231
    x261 = x259 + x260
    x262 = 2 * x137
    x263 = x206 * x233 + x206 * x262
    x264 = x3 * (x232 + x237)
    x265 = x206 * x236
    x266 = x264 + x265
    x267 = x206 * x240
    x268 = x3 * (4 * x168 + 2 * x211 + 2 * x234 + 2 * x235)
    x269 = x267 + x268
    x270 = x206 * x243
    x271 = 3 * x238 + 3 * x239
    x272 = x3 * (2 * x195 + 2 * x213 + x271)
    x273 = x270 + x272
    x274 = x206 * x246 + x3 * (2 * x216 + 2 * x217 + 4 * x241 + 4 * x242)
    x275 = 2 * x165
    x276 = x219 * x248 + x219 * x275
    x277 = x3 * (x247 + x252)
    x278 = x219 * x251
    x279 = x277 + x278
    x280 = x219 * x255
    x281 = x3 * (4 * x190 + 2 * x224 + 2 * x249 + 2 * x250)
    x282 = x280 + x281
    x283 = x219 * x258
    x284 = 3 * x253 + 3 * x254
    x285 = x3 * (2 * x203 + 2 * x226 + x284)
    x286 = x283 + x285
    x287 = x219 * x261 + x3 * (2 * x229 + 2 * x230 + 4 * x256 + 4 * x257)
    x288 = x206 * x263 + x3 * (x194 + 3 * x232)
    x289 = x206 * x266 + x3 * (3 * x234 + 3 * x235 + x263)
    x290 = x206 * x269 + x3 * (2 * x264 + 2 * x265 + x271)
    x291 = x206 * x273 + x3 * (3 * x241 + 3 * x242 + 3 * x267 + 3 * x268)
    x292 = x206 * x274 + x3 * (3 * x244 + 3 * x245 + 4 * x270 + 4 * x272)
    x293 = x219 * x276 + x3 * (x202 + 3 * x247)
    x294 = x219 * x279 + x3 * (3 * x249 + 3 * x250 + x276)
    x295 = x219 * x282 + x3 * (2 * x277 + 2 * x278 + x284)
    x296 = x219 * x286 + x3 * (3 * x256 + 3 * x257 + 3 * x280 + 3 * x281)
    x297 = x219 * x287 + x3 * (3 * x259 + 3 * x260 + 4 * x283 + 4 * x285)
    x298 = -x109 - C[1]
    x299 = x84 + x96
    x300 = x2 * x95 + x3 * (4 * x117 + 4 * x118 + 2 * x75 + 2 * x76)
    x301 = x108 * (x2 * x299 + x3 * (4 * x115 + 4 * x127 + x300 + 3 * x70 + 3 * x78))
    x302 = x189 * x298
    x303 = x137 + x302
    x304 = x115 + x127
    x305 = x119 * x2
    x306 = 3 * x143 + 3 * x144
    x307 = x3 * (x306 + 2 * x60 + 2 * x88)
    x308 = x305 + x307
    x309 = x2 * x304 + x3 * (3 * x151 + 3 * x159 + x308 + 3 * x90 + 3 * x92)
    x310 = x136 * x298
    x311 = x3 * (x189 + x310)
    x312 = x110 * x303
    x313 = x311 + x312
    x314 = x151 + x159
    x315 = x145 * x2
    x316 = x3 * (2 * x116 + 2 * x179 + 2 * x181 + 4 * x47)
    x317 = x315 + x316
    x318 = x2 * x314 + x3 * (x126 + x186 + x317)
    x319 = x194 + 2 * x302
    x320 = x3 * (x138 + x319)
    x321 = x110 * x313
    x322 = x320 + x321
    x323 = x183 + x185
    x324 = x199 + x200
    x325 = x3 * (x142 + x198)
    x326 = x182 * x2
    x327 = x325 + x326
    x328 = x2 * x323 + x3 * (3 * x147 + 3 * x149 + x324 + x327)
    x329 = 3 * x311
    x330 = x3 * (x171 + 3 * x312 + x329)
    x331 = x110 * x322
    x332 = x330 + x331
    x333 = x13 + x198
    x334 = 2 * x13 * x2 + x2 * x333
    x335 = x2 * x324 + x3 * (3 * x174 + 3 * x176 + x334)
    x336 = x208 * x298
    x337 = x137 + x336
    x338 = x206 * x303
    x339 = x311 + x338
    x340 = x206 * x313
    x341 = x320 + x340
    x342 = x206 * x322
    x343 = x330 + x342
    x344 = 4 * x320
    x345 = x3 * (x197 + 4 * x321 + x344)
    x346 = x206 * x332
    x347 = x345 + x346
    x348 = x3 * (x208 + x310)
    x349 = x206 * x337
    x350 = x348 + x349
    x351 = x3 * (x194 + x209 + x302 + x336)
    x352 = x206 * x339
    x353 = x351 + x352
    x354 = x206 * x341
    x355 = 2 * x338
    x356 = x3 * (x212 + x312 + x329 + x355)
    x357 = x354 + x356
    x358 = x206 * x343
    x359 = x3 * (x214 + x321 + 3 * x340 + x344)
    x360 = x358 + x359
    x361 = x206 * x347
    x362 = x3 * (x218 + 5 * x330 + x331 + 4 * x342)
    x363 = x361 + x362
    x364 = x194 + 2 * x336
    x365 = x3 * (x232 + x364)
    x366 = x206 * x350
    x367 = x365 + x366
    x368 = x206 * x353
    x369 = 2 * x311
    x370 = x3 * (x236 + x350 + x355 + x369)
    x371 = x368 + x370
    x372 = x206 * x357
    x373 = 2 * x320
    x374 = 2 * x340
    x375 = 2 * x351 + 2 * x352
    x376 = x3 * (x240 + x373 + x374 + x375)
    x377 = x372 + x376
    x378 = x206 * x360
    x379 = 2 * x330
    x380 = 2 * x342
    x381 = 3 * x354 + 3 * x356
    x382 = x3 * (x243 + x379 + x380 + x381)
    x383 = x378 + x382
    x384 = x206 * x363
    x385 = 2 * x345 + 2 * x346
    x386 = x3 * (x246 + 4 * x358 + 4 * x359 + x385)
    x387 = x384 + x386
    x388 = x206 * x367 + x3 * (x263 + 3 * x348 + 3 * x349)
    x389 = x206 * x371 + x3 * (x266 + 3 * x351 + 3 * x352 + x367)
    x390 = 2 * x368 + 2 * x370
    x391 = x206 * x377 + x3 * (x269 + x381 + x390)
    x392 = x206 * x383 + x3 * (x273 + 3 * x358 + 3 * x359 + 3 * x372 + 3 * x376)
    x393 = x107 * x6
    x394 = x393 * (x206 * x387 + x3 * (x274 + 3 * x361 + 3 * x362 + 4 * x378 + 4 * x382))
    x395 = x10 * x393
    x396 = numpy.pi * x0 * x105 * x6
    x397 = x10 * x396
    x398 = -x134 - C[2]
    x399 = x108 * x398
    x400 = x188 * x398
    x401 = x165 + x400
    x402 = x163 * x398
    x403 = x3 * (x188 + x402)
    x404 = x135 * x401
    x405 = x403 + x404
    x406 = x202 + 2 * x400
    x407 = x3 * (x166 + x406)
    x408 = x135 * x405
    x409 = x407 + x408
    x410 = 3 * x403
    x411 = x3 * (x193 + 3 * x404 + x410)
    x412 = x135 * x409
    x413 = x411 + x412
    x414 = x221 * x398
    x415 = x165 + x414
    x416 = x219 * x401
    x417 = x403 + x416
    x418 = x219 * x405
    x419 = x407 + x418
    x420 = x219 * x409
    x421 = x411 + x420
    x422 = 4 * x407
    x423 = x3 * (x205 + 4 * x408 + x422)
    x424 = x219 * x413
    x425 = x423 + x424
    x426 = x3 * (x221 + x402)
    x427 = x219 * x415
    x428 = x426 + x427
    x429 = x3 * (x202 + x222 + x400 + x414)
    x430 = x219 * x417
    x431 = x429 + x430
    x432 = x219 * x419
    x433 = 2 * x416
    x434 = x3 * (x225 + x404 + x410 + x433)
    x435 = x432 + x434
    x436 = x219 * x421
    x437 = x3 * (x227 + x408 + 3 * x418 + x422)
    x438 = x436 + x437
    x439 = x219 * x425
    x440 = x3 * (x231 + 5 * x411 + x412 + 4 * x420)
    x441 = x439 + x440
    x442 = x202 + 2 * x414
    x443 = x3 * (x247 + x442)
    x444 = x219 * x428
    x445 = x443 + x444
    x446 = x219 * x431
    x447 = 2 * x403
    x448 = x3 * (x251 + x428 + x433 + x447)
    x449 = x446 + x448
    x450 = x219 * x435
    x451 = 2 * x407
    x452 = 2 * x418
    x453 = 2 * x429 + 2 * x430
    x454 = x3 * (x255 + x451 + x452 + x453)
    x455 = x450 + x454
    x456 = x219 * x438
    x457 = 2 * x411
    x458 = 2 * x420
    x459 = 3 * x432 + 3 * x434
    x460 = x3 * (x258 + x457 + x458 + x459)
    x461 = x456 + x460
    x462 = x219 * x441
    x463 = 2 * x423 + 2 * x424
    x464 = x3 * (x261 + 4 * x436 + 4 * x437 + x463)
    x465 = x462 + x464
    x466 = x219 * x445 + x3 * (x276 + 3 * x426 + 3 * x427)
    x467 = x219 * x449 + x3 * (x279 + 3 * x429 + 3 * x430 + x445)
    x468 = 2 * x446 + 2 * x448
    x469 = x219 * x455 + x3 * (x282 + x459 + x468)
    x470 = x219 * x461 + x3 * (x286 + 3 * x436 + 3 * x437 + 3 * x450 + 3 * x454)
    x471 = x396 * (x219 * x465 + x3 * (x287 + 3 * x439 + 3 * x440 + 4 * x456 + 4 * x460))
    x472 = x136 * x298 ** 2
    x473 = x137 + x472
    x474 = x2 * x300 + x3 * (4 * x305 + 4 * x307 + 3 * x93 + 3 * x94)
    x475 = x262 * x298
    x476 = x110 * x473
    x477 = x475 + x476
    x478 = x2 * x308 + x3 * (3 * x117 + 3 * x118 + 3 * x315 + 3 * x316)
    x479 = x3 * (x319 + x472)
    x480 = x110 * x477
    x481 = x479 + x480
    x482 = x2 * x317 + x3 * (x306 + 2 * x325 + 2 * x326)
    x483 = x110 * x481
    x484 = 4 * x137 * x298
    x485 = x369 + x484
    x486 = x3 * (2 * x312 + 2 * x476 + x485)
    x487 = x483 + x486
    x488 = x2 * x327 + x3 * (3 * x179 + 3 * x181 + x334)
    x489 = x373 + 3 * x479
    x490 = x3 * (2 * x321 + 3 * x480 + x489)
    x491 = x110 * x487
    x492 = x490 + x491
    x493 = x2 * x334 + x3 * (3 * x198 + x27)
    x494 = x206 * x473
    x495 = x475 + x494
    x496 = x206 * x477
    x497 = x479 + x496
    x498 = x206 * x481
    x499 = x486 + x498
    x500 = x206 * x487
    x501 = x490 + x500
    x502 = x206 * x492
    x503 = x379 + 4 * x486
    x504 = x3 * (2 * x331 + 4 * x483 + x503)
    x505 = x502 + x504
    x506 = x3 * (x364 + x472)
    x507 = x206 * x495
    x508 = x506 + x507
    x509 = x206 * x497
    x510 = x3 * (x355 + x476 + x485 + x494)
    x511 = x509 + x510
    x512 = x206 * x499
    x513 = 2 * x496
    x514 = x3 * (x374 + x480 + x489 + x513)
    x515 = x512 + x514
    x516 = x206 * x501
    x517 = x3 * (x380 + x483 + 3 * x498 + x503)
    x518 = x516 + x517
    x519 = x206 * x505
    x520 = x3 * (x385 + 5 * x490 + x491 + 4 * x500)
    x521 = x519 + x520
    x522 = x206 * x508 + x3 * (2 * x348 + 2 * x349 + x484 + 2 * x494)
    x523 = x206 * x511
    x524 = x3 * (x375 + 2 * x479 + x508 + x513)
    x525 = x523 + x524
    x526 = x206 * x515
    x527 = x3 * (2 * x354 + 2 * x356 + 2 * x486 + 2 * x498 + 2 * x509 + 2 * x510)
    x528 = x526 + x527
    x529 = x206 * x518
    x530 = 3 * x512 + 3 * x514
    x531 = x3 * (2 * x358 + 2 * x359 + 2 * x490 + 2 * x500 + x530)
    x532 = x529 + x531
    x533 = x206 * x521 + x3 * (
        2 * x361 + 2 * x362 + 2 * x502 + 2 * x504 + 4 * x516 + 4 * x517
    )
    x534 = x393 * x533
    x535 = x2 * x393
    x536 = x206 * x522 + x3 * (2 * x365 + 2 * x366 + 3 * x506 + 3 * x507)
    x537 = x206 * x525 + x3 * (x390 + 3 * x509 + 3 * x510 + x522)
    x538 = x206 * x528 + x3 * (2 * x372 + 2 * x376 + 2 * x523 + 2 * x524 + x530)
    x539 = x393 * (
        x206 * x532
        + x3 * (2 * x378 + 2 * x382 + 3 * x516 + 3 * x517 + 3 * x526 + 3 * x527)
    )
    x540 = x393 * x4
    x541 = x298 * x396
    x542 = x163 * x398 ** 2
    x543 = x165 + x542
    x544 = x275 * x398
    x545 = x135 * x543
    x546 = x544 + x545
    x547 = x3 * (x406 + x542)
    x548 = x135 * x546
    x549 = x547 + x548
    x550 = x135 * x549
    x551 = 4 * x165 * x398
    x552 = x447 + x551
    x553 = x3 * (2 * x404 + 2 * x545 + x552)
    x554 = x550 + x553
    x555 = x451 + 3 * x547
    x556 = x3 * (2 * x408 + 3 * x548 + x555)
    x557 = x135 * x554
    x558 = x556 + x557
    x559 = x219 * x543
    x560 = x544 + x559
    x561 = x219 * x546
    x562 = x547 + x561
    x563 = x219 * x549
    x564 = x553 + x563
    x565 = x219 * x554
    x566 = x556 + x565
    x567 = x219 * x558
    x568 = x457 + 4 * x553
    x569 = x3 * (2 * x412 + 4 * x550 + x568)
    x570 = x567 + x569
    x571 = x3 * (x442 + x542)
    x572 = x219 * x560
    x573 = x571 + x572
    x574 = x219 * x562
    x575 = x3 * (x433 + x545 + x552 + x559)
    x576 = x574 + x575
    x577 = x219 * x564
    x578 = 2 * x561
    x579 = x3 * (x452 + x548 + x555 + x578)
    x580 = x577 + x579
    x581 = x219 * x566
    x582 = x3 * (x458 + x550 + 3 * x563 + x568)
    x583 = x581 + x582
    x584 = x219 * x570
    x585 = x3 * (x463 + 5 * x556 + x557 + 4 * x565)
    x586 = x584 + x585
    x587 = x2 * x396
    x588 = x219 * x573 + x3 * (2 * x426 + 2 * x427 + x551 + 2 * x559)
    x589 = x219 * x576
    x590 = x3 * (x453 + 2 * x547 + x573 + x578)
    x591 = x589 + x590
    x592 = x219 * x580
    x593 = x3 * (2 * x432 + 2 * x434 + 2 * x553 + 2 * x563 + 2 * x574 + 2 * x575)
    x594 = x592 + x593
    x595 = x219 * x583
    x596 = 3 * x577 + 3 * x579
    x597 = x3 * (2 * x436 + 2 * x437 + 2 * x556 + 2 * x565 + x596)
    x598 = x595 + x597
    x599 = x219 * x586 + x3 * (
        2 * x439 + 2 * x440 + 2 * x567 + 2 * x569 + 4 * x581 + 4 * x582
    )
    x600 = x396 * x599
    x601 = x396 * x4
    x602 = x219 * x588 + x3 * (2 * x443 + 2 * x444 + 3 * x571 + 3 * x572)
    x603 = x219 * x591 + x3 * (x468 + 3 * x574 + 3 * x575 + x588)
    x604 = x219 * x594 + x3 * (2 * x450 + 2 * x454 + 2 * x589 + 2 * x590 + x596)
    x605 = x396 * (
        x219 * x598
        + x3 * (2 * x456 + 2 * x460 + 3 * x581 + 3 * x582 + 3 * x592 + 3 * x593)
    )

    # 1350 item(s)
    return numpy.array(
        [
            x108
            * (
                x2 * x82
                + x3 * (4 * x104 + 3 * x57 + 3 * x67 + 2 * x84 + 4 * x86 + 2 * x96)
            ),
            x110 * x133,
            x133 * x135,
            x139 * x162 * x163,
            x110 * x162 * x164,
            x136 * x162 * x167,
            x163 * x171 * x187,
            x139 * x187 * x188,
            x167 * x187 * x189,
            x136 * x187 * x193,
            x163 * x197 * x201,
            x171 * x188 * x201,
            x139 * x167 * x201,
            x189 * x193 * x201,
            x136 * x201 * x205,
            x206 * x207,
            x111 * x163 * x210,
            x111 * x164 * x206,
            x140 * x163 * x212,
            x140 * x188 * x210,
            x140 * x167 * x208,
            x163 * x172 * x214,
            x172 * x188 * x212,
            x167 * x172 * x210,
            x172 * x193 * x208,
            x163 * x177 * x218,
            x177 * x188 * x214,
            x167 * x177 * x212,
            x177 * x193 * x210,
            x177 * x205 * x208,
            x207 * x219,
            x110 * x111 * x220,
            x111 * x136 * x223,
            x139 * x140 * x221,
            x140 * x189 * x223,
            x136 * x140 * x225,
            x171 * x172 * x221,
            x139 * x172 * x223,
            x172 * x189 * x225,
            x136 * x172 * x227,
            x177 * x197 * x221,
            x171 * x177 * x223,
            x139 * x177 * x225,
            x177 * x189 * x227,
            x136 * x177 * x231,
            x163 * x233 * x68,
            x163 * x236 * x85,
            x188 * x233 * x85,
            x112 * x163 * x240,
            x112 * x188 * x236,
            x112 * x167 * x233,
            x160 * x163 * x243,
            x160 * x188 * x240,
            x160 * x167 * x236,
            x160 * x193 * x233,
            x156 * x163 * x246,
            x156 * x188 * x243,
            x156 * x167 * x240,
            x156 * x193 * x236,
            x156 * x205 * x233,
            x206 * x220 * x68,
            x210 * x221 * x85,
            x208 * x223 * x85,
            x112 * x212 * x221,
            x112 * x210 * x223,
            x112 * x208 * x225,
            x160 * x214 * x221,
            x160 * x212 * x223,
            x160 * x210 * x225,
            x160 * x208 * x227,
            x156 * x218 * x221,
            x156 * x214 * x223,
            x156 * x212 * x225,
            x156 * x210 * x227,
            x156 * x208 * x231,
            x136 * x248 * x68,
            x189 * x248 * x85,
            x136 * x251 * x85,
            x112 * x139 * x248,
            x112 * x189 * x251,
            x112 * x136 * x255,
            x160 * x171 * x248,
            x139 * x160 * x251,
            x160 * x189 * x255,
            x136 * x160 * x258,
            x156 * x197 * x248,
            x156 * x171 * x251,
            x139 * x156 * x255,
            x156 * x189 * x258,
            x136 * x156 * x261,
            x163 * x263 * x56,
            x163 * x266 * x71,
            x188 * x263 * x71,
            x101 * x163 * x269,
            x101 * x188 * x266,
            x101 * x167 * x263,
            x128 * x163 * x273,
            x128 * x188 * x269,
            x128 * x167 * x266,
            x128 * x193 * x263,
            x154 * x163 * x274,
            x154 * x188 * x273,
            x154 * x167 * x269,
            x154 * x193 * x266,
            x154 * x205 * x263,
            x221 * x233 * x56,
            x221 * x236 * x71,
            x223 * x233 * x71,
            x101 * x221 * x240,
            x101 * x223 * x236,
            x101 * x225 * x233,
            x128 * x221 * x243,
            x128 * x223 * x240,
            x128 * x225 * x236,
            x128 * x227 * x233,
            x154 * x221 * x246,
            x154 * x223 * x243,
            x154 * x225 * x240,
            x154 * x227 * x236,
            x154 * x231 * x233,
            x208 * x248 * x56,
            x210 * x248 * x71,
            x208 * x251 * x71,
            x101 * x212 * x248,
            x101 * x210 * x251,
            x101 * x208 * x255,
            x128 * x214 * x248,
            x128 * x212 * x251,
            x128 * x210 * x255,
            x128 * x208 * x258,
            x154 * x218 * x248,
            x154 * x214 * x251,
            x154 * x212 * x255,
            x154 * x210 * x258,
            x154 * x208 * x261,
            x136 * x276 * x56,
            x189 * x276 * x71,
            x136 * x279 * x71,
            x101 * x139 * x276,
            x101 * x189 * x279,
            x101 * x136 * x282,
            x128 * x171 * x276,
            x128 * x139 * x279,
            x128 * x189 * x282,
            x128 * x136 * x286,
            x154 * x197 * x276,
            x154 * x171 * x279,
            x139 * x154 * x282,
            x154 * x189 * x286,
            x136 * x154 * x287,
            x163 * x288 * x42,
            x163 * x289 * x40,
            x188 * x288 * x40,
            x163 * x290 * x34,
            x188 * x289 * x34,
            x167 * x288 * x34,
            x163 * x24 * x291,
            x188 * x24 * x290,
            x167 * x24 * x289,
            x193 * x24 * x288,
            x163 * x22 * x292,
            x188 * x22 * x291,
            x167 * x22 * x290,
            x193 * x22 * x289,
            x205 * x22 * x288,
            x221 * x263 * x42,
            x221 * x266 * x40,
            x223 * x263 * x40,
            x221 * x269 * x34,
            x223 * x266 * x34,
            x225 * x263 * x34,
            x221 * x24 * x273,
            x223 * x24 * x269,
            x225 * x24 * x266,
            x227 * x24 * x263,
            x22 * x221 * x274,
            x22 * x223 * x273,
            x22 * x225 * x269,
            x22 * x227 * x266,
            x22 * x231 * x263,
            x233 * x248 * x42,
            x236 * x248 * x40,
            x233 * x251 * x40,
            x240 * x248 * x34,
            x236 * x251 * x34,
            x233 * x255 * x34,
            x24 * x243 * x248,
            x24 * x240 * x251,
            x236 * x24 * x255,
            x233 * x24 * x258,
            x22 * x246 * x248,
            x22 * x243 * x251,
            x22 * x240 * x255,
            x22 * x236 * x258,
            x22 * x233 * x261,
            x208 * x276 * x42,
            x210 * x276 * x40,
            x208 * x279 * x40,
            x212 * x276 * x34,
            x210 * x279 * x34,
            x208 * x282 * x34,
            x214 * x24 * x276,
            x212 * x24 * x279,
            x210 * x24 * x282,
            x208 * x24 * x286,
            x218 * x22 * x276,
            x214 * x22 * x279,
            x212 * x22 * x282,
            x210 * x22 * x286,
            x208 * x22 * x287,
            x136 * x293 * x42,
            x189 * x293 * x40,
            x136 * x294 * x40,
            x139 * x293 * x34,
            x189 * x294 * x34,
            x136 * x295 * x34,
            x171 * x24 * x293,
            x139 * x24 * x294,
            x189 * x24 * x295,
            x136 * x24 * x296,
            x197 * x22 * x293,
            x171 * x22 * x294,
            x139 * x22 * x295,
            x189 * x22 * x296,
            x136 * x22 * x297,
            x298 * x301,
            x163 * x303 * x309,
            x164 * x298 * x309,
            x163 * x313 * x318,
            x188 * x303 * x318,
            x167 * x310 * x318,
            x163 * x322 * x328,
            x188 * x313 * x328,
            x167 * x303 * x328,
            x193 * x310 * x328,
            x163 * x332 * x335,
            x188 * x322 * x335,
            x167 * x313 * x335,
            x193 * x303 * x335,
            x205 * x310 * x335,
            x163 * x299 * x337,
            x163 * x304 * x339,
            x188 * x304 * x337,
            x163 * x314 * x341,
            x188 * x314 * x339,
            x167 * x314 * x337,
            x163 * x323 * x343,
            x188 * x323 * x341,
            x167 * x323 * x339,
            x193 * x323 * x337,
            x163 * x324 * x347,
            x188 * x324 * x343,
            x167 * x324 * x341,
            x193 * x324 * x339,
            x205 * x324 * x337,
            x220 * x298 * x299,
            x221 * x303 * x304,
            x223 * x304 * x310,
            x221 * x313 * x314,
            x223 * x303 * x314,
            x225 * x310 * x314,
            x221 * x322 * x323,
            x223 * x313 * x323,
            x225 * x303 * x323,
            x227 * x310 * x323,
            x221 * x324 * x332,
            x223 * x322 * x324,
            x225 * x313 * x324,
            x227 * x303 * x324,
            x231 * x310 * x324,
            x163 * x350 * x83,
            x114 * x163 * x353,
            x114 * x188 * x350,
            x158 * x163 * x357,
            x158 * x188 * x353,
            x158 * x167 * x350,
            x163 * x184 * x360,
            x184 * x188 * x357,
            x167 * x184 * x353,
            x184 * x193 * x350,
            x163 * x178 * x363,
            x178 * x188 * x360,
            x167 * x178 * x357,
            x178 * x193 * x353,
            x178 * x205 * x350,
            x221 * x337 * x83,
            x114 * x221 * x339,
            x114 * x223 * x337,
            x158 * x221 * x341,
            x158 * x223 * x339,
            x158 * x225 * x337,
            x184 * x221 * x343,
            x184 * x223 * x341,
            x184 * x225 * x339,
            x184 * x227 * x337,
            x178 * x221 * x347,
            x178 * x223 * x343,
            x178 * x225 * x341,
            x178 * x227 * x339,
            x178 * x231 * x337,
            x248 * x310 * x83,
            x114 * x248 * x303,
            x114 * x251 * x310,
            x158 * x248 * x313,
            x158 * x251 * x303,
            x158 * x255 * x310,
            x184 * x248 * x322,
            x184 * x251 * x313,
            x184 * x255 * x303,
            x184 * x258 * x310,
            x178 * x248 * x332,
            x178 * x251 * x322,
            x178 * x255 * x313,
            x178 * x258 * x303,
            x178 * x261 * x310,
            x163 * x367 * x69,
            x163 * x371 * x91,
            x188 * x367 * x91,
            x124 * x163 * x377,
            x124 * x188 * x371,
            x124 * x167 * x367,
            x148 * x163 * x383,
            x148 * x188 * x377,
            x148 * x167 * x371,
            x148 * x193 * x367,
            x163 * x175 * x387,
            x175 * x188 * x383,
            x167 * x175 * x377,
            x175 * x193 * x371,
            x175 * x205 * x367,
            x221 * x350 * x69,
            x221 * x353 * x91,
            x223 * x350 * x91,
            x124 * x221 * x357,
            x124 * x223 * x353,
            x124 * x225 * x350,
            x148 * x221 * x360,
            x148 * x223 * x357,
            x148 * x225 * x353,
            x148 * x227 * x350,
            x175 * x221 * x363,
            x175 * x223 * x360,
            x175 * x225 * x357,
            x175 * x227 * x353,
            x175 * x231 * x350,
            x248 * x337 * x69,
            x248 * x339 * x91,
            x251 * x337 * x91,
            x124 * x248 * x341,
            x124 * x251 * x339,
            x124 * x255 * x337,
            x148 * x248 * x343,
            x148 * x251 * x341,
            x148 * x255 * x339,
            x148 * x258 * x337,
            x175 * x248 * x347,
            x175 * x251 * x343,
            x175 * x255 * x341,
            x175 * x258 * x339,
            x175 * x261 * x337,
            x276 * x310 * x69,
            x276 * x303 * x91,
            x279 * x310 * x91,
            x124 * x276 * x313,
            x124 * x279 * x303,
            x124 * x282 * x310,
            x148 * x276 * x322,
            x148 * x279 * x313,
            x148 * x282 * x303,
            x148 * x286 * x310,
            x175 * x276 * x332,
            x175 * x279 * x322,
            x175 * x282 * x313,
            x175 * x286 * x303,
            x175 * x287 * x310,
            x163 * x388 * x64,
            x163 * x389 * x44,
            x188 * x388 * x44,
            x163 * x17 * x391,
            x17 * x188 * x389,
            x167 * x17 * x388,
            x15 * x163 * x392,
            x15 * x188 * x391,
            x15 * x167 * x389,
            x15 * x193 * x388,
            x10 * x394,
            x135 * x392 * x395,
            x11 * x167 * x391,
            x11 * x193 * x389,
            x11 * x205 * x388,
            x221 * x367 * x64,
            x221 * x371 * x44,
            x223 * x367 * x44,
            x17 * x221 * x377,
            x17 * x223 * x371,
            x17 * x225 * x367,
            x15 * x221 * x383,
            x15 * x223 * x377,
            x15 * x225 * x371,
            x15 * x227 * x367,
            x219 * x387 * x395,
            x11 * x223 * x383,
            x11 * x225 * x377,
            x11 * x227 * x371,
            x11 * x231 * x367,
            x248 * x350 * x64,
            x248 * x353 * x44,
            x251 * x350 * x44,
            x17 * x248 * x357,
            x17 * x251 * x353,
            x17 * x255 * x350,
            x15 * x248 * x360,
            x15 * x251 * x357,
            x15 * x255 * x353,
            x15 * x258 * x350,
            x11 * x248 * x363,
            x11 * x251 * x360,
            x11 * x255 * x357,
            x11 * x258 * x353,
            x11 * x261 * x350,
            x276 * x337 * x64,
            x276 * x339 * x44,
            x279 * x337 * x44,
            x17 * x276 * x341,
            x17 * x279 * x339,
            x17 * x282 * x337,
            x15 * x276 * x343,
            x15 * x279 * x341,
            x15 * x282 * x339,
            x15 * x286 * x337,
            x11 * x276 * x347,
            x11 * x279 * x343,
            x11 * x282 * x341,
            x11 * x286 * x339,
            x11 * x287 * x337,
            x293 * x310 * x64,
            x293 * x303 * x44,
            x294 * x310 * x44,
            x17 * x293 * x313,
            x17 * x294 * x303,
            x17 * x295 * x310,
            x15 * x293 * x322,
            x15 * x294 * x313,
            x15 * x295 * x303,
            x15 * x296 * x310,
            x11 * x293 * x332,
            x11 * x294 * x322,
            x11 * x295 * x313,
            x11 * x296 * x303,
            x297 * x298 * x397,
            x301 * x398,
            x110 * x309 * x399,
            x136 * x309 * x401,
            x139 * x318 * x402,
            x189 * x318 * x401,
            x136 * x318 * x405,
            x171 * x328 * x402,
            x139 * x328 * x401,
            x189 * x328 * x405,
            x136 * x328 * x409,
            x197 * x335 * x402,
            x171 * x335 * x401,
            x139 * x335 * x405,
            x189 * x335 * x409,
            x136 * x335 * x413,
            x206 * x299 * x399,
            x210 * x304 * x402,
            x208 * x304 * x401,
            x212 * x314 * x402,
            x210 * x314 * x401,
            x208 * x314 * x405,
            x214 * x323 * x402,
            x212 * x323 * x401,
            x210 * x323 * x405,
            x208 * x323 * x409,
            x218 * x324 * x402,
            x214 * x324 * x401,
            x212 * x324 * x405,
            x210 * x324 * x409,
            x208 * x324 * x413,
            x136 * x299 * x415,
            x189 * x304 * x415,
            x136 * x304 * x417,
            x139 * x314 * x415,
            x189 * x314 * x417,
            x136 * x314 * x419,
            x171 * x323 * x415,
            x139 * x323 * x417,
            x189 * x323 * x419,
            x136 * x323 * x421,
            x197 * x324 * x415,
            x171 * x324 * x417,
            x139 * x324 * x419,
            x189 * x324 * x421,
            x136 * x324 * x425,
            x233 * x402 * x83,
            x114 * x236 * x402,
            x114 * x233 * x401,
            x158 * x240 * x402,
            x158 * x236 * x401,
            x158 * x233 * x405,
            x184 * x243 * x402,
            x184 * x240 * x401,
            x184 * x236 * x405,
            x184 * x233 * x409,
            x178 * x246 * x402,
            x178 * x243 * x401,
            x178 * x240 * x405,
            x178 * x236 * x409,
            x178 * x233 * x413,
            x208 * x415 * x83,
            x114 * x210 * x415,
            x114 * x208 * x417,
            x158 * x212 * x415,
            x158 * x210 * x417,
            x158 * x208 * x419,
            x184 * x214 * x415,
            x184 * x212 * x417,
            x184 * x210 * x419,
            x184 * x208 * x421,
            x178 * x218 * x415,
            x178 * x214 * x417,
            x178 * x212 * x419,
            x178 * x210 * x421,
            x178 * x208 * x425,
            x136 * x428 * x83,
            x114 * x189 * x428,
            x114 * x136 * x431,
            x139 * x158 * x428,
            x158 * x189 * x431,
            x136 * x158 * x435,
            x171 * x184 * x428,
            x139 * x184 * x431,
            x184 * x189 * x435,
            x136 * x184 * x438,
            x178 * x197 * x428,
            x171 * x178 * x431,
            x139 * x178 * x435,
            x178 * x189 * x438,
            x136 * x178 * x441,
            x263 * x402 * x69,
            x266 * x402 * x91,
            x263 * x401 * x91,
            x124 * x269 * x402,
            x124 * x266 * x401,
            x124 * x263 * x405,
            x148 * x273 * x402,
            x148 * x269 * x401,
            x148 * x266 * x405,
            x148 * x263 * x409,
            x175 * x274 * x402,
            x175 * x273 * x401,
            x175 * x269 * x405,
            x175 * x266 * x409,
            x175 * x263 * x413,
            x233 * x415 * x69,
            x236 * x415 * x91,
            x233 * x417 * x91,
            x124 * x240 * x415,
            x124 * x236 * x417,
            x124 * x233 * x419,
            x148 * x243 * x415,
            x148 * x240 * x417,
            x148 * x236 * x419,
            x148 * x233 * x421,
            x175 * x246 * x415,
            x175 * x243 * x417,
            x175 * x240 * x419,
            x175 * x236 * x421,
            x175 * x233 * x425,
            x208 * x428 * x69,
            x210 * x428 * x91,
            x208 * x431 * x91,
            x124 * x212 * x428,
            x124 * x210 * x431,
            x124 * x208 * x435,
            x148 * x214 * x428,
            x148 * x212 * x431,
            x148 * x210 * x435,
            x148 * x208 * x438,
            x175 * x218 * x428,
            x175 * x214 * x431,
            x175 * x212 * x435,
            x175 * x210 * x438,
            x175 * x208 * x441,
            x136 * x445 * x69,
            x189 * x445 * x91,
            x136 * x449 * x91,
            x124 * x139 * x445,
            x124 * x189 * x449,
            x124 * x136 * x455,
            x148 * x171 * x445,
            x139 * x148 * x449,
            x148 * x189 * x455,
            x136 * x148 * x461,
            x175 * x197 * x445,
            x171 * x175 * x449,
            x139 * x175 * x455,
            x175 * x189 * x461,
            x136 * x175 * x465,
            x288 * x402 * x64,
            x289 * x402 * x44,
            x288 * x401 * x44,
            x17 * x290 * x402,
            x17 * x289 * x401,
            x17 * x288 * x405,
            x15 * x291 * x402,
            x15 * x290 * x401,
            x15 * x289 * x405,
            x15 * x288 * x409,
            x292 * x395 * x398,
            x11 * x291 * x401,
            x11 * x290 * x405,
            x11 * x289 * x409,
            x11 * x288 * x413,
            x263 * x415 * x64,
            x266 * x415 * x44,
            x263 * x417 * x44,
            x17 * x269 * x415,
            x17 * x266 * x417,
            x17 * x263 * x419,
            x15 * x273 * x415,
            x15 * x269 * x417,
            x15 * x266 * x419,
            x15 * x263 * x421,
            x11 * x274 * x415,
            x11 * x273 * x417,
            x11 * x269 * x419,
            x11 * x266 * x421,
            x11 * x263 * x425,
            x233 * x428 * x64,
            x236 * x428 * x44,
            x233 * x431 * x44,
            x17 * x240 * x428,
            x17 * x236 * x431,
            x17 * x233 * x435,
            x15 * x243 * x428,
            x15 * x240 * x431,
            x15 * x236 * x435,
            x15 * x233 * x438,
            x11 * x246 * x428,
            x11 * x243 * x431,
            x11 * x240 * x435,
            x11 * x236 * x438,
            x11 * x233 * x441,
            x208 * x445 * x64,
            x210 * x44 * x445,
            x208 * x44 * x449,
            x17 * x212 * x445,
            x17 * x210 * x449,
            x17 * x208 * x455,
            x15 * x214 * x445,
            x15 * x212 * x449,
            x15 * x210 * x455,
            x15 * x208 * x461,
            x11 * x218 * x445,
            x11 * x214 * x449,
            x11 * x212 * x455,
            x11 * x210 * x461,
            x206 * x397 * x465,
            x136 * x466 * x64,
            x189 * x44 * x466,
            x136 * x44 * x467,
            x139 * x17 * x466,
            x17 * x189 * x467,
            x136 * x17 * x469,
            x15 * x171 * x466,
            x139 * x15 * x467,
            x15 * x189 * x469,
            x136 * x15 * x470,
            x11 * x197 * x466,
            x11 * x171 * x467,
            x11 * x139 * x469,
            x110 * x397 * x470,
            x10 * x471,
            x163 * x473 * x474,
            x163 * x477 * x478,
            x188 * x473 * x478,
            x163 * x481 * x482,
            x188 * x477 * x482,
            x167 * x473 * x482,
            x163 * x487 * x488,
            x188 * x481 * x488,
            x167 * x477 * x488,
            x193 * x473 * x488,
            x163 * x492 * x493,
            x188 * x487 * x493,
            x167 * x481 * x493,
            x193 * x477 * x493,
            x205 * x473 * x493,
            x163 * x300 * x495,
            x163 * x308 * x497,
            x188 * x308 * x495,
            x163 * x317 * x499,
            x188 * x317 * x497,
            x167 * x317 * x495,
            x163 * x327 * x501,
            x188 * x327 * x499,
            x167 * x327 * x497,
            x193 * x327 * x495,
            x163 * x334 * x505,
            x188 * x334 * x501,
            x167 * x334 * x499,
            x193 * x334 * x497,
            x205 * x334 * x495,
            x221 * x300 * x473,
            x221 * x308 * x477,
            x223 * x308 * x473,
            x221 * x317 * x481,
            x223 * x317 * x477,
            x225 * x317 * x473,
            x221 * x327 * x487,
            x223 * x327 * x481,
            x225 * x327 * x477,
            x227 * x327 * x473,
            x221 * x334 * x492,
            x223 * x334 * x487,
            x225 * x334 * x481,
            x227 * x334 * x477,
            x231 * x334 * x473,
            x163 * x508 * x95,
            x119 * x163 * x511,
            x119 * x188 * x508,
            x145 * x163 * x515,
            x145 * x188 * x511,
            x145 * x167 * x508,
            x163 * x182 * x518,
            x182 * x188 * x515,
            x167 * x182 * x511,
            x182 * x193 * x508,
            x163 * x333 * x521,
            x188 * x333 * x518,
            x167 * x333 * x515,
            x193 * x333 * x511,
            x205 * x333 * x508,
            x221 * x495 * x95,
            x119 * x221 * x497,
            x119 * x223 * x495,
            x145 * x221 * x499,
            x145 * x223 * x497,
            x145 * x225 * x495,
            x182 * x221 * x501,
            x182 * x223 * x499,
            x182 * x225 * x497,
            x182 * x227 * x495,
            x221 * x333 * x505,
            x223 * x333 * x501,
            x225 * x333 * x499,
            x227 * x333 * x497,
            x231 * x333 * x495,
            x248 * x473 * x95,
            x119 * x248 * x477,
            x119 * x251 * x473,
            x145 * x248 * x481,
            x145 * x251 * x477,
            x145 * x255 * x473,
            x182 * x248 * x487,
            x182 * x251 * x481,
            x182 * x255 * x477,
            x182 * x258 * x473,
            x248 * x333 * x492,
            x251 * x333 * x487,
            x255 * x333 * x481,
            x258 * x333 * x477,
            x261 * x333 * x473,
            x163 * x522 * x77,
            x163 * x525 * x89,
            x188 * x522 * x89,
            x122 * x163 * x528,
            x122 * x188 * x525,
            x122 * x167 * x522,
            x163 * x180 * x532,
            x180 * x188 * x528,
            x167 * x180 * x525,
            x180 * x193 * x522,
            x2 * x534,
            x135 * x532 * x535,
            x167 * x173 * x528,
            x173 * x193 * x525,
            x173 * x205 * x522,
            x221 * x508 * x77,
            x221 * x511 * x89,
            x223 * x508 * x89,
            x122 * x221 * x515,
            x122 * x223 * x511,
            x122 * x225 * x508,
            x180 * x221 * x518,
            x180 * x223 * x515,
            x180 * x225 * x511,
            x180 * x227 * x508,
            x219 * x521 * x535,
            x173 * x223 * x518,
            x173 * x225 * x515,
            x173 * x227 * x511,
            x173 * x231 * x508,
            x248 * x495 * x77,
            x248 * x497 * x89,
            x251 * x495 * x89,
            x122 * x248 * x499,
            x122 * x251 * x497,
            x122 * x255 * x495,
            x180 * x248 * x501,
            x180 * x251 * x499,
            x180 * x255 * x497,
            x180 * x258 * x495,
            x173 * x248 * x505,
            x173 * x251 * x501,
            x173 * x255 * x499,
            x173 * x258 * x497,
            x173 * x261 * x495,
            x276 * x473 * x77,
            x276 * x477 * x89,
            x279 * x473 * x89,
            x122 * x276 * x481,
            x122 * x279 * x477,
            x122 * x282 * x473,
            x180 * x276 * x487,
            x180 * x279 * x481,
            x180 * x282 * x477,
            x180 * x286 * x473,
            x173 * x276 * x492,
            x173 * x279 * x487,
            x173 * x282 * x481,
            x173 * x286 * x477,
            x173 * x287 * x473,
            x163 * x536 * x62,
            x163 * x51 * x537,
            x188 * x51 * x536,
            x163 * x49 * x538,
            x188 * x49 * x537,
            x167 * x49 * x536,
            x4 * x539,
            x135 * x538 * x540,
            x167 * x537 * x9,
            x193 * x536 * x9,
            x393
            * (
                x206 * x533
                + x3 * (2 * x384 + 2 * x386 + 3 * x519 + 3 * x520 + 4 * x529 + 4 * x531)
            ),
            x135 * x539,
            x167 * x538 * x8,
            x193 * x537 * x8,
            x205 * x536 * x8,
            x221 * x522 * x62,
            x221 * x51 * x525,
            x223 * x51 * x522,
            x221 * x49 * x528,
            x223 * x49 * x525,
            x225 * x49 * x522,
            x219 * x532 * x540,
            x223 * x528 * x9,
            x225 * x525 * x9,
            x227 * x522 * x9,
            x219 * x534,
            x223 * x532 * x8,
            x225 * x528 * x8,
            x227 * x525 * x8,
            x231 * x522 * x8,
            x248 * x508 * x62,
            x248 * x51 * x511,
            x251 * x508 * x51,
            x248 * x49 * x515,
            x251 * x49 * x511,
            x255 * x49 * x508,
            x248 * x518 * x9,
            x251 * x515 * x9,
            x255 * x511 * x9,
            x258 * x508 * x9,
            x248 * x521 * x8,
            x251 * x518 * x8,
            x255 * x515 * x8,
            x258 * x511 * x8,
            x261 * x508 * x8,
            x276 * x495 * x62,
            x276 * x497 * x51,
            x279 * x495 * x51,
            x276 * x49 * x499,
            x279 * x49 * x497,
            x282 * x49 * x495,
            x276 * x501 * x9,
            x279 * x499 * x9,
            x282 * x497 * x9,
            x286 * x495 * x9,
            x276 * x505 * x8,
            x279 * x501 * x8,
            x282 * x499 * x8,
            x286 * x497 * x8,
            x287 * x495 * x8,
            x293 * x473 * x62,
            x293 * x477 * x51,
            x294 * x473 * x51,
            x293 * x481 * x49,
            x294 * x477 * x49,
            x295 * x473 * x49,
            x293 * x487 * x9,
            x294 * x481 * x9,
            x295 * x477 * x9,
            x296 * x473 * x9,
            x293 * x492 * x8,
            x294 * x487 * x8,
            x295 * x481 * x8,
            x296 * x477 * x8,
            x297 * x473 * x8,
            x298 * x399 * x474,
            x303 * x402 * x478,
            x310 * x401 * x478,
            x313 * x402 * x482,
            x303 * x401 * x482,
            x310 * x405 * x482,
            x322 * x402 * x488,
            x313 * x401 * x488,
            x303 * x405 * x488,
            x310 * x409 * x488,
            x332 * x402 * x493,
            x322 * x401 * x493,
            x313 * x405 * x493,
            x303 * x409 * x493,
            x310 * x413 * x493,
            x300 * x337 * x402,
            x308 * x339 * x402,
            x308 * x337 * x401,
            x317 * x341 * x402,
            x317 * x339 * x401,
            x317 * x337 * x405,
            x327 * x343 * x402,
            x327 * x341 * x401,
            x327 * x339 * x405,
            x327 * x337 * x409,
            x334 * x347 * x402,
            x334 * x343 * x401,
            x334 * x341 * x405,
            x334 * x339 * x409,
            x334 * x337 * x413,
            x300 * x310 * x415,
            x303 * x308 * x415,
            x308 * x310 * x417,
            x313 * x317 * x415,
            x303 * x317 * x417,
            x310 * x317 * x419,
            x322 * x327 * x415,
            x313 * x327 * x417,
            x303 * x327 * x419,
            x310 * x327 * x421,
            x332 * x334 * x415,
            x322 * x334 * x417,
            x313 * x334 * x419,
            x303 * x334 * x421,
            x310 * x334 * x425,
            x350 * x402 * x95,
            x119 * x353 * x402,
            x119 * x350 * x401,
            x145 * x357 * x402,
            x145 * x353 * x401,
            x145 * x350 * x405,
            x182 * x360 * x402,
            x182 * x357 * x401,
            x182 * x353 * x405,
            x182 * x350 * x409,
            x333 * x363 * x402,
            x333 * x360 * x401,
            x333 * x357 * x405,
            x333 * x353 * x409,
            x333 * x350 * x413,
            x337 * x415 * x95,
            x119 * x339 * x415,
            x119 * x337 * x417,
            x145 * x341 * x415,
            x145 * x339 * x417,
            x145 * x337 * x419,
            x182 * x343 * x415,
            x182 * x341 * x417,
            x182 * x339 * x419,
            x182 * x337 * x421,
            x333 * x347 * x415,
            x333 * x343 * x417,
            x333 * x341 * x419,
            x333 * x339 * x421,
            x333 * x337 * x425,
            x310 * x428 * x95,
            x119 * x303 * x428,
            x119 * x310 * x431,
            x145 * x313 * x428,
            x145 * x303 * x431,
            x145 * x310 * x435,
            x182 * x322 * x428,
            x182 * x313 * x431,
            x182 * x303 * x435,
            x182 * x310 * x438,
            x332 * x333 * x428,
            x322 * x333 * x431,
            x313 * x333 * x435,
            x303 * x333 * x438,
            x310 * x333 * x441,
            x367 * x402 * x77,
            x371 * x402 * x89,
            x367 * x401 * x89,
            x122 * x377 * x402,
            x122 * x371 * x401,
            x122 * x367 * x405,
            x180 * x383 * x402,
            x180 * x377 * x401,
            x180 * x371 * x405,
            x180 * x367 * x409,
            x387 * x398 * x535,
            x173 * x383 * x401,
            x173 * x377 * x405,
            x173 * x371 * x409,
            x173 * x367 * x413,
            x350 * x415 * x77,
            x353 * x415 * x89,
            x350 * x417 * x89,
            x122 * x357 * x415,
            x122 * x353 * x417,
            x122 * x350 * x419,
            x180 * x360 * x415,
            x180 * x357 * x417,
            x180 * x353 * x419,
            x180 * x350 * x421,
            x173 * x363 * x415,
            x173 * x360 * x417,
            x173 * x357 * x419,
            x173 * x353 * x421,
            x173 * x350 * x425,
            x337 * x428 * x77,
            x339 * x428 * x89,
            x337 * x431 * x89,
            x122 * x341 * x428,
            x122 * x339 * x431,
            x122 * x337 * x435,
            x180 * x343 * x428,
            x180 * x341 * x431,
            x180 * x339 * x435,
            x180 * x337 * x438,
            x173 * x347 * x428,
            x173 * x343 * x431,
            x173 * x341 * x435,
            x173 * x339 * x438,
            x173 * x337 * x441,
            x310 * x445 * x77,
            x303 * x445 * x89,
            x310 * x449 * x89,
            x122 * x313 * x445,
            x122 * x303 * x449,
            x122 * x310 * x455,
            x180 * x322 * x445,
            x180 * x313 * x449,
            x180 * x303 * x455,
            x180 * x310 * x461,
            x173 * x332 * x445,
            x173 * x322 * x449,
            x173 * x313 * x455,
            x173 * x303 * x461,
            x2 * x465 * x541,
            x388 * x402 * x62,
            x389 * x402 * x51,
            x388 * x401 * x51,
            x391 * x402 * x49,
            x389 * x401 * x49,
            x388 * x405 * x49,
            x392 * x398 * x540,
            x391 * x401 * x9,
            x389 * x405 * x9,
            x388 * x409 * x9,
            x394 * x398,
            x392 * x401 * x8,
            x391 * x405 * x8,
            x389 * x409 * x8,
            x388 * x413 * x8,
            x367 * x415 * x62,
            x371 * x415 * x51,
            x367 * x417 * x51,
            x377 * x415 * x49,
            x371 * x417 * x49,
            x367 * x419 * x49,
            x383 * x415 * x9,
            x377 * x417 * x9,
            x371 * x419 * x9,
            x367 * x421 * x9,
            x387 * x415 * x8,
            x383 * x417 * x8,
            x377 * x419 * x8,
            x371 * x421 * x8,
            x367 * x425 * x8,
            x350 * x428 * x62,
            x353 * x428 * x51,
            x350 * x431 * x51,
            x357 * x428 * x49,
            x353 * x431 * x49,
            x350 * x435 * x49,
            x360 * x428 * x9,
            x357 * x431 * x9,
            x353 * x435 * x9,
            x350 * x438 * x9,
            x363 * x428 * x8,
            x360 * x431 * x8,
            x357 * x435 * x8,
            x353 * x438 * x8,
            x350 * x441 * x8,
            x337 * x445 * x62,
            x339 * x445 * x51,
            x337 * x449 * x51,
            x341 * x445 * x49,
            x339 * x449 * x49,
            x337 * x455 * x49,
            x343 * x445 * x9,
            x341 * x449 * x9,
            x339 * x455 * x9,
            x337 * x461 * x9,
            x347 * x445 * x8,
            x343 * x449 * x8,
            x341 * x455 * x8,
            x339 * x461 * x8,
            x337 * x465 * x8,
            x310 * x466 * x62,
            x303 * x466 * x51,
            x310 * x467 * x51,
            x313 * x466 * x49,
            x303 * x467 * x49,
            x310 * x469 * x49,
            x322 * x466 * x9,
            x313 * x467 * x9,
            x303 * x469 * x9,
            x4 * x470 * x541,
            x332 * x466 * x8,
            x322 * x467 * x8,
            x313 * x469 * x8,
            x303 * x470 * x8,
            x298 * x471,
            x136 * x474 * x543,
            x189 * x478 * x543,
            x136 * x478 * x546,
            x139 * x482 * x543,
            x189 * x482 * x546,
            x136 * x482 * x549,
            x171 * x488 * x543,
            x139 * x488 * x546,
            x189 * x488 * x549,
            x136 * x488 * x554,
            x197 * x493 * x543,
            x171 * x493 * x546,
            x139 * x493 * x549,
            x189 * x493 * x554,
            x136 * x493 * x558,
            x208 * x300 * x543,
            x210 * x308 * x543,
            x208 * x308 * x546,
            x212 * x317 * x543,
            x210 * x317 * x546,
            x208 * x317 * x549,
            x214 * x327 * x543,
            x212 * x327 * x546,
            x210 * x327 * x549,
            x208 * x327 * x554,
            x218 * x334 * x543,
            x214 * x334 * x546,
            x212 * x334 * x549,
            x210 * x334 * x554,
            x208 * x334 * x558,
            x136 * x300 * x560,
            x189 * x308 * x560,
            x136 * x308 * x562,
            x139 * x317 * x560,
            x189 * x317 * x562,
            x136 * x317 * x564,
            x171 * x327 * x560,
            x139 * x327 * x562,
            x189 * x327 * x564,
            x136 * x327 * x566,
            x197 * x334 * x560,
            x171 * x334 * x562,
            x139 * x334 * x564,
            x189 * x334 * x566,
            x136 * x334 * x570,
            x233 * x543 * x95,
            x119 * x236 * x543,
            x119 * x233 * x546,
            x145 * x240 * x543,
            x145 * x236 * x546,
            x145 * x233 * x549,
            x182 * x243 * x543,
            x182 * x240 * x546,
            x182 * x236 * x549,
            x182 * x233 * x554,
            x246 * x333 * x543,
            x243 * x333 * x546,
            x240 * x333 * x549,
            x236 * x333 * x554,
            x233 * x333 * x558,
            x208 * x560 * x95,
            x119 * x210 * x560,
            x119 * x208 * x562,
            x145 * x212 * x560,
            x145 * x210 * x562,
            x145 * x208 * x564,
            x182 * x214 * x560,
            x182 * x212 * x562,
            x182 * x210 * x564,
            x182 * x208 * x566,
            x218 * x333 * x560,
            x214 * x333 * x562,
            x212 * x333 * x564,
            x210 * x333 * x566,
            x208 * x333 * x570,
            x136 * x573 * x95,
            x119 * x189 * x573,
            x119 * x136 * x576,
            x139 * x145 * x573,
            x145 * x189 * x576,
            x136 * x145 * x580,
            x171 * x182 * x573,
            x139 * x182 * x576,
            x182 * x189 * x580,
            x136 * x182 * x583,
            x197 * x333 * x573,
            x171 * x333 * x576,
            x139 * x333 * x580,
            x189 * x333 * x583,
            x136 * x333 * x586,
            x263 * x543 * x77,
            x266 * x543 * x89,
            x263 * x546 * x89,
            x122 * x269 * x543,
            x122 * x266 * x546,
            x122 * x263 * x549,
            x180 * x273 * x543,
            x180 * x269 * x546,
            x180 * x266 * x549,
            x180 * x263 * x554,
            x173 * x274 * x543,
            x173 * x273 * x546,
            x173 * x269 * x549,
            x173 * x266 * x554,
            x173 * x263 * x558,
            x233 * x560 * x77,
            x236 * x560 * x89,
            x233 * x562 * x89,
            x122 * x240 * x560,
            x122 * x236 * x562,
            x122 * x233 * x564,
            x180 * x243 * x560,
            x180 * x240 * x562,
            x180 * x236 * x564,
            x180 * x233 * x566,
            x173 * x246 * x560,
            x173 * x243 * x562,
            x173 * x240 * x564,
            x173 * x236 * x566,
            x173 * x233 * x570,
            x208 * x573 * x77,
            x210 * x573 * x89,
            x208 * x576 * x89,
            x122 * x212 * x573,
            x122 * x210 * x576,
            x122 * x208 * x580,
            x180 * x214 * x573,
            x180 * x212 * x576,
            x180 * x210 * x580,
            x180 * x208 * x583,
            x173 * x218 * x573,
            x173 * x214 * x576,
            x173 * x212 * x580,
            x173 * x210 * x583,
            x206 * x586 * x587,
            x136 * x588 * x77,
            x189 * x588 * x89,
            x136 * x591 * x89,
            x122 * x139 * x588,
            x122 * x189 * x591,
            x122 * x136 * x594,
            x171 * x180 * x588,
            x139 * x180 * x591,
            x180 * x189 * x594,
            x136 * x180 * x598,
            x173 * x197 * x588,
            x171 * x173 * x591,
            x139 * x173 * x594,
            x110 * x587 * x598,
            x2 * x600,
            x288 * x543 * x62,
            x289 * x51 * x543,
            x288 * x51 * x546,
            x290 * x49 * x543,
            x289 * x49 * x546,
            x288 * x49 * x549,
            x291 * x543 * x9,
            x290 * x546 * x9,
            x289 * x549 * x9,
            x288 * x554 * x9,
            x292 * x543 * x8,
            x291 * x546 * x8,
            x290 * x549 * x8,
            x289 * x554 * x8,
            x288 * x558 * x8,
            x263 * x560 * x62,
            x266 * x51 * x560,
            x263 * x51 * x562,
            x269 * x49 * x560,
            x266 * x49 * x562,
            x263 * x49 * x564,
            x273 * x560 * x9,
            x269 * x562 * x9,
            x266 * x564 * x9,
            x263 * x566 * x9,
            x274 * x560 * x8,
            x273 * x562 * x8,
            x269 * x564 * x8,
            x266 * x566 * x8,
            x263 * x570 * x8,
            x233 * x573 * x62,
            x236 * x51 * x573,
            x233 * x51 * x576,
            x240 * x49 * x573,
            x236 * x49 * x576,
            x233 * x49 * x580,
            x243 * x573 * x9,
            x240 * x576 * x9,
            x236 * x580 * x9,
            x233 * x583 * x9,
            x246 * x573 * x8,
            x243 * x576 * x8,
            x240 * x580 * x8,
            x236 * x583 * x8,
            x233 * x586 * x8,
            x208 * x588 * x62,
            x210 * x51 * x588,
            x208 * x51 * x591,
            x212 * x49 * x588,
            x210 * x49 * x591,
            x208 * x49 * x594,
            x214 * x588 * x9,
            x212 * x591 * x9,
            x210 * x594 * x9,
            x206 * x598 * x601,
            x218 * x588 * x8,
            x214 * x591 * x8,
            x212 * x594 * x8,
            x210 * x598 * x8,
            x206 * x600,
            x136 * x602 * x62,
            x189 * x51 * x602,
            x136 * x51 * x603,
            x139 * x49 * x602,
            x189 * x49 * x603,
            x136 * x49 * x604,
            x171 * x602 * x9,
            x139 * x603 * x9,
            x110 * x601 * x604,
            x4 * x605,
            x197 * x602 * x8,
            x171 * x603 * x8,
            x139 * x604 * x8,
            x110 * x605,
            x396
            * (
                x219 * x599
                + x3 * (2 * x462 + 2 * x464 + 3 * x584 + 3 * x585 + 4 * x595 + 4 * x597)
            ),
        ]
    )
