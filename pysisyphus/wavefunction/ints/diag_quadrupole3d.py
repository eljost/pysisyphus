import numpy

"""

    Diagonal of the quadrupole moment matrix with operators x², y², z².

    for rr in (xx, yy, zz):
        for bf_a in basis_functions_a:
            for bf_b in basis_functions_b:
                    quadrupole_integrals(bf_a, bf_b, rr)

"""


def diag_quadrupole3d_00(a, A, b, B, C):
    """Cartesian 3D (ss) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = a * b * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x5 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x6 = x4 * x5
    x7 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x8 = numpy.pi * x0 * x7
    x9 = x2 * x5
    x10 = x5 * x7

    # 3 item(s)
    return numpy.array(
        [
            x2 * x8 * (x3 * x6 + x6 * (x0 * (a * A[0] + b * B[0]) - C[0]) ** 2),
            x4 * x8 * (x3 * x9 + x9 * (x0 * (a * A[1] + b * B[1]) - C[1]) ** 2),
            numpy.pi
            * x0
            * x2
            * x4
            * (x10 * x3 + x10 * (x0 * (a * A[2] + b * B[2]) - C[2]) ** 2),
        ]
    )


def diag_quadrupole3d_01(a, A, b, B, C):
    """Cartesian 3D (sp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - C[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = a * b * x0
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x7 = x5 * x6
    x8 = x3 * x7
    x9 = -x1 - B[0]
    x10 = x2 ** 2 * x7 + x8
    x11 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x12 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x13 = numpy.pi * x0 * x12
    x14 = x11 * x13
    x15 = -x0 * (a * A[1] + b * B[1])
    x16 = -x15 - B[1]
    x17 = x10 * x14
    x18 = -x0 * (a * A[2] + b * B[2])
    x19 = -x18 - B[2]
    x20 = x11 * x6
    x21 = x20 * x3
    x22 = -x15 - C[1]
    x23 = x20 * x22 ** 2 + x21
    x24 = x13 * x5
    x25 = x23 * x24
    x26 = x12 * x6
    x27 = x26 * x3
    x28 = -x18 - C[2]
    x29 = x26 * x28 ** 2 + x27
    x30 = numpy.pi * x0 * x11 * x5
    x31 = x29 * x30

    # 9 item(s)
    return numpy.array(
        [
            x14 * (x10 * x9 + 2 * x2 * x8),
            x16 * x17,
            x17 * x19,
            x25 * x9,
            x24 * (x16 * x23 + 2 * x21 * x22),
            x19 * x25,
            x31 * x9,
            x16 * x31,
            x30 * (x19 * x29 + 2 * x27 * x28),
        ]
    )


def diag_quadrupole3d_02(a, A, b, B, C):
    """Cartesian 3D (sd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x11 = 2 * x8
    x12 = x6 + x9
    x13 = x10 * x12 + x11 * x6
    x14 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x15 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x16 = numpy.pi * x1 * x15
    x17 = x14 * x16
    x18 = -x1 * (a * A[1] + b * B[1])
    x19 = -x18 - B[1]
    x20 = x13 * x17
    x21 = -x1 * (a * A[2] + b * B[2])
    x22 = -x21 - B[2]
    x23 = x14 * x4
    x24 = x0 * x23
    x25 = x19 ** 2 * x23 + x24
    x26 = x15 * x4
    x27 = x0 * x26
    x28 = x22 ** 2 * x26 + x27
    x29 = -x18 - C[1]
    x30 = x23 * x29 ** 2
    x31 = x24 + x30
    x32 = x10 ** 2 * x5 + x6
    x33 = 2 * x29
    x34 = x19 * x31 + x24 * x33
    x35 = x16 * x3
    x36 = x34 * x35
    x37 = -x21 - C[2]
    x38 = x26 * x37 ** 2
    x39 = x27 + x38
    x40 = numpy.pi * x1 * x14 * x3
    x41 = 2 * x37
    x42 = x22 * x39 + x27 * x41
    x43 = x40 * x42

    # 18 item(s)
    return numpy.array(
        [
            x17 * (x0 * (x10 * x11 * x5 + 3 * x6 + x9) + x10 * x13),
            x19 * x20,
            x20 * x22,
            x12 * x25 * x26,
            x12 * x17 * x19 * x22,
            x12 * x23 * x28,
            x26 * x31 * x32,
            x10 * x36,
            x10 * x22 * x31 * x35,
            x35 * (x0 * (x19 * x23 * x33 + 3 * x24 + x30) + x19 * x34),
            x22 * x36,
            x28 * x31 * x5,
            x23 * x32 * x39,
            x10 * x19 * x39 * x40,
            x10 * x43,
            x25 * x39 * x5,
            x19 * x43,
            x40 * (x0 * (x22 * x26 * x41 + 3 * x27 + x38) + x22 * x42),
        ]
    )


def diag_quadrupole3d_03(a, A, b, B, C):
    """Cartesian 3D (sf) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - B[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = a * b * x0
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x7 = x5 * x6
    x8 = x3 * x7
    x9 = -x1 - C[0]
    x10 = x7 * x9 ** 2
    x11 = 2 * x2
    x12 = x7 * x9
    x13 = 2 * x3
    x14 = x10 + x8
    x15 = x14 * x2
    x16 = x12 * x13 + x15
    x17 = x16 * x2 + x3 * (x10 + x11 * x12 + 3 * x8)
    x18 = x2 * x5
    x19 = x18 * x6
    x20 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x21 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x22 = numpy.pi * x0 * x21
    x23 = x20 * x22
    x24 = -x0 * (a * A[1] + b * B[1])
    x25 = -x24 - B[1]
    x26 = x17 * x23
    x27 = -x0 * (a * A[2] + b * B[2])
    x28 = -x27 - B[2]
    x29 = x20 * x6
    x30 = x29 * x3
    x31 = x25 ** 2 * x29 + x30
    x32 = x21 * x6
    x33 = x3 * x32
    x34 = x28 ** 2 * x32 + x33
    x35 = x25 * x29
    x36 = x13 * x35 + x25 * x31
    x37 = x28 * x32
    x38 = x13 * x37 + x28 * x34
    x39 = -x24 - C[1]
    x40 = x29 * x39 ** 2
    x41 = x30 + x40
    x42 = x2 ** 2 * x7 + x8
    x43 = x13 * x19 + x2 * x42
    x44 = x29 * x39
    x45 = x25 * x41
    x46 = x13 * x44 + x45
    x47 = 2 * x25
    x48 = x25 * x46 + x3 * (3 * x30 + x40 + x44 * x47)
    x49 = x18 * x22
    x50 = x22 * x5
    x51 = -x27 - C[2]
    x52 = x32 * x51 ** 2
    x53 = x33 + x52
    x54 = x32 * x51
    x55 = x28 * x53
    x56 = x13 * x54 + x55
    x57 = numpy.pi * x0 * x20
    x58 = x18 * x57
    x59 = 2 * x28
    x60 = x28 * x56 + x3 * (3 * x33 + x52 + x54 * x59)
    x61 = x5 * x57

    # 30 item(s)
    return numpy.array(
        [
            x23
            * (
                x17 * x2
                + x3
                * (x11 * (x19 * x9 + x8) + 4 * x12 * x3 + x13 * (x12 + x19) + 2 * x15)
            ),
            x25 * x26,
            x26 * x28,
            x16 * x31 * x32,
            x16 * x23 * x25 * x28,
            x16 * x29 * x34,
            x14 * x32 * x36,
            x14 * x31 * x37,
            x14 * x34 * x35,
            x14 * x29 * x38,
            x32 * x41 * x43,
            x32 * x42 * x46,
            x37 * x41 * x42,
            x48 * x49,
            x28 * x46 * x49,
            x19 * x34 * x41,
            x50
            * (
                x25 * x48
                + x3
                * (x13 * (x35 + x44) + 4 * x30 * x39 + 2 * x45 + x47 * (x30 + x35 * x39))
            ),
            x28 * x48 * x50,
            x34 * x46 * x7,
            x38 * x41 * x7,
            x29 * x43 * x53,
            x35 * x42 * x53,
            x29 * x42 * x56,
            x19 * x31 * x53,
            x25 * x56 * x58,
            x58 * x60,
            x36 * x53 * x7,
            x31 * x56 * x7,
            x25 * x60 * x61,
            x61
            * (
                x28 * x60
                + x3
                * (x13 * (x37 + x54) + 4 * x33 * x51 + 2 * x55 + x59 * (x33 + x37 * x51))
            ),
        ]
    )


def diag_quadrupole3d_04(a, A, b, B, C):
    """Cartesian 3D (sg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x11 = 2 * x3
    x12 = -x2 - C[0]
    x13 = x12 * x7
    x14 = x10 + x11 * x13
    x15 = 2 * x0
    x16 = x12 ** 2 * x7
    x17 = x0 * (x14 + x16)
    x18 = x16 + x9
    x19 = x18 * x3
    x20 = x13 * x15 + x19
    x21 = x20 * x3
    x22 = x3 * x7
    x23 = x0 * (x13 + x22)
    x24 = x3 * (x12 * x22 + x9)
    x25 = x17 + x21
    x26 = x0 * (4 * x0 * x13 + 2 * x19 + 2 * x23 + 2 * x24) + x25 * x3
    x27 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x28 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x29 = numpy.pi * x1 * x28
    x30 = x27 * x29
    x31 = -x1 * (a * A[1] + b * B[1])
    x32 = -x31 - B[1]
    x33 = x26 * x30
    x34 = -x1 * (a * A[2] + b * B[2])
    x35 = -x34 - B[2]
    x36 = x27 * x6
    x37 = x0 * x36
    x38 = x32 ** 2 * x36
    x39 = x37 + x38
    x40 = x28 * x6
    x41 = x0 * x40
    x42 = x35 ** 2 * x40
    x43 = x41 + x42
    x44 = x32 * x36
    x45 = x15 * x44 + x32 * x39
    x46 = x35 * x40
    x47 = x15 * x46 + x35 * x43
    x48 = 3 * x37
    x49 = x0 * (3 * x38 + x48) + x32 * x45
    x50 = 3 * x41
    x51 = x0 * (3 * x42 + x50) + x35 * x47
    x52 = -x31 - C[1]
    x53 = x36 * x52 ** 2
    x54 = x37 + x53
    x55 = x8 + x9
    x56 = x15 * x22 + x3 * x55
    x57 = x0 * (x10 + 3 * x8) + x3 * x56
    x58 = x36 * x52
    x59 = x32 * x54
    x60 = x15 * x58 + x59
    x61 = 2 * x32
    x62 = x48 + x58 * x61
    x63 = x0 * (x53 + x62)
    x64 = x32 * x60
    x65 = x63 + x64
    x66 = x0 * (x44 + x58)
    x67 = x32 * (x37 + x44 * x52)
    x68 = x0 * (4 * x37 * x52 + 2 * x59 + 2 * x66 + 2 * x67) + x32 * x65
    x69 = x29 * x5
    x70 = x68 * x69
    x71 = -x34 - C[2]
    x72 = x40 * x71 ** 2
    x73 = x41 + x72
    x74 = x40 * x71
    x75 = x35 * x73
    x76 = x15 * x74 + x75
    x77 = 2 * x35
    x78 = x50 + x74 * x77
    x79 = x0 * (x72 + x78)
    x80 = x35 * x76
    x81 = x79 + x80
    x82 = numpy.pi * x1 * x27 * x5
    x83 = x0 * (x46 + x74)
    x84 = x35 * (x41 + x46 * x71)
    x85 = x0 * (4 * x41 * x71 + 2 * x75 + 2 * x83 + 2 * x84) + x35 * x81
    x86 = x82 * x85

    # 45 item(s)
    return numpy.array(
        [
            x30
            * (
                x0 * (x11 * (x23 + x24) + x15 * (x14 + x8) + 3 * x17 + 3 * x21) + x26 * x3
            ),
            x32 * x33,
            x33 * x35,
            x25 * x39 * x40,
            x25 * x30 * x32 * x35,
            x25 * x36 * x43,
            x20 * x40 * x45,
            x20 * x39 * x46,
            x20 * x43 * x44,
            x20 * x36 * x47,
            x18 * x40 * x49,
            x18 * x45 * x46,
            x18 * x39 * x43,
            x18 * x44 * x47,
            x18 * x36 * x51,
            x40 * x54 * x57,
            x40 * x56 * x60,
            x46 * x54 * x56,
            x40 * x55 * x65,
            x46 * x55 * x60,
            x43 * x54 * x55,
            x3 * x70,
            x3 * x35 * x65 * x69,
            x22 * x43 * x60,
            x22 * x47 * x54,
            x69
            * (
                x0 * (x15 * (x38 + x62) + x61 * (x66 + x67) + 3 * x63 + 3 * x64)
                + x32 * x68
            ),
            x35 * x70,
            x43 * x65 * x7,
            x47 * x60 * x7,
            x51 * x54 * x7,
            x36 * x57 * x73,
            x44 * x56 * x73,
            x36 * x56 * x76,
            x39 * x55 * x73,
            x44 * x55 * x76,
            x36 * x55 * x81,
            x22 * x45 * x73,
            x22 * x39 * x76,
            x3 * x32 * x81 * x82,
            x3 * x86,
            x49 * x7 * x73,
            x45 * x7 * x76,
            x39 * x7 * x81,
            x32 * x86,
            x82
            * (
                x0 * (x15 * (x42 + x78) + x77 * (x83 + x84) + 3 * x79 + 3 * x80)
                + x35 * x85
            ),
        ]
    )


def diag_quadrupole3d_10(a, A, b, B, C):
    """Cartesian 3D (ps) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - C[0]
    x3 = (2 * a + 2 * b) ** (-1.0)
    x4 = a * b * x0
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(numpy.pi) * numpy.sqrt(x0)
    x7 = x5 * x6
    x8 = x3 * x7
    x9 = -x1 - A[0]
    x10 = x2 ** 2 * x7 + x8
    x11 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x12 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x13 = numpy.pi * x0 * x12
    x14 = x11 * x13
    x15 = -x0 * (a * A[1] + b * B[1])
    x16 = -x15 - A[1]
    x17 = x10 * x14
    x18 = -x0 * (a * A[2] + b * B[2])
    x19 = -x18 - A[2]
    x20 = x11 * x6
    x21 = x20 * x3
    x22 = -x15 - C[1]
    x23 = x20 * x22 ** 2 + x21
    x24 = x13 * x5
    x25 = x23 * x24
    x26 = x12 * x6
    x27 = x26 * x3
    x28 = -x18 - C[2]
    x29 = x26 * x28 ** 2 + x27
    x30 = numpy.pi * x0 * x11 * x5
    x31 = x29 * x30

    # 9 item(s)
    return numpy.array(
        [
            x14 * (x10 * x9 + 2 * x2 * x8),
            x16 * x17,
            x17 * x19,
            x25 * x9,
            x24 * (x16 * x23 + 2 * x21 * x22),
            x19 * x25,
            x31 * x9,
            x16 * x31,
            x30 * (x19 * x29 + 2 * x27 * x28),
        ]
    )


def diag_quadrupole3d_11(a, A, b, B, C):
    """Cartesian 3D (pp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x12 = 2 * x8
    x13 = -x7 - A[0]
    x14 = x12 * x6
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
    x40 = x30 * x39 ** 2
    x41 = x29 + x40
    x42 = x11 * x13 + x6
    x43 = 2 * x39
    x44 = x29 * x43
    x45 = x22 * x41 + x44
    x46 = x19 * x3
    x47 = x45 * x46
    x48 = x41 * x46
    x49 = x46 * (x26 * x41 + x44)
    x50 = -x24 - C[2]
    x51 = x33 * x50 ** 2
    x52 = x36 + x51
    x53 = numpy.pi * x1 * x17 * x3
    x54 = x52 * x53
    x55 = 2 * x50
    x56 = x36 * x55
    x57 = x25 * x52 + x56
    x58 = x53 * x57
    x59 = x53 * (x35 * x52 + x56)

    # 27 item(s)
    return numpy.array(
        [
            x20 * (x0 * (x11 * x12 + 3 * x6 + x9) + x13 * x16),
            x22 * x23,
            x23 * x25,
            x26 * x27,
            x15 * x32 * x33,
            x25 * x26 * x34,
            x27 * x35,
            x22 * x34 * x35,
            x15 * x30 * x38,
            x33 * x41 * x42,
            x13 * x47,
            x13 * x25 * x48,
            x10 * x49,
            x46 * (x0 * (3 * x29 + x31 * x43 + x40) + x26 * x45),
            x25 * x49,
            x10 * x35 * x48,
            x35 * x47,
            x38 * x41 * x5,
            x30 * x42 * x52,
            x13 * x22 * x54,
            x13 * x58,
            x10 * x26 * x54,
            x32 * x5 * x52,
            x26 * x58,
            x10 * x59,
            x22 * x59,
            x53 * (x0 * (3 * x36 + x37 * x55 + x51) + x35 * x57),
        ]
    )


def diag_quadrupole3d_12(a, A, b, B, C):
    """Cartesian 3D (pd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x9 = -x1 - C[0]
    x10 = x7 * x9 ** 2
    x11 = -x1 - B[0]
    x12 = 2 * x11
    x13 = x7 * x9
    x14 = x3 * (x10 + x12 * x13 + 3 * x8)
    x15 = 2 * x3
    x16 = x13 * x15
    x17 = x10 + x8
    x18 = x11 * x17
    x19 = x16 + x18
    x20 = x11 * x19 + x14
    x21 = x11 * x5
    x22 = x21 * x6
    x23 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x24 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x25 = numpy.pi * x0 * x24
    x26 = x23 * x25
    x27 = -x0 * (a * A[1] + b * B[1])
    x28 = -x27 - B[1]
    x29 = x26 * (x14 + x19 * x2)
    x30 = -x0 * (a * A[2] + b * B[2])
    x31 = -x30 - B[2]
    x32 = x16 + x17 * x2
    x33 = x23 * x6
    x34 = x3 * x33
    x35 = x28 ** 2 * x33 + x34
    x36 = x24 * x6
    x37 = x35 * x36
    x38 = x26 * x31
    x39 = x3 * x36
    x40 = x31 ** 2 * x36 + x39
    x41 = x33 * x40
    x42 = -x27 - A[1]
    x43 = x20 * x26
    x44 = x28 * x33
    x45 = x34 + x42 * x44
    x46 = x15 * x44 + x35 * x42
    x47 = x31 * x36
    x48 = -x30 - A[2]
    x49 = x39 + x47 * x48
    x50 = x15 * x47 + x40 * x48
    x51 = x11 ** 2 * x7 + x8
    x52 = x15 * x22 + x2 * x51
    x53 = -x27 - C[1]
    x54 = x33 * x53 ** 2
    x55 = x34 + x54
    x56 = x36 * x55
    x57 = x33 * x53
    x58 = x15 * x57
    x59 = x28 * x55
    x60 = x58 + x59
    x61 = x2 * x22 + x8
    x62 = 2 * x28
    x63 = x3 * (3 * x34 + x54 + x57 * x62)
    x64 = x28 * x60 + x63
    x65 = x25 * x5
    x66 = x64 * x65
    x67 = x31 * x65
    x68 = x40 * x7
    x69 = x42 * x55 + x58
    x70 = x42 * x60 + x63
    x71 = x21 * x25
    x72 = -x30 - C[2]
    x73 = x36 * x72 ** 2
    x74 = x39 + x73
    x75 = x33 * x74
    x76 = x36 * x72
    x77 = x15 * x76
    x78 = x31 * x74
    x79 = x77 + x78
    x80 = x7 * x74
    x81 = numpy.pi * x0 * x23
    x82 = x5 * x81
    x83 = x28 * x82
    x84 = 2 * x31
    x85 = x3 * (3 * x39 + x73 + x76 * x84)
    x86 = x31 * x79 + x85
    x87 = x82 * x86
    x88 = x21 * x81
    x89 = x48 * x74 + x77
    x90 = x48 * x79 + x85

    # 54 item(s)
    return numpy.array(
        [
            x26
            * (
                x2 * x20
                + x3
                * (x12 * (x22 * x9 + x8) + 4 * x13 * x3 + x15 * (x13 + x22) + 2 * x18)
            ),
            x28 * x29,
            x29 * x31,
            x32 * x37,
            x28 * x32 * x38,
            x32 * x41,
            x42 * x43,
            x19 * x36 * x45,
            x19 * x38 * x42,
            x17 * x36 * x46,
            x17 * x45 * x47,
            x17 * x41 * x42,
            x43 * x48,
            x19 * x26 * x28 * x48,
            x19 * x33 * x49,
            x17 * x37 * x48,
            x17 * x44 * x49,
            x17 * x33 * x50,
            x52 * x56,
            x36 * x60 * x61,
            x47 * x55 * x61,
            x2 * x66,
            x2 * x60 * x67,
            x2 * x55 * x68,
            x36 * x51 * x69,
            x70 * x71,
            x31 * x69 * x71,
            x65
            * (
                x3
                * (x15 * (x44 + x57) + 4 * x34 * x53 + 2 * x59 + x62 * (x34 + x44 * x53))
                + x42 * x64
            ),
            x67 * x70,
            x68 * x69,
            x48 * x51 * x56,
            x48 * x60 * x71,
            x22 * x49 * x55,
            x48 * x66,
            x49 * x60 * x7,
            x50 * x55 * x7,
            x52 * x75,
            x44 * x61 * x74,
            x33 * x61 * x79,
            x2 * x35 * x80,
            x2 * x79 * x83,
            x2 * x87,
            x42 * x51 * x75,
            x22 * x45 * x74,
            x42 * x79 * x88,
            x46 * x80,
            x45 * x7 * x79,
            x42 * x87,
            x33 * x51 * x89,
            x28 * x88 * x89,
            x88 * x90,
            x35 * x7 * x89,
            x83 * x90,
            x82
            * (
                x3
                * (x15 * (x47 + x76) + 4 * x39 * x72 + 2 * x78 + x84 * (x39 + x47 * x72))
                + x48 * x86
            ),
        ]
    )


def diag_quadrupole3d_13(a, A, b, B, C):
    """Cartesian 3D (pf) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x11 = 2 * x3
    x12 = -x2 - C[0]
    x13 = x12 * x7
    x14 = x10 + x11 * x13
    x15 = 2 * x0
    x16 = x12 ** 2 * x7
    x17 = x0 * (x14 + x16)
    x18 = x13 * x15
    x19 = x16 + x9
    x20 = x19 * x3
    x21 = x18 + x20
    x22 = x21 * x3
    x23 = x3 * x7
    x24 = x0 * (x13 + x23)
    x25 = x3 * (x12 * x23 + x9)
    x26 = -x2 - A[0]
    x27 = x17 + x22
    x28 = x0 * (4 * x0 * x13 + 2 * x20 + 2 * x24 + 2 * x25)
    x29 = x27 * x3 + x28
    x30 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x31 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x32 = numpy.pi * x1 * x31
    x33 = x30 * x32
    x34 = -x1 * (a * A[1] + b * B[1])
    x35 = -x34 - B[1]
    x36 = x33 * (x26 * x27 + x28)
    x37 = -x1 * (a * A[2] + b * B[2])
    x38 = -x37 - B[2]
    x39 = x17 + x21 * x26
    x40 = x30 * x6
    x41 = x0 * x40
    x42 = x35 ** 2 * x40
    x43 = x41 + x42
    x44 = x31 * x6
    x45 = x43 * x44
    x46 = x33 * x38
    x47 = x0 * x44
    x48 = x38 ** 2 * x44
    x49 = x47 + x48
    x50 = x40 * x49
    x51 = x18 + x19 * x26
    x52 = x35 * x40
    x53 = x15 * x52
    x54 = x35 * x43 + x53
    x55 = x44 * x54
    x56 = x38 * x44
    x57 = x15 * x56
    x58 = x38 * x49 + x57
    x59 = x40 * x58
    x60 = -x34 - A[1]
    x61 = x29 * x33
    x62 = x41 + x52 * x60
    x63 = x43 * x60 + x53
    x64 = 3 * x41
    x65 = x0 * (3 * x42 + x64) + x54 * x60
    x66 = -x37 - A[2]
    x67 = x47 + x56 * x66
    x68 = x49 * x66 + x57
    x69 = 3 * x47
    x70 = x0 * (3 * x48 + x69) + x58 * x66
    x71 = x15 * x23
    x72 = x8 + x9
    x73 = x3 * x72 + x71
    x74 = x0 * (x10 + 3 * x8) + x26 * x73
    x75 = -x34 - C[1]
    x76 = x40 * x75 ** 2
    x77 = x41 + x76
    x78 = x44 * x77
    x79 = x26 * x72 + x71
    x80 = x40 * x75
    x81 = x15 * x80
    x82 = x35 * x77
    x83 = x81 + x82
    x84 = x44 * x83
    x85 = 2 * x35
    x86 = x64 + x80 * x85
    x87 = x0 * (x76 + x86)
    x88 = x35 * x83
    x89 = x87 + x88
    x90 = x23 * x26 + x9
    x91 = x0 * (x52 + x80)
    x92 = x35 * (x41 + x52 * x75)
    x93 = x0 * (4 * x41 * x75 + 2 * x82 + 2 * x91 + 2 * x92)
    x94 = x35 * x89 + x93
    x95 = x32 * x5
    x96 = x94 * x95
    x97 = x38 * x95
    x98 = x49 * x7
    x99 = x58 * x7
    x100 = x60 * x77 + x81
    x101 = x60 * x83 + x87
    x102 = x95 * (x60 * x89 + x93)
    x103 = -x37 - C[2]
    x104 = x103 ** 2 * x44
    x105 = x104 + x47
    x106 = x105 * x40
    x107 = x103 * x44
    x108 = x107 * x15
    x109 = x105 * x38
    x110 = x108 + x109
    x111 = x110 * x40
    x112 = 2 * x38
    x113 = x107 * x112 + x69
    x114 = x0 * (x104 + x113)
    x115 = x110 * x38
    x116 = x114 + x115
    x117 = x105 * x7
    x118 = x110 * x7
    x119 = numpy.pi * x1 * x30 * x5
    x120 = x116 * x119
    x121 = x0 * (x107 + x56)
    x122 = x38 * (x103 * x56 + x47)
    x123 = x0 * (4 * x103 * x47 + 2 * x109 + 2 * x121 + 2 * x122)
    x124 = x116 * x38 + x123
    x125 = x119 * x124
    x126 = x105 * x66 + x108
    x127 = x110 * x66 + x114
    x128 = x119 * (x116 * x66 + x123)

    # 90 item(s)
    return numpy.array(
        [
            x33
            * (
                x0 * (x11 * (x24 + x25) + x15 * (x14 + x8) + 3 * x17 + 3 * x22)
                + x26 * x29
            ),
            x35 * x36,
            x36 * x38,
            x39 * x45,
            x35 * x39 * x46,
            x39 * x50,
            x51 * x55,
            x43 * x51 * x56,
            x49 * x51 * x52,
            x51 * x59,
            x60 * x61,
            x27 * x44 * x62,
            x27 * x46 * x60,
            x21 * x44 * x63,
            x21 * x56 * x62,
            x21 * x50 * x60,
            x19 * x44 * x65,
            x19 * x56 * x63,
            x19 * x49 * x62,
            x19 * x59 * x60,
            x61 * x66,
            x27 * x33 * x35 * x66,
            x27 * x40 * x67,
            x21 * x45 * x66,
            x21 * x52 * x67,
            x21 * x40 * x68,
            x19 * x55 * x66,
            x19 * x43 * x67,
            x19 * x52 * x68,
            x19 * x40 * x70,
            x74 * x78,
            x79 * x84,
            x56 * x77 * x79,
            x44 * x89 * x90,
            x56 * x83 * x90,
            x49 * x77 * x90,
            x26 * x96,
            x26 * x89 * x97,
            x26 * x83 * x98,
            x26 * x77 * x99,
            x100 * x44 * x73,
            x101 * x44 * x72,
            x100 * x56 * x72,
            x102 * x3,
            x101 * x3 * x97,
            x100 * x23 * x49,
            x95
            * (
                x0 * (x15 * (x42 + x86) + x85 * (x91 + x92) + 3 * x87 + 3 * x88)
                + x60 * x94
            ),
            x102 * x38,
            x101 * x98,
            x100 * x99,
            x66 * x73 * x78,
            x66 * x72 * x84,
            x67 * x72 * x77,
            x3 * x66 * x89 * x95,
            x23 * x67 * x83,
            x23 * x68 * x77,
            x66 * x96,
            x67 * x7 * x89,
            x68 * x7 * x83,
            x7 * x70 * x77,
            x106 * x74,
            x105 * x52 * x79,
            x111 * x79,
            x105 * x43 * x90,
            x110 * x52 * x90,
            x116 * x40 * x90,
            x117 * x26 * x54,
            x118 * x26 * x43,
            x120 * x26 * x35,
            x125 * x26,
            x106 * x60 * x73,
            x105 * x62 * x72,
            x111 * x60 * x72,
            x105 * x23 * x63,
            x110 * x23 * x62,
            x120 * x3 * x60,
            x117 * x65,
            x118 * x63,
            x116 * x62 * x7,
            x125 * x60,
            x126 * x40 * x73,
            x126 * x52 * x72,
            x127 * x40 * x72,
            x126 * x23 * x43,
            x119 * x127 * x3 * x35,
            x128 * x3,
            x126 * x54 * x7,
            x127 * x43 * x7,
            x128 * x35,
            x119
            * (
                x0 * (x112 * (x121 + x122) + 3 * x114 + 3 * x115 + x15 * (x113 + x48))
                + x124 * x66
            ),
        ]
    )


def diag_quadrupole3d_14(a, A, b, B, C):
    """Cartesian 3D (pg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x12 = 2 * x4
    x13 = -x1 - C[0]
    x14 = x13 * x8
    x15 = x11 + x12 * x14
    x16 = x3 * (x15 + x9)
    x17 = x13 ** 2 * x8
    x18 = x3 * (x15 + x17)
    x19 = 2 * x3
    x20 = x14 * x19
    x21 = x10 + x17
    x22 = x21 * x4
    x23 = x20 + x22
    x24 = x23 * x4
    x25 = x4 * x8
    x26 = x3 * (x14 + x25)
    x27 = x4 * (x10 + x13 * x25)
    x28 = x4 * (x26 + x27)
    x29 = x3 * (2 * x16 + 3 * x18 + 3 * x24 + 2 * x28)
    x30 = x18 + x24
    x31 = x30 * x4
    x32 = x3 * (4 * x10 * x13 + 2 * x22 + 2 * x26 + 2 * x27)
    x33 = x31 + x32
    x34 = x29 + x33 * x4
    x35 = x19 * x25
    x36 = x10 + x9
    x37 = x36 * x4
    x38 = x35 + x37
    x39 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x40 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x41 = numpy.pi * x0 * x40
    x42 = x39 * x41
    x43 = -x0 * (a * A[1] + b * B[1])
    x44 = -x43 - B[1]
    x45 = x42 * (x2 * x33 + x29)
    x46 = -x0 * (a * A[2] + b * B[2])
    x47 = -x46 - B[2]
    x48 = x2 * x30 + x32
    x49 = x39 * x7
    x50 = x3 * x49
    x51 = x44 ** 2 * x49
    x52 = x50 + x51
    x53 = x40 * x7
    x54 = x52 * x53
    x55 = x42 * x47
    x56 = x3 * x53
    x57 = x47 ** 2 * x53
    x58 = x56 + x57
    x59 = x49 * x58
    x60 = x18 + x2 * x23
    x61 = x44 * x49
    x62 = x19 * x61
    x63 = x44 * x52
    x64 = x62 + x63
    x65 = x53 * x64
    x66 = x47 * x53
    x67 = x19 * x66
    x68 = x47 * x58
    x69 = x67 + x68
    x70 = x49 * x69
    x71 = x2 * x21 + x20
    x72 = 3 * x50
    x73 = x3 * (3 * x51 + x72)
    x74 = x44 * x64 + x73
    x75 = x53 * x74
    x76 = 3 * x56
    x77 = x3 * (3 * x57 + x76)
    x78 = x47 * x69 + x77
    x79 = x49 * x78
    x80 = -x43 - A[1]
    x81 = x34 * x42
    x82 = x50 + x61 * x80
    x83 = x52 * x80 + x62
    x84 = x64 * x80 + x73
    x85 = x3 * (8 * x44 * x50 + 4 * x63) + x74 * x80
    x86 = -x46 - A[2]
    x87 = x56 + x66 * x86
    x88 = x58 * x86 + x67
    x89 = x69 * x86 + x77
    x90 = x3 * (8 * x47 * x56 + 4 * x68) + x78 * x86
    x91 = x3 * (x11 + 3 * x9)
    x92 = x38 * x4 + x91
    x93 = x2 * x92 + x3 * (8 * x10 * x4 + 4 * x37)
    x94 = -x43 - C[1]
    x95 = x49 * x94 ** 2
    x96 = x50 + x95
    x97 = x53 * x96
    x98 = x2 * x38 + x91
    x99 = x49 * x94
    x100 = x19 * x99
    x101 = x44 * x96
    x102 = x100 + x101
    x103 = x102 * x53
    x104 = x2 * x36 + x35
    x105 = 2 * x44
    x106 = x105 * x99 + x72
    x107 = x3 * (x106 + x95)
    x108 = x102 * x44
    x109 = x107 + x108
    x110 = x109 * x53
    x111 = x109 * x44
    x112 = x3 * (x61 + x99)
    x113 = x44 * (x50 + x61 * x94)
    x114 = x3 * (2 * x101 + 2 * x112 + 2 * x113 + 4 * x50 * x94)
    x115 = x111 + x114
    x116 = x10 + x2 * x25
    x117 = x3 * (x106 + x51)
    x118 = x44 * (x112 + x113)
    x119 = x3 * (3 * x107 + 3 * x108 + 2 * x117 + 2 * x118)
    x120 = x115 * x44 + x119
    x121 = x41 * x6
    x122 = x120 * x121
    x123 = x121 * x47
    x124 = x58 * x8
    x125 = x69 * x8
    x126 = x78 * x8
    x127 = x100 + x80 * x96
    x128 = x102 * x80 + x107
    x129 = x109 * x80 + x114
    x130 = x121 * (x115 * x80 + x119)
    x131 = -x46 - C[2]
    x132 = x131 ** 2 * x53
    x133 = x132 + x56
    x134 = x133 * x49
    x135 = x131 * x53
    x136 = x135 * x19
    x137 = x133 * x47
    x138 = x136 + x137
    x139 = x138 * x49
    x140 = 2 * x47
    x141 = x135 * x140 + x76
    x142 = x3 * (x132 + x141)
    x143 = x138 * x47
    x144 = x142 + x143
    x145 = x144 * x49
    x146 = x144 * x47
    x147 = x3 * (x135 + x66)
    x148 = x47 * (x131 * x66 + x56)
    x149 = x3 * (4 * x131 * x56 + 2 * x137 + 2 * x147 + 2 * x148)
    x150 = x146 + x149
    x151 = x133 * x8
    x152 = x138 * x8
    x153 = x144 * x8
    x154 = numpy.pi * x0 * x39 * x6
    x155 = x150 * x154
    x156 = x3 * (x141 + x57)
    x157 = x47 * (x147 + x148)
    x158 = x3 * (3 * x142 + 3 * x143 + 2 * x156 + 2 * x157)
    x159 = x150 * x47 + x158
    x160 = x154 * x159
    x161 = x133 * x86 + x136
    x162 = x138 * x86 + x142
    x163 = x144 * x86 + x149
    x164 = x154 * (x150 * x86 + x158)

    # 135 item(s)
    return numpy.array(
        [
            x42
            * (
                x2 * x34
                + x3
                * (
                    x12 * (x16 + x28)
                    + x19 * (3 * x26 + 3 * x27 + x38)
                    + 4 * x31
                    + 4 * x32
                )
            ),
            x44 * x45,
            x45 * x47,
            x48 * x54,
            x44 * x48 * x55,
            x48 * x59,
            x60 * x65,
            x52 * x60 * x66,
            x58 * x60 * x61,
            x60 * x70,
            x71 * x75,
            x64 * x66 * x71,
            x52 * x58 * x71,
            x61 * x69 * x71,
            x71 * x79,
            x80 * x81,
            x33 * x53 * x82,
            x33 * x55 * x80,
            x30 * x53 * x83,
            x30 * x66 * x82,
            x30 * x59 * x80,
            x23 * x53 * x84,
            x23 * x66 * x83,
            x23 * x58 * x82,
            x23 * x70 * x80,
            x21 * x53 * x85,
            x21 * x66 * x84,
            x21 * x58 * x83,
            x21 * x69 * x82,
            x21 * x79 * x80,
            x81 * x86,
            x33 * x42 * x44 * x86,
            x33 * x49 * x87,
            x30 * x54 * x86,
            x30 * x61 * x87,
            x30 * x49 * x88,
            x23 * x65 * x86,
            x23 * x52 * x87,
            x23 * x61 * x88,
            x23 * x49 * x89,
            x21 * x75 * x86,
            x21 * x64 * x87,
            x21 * x52 * x88,
            x21 * x61 * x89,
            x21 * x49 * x90,
            x93 * x97,
            x103 * x98,
            x66 * x96 * x98,
            x104 * x110,
            x102 * x104 * x66,
            x104 * x58 * x96,
            x115 * x116 * x53,
            x109 * x116 * x66,
            x102 * x116 * x58,
            x116 * x69 * x96,
            x122 * x2,
            x115 * x123 * x2,
            x109 * x124 * x2,
            x102 * x125 * x2,
            x126 * x2 * x96,
            x127 * x53 * x92,
            x128 * x38 * x53,
            x127 * x38 * x66,
            x129 * x36 * x53,
            x128 * x36 * x66,
            x127 * x36 * x58,
            x130 * x4,
            x123 * x129 * x4,
            x128 * x25 * x58,
            x127 * x25 * x69,
            x121
            * (
                x120 * x80
                + x3
                * (
                    x105 * (x117 + x118)
                    + 4 * x111
                    + 4 * x114
                    + x19 * (3 * x112 + 3 * x113 + x64)
                )
            ),
            x130 * x47,
            x124 * x129,
            x125 * x128,
            x126 * x127,
            x86 * x92 * x97,
            x103 * x38 * x86,
            x38 * x87 * x96,
            x110 * x36 * x86,
            x102 * x36 * x87,
            x36 * x88 * x96,
            x115 * x121 * x4 * x86,
            x109 * x25 * x87,
            x102 * x25 * x88,
            x25 * x89 * x96,
            x122 * x86,
            x115 * x8 * x87,
            x109 * x8 * x88,
            x102 * x8 * x89,
            x8 * x90 * x96,
            x134 * x93,
            x133 * x61 * x98,
            x139 * x98,
            x104 * x133 * x52,
            x104 * x138 * x61,
            x104 * x145,
            x116 * x133 * x64,
            x116 * x138 * x52,
            x116 * x144 * x61,
            x116 * x150 * x49,
            x151 * x2 * x74,
            x152 * x2 * x64,
            x153 * x2 * x52,
            x155 * x2 * x44,
            x160 * x2,
            x134 * x80 * x92,
            x133 * x38 * x82,
            x139 * x38 * x80,
            x133 * x36 * x83,
            x138 * x36 * x82,
            x145 * x36 * x80,
            x133 * x25 * x84,
            x138 * x25 * x83,
            x144 * x25 * x82,
            x155 * x4 * x80,
            x151 * x85,
            x152 * x84,
            x153 * x83,
            x150 * x8 * x82,
            x160 * x80,
            x161 * x49 * x92,
            x161 * x38 * x61,
            x162 * x38 * x49,
            x161 * x36 * x52,
            x162 * x36 * x61,
            x163 * x36 * x49,
            x161 * x25 * x64,
            x162 * x25 * x52,
            x154 * x163 * x4 * x44,
            x164 * x4,
            x161 * x74 * x8,
            x162 * x64 * x8,
            x163 * x52 * x8,
            x164 * x44,
            x154
            * (
                x159 * x86
                + x3
                * (
                    x140 * (x156 + x157)
                    + 4 * x146
                    + 4 * x149
                    + x19 * (3 * x147 + 3 * x148 + x69)
                )
            ),
        ]
    )


def diag_quadrupole3d_20(a, A, b, B, C):
    """Cartesian 3D (ds) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x11 = 2 * x8
    x12 = x6 + x9
    x13 = x10 * x12 + x11 * x6
    x14 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x15 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x16 = numpy.pi * x1 * x15
    x17 = x14 * x16
    x18 = -x1 * (a * A[1] + b * B[1])
    x19 = -x18 - A[1]
    x20 = x13 * x17
    x21 = -x1 * (a * A[2] + b * B[2])
    x22 = -x21 - A[2]
    x23 = x14 * x4
    x24 = x0 * x23
    x25 = x19 ** 2 * x23 + x24
    x26 = x15 * x4
    x27 = x0 * x26
    x28 = x22 ** 2 * x26 + x27
    x29 = -x18 - C[1]
    x30 = x23 * x29 ** 2
    x31 = x24 + x30
    x32 = x10 ** 2 * x5 + x6
    x33 = 2 * x29
    x34 = x19 * x31 + x24 * x33
    x35 = x16 * x3
    x36 = x34 * x35
    x37 = -x21 - C[2]
    x38 = x26 * x37 ** 2
    x39 = x27 + x38
    x40 = numpy.pi * x1 * x14 * x3
    x41 = 2 * x37
    x42 = x22 * x39 + x27 * x41
    x43 = x40 * x42

    # 18 item(s)
    return numpy.array(
        [
            x17 * (x0 * (x10 * x11 * x5 + 3 * x6 + x9) + x10 * x13),
            x19 * x20,
            x20 * x22,
            x12 * x25 * x26,
            x12 * x17 * x19 * x22,
            x12 * x23 * x28,
            x26 * x31 * x32,
            x10 * x36,
            x10 * x22 * x31 * x35,
            x35 * (x0 * (x19 * x23 * x33 + 3 * x24 + x30) + x19 * x34),
            x22 * x36,
            x28 * x31 * x5,
            x23 * x32 * x39,
            x10 * x19 * x39 * x40,
            x10 * x43,
            x25 * x39 * x5,
            x19 * x43,
            x40 * (x0 * (x22 * x26 * x41 + 3 * x27 + x38) + x22 * x42),
        ]
    )


def diag_quadrupole3d_21(a, A, b, B, C):
    """Cartesian 3D (dp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x14 = x12 + 3 * x13
    x15 = 2 * x3
    x16 = x4 * x9
    x17 = x15 * x16
    x18 = x12 + x13
    x19 = x18 * x5
    x20 = x17 + x19
    x21 = x2 * x20 + x3 * (2 * x11 + x14)
    x22 = x18 * x2
    x23 = 2 * x2
    x24 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x25 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x26 = numpy.pi * x0 * x25
    x27 = x24 * x26
    x28 = -x0 * (a * A[1] + b * B[1])
    x29 = -x28 - B[1]
    x30 = x17 + x22
    x31 = x27 * (x2 * x30 + x3 * (x14 + x16 * x23))
    x32 = -x0 * (a * A[2] + b * B[2])
    x33 = -x32 - B[2]
    x34 = -x28 - A[1]
    x35 = x21 * x27
    x36 = x3 * x8
    x37 = x24 * x36
    x38 = x24 * x8
    x39 = x34 * x38
    x40 = x29 * x39 + x37
    x41 = x25 * x8
    x42 = x27 * x30
    x43 = -x32 - A[2]
    x44 = x25 * x36
    x45 = x41 * x43
    x46 = x33 * x45 + x44
    x47 = x34 ** 2 * x38 + x37
    x48 = x29 * x38
    x49 = x3 * (x39 + x48) + x34 * x40
    x50 = x33 * x41
    x51 = x41 * x43 ** 2 + x44
    x52 = x3 * (x45 + x50) + x43 * x46
    x53 = -x28 - C[1]
    x54 = x38 * x53 ** 2
    x55 = x37 + x54
    x56 = x2 * x9
    x57 = x10 * x2 + x13
    x58 = x2 * x57 + x3 * (x10 + x56)
    x59 = x38 * x53
    x60 = x15 * x59
    x61 = x29 * x55
    x62 = x60 + x61
    x63 = x13 + x2 ** 2 * x9
    x64 = x34 * x55
    x65 = x60 + x64
    x66 = x48 * x53
    x67 = 3 * x37 + x54
    x68 = x3 * (2 * x66 + x67) + x34 * x62
    x69 = x26 * x7
    x70 = x68 * x69
    x71 = x2 * x69
    x72 = 2 * x34
    x73 = x69 * (x3 * (x59 * x72 + x67) + x34 * x65)
    x74 = -x32 - C[2]
    x75 = x41 * x74 ** 2
    x76 = x44 + x75
    x77 = x41 * x74
    x78 = x15 * x77
    x79 = x33 * x76
    x80 = x78 + x79
    x81 = numpy.pi * x0 * x24 * x7
    x82 = x2 * x81
    x83 = x43 * x76
    x84 = x78 + x83
    x85 = x50 * x74
    x86 = 3 * x44 + x75
    x87 = x3 * (2 * x85 + x86) + x43 * x80
    x88 = x81 * x87
    x89 = 2 * x43
    x90 = x81 * (x3 * (x77 * x89 + x86) + x43 * x84)

    # 54 item(s)
    return numpy.array(
        [
            x27
            * (
                x2 * x21
                + x3 * (4 * x13 * x4 + x15 * (x10 + x16) + x19 + x22 + x23 * (x11 + x13))
            ),
            x29 * x31,
            x31 * x33,
            x34 * x35,
            x30 * x40 * x41,
            x33 * x34 * x42,
            x35 * x43,
            x29 * x42 * x43,
            x30 * x38 * x46,
            x20 * x41 * x47,
            x18 * x41 * x49,
            x18 * x47 * x50,
            x20 * x27 * x34 * x43,
            x18 * x40 * x45,
            x18 * x39 * x46,
            x20 * x38 * x51,
            x18 * x48 * x51,
            x18 * x38 * x52,
            x41 * x55 * x58,
            x41 * x62 * x63,
            x50 * x55 * x63,
            x41 * x57 * x65,
            x2 * x70,
            x33 * x65 * x71,
            x45 * x55 * x57,
            x43 * x62 * x71,
            x46 * x55 * x56,
            x5 * x73,
            x69
            * (
                x3 * (x15 * (x48 + x59) + 4 * x37 * x53 + x61 + x64 + x72 * (x37 + x66))
                + x34 * x68
            ),
            x33 * x73,
            x43 * x5 * x65 * x69,
            x43 * x70,
            x46 * x65 * x9,
            x10 * x51 * x55,
            x51 * x62 * x9,
            x52 * x55 * x9,
            x38 * x58 * x76,
            x48 * x63 * x76,
            x38 * x63 * x80,
            x39 * x57 * x76,
            x40 * x56 * x76,
            x34 * x80 * x82,
            x38 * x57 * x84,
            x29 * x82 * x84,
            x2 * x88,
            x10 * x47 * x76,
            x49 * x76 * x9,
            x47 * x80 * x9,
            x34 * x5 * x81 * x84,
            x40 * x84 * x9,
            x34 * x88,
            x5 * x90,
            x29 * x90,
            x81
            * (
                x3 * (x15 * (x50 + x77) + 4 * x44 * x74 + x79 + x83 + x89 * (x44 + x85))
                + x43 * x87
            ),
        ]
    )


def diag_quadrupole3d_22(a, A, b, B, C):
    """Cartesian 3D (dd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x17 = 2 * x3
    x18 = x4 * x8
    x19 = x17 * x18
    x20 = x10 + x9
    x21 = x12 * x20
    x22 = x19 + x21
    x23 = x12 * x22
    x24 = x16 + x23
    x25 = x10 + x14
    x26 = x12 * x25
    x27 = x3 * (x13 + x18)
    x28 = 4 * x10 * x4 + 2 * x27
    x29 = x2 * x24 + x3 * (2 * x21 + 2 * x26 + x28)
    x30 = x12 ** 2 * x8
    x31 = x2 * x22
    x32 = 2 * x2
    x33 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x34 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x35 = numpy.pi * x0 * x34
    x36 = x33 * x35
    x37 = -x0 * (a * A[1] + b * B[1])
    x38 = -x37 - B[1]
    x39 = x16 + x31
    x40 = x2 * x20
    x41 = x36 * (x2 * x39 + x3 * (x21 + x25 * x32 + x28 + x40))
    x42 = -x0 * (a * A[2] + b * B[2])
    x43 = -x42 - B[2]
    x44 = x33 * x7
    x45 = x3 * x44
    x46 = x38 ** 2 * x44
    x47 = x45 + x46
    x48 = x19 + x40
    x49 = x2 * x48 + x3 * (x11 + x18 * x32 + x9)
    x50 = x34 * x7
    x51 = x36 * x43
    x52 = x3 * x50
    x53 = x43 ** 2 * x50
    x54 = x52 + x53
    x55 = -x37 - A[1]
    x56 = x29 * x36
    x57 = x44 * x55
    x58 = x38 * x57 + x45
    x59 = x38 * x44
    x60 = x17 * x59 + x47 * x55
    x61 = x43 * x50
    x62 = -x42 - A[2]
    x63 = x36 * x62
    x64 = x50 * x62
    x65 = x43 * x64 + x52
    x66 = x17 * x61 + x54 * x62
    x67 = x44 * x55 ** 2 + x45
    x68 = x3 * (x57 + x59) + x55 * x58
    x69 = 2 * x55
    x70 = 3 * x45
    x71 = x46 + x70
    x72 = x3 * (x59 * x69 + x71) + x55 * x60
    x73 = x50 * x62 ** 2 + x52
    x74 = x3 * (x61 + x64) + x62 * x65
    x75 = 2 * x62
    x76 = 3 * x52
    x77 = x53 + x76
    x78 = x3 * (x61 * x75 + x77) + x62 * x66
    x79 = -x37 - C[1]
    x80 = x44 * x79 ** 2
    x81 = x45 + x80
    x82 = x10 + x30
    x83 = x13 * x17 + x2 * x82
    x84 = x2 * x83 + x3 * (x11 + x13 * x32 + x30)
    x85 = x44 * x79
    x86 = x17 * x85
    x87 = x38 * x81
    x88 = x86 + x87
    x89 = x2 * x8
    x90 = x10 + x13 * x2
    x91 = x2 * x90 + x3 * (x13 + x89)
    x92 = x59 * x79
    x93 = 2 * x92
    x94 = x70 + x80
    x95 = x3 * (x93 + x94)
    x96 = x38 * x88
    x97 = x95 + x96
    x98 = x10 + x2 ** 2 * x8
    x99 = x55 * x81
    x100 = x86 + x99
    x101 = x55 * x88
    x102 = x101 + x95
    x103 = x45 + x92
    x104 = x103 * x38
    x105 = x3 * (x59 + x85)
    x106 = 2 * x105 + 4 * x45 * x79
    x107 = x3 * (2 * x104 + x106 + 2 * x87) + x55 * x97
    x108 = x35 * x6
    x109 = x107 * x108
    x110 = x108 * x2
    x111 = x100 * x55 + x3 * (x69 * x85 + x94)
    x112 = x108 * (x102 * x55 + x3 * (x103 * x69 + x106 + x87 + x99))
    x113 = x108 * x12
    x114 = -x42 - C[2]
    x115 = x114 ** 2 * x50
    x116 = x115 + x52
    x117 = x114 * x50
    x118 = x117 * x17
    x119 = x116 * x43
    x120 = x118 + x119
    x121 = x114 * x61
    x122 = 2 * x121
    x123 = x115 + x76
    x124 = x3 * (x122 + x123)
    x125 = x120 * x43
    x126 = x124 + x125
    x127 = numpy.pi * x0 * x33 * x6
    x128 = x127 * x2
    x129 = x116 * x62
    x130 = x118 + x129
    x131 = x120 * x62
    x132 = x124 + x131
    x133 = x121 + x52
    x134 = x133 * x43
    x135 = x3 * (x117 + x61)
    x136 = 4 * x114 * x52 + 2 * x135
    x137 = x126 * x62 + x3 * (2 * x119 + 2 * x134 + x136)
    x138 = x127 * x137
    x139 = x12 * x127
    x140 = x130 * x62 + x3 * (x117 * x75 + x123)
    x141 = x127 * (x132 * x62 + x3 * (x119 + x129 + x133 * x75 + x136))

    # 108 item(s)
    return numpy.array(
        [
            x36
            * (
                x2 * x29
                + x3 * (3 * x16 + x17 * (x15 + x30) + x23 + 2 * x31 + x32 * (x26 + x27))
            ),
            x38 * x41,
            x41 * x43,
            x47 * x49 * x50,
            x38 * x49 * x51,
            x44 * x49 * x54,
            x55 * x56,
            x39 * x50 * x58,
            x39 * x51 * x55,
            x48 * x50 * x60,
            x48 * x58 * x61,
            x48 * x54 * x57,
            x56 * x62,
            x38 * x39 * x63,
            x39 * x44 * x65,
            x47 * x48 * x64,
            x48 * x59 * x65,
            x44 * x48 * x66,
            x24 * x50 * x67,
            x22 * x50 * x68,
            x22 * x61 * x67,
            x20 * x50 * x72,
            x20 * x61 * x68,
            x20 * x54 * x67,
            x24 * x55 * x63,
            x22 * x58 * x64,
            x22 * x57 * x65,
            x20 * x60 * x64,
            x20 * x58 * x65,
            x20 * x57 * x66,
            x24 * x44 * x73,
            x22 * x59 * x73,
            x22 * x44 * x74,
            x20 * x47 * x73,
            x20 * x59 * x74,
            x20 * x44 * x78,
            x50 * x81 * x84,
            x50 * x88 * x91,
            x61 * x81 * x91,
            x50 * x97 * x98,
            x61 * x88 * x98,
            x54 * x81 * x98,
            x100 * x50 * x83,
            x102 * x50 * x90,
            x100 * x61 * x90,
            x109 * x2,
            x102 * x110 * x43,
            x100 * x54 * x89,
            x64 * x81 * x83,
            x64 * x88 * x90,
            x65 * x81 * x90,
            x110 * x62 * x97,
            x65 * x88 * x89,
            x66 * x81 * x89,
            x111 * x50 * x82,
            x112 * x12,
            x111 * x113 * x43,
            x108
            * (
                x107 * x55
                + x3
                * (2 * x101 + x17 * (x71 + x93) + x69 * (x104 + x105) + 3 * x95 + x96)
            ),
            x112 * x43,
            x111 * x54 * x8,
            x100 * x64 * x82,
            x102 * x113 * x62,
            x100 * x13 * x65,
            x109 * x62,
            x102 * x65 * x8,
            x100 * x66 * x8,
            x73 * x81 * x82,
            x13 * x73 * x88,
            x13 * x74 * x81,
            x73 * x8 * x97,
            x74 * x8 * x88,
            x78 * x8 * x81,
            x116 * x44 * x84,
            x116 * x59 * x91,
            x120 * x44 * x91,
            x116 * x47 * x98,
            x120 * x59 * x98,
            x126 * x44 * x98,
            x116 * x57 * x83,
            x116 * x58 * x90,
            x120 * x57 * x90,
            x116 * x60 * x89,
            x120 * x58 * x89,
            x126 * x128 * x55,
            x130 * x44 * x83,
            x130 * x59 * x90,
            x132 * x44 * x90,
            x130 * x47 * x89,
            x128 * x132 * x38,
            x138 * x2,
            x116 * x67 * x82,
            x116 * x13 * x68,
            x120 * x13 * x67,
            x116 * x72 * x8,
            x120 * x68 * x8,
            x126 * x67 * x8,
            x130 * x57 * x82,
            x13 * x130 * x58,
            x132 * x139 * x55,
            x130 * x60 * x8,
            x132 * x58 * x8,
            x138 * x55,
            x140 * x44 * x82,
            x139 * x140 * x38,
            x12 * x141,
            x140 * x47 * x8,
            x141 * x38,
            x127
            * (
                x137 * x62
                + x3
                * (3 * x124 + x125 + 2 * x131 + x17 * (x122 + x77) + x75 * (x134 + x135))
            ),
        ]
    )


def diag_quadrupole3d_23(a, A, b, B, C):
    """Cartesian 3D (df) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x30 = x3 * (x21 + x28)
    x31 = 2 * x29 + 3 * x30
    x32 = x25 + x30
    x33 = x32 * x4
    x34 = 4 * x10 * x13 + 2 * x12
    x35 = x3 * (2 * x16 + 2 * x23 + x34)
    x36 = x33 + x35
    x37 = x2 * x36 + x3 * (2 * x18 + 3 * x25 + x31)
    x38 = x2 * x32
    x39 = x19 * x9
    x40 = x13 + x26
    x41 = x4 * x40
    x42 = x39 + x41
    x43 = 2 * x2
    x44 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x45 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x46 = numpy.pi * x0 * x45
    x47 = x44 * x46
    x48 = -x0 * (a * A[1] + b * B[1])
    x49 = -x48 - B[1]
    x50 = x35 + x38
    x51 = x2 * x24
    x52 = x47 * (x2 * x50 + x3 * (x17 * x43 + x25 + x31 + 2 * x51))
    x53 = -x0 * (a * A[2] + b * B[2])
    x54 = -x53 - B[2]
    x55 = x44 * x7
    x56 = x3 * x55
    x57 = x49 ** 2 * x55
    x58 = x56 + x57
    x59 = x30 + x51
    x60 = x2 * x22
    x61 = x2 * x59 + x3 * (x15 * x43 + x23 + x34 + x60)
    x62 = x45 * x7
    x63 = x47 * x54
    x64 = x3 * x62
    x65 = x54 ** 2 * x62
    x66 = x64 + x65
    x67 = x49 * x55
    x68 = x19 * x67
    x69 = x49 * x58
    x70 = x68 + x69
    x71 = x20 + x60
    x72 = x2 * x71 + x3 * (x11 * x43 + x21 + x27)
    x73 = x54 * x62
    x74 = x19 * x73
    x75 = x54 * x66
    x76 = x74 + x75
    x77 = -x48 - A[1]
    x78 = x37 * x47
    x79 = x55 * x77
    x80 = x49 * x79 + x56
    x81 = x58 * x77
    x82 = x68 + x81
    x83 = 3 * x56
    x84 = x3 * (3 * x57 + x83) + x70 * x77
    x85 = -x53 - A[2]
    x86 = x47 * x85
    x87 = x62 * x85
    x88 = x54 * x87 + x64
    x89 = x66 * x85
    x90 = x74 + x89
    x91 = 3 * x64
    x92 = x3 * (3 * x65 + x91) + x76 * x85
    x93 = x55 * x77 ** 2 + x56
    x94 = x3 * (x67 + x79) + x77 * x80
    x95 = 2 * x77
    x96 = x57 + x83
    x97 = x3 * (x67 * x95 + x96) + x77 * x82
    x98 = x3 * (8 * x49 * x56 + x69 + 3 * x81) + x77 * x84
    x99 = x62 * x85 ** 2 + x64
    x100 = x3 * (x73 + x87) + x85 * x88
    x101 = 2 * x85
    x102 = x65 + x91
    x103 = x3 * (x101 * x73 + x102) + x85 * x90
    x104 = x3 * (8 * x54 * x64 + x75 + 3 * x89) + x85 * x92
    x105 = -x48 - C[1]
    x106 = x105 ** 2 * x55
    x107 = x106 + x56
    x108 = x2 * x40
    x109 = x2 * x42 + x3 * (3 * x26 + x27)
    x110 = x109 * x2 + x3 * (3 * x108 + 8 * x13 * x4 + x41)
    x111 = x105 * x55
    x112 = x111 * x19
    x113 = x107 * x49
    x114 = x112 + x113
    x115 = x108 + x39
    x116 = x115 * x2 + x3 * (x26 + x27 + x43 * x9)
    x117 = x105 * x67
    x118 = 2 * x117
    x119 = x106 + x83
    x120 = x3 * (x118 + x119)
    x121 = x114 * x49
    x122 = x120 + x121
    x123 = x2 * x8
    x124 = x13 + x2 * x9
    x125 = x124 * x2 + x3 * (x123 + x9)
    x126 = x122 * x49
    x127 = x117 + x56
    x128 = x127 * x49
    x129 = x3 * (x111 + x67)
    x130 = 4 * x105 * x56 + 2 * x129
    x131 = x3 * (2 * x113 + 2 * x128 + x130)
    x132 = x126 + x131
    x133 = x13 + x2 ** 2 * x8
    x134 = x107 * x77
    x135 = x112 + x134
    x136 = x114 * x77
    x137 = x120 + x136
    x138 = x122 * x77
    x139 = x131 + x138
    x140 = x128 + x129
    x141 = x140 * x49
    x142 = x3 * (x118 + x96)
    x143 = 3 * x120 + 2 * x142
    x144 = x132 * x77 + x3 * (3 * x121 + 2 * x141 + x143)
    x145 = x46 * x6
    x146 = x144 * x145
    x147 = x145 * x2
    x148 = x135 * x77 + x3 * (x111 * x95 + x119)
    x149 = x137 * x77 + x3 * (x113 + x127 * x95 + x130 + x134)
    x150 = x145 * (x139 * x77 + x3 * (x121 + 2 * x136 + x140 * x95 + x143))
    x151 = x145 * x4
    x152 = -x53 - C[2]
    x153 = x152 ** 2 * x62
    x154 = x153 + x64
    x155 = x152 * x62
    x156 = x155 * x19
    x157 = x154 * x54
    x158 = x156 + x157
    x159 = x152 * x73
    x160 = 2 * x159
    x161 = x153 + x91
    x162 = x3 * (x160 + x161)
    x163 = x158 * x54
    x164 = x162 + x163
    x165 = x164 * x54
    x166 = x159 + x64
    x167 = x166 * x54
    x168 = x3 * (x155 + x73)
    x169 = 4 * x152 * x64 + 2 * x168
    x170 = x3 * (2 * x157 + 2 * x167 + x169)
    x171 = x165 + x170
    x172 = numpy.pi * x0 * x44 * x6
    x173 = x172 * x2
    x174 = x154 * x85
    x175 = x156 + x174
    x176 = x158 * x85
    x177 = x162 + x176
    x178 = x164 * x85
    x179 = x170 + x178
    x180 = x167 + x168
    x181 = x180 * x54
    x182 = x3 * (x102 + x160)
    x183 = 3 * x162 + 2 * x182
    x184 = x171 * x85 + x3 * (3 * x163 + 2 * x181 + x183)
    x185 = x172 * x184
    x186 = x172 * x4
    x187 = x175 * x85 + x3 * (x101 * x155 + x161)
    x188 = x177 * x85 + x3 * (x101 * x166 + x157 + x169 + x174)
    x189 = x172 * (x179 * x85 + x3 * (x101 * x180 + x163 + 2 * x176 + x183))

    # 180 item(s)
    return numpy.array(
        [
            x47
            * (
                x2 * x37
                + x3
                * (
                    x19 * (3 * x12 + 3 * x16 + x42)
                    + x33
                    + 4 * x35
                    + 3 * x38
                    + x43 * (x18 + x29)
                )
            ),
            x49 * x52,
            x52 * x54,
            x58 * x61 * x62,
            x49 * x61 * x63,
            x55 * x61 * x66,
            x62 * x70 * x72,
            x58 * x72 * x73,
            x66 * x67 * x72,
            x55 * x72 * x76,
            x77 * x78,
            x50 * x62 * x80,
            x50 * x63 * x77,
            x59 * x62 * x82,
            x59 * x73 * x80,
            x59 * x66 * x79,
            x62 * x71 * x84,
            x71 * x73 * x82,
            x66 * x71 * x80,
            x71 * x76 * x79,
            x78 * x85,
            x49 * x50 * x86,
            x50 * x55 * x88,
            x58 * x59 * x87,
            x59 * x67 * x88,
            x55 * x59 * x90,
            x70 * x71 * x87,
            x58 * x71 * x88,
            x67 * x71 * x90,
            x55 * x71 * x92,
            x36 * x62 * x93,
            x32 * x62 * x94,
            x32 * x73 * x93,
            x24 * x62 * x97,
            x24 * x73 * x94,
            x24 * x66 * x93,
            x22 * x62 * x98,
            x22 * x73 * x97,
            x22 * x66 * x94,
            x22 * x76 * x93,
            x36 * x77 * x86,
            x32 * x80 * x87,
            x32 * x79 * x88,
            x24 * x82 * x87,
            x24 * x80 * x88,
            x24 * x79 * x90,
            x22 * x84 * x87,
            x22 * x82 * x88,
            x22 * x80 * x90,
            x22 * x79 * x92,
            x36 * x55 * x99,
            x32 * x67 * x99,
            x100 * x32 * x55,
            x24 * x58 * x99,
            x100 * x24 * x67,
            x103 * x24 * x55,
            x22 * x70 * x99,
            x100 * x22 * x58,
            x103 * x22 * x67,
            x104 * x22 * x55,
            x107 * x110 * x62,
            x114 * x116 * x62,
            x107 * x116 * x73,
            x122 * x125 * x62,
            x114 * x125 * x73,
            x107 * x125 * x66,
            x132 * x133 * x62,
            x122 * x133 * x73,
            x114 * x133 * x66,
            x107 * x133 * x76,
            x109 * x135 * x62,
            x115 * x137 * x62,
            x115 * x135 * x73,
            x124 * x139 * x62,
            x124 * x137 * x73,
            x124 * x135 * x66,
            x146 * x2,
            x139 * x147 * x54,
            x123 * x137 * x66,
            x123 * x135 * x76,
            x107 * x109 * x87,
            x114 * x115 * x87,
            x107 * x115 * x88,
            x122 * x124 * x87,
            x114 * x124 * x88,
            x107 * x124 * x90,
            x132 * x147 * x85,
            x122 * x123 * x88,
            x114 * x123 * x90,
            x107 * x123 * x92,
            x148 * x42 * x62,
            x149 * x40 * x62,
            x148 * x40 * x73,
            x150 * x4,
            x149 * x151 * x54,
            x148 * x66 * x9,
            x145
            * (
                x144 * x77
                + x3
                * (
                    x126
                    + 4 * x131
                    + 3 * x138
                    + x19 * (3 * x128 + 3 * x129 + x70)
                    + x95 * (x141 + x142)
                )
            ),
            x150 * x54,
            x149 * x66 * x8,
            x148 * x76 * x8,
            x135 * x42 * x87,
            x137 * x40 * x87,
            x135 * x40 * x88,
            x139 * x151 * x85,
            x137 * x88 * x9,
            x135 * x9 * x90,
            x146 * x85,
            x139 * x8 * x88,
            x137 * x8 * x90,
            x135 * x8 * x92,
            x107 * x42 * x99,
            x114 * x40 * x99,
            x100 * x107 * x40,
            x122 * x9 * x99,
            x100 * x114 * x9,
            x103 * x107 * x9,
            x132 * x8 * x99,
            x100 * x122 * x8,
            x103 * x114 * x8,
            x104 * x107 * x8,
            x110 * x154 * x55,
            x116 * x154 * x67,
            x116 * x158 * x55,
            x125 * x154 * x58,
            x125 * x158 * x67,
            x125 * x164 * x55,
            x133 * x154 * x70,
            x133 * x158 * x58,
            x133 * x164 * x67,
            x133 * x171 * x55,
            x109 * x154 * x79,
            x115 * x154 * x80,
            x115 * x158 * x79,
            x124 * x154 * x82,
            x124 * x158 * x80,
            x124 * x164 * x79,
            x123 * x154 * x84,
            x123 * x158 * x82,
            x123 * x164 * x80,
            x171 * x173 * x77,
            x109 * x175 * x55,
            x115 * x175 * x67,
            x115 * x177 * x55,
            x124 * x175 * x58,
            x124 * x177 * x67,
            x124 * x179 * x55,
            x123 * x175 * x70,
            x123 * x177 * x58,
            x173 * x179 * x49,
            x185 * x2,
            x154 * x42 * x93,
            x154 * x40 * x94,
            x158 * x40 * x93,
            x154 * x9 * x97,
            x158 * x9 * x94,
            x164 * x9 * x93,
            x154 * x8 * x98,
            x158 * x8 * x97,
            x164 * x8 * x94,
            x171 * x8 * x93,
            x175 * x42 * x79,
            x175 * x40 * x80,
            x177 * x40 * x79,
            x175 * x82 * x9,
            x177 * x80 * x9,
            x179 * x186 * x77,
            x175 * x8 * x84,
            x177 * x8 * x82,
            x179 * x8 * x80,
            x185 * x77,
            x187 * x42 * x55,
            x187 * x40 * x67,
            x188 * x40 * x55,
            x187 * x58 * x9,
            x186 * x188 * x49,
            x189 * x4,
            x187 * x70 * x8,
            x188 * x58 * x8,
            x189 * x49,
            x172
            * (
                x184 * x85
                + x3
                * (
                    x101 * (x181 + x182)
                    + x165
                    + 4 * x170
                    + 3 * x178
                    + x19 * (3 * x167 + 3 * x168 + x76)
                )
            ),
        ]
    )


def diag_quadrupole3d_24(a, A, b, B, C):
    """Cartesian 3D (dg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x20 = 2 * x3
    x21 = x12 * x20
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
    x36 = 4 * x12 * x3 + 2 * x13
    x37 = x3 * (2 * x17 + 2 * x24 + x36)
    x38 = x35 + x37
    x39 = x38 * x4
    x40 = x33 + x39
    x41 = x19 + x30
    x42 = x4 * x41
    x43 = x20 * x9
    x44 = x14 + x27
    x45 = x4 * x44
    x46 = x43 + x45
    x47 = x3 * (3 * x13 + 3 * x17 + x46)
    x48 = 4 * x37 + 2 * x47
    x49 = x2 * x40 + x3 * (4 * x35 + 2 * x42 + x48)
    x50 = x3 * (3 * x27 + x28)
    x51 = x4 * x46
    x52 = x50 + x51
    x53 = x2 * x38
    x54 = 2 * x2
    x55 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x56 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x57 = numpy.pi * x0 * x56
    x58 = x55 * x57
    x59 = -x0 * (a * A[1] + b * B[1])
    x60 = -x59 - B[1]
    x61 = x33 + x53
    x62 = x2 * x34
    x63 = x58 * (x2 * x61 + x3 * (x35 + x41 * x54 + x48 + 3 * x62))
    x64 = -x0 * (a * A[2] + b * B[2])
    x65 = -x64 - B[2]
    x66 = x5 * x55
    x67 = x3 * x66
    x68 = x60 ** 2 * x66
    x69 = x67 + x68
    x70 = x37 + x62
    x71 = x2 * x25
    x72 = x2 * x70 + x3 * (x18 * x54 + x26 + x32 + 2 * x71)
    x73 = x5 * x56
    x74 = x58 * x65
    x75 = x3 * x73
    x76 = x65 ** 2 * x73
    x77 = x75 + x76
    x78 = x60 * x66
    x79 = x20 * x78
    x80 = x60 * x69
    x81 = x79 + x80
    x82 = x31 + x71
    x83 = x2 * x23
    x84 = x2 * x82 + x3 * (x16 * x54 + x24 + x36 + x83)
    x85 = x65 * x73
    x86 = x20 * x85
    x87 = x65 * x77
    x88 = x86 + x87
    x89 = 3 * x67
    x90 = x3 * (3 * x68 + x89)
    x91 = x60 * x81
    x92 = x90 + x91
    x93 = x21 + x83
    x94 = x2 * x93 + x3 * (x12 * x54 + x22 + x28)
    x95 = 3 * x75
    x96 = x3 * (3 * x76 + x95)
    x97 = x65 * x88
    x98 = x96 + x97
    x99 = -x59 - A[1]
    x100 = x49 * x58
    x101 = x66 * x99
    x102 = x101 * x60 + x67
    x103 = x69 * x99
    x104 = x103 + x79
    x105 = x81 * x99
    x106 = x105 + x90
    x107 = 8 * x60 * x67
    x108 = x3 * (x107 + 4 * x80) + x92 * x99
    x109 = -x64 - A[2]
    x110 = x109 * x58
    x111 = x109 * x73
    x112 = x111 * x65 + x75
    x113 = x109 * x77
    x114 = x113 + x86
    x115 = x109 * x88
    x116 = x115 + x96
    x117 = 8 * x65 * x75
    x118 = x109 * x98 + x3 * (x117 + 4 * x87)
    x119 = x66 * x99 ** 2 + x67
    x120 = x102 * x99 + x3 * (x101 + x78)
    x121 = 2 * x99
    x122 = x68 + x89
    x123 = x104 * x99 + x3 * (x121 * x78 + x122)
    x124 = x106 * x99 + x3 * (3 * x103 + x107 + x80)
    x125 = x108 * x99 + x3 * (4 * x105 + 5 * x90 + x91)
    x126 = x109 ** 2 * x73 + x75
    x127 = x109 * x112 + x3 * (x111 + x85)
    x128 = 2 * x109
    x129 = x76 + x95
    x130 = x109 * x114 + x3 * (x128 * x85 + x129)
    x131 = x109 * x116 + x3 * (3 * x113 + x117 + x87)
    x132 = x109 * x118 + x3 * (4 * x115 + 5 * x96 + x97)
    x133 = -x59 - C[1]
    x134 = x133 ** 2 * x66
    x135 = x134 + x67
    x136 = x2 * x46
    x137 = 8 * x3 * x9
    x138 = x2 * x52 + x3 * (x137 + 4 * x45)
    x139 = x138 * x2 + x3 * (4 * x136 + 5 * x50 + x51)
    x140 = x133 * x66
    x141 = x140 * x20
    x142 = x135 * x60
    x143 = x141 + x142
    x144 = x2 * x44
    x145 = x136 + x50
    x146 = x145 * x2 + x3 * (x137 + 3 * x144 + x45)
    x147 = x133 * x78
    x148 = 2 * x147
    x149 = x134 + x89
    x150 = x3 * (x148 + x149)
    x151 = x143 * x60
    x152 = x150 + x151
    x153 = x144 + x43
    x154 = x153 * x2 + x3 * (x27 + x28 + x54 * x9)
    x155 = x152 * x60
    x156 = x147 + x67
    x157 = x156 * x60
    x158 = x3 * (x140 + x78)
    x159 = 4 * x133 * x67 + 2 * x158
    x160 = x3 * (2 * x142 + 2 * x157 + x159)
    x161 = x155 + x160
    x162 = x11 * x2
    x163 = x14 + x2 * x9
    x164 = x163 * x2 + x3 * (x162 + x9)
    x165 = x157 + x158
    x166 = x165 * x60
    x167 = x3 * (x122 + x148)
    x168 = 3 * x150 + 2 * x167
    x169 = x3 * (3 * x151 + 2 * x166 + x168)
    x170 = x161 * x60
    x171 = x169 + x170
    x172 = x11 * x2 ** 2 + x14
    x173 = x135 * x99
    x174 = x141 + x173
    x175 = x143 * x99
    x176 = x150 + x175
    x177 = x152 * x99
    x178 = x160 + x177
    x179 = x161 * x99
    x180 = x169 + x179
    x181 = x166 + x167
    x182 = x181 * x60
    x183 = x3 * (3 * x157 + 3 * x158 + x81)
    x184 = 4 * x160 + 2 * x183
    x185 = x171 * x99 + x3 * (4 * x155 + 2 * x182 + x184)
    x186 = x57 * x7
    x187 = x185 * x186
    x188 = x186 * x65
    x189 = x174 * x99 + x3 * (x121 * x140 + x149)
    x190 = x176 * x99 + x3 * (x121 * x156 + x142 + x159 + x173)
    x191 = x178 * x99 + x3 * (x121 * x165 + x151 + x168 + 2 * x175)
    x192 = x180 * x99 + x3 * (x121 * x181 + x155 + 3 * x177 + x184)
    x193 = x57 * x8
    x194 = -x64 - C[2]
    x195 = x194 ** 2 * x73
    x196 = x195 + x75
    x197 = x194 * x73
    x198 = x197 * x20
    x199 = x196 * x65
    x200 = x198 + x199
    x201 = x194 * x85
    x202 = 2 * x201
    x203 = x195 + x95
    x204 = x3 * (x202 + x203)
    x205 = x200 * x65
    x206 = x204 + x205
    x207 = x206 * x65
    x208 = x201 + x75
    x209 = x208 * x65
    x210 = x3 * (x197 + x85)
    x211 = 4 * x194 * x75 + 2 * x210
    x212 = x3 * (2 * x199 + 2 * x209 + x211)
    x213 = x207 + x212
    x214 = x209 + x210
    x215 = x214 * x65
    x216 = x3 * (x129 + x202)
    x217 = 3 * x204 + 2 * x216
    x218 = x3 * (3 * x205 + 2 * x215 + x217)
    x219 = x213 * x65
    x220 = x218 + x219
    x221 = numpy.pi * x0 * x55
    x222 = x221 * x7
    x223 = x109 * x196
    x224 = x198 + x223
    x225 = x109 * x200
    x226 = x204 + x225
    x227 = x109 * x206
    x228 = x212 + x227
    x229 = x109 * x213
    x230 = x218 + x229
    x231 = x222 * x60
    x232 = x215 + x216
    x233 = x232 * x65
    x234 = x3 * (3 * x209 + 3 * x210 + x88)
    x235 = 4 * x212 + 2 * x234
    x236 = x109 * x220 + x3 * (4 * x207 + 2 * x233 + x235)
    x237 = x222 * x236
    x238 = x221 * x8
    x239 = x109 * x224 + x3 * (x128 * x197 + x203)
    x240 = x109 * x226 + x3 * (x128 * x208 + x199 + x211 + x223)
    x241 = x109 * x228 + x3 * (x128 * x214 + x205 + x217 + 2 * x225)
    x242 = x109 * x230 + x3 * (x128 * x232 + x207 + 3 * x227 + x235)

    # 270 item(s)
    return numpy.array(
        [
            x58
            * (
                x2 * x49
                + x3
                * (
                    x20 * (4 * x19 + 4 * x30 + x52)
                    + 5 * x33
                    + x39
                    + 4 * x53
                    + x54 * (x42 + x47)
                )
            ),
            x60 * x63,
            x63 * x65,
            x69 * x72 * x73,
            x60 * x72 * x74,
            x66 * x72 * x77,
            x73 * x81 * x84,
            x69 * x84 * x85,
            x77 * x78 * x84,
            x66 * x84 * x88,
            x73 * x92 * x94,
            x81 * x85 * x94,
            x69 * x77 * x94,
            x78 * x88 * x94,
            x66 * x94 * x98,
            x100 * x99,
            x102 * x61 * x73,
            x61 * x74 * x99,
            x104 * x70 * x73,
            x102 * x70 * x85,
            x101 * x70 * x77,
            x106 * x73 * x82,
            x104 * x82 * x85,
            x102 * x77 * x82,
            x101 * x82 * x88,
            x108 * x73 * x93,
            x106 * x85 * x93,
            x104 * x77 * x93,
            x102 * x88 * x93,
            x101 * x93 * x98,
            x100 * x109,
            x110 * x60 * x61,
            x112 * x61 * x66,
            x111 * x69 * x70,
            x112 * x70 * x78,
            x114 * x66 * x70,
            x111 * x81 * x82,
            x112 * x69 * x82,
            x114 * x78 * x82,
            x116 * x66 * x82,
            x111 * x92 * x93,
            x112 * x81 * x93,
            x114 * x69 * x93,
            x116 * x78 * x93,
            x118 * x66 * x93,
            x119 * x40 * x73,
            x120 * x38 * x73,
            x119 * x38 * x85,
            x123 * x34 * x73,
            x120 * x34 * x85,
            x119 * x34 * x77,
            x124 * x25 * x73,
            x123 * x25 * x85,
            x120 * x25 * x77,
            x119 * x25 * x88,
            x125 * x23 * x73,
            x124 * x23 * x85,
            x123 * x23 * x77,
            x120 * x23 * x88,
            x119 * x23 * x98,
            x110 * x40 * x99,
            x102 * x111 * x38,
            x101 * x112 * x38,
            x104 * x111 * x34,
            x102 * x112 * x34,
            x101 * x114 * x34,
            x106 * x111 * x25,
            x104 * x112 * x25,
            x102 * x114 * x25,
            x101 * x116 * x25,
            x108 * x111 * x23,
            x106 * x112 * x23,
            x104 * x114 * x23,
            x102 * x116 * x23,
            x101 * x118 * x23,
            x126 * x40 * x66,
            x126 * x38 * x78,
            x127 * x38 * x66,
            x126 * x34 * x69,
            x127 * x34 * x78,
            x130 * x34 * x66,
            x126 * x25 * x81,
            x127 * x25 * x69,
            x130 * x25 * x78,
            x131 * x25 * x66,
            x126 * x23 * x92,
            x127 * x23 * x81,
            x130 * x23 * x69,
            x131 * x23 * x78,
            x132 * x23 * x66,
            x135 * x139 * x73,
            x143 * x146 * x73,
            x135 * x146 * x85,
            x152 * x154 * x73,
            x143 * x154 * x85,
            x135 * x154 * x77,
            x161 * x164 * x73,
            x152 * x164 * x85,
            x143 * x164 * x77,
            x135 * x164 * x88,
            x171 * x172 * x73,
            x161 * x172 * x85,
            x152 * x172 * x77,
            x143 * x172 * x88,
            x135 * x172 * x98,
            x138 * x174 * x73,
            x145 * x176 * x73,
            x145 * x174 * x85,
            x153 * x178 * x73,
            x153 * x176 * x85,
            x153 * x174 * x77,
            x163 * x180 * x73,
            x163 * x178 * x85,
            x163 * x176 * x77,
            x163 * x174 * x88,
            x187 * x2,
            x180 * x188 * x2,
            x162 * x178 * x77,
            x162 * x176 * x88,
            x162 * x174 * x98,
            x111 * x135 * x138,
            x111 * x143 * x145,
            x112 * x135 * x145,
            x111 * x152 * x153,
            x112 * x143 * x153,
            x114 * x135 * x153,
            x111 * x161 * x163,
            x112 * x152 * x163,
            x114 * x143 * x163,
            x116 * x135 * x163,
            x109 * x171 * x186 * x2,
            x112 * x161 * x162,
            x114 * x152 * x162,
            x116 * x143 * x162,
            x118 * x135 * x162,
            x189 * x52 * x73,
            x190 * x46 * x73,
            x189 * x46 * x85,
            x191 * x44 * x73,
            x190 * x44 * x85,
            x189 * x44 * x77,
            x192 * x193,
            x191 * x193 * x65,
            x190 * x77 * x9,
            x189 * x88 * x9,
            x186
            * (
                x185 * x99
                + x3
                * (
                    x121 * (x182 + x183)
                    + 5 * x169
                    + x170
                    + 4 * x179
                    + x20 * (4 * x166 + 4 * x167 + x92)
                )
            ),
            x188 * x192,
            x11 * x191 * x77,
            x11 * x190 * x88,
            x11 * x189 * x98,
            x111 * x174 * x52,
            x111 * x176 * x46,
            x112 * x174 * x46,
            x111 * x178 * x44,
            x112 * x176 * x44,
            x114 * x174 * x44,
            x109 * x180 * x193,
            x112 * x178 * x9,
            x114 * x176 * x9,
            x116 * x174 * x9,
            x109 * x187,
            x11 * x112 * x180,
            x11 * x114 * x178,
            x11 * x116 * x176,
            x11 * x118 * x174,
            x126 * x135 * x52,
            x126 * x143 * x46,
            x127 * x135 * x46,
            x126 * x152 * x44,
            x127 * x143 * x44,
            x130 * x135 * x44,
            x126 * x161 * x9,
            x127 * x152 * x9,
            x130 * x143 * x9,
            x131 * x135 * x9,
            x11 * x126 * x171,
            x11 * x127 * x161,
            x11 * x130 * x152,
            x11 * x131 * x143,
            x11 * x132 * x135,
            x139 * x196 * x66,
            x146 * x196 * x78,
            x146 * x200 * x66,
            x154 * x196 * x69,
            x154 * x200 * x78,
            x154 * x206 * x66,
            x164 * x196 * x81,
            x164 * x200 * x69,
            x164 * x206 * x78,
            x164 * x213 * x66,
            x172 * x196 * x92,
            x172 * x200 * x81,
            x172 * x206 * x69,
            x172 * x213 * x78,
            x172 * x220 * x66,
            x101 * x138 * x196,
            x102 * x145 * x196,
            x101 * x145 * x200,
            x104 * x153 * x196,
            x102 * x153 * x200,
            x101 * x153 * x206,
            x106 * x163 * x196,
            x104 * x163 * x200,
            x102 * x163 * x206,
            x101 * x163 * x213,
            x108 * x162 * x196,
            x106 * x162 * x200,
            x104 * x162 * x206,
            x102 * x162 * x213,
            x2 * x220 * x222 * x99,
            x138 * x224 * x66,
            x145 * x224 * x78,
            x145 * x226 * x66,
            x153 * x224 * x69,
            x153 * x226 * x78,
            x153 * x228 * x66,
            x163 * x224 * x81,
            x163 * x226 * x69,
            x163 * x228 * x78,
            x163 * x230 * x66,
            x162 * x224 * x92,
            x162 * x226 * x81,
            x162 * x228 * x69,
            x2 * x230 * x231,
            x2 * x237,
            x119 * x196 * x52,
            x120 * x196 * x46,
            x119 * x200 * x46,
            x123 * x196 * x44,
            x120 * x200 * x44,
            x119 * x206 * x44,
            x124 * x196 * x9,
            x123 * x200 * x9,
            x120 * x206 * x9,
            x119 * x213 * x9,
            x11 * x125 * x196,
            x11 * x124 * x200,
            x11 * x123 * x206,
            x11 * x120 * x213,
            x11 * x119 * x220,
            x101 * x224 * x52,
            x102 * x224 * x46,
            x101 * x226 * x46,
            x104 * x224 * x44,
            x102 * x226 * x44,
            x101 * x228 * x44,
            x106 * x224 * x9,
            x104 * x226 * x9,
            x102 * x228 * x9,
            x230 * x238 * x99,
            x108 * x11 * x224,
            x106 * x11 * x226,
            x104 * x11 * x228,
            x102 * x11 * x230,
            x237 * x99,
            x239 * x52 * x66,
            x239 * x46 * x78,
            x240 * x46 * x66,
            x239 * x44 * x69,
            x240 * x44 * x78,
            x241 * x44 * x66,
            x239 * x81 * x9,
            x240 * x69 * x9,
            x238 * x241 * x60,
            x238 * x242,
            x11 * x239 * x92,
            x11 * x240 * x81,
            x11 * x241 * x69,
            x231 * x242,
            x222
            * (
                x109 * x236
                + x3
                * (
                    x128 * (x233 + x234)
                    + x20 * (4 * x215 + 4 * x216 + x98)
                    + 5 * x218
                    + x219
                    + 4 * x229
                )
            ),
        ]
    )


def diag_quadrupole3d_30(a, A, b, B, C):
    """Cartesian 3D (fs) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x9 = -x1 - C[0]
    x10 = x7 * x9 ** 2
    x11 = 2 * x2
    x12 = x7 * x9
    x13 = 2 * x3
    x14 = x10 + x8
    x15 = x14 * x2
    x16 = x12 * x13 + x15
    x17 = x16 * x2 + x3 * (x10 + x11 * x12 + 3 * x8)
    x18 = x2 * x5
    x19 = x18 * x6
    x20 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x21 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x22 = numpy.pi * x0 * x21
    x23 = x20 * x22
    x24 = -x0 * (a * A[1] + b * B[1])
    x25 = -x24 - A[1]
    x26 = x17 * x23
    x27 = -x0 * (a * A[2] + b * B[2])
    x28 = -x27 - A[2]
    x29 = x20 * x6
    x30 = x29 * x3
    x31 = x25 ** 2 * x29 + x30
    x32 = x21 * x6
    x33 = x3 * x32
    x34 = x28 ** 2 * x32 + x33
    x35 = x25 * x29
    x36 = x13 * x35 + x25 * x31
    x37 = x28 * x32
    x38 = x13 * x37 + x28 * x34
    x39 = -x24 - C[1]
    x40 = x29 * x39 ** 2
    x41 = x30 + x40
    x42 = x2 ** 2 * x7 + x8
    x43 = x13 * x19 + x2 * x42
    x44 = x29 * x39
    x45 = x25 * x41
    x46 = x13 * x44 + x45
    x47 = 2 * x25
    x48 = x25 * x46 + x3 * (3 * x30 + x40 + x44 * x47)
    x49 = x18 * x22
    x50 = x22 * x5
    x51 = -x27 - C[2]
    x52 = x32 * x51 ** 2
    x53 = x33 + x52
    x54 = x32 * x51
    x55 = x28 * x53
    x56 = x13 * x54 + x55
    x57 = numpy.pi * x0 * x20
    x58 = x18 * x57
    x59 = 2 * x28
    x60 = x28 * x56 + x3 * (3 * x33 + x52 + x54 * x59)
    x61 = x5 * x57

    # 30 item(s)
    return numpy.array(
        [
            x23
            * (
                x17 * x2
                + x3
                * (x11 * (x19 * x9 + x8) + 4 * x12 * x3 + x13 * (x12 + x19) + 2 * x15)
            ),
            x25 * x26,
            x26 * x28,
            x16 * x31 * x32,
            x16 * x23 * x25 * x28,
            x16 * x29 * x34,
            x14 * x32 * x36,
            x14 * x31 * x37,
            x14 * x34 * x35,
            x14 * x29 * x38,
            x32 * x41 * x43,
            x32 * x42 * x46,
            x37 * x41 * x42,
            x48 * x49,
            x28 * x46 * x49,
            x19 * x34 * x41,
            x50
            * (
                x25 * x48
                + x3
                * (x13 * (x35 + x44) + 4 * x30 * x39 + 2 * x45 + x47 * (x30 + x35 * x39))
            ),
            x28 * x48 * x50,
            x34 * x46 * x7,
            x38 * x41 * x7,
            x29 * x43 * x53,
            x35 * x42 * x53,
            x29 * x42 * x56,
            x19 * x31 * x53,
            x25 * x56 * x58,
            x58 * x60,
            x36 * x53 * x7,
            x31 * x56 * x7,
            x25 * x60 * x61,
            x61
            * (
                x28 * x60
                + x3
                * (x13 * (x37 + x54) + 4 * x33 * x51 + 2 * x55 + x59 * (x33 + x37 * x51))
            ),
        ]
    )


def diag_quadrupole3d_31(a, A, b, B, C):
    """Cartesian 3D (fp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x17 = 2 * x3
    x18 = x4 * x9
    x19 = x17 * x18
    x20 = x12 + x13
    x21 = x20 * x5
    x22 = x19 + x21
    x23 = x2 * x22
    x24 = x16 + x23
    x25 = 4 * x13 * x4
    x26 = x3 * (x10 + x18)
    x27 = x2 * x20
    x28 = x2 * (x11 + x13)
    x29 = x2 * x24 + x3 * (x21 + x25 + 2 * x26 + x27 + 2 * x28)
    x30 = x10 * x2
    x31 = x18 * x2
    x32 = 2 * x2
    x33 = x19 + x27
    x34 = x2 * x33 + x3 * (x15 + x18 * x32)
    x35 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x36 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x37 = numpy.pi * x0 * x36
    x38 = x35 * x37
    x39 = -x0 * (a * A[1] + b * B[1])
    x40 = -x39 - B[1]
    x41 = x2 * x9
    x42 = x38 * (x2 * x34 + x3 * (x17 * (x18 + x41) + x25 + 2 * x27 + x32 * (x13 + x31)))
    x43 = -x0 * (a * A[2] + b * B[2])
    x44 = -x43 - B[2]
    x45 = -x39 - A[1]
    x46 = x29 * x38
    x47 = x3 * x8
    x48 = x35 * x47
    x49 = x35 * x8
    x50 = x45 * x49
    x51 = x40 * x50
    x52 = x48 + x51
    x53 = x36 * x8
    x54 = x34 * x38
    x55 = -x43 - A[2]
    x56 = x36 * x47
    x57 = x53 * x55
    x58 = x44 * x57
    x59 = x56 + x58
    x60 = x45 ** 2 * x49
    x61 = x48 + x60
    x62 = x40 * x49
    x63 = x3 * (x50 + x62) + x45 * x52
    x64 = x44 * x53
    x65 = x53 * x55 ** 2
    x66 = x56 + x65
    x67 = x3 * (x57 + x64) + x55 * x59
    x68 = x17 * x50 + x45 * x61
    x69 = 3 * x48
    x70 = 2 * x45
    x71 = x3 * (x60 + x62 * x70 + x69) + x45 * x63
    x72 = x17 * x57 + x55 * x66
    x73 = 3 * x56
    x74 = 2 * x55
    x75 = x3 * (x64 * x74 + x65 + x73) + x55 * x67
    x76 = -x39 - C[1]
    x77 = x49 * x76 ** 2
    x78 = x48 + x77
    x79 = x2 ** 2 * x9
    x80 = x13 + x30
    x81 = x2 * x80 + x3 * (x10 + x41)
    x82 = x2 * x81 + x3 * (x10 * x32 + x14 + x79)
    x83 = x49 * x76
    x84 = x17 * x83
    x85 = x40 * x78
    x86 = x84 + x85
    x87 = x13 + x79
    x88 = x17 * x41 + x2 * x87
    x89 = x45 * x78
    x90 = x84 + x89
    x91 = x62 * x76
    x92 = x69 + x77
    x93 = x3 * (2 * x91 + x92)
    x94 = x45 * x86
    x95 = x93 + x94
    x96 = x3 * (x70 * x83 + x92) + x45 * x90
    x97 = 4 * x48 * x76
    x98 = x3 * (x62 + x83)
    x99 = x45 * (x48 + x91)
    x100 = x3 * (x85 + x89 + x97 + 2 * x98 + 2 * x99) + x45 * x95
    x101 = x37 * x7
    x102 = x100 * x101
    x103 = x101 * x2
    x104 = x50 * x76
    x105 = x101 * (
        x3 * (x17 * (x50 + x83) + x70 * (x104 + x48) + 2 * x89 + x97) + x45 * x96
    )
    x106 = -x43 - C[2]
    x107 = x106 ** 2 * x53
    x108 = x107 + x56
    x109 = x106 * x53
    x110 = x109 * x17
    x111 = x108 * x44
    x112 = x110 + x111
    x113 = x108 * x55
    x114 = x110 + x113
    x115 = x106 * x64
    x116 = x107 + x73
    x117 = x3 * (2 * x115 + x116)
    x118 = x112 * x55
    x119 = x117 + x118
    x120 = numpy.pi * x0 * x35 * x7
    x121 = x120 * x2
    x122 = x114 * x55 + x3 * (x109 * x74 + x116)
    x123 = 4 * x106 * x56
    x124 = x3 * (x109 + x64)
    x125 = x55 * (x115 + x56)
    x126 = x119 * x55 + x3 * (x111 + x113 + x123 + 2 * x124 + 2 * x125)
    x127 = x120 * x126
    x128 = x106 * x57
    x129 = x120 * (
        x122 * x55 + x3 * (2 * x113 + x123 + x17 * (x109 + x57) + x74 * (x128 + x56))
    )

    # 90 item(s)
    return numpy.array(
        [
            x38
            * (
                x2 * x29
                + x3
                * (
                    2 * x16
                    + x17 * (x11 + x14 + x30 + x31)
                    + 2 * x23
                    + x32 * (x26 + x28)
                    + x34
                )
            ),
            x40 * x42,
            x42 * x44,
            x45 * x46,
            x34 * x52 * x53,
            x44 * x45 * x54,
            x46 * x55,
            x40 * x54 * x55,
            x34 * x49 * x59,
            x24 * x53 * x61,
            x33 * x53 * x63,
            x33 * x61 * x64,
            x24 * x38 * x45 * x55,
            x33 * x52 * x57,
            x33 * x50 * x59,
            x24 * x49 * x66,
            x33 * x62 * x66,
            x33 * x49 * x67,
            x22 * x53 * x68,
            x20 * x53 * x71,
            x20 * x64 * x68,
            x22 * x57 * x61,
            x20 * x57 * x63,
            x20 * x59 * x61,
            x22 * x50 * x66,
            x20 * x52 * x66,
            x20 * x50 * x67,
            x22 * x49 * x72,
            x20 * x62 * x72,
            x20 * x49 * x75,
            x53 * x78 * x82,
            x53 * x86 * x88,
            x64 * x78 * x88,
            x53 * x81 * x90,
            x53 * x87 * x95,
            x64 * x87 * x90,
            x57 * x78 * x81,
            x57 * x86 * x87,
            x59 * x78 * x87,
            x53 * x80 * x96,
            x102 * x2,
            x103 * x44 * x96,
            x57 * x80 * x90,
            x103 * x55 * x95,
            x41 * x59 * x90,
            x66 * x78 * x80,
            x41 * x66 * x86,
            x41 * x67 * x78,
            x105 * x5,
            x101
            * (
                x100 * x45
                + x3
                * (
                    x17 * (x104 + x51 + x69 + x91)
                    + x70 * (x98 + x99)
                    + 2 * x93
                    + 2 * x94
                    + x96
                )
            ),
            x105 * x44,
            x101 * x5 * x55 * x96,
            x102 * x55,
            x59 * x9 * x96,
            x10 * x66 * x90,
            x66 * x9 * x95,
            x67 * x9 * x90,
            x10 * x72 * x78,
            x72 * x86 * x9,
            x75 * x78 * x9,
            x108 * x49 * x82,
            x108 * x62 * x88,
            x112 * x49 * x88,
            x108 * x50 * x81,
            x108 * x52 * x87,
            x112 * x50 * x87,
            x114 * x49 * x81,
            x114 * x62 * x87,
            x119 * x49 * x87,
            x108 * x61 * x80,
            x108 * x41 * x63,
            x112 * x41 * x61,
            x114 * x50 * x80,
            x114 * x41 * x52,
            x119 * x121 * x45,
            x122 * x49 * x80,
            x121 * x122 * x40,
            x127 * x2,
            x10 * x108 * x68,
            x108 * x71 * x9,
            x112 * x68 * x9,
            x10 * x114 * x61,
            x114 * x63 * x9,
            x119 * x61 * x9,
            x120 * x122 * x45 * x5,
            x122 * x52 * x9,
            x127 * x45,
            x129 * x5,
            x129 * x40,
            x120
            * (
                x126 * x55
                + x3
                * (
                    2 * x117
                    + 2 * x118
                    + x122
                    + x17 * (x115 + x128 + x58 + x73)
                    + x74 * (x124 + x125)
                )
            ),
        ]
    )


def diag_quadrupole3d_32(a, A, b, B, C):
    """Cartesian 3D (fd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x18 = 2 * x3
    x19 = x4 * x8
    x20 = x18 * x19
    x21 = x10 + x9
    x22 = x12 * x21
    x23 = x20 + x22
    x24 = x12 * x23
    x25 = x17 + x24
    x26 = x2 * x25
    x27 = x10 + x15
    x28 = x12 * x27
    x29 = x3 * (x14 + x19)
    x30 = 4 * x3
    x31 = x19 * x30
    x32 = 2 * x29 + x31
    x33 = x3 * (2 * x22 + 2 * x28 + x32)
    x34 = x26 + x33
    x35 = x12 ** 2 * x8
    x36 = x3 * (x16 + x35)
    x37 = x2 * x23
    x38 = 2 * x37
    x39 = x2 * (x28 + x29)
    x40 = x2 * x34 + x3 * (3 * x17 + x24 + 2 * x36 + x38 + 2 * x39)
    x41 = x17 + x37
    x42 = x2 * x41
    x43 = 2 * x2
    x44 = x2 * x21
    x45 = x2 * x27
    x46 = 2 * x45
    x47 = x3 * (x22 + x32 + x44 + x46)
    x48 = x10 + x35
    x49 = x2 * x48
    x50 = x14 * x18 + x49
    x51 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x52 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x53 = numpy.pi * x0 * x52
    x54 = x51 * x53
    x55 = -x0 * (a * A[1] + b * B[1])
    x56 = -x55 - B[1]
    x57 = x42 + x47
    x58 = x14 * x2
    x59 = x19 * x2
    x60 = x20 + x44
    x61 = x2 * x60 + x3 * (x11 + x19 * x43 + x9)
    x62 = x54 * (
        x2 * x57
        + x3 * (2 * x17 + x18 * (x11 + x15 + x58 + x59) + x38 + x43 * (x29 + x45) + x61)
    )
    x63 = -x0 * (a * A[2] + b * B[2])
    x64 = -x63 - B[2]
    x65 = x51 * x7
    x66 = x3 * x65
    x67 = x56 ** 2 * x65
    x68 = x66 + x67
    x69 = x2 * x8
    x70 = x2 * x61 + x3 * (x18 * (x19 + x69) + x31 + x43 * (x10 + x59) + 2 * x44)
    x71 = x52 * x7
    x72 = x54 * x64
    x73 = x3 * x71
    x74 = x64 ** 2 * x71
    x75 = x73 + x74
    x76 = -x55 - A[1]
    x77 = x40 * x54
    x78 = x65 * x76
    x79 = x56 * x78
    x80 = x66 + x79
    x81 = x56 * x65
    x82 = x68 * x76
    x83 = x18 * x81 + x82
    x84 = x64 * x71
    x85 = -x63 - A[2]
    x86 = x54 * x85
    x87 = x71 * x85
    x88 = x64 * x87
    x89 = x73 + x88
    x90 = x75 * x85
    x91 = x18 * x84 + x90
    x92 = x65 * x76 ** 2
    x93 = x66 + x92
    x94 = x3 * (x78 + x81)
    x95 = x76 * x80
    x96 = x94 + x95
    x97 = 3 * x66
    x98 = 2 * x76
    x99 = x81 * x98 + x97
    x100 = x3 * (x67 + x99) + x76 * x83
    x101 = x71 * x85 ** 2
    x102 = x101 + x73
    x103 = x3 * (x84 + x87)
    x104 = x85 * x89
    x105 = x103 + x104
    x106 = 3 * x73
    x107 = 2 * x85
    x108 = x106 + x107 * x84
    x109 = x3 * (x108 + x74) + x85 * x91
    x110 = x18 * x78 + x76 * x93
    x111 = x3 * (x92 + x99) + x76 * x96
    x112 = 4 * x66
    x113 = x100 * x76 + x3 * (x112 * x56 + 2 * x82 + 2 * x94 + 2 * x95)
    x114 = x102 * x85 + x18 * x87
    x115 = x105 * x85 + x3 * (x101 + x108)
    x116 = 4 * x73
    x117 = x109 * x85 + x3 * (2 * x103 + 2 * x104 + x116 * x64 + 2 * x90)
    x118 = -x55 - C[1]
    x119 = x118 ** 2 * x65
    x120 = x119 + x66
    x121 = x11 + x14 * x43
    x122 = x2 * x50 + x3 * (x121 + x35)
    x123 = x3 * (x14 + x69)
    x124 = x10 + x58
    x125 = x124 * x2
    x126 = x122 * x2 + x3 * (2 * x123 + 2 * x125 + x14 * x30 + 2 * x49)
    x127 = x118 * x65
    x128 = x127 * x18
    x129 = x120 * x56
    x130 = x128 + x129
    x131 = x2 ** 2 * x8
    x132 = x123 + x125
    x133 = x132 * x2 + x3 * (x121 + x131)
    x134 = x118 * x81
    x135 = 2 * x134
    x136 = x119 + x97
    x137 = x3 * (x135 + x136)
    x138 = x130 * x56
    x139 = x137 + x138
    x140 = x10 + x131
    x141 = x140 * x2 + x18 * x69
    x142 = x120 * x76
    x143 = x128 + x142
    x144 = x130 * x76
    x145 = x137 + x144
    x146 = x139 * x76
    x147 = x134 + x66
    x148 = x147 * x56
    x149 = x3 * (x127 + x81)
    x150 = x112 * x118
    x151 = 2 * x149 + x150
    x152 = x3 * (2 * x129 + 2 * x148 + x151)
    x153 = x146 + x152
    x154 = x143 * x76 + x3 * (x127 * x98 + x136)
    x155 = x145 * x76
    x156 = x147 * x76
    x157 = 2 * x156
    x158 = x3 * (x129 + x142 + x151 + x157)
    x159 = x155 + x158
    x160 = x3 * (x135 + x67 + x97)
    x161 = 2 * x144
    x162 = x76 * (x148 + x149)
    x163 = x153 * x76 + x3 * (3 * x137 + x138 + 2 * x160 + x161 + 2 * x162)
    x164 = x53 * x6
    x165 = x163 * x164
    x166 = x164 * x64
    x167 = x118 * x78
    x168 = x154 * x76 + x3 * (2 * x142 + x150 + x18 * (x127 + x78) + x98 * (x167 + x66))
    x169 = x159 * x76 + x3 * (
        2 * x137 + x154 + x161 + x18 * (x134 + x167 + x79 + x97) + x98 * (x149 + x156)
    )
    x170 = x13 * x53
    x171 = -x63 - C[2]
    x172 = x171 ** 2 * x71
    x173 = x172 + x73
    x174 = x171 * x71
    x175 = x174 * x18
    x176 = x173 * x64
    x177 = x175 + x176
    x178 = x171 * x84
    x179 = 2 * x178
    x180 = x106 + x172
    x181 = x3 * (x179 + x180)
    x182 = x177 * x64
    x183 = x181 + x182
    x184 = x173 * x85
    x185 = x175 + x184
    x186 = x177 * x85
    x187 = x181 + x186
    x188 = x183 * x85
    x189 = x178 + x73
    x190 = x189 * x64
    x191 = x3 * (x174 + x84)
    x192 = x116 * x171
    x193 = 2 * x191 + x192
    x194 = x3 * (2 * x176 + 2 * x190 + x193)
    x195 = x188 + x194
    x196 = numpy.pi * x0 * x51
    x197 = x196 * x6
    x198 = x185 * x85 + x3 * (x107 * x174 + x180)
    x199 = x187 * x85
    x200 = x189 * x85
    x201 = 2 * x200
    x202 = x3 * (x176 + x184 + x193 + x201)
    x203 = x199 + x202
    x204 = x197 * x56
    x205 = x3 * (x106 + x179 + x74)
    x206 = 2 * x186
    x207 = x85 * (x190 + x191)
    x208 = x195 * x85 + x3 * (3 * x181 + x182 + 2 * x205 + x206 + 2 * x207)
    x209 = x197 * x208
    x210 = x13 * x196
    x211 = x171 * x87
    x212 = x198 * x85 + x3 * (x107 * (x211 + x73) + x18 * (x174 + x87) + 2 * x184 + x192)
    x213 = x203 * x85 + x3 * (
        x107 * (x191 + x200) + x18 * (x106 + x178 + x211 + x88) + 2 * x181 + x198 + x206
    )

    # 180 item(s)
    return numpy.array(
        [
            x54
            * (
                x2 * x40
                + x3
                * (
                    x18 * (x28 + 3 * x29 + x46 + x50)
                    + 2 * x26
                    + 2 * x33
                    + 2 * x42
                    + x43 * (x36 + x39)
                    + 2 * x47
                )
            ),
            x56 * x62,
            x62 * x64,
            x68 * x70 * x71,
            x56 * x70 * x72,
            x65 * x70 * x75,
            x76 * x77,
            x57 * x71 * x80,
            x57 * x72 * x76,
            x61 * x71 * x83,
            x61 * x80 * x84,
            x61 * x75 * x78,
            x77 * x85,
            x56 * x57 * x86,
            x57 * x65 * x89,
            x61 * x68 * x87,
            x61 * x81 * x89,
            x61 * x65 * x91,
            x34 * x71 * x93,
            x41 * x71 * x96,
            x41 * x84 * x93,
            x100 * x60 * x71,
            x60 * x84 * x96,
            x60 * x75 * x93,
            x34 * x76 * x86,
            x41 * x80 * x87,
            x41 * x78 * x89,
            x60 * x83 * x87,
            x60 * x80 * x89,
            x60 * x78 * x91,
            x102 * x34 * x65,
            x102 * x41 * x81,
            x105 * x41 * x65,
            x102 * x60 * x68,
            x105 * x60 * x81,
            x109 * x60 * x65,
            x110 * x25 * x71,
            x111 * x23 * x71,
            x110 * x23 * x84,
            x113 * x21 * x71,
            x111 * x21 * x84,
            x110 * x21 * x75,
            x25 * x87 * x93,
            x23 * x87 * x96,
            x23 * x89 * x93,
            x100 * x21 * x87,
            x21 * x89 * x96,
            x21 * x91 * x93,
            x102 * x25 * x78,
            x102 * x23 * x80,
            x105 * x23 * x78,
            x102 * x21 * x83,
            x105 * x21 * x80,
            x109 * x21 * x78,
            x114 * x25 * x65,
            x114 * x23 * x81,
            x115 * x23 * x65,
            x114 * x21 * x68,
            x115 * x21 * x81,
            x117 * x21 * x65,
            x120 * x126 * x71,
            x130 * x133 * x71,
            x120 * x133 * x84,
            x139 * x141 * x71,
            x130 * x141 * x84,
            x120 * x141 * x75,
            x122 * x143 * x71,
            x132 * x145 * x71,
            x132 * x143 * x84,
            x140 * x153 * x71,
            x140 * x145 * x84,
            x140 * x143 * x75,
            x120 * x122 * x87,
            x130 * x132 * x87,
            x120 * x132 * x89,
            x139 * x140 * x87,
            x130 * x140 * x89,
            x120 * x140 * x91,
            x154 * x50 * x71,
            x124 * x159 * x71,
            x124 * x154 * x84,
            x165 * x2,
            x159 * x166 * x2,
            x154 * x69 * x75,
            x143 * x50 * x87,
            x124 * x145 * x87,
            x124 * x143 * x89,
            x153 * x164 * x2 * x85,
            x145 * x69 * x89,
            x143 * x69 * x91,
            x102 * x120 * x50,
            x102 * x124 * x130,
            x105 * x120 * x124,
            x102 * x139 * x69,
            x105 * x130 * x69,
            x109 * x120 * x69,
            x168 * x48 * x71,
            x169 * x170,
            x168 * x170 * x64,
            x164
            * (
                x163 * x76
                + x3
                * (
                    2 * x146
                    + 2 * x152
                    + 2 * x155
                    + 2 * x158
                    + x18 * (x148 + 3 * x149 + x157 + x83)
                    + x98 * (x160 + x162)
                )
            ),
            x166 * x169,
            x168 * x75 * x8,
            x154 * x48 * x87,
            x159 * x170 * x85,
            x14 * x154 * x89,
            x165 * x85,
            x159 * x8 * x89,
            x154 * x8 * x91,
            x102 * x143 * x48,
            x102 * x14 * x145,
            x105 * x14 * x143,
            x102 * x153 * x8,
            x105 * x145 * x8,
            x109 * x143 * x8,
            x114 * x120 * x48,
            x114 * x130 * x14,
            x115 * x120 * x14,
            x114 * x139 * x8,
            x115 * x130 * x8,
            x117 * x120 * x8,
            x126 * x173 * x65,
            x133 * x173 * x81,
            x133 * x177 * x65,
            x141 * x173 * x68,
            x141 * x177 * x81,
            x141 * x183 * x65,
            x122 * x173 * x78,
            x132 * x173 * x80,
            x132 * x177 * x78,
            x140 * x173 * x83,
            x140 * x177 * x80,
            x140 * x183 * x78,
            x122 * x185 * x65,
            x132 * x185 * x81,
            x132 * x187 * x65,
            x140 * x185 * x68,
            x140 * x187 * x81,
            x140 * x195 * x65,
            x173 * x50 * x93,
            x124 * x173 * x96,
            x124 * x177 * x93,
            x100 * x173 * x69,
            x177 * x69 * x96,
            x183 * x69 * x93,
            x185 * x50 * x78,
            x124 * x185 * x80,
            x124 * x187 * x78,
            x185 * x69 * x83,
            x187 * x69 * x80,
            x195 * x197 * x2 * x76,
            x198 * x50 * x65,
            x124 * x198 * x81,
            x124 * x203 * x65,
            x198 * x68 * x69,
            x2 * x203 * x204,
            x2 * x209,
            x110 * x173 * x48,
            x111 * x14 * x173,
            x110 * x14 * x177,
            x113 * x173 * x8,
            x111 * x177 * x8,
            x110 * x183 * x8,
            x185 * x48 * x93,
            x14 * x185 * x96,
            x14 * x187 * x93,
            x100 * x185 * x8,
            x187 * x8 * x96,
            x195 * x8 * x93,
            x198 * x48 * x78,
            x14 * x198 * x80,
            x203 * x210 * x76,
            x198 * x8 * x83,
            x203 * x8 * x80,
            x209 * x76,
            x212 * x48 * x65,
            x210 * x212 * x56,
            x210 * x213,
            x212 * x68 * x8,
            x204 * x213,
            x197
            * (
                x208 * x85
                + x3
                * (
                    x107 * (x205 + x207)
                    + x18 * (x190 + 3 * x191 + x201 + x91)
                    + 2 * x188
                    + 2 * x194
                    + 2 * x199
                    + 2 * x202
                )
            ),
        ]
    )


def diag_quadrupole3d_33(a, A, b, B, C):
    """Cartesian 3D (ff) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x30 = x3 * (x21 + x28)
    x31 = 2 * x29 + 3 * x30
    x32 = x3 * (2 * x18 + 3 * x25 + x31)
    x33 = x25 + x30
    x34 = x33 * x4
    x35 = 4 * x3
    x36 = x11 * x35
    x37 = 2 * x12 + x36
    x38 = x3 * (2 * x16 + 2 * x23 + x37)
    x39 = x34 + x38
    x40 = x2 * x39
    x41 = x32 + x40
    x42 = x2 * x33
    x43 = 3 * x12
    x44 = x19 * x9
    x45 = x13 + x26
    x46 = x4 * x45
    x47 = x44 + x46
    x48 = x3 * (3 * x16 + x43 + x47)
    x49 = x2 * (x18 + x29)
    x50 = x2 * x41 + x3 * (x34 + 4 * x38 + 3 * x42 + 2 * x48 + 2 * x49)
    x51 = x38 + x42
    x52 = x2 * x51
    x53 = 2 * x2
    x54 = x17 * x2
    x55 = x3 * (3 * x26 + x27)
    x56 = x2 * x47
    x57 = x55 + x56
    x58 = x2 * x24
    x59 = 2 * x58
    x60 = x3 * (x25 + x31 + 2 * x54 + x59)
    x61 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x62 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x63 = numpy.pi * x0 * x62
    x64 = x61 * x63
    x65 = -x0 * (a * A[1] + b * B[1])
    x66 = -x65 - B[1]
    x67 = x52 + x60
    x68 = x30 + x58
    x69 = x2 * x68
    x70 = x2 * x22
    x71 = x15 * x2
    x72 = 2 * x71
    x73 = x3 * (x23 + x37 + x70 + x72)
    x74 = x2 * x45
    x75 = x44 + x74
    x76 = x64 * (
        x2 * x67
        + x3
        * (
            x19 * (x16 + x43 + x72 + x75)
            + 2 * x38
            + 2 * x42
            + x53 * (x29 + x54)
            + 2 * x69
            + 2 * x73
        )
    )
    x77 = -x0 * (a * A[2] + b * B[2])
    x78 = -x77 - B[2]
    x79 = x61 * x7
    x80 = x3 * x79
    x81 = x66 ** 2 * x79
    x82 = x80 + x81
    x83 = x69 + x73
    x84 = x2 * x9
    x85 = x11 * x2
    x86 = x20 + x70
    x87 = x2 * x86 + x3 * (x11 * x53 + x21 + x27)
    x88 = x2 * x83 + x3 * (
        x19 * (x14 + x27 + x84 + x85) + 2 * x30 + x53 * (x12 + x71) + x59 + x87
    )
    x89 = x62 * x7
    x90 = x64 * x78
    x91 = x3 * x89
    x92 = x78 ** 2 * x89
    x93 = x91 + x92
    x94 = x66 * x79
    x95 = x19 * x94
    x96 = x66 * x82
    x97 = x95 + x96
    x98 = x2 * x8
    x99 = x2 * x87 + x3 * (x19 * (x11 + x98) + x36 + x53 * (x13 + x85) + 2 * x70)
    x100 = x78 * x89
    x101 = x100 * x19
    x102 = x78 * x93
    x103 = x101 + x102
    x104 = -x65 - A[1]
    x105 = x50 * x64
    x106 = x104 * x79
    x107 = x106 * x66
    x108 = x107 + x80
    x109 = x104 * x82
    x110 = x109 + x95
    x111 = 3 * x80
    x112 = x3 * (x111 + 3 * x81)
    x113 = x104 * x97
    x114 = x112 + x113
    x115 = -x77 - A[2]
    x116 = x115 * x64
    x117 = x115 * x89
    x118 = x117 * x78
    x119 = x118 + x91
    x120 = x115 * x93
    x121 = x101 + x120
    x122 = 3 * x91
    x123 = x3 * (x122 + 3 * x92)
    x124 = x103 * x115
    x125 = x123 + x124
    x126 = x104 ** 2 * x79
    x127 = x126 + x80
    x128 = x3 * (x106 + x94)
    x129 = x104 * x108
    x130 = x128 + x129
    x131 = 2 * x104
    x132 = x111 + x131 * x94
    x133 = x3 * (x132 + x81)
    x134 = x104 * x110
    x135 = x133 + x134
    x136 = x66 * x80
    x137 = x104 * x114 + x3 * (3 * x109 + 8 * x136 + x96)
    x138 = x115 ** 2 * x89
    x139 = x138 + x91
    x140 = x3 * (x100 + x117)
    x141 = x115 * x119
    x142 = x140 + x141
    x143 = 2 * x115
    x144 = x100 * x143 + x122
    x145 = x3 * (x144 + x92)
    x146 = x115 * x121
    x147 = x145 + x146
    x148 = x78 * x91
    x149 = x115 * x125 + x3 * (x102 + 3 * x120 + 8 * x148)
    x150 = x104 * x127 + x106 * x19
    x151 = x104 * x130 + x3 * (x126 + x132)
    x152 = x104 * x135 + x3 * (2 * x109 + 2 * x128 + 2 * x129 + 4 * x136)
    x153 = x104 * x137 + x3 * (2 * x112 + 2 * x113 + 3 * x133 + 3 * x134)
    x154 = x115 * x139 + x117 * x19
    x155 = x115 * x142 + x3 * (x138 + x144)
    x156 = x115 * x147 + x3 * (2 * x120 + 2 * x140 + 2 * x141 + 4 * x148)
    x157 = x115 * x149 + x3 * (2 * x123 + 2 * x124 + 3 * x145 + 3 * x146)
    x158 = -x65 - C[1]
    x159 = x158 ** 2 * x79
    x160 = x159 + x80
    x161 = x2 * x57 + x3 * (8 * x3 * x9 + x46 + 3 * x74)
    x162 = x27 + x53 * x9
    x163 = x3 * (x162 + x26)
    x164 = x2 * x75
    x165 = x161 * x2 + x3 * (3 * x163 + 3 * x164 + 2 * x55 + 2 * x56)
    x166 = x158 * x79
    x167 = x166 * x19
    x168 = x160 * x66
    x169 = x167 + x168
    x170 = x163 + x164
    x171 = x3 * (x9 + x98)
    x172 = x13 + x84
    x173 = x172 * x2
    x174 = x170 * x2 + x3 * (2 * x171 + 2 * x173 + x35 * x9 + 2 * x74)
    x175 = x158 * x94
    x176 = x111 + 2 * x175
    x177 = x3 * (x159 + x176)
    x178 = x169 * x66
    x179 = x177 + x178
    x180 = x2 ** 2 * x8
    x181 = x171 + x173
    x182 = x181 * x2 + x3 * (x162 + x180)
    x183 = x179 * x66
    x184 = x175 + x80
    x185 = x184 * x66
    x186 = x3 * (x166 + x94)
    x187 = 4 * x158 * x80
    x188 = 2 * x186 + x187
    x189 = x3 * (2 * x168 + 2 * x185 + x188)
    x190 = x183 + x189
    x191 = x13 + x180
    x192 = x19 * x98 + x191 * x2
    x193 = x104 * x160
    x194 = x167 + x193
    x195 = x104 * x169
    x196 = x177 + x195
    x197 = x104 * x179
    x198 = x189 + x197
    x199 = x185 + x186
    x200 = x199 * x66
    x201 = x3 * (x176 + x81)
    x202 = 3 * x177 + 2 * x201
    x203 = x3 * (3 * x178 + 2 * x200 + x202)
    x204 = x104 * x190
    x205 = x203 + x204
    x206 = x104 * x194 + x3 * (x111 + x131 * x166 + x159)
    x207 = x104 * x196
    x208 = x104 * x184
    x209 = 2 * x208
    x210 = x3 * (x168 + x188 + x193 + x209)
    x211 = x207 + x210
    x212 = x104 * x198
    x213 = x104 * x199
    x214 = 2 * x195
    x215 = x3 * (x178 + x202 + 2 * x213 + x214)
    x216 = x212 + x215
    x217 = 3 * x186
    x218 = x3 * (3 * x185 + x217 + x97)
    x219 = x104 * (x200 + x201)
    x220 = x104 * x205 + x3 * (x183 + 4 * x189 + 3 * x197 + 2 * x218 + 2 * x219)
    x221 = x6 * x63
    x222 = x220 * x221
    x223 = x2 * x221
    x224 = x106 * x158
    x225 = x104 * x206 + x3 * (
        x131 * (x224 + x80) + x187 + x19 * (x106 + x166) + 2 * x193
    )
    x226 = x104 * x211 + x3 * (
        x131 * (x186 + x208) + 2 * x177 + x19 * (x107 + x111 + x175 + x224) + x206 + x214
    )
    x227 = x221 * (
        x104 * x216
        + x3
        * (
            x131 * (x201 + x213)
            + 2 * x189
            + x19 * (x110 + x185 + x209 + x217)
            + 2 * x197
            + 2 * x207
            + 2 * x210
        )
    )
    x228 = x221 * x4
    x229 = -x77 - C[2]
    x230 = x229 ** 2 * x89
    x231 = x230 + x91
    x232 = x229 * x89
    x233 = x19 * x232
    x234 = x231 * x78
    x235 = x233 + x234
    x236 = x100 * x229
    x237 = x122 + 2 * x236
    x238 = x3 * (x230 + x237)
    x239 = x235 * x78
    x240 = x238 + x239
    x241 = x240 * x78
    x242 = x236 + x91
    x243 = x242 * x78
    x244 = x3 * (x100 + x232)
    x245 = 4 * x229 * x91
    x246 = 2 * x244 + x245
    x247 = x3 * (2 * x234 + 2 * x243 + x246)
    x248 = x241 + x247
    x249 = x115 * x231
    x250 = x233 + x249
    x251 = x115 * x235
    x252 = x238 + x251
    x253 = x115 * x240
    x254 = x247 + x253
    x255 = x243 + x244
    x256 = x255 * x78
    x257 = x3 * (x237 + x92)
    x258 = 3 * x238 + 2 * x257
    x259 = x3 * (3 * x239 + 2 * x256 + x258)
    x260 = x115 * x248
    x261 = x259 + x260
    x262 = numpy.pi * x0 * x6 * x61
    x263 = x2 * x262
    x264 = x115 * x250 + x3 * (x122 + x143 * x232 + x230)
    x265 = x115 * x252
    x266 = x115 * x242
    x267 = 2 * x266
    x268 = x3 * (x234 + x246 + x249 + x267)
    x269 = x265 + x268
    x270 = x115 * x254
    x271 = x115 * x255
    x272 = 2 * x251
    x273 = x3 * (x239 + x258 + 2 * x271 + x272)
    x274 = x270 + x273
    x275 = 3 * x244
    x276 = x3 * (x103 + 3 * x243 + x275)
    x277 = x115 * (x256 + x257)
    x278 = x115 * x261 + x3 * (x241 + 4 * x247 + 3 * x253 + 2 * x276 + 2 * x277)
    x279 = x262 * x278
    x280 = x262 * x4
    x281 = x117 * x229
    x282 = x115 * x264 + x3 * (
        x143 * (x281 + x91) + x19 * (x117 + x232) + x245 + 2 * x249
    )
    x283 = x115 * x269 + x3 * (
        x143 * (x244 + x266) + x19 * (x118 + x122 + x236 + x281) + 2 * x238 + x264 + x272
    )
    x284 = x262 * (
        x115 * x274
        + x3
        * (
            x143 * (x257 + x271)
            + x19 * (x121 + x243 + x267 + x275)
            + 2 * x247
            + 2 * x253
            + 2 * x265
            + 2 * x268
        )
    )

    # 300 item(s)
    return numpy.array(
        [
            x64
            * (
                x2 * x50
                + x3
                * (
                    x19 * (x18 + 4 * x29 + 3 * x54 + x57)
                    + 2 * x32
                    + 2 * x40
                    + 3 * x52
                    + x53 * (x48 + x49)
                    + 3 * x60
                )
            ),
            x66 * x76,
            x76 * x78,
            x82 * x88 * x89,
            x66 * x88 * x90,
            x79 * x88 * x93,
            x89 * x97 * x99,
            x100 * x82 * x99,
            x93 * x94 * x99,
            x103 * x79 * x99,
            x104 * x105,
            x108 * x67 * x89,
            x104 * x67 * x90,
            x110 * x83 * x89,
            x100 * x108 * x83,
            x106 * x83 * x93,
            x114 * x87 * x89,
            x100 * x110 * x87,
            x108 * x87 * x93,
            x103 * x106 * x87,
            x105 * x115,
            x116 * x66 * x67,
            x119 * x67 * x79,
            x117 * x82 * x83,
            x119 * x83 * x94,
            x121 * x79 * x83,
            x117 * x87 * x97,
            x119 * x82 * x87,
            x121 * x87 * x94,
            x125 * x79 * x87,
            x127 * x41 * x89,
            x130 * x51 * x89,
            x100 * x127 * x51,
            x135 * x68 * x89,
            x100 * x130 * x68,
            x127 * x68 * x93,
            x137 * x86 * x89,
            x100 * x135 * x86,
            x130 * x86 * x93,
            x103 * x127 * x86,
            x104 * x116 * x41,
            x108 * x117 * x51,
            x106 * x119 * x51,
            x110 * x117 * x68,
            x108 * x119 * x68,
            x106 * x121 * x68,
            x114 * x117 * x86,
            x110 * x119 * x86,
            x108 * x121 * x86,
            x106 * x125 * x86,
            x139 * x41 * x79,
            x139 * x51 * x94,
            x142 * x51 * x79,
            x139 * x68 * x82,
            x142 * x68 * x94,
            x147 * x68 * x79,
            x139 * x86 * x97,
            x142 * x82 * x86,
            x147 * x86 * x94,
            x149 * x79 * x86,
            x150 * x39 * x89,
            x151 * x33 * x89,
            x100 * x150 * x33,
            x152 * x24 * x89,
            x100 * x151 * x24,
            x150 * x24 * x93,
            x153 * x22 * x89,
            x100 * x152 * x22,
            x151 * x22 * x93,
            x103 * x150 * x22,
            x117 * x127 * x39,
            x117 * x130 * x33,
            x119 * x127 * x33,
            x117 * x135 * x24,
            x119 * x130 * x24,
            x121 * x127 * x24,
            x117 * x137 * x22,
            x119 * x135 * x22,
            x121 * x130 * x22,
            x125 * x127 * x22,
            x106 * x139 * x39,
            x108 * x139 * x33,
            x106 * x142 * x33,
            x110 * x139 * x24,
            x108 * x142 * x24,
            x106 * x147 * x24,
            x114 * x139 * x22,
            x110 * x142 * x22,
            x108 * x147 * x22,
            x106 * x149 * x22,
            x154 * x39 * x79,
            x154 * x33 * x94,
            x155 * x33 * x79,
            x154 * x24 * x82,
            x155 * x24 * x94,
            x156 * x24 * x79,
            x154 * x22 * x97,
            x155 * x22 * x82,
            x156 * x22 * x94,
            x157 * x22 * x79,
            x160 * x165 * x89,
            x169 * x174 * x89,
            x100 * x160 * x174,
            x179 * x182 * x89,
            x100 * x169 * x182,
            x160 * x182 * x93,
            x190 * x192 * x89,
            x100 * x179 * x192,
            x169 * x192 * x93,
            x103 * x160 * x192,
            x161 * x194 * x89,
            x170 * x196 * x89,
            x100 * x170 * x194,
            x181 * x198 * x89,
            x100 * x181 * x196,
            x181 * x194 * x93,
            x191 * x205 * x89,
            x100 * x191 * x198,
            x191 * x196 * x93,
            x103 * x191 * x194,
            x117 * x160 * x161,
            x117 * x169 * x170,
            x119 * x160 * x170,
            x117 * x179 * x181,
            x119 * x169 * x181,
            x121 * x160 * x181,
            x117 * x190 * x191,
            x119 * x179 * x191,
            x121 * x169 * x191,
            x125 * x160 * x191,
            x206 * x57 * x89,
            x211 * x75 * x89,
            x100 * x206 * x75,
            x172 * x216 * x89,
            x100 * x172 * x211,
            x172 * x206 * x93,
            x2 * x222,
            x216 * x223 * x78,
            x211 * x93 * x98,
            x103 * x206 * x98,
            x117 * x194 * x57,
            x117 * x196 * x75,
            x119 * x194 * x75,
            x117 * x172 * x198,
            x119 * x172 * x196,
            x121 * x172 * x194,
            x115 * x205 * x223,
            x119 * x198 * x98,
            x121 * x196 * x98,
            x125 * x194 * x98,
            x139 * x160 * x57,
            x139 * x169 * x75,
            x142 * x160 * x75,
            x139 * x172 * x179,
            x142 * x169 * x172,
            x147 * x160 * x172,
            x139 * x190 * x98,
            x142 * x179 * x98,
            x147 * x169 * x98,
            x149 * x160 * x98,
            x225 * x47 * x89,
            x226 * x45 * x89,
            x100 * x225 * x45,
            x227 * x4,
            x226 * x228 * x78,
            x225 * x9 * x93,
            x221
            * (
                x104 * x220
                + x3
                * (
                    x131 * (x218 + x219)
                    + x19 * (x114 + x200 + 4 * x201 + 3 * x213)
                    + 2 * x203
                    + 2 * x204
                    + 3 * x212
                    + 3 * x215
                )
            ),
            x227 * x78,
            x226 * x8 * x93,
            x103 * x225 * x8,
            x117 * x206 * x47,
            x117 * x211 * x45,
            x119 * x206 * x45,
            x115 * x216 * x228,
            x119 * x211 * x9,
            x121 * x206 * x9,
            x115 * x222,
            x119 * x216 * x8,
            x121 * x211 * x8,
            x125 * x206 * x8,
            x139 * x194 * x47,
            x139 * x196 * x45,
            x142 * x194 * x45,
            x139 * x198 * x9,
            x142 * x196 * x9,
            x147 * x194 * x9,
            x139 * x205 * x8,
            x142 * x198 * x8,
            x147 * x196 * x8,
            x149 * x194 * x8,
            x154 * x160 * x47,
            x154 * x169 * x45,
            x155 * x160 * x45,
            x154 * x179 * x9,
            x155 * x169 * x9,
            x156 * x160 * x9,
            x154 * x190 * x8,
            x155 * x179 * x8,
            x156 * x169 * x8,
            x157 * x160 * x8,
            x165 * x231 * x79,
            x174 * x231 * x94,
            x174 * x235 * x79,
            x182 * x231 * x82,
            x182 * x235 * x94,
            x182 * x240 * x79,
            x192 * x231 * x97,
            x192 * x235 * x82,
            x192 * x240 * x94,
            x192 * x248 * x79,
            x106 * x161 * x231,
            x108 * x170 * x231,
            x106 * x170 * x235,
            x110 * x181 * x231,
            x108 * x181 * x235,
            x106 * x181 * x240,
            x114 * x191 * x231,
            x110 * x191 * x235,
            x108 * x191 * x240,
            x106 * x191 * x248,
            x161 * x250 * x79,
            x170 * x250 * x94,
            x170 * x252 * x79,
            x181 * x250 * x82,
            x181 * x252 * x94,
            x181 * x254 * x79,
            x191 * x250 * x97,
            x191 * x252 * x82,
            x191 * x254 * x94,
            x191 * x261 * x79,
            x127 * x231 * x57,
            x130 * x231 * x75,
            x127 * x235 * x75,
            x135 * x172 * x231,
            x130 * x172 * x235,
            x127 * x172 * x240,
            x137 * x231 * x98,
            x135 * x235 * x98,
            x130 * x240 * x98,
            x127 * x248 * x98,
            x106 * x250 * x57,
            x108 * x250 * x75,
            x106 * x252 * x75,
            x110 * x172 * x250,
            x108 * x172 * x252,
            x106 * x172 * x254,
            x114 * x250 * x98,
            x110 * x252 * x98,
            x108 * x254 * x98,
            x104 * x261 * x263,
            x264 * x57 * x79,
            x264 * x75 * x94,
            x269 * x75 * x79,
            x172 * x264 * x82,
            x172 * x269 * x94,
            x172 * x274 * x79,
            x264 * x97 * x98,
            x269 * x82 * x98,
            x263 * x274 * x66,
            x2 * x279,
            x150 * x231 * x47,
            x151 * x231 * x45,
            x150 * x235 * x45,
            x152 * x231 * x9,
            x151 * x235 * x9,
            x150 * x240 * x9,
            x153 * x231 * x8,
            x152 * x235 * x8,
            x151 * x240 * x8,
            x150 * x248 * x8,
            x127 * x250 * x47,
            x130 * x250 * x45,
            x127 * x252 * x45,
            x135 * x250 * x9,
            x130 * x252 * x9,
            x127 * x254 * x9,
            x137 * x250 * x8,
            x135 * x252 * x8,
            x130 * x254 * x8,
            x127 * x261 * x8,
            x106 * x264 * x47,
            x108 * x264 * x45,
            x106 * x269 * x45,
            x110 * x264 * x9,
            x108 * x269 * x9,
            x104 * x274 * x280,
            x114 * x264 * x8,
            x110 * x269 * x8,
            x108 * x274 * x8,
            x104 * x279,
            x282 * x47 * x79,
            x282 * x45 * x94,
            x283 * x45 * x79,
            x282 * x82 * x9,
            x280 * x283 * x66,
            x284 * x4,
            x282 * x8 * x97,
            x283 * x8 * x82,
            x284 * x66,
            x262
            * (
                x115 * x278
                + x3
                * (
                    x143 * (x276 + x277)
                    + x19 * (x125 + x256 + 4 * x257 + 3 * x271)
                    + 2 * x259
                    + 2 * x260
                    + 3 * x270
                    + 3 * x273
                )
            ),
        ]
    )


def diag_quadrupole3d_34(a, A, b, B, C):
    """Cartesian 3D (fg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x30 = x3 * (x21 + x28)
    x31 = 2 * x29 + 3 * x30
    x32 = x3 * (2 * x18 + 3 * x25 + x31)
    x33 = x25 + x30
    x34 = x33 * x4
    x35 = 4 * x13
    x36 = x10 * x35
    x37 = 2 * x12 + x36
    x38 = x3 * (2 * x16 + 2 * x23 + x37)
    x39 = x34 + x38
    x40 = x39 * x4
    x41 = x32 + x40
    x42 = x2 * x41
    x43 = x18 + x29
    x44 = x4 * x43
    x45 = 3 * x12
    x46 = x19 * x9
    x47 = x13 + x26
    x48 = x4 * x47
    x49 = x46 + x48
    x50 = x3 * (3 * x16 + x45 + x49)
    x51 = 4 * x38 + 2 * x50
    x52 = x3 * (4 * x34 + 2 * x44 + x51)
    x53 = x42 + x52
    x54 = 4 * x29
    x55 = x3 * (3 * x26 + x27)
    x56 = x4 * x49
    x57 = x55 + x56
    x58 = x3 * (4 * x18 + x54 + x57)
    x59 = x2 * x39
    x60 = x2 * (x44 + x50)
    x61 = x2 * x53 + x3 * (5 * x32 + x40 + 2 * x58 + 4 * x59 + 2 * x60)
    x62 = 2 * x2
    x63 = x32 + x59
    x64 = x2 * x63
    x65 = x2 * x43
    x66 = 8 * x13 * x4
    x67 = x3 * (4 * x48 + x66)
    x68 = x2 * x57
    x69 = x67 + x68
    x70 = x2 * x33
    x71 = x3 * (x34 + x51 + 2 * x65 + 3 * x70)
    x72 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x73 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x74 = numpy.pi * x0 * x73
    x75 = x72 * x74
    x76 = -x0 * (a * A[1] + b * B[1])
    x77 = -x76 - B[1]
    x78 = x64 + x71
    x79 = x38 + x70
    x80 = x2 * x79
    x81 = x17 * x2
    x82 = x2 * x49
    x83 = x55 + x82
    x84 = x2 * x24
    x85 = 2 * x84
    x86 = x3 * (x25 + x31 + 2 * x81 + x85)
    x87 = x75 * (
        x2 * x78
        + x3
        * (
            x19 * (x18 + x54 + 3 * x81 + x83)
            + 2 * x32
            + 2 * x59
            + x62 * (x50 + x65)
            + 3 * x80
            + 3 * x86
        )
    )
    x88 = -x0 * (a * A[2] + b * B[2])
    x89 = -x88 - B[2]
    x90 = x7 * x72
    x91 = x3 * x90
    x92 = x77 ** 2 * x90
    x93 = x91 + x92
    x94 = x80 + x86
    x95 = x30 + x84
    x96 = x2 * x95
    x97 = x2 * x22
    x98 = x15 * x2
    x99 = 2 * x98
    x100 = x3 * (x23 + x37 + x97 + x99)
    x101 = x2 * x47
    x102 = x101 + x46
    x103 = x2 * x94 + x3 * (
        2 * x100
        + x19 * (x102 + x16 + x45 + x99)
        + 2 * x38
        + x62 * (x29 + x81)
        + 2 * x70
        + 2 * x96
    )
    x104 = x7 * x73
    x105 = x75 * x89
    x106 = x104 * x3
    x107 = x104 * x89 ** 2
    x108 = x106 + x107
    x109 = x77 * x90
    x110 = x109 * x19
    x111 = x77 * x93
    x112 = x110 + x111
    x113 = x100 + x96
    x114 = x2 * x9
    x115 = x11 * x2
    x116 = x20 + x97
    x117 = x116 * x2 + x3 * (x11 * x62 + x21 + x27)
    x118 = x113 * x2 + x3 * (
        x117 + x19 * (x114 + x115 + x14 + x27) + 2 * x30 + x62 * (x12 + x98) + x85
    )
    x119 = x104 * x89
    x120 = x119 * x19
    x121 = x108 * x89
    x122 = x120 + x121
    x123 = 3 * x91
    x124 = x3 * (x123 + 3 * x92)
    x125 = x112 * x77
    x126 = x124 + x125
    x127 = x2 * x8
    x128 = x117 * x2 + x3 * (x19 * (x11 + x127) + x36 + x62 * (x115 + x13) + 2 * x97)
    x129 = 3 * x106
    x130 = x3 * (3 * x107 + x129)
    x131 = x122 * x89
    x132 = x130 + x131
    x133 = -x76 - A[1]
    x134 = x61 * x75
    x135 = x133 * x90
    x136 = x135 * x77
    x137 = x136 + x91
    x138 = x133 * x93
    x139 = x110 + x138
    x140 = x112 * x133
    x141 = x124 + x140
    x142 = x77 * x91
    x143 = 8 * x142
    x144 = x3 * (4 * x111 + x143)
    x145 = x126 * x133
    x146 = x144 + x145
    x147 = -x88 - A[2]
    x148 = x147 * x75
    x149 = x104 * x147
    x150 = x149 * x89
    x151 = x106 + x150
    x152 = x108 * x147
    x153 = x120 + x152
    x154 = x122 * x147
    x155 = x130 + x154
    x156 = x106 * x89
    x157 = 8 * x156
    x158 = x3 * (4 * x121 + x157)
    x159 = x132 * x147
    x160 = x158 + x159
    x161 = x133 ** 2 * x90
    x162 = x161 + x91
    x163 = x3 * (x109 + x135)
    x164 = x133 * x137
    x165 = x163 + x164
    x166 = 2 * x133
    x167 = x109 * x166 + x123
    x168 = x3 * (x167 + x92)
    x169 = x133 * x139
    x170 = x168 + x169
    x171 = x3 * (x111 + 3 * x138 + x143)
    x172 = x133 * x141
    x173 = x171 + x172
    x174 = x133 * x146 + x3 * (5 * x124 + x125 + 4 * x140)
    x175 = x104 * x147 ** 2
    x176 = x106 + x175
    x177 = x3 * (x119 + x149)
    x178 = x147 * x151
    x179 = x177 + x178
    x180 = 2 * x147
    x181 = x119 * x180 + x129
    x182 = x3 * (x107 + x181)
    x183 = x147 * x153
    x184 = x182 + x183
    x185 = x3 * (x121 + 3 * x152 + x157)
    x186 = x147 * x155
    x187 = x185 + x186
    x188 = x147 * x160 + x3 * (5 * x130 + x131 + 4 * x154)
    x189 = x133 * x162 + x135 * x19
    x190 = x133 * x165 + x3 * (x161 + x167)
    x191 = x133 * x170 + x3 * (2 * x138 + 4 * x142 + 2 * x163 + 2 * x164)
    x192 = x133 * x173 + x3 * (2 * x124 + 2 * x140 + 3 * x168 + 3 * x169)
    x193 = x133 * x174 + x3 * (2 * x144 + 2 * x145 + 4 * x171 + 4 * x172)
    x194 = x147 * x176 + x149 * x19
    x195 = x147 * x179 + x3 * (x175 + x181)
    x196 = x147 * x184 + x3 * (2 * x152 + 4 * x156 + 2 * x177 + 2 * x178)
    x197 = x147 * x187 + x3 * (2 * x130 + 2 * x154 + 3 * x182 + 3 * x183)
    x198 = x147 * x188 + x3 * (2 * x158 + 2 * x159 + 4 * x185 + 4 * x186)
    x199 = -x76 - C[1]
    x200 = x199 ** 2 * x90
    x201 = x200 + x91
    x202 = x2 * x69 + x3 * (5 * x55 + x56 + 4 * x82)
    x203 = x3 * (3 * x101 + x48 + x66)
    x204 = x2 * x83
    x205 = x2 * x202 + x3 * (4 * x203 + 4 * x204 + 2 * x67 + 2 * x68)
    x206 = x199 * x90
    x207 = x19 * x206
    x208 = x201 * x77
    x209 = x207 + x208
    x210 = x203 + x204
    x211 = x27 + x62 * x9
    x212 = x3 * (x211 + x26)
    x213 = x102 * x2
    x214 = x2 * x210 + x3 * (3 * x212 + 3 * x213 + 2 * x55 + 2 * x82)
    x215 = x109 * x199
    x216 = x123 + 2 * x215
    x217 = x3 * (x200 + x216)
    x218 = x209 * x77
    x219 = x217 + x218
    x220 = x212 + x213
    x221 = x3 * (x127 + x9)
    x222 = x114 + x13
    x223 = x2 * x222
    x224 = x2 * x220 + x3 * (2 * x101 + 2 * x221 + 2 * x223 + x35 * x4)
    x225 = x219 * x77
    x226 = x215 + x91
    x227 = x226 * x77
    x228 = x3 * (x109 + x206)
    x229 = 4 * x199 * x91
    x230 = 2 * x228 + x229
    x231 = x3 * (2 * x208 + 2 * x227 + x230)
    x232 = x225 + x231
    x233 = x2 ** 2 * x8
    x234 = x221 + x223
    x235 = x2 * x234 + x3 * (x211 + x233)
    x236 = x227 + x228
    x237 = x236 * x77
    x238 = x3 * (x216 + x92)
    x239 = 3 * x217 + 2 * x238
    x240 = x3 * (3 * x218 + 2 * x237 + x239)
    x241 = x232 * x77
    x242 = x240 + x241
    x243 = x13 + x233
    x244 = x127 * x19 + x2 * x243
    x245 = x133 * x201
    x246 = x207 + x245
    x247 = x133 * x209
    x248 = x217 + x247
    x249 = x133 * x219
    x250 = x231 + x249
    x251 = x133 * x232
    x252 = x240 + x251
    x253 = x133 * x242
    x254 = x237 + x238
    x255 = x254 * x77
    x256 = 3 * x228
    x257 = x3 * (x112 + 3 * x227 + x256)
    x258 = 4 * x231 + 2 * x257
    x259 = x3 * (4 * x225 + 2 * x255 + x258)
    x260 = x253 + x259
    x261 = x133 * x246 + x3 * (x123 + x166 * x206 + x200)
    x262 = x133 * x248
    x263 = x133 * x226
    x264 = 2 * x263
    x265 = x3 * (x208 + x230 + x245 + x264)
    x266 = x262 + x265
    x267 = x133 * x250
    x268 = x133 * x236
    x269 = 2 * x247
    x270 = x3 * (x218 + x239 + 2 * x268 + x269)
    x271 = x267 + x270
    x272 = x133 * x252
    x273 = x133 * x254
    x274 = x3 * (x225 + 3 * x249 + x258 + 2 * x273)
    x275 = x272 + x274
    x276 = 4 * x238
    x277 = x3 * (x126 + 4 * x237 + x276)
    x278 = x133 * (x255 + x257)
    x279 = x133 * x260 + x3 * (5 * x240 + x241 + 4 * x251 + 2 * x277 + 2 * x278)
    x280 = x6 * x74
    x281 = x279 * x280
    x282 = x2 * x280
    x283 = x135 * x199
    x284 = x133 * x261 + x3 * (
        x166 * (x283 + x91) + x19 * (x135 + x206) + x229 + 2 * x245
    )
    x285 = x133 * x266 + x3 * (
        x166 * (x228 + x263) + x19 * (x123 + x136 + x215 + x283) + 2 * x217 + x261 + x269
    )
    x286 = x133 * x271 + x3 * (
        x166 * (x238 + x268)
        + x19 * (x139 + x227 + x256 + x264)
        + 2 * x231
        + 2 * x249
        + 2 * x262
        + 2 * x265
    )
    x287 = x280 * (
        x133 * x275
        + x3
        * (
            x166 * (x257 + x273)
            + x19 * (x141 + x237 + 3 * x268 + x276)
            + 2 * x240
            + 2 * x251
            + 3 * x267
            + 3 * x270
        )
    )
    x288 = x280 * x4
    x289 = -x88 - C[2]
    x290 = x104 * x289 ** 2
    x291 = x106 + x290
    x292 = x104 * x289
    x293 = x19 * x292
    x294 = x291 * x89
    x295 = x293 + x294
    x296 = x119 * x289
    x297 = x129 + 2 * x296
    x298 = x3 * (x290 + x297)
    x299 = x295 * x89
    x300 = x298 + x299
    x301 = x300 * x89
    x302 = x106 + x296
    x303 = x302 * x89
    x304 = x3 * (x119 + x292)
    x305 = 4 * x106 * x289
    x306 = 2 * x304 + x305
    x307 = x3 * (2 * x294 + 2 * x303 + x306)
    x308 = x301 + x307
    x309 = x303 + x304
    x310 = x309 * x89
    x311 = x3 * (x107 + x297)
    x312 = 3 * x298 + 2 * x311
    x313 = x3 * (3 * x299 + 2 * x310 + x312)
    x314 = x308 * x89
    x315 = x313 + x314
    x316 = x147 * x291
    x317 = x293 + x316
    x318 = x147 * x295
    x319 = x298 + x318
    x320 = x147 * x300
    x321 = x307 + x320
    x322 = x147 * x308
    x323 = x313 + x322
    x324 = x147 * x315
    x325 = x310 + x311
    x326 = x325 * x89
    x327 = 3 * x304
    x328 = x3 * (x122 + 3 * x303 + x327)
    x329 = 4 * x307 + 2 * x328
    x330 = x3 * (4 * x301 + 2 * x326 + x329)
    x331 = x324 + x330
    x332 = numpy.pi * x0 * x6 * x72
    x333 = x2 * x332
    x334 = x147 * x317 + x3 * (x129 + x180 * x292 + x290)
    x335 = x147 * x319
    x336 = x147 * x302
    x337 = 2 * x336
    x338 = x3 * (x294 + x306 + x316 + x337)
    x339 = x335 + x338
    x340 = x147 * x321
    x341 = x147 * x309
    x342 = 2 * x318
    x343 = x3 * (x299 + x312 + 2 * x341 + x342)
    x344 = x340 + x343
    x345 = x147 * x323
    x346 = x147 * x325
    x347 = x3 * (x301 + 3 * x320 + x329 + 2 * x346)
    x348 = x345 + x347
    x349 = 4 * x311
    x350 = x3 * (x132 + 4 * x310 + x349)
    x351 = x147 * (x326 + x328)
    x352 = x147 * x331 + x3 * (5 * x313 + x314 + 4 * x322 + 2 * x350 + 2 * x351)
    x353 = x332 * x352
    x354 = x332 * x4
    x355 = x149 * x289
    x356 = x147 * x334 + x3 * (
        x180 * (x106 + x355) + x19 * (x149 + x292) + x305 + 2 * x316
    )
    x357 = x147 * x339 + x3 * (
        x180 * (x304 + x336) + x19 * (x129 + x150 + x296 + x355) + 2 * x298 + x334 + x342
    )
    x358 = x147 * x344 + x3 * (
        x180 * (x311 + x341)
        + x19 * (x153 + x303 + x327 + x337)
        + 2 * x307
        + 2 * x320
        + 2 * x335
        + 2 * x338
    )
    x359 = x332 * (
        x147 * x348
        + x3
        * (
            x180 * (x328 + x346)
            + x19 * (x155 + x310 + 3 * x341 + x349)
            + 2 * x313
            + 2 * x322
            + 3 * x340
            + 3 * x343
        )
    )

    # 450 item(s)
    return numpy.array(
        [
            x75
            * (
                x2 * x61
                + x3
                * (
                    x19 * (x44 + 5 * x50 + 4 * x65 + x69)
                    + 2 * x42
                    + 2 * x52
                    + x62 * (x58 + x60)
                    + 4 * x64
                    + 4 * x71
                )
            ),
            x77 * x87,
            x87 * x89,
            x103 * x104 * x93,
            x103 * x105 * x77,
            x103 * x108 * x90,
            x104 * x112 * x118,
            x118 * x119 * x93,
            x108 * x109 * x118,
            x118 * x122 * x90,
            x104 * x126 * x128,
            x112 * x119 * x128,
            x108 * x128 * x93,
            x109 * x122 * x128,
            x128 * x132 * x90,
            x133 * x134,
            x104 * x137 * x78,
            x105 * x133 * x78,
            x104 * x139 * x94,
            x119 * x137 * x94,
            x108 * x135 * x94,
            x104 * x113 * x141,
            x113 * x119 * x139,
            x108 * x113 * x137,
            x113 * x122 * x135,
            x104 * x117 * x146,
            x117 * x119 * x141,
            x108 * x117 * x139,
            x117 * x122 * x137,
            x117 * x132 * x135,
            x134 * x147,
            x148 * x77 * x78,
            x151 * x78 * x90,
            x149 * x93 * x94,
            x109 * x151 * x94,
            x153 * x90 * x94,
            x112 * x113 * x149,
            x113 * x151 * x93,
            x109 * x113 * x153,
            x113 * x155 * x90,
            x117 * x126 * x149,
            x112 * x117 * x151,
            x117 * x153 * x93,
            x109 * x117 * x155,
            x117 * x160 * x90,
            x104 * x162 * x53,
            x104 * x165 * x63,
            x119 * x162 * x63,
            x104 * x170 * x79,
            x119 * x165 * x79,
            x108 * x162 * x79,
            x104 * x173 * x95,
            x119 * x170 * x95,
            x108 * x165 * x95,
            x122 * x162 * x95,
            x104 * x116 * x174,
            x116 * x119 * x173,
            x108 * x116 * x170,
            x116 * x122 * x165,
            x116 * x132 * x162,
            x133 * x148 * x53,
            x137 * x149 * x63,
            x135 * x151 * x63,
            x139 * x149 * x79,
            x137 * x151 * x79,
            x135 * x153 * x79,
            x141 * x149 * x95,
            x139 * x151 * x95,
            x137 * x153 * x95,
            x135 * x155 * x95,
            x116 * x146 * x149,
            x116 * x141 * x151,
            x116 * x139 * x153,
            x116 * x137 * x155,
            x116 * x135 * x160,
            x176 * x53 * x90,
            x109 * x176 * x63,
            x179 * x63 * x90,
            x176 * x79 * x93,
            x109 * x179 * x79,
            x184 * x79 * x90,
            x112 * x176 * x95,
            x179 * x93 * x95,
            x109 * x184 * x95,
            x187 * x90 * x95,
            x116 * x126 * x176,
            x112 * x116 * x179,
            x116 * x184 * x93,
            x109 * x116 * x187,
            x116 * x188 * x90,
            x104 * x189 * x41,
            x104 * x190 * x39,
            x119 * x189 * x39,
            x104 * x191 * x33,
            x119 * x190 * x33,
            x108 * x189 * x33,
            x104 * x192 * x24,
            x119 * x191 * x24,
            x108 * x190 * x24,
            x122 * x189 * x24,
            x104 * x193 * x22,
            x119 * x192 * x22,
            x108 * x191 * x22,
            x122 * x190 * x22,
            x132 * x189 * x22,
            x149 * x162 * x41,
            x149 * x165 * x39,
            x151 * x162 * x39,
            x149 * x170 * x33,
            x151 * x165 * x33,
            x153 * x162 * x33,
            x149 * x173 * x24,
            x151 * x170 * x24,
            x153 * x165 * x24,
            x155 * x162 * x24,
            x149 * x174 * x22,
            x151 * x173 * x22,
            x153 * x170 * x22,
            x155 * x165 * x22,
            x160 * x162 * x22,
            x135 * x176 * x41,
            x137 * x176 * x39,
            x135 * x179 * x39,
            x139 * x176 * x33,
            x137 * x179 * x33,
            x135 * x184 * x33,
            x141 * x176 * x24,
            x139 * x179 * x24,
            x137 * x184 * x24,
            x135 * x187 * x24,
            x146 * x176 * x22,
            x141 * x179 * x22,
            x139 * x184 * x22,
            x137 * x187 * x22,
            x135 * x188 * x22,
            x194 * x41 * x90,
            x109 * x194 * x39,
            x195 * x39 * x90,
            x194 * x33 * x93,
            x109 * x195 * x33,
            x196 * x33 * x90,
            x112 * x194 * x24,
            x195 * x24 * x93,
            x109 * x196 * x24,
            x197 * x24 * x90,
            x126 * x194 * x22,
            x112 * x195 * x22,
            x196 * x22 * x93,
            x109 * x197 * x22,
            x198 * x22 * x90,
            x104 * x201 * x205,
            x104 * x209 * x214,
            x119 * x201 * x214,
            x104 * x219 * x224,
            x119 * x209 * x224,
            x108 * x201 * x224,
            x104 * x232 * x235,
            x119 * x219 * x235,
            x108 * x209 * x235,
            x122 * x201 * x235,
            x104 * x242 * x244,
            x119 * x232 * x244,
            x108 * x219 * x244,
            x122 * x209 * x244,
            x132 * x201 * x244,
            x104 * x202 * x246,
            x104 * x210 * x248,
            x119 * x210 * x246,
            x104 * x220 * x250,
            x119 * x220 * x248,
            x108 * x220 * x246,
            x104 * x234 * x252,
            x119 * x234 * x250,
            x108 * x234 * x248,
            x122 * x234 * x246,
            x104 * x243 * x260,
            x119 * x243 * x252,
            x108 * x243 * x250,
            x122 * x243 * x248,
            x132 * x243 * x246,
            x149 * x201 * x202,
            x149 * x209 * x210,
            x151 * x201 * x210,
            x149 * x219 * x220,
            x151 * x209 * x220,
            x153 * x201 * x220,
            x149 * x232 * x234,
            x151 * x219 * x234,
            x153 * x209 * x234,
            x155 * x201 * x234,
            x149 * x242 * x243,
            x151 * x232 * x243,
            x153 * x219 * x243,
            x155 * x209 * x243,
            x160 * x201 * x243,
            x104 * x261 * x69,
            x104 * x266 * x83,
            x119 * x261 * x83,
            x102 * x104 * x271,
            x102 * x119 * x266,
            x102 * x108 * x261,
            x104 * x222 * x275,
            x119 * x222 * x271,
            x108 * x222 * x266,
            x122 * x222 * x261,
            x2 * x281,
            x275 * x282 * x89,
            x108 * x127 * x271,
            x122 * x127 * x266,
            x127 * x132 * x261,
            x149 * x246 * x69,
            x149 * x248 * x83,
            x151 * x246 * x83,
            x102 * x149 * x250,
            x102 * x151 * x248,
            x102 * x153 * x246,
            x149 * x222 * x252,
            x151 * x222 * x250,
            x153 * x222 * x248,
            x155 * x222 * x246,
            x147 * x260 * x282,
            x127 * x151 * x252,
            x127 * x153 * x250,
            x127 * x155 * x248,
            x127 * x160 * x246,
            x176 * x201 * x69,
            x176 * x209 * x83,
            x179 * x201 * x83,
            x102 * x176 * x219,
            x102 * x179 * x209,
            x102 * x184 * x201,
            x176 * x222 * x232,
            x179 * x219 * x222,
            x184 * x209 * x222,
            x187 * x201 * x222,
            x127 * x176 * x242,
            x127 * x179 * x232,
            x127 * x184 * x219,
            x127 * x187 * x209,
            x127 * x188 * x201,
            x104 * x284 * x57,
            x104 * x285 * x49,
            x119 * x284 * x49,
            x104 * x286 * x47,
            x119 * x285 * x47,
            x108 * x284 * x47,
            x287 * x4,
            x286 * x288 * x89,
            x108 * x285 * x9,
            x122 * x284 * x9,
            x280
            * (
                x133 * x279
                + x3
                * (
                    x166 * (x277 + x278)
                    + x19 * (x146 + x255 + 5 * x257 + 4 * x273)
                    + 2 * x253
                    + 2 * x259
                    + 4 * x272
                    + 4 * x274
                )
            ),
            x287 * x89,
            x108 * x286 * x8,
            x122 * x285 * x8,
            x132 * x284 * x8,
            x149 * x261 * x57,
            x149 * x266 * x49,
            x151 * x261 * x49,
            x149 * x271 * x47,
            x151 * x266 * x47,
            x153 * x261 * x47,
            x147 * x275 * x288,
            x151 * x271 * x9,
            x153 * x266 * x9,
            x155 * x261 * x9,
            x147 * x281,
            x151 * x275 * x8,
            x153 * x271 * x8,
            x155 * x266 * x8,
            x160 * x261 * x8,
            x176 * x246 * x57,
            x176 * x248 * x49,
            x179 * x246 * x49,
            x176 * x250 * x47,
            x179 * x248 * x47,
            x184 * x246 * x47,
            x176 * x252 * x9,
            x179 * x250 * x9,
            x184 * x248 * x9,
            x187 * x246 * x9,
            x176 * x260 * x8,
            x179 * x252 * x8,
            x184 * x250 * x8,
            x187 * x248 * x8,
            x188 * x246 * x8,
            x194 * x201 * x57,
            x194 * x209 * x49,
            x195 * x201 * x49,
            x194 * x219 * x47,
            x195 * x209 * x47,
            x196 * x201 * x47,
            x194 * x232 * x9,
            x195 * x219 * x9,
            x196 * x209 * x9,
            x197 * x201 * x9,
            x194 * x242 * x8,
            x195 * x232 * x8,
            x196 * x219 * x8,
            x197 * x209 * x8,
            x198 * x201 * x8,
            x205 * x291 * x90,
            x109 * x214 * x291,
            x214 * x295 * x90,
            x224 * x291 * x93,
            x109 * x224 * x295,
            x224 * x300 * x90,
            x112 * x235 * x291,
            x235 * x295 * x93,
            x109 * x235 * x300,
            x235 * x308 * x90,
            x126 * x244 * x291,
            x112 * x244 * x295,
            x244 * x300 * x93,
            x109 * x244 * x308,
            x244 * x315 * x90,
            x135 * x202 * x291,
            x137 * x210 * x291,
            x135 * x210 * x295,
            x139 * x220 * x291,
            x137 * x220 * x295,
            x135 * x220 * x300,
            x141 * x234 * x291,
            x139 * x234 * x295,
            x137 * x234 * x300,
            x135 * x234 * x308,
            x146 * x243 * x291,
            x141 * x243 * x295,
            x139 * x243 * x300,
            x137 * x243 * x308,
            x135 * x243 * x315,
            x202 * x317 * x90,
            x109 * x210 * x317,
            x210 * x319 * x90,
            x220 * x317 * x93,
            x109 * x220 * x319,
            x220 * x321 * x90,
            x112 * x234 * x317,
            x234 * x319 * x93,
            x109 * x234 * x321,
            x234 * x323 * x90,
            x126 * x243 * x317,
            x112 * x243 * x319,
            x243 * x321 * x93,
            x109 * x243 * x323,
            x243 * x331 * x90,
            x162 * x291 * x69,
            x165 * x291 * x83,
            x162 * x295 * x83,
            x102 * x170 * x291,
            x102 * x165 * x295,
            x102 * x162 * x300,
            x173 * x222 * x291,
            x170 * x222 * x295,
            x165 * x222 * x300,
            x162 * x222 * x308,
            x127 * x174 * x291,
            x127 * x173 * x295,
            x127 * x170 * x300,
            x127 * x165 * x308,
            x127 * x162 * x315,
            x135 * x317 * x69,
            x137 * x317 * x83,
            x135 * x319 * x83,
            x102 * x139 * x317,
            x102 * x137 * x319,
            x102 * x135 * x321,
            x141 * x222 * x317,
            x139 * x222 * x319,
            x137 * x222 * x321,
            x135 * x222 * x323,
            x127 * x146 * x317,
            x127 * x141 * x319,
            x127 * x139 * x321,
            x127 * x137 * x323,
            x133 * x331 * x333,
            x334 * x69 * x90,
            x109 * x334 * x83,
            x339 * x83 * x90,
            x102 * x334 * x93,
            x102 * x109 * x339,
            x102 * x344 * x90,
            x112 * x222 * x334,
            x222 * x339 * x93,
            x109 * x222 * x344,
            x222 * x348 * x90,
            x126 * x127 * x334,
            x112 * x127 * x339,
            x127 * x344 * x93,
            x333 * x348 * x77,
            x2 * x353,
            x189 * x291 * x57,
            x190 * x291 * x49,
            x189 * x295 * x49,
            x191 * x291 * x47,
            x190 * x295 * x47,
            x189 * x300 * x47,
            x192 * x291 * x9,
            x191 * x295 * x9,
            x190 * x300 * x9,
            x189 * x308 * x9,
            x193 * x291 * x8,
            x192 * x295 * x8,
            x191 * x300 * x8,
            x190 * x308 * x8,
            x189 * x315 * x8,
            x162 * x317 * x57,
            x165 * x317 * x49,
            x162 * x319 * x49,
            x170 * x317 * x47,
            x165 * x319 * x47,
            x162 * x321 * x47,
            x173 * x317 * x9,
            x170 * x319 * x9,
            x165 * x321 * x9,
            x162 * x323 * x9,
            x174 * x317 * x8,
            x173 * x319 * x8,
            x170 * x321 * x8,
            x165 * x323 * x8,
            x162 * x331 * x8,
            x135 * x334 * x57,
            x137 * x334 * x49,
            x135 * x339 * x49,
            x139 * x334 * x47,
            x137 * x339 * x47,
            x135 * x344 * x47,
            x141 * x334 * x9,
            x139 * x339 * x9,
            x137 * x344 * x9,
            x133 * x348 * x354,
            x146 * x334 * x8,
            x141 * x339 * x8,
            x139 * x344 * x8,
            x137 * x348 * x8,
            x133 * x353,
            x356 * x57 * x90,
            x109 * x356 * x49,
            x357 * x49 * x90,
            x356 * x47 * x93,
            x109 * x357 * x47,
            x358 * x47 * x90,
            x112 * x356 * x9,
            x357 * x9 * x93,
            x354 * x358 * x77,
            x359 * x4,
            x126 * x356 * x8,
            x112 * x357 * x8,
            x358 * x8 * x93,
            x359 * x77,
            x332
            * (
                x147 * x352
                + x3
                * (
                    x180 * (x350 + x351)
                    + x19 * (x160 + x326 + 5 * x328 + 4 * x346)
                    + 2 * x324
                    + 2 * x330
                    + 4 * x345
                    + 4 * x347
                )
            ),
        ]
    )


def diag_quadrupole3d_40(a, A, b, B, C):
    """Cartesian 3D (gs) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x11 = 2 * x3
    x12 = -x2 - C[0]
    x13 = x12 * x7
    x14 = x10 + x11 * x13
    x15 = 2 * x0
    x16 = x12 ** 2 * x7
    x17 = x0 * (x14 + x16)
    x18 = x16 + x9
    x19 = x18 * x3
    x20 = x13 * x15 + x19
    x21 = x20 * x3
    x22 = x3 * x7
    x23 = x0 * (x13 + x22)
    x24 = x3 * (x12 * x22 + x9)
    x25 = x17 + x21
    x26 = x0 * (4 * x0 * x13 + 2 * x19 + 2 * x23 + 2 * x24) + x25 * x3
    x27 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x28 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x29 = numpy.pi * x1 * x28
    x30 = x27 * x29
    x31 = -x1 * (a * A[1] + b * B[1])
    x32 = -x31 - A[1]
    x33 = x26 * x30
    x34 = -x1 * (a * A[2] + b * B[2])
    x35 = -x34 - A[2]
    x36 = x27 * x6
    x37 = x0 * x36
    x38 = x32 ** 2 * x36
    x39 = x37 + x38
    x40 = x28 * x6
    x41 = x0 * x40
    x42 = x35 ** 2 * x40
    x43 = x41 + x42
    x44 = x32 * x36
    x45 = x15 * x44 + x32 * x39
    x46 = x35 * x40
    x47 = x15 * x46 + x35 * x43
    x48 = 3 * x37
    x49 = x0 * (3 * x38 + x48) + x32 * x45
    x50 = 3 * x41
    x51 = x0 * (3 * x42 + x50) + x35 * x47
    x52 = -x31 - C[1]
    x53 = x36 * x52 ** 2
    x54 = x37 + x53
    x55 = x8 + x9
    x56 = x15 * x22 + x3 * x55
    x57 = x0 * (x10 + 3 * x8) + x3 * x56
    x58 = x36 * x52
    x59 = x32 * x54
    x60 = x15 * x58 + x59
    x61 = 2 * x32
    x62 = x48 + x58 * x61
    x63 = x0 * (x53 + x62)
    x64 = x32 * x60
    x65 = x63 + x64
    x66 = x0 * (x44 + x58)
    x67 = x32 * (x37 + x44 * x52)
    x68 = x0 * (4 * x37 * x52 + 2 * x59 + 2 * x66 + 2 * x67) + x32 * x65
    x69 = x29 * x5
    x70 = x68 * x69
    x71 = -x34 - C[2]
    x72 = x40 * x71 ** 2
    x73 = x41 + x72
    x74 = x40 * x71
    x75 = x35 * x73
    x76 = x15 * x74 + x75
    x77 = 2 * x35
    x78 = x50 + x74 * x77
    x79 = x0 * (x72 + x78)
    x80 = x35 * x76
    x81 = x79 + x80
    x82 = numpy.pi * x1 * x27 * x5
    x83 = x0 * (x46 + x74)
    x84 = x35 * (x41 + x46 * x71)
    x85 = x0 * (4 * x41 * x71 + 2 * x75 + 2 * x83 + 2 * x84) + x35 * x81
    x86 = x82 * x85

    # 45 item(s)
    return numpy.array(
        [
            x30
            * (
                x0 * (x11 * (x23 + x24) + x15 * (x14 + x8) + 3 * x17 + 3 * x21) + x26 * x3
            ),
            x32 * x33,
            x33 * x35,
            x25 * x39 * x40,
            x25 * x30 * x32 * x35,
            x25 * x36 * x43,
            x20 * x40 * x45,
            x20 * x39 * x46,
            x20 * x43 * x44,
            x20 * x36 * x47,
            x18 * x40 * x49,
            x18 * x45 * x46,
            x18 * x39 * x43,
            x18 * x44 * x47,
            x18 * x36 * x51,
            x40 * x54 * x57,
            x40 * x56 * x60,
            x46 * x54 * x56,
            x40 * x55 * x65,
            x46 * x55 * x60,
            x43 * x54 * x55,
            x3 * x70,
            x3 * x35 * x65 * x69,
            x22 * x43 * x60,
            x22 * x47 * x54,
            x69
            * (
                x0 * (x15 * (x38 + x62) + x61 * (x66 + x67) + 3 * x63 + 3 * x64)
                + x32 * x68
            ),
            x35 * x70,
            x43 * x65 * x7,
            x47 * x60 * x7,
            x51 * x54 * x7,
            x36 * x57 * x73,
            x44 * x56 * x73,
            x36 * x56 * x76,
            x39 * x55 * x73,
            x44 * x55 * x76,
            x36 * x55 * x81,
            x22 * x45 * x73,
            x22 * x39 * x76,
            x3 * x32 * x81 * x82,
            x3 * x86,
            x49 * x7 * x73,
            x45 * x7 * x76,
            x39 * x7 * x81,
            x32 * x86,
            x82
            * (
                x0 * (x15 * (x42 + x78) + x77 * (x83 + x84) + 3 * x79 + 3 * x80)
                + x35 * x85
            ),
        ]
    )


def diag_quadrupole3d_41(a, A, b, B, C):
    """Cartesian 3D (gp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x17 = 2 * x3
    x18 = x10 * x17
    x19 = x12 + x13
    x20 = x19 * x4
    x21 = x18 + x20
    x22 = x2 * x21
    x23 = x16 + x22
    x24 = x2 * x23
    x25 = x19 * x2
    x26 = 4 * x10 * x3
    x27 = x4 * x9
    x28 = x3 * (x10 + x27)
    x29 = x2 * (x11 + x13)
    x30 = 2 * x28 + 2 * x29
    x31 = x3 * (x20 + x25 + x26 + x30)
    x32 = x24 + x31
    x33 = x2 * x7
    x34 = x33 * x8
    x35 = x34 * x4
    x36 = x34 * x5
    x37 = x3 * (x11 + x14 + x35 + x36)
    x38 = x2 * (x28 + x29)
    x39 = 2 * x2
    x40 = x10 * x39
    x41 = x3 * (x15 + x40)
    x42 = x18 + x25
    x43 = x2 * x42
    x44 = x41 + x43
    x45 = x2 * x32 + x3 * (2 * x16 + 2 * x22 + 2 * x37 + 2 * x38 + x44)
    x46 = x3 * (x10 + x34)
    x47 = x2 * (x13 + x36)
    x48 = x46 + x47
    x49 = x3 * (x27 + x34)
    x50 = x13 + x35
    x51 = x2 * x50
    x52 = x49 + x51
    x53 = x2 * x44 + x3 * (2 * x25 + x26 + 2 * x46 + 2 * x47)
    x54 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x55 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x56 = numpy.pi * x0 * x55
    x57 = x54 * x56
    x58 = -x0 * (a * A[1] + b * B[1])
    x59 = -x58 - B[1]
    x60 = x2 ** 2 * x9
    x61 = x14 + x60
    x62 = x57 * (x2 * x53 + x3 * (x17 * (x40 + x61) + x39 * x48 + 3 * x41 + 3 * x43))
    x63 = -x0 * (a * A[2] + b * B[2])
    x64 = -x63 - B[2]
    x65 = -x58 - A[1]
    x66 = x45 * x57
    x67 = x3 * x8
    x68 = x54 * x67
    x69 = x54 * x8
    x70 = x65 * x69
    x71 = x59 * x70
    x72 = x68 + x71
    x73 = x55 * x8
    x74 = x53 * x57
    x75 = -x63 - A[2]
    x76 = x55 * x67
    x77 = x73 * x75
    x78 = x64 * x77
    x79 = x76 + x78
    x80 = x65 ** 2 * x69
    x81 = x68 + x80
    x82 = x59 * x69
    x83 = x3 * (x70 + x82)
    x84 = x65 * x72
    x85 = x83 + x84
    x86 = x64 * x73
    x87 = x73 * x75 ** 2
    x88 = x76 + x87
    x89 = x3 * (x77 + x86)
    x90 = x75 * x79
    x91 = x89 + x90
    x92 = x17 * x70 + x65 * x81
    x93 = 2 * x65
    x94 = 3 * x68
    x95 = x80 + x94
    x96 = x3 * (x82 * x93 + x95) + x65 * x85
    x97 = x17 * x77 + x75 * x88
    x98 = 2 * x75
    x99 = 3 * x76
    x100 = x87 + x99
    x101 = x3 * (x100 + x86 * x98) + x75 * x91
    x102 = x3 * (3 * x80 + x94) + x65 * x92
    x103 = x3 * (3 * x83 + 3 * x84 + x92) + x65 * x96
    x104 = x3 * (3 * x87 + x99) + x75 * x97
    x105 = x101 * x75 + x3 * (3 * x89 + 3 * x90 + x97)
    x106 = -x58 - C[1]
    x107 = x106 ** 2 * x69
    x108 = x107 + x68
    x109 = x13 + x60
    x110 = x109 * x2 + x17 * x34
    x111 = x2 * x52 + x3 * (x27 * x39 + x61)
    x112 = x111 * x2 + x3 * (x110 + 3 * x49 + 3 * x51)
    x113 = x106 * x69
    x114 = x113 * x17
    x115 = x108 * x59
    x116 = x114 + x115
    x117 = x110 * x2 + x3 * (x14 + 3 * x60)
    x118 = x108 * x65
    x119 = x114 + x118
    x120 = x106 * x82
    x121 = x107 + x94
    x122 = x3 * (2 * x120 + x121)
    x123 = x116 * x65
    x124 = x122 + x123
    x125 = x113 * x93
    x126 = x3 * (x121 + x125)
    x127 = x119 * x65
    x128 = x126 + x127
    x129 = x124 * x65
    x130 = 4 * x106 * x68
    x131 = x3 * (x113 + x82)
    x132 = x65 * (x120 + x68)
    x133 = 2 * x131 + 2 * x132
    x134 = x3 * (x115 + x118 + x130 + x133)
    x135 = x129 + x134
    x136 = x3 * (x113 + x70)
    x137 = x106 * x70
    x138 = x65 * (x137 + x68)
    x139 = x128 * x65 + x3 * (2 * x118 + x130 + 2 * x136 + 2 * x138)
    x140 = x3 * (x120 + x137 + x71 + x94)
    x141 = x65 * (x131 + x132)
    x142 = x135 * x65 + x3 * (2 * x122 + 2 * x123 + x128 + 2 * x140 + 2 * x141)
    x143 = x33 * x56
    x144 = x136 + x138
    x145 = x56 * x7
    x146 = x145 * (
        x139 * x65 + x3 * (3 * x126 + 3 * x127 + x144 * x93 + x17 * (x125 + x95))
    )
    x147 = x145 * x75
    x148 = -x63 - C[2]
    x149 = x148 ** 2 * x73
    x150 = x149 + x76
    x151 = x148 * x73
    x152 = x151 * x17
    x153 = x150 * x64
    x154 = x152 + x153
    x155 = x150 * x75
    x156 = x152 + x155
    x157 = x148 * x86
    x158 = x149 + x99
    x159 = x3 * (2 * x157 + x158)
    x160 = x154 * x75
    x161 = x159 + x160
    x162 = x151 * x98
    x163 = x3 * (x158 + x162)
    x164 = x156 * x75
    x165 = x163 + x164
    x166 = x161 * x75
    x167 = 4 * x148 * x76
    x168 = x3 * (x151 + x86)
    x169 = x75 * (x157 + x76)
    x170 = 2 * x168 + 2 * x169
    x171 = x3 * (x153 + x155 + x167 + x170)
    x172 = x166 + x171
    x173 = numpy.pi * x0 * x54
    x174 = x173 * x33
    x175 = x3 * (x151 + x77)
    x176 = x148 * x77
    x177 = x75 * (x176 + x76)
    x178 = x165 * x75 + x3 * (2 * x155 + x167 + 2 * x175 + 2 * x177)
    x179 = x3 * (x157 + x176 + x78 + x99)
    x180 = x75 * (x168 + x169)
    x181 = x172 * x75 + x3 * (2 * x159 + 2 * x160 + x165 + 2 * x179 + 2 * x180)
    x182 = x173 * x7
    x183 = x182 * x65
    x184 = x175 + x177
    x185 = x182 * (
        x178 * x75 + x3 * (3 * x163 + 3 * x164 + x17 * (x100 + x162) + x184 * x98)
    )

    # 135 item(s)
    return numpy.array(
        [
            x57
            * (
                x2 * x45
                + x3
                * (x17 * (x30 + x48 + x52) + 3 * x24 + 3 * x31 + x39 * (x37 + x38) + x53)
            ),
            x59 * x62,
            x62 * x64,
            x65 * x66,
            x53 * x72 * x73,
            x64 * x65 * x74,
            x66 * x75,
            x59 * x74 * x75,
            x53 * x69 * x79,
            x32 * x73 * x81,
            x44 * x73 * x85,
            x44 * x81 * x86,
            x32 * x57 * x65 * x75,
            x44 * x72 * x77,
            x44 * x70 * x79,
            x32 * x69 * x88,
            x44 * x82 * x88,
            x44 * x69 * x91,
            x23 * x73 * x92,
            x42 * x73 * x96,
            x42 * x86 * x92,
            x23 * x77 * x81,
            x42 * x77 * x85,
            x42 * x79 * x81,
            x23 * x70 * x88,
            x42 * x72 * x88,
            x42 * x70 * x91,
            x23 * x69 * x97,
            x42 * x82 * x97,
            x101 * x42 * x69,
            x102 * x21 * x73,
            x103 * x19 * x73,
            x102 * x19 * x86,
            x21 * x77 * x92,
            x19 * x77 * x96,
            x19 * x79 * x92,
            x21 * x81 * x88,
            x19 * x85 * x88,
            x19 * x81 * x91,
            x21 * x70 * x97,
            x19 * x72 * x97,
            x101 * x19 * x70,
            x104 * x21 * x69,
            x104 * x19 * x82,
            x105 * x19 * x69,
            x108 * x112 * x73,
            x116 * x117 * x73,
            x108 * x117 * x86,
            x111 * x119 * x73,
            x110 * x124 * x73,
            x110 * x119 * x86,
            x108 * x111 * x77,
            x110 * x116 * x77,
            x108 * x110 * x79,
            x128 * x52 * x73,
            x109 * x135 * x73,
            x109 * x128 * x86,
            x119 * x52 * x77,
            x109 * x124 * x77,
            x109 * x119 * x79,
            x108 * x52 * x88,
            x109 * x116 * x88,
            x108 * x109 * x91,
            x139 * x50 * x73,
            x142 * x143,
            x139 * x143 * x64,
            x128 * x50 * x77,
            x135 * x143 * x75,
            x128 * x34 * x79,
            x119 * x50 * x88,
            x124 * x34 * x88,
            x119 * x34 * x91,
            x108 * x50 * x97,
            x116 * x34 * x97,
            x101 * x108 * x34,
            x146 * x4,
            x145
            * (
                x142 * x65
                + x3
                * (
                    3 * x129
                    + 3 * x134
                    + x139
                    + x17 * (x133 + x144 + x85)
                    + x93 * (x140 + x141)
                )
            ),
            x146 * x64,
            x139 * x147 * x4,
            x142 * x147,
            x139 * x79 * x9,
            x128 * x27 * x88,
            x135 * x88 * x9,
            x128 * x9 * x91,
            x119 * x27 * x97,
            x124 * x9 * x97,
            x101 * x119 * x9,
            x104 * x108 * x27,
            x104 * x116 * x9,
            x105 * x108 * x9,
            x112 * x150 * x69,
            x117 * x150 * x82,
            x117 * x154 * x69,
            x111 * x150 * x70,
            x110 * x150 * x72,
            x110 * x154 * x70,
            x111 * x156 * x69,
            x110 * x156 * x82,
            x110 * x161 * x69,
            x150 * x52 * x81,
            x109 * x150 * x85,
            x109 * x154 * x81,
            x156 * x52 * x70,
            x109 * x156 * x72,
            x109 * x161 * x70,
            x165 * x52 * x69,
            x109 * x165 * x82,
            x109 * x172 * x69,
            x150 * x50 * x92,
            x150 * x34 * x96,
            x154 * x34 * x92,
            x156 * x50 * x81,
            x156 * x34 * x85,
            x161 * x34 * x81,
            x165 * x50 * x70,
            x165 * x34 * x72,
            x172 * x174 * x65,
            x178 * x50 * x69,
            x174 * x178 * x59,
            x174 * x181,
            x102 * x150 * x27,
            x103 * x150 * x9,
            x102 * x154 * x9,
            x156 * x27 * x92,
            x156 * x9 * x96,
            x161 * x9 * x92,
            x165 * x27 * x81,
            x165 * x85 * x9,
            x172 * x81 * x9,
            x178 * x183 * x4,
            x178 * x72 * x9,
            x181 * x183,
            x185 * x4,
            x185 * x59,
            x182
            * (
                x181 * x75
                + x3
                * (
                    3 * x166
                    + x17 * (x170 + x184 + x91)
                    + 3 * x171
                    + x178
                    + x98 * (x179 + x180)
                )
            ),
        ]
    )


def diag_quadrupole3d_42(a, A, b, B, C):
    """Cartesian 3D (gd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x18 = 2 * x3
    x19 = x4 * x9
    x20 = x18 * x19
    x21 = x13 + x14
    x22 = x21 * x5
    x23 = x20 + x22
    x24 = x23 * x5
    x25 = x17 + x24
    x26 = x2 * x25
    x27 = x11 + x14
    x28 = x27 * x5
    x29 = x3 * (x10 + x19)
    x30 = 2 * x29
    x31 = 4 * x14
    x32 = x31 * x4
    x33 = x30 + x32
    x34 = x3 * (2 * x22 + 2 * x28 + x33)
    x35 = x26 + x34
    x36 = x2 * x35
    x37 = x2 * x23
    x38 = 2 * x37
    x39 = x5 ** 2 * x9
    x40 = x15 + x39
    x41 = x3 * (x12 + x40)
    x42 = x2 * (x28 + x29)
    x43 = 2 * x41 + 2 * x42
    x44 = x3 * (3 * x17 + x24 + x38 + x43)
    x45 = x36 + x44
    x46 = x17 + x37
    x47 = x2 * x46
    x48 = x2 * (x41 + x42)
    x49 = x2 * x21
    x50 = x2 * x27
    x51 = 2 * x50
    x52 = x3 * (x22 + x33 + x49 + x51)
    x53 = x14 + x39
    x54 = x2 * x53
    x55 = x10 * x18 + x54
    x56 = x3 * (x28 + 3 * x29 + x51 + x55)
    x57 = x2 * x45 + x3 * (2 * x26 + 2 * x34 + 2 * x47 + 2 * x48 + 2 * x52 + 2 * x56)
    x58 = x47 + x52
    x59 = x2 * x58
    x60 = 2 * x2
    x61 = x19 * x60
    x62 = x3 * (x16 + x61)
    x63 = x20 + x49
    x64 = x2 * x63
    x65 = x62 + x64
    x66 = x10 * x2
    x67 = x19 * x2
    x68 = x3 * (x11 + x15 + x66 + x67)
    x69 = x2 * (x29 + x50)
    x70 = 2 * x68 + 2 * x69
    x71 = x3 * (2 * x17 + x38 + x65 + x70)
    x72 = x10 * x60
    x73 = x3 * (x40 + x72)
    x74 = x2 * x55
    x75 = x73 + x74
    x76 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x77 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x78 = numpy.pi * x0 * x77
    x79 = x76 * x78
    x80 = -x0 * (a * A[1] + b * B[1])
    x81 = -x80 - B[1]
    x82 = x59 + x71
    x83 = x2 * x9
    x84 = x3 * (x19 + x83)
    x85 = x2 * (x14 + x67)
    x86 = x84 + x85
    x87 = x3 * (x10 + x83)
    x88 = x14 + x66
    x89 = x2 * x88
    x90 = x87 + x89
    x91 = x2 * x65 + x3 * (x32 + 2 * x49 + 2 * x84 + 2 * x85)
    x92 = x79 * (
        x2 * x82
        + x3
        * (x18 * (x30 + x51 + x86 + x90) + 3 * x47 + 3 * x52 + x60 * (x68 + x69) + x91)
    )
    x93 = -x0 * (a * A[2] + b * B[2])
    x94 = -x93 - B[2]
    x95 = x76 * x8
    x96 = x3 * x95
    x97 = x81 ** 2 * x95
    x98 = x96 + x97
    x99 = x2 ** 2 * x9
    x100 = x15 + x99
    x101 = x2 * x91 + x3 * (x18 * (x100 + x61) + x60 * x86 + 3 * x62 + 3 * x64)
    x102 = x77 * x8
    x103 = x79 * x94
    x104 = x102 * x3
    x105 = x102 * x94 ** 2
    x106 = x104 + x105
    x107 = -x80 - A[1]
    x108 = x57 * x79
    x109 = x107 * x95
    x110 = x109 * x81
    x111 = x110 + x96
    x112 = x81 * x95
    x113 = x107 * x98
    x114 = x112 * x18 + x113
    x115 = x102 * x94
    x116 = -x93 - A[2]
    x117 = x116 * x79
    x118 = x102 * x116
    x119 = x118 * x94
    x120 = x104 + x119
    x121 = x106 * x116
    x122 = x115 * x18 + x121
    x123 = x107 ** 2 * x95
    x124 = x123 + x96
    x125 = x3 * (x109 + x112)
    x126 = x107 * x111
    x127 = x125 + x126
    x128 = 3 * x96
    x129 = 2 * x107
    x130 = x112 * x129 + x128
    x131 = x3 * (x130 + x97)
    x132 = x107 * x114
    x133 = x131 + x132
    x134 = x102 * x116 ** 2
    x135 = x104 + x134
    x136 = x3 * (x115 + x118)
    x137 = x116 * x120
    x138 = x136 + x137
    x139 = 3 * x104
    x140 = 2 * x116
    x141 = x115 * x140 + x139
    x142 = x3 * (x105 + x141)
    x143 = x116 * x122
    x144 = x142 + x143
    x145 = x107 * x124 + x109 * x18
    x146 = x3 * (x123 + x130)
    x147 = x107 * x127
    x148 = x146 + x147
    x149 = 4 * x96
    x150 = x107 * x133 + x3 * (2 * x113 + 2 * x125 + 2 * x126 + x149 * x81)
    x151 = x116 * x135 + x118 * x18
    x152 = x3 * (x134 + x141)
    x153 = x116 * x138
    x154 = x152 + x153
    x155 = 4 * x104
    x156 = x116 * x144 + x3 * (2 * x121 + 2 * x136 + 2 * x137 + x155 * x94)
    x157 = x107 * x145 + x3 * (3 * x123 + x128)
    x158 = x107 * x148 + x3 * (3 * x125 + 3 * x126 + x145)
    x159 = x107 * x150 + x3 * (3 * x131 + 3 * x132 + 2 * x146 + 2 * x147)
    x160 = x116 * x151 + x3 * (3 * x134 + x139)
    x161 = x116 * x154 + x3 * (3 * x136 + 3 * x137 + x151)
    x162 = x116 * x156 + x3 * (3 * x142 + 3 * x143 + 2 * x152 + 2 * x153)
    x163 = -x80 - C[1]
    x164 = x163 ** 2 * x95
    x165 = x164 + x96
    x166 = x3 * (x100 + x72)
    x167 = x2 * x90
    x168 = x2 * x75 + x3 * (x31 * x5 + 2 * x54 + 2 * x87 + 2 * x89)
    x169 = x168 * x2 + x3 * (2 * x166 + 2 * x167 + 3 * x73 + 3 * x74)
    x170 = x163 * x95
    x171 = x170 * x18
    x172 = x165 * x81
    x173 = x171 + x172
    x174 = x14 + x99
    x175 = x174 * x2 + x18 * x83
    x176 = x166 + x167
    x177 = x176 * x2 + x3 * (x175 + 3 * x87 + 3 * x89)
    x178 = x112 * x163
    x179 = 2 * x178
    x180 = x128 + x164
    x181 = x3 * (x179 + x180)
    x182 = x173 * x81
    x183 = x181 + x182
    x184 = x175 * x2 + x3 * (x15 + 3 * x99)
    x185 = x107 * x165
    x186 = x171 + x185
    x187 = x107 * x173
    x188 = x181 + x187
    x189 = x107 * x183
    x190 = x178 + x96
    x191 = x190 * x81
    x192 = x3 * (x112 + x170)
    x193 = 2 * x192
    x194 = x149 * x163
    x195 = x193 + x194
    x196 = x3 * (2 * x172 + 2 * x191 + x195)
    x197 = x189 + x196
    x198 = x129 * x170
    x199 = x3 * (x180 + x198)
    x200 = x107 * x186
    x201 = x199 + x200
    x202 = x107 * x188
    x203 = x107 * x190
    x204 = 2 * x203
    x205 = x3 * (x172 + x185 + x195 + x204)
    x206 = x202 + x205
    x207 = x107 * x197
    x208 = 2 * x187
    x209 = x3 * (x128 + x179 + x97)
    x210 = x107 * (x191 + x192)
    x211 = 2 * x209 + 2 * x210
    x212 = x3 * (3 * x181 + x182 + x208 + x211)
    x213 = x207 + x212
    x214 = x3 * (x109 + x170)
    x215 = x109 * x163
    x216 = x107 * (x215 + x96)
    x217 = x107 * x201 + x3 * (2 * x185 + x194 + 2 * x214 + 2 * x216)
    x218 = x107 * x206
    x219 = x3 * (x110 + x128 + x178 + x215)
    x220 = x107 * (x192 + x203)
    x221 = 2 * x219 + 2 * x220
    x222 = x3 * (2 * x181 + x201 + x208 + x221)
    x223 = x218 + x222
    x224 = x107 * (x209 + x210)
    x225 = x3 * (x114 + x191 + 3 * x192 + x204)
    x226 = x107 * x213 + x3 * (
        2 * x189 + 2 * x196 + 2 * x202 + 2 * x205 + 2 * x224 + 2 * x225
    )
    x227 = x7 * x78
    x228 = x226 * x227
    x229 = x2 * x227
    x230 = x214 + x216
    x231 = x107 * x217 + x3 * (
        x129 * x230 + x18 * (x123 + x128 + x198) + 3 * x199 + 3 * x200
    )
    x232 = x227 * (
        x107 * x223
        + x3
        * (
            x129 * (x219 + x220)
            + x18 * (x127 + x193 + x204 + x230)
            + 3 * x202
            + 3 * x205
            + x217
        )
    )
    x233 = x227 * x5
    x234 = -x93 - C[2]
    x235 = x102 * x234 ** 2
    x236 = x104 + x235
    x237 = x102 * x234
    x238 = x18 * x237
    x239 = x236 * x94
    x240 = x238 + x239
    x241 = x115 * x234
    x242 = 2 * x241
    x243 = x139 + x235
    x244 = x3 * (x242 + x243)
    x245 = x240 * x94
    x246 = x244 + x245
    x247 = x116 * x236
    x248 = x238 + x247
    x249 = x116 * x240
    x250 = x244 + x249
    x251 = x116 * x246
    x252 = x104 + x241
    x253 = x252 * x94
    x254 = x3 * (x115 + x237)
    x255 = 2 * x254
    x256 = x155 * x234
    x257 = x255 + x256
    x258 = x3 * (2 * x239 + 2 * x253 + x257)
    x259 = x251 + x258
    x260 = x140 * x237
    x261 = x3 * (x243 + x260)
    x262 = x116 * x248
    x263 = x261 + x262
    x264 = x116 * x250
    x265 = x116 * x252
    x266 = 2 * x265
    x267 = x3 * (x239 + x247 + x257 + x266)
    x268 = x264 + x267
    x269 = x116 * x259
    x270 = 2 * x249
    x271 = x3 * (x105 + x139 + x242)
    x272 = x116 * (x253 + x254)
    x273 = 2 * x271 + 2 * x272
    x274 = x3 * (3 * x244 + x245 + x270 + x273)
    x275 = x269 + x274
    x276 = numpy.pi * x0 * x7 * x76
    x277 = x2 * x276
    x278 = x3 * (x118 + x237)
    x279 = x118 * x234
    x280 = x116 * (x104 + x279)
    x281 = x116 * x263 + x3 * (2 * x247 + x256 + 2 * x278 + 2 * x280)
    x282 = x116 * x268
    x283 = x3 * (x119 + x139 + x241 + x279)
    x284 = x116 * (x254 + x265)
    x285 = 2 * x283 + 2 * x284
    x286 = x3 * (2 * x244 + x263 + x270 + x285)
    x287 = x282 + x286
    x288 = x116 * (x271 + x272)
    x289 = x3 * (x122 + x253 + 3 * x254 + x266)
    x290 = x116 * x275 + x3 * (
        2 * x251 + 2 * x258 + 2 * x264 + 2 * x267 + 2 * x288 + 2 * x289
    )
    x291 = x276 * x290
    x292 = x276 * x5
    x293 = x278 + x280
    x294 = x116 * x281 + x3 * (
        x140 * x293 + x18 * (x134 + x139 + x260) + 3 * x261 + 3 * x262
    )
    x295 = x276 * (
        x116 * x287
        + x3
        * (
            x140 * (x283 + x284)
            + x18 * (x138 + x255 + x266 + x293)
            + 3 * x264
            + 3 * x267
            + x281
        )
    )

    # 270 item(s)
    return numpy.array(
        [
            x79
            * (
                x2 * x57
                + x3
                * (
                    x18 * (x43 + x70 + x75)
                    + 3 * x36
                    + 3 * x44
                    + 2 * x59
                    + x60 * (x48 + x56)
                    + 2 * x71
                )
            ),
            x81 * x92,
            x92 * x94,
            x101 * x102 * x98,
            x101 * x103 * x81,
            x101 * x106 * x95,
            x107 * x108,
            x102 * x111 * x82,
            x103 * x107 * x82,
            x102 * x114 * x91,
            x111 * x115 * x91,
            x106 * x109 * x91,
            x108 * x116,
            x117 * x81 * x82,
            x120 * x82 * x95,
            x118 * x91 * x98,
            x112 * x120 * x91,
            x122 * x91 * x95,
            x102 * x124 * x45,
            x102 * x127 * x58,
            x115 * x124 * x58,
            x102 * x133 * x65,
            x115 * x127 * x65,
            x106 * x124 * x65,
            x107 * x117 * x45,
            x111 * x118 * x58,
            x109 * x120 * x58,
            x114 * x118 * x65,
            x111 * x120 * x65,
            x109 * x122 * x65,
            x135 * x45 * x95,
            x112 * x135 * x58,
            x138 * x58 * x95,
            x135 * x65 * x98,
            x112 * x138 * x65,
            x144 * x65 * x95,
            x102 * x145 * x35,
            x102 * x148 * x46,
            x115 * x145 * x46,
            x102 * x150 * x63,
            x115 * x148 * x63,
            x106 * x145 * x63,
            x118 * x124 * x35,
            x118 * x127 * x46,
            x120 * x124 * x46,
            x118 * x133 * x63,
            x120 * x127 * x63,
            x122 * x124 * x63,
            x109 * x135 * x35,
            x111 * x135 * x46,
            x109 * x138 * x46,
            x114 * x135 * x63,
            x111 * x138 * x63,
            x109 * x144 * x63,
            x151 * x35 * x95,
            x112 * x151 * x46,
            x154 * x46 * x95,
            x151 * x63 * x98,
            x112 * x154 * x63,
            x156 * x63 * x95,
            x102 * x157 * x25,
            x102 * x158 * x23,
            x115 * x157 * x23,
            x102 * x159 * x21,
            x115 * x158 * x21,
            x106 * x157 * x21,
            x118 * x145 * x25,
            x118 * x148 * x23,
            x120 * x145 * x23,
            x118 * x150 * x21,
            x120 * x148 * x21,
            x122 * x145 * x21,
            x124 * x135 * x25,
            x127 * x135 * x23,
            x124 * x138 * x23,
            x133 * x135 * x21,
            x127 * x138 * x21,
            x124 * x144 * x21,
            x109 * x151 * x25,
            x111 * x151 * x23,
            x109 * x154 * x23,
            x114 * x151 * x21,
            x111 * x154 * x21,
            x109 * x156 * x21,
            x160 * x25 * x95,
            x112 * x160 * x23,
            x161 * x23 * x95,
            x160 * x21 * x98,
            x112 * x161 * x21,
            x162 * x21 * x95,
            x102 * x165 * x169,
            x102 * x173 * x177,
            x115 * x165 * x177,
            x102 * x183 * x184,
            x115 * x173 * x184,
            x106 * x165 * x184,
            x102 * x168 * x186,
            x102 * x176 * x188,
            x115 * x176 * x186,
            x102 * x175 * x197,
            x115 * x175 * x188,
            x106 * x175 * x186,
            x118 * x165 * x168,
            x118 * x173 * x176,
            x120 * x165 * x176,
            x118 * x175 * x183,
            x120 * x173 * x175,
            x122 * x165 * x175,
            x102 * x201 * x75,
            x102 * x206 * x90,
            x115 * x201 * x90,
            x102 * x174 * x213,
            x115 * x174 * x206,
            x106 * x174 * x201,
            x118 * x186 * x75,
            x118 * x188 * x90,
            x120 * x186 * x90,
            x118 * x174 * x197,
            x120 * x174 * x188,
            x122 * x174 * x186,
            x135 * x165 * x75,
            x135 * x173 * x90,
            x138 * x165 * x90,
            x135 * x174 * x183,
            x138 * x173 * x174,
            x144 * x165 * x174,
            x102 * x217 * x55,
            x102 * x223 * x88,
            x115 * x217 * x88,
            x2 * x228,
            x223 * x229 * x94,
            x106 * x217 * x83,
            x118 * x201 * x55,
            x118 * x206 * x88,
            x120 * x201 * x88,
            x116 * x213 * x229,
            x120 * x206 * x83,
            x122 * x201 * x83,
            x135 * x186 * x55,
            x135 * x188 * x88,
            x138 * x186 * x88,
            x135 * x197 * x83,
            x138 * x188 * x83,
            x144 * x186 * x83,
            x151 * x165 * x55,
            x151 * x173 * x88,
            x154 * x165 * x88,
            x151 * x183 * x83,
            x154 * x173 * x83,
            x156 * x165 * x83,
            x102 * x231 * x53,
            x232 * x5,
            x231 * x233 * x94,
            x227
            * (
                x107 * x226
                + x3
                * (
                    x129 * (x224 + x225)
                    + x18 * (x133 + x211 + x221)
                    + 3 * x207
                    + 3 * x212
                    + 2 * x218
                    + 2 * x222
                )
            ),
            x232 * x94,
            x106 * x231 * x9,
            x118 * x217 * x53,
            x116 * x223 * x233,
            x10 * x120 * x217,
            x116 * x228,
            x120 * x223 * x9,
            x122 * x217 * x9,
            x135 * x201 * x53,
            x10 * x135 * x206,
            x10 * x138 * x201,
            x135 * x213 * x9,
            x138 * x206 * x9,
            x144 * x201 * x9,
            x151 * x186 * x53,
            x10 * x151 * x188,
            x10 * x154 * x186,
            x151 * x197 * x9,
            x154 * x188 * x9,
            x156 * x186 * x9,
            x160 * x165 * x53,
            x10 * x160 * x173,
            x10 * x161 * x165,
            x160 * x183 * x9,
            x161 * x173 * x9,
            x162 * x165 * x9,
            x169 * x236 * x95,
            x112 * x177 * x236,
            x177 * x240 * x95,
            x184 * x236 * x98,
            x112 * x184 * x240,
            x184 * x246 * x95,
            x109 * x168 * x236,
            x111 * x176 * x236,
            x109 * x176 * x240,
            x114 * x175 * x236,
            x111 * x175 * x240,
            x109 * x175 * x246,
            x168 * x248 * x95,
            x112 * x176 * x248,
            x176 * x250 * x95,
            x175 * x248 * x98,
            x112 * x175 * x250,
            x175 * x259 * x95,
            x124 * x236 * x75,
            x127 * x236 * x90,
            x124 * x240 * x90,
            x133 * x174 * x236,
            x127 * x174 * x240,
            x124 * x174 * x246,
            x109 * x248 * x75,
            x111 * x248 * x90,
            x109 * x250 * x90,
            x114 * x174 * x248,
            x111 * x174 * x250,
            x109 * x174 * x259,
            x263 * x75 * x95,
            x112 * x263 * x90,
            x268 * x90 * x95,
            x174 * x263 * x98,
            x112 * x174 * x268,
            x174 * x275 * x95,
            x145 * x236 * x55,
            x148 * x236 * x88,
            x145 * x240 * x88,
            x150 * x236 * x83,
            x148 * x240 * x83,
            x145 * x246 * x83,
            x124 * x248 * x55,
            x127 * x248 * x88,
            x124 * x250 * x88,
            x133 * x248 * x83,
            x127 * x250 * x83,
            x124 * x259 * x83,
            x109 * x263 * x55,
            x111 * x263 * x88,
            x109 * x268 * x88,
            x114 * x263 * x83,
            x111 * x268 * x83,
            x107 * x275 * x277,
            x281 * x55 * x95,
            x112 * x281 * x88,
            x287 * x88 * x95,
            x281 * x83 * x98,
            x277 * x287 * x81,
            x2 * x291,
            x157 * x236 * x53,
            x10 * x158 * x236,
            x10 * x157 * x240,
            x159 * x236 * x9,
            x158 * x240 * x9,
            x157 * x246 * x9,
            x145 * x248 * x53,
            x10 * x148 * x248,
            x10 * x145 * x250,
            x150 * x248 * x9,
            x148 * x250 * x9,
            x145 * x259 * x9,
            x124 * x263 * x53,
            x10 * x127 * x263,
            x10 * x124 * x268,
            x133 * x263 * x9,
            x127 * x268 * x9,
            x124 * x275 * x9,
            x109 * x281 * x53,
            x10 * x111 * x281,
            x107 * x287 * x292,
            x114 * x281 * x9,
            x111 * x287 * x9,
            x107 * x291,
            x294 * x53 * x95,
            x292 * x294 * x81,
            x295 * x5,
            x294 * x9 * x98,
            x295 * x81,
            x276
            * (
                x116 * x290
                + x3
                * (
                    x140 * (x288 + x289)
                    + x18 * (x144 + x273 + x285)
                    + 3 * x269
                    + 3 * x274
                    + 2 * x282
                    + 2 * x286
                )
            ),
        ]
    )


def diag_quadrupole3d_43(a, A, b, B, C):
    """Cartesian 3D (gf) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x37 = 4 * x10 * x13
    x38 = x36 + x37
    x39 = x3 * (2 * x16 + 2 * x23 + x38)
    x40 = x35 + x39
    x41 = x2 * x40
    x42 = x33 + x41
    x43 = x2 * x42
    x44 = x2 * x34
    x45 = 3 * x12
    x46 = x19 * x9
    x47 = x13 + x26
    x48 = x4 * x47
    x49 = x46 + x48
    x50 = x3 * (3 * x16 + x45 + x49)
    x51 = x2 * (x18 + x29)
    x52 = 2 * x50 + 2 * x51
    x53 = x3 * (x35 + 4 * x39 + 3 * x44 + x52)
    x54 = x43 + x53
    x55 = x17 * x2
    x56 = x3 * (3 * x26 + x27)
    x57 = x2 * x49
    x58 = x56 + x57
    x59 = x3 * (x18 + 4 * x29 + 3 * x55 + x58)
    x60 = x2 * (x50 + x51)
    x61 = 2 * x55
    x62 = x2 * x24
    x63 = 2 * x62
    x64 = x3 * (x25 + x32 + x61 + x63)
    x65 = x39 + x44
    x66 = x2 * x65
    x67 = 3 * x64 + 3 * x66
    x68 = x2 * x54 + x3 * (2 * x33 + 2 * x41 + 2 * x59 + 2 * x60 + x67)
    x69 = x64 + x66
    x70 = x2 * x69
    x71 = 2 * x2
    x72 = x15 * x2
    x73 = 2 * x72
    x74 = x2 * x47
    x75 = x46 + x74
    x76 = x3 * (x16 + x45 + x73 + x75)
    x77 = x2 * (x29 + x55)
    x78 = x13 * x4
    x79 = x3 * (x48 + 3 * x74 + 8 * x78)
    x80 = x2 * x58
    x81 = x79 + x80
    x82 = x31 + x62
    x83 = x2 * x82
    x84 = x2 * x22
    x85 = x3 * (x23 + x38 + x73 + x84)
    x86 = x3 * (2 * x39 + 2 * x44 + 2 * x76 + 2 * x77 + 2 * x83 + 2 * x85)
    x87 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x88 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x89 = numpy.pi * x0 * x88
    x90 = x87 * x89
    x91 = -x0 * (a * A[1] + b * B[1])
    x92 = -x91 - B[1]
    x93 = x70 + x86
    x94 = x27 + x71 * x9
    x95 = x3 * (x26 + x94)
    x96 = x2 * x75
    x97 = x95 + x96
    x98 = x2 * x9
    x99 = x11 * x2
    x100 = x3 * (x14 + x27 + x98 + x99)
    x101 = x2 * (x12 + x72)
    x102 = 2 * x100 + 2 * x101
    x103 = x11 * x71 + x27
    x104 = x3 * (x103 + x21)
    x105 = x20 + x84
    x106 = x105 * x2
    x107 = x104 + x106
    x108 = x3 * (x102 + x107 + 2 * x31 + x63)
    x109 = x83 + x85
    x110 = x109 * x2
    x111 = x90 * (
        x2 * x93
        + x3
        * (2 * x108 + 2 * x110 + x19 * (x102 + x30 + x61 + x97) + x67 + x71 * (x76 + x77))
    )
    x112 = -x0 * (a * A[2] + b * B[2])
    x113 = -x112 - B[2]
    x114 = x7 * x87
    x115 = x114 * x3
    x116 = x114 * x92 ** 2
    x117 = x115 + x116
    x118 = x108 + x110
    x119 = x2 * x8
    x120 = x3 * (x11 + x119)
    x121 = x2 * (x13 + x99)
    x122 = x120 + x121
    x123 = x3 * (x119 + x9)
    x124 = x13 + x98
    x125 = x124 * x2
    x126 = x123 + x125
    x127 = x107 * x2 + x3 * (2 * x120 + 2 * x121 + x37 + 2 * x84)
    x128 = x118 * x2 + x3 * (
        x127 + x19 * (x122 + x126 + x36 + x73) + x71 * (x100 + x101) + 3 * x83 + 3 * x85
    )
    x129 = x7 * x88
    x130 = x113 * x90
    x131 = x129 * x3
    x132 = x113 ** 2 * x129
    x133 = x131 + x132
    x134 = x114 * x92
    x135 = x134 * x19
    x136 = x117 * x92
    x137 = x135 + x136
    x138 = x2 ** 2 * x8
    x139 = x127 * x2 + x3 * (3 * x104 + 3 * x106 + x122 * x71 + x19 * (x103 + x138))
    x140 = x113 * x129
    x141 = x140 * x19
    x142 = x113 * x133
    x143 = x141 + x142
    x144 = -x91 - A[1]
    x145 = x68 * x90
    x146 = x114 * x144
    x147 = x146 * x92
    x148 = x115 + x147
    x149 = x117 * x144
    x150 = x135 + x149
    x151 = 3 * x115
    x152 = x3 * (3 * x116 + x151)
    x153 = x137 * x144
    x154 = x152 + x153
    x155 = -x112 - A[2]
    x156 = x155 * x90
    x157 = x129 * x155
    x158 = x113 * x157
    x159 = x131 + x158
    x160 = x133 * x155
    x161 = x141 + x160
    x162 = 3 * x131
    x163 = x3 * (3 * x132 + x162)
    x164 = x143 * x155
    x165 = x163 + x164
    x166 = x114 * x144 ** 2
    x167 = x115 + x166
    x168 = x3 * (x134 + x146)
    x169 = x144 * x148
    x170 = x168 + x169
    x171 = 2 * x144
    x172 = x134 * x171 + x151
    x173 = x3 * (x116 + x172)
    x174 = x144 * x150
    x175 = x173 + x174
    x176 = x115 * x92
    x177 = x3 * (x136 + 3 * x149 + 8 * x176)
    x178 = x144 * x154
    x179 = x177 + x178
    x180 = x129 * x155 ** 2
    x181 = x131 + x180
    x182 = x3 * (x140 + x157)
    x183 = x155 * x159
    x184 = x182 + x183
    x185 = 2 * x155
    x186 = x140 * x185 + x162
    x187 = x3 * (x132 + x186)
    x188 = x155 * x161
    x189 = x187 + x188
    x190 = x113 * x131
    x191 = x3 * (x142 + 3 * x160 + 8 * x190)
    x192 = x155 * x165
    x193 = x191 + x192
    x194 = x144 * x167 + x146 * x19
    x195 = x3 * (x166 + x172)
    x196 = x144 * x170
    x197 = x195 + x196
    x198 = x144 * x175
    x199 = x3 * (2 * x149 + 2 * x168 + 2 * x169 + 4 * x176)
    x200 = x198 + x199
    x201 = 3 * x173 + 3 * x174
    x202 = x144 * x179 + x3 * (2 * x152 + 2 * x153 + x201)
    x203 = x155 * x181 + x157 * x19
    x204 = x3 * (x180 + x186)
    x205 = x155 * x184
    x206 = x204 + x205
    x207 = x155 * x189
    x208 = x3 * (2 * x160 + 2 * x182 + 2 * x183 + 4 * x190)
    x209 = x207 + x208
    x210 = 3 * x187 + 3 * x188
    x211 = x155 * x193 + x3 * (2 * x163 + 2 * x164 + x210)
    x212 = x144 * x194 + x3 * (x151 + 3 * x166)
    x213 = x144 * x197 + x3 * (3 * x168 + 3 * x169 + x194)
    x214 = x144 * x200 + x3 * (2 * x195 + 2 * x196 + x201)
    x215 = x144 * x202 + x3 * (3 * x177 + 3 * x178 + 3 * x198 + 3 * x199)
    x216 = x155 * x203 + x3 * (x162 + 3 * x180)
    x217 = x155 * x206 + x3 * (3 * x182 + 3 * x183 + x203)
    x218 = x155 * x209 + x3 * (2 * x204 + 2 * x205 + x210)
    x219 = x155 * x211 + x3 * (3 * x191 + 3 * x192 + 3 * x207 + 3 * x208)
    x220 = -x91 - C[1]
    x221 = x114 * x220 ** 2
    x222 = x115 + x221
    x223 = 3 * x95 + 3 * x96
    x224 = x2 * x81 + x3 * (x223 + 2 * x56 + 2 * x57)
    x225 = x2 * x97
    x226 = x3 * (2 * x123 + 2 * x125 + 2 * x74 + 4 * x78)
    x227 = x2 * x224 + x3 * (3 * x225 + 3 * x226 + 3 * x79 + 3 * x80)
    x228 = x114 * x220
    x229 = x19 * x228
    x230 = x222 * x92
    x231 = x229 + x230
    x232 = x3 * (x138 + x94)
    x233 = x126 * x2
    x234 = x225 + x226
    x235 = x2 * x234 + x3 * (x223 + 2 * x232 + 2 * x233)
    x236 = x134 * x220
    x237 = x151 + 2 * x236
    x238 = x3 * (x221 + x237)
    x239 = x231 * x92
    x240 = x238 + x239
    x241 = x13 + x138
    x242 = x119 * x19 + x2 * x241
    x243 = x232 + x233
    x244 = x2 * x243 + x3 * (3 * x123 + 3 * x125 + x242)
    x245 = x240 * x92
    x246 = x115 + x236
    x247 = x246 * x92
    x248 = x3 * (x134 + x228)
    x249 = 2 * x248
    x250 = 4 * x115 * x220
    x251 = x249 + x250
    x252 = x3 * (2 * x230 + 2 * x247 + x251)
    x253 = x245 + x252
    x254 = x2 * x242 + x3 * (3 * x138 + x27)
    x255 = x144 * x222
    x256 = x229 + x255
    x257 = x144 * x231
    x258 = x238 + x257
    x259 = x144 * x240
    x260 = x252 + x259
    x261 = x247 + x248
    x262 = x261 * x92
    x263 = x3 * (x116 + x237)
    x264 = 2 * x263
    x265 = 3 * x238 + x264
    x266 = x3 * (3 * x239 + 2 * x262 + x265)
    x267 = x144 * x253
    x268 = x266 + x267
    x269 = x151 + x171 * x228
    x270 = x3 * (x221 + x269)
    x271 = x144 * x256
    x272 = x270 + x271
    x273 = x144 * x258
    x274 = x144 * x246
    x275 = 2 * x274
    x276 = x3 * (x230 + x251 + x255 + x275)
    x277 = x273 + x276
    x278 = x144 * x260
    x279 = x144 * x261
    x280 = 2 * x279
    x281 = 2 * x257
    x282 = x3 * (x239 + x265 + x280 + x281)
    x283 = x278 + x282
    x284 = x144 * x268
    x285 = 3 * x248
    x286 = x3 * (x137 + 3 * x247 + x285)
    x287 = x144 * (x262 + x263)
    x288 = 2 * x286 + 2 * x287
    x289 = x3 * (x245 + 4 * x252 + 3 * x259 + x288)
    x290 = x284 + x289
    x291 = x3 * (x146 + x228)
    x292 = x146 * x220
    x293 = x144 * (x115 + x292)
    x294 = x144 * x272 + x3 * (x250 + 2 * x255 + 2 * x291 + 2 * x293)
    x295 = x144 * x277
    x296 = x3 * (x147 + x151 + x236 + x292)
    x297 = x144 * (x248 + x274)
    x298 = 2 * x296 + 2 * x297
    x299 = x3 * (2 * x238 + x272 + x281 + x298)
    x300 = x295 + x299
    x301 = x144 * x283
    x302 = x144 * (x263 + x279)
    x303 = x3 * (x150 + x247 + x275 + x285)
    x304 = x3 * (2 * x252 + 2 * x259 + 2 * x273 + 2 * x276 + 2 * x302 + 2 * x303)
    x305 = x301 + x304
    x306 = x3 * (x154 + x262 + 4 * x263 + 3 * x279)
    x307 = x144 * (x286 + x287)
    x308 = 3 * x278 + 3 * x282
    x309 = x144 * x290 + x3 * (2 * x266 + 2 * x267 + 2 * x306 + 2 * x307 + x308)
    x310 = x6 * x89
    x311 = x309 * x310
    x312 = x2 * x310
    x313 = x291 + x293
    x314 = x144 * x294 + x3 * (x171 * x313 + x19 * (x166 + x269) + 3 * x270 + 3 * x271)
    x315 = x144 * x300 + x3 * (
        x171 * (x296 + x297)
        + x19 * (x170 + x249 + x275 + x313)
        + 3 * x273
        + 3 * x276
        + x294
    )
    x316 = x310 * (
        x144 * x305
        + x3
        * (
            x171 * (x302 + x303)
            + x19 * (x175 + x264 + x280 + x298)
            + 2 * x295
            + 2 * x299
            + x308
        )
    )
    x317 = x310 * x4
    x318 = -x112 - C[2]
    x319 = x129 * x318 ** 2
    x320 = x131 + x319
    x321 = x129 * x318
    x322 = x19 * x321
    x323 = x113 * x320
    x324 = x322 + x323
    x325 = x140 * x318
    x326 = x162 + 2 * x325
    x327 = x3 * (x319 + x326)
    x328 = x113 * x324
    x329 = x327 + x328
    x330 = x113 * x329
    x331 = x131 + x325
    x332 = x113 * x331
    x333 = x3 * (x140 + x321)
    x334 = 2 * x333
    x335 = 4 * x131 * x318
    x336 = x334 + x335
    x337 = x3 * (2 * x323 + 2 * x332 + x336)
    x338 = x330 + x337
    x339 = x155 * x320
    x340 = x322 + x339
    x341 = x155 * x324
    x342 = x327 + x341
    x343 = x155 * x329
    x344 = x337 + x343
    x345 = x332 + x333
    x346 = x113 * x345
    x347 = x3 * (x132 + x326)
    x348 = 2 * x347
    x349 = 3 * x327 + x348
    x350 = x3 * (3 * x328 + 2 * x346 + x349)
    x351 = x155 * x338
    x352 = x350 + x351
    x353 = x162 + x185 * x321
    x354 = x3 * (x319 + x353)
    x355 = x155 * x340
    x356 = x354 + x355
    x357 = x155 * x342
    x358 = x155 * x331
    x359 = 2 * x358
    x360 = x3 * (x323 + x336 + x339 + x359)
    x361 = x357 + x360
    x362 = x155 * x344
    x363 = x155 * x345
    x364 = 2 * x363
    x365 = 2 * x341
    x366 = x3 * (x328 + x349 + x364 + x365)
    x367 = x362 + x366
    x368 = x155 * x352
    x369 = 3 * x333
    x370 = x3 * (x143 + 3 * x332 + x369)
    x371 = x155 * (x346 + x347)
    x372 = 2 * x370 + 2 * x371
    x373 = x3 * (x330 + 4 * x337 + 3 * x343 + x372)
    x374 = x368 + x373
    x375 = numpy.pi * x0 * x6 * x87
    x376 = x2 * x375
    x377 = x3 * (x157 + x321)
    x378 = x157 * x318
    x379 = x155 * (x131 + x378)
    x380 = x155 * x356 + x3 * (x335 + 2 * x339 + 2 * x377 + 2 * x379)
    x381 = x155 * x361
    x382 = x3 * (x158 + x162 + x325 + x378)
    x383 = x155 * (x333 + x358)
    x384 = 2 * x382 + 2 * x383
    x385 = x3 * (2 * x327 + x356 + x365 + x384)
    x386 = x381 + x385
    x387 = x155 * x367
    x388 = x155 * (x347 + x363)
    x389 = x3 * (x161 + x332 + x359 + x369)
    x390 = x3 * (2 * x337 + 2 * x343 + 2 * x357 + 2 * x360 + 2 * x388 + 2 * x389)
    x391 = x387 + x390
    x392 = x3 * (x165 + x346 + 4 * x347 + 3 * x363)
    x393 = x155 * (x370 + x371)
    x394 = 3 * x362 + 3 * x366
    x395 = x155 * x374 + x3 * (2 * x350 + 2 * x351 + 2 * x392 + 2 * x393 + x394)
    x396 = x375 * x395
    x397 = x375 * x4
    x398 = x377 + x379
    x399 = x155 * x380 + x3 * (x185 * x398 + x19 * (x180 + x353) + 3 * x354 + 3 * x355)
    x400 = x155 * x386 + x3 * (
        x185 * (x382 + x383)
        + x19 * (x184 + x334 + x359 + x398)
        + 3 * x357
        + 3 * x360
        + x380
    )
    x401 = x375 * (
        x155 * x391
        + x3
        * (
            x185 * (x388 + x389)
            + x19 * (x189 + x348 + x364 + x384)
            + 2 * x381
            + 2 * x385
            + x394
        )
    )

    # 450 item(s)
    return numpy.array(
        [
            x90
            * (
                x2 * x68
                + x3
                * (
                    x19 * (x52 + 3 * x76 + 3 * x77 + x81)
                    + 3 * x43
                    + 3 * x53
                    + 3 * x70
                    + x71 * (x59 + x60)
                    + 3 * x86
                )
            ),
            x111 * x92,
            x111 * x113,
            x117 * x128 * x129,
            x128 * x130 * x92,
            x114 * x128 * x133,
            x129 * x137 * x139,
            x117 * x139 * x140,
            x133 * x134 * x139,
            x114 * x139 * x143,
            x144 * x145,
            x129 * x148 * x93,
            x130 * x144 * x93,
            x118 * x129 * x150,
            x118 * x140 * x148,
            x118 * x133 * x146,
            x127 * x129 * x154,
            x127 * x140 * x150,
            x127 * x133 * x148,
            x127 * x143 * x146,
            x145 * x155,
            x156 * x92 * x93,
            x114 * x159 * x93,
            x117 * x118 * x157,
            x118 * x134 * x159,
            x114 * x118 * x161,
            x127 * x137 * x157,
            x117 * x127 * x159,
            x127 * x134 * x161,
            x114 * x127 * x165,
            x129 * x167 * x54,
            x129 * x170 * x69,
            x140 * x167 * x69,
            x109 * x129 * x175,
            x109 * x140 * x170,
            x109 * x133 * x167,
            x107 * x129 * x179,
            x107 * x140 * x175,
            x107 * x133 * x170,
            x107 * x143 * x167,
            x144 * x156 * x54,
            x148 * x157 * x69,
            x146 * x159 * x69,
            x109 * x150 * x157,
            x109 * x148 * x159,
            x109 * x146 * x161,
            x107 * x154 * x157,
            x107 * x150 * x159,
            x107 * x148 * x161,
            x107 * x146 * x165,
            x114 * x181 * x54,
            x134 * x181 * x69,
            x114 * x184 * x69,
            x109 * x117 * x181,
            x109 * x134 * x184,
            x109 * x114 * x189,
            x107 * x137 * x181,
            x107 * x117 * x184,
            x107 * x134 * x189,
            x107 * x114 * x193,
            x129 * x194 * x42,
            x129 * x197 * x65,
            x140 * x194 * x65,
            x129 * x200 * x82,
            x140 * x197 * x82,
            x133 * x194 * x82,
            x105 * x129 * x202,
            x105 * x140 * x200,
            x105 * x133 * x197,
            x105 * x143 * x194,
            x157 * x167 * x42,
            x157 * x170 * x65,
            x159 * x167 * x65,
            x157 * x175 * x82,
            x159 * x170 * x82,
            x161 * x167 * x82,
            x105 * x157 * x179,
            x105 * x159 * x175,
            x105 * x161 * x170,
            x105 * x165 * x167,
            x146 * x181 * x42,
            x148 * x181 * x65,
            x146 * x184 * x65,
            x150 * x181 * x82,
            x148 * x184 * x82,
            x146 * x189 * x82,
            x105 * x154 * x181,
            x105 * x150 * x184,
            x105 * x148 * x189,
            x105 * x146 * x193,
            x114 * x203 * x42,
            x134 * x203 * x65,
            x114 * x206 * x65,
            x117 * x203 * x82,
            x134 * x206 * x82,
            x114 * x209 * x82,
            x105 * x137 * x203,
            x105 * x117 * x206,
            x105 * x134 * x209,
            x105 * x114 * x211,
            x129 * x212 * x40,
            x129 * x213 * x34,
            x140 * x212 * x34,
            x129 * x214 * x24,
            x140 * x213 * x24,
            x133 * x212 * x24,
            x129 * x215 * x22,
            x140 * x214 * x22,
            x133 * x213 * x22,
            x143 * x212 * x22,
            x157 * x194 * x40,
            x157 * x197 * x34,
            x159 * x194 * x34,
            x157 * x200 * x24,
            x159 * x197 * x24,
            x161 * x194 * x24,
            x157 * x202 * x22,
            x159 * x200 * x22,
            x161 * x197 * x22,
            x165 * x194 * x22,
            x167 * x181 * x40,
            x170 * x181 * x34,
            x167 * x184 * x34,
            x175 * x181 * x24,
            x170 * x184 * x24,
            x167 * x189 * x24,
            x179 * x181 * x22,
            x175 * x184 * x22,
            x170 * x189 * x22,
            x167 * x193 * x22,
            x146 * x203 * x40,
            x148 * x203 * x34,
            x146 * x206 * x34,
            x150 * x203 * x24,
            x148 * x206 * x24,
            x146 * x209 * x24,
            x154 * x203 * x22,
            x150 * x206 * x22,
            x148 * x209 * x22,
            x146 * x211 * x22,
            x114 * x216 * x40,
            x134 * x216 * x34,
            x114 * x217 * x34,
            x117 * x216 * x24,
            x134 * x217 * x24,
            x114 * x218 * x24,
            x137 * x216 * x22,
            x117 * x217 * x22,
            x134 * x218 * x22,
            x114 * x219 * x22,
            x129 * x222 * x227,
            x129 * x231 * x235,
            x140 * x222 * x235,
            x129 * x240 * x244,
            x140 * x231 * x244,
            x133 * x222 * x244,
            x129 * x253 * x254,
            x140 * x240 * x254,
            x133 * x231 * x254,
            x143 * x222 * x254,
            x129 * x224 * x256,
            x129 * x234 * x258,
            x140 * x234 * x256,
            x129 * x243 * x260,
            x140 * x243 * x258,
            x133 * x243 * x256,
            x129 * x242 * x268,
            x140 * x242 * x260,
            x133 * x242 * x258,
            x143 * x242 * x256,
            x157 * x222 * x224,
            x157 * x231 * x234,
            x159 * x222 * x234,
            x157 * x240 * x243,
            x159 * x231 * x243,
            x161 * x222 * x243,
            x157 * x242 * x253,
            x159 * x240 * x242,
            x161 * x231 * x242,
            x165 * x222 * x242,
            x129 * x272 * x81,
            x129 * x277 * x97,
            x140 * x272 * x97,
            x126 * x129 * x283,
            x126 * x140 * x277,
            x126 * x133 * x272,
            x129 * x241 * x290,
            x140 * x241 * x283,
            x133 * x241 * x277,
            x143 * x241 * x272,
            x157 * x256 * x81,
            x157 * x258 * x97,
            x159 * x256 * x97,
            x126 * x157 * x260,
            x126 * x159 * x258,
            x126 * x161 * x256,
            x157 * x241 * x268,
            x159 * x241 * x260,
            x161 * x241 * x258,
            x165 * x241 * x256,
            x181 * x222 * x81,
            x181 * x231 * x97,
            x184 * x222 * x97,
            x126 * x181 * x240,
            x126 * x184 * x231,
            x126 * x189 * x222,
            x181 * x241 * x253,
            x184 * x240 * x241,
            x189 * x231 * x241,
            x193 * x222 * x241,
            x129 * x294 * x58,
            x129 * x300 * x75,
            x140 * x294 * x75,
            x124 * x129 * x305,
            x124 * x140 * x300,
            x124 * x133 * x294,
            x2 * x311,
            x113 * x305 * x312,
            x119 * x133 * x300,
            x119 * x143 * x294,
            x157 * x272 * x58,
            x157 * x277 * x75,
            x159 * x272 * x75,
            x124 * x157 * x283,
            x124 * x159 * x277,
            x124 * x161 * x272,
            x155 * x290 * x312,
            x119 * x159 * x283,
            x119 * x161 * x277,
            x119 * x165 * x272,
            x181 * x256 * x58,
            x181 * x258 * x75,
            x184 * x256 * x75,
            x124 * x181 * x260,
            x124 * x184 * x258,
            x124 * x189 * x256,
            x119 * x181 * x268,
            x119 * x184 * x260,
            x119 * x189 * x258,
            x119 * x193 * x256,
            x203 * x222 * x58,
            x203 * x231 * x75,
            x206 * x222 * x75,
            x124 * x203 * x240,
            x124 * x206 * x231,
            x124 * x209 * x222,
            x119 * x203 * x253,
            x119 * x206 * x240,
            x119 * x209 * x231,
            x119 * x211 * x222,
            x129 * x314 * x49,
            x129 * x315 * x47,
            x140 * x314 * x47,
            x316 * x4,
            x113 * x315 * x317,
            x133 * x314 * x9,
            x310
            * (
                x144 * x309
                + x3
                * (
                    x171 * (x306 + x307)
                    + x19 * (x179 + x288 + 3 * x302 + 3 * x303)
                    + 3 * x284
                    + 3 * x289
                    + 3 * x301
                    + 3 * x304
                )
            ),
            x113 * x316,
            x133 * x315 * x8,
            x143 * x314 * x8,
            x157 * x294 * x49,
            x157 * x300 * x47,
            x159 * x294 * x47,
            x155 * x305 * x317,
            x159 * x300 * x9,
            x161 * x294 * x9,
            x155 * x311,
            x159 * x305 * x8,
            x161 * x300 * x8,
            x165 * x294 * x8,
            x181 * x272 * x49,
            x181 * x277 * x47,
            x184 * x272 * x47,
            x181 * x283 * x9,
            x184 * x277 * x9,
            x189 * x272 * x9,
            x181 * x290 * x8,
            x184 * x283 * x8,
            x189 * x277 * x8,
            x193 * x272 * x8,
            x203 * x256 * x49,
            x203 * x258 * x47,
            x206 * x256 * x47,
            x203 * x260 * x9,
            x206 * x258 * x9,
            x209 * x256 * x9,
            x203 * x268 * x8,
            x206 * x260 * x8,
            x209 * x258 * x8,
            x211 * x256 * x8,
            x216 * x222 * x49,
            x216 * x231 * x47,
            x217 * x222 * x47,
            x216 * x240 * x9,
            x217 * x231 * x9,
            x218 * x222 * x9,
            x216 * x253 * x8,
            x217 * x240 * x8,
            x218 * x231 * x8,
            x219 * x222 * x8,
            x114 * x227 * x320,
            x134 * x235 * x320,
            x114 * x235 * x324,
            x117 * x244 * x320,
            x134 * x244 * x324,
            x114 * x244 * x329,
            x137 * x254 * x320,
            x117 * x254 * x324,
            x134 * x254 * x329,
            x114 * x254 * x338,
            x146 * x224 * x320,
            x148 * x234 * x320,
            x146 * x234 * x324,
            x150 * x243 * x320,
            x148 * x243 * x324,
            x146 * x243 * x329,
            x154 * x242 * x320,
            x150 * x242 * x324,
            x148 * x242 * x329,
            x146 * x242 * x338,
            x114 * x224 * x340,
            x134 * x234 * x340,
            x114 * x234 * x342,
            x117 * x243 * x340,
            x134 * x243 * x342,
            x114 * x243 * x344,
            x137 * x242 * x340,
            x117 * x242 * x342,
            x134 * x242 * x344,
            x114 * x242 * x352,
            x167 * x320 * x81,
            x170 * x320 * x97,
            x167 * x324 * x97,
            x126 * x175 * x320,
            x126 * x170 * x324,
            x126 * x167 * x329,
            x179 * x241 * x320,
            x175 * x241 * x324,
            x170 * x241 * x329,
            x167 * x241 * x338,
            x146 * x340 * x81,
            x148 * x340 * x97,
            x146 * x342 * x97,
            x126 * x150 * x340,
            x126 * x148 * x342,
            x126 * x146 * x344,
            x154 * x241 * x340,
            x150 * x241 * x342,
            x148 * x241 * x344,
            x146 * x241 * x352,
            x114 * x356 * x81,
            x134 * x356 * x97,
            x114 * x361 * x97,
            x117 * x126 * x356,
            x126 * x134 * x361,
            x114 * x126 * x367,
            x137 * x241 * x356,
            x117 * x241 * x361,
            x134 * x241 * x367,
            x114 * x241 * x374,
            x194 * x320 * x58,
            x197 * x320 * x75,
            x194 * x324 * x75,
            x124 * x200 * x320,
            x124 * x197 * x324,
            x124 * x194 * x329,
            x119 * x202 * x320,
            x119 * x200 * x324,
            x119 * x197 * x329,
            x119 * x194 * x338,
            x167 * x340 * x58,
            x170 * x340 * x75,
            x167 * x342 * x75,
            x124 * x175 * x340,
            x124 * x170 * x342,
            x124 * x167 * x344,
            x119 * x179 * x340,
            x119 * x175 * x342,
            x119 * x170 * x344,
            x119 * x167 * x352,
            x146 * x356 * x58,
            x148 * x356 * x75,
            x146 * x361 * x75,
            x124 * x150 * x356,
            x124 * x148 * x361,
            x124 * x146 * x367,
            x119 * x154 * x356,
            x119 * x150 * x361,
            x119 * x148 * x367,
            x144 * x374 * x376,
            x114 * x380 * x58,
            x134 * x380 * x75,
            x114 * x386 * x75,
            x117 * x124 * x380,
            x124 * x134 * x386,
            x114 * x124 * x391,
            x119 * x137 * x380,
            x117 * x119 * x386,
            x376 * x391 * x92,
            x2 * x396,
            x212 * x320 * x49,
            x213 * x320 * x47,
            x212 * x324 * x47,
            x214 * x320 * x9,
            x213 * x324 * x9,
            x212 * x329 * x9,
            x215 * x320 * x8,
            x214 * x324 * x8,
            x213 * x329 * x8,
            x212 * x338 * x8,
            x194 * x340 * x49,
            x197 * x340 * x47,
            x194 * x342 * x47,
            x200 * x340 * x9,
            x197 * x342 * x9,
            x194 * x344 * x9,
            x202 * x340 * x8,
            x200 * x342 * x8,
            x197 * x344 * x8,
            x194 * x352 * x8,
            x167 * x356 * x49,
            x170 * x356 * x47,
            x167 * x361 * x47,
            x175 * x356 * x9,
            x170 * x361 * x9,
            x167 * x367 * x9,
            x179 * x356 * x8,
            x175 * x361 * x8,
            x170 * x367 * x8,
            x167 * x374 * x8,
            x146 * x380 * x49,
            x148 * x380 * x47,
            x146 * x386 * x47,
            x150 * x380 * x9,
            x148 * x386 * x9,
            x144 * x391 * x397,
            x154 * x380 * x8,
            x150 * x386 * x8,
            x148 * x391 * x8,
            x144 * x396,
            x114 * x399 * x49,
            x134 * x399 * x47,
            x114 * x400 * x47,
            x117 * x399 * x9,
            x397 * x400 * x92,
            x4 * x401,
            x137 * x399 * x8,
            x117 * x400 * x8,
            x401 * x92,
            x375
            * (
                x155 * x395
                + x3
                * (
                    x185 * (x392 + x393)
                    + x19 * (x193 + x372 + 3 * x388 + 3 * x389)
                    + 3 * x368
                    + 3 * x373
                    + 3 * x387
                    + 3 * x390
                )
            ),
        ]
    )


def diag_quadrupole3d_44(a, A, b, B, C):
    """Cartesian 3D (gg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

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
    x37 = 4 * x10 * x13
    x38 = x36 + x37
    x39 = x3 * (2 * x16 + 2 * x23 + x38)
    x40 = x35 + x39
    x41 = x4 * x40
    x42 = x33 + x41
    x43 = x2 * x42
    x44 = x18 + x29
    x45 = x4 * x44
    x46 = 3 * x12
    x47 = x19 * x9
    x48 = x13 + x26
    x49 = x4 * x48
    x50 = x47 + x49
    x51 = x3 * (3 * x16 + x46 + x50)
    x52 = 2 * x51
    x53 = 4 * x39 + x52
    x54 = x3 * (4 * x35 + 2 * x45 + x53)
    x55 = x43 + x54
    x56 = x2 * x55
    x57 = x2 * x40
    x58 = 4 * x29
    x59 = x3 * (3 * x26 + x27)
    x60 = x4 * x50
    x61 = x59 + x60
    x62 = x3 * (4 * x18 + x58 + x61)
    x63 = x2 * (x45 + x51)
    x64 = 2 * x62 + 2 * x63
    x65 = x3 * (5 * x33 + x41 + 4 * x57 + x64)
    x66 = x56 + x65
    x67 = x2 * (x62 + x63)
    x68 = x33 + x57
    x69 = x2 * x68
    x70 = x2 * x44
    x71 = x13 * x4
    x72 = 8 * x71
    x73 = x3 * (4 * x49 + x72)
    x74 = x2 * x61
    x75 = x73 + x74
    x76 = x3 * (x45 + 5 * x51 + 4 * x70 + x75)
    x77 = 2 * x70
    x78 = x2 * x34
    x79 = x3 * (x35 + x53 + x77 + 3 * x78)
    x80 = x2 * x66 + x3 * (2 * x43 + 2 * x54 + 2 * x67 + 4 * x69 + 2 * x76 + 4 * x79)
    x81 = 2 * x2
    x82 = x69 + x79
    x83 = x2 * x82
    x84 = x17 * x2
    x85 = x2 * x50
    x86 = x59 + x85
    x87 = x3 * (x18 + x58 + 3 * x84 + x86)
    x88 = x2 * (x51 + x70)
    x89 = x3 * (5 * x59 + x60 + 4 * x85)
    x90 = x2 * x75
    x91 = x89 + x90
    x92 = 2 * x84
    x93 = x2 * x24
    x94 = 2 * x93
    x95 = x3 * (x25 + x32 + x92 + x94)
    x96 = x39 + x78
    x97 = x2 * x96
    x98 = 3 * x95 + 3 * x97
    x99 = x3 * (2 * x33 + 2 * x57 + 2 * x87 + 2 * x88 + x98)
    x100 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x101 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x102 = numpy.pi * x0 * x101
    x103 = x100 * x102
    x104 = -x0 * (a * A[1] + b * B[1])
    x105 = -x104 - B[1]
    x106 = x83 + x99
    x107 = x95 + x97
    x108 = x107 * x2
    x109 = x15 * x2
    x110 = 2 * x109
    x111 = x2 * x48
    x112 = x111 + x47
    x113 = x3 * (x110 + x112 + x16 + x46)
    x114 = x2 * (x29 + x84)
    x115 = x3 * (3 * x111 + x49 + x72)
    x116 = x2 * x86
    x117 = x115 + x116
    x118 = x31 + x93
    x119 = x118 * x2
    x120 = x2 * x22
    x121 = x3 * (x110 + x120 + x23 + x38)
    x122 = x3 * (2 * x113 + 2 * x114 + 2 * x119 + 2 * x121 + 2 * x39 + 2 * x78)
    x123 = x103 * (
        x106 * x2
        + x3
        * (
            3 * x108
            + 3 * x122
            + x19 * (3 * x113 + 3 * x114 + x117 + x52 + x77)
            + 3 * x69
            + 3 * x79
            + x81 * (x87 + x88)
        )
    )
    x124 = -x0 * (a * A[2] + b * B[2])
    x125 = -x124 - B[2]
    x126 = x100 * x7
    x127 = x126 * x3
    x128 = x105 ** 2 * x126
    x129 = x127 + x128
    x130 = x108 + x122
    x131 = x27 + x81 * x9
    x132 = x3 * (x131 + x26)
    x133 = x112 * x2
    x134 = x132 + x133
    x135 = x2 * x9
    x136 = x11 * x2
    x137 = x3 * (x135 + x136 + x14 + x27)
    x138 = x2 * (x109 + x12)
    x139 = 2 * x137 + 2 * x138
    x140 = x11 * x81 + x27
    x141 = x3 * (x140 + x21)
    x142 = x120 + x20
    x143 = x142 * x2
    x144 = x141 + x143
    x145 = x3 * (x139 + x144 + 2 * x31 + x94)
    x146 = x119 + x121
    x147 = x146 * x2
    x148 = x130 * x2 + x3 * (
        2 * x145 + 2 * x147 + x19 * (x134 + x139 + x30 + x92) + x81 * (x113 + x114) + x98
    )
    x149 = x101 * x7
    x150 = x103 * x125
    x151 = x149 * x3
    x152 = x125 ** 2 * x149
    x153 = x151 + x152
    x154 = x105 * x126
    x155 = x154 * x19
    x156 = x105 * x129
    x157 = x155 + x156
    x158 = x145 + x147
    x159 = x2 * x8
    x160 = x3 * (x11 + x159)
    x161 = x2 * (x13 + x136)
    x162 = x160 + x161
    x163 = x3 * (x159 + x9)
    x164 = x13 + x135
    x165 = x164 * x2
    x166 = x163 + x165
    x167 = x144 * x2 + x3 * (2 * x120 + 2 * x160 + 2 * x161 + x37)
    x168 = x158 * x2 + x3 * (
        3 * x119
        + 3 * x121
        + x167
        + x19 * (x110 + x162 + x166 + x36)
        + x81 * (x137 + x138)
    )
    x169 = x125 * x149
    x170 = x169 * x19
    x171 = x125 * x153
    x172 = x170 + x171
    x173 = 3 * x127
    x174 = x3 * (3 * x128 + x173)
    x175 = x105 * x157
    x176 = x174 + x175
    x177 = x2 ** 2 * x8
    x178 = x167 * x2 + x3 * (3 * x141 + 3 * x143 + x162 * x81 + x19 * (x140 + x177))
    x179 = 3 * x151
    x180 = x3 * (3 * x152 + x179)
    x181 = x125 * x172
    x182 = x180 + x181
    x183 = -x104 - A[1]
    x184 = x103 * x80
    x185 = x126 * x183
    x186 = x105 * x185
    x187 = x127 + x186
    x188 = x129 * x183
    x189 = x155 + x188
    x190 = x157 * x183
    x191 = x174 + x190
    x192 = x105 * x127
    x193 = 8 * x192
    x194 = x3 * (4 * x156 + x193)
    x195 = x176 * x183
    x196 = x194 + x195
    x197 = -x124 - A[2]
    x198 = x103 * x197
    x199 = x149 * x197
    x200 = x125 * x199
    x201 = x151 + x200
    x202 = x153 * x197
    x203 = x170 + x202
    x204 = x172 * x197
    x205 = x180 + x204
    x206 = x125 * x151
    x207 = 8 * x206
    x208 = x3 * (4 * x171 + x207)
    x209 = x182 * x197
    x210 = x208 + x209
    x211 = x126 * x183 ** 2
    x212 = x127 + x211
    x213 = x3 * (x154 + x185)
    x214 = x183 * x187
    x215 = x213 + x214
    x216 = 2 * x183
    x217 = x154 * x216 + x173
    x218 = x3 * (x128 + x217)
    x219 = x183 * x189
    x220 = x218 + x219
    x221 = x3 * (x156 + 3 * x188 + x193)
    x222 = x183 * x191
    x223 = x221 + x222
    x224 = x3 * (5 * x174 + x175 + 4 * x190)
    x225 = x183 * x196
    x226 = x224 + x225
    x227 = x149 * x197 ** 2
    x228 = x151 + x227
    x229 = x3 * (x169 + x199)
    x230 = x197 * x201
    x231 = x229 + x230
    x232 = 2 * x197
    x233 = x169 * x232 + x179
    x234 = x3 * (x152 + x233)
    x235 = x197 * x203
    x236 = x234 + x235
    x237 = x3 * (x171 + 3 * x202 + x207)
    x238 = x197 * x205
    x239 = x237 + x238
    x240 = x3 * (5 * x180 + x181 + 4 * x204)
    x241 = x197 * x210
    x242 = x240 + x241
    x243 = x183 * x212 + x185 * x19
    x244 = x3 * (x211 + x217)
    x245 = x183 * x215
    x246 = x244 + x245
    x247 = x183 * x220
    x248 = x3 * (2 * x188 + 4 * x192 + 2 * x213 + 2 * x214)
    x249 = x247 + x248
    x250 = x183 * x223
    x251 = 3 * x218 + 3 * x219
    x252 = x3 * (2 * x174 + 2 * x190 + x251)
    x253 = x250 + x252
    x254 = x183 * x226 + x3 * (2 * x194 + 2 * x195 + 4 * x221 + 4 * x222)
    x255 = x19 * x199 + x197 * x228
    x256 = x3 * (x227 + x233)
    x257 = x197 * x231
    x258 = x256 + x257
    x259 = x197 * x236
    x260 = x3 * (2 * x202 + 4 * x206 + 2 * x229 + 2 * x230)
    x261 = x259 + x260
    x262 = x197 * x239
    x263 = 3 * x234 + 3 * x235
    x264 = x3 * (2 * x180 + 2 * x204 + x263)
    x265 = x262 + x264
    x266 = x197 * x242 + x3 * (2 * x208 + 2 * x209 + 4 * x237 + 4 * x238)
    x267 = x183 * x243 + x3 * (x173 + 3 * x211)
    x268 = x183 * x246 + x3 * (3 * x213 + 3 * x214 + x243)
    x269 = x183 * x249 + x3 * (2 * x244 + 2 * x245 + x251)
    x270 = x183 * x253 + x3 * (3 * x221 + 3 * x222 + 3 * x247 + 3 * x248)
    x271 = x183 * x254 + x3 * (3 * x224 + 3 * x225 + 4 * x250 + 4 * x252)
    x272 = x197 * x255 + x3 * (x179 + 3 * x227)
    x273 = x197 * x258 + x3 * (3 * x229 + 3 * x230 + x255)
    x274 = x197 * x261 + x3 * (2 * x256 + 2 * x257 + x263)
    x275 = x197 * x265 + x3 * (3 * x237 + 3 * x238 + 3 * x259 + 3 * x260)
    x276 = x197 * x266 + x3 * (3 * x240 + 3 * x241 + 4 * x262 + 4 * x264)
    x277 = -x104 - C[1]
    x278 = x126 * x277 ** 2
    x279 = x127 + x278
    x280 = x2 * x91 + x3 * (4 * x115 + 4 * x116 + 2 * x73 + 2 * x74)
    x281 = x117 * x2
    x282 = 3 * x132 + 3 * x133
    x283 = x3 * (x282 + 2 * x59 + 2 * x85)
    x284 = x2 * x280 + x3 * (4 * x281 + 4 * x283 + 3 * x89 + 3 * x90)
    x285 = x126 * x277
    x286 = x19 * x285
    x287 = x105 * x279
    x288 = x286 + x287
    x289 = x281 + x283
    x290 = x134 * x2
    x291 = x3 * (2 * x111 + 2 * x163 + 2 * x165 + 4 * x71)
    x292 = x2 * x289 + x3 * (3 * x115 + 3 * x116 + 3 * x290 + 3 * x291)
    x293 = x154 * x277
    x294 = x173 + 2 * x293
    x295 = x3 * (x278 + x294)
    x296 = x105 * x288
    x297 = x295 + x296
    x298 = x3 * (x131 + x177)
    x299 = x166 * x2
    x300 = x290 + x291
    x301 = x2 * x300 + x3 * (x282 + 2 * x298 + 2 * x299)
    x302 = x105 * x297
    x303 = x127 + x293
    x304 = x105 * x303
    x305 = x3 * (x154 + x285)
    x306 = 2 * x305
    x307 = 4 * x127 * x277
    x308 = x306 + x307
    x309 = x3 * (2 * x287 + 2 * x304 + x308)
    x310 = x302 + x309
    x311 = x13 + x177
    x312 = x159 * x19 + x2 * x311
    x313 = x298 + x299
    x314 = x2 * x313 + x3 * (3 * x163 + 3 * x165 + x312)
    x315 = x304 + x305
    x316 = x105 * x315
    x317 = x3 * (x128 + x294)
    x318 = 2 * x317
    x319 = 3 * x295 + x318
    x320 = x3 * (3 * x296 + 2 * x316 + x319)
    x321 = x105 * x310
    x322 = x320 + x321
    x323 = x2 * x312 + x3 * (3 * x177 + x27)
    x324 = x183 * x279
    x325 = x286 + x324
    x326 = x183 * x288
    x327 = x295 + x326
    x328 = x183 * x297
    x329 = x309 + x328
    x330 = x183 * x310
    x331 = x320 + x330
    x332 = x183 * x322
    x333 = x316 + x317
    x334 = x105 * x333
    x335 = 3 * x305
    x336 = x3 * (x157 + 3 * x304 + x335)
    x337 = 2 * x336
    x338 = 4 * x309 + x337
    x339 = x3 * (4 * x302 + 2 * x334 + x338)
    x340 = x332 + x339
    x341 = x173 + x216 * x285
    x342 = x3 * (x278 + x341)
    x343 = x183 * x325
    x344 = x342 + x343
    x345 = x183 * x327
    x346 = x183 * x303
    x347 = 2 * x346
    x348 = x3 * (x287 + x308 + x324 + x347)
    x349 = x345 + x348
    x350 = x183 * x329
    x351 = x183 * x315
    x352 = 2 * x351
    x353 = 2 * x326
    x354 = x3 * (x296 + x319 + x352 + x353)
    x355 = x350 + x354
    x356 = x183 * x331
    x357 = x183 * x333
    x358 = 2 * x357
    x359 = x3 * (x302 + 3 * x328 + x338 + x358)
    x360 = x356 + x359
    x361 = x183 * x340
    x362 = 4 * x317
    x363 = x3 * (x176 + 4 * x316 + x362)
    x364 = x183 * (x334 + x336)
    x365 = 2 * x363 + 2 * x364
    x366 = x3 * (5 * x320 + x321 + 4 * x330 + x365)
    x367 = x361 + x366
    x368 = x3 * (x185 + x285)
    x369 = x185 * x277
    x370 = x183 * (x127 + x369)
    x371 = x183 * x344 + x3 * (x307 + 2 * x324 + 2 * x368 + 2 * x370)
    x372 = x183 * x349
    x373 = x3 * (x173 + x186 + x293 + x369)
    x374 = x183 * (x305 + x346)
    x375 = 2 * x373 + 2 * x374
    x376 = x3 * (2 * x295 + x344 + x353 + x375)
    x377 = x372 + x376
    x378 = x183 * x355
    x379 = x183 * (x317 + x351)
    x380 = x3 * (x189 + x304 + x335 + x347)
    x381 = x3 * (2 * x309 + 2 * x328 + 2 * x345 + 2 * x348 + 2 * x379 + 2 * x380)
    x382 = x378 + x381
    x383 = x183 * x360
    x384 = x3 * (x191 + x316 + 3 * x351 + x362)
    x385 = x183 * (x336 + x357)
    x386 = 3 * x350 + 3 * x354
    x387 = x3 * (2 * x320 + 2 * x330 + 2 * x384 + 2 * x385 + x386)
    x388 = x383 + x387
    x389 = x183 * (x363 + x364)
    x390 = x3 * (x196 + x334 + 5 * x336 + 4 * x357)
    x391 = x183 * x367 + x3 * (
        2 * x332 + 2 * x339 + 4 * x356 + 4 * x359 + 2 * x389 + 2 * x390
    )
    x392 = x102 * x6
    x393 = x391 * x392
    x394 = x2 * x392
    x395 = x368 + x370
    x396 = x183 * x371 + x3 * (x19 * (x211 + x341) + x216 * x395 + 3 * x342 + 3 * x343)
    x397 = x183 * x377 + x3 * (
        x19 * (x215 + x306 + x347 + x395)
        + x216 * (x373 + x374)
        + 3 * x345
        + 3 * x348
        + x371
    )
    x398 = x183 * x382 + x3 * (
        x19 * (x220 + x318 + x352 + x375)
        + x216 * (x379 + x380)
        + 2 * x372
        + 2 * x376
        + x386
    )
    x399 = x392 * (
        x183 * x388
        + x3
        * (
            x19 * (x223 + x337 + x358 + 3 * x379 + 3 * x380)
            + x216 * (x384 + x385)
            + 3 * x356
            + 3 * x359
            + 3 * x378
            + 3 * x381
        )
    )
    x400 = x392 * x4
    x401 = -x124 - C[2]
    x402 = x149 * x401 ** 2
    x403 = x151 + x402
    x404 = x149 * x401
    x405 = x19 * x404
    x406 = x125 * x403
    x407 = x405 + x406
    x408 = x169 * x401
    x409 = x179 + 2 * x408
    x410 = x3 * (x402 + x409)
    x411 = x125 * x407
    x412 = x410 + x411
    x413 = x125 * x412
    x414 = x151 + x408
    x415 = x125 * x414
    x416 = x3 * (x169 + x404)
    x417 = 2 * x416
    x418 = 4 * x151 * x401
    x419 = x417 + x418
    x420 = x3 * (2 * x406 + 2 * x415 + x419)
    x421 = x413 + x420
    x422 = x415 + x416
    x423 = x125 * x422
    x424 = x3 * (x152 + x409)
    x425 = 2 * x424
    x426 = 3 * x410 + x425
    x427 = x3 * (3 * x411 + 2 * x423 + x426)
    x428 = x125 * x421
    x429 = x427 + x428
    x430 = x197 * x403
    x431 = x405 + x430
    x432 = x197 * x407
    x433 = x410 + x432
    x434 = x197 * x412
    x435 = x420 + x434
    x436 = x197 * x421
    x437 = x427 + x436
    x438 = x197 * x429
    x439 = x423 + x424
    x440 = x125 * x439
    x441 = 3 * x416
    x442 = x3 * (x172 + 3 * x415 + x441)
    x443 = 2 * x442
    x444 = 4 * x420 + x443
    x445 = x3 * (4 * x413 + 2 * x440 + x444)
    x446 = x438 + x445
    x447 = x179 + x232 * x404
    x448 = x3 * (x402 + x447)
    x449 = x197 * x431
    x450 = x448 + x449
    x451 = x197 * x433
    x452 = x197 * x414
    x453 = 2 * x452
    x454 = x3 * (x406 + x419 + x430 + x453)
    x455 = x451 + x454
    x456 = x197 * x435
    x457 = x197 * x422
    x458 = 2 * x457
    x459 = 2 * x432
    x460 = x3 * (x411 + x426 + x458 + x459)
    x461 = x456 + x460
    x462 = x197 * x437
    x463 = x197 * x439
    x464 = 2 * x463
    x465 = x3 * (x413 + 3 * x434 + x444 + x464)
    x466 = x462 + x465
    x467 = x197 * x446
    x468 = 4 * x424
    x469 = x3 * (x182 + 4 * x423 + x468)
    x470 = x197 * (x440 + x442)
    x471 = 2 * x469 + 2 * x470
    x472 = x3 * (5 * x427 + x428 + 4 * x436 + x471)
    x473 = x467 + x472
    x474 = numpy.pi * x0 * x100 * x6
    x475 = x2 * x474
    x476 = x3 * (x199 + x404)
    x477 = x199 * x401
    x478 = x197 * (x151 + x477)
    x479 = x197 * x450 + x3 * (x418 + 2 * x430 + 2 * x476 + 2 * x478)
    x480 = x197 * x455
    x481 = x3 * (x179 + x200 + x408 + x477)
    x482 = x197 * (x416 + x452)
    x483 = 2 * x481 + 2 * x482
    x484 = x3 * (2 * x410 + x450 + x459 + x483)
    x485 = x480 + x484
    x486 = x197 * x461
    x487 = x197 * (x424 + x457)
    x488 = x3 * (x203 + x415 + x441 + x453)
    x489 = x3 * (2 * x420 + 2 * x434 + 2 * x451 + 2 * x454 + 2 * x487 + 2 * x488)
    x490 = x486 + x489
    x491 = x197 * x466
    x492 = x3 * (x205 + x423 + 3 * x457 + x468)
    x493 = x197 * (x442 + x463)
    x494 = 3 * x456 + 3 * x460
    x495 = x3 * (2 * x427 + 2 * x436 + 2 * x492 + 2 * x493 + x494)
    x496 = x491 + x495
    x497 = x197 * (x469 + x470)
    x498 = x3 * (x210 + x440 + 5 * x442 + 4 * x463)
    x499 = x197 * x473 + x3 * (
        2 * x438 + 2 * x445 + 4 * x462 + 4 * x465 + 2 * x497 + 2 * x498
    )
    x500 = x474 * x499
    x501 = x4 * x474
    x502 = x476 + x478
    x503 = x197 * x479 + x3 * (x19 * (x227 + x447) + x232 * x502 + 3 * x448 + 3 * x449)
    x504 = x197 * x485 + x3 * (
        x19 * (x231 + x417 + x453 + x502)
        + x232 * (x481 + x482)
        + 3 * x451
        + 3 * x454
        + x479
    )
    x505 = x197 * x490 + x3 * (
        x19 * (x236 + x425 + x458 + x483)
        + x232 * (x487 + x488)
        + 2 * x480
        + 2 * x484
        + x494
    )
    x506 = x474 * (
        x197 * x496
        + x3
        * (
            x19 * (x239 + x443 + x464 + 3 * x487 + 3 * x488)
            + x232 * (x492 + x493)
            + 3 * x462
            + 3 * x465
            + 3 * x486
            + 3 * x489
        )
    )

    # 675 item(s)
    return numpy.array(
        [
            x103
            * (
                x2 * x80
                + x3
                * (
                    x19 * (x64 + 4 * x87 + 4 * x88 + x91)
                    + 3 * x56
                    + 3 * x65
                    + x81 * (x67 + x76)
                    + 4 * x83
                    + 4 * x99
                )
            ),
            x105 * x123,
            x123 * x125,
            x129 * x148 * x149,
            x105 * x148 * x150,
            x126 * x148 * x153,
            x149 * x157 * x168,
            x129 * x168 * x169,
            x153 * x154 * x168,
            x126 * x168 * x172,
            x149 * x176 * x178,
            x157 * x169 * x178,
            x129 * x153 * x178,
            x154 * x172 * x178,
            x126 * x178 * x182,
            x183 * x184,
            x106 * x149 * x187,
            x106 * x150 * x183,
            x130 * x149 * x189,
            x130 * x169 * x187,
            x130 * x153 * x185,
            x149 * x158 * x191,
            x158 * x169 * x189,
            x153 * x158 * x187,
            x158 * x172 * x185,
            x149 * x167 * x196,
            x167 * x169 * x191,
            x153 * x167 * x189,
            x167 * x172 * x187,
            x167 * x182 * x185,
            x184 * x197,
            x105 * x106 * x198,
            x106 * x126 * x201,
            x129 * x130 * x199,
            x130 * x154 * x201,
            x126 * x130 * x203,
            x157 * x158 * x199,
            x129 * x158 * x201,
            x154 * x158 * x203,
            x126 * x158 * x205,
            x167 * x176 * x199,
            x157 * x167 * x201,
            x129 * x167 * x203,
            x154 * x167 * x205,
            x126 * x167 * x210,
            x149 * x212 * x66,
            x149 * x215 * x82,
            x169 * x212 * x82,
            x107 * x149 * x220,
            x107 * x169 * x215,
            x107 * x153 * x212,
            x146 * x149 * x223,
            x146 * x169 * x220,
            x146 * x153 * x215,
            x146 * x172 * x212,
            x144 * x149 * x226,
            x144 * x169 * x223,
            x144 * x153 * x220,
            x144 * x172 * x215,
            x144 * x182 * x212,
            x183 * x198 * x66,
            x187 * x199 * x82,
            x185 * x201 * x82,
            x107 * x189 * x199,
            x107 * x187 * x201,
            x107 * x185 * x203,
            x146 * x191 * x199,
            x146 * x189 * x201,
            x146 * x187 * x203,
            x146 * x185 * x205,
            x144 * x196 * x199,
            x144 * x191 * x201,
            x144 * x189 * x203,
            x144 * x187 * x205,
            x144 * x185 * x210,
            x126 * x228 * x66,
            x154 * x228 * x82,
            x126 * x231 * x82,
            x107 * x129 * x228,
            x107 * x154 * x231,
            x107 * x126 * x236,
            x146 * x157 * x228,
            x129 * x146 * x231,
            x146 * x154 * x236,
            x126 * x146 * x239,
            x144 * x176 * x228,
            x144 * x157 * x231,
            x129 * x144 * x236,
            x144 * x154 * x239,
            x126 * x144 * x242,
            x149 * x243 * x55,
            x149 * x246 * x68,
            x169 * x243 * x68,
            x149 * x249 * x96,
            x169 * x246 * x96,
            x153 * x243 * x96,
            x118 * x149 * x253,
            x118 * x169 * x249,
            x118 * x153 * x246,
            x118 * x172 * x243,
            x142 * x149 * x254,
            x142 * x169 * x253,
            x142 * x153 * x249,
            x142 * x172 * x246,
            x142 * x182 * x243,
            x199 * x212 * x55,
            x199 * x215 * x68,
            x201 * x212 * x68,
            x199 * x220 * x96,
            x201 * x215 * x96,
            x203 * x212 * x96,
            x118 * x199 * x223,
            x118 * x201 * x220,
            x118 * x203 * x215,
            x118 * x205 * x212,
            x142 * x199 * x226,
            x142 * x201 * x223,
            x142 * x203 * x220,
            x142 * x205 * x215,
            x142 * x210 * x212,
            x185 * x228 * x55,
            x187 * x228 * x68,
            x185 * x231 * x68,
            x189 * x228 * x96,
            x187 * x231 * x96,
            x185 * x236 * x96,
            x118 * x191 * x228,
            x118 * x189 * x231,
            x118 * x187 * x236,
            x118 * x185 * x239,
            x142 * x196 * x228,
            x142 * x191 * x231,
            x142 * x189 * x236,
            x142 * x187 * x239,
            x142 * x185 * x242,
            x126 * x255 * x55,
            x154 * x255 * x68,
            x126 * x258 * x68,
            x129 * x255 * x96,
            x154 * x258 * x96,
            x126 * x261 * x96,
            x118 * x157 * x255,
            x118 * x129 * x258,
            x118 * x154 * x261,
            x118 * x126 * x265,
            x142 * x176 * x255,
            x142 * x157 * x258,
            x129 * x142 * x261,
            x142 * x154 * x265,
            x126 * x142 * x266,
            x149 * x267 * x42,
            x149 * x268 * x40,
            x169 * x267 * x40,
            x149 * x269 * x34,
            x169 * x268 * x34,
            x153 * x267 * x34,
            x149 * x24 * x270,
            x169 * x24 * x269,
            x153 * x24 * x268,
            x172 * x24 * x267,
            x149 * x22 * x271,
            x169 * x22 * x270,
            x153 * x22 * x269,
            x172 * x22 * x268,
            x182 * x22 * x267,
            x199 * x243 * x42,
            x199 * x246 * x40,
            x201 * x243 * x40,
            x199 * x249 * x34,
            x201 * x246 * x34,
            x203 * x243 * x34,
            x199 * x24 * x253,
            x201 * x24 * x249,
            x203 * x24 * x246,
            x205 * x24 * x243,
            x199 * x22 * x254,
            x201 * x22 * x253,
            x203 * x22 * x249,
            x205 * x22 * x246,
            x210 * x22 * x243,
            x212 * x228 * x42,
            x215 * x228 * x40,
            x212 * x231 * x40,
            x220 * x228 * x34,
            x215 * x231 * x34,
            x212 * x236 * x34,
            x223 * x228 * x24,
            x220 * x231 * x24,
            x215 * x236 * x24,
            x212 * x239 * x24,
            x22 * x226 * x228,
            x22 * x223 * x231,
            x22 * x220 * x236,
            x215 * x22 * x239,
            x212 * x22 * x242,
            x185 * x255 * x42,
            x187 * x255 * x40,
            x185 * x258 * x40,
            x189 * x255 * x34,
            x187 * x258 * x34,
            x185 * x261 * x34,
            x191 * x24 * x255,
            x189 * x24 * x258,
            x187 * x24 * x261,
            x185 * x24 * x265,
            x196 * x22 * x255,
            x191 * x22 * x258,
            x189 * x22 * x261,
            x187 * x22 * x265,
            x185 * x22 * x266,
            x126 * x272 * x42,
            x154 * x272 * x40,
            x126 * x273 * x40,
            x129 * x272 * x34,
            x154 * x273 * x34,
            x126 * x274 * x34,
            x157 * x24 * x272,
            x129 * x24 * x273,
            x154 * x24 * x274,
            x126 * x24 * x275,
            x176 * x22 * x272,
            x157 * x22 * x273,
            x129 * x22 * x274,
            x154 * x22 * x275,
            x126 * x22 * x276,
            x149 * x279 * x284,
            x149 * x288 * x292,
            x169 * x279 * x292,
            x149 * x297 * x301,
            x169 * x288 * x301,
            x153 * x279 * x301,
            x149 * x310 * x314,
            x169 * x297 * x314,
            x153 * x288 * x314,
            x172 * x279 * x314,
            x149 * x322 * x323,
            x169 * x310 * x323,
            x153 * x297 * x323,
            x172 * x288 * x323,
            x182 * x279 * x323,
            x149 * x280 * x325,
            x149 * x289 * x327,
            x169 * x289 * x325,
            x149 * x300 * x329,
            x169 * x300 * x327,
            x153 * x300 * x325,
            x149 * x313 * x331,
            x169 * x313 * x329,
            x153 * x313 * x327,
            x172 * x313 * x325,
            x149 * x312 * x340,
            x169 * x312 * x331,
            x153 * x312 * x329,
            x172 * x312 * x327,
            x182 * x312 * x325,
            x199 * x279 * x280,
            x199 * x288 * x289,
            x201 * x279 * x289,
            x199 * x297 * x300,
            x201 * x288 * x300,
            x203 * x279 * x300,
            x199 * x310 * x313,
            x201 * x297 * x313,
            x203 * x288 * x313,
            x205 * x279 * x313,
            x199 * x312 * x322,
            x201 * x310 * x312,
            x203 * x297 * x312,
            x205 * x288 * x312,
            x210 * x279 * x312,
            x149 * x344 * x91,
            x117 * x149 * x349,
            x117 * x169 * x344,
            x134 * x149 * x355,
            x134 * x169 * x349,
            x134 * x153 * x344,
            x149 * x166 * x360,
            x166 * x169 * x355,
            x153 * x166 * x349,
            x166 * x172 * x344,
            x149 * x311 * x367,
            x169 * x311 * x360,
            x153 * x311 * x355,
            x172 * x311 * x349,
            x182 * x311 * x344,
            x199 * x325 * x91,
            x117 * x199 * x327,
            x117 * x201 * x325,
            x134 * x199 * x329,
            x134 * x201 * x327,
            x134 * x203 * x325,
            x166 * x199 * x331,
            x166 * x201 * x329,
            x166 * x203 * x327,
            x166 * x205 * x325,
            x199 * x311 * x340,
            x201 * x311 * x331,
            x203 * x311 * x329,
            x205 * x311 * x327,
            x210 * x311 * x325,
            x228 * x279 * x91,
            x117 * x228 * x288,
            x117 * x231 * x279,
            x134 * x228 * x297,
            x134 * x231 * x288,
            x134 * x236 * x279,
            x166 * x228 * x310,
            x166 * x231 * x297,
            x166 * x236 * x288,
            x166 * x239 * x279,
            x228 * x311 * x322,
            x231 * x310 * x311,
            x236 * x297 * x311,
            x239 * x288 * x311,
            x242 * x279 * x311,
            x149 * x371 * x75,
            x149 * x377 * x86,
            x169 * x371 * x86,
            x112 * x149 * x382,
            x112 * x169 * x377,
            x112 * x153 * x371,
            x149 * x164 * x388,
            x164 * x169 * x382,
            x153 * x164 * x377,
            x164 * x172 * x371,
            x2 * x393,
            x125 * x388 * x394,
            x153 * x159 * x382,
            x159 * x172 * x377,
            x159 * x182 * x371,
            x199 * x344 * x75,
            x199 * x349 * x86,
            x201 * x344 * x86,
            x112 * x199 * x355,
            x112 * x201 * x349,
            x112 * x203 * x344,
            x164 * x199 * x360,
            x164 * x201 * x355,
            x164 * x203 * x349,
            x164 * x205 * x344,
            x197 * x367 * x394,
            x159 * x201 * x360,
            x159 * x203 * x355,
            x159 * x205 * x349,
            x159 * x210 * x344,
            x228 * x325 * x75,
            x228 * x327 * x86,
            x231 * x325 * x86,
            x112 * x228 * x329,
            x112 * x231 * x327,
            x112 * x236 * x325,
            x164 * x228 * x331,
            x164 * x231 * x329,
            x164 * x236 * x327,
            x164 * x239 * x325,
            x159 * x228 * x340,
            x159 * x231 * x331,
            x159 * x236 * x329,
            x159 * x239 * x327,
            x159 * x242 * x325,
            x255 * x279 * x75,
            x255 * x288 * x86,
            x258 * x279 * x86,
            x112 * x255 * x297,
            x112 * x258 * x288,
            x112 * x261 * x279,
            x164 * x255 * x310,
            x164 * x258 * x297,
            x164 * x261 * x288,
            x164 * x265 * x279,
            x159 * x255 * x322,
            x159 * x258 * x310,
            x159 * x261 * x297,
            x159 * x265 * x288,
            x159 * x266 * x279,
            x149 * x396 * x61,
            x149 * x397 * x50,
            x169 * x396 * x50,
            x149 * x398 * x48,
            x169 * x397 * x48,
            x153 * x396 * x48,
            x399 * x4,
            x125 * x398 * x400,
            x153 * x397 * x9,
            x172 * x396 * x9,
            x392
            * (
                x183 * x391
                + x3
                * (
                    x19 * (x226 + x365 + 4 * x384 + 4 * x385)
                    + x216 * (x389 + x390)
                    + 3 * x361
                    + 3 * x366
                    + 4 * x383
                    + 4 * x387
                )
            ),
            x125 * x399,
            x153 * x398 * x8,
            x172 * x397 * x8,
            x182 * x396 * x8,
            x199 * x371 * x61,
            x199 * x377 * x50,
            x201 * x371 * x50,
            x199 * x382 * x48,
            x201 * x377 * x48,
            x203 * x371 * x48,
            x197 * x388 * x400,
            x201 * x382 * x9,
            x203 * x377 * x9,
            x205 * x371 * x9,
            x197 * x393,
            x201 * x388 * x8,
            x203 * x382 * x8,
            x205 * x377 * x8,
            x210 * x371 * x8,
            x228 * x344 * x61,
            x228 * x349 * x50,
            x231 * x344 * x50,
            x228 * x355 * x48,
            x231 * x349 * x48,
            x236 * x344 * x48,
            x228 * x360 * x9,
            x231 * x355 * x9,
            x236 * x349 * x9,
            x239 * x344 * x9,
            x228 * x367 * x8,
            x231 * x360 * x8,
            x236 * x355 * x8,
            x239 * x349 * x8,
            x242 * x344 * x8,
            x255 * x325 * x61,
            x255 * x327 * x50,
            x258 * x325 * x50,
            x255 * x329 * x48,
            x258 * x327 * x48,
            x261 * x325 * x48,
            x255 * x331 * x9,
            x258 * x329 * x9,
            x261 * x327 * x9,
            x265 * x325 * x9,
            x255 * x340 * x8,
            x258 * x331 * x8,
            x261 * x329 * x8,
            x265 * x327 * x8,
            x266 * x325 * x8,
            x272 * x279 * x61,
            x272 * x288 * x50,
            x273 * x279 * x50,
            x272 * x297 * x48,
            x273 * x288 * x48,
            x274 * x279 * x48,
            x272 * x310 * x9,
            x273 * x297 * x9,
            x274 * x288 * x9,
            x275 * x279 * x9,
            x272 * x322 * x8,
            x273 * x310 * x8,
            x274 * x297 * x8,
            x275 * x288 * x8,
            x276 * x279 * x8,
            x126 * x284 * x403,
            x154 * x292 * x403,
            x126 * x292 * x407,
            x129 * x301 * x403,
            x154 * x301 * x407,
            x126 * x301 * x412,
            x157 * x314 * x403,
            x129 * x314 * x407,
            x154 * x314 * x412,
            x126 * x314 * x421,
            x176 * x323 * x403,
            x157 * x323 * x407,
            x129 * x323 * x412,
            x154 * x323 * x421,
            x126 * x323 * x429,
            x185 * x280 * x403,
            x187 * x289 * x403,
            x185 * x289 * x407,
            x189 * x300 * x403,
            x187 * x300 * x407,
            x185 * x300 * x412,
            x191 * x313 * x403,
            x189 * x313 * x407,
            x187 * x313 * x412,
            x185 * x313 * x421,
            x196 * x312 * x403,
            x191 * x312 * x407,
            x189 * x312 * x412,
            x187 * x312 * x421,
            x185 * x312 * x429,
            x126 * x280 * x431,
            x154 * x289 * x431,
            x126 * x289 * x433,
            x129 * x300 * x431,
            x154 * x300 * x433,
            x126 * x300 * x435,
            x157 * x313 * x431,
            x129 * x313 * x433,
            x154 * x313 * x435,
            x126 * x313 * x437,
            x176 * x312 * x431,
            x157 * x312 * x433,
            x129 * x312 * x435,
            x154 * x312 * x437,
            x126 * x312 * x446,
            x212 * x403 * x91,
            x117 * x215 * x403,
            x117 * x212 * x407,
            x134 * x220 * x403,
            x134 * x215 * x407,
            x134 * x212 * x412,
            x166 * x223 * x403,
            x166 * x220 * x407,
            x166 * x215 * x412,
            x166 * x212 * x421,
            x226 * x311 * x403,
            x223 * x311 * x407,
            x220 * x311 * x412,
            x215 * x311 * x421,
            x212 * x311 * x429,
            x185 * x431 * x91,
            x117 * x187 * x431,
            x117 * x185 * x433,
            x134 * x189 * x431,
            x134 * x187 * x433,
            x134 * x185 * x435,
            x166 * x191 * x431,
            x166 * x189 * x433,
            x166 * x187 * x435,
            x166 * x185 * x437,
            x196 * x311 * x431,
            x191 * x311 * x433,
            x189 * x311 * x435,
            x187 * x311 * x437,
            x185 * x311 * x446,
            x126 * x450 * x91,
            x117 * x154 * x450,
            x117 * x126 * x455,
            x129 * x134 * x450,
            x134 * x154 * x455,
            x126 * x134 * x461,
            x157 * x166 * x450,
            x129 * x166 * x455,
            x154 * x166 * x461,
            x126 * x166 * x466,
            x176 * x311 * x450,
            x157 * x311 * x455,
            x129 * x311 * x461,
            x154 * x311 * x466,
            x126 * x311 * x473,
            x243 * x403 * x75,
            x246 * x403 * x86,
            x243 * x407 * x86,
            x112 * x249 * x403,
            x112 * x246 * x407,
            x112 * x243 * x412,
            x164 * x253 * x403,
            x164 * x249 * x407,
            x164 * x246 * x412,
            x164 * x243 * x421,
            x159 * x254 * x403,
            x159 * x253 * x407,
            x159 * x249 * x412,
            x159 * x246 * x421,
            x159 * x243 * x429,
            x212 * x431 * x75,
            x215 * x431 * x86,
            x212 * x433 * x86,
            x112 * x220 * x431,
            x112 * x215 * x433,
            x112 * x212 * x435,
            x164 * x223 * x431,
            x164 * x220 * x433,
            x164 * x215 * x435,
            x164 * x212 * x437,
            x159 * x226 * x431,
            x159 * x223 * x433,
            x159 * x220 * x435,
            x159 * x215 * x437,
            x159 * x212 * x446,
            x185 * x450 * x75,
            x187 * x450 * x86,
            x185 * x455 * x86,
            x112 * x189 * x450,
            x112 * x187 * x455,
            x112 * x185 * x461,
            x164 * x191 * x450,
            x164 * x189 * x455,
            x164 * x187 * x461,
            x164 * x185 * x466,
            x159 * x196 * x450,
            x159 * x191 * x455,
            x159 * x189 * x461,
            x159 * x187 * x466,
            x183 * x473 * x475,
            x126 * x479 * x75,
            x154 * x479 * x86,
            x126 * x485 * x86,
            x112 * x129 * x479,
            x112 * x154 * x485,
            x112 * x126 * x490,
            x157 * x164 * x479,
            x129 * x164 * x485,
            x154 * x164 * x490,
            x126 * x164 * x496,
            x159 * x176 * x479,
            x157 * x159 * x485,
            x129 * x159 * x490,
            x105 * x475 * x496,
            x2 * x500,
            x267 * x403 * x61,
            x268 * x403 * x50,
            x267 * x407 * x50,
            x269 * x403 * x48,
            x268 * x407 * x48,
            x267 * x412 * x48,
            x270 * x403 * x9,
            x269 * x407 * x9,
            x268 * x412 * x9,
            x267 * x421 * x9,
            x271 * x403 * x8,
            x270 * x407 * x8,
            x269 * x412 * x8,
            x268 * x421 * x8,
            x267 * x429 * x8,
            x243 * x431 * x61,
            x246 * x431 * x50,
            x243 * x433 * x50,
            x249 * x431 * x48,
            x246 * x433 * x48,
            x243 * x435 * x48,
            x253 * x431 * x9,
            x249 * x433 * x9,
            x246 * x435 * x9,
            x243 * x437 * x9,
            x254 * x431 * x8,
            x253 * x433 * x8,
            x249 * x435 * x8,
            x246 * x437 * x8,
            x243 * x446 * x8,
            x212 * x450 * x61,
            x215 * x450 * x50,
            x212 * x455 * x50,
            x220 * x450 * x48,
            x215 * x455 * x48,
            x212 * x461 * x48,
            x223 * x450 * x9,
            x220 * x455 * x9,
            x215 * x461 * x9,
            x212 * x466 * x9,
            x226 * x450 * x8,
            x223 * x455 * x8,
            x220 * x461 * x8,
            x215 * x466 * x8,
            x212 * x473 * x8,
            x185 * x479 * x61,
            x187 * x479 * x50,
            x185 * x485 * x50,
            x189 * x479 * x48,
            x187 * x48 * x485,
            x185 * x48 * x490,
            x191 * x479 * x9,
            x189 * x485 * x9,
            x187 * x490 * x9,
            x183 * x496 * x501,
            x196 * x479 * x8,
            x191 * x485 * x8,
            x189 * x490 * x8,
            x187 * x496 * x8,
            x183 * x500,
            x126 * x503 * x61,
            x154 * x50 * x503,
            x126 * x50 * x504,
            x129 * x48 * x503,
            x154 * x48 * x504,
            x126 * x48 * x505,
            x157 * x503 * x9,
            x129 * x504 * x9,
            x105 * x501 * x505,
            x4 * x506,
            x176 * x503 * x8,
            x157 * x504 * x8,
            x129 * x505 * x8,
            x105 * x506,
            x474
            * (
                x197 * x499
                + x3
                * (
                    x19 * (x242 + x471 + 4 * x492 + 4 * x493)
                    + x232 * (x497 + x498)
                    + 3 * x467
                    + 3 * x472
                    + 4 * x491
                    + 4 * x495
                )
            ),
        ]
    )
