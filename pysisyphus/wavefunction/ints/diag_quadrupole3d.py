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
    x3 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x4 = 1.77245385090552 * numpy.sqrt(x0)
    x5 = x3 * x4
    x6 = (2.0 * a + 2.0 * b) ** (-1.0)
    x7 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x8 = 3.14159265358979 * x0 * x7
    x9 = x2 * x4
    x10 = x4 * x7

    # 3 item(s)
    return numpy.array(
        [
            x2 * x8 * (x5 * x6 + x5 * (x0 * (a * A[0] + b * B[0]) - C[0]) ** 2),
            x3 * x8 * (x6 * x9 + x9 * (x0 * (a * A[1] + b * B[1]) - C[1]) ** 2),
            3.14159265358979
            * x0
            * x2
            * x3
            * (x10 * x6 + x10 * (x0 * (a * A[2] + b * B[2]) - C[2]) ** 2),
        ]
    )


def diag_quadrupole3d_01(a, A, b, B, C):
    """Cartesian 3D (sp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = -x0 * (a * A[0] + b * B[0])
    x2 = -x1 - B[0]
    x3 = -x1 - C[0]
    x4 = a * b * x0
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x0)
    x7 = x5 * x6
    x8 = (2.0 * a + 2.0 * b) ** (-1.0)
    x9 = x7 * x8
    x10 = x3**2 * x7 + x9
    x11 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x12 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x13 = 3.14159265358979 * x0 * x12
    x14 = x11 * x13
    x15 = -x0 * (a * A[1] + b * B[1])
    x16 = -x15 - B[1]
    x17 = x10 * x14
    x18 = -x0 * (a * A[2] + b * B[2])
    x19 = -x18 - B[2]
    x20 = -x15 - C[1]
    x21 = x11 * x6
    x22 = x21 * x8
    x23 = x20**2 * x21 + x22
    x24 = x13 * x5
    x25 = x23 * x24
    x26 = -x18 - C[2]
    x27 = x12 * x6
    x28 = x27 * x8
    x29 = x26**2 * x27 + x28
    x30 = 3.14159265358979 * x0 * x11 * x5
    x31 = x29 * x30

    # 9 item(s)
    return numpy.array(
        [
            x14 * (x10 * x2 + 2.0 * x3 * x9),
            x16 * x17,
            x17 * x19,
            x2 * x25,
            x24 * (x16 * x23 + 2.0 * x20 * x22),
            x19 * x25,
            x2 * x31,
            x16 * x31,
            x30 * (x19 * x29 + 2.0 * x26 * x28),
        ]
    )


def diag_quadrupole3d_02(a, A, b, B, C):
    """Cartesian 3D (sd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - C[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3**2 * x7
    x9 = x0 * x7
    x10 = -x2 - B[0]
    x11 = 2.0 * x3
    x12 = x8 + x9
    x13 = x10 * x12 + x11 * x9
    x14 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x15 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x16 = 3.14159265358979 * x1 * x15
    x17 = x14 * x16
    x18 = -x1 * (a * A[1] + b * B[1])
    x19 = -x18 - B[1]
    x20 = x13 * x17
    x21 = -x1 * (a * A[2] + b * B[2])
    x22 = -x21 - B[2]
    x23 = x14 * x6
    x24 = x0 * x23
    x25 = x19**2 * x23 + x24
    x26 = x15 * x6
    x27 = x0 * x26
    x28 = x22**2 * x26 + x27
    x29 = -x18 - C[1]
    x30 = x23 * x29**2
    x31 = x24 + x30
    x32 = x10**2 * x7 + x9
    x33 = 2.0 * x29
    x34 = x19 * x31 + x24 * x33
    x35 = x16 * x5
    x36 = x34 * x35
    x37 = -x21 - C[2]
    x38 = x26 * x37**2
    x39 = x27 + x38
    x40 = 3.14159265358979 * x1 * x14 * x5
    x41 = 2.0 * x37
    x42 = x22 * x39 + x27 * x41
    x43 = x40 * x42

    # 18 item(s)
    return numpy.array(
        [
            x17 * (x0 * (x10 * x11 * x7 + x8 + 3.0 * x9) + x10 * x13),
            x19 * x20,
            x20 * x22,
            x12 * x25 * x26,
            x12 * x17 * x19 * x22,
            x12 * x23 * x28,
            x26 * x31 * x32,
            x10 * x36,
            x10 * x22 * x31 * x35,
            x35 * (x0 * (x19 * x23 * x33 + 3.0 * x24 + x30) + x19 * x34),
            x22 * x36,
            x28 * x31 * x7,
            x23 * x32 * x39,
            x10 * x19 * x39 * x40,
            x10 * x43,
            x25 * x39 * x7,
            x19 * x43,
            x40 * (x0 * (x22 * x26 * x41 + 3.0 * x27 + x38) + x22 * x42),
        ]
    )


def diag_quadrupole3d_03(a, A, b, B, C):
    """Cartesian 3D (sf) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = 1.77245385090552 * numpy.sqrt(x1)
    x3 = -x1 * (a * A[0] + b * B[0])
    x4 = -x3 - B[0]
    x5 = a * b * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = x4 * x6
    x8 = x2 * x7
    x9 = -x3 - C[0]
    x10 = x2 * x6
    x11 = x10 * x9
    x12 = 2.0 * x0
    x13 = x10 * x9**2
    x14 = x0 * x10
    x15 = x13 + x14
    x16 = x15 * x4
    x17 = 2.0 * x4
    x18 = x11 * x12 + x16
    x19 = x0 * (x11 * x17 + x13 + 3.0 * x14) + x18 * x4
    x20 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x21 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x22 = 3.14159265358979 * x1 * x21
    x23 = x20 * x22
    x24 = -x1 * (a * A[1] + b * B[1])
    x25 = -x24 - B[1]
    x26 = x19 * x23
    x27 = -x1 * (a * A[2] + b * B[2])
    x28 = -x27 - B[2]
    x29 = x2 * x20
    x30 = x0 * x29
    x31 = x25**2 * x29 + x30
    x32 = x2 * x21
    x33 = x0 * x32
    x34 = x28**2 * x32 + x33
    x35 = x25 * x29
    x36 = x12 * x35 + x25 * x31
    x37 = x28 * x32
    x38 = x12 * x37 + x28 * x34
    x39 = -x24 - C[1]
    x40 = x29 * x39**2
    x41 = x30 + x40
    x42 = x10 * x4**2 + x14
    x43 = x12 * x8 + x4 * x42
    x44 = x25 * x41
    x45 = x29 * x39
    x46 = x12 * x45 + x44
    x47 = 2.0 * x25
    x48 = x0 * (3.0 * x30 + x40 + x45 * x47) + x25 * x46
    x49 = x22 * x7
    x50 = x22 * x6
    x51 = -x27 - C[2]
    x52 = x32 * x51**2
    x53 = x33 + x52
    x54 = x28 * x53
    x55 = x32 * x51
    x56 = x12 * x55 + x54
    x57 = 3.14159265358979 * x1 * x20
    x58 = x57 * x7
    x59 = 2.0 * x28
    x60 = x0 * (3.0 * x33 + x52 + x55 * x59) + x28 * x56
    x61 = x57 * x6

    # 30 item(s)
    return numpy.array(
        [
            x23
            * (
                x0
                * (4.0 * x0 * x11 + x12 * (x11 + x8) + 2.0 * x16 + x17 * (x14 + x8 * x9))
                + x19 * x4
            ),
            x25 * x26,
            x26 * x28,
            x18 * x31 * x32,
            x18 * x23 * x25 * x28,
            x18 * x29 * x34,
            x15 * x32 * x36,
            x15 * x31 * x37,
            x15 * x34 * x35,
            x15 * x29 * x38,
            x32 * x41 * x43,
            x32 * x42 * x46,
            x37 * x41 * x42,
            x48 * x49,
            x28 * x46 * x49,
            x34 * x41 * x8,
            x50
            * (
                x0
                * (
                    x12 * (x35 + x45)
                    + 4.0 * x30 * x39
                    + 2.0 * x44
                    + x47 * (x30 + x35 * x39)
                )
                + x25 * x48
            ),
            x28 * x48 * x50,
            x10 * x34 * x46,
            x10 * x38 * x41,
            x29 * x43 * x53,
            x35 * x42 * x53,
            x29 * x42 * x56,
            x31 * x53 * x8,
            x25 * x56 * x58,
            x58 * x60,
            x10 * x36 * x53,
            x10 * x31 * x56,
            x25 * x60 * x61,
            x61
            * (
                x0
                * (
                    x12 * (x37 + x55)
                    + 4.0 * x33 * x51
                    + 2.0 * x54
                    + x59 * (x33 + x37 * x51)
                )
                + x28 * x60
            ),
        ]
    )


def diag_quadrupole3d_04(a, A, b, B, C):
    """Cartesian 3D (sg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3**2 * x7
    x9 = x0 * x7
    x10 = 3.0 * x9
    x11 = 2.0 * x3
    x12 = -x2 - C[0]
    x13 = x12 * x7
    x14 = x10 + x11 * x13
    x15 = 2.0 * x0
    x16 = x3 * x7
    x17 = x0 * (x13 + x16)
    x18 = x3 * (x12 * x16 + x9)
    x19 = x12**2 * x7
    x20 = x0 * (x14 + x19)
    x21 = x19 + x9
    x22 = x21 * x3
    x23 = x13 * x15 + x22
    x24 = x23 * x3
    x25 = x20 + x24
    x26 = x0 * (4.0 * x0 * x13 + 2.0 * x17 + 2.0 * x18 + 2.0 * x22) + x25 * x3
    x27 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x28 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x29 = 3.14159265358979 * x1 * x28
    x30 = x27 * x29
    x31 = -x1 * (a * A[1] + b * B[1])
    x32 = -x31 - B[1]
    x33 = x26 * x30
    x34 = -x1 * (a * A[2] + b * B[2])
    x35 = -x34 - B[2]
    x36 = x27 * x6
    x37 = x32**2 * x36
    x38 = x0 * x36
    x39 = x37 + x38
    x40 = x28 * x6
    x41 = x35**2 * x40
    x42 = x0 * x40
    x43 = x41 + x42
    x44 = x32 * x36
    x45 = x15 * x44 + x32 * x39
    x46 = x35 * x40
    x47 = x15 * x46 + x35 * x43
    x48 = 3.0 * x38
    x49 = x0 * (3.0 * x37 + x48) + x32 * x45
    x50 = 3.0 * x42
    x51 = x0 * (3.0 * x41 + x50) + x35 * x47
    x52 = -x31 - C[1]
    x53 = x36 * x52**2
    x54 = x38 + x53
    x55 = x8 + x9
    x56 = x15 * x16 + x3 * x55
    x57 = x0 * (x10 + 3.0 * x8) + x3 * x56
    x58 = x32 * x54
    x59 = x36 * x52
    x60 = x15 * x59 + x58
    x61 = 2.0 * x32
    x62 = x48 + x59 * x61
    x63 = x0 * (x53 + x62)
    x64 = x32 * x60
    x65 = x63 + x64
    x66 = x0 * (x44 + x59)
    x67 = x32 * (x38 + x44 * x52)
    x68 = x0 * (4.0 * x38 * x52 + 2.0 * x58 + 2.0 * x66 + 2.0 * x67) + x32 * x65
    x69 = x29 * x5
    x70 = x68 * x69
    x71 = -x34 - C[2]
    x72 = x40 * x71**2
    x73 = x42 + x72
    x74 = x35 * x73
    x75 = x40 * x71
    x76 = x15 * x75 + x74
    x77 = 2.0 * x35
    x78 = x50 + x75 * x77
    x79 = x0 * (x72 + x78)
    x80 = x35 * x76
    x81 = x79 + x80
    x82 = 3.14159265358979 * x1 * x27 * x5
    x83 = x0 * (x46 + x75)
    x84 = x35 * (x42 + x46 * x71)
    x85 = x0 * (4.0 * x42 * x71 + 2.0 * x74 + 2.0 * x83 + 2.0 * x84) + x35 * x81
    x86 = x82 * x85

    # 45 item(s)
    return numpy.array(
        [
            x30
            * (
                x0 * (x11 * (x17 + x18) + x15 * (x14 + x8) + 3.0 * x20 + 3.0 * x24)
                + x26 * x3
            ),
            x32 * x33,
            x33 * x35,
            x25 * x39 * x40,
            x25 * x30 * x32 * x35,
            x25 * x36 * x43,
            x23 * x40 * x45,
            x23 * x39 * x46,
            x23 * x43 * x44,
            x23 * x36 * x47,
            x21 * x40 * x49,
            x21 * x45 * x46,
            x21 * x39 * x43,
            x21 * x44 * x47,
            x21 * x36 * x51,
            x40 * x54 * x57,
            x40 * x56 * x60,
            x46 * x54 * x56,
            x40 * x55 * x65,
            x46 * x55 * x60,
            x43 * x54 * x55,
            x3 * x70,
            x3 * x35 * x65 * x69,
            x16 * x43 * x60,
            x16 * x47 * x54,
            x69
            * (
                x0 * (x15 * (x37 + x62) + x61 * (x66 + x67) + 3.0 * x63 + 3.0 * x64)
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
            x16 * x45 * x73,
            x16 * x39 * x76,
            x3 * x32 * x81 * x82,
            x3 * x86,
            x49 * x7 * x73,
            x45 * x7 * x76,
            x39 * x7 * x81,
            x32 * x86,
            x82
            * (
                x0 * (x15 * (x41 + x78) + x77 * (x83 + x84) + 3.0 * x79 + 3.0 * x80)
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
    x2 = -x1 - A[0]
    x3 = -x1 - C[0]
    x4 = a * b * x0
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x0)
    x7 = x5 * x6
    x8 = (2.0 * a + 2.0 * b) ** (-1.0)
    x9 = x7 * x8
    x10 = x3**2 * x7 + x9
    x11 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x12 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x13 = 3.14159265358979 * x0 * x12
    x14 = x11 * x13
    x15 = -x0 * (a * A[1] + b * B[1])
    x16 = -x15 - A[1]
    x17 = x10 * x14
    x18 = -x0 * (a * A[2] + b * B[2])
    x19 = -x18 - A[2]
    x20 = -x15 - C[1]
    x21 = x11 * x6
    x22 = x21 * x8
    x23 = x20**2 * x21 + x22
    x24 = x13 * x5
    x25 = x23 * x24
    x26 = -x18 - C[2]
    x27 = x12 * x6
    x28 = x27 * x8
    x29 = x26**2 * x27 + x28
    x30 = 3.14159265358979 * x0 * x11 * x5
    x31 = x29 * x30

    # 9 item(s)
    return numpy.array(
        [
            x14 * (x10 * x2 + 2.0 * x3 * x9),
            x16 * x17,
            x17 * x19,
            x2 * x25,
            x24 * (x16 * x23 + 2.0 * x20 * x22),
            x19 * x25,
            x2 * x31,
            x16 * x31,
            x30 * (x19 * x29 + 2.0 * x26 * x28),
        ]
    )


def diag_quadrupole3d_11(a, A, b, B, C):
    """Cartesian 3D (pp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - C[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3**2 * x7
    x9 = x0 * x7
    x10 = -x2 - B[0]
    x11 = x10 * x7
    x12 = 2.0 * x3
    x13 = -x2 - A[0]
    x14 = x8 + x9
    x15 = x12 * x9
    x16 = x10 * x14 + x15
    x17 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x18 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x19 = 3.14159265358979 * x1 * x18
    x20 = x17 * x19
    x21 = -x1 * (a * A[1] + b * B[1])
    x22 = -x21 - B[1]
    x23 = x20 * (x13 * x14 + x15)
    x24 = -x1 * (a * A[2] + b * B[2])
    x25 = -x24 - B[2]
    x26 = -x21 - A[1]
    x27 = x16 * x20
    x28 = x0 * x6
    x29 = x17 * x28
    x30 = x17 * x6
    x31 = x22 * x30
    x32 = x26 * x31 + x29
    x33 = x18 * x6
    x34 = x14 * x20
    x35 = -x24 - A[2]
    x36 = x18 * x28
    x37 = x25 * x33
    x38 = x35 * x37 + x36
    x39 = -x21 - C[1]
    x40 = x30 * x39**2
    x41 = x29 + x40
    x42 = x11 * x13 + x9
    x43 = 2.0 * x39
    x44 = x29 * x43
    x45 = x22 * x41 + x44
    x46 = x19 * x5
    x47 = x45 * x46
    x48 = x41 * x46
    x49 = x46 * (x26 * x41 + x44)
    x50 = -x24 - C[2]
    x51 = x33 * x50**2
    x52 = x36 + x51
    x53 = 3.14159265358979 * x1 * x17 * x5
    x54 = x52 * x53
    x55 = 2.0 * x50
    x56 = x36 * x55
    x57 = x25 * x52 + x56
    x58 = x53 * x57
    x59 = x53 * (x35 * x52 + x56)

    # 27 item(s)
    return numpy.array(
        [
            x20 * (x0 * (x11 * x12 + x8 + 3.0 * x9) + x13 * x16),
            x22 * x23,
            x23 * x25,
            x26 * x27,
            x14 * x32 * x33,
            x25 * x26 * x34,
            x27 * x35,
            x22 * x34 * x35,
            x14 * x30 * x38,
            x33 * x41 * x42,
            x13 * x47,
            x13 * x25 * x48,
            x10 * x49,
            x46 * (x0 * (3.0 * x29 + x31 * x43 + x40) + x26 * x45),
            x25 * x49,
            x10 * x35 * x48,
            x35 * x47,
            x38 * x41 * x7,
            x30 * x42 * x52,
            x13 * x22 * x54,
            x13 * x58,
            x10 * x26 * x54,
            x32 * x52 * x7,
            x26 * x58,
            x10 * x59,
            x22 * x59,
            x53 * (x0 * (3.0 * x36 + x37 * x55 + x51) + x35 * x57),
        ]
    )


def diag_quadrupole3d_12(a, A, b, B, C):
    """Cartesian 3D (pd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = 1.77245385090552 * numpy.sqrt(x1)
    x3 = -x1 * (a * A[0] + b * B[0])
    x4 = -x3 - B[0]
    x5 = a * b * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = x4 * x6
    x8 = x2 * x7
    x9 = -x3 - C[0]
    x10 = x2 * x6
    x11 = x10 * x9
    x12 = 2.0 * x0
    x13 = x10 * x9**2
    x14 = x0 * x10
    x15 = x13 + x14
    x16 = x15 * x4
    x17 = 2.0 * x4
    x18 = -x3 - A[0]
    x19 = x0 * (x11 * x17 + x13 + 3.0 * x14)
    x20 = x11 * x12
    x21 = x16 + x20
    x22 = x19 + x21 * x4
    x23 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x24 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x25 = 3.14159265358979 * x1 * x24
    x26 = x23 * x25
    x27 = -x1 * (a * A[1] + b * B[1])
    x28 = -x27 - B[1]
    x29 = x26 * (x18 * x21 + x19)
    x30 = -x1 * (a * A[2] + b * B[2])
    x31 = -x30 - B[2]
    x32 = x15 * x18 + x20
    x33 = x2 * x23
    x34 = x0 * x33
    x35 = x28**2 * x33 + x34
    x36 = x2 * x24
    x37 = x35 * x36
    x38 = x26 * x31
    x39 = x0 * x36
    x40 = x31**2 * x36 + x39
    x41 = x33 * x40
    x42 = -x27 - A[1]
    x43 = x22 * x26
    x44 = x28 * x33
    x45 = x34 + x42 * x44
    x46 = x12 * x44 + x35 * x42
    x47 = x31 * x36
    x48 = -x30 - A[2]
    x49 = x39 + x47 * x48
    x50 = x12 * x47 + x40 * x48
    x51 = x10 * x4**2 + x14
    x52 = x12 * x8 + x18 * x51
    x53 = -x27 - C[1]
    x54 = x33 * x53**2
    x55 = x34 + x54
    x56 = x36 * x55
    x57 = x28 * x55
    x58 = x33 * x53
    x59 = x12 * x58
    x60 = x57 + x59
    x61 = x14 + x18 * x8
    x62 = 2.0 * x28
    x63 = x0 * (3.0 * x34 + x54 + x58 * x62)
    x64 = x28 * x60 + x63
    x65 = x25 * x6
    x66 = x64 * x65
    x67 = x31 * x65
    x68 = x10 * x40
    x69 = x42 * x55 + x59
    x70 = x42 * x60 + x63
    x71 = x25 * x7
    x72 = -x30 - C[2]
    x73 = x36 * x72**2
    x74 = x39 + x73
    x75 = x33 * x74
    x76 = x31 * x74
    x77 = x36 * x72
    x78 = x12 * x77
    x79 = x76 + x78
    x80 = x10 * x74
    x81 = 3.14159265358979 * x1 * x23
    x82 = x6 * x81
    x83 = x28 * x82
    x84 = 2.0 * x31
    x85 = x0 * (3.0 * x39 + x73 + x77 * x84)
    x86 = x31 * x79 + x85
    x87 = x82 * x86
    x88 = x7 * x81
    x89 = x48 * x74 + x78
    x90 = x48 * x79 + x85

    # 54 item(s)
    return numpy.array(
        [
            x26
            * (
                x0
                * (4.0 * x0 * x11 + x12 * (x11 + x8) + 2.0 * x16 + x17 * (x14 + x8 * x9))
                + x18 * x22
            ),
            x28 * x29,
            x29 * x31,
            x32 * x37,
            x28 * x32 * x38,
            x32 * x41,
            x42 * x43,
            x21 * x36 * x45,
            x21 * x38 * x42,
            x15 * x36 * x46,
            x15 * x45 * x47,
            x15 * x41 * x42,
            x43 * x48,
            x21 * x26 * x28 * x48,
            x21 * x33 * x49,
            x15 * x37 * x48,
            x15 * x44 * x49,
            x15 * x33 * x50,
            x52 * x56,
            x36 * x60 * x61,
            x47 * x55 * x61,
            x18 * x66,
            x18 * x60 * x67,
            x18 * x55 * x68,
            x36 * x51 * x69,
            x70 * x71,
            x31 * x69 * x71,
            x65
            * (
                x0
                * (
                    x12 * (x44 + x58)
                    + 4.0 * x34 * x53
                    + 2.0 * x57
                    + x62 * (x34 + x44 * x53)
                )
                + x42 * x64
            ),
            x67 * x70,
            x68 * x69,
            x48 * x51 * x56,
            x48 * x60 * x71,
            x49 * x55 * x8,
            x48 * x66,
            x10 * x49 * x60,
            x10 * x50 * x55,
            x52 * x75,
            x44 * x61 * x74,
            x33 * x61 * x79,
            x18 * x35 * x80,
            x18 * x79 * x83,
            x18 * x87,
            x42 * x51 * x75,
            x45 * x74 * x8,
            x42 * x79 * x88,
            x46 * x80,
            x10 * x45 * x79,
            x42 * x87,
            x33 * x51 * x89,
            x28 * x88 * x89,
            x88 * x90,
            x10 * x35 * x89,
            x83 * x90,
            x82
            * (
                x0
                * (
                    x12 * (x47 + x77)
                    + 4.0 * x39 * x72
                    + 2.0 * x76
                    + x84 * (x39 + x47 * x72)
                )
                + x48 * x86
            ),
        ]
    )


def diag_quadrupole3d_13(a, A, b, B, C):
    """Cartesian 3D (pf) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3**2 * x7
    x9 = x0 * x7
    x10 = 3.0 * x9
    x11 = 2.0 * x3
    x12 = -x2 - C[0]
    x13 = x12 * x7
    x14 = x10 + x11 * x13
    x15 = 2.0 * x0
    x16 = x3 * x7
    x17 = x0 * (x13 + x16)
    x18 = x3 * (x12 * x16 + x9)
    x19 = x12**2 * x7
    x20 = x0 * (x14 + x19)
    x21 = x19 + x9
    x22 = x21 * x3
    x23 = x13 * x15
    x24 = x22 + x23
    x25 = x24 * x3
    x26 = -x2 - A[0]
    x27 = x0 * (4.0 * x0 * x13 + 2.0 * x17 + 2.0 * x18 + 2.0 * x22)
    x28 = x20 + x25
    x29 = x27 + x28 * x3
    x30 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x31 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x32 = 3.14159265358979 * x1 * x31
    x33 = x30 * x32
    x34 = -x1 * (a * A[1] + b * B[1])
    x35 = -x34 - B[1]
    x36 = x33 * (x26 * x28 + x27)
    x37 = -x1 * (a * A[2] + b * B[2])
    x38 = -x37 - B[2]
    x39 = x20 + x24 * x26
    x40 = x30 * x6
    x41 = x35**2 * x40
    x42 = x0 * x40
    x43 = x41 + x42
    x44 = x31 * x6
    x45 = x43 * x44
    x46 = x33 * x38
    x47 = x38**2 * x44
    x48 = x0 * x44
    x49 = x47 + x48
    x50 = x40 * x49
    x51 = x21 * x26 + x23
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
    x62 = x42 + x52 * x60
    x63 = x43 * x60 + x53
    x64 = 3.0 * x42
    x65 = x0 * (3.0 * x41 + x64) + x54 * x60
    x66 = -x37 - A[2]
    x67 = x48 + x56 * x66
    x68 = x49 * x66 + x57
    x69 = 3.0 * x48
    x70 = x0 * (3.0 * x47 + x69) + x58 * x66
    x71 = x8 + x9
    x72 = x15 * x16
    x73 = x3 * x71 + x72
    x74 = x0 * (x10 + 3.0 * x8) + x26 * x73
    x75 = -x34 - C[1]
    x76 = x40 * x75**2
    x77 = x42 + x76
    x78 = x44 * x77
    x79 = x26 * x71 + x72
    x80 = x35 * x77
    x81 = x40 * x75
    x82 = x15 * x81
    x83 = x80 + x82
    x84 = x44 * x83
    x85 = 2.0 * x35
    x86 = x64 + x81 * x85
    x87 = x0 * (x76 + x86)
    x88 = x35 * x83
    x89 = x87 + x88
    x90 = x16 * x26 + x9
    x91 = x0 * (x52 + x81)
    x92 = x35 * (x42 + x52 * x75)
    x93 = x0 * (4.0 * x42 * x75 + 2.0 * x80 + 2.0 * x91 + 2.0 * x92)
    x94 = x35 * x89 + x93
    x95 = x32 * x5
    x96 = x94 * x95
    x97 = x38 * x95
    x98 = x49 * x7
    x99 = x58 * x7
    x100 = x60 * x77 + x82
    x101 = x60 * x83 + x87
    x102 = x95 * (x60 * x89 + x93)
    x103 = -x37 - C[2]
    x104 = x103**2 * x44
    x105 = x104 + x48
    x106 = x105 * x40
    x107 = x105 * x38
    x108 = x103 * x44
    x109 = x108 * x15
    x110 = x107 + x109
    x111 = x110 * x40
    x112 = 2.0 * x38
    x113 = x108 * x112 + x69
    x114 = x0 * (x104 + x113)
    x115 = x110 * x38
    x116 = x114 + x115
    x117 = x105 * x7
    x118 = x110 * x7
    x119 = 3.14159265358979 * x1 * x30 * x5
    x120 = x116 * x119
    x121 = x0 * (x108 + x56)
    x122 = x38 * (x103 * x56 + x48)
    x123 = x0 * (4.0 * x103 * x48 + 2.0 * x107 + 2.0 * x121 + 2.0 * x122)
    x124 = x116 * x38 + x123
    x125 = x119 * x124
    x126 = x105 * x66 + x109
    x127 = x110 * x66 + x114
    x128 = x119 * (x116 * x66 + x123)

    # 90 item(s)
    return numpy.array(
        [
            x33
            * (
                x0 * (x11 * (x17 + x18) + x15 * (x14 + x8) + 3.0 * x20 + 3.0 * x25)
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
            x28 * x44 * x62,
            x28 * x46 * x60,
            x24 * x44 * x63,
            x24 * x56 * x62,
            x24 * x50 * x60,
            x21 * x44 * x65,
            x21 * x56 * x63,
            x21 * x49 * x62,
            x21 * x59 * x60,
            x61 * x66,
            x28 * x33 * x35 * x66,
            x28 * x40 * x67,
            x24 * x45 * x66,
            x24 * x52 * x67,
            x24 * x40 * x68,
            x21 * x55 * x66,
            x21 * x43 * x67,
            x21 * x52 * x68,
            x21 * x40 * x70,
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
            x101 * x44 * x71,
            x100 * x56 * x71,
            x102 * x3,
            x101 * x3 * x97,
            x100 * x16 * x49,
            x95
            * (
                x0 * (x15 * (x41 + x86) + x85 * (x91 + x92) + 3.0 * x87 + 3.0 * x88)
                + x60 * x94
            ),
            x102 * x38,
            x101 * x98,
            x100 * x99,
            x66 * x73 * x78,
            x66 * x71 * x84,
            x67 * x71 * x77,
            x3 * x66 * x89 * x95,
            x16 * x67 * x83,
            x16 * x68 * x77,
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
            x105 * x62 * x71,
            x111 * x60 * x71,
            x105 * x16 * x63,
            x110 * x16 * x62,
            x120 * x3 * x60,
            x117 * x65,
            x118 * x63,
            x116 * x62 * x7,
            x125 * x60,
            x126 * x40 * x73,
            x126 * x52 * x71,
            x127 * x40 * x71,
            x126 * x16 * x43,
            x119 * x127 * x3 * x35,
            x128 * x3,
            x126 * x54 * x7,
            x127 * x43 * x7,
            x128 * x35,
            x119
            * (
                x0 * (x112 * (x121 + x122) + 3.0 * x114 + 3.0 * x115 + x15 * (x113 + x47))
                + x124 * x66
            ),
        ]
    )


def diag_quadrupole3d_14(a, A, b, B, C):
    """Cartesian 3D (pg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3 * x7
    x9 = -x2 - C[0]
    x10 = x7 * x9
    x11 = x0 * (x10 + x8)
    x12 = x0 * x7
    x13 = x3 * (x12 + x8 * x9)
    x14 = x3**2 * x7
    x15 = x12 + x14
    x16 = x15 * x3
    x17 = 2.0 * x0
    x18 = x17 * x8
    x19 = x16 + x18
    x20 = 3.0 * x12
    x21 = 2.0 * x3
    x22 = x10 * x21 + x20
    x23 = x0 * (x14 + x22)
    x24 = x3 * (x11 + x13)
    x25 = x7 * x9**2
    x26 = x12 + x25
    x27 = x26 * x3
    x28 = x0 * (2.0 * x11 + 4.0 * x12 * x9 + 2.0 * x13 + 2.0 * x27)
    x29 = x0 * (x22 + x25)
    x30 = x10 * x17
    x31 = x27 + x30
    x32 = x3 * x31
    x33 = x29 + x32
    x34 = x3 * x33
    x35 = -x2 - A[0]
    x36 = x0 * (2.0 * x23 + 2.0 * x24 + 3.0 * x29 + 3.0 * x32)
    x37 = x28 + x34
    x38 = x3 * x37 + x36
    x39 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x40 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x41 = 3.14159265358979 * x1 * x40
    x42 = x39 * x41
    x43 = -x1 * (a * A[1] + b * B[1])
    x44 = -x43 - B[1]
    x45 = x42 * (x35 * x37 + x36)
    x46 = -x1 * (a * A[2] + b * B[2])
    x47 = -x46 - B[2]
    x48 = x28 + x33 * x35
    x49 = x39 * x6
    x50 = x44**2 * x49
    x51 = x0 * x49
    x52 = x50 + x51
    x53 = x40 * x6
    x54 = x52 * x53
    x55 = x42 * x47
    x56 = x47**2 * x53
    x57 = x0 * x53
    x58 = x56 + x57
    x59 = x49 * x58
    x60 = x29 + x31 * x35
    x61 = x44 * x52
    x62 = x44 * x49
    x63 = x17 * x62
    x64 = x61 + x63
    x65 = x53 * x64
    x66 = x47 * x53
    x67 = x47 * x58
    x68 = x17 * x66
    x69 = x67 + x68
    x70 = x49 * x69
    x71 = x26 * x35 + x30
    x72 = 3.0 * x51
    x73 = x0 * (3.0 * x50 + x72)
    x74 = x44 * x64 + x73
    x75 = x53 * x74
    x76 = 3.0 * x57
    x77 = x0 * (3.0 * x56 + x76)
    x78 = x47 * x69 + x77
    x79 = x49 * x78
    x80 = -x43 - A[1]
    x81 = x38 * x42
    x82 = x51 + x62 * x80
    x83 = x52 * x80 + x63
    x84 = x64 * x80 + x73
    x85 = x0 * (8.0 * x44 * x51 + 4.0 * x61) + x74 * x80
    x86 = -x46 - A[2]
    x87 = x57 + x66 * x86
    x88 = x58 * x86 + x68
    x89 = x69 * x86 + x77
    x90 = x0 * (8.0 * x47 * x57 + 4.0 * x67) + x78 * x86
    x91 = x0 * (3.0 * x14 + x20)
    x92 = x19 * x3 + x91
    x93 = x0 * (8.0 * x12 * x3 + 4.0 * x16) + x35 * x92
    x94 = -x43 - C[1]
    x95 = x49 * x94**2
    x96 = x51 + x95
    x97 = x53 * x96
    x98 = x19 * x35 + x91
    x99 = x44 * x96
    x100 = x49 * x94
    x101 = x100 * x17
    x102 = x101 + x99
    x103 = x102 * x53
    x104 = x15 * x35 + x18
    x105 = 2.0 * x44
    x106 = x100 * x105 + x72
    x107 = x0 * (x106 + x95)
    x108 = x102 * x44
    x109 = x107 + x108
    x110 = x109 * x53
    x111 = x0 * (x100 + x62)
    x112 = x44 * (x51 + x62 * x94)
    x113 = x0 * (2.0 * x111 + 2.0 * x112 + 4.0 * x51 * x94 + 2.0 * x99)
    x114 = x109 * x44
    x115 = x113 + x114
    x116 = x12 + x35 * x8
    x117 = x0 * (x106 + x50)
    x118 = x44 * (x111 + x112)
    x119 = x0 * (3.0 * x107 + 3.0 * x108 + 2.0 * x117 + 2.0 * x118)
    x120 = x115 * x44 + x119
    x121 = x41 * x5
    x122 = x120 * x121
    x123 = x121 * x47
    x124 = x58 * x7
    x125 = x69 * x7
    x126 = x7 * x78
    x127 = x101 + x80 * x96
    x128 = x102 * x80 + x107
    x129 = x109 * x80 + x113
    x130 = x121 * (x115 * x80 + x119)
    x131 = -x46 - C[2]
    x132 = x131**2 * x53
    x133 = x132 + x57
    x134 = x133 * x49
    x135 = x133 * x47
    x136 = x131 * x53
    x137 = x136 * x17
    x138 = x135 + x137
    x139 = x138 * x49
    x140 = 2.0 * x47
    x141 = x136 * x140 + x76
    x142 = x0 * (x132 + x141)
    x143 = x138 * x47
    x144 = x142 + x143
    x145 = x144 * x49
    x146 = x0 * (x136 + x66)
    x147 = x47 * (x131 * x66 + x57)
    x148 = x0 * (4.0 * x131 * x57 + 2.0 * x135 + 2.0 * x146 + 2.0 * x147)
    x149 = x144 * x47
    x150 = x148 + x149
    x151 = x133 * x7
    x152 = x138 * x7
    x153 = x144 * x7
    x154 = 3.14159265358979 * x1 * x39 * x5
    x155 = x150 * x154
    x156 = x0 * (x141 + x56)
    x157 = x47 * (x146 + x147)
    x158 = x0 * (3.0 * x142 + 3.0 * x143 + 2.0 * x156 + 2.0 * x157)
    x159 = x150 * x47 + x158
    x160 = x154 * x159
    x161 = x133 * x86 + x137
    x162 = x138 * x86 + x142
    x163 = x144 * x86 + x148
    x164 = x154 * (x150 * x86 + x158)

    # 135 item(s)
    return numpy.array(
        [
            x42
            * (
                x0
                * (
                    x17 * (3.0 * x11 + 3.0 * x13 + x19)
                    + x21 * (x23 + x24)
                    + 4.0 * x28
                    + 4.0 * x34
                )
                + x35 * x38
            ),
            x44 * x45,
            x45 * x47,
            x48 * x54,
            x44 * x48 * x55,
            x48 * x59,
            x60 * x65,
            x52 * x60 * x66,
            x58 * x60 * x62,
            x60 * x70,
            x71 * x75,
            x64 * x66 * x71,
            x52 * x58 * x71,
            x62 * x69 * x71,
            x71 * x79,
            x80 * x81,
            x37 * x53 * x82,
            x37 * x55 * x80,
            x33 * x53 * x83,
            x33 * x66 * x82,
            x33 * x59 * x80,
            x31 * x53 * x84,
            x31 * x66 * x83,
            x31 * x58 * x82,
            x31 * x70 * x80,
            x26 * x53 * x85,
            x26 * x66 * x84,
            x26 * x58 * x83,
            x26 * x69 * x82,
            x26 * x79 * x80,
            x81 * x86,
            x37 * x42 * x44 * x86,
            x37 * x49 * x87,
            x33 * x54 * x86,
            x33 * x62 * x87,
            x33 * x49 * x88,
            x31 * x65 * x86,
            x31 * x52 * x87,
            x31 * x62 * x88,
            x31 * x49 * x89,
            x26 * x75 * x86,
            x26 * x64 * x87,
            x26 * x52 * x88,
            x26 * x62 * x89,
            x26 * x49 * x90,
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
            x122 * x35,
            x115 * x123 * x35,
            x109 * x124 * x35,
            x102 * x125 * x35,
            x126 * x35 * x96,
            x127 * x53 * x92,
            x128 * x19 * x53,
            x127 * x19 * x66,
            x129 * x15 * x53,
            x128 * x15 * x66,
            x127 * x15 * x58,
            x130 * x3,
            x123 * x129 * x3,
            x128 * x58 * x8,
            x127 * x69 * x8,
            x121
            * (
                x0
                * (
                    x105 * (x117 + x118)
                    + 4.0 * x113
                    + 4.0 * x114
                    + x17 * (3.0 * x111 + 3.0 * x112 + x64)
                )
                + x120 * x80
            ),
            x130 * x47,
            x124 * x129,
            x125 * x128,
            x126 * x127,
            x86 * x92 * x97,
            x103 * x19 * x86,
            x19 * x87 * x96,
            x110 * x15 * x86,
            x102 * x15 * x87,
            x15 * x88 * x96,
            x115 * x121 * x3 * x86,
            x109 * x8 * x87,
            x102 * x8 * x88,
            x8 * x89 * x96,
            x122 * x86,
            x115 * x7 * x87,
            x109 * x7 * x88,
            x102 * x7 * x89,
            x7 * x90 * x96,
            x134 * x93,
            x133 * x62 * x98,
            x139 * x98,
            x104 * x133 * x52,
            x104 * x138 * x62,
            x104 * x145,
            x116 * x133 * x64,
            x116 * x138 * x52,
            x116 * x144 * x62,
            x116 * x150 * x49,
            x151 * x35 * x74,
            x152 * x35 * x64,
            x153 * x35 * x52,
            x155 * x35 * x44,
            x160 * x35,
            x134 * x80 * x92,
            x133 * x19 * x82,
            x139 * x19 * x80,
            x133 * x15 * x83,
            x138 * x15 * x82,
            x145 * x15 * x80,
            x133 * x8 * x84,
            x138 * x8 * x83,
            x144 * x8 * x82,
            x155 * x3 * x80,
            x151 * x85,
            x152 * x84,
            x153 * x83,
            x150 * x7 * x82,
            x160 * x80,
            x161 * x49 * x92,
            x161 * x19 * x62,
            x162 * x19 * x49,
            x15 * x161 * x52,
            x15 * x162 * x62,
            x15 * x163 * x49,
            x161 * x64 * x8,
            x162 * x52 * x8,
            x154 * x163 * x3 * x44,
            x164 * x3,
            x161 * x7 * x74,
            x162 * x64 * x7,
            x163 * x52 * x7,
            x164 * x44,
            x154
            * (
                x0
                * (
                    x140 * (x156 + x157)
                    + 4.0 * x148
                    + 4.0 * x149
                    + x17 * (3.0 * x146 + 3.0 * x147 + x69)
                )
                + x159 * x86
            ),
        ]
    )


def diag_quadrupole3d_20(a, A, b, B, C):
    """Cartesian 3D (ds) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - C[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3**2 * x7
    x9 = x0 * x7
    x10 = -x2 - A[0]
    x11 = 2.0 * x3
    x12 = x8 + x9
    x13 = x10 * x12 + x11 * x9
    x14 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x15 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x16 = 3.14159265358979 * x1 * x15
    x17 = x14 * x16
    x18 = -x1 * (a * A[1] + b * B[1])
    x19 = -x18 - A[1]
    x20 = x13 * x17
    x21 = -x1 * (a * A[2] + b * B[2])
    x22 = -x21 - A[2]
    x23 = x14 * x6
    x24 = x0 * x23
    x25 = x19**2 * x23 + x24
    x26 = x15 * x6
    x27 = x0 * x26
    x28 = x22**2 * x26 + x27
    x29 = -x18 - C[1]
    x30 = x23 * x29**2
    x31 = x24 + x30
    x32 = x10**2 * x7 + x9
    x33 = 2.0 * x29
    x34 = x19 * x31 + x24 * x33
    x35 = x16 * x5
    x36 = x34 * x35
    x37 = -x21 - C[2]
    x38 = x26 * x37**2
    x39 = x27 + x38
    x40 = 3.14159265358979 * x1 * x14 * x5
    x41 = 2.0 * x37
    x42 = x22 * x39 + x27 * x41
    x43 = x40 * x42

    # 18 item(s)
    return numpy.array(
        [
            x17 * (x0 * (x10 * x11 * x7 + x8 + 3.0 * x9) + x10 * x13),
            x19 * x20,
            x20 * x22,
            x12 * x25 * x26,
            x12 * x17 * x19 * x22,
            x12 * x23 * x28,
            x26 * x31 * x32,
            x10 * x36,
            x10 * x22 * x31 * x35,
            x35 * (x0 * (x19 * x23 * x33 + 3.0 * x24 + x30) + x19 * x34),
            x22 * x36,
            x28 * x31 * x7,
            x23 * x32 * x39,
            x10 * x19 * x39 * x40,
            x10 * x43,
            x25 * x39 * x7,
            x19 * x43,
            x40 * (x0 * (x22 * x26 * x41 + 3.0 * x27 + x38) + x22 * x42),
        ]
    )


def diag_quadrupole3d_21(a, A, b, B, C):
    """Cartesian 3D (dp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - A[0]
    x4 = -x2 - C[0]
    x5 = a * b * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = 1.77245385090552 * numpy.sqrt(x1)
    x8 = x6 * x7
    x9 = x4**2 * x8
    x10 = x0 * x8
    x11 = x10 + x9
    x12 = x11 * x3
    x13 = -x2 - B[0]
    x14 = x11 * x13
    x15 = x13 * x8
    x16 = x4 * x8
    x17 = 2.0 * x0
    x18 = x15 * x4
    x19 = 2.0 * x3
    x20 = 3.0 * x10 + x9
    x21 = x16 * x17
    x22 = x14 + x21
    x23 = x0 * (2.0 * x18 + x20) + x22 * x3
    x24 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x25 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x26 = 3.14159265358979 * x1 * x25
    x27 = x24 * x26
    x28 = -x1 * (a * A[1] + b * B[1])
    x29 = -x28 - B[1]
    x30 = x12 + x21
    x31 = x27 * (x0 * (x16 * x19 + x20) + x3 * x30)
    x32 = -x1 * (a * A[2] + b * B[2])
    x33 = -x32 - B[2]
    x34 = -x28 - A[1]
    x35 = x23 * x27
    x36 = x0 * x7
    x37 = x24 * x36
    x38 = x24 * x7
    x39 = x34 * x38
    x40 = x29 * x39 + x37
    x41 = x25 * x7
    x42 = x27 * x30
    x43 = -x32 - A[2]
    x44 = x25 * x36
    x45 = x41 * x43
    x46 = x33 * x45 + x44
    x47 = x34**2 * x38 + x37
    x48 = x29 * x38
    x49 = x0 * (x39 + x48) + x34 * x40
    x50 = x33 * x41
    x51 = x41 * x43**2 + x44
    x52 = x0 * (x45 + x50) + x43 * x46
    x53 = -x28 - C[1]
    x54 = x38 * x53**2
    x55 = x37 + x54
    x56 = x3 * x8
    x57 = x10 + x15 * x3
    x58 = x0 * (x15 + x56) + x3 * x57
    x59 = x29 * x55
    x60 = x38 * x53
    x61 = x17 * x60
    x62 = x59 + x61
    x63 = x10 + x3**2 * x8
    x64 = x34 * x55
    x65 = x61 + x64
    x66 = x48 * x53
    x67 = 3.0 * x37 + x54
    x68 = x0 * (2.0 * x66 + x67) + x34 * x62
    x69 = x26 * x6
    x70 = x68 * x69
    x71 = x3 * x69
    x72 = 2.0 * x34
    x73 = x69 * (x0 * (x60 * x72 + x67) + x34 * x65)
    x74 = -x32 - C[2]
    x75 = x41 * x74**2
    x76 = x44 + x75
    x77 = x33 * x76
    x78 = x41 * x74
    x79 = x17 * x78
    x80 = x77 + x79
    x81 = 3.14159265358979 * x1 * x24 * x6
    x82 = x3 * x81
    x83 = x43 * x76
    x84 = x79 + x83
    x85 = x50 * x74
    x86 = 3.0 * x44 + x75
    x87 = x0 * (2.0 * x85 + x86) + x43 * x80
    x88 = x81 * x87
    x89 = 2.0 * x43
    x90 = x81 * (x0 * (x78 * x89 + x86) + x43 * x84)

    # 54 item(s)
    return numpy.array(
        [
            x27
            * (
                x0 * (4.0 * x10 * x4 + x12 + x14 + x17 * (x15 + x16) + x19 * (x10 + x18))
                + x23 * x3
            ),
            x29 * x31,
            x31 * x33,
            x34 * x35,
            x30 * x40 * x41,
            x33 * x34 * x42,
            x35 * x43,
            x29 * x42 * x43,
            x30 * x38 * x46,
            x22 * x41 * x47,
            x11 * x41 * x49,
            x11 * x47 * x50,
            x22 * x27 * x34 * x43,
            x11 * x40 * x45,
            x11 * x39 * x46,
            x22 * x38 * x51,
            x11 * x48 * x51,
            x11 * x38 * x52,
            x41 * x55 * x58,
            x41 * x62 * x63,
            x50 * x55 * x63,
            x41 * x57 * x65,
            x3 * x70,
            x33 * x65 * x71,
            x45 * x55 * x57,
            x43 * x62 * x71,
            x46 * x55 * x56,
            x13 * x73,
            x69
            * (
                x0 * (x17 * (x48 + x60) + 4.0 * x37 * x53 + x59 + x64 + x72 * (x37 + x66))
                + x34 * x68
            ),
            x33 * x73,
            x13 * x43 * x65 * x69,
            x43 * x70,
            x46 * x65 * x8,
            x15 * x51 * x55,
            x51 * x62 * x8,
            x52 * x55 * x8,
            x38 * x58 * x76,
            x48 * x63 * x76,
            x38 * x63 * x80,
            x39 * x57 * x76,
            x40 * x56 * x76,
            x34 * x80 * x82,
            x38 * x57 * x84,
            x29 * x82 * x84,
            x3 * x88,
            x15 * x47 * x76,
            x49 * x76 * x8,
            x47 * x8 * x80,
            x13 * x34 * x81 * x84,
            x40 * x8 * x84,
            x34 * x88,
            x13 * x90,
            x29 * x90,
            x81
            * (
                x0 * (x17 * (x50 + x78) + 4.0 * x44 * x74 + x77 + x83 + x89 * (x44 + x85))
                + x43 * x87
            ),
        ]
    )


def diag_quadrupole3d_22(a, A, b, B, C):
    """Cartesian 3D (dd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = -x2 - C[0]
    x5 = a * b * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = 1.77245385090552 * numpy.sqrt(x1)
    x8 = x6 * x7
    x9 = x4**2 * x8
    x10 = x0 * x8
    x11 = x10 + x9
    x12 = x11 * x3
    x13 = 2.0 * x0
    x14 = x4 * x8
    x15 = x13 * x14
    x16 = x12 + x15
    x17 = x16 * x3
    x18 = x3**2 * x8
    x19 = 3.0 * x10
    x20 = x3 * x8
    x21 = x20 * x4
    x22 = x19 + 2.0 * x21
    x23 = x0 * (x14 + x20)
    x24 = x10 + x21
    x25 = x24 * x3
    x26 = -x2 - A[0]
    x27 = 2.0 * x26
    x28 = x16 * x26
    x29 = x0 * (x22 + x9)
    x30 = 4.0 * x10 * x4 + 2.0 * x23
    x31 = x17 + x29
    x32 = x0 * (2.0 * x12 + 2.0 * x25 + x30) + x26 * x31
    x33 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x34 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x35 = 3.14159265358979 * x1 * x34
    x36 = x33 * x35
    x37 = -x1 * (a * A[1] + b * B[1])
    x38 = -x37 - B[1]
    x39 = x11 * x26
    x40 = x28 + x29
    x41 = x36 * (x0 * (x12 + x24 * x27 + x30 + x39) + x26 * x40)
    x42 = -x1 * (a * A[2] + b * B[2])
    x43 = -x42 - B[2]
    x44 = x33 * x7
    x45 = x38**2 * x44
    x46 = x0 * x44
    x47 = x45 + x46
    x48 = x15 + x39
    x49 = x0 * (x14 * x27 + x19 + x9) + x26 * x48
    x50 = x34 * x7
    x51 = x36 * x43
    x52 = x43**2 * x50
    x53 = x0 * x50
    x54 = x52 + x53
    x55 = -x37 - A[1]
    x56 = x32 * x36
    x57 = x44 * x55
    x58 = x38 * x57 + x46
    x59 = x38 * x44
    x60 = x13 * x59 + x47 * x55
    x61 = x43 * x50
    x62 = -x42 - A[2]
    x63 = x36 * x62
    x64 = x50 * x62
    x65 = x43 * x64 + x53
    x66 = x13 * x61 + x54 * x62
    x67 = x44 * x55**2 + x46
    x68 = x0 * (x57 + x59) + x55 * x58
    x69 = 2.0 * x55
    x70 = 3.0 * x46
    x71 = x45 + x70
    x72 = x0 * (x59 * x69 + x71) + x55 * x60
    x73 = x50 * x62**2 + x53
    x74 = x0 * (x61 + x64) + x62 * x65
    x75 = 2.0 * x62
    x76 = 3.0 * x53
    x77 = x52 + x76
    x78 = x0 * (x61 * x75 + x77) + x62 * x66
    x79 = -x37 - C[1]
    x80 = x44 * x79**2
    x81 = x46 + x80
    x82 = x10 + x18
    x83 = x13 * x20 + x26 * x82
    x84 = x0 * (x18 + x19 + x20 * x27) + x26 * x83
    x85 = x38 * x81
    x86 = x44 * x79
    x87 = x13 * x86
    x88 = x85 + x87
    x89 = x26 * x8
    x90 = x10 + x20 * x26
    x91 = x0 * (x20 + x89) + x26 * x90
    x92 = x59 * x79
    x93 = 2.0 * x92
    x94 = x70 + x80
    x95 = x0 * (x93 + x94)
    x96 = x38 * x88
    x97 = x95 + x96
    x98 = x10 + x26**2 * x8
    x99 = x55 * x81
    x100 = x87 + x99
    x101 = x55 * x88
    x102 = x101 + x95
    x103 = x46 + x92
    x104 = x103 * x38
    x105 = x0 * (x59 + x86)
    x106 = 2.0 * x105 + 4.0 * x46 * x79
    x107 = x0 * (2.0 * x104 + x106 + 2.0 * x85) + x55 * x97
    x108 = x35 * x6
    x109 = x107 * x108
    x110 = x108 * x26
    x111 = x0 * (x69 * x86 + x94) + x100 * x55
    x112 = x108 * (x0 * (x103 * x69 + x106 + x85 + x99) + x102 * x55)
    x113 = x108 * x3
    x114 = -x42 - C[2]
    x115 = x114**2 * x50
    x116 = x115 + x53
    x117 = x116 * x43
    x118 = x114 * x50
    x119 = x118 * x13
    x120 = x117 + x119
    x121 = x114 * x61
    x122 = 2.0 * x121
    x123 = x115 + x76
    x124 = x0 * (x122 + x123)
    x125 = x120 * x43
    x126 = x124 + x125
    x127 = 3.14159265358979 * x1 * x33 * x6
    x128 = x127 * x26
    x129 = x116 * x62
    x130 = x119 + x129
    x131 = x120 * x62
    x132 = x124 + x131
    x133 = x121 + x53
    x134 = x133 * x43
    x135 = x0 * (x118 + x61)
    x136 = 4.0 * x114 * x53 + 2.0 * x135
    x137 = x0 * (2.0 * x117 + 2.0 * x134 + x136) + x126 * x62
    x138 = x127 * x137
    x139 = x127 * x3
    x140 = x0 * (x118 * x75 + x123) + x130 * x62
    x141 = x127 * (x0 * (x117 + x129 + x133 * x75 + x136) + x132 * x62)

    # 108 item(s)
    return numpy.array(
        [
            x36
            * (
                x0 * (x13 * (x18 + x22) + x17 + x27 * (x23 + x25) + 2.0 * x28 + 3.0 * x29)
                + x26 * x32
            ),
            x38 * x41,
            x41 * x43,
            x47 * x49 * x50,
            x38 * x49 * x51,
            x44 * x49 * x54,
            x55 * x56,
            x40 * x50 * x58,
            x40 * x51 * x55,
            x48 * x50 * x60,
            x48 * x58 * x61,
            x48 * x54 * x57,
            x56 * x62,
            x38 * x40 * x63,
            x40 * x44 * x65,
            x47 * x48 * x64,
            x48 * x59 * x65,
            x44 * x48 * x66,
            x31 * x50 * x67,
            x16 * x50 * x68,
            x16 * x61 * x67,
            x11 * x50 * x72,
            x11 * x61 * x68,
            x11 * x54 * x67,
            x31 * x55 * x63,
            x16 * x58 * x64,
            x16 * x57 * x65,
            x11 * x60 * x64,
            x11 * x58 * x65,
            x11 * x57 * x66,
            x31 * x44 * x73,
            x16 * x59 * x73,
            x16 * x44 * x74,
            x11 * x47 * x73,
            x11 * x59 * x74,
            x11 * x44 * x78,
            x50 * x81 * x84,
            x50 * x88 * x91,
            x61 * x81 * x91,
            x50 * x97 * x98,
            x61 * x88 * x98,
            x54 * x81 * x98,
            x100 * x50 * x83,
            x102 * x50 * x90,
            x100 * x61 * x90,
            x109 * x26,
            x102 * x110 * x43,
            x100 * x54 * x89,
            x64 * x81 * x83,
            x64 * x88 * x90,
            x65 * x81 * x90,
            x110 * x62 * x97,
            x65 * x88 * x89,
            x66 * x81 * x89,
            x111 * x50 * x82,
            x112 * x3,
            x111 * x113 * x43,
            x108
            * (
                x0
                * (2.0 * x101 + x13 * (x71 + x93) + x69 * (x104 + x105) + 3.0 * x95 + x96)
                + x107 * x55
            ),
            x112 * x43,
            x111 * x54 * x8,
            x100 * x64 * x82,
            x102 * x113 * x62,
            x100 * x20 * x65,
            x109 * x62,
            x102 * x65 * x8,
            x100 * x66 * x8,
            x73 * x81 * x82,
            x20 * x73 * x88,
            x20 * x74 * x81,
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
            x138 * x26,
            x116 * x67 * x82,
            x116 * x20 * x68,
            x120 * x20 * x67,
            x116 * x72 * x8,
            x120 * x68 * x8,
            x126 * x67 * x8,
            x130 * x57 * x82,
            x130 * x20 * x58,
            x132 * x139 * x55,
            x130 * x60 * x8,
            x132 * x58 * x8,
            x138 * x55,
            x140 * x44 * x82,
            x139 * x140 * x38,
            x141 * x3,
            x140 * x47 * x8,
            x141 * x38,
            x127
            * (
                x0
                * (
                    3.0 * x124
                    + x125
                    + x13 * (x122 + x77)
                    + 2.0 * x131
                    + x75 * (x134 + x135)
                )
                + x137 * x62
            ),
        ]
    )


def diag_quadrupole3d_23(a, A, b, B, C):
    """Cartesian 3D (df) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = -x2 - C[0]
    x5 = a * b * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = 1.77245385090552 * numpy.sqrt(x1)
    x8 = x6 * x7
    x9 = x4**2 * x8
    x10 = x0 * x8
    x11 = 3.0 * x10
    x12 = x3 * x8
    x13 = x12 * x4
    x14 = x11 + 2.0 * x13
    x15 = x0 * (x14 + x9)
    x16 = x10 + x9
    x17 = x16 * x3
    x18 = 2.0 * x0
    x19 = x4 * x8
    x20 = x18 * x19
    x21 = x17 + x20
    x22 = x21 * x3
    x23 = x15 + x22
    x24 = x23 * x3
    x25 = x0 * (x12 + x19)
    x26 = x10 + x13
    x27 = x26 * x3
    x28 = x3**2 * x8
    x29 = x10 + x28
    x30 = x29 * x3
    x31 = x12 * x18
    x32 = x30 + x31
    x33 = x0 * (x14 + x28)
    x34 = x25 + x27
    x35 = x3 * x34
    x36 = -x2 - A[0]
    x37 = 2.0 * x36
    x38 = x23 * x36
    x39 = 4.0 * x10 * x4 + 2.0 * x25
    x40 = x0 * (2.0 * x17 + 2.0 * x27 + x39)
    x41 = 3.0 * x15 + 2.0 * x33
    x42 = x24 + x40
    x43 = x0 * (3.0 * x22 + 2.0 * x35 + x41) + x36 * x42
    x44 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x45 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x46 = 3.14159265358979 * x1 * x45
    x47 = x44 * x46
    x48 = -x1 * (a * A[1] + b * B[1])
    x49 = -x48 - B[1]
    x50 = x21 * x36
    x51 = x38 + x40
    x52 = x47 * (x0 * (x22 + x34 * x37 + x41 + 2.0 * x50) + x36 * x51)
    x53 = -x1 * (a * A[2] + b * B[2])
    x54 = -x53 - B[2]
    x55 = x44 * x7
    x56 = x49**2 * x55
    x57 = x0 * x55
    x58 = x56 + x57
    x59 = x16 * x36
    x60 = x15 + x50
    x61 = x0 * (x17 + x26 * x37 + x39 + x59) + x36 * x60
    x62 = x45 * x7
    x63 = x47 * x54
    x64 = x54**2 * x62
    x65 = x0 * x62
    x66 = x64 + x65
    x67 = x49 * x58
    x68 = x49 * x55
    x69 = x18 * x68
    x70 = x67 + x69
    x71 = x20 + x59
    x72 = x0 * (x11 + x19 * x37 + x9) + x36 * x71
    x73 = x54 * x62
    x74 = x54 * x66
    x75 = x18 * x73
    x76 = x74 + x75
    x77 = -x48 - A[1]
    x78 = x43 * x47
    x79 = x55 * x77
    x80 = x49 * x79 + x57
    x81 = x58 * x77
    x82 = x69 + x81
    x83 = 3.0 * x57
    x84 = x0 * (3.0 * x56 + x83) + x70 * x77
    x85 = -x53 - A[2]
    x86 = x47 * x85
    x87 = x62 * x85
    x88 = x54 * x87 + x65
    x89 = x66 * x85
    x90 = x75 + x89
    x91 = 3.0 * x65
    x92 = x0 * (3.0 * x64 + x91) + x76 * x85
    x93 = x55 * x77**2 + x57
    x94 = x0 * (x68 + x79) + x77 * x80
    x95 = 2.0 * x77
    x96 = x56 + x83
    x97 = x0 * (x68 * x95 + x96) + x77 * x82
    x98 = x0 * (8.0 * x49 * x57 + x67 + 3.0 * x81) + x77 * x84
    x99 = x62 * x85**2 + x65
    x100 = x0 * (x73 + x87) + x85 * x88
    x101 = 2.0 * x85
    x102 = x64 + x91
    x103 = x0 * (x101 * x73 + x102) + x85 * x90
    x104 = x0 * (8.0 * x54 * x65 + x74 + 3.0 * x89) + x85 * x92
    x105 = -x48 - C[1]
    x106 = x105**2 * x55
    x107 = x106 + x57
    x108 = x29 * x36
    x109 = x0 * (x11 + 3.0 * x28) + x32 * x36
    x110 = x0 * (8.0 * x10 * x3 + 3.0 * x108 + x30) + x109 * x36
    x111 = x107 * x49
    x112 = x105 * x55
    x113 = x112 * x18
    x114 = x111 + x113
    x115 = x108 + x31
    x116 = x0 * (x11 + x12 * x37 + x28) + x115 * x36
    x117 = x105 * x68
    x118 = 2.0 * x117
    x119 = x106 + x83
    x120 = x0 * (x118 + x119)
    x121 = x114 * x49
    x122 = x120 + x121
    x123 = x36 * x8
    x124 = x10 + x12 * x36
    x125 = x0 * (x12 + x123) + x124 * x36
    x126 = x117 + x57
    x127 = x126 * x49
    x128 = x0 * (x112 + x68)
    x129 = 4.0 * x105 * x57 + 2.0 * x128
    x130 = x0 * (2.0 * x111 + 2.0 * x127 + x129)
    x131 = x122 * x49
    x132 = x130 + x131
    x133 = x10 + x36**2 * x8
    x134 = x107 * x77
    x135 = x113 + x134
    x136 = x114 * x77
    x137 = x120 + x136
    x138 = x122 * x77
    x139 = x130 + x138
    x140 = x127 + x128
    x141 = x140 * x49
    x142 = x0 * (x118 + x96)
    x143 = 3.0 * x120 + 2.0 * x142
    x144 = x0 * (3.0 * x121 + 2.0 * x141 + x143) + x132 * x77
    x145 = x46 * x6
    x146 = x144 * x145
    x147 = x145 * x36
    x148 = x0 * (x112 * x95 + x119) + x135 * x77
    x149 = x0 * (x111 + x126 * x95 + x129 + x134) + x137 * x77
    x150 = x145 * (x0 * (x121 + 2.0 * x136 + x140 * x95 + x143) + x139 * x77)
    x151 = x145 * x3
    x152 = -x53 - C[2]
    x153 = x152**2 * x62
    x154 = x153 + x65
    x155 = x154 * x54
    x156 = x152 * x62
    x157 = x156 * x18
    x158 = x155 + x157
    x159 = x152 * x73
    x160 = 2.0 * x159
    x161 = x153 + x91
    x162 = x0 * (x160 + x161)
    x163 = x158 * x54
    x164 = x162 + x163
    x165 = x159 + x65
    x166 = x165 * x54
    x167 = x0 * (x156 + x73)
    x168 = 4.0 * x152 * x65 + 2.0 * x167
    x169 = x0 * (2.0 * x155 + 2.0 * x166 + x168)
    x170 = x164 * x54
    x171 = x169 + x170
    x172 = 3.14159265358979 * x1 * x44 * x6
    x173 = x172 * x36
    x174 = x154 * x85
    x175 = x157 + x174
    x176 = x158 * x85
    x177 = x162 + x176
    x178 = x164 * x85
    x179 = x169 + x178
    x180 = x166 + x167
    x181 = x180 * x54
    x182 = x0 * (x102 + x160)
    x183 = 3.0 * x162 + 2.0 * x182
    x184 = x0 * (3.0 * x163 + 2.0 * x181 + x183) + x171 * x85
    x185 = x172 * x184
    x186 = x172 * x3
    x187 = x0 * (x101 * x156 + x161) + x175 * x85
    x188 = x0 * (x101 * x165 + x155 + x168 + x174) + x177 * x85
    x189 = x172 * (x0 * (x101 * x180 + x163 + 2.0 * x176 + x183) + x179 * x85)

    # 180 item(s)
    return numpy.array(
        [
            x47
            * (
                x0
                * (
                    x18 * (3.0 * x25 + 3.0 * x27 + x32)
                    + x24
                    + x37 * (x33 + x35)
                    + 3.0 * x38
                    + 4.0 * x40
                )
                + x36 * x43
            ),
            x49 * x52,
            x52 * x54,
            x58 * x61 * x62,
            x49 * x61 * x63,
            x55 * x61 * x66,
            x62 * x70 * x72,
            x58 * x72 * x73,
            x66 * x68 * x72,
            x55 * x72 * x76,
            x77 * x78,
            x51 * x62 * x80,
            x51 * x63 * x77,
            x60 * x62 * x82,
            x60 * x73 * x80,
            x60 * x66 * x79,
            x62 * x71 * x84,
            x71 * x73 * x82,
            x66 * x71 * x80,
            x71 * x76 * x79,
            x78 * x85,
            x49 * x51 * x86,
            x51 * x55 * x88,
            x58 * x60 * x87,
            x60 * x68 * x88,
            x55 * x60 * x90,
            x70 * x71 * x87,
            x58 * x71 * x88,
            x68 * x71 * x90,
            x55 * x71 * x92,
            x42 * x62 * x93,
            x23 * x62 * x94,
            x23 * x73 * x93,
            x21 * x62 * x97,
            x21 * x73 * x94,
            x21 * x66 * x93,
            x16 * x62 * x98,
            x16 * x73 * x97,
            x16 * x66 * x94,
            x16 * x76 * x93,
            x42 * x77 * x86,
            x23 * x80 * x87,
            x23 * x79 * x88,
            x21 * x82 * x87,
            x21 * x80 * x88,
            x21 * x79 * x90,
            x16 * x84 * x87,
            x16 * x82 * x88,
            x16 * x80 * x90,
            x16 * x79 * x92,
            x42 * x55 * x99,
            x23 * x68 * x99,
            x100 * x23 * x55,
            x21 * x58 * x99,
            x100 * x21 * x68,
            x103 * x21 * x55,
            x16 * x70 * x99,
            x100 * x16 * x58,
            x103 * x16 * x68,
            x104 * x16 * x55,
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
            x146 * x36,
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
            x148 * x32 * x62,
            x149 * x29 * x62,
            x148 * x29 * x73,
            x150 * x3,
            x149 * x151 * x54,
            x12 * x148 * x66,
            x145
            * (
                x0
                * (
                    4.0 * x130
                    + x131
                    + 3.0 * x138
                    + x18 * (3.0 * x127 + 3.0 * x128 + x70)
                    + x95 * (x141 + x142)
                )
                + x144 * x77
            ),
            x150 * x54,
            x149 * x66 * x8,
            x148 * x76 * x8,
            x135 * x32 * x87,
            x137 * x29 * x87,
            x135 * x29 * x88,
            x139 * x151 * x85,
            x12 * x137 * x88,
            x12 * x135 * x90,
            x146 * x85,
            x139 * x8 * x88,
            x137 * x8 * x90,
            x135 * x8 * x92,
            x107 * x32 * x99,
            x114 * x29 * x99,
            x100 * x107 * x29,
            x12 * x122 * x99,
            x100 * x114 * x12,
            x103 * x107 * x12,
            x132 * x8 * x99,
            x100 * x122 * x8,
            x103 * x114 * x8,
            x104 * x107 * x8,
            x110 * x154 * x55,
            x116 * x154 * x68,
            x116 * x158 * x55,
            x125 * x154 * x58,
            x125 * x158 * x68,
            x125 * x164 * x55,
            x133 * x154 * x70,
            x133 * x158 * x58,
            x133 * x164 * x68,
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
            x115 * x175 * x68,
            x115 * x177 * x55,
            x124 * x175 * x58,
            x124 * x177 * x68,
            x124 * x179 * x55,
            x123 * x175 * x70,
            x123 * x177 * x58,
            x173 * x179 * x49,
            x185 * x36,
            x154 * x32 * x93,
            x154 * x29 * x94,
            x158 * x29 * x93,
            x12 * x154 * x97,
            x12 * x158 * x94,
            x12 * x164 * x93,
            x154 * x8 * x98,
            x158 * x8 * x97,
            x164 * x8 * x94,
            x171 * x8 * x93,
            x175 * x32 * x79,
            x175 * x29 * x80,
            x177 * x29 * x79,
            x12 * x175 * x82,
            x12 * x177 * x80,
            x179 * x186 * x77,
            x175 * x8 * x84,
            x177 * x8 * x82,
            x179 * x8 * x80,
            x185 * x77,
            x187 * x32 * x55,
            x187 * x29 * x68,
            x188 * x29 * x55,
            x12 * x187 * x58,
            x186 * x188 * x49,
            x189 * x3,
            x187 * x70 * x8,
            x188 * x58 * x8,
            x189 * x49,
            x172
            * (
                x0
                * (
                    x101 * (x181 + x182)
                    + 4.0 * x169
                    + x170
                    + 3.0 * x178
                    + x18 * (3.0 * x166 + 3.0 * x167 + x76)
                )
                + x184 * x85
            ),
        ]
    )


def diag_quadrupole3d_24(a, A, b, B, C):
    """Cartesian 3D (dg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = -x2 - C[0]
    x5 = a * b * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = 1.77245385090552 * numpy.sqrt(x1)
    x8 = x6 * x7
    x9 = x4**2 * x8
    x10 = x0 * x8
    x11 = x10 + x9
    x12 = x11 * x3
    x13 = x3 * x6
    x14 = x13 * x7
    x15 = x14 * x4
    x16 = x10 + x15
    x17 = x16 * x3
    x18 = x4 * x8
    x19 = x0 * (x14 + x18)
    x20 = 4.0 * x0 * x18 + 2.0 * x19
    x21 = x0 * (2.0 * x12 + 2.0 * x17 + x20)
    x22 = 3.0 * x10
    x23 = 2.0 * x15 + x22
    x24 = x0 * (x23 + x9)
    x25 = 2.0 * x0
    x26 = x18 * x25
    x27 = x12 + x26
    x28 = x27 * x3
    x29 = x24 + x28
    x30 = x29 * x3
    x31 = x21 + x30
    x32 = x3 * x31
    x33 = x3**2 * x8
    x34 = x0 * (x23 + x33)
    x35 = x17 + x19
    x36 = x3 * x35
    x37 = x0 * (x22 + 3.0 * x33)
    x38 = x10 + x33
    x39 = x3 * x38
    x40 = x14 * x25
    x41 = x39 + x40
    x42 = x3 * x41
    x43 = x37 + x42
    x44 = x0 * (3.0 * x17 + 3.0 * x19 + x41)
    x45 = x34 + x36
    x46 = x3 * x45
    x47 = -x2 - A[0]
    x48 = 2.0 * x47
    x49 = x31 * x47
    x50 = 3.0 * x24 + 2.0 * x34
    x51 = x0 * (3.0 * x28 + 2.0 * x36 + x50)
    x52 = 4.0 * x21 + 2.0 * x44
    x53 = x32 + x51
    x54 = x0 * (4.0 * x30 + 2.0 * x46 + x52) + x47 * x53
    x55 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x56 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x57 = 3.14159265358979 * x1 * x56
    x58 = x55 * x57
    x59 = -x1 * (a * A[1] + b * B[1])
    x60 = -x59 - B[1]
    x61 = x29 * x47
    x62 = x49 + x51
    x63 = x58 * (x0 * (x30 + x45 * x48 + x52 + 3.0 * x61) + x47 * x62)
    x64 = -x1 * (a * A[2] + b * B[2])
    x65 = -x64 - B[2]
    x66 = x55 * x7
    x67 = x60**2 * x66
    x68 = x0 * x66
    x69 = x67 + x68
    x70 = x27 * x47
    x71 = x21 + x61
    x72 = x0 * (x28 + x35 * x48 + x50 + 2.0 * x70) + x47 * x71
    x73 = x56 * x7
    x74 = x58 * x65
    x75 = x65**2 * x73
    x76 = x0 * x73
    x77 = x75 + x76
    x78 = x60 * x69
    x79 = x60 * x66
    x80 = x25 * x79
    x81 = x78 + x80
    x82 = x11 * x47
    x83 = x24 + x70
    x84 = x0 * (x12 + x16 * x48 + x20 + x82) + x47 * x83
    x85 = x65 * x73
    x86 = x65 * x77
    x87 = x25 * x85
    x88 = x86 + x87
    x89 = 3.0 * x68
    x90 = x0 * (3.0 * x67 + x89)
    x91 = x60 * x81
    x92 = x90 + x91
    x93 = x26 + x82
    x94 = x0 * (x18 * x48 + x22 + x9) + x47 * x93
    x95 = 3.0 * x76
    x96 = x0 * (3.0 * x75 + x95)
    x97 = x65 * x88
    x98 = x96 + x97
    x99 = -x59 - A[1]
    x100 = x54 * x58
    x101 = x66 * x99
    x102 = x101 * x60 + x68
    x103 = x69 * x99
    x104 = x103 + x80
    x105 = x81 * x99
    x106 = x105 + x90
    x107 = 8.0 * x60 * x68
    x108 = x0 * (x107 + 4.0 * x78) + x92 * x99
    x109 = -x64 - A[2]
    x110 = x109 * x58
    x111 = x109 * x73
    x112 = x111 * x65 + x76
    x113 = x109 * x77
    x114 = x113 + x87
    x115 = x109 * x88
    x116 = x115 + x96
    x117 = 8.0 * x65 * x76
    x118 = x0 * (x117 + 4.0 * x86) + x109 * x98
    x119 = x66 * x99**2 + x68
    x120 = x0 * (x101 + x79) + x102 * x99
    x121 = 2.0 * x99
    x122 = x67 + x89
    x123 = x0 * (x121 * x79 + x122) + x104 * x99
    x124 = x0 * (3.0 * x103 + x107 + x78) + x106 * x99
    x125 = x0 * (4.0 * x105 + 5.0 * x90 + x91) + x108 * x99
    x126 = x109**2 * x73 + x76
    x127 = x0 * (x111 + x85) + x109 * x112
    x128 = 2.0 * x109
    x129 = x75 + x95
    x130 = x0 * (x128 * x85 + x129) + x109 * x114
    x131 = x0 * (3.0 * x113 + x117 + x86) + x109 * x116
    x132 = x0 * (4.0 * x115 + 5.0 * x96 + x97) + x109 * x118
    x133 = -x59 - C[1]
    x134 = x133**2 * x66
    x135 = x134 + x68
    x136 = x41 * x47
    x137 = 8.0 * x0 * x14
    x138 = x0 * (x137 + 4.0 * x39) + x43 * x47
    x139 = x0 * (4.0 * x136 + 5.0 * x37 + x42) + x138 * x47
    x140 = x135 * x60
    x141 = x133 * x66
    x142 = x141 * x25
    x143 = x140 + x142
    x144 = x38 * x47
    x145 = x136 + x37
    x146 = x0 * (x137 + 3.0 * x144 + x39) + x145 * x47
    x147 = x133 * x79
    x148 = 2.0 * x147
    x149 = x134 + x89
    x150 = x0 * (x148 + x149)
    x151 = x143 * x60
    x152 = x150 + x151
    x153 = x144 + x40
    x154 = x0 * (x14 * x48 + x22 + x33) + x153 * x47
    x155 = x147 + x68
    x156 = x155 * x60
    x157 = x0 * (x141 + x79)
    x158 = 4.0 * x133 * x68 + 2.0 * x157
    x159 = x0 * (2.0 * x140 + 2.0 * x156 + x158)
    x160 = x152 * x60
    x161 = x159 + x160
    x162 = x47 * x8
    x163 = x10 + x14 * x47
    x164 = x0 * (x14 + x162) + x163 * x47
    x165 = x156 + x157
    x166 = x165 * x60
    x167 = x0 * (x122 + x148)
    x168 = 3.0 * x150 + 2.0 * x167
    x169 = x0 * (3.0 * x151 + 2.0 * x166 + x168)
    x170 = x161 * x60
    x171 = x169 + x170
    x172 = x10 + x47**2 * x8
    x173 = x135 * x99
    x174 = x142 + x173
    x175 = x143 * x99
    x176 = x150 + x175
    x177 = x152 * x99
    x178 = x159 + x177
    x179 = x161 * x99
    x180 = x169 + x179
    x181 = x166 + x167
    x182 = x181 * x60
    x183 = x0 * (3.0 * x156 + 3.0 * x157 + x81)
    x184 = 4.0 * x159 + 2.0 * x183
    x185 = x0 * (4.0 * x160 + 2.0 * x182 + x184) + x171 * x99
    x186 = x57 * x6
    x187 = x185 * x186
    x188 = x186 * x65
    x189 = x0 * (x121 * x141 + x149) + x174 * x99
    x190 = x0 * (x121 * x155 + x140 + x158 + x173) + x176 * x99
    x191 = x0 * (x121 * x165 + x151 + x168 + 2.0 * x175) + x178 * x99
    x192 = x0 * (x121 * x181 + x160 + 3.0 * x177 + x184) + x180 * x99
    x193 = x13 * x57
    x194 = -x64 - C[2]
    x195 = x194**2 * x73
    x196 = x195 + x76
    x197 = x196 * x65
    x198 = x194 * x73
    x199 = x198 * x25
    x200 = x197 + x199
    x201 = x194 * x85
    x202 = 2.0 * x201
    x203 = x195 + x95
    x204 = x0 * (x202 + x203)
    x205 = x200 * x65
    x206 = x204 + x205
    x207 = x201 + x76
    x208 = x207 * x65
    x209 = x0 * (x198 + x85)
    x210 = 4.0 * x194 * x76 + 2.0 * x209
    x211 = x0 * (2.0 * x197 + 2.0 * x208 + x210)
    x212 = x206 * x65
    x213 = x211 + x212
    x214 = x208 + x209
    x215 = x214 * x65
    x216 = x0 * (x129 + x202)
    x217 = 3.0 * x204 + 2.0 * x216
    x218 = x0 * (3.0 * x205 + 2.0 * x215 + x217)
    x219 = x213 * x65
    x220 = x218 + x219
    x221 = 3.14159265358979 * x1 * x55
    x222 = x221 * x6
    x223 = x109 * x196
    x224 = x199 + x223
    x225 = x109 * x200
    x226 = x204 + x225
    x227 = x109 * x206
    x228 = x211 + x227
    x229 = x109 * x213
    x230 = x218 + x229
    x231 = x222 * x60
    x232 = x215 + x216
    x233 = x232 * x65
    x234 = x0 * (3.0 * x208 + 3.0 * x209 + x88)
    x235 = 4.0 * x211 + 2.0 * x234
    x236 = x0 * (4.0 * x212 + 2.0 * x233 + x235) + x109 * x220
    x237 = x222 * x236
    x238 = x13 * x221
    x239 = x0 * (x128 * x198 + x203) + x109 * x224
    x240 = x0 * (x128 * x207 + x197 + x210 + x223) + x109 * x226
    x241 = x0 * (x128 * x214 + x205 + x217 + 2.0 * x225) + x109 * x228
    x242 = x0 * (x128 * x232 + x212 + 3.0 * x227 + x235) + x109 * x230

    # 270 item(s)
    return numpy.array(
        [
            x58
            * (
                x0
                * (
                    x25 * (4.0 * x34 + 4.0 * x36 + x43)
                    + x32
                    + x48 * (x44 + x46)
                    + 4.0 * x49
                    + 5.0 * x51
                )
                + x47 * x54
            ),
            x60 * x63,
            x63 * x65,
            x69 * x72 * x73,
            x60 * x72 * x74,
            x66 * x72 * x77,
            x73 * x81 * x84,
            x69 * x84 * x85,
            x77 * x79 * x84,
            x66 * x84 * x88,
            x73 * x92 * x94,
            x81 * x85 * x94,
            x69 * x77 * x94,
            x79 * x88 * x94,
            x66 * x94 * x98,
            x100 * x99,
            x102 * x62 * x73,
            x62 * x74 * x99,
            x104 * x71 * x73,
            x102 * x71 * x85,
            x101 * x71 * x77,
            x106 * x73 * x83,
            x104 * x83 * x85,
            x102 * x77 * x83,
            x101 * x83 * x88,
            x108 * x73 * x93,
            x106 * x85 * x93,
            x104 * x77 * x93,
            x102 * x88 * x93,
            x101 * x93 * x98,
            x100 * x109,
            x110 * x60 * x62,
            x112 * x62 * x66,
            x111 * x69 * x71,
            x112 * x71 * x79,
            x114 * x66 * x71,
            x111 * x81 * x83,
            x112 * x69 * x83,
            x114 * x79 * x83,
            x116 * x66 * x83,
            x111 * x92 * x93,
            x112 * x81 * x93,
            x114 * x69 * x93,
            x116 * x79 * x93,
            x118 * x66 * x93,
            x119 * x53 * x73,
            x120 * x31 * x73,
            x119 * x31 * x85,
            x123 * x29 * x73,
            x120 * x29 * x85,
            x119 * x29 * x77,
            x124 * x27 * x73,
            x123 * x27 * x85,
            x120 * x27 * x77,
            x119 * x27 * x88,
            x11 * x125 * x73,
            x11 * x124 * x85,
            x11 * x123 * x77,
            x11 * x120 * x88,
            x11 * x119 * x98,
            x110 * x53 * x99,
            x102 * x111 * x31,
            x101 * x112 * x31,
            x104 * x111 * x29,
            x102 * x112 * x29,
            x101 * x114 * x29,
            x106 * x111 * x27,
            x104 * x112 * x27,
            x102 * x114 * x27,
            x101 * x116 * x27,
            x108 * x11 * x111,
            x106 * x11 * x112,
            x104 * x11 * x114,
            x102 * x11 * x116,
            x101 * x11 * x118,
            x126 * x53 * x66,
            x126 * x31 * x79,
            x127 * x31 * x66,
            x126 * x29 * x69,
            x127 * x29 * x79,
            x130 * x29 * x66,
            x126 * x27 * x81,
            x127 * x27 * x69,
            x130 * x27 * x79,
            x131 * x27 * x66,
            x11 * x126 * x92,
            x11 * x127 * x81,
            x11 * x130 * x69,
            x11 * x131 * x79,
            x11 * x132 * x66,
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
            x187 * x47,
            x180 * x188 * x47,
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
            x109 * x171 * x186 * x47,
            x112 * x161 * x162,
            x114 * x152 * x162,
            x116 * x143 * x162,
            x118 * x135 * x162,
            x189 * x43 * x73,
            x190 * x41 * x73,
            x189 * x41 * x85,
            x191 * x38 * x73,
            x190 * x38 * x85,
            x189 * x38 * x77,
            x192 * x193,
            x191 * x193 * x65,
            x14 * x190 * x77,
            x14 * x189 * x88,
            x186
            * (
                x0
                * (
                    x121 * (x182 + x183)
                    + 5.0 * x169
                    + x170
                    + 4.0 * x179
                    + x25 * (4.0 * x166 + 4.0 * x167 + x92)
                )
                + x185 * x99
            ),
            x188 * x192,
            x191 * x77 * x8,
            x190 * x8 * x88,
            x189 * x8 * x98,
            x111 * x174 * x43,
            x111 * x176 * x41,
            x112 * x174 * x41,
            x111 * x178 * x38,
            x112 * x176 * x38,
            x114 * x174 * x38,
            x109 * x180 * x193,
            x112 * x14 * x178,
            x114 * x14 * x176,
            x116 * x14 * x174,
            x109 * x187,
            x112 * x180 * x8,
            x114 * x178 * x8,
            x116 * x176 * x8,
            x118 * x174 * x8,
            x126 * x135 * x43,
            x126 * x143 * x41,
            x127 * x135 * x41,
            x126 * x152 * x38,
            x127 * x143 * x38,
            x130 * x135 * x38,
            x126 * x14 * x161,
            x127 * x14 * x152,
            x130 * x14 * x143,
            x131 * x135 * x14,
            x126 * x171 * x8,
            x127 * x161 * x8,
            x130 * x152 * x8,
            x131 * x143 * x8,
            x132 * x135 * x8,
            x139 * x196 * x66,
            x146 * x196 * x79,
            x146 * x200 * x66,
            x154 * x196 * x69,
            x154 * x200 * x79,
            x154 * x206 * x66,
            x164 * x196 * x81,
            x164 * x200 * x69,
            x164 * x206 * x79,
            x164 * x213 * x66,
            x172 * x196 * x92,
            x172 * x200 * x81,
            x172 * x206 * x69,
            x172 * x213 * x79,
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
            x220 * x222 * x47 * x99,
            x138 * x224 * x66,
            x145 * x224 * x79,
            x145 * x226 * x66,
            x153 * x224 * x69,
            x153 * x226 * x79,
            x153 * x228 * x66,
            x163 * x224 * x81,
            x163 * x226 * x69,
            x163 * x228 * x79,
            x163 * x230 * x66,
            x162 * x224 * x92,
            x162 * x226 * x81,
            x162 * x228 * x69,
            x230 * x231 * x47,
            x237 * x47,
            x119 * x196 * x43,
            x120 * x196 * x41,
            x119 * x200 * x41,
            x123 * x196 * x38,
            x120 * x200 * x38,
            x119 * x206 * x38,
            x124 * x14 * x196,
            x123 * x14 * x200,
            x120 * x14 * x206,
            x119 * x14 * x213,
            x125 * x196 * x8,
            x124 * x200 * x8,
            x123 * x206 * x8,
            x120 * x213 * x8,
            x119 * x220 * x8,
            x101 * x224 * x43,
            x102 * x224 * x41,
            x101 * x226 * x41,
            x104 * x224 * x38,
            x102 * x226 * x38,
            x101 * x228 * x38,
            x106 * x14 * x224,
            x104 * x14 * x226,
            x102 * x14 * x228,
            x230 * x238 * x99,
            x108 * x224 * x8,
            x106 * x226 * x8,
            x104 * x228 * x8,
            x102 * x230 * x8,
            x237 * x99,
            x239 * x43 * x66,
            x239 * x41 * x79,
            x240 * x41 * x66,
            x239 * x38 * x69,
            x240 * x38 * x79,
            x241 * x38 * x66,
            x14 * x239 * x81,
            x14 * x240 * x69,
            x238 * x241 * x60,
            x238 * x242,
            x239 * x8 * x92,
            x240 * x8 * x81,
            x241 * x69 * x8,
            x231 * x242,
            x222
            * (
                x0
                * (
                    x128 * (x233 + x234)
                    + 5.0 * x218
                    + x219
                    + 4.0 * x229
                    + x25 * (4.0 * x215 + 4.0 * x216 + x98)
                )
                + x109 * x236
            ),
        ]
    )


def diag_quadrupole3d_30(a, A, b, B, C):
    """Cartesian 3D (fs) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = 1.77245385090552 * numpy.sqrt(x1)
    x3 = -x1 * (a * A[0] + b * B[0])
    x4 = -x3 - A[0]
    x5 = a * b * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = x4 * x6
    x8 = x2 * x7
    x9 = -x3 - C[0]
    x10 = x2 * x6
    x11 = x10 * x9
    x12 = 2.0 * x0
    x13 = x10 * x9**2
    x14 = x0 * x10
    x15 = x13 + x14
    x16 = x15 * x4
    x17 = 2.0 * x4
    x18 = x11 * x12 + x16
    x19 = x0 * (x11 * x17 + x13 + 3.0 * x14) + x18 * x4
    x20 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x21 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x22 = 3.14159265358979 * x1 * x21
    x23 = x20 * x22
    x24 = -x1 * (a * A[1] + b * B[1])
    x25 = -x24 - A[1]
    x26 = x19 * x23
    x27 = -x1 * (a * A[2] + b * B[2])
    x28 = -x27 - A[2]
    x29 = x2 * x20
    x30 = x0 * x29
    x31 = x25**2 * x29 + x30
    x32 = x2 * x21
    x33 = x0 * x32
    x34 = x28**2 * x32 + x33
    x35 = x25 * x29
    x36 = x12 * x35 + x25 * x31
    x37 = x28 * x32
    x38 = x12 * x37 + x28 * x34
    x39 = -x24 - C[1]
    x40 = x29 * x39**2
    x41 = x30 + x40
    x42 = x10 * x4**2 + x14
    x43 = x12 * x8 + x4 * x42
    x44 = x25 * x41
    x45 = x29 * x39
    x46 = x12 * x45 + x44
    x47 = 2.0 * x25
    x48 = x0 * (3.0 * x30 + x40 + x45 * x47) + x25 * x46
    x49 = x22 * x7
    x50 = x22 * x6
    x51 = -x27 - C[2]
    x52 = x32 * x51**2
    x53 = x33 + x52
    x54 = x28 * x53
    x55 = x32 * x51
    x56 = x12 * x55 + x54
    x57 = 3.14159265358979 * x1 * x20
    x58 = x57 * x7
    x59 = 2.0 * x28
    x60 = x0 * (3.0 * x33 + x52 + x55 * x59) + x28 * x56
    x61 = x57 * x6

    # 30 item(s)
    return numpy.array(
        [
            x23
            * (
                x0
                * (4.0 * x0 * x11 + x12 * (x11 + x8) + 2.0 * x16 + x17 * (x14 + x8 * x9))
                + x19 * x4
            ),
            x25 * x26,
            x26 * x28,
            x18 * x31 * x32,
            x18 * x23 * x25 * x28,
            x18 * x29 * x34,
            x15 * x32 * x36,
            x15 * x31 * x37,
            x15 * x34 * x35,
            x15 * x29 * x38,
            x32 * x41 * x43,
            x32 * x42 * x46,
            x37 * x41 * x42,
            x48 * x49,
            x28 * x46 * x49,
            x34 * x41 * x8,
            x50
            * (
                x0
                * (
                    x12 * (x35 + x45)
                    + 4.0 * x30 * x39
                    + 2.0 * x44
                    + x47 * (x30 + x35 * x39)
                )
                + x25 * x48
            ),
            x28 * x48 * x50,
            x10 * x34 * x46,
            x10 * x38 * x41,
            x29 * x43 * x53,
            x35 * x42 * x53,
            x29 * x42 * x56,
            x31 * x53 * x8,
            x25 * x56 * x58,
            x58 * x60,
            x10 * x36 * x53,
            x10 * x31 * x56,
            x25 * x60 * x61,
            x61
            * (
                x0
                * (
                    x12 * (x37 + x55)
                    + 4.0 * x33 * x51
                    + 2.0 * x54
                    + x59 * (x33 + x37 * x51)
                )
                + x28 * x60
            ),
        ]
    )


def diag_quadrupole3d_31(a, A, b, B, C):
    """Cartesian 3D (fp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - C[0]
    x4 = -x2 - B[0]
    x5 = a * b * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = 1.77245385090552 * numpy.sqrt(x1)
    x8 = x6 * x7
    x9 = x4 * x8
    x10 = x3 * x9
    x11 = x3**2 * x8
    x12 = x0 * x8
    x13 = 3.0 * x12
    x14 = x11 + x13
    x15 = x0 * (2.0 * x10 + x14)
    x16 = -x2 - A[0]
    x17 = x16 * x9
    x18 = x3 * x8
    x19 = x16 * x18
    x20 = 2.0 * x0
    x21 = x0 * (x18 + x9)
    x22 = x16 * (x10 + x12)
    x23 = 2.0 * x16
    x24 = x11 + x12
    x25 = x24 * x4
    x26 = x18 * x20
    x27 = x25 + x26
    x28 = x16 * x27
    x29 = x16 * x24
    x30 = x26 + x29
    x31 = x0 * (x14 + x18 * x23) + x16 * x30
    x32 = 4.0 * x12 * x3
    x33 = x15 + x28
    x34 = x0 * (2.0 * x21 + 2.0 * x22 + x25 + x29 + x32) + x16 * x33
    x35 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x36 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x37 = 3.14159265358979 * x1 * x36
    x38 = x35 * x37
    x39 = -x1 * (a * A[1] + b * B[1])
    x40 = -x39 - B[1]
    x41 = x16 * x8
    x42 = x38 * (
        x0 * (x20 * (x18 + x41) + x23 * (x12 + x19) + 2.0 * x29 + x32) + x16 * x31
    )
    x43 = -x1 * (a * A[2] + b * B[2])
    x44 = -x43 - B[2]
    x45 = -x39 - A[1]
    x46 = x34 * x38
    x47 = x0 * x7
    x48 = x35 * x47
    x49 = x35 * x7
    x50 = x45 * x49
    x51 = x40 * x50
    x52 = x48 + x51
    x53 = x36 * x7
    x54 = x31 * x38
    x55 = -x43 - A[2]
    x56 = x36 * x47
    x57 = x53 * x55
    x58 = x44 * x57
    x59 = x56 + x58
    x60 = x45**2 * x49
    x61 = x48 + x60
    x62 = x40 * x49
    x63 = x0 * (x50 + x62) + x45 * x52
    x64 = x44 * x53
    x65 = x53 * x55**2
    x66 = x56 + x65
    x67 = x0 * (x57 + x64) + x55 * x59
    x68 = x20 * x50 + x45 * x61
    x69 = 3.0 * x48
    x70 = 2.0 * x45
    x71 = x0 * (x60 + x62 * x70 + x69) + x45 * x63
    x72 = x20 * x57 + x55 * x66
    x73 = 3.0 * x56
    x74 = 2.0 * x55
    x75 = x0 * (x64 * x74 + x65 + x73) + x55 * x67
    x76 = -x39 - C[1]
    x77 = x49 * x76**2
    x78 = x48 + x77
    x79 = x16**2 * x8
    x80 = x12 + x17
    x81 = x0 * (x41 + x9) + x16 * x80
    x82 = x0 * (x13 + x23 * x9 + x79) + x16 * x81
    x83 = x40 * x78
    x84 = x49 * x76
    x85 = x20 * x84
    x86 = x83 + x85
    x87 = x12 + x79
    x88 = x16 * x87 + x20 * x41
    x89 = x45 * x78
    x90 = x85 + x89
    x91 = x62 * x76
    x92 = x69 + x77
    x93 = x0 * (2.0 * x91 + x92)
    x94 = x45 * x86
    x95 = x93 + x94
    x96 = x0 * (x70 * x84 + x92) + x45 * x90
    x97 = x0 * (x62 + x84)
    x98 = x45 * (x48 + x91)
    x99 = 4.0 * x48 * x76
    x100 = x0 * (x83 + x89 + 2.0 * x97 + 2.0 * x98 + x99) + x45 * x95
    x101 = x37 * x6
    x102 = x100 * x101
    x103 = x101 * x16
    x104 = x50 * x76
    x105 = x101 * (
        x0 * (x20 * (x50 + x84) + x70 * (x104 + x48) + 2.0 * x89 + x99) + x45 * x96
    )
    x106 = -x43 - C[2]
    x107 = x106**2 * x53
    x108 = x107 + x56
    x109 = x108 * x44
    x110 = x106 * x53
    x111 = x110 * x20
    x112 = x109 + x111
    x113 = x108 * x55
    x114 = x111 + x113
    x115 = x106 * x64
    x116 = x107 + x73
    x117 = x0 * (2.0 * x115 + x116)
    x118 = x112 * x55
    x119 = x117 + x118
    x120 = 3.14159265358979 * x1 * x35 * x6
    x121 = x120 * x16
    x122 = x0 * (x110 * x74 + x116) + x114 * x55
    x123 = x0 * (x110 + x64)
    x124 = x55 * (x115 + x56)
    x125 = 4.0 * x106 * x56
    x126 = x0 * (x109 + x113 + 2.0 * x123 + 2.0 * x124 + x125) + x119 * x55
    x127 = x120 * x126
    x128 = x106 * x57
    x129 = x120 * (
        x0 * (2.0 * x113 + x125 + x20 * (x110 + x57) + x74 * (x128 + x56)) + x122 * x55
    )

    # 90 item(s)
    return numpy.array(
        [
            x38
            * (
                x0
                * (
                    2.0 * x15
                    + x20 * (x10 + x13 + x17 + x19)
                    + x23 * (x21 + x22)
                    + 2.0 * x28
                    + x31
                )
                + x16 * x34
            ),
            x40 * x42,
            x42 * x44,
            x45 * x46,
            x31 * x52 * x53,
            x44 * x45 * x54,
            x46 * x55,
            x40 * x54 * x55,
            x31 * x49 * x59,
            x33 * x53 * x61,
            x30 * x53 * x63,
            x30 * x61 * x64,
            x33 * x38 * x45 * x55,
            x30 * x52 * x57,
            x30 * x50 * x59,
            x33 * x49 * x66,
            x30 * x62 * x66,
            x30 * x49 * x67,
            x27 * x53 * x68,
            x24 * x53 * x71,
            x24 * x64 * x68,
            x27 * x57 * x61,
            x24 * x57 * x63,
            x24 * x59 * x61,
            x27 * x50 * x66,
            x24 * x52 * x66,
            x24 * x50 * x67,
            x27 * x49 * x72,
            x24 * x62 * x72,
            x24 * x49 * x75,
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
            x102 * x16,
            x103 * x44 * x96,
            x57 * x80 * x90,
            x103 * x55 * x95,
            x41 * x59 * x90,
            x66 * x78 * x80,
            x41 * x66 * x86,
            x41 * x67 * x78,
            x105 * x4,
            x101
            * (
                x0
                * (
                    x20 * (x104 + x51 + x69 + x91)
                    + x70 * (x97 + x98)
                    + 2.0 * x93
                    + 2.0 * x94
                    + x96
                )
                + x100 * x45
            ),
            x105 * x44,
            x101 * x4 * x55 * x96,
            x102 * x55,
            x59 * x8 * x96,
            x66 * x9 * x90,
            x66 * x8 * x95,
            x67 * x8 * x90,
            x72 * x78 * x9,
            x72 * x8 * x86,
            x75 * x78 * x8,
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
            x127 * x16,
            x108 * x68 * x9,
            x108 * x71 * x8,
            x112 * x68 * x8,
            x114 * x61 * x9,
            x114 * x63 * x8,
            x119 * x61 * x8,
            x120 * x122 * x4 * x45,
            x122 * x52 * x8,
            x127 * x45,
            x129 * x4,
            x129 * x40,
            x120
            * (
                x0
                * (
                    2.0 * x117
                    + 2.0 * x118
                    + x122
                    + x20 * (x115 + x128 + x58 + x73)
                    + x74 * (x123 + x124)
                )
                + x126 * x55
            ),
        ]
    )


def diag_quadrupole3d_32(a, A, b, B, C):
    """Cartesian 3D (fd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = -x2 - C[0]
    x5 = a * b * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = 1.77245385090552 * numpy.sqrt(x1)
    x8 = x6 * x7
    x9 = x4**2 * x8
    x10 = x0 * x8
    x11 = x10 + x9
    x12 = x11 * x3
    x13 = x3 * x6
    x14 = x13 * x7
    x15 = x14 * x4
    x16 = x10 + x15
    x17 = x16 * x3
    x18 = x4 * x8
    x19 = x0 * (x14 + x18)
    x20 = 4.0 * x0
    x21 = x18 * x20
    x22 = 2.0 * x19 + x21
    x23 = x0 * (2.0 * x12 + 2.0 * x17 + x22)
    x24 = -x2 - A[0]
    x25 = x16 * x24
    x26 = 2.0 * x25
    x27 = x3**2 * x8
    x28 = x10 + x27
    x29 = x24 * x28
    x30 = 2.0 * x0
    x31 = x14 * x30 + x29
    x32 = x11 * x24
    x33 = x0 * (x12 + x22 + x26 + x32)
    x34 = 3.0 * x10
    x35 = 2.0 * x15 + x34
    x36 = x0 * (x27 + x35)
    x37 = x24 * (x17 + x19)
    x38 = 2.0 * x24
    x39 = x0 * (x35 + x9)
    x40 = x18 * x30
    x41 = x12 + x40
    x42 = x24 * x41
    x43 = x39 + x42
    x44 = x24 * x43
    x45 = x3 * x41
    x46 = x39 + x45
    x47 = x24 * x46
    x48 = 2.0 * x42
    x49 = x23 + x47
    x50 = x0 * (2.0 * x36 + 2.0 * x37 + 3.0 * x39 + x45 + x48) + x24 * x49
    x51 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x52 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x53 = 3.14159265358979 * x1 * x52
    x54 = x51 * x53
    x55 = -x1 * (a * A[1] + b * B[1])
    x56 = -x55 - B[1]
    x57 = x14 * x24
    x58 = x18 * x24
    x59 = x32 + x40
    x60 = x0 * (x18 * x38 + x34 + x9) + x24 * x59
    x61 = x33 + x44
    x62 = x54 * (
        x0 * (x30 * (x15 + x34 + x57 + x58) + x38 * (x19 + x25) + 2.0 * x39 + x48 + x60)
        + x24 * x61
    )
    x63 = -x1 * (a * A[2] + b * B[2])
    x64 = -x63 - B[2]
    x65 = x51 * x7
    x66 = x56**2 * x65
    x67 = x0 * x65
    x68 = x66 + x67
    x69 = x24 * x8
    x70 = x0 * (x21 + x30 * (x18 + x69) + 2.0 * x32 + x38 * (x10 + x58)) + x24 * x60
    x71 = x52 * x7
    x72 = x54 * x64
    x73 = x64**2 * x71
    x74 = x0 * x71
    x75 = x73 + x74
    x76 = -x55 - A[1]
    x77 = x50 * x54
    x78 = x65 * x76
    x79 = x56 * x78
    x80 = x67 + x79
    x81 = x68 * x76
    x82 = x56 * x65
    x83 = x30 * x82 + x81
    x84 = x64 * x71
    x85 = -x63 - A[2]
    x86 = x54 * x85
    x87 = x71 * x85
    x88 = x64 * x87
    x89 = x74 + x88
    x90 = x75 * x85
    x91 = x30 * x84 + x90
    x92 = x65 * x76**2
    x93 = x67 + x92
    x94 = x0 * (x78 + x82)
    x95 = x76 * x80
    x96 = x94 + x95
    x97 = 3.0 * x67
    x98 = 2.0 * x76
    x99 = x82 * x98 + x97
    x100 = x0 * (x66 + x99) + x76 * x83
    x101 = x71 * x85**2
    x102 = x101 + x74
    x103 = x0 * (x84 + x87)
    x104 = x85 * x89
    x105 = x103 + x104
    x106 = 3.0 * x74
    x107 = 2.0 * x85
    x108 = x106 + x107 * x84
    x109 = x0 * (x108 + x73) + x85 * x91
    x110 = x30 * x78 + x76 * x93
    x111 = x0 * (x92 + x99) + x76 * x96
    x112 = 4.0 * x67
    x113 = x0 * (x112 * x56 + 2.0 * x81 + 2.0 * x94 + 2.0 * x95) + x100 * x76
    x114 = x102 * x85 + x30 * x87
    x115 = x0 * (x101 + x108) + x105 * x85
    x116 = 4.0 * x74
    x117 = x0 * (2.0 * x103 + 2.0 * x104 + x116 * x64 + 2.0 * x90) + x109 * x85
    x118 = -x55 - C[1]
    x119 = x118**2 * x65
    x120 = x119 + x67
    x121 = x0 * (x14 + x69)
    x122 = x10 + x57
    x123 = x122 * x24
    x124 = x14 * x38 + x34
    x125 = x0 * (x124 + x27) + x24 * x31
    x126 = x0 * (2.0 * x121 + 2.0 * x123 + x14 * x20 + 2.0 * x29) + x125 * x24
    x127 = x120 * x56
    x128 = x118 * x65
    x129 = x128 * x30
    x130 = x127 + x129
    x131 = x24**2 * x8
    x132 = x121 + x123
    x133 = x0 * (x124 + x131) + x132 * x24
    x134 = x118 * x82
    x135 = 2.0 * x134
    x136 = x119 + x97
    x137 = x0 * (x135 + x136)
    x138 = x130 * x56
    x139 = x137 + x138
    x140 = x10 + x131
    x141 = x140 * x24 + x30 * x69
    x142 = x120 * x76
    x143 = x129 + x142
    x144 = x130 * x76
    x145 = x137 + x144
    x146 = x134 + x67
    x147 = x146 * x56
    x148 = x0 * (x128 + x82)
    x149 = x112 * x118
    x150 = 2.0 * x148 + x149
    x151 = x0 * (2.0 * x127 + 2.0 * x147 + x150)
    x152 = x139 * x76
    x153 = x151 + x152
    x154 = x0 * (x128 * x98 + x136) + x143 * x76
    x155 = x146 * x76
    x156 = 2.0 * x155
    x157 = x0 * (x127 + x142 + x150 + x156)
    x158 = x145 * x76
    x159 = x157 + x158
    x160 = x0 * (x135 + x66 + x97)
    x161 = x76 * (x147 + x148)
    x162 = 2.0 * x144
    x163 = x0 * (3.0 * x137 + x138 + 2.0 * x160 + 2.0 * x161 + x162) + x153 * x76
    x164 = x53 * x6
    x165 = x163 * x164
    x166 = x164 * x64
    x167 = x118 * x78
    x168 = x0 * (2.0 * x142 + x149 + x30 * (x128 + x78) + x98 * (x167 + x67)) + x154 * x76
    x169 = (
        x0
        * (
            2.0 * x137
            + x154
            + x162
            + x30 * (x134 + x167 + x79 + x97)
            + x98 * (x148 + x155)
        )
        + x159 * x76
    )
    x170 = x13 * x53
    x171 = -x63 - C[2]
    x172 = x171**2 * x71
    x173 = x172 + x74
    x174 = x173 * x64
    x175 = x171 * x71
    x176 = x175 * x30
    x177 = x174 + x176
    x178 = x171 * x84
    x179 = 2.0 * x178
    x180 = x106 + x172
    x181 = x0 * (x179 + x180)
    x182 = x177 * x64
    x183 = x181 + x182
    x184 = x173 * x85
    x185 = x176 + x184
    x186 = x177 * x85
    x187 = x181 + x186
    x188 = x178 + x74
    x189 = x188 * x64
    x190 = x0 * (x175 + x84)
    x191 = x116 * x171
    x192 = 2.0 * x190 + x191
    x193 = x0 * (2.0 * x174 + 2.0 * x189 + x192)
    x194 = x183 * x85
    x195 = x193 + x194
    x196 = 3.14159265358979 * x1 * x51
    x197 = x196 * x6
    x198 = x0 * (x107 * x175 + x180) + x185 * x85
    x199 = x188 * x85
    x200 = 2.0 * x199
    x201 = x0 * (x174 + x184 + x192 + x200)
    x202 = x187 * x85
    x203 = x201 + x202
    x204 = x197 * x56
    x205 = x0 * (x106 + x179 + x73)
    x206 = x85 * (x189 + x190)
    x207 = 2.0 * x186
    x208 = x0 * (3.0 * x181 + x182 + 2.0 * x205 + 2.0 * x206 + x207) + x195 * x85
    x209 = x197 * x208
    x210 = x13 * x196
    x211 = x171 * x87
    x212 = (
        x0 * (x107 * (x211 + x74) + 2.0 * x184 + x191 + x30 * (x175 + x87)) + x198 * x85
    )
    x213 = (
        x0
        * (
            x107 * (x190 + x199)
            + 2.0 * x181
            + x198
            + x207
            + x30 * (x106 + x178 + x211 + x88)
        )
        + x203 * x85
    )

    # 180 item(s)
    return numpy.array(
        [
            x54
            * (
                x0
                * (
                    2.0 * x23
                    + x30 * (x17 + 3.0 * x19 + x26 + x31)
                    + 2.0 * x33
                    + x38 * (x36 + x37)
                    + 2.0 * x44
                    + 2.0 * x47
                )
                + x24 * x50
            ),
            x56 * x62,
            x62 * x64,
            x68 * x70 * x71,
            x56 * x70 * x72,
            x65 * x70 * x75,
            x76 * x77,
            x61 * x71 * x80,
            x61 * x72 * x76,
            x60 * x71 * x83,
            x60 * x80 * x84,
            x60 * x75 * x78,
            x77 * x85,
            x56 * x61 * x86,
            x61 * x65 * x89,
            x60 * x68 * x87,
            x60 * x82 * x89,
            x60 * x65 * x91,
            x49 * x71 * x93,
            x43 * x71 * x96,
            x43 * x84 * x93,
            x100 * x59 * x71,
            x59 * x84 * x96,
            x59 * x75 * x93,
            x49 * x76 * x86,
            x43 * x80 * x87,
            x43 * x78 * x89,
            x59 * x83 * x87,
            x59 * x80 * x89,
            x59 * x78 * x91,
            x102 * x49 * x65,
            x102 * x43 * x82,
            x105 * x43 * x65,
            x102 * x59 * x68,
            x105 * x59 * x82,
            x109 * x59 * x65,
            x110 * x46 * x71,
            x111 * x41 * x71,
            x110 * x41 * x84,
            x11 * x113 * x71,
            x11 * x111 * x84,
            x11 * x110 * x75,
            x46 * x87 * x93,
            x41 * x87 * x96,
            x41 * x89 * x93,
            x100 * x11 * x87,
            x11 * x89 * x96,
            x11 * x91 * x93,
            x102 * x46 * x78,
            x102 * x41 * x80,
            x105 * x41 * x78,
            x102 * x11 * x83,
            x105 * x11 * x80,
            x109 * x11 * x78,
            x114 * x46 * x65,
            x114 * x41 * x82,
            x115 * x41 * x65,
            x11 * x114 * x68,
            x11 * x115 * x82,
            x11 * x117 * x65,
            x120 * x126 * x71,
            x130 * x133 * x71,
            x120 * x133 * x84,
            x139 * x141 * x71,
            x130 * x141 * x84,
            x120 * x141 * x75,
            x125 * x143 * x71,
            x132 * x145 * x71,
            x132 * x143 * x84,
            x140 * x153 * x71,
            x140 * x145 * x84,
            x140 * x143 * x75,
            x120 * x125 * x87,
            x130 * x132 * x87,
            x120 * x132 * x89,
            x139 * x140 * x87,
            x130 * x140 * x89,
            x120 * x140 * x91,
            x154 * x31 * x71,
            x122 * x159 * x71,
            x122 * x154 * x84,
            x165 * x24,
            x159 * x166 * x24,
            x154 * x69 * x75,
            x143 * x31 * x87,
            x122 * x145 * x87,
            x122 * x143 * x89,
            x153 * x164 * x24 * x85,
            x145 * x69 * x89,
            x143 * x69 * x91,
            x102 * x120 * x31,
            x102 * x122 * x130,
            x105 * x120 * x122,
            x102 * x139 * x69,
            x105 * x130 * x69,
            x109 * x120 * x69,
            x168 * x28 * x71,
            x169 * x170,
            x168 * x170 * x64,
            x164
            * (
                x0
                * (
                    2.0 * x151
                    + 2.0 * x152
                    + 2.0 * x157
                    + 2.0 * x158
                    + x30 * (x147 + 3.0 * x148 + x156 + x83)
                    + x98 * (x160 + x161)
                )
                + x163 * x76
            ),
            x166 * x169,
            x168 * x75 * x8,
            x154 * x28 * x87,
            x159 * x170 * x85,
            x14 * x154 * x89,
            x165 * x85,
            x159 * x8 * x89,
            x154 * x8 * x91,
            x102 * x143 * x28,
            x102 * x14 * x145,
            x105 * x14 * x143,
            x102 * x153 * x8,
            x105 * x145 * x8,
            x109 * x143 * x8,
            x114 * x120 * x28,
            x114 * x130 * x14,
            x115 * x120 * x14,
            x114 * x139 * x8,
            x115 * x130 * x8,
            x117 * x120 * x8,
            x126 * x173 * x65,
            x133 * x173 * x82,
            x133 * x177 * x65,
            x141 * x173 * x68,
            x141 * x177 * x82,
            x141 * x183 * x65,
            x125 * x173 * x78,
            x132 * x173 * x80,
            x132 * x177 * x78,
            x140 * x173 * x83,
            x140 * x177 * x80,
            x140 * x183 * x78,
            x125 * x185 * x65,
            x132 * x185 * x82,
            x132 * x187 * x65,
            x140 * x185 * x68,
            x140 * x187 * x82,
            x140 * x195 * x65,
            x173 * x31 * x93,
            x122 * x173 * x96,
            x122 * x177 * x93,
            x100 * x173 * x69,
            x177 * x69 * x96,
            x183 * x69 * x93,
            x185 * x31 * x78,
            x122 * x185 * x80,
            x122 * x187 * x78,
            x185 * x69 * x83,
            x187 * x69 * x80,
            x195 * x197 * x24 * x76,
            x198 * x31 * x65,
            x122 * x198 * x82,
            x122 * x203 * x65,
            x198 * x68 * x69,
            x203 * x204 * x24,
            x209 * x24,
            x110 * x173 * x28,
            x111 * x14 * x173,
            x110 * x14 * x177,
            x113 * x173 * x8,
            x111 * x177 * x8,
            x110 * x183 * x8,
            x185 * x28 * x93,
            x14 * x185 * x96,
            x14 * x187 * x93,
            x100 * x185 * x8,
            x187 * x8 * x96,
            x195 * x8 * x93,
            x198 * x28 * x78,
            x14 * x198 * x80,
            x203 * x210 * x76,
            x198 * x8 * x83,
            x203 * x8 * x80,
            x209 * x76,
            x212 * x28 * x65,
            x210 * x212 * x56,
            x210 * x213,
            x212 * x68 * x8,
            x204 * x213,
            x197
            * (
                x0
                * (
                    x107 * (x205 + x206)
                    + 2.0 * x193
                    + 2.0 * x194
                    + 2.0 * x201
                    + 2.0 * x202
                    + x30 * (x189 + 3.0 * x190 + x200 + x91)
                )
                + x208 * x85
            ),
        ]
    )


def diag_quadrupole3d_33(a, A, b, B, C):
    """Cartesian 3D (ff) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3 * x7
    x9 = -x2 - C[0]
    x10 = x7 * x9
    x11 = x0 * (x10 + x8)
    x12 = x0 * x7
    x13 = x8 * x9
    x14 = x12 + x13
    x15 = x14 * x3
    x16 = x11 + x15
    x17 = x16 * x3
    x18 = x7 * x9**2
    x19 = x12 + x18
    x20 = x19 * x3
    x21 = 2.0 * x0
    x22 = x10 * x21
    x23 = x20 + x22
    x24 = x23 * x3
    x25 = x3**2 * x7
    x26 = 3.0 * x12
    x27 = 2.0 * x13 + x26
    x28 = x0 * (x25 + x27)
    x29 = x0 * (x18 + x27)
    x30 = 2.0 * x28 + 3.0 * x29
    x31 = x0 * (2.0 * x17 + 3.0 * x24 + x30)
    x32 = -x2 - A[0]
    x33 = x16 * x32
    x34 = x0 * (3.0 * x25 + x26)
    x35 = x12 + x25
    x36 = x3 * x35
    x37 = x21 * x8
    x38 = x36 + x37
    x39 = x32 * x38
    x40 = x34 + x39
    x41 = 3.0 * x11
    x42 = x0 * (3.0 * x15 + x38 + x41)
    x43 = x32 * (x17 + x28)
    x44 = 2.0 * x32
    x45 = 4.0 * x0
    x46 = x10 * x45
    x47 = 2.0 * x11 + x46
    x48 = x0 * (2.0 * x15 + 2.0 * x20 + x47)
    x49 = x24 + x29
    x50 = x3 * x49
    x51 = x48 + x50
    x52 = x32 * x51
    x53 = x23 * x32
    x54 = 2.0 * x53
    x55 = x0 * (x24 + x30 + 2.0 * x33 + x54)
    x56 = x32 * x49
    x57 = x48 + x56
    x58 = x32 * x57
    x59 = x31 + x52
    x60 = x0 * (2.0 * x42 + 2.0 * x43 + 4.0 * x48 + x50 + 3.0 * x56) + x32 * x59
    x61 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x62 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x63 = 3.14159265358979 * x1 * x62
    x64 = x61 * x63
    x65 = -x1 * (a * A[1] + b * B[1])
    x66 = -x65 - B[1]
    x67 = x14 * x32
    x68 = 2.0 * x67
    x69 = x32 * x35
    x70 = x37 + x69
    x71 = x19 * x32
    x72 = x0 * (x20 + x47 + x68 + x71)
    x73 = x29 + x53
    x74 = x32 * x73
    x75 = x55 + x58
    x76 = x64 * (
        x0
        * (
            x21 * (x15 + x41 + x68 + x70)
            + x44 * (x28 + x33)
            + 2.0 * x48
            + 2.0 * x56
            + 2.0 * x72
            + 2.0 * x74
        )
        + x32 * x75
    )
    x77 = -x1 * (a * A[2] + b * B[2])
    x78 = -x77 - B[2]
    x79 = x6 * x61
    x80 = x66**2 * x79
    x81 = x0 * x79
    x82 = x80 + x81
    x83 = x32 * x8
    x84 = x10 * x32
    x85 = x22 + x71
    x86 = x0 * (x10 * x44 + x18 + x26) + x32 * x85
    x87 = x72 + x74
    x88 = (
        x0 * (x21 * (x13 + x26 + x83 + x84) + 2.0 * x29 + x44 * (x11 + x67) + x54 + x86)
        + x32 * x87
    )
    x89 = x6 * x62
    x90 = x64 * x78
    x91 = x78**2 * x89
    x92 = x0 * x89
    x93 = x91 + x92
    x94 = x66 * x82
    x95 = x66 * x79
    x96 = x21 * x95
    x97 = x94 + x96
    x98 = x32 * x7
    x99 = x0 * (x21 * (x10 + x98) + x44 * (x12 + x84) + x46 + 2.0 * x71) + x32 * x86
    x100 = x78 * x89
    x101 = x78 * x93
    x102 = x100 * x21
    x103 = x101 + x102
    x104 = -x65 - A[1]
    x105 = x60 * x64
    x106 = x104 * x79
    x107 = x106 * x66
    x108 = x107 + x81
    x109 = x104 * x82
    x110 = x109 + x96
    x111 = 3.0 * x81
    x112 = x0 * (x111 + 3.0 * x80)
    x113 = x104 * x97
    x114 = x112 + x113
    x115 = -x77 - A[2]
    x116 = x115 * x64
    x117 = x115 * x89
    x118 = x117 * x78
    x119 = x118 + x92
    x120 = x115 * x93
    x121 = x102 + x120
    x122 = 3.0 * x92
    x123 = x0 * (x122 + 3.0 * x91)
    x124 = x103 * x115
    x125 = x123 + x124
    x126 = x104**2 * x79
    x127 = x126 + x81
    x128 = x0 * (x106 + x95)
    x129 = x104 * x108
    x130 = x128 + x129
    x131 = 2.0 * x104
    x132 = x111 + x131 * x95
    x133 = x0 * (x132 + x80)
    x134 = x104 * x110
    x135 = x133 + x134
    x136 = x66 * x81
    x137 = x0 * (3.0 * x109 + 8.0 * x136 + x94) + x104 * x114
    x138 = x115**2 * x89
    x139 = x138 + x92
    x140 = x0 * (x100 + x117)
    x141 = x115 * x119
    x142 = x140 + x141
    x143 = 2.0 * x115
    x144 = x100 * x143 + x122
    x145 = x0 * (x144 + x91)
    x146 = x115 * x121
    x147 = x145 + x146
    x148 = x78 * x92
    x149 = x0 * (x101 + 3.0 * x120 + 8.0 * x148) + x115 * x125
    x150 = x104 * x127 + x106 * x21
    x151 = x0 * (x126 + x132) + x104 * x130
    x152 = x0 * (2.0 * x109 + 2.0 * x128 + 2.0 * x129 + 4.0 * x136) + x104 * x135
    x153 = x0 * (2.0 * x112 + 2.0 * x113 + 3.0 * x133 + 3.0 * x134) + x104 * x137
    x154 = x115 * x139 + x117 * x21
    x155 = x0 * (x138 + x144) + x115 * x142
    x156 = x0 * (2.0 * x120 + 2.0 * x140 + 2.0 * x141 + 4.0 * x148) + x115 * x147
    x157 = x0 * (2.0 * x123 + 2.0 * x124 + 3.0 * x145 + 3.0 * x146) + x115 * x149
    x158 = -x65 - C[1]
    x159 = x158**2 * x79
    x160 = x159 + x81
    x161 = x26 + x44 * x8
    x162 = x0 * (x161 + x25)
    x163 = x32 * x70
    x164 = x0 * (8.0 * x0 * x8 + x36 + 3.0 * x69) + x32 * x40
    x165 = x0 * (3.0 * x162 + 3.0 * x163 + 2.0 * x34 + 2.0 * x39) + x164 * x32
    x166 = x160 * x66
    x167 = x158 * x79
    x168 = x167 * x21
    x169 = x166 + x168
    x170 = x0 * (x8 + x98)
    x171 = x12 + x83
    x172 = x171 * x32
    x173 = x162 + x163
    x174 = x0 * (2.0 * x170 + 2.0 * x172 + x45 * x8 + 2.0 * x69) + x173 * x32
    x175 = x158 * x95
    x176 = x111 + 2.0 * x175
    x177 = x0 * (x159 + x176)
    x178 = x169 * x66
    x179 = x177 + x178
    x180 = x32**2 * x7
    x181 = x170 + x172
    x182 = x0 * (x161 + x180) + x181 * x32
    x183 = x175 + x81
    x184 = x183 * x66
    x185 = x0 * (x167 + x95)
    x186 = 4.0 * x158 * x81
    x187 = 2.0 * x185 + x186
    x188 = x0 * (2.0 * x166 + 2.0 * x184 + x187)
    x189 = x179 * x66
    x190 = x188 + x189
    x191 = x12 + x180
    x192 = x191 * x32 + x21 * x98
    x193 = x104 * x160
    x194 = x168 + x193
    x195 = x104 * x169
    x196 = x177 + x195
    x197 = x104 * x179
    x198 = x188 + x197
    x199 = x184 + x185
    x200 = x199 * x66
    x201 = x0 * (x176 + x80)
    x202 = 3.0 * x177 + 2.0 * x201
    x203 = x0 * (3.0 * x178 + 2.0 * x200 + x202)
    x204 = x104 * x190
    x205 = x203 + x204
    x206 = x0 * (x111 + x131 * x167 + x159) + x104 * x194
    x207 = x104 * x183
    x208 = 2.0 * x207
    x209 = x0 * (x166 + x187 + x193 + x208)
    x210 = x104 * x196
    x211 = x209 + x210
    x212 = x104 * x199
    x213 = 2.0 * x195
    x214 = x0 * (x178 + x202 + 2.0 * x212 + x213)
    x215 = x104 * x198
    x216 = x214 + x215
    x217 = 3.0 * x185
    x218 = x0 * (3.0 * x184 + x217 + x97)
    x219 = x104 * (x200 + x201)
    x220 = x0 * (4.0 * x188 + x189 + 3.0 * x197 + 2.0 * x218 + 2.0 * x219) + x104 * x205
    x221 = x5 * x63
    x222 = x220 * x221
    x223 = x221 * x32
    x224 = x106 * x158
    x225 = (
        x0 * (x131 * (x224 + x81) + x186 + 2.0 * x193 + x21 * (x106 + x167)) + x104 * x206
    )
    x226 = (
        x0
        * (
            x131 * (x185 + x207)
            + 2.0 * x177
            + x206
            + x21 * (x107 + x111 + x175 + x224)
            + x213
        )
        + x104 * x211
    )
    x227 = x221 * (
        x0
        * (
            x131 * (x201 + x212)
            + 2.0 * x188
            + 2.0 * x197
            + 2.0 * x209
            + x21 * (x110 + x184 + x208 + x217)
            + 2.0 * x210
        )
        + x104 * x216
    )
    x228 = x221 * x3
    x229 = -x77 - C[2]
    x230 = x229**2 * x89
    x231 = x230 + x92
    x232 = x231 * x78
    x233 = x229 * x89
    x234 = x21 * x233
    x235 = x232 + x234
    x236 = x100 * x229
    x237 = x122 + 2.0 * x236
    x238 = x0 * (x230 + x237)
    x239 = x235 * x78
    x240 = x238 + x239
    x241 = x236 + x92
    x242 = x241 * x78
    x243 = x0 * (x100 + x233)
    x244 = 4.0 * x229 * x92
    x245 = 2.0 * x243 + x244
    x246 = x0 * (2.0 * x232 + 2.0 * x242 + x245)
    x247 = x240 * x78
    x248 = x246 + x247
    x249 = x115 * x231
    x250 = x234 + x249
    x251 = x115 * x235
    x252 = x238 + x251
    x253 = x115 * x240
    x254 = x246 + x253
    x255 = x242 + x243
    x256 = x255 * x78
    x257 = x0 * (x237 + x91)
    x258 = 3.0 * x238 + 2.0 * x257
    x259 = x0 * (3.0 * x239 + 2.0 * x256 + x258)
    x260 = x115 * x248
    x261 = x259 + x260
    x262 = 3.14159265358979 * x1 * x5 * x61
    x263 = x262 * x32
    x264 = x0 * (x122 + x143 * x233 + x230) + x115 * x250
    x265 = x115 * x241
    x266 = 2.0 * x265
    x267 = x0 * (x232 + x245 + x249 + x266)
    x268 = x115 * x252
    x269 = x267 + x268
    x270 = x115 * x255
    x271 = 2.0 * x251
    x272 = x0 * (x239 + x258 + 2.0 * x270 + x271)
    x273 = x115 * x254
    x274 = x272 + x273
    x275 = 3.0 * x243
    x276 = x0 * (x103 + 3.0 * x242 + x275)
    x277 = x115 * (x256 + x257)
    x278 = x0 * (4.0 * x246 + x247 + 3.0 * x253 + 2.0 * x276 + 2.0 * x277) + x115 * x261
    x279 = x262 * x278
    x280 = x262 * x3
    x281 = x117 * x229
    x282 = (
        x0 * (x143 * (x281 + x92) + x21 * (x117 + x233) + x244 + 2.0 * x249) + x115 * x264
    )
    x283 = (
        x0
        * (
            x143 * (x243 + x265)
            + x21 * (x118 + x122 + x236 + x281)
            + 2.0 * x238
            + x264
            + x271
        )
        + x115 * x269
    )
    x284 = x262 * (
        x0
        * (
            x143 * (x257 + x270)
            + x21 * (x121 + x242 + x266 + x275)
            + 2.0 * x246
            + 2.0 * x253
            + 2.0 * x267
            + 2.0 * x268
        )
        + x115 * x274
    )

    # 300 item(s)
    return numpy.array(
        [
            x64
            * (
                x0
                * (
                    x21 * (x17 + 4.0 * x28 + 3.0 * x33 + x40)
                    + 2.0 * x31
                    + x44 * (x42 + x43)
                    + 2.0 * x52
                    + 3.0 * x55
                    + 3.0 * x58
                )
                + x32 * x60
            ),
            x66 * x76,
            x76 * x78,
            x82 * x88 * x89,
            x66 * x88 * x90,
            x79 * x88 * x93,
            x89 * x97 * x99,
            x100 * x82 * x99,
            x93 * x95 * x99,
            x103 * x79 * x99,
            x104 * x105,
            x108 * x75 * x89,
            x104 * x75 * x90,
            x110 * x87 * x89,
            x100 * x108 * x87,
            x106 * x87 * x93,
            x114 * x86 * x89,
            x100 * x110 * x86,
            x108 * x86 * x93,
            x103 * x106 * x86,
            x105 * x115,
            x116 * x66 * x75,
            x119 * x75 * x79,
            x117 * x82 * x87,
            x119 * x87 * x95,
            x121 * x79 * x87,
            x117 * x86 * x97,
            x119 * x82 * x86,
            x121 * x86 * x95,
            x125 * x79 * x86,
            x127 * x59 * x89,
            x130 * x57 * x89,
            x100 * x127 * x57,
            x135 * x73 * x89,
            x100 * x130 * x73,
            x127 * x73 * x93,
            x137 * x85 * x89,
            x100 * x135 * x85,
            x130 * x85 * x93,
            x103 * x127 * x85,
            x104 * x116 * x59,
            x108 * x117 * x57,
            x106 * x119 * x57,
            x110 * x117 * x73,
            x108 * x119 * x73,
            x106 * x121 * x73,
            x114 * x117 * x85,
            x110 * x119 * x85,
            x108 * x121 * x85,
            x106 * x125 * x85,
            x139 * x59 * x79,
            x139 * x57 * x95,
            x142 * x57 * x79,
            x139 * x73 * x82,
            x142 * x73 * x95,
            x147 * x73 * x79,
            x139 * x85 * x97,
            x142 * x82 * x85,
            x147 * x85 * x95,
            x149 * x79 * x85,
            x150 * x51 * x89,
            x151 * x49 * x89,
            x100 * x150 * x49,
            x152 * x23 * x89,
            x100 * x151 * x23,
            x150 * x23 * x93,
            x153 * x19 * x89,
            x100 * x152 * x19,
            x151 * x19 * x93,
            x103 * x150 * x19,
            x117 * x127 * x51,
            x117 * x130 * x49,
            x119 * x127 * x49,
            x117 * x135 * x23,
            x119 * x130 * x23,
            x121 * x127 * x23,
            x117 * x137 * x19,
            x119 * x135 * x19,
            x121 * x130 * x19,
            x125 * x127 * x19,
            x106 * x139 * x51,
            x108 * x139 * x49,
            x106 * x142 * x49,
            x110 * x139 * x23,
            x108 * x142 * x23,
            x106 * x147 * x23,
            x114 * x139 * x19,
            x110 * x142 * x19,
            x108 * x147 * x19,
            x106 * x149 * x19,
            x154 * x51 * x79,
            x154 * x49 * x95,
            x155 * x49 * x79,
            x154 * x23 * x82,
            x155 * x23 * x95,
            x156 * x23 * x79,
            x154 * x19 * x97,
            x155 * x19 * x82,
            x156 * x19 * x95,
            x157 * x19 * x79,
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
            x164 * x194 * x89,
            x173 * x196 * x89,
            x100 * x173 * x194,
            x181 * x198 * x89,
            x100 * x181 * x196,
            x181 * x194 * x93,
            x191 * x205 * x89,
            x100 * x191 * x198,
            x191 * x196 * x93,
            x103 * x191 * x194,
            x117 * x160 * x164,
            x117 * x169 * x173,
            x119 * x160 * x173,
            x117 * x179 * x181,
            x119 * x169 * x181,
            x121 * x160 * x181,
            x117 * x190 * x191,
            x119 * x179 * x191,
            x121 * x169 * x191,
            x125 * x160 * x191,
            x206 * x40 * x89,
            x211 * x70 * x89,
            x100 * x206 * x70,
            x171 * x216 * x89,
            x100 * x171 * x211,
            x171 * x206 * x93,
            x222 * x32,
            x216 * x223 * x78,
            x211 * x93 * x98,
            x103 * x206 * x98,
            x117 * x194 * x40,
            x117 * x196 * x70,
            x119 * x194 * x70,
            x117 * x171 * x198,
            x119 * x171 * x196,
            x121 * x171 * x194,
            x115 * x205 * x223,
            x119 * x198 * x98,
            x121 * x196 * x98,
            x125 * x194 * x98,
            x139 * x160 * x40,
            x139 * x169 * x70,
            x142 * x160 * x70,
            x139 * x171 * x179,
            x142 * x169 * x171,
            x147 * x160 * x171,
            x139 * x190 * x98,
            x142 * x179 * x98,
            x147 * x169 * x98,
            x149 * x160 * x98,
            x225 * x38 * x89,
            x226 * x35 * x89,
            x100 * x225 * x35,
            x227 * x3,
            x226 * x228 * x78,
            x225 * x8 * x93,
            x221
            * (
                x0
                * (
                    x131 * (x218 + x219)
                    + 2.0 * x203
                    + 2.0 * x204
                    + x21 * (x114 + x200 + 4.0 * x201 + 3.0 * x212)
                    + 3.0 * x214
                    + 3.0 * x215
                )
                + x104 * x220
            ),
            x227 * x78,
            x226 * x7 * x93,
            x103 * x225 * x7,
            x117 * x206 * x38,
            x117 * x211 * x35,
            x119 * x206 * x35,
            x115 * x216 * x228,
            x119 * x211 * x8,
            x121 * x206 * x8,
            x115 * x222,
            x119 * x216 * x7,
            x121 * x211 * x7,
            x125 * x206 * x7,
            x139 * x194 * x38,
            x139 * x196 * x35,
            x142 * x194 * x35,
            x139 * x198 * x8,
            x142 * x196 * x8,
            x147 * x194 * x8,
            x139 * x205 * x7,
            x142 * x198 * x7,
            x147 * x196 * x7,
            x149 * x194 * x7,
            x154 * x160 * x38,
            x154 * x169 * x35,
            x155 * x160 * x35,
            x154 * x179 * x8,
            x155 * x169 * x8,
            x156 * x160 * x8,
            x154 * x190 * x7,
            x155 * x179 * x7,
            x156 * x169 * x7,
            x157 * x160 * x7,
            x165 * x231 * x79,
            x174 * x231 * x95,
            x174 * x235 * x79,
            x182 * x231 * x82,
            x182 * x235 * x95,
            x182 * x240 * x79,
            x192 * x231 * x97,
            x192 * x235 * x82,
            x192 * x240 * x95,
            x192 * x248 * x79,
            x106 * x164 * x231,
            x108 * x173 * x231,
            x106 * x173 * x235,
            x110 * x181 * x231,
            x108 * x181 * x235,
            x106 * x181 * x240,
            x114 * x191 * x231,
            x110 * x191 * x235,
            x108 * x191 * x240,
            x106 * x191 * x248,
            x164 * x250 * x79,
            x173 * x250 * x95,
            x173 * x252 * x79,
            x181 * x250 * x82,
            x181 * x252 * x95,
            x181 * x254 * x79,
            x191 * x250 * x97,
            x191 * x252 * x82,
            x191 * x254 * x95,
            x191 * x261 * x79,
            x127 * x231 * x40,
            x130 * x231 * x70,
            x127 * x235 * x70,
            x135 * x171 * x231,
            x130 * x171 * x235,
            x127 * x171 * x240,
            x137 * x231 * x98,
            x135 * x235 * x98,
            x130 * x240 * x98,
            x127 * x248 * x98,
            x106 * x250 * x40,
            x108 * x250 * x70,
            x106 * x252 * x70,
            x110 * x171 * x250,
            x108 * x171 * x252,
            x106 * x171 * x254,
            x114 * x250 * x98,
            x110 * x252 * x98,
            x108 * x254 * x98,
            x104 * x261 * x263,
            x264 * x40 * x79,
            x264 * x70 * x95,
            x269 * x70 * x79,
            x171 * x264 * x82,
            x171 * x269 * x95,
            x171 * x274 * x79,
            x264 * x97 * x98,
            x269 * x82 * x98,
            x263 * x274 * x66,
            x279 * x32,
            x150 * x231 * x38,
            x151 * x231 * x35,
            x150 * x235 * x35,
            x152 * x231 * x8,
            x151 * x235 * x8,
            x150 * x240 * x8,
            x153 * x231 * x7,
            x152 * x235 * x7,
            x151 * x240 * x7,
            x150 * x248 * x7,
            x127 * x250 * x38,
            x130 * x250 * x35,
            x127 * x252 * x35,
            x135 * x250 * x8,
            x130 * x252 * x8,
            x127 * x254 * x8,
            x137 * x250 * x7,
            x135 * x252 * x7,
            x130 * x254 * x7,
            x127 * x261 * x7,
            x106 * x264 * x38,
            x108 * x264 * x35,
            x106 * x269 * x35,
            x110 * x264 * x8,
            x108 * x269 * x8,
            x104 * x274 * x280,
            x114 * x264 * x7,
            x110 * x269 * x7,
            x108 * x274 * x7,
            x104 * x279,
            x282 * x38 * x79,
            x282 * x35 * x95,
            x283 * x35 * x79,
            x282 * x8 * x82,
            x280 * x283 * x66,
            x284 * x3,
            x282 * x7 * x97,
            x283 * x7 * x82,
            x284 * x66,
            x262
            * (
                x0
                * (
                    x143 * (x276 + x277)
                    + x21 * (x125 + x256 + 4.0 * x257 + 3.0 * x270)
                    + 2.0 * x259
                    + 2.0 * x260
                    + 3.0 * x272
                    + 3.0 * x273
                )
                + x115 * x278
            ),
        ]
    )


def diag_quadrupole3d_34(a, A, b, B, C):
    """Cartesian 3D (fg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3**2 * x7
    x9 = x0 * x7
    x10 = 3.0 * x9
    x11 = -x2 - C[0]
    x12 = x3 * x7
    x13 = x11 * x12
    x14 = x10 + 2.0 * x13
    x15 = x0 * (x14 + x8)
    x16 = x11 * x7
    x17 = x0 * (x12 + x16)
    x18 = x13 + x9
    x19 = x18 * x3
    x20 = x17 + x19
    x21 = x20 * x3
    x22 = x15 + x21
    x23 = x22 * x3
    x24 = x11**2 * x7
    x25 = x0 * (x14 + x24)
    x26 = x24 + x9
    x27 = x26 * x3
    x28 = 2.0 * x0
    x29 = x16 * x28
    x30 = x27 + x29
    x31 = x3 * x30
    x32 = x25 + x31
    x33 = x3 * x32
    x34 = 3.0 * x17
    x35 = x8 + x9
    x36 = x3 * x35
    x37 = x12 * x28
    x38 = x36 + x37
    x39 = x0 * (3.0 * x19 + x34 + x38)
    x40 = 4.0 * x9
    x41 = x11 * x40
    x42 = 2.0 * x17 + x41
    x43 = x0 * (2.0 * x19 + 2.0 * x27 + x42)
    x44 = 2.0 * x39 + 4.0 * x43
    x45 = x0 * (2.0 * x23 + 4.0 * x33 + x44)
    x46 = -x2 - A[0]
    x47 = x22 * x46
    x48 = 8.0 * x3 * x9
    x49 = x0 * (4.0 * x36 + x48)
    x50 = x0 * (x10 + 3.0 * x8)
    x51 = x3 * x38
    x52 = x50 + x51
    x53 = x46 * x52
    x54 = x49 + x53
    x55 = 4.0 * x15
    x56 = x0 * (4.0 * x21 + x52 + x55)
    x57 = x46 * (x23 + x39)
    x58 = 2.0 * x46
    x59 = 2.0 * x15 + 3.0 * x25
    x60 = x0 * (2.0 * x21 + 3.0 * x31 + x59)
    x61 = x33 + x43
    x62 = x3 * x61
    x63 = x60 + x62
    x64 = x46 * x63
    x65 = x32 * x46
    x66 = x0 * (x33 + x44 + 2.0 * x47 + 3.0 * x65)
    x67 = x46 * x61
    x68 = x60 + x67
    x69 = x46 * x68
    x70 = x45 + x64
    x71 = x0 * (2.0 * x56 + 2.0 * x57 + 5.0 * x60 + x62 + 4.0 * x67) + x46 * x70
    x72 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x73 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x74 = 3.14159265358979 * x1 * x73
    x75 = x72 * x74
    x76 = -x1 * (a * A[1] + b * B[1])
    x77 = -x76 - B[1]
    x78 = x20 * x46
    x79 = x38 * x46
    x80 = x50 + x79
    x81 = x30 * x46
    x82 = 2.0 * x81
    x83 = x0 * (x31 + x59 + 2.0 * x78 + x82)
    x84 = x43 + x65
    x85 = x46 * x84
    x86 = x66 + x69
    x87 = x75 * (
        x0
        * (
            x28 * (x21 + x55 + 3.0 * x78 + x80)
            + x58 * (x39 + x47)
            + 2.0 * x60
            + 2.0 * x67
            + 3.0 * x83
            + 3.0 * x85
        )
        + x46 * x86
    )
    x88 = -x1 * (a * A[2] + b * B[2])
    x89 = -x88 - B[2]
    x90 = x6 * x72
    x91 = x77**2 * x90
    x92 = x0 * x90
    x93 = x91 + x92
    x94 = x18 * x46
    x95 = 2.0 * x94
    x96 = x35 * x46
    x97 = x37 + x96
    x98 = x26 * x46
    x99 = x0 * (x27 + x42 + x95 + x98)
    x100 = x25 + x81
    x101 = x100 * x46
    x102 = x83 + x85
    x103 = (
        x0
        * (
            2.0 * x101
            + x28 * (x19 + x34 + x95 + x97)
            + 2.0 * x43
            + x58 * (x15 + x78)
            + 2.0 * x65
            + 2.0 * x99
        )
        + x102 * x46
    )
    x104 = x6 * x73
    x105 = x75 * x89
    x106 = x104 * x89**2
    x107 = x0 * x104
    x108 = x106 + x107
    x109 = x77 * x93
    x110 = x77 * x90
    x111 = x110 * x28
    x112 = x109 + x111
    x113 = x12 * x46
    x114 = x16 * x46
    x115 = x29 + x98
    x116 = x0 * (x10 + x16 * x58 + x24) + x115 * x46
    x117 = x101 + x99
    x118 = (
        x0
        * (x116 + 2.0 * x25 + x28 * (x10 + x113 + x114 + x13) + x58 * (x17 + x94) + x82)
        + x117 * x46
    )
    x119 = x104 * x89
    x120 = x108 * x89
    x121 = x119 * x28
    x122 = x120 + x121
    x123 = 3.0 * x92
    x124 = x0 * (x123 + 3.0 * x91)
    x125 = x112 * x77
    x126 = x124 + x125
    x127 = x46 * x7
    x128 = x0 * (x28 * (x127 + x16) + x41 + x58 * (x114 + x9) + 2.0 * x98) + x116 * x46
    x129 = 3.0 * x107
    x130 = x0 * (3.0 * x106 + x129)
    x131 = x122 * x89
    x132 = x130 + x131
    x133 = -x76 - A[1]
    x134 = x71 * x75
    x135 = x133 * x90
    x136 = x135 * x77
    x137 = x136 + x92
    x138 = x133 * x93
    x139 = x111 + x138
    x140 = x112 * x133
    x141 = x124 + x140
    x142 = x77 * x92
    x143 = 8.0 * x142
    x144 = x0 * (4.0 * x109 + x143)
    x145 = x126 * x133
    x146 = x144 + x145
    x147 = -x88 - A[2]
    x148 = x147 * x75
    x149 = x104 * x147
    x150 = x149 * x89
    x151 = x107 + x150
    x152 = x108 * x147
    x153 = x121 + x152
    x154 = x122 * x147
    x155 = x130 + x154
    x156 = x107 * x89
    x157 = 8.0 * x156
    x158 = x0 * (4.0 * x120 + x157)
    x159 = x132 * x147
    x160 = x158 + x159
    x161 = x133**2 * x90
    x162 = x161 + x92
    x163 = x0 * (x110 + x135)
    x164 = x133 * x137
    x165 = x163 + x164
    x166 = 2.0 * x133
    x167 = x110 * x166 + x123
    x168 = x0 * (x167 + x91)
    x169 = x133 * x139
    x170 = x168 + x169
    x171 = x0 * (x109 + 3.0 * x138 + x143)
    x172 = x133 * x141
    x173 = x171 + x172
    x174 = x0 * (5.0 * x124 + x125 + 4.0 * x140) + x133 * x146
    x175 = x104 * x147**2
    x176 = x107 + x175
    x177 = x0 * (x119 + x149)
    x178 = x147 * x151
    x179 = x177 + x178
    x180 = 2.0 * x147
    x181 = x119 * x180 + x129
    x182 = x0 * (x106 + x181)
    x183 = x147 * x153
    x184 = x182 + x183
    x185 = x0 * (x120 + 3.0 * x152 + x157)
    x186 = x147 * x155
    x187 = x185 + x186
    x188 = x0 * (5.0 * x130 + x131 + 4.0 * x154) + x147 * x160
    x189 = x133 * x162 + x135 * x28
    x190 = x0 * (x161 + x167) + x133 * x165
    x191 = x0 * (2.0 * x138 + 4.0 * x142 + 2.0 * x163 + 2.0 * x164) + x133 * x170
    x192 = x0 * (2.0 * x124 + 2.0 * x140 + 3.0 * x168 + 3.0 * x169) + x133 * x173
    x193 = x0 * (2.0 * x144 + 2.0 * x145 + 4.0 * x171 + 4.0 * x172) + x133 * x174
    x194 = x147 * x176 + x149 * x28
    x195 = x0 * (x175 + x181) + x147 * x179
    x196 = x0 * (2.0 * x152 + 4.0 * x156 + 2.0 * x177 + 2.0 * x178) + x147 * x184
    x197 = x0 * (2.0 * x130 + 2.0 * x154 + 3.0 * x182 + 3.0 * x183) + x147 * x187
    x198 = x0 * (2.0 * x158 + 2.0 * x159 + 4.0 * x185 + 4.0 * x186) + x147 * x188
    x199 = -x76 - C[1]
    x200 = x199**2 * x90
    x201 = x200 + x92
    x202 = x0 * (x36 + x48 + 3.0 * x96)
    x203 = x46 * x80
    x204 = x0 * (5.0 * x50 + x51 + 4.0 * x79) + x46 * x54
    x205 = x0 * (4.0 * x202 + 4.0 * x203 + 2.0 * x49 + 2.0 * x53) + x204 * x46
    x206 = x201 * x77
    x207 = x199 * x90
    x208 = x207 * x28
    x209 = x206 + x208
    x210 = x10 + x12 * x58
    x211 = x0 * (x210 + x8)
    x212 = x46 * x97
    x213 = x202 + x203
    x214 = x0 * (3.0 * x211 + 3.0 * x212 + 2.0 * x50 + 2.0 * x79) + x213 * x46
    x215 = x110 * x199
    x216 = x123 + 2.0 * x215
    x217 = x0 * (x200 + x216)
    x218 = x209 * x77
    x219 = x217 + x218
    x220 = x0 * (x12 + x127)
    x221 = x113 + x9
    x222 = x221 * x46
    x223 = x211 + x212
    x224 = x0 * (2.0 * x220 + 2.0 * x222 + x3 * x40 + 2.0 * x96) + x223 * x46
    x225 = x215 + x92
    x226 = x225 * x77
    x227 = x0 * (x110 + x207)
    x228 = 4.0 * x199 * x92
    x229 = 2.0 * x227 + x228
    x230 = x0 * (2.0 * x206 + 2.0 * x226 + x229)
    x231 = x219 * x77
    x232 = x230 + x231
    x233 = x46**2 * x7
    x234 = x220 + x222
    x235 = x0 * (x210 + x233) + x234 * x46
    x236 = x226 + x227
    x237 = x236 * x77
    x238 = x0 * (x216 + x91)
    x239 = 3.0 * x217 + 2.0 * x238
    x240 = x0 * (3.0 * x218 + 2.0 * x237 + x239)
    x241 = x232 * x77
    x242 = x240 + x241
    x243 = x233 + x9
    x244 = x127 * x28 + x243 * x46
    x245 = x133 * x201
    x246 = x208 + x245
    x247 = x133 * x209
    x248 = x217 + x247
    x249 = x133 * x219
    x250 = x230 + x249
    x251 = x133 * x232
    x252 = x240 + x251
    x253 = x237 + x238
    x254 = x253 * x77
    x255 = 3.0 * x227
    x256 = x0 * (x112 + 3.0 * x226 + x255)
    x257 = 4.0 * x230 + 2.0 * x256
    x258 = x0 * (4.0 * x231 + 2.0 * x254 + x257)
    x259 = x133 * x242
    x260 = x258 + x259
    x261 = x0 * (x123 + x166 * x207 + x200) + x133 * x246
    x262 = x133 * x225
    x263 = 2.0 * x262
    x264 = x0 * (x206 + x229 + x245 + x263)
    x265 = x133 * x248
    x266 = x264 + x265
    x267 = x133 * x236
    x268 = 2.0 * x247
    x269 = x0 * (x218 + x239 + 2.0 * x267 + x268)
    x270 = x133 * x250
    x271 = x269 + x270
    x272 = x133 * x253
    x273 = x0 * (x231 + 3.0 * x249 + x257 + 2.0 * x272)
    x274 = x133 * x252
    x275 = x273 + x274
    x276 = 4.0 * x238
    x277 = x0 * (x126 + 4.0 * x237 + x276)
    x278 = x133 * (x254 + x256)
    x279 = x0 * (5.0 * x240 + x241 + 4.0 * x251 + 2.0 * x277 + 2.0 * x278) + x133 * x260
    x280 = x5 * x74
    x281 = x279 * x280
    x282 = x280 * x46
    x283 = x135 * x199
    x284 = (
        x0 * (x166 * (x283 + x92) + x228 + 2.0 * x245 + x28 * (x135 + x207)) + x133 * x261
    )
    x285 = (
        x0
        * (
            x166 * (x227 + x262)
            + 2.0 * x217
            + x261
            + x268
            + x28 * (x123 + x136 + x215 + x283)
        )
        + x133 * x266
    )
    x286 = (
        x0
        * (
            x166 * (x238 + x267)
            + 2.0 * x230
            + 2.0 * x249
            + 2.0 * x264
            + 2.0 * x265
            + x28 * (x139 + x226 + x255 + x263)
        )
        + x133 * x271
    )
    x287 = x280 * (
        x0
        * (
            x166 * (x256 + x272)
            + 2.0 * x240
            + 2.0 * x251
            + 3.0 * x269
            + 3.0 * x270
            + x28 * (x141 + x237 + 3.0 * x267 + x276)
        )
        + x133 * x275
    )
    x288 = x280 * x3
    x289 = -x88 - C[2]
    x290 = x104 * x289**2
    x291 = x107 + x290
    x292 = x291 * x89
    x293 = x104 * x289
    x294 = x28 * x293
    x295 = x292 + x294
    x296 = x119 * x289
    x297 = x129 + 2.0 * x296
    x298 = x0 * (x290 + x297)
    x299 = x295 * x89
    x300 = x298 + x299
    x301 = x107 + x296
    x302 = x301 * x89
    x303 = x0 * (x119 + x293)
    x304 = 4.0 * x107 * x289
    x305 = 2.0 * x303 + x304
    x306 = x0 * (2.0 * x292 + 2.0 * x302 + x305)
    x307 = x300 * x89
    x308 = x306 + x307
    x309 = x302 + x303
    x310 = x309 * x89
    x311 = x0 * (x106 + x297)
    x312 = 3.0 * x298 + 2.0 * x311
    x313 = x0 * (3.0 * x299 + 2.0 * x310 + x312)
    x314 = x308 * x89
    x315 = x313 + x314
    x316 = x147 * x291
    x317 = x294 + x316
    x318 = x147 * x295
    x319 = x298 + x318
    x320 = x147 * x300
    x321 = x306 + x320
    x322 = x147 * x308
    x323 = x313 + x322
    x324 = x310 + x311
    x325 = x324 * x89
    x326 = 3.0 * x303
    x327 = x0 * (x122 + 3.0 * x302 + x326)
    x328 = 4.0 * x306 + 2.0 * x327
    x329 = x0 * (4.0 * x307 + 2.0 * x325 + x328)
    x330 = x147 * x315
    x331 = x329 + x330
    x332 = 3.14159265358979 * x1 * x5 * x72
    x333 = x332 * x46
    x334 = x0 * (x129 + x180 * x293 + x290) + x147 * x317
    x335 = x147 * x301
    x336 = 2.0 * x335
    x337 = x0 * (x292 + x305 + x316 + x336)
    x338 = x147 * x319
    x339 = x337 + x338
    x340 = x147 * x309
    x341 = 2.0 * x318
    x342 = x0 * (x299 + x312 + 2.0 * x340 + x341)
    x343 = x147 * x321
    x344 = x342 + x343
    x345 = x147 * x324
    x346 = x0 * (x307 + 3.0 * x320 + x328 + 2.0 * x345)
    x347 = x147 * x323
    x348 = x346 + x347
    x349 = 4.0 * x311
    x350 = x0 * (x132 + 4.0 * x310 + x349)
    x351 = x147 * (x325 + x327)
    x352 = x0 * (5.0 * x313 + x314 + 4.0 * x322 + 2.0 * x350 + 2.0 * x351) + x147 * x331
    x353 = x332 * x352
    x354 = x3 * x332
    x355 = x149 * x289
    x356 = (
        x0 * (x180 * (x107 + x355) + x28 * (x149 + x293) + x304 + 2.0 * x316)
        + x147 * x334
    )
    x357 = (
        x0
        * (
            x180 * (x303 + x335)
            + x28 * (x129 + x150 + x296 + x355)
            + 2.0 * x298
            + x334
            + x341
        )
        + x147 * x339
    )
    x358 = (
        x0
        * (
            x180 * (x311 + x340)
            + x28 * (x153 + x302 + x326 + x336)
            + 2.0 * x306
            + 2.0 * x320
            + 2.0 * x337
            + 2.0 * x338
        )
        + x147 * x344
    )
    x359 = x332 * (
        x0
        * (
            x180 * (x327 + x345)
            + x28 * (x155 + x310 + 3.0 * x340 + x349)
            + 2.0 * x313
            + 2.0 * x322
            + 3.0 * x342
            + 3.0 * x343
        )
        + x147 * x348
    )

    # 450 item(s)
    return numpy.array(
        [
            x75
            * (
                x0
                * (
                    x28 * (x23 + 5.0 * x39 + 4.0 * x47 + x54)
                    + 2.0 * x45
                    + x58 * (x56 + x57)
                    + 2.0 * x64
                    + 4.0 * x66
                    + 4.0 * x69
                )
                + x46 * x71
            ),
            x77 * x87,
            x87 * x89,
            x103 * x104 * x93,
            x103 * x105 * x77,
            x103 * x108 * x90,
            x104 * x112 * x118,
            x118 * x119 * x93,
            x108 * x110 * x118,
            x118 * x122 * x90,
            x104 * x126 * x128,
            x112 * x119 * x128,
            x108 * x128 * x93,
            x110 * x122 * x128,
            x128 * x132 * x90,
            x133 * x134,
            x104 * x137 * x86,
            x105 * x133 * x86,
            x102 * x104 * x139,
            x102 * x119 * x137,
            x102 * x108 * x135,
            x104 * x117 * x141,
            x117 * x119 * x139,
            x108 * x117 * x137,
            x117 * x122 * x135,
            x104 * x116 * x146,
            x116 * x119 * x141,
            x108 * x116 * x139,
            x116 * x122 * x137,
            x116 * x132 * x135,
            x134 * x147,
            x148 * x77 * x86,
            x151 * x86 * x90,
            x102 * x149 * x93,
            x102 * x110 * x151,
            x102 * x153 * x90,
            x112 * x117 * x149,
            x117 * x151 * x93,
            x110 * x117 * x153,
            x117 * x155 * x90,
            x116 * x126 * x149,
            x112 * x116 * x151,
            x116 * x153 * x93,
            x110 * x116 * x155,
            x116 * x160 * x90,
            x104 * x162 * x70,
            x104 * x165 * x68,
            x119 * x162 * x68,
            x104 * x170 * x84,
            x119 * x165 * x84,
            x108 * x162 * x84,
            x100 * x104 * x173,
            x100 * x119 * x170,
            x100 * x108 * x165,
            x100 * x122 * x162,
            x104 * x115 * x174,
            x115 * x119 * x173,
            x108 * x115 * x170,
            x115 * x122 * x165,
            x115 * x132 * x162,
            x133 * x148 * x70,
            x137 * x149 * x68,
            x135 * x151 * x68,
            x139 * x149 * x84,
            x137 * x151 * x84,
            x135 * x153 * x84,
            x100 * x141 * x149,
            x100 * x139 * x151,
            x100 * x137 * x153,
            x100 * x135 * x155,
            x115 * x146 * x149,
            x115 * x141 * x151,
            x115 * x139 * x153,
            x115 * x137 * x155,
            x115 * x135 * x160,
            x176 * x70 * x90,
            x110 * x176 * x68,
            x179 * x68 * x90,
            x176 * x84 * x93,
            x110 * x179 * x84,
            x184 * x84 * x90,
            x100 * x112 * x176,
            x100 * x179 * x93,
            x100 * x110 * x184,
            x100 * x187 * x90,
            x115 * x126 * x176,
            x112 * x115 * x179,
            x115 * x184 * x93,
            x110 * x115 * x187,
            x115 * x188 * x90,
            x104 * x189 * x63,
            x104 * x190 * x61,
            x119 * x189 * x61,
            x104 * x191 * x32,
            x119 * x190 * x32,
            x108 * x189 * x32,
            x104 * x192 * x30,
            x119 * x191 * x30,
            x108 * x190 * x30,
            x122 * x189 * x30,
            x104 * x193 * x26,
            x119 * x192 * x26,
            x108 * x191 * x26,
            x122 * x190 * x26,
            x132 * x189 * x26,
            x149 * x162 * x63,
            x149 * x165 * x61,
            x151 * x162 * x61,
            x149 * x170 * x32,
            x151 * x165 * x32,
            x153 * x162 * x32,
            x149 * x173 * x30,
            x151 * x170 * x30,
            x153 * x165 * x30,
            x155 * x162 * x30,
            x149 * x174 * x26,
            x151 * x173 * x26,
            x153 * x170 * x26,
            x155 * x165 * x26,
            x160 * x162 * x26,
            x135 * x176 * x63,
            x137 * x176 * x61,
            x135 * x179 * x61,
            x139 * x176 * x32,
            x137 * x179 * x32,
            x135 * x184 * x32,
            x141 * x176 * x30,
            x139 * x179 * x30,
            x137 * x184 * x30,
            x135 * x187 * x30,
            x146 * x176 * x26,
            x141 * x179 * x26,
            x139 * x184 * x26,
            x137 * x187 * x26,
            x135 * x188 * x26,
            x194 * x63 * x90,
            x110 * x194 * x61,
            x195 * x61 * x90,
            x194 * x32 * x93,
            x110 * x195 * x32,
            x196 * x32 * x90,
            x112 * x194 * x30,
            x195 * x30 * x93,
            x110 * x196 * x30,
            x197 * x30 * x90,
            x126 * x194 * x26,
            x112 * x195 * x26,
            x196 * x26 * x93,
            x110 * x197 * x26,
            x198 * x26 * x90,
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
            x104 * x204 * x246,
            x104 * x213 * x248,
            x119 * x213 * x246,
            x104 * x223 * x250,
            x119 * x223 * x248,
            x108 * x223 * x246,
            x104 * x234 * x252,
            x119 * x234 * x250,
            x108 * x234 * x248,
            x122 * x234 * x246,
            x104 * x243 * x260,
            x119 * x243 * x252,
            x108 * x243 * x250,
            x122 * x243 * x248,
            x132 * x243 * x246,
            x149 * x201 * x204,
            x149 * x209 * x213,
            x151 * x201 * x213,
            x149 * x219 * x223,
            x151 * x209 * x223,
            x153 * x201 * x223,
            x149 * x232 * x234,
            x151 * x219 * x234,
            x153 * x209 * x234,
            x155 * x201 * x234,
            x149 * x242 * x243,
            x151 * x232 * x243,
            x153 * x219 * x243,
            x155 * x209 * x243,
            x160 * x201 * x243,
            x104 * x261 * x54,
            x104 * x266 * x80,
            x119 * x261 * x80,
            x104 * x271 * x97,
            x119 * x266 * x97,
            x108 * x261 * x97,
            x104 * x221 * x275,
            x119 * x221 * x271,
            x108 * x221 * x266,
            x122 * x221 * x261,
            x281 * x46,
            x275 * x282 * x89,
            x108 * x127 * x271,
            x122 * x127 * x266,
            x127 * x132 * x261,
            x149 * x246 * x54,
            x149 * x248 * x80,
            x151 * x246 * x80,
            x149 * x250 * x97,
            x151 * x248 * x97,
            x153 * x246 * x97,
            x149 * x221 * x252,
            x151 * x221 * x250,
            x153 * x221 * x248,
            x155 * x221 * x246,
            x147 * x260 * x282,
            x127 * x151 * x252,
            x127 * x153 * x250,
            x127 * x155 * x248,
            x127 * x160 * x246,
            x176 * x201 * x54,
            x176 * x209 * x80,
            x179 * x201 * x80,
            x176 * x219 * x97,
            x179 * x209 * x97,
            x184 * x201 * x97,
            x176 * x221 * x232,
            x179 * x219 * x221,
            x184 * x209 * x221,
            x187 * x201 * x221,
            x127 * x176 * x242,
            x127 * x179 * x232,
            x127 * x184 * x219,
            x127 * x187 * x209,
            x127 * x188 * x201,
            x104 * x284 * x52,
            x104 * x285 * x38,
            x119 * x284 * x38,
            x104 * x286 * x35,
            x119 * x285 * x35,
            x108 * x284 * x35,
            x287 * x3,
            x286 * x288 * x89,
            x108 * x12 * x285,
            x12 * x122 * x284,
            x280
            * (
                x0
                * (
                    x166 * (x277 + x278)
                    + 2.0 * x258
                    + 2.0 * x259
                    + 4.0 * x273
                    + 4.0 * x274
                    + x28 * (x146 + x254 + 5.0 * x256 + 4.0 * x272)
                )
                + x133 * x279
            ),
            x287 * x89,
            x108 * x286 * x7,
            x122 * x285 * x7,
            x132 * x284 * x7,
            x149 * x261 * x52,
            x149 * x266 * x38,
            x151 * x261 * x38,
            x149 * x271 * x35,
            x151 * x266 * x35,
            x153 * x261 * x35,
            x147 * x275 * x288,
            x12 * x151 * x271,
            x12 * x153 * x266,
            x12 * x155 * x261,
            x147 * x281,
            x151 * x275 * x7,
            x153 * x271 * x7,
            x155 * x266 * x7,
            x160 * x261 * x7,
            x176 * x246 * x52,
            x176 * x248 * x38,
            x179 * x246 * x38,
            x176 * x250 * x35,
            x179 * x248 * x35,
            x184 * x246 * x35,
            x12 * x176 * x252,
            x12 * x179 * x250,
            x12 * x184 * x248,
            x12 * x187 * x246,
            x176 * x260 * x7,
            x179 * x252 * x7,
            x184 * x250 * x7,
            x187 * x248 * x7,
            x188 * x246 * x7,
            x194 * x201 * x52,
            x194 * x209 * x38,
            x195 * x201 * x38,
            x194 * x219 * x35,
            x195 * x209 * x35,
            x196 * x201 * x35,
            x12 * x194 * x232,
            x12 * x195 * x219,
            x12 * x196 * x209,
            x12 * x197 * x201,
            x194 * x242 * x7,
            x195 * x232 * x7,
            x196 * x219 * x7,
            x197 * x209 * x7,
            x198 * x201 * x7,
            x205 * x291 * x90,
            x110 * x214 * x291,
            x214 * x295 * x90,
            x224 * x291 * x93,
            x110 * x224 * x295,
            x224 * x300 * x90,
            x112 * x235 * x291,
            x235 * x295 * x93,
            x110 * x235 * x300,
            x235 * x308 * x90,
            x126 * x244 * x291,
            x112 * x244 * x295,
            x244 * x300 * x93,
            x110 * x244 * x308,
            x244 * x315 * x90,
            x135 * x204 * x291,
            x137 * x213 * x291,
            x135 * x213 * x295,
            x139 * x223 * x291,
            x137 * x223 * x295,
            x135 * x223 * x300,
            x141 * x234 * x291,
            x139 * x234 * x295,
            x137 * x234 * x300,
            x135 * x234 * x308,
            x146 * x243 * x291,
            x141 * x243 * x295,
            x139 * x243 * x300,
            x137 * x243 * x308,
            x135 * x243 * x315,
            x204 * x317 * x90,
            x110 * x213 * x317,
            x213 * x319 * x90,
            x223 * x317 * x93,
            x110 * x223 * x319,
            x223 * x321 * x90,
            x112 * x234 * x317,
            x234 * x319 * x93,
            x110 * x234 * x321,
            x234 * x323 * x90,
            x126 * x243 * x317,
            x112 * x243 * x319,
            x243 * x321 * x93,
            x110 * x243 * x323,
            x243 * x331 * x90,
            x162 * x291 * x54,
            x165 * x291 * x80,
            x162 * x295 * x80,
            x170 * x291 * x97,
            x165 * x295 * x97,
            x162 * x300 * x97,
            x173 * x221 * x291,
            x170 * x221 * x295,
            x165 * x221 * x300,
            x162 * x221 * x308,
            x127 * x174 * x291,
            x127 * x173 * x295,
            x127 * x170 * x300,
            x127 * x165 * x308,
            x127 * x162 * x315,
            x135 * x317 * x54,
            x137 * x317 * x80,
            x135 * x319 * x80,
            x139 * x317 * x97,
            x137 * x319 * x97,
            x135 * x321 * x97,
            x141 * x221 * x317,
            x139 * x221 * x319,
            x137 * x221 * x321,
            x135 * x221 * x323,
            x127 * x146 * x317,
            x127 * x141 * x319,
            x127 * x139 * x321,
            x127 * x137 * x323,
            x133 * x331 * x333,
            x334 * x54 * x90,
            x110 * x334 * x80,
            x339 * x80 * x90,
            x334 * x93 * x97,
            x110 * x339 * x97,
            x344 * x90 * x97,
            x112 * x221 * x334,
            x221 * x339 * x93,
            x110 * x221 * x344,
            x221 * x348 * x90,
            x126 * x127 * x334,
            x112 * x127 * x339,
            x127 * x344 * x93,
            x333 * x348 * x77,
            x353 * x46,
            x189 * x291 * x52,
            x190 * x291 * x38,
            x189 * x295 * x38,
            x191 * x291 * x35,
            x190 * x295 * x35,
            x189 * x300 * x35,
            x12 * x192 * x291,
            x12 * x191 * x295,
            x12 * x190 * x300,
            x12 * x189 * x308,
            x193 * x291 * x7,
            x192 * x295 * x7,
            x191 * x300 * x7,
            x190 * x308 * x7,
            x189 * x315 * x7,
            x162 * x317 * x52,
            x165 * x317 * x38,
            x162 * x319 * x38,
            x170 * x317 * x35,
            x165 * x319 * x35,
            x162 * x321 * x35,
            x12 * x173 * x317,
            x12 * x170 * x319,
            x12 * x165 * x321,
            x12 * x162 * x323,
            x174 * x317 * x7,
            x173 * x319 * x7,
            x170 * x321 * x7,
            x165 * x323 * x7,
            x162 * x331 * x7,
            x135 * x334 * x52,
            x137 * x334 * x38,
            x135 * x339 * x38,
            x139 * x334 * x35,
            x137 * x339 * x35,
            x135 * x344 * x35,
            x12 * x141 * x334,
            x12 * x139 * x339,
            x12 * x137 * x344,
            x133 * x348 * x354,
            x146 * x334 * x7,
            x141 * x339 * x7,
            x139 * x344 * x7,
            x137 * x348 * x7,
            x133 * x353,
            x356 * x52 * x90,
            x110 * x356 * x38,
            x357 * x38 * x90,
            x35 * x356 * x93,
            x110 * x35 * x357,
            x35 * x358 * x90,
            x112 * x12 * x356,
            x12 * x357 * x93,
            x354 * x358 * x77,
            x3 * x359,
            x126 * x356 * x7,
            x112 * x357 * x7,
            x358 * x7 * x93,
            x359 * x77,
            x332
            * (
                x0
                * (
                    x180 * (x350 + x351)
                    + x28 * (x160 + x325 + 5.0 * x327 + 4.0 * x345)
                    + 2.0 * x329
                    + 2.0 * x330
                    + 4.0 * x346
                    + 4.0 * x347
                )
                + x147 * x352
            ),
        ]
    )


def diag_quadrupole3d_40(a, A, b, B, C):
    """Cartesian 3D (gs) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - A[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3**2 * x7
    x9 = x0 * x7
    x10 = 3.0 * x9
    x11 = 2.0 * x3
    x12 = -x2 - C[0]
    x13 = x12 * x7
    x14 = x10 + x11 * x13
    x15 = 2.0 * x0
    x16 = x3 * x7
    x17 = x0 * (x13 + x16)
    x18 = x3 * (x12 * x16 + x9)
    x19 = x12**2 * x7
    x20 = x0 * (x14 + x19)
    x21 = x19 + x9
    x22 = x21 * x3
    x23 = x13 * x15 + x22
    x24 = x23 * x3
    x25 = x20 + x24
    x26 = x0 * (4.0 * x0 * x13 + 2.0 * x17 + 2.0 * x18 + 2.0 * x22) + x25 * x3
    x27 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x28 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x29 = 3.14159265358979 * x1 * x28
    x30 = x27 * x29
    x31 = -x1 * (a * A[1] + b * B[1])
    x32 = -x31 - A[1]
    x33 = x26 * x30
    x34 = -x1 * (a * A[2] + b * B[2])
    x35 = -x34 - A[2]
    x36 = x27 * x6
    x37 = x32**2 * x36
    x38 = x0 * x36
    x39 = x37 + x38
    x40 = x28 * x6
    x41 = x35**2 * x40
    x42 = x0 * x40
    x43 = x41 + x42
    x44 = x32 * x36
    x45 = x15 * x44 + x32 * x39
    x46 = x35 * x40
    x47 = x15 * x46 + x35 * x43
    x48 = 3.0 * x38
    x49 = x0 * (3.0 * x37 + x48) + x32 * x45
    x50 = 3.0 * x42
    x51 = x0 * (3.0 * x41 + x50) + x35 * x47
    x52 = -x31 - C[1]
    x53 = x36 * x52**2
    x54 = x38 + x53
    x55 = x8 + x9
    x56 = x15 * x16 + x3 * x55
    x57 = x0 * (x10 + 3.0 * x8) + x3 * x56
    x58 = x32 * x54
    x59 = x36 * x52
    x60 = x15 * x59 + x58
    x61 = 2.0 * x32
    x62 = x48 + x59 * x61
    x63 = x0 * (x53 + x62)
    x64 = x32 * x60
    x65 = x63 + x64
    x66 = x0 * (x44 + x59)
    x67 = x32 * (x38 + x44 * x52)
    x68 = x0 * (4.0 * x38 * x52 + 2.0 * x58 + 2.0 * x66 + 2.0 * x67) + x32 * x65
    x69 = x29 * x5
    x70 = x68 * x69
    x71 = -x34 - C[2]
    x72 = x40 * x71**2
    x73 = x42 + x72
    x74 = x35 * x73
    x75 = x40 * x71
    x76 = x15 * x75 + x74
    x77 = 2.0 * x35
    x78 = x50 + x75 * x77
    x79 = x0 * (x72 + x78)
    x80 = x35 * x76
    x81 = x79 + x80
    x82 = 3.14159265358979 * x1 * x27 * x5
    x83 = x0 * (x46 + x75)
    x84 = x35 * (x42 + x46 * x71)
    x85 = x0 * (4.0 * x42 * x71 + 2.0 * x74 + 2.0 * x83 + 2.0 * x84) + x35 * x81
    x86 = x82 * x85

    # 45 item(s)
    return numpy.array(
        [
            x30
            * (
                x0 * (x11 * (x17 + x18) + x15 * (x14 + x8) + 3.0 * x20 + 3.0 * x24)
                + x26 * x3
            ),
            x32 * x33,
            x33 * x35,
            x25 * x39 * x40,
            x25 * x30 * x32 * x35,
            x25 * x36 * x43,
            x23 * x40 * x45,
            x23 * x39 * x46,
            x23 * x43 * x44,
            x23 * x36 * x47,
            x21 * x40 * x49,
            x21 * x45 * x46,
            x21 * x39 * x43,
            x21 * x44 * x47,
            x21 * x36 * x51,
            x40 * x54 * x57,
            x40 * x56 * x60,
            x46 * x54 * x56,
            x40 * x55 * x65,
            x46 * x55 * x60,
            x43 * x54 * x55,
            x3 * x70,
            x3 * x35 * x65 * x69,
            x16 * x43 * x60,
            x16 * x47 * x54,
            x69
            * (
                x0 * (x15 * (x37 + x62) + x61 * (x66 + x67) + 3.0 * x63 + 3.0 * x64)
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
            x16 * x45 * x73,
            x16 * x39 * x76,
            x3 * x32 * x81 * x82,
            x3 * x86,
            x49 * x7 * x73,
            x45 * x7 * x76,
            x39 * x7 * x81,
            x32 * x86,
            x82
            * (
                x0 * (x15 * (x41 + x78) + x77 * (x83 + x84) + 3.0 * x79 + 3.0 * x80)
                + x35 * x85
            ),
        ]
    )


def diag_quadrupole3d_41(a, A, b, B, C):
    """Cartesian 3D (gp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = 1.77245385090552 * numpy.sqrt(x1)
    x3 = -x1 * (a * A[0] + b * B[0])
    x4 = -x3 - A[0]
    x5 = a * b * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = x4 * x6
    x8 = x2 * x7
    x9 = -x3 - C[0]
    x10 = x2 * x6
    x11 = x10 * x9
    x12 = x0 * (x11 + x8)
    x13 = x0 * x10
    x14 = x8 * x9
    x15 = x4 * (x13 + x14)
    x16 = x12 + x15
    x17 = -x3 - B[0]
    x18 = x10 * x17
    x19 = x0 * (x18 + x8)
    x20 = x17 * x8
    x21 = x13 + x20
    x22 = x21 * x4
    x23 = x19 + x22
    x24 = x0 * (x11 + x18)
    x25 = x11 * x17
    x26 = x4 * (x13 + x25)
    x27 = 2.0 * x24 + 2.0 * x26
    x28 = 2.0 * x0
    x29 = 3.0 * x13
    x30 = x0 * (x14 + x20 + x25 + x29)
    x31 = x4 * (x24 + x26)
    x32 = 2.0 * x4
    x33 = x10 * x9**2
    x34 = x13 + x33
    x35 = x34 * x4
    x36 = x17 * x34
    x37 = 4.0 * x0 * x11
    x38 = x0 * (x27 + x35 + x36 + x37)
    x39 = x29 + x33
    x40 = x0 * (2.0 * x25 + x39)
    x41 = x11 * x28
    x42 = x36 + x41
    x43 = x4 * x42
    x44 = x40 + x43
    x45 = x4 * x44
    x46 = x11 * x32
    x47 = x0 * (x39 + x46)
    x48 = x35 + x41
    x49 = x4 * x48
    x50 = x47 + x49
    x51 = x0 * (2.0 * x12 + 2.0 * x15 + 2.0 * x35 + x37) + x4 * x50
    x52 = x38 + x45
    x53 = x0 * (2.0 * x30 + 2.0 * x31 + 2.0 * x40 + 2.0 * x43 + x50) + x4 * x52
    x54 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x55 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x56 = 3.14159265358979 * x1 * x55
    x57 = x54 * x56
    x58 = -x1 * (a * A[1] + b * B[1])
    x59 = -x58 - B[1]
    x60 = x10 * x4**2
    x61 = x29 + x60
    x62 = x57 * (x0 * (x16 * x32 + x28 * (x46 + x61) + 3.0 * x47 + 3.0 * x49) + x4 * x51)
    x63 = -x1 * (a * A[2] + b * B[2])
    x64 = -x63 - B[2]
    x65 = -x58 - A[1]
    x66 = x53 * x57
    x67 = x0 * x2
    x68 = x54 * x67
    x69 = x2 * x54
    x70 = x65 * x69
    x71 = x59 * x70
    x72 = x68 + x71
    x73 = x2 * x55
    x74 = x51 * x57
    x75 = -x63 - A[2]
    x76 = x55 * x67
    x77 = x73 * x75
    x78 = x64 * x77
    x79 = x76 + x78
    x80 = x65**2 * x69
    x81 = x68 + x80
    x82 = x59 * x69
    x83 = x0 * (x70 + x82)
    x84 = x65 * x72
    x85 = x83 + x84
    x86 = x64 * x73
    x87 = x73 * x75**2
    x88 = x76 + x87
    x89 = x0 * (x77 + x86)
    x90 = x75 * x79
    x91 = x89 + x90
    x92 = x28 * x70 + x65 * x81
    x93 = 2.0 * x65
    x94 = 3.0 * x68
    x95 = x80 + x94
    x96 = x0 * (x82 * x93 + x95) + x65 * x85
    x97 = x28 * x77 + x75 * x88
    x98 = 2.0 * x75
    x99 = 3.0 * x76
    x100 = x87 + x99
    x101 = x0 * (x100 + x86 * x98) + x75 * x91
    x102 = x0 * (3.0 * x80 + x94) + x65 * x92
    x103 = x0 * (3.0 * x83 + 3.0 * x84 + x92) + x65 * x96
    x104 = x0 * (3.0 * x87 + x99) + x75 * x97
    x105 = x0 * (3.0 * x89 + 3.0 * x90 + x97) + x101 * x75
    x106 = -x58 - C[1]
    x107 = x106**2 * x69
    x108 = x107 + x68
    x109 = x13 + x60
    x110 = x109 * x4 + x28 * x8
    x111 = x0 * (x18 * x32 + x61) + x23 * x4
    x112 = x0 * (x110 + 3.0 * x19 + 3.0 * x22) + x111 * x4
    x113 = x108 * x59
    x114 = x106 * x69
    x115 = x114 * x28
    x116 = x113 + x115
    x117 = x0 * (x29 + 3.0 * x60) + x110 * x4
    x118 = x108 * x65
    x119 = x115 + x118
    x120 = x106 * x82
    x121 = x107 + x94
    x122 = x0 * (2.0 * x120 + x121)
    x123 = x116 * x65
    x124 = x122 + x123
    x125 = x114 * x93
    x126 = x0 * (x121 + x125)
    x127 = x119 * x65
    x128 = x126 + x127
    x129 = 4.0 * x106 * x68
    x130 = x0 * (x114 + x82)
    x131 = x65 * (x120 + x68)
    x132 = 2.0 * x130 + 2.0 * x131
    x133 = x0 * (x113 + x118 + x129 + x132)
    x134 = x124 * x65
    x135 = x133 + x134
    x136 = x0 * (x114 + x70)
    x137 = x106 * x70
    x138 = x65 * (x137 + x68)
    x139 = x0 * (2.0 * x118 + x129 + 2.0 * x136 + 2.0 * x138) + x128 * x65
    x140 = x0 * (x120 + x137 + x71 + x94)
    x141 = x65 * (x130 + x131)
    x142 = x0 * (2.0 * x122 + 2.0 * x123 + x128 + 2.0 * x140 + 2.0 * x141) + x135 * x65
    x143 = x56 * x7
    x144 = x136 + x138
    x145 = x56 * x6
    x146 = x145 * (
        x0 * (3.0 * x126 + 3.0 * x127 + x144 * x93 + x28 * (x125 + x95)) + x139 * x65
    )
    x147 = x145 * x75
    x148 = -x63 - C[2]
    x149 = x148**2 * x73
    x150 = x149 + x76
    x151 = x150 * x64
    x152 = x148 * x73
    x153 = x152 * x28
    x154 = x151 + x153
    x155 = x150 * x75
    x156 = x153 + x155
    x157 = x148 * x86
    x158 = x149 + x99
    x159 = x0 * (2.0 * x157 + x158)
    x160 = x154 * x75
    x161 = x159 + x160
    x162 = x152 * x98
    x163 = x0 * (x158 + x162)
    x164 = x156 * x75
    x165 = x163 + x164
    x166 = 4.0 * x148 * x76
    x167 = x0 * (x152 + x86)
    x168 = x75 * (x157 + x76)
    x169 = 2.0 * x167 + 2.0 * x168
    x170 = x0 * (x151 + x155 + x166 + x169)
    x171 = x161 * x75
    x172 = x170 + x171
    x173 = 3.14159265358979 * x1 * x54
    x174 = x173 * x7
    x175 = x0 * (x152 + x77)
    x176 = x148 * x77
    x177 = x75 * (x176 + x76)
    x178 = x0 * (2.0 * x155 + x166 + 2.0 * x175 + 2.0 * x177) + x165 * x75
    x179 = x0 * (x157 + x176 + x78 + x99)
    x180 = x75 * (x167 + x168)
    x181 = x0 * (2.0 * x159 + 2.0 * x160 + x165 + 2.0 * x179 + 2.0 * x180) + x172 * x75
    x182 = x173 * x6
    x183 = x182 * x65
    x184 = x175 + x177
    x185 = x182 * (
        x0 * (3.0 * x163 + 3.0 * x164 + x184 * x98 + x28 * (x100 + x162)) + x178 * x75
    )

    # 135 item(s)
    return numpy.array(
        [
            x57
            * (
                x0
                * (
                    x28 * (x16 + x23 + x27)
                    + x32 * (x30 + x31)
                    + 3.0 * x38
                    + 3.0 * x45
                    + x51
                )
                + x4 * x53
            ),
            x59 * x62,
            x62 * x64,
            x65 * x66,
            x51 * x72 * x73,
            x64 * x65 * x74,
            x66 * x75,
            x59 * x74 * x75,
            x51 * x69 * x79,
            x52 * x73 * x81,
            x50 * x73 * x85,
            x50 * x81 * x86,
            x52 * x57 * x65 * x75,
            x50 * x72 * x77,
            x50 * x70 * x79,
            x52 * x69 * x88,
            x50 * x82 * x88,
            x50 * x69 * x91,
            x44 * x73 * x92,
            x48 * x73 * x96,
            x48 * x86 * x92,
            x44 * x77 * x81,
            x48 * x77 * x85,
            x48 * x79 * x81,
            x44 * x70 * x88,
            x48 * x72 * x88,
            x48 * x70 * x91,
            x44 * x69 * x97,
            x48 * x82 * x97,
            x101 * x48 * x69,
            x102 * x42 * x73,
            x103 * x34 * x73,
            x102 * x34 * x86,
            x42 * x77 * x92,
            x34 * x77 * x96,
            x34 * x79 * x92,
            x42 * x81 * x88,
            x34 * x85 * x88,
            x34 * x81 * x91,
            x42 * x70 * x97,
            x34 * x72 * x97,
            x101 * x34 * x70,
            x104 * x42 * x69,
            x104 * x34 * x82,
            x105 * x34 * x69,
            x108 * x112 * x73,
            x116 * x117 * x73,
            x108 * x117 * x86,
            x111 * x119 * x73,
            x110 * x124 * x73,
            x110 * x119 * x86,
            x108 * x111 * x77,
            x110 * x116 * x77,
            x108 * x110 * x79,
            x128 * x23 * x73,
            x109 * x135 * x73,
            x109 * x128 * x86,
            x119 * x23 * x77,
            x109 * x124 * x77,
            x109 * x119 * x79,
            x108 * x23 * x88,
            x109 * x116 * x88,
            x108 * x109 * x91,
            x139 * x21 * x73,
            x142 * x143,
            x139 * x143 * x64,
            x128 * x21 * x77,
            x135 * x143 * x75,
            x128 * x79 * x8,
            x119 * x21 * x88,
            x124 * x8 * x88,
            x119 * x8 * x91,
            x108 * x21 * x97,
            x116 * x8 * x97,
            x101 * x108 * x8,
            x146 * x17,
            x145
            * (
                x0
                * (
                    3.0 * x133
                    + 3.0 * x134
                    + x139
                    + x28 * (x132 + x144 + x85)
                    + x93 * (x140 + x141)
                )
                + x142 * x65
            ),
            x146 * x64,
            x139 * x147 * x17,
            x142 * x147,
            x10 * x139 * x79,
            x128 * x18 * x88,
            x10 * x135 * x88,
            x10 * x128 * x91,
            x119 * x18 * x97,
            x10 * x124 * x97,
            x10 * x101 * x119,
            x104 * x108 * x18,
            x10 * x104 * x116,
            x10 * x105 * x108,
            x112 * x150 * x69,
            x117 * x150 * x82,
            x117 * x154 * x69,
            x111 * x150 * x70,
            x110 * x150 * x72,
            x110 * x154 * x70,
            x111 * x156 * x69,
            x110 * x156 * x82,
            x110 * x161 * x69,
            x150 * x23 * x81,
            x109 * x150 * x85,
            x109 * x154 * x81,
            x156 * x23 * x70,
            x109 * x156 * x72,
            x109 * x161 * x70,
            x165 * x23 * x69,
            x109 * x165 * x82,
            x109 * x172 * x69,
            x150 * x21 * x92,
            x150 * x8 * x96,
            x154 * x8 * x92,
            x156 * x21 * x81,
            x156 * x8 * x85,
            x161 * x8 * x81,
            x165 * x21 * x70,
            x165 * x72 * x8,
            x172 * x174 * x65,
            x178 * x21 * x69,
            x174 * x178 * x59,
            x174 * x181,
            x102 * x150 * x18,
            x10 * x103 * x150,
            x10 * x102 * x154,
            x156 * x18 * x92,
            x10 * x156 * x96,
            x10 * x161 * x92,
            x165 * x18 * x81,
            x10 * x165 * x85,
            x10 * x172 * x81,
            x17 * x178 * x183,
            x10 * x178 * x72,
            x181 * x183,
            x17 * x185,
            x185 * x59,
            x182
            * (
                x0
                * (
                    3.0 * x170
                    + 3.0 * x171
                    + x178
                    + x28 * (x169 + x184 + x91)
                    + x98 * (x179 + x180)
                )
                + x181 * x75
            ),
        ]
    )


def diag_quadrupole3d_42(a, A, b, B, C):
    """Cartesian 3D (gd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - A[0]
    x4 = 2.0 * x3
    x5 = -x2 - B[0]
    x6 = a * b * x1
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = 1.77245385090552 * numpy.sqrt(x1)
    x9 = x7 * x8
    x10 = x5 * x9
    x11 = x10 * x4
    x12 = x5**2 * x9
    x13 = x0 * x9
    x14 = 3.0 * x13
    x15 = x12 + x14
    x16 = x0 * (x11 + x15)
    x17 = x12 + x13
    x18 = x17 * x3
    x19 = 2.0 * x0
    x20 = x10 * x19 + x18
    x21 = x20 * x3
    x22 = x16 + x21
    x23 = -x2 - C[0]
    x24 = x10 * x23
    x25 = 2.0 * x24
    x26 = x0 * (x15 + x25)
    x27 = x23 * x9
    x28 = x0 * (x10 + x27)
    x29 = x13 + x24
    x30 = x29 * x5
    x31 = x3 * (x28 + x30)
    x32 = 2.0 * x26 + 2.0 * x31
    x33 = x10 * x3
    x34 = x27 * x3
    x35 = x0 * (x14 + x24 + x33 + x34)
    x36 = x29 * x3
    x37 = x3 * (x28 + x36)
    x38 = 2.0 * x35 + 2.0 * x37
    x39 = x23**2 * x9
    x40 = x13 + x39
    x41 = x40 * x5
    x42 = x19 * x27
    x43 = x41 + x42
    x44 = x3 * x43
    x45 = 2.0 * x44
    x46 = x14 + x39
    x47 = x0 * (x25 + x46)
    x48 = x27 * x4
    x49 = x0 * (x46 + x48)
    x50 = x3 * x40
    x51 = x42 + x50
    x52 = x3 * x51
    x53 = x49 + x52
    x54 = x0 * (x38 + x45 + 2.0 * x47 + x53)
    x55 = 2.0 * x36
    x56 = x0 * (x20 + 3.0 * x28 + x30 + x55)
    x57 = x3 * (x26 + x31)
    x58 = 2.0 * x28
    x59 = 4.0 * x13
    x60 = x23 * x59
    x61 = x58 + x60
    x62 = x0 * (x41 + x50 + x55 + x61)
    x63 = x44 + x47
    x64 = x3 * x63
    x65 = x62 + x64
    x66 = x3 * x65
    x67 = x43 * x5
    x68 = x0 * (x32 + x45 + 3.0 * x47 + x67)
    x69 = x0 * (2.0 * x30 + 2.0 * x41 + x61)
    x70 = x47 + x67
    x71 = x3 * x70
    x72 = x69 + x71
    x73 = x3 * x72
    x74 = x68 + x73
    x75 = (
        x0 * (2.0 * x56 + 2.0 * x57 + 2.0 * x62 + 2.0 * x64 + 2.0 * x69 + 2.0 * x71)
        + x3 * x74
    )
    x76 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x77 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x78 = 3.14159265358979 * x1 * x77
    x79 = x76 * x78
    x80 = -x1 * (a * A[1] + b * B[1])
    x81 = -x80 - B[1]
    x82 = x3 * x9
    x83 = x0 * (x27 + x82)
    x84 = x3 * (x13 + x34)
    x85 = x83 + x84
    x86 = x0 * (x10 + x82)
    x87 = x13 + x33
    x88 = x3 * x87
    x89 = x86 + x88
    x90 = x0 * (2.0 * x50 + x60 + 2.0 * x83 + 2.0 * x84) + x3 * x53
    x91 = x54 + x66
    x92 = x79 * (
        x0
        * (x19 * (x55 + x58 + x85 + x89) + x4 * (x35 + x37) + 3.0 * x62 + 3.0 * x64 + x90)
        + x3 * x91
    )
    x93 = -x1 * (a * A[2] + b * B[2])
    x94 = -x93 - B[2]
    x95 = x76 * x8
    x96 = x81**2 * x95
    x97 = x0 * x95
    x98 = x96 + x97
    x99 = x3**2 * x9
    x100 = x14 + x99
    x101 = x0 * (x19 * (x100 + x48) + x4 * x85 + 3.0 * x49 + 3.0 * x52) + x3 * x90
    x102 = x77 * x8
    x103 = x79 * x94
    x104 = x102 * x94**2
    x105 = x0 * x102
    x106 = x104 + x105
    x107 = -x80 - A[1]
    x108 = x75 * x79
    x109 = x107 * x95
    x110 = x109 * x81
    x111 = x110 + x97
    x112 = x107 * x98
    x113 = x81 * x95
    x114 = x112 + x113 * x19
    x115 = x102 * x94
    x116 = -x93 - A[2]
    x117 = x116 * x79
    x118 = x102 * x116
    x119 = x118 * x94
    x120 = x105 + x119
    x121 = x106 * x116
    x122 = x115 * x19 + x121
    x123 = x107**2 * x95
    x124 = x123 + x97
    x125 = x0 * (x109 + x113)
    x126 = x107 * x111
    x127 = x125 + x126
    x128 = 3.0 * x97
    x129 = 2.0 * x107
    x130 = x113 * x129 + x128
    x131 = x0 * (x130 + x96)
    x132 = x107 * x114
    x133 = x131 + x132
    x134 = x102 * x116**2
    x135 = x105 + x134
    x136 = x0 * (x115 + x118)
    x137 = x116 * x120
    x138 = x136 + x137
    x139 = 3.0 * x105
    x140 = 2.0 * x116
    x141 = x115 * x140 + x139
    x142 = x0 * (x104 + x141)
    x143 = x116 * x122
    x144 = x142 + x143
    x145 = x107 * x124 + x109 * x19
    x146 = x0 * (x123 + x130)
    x147 = x107 * x127
    x148 = x146 + x147
    x149 = 4.0 * x97
    x150 = x0 * (2.0 * x112 + 2.0 * x125 + 2.0 * x126 + x149 * x81) + x107 * x133
    x151 = x116 * x135 + x118 * x19
    x152 = x0 * (x134 + x141)
    x153 = x116 * x138
    x154 = x152 + x153
    x155 = 4.0 * x105
    x156 = x0 * (2.0 * x121 + 2.0 * x136 + 2.0 * x137 + x155 * x94) + x116 * x144
    x157 = x0 * (3.0 * x123 + x128) + x107 * x145
    x158 = x0 * (3.0 * x125 + 3.0 * x126 + x145) + x107 * x148
    x159 = x0 * (3.0 * x131 + 3.0 * x132 + 2.0 * x146 + 2.0 * x147) + x107 * x150
    x160 = x0 * (3.0 * x134 + x139) + x116 * x151
    x161 = x0 * (3.0 * x136 + 3.0 * x137 + x151) + x116 * x154
    x162 = x0 * (3.0 * x142 + 3.0 * x143 + 2.0 * x152 + 2.0 * x153) + x116 * x156
    x163 = -x80 - C[1]
    x164 = x163**2 * x95
    x165 = x164 + x97
    x166 = x0 * (x100 + x11)
    x167 = x3 * x89
    x168 = x0 * (2.0 * x18 + x5 * x59 + 2.0 * x86 + 2.0 * x88) + x22 * x3
    x169 = x0 * (3.0 * x16 + 2.0 * x166 + 2.0 * x167 + 3.0 * x21) + x168 * x3
    x170 = x165 * x81
    x171 = x163 * x95
    x172 = x171 * x19
    x173 = x170 + x172
    x174 = x13 + x99
    x175 = x174 * x3 + x19 * x82
    x176 = x166 + x167
    x177 = x0 * (x175 + 3.0 * x86 + 3.0 * x88) + x176 * x3
    x178 = x113 * x163
    x179 = 2.0 * x178
    x180 = x128 + x164
    x181 = x0 * (x179 + x180)
    x182 = x173 * x81
    x183 = x181 + x182
    x184 = x0 * (x14 + 3.0 * x99) + x175 * x3
    x185 = x107 * x165
    x186 = x172 + x185
    x187 = x107 * x173
    x188 = x181 + x187
    x189 = x178 + x97
    x190 = x189 * x81
    x191 = x0 * (x113 + x171)
    x192 = 2.0 * x191
    x193 = x149 * x163
    x194 = x192 + x193
    x195 = x0 * (2.0 * x170 + 2.0 * x190 + x194)
    x196 = x107 * x183
    x197 = x195 + x196
    x198 = x129 * x171
    x199 = x0 * (x180 + x198)
    x200 = x107 * x186
    x201 = x199 + x200
    x202 = x107 * x189
    x203 = 2.0 * x202
    x204 = x0 * (x170 + x185 + x194 + x203)
    x205 = x107 * x188
    x206 = x204 + x205
    x207 = 2.0 * x187
    x208 = x0 * (x128 + x179 + x96)
    x209 = x107 * (x190 + x191)
    x210 = 2.0 * x208 + 2.0 * x209
    x211 = x0 * (3.0 * x181 + x182 + x207 + x210)
    x212 = x107 * x197
    x213 = x211 + x212
    x214 = x0 * (x109 + x171)
    x215 = x109 * x163
    x216 = x107 * (x215 + x97)
    x217 = x0 * (2.0 * x185 + x193 + 2.0 * x214 + 2.0 * x216) + x107 * x201
    x218 = x0 * (x110 + x128 + x178 + x215)
    x219 = x107 * (x191 + x202)
    x220 = 2.0 * x218 + 2.0 * x219
    x221 = x0 * (2.0 * x181 + x201 + x207 + x220)
    x222 = x107 * x206
    x223 = x221 + x222
    x224 = x0 * (x114 + x190 + 3.0 * x191 + x203)
    x225 = x107 * (x208 + x209)
    x226 = (
        x0 * (2.0 * x195 + 2.0 * x196 + 2.0 * x204 + 2.0 * x205 + 2.0 * x224 + 2.0 * x225)
        + x107 * x213
    )
    x227 = x7 * x78
    x228 = x226 * x227
    x229 = x227 * x3
    x230 = x214 + x216
    x231 = (
        x0 * (x129 * x230 + x19 * (x123 + x128 + x198) + 3.0 * x199 + 3.0 * x200)
        + x107 * x217
    )
    x232 = x227 * (
        x0
        * (
            x129 * (x218 + x219)
            + x19 * (x127 + x192 + x203 + x230)
            + 3.0 * x204
            + 3.0 * x205
            + x217
        )
        + x107 * x223
    )
    x233 = x227 * x5
    x234 = -x93 - C[2]
    x235 = x102 * x234**2
    x236 = x105 + x235
    x237 = x236 * x94
    x238 = x102 * x234
    x239 = x19 * x238
    x240 = x237 + x239
    x241 = x115 * x234
    x242 = 2.0 * x241
    x243 = x139 + x235
    x244 = x0 * (x242 + x243)
    x245 = x240 * x94
    x246 = x244 + x245
    x247 = x116 * x236
    x248 = x239 + x247
    x249 = x116 * x240
    x250 = x244 + x249
    x251 = x105 + x241
    x252 = x251 * x94
    x253 = x0 * (x115 + x238)
    x254 = 2.0 * x253
    x255 = x155 * x234
    x256 = x254 + x255
    x257 = x0 * (2.0 * x237 + 2.0 * x252 + x256)
    x258 = x116 * x246
    x259 = x257 + x258
    x260 = x140 * x238
    x261 = x0 * (x243 + x260)
    x262 = x116 * x248
    x263 = x261 + x262
    x264 = x116 * x251
    x265 = 2.0 * x264
    x266 = x0 * (x237 + x247 + x256 + x265)
    x267 = x116 * x250
    x268 = x266 + x267
    x269 = 2.0 * x249
    x270 = x0 * (x104 + x139 + x242)
    x271 = x116 * (x252 + x253)
    x272 = 2.0 * x270 + 2.0 * x271
    x273 = x0 * (3.0 * x244 + x245 + x269 + x272)
    x274 = x116 * x259
    x275 = x273 + x274
    x276 = 3.14159265358979 * x1 * x7 * x76
    x277 = x276 * x3
    x278 = x0 * (x118 + x238)
    x279 = x118 * x234
    x280 = x116 * (x105 + x279)
    x281 = x0 * (2.0 * x247 + x255 + 2.0 * x278 + 2.0 * x280) + x116 * x263
    x282 = x0 * (x119 + x139 + x241 + x279)
    x283 = x116 * (x253 + x264)
    x284 = 2.0 * x282 + 2.0 * x283
    x285 = x0 * (2.0 * x244 + x263 + x269 + x284)
    x286 = x116 * x268
    x287 = x285 + x286
    x288 = x0 * (x122 + x252 + 3.0 * x253 + x265)
    x289 = x116 * (x270 + x271)
    x290 = (
        x0 * (2.0 * x257 + 2.0 * x258 + 2.0 * x266 + 2.0 * x267 + 2.0 * x288 + 2.0 * x289)
        + x116 * x275
    )
    x291 = x276 * x290
    x292 = x276 * x5
    x293 = x278 + x280
    x294 = (
        x0 * (x140 * x293 + x19 * (x134 + x139 + x260) + 3.0 * x261 + 3.0 * x262)
        + x116 * x281
    )
    x295 = x276 * (
        x0
        * (
            x140 * (x282 + x283)
            + x19 * (x138 + x254 + x265 + x293)
            + 3.0 * x266
            + 3.0 * x267
            + x281
        )
        + x116 * x287
    )

    # 270 item(s)
    return numpy.array(
        [
            x79
            * (
                x0
                * (
                    x19 * (x22 + x32 + x38)
                    + x4 * (x56 + x57)
                    + 2.0 * x54
                    + 2.0 * x66
                    + 3.0 * x68
                    + 3.0 * x73
                )
                + x3 * x75
            ),
            x81 * x92,
            x92 * x94,
            x101 * x102 * x98,
            x101 * x103 * x81,
            x101 * x106 * x95,
            x107 * x108,
            x102 * x111 * x91,
            x103 * x107 * x91,
            x102 * x114 * x90,
            x111 * x115 * x90,
            x106 * x109 * x90,
            x108 * x116,
            x117 * x81 * x91,
            x120 * x91 * x95,
            x118 * x90 * x98,
            x113 * x120 * x90,
            x122 * x90 * x95,
            x102 * x124 * x74,
            x102 * x127 * x65,
            x115 * x124 * x65,
            x102 * x133 * x53,
            x115 * x127 * x53,
            x106 * x124 * x53,
            x107 * x117 * x74,
            x111 * x118 * x65,
            x109 * x120 * x65,
            x114 * x118 * x53,
            x111 * x120 * x53,
            x109 * x122 * x53,
            x135 * x74 * x95,
            x113 * x135 * x65,
            x138 * x65 * x95,
            x135 * x53 * x98,
            x113 * x138 * x53,
            x144 * x53 * x95,
            x102 * x145 * x72,
            x102 * x148 * x63,
            x115 * x145 * x63,
            x102 * x150 * x51,
            x115 * x148 * x51,
            x106 * x145 * x51,
            x118 * x124 * x72,
            x118 * x127 * x63,
            x120 * x124 * x63,
            x118 * x133 * x51,
            x120 * x127 * x51,
            x122 * x124 * x51,
            x109 * x135 * x72,
            x111 * x135 * x63,
            x109 * x138 * x63,
            x114 * x135 * x51,
            x111 * x138 * x51,
            x109 * x144 * x51,
            x151 * x72 * x95,
            x113 * x151 * x63,
            x154 * x63 * x95,
            x151 * x51 * x98,
            x113 * x154 * x51,
            x156 * x51 * x95,
            x102 * x157 * x70,
            x102 * x158 * x43,
            x115 * x157 * x43,
            x102 * x159 * x40,
            x115 * x158 * x40,
            x106 * x157 * x40,
            x118 * x145 * x70,
            x118 * x148 * x43,
            x120 * x145 * x43,
            x118 * x150 * x40,
            x120 * x148 * x40,
            x122 * x145 * x40,
            x124 * x135 * x70,
            x127 * x135 * x43,
            x124 * x138 * x43,
            x133 * x135 * x40,
            x127 * x138 * x40,
            x124 * x144 * x40,
            x109 * x151 * x70,
            x111 * x151 * x43,
            x109 * x154 * x43,
            x114 * x151 * x40,
            x111 * x154 * x40,
            x109 * x156 * x40,
            x160 * x70 * x95,
            x113 * x160 * x43,
            x161 * x43 * x95,
            x160 * x40 * x98,
            x113 * x161 * x40,
            x162 * x40 * x95,
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
            x102 * x201 * x22,
            x102 * x206 * x89,
            x115 * x201 * x89,
            x102 * x174 * x213,
            x115 * x174 * x206,
            x106 * x174 * x201,
            x118 * x186 * x22,
            x118 * x188 * x89,
            x120 * x186 * x89,
            x118 * x174 * x197,
            x120 * x174 * x188,
            x122 * x174 * x186,
            x135 * x165 * x22,
            x135 * x173 * x89,
            x138 * x165 * x89,
            x135 * x174 * x183,
            x138 * x173 * x174,
            x144 * x165 * x174,
            x102 * x20 * x217,
            x102 * x223 * x87,
            x115 * x217 * x87,
            x228 * x3,
            x223 * x229 * x94,
            x106 * x217 * x82,
            x118 * x20 * x201,
            x118 * x206 * x87,
            x120 * x201 * x87,
            x116 * x213 * x229,
            x120 * x206 * x82,
            x122 * x201 * x82,
            x135 * x186 * x20,
            x135 * x188 * x87,
            x138 * x186 * x87,
            x135 * x197 * x82,
            x138 * x188 * x82,
            x144 * x186 * x82,
            x151 * x165 * x20,
            x151 * x173 * x87,
            x154 * x165 * x87,
            x151 * x183 * x82,
            x154 * x173 * x82,
            x156 * x165 * x82,
            x102 * x17 * x231,
            x232 * x5,
            x231 * x233 * x94,
            x227
            * (
                x0
                * (
                    x129 * (x224 + x225)
                    + x19 * (x133 + x210 + x220)
                    + 3.0 * x211
                    + 3.0 * x212
                    + 2.0 * x221
                    + 2.0 * x222
                )
                + x107 * x226
            ),
            x232 * x94,
            x106 * x231 * x9,
            x118 * x17 * x217,
            x116 * x223 * x233,
            x10 * x120 * x217,
            x116 * x228,
            x120 * x223 * x9,
            x122 * x217 * x9,
            x135 * x17 * x201,
            x10 * x135 * x206,
            x10 * x138 * x201,
            x135 * x213 * x9,
            x138 * x206 * x9,
            x144 * x201 * x9,
            x151 * x17 * x186,
            x10 * x151 * x188,
            x10 * x154 * x186,
            x151 * x197 * x9,
            x154 * x188 * x9,
            x156 * x186 * x9,
            x160 * x165 * x17,
            x10 * x160 * x173,
            x10 * x161 * x165,
            x160 * x183 * x9,
            x161 * x173 * x9,
            x162 * x165 * x9,
            x169 * x236 * x95,
            x113 * x177 * x236,
            x177 * x240 * x95,
            x184 * x236 * x98,
            x113 * x184 * x240,
            x184 * x246 * x95,
            x109 * x168 * x236,
            x111 * x176 * x236,
            x109 * x176 * x240,
            x114 * x175 * x236,
            x111 * x175 * x240,
            x109 * x175 * x246,
            x168 * x248 * x95,
            x113 * x176 * x248,
            x176 * x250 * x95,
            x175 * x248 * x98,
            x113 * x175 * x250,
            x175 * x259 * x95,
            x124 * x22 * x236,
            x127 * x236 * x89,
            x124 * x240 * x89,
            x133 * x174 * x236,
            x127 * x174 * x240,
            x124 * x174 * x246,
            x109 * x22 * x248,
            x111 * x248 * x89,
            x109 * x250 * x89,
            x114 * x174 * x248,
            x111 * x174 * x250,
            x109 * x174 * x259,
            x22 * x263 * x95,
            x113 * x263 * x89,
            x268 * x89 * x95,
            x174 * x263 * x98,
            x113 * x174 * x268,
            x174 * x275 * x95,
            x145 * x20 * x236,
            x148 * x236 * x87,
            x145 * x240 * x87,
            x150 * x236 * x82,
            x148 * x240 * x82,
            x145 * x246 * x82,
            x124 * x20 * x248,
            x127 * x248 * x87,
            x124 * x250 * x87,
            x133 * x248 * x82,
            x127 * x250 * x82,
            x124 * x259 * x82,
            x109 * x20 * x263,
            x111 * x263 * x87,
            x109 * x268 * x87,
            x114 * x263 * x82,
            x111 * x268 * x82,
            x107 * x275 * x277,
            x20 * x281 * x95,
            x113 * x281 * x87,
            x287 * x87 * x95,
            x281 * x82 * x98,
            x277 * x287 * x81,
            x291 * x3,
            x157 * x17 * x236,
            x10 * x158 * x236,
            x10 * x157 * x240,
            x159 * x236 * x9,
            x158 * x240 * x9,
            x157 * x246 * x9,
            x145 * x17 * x248,
            x10 * x148 * x248,
            x10 * x145 * x250,
            x150 * x248 * x9,
            x148 * x250 * x9,
            x145 * x259 * x9,
            x124 * x17 * x263,
            x10 * x127 * x263,
            x10 * x124 * x268,
            x133 * x263 * x9,
            x127 * x268 * x9,
            x124 * x275 * x9,
            x109 * x17 * x281,
            x10 * x111 * x281,
            x107 * x287 * x292,
            x114 * x281 * x9,
            x111 * x287 * x9,
            x107 * x291,
            x17 * x294 * x95,
            x292 * x294 * x81,
            x295 * x5,
            x294 * x9 * x98,
            x295 * x81,
            x276
            * (
                x0
                * (
                    x140 * (x288 + x289)
                    + x19 * (x144 + x272 + x284)
                    + 3.0 * x273
                    + 3.0 * x274
                    + 2.0 * x285
                    + 2.0 * x286
                )
                + x116 * x290
            ),
        ]
    )


def diag_quadrupole3d_43(a, A, b, B, C):
    """Cartesian 3D (gf) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x0 * x7
    x9 = -x2 - C[0]
    x10 = x3 * x7
    x11 = x10 * x9
    x12 = x11 + x8
    x13 = x12 * x3
    x14 = x7 * x9
    x15 = x0 * (x10 + x14)
    x16 = 3.0 * x15
    x17 = -x2 - A[0]
    x18 = x12 * x17
    x19 = 2.0 * x18
    x20 = x3**2 * x7
    x21 = x20 + x8
    x22 = x17 * x21
    x23 = 2.0 * x0
    x24 = x10 * x23
    x25 = x22 + x24
    x26 = x0 * (x13 + x16 + x19 + x25)
    x27 = 3.0 * x8
    x28 = 2.0 * x11 + x27
    x29 = x0 * (x20 + x28)
    x30 = x13 + x15
    x31 = x17 * x30
    x32 = x17 * (x29 + x31)
    x33 = x21 * x3
    x34 = x3 * x8
    x35 = x0 * (3.0 * x22 + x33 + 8.0 * x34)
    x36 = x0 * (3.0 * x20 + x27)
    x37 = x24 + x33
    x38 = x17 * x37
    x39 = x36 + x38
    x40 = x17 * x39
    x41 = x35 + x40
    x42 = x0 * (3.0 * x13 + x16 + x37)
    x43 = x3 * x30
    x44 = x17 * (x29 + x43)
    x45 = 2.0 * x42 + 2.0 * x44
    x46 = x0 * (4.0 * x29 + 3.0 * x31 + x39 + x43)
    x47 = x17 * (x42 + x44)
    x48 = 2.0 * x17
    x49 = x7 * x9**2
    x50 = x0 * (x28 + x49)
    x51 = x49 + x8
    x52 = x3 * x51
    x53 = x14 * x23
    x54 = x52 + x53
    x55 = x3 * x54
    x56 = x50 + x55
    x57 = x3 * x56
    x58 = x17 * x56
    x59 = 2.0 * x15
    x60 = 4.0 * x8 * x9
    x61 = x59 + x60
    x62 = x0 * (2.0 * x13 + 2.0 * x52 + x61)
    x63 = x0 * (x45 + x57 + 3.0 * x58 + 4.0 * x62)
    x64 = x17 * x51
    x65 = x0 * (x19 + x52 + x61 + x64)
    x66 = x17 * x54
    x67 = x50 + x66
    x68 = x17 * x67
    x69 = x0 * (2.0 * x26 + 2.0 * x32 + 2.0 * x58 + 2.0 * x62 + 2.0 * x65 + 2.0 * x68)
    x70 = 2.0 * x29
    x71 = 3.0 * x50 + x70
    x72 = x0 * (2.0 * x43 + 3.0 * x55 + x71)
    x73 = x57 + x62
    x74 = x17 * x73
    x75 = x72 + x74
    x76 = x17 * x75
    x77 = 2.0 * x31
    x78 = 2.0 * x66
    x79 = x0 * (x55 + x71 + x77 + x78)
    x80 = x58 + x62
    x81 = x17 * x80
    x82 = x79 + x81
    x83 = x17 * x82
    x84 = 3.0 * x79 + 3.0 * x81
    x85 = x63 + x76
    x86 = x0 * (2.0 * x46 + 2.0 * x47 + 2.0 * x72 + 2.0 * x74 + x84) + x17 * x85
    x87 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x88 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x89 = 3.14159265358979 * x1 * x88
    x90 = x87 * x89
    x91 = -x1 * (a * A[1] + b * B[1])
    x92 = -x91 - B[1]
    x93 = x10 * x48 + x27
    x94 = x0 * (x20 + x93)
    x95 = x17 * x25
    x96 = x94 + x95
    x97 = x10 * x17
    x98 = x14 * x17
    x99 = x0 * (x11 + x27 + x97 + x98)
    x100 = x17 * (x15 + x18)
    x101 = 2.0 * x100 + 2.0 * x99
    x102 = x14 * x48 + x27
    x103 = x0 * (x102 + x49)
    x104 = x53 + x64
    x105 = x104 * x17
    x106 = x103 + x105
    x107 = x0 * (x101 + x106 + 2.0 * x50 + x78)
    x108 = x65 + x68
    x109 = x108 * x17
    x110 = x69 + x83
    x111 = x90 * (
        x0
        * (
            2.0 * x107
            + 2.0 * x109
            + x23 * (x101 + x70 + x77 + x96)
            + x48 * (x26 + x32)
            + x84
        )
        + x110 * x17
    )
    x112 = -x1 * (a * A[2] + b * B[2])
    x113 = -x112 - B[2]
    x114 = x6 * x87
    x115 = x114 * x92**2
    x116 = x0 * x114
    x117 = x115 + x116
    x118 = x17 * x7
    x119 = x0 * (x118 + x14)
    x120 = x17 * (x8 + x98)
    x121 = x119 + x120
    x122 = x0 * (x10 + x118)
    x123 = x8 + x97
    x124 = x123 * x17
    x125 = x122 + x124
    x126 = x0 * (2.0 * x119 + 2.0 * x120 + x60 + 2.0 * x64) + x106 * x17
    x127 = x107 + x109
    x128 = (
        x0
        * (
            x126
            + x23 * (x121 + x125 + x19 + x59)
            + x48 * (x100 + x99)
            + 3.0 * x65
            + 3.0 * x68
        )
        + x127 * x17
    )
    x129 = x6 * x88
    x130 = x113 * x90
    x131 = x113**2 * x129
    x132 = x0 * x129
    x133 = x131 + x132
    x134 = x117 * x92
    x135 = x114 * x92
    x136 = x135 * x23
    x137 = x134 + x136
    x138 = x17**2 * x7
    x139 = x0 * (3.0 * x103 + 3.0 * x105 + x121 * x48 + x23 * (x102 + x138)) + x126 * x17
    x140 = x113 * x129
    x141 = x113 * x133
    x142 = x140 * x23
    x143 = x141 + x142
    x144 = -x91 - A[1]
    x145 = x86 * x90
    x146 = x114 * x144
    x147 = x146 * x92
    x148 = x116 + x147
    x149 = x117 * x144
    x150 = x136 + x149
    x151 = 3.0 * x116
    x152 = x0 * (3.0 * x115 + x151)
    x153 = x137 * x144
    x154 = x152 + x153
    x155 = -x112 - A[2]
    x156 = x155 * x90
    x157 = x129 * x155
    x158 = x113 * x157
    x159 = x132 + x158
    x160 = x133 * x155
    x161 = x142 + x160
    x162 = 3.0 * x132
    x163 = x0 * (3.0 * x131 + x162)
    x164 = x143 * x155
    x165 = x163 + x164
    x166 = x114 * x144**2
    x167 = x116 + x166
    x168 = x0 * (x135 + x146)
    x169 = x144 * x148
    x170 = x168 + x169
    x171 = 2.0 * x144
    x172 = x135 * x171 + x151
    x173 = x0 * (x115 + x172)
    x174 = x144 * x150
    x175 = x173 + x174
    x176 = x116 * x92
    x177 = x0 * (x134 + 3.0 * x149 + 8.0 * x176)
    x178 = x144 * x154
    x179 = x177 + x178
    x180 = x129 * x155**2
    x181 = x132 + x180
    x182 = x0 * (x140 + x157)
    x183 = x155 * x159
    x184 = x182 + x183
    x185 = 2.0 * x155
    x186 = x140 * x185 + x162
    x187 = x0 * (x131 + x186)
    x188 = x155 * x161
    x189 = x187 + x188
    x190 = x113 * x132
    x191 = x0 * (x141 + 3.0 * x160 + 8.0 * x190)
    x192 = x155 * x165
    x193 = x191 + x192
    x194 = x144 * x167 + x146 * x23
    x195 = x0 * (x166 + x172)
    x196 = x144 * x170
    x197 = x195 + x196
    x198 = x0 * (2.0 * x149 + 2.0 * x168 + 2.0 * x169 + 4.0 * x176)
    x199 = x144 * x175
    x200 = x198 + x199
    x201 = 3.0 * x173 + 3.0 * x174
    x202 = x0 * (2.0 * x152 + 2.0 * x153 + x201) + x144 * x179
    x203 = x155 * x181 + x157 * x23
    x204 = x0 * (x180 + x186)
    x205 = x155 * x184
    x206 = x204 + x205
    x207 = x0 * (2.0 * x160 + 2.0 * x182 + 2.0 * x183 + 4.0 * x190)
    x208 = x155 * x189
    x209 = x207 + x208
    x210 = 3.0 * x187 + 3.0 * x188
    x211 = x0 * (2.0 * x163 + 2.0 * x164 + x210) + x155 * x193
    x212 = x0 * (x151 + 3.0 * x166) + x144 * x194
    x213 = x0 * (3.0 * x168 + 3.0 * x169 + x194) + x144 * x197
    x214 = x0 * (2.0 * x195 + 2.0 * x196 + x201) + x144 * x200
    x215 = x0 * (3.0 * x177 + 3.0 * x178 + 3.0 * x198 + 3.0 * x199) + x144 * x202
    x216 = x0 * (x162 + 3.0 * x180) + x155 * x203
    x217 = x0 * (3.0 * x182 + 3.0 * x183 + x203) + x155 * x206
    x218 = x0 * (2.0 * x204 + 2.0 * x205 + x210) + x155 * x209
    x219 = x0 * (3.0 * x191 + 3.0 * x192 + 3.0 * x207 + 3.0 * x208) + x155 * x211
    x220 = -x91 - C[1]
    x221 = x114 * x220**2
    x222 = x116 + x221
    x223 = x0 * (2.0 * x122 + 2.0 * x124 + 2.0 * x22 + 4.0 * x34)
    x224 = x17 * x96
    x225 = 3.0 * x94 + 3.0 * x95
    x226 = x0 * (x225 + 2.0 * x36 + 2.0 * x38) + x17 * x41
    x227 = x0 * (3.0 * x223 + 3.0 * x224 + 3.0 * x35 + 3.0 * x40) + x17 * x226
    x228 = x222 * x92
    x229 = x114 * x220
    x230 = x229 * x23
    x231 = x228 + x230
    x232 = x0 * (x138 + x93)
    x233 = x125 * x17
    x234 = x223 + x224
    x235 = x0 * (x225 + 2.0 * x232 + 2.0 * x233) + x17 * x234
    x236 = x135 * x220
    x237 = x151 + 2.0 * x236
    x238 = x0 * (x221 + x237)
    x239 = x231 * x92
    x240 = x238 + x239
    x241 = x138 + x8
    x242 = x118 * x23 + x17 * x241
    x243 = x232 + x233
    x244 = x0 * (3.0 * x122 + 3.0 * x124 + x242) + x17 * x243
    x245 = x116 + x236
    x246 = x245 * x92
    x247 = x0 * (x135 + x229)
    x248 = 2.0 * x247
    x249 = 4.0 * x116 * x220
    x250 = x248 + x249
    x251 = x0 * (2.0 * x228 + 2.0 * x246 + x250)
    x252 = x240 * x92
    x253 = x251 + x252
    x254 = x0 * (3.0 * x138 + x27) + x17 * x242
    x255 = x144 * x222
    x256 = x230 + x255
    x257 = x144 * x231
    x258 = x238 + x257
    x259 = x144 * x240
    x260 = x251 + x259
    x261 = x246 + x247
    x262 = x261 * x92
    x263 = x0 * (x115 + x237)
    x264 = 2.0 * x263
    x265 = 3.0 * x238 + x264
    x266 = x0 * (3.0 * x239 + 2.0 * x262 + x265)
    x267 = x144 * x253
    x268 = x266 + x267
    x269 = x151 + x171 * x229
    x270 = x0 * (x221 + x269)
    x271 = x144 * x256
    x272 = x270 + x271
    x273 = x144 * x245
    x274 = 2.0 * x273
    x275 = x0 * (x228 + x250 + x255 + x274)
    x276 = x144 * x258
    x277 = x275 + x276
    x278 = x144 * x261
    x279 = 2.0 * x278
    x280 = 2.0 * x257
    x281 = x0 * (x239 + x265 + x279 + x280)
    x282 = x144 * x260
    x283 = x281 + x282
    x284 = 3.0 * x247
    x285 = x0 * (x137 + 3.0 * x246 + x284)
    x286 = x144 * (x262 + x263)
    x287 = 2.0 * x285 + 2.0 * x286
    x288 = x0 * (4.0 * x251 + x252 + 3.0 * x259 + x287)
    x289 = x144 * x268
    x290 = x288 + x289
    x291 = x0 * (x146 + x229)
    x292 = x146 * x220
    x293 = x144 * (x116 + x292)
    x294 = x0 * (x249 + 2.0 * x255 + 2.0 * x291 + 2.0 * x293) + x144 * x272
    x295 = x0 * (x147 + x151 + x236 + x292)
    x296 = x144 * (x247 + x273)
    x297 = 2.0 * x295 + 2.0 * x296
    x298 = x0 * (2.0 * x238 + x272 + x280 + x297)
    x299 = x144 * x277
    x300 = x298 + x299
    x301 = x0 * (x150 + x246 + x274 + x284)
    x302 = x144 * (x263 + x278)
    x303 = x0 * (
        2.0 * x251 + 2.0 * x259 + 2.0 * x275 + 2.0 * x276 + 2.0 * x301 + 2.0 * x302
    )
    x304 = x144 * x283
    x305 = x303 + x304
    x306 = x0 * (x154 + x262 + 4.0 * x263 + 3.0 * x278)
    x307 = x144 * (x285 + x286)
    x308 = 3.0 * x281 + 3.0 * x282
    x309 = x0 * (2.0 * x266 + 2.0 * x267 + 2.0 * x306 + 2.0 * x307 + x308) + x144 * x290
    x310 = x5 * x89
    x311 = x309 * x310
    x312 = x17 * x310
    x313 = x291 + x293
    x314 = (
        x0 * (x171 * x313 + x23 * (x166 + x269) + 3.0 * x270 + 3.0 * x271) + x144 * x294
    )
    x315 = (
        x0
        * (
            x171 * (x295 + x296)
            + x23 * (x170 + x248 + x274 + x313)
            + 3.0 * x275
            + 3.0 * x276
            + x294
        )
        + x144 * x300
    )
    x316 = x310 * (
        x0
        * (
            x171 * (x301 + x302)
            + x23 * (x175 + x264 + x279 + x297)
            + 2.0 * x298
            + 2.0 * x299
            + x308
        )
        + x144 * x305
    )
    x317 = x3 * x310
    x318 = -x112 - C[2]
    x319 = x129 * x318**2
    x320 = x132 + x319
    x321 = x113 * x320
    x322 = x129 * x318
    x323 = x23 * x322
    x324 = x321 + x323
    x325 = x140 * x318
    x326 = x162 + 2.0 * x325
    x327 = x0 * (x319 + x326)
    x328 = x113 * x324
    x329 = x327 + x328
    x330 = x132 + x325
    x331 = x113 * x330
    x332 = x0 * (x140 + x322)
    x333 = 2.0 * x332
    x334 = 4.0 * x132 * x318
    x335 = x333 + x334
    x336 = x0 * (2.0 * x321 + 2.0 * x331 + x335)
    x337 = x113 * x329
    x338 = x336 + x337
    x339 = x155 * x320
    x340 = x323 + x339
    x341 = x155 * x324
    x342 = x327 + x341
    x343 = x155 * x329
    x344 = x336 + x343
    x345 = x331 + x332
    x346 = x113 * x345
    x347 = x0 * (x131 + x326)
    x348 = 2.0 * x347
    x349 = 3.0 * x327 + x348
    x350 = x0 * (3.0 * x328 + 2.0 * x346 + x349)
    x351 = x155 * x338
    x352 = x350 + x351
    x353 = x162 + x185 * x322
    x354 = x0 * (x319 + x353)
    x355 = x155 * x340
    x356 = x354 + x355
    x357 = x155 * x330
    x358 = 2.0 * x357
    x359 = x0 * (x321 + x335 + x339 + x358)
    x360 = x155 * x342
    x361 = x359 + x360
    x362 = x155 * x345
    x363 = 2.0 * x362
    x364 = 2.0 * x341
    x365 = x0 * (x328 + x349 + x363 + x364)
    x366 = x155 * x344
    x367 = x365 + x366
    x368 = 3.0 * x332
    x369 = x0 * (x143 + 3.0 * x331 + x368)
    x370 = x155 * (x346 + x347)
    x371 = 2.0 * x369 + 2.0 * x370
    x372 = x0 * (4.0 * x336 + x337 + 3.0 * x343 + x371)
    x373 = x155 * x352
    x374 = x372 + x373
    x375 = 3.14159265358979 * x1 * x5 * x87
    x376 = x17 * x375
    x377 = x0 * (x157 + x322)
    x378 = x157 * x318
    x379 = x155 * (x132 + x378)
    x380 = x0 * (x334 + 2.0 * x339 + 2.0 * x377 + 2.0 * x379) + x155 * x356
    x381 = x0 * (x158 + x162 + x325 + x378)
    x382 = x155 * (x332 + x357)
    x383 = 2.0 * x381 + 2.0 * x382
    x384 = x0 * (2.0 * x327 + x356 + x364 + x383)
    x385 = x155 * x361
    x386 = x384 + x385
    x387 = x0 * (x161 + x331 + x358 + x368)
    x388 = x155 * (x347 + x362)
    x389 = x0 * (
        2.0 * x336 + 2.0 * x343 + 2.0 * x359 + 2.0 * x360 + 2.0 * x387 + 2.0 * x388
    )
    x390 = x155 * x367
    x391 = x389 + x390
    x392 = x0 * (x165 + x346 + 4.0 * x347 + 3.0 * x362)
    x393 = x155 * (x369 + x370)
    x394 = 3.0 * x365 + 3.0 * x366
    x395 = x0 * (2.0 * x350 + 2.0 * x351 + 2.0 * x392 + 2.0 * x393 + x394) + x155 * x374
    x396 = x375 * x395
    x397 = x3 * x375
    x398 = x377 + x379
    x399 = (
        x0 * (x185 * x398 + x23 * (x180 + x353) + 3.0 * x354 + 3.0 * x355) + x155 * x380
    )
    x400 = (
        x0
        * (
            x185 * (x381 + x382)
            + x23 * (x184 + x333 + x358 + x398)
            + 3.0 * x359
            + 3.0 * x360
            + x380
        )
        + x155 * x386
    )
    x401 = x375 * (
        x0
        * (
            x185 * (x387 + x388)
            + x23 * (x189 + x348 + x363 + x383)
            + 2.0 * x384
            + 2.0 * x385
            + x394
        )
        + x155 * x391
    )

    # 450 item(s)
    return numpy.array(
        [
            x90
            * (
                x0
                * (
                    x23 * (3.0 * x26 + 3.0 * x32 + x41 + x45)
                    + x48 * (x46 + x47)
                    + 3.0 * x63
                    + 3.0 * x69
                    + 3.0 * x76
                    + 3.0 * x83
                )
                + x17 * x86
            ),
            x111 * x92,
            x111 * x113,
            x117 * x128 * x129,
            x128 * x130 * x92,
            x114 * x128 * x133,
            x129 * x137 * x139,
            x117 * x139 * x140,
            x133 * x135 * x139,
            x114 * x139 * x143,
            x144 * x145,
            x110 * x129 * x148,
            x110 * x130 * x144,
            x127 * x129 * x150,
            x127 * x140 * x148,
            x127 * x133 * x146,
            x126 * x129 * x154,
            x126 * x140 * x150,
            x126 * x133 * x148,
            x126 * x143 * x146,
            x145 * x155,
            x110 * x156 * x92,
            x110 * x114 * x159,
            x117 * x127 * x157,
            x127 * x135 * x159,
            x114 * x127 * x161,
            x126 * x137 * x157,
            x117 * x126 * x159,
            x126 * x135 * x161,
            x114 * x126 * x165,
            x129 * x167 * x85,
            x129 * x170 * x82,
            x140 * x167 * x82,
            x108 * x129 * x175,
            x108 * x140 * x170,
            x108 * x133 * x167,
            x106 * x129 * x179,
            x106 * x140 * x175,
            x106 * x133 * x170,
            x106 * x143 * x167,
            x144 * x156 * x85,
            x148 * x157 * x82,
            x146 * x159 * x82,
            x108 * x150 * x157,
            x108 * x148 * x159,
            x108 * x146 * x161,
            x106 * x154 * x157,
            x106 * x150 * x159,
            x106 * x148 * x161,
            x106 * x146 * x165,
            x114 * x181 * x85,
            x135 * x181 * x82,
            x114 * x184 * x82,
            x108 * x117 * x181,
            x108 * x135 * x184,
            x108 * x114 * x189,
            x106 * x137 * x181,
            x106 * x117 * x184,
            x106 * x135 * x189,
            x106 * x114 * x193,
            x129 * x194 * x75,
            x129 * x197 * x80,
            x140 * x194 * x80,
            x129 * x200 * x67,
            x140 * x197 * x67,
            x133 * x194 * x67,
            x104 * x129 * x202,
            x104 * x140 * x200,
            x104 * x133 * x197,
            x104 * x143 * x194,
            x157 * x167 * x75,
            x157 * x170 * x80,
            x159 * x167 * x80,
            x157 * x175 * x67,
            x159 * x170 * x67,
            x161 * x167 * x67,
            x104 * x157 * x179,
            x104 * x159 * x175,
            x104 * x161 * x170,
            x104 * x165 * x167,
            x146 * x181 * x75,
            x148 * x181 * x80,
            x146 * x184 * x80,
            x150 * x181 * x67,
            x148 * x184 * x67,
            x146 * x189 * x67,
            x104 * x154 * x181,
            x104 * x150 * x184,
            x104 * x148 * x189,
            x104 * x146 * x193,
            x114 * x203 * x75,
            x135 * x203 * x80,
            x114 * x206 * x80,
            x117 * x203 * x67,
            x135 * x206 * x67,
            x114 * x209 * x67,
            x104 * x137 * x203,
            x104 * x117 * x206,
            x104 * x135 * x209,
            x104 * x114 * x211,
            x129 * x212 * x73,
            x129 * x213 * x56,
            x140 * x212 * x56,
            x129 * x214 * x54,
            x140 * x213 * x54,
            x133 * x212 * x54,
            x129 * x215 * x51,
            x140 * x214 * x51,
            x133 * x213 * x51,
            x143 * x212 * x51,
            x157 * x194 * x73,
            x157 * x197 * x56,
            x159 * x194 * x56,
            x157 * x200 * x54,
            x159 * x197 * x54,
            x161 * x194 * x54,
            x157 * x202 * x51,
            x159 * x200 * x51,
            x161 * x197 * x51,
            x165 * x194 * x51,
            x167 * x181 * x73,
            x170 * x181 * x56,
            x167 * x184 * x56,
            x175 * x181 * x54,
            x170 * x184 * x54,
            x167 * x189 * x54,
            x179 * x181 * x51,
            x175 * x184 * x51,
            x170 * x189 * x51,
            x167 * x193 * x51,
            x146 * x203 * x73,
            x148 * x203 * x56,
            x146 * x206 * x56,
            x150 * x203 * x54,
            x148 * x206 * x54,
            x146 * x209 * x54,
            x154 * x203 * x51,
            x150 * x206 * x51,
            x148 * x209 * x51,
            x146 * x211 * x51,
            x114 * x216 * x73,
            x135 * x216 * x56,
            x114 * x217 * x56,
            x117 * x216 * x54,
            x135 * x217 * x54,
            x114 * x218 * x54,
            x137 * x216 * x51,
            x117 * x217 * x51,
            x135 * x218 * x51,
            x114 * x219 * x51,
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
            x129 * x226 * x256,
            x129 * x234 * x258,
            x140 * x234 * x256,
            x129 * x243 * x260,
            x140 * x243 * x258,
            x133 * x243 * x256,
            x129 * x242 * x268,
            x140 * x242 * x260,
            x133 * x242 * x258,
            x143 * x242 * x256,
            x157 * x222 * x226,
            x157 * x231 * x234,
            x159 * x222 * x234,
            x157 * x240 * x243,
            x159 * x231 * x243,
            x161 * x222 * x243,
            x157 * x242 * x253,
            x159 * x240 * x242,
            x161 * x231 * x242,
            x165 * x222 * x242,
            x129 * x272 * x41,
            x129 * x277 * x96,
            x140 * x272 * x96,
            x125 * x129 * x283,
            x125 * x140 * x277,
            x125 * x133 * x272,
            x129 * x241 * x290,
            x140 * x241 * x283,
            x133 * x241 * x277,
            x143 * x241 * x272,
            x157 * x256 * x41,
            x157 * x258 * x96,
            x159 * x256 * x96,
            x125 * x157 * x260,
            x125 * x159 * x258,
            x125 * x161 * x256,
            x157 * x241 * x268,
            x159 * x241 * x260,
            x161 * x241 * x258,
            x165 * x241 * x256,
            x181 * x222 * x41,
            x181 * x231 * x96,
            x184 * x222 * x96,
            x125 * x181 * x240,
            x125 * x184 * x231,
            x125 * x189 * x222,
            x181 * x241 * x253,
            x184 * x240 * x241,
            x189 * x231 * x241,
            x193 * x222 * x241,
            x129 * x294 * x39,
            x129 * x25 * x300,
            x140 * x25 * x294,
            x123 * x129 * x305,
            x123 * x140 * x300,
            x123 * x133 * x294,
            x17 * x311,
            x113 * x305 * x312,
            x118 * x133 * x300,
            x118 * x143 * x294,
            x157 * x272 * x39,
            x157 * x25 * x277,
            x159 * x25 * x272,
            x123 * x157 * x283,
            x123 * x159 * x277,
            x123 * x161 * x272,
            x155 * x290 * x312,
            x118 * x159 * x283,
            x118 * x161 * x277,
            x118 * x165 * x272,
            x181 * x256 * x39,
            x181 * x25 * x258,
            x184 * x25 * x256,
            x123 * x181 * x260,
            x123 * x184 * x258,
            x123 * x189 * x256,
            x118 * x181 * x268,
            x118 * x184 * x260,
            x118 * x189 * x258,
            x118 * x193 * x256,
            x203 * x222 * x39,
            x203 * x231 * x25,
            x206 * x222 * x25,
            x123 * x203 * x240,
            x123 * x206 * x231,
            x123 * x209 * x222,
            x118 * x203 * x253,
            x118 * x206 * x240,
            x118 * x209 * x231,
            x118 * x211 * x222,
            x129 * x314 * x37,
            x129 * x21 * x315,
            x140 * x21 * x314,
            x3 * x316,
            x113 * x315 * x317,
            x10 * x133 * x314,
            x310
            * (
                x0
                * (
                    x171 * (x306 + x307)
                    + x23 * (x179 + x287 + 3.0 * x301 + 3.0 * x302)
                    + 3.0 * x288
                    + 3.0 * x289
                    + 3.0 * x303
                    + 3.0 * x304
                )
                + x144 * x309
            ),
            x113 * x316,
            x133 * x315 * x7,
            x143 * x314 * x7,
            x157 * x294 * x37,
            x157 * x21 * x300,
            x159 * x21 * x294,
            x155 * x305 * x317,
            x10 * x159 * x300,
            x10 * x161 * x294,
            x155 * x311,
            x159 * x305 * x7,
            x161 * x300 * x7,
            x165 * x294 * x7,
            x181 * x272 * x37,
            x181 * x21 * x277,
            x184 * x21 * x272,
            x10 * x181 * x283,
            x10 * x184 * x277,
            x10 * x189 * x272,
            x181 * x290 * x7,
            x184 * x283 * x7,
            x189 * x277 * x7,
            x193 * x272 * x7,
            x203 * x256 * x37,
            x203 * x21 * x258,
            x206 * x21 * x256,
            x10 * x203 * x260,
            x10 * x206 * x258,
            x10 * x209 * x256,
            x203 * x268 * x7,
            x206 * x260 * x7,
            x209 * x258 * x7,
            x211 * x256 * x7,
            x216 * x222 * x37,
            x21 * x216 * x231,
            x21 * x217 * x222,
            x10 * x216 * x240,
            x10 * x217 * x231,
            x10 * x218 * x222,
            x216 * x253 * x7,
            x217 * x240 * x7,
            x218 * x231 * x7,
            x219 * x222 * x7,
            x114 * x227 * x320,
            x135 * x235 * x320,
            x114 * x235 * x324,
            x117 * x244 * x320,
            x135 * x244 * x324,
            x114 * x244 * x329,
            x137 * x254 * x320,
            x117 * x254 * x324,
            x135 * x254 * x329,
            x114 * x254 * x338,
            x146 * x226 * x320,
            x148 * x234 * x320,
            x146 * x234 * x324,
            x150 * x243 * x320,
            x148 * x243 * x324,
            x146 * x243 * x329,
            x154 * x242 * x320,
            x150 * x242 * x324,
            x148 * x242 * x329,
            x146 * x242 * x338,
            x114 * x226 * x340,
            x135 * x234 * x340,
            x114 * x234 * x342,
            x117 * x243 * x340,
            x135 * x243 * x342,
            x114 * x243 * x344,
            x137 * x242 * x340,
            x117 * x242 * x342,
            x135 * x242 * x344,
            x114 * x242 * x352,
            x167 * x320 * x41,
            x170 * x320 * x96,
            x167 * x324 * x96,
            x125 * x175 * x320,
            x125 * x170 * x324,
            x125 * x167 * x329,
            x179 * x241 * x320,
            x175 * x241 * x324,
            x170 * x241 * x329,
            x167 * x241 * x338,
            x146 * x340 * x41,
            x148 * x340 * x96,
            x146 * x342 * x96,
            x125 * x150 * x340,
            x125 * x148 * x342,
            x125 * x146 * x344,
            x154 * x241 * x340,
            x150 * x241 * x342,
            x148 * x241 * x344,
            x146 * x241 * x352,
            x114 * x356 * x41,
            x135 * x356 * x96,
            x114 * x361 * x96,
            x117 * x125 * x356,
            x125 * x135 * x361,
            x114 * x125 * x367,
            x137 * x241 * x356,
            x117 * x241 * x361,
            x135 * x241 * x367,
            x114 * x241 * x374,
            x194 * x320 * x39,
            x197 * x25 * x320,
            x194 * x25 * x324,
            x123 * x200 * x320,
            x123 * x197 * x324,
            x123 * x194 * x329,
            x118 * x202 * x320,
            x118 * x200 * x324,
            x118 * x197 * x329,
            x118 * x194 * x338,
            x167 * x340 * x39,
            x170 * x25 * x340,
            x167 * x25 * x342,
            x123 * x175 * x340,
            x123 * x170 * x342,
            x123 * x167 * x344,
            x118 * x179 * x340,
            x118 * x175 * x342,
            x118 * x170 * x344,
            x118 * x167 * x352,
            x146 * x356 * x39,
            x148 * x25 * x356,
            x146 * x25 * x361,
            x123 * x150 * x356,
            x123 * x148 * x361,
            x123 * x146 * x367,
            x118 * x154 * x356,
            x118 * x150 * x361,
            x118 * x148 * x367,
            x144 * x374 * x376,
            x114 * x380 * x39,
            x135 * x25 * x380,
            x114 * x25 * x386,
            x117 * x123 * x380,
            x123 * x135 * x386,
            x114 * x123 * x391,
            x118 * x137 * x380,
            x117 * x118 * x386,
            x376 * x391 * x92,
            x17 * x396,
            x212 * x320 * x37,
            x21 * x213 * x320,
            x21 * x212 * x324,
            x10 * x214 * x320,
            x10 * x213 * x324,
            x10 * x212 * x329,
            x215 * x320 * x7,
            x214 * x324 * x7,
            x213 * x329 * x7,
            x212 * x338 * x7,
            x194 * x340 * x37,
            x197 * x21 * x340,
            x194 * x21 * x342,
            x10 * x200 * x340,
            x10 * x197 * x342,
            x10 * x194 * x344,
            x202 * x340 * x7,
            x200 * x342 * x7,
            x197 * x344 * x7,
            x194 * x352 * x7,
            x167 * x356 * x37,
            x170 * x21 * x356,
            x167 * x21 * x361,
            x10 * x175 * x356,
            x10 * x170 * x361,
            x10 * x167 * x367,
            x179 * x356 * x7,
            x175 * x361 * x7,
            x170 * x367 * x7,
            x167 * x374 * x7,
            x146 * x37 * x380,
            x148 * x21 * x380,
            x146 * x21 * x386,
            x10 * x150 * x380,
            x10 * x148 * x386,
            x144 * x391 * x397,
            x154 * x380 * x7,
            x150 * x386 * x7,
            x148 * x391 * x7,
            x144 * x396,
            x114 * x37 * x399,
            x135 * x21 * x399,
            x114 * x21 * x400,
            x10 * x117 * x399,
            x397 * x400 * x92,
            x3 * x401,
            x137 * x399 * x7,
            x117 * x400 * x7,
            x401 * x92,
            x375
            * (
                x0
                * (
                    x185 * (x392 + x393)
                    + x23 * (x193 + x371 + 3.0 * x387 + 3.0 * x388)
                    + 3.0 * x372
                    + 3.0 * x373
                    + 3.0 * x389
                    + 3.0 * x390
                )
                + x155 * x395
            ),
        ]
    )


def diag_quadrupole3d_44(a, A, b, B, C):
    """Cartesian 3D (gg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - B[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x3 * x7
    x9 = -x2 - C[0]
    x10 = x7 * x9
    x11 = x0 * (x10 + x8)
    x12 = x0 * x7
    x13 = x8 * x9
    x14 = x12 + x13
    x15 = x14 * x3
    x16 = x11 + x15
    x17 = x16 * x3
    x18 = x3**2 * x7
    x19 = 3.0 * x12
    x20 = 2.0 * x13 + x19
    x21 = x0 * (x18 + x20)
    x22 = 4.0 * x21
    x23 = -x2 - A[0]
    x24 = x16 * x23
    x25 = x0 * (3.0 * x18 + x19)
    x26 = x12 + x18
    x27 = x26 * x3
    x28 = 2.0 * x0
    x29 = x28 * x8
    x30 = x27 + x29
    x31 = x23 * x30
    x32 = x25 + x31
    x33 = x0 * (x17 + x22 + 3.0 * x24 + x32)
    x34 = 3.0 * x11
    x35 = x0 * (3.0 * x15 + x30 + x34)
    x36 = x17 + x21
    x37 = x23 * x36
    x38 = x23 * (x35 + x37)
    x39 = x3 * x30
    x40 = x0 * (5.0 * x25 + 4.0 * x31 + x39)
    x41 = x12 * x3
    x42 = 8.0 * x41
    x43 = x0 * (4.0 * x27 + x42)
    x44 = x25 + x39
    x45 = x23 * x44
    x46 = x43 + x45
    x47 = x23 * x46
    x48 = x40 + x47
    x49 = x0 * (4.0 * x17 + x22 + x44)
    x50 = x3 * x36
    x51 = x23 * (x35 + x50)
    x52 = 2.0 * x49 + 2.0 * x51
    x53 = x0 * (5.0 * x35 + 4.0 * x37 + x46 + x50)
    x54 = x23 * (x49 + x51)
    x55 = 2.0 * x23
    x56 = x7 * x9**2
    x57 = x12 + x56
    x58 = x3 * x57
    x59 = 2.0 * x11
    x60 = 4.0 * x12 * x9
    x61 = x59 + x60
    x62 = x0 * (2.0 * x15 + 2.0 * x58 + x61)
    x63 = x0 * (x20 + x56)
    x64 = x10 * x28
    x65 = x58 + x64
    x66 = x3 * x65
    x67 = x63 + x66
    x68 = x3 * x67
    x69 = x62 + x68
    x70 = x3 * x69
    x71 = x23 * x69
    x72 = 2.0 * x21
    x73 = 3.0 * x63 + x72
    x74 = x0 * (2.0 * x17 + 3.0 * x66 + x73)
    x75 = x0 * (x52 + x70 + 4.0 * x71 + 5.0 * x74)
    x76 = 2.0 * x35
    x77 = 4.0 * x62 + x76
    x78 = x0 * (2.0 * x50 + 4.0 * x68 + x77)
    x79 = x70 + x74
    x80 = x23 * x79
    x81 = x78 + x80
    x82 = x23 * x81
    x83 = 2.0 * x24
    x84 = x23 * x65
    x85 = 2.0 * x84
    x86 = x0 * (x66 + x73 + x83 + x85)
    x87 = x23 * x67
    x88 = x62 + x87
    x89 = x23 * x88
    x90 = 3.0 * x86 + 3.0 * x89
    x91 = x0 * (2.0 * x33 + 2.0 * x38 + 2.0 * x71 + 2.0 * x74 + x90)
    x92 = 2.0 * x37
    x93 = x0 * (x68 + x77 + 3.0 * x87 + x92)
    x94 = x71 + x74
    x95 = x23 * x94
    x96 = x93 + x95
    x97 = x23 * x96
    x98 = x75 + x82
    x99 = (
        x0 * (2.0 * x53 + 2.0 * x54 + 2.0 * x78 + 2.0 * x80 + 4.0 * x93 + 4.0 * x95)
        + x23 * x98
    )
    x100 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x101 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x102 = 3.14159265358979 * x1 * x101
    x103 = x100 * x102
    x104 = -x1 * (a * A[1] + b * B[1])
    x105 = -x104 - B[1]
    x106 = x14 * x23
    x107 = 2.0 * x106
    x108 = x23 * x26
    x109 = x108 + x29
    x110 = x0 * (x107 + x109 + x15 + x34)
    x111 = x23 * (x21 + x24)
    x112 = x0 * (3.0 * x108 + x27 + x42)
    x113 = x23 * x32
    x114 = x112 + x113
    x115 = x23 * x57
    x116 = x0 * (x107 + x115 + x58 + x61)
    x117 = x63 + x84
    x118 = x117 * x23
    x119 = x0 * (
        2.0 * x110 + 2.0 * x111 + 2.0 * x116 + 2.0 * x118 + 2.0 * x62 + 2.0 * x87
    )
    x120 = x86 + x89
    x121 = x120 * x23
    x122 = x91 + x97
    x123 = x103 * (
        x0
        * (
            3.0 * x119
            + 3.0 * x121
            + x28 * (3.0 * x110 + 3.0 * x111 + x114 + x76 + x92)
            + x55 * (x33 + x38)
            + 3.0 * x93
            + 3.0 * x95
        )
        + x122 * x23
    )
    x124 = -x1 * (a * A[2] + b * B[2])
    x125 = -x124 - B[2]
    x126 = x100 * x6
    x127 = x105**2 * x126
    x128 = x0 * x126
    x129 = x127 + x128
    x130 = x19 + x55 * x8
    x131 = x0 * (x130 + x18)
    x132 = x109 * x23
    x133 = x131 + x132
    x134 = x23 * x8
    x135 = x10 * x23
    x136 = x0 * (x13 + x134 + x135 + x19)
    x137 = x23 * (x106 + x11)
    x138 = 2.0 * x136 + 2.0 * x137
    x139 = x10 * x55 + x19
    x140 = x0 * (x139 + x56)
    x141 = x115 + x64
    x142 = x141 * x23
    x143 = x140 + x142
    x144 = x0 * (x138 + x143 + 2.0 * x63 + x85)
    x145 = x116 + x118
    x146 = x145 * x23
    x147 = x119 + x121
    x148 = (
        x0
        * (
            2.0 * x144
            + 2.0 * x146
            + x28 * (x133 + x138 + x72 + x83)
            + x55 * (x110 + x111)
            + x90
        )
        + x147 * x23
    )
    x149 = x101 * x6
    x150 = x103 * x125
    x151 = x125**2 * x149
    x152 = x0 * x149
    x153 = x151 + x152
    x154 = x105 * x129
    x155 = x105 * x126
    x156 = x155 * x28
    x157 = x154 + x156
    x158 = x23 * x7
    x159 = x0 * (x10 + x158)
    x160 = x23 * (x12 + x135)
    x161 = x159 + x160
    x162 = x0 * (x158 + x8)
    x163 = x12 + x134
    x164 = x163 * x23
    x165 = x162 + x164
    x166 = x0 * (2.0 * x115 + 2.0 * x159 + 2.0 * x160 + x60) + x143 * x23
    x167 = x144 + x146
    x168 = (
        x0
        * (
            3.0 * x116
            + 3.0 * x118
            + x166
            + x28 * (x107 + x161 + x165 + x59)
            + x55 * (x136 + x137)
        )
        + x167 * x23
    )
    x169 = x125 * x149
    x170 = x125 * x153
    x171 = x169 * x28
    x172 = x170 + x171
    x173 = 3.0 * x128
    x174 = x0 * (3.0 * x127 + x173)
    x175 = x105 * x157
    x176 = x174 + x175
    x177 = x23**2 * x7
    x178 = x0 * (3.0 * x140 + 3.0 * x142 + x161 * x55 + x28 * (x139 + x177)) + x166 * x23
    x179 = 3.0 * x152
    x180 = x0 * (3.0 * x151 + x179)
    x181 = x125 * x172
    x182 = x180 + x181
    x183 = -x104 - A[1]
    x184 = x103 * x99
    x185 = x126 * x183
    x186 = x105 * x185
    x187 = x128 + x186
    x188 = x129 * x183
    x189 = x156 + x188
    x190 = x157 * x183
    x191 = x174 + x190
    x192 = x105 * x128
    x193 = 8.0 * x192
    x194 = x0 * (4.0 * x154 + x193)
    x195 = x176 * x183
    x196 = x194 + x195
    x197 = -x124 - A[2]
    x198 = x103 * x197
    x199 = x149 * x197
    x200 = x125 * x199
    x201 = x152 + x200
    x202 = x153 * x197
    x203 = x171 + x202
    x204 = x172 * x197
    x205 = x180 + x204
    x206 = x125 * x152
    x207 = 8.0 * x206
    x208 = x0 * (4.0 * x170 + x207)
    x209 = x182 * x197
    x210 = x208 + x209
    x211 = x126 * x183**2
    x212 = x128 + x211
    x213 = x0 * (x155 + x185)
    x214 = x183 * x187
    x215 = x213 + x214
    x216 = 2.0 * x183
    x217 = x155 * x216 + x173
    x218 = x0 * (x127 + x217)
    x219 = x183 * x189
    x220 = x218 + x219
    x221 = x0 * (x154 + 3.0 * x188 + x193)
    x222 = x183 * x191
    x223 = x221 + x222
    x224 = x0 * (5.0 * x174 + x175 + 4.0 * x190)
    x225 = x183 * x196
    x226 = x224 + x225
    x227 = x149 * x197**2
    x228 = x152 + x227
    x229 = x0 * (x169 + x199)
    x230 = x197 * x201
    x231 = x229 + x230
    x232 = 2.0 * x197
    x233 = x169 * x232 + x179
    x234 = x0 * (x151 + x233)
    x235 = x197 * x203
    x236 = x234 + x235
    x237 = x0 * (x170 + 3.0 * x202 + x207)
    x238 = x197 * x205
    x239 = x237 + x238
    x240 = x0 * (5.0 * x180 + x181 + 4.0 * x204)
    x241 = x197 * x210
    x242 = x240 + x241
    x243 = x183 * x212 + x185 * x28
    x244 = x0 * (x211 + x217)
    x245 = x183 * x215
    x246 = x244 + x245
    x247 = x0 * (2.0 * x188 + 4.0 * x192 + 2.0 * x213 + 2.0 * x214)
    x248 = x183 * x220
    x249 = x247 + x248
    x250 = 3.0 * x218 + 3.0 * x219
    x251 = x0 * (2.0 * x174 + 2.0 * x190 + x250)
    x252 = x183 * x223
    x253 = x251 + x252
    x254 = x0 * (2.0 * x194 + 2.0 * x195 + 4.0 * x221 + 4.0 * x222) + x183 * x226
    x255 = x197 * x228 + x199 * x28
    x256 = x0 * (x227 + x233)
    x257 = x197 * x231
    x258 = x256 + x257
    x259 = x0 * (2.0 * x202 + 4.0 * x206 + 2.0 * x229 + 2.0 * x230)
    x260 = x197 * x236
    x261 = x259 + x260
    x262 = 3.0 * x234 + 3.0 * x235
    x263 = x0 * (2.0 * x180 + 2.0 * x204 + x262)
    x264 = x197 * x239
    x265 = x263 + x264
    x266 = x0 * (2.0 * x208 + 2.0 * x209 + 4.0 * x237 + 4.0 * x238) + x197 * x242
    x267 = x0 * (x173 + 3.0 * x211) + x183 * x243
    x268 = x0 * (3.0 * x213 + 3.0 * x214 + x243) + x183 * x246
    x269 = x0 * (2.0 * x244 + 2.0 * x245 + x250) + x183 * x249
    x270 = x0 * (3.0 * x221 + 3.0 * x222 + 3.0 * x247 + 3.0 * x248) + x183 * x253
    x271 = x0 * (3.0 * x224 + 3.0 * x225 + 4.0 * x251 + 4.0 * x252) + x183 * x254
    x272 = x0 * (x179 + 3.0 * x227) + x197 * x255
    x273 = x0 * (3.0 * x229 + 3.0 * x230 + x255) + x197 * x258
    x274 = x0 * (2.0 * x256 + 2.0 * x257 + x262) + x197 * x261
    x275 = x0 * (3.0 * x237 + 3.0 * x238 + 3.0 * x259 + 3.0 * x260) + x197 * x265
    x276 = x0 * (3.0 * x240 + 3.0 * x241 + 4.0 * x263 + 4.0 * x264) + x197 * x266
    x277 = -x104 - C[1]
    x278 = x126 * x277**2
    x279 = x128 + x278
    x280 = 3.0 * x131 + 3.0 * x132
    x281 = x0 * (2.0 * x25 + x280 + 2.0 * x31)
    x282 = x114 * x23
    x283 = x0 * (4.0 * x112 + 4.0 * x113 + 2.0 * x43 + 2.0 * x45) + x23 * x48
    x284 = x0 * (4.0 * x281 + 4.0 * x282 + 3.0 * x40 + 3.0 * x47) + x23 * x283
    x285 = x105 * x279
    x286 = x126 * x277
    x287 = x28 * x286
    x288 = x285 + x287
    x289 = x0 * (2.0 * x108 + 2.0 * x162 + 2.0 * x164 + 4.0 * x41)
    x290 = x133 * x23
    x291 = x281 + x282
    x292 = x0 * (3.0 * x112 + 3.0 * x113 + 3.0 * x289 + 3.0 * x290) + x23 * x291
    x293 = x155 * x277
    x294 = x173 + 2.0 * x293
    x295 = x0 * (x278 + x294)
    x296 = x105 * x288
    x297 = x295 + x296
    x298 = x0 * (x130 + x177)
    x299 = x165 * x23
    x300 = x289 + x290
    x301 = x0 * (x280 + 2.0 * x298 + 2.0 * x299) + x23 * x300
    x302 = x128 + x293
    x303 = x105 * x302
    x304 = x0 * (x155 + x286)
    x305 = 2.0 * x304
    x306 = 4.0 * x128 * x277
    x307 = x305 + x306
    x308 = x0 * (2.0 * x285 + 2.0 * x303 + x307)
    x309 = x105 * x297
    x310 = x308 + x309
    x311 = x12 + x177
    x312 = x158 * x28 + x23 * x311
    x313 = x298 + x299
    x314 = x0 * (3.0 * x162 + 3.0 * x164 + x312) + x23 * x313
    x315 = x303 + x304
    x316 = x105 * x315
    x317 = x0 * (x127 + x294)
    x318 = 2.0 * x317
    x319 = 3.0 * x295 + x318
    x320 = x0 * (3.0 * x296 + 2.0 * x316 + x319)
    x321 = x105 * x310
    x322 = x320 + x321
    x323 = x0 * (3.0 * x177 + x19) + x23 * x312
    x324 = x183 * x279
    x325 = x287 + x324
    x326 = x183 * x288
    x327 = x295 + x326
    x328 = x183 * x297
    x329 = x308 + x328
    x330 = x183 * x310
    x331 = x320 + x330
    x332 = x316 + x317
    x333 = x105 * x332
    x334 = 3.0 * x304
    x335 = x0 * (x157 + 3.0 * x303 + x334)
    x336 = 2.0 * x335
    x337 = 4.0 * x308 + x336
    x338 = x0 * (4.0 * x309 + 2.0 * x333 + x337)
    x339 = x183 * x322
    x340 = x338 + x339
    x341 = x173 + x216 * x286
    x342 = x0 * (x278 + x341)
    x343 = x183 * x325
    x344 = x342 + x343
    x345 = x183 * x302
    x346 = 2.0 * x345
    x347 = x0 * (x285 + x307 + x324 + x346)
    x348 = x183 * x327
    x349 = x347 + x348
    x350 = x183 * x315
    x351 = 2.0 * x350
    x352 = 2.0 * x326
    x353 = x0 * (x296 + x319 + x351 + x352)
    x354 = x183 * x329
    x355 = x353 + x354
    x356 = x183 * x332
    x357 = 2.0 * x356
    x358 = x0 * (x309 + 3.0 * x328 + x337 + x357)
    x359 = x183 * x331
    x360 = x358 + x359
    x361 = 4.0 * x317
    x362 = x0 * (x176 + 4.0 * x316 + x361)
    x363 = x183 * (x333 + x335)
    x364 = 2.0 * x362 + 2.0 * x363
    x365 = x0 * (5.0 * x320 + x321 + 4.0 * x330 + x364)
    x366 = x183 * x340
    x367 = x365 + x366
    x368 = x0 * (x185 + x286)
    x369 = x185 * x277
    x370 = x183 * (x128 + x369)
    x371 = x0 * (x306 + 2.0 * x324 + 2.0 * x368 + 2.0 * x370) + x183 * x344
    x372 = x0 * (x173 + x186 + x293 + x369)
    x373 = x183 * (x304 + x345)
    x374 = 2.0 * x372 + 2.0 * x373
    x375 = x0 * (2.0 * x295 + x344 + x352 + x374)
    x376 = x183 * x349
    x377 = x375 + x376
    x378 = x0 * (x189 + x303 + x334 + x346)
    x379 = x183 * (x317 + x350)
    x380 = x0 * (
        2.0 * x308 + 2.0 * x328 + 2.0 * x347 + 2.0 * x348 + 2.0 * x378 + 2.0 * x379
    )
    x381 = x183 * x355
    x382 = x380 + x381
    x383 = x0 * (x191 + x316 + 3.0 * x350 + x361)
    x384 = x183 * (x335 + x356)
    x385 = 3.0 * x353 + 3.0 * x354
    x386 = x0 * (2.0 * x320 + 2.0 * x330 + 2.0 * x383 + 2.0 * x384 + x385)
    x387 = x183 * x360
    x388 = x386 + x387
    x389 = x0 * (x196 + x333 + 5.0 * x335 + 4.0 * x356)
    x390 = x183 * (x362 + x363)
    x391 = (
        x0 * (2.0 * x338 + 2.0 * x339 + 4.0 * x358 + 4.0 * x359 + 2.0 * x389 + 2.0 * x390)
        + x183 * x367
    )
    x392 = x102 * x5
    x393 = x391 * x392
    x394 = x23 * x392
    x395 = x368 + x370
    x396 = (
        x0 * (x216 * x395 + x28 * (x211 + x341) + 3.0 * x342 + 3.0 * x343) + x183 * x371
    )
    x397 = (
        x0
        * (
            x216 * (x372 + x373)
            + x28 * (x215 + x305 + x346 + x395)
            + 3.0 * x347
            + 3.0 * x348
            + x371
        )
        + x183 * x377
    )
    x398 = (
        x0
        * (
            x216 * (x378 + x379)
            + x28 * (x220 + x318 + x351 + x374)
            + 2.0 * x375
            + 2.0 * x376
            + x385
        )
        + x183 * x382
    )
    x399 = x392 * (
        x0
        * (
            x216 * (x383 + x384)
            + x28 * (x223 + x336 + x357 + 3.0 * x378 + 3.0 * x379)
            + 3.0 * x358
            + 3.0 * x359
            + 3.0 * x380
            + 3.0 * x381
        )
        + x183 * x388
    )
    x400 = x3 * x392
    x401 = -x124 - C[2]
    x402 = x149 * x401**2
    x403 = x152 + x402
    x404 = x125 * x403
    x405 = x149 * x401
    x406 = x28 * x405
    x407 = x404 + x406
    x408 = x169 * x401
    x409 = x179 + 2.0 * x408
    x410 = x0 * (x402 + x409)
    x411 = x125 * x407
    x412 = x410 + x411
    x413 = x152 + x408
    x414 = x125 * x413
    x415 = x0 * (x169 + x405)
    x416 = 2.0 * x415
    x417 = 4.0 * x152 * x401
    x418 = x416 + x417
    x419 = x0 * (2.0 * x404 + 2.0 * x414 + x418)
    x420 = x125 * x412
    x421 = x419 + x420
    x422 = x414 + x415
    x423 = x125 * x422
    x424 = x0 * (x151 + x409)
    x425 = 2.0 * x424
    x426 = 3.0 * x410 + x425
    x427 = x0 * (3.0 * x411 + 2.0 * x423 + x426)
    x428 = x125 * x421
    x429 = x427 + x428
    x430 = x197 * x403
    x431 = x406 + x430
    x432 = x197 * x407
    x433 = x410 + x432
    x434 = x197 * x412
    x435 = x419 + x434
    x436 = x197 * x421
    x437 = x427 + x436
    x438 = x423 + x424
    x439 = x125 * x438
    x440 = 3.0 * x415
    x441 = x0 * (x172 + 3.0 * x414 + x440)
    x442 = 2.0 * x441
    x443 = 4.0 * x419 + x442
    x444 = x0 * (4.0 * x420 + 2.0 * x439 + x443)
    x445 = x197 * x429
    x446 = x444 + x445
    x447 = x179 + x232 * x405
    x448 = x0 * (x402 + x447)
    x449 = x197 * x431
    x450 = x448 + x449
    x451 = x197 * x413
    x452 = 2.0 * x451
    x453 = x0 * (x404 + x418 + x430 + x452)
    x454 = x197 * x433
    x455 = x453 + x454
    x456 = x197 * x422
    x457 = 2.0 * x456
    x458 = 2.0 * x432
    x459 = x0 * (x411 + x426 + x457 + x458)
    x460 = x197 * x435
    x461 = x459 + x460
    x462 = x197 * x438
    x463 = 2.0 * x462
    x464 = x0 * (x420 + 3.0 * x434 + x443 + x463)
    x465 = x197 * x437
    x466 = x464 + x465
    x467 = 4.0 * x424
    x468 = x0 * (x182 + 4.0 * x423 + x467)
    x469 = x197 * (x439 + x441)
    x470 = 2.0 * x468 + 2.0 * x469
    x471 = x0 * (5.0 * x427 + x428 + 4.0 * x436 + x470)
    x472 = x197 * x446
    x473 = x471 + x472
    x474 = 3.14159265358979 * x1 * x100 * x5
    x475 = x23 * x474
    x476 = x0 * (x199 + x405)
    x477 = x199 * x401
    x478 = x197 * (x152 + x477)
    x479 = x0 * (x417 + 2.0 * x430 + 2.0 * x476 + 2.0 * x478) + x197 * x450
    x480 = x0 * (x179 + x200 + x408 + x477)
    x481 = x197 * (x415 + x451)
    x482 = 2.0 * x480 + 2.0 * x481
    x483 = x0 * (2.0 * x410 + x450 + x458 + x482)
    x484 = x197 * x455
    x485 = x483 + x484
    x486 = x0 * (x203 + x414 + x440 + x452)
    x487 = x197 * (x424 + x456)
    x488 = x0 * (
        2.0 * x419 + 2.0 * x434 + 2.0 * x453 + 2.0 * x454 + 2.0 * x486 + 2.0 * x487
    )
    x489 = x197 * x461
    x490 = x488 + x489
    x491 = x0 * (x205 + x423 + 3.0 * x456 + x467)
    x492 = x197 * (x441 + x462)
    x493 = 3.0 * x459 + 3.0 * x460
    x494 = x0 * (2.0 * x427 + 2.0 * x436 + 2.0 * x491 + 2.0 * x492 + x493)
    x495 = x197 * x466
    x496 = x494 + x495
    x497 = x0 * (x210 + x439 + 5.0 * x441 + 4.0 * x462)
    x498 = x197 * (x468 + x469)
    x499 = (
        x0 * (2.0 * x444 + 2.0 * x445 + 4.0 * x464 + 4.0 * x465 + 2.0 * x497 + 2.0 * x498)
        + x197 * x473
    )
    x500 = x474 * x499
    x501 = x3 * x474
    x502 = x476 + x478
    x503 = (
        x0 * (x232 * x502 + x28 * (x227 + x447) + 3.0 * x448 + 3.0 * x449) + x197 * x479
    )
    x504 = (
        x0
        * (
            x232 * (x480 + x481)
            + x28 * (x231 + x416 + x452 + x502)
            + 3.0 * x453
            + 3.0 * x454
            + x479
        )
        + x197 * x485
    )
    x505 = (
        x0
        * (
            x232 * (x486 + x487)
            + x28 * (x236 + x425 + x457 + x482)
            + 2.0 * x483
            + 2.0 * x484
            + x493
        )
        + x197 * x490
    )
    x506 = x474 * (
        x0
        * (
            x232 * (x491 + x492)
            + x28 * (x239 + x442 + x463 + 3.0 * x486 + 3.0 * x487)
            + 3.0 * x464
            + 3.0 * x465
            + 3.0 * x488
            + 3.0 * x489
        )
        + x197 * x496
    )

    # 675 item(s)
    return numpy.array(
        [
            x103
            * (
                x0
                * (
                    x28 * (4.0 * x33 + 4.0 * x38 + x48 + x52)
                    + x55 * (x53 + x54)
                    + 3.0 * x75
                    + 3.0 * x82
                    + 4.0 * x91
                    + 4.0 * x97
                )
                + x23 * x99
            ),
            x105 * x123,
            x123 * x125,
            x129 * x148 * x149,
            x105 * x148 * x150,
            x126 * x148 * x153,
            x149 * x157 * x168,
            x129 * x168 * x169,
            x153 * x155 * x168,
            x126 * x168 * x172,
            x149 * x176 * x178,
            x157 * x169 * x178,
            x129 * x153 * x178,
            x155 * x172 * x178,
            x126 * x178 * x182,
            x183 * x184,
            x122 * x149 * x187,
            x122 * x150 * x183,
            x147 * x149 * x189,
            x147 * x169 * x187,
            x147 * x153 * x185,
            x149 * x167 * x191,
            x167 * x169 * x189,
            x153 * x167 * x187,
            x167 * x172 * x185,
            x149 * x166 * x196,
            x166 * x169 * x191,
            x153 * x166 * x189,
            x166 * x172 * x187,
            x166 * x182 * x185,
            x184 * x197,
            x105 * x122 * x198,
            x122 * x126 * x201,
            x129 * x147 * x199,
            x147 * x155 * x201,
            x126 * x147 * x203,
            x157 * x167 * x199,
            x129 * x167 * x201,
            x155 * x167 * x203,
            x126 * x167 * x205,
            x166 * x176 * x199,
            x157 * x166 * x201,
            x129 * x166 * x203,
            x155 * x166 * x205,
            x126 * x166 * x210,
            x149 * x212 * x98,
            x149 * x215 * x96,
            x169 * x212 * x96,
            x120 * x149 * x220,
            x120 * x169 * x215,
            x120 * x153 * x212,
            x145 * x149 * x223,
            x145 * x169 * x220,
            x145 * x153 * x215,
            x145 * x172 * x212,
            x143 * x149 * x226,
            x143 * x169 * x223,
            x143 * x153 * x220,
            x143 * x172 * x215,
            x143 * x182 * x212,
            x183 * x198 * x98,
            x187 * x199 * x96,
            x185 * x201 * x96,
            x120 * x189 * x199,
            x120 * x187 * x201,
            x120 * x185 * x203,
            x145 * x191 * x199,
            x145 * x189 * x201,
            x145 * x187 * x203,
            x145 * x185 * x205,
            x143 * x196 * x199,
            x143 * x191 * x201,
            x143 * x189 * x203,
            x143 * x187 * x205,
            x143 * x185 * x210,
            x126 * x228 * x98,
            x155 * x228 * x96,
            x126 * x231 * x96,
            x120 * x129 * x228,
            x120 * x155 * x231,
            x120 * x126 * x236,
            x145 * x157 * x228,
            x129 * x145 * x231,
            x145 * x155 * x236,
            x126 * x145 * x239,
            x143 * x176 * x228,
            x143 * x157 * x231,
            x129 * x143 * x236,
            x143 * x155 * x239,
            x126 * x143 * x242,
            x149 * x243 * x81,
            x149 * x246 * x94,
            x169 * x243 * x94,
            x149 * x249 * x88,
            x169 * x246 * x88,
            x153 * x243 * x88,
            x117 * x149 * x253,
            x117 * x169 * x249,
            x117 * x153 * x246,
            x117 * x172 * x243,
            x141 * x149 * x254,
            x141 * x169 * x253,
            x141 * x153 * x249,
            x141 * x172 * x246,
            x141 * x182 * x243,
            x199 * x212 * x81,
            x199 * x215 * x94,
            x201 * x212 * x94,
            x199 * x220 * x88,
            x201 * x215 * x88,
            x203 * x212 * x88,
            x117 * x199 * x223,
            x117 * x201 * x220,
            x117 * x203 * x215,
            x117 * x205 * x212,
            x141 * x199 * x226,
            x141 * x201 * x223,
            x141 * x203 * x220,
            x141 * x205 * x215,
            x141 * x210 * x212,
            x185 * x228 * x81,
            x187 * x228 * x94,
            x185 * x231 * x94,
            x189 * x228 * x88,
            x187 * x231 * x88,
            x185 * x236 * x88,
            x117 * x191 * x228,
            x117 * x189 * x231,
            x117 * x187 * x236,
            x117 * x185 * x239,
            x141 * x196 * x228,
            x141 * x191 * x231,
            x141 * x189 * x236,
            x141 * x187 * x239,
            x141 * x185 * x242,
            x126 * x255 * x81,
            x155 * x255 * x94,
            x126 * x258 * x94,
            x129 * x255 * x88,
            x155 * x258 * x88,
            x126 * x261 * x88,
            x117 * x157 * x255,
            x117 * x129 * x258,
            x117 * x155 * x261,
            x117 * x126 * x265,
            x141 * x176 * x255,
            x141 * x157 * x258,
            x129 * x141 * x261,
            x141 * x155 * x265,
            x126 * x141 * x266,
            x149 * x267 * x79,
            x149 * x268 * x69,
            x169 * x267 * x69,
            x149 * x269 * x67,
            x169 * x268 * x67,
            x153 * x267 * x67,
            x149 * x270 * x65,
            x169 * x269 * x65,
            x153 * x268 * x65,
            x172 * x267 * x65,
            x149 * x271 * x57,
            x169 * x270 * x57,
            x153 * x269 * x57,
            x172 * x268 * x57,
            x182 * x267 * x57,
            x199 * x243 * x79,
            x199 * x246 * x69,
            x201 * x243 * x69,
            x199 * x249 * x67,
            x201 * x246 * x67,
            x203 * x243 * x67,
            x199 * x253 * x65,
            x201 * x249 * x65,
            x203 * x246 * x65,
            x205 * x243 * x65,
            x199 * x254 * x57,
            x201 * x253 * x57,
            x203 * x249 * x57,
            x205 * x246 * x57,
            x210 * x243 * x57,
            x212 * x228 * x79,
            x215 * x228 * x69,
            x212 * x231 * x69,
            x220 * x228 * x67,
            x215 * x231 * x67,
            x212 * x236 * x67,
            x223 * x228 * x65,
            x220 * x231 * x65,
            x215 * x236 * x65,
            x212 * x239 * x65,
            x226 * x228 * x57,
            x223 * x231 * x57,
            x220 * x236 * x57,
            x215 * x239 * x57,
            x212 * x242 * x57,
            x185 * x255 * x79,
            x187 * x255 * x69,
            x185 * x258 * x69,
            x189 * x255 * x67,
            x187 * x258 * x67,
            x185 * x261 * x67,
            x191 * x255 * x65,
            x189 * x258 * x65,
            x187 * x261 * x65,
            x185 * x265 * x65,
            x196 * x255 * x57,
            x191 * x258 * x57,
            x189 * x261 * x57,
            x187 * x265 * x57,
            x185 * x266 * x57,
            x126 * x272 * x79,
            x155 * x272 * x69,
            x126 * x273 * x69,
            x129 * x272 * x67,
            x155 * x273 * x67,
            x126 * x274 * x67,
            x157 * x272 * x65,
            x129 * x273 * x65,
            x155 * x274 * x65,
            x126 * x275 * x65,
            x176 * x272 * x57,
            x157 * x273 * x57,
            x129 * x274 * x57,
            x155 * x275 * x57,
            x126 * x276 * x57,
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
            x149 * x283 * x325,
            x149 * x291 * x327,
            x169 * x291 * x325,
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
            x199 * x279 * x283,
            x199 * x288 * x291,
            x201 * x279 * x291,
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
            x149 * x344 * x48,
            x114 * x149 * x349,
            x114 * x169 * x344,
            x133 * x149 * x355,
            x133 * x169 * x349,
            x133 * x153 * x344,
            x149 * x165 * x360,
            x165 * x169 * x355,
            x153 * x165 * x349,
            x165 * x172 * x344,
            x149 * x311 * x367,
            x169 * x311 * x360,
            x153 * x311 * x355,
            x172 * x311 * x349,
            x182 * x311 * x344,
            x199 * x325 * x48,
            x114 * x199 * x327,
            x114 * x201 * x325,
            x133 * x199 * x329,
            x133 * x201 * x327,
            x133 * x203 * x325,
            x165 * x199 * x331,
            x165 * x201 * x329,
            x165 * x203 * x327,
            x165 * x205 * x325,
            x199 * x311 * x340,
            x201 * x311 * x331,
            x203 * x311 * x329,
            x205 * x311 * x327,
            x210 * x311 * x325,
            x228 * x279 * x48,
            x114 * x228 * x288,
            x114 * x231 * x279,
            x133 * x228 * x297,
            x133 * x231 * x288,
            x133 * x236 * x279,
            x165 * x228 * x310,
            x165 * x231 * x297,
            x165 * x236 * x288,
            x165 * x239 * x279,
            x228 * x311 * x322,
            x231 * x310 * x311,
            x236 * x297 * x311,
            x239 * x288 * x311,
            x242 * x279 * x311,
            x149 * x371 * x46,
            x149 * x32 * x377,
            x169 * x32 * x371,
            x109 * x149 * x382,
            x109 * x169 * x377,
            x109 * x153 * x371,
            x149 * x163 * x388,
            x163 * x169 * x382,
            x153 * x163 * x377,
            x163 * x172 * x371,
            x23 * x393,
            x125 * x388 * x394,
            x153 * x158 * x382,
            x158 * x172 * x377,
            x158 * x182 * x371,
            x199 * x344 * x46,
            x199 * x32 * x349,
            x201 * x32 * x344,
            x109 * x199 * x355,
            x109 * x201 * x349,
            x109 * x203 * x344,
            x163 * x199 * x360,
            x163 * x201 * x355,
            x163 * x203 * x349,
            x163 * x205 * x344,
            x197 * x367 * x394,
            x158 * x201 * x360,
            x158 * x203 * x355,
            x158 * x205 * x349,
            x158 * x210 * x344,
            x228 * x325 * x46,
            x228 * x32 * x327,
            x231 * x32 * x325,
            x109 * x228 * x329,
            x109 * x231 * x327,
            x109 * x236 * x325,
            x163 * x228 * x331,
            x163 * x231 * x329,
            x163 * x236 * x327,
            x163 * x239 * x325,
            x158 * x228 * x340,
            x158 * x231 * x331,
            x158 * x236 * x329,
            x158 * x239 * x327,
            x158 * x242 * x325,
            x255 * x279 * x46,
            x255 * x288 * x32,
            x258 * x279 * x32,
            x109 * x255 * x297,
            x109 * x258 * x288,
            x109 * x261 * x279,
            x163 * x255 * x310,
            x163 * x258 * x297,
            x163 * x261 * x288,
            x163 * x265 * x279,
            x158 * x255 * x322,
            x158 * x258 * x310,
            x158 * x261 * x297,
            x158 * x265 * x288,
            x158 * x266 * x279,
            x149 * x396 * x44,
            x149 * x30 * x397,
            x169 * x30 * x396,
            x149 * x26 * x398,
            x169 * x26 * x397,
            x153 * x26 * x396,
            x3 * x399,
            x125 * x398 * x400,
            x153 * x397 * x8,
            x172 * x396 * x8,
            x392
            * (
                x0
                * (
                    x216 * (x389 + x390)
                    + x28 * (x226 + x364 + 4.0 * x383 + 4.0 * x384)
                    + 3.0 * x365
                    + 3.0 * x366
                    + 4.0 * x386
                    + 4.0 * x387
                )
                + x183 * x391
            ),
            x125 * x399,
            x153 * x398 * x7,
            x172 * x397 * x7,
            x182 * x396 * x7,
            x199 * x371 * x44,
            x199 * x30 * x377,
            x201 * x30 * x371,
            x199 * x26 * x382,
            x201 * x26 * x377,
            x203 * x26 * x371,
            x197 * x388 * x400,
            x201 * x382 * x8,
            x203 * x377 * x8,
            x205 * x371 * x8,
            x197 * x393,
            x201 * x388 * x7,
            x203 * x382 * x7,
            x205 * x377 * x7,
            x210 * x371 * x7,
            x228 * x344 * x44,
            x228 * x30 * x349,
            x231 * x30 * x344,
            x228 * x26 * x355,
            x231 * x26 * x349,
            x236 * x26 * x344,
            x228 * x360 * x8,
            x231 * x355 * x8,
            x236 * x349 * x8,
            x239 * x344 * x8,
            x228 * x367 * x7,
            x231 * x360 * x7,
            x236 * x355 * x7,
            x239 * x349 * x7,
            x242 * x344 * x7,
            x255 * x325 * x44,
            x255 * x30 * x327,
            x258 * x30 * x325,
            x255 * x26 * x329,
            x258 * x26 * x327,
            x26 * x261 * x325,
            x255 * x331 * x8,
            x258 * x329 * x8,
            x261 * x327 * x8,
            x265 * x325 * x8,
            x255 * x340 * x7,
            x258 * x331 * x7,
            x261 * x329 * x7,
            x265 * x327 * x7,
            x266 * x325 * x7,
            x272 * x279 * x44,
            x272 * x288 * x30,
            x273 * x279 * x30,
            x26 * x272 * x297,
            x26 * x273 * x288,
            x26 * x274 * x279,
            x272 * x310 * x8,
            x273 * x297 * x8,
            x274 * x288 * x8,
            x275 * x279 * x8,
            x272 * x322 * x7,
            x273 * x310 * x7,
            x274 * x297 * x7,
            x275 * x288 * x7,
            x276 * x279 * x7,
            x126 * x284 * x403,
            x155 * x292 * x403,
            x126 * x292 * x407,
            x129 * x301 * x403,
            x155 * x301 * x407,
            x126 * x301 * x412,
            x157 * x314 * x403,
            x129 * x314 * x407,
            x155 * x314 * x412,
            x126 * x314 * x421,
            x176 * x323 * x403,
            x157 * x323 * x407,
            x129 * x323 * x412,
            x155 * x323 * x421,
            x126 * x323 * x429,
            x185 * x283 * x403,
            x187 * x291 * x403,
            x185 * x291 * x407,
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
            x126 * x283 * x431,
            x155 * x291 * x431,
            x126 * x291 * x433,
            x129 * x300 * x431,
            x155 * x300 * x433,
            x126 * x300 * x435,
            x157 * x313 * x431,
            x129 * x313 * x433,
            x155 * x313 * x435,
            x126 * x313 * x437,
            x176 * x312 * x431,
            x157 * x312 * x433,
            x129 * x312 * x435,
            x155 * x312 * x437,
            x126 * x312 * x446,
            x212 * x403 * x48,
            x114 * x215 * x403,
            x114 * x212 * x407,
            x133 * x220 * x403,
            x133 * x215 * x407,
            x133 * x212 * x412,
            x165 * x223 * x403,
            x165 * x220 * x407,
            x165 * x215 * x412,
            x165 * x212 * x421,
            x226 * x311 * x403,
            x223 * x311 * x407,
            x220 * x311 * x412,
            x215 * x311 * x421,
            x212 * x311 * x429,
            x185 * x431 * x48,
            x114 * x187 * x431,
            x114 * x185 * x433,
            x133 * x189 * x431,
            x133 * x187 * x433,
            x133 * x185 * x435,
            x165 * x191 * x431,
            x165 * x189 * x433,
            x165 * x187 * x435,
            x165 * x185 * x437,
            x196 * x311 * x431,
            x191 * x311 * x433,
            x189 * x311 * x435,
            x187 * x311 * x437,
            x185 * x311 * x446,
            x126 * x450 * x48,
            x114 * x155 * x450,
            x114 * x126 * x455,
            x129 * x133 * x450,
            x133 * x155 * x455,
            x126 * x133 * x461,
            x157 * x165 * x450,
            x129 * x165 * x455,
            x155 * x165 * x461,
            x126 * x165 * x466,
            x176 * x311 * x450,
            x157 * x311 * x455,
            x129 * x311 * x461,
            x155 * x311 * x466,
            x126 * x311 * x473,
            x243 * x403 * x46,
            x246 * x32 * x403,
            x243 * x32 * x407,
            x109 * x249 * x403,
            x109 * x246 * x407,
            x109 * x243 * x412,
            x163 * x253 * x403,
            x163 * x249 * x407,
            x163 * x246 * x412,
            x163 * x243 * x421,
            x158 * x254 * x403,
            x158 * x253 * x407,
            x158 * x249 * x412,
            x158 * x246 * x421,
            x158 * x243 * x429,
            x212 * x431 * x46,
            x215 * x32 * x431,
            x212 * x32 * x433,
            x109 * x220 * x431,
            x109 * x215 * x433,
            x109 * x212 * x435,
            x163 * x223 * x431,
            x163 * x220 * x433,
            x163 * x215 * x435,
            x163 * x212 * x437,
            x158 * x226 * x431,
            x158 * x223 * x433,
            x158 * x220 * x435,
            x158 * x215 * x437,
            x158 * x212 * x446,
            x185 * x450 * x46,
            x187 * x32 * x450,
            x185 * x32 * x455,
            x109 * x189 * x450,
            x109 * x187 * x455,
            x109 * x185 * x461,
            x163 * x191 * x450,
            x163 * x189 * x455,
            x163 * x187 * x461,
            x163 * x185 * x466,
            x158 * x196 * x450,
            x158 * x191 * x455,
            x158 * x189 * x461,
            x158 * x187 * x466,
            x183 * x473 * x475,
            x126 * x46 * x479,
            x155 * x32 * x479,
            x126 * x32 * x485,
            x109 * x129 * x479,
            x109 * x155 * x485,
            x109 * x126 * x490,
            x157 * x163 * x479,
            x129 * x163 * x485,
            x155 * x163 * x490,
            x126 * x163 * x496,
            x158 * x176 * x479,
            x157 * x158 * x485,
            x129 * x158 * x490,
            x105 * x475 * x496,
            x23 * x500,
            x267 * x403 * x44,
            x268 * x30 * x403,
            x267 * x30 * x407,
            x26 * x269 * x403,
            x26 * x268 * x407,
            x26 * x267 * x412,
            x270 * x403 * x8,
            x269 * x407 * x8,
            x268 * x412 * x8,
            x267 * x421 * x8,
            x271 * x403 * x7,
            x270 * x407 * x7,
            x269 * x412 * x7,
            x268 * x421 * x7,
            x267 * x429 * x7,
            x243 * x431 * x44,
            x246 * x30 * x431,
            x243 * x30 * x433,
            x249 * x26 * x431,
            x246 * x26 * x433,
            x243 * x26 * x435,
            x253 * x431 * x8,
            x249 * x433 * x8,
            x246 * x435 * x8,
            x243 * x437 * x8,
            x254 * x431 * x7,
            x253 * x433 * x7,
            x249 * x435 * x7,
            x246 * x437 * x7,
            x243 * x446 * x7,
            x212 * x44 * x450,
            x215 * x30 * x450,
            x212 * x30 * x455,
            x220 * x26 * x450,
            x215 * x26 * x455,
            x212 * x26 * x461,
            x223 * x450 * x8,
            x220 * x455 * x8,
            x215 * x461 * x8,
            x212 * x466 * x8,
            x226 * x450 * x7,
            x223 * x455 * x7,
            x220 * x461 * x7,
            x215 * x466 * x7,
            x212 * x473 * x7,
            x185 * x44 * x479,
            x187 * x30 * x479,
            x185 * x30 * x485,
            x189 * x26 * x479,
            x187 * x26 * x485,
            x185 * x26 * x490,
            x191 * x479 * x8,
            x189 * x485 * x8,
            x187 * x490 * x8,
            x183 * x496 * x501,
            x196 * x479 * x7,
            x191 * x485 * x7,
            x189 * x490 * x7,
            x187 * x496 * x7,
            x183 * x500,
            x126 * x44 * x503,
            x155 * x30 * x503,
            x126 * x30 * x504,
            x129 * x26 * x503,
            x155 * x26 * x504,
            x126 * x26 * x505,
            x157 * x503 * x8,
            x129 * x504 * x8,
            x105 * x501 * x505,
            x3 * x506,
            x176 * x503 * x7,
            x157 * x504 * x7,
            x129 * x505 * x7,
            x105 * x506,
            x474
            * (
                x0
                * (
                    x232 * (x497 + x498)
                    + x28 * (x242 + x470 + 4.0 * x491 + 4.0 * x492)
                    + 3.0 * x471
                    + 3.0 * x472
                    + 4.0 * x494
                    + 4.0 * x495
                )
                + x197 * x499
            ),
        ]
    )
