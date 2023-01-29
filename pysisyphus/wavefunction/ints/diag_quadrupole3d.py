"""

        Diagonal of the quadrupole moment matrix with operators x², y², z².

        for rr in (xx, yy, zz):
            for bf_a in basis_functions_a:
                for bf_b in basis_functions_b:
                        quadrupole_integrals(bf_a, bf_b, rr)

"""


import numpy


_L_MAX = 4


def diag_quadrupole3d_00(ax, da, A, bx, db, B, C):
    """Cartesian 3D (ss) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 1, 1), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = ax * bx * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x4 = 1.77245385090552 * numpy.sqrt(x0)
    x5 = x3 * x4
    x6 = 0.5 / (ax + bx)
    x7 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x8 = 0.564189583547756
    x9 = numpy.sqrt(ax**1.5)
    x10 = numpy.sqrt(bx**1.5)
    x11 = 2.82842712474619 * da * db * x0 * x10 * x7 * x8 * x9
    x12 = x2 * x4
    x13 = x4 * x7

    # 3 item(s)
    result[0, 0, 0] = numpy.sum(
        x11 * x2 * x5 * (x6 + (-x0 * (ax * A[0] + bx * B[0]) + C[0]) ** 2)
    )
    result[1, 0, 0] = numpy.sum(
        x11 * x12 * x3 * (x6 + (-x0 * (ax * A[1] + bx * B[1]) + C[1]) ** 2)
    )
    result[2, 0, 0] = numpy.sum(
        2.82842712474619
        * da
        * db
        * x0
        * x10
        * x13
        * x2
        * x3
        * x8
        * x9
        * (x6 + (-x0 * (ax * A[2] + bx * B[2]) + C[2]) ** 2)
    )
    return result


def diag_quadrupole3d_01(ax, da, A, bx, db, B, C):
    """Cartesian 3D (sp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 1, 3), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = -x0 * (ax * A[0] + bx * B[0])
    x2 = -x1 - B[0]
    x3 = -x1 - C[0]
    x4 = ax * bx * x0
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x0)
    x7 = x5 * x6
    x8 = 0.5 / (ax + bx)
    x9 = x7 * x8
    x10 = x3**2 * x7 + x9
    x11 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x12 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x13 = 0.564189583547756
    x14 = numpy.sqrt(ax**1.5)
    x15 = numpy.sqrt(bx**2.5)
    x16 = 5.65685424949238 * da * db * x0 * x12 * x13 * x14 * x15
    x17 = x11 * x16
    x18 = -x0 * (ax * A[1] + bx * B[1])
    x19 = -x18 - B[1]
    x20 = x10 * x17
    x21 = -x0 * (ax * A[2] + bx * B[2])
    x22 = -x21 - B[2]
    x23 = -x18 - C[1]
    x24 = x11 * x6
    x25 = x24 * x8
    x26 = x23**2 * x24 + x25
    x27 = x16 * x5
    x28 = x26 * x27
    x29 = -x21 - C[2]
    x30 = x12 * x6
    x31 = x30 * x8
    x32 = x29**2 * x30 + x31
    x33 = 5.65685424949238 * da * db * x0 * x11 * x13 * x14 * x15 * x5
    x34 = x32 * x33

    # 9 item(s)
    result[0, 0, 0] = numpy.sum(x17 * (x10 * x2 + 2.0 * x3 * x9))
    result[0, 0, 1] = numpy.sum(x19 * x20)
    result[0, 0, 2] = numpy.sum(x20 * x22)
    result[1, 0, 0] = numpy.sum(x2 * x28)
    result[1, 0, 1] = numpy.sum(x27 * (x19 * x26 + 2.0 * x23 * x25))
    result[1, 0, 2] = numpy.sum(x22 * x28)
    result[2, 0, 0] = numpy.sum(x2 * x34)
    result[2, 0, 1] = numpy.sum(x19 * x34)
    result[2, 0, 2] = numpy.sum(x33 * (x22 * x32 + 2.0 * x29 * x31))
    return result


def diag_quadrupole3d_02(ax, da, A, bx, db, B, C):
    """Cartesian 3D (sd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 1, 6), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - C[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = -x2 - B[0]
    x12 = 2.0 * x3
    x13 = x10 + x9
    x14 = x10 * x12 + x11 * x13
    x15 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x16 = numpy.sqrt(ax**1.5)
    x17 = numpy.sqrt(bx**3.5)
    x18 = 6.53197264742181 * da * db * x16 * x17
    x19 = x15 * x18
    x20 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x21 = 0.564189583547756 * x1
    x22 = x20 * x21
    x23 = -x1 * (ax * A[1] + bx * B[1])
    x24 = -x23 - B[1]
    x25 = 11.3137084989848 * da * db * x16 * x17
    x26 = x15 * x25
    x27 = x14 * x22 * x26
    x28 = -x1 * (ax * A[2] + bx * B[2])
    x29 = -x28 - B[2]
    x30 = x20 * x7
    x31 = x0 * x30
    x32 = x24**2 * x30 + x31
    x33 = 0.318309886183791 * x6
    x34 = x13 * x33
    x35 = x26 * x29
    x36 = x15 * x7
    x37 = x0 * x36
    x38 = x18 * (x29**2 * x36 + x37)
    x39 = x10 + x11**2 * x8
    x40 = -x23 - C[1]
    x41 = x30 * x40**2
    x42 = x31 + x41
    x43 = x33 * x42
    x44 = 2.0 * x40
    x45 = x24 * x42 + x31 * x44
    x46 = x21 * x5
    x47 = x26 * x45 * x46
    x48 = -x28 - C[2]
    x49 = x36 * x48**2
    x50 = x37 + x49
    x51 = x18 * x33 * x50
    x52 = x22 * x5
    x53 = 2.0 * x48
    x54 = x29 * x50 + x37 * x53
    x55 = x25 * x52 * x54

    # 18 item(s)
    result[0, 0, 0] = numpy.sum(
        x19 * x22 * (x0 * (3.0 * x10 + x11 * x12 * x8 + x9) + x11 * x14)
    )
    result[0, 0, 1] = numpy.sum(x24 * x27)
    result[0, 0, 2] = numpy.sum(x27 * x29)
    result[0, 0, 3] = numpy.sum(x19 * x32 * x34)
    result[0, 0, 4] = numpy.sum(x13 * x22 * x24 * x35)
    result[0, 0, 5] = numpy.sum(x20 * x34 * x38)
    result[1, 0, 0] = numpy.sum(x19 * x39 * x43)
    result[1, 0, 1] = numpy.sum(x11 * x47)
    result[1, 0, 2] = numpy.sum(x11 * x35 * x42 * x46)
    result[1, 0, 3] = numpy.sum(
        x19 * x46 * (x0 * (x24 * x30 * x44 + 3.0 * x31 + x41) + x24 * x45)
    )
    result[1, 0, 4] = numpy.sum(x29 * x47)
    result[1, 0, 5] = numpy.sum(x38 * x43 * x5)
    result[2, 0, 0] = numpy.sum(x20 * x39 * x51)
    result[2, 0, 1] = numpy.sum(x11 * x24 * x25 * x50 * x52)
    result[2, 0, 2] = numpy.sum(x11 * x55)
    result[2, 0, 3] = numpy.sum(x32 * x5 * x51)
    result[2, 0, 4] = numpy.sum(x24 * x55)
    result[2, 0, 5] = numpy.sum(
        x18 * x52 * (x0 * (x29 * x36 * x53 + 3.0 * x37 + x49) + x29 * x54)
    )
    return result


def diag_quadrupole3d_03(ax, da, A, bx, db, B, C):
    """Cartesian 3D (sf) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 1, 10), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = numpy.sqrt(x1)
    x3 = 1.77245385090552 * x2
    x4 = -x1 * (ax * A[0] + bx * B[0])
    x5 = -x4 - B[0]
    x6 = ax * bx * x1
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x4 - C[0]
    x11 = x3 * x7
    x12 = x10 * x11
    x13 = 2.0 * x0
    x14 = x10**2 * x11
    x15 = x0 * x11
    x16 = x14 + x15
    x17 = x16 * x5
    x18 = 2.0 * x5
    x19 = x12 * x13 + x17
    x20 = x0 * (x12 * x18 + x14 + 3.0 * x15) + x19 * x5
    x21 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x22 = da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**4.5)
    x23 = x21 * x22
    x24 = 5.84237394672177 * x23
    x25 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x26 = 0.564189583547756 * x1
    x27 = x25 * x26
    x28 = -x1 * (ax * A[1] + bx * B[1])
    x29 = -x28 - B[1]
    x30 = x27 * x29
    x31 = 13.0639452948436 * x23
    x32 = x20 * x31
    x33 = -x1 * (ax * A[2] + bx * B[2])
    x34 = -x33 - B[2]
    x35 = 0.318309886183791 * x2
    x36 = x19 * x35
    x37 = x25 * x3
    x38 = x0 * x37
    x39 = x29**2 * x37 + x38
    x40 = x31 * x39
    x41 = 22.6274169979695 * x30
    x42 = x23 * x34
    x43 = x22 * x25
    x44 = x21 * x3
    x45 = x0 * x44
    x46 = x34**2 * x44 + x45
    x47 = 13.0639452948436 * x46
    x48 = x43 * x47
    x49 = x29 * x37
    x50 = x13 * x49 + x29 * x39
    x51 = x16 * x35
    x52 = x34 * x44
    x53 = 5.84237394672177 * x13 * x52 + 5.84237394672177 * x34 * x46
    x54 = x11 * x5**2 + x15
    x55 = x13 * x9 + x5 * x54
    x56 = -x28 - C[1]
    x57 = x37 * x56**2
    x58 = x38 + x57
    x59 = x35 * x58
    x60 = x29 * x58
    x61 = x37 * x56
    x62 = x13 * x61 + x60
    x63 = x35 * x62
    x64 = x31 * x54
    x65 = 2.0 * x29
    x66 = x0 * (3.0 * x38 + x57 + x61 * x65) + x29 * x62
    x67 = x31 * x66
    x68 = x26 * x8
    x69 = x22 * x8
    x70 = x26 * x7
    x71 = x22 * x7
    x72 = x35 * x43
    x73 = -x33 - C[2]
    x74 = x44 * x73**2
    x75 = x45 + x74
    x76 = 5.84237394672177 * x75
    x77 = x54 * x72
    x78 = 13.0639452948436 * x75
    x79 = x34 * x75
    x80 = x44 * x73
    x81 = x13 * x80 + x79
    x82 = 13.0639452948436 * x81
    x83 = 2.0 * x34
    x84 = x0 * (3.0 * x45 + x74 + x80 * x83) + x34 * x81
    x85 = 13.0639452948436 * x84
    x86 = x35 * x71
    x87 = x27 * x71

    # 30 item(s)
    result[0, 0, 0] = numpy.sum(
        x24
        * x27
        * (
            x0 * (4.0 * x0 * x12 + x13 * (x12 + x9) + 2.0 * x17 + x18 * (x10 * x9 + x15))
            + x20 * x5
        )
    )
    result[0, 0, 1] = numpy.sum(x30 * x32)
    result[0, 0, 2] = numpy.sum(x27 * x32 * x34)
    result[0, 0, 3] = numpy.sum(x36 * x40)
    result[0, 0, 4] = numpy.sum(x19 * x41 * x42)
    result[0, 0, 5] = numpy.sum(x36 * x48)
    result[0, 0, 6] = numpy.sum(x24 * x50 * x51)
    result[0, 0, 7] = numpy.sum(x34 * x40 * x51)
    result[0, 0, 8] = numpy.sum(x29 * x48 * x51)
    result[0, 0, 9] = numpy.sum(x43 * x51 * x53)
    result[1, 0, 0] = numpy.sum(x24 * x55 * x59)
    result[1, 0, 1] = numpy.sum(x63 * x64)
    result[1, 0, 2] = numpy.sum(x34 * x59 * x64)
    result[1, 0, 3] = numpy.sum(x67 * x68)
    result[1, 0, 4] = numpy.sum(22.6274169979695 * x42 * x62 * x68)
    result[1, 0, 5] = numpy.sum(x47 * x59 * x69)
    result[1, 0, 6] = numpy.sum(
        x24
        * x70
        * (
            x0
            * (x13 * (x49 + x61) + 4.0 * x38 * x56 + 2.0 * x60 + x65 * (x38 + x49 * x56))
            + x29 * x66
        )
    )
    result[1, 0, 7] = numpy.sum(x34 * x67 * x70)
    result[1, 0, 8] = numpy.sum(x47 * x63 * x71)
    result[1, 0, 9] = numpy.sum(x53 * x59 * x71)
    result[2, 0, 0] = numpy.sum(x55 * x72 * x76)
    result[2, 0, 1] = numpy.sum(x29 * x77 * x78)
    result[2, 0, 2] = numpy.sum(x77 * x82)
    result[2, 0, 3] = numpy.sum(x35 * x39 * x69 * x78)
    result[2, 0, 4] = numpy.sum(x41 * x69 * x81)
    result[2, 0, 5] = numpy.sum(x27 * x69 * x85)
    result[2, 0, 6] = numpy.sum(x50 * x76 * x86)
    result[2, 0, 7] = numpy.sum(x39 * x82 * x86)
    result[2, 0, 8] = numpy.sum(x29 * x85 * x87)
    result[2, 0, 9] = numpy.sum(
        5.84237394672177
        * x87
        * (
            x0
            * (x13 * (x52 + x80) + 4.0 * x45 * x73 + 2.0 * x79 + x83 * (x45 + x52 * x73))
            + x34 * x84
        )
    )
    return result


def diag_quadrupole3d_04(ax, da, A, bx, db, B, C):
    """Cartesian 3D (sg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 1, 15), dtype=float)

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
    x12 = 2.0 * x3
    x13 = -x2 - C[0]
    x14 = x13 * x8
    x15 = x11 + x12 * x14
    x16 = 2.0 * x0
    x17 = x3 * x8
    x18 = x0 * (x14 + x17)
    x19 = x3 * (x10 + x13 * x17)
    x20 = x13**2 * x8
    x21 = x0 * (x15 + x20)
    x22 = x10 + x20
    x23 = x22 * x3
    x24 = x14 * x16 + x23
    x25 = x24 * x3
    x26 = x21 + x25
    x27 = 2.0 * x0 * (2.0 * x0 * x14 + x18 + x19 + x23) + x26 * x3
    x28 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x29 = 4.41641957979107 * x28
    x30 = 0.564189583547756
    x31 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x32 = da * db * numpy.sqrt(ax**1.5) * numpy.sqrt(bx**5.5)
    x33 = x31 * x32
    x34 = x1 * x30 * x33
    x35 = -x1 * (ax * A[1] + bx * B[1])
    x36 = -x35 - B[1]
    x37 = 11.6847478934435 * x36
    x38 = x28 * x37
    x39 = x27 * x34
    x40 = -x1 * (ax * A[2] + bx * B[2])
    x41 = -x40 - B[2]
    x42 = 11.6847478934435 * x41
    x43 = x28 * x7
    x44 = x36**2 * x43
    x45 = x0 * x43
    x46 = x44 + x45
    x47 = 15.084944665313 * x46
    x48 = 0.318309886183791 * x6
    x49 = x26 * x48
    x50 = 26.1278905896872 * x41
    x51 = x28 * x32
    x52 = x31 * x7
    x53 = x41**2 * x52
    x54 = x0 * x52
    x55 = x53 + x54
    x56 = 15.084944665313 * x55
    x57 = x36 * x43
    x58 = x16 * x57 + x36 * x46
    x59 = x33 * x48
    x60 = 11.6847478934435 * x24
    x61 = 26.1278905896872 * x24
    x62 = x48 * x51
    x63 = x41 * x52
    x64 = x16 * x63 + x41 * x55
    x65 = x48 * x64
    x66 = 3.0 * x45
    x67 = 4.41641957979107 * x0 * (3.0 * x44 + x66) + 4.41641957979107 * x36 * x58
    x68 = x22 * x32
    x69 = x48 * x68
    x70 = x31 * x69
    x71 = 0.179587122125167
    x72 = x55 * x71
    x73 = 3.0 * x54
    x74 = x0 * (3.0 * x53 + x73) + x41 * x64
    x75 = x10 + x9
    x76 = x16 * x17 + x3 * x75
    x77 = x0 * (x11 + 3.0 * x9) + x3 * x76
    x78 = -x35 - C[1]
    x79 = x43 * x78**2
    x80 = x45 + x79
    x81 = 4.41641957979107 * x80
    x82 = x36 * x80
    x83 = x43 * x78
    x84 = x16 * x83 + x82
    x85 = 11.6847478934435 * x84
    x86 = x59 * x76
    x87 = x32 * x75
    x88 = 15.084944665313 * x87
    x89 = 2.0 * x36
    x90 = x66 + x83 * x89
    x91 = x0 * (x79 + x90)
    x92 = x36 * x84
    x93 = x91 + x92
    x94 = x48 * x93
    x95 = 26.1278905896872 * x84
    x96 = 11.6847478934435 * x3
    x97 = x0 * (x57 + x83)
    x98 = x36 * (x45 + x57 * x78)
    x99 = 2.0 * x0 * (2.0 * x45 * x78 + x82 + x97 + x98) + x36 * x93
    x100 = x34 * x5
    x101 = x100 * x99
    x102 = x32 * x5
    x103 = x102 * x48
    x104 = x103 * x3
    x105 = x102 * x65
    x106 = -x40 - C[2]
    x107 = x106**2 * x52
    x108 = x107 + x54
    x109 = x108 * x48
    x110 = x37 * x51
    x111 = x108 * x41
    x112 = x106 * x52
    x113 = x111 + x112 * x16
    x114 = 11.6847478934435 * x113
    x115 = x28 * x48
    x116 = 26.1278905896872 * x113
    x117 = 2.0 * x41
    x118 = x112 * x117 + x73
    x119 = x0 * (x107 + x118)
    x120 = x113 * x41
    x121 = x119 + x120
    x122 = x103 * x108
    x123 = x1 * x30
    x124 = x0 * (x112 + x63)
    x125 = x41 * (x106 * x63 + x54)
    x126 = 2.0 * x0 * (2.0 * x106 * x54 + x111 + x124 + x125) + x121 * x41
    x127 = x123 * x126 * x5

    # 45 item(s)
    result[0, 0, 0] = numpy.sum(
        x29
        * x34
        * (x0 * (x12 * (x18 + x19) + x16 * (x15 + x9) + 3.0 * x21 + 3.0 * x25) + x27 * x3)
    )
    result[0, 0, 1] = numpy.sum(x38 * x39)
    result[0, 0, 2] = numpy.sum(x28 * x39 * x42)
    result[0, 0, 3] = numpy.sum(x33 * x47 * x49)
    result[0, 0, 4] = numpy.sum(x26 * x28 * x34 * x36 * x50)
    result[0, 0, 5] = numpy.sum(x49 * x51 * x56)
    result[0, 0, 6] = numpy.sum(x58 * x59 * x60)
    result[0, 0, 7] = numpy.sum(x41 * x46 * x59 * x61)
    result[0, 0, 8] = numpy.sum(x36 * x55 * x61 * x62)
    result[0, 0, 9] = numpy.sum(x51 * x60 * x65)
    result[0, 0, 10] = numpy.sum(x67 * x70)
    result[0, 0, 11] = numpy.sum(x42 * x58 * x70)
    result[0, 0, 12] = numpy.sum(x47 * x68 * x72)
    result[0, 0, 13] = numpy.sum(x38 * x65 * x68)
    result[0, 0, 14] = numpy.sum(x29 * x69 * x74)
    result[1, 0, 0] = numpy.sum(x59 * x77 * x81)
    result[1, 0, 1] = numpy.sum(x85 * x86)
    result[1, 0, 2] = numpy.sum(x42 * x80 * x86)
    result[1, 0, 3] = numpy.sum(x31 * x88 * x94)
    result[1, 0, 4] = numpy.sum(x31 * x41 * x48 * x87 * x95)
    result[1, 0, 5] = numpy.sum(x72 * x80 * x88)
    result[1, 0, 6] = numpy.sum(x101 * x96)
    result[1, 0, 7] = numpy.sum(x100 * x3 * x50 * x93)
    result[1, 0, 8] = numpy.sum(x104 * x55 * x95)
    result[1, 0, 9] = numpy.sum(x105 * x80 * x96)
    result[1, 0, 10] = numpy.sum(
        4.41641957979107
        * x100
        * (
            x0 * (x16 * (x44 + x90) + x89 * (x97 + x98) + 3.0 * x91 + 3.0 * x92)
            + x36 * x99
        )
    )
    result[1, 0, 11] = numpy.sum(x101 * x42)
    result[1, 0, 12] = numpy.sum(x102 * x56 * x94)
    result[1, 0, 13] = numpy.sum(x105 * x85)
    result[1, 0, 14] = numpy.sum(x103 * x74 * x81)
    result[2, 0, 0] = numpy.sum(x109 * x29 * x32 * x77)
    result[2, 0, 1] = numpy.sum(x109 * x110 * x76)
    result[2, 0, 2] = numpy.sum(x114 * x62 * x76)
    result[2, 0, 3] = numpy.sum(x108 * x47 * x71 * x87)
    result[2, 0, 4] = numpy.sum(x115 * x116 * x36 * x87)
    result[2, 0, 5] = numpy.sum(x115 * x121 * x88)
    result[2, 0, 6] = numpy.sum(x122 * x58 * x96)
    result[2, 0, 7] = numpy.sum(x104 * x116 * x46)
    result[2, 0, 8] = numpy.sum(26.1278905896872 * x121 * x123 * x3 * x36 * x5 * x51)
    result[2, 0, 9] = numpy.sum(x127 * x51 * x96)
    result[2, 0, 10] = numpy.sum(x122 * x67)
    result[2, 0, 11] = numpy.sum(x103 * x114 * x58)
    result[2, 0, 12] = numpy.sum(x103 * x121 * x47)
    result[2, 0, 13] = numpy.sum(x110 * x127)
    result[2, 0, 14] = numpy.sum(
        x102
        * x123
        * x29
        * (
            x0 * (x117 * (x124 + x125) + 3.0 * x119 + 3.0 * x120 + x16 * (x118 + x53))
            + x126 * x41
        )
    )
    return result


def diag_quadrupole3d_10(ax, da, A, bx, db, B, C):
    """Cartesian 3D (ps) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 3, 1), dtype=float)

    x0 = (ax + bx) ** (-1.0)
    x1 = -x0 * (ax * A[0] + bx * B[0])
    x2 = -x1 - A[0]
    x3 = -x1 - C[0]
    x4 = ax * bx * x0
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x0)
    x7 = x5 * x6
    x8 = 0.5 / (ax + bx)
    x9 = x7 * x8
    x10 = x3**2 * x7 + x9
    x11 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x12 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x13 = 0.564189583547756
    x14 = numpy.sqrt(ax**2.5)
    x15 = numpy.sqrt(bx**1.5)
    x16 = 5.65685424949238 * da * db * x0 * x12 * x13 * x14 * x15
    x17 = x11 * x16
    x18 = -x0 * (ax * A[1] + bx * B[1])
    x19 = -x18 - A[1]
    x20 = x10 * x17
    x21 = -x0 * (ax * A[2] + bx * B[2])
    x22 = -x21 - A[2]
    x23 = -x18 - C[1]
    x24 = x11 * x6
    x25 = x24 * x8
    x26 = x23**2 * x24 + x25
    x27 = x16 * x5
    x28 = x26 * x27
    x29 = -x21 - C[2]
    x30 = x12 * x6
    x31 = x30 * x8
    x32 = x29**2 * x30 + x31
    x33 = 5.65685424949238 * da * db * x0 * x11 * x13 * x14 * x15 * x5
    x34 = x32 * x33

    # 9 item(s)
    result[0, 0, 0] = numpy.sum(x17 * (x10 * x2 + 2.0 * x3 * x9))
    result[0, 1, 0] = numpy.sum(x19 * x20)
    result[0, 2, 0] = numpy.sum(x20 * x22)
    result[1, 0, 0] = numpy.sum(x2 * x28)
    result[1, 1, 0] = numpy.sum(x27 * (x19 * x26 + 2.0 * x23 * x25))
    result[1, 2, 0] = numpy.sum(x22 * x28)
    result[2, 0, 0] = numpy.sum(x2 * x34)
    result[2, 1, 0] = numpy.sum(x19 * x34)
    result[2, 2, 0] = numpy.sum(x33 * (x22 * x32 + 2.0 * x29 * x31))
    return result


def diag_quadrupole3d_11(ax, da, A, bx, db, B, C):
    """Cartesian 3D (pp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 3, 3), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - C[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = -x2 - B[0]
    x12 = x11 * x8
    x13 = 2.0 * x3
    x14 = -x2 - A[0]
    x15 = x10 + x9
    x16 = x10 * x13
    x17 = x11 * x15 + x16
    x18 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x19 = 0.564189583547756
    x20 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x21 = 11.3137084989848 * da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**2.5)
    x22 = x20 * x21
    x23 = x1 * x19 * x22
    x24 = x18 * x23
    x25 = -x1 * (ax * A[1] + bx * B[1])
    x26 = -x25 - B[1]
    x27 = x24 * (x14 * x15 + x16)
    x28 = -x1 * (ax * A[2] + bx * B[2])
    x29 = -x28 - B[2]
    x30 = -x25 - A[1]
    x31 = x17 * x24
    x32 = x0 * x7
    x33 = x18 * x32
    x34 = x18 * x7
    x35 = x26 * x34
    x36 = x30 * x35 + x33
    x37 = 0.318309886183791 * x6
    x38 = x15 * x37
    x39 = x15 * x24
    x40 = -x28 - A[2]
    x41 = x20 * x32
    x42 = x20 * x7
    x43 = x29 * x42
    x44 = x21 * (x40 * x43 + x41)
    x45 = x10 + x12 * x14
    x46 = -x25 - C[1]
    x47 = x34 * x46**2
    x48 = x33 + x47
    x49 = x37 * x48
    x50 = 2.0 * x46
    x51 = x33 * x50
    x52 = x26 * x48 + x51
    x53 = x23 * x5
    x54 = x52 * x53
    x55 = x48 * x53
    x56 = x53 * (x30 * x48 + x51)
    x57 = x18 * x21
    x58 = -x28 - C[2]
    x59 = x42 * x58**2
    x60 = x41 + x59
    x61 = x37 * x60
    x62 = x1 * x19 * x5 * x57
    x63 = x60 * x62
    x64 = 2.0 * x58
    x65 = x41 * x64
    x66 = x29 * x60 + x65
    x67 = x62 * x66
    x68 = x62 * (x40 * x60 + x65)

    # 27 item(s)
    result[0, 0, 0] = numpy.sum(x24 * (x0 * (3.0 * x10 + x12 * x13 + x9) + x14 * x17))
    result[0, 0, 1] = numpy.sum(x26 * x27)
    result[0, 0, 2] = numpy.sum(x27 * x29)
    result[0, 1, 0] = numpy.sum(x30 * x31)
    result[0, 1, 1] = numpy.sum(x22 * x36 * x38)
    result[0, 1, 2] = numpy.sum(x29 * x30 * x39)
    result[0, 2, 0] = numpy.sum(x31 * x40)
    result[0, 2, 1] = numpy.sum(x26 * x39 * x40)
    result[0, 2, 2] = numpy.sum(x18 * x38 * x44)
    result[1, 0, 0] = numpy.sum(x22 * x45 * x49)
    result[1, 0, 1] = numpy.sum(x14 * x54)
    result[1, 0, 2] = numpy.sum(x14 * x29 * x55)
    result[1, 1, 0] = numpy.sum(x11 * x56)
    result[1, 1, 1] = numpy.sum(x53 * (x0 * (3.0 * x33 + x35 * x50 + x47) + x30 * x52))
    result[1, 1, 2] = numpy.sum(x29 * x56)
    result[1, 2, 0] = numpy.sum(x11 * x40 * x55)
    result[1, 2, 1] = numpy.sum(x40 * x54)
    result[1, 2, 2] = numpy.sum(x44 * x49 * x5)
    result[2, 0, 0] = numpy.sum(x45 * x57 * x61)
    result[2, 0, 1] = numpy.sum(x14 * x26 * x63)
    result[2, 0, 2] = numpy.sum(x14 * x67)
    result[2, 1, 0] = numpy.sum(x11 * x30 * x63)
    result[2, 1, 1] = numpy.sum(x21 * x36 * x5 * x61)
    result[2, 1, 2] = numpy.sum(x30 * x67)
    result[2, 2, 0] = numpy.sum(x11 * x68)
    result[2, 2, 1] = numpy.sum(x26 * x68)
    result[2, 2, 2] = numpy.sum(x62 * (x0 * (3.0 * x41 + x43 * x64 + x59) + x40 * x66))
    return result


def diag_quadrupole3d_12(ax, da, A, bx, db, B, C):
    """Cartesian 3D (pd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 3, 6), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = numpy.sqrt(x1)
    x3 = 1.77245385090552 * x2
    x4 = -x1 * (ax * A[0] + bx * B[0])
    x5 = -x4 - B[0]
    x6 = ax * bx * x1
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x4 - C[0]
    x11 = x3 * x7
    x12 = x10 * x11
    x13 = 2.0 * x0
    x14 = x10**2 * x11
    x15 = x0 * x11
    x16 = x14 + x15
    x17 = x16 * x5
    x18 = 2.0 * x5
    x19 = -x4 - A[0]
    x20 = x0 * (x12 * x18 + x14 + 3.0 * x15)
    x21 = x12 * x13
    x22 = x17 + x21
    x23 = x20 + x22 * x5
    x24 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x25 = da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**3.5)
    x26 = x24 * x25
    x27 = 13.0639452948436 * x26
    x28 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x29 = 0.564189583547756 * x1
    x30 = x28 * x29
    x31 = x27 * x30
    x32 = -x1 * (ax * A[1] + bx * B[1])
    x33 = -x32 - B[1]
    x34 = x30 * x33
    x35 = 22.6274169979695 * x26
    x36 = x35 * (x19 * x22 + x20)
    x37 = -x1 * (ax * A[2] + bx * B[2])
    x38 = -x37 - B[2]
    x39 = x30 * x38
    x40 = x16 * x19 + x21
    x41 = 0.318309886183791 * x2
    x42 = x40 * x41
    x43 = x28 * x3
    x44 = x0 * x43
    x45 = x33**2 * x43 + x44
    x46 = x27 * x45
    x47 = x34 * x35
    x48 = x24 * x3
    x49 = x0 * x48
    x50 = x38**2 * x48 + x49
    x51 = x25 * x28
    x52 = 13.0639452948436 * x51
    x53 = x50 * x52
    x54 = -x32 - A[1]
    x55 = x23 * x31
    x56 = 22.6274169979695 * x41
    x57 = x22 * x56
    x58 = x33 * x43
    x59 = x44 + x54 * x58
    x60 = x26 * x59
    x61 = x13 * x58 + x45 * x54
    x62 = x16 * x41
    x63 = x16 * x56
    x64 = -x37 - A[2]
    x65 = x38 * x48
    x66 = x49 + x64 * x65
    x67 = x51 * x66
    x68 = x13 * x65 + x50 * x64
    x69 = x11 * x5**2 + x15
    x70 = x13 * x9 + x19 * x69
    x71 = -x32 - C[1]
    x72 = x43 * x71**2
    x73 = x44 + x72
    x74 = x27 * x41
    x75 = x73 * x74
    x76 = x33 * x73
    x77 = x43 * x71
    x78 = x13 * x77
    x79 = x76 + x78
    x80 = x56 * x79
    x81 = x15 + x19 * x9
    x82 = x26 * x81
    x83 = x56 * x73
    x84 = 2.0 * x33
    x85 = x0 * (3.0 * x44 + x72 + x77 * x84)
    x86 = x33 * x79 + x85
    x87 = x29 * x7
    x88 = x27 * x87
    x89 = x86 * x88
    x90 = x38 * x87
    x91 = x35 * x79
    x92 = x25 * x7
    x93 = 13.0639452948436 * x92
    x94 = x41 * x93
    x95 = x50 * x94
    x96 = x54 * x73 + x78
    x97 = x35 * (x54 * x79 + x85)
    x98 = x29 * x8
    x99 = x25 * x8
    x100 = -x37 - C[2]
    x101 = x100**2 * x48
    x102 = x101 + x49
    x103 = x41 * x52
    x104 = x102 * x103
    x105 = x51 * x81
    x106 = x102 * x56
    x107 = x102 * x38
    x108 = x100 * x48
    x109 = x108 * x13
    x110 = x107 + x109
    x111 = x110 * x56
    x112 = x102 * x94
    x113 = x34 * x92
    x114 = 22.6274169979695 * x110
    x115 = 2.0 * x38
    x116 = x0 * (x101 + x108 * x115 + 3.0 * x49)
    x117 = x110 * x38 + x116
    x118 = x30 * x93
    x119 = x117 * x118
    x120 = x30 * x99
    x121 = x102 * x64 + x109
    x122 = 22.6274169979695 * x110 * x64 + 22.6274169979695 * x116

    # 54 item(s)
    result[0, 0, 0] = numpy.sum(
        x31
        * (
            x0 * (4.0 * x0 * x12 + x13 * (x12 + x9) + 2.0 * x17 + x18 * (x10 * x9 + x15))
            + x19 * x23
        )
    )
    result[0, 0, 1] = numpy.sum(x34 * x36)
    result[0, 0, 2] = numpy.sum(x36 * x39)
    result[0, 0, 3] = numpy.sum(x42 * x46)
    result[0, 0, 4] = numpy.sum(x38 * x40 * x47)
    result[0, 0, 5] = numpy.sum(x42 * x53)
    result[0, 1, 0] = numpy.sum(x54 * x55)
    result[0, 1, 1] = numpy.sum(x57 * x60)
    result[0, 1, 2] = numpy.sum(x22 * x35 * x39 * x54)
    result[0, 1, 3] = numpy.sum(x27 * x61 * x62)
    result[0, 1, 4] = numpy.sum(x38 * x60 * x63)
    result[0, 1, 5] = numpy.sum(x53 * x54 * x62)
    result[0, 2, 0] = numpy.sum(x55 * x64)
    result[0, 2, 1] = numpy.sum(x22 * x47 * x64)
    result[0, 2, 2] = numpy.sum(x57 * x67)
    result[0, 2, 3] = numpy.sum(x46 * x62 * x64)
    result[0, 2, 4] = numpy.sum(x33 * x63 * x67)
    result[0, 2, 5] = numpy.sum(x52 * x62 * x68)
    result[1, 0, 0] = numpy.sum(x70 * x75)
    result[1, 0, 1] = numpy.sum(x80 * x82)
    result[1, 0, 2] = numpy.sum(x38 * x82 * x83)
    result[1, 0, 3] = numpy.sum(x19 * x89)
    result[1, 0, 4] = numpy.sum(x19 * x90 * x91)
    result[1, 0, 5] = numpy.sum(x19 * x73 * x95)
    result[1, 1, 0] = numpy.sum(x69 * x74 * x96)
    result[1, 1, 1] = numpy.sum(x97 * x98)
    result[1, 1, 2] = numpy.sum(x35 * x38 * x96 * x98)
    result[1, 1, 3] = numpy.sum(
        x88
        * (
            x0
            * (x13 * (x58 + x77) + 4.0 * x44 * x71 + 2.0 * x76 + x84 * (x44 + x58 * x71))
            + x54 * x86
        )
    )
    result[1, 1, 4] = numpy.sum(x90 * x97)
    result[1, 1, 5] = numpy.sum(x95 * x96)
    result[1, 2, 0] = numpy.sum(x64 * x69 * x75)
    result[1, 2, 1] = numpy.sum(x64 * x91 * x98)
    result[1, 2, 2] = numpy.sum(x66 * x83 * x99)
    result[1, 2, 3] = numpy.sum(x64 * x89)
    result[1, 2, 4] = numpy.sum(x66 * x80 * x92)
    result[1, 2, 5] = numpy.sum(x68 * x73 * x94)
    result[2, 0, 0] = numpy.sum(x104 * x70)
    result[2, 0, 1] = numpy.sum(x105 * x106 * x33)
    result[2, 0, 2] = numpy.sum(x105 * x111)
    result[2, 0, 3] = numpy.sum(x112 * x19 * x45)
    result[2, 0, 4] = numpy.sum(x113 * x114 * x19)
    result[2, 0, 5] = numpy.sum(x119 * x19)
    result[2, 1, 0] = numpy.sum(x104 * x54 * x69)
    result[2, 1, 1] = numpy.sum(x106 * x59 * x99)
    result[2, 1, 2] = numpy.sum(x114 * x120 * x54)
    result[2, 1, 3] = numpy.sum(x112 * x61)
    result[2, 1, 4] = numpy.sum(x111 * x59 * x92)
    result[2, 1, 5] = numpy.sum(x119 * x54)
    result[2, 2, 0] = numpy.sum(x103 * x121 * x69)
    result[2, 2, 1] = numpy.sum(22.6274169979695 * x121 * x34 * x99)
    result[2, 2, 2] = numpy.sum(x120 * x122)
    result[2, 2, 3] = numpy.sum(x121 * x45 * x94)
    result[2, 2, 4] = numpy.sum(x113 * x122)
    result[2, 2, 5] = numpy.sum(
        x118
        * (
            x0
            * (
                4.0 * x100 * x49
                + 2.0 * x107
                + x115 * (x100 * x65 + x49)
                + x13 * (x108 + x65)
            )
            + x117 * x64
        )
    )
    return result


def diag_quadrupole3d_13(ax, da, A, bx, db, B, C):
    """Cartesian 3D (pf) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 3, 10), dtype=float)

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
    x12 = 2.0 * x3
    x13 = -x2 - C[0]
    x14 = x13 * x8
    x15 = x11 + x12 * x14
    x16 = 2.0 * x0
    x17 = x3 * x8
    x18 = x0 * (x14 + x17)
    x19 = x3 * (x10 + x13 * x17)
    x20 = x13**2 * x8
    x21 = x0 * (x15 + x20)
    x22 = x10 + x20
    x23 = x22 * x3
    x24 = x14 * x16
    x25 = x23 + x24
    x26 = x25 * x3
    x27 = -x2 - A[0]
    x28 = 2.0 * x0 * (2.0 * x0 * x14 + x18 + x19 + x23)
    x29 = x21 + x26
    x30 = x28 + x29 * x3
    x31 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x32 = da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**4.5)
    x33 = 11.6847478934435 * x32
    x34 = x31 * x33
    x35 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x36 = 0.564189583547756 * x1
    x37 = x35 * x36
    x38 = x34 * x37
    x39 = -x1 * (ax * A[1] + bx * B[1])
    x40 = -x39 - B[1]
    x41 = x37 * x40
    x42 = 26.1278905896872 * x32
    x43 = x31 * x42
    x44 = x43 * (x27 * x29 + x28)
    x45 = -x1 * (ax * A[2] + bx * B[2])
    x46 = -x45 - B[2]
    x47 = x37 * x46
    x48 = x21 + x25 * x27
    x49 = 0.318309886183791 * x6
    x50 = x48 * x49
    x51 = x35 * x7
    x52 = x40**2 * x51
    x53 = x0 * x51
    x54 = x52 + x53
    x55 = x43 * x54
    x56 = 45.2548339959391 * x32
    x57 = x31 * x46
    x58 = x56 * x57
    x59 = x31 * x7
    x60 = x46**2 * x59
    x61 = x0 * x59
    x62 = x60 + x61
    x63 = x35 * x42
    x64 = x62 * x63
    x65 = x49 * (x22 * x27 + x24)
    x66 = x40 * x51
    x67 = x16 * x66
    x68 = x40 * x54 + x67
    x69 = x34 * x68
    x70 = x46 * x59
    x71 = x16 * x70
    x72 = x46 * x62 + x71
    x73 = x33 * x35
    x74 = x72 * x73
    x75 = -x39 - A[1]
    x76 = x30 * x38
    x77 = x53 + x66 * x75
    x78 = x43 * x49
    x79 = x29 * x43
    x80 = x54 * x75 + x67
    x81 = x78 * x80
    x82 = x25 * x49
    x83 = x56 * x77
    x84 = 3.0 * x53
    x85 = x0 * (3.0 * x52 + x84) + x68 * x75
    x86 = x22 * x49
    x87 = 0.179587122125167 * x42
    x88 = x22 * x87
    x89 = -x45 - A[2]
    x90 = x61 + x70 * x89
    x91 = x49 * x63
    x92 = x56 * x90
    x93 = x35 * x40
    x94 = x62 * x89 + x71
    x95 = 3.0 * x61
    x96 = x0 * (3.0 * x60 + x95) + x72 * x89
    x97 = x10 + x9
    x98 = x16 * x17
    x99 = x3 * x97 + x98
    x100 = x0 * (x11 + 3.0 * x9) + x27 * x99
    x101 = -x39 - C[1]
    x102 = x101**2 * x51
    x103 = x102 + x53
    x104 = x34 * x49
    x105 = x103 * x104
    x106 = x103 * x40
    x107 = x101 * x51
    x108 = x107 * x16
    x109 = x106 + x108
    x110 = x27 * x97 + x98
    x111 = x110 * x78
    x112 = 2.0 * x40
    x113 = x107 * x112 + x84
    x114 = x0 * (x102 + x113)
    x115 = x109 * x40
    x116 = x114 + x115
    x117 = x10 + x17 * x27
    x118 = x117 * x49
    x119 = x103 * x87
    x120 = x0 * (x107 + x66)
    x121 = x40 * (x101 * x66 + x53)
    x122 = 2.0 * x0 * (2.0 * x101 * x53 + x106 + x120 + x121)
    x123 = x116 * x40 + x122
    x124 = x36 * x5
    x125 = x124 * x34
    x126 = x123 * x125
    x127 = x124 * x43
    x128 = x116 * x127
    x129 = x49 * x5
    x130 = x129 * x42
    x131 = x130 * x62
    x132 = x129 * x33
    x133 = x132 * x72
    x134 = x103 * x75 + x108
    x135 = x109 * x75 + x114
    x136 = x78 * x97
    x137 = x127 * (x116 * x75 + x122)
    x138 = x129 * x3
    x139 = x130 * x94
    x140 = -x45 - C[2]
    x141 = x140**2 * x59
    x142 = x141 + x61
    x143 = x49 * x73
    x144 = x142 * x143
    x145 = x110 * x91
    x146 = x142 * x46
    x147 = x140 * x59
    x148 = x147 * x16
    x149 = x146 + x148
    x150 = x142 * x87
    x151 = 2.0 * x46
    x152 = x147 * x151 + x95
    x153 = x0 * (x141 + x152)
    x154 = x149 * x46
    x155 = x153 + x154
    x156 = x132 * x142
    x157 = x130 * x149
    x158 = x37 * x5
    x159 = x158 * x42
    x160 = x155 * x159
    x161 = x0 * (x147 + x70)
    x162 = x46 * (x140 * x70 + x61)
    x163 = 2.0 * x0 * (2.0 * x140 * x61 + x146 + x161 + x162)
    x164 = x155 * x46 + x163
    x165 = x158 * x33
    x166 = x164 * x165
    x167 = x91 * x97
    x168 = x142 * x89 + x148
    x169 = x149 * x89 + x153
    x170 = x130 * x54
    x171 = x159 * (x155 * x89 + x163)

    # 90 item(s)
    result[0, 0, 0] = numpy.sum(
        x38
        * (
            x0 * (x12 * (x18 + x19) + x16 * (x15 + x9) + 3.0 * x21 + 3.0 * x26)
            + x27 * x30
        )
    )
    result[0, 0, 1] = numpy.sum(x41 * x44)
    result[0, 0, 2] = numpy.sum(x44 * x47)
    result[0, 0, 3] = numpy.sum(x50 * x55)
    result[0, 0, 4] = numpy.sum(x41 * x48 * x58)
    result[0, 0, 5] = numpy.sum(x50 * x64)
    result[0, 0, 6] = numpy.sum(x65 * x69)
    result[0, 0, 7] = numpy.sum(x46 * x55 * x65)
    result[0, 0, 8] = numpy.sum(x40 * x64 * x65)
    result[0, 0, 9] = numpy.sum(x65 * x74)
    result[0, 1, 0] = numpy.sum(x75 * x76)
    result[0, 1, 1] = numpy.sum(x29 * x77 * x78)
    result[0, 1, 2] = numpy.sum(x47 * x75 * x79)
    result[0, 1, 3] = numpy.sum(x25 * x81)
    result[0, 1, 4] = numpy.sum(x57 * x82 * x83)
    result[0, 1, 5] = numpy.sum(x64 * x75 * x82)
    result[0, 1, 6] = numpy.sum(x34 * x85 * x86)
    result[0, 1, 7] = numpy.sum(x22 * x46 * x81)
    result[0, 1, 8] = numpy.sum(x62 * x77 * x88)
    result[0, 1, 9] = numpy.sum(x74 * x75 * x86)
    result[0, 2, 0] = numpy.sum(x76 * x89)
    result[0, 2, 1] = numpy.sum(x41 * x79 * x89)
    result[0, 2, 2] = numpy.sum(x29 * x90 * x91)
    result[0, 2, 3] = numpy.sum(x55 * x82 * x89)
    result[0, 2, 4] = numpy.sum(x82 * x92 * x93)
    result[0, 2, 5] = numpy.sum(x25 * x91 * x94)
    result[0, 2, 6] = numpy.sum(x69 * x86 * x89)
    result[0, 2, 7] = numpy.sum(x54 * x88 * x90)
    result[0, 2, 8] = numpy.sum(x40 * x63 * x86 * x94)
    result[0, 2, 9] = numpy.sum(x73 * x86 * x96)
    result[1, 0, 0] = numpy.sum(x100 * x105)
    result[1, 0, 1] = numpy.sum(x109 * x111)
    result[1, 0, 2] = numpy.sum(x103 * x111 * x46)
    result[1, 0, 3] = numpy.sum(x116 * x117 * x78)
    result[1, 0, 4] = numpy.sum(x109 * x118 * x58)
    result[1, 0, 5] = numpy.sum(x117 * x119 * x62)
    result[1, 0, 6] = numpy.sum(x126 * x27)
    result[1, 0, 7] = numpy.sum(x128 * x27 * x46)
    result[1, 0, 8] = numpy.sum(x109 * x131 * x27)
    result[1, 0, 9] = numpy.sum(x103 * x133 * x27)
    result[1, 1, 0] = numpy.sum(x104 * x134 * x99)
    result[1, 1, 1] = numpy.sum(x135 * x136)
    result[1, 1, 2] = numpy.sum(x134 * x136 * x46)
    result[1, 1, 3] = numpy.sum(x137 * x3)
    result[1, 1, 4] = numpy.sum(x124 * x135 * x3 * x58)
    result[1, 1, 5] = numpy.sum(x131 * x134 * x3)
    result[1, 1, 6] = numpy.sum(
        x125
        * (
            x0 * (x112 * (x120 + x121) + 3.0 * x114 + 3.0 * x115 + x16 * (x113 + x52))
            + x123 * x75
        )
    )
    result[1, 1, 7] = numpy.sum(x137 * x46)
    result[1, 1, 8] = numpy.sum(x131 * x135)
    result[1, 1, 9] = numpy.sum(x133 * x134)
    result[1, 2, 0] = numpy.sum(x105 * x89 * x99)
    result[1, 2, 1] = numpy.sum(x109 * x136 * x89)
    result[1, 2, 2] = numpy.sum(x119 * x90 * x97)
    result[1, 2, 3] = numpy.sum(x128 * x3 * x89)
    result[1, 2, 4] = numpy.sum(x109 * x138 * x92)
    result[1, 2, 5] = numpy.sum(x103 * x139 * x3)
    result[1, 2, 6] = numpy.sum(x126 * x89)
    result[1, 2, 7] = numpy.sum(x116 * x130 * x90)
    result[1, 2, 8] = numpy.sum(x109 * x139)
    result[1, 2, 9] = numpy.sum(x103 * x132 * x96)
    result[2, 0, 0] = numpy.sum(x100 * x144)
    result[2, 0, 1] = numpy.sum(x142 * x145 * x40)
    result[2, 0, 2] = numpy.sum(x145 * x149)
    result[2, 0, 3] = numpy.sum(x117 * x150 * x54)
    result[2, 0, 4] = numpy.sum(x118 * x149 * x56 * x93)
    result[2, 0, 5] = numpy.sum(x117 * x155 * x91)
    result[2, 0, 6] = numpy.sum(x156 * x27 * x68)
    result[2, 0, 7] = numpy.sum(x157 * x27 * x54)
    result[2, 0, 8] = numpy.sum(x160 * x27 * x40)
    result[2, 0, 9] = numpy.sum(x166 * x27)
    result[2, 1, 0] = numpy.sum(x144 * x75 * x99)
    result[2, 1, 1] = numpy.sum(x150 * x77 * x97)
    result[2, 1, 2] = numpy.sum(x149 * x167 * x75)
    result[2, 1, 3] = numpy.sum(x130 * x142 * x3 * x80)
    result[2, 1, 4] = numpy.sum(x138 * x149 * x83)
    result[2, 1, 5] = numpy.sum(x160 * x3 * x75)
    result[2, 1, 6] = numpy.sum(x156 * x85)
    result[2, 1, 7] = numpy.sum(x157 * x80)
    result[2, 1, 8] = numpy.sum(x130 * x155 * x77)
    result[2, 1, 9] = numpy.sum(x166 * x75)
    result[2, 2, 0] = numpy.sum(x143 * x168 * x99)
    result[2, 2, 1] = numpy.sum(x167 * x168 * x40)
    result[2, 2, 2] = numpy.sum(x167 * x169)
    result[2, 2, 3] = numpy.sum(x168 * x170 * x3)
    result[2, 2, 4] = numpy.sum(x158 * x169 * x3 * x40 * x56)
    result[2, 2, 5] = numpy.sum(x171 * x3)
    result[2, 2, 6] = numpy.sum(x132 * x168 * x68)
    result[2, 2, 7] = numpy.sum(x169 * x170)
    result[2, 2, 8] = numpy.sum(x171 * x40)
    result[2, 2, 9] = numpy.sum(
        x165
        * (
            x0 * (x151 * (x161 + x162) + 3.0 * x153 + 3.0 * x154 + x16 * (x152 + x60))
            + x164 * x89
        )
    )
    return result


def diag_quadrupole3d_14(ax, da, A, bx, db, B, C):
    """Cartesian 3D (pg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 3, 15), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x2 - C[0]
    x11 = x10 * x8
    x12 = x0 * (x11 + x9)
    x13 = x0 * x8
    x14 = x3 * (x10 * x9 + x13)
    x15 = x3**2 * x8
    x16 = x13 + x15
    x17 = x16 * x3
    x18 = 2.0 * x0
    x19 = x18 * x9
    x20 = x17 + x19
    x21 = 3.0 * x13
    x22 = 2.0 * x3
    x23 = x11 * x22 + x21
    x24 = x0 * (x15 + x23)
    x25 = x3 * (x12 + x14)
    x26 = x10**2 * x8
    x27 = x13 + x26
    x28 = x27 * x3
    x29 = 2.0 * x0 * (2.0 * x10 * x13 + x12 + x14 + x28)
    x30 = x0 * (x23 + x26)
    x31 = x11 * x18
    x32 = x28 + x31
    x33 = x3 * x32
    x34 = x30 + x33
    x35 = x3 * x34
    x36 = -x2 - A[0]
    x37 = x0 * (2.0 * x24 + 2.0 * x25 + 3.0 * x30 + 3.0 * x33)
    x38 = x29 + x35
    x39 = x3 * x38 + x37
    x40 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x41 = da * db * numpy.sqrt(ax**2.5) * numpy.sqrt(bx**5.5)
    x42 = x40 * x41
    x43 = 8.83283915958214 * x42
    x44 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x45 = 0.564189583547756 * x1
    x46 = x44 * x45
    x47 = x43 * x46
    x48 = -x1 * (ax * A[1] + bx * B[1])
    x49 = -x48 - B[1]
    x50 = x46 * x49
    x51 = 23.3694957868871 * x42
    x52 = x51 * (x36 * x38 + x37)
    x53 = -x1 * (ax * A[2] + bx * B[2])
    x54 = -x53 - B[2]
    x55 = x46 * x54
    x56 = x29 + x34 * x36
    x57 = x42 * x56
    x58 = x44 * x7
    x59 = x49**2 * x58
    x60 = x0 * x58
    x61 = x59 + x60
    x62 = 30.169889330626 * x61
    x63 = 0.318309886183791 * x6
    x64 = x62 * x63
    x65 = x40 * x7
    x66 = x54**2 * x65
    x67 = x0 * x65
    x68 = x66 + x67
    x69 = 30.169889330626 * x63
    x70 = x41 * x44
    x71 = x69 * x70
    x72 = x68 * x71
    x73 = x49 * x61
    x74 = x49 * x58
    x75 = x18 * x74
    x76 = x73 + x75
    x77 = 23.3694957868871 * x76
    x78 = x63 * (x30 + x32 * x36)
    x79 = x42 * x78
    x80 = 52.2557811793745 * x61
    x81 = x70 * x78
    x82 = 52.2557811793745 * x49
    x83 = x54 * x68
    x84 = x54 * x65
    x85 = x18 * x84
    x86 = x83 + x85
    x87 = 23.3694957868871 * x86
    x88 = x27 * x36 + x31
    x89 = x63 * x88
    x90 = 3.0 * x60
    x91 = x0 * (3.0 * x59 + x90)
    x92 = x49 * x76 + x91
    x93 = x43 * x92
    x94 = x42 * x54
    x95 = 0.179587122125167 * x41
    x96 = x68 * x95
    x97 = 3.0 * x67
    x98 = x0 * (3.0 * x66 + x97)
    x99 = x54 * x86 + x98
    x100 = 8.83283915958214 * x70
    x101 = x100 * x99
    x102 = -x48 - A[1]
    x103 = x39 * x47
    x104 = x102 * x74 + x60
    x105 = 23.3694957868871 * x63
    x106 = x105 * x42
    x107 = x38 * x51
    x108 = x102 * x61 + x75
    x109 = x42 * x69
    x110 = 52.2557811793745 * x104
    x111 = x42 * x63
    x112 = x111 * x54
    x113 = x102 * x76 + x91
    x114 = x106 * x113
    x115 = 52.2557811793745 * x32
    x116 = x63 * x70
    x117 = 4.0 * x0 * (2.0 * x49 * x60 + x73) + x102 * x92
    x118 = x27 * x63
    x119 = 30.169889330626 * x96
    x120 = x27 * x95
    x121 = -x53 - A[2]
    x122 = x121 * x84 + x67
    x123 = x105 * x70
    x124 = x116 * x82
    x125 = x121 * x68 + x85
    x126 = x122 * x95
    x127 = x121 * x86 + x98
    x128 = x123 * x127
    x129 = 4.0 * x0 * (2.0 * x54 * x67 + x83) + x121 * x99
    x130 = x0 * (3.0 * x15 + x21)
    x131 = x130 + x20 * x3
    x132 = 4.0 * x0 * (2.0 * x13 * x3 + x17) + x131 * x36
    x133 = -x48 - C[1]
    x134 = x133**2 * x58
    x135 = x134 + x60
    x136 = x43 * x63
    x137 = x135 * x136
    x138 = x135 * x49
    x139 = x133 * x58
    x140 = x139 * x18
    x141 = x138 + x140
    x142 = x130 + x20 * x36
    x143 = x106 * x142
    x144 = 2.0 * x49
    x145 = x139 * x144 + x90
    x146 = x0 * (x134 + x145)
    x147 = x141 * x49
    x148 = x146 + x147
    x149 = x16 * x36 + x19
    x150 = 52.2557811793745 * x141
    x151 = x0 * (x139 + x74)
    x152 = x49 * (x133 * x74 + x60)
    x153 = 2.0 * x0 * (2.0 * x133 * x60 + x138 + x151 + x152)
    x154 = x148 * x49
    x155 = x153 + x154
    x156 = x13 + x36 * x9
    x157 = 52.2557811793745 * x112
    x158 = x135 * x95
    x159 = x0 * (x145 + x59)
    x160 = x49 * (x151 + x152)
    x161 = x0 * (3.0 * x146 + 3.0 * x147 + 2.0 * x159 + 2.0 * x160)
    x162 = x155 * x49 + x161
    x163 = x45 * x5
    x164 = x163 * x43
    x165 = x162 * x164
    x166 = x163 * x51
    x167 = x155 * x166
    x168 = x41 * x5
    x169 = x168 * x69
    x170 = x169 * x68
    x171 = x168 * x63
    x172 = x171 * x87
    x173 = 8.83283915958214 * x171
    x174 = x173 * x99
    x175 = x102 * x135 + x140
    x176 = x102 * x141 + x146
    x177 = 23.3694957868871 * x20
    x178 = x111 * x177
    x179 = x102 * x148 + x153
    x180 = 30.169889330626 * x16
    x181 = x111 * x180
    x182 = x166 * (x102 * x155 + x161)
    x183 = x171 * x3
    x184 = 52.2557811793745 * x183
    x185 = x105 * x168
    x186 = x127 * x185
    x187 = -x53 - C[2]
    x188 = x187**2 * x65
    x189 = x188 + x67
    x190 = x100 * x63
    x191 = x189 * x190
    x192 = x123 * x142
    x193 = x189 * x54
    x194 = x187 * x65
    x195 = x18 * x194
    x196 = x193 + x195
    x197 = x189 * x95
    x198 = 2.0 * x54
    x199 = x194 * x198 + x97
    x200 = x0 * (x188 + x199)
    x201 = x196 * x54
    x202 = x200 + x201
    x203 = x196 * x95
    x204 = x0 * (x194 + x84)
    x205 = x54 * (x187 * x84 + x67)
    x206 = 2.0 * x0 * (2.0 * x187 * x67 + x193 + x204 + x205)
    x207 = x202 * x54
    x208 = x206 + x207
    x209 = x173 * x189
    x210 = x171 * x77
    x211 = x168 * x64
    x212 = x168 * x46
    x213 = 23.3694957868871 * x212
    x214 = x208 * x213
    x215 = x0 * (x199 + x66)
    x216 = x54 * (x204 + x205)
    x217 = x0 * (3.0 * x200 + 3.0 * x201 + 2.0 * x215 + 2.0 * x216)
    x218 = x208 * x54 + x217
    x219 = 8.83283915958214 * x212
    x220 = x218 * x219
    x221 = x116 * x177
    x222 = x116 * x180
    x223 = x113 * x185
    x224 = x121 * x189 + x195
    x225 = x121 * x196 + x200
    x226 = x121 * x202 + x206
    x227 = x213 * (x121 * x208 + x217)

    # 135 item(s)
    result[0, 0, 0] = numpy.sum(
        x47
        * (
            x0
            * (
                x18 * (3.0 * x12 + 3.0 * x14 + x20)
                + x22 * (x24 + x25)
                + 4.0 * x29
                + 4.0 * x35
            )
            + x36 * x39
        )
    )
    result[0, 0, 1] = numpy.sum(x50 * x52)
    result[0, 0, 2] = numpy.sum(x52 * x55)
    result[0, 0, 3] = numpy.sum(x57 * x64)
    result[0, 0, 4] = numpy.sum(52.2557811793745 * x50 * x54 * x57)
    result[0, 0, 5] = numpy.sum(x56 * x72)
    result[0, 0, 6] = numpy.sum(x77 * x79)
    result[0, 0, 7] = numpy.sum(x54 * x79 * x80)
    result[0, 0, 8] = numpy.sum(x68 * x81 * x82)
    result[0, 0, 9] = numpy.sum(x81 * x87)
    result[0, 0, 10] = numpy.sum(x89 * x93)
    result[0, 0, 11] = numpy.sum(x77 * x89 * x94)
    result[0, 0, 12] = numpy.sum(x62 * x88 * x96)
    result[0, 0, 13] = numpy.sum(x49 * x70 * x87 * x89)
    result[0, 0, 14] = numpy.sum(x101 * x89)
    result[0, 1, 0] = numpy.sum(x102 * x103)
    result[0, 1, 1] = numpy.sum(x104 * x106 * x38)
    result[0, 1, 2] = numpy.sum(x102 * x107 * x55)
    result[0, 1, 3] = numpy.sum(x108 * x109 * x34)
    result[0, 1, 4] = numpy.sum(x110 * x112 * x34)
    result[0, 1, 5] = numpy.sum(x102 * x34 * x72)
    result[0, 1, 6] = numpy.sum(x114 * x32)
    result[0, 1, 7] = numpy.sum(x108 * x112 * x115)
    result[0, 1, 8] = numpy.sum(x104 * x115 * x96)
    result[0, 1, 9] = numpy.sum(x102 * x116 * x32 * x87)
    result[0, 1, 10] = numpy.sum(x117 * x118 * x43)
    result[0, 1, 11] = numpy.sum(x114 * x27 * x54)
    result[0, 1, 12] = numpy.sum(x108 * x119 * x27)
    result[0, 1, 13] = numpy.sum(x104 * x120 * x87)
    result[0, 1, 14] = numpy.sum(x101 * x102 * x118)
    result[0, 2, 0] = numpy.sum(x103 * x121)
    result[0, 2, 1] = numpy.sum(x107 * x121 * x50)
    result[0, 2, 2] = numpy.sum(x122 * x123 * x38)
    result[0, 2, 3] = numpy.sum(x121 * x34 * x42 * x64)
    result[0, 2, 4] = numpy.sum(x122 * x124 * x34)
    result[0, 2, 5] = numpy.sum(x125 * x34 * x71)
    result[0, 2, 6] = numpy.sum(x111 * x121 * x32 * x77)
    result[0, 2, 7] = numpy.sum(x115 * x126 * x61)
    result[0, 2, 8] = numpy.sum(x115 * x116 * x125 * x49)
    result[0, 2, 9] = numpy.sum(x128 * x32)
    result[0, 2, 10] = numpy.sum(x118 * x121 * x93)
    result[0, 2, 11] = numpy.sum(x120 * x122 * x77)
    result[0, 2, 12] = numpy.sum(x120 * x125 * x62)
    result[0, 2, 13] = numpy.sum(x128 * x27 * x49)
    result[0, 2, 14] = numpy.sum(x100 * x118 * x129)
    result[1, 0, 0] = numpy.sum(x132 * x137)
    result[1, 0, 1] = numpy.sum(x141 * x143)
    result[1, 0, 2] = numpy.sum(x135 * x143 * x54)
    result[1, 0, 3] = numpy.sum(x109 * x148 * x149)
    result[1, 0, 4] = numpy.sum(x112 * x149 * x150)
    result[1, 0, 5] = numpy.sum(x119 * x135 * x149)
    result[1, 0, 6] = numpy.sum(x106 * x155 * x156)
    result[1, 0, 7] = numpy.sum(x148 * x156 * x157)
    result[1, 0, 8] = numpy.sum(x150 * x156 * x96)
    result[1, 0, 9] = numpy.sum(x156 * x158 * x87)
    result[1, 0, 10] = numpy.sum(x165 * x36)
    result[1, 0, 11] = numpy.sum(x167 * x36 * x54)
    result[1, 0, 12] = numpy.sum(x148 * x170 * x36)
    result[1, 0, 13] = numpy.sum(x141 * x172 * x36)
    result[1, 0, 14] = numpy.sum(x135 * x174 * x36)
    result[1, 1, 0] = numpy.sum(x131 * x136 * x175)
    result[1, 1, 1] = numpy.sum(x176 * x178)
    result[1, 1, 2] = numpy.sum(x175 * x178 * x54)
    result[1, 1, 3] = numpy.sum(x179 * x181)
    result[1, 1, 4] = numpy.sum(x157 * x16 * x176)
    result[1, 1, 5] = numpy.sum(x119 * x16 * x175)
    result[1, 1, 6] = numpy.sum(x182 * x3)
    result[1, 1, 7] = numpy.sum(52.2557811793745 * x163 * x179 * x3 * x94)
    result[1, 1, 8] = numpy.sum(x176 * x184 * x68)
    result[1, 1, 9] = numpy.sum(x172 * x175 * x3)
    result[1, 1, 10] = numpy.sum(
        x164
        * (
            x0
            * (
                x144 * (x159 + x160)
                + 4.0 * x153
                + 4.0 * x154
                + x18 * (3.0 * x151 + 3.0 * x152 + x76)
            )
            + x102 * x162
        )
    )
    result[1, 1, 11] = numpy.sum(x182 * x54)
    result[1, 1, 12] = numpy.sum(x170 * x179)
    result[1, 1, 13] = numpy.sum(x172 * x176)
    result[1, 1, 14] = numpy.sum(x174 * x175)
    result[1, 2, 0] = numpy.sum(x121 * x131 * x137)
    result[1, 2, 1] = numpy.sum(x121 * x141 * x178)
    result[1, 2, 2] = numpy.sum(x126 * x135 * x177)
    result[1, 2, 3] = numpy.sum(x121 * x148 * x181)
    result[1, 2, 4] = numpy.sum(x126 * x150 * x16)
    result[1, 2, 5] = numpy.sum(x125 * x158 * x180)
    result[1, 2, 6] = numpy.sum(x121 * x167 * x3)
    result[1, 2, 7] = numpy.sum(x122 * x148 * x184)
    result[1, 2, 8] = numpy.sum(x125 * x150 * x183)
    result[1, 2, 9] = numpy.sum(x135 * x186 * x3)
    result[1, 2, 10] = numpy.sum(x121 * x165)
    result[1, 2, 11] = numpy.sum(x122 * x155 * x185)
    result[1, 2, 12] = numpy.sum(x125 * x148 * x169)
    result[1, 2, 13] = numpy.sum(x141 * x186)
    result[1, 2, 14] = numpy.sum(x129 * x135 * x173)
    result[2, 0, 0] = numpy.sum(x132 * x191)
    result[2, 0, 1] = numpy.sum(x189 * x192 * x49)
    result[2, 0, 2] = numpy.sum(x192 * x196)
    result[2, 0, 3] = numpy.sum(x149 * x197 * x62)
    result[2, 0, 4] = numpy.sum(x124 * x149 * x196)
    result[2, 0, 5] = numpy.sum(x149 * x202 * x71)
    result[2, 0, 6] = numpy.sum(x156 * x197 * x77)
    result[2, 0, 7] = numpy.sum(x156 * x203 * x80)
    result[2, 0, 8] = numpy.sum(x124 * x156 * x202)
    result[2, 0, 9] = numpy.sum(x123 * x156 * x208)
    result[2, 0, 10] = numpy.sum(x209 * x36 * x92)
    result[2, 0, 11] = numpy.sum(x196 * x210 * x36)
    result[2, 0, 12] = numpy.sum(x202 * x211 * x36)
    result[2, 0, 13] = numpy.sum(x214 * x36 * x49)
    result[2, 0, 14] = numpy.sum(x220 * x36)
    result[2, 1, 0] = numpy.sum(x102 * x131 * x191)
    result[2, 1, 1] = numpy.sum(x104 * x177 * x197)
    result[2, 1, 2] = numpy.sum(x102 * x196 * x221)
    result[2, 1, 3] = numpy.sum(x108 * x180 * x197)
    result[2, 1, 4] = numpy.sum(x110 * x16 * x203)
    result[2, 1, 5] = numpy.sum(x102 * x202 * x222)
    result[2, 1, 6] = numpy.sum(x189 * x223 * x3)
    result[2, 1, 7] = numpy.sum(x108 * x184 * x196)
    result[2, 1, 8] = numpy.sum(x110 * x183 * x202)
    result[2, 1, 9] = numpy.sum(x102 * x214 * x3)
    result[2, 1, 10] = numpy.sum(x117 * x209)
    result[2, 1, 11] = numpy.sum(x196 * x223)
    result[2, 1, 12] = numpy.sum(x108 * x169 * x202)
    result[2, 1, 13] = numpy.sum(x104 * x185 * x208)
    result[2, 1, 14] = numpy.sum(x102 * x220)
    result[2, 2, 0] = numpy.sum(x131 * x190 * x224)
    result[2, 2, 1] = numpy.sum(x221 * x224 * x49)
    result[2, 2, 2] = numpy.sum(x221 * x225)
    result[2, 2, 3] = numpy.sum(x16 * x224 * x62 * x95)
    result[2, 2, 4] = numpy.sum(x124 * x16 * x225)
    result[2, 2, 5] = numpy.sum(x222 * x226)
    result[2, 2, 6] = numpy.sum(x210 * x224 * x3)
    result[2, 2, 7] = numpy.sum(x183 * x225 * x80)
    result[2, 2, 8] = numpy.sum(x212 * x226 * x3 * x82)
    result[2, 2, 9] = numpy.sum(x227 * x3)
    result[2, 2, 10] = numpy.sum(x173 * x224 * x92)
    result[2, 2, 11] = numpy.sum(x210 * x225)
    result[2, 2, 12] = numpy.sum(x211 * x226)
    result[2, 2, 13] = numpy.sum(x227 * x49)
    result[2, 2, 14] = numpy.sum(
        x219
        * (
            x0
            * (
                x18 * (3.0 * x204 + 3.0 * x205 + x86)
                + x198 * (x215 + x216)
                + 4.0 * x206
                + 4.0 * x207
            )
            + x121 * x218
        )
    )
    return result


def diag_quadrupole3d_20(ax, da, A, bx, db, B, C):
    """Cartesian 3D (ds) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 6, 1), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - C[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3**2 * x8
    x10 = x0 * x8
    x11 = -x2 - A[0]
    x12 = 2.0 * x3
    x13 = x10 + x9
    x14 = x10 * x12 + x11 * x13
    x15 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x16 = numpy.sqrt(ax**3.5)
    x17 = numpy.sqrt(bx**1.5)
    x18 = 6.53197264742181 * da * db * x16 * x17
    x19 = x15 * x18
    x20 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x21 = 0.564189583547756 * x1
    x22 = x20 * x21
    x23 = -x1 * (ax * A[1] + bx * B[1])
    x24 = -x23 - A[1]
    x25 = 11.3137084989848 * da * db * x16 * x17
    x26 = x15 * x25
    x27 = x14 * x22 * x26
    x28 = -x1 * (ax * A[2] + bx * B[2])
    x29 = -x28 - A[2]
    x30 = x20 * x7
    x31 = x0 * x30
    x32 = x24**2 * x30 + x31
    x33 = 0.318309886183791 * x6
    x34 = x13 * x33
    x35 = x26 * x29
    x36 = x15 * x7
    x37 = x0 * x36
    x38 = x18 * (x29**2 * x36 + x37)
    x39 = x10 + x11**2 * x8
    x40 = -x23 - C[1]
    x41 = x30 * x40**2
    x42 = x31 + x41
    x43 = x33 * x42
    x44 = 2.0 * x40
    x45 = x24 * x42 + x31 * x44
    x46 = x21 * x5
    x47 = x26 * x45 * x46
    x48 = -x28 - C[2]
    x49 = x36 * x48**2
    x50 = x37 + x49
    x51 = x18 * x33 * x50
    x52 = x22 * x5
    x53 = 2.0 * x48
    x54 = x29 * x50 + x37 * x53
    x55 = x25 * x52 * x54

    # 18 item(s)
    result[0, 0, 0] = numpy.sum(
        x19 * x22 * (x0 * (3.0 * x10 + x11 * x12 * x8 + x9) + x11 * x14)
    )
    result[0, 1, 0] = numpy.sum(x24 * x27)
    result[0, 2, 0] = numpy.sum(x27 * x29)
    result[0, 3, 0] = numpy.sum(x19 * x32 * x34)
    result[0, 4, 0] = numpy.sum(x13 * x22 * x24 * x35)
    result[0, 5, 0] = numpy.sum(x20 * x34 * x38)
    result[1, 0, 0] = numpy.sum(x19 * x39 * x43)
    result[1, 1, 0] = numpy.sum(x11 * x47)
    result[1, 2, 0] = numpy.sum(x11 * x35 * x42 * x46)
    result[1, 3, 0] = numpy.sum(
        x19 * x46 * (x0 * (x24 * x30 * x44 + 3.0 * x31 + x41) + x24 * x45)
    )
    result[1, 4, 0] = numpy.sum(x29 * x47)
    result[1, 5, 0] = numpy.sum(x38 * x43 * x5)
    result[2, 0, 0] = numpy.sum(x20 * x39 * x51)
    result[2, 1, 0] = numpy.sum(x11 * x24 * x25 * x50 * x52)
    result[2, 2, 0] = numpy.sum(x11 * x55)
    result[2, 3, 0] = numpy.sum(x32 * x5 * x51)
    result[2, 4, 0] = numpy.sum(x24 * x55)
    result[2, 5, 0] = numpy.sum(
        x18 * x52 * (x0 * (x29 * x36 * x53 + 3.0 * x37 + x49) + x29 * x54)
    )
    return result


def diag_quadrupole3d_21(ax, da, A, bx, db, B, C):
    """Cartesian 3D (dp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 6, 3), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = -x2 - C[0]
    x5 = ax * bx * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(x1)
    x8 = 1.77245385090552 * x7
    x9 = x6 * x8
    x10 = x4**2 * x9
    x11 = x0 * x9
    x12 = x10 + x11
    x13 = x12 * x3
    x14 = -x2 - B[0]
    x15 = x12 * x14
    x16 = x14 * x9
    x17 = x4 * x9
    x18 = 2.0 * x0
    x19 = x16 * x4
    x20 = 2.0 * x3
    x21 = x10 + 3.0 * x11
    x22 = x17 * x18
    x23 = x15 + x22
    x24 = x0 * (2.0 * x19 + x21) + x23 * x3
    x25 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x26 = da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**2.5)
    x27 = x25 * x26
    x28 = 13.0639452948436 * x27
    x29 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x30 = 0.564189583547756 * x1
    x31 = x29 * x30
    x32 = x28 * x31
    x33 = -x1 * (ax * A[1] + bx * B[1])
    x34 = -x33 - B[1]
    x35 = x13 + x22
    x36 = x32 * (x0 * (x17 * x20 + x21) + x3 * x35)
    x37 = -x1 * (ax * A[2] + bx * B[2])
    x38 = -x37 - B[2]
    x39 = -x33 - A[1]
    x40 = x31 * x39
    x41 = 22.6274169979695 * x27
    x42 = x24 * x41
    x43 = x35 * x41
    x44 = x0 * x8
    x45 = x29 * x44
    x46 = x29 * x8
    x47 = x39 * x46
    x48 = x34 * x47 + x45
    x49 = 0.318309886183791 * x7
    x50 = x48 * x49
    x51 = -x37 - A[2]
    x52 = x31 * x51
    x53 = x26 * x49
    x54 = x29 * x53
    x55 = x25 * x44
    x56 = x25 * x8
    x57 = x51 * x56
    x58 = x38 * x57 + x55
    x59 = 22.6274169979695 * x58
    x60 = x54 * x59
    x61 = x39**2 * x46 + x45
    x62 = x28 * x49
    x63 = x61 * x62
    x64 = x34 * x46
    x65 = x0 * (x47 + x64) + x39 * x48
    x66 = x41 * x51
    x67 = x51**2 * x56 + x55
    x68 = 13.0639452948436 * x54
    x69 = x67 * x68
    x70 = x38 * x56
    x71 = x0 * (x57 + x70) + x51 * x58
    x72 = x11 + x16 * x3
    x73 = x0 * (x16 + x3 * x9) + x3 * x72
    x74 = -x33 - C[1]
    x75 = x46 * x74**2
    x76 = x45 + x75
    x77 = x62 * x76
    x78 = x34 * x76
    x79 = x46 * x74
    x80 = x18 * x79
    x81 = x78 + x80
    x82 = x11 + x3**2 * x9
    x83 = x39 * x76
    x84 = x80 + x83
    x85 = x41 * x84
    x86 = x49 * x72
    x87 = x64 * x74
    x88 = 3.0 * x45 + x75
    x89 = x0 * (2.0 * x87 + x88) + x39 * x81
    x90 = x30 * x6
    x91 = x89 * x90
    x92 = x3 * x90
    x93 = x53 * x6
    x94 = x59 * x93
    x95 = 2.0 * x39
    x96 = x28 * x90
    x97 = x96 * (x0 * (x79 * x95 + x88) + x39 * x84)
    x98 = 13.0639452948436 * x93
    x99 = x67 * x98
    x100 = -x37 - C[2]
    x101 = x100**2 * x56
    x102 = x101 + x55
    x103 = x102 * x68
    x104 = x102 * x38
    x105 = x100 * x56
    x106 = x105 * x18
    x107 = x104 + x106
    x108 = x54 * x72
    x109 = 22.6274169979695 * x102
    x110 = x26 * x6
    x111 = x110 * x50
    x112 = x110 * x31
    x113 = x112 * x3
    x114 = x102 * x51
    x115 = x106 + x114
    x116 = 22.6274169979695 * x115
    x117 = x100 * x70
    x118 = x101 + 3.0 * x55
    x119 = x0 * (2.0 * x117 + x118) + x107 * x51
    x120 = 22.6274169979695 * x119
    x121 = x102 * x98
    x122 = x112 * x39
    x123 = 2.0 * x51
    x124 = 13.0639452948436 * x112
    x125 = x124 * (x0 * (x105 * x123 + x118) + x115 * x51)

    # 54 item(s)
    result[0, 0, 0] = numpy.sum(
        x32
        * (
            x0 * (4.0 * x11 * x4 + x13 + x15 + x18 * (x16 + x17) + x20 * (x11 + x19))
            + x24 * x3
        )
    )
    result[0, 0, 1] = numpy.sum(x34 * x36)
    result[0, 0, 2] = numpy.sum(x36 * x38)
    result[0, 1, 0] = numpy.sum(x40 * x42)
    result[0, 1, 1] = numpy.sum(x43 * x50)
    result[0, 1, 2] = numpy.sum(x38 * x40 * x43)
    result[0, 2, 0] = numpy.sum(x42 * x52)
    result[0, 2, 1] = numpy.sum(x34 * x43 * x52)
    result[0, 2, 2] = numpy.sum(x35 * x60)
    result[0, 3, 0] = numpy.sum(x23 * x63)
    result[0, 3, 1] = numpy.sum(x12 * x62 * x65)
    result[0, 3, 2] = numpy.sum(x12 * x38 * x63)
    result[0, 4, 0] = numpy.sum(x23 * x40 * x66)
    result[0, 4, 1] = numpy.sum(x12 * x50 * x66)
    result[0, 4, 2] = numpy.sum(x12 * x39 * x60)
    result[0, 5, 0] = numpy.sum(x23 * x69)
    result[0, 5, 1] = numpy.sum(x12 * x34 * x69)
    result[0, 5, 2] = numpy.sum(x12 * x68 * x71)
    result[1, 0, 0] = numpy.sum(x73 * x77)
    result[1, 0, 1] = numpy.sum(x62 * x81 * x82)
    result[1, 0, 2] = numpy.sum(x38 * x77 * x82)
    result[1, 1, 0] = numpy.sum(x85 * x86)
    result[1, 1, 1] = numpy.sum(x3 * x41 * x91)
    result[1, 1, 2] = numpy.sum(x38 * x85 * x92)
    result[1, 2, 0] = numpy.sum(x66 * x76 * x86)
    result[1, 2, 1] = numpy.sum(x66 * x81 * x92)
    result[1, 2, 2] = numpy.sum(x3 * x76 * x94)
    result[1, 3, 0] = numpy.sum(x14 * x97)
    result[1, 3, 1] = numpy.sum(
        x96
        * (
            x0 * (x18 * (x64 + x79) + 4.0 * x45 * x74 + x78 + x83 + x95 * (x45 + x87))
            + x39 * x89
        )
    )
    result[1, 3, 2] = numpy.sum(x38 * x97)
    result[1, 4, 0] = numpy.sum(x14 * x51 * x85 * x90)
    result[1, 4, 1] = numpy.sum(x66 * x91)
    result[1, 4, 2] = numpy.sum(x84 * x94)
    result[1, 5, 0] = numpy.sum(x14 * x76 * x99)
    result[1, 5, 1] = numpy.sum(x81 * x99)
    result[1, 5, 2] = numpy.sum(x71 * x76 * x98)
    result[2, 0, 0] = numpy.sum(x103 * x73)
    result[2, 0, 1] = numpy.sum(x103 * x34 * x82)
    result[2, 0, 2] = numpy.sum(x107 * x68 * x82)
    result[2, 1, 0] = numpy.sum(x108 * x109 * x39)
    result[2, 1, 1] = numpy.sum(x109 * x111 * x3)
    result[2, 1, 2] = numpy.sum(22.6274169979695 * x107 * x113 * x39)
    result[2, 2, 0] = numpy.sum(x108 * x116)
    result[2, 2, 1] = numpy.sum(x113 * x116 * x34)
    result[2, 2, 2] = numpy.sum(x113 * x120)
    result[2, 3, 0] = numpy.sum(x121 * x14 * x61)
    result[2, 3, 1] = numpy.sum(x121 * x65)
    result[2, 3, 2] = numpy.sum(x107 * x61 * x98)
    result[2, 4, 0] = numpy.sum(x116 * x122 * x14)
    result[2, 4, 1] = numpy.sum(x111 * x116)
    result[2, 4, 2] = numpy.sum(x120 * x122)
    result[2, 5, 0] = numpy.sum(x125 * x14)
    result[2, 5, 1] = numpy.sum(x125 * x34)
    result[2, 5, 2] = numpy.sum(
        x124
        * (
            x0
            * (4.0 * x100 * x55 + x104 + x114 + x123 * (x117 + x55) + x18 * (x105 + x70))
            + x119 * x51
        )
    )
    return result


def diag_quadrupole3d_22(ax, da, A, bx, db, B, C):
    """Cartesian 3D (dd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 6, 6), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = -x2 - C[0]
    x5 = ax * bx * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(x1)
    x8 = 1.77245385090552 * x7
    x9 = x6 * x8
    x10 = x4**2 * x9
    x11 = x0 * x9
    x12 = x10 + x11
    x13 = x12 * x3
    x14 = 2.0 * x0
    x15 = x4 * x9
    x16 = x14 * x15
    x17 = x13 + x16
    x18 = x17 * x3
    x19 = x3**2 * x9
    x20 = 3.0 * x11
    x21 = x3 * x9
    x22 = x21 * x4
    x23 = x20 + 2.0 * x22
    x24 = x0 * (x15 + x21)
    x25 = x11 + x22
    x26 = x25 * x3
    x27 = -x2 - A[0]
    x28 = 2.0 * x27
    x29 = x17 * x27
    x30 = x0 * (x10 + x23)
    x31 = 4.0 * x11 * x4 + 2.0 * x24
    x32 = x18 + x30
    x33 = x0 * (2.0 * x13 + 2.0 * x26 + x31) + x27 * x32
    x34 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x35 = da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**3.5)
    x36 = x34 * x35
    x37 = 15.084944665313 * x36
    x38 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x39 = 0.564189583547756 * x1
    x40 = x38 * x39
    x41 = -x1 * (ax * A[1] + bx * B[1])
    x42 = -x41 - B[1]
    x43 = 26.1278905896872 * x42
    x44 = x12 * x27
    x45 = x29 + x30
    x46 = x36 * x40
    x47 = x46 * (x0 * (x13 + x25 * x28 + x31 + x44) + x27 * x45)
    x48 = -x1 * (ax * A[2] + bx * B[2])
    x49 = -x48 - B[2]
    x50 = 26.1278905896872 * x49
    x51 = x16 + x44
    x52 = x0 * (x10 + x15 * x28 + x20) + x27 * x51
    x53 = x36 * x52
    x54 = x38 * x8
    x55 = x42**2 * x54
    x56 = x0 * x54
    x57 = x55 + x56
    x58 = 15.084944665313 * x57
    x59 = 0.318309886183791 * x7
    x60 = x58 * x59
    x61 = x40 * x43
    x62 = x34 * x8
    x63 = x49**2 * x62
    x64 = x0 * x62
    x65 = x63 + x64
    x66 = 15.084944665313 * x65
    x67 = x35 * x59
    x68 = x38 * x67
    x69 = -x41 - A[1]
    x70 = 26.1278905896872 * x69
    x71 = x33 * x46
    x72 = x54 * x69
    x73 = x42 * x72 + x56
    x74 = 45.2548339959391 * x73
    x75 = x36 * x59
    x76 = x74 * x75
    x77 = 45.2548339959391 * x49
    x78 = x45 * x46
    x79 = x42 * x54
    x80 = x14 * x79 + x57 * x69
    x81 = 26.1278905896872 * x51
    x82 = x75 * x81
    x83 = x68 * x81
    x84 = -x48 - A[2]
    x85 = 26.1278905896872 * x84
    x86 = 45.2548339959391 * x84
    x87 = x62 * x84
    x88 = x49 * x87 + x64
    x89 = 45.2548339959391 * x88
    x90 = x68 * x89
    x91 = x49 * x62
    x92 = x14 * x91 + x65 * x84
    x93 = x54 * x69**2 + x56
    x94 = 15.084944665313 * x93
    x95 = x0 * (x72 + x79) + x69 * x73
    x96 = 26.1278905896872 * x17
    x97 = x75 * x96
    x98 = 2.0 * x69
    x99 = 3.0 * x56
    x100 = x55 + x99
    x101 = x0 * (x100 + x79 * x98) + x69 * x80
    x102 = x37 * x59
    x103 = 26.1278905896872 * x75
    x104 = x103 * x12
    x105 = 0.179587122125167 * x35
    x106 = x105 * x12
    x107 = x12 * x68
    x108 = x62 * x84**2 + x64
    x109 = 15.084944665313 * x108
    x110 = x68 * x96
    x111 = x0 * (x87 + x91) + x84 * x88
    x112 = 2.0 * x84
    x113 = 3.0 * x64
    x114 = x113 + x63
    x115 = 15.084944665313 * x0 * (x112 * x91 + x114) + 15.084944665313 * x84 * x92
    x116 = -x41 - C[1]
    x117 = x116**2 * x54
    x118 = x117 + x56
    x119 = x11 + x19
    x120 = x119 * x27 + x14 * x21
    x121 = x0 * (x19 + x20 + x21 * x28) + x120 * x27
    x122 = x118 * x42
    x123 = x116 * x54
    x124 = x123 * x14
    x125 = x122 + x124
    x126 = x11 + x21 * x27
    x127 = x0 * (x21 + x27 * x9) + x126 * x27
    x128 = x103 * x127
    x129 = x116 * x79
    x130 = 2.0 * x129
    x131 = x117 + x99
    x132 = x0 * (x130 + x131)
    x133 = x125 * x42
    x134 = x132 + x133
    x135 = x11 + x27**2 * x9
    x136 = x105 * x118
    x137 = x118 * x69
    x138 = x124 + x137
    x139 = x103 * x120
    x140 = x125 * x69
    x141 = x132 + x140
    x142 = 45.2548339959391 * x126
    x143 = x142 * x75
    x144 = 26.1278905896872 * x27
    x145 = x129 + x56
    x146 = x145 * x42
    x147 = x0 * (x123 + x79)
    x148 = 4.0 * x116 * x56 + 2.0 * x147
    x149 = x0 * (2.0 * x122 + 2.0 * x146 + x148) + x134 * x69
    x150 = x39 * x6
    x151 = x150 * x36
    x152 = x149 * x151
    x153 = x151 * x27
    x154 = x6 * x67
    x155 = 26.1278905896872 * x154
    x156 = x138 * x155
    x157 = x154 * x89
    x158 = x118 * x155
    x159 = x0 * (x123 * x98 + x131) + x138 * x69
    x160 = x151 * (x0 * (x122 + x137 + x145 * x98 + x148) + x141 * x69)
    x161 = 26.1278905896872 * x3
    x162 = x151 * x3
    x163 = x125 * x155
    x164 = -x48 - C[2]
    x165 = x164**2 * x62
    x166 = x165 + x64
    x167 = 15.084944665313 * x68
    x168 = 26.1278905896872 * x68
    x169 = x127 * x168
    x170 = x166 * x49
    x171 = x164 * x62
    x172 = x14 * x171
    x173 = x170 + x172
    x174 = x105 * x166
    x175 = x164 * x91
    x176 = 2.0 * x175
    x177 = x113 + x165
    x178 = x0 * (x176 + x177)
    x179 = x173 * x49
    x180 = x178 + x179
    x181 = x120 * x168
    x182 = x142 * x68
    x183 = x155 * x80
    x184 = x154 * x74
    x185 = x35 * x6
    x186 = x185 * x40
    x187 = x186 * x27
    x188 = x166 * x84
    x189 = x172 + x188
    x190 = x173 * x84
    x191 = x178 + x190
    x192 = 45.2548339959391 * x191
    x193 = x175 + x64
    x194 = x193 * x49
    x195 = x0 * (x171 + x91)
    x196 = 4.0 * x164 * x64 + 2.0 * x195
    x197 = x0 * (2.0 * x170 + 2.0 * x194 + x196) + x180 * x84
    x198 = x186 * x197
    x199 = x155 * x95
    x200 = x0 * (x112 * x171 + x177) + x189 * x84
    x201 = x185 * x200
    x202 = x186 * (x0 * (x112 * x193 + x170 + x188 + x196) + x191 * x84)

    # 108 item(s)
    result[0, 0, 0] = numpy.sum(
        x37
        * x40
        * (
            x0 * (x14 * (x19 + x23) + x18 + x28 * (x24 + x26) + 2.0 * x29 + 3.0 * x30)
            + x27 * x33
        )
    )
    result[0, 0, 1] = numpy.sum(x43 * x47)
    result[0, 0, 2] = numpy.sum(x47 * x50)
    result[0, 0, 3] = numpy.sum(x53 * x60)
    result[0, 0, 4] = numpy.sum(x49 * x53 * x61)
    result[0, 0, 5] = numpy.sum(x52 * x66 * x68)
    result[0, 1, 0] = numpy.sum(x70 * x71)
    result[0, 1, 1] = numpy.sum(x45 * x76)
    result[0, 1, 2] = numpy.sum(x69 * x77 * x78)
    result[0, 1, 3] = numpy.sum(x80 * x82)
    result[0, 1, 4] = numpy.sum(x49 * x51 * x76)
    result[0, 1, 5] = numpy.sum(x65 * x69 * x83)
    result[0, 2, 0] = numpy.sum(x71 * x85)
    result[0, 2, 1] = numpy.sum(x42 * x78 * x86)
    result[0, 2, 2] = numpy.sum(x45 * x90)
    result[0, 2, 3] = numpy.sum(x57 * x82 * x84)
    result[0, 2, 4] = numpy.sum(x42 * x51 * x90)
    result[0, 2, 5] = numpy.sum(x83 * x92)
    result[0, 3, 0] = numpy.sum(x32 * x75 * x94)
    result[0, 3, 1] = numpy.sum(x95 * x97)
    result[0, 3, 2] = numpy.sum(x49 * x93 * x97)
    result[0, 3, 3] = numpy.sum(x101 * x102 * x12)
    result[0, 3, 4] = numpy.sum(x104 * x49 * x95)
    result[0, 3, 5] = numpy.sum(x106 * x66 * x93)
    result[0, 4, 0] = numpy.sum(x32 * x46 * x70 * x84)
    result[0, 4, 1] = numpy.sum(x17 * x76 * x84)
    result[0, 4, 2] = numpy.sum(x17 * x69 * x90)
    result[0, 4, 3] = numpy.sum(x104 * x80 * x84)
    result[0, 4, 4] = numpy.sum(x106 * x73 * x89)
    result[0, 4, 5] = numpy.sum(x107 * x70 * x92)
    result[0, 5, 0] = numpy.sum(x109 * x32 * x68)
    result[0, 5, 1] = numpy.sum(x108 * x110 * x42)
    result[0, 5, 2] = numpy.sum(x110 * x111)
    result[0, 5, 3] = numpy.sum(x106 * x109 * x57)
    result[0, 5, 4] = numpy.sum(x107 * x111 * x43)
    result[0, 5, 5] = numpy.sum(x107 * x115)
    result[1, 0, 0] = numpy.sum(x102 * x118 * x121)
    result[1, 0, 1] = numpy.sum(x125 * x128)
    result[1, 0, 2] = numpy.sum(x118 * x128 * x49)
    result[1, 0, 3] = numpy.sum(x102 * x134 * x135)
    result[1, 0, 4] = numpy.sum(x103 * x125 * x135 * x49)
    result[1, 0, 5] = numpy.sum(x135 * x136 * x66)
    result[1, 1, 0] = numpy.sum(x138 * x139)
    result[1, 1, 1] = numpy.sum(x141 * x143)
    result[1, 1, 2] = numpy.sum(x138 * x143 * x49)
    result[1, 1, 3] = numpy.sum(x144 * x152)
    result[1, 1, 4] = numpy.sum(x141 * x153 * x77)
    result[1, 1, 5] = numpy.sum(x156 * x27 * x65)
    result[1, 2, 0] = numpy.sum(x118 * x139 * x84)
    result[1, 2, 1] = numpy.sum(x125 * x143 * x84)
    result[1, 2, 2] = numpy.sum(x126 * x136 * x89)
    result[1, 2, 3] = numpy.sum(x134 * x153 * x85)
    result[1, 2, 4] = numpy.sum(x125 * x157 * x27)
    result[1, 2, 5] = numpy.sum(x158 * x27 * x92)
    result[1, 3, 0] = numpy.sum(x102 * x119 * x159)
    result[1, 3, 1] = numpy.sum(x160 * x161)
    result[1, 3, 2] = numpy.sum(x159 * x162 * x50)
    result[1, 3, 3] = numpy.sum(
        x150
        * x37
        * (
            x0
            * (3.0 * x132 + x133 + x14 * (x100 + x130) + 2.0 * x140 + x98 * (x146 + x147))
            + x149 * x69
        )
    )
    result[1, 3, 4] = numpy.sum(x160 * x50)
    result[1, 3, 5] = numpy.sum(x154 * x159 * x66)
    result[1, 4, 0] = numpy.sum(x103 * x119 * x138 * x84)
    result[1, 4, 1] = numpy.sum(x141 * x162 * x86)
    result[1, 4, 2] = numpy.sum(x138 * x157 * x3)
    result[1, 4, 3] = numpy.sum(x152 * x85)
    result[1, 4, 4] = numpy.sum(x141 * x157)
    result[1, 4, 5] = numpy.sum(x156 * x92)
    result[1, 5, 0] = numpy.sum(x109 * x119 * x136)
    result[1, 5, 1] = numpy.sum(x108 * x163 * x3)
    result[1, 5, 2] = numpy.sum(x111 * x158 * x3)
    result[1, 5, 3] = numpy.sum(x109 * x134 * x154)
    result[1, 5, 4] = numpy.sum(x111 * x163)
    result[1, 5, 5] = numpy.sum(x115 * x118 * x154)
    result[2, 0, 0] = numpy.sum(x121 * x166 * x167)
    result[2, 0, 1] = numpy.sum(x166 * x169 * x42)
    result[2, 0, 2] = numpy.sum(x169 * x173)
    result[2, 0, 3] = numpy.sum(x135 * x174 * x58)
    result[2, 0, 4] = numpy.sum(x135 * x168 * x173 * x42)
    result[2, 0, 5] = numpy.sum(x135 * x167 * x180)
    result[2, 1, 0] = numpy.sum(x166 * x181 * x69)
    result[2, 1, 1] = numpy.sum(x126 * x174 * x74)
    result[2, 1, 2] = numpy.sum(x173 * x182 * x69)
    result[2, 1, 3] = numpy.sum(x166 * x183 * x27)
    result[2, 1, 4] = numpy.sum(x173 * x184 * x27)
    result[2, 1, 5] = numpy.sum(x180 * x187 * x70)
    result[2, 2, 0] = numpy.sum(x181 * x189)
    result[2, 2, 1] = numpy.sum(x182 * x189 * x42)
    result[2, 2, 2] = numpy.sum(x182 * x191)
    result[2, 2, 3] = numpy.sum(x155 * x189 * x27 * x57)
    result[2, 2, 4] = numpy.sum(x187 * x192 * x42)
    result[2, 2, 5] = numpy.sum(x144 * x198)
    result[2, 3, 0] = numpy.sum(x119 * x174 * x94)
    result[2, 3, 1] = numpy.sum(x166 * x199 * x3)
    result[2, 3, 2] = numpy.sum(x155 * x173 * x3 * x93)
    result[2, 3, 3] = numpy.sum(15.084944665313 * x101 * x154 * x166)
    result[2, 3, 4] = numpy.sum(x173 * x199)
    result[2, 3, 5] = numpy.sum(x154 * x180 * x94)
    result[2, 4, 0] = numpy.sum(x119 * x168 * x189 * x69)
    result[2, 4, 1] = numpy.sum(x184 * x189 * x3)
    result[2, 4, 2] = numpy.sum(x186 * x192 * x3 * x69)
    result[2, 4, 3] = numpy.sum(x183 * x189)
    result[2, 4, 4] = numpy.sum(x184 * x191)
    result[2, 4, 5] = numpy.sum(x198 * x70)
    result[2, 5, 0] = numpy.sum(x119 * x167 * x200)
    result[2, 5, 1] = numpy.sum(x201 * x3 * x61)
    result[2, 5, 2] = numpy.sum(x161 * x202)
    result[2, 5, 3] = numpy.sum(x201 * x60)
    result[2, 5, 4] = numpy.sum(x202 * x43)
    result[2, 5, 5] = numpy.sum(
        15.084944665313
        * x186
        * (
            x0
            * (
                x112 * (x194 + x195)
                + x14 * (x114 + x176)
                + 3.0 * x178
                + x179
                + 2.0 * x190
            )
            + x197 * x84
        )
    )
    return result


def diag_quadrupole3d_23(ax, da, A, bx, db, B, C):
    """Cartesian 3D (df) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 6, 10), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = -x2 - C[0]
    x5 = ax * bx * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(x1)
    x8 = 1.77245385090552 * x7
    x9 = x6 * x8
    x10 = x4**2 * x9
    x11 = x0 * x9
    x12 = 3.0 * x11
    x13 = x3 * x9
    x14 = x13 * x4
    x15 = x12 + 2.0 * x14
    x16 = x0 * (x10 + x15)
    x17 = x10 + x11
    x18 = x17 * x3
    x19 = 2.0 * x0
    x20 = x4 * x9
    x21 = x19 * x20
    x22 = x18 + x21
    x23 = x22 * x3
    x24 = x16 + x23
    x25 = x24 * x3
    x26 = x0 * (x13 + x20)
    x27 = x11 + x14
    x28 = x27 * x3
    x29 = x3**2 * x9
    x30 = x11 + x29
    x31 = x3 * x30
    x32 = x13 * x19
    x33 = x31 + x32
    x34 = x0 * (x15 + x29)
    x35 = x26 + x28
    x36 = x3 * x35
    x37 = -x2 - A[0]
    x38 = 2.0 * x37
    x39 = x24 * x37
    x40 = 4.0 * x11 * x4 + 2.0 * x26
    x41 = x0 * (2.0 * x18 + 2.0 * x28 + x40)
    x42 = 3.0 * x16 + 2.0 * x34
    x43 = x25 + x41
    x44 = x0 * (3.0 * x23 + 2.0 * x36 + x42) + x37 * x43
    x45 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x46 = da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**4.5)
    x47 = x45 * x46
    x48 = 13.4923846833851 * x47
    x49 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x50 = 0.564189583547756 * x1
    x51 = x49 * x50
    x52 = -x1 * (ax * A[1] + bx * B[1])
    x53 = -x52 - B[1]
    x54 = 30.169889330626 * x53
    x55 = x22 * x37
    x56 = x39 + x41
    x57 = x47 * x51
    x58 = x57 * (x0 * (x23 + x35 * x38 + x42 + 2.0 * x55) + x37 * x56)
    x59 = -x1 * (ax * A[2] + bx * B[2])
    x60 = -x59 - B[2]
    x61 = 30.169889330626 * x60
    x62 = x17 * x37
    x63 = x16 + x55
    x64 = x0 * (x18 + x27 * x38 + x40 + x62) + x37 * x63
    x65 = x47 * x64
    x66 = x49 * x8
    x67 = x53**2 * x66
    x68 = x0 * x66
    x69 = x67 + x68
    x70 = 30.169889330626 * x69
    x71 = 0.318309886183791 * x7
    x72 = x70 * x71
    x73 = 52.2557811793745 * x60
    x74 = x45 * x8
    x75 = x60**2 * x74
    x76 = x0 * x74
    x77 = x75 + x76
    x78 = 30.169889330626 * x71
    x79 = x46 * x49
    x80 = x78 * x79
    x81 = x77 * x80
    x82 = x53 * x69
    x83 = x53 * x66
    x84 = x19 * x83
    x85 = x82 + x84
    x86 = 13.4923846833851 * x85
    x87 = x21 + x62
    x88 = x0 * (x10 + x12 + x20 * x38) + x37 * x87
    x89 = x71 * x88
    x90 = x47 * x60
    x91 = x60 * x77
    x92 = x60 * x74
    x93 = x19 * x92
    x94 = x91 + x93
    x95 = 13.4923846833851 * x94
    x96 = -x52 - A[1]
    x97 = 23.3694957868871 * x96
    x98 = x44 * x57
    x99 = x66 * x96
    x100 = x53 * x99 + x68
    x101 = 52.2557811793745 * x100
    x102 = x47 * x71
    x103 = x101 * x102
    x104 = 52.2557811793745 * x56
    x105 = x104 * x57
    x106 = 52.2557811793745 * x63
    x107 = x69 * x96
    x108 = x107 + x84
    x109 = x102 * x108
    x110 = 90.5096679918781 * x100
    x111 = x102 * x60
    x112 = x71 * x79
    x113 = x106 * x112
    x114 = 3.0 * x68
    x115 = x0 * (x114 + 3.0 * x67) + x85 * x96
    x116 = 23.3694957868871 * x87
    x117 = x102 * x116
    x118 = 52.2557811793745 * x87
    x119 = 0.179587122125167 * x46
    x120 = x119 * x77
    x121 = x112 * x116
    x122 = -x59 - A[2]
    x123 = 23.3694957868871 * x122
    x124 = x122 * x74
    x125 = x124 * x60 + x76
    x126 = x112 * x125
    x127 = x102 * x122
    x128 = 90.5096679918781 * x53
    x129 = x122 * x77
    x130 = x129 + x93
    x131 = x119 * x125
    x132 = 3.0 * x76
    x133 = x0 * (x132 + 3.0 * x75) + x122 * x94
    x134 = x66 * x96**2 + x68
    x135 = x48 * x71
    x136 = x0 * (x83 + x99) + x100 * x96
    x137 = x47 * x78
    x138 = x137 * x24
    x139 = 2.0 * x96
    x140 = x114 + x67
    x141 = x0 * (x139 * x83 + x140) + x108 * x96
    x142 = x137 * x141
    x143 = 52.2557811793745 * x22
    x144 = 30.169889330626 * x120
    x145 = x0 * (3.0 * x107 + 8.0 * x53 * x68 + x82) + x115 * x96
    x146 = x119 * x17
    x147 = 52.2557811793745 * x96
    x148 = 90.5096679918781 * x131
    x149 = 52.2557811793745 * x112
    x150 = x149 * x22
    x151 = 52.2557811793745 * x131
    x152 = x122**2 * x74 + x76
    x153 = 13.4923846833851 * x112
    x154 = x24 * x80
    x155 = x0 * (x124 + x92) + x122 * x125
    x156 = x119 * x152
    x157 = 2.0 * x122
    x158 = x132 + x75
    x159 = x0 * (x157 * x92 + x158) + x122 * x130
    x160 = x159 * x80
    x161 = x0 * (3.0 * x129 + 8.0 * x60 * x76 + x91) + x122 * x133
    x162 = -x52 - C[1]
    x163 = x162**2 * x66
    x164 = x163 + x68
    x165 = x30 * x37
    x166 = x0 * (x12 + 3.0 * x29) + x33 * x37
    x167 = x0 * (8.0 * x11 * x3 + 3.0 * x165 + x31) + x166 * x37
    x168 = x164 * x53
    x169 = x162 * x66
    x170 = x169 * x19
    x171 = x168 + x170
    x172 = x165 + x32
    x173 = x0 * (x12 + x13 * x38 + x29) + x172 * x37
    x174 = x137 * x173
    x175 = x11 + x13 * x37
    x176 = x0 * (x13 + x37 * x9) + x175 * x37
    x177 = x162 * x83
    x178 = 2.0 * x177
    x179 = x114 + x163
    x180 = x0 * (x178 + x179)
    x181 = x171 * x53
    x182 = x180 + x181
    x183 = x137 * x182
    x184 = 52.2557811793745 * x171
    x185 = x177 + x68
    x186 = x185 * x53
    x187 = x0 * (x169 + x83)
    x188 = 4.0 * x162 * x68 + 2.0 * x187
    x189 = x0 * (2.0 * x168 + 2.0 * x186 + x188)
    x190 = x182 * x53
    x191 = x189 + x190
    x192 = x11 + x37**2 * x9
    x193 = x119 * x164
    x194 = x164 * x96
    x195 = x170 + x194
    x196 = 23.3694957868871 * x195
    x197 = x102 * x166
    x198 = x171 * x96
    x199 = x180 + x198
    x200 = 52.2557811793745 * x199
    x201 = x102 * x172
    x202 = x182 * x96
    x203 = x189 + x202
    x204 = 52.2557811793745 * x175
    x205 = x102 * x204
    x206 = 90.5096679918781 * x199
    x207 = 23.3694957868871 * x37
    x208 = x186 + x187
    x209 = x208 * x53
    x210 = x0 * (x140 + x178)
    x211 = 3.0 * x180 + 2.0 * x210
    x212 = x0 * (3.0 * x181 + 2.0 * x209 + x211) + x191 * x96
    x213 = x50 * x6
    x214 = x213 * x47
    x215 = x212 * x214
    x216 = 52.2557811793745 * x203
    x217 = x213 * x90
    x218 = x46 * x6
    x219 = x218 * x71
    x220 = x200 * x219
    x221 = x196 * x219
    x222 = 52.2557811793745 * x219
    x223 = x125 * x222
    x224 = x130 * x222
    x225 = 13.4923846833851 * x33
    x226 = x0 * (x139 * x169 + x179) + x195 * x96
    x227 = x102 * x226
    x228 = x0 * (x139 * x185 + x168 + x188 + x194) + x199 * x96
    x229 = 30.169889330626 * x30
    x230 = x214 * (x0 * (x139 * x208 + x181 + 2.0 * x198 + x211) + x203 * x96)
    x231 = 30.169889330626 * x3
    x232 = 52.2557811793745 * x3
    x233 = x218 * x78
    x234 = x233 * x77
    x235 = x219 * x3
    x236 = x182 * x233
    x237 = x222 * x3
    x238 = x159 * x233
    x239 = 13.4923846833851 * x219
    x240 = -x59 - C[2]
    x241 = x240**2 * x74
    x242 = x241 + x76
    x243 = x173 * x80
    x244 = x242 * x60
    x245 = x240 * x74
    x246 = x19 * x245
    x247 = x244 + x246
    x248 = x119 * x242
    x249 = x240 * x92
    x250 = 2.0 * x249
    x251 = x132 + x241
    x252 = x0 * (x250 + x251)
    x253 = x247 * x60
    x254 = x252 + x253
    x255 = x254 * x80
    x256 = x119 * x247
    x257 = x249 + x76
    x258 = x257 * x60
    x259 = x0 * (x245 + x92)
    x260 = 4.0 * x240 * x76 + 2.0 * x259
    x261 = x0 * (2.0 * x244 + 2.0 * x258 + x260)
    x262 = x254 * x60
    x263 = x261 + x262
    x264 = x112 * x166
    x265 = x149 * x172
    x266 = x112 * x204
    x267 = x115 * x219
    x268 = x108 * x222
    x269 = x101 * x219
    x270 = x218 * x51
    x271 = x270 * x37
    x272 = x122 * x242
    x273 = x246 + x272
    x274 = 23.3694957868871 * x273
    x275 = x122 * x247
    x276 = x252 + x275
    x277 = x119 * x273
    x278 = x122 * x254
    x279 = x261 + x278
    x280 = x258 + x259
    x281 = x280 * x60
    x282 = x0 * (x158 + x250)
    x283 = 3.0 * x252 + 2.0 * x282
    x284 = x0 * (3.0 * x253 + 2.0 * x281 + x283) + x122 * x263
    x285 = x270 * x284
    x286 = x141 * x233
    x287 = x233 * x254
    x288 = x0 * (x157 * x245 + x251) + x122 * x273
    x289 = x112 * x288
    x290 = x0 * (x157 * x257 + x244 + x260 + x272) + x122 * x276
    x291 = x218 * x72
    x292 = x270 * (x0 * (x157 * x280 + x253 + 2.0 * x275 + x283) + x122 * x279)

    # 180 item(s)
    result[0, 0, 0] = numpy.sum(
        x48
        * x51
        * (
            x0
            * (
                x19 * (3.0 * x26 + 3.0 * x28 + x33)
                + x25
                + x38 * (x34 + x36)
                + 3.0 * x39
                + 4.0 * x41
            )
            + x37 * x44
        )
    )
    result[0, 0, 1] = numpy.sum(x54 * x58)
    result[0, 0, 2] = numpy.sum(x58 * x61)
    result[0, 0, 3] = numpy.sum(x65 * x72)
    result[0, 0, 4] = numpy.sum(x51 * x53 * x65 * x73)
    result[0, 0, 5] = numpy.sum(x64 * x81)
    result[0, 0, 6] = numpy.sum(x47 * x86 * x89)
    result[0, 0, 7] = numpy.sum(x72 * x88 * x90)
    result[0, 0, 8] = numpy.sum(x53 * x81 * x88)
    result[0, 0, 9] = numpy.sum(x79 * x89 * x95)
    result[0, 1, 0] = numpy.sum(x97 * x98)
    result[0, 1, 1] = numpy.sum(x103 * x56)
    result[0, 1, 2] = numpy.sum(x105 * x60 * x96)
    result[0, 1, 3] = numpy.sum(x106 * x109)
    result[0, 1, 4] = numpy.sum(x110 * x111 * x63)
    result[0, 1, 5] = numpy.sum(x113 * x77 * x96)
    result[0, 1, 6] = numpy.sum(x115 * x117)
    result[0, 1, 7] = numpy.sum(x109 * x118 * x60)
    result[0, 1, 8] = numpy.sum(x101 * x120 * x87)
    result[0, 1, 9] = numpy.sum(x121 * x94 * x96)
    result[0, 2, 0] = numpy.sum(x123 * x98)
    result[0, 2, 1] = numpy.sum(x105 * x122 * x53)
    result[0, 2, 2] = numpy.sum(x104 * x126)
    result[0, 2, 3] = numpy.sum(x106 * x127 * x69)
    result[0, 2, 4] = numpy.sum(x126 * x128 * x63)
    result[0, 2, 5] = numpy.sum(x113 * x130)
    result[0, 2, 6] = numpy.sum(x117 * x122 * x85)
    result[0, 2, 7] = numpy.sum(x118 * x131 * x69)
    result[0, 2, 8] = numpy.sum(x112 * x118 * x130 * x53)
    result[0, 2, 9] = numpy.sum(x121 * x133)
    result[0, 3, 0] = numpy.sum(x134 * x135 * x43)
    result[0, 3, 1] = numpy.sum(x136 * x138)
    result[0, 3, 2] = numpy.sum(x134 * x138 * x60)
    result[0, 3, 3] = numpy.sum(x142 * x22)
    result[0, 3, 4] = numpy.sum(x111 * x136 * x143)
    result[0, 3, 5] = numpy.sum(x134 * x144 * x22)
    result[0, 3, 6] = numpy.sum(x135 * x145 * x17)
    result[0, 3, 7] = numpy.sum(x142 * x17 * x60)
    result[0, 3, 8] = numpy.sum(x136 * x144 * x17)
    result[0, 3, 9] = numpy.sum(x134 * x146 * x95)
    result[0, 4, 0] = numpy.sum(x122 * x43 * x57 * x97)
    result[0, 4, 1] = numpy.sum(x103 * x122 * x24)
    result[0, 4, 2] = numpy.sum(x126 * x147 * x24)
    result[0, 4, 3] = numpy.sum(x109 * x122 * x143)
    result[0, 4, 4] = numpy.sum(x100 * x148 * x22)
    result[0, 4, 5] = numpy.sum(x130 * x150 * x96)
    result[0, 4, 6] = numpy.sum(x102 * x115 * x123 * x17)
    result[0, 4, 7] = numpy.sum(x108 * x151 * x17)
    result[0, 4, 8] = numpy.sum(x101 * x130 * x146)
    result[0, 4, 9] = numpy.sum(x112 * x133 * x17 * x97)
    result[0, 5, 0] = numpy.sum(x152 * x153 * x43)
    result[0, 5, 1] = numpy.sum(x152 * x154 * x53)
    result[0, 5, 2] = numpy.sum(x154 * x155)
    result[0, 5, 3] = numpy.sum(x156 * x22 * x70)
    result[0, 5, 4] = numpy.sum(x150 * x155 * x53)
    result[0, 5, 5] = numpy.sum(x160 * x22)
    result[0, 5, 6] = numpy.sum(x146 * x152 * x86)
    result[0, 5, 7] = numpy.sum(x146 * x155 * x70)
    result[0, 5, 8] = numpy.sum(x160 * x17 * x53)
    result[0, 5, 9] = numpy.sum(x153 * x161 * x17)
    result[1, 0, 0] = numpy.sum(x135 * x164 * x167)
    result[1, 0, 1] = numpy.sum(x171 * x174)
    result[1, 0, 2] = numpy.sum(x164 * x174 * x60)
    result[1, 0, 3] = numpy.sum(x176 * x183)
    result[1, 0, 4] = numpy.sum(x111 * x176 * x184)
    result[1, 0, 5] = numpy.sum(x144 * x164 * x176)
    result[1, 0, 6] = numpy.sum(x135 * x191 * x192)
    result[1, 0, 7] = numpy.sum(x183 * x192 * x60)
    result[1, 0, 8] = numpy.sum(x144 * x171 * x192)
    result[1, 0, 9] = numpy.sum(x192 * x193 * x95)
    result[1, 1, 0] = numpy.sum(x196 * x197)
    result[1, 1, 1] = numpy.sum(x200 * x201)
    result[1, 1, 2] = numpy.sum(x195 * x201 * x73)
    result[1, 1, 3] = numpy.sum(x203 * x205)
    result[1, 1, 4] = numpy.sum(x111 * x175 * x206)
    result[1, 1, 5] = numpy.sum(x120 * x195 * x204)
    result[1, 1, 6] = numpy.sum(x207 * x215)
    result[1, 1, 7] = numpy.sum(x216 * x217 * x37)
    result[1, 1, 8] = numpy.sum(x220 * x37 * x77)
    result[1, 1, 9] = numpy.sum(x221 * x37 * x94)
    result[1, 2, 0] = numpy.sum(x123 * x164 * x197)
    result[1, 2, 1] = numpy.sum(x122 * x184 * x201)
    result[1, 2, 2] = numpy.sum(x151 * x164 * x172)
    result[1, 2, 3] = numpy.sum(x122 * x182 * x205)
    result[1, 2, 4] = numpy.sum(x148 * x171 * x175)
    result[1, 2, 5] = numpy.sum(x130 * x193 * x204)
    result[1, 2, 6] = numpy.sum(x123 * x191 * x214 * x37)
    result[1, 2, 7] = numpy.sum(x182 * x223 * x37)
    result[1, 2, 8] = numpy.sum(x171 * x224 * x37)
    result[1, 2, 9] = numpy.sum(x133 * x164 * x207 * x219)
    result[1, 3, 0] = numpy.sum(x225 * x227)
    result[1, 3, 1] = numpy.sum(x102 * x228 * x229)
    result[1, 3, 2] = numpy.sum(x227 * x229 * x60)
    result[1, 3, 3] = numpy.sum(x230 * x231)
    result[1, 3, 4] = numpy.sum(x217 * x228 * x232)
    result[1, 3, 5] = numpy.sum(x226 * x234 * x3)
    result[1, 3, 6] = numpy.sum(
        x213
        * x48
        * (
            x0
            * (
                x139 * (x209 + x210)
                + 4.0 * x189
                + x19 * (3.0 * x186 + 3.0 * x187 + x85)
                + x190
                + 3.0 * x202
            )
            + x212 * x96
        )
    )
    result[1, 3, 7] = numpy.sum(x230 * x61)
    result[1, 3, 8] = numpy.sum(x228 * x234)
    result[1, 3, 9] = numpy.sum(x219 * x226 * x95)
    result[1, 4, 0] = numpy.sum(x127 * x196 * x33)
    result[1, 4, 1] = numpy.sum(x127 * x200 * x30)
    result[1, 4, 2] = numpy.sum(x151 * x195 * x30)
    result[1, 4, 3] = numpy.sum(x122 * x214 * x216 * x3)
    result[1, 4, 4] = numpy.sum(x125 * x206 * x235)
    result[1, 4, 5] = numpy.sum(x195 * x224 * x3)
    result[1, 4, 6] = numpy.sum(x123 * x215)
    result[1, 4, 7] = numpy.sum(x203 * x223)
    result[1, 4, 8] = numpy.sum(x130 * x220)
    result[1, 4, 9] = numpy.sum(x133 * x221)
    result[1, 5, 0] = numpy.sum(x156 * x164 * x225)
    result[1, 5, 1] = numpy.sum(x156 * x171 * x229)
    result[1, 5, 2] = numpy.sum(x155 * x193 * x229)
    result[1, 5, 3] = numpy.sum(x152 * x236 * x3)
    result[1, 5, 4] = numpy.sum(x155 * x171 * x237)
    result[1, 5, 5] = numpy.sum(x164 * x238 * x3)
    result[1, 5, 6] = numpy.sum(x152 * x191 * x239)
    result[1, 5, 7] = numpy.sum(x155 * x236)
    result[1, 5, 8] = numpy.sum(x171 * x238)
    result[1, 5, 9] = numpy.sum(x161 * x164 * x239)
    result[2, 0, 0] = numpy.sum(x153 * x167 * x242)
    result[2, 0, 1] = numpy.sum(x242 * x243 * x53)
    result[2, 0, 2] = numpy.sum(x243 * x247)
    result[2, 0, 3] = numpy.sum(x176 * x248 * x70)
    result[2, 0, 4] = numpy.sum(x149 * x176 * x247 * x53)
    result[2, 0, 5] = numpy.sum(x176 * x255)
    result[2, 0, 6] = numpy.sum(x192 * x248 * x86)
    result[2, 0, 7] = numpy.sum(x192 * x256 * x70)
    result[2, 0, 8] = numpy.sum(x192 * x255 * x53)
    result[2, 0, 9] = numpy.sum(x153 * x192 * x263)
    result[2, 1, 0] = numpy.sum(x242 * x264 * x97)
    result[2, 1, 1] = numpy.sum(x101 * x172 * x248)
    result[2, 1, 2] = numpy.sum(x247 * x265 * x96)
    result[2, 1, 3] = numpy.sum(x108 * x204 * x248)
    result[2, 1, 4] = numpy.sum(x110 * x175 * x256)
    result[2, 1, 5] = numpy.sum(x254 * x266 * x96)
    result[2, 1, 6] = numpy.sum(x207 * x242 * x267)
    result[2, 1, 7] = numpy.sum(x247 * x268 * x37)
    result[2, 1, 8] = numpy.sum(x254 * x269 * x37)
    result[2, 1, 9] = numpy.sum(x263 * x271 * x97)
    result[2, 2, 0] = numpy.sum(x264 * x274)
    result[2, 2, 1] = numpy.sum(x265 * x273 * x53)
    result[2, 2, 2] = numpy.sum(x265 * x276)
    result[2, 2, 3] = numpy.sum(x204 * x277 * x69)
    result[2, 2, 4] = numpy.sum(x112 * x128 * x175 * x276)
    result[2, 2, 5] = numpy.sum(x266 * x279)
    result[2, 2, 6] = numpy.sum(x219 * x274 * x37 * x85)
    result[2, 2, 7] = numpy.sum(x222 * x276 * x37 * x69)
    result[2, 2, 8] = numpy.sum(52.2557811793745 * x271 * x279 * x53)
    result[2, 2, 9] = numpy.sum(x207 * x285)
    result[2, 3, 0] = numpy.sum(x134 * x225 * x248)
    result[2, 3, 1] = numpy.sum(x136 * x229 * x248)
    result[2, 3, 2] = numpy.sum(x134 * x229 * x256)
    result[2, 3, 3] = numpy.sum(x242 * x286 * x3)
    result[2, 3, 4] = numpy.sum(x136 * x237 * x247)
    result[2, 3, 5] = numpy.sum(x134 * x287 * x3)
    result[2, 3, 6] = numpy.sum(x145 * x239 * x242)
    result[2, 3, 7] = numpy.sum(x247 * x286)
    result[2, 3, 8] = numpy.sum(x136 * x287)
    result[2, 3, 9] = numpy.sum(x134 * x239 * x263)
    result[2, 4, 0] = numpy.sum(x112 * x274 * x33 * x96)
    result[2, 4, 1] = numpy.sum(x101 * x277 * x30)
    result[2, 4, 2] = numpy.sum(x149 * x276 * x30 * x96)
    result[2, 4, 3] = numpy.sum(x268 * x273 * x3)
    result[2, 4, 4] = numpy.sum(x110 * x235 * x276)
    result[2, 4, 5] = numpy.sum(x147 * x270 * x279 * x3)
    result[2, 4, 6] = numpy.sum(x267 * x274)
    result[2, 4, 7] = numpy.sum(x268 * x276)
    result[2, 4, 8] = numpy.sum(x269 * x279)
    result[2, 4, 9] = numpy.sum(x285 * x97)
    result[2, 5, 0] = numpy.sum(x225 * x289)
    result[2, 5, 1] = numpy.sum(x229 * x289 * x53)
    result[2, 5, 2] = numpy.sum(x112 * x229 * x290)
    result[2, 5, 3] = numpy.sum(x288 * x291 * x3)
    result[2, 5, 4] = numpy.sum(x232 * x270 * x290 * x53)
    result[2, 5, 5] = numpy.sum(x231 * x292)
    result[2, 5, 6] = numpy.sum(x219 * x288 * x86)
    result[2, 5, 7] = numpy.sum(x290 * x291)
    result[2, 5, 8] = numpy.sum(x292 * x54)
    result[2, 5, 9] = numpy.sum(
        13.4923846833851
        * x270
        * (
            x0
            * (
                x157 * (x281 + x282)
                + x19 * (3.0 * x258 + 3.0 * x259 + x94)
                + 4.0 * x261
                + x262
                + 3.0 * x278
            )
            + x122 * x284
        )
    )
    return result


def diag_quadrupole3d_24(ax, da, A, bx, db, B, C):
    """Cartesian 3D (dg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 6, 15), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = -x2 - C[0]
    x5 = ax * bx * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(x1)
    x8 = 1.77245385090552 * x7
    x9 = x6 * x8
    x10 = x4**2 * x9
    x11 = x0 * x9
    x12 = x10 + x11
    x13 = x12 * x3
    x14 = x3 * x6
    x15 = x14 * x8
    x16 = x15 * x4
    x17 = x11 + x16
    x18 = x17 * x3
    x19 = x4 * x9
    x20 = x0 * (x15 + x19)
    x21 = 4.0 * x0 * x19 + 2.0 * x20
    x22 = x0 * (2.0 * x13 + 2.0 * x18 + x21)
    x23 = 3.0 * x11
    x24 = 2.0 * x16 + x23
    x25 = x0 * (x10 + x24)
    x26 = 2.0 * x0
    x27 = x19 * x26
    x28 = x13 + x27
    x29 = x28 * x3
    x30 = x25 + x29
    x31 = x3 * x30
    x32 = x22 + x31
    x33 = x3 * x32
    x34 = x3**2 * x9
    x35 = x0 * (x24 + x34)
    x36 = x18 + x20
    x37 = x3 * x36
    x38 = x0 * (x23 + 3.0 * x34)
    x39 = x11 + x34
    x40 = x3 * x39
    x41 = x15 * x26
    x42 = x40 + x41
    x43 = x3 * x42
    x44 = x38 + x43
    x45 = x0 * (3.0 * x18 + 3.0 * x20 + x42)
    x46 = x35 + x37
    x47 = x3 * x46
    x48 = -x2 - A[0]
    x49 = 2.0 * x48
    x50 = x32 * x48
    x51 = 3.0 * x25 + 2.0 * x35
    x52 = x0 * (3.0 * x29 + 2.0 * x37 + x51)
    x53 = 4.0 * x22 + 2.0 * x45
    x54 = x33 + x52
    x55 = x0 * (4.0 * x31 + 2.0 * x47 + x53) + x48 * x54
    x56 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x57 = da * db * numpy.sqrt(ax**3.5) * numpy.sqrt(bx**5.5)
    x58 = x56 * x57
    x59 = 10.1992841329868 * x58
    x60 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x61 = 0.564189583547756 * x1
    x62 = x60 * x61
    x63 = -x1 * (ax * A[1] + bx * B[1])
    x64 = -x63 - B[1]
    x65 = 26.9847693667702 * x64
    x66 = x30 * x48
    x67 = x50 + x52
    x68 = x58 * x62
    x69 = x68 * (x0 * (x31 + x46 * x49 + x53 + 3.0 * x66) + x48 * x67)
    x70 = -x1 * (ax * A[2] + bx * B[2])
    x71 = -x70 - B[2]
    x72 = 26.9847693667702 * x71
    x73 = x28 * x48
    x74 = x22 + x66
    x75 = x0 * (x29 + x36 * x49 + x51 + 2.0 * x73) + x48 * x74
    x76 = x58 * x75
    x77 = x60 * x8
    x78 = x64**2 * x77
    x79 = x0 * x77
    x80 = x78 + x79
    x81 = 34.8371874529163 * x80
    x82 = 0.318309886183791 * x7
    x83 = x81 * x82
    x84 = 60.3397786612521 * x64
    x85 = x56 * x8
    x86 = x71**2 * x85
    x87 = x0 * x85
    x88 = x86 + x87
    x89 = 34.8371874529163 * x82
    x90 = x57 * x60
    x91 = x89 * x90
    x92 = x64 * x80
    x93 = x64 * x77
    x94 = x26 * x93
    x95 = x92 + x94
    x96 = 26.9847693667702 * x95
    x97 = x12 * x48
    x98 = x25 + x73
    x99 = x82 * (x0 * (x13 + x17 * x49 + x21 + x97) + x48 * x98)
    x100 = x58 * x99
    x101 = 60.3397786612521 * x80
    x102 = x71 * x88
    x103 = x71 * x85
    x104 = x103 * x26
    x105 = x102 + x104
    x106 = 26.9847693667702 * x90
    x107 = x105 * x106
    x108 = 3.0 * x79
    x109 = x0 * (x108 + 3.0 * x78)
    x110 = x64 * x95
    x111 = x109 + x110
    x112 = 10.1992841329868 * x111
    x113 = x27 + x97
    x114 = x0 * (x10 + x19 * x49 + x23) + x113 * x48
    x115 = x114 * x82
    x116 = x115 * x58
    x117 = 0.179587122125167 * x57
    x118 = x117 * x88
    x119 = 3.0 * x87
    x120 = x0 * (x119 + 3.0 * x86)
    x121 = x105 * x71
    x122 = x120 + x121
    x123 = 10.1992841329868 * x122
    x124 = -x63 - A[1]
    x125 = 17.6656783191643 * x124
    x126 = x55 * x68
    x127 = x124 * x77
    x128 = x127 * x64 + x79
    x129 = 46.7389915737742 * x128
    x130 = x58 * x82
    x131 = x129 * x130
    x132 = 46.7389915737742 * x67
    x133 = x132 * x68
    x134 = x124 * x80
    x135 = x134 + x94
    x136 = 60.3397786612521 * x135
    x137 = x130 * x74
    x138 = 104.511562358749 * x128
    x139 = x82 * x90
    x140 = x139 * x74
    x141 = 60.3397786612521 * x88
    x142 = 46.7389915737742 * x98
    x143 = x124 * x95
    x144 = x109 + x143
    x145 = x130 * x144
    x146 = 104.511562358749 * x98
    x147 = x130 * x71
    x148 = x139 * x142
    x149 = 8.0 * x64 * x79
    x150 = x0 * (x149 + 4.0 * x92) + x111 * x124
    x151 = 17.6656783191643 * x113
    x152 = x130 * x151
    x153 = 46.7389915737742 * x145
    x154 = 60.3397786612521 * x113
    x155 = x105 * x117
    x156 = x139 * x151
    x157 = -x70 - A[2]
    x158 = 17.6656783191643 * x157
    x159 = x157 * x85
    x160 = x159 * x71 + x87
    x161 = x139 * x160
    x162 = 60.3397786612521 * x157
    x163 = 104.511562358749 * x64
    x164 = x157 * x88
    x165 = x104 + x164
    x166 = 60.3397786612521 * x165
    x167 = x130 * x157
    x168 = x117 * x160
    x169 = x139 * x64
    x170 = x105 * x157
    x171 = x120 + x170
    x172 = 46.7389915737742 * x168
    x173 = x117 * x165
    x174 = 46.7389915737742 * x139
    x175 = x171 * x174
    x176 = 8.0 * x71 * x87
    x177 = x0 * (4.0 * x102 + x176) + x122 * x157
    x178 = x124**2 * x77 + x79
    x179 = x59 * x82
    x180 = x0 * (x127 + x93) + x124 * x128
    x181 = 26.9847693667702 * x130
    x182 = x181 * x32
    x183 = 2.0 * x124
    x184 = x108 + x78
    x185 = x0 * (x183 * x93 + x184) + x124 * x135
    x186 = x58 * x89
    x187 = 60.3397786612521 * x180
    x188 = 34.8371874529163 * x118
    x189 = x0 * (3.0 * x134 + x149 + x92) + x124 * x144
    x190 = x181 * x189
    x191 = 60.3397786612521 * x28
    x192 = 26.9847693667702 * x155
    x193 = x0 * (5.0 * x109 + x110 + 4.0 * x143) + x124 * x150
    x194 = x117 * x12
    x195 = 46.7389915737742 * x124
    x196 = 104.511562358749 * x168
    x197 = x139 * x30
    x198 = 104.511562358749 * x173
    x199 = 60.3397786612521 * x173
    x200 = x157**2 * x85 + x87
    x201 = 10.1992841329868 * x139
    x202 = x106 * x82
    x203 = x202 * x32
    x204 = x0 * (x103 + x159) + x157 * x160
    x205 = x117 * x200
    x206 = 2.0 * x157
    x207 = x119 + x86
    x208 = x0 * (x103 * x206 + x207) + x157 * x165
    x209 = x117 * x204
    x210 = x0 * (x102 + 3.0 * x164 + x176) + x157 * x171
    x211 = x202 * x210
    x212 = x0 * (5.0 * x120 + x121 + 4.0 * x170) + x157 * x177
    x213 = -x63 - C[1]
    x214 = x213**2 * x77
    x215 = x214 + x79
    x216 = x42 * x48
    x217 = 8.0 * x0 * x15
    x218 = x0 * (x217 + 4.0 * x40) + x44 * x48
    x219 = x0 * (4.0 * x216 + 5.0 * x38 + x43) + x218 * x48
    x220 = x215 * x64
    x221 = x213 * x77
    x222 = x221 * x26
    x223 = x220 + x222
    x224 = x39 * x48
    x225 = x216 + x38
    x226 = x0 * (x217 + 3.0 * x224 + x40) + x225 * x48
    x227 = x181 * x226
    x228 = x213 * x93
    x229 = 2.0 * x228
    x230 = x108 + x214
    x231 = x0 * (x229 + x230)
    x232 = x223 * x64
    x233 = x231 + x232
    x234 = x224 + x41
    x235 = x0 * (x15 * x49 + x23 + x34) + x234 * x48
    x236 = 60.3397786612521 * x223
    x237 = x11 + x15 * x48
    x238 = x0 * (x15 + x48 * x9) + x237 * x48
    x239 = x228 + x79
    x240 = x239 * x64
    x241 = x0 * (x221 + x93)
    x242 = 4.0 * x213 * x79 + 2.0 * x241
    x243 = x0 * (2.0 * x220 + 2.0 * x240 + x242)
    x244 = x233 * x64
    x245 = x243 + x244
    x246 = x181 * x245
    x247 = 60.3397786612521 * x233
    x248 = x240 + x241
    x249 = x248 * x64
    x250 = x0 * (x184 + x229)
    x251 = 3.0 * x231 + 2.0 * x250
    x252 = x0 * (3.0 * x232 + 2.0 * x249 + x251)
    x253 = x245 * x64
    x254 = x252 + x253
    x255 = x11 + x48**2 * x9
    x256 = x117 * x215
    x257 = x124 * x215
    x258 = x222 + x257
    x259 = 17.6656783191643 * x258
    x260 = x130 * x218
    x261 = x124 * x223
    x262 = x231 + x261
    x263 = 46.7389915737742 * x262
    x264 = x130 * x225
    x265 = 46.7389915737742 * x264
    x266 = 60.3397786612521 * x234
    x267 = x124 * x233
    x268 = x243 + x267
    x269 = x130 * x268
    x270 = 104.511562358749 * x262
    x271 = x124 * x245
    x272 = x252 + x271
    x273 = 46.7389915737742 * x237
    x274 = x130 * x273
    x275 = 104.511562358749 * x237
    x276 = 17.6656783191643 * x48
    x277 = x249 + x250
    x278 = x277 * x64
    x279 = x0 * (3.0 * x240 + 3.0 * x241 + x95)
    x280 = 4.0 * x243 + 2.0 * x279
    x281 = x0 * (4.0 * x244 + 2.0 * x278 + x280) + x124 * x254
    x282 = x6 * x61
    x283 = x282 * x58
    x284 = x281 * x283
    x285 = x283 * x48
    x286 = 46.7389915737742 * x272
    x287 = x57 * x6
    x288 = x287 * x82
    x289 = x268 * x288
    x290 = x263 * x288
    x291 = x259 * x288
    x292 = 46.7389915737742 * x288
    x293 = x160 * x292
    x294 = x288 * x48
    x295 = x292 * x48
    x296 = 10.1992841329868 * x44
    x297 = x0 * (x183 * x221 + x230) + x124 * x258
    x298 = x130 * x297
    x299 = 26.9847693667702 * x42
    x300 = x0 * (x183 * x239 + x220 + x242 + x257) + x124 * x262
    x301 = x130 * x300
    x302 = x0 * (x183 * x248 + x232 + x251 + 2.0 * x261) + x124 * x268
    x303 = 34.8371874529163 * x39
    x304 = 60.3397786612521 * x71
    x305 = x0 * (x183 * x277 + x244 + 3.0 * x267 + x280) + x124 * x272
    x306 = 26.9847693667702 * x14
    x307 = x58 * x61
    x308 = x14 * x307
    x309 = x57 * x82
    x310 = x14 * x309
    x311 = 60.3397786612521 * x310
    x312 = x306 * x309
    x313 = x287 * x89
    x314 = 26.9847693667702 * x288
    x315 = 104.511562358749 * x310
    x316 = 46.7389915737742 * x310
    x317 = 10.1992841329868 * x288
    x318 = -x70 - C[2]
    x319 = x318**2 * x85
    x320 = x319 + x87
    x321 = x202 * x226
    x322 = x320 * x71
    x323 = x318 * x85
    x324 = x26 * x323
    x325 = x322 + x324
    x326 = x117 * x320
    x327 = x139 * x84
    x328 = x103 * x318
    x329 = 2.0 * x328
    x330 = x119 + x319
    x331 = x0 * (x329 + x330)
    x332 = x325 * x71
    x333 = x331 + x332
    x334 = x117 * x325
    x335 = x328 + x87
    x336 = x335 * x71
    x337 = x0 * (x103 + x323)
    x338 = 4.0 * x318 * x87 + 2.0 * x337
    x339 = x0 * (2.0 * x322 + 2.0 * x336 + x338)
    x340 = x333 * x71
    x341 = x339 + x340
    x342 = x202 * x341
    x343 = x117 * x333
    x344 = x336 + x337
    x345 = x344 * x71
    x346 = x0 * (x207 + x329)
    x347 = 3.0 * x331 + 2.0 * x346
    x348 = x0 * (3.0 * x332 + 2.0 * x345 + x347)
    x349 = x341 * x71
    x350 = x348 + x349
    x351 = x139 * x218
    x352 = x174 * x225
    x353 = x139 * x266
    x354 = x139 * x273
    x355 = x150 * x288
    x356 = x144 * x292
    x357 = x136 * x288
    x358 = x129 * x288
    x359 = x287 * x62
    x360 = x359 * x48
    x361 = x157 * x320
    x362 = x324 + x361
    x363 = 17.6656783191643 * x362
    x364 = x157 * x325
    x365 = x331 + x364
    x366 = x117 * x362
    x367 = x157 * x333
    x368 = x339 + x367
    x369 = x117 * x365
    x370 = x157 * x341
    x371 = x348 + x370
    x372 = 60.3397786612521 * x368
    x373 = x345 + x346
    x374 = x373 * x71
    x375 = x0 * (x105 + 3.0 * x336 + 3.0 * x337)
    x376 = 4.0 * x339 + 2.0 * x375
    x377 = x0 * (4.0 * x340 + 2.0 * x374 + x376) + x157 * x350
    x378 = x359 * x377
    x379 = x124 * x139
    x380 = x57 * x62
    x381 = x14 * x380
    x382 = x0 * (x206 * x323 + x330) + x157 * x362
    x383 = x139 * x382
    x384 = x0 * (x206 * x335 + x322 + x338 + x361) + x157 * x365
    x385 = x139 * x384
    x386 = x0 * (x206 * x344 + x332 + x347 + 2.0 * x364) + x157 * x368
    x387 = x0 * (x206 * x373 + x340 + 3.0 * x367 + x376) + x157 * x371

    # 270 item(s)
    result[0, 0, 0] = numpy.sum(
        x59
        * x62
        * (
            x0
            * (
                x26 * (4.0 * x35 + 4.0 * x37 + x44)
                + x33
                + x49 * (x45 + x47)
                + 4.0 * x50
                + 5.0 * x52
            )
            + x48 * x55
        )
    )
    result[0, 0, 1] = numpy.sum(x65 * x69)
    result[0, 0, 2] = numpy.sum(x69 * x72)
    result[0, 0, 3] = numpy.sum(x76 * x83)
    result[0, 0, 4] = numpy.sum(x62 * x71 * x76 * x84)
    result[0, 0, 5] = numpy.sum(x75 * x88 * x91)
    result[0, 0, 6] = numpy.sum(x100 * x96)
    result[0, 0, 7] = numpy.sum(x100 * x101 * x71)
    result[0, 0, 8] = numpy.sum(x84 * x88 * x90 * x99)
    result[0, 0, 9] = numpy.sum(x107 * x99)
    result[0, 0, 10] = numpy.sum(x112 * x116)
    result[0, 0, 11] = numpy.sum(x116 * x71 * x96)
    result[0, 0, 12] = numpy.sum(x114 * x118 * x81)
    result[0, 0, 13] = numpy.sum(x107 * x115 * x64)
    result[0, 0, 14] = numpy.sum(x115 * x123 * x90)
    result[0, 1, 0] = numpy.sum(x125 * x126)
    result[0, 1, 1] = numpy.sum(x131 * x67)
    result[0, 1, 2] = numpy.sum(x124 * x133 * x71)
    result[0, 1, 3] = numpy.sum(x136 * x137)
    result[0, 1, 4] = numpy.sum(x137 * x138 * x71)
    result[0, 1, 5] = numpy.sum(x124 * x140 * x141)
    result[0, 1, 6] = numpy.sum(x142 * x145)
    result[0, 1, 7] = numpy.sum(x135 * x146 * x147)
    result[0, 1, 8] = numpy.sum(x118 * x128 * x146)
    result[0, 1, 9] = numpy.sum(x105 * x124 * x148)
    result[0, 1, 10] = numpy.sum(x150 * x152)
    result[0, 1, 11] = numpy.sum(x113 * x153 * x71)
    result[0, 1, 12] = numpy.sum(x118 * x135 * x154)
    result[0, 1, 13] = numpy.sum(x113 * x129 * x155)
    result[0, 1, 14] = numpy.sum(x122 * x124 * x156)
    result[0, 2, 0] = numpy.sum(x126 * x158)
    result[0, 2, 1] = numpy.sum(x133 * x157 * x64)
    result[0, 2, 2] = numpy.sum(x132 * x161)
    result[0, 2, 3] = numpy.sum(x137 * x162 * x80)
    result[0, 2, 4] = numpy.sum(x161 * x163 * x74)
    result[0, 2, 5] = numpy.sum(x140 * x166)
    result[0, 2, 6] = numpy.sum(x142 * x167 * x95)
    result[0, 2, 7] = numpy.sum(x146 * x168 * x80)
    result[0, 2, 8] = numpy.sum(x146 * x165 * x169)
    result[0, 2, 9] = numpy.sum(x148 * x171)
    result[0, 2, 10] = numpy.sum(x111 * x152 * x157)
    result[0, 2, 11] = numpy.sum(x113 * x172 * x95)
    result[0, 2, 12] = numpy.sum(x154 * x173 * x80)
    result[0, 2, 13] = numpy.sum(x113 * x175 * x64)
    result[0, 2, 14] = numpy.sum(x156 * x177)
    result[0, 3, 0] = numpy.sum(x178 * x179 * x54)
    result[0, 3, 1] = numpy.sum(x180 * x182)
    result[0, 3, 2] = numpy.sum(x178 * x182 * x71)
    result[0, 3, 3] = numpy.sum(x185 * x186 * x30)
    result[0, 3, 4] = numpy.sum(x147 * x187 * x30)
    result[0, 3, 5] = numpy.sum(x178 * x188 * x30)
    result[0, 3, 6] = numpy.sum(x190 * x28)
    result[0, 3, 7] = numpy.sum(x147 * x185 * x191)
    result[0, 3, 8] = numpy.sum(x118 * x180 * x191)
    result[0, 3, 9] = numpy.sum(x178 * x192 * x28)
    result[0, 3, 10] = numpy.sum(x12 * x179 * x193)
    result[0, 3, 11] = numpy.sum(x12 * x190 * x71)
    result[0, 3, 12] = numpy.sum(x12 * x185 * x188)
    result[0, 3, 13] = numpy.sum(x12 * x180 * x192)
    result[0, 3, 14] = numpy.sum(x123 * x178 * x194)
    result[0, 4, 0] = numpy.sum(x125 * x157 * x54 * x68)
    result[0, 4, 1] = numpy.sum(x131 * x157 * x32)
    result[0, 4, 2] = numpy.sum(x161 * x195 * x32)
    result[0, 4, 3] = numpy.sum(x136 * x167 * x30)
    result[0, 4, 4] = numpy.sum(x128 * x196 * x30)
    result[0, 4, 5] = numpy.sum(x124 * x166 * x197)
    result[0, 4, 6] = numpy.sum(x153 * x157 * x28)
    result[0, 4, 7] = numpy.sum(x135 * x196 * x28)
    result[0, 4, 8] = numpy.sum(x128 * x198 * x28)
    result[0, 4, 9] = numpy.sum(x124 * x175 * x28)
    result[0, 4, 10] = numpy.sum(x12 * x130 * x150 * x158)
    result[0, 4, 11] = numpy.sum(x12 * x144 * x172)
    result[0, 4, 12] = numpy.sum(x12 * x135 * x199)
    result[0, 4, 13] = numpy.sum(x129 * x171 * x194)
    result[0, 4, 14] = numpy.sum(x12 * x125 * x139 * x177)
    result[0, 5, 0] = numpy.sum(x200 * x201 * x54)
    result[0, 5, 1] = numpy.sum(x200 * x203 * x64)
    result[0, 5, 2] = numpy.sum(x203 * x204)
    result[0, 5, 3] = numpy.sum(x205 * x30 * x81)
    result[0, 5, 4] = numpy.sum(x197 * x204 * x84)
    result[0, 5, 5] = numpy.sum(x208 * x30 * x91)
    result[0, 5, 6] = numpy.sum(x205 * x28 * x96)
    result[0, 5, 7] = numpy.sum(x191 * x209 * x80)
    result[0, 5, 8] = numpy.sum(x169 * x191 * x208)
    result[0, 5, 9] = numpy.sum(x211 * x28)
    result[0, 5, 10] = numpy.sum(x112 * x194 * x200)
    result[0, 5, 11] = numpy.sum(x194 * x204 * x96)
    result[0, 5, 12] = numpy.sum(x194 * x208 * x81)
    result[0, 5, 13] = numpy.sum(x12 * x211 * x64)
    result[0, 5, 14] = numpy.sum(x12 * x201 * x212)
    result[1, 0, 0] = numpy.sum(x179 * x215 * x219)
    result[1, 0, 1] = numpy.sum(x223 * x227)
    result[1, 0, 2] = numpy.sum(x215 * x227 * x71)
    result[1, 0, 3] = numpy.sum(x186 * x233 * x235)
    result[1, 0, 4] = numpy.sum(x147 * x235 * x236)
    result[1, 0, 5] = numpy.sum(x188 * x215 * x235)
    result[1, 0, 6] = numpy.sum(x238 * x246)
    result[1, 0, 7] = numpy.sum(x147 * x238 * x247)
    result[1, 0, 8] = numpy.sum(x118 * x236 * x238)
    result[1, 0, 9] = numpy.sum(x192 * x215 * x238)
    result[1, 0, 10] = numpy.sum(x179 * x254 * x255)
    result[1, 0, 11] = numpy.sum(x246 * x255 * x71)
    result[1, 0, 12] = numpy.sum(x188 * x233 * x255)
    result[1, 0, 13] = numpy.sum(x192 * x223 * x255)
    result[1, 0, 14] = numpy.sum(x123 * x255 * x256)
    result[1, 1, 0] = numpy.sum(x259 * x260)
    result[1, 1, 1] = numpy.sum(x263 * x264)
    result[1, 1, 2] = numpy.sum(x258 * x265 * x71)
    result[1, 1, 3] = numpy.sum(x266 * x269)
    result[1, 1, 4] = numpy.sum(x147 * x234 * x270)
    result[1, 1, 5] = numpy.sum(x118 * x258 * x266)
    result[1, 1, 6] = numpy.sum(x272 * x274)
    result[1, 1, 7] = numpy.sum(x269 * x275 * x71)
    result[1, 1, 8] = numpy.sum(x118 * x262 * x275)
    result[1, 1, 9] = numpy.sum(x155 * x258 * x273)
    result[1, 1, 10] = numpy.sum(x276 * x284)
    result[1, 1, 11] = numpy.sum(x285 * x286 * x71)
    result[1, 1, 12] = numpy.sum(x141 * x289 * x48)
    result[1, 1, 13] = numpy.sum(x105 * x290 * x48)
    result[1, 1, 14] = numpy.sum(x122 * x291 * x48)
    result[1, 2, 0] = numpy.sum(x158 * x215 * x260)
    result[1, 2, 1] = numpy.sum(x157 * x223 * x265)
    result[1, 2, 2] = numpy.sum(x172 * x215 * x225)
    result[1, 2, 3] = numpy.sum(x167 * x233 * x266)
    result[1, 2, 4] = numpy.sum(x196 * x223 * x234)
    result[1, 2, 5] = numpy.sum(x199 * x215 * x234)
    result[1, 2, 6] = numpy.sum(x157 * x245 * x274)
    result[1, 2, 7] = numpy.sum(x196 * x233 * x237)
    result[1, 2, 8] = numpy.sum(x198 * x223 * x237)
    result[1, 2, 9] = numpy.sum(x171 * x256 * x273)
    result[1, 2, 10] = numpy.sum(x158 * x254 * x285)
    result[1, 2, 11] = numpy.sum(x245 * x293 * x48)
    result[1, 2, 12] = numpy.sum(x166 * x233 * x294)
    result[1, 2, 13] = numpy.sum(x171 * x223 * x295)
    result[1, 2, 14] = numpy.sum(x177 * x215 * x276 * x288)
    result[1, 3, 0] = numpy.sum(x296 * x298)
    result[1, 3, 1] = numpy.sum(x299 * x301)
    result[1, 3, 2] = numpy.sum(x298 * x299 * x71)
    result[1, 3, 3] = numpy.sum(x130 * x302 * x303)
    result[1, 3, 4] = numpy.sum(x301 * x304 * x39)
    result[1, 3, 5] = numpy.sum(x188 * x297 * x39)
    result[1, 3, 6] = numpy.sum(x305 * x306 * x307)
    result[1, 3, 7] = numpy.sum(x302 * x304 * x308)
    result[1, 3, 8] = numpy.sum(x300 * x311 * x88)
    result[1, 3, 9] = numpy.sum(x105 * x297 * x312)
    result[1, 3, 10] = numpy.sum(
        x282
        * x59
        * (
            x0
            * (
                x183 * (x278 + x279)
                + 5.0 * x252
                + x253
                + x26 * (x111 + 4.0 * x249 + 4.0 * x250)
                + 4.0 * x271
            )
            + x124 * x281
        )
    )
    result[1, 3, 11] = numpy.sum(x283 * x305 * x72)
    result[1, 3, 12] = numpy.sum(x302 * x313 * x88)
    result[1, 3, 13] = numpy.sum(x105 * x300 * x314)
    result[1, 3, 14] = numpy.sum(x123 * x288 * x297)
    result[1, 4, 0] = numpy.sum(x167 * x259 * x44)
    result[1, 4, 1] = numpy.sum(x167 * x263 * x42)
    result[1, 4, 2] = numpy.sum(x172 * x258 * x42)
    result[1, 4, 3] = numpy.sum(x162 * x269 * x39)
    result[1, 4, 4] = numpy.sum(x196 * x262 * x39)
    result[1, 4, 5] = numpy.sum(x199 * x258 * x39)
    result[1, 4, 6] = numpy.sum(x157 * x286 * x308)
    result[1, 4, 7] = numpy.sum(x160 * x268 * x315)
    result[1, 4, 8] = numpy.sum(x165 * x270 * x310)
    result[1, 4, 9] = numpy.sum(x171 * x258 * x316)
    result[1, 4, 10] = numpy.sum(x158 * x284)
    result[1, 4, 11] = numpy.sum(x272 * x293)
    result[1, 4, 12] = numpy.sum(x166 * x289)
    result[1, 4, 13] = numpy.sum(x171 * x290)
    result[1, 4, 14] = numpy.sum(x177 * x291)
    result[1, 5, 0] = numpy.sum(x205 * x215 * x296)
    result[1, 5, 1] = numpy.sum(x205 * x223 * x299)
    result[1, 5, 2] = numpy.sum(x209 * x215 * x299)
    result[1, 5, 3] = numpy.sum(x205 * x233 * x303)
    result[1, 5, 4] = numpy.sum(x209 * x236 * x39)
    result[1, 5, 5] = numpy.sum(x208 * x256 * x303)
    result[1, 5, 6] = numpy.sum(x200 * x245 * x312)
    result[1, 5, 7] = numpy.sum(x204 * x247 * x310)
    result[1, 5, 8] = numpy.sum(x208 * x236 * x310)
    result[1, 5, 9] = numpy.sum(x210 * x215 * x312)
    result[1, 5, 10] = numpy.sum(x200 * x254 * x317)
    result[1, 5, 11] = numpy.sum(x204 * x245 * x314)
    result[1, 5, 12] = numpy.sum(x208 * x233 * x313)
    result[1, 5, 13] = numpy.sum(x210 * x223 * x314)
    result[1, 5, 14] = numpy.sum(x212 * x215 * x317)
    result[2, 0, 0] = numpy.sum(x201 * x219 * x320)
    result[2, 0, 1] = numpy.sum(x320 * x321 * x64)
    result[2, 0, 2] = numpy.sum(x321 * x325)
    result[2, 0, 3] = numpy.sum(x235 * x326 * x81)
    result[2, 0, 4] = numpy.sum(x235 * x325 * x327)
    result[2, 0, 5] = numpy.sum(x235 * x333 * x91)
    result[2, 0, 6] = numpy.sum(x238 * x326 * x96)
    result[2, 0, 7] = numpy.sum(x101 * x238 * x334)
    result[2, 0, 8] = numpy.sum(x238 * x327 * x333)
    result[2, 0, 9] = numpy.sum(x238 * x342)
    result[2, 0, 10] = numpy.sum(x112 * x255 * x326)
    result[2, 0, 11] = numpy.sum(x255 * x334 * x96)
    result[2, 0, 12] = numpy.sum(x255 * x343 * x81)
    result[2, 0, 13] = numpy.sum(x255 * x342 * x64)
    result[2, 0, 14] = numpy.sum(x201 * x255 * x350)
    result[2, 1, 0] = numpy.sum(x125 * x320 * x351)
    result[2, 1, 1] = numpy.sum(x129 * x225 * x326)
    result[2, 1, 2] = numpy.sum(x124 * x325 * x352)
    result[2, 1, 3] = numpy.sum(x135 * x266 * x326)
    result[2, 1, 4] = numpy.sum(x138 * x234 * x334)
    result[2, 1, 5] = numpy.sum(x124 * x333 * x353)
    result[2, 1, 6] = numpy.sum(x144 * x273 * x326)
    result[2, 1, 7] = numpy.sum(x135 * x275 * x334)
    result[2, 1, 8] = numpy.sum(x128 * x275 * x343)
    result[2, 1, 9] = numpy.sum(x124 * x341 * x354)
    result[2, 1, 10] = numpy.sum(x276 * x320 * x355)
    result[2, 1, 11] = numpy.sum(x325 * x356 * x48)
    result[2, 1, 12] = numpy.sum(x333 * x357 * x48)
    result[2, 1, 13] = numpy.sum(x341 * x358 * x48)
    result[2, 1, 14] = numpy.sum(x125 * x350 * x360)
    result[2, 2, 0] = numpy.sum(x351 * x363)
    result[2, 2, 1] = numpy.sum(x352 * x362 * x64)
    result[2, 2, 2] = numpy.sum(x352 * x365)
    result[2, 2, 3] = numpy.sum(x266 * x366 * x80)
    result[2, 2, 4] = numpy.sum(x139 * x163 * x234 * x365)
    result[2, 2, 5] = numpy.sum(x353 * x368)
    result[2, 2, 6] = numpy.sum(x273 * x366 * x95)
    result[2, 2, 7] = numpy.sum(x275 * x369 * x80)
    result[2, 2, 8] = numpy.sum(x169 * x275 * x368)
    result[2, 2, 9] = numpy.sum(x354 * x371)
    result[2, 2, 10] = numpy.sum(x111 * x294 * x363)
    result[2, 2, 11] = numpy.sum(x295 * x365 * x95)
    result[2, 2, 12] = numpy.sum(x294 * x372 * x80)
    result[2, 2, 13] = numpy.sum(46.7389915737742 * x360 * x371 * x64)
    result[2, 2, 14] = numpy.sum(x276 * x378)
    result[2, 3, 0] = numpy.sum(x178 * x296 * x326)
    result[2, 3, 1] = numpy.sum(x180 * x299 * x326)
    result[2, 3, 2] = numpy.sum(x178 * x299 * x334)
    result[2, 3, 3] = numpy.sum(x185 * x303 * x326)
    result[2, 3, 4] = numpy.sum(x187 * x334 * x39)
    result[2, 3, 5] = numpy.sum(x178 * x303 * x343)
    result[2, 3, 6] = numpy.sum(x189 * x312 * x320)
    result[2, 3, 7] = numpy.sum(x185 * x311 * x325)
    result[2, 3, 8] = numpy.sum(x187 * x310 * x333)
    result[2, 3, 9] = numpy.sum(x178 * x312 * x341)
    result[2, 3, 10] = numpy.sum(x193 * x317 * x320)
    result[2, 3, 11] = numpy.sum(x189 * x314 * x325)
    result[2, 3, 12] = numpy.sum(x185 * x313 * x333)
    result[2, 3, 13] = numpy.sum(x180 * x314 * x341)
    result[2, 3, 14] = numpy.sum(x178 * x317 * x350)
    result[2, 4, 0] = numpy.sum(x363 * x379 * x44)
    result[2, 4, 1] = numpy.sum(x129 * x366 * x42)
    result[2, 4, 2] = numpy.sum(x124 * x174 * x365 * x42)
    result[2, 4, 3] = numpy.sum(x136 * x366 * x39)
    result[2, 4, 4] = numpy.sum(x138 * x369 * x39)
    result[2, 4, 5] = numpy.sum(x372 * x379 * x39)
    result[2, 4, 6] = numpy.sum(x144 * x316 * x362)
    result[2, 4, 7] = numpy.sum(x135 * x315 * x365)
    result[2, 4, 8] = numpy.sum(x138 * x310 * x368)
    result[2, 4, 9] = numpy.sum(x195 * x371 * x381)
    result[2, 4, 10] = numpy.sum(x355 * x363)
    result[2, 4, 11] = numpy.sum(x356 * x365)
    result[2, 4, 12] = numpy.sum(x357 * x368)
    result[2, 4, 13] = numpy.sum(x358 * x371)
    result[2, 4, 14] = numpy.sum(x125 * x378)
    result[2, 5, 0] = numpy.sum(x296 * x383)
    result[2, 5, 1] = numpy.sum(x299 * x383 * x64)
    result[2, 5, 2] = numpy.sum(x299 * x385)
    result[2, 5, 3] = numpy.sum(x117 * x382 * x39 * x81)
    result[2, 5, 4] = numpy.sum(x385 * x39 * x84)
    result[2, 5, 5] = numpy.sum(x139 * x303 * x386)
    result[2, 5, 6] = numpy.sum(x310 * x382 * x96)
    result[2, 5, 7] = numpy.sum(x101 * x310 * x384)
    result[2, 5, 8] = numpy.sum(x381 * x386 * x84)
    result[2, 5, 9] = numpy.sum(x306 * x380 * x387)
    result[2, 5, 10] = numpy.sum(x112 * x288 * x382)
    result[2, 5, 11] = numpy.sum(x288 * x384 * x96)
    result[2, 5, 12] = numpy.sum(x287 * x386 * x83)
    result[2, 5, 13] = numpy.sum(x359 * x387 * x65)
    result[2, 5, 14] = numpy.sum(
        10.1992841329868
        * x359
        * (
            x0
            * (
                x206 * (x374 + x375)
                + x26 * (x122 + 4.0 * x345 + 4.0 * x346)
                + 5.0 * x348
                + x349
                + 4.0 * x370
            )
            + x157 * x377
        )
    )
    return result


def diag_quadrupole3d_30(ax, da, A, bx, db, B, C):
    """Cartesian 3D (fs) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 10, 1), dtype=float)

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
    x10 = -x4 - C[0]
    x11 = x3 * x7
    x12 = x10 * x11
    x13 = 2.0 * x0
    x14 = x10**2 * x11
    x15 = x0 * x11
    x16 = x14 + x15
    x17 = x16 * x5
    x18 = 2.0 * x5
    x19 = x12 * x13 + x17
    x20 = x0 * (x12 * x18 + x14 + 3.0 * x15) + x19 * x5
    x21 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x22 = da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**1.5)
    x23 = x21 * x22
    x24 = 5.84237394672177 * x23
    x25 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x26 = 0.564189583547756 * x1
    x27 = x25 * x26
    x28 = -x1 * (ax * A[1] + bx * B[1])
    x29 = -x28 - A[1]
    x30 = x27 * x29
    x31 = 13.0639452948436 * x23
    x32 = x20 * x31
    x33 = -x1 * (ax * A[2] + bx * B[2])
    x34 = -x33 - A[2]
    x35 = 0.318309886183791 * x2
    x36 = x19 * x35
    x37 = x25 * x3
    x38 = x0 * x37
    x39 = x29**2 * x37 + x38
    x40 = x31 * x39
    x41 = 22.6274169979695 * x30
    x42 = x23 * x34
    x43 = x22 * x25
    x44 = x21 * x3
    x45 = x0 * x44
    x46 = x34**2 * x44 + x45
    x47 = 13.0639452948436 * x46
    x48 = x43 * x47
    x49 = x29 * x37
    x50 = x13 * x49 + x29 * x39
    x51 = x16 * x35
    x52 = x34 * x44
    x53 = 5.84237394672177 * x13 * x52 + 5.84237394672177 * x34 * x46
    x54 = x11 * x5**2 + x15
    x55 = x13 * x9 + x5 * x54
    x56 = -x28 - C[1]
    x57 = x37 * x56**2
    x58 = x38 + x57
    x59 = x35 * x58
    x60 = x29 * x58
    x61 = x37 * x56
    x62 = x13 * x61 + x60
    x63 = x35 * x62
    x64 = x31 * x54
    x65 = 2.0 * x29
    x66 = x0 * (3.0 * x38 + x57 + x61 * x65) + x29 * x62
    x67 = x31 * x66
    x68 = x26 * x8
    x69 = x22 * x8
    x70 = x26 * x7
    x71 = x22 * x7
    x72 = x35 * x43
    x73 = -x33 - C[2]
    x74 = x44 * x73**2
    x75 = x45 + x74
    x76 = 5.84237394672177 * x75
    x77 = x54 * x72
    x78 = 13.0639452948436 * x75
    x79 = x34 * x75
    x80 = x44 * x73
    x81 = x13 * x80 + x79
    x82 = 13.0639452948436 * x81
    x83 = 2.0 * x34
    x84 = x0 * (3.0 * x45 + x74 + x80 * x83) + x34 * x81
    x85 = 13.0639452948436 * x84
    x86 = x35 * x71
    x87 = x27 * x71

    # 30 item(s)
    result[0, 0, 0] = numpy.sum(
        x24
        * x27
        * (
            x0 * (4.0 * x0 * x12 + x13 * (x12 + x9) + 2.0 * x17 + x18 * (x10 * x9 + x15))
            + x20 * x5
        )
    )
    result[0, 1, 0] = numpy.sum(x30 * x32)
    result[0, 2, 0] = numpy.sum(x27 * x32 * x34)
    result[0, 3, 0] = numpy.sum(x36 * x40)
    result[0, 4, 0] = numpy.sum(x19 * x41 * x42)
    result[0, 5, 0] = numpy.sum(x36 * x48)
    result[0, 6, 0] = numpy.sum(x24 * x50 * x51)
    result[0, 7, 0] = numpy.sum(x34 * x40 * x51)
    result[0, 8, 0] = numpy.sum(x29 * x48 * x51)
    result[0, 9, 0] = numpy.sum(x43 * x51 * x53)
    result[1, 0, 0] = numpy.sum(x24 * x55 * x59)
    result[1, 1, 0] = numpy.sum(x63 * x64)
    result[1, 2, 0] = numpy.sum(x34 * x59 * x64)
    result[1, 3, 0] = numpy.sum(x67 * x68)
    result[1, 4, 0] = numpy.sum(22.6274169979695 * x42 * x62 * x68)
    result[1, 5, 0] = numpy.sum(x47 * x59 * x69)
    result[1, 6, 0] = numpy.sum(
        x24
        * x70
        * (
            x0
            * (x13 * (x49 + x61) + 4.0 * x38 * x56 + 2.0 * x60 + x65 * (x38 + x49 * x56))
            + x29 * x66
        )
    )
    result[1, 7, 0] = numpy.sum(x34 * x67 * x70)
    result[1, 8, 0] = numpy.sum(x47 * x63 * x71)
    result[1, 9, 0] = numpy.sum(x53 * x59 * x71)
    result[2, 0, 0] = numpy.sum(x55 * x72 * x76)
    result[2, 1, 0] = numpy.sum(x29 * x77 * x78)
    result[2, 2, 0] = numpy.sum(x77 * x82)
    result[2, 3, 0] = numpy.sum(x35 * x39 * x69 * x78)
    result[2, 4, 0] = numpy.sum(x41 * x69 * x81)
    result[2, 5, 0] = numpy.sum(x27 * x69 * x85)
    result[2, 6, 0] = numpy.sum(x50 * x76 * x86)
    result[2, 7, 0] = numpy.sum(x39 * x82 * x86)
    result[2, 8, 0] = numpy.sum(x29 * x85 * x87)
    result[2, 9, 0] = numpy.sum(
        5.84237394672177
        * x87
        * (
            x0
            * (x13 * (x52 + x80) + 4.0 * x45 * x73 + 2.0 * x79 + x83 * (x45 + x52 * x73))
            + x34 * x84
        )
    )
    return result


def diag_quadrupole3d_31(ax, da, A, bx, db, B, C):
    """Cartesian 3D (fp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 10, 3), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - C[0]
    x4 = -x2 - B[0]
    x5 = ax * bx * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(x1)
    x8 = 1.77245385090552 * x7
    x9 = x6 * x8
    x10 = x4 * x9
    x11 = x10 * x3
    x12 = x3**2 * x9
    x13 = x0 * x9
    x14 = 3.0 * x13
    x15 = x12 + x14
    x16 = x0 * (2.0 * x11 + x15)
    x17 = -x2 - A[0]
    x18 = x10 * x17
    x19 = x3 * x9
    x20 = x17 * x19
    x21 = 2.0 * x0
    x22 = x0 * (x10 + x19)
    x23 = x17 * (x11 + x13)
    x24 = 2.0 * x17
    x25 = x12 + x13
    x26 = x25 * x4
    x27 = x19 * x21
    x28 = x26 + x27
    x29 = x17 * x28
    x30 = x17 * x25
    x31 = x27 + x30
    x32 = x0 * (x15 + x19 * x24) + x17 * x31
    x33 = 4.0 * x13 * x3
    x34 = x16 + x29
    x35 = x0 * (2.0 * x22 + 2.0 * x23 + x26 + x30 + x33) + x17 * x34
    x36 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x37 = da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**2.5)
    x38 = 11.6847478934435 * x37
    x39 = x36 * x38
    x40 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x41 = 0.564189583547756 * x1
    x42 = x40 * x41
    x43 = x39 * x42
    x44 = -x1 * (ax * A[1] + bx * B[1])
    x45 = -x44 - B[1]
    x46 = x17 * x9
    x47 = x43 * (
        x0 * (x21 * (x19 + x46) + x24 * (x13 + x20) + 2.0 * x30 + x33) + x17 * x32
    )
    x48 = -x1 * (ax * A[2] + bx * B[2])
    x49 = -x48 - B[2]
    x50 = -x44 - A[1]
    x51 = x42 * x50
    x52 = 26.1278905896872 * x37
    x53 = x36 * x52
    x54 = x35 * x53
    x55 = x0 * x8
    x56 = x40 * x55
    x57 = x40 * x8
    x58 = x50 * x57
    x59 = x45 * x58
    x60 = x56 + x59
    x61 = 0.318309886183791 * x7
    x62 = x60 * x61
    x63 = x32 * x53
    x64 = -x48 - A[2]
    x65 = x42 * x64
    x66 = x36 * x55
    x67 = x36 * x8
    x68 = x64 * x67
    x69 = x49 * x68
    x70 = x66 + x69
    x71 = x40 * x61
    x72 = x52 * x71
    x73 = x50**2 * x57
    x74 = x56 + x73
    x75 = x53 * x61
    x76 = x74 * x75
    x77 = x45 * x57
    x78 = x0 * (x58 + x77) + x50 * x60
    x79 = x75 * x78
    x80 = 45.2548339959391 * x37
    x81 = x36 * x64
    x82 = x80 * x81
    x83 = x31 * x80
    x84 = x50 * x71
    x85 = x64**2 * x67
    x86 = x66 + x85
    x87 = x72 * x86
    x88 = x49 * x67
    x89 = x0 * (x68 + x88) + x64 * x70
    x90 = x72 * x89
    x91 = x21 * x58 + x50 * x74
    x92 = x39 * x61
    x93 = x91 * x92
    x94 = 3.0 * x56
    x95 = 2.0 * x50
    x96 = x0 * (x73 + x77 * x95 + x94) + x50 * x78
    x97 = 0.179587122125167 * x52
    x98 = x25 * x97
    x99 = x21 * x68 + x64 * x86
    x100 = x38 * x71
    x101 = x100 * x99
    x102 = 3.0 * x66
    x103 = 2.0 * x64
    x104 = x0 * (x102 + x103 * x88 + x85) + x64 * x89
    x105 = x17**2 * x9
    x106 = x13 + x18
    x107 = x0 * (x10 + x46) + x106 * x17
    x108 = x0 * (x10 * x24 + x105 + x14) + x107 * x17
    x109 = -x44 - C[1]
    x110 = x109**2 * x57
    x111 = x110 + x56
    x112 = x111 * x92
    x113 = x111 * x45
    x114 = x109 * x57
    x115 = x114 * x21
    x116 = x113 + x115
    x117 = x105 + x13
    x118 = x117 * x17 + x21 * x46
    x119 = x111 * x50
    x120 = x115 + x119
    x121 = x120 * x75
    x122 = x109 * x77
    x123 = x110 + x94
    x124 = x0 * (2.0 * x122 + x123)
    x125 = x116 * x50
    x126 = x124 + x125
    x127 = x117 * x75
    x128 = x111 * x97
    x129 = x0 * (x114 * x95 + x123) + x120 * x50
    x130 = x41 * x6
    x131 = x130 * x17
    x132 = x0 * (x114 + x77)
    x133 = x50 * (x122 + x56)
    x134 = 4.0 * x109 * x56
    x135 = x0 * (x113 + x119 + 2.0 * x132 + 2.0 * x133 + x134) + x126 * x50
    x136 = x135 * x53
    x137 = x129 * x53
    x138 = x120 * x80
    x139 = x6 * x61
    x140 = x139 * x52
    x141 = x140 * x86
    x142 = x140 * x89
    x143 = x109 * x58
    x144 = x130 * x39
    x145 = x144 * (
        x0 * (2.0 * x119 + x134 + x21 * (x114 + x58) + x95 * (x143 + x56)) + x129 * x50
    )
    x146 = x130 * x64
    x147 = x139 * x38
    x148 = x147 * x99
    x149 = -x48 - C[2]
    x150 = x149**2 * x67
    x151 = x150 + x66
    x152 = x100 * x151
    x153 = x151 * x49
    x154 = x149 * x67
    x155 = x154 * x21
    x156 = x153 + x155
    x157 = x107 * x72
    x158 = x151 * x97
    x159 = x117 * x72
    x160 = x151 * x64
    x161 = x155 + x160
    x162 = x149 * x88
    x163 = x102 + x150
    x164 = x0 * (2.0 * x162 + x163)
    x165 = x156 * x64
    x166 = x164 + x165
    x167 = x140 * x78
    x168 = x140 * x74
    x169 = x161 * x80
    x170 = x6 * x62
    x171 = x42 * x6
    x172 = x17 * x171
    x173 = x0 * (x103 * x154 + x163) + x161 * x64
    x174 = x173 * x52
    x175 = x0 * (x154 + x88)
    x176 = x64 * (x162 + x66)
    x177 = 4.0 * x149 * x66
    x178 = x0 * (x153 + x160 + 2.0 * x175 + 2.0 * x176 + x177) + x166 * x64
    x179 = x178 * x52
    x180 = x147 * x151
    x181 = x171 * x50
    x182 = x149 * x68
    x183 = x171 * x38
    x184 = x183 * (
        x0 * (x103 * (x182 + x66) + 2.0 * x160 + x177 + x21 * (x154 + x68)) + x173 * x64
    )

    # 90 item(s)
    result[0, 0, 0] = numpy.sum(
        x43
        * (
            x0
            * (
                2.0 * x16
                + x21 * (x11 + x14 + x18 + x20)
                + x24 * (x22 + x23)
                + 2.0 * x29
                + x32
            )
            + x17 * x35
        )
    )
    result[0, 0, 1] = numpy.sum(x45 * x47)
    result[0, 0, 2] = numpy.sum(x47 * x49)
    result[0, 1, 0] = numpy.sum(x51 * x54)
    result[0, 1, 1] = numpy.sum(x62 * x63)
    result[0, 1, 2] = numpy.sum(x49 * x51 * x63)
    result[0, 2, 0] = numpy.sum(x54 * x65)
    result[0, 2, 1] = numpy.sum(x45 * x63 * x65)
    result[0, 2, 2] = numpy.sum(x32 * x70 * x72)
    result[0, 3, 0] = numpy.sum(x34 * x76)
    result[0, 3, 1] = numpy.sum(x31 * x79)
    result[0, 3, 2] = numpy.sum(x31 * x49 * x76)
    result[0, 4, 0] = numpy.sum(x34 * x51 * x82)
    result[0, 4, 1] = numpy.sum(x62 * x81 * x83)
    result[0, 4, 2] = numpy.sum(x70 * x83 * x84)
    result[0, 5, 0] = numpy.sum(x34 * x87)
    result[0, 5, 1] = numpy.sum(x31 * x45 * x87)
    result[0, 5, 2] = numpy.sum(x31 * x90)
    result[0, 6, 0] = numpy.sum(x28 * x93)
    result[0, 6, 1] = numpy.sum(x25 * x92 * x96)
    result[0, 6, 2] = numpy.sum(x25 * x49 * x93)
    result[0, 7, 0] = numpy.sum(x28 * x64 * x76)
    result[0, 7, 1] = numpy.sum(x25 * x64 * x79)
    result[0, 7, 2] = numpy.sum(x70 * x74 * x98)
    result[0, 8, 0] = numpy.sum(x28 * x50 * x87)
    result[0, 8, 1] = numpy.sum(x60 * x86 * x98)
    result[0, 8, 2] = numpy.sum(x25 * x50 * x90)
    result[0, 9, 0] = numpy.sum(x101 * x28)
    result[0, 9, 1] = numpy.sum(x101 * x25 * x45)
    result[0, 9, 2] = numpy.sum(x100 * x104 * x25)
    result[1, 0, 0] = numpy.sum(x108 * x112)
    result[1, 0, 1] = numpy.sum(x116 * x118 * x92)
    result[1, 0, 2] = numpy.sum(x112 * x118 * x49)
    result[1, 1, 0] = numpy.sum(x107 * x121)
    result[1, 1, 1] = numpy.sum(x126 * x127)
    result[1, 1, 2] = numpy.sum(x117 * x121 * x49)
    result[1, 2, 0] = numpy.sum(x107 * x111 * x64 * x75)
    result[1, 2, 1] = numpy.sum(x116 * x127 * x64)
    result[1, 2, 2] = numpy.sum(x117 * x128 * x70)
    result[1, 3, 0] = numpy.sum(x106 * x129 * x75)
    result[1, 3, 1] = numpy.sum(x131 * x136)
    result[1, 3, 2] = numpy.sum(x131 * x137 * x49)
    result[1, 4, 0] = numpy.sum(x106 * x138 * x61 * x81)
    result[1, 4, 1] = numpy.sum(x126 * x131 * x82)
    result[1, 4, 2] = numpy.sum(x138 * x139 * x17 * x70)
    result[1, 5, 0] = numpy.sum(x106 * x128 * x86)
    result[1, 5, 1] = numpy.sum(x116 * x141 * x17)
    result[1, 5, 2] = numpy.sum(x111 * x142 * x17)
    result[1, 6, 0] = numpy.sum(x145 * x4)
    result[1, 6, 1] = numpy.sum(
        x144
        * (
            x0
            * (
                2.0 * x124
                + 2.0 * x125
                + x129
                + x21 * (x122 + x143 + x59 + x94)
                + x95 * (x132 + x133)
            )
            + x135 * x50
        )
    )
    result[1, 6, 2] = numpy.sum(x145 * x49)
    result[1, 7, 0] = numpy.sum(x137 * x146 * x4)
    result[1, 7, 1] = numpy.sum(x136 * x146)
    result[1, 7, 2] = numpy.sum(x129 * x140 * x70)
    result[1, 8, 0] = numpy.sum(x120 * x141 * x4)
    result[1, 8, 1] = numpy.sum(x126 * x141)
    result[1, 8, 2] = numpy.sum(x120 * x142)
    result[1, 9, 0] = numpy.sum(x111 * x148 * x4)
    result[1, 9, 1] = numpy.sum(x116 * x148)
    result[1, 9, 2] = numpy.sum(x104 * x111 * x147)
    result[2, 0, 0] = numpy.sum(x108 * x152)
    result[2, 0, 1] = numpy.sum(x118 * x152 * x45)
    result[2, 0, 2] = numpy.sum(x100 * x118 * x156)
    result[2, 1, 0] = numpy.sum(x151 * x157 * x50)
    result[2, 1, 1] = numpy.sum(x117 * x158 * x60)
    result[2, 1, 2] = numpy.sum(x156 * x159 * x50)
    result[2, 2, 0] = numpy.sum(x157 * x161)
    result[2, 2, 1] = numpy.sum(x159 * x161 * x45)
    result[2, 2, 2] = numpy.sum(x159 * x166)
    result[2, 3, 0] = numpy.sum(x106 * x158 * x74)
    result[2, 3, 1] = numpy.sum(x151 * x167 * x17)
    result[2, 3, 2] = numpy.sum(x156 * x168 * x17)
    result[2, 4, 0] = numpy.sum(x106 * x169 * x84)
    result[2, 4, 1] = numpy.sum(x169 * x17 * x170)
    result[2, 4, 2] = numpy.sum(x166 * x172 * x50 * x80)
    result[2, 5, 0] = numpy.sum(x106 * x173 * x72)
    result[2, 5, 1] = numpy.sum(x172 * x174 * x45)
    result[2, 5, 2] = numpy.sum(x172 * x179)
    result[2, 6, 0] = numpy.sum(x180 * x4 * x91)
    result[2, 6, 1] = numpy.sum(x180 * x96)
    result[2, 6, 2] = numpy.sum(x147 * x156 * x91)
    result[2, 7, 0] = numpy.sum(x161 * x168 * x4)
    result[2, 7, 1] = numpy.sum(x161 * x167)
    result[2, 7, 2] = numpy.sum(x166 * x168)
    result[2, 8, 0] = numpy.sum(x174 * x181 * x4)
    result[2, 8, 1] = numpy.sum(x170 * x174)
    result[2, 8, 2] = numpy.sum(x179 * x181)
    result[2, 9, 0] = numpy.sum(x184 * x4)
    result[2, 9, 1] = numpy.sum(x184 * x45)
    result[2, 9, 2] = numpy.sum(
        x183
        * (
            x0
            * (
                x103 * (x175 + x176)
                + 2.0 * x164
                + 2.0 * x165
                + x173
                + x21 * (x102 + x162 + x182 + x69)
            )
            + x178 * x64
        )
    )
    return result


def diag_quadrupole3d_32(ax, da, A, bx, db, B, C):
    """Cartesian 3D (fd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 10, 6), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = -x2 - C[0]
    x5 = ax * bx * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = numpy.sqrt(x1)
    x8 = 1.77245385090552 * x7
    x9 = x6 * x8
    x10 = x4**2 * x9
    x11 = x0 * x9
    x12 = x10 + x11
    x13 = x12 * x3
    x14 = x3 * x6
    x15 = x14 * x8
    x16 = x15 * x4
    x17 = x11 + x16
    x18 = x17 * x3
    x19 = x4 * x9
    x20 = x0 * (x15 + x19)
    x21 = 4.0 * x0
    x22 = x19 * x21
    x23 = 2.0 * x20 + x22
    x24 = x0 * (2.0 * x13 + 2.0 * x18 + x23)
    x25 = -x2 - A[0]
    x26 = x17 * x25
    x27 = 2.0 * x26
    x28 = x3**2 * x9
    x29 = x11 + x28
    x30 = x25 * x29
    x31 = 2.0 * x0
    x32 = x15 * x31 + x30
    x33 = x12 * x25
    x34 = x0 * (x13 + x23 + x27 + x33)
    x35 = 3.0 * x11
    x36 = 2.0 * x16 + x35
    x37 = x0 * (x28 + x36)
    x38 = x25 * (x18 + x20)
    x39 = 2.0 * x25
    x40 = x0 * (x10 + x36)
    x41 = x19 * x31
    x42 = x13 + x41
    x43 = x25 * x42
    x44 = x40 + x43
    x45 = x25 * x44
    x46 = x3 * x42
    x47 = x40 + x46
    x48 = x25 * x47
    x49 = 2.0 * x43
    x50 = x24 + x48
    x51 = x0 * (2.0 * x37 + 2.0 * x38 + 3.0 * x40 + x46 + x49) + x25 * x50
    x52 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x53 = da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**3.5)
    x54 = x52 * x53
    x55 = 13.4923846833851 * x54
    x56 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x57 = 0.564189583547756 * x1
    x58 = x56 * x57
    x59 = -x1 * (ax * A[1] + bx * B[1])
    x60 = -x59 - B[1]
    x61 = 23.3694957868871 * x60
    x62 = x15 * x25
    x63 = x19 * x25
    x64 = x33 + x41
    x65 = x0 * (x10 + x19 * x39 + x35) + x25 * x64
    x66 = x34 + x45
    x67 = x54 * x58
    x68 = x67 * (
        x0 * (x31 * (x16 + x35 + x62 + x63) + x39 * (x20 + x26) + 2.0 * x40 + x49 + x65)
        + x25 * x66
    )
    x69 = -x1 * (ax * A[2] + bx * B[2])
    x70 = -x69 - B[2]
    x71 = 23.3694957868871 * x70
    x72 = x25 * x9
    x73 = x0 * (x22 + x31 * (x19 + x72) + 2.0 * x33 + x39 * (x11 + x63)) + x25 * x65
    x74 = x54 * x73
    x75 = x56 * x8
    x76 = x60**2 * x75
    x77 = x0 * x75
    x78 = x76 + x77
    x79 = 13.4923846833851 * x78
    x80 = 0.318309886183791 * x7
    x81 = x79 * x80
    x82 = x52 * x8
    x83 = x70**2 * x82
    x84 = x0 * x82
    x85 = x83 + x84
    x86 = 13.4923846833851 * x85
    x87 = x53 * x80
    x88 = x56 * x87
    x89 = -x59 - A[1]
    x90 = 30.169889330626 * x89
    x91 = x51 * x67
    x92 = x75 * x89
    x93 = x60 * x92
    x94 = x77 + x93
    x95 = 52.2557811793745 * x94
    x96 = x54 * x80
    x97 = x95 * x96
    x98 = 52.2557811793745 * x70
    x99 = x66 * x67
    x100 = x78 * x89
    x101 = x60 * x75
    x102 = x100 + x101 * x31
    x103 = 30.169889330626 * x96
    x104 = x103 * x65
    x105 = x65 * x88
    x106 = -x69 - A[2]
    x107 = 30.169889330626 * x106
    x108 = 52.2557811793745 * x106
    x109 = x106 * x82
    x110 = x109 * x70
    x111 = x110 + x84
    x112 = 52.2557811793745 * x111
    x113 = x112 * x88
    x114 = x106 * x85
    x115 = x70 * x82
    x116 = x114 + x115 * x31
    x117 = 30.169889330626 * x116
    x118 = x75 * x89**2
    x119 = x118 + x77
    x120 = 30.169889330626 * x119
    x121 = x120 * x96
    x122 = x0 * (x101 + x92)
    x123 = x89 * x94
    x124 = x122 + x123
    x125 = 52.2557811793745 * x96
    x126 = x125 * x44
    x127 = 3.0 * x77
    x128 = 2.0 * x89
    x129 = x101 * x128 + x127
    x130 = x0 * (x129 + x76) + x102 * x89
    x131 = x103 * x130
    x132 = x125 * x64
    x133 = 0.179587122125167 * x53
    x134 = x133 * x64
    x135 = 90.5096679918781 * x94
    x136 = x106 * x96
    x137 = 90.5096679918781 * x111
    x138 = x44 * x88
    x139 = 52.2557811793745 * x88
    x140 = x139 * x89
    x141 = x106**2 * x82
    x142 = x141 + x84
    x143 = 30.169889330626 * x142
    x144 = x143 * x88
    x145 = 52.2557811793745 * x142
    x146 = x0 * (x109 + x115)
    x147 = x106 * x111
    x148 = x146 + x147
    x149 = 52.2557811793745 * x148
    x150 = x149 * x88
    x151 = 3.0 * x84
    x152 = 2.0 * x106
    x153 = x115 * x152 + x151
    x154 = x0 * (x153 + x83) + x106 * x116
    x155 = 30.169889330626 * x88
    x156 = x154 * x155
    x157 = x119 * x89 + x31 * x92
    x158 = 13.4923846833851 * x157
    x159 = x0 * (x118 + x129) + x124 * x89
    x160 = 23.3694957868871 * x42
    x161 = x160 * x96
    x162 = 4.0 * x77
    x163 = x0 * (2.0 * x100 + 2.0 * x122 + 2.0 * x123 + x162 * x60) + x130 * x89
    x164 = x55 * x80
    x165 = x12 * x133
    x166 = x106 * x125
    x167 = x133 * x42
    x168 = x106 * x142 + x109 * x31
    x169 = 13.4923846833851 * x168
    x170 = x160 * x88
    x171 = x0 * (x141 + x153) + x106 * x148
    x172 = x12 * x88
    x173 = 4.0 * x84
    x174 = (
        13.4923846833851 * x0 * (2.0 * x114 + 2.0 * x146 + 2.0 * x147 + x173 * x70)
        + 13.4923846833851 * x106 * x154
    )
    x175 = -x59 - C[1]
    x176 = x175**2 * x75
    x177 = x176 + x77
    x178 = x0 * (x15 + x72)
    x179 = x11 + x62
    x180 = x179 * x25
    x181 = x15 * x39 + x35
    x182 = x0 * (x181 + x28) + x25 * x32
    x183 = x0 * (x15 * x21 + 2.0 * x178 + 2.0 * x180 + 2.0 * x30) + x182 * x25
    x184 = x177 * x60
    x185 = x175 * x75
    x186 = x185 * x31
    x187 = x184 + x186
    x188 = 23.3694957868871 * x187
    x189 = x25**2 * x9
    x190 = x178 + x180
    x191 = x0 * (x181 + x189) + x190 * x25
    x192 = x191 * x96
    x193 = x101 * x175
    x194 = 2.0 * x193
    x195 = x127 + x176
    x196 = x0 * (x194 + x195)
    x197 = x187 * x60
    x198 = x196 + x197
    x199 = x11 + x189
    x200 = x199 * x25 + x31 * x72
    x201 = x133 * x177
    x202 = x177 * x89
    x203 = x186 + x202
    x204 = 30.169889330626 * x203
    x205 = x187 * x89
    x206 = x196 + x205
    x207 = x125 * x190
    x208 = x193 + x77
    x209 = x208 * x60
    x210 = x0 * (x101 + x185)
    x211 = x162 * x175
    x212 = 2.0 * x210 + x211
    x213 = x0 * (2.0 * x184 + 2.0 * x209 + x212)
    x214 = x198 * x89
    x215 = x213 + x214
    x216 = x103 * x199
    x217 = x133 * x199
    x218 = x0 * (x128 * x185 + x195) + x203 * x89
    x219 = x103 * x218
    x220 = x208 * x89
    x221 = 2.0 * x220
    x222 = x0 * (x184 + x202 + x212 + x221)
    x223 = x206 * x89
    x224 = x222 + x223
    x225 = x125 * x179
    x226 = 30.169889330626 * x25
    x227 = x0 * (x127 + x194 + x76)
    x228 = x89 * (x209 + x210)
    x229 = 2.0 * x205
    x230 = x0 * (3.0 * x196 + x197 + 2.0 * x227 + 2.0 * x228 + x229) + x215 * x89
    x231 = x57 * x6
    x232 = x231 * x54
    x233 = x230 * x232
    x234 = x232 * x25
    x235 = x6 * x87
    x236 = x218 * x235
    x237 = 90.5096679918781 * x179
    x238 = x133 * x179
    x239 = x206 * x235
    x240 = 52.2557811793745 * x235
    x241 = x240 * x25
    x242 = 52.2557811793745 * x238
    x243 = x143 * x235
    x244 = x235 * x25
    x245 = x154 * x235
    x246 = x175 * x92
    x247 = (
        x0 * (x128 * (x246 + x77) + 2.0 * x202 + x211 + x31 * (x185 + x92)) + x218 * x89
    )
    x248 = (
        x0
        * (
            x128 * (x210 + x220)
            + 2.0 * x196
            + x218
            + x229
            + x31 * (x127 + x193 + x246 + x93)
        )
        + x224 * x89
    )
    x249 = 23.3694957868871 * x14
    x250 = x54 * x57
    x251 = x14 * x250
    x252 = x14 * x87
    x253 = x133 * x29
    x254 = x249 * x87
    x255 = -x69 - C[2]
    x256 = x255**2 * x82
    x257 = x256 + x84
    x258 = 13.4923846833851 * x88
    x259 = x191 * x88
    x260 = x257 * x70
    x261 = x255 * x82
    x262 = x261 * x31
    x263 = x260 + x262
    x264 = 23.3694957868871 * x263
    x265 = x133 * x257
    x266 = x115 * x255
    x267 = 2.0 * x266
    x268 = x151 + x256
    x269 = x0 * (x267 + x268)
    x270 = x263 * x70
    x271 = x269 + x270
    x272 = x155 * x182
    x273 = x139 * x190
    x274 = 30.169889330626 * x217
    x275 = x155 * x199
    x276 = x106 * x257
    x277 = x262 + x276
    x278 = x106 * x263
    x279 = x269 + x278
    x280 = x266 + x84
    x281 = x280 * x70
    x282 = x0 * (x115 + x261)
    x283 = x173 * x255
    x284 = 2.0 * x282 + x283
    x285 = x0 * (2.0 * x260 + 2.0 * x281 + x284)
    x286 = x106 * x271
    x287 = x285 + x286
    x288 = x235 * x257
    x289 = x124 * x240
    x290 = x235 * x271
    x291 = 52.2557811793745 * x89
    x292 = x53 * x6
    x293 = x292 * x58
    x294 = x25 * x293
    x295 = x0 * (x152 * x261 + x268) + x106 * x277
    x296 = x155 * x295
    x297 = x139 * x179
    x298 = x106 * x280
    x299 = 2.0 * x298
    x300 = x0 * (x260 + x276 + x284 + x299)
    x301 = x106 * x279
    x302 = x300 + x301
    x303 = 30.169889330626 * x235
    x304 = x295 * x303
    x305 = x0 * (x151 + x267 + x83)
    x306 = x106 * (x281 + x282)
    x307 = 2.0 * x278
    x308 = x0 * (3.0 * x269 + x270 + 2.0 * x305 + 2.0 * x306 + x307) + x106 * x287
    x309 = x293 * x308
    x310 = 52.2557811793745 * x252
    x311 = x53 * x58
    x312 = x14 * x311
    x313 = x109 * x255
    x314 = (
        x0 * (x152 * (x313 + x84) + 2.0 * x276 + x283 + x31 * (x109 + x261)) + x106 * x295
    )
    x315 = (
        x0
        * (
            x152 * (x282 + x298)
            + 2.0 * x269
            + x295
            + x307
            + x31 * (x110 + x151 + x266 + x313)
        )
        + x106 * x302
    )

    # 180 item(s)
    result[0, 0, 0] = numpy.sum(
        x55
        * x58
        * (
            x0
            * (
                2.0 * x24
                + x31 * (x18 + 3.0 * x20 + x27 + x32)
                + 2.0 * x34
                + x39 * (x37 + x38)
                + 2.0 * x45
                + 2.0 * x48
            )
            + x25 * x51
        )
    )
    result[0, 0, 1] = numpy.sum(x61 * x68)
    result[0, 0, 2] = numpy.sum(x68 * x71)
    result[0, 0, 3] = numpy.sum(x74 * x81)
    result[0, 0, 4] = numpy.sum(x58 * x61 * x70 * x74)
    result[0, 0, 5] = numpy.sum(x73 * x86 * x88)
    result[0, 1, 0] = numpy.sum(x90 * x91)
    result[0, 1, 1] = numpy.sum(x66 * x97)
    result[0, 1, 2] = numpy.sum(x89 * x98 * x99)
    result[0, 1, 3] = numpy.sum(x102 * x104)
    result[0, 1, 4] = numpy.sum(x65 * x70 * x97)
    result[0, 1, 5] = numpy.sum(x105 * x85 * x90)
    result[0, 2, 0] = numpy.sum(x107 * x91)
    result[0, 2, 1] = numpy.sum(x108 * x60 * x99)
    result[0, 2, 2] = numpy.sum(x113 * x66)
    result[0, 2, 3] = numpy.sum(x104 * x106 * x78)
    result[0, 2, 4] = numpy.sum(x113 * x60 * x65)
    result[0, 2, 5] = numpy.sum(x105 * x117)
    result[0, 3, 0] = numpy.sum(x121 * x50)
    result[0, 3, 1] = numpy.sum(x124 * x126)
    result[0, 3, 2] = numpy.sum(x119 * x126 * x70)
    result[0, 3, 3] = numpy.sum(x131 * x64)
    result[0, 3, 4] = numpy.sum(x124 * x132 * x70)
    result[0, 3, 5] = numpy.sum(x120 * x134 * x85)
    result[0, 4, 0] = numpy.sum(x108 * x50 * x67 * x89)
    result[0, 4, 1] = numpy.sum(x135 * x136 * x44)
    result[0, 4, 2] = numpy.sum(x137 * x138 * x89)
    result[0, 4, 3] = numpy.sum(x102 * x106 * x132)
    result[0, 4, 4] = numpy.sum(x134 * x137 * x94)
    result[0, 4, 5] = numpy.sum(x116 * x140 * x64)
    result[0, 5, 0] = numpy.sum(x144 * x50)
    result[0, 5, 1] = numpy.sum(x138 * x145 * x60)
    result[0, 5, 2] = numpy.sum(x138 * x149)
    result[0, 5, 3] = numpy.sum(x134 * x143 * x78)
    result[0, 5, 4] = numpy.sum(x150 * x60 * x64)
    result[0, 5, 5] = numpy.sum(x156 * x64)
    result[0, 6, 0] = numpy.sum(x158 * x47 * x96)
    result[0, 6, 1] = numpy.sum(x159 * x161)
    result[0, 6, 2] = numpy.sum(x157 * x161 * x70)
    result[0, 6, 3] = numpy.sum(x12 * x163 * x164)
    result[0, 6, 4] = numpy.sum(x12 * x159 * x71 * x96)
    result[0, 6, 5] = numpy.sum(x157 * x165 * x86)
    result[0, 7, 0] = numpy.sum(x106 * x121 * x47)
    result[0, 7, 1] = numpy.sum(x124 * x166 * x42)
    result[0, 7, 2] = numpy.sum(x112 * x119 * x167)
    result[0, 7, 3] = numpy.sum(x106 * x12 * x131)
    result[0, 7, 4] = numpy.sum(x112 * x124 * x165)
    result[0, 7, 5] = numpy.sum(x116 * x120 * x165)
    result[0, 8, 0] = numpy.sum(x144 * x47 * x89)
    result[0, 8, 1] = numpy.sum(x142 * x167 * x95)
    result[0, 8, 2] = numpy.sum(x150 * x42 * x89)
    result[0, 8, 3] = numpy.sum(x102 * x143 * x165)
    result[0, 8, 4] = numpy.sum(x148 * x165 * x95)
    result[0, 8, 5] = numpy.sum(x12 * x156 * x89)
    result[0, 9, 0] = numpy.sum(x169 * x47 * x88)
    result[0, 9, 1] = numpy.sum(x168 * x170 * x60)
    result[0, 9, 2] = numpy.sum(x170 * x171)
    result[0, 9, 3] = numpy.sum(x165 * x169 * x78)
    result[0, 9, 4] = numpy.sum(x171 * x172 * x61)
    result[0, 9, 5] = numpy.sum(x172 * x174)
    result[1, 0, 0] = numpy.sum(x164 * x177 * x183)
    result[1, 0, 1] = numpy.sum(x188 * x192)
    result[1, 0, 2] = numpy.sum(x177 * x192 * x71)
    result[1, 0, 3] = numpy.sum(x164 * x198 * x200)
    result[1, 0, 4] = numpy.sum(x188 * x200 * x70 * x96)
    result[1, 0, 5] = numpy.sum(x200 * x201 * x86)
    result[1, 1, 0] = numpy.sum(x182 * x204 * x96)
    result[1, 1, 1] = numpy.sum(x206 * x207)
    result[1, 1, 2] = numpy.sum(x203 * x207 * x70)
    result[1, 1, 3] = numpy.sum(x215 * x216)
    result[1, 1, 4] = numpy.sum(x125 * x199 * x206 * x70)
    result[1, 1, 5] = numpy.sum(x204 * x217 * x85)
    result[1, 2, 0] = numpy.sum(x103 * x106 * x177 * x182)
    result[1, 2, 1] = numpy.sum(x106 * x187 * x207)
    result[1, 2, 2] = numpy.sum(x112 * x190 * x201)
    result[1, 2, 3] = numpy.sum(x106 * x198 * x216)
    result[1, 2, 4] = numpy.sum(x112 * x187 * x217)
    result[1, 2, 5] = numpy.sum(x117 * x199 * x201)
    result[1, 3, 0] = numpy.sum(x219 * x32)
    result[1, 3, 1] = numpy.sum(x224 * x225)
    result[1, 3, 2] = numpy.sum(x218 * x225 * x70)
    result[1, 3, 3] = numpy.sum(x226 * x233)
    result[1, 3, 4] = numpy.sum(x224 * x234 * x98)
    result[1, 3, 5] = numpy.sum(x226 * x236 * x85)
    result[1, 4, 0] = numpy.sum(x166 * x203 * x32)
    result[1, 4, 1] = numpy.sum(x136 * x206 * x237)
    result[1, 4, 2] = numpy.sum(x137 * x203 * x238)
    result[1, 4, 3] = numpy.sum(x108 * x215 * x234)
    result[1, 4, 4] = numpy.sum(x137 * x239 * x25)
    result[1, 4, 5] = numpy.sum(x116 * x203 * x241)
    result[1, 5, 0] = numpy.sum(x143 * x201 * x32)
    result[1, 5, 1] = numpy.sum(x142 * x187 * x242)
    result[1, 5, 2] = numpy.sum(x149 * x179 * x201)
    result[1, 5, 3] = numpy.sum(x198 * x243 * x25)
    result[1, 5, 4] = numpy.sum(x149 * x187 * x244)
    result[1, 5, 5] = numpy.sum(x177 * x226 * x245)
    result[1, 6, 0] = numpy.sum(x164 * x247 * x29)
    result[1, 6, 1] = numpy.sum(x248 * x249 * x250)
    result[1, 6, 2] = numpy.sum(x247 * x251 * x71)
    result[1, 6, 3] = numpy.sum(
        x231
        * x55
        * (
            x0
            * (
                x128 * (x227 + x228)
                + 2.0 * x213
                + 2.0 * x214
                + 2.0 * x222
                + 2.0 * x223
                + x31 * (x102 + x209 + 3.0 * x210 + x221)
            )
            + x230 * x89
        )
    )
    result[1, 6, 4] = numpy.sum(x232 * x248 * x71)
    result[1, 6, 5] = numpy.sum(x235 * x247 * x86)
    result[1, 7, 0] = numpy.sum(x106 * x219 * x29)
    result[1, 7, 1] = numpy.sum(x108 * x224 * x251)
    result[1, 7, 2] = numpy.sum(x112 * x218 * x252)
    result[1, 7, 3] = numpy.sum(x107 * x233)
    result[1, 7, 4] = numpy.sum(x112 * x224 * x235)
    result[1, 7, 5] = numpy.sum(x117 * x236)
    result[1, 8, 0] = numpy.sum(x143 * x203 * x253)
    result[1, 8, 1] = numpy.sum(x145 * x206 * x252)
    result[1, 8, 2] = numpy.sum(x149 * x203 * x252)
    result[1, 8, 3] = numpy.sum(x215 * x243)
    result[1, 8, 4] = numpy.sum(x149 * x239)
    result[1, 8, 5] = numpy.sum(x204 * x245)
    result[1, 9, 0] = numpy.sum(x169 * x201 * x29)
    result[1, 9, 1] = numpy.sum(x168 * x188 * x252)
    result[1, 9, 2] = numpy.sum(x171 * x177 * x254)
    result[1, 9, 3] = numpy.sum(x169 * x198 * x235)
    result[1, 9, 4] = numpy.sum(x171 * x188 * x235)
    result[1, 9, 5] = numpy.sum(x174 * x177 * x235)
    result[2, 0, 0] = numpy.sum(x183 * x257 * x258)
    result[2, 0, 1] = numpy.sum(x257 * x259 * x61)
    result[2, 0, 2] = numpy.sum(x259 * x264)
    result[2, 0, 3] = numpy.sum(x200 * x265 * x79)
    result[2, 0, 4] = numpy.sum(x200 * x264 * x60 * x88)
    result[2, 0, 5] = numpy.sum(x200 * x258 * x271)
    result[2, 1, 0] = numpy.sum(x257 * x272 * x89)
    result[2, 1, 1] = numpy.sum(x190 * x265 * x95)
    result[2, 1, 2] = numpy.sum(x263 * x273 * x89)
    result[2, 1, 3] = numpy.sum(x102 * x257 * x274)
    result[2, 1, 4] = numpy.sum(x217 * x263 * x95)
    result[2, 1, 5] = numpy.sum(x271 * x275 * x89)
    result[2, 2, 0] = numpy.sum(x272 * x277)
    result[2, 2, 1] = numpy.sum(x273 * x277 * x60)
    result[2, 2, 2] = numpy.sum(x273 * x279)
    result[2, 2, 3] = numpy.sum(x274 * x277 * x78)
    result[2, 2, 4] = numpy.sum(x139 * x199 * x279 * x60)
    result[2, 2, 5] = numpy.sum(x275 * x287)
    result[2, 3, 0] = numpy.sum(x120 * x265 * x32)
    result[2, 3, 1] = numpy.sum(x124 * x242 * x257)
    result[2, 3, 2] = numpy.sum(x119 * x242 * x263)
    result[2, 3, 3] = numpy.sum(x130 * x226 * x288)
    result[2, 3, 4] = numpy.sum(x25 * x263 * x289)
    result[2, 3, 5] = numpy.sum(x120 * x25 * x290)
    result[2, 4, 0] = numpy.sum(x140 * x277 * x32)
    result[2, 4, 1] = numpy.sum(x135 * x238 * x277)
    result[2, 4, 2] = numpy.sum(x237 * x279 * x88 * x89)
    result[2, 4, 3] = numpy.sum(x102 * x241 * x277)
    result[2, 4, 4] = numpy.sum(x135 * x244 * x279)
    result[2, 4, 5] = numpy.sum(x287 * x291 * x294)
    result[2, 5, 0] = numpy.sum(x296 * x32)
    result[2, 5, 1] = numpy.sum(x295 * x297 * x60)
    result[2, 5, 2] = numpy.sum(x297 * x302)
    result[2, 5, 3] = numpy.sum(x25 * x304 * x78)
    result[2, 5, 4] = numpy.sum(52.2557811793745 * x294 * x302 * x60)
    result[2, 5, 5] = numpy.sum(x226 * x309)
    result[2, 6, 0] = numpy.sum(x158 * x253 * x257)
    result[2, 6, 1] = numpy.sum(x159 * x254 * x257)
    result[2, 6, 2] = numpy.sum(x157 * x252 * x264)
    result[2, 6, 3] = numpy.sum(13.4923846833851 * x163 * x288)
    result[2, 6, 4] = numpy.sum(x159 * x235 * x264)
    result[2, 6, 5] = numpy.sum(x158 * x290)
    result[2, 7, 0] = numpy.sum(x120 * x253 * x277)
    result[2, 7, 1] = numpy.sum(x124 * x277 * x310)
    result[2, 7, 2] = numpy.sum(x119 * x279 * x310)
    result[2, 7, 3] = numpy.sum(x130 * x277 * x303)
    result[2, 7, 4] = numpy.sum(x279 * x289)
    result[2, 7, 5] = numpy.sum(x120 * x235 * x287)
    result[2, 8, 0] = numpy.sum(x29 * x296 * x89)
    result[2, 8, 1] = numpy.sum(x252 * x295 * x95)
    result[2, 8, 2] = numpy.sum(x291 * x302 * x312)
    result[2, 8, 3] = numpy.sum(x102 * x304)
    result[2, 8, 4] = numpy.sum(x235 * x302 * x95)
    result[2, 8, 5] = numpy.sum(x309 * x90)
    result[2, 9, 0] = numpy.sum(x258 * x29 * x314)
    result[2, 9, 1] = numpy.sum(x312 * x314 * x61)
    result[2, 9, 2] = numpy.sum(x249 * x311 * x315)
    result[2, 9, 3] = numpy.sum(x292 * x314 * x81)
    result[2, 9, 4] = numpy.sum(x293 * x315 * x61)
    result[2, 9, 5] = numpy.sum(
        13.4923846833851
        * x293
        * (
            x0
            * (
                x152 * (x305 + x306)
                + 2.0 * x285
                + 2.0 * x286
                + 2.0 * x300
                + 2.0 * x301
                + x31 * (x116 + x281 + 3.0 * x282 + x299)
            )
            + x106 * x308
        )
    )
    return result


def diag_quadrupole3d_33(ax, da, A, bx, db, B, C):
    """Cartesian 3D (ff) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 10, 10), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x2 - C[0]
    x11 = x10 * x8
    x12 = x0 * (x11 + x9)
    x13 = x0 * x8
    x14 = x10 * x9
    x15 = x13 + x14
    x16 = x15 * x3
    x17 = x12 + x16
    x18 = x17 * x3
    x19 = x10**2 * x8
    x20 = x13 + x19
    x21 = x20 * x3
    x22 = 2.0 * x0
    x23 = x11 * x22
    x24 = x21 + x23
    x25 = x24 * x3
    x26 = x3**2 * x8
    x27 = 3.0 * x13
    x28 = 2.0 * x14 + x27
    x29 = x0 * (x26 + x28)
    x30 = x0 * (x19 + x28)
    x31 = 2.0 * x29 + 3.0 * x30
    x32 = x0 * (2.0 * x18 + 3.0 * x25 + x31)
    x33 = -x2 - A[0]
    x34 = x17 * x33
    x35 = x0 * (3.0 * x26 + x27)
    x36 = x13 + x26
    x37 = x3 * x36
    x38 = x22 * x9
    x39 = x37 + x38
    x40 = x33 * x39
    x41 = x35 + x40
    x42 = 3.0 * x12
    x43 = x0 * (3.0 * x16 + x39 + x42)
    x44 = x33 * (x18 + x29)
    x45 = 2.0 * x33
    x46 = 4.0 * x0
    x47 = x11 * x46
    x48 = 2.0 * x12 + x47
    x49 = x0 * (2.0 * x16 + 2.0 * x21 + x48)
    x50 = x25 + x30
    x51 = x3 * x50
    x52 = x49 + x51
    x53 = x33 * x52
    x54 = x24 * x33
    x55 = 2.0 * x54
    x56 = x0 * (x25 + x31 + 2.0 * x34 + x55)
    x57 = x33 * x50
    x58 = x49 + x57
    x59 = x33 * x58
    x60 = x32 + x53
    x61 = x0 * (2.0 * x43 + 2.0 * x44 + 4.0 * x49 + x51 + 3.0 * x57) + x33 * x60
    x62 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x63 = da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**4.5)
    x64 = x62 * x63
    x65 = 12.0679557322504 * x64
    x66 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x67 = 0.564189583547756 * x1
    x68 = x66 * x67
    x69 = -x1 * (ax * A[1] + bx * B[1])
    x70 = -x69 - B[1]
    x71 = x68 * x70
    x72 = x15 * x33
    x73 = 2.0 * x72
    x74 = x33 * x36
    x75 = x38 + x74
    x76 = x20 * x33
    x77 = x0 * (x21 + x48 + x73 + x76)
    x78 = x30 + x54
    x79 = x33 * x78
    x80 = x56 + x59
    x81 = 26.9847693667702 * x64
    x82 = x81 * (
        x0
        * (
            x22 * (x16 + x42 + x73 + x75)
            + x45 * (x29 + x34)
            + 2.0 * x49
            + 2.0 * x57
            + 2.0 * x77
            + 2.0 * x79
        )
        + x33 * x80
    )
    x83 = -x1 * (ax * A[2] + bx * B[2])
    x84 = -x83 - B[2]
    x85 = x68 * x84
    x86 = x33 * x9
    x87 = x11 * x33
    x88 = x23 + x76
    x89 = x0 * (x11 * x45 + x19 + x27) + x33 * x88
    x90 = x77 + x79
    x91 = (
        x0 * (x22 * (x14 + x27 + x86 + x87) + 2.0 * x30 + x45 * (x12 + x72) + x55 + x89)
        + x33 * x90
    )
    x92 = x64 * x91
    x93 = x66 * x7
    x94 = x70**2 * x93
    x95 = x0 * x93
    x96 = x94 + x95
    x97 = 26.9847693667702 * x96
    x98 = 0.318309886183791 * x6
    x99 = x97 * x98
    x100 = x62 * x7
    x101 = x100 * x84**2
    x102 = x0 * x100
    x103 = x101 + x102
    x104 = 26.9847693667702 * x103
    x105 = x63 * x98
    x106 = x105 * x66
    x107 = x104 * x106
    x108 = x33 * x8
    x109 = x0 * (x22 * (x108 + x11) + x45 * (x13 + x87) + x47 + 2.0 * x76) + x33 * x89
    x110 = x70 * x96
    x111 = x70 * x93
    x112 = x111 * x22
    x113 = x110 + x112
    x114 = 12.0679557322504 * x113
    x115 = x64 * x98
    x116 = x64 * x84
    x117 = x103 * x84
    x118 = x100 * x84
    x119 = x118 * x22
    x120 = x117 + x119
    x121 = 12.0679557322504 * x120
    x122 = -x69 - A[1]
    x123 = x122 * x68
    x124 = x61 * x81
    x125 = 60.3397786612521 * x80
    x126 = x122 * x93
    x127 = x126 * x70
    x128 = x127 + x95
    x129 = x115 * x128
    x130 = x125 * x64
    x131 = x122 * x96
    x132 = x112 + x131
    x133 = 60.3397786612521 * x132
    x134 = x115 * x90
    x135 = 104.511562358749 * x90
    x136 = 60.3397786612521 * x103
    x137 = x106 * x90
    x138 = 3.0 * x95
    x139 = x0 * (x138 + 3.0 * x94)
    x140 = x113 * x122
    x141 = x139 + x140
    x142 = 26.9847693667702 * x141
    x143 = x115 * x89
    x144 = 60.3397786612521 * x89
    x145 = x115 * x84
    x146 = 0.179587122125167 * x63
    x147 = x128 * x146
    x148 = 26.9847693667702 * x106
    x149 = x148 * x89
    x150 = -x83 - A[2]
    x151 = x100 * x150
    x152 = x151 * x84
    x153 = x102 + x152
    x154 = x106 * x153
    x155 = 60.3397786612521 * x96
    x156 = x103 * x150
    x157 = x119 + x156
    x158 = 60.3397786612521 * x157
    x159 = 26.9847693667702 * x113
    x160 = x146 * x153
    x161 = x106 * x70
    x162 = 3.0 * x102
    x163 = x0 * (3.0 * x101 + x162)
    x164 = x120 * x150
    x165 = x163 + x164
    x166 = x122**2 * x93
    x167 = x166 + x95
    x168 = 26.9847693667702 * x167
    x169 = x115 * x168
    x170 = x0 * (x111 + x126)
    x171 = x122 * x128
    x172 = x170 + x171
    x173 = 60.3397786612521 * x115
    x174 = x172 * x173
    x175 = 60.3397786612521 * x167
    x176 = 2.0 * x122
    x177 = x111 * x176 + x138
    x178 = x0 * (x177 + x94)
    x179 = x122 * x132
    x180 = x178 + x179
    x181 = x173 * x180
    x182 = 104.511562358749 * x78
    x183 = x136 * x146
    x184 = x70 * x95
    x185 = x0 * (x110 + 3.0 * x131 + 8.0 * x184) + x122 * x141
    x186 = 26.9847693667702 * x115
    x187 = x185 * x186
    x188 = x146 * x88
    x189 = x150 * x64
    x190 = 46.7389915737742 * x189
    x191 = 104.511562358749 * x58
    x192 = x115 * x150
    x193 = 181.019335983756 * x147
    x194 = 104.511562358749 * x157
    x195 = x106 * x122
    x196 = 46.7389915737742 * x88
    x197 = 104.511562358749 * x88
    x198 = x100 * x150**2
    x199 = x102 + x198
    x200 = x148 * x199
    x201 = 60.3397786612521 * x106
    x202 = x201 * x58
    x203 = x0 * (x118 + x151)
    x204 = x150 * x153
    x205 = x203 + x204
    x206 = x146 * x199
    x207 = 104.511562358749 * x205
    x208 = 2.0 * x150
    x209 = x118 * x208 + x162
    x210 = x0 * (x101 + x209)
    x211 = x150 * x157
    x212 = x210 + x211
    x213 = x201 * x212
    x214 = x102 * x84
    x215 = x0 * (x117 + 3.0 * x156 + 8.0 * x214) + x150 * x165
    x216 = x148 * x215
    x217 = x122 * x167 + x126 * x22
    x218 = x65 * x98
    x219 = x0 * (x166 + x177) + x122 * x172
    x220 = x186 * x50
    x221 = 2.0 * x0 * (x131 + x170 + x171 + 2.0 * x184) + x122 * x180
    x222 = x186 * x221
    x223 = 46.7389915737742 * x24
    x224 = x146 * x217
    x225 = x0 * (2.0 * x139 + 2.0 * x140 + 3.0 * x178 + 3.0 * x179) + x122 * x185
    x226 = x146 * x20
    x227 = 104.511562358749 * x24
    x228 = x146 * x24
    x229 = 60.3397786612521 * x20
    x230 = 60.3397786612521 * x226
    x231 = 60.3397786612521 * x147
    x232 = 60.3397786612521 * x206
    x233 = 26.9847693667702 * x206
    x234 = x150 * x199 + x151 * x22
    x235 = 12.0679557322504 * x234
    x236 = x148 * x50
    x237 = x0 * (x198 + x209) + x150 * x205
    x238 = 2.0 * x0 * (x156 + x203 + x204 + 2.0 * x214) + x150 * x212
    x239 = x148 * x238
    x240 = (
        12.0679557322504 * x0 * (2.0 * x163 + 2.0 * x164 + 3.0 * x210 + 3.0 * x211)
        + 12.0679557322504 * x150 * x215
    )
    x241 = -x69 - C[1]
    x242 = x241**2 * x93
    x243 = x242 + x95
    x244 = x27 + x45 * x9
    x245 = x0 * (x244 + x26)
    x246 = x33 * x75
    x247 = x0 * (8.0 * x0 * x9 + x37 + 3.0 * x74) + x33 * x41
    x248 = x0 * (3.0 * x245 + 3.0 * x246 + 2.0 * x35 + 2.0 * x40) + x247 * x33
    x249 = x243 * x70
    x250 = x241 * x93
    x251 = x22 * x250
    x252 = x249 + x251
    x253 = x0 * (x108 + x9)
    x254 = x13 + x86
    x255 = x254 * x33
    x256 = x245 + x246
    x257 = x0 * (2.0 * x253 + 2.0 * x255 + x46 * x9 + 2.0 * x74) + x256 * x33
    x258 = x186 * x257
    x259 = x33**2 * x8
    x260 = x253 + x255
    x261 = x0 * (x244 + x259) + x260 * x33
    x262 = x111 * x241
    x263 = x138 + 2.0 * x262
    x264 = x0 * (x242 + x263)
    x265 = x252 * x70
    x266 = x264 + x265
    x267 = x186 * x266
    x268 = 46.7389915737742 * x252
    x269 = x104 * x146
    x270 = x262 + x95
    x271 = x270 * x70
    x272 = x0 * (x111 + x250)
    x273 = 4.0 * x241 * x95
    x274 = 2.0 * x272 + x273
    x275 = x0 * (2.0 * x249 + 2.0 * x271 + x274)
    x276 = x266 * x70
    x277 = x275 + x276
    x278 = x13 + x259
    x279 = x108 * x22 + x278 * x33
    x280 = x146 * x243
    x281 = x122 * x243
    x282 = x251 + x281
    x283 = 26.9847693667702 * x282
    x284 = x122 * x252
    x285 = x264 + x284
    x286 = x173 * x256
    x287 = x122 * x266
    x288 = x275 + x287
    x289 = x173 * x288
    x290 = 104.511562358749 * x285
    x291 = x271 + x272
    x292 = x291 * x70
    x293 = x0 * (x263 + x94)
    x294 = 3.0 * x264 + 2.0 * x293
    x295 = x0 * (3.0 * x265 + 2.0 * x292 + x294)
    x296 = x122 * x277
    x297 = x295 + x296
    x298 = x186 * x278
    x299 = x146 * x278
    x300 = 60.3397786612521 * x160
    x301 = 104.511562358749 * x160
    x302 = 26.9847693667702 * x280
    x303 = x0 * (x138 + x176 * x250 + x242) + x122 * x282
    x304 = x186 * x303
    x305 = x122 * x270
    x306 = 2.0 * x305
    x307 = x0 * (x249 + x274 + x281 + x306)
    x308 = x122 * x285
    x309 = x307 + x308
    x310 = x173 * x75
    x311 = x122 * x291
    x312 = 2.0 * x284
    x313 = x0 * (x265 + x294 + 2.0 * x311 + x312)
    x314 = x122 * x288
    x315 = x313 + x314
    x316 = 104.511562358749 * x254
    x317 = 3.0 * x272
    x318 = x0 * (x113 + 3.0 * x271 + x317)
    x319 = x122 * (x292 + x293)
    x320 = x0 * (4.0 * x275 + x276 + 3.0 * x287 + 2.0 * x318 + 2.0 * x319) + x122 * x297
    x321 = x5 * x67
    x322 = x321 * x81
    x323 = x320 * x322
    x324 = x321 * x33
    x325 = 60.3397786612521 * x315
    x326 = x105 * x5
    x327 = x309 * x326
    x328 = 26.9847693667702 * x326
    x329 = x303 * x328
    x330 = 46.7389915737742 * x282
    x331 = x146 * x282
    x332 = x326 * x33
    x333 = 104.511562358749 * x153
    x334 = 60.3397786612521 * x280
    x335 = x146 * x252
    x336 = x199 * x328
    x337 = 60.3397786612521 * x326
    x338 = x205 * x337
    x339 = x212 * x337
    x340 = 12.0679557322504 * x39
    x341 = x126 * x241
    x342 = (
        x0 * (x176 * (x341 + x95) + x22 * (x126 + x250) + x273 + 2.0 * x281) + x122 * x303
    )
    x343 = x115 * x342
    x344 = (
        x0
        * (
            x176 * (x272 + x305)
            + x22 * (x127 + x138 + x262 + x341)
            + 2.0 * x264
            + x303
            + x312
        )
        + x122 * x309
    )
    x345 = 26.9847693667702 * x36
    x346 = x322 * (
        x0
        * (
            x176 * (x293 + x311)
            + x22 * (x132 + x271 + x306 + x317)
            + 2.0 * x275
            + 2.0 * x287
            + 2.0 * x307
            + 2.0 * x308
        )
        + x122 * x315
    )
    x347 = x3 * x321
    x348 = x104 * x326
    x349 = 60.3397786612521 * x36
    x350 = x3 * x326
    x351 = x266 * x328
    x352 = x238 * x328
    x353 = -x83 - C[2]
    x354 = x100 * x353**2
    x355 = x102 + x354
    x356 = 12.0679557322504 * x106
    x357 = x148 * x257
    x358 = x355 * x84
    x359 = x100 * x353
    x360 = x22 * x359
    x361 = x358 + x360
    x362 = x146 * x355
    x363 = 46.7389915737742 * x361
    x364 = x118 * x353
    x365 = x162 + 2.0 * x364
    x366 = x0 * (x354 + x365)
    x367 = x361 * x84
    x368 = x366 + x367
    x369 = x148 * x368
    x370 = x146 * x361
    x371 = x102 + x364
    x372 = x371 * x84
    x373 = x0 * (x118 + x359)
    x374 = 4.0 * x102 * x353
    x375 = 2.0 * x373 + x374
    x376 = x0 * (2.0 * x358 + 2.0 * x372 + x375)
    x377 = x368 * x84
    x378 = x376 + x377
    x379 = x148 * x247
    x380 = x201 * x256
    x381 = 104.511562358749 * x147
    x382 = x201 * x260
    x383 = x148 * x278
    x384 = x150 * x355
    x385 = x360 + x384
    x386 = x150 * x361
    x387 = x366 + x386
    x388 = x146 * x385
    x389 = 104.511562358749 * x387
    x390 = x150 * x368
    x391 = x376 + x390
    x392 = x372 + x373
    x393 = x392 * x84
    x394 = x0 * (x101 + x365)
    x395 = 3.0 * x366 + 2.0 * x394
    x396 = x0 * (3.0 * x367 + 2.0 * x393 + x395)
    x397 = x150 * x378
    x398 = x396 + x397
    x399 = 60.3397786612521 * x362
    x400 = x146 * x254
    x401 = x185 * x328
    x402 = x180 * x337
    x403 = x172 * x337
    x404 = x168 * x326
    x405 = 46.7389915737742 * x385
    x406 = 104.511562358749 * x128
    x407 = x5 * x63
    x408 = x407 * x68
    x409 = x33 * x408
    x410 = x0 * (x162 + x208 * x359 + x354) + x150 * x385
    x411 = x148 * x410
    x412 = x201 * x75
    x413 = x150 * x371
    x414 = 2.0 * x413
    x415 = x0 * (x358 + x375 + x384 + x414)
    x416 = x150 * x387
    x417 = x415 + x416
    x418 = x150 * x392
    x419 = 2.0 * x386
    x420 = x0 * (x367 + x395 + 2.0 * x418 + x419)
    x421 = x150 * x391
    x422 = x420 + x421
    x423 = x326 * x410
    x424 = x326 * x417
    x425 = 60.3397786612521 * x422
    x426 = 3.0 * x373
    x427 = x0 * (x120 + 3.0 * x372 + x426)
    x428 = x150 * (x393 + x394)
    x429 = x0 * (4.0 * x376 + x377 + 3.0 * x390 + 2.0 * x427 + 2.0 * x428) + x150 * x398
    x430 = 26.9847693667702 * x408
    x431 = x429 * x430
    x432 = x221 * x328
    x433 = x328 * x368
    x434 = 12.0679557322504 * x326
    x435 = x3 * x408
    x436 = x151 * x353
    x437 = (
        x0 * (x208 * (x102 + x436) + x22 * (x151 + x359) + x374 + 2.0 * x384)
        + x150 * x410
    )
    x438 = x106 * x437
    x439 = (
        x0
        * (
            x208 * (x373 + x413)
            + x22 * (x152 + x162 + x364 + x436)
            + 2.0 * x366
            + x410
            + x419
        )
        + x150 * x417
    )
    x440 = x407 * x99
    x441 = x430 * (
        x0
        * (
            x208 * (x394 + x418)
            + x22 * (x157 + x372 + x414 + x426)
            + 2.0 * x376
            + 2.0 * x390
            + 2.0 * x415
            + 2.0 * x416
        )
        + x150 * x422
    )

    # 300 item(s)
    result[0, 0, 0] = numpy.sum(
        x65
        * x68
        * (
            x0
            * (
                x22 * (x18 + 4.0 * x29 + 3.0 * x34 + x41)
                + 2.0 * x32
                + x45 * (x43 + x44)
                + 2.0 * x53
                + 3.0 * x56
                + 3.0 * x59
            )
            + x33 * x61
        )
    )
    result[0, 0, 1] = numpy.sum(x71 * x82)
    result[0, 0, 2] = numpy.sum(x82 * x85)
    result[0, 0, 3] = numpy.sum(x92 * x99)
    result[0, 0, 4] = numpy.sum(46.7389915737742 * x71 * x84 * x92)
    result[0, 0, 5] = numpy.sum(x107 * x91)
    result[0, 0, 6] = numpy.sum(x109 * x114 * x115)
    result[0, 0, 7] = numpy.sum(x109 * x116 * x99)
    result[0, 0, 8] = numpy.sum(x107 * x109 * x70)
    result[0, 0, 9] = numpy.sum(x106 * x109 * x121)
    result[0, 1, 0] = numpy.sum(x123 * x124)
    result[0, 1, 1] = numpy.sum(x125 * x129)
    result[0, 1, 2] = numpy.sum(x122 * x130 * x85)
    result[0, 1, 3] = numpy.sum(x133 * x134)
    result[0, 1, 4] = numpy.sum(x129 * x135 * x84)
    result[0, 1, 5] = numpy.sum(x122 * x136 * x137)
    result[0, 1, 6] = numpy.sum(x142 * x143)
    result[0, 1, 7] = numpy.sum(x132 * x144 * x145)
    result[0, 1, 8] = numpy.sum(x136 * x147 * x89)
    result[0, 1, 9] = numpy.sum(x120 * x122 * x149)
    result[0, 2, 0] = numpy.sum(x124 * x150 * x68)
    result[0, 2, 1] = numpy.sum(x130 * x150 * x71)
    result[0, 2, 2] = numpy.sum(x125 * x154)
    result[0, 2, 3] = numpy.sum(x134 * x150 * x155)
    result[0, 2, 4] = numpy.sum(x135 * x154 * x70)
    result[0, 2, 5] = numpy.sum(x137 * x158)
    result[0, 2, 6] = numpy.sum(x143 * x150 * x159)
    result[0, 2, 7] = numpy.sum(x144 * x160 * x96)
    result[0, 2, 8] = numpy.sum(x144 * x157 * x161)
    result[0, 2, 9] = numpy.sum(x149 * x165)
    result[0, 3, 0] = numpy.sum(x169 * x60)
    result[0, 3, 1] = numpy.sum(x174 * x58)
    result[0, 3, 2] = numpy.sum(x145 * x175 * x58)
    result[0, 3, 3] = numpy.sum(x181 * x78)
    result[0, 3, 4] = numpy.sum(x145 * x172 * x182)
    result[0, 3, 5] = numpy.sum(x167 * x183 * x78)
    result[0, 3, 6] = numpy.sum(x187 * x88)
    result[0, 3, 7] = numpy.sum(x181 * x84 * x88)
    result[0, 3, 8] = numpy.sum(x172 * x183 * x88)
    result[0, 3, 9] = numpy.sum(x120 * x168 * x188)
    result[0, 4, 0] = numpy.sum(x123 * x190 * x60)
    result[0, 4, 1] = numpy.sum(x129 * x150 * x191)
    result[0, 4, 2] = numpy.sum(x122 * x154 * x191)
    result[0, 4, 3] = numpy.sum(x132 * x182 * x192)
    result[0, 4, 4] = numpy.sum(x153 * x193 * x78)
    result[0, 4, 5] = numpy.sum(x194 * x195 * x78)
    result[0, 4, 6] = numpy.sum(x141 * x192 * x196)
    result[0, 4, 7] = numpy.sum(x132 * x160 * x197)
    result[0, 4, 8] = numpy.sum(x147 * x157 * x197)
    result[0, 4, 9] = numpy.sum(x165 * x195 * x196)
    result[0, 5, 0] = numpy.sum(x200 * x60)
    result[0, 5, 1] = numpy.sum(x199 * x202 * x70)
    result[0, 5, 2] = numpy.sum(x202 * x205)
    result[0, 5, 3] = numpy.sum(x155 * x206 * x78)
    result[0, 5, 4] = numpy.sum(x161 * x207 * x78)
    result[0, 5, 5] = numpy.sum(x213 * x78)
    result[0, 5, 6] = numpy.sum(x159 * x188 * x199)
    result[0, 5, 7] = numpy.sum(x155 * x188 * x205)
    result[0, 5, 8] = numpy.sum(x213 * x70 * x88)
    result[0, 5, 9] = numpy.sum(x216 * x88)
    result[0, 6, 0] = numpy.sum(x217 * x218 * x52)
    result[0, 6, 1] = numpy.sum(x219 * x220)
    result[0, 6, 2] = numpy.sum(x217 * x220 * x84)
    result[0, 6, 3] = numpy.sum(x222 * x24)
    result[0, 6, 4] = numpy.sum(x145 * x219 * x223)
    result[0, 6, 5] = numpy.sum(x104 * x224 * x24)
    result[0, 6, 6] = numpy.sum(x20 * x218 * x225)
    result[0, 6, 7] = numpy.sum(x20 * x222 * x84)
    result[0, 6, 8] = numpy.sum(x104 * x219 * x226)
    result[0, 6, 9] = numpy.sum(x121 * x20 * x224)
    result[0, 7, 0] = numpy.sum(x150 * x169 * x52)
    result[0, 7, 1] = numpy.sum(x150 * x174 * x50)
    result[0, 7, 2] = numpy.sum(x160 * x175 * x50)
    result[0, 7, 3] = numpy.sum(x150 * x181 * x24)
    result[0, 7, 4] = numpy.sum(x160 * x172 * x227)
    result[0, 7, 5] = numpy.sum(x157 * x175 * x228)
    result[0, 7, 6] = numpy.sum(x150 * x187 * x20)
    result[0, 7, 7] = numpy.sum(x160 * x180 * x229)
    result[0, 7, 8] = numpy.sum(x157 * x172 * x230)
    result[0, 7, 9] = numpy.sum(x165 * x168 * x226)
    result[0, 8, 0] = numpy.sum(x122 * x200 * x52)
    result[0, 8, 1] = numpy.sum(x199 * x231 * x50)
    result[0, 8, 2] = numpy.sum(x122 * x201 * x205 * x50)
    result[0, 8, 3] = numpy.sum(x132 * x232 * x24)
    result[0, 8, 4] = numpy.sum(x147 * x205 * x227)
    result[0, 8, 5] = numpy.sum(x122 * x213 * x24)
    result[0, 8, 6] = numpy.sum(x141 * x20 * x233)
    result[0, 8, 7] = numpy.sum(x132 * x205 * x230)
    result[0, 8, 8] = numpy.sum(x147 * x212 * x229)
    result[0, 8, 9] = numpy.sum(x122 * x20 * x216)
    result[0, 9, 0] = numpy.sum(x106 * x235 * x52)
    result[0, 9, 1] = numpy.sum(x234 * x236 * x70)
    result[0, 9, 2] = numpy.sum(x236 * x237)
    result[0, 9, 3] = numpy.sum(x228 * x234 * x97)
    result[0, 9, 4] = numpy.sum(x161 * x223 * x237)
    result[0, 9, 5] = numpy.sum(x239 * x24)
    result[0, 9, 6] = numpy.sum(x113 * x226 * x235)
    result[0, 9, 7] = numpy.sum(x226 * x237 * x97)
    result[0, 9, 8] = numpy.sum(x20 * x239 * x70)
    result[0, 9, 9] = numpy.sum(x106 * x20 * x240)
    result[1, 0, 0] = numpy.sum(x218 * x243 * x248)
    result[1, 0, 1] = numpy.sum(x252 * x258)
    result[1, 0, 2] = numpy.sum(x243 * x258 * x84)
    result[1, 0, 3] = numpy.sum(x261 * x267)
    result[1, 0, 4] = numpy.sum(x145 * x261 * x268)
    result[1, 0, 5] = numpy.sum(x243 * x261 * x269)
    result[1, 0, 6] = numpy.sum(x218 * x277 * x279)
    result[1, 0, 7] = numpy.sum(x267 * x279 * x84)
    result[1, 0, 8] = numpy.sum(x252 * x269 * x279)
    result[1, 0, 9] = numpy.sum(x121 * x279 * x280)
    result[1, 1, 0] = numpy.sum(x115 * x247 * x283)
    result[1, 1, 1] = numpy.sum(x285 * x286)
    result[1, 1, 2] = numpy.sum(x282 * x286 * x84)
    result[1, 1, 3] = numpy.sum(x260 * x289)
    result[1, 1, 4] = numpy.sum(x145 * x260 * x290)
    result[1, 1, 5] = numpy.sum(x183 * x260 * x282)
    result[1, 1, 6] = numpy.sum(x297 * x298)
    result[1, 1, 7] = numpy.sum(x278 * x289 * x84)
    result[1, 1, 8] = numpy.sum(x183 * x278 * x285)
    result[1, 1, 9] = numpy.sum(x120 * x283 * x299)
    result[1, 2, 0] = numpy.sum(x150 * x186 * x243 * x247)
    result[1, 2, 1] = numpy.sum(x150 * x252 * x286)
    result[1, 2, 2] = numpy.sum(x243 * x256 * x300)
    result[1, 2, 3] = numpy.sum(x150 * x173 * x260 * x266)
    result[1, 2, 4] = numpy.sum(x252 * x260 * x301)
    result[1, 2, 5] = numpy.sum(x158 * x260 * x280)
    result[1, 2, 6] = numpy.sum(x150 * x277 * x298)
    result[1, 2, 7] = numpy.sum(x266 * x278 * x300)
    result[1, 2, 8] = numpy.sum(x158 * x252 * x299)
    result[1, 2, 9] = numpy.sum(x165 * x278 * x302)
    result[1, 3, 0] = numpy.sum(x304 * x41)
    result[1, 3, 1] = numpy.sum(x309 * x310)
    result[1, 3, 2] = numpy.sum(x303 * x310 * x84)
    result[1, 3, 3] = numpy.sum(x173 * x254 * x315)
    result[1, 3, 4] = numpy.sum(x145 * x309 * x316)
    result[1, 3, 5] = numpy.sum(x183 * x254 * x303)
    result[1, 3, 6] = numpy.sum(x323 * x33)
    result[1, 3, 7] = numpy.sum(x116 * x324 * x325)
    result[1, 3, 8] = numpy.sum(x136 * x327 * x33)
    result[1, 3, 9] = numpy.sum(x120 * x329 * x33)
    result[1, 4, 0] = numpy.sum(x192 * x330 * x41)
    result[1, 4, 1] = numpy.sum(x192 * x290 * x75)
    result[1, 4, 2] = numpy.sum(x282 * x301 * x75)
    result[1, 4, 3] = numpy.sum(x192 * x288 * x316)
    result[1, 4, 4] = numpy.sum(181.019335983756 * x160 * x254 * x285)
    result[1, 4, 5] = numpy.sum(x194 * x254 * x331)
    result[1, 4, 6] = numpy.sum(x190 * x297 * x324)
    result[1, 4, 7] = numpy.sum(x288 * x332 * x333)
    result[1, 4, 8] = numpy.sum(x194 * x285 * x332)
    result[1, 4, 9] = numpy.sum(x165 * x330 * x332)
    result[1, 5, 0] = numpy.sum(x233 * x243 * x41)
    result[1, 5, 1] = numpy.sum(x232 * x252 * x75)
    result[1, 5, 2] = numpy.sum(x205 * x334 * x75)
    result[1, 5, 3] = numpy.sum(x232 * x254 * x266)
    result[1, 5, 4] = numpy.sum(x207 * x254 * x335)
    result[1, 5, 5] = numpy.sum(x212 * x254 * x334)
    result[1, 5, 6] = numpy.sum(x277 * x33 * x336)
    result[1, 5, 7] = numpy.sum(x266 * x33 * x338)
    result[1, 5, 8] = numpy.sum(x252 * x33 * x339)
    result[1, 5, 9] = numpy.sum(x215 * x243 * x328 * x33)
    result[1, 6, 0] = numpy.sum(x340 * x343)
    result[1, 6, 1] = numpy.sum(x115 * x344 * x345)
    result[1, 6, 2] = numpy.sum(x343 * x345 * x84)
    result[1, 6, 3] = numpy.sum(x3 * x346)
    result[1, 6, 4] = numpy.sum(46.7389915737742 * x116 * x344 * x347)
    result[1, 6, 5] = numpy.sum(x3 * x342 * x348)
    result[1, 6, 6] = numpy.sum(
        x321
        * x65
        * (
            x0
            * (
                x176 * (x318 + x319)
                + x22 * (x141 + x292 + 4.0 * x293 + 3.0 * x311)
                + 2.0 * x295
                + 2.0 * x296
                + 3.0 * x313
                + 3.0 * x314
            )
            + x122 * x320
        )
    )
    result[1, 6, 7] = numpy.sum(x346 * x84)
    result[1, 6, 8] = numpy.sum(x344 * x348)
    result[1, 6, 9] = numpy.sum(x121 * x326 * x342)
    result[1, 7, 0] = numpy.sum(x150 * x304 * x39)
    result[1, 7, 1] = numpy.sum(x192 * x309 * x349)
    result[1, 7, 2] = numpy.sum(x300 * x303 * x36)
    result[1, 7, 3] = numpy.sum(x189 * x325 * x347)
    result[1, 7, 4] = numpy.sum(x3 * x327 * x333)
    result[1, 7, 5] = numpy.sum(x158 * x303 * x350)
    result[1, 7, 6] = numpy.sum(x150 * x323)
    result[1, 7, 7] = numpy.sum(x153 * x315 * x337)
    result[1, 7, 8] = numpy.sum(x158 * x327)
    result[1, 7, 9] = numpy.sum(x165 * x329)
    result[1, 8, 0] = numpy.sum(x233 * x282 * x39)
    result[1, 8, 1] = numpy.sum(x232 * x285 * x36)
    result[1, 8, 2] = numpy.sum(x205 * x331 * x349)
    result[1, 8, 3] = numpy.sum(x199 * x288 * x3 * x337)
    result[1, 8, 4] = numpy.sum(x207 * x285 * x350)
    result[1, 8, 5] = numpy.sum(x282 * x3 * x339)
    result[1, 8, 6] = numpy.sum(x297 * x336)
    result[1, 8, 7] = numpy.sum(x288 * x338)
    result[1, 8, 8] = numpy.sum(x285 * x339)
    result[1, 8, 9] = numpy.sum(x215 * x283 * x326)
    result[1, 9, 0] = numpy.sum(x235 * x280 * x39)
    result[1, 9, 1] = numpy.sum(x234 * x335 * x345)
    result[1, 9, 2] = numpy.sum(x237 * x302 * x36)
    result[1, 9, 3] = numpy.sum(x234 * x3 * x351)
    result[1, 9, 4] = numpy.sum(x237 * x268 * x350)
    result[1, 9, 5] = numpy.sum(x243 * x3 * x352)
    result[1, 9, 6] = numpy.sum(x235 * x277 * x326)
    result[1, 9, 7] = numpy.sum(x237 * x351)
    result[1, 9, 8] = numpy.sum(x252 * x352)
    result[1, 9, 9] = numpy.sum(x240 * x243 * x326)
    result[2, 0, 0] = numpy.sum(x248 * x355 * x356)
    result[2, 0, 1] = numpy.sum(x355 * x357 * x70)
    result[2, 0, 2] = numpy.sum(x357 * x361)
    result[2, 0, 3] = numpy.sum(x261 * x362 * x97)
    result[2, 0, 4] = numpy.sum(x161 * x261 * x363)
    result[2, 0, 5] = numpy.sum(x261 * x369)
    result[2, 0, 6] = numpy.sum(x114 * x279 * x362)
    result[2, 0, 7] = numpy.sum(x279 * x370 * x97)
    result[2, 0, 8] = numpy.sum(x279 * x369 * x70)
    result[2, 0, 9] = numpy.sum(x279 * x356 * x378)
    result[2, 1, 0] = numpy.sum(x122 * x355 * x379)
    result[2, 1, 1] = numpy.sum(x231 * x256 * x355)
    result[2, 1, 2] = numpy.sum(x122 * x361 * x380)
    result[2, 1, 3] = numpy.sum(x133 * x260 * x362)
    result[2, 1, 4] = numpy.sum(x260 * x361 * x381)
    result[2, 1, 5] = numpy.sum(x122 * x368 * x382)
    result[2, 1, 6] = numpy.sum(x142 * x299 * x355)
    result[2, 1, 7] = numpy.sum(x133 * x299 * x361)
    result[2, 1, 8] = numpy.sum(x231 * x278 * x368)
    result[2, 1, 9] = numpy.sum(x122 * x378 * x383)
    result[2, 2, 0] = numpy.sum(x379 * x385)
    result[2, 2, 1] = numpy.sum(x380 * x385 * x70)
    result[2, 2, 2] = numpy.sum(x380 * x387)
    result[2, 2, 3] = numpy.sum(x155 * x260 * x388)
    result[2, 2, 4] = numpy.sum(x161 * x260 * x389)
    result[2, 2, 5] = numpy.sum(x382 * x391)
    result[2, 2, 6] = numpy.sum(x159 * x299 * x385)
    result[2, 2, 7] = numpy.sum(x155 * x299 * x387)
    result[2, 2, 8] = numpy.sum(x201 * x278 * x391 * x70)
    result[2, 2, 9] = numpy.sum(x383 * x398)
    result[2, 3, 0] = numpy.sum(x168 * x362 * x41)
    result[2, 3, 1] = numpy.sum(x172 * x399 * x75)
    result[2, 3, 2] = numpy.sum(x175 * x370 * x75)
    result[2, 3, 3] = numpy.sum(x180 * x254 * x399)
    result[2, 3, 4] = numpy.sum(x172 * x316 * x370)
    result[2, 3, 5] = numpy.sum(x175 * x368 * x400)
    result[2, 3, 6] = numpy.sum(x33 * x355 * x401)
    result[2, 3, 7] = numpy.sum(x33 * x361 * x402)
    result[2, 3, 8] = numpy.sum(x33 * x368 * x403)
    result[2, 3, 9] = numpy.sum(x33 * x378 * x404)
    result[2, 4, 0] = numpy.sum(x195 * x405 * x41)
    result[2, 4, 1] = numpy.sum(x381 * x385 * x75)
    result[2, 4, 2] = numpy.sum(x195 * x389 * x75)
    result[2, 4, 3] = numpy.sum(x132 * x316 * x388)
    result[2, 4, 4] = numpy.sum(x193 * x254 * x387)
    result[2, 4, 5] = numpy.sum(x195 * x316 * x391)
    result[2, 4, 6] = numpy.sum(x141 * x332 * x405)
    result[2, 4, 7] = numpy.sum(x132 * x332 * x389)
    result[2, 4, 8] = numpy.sum(x332 * x391 * x406)
    result[2, 4, 9] = numpy.sum(46.7389915737742 * x122 * x398 * x409)
    result[2, 5, 0] = numpy.sum(x41 * x411)
    result[2, 5, 1] = numpy.sum(x410 * x412 * x70)
    result[2, 5, 2] = numpy.sum(x412 * x417)
    result[2, 5, 3] = numpy.sum(x155 * x400 * x410)
    result[2, 5, 4] = numpy.sum(x161 * x316 * x417)
    result[2, 5, 5] = numpy.sum(x201 * x254 * x422)
    result[2, 5, 6] = numpy.sum(x159 * x33 * x423)
    result[2, 5, 7] = numpy.sum(x155 * x33 * x424)
    result[2, 5, 8] = numpy.sum(x409 * x425 * x70)
    result[2, 5, 9] = numpy.sum(x33 * x431)
    result[2, 6, 0] = numpy.sum(x224 * x340 * x355)
    result[2, 6, 1] = numpy.sum(x219 * x345 * x362)
    result[2, 6, 2] = numpy.sum(x224 * x345 * x361)
    result[2, 6, 3] = numpy.sum(x3 * x355 * x432)
    result[2, 6, 4] = numpy.sum(x219 * x350 * x363)
    result[2, 6, 5] = numpy.sum(x217 * x3 * x433)
    result[2, 6, 6] = numpy.sum(x225 * x355 * x434)
    result[2, 6, 7] = numpy.sum(x361 * x432)
    result[2, 6, 8] = numpy.sum(x219 * x433)
    result[2, 6, 9] = numpy.sum(x217 * x378 * x434)
    result[2, 7, 0] = numpy.sum(x168 * x388 * x39)
    result[2, 7, 1] = numpy.sum(x172 * x349 * x388)
    result[2, 7, 2] = numpy.sum(x146 * x175 * x36 * x387)
    result[2, 7, 3] = numpy.sum(x3 * x385 * x402)
    result[2, 7, 4] = numpy.sum(x172 * x350 * x389)
    result[2, 7, 5] = numpy.sum(x175 * x350 * x391)
    result[2, 7, 6] = numpy.sum(x385 * x401)
    result[2, 7, 7] = numpy.sum(x387 * x402)
    result[2, 7, 8] = numpy.sum(x391 * x403)
    result[2, 7, 9] = numpy.sum(x398 * x404)
    result[2, 8, 0] = numpy.sum(x122 * x39 * x411)
    result[2, 8, 1] = numpy.sum(x231 * x36 * x410)
    result[2, 8, 2] = numpy.sum(x195 * x349 * x417)
    result[2, 8, 3] = numpy.sum(x133 * x3 * x423)
    result[2, 8, 4] = numpy.sum(x3 * x406 * x424)
    result[2, 8, 5] = numpy.sum(x122 * x425 * x435)
    result[2, 8, 6] = numpy.sum(x142 * x423)
    result[2, 8, 7] = numpy.sum(x133 * x424)
    result[2, 8, 8] = numpy.sum(x128 * x337 * x422)
    result[2, 8, 9] = numpy.sum(x122 * x431)
    result[2, 9, 0] = numpy.sum(x340 * x438)
    result[2, 9, 1] = numpy.sum(x345 * x438 * x70)
    result[2, 9, 2] = numpy.sum(x106 * x345 * x439)
    result[2, 9, 3] = numpy.sum(x3 * x437 * x440)
    result[2, 9, 4] = numpy.sum(46.7389915737742 * x435 * x439 * x70)
    result[2, 9, 5] = numpy.sum(x3 * x441)
    result[2, 9, 6] = numpy.sum(x114 * x326 * x437)
    result[2, 9, 7] = numpy.sum(x439 * x440)
    result[2, 9, 8] = numpy.sum(x441 * x70)
    result[2, 9, 9] = numpy.sum(
        12.0679557322504
        * x408
        * (
            x0
            * (
                x208 * (x427 + x428)
                + x22 * (x165 + x393 + 4.0 * x394 + 3.0 * x418)
                + 2.0 * x396
                + 2.0 * x397
                + 3.0 * x420
                + 3.0 * x421
            )
            + x150 * x429
        )
    )
    return result


def diag_quadrupole3d_34(ax, da, A, bx, db, B, C):
    """Cartesian 3D (fg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 10, 15), dtype=float)

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
    x12 = -x2 - C[0]
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
    x23 = x16 + x22
    x24 = x23 * x3
    x25 = x12**2 * x8
    x26 = x0 * (x15 + x25)
    x27 = x10 + x25
    x28 = x27 * x3
    x29 = 2.0 * x0
    x30 = x17 * x29
    x31 = x28 + x30
    x32 = x3 * x31
    x33 = x26 + x32
    x34 = x3 * x33
    x35 = 3.0 * x18
    x36 = x10 + x9
    x37 = x3 * x36
    x38 = x13 * x29
    x39 = x37 + x38
    x40 = x0 * (3.0 * x20 + x35 + x39)
    x41 = 4.0 * x10
    x42 = x12 * x41
    x43 = 2.0 * x18 + x42
    x44 = x0 * (2.0 * x20 + 2.0 * x28 + x43)
    x45 = 2.0 * x40 + 4.0 * x44
    x46 = x0 * (2.0 * x24 + 4.0 * x34 + x45)
    x47 = -x2 - A[0]
    x48 = x23 * x47
    x49 = 8.0 * x10 * x3
    x50 = x0 * (4.0 * x37 + x49)
    x51 = x0 * (x11 + 3.0 * x9)
    x52 = x3 * x39
    x53 = x51 + x52
    x54 = x47 * x53
    x55 = x50 + x54
    x56 = 4.0 * x16
    x57 = x0 * (4.0 * x22 + x53 + x56)
    x58 = x47 * (x24 + x40)
    x59 = 2.0 * x47
    x60 = 2.0 * x16 + 3.0 * x26
    x61 = x0 * (2.0 * x22 + 3.0 * x32 + x60)
    x62 = x34 + x44
    x63 = x3 * x62
    x64 = x61 + x63
    x65 = x47 * x64
    x66 = x33 * x47
    x67 = x0 * (x34 + x45 + 2.0 * x48 + 3.0 * x66)
    x68 = x47 * x62
    x69 = x61 + x68
    x70 = x47 * x69
    x71 = x46 + x65
    x72 = x0 * (2.0 * x57 + 2.0 * x58 + 5.0 * x61 + x63 + 4.0 * x68) + x47 * x71
    x73 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x74 = da * db * numpy.sqrt(ax**4.5) * numpy.sqrt(bx**5.5)
    x75 = x73 * x74
    x76 = 9.12251705727742 * x75
    x77 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x78 = 0.564189583547756 * x1
    x79 = x77 * x78
    x80 = -x1 * (ax * A[1] + bx * B[1])
    x81 = -x80 - B[1]
    x82 = 24.1359114645008 * x81
    x83 = x21 * x47
    x84 = x39 * x47
    x85 = x51 + x84
    x86 = x31 * x47
    x87 = 2.0 * x86
    x88 = x0 * (x32 + x60 + 2.0 * x83 + x87)
    x89 = x44 + x66
    x90 = x47 * x89
    x91 = x67 + x70
    x92 = x75 * x79
    x93 = x92 * (
        x0
        * (
            x29 * (x22 + x56 + 3.0 * x83 + x85)
            + x59 * (x40 + x48)
            + 2.0 * x61
            + 2.0 * x68
            + 3.0 * x88
            + 3.0 * x90
        )
        + x47 * x91
    )
    x94 = -x1 * (ax * A[2] + bx * B[2])
    x95 = -x94 - B[2]
    x96 = 24.1359114645008 * x95
    x97 = x19 * x47
    x98 = 2.0 * x97
    x99 = x36 * x47
    x100 = x38 + x99
    x101 = x27 * x47
    x102 = x0 * (x101 + x28 + x43 + x98)
    x103 = x26 + x86
    x104 = x103 * x47
    x105 = x88 + x90
    x106 = (
        x0
        * (
            2.0 * x102
            + 2.0 * x104
            + x29 * (x100 + x20 + x35 + x98)
            + 2.0 * x44
            + x59 * (x16 + x83)
            + 2.0 * x66
        )
        + x105 * x47
    )
    x107 = x106 * x75
    x108 = x7 * x77
    x109 = x108 * x81**2
    x110 = x0 * x108
    x111 = x109 + x110
    x112 = 31.1593277158494 * x111
    x113 = 0.318309886183791 * x6
    x114 = x112 * x113
    x115 = 53.9695387335403 * x95
    x116 = x79 * x81
    x117 = x7 * x73
    x118 = x117 * x95**2
    x119 = x0 * x117
    x120 = x118 + x119
    x121 = 31.1593277158494 * x113
    x122 = x74 * x77
    x123 = x121 * x122
    x124 = x111 * x81
    x125 = x108 * x81
    x126 = x125 * x29
    x127 = x124 + x126
    x128 = 24.1359114645008 * x127
    x129 = x13 * x47
    x130 = x17 * x47
    x131 = x101 + x30
    x132 = x0 * (x11 + x17 * x59 + x25) + x131 * x47
    x133 = x102 + x104
    x134 = x113 * (
        x0
        * (x132 + 2.0 * x26 + x29 * (x11 + x129 + x130 + x14) + x59 * (x18 + x97) + x87)
        + x133 * x47
    )
    x135 = x134 * x75
    x136 = 53.9695387335403 * x111
    x137 = x122 * x134
    x138 = 53.9695387335403 * x120
    x139 = x120 * x95
    x140 = x117 * x95
    x141 = x140 * x29
    x142 = x139 + x141
    x143 = 24.1359114645008 * x142
    x144 = 3.0 * x110
    x145 = x0 * (3.0 * x109 + x144)
    x146 = x127 * x81
    x147 = x145 + x146
    x148 = 9.12251705727742 * x147
    x149 = x47 * x8
    x150 = x0 * (2.0 * x101 + x29 * (x149 + x17) + x42 + x59 * (x10 + x130)) + x132 * x47
    x151 = x113 * x150
    x152 = x151 * x75
    x153 = 0.179587122125167 * x74
    x154 = x120 * x153
    x155 = x122 * x151
    x156 = 3.0 * x119
    x157 = x0 * (3.0 * x118 + x156)
    x158 = x142 * x95
    x159 = x157 + x158
    x160 = 9.12251705727742 * x159
    x161 = -x80 - A[1]
    x162 = 20.3985682659737 * x161
    x163 = x72 * x92
    x164 = 53.9695387335403 * x91
    x165 = x108 * x161
    x166 = x165 * x81
    x167 = x110 + x166
    x168 = x113 * x75
    x169 = x167 * x168
    x170 = x164 * x92
    x171 = x111 * x161
    x172 = x126 + x171
    x173 = 69.6743749058326 * x172
    x174 = x105 * x168
    x175 = 120.679557322504 * x105
    x176 = x113 * x122
    x177 = 69.6743749058326 * x176
    x178 = x105 * x177
    x179 = x127 * x161
    x180 = x145 + x179
    x181 = 53.9695387335403 * x180
    x182 = x168 * x181
    x183 = 120.679557322504 * x133
    x184 = x168 * x95
    x185 = 53.9695387335403 * x142
    x186 = x133 * x176
    x187 = x110 * x81
    x188 = 8.0 * x187
    x189 = x0 * (4.0 * x124 + x188)
    x190 = x147 * x161
    x191 = x189 + x190
    x192 = 20.3985682659737 * x191
    x193 = x132 * x168
    x194 = 69.6743749058326 * x132
    x195 = x153 * x167
    x196 = x132 * x176
    x197 = -x94 - A[2]
    x198 = 20.3985682659737 * x197
    x199 = x117 * x197
    x200 = x199 * x95
    x201 = x119 + x200
    x202 = x176 * x201
    x203 = 69.6743749058326 * x111
    x204 = x120 * x197
    x205 = x141 + x204
    x206 = 53.9695387335403 * x127
    x207 = x168 * x197
    x208 = x153 * x201
    x209 = x176 * x205
    x210 = x142 * x197
    x211 = x157 + x210
    x212 = 53.9695387335403 * x211
    x213 = 20.3985682659737 * x147
    x214 = x153 * x205
    x215 = x119 * x95
    x216 = 8.0 * x215
    x217 = x0 * (4.0 * x139 + x216)
    x218 = x159 * x197
    x219 = x217 + x218
    x220 = 20.3985682659737 * x219
    x221 = x108 * x161**2
    x222 = x110 + x221
    x223 = 20.3985682659737 * x222
    x224 = x168 * x223
    x225 = x0 * (x125 + x165)
    x226 = x161 * x167
    x227 = x225 + x226
    x228 = 53.9695387335403 * x168
    x229 = x228 * x69
    x230 = 2.0 * x161
    x231 = x125 * x230 + x144
    x232 = x0 * (x109 + x231)
    x233 = x161 * x172
    x234 = x232 + x233
    x235 = 69.6743749058326 * x168
    x236 = x234 * x235
    x237 = 120.679557322504 * x227
    x238 = 69.6743749058326 * x154
    x239 = x0 * (x124 + 3.0 * x171 + x188)
    x240 = x161 * x180
    x241 = x239 + x240
    x242 = x228 * x241
    x243 = 120.679557322504 * x103
    x244 = x153 * x185
    x245 = 20.3985682659737 * x131
    x246 = x0 * (5.0 * x145 + x146 + 4.0 * x179) + x161 * x191
    x247 = x168 * x246
    x248 = x131 * x153
    x249 = 93.4779831475484 * x69
    x250 = 120.679557322504 * x172
    x251 = 209.023124717498 * x195
    x252 = 120.679557322504 * x161
    x253 = 93.4779831475484 * x103
    x254 = 209.023124717498 * x208
    x255 = x161 * x176
    x256 = 35.3313566383285 * x131
    x257 = 93.4779831475484 * x131
    x258 = 120.679557322504 * x214
    x259 = x117 * x197**2
    x260 = x119 + x259
    x261 = 20.3985682659737 * x260
    x262 = x176 * x261
    x263 = 53.9695387335403 * x176
    x264 = x263 * x69
    x265 = x0 * (x140 + x199)
    x266 = x197 * x201
    x267 = x265 + x266
    x268 = x153 * x260
    x269 = x176 * x81
    x270 = 120.679557322504 * x267
    x271 = 2.0 * x197
    x272 = x140 * x271 + x156
    x273 = x0 * (x118 + x272)
    x274 = x197 * x205
    x275 = x273 + x274
    x276 = x177 * x275
    x277 = x153 * x267
    x278 = x0 * (x139 + 3.0 * x204 + x216)
    x279 = x197 * x211
    x280 = x278 + x279
    x281 = x263 * x280
    x282 = x0 * (5.0 * x157 + x158 + 4.0 * x210) + x197 * x219
    x283 = x176 * x282
    x284 = x161 * x222 + x165 * x29
    x285 = x113 * x76
    x286 = x0 * (x221 + x231) + x161 * x227
    x287 = 24.1359114645008 * x168
    x288 = x287 * x62
    x289 = 2.0 * x0 * (x171 + 2.0 * x187 + x225 + x226) + x161 * x234
    x290 = x121 * x75
    x291 = 53.9695387335403 * x286
    x292 = 31.1593277158494 * x154
    x293 = x0 * (2.0 * x145 + 2.0 * x179 + 3.0 * x232 + 3.0 * x233) + x161 * x241
    x294 = x287 * x293
    x295 = 53.9695387335403 * x31
    x296 = x153 * x284
    x297 = 2.0 * x0 * (x189 + x190 + 2.0 * x239 + 2.0 * x240) + x161 * x246
    x298 = x153 * x27
    x299 = x197 * x228
    x300 = 53.9695387335403 * x208
    x301 = 120.679557322504 * x208
    x302 = 69.6743749058326 * x214
    x303 = 120.679557322504 * x31
    x304 = x153 * x295
    x305 = 53.9695387335403 * x195
    x306 = 120.679557322504 * x195
    x307 = 20.3985682659737 * x268
    x308 = 53.9695387335403 * x277
    x309 = x197 * x260 + x199 * x29
    x310 = 9.12251705727742 * x176
    x311 = 24.1359114645008 * x176
    x312 = x311 * x62
    x313 = x0 * (x259 + x272) + x197 * x267
    x314 = x153 * x309
    x315 = x263 * x81
    x316 = 2.0 * x0 * (x204 + 2.0 * x215 + x265 + x266) + x197 * x275
    x317 = x0 * (2.0 * x157 + 2.0 * x210 + 3.0 * x273 + 3.0 * x274) + x197 * x280
    x318 = x311 * x317
    x319 = 2.0 * x0 * (x217 + x218 + 2.0 * x278 + 2.0 * x279) + x197 * x282
    x320 = -x80 - C[1]
    x321 = x108 * x320**2
    x322 = x110 + x321
    x323 = x0 * (x37 + x49 + 3.0 * x99)
    x324 = x47 * x85
    x325 = x0 * (5.0 * x51 + x52 + 4.0 * x84) + x47 * x55
    x326 = 2.0 * x0 * (2.0 * x323 + 2.0 * x324 + x50 + x54) + x325 * x47
    x327 = x322 * x81
    x328 = x108 * x320
    x329 = x29 * x328
    x330 = x327 + x329
    x331 = x11 + x13 * x59
    x332 = x0 * (x331 + x9)
    x333 = x100 * x47
    x334 = x323 + x324
    x335 = x0 * (3.0 * x332 + 3.0 * x333 + 2.0 * x51 + 2.0 * x84) + x334 * x47
    x336 = x287 * x335
    x337 = x125 * x320
    x338 = x144 + 2.0 * x337
    x339 = x0 * (x321 + x338)
    x340 = x330 * x81
    x341 = x339 + x340
    x342 = x0 * (x13 + x149)
    x343 = x10 + x129
    x344 = x343 * x47
    x345 = x332 + x333
    x346 = x0 * (x3 * x41 + 2.0 * x342 + 2.0 * x344 + 2.0 * x99) + x345 * x47
    x347 = 53.9695387335403 * x330
    x348 = x47**2 * x8
    x349 = x342 + x344
    x350 = x0 * (x331 + x348) + x349 * x47
    x351 = x110 + x337
    x352 = x351 * x81
    x353 = x0 * (x125 + x328)
    x354 = 4.0 * x110 * x320
    x355 = 2.0 * x353 + x354
    x356 = x0 * (2.0 * x327 + 2.0 * x352 + x355)
    x357 = x341 * x81
    x358 = x356 + x357
    x359 = x287 * x358
    x360 = x228 * x95
    x361 = x143 * x153
    x362 = x352 + x353
    x363 = x362 * x81
    x364 = x0 * (x109 + x338)
    x365 = 3.0 * x339 + 2.0 * x364
    x366 = x0 * (3.0 * x340 + 2.0 * x363 + x365)
    x367 = x358 * x81
    x368 = x366 + x367
    x369 = x10 + x348
    x370 = x149 * x29 + x369 * x47
    x371 = x153 * x322
    x372 = x161 * x322
    x373 = x329 + x372
    x374 = 20.3985682659737 * x373
    x375 = x168 * x325
    x376 = x161 * x330
    x377 = x339 + x376
    x378 = x228 * x334
    x379 = x161 * x341
    x380 = x356 + x379
    x381 = x235 * x345
    x382 = 120.679557322504 * x184
    x383 = x161 * x358
    x384 = x366 + x383
    x385 = x228 * x384
    x386 = 120.679557322504 * x349
    x387 = x363 + x364
    x388 = x387 * x81
    x389 = 3.0 * x353
    x390 = x0 * (x127 + 3.0 * x352 + x389)
    x391 = 4.0 * x356 + 2.0 * x390
    x392 = x0 * (4.0 * x357 + 2.0 * x388 + x391)
    x393 = x161 * x368
    x394 = x392 + x393
    x395 = 20.3985682659737 * x168
    x396 = x369 * x395
    x397 = x153 * x369
    x398 = x0 * (x144 + x230 * x328 + x321) + x161 * x373
    x399 = x395 * x398
    x400 = x161 * x351
    x401 = 2.0 * x400
    x402 = x0 * (x327 + x355 + x372 + x401)
    x403 = x161 * x377
    x404 = x402 + x403
    x405 = 53.9695387335403 * x85
    x406 = x168 * x405
    x407 = 69.6743749058326 * x100
    x408 = x161 * x362
    x409 = 2.0 * x376
    x410 = x0 * (x340 + x365 + 2.0 * x408 + x409)
    x411 = x161 * x380
    x412 = x410 + x411
    x413 = x168 * x412
    x414 = x161 * x387
    x415 = x0 * (x357 + 3.0 * x379 + x391 + 2.0 * x414)
    x416 = x161 * x384
    x417 = x415 + x416
    x418 = 53.9695387335403 * x343
    x419 = 120.679557322504 * x343
    x420 = 20.3985682659737 * x47
    x421 = 4.0 * x364
    x422 = x0 * (x147 + 4.0 * x363 + x421)
    x423 = x161 * (x388 + x390)
    x424 = x0 * (5.0 * x366 + x367 + 4.0 * x383 + 2.0 * x422 + 2.0 * x423) + x161 * x394
    x425 = x5 * x78
    x426 = x425 * x75
    x427 = x424 * x426
    x428 = x115 * x426
    x429 = x5 * x74
    x430 = x113 * x429
    x431 = 69.6743749058326 * x430
    x432 = x412 * x431
    x433 = x404 * x430
    x434 = x398 * x430
    x435 = 35.3313566383285 * x373
    x436 = 93.4779831475484 * x207
    x437 = 93.4779831475484 * x373
    x438 = 120.679557322504 * x380
    x439 = 209.023124717498 * x343
    x440 = x153 * x343
    x441 = 35.3313566383285 * x47
    x442 = x197 * x426
    x443 = x430 * x47
    x444 = 93.4779831475484 * x443
    x445 = 120.679557322504 * x277
    x446 = x153 * x275
    x447 = x261 * x430
    x448 = 53.9695387335403 * x430
    x449 = x267 * x448
    x450 = x275 * x431
    x451 = x282 * x430
    x452 = 9.12251705727742 * x53
    x453 = x165 * x320
    x454 = (
        x0 * (x230 * (x110 + x453) + x29 * (x165 + x328) + x354 + 2.0 * x372)
        + x161 * x398
    )
    x455 = x168 * x454
    x456 = (
        x0
        * (
            x230 * (x353 + x400)
            + x29 * (x144 + x166 + x337 + x453)
            + 2.0 * x339
            + x398
            + x409
        )
        + x161 * x404
    )
    x457 = 24.1359114645008 * x39
    x458 = (
        x0
        * (
            x230 * (x364 + x408)
            + x29 * (x172 + x352 + x389 + x401)
            + 2.0 * x356
            + 2.0 * x379
            + 2.0 * x402
            + 2.0 * x403
        )
        + x161 * x412
    )
    x459 = 31.1593277158494 * x36
    x460 = x426 * (
        x0
        * (
            x230 * (x390 + x414)
            + x29 * (x180 + x363 + 3.0 * x408 + x421)
            + 2.0 * x366
            + 2.0 * x383
            + 3.0 * x410
            + 3.0 * x411
        )
        + x161 * x417
    )
    x461 = 24.1359114645008 * x3
    x462 = x430 * x456
    x463 = x430 * x454
    x464 = x121 * x429
    x465 = 53.9695387335403 * x39
    x466 = 69.6743749058326 * x36
    x467 = 53.9695387335403 * x3
    x468 = 120.679557322504 * x3
    x469 = x430 * x468
    x470 = x3 * x448
    x471 = x3 * x430
    x472 = x280 * x448
    x473 = x153 * x36
    x474 = 24.1359114645008 * x430
    x475 = x358 * x474
    x476 = x317 * x474
    x477 = 9.12251705727742 * x430
    x478 = -x94 - C[2]
    x479 = x117 * x478**2
    x480 = x119 + x479
    x481 = x311 * x335
    x482 = x480 * x95
    x483 = x117 * x478
    x484 = x29 * x483
    x485 = x482 + x484
    x486 = x153 * x480
    x487 = x140 * x478
    x488 = x156 + 2.0 * x487
    x489 = x0 * (x479 + x488)
    x490 = x485 * x95
    x491 = x489 + x490
    x492 = x153 * x485
    x493 = x119 + x487
    x494 = x493 * x95
    x495 = x0 * (x140 + x483)
    x496 = 4.0 * x119 * x478
    x497 = 2.0 * x495 + x496
    x498 = x0 * (2.0 * x482 + 2.0 * x494 + x497)
    x499 = x491 * x95
    x500 = x498 + x499
    x501 = x311 * x500
    x502 = x153 * x491
    x503 = x494 + x495
    x504 = x503 * x95
    x505 = x0 * (x118 + x488)
    x506 = 3.0 * x489 + 2.0 * x505
    x507 = x0 * (3.0 * x490 + 2.0 * x504 + x506)
    x508 = x500 * x95
    x509 = x507 + x508
    x510 = 20.3985682659737 * x176
    x511 = x325 * x510
    x512 = x263 * x334
    x513 = x177 * x345
    x514 = x263 * x349
    x515 = x369 * x510
    x516 = x197 * x480
    x517 = x484 + x516
    x518 = x197 * x485
    x519 = x489 + x518
    x520 = x153 * x517
    x521 = 120.679557322504 * x269
    x522 = x197 * x491
    x523 = x498 + x522
    x524 = x153 * x519
    x525 = x197 * x500
    x526 = x507 + x525
    x527 = x504 + x505
    x528 = x527 * x95
    x529 = 3.0 * x495
    x530 = x0 * (x142 + 3.0 * x494 + x529)
    x531 = 4.0 * x498 + 2.0 * x530
    x532 = x0 * (4.0 * x499 + 2.0 * x528 + x531)
    x533 = x197 * x509
    x534 = x532 + x533
    x535 = 53.9695387335403 * x440
    x536 = 20.3985682659737 * x246
    x537 = x241 * x448
    x538 = x234 * x431
    x539 = x227 * x448
    x540 = x223 * x430
    x541 = 93.4779831475484 * x517
    x542 = 93.4779831475484 * x255
    x543 = x430 * x517
    x544 = x429 * x79
    x545 = x161 * x544
    x546 = x0 * (x156 + x271 * x483 + x479) + x197 * x517
    x547 = x510 * x546
    x548 = x176 * x405
    x549 = x197 * x493
    x550 = 2.0 * x549
    x551 = x0 * (x482 + x497 + x516 + x550)
    x552 = x197 * x519
    x553 = x551 + x552
    x554 = x197 * x503
    x555 = 2.0 * x518
    x556 = x0 * (x490 + x506 + 2.0 * x554 + x555)
    x557 = x197 * x523
    x558 = x556 + x557
    x559 = x176 * x558
    x560 = x197 * x527
    x561 = x0 * (x499 + 3.0 * x522 + x531 + 2.0 * x560)
    x562 = x197 * x526
    x563 = x561 + x562
    x564 = x430 * x546
    x565 = x430 * x553
    x566 = x430 * x558
    x567 = 4.0 * x505
    x568 = x0 * (x159 + 4.0 * x504 + x567)
    x569 = x197 * (x528 + x530)
    x570 = x0 * (5.0 * x507 + x508 + 4.0 * x525 + 2.0 * x568 + 2.0 * x569) + x197 * x534
    x571 = x544 * x570
    x572 = x293 * x474
    x573 = x474 * x500
    x574 = x199 * x478
    x575 = (
        x0 * (x271 * (x119 + x574) + x29 * (x199 + x483) + x496 + 2.0 * x516)
        + x197 * x546
    )
    x576 = x176 * x575
    x577 = (
        x0
        * (
            x271 * (x495 + x549)
            + x29 * (x156 + x200 + x487 + x574)
            + 2.0 * x489
            + x546
            + x555
        )
        + x197 * x553
    )
    x578 = (
        x0
        * (
            x271 * (x505 + x554)
            + x29 * (x205 + x494 + x529 + x550)
            + 2.0 * x498
            + 2.0 * x522
            + 2.0 * x551
            + 2.0 * x552
        )
        + x197 * x558
    )
    x579 = x430 * x575
    x580 = x430 * x577
    x581 = x429 * x578
    x582 = x544 * (
        x0
        * (
            x271 * (x530 + x560)
            + x29 * (x211 + x504 + 3.0 * x554 + x567)
            + 2.0 * x507
            + 2.0 * x525
            + 3.0 * x556
            + 3.0 * x557
        )
        + x197 * x563
    )

    # 450 item(s)
    result[0, 0, 0] = numpy.sum(
        x76
        * x79
        * (
            x0
            * (
                x29 * (x24 + 5.0 * x40 + 4.0 * x48 + x55)
                + 2.0 * x46
                + x59 * (x57 + x58)
                + 2.0 * x65
                + 4.0 * x67
                + 4.0 * x70
            )
            + x47 * x72
        )
    )
    result[0, 0, 1] = numpy.sum(x82 * x93)
    result[0, 0, 2] = numpy.sum(x93 * x96)
    result[0, 0, 3] = numpy.sum(x107 * x114)
    result[0, 0, 4] = numpy.sum(x107 * x115 * x116)
    result[0, 0, 5] = numpy.sum(x106 * x120 * x123)
    result[0, 0, 6] = numpy.sum(x128 * x135)
    result[0, 0, 7] = numpy.sum(x135 * x136 * x95)
    result[0, 0, 8] = numpy.sum(x137 * x138 * x81)
    result[0, 0, 9] = numpy.sum(x137 * x143)
    result[0, 0, 10] = numpy.sum(x148 * x152)
    result[0, 0, 11] = numpy.sum(x128 * x152 * x95)
    result[0, 0, 12] = numpy.sum(x112 * x150 * x154)
    result[0, 0, 13] = numpy.sum(x143 * x155 * x81)
    result[0, 0, 14] = numpy.sum(x155 * x160)
    result[0, 1, 0] = numpy.sum(x162 * x163)
    result[0, 1, 1] = numpy.sum(x164 * x169)
    result[0, 1, 2] = numpy.sum(x161 * x170 * x95)
    result[0, 1, 3] = numpy.sum(x173 * x174)
    result[0, 1, 4] = numpy.sum(x169 * x175 * x95)
    result[0, 1, 5] = numpy.sum(x120 * x161 * x178)
    result[0, 1, 6] = numpy.sum(x133 * x182)
    result[0, 1, 7] = numpy.sum(x172 * x183 * x184)
    result[0, 1, 8] = numpy.sum(x154 * x167 * x183)
    result[0, 1, 9] = numpy.sum(x161 * x185 * x186)
    result[0, 1, 10] = numpy.sum(x192 * x193)
    result[0, 1, 11] = numpy.sum(x132 * x182 * x95)
    result[0, 1, 12] = numpy.sum(x154 * x172 * x194)
    result[0, 1, 13] = numpy.sum(x132 * x185 * x195)
    result[0, 1, 14] = numpy.sum(x159 * x162 * x196)
    result[0, 2, 0] = numpy.sum(x163 * x198)
    result[0, 2, 1] = numpy.sum(x170 * x197 * x81)
    result[0, 2, 2] = numpy.sum(x164 * x202)
    result[0, 2, 3] = numpy.sum(x174 * x197 * x203)
    result[0, 2, 4] = numpy.sum(x175 * x202 * x81)
    result[0, 2, 5] = numpy.sum(x178 * x205)
    result[0, 2, 6] = numpy.sum(x133 * x206 * x207)
    result[0, 2, 7] = numpy.sum(x111 * x183 * x208)
    result[0, 2, 8] = numpy.sum(x183 * x209 * x81)
    result[0, 2, 9] = numpy.sum(x186 * x212)
    result[0, 2, 10] = numpy.sum(x193 * x197 * x213)
    result[0, 2, 11] = numpy.sum(x132 * x206 * x208)
    result[0, 2, 12] = numpy.sum(x111 * x194 * x214)
    result[0, 2, 13] = numpy.sum(x196 * x212 * x81)
    result[0, 2, 14] = numpy.sum(x196 * x220)
    result[0, 3, 0] = numpy.sum(x224 * x71)
    result[0, 3, 1] = numpy.sum(x227 * x229)
    result[0, 3, 2] = numpy.sum(x222 * x229 * x95)
    result[0, 3, 3] = numpy.sum(x236 * x89)
    result[0, 3, 4] = numpy.sum(x184 * x237 * x89)
    result[0, 3, 5] = numpy.sum(x222 * x238 * x89)
    result[0, 3, 6] = numpy.sum(x103 * x242)
    result[0, 3, 7] = numpy.sum(x184 * x234 * x243)
    result[0, 3, 8] = numpy.sum(x154 * x227 * x243)
    result[0, 3, 9] = numpy.sum(x103 * x222 * x244)
    result[0, 3, 10] = numpy.sum(x245 * x247)
    result[0, 3, 11] = numpy.sum(x131 * x242 * x95)
    result[0, 3, 12] = numpy.sum(x131 * x234 * x238)
    result[0, 3, 13] = numpy.sum(x131 * x227 * x244)
    result[0, 3, 14] = numpy.sum(x159 * x223 * x248)
    result[0, 4, 0] = numpy.sum(35.3313566383285 * x161 * x197 * x71 * x92)
    result[0, 4, 1] = numpy.sum(x169 * x197 * x249)
    result[0, 4, 2] = numpy.sum(x161 * x202 * x249)
    result[0, 4, 3] = numpy.sum(x207 * x250 * x89)
    result[0, 4, 4] = numpy.sum(x201 * x251 * x89)
    result[0, 4, 5] = numpy.sum(x209 * x252 * x89)
    result[0, 4, 6] = numpy.sum(x180 * x207 * x253)
    result[0, 4, 7] = numpy.sum(x103 * x172 * x254)
    result[0, 4, 8] = numpy.sum(x103 * x205 * x251)
    result[0, 4, 9] = numpy.sum(x211 * x253 * x255)
    result[0, 4, 10] = numpy.sum(x191 * x207 * x256)
    result[0, 4, 11] = numpy.sum(x180 * x208 * x257)
    result[0, 4, 12] = numpy.sum(x131 * x172 * x258)
    result[0, 4, 13] = numpy.sum(x195 * x211 * x257)
    result[0, 4, 14] = numpy.sum(x219 * x255 * x256)
    result[0, 5, 0] = numpy.sum(x262 * x71)
    result[0, 5, 1] = numpy.sum(x260 * x264 * x81)
    result[0, 5, 2] = numpy.sum(x264 * x267)
    result[0, 5, 3] = numpy.sum(x203 * x268 * x89)
    result[0, 5, 4] = numpy.sum(x269 * x270 * x89)
    result[0, 5, 5] = numpy.sum(x276 * x89)
    result[0, 5, 6] = numpy.sum(x103 * x206 * x268)
    result[0, 5, 7] = numpy.sum(x111 * x243 * x277)
    result[0, 5, 8] = numpy.sum(x243 * x269 * x275)
    result[0, 5, 9] = numpy.sum(x103 * x281)
    result[0, 5, 10] = numpy.sum(x213 * x248 * x260)
    result[0, 5, 11] = numpy.sum(x206 * x248 * x267)
    result[0, 5, 12] = numpy.sum(x203 * x248 * x275)
    result[0, 5, 13] = numpy.sum(x131 * x281 * x81)
    result[0, 5, 14] = numpy.sum(x245 * x283)
    result[0, 6, 0] = numpy.sum(x284 * x285 * x64)
    result[0, 6, 1] = numpy.sum(x286 * x288)
    result[0, 6, 2] = numpy.sum(x284 * x288 * x95)
    result[0, 6, 3] = numpy.sum(x289 * x290 * x33)
    result[0, 6, 4] = numpy.sum(x184 * x291 * x33)
    result[0, 6, 5] = numpy.sum(x284 * x292 * x33)
    result[0, 6, 6] = numpy.sum(x294 * x31)
    result[0, 6, 7] = numpy.sum(x184 * x289 * x295)
    result[0, 6, 8] = numpy.sum(x154 * x286 * x295)
    result[0, 6, 9] = numpy.sum(x143 * x296 * x31)
    result[0, 6, 10] = numpy.sum(x27 * x285 * x297)
    result[0, 6, 11] = numpy.sum(x27 * x294 * x95)
    result[0, 6, 12] = numpy.sum(x27 * x289 * x292)
    result[0, 6, 13] = numpy.sum(x143 * x286 * x298)
    result[0, 6, 14] = numpy.sum(x160 * x27 * x296)
    result[0, 7, 0] = numpy.sum(x197 * x224 * x64)
    result[0, 7, 1] = numpy.sum(x227 * x299 * x62)
    result[0, 7, 2] = numpy.sum(x222 * x300 * x62)
    result[0, 7, 3] = numpy.sum(x197 * x236 * x33)
    result[0, 7, 4] = numpy.sum(x227 * x301 * x33)
    result[0, 7, 5] = numpy.sum(x222 * x302 * x33)
    result[0, 7, 6] = numpy.sum(x207 * x241 * x295)
    result[0, 7, 7] = numpy.sum(x234 * x301 * x31)
    result[0, 7, 8] = numpy.sum(x214 * x227 * x303)
    result[0, 7, 9] = numpy.sum(x211 * x222 * x304)
    result[0, 7, 10] = numpy.sum(x198 * x247 * x27)
    result[0, 7, 11] = numpy.sum(x241 * x27 * x300)
    result[0, 7, 12] = numpy.sum(x234 * x27 * x302)
    result[0, 7, 13] = numpy.sum(x212 * x227 * x298)
    result[0, 7, 14] = numpy.sum(x219 * x223 * x298)
    result[0, 8, 0] = numpy.sum(x161 * x262 * x64)
    result[0, 8, 1] = numpy.sum(x260 * x305 * x62)
    result[0, 8, 2] = numpy.sum(x161 * x263 * x267 * x62)
    result[0, 8, 3] = numpy.sum(x173 * x268 * x33)
    result[0, 8, 4] = numpy.sum(x267 * x306 * x33)
    result[0, 8, 5] = numpy.sum(x161 * x276 * x33)
    result[0, 8, 6] = numpy.sum(x180 * x268 * x295)
    result[0, 8, 7] = numpy.sum(x172 * x277 * x303)
    result[0, 8, 8] = numpy.sum(x195 * x275 * x303)
    result[0, 8, 9] = numpy.sum(x255 * x280 * x295)
    result[0, 8, 10] = numpy.sum(x191 * x27 * x307)
    result[0, 8, 11] = numpy.sum(x180 * x27 * x308)
    result[0, 8, 12] = numpy.sum(x173 * x275 * x298)
    result[0, 8, 13] = numpy.sum(x27 * x280 * x305)
    result[0, 8, 14] = numpy.sum(x162 * x27 * x283)
    result[0, 9, 0] = numpy.sum(x309 * x310 * x64)
    result[0, 9, 1] = numpy.sum(x309 * x312 * x81)
    result[0, 9, 2] = numpy.sum(x312 * x313)
    result[0, 9, 3] = numpy.sum(x112 * x314 * x33)
    result[0, 9, 4] = numpy.sum(x313 * x315 * x33)
    result[0, 9, 5] = numpy.sum(x123 * x316 * x33)
    result[0, 9, 6] = numpy.sum(x128 * x31 * x314)
    result[0, 9, 7] = numpy.sum(x111 * x304 * x313)
    result[0, 9, 8] = numpy.sum(x269 * x295 * x316)
    result[0, 9, 9] = numpy.sum(x31 * x318)
    result[0, 9, 10] = numpy.sum(x148 * x298 * x309)
    result[0, 9, 11] = numpy.sum(x128 * x298 * x313)
    result[0, 9, 12] = numpy.sum(x112 * x298 * x316)
    result[0, 9, 13] = numpy.sum(x27 * x318 * x81)
    result[0, 9, 14] = numpy.sum(x27 * x310 * x319)
    result[1, 0, 0] = numpy.sum(x285 * x322 * x326)
    result[1, 0, 1] = numpy.sum(x330 * x336)
    result[1, 0, 2] = numpy.sum(x322 * x336 * x95)
    result[1, 0, 3] = numpy.sum(x290 * x341 * x346)
    result[1, 0, 4] = numpy.sum(x184 * x346 * x347)
    result[1, 0, 5] = numpy.sum(x292 * x322 * x346)
    result[1, 0, 6] = numpy.sum(x350 * x359)
    result[1, 0, 7] = numpy.sum(x341 * x350 * x360)
    result[1, 0, 8] = numpy.sum(x154 * x347 * x350)
    result[1, 0, 9] = numpy.sum(x322 * x350 * x361)
    result[1, 0, 10] = numpy.sum(x285 * x368 * x370)
    result[1, 0, 11] = numpy.sum(x359 * x370 * x95)
    result[1, 0, 12] = numpy.sum(x292 * x341 * x370)
    result[1, 0, 13] = numpy.sum(x330 * x361 * x370)
    result[1, 0, 14] = numpy.sum(x160 * x370 * x371)
    result[1, 1, 0] = numpy.sum(x374 * x375)
    result[1, 1, 1] = numpy.sum(x377 * x378)
    result[1, 1, 2] = numpy.sum(x373 * x378 * x95)
    result[1, 1, 3] = numpy.sum(x380 * x381)
    result[1, 1, 4] = numpy.sum(x345 * x377 * x382)
    result[1, 1, 5] = numpy.sum(x238 * x345 * x373)
    result[1, 1, 6] = numpy.sum(x349 * x385)
    result[1, 1, 7] = numpy.sum(x184 * x380 * x386)
    result[1, 1, 8] = numpy.sum(x154 * x377 * x386)
    result[1, 1, 9] = numpy.sum(x244 * x349 * x373)
    result[1, 1, 10] = numpy.sum(x394 * x396)
    result[1, 1, 11] = numpy.sum(x369 * x385 * x95)
    result[1, 1, 12] = numpy.sum(x238 * x369 * x380)
    result[1, 1, 13] = numpy.sum(x244 * x369 * x377)
    result[1, 1, 14] = numpy.sum(x159 * x374 * x397)
    result[1, 2, 0] = numpy.sum(x198 * x322 * x375)
    result[1, 2, 1] = numpy.sum(x207 * x334 * x347)
    result[1, 2, 2] = numpy.sum(x300 * x322 * x334)
    result[1, 2, 3] = numpy.sum(x197 * x341 * x381)
    result[1, 2, 4] = numpy.sum(x301 * x330 * x345)
    result[1, 2, 5] = numpy.sum(x302 * x322 * x345)
    result[1, 2, 6] = numpy.sum(x299 * x349 * x358)
    result[1, 2, 7] = numpy.sum(x301 * x341 * x349)
    result[1, 2, 8] = numpy.sum(x214 * x330 * x386)
    result[1, 2, 9] = numpy.sum(x212 * x349 * x371)
    result[1, 2, 10] = numpy.sum(x197 * x368 * x396)
    result[1, 2, 11] = numpy.sum(x300 * x358 * x369)
    result[1, 2, 12] = numpy.sum(x302 * x341 * x369)
    result[1, 2, 13] = numpy.sum(x212 * x330 * x397)
    result[1, 2, 14] = numpy.sum(x220 * x369 * x371)
    result[1, 3, 0] = numpy.sum(x399 * x55)
    result[1, 3, 1] = numpy.sum(x404 * x406)
    result[1, 3, 2] = numpy.sum(x398 * x406 * x95)
    result[1, 3, 3] = numpy.sum(x407 * x413)
    result[1, 3, 4] = numpy.sum(x100 * x382 * x404)
    result[1, 3, 5] = numpy.sum(x100 * x238 * x398)
    result[1, 3, 6] = numpy.sum(x168 * x417 * x418)
    result[1, 3, 7] = numpy.sum(x413 * x419 * x95)
    result[1, 3, 8] = numpy.sum(x154 * x404 * x419)
    result[1, 3, 9] = numpy.sum(x244 * x343 * x398)
    result[1, 3, 10] = numpy.sum(x420 * x427)
    result[1, 3, 11] = numpy.sum(x417 * x428 * x47)
    result[1, 3, 12] = numpy.sum(x120 * x432 * x47)
    result[1, 3, 13] = numpy.sum(x185 * x433 * x47)
    result[1, 3, 14] = numpy.sum(x159 * x420 * x434)
    result[1, 4, 0] = numpy.sum(x207 * x435 * x55)
    result[1, 4, 1] = numpy.sum(x377 * x436 * x85)
    result[1, 4, 2] = numpy.sum(x208 * x437 * x85)
    result[1, 4, 3] = numpy.sum(x100 * x207 * x438)
    result[1, 4, 4] = numpy.sum(x100 * x254 * x377)
    result[1, 4, 5] = numpy.sum(x100 * x258 * x373)
    result[1, 4, 6] = numpy.sum(x343 * x384 * x436)
    result[1, 4, 7] = numpy.sum(x254 * x343 * x380)
    result[1, 4, 8] = numpy.sum(x214 * x377 * x439)
    result[1, 4, 9] = numpy.sum(x211 * x437 * x440)
    result[1, 4, 10] = numpy.sum(x394 * x441 * x442)
    result[1, 4, 11] = numpy.sum(x201 * x384 * x444)
    result[1, 4, 12] = numpy.sum(x205 * x438 * x443)
    result[1, 4, 13] = numpy.sum(x211 * x377 * x444)
    result[1, 4, 14] = numpy.sum(x219 * x435 * x443)
    result[1, 5, 0] = numpy.sum(x307 * x322 * x55)
    result[1, 5, 1] = numpy.sum(x268 * x347 * x85)
    result[1, 5, 2] = numpy.sum(x308 * x322 * x85)
    result[1, 5, 3] = numpy.sum(x268 * x341 * x407)
    result[1, 5, 4] = numpy.sum(x100 * x330 * x445)
    result[1, 5, 5] = numpy.sum(x275 * x371 * x407)
    result[1, 5, 6] = numpy.sum(x268 * x358 * x418)
    result[1, 5, 7] = numpy.sum(x277 * x341 * x419)
    result[1, 5, 8] = numpy.sum(x330 * x419 * x446)
    result[1, 5, 9] = numpy.sum(x280 * x371 * x418)
    result[1, 5, 10] = numpy.sum(x368 * x447 * x47)
    result[1, 5, 11] = numpy.sum(x358 * x449 * x47)
    result[1, 5, 12] = numpy.sum(x341 * x450 * x47)
    result[1, 5, 13] = numpy.sum(x280 * x347 * x443)
    result[1, 5, 14] = numpy.sum(x322 * x420 * x451)
    result[1, 6, 0] = numpy.sum(x452 * x455)
    result[1, 6, 1] = numpy.sum(x168 * x456 * x457)
    result[1, 6, 2] = numpy.sum(x455 * x457 * x95)
    result[1, 6, 3] = numpy.sum(x168 * x458 * x459)
    result[1, 6, 4] = numpy.sum(x36 * x360 * x456)
    result[1, 6, 5] = numpy.sum(x292 * x36 * x454)
    result[1, 6, 6] = numpy.sum(x460 * x461)
    result[1, 6, 7] = numpy.sum(x3 * x428 * x458)
    result[1, 6, 8] = numpy.sum(x138 * x3 * x462)
    result[1, 6, 9] = numpy.sum(x143 * x3 * x463)
    result[1, 6, 10] = numpy.sum(
        x425
        * x76
        * (
            x0
            * (
                x230 * (x422 + x423)
                + x29 * (x191 + x388 + 5.0 * x390 + 4.0 * x414)
                + 2.0 * x392
                + 2.0 * x393
                + 4.0 * x415
                + 4.0 * x416
            )
            + x161 * x424
        )
    )
    result[1, 6, 11] = numpy.sum(x460 * x96)
    result[1, 6, 12] = numpy.sum(x120 * x458 * x464)
    result[1, 6, 13] = numpy.sum(x143 * x462)
    result[1, 6, 14] = numpy.sum(x160 * x463)
    result[1, 7, 0] = numpy.sum(x197 * x399 * x53)
    result[1, 7, 1] = numpy.sum(x207 * x404 * x465)
    result[1, 7, 2] = numpy.sum(x300 * x39 * x398)
    result[1, 7, 3] = numpy.sum(x197 * x413 * x466)
    result[1, 7, 4] = numpy.sum(x301 * x36 * x404)
    result[1, 7, 5] = numpy.sum(x302 * x36 * x398)
    result[1, 7, 6] = numpy.sum(x417 * x442 * x467)
    result[1, 7, 7] = numpy.sum(x201 * x412 * x469)
    result[1, 7, 8] = numpy.sum(x205 * x433 * x468)
    result[1, 7, 9] = numpy.sum(x212 * x3 * x434)
    result[1, 7, 10] = numpy.sum(x198 * x427)
    result[1, 7, 11] = numpy.sum(x201 * x417 * x448)
    result[1, 7, 12] = numpy.sum(x205 * x432)
    result[1, 7, 13] = numpy.sum(x212 * x433)
    result[1, 7, 14] = numpy.sum(x220 * x434)
    result[1, 8, 0] = numpy.sum(x307 * x373 * x53)
    result[1, 8, 1] = numpy.sum(x268 * x377 * x465)
    result[1, 8, 2] = numpy.sum(x308 * x373 * x39)
    result[1, 8, 3] = numpy.sum(x268 * x380 * x466)
    result[1, 8, 4] = numpy.sum(x36 * x377 * x445)
    result[1, 8, 5] = numpy.sum(x373 * x446 * x466)
    result[1, 8, 6] = numpy.sum(x260 * x384 * x470)
    result[1, 8, 7] = numpy.sum(x270 * x380 * x471)
    result[1, 8, 8] = numpy.sum(x275 * x377 * x469)
    result[1, 8, 9] = numpy.sum(x3 * x373 * x472)
    result[1, 8, 10] = numpy.sum(x394 * x447)
    result[1, 8, 11] = numpy.sum(x384 * x449)
    result[1, 8, 12] = numpy.sum(x380 * x450)
    result[1, 8, 13] = numpy.sum(x377 * x472)
    result[1, 8, 14] = numpy.sum(x374 * x451)
    result[1, 9, 0] = numpy.sum(x314 * x322 * x452)
    result[1, 9, 1] = numpy.sum(x314 * x330 * x457)
    result[1, 9, 2] = numpy.sum(x313 * x371 * x457)
    result[1, 9, 3] = numpy.sum(x314 * x341 * x459)
    result[1, 9, 4] = numpy.sum(x313 * x347 * x473)
    result[1, 9, 5] = numpy.sum(x316 * x371 * x459)
    result[1, 9, 6] = numpy.sum(x3 * x309 * x475)
    result[1, 9, 7] = numpy.sum(x313 * x341 * x470)
    result[1, 9, 8] = numpy.sum(x316 * x347 * x471)
    result[1, 9, 9] = numpy.sum(x3 * x322 * x476)
    result[1, 9, 10] = numpy.sum(x309 * x368 * x477)
    result[1, 9, 11] = numpy.sum(x313 * x475)
    result[1, 9, 12] = numpy.sum(x316 * x341 * x464)
    result[1, 9, 13] = numpy.sum(x330 * x476)
    result[1, 9, 14] = numpy.sum(x319 * x322 * x477)
    result[2, 0, 0] = numpy.sum(x310 * x326 * x480)
    result[2, 0, 1] = numpy.sum(x480 * x481 * x81)
    result[2, 0, 2] = numpy.sum(x481 * x485)
    result[2, 0, 3] = numpy.sum(x112 * x346 * x486)
    result[2, 0, 4] = numpy.sum(x315 * x346 * x485)
    result[2, 0, 5] = numpy.sum(x123 * x346 * x491)
    result[2, 0, 6] = numpy.sum(x128 * x350 * x486)
    result[2, 0, 7] = numpy.sum(x136 * x350 * x492)
    result[2, 0, 8] = numpy.sum(x315 * x350 * x491)
    result[2, 0, 9] = numpy.sum(x350 * x501)
    result[2, 0, 10] = numpy.sum(x148 * x370 * x486)
    result[2, 0, 11] = numpy.sum(x128 * x370 * x492)
    result[2, 0, 12] = numpy.sum(x112 * x370 * x502)
    result[2, 0, 13] = numpy.sum(x370 * x501 * x81)
    result[2, 0, 14] = numpy.sum(x310 * x370 * x509)
    result[2, 1, 0] = numpy.sum(x161 * x480 * x511)
    result[2, 1, 1] = numpy.sum(x305 * x334 * x480)
    result[2, 1, 2] = numpy.sum(x161 * x485 * x512)
    result[2, 1, 3] = numpy.sum(x173 * x345 * x486)
    result[2, 1, 4] = numpy.sum(x306 * x345 * x485)
    result[2, 1, 5] = numpy.sum(x161 * x491 * x513)
    result[2, 1, 6] = numpy.sum(x181 * x349 * x486)
    result[2, 1, 7] = numpy.sum(x172 * x386 * x492)
    result[2, 1, 8] = numpy.sum(x306 * x349 * x491)
    result[2, 1, 9] = numpy.sum(x161 * x500 * x514)
    result[2, 1, 10] = numpy.sum(x192 * x397 * x480)
    result[2, 1, 11] = numpy.sum(x181 * x397 * x485)
    result[2, 1, 12] = numpy.sum(x173 * x397 * x491)
    result[2, 1, 13] = numpy.sum(x305 * x369 * x500)
    result[2, 1, 14] = numpy.sum(x161 * x509 * x515)
    result[2, 2, 0] = numpy.sum(x511 * x517)
    result[2, 2, 1] = numpy.sum(x512 * x517 * x81)
    result[2, 2, 2] = numpy.sum(x512 * x519)
    result[2, 2, 3] = numpy.sum(x203 * x345 * x520)
    result[2, 2, 4] = numpy.sum(x345 * x519 * x521)
    result[2, 2, 5] = numpy.sum(x513 * x523)
    result[2, 2, 6] = numpy.sum(x206 * x349 * x520)
    result[2, 2, 7] = numpy.sum(x111 * x386 * x524)
    result[2, 2, 8] = numpy.sum(x269 * x386 * x523)
    result[2, 2, 9] = numpy.sum(x514 * x526)
    result[2, 2, 10] = numpy.sum(x213 * x397 * x517)
    result[2, 2, 11] = numpy.sum(x206 * x397 * x519)
    result[2, 2, 12] = numpy.sum(x203 * x397 * x523)
    result[2, 2, 13] = numpy.sum(x315 * x369 * x526)
    result[2, 2, 14] = numpy.sum(x515 * x534)
    result[2, 3, 0] = numpy.sum(x223 * x486 * x55)
    result[2, 3, 1] = numpy.sum(x227 * x405 * x486)
    result[2, 3, 2] = numpy.sum(x222 * x405 * x492)
    result[2, 3, 3] = numpy.sum(x234 * x407 * x486)
    result[2, 3, 4] = numpy.sum(x100 * x237 * x492)
    result[2, 3, 5] = numpy.sum(x222 * x407 * x502)
    result[2, 3, 6] = numpy.sum(x241 * x480 * x535)
    result[2, 3, 7] = numpy.sum(x234 * x419 * x492)
    result[2, 3, 8] = numpy.sum(x227 * x419 * x502)
    result[2, 3, 9] = numpy.sum(x222 * x500 * x535)
    result[2, 3, 10] = numpy.sum(x443 * x480 * x536)
    result[2, 3, 11] = numpy.sum(x47 * x485 * x537)
    result[2, 3, 12] = numpy.sum(x47 * x491 * x538)
    result[2, 3, 13] = numpy.sum(x47 * x500 * x539)
    result[2, 3, 14] = numpy.sum(x47 * x509 * x540)
    result[2, 4, 0] = numpy.sum(35.3313566383285 * x255 * x517 * x55)
    result[2, 4, 1] = numpy.sum(x195 * x541 * x85)
    result[2, 4, 2] = numpy.sum(x519 * x542 * x85)
    result[2, 4, 3] = numpy.sum(x100 * x250 * x520)
    result[2, 4, 4] = numpy.sum(x100 * x251 * x519)
    result[2, 4, 5] = numpy.sum(x100 * x176 * x252 * x523)
    result[2, 4, 6] = numpy.sum(x180 * x440 * x541)
    result[2, 4, 7] = numpy.sum(x172 * x439 * x524)
    result[2, 4, 8] = numpy.sum(x251 * x343 * x523)
    result[2, 4, 9] = numpy.sum(x343 * x526 * x542)
    result[2, 4, 10] = numpy.sum(x191 * x441 * x543)
    result[2, 4, 11] = numpy.sum(x180 * x444 * x519)
    result[2, 4, 12] = numpy.sum(x250 * x443 * x523)
    result[2, 4, 13] = numpy.sum(x167 * x444 * x526)
    result[2, 4, 14] = numpy.sum(x441 * x534 * x545)
    result[2, 5, 0] = numpy.sum(x547 * x55)
    result[2, 5, 1] = numpy.sum(x546 * x548 * x81)
    result[2, 5, 2] = numpy.sum(x548 * x553)
    result[2, 5, 3] = numpy.sum(x100 * x153 * x203 * x546)
    result[2, 5, 4] = numpy.sum(x100 * x521 * x553)
    result[2, 5, 5] = numpy.sum(x407 * x559)
    result[2, 5, 6] = numpy.sum(x206 * x440 * x546)
    result[2, 5, 7] = numpy.sum(x111 * x153 * x419 * x553)
    result[2, 5, 8] = numpy.sum(x419 * x559 * x81)
    result[2, 5, 9] = numpy.sum(x176 * x418 * x563)
    result[2, 5, 10] = numpy.sum(x213 * x47 * x564)
    result[2, 5, 11] = numpy.sum(x206 * x47 * x565)
    result[2, 5, 12] = numpy.sum(x203 * x47 * x566)
    result[2, 5, 13] = numpy.sum(53.9695387335403 * x47 * x544 * x563 * x81)
    result[2, 5, 14] = numpy.sum(x420 * x571)
    result[2, 6, 0] = numpy.sum(x296 * x452 * x480)
    result[2, 6, 1] = numpy.sum(x286 * x457 * x486)
    result[2, 6, 2] = numpy.sum(x296 * x457 * x485)
    result[2, 6, 3] = numpy.sum(x289 * x459 * x486)
    result[2, 6, 4] = numpy.sum(x291 * x473 * x485)
    result[2, 6, 5] = numpy.sum(x296 * x459 * x491)
    result[2, 6, 6] = numpy.sum(x3 * x480 * x572)
    result[2, 6, 7] = numpy.sum(x289 * x470 * x485)
    result[2, 6, 8] = numpy.sum(x291 * x471 * x491)
    result[2, 6, 9] = numpy.sum(x284 * x3 * x573)
    result[2, 6, 10] = numpy.sum(x297 * x477 * x480)
    result[2, 6, 11] = numpy.sum(x485 * x572)
    result[2, 6, 12] = numpy.sum(x289 * x464 * x491)
    result[2, 6, 13] = numpy.sum(x286 * x573)
    result[2, 6, 14] = numpy.sum(x284 * x477 * x509)
    result[2, 7, 0] = numpy.sum(x223 * x520 * x53)
    result[2, 7, 1] = numpy.sum(x227 * x465 * x520)
    result[2, 7, 2] = numpy.sum(x222 * x465 * x524)
    result[2, 7, 3] = numpy.sum(x234 * x466 * x520)
    result[2, 7, 4] = numpy.sum(x237 * x473 * x519)
    result[2, 7, 5] = numpy.sum(x153 * x222 * x466 * x523)
    result[2, 7, 6] = numpy.sum(x3 * x517 * x537)
    result[2, 7, 7] = numpy.sum(x234 * x469 * x519)
    result[2, 7, 8] = numpy.sum(x237 * x471 * x523)
    result[2, 7, 9] = numpy.sum(x222 * x470 * x526)
    result[2, 7, 10] = numpy.sum(x536 * x543)
    result[2, 7, 11] = numpy.sum(x519 * x537)
    result[2, 7, 12] = numpy.sum(x523 * x538)
    result[2, 7, 13] = numpy.sum(x526 * x539)
    result[2, 7, 14] = numpy.sum(x534 * x540)
    result[2, 8, 0] = numpy.sum(x161 * x53 * x547)
    result[2, 8, 1] = numpy.sum(x305 * x39 * x546)
    result[2, 8, 2] = numpy.sum(x255 * x465 * x553)
    result[2, 8, 3] = numpy.sum(x173 * x473 * x546)
    result[2, 8, 4] = numpy.sum(x306 * x36 * x553)
    result[2, 8, 5] = numpy.sum(x161 * x466 * x559)
    result[2, 8, 6] = numpy.sum(x181 * x3 * x564)
    result[2, 8, 7] = numpy.sum(x172 * x468 * x565)
    result[2, 8, 8] = numpy.sum(x167 * x468 * x566)
    result[2, 8, 9] = numpy.sum(x467 * x545 * x563)
    result[2, 8, 10] = numpy.sum(x192 * x564)
    result[2, 8, 11] = numpy.sum(x181 * x565)
    result[2, 8, 12] = numpy.sum(x173 * x566)
    result[2, 8, 13] = numpy.sum(x167 * x448 * x563)
    result[2, 8, 14] = numpy.sum(x162 * x571)
    result[2, 9, 0] = numpy.sum(x452 * x576)
    result[2, 9, 1] = numpy.sum(x457 * x576 * x81)
    result[2, 9, 2] = numpy.sum(x176 * x457 * x577)
    result[2, 9, 3] = numpy.sum(x112 * x473 * x575)
    result[2, 9, 4] = numpy.sum(x315 * x36 * x577)
    result[2, 9, 5] = numpy.sum(x176 * x459 * x578)
    result[2, 9, 6] = numpy.sum(x128 * x3 * x579)
    result[2, 9, 7] = numpy.sum(x136 * x3 * x580)
    result[2, 9, 8] = numpy.sum(x116 * x467 * x581)
    result[2, 9, 9] = numpy.sum(x461 * x582)
    result[2, 9, 10] = numpy.sum(x148 * x579)
    result[2, 9, 11] = numpy.sum(x128 * x580)
    result[2, 9, 12] = numpy.sum(x114 * x581)
    result[2, 9, 13] = numpy.sum(x582 * x82)
    result[2, 9, 14] = numpy.sum(
        9.12251705727742
        * x544
        * (
            x0
            * (
                x271 * (x568 + x569)
                + x29 * (x219 + x528 + 5.0 * x530 + 4.0 * x560)
                + 2.0 * x532
                + 2.0 * x533
                + 4.0 * x561
                + 4.0 * x562
            )
            + x197 * x570
        )
    )
    return result


def diag_quadrupole3d_40(ax, da, A, bx, db, B, C):
    """Cartesian 3D (gs) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 15, 1), dtype=float)

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
    x12 = 2.0 * x3
    x13 = -x2 - C[0]
    x14 = x13 * x8
    x15 = x11 + x12 * x14
    x16 = 2.0 * x0
    x17 = x3 * x8
    x18 = x0 * (x14 + x17)
    x19 = x3 * (x10 + x13 * x17)
    x20 = x13**2 * x8
    x21 = x0 * (x15 + x20)
    x22 = x10 + x20
    x23 = x22 * x3
    x24 = x14 * x16 + x23
    x25 = x24 * x3
    x26 = x21 + x25
    x27 = 2.0 * x0 * (2.0 * x0 * x14 + x18 + x19 + x23) + x26 * x3
    x28 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x29 = 4.41641957979107 * x28
    x30 = 0.564189583547756
    x31 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x32 = da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**1.5)
    x33 = x31 * x32
    x34 = x1 * x30 * x33
    x35 = -x1 * (ax * A[1] + bx * B[1])
    x36 = -x35 - A[1]
    x37 = 11.6847478934435 * x36
    x38 = x28 * x37
    x39 = x27 * x34
    x40 = -x1 * (ax * A[2] + bx * B[2])
    x41 = -x40 - A[2]
    x42 = 11.6847478934435 * x41
    x43 = x28 * x7
    x44 = x36**2 * x43
    x45 = x0 * x43
    x46 = x44 + x45
    x47 = 15.084944665313 * x46
    x48 = 0.318309886183791 * x6
    x49 = x26 * x48
    x50 = 26.1278905896872 * x41
    x51 = x28 * x32
    x52 = x31 * x7
    x53 = x41**2 * x52
    x54 = x0 * x52
    x55 = x53 + x54
    x56 = 15.084944665313 * x55
    x57 = x36 * x43
    x58 = x16 * x57 + x36 * x46
    x59 = x33 * x48
    x60 = 11.6847478934435 * x24
    x61 = 26.1278905896872 * x24
    x62 = x48 * x51
    x63 = x41 * x52
    x64 = x16 * x63 + x41 * x55
    x65 = x48 * x64
    x66 = 3.0 * x45
    x67 = 4.41641957979107 * x0 * (3.0 * x44 + x66) + 4.41641957979107 * x36 * x58
    x68 = x22 * x32
    x69 = x48 * x68
    x70 = x31 * x69
    x71 = 0.179587122125167
    x72 = x55 * x71
    x73 = 3.0 * x54
    x74 = x0 * (3.0 * x53 + x73) + x41 * x64
    x75 = x10 + x9
    x76 = x16 * x17 + x3 * x75
    x77 = x0 * (x11 + 3.0 * x9) + x3 * x76
    x78 = -x35 - C[1]
    x79 = x43 * x78**2
    x80 = x45 + x79
    x81 = 4.41641957979107 * x80
    x82 = x36 * x80
    x83 = x43 * x78
    x84 = x16 * x83 + x82
    x85 = 11.6847478934435 * x84
    x86 = x59 * x76
    x87 = x32 * x75
    x88 = 15.084944665313 * x87
    x89 = 2.0 * x36
    x90 = x66 + x83 * x89
    x91 = x0 * (x79 + x90)
    x92 = x36 * x84
    x93 = x91 + x92
    x94 = x48 * x93
    x95 = 26.1278905896872 * x84
    x96 = 11.6847478934435 * x3
    x97 = x0 * (x57 + x83)
    x98 = x36 * (x45 + x57 * x78)
    x99 = 2.0 * x0 * (2.0 * x45 * x78 + x82 + x97 + x98) + x36 * x93
    x100 = x34 * x5
    x101 = x100 * x99
    x102 = x32 * x5
    x103 = x102 * x48
    x104 = x103 * x3
    x105 = x102 * x65
    x106 = -x40 - C[2]
    x107 = x106**2 * x52
    x108 = x107 + x54
    x109 = x108 * x48
    x110 = x37 * x51
    x111 = x108 * x41
    x112 = x106 * x52
    x113 = x111 + x112 * x16
    x114 = 11.6847478934435 * x113
    x115 = x28 * x48
    x116 = 26.1278905896872 * x113
    x117 = 2.0 * x41
    x118 = x112 * x117 + x73
    x119 = x0 * (x107 + x118)
    x120 = x113 * x41
    x121 = x119 + x120
    x122 = x103 * x108
    x123 = x1 * x30
    x124 = x0 * (x112 + x63)
    x125 = x41 * (x106 * x63 + x54)
    x126 = 2.0 * x0 * (2.0 * x106 * x54 + x111 + x124 + x125) + x121 * x41
    x127 = x123 * x126 * x5

    # 45 item(s)
    result[0, 0, 0] = numpy.sum(
        x29
        * x34
        * (x0 * (x12 * (x18 + x19) + x16 * (x15 + x9) + 3.0 * x21 + 3.0 * x25) + x27 * x3)
    )
    result[0, 1, 0] = numpy.sum(x38 * x39)
    result[0, 2, 0] = numpy.sum(x28 * x39 * x42)
    result[0, 3, 0] = numpy.sum(x33 * x47 * x49)
    result[0, 4, 0] = numpy.sum(x26 * x28 * x34 * x36 * x50)
    result[0, 5, 0] = numpy.sum(x49 * x51 * x56)
    result[0, 6, 0] = numpy.sum(x58 * x59 * x60)
    result[0, 7, 0] = numpy.sum(x41 * x46 * x59 * x61)
    result[0, 8, 0] = numpy.sum(x36 * x55 * x61 * x62)
    result[0, 9, 0] = numpy.sum(x51 * x60 * x65)
    result[0, 10, 0] = numpy.sum(x67 * x70)
    result[0, 11, 0] = numpy.sum(x42 * x58 * x70)
    result[0, 12, 0] = numpy.sum(x47 * x68 * x72)
    result[0, 13, 0] = numpy.sum(x38 * x65 * x68)
    result[0, 14, 0] = numpy.sum(x29 * x69 * x74)
    result[1, 0, 0] = numpy.sum(x59 * x77 * x81)
    result[1, 1, 0] = numpy.sum(x85 * x86)
    result[1, 2, 0] = numpy.sum(x42 * x80 * x86)
    result[1, 3, 0] = numpy.sum(x31 * x88 * x94)
    result[1, 4, 0] = numpy.sum(x31 * x41 * x48 * x87 * x95)
    result[1, 5, 0] = numpy.sum(x72 * x80 * x88)
    result[1, 6, 0] = numpy.sum(x101 * x96)
    result[1, 7, 0] = numpy.sum(x100 * x3 * x50 * x93)
    result[1, 8, 0] = numpy.sum(x104 * x55 * x95)
    result[1, 9, 0] = numpy.sum(x105 * x80 * x96)
    result[1, 10, 0] = numpy.sum(
        4.41641957979107
        * x100
        * (
            x0 * (x16 * (x44 + x90) + x89 * (x97 + x98) + 3.0 * x91 + 3.0 * x92)
            + x36 * x99
        )
    )
    result[1, 11, 0] = numpy.sum(x101 * x42)
    result[1, 12, 0] = numpy.sum(x102 * x56 * x94)
    result[1, 13, 0] = numpy.sum(x105 * x85)
    result[1, 14, 0] = numpy.sum(x103 * x74 * x81)
    result[2, 0, 0] = numpy.sum(x109 * x29 * x32 * x77)
    result[2, 1, 0] = numpy.sum(x109 * x110 * x76)
    result[2, 2, 0] = numpy.sum(x114 * x62 * x76)
    result[2, 3, 0] = numpy.sum(x108 * x47 * x71 * x87)
    result[2, 4, 0] = numpy.sum(x115 * x116 * x36 * x87)
    result[2, 5, 0] = numpy.sum(x115 * x121 * x88)
    result[2, 6, 0] = numpy.sum(x122 * x58 * x96)
    result[2, 7, 0] = numpy.sum(x104 * x116 * x46)
    result[2, 8, 0] = numpy.sum(26.1278905896872 * x121 * x123 * x3 * x36 * x5 * x51)
    result[2, 9, 0] = numpy.sum(x127 * x51 * x96)
    result[2, 10, 0] = numpy.sum(x122 * x67)
    result[2, 11, 0] = numpy.sum(x103 * x114 * x58)
    result[2, 12, 0] = numpy.sum(x103 * x121 * x47)
    result[2, 13, 0] = numpy.sum(x110 * x127)
    result[2, 14, 0] = numpy.sum(
        x102
        * x123
        * x29
        * (
            x0 * (x117 * (x124 + x125) + 3.0 * x119 + 3.0 * x120 + x16 * (x118 + x53))
            + x126 * x41
        )
    )
    return result


def diag_quadrupole3d_41(ax, da, A, bx, db, B, C):
    """Cartesian 3D (gp) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 15, 3), dtype=float)

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
    x10 = -x4 - C[0]
    x11 = x3 * x7
    x12 = x10 * x11
    x13 = x0 * (x12 + x9)
    x14 = x0 * x11
    x15 = x10 * x9
    x16 = x5 * (x14 + x15)
    x17 = x13 + x16
    x18 = -x4 - B[0]
    x19 = x11 * x18
    x20 = x0 * (x19 + x9)
    x21 = x18 * x9
    x22 = x14 + x21
    x23 = x22 * x5
    x24 = x20 + x23
    x25 = x0 * (x12 + x19)
    x26 = x12 * x18
    x27 = x5 * (x14 + x26)
    x28 = 2.0 * x25 + 2.0 * x27
    x29 = 2.0 * x0
    x30 = 3.0 * x14
    x31 = x0 * (x15 + x21 + x26 + x30)
    x32 = x5 * (x25 + x27)
    x33 = 2.0 * x5
    x34 = x10**2 * x11
    x35 = x14 + x34
    x36 = x35 * x5
    x37 = x18 * x35
    x38 = 4.0 * x0 * x12
    x39 = x0 * (x28 + x36 + x37 + x38)
    x40 = x30 + x34
    x41 = x0 * (2.0 * x26 + x40)
    x42 = x12 * x29
    x43 = x37 + x42
    x44 = x43 * x5
    x45 = x41 + x44
    x46 = x45 * x5
    x47 = x12 * x33
    x48 = x0 * (x40 + x47)
    x49 = x36 + x42
    x50 = x49 * x5
    x51 = x48 + x50
    x52 = x0 * (2.0 * x13 + 2.0 * x16 + 2.0 * x36 + x38) + x5 * x51
    x53 = x39 + x46
    x54 = x0 * (2.0 * x31 + 2.0 * x32 + 2.0 * x41 + 2.0 * x44 + x51) + x5 * x53
    x55 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x56 = da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**2.5)
    x57 = x55 * x56
    x58 = 8.83283915958214 * x57
    x59 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x60 = 0.564189583547756 * x1
    x61 = x59 * x60
    x62 = x58 * x61
    x63 = -x1 * (ax * A[1] + bx * B[1])
    x64 = -x63 - B[1]
    x65 = x11 * x5**2
    x66 = x30 + x65
    x67 = x62 * (x0 * (x17 * x33 + x29 * (x47 + x66) + 3.0 * x48 + 3.0 * x50) + x5 * x52)
    x68 = -x1 * (ax * A[2] + bx * B[2])
    x69 = -x68 - B[2]
    x70 = -x63 - A[1]
    x71 = x61 * x70
    x72 = 23.3694957868871 * x57
    x73 = x54 * x72
    x74 = x52 * x57
    x75 = x0 * x3
    x76 = x59 * x75
    x77 = x3 * x59
    x78 = x70 * x77
    x79 = x64 * x78
    x80 = x76 + x79
    x81 = 23.3694957868871 * x80
    x82 = 0.318309886183791 * x2
    x83 = x81 * x82
    x84 = 23.3694957868871 * x74
    x85 = -x68 - A[2]
    x86 = x61 * x85
    x87 = x55 * x75
    x88 = x3 * x55
    x89 = x85 * x88
    x90 = x69 * x89
    x91 = x87 + x90
    x92 = 23.3694957868871 * x82
    x93 = x56 * x59
    x94 = x92 * x93
    x95 = x53 * x57
    x96 = x70**2 * x77
    x97 = x76 + x96
    x98 = 30.169889330626 * x82
    x99 = x97 * x98
    x100 = x51 * x98
    x101 = x64 * x77
    x102 = x0 * (x101 + x78)
    x103 = x70 * x80
    x104 = x102 + x103
    x105 = x104 * x57
    x106 = x51 * x57
    x107 = 52.2557811793745 * x85
    x108 = 52.2557811793745 * x80
    x109 = x82 * x85
    x110 = x70 * x93
    x111 = 52.2557811793745 * x82
    x112 = x111 * x91
    x113 = x85**2 * x88
    x114 = x113 + x87
    x115 = x114 * x98
    x116 = x100 * x93
    x117 = x69 * x88
    x118 = x0 * (x117 + x89)
    x119 = x85 * x91
    x120 = x118 + x119
    x121 = x29 * x78 + x70 * x97
    x122 = 23.3694957868871 * x121
    x123 = x122 * x82
    x124 = x123 * x57
    x125 = 2.0 * x70
    x126 = 3.0 * x76
    x127 = x126 + x96
    x128 = x0 * (x101 * x125 + x127) + x104 * x70
    x129 = x57 * x92
    x130 = 52.2557811793745 * x97
    x131 = x109 * x57
    x132 = 0.179587122125167 * x56
    x133 = x132 * x91
    x134 = x110 * x111
    x135 = x114 * x132
    x136 = x114 * x85 + x29 * x89
    x137 = x136 * x94
    x138 = 2.0 * x85
    x139 = 3.0 * x87
    x140 = x113 + x139
    x141 = x0 * (x117 * x138 + x140) + x120 * x85
    x142 = x0 * (x126 + 3.0 * x96) + x121 * x70
    x143 = x58 * x82
    x144 = x142 * x143
    x145 = x0 * (3.0 * x102 + 3.0 * x103 + x121) + x128 * x70
    x146 = 23.3694957868871 * x35
    x147 = 30.169889330626 * x135
    x148 = 30.169889330626 * x132
    x149 = x120 * x148
    x150 = x132 * x136
    x151 = x141 * x82
    x152 = x0 * (3.0 * x113 + x139) + x136 * x85
    x153 = 8.83283915958214 * x82
    x154 = x153 * x93
    x155 = x152 * x154
    x156 = x0 * (3.0 * x118 + 3.0 * x119 + x136) + x141 * x85
    x157 = x14 + x65
    x158 = x157 * x5 + x29 * x9
    x159 = x0 * (x19 * x33 + x66) + x24 * x5
    x160 = x0 * (x158 + 3.0 * x20 + 3.0 * x23) + x159 * x5
    x161 = -x63 - C[1]
    x162 = x161**2 * x77
    x163 = x162 + x76
    x164 = x143 * x163
    x165 = x163 * x64
    x166 = x161 * x77
    x167 = x166 * x29
    x168 = x165 + x167
    x169 = x0 * (x30 + 3.0 * x65) + x158 * x5
    x170 = x163 * x70
    x171 = x167 + x170
    x172 = x129 * x171
    x173 = x101 * x161
    x174 = x126 + x162
    x175 = x0 * (2.0 * x173 + x174)
    x176 = x168 * x70
    x177 = x175 + x176
    x178 = x129 * x158
    x179 = 23.3694957868871 * x163
    x180 = x125 * x166
    x181 = x0 * (x174 + x180)
    x182 = x171 * x70
    x183 = x181 + x182
    x184 = x57 * x98
    x185 = x183 * x184
    x186 = 4.0 * x161 * x76
    x187 = x0 * (x101 + x166)
    x188 = x70 * (x173 + x76)
    x189 = 2.0 * x187 + 2.0 * x188
    x190 = x0 * (x165 + x170 + x186 + x189)
    x191 = x177 * x70
    x192 = x190 + x191
    x193 = 52.2557811793745 * x171
    x194 = 52.2557811793745 * x131
    x195 = x0 * (x166 + x78)
    x196 = x161 * x78
    x197 = x70 * (x196 + x76)
    x198 = x0 * (2.0 * x170 + x186 + 2.0 * x195 + 2.0 * x197) + x183 * x70
    x199 = x0 * (x126 + x173 + x196 + x79)
    x200 = x70 * (x187 + x188)
    x201 = x0 * (2.0 * x175 + 2.0 * x176 + x183 + 2.0 * x199 + 2.0 * x200) + x192 * x70
    x202 = x201 * x72
    x203 = x60 * x8
    x204 = x198 * x72
    x205 = x56 * x8
    x206 = x111 * x205
    x207 = x205 * x82
    x208 = x205 * x92
    x209 = x195 + x197
    x210 = x60 * x7
    x211 = x210 * x58
    x212 = x211 * (
        x0 * (x125 * x209 + 3.0 * x181 + 3.0 * x182 + x29 * (x127 + x180)) + x198 * x70
    )
    x213 = x210 * x85
    x214 = x56 * x7
    x215 = x214 * x92
    x216 = x115 * x214
    x217 = x214 * x98
    x218 = x136 * x215
    x219 = x153 * x214
    x220 = x152 * x219
    x221 = -x68 - C[2]
    x222 = x221**2 * x88
    x223 = x222 + x87
    x224 = x154 * x223
    x225 = x223 * x69
    x226 = x221 * x88
    x227 = x226 * x29
    x228 = x225 + x227
    x229 = x159 * x94
    x230 = x132 * x223
    x231 = x158 * x94
    x232 = x223 * x85
    x233 = x227 + x232
    x234 = x117 * x221
    x235 = x139 + x222
    x236 = x0 * (2.0 * x234 + x235)
    x237 = x228 * x85
    x238 = x236 + x237
    x239 = x148 * x223
    x240 = x132 * x233
    x241 = x138 * x226
    x242 = x0 * (x235 + x241)
    x243 = x233 * x85
    x244 = x242 + x243
    x245 = x93 * x98
    x246 = x244 * x245
    x247 = 4.0 * x221 * x87
    x248 = x0 * (x117 + x226)
    x249 = x85 * (x234 + x87)
    x250 = 2.0 * x248 + 2.0 * x249
    x251 = x0 * (x225 + x232 + x247 + x250)
    x252 = x238 * x85
    x253 = x251 + x252
    x254 = x0 * (x226 + x89)
    x255 = x221 * x89
    x256 = x85 * (x255 + x87)
    x257 = x0 * (2.0 * x232 + x247 + 2.0 * x254 + 2.0 * x256) + x244 * x85
    x258 = x205 * x61
    x259 = x0 * (x139 + x234 + x255 + x90)
    x260 = x85 * (x248 + x249)
    x261 = x0 * (2.0 * x236 + 2.0 * x237 + x244 + 2.0 * x259 + 2.0 * x260) + x253 * x85
    x262 = 23.3694957868871 * x261
    x263 = x219 * x223
    x264 = x123 * x214
    x265 = x214 * x99
    x266 = x214 * x257
    x267 = x214 * x61
    x268 = x254 + x256
    x269 = 8.83283915958214 * x267
    x270 = x269 * (
        x0 * (x138 * x268 + 3.0 * x242 + 3.0 * x243 + x29 * (x140 + x241)) + x257 * x85
    )

    # 135 item(s)
    result[0, 0, 0] = numpy.sum(
        x62
        * (
            x0
            * (x29 * (x17 + x24 + x28) + x33 * (x31 + x32) + 3.0 * x39 + 3.0 * x46 + x52)
            + x5 * x54
        )
    )
    result[0, 0, 1] = numpy.sum(x64 * x67)
    result[0, 0, 2] = numpy.sum(x67 * x69)
    result[0, 1, 0] = numpy.sum(x71 * x73)
    result[0, 1, 1] = numpy.sum(x74 * x83)
    result[0, 1, 2] = numpy.sum(x69 * x71 * x84)
    result[0, 2, 0] = numpy.sum(x73 * x86)
    result[0, 2, 1] = numpy.sum(x64 * x84 * x86)
    result[0, 2, 2] = numpy.sum(x52 * x91 * x94)
    result[0, 3, 0] = numpy.sum(x95 * x99)
    result[0, 3, 1] = numpy.sum(x100 * x105)
    result[0, 3, 2] = numpy.sum(x106 * x69 * x99)
    result[0, 4, 0] = numpy.sum(x107 * x71 * x95)
    result[0, 4, 1] = numpy.sum(x106 * x108 * x109)
    result[0, 4, 2] = numpy.sum(x110 * x112 * x51)
    result[0, 5, 0] = numpy.sum(x115 * x53 * x93)
    result[0, 5, 1] = numpy.sum(x114 * x116 * x64)
    result[0, 5, 2] = numpy.sum(x116 * x120)
    result[0, 6, 0] = numpy.sum(x124 * x45)
    result[0, 6, 1] = numpy.sum(x128 * x129 * x49)
    result[0, 6, 2] = numpy.sum(x124 * x49 * x69)
    result[0, 7, 0] = numpy.sum(x130 * x131 * x45)
    result[0, 7, 1] = numpy.sum(52.2557811793745 * x105 * x109 * x49)
    result[0, 7, 2] = numpy.sum(x130 * x133 * x49)
    result[0, 8, 0] = numpy.sum(x114 * x134 * x45)
    result[0, 8, 1] = numpy.sum(x108 * x135 * x49)
    result[0, 8, 2] = numpy.sum(x120 * x134 * x49)
    result[0, 9, 0] = numpy.sum(x137 * x45)
    result[0, 9, 1] = numpy.sum(x137 * x49 * x64)
    result[0, 9, 2] = numpy.sum(x141 * x49 * x94)
    result[0, 10, 0] = numpy.sum(x144 * x43)
    result[0, 10, 1] = numpy.sum(x143 * x145 * x35)
    result[0, 10, 2] = numpy.sum(x144 * x35 * x69)
    result[0, 11, 0] = numpy.sum(x124 * x43 * x85)
    result[0, 11, 1] = numpy.sum(x128 * x131 * x146)
    result[0, 11, 2] = numpy.sum(x121 * x133 * x146)
    result[0, 12, 0] = numpy.sum(x147 * x43 * x97)
    result[0, 12, 1] = numpy.sum(x104 * x147 * x35)
    result[0, 12, 2] = numpy.sum(x149 * x35 * x97)
    result[0, 13, 0] = numpy.sum(x137 * x43 * x70)
    result[0, 13, 1] = numpy.sum(x146 * x150 * x80)
    result[0, 13, 2] = numpy.sum(x110 * x146 * x151)
    result[0, 14, 0] = numpy.sum(x155 * x43)
    result[0, 14, 1] = numpy.sum(x155 * x35 * x64)
    result[0, 14, 2] = numpy.sum(x154 * x156 * x35)
    result[1, 0, 0] = numpy.sum(x160 * x164)
    result[1, 0, 1] = numpy.sum(x143 * x168 * x169)
    result[1, 0, 2] = numpy.sum(x164 * x169 * x69)
    result[1, 1, 0] = numpy.sum(x159 * x172)
    result[1, 1, 1] = numpy.sum(x177 * x178)
    result[1, 1, 2] = numpy.sum(x158 * x172 * x69)
    result[1, 2, 0] = numpy.sum(x131 * x159 * x179)
    result[1, 2, 1] = numpy.sum(x168 * x178 * x85)
    result[1, 2, 2] = numpy.sum(x133 * x158 * x179)
    result[1, 3, 0] = numpy.sum(x185 * x24)
    result[1, 3, 1] = numpy.sum(x157 * x184 * x192)
    result[1, 3, 2] = numpy.sum(x157 * x185 * x69)
    result[1, 4, 0] = numpy.sum(x131 * x193 * x24)
    result[1, 4, 1] = numpy.sum(x157 * x177 * x194)
    result[1, 4, 2] = numpy.sum(x133 * x157 * x193)
    result[1, 5, 0] = numpy.sum(x147 * x163 * x24)
    result[1, 5, 1] = numpy.sum(x147 * x157 * x168)
    result[1, 5, 2] = numpy.sum(x149 * x157 * x163)
    result[1, 6, 0] = numpy.sum(x129 * x198 * x22)
    result[1, 6, 1] = numpy.sum(x202 * x203)
    result[1, 6, 2] = numpy.sum(x203 * x204 * x69)
    result[1, 7, 0] = numpy.sum(x183 * x194 * x22)
    result[1, 7, 1] = numpy.sum(x107 * x192 * x203 * x57)
    result[1, 7, 2] = numpy.sum(x112 * x183 * x205)
    result[1, 8, 0] = numpy.sum(x135 * x193 * x22)
    result[1, 8, 1] = numpy.sum(x114 * x177 * x206)
    result[1, 8, 2] = numpy.sum(x120 * x193 * x207)
    result[1, 9, 0] = numpy.sum(x150 * x179 * x22)
    result[1, 9, 1] = numpy.sum(x136 * x168 * x208)
    result[1, 9, 2] = numpy.sum(x151 * x179 * x205)
    result[1, 10, 0] = numpy.sum(x18 * x212)
    result[1, 10, 1] = numpy.sum(
        x211
        * (
            x0
            * (
                x125 * (x199 + x200)
                + 3.0 * x190
                + 3.0 * x191
                + x198
                + x29 * (x104 + x189 + x209)
            )
            + x201 * x70
        )
    )
    result[1, 10, 2] = numpy.sum(x212 * x69)
    result[1, 11, 0] = numpy.sum(x18 * x204 * x213)
    result[1, 11, 1] = numpy.sum(x202 * x213)
    result[1, 11, 2] = numpy.sum(x198 * x215 * x91)
    result[1, 12, 0] = numpy.sum(x18 * x183 * x216)
    result[1, 12, 1] = numpy.sum(x192 * x216)
    result[1, 12, 2] = numpy.sum(x120 * x183 * x217)
    result[1, 13, 0] = numpy.sum(x171 * x18 * x218)
    result[1, 13, 1] = numpy.sum(x177 * x218)
    result[1, 13, 2] = numpy.sum(x141 * x171 * x215)
    result[1, 14, 0] = numpy.sum(x163 * x18 * x220)
    result[1, 14, 1] = numpy.sum(x168 * x220)
    result[1, 14, 2] = numpy.sum(x156 * x163 * x219)
    result[2, 0, 0] = numpy.sum(x160 * x224)
    result[2, 0, 1] = numpy.sum(x169 * x224 * x64)
    result[2, 0, 2] = numpy.sum(x154 * x169 * x228)
    result[2, 1, 0] = numpy.sum(x223 * x229 * x70)
    result[2, 1, 1] = numpy.sum(x158 * x230 * x81)
    result[2, 1, 2] = numpy.sum(x228 * x231 * x70)
    result[2, 2, 0] = numpy.sum(x229 * x233)
    result[2, 2, 1] = numpy.sum(x231 * x233 * x64)
    result[2, 2, 2] = numpy.sum(x231 * x238)
    result[2, 3, 0] = numpy.sum(x239 * x24 * x97)
    result[2, 3, 1] = numpy.sum(x104 * x157 * x239)
    result[2, 3, 2] = numpy.sum(x148 * x157 * x228 * x97)
    result[2, 4, 0] = numpy.sum(x134 * x233 * x24)
    result[2, 4, 1] = numpy.sum(x108 * x157 * x240)
    result[2, 4, 2] = numpy.sum(x134 * x157 * x238)
    result[2, 5, 0] = numpy.sum(x24 * x246)
    result[2, 5, 1] = numpy.sum(x157 * x246 * x64)
    result[2, 5, 2] = numpy.sum(x157 * x245 * x253)
    result[2, 6, 0] = numpy.sum(x122 * x22 * x230)
    result[2, 6, 1] = numpy.sum(x128 * x208 * x223)
    result[2, 6, 2] = numpy.sum(x123 * x205 * x228)
    result[2, 7, 0] = numpy.sum(x130 * x22 * x240)
    result[2, 7, 1] = numpy.sum(x104 * x206 * x233)
    result[2, 7, 2] = numpy.sum(x130 * x207 * x238)
    result[2, 8, 0] = numpy.sum(x134 * x22 * x244)
    result[2, 8, 1] = numpy.sum(x108 * x207 * x244)
    result[2, 8, 2] = numpy.sum(52.2557811793745 * x205 * x253 * x71)
    result[2, 9, 0] = numpy.sum(x22 * x257 * x94)
    result[2, 9, 1] = numpy.sum(23.3694957868871 * x257 * x258 * x64)
    result[2, 9, 2] = numpy.sum(x258 * x262)
    result[2, 10, 0] = numpy.sum(x142 * x18 * x263)
    result[2, 10, 1] = numpy.sum(x145 * x263)
    result[2, 10, 2] = numpy.sum(x142 * x219 * x228)
    result[2, 11, 0] = numpy.sum(x18 * x233 * x264)
    result[2, 11, 1] = numpy.sum(x128 * x215 * x233)
    result[2, 11, 2] = numpy.sum(x238 * x264)
    result[2, 12, 0] = numpy.sum(x18 * x244 * x265)
    result[2, 12, 1] = numpy.sum(x104 * x217 * x244)
    result[2, 12, 2] = numpy.sum(x253 * x265)
    result[2, 13, 0] = numpy.sum(23.3694957868871 * x18 * x266 * x71)
    result[2, 13, 1] = numpy.sum(x266 * x83)
    result[2, 13, 2] = numpy.sum(x262 * x267 * x70)
    result[2, 14, 0] = numpy.sum(x18 * x270)
    result[2, 14, 1] = numpy.sum(x270 * x64)
    result[2, 14, 2] = numpy.sum(
        x269
        * (
            x0
            * (
                x138 * (x259 + x260)
                + 3.0 * x251
                + 3.0 * x252
                + x257
                + x29 * (x120 + x250 + x268)
            )
            + x261 * x85
        )
    )
    return result


def diag_quadrupole3d_42(ax, da, A, bx, db, B, C):
    """Cartesian 3D (gd) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 15, 6), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - A[0]
    x4 = 2.0 * x3
    x5 = -x2 - B[0]
    x6 = ax * bx * x1
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = numpy.sqrt(x1)
    x9 = 1.77245385090552 * x8
    x10 = x7 * x9
    x11 = x10 * x5
    x12 = x11 * x4
    x13 = x10 * x5**2
    x14 = x0 * x10
    x15 = 3.0 * x14
    x16 = x13 + x15
    x17 = x0 * (x12 + x16)
    x18 = x13 + x14
    x19 = x18 * x3
    x20 = 2.0 * x0
    x21 = x11 * x20 + x19
    x22 = x21 * x3
    x23 = x17 + x22
    x24 = -x2 - C[0]
    x25 = x11 * x24
    x26 = 2.0 * x25
    x27 = x0 * (x16 + x26)
    x28 = x10 * x24
    x29 = x0 * (x11 + x28)
    x30 = x14 + x25
    x31 = x30 * x5
    x32 = x3 * (x29 + x31)
    x33 = 2.0 * x27 + 2.0 * x32
    x34 = x11 * x3
    x35 = x28 * x3
    x36 = x0 * (x15 + x25 + x34 + x35)
    x37 = x3 * x30
    x38 = x3 * (x29 + x37)
    x39 = 2.0 * x36 + 2.0 * x38
    x40 = x10 * x24**2
    x41 = x14 + x40
    x42 = x41 * x5
    x43 = x20 * x28
    x44 = x42 + x43
    x45 = x3 * x44
    x46 = 2.0 * x45
    x47 = x15 + x40
    x48 = x0 * (x26 + x47)
    x49 = x28 * x4
    x50 = x0 * (x47 + x49)
    x51 = x3 * x41
    x52 = x43 + x51
    x53 = x3 * x52
    x54 = x50 + x53
    x55 = x0 * (x39 + x46 + 2.0 * x48 + x54)
    x56 = 2.0 * x37
    x57 = x0 * (x21 + 3.0 * x29 + x31 + x56)
    x58 = x3 * (x27 + x32)
    x59 = 2.0 * x29
    x60 = 4.0 * x14
    x61 = x24 * x60
    x62 = x59 + x61
    x63 = x0 * (x42 + x51 + x56 + x62)
    x64 = x45 + x48
    x65 = x3 * x64
    x66 = x63 + x65
    x67 = x3 * x66
    x68 = x44 * x5
    x69 = x0 * (x33 + x46 + 3.0 * x48 + x68)
    x70 = x0 * (2.0 * x31 + 2.0 * x42 + x62)
    x71 = x48 + x68
    x72 = x3 * x71
    x73 = x70 + x72
    x74 = x3 * x73
    x75 = x69 + x74
    x76 = 2.0 * x0 * (x57 + x58 + x63 + x65 + x70 + x72) + x3 * x75
    x77 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x78 = da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**3.5)
    x79 = x77 * x78
    x80 = 10.1992841329868 * x79
    x81 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x82 = 0.564189583547756 * x1
    x83 = x81 * x82
    x84 = -x1 * (ax * A[1] + bx * B[1])
    x85 = -x84 - B[1]
    x86 = 17.6656783191643 * x85
    x87 = x10 * x3
    x88 = x0 * (x28 + x87)
    x89 = x3 * (x14 + x35)
    x90 = x88 + x89
    x91 = x0 * (x11 + x87)
    x92 = x14 + x34
    x93 = x3 * x92
    x94 = x91 + x93
    x95 = x0 * (2.0 * x51 + x61 + 2.0 * x88 + 2.0 * x89) + x3 * x54
    x96 = x55 + x67
    x97 = x79 * x83
    x98 = x97 * (
        x0
        * (x20 * (x56 + x59 + x90 + x94) + x4 * (x36 + x38) + 3.0 * x63 + 3.0 * x65 + x95)
        + x3 * x96
    )
    x99 = -x1 * (ax * A[2] + bx * B[2])
    x100 = -x99 - B[2]
    x101 = 17.6656783191643 * x100
    x102 = x10 * x3**2
    x103 = x102 + x15
    x104 = x0 * (x20 * (x103 + x49) + x4 * x90 + 3.0 * x50 + 3.0 * x53) + x3 * x95
    x105 = x104 * x79
    x106 = x81 * x9
    x107 = x106 * x85**2
    x108 = x0 * x106
    x109 = x107 + x108
    x110 = 10.1992841329868 * x109
    x111 = 0.318309886183791 * x8
    x112 = x110 * x111
    x113 = x83 * x86
    x114 = x77 * x9
    x115 = x100**2 * x114
    x116 = x0 * x114
    x117 = x115 + x116
    x118 = 10.1992841329868 * x117
    x119 = x111 * x78
    x120 = x119 * x81
    x121 = -x84 - A[1]
    x122 = 26.9847693667702 * x121
    x123 = x76 * x97
    x124 = x106 * x121
    x125 = x124 * x85
    x126 = x108 + x125
    x127 = 46.7389915737742 * x126
    x128 = x111 * x79
    x129 = x127 * x128
    x130 = 46.7389915737742 * x100
    x131 = x96 * x97
    x132 = x109 * x121
    x133 = x106 * x85
    x134 = x132 + x133 * x20
    x135 = 26.9847693667702 * x128
    x136 = x135 * x95
    x137 = x120 * x95
    x138 = -x99 - A[2]
    x139 = 26.9847693667702 * x138
    x140 = x114 * x138
    x141 = x100 * x140
    x142 = x116 + x141
    x143 = 46.7389915737742 * x142
    x144 = x120 * x143
    x145 = x117 * x138
    x146 = x100 * x114
    x147 = x145 + x146 * x20
    x148 = 26.9847693667702 * x147
    x149 = x106 * x121**2
    x150 = x108 + x149
    x151 = 34.8371874529163 * x150
    x152 = x0 * (x124 + x133)
    x153 = x121 * x126
    x154 = x152 + x153
    x155 = 60.3397786612521 * x154
    x156 = x128 * x66
    x157 = 60.3397786612521 * x100
    x158 = 3.0 * x108
    x159 = 2.0 * x121
    x160 = x133 * x159 + x158
    x161 = x0 * (x107 + x160)
    x162 = x121 * x134
    x163 = x161 + x162
    x164 = 34.8371874529163 * x128
    x165 = x128 * x54
    x166 = 0.179587122125167 * x78
    x167 = x166 * x54
    x168 = 60.3397786612521 * x138
    x169 = 104.511562358749 * x126
    x170 = 104.511562358749 * x142
    x171 = x120 * x66
    x172 = x120 * x54
    x173 = 60.3397786612521 * x121
    x174 = x114 * x138**2
    x175 = x116 + x174
    x176 = 34.8371874529163 * x175
    x177 = 60.3397786612521 * x175
    x178 = x0 * (x140 + x146)
    x179 = x138 * x142
    x180 = x178 + x179
    x181 = 60.3397786612521 * x180
    x182 = 3.0 * x116
    x183 = 2.0 * x138
    x184 = x146 * x183 + x182
    x185 = x0 * (x115 + x184)
    x186 = x138 * x147
    x187 = x185 + x186
    x188 = 34.8371874529163 * x187
    x189 = x121 * x150 + x124 * x20
    x190 = 26.9847693667702 * x189
    x191 = x128 * x73
    x192 = x0 * (x149 + x160)
    x193 = x121 * x154
    x194 = x192 + x193
    x195 = 46.7389915737742 * x64
    x196 = x128 * x195
    x197 = 4.0 * x108
    x198 = x0 * (2.0 * x132 + 2.0 * x152 + 2.0 * x153 + x197 * x85) + x121 * x163
    x199 = x135 * x198
    x200 = x128 * x52
    x201 = x166 * x52
    x202 = 104.511562358749 * x64
    x203 = x128 * x138
    x204 = x166 * x64
    x205 = 60.3397786612521 * x201
    x206 = x120 * x73
    x207 = x120 * x121
    x208 = x120 * x52
    x209 = x138 * x175 + x140 * x20
    x210 = 26.9847693667702 * x209
    x211 = x120 * x195
    x212 = x0 * (x174 + x184)
    x213 = x138 * x180
    x214 = x212 + x213
    x215 = 46.7389915737742 * x214
    x216 = 4.0 * x116
    x217 = x0 * (x100 * x216 + 2.0 * x145 + 2.0 * x178 + 2.0 * x179) + x138 * x187
    x218 = 26.9847693667702 * x120
    x219 = x217 * x218
    x220 = x0 * (3.0 * x149 + x158) + x121 * x189
    x221 = 10.1992841329868 * x220
    x222 = x128 * x71
    x223 = x0 * (3.0 * x152 + 3.0 * x153 + x189) + x121 * x194
    x224 = 17.6656783191643 * x44
    x225 = x128 * x224
    x226 = x0 * (3.0 * x161 + 3.0 * x162 + 2.0 * x192 + 2.0 * x193) + x121 * x198
    x227 = x111 * x80
    x228 = x166 * x41
    x229 = x166 * x44
    x230 = x166 * x175
    x231 = 60.3397786612521 * x229
    x232 = x120 * x71
    x233 = x0 * (3.0 * x174 + x182) + x138 * x209
    x234 = 10.1992841329868 * x233
    x235 = x120 * x224
    x236 = x0 * (3.0 * x178 + 3.0 * x179 + x209) + x138 * x214
    x237 = x120 * x41
    x238 = (
        10.1992841329868 * x0 * (3.0 * x185 + 3.0 * x186 + 2.0 * x212 + 2.0 * x213)
        + 10.1992841329868 * x138 * x217
    )
    x239 = -x84 - C[1]
    x240 = x106 * x239**2
    x241 = x108 + x240
    x242 = x0 * (x103 + x12)
    x243 = x3 * x94
    x244 = x0 * (2.0 * x19 + x5 * x60 + 2.0 * x91 + 2.0 * x93) + x23 * x3
    x245 = x0 * (3.0 * x17 + 3.0 * x22 + 2.0 * x242 + 2.0 * x243) + x244 * x3
    x246 = x241 * x85
    x247 = x106 * x239
    x248 = x20 * x247
    x249 = x246 + x248
    x250 = 17.6656783191643 * x249
    x251 = x102 + x14
    x252 = x20 * x87 + x251 * x3
    x253 = x242 + x243
    x254 = x0 * (x252 + 3.0 * x91 + 3.0 * x93) + x253 * x3
    x255 = x128 * x254
    x256 = x133 * x239
    x257 = 2.0 * x256
    x258 = x158 + x240
    x259 = x0 * (x257 + x258)
    x260 = x249 * x85
    x261 = x259 + x260
    x262 = x0 * (3.0 * x102 + x15) + x252 * x3
    x263 = x100 * x128
    x264 = x166 * x241
    x265 = x121 * x241
    x266 = x248 + x265
    x267 = 26.9847693667702 * x266
    x268 = x121 * x249
    x269 = x259 + x268
    x270 = 46.7389915737742 * x269
    x271 = x128 * x253
    x272 = x108 + x256
    x273 = x272 * x85
    x274 = x0 * (x133 + x247)
    x275 = 2.0 * x274
    x276 = x197 * x239
    x277 = x275 + x276
    x278 = x0 * (2.0 * x246 + 2.0 * x273 + x277)
    x279 = x121 * x261
    x280 = x278 + x279
    x281 = x135 * x252
    x282 = x166 * x252
    x283 = 46.7389915737742 * x249
    x284 = x159 * x247
    x285 = x0 * (x258 + x284)
    x286 = x121 * x266
    x287 = x285 + x286
    x288 = 34.8371874529163 * x287
    x289 = x128 * x23
    x290 = x121 * x272
    x291 = 2.0 * x290
    x292 = x0 * (x246 + x265 + x277 + x291)
    x293 = x121 * x269
    x294 = x292 + x293
    x295 = 60.3397786612521 * x94
    x296 = x128 * x295
    x297 = 2.0 * x268
    x298 = x0 * (x107 + x158 + x257)
    x299 = x121 * (x273 + x274)
    x300 = 2.0 * x298 + 2.0 * x299
    x301 = x0 * (3.0 * x259 + x260 + x297 + x300)
    x302 = x121 * x280
    x303 = x301 + x302
    x304 = x128 * x251
    x305 = x166 * x251
    x306 = 60.3397786612521 * x266
    x307 = 104.511562358749 * x269
    x308 = x166 * x170
    x309 = x0 * (x124 + x247)
    x310 = x124 * x239
    x311 = x121 * (x108 + x310)
    x312 = x0 * (2.0 * x265 + x276 + 2.0 * x309 + 2.0 * x311) + x121 * x287
    x313 = x135 * x312
    x314 = x0 * (x125 + x158 + x256 + x310)
    x315 = x121 * (x274 + x290)
    x316 = 2.0 * x314 + 2.0 * x315
    x317 = x0 * (2.0 * x259 + x287 + x297 + x316)
    x318 = x121 * x294
    x319 = x317 + x318
    x320 = 46.7389915737742 * x92
    x321 = x128 * x320
    x322 = 26.9847693667702 * x3
    x323 = x0 * (x134 + x273 + 3.0 * x274 + x291)
    x324 = x121 * (x298 + x299)
    x325 = 2.0 * x0 * (x278 + x279 + x292 + x293 + x323 + x324) + x121 * x303
    x326 = x7 * x82
    x327 = x326 * x79
    x328 = x325 * x327
    x329 = x3 * x327
    x330 = x119 * x7
    x331 = x312 * x330
    x332 = 104.511562358749 * x92
    x333 = x294 * x330
    x334 = 60.3397786612521 * x3
    x335 = x287 * x330
    x336 = x166 * x266
    x337 = x280 * x330
    x338 = x3 * x330
    x339 = x187 * x330
    x340 = x166 * x320
    x341 = x261 * x330
    x342 = x214 * x330
    x343 = x217 * x330
    x344 = x309 + x311
    x345 = (
        x0 * (x159 * x344 + x20 * (x149 + x158 + x284) + 3.0 * x285 + 3.0 * x286)
        + x121 * x312
    )
    x346 = x327 * (
        x0
        * (
            x159 * (x314 + x315)
            + x20 * (x154 + x275 + x291 + x344)
            + 3.0 * x292
            + 3.0 * x293
            + x312
        )
        + x121 * x319
    )
    x347 = 17.6656783191643 * x5
    x348 = 46.7389915737742 * x5
    x349 = x143 * x330
    x350 = x166 * x18
    x351 = x330 * x5
    x352 = x250 * x330
    x353 = x241 * x330
    x354 = -x99 - C[2]
    x355 = x114 * x354**2
    x356 = x116 + x355
    x357 = 10.1992841329868 * x120
    x358 = x120 * x254
    x359 = x100 * x356
    x360 = x114 * x354
    x361 = x20 * x360
    x362 = x359 + x361
    x363 = 17.6656783191643 * x362
    x364 = x166 * x356
    x365 = x120 * x85
    x366 = x146 * x354
    x367 = 2.0 * x366
    x368 = x182 + x355
    x369 = x0 * (x367 + x368)
    x370 = x100 * x362
    x371 = x369 + x370
    x372 = x218 * x244
    x373 = x120 * x253
    x374 = 46.7389915737742 * x373
    x375 = 26.9847693667702 * x282
    x376 = x218 * x252
    x377 = x138 * x356
    x378 = x361 + x377
    x379 = x138 * x362
    x380 = x369 + x379
    x381 = 46.7389915737742 * x380
    x382 = x116 + x366
    x383 = x100 * x382
    x384 = x0 * (x146 + x360)
    x385 = 2.0 * x384
    x386 = x216 * x354
    x387 = x385 + x386
    x388 = x0 * (2.0 * x359 + 2.0 * x383 + x387)
    x389 = x138 * x371
    x390 = x388 + x389
    x391 = x150 * x166
    x392 = 34.8371874529163 * x305
    x393 = 60.3397786612521 * x378
    x394 = x166 * x378
    x395 = 104.511562358749 * x380
    x396 = x120 * x173
    x397 = x183 * x360
    x398 = x0 * (x368 + x397)
    x399 = x138 * x378
    x400 = x398 + x399
    x401 = 34.8371874529163 * x120
    x402 = x120 * x295
    x403 = x138 * x382
    x404 = 2.0 * x403
    x405 = x0 * (x359 + x377 + x387 + x404)
    x406 = x138 * x380
    x407 = x405 + x406
    x408 = 2.0 * x379
    x409 = x0 * (x115 + x182 + x367)
    x410 = x138 * (x383 + x384)
    x411 = 2.0 * x409 + 2.0 * x410
    x412 = x0 * (3.0 * x369 + x370 + x408 + x411)
    x413 = x138 * x390
    x414 = x412 + x413
    x415 = x330 * x356
    x416 = x194 * x330
    x417 = 46.7389915737742 * x3
    x418 = x330 * x371
    x419 = x163 * x330
    x420 = x330 * x390
    x421 = x330 * x407
    x422 = x7 * x78
    x423 = x422 * x83
    x424 = x0 * (x140 + x360)
    x425 = x140 * x354
    x426 = x138 * (x116 + x425)
    x427 = x0 * (2.0 * x377 + x386 + 2.0 * x424 + 2.0 * x426) + x138 * x400
    x428 = x218 * x427
    x429 = x120 * x320
    x430 = x0 * (x141 + x182 + x366 + x425)
    x431 = x138 * (x384 + x403)
    x432 = 2.0 * x430 + 2.0 * x431
    x433 = x0 * (2.0 * x369 + x400 + x408 + x432)
    x434 = x138 * x407
    x435 = x433 + x434
    x436 = 26.9847693667702 * x330
    x437 = x427 * x436
    x438 = x423 * x435
    x439 = x0 * (x147 + x383 + 3.0 * x384 + x404)
    x440 = x138 * (x409 + x410)
    x441 = 2.0 * x0 * (x388 + x389 + x405 + x406 + x439 + x440) + x138 * x414
    x442 = x423 * x441
    x443 = x330 * x363
    x444 = x127 * x330
    x445 = x424 + x426
    x446 = (
        x0 * (x183 * x445 + x20 * (x174 + x182 + x397) + 3.0 * x398 + 3.0 * x399)
        + x138 * x427
    )
    x447 = x422 * x446
    x448 = x423 * (
        x0
        * (
            x183 * (x430 + x431)
            + x20 * (x180 + x385 + x404 + x445)
            + 3.0 * x405
            + 3.0 * x406
            + x427
        )
        + x138 * x435
    )

    # 270 item(s)
    result[0, 0, 0] = numpy.sum(
        x80
        * x83
        * (
            x0
            * (
                x20 * (x23 + x33 + x39)
                + x4 * (x57 + x58)
                + 2.0 * x55
                + 2.0 * x67
                + 3.0 * x69
                + 3.0 * x74
            )
            + x3 * x76
        )
    )
    result[0, 0, 1] = numpy.sum(x86 * x98)
    result[0, 0, 2] = numpy.sum(x101 * x98)
    result[0, 0, 3] = numpy.sum(x105 * x112)
    result[0, 0, 4] = numpy.sum(x100 * x105 * x113)
    result[0, 0, 5] = numpy.sum(x104 * x118 * x120)
    result[0, 1, 0] = numpy.sum(x122 * x123)
    result[0, 1, 1] = numpy.sum(x129 * x96)
    result[0, 1, 2] = numpy.sum(x121 * x130 * x131)
    result[0, 1, 3] = numpy.sum(x134 * x136)
    result[0, 1, 4] = numpy.sum(x100 * x129 * x95)
    result[0, 1, 5] = numpy.sum(x117 * x122 * x137)
    result[0, 2, 0] = numpy.sum(x123 * x139)
    result[0, 2, 1] = numpy.sum(46.7389915737742 * x131 * x138 * x85)
    result[0, 2, 2] = numpy.sum(x144 * x96)
    result[0, 2, 3] = numpy.sum(x109 * x136 * x138)
    result[0, 2, 4] = numpy.sum(x144 * x85 * x95)
    result[0, 2, 5] = numpy.sum(x137 * x148)
    result[0, 3, 0] = numpy.sum(x128 * x151 * x75)
    result[0, 3, 1] = numpy.sum(x155 * x156)
    result[0, 3, 2] = numpy.sum(x150 * x156 * x157)
    result[0, 3, 3] = numpy.sum(x163 * x164 * x54)
    result[0, 3, 4] = numpy.sum(x100 * x155 * x165)
    result[0, 3, 5] = numpy.sum(x117 * x151 * x167)
    result[0, 4, 0] = numpy.sum(x121 * x168 * x75 * x97)
    result[0, 4, 1] = numpy.sum(x138 * x156 * x169)
    result[0, 4, 2] = numpy.sum(x121 * x170 * x171)
    result[0, 4, 3] = numpy.sum(x134 * x165 * x168)
    result[0, 4, 4] = numpy.sum(x126 * x167 * x170)
    result[0, 4, 5] = numpy.sum(x147 * x172 * x173)
    result[0, 5, 0] = numpy.sum(x120 * x176 * x75)
    result[0, 5, 1] = numpy.sum(x171 * x177 * x85)
    result[0, 5, 2] = numpy.sum(x171 * x181)
    result[0, 5, 3] = numpy.sum(x109 * x167 * x176)
    result[0, 5, 4] = numpy.sum(x172 * x181 * x85)
    result[0, 5, 5] = numpy.sum(x172 * x188)
    result[0, 6, 0] = numpy.sum(x190 * x191)
    result[0, 6, 1] = numpy.sum(x194 * x196)
    result[0, 6, 2] = numpy.sum(x100 * x189 * x196)
    result[0, 6, 3] = numpy.sum(x199 * x52)
    result[0, 6, 4] = numpy.sum(x130 * x194 * x200)
    result[0, 6, 5] = numpy.sum(x117 * x190 * x201)
    result[0, 7, 0] = numpy.sum(x150 * x168 * x191)
    result[0, 7, 1] = numpy.sum(x154 * x202 * x203)
    result[0, 7, 2] = numpy.sum(x150 * x170 * x204)
    result[0, 7, 3] = numpy.sum(x163 * x168 * x200)
    result[0, 7, 4] = numpy.sum(x154 * x170 * x201)
    result[0, 7, 5] = numpy.sum(x147 * x150 * x205)
    result[0, 8, 0] = numpy.sum(x173 * x175 * x206)
    result[0, 8, 1] = numpy.sum(x169 * x175 * x204)
    result[0, 8, 2] = numpy.sum(x180 * x202 * x207)
    result[0, 8, 3] = numpy.sum(x134 * x175 * x205)
    result[0, 8, 4] = numpy.sum(x169 * x180 * x201)
    result[0, 8, 5] = numpy.sum(x173 * x187 * x208)
    result[0, 9, 0] = numpy.sum(x206 * x210)
    result[0, 9, 1] = numpy.sum(x209 * x211 * x85)
    result[0, 9, 2] = numpy.sum(x211 * x214)
    result[0, 9, 3] = numpy.sum(x109 * x201 * x210)
    result[0, 9, 4] = numpy.sum(x208 * x215 * x85)
    result[0, 9, 5] = numpy.sum(x219 * x52)
    result[0, 10, 0] = numpy.sum(x221 * x222)
    result[0, 10, 1] = numpy.sum(x223 * x225)
    result[0, 10, 2] = numpy.sum(x100 * x220 * x225)
    result[0, 10, 3] = numpy.sum(x226 * x227 * x41)
    result[0, 10, 4] = numpy.sum(x101 * x128 * x223 * x41)
    result[0, 10, 5] = numpy.sum(x118 * x220 * x228)
    result[0, 11, 0] = numpy.sum(x138 * x190 * x222)
    result[0, 11, 1] = numpy.sum(46.7389915737742 * x194 * x203 * x44)
    result[0, 11, 2] = numpy.sum(x143 * x189 * x229)
    result[0, 11, 3] = numpy.sum(x138 * x199 * x41)
    result[0, 11, 4] = numpy.sum(x143 * x194 * x228)
    result[0, 11, 5] = numpy.sum(x147 * x190 * x228)
    result[0, 12, 0] = numpy.sum(x151 * x230 * x71)
    result[0, 12, 1] = numpy.sum(x154 * x175 * x231)
    result[0, 12, 2] = numpy.sum(x150 * x180 * x231)
    result[0, 12, 3] = numpy.sum(x163 * x176 * x228)
    result[0, 12, 4] = numpy.sum(x154 * x181 * x228)
    result[0, 12, 5] = numpy.sum(x151 * x187 * x228)
    result[0, 13, 0] = numpy.sum(x121 * x210 * x232)
    result[0, 13, 1] = numpy.sum(x127 * x209 * x229)
    result[0, 13, 2] = numpy.sum(x207 * x215 * x44)
    result[0, 13, 3] = numpy.sum(x134 * x210 * x228)
    result[0, 13, 4] = numpy.sum(x127 * x214 * x228)
    result[0, 13, 5] = numpy.sum(x121 * x219 * x41)
    result[0, 14, 0] = numpy.sum(x232 * x234)
    result[0, 14, 1] = numpy.sum(x233 * x235 * x85)
    result[0, 14, 2] = numpy.sum(x235 * x236)
    result[0, 14, 3] = numpy.sum(x109 * x228 * x234)
    result[0, 14, 4] = numpy.sum(x236 * x237 * x86)
    result[0, 14, 5] = numpy.sum(x237 * x238)
    result[1, 0, 0] = numpy.sum(x227 * x241 * x245)
    result[1, 0, 1] = numpy.sum(x250 * x255)
    result[1, 0, 2] = numpy.sum(x101 * x241 * x255)
    result[1, 0, 3] = numpy.sum(x227 * x261 * x262)
    result[1, 0, 4] = numpy.sum(x250 * x262 * x263)
    result[1, 0, 5] = numpy.sum(x118 * x262 * x264)
    result[1, 1, 0] = numpy.sum(x128 * x244 * x267)
    result[1, 1, 1] = numpy.sum(x270 * x271)
    result[1, 1, 2] = numpy.sum(x130 * x266 * x271)
    result[1, 1, 3] = numpy.sum(x280 * x281)
    result[1, 1, 4] = numpy.sum(x252 * x263 * x270)
    result[1, 1, 5] = numpy.sum(x117 * x267 * x282)
    result[1, 2, 0] = numpy.sum(x135 * x138 * x241 * x244)
    result[1, 2, 1] = numpy.sum(x138 * x271 * x283)
    result[1, 2, 2] = numpy.sum(x143 * x253 * x264)
    result[1, 2, 3] = numpy.sum(x138 * x261 * x281)
    result[1, 2, 4] = numpy.sum(x143 * x249 * x282)
    result[1, 2, 5] = numpy.sum(x148 * x252 * x264)
    result[1, 3, 0] = numpy.sum(x288 * x289)
    result[1, 3, 1] = numpy.sum(x294 * x296)
    result[1, 3, 2] = numpy.sum(x100 * x287 * x296)
    result[1, 3, 3] = numpy.sum(x164 * x251 * x303)
    result[1, 3, 4] = numpy.sum(x157 * x294 * x304)
    result[1, 3, 5] = numpy.sum(x117 * x288 * x305)
    result[1, 4, 0] = numpy.sum(x138 * x289 * x306)
    result[1, 4, 1] = numpy.sum(x203 * x307 * x94)
    result[1, 4, 2] = numpy.sum(x266 * x308 * x94)
    result[1, 4, 3] = numpy.sum(x168 * x280 * x304)
    result[1, 4, 4] = numpy.sum(x170 * x269 * x305)
    result[1, 4, 5] = numpy.sum(x147 * x305 * x306)
    result[1, 5, 0] = numpy.sum(x176 * x23 * x264)
    result[1, 5, 1] = numpy.sum(x230 * x249 * x295)
    result[1, 5, 2] = numpy.sum(x181 * x264 * x94)
    result[1, 5, 3] = numpy.sum(x176 * x261 * x305)
    result[1, 5, 4] = numpy.sum(x181 * x249 * x305)
    result[1, 5, 5] = numpy.sum(x188 * x251 * x264)
    result[1, 6, 0] = numpy.sum(x21 * x313)
    result[1, 6, 1] = numpy.sum(x319 * x321)
    result[1, 6, 2] = numpy.sum(x100 * x312 * x321)
    result[1, 6, 3] = numpy.sum(x322 * x328)
    result[1, 6, 4] = numpy.sum(x130 * x319 * x329)
    result[1, 6, 5] = numpy.sum(x117 * x322 * x331)
    result[1, 7, 0] = numpy.sum(x128 * x168 * x21 * x287)
    result[1, 7, 1] = numpy.sum(x203 * x294 * x332)
    result[1, 7, 2] = numpy.sum(x287 * x308 * x92)
    result[1, 7, 3] = numpy.sum(x168 * x303 * x329)
    result[1, 7, 4] = numpy.sum(x170 * x3 * x333)
    result[1, 7, 5] = numpy.sum(x147 * x334 * x335)
    result[1, 8, 0] = numpy.sum(x21 * x230 * x306)
    result[1, 8, 1] = numpy.sum(x230 * x269 * x332)
    result[1, 8, 2] = numpy.sum(x180 * x332 * x336)
    result[1, 8, 3] = numpy.sum(x175 * x334 * x337)
    result[1, 8, 4] = numpy.sum(x180 * x307 * x338)
    result[1, 8, 5] = numpy.sum(x3 * x306 * x339)
    result[1, 9, 0] = numpy.sum(x21 * x210 * x264)
    result[1, 9, 1] = numpy.sum(x209 * x249 * x340)
    result[1, 9, 2] = numpy.sum(x214 * x264 * x320)
    result[1, 9, 3] = numpy.sum(x210 * x3 * x341)
    result[1, 9, 4] = numpy.sum(x283 * x3 * x342)
    result[1, 9, 5] = numpy.sum(x241 * x322 * x343)
    result[1, 10, 0] = numpy.sum(x18 * x227 * x345)
    result[1, 10, 1] = numpy.sum(x346 * x347)
    result[1, 10, 2] = numpy.sum(x101 * x327 * x345 * x5)
    result[1, 10, 3] = numpy.sum(
        x326
        * x80
        * (
            x0
            * (
                x159 * (x323 + x324)
                + x20 * (x163 + x300 + x316)
                + 3.0 * x301
                + 3.0 * x302
                + 2.0 * x317
                + 2.0 * x318
            )
            + x121 * x325
        )
    )
    result[1, 10, 4] = numpy.sum(x101 * x346)
    result[1, 10, 5] = numpy.sum(x118 * x330 * x345)
    result[1, 11, 0] = numpy.sum(x138 * x18 * x313)
    result[1, 11, 1] = numpy.sum(x138 * x319 * x327 * x348)
    result[1, 11, 2] = numpy.sum(x312 * x349 * x5)
    result[1, 11, 3] = numpy.sum(x139 * x328)
    result[1, 11, 4] = numpy.sum(x319 * x349)
    result[1, 11, 5] = numpy.sum(x148 * x331)
    result[1, 12, 0] = numpy.sum(x176 * x287 * x350)
    result[1, 12, 1] = numpy.sum(x177 * x333 * x5)
    result[1, 12, 2] = numpy.sum(x181 * x335 * x5)
    result[1, 12, 3] = numpy.sum(x176 * x303 * x330)
    result[1, 12, 4] = numpy.sum(x181 * x333)
    result[1, 12, 5] = numpy.sum(x288 * x339)
    result[1, 13, 0] = numpy.sum(x18 * x210 * x336)
    result[1, 13, 1] = numpy.sum(x209 * x270 * x351)
    result[1, 13, 2] = numpy.sum(x266 * x342 * x348)
    result[1, 13, 3] = numpy.sum(x210 * x337)
    result[1, 13, 4] = numpy.sum(x270 * x342)
    result[1, 13, 5] = numpy.sum(x267 * x343)
    result[1, 14, 0] = numpy.sum(x18 * x234 * x264)
    result[1, 14, 1] = numpy.sum(x233 * x352 * x5)
    result[1, 14, 2] = numpy.sum(x236 * x347 * x353)
    result[1, 14, 3] = numpy.sum(x234 * x341)
    result[1, 14, 4] = numpy.sum(x236 * x352)
    result[1, 14, 5] = numpy.sum(x238 * x353)
    result[2, 0, 0] = numpy.sum(x245 * x356 * x357)
    result[2, 0, 1] = numpy.sum(x356 * x358 * x86)
    result[2, 0, 2] = numpy.sum(x358 * x363)
    result[2, 0, 3] = numpy.sum(x110 * x262 * x364)
    result[2, 0, 4] = numpy.sum(x262 * x363 * x365)
    result[2, 0, 5] = numpy.sum(x262 * x357 * x371)
    result[2, 1, 0] = numpy.sum(x121 * x356 * x372)
    result[2, 1, 1] = numpy.sum(x127 * x253 * x364)
    result[2, 1, 2] = numpy.sum(x121 * x362 * x374)
    result[2, 1, 3] = numpy.sum(x134 * x356 * x375)
    result[2, 1, 4] = numpy.sum(x127 * x282 * x362)
    result[2, 1, 5] = numpy.sum(x121 * x371 * x376)
    result[2, 2, 0] = numpy.sum(x372 * x378)
    result[2, 2, 1] = numpy.sum(x374 * x378 * x85)
    result[2, 2, 2] = numpy.sum(x373 * x381)
    result[2, 2, 3] = numpy.sum(x109 * x375 * x378)
    result[2, 2, 4] = numpy.sum(x252 * x365 * x381)
    result[2, 2, 5] = numpy.sum(x376 * x390)
    result[2, 3, 0] = numpy.sum(x151 * x23 * x364)
    result[2, 3, 1] = numpy.sum(x154 * x295 * x364)
    result[2, 3, 2] = numpy.sum(x295 * x362 * x391)
    result[2, 3, 3] = numpy.sum(x163 * x356 * x392)
    result[2, 3, 4] = numpy.sum(x155 * x305 * x362)
    result[2, 3, 5] = numpy.sum(x151 * x305 * x371)
    result[2, 4, 0] = numpy.sum(x207 * x23 * x393)
    result[2, 4, 1] = numpy.sum(x169 * x394 * x94)
    result[2, 4, 2] = numpy.sum(x207 * x395 * x94)
    result[2, 4, 3] = numpy.sum(x134 * x305 * x393)
    result[2, 4, 4] = numpy.sum(x169 * x305 * x380)
    result[2, 4, 5] = numpy.sum(x251 * x390 * x396)
    result[2, 5, 0] = numpy.sum(x23 * x400 * x401)
    result[2, 5, 1] = numpy.sum(x400 * x402 * x85)
    result[2, 5, 2] = numpy.sum(x402 * x407)
    result[2, 5, 3] = numpy.sum(x109 * x392 * x400)
    result[2, 5, 4] = numpy.sum(60.3397786612521 * x251 * x365 * x407)
    result[2, 5, 5] = numpy.sum(x251 * x401 * x414)
    result[2, 6, 0] = numpy.sum(x190 * x21 * x364)
    result[2, 6, 1] = numpy.sum(x194 * x340 * x356)
    result[2, 6, 2] = numpy.sum(x189 * x340 * x362)
    result[2, 6, 3] = numpy.sum(x198 * x322 * x415)
    result[2, 6, 4] = numpy.sum(x362 * x416 * x417)
    result[2, 6, 5] = numpy.sum(x190 * x3 * x418)
    result[2, 7, 0] = numpy.sum(x21 * x391 * x393)
    result[2, 7, 1] = numpy.sum(x154 * x332 * x394)
    result[2, 7, 2] = numpy.sum(x332 * x380 * x391)
    result[2, 7, 3] = numpy.sum(x3 * x393 * x419)
    result[2, 7, 4] = numpy.sum(x154 * x338 * x395)
    result[2, 7, 5] = numpy.sum(x150 * x334 * x420)
    result[2, 8, 0] = numpy.sum(x21 * x396 * x400)
    result[2, 8, 1] = numpy.sum(x166 * x169 * x400 * x92)
    result[2, 8, 2] = numpy.sum(x207 * x332 * x407)
    result[2, 8, 3] = numpy.sum(x134 * x330 * x334 * x400)
    result[2, 8, 4] = numpy.sum(x169 * x3 * x421)
    result[2, 8, 5] = numpy.sum(x173 * x3 * x414 * x423)
    result[2, 9, 0] = numpy.sum(x21 * x428)
    result[2, 9, 1] = numpy.sum(x427 * x429 * x85)
    result[2, 9, 2] = numpy.sum(x429 * x435)
    result[2, 9, 3] = numpy.sum(x109 * x3 * x437)
    result[2, 9, 4] = numpy.sum(x417 * x438 * x85)
    result[2, 9, 5] = numpy.sum(x322 * x442)
    result[2, 10, 0] = numpy.sum(x221 * x350 * x356)
    result[2, 10, 1] = numpy.sum(x223 * x347 * x415)
    result[2, 10, 2] = numpy.sum(x220 * x443 * x5)
    result[2, 10, 3] = numpy.sum(10.1992841329868 * x226 * x415)
    result[2, 10, 4] = numpy.sum(x223 * x443)
    result[2, 10, 5] = numpy.sum(x221 * x418)
    result[2, 11, 0] = numpy.sum(x190 * x350 * x378)
    result[2, 11, 1] = numpy.sum(x348 * x378 * x416)
    result[2, 11, 2] = numpy.sum(x189 * x351 * x381)
    result[2, 11, 3] = numpy.sum(x198 * x378 * x436)
    result[2, 11, 4] = numpy.sum(x381 * x416)
    result[2, 11, 5] = numpy.sum(x190 * x420)
    result[2, 12, 0] = numpy.sum(x151 * x350 * x400)
    result[2, 12, 1] = numpy.sum(x155 * x351 * x400)
    result[2, 12, 2] = numpy.sum(60.3397786612521 * x150 * x421 * x5)
    result[2, 12, 3] = numpy.sum(34.8371874529163 * x400 * x419)
    result[2, 12, 4] = numpy.sum(x155 * x421)
    result[2, 12, 5] = numpy.sum(x151 * x330 * x414)
    result[2, 13, 0] = numpy.sum(x121 * x18 * x428)
    result[2, 13, 1] = numpy.sum(x427 * x444 * x5)
    result[2, 13, 2] = numpy.sum(x121 * x348 * x438)
    result[2, 13, 3] = numpy.sum(x134 * x437)
    result[2, 13, 4] = numpy.sum(x435 * x444)
    result[2, 13, 5] = numpy.sum(x122 * x442)
    result[2, 14, 0] = numpy.sum(x18 * x357 * x446)
    result[2, 14, 1] = numpy.sum(x113 * x447 * x5)
    result[2, 14, 2] = numpy.sum(x347 * x448)
    result[2, 14, 3] = numpy.sum(x112 * x447)
    result[2, 14, 4] = numpy.sum(x448 * x86)
    result[2, 14, 5] = numpy.sum(
        10.1992841329868
        * x423
        * (
            x0
            * (
                x183 * (x439 + x440)
                + x20 * (x187 + x411 + x432)
                + 3.0 * x412
                + 3.0 * x413
                + 2.0 * x433
                + 2.0 * x434
            )
            + x138 * x441
        )
    )
    return result


def diag_quadrupole3d_43(ax, da, A, bx, db, B, C):
    """Cartesian 3D (gf) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 15, 10), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x0 * x8
    x10 = -x2 - C[0]
    x11 = x3 * x8
    x12 = x10 * x11
    x13 = x12 + x9
    x14 = x13 * x3
    x15 = x10 * x8
    x16 = x0 * (x11 + x15)
    x17 = 3.0 * x16
    x18 = -x2 - A[0]
    x19 = x13 * x18
    x20 = 2.0 * x19
    x21 = x3**2 * x8
    x22 = x21 + x9
    x23 = x18 * x22
    x24 = 2.0 * x0
    x25 = x11 * x24
    x26 = x23 + x25
    x27 = x0 * (x14 + x17 + x20 + x26)
    x28 = 3.0 * x9
    x29 = 2.0 * x12 + x28
    x30 = x0 * (x21 + x29)
    x31 = x14 + x16
    x32 = x18 * x31
    x33 = x18 * (x30 + x32)
    x34 = x22 * x3
    x35 = x3 * x9
    x36 = x0 * (3.0 * x23 + x34 + 8.0 * x35)
    x37 = x0 * (3.0 * x21 + x28)
    x38 = x25 + x34
    x39 = x18 * x38
    x40 = x37 + x39
    x41 = x18 * x40
    x42 = x36 + x41
    x43 = x0 * (3.0 * x14 + x17 + x38)
    x44 = x3 * x31
    x45 = x18 * (x30 + x44)
    x46 = 2.0 * x43 + 2.0 * x45
    x47 = x0 * (4.0 * x30 + 3.0 * x32 + x40 + x44)
    x48 = x18 * (x43 + x45)
    x49 = 2.0 * x18
    x50 = x10**2 * x8
    x51 = x0 * (x29 + x50)
    x52 = x50 + x9
    x53 = x3 * x52
    x54 = x15 * x24
    x55 = x53 + x54
    x56 = x3 * x55
    x57 = x51 + x56
    x58 = x3 * x57
    x59 = x18 * x57
    x60 = 2.0 * x16
    x61 = 4.0 * x10 * x9
    x62 = x60 + x61
    x63 = x0 * (2.0 * x14 + 2.0 * x53 + x62)
    x64 = x0 * (x46 + x58 + 3.0 * x59 + 4.0 * x63)
    x65 = x18 * x52
    x66 = x0 * (x20 + x53 + x62 + x65)
    x67 = x18 * x55
    x68 = x51 + x67
    x69 = x18 * x68
    x70 = 2.0 * x0 * (x27 + x33 + x59 + x63 + x66 + x69)
    x71 = 2.0 * x30
    x72 = 3.0 * x51 + x71
    x73 = x0 * (2.0 * x44 + 3.0 * x56 + x72)
    x74 = x58 + x63
    x75 = x18 * x74
    x76 = x73 + x75
    x77 = x18 * x76
    x78 = 2.0 * x32
    x79 = 2.0 * x67
    x80 = x0 * (x56 + x72 + x78 + x79)
    x81 = x59 + x63
    x82 = x18 * x81
    x83 = x80 + x82
    x84 = x18 * x83
    x85 = 3.0 * x80 + 3.0 * x82
    x86 = x64 + x77
    x87 = x0 * (2.0 * x47 + 2.0 * x48 + 2.0 * x73 + 2.0 * x75 + x85) + x18 * x86
    x88 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x89 = da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**4.5)
    x90 = x88 * x89
    x91 = 9.12251705727742 * x90
    x92 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x93 = 0.564189583547756 * x1
    x94 = x92 * x93
    x95 = -x1 * (ax * A[1] + bx * B[1])
    x96 = -x95 - B[1]
    x97 = 20.3985682659737 * x96
    x98 = x11 * x49 + x28
    x99 = x0 * (x21 + x98)
    x100 = x18 * x26
    x101 = x100 + x99
    x102 = x11 * x18
    x103 = x15 * x18
    x104 = x0 * (x102 + x103 + x12 + x28)
    x105 = x18 * (x16 + x19)
    x106 = 2.0 * x104 + 2.0 * x105
    x107 = x15 * x49 + x28
    x108 = x0 * (x107 + x50)
    x109 = x54 + x65
    x110 = x109 * x18
    x111 = x108 + x110
    x112 = x0 * (x106 + x111 + 2.0 * x51 + x79)
    x113 = x66 + x69
    x114 = x113 * x18
    x115 = x70 + x84
    x116 = x90 * x94
    x117 = x116 * (
        x0
        * (
            2.0 * x112
            + 2.0 * x114
            + x24 * (x101 + x106 + x71 + x78)
            + x49 * (x27 + x33)
            + x85
        )
        + x115 * x18
    )
    x118 = -x1 * (ax * A[2] + bx * B[2])
    x119 = -x118 - B[2]
    x120 = 20.3985682659737 * x119
    x121 = x18 * x8
    x122 = x0 * (x121 + x15)
    x123 = x18 * (x103 + x9)
    x124 = x122 + x123
    x125 = x0 * (x11 + x121)
    x126 = x102 + x9
    x127 = x126 * x18
    x128 = x125 + x127
    x129 = x0 * (2.0 * x122 + 2.0 * x123 + x61 + 2.0 * x65) + x111 * x18
    x130 = x112 + x114
    x131 = (
        x0
        * (
            x129
            + x24 * (x124 + x128 + x20 + x60)
            + x49 * (x104 + x105)
            + 3.0 * x66
            + 3.0 * x69
        )
        + x130 * x18
    )
    x132 = x131 * x90
    x133 = x7 * x92
    x134 = x133 * x96**2
    x135 = x0 * x133
    x136 = x134 + x135
    x137 = 20.3985682659737 * x136
    x138 = 0.318309886183791 * x6
    x139 = x137 * x138
    x140 = 35.3313566383285 * x96
    x141 = x7 * x88
    x142 = x119**2 * x141
    x143 = x0 * x141
    x144 = x142 + x143
    x145 = 20.3985682659737 * x138
    x146 = x89 * x92
    x147 = x145 * x146
    x148 = x144 * x147
    x149 = x136 * x96
    x150 = x133 * x96
    x151 = x150 * x24
    x152 = x149 + x151
    x153 = 9.12251705727742 * x152
    x154 = x18**2 * x8
    x155 = x0 * (3.0 * x108 + 3.0 * x110 + x124 * x49 + x24 * (x107 + x154)) + x129 * x18
    x156 = x138 * x155
    x157 = x119 * x90
    x158 = x119 * x144
    x159 = x119 * x141
    x160 = x159 * x24
    x161 = x158 + x160
    x162 = 9.12251705727742 * x161
    x163 = -x95 - A[1]
    x164 = 24.1359114645008 * x163
    x165 = x116 * x87
    x166 = x133 * x163
    x167 = x166 * x96
    x168 = x135 + x167
    x169 = 53.9695387335403 * x168
    x170 = x138 * x90
    x171 = 53.9695387335403 * x115
    x172 = x116 * x171
    x173 = x136 * x163
    x174 = x151 + x173
    x175 = 53.9695387335403 * x174
    x176 = x130 * x170
    x177 = 93.4779831475484 * x168
    x178 = x138 * x146
    x179 = x130 * x178
    x180 = 53.9695387335403 * x163
    x181 = 3.0 * x135
    x182 = x0 * (3.0 * x134 + x181)
    x183 = x152 * x163
    x184 = x182 + x183
    x185 = 24.1359114645008 * x170
    x186 = x129 * x185
    x187 = 53.9695387335403 * x129
    x188 = x119 * x170
    x189 = 0.179587122125167 * x89
    x190 = x144 * x189
    x191 = x129 * x178
    x192 = -x118 - A[2]
    x193 = 24.1359114645008 * x192
    x194 = x141 * x192
    x195 = x119 * x194
    x196 = x143 + x195
    x197 = x178 * x196
    x198 = 53.9695387335403 * x136
    x199 = 93.4779831475484 * x96
    x200 = x144 * x192
    x201 = x160 + x200
    x202 = 53.9695387335403 * x201
    x203 = x189 * x196
    x204 = x178 * x96
    x205 = 3.0 * x143
    x206 = x0 * (3.0 * x142 + x205)
    x207 = x161 * x192
    x208 = x206 + x207
    x209 = 24.1359114645008 * x208
    x210 = x133 * x163**2
    x211 = x135 + x210
    x212 = 31.1593277158494 * x211
    x213 = 69.6743749058326 * x83
    x214 = x0 * (x150 + x166)
    x215 = x163 * x168
    x216 = x214 + x215
    x217 = x170 * x216
    x218 = 2.0 * x163
    x219 = x150 * x218 + x181
    x220 = x0 * (x134 + x219)
    x221 = x163 * x174
    x222 = x220 + x221
    x223 = 69.6743749058326 * x222
    x224 = x113 * x170
    x225 = 120.679557322504 * x216
    x226 = 69.6743749058326 * x190
    x227 = x135 * x96
    x228 = x0 * (x149 + 3.0 * x173 + 8.0 * x227)
    x229 = x163 * x184
    x230 = x228 + x229
    x231 = 31.1593277158494 * x230
    x232 = x111 * x170
    x233 = x111 * x189
    x234 = 120.679557322504 * x168
    x235 = x170 * x192
    x236 = 120.679557322504 * x163
    x237 = 120.679557322504 * x192
    x238 = 209.023124717498 * x203
    x239 = 120.679557322504 * x201
    x240 = x113 * x178
    x241 = 53.9695387335403 * x184
    x242 = 120.679557322504 * x203
    x243 = x111 * x178
    x244 = x141 * x192**2
    x245 = x143 + x244
    x246 = 31.1593277158494 * x245
    x247 = x178 * x213
    x248 = x0 * (x159 + x194)
    x249 = x192 * x196
    x250 = x248 + x249
    x251 = 69.6743749058326 * x136
    x252 = x189 * x245
    x253 = 2.0 * x192
    x254 = x159 * x253 + x205
    x255 = x0 * (x142 + x254)
    x256 = x192 * x201
    x257 = x255 + x256
    x258 = 69.6743749058326 * x257
    x259 = 31.1593277158494 * x152
    x260 = x119 * x143
    x261 = x0 * (x158 + 3.0 * x200 + 8.0 * x260)
    x262 = x192 * x208
    x263 = x261 + x262
    x264 = 31.1593277158494 * x263
    x265 = x163 * x211 + x166 * x24
    x266 = 24.1359114645008 * x265
    x267 = x170 * x76
    x268 = 53.9695387335403 * x81
    x269 = x0 * (x210 + x219)
    x270 = x163 * x216
    x271 = x269 + x270
    x272 = x170 * x271
    x273 = 53.9695387335403 * x265
    x274 = 2.0 * x0 * (x173 + x214 + x215 + 2.0 * x227)
    x275 = x163 * x222
    x276 = x274 + x275
    x277 = 53.9695387335403 * x276
    x278 = x170 * x68
    x279 = 93.4779831475484 * x271
    x280 = 53.9695387335403 * x190
    x281 = 3.0 * x220 + 3.0 * x221
    x282 = x0 * (2.0 * x182 + 2.0 * x183 + x281) + x163 * x230
    x283 = x185 * x282
    x284 = x109 * x189
    x285 = 53.9695387335403 * x192
    x286 = x189 * x68
    x287 = 53.9695387335403 * x170
    x288 = x192 * x287
    x289 = 53.9695387335403 * x284
    x290 = x178 * x76
    x291 = 120.679557322504 * x252
    x292 = x178 * x236
    x293 = 209.023124717498 * x250
    x294 = x178 * x68
    x295 = 53.9695387335403 * x252
    x296 = 120.679557322504 * x284
    x297 = 53.9695387335403 * x178
    x298 = x163 * x297
    x299 = x192 * x245 + x194 * x24
    x300 = 24.1359114645008 * x299
    x301 = 53.9695387335403 * x299
    x302 = x0 * (x244 + x254)
    x303 = x192 * x250
    x304 = x302 + x303
    x305 = x178 * x304
    x306 = 2.0 * x0 * (x200 + x248 + x249 + 2.0 * x260)
    x307 = x192 * x257
    x308 = x306 + x307
    x309 = 53.9695387335403 * x308
    x310 = 3.0 * x255 + 3.0 * x256
    x311 = x0 * (2.0 * x206 + 2.0 * x207 + x310) + x192 * x263
    x312 = 24.1359114645008 * x178
    x313 = x311 * x312
    x314 = x0 * (x181 + 3.0 * x210) + x163 * x265
    x315 = 9.12251705727742 * x314
    x316 = x170 * x74
    x317 = x0 * (3.0 * x214 + 3.0 * x215 + x265) + x163 * x271
    x318 = x145 * x90
    x319 = x318 * x57
    x320 = x0 * (2.0 * x269 + 2.0 * x270 + x281) + x163 * x276
    x321 = x318 * x320
    x322 = 35.3313566383285 * x55
    x323 = 20.3985682659737 * x190
    x324 = 3.0 * x0 * (x228 + x229 + x274 + x275) + x163 * x282
    x325 = x138 * x91
    x326 = x189 * x52
    x327 = 93.4779831475484 * x55
    x328 = x189 * x55
    x329 = 53.9695387335403 * x203
    x330 = 53.9695387335403 * x326
    x331 = 69.6743749058326 * x57
    x332 = x189 * x250
    x333 = 31.1593277158494 * x252
    x334 = x178 * x74
    x335 = x169 * x189
    x336 = x189 * x304
    x337 = x0 * (x205 + 3.0 * x244) + x192 * x299
    x338 = 9.12251705727742 * x337
    x339 = x147 * x57
    x340 = x0 * (3.0 * x248 + 3.0 * x249 + x299) + x192 * x304
    x341 = x0 * (2.0 * x302 + 2.0 * x303 + x310) + x192 * x308
    x342 = x147 * x341
    x343 = (
        27.3675511718323 * x0 * (x261 + x262 + x306 + x307)
        + 9.12251705727742 * x192 * x311
    )
    x344 = -x95 - C[1]
    x345 = x133 * x344**2
    x346 = x135 + x345
    x347 = 2.0 * x0 * (x125 + x127 + x23 + 2.0 * x35)
    x348 = x101 * x18
    x349 = 3.0 * x100 + 3.0 * x99
    x350 = x0 * (x349 + 2.0 * x37 + 2.0 * x39) + x18 * x42
    x351 = 3.0 * x0 * (x347 + x348 + x36 + x41) + x18 * x350
    x352 = x346 * x96
    x353 = x133 * x344
    x354 = x24 * x353
    x355 = x352 + x354
    x356 = x0 * (x154 + x98)
    x357 = x128 * x18
    x358 = x347 + x348
    x359 = x0 * (x349 + 2.0 * x356 + 2.0 * x357) + x18 * x358
    x360 = x318 * x359
    x361 = x154 + x9
    x362 = x121 * x24 + x18 * x361
    x363 = x356 + x357
    x364 = x0 * (3.0 * x125 + 3.0 * x127 + x362) + x18 * x363
    x365 = x150 * x344
    x366 = x181 + 2.0 * x365
    x367 = x0 * (x345 + x366)
    x368 = x355 * x96
    x369 = x367 + x368
    x370 = x318 * x369
    x371 = 35.3313566383285 * x355
    x372 = x135 + x365
    x373 = x372 * x96
    x374 = x0 * (x150 + x353)
    x375 = 2.0 * x374
    x376 = 4.0 * x135 * x344
    x377 = x375 + x376
    x378 = x0 * (2.0 * x352 + 2.0 * x373 + x377)
    x379 = x369 * x96
    x380 = x378 + x379
    x381 = x0 * (3.0 * x154 + x28) + x18 * x362
    x382 = x189 * x346
    x383 = x163 * x346
    x384 = x354 + x383
    x385 = 24.1359114645008 * x384
    x386 = x163 * x355
    x387 = x367 + x386
    x388 = x287 * x358
    x389 = 53.9695387335403 * x384
    x390 = x163 * x369
    x391 = x378 + x390
    x392 = x287 * x391
    x393 = 93.4779831475484 * x188
    x394 = x373 + x374
    x395 = x394 * x96
    x396 = x0 * (x134 + x366)
    x397 = 2.0 * x396
    x398 = 3.0 * x367 + x397
    x399 = x0 * (3.0 * x368 + 2.0 * x395 + x398)
    x400 = x163 * x380
    x401 = x399 + x400
    x402 = x185 * x362
    x403 = x189 * x362
    x404 = 93.4779831475484 * x355
    x405 = x181 + x218 * x353
    x406 = x0 * (x345 + x405)
    x407 = x163 * x384
    x408 = x406 + x407
    x409 = 31.1593277158494 * x408
    x410 = x170 * x42
    x411 = x163 * x372
    x412 = 2.0 * x411
    x413 = x0 * (x352 + x377 + x383 + x412)
    x414 = x163 * x387
    x415 = x413 + x414
    x416 = 69.6743749058326 * x101
    x417 = x170 * x416
    x418 = x163 * x394
    x419 = 2.0 * x418
    x420 = 2.0 * x386
    x421 = x0 * (x368 + x398 + x419 + x420)
    x422 = x163 * x391
    x423 = x421 + x422
    x424 = 69.6743749058326 * x423
    x425 = x128 * x170
    x426 = 120.679557322504 * x128
    x427 = 3.0 * x374
    x428 = x0 * (x152 + 3.0 * x373 + x427)
    x429 = x163 * (x395 + x396)
    x430 = 2.0 * x428 + 2.0 * x429
    x431 = x0 * (4.0 * x378 + x379 + 3.0 * x390 + x430)
    x432 = x163 * x401
    x433 = x431 + x432
    x434 = 31.1593277158494 * x361
    x435 = x189 * x361
    x436 = 120.679557322504 * x235
    x437 = x189 * x239
    x438 = 69.6743749058326 * x369
    x439 = x0 * (x166 + x353)
    x440 = x166 * x344
    x441 = x163 * (x135 + x440)
    x442 = x0 * (x376 + 2.0 * x383 + 2.0 * x439 + 2.0 * x441) + x163 * x408
    x443 = x185 * x442
    x444 = x0 * (x167 + x181 + x365 + x440)
    x445 = x163 * (x374 + x411)
    x446 = 2.0 * x444 + 2.0 * x445
    x447 = x0 * (2.0 * x367 + x408 + x420 + x446)
    x448 = x163 * x415
    x449 = x447 + x448
    x450 = 53.9695387335403 * x26
    x451 = x170 * x450
    x452 = x0 * (x174 + x373 + x412 + x427)
    x453 = x163 * (x396 + x418)
    x454 = 2.0 * x0 * (x378 + x390 + x413 + x414 + x452 + x453)
    x455 = x163 * x423
    x456 = x454 + x455
    x457 = 24.1359114645008 * x18
    x458 = x0 * (x184 + x395 + 4.0 * x396 + 3.0 * x418)
    x459 = x163 * (x428 + x429)
    x460 = 3.0 * x421 + 3.0 * x422
    x461 = x0 * (2.0 * x399 + 2.0 * x400 + 2.0 * x458 + 2.0 * x459 + x460) + x163 * x433
    x462 = x5 * x93
    x463 = x462 * x90
    x464 = x461 * x463
    x465 = x157 * x462
    x466 = 53.9695387335403 * x18
    x467 = x5 * x89
    x468 = x138 * x467
    x469 = 53.9695387335403 * x468
    x470 = x18 * x469
    x471 = x442 * x468
    x472 = x285 * x463
    x473 = 120.679557322504 * x18
    x474 = x468 * x473
    x475 = x415 * x468
    x476 = 120.679557322504 * x384
    x477 = x126 * x189
    x478 = x250 * x468
    x479 = x387 * x468
    x480 = x263 * x468
    x481 = x189 * x301
    x482 = x300 * x468
    x483 = x304 * x469
    x484 = x18 * x468
    x485 = x311 * x468
    x486 = x439 + x441
    x487 = (
        x0 * (x218 * x486 + x24 * (x210 + x405) + 3.0 * x406 + 3.0 * x407) + x163 * x442
    )
    x488 = (
        x0
        * (
            x218 * (x444 + x445)
            + x24 * (x216 + x375 + x412 + x486)
            + 3.0 * x413
            + 3.0 * x414
            + x442
        )
        + x163 * x449
    )
    x489 = 20.3985682659737 * x22
    x490 = x170 * x489
    x491 = x463 * (
        x0
        * (
            x218 * (x452 + x453)
            + x24 * (x222 + x397 + x419 + x446)
            + 2.0 * x447
            + 2.0 * x448
            + x460
        )
        + x163 * x456
    )
    x492 = 20.3985682659737 * x3
    x493 = x145 * x467
    x494 = x144 * x493
    x495 = x449 * x468
    x496 = 93.4779831475484 * x3
    x497 = 69.6743749058326 * x22
    x498 = x3 * x468
    x499 = x189 * x38
    x500 = x189 * x489
    x501 = x369 * x493
    x502 = x341 * x493
    x503 = -x118 - C[2]
    x504 = x141 * x503**2
    x505 = x143 + x504
    x506 = 9.12251705727742 * x178
    x507 = x147 * x359
    x508 = x119 * x505
    x509 = x141 * x503
    x510 = x24 * x509
    x511 = x508 + x510
    x512 = x189 * x505
    x513 = 35.3313566383285 * x511
    x514 = x159 * x503
    x515 = x205 + 2.0 * x514
    x516 = x0 * (x504 + x515)
    x517 = x119 * x511
    x518 = x516 + x517
    x519 = x147 * x518
    x520 = x189 * x511
    x521 = x143 + x514
    x522 = x119 * x521
    x523 = x0 * (x159 + x509)
    x524 = 2.0 * x523
    x525 = 4.0 * x143 * x503
    x526 = x524 + x525
    x527 = x0 * (2.0 * x508 + 2.0 * x522 + x526)
    x528 = x119 * x518
    x529 = x527 + x528
    x530 = x312 * x350
    x531 = x297 * x358
    x532 = x297 * x363
    x533 = 24.1359114645008 * x403
    x534 = x312 * x362
    x535 = x192 * x505
    x536 = x510 + x535
    x537 = x192 * x511
    x538 = x516 + x537
    x539 = x189 * x536
    x540 = x178 * x199
    x541 = x192 * x518
    x542 = x527 + x541
    x543 = x522 + x523
    x544 = x119 * x543
    x545 = x0 * (x142 + x515)
    x546 = 2.0 * x545
    x547 = 3.0 * x516 + x546
    x548 = x0 * (3.0 * x517 + 2.0 * x544 + x547)
    x549 = x192 * x529
    x550 = x548 + x549
    x551 = 69.6743749058326 * x518
    x552 = x128 * x189
    x553 = 120.679557322504 * x539
    x554 = 209.023124717498 * x538
    x555 = x128 * x178
    x556 = 120.679557322504 * x435
    x557 = x205 + x253 * x509
    x558 = x0 * (x504 + x557)
    x559 = x192 * x536
    x560 = x558 + x559
    x561 = x178 * x560
    x562 = x192 * x521
    x563 = 2.0 * x562
    x564 = x0 * (x508 + x526 + x535 + x563)
    x565 = x192 * x538
    x566 = x564 + x565
    x567 = x178 * x566
    x568 = x192 * x543
    x569 = 2.0 * x568
    x570 = 2.0 * x537
    x571 = x0 * (x517 + x547 + x569 + x570)
    x572 = x192 * x542
    x573 = x571 + x572
    x574 = 69.6743749058326 * x573
    x575 = 3.0 * x523
    x576 = x0 * (x161 + 3.0 * x522 + x575)
    x577 = x192 * (x544 + x545)
    x578 = 2.0 * x576 + 2.0 * x577
    x579 = x0 * (4.0 * x527 + x528 + 3.0 * x541 + x578)
    x580 = x192 * x550
    x581 = x579 + x580
    x582 = x468 * x505
    x583 = x277 * x468
    x584 = x271 * x469
    x585 = x468 * x529
    x586 = 53.9695387335403 * x539
    x587 = 120.679557322504 * x211
    x588 = x189 * x538
    x589 = 120.679557322504 * x477
    x590 = x216 * x468
    x591 = x189 * x560
    x592 = x468 * x560
    x593 = x468 * x566
    x594 = x467 * x94
    x595 = x180 * x594
    x596 = x0 * (x194 + x509)
    x597 = x194 * x503
    x598 = x192 * (x143 + x597)
    x599 = x0 * (x525 + 2.0 * x535 + 2.0 * x596 + 2.0 * x598) + x192 * x560
    x600 = x312 * x599
    x601 = x178 * x450
    x602 = x0 * (x195 + x205 + x514 + x597)
    x603 = x192 * (x523 + x562)
    x604 = 2.0 * x602 + 2.0 * x603
    x605 = x0 * (2.0 * x516 + x560 + x570 + x604)
    x606 = x192 * x566
    x607 = x605 + x606
    x608 = x0 * (x201 + x522 + x563 + x575)
    x609 = x192 * (x545 + x568)
    x610 = 2.0 * x0 * (x527 + x541 + x564 + x565 + x608 + x609)
    x611 = x192 * x573
    x612 = x610 + x611
    x613 = 24.1359114645008 * x468
    x614 = x599 * x613
    x615 = x468 * x607
    x616 = x0 * (x208 + x544 + 4.0 * x545 + 3.0 * x568)
    x617 = x192 * (x576 + x577)
    x618 = 3.0 * x571 + 3.0 * x572
    x619 = x0 * (2.0 * x548 + 2.0 * x549 + 2.0 * x616 + 2.0 * x617 + x618) + x192 * x581
    x620 = x594 * x619
    x621 = x320 * x493
    x622 = x493 * x518
    x623 = x596 + x598
    x624 = (
        x0 * (x24 * (x244 + x557) + x253 * x623 + 3.0 * x558 + 3.0 * x559) + x192 * x599
    )
    x625 = x178 * x489
    x626 = (
        x0
        * (
            x24 * (x250 + x524 + x563 + x623)
            + x253 * (x602 + x603)
            + 3.0 * x564
            + 3.0 * x565
            + x599
        )
        + x192 * x607
    )
    x627 = x139 * x467
    x628 = x594 * (
        x0
        * (
            x24 * (x257 + x546 + x569 + x604)
            + x253 * (x608 + x609)
            + 2.0 * x605
            + 2.0 * x606
            + x618
        )
        + x192 * x612
    )

    # 450 item(s)
    result[0, 0, 0] = numpy.sum(
        x91
        * x94
        * (
            x0
            * (
                x24 * (3.0 * x27 + 3.0 * x33 + x42 + x46)
                + x49 * (x47 + x48)
                + 3.0 * x64
                + 3.0 * x70
                + 3.0 * x77
                + 3.0 * x84
            )
            + x18 * x87
        )
    )
    result[0, 0, 1] = numpy.sum(x117 * x97)
    result[0, 0, 2] = numpy.sum(x117 * x120)
    result[0, 0, 3] = numpy.sum(x132 * x139)
    result[0, 0, 4] = numpy.sum(x119 * x132 * x140 * x94)
    result[0, 0, 5] = numpy.sum(x131 * x148)
    result[0, 0, 6] = numpy.sum(x153 * x156 * x90)
    result[0, 0, 7] = numpy.sum(x139 * x155 * x157)
    result[0, 0, 8] = numpy.sum(x148 * x155 * x96)
    result[0, 0, 9] = numpy.sum(x146 * x156 * x162)
    result[0, 1, 0] = numpy.sum(x164 * x165)
    result[0, 1, 1] = numpy.sum(x115 * x169 * x170)
    result[0, 1, 2] = numpy.sum(x119 * x163 * x172)
    result[0, 1, 3] = numpy.sum(x175 * x176)
    result[0, 1, 4] = numpy.sum(x119 * x176 * x177)
    result[0, 1, 5] = numpy.sum(x144 * x179 * x180)
    result[0, 1, 6] = numpy.sum(x184 * x186)
    result[0, 1, 7] = numpy.sum(x174 * x187 * x188)
    result[0, 1, 8] = numpy.sum(x129 * x169 * x190)
    result[0, 1, 9] = numpy.sum(x161 * x164 * x191)
    result[0, 2, 0] = numpy.sum(x165 * x193)
    result[0, 2, 1] = numpy.sum(x172 * x192 * x96)
    result[0, 2, 2] = numpy.sum(x171 * x197)
    result[0, 2, 3] = numpy.sum(x176 * x192 * x198)
    result[0, 2, 4] = numpy.sum(x130 * x197 * x199)
    result[0, 2, 5] = numpy.sum(x179 * x202)
    result[0, 2, 6] = numpy.sum(x152 * x186 * x192)
    result[0, 2, 7] = numpy.sum(x136 * x187 * x203)
    result[0, 2, 8] = numpy.sum(x187 * x201 * x204)
    result[0, 2, 9] = numpy.sum(x191 * x209)
    result[0, 3, 0] = numpy.sum(x170 * x212 * x86)
    result[0, 3, 1] = numpy.sum(x213 * x217)
    result[0, 3, 2] = numpy.sum(x188 * x211 * x213)
    result[0, 3, 3] = numpy.sum(x223 * x224)
    result[0, 3, 4] = numpy.sum(x119 * x224 * x225)
    result[0, 3, 5] = numpy.sum(x113 * x211 * x226)
    result[0, 3, 6] = numpy.sum(x231 * x232)
    result[0, 3, 7] = numpy.sum(x119 * x223 * x232)
    result[0, 3, 8] = numpy.sum(x111 * x216 * x226)
    result[0, 3, 9] = numpy.sum(x161 * x212 * x233)
    result[0, 4, 0] = numpy.sum(x116 * x180 * x192 * x86)
    result[0, 4, 1] = numpy.sum(x234 * x235 * x83)
    result[0, 4, 2] = numpy.sum(x197 * x236 * x83)
    result[0, 4, 3] = numpy.sum(x174 * x224 * x237)
    result[0, 4, 4] = numpy.sum(x113 * x168 * x238)
    result[0, 4, 5] = numpy.sum(x163 * x239 * x240)
    result[0, 4, 6] = numpy.sum(x192 * x232 * x241)
    result[0, 4, 7] = numpy.sum(x111 * x174 * x242)
    result[0, 4, 8] = numpy.sum(x168 * x233 * x239)
    result[0, 4, 9] = numpy.sum(x180 * x208 * x243)
    result[0, 5, 0] = numpy.sum(x178 * x246 * x86)
    result[0, 5, 1] = numpy.sum(x245 * x247 * x96)
    result[0, 5, 2] = numpy.sum(x247 * x250)
    result[0, 5, 3] = numpy.sum(x113 * x251 * x252)
    result[0, 5, 4] = numpy.sum(120.679557322504 * x240 * x250 * x96)
    result[0, 5, 5] = numpy.sum(x240 * x258)
    result[0, 5, 6] = numpy.sum(x233 * x245 * x259)
    result[0, 5, 7] = numpy.sum(x233 * x250 * x251)
    result[0, 5, 8] = numpy.sum(x243 * x258 * x96)
    result[0, 5, 9] = numpy.sum(x243 * x264)
    result[0, 6, 0] = numpy.sum(x266 * x267)
    result[0, 6, 1] = numpy.sum(x268 * x272)
    result[0, 6, 2] = numpy.sum(x188 * x273 * x81)
    result[0, 6, 3] = numpy.sum(x277 * x278)
    result[0, 6, 4] = numpy.sum(x119 * x278 * x279)
    result[0, 6, 5] = numpy.sum(x265 * x280 * x68)
    result[0, 6, 6] = numpy.sum(x109 * x283)
    result[0, 6, 7] = numpy.sum(x109 * x188 * x277)
    result[0, 6, 8] = numpy.sum(x109 * x271 * x280)
    result[0, 6, 9] = numpy.sum(x161 * x266 * x284)
    result[0, 7, 0] = numpy.sum(x211 * x267 * x285)
    result[0, 7, 1] = numpy.sum(x217 * x237 * x81)
    result[0, 7, 2] = numpy.sum(x211 * x242 * x81)
    result[0, 7, 3] = numpy.sum(x222 * x237 * x278)
    result[0, 7, 4] = numpy.sum(x216 * x238 * x68)
    result[0, 7, 5] = numpy.sum(x211 * x239 * x286)
    result[0, 7, 6] = numpy.sum(x109 * x230 * x288)
    result[0, 7, 7] = numpy.sum(x109 * x222 * x242)
    result[0, 7, 8] = numpy.sum(x216 * x239 * x284)
    result[0, 7, 9] = numpy.sum(x208 * x211 * x289)
    result[0, 8, 0] = numpy.sum(x180 * x245 * x290)
    result[0, 8, 1] = numpy.sum(x168 * x291 * x81)
    result[0, 8, 2] = numpy.sum(x250 * x292 * x81)
    result[0, 8, 3] = numpy.sum(x174 * x291 * x68)
    result[0, 8, 4] = numpy.sum(x168 * x286 * x293)
    result[0, 8, 5] = numpy.sum(x236 * x257 * x294)
    result[0, 8, 6] = numpy.sum(x109 * x184 * x295)
    result[0, 8, 7] = numpy.sum(x174 * x250 * x296)
    result[0, 8, 8] = numpy.sum(x168 * x257 * x296)
    result[0, 8, 9] = numpy.sum(x109 * x263 * x298)
    result[0, 9, 0] = numpy.sum(x290 * x300)
    result[0, 9, 1] = numpy.sum(x204 * x301 * x81)
    result[0, 9, 2] = numpy.sum(x268 * x305)
    result[0, 9, 3] = numpy.sum(x136 * x286 * x301)
    result[0, 9, 4] = numpy.sum(x199 * x305 * x68)
    result[0, 9, 5] = numpy.sum(x294 * x309)
    result[0, 9, 6] = numpy.sum(x152 * x284 * x300)
    result[0, 9, 7] = numpy.sum(x136 * x289 * x304)
    result[0, 9, 8] = numpy.sum(x109 * x204 * x309)
    result[0, 9, 9] = numpy.sum(x109 * x313)
    result[0, 10, 0] = numpy.sum(x315 * x316)
    result[0, 10, 1] = numpy.sum(x317 * x319)
    result[0, 10, 2] = numpy.sum(x119 * x314 * x319)
    result[0, 10, 3] = numpy.sum(x321 * x55)
    result[0, 10, 4] = numpy.sum(x188 * x317 * x322)
    result[0, 10, 5] = numpy.sum(x314 * x323 * x55)
    result[0, 10, 6] = numpy.sum(x324 * x325 * x52)
    result[0, 10, 7] = numpy.sum(x119 * x321 * x52)
    result[0, 10, 8] = numpy.sum(x317 * x323 * x52)
    result[0, 10, 9] = numpy.sum(x162 * x314 * x326)
    result[0, 11, 0] = numpy.sum(x192 * x266 * x316)
    result[0, 11, 1] = numpy.sum(x272 * x285 * x57)
    result[0, 11, 2] = numpy.sum(x203 * x273 * x57)
    result[0, 11, 3] = numpy.sum(x235 * x277 * x55)
    result[0, 11, 4] = numpy.sum(x203 * x271 * x327)
    result[0, 11, 5] = numpy.sum(x201 * x273 * x328)
    result[0, 11, 6] = numpy.sum(x192 * x283 * x52)
    result[0, 11, 7] = numpy.sum(x276 * x329 * x52)
    result[0, 11, 8] = numpy.sum(x201 * x271 * x330)
    result[0, 11, 9] = numpy.sum(x208 * x266 * x326)
    result[0, 12, 0] = numpy.sum(x212 * x252 * x74)
    result[0, 12, 1] = numpy.sum(x216 * x252 * x331)
    result[0, 12, 2] = numpy.sum(x211 * x331 * x332)
    result[0, 12, 3] = numpy.sum(x223 * x252 * x55)
    result[0, 12, 4] = numpy.sum(x225 * x250 * x328)
    result[0, 12, 5] = numpy.sum(x211 * x258 * x328)
    result[0, 12, 6] = numpy.sum(x230 * x333 * x52)
    result[0, 12, 7] = numpy.sum(x223 * x250 * x326)
    result[0, 12, 8] = numpy.sum(x216 * x258 * x326)
    result[0, 12, 9] = numpy.sum(x212 * x263 * x326)
    result[0, 13, 0] = numpy.sum(x163 * x300 * x334)
    result[0, 13, 1] = numpy.sum(x299 * x335 * x57)
    result[0, 13, 2] = numpy.sum(x180 * x305 * x57)
    result[0, 13, 3] = numpy.sum(x174 * x301 * x328)
    result[0, 13, 4] = numpy.sum(x168 * x327 * x336)
    result[0, 13, 5] = numpy.sum(x163 * x178 * x309 * x55)
    result[0, 13, 6] = numpy.sum(x184 * x300 * x326)
    result[0, 13, 7] = numpy.sum(x174 * x304 * x330)
    result[0, 13, 8] = numpy.sum(x169 * x308 * x326)
    result[0, 13, 9] = numpy.sum(x163 * x313 * x52)
    result[0, 14, 0] = numpy.sum(x334 * x338)
    result[0, 14, 1] = numpy.sum(x337 * x339 * x96)
    result[0, 14, 2] = numpy.sum(x339 * x340)
    result[0, 14, 3] = numpy.sum(x137 * x328 * x337)
    result[0, 14, 4] = numpy.sum(x204 * x322 * x340)
    result[0, 14, 5] = numpy.sum(x342 * x55)
    result[0, 14, 6] = numpy.sum(x152 * x326 * x338)
    result[0, 14, 7] = numpy.sum(x137 * x326 * x340)
    result[0, 14, 8] = numpy.sum(x342 * x52 * x96)
    result[0, 14, 9] = numpy.sum(x178 * x343 * x52)
    result[1, 0, 0] = numpy.sum(x325 * x346 * x351)
    result[1, 0, 1] = numpy.sum(x355 * x360)
    result[1, 0, 2] = numpy.sum(x119 * x346 * x360)
    result[1, 0, 3] = numpy.sum(x364 * x370)
    result[1, 0, 4] = numpy.sum(x188 * x364 * x371)
    result[1, 0, 5] = numpy.sum(x323 * x346 * x364)
    result[1, 0, 6] = numpy.sum(x325 * x380 * x381)
    result[1, 0, 7] = numpy.sum(x119 * x370 * x381)
    result[1, 0, 8] = numpy.sum(x323 * x355 * x381)
    result[1, 0, 9] = numpy.sum(x162 * x381 * x382)
    result[1, 1, 0] = numpy.sum(x170 * x350 * x385)
    result[1, 1, 1] = numpy.sum(x387 * x388)
    result[1, 1, 2] = numpy.sum(x188 * x358 * x389)
    result[1, 1, 3] = numpy.sum(x363 * x392)
    result[1, 1, 4] = numpy.sum(x363 * x387 * x393)
    result[1, 1, 5] = numpy.sum(x280 * x363 * x384)
    result[1, 1, 6] = numpy.sum(x401 * x402)
    result[1, 1, 7] = numpy.sum(x119 * x362 * x392)
    result[1, 1, 8] = numpy.sum(x280 * x362 * x387)
    result[1, 1, 9] = numpy.sum(x161 * x385 * x403)
    result[1, 2, 0] = numpy.sum(x185 * x192 * x346 * x350)
    result[1, 2, 1] = numpy.sum(x192 * x355 * x388)
    result[1, 2, 2] = numpy.sum(x329 * x346 * x358)
    result[1, 2, 3] = numpy.sum(x288 * x363 * x369)
    result[1, 2, 4] = numpy.sum(x203 * x363 * x404)
    result[1, 2, 5] = numpy.sum(x202 * x363 * x382)
    result[1, 2, 6] = numpy.sum(x192 * x380 * x402)
    result[1, 2, 7] = numpy.sum(x329 * x362 * x369)
    result[1, 2, 8] = numpy.sum(x202 * x355 * x403)
    result[1, 2, 9] = numpy.sum(x209 * x362 * x382)
    result[1, 3, 0] = numpy.sum(x409 * x410)
    result[1, 3, 1] = numpy.sum(x415 * x417)
    result[1, 3, 2] = numpy.sum(x119 * x408 * x417)
    result[1, 3, 3] = numpy.sum(x424 * x425)
    result[1, 3, 4] = numpy.sum(x188 * x415 * x426)
    result[1, 3, 5] = numpy.sum(x128 * x226 * x408)
    result[1, 3, 6] = numpy.sum(x170 * x433 * x434)
    result[1, 3, 7] = numpy.sum(x188 * x361 * x424)
    result[1, 3, 8] = numpy.sum(x226 * x361 * x415)
    result[1, 3, 9] = numpy.sum(x161 * x409 * x435)
    result[1, 4, 0] = numpy.sum(x192 * x389 * x410)
    result[1, 4, 1] = numpy.sum(x101 * x387 * x436)
    result[1, 4, 2] = numpy.sum(x101 * x242 * x384)
    result[1, 4, 3] = numpy.sum(x237 * x391 * x425)
    result[1, 4, 4] = numpy.sum(x128 * x238 * x387)
    result[1, 4, 5] = numpy.sum(x128 * x384 * x437)
    result[1, 4, 6] = numpy.sum(x288 * x361 * x401)
    result[1, 4, 7] = numpy.sum(x242 * x361 * x391)
    result[1, 4, 8] = numpy.sum(x239 * x387 * x435)
    result[1, 4, 9] = numpy.sum(x208 * x389 * x435)
    result[1, 5, 0] = numpy.sum(x333 * x346 * x42)
    result[1, 5, 1] = numpy.sum(x252 * x355 * x416)
    result[1, 5, 2] = numpy.sum(x332 * x346 * x416)
    result[1, 5, 3] = numpy.sum(x128 * x252 * x438)
    result[1, 5, 4] = numpy.sum(x332 * x355 * x426)
    result[1, 5, 5] = numpy.sum(x128 * x258 * x382)
    result[1, 5, 6] = numpy.sum(x333 * x361 * x380)
    result[1, 5, 7] = numpy.sum(x332 * x361 * x438)
    result[1, 5, 8] = numpy.sum(x258 * x355 * x435)
    result[1, 5, 9] = numpy.sum(x264 * x361 * x382)
    result[1, 6, 0] = numpy.sum(x40 * x443)
    result[1, 6, 1] = numpy.sum(x449 * x451)
    result[1, 6, 2] = numpy.sum(x119 * x442 * x451)
    result[1, 6, 3] = numpy.sum(x126 * x287 * x456)
    result[1, 6, 4] = numpy.sum(x126 * x393 * x449)
    result[1, 6, 5] = numpy.sum(x126 * x280 * x442)
    result[1, 6, 6] = numpy.sum(x457 * x464)
    result[1, 6, 7] = numpy.sum(x456 * x465 * x466)
    result[1, 6, 8] = numpy.sum(x144 * x449 * x470)
    result[1, 6, 9] = numpy.sum(x161 * x457 * x471)
    result[1, 7, 0] = numpy.sum(x288 * x40 * x408)
    result[1, 7, 1] = numpy.sum(x26 * x415 * x436)
    result[1, 7, 2] = numpy.sum(x242 * x26 * x408)
    result[1, 7, 3] = numpy.sum(x126 * x423 * x436)
    result[1, 7, 4] = numpy.sum(x126 * x238 * x415)
    result[1, 7, 5] = numpy.sum(x126 * x408 * x437)
    result[1, 7, 6] = numpy.sum(x18 * x433 * x472)
    result[1, 7, 7] = numpy.sum(x196 * x423 * x474)
    result[1, 7, 8] = numpy.sum(x18 * x239 * x475)
    result[1, 7, 9] = numpy.sum(x208 * x408 * x470)
    result[1, 8, 0] = numpy.sum(x295 * x384 * x40)
    result[1, 8, 1] = numpy.sum(x26 * x291 * x387)
    result[1, 8, 2] = numpy.sum(x26 * x332 * x476)
    result[1, 8, 3] = numpy.sum(x126 * x291 * x391)
    result[1, 8, 4] = numpy.sum(x293 * x387 * x477)
    result[1, 8, 5] = numpy.sum(x257 * x476 * x477)
    result[1, 8, 6] = numpy.sum(x245 * x401 * x470)
    result[1, 8, 7] = numpy.sum(x391 * x473 * x478)
    result[1, 8, 8] = numpy.sum(x257 * x473 * x479)
    result[1, 8, 9] = numpy.sum(x18 * x389 * x480)
    result[1, 9, 0] = numpy.sum(x300 * x382 * x40)
    result[1, 9, 1] = numpy.sum(x26 * x355 * x481)
    result[1, 9, 2] = numpy.sum(x336 * x346 * x450)
    result[1, 9, 3] = numpy.sum(x301 * x369 * x477)
    result[1, 9, 4] = numpy.sum(x126 * x336 * x404)
    result[1, 9, 5] = numpy.sum(x126 * x309 * x382)
    result[1, 9, 6] = numpy.sum(x18 * x380 * x482)
    result[1, 9, 7] = numpy.sum(x18 * x369 * x483)
    result[1, 9, 8] = numpy.sum(x309 * x355 * x484)
    result[1, 9, 9] = numpy.sum(x346 * x457 * x485)
    result[1, 10, 0] = numpy.sum(x325 * x38 * x487)
    result[1, 10, 1] = numpy.sum(x488 * x490)
    result[1, 10, 2] = numpy.sum(x119 * x487 * x490)
    result[1, 10, 3] = numpy.sum(x491 * x492)
    result[1, 10, 4] = numpy.sum(35.3313566383285 * x3 * x465 * x488)
    result[1, 10, 5] = numpy.sum(x3 * x487 * x494)
    result[1, 10, 6] = numpy.sum(
        x462
        * x91
        * (
            x0
            * (
                x218 * (x458 + x459)
                + x24 * (x230 + x430 + 3.0 * x452 + 3.0 * x453)
                + 3.0 * x431
                + 3.0 * x432
                + 3.0 * x454
                + 3.0 * x455
            )
            + x163 * x461
        )
    )
    result[1, 10, 7] = numpy.sum(x120 * x491)
    result[1, 10, 8] = numpy.sum(x488 * x494)
    result[1, 10, 9] = numpy.sum(x162 * x468 * x487)
    result[1, 11, 0] = numpy.sum(x192 * x38 * x443)
    result[1, 11, 1] = numpy.sum(x22 * x288 * x449)
    result[1, 11, 2] = numpy.sum(x22 * x329 * x442)
    result[1, 11, 3] = numpy.sum(x3 * x456 * x472)
    result[1, 11, 4] = numpy.sum(x196 * x495 * x496)
    result[1, 11, 5] = numpy.sum(x202 * x3 * x471)
    result[1, 11, 6] = numpy.sum(x193 * x464)
    result[1, 11, 7] = numpy.sum(x196 * x456 * x469)
    result[1, 11, 8] = numpy.sum(x202 * x495)
    result[1, 11, 9] = numpy.sum(x209 * x471)
    result[1, 12, 0] = numpy.sum(x333 * x38 * x408)
    result[1, 12, 1] = numpy.sum(x252 * x415 * x497)
    result[1, 12, 2] = numpy.sum(x332 * x408 * x497)
    result[1, 12, 3] = numpy.sum(x245 * x424 * x498)
    result[1, 12, 4] = numpy.sum(120.679557322504 * x3 * x415 * x478)
    result[1, 12, 5] = numpy.sum(x258 * x408 * x498)
    result[1, 12, 6] = numpy.sum(x246 * x433 * x468)
    result[1, 12, 7] = numpy.sum(x424 * x478)
    result[1, 12, 8] = numpy.sum(x258 * x475)
    result[1, 12, 9] = numpy.sum(x409 * x480)
    result[1, 13, 0] = numpy.sum(x300 * x384 * x499)
    result[1, 13, 1] = numpy.sum(x22 * x387 * x481)
    result[1, 13, 2] = numpy.sum(x22 * x336 * x389)
    result[1, 13, 3] = numpy.sum(x301 * x391 * x498)
    result[1, 13, 4] = numpy.sum(x304 * x479 * x496)
    result[1, 13, 5] = numpy.sum(x308 * x389 * x498)
    result[1, 13, 6] = numpy.sum(x401 * x482)
    result[1, 13, 7] = numpy.sum(x391 * x483)
    result[1, 13, 8] = numpy.sum(x309 * x479)
    result[1, 13, 9] = numpy.sum(x385 * x485)
    result[1, 14, 0] = numpy.sum(x338 * x38 * x382)
    result[1, 14, 1] = numpy.sum(x337 * x355 * x500)
    result[1, 14, 2] = numpy.sum(x340 * x382 * x489)
    result[1, 14, 3] = numpy.sum(x3 * x337 * x501)
    result[1, 14, 4] = numpy.sum(x340 * x371 * x498)
    result[1, 14, 5] = numpy.sum(x3 * x346 * x502)
    result[1, 14, 6] = numpy.sum(x338 * x380 * x468)
    result[1, 14, 7] = numpy.sum(x340 * x501)
    result[1, 14, 8] = numpy.sum(x355 * x502)
    result[1, 14, 9] = numpy.sum(x343 * x346 * x468)
    result[2, 0, 0] = numpy.sum(x351 * x505 * x506)
    result[2, 0, 1] = numpy.sum(x505 * x507 * x96)
    result[2, 0, 2] = numpy.sum(x507 * x511)
    result[2, 0, 3] = numpy.sum(x137 * x364 * x512)
    result[2, 0, 4] = numpy.sum(x204 * x364 * x513)
    result[2, 0, 5] = numpy.sum(x364 * x519)
    result[2, 0, 6] = numpy.sum(x153 * x381 * x512)
    result[2, 0, 7] = numpy.sum(x137 * x381 * x520)
    result[2, 0, 8] = numpy.sum(x381 * x519 * x96)
    result[2, 0, 9] = numpy.sum(x381 * x506 * x529)
    result[2, 1, 0] = numpy.sum(x163 * x505 * x530)
    result[2, 1, 1] = numpy.sum(x335 * x358 * x505)
    result[2, 1, 2] = numpy.sum(x163 * x511 * x531)
    result[2, 1, 3] = numpy.sum(x175 * x363 * x512)
    result[2, 1, 4] = numpy.sum(x177 * x363 * x520)
    result[2, 1, 5] = numpy.sum(x163 * x518 * x532)
    result[2, 1, 6] = numpy.sum(x184 * x505 * x533)
    result[2, 1, 7] = numpy.sum(x175 * x403 * x511)
    result[2, 1, 8] = numpy.sum(x335 * x362 * x518)
    result[2, 1, 9] = numpy.sum(x163 * x529 * x534)
    result[2, 2, 0] = numpy.sum(x530 * x536)
    result[2, 2, 1] = numpy.sum(x531 * x536 * x96)
    result[2, 2, 2] = numpy.sum(x531 * x538)
    result[2, 2, 3] = numpy.sum(x198 * x363 * x539)
    result[2, 2, 4] = numpy.sum(x363 * x538 * x540)
    result[2, 2, 5] = numpy.sum(x532 * x542)
    result[2, 2, 6] = numpy.sum(x152 * x533 * x536)
    result[2, 2, 7] = numpy.sum(x198 * x403 * x538)
    result[2, 2, 8] = numpy.sum(x297 * x362 * x542 * x96)
    result[2, 2, 9] = numpy.sum(x534 * x550)
    result[2, 3, 0] = numpy.sum(x212 * x42 * x512)
    result[2, 3, 1] = numpy.sum(x216 * x416 * x512)
    result[2, 3, 2] = numpy.sum(x211 * x416 * x520)
    result[2, 3, 3] = numpy.sum(x128 * x223 * x512)
    result[2, 3, 4] = numpy.sum(x128 * x225 * x520)
    result[2, 3, 5] = numpy.sum(x211 * x551 * x552)
    result[2, 3, 6] = numpy.sum(x231 * x435 * x505)
    result[2, 3, 7] = numpy.sum(x223 * x435 * x511)
    result[2, 3, 8] = numpy.sum(x216 * x435 * x551)
    result[2, 3, 9] = numpy.sum(x212 * x435 * x529)
    result[2, 4, 0] = numpy.sum(x298 * x42 * x536)
    result[2, 4, 1] = numpy.sum(x101 * x168 * x553)
    result[2, 4, 2] = numpy.sum(x101 * x292 * x538)
    result[2, 4, 3] = numpy.sum(x128 * x174 * x553)
    result[2, 4, 4] = numpy.sum(x168 * x552 * x554)
    result[2, 4, 5] = numpy.sum(x236 * x542 * x555)
    result[2, 4, 6] = numpy.sum(x241 * x435 * x536)
    result[2, 4, 7] = numpy.sum(x174 * x538 * x556)
    result[2, 4, 8] = numpy.sum(x168 * x542 * x556)
    result[2, 4, 9] = numpy.sum(x298 * x361 * x550)
    result[2, 5, 0] = numpy.sum(31.1593277158494 * x42 * x561)
    result[2, 5, 1] = numpy.sum(x416 * x561 * x96)
    result[2, 5, 2] = numpy.sum(x416 * x567)
    result[2, 5, 3] = numpy.sum(x251 * x552 * x560)
    result[2, 5, 4] = numpy.sum(x426 * x567 * x96)
    result[2, 5, 5] = numpy.sum(x555 * x574)
    result[2, 5, 6] = numpy.sum(x259 * x435 * x560)
    result[2, 5, 7] = numpy.sum(x251 * x435 * x566)
    result[2, 5, 8] = numpy.sum(x204 * x361 * x574)
    result[2, 5, 9] = numpy.sum(x178 * x434 * x581)
    result[2, 6, 0] = numpy.sum(x266 * x40 * x512)
    result[2, 6, 1] = numpy.sum(x271 * x450 * x512)
    result[2, 6, 2] = numpy.sum(x26 * x273 * x520)
    result[2, 6, 3] = numpy.sum(x277 * x477 * x505)
    result[2, 6, 4] = numpy.sum(x279 * x477 * x511)
    result[2, 6, 5] = numpy.sum(x273 * x477 * x518)
    result[2, 6, 6] = numpy.sum(x282 * x457 * x582)
    result[2, 6, 7] = numpy.sum(x18 * x511 * x583)
    result[2, 6, 8] = numpy.sum(x18 * x518 * x584)
    result[2, 6, 9] = numpy.sum(x18 * x266 * x585)
    result[2, 7, 0] = numpy.sum(x211 * x40 * x586)
    result[2, 7, 1] = numpy.sum(x216 * x26 * x553)
    result[2, 7, 2] = numpy.sum(x26 * x587 * x588)
    result[2, 7, 3] = numpy.sum(x222 * x536 * x589)
    result[2, 7, 4] = numpy.sum(x216 * x477 * x554)
    result[2, 7, 5] = numpy.sum(x477 * x542 * x587)
    result[2, 7, 6] = numpy.sum(x230 * x470 * x536)
    result[2, 7, 7] = numpy.sum(x222 * x474 * x538)
    result[2, 7, 8] = numpy.sum(x473 * x542 * x590)
    result[2, 7, 9] = numpy.sum(x211 * x470 * x550)
    result[2, 8, 0] = numpy.sum(x298 * x40 * x560)
    result[2, 8, 1] = numpy.sum(x234 * x26 * x591)
    result[2, 8, 2] = numpy.sum(x236 * x26 * x567)
    result[2, 8, 3] = numpy.sum(x174 * x560 * x589)
    result[2, 8, 4] = numpy.sum(209.023124717498 * x168 * x477 * x566)
    result[2, 8, 5] = numpy.sum(x126 * x292 * x573)
    result[2, 8, 6] = numpy.sum(x18 * x241 * x592)
    result[2, 8, 7] = numpy.sum(x174 * x473 * x593)
    result[2, 8, 8] = numpy.sum(x234 * x484 * x573)
    result[2, 8, 9] = numpy.sum(x18 * x581 * x595)
    result[2, 9, 0] = numpy.sum(x40 * x600)
    result[2, 9, 1] = numpy.sum(x599 * x601 * x96)
    result[2, 9, 2] = numpy.sum(x601 * x607)
    result[2, 9, 3] = numpy.sum(x198 * x477 * x599)
    result[2, 9, 4] = numpy.sum(x126 * x540 * x607)
    result[2, 9, 5] = numpy.sum(x126 * x297 * x612)
    result[2, 9, 6] = numpy.sum(x152 * x18 * x614)
    result[2, 9, 7] = numpy.sum(x18 * x198 * x615)
    result[2, 9, 8] = numpy.sum(x466 * x594 * x612 * x96)
    result[2, 9, 9] = numpy.sum(x457 * x620)
    result[2, 10, 0] = numpy.sum(x315 * x499 * x505)
    result[2, 10, 1] = numpy.sum(x317 * x500 * x505)
    result[2, 10, 2] = numpy.sum(x314 * x500 * x511)
    result[2, 10, 3] = numpy.sum(x3 * x505 * x621)
    result[2, 10, 4] = numpy.sum(x317 * x498 * x513)
    result[2, 10, 5] = numpy.sum(x3 * x314 * x622)
    result[2, 10, 6] = numpy.sum(9.12251705727742 * x324 * x582)
    result[2, 10, 7] = numpy.sum(x511 * x621)
    result[2, 10, 8] = numpy.sum(x317 * x622)
    result[2, 10, 9] = numpy.sum(x315 * x585)
    result[2, 11, 0] = numpy.sum(x266 * x499 * x536)
    result[2, 11, 1] = numpy.sum(x22 * x271 * x586)
    result[2, 11, 2] = numpy.sum(x22 * x273 * x588)
    result[2, 11, 3] = numpy.sum(x3 * x536 * x583)
    result[2, 11, 4] = numpy.sum(x279 * x498 * x538)
    result[2, 11, 5] = numpy.sum(x273 * x498 * x542)
    result[2, 11, 6] = numpy.sum(x282 * x536 * x613)
    result[2, 11, 7] = numpy.sum(x538 * x583)
    result[2, 11, 8] = numpy.sum(x542 * x584)
    result[2, 11, 9] = numpy.sum(x266 * x468 * x550)
    result[2, 12, 0] = numpy.sum(x212 * x499 * x560)
    result[2, 12, 1] = numpy.sum(x216 * x497 * x591)
    result[2, 12, 2] = numpy.sum(x189 * x211 * x497 * x566)
    result[2, 12, 3] = numpy.sum(x223 * x3 * x592)
    result[2, 12, 4] = numpy.sum(x225 * x3 * x593)
    result[2, 12, 5] = numpy.sum(x211 * x498 * x574)
    result[2, 12, 6] = numpy.sum(x231 * x592)
    result[2, 12, 7] = numpy.sum(x223 * x593)
    result[2, 12, 8] = numpy.sum(x574 * x590)
    result[2, 12, 9] = numpy.sum(x212 * x468 * x581)
    result[2, 13, 0] = numpy.sum(x163 * x38 * x600)
    result[2, 13, 1] = numpy.sum(x22 * x335 * x599)
    result[2, 13, 2] = numpy.sum(x22 * x298 * x607)
    result[2, 13, 3] = numpy.sum(x175 * x498 * x599)
    result[2, 13, 4] = numpy.sum(x177 * x3 * x615)
    result[2, 13, 5] = numpy.sum(x3 * x595 * x612)
    result[2, 13, 6] = numpy.sum(x184 * x614)
    result[2, 13, 7] = numpy.sum(x175 * x615)
    result[2, 13, 8] = numpy.sum(x169 * x468 * x612)
    result[2, 13, 9] = numpy.sum(x164 * x620)
    result[2, 14, 0] = numpy.sum(x38 * x506 * x624)
    result[2, 14, 1] = numpy.sum(x624 * x625 * x96)
    result[2, 14, 2] = numpy.sum(x625 * x626)
    result[2, 14, 3] = numpy.sum(x3 * x624 * x627)
    result[2, 14, 4] = numpy.sum(x140 * x3 * x594 * x626)
    result[2, 14, 5] = numpy.sum(x492 * x628)
    result[2, 14, 6] = numpy.sum(x153 * x468 * x624)
    result[2, 14, 7] = numpy.sum(x626 * x627)
    result[2, 14, 8] = numpy.sum(x628 * x97)
    result[2, 14, 9] = numpy.sum(
        9.12251705727742
        * x594
        * (
            x0
            * (
                x24 * (x263 + x578 + 3.0 * x608 + 3.0 * x609)
                + x253 * (x616 + x617)
                + 3.0 * x579
                + 3.0 * x580
                + 3.0 * x610
                + 3.0 * x611
            )
            + x192 * x619
        )
    )
    return result


def diag_quadrupole3d_44(ax, da, A, bx, db, B, C):
    """Cartesian 3D (gg) quadrupole moment integrals
    for operators x², y² and z². The origin is at C.

    Generated code; DO NOT modify by hand!"""

    result = numpy.zeros((3, 15, 15), dtype=float)

    x0 = 0.5 / (ax + bx)
    x1 = (ax + bx) ** (-1.0)
    x2 = -x1 * (ax * A[0] + bx * B[0])
    x3 = -x2 - B[0]
    x4 = ax * bx * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = numpy.sqrt(x1)
    x7 = 1.77245385090552 * x6
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x2 - C[0]
    x11 = x10 * x8
    x12 = x0 * (x11 + x9)
    x13 = x0 * x8
    x14 = x10 * x9
    x15 = x13 + x14
    x16 = x15 * x3
    x17 = x12 + x16
    x18 = x17 * x3
    x19 = x3**2 * x8
    x20 = 3.0 * x13
    x21 = 2.0 * x14 + x20
    x22 = x0 * (x19 + x21)
    x23 = 4.0 * x22
    x24 = -x2 - A[0]
    x25 = x17 * x24
    x26 = x0 * (3.0 * x19 + x20)
    x27 = x13 + x19
    x28 = x27 * x3
    x29 = 2.0 * x0
    x30 = x29 * x9
    x31 = x28 + x30
    x32 = x24 * x31
    x33 = x26 + x32
    x34 = x0 * (x18 + x23 + 3.0 * x25 + x33)
    x35 = 3.0 * x12
    x36 = x0 * (3.0 * x16 + x31 + x35)
    x37 = x18 + x22
    x38 = x24 * x37
    x39 = x24 * (x36 + x38)
    x40 = x3 * x31
    x41 = x0 * (5.0 * x26 + 4.0 * x32 + x40)
    x42 = x13 * x3
    x43 = 8.0 * x42
    x44 = x0 * (4.0 * x28 + x43)
    x45 = x26 + x40
    x46 = x24 * x45
    x47 = x44 + x46
    x48 = x24 * x47
    x49 = x41 + x48
    x50 = x0 * (4.0 * x18 + x23 + x45)
    x51 = x3 * x37
    x52 = x24 * (x36 + x51)
    x53 = 2.0 * x50 + 2.0 * x52
    x54 = x0 * (5.0 * x36 + 4.0 * x38 + x47 + x51)
    x55 = x24 * (x50 + x52)
    x56 = 2.0 * x24
    x57 = x10**2 * x8
    x58 = x13 + x57
    x59 = x3 * x58
    x60 = 2.0 * x12
    x61 = 4.0 * x10 * x13
    x62 = x60 + x61
    x63 = x0 * (2.0 * x16 + 2.0 * x59 + x62)
    x64 = x0 * (x21 + x57)
    x65 = x11 * x29
    x66 = x59 + x65
    x67 = x3 * x66
    x68 = x64 + x67
    x69 = x3 * x68
    x70 = x63 + x69
    x71 = x3 * x70
    x72 = x24 * x70
    x73 = 2.0 * x22
    x74 = 3.0 * x64 + x73
    x75 = x0 * (2.0 * x18 + 3.0 * x67 + x74)
    x76 = x0 * (x53 + x71 + 4.0 * x72 + 5.0 * x75)
    x77 = 2.0 * x36
    x78 = 4.0 * x63 + x77
    x79 = x0 * (2.0 * x51 + 4.0 * x69 + x78)
    x80 = x71 + x75
    x81 = x24 * x80
    x82 = x79 + x81
    x83 = x24 * x82
    x84 = 2.0 * x25
    x85 = x24 * x66
    x86 = 2.0 * x85
    x87 = x0 * (x67 + x74 + x84 + x86)
    x88 = x24 * x68
    x89 = x63 + x88
    x90 = x24 * x89
    x91 = 3.0 * x87 + 3.0 * x90
    x92 = x0 * (2.0 * x34 + 2.0 * x39 + 2.0 * x72 + 2.0 * x75 + x91)
    x93 = 2.0 * x38
    x94 = x0 * (x69 + x78 + 3.0 * x88 + x93)
    x95 = x72 + x75
    x96 = x24 * x95
    x97 = x94 + x96
    x98 = x24 * x97
    x99 = x76 + x83
    x100 = 2.0 * x0 * (x54 + x55 + x79 + x81 + 2.0 * x94 + 2.0 * x96) + x24 * x99
    x101 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x102 = da * db * numpy.sqrt(ax**5.5) * numpy.sqrt(bx**5.5)
    x103 = x101 * x102
    x104 = 6.89597470414309 * x103
    x105 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x106 = 0.564189583547756 * x1
    x107 = x105 * x106
    x108 = -x1 * (ax * A[1] + bx * B[1])
    x109 = -x108 - B[1]
    x110 = x107 * x109
    x111 = x15 * x24
    x112 = 2.0 * x111
    x113 = x24 * x27
    x114 = x113 + x30
    x115 = x0 * (x112 + x114 + x16 + x35)
    x116 = x24 * (x22 + x25)
    x117 = x0 * (3.0 * x113 + x28 + x43)
    x118 = x24 * x33
    x119 = x117 + x118
    x120 = x24 * x58
    x121 = x0 * (x112 + x120 + x59 + x62)
    x122 = x64 + x85
    x123 = x122 * x24
    x124 = 2.0 * x0 * (x115 + x116 + x121 + x123 + x63 + x88)
    x125 = x87 + x90
    x126 = x125 * x24
    x127 = x92 + x98
    x128 = 18.2450341145548 * x103
    x129 = x128 * (
        x0
        * (
            3.0 * x124
            + 3.0 * x126
            + x29 * (3.0 * x115 + 3.0 * x116 + x119 + x77 + x93)
            + x56 * (x34 + x39)
            + 3.0 * x94
            + 3.0 * x96
        )
        + x127 * x24
    )
    x130 = -x1 * (ax * A[2] + bx * B[2])
    x131 = -x130 - B[2]
    x132 = x107 * x131
    x133 = x20 + x56 * x9
    x134 = x0 * (x133 + x19)
    x135 = x114 * x24
    x136 = x134 + x135
    x137 = x24 * x9
    x138 = x11 * x24
    x139 = x0 * (x137 + x138 + x14 + x20)
    x140 = x24 * (x111 + x12)
    x141 = 2.0 * x139 + 2.0 * x140
    x142 = x11 * x56 + x20
    x143 = x0 * (x142 + x57)
    x144 = x120 + x65
    x145 = x144 * x24
    x146 = x143 + x145
    x147 = x0 * (x141 + x146 + 2.0 * x64 + x86)
    x148 = x121 + x123
    x149 = x148 * x24
    x150 = x124 + x126
    x151 = (
        x0
        * (
            2.0 * x147
            + 2.0 * x149
            + x29 * (x136 + x141 + x73 + x84)
            + x56 * (x115 + x116)
            + x91
        )
        + x150 * x24
    )
    x152 = x103 * x151
    x153 = x105 * x7
    x154 = x109**2 * x153
    x155 = x0 * x153
    x156 = x154 + x155
    x157 = 23.5542377588857 * x156
    x158 = 0.318309886183791 * x6
    x159 = x157 * x158
    x160 = 40.7971365319473 * x131
    x161 = x101 * x7
    x162 = x131**2 * x161
    x163 = x0 * x161
    x164 = x162 + x163
    x165 = 23.5542377588857 * x158
    x166 = x102 * x105
    x167 = x165 * x166
    x168 = x109 * x156
    x169 = x109 * x153
    x170 = x169 * x29
    x171 = x168 + x170
    x172 = 18.2450341145548 * x171
    x173 = x24 * x8
    x174 = x0 * (x11 + x173)
    x175 = x24 * (x13 + x138)
    x176 = x174 + x175
    x177 = x0 * (x173 + x9)
    x178 = x13 + x137
    x179 = x178 * x24
    x180 = x177 + x179
    x181 = x0 * (2.0 * x120 + 2.0 * x174 + 2.0 * x175 + x61) + x146 * x24
    x182 = x147 + x149
    x183 = x158 * (
        x0
        * (
            3.0 * x121
            + 3.0 * x123
            + x181
            + x29 * (x112 + x176 + x180 + x60)
            + x56 * (x139 + x140)
        )
        + x182 * x24
    )
    x184 = x103 * x183
    x185 = 40.7971365319473 * x109
    x186 = x131 * x164
    x187 = x131 * x161
    x188 = x187 * x29
    x189 = x186 + x188
    x190 = 18.2450341145548 * x166
    x191 = x189 * x190
    x192 = 3.0 * x155
    x193 = x0 * (3.0 * x154 + x192)
    x194 = x109 * x171
    x195 = x193 + x194
    x196 = 6.89597470414309 * x195
    x197 = x24**2 * x8
    x198 = x0 * (3.0 * x143 + 3.0 * x145 + x176 * x56 + x29 * (x142 + x197)) + x181 * x24
    x199 = x158 * x198
    x200 = x103 * x199
    x201 = 0.179587122125167 * x102
    x202 = x164 * x201
    x203 = 3.0 * x163
    x204 = x0 * (3.0 * x162 + x203)
    x205 = x131 * x189
    x206 = x204 + x205
    x207 = 6.89597470414309 * x206
    x208 = -x108 - A[1]
    x209 = x107 * x208
    x210 = x100 * x128
    x211 = x153 * x208
    x212 = x109 * x211
    x213 = x155 + x212
    x214 = 48.2718229290016 * x213
    x215 = x103 * x158
    x216 = 48.2718229290016 * x127
    x217 = x103 * x216
    x218 = x156 * x208
    x219 = x170 + x218
    x220 = 62.3186554316989 * x219
    x221 = x150 * x215
    x222 = 107.939077467081 * x213
    x223 = x158 * x166
    x224 = 62.3186554316989 * x223
    x225 = x150 * x224
    x226 = x171 * x208
    x227 = x193 + x226
    x228 = 48.2718229290016 * x227
    x229 = x215 * x228
    x230 = 107.939077467081 * x182
    x231 = x131 * x215
    x232 = x182 * x223
    x233 = 48.2718229290016 * x208
    x234 = x109 * x155
    x235 = 8.0 * x234
    x236 = x0 * (4.0 * x168 + x235)
    x237 = x195 * x208
    x238 = x236 + x237
    x239 = 18.2450341145548 * x215
    x240 = x181 * x239
    x241 = 62.3186554316989 * x181
    x242 = x189 * x201
    x243 = x158 * x190
    x244 = x181 * x243
    x245 = -x130 - A[2]
    x246 = x161 * x245
    x247 = x131 * x246
    x248 = x163 + x247
    x249 = x223 * x248
    x250 = 107.939077467081 * x249
    x251 = x164 * x245
    x252 = x188 + x251
    x253 = 48.2718229290016 * x171
    x254 = x215 * x245
    x255 = x201 * x248
    x256 = x109 * x223
    x257 = x189 * x245
    x258 = x204 + x257
    x259 = 48.2718229290016 * x258
    x260 = x201 * x252
    x261 = x131 * x163
    x262 = 8.0 * x261
    x263 = x0 * (4.0 * x186 + x262)
    x264 = x206 * x245
    x265 = x263 + x264
    x266 = x153 * x208**2
    x267 = x155 + x266
    x268 = 23.5542377588857 * x267
    x269 = x0 * (x169 + x211)
    x270 = x208 * x213
    x271 = x269 + x270
    x272 = 62.3186554316989 * x215
    x273 = x272 * x97
    x274 = 2.0 * x208
    x275 = x169 * x274 + x192
    x276 = x0 * (x154 + x275)
    x277 = x208 * x219
    x278 = x276 + x277
    x279 = 80.4530382150027 * x278
    x280 = x125 * x215
    x281 = 139.348749811665 * x271
    x282 = 80.4530382150027 * x202
    x283 = x0 * (x168 + 3.0 * x218 + x235)
    x284 = x208 * x227
    x285 = x283 + x284
    x286 = 62.3186554316989 * x285
    x287 = x215 * x286
    x288 = 139.348749811665 * x148
    x289 = 62.3186554316989 * x242
    x290 = x0 * (5.0 * x193 + x194 + 4.0 * x226)
    x291 = x208 * x238
    x292 = x290 + x291
    x293 = 23.5542377588857 * x292
    x294 = x146 * x215
    x295 = x146 * x201
    x296 = 40.7971365319473 * x103 * x245
    x297 = 139.348749811665 * x219
    x298 = 241.359114645008 * x255
    x299 = x125 * x223
    x300 = 139.348749811665 * x299
    x301 = 107.939077467081 * x227
    x302 = 241.359114645008 * x260
    x303 = 107.939077467081 * x258
    x304 = x148 * x223
    x305 = 40.7971365319473 * x238
    x306 = 107.939077467081 * x255
    x307 = 139.348749811665 * x260
    x308 = 40.7971365319473 * x265
    x309 = x208 * x223
    x310 = x161 * x245**2
    x311 = x163 + x310
    x312 = x224 * x97
    x313 = x0 * (x187 + x246)
    x314 = x245 * x248
    x315 = x313 + x314
    x316 = x201 * x311
    x317 = 80.4530382150027 * x156
    x318 = 2.0 * x245
    x319 = x187 * x318 + x203
    x320 = x0 * (x162 + x319)
    x321 = x245 * x252
    x322 = x320 + x321
    x323 = 80.4530382150027 * x322
    x324 = 62.3186554316989 * x171
    x325 = x201 * x315
    x326 = x0 * (x186 + 3.0 * x251 + x262)
    x327 = x245 * x258
    x328 = x326 + x327
    x329 = 62.3186554316989 * x328
    x330 = 23.5542377588857 * x195
    x331 = x0 * (5.0 * x204 + x205 + 4.0 * x257)
    x332 = x245 * x265
    x333 = x331 + x332
    x334 = x208 * x267 + x211 * x29
    x335 = 18.2450341145548 * x334
    x336 = x215 * x82
    x337 = 48.2718229290016 * x95
    x338 = x0 * (x266 + x275)
    x339 = x208 * x271
    x340 = x338 + x339
    x341 = x215 * x340
    x342 = 48.2718229290016 * x334
    x343 = 62.3186554316989 * x89
    x344 = 2.0 * x0 * (x218 + 2.0 * x234 + x269 + x270)
    x345 = x208 * x278
    x346 = x344 + x345
    x347 = x215 * x346
    x348 = 107.939077467081 * x340
    x349 = 62.3186554316989 * x202
    x350 = 3.0 * x276 + 3.0 * x277
    x351 = x0 * (2.0 * x193 + 2.0 * x226 + x350)
    x352 = x208 * x285
    x353 = x351 + x352
    x354 = 48.2718229290016 * x353
    x355 = x215 * x354
    x356 = 107.939077467081 * x122
    x357 = 48.2718229290016 * x242
    x358 = 2.0 * x0 * (x236 + x237 + 2.0 * x283 + 2.0 * x284) + x208 * x292
    x359 = x239 * x358
    x360 = x144 * x201
    x361 = 40.7971365319473 * x267
    x362 = 107.939077467081 * x95
    x363 = 139.348749811665 * x89
    x364 = x122 * x201
    x365 = 40.7971365319473 * x144
    x366 = x223 * x82
    x367 = 40.7971365319473 * x311
    x368 = 241.359114645008 * x325
    x369 = 241.359114645008 * x322
    x370 = 40.7971365319473 * x316
    x371 = 107.939077467081 * x325
    x372 = x245 * x311 + x246 * x29
    x373 = 18.2450341145548 * x372
    x374 = x223 * x337
    x375 = x0 * (x310 + x319)
    x376 = x245 * x315
    x377 = x375 + x376
    x378 = x156 * x201
    x379 = 107.939077467081 * x377
    x380 = 2.0 * x0 * (x251 + 2.0 * x261 + x313 + x314)
    x381 = x245 * x322
    x382 = x380 + x381
    x383 = x223 * x382
    x384 = 3.0 * x320 + 3.0 * x321
    x385 = x0 * (2.0 * x204 + 2.0 * x257 + x384)
    x386 = x245 * x328
    x387 = x385 + x386
    x388 = 48.2718229290016 * x387
    x389 = x223 * x388
    x390 = 62.3186554316989 * x382
    x391 = 2.0 * x0 * (x263 + x264 + 2.0 * x326 + 2.0 * x327) + x245 * x333
    x392 = x243 * x391
    x393 = x0 * (x192 + 3.0 * x266) + x208 * x334
    x394 = x104 * x158
    x395 = x0 * (3.0 * x269 + 3.0 * x270 + x334) + x208 * x340
    x396 = x239 * x70
    x397 = x0 * (2.0 * x338 + 2.0 * x339 + x350) + x208 * x346
    x398 = x103 * x165
    x399 = x160 * x215
    x400 = 23.5542377588857 * x202
    x401 = 3.0 * x0 * (x283 + x284 + x344 + x345) + x208 * x353
    x402 = x239 * x401
    x403 = 40.7971365319473 * x66
    x404 = 18.2450341145548 * x242
    x405 = x0 * (3.0 * x290 + 3.0 * x291 + 4.0 * x351 + 4.0 * x352) + x208 * x358
    x406 = x201 * x58
    x407 = 48.2718229290016 * x245
    x408 = x245 * x272
    x409 = 62.3186554316989 * x260
    x410 = 107.939077467081 * x66
    x411 = x201 * x66
    x412 = 48.2718229290016 * x255
    x413 = 48.2718229290016 * x406
    x414 = 62.3186554316989 * x70
    x415 = 80.4530382150027 * x68
    x416 = 139.348749811665 * x325
    x417 = x201 * x322
    x418 = 23.5542377588857 * x316
    x419 = x201 * x372
    x420 = 48.2718229290016 * x223
    x421 = 62.3186554316989 * x419
    x422 = x201 * x377
    x423 = x0 * (x203 + 3.0 * x310) + x245 * x372
    x424 = 6.89597470414309 * x223
    x425 = x243 * x70
    x426 = x0 * (3.0 * x313 + 3.0 * x314 + x372) + x245 * x377
    x427 = x201 * x423
    x428 = x185 * x223
    x429 = x0 * (2.0 * x375 + 2.0 * x376 + x384) + x245 * x382
    x430 = 3.0 * x0 * (x326 + x327 + x380 + x381) + x245 * x387
    x431 = x243 * x430
    x432 = x0 * (3.0 * x331 + 3.0 * x332 + 4.0 * x385 + 4.0 * x386) + x245 * x391
    x433 = -x108 - C[1]
    x434 = x153 * x433**2
    x435 = x155 + x434
    x436 = 3.0 * x134 + 3.0 * x135
    x437 = x0 * (2.0 * x26 + 2.0 * x32 + x436)
    x438 = x119 * x24
    x439 = 2.0 * x0 * (2.0 * x117 + 2.0 * x118 + x44 + x46) + x24 * x49
    x440 = x0 * (3.0 * x41 + 4.0 * x437 + 4.0 * x438 + 3.0 * x48) + x24 * x439
    x441 = x109 * x435
    x442 = x153 * x433
    x443 = x29 * x442
    x444 = x441 + x443
    x445 = 2.0 * x0 * (x113 + x177 + x179 + 2.0 * x42)
    x446 = x136 * x24
    x447 = x437 + x438
    x448 = 3.0 * x0 * (x117 + x118 + x445 + x446) + x24 * x447
    x449 = x239 * x448
    x450 = x169 * x433
    x451 = x192 + 2.0 * x450
    x452 = x0 * (x434 + x451)
    x453 = x109 * x444
    x454 = x452 + x453
    x455 = x0 * (x133 + x197)
    x456 = x180 * x24
    x457 = x445 + x446
    x458 = x0 * (x436 + 2.0 * x455 + 2.0 * x456) + x24 * x457
    x459 = 40.7971365319473 * x444
    x460 = x13 + x197
    x461 = x173 * x29 + x24 * x460
    x462 = x455 + x456
    x463 = x0 * (3.0 * x177 + 3.0 * x179 + x461) + x24 * x462
    x464 = x155 + x450
    x465 = x109 * x464
    x466 = x0 * (x169 + x442)
    x467 = 2.0 * x466
    x468 = 4.0 * x155 * x433
    x469 = x467 + x468
    x470 = x0 * (2.0 * x441 + 2.0 * x465 + x469)
    x471 = x109 * x454
    x472 = x470 + x471
    x473 = x239 * x472
    x474 = x465 + x466
    x475 = x109 * x474
    x476 = x0 * (x154 + x451)
    x477 = 2.0 * x476
    x478 = 3.0 * x452 + x477
    x479 = x0 * (3.0 * x453 + 2.0 * x475 + x478)
    x480 = x109 * x472
    x481 = x479 + x480
    x482 = x0 * (3.0 * x197 + x20) + x24 * x461
    x483 = x201 * x435
    x484 = x208 * x435
    x485 = x443 + x484
    x486 = 18.2450341145548 * x485
    x487 = x208 * x444
    x488 = x452 + x487
    x489 = 48.2718229290016 * x215
    x490 = x447 * x489
    x491 = x208 * x454
    x492 = x470 + x491
    x493 = x272 * x457
    x494 = 107.939077467081 * x488
    x495 = x208 * x472
    x496 = x479 + x495
    x497 = x489 * x496
    x498 = 107.939077467081 * x462
    x499 = x475 + x476
    x500 = x109 * x499
    x501 = 3.0 * x466
    x502 = x0 * (x171 + 3.0 * x465 + x501)
    x503 = 2.0 * x502
    x504 = 4.0 * x470 + x503
    x505 = x0 * (4.0 * x471 + 2.0 * x500 + x504)
    x506 = x208 * x481
    x507 = x505 + x506
    x508 = x239 * x461
    x509 = x201 * x461
    x510 = 18.2450341145548 * x483
    x511 = x192 + x274 * x442
    x512 = x0 * (x434 + x511)
    x513 = x208 * x485
    x514 = x512 + x513
    x515 = 23.5542377588857 * x514
    x516 = x215 * x49
    x517 = x208 * x464
    x518 = 2.0 * x517
    x519 = x0 * (x441 + x469 + x484 + x518)
    x520 = x208 * x488
    x521 = x519 + x520
    x522 = 62.3186554316989 * x119
    x523 = x215 * x522
    x524 = 80.4530382150027 * x136
    x525 = x208 * x474
    x526 = 2.0 * x525
    x527 = 2.0 * x487
    x528 = x0 * (x453 + x478 + x526 + x527)
    x529 = x208 * x492
    x530 = x528 + x529
    x531 = x215 * x530
    x532 = 139.348749811665 * x136
    x533 = x208 * x499
    x534 = 2.0 * x533
    x535 = x0 * (x471 + 3.0 * x491 + x504 + x534)
    x536 = x208 * x496
    x537 = x535 + x536
    x538 = x272 * x537
    x539 = 139.348749811665 * x180
    x540 = 4.0 * x476
    x541 = x0 * (x195 + 4.0 * x475 + x540)
    x542 = x208 * (x500 + x502)
    x543 = 2.0 * x541 + 2.0 * x542
    x544 = x0 * (5.0 * x479 + x480 + 4.0 * x495 + x543)
    x545 = x208 * x507
    x546 = x544 + x545
    x547 = x201 * x460
    x548 = 40.7971365319473 * x485
    x549 = 107.939077467081 * x496
    x550 = x201 * x303
    x551 = 40.7971365319473 * x254
    x552 = 62.3186554316989 * x472
    x553 = 80.4530382150027 * x417
    x554 = 23.5542377588857 * x483
    x555 = x0 * (x211 + x442)
    x556 = x211 * x433
    x557 = x208 * (x155 + x556)
    x558 = x0 * (x468 + 2.0 * x484 + 2.0 * x555 + 2.0 * x557) + x208 * x514
    x559 = x239 * x558
    x560 = x0 * (x192 + x212 + x450 + x556)
    x561 = x208 * (x466 + x517)
    x562 = 2.0 * x560 + 2.0 * x561
    x563 = x0 * (2.0 * x452 + x514 + x527 + x562)
    x564 = x208 * x521
    x565 = x563 + x564
    x566 = 48.2718229290016 * x33
    x567 = x215 * x566
    x568 = 62.3186554316989 * x114
    x569 = x0 * (x219 + x465 + x501 + x518)
    x570 = x208 * (x476 + x525)
    x571 = 2.0 * x0 * (x470 + x491 + x519 + x520 + x569 + x570)
    x572 = x208 * x530
    x573 = x571 + x572
    x574 = x215 * x573
    x575 = 107.939077467081 * x114
    x576 = x0 * (x227 + x475 + 3.0 * x525 + x540)
    x577 = x208 * (x502 + x533)
    x578 = 3.0 * x528 + 3.0 * x529
    x579 = x0 * (2.0 * x479 + 2.0 * x495 + 2.0 * x576 + 2.0 * x577 + x578)
    x580 = x208 * x537
    x581 = x579 + x580
    x582 = 48.2718229290016 * x178
    x583 = 107.939077467081 * x178
    x584 = x0 * (x238 + x500 + 5.0 * x502 + 4.0 * x533)
    x585 = x208 * (x541 + x542)
    x586 = 2.0 * x0 * (x505 + x506 + 2.0 * x535 + 2.0 * x536 + x584 + x585) + x208 * x546
    x587 = x106 * x5
    x588 = x128 * x587
    x589 = x586 * x588
    x590 = x24 * x587
    x591 = x103 * x581
    x592 = x102 * x5
    x593 = x158 * x592
    x594 = 62.3186554316989 * x593
    x595 = x573 * x594
    x596 = 48.2718229290016 * x593
    x597 = 18.2450341145548 * x593
    x598 = x558 * x597
    x599 = 107.939077467081 * x33
    x600 = 139.348749811665 * x114
    x601 = x24 * x593
    x602 = 107.939077467081 * x248
    x603 = x530 * x593
    x604 = 139.348749811665 * x24
    x605 = x521 * x593
    x606 = x178 * x201
    x607 = x201 * x583
    x608 = x507 * x593
    x609 = x492 * x593
    x610 = x333 * x593
    x611 = x377 * x596
    x612 = x388 * x593
    x613 = 6.89597470414309 * x45
    x614 = x555 + x557
    x615 = (
        x0 * (x274 * x614 + x29 * (x266 + x511) + 3.0 * x512 + 3.0 * x513) + x208 * x558
    )
    x616 = x215 * x615
    x617 = 18.2450341145548 * x31
    x618 = (
        x0
        * (
            x274 * (x560 + x561)
            + x29 * (x271 + x467 + x518 + x614)
            + 3.0 * x519
            + 3.0 * x520
            + x558
        )
        + x208 * x565
    )
    x619 = x215 * x618
    x620 = (
        x0
        * (
            x274 * (x569 + x570)
            + x29 * (x278 + x477 + x526 + x562)
            + 2.0 * x563
            + 2.0 * x564
            + x578
        )
        + x208 * x573
    )
    x621 = 23.5542377588857 * x27
    x622 = x588 * (
        x0
        * (
            x274 * (x576 + x577)
            + x29 * (x285 + x503 + x534 + 3.0 * x569 + 3.0 * x570)
            + 3.0 * x535
            + 3.0 * x536
            + 3.0 * x571
            + 3.0 * x572
        )
        + x208 * x581
    )
    x623 = x3 * x587
    x624 = x3 * x593
    x625 = 40.7971365319473 * x624
    x626 = x189 * x597
    x627 = x165 * x592
    x628 = 48.2718229290016 * x31
    x629 = x565 * x593
    x630 = 107.939077467081 * x3
    x631 = 62.3186554316989 * x31
    x632 = 80.4530382150027 * x27
    x633 = x537 * x594
    x634 = 139.348749811665 * x3
    x635 = x201 * x485
    x636 = x201 * x27
    x637 = x472 * x597
    x638 = x430 * x597
    x639 = 6.89597470414309 * x593
    x640 = -x130 - C[2]
    x641 = x161 * x640**2
    x642 = x163 + x641
    x643 = x243 * x448
    x644 = x131 * x642
    x645 = x161 * x640
    x646 = x29 * x645
    x647 = x644 + x646
    x648 = x201 * x642
    x649 = 40.7971365319473 * x647
    x650 = x187 * x640
    x651 = x203 + 2.0 * x650
    x652 = x0 * (x641 + x651)
    x653 = x131 * x647
    x654 = x652 + x653
    x655 = x163 + x650
    x656 = x131 * x655
    x657 = x0 * (x187 + x645)
    x658 = 2.0 * x657
    x659 = 4.0 * x163 * x640
    x660 = x658 + x659
    x661 = x0 * (2.0 * x644 + 2.0 * x656 + x660)
    x662 = x131 * x654
    x663 = x661 + x662
    x664 = x243 * x663
    x665 = x201 * x482
    x666 = x656 + x657
    x667 = x131 * x666
    x668 = x0 * (x162 + x651)
    x669 = 2.0 * x668
    x670 = 3.0 * x652 + x669
    x671 = x0 * (3.0 * x653 + 2.0 * x667 + x670)
    x672 = x131 * x663
    x673 = x671 + x672
    x674 = x243 * x439
    x675 = x420 * x447
    x676 = x201 * x647
    x677 = x224 * x457
    x678 = x201 * x462
    x679 = x420 * x462
    x680 = 18.2450341145548 * x509
    x681 = x243 * x461
    x682 = x245 * x642
    x683 = x646 + x682
    x684 = x245 * x647
    x685 = x652 + x684
    x686 = 62.3186554316989 * x378
    x687 = 107.939077467081 * x685
    x688 = x245 * x654
    x689 = x661 + x688
    x690 = x245 * x663
    x691 = x671 + x690
    x692 = x667 + x668
    x693 = x131 * x692
    x694 = 3.0 * x657
    x695 = x0 * (x189 + 3.0 * x656 + x694)
    x696 = 2.0 * x695
    x697 = 4.0 * x661 + x696
    x698 = x0 * (4.0 * x662 + 2.0 * x693 + x697)
    x699 = x245 * x673
    x700 = x698 + x699
    x701 = x201 * x654
    x702 = 62.3186554316989 * x663
    x703 = x180 * x201
    x704 = 40.7971365319473 * x309
    x705 = x201 * x683
    x706 = 241.359114645008 * x685
    x707 = x136 * x201
    x708 = 241.359114645008 * x689
    x709 = 107.939077467081 * x691
    x710 = x203 + x318 * x645
    x711 = x0 * (x641 + x710)
    x712 = x245 * x683
    x713 = x711 + x712
    x714 = x223 * x522
    x715 = x245 * x655
    x716 = 2.0 * x715
    x717 = x0 * (x644 + x660 + x682 + x716)
    x718 = x245 * x685
    x719 = x717 + x718
    x720 = x245 * x666
    x721 = 2.0 * x720
    x722 = 2.0 * x684
    x723 = x0 * (x653 + x670 + x721 + x722)
    x724 = x245 * x689
    x725 = x723 + x724
    x726 = x223 * x725
    x727 = x245 * x692
    x728 = 2.0 * x727
    x729 = x0 * (x662 + 3.0 * x688 + x697 + x728)
    x730 = x245 * x691
    x731 = x729 + x730
    x732 = x224 * x731
    x733 = 4.0 * x668
    x734 = x0 * (x206 + 4.0 * x667 + x733)
    x735 = x245 * (x693 + x695)
    x736 = 2.0 * x734 + 2.0 * x735
    x737 = x0 * (5.0 * x671 + x672 + 4.0 * x690 + x736)
    x738 = x245 * x700
    x739 = x737 + x738
    x740 = x358 * x597
    x741 = x354 * x593
    x742 = x346 * x594
    x743 = x340 * x596
    x744 = x335 * x593
    x745 = x201 * x267
    x746 = x114 * x201
    x747 = x201 * x713
    x748 = 241.359114645008 * x719
    x749 = x593 * x713
    x750 = x593 * x719
    x751 = x593 * x725
    x752 = x107 * x592
    x753 = x24 * x752
    x754 = x0 * (x246 + x645)
    x755 = x246 * x640
    x756 = x245 * (x163 + x755)
    x757 = x0 * (x659 + 2.0 * x682 + 2.0 * x754 + 2.0 * x756) + x245 * x713
    x758 = x243 * x757
    x759 = x223 * x566
    x760 = x0 * (x203 + x247 + x650 + x755)
    x761 = x245 * (x657 + x715)
    x762 = 2.0 * x760 + 2.0 * x761
    x763 = x0 * (2.0 * x652 + x713 + x722 + x762)
    x764 = x245 * x719
    x765 = x763 + x764
    x766 = x0 * (x252 + x656 + x694 + x716)
    x767 = x245 * (x668 + x720)
    x768 = 2.0 * x0 * (x661 + x688 + x717 + x718 + x766 + x767)
    x769 = x245 * x725
    x770 = x768 + x769
    x771 = x223 * x770
    x772 = x0 * (x258 + x667 + 3.0 * x720 + x733)
    x773 = x245 * (x695 + x727)
    x774 = 3.0 * x723 + 3.0 * x724
    x775 = x0 * (2.0 * x671 + 2.0 * x690 + 2.0 * x772 + 2.0 * x773 + x774)
    x776 = x245 * x731
    x777 = x775 + x776
    x778 = x597 * x757
    x779 = x593 * x765
    x780 = x0 * (x265 + x693 + 5.0 * x695 + 4.0 * x727)
    x781 = x245 * (x734 + x735)
    x782 = 2.0 * x0 * (x698 + x699 + 2.0 * x729 + 2.0 * x730 + x780 + x781) + x245 * x739
    x783 = 18.2450341145548 * x752
    x784 = x782 * x783
    x785 = x401 * x597
    x786 = x597 * x663
    x787 = x201 * x31
    x788 = 62.3186554316989 * x636
    x789 = x594 * x731
    x790 = x593 * x770
    x791 = x754 + x756
    x792 = (
        x0 * (x29 * (x310 + x710) + x318 * x791 + 3.0 * x711 + 3.0 * x712) + x245 * x757
    )
    x793 = x223 * x792
    x794 = (
        x0
        * (
            x29 * (x315 + x658 + x716 + x791)
            + x318 * (x760 + x761)
            + 3.0 * x717
            + 3.0 * x718
            + x757
        )
        + x245 * x765
    )
    x795 = x223 * x794
    x796 = (
        x0
        * (
            x29 * (x322 + x669 + x721 + x762)
            + x318 * (x766 + x767)
            + 2.0 * x763
            + 2.0 * x764
            + x774
        )
        + x245 * x770
    )
    x797 = x593 * x792
    x798 = x593 * x794
    x799 = 40.7971365319473 * x3
    x800 = x592 * x796
    x801 = x783 * (
        x0
        * (
            x29 * (x328 + x696 + x728 + 3.0 * x766 + 3.0 * x767)
            + x318 * (x772 + x773)
            + 3.0 * x729
            + 3.0 * x730
            + 3.0 * x768
            + 3.0 * x769
        )
        + x245 * x777
    )

    # 675 item(s)
    result[0, 0, 0] = numpy.sum(
        x104
        * x107
        * (
            x0
            * (
                x29 * (4.0 * x34 + 4.0 * x39 + x49 + x53)
                + x56 * (x54 + x55)
                + 3.0 * x76
                + 3.0 * x83
                + 4.0 * x92
                + 4.0 * x98
            )
            + x100 * x24
        )
    )
    result[0, 0, 1] = numpy.sum(x110 * x129)
    result[0, 0, 2] = numpy.sum(x129 * x132)
    result[0, 0, 3] = numpy.sum(x152 * x159)
    result[0, 0, 4] = numpy.sum(x110 * x152 * x160)
    result[0, 0, 5] = numpy.sum(x151 * x164 * x167)
    result[0, 0, 6] = numpy.sum(x172 * x184)
    result[0, 0, 7] = numpy.sum(x156 * x160 * x184)
    result[0, 0, 8] = numpy.sum(x164 * x166 * x183 * x185)
    result[0, 0, 9] = numpy.sum(x183 * x191)
    result[0, 0, 10] = numpy.sum(x196 * x200)
    result[0, 0, 11] = numpy.sum(x131 * x172 * x200)
    result[0, 0, 12] = numpy.sum(x157 * x198 * x202)
    result[0, 0, 13] = numpy.sum(x109 * x191 * x199)
    result[0, 0, 14] = numpy.sum(x166 * x199 * x207)
    result[0, 1, 0] = numpy.sum(x209 * x210)
    result[0, 1, 1] = numpy.sum(x127 * x214 * x215)
    result[0, 1, 2] = numpy.sum(x132 * x208 * x217)
    result[0, 1, 3] = numpy.sum(x220 * x221)
    result[0, 1, 4] = numpy.sum(x131 * x221 * x222)
    result[0, 1, 5] = numpy.sum(x164 * x208 * x225)
    result[0, 1, 6] = numpy.sum(x182 * x229)
    result[0, 1, 7] = numpy.sum(x219 * x230 * x231)
    result[0, 1, 8] = numpy.sum(x202 * x213 * x230)
    result[0, 1, 9] = numpy.sum(x189 * x232 * x233)
    result[0, 1, 10] = numpy.sum(x238 * x240)
    result[0, 1, 11] = numpy.sum(x131 * x181 * x229)
    result[0, 1, 12] = numpy.sum(x202 * x219 * x241)
    result[0, 1, 13] = numpy.sum(x181 * x214 * x242)
    result[0, 1, 14] = numpy.sum(x206 * x208 * x244)
    result[0, 2, 0] = numpy.sum(x107 * x210 * x245)
    result[0, 2, 1] = numpy.sum(x110 * x217 * x245)
    result[0, 2, 2] = numpy.sum(x216 * x249)
    result[0, 2, 3] = numpy.sum(62.3186554316989 * x156 * x221 * x245)
    result[0, 2, 4] = numpy.sum(x109 * x150 * x250)
    result[0, 2, 5] = numpy.sum(x225 * x252)
    result[0, 2, 6] = numpy.sum(x182 * x253 * x254)
    result[0, 2, 7] = numpy.sum(x156 * x230 * x255)
    result[0, 2, 8] = numpy.sum(x230 * x252 * x256)
    result[0, 2, 9] = numpy.sum(x232 * x259)
    result[0, 2, 10] = numpy.sum(x195 * x240 * x245)
    result[0, 2, 11] = numpy.sum(x181 * x253 * x255)
    result[0, 2, 12] = numpy.sum(x156 * x241 * x260)
    result[0, 2, 13] = numpy.sum(x181 * x256 * x259)
    result[0, 2, 14] = numpy.sum(x244 * x265)
    result[0, 3, 0] = numpy.sum(x215 * x268 * x99)
    result[0, 3, 1] = numpy.sum(x271 * x273)
    result[0, 3, 2] = numpy.sum(x131 * x267 * x273)
    result[0, 3, 3] = numpy.sum(x279 * x280)
    result[0, 3, 4] = numpy.sum(x131 * x280 * x281)
    result[0, 3, 5] = numpy.sum(x125 * x267 * x282)
    result[0, 3, 6] = numpy.sum(x148 * x287)
    result[0, 3, 7] = numpy.sum(x231 * x278 * x288)
    result[0, 3, 8] = numpy.sum(x202 * x271 * x288)
    result[0, 3, 9] = numpy.sum(x148 * x267 * x289)
    result[0, 3, 10] = numpy.sum(x293 * x294)
    result[0, 3, 11] = numpy.sum(x131 * x146 * x287)
    result[0, 3, 12] = numpy.sum(x146 * x278 * x282)
    result[0, 3, 13] = numpy.sum(x146 * x271 * x289)
    result[0, 3, 14] = numpy.sum(x206 * x268 * x295)
    result[0, 4, 0] = numpy.sum(x209 * x296 * x99)
    result[0, 4, 1] = numpy.sum(x222 * x254 * x97)
    result[0, 4, 2] = numpy.sum(x208 * x250 * x97)
    result[0, 4, 3] = numpy.sum(x245 * x280 * x297)
    result[0, 4, 4] = numpy.sum(x125 * x213 * x298)
    result[0, 4, 5] = numpy.sum(x208 * x252 * x300)
    result[0, 4, 6] = numpy.sum(x148 * x254 * x301)
    result[0, 4, 7] = numpy.sum(x148 * x219 * x298)
    result[0, 4, 8] = numpy.sum(x148 * x213 * x302)
    result[0, 4, 9] = numpy.sum(x208 * x303 * x304)
    result[0, 4, 10] = numpy.sum(x245 * x294 * x305)
    result[0, 4, 11] = numpy.sum(x146 * x227 * x306)
    result[0, 4, 12] = numpy.sum(x146 * x219 * x307)
    result[0, 4, 13] = numpy.sum(x213 * x295 * x303)
    result[0, 4, 14] = numpy.sum(x146 * x308 * x309)
    result[0, 5, 0] = numpy.sum(x167 * x311 * x99)
    result[0, 5, 1] = numpy.sum(x109 * x311 * x312)
    result[0, 5, 2] = numpy.sum(x312 * x315)
    result[0, 5, 3] = numpy.sum(x125 * x316 * x317)
    result[0, 5, 4] = numpy.sum(x109 * x300 * x315)
    result[0, 5, 5] = numpy.sum(x299 * x323)
    result[0, 5, 6] = numpy.sum(x148 * x316 * x324)
    result[0, 5, 7] = numpy.sum(x156 * x288 * x325)
    result[0, 5, 8] = numpy.sum(x256 * x288 * x322)
    result[0, 5, 9] = numpy.sum(x304 * x329)
    result[0, 5, 10] = numpy.sum(x295 * x311 * x330)
    result[0, 5, 11] = numpy.sum(x295 * x315 * x324)
    result[0, 5, 12] = numpy.sum(x295 * x317 * x322)
    result[0, 5, 13] = numpy.sum(x146 * x256 * x329)
    result[0, 5, 14] = numpy.sum(x146 * x167 * x333)
    result[0, 6, 0] = numpy.sum(x335 * x336)
    result[0, 6, 1] = numpy.sum(x337 * x341)
    result[0, 6, 2] = numpy.sum(x231 * x342 * x95)
    result[0, 6, 3] = numpy.sum(x343 * x347)
    result[0, 6, 4] = numpy.sum(x231 * x348 * x89)
    result[0, 6, 5] = numpy.sum(x334 * x349 * x89)
    result[0, 6, 6] = numpy.sum(x122 * x355)
    result[0, 6, 7] = numpy.sum(x131 * x347 * x356)
    result[0, 6, 8] = numpy.sum(x202 * x340 * x356)
    result[0, 6, 9] = numpy.sum(x122 * x334 * x357)
    result[0, 6, 10] = numpy.sum(x144 * x359)
    result[0, 6, 11] = numpy.sum(x131 * x144 * x355)
    result[0, 6, 12] = numpy.sum(x144 * x346 * x349)
    result[0, 6, 13] = numpy.sum(x144 * x340 * x357)
    result[0, 6, 14] = numpy.sum(x206 * x335 * x360)
    result[0, 7, 0] = numpy.sum(x245 * x336 * x361)
    result[0, 7, 1] = numpy.sum(x254 * x271 * x362)
    result[0, 7, 2] = numpy.sum(x267 * x306 * x95)
    result[0, 7, 3] = numpy.sum(x254 * x278 * x363)
    result[0, 7, 4] = numpy.sum(x271 * x298 * x89)
    result[0, 7, 5] = numpy.sum(x267 * x307 * x89)
    result[0, 7, 6] = numpy.sum(x254 * x285 * x356)
    result[0, 7, 7] = numpy.sum(x122 * x278 * x298)
    result[0, 7, 8] = numpy.sum(x122 * x271 * x302)
    result[0, 7, 9] = numpy.sum(x267 * x303 * x364)
    result[0, 7, 10] = numpy.sum(x254 * x292 * x365)
    result[0, 7, 11] = numpy.sum(x144 * x285 * x306)
    result[0, 7, 12] = numpy.sum(x144 * x278 * x307)
    result[0, 7, 13] = numpy.sum(x271 * x303 * x360)
    result[0, 7, 14] = numpy.sum(x267 * x308 * x360)
    result[0, 8, 0] = numpy.sum(x208 * x366 * x367)
    result[0, 8, 1] = numpy.sum(x222 * x316 * x95)
    result[0, 8, 2] = numpy.sum(x309 * x315 * x362)
    result[0, 8, 3] = numpy.sum(x297 * x316 * x89)
    result[0, 8, 4] = numpy.sum(x213 * x368 * x89)
    result[0, 8, 5] = numpy.sum(x309 * x322 * x363)
    result[0, 8, 6] = numpy.sum(x227 * x316 * x356)
    result[0, 8, 7] = numpy.sum(x122 * x219 * x368)
    result[0, 8, 8] = numpy.sum(x213 * x364 * x369)
    result[0, 8, 9] = numpy.sum(x309 * x328 * x356)
    result[0, 8, 10] = numpy.sum(x144 * x238 * x370)
    result[0, 8, 11] = numpy.sum(x144 * x227 * x371)
    result[0, 8, 12] = numpy.sum(x297 * x322 * x360)
    result[0, 8, 13] = numpy.sum(x222 * x328 * x360)
    result[0, 8, 14] = numpy.sum(x309 * x333 * x365)
    result[0, 9, 0] = numpy.sum(x366 * x373)
    result[0, 9, 1] = numpy.sum(x109 * x372 * x374)
    result[0, 9, 2] = numpy.sum(x374 * x377)
    result[0, 9, 3] = numpy.sum(x343 * x372 * x378)
    result[0, 9, 4] = numpy.sum(x256 * x379 * x89)
    result[0, 9, 5] = numpy.sum(x343 * x383)
    result[0, 9, 6] = numpy.sum(x253 * x364 * x372)
    result[0, 9, 7] = numpy.sum(x356 * x377 * x378)
    result[0, 9, 8] = numpy.sum(x109 * x356 * x383)
    result[0, 9, 9] = numpy.sum(x122 * x389)
    result[0, 9, 10] = numpy.sum(x195 * x360 * x373)
    result[0, 9, 11] = numpy.sum(x253 * x360 * x377)
    result[0, 9, 12] = numpy.sum(x156 * x360 * x390)
    result[0, 9, 13] = numpy.sum(x109 * x144 * x389)
    result[0, 9, 14] = numpy.sum(x144 * x392)
    result[0, 10, 0] = numpy.sum(x393 * x394 * x80)
    result[0, 10, 1] = numpy.sum(x395 * x396)
    result[0, 10, 2] = numpy.sum(x131 * x393 * x396)
    result[0, 10, 3] = numpy.sum(x397 * x398 * x68)
    result[0, 10, 4] = numpy.sum(x395 * x399 * x68)
    result[0, 10, 5] = numpy.sum(x393 * x400 * x68)
    result[0, 10, 6] = numpy.sum(x402 * x66)
    result[0, 10, 7] = numpy.sum(x231 * x397 * x403)
    result[0, 10, 8] = numpy.sum(x202 * x395 * x403)
    result[0, 10, 9] = numpy.sum(x393 * x404 * x66)
    result[0, 10, 10] = numpy.sum(x394 * x405 * x58)
    result[0, 10, 11] = numpy.sum(x131 * x402 * x58)
    result[0, 10, 12] = numpy.sum(x397 * x400 * x58)
    result[0, 10, 13] = numpy.sum(x395 * x404 * x58)
    result[0, 10, 14] = numpy.sum(x207 * x393 * x406)
    result[0, 11, 0] = numpy.sum(x254 * x335 * x80)
    result[0, 11, 1] = numpy.sum(x341 * x407 * x70)
    result[0, 11, 2] = numpy.sum(x255 * x342 * x70)
    result[0, 11, 3] = numpy.sum(x346 * x408 * x68)
    result[0, 11, 4] = numpy.sum(x306 * x340 * x68)
    result[0, 11, 5] = numpy.sum(x334 * x409 * x68)
    result[0, 11, 6] = numpy.sum(x245 * x355 * x66)
    result[0, 11, 7] = numpy.sum(x306 * x346 * x66)
    result[0, 11, 8] = numpy.sum(x260 * x340 * x410)
    result[0, 11, 9] = numpy.sum(x258 * x342 * x411)
    result[0, 11, 10] = numpy.sum(x245 * x359 * x58)
    result[0, 11, 11] = numpy.sum(x353 * x412 * x58)
    result[0, 11, 12] = numpy.sum(x346 * x409 * x58)
    result[0, 11, 13] = numpy.sum(x258 * x340 * x413)
    result[0, 11, 14] = numpy.sum(x265 * x335 * x406)
    result[0, 12, 0] = numpy.sum(x268 * x316 * x80)
    result[0, 12, 1] = numpy.sum(x271 * x316 * x414)
    result[0, 12, 2] = numpy.sum(x267 * x325 * x414)
    result[0, 12, 3] = numpy.sum(x278 * x316 * x415)
    result[0, 12, 4] = numpy.sum(x271 * x416 * x68)
    result[0, 12, 5] = numpy.sum(x267 * x415 * x417)
    result[0, 12, 6] = numpy.sum(x286 * x316 * x66)
    result[0, 12, 7] = numpy.sum(x278 * x416 * x66)
    result[0, 12, 8] = numpy.sum(x281 * x322 * x411)
    result[0, 12, 9] = numpy.sum(x267 * x329 * x411)
    result[0, 12, 10] = numpy.sum(x292 * x418 * x58)
    result[0, 12, 11] = numpy.sum(x286 * x325 * x58)
    result[0, 12, 12] = numpy.sum(x279 * x322 * x406)
    result[0, 12, 13] = numpy.sum(x271 * x329 * x406)
    result[0, 12, 14] = numpy.sum(x268 * x333 * x406)
    result[0, 13, 0] = numpy.sum(x309 * x373 * x80)
    result[0, 13, 1] = numpy.sum(x214 * x419 * x70)
    result[0, 13, 2] = numpy.sum(x208 * x377 * x420 * x70)
    result[0, 13, 3] = numpy.sum(x219 * x421 * x68)
    result[0, 13, 4] = numpy.sum(x222 * x422 * x68)
    result[0, 13, 5] = numpy.sum(x309 * x390 * x68)
    result[0, 13, 6] = numpy.sum(x228 * x372 * x411)
    result[0, 13, 7] = numpy.sum(x219 * x410 * x422)
    result[0, 13, 8] = numpy.sum(x222 * x382 * x411)
    result[0, 13, 9] = numpy.sum(x208 * x389 * x66)
    result[0, 13, 10] = numpy.sum(x238 * x373 * x406)
    result[0, 13, 11] = numpy.sum(x227 * x377 * x413)
    result[0, 13, 12] = numpy.sum(x219 * x390 * x406)
    result[0, 13, 13] = numpy.sum(x214 * x387 * x406)
    result[0, 13, 14] = numpy.sum(x208 * x392 * x58)
    result[0, 14, 0] = numpy.sum(x423 * x424 * x80)
    result[0, 14, 1] = numpy.sum(x109 * x423 * x425)
    result[0, 14, 2] = numpy.sum(x425 * x426)
    result[0, 14, 3] = numpy.sum(x157 * x427 * x68)
    result[0, 14, 4] = numpy.sum(x426 * x428 * x68)
    result[0, 14, 5] = numpy.sum(x167 * x429 * x68)
    result[0, 14, 6] = numpy.sum(x172 * x411 * x423)
    result[0, 14, 7] = numpy.sum(x378 * x403 * x426)
    result[0, 14, 8] = numpy.sum(x256 * x403 * x429)
    result[0, 14, 9] = numpy.sum(x431 * x66)
    result[0, 14, 10] = numpy.sum(x196 * x406 * x423)
    result[0, 14, 11] = numpy.sum(x172 * x406 * x426)
    result[0, 14, 12] = numpy.sum(x157 * x406 * x429)
    result[0, 14, 13] = numpy.sum(x109 * x431 * x58)
    result[0, 14, 14] = numpy.sum(x424 * x432 * x58)
    result[1, 0, 0] = numpy.sum(x394 * x435 * x440)
    result[1, 0, 1] = numpy.sum(x444 * x449)
    result[1, 0, 2] = numpy.sum(x131 * x435 * x449)
    result[1, 0, 3] = numpy.sum(x398 * x454 * x458)
    result[1, 0, 4] = numpy.sum(x231 * x458 * x459)
    result[1, 0, 5] = numpy.sum(x400 * x435 * x458)
    result[1, 0, 6] = numpy.sum(x463 * x473)
    result[1, 0, 7] = numpy.sum(x399 * x454 * x463)
    result[1, 0, 8] = numpy.sum(x202 * x459 * x463)
    result[1, 0, 9] = numpy.sum(x404 * x435 * x463)
    result[1, 0, 10] = numpy.sum(x394 * x481 * x482)
    result[1, 0, 11] = numpy.sum(x131 * x473 * x482)
    result[1, 0, 12] = numpy.sum(x400 * x454 * x482)
    result[1, 0, 13] = numpy.sum(x404 * x444 * x482)
    result[1, 0, 14] = numpy.sum(x207 * x482 * x483)
    result[1, 1, 0] = numpy.sum(x215 * x439 * x486)
    result[1, 1, 1] = numpy.sum(x488 * x490)
    result[1, 1, 2] = numpy.sum(x131 * x485 * x490)
    result[1, 1, 3] = numpy.sum(x492 * x493)
    result[1, 1, 4] = numpy.sum(x231 * x457 * x494)
    result[1, 1, 5] = numpy.sum(x349 * x457 * x485)
    result[1, 1, 6] = numpy.sum(x462 * x497)
    result[1, 1, 7] = numpy.sum(x231 * x492 * x498)
    result[1, 1, 8] = numpy.sum(x202 * x488 * x498)
    result[1, 1, 9] = numpy.sum(x357 * x462 * x485)
    result[1, 1, 10] = numpy.sum(x507 * x508)
    result[1, 1, 11] = numpy.sum(x131 * x461 * x497)
    result[1, 1, 12] = numpy.sum(x349 * x461 * x492)
    result[1, 1, 13] = numpy.sum(x357 * x461 * x488)
    result[1, 1, 14] = numpy.sum(x206 * x486 * x509)
    result[1, 2, 0] = numpy.sum(x239 * x245 * x435 * x439)
    result[1, 2, 1] = numpy.sum(x245 * x444 * x490)
    result[1, 2, 2] = numpy.sum(x412 * x435 * x447)
    result[1, 2, 3] = numpy.sum(x245 * x454 * x493)
    result[1, 2, 4] = numpy.sum(x306 * x444 * x457)
    result[1, 2, 5] = numpy.sum(x409 * x435 * x457)
    result[1, 2, 6] = numpy.sum(x245 * x462 * x472 * x489)
    result[1, 2, 7] = numpy.sum(x306 * x454 * x462)
    result[1, 2, 8] = numpy.sum(x260 * x444 * x498)
    result[1, 2, 9] = numpy.sum(x259 * x462 * x483)
    result[1, 2, 10] = numpy.sum(x245 * x481 * x508)
    result[1, 2, 11] = numpy.sum(x412 * x461 * x472)
    result[1, 2, 12] = numpy.sum(x409 * x454 * x461)
    result[1, 2, 13] = numpy.sum(x259 * x444 * x509)
    result[1, 2, 14] = numpy.sum(x265 * x461 * x510)
    result[1, 3, 0] = numpy.sum(x515 * x516)
    result[1, 3, 1] = numpy.sum(x521 * x523)
    result[1, 3, 2] = numpy.sum(x131 * x514 * x523)
    result[1, 3, 3] = numpy.sum(x524 * x531)
    result[1, 3, 4] = numpy.sum(x231 * x521 * x532)
    result[1, 3, 5] = numpy.sum(x136 * x282 * x514)
    result[1, 3, 6] = numpy.sum(x180 * x538)
    result[1, 3, 7] = numpy.sum(x131 * x531 * x539)
    result[1, 3, 8] = numpy.sum(x202 * x521 * x539)
    result[1, 3, 9] = numpy.sum(x180 * x289 * x514)
    result[1, 3, 10] = numpy.sum(x398 * x460 * x546)
    result[1, 3, 11] = numpy.sum(x131 * x460 * x538)
    result[1, 3, 12] = numpy.sum(x282 * x460 * x530)
    result[1, 3, 13] = numpy.sum(x289 * x460 * x521)
    result[1, 3, 14] = numpy.sum(x206 * x515 * x547)
    result[1, 4, 0] = numpy.sum(x245 * x516 * x548)
    result[1, 4, 1] = numpy.sum(x119 * x254 * x494)
    result[1, 4, 2] = numpy.sum(x119 * x306 * x485)
    result[1, 4, 3] = numpy.sum(x254 * x492 * x532)
    result[1, 4, 4] = numpy.sum(x136 * x298 * x488)
    result[1, 4, 5] = numpy.sum(x136 * x307 * x485)
    result[1, 4, 6] = numpy.sum(x180 * x254 * x549)
    result[1, 4, 7] = numpy.sum(x180 * x298 * x492)
    result[1, 4, 8] = numpy.sum(x180 * x302 * x488)
    result[1, 4, 9] = numpy.sum(x180 * x485 * x550)
    result[1, 4, 10] = numpy.sum(x460 * x507 * x551)
    result[1, 4, 11] = numpy.sum(x306 * x460 * x496)
    result[1, 4, 12] = numpy.sum(x307 * x460 * x492)
    result[1, 4, 13] = numpy.sum(x303 * x488 * x547)
    result[1, 4, 14] = numpy.sum(x308 * x485 * x547)
    result[1, 5, 0] = numpy.sum(x418 * x435 * x49)
    result[1, 5, 1] = numpy.sum(x316 * x444 * x522)
    result[1, 5, 2] = numpy.sum(x325 * x435 * x522)
    result[1, 5, 3] = numpy.sum(x316 * x454 * x524)
    result[1, 5, 4] = numpy.sum(x136 * x416 * x444)
    result[1, 5, 5] = numpy.sum(x417 * x435 * x524)
    result[1, 5, 6] = numpy.sum(x180 * x316 * x552)
    result[1, 5, 7] = numpy.sum(x180 * x416 * x454)
    result[1, 5, 8] = numpy.sum(x417 * x444 * x539)
    result[1, 5, 9] = numpy.sum(x180 * x329 * x483)
    result[1, 5, 10] = numpy.sum(x418 * x460 * x481)
    result[1, 5, 11] = numpy.sum(x325 * x460 * x552)
    result[1, 5, 12] = numpy.sum(x454 * x460 * x553)
    result[1, 5, 13] = numpy.sum(x329 * x444 * x547)
    result[1, 5, 14] = numpy.sum(x333 * x460 * x554)
    result[1, 6, 0] = numpy.sum(x47 * x559)
    result[1, 6, 1] = numpy.sum(x565 * x567)
    result[1, 6, 2] = numpy.sum(x131 * x558 * x567)
    result[1, 6, 3] = numpy.sum(x568 * x574)
    result[1, 6, 4] = numpy.sum(x231 * x565 * x575)
    result[1, 6, 5] = numpy.sum(x114 * x349 * x558)
    result[1, 6, 6] = numpy.sum(x215 * x581 * x582)
    result[1, 6, 7] = numpy.sum(x131 * x574 * x583)
    result[1, 6, 8] = numpy.sum(x202 * x565 * x583)
    result[1, 6, 9] = numpy.sum(x178 * x357 * x558)
    result[1, 6, 10] = numpy.sum(x24 * x589)
    result[1, 6, 11] = numpy.sum(48.2718229290016 * x131 * x590 * x591)
    result[1, 6, 12] = numpy.sum(x164 * x24 * x595)
    result[1, 6, 13] = numpy.sum(x189 * x24 * x565 * x596)
    result[1, 6, 14] = numpy.sum(x206 * x24 * x598)
    result[1, 7, 0] = numpy.sum(x47 * x514 * x551)
    result[1, 7, 1] = numpy.sum(x254 * x521 * x599)
    result[1, 7, 2] = numpy.sum(x306 * x33 * x514)
    result[1, 7, 3] = numpy.sum(x245 * x531 * x600)
    result[1, 7, 4] = numpy.sum(x114 * x298 * x521)
    result[1, 7, 5] = numpy.sum(x114 * x307 * x514)
    result[1, 7, 6] = numpy.sum(x254 * x537 * x583)
    result[1, 7, 7] = numpy.sum(x178 * x298 * x530)
    result[1, 7, 8] = numpy.sum(x178 * x302 * x521)
    result[1, 7, 9] = numpy.sum(x178 * x514 * x550)
    result[1, 7, 10] = numpy.sum(x296 * x546 * x590)
    result[1, 7, 11] = numpy.sum(x537 * x601 * x602)
    result[1, 7, 12] = numpy.sum(x252 * x603 * x604)
    result[1, 7, 13] = numpy.sum(x24 * x303 * x605)
    result[1, 7, 14] = numpy.sum(x308 * x514 * x601)
    result[1, 8, 0] = numpy.sum(x370 * x47 * x485)
    result[1, 8, 1] = numpy.sum(x316 * x33 * x494)
    result[1, 8, 2] = numpy.sum(x33 * x371 * x485)
    result[1, 8, 3] = numpy.sum(x316 * x492 * x600)
    result[1, 8, 4] = numpy.sum(x114 * x368 * x488)
    result[1, 8, 5] = numpy.sum(x417 * x485 * x600)
    result[1, 8, 6] = numpy.sum(x316 * x496 * x583)
    result[1, 8, 7] = numpy.sum(x178 * x368 * x492)
    result[1, 8, 8] = numpy.sum(x369 * x488 * x606)
    result[1, 8, 9] = numpy.sum(x328 * x485 * x607)
    result[1, 8, 10] = numpy.sum(x24 * x367 * x608)
    result[1, 8, 11] = numpy.sum(x315 * x549 * x601)
    result[1, 8, 12] = numpy.sum(x322 * x604 * x609)
    result[1, 8, 13] = numpy.sum(x328 * x494 * x601)
    result[1, 8, 14] = numpy.sum(x24 * x548 * x610)
    result[1, 9, 0] = numpy.sum(x373 * x47 * x483)
    result[1, 9, 1] = numpy.sum(x419 * x444 * x566)
    result[1, 9, 2] = numpy.sum(x422 * x435 * x566)
    result[1, 9, 3] = numpy.sum(x114 * x421 * x454)
    result[1, 9, 4] = numpy.sum(x422 * x444 * x575)
    result[1, 9, 5] = numpy.sum(x114 * x390 * x483)
    result[1, 9, 6] = numpy.sum(x419 * x472 * x582)
    result[1, 9, 7] = numpy.sum(x422 * x454 * x583)
    result[1, 9, 8] = numpy.sum(x382 * x444 * x607)
    result[1, 9, 9] = numpy.sum(x387 * x483 * x582)
    result[1, 9, 10] = numpy.sum(x373 * x481 * x601)
    result[1, 9, 11] = numpy.sum(x24 * x472 * x611)
    result[1, 9, 12] = numpy.sum(x390 * x454 * x601)
    result[1, 9, 13] = numpy.sum(x24 * x444 * x612)
    result[1, 9, 14] = numpy.sum(x24 * x391 * x435 * x597)
    result[1, 10, 0] = numpy.sum(x613 * x616)
    result[1, 10, 1] = numpy.sum(x617 * x619)
    result[1, 10, 2] = numpy.sum(x131 * x616 * x617)
    result[1, 10, 3] = numpy.sum(x215 * x620 * x621)
    result[1, 10, 4] = numpy.sum(x160 * x27 * x619)
    result[1, 10, 5] = numpy.sum(x27 * x400 * x615)
    result[1, 10, 6] = numpy.sum(x3 * x622)
    result[1, 10, 7] = numpy.sum(x103 * x160 * x620 * x623)
    result[1, 10, 8] = numpy.sum(x164 * x618 * x625)
    result[1, 10, 9] = numpy.sum(x3 * x615 * x626)
    result[1, 10, 10] = numpy.sum(
        x104
        * x587
        * (
            x0
            * (
                x274 * (x584 + x585)
                + x29 * (x292 + x543 + 4.0 * x576 + 4.0 * x577)
                + 3.0 * x544
                + 3.0 * x545
                + 4.0 * x579
                + 4.0 * x580
            )
            + x208 * x586
        )
    )
    result[1, 10, 11] = numpy.sum(x131 * x622)
    result[1, 10, 12] = numpy.sum(x164 * x620 * x627)
    result[1, 10, 13] = numpy.sum(x618 * x626)
    result[1, 10, 14] = numpy.sum(x207 * x593 * x615)
    result[1, 11, 0] = numpy.sum(x245 * x45 * x559)
    result[1, 11, 1] = numpy.sum(x254 * x565 * x628)
    result[1, 11, 2] = numpy.sum(x31 * x412 * x558)
    result[1, 11, 3] = numpy.sum(x27 * x408 * x573)
    result[1, 11, 4] = numpy.sum(x27 * x306 * x565)
    result[1, 11, 5] = numpy.sum(x27 * x409 * x558)
    result[1, 11, 6] = numpy.sum(x407 * x591 * x623)
    result[1, 11, 7] = numpy.sum(x573 * x602 * x624)
    result[1, 11, 8] = numpy.sum(x252 * x629 * x630)
    result[1, 11, 9] = numpy.sum(x259 * x558 * x624)
    result[1, 11, 10] = numpy.sum(x245 * x589)
    result[1, 11, 11] = numpy.sum(x248 * x581 * x596)
    result[1, 11, 12] = numpy.sum(x252 * x595)
    result[1, 11, 13] = numpy.sum(x259 * x629)
    result[1, 11, 14] = numpy.sum(x265 * x598)
    result[1, 12, 0] = numpy.sum(x418 * x45 * x514)
    result[1, 12, 1] = numpy.sum(x316 * x521 * x631)
    result[1, 12, 2] = numpy.sum(x325 * x514 * x631)
    result[1, 12, 3] = numpy.sum(x316 * x530 * x632)
    result[1, 12, 4] = numpy.sum(x27 * x416 * x521)
    result[1, 12, 5] = numpy.sum(x27 * x514 * x553)
    result[1, 12, 6] = numpy.sum(x3 * x311 * x633)
    result[1, 12, 7] = numpy.sum(x315 * x603 * x634)
    result[1, 12, 8] = numpy.sum(x322 * x605 * x634)
    result[1, 12, 9] = numpy.sum(x329 * x514 * x624)
    result[1, 12, 10] = numpy.sum(x311 * x546 * x627)
    result[1, 12, 11] = numpy.sum(x315 * x633)
    result[1, 12, 12] = numpy.sum(x323 * x603)
    result[1, 12, 13] = numpy.sum(x329 * x605)
    result[1, 12, 14] = numpy.sum(x515 * x610)
    result[1, 13, 0] = numpy.sum(x373 * x45 * x635)
    result[1, 13, 1] = numpy.sum(x419 * x488 * x628)
    result[1, 13, 2] = numpy.sum(x422 * x485 * x628)
    result[1, 13, 3] = numpy.sum(x27 * x421 * x492)
    result[1, 13, 4] = numpy.sum(x27 * x422 * x494)
    result[1, 13, 5] = numpy.sum(x27 * x390 * x635)
    result[1, 13, 6] = numpy.sum(x3 * x372 * x496 * x596)
    result[1, 13, 7] = numpy.sum(x3 * x379 * x609)
    result[1, 13, 8] = numpy.sum(x382 * x494 * x624)
    result[1, 13, 9] = numpy.sum(x3 * x485 * x612)
    result[1, 13, 10] = numpy.sum(x373 * x608)
    result[1, 13, 11] = numpy.sum(x496 * x611)
    result[1, 13, 12] = numpy.sum(x390 * x609)
    result[1, 13, 13] = numpy.sum(x488 * x612)
    result[1, 13, 14] = numpy.sum(x391 * x486 * x593)
    result[1, 14, 0] = numpy.sum(x427 * x435 * x613)
    result[1, 14, 1] = numpy.sum(x427 * x444 * x617)
    result[1, 14, 2] = numpy.sum(x31 * x426 * x510)
    result[1, 14, 3] = numpy.sum(x427 * x454 * x621)
    result[1, 14, 4] = numpy.sum(x426 * x459 * x636)
    result[1, 14, 5] = numpy.sum(x27 * x429 * x554)
    result[1, 14, 6] = numpy.sum(x3 * x423 * x637)
    result[1, 14, 7] = numpy.sum(x426 * x454 * x625)
    result[1, 14, 8] = numpy.sum(x429 * x459 * x624)
    result[1, 14, 9] = numpy.sum(x3 * x435 * x638)
    result[1, 14, 10] = numpy.sum(x423 * x481 * x639)
    result[1, 14, 11] = numpy.sum(x426 * x637)
    result[1, 14, 12] = numpy.sum(x429 * x454 * x627)
    result[1, 14, 13] = numpy.sum(x444 * x638)
    result[1, 14, 14] = numpy.sum(x432 * x435 * x639)
    result[2, 0, 0] = numpy.sum(x424 * x440 * x642)
    result[2, 0, 1] = numpy.sum(x109 * x642 * x643)
    result[2, 0, 2] = numpy.sum(x643 * x647)
    result[2, 0, 3] = numpy.sum(x157 * x458 * x648)
    result[2, 0, 4] = numpy.sum(x256 * x458 * x649)
    result[2, 0, 5] = numpy.sum(x167 * x458 * x654)
    result[2, 0, 6] = numpy.sum(x172 * x463 * x648)
    result[2, 0, 7] = numpy.sum(x378 * x463 * x649)
    result[2, 0, 8] = numpy.sum(x428 * x463 * x654)
    result[2, 0, 9] = numpy.sum(x463 * x664)
    result[2, 0, 10] = numpy.sum(x196 * x482 * x648)
    result[2, 0, 11] = numpy.sum(x172 * x647 * x665)
    result[2, 0, 12] = numpy.sum(x157 * x654 * x665)
    result[2, 0, 13] = numpy.sum(x109 * x482 * x664)
    result[2, 0, 14] = numpy.sum(x424 * x482 * x673)
    result[2, 1, 0] = numpy.sum(x208 * x642 * x674)
    result[2, 1, 1] = numpy.sum(x214 * x447 * x648)
    result[2, 1, 2] = numpy.sum(x208 * x647 * x675)
    result[2, 1, 3] = numpy.sum(x220 * x457 * x648)
    result[2, 1, 4] = numpy.sum(x222 * x457 * x676)
    result[2, 1, 5] = numpy.sum(x208 * x654 * x677)
    result[2, 1, 6] = numpy.sum(x228 * x462 * x648)
    result[2, 1, 7] = numpy.sum(x219 * x498 * x676)
    result[2, 1, 8] = numpy.sum(x222 * x654 * x678)
    result[2, 1, 9] = numpy.sum(x208 * x663 * x679)
    result[2, 1, 10] = numpy.sum(x238 * x642 * x680)
    result[2, 1, 11] = numpy.sum(x228 * x509 * x647)
    result[2, 1, 12] = numpy.sum(x220 * x509 * x654)
    result[2, 1, 13] = numpy.sum(x214 * x509 * x663)
    result[2, 1, 14] = numpy.sum(x208 * x673 * x681)
    result[2, 2, 0] = numpy.sum(x674 * x683)
    result[2, 2, 1] = numpy.sum(x109 * x675 * x683)
    result[2, 2, 2] = numpy.sum(x675 * x685)
    result[2, 2, 3] = numpy.sum(x457 * x683 * x686)
    result[2, 2, 4] = numpy.sum(x256 * x457 * x687)
    result[2, 2, 5] = numpy.sum(x677 * x689)
    result[2, 2, 6] = numpy.sum(x253 * x678 * x683)
    result[2, 2, 7] = numpy.sum(x378 * x498 * x685)
    result[2, 2, 8] = numpy.sum(x256 * x498 * x689)
    result[2, 2, 9] = numpy.sum(x679 * x691)
    result[2, 2, 10] = numpy.sum(x195 * x680 * x683)
    result[2, 2, 11] = numpy.sum(x253 * x509 * x685)
    result[2, 2, 12] = numpy.sum(x461 * x686 * x689)
    result[2, 2, 13] = numpy.sum(x109 * x420 * x461 * x691)
    result[2, 2, 14] = numpy.sum(x681 * x700)
    result[2, 3, 0] = numpy.sum(x268 * x49 * x648)
    result[2, 3, 1] = numpy.sum(x271 * x522 * x648)
    result[2, 3, 2] = numpy.sum(x267 * x522 * x676)
    result[2, 3, 3] = numpy.sum(x136 * x279 * x648)
    result[2, 3, 4] = numpy.sum(x136 * x281 * x676)
    result[2, 3, 5] = numpy.sum(x267 * x524 * x701)
    result[2, 3, 6] = numpy.sum(x180 * x286 * x648)
    result[2, 3, 7] = numpy.sum(x278 * x539 * x676)
    result[2, 3, 8] = numpy.sum(x180 * x281 * x701)
    result[2, 3, 9] = numpy.sum(x267 * x702 * x703)
    result[2, 3, 10] = numpy.sum(x293 * x547 * x642)
    result[2, 3, 11] = numpy.sum(x286 * x547 * x647)
    result[2, 3, 12] = numpy.sum(x279 * x547 * x654)
    result[2, 3, 13] = numpy.sum(x271 * x547 * x702)
    result[2, 3, 14] = numpy.sum(x268 * x547 * x673)
    result[2, 4, 0] = numpy.sum(x49 * x683 * x704)
    result[2, 4, 1] = numpy.sum(x119 * x222 * x705)
    result[2, 4, 2] = numpy.sum(x119 * x309 * x687)
    result[2, 4, 3] = numpy.sum(x136 * x297 * x705)
    result[2, 4, 4] = numpy.sum(x213 * x706 * x707)
    result[2, 4, 5] = numpy.sum(x309 * x532 * x689)
    result[2, 4, 6] = numpy.sum(x301 * x683 * x703)
    result[2, 4, 7] = numpy.sum(x219 * x703 * x706)
    result[2, 4, 8] = numpy.sum(x213 * x703 * x708)
    result[2, 4, 9] = numpy.sum(x180 * x309 * x709)
    result[2, 4, 10] = numpy.sum(x305 * x547 * x683)
    result[2, 4, 11] = numpy.sum(x301 * x547 * x685)
    result[2, 4, 12] = numpy.sum(x297 * x547 * x689)
    result[2, 4, 13] = numpy.sum(x222 * x547 * x691)
    result[2, 4, 14] = numpy.sum(x460 * x700 * x704)
    result[2, 5, 0] = numpy.sum(x167 * x49 * x713)
    result[2, 5, 1] = numpy.sum(x109 * x713 * x714)
    result[2, 5, 2] = numpy.sum(x714 * x719)
    result[2, 5, 3] = numpy.sum(x317 * x707 * x713)
    result[2, 5, 4] = numpy.sum(x256 * x532 * x719)
    result[2, 5, 5] = numpy.sum(x524 * x726)
    result[2, 5, 6] = numpy.sum(x324 * x703 * x713)
    result[2, 5, 7] = numpy.sum(x378 * x539 * x719)
    result[2, 5, 8] = numpy.sum(x109 * x539 * x726)
    result[2, 5, 9] = numpy.sum(x180 * x732)
    result[2, 5, 10] = numpy.sum(x330 * x547 * x713)
    result[2, 5, 11] = numpy.sum(x324 * x547 * x719)
    result[2, 5, 12] = numpy.sum(x317 * x547 * x725)
    result[2, 5, 13] = numpy.sum(x109 * x460 * x732)
    result[2, 5, 14] = numpy.sum(x167 * x460 * x739)
    result[2, 6, 0] = numpy.sum(x335 * x47 * x648)
    result[2, 6, 1] = numpy.sum(x340 * x566 * x648)
    result[2, 6, 2] = numpy.sum(x33 * x342 * x676)
    result[2, 6, 3] = numpy.sum(x346 * x568 * x648)
    result[2, 6, 4] = numpy.sum(x340 * x575 * x676)
    result[2, 6, 5] = numpy.sum(x334 * x568 * x701)
    result[2, 6, 6] = numpy.sum(x354 * x606 * x642)
    result[2, 6, 7] = numpy.sum(x346 * x607 * x647)
    result[2, 6, 8] = numpy.sum(x340 * x607 * x654)
    result[2, 6, 9] = numpy.sum(x342 * x606 * x663)
    result[2, 6, 10] = numpy.sum(x24 * x642 * x740)
    result[2, 6, 11] = numpy.sum(x24 * x647 * x741)
    result[2, 6, 12] = numpy.sum(x24 * x654 * x742)
    result[2, 6, 13] = numpy.sum(x24 * x663 * x743)
    result[2, 6, 14] = numpy.sum(x24 * x673 * x744)
    result[2, 7, 0] = numpy.sum(x361 * x47 * x705)
    result[2, 7, 1] = numpy.sum(x271 * x599 * x705)
    result[2, 7, 2] = numpy.sum(x599 * x685 * x745)
    result[2, 7, 3] = numpy.sum(x278 * x600 * x705)
    result[2, 7, 4] = numpy.sum(x271 * x706 * x746)
    result[2, 7, 5] = numpy.sum(x600 * x689 * x745)
    result[2, 7, 6] = numpy.sum(x285 * x607 * x683)
    result[2, 7, 7] = numpy.sum(x278 * x606 * x706)
    result[2, 7, 8] = numpy.sum(x271 * x606 * x708)
    result[2, 7, 9] = numpy.sum(x267 * x607 * x691)
    result[2, 7, 10] = numpy.sum(40.7971365319473 * x292 * x601 * x683)
    result[2, 7, 11] = numpy.sum(x285 * x601 * x687)
    result[2, 7, 12] = numpy.sum(139.348749811665 * x278 * x601 * x689)
    result[2, 7, 13] = numpy.sum(x271 * x601 * x709)
    result[2, 7, 14] = numpy.sum(x361 * x601 * x700)
    result[2, 8, 0] = numpy.sum(x47 * x704 * x713)
    result[2, 8, 1] = numpy.sum(x222 * x33 * x747)
    result[2, 8, 2] = numpy.sum(x309 * x599 * x719)
    result[2, 8, 3] = numpy.sum(x297 * x713 * x746)
    result[2, 8, 4] = numpy.sum(x213 * x746 * x748)
    result[2, 8, 5] = numpy.sum(x208 * x600 * x726)
    result[2, 8, 6] = numpy.sum(x227 * x607 * x713)
    result[2, 8, 7] = numpy.sum(x219 * x606 * x748)
    result[2, 8, 8] = numpy.sum(241.359114645008 * x213 * x606 * x725)
    result[2, 8, 9] = numpy.sum(x309 * x583 * x731)
    result[2, 8, 10] = numpy.sum(x24 * x305 * x749)
    result[2, 8, 11] = numpy.sum(x24 * x301 * x750)
    result[2, 8, 12] = numpy.sum(x24 * x297 * x751)
    result[2, 8, 13] = numpy.sum(x222 * x601 * x731)
    result[2, 8, 14] = numpy.sum(40.7971365319473 * x208 * x739 * x753)
    result[2, 9, 0] = numpy.sum(x47 * x758)
    result[2, 9, 1] = numpy.sum(x109 * x757 * x759)
    result[2, 9, 2] = numpy.sum(x759 * x765)
    result[2, 9, 3] = numpy.sum(x114 * x686 * x757)
    result[2, 9, 4] = numpy.sum(x256 * x575 * x765)
    result[2, 9, 5] = numpy.sum(x568 * x771)
    result[2, 9, 6] = numpy.sum(x253 * x606 * x757)
    result[2, 9, 7] = numpy.sum(x378 * x583 * x765)
    result[2, 9, 8] = numpy.sum(x109 * x583 * x771)
    result[2, 9, 9] = numpy.sum(x223 * x582 * x777)
    result[2, 9, 10] = numpy.sum(x195 * x24 * x778)
    result[2, 9, 11] = numpy.sum(x24 * x253 * x779)
    result[2, 9, 12] = numpy.sum(x156 * x24 * x594 * x770)
    result[2, 9, 13] = numpy.sum(48.2718229290016 * x109 * x753 * x777)
    result[2, 9, 14] = numpy.sum(x24 * x784)
    result[2, 10, 0] = numpy.sum(x393 * x613 * x648)
    result[2, 10, 1] = numpy.sum(x395 * x617 * x648)
    result[2, 10, 2] = numpy.sum(x393 * x617 * x676)
    result[2, 10, 3] = numpy.sum(x397 * x621 * x648)
    result[2, 10, 4] = numpy.sum(x395 * x636 * x649)
    result[2, 10, 5] = numpy.sum(x393 * x621 * x701)
    result[2, 10, 6] = numpy.sum(x3 * x642 * x785)
    result[2, 10, 7] = numpy.sum(x397 * x624 * x649)
    result[2, 10, 8] = numpy.sum(x395 * x625 * x654)
    result[2, 10, 9] = numpy.sum(x3 * x393 * x786)
    result[2, 10, 10] = numpy.sum(x405 * x639 * x642)
    result[2, 10, 11] = numpy.sum(x647 * x785)
    result[2, 10, 12] = numpy.sum(x397 * x627 * x654)
    result[2, 10, 13] = numpy.sum(x395 * x786)
    result[2, 10, 14] = numpy.sum(x393 * x639 * x673)
    result[2, 11, 0] = numpy.sum(x335 * x45 * x705)
    result[2, 11, 1] = numpy.sum(x340 * x628 * x705)
    result[2, 11, 2] = numpy.sum(x342 * x685 * x787)
    result[2, 11, 3] = numpy.sum(x346 * x683 * x788)
    result[2, 11, 4] = numpy.sum(x348 * x636 * x685)
    result[2, 11, 5] = numpy.sum(x334 * x689 * x788)
    result[2, 11, 6] = numpy.sum(x3 * x683 * x741)
    result[2, 11, 7] = numpy.sum(x346 * x624 * x687)
    result[2, 11, 8] = numpy.sum(x348 * x624 * x689)
    result[2, 11, 9] = numpy.sum(x342 * x624 * x691)
    result[2, 11, 10] = numpy.sum(x683 * x740)
    result[2, 11, 11] = numpy.sum(x685 * x741)
    result[2, 11, 12] = numpy.sum(x689 * x742)
    result[2, 11, 13] = numpy.sum(x691 * x743)
    result[2, 11, 14] = numpy.sum(x700 * x744)
    result[2, 12, 0] = numpy.sum(x268 * x45 * x747)
    result[2, 12, 1] = numpy.sum(x271 * x631 * x747)
    result[2, 12, 2] = numpy.sum(x631 * x719 * x745)
    result[2, 12, 3] = numpy.sum(x279 * x636 * x713)
    result[2, 12, 4] = numpy.sum(x281 * x636 * x719)
    result[2, 12, 5] = numpy.sum(x632 * x725 * x745)
    result[2, 12, 6] = numpy.sum(x286 * x3 * x749)
    result[2, 12, 7] = numpy.sum(x278 * x634 * x750)
    result[2, 12, 8] = numpy.sum(x281 * x3 * x751)
    result[2, 12, 9] = numpy.sum(x267 * x3 * x789)
    result[2, 12, 10] = numpy.sum(x293 * x749)
    result[2, 12, 11] = numpy.sum(x286 * x750)
    result[2, 12, 12] = numpy.sum(x279 * x751)
    result[2, 12, 13] = numpy.sum(x271 * x789)
    result[2, 12, 14] = numpy.sum(x268 * x593 * x739)
    result[2, 13, 0] = numpy.sum(x208 * x45 * x758)
    result[2, 13, 1] = numpy.sum(x214 * x757 * x787)
    result[2, 13, 2] = numpy.sum(x309 * x628 * x765)
    result[2, 13, 3] = numpy.sum(x220 * x636 * x757)
    result[2, 13, 4] = numpy.sum(x222 * x636 * x765)
    result[2, 13, 5] = numpy.sum(x208 * x224 * x27 * x770)
    result[2, 13, 6] = numpy.sum(x228 * x624 * x757)
    result[2, 13, 7] = numpy.sum(x219 * x630 * x779)
    result[2, 13, 8] = numpy.sum(x222 * x3 * x790)
    result[2, 13, 9] = numpy.sum(x233 * x3 * x752 * x777)
    result[2, 13, 10] = numpy.sum(x238 * x778)
    result[2, 13, 11] = numpy.sum(x228 * x779)
    result[2, 13, 12] = numpy.sum(x220 * x790)
    result[2, 13, 13] = numpy.sum(x214 * x593 * x777)
    result[2, 13, 14] = numpy.sum(x208 * x784)
    result[2, 14, 0] = numpy.sum(x613 * x793)
    result[2, 14, 1] = numpy.sum(x109 * x617 * x793)
    result[2, 14, 2] = numpy.sum(x617 * x795)
    result[2, 14, 3] = numpy.sum(x157 * x636 * x792)
    result[2, 14, 4] = numpy.sum(x185 * x27 * x795)
    result[2, 14, 5] = numpy.sum(x223 * x621 * x796)
    result[2, 14, 6] = numpy.sum(x172 * x3 * x797)
    result[2, 14, 7] = numpy.sum(x156 * x798 * x799)
    result[2, 14, 8] = numpy.sum(x110 * x799 * x800)
    result[2, 14, 9] = numpy.sum(x3 * x801)
    result[2, 14, 10] = numpy.sum(x196 * x797)
    result[2, 14, 11] = numpy.sum(x172 * x798)
    result[2, 14, 12] = numpy.sum(x159 * x800)
    result[2, 14, 13] = numpy.sum(x109 * x801)
    result[2, 14, 14] = numpy.sum(
        6.89597470414309
        * x752
        * (
            x0
            * (
                x29 * (x333 + x736 + 4.0 * x772 + 4.0 * x773)
                + x318 * (x780 + x781)
                + 3.0 * x737
                + 3.0 * x738
                + 4.0 * x775
                + 4.0 * x776
            )
            + x245 * x782
        )
    )
    return result


diag_quadrupole3d = {
    (0, 0): diag_quadrupole3d_00,
    (0, 1): diag_quadrupole3d_01,
    (0, 2): diag_quadrupole3d_02,
    (0, 3): diag_quadrupole3d_03,
    (0, 4): diag_quadrupole3d_04,
    (1, 0): diag_quadrupole3d_10,
    (1, 1): diag_quadrupole3d_11,
    (1, 2): diag_quadrupole3d_12,
    (1, 3): diag_quadrupole3d_13,
    (1, 4): diag_quadrupole3d_14,
    (2, 0): diag_quadrupole3d_20,
    (2, 1): diag_quadrupole3d_21,
    (2, 2): diag_quadrupole3d_22,
    (2, 3): diag_quadrupole3d_23,
    (2, 4): diag_quadrupole3d_24,
    (3, 0): diag_quadrupole3d_30,
    (3, 1): diag_quadrupole3d_31,
    (3, 2): diag_quadrupole3d_32,
    (3, 3): diag_quadrupole3d_33,
    (3, 4): diag_quadrupole3d_34,
    (4, 0): diag_quadrupole3d_40,
    (4, 1): diag_quadrupole3d_41,
    (4, 2): diag_quadrupole3d_42,
    (4, 3): diag_quadrupole3d_43,
    (4, 4): diag_quadrupole3d_44,
}
