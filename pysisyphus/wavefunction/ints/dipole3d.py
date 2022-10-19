import numpy

"""

        Dipole integrals are given in the order:
        for bf_a in basis_functions_a:
            for bf_b in basis_functions_b:
                for cart_dir in (x, y, z):
                    dipole_integrals(bf_a, bf_b, cart_dir)

        So for <s_a|μ|s_b> it will be:

            <s_a|x|s_b>
            <s_a|y|s_b>
            <s_a|z|s_b>

"""


def dipole3d_00(a, A, b, B, C):
    """Cartesian 3D (ss) dipole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = a * b * x0
    x2 = (
        5.56832799683171
        * x0**1.5
        * numpy.exp(-x1 * (A[0] - B[0]) ** 2)
        * numpy.exp(-x1 * (A[1] - B[1]) ** 2)
        * numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    )

    # 3 item(s)
    return numpy.array(
        [
            x2 * (x0 * (a * A[0] + b * B[0]) - C[0]),
            x2 * (x0 * (a * A[1] + b * B[1]) - C[1]),
            x2 * (x0 * (a * A[2] + b * B[2]) - C[2]),
        ]
    )


def dipole3d_01(a, A, b, B, C):
    """Cartesian 3D (sp) dipole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = a * b * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(x0)
    x5 = 1.77245385090552 * x4
    x6 = x5 / (2.0 * a + 2.0 * b)
    x7 = -x0 * (a * A[0] + b * B[0])
    x8 = -x7 - B[0]
    x9 = x3 * (-x7 - C[0])
    x10 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x11 = 3.14159265358979 * x0 * x10
    x12 = 5.56832799683171 * x4
    x13 = x12 * x9
    x14 = -x0 * (a * A[1] + b * B[1])
    x15 = x2 * (-x14 - B[1])
    x16 = x0 * x10 * x15
    x17 = -x0 * (a * A[2] + b * B[2])
    x18 = x10 * (-x17 - B[2])
    x19 = x0 * x2
    x20 = -x14 - C[1]
    x21 = x19 * x3
    x22 = x12 * x20 * x21
    x23 = x10 * x8
    x24 = -x17 - C[2]
    x25 = x12 * x24

    # 9 item(s)
    return numpy.array(
        [
            x11 * x2 * (x3 * x6 + x5 * x8 * x9),
            x13 * x16,
            x13 * x18 * x19,
            x22 * x23,
            x11 * x3 * (x15 * x20 * x5 + x2 * x6),
            x18 * x22,
            x21 * x23 * x25,
            x16 * x25 * x3,
            3.14159265358979 * x21 * (x10 * x6 + x18 * x24 * x5),
        ]
    )


def dipole3d_02(a, A, b, B, C):
    """Cartesian 3D (sd) dipole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = numpy.sqrt(x1)
    x3 = 1.77245385090552 * x2
    x4 = -x1 * (a * A[0] + b * B[0])
    x5 = -x4 - B[0]
    x6 = a * b * x1
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x4 - C[0]
    x11 = x3 * x7
    x12 = x0 * x11
    x13 = x10 * x9 + x12
    x14 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x15 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x16 = 3.14159265358979 * x1 * x15
    x17 = x14 * x16
    x18 = -x1 * (a * A[1] + b * B[1])
    x19 = -x18 - B[1]
    x20 = x13 * x17
    x21 = -x1 * (a * A[2] + b * B[2])
    x22 = -x21 - B[2]
    x23 = x14 * x3
    x24 = x0 * x23
    x25 = x16 * x7
    x26 = x25 * (x19**2 * x23 + x24)
    x27 = 5.56832799683171
    x28 = x1 * x14
    x29 = x15 * x2 * x22 * x27 * x28
    x30 = x15 * x3
    x31 = x0 * x30
    x32 = 3.14159265358979 * x28
    x33 = x32 * x7
    x34 = x33 * (x22**2 * x30 + x31)
    x35 = -x18 - C[1]
    x36 = x17 * (x11 * x5**2 + x12)
    x37 = x19 * x23
    x38 = x24 + x35 * x37
    x39 = -x21 - C[2]
    x40 = x22 * x30
    x41 = x31 + x39 * x40

    # 18 item(s)
    return numpy.array(
        [
            x17 * (x0 * (x10 * x11 + x9) + x13 * x5),
            x19 * x20,
            x20 * x22,
            x10 * x26,
            x10 * x19 * x29 * x7,
            x10 * x34,
            x35 * x36,
            x16 * x38 * x8,
            x29 * x35 * x8,
            x25 * (x0 * (x23 * x35 + x37) + x19 * x38),
            x22 * x25 * x38,
            x34 * x35,
            x36 * x39,
            x15 * x19 * x2 * x27 * x28 * x39 * x8,
            x32 * x41 * x8,
            x26 * x39,
            x19 * x33 * x41,
            x33 * (x0 * (x30 * x39 + x40) + x22 * x41),
        ]
    )


def dipole3d_03(a, A, b, B, C):
    """Cartesian 3D (sf) dipole moment integrals.
    The origin is at C.

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
    x10 = -x2 - C[0]
    x11 = x3 * x7
    x12 = x10 * x11
    x13 = x12 + x9
    x14 = x0 * (x10 * x7 + x11) + x13 * x3
    x15 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x16 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x17 = 3.14159265358979 * x1 * x16
    x18 = x15 * x17
    x19 = -x1 * (a * A[1] + b * B[1])
    x20 = -x19 - B[1]
    x21 = x14 * x18
    x22 = -x1 * (a * A[2] + b * B[2])
    x23 = -x22 - B[2]
    x24 = x15 * x6
    x25 = x20**2 * x24
    x26 = x0 * x24
    x27 = x25 + x26
    x28 = x16 * x6
    x29 = x18 * x23
    x30 = x23**2 * x28
    x31 = x0 * x28
    x32 = x30 + x31
    x33 = x17 * x5
    x34 = x33 * (2.0 * x20 * x26 + x20 * x27)
    x35 = x23 * x33
    x36 = 3.14159265358979 * x1 * x15 * x5
    x37 = x32 * x36
    x38 = x36 * (2.0 * x23 * x31 + x23 * x32)
    x39 = -x19 - C[1]
    x40 = x8 + x9
    x41 = x18 * (2.0 * x0 * x11 + x3 * x40)
    x42 = x20 * x24
    x43 = x39 * x42
    x44 = x26 + x43
    x45 = x0 * (x24 * x39 + x42) + x20 * x44
    x46 = x33 * x45
    x47 = -x22 - C[2]
    x48 = x23 * x28
    x49 = x47 * x48
    x50 = x31 + x49
    x51 = x0 * (x28 * x47 + x48) + x23 * x50
    x52 = x36 * x51

    # 30 item(s)
    return numpy.array(
        [
            x18 * (x0 * (2.0 * x12 + x8 + 3.0 * x9) + x14 * x3),
            x20 * x21,
            x21 * x23,
            x13 * x27 * x28,
            x13 * x20 * x29,
            x13 * x24 * x32,
            x10 * x34,
            x10 * x27 * x35,
            x10 * x20 * x37,
            x10 * x38,
            x39 * x41,
            x28 * x40 * x44,
            x29 * x39 * x40,
            x3 * x46,
            x3 * x35 * x44,
            x3 * x37 * x39,
            x33 * (x0 * (x25 + 3.0 * x26 + 2.0 * x43) + x20 * x45),
            x23 * x46,
            x32 * x44 * x7,
            x38 * x39,
            x41 * x47,
            x18 * x20 * x40 * x47,
            x24 * x40 * x50,
            x27 * x3 * x33 * x47,
            x20 * x3 * x36 * x50,
            x3 * x52,
            x34 * x47,
            x27 * x50 * x7,
            x20 * x52,
            x36 * (x0 * (x30 + 3.0 * x31 + 2.0 * x49) + x23 * x51),
        ]
    )


def dipole3d_04(a, A, b, B, C):
    """Cartesian 3D (sg) dipole moment integrals.
    The origin is at C.

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
    x16 = x3**2 * x7
    x17 = x12 + x16
    x18 = 2.0 * x12 * x3 + x17 * x3
    x19 = 3.0 * x12
    x20 = x11 + x15
    x21 = x0 * (2.0 * x13 + x16 + x19) + x20 * x3
    x22 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x23 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x24 = 3.14159265358979 * x1 * x23
    x25 = x22 * x24
    x26 = -x1 * (a * A[1] + b * B[1])
    x27 = -x26 - B[1]
    x28 = x21 * x25
    x29 = -x1 * (a * A[2] + b * B[2])
    x30 = -x29 - B[2]
    x31 = x22 * x6
    x32 = x27**2 * x31
    x33 = x0 * x31
    x34 = x32 + x33
    x35 = x23 * x6
    x36 = x25 * x30
    x37 = x30**2 * x35
    x38 = x0 * x35
    x39 = x37 + x38
    x40 = 2.0 * x27 * x33 + x27 * x34
    x41 = x30 * x35
    x42 = x27 * x31
    x43 = 2.0 * x30 * x38 + x30 * x39
    x44 = 3.0 * x33
    x45 = x24 * x5
    x46 = x45 * (x0 * (3.0 * x32 + x44) + x27 * x40)
    x47 = x30 * x45
    x48 = 3.14159265358979 * x1 * x22 * x5
    x49 = x43 * x48
    x50 = 3.0 * x38
    x51 = x48 * (x0 * (3.0 * x37 + x50) + x30 * x43)
    x52 = -x26 - C[1]
    x53 = x25 * (x0 * (3.0 * x16 + x19) + x18 * x3)
    x54 = x42 * x52
    x55 = x33 + x54
    x56 = x31 * x52
    x57 = x0 * (x42 + x56)
    x58 = x27 * x55
    x59 = x57 + x58
    x60 = x0 * (x32 + x44 + 2.0 * x54) + x27 * x59
    x61 = x45 * x60
    x62 = -x29 - C[2]
    x63 = x41 * x62
    x64 = x38 + x63
    x65 = x35 * x62
    x66 = x0 * (x41 + x65)
    x67 = x30 * x64
    x68 = x66 + x67
    x69 = x0 * (x37 + x50 + 2.0 * x63) + x30 * x68
    x70 = x48 * x69

    # 45 item(s)
    return numpy.array(
        [
            x25 * (x0 * (3.0 * x11 + 3.0 * x15 + x18) + x21 * x3),
            x27 * x28,
            x28 * x30,
            x20 * x34 * x35,
            x20 * x27 * x36,
            x20 * x31 * x39,
            x14 * x35 * x40,
            x14 * x34 * x41,
            x14 * x39 * x42,
            x14 * x31 * x43,
            x46 * x9,
            x40 * x47 * x9,
            x10 * x34 * x39,
            x27 * x49 * x9,
            x51 * x9,
            x52 * x53,
            x18 * x35 * x55,
            x18 * x36 * x52,
            x17 * x35 * x59,
            x17 * x41 * x55,
            x17 * x39 * x56,
            x3 * x61,
            x3 * x47 * x59,
            x39 * x55 * x8,
            x3 * x49 * x52,
            x45 * (x0 * (x40 + 3.0 * x57 + 3.0 * x58) + x27 * x60),
            x30 * x61,
            x39 * x59 * x7,
            x43 * x55 * x7,
            x51 * x52,
            x53 * x62,
            x18 * x25 * x27 * x62,
            x18 * x31 * x64,
            x17 * x34 * x65,
            x17 * x42 * x64,
            x17 * x31 * x68,
            x3 * x40 * x45 * x62,
            x34 * x64 * x8,
            x27 * x3 * x48 * x68,
            x3 * x70,
            x46 * x62,
            x40 * x64 * x7,
            x34 * x68 * x7,
            x27 * x70,
            x48 * (x0 * (x43 + 3.0 * x66 + 3.0 * x67) + x30 * x69),
        ]
    )


def dipole3d_10(a, A, b, B, C):
    """Cartesian 3D (ps) dipole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (a + b) ** (-1.0)
    x1 = a * b * x0
    x2 = numpy.exp(-x1 * (A[1] - B[1]) ** 2)
    x3 = numpy.exp(-x1 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(x0)
    x5 = 1.77245385090552 * x4
    x6 = x5 / (2.0 * a + 2.0 * b)
    x7 = -x0 * (a * A[0] + b * B[0])
    x8 = -x7 - A[0]
    x9 = x3 * (-x7 - C[0])
    x10 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x11 = 3.14159265358979 * x0 * x10
    x12 = 5.56832799683171 * x4
    x13 = x12 * x9
    x14 = -x0 * (a * A[1] + b * B[1])
    x15 = x2 * (-x14 - A[1])
    x16 = x0 * x10 * x15
    x17 = -x0 * (a * A[2] + b * B[2])
    x18 = x10 * (-x17 - A[2])
    x19 = x0 * x2
    x20 = -x14 - C[1]
    x21 = x19 * x3
    x22 = x12 * x20 * x21
    x23 = x10 * x8
    x24 = -x17 - C[2]
    x25 = x12 * x24

    # 9 item(s)
    return numpy.array(
        [
            x11 * x2 * (x3 * x6 + x5 * x8 * x9),
            x13 * x16,
            x13 * x18 * x19,
            x22 * x23,
            x11 * x3 * (x15 * x20 * x5 + x2 * x6),
            x18 * x22,
            x21 * x23 * x25,
            x16 * x25 * x3,
            3.14159265358979 * x21 * (x10 * x6 + x18 * x24 * x5),
        ]
    )


def dipole3d_11(a, A, b, B, C):
    """Cartesian 3D (pp) dipole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = numpy.sqrt(x1)
    x3 = 1.77245385090552 * x2
    x4 = -x1 * (a * A[0] + b * B[0])
    x5 = a * b * x1
    x6 = numpy.exp(-x5 * (A[0] - B[0]) ** 2)
    x7 = x6 * (-x4 - B[0])
    x8 = x3 * x7
    x9 = -x4 - C[0]
    x10 = x3 * x6
    x11 = x10 * x9
    x12 = -x4 - A[0]
    x13 = x0 * x10
    x14 = x13 + x8 * x9
    x15 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x16 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x17 = 3.14159265358979 * x1 * x16
    x18 = x15 * x17
    x19 = -x1 * (a * A[1] + b * B[1])
    x20 = -x19 - B[1]
    x21 = x18 * (x11 * x12 + x13)
    x22 = -x1 * (a * A[2] + b * B[2])
    x23 = -x22 - B[2]
    x24 = -x19 - A[1]
    x25 = x14 * x18
    x26 = x0 * x3
    x27 = x15 * x26
    x28 = x15 * x3
    x29 = x20 * x28
    x30 = x17 * x6
    x31 = x30 * (x24 * x29 + x27)
    x32 = x1 * x15
    x33 = 5.56832799683171 * x16 * x2 * x32
    x34 = x33 * x6
    x35 = x34 * x9
    x36 = -x22 - A[2]
    x37 = x16 * x26
    x38 = x16 * x3
    x39 = x23 * x38
    x40 = 3.14159265358979 * x32
    x41 = x40 * x6
    x42 = x41 * (x36 * x39 + x37)
    x43 = -x19 - C[1]
    x44 = x18 * (x12 * x8 + x13)
    x45 = x27 + x29 * x43
    x46 = x30 * x45
    x47 = x12 * x34
    x48 = x28 * x43
    x49 = x24 * x48 + x27
    x50 = x33 * x7
    x51 = -x22 - C[2]
    x52 = x37 + x39 * x51
    x53 = x41 * x52
    x54 = x38 * x51
    x55 = x36 * x54 + x37

    # 27 item(s)
    return numpy.array(
        [
            x18 * (x0 * (x11 + x8) + x12 * x14),
            x20 * x21,
            x21 * x23,
            x24 * x25,
            x31 * x9,
            x23 * x24 * x35,
            x25 * x36,
            x20 * x35 * x36,
            x42 * x9,
            x43 * x44,
            x12 * x46,
            x23 * x43 * x47,
            x17 * x49 * x7,
            x30 * (x0 * (x29 + x48) + x24 * x45),
            x23 * x30 * x49,
            x36 * x43 * x50,
            x36 * x46,
            x42 * x43,
            x44 * x51,
            x20 * x47 * x51,
            x12 * x53,
            x24 * x50 * x51,
            x31 * x51,
            x24 * x53,
            x40 * x55 * x7,
            x20 * x41 * x55,
            x41 * (x0 * (x39 + x54) + x36 * x52),
        ]
    )


def dipole3d_12(a, A, b, B, C):
    """Cartesian 3D (pd) dipole moment integrals.
    The origin is at C.

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
    x10 = -x2 - C[0]
    x11 = x3 * x7
    x12 = x10 * x11
    x13 = -x2 - A[0]
    x14 = x10 * x7
    x15 = x0 * (x11 + x14)
    x16 = x12 + x9
    x17 = x15 + x16 * x3
    x18 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x19 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x20 = 3.14159265358979 * x1 * x19
    x21 = x18 * x20
    x22 = -x1 * (a * A[1] + b * B[1])
    x23 = -x22 - B[1]
    x24 = x21 * (x13 * x16 + x15)
    x25 = -x1 * (a * A[2] + b * B[2])
    x26 = -x25 - B[2]
    x27 = x18 * x6
    x28 = x23**2 * x27
    x29 = x0 * x27
    x30 = x28 + x29
    x31 = x13 * x14 + x9
    x32 = x19 * x6
    x33 = x21 * x26
    x34 = x26**2 * x32
    x35 = x0 * x32
    x36 = x34 + x35
    x37 = -x22 - A[1]
    x38 = x17 * x21
    x39 = x23 * x27
    x40 = x29 + x37 * x39
    x41 = x20 * x5
    x42 = x41 * (2.0 * x23 * x29 + x30 * x37)
    x43 = x10 * x41
    x44 = 3.14159265358979 * x1 * x18 * x5
    x45 = x10 * x44
    x46 = -x25 - A[2]
    x47 = x21 * x46
    x48 = x26 * x32
    x49 = x35 + x46 * x48
    x50 = x44 * (2.0 * x26 * x35 + x36 * x46)
    x51 = -x22 - C[1]
    x52 = x8 + x9
    x53 = x21 * (2.0 * x0 * x11 + x13 * x52)
    x54 = x39 * x51
    x55 = x29 + x54
    x56 = x11 * x13 + x9
    x57 = x27 * x51
    x58 = x0 * (x39 + x57)
    x59 = x23 * x55 + x58
    x60 = x41 * x59
    x61 = x26 * x41
    x62 = x44 * x51
    x63 = x29 + x37 * x57
    x64 = x41 * (x37 * x55 + x58)
    x65 = x3 * x41
    x66 = -x25 - C[2]
    x67 = x21 * x66
    x68 = x48 * x66
    x69 = x35 + x68
    x70 = x44 * x69
    x71 = x32 * x66
    x72 = x0 * (x48 + x71)
    x73 = x26 * x69 + x72
    x74 = x44 * x73
    x75 = x35 + x46 * x71
    x76 = x44 * (x46 * x69 + x72)

    # 54 item(s)
    return numpy.array(
        [
            x21 * (x0 * (2.0 * x12 + x8 + 3.0 * x9) + x13 * x17),
            x23 * x24,
            x24 * x26,
            x30 * x31 * x32,
            x23 * x31 * x33,
            x27 * x31 * x36,
            x37 * x38,
            x16 * x32 * x40,
            x16 * x33 * x37,
            x10 * x42,
            x26 * x40 * x43,
            x36 * x37 * x45,
            x38 * x46,
            x16 * x23 * x47,
            x16 * x27 * x49,
            x30 * x43 * x46,
            x23 * x45 * x49,
            x10 * x50,
            x51 * x53,
            x32 * x55 * x56,
            x33 * x51 * x56,
            x13 * x60,
            x13 * x55 * x61,
            x13 * x36 * x62,
            x32 * x52 * x63,
            x3 * x64,
            x3 * x61 * x63,
            x41 * (x0 * (x28 + 3.0 * x29 + 2.0 * x54) + x37 * x59),
            x26 * x64,
            x36 * x63 * x7,
            x47 * x51 * x52,
            x46 * x55 * x65,
            x3 * x49 * x62,
            x46 * x60,
            x49 * x55 * x7,
            x50 * x51,
            x53 * x66,
            x23 * x56 * x67,
            x27 * x56 * x69,
            x13 * x30 * x41 * x66,
            x13 * x23 * x70,
            x13 * x74,
            x37 * x52 * x67,
            x40 * x65 * x66,
            x3 * x37 * x70,
            x42 * x66,
            x40 * x69 * x7,
            x37 * x74,
            x27 * x52 * x75,
            x23 * x3 * x44 * x75,
            x3 * x76,
            x30 * x7 * x75,
            x23 * x76,
            x44 * (x0 * (x34 + 3.0 * x35 + 2.0 * x68) + x46 * x73),
        ]
    )


def dipole3d_13(a, A, b, B, C):
    """Cartesian 3D (pf) dipole moment integrals.
    The origin is at C.

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
    x16 = x3**2 * x7
    x17 = x12 + x16
    x18 = 2.0 * x12 * x3
    x19 = x17 * x3 + x18
    x20 = -x2 - A[0]
    x21 = 3.0 * x12
    x22 = x0 * (2.0 * x13 + x16 + x21)
    x23 = x11 + x15
    x24 = x22 + x23 * x3
    x25 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x26 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x27 = 3.14159265358979 * x1 * x26
    x28 = x25 * x27
    x29 = -x1 * (a * A[1] + b * B[1])
    x30 = -x29 - B[1]
    x31 = x28 * (x20 * x23 + x22)
    x32 = -x1 * (a * A[2] + b * B[2])
    x33 = -x32 - B[2]
    x34 = x11 + x14 * x20
    x35 = x25 * x6
    x36 = x30**2 * x35
    x37 = x0 * x35
    x38 = x36 + x37
    x39 = x26 * x6
    x40 = x38 * x39
    x41 = x28 * x33
    x42 = x33**2 * x39
    x43 = x0 * x39
    x44 = x42 + x43
    x45 = x35 * x44
    x46 = 2.0 * x30 * x37
    x47 = x30 * x38 + x46
    x48 = x10 * x20 + x12
    x49 = x33 * x39
    x50 = x30 * x35
    x51 = 2.0 * x33 * x43
    x52 = x33 * x44 + x51
    x53 = -x29 - A[1]
    x54 = x24 * x28
    x55 = x37 + x50 * x53
    x56 = x38 * x53 + x46
    x57 = 3.0 * x37
    x58 = x27 * x5
    x59 = x58 * (x0 * (3.0 * x36 + x57) + x47 * x53)
    x60 = x58 * x9
    x61 = 3.14159265358979 * x1 * x25 * x5
    x62 = x61 * x9
    x63 = -x32 - A[2]
    x64 = x28 * x63
    x65 = x43 + x49 * x63
    x66 = x44 * x63 + x51
    x67 = 3.0 * x43
    x68 = x61 * (x0 * (3.0 * x42 + x67) + x52 * x63)
    x69 = -x29 - C[1]
    x70 = x28 * (x0 * (3.0 * x16 + x21) + x19 * x20)
    x71 = x17 * x20 + x18
    x72 = x50 * x69
    x73 = x37 + x72
    x74 = x39 * x73
    x75 = x35 * x69
    x76 = x0 * (x50 + x75)
    x77 = x30 * x73
    x78 = x76 + x77
    x79 = x12 + x20 * x8
    x80 = x0 * (x36 + x57 + 2.0 * x72)
    x81 = x30 * x78 + x80
    x82 = x58 * x81
    x83 = x33 * x58
    x84 = x44 * x7
    x85 = x61 * x69
    x86 = x37 + x53 * x75
    x87 = x53 * x73 + x76
    x88 = x58 * (x53 * x78 + x80)
    x89 = x3 * x58
    x90 = -x32 - C[2]
    x91 = x28 * x90
    x92 = x49 * x90
    x93 = x43 + x92
    x94 = x35 * x93
    x95 = x39 * x90
    x96 = x0 * (x49 + x95)
    x97 = x33 * x93
    x98 = x96 + x97
    x99 = x7 * x93
    x100 = x61 * x98
    x101 = x0 * (x42 + x67 + 2.0 * x92)
    x102 = x101 + x33 * x98
    x103 = x102 * x61
    x104 = x43 + x63 * x95
    x105 = x63 * x93 + x96
    x106 = x61 * (x101 + x63 * x98)

    # 90 item(s)
    return numpy.array(
        [
            x28 * (x0 * (3.0 * x11 + 3.0 * x15 + x19) + x20 * x24),
            x30 * x31,
            x31 * x33,
            x34 * x40,
            x30 * x34 * x41,
            x34 * x45,
            x39 * x47 * x48,
            x38 * x48 * x49,
            x44 * x48 * x50,
            x35 * x48 * x52,
            x53 * x54,
            x23 * x39 * x55,
            x23 * x41 * x53,
            x14 * x39 * x56,
            x14 * x49 * x55,
            x14 * x45 * x53,
            x59 * x9,
            x33 * x56 * x60,
            x10 * x44 * x55,
            x52 * x53 * x62,
            x54 * x63,
            x23 * x30 * x64,
            x23 * x35 * x65,
            x14 * x40 * x63,
            x14 * x50 * x65,
            x14 * x35 * x66,
            x47 * x60 * x63,
            x10 * x38 * x65,
            x30 * x62 * x66,
            x68 * x9,
            x69 * x70,
            x71 * x74,
            x41 * x69 * x71,
            x39 * x78 * x79,
            x49 * x73 * x79,
            x44 * x75 * x79,
            x20 * x82,
            x20 * x78 * x83,
            x20 * x73 * x84,
            x20 * x52 * x85,
            x19 * x39 * x86,
            x17 * x39 * x87,
            x17 * x49 * x86,
            x3 * x88,
            x3 * x83 * x87,
            x44 * x8 * x86,
            x58 * (x0 * (x47 + 3.0 * x76 + 3.0 * x77) + x53 * x81),
            x33 * x88,
            x84 * x87,
            x52 * x7 * x86,
            x19 * x64 * x69,
            x17 * x63 * x74,
            x17 * x65 * x75,
            x63 * x78 * x89,
            x65 * x73 * x8,
            x3 * x66 * x85,
            x63 * x82,
            x65 * x7 * x78,
            x66 * x7 * x73,
            x68 * x69,
            x70 * x90,
            x30 * x71 * x91,
            x71 * x94,
            x38 * x79 * x95,
            x50 * x79 * x93,
            x35 * x79 * x98,
            x20 * x47 * x58 * x90,
            x20 * x38 * x99,
            x100 * x20 * x30,
            x103 * x20,
            x19 * x53 * x91,
            x17 * x55 * x95,
            x17 * x53 * x94,
            x56 * x89 * x90,
            x55 * x8 * x93,
            x100 * x3 * x53,
            x59 * x90,
            x56 * x99,
            x55 * x7 * x98,
            x103 * x53,
            x104 * x19 * x35,
            x104 * x17 * x50,
            x105 * x17 * x35,
            x104 * x38 * x8,
            x105 * x3 * x30 * x61,
            x106 * x3,
            x104 * x47 * x7,
            x105 * x38 * x7,
            x106 * x30,
            x61 * (x0 * (x52 + 3.0 * x96 + 3.0 * x97) + x102 * x63),
        ]
    )


def dipole3d_14(a, A, b, B, C):
    """Cartesian 3D (pg) dipole moment integrals.
    The origin is at C.

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
    x14 = x0 * (x10 + 2.0 * x13 + x8)
    x15 = x11 * x7
    x16 = x0 * (x12 + x15)
    x17 = x13 + x9
    x18 = x17 * x3
    x19 = x16 + x18
    x20 = x19 * x3
    x21 = x0 * (x10 + 3.0 * x8)
    x22 = x8 + x9
    x23 = x22 * x3
    x24 = x3 * x9
    x25 = 2.0 * x24
    x26 = x23 + x25
    x27 = x21 + x26 * x3
    x28 = -x2 - A[0]
    x29 = x0 * (3.0 * x16 + 3.0 * x18 + x26)
    x30 = x14 + x20
    x31 = x29 + x3 * x30
    x32 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x33 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x34 = 3.14159265358979 * x1 * x33
    x35 = x32 * x34
    x36 = -x1 * (a * A[1] + b * B[1])
    x37 = -x36 - B[1]
    x38 = x35 * (x28 * x30 + x29)
    x39 = -x1 * (a * A[2] + b * B[2])
    x40 = -x39 - B[2]
    x41 = x14 + x19 * x28
    x42 = x32 * x6
    x43 = x37**2 * x42
    x44 = x0 * x42
    x45 = x43 + x44
    x46 = x33 * x6
    x47 = x45 * x46
    x48 = x35 * x40
    x49 = x40**2 * x46
    x50 = x0 * x46
    x51 = x49 + x50
    x52 = x42 * x51
    x53 = x16 + x17 * x28
    x54 = x37 * x45
    x55 = x37 * x44
    x56 = 2.0 * x55
    x57 = x54 + x56
    x58 = x46 * x57
    x59 = x40 * x46
    x60 = x37 * x42
    x61 = x40 * x51
    x62 = x40 * x50
    x63 = 2.0 * x62
    x64 = x61 + x63
    x65 = x42 * x64
    x66 = 3.0 * x44
    x67 = x0 * (3.0 * x43 + x66)
    x68 = x37 * x57 + x67
    x69 = x15 * x28 + x9
    x70 = 3.0 * x50
    x71 = x0 * (3.0 * x49 + x70)
    x72 = x40 * x64 + x71
    x73 = -x36 - A[1]
    x74 = x31 * x35
    x75 = x44 + x60 * x73
    x76 = x45 * x73 + x56
    x77 = x57 * x73 + x67
    x78 = x34 * x5
    x79 = x78 * (x0 * (4.0 * x54 + 8.0 * x55) + x68 * x73)
    x80 = x11 * x78
    x81 = 3.14159265358979 * x1 * x32 * x5
    x82 = x11 * x81
    x83 = -x39 - A[2]
    x84 = x35 * x83
    x85 = x50 + x59 * x83
    x86 = x51 * x83 + x63
    x87 = x64 * x83 + x71
    x88 = x81 * (x0 * (4.0 * x61 + 8.0 * x62) + x72 * x83)
    x89 = -x36 - C[1]
    x90 = x35 * (x0 * (4.0 * x23 + 8.0 * x24) + x27 * x28)
    x91 = x21 + x26 * x28
    x92 = x60 * x89
    x93 = x44 + x92
    x94 = x46 * x93
    x95 = x22 * x28 + x25
    x96 = x42 * x89
    x97 = x0 * (x60 + x96)
    x98 = x37 * x93
    x99 = x97 + x98
    x100 = x46 * x99
    x101 = x0 * (x43 + x66 + 2.0 * x92)
    x102 = x37 * x99
    x103 = x101 + x102
    x104 = x12 * x28 + x9
    x105 = x0 * (x57 + 3.0 * x97 + 3.0 * x98)
    x106 = x103 * x37 + x105
    x107 = x106 * x78
    x108 = x40 * x78
    x109 = x51 * x7
    x110 = x64 * x7
    x111 = x81 * x89
    x112 = x44 + x73 * x96
    x113 = x73 * x93 + x97
    x114 = x101 + x73 * x99
    x115 = x78 * (x103 * x73 + x105)
    x116 = x3 * x78
    x117 = -x39 - C[2]
    x118 = x117 * x35
    x119 = x117 * x59
    x120 = x119 + x50
    x121 = x120 * x42
    x122 = x117 * x46
    x123 = x0 * (x122 + x59)
    x124 = x120 * x40
    x125 = x123 + x124
    x126 = x125 * x42
    x127 = x0 * (2.0 * x119 + x49 + x70)
    x128 = x125 * x40
    x129 = x127 + x128
    x130 = x120 * x7
    x131 = x125 * x7
    x132 = x129 * x81
    x133 = x0 * (3.0 * x123 + 3.0 * x124 + x64)
    x134 = x129 * x40 + x133
    x135 = x134 * x81
    x136 = x122 * x83 + x50
    x137 = x120 * x83 + x123
    x138 = x125 * x83 + x127
    x139 = x81 * (x129 * x83 + x133)

    # 135 item(s)
    return numpy.array(
        [
            x35 * (x0 * (4.0 * x14 + 4.0 * x20 + x27) + x28 * x31),
            x37 * x38,
            x38 * x40,
            x41 * x47,
            x37 * x41 * x48,
            x41 * x52,
            x53 * x58,
            x45 * x53 * x59,
            x51 * x53 * x60,
            x53 * x65,
            x46 * x68 * x69,
            x57 * x59 * x69,
            x45 * x51 * x69,
            x60 * x64 * x69,
            x42 * x69 * x72,
            x73 * x74,
            x30 * x46 * x75,
            x30 * x48 * x73,
            x19 * x46 * x76,
            x19 * x59 * x75,
            x19 * x52 * x73,
            x17 * x46 * x77,
            x17 * x59 * x76,
            x17 * x51 * x75,
            x17 * x65 * x73,
            x11 * x79,
            x40 * x77 * x80,
            x15 * x51 * x76,
            x15 * x64 * x75,
            x72 * x73 * x82,
            x74 * x83,
            x30 * x37 * x84,
            x30 * x42 * x85,
            x19 * x47 * x83,
            x19 * x60 * x85,
            x19 * x42 * x86,
            x17 * x58 * x83,
            x17 * x45 * x85,
            x17 * x60 * x86,
            x17 * x42 * x87,
            x68 * x80 * x83,
            x15 * x57 * x85,
            x15 * x45 * x86,
            x37 * x82 * x87,
            x11 * x88,
            x89 * x90,
            x91 * x94,
            x48 * x89 * x91,
            x100 * x95,
            x59 * x93 * x95,
            x51 * x95 * x96,
            x103 * x104 * x46,
            x104 * x59 * x99,
            x104 * x51 * x93,
            x104 * x64 * x96,
            x107 * x28,
            x103 * x108 * x28,
            x109 * x28 * x99,
            x110 * x28 * x93,
            x111 * x28 * x72,
            x112 * x27 * x46,
            x113 * x26 * x46,
            x112 * x26 * x59,
            x114 * x22 * x46,
            x113 * x22 * x59,
            x112 * x22 * x51,
            x115 * x3,
            x108 * x114 * x3,
            x113 * x12 * x51,
            x112 * x12 * x64,
            x78 * (x0 * (4.0 * x101 + 4.0 * x102 + x68) + x106 * x73),
            x115 * x40,
            x109 * x114,
            x110 * x113,
            x112 * x7 * x72,
            x27 * x84 * x89,
            x26 * x83 * x94,
            x26 * x85 * x96,
            x100 * x22 * x83,
            x22 * x85 * x93,
            x22 * x86 * x96,
            x103 * x116 * x83,
            x12 * x85 * x99,
            x12 * x86 * x93,
            x111 * x3 * x87,
            x107 * x83,
            x103 * x7 * x85,
            x7 * x86 * x99,
            x7 * x87 * x93,
            x88 * x89,
            x117 * x90,
            x118 * x37 * x91,
            x121 * x91,
            x122 * x45 * x95,
            x120 * x60 * x95,
            x126 * x95,
            x104 * x122 * x57,
            x104 * x120 * x45,
            x104 * x125 * x60,
            x104 * x129 * x42,
            x117 * x28 * x68 * x78,
            x130 * x28 * x57,
            x131 * x28 * x45,
            x132 * x28 * x37,
            x135 * x28,
            x118 * x27 * x73,
            x122 * x26 * x75,
            x121 * x26 * x73,
            x122 * x22 * x76,
            x120 * x22 * x75,
            x126 * x22 * x73,
            x116 * x117 * x77,
            x12 * x120 * x76,
            x12 * x125 * x75,
            x132 * x3 * x73,
            x117 * x79,
            x130 * x77,
            x131 * x76,
            x129 * x7 * x75,
            x135 * x73,
            x136 * x27 * x42,
            x136 * x26 * x60,
            x137 * x26 * x42,
            x136 * x22 * x45,
            x137 * x22 * x60,
            x138 * x22 * x42,
            x12 * x136 * x57,
            x12 * x137 * x45,
            x138 * x3 * x37 * x81,
            x139 * x3,
            x136 * x68 * x7,
            x137 * x57 * x7,
            x138 * x45 * x7,
            x139 * x37,
            x81 * (x0 * (4.0 * x127 + 4.0 * x128 + x72) + x134 * x83),
        ]
    )


def dipole3d_20(a, A, b, B, C):
    """Cartesian 3D (ds) dipole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = numpy.sqrt(x1)
    x3 = 1.77245385090552 * x2
    x4 = -x1 * (a * A[0] + b * B[0])
    x5 = -x4 - A[0]
    x6 = a * b * x1
    x7 = numpy.exp(-x6 * (A[0] - B[0]) ** 2)
    x8 = x5 * x7
    x9 = x3 * x8
    x10 = -x4 - C[0]
    x11 = x3 * x7
    x12 = x0 * x11
    x13 = x10 * x9 + x12
    x14 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x15 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x16 = 3.14159265358979 * x1 * x15
    x17 = x14 * x16
    x18 = -x1 * (a * A[1] + b * B[1])
    x19 = -x18 - A[1]
    x20 = x13 * x17
    x21 = -x1 * (a * A[2] + b * B[2])
    x22 = -x21 - A[2]
    x23 = x14 * x3
    x24 = x0 * x23
    x25 = x16 * x7
    x26 = x25 * (x19**2 * x23 + x24)
    x27 = 5.56832799683171
    x28 = x1 * x14
    x29 = x15 * x2 * x22 * x27 * x28
    x30 = x15 * x3
    x31 = x0 * x30
    x32 = 3.14159265358979 * x28
    x33 = x32 * x7
    x34 = x33 * (x22**2 * x30 + x31)
    x35 = -x18 - C[1]
    x36 = x17 * (x11 * x5**2 + x12)
    x37 = x19 * x23
    x38 = x24 + x35 * x37
    x39 = -x21 - C[2]
    x40 = x22 * x30
    x41 = x31 + x39 * x40

    # 18 item(s)
    return numpy.array(
        [
            x17 * (x0 * (x10 * x11 + x9) + x13 * x5),
            x19 * x20,
            x20 * x22,
            x10 * x26,
            x10 * x19 * x29 * x7,
            x10 * x34,
            x35 * x36,
            x16 * x38 * x8,
            x29 * x35 * x8,
            x25 * (x0 * (x23 * x35 + x37) + x19 * x38),
            x22 * x25 * x38,
            x34 * x35,
            x36 * x39,
            x15 * x19 * x2 * x27 * x28 * x39 * x8,
            x32 * x41 * x8,
            x26 * x39,
            x19 * x33 * x41,
            x33 * (x0 * (x30 * x39 + x40) + x22 * x41),
        ]
    )


def dipole3d_21(a, A, b, B, C):
    """Cartesian 3D (dp) dipole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = a * b * x1
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = 1.77245385090552 * numpy.sqrt(x1)
    x5 = x3 * x4
    x6 = x0 * x5
    x7 = -x1 * (a * A[0] + b * B[0])
    x8 = -x7 - A[0]
    x9 = x3 * (-x7 - B[0])
    x10 = x4 * x9
    x11 = x10 * x8
    x12 = -x7 - C[0]
    x13 = x12 * x5
    x14 = x13 * x8
    x15 = x10 * x12
    x16 = x15 + x6
    x17 = x0 * (x10 + x13) + x16 * x8
    x18 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x19 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x20 = 3.14159265358979 * x1 * x19
    x21 = x18 * x20
    x22 = -x1 * (a * A[1] + b * B[1])
    x23 = -x22 - B[1]
    x24 = x5 * x8
    x25 = x14 + x6
    x26 = x21 * (x0 * (x13 + x24) + x25 * x8)
    x27 = -x1 * (a * A[2] + b * B[2])
    x28 = -x27 - B[2]
    x29 = -x22 - A[1]
    x30 = x17 * x21
    x31 = x0 * x4
    x32 = x18 * x31
    x33 = x18 * x4
    x34 = x29 * x33
    x35 = x23 * x34
    x36 = x32 + x35
    x37 = x19 * x4
    x38 = x21 * x25
    x39 = -x27 - A[2]
    x40 = x19 * x31
    x41 = x37 * x39
    x42 = x28 * x41
    x43 = x40 + x42
    x44 = x29**2 * x33 + x32
    x45 = x23 * x33
    x46 = x20 * x3
    x47 = x46 * (x0 * (x34 + x45) + x29 * x36)
    x48 = x28 * x46
    x49 = x21 * x39
    x50 = x39 * x46
    x51 = 3.14159265358979 * x1 * x18
    x52 = x3 * x51
    x53 = x43 * x52
    x54 = x37 * x39**2 + x40
    x55 = x23 * x52
    x56 = x28 * x37
    x57 = x52 * (x0 * (x41 + x56) + x39 * x43)
    x58 = -x22 - C[1]
    x59 = x11 + x6
    x60 = x21 * (x0 * (x10 + x24) + x59 * x8)
    x61 = x45 * x58
    x62 = x32 + x61
    x63 = x5 * x8**2 + x6
    x64 = x21 * x63
    x65 = x34 * x58
    x66 = x32 + x65
    x67 = x33 * x58
    x68 = x0 * (x45 + x67) + x29 * x62
    x69 = x46 * x68
    x70 = x0 * (x34 + x67) + x29 * x66
    x71 = x20 * x9
    x72 = x51 * x9
    x73 = -x27 - C[2]
    x74 = x56 * x73
    x75 = x40 + x74
    x76 = x41 * x73
    x77 = x40 + x76
    x78 = x37 * x73
    x79 = x0 * (x56 + x78) + x39 * x75
    x80 = x52 * x79
    x81 = x0 * (x41 + x78) + x39 * x77

    # 54 item(s)
    return numpy.array(
        [
            x21 * (x0 * (x11 + x14 + x15 + 3.0 * x6) + x17 * x8),
            x23 * x26,
            x26 * x28,
            x29 * x30,
            x25 * x36 * x37,
            x28 * x29 * x38,
            x30 * x39,
            x23 * x38 * x39,
            x25 * x33 * x43,
            x16 * x37 * x44,
            x12 * x47,
            x12 * x44 * x48,
            x16 * x29 * x49,
            x12 * x36 * x50,
            x12 * x29 * x53,
            x16 * x33 * x54,
            x12 * x54 * x55,
            x12 * x57,
            x58 * x60,
            x37 * x62 * x63,
            x28 * x58 * x64,
            x37 * x59 * x66,
            x69 * x8,
            x48 * x66 * x8,
            x49 * x58 * x59,
            x50 * x62 * x8,
            x53 * x58 * x8,
            x70 * x71,
            x46 * (x0 * (3.0 * x32 + x35 + x61 + x65) + x29 * x68),
            x48 * x70,
            x39 * x66 * x71,
            x39 * x69,
            x43 * x5 * x66,
            x54 * x58 * x72,
            x5 * x54 * x62,
            x57 * x58,
            x60 * x73,
            x23 * x64 * x73,
            x33 * x63 * x75,
            x21 * x29 * x59 * x73,
            x36 * x46 * x73 * x8,
            x29 * x52 * x75 * x8,
            x33 * x59 * x77,
            x55 * x77 * x8,
            x8 * x80,
            x44 * x71 * x73,
            x47 * x73,
            x44 * x5 * x75,
            x29 * x72 * x77,
            x36 * x5 * x77,
            x29 * x80,
            x72 * x81,
            x55 * x81,
            x52 * (x0 * (3.0 * x40 + x42 + x74 + x76) + x39 * x79),
        ]
    )


def dipole3d_22(a, A, b, B, C):
    """Cartesian 3D (dd) dipole moment integrals.
    The origin is at C.

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
    x14 = -x2 - A[0]
    x15 = x12 * x14
    x16 = x7 * x9
    x17 = x0 * (x10 + x16)
    x18 = x3**2 * x7
    x19 = x18 + x8
    x20 = x14 * x19 + 2.0 * x3 * x8
    x21 = 3.0 * x8
    x22 = x18 + x21
    x23 = x13 + x17
    x24 = x0 * (2.0 * x11 + x22) + x14 * x23
    x25 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x26 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x27 = 3.14159265358979 * x1 * x26
    x28 = x25 * x27
    x29 = -x1 * (a * A[1] + b * B[1])
    x30 = -x29 - B[1]
    x31 = x10 * x14
    x32 = x14 * x16
    x33 = x15 + x17
    x34 = x28 * (x0 * (x11 + x21 + x31 + x32) + x14 * x33)
    x35 = -x1 * (a * A[2] + b * B[2])
    x36 = -x35 - B[2]
    x37 = x25 * x6
    x38 = x30**2 * x37
    x39 = x0 * x37
    x40 = x38 + x39
    x41 = x14 * x7
    x42 = x32 + x8
    x43 = x0 * (x16 + x41) + x14 * x42
    x44 = x26 * x6
    x45 = x28 * x36
    x46 = x36**2 * x44
    x47 = x0 * x44
    x48 = x46 + x47
    x49 = -x29 - A[1]
    x50 = x24 * x28
    x51 = x37 * x49
    x52 = x30 * x51
    x53 = x39 + x52
    x54 = 2.0 * x30 * x39 + x40 * x49
    x55 = x36 * x44
    x56 = -x35 - A[2]
    x57 = x28 * x56
    x58 = x44 * x56
    x59 = x36 * x58
    x60 = x47 + x59
    x61 = x30 * x37
    x62 = 2.0 * x36 * x47 + x48 * x56
    x63 = x37 * x49**2 + x39
    x64 = x0 * (x51 + x61) + x49 * x53
    x65 = 3.0 * x39
    x66 = x38 + x65
    x67 = x27 * x5
    x68 = x67 * (x0 * (2.0 * x52 + x66) + x49 * x54)
    x69 = x67 * x9
    x70 = 3.14159265358979 * x1 * x25 * x5
    x71 = x70 * x9
    x72 = x44 * x56**2 + x47
    x73 = x0 * (x55 + x58) + x56 * x60
    x74 = 3.0 * x47
    x75 = x46 + x74
    x76 = x70 * (x0 * (2.0 * x59 + x75) + x56 * x62)
    x77 = -x29 - C[1]
    x78 = x28 * (x0 * (x22 + 2.0 * x31) + x14 * x20)
    x79 = x61 * x77
    x80 = x39 + x79
    x81 = x31 + x8
    x82 = x0 * (x10 + x41) + x14 * x81
    x83 = x37 * x77
    x84 = x0 * (x61 + x83)
    x85 = x30 * x80
    x86 = x84 + x85
    x87 = x14**2 * x7 + x8
    x88 = x51 * x77
    x89 = x39 + x88
    x90 = x49 * x80
    x91 = x84 + x90
    x92 = x0 * (x66 + 2.0 * x79) + x49 * x86
    x93 = x67 * x92
    x94 = x14 * x67
    x95 = x70 * x77
    x96 = x0 * (x51 + x83) + x49 * x89
    x97 = x67 * (x0 * (x52 + x65 + x79 + x88) + x49 * x91)
    x98 = x3 * x67
    x99 = -x35 - C[2]
    x100 = x28 * x99
    x101 = x55 * x99
    x102 = x101 + x47
    x103 = x44 * x99
    x104 = x0 * (x103 + x55)
    x105 = x102 * x36
    x106 = x104 + x105
    x107 = x14 * x70
    x108 = x58 * x99
    x109 = x108 + x47
    x110 = x102 * x56
    x111 = x104 + x110
    x112 = x0 * (2.0 * x101 + x75) + x106 * x56
    x113 = x112 * x70
    x114 = x3 * x70
    x115 = x0 * (x103 + x58) + x109 * x56
    x116 = x70 * (x0 * (x101 + x108 + x59 + x74) + x111 * x56)

    # 108 item(s)
    return numpy.array(
        [
            x28 * (x0 * (x13 + 2.0 * x15 + 3.0 * x17 + x20) + x14 * x24),
            x30 * x34,
            x34 * x36,
            x40 * x43 * x44,
            x30 * x43 * x45,
            x37 * x43 * x48,
            x49 * x50,
            x33 * x44 * x53,
            x33 * x45 * x49,
            x42 * x44 * x54,
            x42 * x53 * x55,
            x42 * x48 * x51,
            x50 * x56,
            x30 * x33 * x57,
            x33 * x37 * x60,
            x40 * x42 * x58,
            x42 * x60 * x61,
            x37 * x42 * x62,
            x23 * x44 * x63,
            x12 * x44 * x64,
            x12 * x55 * x63,
            x68 * x9,
            x36 * x64 * x69,
            x16 * x48 * x63,
            x23 * x49 * x57,
            x12 * x53 * x58,
            x12 * x51 * x60,
            x54 * x56 * x69,
            x16 * x53 * x60,
            x49 * x62 * x71,
            x23 * x37 * x72,
            x12 * x61 * x72,
            x12 * x37 * x73,
            x16 * x40 * x72,
            x30 * x71 * x73,
            x76 * x9,
            x77 * x78,
            x44 * x80 * x82,
            x45 * x77 * x82,
            x44 * x86 * x87,
            x55 * x80 * x87,
            x48 * x83 * x87,
            x20 * x44 * x89,
            x44 * x81 * x91,
            x55 * x81 * x89,
            x14 * x93,
            x36 * x91 * x94,
            x41 * x48 * x89,
            x20 * x57 * x77,
            x58 * x80 * x81,
            x60 * x81 * x83,
            x56 * x86 * x94,
            x41 * x60 * x80,
            x14 * x62 * x95,
            x19 * x44 * x96,
            x3 * x97,
            x36 * x96 * x98,
            x67 * (x0 * (x54 + 3.0 * x84 + x85 + 2.0 * x90) + x49 * x92),
            x36 * x97,
            x48 * x7 * x96,
            x19 * x58 * x89,
            x56 * x91 * x98,
            x10 * x60 * x89,
            x56 * x93,
            x60 * x7 * x91,
            x62 * x7 * x89,
            x19 * x72 * x83,
            x10 * x72 * x80,
            x3 * x73 * x95,
            x7 * x72 * x86,
            x7 * x73 * x80,
            x76 * x77,
            x78 * x99,
            x100 * x30 * x82,
            x102 * x37 * x82,
            x103 * x40 * x87,
            x102 * x61 * x87,
            x106 * x37 * x87,
            x100 * x20 * x49,
            x103 * x53 * x81,
            x102 * x51 * x81,
            x54 * x94 * x99,
            x102 * x41 * x53,
            x106 * x107 * x49,
            x109 * x20 * x37,
            x109 * x61 * x81,
            x111 * x37 * x81,
            x109 * x40 * x41,
            x107 * x111 * x30,
            x113 * x14,
            x103 * x19 * x63,
            x64 * x98 * x99,
            x10 * x102 * x63,
            x68 * x99,
            x102 * x64 * x7,
            x106 * x63 * x7,
            x109 * x19 * x51,
            x10 * x109 * x53,
            x111 * x114 * x49,
            x109 * x54 * x7,
            x111 * x53 * x7,
            x113 * x49,
            x115 * x19 * x37,
            x114 * x115 * x30,
            x116 * x3,
            x115 * x40 * x7,
            x116 * x30,
            x70 * (x0 * (3.0 * x104 + x105 + 2.0 * x110 + x62) + x112 * x56),
        ]
    )


def dipole3d_23(a, A, b, B, C):
    """Cartesian 3D (df) dipole moment integrals.
    The origin is at C.

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
    x18 = -x2 - A[0]
    x19 = x16 * x18
    x20 = 3.0 * x12
    x21 = x3**2 * x7
    x22 = x20 + x21
    x23 = x0 * (2.0 * x13 + x22)
    x24 = x12 + x21
    x25 = x24 * x3
    x26 = x12 * x3
    x27 = 2.0 * x26
    x28 = x25 + x27
    x29 = x0 * (x20 + 3.0 * x21) + x18 * x28
    x30 = 3.0 * x11
    x31 = x17 + x23
    x32 = x0 * (3.0 * x15 + x28 + x30) + x18 * x31
    x33 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x34 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x35 = 3.14159265358979 * x1 * x34
    x36 = x33 * x35
    x37 = -x1 * (a * A[1] + b * B[1])
    x38 = -x37 - B[1]
    x39 = x14 * x18
    x40 = x18 * x24
    x41 = x27 + x40
    x42 = x19 + x23
    x43 = x36 * (x0 * (x15 + x30 + 2.0 * x39 + x41) + x18 * x42)
    x44 = -x1 * (a * A[2] + b * B[2])
    x45 = -x44 - B[2]
    x46 = x33 * x6
    x47 = x38**2 * x46
    x48 = x0 * x46
    x49 = x47 + x48
    x50 = x18 * x8
    x51 = x10 * x18
    x52 = x11 + x39
    x53 = x0 * (x13 + x20 + x50 + x51) + x18 * x52
    x54 = x34 * x6
    x55 = x36 * x45
    x56 = x45**2 * x54
    x57 = x0 * x54
    x58 = x56 + x57
    x59 = x38 * x49
    x60 = x38 * x48
    x61 = 2.0 * x60
    x62 = x59 + x61
    x63 = x18 * x7
    x64 = x12 + x51
    x65 = x0 * (x10 + x63) + x18 * x64
    x66 = x45 * x54
    x67 = x38 * x46
    x68 = x45 * x58
    x69 = x45 * x57
    x70 = 2.0 * x69
    x71 = x68 + x70
    x72 = -x37 - A[1]
    x73 = x32 * x36
    x74 = x46 * x72
    x75 = x38 * x74
    x76 = x48 + x75
    x77 = x49 * x72
    x78 = x61 + x77
    x79 = 3.0 * x48
    x80 = x0 * (3.0 * x47 + x79) + x62 * x72
    x81 = -x44 - A[2]
    x82 = x36 * x81
    x83 = x54 * x81
    x84 = x45 * x83
    x85 = x57 + x84
    x86 = x58 * x81
    x87 = x70 + x86
    x88 = 3.0 * x57
    x89 = x0 * (3.0 * x56 + x88) + x71 * x81
    x90 = x46 * x72**2 + x48
    x91 = x0 * (x67 + x74) + x72 * x76
    x92 = x47 + x79
    x93 = x0 * (2.0 * x75 + x92) + x72 * x78
    x94 = x35 * x5
    x95 = x94 * (x0 * (x59 + 8.0 * x60 + 3.0 * x77) + x72 * x80)
    x96 = x9 * x94
    x97 = 3.14159265358979 * x1 * x33 * x5
    x98 = x9 * x97
    x99 = x54 * x81**2 + x57
    x100 = x0 * (x66 + x83) + x81 * x85
    x101 = x56 + x88
    x102 = x0 * (x101 + 2.0 * x84) + x81 * x87
    x103 = x97 * (x0 * (x68 + 8.0 * x69 + 3.0 * x86) + x81 * x89)
    x104 = -x37 - C[1]
    x105 = x36 * (x0 * (x25 + 8.0 * x26 + 3.0 * x40) + x18 * x29)
    x106 = x104 * x67
    x107 = x106 + x48
    x108 = x0 * (x22 + 2.0 * x50) + x18 * x41
    x109 = x104 * x46
    x110 = x0 * (x109 + x67)
    x111 = x107 * x38
    x112 = x110 + x111
    x113 = x12 + x50
    x114 = x0 * (x63 + x8) + x113 * x18
    x115 = x0 * (2.0 * x106 + x92)
    x116 = x112 * x38
    x117 = x115 + x116
    x118 = x12 + x18**2 * x7
    x119 = x104 * x74
    x120 = x119 + x48
    x121 = x107 * x72
    x122 = x110 + x121
    x123 = x112 * x72
    x124 = x115 + x123
    x125 = 3.0 * x110
    x126 = x0 * (3.0 * x111 + x125 + x62) + x117 * x72
    x127 = x126 * x94
    x128 = x18 * x94
    x129 = x104 * x97
    x130 = x0 * (x109 + x74) + x120 * x72
    x131 = x0 * (x106 + x119 + x75 + x79) + x122 * x72
    x132 = x94 * (x0 * (x111 + 2.0 * x121 + x125 + x78) + x124 * x72)
    x133 = x3 * x94
    x134 = -x44 - C[2]
    x135 = x134 * x36
    x136 = x134 * x66
    x137 = x136 + x57
    x138 = x134 * x54
    x139 = x0 * (x138 + x66)
    x140 = x137 * x45
    x141 = x139 + x140
    x142 = x0 * (x101 + 2.0 * x136)
    x143 = x141 * x45
    x144 = x142 + x143
    x145 = x18 * x97
    x146 = x134 * x83
    x147 = x146 + x57
    x148 = x137 * x81
    x149 = x139 + x148
    x150 = x141 * x81
    x151 = x142 + x150
    x152 = 3.0 * x139
    x153 = x0 * (3.0 * x140 + x152 + x71) + x144 * x81
    x154 = x153 * x97
    x155 = x3 * x97
    x156 = x0 * (x138 + x83) + x147 * x81
    x157 = x0 * (x136 + x146 + x84 + x88) + x149 * x81
    x158 = x97 * (x0 * (x140 + 2.0 * x148 + x152 + x87) + x151 * x81)

    # 180 item(s)
    return numpy.array(
        [
            x36 * (x0 * (x17 + 3.0 * x19 + 4.0 * x23 + x29) + x18 * x32),
            x38 * x43,
            x43 * x45,
            x49 * x53 * x54,
            x38 * x53 * x55,
            x46 * x53 * x58,
            x54 * x62 * x65,
            x49 * x65 * x66,
            x58 * x65 * x67,
            x46 * x65 * x71,
            x72 * x73,
            x42 * x54 * x76,
            x42 * x55 * x72,
            x52 * x54 * x78,
            x52 * x66 * x76,
            x52 * x58 * x74,
            x54 * x64 * x80,
            x64 * x66 * x78,
            x58 * x64 * x76,
            x64 * x71 * x74,
            x73 * x81,
            x38 * x42 * x82,
            x42 * x46 * x85,
            x49 * x52 * x83,
            x52 * x67 * x85,
            x46 * x52 * x87,
            x62 * x64 * x83,
            x49 * x64 * x85,
            x64 * x67 * x87,
            x46 * x64 * x89,
            x31 * x54 * x90,
            x16 * x54 * x91,
            x16 * x66 * x90,
            x14 * x54 * x93,
            x14 * x66 * x91,
            x14 * x58 * x90,
            x9 * x95,
            x45 * x93 * x96,
            x10 * x58 * x91,
            x10 * x71 * x90,
            x31 * x72 * x82,
            x16 * x76 * x83,
            x16 * x74 * x85,
            x14 * x78 * x83,
            x14 * x76 * x85,
            x14 * x74 * x87,
            x80 * x81 * x96,
            x10 * x78 * x85,
            x10 * x76 * x87,
            x72 * x89 * x98,
            x31 * x46 * x99,
            x16 * x67 * x99,
            x100 * x16 * x46,
            x14 * x49 * x99,
            x100 * x14 * x67,
            x102 * x14 * x46,
            x10 * x62 * x99,
            x10 * x100 * x49,
            x102 * x38 * x98,
            x103 * x9,
            x104 * x105,
            x107 * x108 * x54,
            x104 * x108 * x55,
            x112 * x114 * x54,
            x107 * x114 * x66,
            x109 * x114 * x58,
            x117 * x118 * x54,
            x112 * x118 * x66,
            x107 * x118 * x58,
            x109 * x118 * x71,
            x120 * x29 * x54,
            x122 * x41 * x54,
            x120 * x41 * x66,
            x113 * x124 * x54,
            x113 * x122 * x66,
            x113 * x120 * x58,
            x127 * x18,
            x124 * x128 * x45,
            x122 * x58 * x63,
            x120 * x63 * x71,
            x104 * x29 * x82,
            x107 * x41 * x83,
            x109 * x41 * x85,
            x112 * x113 * x83,
            x107 * x113 * x85,
            x109 * x113 * x87,
            x117 * x128 * x81,
            x112 * x63 * x85,
            x107 * x63 * x87,
            x129 * x18 * x89,
            x130 * x28 * x54,
            x131 * x24 * x54,
            x130 * x24 * x66,
            x132 * x3,
            x131 * x133 * x45,
            x130 * x58 * x8,
            x94 * (x0 * (4.0 * x115 + x116 + 3.0 * x123 + x80) + x126 * x72),
            x132 * x45,
            x131 * x58 * x7,
            x130 * x7 * x71,
            x120 * x28 * x83,
            x122 * x24 * x83,
            x120 * x24 * x85,
            x124 * x133 * x81,
            x122 * x8 * x85,
            x120 * x8 * x87,
            x127 * x81,
            x124 * x7 * x85,
            x122 * x7 * x87,
            x120 * x7 * x89,
            x109 * x28 * x99,
            x107 * x24 * x99,
            x100 * x109 * x24,
            x112 * x8 * x99,
            x100 * x107 * x8,
            x102 * x129 * x3,
            x117 * x7 * x99,
            x100 * x112 * x7,
            x102 * x107 * x7,
            x103 * x104,
            x105 * x134,
            x108 * x135 * x38,
            x108 * x137 * x46,
            x114 * x138 * x49,
            x114 * x137 * x67,
            x114 * x141 * x46,
            x118 * x138 * x62,
            x118 * x137 * x49,
            x118 * x141 * x67,
            x118 * x144 * x46,
            x135 * x29 * x72,
            x138 * x41 * x76,
            x137 * x41 * x74,
            x113 * x138 * x78,
            x113 * x137 * x76,
            x113 * x141 * x74,
            x128 * x134 * x80,
            x137 * x63 * x78,
            x141 * x63 * x76,
            x144 * x145 * x72,
            x147 * x29 * x46,
            x147 * x41 * x67,
            x149 * x41 * x46,
            x113 * x147 * x49,
            x113 * x149 * x67,
            x113 * x151 * x46,
            x147 * x62 * x63,
            x149 * x49 * x63,
            x145 * x151 * x38,
            x154 * x18,
            x138 * x28 * x90,
            x138 * x24 * x91,
            x137 * x24 * x90,
            x133 * x134 * x93,
            x137 * x8 * x91,
            x141 * x8 * x90,
            x134 * x95,
            x137 * x7 * x93,
            x141 * x7 * x91,
            x144 * x7 * x90,
            x147 * x28 * x74,
            x147 * x24 * x76,
            x149 * x24 * x74,
            x147 * x78 * x8,
            x149 * x76 * x8,
            x151 * x155 * x72,
            x147 * x7 * x80,
            x149 * x7 * x78,
            x151 * x7 * x76,
            x154 * x72,
            x156 * x28 * x46,
            x156 * x24 * x67,
            x157 * x24 * x46,
            x156 * x49 * x8,
            x155 * x157 * x38,
            x158 * x3,
            x156 * x62 * x7,
            x157 * x49 * x7,
            x158 * x38,
            x97 * (x0 * (4.0 * x142 + x143 + 3.0 * x150 + x89) + x153 * x81),
        ]
    )


def dipole3d_24(a, A, b, B, C):
    """Cartesian 3D (dg) dipole moment integrals.
    The origin is at C.

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
    x9 = x3 * x8
    x10 = x4 * x9
    x11 = x3**2 * x8
    x12 = x0 * x8
    x13 = 3.0 * x12
    x14 = x11 + x13
    x15 = x0 * (2.0 * x10 + x14)
    x16 = x4 * x8
    x17 = x0 * (x16 + x9)
    x18 = x10 + x12
    x19 = x18 * x3
    x20 = x17 + x19
    x21 = x20 * x3
    x22 = x15 + x21
    x23 = x22 * x3
    x24 = -x2 - A[0]
    x25 = x22 * x24
    x26 = 3.0 * x17
    x27 = x11 + x12
    x28 = x27 * x3
    x29 = x12 * x3
    x30 = 2.0 * x29
    x31 = x28 + x30
    x32 = x0 * (3.0 * x19 + x26 + x31)
    x33 = 8.0 * x29
    x34 = x0 * (3.0 * x11 + x13)
    x35 = x3 * x31
    x36 = x34 + x35
    x37 = x0 * (4.0 * x28 + x33) + x24 * x36
    x38 = 4.0 * x15
    x39 = x23 + x32
    x40 = x0 * (4.0 * x21 + x36 + x38) + x24 * x39
    x41 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x42 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x43 = 3.14159265358979 * x1 * x42
    x44 = x41 * x43
    x45 = -x1 * (a * A[1] + b * B[1])
    x46 = -x45 - B[1]
    x47 = x20 * x24
    x48 = x24 * x31
    x49 = x34 + x48
    x50 = x25 + x32
    x51 = x44 * (x0 * (x21 + x38 + 3.0 * x47 + x49) + x24 * x50)
    x52 = -x1 * (a * A[2] + b * B[2])
    x53 = -x52 - B[2]
    x54 = x41 * x7
    x55 = x46**2 * x54
    x56 = x0 * x54
    x57 = x55 + x56
    x58 = x18 * x24
    x59 = x24 * x27
    x60 = x30 + x59
    x61 = x15 + x47
    x62 = x0 * (x19 + x26 + 2.0 * x58 + x60) + x24 * x61
    x63 = x42 * x7
    x64 = x44 * x53
    x65 = x53**2 * x63
    x66 = x0 * x63
    x67 = x65 + x66
    x68 = x46 * x57
    x69 = x46 * x56
    x70 = 2.0 * x69
    x71 = x68 + x70
    x72 = x24 * x9
    x73 = x16 * x24
    x74 = x17 + x58
    x75 = x0 * (x10 + x13 + x72 + x73) + x24 * x74
    x76 = x53 * x63
    x77 = x46 * x54
    x78 = x53 * x67
    x79 = x53 * x66
    x80 = 2.0 * x79
    x81 = x78 + x80
    x82 = 3.0 * x56
    x83 = x0 * (3.0 * x55 + x82)
    x84 = x46 * x71
    x85 = x83 + x84
    x86 = x24 * x8
    x87 = x12 + x73
    x88 = x0 * (x16 + x86) + x24 * x87
    x89 = 3.0 * x66
    x90 = x0 * (3.0 * x65 + x89)
    x91 = x53 * x81
    x92 = x90 + x91
    x93 = -x45 - A[1]
    x94 = x40 * x44
    x95 = x54 * x93
    x96 = x46 * x95
    x97 = x56 + x96
    x98 = x57 * x93
    x99 = x70 + x98
    x100 = x71 * x93
    x101 = x100 + x83
    x102 = 8.0 * x69
    x103 = x0 * (x102 + 4.0 * x68) + x85 * x93
    x104 = -x52 - A[2]
    x105 = x104 * x44
    x106 = x104 * x63
    x107 = x106 * x53
    x108 = x107 + x66
    x109 = x104 * x67
    x110 = x109 + x80
    x111 = x104 * x81
    x112 = x111 + x90
    x113 = 8.0 * x79
    x114 = x0 * (x113 + 4.0 * x78) + x104 * x92
    x115 = x54 * x93**2 + x56
    x116 = x0 * (x77 + x95) + x93 * x97
    x117 = x55 + x82
    x118 = x0 * (x117 + 2.0 * x96) + x93 * x99
    x119 = x0 * (x102 + x68 + 3.0 * x98) + x101 * x93
    x120 = x43 * x6
    x121 = x120 * (x0 * (4.0 * x100 + 5.0 * x83 + x84) + x103 * x93)
    x122 = x120 * x4
    x123 = 3.14159265358979 * x1 * x41 * x6
    x124 = x123 * x4
    x125 = x104**2 * x63 + x66
    x126 = x0 * (x106 + x76) + x104 * x108
    x127 = x65 + x89
    x128 = x0 * (2.0 * x107 + x127) + x104 * x110
    x129 = x0 * (3.0 * x109 + x113 + x78) + x104 * x112
    x130 = x123 * (x0 * (4.0 * x111 + 5.0 * x90 + x91) + x104 * x114)
    x131 = -x45 - C[1]
    x132 = x44 * (x0 * (5.0 * x34 + x35 + 4.0 * x48) + x24 * x37)
    x133 = x131 * x77
    x134 = x133 + x56
    x135 = x0 * (x28 + x33 + 3.0 * x59) + x24 * x49
    x136 = x131 * x54
    x137 = x0 * (x136 + x77)
    x138 = x134 * x46
    x139 = x137 + x138
    x140 = x0 * (x14 + 2.0 * x72) + x24 * x60
    x141 = x0 * (x117 + 2.0 * x133)
    x142 = x139 * x46
    x143 = x141 + x142
    x144 = x12 + x72
    x145 = x0 * (x86 + x9) + x144 * x24
    x146 = 3.0 * x137
    x147 = x0 * (3.0 * x138 + x146 + x71)
    x148 = x143 * x46
    x149 = x147 + x148
    x150 = x12 + x24**2 * x8
    x151 = x131 * x95
    x152 = x151 + x56
    x153 = x134 * x93
    x154 = x137 + x153
    x155 = x139 * x93
    x156 = x141 + x155
    x157 = x143 * x93
    x158 = x147 + x157
    x159 = 4.0 * x141
    x160 = x0 * (4.0 * x142 + x159 + x85) + x149 * x93
    x161 = x120 * x160
    x162 = x120 * x24
    x163 = x123 * x131
    x164 = x0 * (x136 + x95) + x152 * x93
    x165 = x0 * (x133 + x151 + x82 + x96) + x154 * x93
    x166 = x0 * (x138 + x146 + 2.0 * x153 + x99) + x156 * x93
    x167 = x120 * (x0 * (x101 + x142 + 3.0 * x155 + x159) + x158 * x93)
    x168 = x120 * x3
    x169 = -x52 - C[2]
    x170 = x169 * x44
    x171 = x169 * x76
    x172 = x171 + x66
    x173 = x169 * x63
    x174 = x0 * (x173 + x76)
    x175 = x172 * x53
    x176 = x174 + x175
    x177 = x0 * (x127 + 2.0 * x171)
    x178 = x176 * x53
    x179 = x177 + x178
    x180 = 3.0 * x174
    x181 = x0 * (3.0 * x175 + x180 + x81)
    x182 = x179 * x53
    x183 = x181 + x182
    x184 = x123 * x24
    x185 = x106 * x169
    x186 = x185 + x66
    x187 = x104 * x172
    x188 = x174 + x187
    x189 = x104 * x176
    x190 = x177 + x189
    x191 = x104 * x179
    x192 = x181 + x191
    x193 = 4.0 * x177
    x194 = x0 * (4.0 * x178 + x193 + x92) + x104 * x183
    x195 = x123 * x194
    x196 = x123 * x3
    x197 = x0 * (x106 + x173) + x104 * x186
    x198 = x0 * (x107 + x171 + x185 + x89) + x104 * x188
    x199 = x0 * (x110 + x175 + x180 + 2.0 * x187) + x104 * x190
    x200 = x123 * (x0 * (x112 + x178 + 3.0 * x189 + x193) + x104 * x192)

    # 270 item(s)
    return numpy.array(
        [
            x44 * (x0 * (x23 + 4.0 * x25 + 5.0 * x32 + x37) + x24 * x40),
            x46 * x51,
            x51 * x53,
            x57 * x62 * x63,
            x46 * x62 * x64,
            x54 * x62 * x67,
            x63 * x71 * x75,
            x57 * x75 * x76,
            x67 * x75 * x77,
            x54 * x75 * x81,
            x63 * x85 * x88,
            x71 * x76 * x88,
            x57 * x67 * x88,
            x77 * x81 * x88,
            x54 * x88 * x92,
            x93 * x94,
            x50 * x63 * x97,
            x50 * x64 * x93,
            x61 * x63 * x99,
            x61 * x76 * x97,
            x61 * x67 * x95,
            x101 * x63 * x74,
            x74 * x76 * x99,
            x67 * x74 * x97,
            x74 * x81 * x95,
            x103 * x63 * x87,
            x101 * x76 * x87,
            x67 * x87 * x99,
            x81 * x87 * x97,
            x87 * x92 * x95,
            x104 * x94,
            x105 * x46 * x50,
            x108 * x50 * x54,
            x106 * x57 * x61,
            x108 * x61 * x77,
            x110 * x54 * x61,
            x106 * x71 * x74,
            x108 * x57 * x74,
            x110 * x74 * x77,
            x112 * x54 * x74,
            x106 * x85 * x87,
            x108 * x71 * x87,
            x110 * x57 * x87,
            x112 * x77 * x87,
            x114 * x54 * x87,
            x115 * x39 * x63,
            x116 * x22 * x63,
            x115 * x22 * x76,
            x118 * x20 * x63,
            x116 * x20 * x76,
            x115 * x20 * x67,
            x119 * x18 * x63,
            x118 * x18 * x76,
            x116 * x18 * x67,
            x115 * x18 * x81,
            x121 * x4,
            x119 * x122 * x53,
            x118 * x16 * x67,
            x116 * x16 * x81,
            x115 * x16 * x92,
            x105 * x39 * x93,
            x106 * x22 * x97,
            x108 * x22 * x95,
            x106 * x20 * x99,
            x108 * x20 * x97,
            x110 * x20 * x95,
            x101 * x106 * x18,
            x108 * x18 * x99,
            x110 * x18 * x97,
            x112 * x18 * x95,
            x103 * x104 * x122,
            x101 * x108 * x16,
            x110 * x16 * x99,
            x112 * x16 * x97,
            x114 * x124 * x93,
            x125 * x39 * x54,
            x125 * x22 * x77,
            x126 * x22 * x54,
            x125 * x20 * x57,
            x126 * x20 * x77,
            x128 * x20 * x54,
            x125 * x18 * x71,
            x126 * x18 * x57,
            x128 * x18 * x77,
            x129 * x18 * x54,
            x125 * x16 * x85,
            x126 * x16 * x71,
            x128 * x16 * x57,
            x124 * x129 * x46,
            x130 * x4,
            x131 * x132,
            x134 * x135 * x63,
            x131 * x135 * x64,
            x139 * x140 * x63,
            x134 * x140 * x76,
            x136 * x140 * x67,
            x143 * x145 * x63,
            x139 * x145 * x76,
            x134 * x145 * x67,
            x136 * x145 * x81,
            x149 * x150 * x63,
            x143 * x150 * x76,
            x139 * x150 * x67,
            x134 * x150 * x81,
            x136 * x150 * x92,
            x152 * x37 * x63,
            x154 * x49 * x63,
            x152 * x49 * x76,
            x156 * x60 * x63,
            x154 * x60 * x76,
            x152 * x60 * x67,
            x144 * x158 * x63,
            x144 * x156 * x76,
            x144 * x154 * x67,
            x144 * x152 * x81,
            x161 * x24,
            x158 * x162 * x53,
            x156 * x67 * x86,
            x154 * x81 * x86,
            x152 * x86 * x92,
            x105 * x131 * x37,
            x106 * x134 * x49,
            x108 * x136 * x49,
            x106 * x139 * x60,
            x108 * x134 * x60,
            x110 * x136 * x60,
            x106 * x143 * x144,
            x108 * x139 * x144,
            x110 * x134 * x144,
            x112 * x136 * x144,
            x104 * x149 * x162,
            x108 * x143 * x86,
            x110 * x139 * x86,
            x112 * x134 * x86,
            x114 * x163 * x24,
            x164 * x36 * x63,
            x165 * x31 * x63,
            x164 * x31 * x76,
            x166 * x27 * x63,
            x165 * x27 * x76,
            x164 * x27 * x67,
            x167 * x3,
            x166 * x168 * x53,
            x165 * x67 * x9,
            x164 * x81 * x9,
            x120 * (x0 * (x103 + 5.0 * x147 + x148 + 4.0 * x157) + x160 * x93),
            x167 * x53,
            x166 * x67 * x8,
            x165 * x8 * x81,
            x164 * x8 * x92,
            x106 * x152 * x36,
            x106 * x154 * x31,
            x108 * x152 * x31,
            x106 * x156 * x27,
            x108 * x154 * x27,
            x110 * x152 * x27,
            x104 * x158 * x168,
            x108 * x156 * x9,
            x110 * x154 * x9,
            x112 * x152 * x9,
            x104 * x161,
            x108 * x158 * x8,
            x110 * x156 * x8,
            x112 * x154 * x8,
            x114 * x152 * x8,
            x125 * x136 * x36,
            x125 * x134 * x31,
            x126 * x136 * x31,
            x125 * x139 * x27,
            x126 * x134 * x27,
            x128 * x136 * x27,
            x125 * x143 * x9,
            x126 * x139 * x9,
            x128 * x134 * x9,
            x129 * x163 * x3,
            x125 * x149 * x8,
            x126 * x143 * x8,
            x128 * x139 * x8,
            x129 * x134 * x8,
            x130 * x131,
            x132 * x169,
            x135 * x170 * x46,
            x135 * x172 * x54,
            x140 * x173 * x57,
            x140 * x172 * x77,
            x140 * x176 * x54,
            x145 * x173 * x71,
            x145 * x172 * x57,
            x145 * x176 * x77,
            x145 * x179 * x54,
            x150 * x173 * x85,
            x150 * x172 * x71,
            x150 * x176 * x57,
            x150 * x179 * x77,
            x150 * x183 * x54,
            x170 * x37 * x93,
            x173 * x49 * x97,
            x172 * x49 * x95,
            x173 * x60 * x99,
            x172 * x60 * x97,
            x176 * x60 * x95,
            x101 * x144 * x173,
            x144 * x172 * x99,
            x144 * x176 * x97,
            x144 * x179 * x95,
            x103 * x162 * x169,
            x101 * x172 * x86,
            x176 * x86 * x99,
            x179 * x86 * x97,
            x183 * x184 * x93,
            x186 * x37 * x54,
            x186 * x49 * x77,
            x188 * x49 * x54,
            x186 * x57 * x60,
            x188 * x60 * x77,
            x190 * x54 * x60,
            x144 * x186 * x71,
            x144 * x188 * x57,
            x144 * x190 * x77,
            x144 * x192 * x54,
            x186 * x85 * x86,
            x188 * x71 * x86,
            x190 * x57 * x86,
            x184 * x192 * x46,
            x195 * x24,
            x115 * x173 * x36,
            x116 * x173 * x31,
            x115 * x172 * x31,
            x118 * x173 * x27,
            x116 * x172 * x27,
            x115 * x176 * x27,
            x119 * x168 * x169,
            x118 * x172 * x9,
            x116 * x176 * x9,
            x115 * x179 * x9,
            x121 * x169,
            x119 * x172 * x8,
            x118 * x176 * x8,
            x116 * x179 * x8,
            x115 * x183 * x8,
            x186 * x36 * x95,
            x186 * x31 * x97,
            x188 * x31 * x95,
            x186 * x27 * x99,
            x188 * x27 * x97,
            x190 * x27 * x95,
            x101 * x186 * x9,
            x188 * x9 * x99,
            x190 * x9 * x97,
            x192 * x196 * x93,
            x103 * x186 * x8,
            x101 * x188 * x8,
            x190 * x8 * x99,
            x192 * x8 * x97,
            x195 * x93,
            x197 * x36 * x54,
            x197 * x31 * x77,
            x198 * x31 * x54,
            x197 * x27 * x57,
            x198 * x27 * x77,
            x199 * x27 * x54,
            x197 * x71 * x9,
            x198 * x57 * x9,
            x196 * x199 * x46,
            x200 * x3,
            x197 * x8 * x85,
            x198 * x71 * x8,
            x199 * x57 * x8,
            x200 * x46,
            x123 * (x0 * (x114 + 5.0 * x181 + x182 + 4.0 * x191) + x104 * x194),
        ]
    )


def dipole3d_30(a, A, b, B, C):
    """Cartesian 3D (fs) dipole moment integrals.
    The origin is at C.

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
    x10 = -x2 - C[0]
    x11 = x3 * x7
    x12 = x10 * x11
    x13 = x12 + x9
    x14 = x0 * (x10 * x7 + x11) + x13 * x3
    x15 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x16 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x17 = 3.14159265358979 * x1 * x16
    x18 = x15 * x17
    x19 = -x1 * (a * A[1] + b * B[1])
    x20 = -x19 - A[1]
    x21 = x14 * x18
    x22 = -x1 * (a * A[2] + b * B[2])
    x23 = -x22 - A[2]
    x24 = x15 * x6
    x25 = x20**2 * x24
    x26 = x0 * x24
    x27 = x25 + x26
    x28 = x16 * x6
    x29 = x18 * x23
    x30 = x23**2 * x28
    x31 = x0 * x28
    x32 = x30 + x31
    x33 = x17 * x5
    x34 = x33 * (2.0 * x20 * x26 + x20 * x27)
    x35 = x23 * x33
    x36 = 3.14159265358979 * x1 * x15 * x5
    x37 = x32 * x36
    x38 = x36 * (2.0 * x23 * x31 + x23 * x32)
    x39 = -x19 - C[1]
    x40 = x8 + x9
    x41 = x18 * (2.0 * x0 * x11 + x3 * x40)
    x42 = x20 * x24
    x43 = x39 * x42
    x44 = x26 + x43
    x45 = x0 * (x24 * x39 + x42) + x20 * x44
    x46 = x33 * x45
    x47 = -x22 - C[2]
    x48 = x23 * x28
    x49 = x47 * x48
    x50 = x31 + x49
    x51 = x0 * (x28 * x47 + x48) + x23 * x50
    x52 = x36 * x51

    # 30 item(s)
    return numpy.array(
        [
            x18 * (x0 * (2.0 * x12 + x8 + 3.0 * x9) + x14 * x3),
            x20 * x21,
            x21 * x23,
            x13 * x27 * x28,
            x13 * x20 * x29,
            x13 * x24 * x32,
            x10 * x34,
            x10 * x27 * x35,
            x10 * x20 * x37,
            x10 * x38,
            x39 * x41,
            x28 * x40 * x44,
            x29 * x39 * x40,
            x3 * x46,
            x3 * x35 * x44,
            x3 * x37 * x39,
            x33 * (x0 * (x25 + 3.0 * x26 + 2.0 * x43) + x20 * x45),
            x23 * x46,
            x32 * x44 * x7,
            x38 * x39,
            x41 * x47,
            x18 * x20 * x40 * x47,
            x24 * x40 * x50,
            x27 * x3 * x33 * x47,
            x20 * x3 * x36 * x50,
            x3 * x52,
            x34 * x47,
            x27 * x50 * x7,
            x20 * x52,
            x36 * (x0 * (x30 + 3.0 * x31 + 2.0 * x49) + x23 * x51),
        ]
    )


def dipole3d_31(a, A, b, B, C):
    """Cartesian 3D (fp) dipole moment integrals.
    The origin is at C.

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
    x12 = -x2 - A[0]
    x13 = x0 * x7
    x14 = x8 * x9
    x15 = x13 + x14
    x16 = x12 * x15
    x17 = x12 * x5
    x18 = x17 * x6
    x19 = x18 * x9
    x20 = x13 + x19
    x21 = x0 * (x10 + x18) + x12 * x20
    x22 = x18 * x3
    x23 = x13 + x22
    x24 = x0 * (x18 + x8) + x12 * x23
    x25 = 3.0 * x13
    x26 = x11 + x16
    x27 = x0 * (x14 + x19 + x22 + x25) + x12 * x26
    x28 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x29 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x30 = 3.14159265358979 * x1 * x29
    x31 = x28 * x30
    x32 = -x1 * (a * A[1] + b * B[1])
    x33 = -x32 - B[1]
    x34 = x12**2 * x7
    x35 = x25 + x34
    x36 = x31 * (x0 * (2.0 * x19 + x35) + x12 * x21)
    x37 = -x1 * (a * A[2] + b * B[2])
    x38 = -x37 - B[2]
    x39 = -x32 - A[1]
    x40 = x27 * x31
    x41 = x0 * x6
    x42 = x28 * x41
    x43 = x28 * x6
    x44 = x39 * x43
    x45 = x33 * x44
    x46 = x42 + x45
    x47 = x29 * x6
    x48 = x21 * x31
    x49 = -x37 - A[2]
    x50 = x29 * x41
    x51 = x47 * x49
    x52 = x38 * x51
    x53 = x50 + x52
    x54 = x39**2 * x43
    x55 = x42 + x54
    x56 = x33 * x43
    x57 = x0 * (x44 + x56) + x39 * x46
    x58 = x38 * x47
    x59 = x31 * x49
    x60 = x47 * x49**2
    x61 = x50 + x60
    x62 = x0 * (x51 + x58) + x49 * x53
    x63 = 2.0 * x39 * x42 + x39 * x55
    x64 = 3.0 * x42
    x65 = x54 + x64
    x66 = x30 * x5
    x67 = x66 * (x0 * (2.0 * x45 + x65) + x39 * x57)
    x68 = x63 * x66
    x69 = x49 * x66
    x70 = 3.14159265358979 * x1 * x28
    x71 = x5 * x70
    x72 = x39 * x71
    x73 = 2.0 * x49 * x50 + x49 * x61
    x74 = x71 * x73
    x75 = 3.0 * x50
    x76 = x60 + x75
    x77 = x71 * (x0 * (2.0 * x52 + x76) + x49 * x62)
    x78 = -x32 - C[1]
    x79 = x31 * (x0 * (2.0 * x22 + x35) + x12 * x24)
    x80 = x56 * x78
    x81 = x42 + x80
    x82 = x13 + x34
    x83 = 2.0 * x0 * x18 + x12 * x82
    x84 = x31 * x83
    x85 = x44 * x78
    x86 = x42 + x85
    x87 = x43 * x78
    x88 = x0 * (x56 + x87)
    x89 = x39 * x81
    x90 = x88 + x89
    x91 = x0 * (x44 + x87) + x39 * x86
    x92 = x0 * (x45 + x64 + x80 + x85) + x39 * x90
    x93 = x17 * x30
    x94 = x17 * x70
    x95 = x66 * (x0 * (x65 + 2.0 * x85) + x39 * x91)
    x96 = -x37 - C[2]
    x97 = x58 * x96
    x98 = x50 + x97
    x99 = x47 * x96
    x100 = x51 * x96
    x101 = x100 + x50
    x102 = x0 * (x58 + x99)
    x103 = x49 * x98
    x104 = x102 + x103
    x105 = x0 * (x51 + x99) + x101 * x49
    x106 = x0 * (x100 + x52 + x75 + x97) + x104 * x49
    x107 = x71 * (x0 * (2.0 * x100 + x76) + x105 * x49)

    # 90 item(s)
    return numpy.array(
        [
            x31 * (x0 * (2.0 * x11 + 2.0 * x16 + x21 + x24) + x12 * x27),
            x33 * x36,
            x36 * x38,
            x39 * x40,
            x21 * x46 * x47,
            x38 * x39 * x48,
            x40 * x49,
            x33 * x48 * x49,
            x21 * x43 * x53,
            x26 * x47 * x55,
            x20 * x47 * x57,
            x20 * x55 * x58,
            x26 * x39 * x59,
            x20 * x46 * x51,
            x20 * x44 * x53,
            x26 * x43 * x61,
            x20 * x56 * x61,
            x20 * x43 * x62,
            x15 * x47 * x63,
            x67 * x9,
            x38 * x68 * x9,
            x15 * x51 * x55,
            x57 * x69 * x9,
            x10 * x53 * x55,
            x15 * x44 * x61,
            x10 * x46 * x61,
            x62 * x72 * x9,
            x15 * x43 * x73,
            x33 * x74 * x9,
            x77 * x9,
            x78 * x79,
            x47 * x81 * x83,
            x38 * x78 * x84,
            x24 * x47 * x86,
            x47 * x82 * x90,
            x58 * x82 * x86,
            x24 * x59 * x78,
            x51 * x81 * x82,
            x53 * x82 * x87,
            x23 * x47 * x91,
            x92 * x93,
            x38 * x91 * x93,
            x23 * x51 * x86,
            x49 * x90 * x93,
            x18 * x53 * x86,
            x23 * x61 * x87,
            x18 * x61 * x81,
            x62 * x78 * x94,
            x3 * x95,
            x66 * (x0 * (x57 + 2.0 * x88 + 2.0 * x89 + x91) + x39 * x92),
            x38 * x95,
            x3 * x69 * x91,
            x69 * x92,
            x53 * x7 * x91,
            x61 * x8 * x86,
            x61 * x7 * x90,
            x62 * x7 * x86,
            x3 * x74 * x78,
            x7 * x73 * x81,
            x77 * x78,
            x79 * x96,
            x33 * x84 * x96,
            x43 * x83 * x98,
            x24 * x31 * x39 * x96,
            x46 * x82 * x99,
            x44 * x82 * x98,
            x101 * x24 * x43,
            x101 * x56 * x82,
            x104 * x43 * x82,
            x23 * x55 * x99,
            x57 * x93 * x96,
            x18 * x55 * x98,
            x101 * x23 * x44,
            x101 * x18 * x46,
            x104 * x39 * x94,
            x105 * x23 * x43,
            x105 * x33 * x94,
            x106 * x94,
            x3 * x68 * x96,
            x67 * x96,
            x63 * x7 * x98,
            x101 * x55 * x8,
            x101 * x57 * x7,
            x104 * x55 * x7,
            x105 * x3 * x72,
            x105 * x46 * x7,
            x106 * x72,
            x107 * x3,
            x107 * x33,
            x71 * (x0 * (2.0 * x102 + 2.0 * x103 + x105 + x62) + x106 * x49),
        ]
    )


def dipole3d_32(a, A, b, B, C):
    """Cartesian 3D (fd) dipole moment integrals.
    The origin is at C.

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
    x11 = x4**2 * x8
    x12 = x0 * x8
    x13 = 3.0 * x12
    x14 = x11 + x13
    x15 = x0 * (2.0 * x10 + x14)
    x16 = -x2 - A[0]
    x17 = x16 * x9
    x18 = x3 * x8
    x19 = x16 * x18
    x20 = x0 * (x10 + x13 + x17 + x19)
    x21 = x0 * (x18 + x9)
    x22 = x10 + x12
    x23 = x16 * x22
    x24 = x21 + x23
    x25 = x16 * x24
    x26 = x22 * x4
    x27 = x21 + x26
    x28 = x16 * x27
    x29 = 2.0 * x17
    x30 = x11 + x12
    x31 = x16 * x30
    x32 = x12 * x4
    x33 = x31 + 2.0 * x32
    x34 = x0 * (x14 + x29) + x16 * x33
    x35 = 2.0 * x23
    x36 = x15 + x28
    x37 = x0 * (3.0 * x21 + x26 + x33 + x35) + x16 * x36
    x38 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x39 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x40 = 3.14159265358979 * x1 * x39
    x41 = x38 * x40
    x42 = -x1 * (a * A[1] + b * B[1])
    x43 = -x42 - B[1]
    x44 = x16 * x8
    x45 = x12 + x19
    x46 = x0 * (x18 + x44) + x16 * x45
    x47 = x0 * (x44 + x9)
    x48 = x12 + x17
    x49 = x16 * x48
    x50 = x47 + x49
    x51 = x20 + x25
    x52 = x41 * (x0 * (2.0 * x21 + x35 + x46 + x50) + x16 * x51)
    x53 = -x1 * (a * A[2] + b * B[2])
    x54 = -x53 - B[2]
    x55 = x38 * x7
    x56 = x43**2 * x55
    x57 = x0 * x55
    x58 = x56 + x57
    x59 = x16**2 * x8
    x60 = x13 + x59
    x61 = x0 * (2.0 * x19 + x60) + x16 * x46
    x62 = x39 * x7
    x63 = x41 * x54
    x64 = x54**2 * x62
    x65 = x0 * x62
    x66 = x64 + x65
    x67 = -x42 - A[1]
    x68 = x37 * x41
    x69 = x55 * x67
    x70 = x43 * x69
    x71 = x57 + x70
    x72 = x58 * x67
    x73 = 2.0 * x57
    x74 = x43 * x73 + x72
    x75 = x54 * x62
    x76 = -x53 - A[2]
    x77 = x41 * x76
    x78 = x62 * x76
    x79 = x54 * x78
    x80 = x65 + x79
    x81 = x43 * x55
    x82 = x66 * x76
    x83 = 2.0 * x65
    x84 = x54 * x83 + x82
    x85 = x55 * x67**2
    x86 = x57 + x85
    x87 = x0 * (x69 + x81)
    x88 = x67 * x71
    x89 = x87 + x88
    x90 = 3.0 * x57
    x91 = 2.0 * x70 + x90
    x92 = x0 * (x56 + x91) + x67 * x74
    x93 = x62 * x76**2
    x94 = x65 + x93
    x95 = x0 * (x75 + x78)
    x96 = x76 * x80
    x97 = x95 + x96
    x98 = 3.0 * x65
    x99 = 2.0 * x79 + x98
    x100 = x0 * (x64 + x99) + x76 * x84
    x101 = x67 * x73 + x67 * x86
    x102 = x0 * (x85 + x91) + x67 * x89
    x103 = x40 * x6
    x104 = x103 * (x0 * (4.0 * x43 * x57 + 2.0 * x72 + 2.0 * x87 + 2.0 * x88) + x67 * x92)
    x105 = x103 * x3
    x106 = 3.14159265358979 * x1 * x38 * x6
    x107 = x106 * x3
    x108 = x76 * x83 + x76 * x94
    x109 = x0 * (x93 + x99) + x76 * x97
    x110 = x106 * (
        x0 * (4.0 * x54 * x65 + 2.0 * x82 + 2.0 * x95 + 2.0 * x96) + x100 * x76
    )
    x111 = -x42 - C[1]
    x112 = x41 * (x0 * (2.0 * x31 + 4.0 * x32 + 2.0 * x47 + 2.0 * x49) + x16 * x34)
    x113 = x111 * x81
    x114 = x113 + x57
    x115 = x0 * (x29 + x60) + x16 * x50
    x116 = x111 * x55
    x117 = x0 * (x116 + x81)
    x118 = x114 * x43
    x119 = x117 + x118
    x120 = x12 + x59
    x121 = 2.0 * x12 * x16 + x120 * x16
    x122 = x111 * x69
    x123 = x122 + x57
    x124 = x114 * x67
    x125 = x117 + x124
    x126 = x0 * (2.0 * x113 + x56 + x90)
    x127 = x119 * x67
    x128 = x126 + x127
    x129 = x0 * (x116 + x69) + x123 * x67
    x130 = x0 * (x113 + x122 + x70 + x90)
    x131 = x125 * x67
    x132 = x130 + x131
    x133 = 2.0 * x124
    x134 = x0 * (3.0 * x117 + x118 + x133 + x74) + x128 * x67
    x135 = x103 * x134
    x136 = x103 * x16
    x137 = x106 * x111
    x138 = x0 * (2.0 * x122 + x85 + x90) + x129 * x67
    x139 = x103 * (x0 * (2.0 * x117 + x129 + x133 + x89) + x132 * x67)
    x140 = x103 * x4
    x141 = -x53 - C[2]
    x142 = x141 * x41
    x143 = x141 * x75
    x144 = x143 + x65
    x145 = x141 * x62
    x146 = x0 * (x145 + x75)
    x147 = x144 * x54
    x148 = x146 + x147
    x149 = x141 * x78
    x150 = x149 + x65
    x151 = x144 * x76
    x152 = x146 + x151
    x153 = x0 * (2.0 * x143 + x64 + x98)
    x154 = x148 * x76
    x155 = x153 + x154
    x156 = x106 * x16
    x157 = x0 * (x145 + x78) + x150 * x76
    x158 = x0 * (x143 + x149 + x79 + x98)
    x159 = x152 * x76
    x160 = x158 + x159
    x161 = 2.0 * x151
    x162 = x0 * (3.0 * x146 + x147 + x161 + x84) + x155 * x76
    x163 = x106 * x162
    x164 = x106 * x4
    x165 = x0 * (2.0 * x149 + x93 + x98) + x157 * x76
    x166 = x106 * (x0 * (2.0 * x146 + x157 + x161 + x97) + x160 * x76)

    # 180 item(s)
    return numpy.array(
        [
            x41
            * (x0 * (2.0 * x15 + 2.0 * x20 + 2.0 * x25 + 2.0 * x28 + x34) + x16 * x37),
            x43 * x52,
            x52 * x54,
            x58 * x61 * x62,
            x43 * x61 * x63,
            x55 * x61 * x66,
            x67 * x68,
            x51 * x62 * x71,
            x51 * x63 * x67,
            x46 * x62 * x74,
            x46 * x71 * x75,
            x46 * x66 * x69,
            x68 * x76,
            x43 * x51 * x77,
            x51 * x55 * x80,
            x46 * x58 * x78,
            x46 * x80 * x81,
            x46 * x55 * x84,
            x36 * x62 * x86,
            x24 * x62 * x89,
            x24 * x75 * x86,
            x45 * x62 * x92,
            x45 * x75 * x89,
            x45 * x66 * x86,
            x36 * x67 * x77,
            x24 * x71 * x78,
            x24 * x69 * x80,
            x45 * x74 * x78,
            x45 * x71 * x80,
            x45 * x69 * x84,
            x36 * x55 * x94,
            x24 * x81 * x94,
            x24 * x55 * x97,
            x45 * x58 * x94,
            x45 * x81 * x97,
            x100 * x45 * x55,
            x101 * x27 * x62,
            x102 * x22 * x62,
            x101 * x22 * x75,
            x104 * x3,
            x102 * x105 * x54,
            x101 * x18 * x66,
            x27 * x78 * x86,
            x22 * x78 * x89,
            x22 * x80 * x86,
            x105 * x76 * x92,
            x18 * x80 * x89,
            x18 * x84 * x86,
            x27 * x69 * x94,
            x22 * x71 * x94,
            x22 * x69 * x97,
            x18 * x74 * x94,
            x18 * x71 * x97,
            x100 * x107 * x67,
            x108 * x27 * x55,
            x108 * x22 * x81,
            x109 * x22 * x55,
            x108 * x18 * x58,
            x107 * x109 * x43,
            x110 * x3,
            x111 * x112,
            x114 * x115 * x62,
            x111 * x115 * x63,
            x119 * x121 * x62,
            x114 * x121 * x75,
            x116 * x121 * x66,
            x123 * x34 * x62,
            x125 * x50 * x62,
            x123 * x50 * x75,
            x120 * x128 * x62,
            x120 * x125 * x75,
            x120 * x123 * x66,
            x111 * x34 * x77,
            x114 * x50 * x78,
            x116 * x50 * x80,
            x119 * x120 * x78,
            x114 * x120 * x80,
            x116 * x120 * x84,
            x129 * x33 * x62,
            x132 * x48 * x62,
            x129 * x48 * x75,
            x135 * x16,
            x132 * x136 * x54,
            x129 * x44 * x66,
            x123 * x33 * x78,
            x125 * x48 * x78,
            x123 * x48 * x80,
            x128 * x136 * x76,
            x125 * x44 * x80,
            x123 * x44 * x84,
            x116 * x33 * x94,
            x114 * x48 * x94,
            x116 * x48 * x97,
            x119 * x44 * x94,
            x114 * x44 * x97,
            x100 * x137 * x16,
            x138 * x30 * x62,
            x139 * x4,
            x138 * x140 * x54,
            x103
            * (
                x0 * (2.0 * x126 + 2.0 * x127 + 2.0 * x130 + 2.0 * x131 + x92)
                + x134 * x67
            ),
            x139 * x54,
            x138 * x66 * x8,
            x129 * x30 * x78,
            x132 * x140 * x76,
            x129 * x80 * x9,
            x135 * x76,
            x132 * x8 * x80,
            x129 * x8 * x84,
            x123 * x30 * x94,
            x125 * x9 * x94,
            x123 * x9 * x97,
            x128 * x8 * x94,
            x125 * x8 * x97,
            x100 * x123 * x8,
            x108 * x116 * x30,
            x108 * x114 * x9,
            x109 * x137 * x4,
            x108 * x119 * x8,
            x109 * x114 * x8,
            x110 * x111,
            x112 * x141,
            x115 * x142 * x43,
            x115 * x144 * x55,
            x121 * x145 * x58,
            x121 * x144 * x81,
            x121 * x148 * x55,
            x142 * x34 * x67,
            x145 * x50 * x71,
            x144 * x50 * x69,
            x120 * x145 * x74,
            x120 * x144 * x71,
            x120 * x148 * x69,
            x150 * x34 * x55,
            x150 * x50 * x81,
            x152 * x50 * x55,
            x120 * x150 * x58,
            x120 * x152 * x81,
            x120 * x155 * x55,
            x145 * x33 * x86,
            x145 * x48 * x89,
            x144 * x48 * x86,
            x136 * x141 * x92,
            x144 * x44 * x89,
            x148 * x44 * x86,
            x150 * x33 * x69,
            x150 * x48 * x71,
            x152 * x48 * x69,
            x150 * x44 * x74,
            x152 * x44 * x71,
            x155 * x156 * x67,
            x157 * x33 * x55,
            x157 * x48 * x81,
            x160 * x48 * x55,
            x157 * x44 * x58,
            x156 * x160 * x43,
            x16 * x163,
            x101 * x145 * x30,
            x102 * x140 * x141,
            x101 * x144 * x9,
            x104 * x141,
            x102 * x144 * x8,
            x101 * x148 * x8,
            x150 * x30 * x86,
            x150 * x89 * x9,
            x152 * x86 * x9,
            x150 * x8 * x92,
            x152 * x8 * x89,
            x155 * x8 * x86,
            x157 * x30 * x69,
            x157 * x71 * x9,
            x160 * x164 * x67,
            x157 * x74 * x8,
            x160 * x71 * x8,
            x163 * x67,
            x165 * x30 * x55,
            x164 * x165 * x43,
            x166 * x4,
            x165 * x58 * x8,
            x166 * x43,
            x106
            * (
                x0 * (x100 + 2.0 * x153 + 2.0 * x154 + 2.0 * x158 + 2.0 * x159)
                + x162 * x76
            ),
        ]
    )


def dipole3d_33(a, A, b, B, C):
    """Cartesian 3D (ff) dipole moment integrals.
    The origin is at C.

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
    x12 = 3.0 * x11
    x13 = x0 * x7
    x14 = x8 * x9
    x15 = x13 + x14
    x16 = x15 * x3
    x17 = x3**2 * x7
    x18 = x13 + x17
    x19 = x18 * x3
    x20 = x13 * x3
    x21 = 2.0 * x20
    x22 = x19 + x21
    x23 = x0 * (x12 + 3.0 * x16 + x22)
    x24 = -x2 - A[0]
    x25 = 3.0 * x13
    x26 = x17 + x25
    x27 = x0 * (2.0 * x14 + x26)
    x28 = x11 + x16
    x29 = x28 * x3
    x30 = x27 + x29
    x31 = x24 * x30
    x32 = x15 * x24
    x33 = 2.0 * x32
    x34 = x18 * x24
    x35 = x21 + x34
    x36 = x0 * (x12 + x16 + x33 + x35)
    x37 = x24 * x28
    x38 = x27 + x37
    x39 = x24 * x38
    x40 = x0 * (3.0 * x17 + x25)
    x41 = x22 * x24
    x42 = x40 + x41
    x43 = x0 * (x19 + 8.0 * x20 + 3.0 * x34) + x24 * x42
    x44 = x23 + x31
    x45 = x0 * (4.0 * x27 + x29 + 3.0 * x37 + x42) + x24 * x44
    x46 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x47 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x48 = 3.14159265358979 * x1 * x47
    x49 = x46 * x48
    x50 = -x1 * (a * A[1] + b * B[1])
    x51 = -x50 - B[1]
    x52 = x24 * x8
    x53 = x10 * x24
    x54 = x0 * (x14 + x25 + x52 + x53)
    x55 = x11 + x32
    x56 = x24 * x55
    x57 = 2.0 * x52
    x58 = x0 * (x26 + x57)
    x59 = x24 * x35
    x60 = x58 + x59
    x61 = x36 + x39
    x62 = x49 * (x0 * (2.0 * x27 + 2.0 * x37 + 2.0 * x54 + 2.0 * x56 + x60) + x24 * x61)
    x63 = -x1 * (a * A[2] + b * B[2])
    x64 = -x63 - B[2]
    x65 = x46 * x6
    x66 = x51**2 * x65
    x67 = x0 * x65
    x68 = x66 + x67
    x69 = x24 * x7
    x70 = x13 + x53
    x71 = x0 * (x10 + x69) + x24 * x70
    x72 = x0 * (x69 + x8)
    x73 = x13 + x52
    x74 = x24 * x73
    x75 = x72 + x74
    x76 = x54 + x56
    x77 = x0 * (2.0 * x11 + x33 + x71 + x75) + x24 * x76
    x78 = x47 * x6
    x79 = x49 * x64
    x80 = x64**2 * x78
    x81 = x0 * x78
    x82 = x80 + x81
    x83 = x51 * x68
    x84 = x51 * x67
    x85 = 2.0 * x84
    x86 = x83 + x85
    x87 = x24**2 * x7
    x88 = x25 + x87
    x89 = x0 * (2.0 * x53 + x88) + x24 * x71
    x90 = x64 * x78
    x91 = x51 * x65
    x92 = x64 * x82
    x93 = x64 * x81
    x94 = 2.0 * x93
    x95 = x92 + x94
    x96 = -x50 - A[1]
    x97 = x45 * x49
    x98 = x65 * x96
    x99 = x51 * x98
    x100 = x67 + x99
    x101 = x68 * x96
    x102 = x101 + x85
    x103 = 3.0 * x67
    x104 = x0 * (x103 + 3.0 * x66)
    x105 = x86 * x96
    x106 = x104 + x105
    x107 = -x63 - A[2]
    x108 = x107 * x49
    x109 = x107 * x78
    x110 = x109 * x64
    x111 = x110 + x81
    x112 = x107 * x82
    x113 = x112 + x94
    x114 = 3.0 * x81
    x115 = x0 * (x114 + 3.0 * x80)
    x116 = x107 * x95
    x117 = x115 + x116
    x118 = x65 * x96**2
    x119 = x118 + x67
    x120 = x0 * (x91 + x98)
    x121 = x100 * x96
    x122 = x120 + x121
    x123 = x103 + 2.0 * x99
    x124 = x0 * (x123 + x66)
    x125 = x102 * x96
    x126 = x124 + x125
    x127 = x0 * (3.0 * x101 + x83 + 8.0 * x84) + x106 * x96
    x128 = x107**2 * x78
    x129 = x128 + x81
    x130 = x0 * (x109 + x90)
    x131 = x107 * x111
    x132 = x130 + x131
    x133 = 2.0 * x110 + x114
    x134 = x0 * (x133 + x80)
    x135 = x107 * x113
    x136 = x134 + x135
    x137 = x0 * (3.0 * x112 + x92 + 8.0 * x93) + x107 * x117
    x138 = x119 * x96 + 2.0 * x67 * x96
    x139 = x0 * (x118 + x123) + x122 * x96
    x140 = x0 * (2.0 * x101 + 2.0 * x120 + 2.0 * x121 + 4.0 * x84) + x126 * x96
    x141 = x48 * x5
    x142 = x141 * (x0 * (2.0 * x104 + 2.0 * x105 + 3.0 * x124 + 3.0 * x125) + x127 * x96)
    x143 = x141 * x9
    x144 = 3.14159265358979 * x1 * x46 * x5
    x145 = x144 * x9
    x146 = x107 * x129 + 2.0 * x107 * x81
    x147 = x0 * (x128 + x133) + x107 * x132
    x148 = x0 * (2.0 * x112 + 2.0 * x130 + 2.0 * x131 + 4.0 * x93) + x107 * x136
    x149 = x144 * (x0 * (2.0 * x115 + 2.0 * x116 + 3.0 * x134 + 3.0 * x135) + x107 * x137)
    x150 = -x50 - C[1]
    x151 = x49 * (x0 * (2.0 * x40 + 2.0 * x41 + 3.0 * x58 + 3.0 * x59) + x24 * x43)
    x152 = x150 * x91
    x153 = x152 + x67
    x154 = x0 * (4.0 * x20 + 2.0 * x34 + 2.0 * x72 + 2.0 * x74) + x24 * x60
    x155 = x150 * x65
    x156 = x0 * (x155 + x91)
    x157 = x153 * x51
    x158 = x156 + x157
    x159 = x0 * (x57 + x88) + x24 * x75
    x160 = x0 * (x103 + 2.0 * x152 + x66)
    x161 = x158 * x51
    x162 = x160 + x161
    x163 = x13 + x87
    x164 = 2.0 * x13 * x24 + x163 * x24
    x165 = x150 * x98
    x166 = x165 + x67
    x167 = x153 * x96
    x168 = x156 + x167
    x169 = x158 * x96
    x170 = x160 + x169
    x171 = 3.0 * x156
    x172 = x0 * (3.0 * x157 + x171 + x86)
    x173 = x162 * x96
    x174 = x172 + x173
    x175 = x0 * (x155 + x98) + x166 * x96
    x176 = x0 * (x103 + x152 + x165 + x99)
    x177 = x168 * x96
    x178 = x176 + x177
    x179 = 2.0 * x167
    x180 = x0 * (x102 + x157 + x171 + x179)
    x181 = x170 * x96
    x182 = x180 + x181
    x183 = x0 * (x106 + 4.0 * x160 + x161 + 3.0 * x169) + x174 * x96
    x184 = x141 * x183
    x185 = x141 * x24
    x186 = x144 * x150
    x187 = x0 * (x103 + x118 + 2.0 * x165) + x175 * x96
    x188 = x0 * (x122 + 2.0 * x156 + x175 + x179) + x178 * x96
    x189 = x141 * (
        x0 * (x126 + 2.0 * x160 + 2.0 * x169 + 2.0 * x176 + 2.0 * x177) + x182 * x96
    )
    x190 = x141 * x3
    x191 = -x63 - C[2]
    x192 = x191 * x49
    x193 = x191 * x90
    x194 = x193 + x81
    x195 = x191 * x78
    x196 = x0 * (x195 + x90)
    x197 = x194 * x64
    x198 = x196 + x197
    x199 = x0 * (x114 + 2.0 * x193 + x80)
    x200 = x198 * x64
    x201 = x199 + x200
    x202 = x109 * x191
    x203 = x202 + x81
    x204 = x107 * x194
    x205 = x196 + x204
    x206 = x107 * x198
    x207 = x199 + x206
    x208 = 3.0 * x196
    x209 = x0 * (3.0 * x197 + x208 + x95)
    x210 = x107 * x201
    x211 = x209 + x210
    x212 = x144 * x24
    x213 = x0 * (x109 + x195) + x107 * x203
    x214 = x0 * (x110 + x114 + x193 + x202)
    x215 = x107 * x205
    x216 = x214 + x215
    x217 = 2.0 * x204
    x218 = x0 * (x113 + x197 + x208 + x217)
    x219 = x107 * x207
    x220 = x218 + x219
    x221 = x0 * (x117 + 4.0 * x199 + x200 + 3.0 * x206) + x107 * x211
    x222 = x144 * x221
    x223 = x144 * x3
    x224 = x0 * (x114 + x128 + 2.0 * x202) + x107 * x213
    x225 = x0 * (x132 + 2.0 * x196 + x213 + x217) + x107 * x216
    x226 = x144 * (
        x0 * (x136 + 2.0 * x199 + 2.0 * x206 + 2.0 * x214 + 2.0 * x215) + x107 * x220
    )

    # 300 item(s)
    return numpy.array(
        [
            x49
            * (x0 * (2.0 * x23 + 2.0 * x31 + 3.0 * x36 + 3.0 * x39 + x43) + x24 * x45),
            x51 * x62,
            x62 * x64,
            x68 * x77 * x78,
            x51 * x77 * x79,
            x65 * x77 * x82,
            x78 * x86 * x89,
            x68 * x89 * x90,
            x82 * x89 * x91,
            x65 * x89 * x95,
            x96 * x97,
            x100 * x61 * x78,
            x61 * x79 * x96,
            x102 * x76 * x78,
            x100 * x76 * x90,
            x76 * x82 * x98,
            x106 * x71 * x78,
            x102 * x71 * x90,
            x100 * x71 * x82,
            x71 * x95 * x98,
            x107 * x97,
            x108 * x51 * x61,
            x111 * x61 * x65,
            x109 * x68 * x76,
            x111 * x76 * x91,
            x113 * x65 * x76,
            x109 * x71 * x86,
            x111 * x68 * x71,
            x113 * x71 * x91,
            x117 * x65 * x71,
            x119 * x44 * x78,
            x122 * x38 * x78,
            x119 * x38 * x90,
            x126 * x55 * x78,
            x122 * x55 * x90,
            x119 * x55 * x82,
            x127 * x70 * x78,
            x126 * x70 * x90,
            x122 * x70 * x82,
            x119 * x70 * x95,
            x108 * x44 * x96,
            x100 * x109 * x38,
            x111 * x38 * x98,
            x102 * x109 * x55,
            x100 * x111 * x55,
            x113 * x55 * x98,
            x106 * x109 * x70,
            x102 * x111 * x70,
            x100 * x113 * x70,
            x117 * x70 * x98,
            x129 * x44 * x65,
            x129 * x38 * x91,
            x132 * x38 * x65,
            x129 * x55 * x68,
            x132 * x55 * x91,
            x136 * x55 * x65,
            x129 * x70 * x86,
            x132 * x68 * x70,
            x136 * x70 * x91,
            x137 * x65 * x70,
            x138 * x30 * x78,
            x139 * x28 * x78,
            x138 * x28 * x90,
            x140 * x15 * x78,
            x139 * x15 * x90,
            x138 * x15 * x82,
            x142 * x9,
            x140 * x143 * x64,
            x10 * x139 * x82,
            x10 * x138 * x95,
            x109 * x119 * x30,
            x109 * x122 * x28,
            x111 * x119 * x28,
            x109 * x126 * x15,
            x111 * x122 * x15,
            x113 * x119 * x15,
            x107 * x127 * x143,
            x10 * x111 * x126,
            x10 * x113 * x122,
            x10 * x117 * x119,
            x129 * x30 * x98,
            x100 * x129 * x28,
            x132 * x28 * x98,
            x102 * x129 * x15,
            x100 * x132 * x15,
            x136 * x15 * x98,
            x10 * x106 * x129,
            x10 * x102 * x132,
            x10 * x100 * x136,
            x137 * x145 * x96,
            x146 * x30 * x65,
            x146 * x28 * x91,
            x147 * x28 * x65,
            x146 * x15 * x68,
            x147 * x15 * x91,
            x148 * x15 * x65,
            x10 * x146 * x86,
            x10 * x147 * x68,
            x145 * x148 * x51,
            x149 * x9,
            x150 * x151,
            x153 * x154 * x78,
            x150 * x154 * x79,
            x158 * x159 * x78,
            x153 * x159 * x90,
            x155 * x159 * x82,
            x162 * x164 * x78,
            x158 * x164 * x90,
            x153 * x164 * x82,
            x155 * x164 * x95,
            x166 * x43 * x78,
            x168 * x60 * x78,
            x166 * x60 * x90,
            x170 * x75 * x78,
            x168 * x75 * x90,
            x166 * x75 * x82,
            x163 * x174 * x78,
            x163 * x170 * x90,
            x163 * x168 * x82,
            x163 * x166 * x95,
            x108 * x150 * x43,
            x109 * x153 * x60,
            x111 * x155 * x60,
            x109 * x158 * x75,
            x111 * x153 * x75,
            x113 * x155 * x75,
            x109 * x162 * x163,
            x111 * x158 * x163,
            x113 * x153 * x163,
            x117 * x155 * x163,
            x175 * x42 * x78,
            x178 * x35 * x78,
            x175 * x35 * x90,
            x182 * x73 * x78,
            x178 * x73 * x90,
            x175 * x73 * x82,
            x184 * x24,
            x182 * x185 * x64,
            x178 * x69 * x82,
            x175 * x69 * x95,
            x109 * x166 * x42,
            x109 * x168 * x35,
            x111 * x166 * x35,
            x109 * x170 * x73,
            x111 * x168 * x73,
            x113 * x166 * x73,
            x107 * x174 * x185,
            x111 * x170 * x69,
            x113 * x168 * x69,
            x117 * x166 * x69,
            x129 * x155 * x42,
            x129 * x153 * x35,
            x132 * x155 * x35,
            x129 * x158 * x73,
            x132 * x153 * x73,
            x136 * x155 * x73,
            x129 * x162 * x69,
            x132 * x158 * x69,
            x136 * x153 * x69,
            x137 * x186 * x24,
            x187 * x22 * x78,
            x18 * x188 * x78,
            x18 * x187 * x90,
            x189 * x3,
            x188 * x190 * x64,
            x187 * x8 * x82,
            x141
            * (
                x0 * (x127 + 2.0 * x172 + 2.0 * x173 + 3.0 * x180 + 3.0 * x181)
                + x183 * x96
            ),
            x189 * x64,
            x188 * x7 * x82,
            x187 * x7 * x95,
            x109 * x175 * x22,
            x109 * x178 * x18,
            x111 * x175 * x18,
            x107 * x182 * x190,
            x111 * x178 * x8,
            x113 * x175 * x8,
            x107 * x184,
            x111 * x182 * x7,
            x113 * x178 * x7,
            x117 * x175 * x7,
            x129 * x166 * x22,
            x129 * x168 * x18,
            x132 * x166 * x18,
            x129 * x170 * x8,
            x132 * x168 * x8,
            x136 * x166 * x8,
            x129 * x174 * x7,
            x132 * x170 * x7,
            x136 * x168 * x7,
            x137 * x166 * x7,
            x146 * x155 * x22,
            x146 * x153 * x18,
            x147 * x155 * x18,
            x146 * x158 * x8,
            x147 * x153 * x8,
            x148 * x186 * x3,
            x146 * x162 * x7,
            x147 * x158 * x7,
            x148 * x153 * x7,
            x149 * x150,
            x151 * x191,
            x154 * x192 * x51,
            x154 * x194 * x65,
            x159 * x195 * x68,
            x159 * x194 * x91,
            x159 * x198 * x65,
            x164 * x195 * x86,
            x164 * x194 * x68,
            x164 * x198 * x91,
            x164 * x201 * x65,
            x192 * x43 * x96,
            x100 * x195 * x60,
            x194 * x60 * x98,
            x102 * x195 * x75,
            x100 * x194 * x75,
            x198 * x75 * x98,
            x106 * x163 * x195,
            x102 * x163 * x194,
            x100 * x163 * x198,
            x163 * x201 * x98,
            x203 * x43 * x65,
            x203 * x60 * x91,
            x205 * x60 * x65,
            x203 * x68 * x75,
            x205 * x75 * x91,
            x207 * x65 * x75,
            x163 * x203 * x86,
            x163 * x205 * x68,
            x163 * x207 * x91,
            x163 * x211 * x65,
            x119 * x195 * x42,
            x122 * x195 * x35,
            x119 * x194 * x35,
            x126 * x195 * x73,
            x122 * x194 * x73,
            x119 * x198 * x73,
            x127 * x185 * x191,
            x126 * x194 * x69,
            x122 * x198 * x69,
            x119 * x201 * x69,
            x203 * x42 * x98,
            x100 * x203 * x35,
            x205 * x35 * x98,
            x102 * x203 * x73,
            x100 * x205 * x73,
            x207 * x73 * x98,
            x106 * x203 * x69,
            x102 * x205 * x69,
            x100 * x207 * x69,
            x211 * x212 * x96,
            x213 * x42 * x65,
            x213 * x35 * x91,
            x216 * x35 * x65,
            x213 * x68 * x73,
            x216 * x73 * x91,
            x220 * x65 * x73,
            x213 * x69 * x86,
            x216 * x68 * x69,
            x212 * x220 * x51,
            x222 * x24,
            x138 * x195 * x22,
            x139 * x18 * x195,
            x138 * x18 * x194,
            x140 * x190 * x191,
            x139 * x194 * x8,
            x138 * x198 * x8,
            x142 * x191,
            x140 * x194 * x7,
            x139 * x198 * x7,
            x138 * x201 * x7,
            x119 * x203 * x22,
            x122 * x18 * x203,
            x119 * x18 * x205,
            x126 * x203 * x8,
            x122 * x205 * x8,
            x119 * x207 * x8,
            x127 * x203 * x7,
            x126 * x205 * x7,
            x122 * x207 * x7,
            x119 * x211 * x7,
            x213 * x22 * x98,
            x100 * x18 * x213,
            x18 * x216 * x98,
            x102 * x213 * x8,
            x100 * x216 * x8,
            x220 * x223 * x96,
            x106 * x213 * x7,
            x102 * x216 * x7,
            x100 * x220 * x7,
            x222 * x96,
            x22 * x224 * x65,
            x18 * x224 * x91,
            x18 * x225 * x65,
            x224 * x68 * x8,
            x223 * x225 * x51,
            x226 * x3,
            x224 * x7 * x86,
            x225 * x68 * x7,
            x226 * x51,
            x144
            * (
                x0 * (x137 + 2.0 * x209 + 2.0 * x210 + 3.0 * x218 + 3.0 * x219)
                + x107 * x221
            ),
        ]
    )


def dipole3d_34(a, A, b, B, C):
    """Cartesian 3D (fg) dipole moment integrals.
    The origin is at C.

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
    x11 = x4**2 * x8
    x12 = x0 * x8
    x13 = 3.0 * x12
    x14 = x11 + x13
    x15 = x0 * (2.0 * x10 + x14)
    x16 = 4.0 * x15
    x17 = x3 * x8
    x18 = x0 * (x17 + x9)
    x19 = x10 + x12
    x20 = x19 * x4
    x21 = x18 + x20
    x22 = x21 * x4
    x23 = x0 * (3.0 * x11 + x13)
    x24 = x11 + x12
    x25 = x24 * x4
    x26 = x12 * x4
    x27 = 2.0 * x26
    x28 = x25 + x27
    x29 = x28 * x4
    x30 = x23 + x29
    x31 = x0 * (x16 + 4.0 * x22 + x30)
    x32 = -x2 - A[0]
    x33 = 3.0 * x18
    x34 = x0 * (3.0 * x20 + x28 + x33)
    x35 = x15 + x22
    x36 = x35 * x4
    x37 = x34 + x36
    x38 = x32 * x37
    x39 = x21 * x32
    x40 = x28 * x32
    x41 = x23 + x40
    x42 = x0 * (x16 + x22 + 3.0 * x39 + x41)
    x43 = x32 * x35
    x44 = x34 + x43
    x45 = x32 * x44
    x46 = 8.0 * x26
    x47 = x0 * (4.0 * x25 + x46)
    x48 = x30 * x32
    x49 = x47 + x48
    x50 = x0 * (5.0 * x23 + x29 + 4.0 * x40) + x32 * x49
    x51 = x31 + x38
    x52 = x0 * (5.0 * x34 + x36 + 4.0 * x43 + x49) + x32 * x51
    x53 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x54 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x55 = 3.14159265358979 * x1 * x54
    x56 = x53 * x55
    x57 = -x1 * (a * A[1] + b * B[1])
    x58 = -x57 - B[1]
    x59 = x19 * x32
    x60 = 2.0 * x59
    x61 = x24 * x32
    x62 = x27 + x61
    x63 = x0 * (x20 + x33 + x60 + x62)
    x64 = x15 + x39
    x65 = x32 * x64
    x66 = x0 * (x25 + x46 + 3.0 * x61)
    x67 = x32 * x41
    x68 = x66 + x67
    x69 = x42 + x45
    x70 = x56 * (x0 * (2.0 * x34 + 2.0 * x43 + 3.0 * x63 + 3.0 * x65 + x68) + x32 * x69)
    x71 = -x1 * (a * A[2] + b * B[2])
    x72 = -x71 - B[2]
    x73 = x53 * x7
    x74 = x58**2 * x73
    x75 = x0 * x73
    x76 = x74 + x75
    x77 = x32 * x9
    x78 = x17 * x32
    x79 = x0 * (x10 + x13 + x77 + x78)
    x80 = x18 + x59
    x81 = x32 * x80
    x82 = 2.0 * x77
    x83 = x0 * (x14 + x82)
    x84 = x32 * x62
    x85 = x83 + x84
    x86 = x63 + x65
    x87 = x0 * (2.0 * x15 + 2.0 * x39 + 2.0 * x79 + 2.0 * x81 + x85) + x32 * x86
    x88 = x54 * x7
    x89 = x56 * x72
    x90 = x72**2 * x88
    x91 = x0 * x88
    x92 = x90 + x91
    x93 = x58 * x76
    x94 = x58 * x75
    x95 = 2.0 * x94
    x96 = x93 + x95
    x97 = x32 * x8
    x98 = x12 + x78
    x99 = x0 * (x17 + x97) + x32 * x98
    x100 = x0 * (x9 + x97)
    x101 = x12 + x77
    x102 = x101 * x32
    x103 = x100 + x102
    x104 = x79 + x81
    x105 = x0 * (x103 + 2.0 * x18 + x60 + x99) + x104 * x32
    x106 = x72 * x88
    x107 = x58 * x73
    x108 = x72 * x92
    x109 = x72 * x91
    x110 = 2.0 * x109
    x111 = x108 + x110
    x112 = 3.0 * x75
    x113 = x0 * (x112 + 3.0 * x74)
    x114 = x58 * x96
    x115 = x113 + x114
    x116 = x32**2 * x8
    x117 = x116 + x13
    x118 = x0 * (x117 + 2.0 * x78) + x32 * x99
    x119 = 3.0 * x91
    x120 = x0 * (x119 + 3.0 * x90)
    x121 = x111 * x72
    x122 = x120 + x121
    x123 = -x57 - A[1]
    x124 = x52 * x56
    x125 = x123 * x73
    x126 = x125 * x58
    x127 = x126 + x75
    x128 = x123 * x76
    x129 = x128 + x95
    x130 = x123 * x96
    x131 = x113 + x130
    x132 = 8.0 * x94
    x133 = x0 * (x132 + 4.0 * x93)
    x134 = x115 * x123
    x135 = x133 + x134
    x136 = -x71 - A[2]
    x137 = x136 * x56
    x138 = x136 * x88
    x139 = x138 * x72
    x140 = x139 + x91
    x141 = x136 * x92
    x142 = x110 + x141
    x143 = x111 * x136
    x144 = x120 + x143
    x145 = 8.0 * x109
    x146 = x0 * (4.0 * x108 + x145)
    x147 = x122 * x136
    x148 = x146 + x147
    x149 = x123**2 * x73
    x150 = x149 + x75
    x151 = x0 * (x107 + x125)
    x152 = x123 * x127
    x153 = x151 + x152
    x154 = x112 + 2.0 * x126
    x155 = x0 * (x154 + x74)
    x156 = x123 * x129
    x157 = x155 + x156
    x158 = x0 * (3.0 * x128 + x132 + x93)
    x159 = x123 * x131
    x160 = x158 + x159
    x161 = x0 * (5.0 * x113 + x114 + 4.0 * x130) + x123 * x135
    x162 = x136**2 * x88
    x163 = x162 + x91
    x164 = x0 * (x106 + x138)
    x165 = x136 * x140
    x166 = x164 + x165
    x167 = x119 + 2.0 * x139
    x168 = x0 * (x167 + x90)
    x169 = x136 * x142
    x170 = x168 + x169
    x171 = x0 * (x108 + 3.0 * x141 + x145)
    x172 = x136 * x144
    x173 = x171 + x172
    x174 = x0 * (5.0 * x120 + x121 + 4.0 * x143) + x136 * x148
    x175 = x123 * x150 + 2.0 * x123 * x75
    x176 = x0 * (x149 + x154) + x123 * x153
    x177 = x0 * (2.0 * x128 + 2.0 * x151 + 2.0 * x152 + 4.0 * x94) + x123 * x157
    x178 = x0 * (2.0 * x113 + 2.0 * x130 + 3.0 * x155 + 3.0 * x156) + x123 * x160
    x179 = x55 * x6
    x180 = x179 * (x0 * (2.0 * x133 + 2.0 * x134 + 4.0 * x158 + 4.0 * x159) + x123 * x161)
    x181 = x179 * x3
    x182 = 3.14159265358979 * x1 * x53 * x6
    x183 = x182 * x3
    x184 = x136 * x163 + 2.0 * x136 * x91
    x185 = x0 * (x162 + x167) + x136 * x166
    x186 = x0 * (4.0 * x109 + 2.0 * x141 + 2.0 * x164 + 2.0 * x165) + x136 * x170
    x187 = x0 * (2.0 * x120 + 2.0 * x143 + 3.0 * x168 + 3.0 * x169) + x136 * x173
    x188 = x182 * (x0 * (2.0 * x146 + 2.0 * x147 + 4.0 * x171 + 4.0 * x172) + x136 * x174)
    x189 = -x57 - C[1]
    x190 = x56 * (x0 * (2.0 * x47 + 2.0 * x48 + 4.0 * x66 + 4.0 * x67) + x32 * x50)
    x191 = x107 * x189
    x192 = x191 + x75
    x193 = x0 * (2.0 * x23 + 2.0 * x40 + 3.0 * x83 + 3.0 * x84) + x32 * x68
    x194 = x189 * x73
    x195 = x0 * (x107 + x194)
    x196 = x192 * x58
    x197 = x195 + x196
    x198 = x0 * (2.0 * x100 + 2.0 * x102 + 4.0 * x26 + 2.0 * x61) + x32 * x85
    x199 = x0 * (x112 + 2.0 * x191 + x74)
    x200 = x197 * x58
    x201 = x199 + x200
    x202 = x0 * (x117 + x82) + x103 * x32
    x203 = 3.0 * x195
    x204 = x0 * (3.0 * x196 + x203 + x96)
    x205 = x201 * x58
    x206 = x204 + x205
    x207 = x116 + x12
    x208 = 2.0 * x12 * x32 + x207 * x32
    x209 = x125 * x189
    x210 = x209 + x75
    x211 = x123 * x192
    x212 = x195 + x211
    x213 = x123 * x197
    x214 = x199 + x213
    x215 = x123 * x201
    x216 = x204 + x215
    x217 = 4.0 * x199
    x218 = x0 * (x115 + 4.0 * x200 + x217)
    x219 = x123 * x206
    x220 = x218 + x219
    x221 = x0 * (x125 + x194) + x123 * x210
    x222 = x0 * (x112 + x126 + x191 + x209)
    x223 = x123 * x212
    x224 = x222 + x223
    x225 = 2.0 * x211
    x226 = x0 * (x129 + x196 + x203 + x225)
    x227 = x123 * x214
    x228 = x226 + x227
    x229 = x0 * (x131 + x200 + 3.0 * x213 + x217)
    x230 = x123 * x216
    x231 = x229 + x230
    x232 = x0 * (x135 + 5.0 * x204 + x205 + 4.0 * x215) + x123 * x220
    x233 = x179 * x232
    x234 = x179 * x32
    x235 = x182 * x189
    x236 = x0 * (x112 + x149 + 2.0 * x209) + x123 * x221
    x237 = x0 * (x153 + 2.0 * x195 + x221 + x225) + x123 * x224
    x238 = x0 * (x157 + 2.0 * x199 + 2.0 * x213 + 2.0 * x222 + 2.0 * x223) + x123 * x228
    x239 = x179 * (
        x0 * (x160 + 2.0 * x204 + 2.0 * x215 + 3.0 * x226 + 3.0 * x227) + x123 * x231
    )
    x240 = x179 * x4
    x241 = -x71 - C[2]
    x242 = x241 * x56
    x243 = x106 * x241
    x244 = x243 + x91
    x245 = x241 * x88
    x246 = x0 * (x106 + x245)
    x247 = x244 * x72
    x248 = x246 + x247
    x249 = x0 * (x119 + 2.0 * x243 + x90)
    x250 = x248 * x72
    x251 = x249 + x250
    x252 = 3.0 * x246
    x253 = x0 * (x111 + 3.0 * x247 + x252)
    x254 = x251 * x72
    x255 = x253 + x254
    x256 = x138 * x241
    x257 = x256 + x91
    x258 = x136 * x244
    x259 = x246 + x258
    x260 = x136 * x248
    x261 = x249 + x260
    x262 = x136 * x251
    x263 = x253 + x262
    x264 = 4.0 * x249
    x265 = x0 * (x122 + 4.0 * x250 + x264)
    x266 = x136 * x255
    x267 = x265 + x266
    x268 = x182 * x32
    x269 = x0 * (x138 + x245) + x136 * x257
    x270 = x0 * (x119 + x139 + x243 + x256)
    x271 = x136 * x259
    x272 = x270 + x271
    x273 = 2.0 * x258
    x274 = x0 * (x142 + x247 + x252 + x273)
    x275 = x136 * x261
    x276 = x274 + x275
    x277 = x0 * (x144 + x250 + 3.0 * x260 + x264)
    x278 = x136 * x263
    x279 = x277 + x278
    x280 = x0 * (x148 + 5.0 * x253 + x254 + 4.0 * x262) + x136 * x267
    x281 = x182 * x280
    x282 = x182 * x4
    x283 = x0 * (x119 + x162 + 2.0 * x256) + x136 * x269
    x284 = x0 * (x166 + 2.0 * x246 + x269 + x273) + x136 * x272
    x285 = x0 * (x170 + 2.0 * x249 + 2.0 * x260 + 2.0 * x270 + 2.0 * x271) + x136 * x276
    x286 = x182 * (
        x0 * (x173 + 2.0 * x253 + 2.0 * x262 + 3.0 * x274 + 3.0 * x275) + x136 * x279
    )

    # 450 item(s)
    return numpy.array(
        [
            x56
            * (x0 * (2.0 * x31 + 2.0 * x38 + 4.0 * x42 + 4.0 * x45 + x50) + x32 * x52),
            x58 * x70,
            x70 * x72,
            x76 * x87 * x88,
            x58 * x87 * x89,
            x73 * x87 * x92,
            x105 * x88 * x96,
            x105 * x106 * x76,
            x105 * x107 * x92,
            x105 * x111 * x73,
            x115 * x118 * x88,
            x106 * x118 * x96,
            x118 * x76 * x92,
            x107 * x111 * x118,
            x118 * x122 * x73,
            x123 * x124,
            x127 * x69 * x88,
            x123 * x69 * x89,
            x129 * x86 * x88,
            x106 * x127 * x86,
            x125 * x86 * x92,
            x104 * x131 * x88,
            x104 * x106 * x129,
            x104 * x127 * x92,
            x104 * x111 * x125,
            x135 * x88 * x99,
            x106 * x131 * x99,
            x129 * x92 * x99,
            x111 * x127 * x99,
            x122 * x125 * x99,
            x124 * x136,
            x137 * x58 * x69,
            x140 * x69 * x73,
            x138 * x76 * x86,
            x107 * x140 * x86,
            x142 * x73 * x86,
            x104 * x138 * x96,
            x104 * x140 * x76,
            x104 * x107 * x142,
            x104 * x144 * x73,
            x115 * x138 * x99,
            x140 * x96 * x99,
            x142 * x76 * x99,
            x107 * x144 * x99,
            x148 * x73 * x99,
            x150 * x51 * x88,
            x153 * x44 * x88,
            x106 * x150 * x44,
            x157 * x64 * x88,
            x106 * x153 * x64,
            x150 * x64 * x92,
            x160 * x80 * x88,
            x106 * x157 * x80,
            x153 * x80 * x92,
            x111 * x150 * x80,
            x161 * x88 * x98,
            x106 * x160 * x98,
            x157 * x92 * x98,
            x111 * x153 * x98,
            x122 * x150 * x98,
            x123 * x137 * x51,
            x127 * x138 * x44,
            x125 * x140 * x44,
            x129 * x138 * x64,
            x127 * x140 * x64,
            x125 * x142 * x64,
            x131 * x138 * x80,
            x129 * x140 * x80,
            x127 * x142 * x80,
            x125 * x144 * x80,
            x135 * x138 * x98,
            x131 * x140 * x98,
            x129 * x142 * x98,
            x127 * x144 * x98,
            x125 * x148 * x98,
            x163 * x51 * x73,
            x107 * x163 * x44,
            x166 * x44 * x73,
            x163 * x64 * x76,
            x107 * x166 * x64,
            x170 * x64 * x73,
            x163 * x80 * x96,
            x166 * x76 * x80,
            x107 * x170 * x80,
            x173 * x73 * x80,
            x115 * x163 * x98,
            x166 * x96 * x98,
            x170 * x76 * x98,
            x107 * x173 * x98,
            x174 * x73 * x98,
            x175 * x37 * x88,
            x176 * x35 * x88,
            x106 * x175 * x35,
            x177 * x21 * x88,
            x106 * x176 * x21,
            x175 * x21 * x92,
            x178 * x19 * x88,
            x106 * x177 * x19,
            x176 * x19 * x92,
            x111 * x175 * x19,
            x180 * x3,
            x178 * x181 * x72,
            x17 * x177 * x92,
            x111 * x17 * x176,
            x122 * x17 * x175,
            x138 * x150 * x37,
            x138 * x153 * x35,
            x140 * x150 * x35,
            x138 * x157 * x21,
            x140 * x153 * x21,
            x142 * x150 * x21,
            x138 * x160 * x19,
            x140 * x157 * x19,
            x142 * x153 * x19,
            x144 * x150 * x19,
            x136 * x161 * x181,
            x140 * x160 * x17,
            x142 * x157 * x17,
            x144 * x153 * x17,
            x148 * x150 * x17,
            x125 * x163 * x37,
            x127 * x163 * x35,
            x125 * x166 * x35,
            x129 * x163 * x21,
            x127 * x166 * x21,
            x125 * x170 * x21,
            x131 * x163 * x19,
            x129 * x166 * x19,
            x127 * x170 * x19,
            x125 * x173 * x19,
            x135 * x163 * x17,
            x131 * x166 * x17,
            x129 * x17 * x170,
            x127 * x17 * x173,
            x123 * x174 * x183,
            x184 * x37 * x73,
            x107 * x184 * x35,
            x185 * x35 * x73,
            x184 * x21 * x76,
            x107 * x185 * x21,
            x186 * x21 * x73,
            x184 * x19 * x96,
            x185 * x19 * x76,
            x107 * x186 * x19,
            x187 * x19 * x73,
            x115 * x17 * x184,
            x17 * x185 * x96,
            x17 * x186 * x76,
            x183 * x187 * x58,
            x188 * x3,
            x189 * x190,
            x192 * x193 * x88,
            x189 * x193 * x89,
            x197 * x198 * x88,
            x106 * x192 * x198,
            x194 * x198 * x92,
            x201 * x202 * x88,
            x106 * x197 * x202,
            x192 * x202 * x92,
            x111 * x194 * x202,
            x206 * x208 * x88,
            x106 * x201 * x208,
            x197 * x208 * x92,
            x111 * x192 * x208,
            x122 * x194 * x208,
            x210 * x50 * x88,
            x212 * x68 * x88,
            x106 * x210 * x68,
            x214 * x85 * x88,
            x106 * x212 * x85,
            x210 * x85 * x92,
            x103 * x216 * x88,
            x103 * x106 * x214,
            x103 * x212 * x92,
            x103 * x111 * x210,
            x207 * x220 * x88,
            x106 * x207 * x216,
            x207 * x214 * x92,
            x111 * x207 * x212,
            x122 * x207 * x210,
            x137 * x189 * x50,
            x138 * x192 * x68,
            x140 * x194 * x68,
            x138 * x197 * x85,
            x140 * x192 * x85,
            x142 * x194 * x85,
            x103 * x138 * x201,
            x103 * x140 * x197,
            x103 * x142 * x192,
            x103 * x144 * x194,
            x138 * x206 * x207,
            x140 * x201 * x207,
            x142 * x197 * x207,
            x144 * x192 * x207,
            x148 * x194 * x207,
            x221 * x49 * x88,
            x224 * x41 * x88,
            x106 * x221 * x41,
            x228 * x62 * x88,
            x106 * x224 * x62,
            x221 * x62 * x92,
            x101 * x231 * x88,
            x101 * x106 * x228,
            x101 * x224 * x92,
            x101 * x111 * x221,
            x233 * x32,
            x231 * x234 * x72,
            x228 * x92 * x97,
            x111 * x224 * x97,
            x122 * x221 * x97,
            x138 * x210 * x49,
            x138 * x212 * x41,
            x140 * x210 * x41,
            x138 * x214 * x62,
            x140 * x212 * x62,
            x142 * x210 * x62,
            x101 * x138 * x216,
            x101 * x140 * x214,
            x101 * x142 * x212,
            x101 * x144 * x210,
            x136 * x220 * x234,
            x140 * x216 * x97,
            x142 * x214 * x97,
            x144 * x212 * x97,
            x148 * x210 * x97,
            x163 * x194 * x49,
            x163 * x192 * x41,
            x166 * x194 * x41,
            x163 * x197 * x62,
            x166 * x192 * x62,
            x170 * x194 * x62,
            x101 * x163 * x201,
            x101 * x166 * x197,
            x101 * x170 * x192,
            x101 * x173 * x194,
            x163 * x206 * x97,
            x166 * x201 * x97,
            x170 * x197 * x97,
            x173 * x192 * x97,
            x174 * x235 * x32,
            x236 * x30 * x88,
            x237 * x28 * x88,
            x106 * x236 * x28,
            x238 * x24 * x88,
            x106 * x237 * x24,
            x236 * x24 * x92,
            x239 * x4,
            x238 * x240 * x72,
            x237 * x9 * x92,
            x111 * x236 * x9,
            x179
            * (
                x0 * (x161 + 2.0 * x218 + 2.0 * x219 + 4.0 * x229 + 4.0 * x230)
                + x123 * x232
            ),
            x239 * x72,
            x238 * x8 * x92,
            x111 * x237 * x8,
            x122 * x236 * x8,
            x138 * x221 * x30,
            x138 * x224 * x28,
            x140 * x221 * x28,
            x138 * x228 * x24,
            x140 * x224 * x24,
            x142 * x221 * x24,
            x136 * x231 * x240,
            x140 * x228 * x9,
            x142 * x224 * x9,
            x144 * x221 * x9,
            x136 * x233,
            x140 * x231 * x8,
            x142 * x228 * x8,
            x144 * x224 * x8,
            x148 * x221 * x8,
            x163 * x210 * x30,
            x163 * x212 * x28,
            x166 * x210 * x28,
            x163 * x214 * x24,
            x166 * x212 * x24,
            x170 * x210 * x24,
            x163 * x216 * x9,
            x166 * x214 * x9,
            x170 * x212 * x9,
            x173 * x210 * x9,
            x163 * x220 * x8,
            x166 * x216 * x8,
            x170 * x214 * x8,
            x173 * x212 * x8,
            x174 * x210 * x8,
            x184 * x194 * x30,
            x184 * x192 * x28,
            x185 * x194 * x28,
            x184 * x197 * x24,
            x185 * x192 * x24,
            x186 * x194 * x24,
            x184 * x201 * x9,
            x185 * x197 * x9,
            x186 * x192 * x9,
            x187 * x235 * x4,
            x184 * x206 * x8,
            x185 * x201 * x8,
            x186 * x197 * x8,
            x187 * x192 * x8,
            x188 * x189,
            x190 * x241,
            x193 * x242 * x58,
            x193 * x244 * x73,
            x198 * x245 * x76,
            x107 * x198 * x244,
            x198 * x248 * x73,
            x202 * x245 * x96,
            x202 * x244 * x76,
            x107 * x202 * x248,
            x202 * x251 * x73,
            x115 * x208 * x245,
            x208 * x244 * x96,
            x208 * x248 * x76,
            x107 * x208 * x251,
            x208 * x255 * x73,
            x123 * x242 * x50,
            x127 * x245 * x68,
            x125 * x244 * x68,
            x129 * x245 * x85,
            x127 * x244 * x85,
            x125 * x248 * x85,
            x103 * x131 * x245,
            x103 * x129 * x244,
            x103 * x127 * x248,
            x103 * x125 * x251,
            x135 * x207 * x245,
            x131 * x207 * x244,
            x129 * x207 * x248,
            x127 * x207 * x251,
            x125 * x207 * x255,
            x257 * x50 * x73,
            x107 * x257 * x68,
            x259 * x68 * x73,
            x257 * x76 * x85,
            x107 * x259 * x85,
            x261 * x73 * x85,
            x103 * x257 * x96,
            x103 * x259 * x76,
            x103 * x107 * x261,
            x103 * x263 * x73,
            x115 * x207 * x257,
            x207 * x259 * x96,
            x207 * x261 * x76,
            x107 * x207 * x263,
            x207 * x267 * x73,
            x150 * x245 * x49,
            x153 * x245 * x41,
            x150 * x244 * x41,
            x157 * x245 * x62,
            x153 * x244 * x62,
            x150 * x248 * x62,
            x101 * x160 * x245,
            x101 * x157 * x244,
            x101 * x153 * x248,
            x101 * x150 * x251,
            x161 * x234 * x241,
            x160 * x244 * x97,
            x157 * x248 * x97,
            x153 * x251 * x97,
            x150 * x255 * x97,
            x125 * x257 * x49,
            x127 * x257 * x41,
            x125 * x259 * x41,
            x129 * x257 * x62,
            x127 * x259 * x62,
            x125 * x261 * x62,
            x101 * x131 * x257,
            x101 * x129 * x259,
            x101 * x127 * x261,
            x101 * x125 * x263,
            x135 * x257 * x97,
            x131 * x259 * x97,
            x129 * x261 * x97,
            x127 * x263 * x97,
            x123 * x267 * x268,
            x269 * x49 * x73,
            x107 * x269 * x41,
            x272 * x41 * x73,
            x269 * x62 * x76,
            x107 * x272 * x62,
            x276 * x62 * x73,
            x101 * x269 * x96,
            x101 * x272 * x76,
            x101 * x107 * x276,
            x101 * x279 * x73,
            x115 * x269 * x97,
            x272 * x96 * x97,
            x276 * x76 * x97,
            x268 * x279 * x58,
            x281 * x32,
            x175 * x245 * x30,
            x176 * x245 * x28,
            x175 * x244 * x28,
            x177 * x24 * x245,
            x176 * x24 * x244,
            x175 * x24 * x248,
            x178 * x240 * x241,
            x177 * x244 * x9,
            x176 * x248 * x9,
            x175 * x251 * x9,
            x180 * x241,
            x178 * x244 * x8,
            x177 * x248 * x8,
            x176 * x251 * x8,
            x175 * x255 * x8,
            x150 * x257 * x30,
            x153 * x257 * x28,
            x150 * x259 * x28,
            x157 * x24 * x257,
            x153 * x24 * x259,
            x150 * x24 * x261,
            x160 * x257 * x9,
            x157 * x259 * x9,
            x153 * x261 * x9,
            x150 * x263 * x9,
            x161 * x257 * x8,
            x160 * x259 * x8,
            x157 * x261 * x8,
            x153 * x263 * x8,
            x150 * x267 * x8,
            x125 * x269 * x30,
            x127 * x269 * x28,
            x125 * x272 * x28,
            x129 * x24 * x269,
            x127 * x24 * x272,
            x125 * x24 * x276,
            x131 * x269 * x9,
            x129 * x272 * x9,
            x127 * x276 * x9,
            x123 * x279 * x282,
            x135 * x269 * x8,
            x131 * x272 * x8,
            x129 * x276 * x8,
            x127 * x279 * x8,
            x123 * x281,
            x283 * x30 * x73,
            x107 * x28 * x283,
            x28 * x284 * x73,
            x24 * x283 * x76,
            x107 * x24 * x284,
            x24 * x285 * x73,
            x283 * x9 * x96,
            x284 * x76 * x9,
            x282 * x285 * x58,
            x286 * x4,
            x115 * x283 * x8,
            x284 * x8 * x96,
            x285 * x76 * x8,
            x286 * x58,
            x182
            * (
                x0 * (x174 + 2.0 * x265 + 2.0 * x266 + 4.0 * x277 + 4.0 * x278)
                + x136 * x280
            ),
        ]
    )


def dipole3d_40(a, A, b, B, C):
    """Cartesian 3D (gs) dipole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - A[0]
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
    x16 = x3**2 * x7
    x17 = x12 + x16
    x18 = 2.0 * x12 * x3 + x17 * x3
    x19 = 3.0 * x12
    x20 = x11 + x15
    x21 = x0 * (2.0 * x13 + x16 + x19) + x20 * x3
    x22 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x23 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x24 = 3.14159265358979 * x1 * x23
    x25 = x22 * x24
    x26 = -x1 * (a * A[1] + b * B[1])
    x27 = -x26 - A[1]
    x28 = x21 * x25
    x29 = -x1 * (a * A[2] + b * B[2])
    x30 = -x29 - A[2]
    x31 = x22 * x6
    x32 = x27**2 * x31
    x33 = x0 * x31
    x34 = x32 + x33
    x35 = x23 * x6
    x36 = x25 * x30
    x37 = x30**2 * x35
    x38 = x0 * x35
    x39 = x37 + x38
    x40 = 2.0 * x27 * x33 + x27 * x34
    x41 = x30 * x35
    x42 = x27 * x31
    x43 = 2.0 * x30 * x38 + x30 * x39
    x44 = 3.0 * x33
    x45 = x24 * x5
    x46 = x45 * (x0 * (3.0 * x32 + x44) + x27 * x40)
    x47 = x30 * x45
    x48 = 3.14159265358979 * x1 * x22 * x5
    x49 = x43 * x48
    x50 = 3.0 * x38
    x51 = x48 * (x0 * (3.0 * x37 + x50) + x30 * x43)
    x52 = -x26 - C[1]
    x53 = x25 * (x0 * (3.0 * x16 + x19) + x18 * x3)
    x54 = x42 * x52
    x55 = x33 + x54
    x56 = x31 * x52
    x57 = x0 * (x42 + x56)
    x58 = x27 * x55
    x59 = x57 + x58
    x60 = x0 * (x32 + x44 + 2.0 * x54) + x27 * x59
    x61 = x45 * x60
    x62 = -x29 - C[2]
    x63 = x41 * x62
    x64 = x38 + x63
    x65 = x35 * x62
    x66 = x0 * (x41 + x65)
    x67 = x30 * x64
    x68 = x66 + x67
    x69 = x0 * (x37 + x50 + 2.0 * x63) + x30 * x68
    x70 = x48 * x69

    # 45 item(s)
    return numpy.array(
        [
            x25 * (x0 * (3.0 * x11 + 3.0 * x15 + x18) + x21 * x3),
            x27 * x28,
            x28 * x30,
            x20 * x34 * x35,
            x20 * x27 * x36,
            x20 * x31 * x39,
            x14 * x35 * x40,
            x14 * x34 * x41,
            x14 * x39 * x42,
            x14 * x31 * x43,
            x46 * x9,
            x40 * x47 * x9,
            x10 * x34 * x39,
            x27 * x49 * x9,
            x51 * x9,
            x52 * x53,
            x18 * x35 * x55,
            x18 * x36 * x52,
            x17 * x35 * x59,
            x17 * x41 * x55,
            x17 * x39 * x56,
            x3 * x61,
            x3 * x47 * x59,
            x39 * x55 * x8,
            x3 * x49 * x52,
            x45 * (x0 * (x40 + 3.0 * x57 + 3.0 * x58) + x27 * x60),
            x30 * x61,
            x39 * x59 * x7,
            x43 * x55 * x7,
            x51 * x52,
            x53 * x62,
            x18 * x25 * x27 * x62,
            x18 * x31 * x64,
            x17 * x34 * x65,
            x17 * x42 * x64,
            x17 * x31 * x68,
            x3 * x40 * x45 * x62,
            x34 * x64 * x8,
            x27 * x3 * x48 * x68,
            x3 * x70,
            x46 * x62,
            x40 * x64 * x7,
            x34 * x68 * x7,
            x27 * x70,
            x48 * (x0 * (x43 + 3.0 * x66 + 3.0 * x67) + x30 * x69),
        ]
    )


def dipole3d_41(a, A, b, B, C):
    """Cartesian 3D (gp) dipole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = a * b * x1
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = 1.77245385090552 * numpy.sqrt(x1)
    x5 = x3 * x4
    x6 = x0 * x5
    x7 = 3.0 * x6
    x8 = -x1 * (a * A[0] + b * B[0])
    x9 = -x8 - B[0]
    x10 = -x8 - A[0]
    x11 = x10 * x5
    x12 = x11 * x9
    x13 = -x8 - C[0]
    x14 = x11 * x13
    x15 = x5 * x9
    x16 = x13 * x15
    x17 = x0 * (x12 + x14 + x16 + x7)
    x18 = x13 * x5
    x19 = x0 * (x15 + x18)
    x20 = x16 + x6
    x21 = x10 * x20
    x22 = x19 + x21
    x23 = x10 * x22
    x24 = x10**2 * x5
    x25 = x24 + x7
    x26 = x0 * (x11 + x18)
    x27 = x14 + x6
    x28 = x10 * x27
    x29 = x26 + x28
    x30 = x0 * (2.0 * x14 + x25) + x10 * x29
    x31 = x0 * (x11 + x15)
    x32 = x12 + x6
    x33 = x10 * x32
    x34 = x31 + x33
    x35 = x0 * (2.0 * x12 + x25) + x10 * x34
    x36 = x17 + x23
    x37 = x0 * (2.0 * x19 + 2.0 * x21 + x29 + x34) + x10 * x36
    x38 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x39 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x40 = 3.14159265358979 * x1 * x39
    x41 = x38 * x40
    x42 = -x1 * (a * A[1] + b * B[1])
    x43 = -x42 - B[1]
    x44 = x24 + x6
    x45 = 2.0 * x0 * x11 + x10 * x44
    x46 = x41 * (x0 * (3.0 * x26 + 3.0 * x28 + x45) + x10 * x30)
    x47 = -x1 * (a * A[2] + b * B[2])
    x48 = -x47 - B[2]
    x49 = -x42 - A[1]
    x50 = x37 * x41
    x51 = x0 * x4
    x52 = x38 * x51
    x53 = x38 * x4
    x54 = x49 * x53
    x55 = x43 * x54
    x56 = x52 + x55
    x57 = x39 * x4
    x58 = x30 * x41
    x59 = -x47 - A[2]
    x60 = x39 * x51
    x61 = x57 * x59
    x62 = x48 * x61
    x63 = x60 + x62
    x64 = x49**2 * x53
    x65 = x52 + x64
    x66 = x43 * x53
    x67 = x0 * (x54 + x66)
    x68 = x49 * x56
    x69 = x67 + x68
    x70 = x48 * x57
    x71 = x41 * x59
    x72 = x57 * x59**2
    x73 = x60 + x72
    x74 = x0 * (x61 + x70)
    x75 = x59 * x63
    x76 = x74 + x75
    x77 = 2.0 * x49 * x52 + x49 * x65
    x78 = 3.0 * x52
    x79 = x64 + x78
    x80 = x0 * (2.0 * x55 + x79) + x49 * x69
    x81 = 2.0 * x59 * x60 + x59 * x73
    x82 = 3.0 * x60
    x83 = x72 + x82
    x84 = x0 * (2.0 * x62 + x83) + x59 * x76
    x85 = x0 * (3.0 * x64 + x78) + x49 * x77
    x86 = x3 * x40
    x87 = x86 * (x0 * (3.0 * x67 + 3.0 * x68 + x77) + x49 * x80)
    x88 = x13 * x86
    x89 = 3.14159265358979 * x1 * x3 * x38
    x90 = x13 * x89
    x91 = x0 * (3.0 * x72 + x82) + x59 * x81
    x92 = x89 * (x0 * (3.0 * x74 + 3.0 * x75 + x81) + x59 * x84)
    x93 = -x42 - C[1]
    x94 = x41 * (x0 * (3.0 * x31 + 3.0 * x33 + x45) + x10 * x35)
    x95 = x66 * x93
    x96 = x52 + x95
    x97 = x0 * (3.0 * x24 + x7) + x10 * x45
    x98 = x41 * x97
    x99 = x54 * x93
    x100 = x52 + x99
    x101 = x53 * x93
    x102 = x0 * (x101 + x66)
    x103 = x49 * x96
    x104 = x102 + x103
    x105 = x0 * (x101 + x54)
    x106 = x100 * x49
    x107 = x105 + x106
    x108 = x0 * (x55 + x78 + x95 + x99)
    x109 = x104 * x49
    x110 = x108 + x109
    x111 = x0 * (x79 + 2.0 * x99) + x107 * x49
    x112 = x0 * (2.0 * x102 + 2.0 * x103 + x107 + x69) + x110 * x49
    x113 = x112 * x86
    x114 = x10 * x86
    x115 = x89 * x93
    x116 = x86 * (x0 * (3.0 * x105 + 3.0 * x106 + x77) + x111 * x49)
    x117 = x86 * x9
    x118 = -x47 - C[2]
    x119 = x118 * x70
    x120 = x119 + x60
    x121 = x118 * x57
    x122 = x118 * x61
    x123 = x122 + x60
    x124 = x0 * (x121 + x70)
    x125 = x120 * x59
    x126 = x124 + x125
    x127 = x0 * (x121 + x61)
    x128 = x123 * x59
    x129 = x127 + x128
    x130 = x0 * (x119 + x122 + x62 + x82)
    x131 = x126 * x59
    x132 = x130 + x131
    x133 = x10 * x89
    x134 = x0 * (2.0 * x122 + x83) + x129 * x59
    x135 = x0 * (2.0 * x124 + 2.0 * x125 + x129 + x76) + x132 * x59
    x136 = x135 * x89
    x137 = x89 * (x0 * (3.0 * x127 + 3.0 * x128 + x81) + x134 * x59)

    # 135 item(s)
    return numpy.array(
        [
            x41 * (x0 * (3.0 * x17 + 3.0 * x23 + x30 + x35) + x10 * x37),
            x43 * x46,
            x46 * x48,
            x49 * x50,
            x30 * x56 * x57,
            x48 * x49 * x58,
            x50 * x59,
            x43 * x58 * x59,
            x30 * x53 * x63,
            x36 * x57 * x65,
            x29 * x57 * x69,
            x29 * x65 * x70,
            x36 * x49 * x71,
            x29 * x56 * x61,
            x29 * x54 * x63,
            x36 * x53 * x73,
            x29 * x66 * x73,
            x29 * x53 * x76,
            x22 * x57 * x77,
            x27 * x57 * x80,
            x27 * x70 * x77,
            x22 * x61 * x65,
            x27 * x61 * x69,
            x27 * x63 * x65,
            x22 * x54 * x73,
            x27 * x56 * x73,
            x27 * x54 * x76,
            x22 * x53 * x81,
            x27 * x66 * x81,
            x27 * x53 * x84,
            x20 * x57 * x85,
            x13 * x87,
            x48 * x85 * x88,
            x20 * x61 * x77,
            x59 * x80 * x88,
            x18 * x63 * x77,
            x20 * x65 * x73,
            x18 * x69 * x73,
            x18 * x65 * x76,
            x20 * x54 * x81,
            x18 * x56 * x81,
            x49 * x84 * x90,
            x20 * x53 * x91,
            x43 * x90 * x91,
            x13 * x92,
            x93 * x94,
            x57 * x96 * x97,
            x48 * x93 * x98,
            x100 * x35 * x57,
            x104 * x45 * x57,
            x100 * x45 * x70,
            x35 * x71 * x93,
            x45 * x61 * x96,
            x101 * x45 * x63,
            x107 * x34 * x57,
            x110 * x44 * x57,
            x107 * x44 * x70,
            x100 * x34 * x61,
            x104 * x44 * x61,
            x100 * x44 * x63,
            x101 * x34 * x73,
            x44 * x73 * x96,
            x101 * x44 * x76,
            x111 * x32 * x57,
            x10 * x113,
            x111 * x114 * x48,
            x107 * x32 * x61,
            x110 * x114 * x59,
            x107 * x11 * x63,
            x100 * x32 * x73,
            x104 * x11 * x73,
            x100 * x11 * x76,
            x101 * x32 * x81,
            x11 * x81 * x96,
            x10 * x115 * x84,
            x116 * x9,
            x86 * (x0 * (3.0 * x108 + 3.0 * x109 + x111 + x80) + x112 * x49),
            x116 * x48,
            x111 * x117 * x59,
            x113 * x59,
            x111 * x5 * x63,
            x107 * x15 * x73,
            x110 * x5 * x73,
            x107 * x5 * x76,
            x100 * x15 * x81,
            x104 * x5 * x81,
            x100 * x5 * x84,
            x115 * x9 * x91,
            x5 * x91 * x96,
            x92 * x93,
            x118 * x94,
            x118 * x43 * x98,
            x120 * x53 * x97,
            x118 * x35 * x41 * x49,
            x121 * x45 * x56,
            x120 * x45 * x54,
            x123 * x35 * x53,
            x123 * x45 * x66,
            x126 * x45 * x53,
            x121 * x34 * x65,
            x121 * x44 * x69,
            x120 * x44 * x65,
            x123 * x34 * x54,
            x123 * x44 * x56,
            x126 * x44 * x54,
            x129 * x34 * x53,
            x129 * x44 * x66,
            x132 * x44 * x53,
            x121 * x32 * x77,
            x114 * x118 * x80,
            x11 * x120 * x77,
            x123 * x32 * x65,
            x11 * x123 * x69,
            x11 * x126 * x65,
            x129 * x32 * x54,
            x11 * x129 * x56,
            x132 * x133 * x49,
            x134 * x32 * x53,
            x133 * x134 * x43,
            x10 * x136,
            x117 * x118 * x85,
            x118 * x87,
            x120 * x5 * x85,
            x123 * x15 * x77,
            x123 * x5 * x80,
            x126 * x5 * x77,
            x129 * x15 * x65,
            x129 * x5 * x69,
            x132 * x5 * x65,
            x134 * x49 * x89 * x9,
            x134 * x5 * x56,
            x136 * x49,
            x137 * x9,
            x137 * x43,
            x89 * (x0 * (3.0 * x130 + 3.0 * x131 + x134 + x84) + x135 * x59),
        ]
    )


def dipole3d_42(a, A, b, B, C):
    """Cartesian 3D (gd) dipole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2.0 * a + 2.0 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = -x1 * (a * A[0] + b * B[0])
    x3 = -x2 - A[0]
    x4 = a * b * x1
    x5 = numpy.exp(-x4 * (A[0] - B[0]) ** 2)
    x6 = 1.77245385090552 * numpy.sqrt(x1)
    x7 = x5 * x6
    x8 = x0 * x7
    x9 = -x2 - C[0]
    x10 = -x2 - B[0]
    x11 = x10 * x7
    x12 = x11 * x9
    x13 = x12 + x8
    x14 = x13 * x3
    x15 = 2.0 * x14
    x16 = x7 * x9
    x17 = x0 * (x11 + x16)
    x18 = x3 * x5
    x19 = x18 * x6
    x20 = x0 * (x11 + x19)
    x21 = x10 * x19
    x22 = x21 + x8
    x23 = x22 * x3
    x24 = x20 + x23
    x25 = x0 * (x16 + x19)
    x26 = x19 * x9
    x27 = x26 + x8
    x28 = x27 * x3
    x29 = x25 + x28
    x30 = x0 * (x15 + 2.0 * x17 + x24 + x29)
    x31 = 3.0 * x8
    x32 = x0 * (x12 + x21 + x26 + x31)
    x33 = x14 + x17
    x34 = x3 * x33
    x35 = x32 + x34
    x36 = x3 * x35
    x37 = x10 * x13
    x38 = x10**2 * x7
    x39 = x38 + x8
    x40 = x3 * x39
    x41 = x0 * x11
    x42 = x40 + 2.0 * x41
    x43 = x0 * (x15 + 3.0 * x17 + x37 + x42)
    x44 = x31 + x38
    x45 = x0 * (2.0 * x12 + x44)
    x46 = x17 + x37
    x47 = x3 * x46
    x48 = x45 + x47
    x49 = x3 * x48
    x50 = 2.0 * x21
    x51 = x0 * (x44 + x50)
    x52 = x3 * x42
    x53 = x51 + x52
    x54 = x0 * (2.0 * x20 + 2.0 * x23 + 2.0 * x40 + 4.0 * x41) + x3 * x53
    x55 = x43 + x49
    x56 = x0 * (2.0 * x32 + 2.0 * x34 + 2.0 * x45 + 2.0 * x47 + x53) + x3 * x55
    x57 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x58 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x59 = 3.14159265358979 * x1 * x58
    x60 = x57 * x59
    x61 = -x1 * (a * A[1] + b * B[1])
    x62 = -x61 - B[1]
    x63 = x3**2 * x7
    x64 = x31 + x63
    x65 = x0 * (2.0 * x26 + x64) + x29 * x3
    x66 = x0 * (x50 + x64)
    x67 = x24 * x3
    x68 = x66 + x67
    x69 = x30 + x36
    x70 = x60 * (x0 * (3.0 * x32 + 3.0 * x34 + x65 + x68) + x3 * x69)
    x71 = -x1 * (a * A[2] + b * B[2])
    x72 = -x71 - B[2]
    x73 = x57 * x6
    x74 = x62**2 * x73
    x75 = x0 * x73
    x76 = x74 + x75
    x77 = x63 + x8
    x78 = 2.0 * x0 * x19 + x3 * x77
    x79 = x0 * (3.0 * x25 + 3.0 * x28 + x78) + x3 * x65
    x80 = x58 * x6
    x81 = x60 * x72
    x82 = x72**2 * x80
    x83 = x0 * x80
    x84 = x82 + x83
    x85 = -x61 - A[1]
    x86 = x56 * x60
    x87 = x73 * x85
    x88 = x62 * x87
    x89 = x75 + x88
    x90 = x76 * x85
    x91 = 2.0 * x75
    x92 = x62 * x91 + x90
    x93 = x72 * x80
    x94 = -x71 - A[2]
    x95 = x60 * x94
    x96 = x80 * x94
    x97 = x72 * x96
    x98 = x83 + x97
    x99 = x62 * x73
    x100 = x84 * x94
    x101 = 2.0 * x83
    x102 = x100 + x101 * x72
    x103 = x73 * x85**2
    x104 = x103 + x75
    x105 = x0 * (x87 + x99)
    x106 = x85 * x89
    x107 = x105 + x106
    x108 = 3.0 * x75
    x109 = x108 + 2.0 * x88
    x110 = x0 * (x109 + x74)
    x111 = x85 * x92
    x112 = x110 + x111
    x113 = x80 * x94**2
    x114 = x113 + x83
    x115 = x0 * (x93 + x96)
    x116 = x94 * x98
    x117 = x115 + x116
    x118 = 3.0 * x83
    x119 = x118 + 2.0 * x97
    x120 = x0 * (x119 + x82)
    x121 = x102 * x94
    x122 = x120 + x121
    x123 = x104 * x85 + x85 * x91
    x124 = x0 * (x103 + x109)
    x125 = x107 * x85
    x126 = x124 + x125
    x127 = x0 * (2.0 * x105 + 2.0 * x106 + 4.0 * x62 * x75 + 2.0 * x90) + x112 * x85
    x128 = x101 * x94 + x114 * x94
    x129 = x0 * (x113 + x119)
    x130 = x117 * x94
    x131 = x129 + x130
    x132 = x0 * (2.0 * x100 + 2.0 * x115 + 2.0 * x116 + 4.0 * x72 * x83) + x122 * x94
    x133 = x0 * (3.0 * x103 + x108) + x123 * x85
    x134 = x0 * (3.0 * x105 + 3.0 * x106 + x123) + x126 * x85
    x135 = x5 * x59
    x136 = x135 * (x0 * (3.0 * x110 + 3.0 * x111 + 2.0 * x124 + 2.0 * x125) + x127 * x85)
    x137 = x135 * x72
    x138 = x135 * x94
    x139 = 3.14159265358979 * x1 * x57
    x140 = x139 * x5
    x141 = x140 * x85
    x142 = x0 * (3.0 * x113 + x118) + x128 * x94
    x143 = x0 * (3.0 * x115 + 3.0 * x116 + x128) + x131 * x94
    x144 = x140 * x143
    x145 = x140 * (x0 * (3.0 * x120 + 3.0 * x121 + 2.0 * x129 + 2.0 * x130) + x132 * x94)
    x146 = -x61 - C[1]
    x147 = x60 * (x0 * (3.0 * x51 + 3.0 * x52 + 2.0 * x66 + 2.0 * x67) + x3 * x54)
    x148 = x146 * x99
    x149 = x148 + x75
    x150 = x0 * (3.0 * x20 + 3.0 * x23 + x78) + x3 * x68
    x151 = x146 * x73
    x152 = x0 * (x151 + x99)
    x153 = x149 * x62
    x154 = x152 + x153
    x155 = x0 * (x31 + 3.0 * x63) + x3 * x78
    x156 = x146 * x87
    x157 = x156 + x75
    x158 = x149 * x85
    x159 = x152 + x158
    x160 = x0 * (x108 + 2.0 * x148 + x74)
    x161 = x154 * x85
    x162 = x160 + x161
    x163 = x0 * (x151 + x87)
    x164 = x157 * x85
    x165 = x163 + x164
    x166 = x0 * (x108 + x148 + x156 + x88)
    x167 = x159 * x85
    x168 = x166 + x167
    x169 = 2.0 * x158
    x170 = x0 * (3.0 * x152 + x153 + x169 + x92)
    x171 = x162 * x85
    x172 = x170 + x171
    x173 = x0 * (x103 + x108 + 2.0 * x156) + x165 * x85
    x174 = x0 * (x107 + 2.0 * x152 + x165 + x169)
    x175 = x168 * x85
    x176 = x174 + x175
    x177 = x0 * (x112 + 2.0 * x160 + 2.0 * x161 + 2.0 * x166 + 2.0 * x167) + x172 * x85
    x178 = x18 * x59
    x179 = x139 * x18
    x180 = x0 * (x123 + 3.0 * x163 + 3.0 * x164) + x173 * x85
    x181 = x135 * (x0 * (x126 + 3.0 * x166 + 3.0 * x167 + x173) + x176 * x85)
    x182 = -x71 - C[2]
    x183 = x182 * x60
    x184 = x182 * x93
    x185 = x184 + x83
    x186 = x182 * x80
    x187 = x0 * (x186 + x93)
    x188 = x185 * x72
    x189 = x187 + x188
    x190 = x182 * x96
    x191 = x190 + x83
    x192 = x185 * x94
    x193 = x187 + x192
    x194 = x0 * (x118 + 2.0 * x184 + x82)
    x195 = x189 * x94
    x196 = x194 + x195
    x197 = x0 * (x186 + x96)
    x198 = x191 * x94
    x199 = x197 + x198
    x200 = x0 * (x118 + x184 + x190 + x97)
    x201 = x193 * x94
    x202 = x200 + x201
    x203 = 2.0 * x192
    x204 = x0 * (x102 + 3.0 * x187 + x188 + x203)
    x205 = x196 * x94
    x206 = x204 + x205
    x207 = x0 * (x113 + x118 + 2.0 * x190) + x199 * x94
    x208 = x0 * (x117 + 2.0 * x187 + x199 + x203)
    x209 = x202 * x94
    x210 = x208 + x209
    x211 = x0 * (x122 + 2.0 * x194 + 2.0 * x195 + 2.0 * x200 + 2.0 * x201) + x206 * x94
    x212 = x0 * (x128 + 3.0 * x197 + 3.0 * x198) + x207 * x94
    x213 = x140 * (x0 * (x131 + 3.0 * x200 + 3.0 * x201 + x207) + x210 * x94)

    # 270 item(s)
    return numpy.array(
        [
            x60 * (x0 * (2.0 * x30 + 2.0 * x36 + 3.0 * x43 + 3.0 * x49 + x54) + x3 * x56),
            x62 * x70,
            x70 * x72,
            x76 * x79 * x80,
            x62 * x79 * x81,
            x73 * x79 * x84,
            x85 * x86,
            x69 * x80 * x89,
            x69 * x81 * x85,
            x65 * x80 * x92,
            x65 * x89 * x93,
            x65 * x84 * x87,
            x86 * x94,
            x62 * x69 * x95,
            x69 * x73 * x98,
            x65 * x76 * x96,
            x65 * x98 * x99,
            x102 * x65 * x73,
            x104 * x55 * x80,
            x107 * x35 * x80,
            x104 * x35 * x93,
            x112 * x29 * x80,
            x107 * x29 * x93,
            x104 * x29 * x84,
            x55 * x85 * x95,
            x35 * x89 * x96,
            x35 * x87 * x98,
            x29 * x92 * x96,
            x29 * x89 * x98,
            x102 * x29 * x87,
            x114 * x55 * x73,
            x114 * x35 * x99,
            x117 * x35 * x73,
            x114 * x29 * x76,
            x117 * x29 * x99,
            x122 * x29 * x73,
            x123 * x48 * x80,
            x126 * x33 * x80,
            x123 * x33 * x93,
            x127 * x27 * x80,
            x126 * x27 * x93,
            x123 * x27 * x84,
            x104 * x48 * x96,
            x107 * x33 * x96,
            x104 * x33 * x98,
            x112 * x27 * x96,
            x107 * x27 * x98,
            x102 * x104 * x27,
            x114 * x48 * x87,
            x114 * x33 * x89,
            x117 * x33 * x87,
            x114 * x27 * x92,
            x117 * x27 * x89,
            x122 * x27 * x87,
            x128 * x48 * x73,
            x128 * x33 * x99,
            x131 * x33 * x73,
            x128 * x27 * x76,
            x131 * x27 * x99,
            x132 * x27 * x73,
            x133 * x46 * x80,
            x13 * x134 * x80,
            x13 * x133 * x93,
            x136 * x9,
            x134 * x137 * x9,
            x133 * x16 * x84,
            x123 * x46 * x96,
            x126 * x13 * x96,
            x123 * x13 * x98,
            x127 * x138 * x9,
            x126 * x16 * x98,
            x102 * x123 * x16,
            x104 * x114 * x46,
            x107 * x114 * x13,
            x104 * x117 * x13,
            x112 * x114 * x16,
            x107 * x117 * x16,
            x104 * x122 * x16,
            x128 * x46 * x87,
            x128 * x13 * x89,
            x13 * x131 * x87,
            x128 * x16 * x92,
            x131 * x16 * x89,
            x132 * x141 * x9,
            x142 * x46 * x73,
            x13 * x142 * x99,
            x13 * x143 * x73,
            x142 * x16 * x76,
            x144 * x62 * x9,
            x145 * x9,
            x146 * x147,
            x149 * x150 * x80,
            x146 * x150 * x81,
            x154 * x155 * x80,
            x149 * x155 * x93,
            x151 * x155 * x84,
            x157 * x54 * x80,
            x159 * x68 * x80,
            x157 * x68 * x93,
            x162 * x78 * x80,
            x159 * x78 * x93,
            x157 * x78 * x84,
            x146 * x54 * x95,
            x149 * x68 * x96,
            x151 * x68 * x98,
            x154 * x78 * x96,
            x149 * x78 * x98,
            x102 * x151 * x78,
            x165 * x53 * x80,
            x168 * x24 * x80,
            x165 * x24 * x93,
            x172 * x77 * x80,
            x168 * x77 * x93,
            x165 * x77 * x84,
            x157 * x53 * x96,
            x159 * x24 * x96,
            x157 * x24 * x98,
            x162 * x77 * x96,
            x159 * x77 * x98,
            x102 * x157 * x77,
            x114 * x151 * x53,
            x114 * x149 * x24,
            x117 * x151 * x24,
            x114 * x154 * x77,
            x117 * x149 * x77,
            x122 * x151 * x77,
            x173 * x42 * x80,
            x176 * x22 * x80,
            x173 * x22 * x93,
            x177 * x178,
            x176 * x178 * x72,
            x173 * x19 * x84,
            x165 * x42 * x96,
            x168 * x22 * x96,
            x165 * x22 * x98,
            x172 * x178 * x94,
            x168 * x19 * x98,
            x102 * x165 * x19,
            x114 * x157 * x42,
            x114 * x159 * x22,
            x117 * x157 * x22,
            x114 * x162 * x19,
            x117 * x159 * x19,
            x122 * x157 * x19,
            x128 * x151 * x42,
            x128 * x149 * x22,
            x131 * x151 * x22,
            x128 * x154 * x19,
            x131 * x149 * x19,
            x132 * x146 * x179,
            x180 * x39 * x80,
            x10 * x181,
            x10 * x137 * x180,
            x135
            * (
                x0 * (x127 + 3.0 * x170 + 3.0 * x171 + 2.0 * x174 + 2.0 * x175)
                + x177 * x85
            ),
            x181 * x72,
            x180 * x7 * x84,
            x173 * x39 * x96,
            x10 * x138 * x176,
            x11 * x173 * x98,
            x138 * x177,
            x176 * x7 * x98,
            x102 * x173 * x7,
            x114 * x165 * x39,
            x11 * x114 * x168,
            x11 * x117 * x165,
            x114 * x172 * x7,
            x117 * x168 * x7,
            x122 * x165 * x7,
            x128 * x157 * x39,
            x11 * x128 * x159,
            x11 * x131 * x157,
            x128 * x162 * x7,
            x131 * x159 * x7,
            x132 * x157 * x7,
            x142 * x151 * x39,
            x11 * x142 * x149,
            x10 * x144 * x146,
            x142 * x154 * x7,
            x143 * x149 * x7,
            x145 * x146,
            x147 * x182,
            x150 * x183 * x62,
            x150 * x185 * x73,
            x155 * x186 * x76,
            x155 * x185 * x99,
            x155 * x189 * x73,
            x183 * x54 * x85,
            x186 * x68 * x89,
            x185 * x68 * x87,
            x186 * x78 * x92,
            x185 * x78 * x89,
            x189 * x78 * x87,
            x191 * x54 * x73,
            x191 * x68 * x99,
            x193 * x68 * x73,
            x191 * x76 * x78,
            x193 * x78 * x99,
            x196 * x73 * x78,
            x104 * x186 * x53,
            x107 * x186 * x24,
            x104 * x185 * x24,
            x112 * x186 * x77,
            x107 * x185 * x77,
            x104 * x189 * x77,
            x191 * x53 * x87,
            x191 * x24 * x89,
            x193 * x24 * x87,
            x191 * x77 * x92,
            x193 * x77 * x89,
            x196 * x77 * x87,
            x199 * x53 * x73,
            x199 * x24 * x99,
            x202 * x24 * x73,
            x199 * x76 * x77,
            x202 * x77 * x99,
            x206 * x73 * x77,
            x123 * x186 * x42,
            x126 * x186 * x22,
            x123 * x185 * x22,
            x127 * x178 * x182,
            x126 * x185 * x19,
            x123 * x189 * x19,
            x104 * x191 * x42,
            x107 * x191 * x22,
            x104 * x193 * x22,
            x112 * x19 * x191,
            x107 * x19 * x193,
            x104 * x19 * x196,
            x199 * x42 * x87,
            x199 * x22 * x89,
            x202 * x22 * x87,
            x19 * x199 * x92,
            x19 * x202 * x89,
            x179 * x206 * x85,
            x207 * x42 * x73,
            x207 * x22 * x99,
            x210 * x22 * x73,
            x19 * x207 * x76,
            x179 * x210 * x62,
            x179 * x211,
            x133 * x186 * x39,
            x10 * x134 * x135 * x182,
            x11 * x133 * x185,
            x136 * x182,
            x134 * x185 * x7,
            x133 * x189 * x7,
            x123 * x191 * x39,
            x11 * x126 * x191,
            x11 * x123 * x193,
            x127 * x191 * x7,
            x126 * x193 * x7,
            x123 * x196 * x7,
            x104 * x199 * x39,
            x107 * x11 * x199,
            x104 * x11 * x202,
            x112 * x199 * x7,
            x107 * x202 * x7,
            x104 * x206 * x7,
            x207 * x39 * x87,
            x11 * x207 * x89,
            x10 * x141 * x210,
            x207 * x7 * x92,
            x210 * x7 * x89,
            x141 * x211,
            x212 * x39 * x73,
            x10 * x140 * x212 * x62,
            x10 * x213,
            x212 * x7 * x76,
            x213 * x62,
            x140
            * (
                x0 * (x132 + 3.0 * x204 + 3.0 * x205 + 2.0 * x208 + 2.0 * x209)
                + x211 * x94
            ),
        ]
    )


def dipole3d_43(a, A, b, B, C):
    """Cartesian 3D (gf) dipole moment integrals.
    The origin is at C.

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
    x18 = -x2 - A[0]
    x19 = x16 * x18
    x20 = 3.0 * x12
    x21 = x3**2 * x7
    x22 = x20 + x21
    x23 = x0 * (2.0 * x13 + x22)
    x24 = x0 * (x20 + 3.0 * x21)
    x25 = x12 + x21
    x26 = x25 * x3
    x27 = x12 * x3
    x28 = 2.0 * x27
    x29 = x26 + x28
    x30 = x18 * x29
    x31 = x24 + x30
    x32 = x0 * (x17 + 3.0 * x19 + 4.0 * x23 + x31)
    x33 = x18 * x8
    x34 = x10 * x18
    x35 = x0 * (x13 + x20 + x33 + x34)
    x36 = x14 * x18
    x37 = x11 + x36
    x38 = x18 * x37
    x39 = 2.0 * x33
    x40 = x0 * (x22 + x39)
    x41 = x18 * x25
    x42 = x28 + x41
    x43 = x18 * x42
    x44 = x40 + x43
    x45 = x0 * (2.0 * x19 + 2.0 * x23 + 2.0 * x35 + 2.0 * x38 + x44)
    x46 = 3.0 * x11
    x47 = x0 * (3.0 * x15 + x29 + x46)
    x48 = x17 + x23
    x49 = x18 * x48
    x50 = x47 + x49
    x51 = x18 * x50
    x52 = 2.0 * x36
    x53 = x0 * (x15 + x42 + x46 + x52)
    x54 = x19 + x23
    x55 = x18 * x54
    x56 = x53 + x55
    x57 = x18 * x56
    x58 = 3.0 * x40 + 3.0 * x43
    x59 = x0 * (x26 + 8.0 * x27 + 3.0 * x41)
    x60 = x18 * x31
    x61 = x59 + x60
    x62 = x0 * (2.0 * x24 + 2.0 * x30 + x58) + x18 * x61
    x63 = 3.0 * x53 + 3.0 * x55
    x64 = x32 + x51
    x65 = x0 * (2.0 * x47 + 2.0 * x49 + x61 + x63) + x18 * x64
    x66 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x67 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x68 = 3.14159265358979 * x1 * x67
    x69 = x66 * x68
    x70 = -x1 * (a * A[1] + b * B[1])
    x71 = -x70 - B[1]
    x72 = x18 * x7
    x73 = x0 * (x72 + x8)
    x74 = x12 + x33
    x75 = x18 * x74
    x76 = x73 + x75
    x77 = x0 * (x10 + x72)
    x78 = x12 + x34
    x79 = x18 * x78
    x80 = x77 + x79
    x81 = x0 * (2.0 * x11 + x52 + x76 + x80)
    x82 = x35 + x38
    x83 = x18 * x82
    x84 = x0 * (4.0 * x27 + 2.0 * x41 + 2.0 * x73 + 2.0 * x75)
    x85 = x18 * x44
    x86 = x84 + x85
    x87 = x45 + x57
    x88 = x69 * (x0 * (x63 + 2.0 * x81 + 2.0 * x83 + x86) + x18 * x87)
    x89 = -x1 * (a * A[2] + b * B[2])
    x90 = -x89 - B[2]
    x91 = x6 * x66
    x92 = x71**2 * x91
    x93 = x0 * x91
    x94 = x92 + x93
    x95 = x18**2 * x7
    x96 = x20 + x95
    x97 = x0 * (2.0 * x34 + x96) + x18 * x80
    x98 = x0 * (x39 + x96)
    x99 = x18 * x76
    x100 = x98 + x99
    x101 = x81 + x83
    x102 = x0 * (x100 + 3.0 * x35 + 3.0 * x38 + x97) + x101 * x18
    x103 = x6 * x67
    x104 = x69 * x90
    x105 = x103 * x90**2
    x106 = x0 * x103
    x107 = x105 + x106
    x108 = x71 * x94
    x109 = x71 * x93
    x110 = 2.0 * x109
    x111 = x108 + x110
    x112 = x12 + x95
    x113 = x112 * x18 + 2.0 * x12 * x18
    x114 = x0 * (x113 + 3.0 * x77 + 3.0 * x79) + x18 * x97
    x115 = x103 * x90
    x116 = x71 * x91
    x117 = x107 * x90
    x118 = x106 * x90
    x119 = 2.0 * x118
    x120 = x117 + x119
    x121 = -x70 - A[1]
    x122 = x65 * x69
    x123 = x121 * x91
    x124 = x123 * x71
    x125 = x124 + x93
    x126 = x121 * x94
    x127 = x110 + x126
    x128 = 3.0 * x93
    x129 = x0 * (x128 + 3.0 * x92)
    x130 = x111 * x121
    x131 = x129 + x130
    x132 = -x89 - A[2]
    x133 = x132 * x69
    x134 = x103 * x132
    x135 = x134 * x90
    x136 = x106 + x135
    x137 = x107 * x132
    x138 = x119 + x137
    x139 = 3.0 * x106
    x140 = x0 * (3.0 * x105 + x139)
    x141 = x120 * x132
    x142 = x140 + x141
    x143 = x121**2 * x91
    x144 = x143 + x93
    x145 = x0 * (x116 + x123)
    x146 = x121 * x125
    x147 = x145 + x146
    x148 = 2.0 * x124 + x128
    x149 = x0 * (x148 + x92)
    x150 = x121 * x127
    x151 = x149 + x150
    x152 = x0 * (x108 + 8.0 * x109 + 3.0 * x126)
    x153 = x121 * x131
    x154 = x152 + x153
    x155 = x103 * x132**2
    x156 = x106 + x155
    x157 = x0 * (x115 + x134)
    x158 = x132 * x136
    x159 = x157 + x158
    x160 = 2.0 * x135 + x139
    x161 = x0 * (x105 + x160)
    x162 = x132 * x138
    x163 = x161 + x162
    x164 = x0 * (x117 + 8.0 * x118 + 3.0 * x137)
    x165 = x132 * x142
    x166 = x164 + x165
    x167 = x121 * x144 + 2.0 * x121 * x93
    x168 = x0 * (x143 + x148)
    x169 = x121 * x147
    x170 = x168 + x169
    x171 = x0 * (4.0 * x109 + 2.0 * x126 + 2.0 * x145 + 2.0 * x146)
    x172 = x121 * x151
    x173 = x171 + x172
    x174 = 3.0 * x149 + 3.0 * x150
    x175 = x0 * (2.0 * x129 + 2.0 * x130 + x174) + x121 * x154
    x176 = 2.0 * x106 * x132 + x132 * x156
    x177 = x0 * (x155 + x160)
    x178 = x132 * x159
    x179 = x177 + x178
    x180 = x0 * (4.0 * x118 + 2.0 * x137 + 2.0 * x157 + 2.0 * x158)
    x181 = x132 * x163
    x182 = x180 + x181
    x183 = 3.0 * x161 + 3.0 * x162
    x184 = x0 * (2.0 * x140 + 2.0 * x141 + x183) + x132 * x166
    x185 = x0 * (x128 + 3.0 * x143) + x121 * x167
    x186 = x0 * (3.0 * x145 + 3.0 * x146 + x167) + x121 * x170
    x187 = x0 * (2.0 * x168 + 2.0 * x169 + x174) + x121 * x173
    x188 = x5 * x68
    x189 = x188 * (x0 * (3.0 * x152 + 3.0 * x153 + 3.0 * x171 + 3.0 * x172) + x121 * x175)
    x190 = x188 * x9
    x191 = 3.14159265358979 * x1 * x5 * x66
    x192 = x191 * x9
    x193 = x0 * (x139 + 3.0 * x155) + x132 * x176
    x194 = x0 * (3.0 * x157 + 3.0 * x158 + x176) + x132 * x179
    x195 = x0 * (2.0 * x177 + 2.0 * x178 + x183) + x132 * x182
    x196 = x191 * (x0 * (3.0 * x164 + 3.0 * x165 + 3.0 * x180 + 3.0 * x181) + x132 * x184)
    x197 = -x70 - C[1]
    x198 = x69 * (x0 * (3.0 * x59 + 3.0 * x60 + 3.0 * x84 + 3.0 * x85) + x18 * x62)
    x199 = x116 * x197
    x200 = x199 + x93
    x201 = x0 * (x58 + 2.0 * x98 + 2.0 * x99) + x18 * x86
    x202 = x197 * x91
    x203 = x0 * (x116 + x202)
    x204 = x200 * x71
    x205 = x203 + x204
    x206 = x0 * (x113 + 3.0 * x73 + 3.0 * x75) + x100 * x18
    x207 = x0 * (x128 + 2.0 * x199 + x92)
    x208 = x205 * x71
    x209 = x207 + x208
    x210 = x0 * (x20 + 3.0 * x95) + x113 * x18
    x211 = x123 * x197
    x212 = x211 + x93
    x213 = x121 * x200
    x214 = x203 + x213
    x215 = x121 * x205
    x216 = x207 + x215
    x217 = 3.0 * x203
    x218 = x0 * (x111 + 3.0 * x204 + x217)
    x219 = x121 * x209
    x220 = x218 + x219
    x221 = x0 * (x123 + x202)
    x222 = x121 * x212
    x223 = x221 + x222
    x224 = x0 * (x124 + x128 + x199 + x211)
    x225 = x121 * x214
    x226 = x224 + x225
    x227 = 2.0 * x213
    x228 = x0 * (x127 + x204 + x217 + x227)
    x229 = x121 * x216
    x230 = x228 + x229
    x231 = x0 * (x131 + 4.0 * x207 + x208 + 3.0 * x215)
    x232 = x121 * x220
    x233 = x231 + x232
    x234 = x0 * (x128 + x143 + 2.0 * x211) + x121 * x223
    x235 = x0 * (x147 + 2.0 * x203 + x223 + x227)
    x236 = x121 * x226
    x237 = x235 + x236
    x238 = x0 * (x151 + 2.0 * x207 + 2.0 * x215 + 2.0 * x224 + 2.0 * x225)
    x239 = x121 * x230
    x240 = x238 + x239
    x241 = 3.0 * x228 + 3.0 * x229
    x242 = x0 * (x154 + 2.0 * x218 + 2.0 * x219 + x241) + x121 * x233
    x243 = x188 * x242
    x244 = x18 * x188
    x245 = x191 * x197
    x246 = x0 * (x167 + 3.0 * x221 + 3.0 * x222) + x121 * x234
    x247 = x0 * (x170 + 3.0 * x224 + 3.0 * x225 + x234) + x121 * x237
    x248 = x188 * (x0 * (x173 + 2.0 * x235 + 2.0 * x236 + x241) + x121 * x240)
    x249 = x188 * x3
    x250 = -x89 - C[2]
    x251 = x250 * x69
    x252 = x115 * x250
    x253 = x106 + x252
    x254 = x103 * x250
    x255 = x0 * (x115 + x254)
    x256 = x253 * x90
    x257 = x255 + x256
    x258 = x0 * (x105 + x139 + 2.0 * x252)
    x259 = x257 * x90
    x260 = x258 + x259
    x261 = x134 * x250
    x262 = x106 + x261
    x263 = x132 * x253
    x264 = x255 + x263
    x265 = x132 * x257
    x266 = x258 + x265
    x267 = 3.0 * x255
    x268 = x0 * (x120 + 3.0 * x256 + x267)
    x269 = x132 * x260
    x270 = x268 + x269
    x271 = x0 * (x134 + x254)
    x272 = x132 * x262
    x273 = x271 + x272
    x274 = x0 * (x135 + x139 + x252 + x261)
    x275 = x132 * x264
    x276 = x274 + x275
    x277 = 2.0 * x263
    x278 = x0 * (x138 + x256 + x267 + x277)
    x279 = x132 * x266
    x280 = x278 + x279
    x281 = x0 * (x142 + 4.0 * x258 + x259 + 3.0 * x265)
    x282 = x132 * x270
    x283 = x281 + x282
    x284 = x18 * x191
    x285 = x0 * (x139 + x155 + 2.0 * x261) + x132 * x273
    x286 = x0 * (x159 + 2.0 * x255 + x273 + x277)
    x287 = x132 * x276
    x288 = x286 + x287
    x289 = x0 * (x163 + 2.0 * x258 + 2.0 * x265 + 2.0 * x274 + 2.0 * x275)
    x290 = x132 * x280
    x291 = x289 + x290
    x292 = 3.0 * x278 + 3.0 * x279
    x293 = x0 * (x166 + 2.0 * x268 + 2.0 * x269 + x292) + x132 * x283
    x294 = x191 * x293
    x295 = x191 * x3
    x296 = x0 * (x176 + 3.0 * x271 + 3.0 * x272) + x132 * x285
    x297 = x0 * (x179 + 3.0 * x274 + 3.0 * x275 + x285) + x132 * x288
    x298 = x191 * (x0 * (x182 + 2.0 * x286 + 2.0 * x287 + x292) + x132 * x291)

    # 450 item(s)
    return numpy.array(
        [
            x69
            * (x0 * (3.0 * x32 + 3.0 * x45 + 3.0 * x51 + 3.0 * x57 + x62) + x18 * x65),
            x71 * x88,
            x88 * x90,
            x102 * x103 * x94,
            x102 * x104 * x71,
            x102 * x107 * x91,
            x103 * x111 * x114,
            x114 * x115 * x94,
            x107 * x114 * x116,
            x114 * x120 * x91,
            x121 * x122,
            x103 * x125 * x87,
            x104 * x121 * x87,
            x101 * x103 * x127,
            x101 * x115 * x125,
            x101 * x107 * x123,
            x103 * x131 * x97,
            x115 * x127 * x97,
            x107 * x125 * x97,
            x120 * x123 * x97,
            x122 * x132,
            x133 * x71 * x87,
            x136 * x87 * x91,
            x101 * x134 * x94,
            x101 * x116 * x136,
            x101 * x138 * x91,
            x111 * x134 * x97,
            x136 * x94 * x97,
            x116 * x138 * x97,
            x142 * x91 * x97,
            x103 * x144 * x64,
            x103 * x147 * x56,
            x115 * x144 * x56,
            x103 * x151 * x82,
            x115 * x147 * x82,
            x107 * x144 * x82,
            x103 * x154 * x80,
            x115 * x151 * x80,
            x107 * x147 * x80,
            x120 * x144 * x80,
            x121 * x133 * x64,
            x125 * x134 * x56,
            x123 * x136 * x56,
            x127 * x134 * x82,
            x125 * x136 * x82,
            x123 * x138 * x82,
            x131 * x134 * x80,
            x127 * x136 * x80,
            x125 * x138 * x80,
            x123 * x142 * x80,
            x156 * x64 * x91,
            x116 * x156 * x56,
            x159 * x56 * x91,
            x156 * x82 * x94,
            x116 * x159 * x82,
            x163 * x82 * x91,
            x111 * x156 * x80,
            x159 * x80 * x94,
            x116 * x163 * x80,
            x166 * x80 * x91,
            x103 * x167 * x50,
            x103 * x170 * x54,
            x115 * x167 * x54,
            x103 * x173 * x37,
            x115 * x170 * x37,
            x107 * x167 * x37,
            x103 * x175 * x78,
            x115 * x173 * x78,
            x107 * x170 * x78,
            x120 * x167 * x78,
            x134 * x144 * x50,
            x134 * x147 * x54,
            x136 * x144 * x54,
            x134 * x151 * x37,
            x136 * x147 * x37,
            x138 * x144 * x37,
            x134 * x154 * x78,
            x136 * x151 * x78,
            x138 * x147 * x78,
            x142 * x144 * x78,
            x123 * x156 * x50,
            x125 * x156 * x54,
            x123 * x159 * x54,
            x127 * x156 * x37,
            x125 * x159 * x37,
            x123 * x163 * x37,
            x131 * x156 * x78,
            x127 * x159 * x78,
            x125 * x163 * x78,
            x123 * x166 * x78,
            x176 * x50 * x91,
            x116 * x176 * x54,
            x179 * x54 * x91,
            x176 * x37 * x94,
            x116 * x179 * x37,
            x182 * x37 * x91,
            x111 * x176 * x78,
            x179 * x78 * x94,
            x116 * x182 * x78,
            x184 * x78 * x91,
            x103 * x185 * x48,
            x103 * x16 * x186,
            x115 * x16 * x185,
            x103 * x14 * x187,
            x115 * x14 * x186,
            x107 * x14 * x185,
            x189 * x9,
            x187 * x190 * x90,
            x10 * x107 * x186,
            x10 * x120 * x185,
            x134 * x167 * x48,
            x134 * x16 * x170,
            x136 * x16 * x167,
            x134 * x14 * x173,
            x136 * x14 * x170,
            x138 * x14 * x167,
            x132 * x175 * x190,
            x10 * x136 * x173,
            x10 * x138 * x170,
            x10 * x142 * x167,
            x144 * x156 * x48,
            x147 * x156 * x16,
            x144 * x159 * x16,
            x14 * x151 * x156,
            x14 * x147 * x159,
            x14 * x144 * x163,
            x10 * x154 * x156,
            x10 * x151 * x159,
            x10 * x147 * x163,
            x10 * x144 * x166,
            x123 * x176 * x48,
            x125 * x16 * x176,
            x123 * x16 * x179,
            x127 * x14 * x176,
            x125 * x14 * x179,
            x123 * x14 * x182,
            x10 * x131 * x176,
            x10 * x127 * x179,
            x10 * x125 * x182,
            x121 * x184 * x192,
            x193 * x48 * x91,
            x116 * x16 * x193,
            x16 * x194 * x91,
            x14 * x193 * x94,
            x116 * x14 * x194,
            x14 * x195 * x91,
            x10 * x111 * x193,
            x10 * x194 * x94,
            x192 * x195 * x71,
            x196 * x9,
            x197 * x198,
            x103 * x200 * x201,
            x104 * x197 * x201,
            x103 * x205 * x206,
            x115 * x200 * x206,
            x107 * x202 * x206,
            x103 * x209 * x210,
            x115 * x205 * x210,
            x107 * x200 * x210,
            x120 * x202 * x210,
            x103 * x212 * x62,
            x103 * x214 * x86,
            x115 * x212 * x86,
            x100 * x103 * x216,
            x100 * x115 * x214,
            x100 * x107 * x212,
            x103 * x113 * x220,
            x113 * x115 * x216,
            x107 * x113 * x214,
            x113 * x120 * x212,
            x133 * x197 * x62,
            x134 * x200 * x86,
            x136 * x202 * x86,
            x100 * x134 * x205,
            x100 * x136 * x200,
            x100 * x138 * x202,
            x113 * x134 * x209,
            x113 * x136 * x205,
            x113 * x138 * x200,
            x113 * x142 * x202,
            x103 * x223 * x61,
            x103 * x226 * x44,
            x115 * x223 * x44,
            x103 * x230 * x76,
            x115 * x226 * x76,
            x107 * x223 * x76,
            x103 * x112 * x233,
            x112 * x115 * x230,
            x107 * x112 * x226,
            x112 * x120 * x223,
            x134 * x212 * x61,
            x134 * x214 * x44,
            x136 * x212 * x44,
            x134 * x216 * x76,
            x136 * x214 * x76,
            x138 * x212 * x76,
            x112 * x134 * x220,
            x112 * x136 * x216,
            x112 * x138 * x214,
            x112 * x142 * x212,
            x156 * x202 * x61,
            x156 * x200 * x44,
            x159 * x202 * x44,
            x156 * x205 * x76,
            x159 * x200 * x76,
            x163 * x202 * x76,
            x112 * x156 * x209,
            x112 * x159 * x205,
            x112 * x163 * x200,
            x112 * x166 * x202,
            x103 * x234 * x31,
            x103 * x237 * x42,
            x115 * x234 * x42,
            x103 * x240 * x74,
            x115 * x237 * x74,
            x107 * x234 * x74,
            x18 * x243,
            x240 * x244 * x90,
            x107 * x237 * x72,
            x120 * x234 * x72,
            x134 * x223 * x31,
            x134 * x226 * x42,
            x136 * x223 * x42,
            x134 * x230 * x74,
            x136 * x226 * x74,
            x138 * x223 * x74,
            x132 * x233 * x244,
            x136 * x230 * x72,
            x138 * x226 * x72,
            x142 * x223 * x72,
            x156 * x212 * x31,
            x156 * x214 * x42,
            x159 * x212 * x42,
            x156 * x216 * x74,
            x159 * x214 * x74,
            x163 * x212 * x74,
            x156 * x220 * x72,
            x159 * x216 * x72,
            x163 * x214 * x72,
            x166 * x212 * x72,
            x176 * x202 * x31,
            x176 * x200 * x42,
            x179 * x202 * x42,
            x176 * x205 * x74,
            x179 * x200 * x74,
            x182 * x202 * x74,
            x176 * x209 * x72,
            x179 * x205 * x72,
            x182 * x200 * x72,
            x18 * x184 * x245,
            x103 * x246 * x29,
            x103 * x247 * x25,
            x115 * x246 * x25,
            x248 * x3,
            x247 * x249 * x90,
            x107 * x246 * x8,
            x188
            * (
                x0 * (x175 + 3.0 * x231 + 3.0 * x232 + 3.0 * x238 + 3.0 * x239)
                + x121 * x242
            ),
            x248 * x90,
            x107 * x247 * x7,
            x120 * x246 * x7,
            x134 * x234 * x29,
            x134 * x237 * x25,
            x136 * x234 * x25,
            x132 * x240 * x249,
            x136 * x237 * x8,
            x138 * x234 * x8,
            x132 * x243,
            x136 * x240 * x7,
            x138 * x237 * x7,
            x142 * x234 * x7,
            x156 * x223 * x29,
            x156 * x226 * x25,
            x159 * x223 * x25,
            x156 * x230 * x8,
            x159 * x226 * x8,
            x163 * x223 * x8,
            x156 * x233 * x7,
            x159 * x230 * x7,
            x163 * x226 * x7,
            x166 * x223 * x7,
            x176 * x212 * x29,
            x176 * x214 * x25,
            x179 * x212 * x25,
            x176 * x216 * x8,
            x179 * x214 * x8,
            x182 * x212 * x8,
            x176 * x220 * x7,
            x179 * x216 * x7,
            x182 * x214 * x7,
            x184 * x212 * x7,
            x193 * x202 * x29,
            x193 * x200 * x25,
            x194 * x202 * x25,
            x193 * x205 * x8,
            x194 * x200 * x8,
            x195 * x245 * x3,
            x193 * x209 * x7,
            x194 * x205 * x7,
            x195 * x200 * x7,
            x196 * x197,
            x198 * x250,
            x201 * x251 * x71,
            x201 * x253 * x91,
            x206 * x254 * x94,
            x116 * x206 * x253,
            x206 * x257 * x91,
            x111 * x210 * x254,
            x210 * x253 * x94,
            x116 * x210 * x257,
            x210 * x260 * x91,
            x121 * x251 * x62,
            x125 * x254 * x86,
            x123 * x253 * x86,
            x100 * x127 * x254,
            x100 * x125 * x253,
            x100 * x123 * x257,
            x113 * x131 * x254,
            x113 * x127 * x253,
            x113 * x125 * x257,
            x113 * x123 * x260,
            x262 * x62 * x91,
            x116 * x262 * x86,
            x264 * x86 * x91,
            x100 * x262 * x94,
            x100 * x116 * x264,
            x100 * x266 * x91,
            x111 * x113 * x262,
            x113 * x264 * x94,
            x113 * x116 * x266,
            x113 * x270 * x91,
            x144 * x254 * x61,
            x147 * x254 * x44,
            x144 * x253 * x44,
            x151 * x254 * x76,
            x147 * x253 * x76,
            x144 * x257 * x76,
            x112 * x154 * x254,
            x112 * x151 * x253,
            x112 * x147 * x257,
            x112 * x144 * x260,
            x123 * x262 * x61,
            x125 * x262 * x44,
            x123 * x264 * x44,
            x127 * x262 * x76,
            x125 * x264 * x76,
            x123 * x266 * x76,
            x112 * x131 * x262,
            x112 * x127 * x264,
            x112 * x125 * x266,
            x112 * x123 * x270,
            x273 * x61 * x91,
            x116 * x273 * x44,
            x276 * x44 * x91,
            x273 * x76 * x94,
            x116 * x276 * x76,
            x280 * x76 * x91,
            x111 * x112 * x273,
            x112 * x276 * x94,
            x112 * x116 * x280,
            x112 * x283 * x91,
            x167 * x254 * x31,
            x170 * x254 * x42,
            x167 * x253 * x42,
            x173 * x254 * x74,
            x170 * x253 * x74,
            x167 * x257 * x74,
            x175 * x244 * x250,
            x173 * x253 * x72,
            x170 * x257 * x72,
            x167 * x260 * x72,
            x144 * x262 * x31,
            x147 * x262 * x42,
            x144 * x264 * x42,
            x151 * x262 * x74,
            x147 * x264 * x74,
            x144 * x266 * x74,
            x154 * x262 * x72,
            x151 * x264 * x72,
            x147 * x266 * x72,
            x144 * x270 * x72,
            x123 * x273 * x31,
            x125 * x273 * x42,
            x123 * x276 * x42,
            x127 * x273 * x74,
            x125 * x276 * x74,
            x123 * x280 * x74,
            x131 * x273 * x72,
            x127 * x276 * x72,
            x125 * x280 * x72,
            x121 * x283 * x284,
            x285 * x31 * x91,
            x116 * x285 * x42,
            x288 * x42 * x91,
            x285 * x74 * x94,
            x116 * x288 * x74,
            x291 * x74 * x91,
            x111 * x285 * x72,
            x288 * x72 * x94,
            x284 * x291 * x71,
            x18 * x294,
            x185 * x254 * x29,
            x186 * x25 * x254,
            x185 * x25 * x253,
            x187 * x249 * x250,
            x186 * x253 * x8,
            x185 * x257 * x8,
            x189 * x250,
            x187 * x253 * x7,
            x186 * x257 * x7,
            x185 * x260 * x7,
            x167 * x262 * x29,
            x170 * x25 * x262,
            x167 * x25 * x264,
            x173 * x262 * x8,
            x170 * x264 * x8,
            x167 * x266 * x8,
            x175 * x262 * x7,
            x173 * x264 * x7,
            x170 * x266 * x7,
            x167 * x270 * x7,
            x144 * x273 * x29,
            x147 * x25 * x273,
            x144 * x25 * x276,
            x151 * x273 * x8,
            x147 * x276 * x8,
            x144 * x280 * x8,
            x154 * x273 * x7,
            x151 * x276 * x7,
            x147 * x280 * x7,
            x144 * x283 * x7,
            x123 * x285 * x29,
            x125 * x25 * x285,
            x123 * x25 * x288,
            x127 * x285 * x8,
            x125 * x288 * x8,
            x121 * x291 * x295,
            x131 * x285 * x7,
            x127 * x288 * x7,
            x125 * x291 * x7,
            x121 * x294,
            x29 * x296 * x91,
            x116 * x25 * x296,
            x25 * x297 * x91,
            x296 * x8 * x94,
            x295 * x297 * x71,
            x298 * x3,
            x111 * x296 * x7,
            x297 * x7 * x94,
            x298 * x71,
            x191
            * (
                x0 * (x184 + 3.0 * x281 + 3.0 * x282 + 3.0 * x289 + 3.0 * x290)
                + x132 * x293
            ),
        ]
    )


def dipole3d_44(a, A, b, B, C):
    """Cartesian 3D (gg) dipole moment integrals.
    The origin is at C.

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
    x9 = x3 * x8
    x10 = x4 * x9
    x11 = x3**2 * x8
    x12 = x0 * x8
    x13 = 3.0 * x12
    x14 = x11 + x13
    x15 = x0 * (2.0 * x10 + x14)
    x16 = x4 * x8
    x17 = x0 * (x16 + x9)
    x18 = x10 + x12
    x19 = x18 * x3
    x20 = x17 + x19
    x21 = x20 * x3
    x22 = x15 + x21
    x23 = x22 * x3
    x24 = -x2 - A[0]
    x25 = x22 * x24
    x26 = 3.0 * x17
    x27 = x11 + x12
    x28 = x27 * x3
    x29 = x12 * x3
    x30 = 2.0 * x29
    x31 = x28 + x30
    x32 = x0 * (3.0 * x19 + x26 + x31)
    x33 = 8.0 * x29
    x34 = x0 * (4.0 * x28 + x33)
    x35 = x0 * (3.0 * x11 + x13)
    x36 = x3 * x31
    x37 = x35 + x36
    x38 = x24 * x37
    x39 = x34 + x38
    x40 = x0 * (x23 + 4.0 * x25 + 5.0 * x32 + x39)
    x41 = 4.0 * x15
    x42 = x0 * (4.0 * x21 + x37 + x41)
    x43 = x23 + x32
    x44 = x24 * x43
    x45 = x42 + x44
    x46 = x24 * x45
    x47 = x24 * x27
    x48 = x0 * (x28 + x33 + 3.0 * x47)
    x49 = x24 * x31
    x50 = x35 + x49
    x51 = x24 * x50
    x52 = x48 + x51
    x53 = x18 * x24
    x54 = 2.0 * x53
    x55 = x30 + x47
    x56 = x0 * (x19 + x26 + x54 + x55)
    x57 = x20 * x24
    x58 = x15 + x57
    x59 = x24 * x58
    x60 = 3.0 * x56 + 3.0 * x59
    x61 = x0 * (2.0 * x25 + 2.0 * x32 + x52 + x60)
    x62 = x0 * (x21 + x41 + x50 + 3.0 * x57)
    x63 = x25 + x32
    x64 = x24 * x63
    x65 = x62 + x64
    x66 = x24 * x65
    x67 = x0 * (5.0 * x35 + x36 + 4.0 * x49)
    x68 = x24 * x39
    x69 = x67 + x68
    x70 = x0 * (2.0 * x34 + 2.0 * x38 + 4.0 * x48 + 4.0 * x51) + x24 * x69
    x71 = x40 + x46
    x72 = x0 * (2.0 * x42 + 2.0 * x44 + 4.0 * x62 + 4.0 * x64 + x69) + x24 * x71
    x73 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x74 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x75 = 3.14159265358979 * x1 * x74
    x76 = x73 * x75
    x77 = -x1 * (a * A[1] + b * B[1])
    x78 = -x77 - B[1]
    x79 = x24 * x9
    x80 = x16 * x24
    x81 = x0 * (x10 + x13 + x79 + x80)
    x82 = x17 + x53
    x83 = x24 * x82
    x84 = 2.0 * x79
    x85 = x0 * (x14 + x84)
    x86 = x24 * x55
    x87 = x85 + x86
    x88 = x0 * (2.0 * x15 + 2.0 * x57 + 2.0 * x81 + 2.0 * x83 + x87)
    x89 = x56 + x59
    x90 = x24 * x89
    x91 = 3.0 * x85 + 3.0 * x86
    x92 = x0 * (2.0 * x35 + 2.0 * x49 + x91)
    x93 = x24 * x52
    x94 = x92 + x93
    x95 = x61 + x66
    x96 = x76 * (x0 * (3.0 * x62 + 3.0 * x64 + 3.0 * x88 + 3.0 * x90 + x94) + x24 * x95)
    x97 = -x1 * (a * A[2] + b * B[2])
    x98 = -x97 - B[2]
    x99 = x7 * x73
    x100 = x78**2 * x99
    x101 = x0 * x99
    x102 = x100 + x101
    x103 = x24 * x8
    x104 = x0 * (x103 + x9)
    x105 = x12 + x79
    x106 = x105 * x24
    x107 = x104 + x106
    x108 = x0 * (x103 + x16)
    x109 = x12 + x80
    x110 = x109 * x24
    x111 = x108 + x110
    x112 = x0 * (x107 + x111 + 2.0 * x17 + x54)
    x113 = x81 + x83
    x114 = x113 * x24
    x115 = x0 * (2.0 * x104 + 2.0 * x106 + 4.0 * x29 + 2.0 * x47)
    x116 = x24 * x87
    x117 = x115 + x116
    x118 = x88 + x90
    x119 = x0 * (2.0 * x112 + 2.0 * x114 + x117 + x60) + x118 * x24
    x120 = x7 * x74
    x121 = x76 * x98
    x122 = x120 * x98**2
    x123 = x0 * x120
    x124 = x122 + x123
    x125 = x102 * x78
    x126 = x101 * x78
    x127 = 2.0 * x126
    x128 = x125 + x127
    x129 = x24**2 * x8
    x130 = x129 + x13
    x131 = x0 * (x130 + 2.0 * x80) + x111 * x24
    x132 = x0 * (x130 + x84)
    x133 = x107 * x24
    x134 = x132 + x133
    x135 = x112 + x114
    x136 = x0 * (x131 + x134 + 3.0 * x81 + 3.0 * x83) + x135 * x24
    x137 = x120 * x98
    x138 = x78 * x99
    x139 = x124 * x98
    x140 = x123 * x98
    x141 = 2.0 * x140
    x142 = x139 + x141
    x143 = 3.0 * x101
    x144 = x0 * (3.0 * x100 + x143)
    x145 = x128 * x78
    x146 = x144 + x145
    x147 = x12 + x129
    x148 = 2.0 * x12 * x24 + x147 * x24
    x149 = x0 * (3.0 * x108 + 3.0 * x110 + x148) + x131 * x24
    x150 = 3.0 * x123
    x151 = x0 * (3.0 * x122 + x150)
    x152 = x142 * x98
    x153 = x151 + x152
    x154 = -x77 - A[1]
    x155 = x72 * x76
    x156 = x154 * x99
    x157 = x156 * x78
    x158 = x101 + x157
    x159 = x102 * x154
    x160 = x127 + x159
    x161 = x128 * x154
    x162 = x144 + x161
    x163 = 8.0 * x126
    x164 = x0 * (4.0 * x125 + x163)
    x165 = x146 * x154
    x166 = x164 + x165
    x167 = -x97 - A[2]
    x168 = x167 * x76
    x169 = x120 * x167
    x170 = x169 * x98
    x171 = x123 + x170
    x172 = x124 * x167
    x173 = x141 + x172
    x174 = x142 * x167
    x175 = x151 + x174
    x176 = 8.0 * x140
    x177 = x0 * (4.0 * x139 + x176)
    x178 = x153 * x167
    x179 = x177 + x178
    x180 = x154**2 * x99
    x181 = x101 + x180
    x182 = x0 * (x138 + x156)
    x183 = x154 * x158
    x184 = x182 + x183
    x185 = x143 + 2.0 * x157
    x186 = x0 * (x100 + x185)
    x187 = x154 * x160
    x188 = x186 + x187
    x189 = x0 * (x125 + 3.0 * x159 + x163)
    x190 = x154 * x162
    x191 = x189 + x190
    x192 = x0 * (5.0 * x144 + x145 + 4.0 * x161)
    x193 = x154 * x166
    x194 = x192 + x193
    x195 = x120 * x167**2
    x196 = x123 + x195
    x197 = x0 * (x137 + x169)
    x198 = x167 * x171
    x199 = x197 + x198
    x200 = x150 + 2.0 * x170
    x201 = x0 * (x122 + x200)
    x202 = x167 * x173
    x203 = x201 + x202
    x204 = x0 * (x139 + 3.0 * x172 + x176)
    x205 = x167 * x175
    x206 = x204 + x205
    x207 = x0 * (5.0 * x151 + x152 + 4.0 * x174)
    x208 = x167 * x179
    x209 = x207 + x208
    x210 = 2.0 * x101 * x154 + x154 * x181
    x211 = x0 * (x180 + x185)
    x212 = x154 * x184
    x213 = x211 + x212
    x214 = x0 * (4.0 * x126 + 2.0 * x159 + 2.0 * x182 + 2.0 * x183)
    x215 = x154 * x188
    x216 = x214 + x215
    x217 = 3.0 * x186 + 3.0 * x187
    x218 = x0 * (2.0 * x144 + 2.0 * x161 + x217)
    x219 = x154 * x191
    x220 = x218 + x219
    x221 = x0 * (2.0 * x164 + 2.0 * x165 + 4.0 * x189 + 4.0 * x190) + x154 * x194
    x222 = 2.0 * x123 * x167 + x167 * x196
    x223 = x0 * (x195 + x200)
    x224 = x167 * x199
    x225 = x223 + x224
    x226 = x0 * (4.0 * x140 + 2.0 * x172 + 2.0 * x197 + 2.0 * x198)
    x227 = x167 * x203
    x228 = x226 + x227
    x229 = 3.0 * x201 + 3.0 * x202
    x230 = x0 * (2.0 * x151 + 2.0 * x174 + x229)
    x231 = x167 * x206
    x232 = x230 + x231
    x233 = x0 * (2.0 * x177 + 2.0 * x178 + 4.0 * x204 + 4.0 * x205) + x167 * x209
    x234 = x0 * (x143 + 3.0 * x180) + x154 * x210
    x235 = x0 * (3.0 * x182 + 3.0 * x183 + x210) + x154 * x213
    x236 = x0 * (2.0 * x211 + 2.0 * x212 + x217) + x154 * x216
    x237 = x0 * (3.0 * x189 + 3.0 * x190 + 3.0 * x214 + 3.0 * x215) + x154 * x220
    x238 = x6 * x75
    x239 = x238 * (x0 * (3.0 * x192 + 3.0 * x193 + 4.0 * x218 + 4.0 * x219) + x154 * x221)
    x240 = x238 * x4
    x241 = 3.14159265358979 * x1 * x6 * x73
    x242 = x241 * x4
    x243 = x0 * (x150 + 3.0 * x195) + x167 * x222
    x244 = x0 * (3.0 * x197 + 3.0 * x198 + x222) + x167 * x225
    x245 = x0 * (2.0 * x223 + 2.0 * x224 + x229) + x167 * x228
    x246 = x0 * (3.0 * x204 + 3.0 * x205 + 3.0 * x226 + 3.0 * x227) + x167 * x232
    x247 = x241 * (x0 * (3.0 * x207 + 3.0 * x208 + 4.0 * x230 + 4.0 * x231) + x167 * x233)
    x248 = -x77 - C[1]
    x249 = x76 * (x0 * (3.0 * x67 + 3.0 * x68 + 4.0 * x92 + 4.0 * x93) + x24 * x70)
    x250 = x138 * x248
    x251 = x101 + x250
    x252 = x0 * (3.0 * x115 + 3.0 * x116 + 3.0 * x48 + 3.0 * x51) + x24 * x94
    x253 = x248 * x99
    x254 = x0 * (x138 + x253)
    x255 = x251 * x78
    x256 = x254 + x255
    x257 = x0 * (2.0 * x132 + 2.0 * x133 + x91) + x117 * x24
    x258 = x0 * (x100 + x143 + 2.0 * x250)
    x259 = x256 * x78
    x260 = x258 + x259
    x261 = x0 * (3.0 * x104 + 3.0 * x106 + x148) + x134 * x24
    x262 = 3.0 * x254
    x263 = x0 * (x128 + 3.0 * x255 + x262)
    x264 = x260 * x78
    x265 = x263 + x264
    x266 = x0 * (3.0 * x129 + x13) + x148 * x24
    x267 = x156 * x248
    x268 = x101 + x267
    x269 = x154 * x251
    x270 = x254 + x269
    x271 = x154 * x256
    x272 = x258 + x271
    x273 = x154 * x260
    x274 = x263 + x273
    x275 = 4.0 * x258
    x276 = x0 * (x146 + 4.0 * x259 + x275)
    x277 = x154 * x265
    x278 = x276 + x277
    x279 = x0 * (x156 + x253)
    x280 = x154 * x268
    x281 = x279 + x280
    x282 = x0 * (x143 + x157 + x250 + x267)
    x283 = x154 * x270
    x284 = x282 + x283
    x285 = 2.0 * x269
    x286 = x0 * (x160 + x255 + x262 + x285)
    x287 = x154 * x272
    x288 = x286 + x287
    x289 = x0 * (x162 + x259 + 3.0 * x271 + x275)
    x290 = x154 * x274
    x291 = x289 + x290
    x292 = x0 * (x166 + 5.0 * x263 + x264 + 4.0 * x273)
    x293 = x154 * x278
    x294 = x292 + x293
    x295 = x0 * (x143 + x180 + 2.0 * x267) + x154 * x281
    x296 = x0 * (x184 + 2.0 * x254 + x281 + x285)
    x297 = x154 * x284
    x298 = x296 + x297
    x299 = x0 * (x188 + 2.0 * x258 + 2.0 * x271 + 2.0 * x282 + 2.0 * x283)
    x300 = x154 * x288
    x301 = x299 + x300
    x302 = 3.0 * x286 + 3.0 * x287
    x303 = x0 * (x191 + 2.0 * x263 + 2.0 * x273 + x302)
    x304 = x154 * x291
    x305 = x303 + x304
    x306 = x0 * (x194 + 2.0 * x276 + 2.0 * x277 + 4.0 * x289 + 4.0 * x290) + x154 * x294
    x307 = x238 * x306
    x308 = x238 * x24
    x309 = x241 * x248
    x310 = x0 * (x210 + 3.0 * x279 + 3.0 * x280) + x154 * x295
    x311 = x0 * (x213 + 3.0 * x282 + 3.0 * x283 + x295) + x154 * x298
    x312 = x0 * (x216 + 2.0 * x296 + 2.0 * x297 + x302) + x154 * x301
    x313 = x238 * (
        x0 * (x220 + 3.0 * x289 + 3.0 * x290 + 3.0 * x299 + 3.0 * x300) + x154 * x305
    )
    x314 = x238 * x3
    x315 = -x97 - C[2]
    x316 = x315 * x76
    x317 = x137 * x315
    x318 = x123 + x317
    x319 = x120 * x315
    x320 = x0 * (x137 + x319)
    x321 = x318 * x98
    x322 = x320 + x321
    x323 = x0 * (x122 + x150 + 2.0 * x317)
    x324 = x322 * x98
    x325 = x323 + x324
    x326 = 3.0 * x320
    x327 = x0 * (x142 + 3.0 * x321 + x326)
    x328 = x325 * x98
    x329 = x327 + x328
    x330 = x169 * x315
    x331 = x123 + x330
    x332 = x167 * x318
    x333 = x320 + x332
    x334 = x167 * x322
    x335 = x323 + x334
    x336 = x167 * x325
    x337 = x327 + x336
    x338 = 4.0 * x323
    x339 = x0 * (x153 + 4.0 * x324 + x338)
    x340 = x167 * x329
    x341 = x339 + x340
    x342 = x0 * (x169 + x319)
    x343 = x167 * x331
    x344 = x342 + x343
    x345 = x0 * (x150 + x170 + x317 + x330)
    x346 = x167 * x333
    x347 = x345 + x346
    x348 = 2.0 * x332
    x349 = x0 * (x173 + x321 + x326 + x348)
    x350 = x167 * x335
    x351 = x349 + x350
    x352 = x0 * (x175 + x324 + 3.0 * x334 + x338)
    x353 = x167 * x337
    x354 = x352 + x353
    x355 = x0 * (x179 + 5.0 * x327 + x328 + 4.0 * x336)
    x356 = x167 * x341
    x357 = x355 + x356
    x358 = x24 * x241
    x359 = x0 * (x150 + x195 + 2.0 * x330) + x167 * x344
    x360 = x0 * (x199 + 2.0 * x320 + x344 + x348)
    x361 = x167 * x347
    x362 = x360 + x361
    x363 = x0 * (x203 + 2.0 * x323 + 2.0 * x334 + 2.0 * x345 + 2.0 * x346)
    x364 = x167 * x351
    x365 = x363 + x364
    x366 = 3.0 * x349 + 3.0 * x350
    x367 = x0 * (x206 + 2.0 * x327 + 2.0 * x336 + x366)
    x368 = x167 * x354
    x369 = x367 + x368
    x370 = x0 * (x209 + 2.0 * x339 + 2.0 * x340 + 4.0 * x352 + 4.0 * x353) + x167 * x357
    x371 = x241 * x370
    x372 = x241 * x3
    x373 = x0 * (x222 + 3.0 * x342 + 3.0 * x343) + x167 * x359
    x374 = x0 * (x225 + 3.0 * x345 + 3.0 * x346 + x359) + x167 * x362
    x375 = x0 * (x228 + 2.0 * x360 + 2.0 * x361 + x366) + x167 * x365
    x376 = x241 * (
        x0 * (x232 + 3.0 * x352 + 3.0 * x353 + 3.0 * x363 + 3.0 * x364) + x167 * x369
    )

    # 675 item(s)
    return numpy.array(
        [
            x76
            * (x0 * (3.0 * x40 + 3.0 * x46 + 4.0 * x61 + 4.0 * x66 + x70) + x24 * x72),
            x78 * x96,
            x96 * x98,
            x102 * x119 * x120,
            x119 * x121 * x78,
            x119 * x124 * x99,
            x120 * x128 * x136,
            x102 * x136 * x137,
            x124 * x136 * x138,
            x136 * x142 * x99,
            x120 * x146 * x149,
            x128 * x137 * x149,
            x102 * x124 * x149,
            x138 * x142 * x149,
            x149 * x153 * x99,
            x154 * x155,
            x120 * x158 * x95,
            x121 * x154 * x95,
            x118 * x120 * x160,
            x118 * x137 * x158,
            x118 * x124 * x156,
            x120 * x135 * x162,
            x135 * x137 * x160,
            x124 * x135 * x158,
            x135 * x142 * x156,
            x120 * x131 * x166,
            x131 * x137 * x162,
            x124 * x131 * x160,
            x131 * x142 * x158,
            x131 * x153 * x156,
            x155 * x167,
            x168 * x78 * x95,
            x171 * x95 * x99,
            x102 * x118 * x169,
            x118 * x138 * x171,
            x118 * x173 * x99,
            x128 * x135 * x169,
            x102 * x135 * x171,
            x135 * x138 * x173,
            x135 * x175 * x99,
            x131 * x146 * x169,
            x128 * x131 * x171,
            x102 * x131 * x173,
            x131 * x138 * x175,
            x131 * x179 * x99,
            x120 * x181 * x71,
            x120 * x184 * x65,
            x137 * x181 * x65,
            x120 * x188 * x89,
            x137 * x184 * x89,
            x124 * x181 * x89,
            x113 * x120 * x191,
            x113 * x137 * x188,
            x113 * x124 * x184,
            x113 * x142 * x181,
            x111 * x120 * x194,
            x111 * x137 * x191,
            x111 * x124 * x188,
            x111 * x142 * x184,
            x111 * x153 * x181,
            x154 * x168 * x71,
            x158 * x169 * x65,
            x156 * x171 * x65,
            x160 * x169 * x89,
            x158 * x171 * x89,
            x156 * x173 * x89,
            x113 * x162 * x169,
            x113 * x160 * x171,
            x113 * x158 * x173,
            x113 * x156 * x175,
            x111 * x166 * x169,
            x111 * x162 * x171,
            x111 * x160 * x173,
            x111 * x158 * x175,
            x111 * x156 * x179,
            x196 * x71 * x99,
            x138 * x196 * x65,
            x199 * x65 * x99,
            x102 * x196 * x89,
            x138 * x199 * x89,
            x203 * x89 * x99,
            x113 * x128 * x196,
            x102 * x113 * x199,
            x113 * x138 * x203,
            x113 * x206 * x99,
            x111 * x146 * x196,
            x111 * x128 * x199,
            x102 * x111 * x203,
            x111 * x138 * x206,
            x111 * x209 * x99,
            x120 * x210 * x45,
            x120 * x213 * x63,
            x137 * x210 * x63,
            x120 * x216 * x58,
            x137 * x213 * x58,
            x124 * x210 * x58,
            x120 * x220 * x82,
            x137 * x216 * x82,
            x124 * x213 * x82,
            x142 * x210 * x82,
            x109 * x120 * x221,
            x109 * x137 * x220,
            x109 * x124 * x216,
            x109 * x142 * x213,
            x109 * x153 * x210,
            x169 * x181 * x45,
            x169 * x184 * x63,
            x171 * x181 * x63,
            x169 * x188 * x58,
            x171 * x184 * x58,
            x173 * x181 * x58,
            x169 * x191 * x82,
            x171 * x188 * x82,
            x173 * x184 * x82,
            x175 * x181 * x82,
            x109 * x169 * x194,
            x109 * x171 * x191,
            x109 * x173 * x188,
            x109 * x175 * x184,
            x109 * x179 * x181,
            x156 * x196 * x45,
            x158 * x196 * x63,
            x156 * x199 * x63,
            x160 * x196 * x58,
            x158 * x199 * x58,
            x156 * x203 * x58,
            x162 * x196 * x82,
            x160 * x199 * x82,
            x158 * x203 * x82,
            x156 * x206 * x82,
            x109 * x166 * x196,
            x109 * x162 * x199,
            x109 * x160 * x203,
            x109 * x158 * x206,
            x109 * x156 * x209,
            x222 * x45 * x99,
            x138 * x222 * x63,
            x225 * x63 * x99,
            x102 * x222 * x58,
            x138 * x225 * x58,
            x228 * x58 * x99,
            x128 * x222 * x82,
            x102 * x225 * x82,
            x138 * x228 * x82,
            x232 * x82 * x99,
            x109 * x146 * x222,
            x109 * x128 * x225,
            x102 * x109 * x228,
            x109 * x138 * x232,
            x109 * x233 * x99,
            x120 * x234 * x43,
            x120 * x22 * x235,
            x137 * x22 * x234,
            x120 * x20 * x236,
            x137 * x20 * x235,
            x124 * x20 * x234,
            x120 * x18 * x237,
            x137 * x18 * x236,
            x124 * x18 * x235,
            x142 * x18 * x234,
            x239 * x4,
            x237 * x240 * x98,
            x124 * x16 * x236,
            x142 * x16 * x235,
            x153 * x16 * x234,
            x169 * x210 * x43,
            x169 * x213 * x22,
            x171 * x210 * x22,
            x169 * x20 * x216,
            x171 * x20 * x213,
            x173 * x20 * x210,
            x169 * x18 * x220,
            x171 * x18 * x216,
            x173 * x18 * x213,
            x175 * x18 * x210,
            x167 * x221 * x240,
            x16 * x171 * x220,
            x16 * x173 * x216,
            x16 * x175 * x213,
            x16 * x179 * x210,
            x181 * x196 * x43,
            x184 * x196 * x22,
            x181 * x199 * x22,
            x188 * x196 * x20,
            x184 * x199 * x20,
            x181 * x20 * x203,
            x18 * x191 * x196,
            x18 * x188 * x199,
            x18 * x184 * x203,
            x18 * x181 * x206,
            x16 * x194 * x196,
            x16 * x191 * x199,
            x16 * x188 * x203,
            x16 * x184 * x206,
            x16 * x181 * x209,
            x156 * x222 * x43,
            x158 * x22 * x222,
            x156 * x22 * x225,
            x160 * x20 * x222,
            x158 * x20 * x225,
            x156 * x20 * x228,
            x162 * x18 * x222,
            x160 * x18 * x225,
            x158 * x18 * x228,
            x156 * x18 * x232,
            x16 * x166 * x222,
            x16 * x162 * x225,
            x16 * x160 * x228,
            x158 * x16 * x232,
            x154 * x233 * x242,
            x243 * x43 * x99,
            x138 * x22 * x243,
            x22 * x244 * x99,
            x102 * x20 * x243,
            x138 * x20 * x244,
            x20 * x245 * x99,
            x128 * x18 * x243,
            x102 * x18 * x244,
            x138 * x18 * x245,
            x18 * x246 * x99,
            x146 * x16 * x243,
            x128 * x16 * x244,
            x102 * x16 * x245,
            x242 * x246 * x78,
            x247 * x4,
            x248 * x249,
            x120 * x251 * x252,
            x121 * x248 * x252,
            x120 * x256 * x257,
            x137 * x251 * x257,
            x124 * x253 * x257,
            x120 * x260 * x261,
            x137 * x256 * x261,
            x124 * x251 * x261,
            x142 * x253 * x261,
            x120 * x265 * x266,
            x137 * x260 * x266,
            x124 * x256 * x266,
            x142 * x251 * x266,
            x153 * x253 * x266,
            x120 * x268 * x70,
            x120 * x270 * x94,
            x137 * x268 * x94,
            x117 * x120 * x272,
            x117 * x137 * x270,
            x117 * x124 * x268,
            x120 * x134 * x274,
            x134 * x137 * x272,
            x124 * x134 * x270,
            x134 * x142 * x268,
            x120 * x148 * x278,
            x137 * x148 * x274,
            x124 * x148 * x272,
            x142 * x148 * x270,
            x148 * x153 * x268,
            x168 * x248 * x70,
            x169 * x251 * x94,
            x171 * x253 * x94,
            x117 * x169 * x256,
            x117 * x171 * x251,
            x117 * x173 * x253,
            x134 * x169 * x260,
            x134 * x171 * x256,
            x134 * x173 * x251,
            x134 * x175 * x253,
            x148 * x169 * x265,
            x148 * x171 * x260,
            x148 * x173 * x256,
            x148 * x175 * x251,
            x148 * x179 * x253,
            x120 * x281 * x69,
            x120 * x284 * x52,
            x137 * x281 * x52,
            x120 * x288 * x87,
            x137 * x284 * x87,
            x124 * x281 * x87,
            x107 * x120 * x291,
            x107 * x137 * x288,
            x107 * x124 * x284,
            x107 * x142 * x281,
            x120 * x147 * x294,
            x137 * x147 * x291,
            x124 * x147 * x288,
            x142 * x147 * x284,
            x147 * x153 * x281,
            x169 * x268 * x69,
            x169 * x270 * x52,
            x171 * x268 * x52,
            x169 * x272 * x87,
            x171 * x270 * x87,
            x173 * x268 * x87,
            x107 * x169 * x274,
            x107 * x171 * x272,
            x107 * x173 * x270,
            x107 * x175 * x268,
            x147 * x169 * x278,
            x147 * x171 * x274,
            x147 * x173 * x272,
            x147 * x175 * x270,
            x147 * x179 * x268,
            x196 * x253 * x69,
            x196 * x251 * x52,
            x199 * x253 * x52,
            x196 * x256 * x87,
            x199 * x251 * x87,
            x203 * x253 * x87,
            x107 * x196 * x260,
            x107 * x199 * x256,
            x107 * x203 * x251,
            x107 * x206 * x253,
            x147 * x196 * x265,
            x147 * x199 * x260,
            x147 * x203 * x256,
            x147 * x206 * x251,
            x147 * x209 * x253,
            x120 * x295 * x39,
            x120 * x298 * x50,
            x137 * x295 * x50,
            x120 * x301 * x55,
            x137 * x298 * x55,
            x124 * x295 * x55,
            x105 * x120 * x305,
            x105 * x137 * x301,
            x105 * x124 * x298,
            x105 * x142 * x295,
            x24 * x307,
            x305 * x308 * x98,
            x103 * x124 * x301,
            x103 * x142 * x298,
            x103 * x153 * x295,
            x169 * x281 * x39,
            x169 * x284 * x50,
            x171 * x281 * x50,
            x169 * x288 * x55,
            x171 * x284 * x55,
            x173 * x281 * x55,
            x105 * x169 * x291,
            x105 * x171 * x288,
            x105 * x173 * x284,
            x105 * x175 * x281,
            x167 * x294 * x308,
            x103 * x171 * x291,
            x103 * x173 * x288,
            x103 * x175 * x284,
            x103 * x179 * x281,
            x196 * x268 * x39,
            x196 * x270 * x50,
            x199 * x268 * x50,
            x196 * x272 * x55,
            x199 * x270 * x55,
            x203 * x268 * x55,
            x105 * x196 * x274,
            x105 * x199 * x272,
            x105 * x203 * x270,
            x105 * x206 * x268,
            x103 * x196 * x278,
            x103 * x199 * x274,
            x103 * x203 * x272,
            x103 * x206 * x270,
            x103 * x209 * x268,
            x222 * x253 * x39,
            x222 * x251 * x50,
            x225 * x253 * x50,
            x222 * x256 * x55,
            x225 * x251 * x55,
            x228 * x253 * x55,
            x105 * x222 * x260,
            x105 * x225 * x256,
            x105 * x228 * x251,
            x105 * x232 * x253,
            x103 * x222 * x265,
            x103 * x225 * x260,
            x103 * x228 * x256,
            x103 * x232 * x251,
            x233 * x24 * x309,
            x120 * x310 * x37,
            x120 * x31 * x311,
            x137 * x31 * x310,
            x120 * x27 * x312,
            x137 * x27 * x311,
            x124 * x27 * x310,
            x3 * x313,
            x312 * x314 * x98,
            x124 * x311 * x9,
            x142 * x310 * x9,
            x238
            * (
                x0 * (x221 + 3.0 * x292 + 3.0 * x293 + 4.0 * x303 + 4.0 * x304)
                + x154 * x306
            ),
            x313 * x98,
            x124 * x312 * x8,
            x142 * x311 * x8,
            x153 * x310 * x8,
            x169 * x295 * x37,
            x169 * x298 * x31,
            x171 * x295 * x31,
            x169 * x27 * x301,
            x171 * x27 * x298,
            x173 * x27 * x295,
            x167 * x305 * x314,
            x171 * x301 * x9,
            x173 * x298 * x9,
            x175 * x295 * x9,
            x167 * x307,
            x171 * x305 * x8,
            x173 * x301 * x8,
            x175 * x298 * x8,
            x179 * x295 * x8,
            x196 * x281 * x37,
            x196 * x284 * x31,
            x199 * x281 * x31,
            x196 * x27 * x288,
            x199 * x27 * x284,
            x203 * x27 * x281,
            x196 * x291 * x9,
            x199 * x288 * x9,
            x203 * x284 * x9,
            x206 * x281 * x9,
            x196 * x294 * x8,
            x199 * x291 * x8,
            x203 * x288 * x8,
            x206 * x284 * x8,
            x209 * x281 * x8,
            x222 * x268 * x37,
            x222 * x270 * x31,
            x225 * x268 * x31,
            x222 * x27 * x272,
            x225 * x27 * x270,
            x228 * x268 * x27,
            x222 * x274 * x9,
            x225 * x272 * x9,
            x228 * x270 * x9,
            x232 * x268 * x9,
            x222 * x278 * x8,
            x225 * x274 * x8,
            x228 * x272 * x8,
            x232 * x270 * x8,
            x233 * x268 * x8,
            x243 * x253 * x37,
            x243 * x251 * x31,
            x244 * x253 * x31,
            x243 * x256 * x27,
            x244 * x251 * x27,
            x245 * x253 * x27,
            x243 * x260 * x9,
            x244 * x256 * x9,
            x245 * x251 * x9,
            x246 * x3 * x309,
            x243 * x265 * x8,
            x244 * x260 * x8,
            x245 * x256 * x8,
            x246 * x251 * x8,
            x247 * x248,
            x249 * x315,
            x252 * x316 * x78,
            x252 * x318 * x99,
            x102 * x257 * x319,
            x138 * x257 * x318,
            x257 * x322 * x99,
            x128 * x261 * x319,
            x102 * x261 * x318,
            x138 * x261 * x322,
            x261 * x325 * x99,
            x146 * x266 * x319,
            x128 * x266 * x318,
            x102 * x266 * x322,
            x138 * x266 * x325,
            x266 * x329 * x99,
            x154 * x316 * x70,
            x158 * x319 * x94,
            x156 * x318 * x94,
            x117 * x160 * x319,
            x117 * x158 * x318,
            x117 * x156 * x322,
            x134 * x162 * x319,
            x134 * x160 * x318,
            x134 * x158 * x322,
            x134 * x156 * x325,
            x148 * x166 * x319,
            x148 * x162 * x318,
            x148 * x160 * x322,
            x148 * x158 * x325,
            x148 * x156 * x329,
            x331 * x70 * x99,
            x138 * x331 * x94,
            x333 * x94 * x99,
            x102 * x117 * x331,
            x117 * x138 * x333,
            x117 * x335 * x99,
            x128 * x134 * x331,
            x102 * x134 * x333,
            x134 * x138 * x335,
            x134 * x337 * x99,
            x146 * x148 * x331,
            x128 * x148 * x333,
            x102 * x148 * x335,
            x138 * x148 * x337,
            x148 * x341 * x99,
            x181 * x319 * x69,
            x184 * x319 * x52,
            x181 * x318 * x52,
            x188 * x319 * x87,
            x184 * x318 * x87,
            x181 * x322 * x87,
            x107 * x191 * x319,
            x107 * x188 * x318,
            x107 * x184 * x322,
            x107 * x181 * x325,
            x147 * x194 * x319,
            x147 * x191 * x318,
            x147 * x188 * x322,
            x147 * x184 * x325,
            x147 * x181 * x329,
            x156 * x331 * x69,
            x158 * x331 * x52,
            x156 * x333 * x52,
            x160 * x331 * x87,
            x158 * x333 * x87,
            x156 * x335 * x87,
            x107 * x162 * x331,
            x107 * x160 * x333,
            x107 * x158 * x335,
            x107 * x156 * x337,
            x147 * x166 * x331,
            x147 * x162 * x333,
            x147 * x160 * x335,
            x147 * x158 * x337,
            x147 * x156 * x341,
            x344 * x69 * x99,
            x138 * x344 * x52,
            x347 * x52 * x99,
            x102 * x344 * x87,
            x138 * x347 * x87,
            x351 * x87 * x99,
            x107 * x128 * x344,
            x102 * x107 * x347,
            x107 * x138 * x351,
            x107 * x354 * x99,
            x146 * x147 * x344,
            x128 * x147 * x347,
            x102 * x147 * x351,
            x138 * x147 * x354,
            x147 * x357 * x99,
            x210 * x319 * x39,
            x213 * x319 * x50,
            x210 * x318 * x50,
            x216 * x319 * x55,
            x213 * x318 * x55,
            x210 * x322 * x55,
            x105 * x220 * x319,
            x105 * x216 * x318,
            x105 * x213 * x322,
            x105 * x210 * x325,
            x221 * x308 * x315,
            x103 * x220 * x318,
            x103 * x216 * x322,
            x103 * x213 * x325,
            x103 * x210 * x329,
            x181 * x331 * x39,
            x184 * x331 * x50,
            x181 * x333 * x50,
            x188 * x331 * x55,
            x184 * x333 * x55,
            x181 * x335 * x55,
            x105 * x191 * x331,
            x105 * x188 * x333,
            x105 * x184 * x335,
            x105 * x181 * x337,
            x103 * x194 * x331,
            x103 * x191 * x333,
            x103 * x188 * x335,
            x103 * x184 * x337,
            x103 * x181 * x341,
            x156 * x344 * x39,
            x158 * x344 * x50,
            x156 * x347 * x50,
            x160 * x344 * x55,
            x158 * x347 * x55,
            x156 * x351 * x55,
            x105 * x162 * x344,
            x105 * x160 * x347,
            x105 * x158 * x351,
            x105 * x156 * x354,
            x103 * x166 * x344,
            x103 * x162 * x347,
            x103 * x160 * x351,
            x103 * x158 * x354,
            x154 * x357 * x358,
            x359 * x39 * x99,
            x138 * x359 * x50,
            x362 * x50 * x99,
            x102 * x359 * x55,
            x138 * x362 * x55,
            x365 * x55 * x99,
            x105 * x128 * x359,
            x102 * x105 * x362,
            x105 * x138 * x365,
            x105 * x369 * x99,
            x103 * x146 * x359,
            x103 * x128 * x362,
            x102 * x103 * x365,
            x358 * x369 * x78,
            x24 * x371,
            x234 * x319 * x37,
            x235 * x31 * x319,
            x234 * x31 * x318,
            x236 * x27 * x319,
            x235 * x27 * x318,
            x234 * x27 * x322,
            x237 * x314 * x315,
            x236 * x318 * x9,
            x235 * x322 * x9,
            x234 * x325 * x9,
            x239 * x315,
            x237 * x318 * x8,
            x236 * x322 * x8,
            x235 * x325 * x8,
            x234 * x329 * x8,
            x210 * x331 * x37,
            x213 * x31 * x331,
            x210 * x31 * x333,
            x216 * x27 * x331,
            x213 * x27 * x333,
            x210 * x27 * x335,
            x220 * x331 * x9,
            x216 * x333 * x9,
            x213 * x335 * x9,
            x210 * x337 * x9,
            x221 * x331 * x8,
            x220 * x333 * x8,
            x216 * x335 * x8,
            x213 * x337 * x8,
            x210 * x341 * x8,
            x181 * x344 * x37,
            x184 * x31 * x344,
            x181 * x31 * x347,
            x188 * x27 * x344,
            x184 * x27 * x347,
            x181 * x27 * x351,
            x191 * x344 * x9,
            x188 * x347 * x9,
            x184 * x351 * x9,
            x181 * x354 * x9,
            x194 * x344 * x8,
            x191 * x347 * x8,
            x188 * x351 * x8,
            x184 * x354 * x8,
            x181 * x357 * x8,
            x156 * x359 * x37,
            x158 * x31 * x359,
            x156 * x31 * x362,
            x160 * x27 * x359,
            x158 * x27 * x362,
            x156 * x27 * x365,
            x162 * x359 * x9,
            x160 * x362 * x9,
            x158 * x365 * x9,
            x154 * x369 * x372,
            x166 * x359 * x8,
            x162 * x362 * x8,
            x160 * x365 * x8,
            x158 * x369 * x8,
            x154 * x371,
            x37 * x373 * x99,
            x138 * x31 * x373,
            x31 * x374 * x99,
            x102 * x27 * x373,
            x138 * x27 * x374,
            x27 * x375 * x99,
            x128 * x373 * x9,
            x102 * x374 * x9,
            x372 * x375 * x78,
            x3 * x376,
            x146 * x373 * x8,
            x128 * x374 * x8,
            x102 * x375 * x8,
            x376 * x78,
            x241
            * (
                x0 * (x233 + 3.0 * x355 + 3.0 * x356 + 4.0 * x367 + 4.0 * x368)
                + x167 * x370
            ),
        ]
    )
