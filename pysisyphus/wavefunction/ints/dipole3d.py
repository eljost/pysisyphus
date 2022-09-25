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
        numpy.pi ** (3 / 2)
        * x0 ** (3 / 2)
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
    x5 = numpy.sqrt(numpy.pi) * x4
    x6 = x5 / (2 * a + 2 * b)
    x7 = -x0 * (a * A[0] + b * B[0])
    x8 = -x7 - B[0]
    x9 = x3 * (-x7 - C[0])
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
            numpy.pi * x21 * (x10 * x6 + x18 * x24 * x5),
        ]
    )


def dipole3d_02(a, A, b, B, C):
    """Cartesian 3D (sd) dipole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = numpy.sqrt(x1)
    x3 = numpy.sqrt(numpy.pi) * x2
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
    x16 = numpy.pi * x1 * x15
    x17 = x14 * x16
    x18 = -x1 * (a * A[1] + b * B[1])
    x19 = -x18 - B[1]
    x20 = x13 * x17
    x21 = -x1 * (a * A[2] + b * B[2])
    x22 = -x21 - B[2]
    x23 = x14 * x3
    x24 = x0 * x23
    x25 = x16 * x7
    x26 = x25 * (x19 ** 2 * x23 + x24)
    x27 = numpy.pi ** (3 / 2)
    x28 = x1 * x14
    x29 = x15 * x2 * x22 * x27 * x28
    x30 = x15 * x3
    x31 = x0 * x30
    x32 = numpy.pi * x28
    x33 = x32 * x7
    x34 = x33 * (x22 ** 2 * x30 + x31)
    x35 = -x18 - C[1]
    x36 = x17 * (x11 * x5 ** 2 + x12)
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
    x10 = -x7 - C[0]
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
    x21 = x14 * x18
    x22 = -x1 * (a * A[2] + b * B[2])
    x23 = -x22 - B[2]
    x24 = x15 * x4
    x25 = x0 * x24
    x26 = x20 ** 2 * x24
    x27 = x25 + x26
    x28 = x16 * x4
    x29 = x18 * x23
    x30 = x0 * x28
    x31 = x23 ** 2 * x28
    x32 = x30 + x31
    x33 = x17 * x3
    x34 = x33 * (2 * x20 * x25 + x20 * x27)
    x35 = x23 * x33
    x36 = numpy.pi * x1 * x15 * x3
    x37 = x32 * x36
    x38 = x36 * (2 * x23 * x30 + x23 * x32)
    x39 = -x19 - C[1]
    x40 = x6 + x9
    x41 = x18 * (2 * x0 * x11 + x40 * x8)
    x42 = x20 * x24
    x43 = x39 * x42
    x44 = x25 + x43
    x45 = x0 * (x24 * x39 + x42) + x20 * x44
    x46 = x33 * x45
    x47 = -x22 - C[2]
    x48 = x23 * x28
    x49 = x47 * x48
    x50 = x30 + x49
    x51 = x0 * (x28 * x47 + x48) + x23 * x50
    x52 = x36 * x51

    # 30 item(s)
    return numpy.array(
        [
            x18 * (x0 * (2 * x12 + 3 * x6 + x9) + x14 * x8),
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
            x46 * x8,
            x35 * x44 * x8,
            x37 * x39 * x8,
            x33 * (x0 * (3 * x25 + x26 + 2 * x43) + x20 * x45),
            x23 * x46,
            x32 * x44 * x5,
            x38 * x39,
            x41 * x47,
            x18 * x20 * x40 * x47,
            x24 * x40 * x50,
            x27 * x33 * x47 * x8,
            x20 * x36 * x50 * x8,
            x52 * x8,
            x34 * x47,
            x27 * x5 * x50,
            x20 * x52,
            x36 * (x0 * (3 * x30 + x31 + 2 * x49) + x23 * x51),
        ]
    )


def dipole3d_04(a, A, b, B, C):
    """Cartesian 3D (sg) dipole moment integrals.
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
    x8 = x3 * x7
    x9 = -x2 - C[0]
    x10 = x7 * x9
    x11 = x0 * (x10 + x8)
    x12 = x0 * x7
    x13 = x8 * x9
    x14 = x12 + x13
    x15 = x14 * x3
    x16 = x3 ** 2 * x7
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
    x28 = x21 * x25
    x29 = -x1 * (a * A[2] + b * B[2])
    x30 = -x29 - B[2]
    x31 = x22 * x6
    x32 = x0 * x31
    x33 = x27 ** 2 * x31
    x34 = x32 + x33
    x35 = x23 * x6
    x36 = x25 * x30
    x37 = x0 * x35
    x38 = x30 ** 2 * x35
    x39 = x37 + x38
    x40 = 2 * x27 * x32 + x27 * x34
    x41 = x30 * x35
    x42 = x27 * x31
    x43 = 2 * x30 * x37 + x30 * x39
    x44 = 3 * x32
    x45 = x24 * x5
    x46 = x45 * (x0 * (3 * x33 + x44) + x27 * x40)
    x47 = x30 * x45
    x48 = numpy.pi * x1 * x22 * x5
    x49 = x43 * x48
    x50 = 3 * x37
    x51 = x48 * (x0 * (3 * x38 + x50) + x30 * x43)
    x52 = -x26 - C[1]
    x53 = x25 * (x0 * (3 * x16 + x19) + x18 * x3)
    x54 = x42 * x52
    x55 = x32 + x54
    x56 = x31 * x52
    x57 = x0 * (x42 + x56)
    x58 = x27 * x55
    x59 = x57 + x58
    x60 = x0 * (x33 + x44 + 2 * x54) + x27 * x59
    x61 = x45 * x60
    x62 = -x29 - C[2]
    x63 = x41 * x62
    x64 = x37 + x63
    x65 = x35 * x62
    x66 = x0 * (x41 + x65)
    x67 = x30 * x64
    x68 = x66 + x67
    x69 = x0 * (x38 + x50 + 2 * x63) + x30 * x68
    x70 = x48 * x69

    # 45 item(s)
    return numpy.array(
        [
            x25 * (x0 * (3 * x11 + 3 * x15 + x18) + x21 * x3),
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
            x45 * (x0 * (x40 + 3 * x57 + 3 * x58) + x27 * x60),
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
            x48 * (x0 * (x43 + 3 * x66 + 3 * x67) + x30 * x69),
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
    x5 = numpy.sqrt(numpy.pi) * x4
    x6 = x5 / (2 * a + 2 * b)
    x7 = -x0 * (a * A[0] + b * B[0])
    x8 = -x7 - A[0]
    x9 = x3 * (-x7 - C[0])
    x10 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x11 = numpy.pi * x0 * x10
    x12 = numpy.pi ** (3 / 2) * x4
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
            numpy.pi * x21 * (x10 * x6 + x18 * x24 * x5),
        ]
    )


def dipole3d_11(a, A, b, B, C):
    """Cartesian 3D (pp) dipole moment integrals.
    The origin is at C.

    Generated code; DO NOT modify by hand!"""

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = numpy.sqrt(x1)
    x3 = numpy.sqrt(numpy.pi) * x2
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
    x17 = numpy.pi * x1 * x16
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
    x33 = numpy.pi ** (3 / 2) * x16 * x2 * x32
    x34 = x33 * x6
    x35 = x34 * x9
    x36 = -x22 - A[2]
    x37 = x16 * x26
    x38 = x16 * x3
    x39 = x23 * x38
    x40 = numpy.pi * x32
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
    x10 = -x7 - C[0]
    x11 = x5 * x8
    x12 = x10 * x11
    x13 = -x7 - A[0]
    x14 = x10 * x5
    x15 = x0 * (x11 + x14)
    x16 = x12 + x6
    x17 = x15 + x16 * x8
    x18 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x19 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x20 = numpy.pi * x1 * x19
    x21 = x18 * x20
    x22 = -x1 * (a * A[1] + b * B[1])
    x23 = -x22 - B[1]
    x24 = x21 * (x13 * x16 + x15)
    x25 = -x1 * (a * A[2] + b * B[2])
    x26 = -x25 - B[2]
    x27 = x18 * x4
    x28 = x0 * x27
    x29 = x23 ** 2 * x27
    x30 = x28 + x29
    x31 = x13 * x14 + x6
    x32 = x19 * x4
    x33 = x21 * x26
    x34 = x0 * x32
    x35 = x26 ** 2 * x32
    x36 = x34 + x35
    x37 = -x22 - A[1]
    x38 = x17 * x21
    x39 = x23 * x27
    x40 = x28 + x37 * x39
    x41 = x20 * x3
    x42 = x41 * (2 * x23 * x28 + x30 * x37)
    x43 = x10 * x41
    x44 = numpy.pi * x1 * x18 * x3
    x45 = x10 * x44
    x46 = -x25 - A[2]
    x47 = x21 * x46
    x48 = x26 * x32
    x49 = x34 + x46 * x48
    x50 = x44 * (2 * x26 * x34 + x36 * x46)
    x51 = -x22 - C[1]
    x52 = x6 + x9
    x53 = x21 * (2 * x0 * x11 + x13 * x52)
    x54 = x39 * x51
    x55 = x28 + x54
    x56 = x11 * x13 + x6
    x57 = x27 * x51
    x58 = x0 * (x39 + x57)
    x59 = x23 * x55 + x58
    x60 = x41 * x59
    x61 = x26 * x41
    x62 = x44 * x51
    x63 = x28 + x37 * x57
    x64 = x41 * (x37 * x55 + x58)
    x65 = x41 * x8
    x66 = -x25 - C[2]
    x67 = x21 * x66
    x68 = x48 * x66
    x69 = x34 + x68
    x70 = x44 * x69
    x71 = x32 * x66
    x72 = x0 * (x48 + x71)
    x73 = x26 * x69 + x72
    x74 = x44 * x73
    x75 = x34 + x46 * x71
    x76 = x44 * (x46 * x69 + x72)

    # 54 item(s)
    return numpy.array(
        [
            x21 * (x0 * (2 * x12 + 3 * x6 + x9) + x13 * x17),
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
            x64 * x8,
            x61 * x63 * x8,
            x41 * (x0 * (3 * x28 + x29 + 2 * x54) + x37 * x59),
            x26 * x64,
            x36 * x5 * x63,
            x47 * x51 * x52,
            x46 * x55 * x65,
            x49 * x62 * x8,
            x46 * x60,
            x49 * x5 * x55,
            x50 * x51,
            x53 * x66,
            x23 * x56 * x67,
            x27 * x56 * x69,
            x13 * x30 * x41 * x66,
            x13 * x23 * x70,
            x13 * x74,
            x37 * x52 * x67,
            x40 * x65 * x66,
            x37 * x70 * x8,
            x42 * x66,
            x40 * x5 * x69,
            x37 * x74,
            x27 * x52 * x75,
            x23 * x44 * x75 * x8,
            x76 * x8,
            x30 * x5 * x75,
            x23 * x76,
            x44 * (x0 * (3 * x34 + x35 + 2 * x68) + x46 * x73),
        ]
    )


def dipole3d_13(a, A, b, B, C):
    """Cartesian 3D (pf) dipole moment integrals.
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
    x8 = x3 * x7
    x9 = -x2 - C[0]
    x10 = x7 * x9
    x11 = x0 * (x10 + x8)
    x12 = x0 * x7
    x13 = x8 * x9
    x14 = x12 + x13
    x15 = x14 * x3
    x16 = 2 * x12 * x3
    x17 = x3 ** 2 * x7
    x18 = x12 + x17
    x19 = x16 + x18 * x3
    x20 = -x2 - A[0]
    x21 = 3 * x12
    x22 = x0 * (2 * x13 + x17 + x21)
    x23 = x11 + x15
    x24 = x22 + x23 * x3
    x25 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x26 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x27 = numpy.pi * x1 * x26
    x28 = x25 * x27
    x29 = -x1 * (a * A[1] + b * B[1])
    x30 = -x29 - B[1]
    x31 = x28 * (x20 * x23 + x22)
    x32 = -x1 * (a * A[2] + b * B[2])
    x33 = -x32 - B[2]
    x34 = x11 + x14 * x20
    x35 = x25 * x6
    x36 = x0 * x35
    x37 = x30 ** 2 * x35
    x38 = x36 + x37
    x39 = x26 * x6
    x40 = x38 * x39
    x41 = x28 * x33
    x42 = x0 * x39
    x43 = x33 ** 2 * x39
    x44 = x42 + x43
    x45 = x35 * x44
    x46 = 2 * x30 * x36
    x47 = x30 * x38 + x46
    x48 = x10 * x20 + x12
    x49 = x33 * x39
    x50 = x30 * x35
    x51 = 2 * x33 * x42
    x52 = x33 * x44 + x51
    x53 = -x29 - A[1]
    x54 = x24 * x28
    x55 = x36 + x50 * x53
    x56 = x38 * x53 + x46
    x57 = 3 * x36
    x58 = x27 * x5
    x59 = x58 * (x0 * (3 * x37 + x57) + x47 * x53)
    x60 = x58 * x9
    x61 = numpy.pi * x1 * x25 * x5
    x62 = x61 * x9
    x63 = -x32 - A[2]
    x64 = x28 * x63
    x65 = x42 + x49 * x63
    x66 = x44 * x63 + x51
    x67 = 3 * x42
    x68 = x61 * (x0 * (3 * x43 + x67) + x52 * x63)
    x69 = -x29 - C[1]
    x70 = x28 * (x0 * (3 * x17 + x21) + x19 * x20)
    x71 = x16 + x18 * x20
    x72 = x50 * x69
    x73 = x36 + x72
    x74 = x39 * x73
    x75 = x35 * x69
    x76 = x0 * (x50 + x75)
    x77 = x30 * x73
    x78 = x76 + x77
    x79 = x12 + x20 * x8
    x80 = x0 * (x37 + x57 + 2 * x72)
    x81 = x30 * x78 + x80
    x82 = x58 * x81
    x83 = x33 * x58
    x84 = x44 * x7
    x85 = x61 * x69
    x86 = x36 + x53 * x75
    x87 = x53 * x73 + x76
    x88 = x58 * (x53 * x78 + x80)
    x89 = x3 * x58
    x90 = -x32 - C[2]
    x91 = x28 * x90
    x92 = x49 * x90
    x93 = x42 + x92
    x94 = x35 * x93
    x95 = x39 * x90
    x96 = x0 * (x49 + x95)
    x97 = x33 * x93
    x98 = x96 + x97
    x99 = x7 * x93
    x100 = x61 * x98
    x101 = x0 * (x43 + x67 + 2 * x92)
    x102 = x101 + x33 * x98
    x103 = x102 * x61
    x104 = x42 + x63 * x95
    x105 = x63 * x93 + x96
    x106 = x61 * (x101 + x63 * x98)

    # 90 item(s)
    return numpy.array(
        [
            x28 * (x0 * (3 * x11 + 3 * x15 + x19) + x20 * x24),
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
            x18 * x39 * x87,
            x18 * x49 * x86,
            x3 * x88,
            x3 * x83 * x87,
            x44 * x8 * x86,
            x58 * (x0 * (x47 + 3 * x76 + 3 * x77) + x53 * x81),
            x33 * x88,
            x84 * x87,
            x52 * x7 * x86,
            x19 * x64 * x69,
            x18 * x63 * x74,
            x18 * x65 * x75,
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
            x18 * x55 * x95,
            x18 * x53 * x94,
            x56 * x89 * x90,
            x55 * x8 * x93,
            x100 * x3 * x53,
            x59 * x90,
            x56 * x99,
            x55 * x7 * x98,
            x103 * x53,
            x104 * x19 * x35,
            x104 * x18 * x50,
            x105 * x18 * x35,
            x104 * x38 * x8,
            x105 * x3 * x30 * x61,
            x106 * x3,
            x104 * x47 * x7,
            x105 * x38 * x7,
            x106 * x30,
            x61 * (x0 * (x52 + 3 * x96 + 3 * x97) + x102 * x63),
        ]
    )


def dipole3d_14(a, A, b, B, C):
    """Cartesian 3D (pg) dipole moment integrals.
    The origin is at C.

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
    x10 = x5 * x9 ** 2
    x11 = -x8 - C[0]
    x12 = x5 * x9
    x13 = x11 * x12
    x14 = x0 * (x10 + 2 * x13 + x7)
    x15 = x11 * x5
    x16 = x0 * (x12 + x15)
    x17 = x13 + x6
    x18 = x17 * x9
    x19 = x16 + x18
    x20 = x19 * x9
    x21 = x0 * (3 * x10 + x7)
    x22 = x6 * x9
    x23 = 2 * x22
    x24 = x10 + x6
    x25 = x24 * x9
    x26 = x23 + x25
    x27 = x21 + x26 * x9
    x28 = -x8 - A[0]
    x29 = x0 * (3 * x16 + 3 * x18 + x26)
    x30 = x14 + x20
    x31 = x29 + x30 * x9
    x32 = numpy.exp(-x2 * (A[1] - B[1]) ** 2)
    x33 = numpy.exp(-x2 * (A[2] - B[2]) ** 2)
    x34 = numpy.pi * x1 * x33
    x35 = x32 * x34
    x36 = -x1 * (a * A[1] + b * B[1])
    x37 = -x36 - B[1]
    x38 = x35 * (x28 * x30 + x29)
    x39 = -x1 * (a * A[2] + b * B[2])
    x40 = -x39 - B[2]
    x41 = x14 + x19 * x28
    x42 = x32 * x4
    x43 = x0 * x42
    x44 = x37 ** 2 * x42
    x45 = x43 + x44
    x46 = x33 * x4
    x47 = x45 * x46
    x48 = x35 * x40
    x49 = x0 * x46
    x50 = x40 ** 2 * x46
    x51 = x49 + x50
    x52 = x42 * x51
    x53 = x16 + x17 * x28
    x54 = x37 * x43
    x55 = 2 * x54
    x56 = x37 * x45
    x57 = x55 + x56
    x58 = x46 * x57
    x59 = x40 * x46
    x60 = x37 * x42
    x61 = x40 * x49
    x62 = 2 * x61
    x63 = x40 * x51
    x64 = x62 + x63
    x65 = x42 * x64
    x66 = 3 * x43
    x67 = x0 * (3 * x44 + x66)
    x68 = x37 * x57 + x67
    x69 = x15 * x28 + x6
    x70 = 3 * x49
    x71 = x0 * (3 * x50 + x70)
    x72 = x40 * x64 + x71
    x73 = -x36 - A[1]
    x74 = x31 * x35
    x75 = x43 + x60 * x73
    x76 = x45 * x73 + x55
    x77 = x57 * x73 + x67
    x78 = x3 * x34
    x79 = x78 * (x0 * (8 * x54 + 4 * x56) + x68 * x73)
    x80 = x11 * x78
    x81 = numpy.pi * x1 * x3 * x32
    x82 = x11 * x81
    x83 = -x39 - A[2]
    x84 = x35 * x83
    x85 = x49 + x59 * x83
    x86 = x51 * x83 + x62
    x87 = x64 * x83 + x71
    x88 = x81 * (x0 * (8 * x61 + 4 * x63) + x72 * x83)
    x89 = -x36 - C[1]
    x90 = x35 * (x0 * (8 * x22 + 4 * x25) + x27 * x28)
    x91 = x21 + x26 * x28
    x92 = x60 * x89
    x93 = x43 + x92
    x94 = x46 * x93
    x95 = x23 + x24 * x28
    x96 = x42 * x89
    x97 = x0 * (x60 + x96)
    x98 = x37 * x93
    x99 = x97 + x98
    x100 = x46 * x99
    x101 = x0 * (x44 + x66 + 2 * x92)
    x102 = x37 * x99
    x103 = x101 + x102
    x104 = x12 * x28 + x6
    x105 = x0 * (x57 + 3 * x97 + 3 * x98)
    x106 = x103 * x37 + x105
    x107 = x106 * x78
    x108 = x40 * x78
    x109 = x5 * x51
    x110 = x5 * x64
    x111 = x81 * x89
    x112 = x43 + x73 * x96
    x113 = x73 * x93 + x97
    x114 = x101 + x73 * x99
    x115 = x78 * (x103 * x73 + x105)
    x116 = x78 * x9
    x117 = -x39 - C[2]
    x118 = x117 * x35
    x119 = x117 * x59
    x120 = x119 + x49
    x121 = x120 * x42
    x122 = x117 * x46
    x123 = x0 * (x122 + x59)
    x124 = x120 * x40
    x125 = x123 + x124
    x126 = x125 * x42
    x127 = x0 * (2 * x119 + x50 + x70)
    x128 = x125 * x40
    x129 = x127 + x128
    x130 = x120 * x5
    x131 = x125 * x5
    x132 = x129 * x81
    x133 = x0 * (3 * x123 + 3 * x124 + x64)
    x134 = x129 * x40 + x133
    x135 = x134 * x81
    x136 = x122 * x83 + x49
    x137 = x120 * x83 + x123
    x138 = x125 * x83 + x127
    x139 = x81 * (x129 * x83 + x133)

    # 135 item(s)
    return numpy.array(
        [
            x35 * (x0 * (4 * x14 + 4 * x20 + x27) + x28 * x31),
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
            x114 * x24 * x46,
            x113 * x24 * x59,
            x112 * x24 * x51,
            x115 * x9,
            x108 * x114 * x9,
            x113 * x12 * x51,
            x112 * x12 * x64,
            x78 * (x0 * (4 * x101 + 4 * x102 + x68) + x106 * x73),
            x115 * x40,
            x109 * x114,
            x110 * x113,
            x112 * x5 * x72,
            x27 * x84 * x89,
            x26 * x83 * x94,
            x26 * x85 * x96,
            x100 * x24 * x83,
            x24 * x85 * x93,
            x24 * x86 * x96,
            x103 * x116 * x83,
            x12 * x85 * x99,
            x12 * x86 * x93,
            x111 * x87 * x9,
            x107 * x83,
            x103 * x5 * x85,
            x5 * x86 * x99,
            x5 * x87 * x93,
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
            x122 * x24 * x76,
            x120 * x24 * x75,
            x126 * x24 * x73,
            x116 * x117 * x77,
            x12 * x120 * x76,
            x12 * x125 * x75,
            x132 * x73 * x9,
            x117 * x79,
            x130 * x77,
            x131 * x76,
            x129 * x5 * x75,
            x135 * x73,
            x136 * x27 * x42,
            x136 * x26 * x60,
            x137 * x26 * x42,
            x136 * x24 * x45,
            x137 * x24 * x60,
            x138 * x24 * x42,
            x12 * x136 * x57,
            x12 * x137 * x45,
            x138 * x37 * x81 * x9,
            x139 * x9,
            x136 * x5 * x68,
            x137 * x5 * x57,
            x138 * x45 * x5,
            x139 * x37,
            x81 * (x0 * (4 * x127 + 4 * x128 + x72) + x134 * x83),
        ]
    )


def dipole3d_20(a, A, b, B, C):
    """Cartesian 3D (ds) dipole moment integrals.
    The origin is at C.

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
    x10 = -x4 - C[0]
    x11 = x3 * x7
    x12 = x0 * x11
    x13 = x10 * x9 + x12
    x14 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x15 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x16 = numpy.pi * x1 * x15
    x17 = x14 * x16
    x18 = -x1 * (a * A[1] + b * B[1])
    x19 = -x18 - A[1]
    x20 = x13 * x17
    x21 = -x1 * (a * A[2] + b * B[2])
    x22 = -x21 - A[2]
    x23 = x14 * x3
    x24 = x0 * x23
    x25 = x16 * x7
    x26 = x25 * (x19 ** 2 * x23 + x24)
    x27 = numpy.pi ** (3 / 2)
    x28 = x1 * x14
    x29 = x15 * x2 * x22 * x27 * x28
    x30 = x15 * x3
    x31 = x0 * x30
    x32 = numpy.pi * x28
    x33 = x32 * x7
    x34 = x33 * (x22 ** 2 * x30 + x31)
    x35 = -x18 - C[1]
    x36 = x17 * (x11 * x5 ** 2 + x12)
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

    x0 = (2 * a + 2 * b) ** (-1.0)
    x1 = (a + b) ** (-1.0)
    x2 = a * b * x1
    x3 = numpy.exp(-x2 * (A[0] - B[0]) ** 2)
    x4 = numpy.sqrt(numpy.pi) * numpy.sqrt(x1)
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
    x20 = numpy.pi * x1 * x19
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
    x44 = x29 ** 2 * x33 + x32
    x45 = x23 * x33
    x46 = x20 * x3
    x47 = x46 * (x0 * (x34 + x45) + x29 * x36)
    x48 = x28 * x46
    x49 = x21 * x39
    x50 = x39 * x46
    x51 = numpy.pi * x1 * x18
    x52 = x3 * x51
    x53 = x43 * x52
    x54 = x37 * x39 ** 2 + x40
    x55 = x23 * x52
    x56 = x28 * x37
    x57 = x52 * (x0 * (x41 + x56) + x39 * x43)
    x58 = -x22 - C[1]
    x59 = x11 + x6
    x60 = x21 * (x0 * (x10 + x24) + x59 * x8)
    x61 = x45 * x58
    x62 = x32 + x61
    x63 = x5 * x8 ** 2 + x6
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
            x21 * (x0 * (x11 + x14 + x15 + 3 * x6) + x17 * x8),
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
            x46 * (x0 * (3 * x32 + x35 + x61 + x65) + x29 * x68),
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
            x52 * (x0 * (3 * x40 + x42 + x74 + x76) + x39 * x79),
        ]
    )


def dipole3d_22(a, A, b, B, C):
    """Cartesian 3D (dd) dipole moment integrals.
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
    x12 = x5 ** 2 * x9
    x13 = x3 * x9
    x14 = 3 * x13
    x15 = x12 + x14
    x16 = x4 * x9
    x17 = x3 * (x10 + x16)
    x18 = x11 + x13
    x19 = x18 * x5
    x20 = x17 + x19
    x21 = x2 * x20 + x3 * (2 * x11 + x15)
    x22 = x18 * x2
    x23 = x12 + x13
    x24 = 2 * x13 * x5 + x2 * x23
    x25 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x26 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x27 = numpy.pi * x0 * x26
    x28 = x25 * x27
    x29 = -x0 * (a * A[1] + b * B[1])
    x30 = -x29 - B[1]
    x31 = x10 * x2
    x32 = x16 * x2
    x33 = x17 + x22
    x34 = x28 * (x2 * x33 + x3 * (x11 + x14 + x31 + x32))
    x35 = -x0 * (a * A[2] + b * B[2])
    x36 = -x35 - B[2]
    x37 = x25 * x8
    x38 = x3 * x37
    x39 = x30 ** 2 * x37
    x40 = x38 + x39
    x41 = x2 * x9
    x42 = x13 + x32
    x43 = x2 * x42 + x3 * (x16 + x41)
    x44 = x26 * x8
    x45 = x28 * x36
    x46 = x3 * x44
    x47 = x36 ** 2 * x44
    x48 = x46 + x47
    x49 = -x29 - A[1]
    x50 = x21 * x28
    x51 = x37 * x49
    x52 = x30 * x51
    x53 = x38 + x52
    x54 = 2 * x30 * x38 + x40 * x49
    x55 = x36 * x44
    x56 = -x35 - A[2]
    x57 = x28 * x56
    x58 = x44 * x56
    x59 = x36 * x58
    x60 = x46 + x59
    x61 = x30 * x37
    x62 = 2 * x36 * x46 + x48 * x56
    x63 = x37 * x49 ** 2 + x38
    x64 = x3 * (x51 + x61) + x49 * x53
    x65 = 3 * x38
    x66 = x39 + x65
    x67 = x27 * x7
    x68 = x67 * (x3 * (2 * x52 + x66) + x49 * x54)
    x69 = x4 * x67
    x70 = numpy.pi * x0 * x25 * x7
    x71 = x4 * x70
    x72 = x44 * x56 ** 2 + x46
    x73 = x3 * (x55 + x58) + x56 * x60
    x74 = 3 * x46
    x75 = x47 + x74
    x76 = x70 * (x3 * (2 * x59 + x75) + x56 * x62)
    x77 = -x29 - C[1]
    x78 = x28 * (x2 * x24 + x3 * (x15 + 2 * x31))
    x79 = x61 * x77
    x80 = x38 + x79
    x81 = x13 + x31
    x82 = x2 * x81 + x3 * (x10 + x41)
    x83 = x37 * x77
    x84 = x3 * (x61 + x83)
    x85 = x30 * x80
    x86 = x84 + x85
    x87 = x13 + x2 ** 2 * x9
    x88 = x51 * x77
    x89 = x38 + x88
    x90 = x49 * x80
    x91 = x84 + x90
    x92 = x3 * (x66 + 2 * x79) + x49 * x86
    x93 = x67 * x92
    x94 = x2 * x67
    x95 = x70 * x77
    x96 = x3 * (x51 + x83) + x49 * x89
    x97 = x67 * (x3 * (x52 + x65 + x79 + x88) + x49 * x91)
    x98 = x5 * x67
    x99 = -x35 - C[2]
    x100 = x28 * x99
    x101 = x55 * x99
    x102 = x101 + x46
    x103 = x44 * x99
    x104 = x3 * (x103 + x55)
    x105 = x102 * x36
    x106 = x104 + x105
    x107 = x2 * x70
    x108 = x58 * x99
    x109 = x108 + x46
    x110 = x102 * x56
    x111 = x104 + x110
    x112 = x106 * x56 + x3 * (2 * x101 + x75)
    x113 = x112 * x70
    x114 = x5 * x70
    x115 = x109 * x56 + x3 * (x103 + x58)
    x116 = x70 * (x111 * x56 + x3 * (x101 + x108 + x59 + x74))

    # 108 item(s)
    return numpy.array(
        [
            x28 * (x2 * x21 + x3 * (3 * x17 + x19 + 2 * x22 + x24)),
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
            x20 * x44 * x63,
            x18 * x44 * x64,
            x18 * x55 * x63,
            x4 * x68,
            x36 * x64 * x69,
            x16 * x48 * x63,
            x20 * x49 * x57,
            x18 * x53 * x58,
            x18 * x51 * x60,
            x54 * x56 * x69,
            x16 * x53 * x60,
            x49 * x62 * x71,
            x20 * x37 * x72,
            x18 * x61 * x72,
            x18 * x37 * x73,
            x16 * x40 * x72,
            x30 * x71 * x73,
            x4 * x76,
            x77 * x78,
            x44 * x80 * x82,
            x45 * x77 * x82,
            x44 * x86 * x87,
            x55 * x80 * x87,
            x48 * x83 * x87,
            x24 * x44 * x89,
            x44 * x81 * x91,
            x55 * x81 * x89,
            x2 * x93,
            x36 * x91 * x94,
            x41 * x48 * x89,
            x24 * x57 * x77,
            x58 * x80 * x81,
            x60 * x81 * x83,
            x56 * x86 * x94,
            x41 * x60 * x80,
            x2 * x62 * x95,
            x23 * x44 * x96,
            x5 * x97,
            x36 * x96 * x98,
            x67 * (x3 * (x54 + 3 * x84 + x85 + 2 * x90) + x49 * x92),
            x36 * x97,
            x48 * x9 * x96,
            x23 * x58 * x89,
            x56 * x91 * x98,
            x10 * x60 * x89,
            x56 * x93,
            x60 * x9 * x91,
            x62 * x89 * x9,
            x23 * x72 * x83,
            x10 * x72 * x80,
            x5 * x73 * x95,
            x72 * x86 * x9,
            x73 * x80 * x9,
            x76 * x77,
            x78 * x99,
            x100 * x30 * x82,
            x102 * x37 * x82,
            x103 * x40 * x87,
            x102 * x61 * x87,
            x106 * x37 * x87,
            x100 * x24 * x49,
            x103 * x53 * x81,
            x102 * x51 * x81,
            x54 * x94 * x99,
            x102 * x41 * x53,
            x106 * x107 * x49,
            x109 * x24 * x37,
            x109 * x61 * x81,
            x111 * x37 * x81,
            x109 * x40 * x41,
            x107 * x111 * x30,
            x113 * x2,
            x103 * x23 * x63,
            x64 * x98 * x99,
            x10 * x102 * x63,
            x68 * x99,
            x102 * x64 * x9,
            x106 * x63 * x9,
            x109 * x23 * x51,
            x10 * x109 * x53,
            x111 * x114 * x49,
            x109 * x54 * x9,
            x111 * x53 * x9,
            x113 * x49,
            x115 * x23 * x37,
            x114 * x115 * x30,
            x116 * x5,
            x115 * x40 * x9,
            x116 * x30,
            x70 * (x112 * x56 + x3 * (3 * x104 + x105 + 2 * x110 + x62)),
        ]
    )


def dipole3d_23(a, A, b, B, C):
    """Cartesian 3D (df) dipole moment integrals.
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
    x13 = 3 * x12
    x14 = x3 * x8
    x15 = x10 * x9
    x16 = x14 + x15
    x17 = x16 * x4
    x18 = x14 * x4
    x19 = 2 * x18
    x20 = x4 ** 2 * x8
    x21 = x14 + x20
    x22 = x21 * x4
    x23 = x19 + x22
    x24 = 3 * x14
    x25 = x20 + x24
    x26 = x3 * (2 * x15 + x25)
    x27 = x12 + x17
    x28 = x27 * x4
    x29 = x26 + x28
    x30 = x2 * x29 + x3 * (x13 + 3 * x17 + x23)
    x31 = x2 * x27
    x32 = x2 * x23 + x3 * (3 * x20 + x24)
    x33 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x34 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x35 = numpy.pi * x0 * x34
    x36 = x33 * x35
    x37 = -x0 * (a * A[1] + b * B[1])
    x38 = -x37 - B[1]
    x39 = x26 + x31
    x40 = x16 * x2
    x41 = x2 * x21
    x42 = x19 + x41
    x43 = x36 * (x2 * x39 + x3 * (x13 + x17 + 2 * x40 + x42))
    x44 = -x0 * (a * A[2] + b * B[2])
    x45 = -x44 - B[2]
    x46 = x33 * x7
    x47 = x3 * x46
    x48 = x38 ** 2 * x46
    x49 = x47 + x48
    x50 = x2 * x9
    x51 = x11 * x2
    x52 = x12 + x40
    x53 = x2 * x52 + x3 * (x15 + x24 + x50 + x51)
    x54 = x34 * x7
    x55 = x36 * x45
    x56 = x3 * x54
    x57 = x45 ** 2 * x54
    x58 = x56 + x57
    x59 = x38 * x47
    x60 = 2 * x59
    x61 = x38 * x49
    x62 = x60 + x61
    x63 = x2 * x8
    x64 = x14 + x51
    x65 = x2 * x64 + x3 * (x11 + x63)
    x66 = x45 * x54
    x67 = x38 * x46
    x68 = x45 * x56
    x69 = 2 * x68
    x70 = x45 * x58
    x71 = x69 + x70
    x72 = -x37 - A[1]
    x73 = x30 * x36
    x74 = x46 * x72
    x75 = x38 * x74
    x76 = x47 + x75
    x77 = x49 * x72
    x78 = x60 + x77
    x79 = 3 * x47
    x80 = x3 * (3 * x48 + x79) + x62 * x72
    x81 = -x44 - A[2]
    x82 = x36 * x81
    x83 = x54 * x81
    x84 = x45 * x83
    x85 = x56 + x84
    x86 = x58 * x81
    x87 = x69 + x86
    x88 = 3 * x56
    x89 = x3 * (3 * x57 + x88) + x71 * x81
    x90 = x46 * x72 ** 2 + x47
    x91 = x3 * (x67 + x74) + x72 * x76
    x92 = x48 + x79
    x93 = x3 * (2 * x75 + x92) + x72 * x78
    x94 = x35 * x6
    x95 = x94 * (x3 * (8 * x59 + x61 + 3 * x77) + x72 * x80)
    x96 = x10 * x94
    x97 = numpy.pi * x0 * x33 * x6
    x98 = x10 * x97
    x99 = x54 * x81 ** 2 + x56
    x100 = x3 * (x66 + x83) + x81 * x85
    x101 = x57 + x88
    x102 = x3 * (x101 + 2 * x84) + x81 * x87
    x103 = x97 * (x3 * (8 * x68 + x70 + 3 * x86) + x81 * x89)
    x104 = -x37 - C[1]
    x105 = x36 * (x2 * x32 + x3 * (8 * x18 + x22 + 3 * x41))
    x106 = x104 * x67
    x107 = x106 + x47
    x108 = x2 * x42 + x3 * (x25 + 2 * x50)
    x109 = x104 * x46
    x110 = x3 * (x109 + x67)
    x111 = x107 * x38
    x112 = x110 + x111
    x113 = x14 + x50
    x114 = x113 * x2 + x3 * (x63 + x9)
    x115 = x3 * (2 * x106 + x92)
    x116 = x112 * x38
    x117 = x115 + x116
    x118 = x14 + x2 ** 2 * x8
    x119 = x104 * x74
    x120 = x119 + x47
    x121 = x107 * x72
    x122 = x110 + x121
    x123 = x112 * x72
    x124 = x115 + x123
    x125 = 3 * x110
    x126 = x117 * x72 + x3 * (3 * x111 + x125 + x62)
    x127 = x126 * x94
    x128 = x2 * x94
    x129 = x104 * x97
    x130 = x120 * x72 + x3 * (x109 + x74)
    x131 = x122 * x72 + x3 * (x106 + x119 + x75 + x79)
    x132 = x94 * (x124 * x72 + x3 * (x111 + 2 * x121 + x125 + x78))
    x133 = x4 * x94
    x134 = -x44 - C[2]
    x135 = x134 * x36
    x136 = x134 * x66
    x137 = x136 + x56
    x138 = x134 * x54
    x139 = x3 * (x138 + x66)
    x140 = x137 * x45
    x141 = x139 + x140
    x142 = x3 * (x101 + 2 * x136)
    x143 = x141 * x45
    x144 = x142 + x143
    x145 = x2 * x97
    x146 = x134 * x83
    x147 = x146 + x56
    x148 = x137 * x81
    x149 = x139 + x148
    x150 = x141 * x81
    x151 = x142 + x150
    x152 = 3 * x139
    x153 = x144 * x81 + x3 * (3 * x140 + x152 + x71)
    x154 = x153 * x97
    x155 = x4 * x97
    x156 = x147 * x81 + x3 * (x138 + x83)
    x157 = x149 * x81 + x3 * (x136 + x146 + x84 + x88)
    x158 = x97 * (x151 * x81 + x3 * (x140 + 2 * x148 + x152 + x87))

    # 180 item(s)
    return numpy.array(
        [
            x36 * (x2 * x30 + x3 * (4 * x26 + x28 + 3 * x31 + x32)),
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
            x39 * x54 * x76,
            x39 * x55 * x72,
            x52 * x54 * x78,
            x52 * x66 * x76,
            x52 * x58 * x74,
            x54 * x64 * x80,
            x64 * x66 * x78,
            x58 * x64 * x76,
            x64 * x71 * x74,
            x73 * x81,
            x38 * x39 * x82,
            x39 * x46 * x85,
            x49 * x52 * x83,
            x52 * x67 * x85,
            x46 * x52 * x87,
            x62 * x64 * x83,
            x49 * x64 * x85,
            x64 * x67 * x87,
            x46 * x64 * x89,
            x29 * x54 * x90,
            x27 * x54 * x91,
            x27 * x66 * x90,
            x16 * x54 * x93,
            x16 * x66 * x91,
            x16 * x58 * x90,
            x10 * x95,
            x45 * x93 * x96,
            x11 * x58 * x91,
            x11 * x71 * x90,
            x29 * x72 * x82,
            x27 * x76 * x83,
            x27 * x74 * x85,
            x16 * x78 * x83,
            x16 * x76 * x85,
            x16 * x74 * x87,
            x80 * x81 * x96,
            x11 * x78 * x85,
            x11 * x76 * x87,
            x72 * x89 * x98,
            x29 * x46 * x99,
            x27 * x67 * x99,
            x100 * x27 * x46,
            x16 * x49 * x99,
            x100 * x16 * x67,
            x102 * x16 * x46,
            x11 * x62 * x99,
            x100 * x11 * x49,
            x102 * x38 * x98,
            x10 * x103,
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
            x120 * x32 * x54,
            x122 * x42 * x54,
            x120 * x42 * x66,
            x113 * x124 * x54,
            x113 * x122 * x66,
            x113 * x120 * x58,
            x127 * x2,
            x124 * x128 * x45,
            x122 * x58 * x63,
            x120 * x63 * x71,
            x104 * x32 * x82,
            x107 * x42 * x83,
            x109 * x42 * x85,
            x112 * x113 * x83,
            x107 * x113 * x85,
            x109 * x113 * x87,
            x117 * x128 * x81,
            x112 * x63 * x85,
            x107 * x63 * x87,
            x129 * x2 * x89,
            x130 * x23 * x54,
            x131 * x21 * x54,
            x130 * x21 * x66,
            x132 * x4,
            x131 * x133 * x45,
            x130 * x58 * x9,
            x94 * (x126 * x72 + x3 * (4 * x115 + x116 + 3 * x123 + x80)),
            x132 * x45,
            x131 * x58 * x8,
            x130 * x71 * x8,
            x120 * x23 * x83,
            x122 * x21 * x83,
            x120 * x21 * x85,
            x124 * x133 * x81,
            x122 * x85 * x9,
            x120 * x87 * x9,
            x127 * x81,
            x124 * x8 * x85,
            x122 * x8 * x87,
            x120 * x8 * x89,
            x109 * x23 * x99,
            x107 * x21 * x99,
            x100 * x109 * x21,
            x112 * x9 * x99,
            x100 * x107 * x9,
            x102 * x129 * x4,
            x117 * x8 * x99,
            x100 * x112 * x8,
            x102 * x107 * x8,
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
            x135 * x32 * x72,
            x138 * x42 * x76,
            x137 * x42 * x74,
            x113 * x138 * x78,
            x113 * x137 * x76,
            x113 * x141 * x74,
            x128 * x134 * x80,
            x137 * x63 * x78,
            x141 * x63 * x76,
            x144 * x145 * x72,
            x147 * x32 * x46,
            x147 * x42 * x67,
            x149 * x42 * x46,
            x113 * x147 * x49,
            x113 * x149 * x67,
            x113 * x151 * x46,
            x147 * x62 * x63,
            x149 * x49 * x63,
            x145 * x151 * x38,
            x154 * x2,
            x138 * x23 * x90,
            x138 * x21 * x91,
            x137 * x21 * x90,
            x133 * x134 * x93,
            x137 * x9 * x91,
            x141 * x9 * x90,
            x134 * x95,
            x137 * x8 * x93,
            x141 * x8 * x91,
            x144 * x8 * x90,
            x147 * x23 * x74,
            x147 * x21 * x76,
            x149 * x21 * x74,
            x147 * x78 * x9,
            x149 * x76 * x9,
            x151 * x155 * x72,
            x147 * x8 * x80,
            x149 * x78 * x8,
            x151 * x76 * x8,
            x154 * x72,
            x156 * x23 * x46,
            x156 * x21 * x67,
            x157 * x21 * x46,
            x156 * x49 * x9,
            x155 * x157 * x38,
            x158 * x4,
            x156 * x62 * x8,
            x157 * x49 * x8,
            x158 * x38,
            x97 * (x153 * x81 + x3 * (4 * x142 + x143 + 3 * x150 + x89)),
        ]
    )


def dipole3d_24(a, A, b, B, C):
    """Cartesian 3D (dg) dipole moment integrals.
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
    x12 = x5 ** 2 * x9
    x13 = x3 * x9
    x14 = 3 * x13
    x15 = x12 + x14
    x16 = x3 * (2 * x11 + x15)
    x17 = 4 * x16
    x18 = x4 * x9
    x19 = x3 * (x10 + x18)
    x20 = x11 + x13
    x21 = x20 * x5
    x22 = x19 + x21
    x23 = x22 * x5
    x24 = x3 * (3 * x12 + x14)
    x25 = x13 * x5
    x26 = 2 * x25
    x27 = x12 + x13
    x28 = x27 * x5
    x29 = x26 + x28
    x30 = x29 * x5
    x31 = x24 + x30
    x32 = 3 * x19
    x33 = x3 * (3 * x21 + x29 + x32)
    x34 = x16 + x23
    x35 = x34 * x5
    x36 = x33 + x35
    x37 = x2 * x36 + x3 * (x17 + 4 * x23 + x31)
    x38 = x2 * x34
    x39 = 8 * x25
    x40 = x2 * x31 + x3 * (4 * x28 + x39)
    x41 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x42 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x43 = numpy.pi * x0 * x42
    x44 = x41 * x43
    x45 = -x0 * (a * A[1] + b * B[1])
    x46 = -x45 - B[1]
    x47 = x33 + x38
    x48 = x2 * x22
    x49 = x2 * x29
    x50 = x24 + x49
    x51 = x44 * (x2 * x47 + x3 * (x17 + x23 + 3 * x48 + x50))
    x52 = -x0 * (a * A[2] + b * B[2])
    x53 = -x52 - B[2]
    x54 = x41 * x8
    x55 = x3 * x54
    x56 = x46 ** 2 * x54
    x57 = x55 + x56
    x58 = x16 + x48
    x59 = x2 * x20
    x60 = x2 * x27
    x61 = x26 + x60
    x62 = x2 * x58 + x3 * (x21 + x32 + 2 * x59 + x61)
    x63 = x42 * x8
    x64 = x44 * x53
    x65 = x3 * x63
    x66 = x53 ** 2 * x63
    x67 = x65 + x66
    x68 = x46 * x55
    x69 = 2 * x68
    x70 = x46 * x57
    x71 = x69 + x70
    x72 = x10 * x2
    x73 = x18 * x2
    x74 = x19 + x59
    x75 = x2 * x74 + x3 * (x11 + x14 + x72 + x73)
    x76 = x53 * x63
    x77 = x46 * x54
    x78 = x53 * x65
    x79 = 2 * x78
    x80 = x53 * x67
    x81 = x79 + x80
    x82 = 3 * x55
    x83 = x3 * (3 * x56 + x82)
    x84 = x46 * x71
    x85 = x83 + x84
    x86 = x2 * x9
    x87 = x13 + x73
    x88 = x2 * x87 + x3 * (x18 + x86)
    x89 = 3 * x65
    x90 = x3 * (3 * x66 + x89)
    x91 = x53 * x81
    x92 = x90 + x91
    x93 = -x45 - A[1]
    x94 = x37 * x44
    x95 = x54 * x93
    x96 = x46 * x95
    x97 = x55 + x96
    x98 = x57 * x93
    x99 = x69 + x98
    x100 = x71 * x93
    x101 = x100 + x83
    x102 = 8 * x68
    x103 = x3 * (x102 + 4 * x70) + x85 * x93
    x104 = -x52 - A[2]
    x105 = x104 * x44
    x106 = x104 * x63
    x107 = x106 * x53
    x108 = x107 + x65
    x109 = x104 * x67
    x110 = x109 + x79
    x111 = x104 * x81
    x112 = x111 + x90
    x113 = 8 * x78
    x114 = x104 * x92 + x3 * (x113 + 4 * x80)
    x115 = x54 * x93 ** 2 + x55
    x116 = x3 * (x77 + x95) + x93 * x97
    x117 = x56 + x82
    x118 = x3 * (x117 + 2 * x96) + x93 * x99
    x119 = x101 * x93 + x3 * (x102 + x70 + 3 * x98)
    x120 = x43 * x7
    x121 = x120 * (x103 * x93 + x3 * (4 * x100 + 5 * x83 + x84))
    x122 = x120 * x4
    x123 = numpy.pi * x0 * x41 * x7
    x124 = x123 * x4
    x125 = x104 ** 2 * x63 + x65
    x126 = x104 * x108 + x3 * (x106 + x76)
    x127 = x66 + x89
    x128 = x104 * x110 + x3 * (2 * x107 + x127)
    x129 = x104 * x112 + x3 * (3 * x109 + x113 + x80)
    x130 = x123 * (x104 * x114 + x3 * (4 * x111 + 5 * x90 + x91))
    x131 = -x45 - C[1]
    x132 = x44 * (x2 * x40 + x3 * (5 * x24 + x30 + 4 * x49))
    x133 = x131 * x77
    x134 = x133 + x55
    x135 = x2 * x50 + x3 * (x28 + x39 + 3 * x60)
    x136 = x131 * x54
    x137 = x3 * (x136 + x77)
    x138 = x134 * x46
    x139 = x137 + x138
    x140 = x2 * x61 + x3 * (x15 + 2 * x72)
    x141 = x3 * (x117 + 2 * x133)
    x142 = x139 * x46
    x143 = x141 + x142
    x144 = x13 + x72
    x145 = x144 * x2 + x3 * (x10 + x86)
    x146 = 3 * x137
    x147 = x3 * (3 * x138 + x146 + x71)
    x148 = x143 * x46
    x149 = x147 + x148
    x150 = x13 + x2 ** 2 * x9
    x151 = x131 * x95
    x152 = x151 + x55
    x153 = x134 * x93
    x154 = x137 + x153
    x155 = x139 * x93
    x156 = x141 + x155
    x157 = x143 * x93
    x158 = x147 + x157
    x159 = 4 * x141
    x160 = x149 * x93 + x3 * (4 * x142 + x159 + x85)
    x161 = x120 * x160
    x162 = x120 * x2
    x163 = x123 * x131
    x164 = x152 * x93 + x3 * (x136 + x95)
    x165 = x154 * x93 + x3 * (x133 + x151 + x82 + x96)
    x166 = x156 * x93 + x3 * (x138 + x146 + 2 * x153 + x99)
    x167 = x120 * (x158 * x93 + x3 * (x101 + x142 + 3 * x155 + x159))
    x168 = x120 * x5
    x169 = -x52 - C[2]
    x170 = x169 * x44
    x171 = x169 * x76
    x172 = x171 + x65
    x173 = x169 * x63
    x174 = x3 * (x173 + x76)
    x175 = x172 * x53
    x176 = x174 + x175
    x177 = x3 * (x127 + 2 * x171)
    x178 = x176 * x53
    x179 = x177 + x178
    x180 = 3 * x174
    x181 = x3 * (3 * x175 + x180 + x81)
    x182 = x179 * x53
    x183 = x181 + x182
    x184 = x123 * x2
    x185 = x106 * x169
    x186 = x185 + x65
    x187 = x104 * x172
    x188 = x174 + x187
    x189 = x104 * x176
    x190 = x177 + x189
    x191 = x104 * x179
    x192 = x181 + x191
    x193 = 4 * x177
    x194 = x104 * x183 + x3 * (4 * x178 + x193 + x92)
    x195 = x123 * x194
    x196 = x123 * x5
    x197 = x104 * x186 + x3 * (x106 + x173)
    x198 = x104 * x188 + x3 * (x107 + x171 + x185 + x89)
    x199 = x104 * x190 + x3 * (x110 + x175 + x180 + 2 * x187)
    x200 = x123 * (x104 * x192 + x3 * (x112 + x178 + 3 * x189 + x193))

    # 270 item(s)
    return numpy.array(
        [
            x44 * (x2 * x37 + x3 * (5 * x33 + x35 + 4 * x38 + x40)),
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
            x47 * x63 * x97,
            x47 * x64 * x93,
            x58 * x63 * x99,
            x58 * x76 * x97,
            x58 * x67 * x95,
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
            x105 * x46 * x47,
            x108 * x47 * x54,
            x106 * x57 * x58,
            x108 * x58 * x77,
            x110 * x54 * x58,
            x106 * x71 * x74,
            x108 * x57 * x74,
            x110 * x74 * x77,
            x112 * x54 * x74,
            x106 * x85 * x87,
            x108 * x71 * x87,
            x110 * x57 * x87,
            x112 * x77 * x87,
            x114 * x54 * x87,
            x115 * x36 * x63,
            x116 * x34 * x63,
            x115 * x34 * x76,
            x118 * x22 * x63,
            x116 * x22 * x76,
            x115 * x22 * x67,
            x119 * x20 * x63,
            x118 * x20 * x76,
            x116 * x20 * x67,
            x115 * x20 * x81,
            x121 * x4,
            x119 * x122 * x53,
            x118 * x18 * x67,
            x116 * x18 * x81,
            x115 * x18 * x92,
            x105 * x36 * x93,
            x106 * x34 * x97,
            x108 * x34 * x95,
            x106 * x22 * x99,
            x108 * x22 * x97,
            x110 * x22 * x95,
            x101 * x106 * x20,
            x108 * x20 * x99,
            x110 * x20 * x97,
            x112 * x20 * x95,
            x103 * x104 * x122,
            x101 * x108 * x18,
            x110 * x18 * x99,
            x112 * x18 * x97,
            x114 * x124 * x93,
            x125 * x36 * x54,
            x125 * x34 * x77,
            x126 * x34 * x54,
            x125 * x22 * x57,
            x126 * x22 * x77,
            x128 * x22 * x54,
            x125 * x20 * x71,
            x126 * x20 * x57,
            x128 * x20 * x77,
            x129 * x20 * x54,
            x125 * x18 * x85,
            x126 * x18 * x71,
            x128 * x18 * x57,
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
            x152 * x40 * x63,
            x154 * x50 * x63,
            x152 * x50 * x76,
            x156 * x61 * x63,
            x154 * x61 * x76,
            x152 * x61 * x67,
            x144 * x158 * x63,
            x144 * x156 * x76,
            x144 * x154 * x67,
            x144 * x152 * x81,
            x161 * x2,
            x158 * x162 * x53,
            x156 * x67 * x86,
            x154 * x81 * x86,
            x152 * x86 * x92,
            x105 * x131 * x40,
            x106 * x134 * x50,
            x108 * x136 * x50,
            x106 * x139 * x61,
            x108 * x134 * x61,
            x110 * x136 * x61,
            x106 * x143 * x144,
            x108 * x139 * x144,
            x110 * x134 * x144,
            x112 * x136 * x144,
            x104 * x149 * x162,
            x108 * x143 * x86,
            x110 * x139 * x86,
            x112 * x134 * x86,
            x114 * x163 * x2,
            x164 * x31 * x63,
            x165 * x29 * x63,
            x164 * x29 * x76,
            x166 * x27 * x63,
            x165 * x27 * x76,
            x164 * x27 * x67,
            x167 * x5,
            x166 * x168 * x53,
            x10 * x165 * x67,
            x10 * x164 * x81,
            x120 * (x160 * x93 + x3 * (x103 + 5 * x147 + x148 + 4 * x157)),
            x167 * x53,
            x166 * x67 * x9,
            x165 * x81 * x9,
            x164 * x9 * x92,
            x106 * x152 * x31,
            x106 * x154 * x29,
            x108 * x152 * x29,
            x106 * x156 * x27,
            x108 * x154 * x27,
            x110 * x152 * x27,
            x104 * x158 * x168,
            x10 * x108 * x156,
            x10 * x110 * x154,
            x10 * x112 * x152,
            x104 * x161,
            x108 * x158 * x9,
            x110 * x156 * x9,
            x112 * x154 * x9,
            x114 * x152 * x9,
            x125 * x136 * x31,
            x125 * x134 * x29,
            x126 * x136 * x29,
            x125 * x139 * x27,
            x126 * x134 * x27,
            x128 * x136 * x27,
            x10 * x125 * x143,
            x10 * x126 * x139,
            x10 * x128 * x134,
            x129 * x163 * x5,
            x125 * x149 * x9,
            x126 * x143 * x9,
            x128 * x139 * x9,
            x129 * x134 * x9,
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
            x170 * x40 * x93,
            x173 * x50 * x97,
            x172 * x50 * x95,
            x173 * x61 * x99,
            x172 * x61 * x97,
            x176 * x61 * x95,
            x101 * x144 * x173,
            x144 * x172 * x99,
            x144 * x176 * x97,
            x144 * x179 * x95,
            x103 * x162 * x169,
            x101 * x172 * x86,
            x176 * x86 * x99,
            x179 * x86 * x97,
            x183 * x184 * x93,
            x186 * x40 * x54,
            x186 * x50 * x77,
            x188 * x50 * x54,
            x186 * x57 * x61,
            x188 * x61 * x77,
            x190 * x54 * x61,
            x144 * x186 * x71,
            x144 * x188 * x57,
            x144 * x190 * x77,
            x144 * x192 * x54,
            x186 * x85 * x86,
            x188 * x71 * x86,
            x190 * x57 * x86,
            x184 * x192 * x46,
            x195 * x2,
            x115 * x173 * x31,
            x116 * x173 * x29,
            x115 * x172 * x29,
            x118 * x173 * x27,
            x116 * x172 * x27,
            x115 * x176 * x27,
            x119 * x168 * x169,
            x10 * x118 * x172,
            x10 * x116 * x176,
            x10 * x115 * x179,
            x121 * x169,
            x119 * x172 * x9,
            x118 * x176 * x9,
            x116 * x179 * x9,
            x115 * x183 * x9,
            x186 * x31 * x95,
            x186 * x29 * x97,
            x188 * x29 * x95,
            x186 * x27 * x99,
            x188 * x27 * x97,
            x190 * x27 * x95,
            x10 * x101 * x186,
            x10 * x188 * x99,
            x10 * x190 * x97,
            x192 * x196 * x93,
            x103 * x186 * x9,
            x101 * x188 * x9,
            x190 * x9 * x99,
            x192 * x9 * x97,
            x195 * x93,
            x197 * x31 * x54,
            x197 * x29 * x77,
            x198 * x29 * x54,
            x197 * x27 * x57,
            x198 * x27 * x77,
            x199 * x27 * x54,
            x10 * x197 * x71,
            x10 * x198 * x57,
            x196 * x199 * x46,
            x200 * x5,
            x197 * x85 * x9,
            x198 * x71 * x9,
            x199 * x57 * x9,
            x200 * x46,
            x123 * (x104 * x194 + x3 * (x114 + 5 * x181 + x182 + 4 * x191)),
        ]
    )


def dipole3d_30(a, A, b, B, C):
    """Cartesian 3D (fs) dipole moment integrals.
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
    x8 = -x7 - A[0]
    x9 = x5 * x8 ** 2
    x10 = -x7 - C[0]
    x11 = x5 * x8
    x12 = x10 * x11
    x13 = x12 + x6
    x14 = x0 * (x10 * x5 + x11) + x13 * x8
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
    x26 = x20 ** 2 * x24
    x27 = x25 + x26
    x28 = x16 * x4
    x29 = x18 * x23
    x30 = x0 * x28
    x31 = x23 ** 2 * x28
    x32 = x30 + x31
    x33 = x17 * x3
    x34 = x33 * (2 * x20 * x25 + x20 * x27)
    x35 = x23 * x33
    x36 = numpy.pi * x1 * x15 * x3
    x37 = x32 * x36
    x38 = x36 * (2 * x23 * x30 + x23 * x32)
    x39 = -x19 - C[1]
    x40 = x6 + x9
    x41 = x18 * (2 * x0 * x11 + x40 * x8)
    x42 = x20 * x24
    x43 = x39 * x42
    x44 = x25 + x43
    x45 = x0 * (x24 * x39 + x42) + x20 * x44
    x46 = x33 * x45
    x47 = -x22 - C[2]
    x48 = x23 * x28
    x49 = x47 * x48
    x50 = x30 + x49
    x51 = x0 * (x28 * x47 + x48) + x23 * x50
    x52 = x36 * x51

    # 30 item(s)
    return numpy.array(
        [
            x18 * (x0 * (2 * x12 + 3 * x6 + x9) + x14 * x8),
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
            x46 * x8,
            x35 * x44 * x8,
            x37 * x39 * x8,
            x33 * (x0 * (3 * x25 + x26 + 2 * x43) + x20 * x45),
            x23 * x46,
            x32 * x44 * x5,
            x38 * x39,
            x41 * x47,
            x18 * x20 * x40 * x47,
            x24 * x40 * x50,
            x27 * x33 * x47 * x8,
            x20 * x36 * x50 * x8,
            x52 * x8,
            x34 * x47,
            x27 * x5 * x50,
            x20 * x52,
            x36 * (x0 * (3 * x30 + x31 + 2 * x49) + x23 * x51),
        ]
    )


def dipole3d_31(a, A, b, B, C):
    """Cartesian 3D (fp) dipole moment integrals.
    The origin is at C.

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
    x11 = x2 * x5
    x12 = x11 * x6
    x13 = x10 * x12
    x14 = -x1 - C[0]
    x15 = x12 * x14
    x16 = x10 * x7
    x17 = x14 * x16
    x18 = x14 * x7
    x19 = x3 * (x16 + x18)
    x20 = x17 + x8
    x21 = x2 * x20
    x22 = x19 + x21
    x23 = x2 * x22 + x3 * (x13 + x15 + x17 + x9)
    x24 = x15 + x8
    x25 = x2 * x24 + x3 * (x12 + x18)
    x26 = x13 + x8
    x27 = x2 * x26 + x3 * (x12 + x16)
    x28 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x29 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x30 = numpy.pi * x0 * x29
    x31 = x28 * x30
    x32 = -x0 * (a * A[1] + b * B[1])
    x33 = -x32 - B[1]
    x34 = x2 ** 2 * x7
    x35 = x34 + x9
    x36 = x31 * (x2 * x25 + x3 * (2 * x15 + x35))
    x37 = -x0 * (a * A[2] + b * B[2])
    x38 = -x37 - B[2]
    x39 = -x32 - A[1]
    x40 = x23 * x31
    x41 = x3 * x6
    x42 = x28 * x41
    x43 = x28 * x6
    x44 = x39 * x43
    x45 = x33 * x44
    x46 = x42 + x45
    x47 = x29 * x6
    x48 = x25 * x31
    x49 = -x37 - A[2]
    x50 = x29 * x41
    x51 = x47 * x49
    x52 = x38 * x51
    x53 = x50 + x52
    x54 = x39 ** 2 * x43
    x55 = x42 + x54
    x56 = x33 * x43
    x57 = x3 * (x44 + x56) + x39 * x46
    x58 = x38 * x47
    x59 = x31 * x49
    x60 = x47 * x49 ** 2
    x61 = x50 + x60
    x62 = x3 * (x51 + x58) + x49 * x53
    x63 = 2 * x39 * x42 + x39 * x55
    x64 = 3 * x42
    x65 = x54 + x64
    x66 = x30 * x5
    x67 = x66 * (x3 * (2 * x45 + x65) + x39 * x57)
    x68 = x63 * x66
    x69 = x49 * x66
    x70 = numpy.pi * x0 * x28
    x71 = x5 * x70
    x72 = x39 * x71
    x73 = 2 * x49 * x50 + x49 * x61
    x74 = x71 * x73
    x75 = 3 * x50
    x76 = x60 + x75
    x77 = x71 * (x3 * (2 * x52 + x76) + x49 * x62)
    x78 = -x32 - C[1]
    x79 = x31 * (x2 * x27 + x3 * (2 * x13 + x35))
    x80 = x56 * x78
    x81 = x42 + x80
    x82 = x34 + x8
    x83 = 2 * x12 * x3 + x2 * x82
    x84 = x31 * x83
    x85 = x44 * x78
    x86 = x42 + x85
    x87 = x43 * x78
    x88 = x3 * (x56 + x87)
    x89 = x39 * x81
    x90 = x88 + x89
    x91 = x3 * (x44 + x87) + x39 * x86
    x92 = x3 * (x45 + x64 + x80 + x85) + x39 * x90
    x93 = x11 * x30
    x94 = x11 * x70
    x95 = x66 * (x3 * (x65 + 2 * x85) + x39 * x91)
    x96 = -x37 - C[2]
    x97 = x58 * x96
    x98 = x50 + x97
    x99 = x47 * x96
    x100 = x51 * x96
    x101 = x100 + x50
    x102 = x3 * (x58 + x99)
    x103 = x49 * x98
    x104 = x102 + x103
    x105 = x101 * x49 + x3 * (x51 + x99)
    x106 = x104 * x49 + x3 * (x100 + x52 + x75 + x97)
    x107 = x71 * (x105 * x49 + x3 * (2 * x100 + x76))

    # 90 item(s)
    return numpy.array(
        [
            x31 * (x2 * x23 + x3 * (2 * x19 + 2 * x21 + x25 + x27)),
            x33 * x36,
            x36 * x38,
            x39 * x40,
            x25 * x46 * x47,
            x38 * x39 * x48,
            x40 * x49,
            x33 * x48 * x49,
            x25 * x43 * x53,
            x22 * x47 * x55,
            x24 * x47 * x57,
            x24 * x55 * x58,
            x22 * x39 * x59,
            x24 * x46 * x51,
            x24 * x44 * x53,
            x22 * x43 * x61,
            x24 * x56 * x61,
            x24 * x43 * x62,
            x20 * x47 * x63,
            x14 * x67,
            x14 * x38 * x68,
            x20 * x51 * x55,
            x14 * x57 * x69,
            x18 * x53 * x55,
            x20 * x44 * x61,
            x18 * x46 * x61,
            x14 * x62 * x72,
            x20 * x43 * x73,
            x14 * x33 * x74,
            x14 * x77,
            x78 * x79,
            x47 * x81 * x83,
            x38 * x78 * x84,
            x27 * x47 * x86,
            x47 * x82 * x90,
            x58 * x82 * x86,
            x27 * x59 * x78,
            x51 * x81 * x82,
            x53 * x82 * x87,
            x26 * x47 * x91,
            x92 * x93,
            x38 * x91 * x93,
            x26 * x51 * x86,
            x49 * x90 * x93,
            x12 * x53 * x86,
            x26 * x61 * x87,
            x12 * x61 * x81,
            x62 * x78 * x94,
            x10 * x95,
            x66 * (x3 * (x57 + 2 * x88 + 2 * x89 + x91) + x39 * x92),
            x38 * x95,
            x10 * x69 * x91,
            x69 * x92,
            x53 * x7 * x91,
            x16 * x61 * x86,
            x61 * x7 * x90,
            x62 * x7 * x86,
            x10 * x74 * x78,
            x7 * x73 * x81,
            x77 * x78,
            x79 * x96,
            x33 * x84 * x96,
            x43 * x83 * x98,
            x27 * x31 * x39 * x96,
            x46 * x82 * x99,
            x44 * x82 * x98,
            x101 * x27 * x43,
            x101 * x56 * x82,
            x104 * x43 * x82,
            x26 * x55 * x99,
            x57 * x93 * x96,
            x12 * x55 * x98,
            x101 * x26 * x44,
            x101 * x12 * x46,
            x104 * x39 * x94,
            x105 * x26 * x43,
            x105 * x33 * x94,
            x106 * x94,
            x10 * x68 * x96,
            x67 * x96,
            x63 * x7 * x98,
            x101 * x16 * x55,
            x101 * x57 * x7,
            x104 * x55 * x7,
            x10 * x105 * x72,
            x105 * x46 * x7,
            x106 * x72,
            x10 * x107,
            x107 * x33,
            x71 * (x106 * x49 + x3 * (2 * x102 + 2 * x103 + x105 + x62)),
        ]
    )


def dipole3d_32(a, A, b, B, C):
    """Cartesian 3D (fd) dipole moment integrals.
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
    x12 = x5 ** 2 * x9
    x13 = x3 * x9
    x14 = 3 * x13
    x15 = x12 + x14
    x16 = x3 * (2 * x11 + x15)
    x17 = x4 * x9
    x18 = x3 * (x10 + x17)
    x19 = x11 + x13
    x20 = x19 * x5
    x21 = x18 + x20
    x22 = x2 * x21
    x23 = x16 + x22
    x24 = x19 * x2
    x25 = 2 * x24
    x26 = x13 * x5
    x27 = x12 + x13
    x28 = x2 * x27
    x29 = 2 * x26 + x28
    x30 = x2 * x23 + x3 * (3 * x18 + x20 + x25 + x29)
    x31 = x10 * x2
    x32 = x17 * x2
    x33 = x3 * (x11 + x14 + x31 + x32)
    x34 = x18 + x24
    x35 = x2 * x34
    x36 = 2 * x31
    x37 = x2 * x29 + x3 * (x15 + x36)
    x38 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x39 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x40 = numpy.pi * x0 * x39
    x41 = x38 * x40
    x42 = -x0 * (a * A[1] + b * B[1])
    x43 = -x42 - B[1]
    x44 = x33 + x35
    x45 = x2 * x9
    x46 = x13 + x32
    x47 = x2 * x46 + x3 * (x17 + x45)
    x48 = x3 * (x10 + x45)
    x49 = x13 + x31
    x50 = x2 * x49
    x51 = x48 + x50
    x52 = x41 * (x2 * x44 + x3 * (2 * x18 + x25 + x47 + x51))
    x53 = -x0 * (a * A[2] + b * B[2])
    x54 = -x53 - B[2]
    x55 = x38 * x8
    x56 = x3 * x55
    x57 = x43 ** 2 * x55
    x58 = x56 + x57
    x59 = x2 ** 2 * x9
    x60 = x14 + x59
    x61 = x2 * x47 + x3 * (2 * x32 + x60)
    x62 = x39 * x8
    x63 = x41 * x54
    x64 = x3 * x62
    x65 = x54 ** 2 * x62
    x66 = x64 + x65
    x67 = -x42 - A[1]
    x68 = x30 * x41
    x69 = x55 * x67
    x70 = x43 * x69
    x71 = x56 + x70
    x72 = 2 * x56
    x73 = x58 * x67
    x74 = x43 * x72 + x73
    x75 = x54 * x62
    x76 = -x53 - A[2]
    x77 = x41 * x76
    x78 = x62 * x76
    x79 = x54 * x78
    x80 = x64 + x79
    x81 = x43 * x55
    x82 = 2 * x64
    x83 = x66 * x76
    x84 = x54 * x82 + x83
    x85 = x55 * x67 ** 2
    x86 = x56 + x85
    x87 = x3 * (x69 + x81)
    x88 = x67 * x71
    x89 = x87 + x88
    x90 = 3 * x56
    x91 = 2 * x70 + x90
    x92 = x3 * (x57 + x91) + x67 * x74
    x93 = x62 * x76 ** 2
    x94 = x64 + x93
    x95 = x3 * (x75 + x78)
    x96 = x76 * x80
    x97 = x95 + x96
    x98 = 3 * x64
    x99 = 2 * x79 + x98
    x100 = x3 * (x65 + x99) + x76 * x84
    x101 = x67 * x72 + x67 * x86
    x102 = x3 * (x85 + x91) + x67 * x89
    x103 = x40 * x7
    x104 = x103 * (x3 * (4 * x43 * x56 + 2 * x73 + 2 * x87 + 2 * x88) + x67 * x92)
    x105 = x103 * x4
    x106 = numpy.pi * x0 * x38 * x7
    x107 = x106 * x4
    x108 = x76 * x82 + x76 * x94
    x109 = x3 * (x93 + x99) + x76 * x97
    x110 = x106 * (x100 * x76 + x3 * (4 * x54 * x64 + 2 * x83 + 2 * x95 + 2 * x96))
    x111 = -x42 - C[1]
    x112 = x41 * (x2 * x37 + x3 * (4 * x26 + 2 * x28 + 2 * x48 + 2 * x50))
    x113 = x111 * x81
    x114 = x113 + x56
    x115 = x2 * x51 + x3 * (x36 + x60)
    x116 = x111 * x55
    x117 = x3 * (x116 + x81)
    x118 = x114 * x43
    x119 = x117 + x118
    x120 = x13 + x59
    x121 = x120 * x2 + 2 * x13 * x2
    x122 = x111 * x69
    x123 = x122 + x56
    x124 = x114 * x67
    x125 = x117 + x124
    x126 = x3 * (2 * x113 + x57 + x90)
    x127 = x119 * x67
    x128 = x126 + x127
    x129 = x123 * x67 + x3 * (x116 + x69)
    x130 = x3 * (x113 + x122 + x70 + x90)
    x131 = x125 * x67
    x132 = x130 + x131
    x133 = 2 * x124
    x134 = x128 * x67 + x3 * (3 * x117 + x118 + x133 + x74)
    x135 = x103 * x134
    x136 = x103 * x2
    x137 = x106 * x111
    x138 = x129 * x67 + x3 * (2 * x122 + x85 + x90)
    x139 = x103 * (x132 * x67 + x3 * (2 * x117 + x129 + x133 + x89))
    x140 = x103 * x5
    x141 = -x53 - C[2]
    x142 = x141 * x41
    x143 = x141 * x75
    x144 = x143 + x64
    x145 = x141 * x62
    x146 = x3 * (x145 + x75)
    x147 = x144 * x54
    x148 = x146 + x147
    x149 = x141 * x78
    x150 = x149 + x64
    x151 = x144 * x76
    x152 = x146 + x151
    x153 = x3 * (2 * x143 + x65 + x98)
    x154 = x148 * x76
    x155 = x153 + x154
    x156 = x106 * x2
    x157 = x150 * x76 + x3 * (x145 + x78)
    x158 = x3 * (x143 + x149 + x79 + x98)
    x159 = x152 * x76
    x160 = x158 + x159
    x161 = 2 * x151
    x162 = x155 * x76 + x3 * (3 * x146 + x147 + x161 + x84)
    x163 = x106 * x162
    x164 = x106 * x5
    x165 = x157 * x76 + x3 * (2 * x149 + x93 + x98)
    x166 = x106 * (x160 * x76 + x3 * (2 * x146 + x157 + x161 + x97))

    # 180 item(s)
    return numpy.array(
        [
            x41 * (x2 * x30 + x3 * (2 * x16 + 2 * x22 + 2 * x33 + 2 * x35 + x37)),
            x43 * x52,
            x52 * x54,
            x58 * x61 * x62,
            x43 * x61 * x63,
            x55 * x61 * x66,
            x67 * x68,
            x44 * x62 * x71,
            x44 * x63 * x67,
            x47 * x62 * x74,
            x47 * x71 * x75,
            x47 * x66 * x69,
            x68 * x76,
            x43 * x44 * x77,
            x44 * x55 * x80,
            x47 * x58 * x78,
            x47 * x80 * x81,
            x47 * x55 * x84,
            x23 * x62 * x86,
            x34 * x62 * x89,
            x34 * x75 * x86,
            x46 * x62 * x92,
            x46 * x75 * x89,
            x46 * x66 * x86,
            x23 * x67 * x77,
            x34 * x71 * x78,
            x34 * x69 * x80,
            x46 * x74 * x78,
            x46 * x71 * x80,
            x46 * x69 * x84,
            x23 * x55 * x94,
            x34 * x81 * x94,
            x34 * x55 * x97,
            x46 * x58 * x94,
            x46 * x81 * x97,
            x100 * x46 * x55,
            x101 * x21 * x62,
            x102 * x19 * x62,
            x101 * x19 * x75,
            x104 * x4,
            x102 * x105 * x54,
            x101 * x17 * x66,
            x21 * x78 * x86,
            x19 * x78 * x89,
            x19 * x80 * x86,
            x105 * x76 * x92,
            x17 * x80 * x89,
            x17 * x84 * x86,
            x21 * x69 * x94,
            x19 * x71 * x94,
            x19 * x69 * x97,
            x17 * x74 * x94,
            x17 * x71 * x97,
            x100 * x107 * x67,
            x108 * x21 * x55,
            x108 * x19 * x81,
            x109 * x19 * x55,
            x108 * x17 * x58,
            x107 * x109 * x43,
            x110 * x4,
            x111 * x112,
            x114 * x115 * x62,
            x111 * x115 * x63,
            x119 * x121 * x62,
            x114 * x121 * x75,
            x116 * x121 * x66,
            x123 * x37 * x62,
            x125 * x51 * x62,
            x123 * x51 * x75,
            x120 * x128 * x62,
            x120 * x125 * x75,
            x120 * x123 * x66,
            x111 * x37 * x77,
            x114 * x51 * x78,
            x116 * x51 * x80,
            x119 * x120 * x78,
            x114 * x120 * x80,
            x116 * x120 * x84,
            x129 * x29 * x62,
            x132 * x49 * x62,
            x129 * x49 * x75,
            x135 * x2,
            x132 * x136 * x54,
            x129 * x45 * x66,
            x123 * x29 * x78,
            x125 * x49 * x78,
            x123 * x49 * x80,
            x128 * x136 * x76,
            x125 * x45 * x80,
            x123 * x45 * x84,
            x116 * x29 * x94,
            x114 * x49 * x94,
            x116 * x49 * x97,
            x119 * x45 * x94,
            x114 * x45 * x97,
            x100 * x137 * x2,
            x138 * x27 * x62,
            x139 * x5,
            x138 * x140 * x54,
            x103 * (x134 * x67 + x3 * (2 * x126 + 2 * x127 + 2 * x130 + 2 * x131 + x92)),
            x139 * x54,
            x138 * x66 * x9,
            x129 * x27 * x78,
            x132 * x140 * x76,
            x10 * x129 * x80,
            x135 * x76,
            x132 * x80 * x9,
            x129 * x84 * x9,
            x123 * x27 * x94,
            x10 * x125 * x94,
            x10 * x123 * x97,
            x128 * x9 * x94,
            x125 * x9 * x97,
            x100 * x123 * x9,
            x108 * x116 * x27,
            x10 * x108 * x114,
            x109 * x137 * x5,
            x108 * x119 * x9,
            x109 * x114 * x9,
            x110 * x111,
            x112 * x141,
            x115 * x142 * x43,
            x115 * x144 * x55,
            x121 * x145 * x58,
            x121 * x144 * x81,
            x121 * x148 * x55,
            x142 * x37 * x67,
            x145 * x51 * x71,
            x144 * x51 * x69,
            x120 * x145 * x74,
            x120 * x144 * x71,
            x120 * x148 * x69,
            x150 * x37 * x55,
            x150 * x51 * x81,
            x152 * x51 * x55,
            x120 * x150 * x58,
            x120 * x152 * x81,
            x120 * x155 * x55,
            x145 * x29 * x86,
            x145 * x49 * x89,
            x144 * x49 * x86,
            x136 * x141 * x92,
            x144 * x45 * x89,
            x148 * x45 * x86,
            x150 * x29 * x69,
            x150 * x49 * x71,
            x152 * x49 * x69,
            x150 * x45 * x74,
            x152 * x45 * x71,
            x155 * x156 * x67,
            x157 * x29 * x55,
            x157 * x49 * x81,
            x160 * x49 * x55,
            x157 * x45 * x58,
            x156 * x160 * x43,
            x163 * x2,
            x101 * x145 * x27,
            x102 * x140 * x141,
            x10 * x101 * x144,
            x104 * x141,
            x102 * x144 * x9,
            x101 * x148 * x9,
            x150 * x27 * x86,
            x10 * x150 * x89,
            x10 * x152 * x86,
            x150 * x9 * x92,
            x152 * x89 * x9,
            x155 * x86 * x9,
            x157 * x27 * x69,
            x10 * x157 * x71,
            x160 * x164 * x67,
            x157 * x74 * x9,
            x160 * x71 * x9,
            x163 * x67,
            x165 * x27 * x55,
            x164 * x165 * x43,
            x166 * x5,
            x165 * x58 * x9,
            x166 * x43,
            x106 * (x162 * x76 + x3 * (x100 + 2 * x153 + 2 * x154 + 2 * x158 + 2 * x159)),
        ]
    )


def dipole3d_33(a, A, b, B, C):
    """Cartesian 3D (ff) dipole moment integrals.
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
    x13 = 3 * x12
    x14 = x3 * x8
    x15 = x10 * x9
    x16 = x14 + x15
    x17 = x16 * x4
    x18 = x14 * x4
    x19 = 2 * x18
    x20 = x4 ** 2 * x8
    x21 = x14 + x20
    x22 = x21 * x4
    x23 = x19 + x22
    x24 = x3 * (x13 + 3 * x17 + x23)
    x25 = 3 * x14
    x26 = x20 + x25
    x27 = x3 * (2 * x15 + x26)
    x28 = x12 + x17
    x29 = x28 * x4
    x30 = x27 + x29
    x31 = x2 * x30
    x32 = x24 + x31
    x33 = x2 * x28
    x34 = x3 * (3 * x20 + x25)
    x35 = x2 * x23
    x36 = x34 + x35
    x37 = x2 * x32 + x3 * (4 * x27 + x29 + 3 * x33 + x36)
    x38 = x16 * x2
    x39 = 2 * x38
    x40 = x2 * x21
    x41 = x19 + x40
    x42 = x3 * (x13 + x17 + x39 + x41)
    x43 = x27 + x33
    x44 = x2 * x43
    x45 = x2 * x36 + x3 * (8 * x18 + x22 + 3 * x40)
    x46 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x47 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x48 = numpy.pi * x0 * x47
    x49 = x46 * x48
    x50 = -x0 * (a * A[1] + b * B[1])
    x51 = -x50 - B[1]
    x52 = x42 + x44
    x53 = x2 * x9
    x54 = x11 * x2
    x55 = x3 * (x15 + x25 + x53 + x54)
    x56 = x12 + x38
    x57 = x2 * x56
    x58 = 2 * x53
    x59 = x3 * (x26 + x58)
    x60 = x2 * x41
    x61 = x59 + x60
    x62 = x49 * (x2 * x52 + x3 * (2 * x27 + 2 * x33 + 2 * x55 + 2 * x57 + x61))
    x63 = -x0 * (a * A[2] + b * B[2])
    x64 = -x63 - B[2]
    x65 = x46 * x7
    x66 = x3 * x65
    x67 = x51 ** 2 * x65
    x68 = x66 + x67
    x69 = x55 + x57
    x70 = x2 * x8
    x71 = x14 + x54
    x72 = x2 * x71 + x3 * (x11 + x70)
    x73 = x3 * (x70 + x9)
    x74 = x14 + x53
    x75 = x2 * x74
    x76 = x73 + x75
    x77 = x2 * x69 + x3 * (2 * x12 + x39 + x72 + x76)
    x78 = x47 * x7
    x79 = x49 * x64
    x80 = x3 * x78
    x81 = x64 ** 2 * x78
    x82 = x80 + x81
    x83 = x51 * x66
    x84 = 2 * x83
    x85 = x51 * x68
    x86 = x84 + x85
    x87 = x2 ** 2 * x8
    x88 = x25 + x87
    x89 = x2 * x72 + x3 * (2 * x54 + x88)
    x90 = x64 * x78
    x91 = x51 * x65
    x92 = x64 * x80
    x93 = 2 * x92
    x94 = x64 * x82
    x95 = x93 + x94
    x96 = -x50 - A[1]
    x97 = x37 * x49
    x98 = x65 * x96
    x99 = x51 * x98
    x100 = x66 + x99
    x101 = x68 * x96
    x102 = x101 + x84
    x103 = 3 * x66
    x104 = x3 * (x103 + 3 * x67)
    x105 = x86 * x96
    x106 = x104 + x105
    x107 = -x63 - A[2]
    x108 = x107 * x49
    x109 = x107 * x78
    x110 = x109 * x64
    x111 = x110 + x80
    x112 = x107 * x82
    x113 = x112 + x93
    x114 = 3 * x80
    x115 = x3 * (x114 + 3 * x81)
    x116 = x107 * x95
    x117 = x115 + x116
    x118 = x65 * x96 ** 2
    x119 = x118 + x66
    x120 = x3 * (x91 + x98)
    x121 = x100 * x96
    x122 = x120 + x121
    x123 = x103 + 2 * x99
    x124 = x3 * (x123 + x67)
    x125 = x102 * x96
    x126 = x124 + x125
    x127 = x106 * x96 + x3 * (3 * x101 + 8 * x83 + x85)
    x128 = x107 ** 2 * x78
    x129 = x128 + x80
    x130 = x3 * (x109 + x90)
    x131 = x107 * x111
    x132 = x130 + x131
    x133 = 2 * x110 + x114
    x134 = x3 * (x133 + x81)
    x135 = x107 * x113
    x136 = x134 + x135
    x137 = x107 * x117 + x3 * (3 * x112 + 8 * x92 + x94)
    x138 = x119 * x96 + 2 * x66 * x96
    x139 = x122 * x96 + x3 * (x118 + x123)
    x140 = x126 * x96 + x3 * (2 * x101 + 2 * x120 + 2 * x121 + 4 * x83)
    x141 = x48 * x6
    x142 = x141 * (x127 * x96 + x3 * (2 * x104 + 2 * x105 + 3 * x124 + 3 * x125))
    x143 = x10 * x141
    x144 = numpy.pi * x0 * x46 * x6
    x145 = x10 * x144
    x146 = x107 * x129 + 2 * x107 * x80
    x147 = x107 * x132 + x3 * (x128 + x133)
    x148 = x107 * x136 + x3 * (2 * x112 + 2 * x130 + 2 * x131 + 4 * x92)
    x149 = x144 * (x107 * x137 + x3 * (2 * x115 + 2 * x116 + 3 * x134 + 3 * x135))
    x150 = -x50 - C[1]
    x151 = x49 * (x2 * x45 + x3 * (2 * x34 + 2 * x35 + 3 * x59 + 3 * x60))
    x152 = x150 * x91
    x153 = x152 + x66
    x154 = x2 * x61 + x3 * (4 * x18 + 2 * x40 + 2 * x73 + 2 * x75)
    x155 = x150 * x65
    x156 = x3 * (x155 + x91)
    x157 = x153 * x51
    x158 = x156 + x157
    x159 = x2 * x76 + x3 * (x58 + x88)
    x160 = x3 * (x103 + 2 * x152 + x67)
    x161 = x158 * x51
    x162 = x160 + x161
    x163 = x14 + x87
    x164 = 2 * x14 * x2 + x163 * x2
    x165 = x150 * x98
    x166 = x165 + x66
    x167 = x153 * x96
    x168 = x156 + x167
    x169 = x158 * x96
    x170 = x160 + x169
    x171 = 3 * x156
    x172 = x3 * (3 * x157 + x171 + x86)
    x173 = x162 * x96
    x174 = x172 + x173
    x175 = x166 * x96 + x3 * (x155 + x98)
    x176 = x3 * (x103 + x152 + x165 + x99)
    x177 = x168 * x96
    x178 = x176 + x177
    x179 = x170 * x96
    x180 = 2 * x167
    x181 = x3 * (x102 + x157 + x171 + x180)
    x182 = x179 + x181
    x183 = x174 * x96 + x3 * (x106 + 4 * x160 + x161 + 3 * x169)
    x184 = x141 * x183
    x185 = x141 * x2
    x186 = x144 * x150
    x187 = x175 * x96 + x3 * (x103 + x118 + 2 * x165)
    x188 = x178 * x96 + x3 * (x122 + 2 * x156 + x175 + x180)
    x189 = x141 * (x182 * x96 + x3 * (x126 + 2 * x160 + 2 * x169 + 2 * x176 + 2 * x177))
    x190 = x141 * x4
    x191 = -x63 - C[2]
    x192 = x191 * x49
    x193 = x191 * x90
    x194 = x193 + x80
    x195 = x191 * x78
    x196 = x3 * (x195 + x90)
    x197 = x194 * x64
    x198 = x196 + x197
    x199 = x3 * (x114 + 2 * x193 + x81)
    x200 = x198 * x64
    x201 = x199 + x200
    x202 = x109 * x191
    x203 = x202 + x80
    x204 = x107 * x194
    x205 = x196 + x204
    x206 = x107 * x198
    x207 = x199 + x206
    x208 = 3 * x196
    x209 = x3 * (3 * x197 + x208 + x95)
    x210 = x107 * x201
    x211 = x209 + x210
    x212 = x144 * x2
    x213 = x107 * x203 + x3 * (x109 + x195)
    x214 = x3 * (x110 + x114 + x193 + x202)
    x215 = x107 * x205
    x216 = x214 + x215
    x217 = x107 * x207
    x218 = 2 * x204
    x219 = x3 * (x113 + x197 + x208 + x218)
    x220 = x217 + x219
    x221 = x107 * x211 + x3 * (x117 + 4 * x199 + x200 + 3 * x206)
    x222 = x144 * x221
    x223 = x144 * x4
    x224 = x107 * x213 + x3 * (x114 + x128 + 2 * x202)
    x225 = x107 * x216 + x3 * (x132 + 2 * x196 + x213 + x218)
    x226 = x144 * (x107 * x220 + x3 * (x136 + 2 * x199 + 2 * x206 + 2 * x214 + 2 * x215))

    # 300 item(s)
    return numpy.array(
        [
            x49 * (x2 * x37 + x3 * (2 * x24 + 2 * x31 + 3 * x42 + 3 * x44 + x45)),
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
            x100 * x52 * x78,
            x52 * x79 * x96,
            x102 * x69 * x78,
            x100 * x69 * x90,
            x69 * x82 * x98,
            x106 * x72 * x78,
            x102 * x72 * x90,
            x100 * x72 * x82,
            x72 * x95 * x98,
            x107 * x97,
            x108 * x51 * x52,
            x111 * x52 * x65,
            x109 * x68 * x69,
            x111 * x69 * x91,
            x113 * x65 * x69,
            x109 * x72 * x86,
            x111 * x68 * x72,
            x113 * x72 * x91,
            x117 * x65 * x72,
            x119 * x32 * x78,
            x122 * x43 * x78,
            x119 * x43 * x90,
            x126 * x56 * x78,
            x122 * x56 * x90,
            x119 * x56 * x82,
            x127 * x71 * x78,
            x126 * x71 * x90,
            x122 * x71 * x82,
            x119 * x71 * x95,
            x108 * x32 * x96,
            x100 * x109 * x43,
            x111 * x43 * x98,
            x102 * x109 * x56,
            x100 * x111 * x56,
            x113 * x56 * x98,
            x106 * x109 * x71,
            x102 * x111 * x71,
            x100 * x113 * x71,
            x117 * x71 * x98,
            x129 * x32 * x65,
            x129 * x43 * x91,
            x132 * x43 * x65,
            x129 * x56 * x68,
            x132 * x56 * x91,
            x136 * x56 * x65,
            x129 * x71 * x86,
            x132 * x68 * x71,
            x136 * x71 * x91,
            x137 * x65 * x71,
            x138 * x30 * x78,
            x139 * x28 * x78,
            x138 * x28 * x90,
            x140 * x16 * x78,
            x139 * x16 * x90,
            x138 * x16 * x82,
            x10 * x142,
            x140 * x143 * x64,
            x11 * x139 * x82,
            x11 * x138 * x95,
            x109 * x119 * x30,
            x109 * x122 * x28,
            x111 * x119 * x28,
            x109 * x126 * x16,
            x111 * x122 * x16,
            x113 * x119 * x16,
            x107 * x127 * x143,
            x11 * x111 * x126,
            x11 * x113 * x122,
            x11 * x117 * x119,
            x129 * x30 * x98,
            x100 * x129 * x28,
            x132 * x28 * x98,
            x102 * x129 * x16,
            x100 * x132 * x16,
            x136 * x16 * x98,
            x106 * x11 * x129,
            x102 * x11 * x132,
            x100 * x11 * x136,
            x137 * x145 * x96,
            x146 * x30 * x65,
            x146 * x28 * x91,
            x147 * x28 * x65,
            x146 * x16 * x68,
            x147 * x16 * x91,
            x148 * x16 * x65,
            x11 * x146 * x86,
            x11 * x147 * x68,
            x145 * x148 * x51,
            x10 * x149,
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
            x166 * x45 * x78,
            x168 * x61 * x78,
            x166 * x61 * x90,
            x170 * x76 * x78,
            x168 * x76 * x90,
            x166 * x76 * x82,
            x163 * x174 * x78,
            x163 * x170 * x90,
            x163 * x168 * x82,
            x163 * x166 * x95,
            x108 * x150 * x45,
            x109 * x153 * x61,
            x111 * x155 * x61,
            x109 * x158 * x76,
            x111 * x153 * x76,
            x113 * x155 * x76,
            x109 * x162 * x163,
            x111 * x158 * x163,
            x113 * x153 * x163,
            x117 * x155 * x163,
            x175 * x36 * x78,
            x178 * x41 * x78,
            x175 * x41 * x90,
            x182 * x74 * x78,
            x178 * x74 * x90,
            x175 * x74 * x82,
            x184 * x2,
            x182 * x185 * x64,
            x178 * x70 * x82,
            x175 * x70 * x95,
            x109 * x166 * x36,
            x109 * x168 * x41,
            x111 * x166 * x41,
            x109 * x170 * x74,
            x111 * x168 * x74,
            x113 * x166 * x74,
            x107 * x174 * x185,
            x111 * x170 * x70,
            x113 * x168 * x70,
            x117 * x166 * x70,
            x129 * x155 * x36,
            x129 * x153 * x41,
            x132 * x155 * x41,
            x129 * x158 * x74,
            x132 * x153 * x74,
            x136 * x155 * x74,
            x129 * x162 * x70,
            x132 * x158 * x70,
            x136 * x153 * x70,
            x137 * x186 * x2,
            x187 * x23 * x78,
            x188 * x21 * x78,
            x187 * x21 * x90,
            x189 * x4,
            x188 * x190 * x64,
            x187 * x82 * x9,
            x141 * (x183 * x96 + x3 * (x127 + 2 * x172 + 2 * x173 + 3 * x179 + 3 * x181)),
            x189 * x64,
            x188 * x8 * x82,
            x187 * x8 * x95,
            x109 * x175 * x23,
            x109 * x178 * x21,
            x111 * x175 * x21,
            x107 * x182 * x190,
            x111 * x178 * x9,
            x113 * x175 * x9,
            x107 * x184,
            x111 * x182 * x8,
            x113 * x178 * x8,
            x117 * x175 * x8,
            x129 * x166 * x23,
            x129 * x168 * x21,
            x132 * x166 * x21,
            x129 * x170 * x9,
            x132 * x168 * x9,
            x136 * x166 * x9,
            x129 * x174 * x8,
            x132 * x170 * x8,
            x136 * x168 * x8,
            x137 * x166 * x8,
            x146 * x155 * x23,
            x146 * x153 * x21,
            x147 * x155 * x21,
            x146 * x158 * x9,
            x147 * x153 * x9,
            x148 * x186 * x4,
            x146 * x162 * x8,
            x147 * x158 * x8,
            x148 * x153 * x8,
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
            x192 * x45 * x96,
            x100 * x195 * x61,
            x194 * x61 * x98,
            x102 * x195 * x76,
            x100 * x194 * x76,
            x198 * x76 * x98,
            x106 * x163 * x195,
            x102 * x163 * x194,
            x100 * x163 * x198,
            x163 * x201 * x98,
            x203 * x45 * x65,
            x203 * x61 * x91,
            x205 * x61 * x65,
            x203 * x68 * x76,
            x205 * x76 * x91,
            x207 * x65 * x76,
            x163 * x203 * x86,
            x163 * x205 * x68,
            x163 * x207 * x91,
            x163 * x211 * x65,
            x119 * x195 * x36,
            x122 * x195 * x41,
            x119 * x194 * x41,
            x126 * x195 * x74,
            x122 * x194 * x74,
            x119 * x198 * x74,
            x127 * x185 * x191,
            x126 * x194 * x70,
            x122 * x198 * x70,
            x119 * x201 * x70,
            x203 * x36 * x98,
            x100 * x203 * x41,
            x205 * x41 * x98,
            x102 * x203 * x74,
            x100 * x205 * x74,
            x207 * x74 * x98,
            x106 * x203 * x70,
            x102 * x205 * x70,
            x100 * x207 * x70,
            x211 * x212 * x96,
            x213 * x36 * x65,
            x213 * x41 * x91,
            x216 * x41 * x65,
            x213 * x68 * x74,
            x216 * x74 * x91,
            x220 * x65 * x74,
            x213 * x70 * x86,
            x216 * x68 * x70,
            x212 * x220 * x51,
            x2 * x222,
            x138 * x195 * x23,
            x139 * x195 * x21,
            x138 * x194 * x21,
            x140 * x190 * x191,
            x139 * x194 * x9,
            x138 * x198 * x9,
            x142 * x191,
            x140 * x194 * x8,
            x139 * x198 * x8,
            x138 * x201 * x8,
            x119 * x203 * x23,
            x122 * x203 * x21,
            x119 * x205 * x21,
            x126 * x203 * x9,
            x122 * x205 * x9,
            x119 * x207 * x9,
            x127 * x203 * x8,
            x126 * x205 * x8,
            x122 * x207 * x8,
            x119 * x211 * x8,
            x213 * x23 * x98,
            x100 * x21 * x213,
            x21 * x216 * x98,
            x102 * x213 * x9,
            x100 * x216 * x9,
            x220 * x223 * x96,
            x106 * x213 * x8,
            x102 * x216 * x8,
            x100 * x220 * x8,
            x222 * x96,
            x224 * x23 * x65,
            x21 * x224 * x91,
            x21 * x225 * x65,
            x224 * x68 * x9,
            x223 * x225 * x51,
            x226 * x4,
            x224 * x8 * x86,
            x225 * x68 * x8,
            x226 * x51,
            x144
            * (x107 * x221 + x3 * (x137 + 2 * x209 + 2 * x210 + 3 * x217 + 3 * x219)),
        ]
    )


def dipole3d_34(a, A, b, B, C):
    """Cartesian 3D (fg) dipole moment integrals.
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
    x12 = x5 ** 2 * x9
    x13 = x3 * x9
    x14 = 3 * x13
    x15 = x12 + x14
    x16 = x3 * (2 * x11 + x15)
    x17 = 4 * x16
    x18 = x4 * x9
    x19 = x3 * (x10 + x18)
    x20 = x11 + x13
    x21 = x20 * x5
    x22 = x19 + x21
    x23 = x22 * x5
    x24 = x3 * (3 * x12 + x14)
    x25 = x13 * x5
    x26 = 2 * x25
    x27 = x12 + x13
    x28 = x27 * x5
    x29 = x26 + x28
    x30 = x29 * x5
    x31 = x24 + x30
    x32 = x3 * (x17 + 4 * x23 + x31)
    x33 = 3 * x19
    x34 = x3 * (3 * x21 + x29 + x33)
    x35 = x16 + x23
    x36 = x35 * x5
    x37 = x34 + x36
    x38 = x2 * x37
    x39 = x32 + x38
    x40 = x2 * x35
    x41 = 8 * x25
    x42 = x3 * (4 * x28 + x41)
    x43 = x2 * x31
    x44 = x42 + x43
    x45 = x2 * x39 + x3 * (5 * x34 + x36 + 4 * x40 + x44)
    x46 = x2 * x22
    x47 = x2 * x29
    x48 = x24 + x47
    x49 = x3 * (x17 + x23 + 3 * x46 + x48)
    x50 = x34 + x40
    x51 = x2 * x50
    x52 = x2 * x44 + x3 * (5 * x24 + x30 + 4 * x47)
    x53 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x54 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x55 = numpy.pi * x0 * x54
    x56 = x53 * x55
    x57 = -x0 * (a * A[1] + b * B[1])
    x58 = -x57 - B[1]
    x59 = x49 + x51
    x60 = x2 * x20
    x61 = 2 * x60
    x62 = x2 * x27
    x63 = x26 + x62
    x64 = x3 * (x21 + x33 + x61 + x63)
    x65 = x16 + x46
    x66 = x2 * x65
    x67 = x3 * (x28 + x41 + 3 * x62)
    x68 = x2 * x48
    x69 = x67 + x68
    x70 = x56 * (x2 * x59 + x3 * (2 * x34 + 2 * x40 + 3 * x64 + 3 * x66 + x69))
    x71 = -x0 * (a * A[2] + b * B[2])
    x72 = -x71 - B[2]
    x73 = x53 * x8
    x74 = x3 * x73
    x75 = x58 ** 2 * x73
    x76 = x74 + x75
    x77 = x64 + x66
    x78 = x10 * x2
    x79 = x18 * x2
    x80 = x3 * (x11 + x14 + x78 + x79)
    x81 = x19 + x60
    x82 = x2 * x81
    x83 = 2 * x78
    x84 = x3 * (x15 + x83)
    x85 = x2 * x63
    x86 = x84 + x85
    x87 = x2 * x77 + x3 * (2 * x16 + 2 * x46 + 2 * x80 + 2 * x82 + x86)
    x88 = x54 * x8
    x89 = x56 * x72
    x90 = x3 * x88
    x91 = x72 ** 2 * x88
    x92 = x90 + x91
    x93 = x58 * x74
    x94 = 2 * x93
    x95 = x58 * x76
    x96 = x94 + x95
    x97 = x80 + x82
    x98 = x2 * x9
    x99 = x13 + x79
    x100 = x2 * x99 + x3 * (x18 + x98)
    x101 = x3 * (x10 + x98)
    x102 = x13 + x78
    x103 = x102 * x2
    x104 = x101 + x103
    x105 = x2 * x97 + x3 * (x100 + x104 + 2 * x19 + x61)
    x106 = x72 * x88
    x107 = x58 * x73
    x108 = x72 * x90
    x109 = 2 * x108
    x110 = x72 * x92
    x111 = x109 + x110
    x112 = 3 * x74
    x113 = x3 * (x112 + 3 * x75)
    x114 = x58 * x96
    x115 = x113 + x114
    x116 = x2 ** 2 * x9
    x117 = x116 + x14
    x118 = x100 * x2 + x3 * (x117 + 2 * x79)
    x119 = 3 * x90
    x120 = x3 * (x119 + 3 * x91)
    x121 = x111 * x72
    x122 = x120 + x121
    x123 = -x57 - A[1]
    x124 = x45 * x56
    x125 = x123 * x73
    x126 = x125 * x58
    x127 = x126 + x74
    x128 = x123 * x76
    x129 = x128 + x94
    x130 = x123 * x96
    x131 = x113 + x130
    x132 = 8 * x93
    x133 = x3 * (x132 + 4 * x95)
    x134 = x115 * x123
    x135 = x133 + x134
    x136 = -x71 - A[2]
    x137 = x136 * x56
    x138 = x136 * x88
    x139 = x138 * x72
    x140 = x139 + x90
    x141 = x136 * x92
    x142 = x109 + x141
    x143 = x111 * x136
    x144 = x120 + x143
    x145 = 8 * x108
    x146 = x3 * (4 * x110 + x145)
    x147 = x122 * x136
    x148 = x146 + x147
    x149 = x123 ** 2 * x73
    x150 = x149 + x74
    x151 = x3 * (x107 + x125)
    x152 = x123 * x127
    x153 = x151 + x152
    x154 = x112 + 2 * x126
    x155 = x3 * (x154 + x75)
    x156 = x123 * x129
    x157 = x155 + x156
    x158 = x3 * (3 * x128 + x132 + x95)
    x159 = x123 * x131
    x160 = x158 + x159
    x161 = x123 * x135 + x3 * (5 * x113 + x114 + 4 * x130)
    x162 = x136 ** 2 * x88
    x163 = x162 + x90
    x164 = x3 * (x106 + x138)
    x165 = x136 * x140
    x166 = x164 + x165
    x167 = x119 + 2 * x139
    x168 = x3 * (x167 + x91)
    x169 = x136 * x142
    x170 = x168 + x169
    x171 = x3 * (x110 + 3 * x141 + x145)
    x172 = x136 * x144
    x173 = x171 + x172
    x174 = x136 * x148 + x3 * (5 * x120 + x121 + 4 * x143)
    x175 = x123 * x150 + 2 * x123 * x74
    x176 = x123 * x153 + x3 * (x149 + x154)
    x177 = x123 * x157 + x3 * (2 * x128 + 2 * x151 + 2 * x152 + 4 * x93)
    x178 = x123 * x160 + x3 * (2 * x113 + 2 * x130 + 3 * x155 + 3 * x156)
    x179 = x55 * x7
    x180 = x179 * (x123 * x161 + x3 * (2 * x133 + 2 * x134 + 4 * x158 + 4 * x159))
    x181 = x179 * x4
    x182 = numpy.pi * x0 * x53 * x7
    x183 = x182 * x4
    x184 = x136 * x163 + 2 * x136 * x90
    x185 = x136 * x166 + x3 * (x162 + x167)
    x186 = x136 * x170 + x3 * (4 * x108 + 2 * x141 + 2 * x164 + 2 * x165)
    x187 = x136 * x173 + x3 * (2 * x120 + 2 * x143 + 3 * x168 + 3 * x169)
    x188 = x182 * (x136 * x174 + x3 * (2 * x146 + 2 * x147 + 4 * x171 + 4 * x172))
    x189 = -x57 - C[1]
    x190 = x56 * (x2 * x52 + x3 * (2 * x42 + 2 * x43 + 4 * x67 + 4 * x68))
    x191 = x107 * x189
    x192 = x191 + x74
    x193 = x2 * x69 + x3 * (2 * x24 + 2 * x47 + 3 * x84 + 3 * x85)
    x194 = x189 * x73
    x195 = x3 * (x107 + x194)
    x196 = x192 * x58
    x197 = x195 + x196
    x198 = x2 * x86 + x3 * (2 * x101 + 2 * x103 + 4 * x25 + 2 * x62)
    x199 = x3 * (x112 + 2 * x191 + x75)
    x200 = x197 * x58
    x201 = x199 + x200
    x202 = x104 * x2 + x3 * (x117 + x83)
    x203 = 3 * x195
    x204 = x3 * (3 * x196 + x203 + x96)
    x205 = x201 * x58
    x206 = x204 + x205
    x207 = x116 + x13
    x208 = 2 * x13 * x2 + x2 * x207
    x209 = x125 * x189
    x210 = x209 + x74
    x211 = x123 * x192
    x212 = x195 + x211
    x213 = x123 * x197
    x214 = x199 + x213
    x215 = x123 * x201
    x216 = x204 + x215
    x217 = 4 * x199
    x218 = x3 * (x115 + 4 * x200 + x217)
    x219 = x123 * x206
    x220 = x218 + x219
    x221 = x123 * x210 + x3 * (x125 + x194)
    x222 = x3 * (x112 + x126 + x191 + x209)
    x223 = x123 * x212
    x224 = x222 + x223
    x225 = x123 * x214
    x226 = 2 * x211
    x227 = x3 * (x129 + x196 + x203 + x226)
    x228 = x225 + x227
    x229 = x123 * x216
    x230 = x3 * (x131 + x200 + 3 * x213 + x217)
    x231 = x229 + x230
    x232 = x123 * x220 + x3 * (x135 + 5 * x204 + x205 + 4 * x215)
    x233 = x179 * x232
    x234 = x179 * x2
    x235 = x182 * x189
    x236 = x123 * x221 + x3 * (x112 + x149 + 2 * x209)
    x237 = x123 * x224 + x3 * (x153 + 2 * x195 + x221 + x226)
    x238 = x123 * x228 + x3 * (x157 + 2 * x199 + 2 * x213 + 2 * x222 + 2 * x223)
    x239 = x179 * (x123 * x231 + x3 * (x160 + 2 * x204 + 2 * x215 + 3 * x225 + 3 * x227))
    x240 = x179 * x5
    x241 = -x71 - C[2]
    x242 = x241 * x56
    x243 = x106 * x241
    x244 = x243 + x90
    x245 = x241 * x88
    x246 = x3 * (x106 + x245)
    x247 = x244 * x72
    x248 = x246 + x247
    x249 = x3 * (x119 + 2 * x243 + x91)
    x250 = x248 * x72
    x251 = x249 + x250
    x252 = 3 * x246
    x253 = x3 * (x111 + 3 * x247 + x252)
    x254 = x251 * x72
    x255 = x253 + x254
    x256 = x138 * x241
    x257 = x256 + x90
    x258 = x136 * x244
    x259 = x246 + x258
    x260 = x136 * x248
    x261 = x249 + x260
    x262 = x136 * x251
    x263 = x253 + x262
    x264 = 4 * x249
    x265 = x3 * (x122 + 4 * x250 + x264)
    x266 = x136 * x255
    x267 = x265 + x266
    x268 = x182 * x2
    x269 = x136 * x257 + x3 * (x138 + x245)
    x270 = x3 * (x119 + x139 + x243 + x256)
    x271 = x136 * x259
    x272 = x270 + x271
    x273 = x136 * x261
    x274 = 2 * x258
    x275 = x3 * (x142 + x247 + x252 + x274)
    x276 = x273 + x275
    x277 = x136 * x263
    x278 = x3 * (x144 + x250 + 3 * x260 + x264)
    x279 = x277 + x278
    x280 = x136 * x267 + x3 * (x148 + 5 * x253 + x254 + 4 * x262)
    x281 = x182 * x280
    x282 = x182 * x5
    x283 = x136 * x269 + x3 * (x119 + x162 + 2 * x256)
    x284 = x136 * x272 + x3 * (x166 + 2 * x246 + x269 + x274)
    x285 = x136 * x276 + x3 * (x170 + 2 * x249 + 2 * x260 + 2 * x270 + 2 * x271)
    x286 = x182 * (x136 * x279 + x3 * (x173 + 2 * x253 + 2 * x262 + 3 * x273 + 3 * x275))

    # 450 item(s)
    return numpy.array(
        [
            x56 * (x2 * x45 + x3 * (2 * x32 + 2 * x38 + 4 * x49 + 4 * x51 + x52)),
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
            x127 * x59 * x88,
            x123 * x59 * x89,
            x129 * x77 * x88,
            x106 * x127 * x77,
            x125 * x77 * x92,
            x131 * x88 * x97,
            x106 * x129 * x97,
            x127 * x92 * x97,
            x111 * x125 * x97,
            x100 * x135 * x88,
            x100 * x106 * x131,
            x100 * x129 * x92,
            x100 * x111 * x127,
            x100 * x122 * x125,
            x124 * x136,
            x137 * x58 * x59,
            x140 * x59 * x73,
            x138 * x76 * x77,
            x107 * x140 * x77,
            x142 * x73 * x77,
            x138 * x96 * x97,
            x140 * x76 * x97,
            x107 * x142 * x97,
            x144 * x73 * x97,
            x100 * x115 * x138,
            x100 * x140 * x96,
            x100 * x142 * x76,
            x100 * x107 * x144,
            x100 * x148 * x73,
            x150 * x39 * x88,
            x153 * x50 * x88,
            x106 * x150 * x50,
            x157 * x65 * x88,
            x106 * x153 * x65,
            x150 * x65 * x92,
            x160 * x81 * x88,
            x106 * x157 * x81,
            x153 * x81 * x92,
            x111 * x150 * x81,
            x161 * x88 * x99,
            x106 * x160 * x99,
            x157 * x92 * x99,
            x111 * x153 * x99,
            x122 * x150 * x99,
            x123 * x137 * x39,
            x127 * x138 * x50,
            x125 * x140 * x50,
            x129 * x138 * x65,
            x127 * x140 * x65,
            x125 * x142 * x65,
            x131 * x138 * x81,
            x129 * x140 * x81,
            x127 * x142 * x81,
            x125 * x144 * x81,
            x135 * x138 * x99,
            x131 * x140 * x99,
            x129 * x142 * x99,
            x127 * x144 * x99,
            x125 * x148 * x99,
            x163 * x39 * x73,
            x107 * x163 * x50,
            x166 * x50 * x73,
            x163 * x65 * x76,
            x107 * x166 * x65,
            x170 * x65 * x73,
            x163 * x81 * x96,
            x166 * x76 * x81,
            x107 * x170 * x81,
            x173 * x73 * x81,
            x115 * x163 * x99,
            x166 * x96 * x99,
            x170 * x76 * x99,
            x107 * x173 * x99,
            x174 * x73 * x99,
            x175 * x37 * x88,
            x176 * x35 * x88,
            x106 * x175 * x35,
            x177 * x22 * x88,
            x106 * x176 * x22,
            x175 * x22 * x92,
            x178 * x20 * x88,
            x106 * x177 * x20,
            x176 * x20 * x92,
            x111 * x175 * x20,
            x180 * x4,
            x178 * x181 * x72,
            x177 * x18 * x92,
            x111 * x176 * x18,
            x122 * x175 * x18,
            x138 * x150 * x37,
            x138 * x153 * x35,
            x140 * x150 * x35,
            x138 * x157 * x22,
            x140 * x153 * x22,
            x142 * x150 * x22,
            x138 * x160 * x20,
            x140 * x157 * x20,
            x142 * x153 * x20,
            x144 * x150 * x20,
            x136 * x161 * x181,
            x140 * x160 * x18,
            x142 * x157 * x18,
            x144 * x153 * x18,
            x148 * x150 * x18,
            x125 * x163 * x37,
            x127 * x163 * x35,
            x125 * x166 * x35,
            x129 * x163 * x22,
            x127 * x166 * x22,
            x125 * x170 * x22,
            x131 * x163 * x20,
            x129 * x166 * x20,
            x127 * x170 * x20,
            x125 * x173 * x20,
            x135 * x163 * x18,
            x131 * x166 * x18,
            x129 * x170 * x18,
            x127 * x173 * x18,
            x123 * x174 * x183,
            x184 * x37 * x73,
            x107 * x184 * x35,
            x185 * x35 * x73,
            x184 * x22 * x76,
            x107 * x185 * x22,
            x186 * x22 * x73,
            x184 * x20 * x96,
            x185 * x20 * x76,
            x107 * x186 * x20,
            x187 * x20 * x73,
            x115 * x18 * x184,
            x18 * x185 * x96,
            x18 * x186 * x76,
            x183 * x187 * x58,
            x188 * x4,
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
            x210 * x52 * x88,
            x212 * x69 * x88,
            x106 * x210 * x69,
            x214 * x86 * x88,
            x106 * x212 * x86,
            x210 * x86 * x92,
            x104 * x216 * x88,
            x104 * x106 * x214,
            x104 * x212 * x92,
            x104 * x111 * x210,
            x207 * x220 * x88,
            x106 * x207 * x216,
            x207 * x214 * x92,
            x111 * x207 * x212,
            x122 * x207 * x210,
            x137 * x189 * x52,
            x138 * x192 * x69,
            x140 * x194 * x69,
            x138 * x197 * x86,
            x140 * x192 * x86,
            x142 * x194 * x86,
            x104 * x138 * x201,
            x104 * x140 * x197,
            x104 * x142 * x192,
            x104 * x144 * x194,
            x138 * x206 * x207,
            x140 * x201 * x207,
            x142 * x197 * x207,
            x144 * x192 * x207,
            x148 * x194 * x207,
            x221 * x44 * x88,
            x224 * x48 * x88,
            x106 * x221 * x48,
            x228 * x63 * x88,
            x106 * x224 * x63,
            x221 * x63 * x92,
            x102 * x231 * x88,
            x102 * x106 * x228,
            x102 * x224 * x92,
            x102 * x111 * x221,
            x2 * x233,
            x231 * x234 * x72,
            x228 * x92 * x98,
            x111 * x224 * x98,
            x122 * x221 * x98,
            x138 * x210 * x44,
            x138 * x212 * x48,
            x140 * x210 * x48,
            x138 * x214 * x63,
            x140 * x212 * x63,
            x142 * x210 * x63,
            x102 * x138 * x216,
            x102 * x140 * x214,
            x102 * x142 * x212,
            x102 * x144 * x210,
            x136 * x220 * x234,
            x140 * x216 * x98,
            x142 * x214 * x98,
            x144 * x212 * x98,
            x148 * x210 * x98,
            x163 * x194 * x44,
            x163 * x192 * x48,
            x166 * x194 * x48,
            x163 * x197 * x63,
            x166 * x192 * x63,
            x170 * x194 * x63,
            x102 * x163 * x201,
            x102 * x166 * x197,
            x102 * x170 * x192,
            x102 * x173 * x194,
            x163 * x206 * x98,
            x166 * x201 * x98,
            x170 * x197 * x98,
            x173 * x192 * x98,
            x174 * x2 * x235,
            x236 * x31 * x88,
            x237 * x29 * x88,
            x106 * x236 * x29,
            x238 * x27 * x88,
            x106 * x237 * x27,
            x236 * x27 * x92,
            x239 * x5,
            x238 * x240 * x72,
            x10 * x237 * x92,
            x10 * x111 * x236,
            x179
            * (x123 * x232 + x3 * (x161 + 2 * x218 + 2 * x219 + 4 * x229 + 4 * x230)),
            x239 * x72,
            x238 * x9 * x92,
            x111 * x237 * x9,
            x122 * x236 * x9,
            x138 * x221 * x31,
            x138 * x224 * x29,
            x140 * x221 * x29,
            x138 * x228 * x27,
            x140 * x224 * x27,
            x142 * x221 * x27,
            x136 * x231 * x240,
            x10 * x140 * x228,
            x10 * x142 * x224,
            x10 * x144 * x221,
            x136 * x233,
            x140 * x231 * x9,
            x142 * x228 * x9,
            x144 * x224 * x9,
            x148 * x221 * x9,
            x163 * x210 * x31,
            x163 * x212 * x29,
            x166 * x210 * x29,
            x163 * x214 * x27,
            x166 * x212 * x27,
            x170 * x210 * x27,
            x10 * x163 * x216,
            x10 * x166 * x214,
            x10 * x170 * x212,
            x10 * x173 * x210,
            x163 * x220 * x9,
            x166 * x216 * x9,
            x170 * x214 * x9,
            x173 * x212 * x9,
            x174 * x210 * x9,
            x184 * x194 * x31,
            x184 * x192 * x29,
            x185 * x194 * x29,
            x184 * x197 * x27,
            x185 * x192 * x27,
            x186 * x194 * x27,
            x10 * x184 * x201,
            x10 * x185 * x197,
            x10 * x186 * x192,
            x187 * x235 * x5,
            x184 * x206 * x9,
            x185 * x201 * x9,
            x186 * x197 * x9,
            x187 * x192 * x9,
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
            x123 * x242 * x52,
            x127 * x245 * x69,
            x125 * x244 * x69,
            x129 * x245 * x86,
            x127 * x244 * x86,
            x125 * x248 * x86,
            x104 * x131 * x245,
            x104 * x129 * x244,
            x104 * x127 * x248,
            x104 * x125 * x251,
            x135 * x207 * x245,
            x131 * x207 * x244,
            x129 * x207 * x248,
            x127 * x207 * x251,
            x125 * x207 * x255,
            x257 * x52 * x73,
            x107 * x257 * x69,
            x259 * x69 * x73,
            x257 * x76 * x86,
            x107 * x259 * x86,
            x261 * x73 * x86,
            x104 * x257 * x96,
            x104 * x259 * x76,
            x104 * x107 * x261,
            x104 * x263 * x73,
            x115 * x207 * x257,
            x207 * x259 * x96,
            x207 * x261 * x76,
            x107 * x207 * x263,
            x207 * x267 * x73,
            x150 * x245 * x44,
            x153 * x245 * x48,
            x150 * x244 * x48,
            x157 * x245 * x63,
            x153 * x244 * x63,
            x150 * x248 * x63,
            x102 * x160 * x245,
            x102 * x157 * x244,
            x102 * x153 * x248,
            x102 * x150 * x251,
            x161 * x234 * x241,
            x160 * x244 * x98,
            x157 * x248 * x98,
            x153 * x251 * x98,
            x150 * x255 * x98,
            x125 * x257 * x44,
            x127 * x257 * x48,
            x125 * x259 * x48,
            x129 * x257 * x63,
            x127 * x259 * x63,
            x125 * x261 * x63,
            x102 * x131 * x257,
            x102 * x129 * x259,
            x102 * x127 * x261,
            x102 * x125 * x263,
            x135 * x257 * x98,
            x131 * x259 * x98,
            x129 * x261 * x98,
            x127 * x263 * x98,
            x123 * x267 * x268,
            x269 * x44 * x73,
            x107 * x269 * x48,
            x272 * x48 * x73,
            x269 * x63 * x76,
            x107 * x272 * x63,
            x276 * x63 * x73,
            x102 * x269 * x96,
            x102 * x272 * x76,
            x102 * x107 * x276,
            x102 * x279 * x73,
            x115 * x269 * x98,
            x272 * x96 * x98,
            x276 * x76 * x98,
            x268 * x279 * x58,
            x2 * x281,
            x175 * x245 * x31,
            x176 * x245 * x29,
            x175 * x244 * x29,
            x177 * x245 * x27,
            x176 * x244 * x27,
            x175 * x248 * x27,
            x178 * x240 * x241,
            x10 * x177 * x244,
            x10 * x176 * x248,
            x10 * x175 * x251,
            x180 * x241,
            x178 * x244 * x9,
            x177 * x248 * x9,
            x176 * x251 * x9,
            x175 * x255 * x9,
            x150 * x257 * x31,
            x153 * x257 * x29,
            x150 * x259 * x29,
            x157 * x257 * x27,
            x153 * x259 * x27,
            x150 * x261 * x27,
            x10 * x160 * x257,
            x10 * x157 * x259,
            x10 * x153 * x261,
            x10 * x150 * x263,
            x161 * x257 * x9,
            x160 * x259 * x9,
            x157 * x261 * x9,
            x153 * x263 * x9,
            x150 * x267 * x9,
            x125 * x269 * x31,
            x127 * x269 * x29,
            x125 * x272 * x29,
            x129 * x269 * x27,
            x127 * x27 * x272,
            x125 * x27 * x276,
            x10 * x131 * x269,
            x10 * x129 * x272,
            x10 * x127 * x276,
            x123 * x279 * x282,
            x135 * x269 * x9,
            x131 * x272 * x9,
            x129 * x276 * x9,
            x127 * x279 * x9,
            x123 * x281,
            x283 * x31 * x73,
            x107 * x283 * x29,
            x284 * x29 * x73,
            x27 * x283 * x76,
            x107 * x27 * x284,
            x27 * x285 * x73,
            x10 * x283 * x96,
            x10 * x284 * x76,
            x282 * x285 * x58,
            x286 * x5,
            x115 * x283 * x9,
            x284 * x9 * x96,
            x285 * x76 * x9,
            x286 * x58,
            x182
            * (x136 * x280 + x3 * (x174 + 2 * x265 + 2 * x266 + 4 * x277 + 4 * x278)),
        ]
    )


def dipole3d_40(a, A, b, B, C):
    """Cartesian 3D (gs) dipole moment integrals.
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
    x8 = x3 * x7
    x9 = -x2 - C[0]
    x10 = x7 * x9
    x11 = x0 * (x10 + x8)
    x12 = x0 * x7
    x13 = x8 * x9
    x14 = x12 + x13
    x15 = x14 * x3
    x16 = x3 ** 2 * x7
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
    x27 = -x26 - A[1]
    x28 = x21 * x25
    x29 = -x1 * (a * A[2] + b * B[2])
    x30 = -x29 - A[2]
    x31 = x22 * x6
    x32 = x0 * x31
    x33 = x27 ** 2 * x31
    x34 = x32 + x33
    x35 = x23 * x6
    x36 = x25 * x30
    x37 = x0 * x35
    x38 = x30 ** 2 * x35
    x39 = x37 + x38
    x40 = 2 * x27 * x32 + x27 * x34
    x41 = x30 * x35
    x42 = x27 * x31
    x43 = 2 * x30 * x37 + x30 * x39
    x44 = 3 * x32
    x45 = x24 * x5
    x46 = x45 * (x0 * (3 * x33 + x44) + x27 * x40)
    x47 = x30 * x45
    x48 = numpy.pi * x1 * x22 * x5
    x49 = x43 * x48
    x50 = 3 * x37
    x51 = x48 * (x0 * (3 * x38 + x50) + x30 * x43)
    x52 = -x26 - C[1]
    x53 = x25 * (x0 * (3 * x16 + x19) + x18 * x3)
    x54 = x42 * x52
    x55 = x32 + x54
    x56 = x31 * x52
    x57 = x0 * (x42 + x56)
    x58 = x27 * x55
    x59 = x57 + x58
    x60 = x0 * (x33 + x44 + 2 * x54) + x27 * x59
    x61 = x45 * x60
    x62 = -x29 - C[2]
    x63 = x41 * x62
    x64 = x37 + x63
    x65 = x35 * x62
    x66 = x0 * (x41 + x65)
    x67 = x30 * x64
    x68 = x66 + x67
    x69 = x0 * (x38 + x50 + 2 * x63) + x30 * x68
    x70 = x48 * x69

    # 45 item(s)
    return numpy.array(
        [
            x25 * (x0 * (3 * x11 + 3 * x15 + x18) + x21 * x3),
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
            x45 * (x0 * (x40 + 3 * x57 + 3 * x58) + x27 * x60),
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
            x48 * (x0 * (x43 + 3 * x66 + 3 * x67) + x30 * x69),
        ]
    )


def dipole3d_41(a, A, b, B, C):
    """Cartesian 3D (gp) dipole moment integrals.
    The origin is at C.

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
    x11 = x2 * x7
    x12 = x10 * x11
    x13 = -x1 - C[0]
    x14 = x11 * x13
    x15 = x10 * x7
    x16 = x13 * x15
    x17 = x3 * (x12 + x14 + x16 + x9)
    x18 = x13 * x7
    x19 = x3 * (x15 + x18)
    x20 = x16 + x8
    x21 = x2 * x20
    x22 = x19 + x21
    x23 = x2 * x22
    x24 = x17 + x23
    x25 = x3 * (x11 + x15)
    x26 = x12 + x8
    x27 = x2 * x26
    x28 = x25 + x27
    x29 = x3 * (x11 + x18)
    x30 = x14 + x8
    x31 = x2 * x30
    x32 = x29 + x31
    x33 = x2 * x24 + x3 * (2 * x19 + 2 * x21 + x28 + x32)
    x34 = x2 ** 2 * x7
    x35 = x34 + x9
    x36 = x2 * x32 + x3 * (2 * x14 + x35)
    x37 = x2 * x28 + x3 * (2 * x12 + x35)
    x38 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x39 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x40 = numpy.pi * x0 * x39
    x41 = x38 * x40
    x42 = -x0 * (a * A[1] + b * B[1])
    x43 = -x42 - B[1]
    x44 = x34 + x8
    x45 = 2 * x11 * x3 + x2 * x44
    x46 = x41 * (x2 * x36 + x3 * (3 * x29 + 3 * x31 + x45))
    x47 = -x0 * (a * A[2] + b * B[2])
    x48 = -x47 - B[2]
    x49 = -x42 - A[1]
    x50 = x33 * x41
    x51 = x3 * x6
    x52 = x38 * x51
    x53 = x38 * x6
    x54 = x49 * x53
    x55 = x43 * x54
    x56 = x52 + x55
    x57 = x39 * x6
    x58 = x36 * x41
    x59 = -x47 - A[2]
    x60 = x39 * x51
    x61 = x57 * x59
    x62 = x48 * x61
    x63 = x60 + x62
    x64 = x49 ** 2 * x53
    x65 = x52 + x64
    x66 = x43 * x53
    x67 = x3 * (x54 + x66)
    x68 = x49 * x56
    x69 = x67 + x68
    x70 = x48 * x57
    x71 = x41 * x59
    x72 = x57 * x59 ** 2
    x73 = x60 + x72
    x74 = x3 * (x61 + x70)
    x75 = x59 * x63
    x76 = x74 + x75
    x77 = 2 * x49 * x52 + x49 * x65
    x78 = 3 * x52
    x79 = x64 + x78
    x80 = x3 * (2 * x55 + x79) + x49 * x69
    x81 = 2 * x59 * x60 + x59 * x73
    x82 = 3 * x60
    x83 = x72 + x82
    x84 = x3 * (2 * x62 + x83) + x59 * x76
    x85 = x3 * (3 * x64 + x78) + x49 * x77
    x86 = x40 * x5
    x87 = x86 * (x3 * (3 * x67 + 3 * x68 + x77) + x49 * x80)
    x88 = x13 * x86
    x89 = numpy.pi * x0 * x38 * x5
    x90 = x13 * x89
    x91 = x3 * (3 * x72 + x82) + x59 * x81
    x92 = x89 * (x3 * (3 * x74 + 3 * x75 + x81) + x59 * x84)
    x93 = -x42 - C[1]
    x94 = x41 * (x2 * x37 + x3 * (3 * x25 + 3 * x27 + x45))
    x95 = x66 * x93
    x96 = x52 + x95
    x97 = x2 * x45 + x3 * (3 * x34 + x9)
    x98 = x41 * x97
    x99 = x54 * x93
    x100 = x52 + x99
    x101 = x53 * x93
    x102 = x3 * (x101 + x66)
    x103 = x49 * x96
    x104 = x102 + x103
    x105 = x3 * (x101 + x54)
    x106 = x100 * x49
    x107 = x105 + x106
    x108 = x3 * (x55 + x78 + x95 + x99)
    x109 = x104 * x49
    x110 = x108 + x109
    x111 = x107 * x49 + x3 * (x79 + 2 * x99)
    x112 = x110 * x49 + x3 * (2 * x102 + 2 * x103 + x107 + x69)
    x113 = x112 * x86
    x114 = x2 * x86
    x115 = x89 * x93
    x116 = x86 * (x111 * x49 + x3 * (3 * x105 + 3 * x106 + x77))
    x117 = x10 * x86
    x118 = -x47 - C[2]
    x119 = x118 * x70
    x120 = x119 + x60
    x121 = x118 * x57
    x122 = x118 * x61
    x123 = x122 + x60
    x124 = x3 * (x121 + x70)
    x125 = x120 * x59
    x126 = x124 + x125
    x127 = x3 * (x121 + x61)
    x128 = x123 * x59
    x129 = x127 + x128
    x130 = x3 * (x119 + x122 + x62 + x82)
    x131 = x126 * x59
    x132 = x130 + x131
    x133 = x2 * x89
    x134 = x129 * x59 + x3 * (2 * x122 + x83)
    x135 = x132 * x59 + x3 * (2 * x124 + 2 * x125 + x129 + x76)
    x136 = x135 * x89
    x137 = x89 * (x134 * x59 + x3 * (3 * x127 + 3 * x128 + x81))

    # 135 item(s)
    return numpy.array(
        [
            x41 * (x2 * x33 + x3 * (3 * x17 + 3 * x23 + x36 + x37)),
            x43 * x46,
            x46 * x48,
            x49 * x50,
            x36 * x56 * x57,
            x48 * x49 * x58,
            x50 * x59,
            x43 * x58 * x59,
            x36 * x53 * x63,
            x24 * x57 * x65,
            x32 * x57 * x69,
            x32 * x65 * x70,
            x24 * x49 * x71,
            x32 * x56 * x61,
            x32 * x54 * x63,
            x24 * x53 * x73,
            x32 * x66 * x73,
            x32 * x53 * x76,
            x22 * x57 * x77,
            x30 * x57 * x80,
            x30 * x70 * x77,
            x22 * x61 * x65,
            x30 * x61 * x69,
            x30 * x63 * x65,
            x22 * x54 * x73,
            x30 * x56 * x73,
            x30 * x54 * x76,
            x22 * x53 * x81,
            x30 * x66 * x81,
            x30 * x53 * x84,
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
            x100 * x37 * x57,
            x104 * x45 * x57,
            x100 * x45 * x70,
            x37 * x71 * x93,
            x45 * x61 * x96,
            x101 * x45 * x63,
            x107 * x28 * x57,
            x110 * x44 * x57,
            x107 * x44 * x70,
            x100 * x28 * x61,
            x104 * x44 * x61,
            x100 * x44 * x63,
            x101 * x28 * x73,
            x44 * x73 * x96,
            x101 * x44 * x76,
            x111 * x26 * x57,
            x113 * x2,
            x111 * x114 * x48,
            x107 * x26 * x61,
            x110 * x114 * x59,
            x107 * x11 * x63,
            x100 * x26 * x73,
            x104 * x11 * x73,
            x100 * x11 * x76,
            x101 * x26 * x81,
            x11 * x81 * x96,
            x115 * x2 * x84,
            x10 * x116,
            x86 * (x112 * x49 + x3 * (3 * x108 + 3 * x109 + x111 + x80)),
            x116 * x48,
            x111 * x117 * x59,
            x113 * x59,
            x111 * x63 * x7,
            x107 * x15 * x73,
            x110 * x7 * x73,
            x107 * x7 * x76,
            x100 * x15 * x81,
            x104 * x7 * x81,
            x100 * x7 * x84,
            x10 * x115 * x91,
            x7 * x91 * x96,
            x92 * x93,
            x118 * x94,
            x118 * x43 * x98,
            x120 * x53 * x97,
            x118 * x37 * x41 * x49,
            x121 * x45 * x56,
            x120 * x45 * x54,
            x123 * x37 * x53,
            x123 * x45 * x66,
            x126 * x45 * x53,
            x121 * x28 * x65,
            x121 * x44 * x69,
            x120 * x44 * x65,
            x123 * x28 * x54,
            x123 * x44 * x56,
            x126 * x44 * x54,
            x129 * x28 * x53,
            x129 * x44 * x66,
            x132 * x44 * x53,
            x121 * x26 * x77,
            x114 * x118 * x80,
            x11 * x120 * x77,
            x123 * x26 * x65,
            x11 * x123 * x69,
            x11 * x126 * x65,
            x129 * x26 * x54,
            x11 * x129 * x56,
            x132 * x133 * x49,
            x134 * x26 * x53,
            x133 * x134 * x43,
            x136 * x2,
            x117 * x118 * x85,
            x118 * x87,
            x120 * x7 * x85,
            x123 * x15 * x77,
            x123 * x7 * x80,
            x126 * x7 * x77,
            x129 * x15 * x65,
            x129 * x69 * x7,
            x132 * x65 * x7,
            x10 * x134 * x49 * x89,
            x134 * x56 * x7,
            x136 * x49,
            x10 * x137,
            x137 * x43,
            x89 * (x135 * x59 + x3 * (3 * x130 + 3 * x131 + x134 + x84)),
        ]
    )


def dipole3d_42(a, A, b, B, C):
    """Cartesian 3D (gd) dipole moment integrals.
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
    x12 = x5 ** 2 * x9
    x13 = x3 * x9
    x14 = 3 * x13
    x15 = x12 + x14
    x16 = x3 * (2 * x11 + x15)
    x17 = x4 * x9
    x18 = x3 * (x10 + x17)
    x19 = x11 + x13
    x20 = x19 * x5
    x21 = x18 + x20
    x22 = x2 * x21
    x23 = x16 + x22
    x24 = x2 * x23
    x25 = x19 * x2
    x26 = 2 * x25
    x27 = x10 * x3
    x28 = x12 + x13
    x29 = x2 * x28
    x30 = 2 * x27 + x29
    x31 = x3 * (3 * x18 + x20 + x26 + x30)
    x32 = x24 + x31
    x33 = x2 * x7
    x34 = x33 * x8
    x35 = x34 * x5
    x36 = x34 * x4
    x37 = x3 * (x11 + x14 + x35 + x36)
    x38 = x18 + x25
    x39 = x2 * x38
    x40 = 2 * x35
    x41 = x3 * (x15 + x40)
    x42 = x2 * x30
    x43 = x41 + x42
    x44 = x2 * x32 + x3 * (2 * x16 + 2 * x22 + 2 * x37 + 2 * x39 + x43)
    x45 = x3 * (x10 + x34)
    x46 = x13 + x35
    x47 = x2 * x46
    x48 = x45 + x47
    x49 = x3 * (x17 + x34)
    x50 = x13 + x36
    x51 = x2 * x50
    x52 = x49 + x51
    x53 = x3 * (2 * x18 + x26 + x48 + x52)
    x54 = x37 + x39
    x55 = x2 * x54
    x56 = x2 * x43 + x3 * (4 * x27 + 2 * x29 + 2 * x45 + 2 * x47)
    x57 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x58 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x59 = numpy.pi * x0 * x58
    x60 = x57 * x59
    x61 = -x0 * (a * A[1] + b * B[1])
    x62 = -x61 - B[1]
    x63 = x53 + x55
    x64 = x2 ** 2 * x9
    x65 = x14 + x64
    x66 = x2 * x52 + x3 * (2 * x36 + x65)
    x67 = x3 * (x40 + x65)
    x68 = x2 * x48
    x69 = x67 + x68
    x70 = x60 * (x2 * x63 + x3 * (3 * x37 + 3 * x39 + x66 + x69))
    x71 = -x0 * (a * A[2] + b * B[2])
    x72 = -x71 - B[2]
    x73 = x57 * x8
    x74 = x3 * x73
    x75 = x62 ** 2 * x73
    x76 = x74 + x75
    x77 = x13 + x64
    x78 = x2 * x77 + 2 * x3 * x34
    x79 = x2 * x66 + x3 * (3 * x49 + 3 * x51 + x78)
    x80 = x58 * x8
    x81 = x60 * x72
    x82 = x3 * x80
    x83 = x72 ** 2 * x80
    x84 = x82 + x83
    x85 = -x61 - A[1]
    x86 = x44 * x60
    x87 = x73 * x85
    x88 = x62 * x87
    x89 = x74 + x88
    x90 = 2 * x74
    x91 = x76 * x85
    x92 = x62 * x90 + x91
    x93 = x72 * x80
    x94 = -x71 - A[2]
    x95 = x60 * x94
    x96 = x80 * x94
    x97 = x72 * x96
    x98 = x82 + x97
    x99 = x62 * x73
    x100 = 2 * x82
    x101 = x84 * x94
    x102 = x100 * x72 + x101
    x103 = x73 * x85 ** 2
    x104 = x103 + x74
    x105 = x3 * (x87 + x99)
    x106 = x85 * x89
    x107 = x105 + x106
    x108 = 3 * x74
    x109 = x108 + 2 * x88
    x110 = x3 * (x109 + x75)
    x111 = x85 * x92
    x112 = x110 + x111
    x113 = x80 * x94 ** 2
    x114 = x113 + x82
    x115 = x3 * (x93 + x96)
    x116 = x94 * x98
    x117 = x115 + x116
    x118 = 3 * x82
    x119 = x118 + 2 * x97
    x120 = x3 * (x119 + x83)
    x121 = x102 * x94
    x122 = x120 + x121
    x123 = x104 * x85 + x85 * x90
    x124 = x3 * (x103 + x109)
    x125 = x107 * x85
    x126 = x124 + x125
    x127 = x112 * x85 + x3 * (2 * x105 + 2 * x106 + 4 * x62 * x74 + 2 * x91)
    x128 = x100 * x94 + x114 * x94
    x129 = x3 * (x113 + x119)
    x130 = x117 * x94
    x131 = x129 + x130
    x132 = x122 * x94 + x3 * (2 * x101 + 2 * x115 + 2 * x116 + 4 * x72 * x82)
    x133 = x123 * x85 + x3 * (3 * x103 + x108)
    x134 = x126 * x85 + x3 * (3 * x105 + 3 * x106 + x123)
    x135 = x59 * x7
    x136 = x135 * (x127 * x85 + x3 * (3 * x110 + 3 * x111 + 2 * x124 + 2 * x125))
    x137 = x135 * x72
    x138 = x135 * x94
    x139 = numpy.pi * x0 * x57
    x140 = x139 * x7
    x141 = x140 * x85
    x142 = x128 * x94 + x3 * (3 * x113 + x118)
    x143 = x131 * x94 + x3 * (3 * x115 + 3 * x116 + x128)
    x144 = x140 * x143
    x145 = x140 * (x132 * x94 + x3 * (3 * x120 + 3 * x121 + 2 * x129 + 2 * x130))
    x146 = -x61 - C[1]
    x147 = x60 * (x2 * x56 + x3 * (3 * x41 + 3 * x42 + 2 * x67 + 2 * x68))
    x148 = x146 * x99
    x149 = x148 + x74
    x150 = x2 * x69 + x3 * (3 * x45 + 3 * x47 + x78)
    x151 = x146 * x73
    x152 = x3 * (x151 + x99)
    x153 = x149 * x62
    x154 = x152 + x153
    x155 = x2 * x78 + x3 * (x14 + 3 * x64)
    x156 = x146 * x87
    x157 = x156 + x74
    x158 = x149 * x85
    x159 = x152 + x158
    x160 = x3 * (x108 + 2 * x148 + x75)
    x161 = x154 * x85
    x162 = x160 + x161
    x163 = x3 * (x151 + x87)
    x164 = x157 * x85
    x165 = x163 + x164
    x166 = x3 * (x108 + x148 + x156 + x88)
    x167 = x159 * x85
    x168 = x166 + x167
    x169 = x162 * x85
    x170 = 2 * x158
    x171 = x3 * (3 * x152 + x153 + x170 + x92)
    x172 = x169 + x171
    x173 = x165 * x85 + x3 * (x103 + x108 + 2 * x156)
    x174 = x168 * x85
    x175 = x3 * (x107 + 2 * x152 + x165 + x170)
    x176 = x174 + x175
    x177 = x172 * x85 + x3 * (x112 + 2 * x160 + 2 * x161 + 2 * x166 + 2 * x167)
    x178 = x33 * x59
    x179 = x139 * x33
    x180 = x173 * x85 + x3 * (x123 + 3 * x163 + 3 * x164)
    x181 = x135 * (x176 * x85 + x3 * (x126 + 3 * x166 + 3 * x167 + x173))
    x182 = -x71 - C[2]
    x183 = x182 * x60
    x184 = x182 * x93
    x185 = x184 + x82
    x186 = x182 * x80
    x187 = x3 * (x186 + x93)
    x188 = x185 * x72
    x189 = x187 + x188
    x190 = x182 * x96
    x191 = x190 + x82
    x192 = x185 * x94
    x193 = x187 + x192
    x194 = x3 * (x118 + 2 * x184 + x83)
    x195 = x189 * x94
    x196 = x194 + x195
    x197 = x3 * (x186 + x96)
    x198 = x191 * x94
    x199 = x197 + x198
    x200 = x3 * (x118 + x184 + x190 + x97)
    x201 = x193 * x94
    x202 = x200 + x201
    x203 = x196 * x94
    x204 = 2 * x192
    x205 = x3 * (x102 + 3 * x187 + x188 + x204)
    x206 = x203 + x205
    x207 = x199 * x94 + x3 * (x113 + x118 + 2 * x190)
    x208 = x202 * x94
    x209 = x3 * (x117 + 2 * x187 + x199 + x204)
    x210 = x208 + x209
    x211 = x206 * x94 + x3 * (x122 + 2 * x194 + 2 * x195 + 2 * x200 + 2 * x201)
    x212 = x207 * x94 + x3 * (x128 + 3 * x197 + 3 * x198)
    x213 = x140 * (x210 * x94 + x3 * (x131 + 3 * x200 + 3 * x201 + x207))

    # 270 item(s)
    return numpy.array(
        [
            x60 * (x2 * x44 + x3 * (3 * x24 + 3 * x31 + 2 * x53 + 2 * x55 + x56)),
            x62 * x70,
            x70 * x72,
            x76 * x79 * x80,
            x62 * x79 * x81,
            x73 * x79 * x84,
            x85 * x86,
            x63 * x80 * x89,
            x63 * x81 * x85,
            x66 * x80 * x92,
            x66 * x89 * x93,
            x66 * x84 * x87,
            x86 * x94,
            x62 * x63 * x95,
            x63 * x73 * x98,
            x66 * x76 * x96,
            x66 * x98 * x99,
            x102 * x66 * x73,
            x104 * x32 * x80,
            x107 * x54 * x80,
            x104 * x54 * x93,
            x112 * x52 * x80,
            x107 * x52 * x93,
            x104 * x52 * x84,
            x32 * x85 * x95,
            x54 * x89 * x96,
            x54 * x87 * x98,
            x52 * x92 * x96,
            x52 * x89 * x98,
            x102 * x52 * x87,
            x114 * x32 * x73,
            x114 * x54 * x99,
            x117 * x54 * x73,
            x114 * x52 * x76,
            x117 * x52 * x99,
            x122 * x52 * x73,
            x123 * x23 * x80,
            x126 * x38 * x80,
            x123 * x38 * x93,
            x127 * x50 * x80,
            x126 * x50 * x93,
            x123 * x50 * x84,
            x104 * x23 * x96,
            x107 * x38 * x96,
            x104 * x38 * x98,
            x112 * x50 * x96,
            x107 * x50 * x98,
            x102 * x104 * x50,
            x114 * x23 * x87,
            x114 * x38 * x89,
            x117 * x38 * x87,
            x114 * x50 * x92,
            x117 * x50 * x89,
            x122 * x50 * x87,
            x128 * x23 * x73,
            x128 * x38 * x99,
            x131 * x38 * x73,
            x128 * x50 * x76,
            x131 * x50 * x99,
            x132 * x50 * x73,
            x133 * x21 * x80,
            x134 * x19 * x80,
            x133 * x19 * x93,
            x136 * x4,
            x134 * x137 * x4,
            x133 * x17 * x84,
            x123 * x21 * x96,
            x126 * x19 * x96,
            x123 * x19 * x98,
            x127 * x138 * x4,
            x126 * x17 * x98,
            x102 * x123 * x17,
            x104 * x114 * x21,
            x107 * x114 * x19,
            x104 * x117 * x19,
            x112 * x114 * x17,
            x107 * x117 * x17,
            x104 * x122 * x17,
            x128 * x21 * x87,
            x128 * x19 * x89,
            x131 * x19 * x87,
            x128 * x17 * x92,
            x131 * x17 * x89,
            x132 * x141 * x4,
            x142 * x21 * x73,
            x142 * x19 * x99,
            x143 * x19 * x73,
            x142 * x17 * x76,
            x144 * x4 * x62,
            x145 * x4,
            x146 * x147,
            x149 * x150 * x80,
            x146 * x150 * x81,
            x154 * x155 * x80,
            x149 * x155 * x93,
            x151 * x155 * x84,
            x157 * x56 * x80,
            x159 * x69 * x80,
            x157 * x69 * x93,
            x162 * x78 * x80,
            x159 * x78 * x93,
            x157 * x78 * x84,
            x146 * x56 * x95,
            x149 * x69 * x96,
            x151 * x69 * x98,
            x154 * x78 * x96,
            x149 * x78 * x98,
            x102 * x151 * x78,
            x165 * x43 * x80,
            x168 * x48 * x80,
            x165 * x48 * x93,
            x172 * x77 * x80,
            x168 * x77 * x93,
            x165 * x77 * x84,
            x157 * x43 * x96,
            x159 * x48 * x96,
            x157 * x48 * x98,
            x162 * x77 * x96,
            x159 * x77 * x98,
            x102 * x157 * x77,
            x114 * x151 * x43,
            x114 * x149 * x48,
            x117 * x151 * x48,
            x114 * x154 * x77,
            x117 * x149 * x77,
            x122 * x151 * x77,
            x173 * x30 * x80,
            x176 * x46 * x80,
            x173 * x46 * x93,
            x177 * x178,
            x176 * x178 * x72,
            x173 * x34 * x84,
            x165 * x30 * x96,
            x168 * x46 * x96,
            x165 * x46 * x98,
            x172 * x178 * x94,
            x168 * x34 * x98,
            x102 * x165 * x34,
            x114 * x157 * x30,
            x114 * x159 * x46,
            x117 * x157 * x46,
            x114 * x162 * x34,
            x117 * x159 * x34,
            x122 * x157 * x34,
            x128 * x151 * x30,
            x128 * x149 * x46,
            x131 * x151 * x46,
            x128 * x154 * x34,
            x131 * x149 * x34,
            x132 * x146 * x179,
            x180 * x28 * x80,
            x181 * x5,
            x137 * x180 * x5,
            x135 * (x177 * x85 + x3 * (x127 + 3 * x169 + 3 * x171 + 2 * x174 + 2 * x175)),
            x181 * x72,
            x180 * x84 * x9,
            x173 * x28 * x96,
            x138 * x176 * x5,
            x10 * x173 * x98,
            x138 * x177,
            x176 * x9 * x98,
            x102 * x173 * x9,
            x114 * x165 * x28,
            x10 * x114 * x168,
            x10 * x117 * x165,
            x114 * x172 * x9,
            x117 * x168 * x9,
            x122 * x165 * x9,
            x128 * x157 * x28,
            x10 * x128 * x159,
            x10 * x131 * x157,
            x128 * x162 * x9,
            x131 * x159 * x9,
            x132 * x157 * x9,
            x142 * x151 * x28,
            x10 * x142 * x149,
            x144 * x146 * x5,
            x142 * x154 * x9,
            x143 * x149 * x9,
            x145 * x146,
            x147 * x182,
            x150 * x183 * x62,
            x150 * x185 * x73,
            x155 * x186 * x76,
            x155 * x185 * x99,
            x155 * x189 * x73,
            x183 * x56 * x85,
            x186 * x69 * x89,
            x185 * x69 * x87,
            x186 * x78 * x92,
            x185 * x78 * x89,
            x189 * x78 * x87,
            x191 * x56 * x73,
            x191 * x69 * x99,
            x193 * x69 * x73,
            x191 * x76 * x78,
            x193 * x78 * x99,
            x196 * x73 * x78,
            x104 * x186 * x43,
            x107 * x186 * x48,
            x104 * x185 * x48,
            x112 * x186 * x77,
            x107 * x185 * x77,
            x104 * x189 * x77,
            x191 * x43 * x87,
            x191 * x48 * x89,
            x193 * x48 * x87,
            x191 * x77 * x92,
            x193 * x77 * x89,
            x196 * x77 * x87,
            x199 * x43 * x73,
            x199 * x48 * x99,
            x202 * x48 * x73,
            x199 * x76 * x77,
            x202 * x77 * x99,
            x206 * x73 * x77,
            x123 * x186 * x30,
            x126 * x186 * x46,
            x123 * x185 * x46,
            x127 * x178 * x182,
            x126 * x185 * x34,
            x123 * x189 * x34,
            x104 * x191 * x30,
            x107 * x191 * x46,
            x104 * x193 * x46,
            x112 * x191 * x34,
            x107 * x193 * x34,
            x104 * x196 * x34,
            x199 * x30 * x87,
            x199 * x46 * x89,
            x202 * x46 * x87,
            x199 * x34 * x92,
            x202 * x34 * x89,
            x179 * x206 * x85,
            x207 * x30 * x73,
            x207 * x46 * x99,
            x210 * x46 * x73,
            x207 * x34 * x76,
            x179 * x210 * x62,
            x179 * x211,
            x133 * x186 * x28,
            x134 * x135 * x182 * x5,
            x10 * x133 * x185,
            x136 * x182,
            x134 * x185 * x9,
            x133 * x189 * x9,
            x123 * x191 * x28,
            x10 * x126 * x191,
            x10 * x123 * x193,
            x127 * x191 * x9,
            x126 * x193 * x9,
            x123 * x196 * x9,
            x104 * x199 * x28,
            x10 * x107 * x199,
            x10 * x104 * x202,
            x112 * x199 * x9,
            x107 * x202 * x9,
            x104 * x206 * x9,
            x207 * x28 * x87,
            x10 * x207 * x89,
            x141 * x210 * x5,
            x207 * x9 * x92,
            x210 * x89 * x9,
            x141 * x211,
            x212 * x28 * x73,
            x140 * x212 * x5 * x62,
            x213 * x5,
            x212 * x76 * x9,
            x213 * x62,
            x140 * (x211 * x94 + x3 * (x132 + 3 * x203 + 3 * x205 + 2 * x208 + 2 * x209)),
        ]
    )


def dipole3d_43(a, A, b, B, C):
    """Cartesian 3D (gf) dipole moment integrals.
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
    x13 = 3 * x12
    x14 = x3 * x8
    x15 = x10 * x9
    x16 = x14 + x15
    x17 = x16 * x4
    x18 = x14 * x4
    x19 = 2 * x18
    x20 = x4 ** 2 * x8
    x21 = x14 + x20
    x22 = x21 * x4
    x23 = x19 + x22
    x24 = x3 * (x13 + 3 * x17 + x23)
    x25 = 3 * x14
    x26 = x20 + x25
    x27 = x3 * (2 * x15 + x26)
    x28 = x12 + x17
    x29 = x28 * x4
    x30 = x27 + x29
    x31 = x2 * x30
    x32 = x24 + x31
    x33 = x2 * x32
    x34 = x2 * x28
    x35 = x3 * (3 * x20 + x25)
    x36 = x2 * x23
    x37 = x35 + x36
    x38 = x3 * (4 * x27 + x29 + 3 * x34 + x37)
    x39 = x33 + x38
    x40 = x2 * x21
    x41 = x3 * (8 * x18 + x22 + 3 * x40)
    x42 = x2 * x37
    x43 = x41 + x42
    x44 = x16 * x2
    x45 = 2 * x44
    x46 = x19 + x40
    x47 = x3 * (x13 + x17 + x45 + x46)
    x48 = x27 + x34
    x49 = x2 * x48
    x50 = 3 * x47 + 3 * x49
    x51 = x2 * x39 + x3 * (2 * x24 + 2 * x31 + x43 + x50)
    x52 = x2 * x9
    x53 = x11 * x2
    x54 = x3 * (x15 + x25 + x52 + x53)
    x55 = x12 + x44
    x56 = x2 * x55
    x57 = 2 * x52
    x58 = x3 * (x26 + x57)
    x59 = x2 * x46
    x60 = x58 + x59
    x61 = x3 * (2 * x27 + 2 * x34 + 2 * x54 + 2 * x56 + x60)
    x62 = x47 + x49
    x63 = x2 * x62
    x64 = 3 * x58 + 3 * x59
    x65 = x2 * x43 + x3 * (2 * x35 + 2 * x36 + x64)
    x66 = numpy.exp(-x5 * (A[1] - B[1]) ** 2)
    x67 = numpy.exp(-x5 * (A[2] - B[2]) ** 2)
    x68 = numpy.pi * x0 * x67
    x69 = x66 * x68
    x70 = -x0 * (a * A[1] + b * B[1])
    x71 = -x70 - B[1]
    x72 = x61 + x63
    x73 = x2 * x8
    x74 = x3 * (x73 + x9)
    x75 = x14 + x52
    x76 = x2 * x75
    x77 = x74 + x76
    x78 = x3 * (x11 + x73)
    x79 = x14 + x53
    x80 = x2 * x79
    x81 = x78 + x80
    x82 = x3 * (2 * x12 + x45 + x77 + x81)
    x83 = x54 + x56
    x84 = x2 * x83
    x85 = x2 * x60
    x86 = x3 * (4 * x18 + 2 * x40 + 2 * x74 + 2 * x76)
    x87 = x85 + x86
    x88 = x69 * (x2 * x72 + x3 * (x50 + 2 * x82 + 2 * x84 + x87))
    x89 = -x0 * (a * A[2] + b * B[2])
    x90 = -x89 - B[2]
    x91 = x66 * x7
    x92 = x3 * x91
    x93 = x71 ** 2 * x91
    x94 = x92 + x93
    x95 = x82 + x84
    x96 = x2 ** 2 * x8
    x97 = x25 + x96
    x98 = x2 * x81 + x3 * (2 * x53 + x97)
    x99 = x3 * (x57 + x97)
    x100 = x2 * x77
    x101 = x100 + x99
    x102 = x2 * x95 + x3 * (x101 + 3 * x54 + 3 * x56 + x98)
    x103 = x67 * x7
    x104 = x69 * x90
    x105 = x103 * x3
    x106 = x103 * x90 ** 2
    x107 = x105 + x106
    x108 = x71 * x92
    x109 = 2 * x108
    x110 = x71 * x94
    x111 = x109 + x110
    x112 = x14 + x96
    x113 = x112 * x2 + 2 * x14 * x2
    x114 = x2 * x98 + x3 * (x113 + 3 * x78 + 3 * x80)
    x115 = x103 * x90
    x116 = x71 * x91
    x117 = x105 * x90
    x118 = 2 * x117
    x119 = x107 * x90
    x120 = x118 + x119
    x121 = -x70 - A[1]
    x122 = x51 * x69
    x123 = x121 * x91
    x124 = x123 * x71
    x125 = x124 + x92
    x126 = x121 * x94
    x127 = x109 + x126
    x128 = 3 * x92
    x129 = x3 * (x128 + 3 * x93)
    x130 = x111 * x121
    x131 = x129 + x130
    x132 = -x89 - A[2]
    x133 = x132 * x69
    x134 = x103 * x132
    x135 = x134 * x90
    x136 = x105 + x135
    x137 = x107 * x132
    x138 = x118 + x137
    x139 = 3 * x105
    x140 = x3 * (3 * x106 + x139)
    x141 = x120 * x132
    x142 = x140 + x141
    x143 = x121 ** 2 * x91
    x144 = x143 + x92
    x145 = x3 * (x116 + x123)
    x146 = x121 * x125
    x147 = x145 + x146
    x148 = 2 * x124 + x128
    x149 = x3 * (x148 + x93)
    x150 = x121 * x127
    x151 = x149 + x150
    x152 = x3 * (8 * x108 + x110 + 3 * x126)
    x153 = x121 * x131
    x154 = x152 + x153
    x155 = x103 * x132 ** 2
    x156 = x105 + x155
    x157 = x3 * (x115 + x134)
    x158 = x132 * x136
    x159 = x157 + x158
    x160 = 2 * x135 + x139
    x161 = x3 * (x106 + x160)
    x162 = x132 * x138
    x163 = x161 + x162
    x164 = x3 * (8 * x117 + x119 + 3 * x137)
    x165 = x132 * x142
    x166 = x164 + x165
    x167 = x121 * x144 + 2 * x121 * x92
    x168 = x3 * (x143 + x148)
    x169 = x121 * x147
    x170 = x168 + x169
    x171 = x121 * x151
    x172 = x3 * (4 * x108 + 2 * x126 + 2 * x145 + 2 * x146)
    x173 = x171 + x172
    x174 = 3 * x149 + 3 * x150
    x175 = x121 * x154 + x3 * (2 * x129 + 2 * x130 + x174)
    x176 = 2 * x105 * x132 + x132 * x156
    x177 = x3 * (x155 + x160)
    x178 = x132 * x159
    x179 = x177 + x178
    x180 = x132 * x163
    x181 = x3 * (4 * x117 + 2 * x137 + 2 * x157 + 2 * x158)
    x182 = x180 + x181
    x183 = 3 * x161 + 3 * x162
    x184 = x132 * x166 + x3 * (2 * x140 + 2 * x141 + x183)
    x185 = x121 * x167 + x3 * (x128 + 3 * x143)
    x186 = x121 * x170 + x3 * (3 * x145 + 3 * x146 + x167)
    x187 = x121 * x173 + x3 * (2 * x168 + 2 * x169 + x174)
    x188 = x6 * x68
    x189 = x188 * (x121 * x175 + x3 * (3 * x152 + 3 * x153 + 3 * x171 + 3 * x172))
    x190 = x10 * x188
    x191 = numpy.pi * x0 * x6 * x66
    x192 = x10 * x191
    x193 = x132 * x176 + x3 * (x139 + 3 * x155)
    x194 = x132 * x179 + x3 * (3 * x157 + 3 * x158 + x176)
    x195 = x132 * x182 + x3 * (2 * x177 + 2 * x178 + x183)
    x196 = x191 * (x132 * x184 + x3 * (3 * x164 + 3 * x165 + 3 * x180 + 3 * x181))
    x197 = -x70 - C[1]
    x198 = x69 * (x2 * x65 + x3 * (3 * x41 + 3 * x42 + 3 * x85 + 3 * x86))
    x199 = x116 * x197
    x200 = x199 + x92
    x201 = x2 * x87 + x3 * (2 * x100 + x64 + 2 * x99)
    x202 = x197 * x91
    x203 = x3 * (x116 + x202)
    x204 = x200 * x71
    x205 = x203 + x204
    x206 = x101 * x2 + x3 * (x113 + 3 * x74 + 3 * x76)
    x207 = x3 * (x128 + 2 * x199 + x93)
    x208 = x205 * x71
    x209 = x207 + x208
    x210 = x113 * x2 + x3 * (x25 + 3 * x96)
    x211 = x123 * x197
    x212 = x211 + x92
    x213 = x121 * x200
    x214 = x203 + x213
    x215 = x121 * x205
    x216 = x207 + x215
    x217 = 3 * x203
    x218 = x3 * (x111 + 3 * x204 + x217)
    x219 = x121 * x209
    x220 = x218 + x219
    x221 = x3 * (x123 + x202)
    x222 = x121 * x212
    x223 = x221 + x222
    x224 = x3 * (x124 + x128 + x199 + x211)
    x225 = x121 * x214
    x226 = x224 + x225
    x227 = x121 * x216
    x228 = 2 * x213
    x229 = x3 * (x127 + x204 + x217 + x228)
    x230 = x227 + x229
    x231 = x121 * x220
    x232 = x3 * (x131 + 4 * x207 + x208 + 3 * x215)
    x233 = x231 + x232
    x234 = x121 * x223 + x3 * (x128 + x143 + 2 * x211)
    x235 = x121 * x226
    x236 = x3 * (x147 + 2 * x203 + x223 + x228)
    x237 = x235 + x236
    x238 = x121 * x230
    x239 = x3 * (x151 + 2 * x207 + 2 * x215 + 2 * x224 + 2 * x225)
    x240 = x238 + x239
    x241 = 3 * x227 + 3 * x229
    x242 = x121 * x233 + x3 * (x154 + 2 * x218 + 2 * x219 + x241)
    x243 = x188 * x242
    x244 = x188 * x2
    x245 = x191 * x197
    x246 = x121 * x234 + x3 * (x167 + 3 * x221 + 3 * x222)
    x247 = x121 * x237 + x3 * (x170 + 3 * x224 + 3 * x225 + x234)
    x248 = x188 * (x121 * x240 + x3 * (x173 + 2 * x235 + 2 * x236 + x241))
    x249 = x188 * x4
    x250 = -x89 - C[2]
    x251 = x250 * x69
    x252 = x115 * x250
    x253 = x105 + x252
    x254 = x103 * x250
    x255 = x3 * (x115 + x254)
    x256 = x253 * x90
    x257 = x255 + x256
    x258 = x3 * (x106 + x139 + 2 * x252)
    x259 = x257 * x90
    x260 = x258 + x259
    x261 = x134 * x250
    x262 = x105 + x261
    x263 = x132 * x253
    x264 = x255 + x263
    x265 = x132 * x257
    x266 = x258 + x265
    x267 = 3 * x255
    x268 = x3 * (x120 + 3 * x256 + x267)
    x269 = x132 * x260
    x270 = x268 + x269
    x271 = x3 * (x134 + x254)
    x272 = x132 * x262
    x273 = x271 + x272
    x274 = x3 * (x135 + x139 + x252 + x261)
    x275 = x132 * x264
    x276 = x274 + x275
    x277 = x132 * x266
    x278 = 2 * x263
    x279 = x3 * (x138 + x256 + x267 + x278)
    x280 = x277 + x279
    x281 = x132 * x270
    x282 = x3 * (x142 + 4 * x258 + x259 + 3 * x265)
    x283 = x281 + x282
    x284 = x191 * x2
    x285 = x132 * x273 + x3 * (x139 + x155 + 2 * x261)
    x286 = x132 * x276
    x287 = x3 * (x159 + 2 * x255 + x273 + x278)
    x288 = x286 + x287
    x289 = x132 * x280
    x290 = x3 * (x163 + 2 * x258 + 2 * x265 + 2 * x274 + 2 * x275)
    x291 = x289 + x290
    x292 = 3 * x277 + 3 * x279
    x293 = x132 * x283 + x3 * (x166 + 2 * x268 + 2 * x269 + x292)
    x294 = x191 * x293
    x295 = x191 * x4
    x296 = x132 * x285 + x3 * (x176 + 3 * x271 + 3 * x272)
    x297 = x132 * x288 + x3 * (x179 + 3 * x274 + 3 * x275 + x285)
    x298 = x191 * (x132 * x291 + x3 * (x182 + 2 * x286 + 2 * x287 + x292))

    # 450 item(s)
    return numpy.array(
        [
            x69 * (x2 * x51 + x3 * (3 * x33 + 3 * x38 + 3 * x61 + 3 * x63 + x65)),
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
            x103 * x125 * x72,
            x104 * x121 * x72,
            x103 * x127 * x95,
            x115 * x125 * x95,
            x107 * x123 * x95,
            x103 * x131 * x98,
            x115 * x127 * x98,
            x107 * x125 * x98,
            x120 * x123 * x98,
            x122 * x132,
            x133 * x71 * x72,
            x136 * x72 * x91,
            x134 * x94 * x95,
            x116 * x136 * x95,
            x138 * x91 * x95,
            x111 * x134 * x98,
            x136 * x94 * x98,
            x116 * x138 * x98,
            x142 * x91 * x98,
            x103 * x144 * x39,
            x103 * x147 * x62,
            x115 * x144 * x62,
            x103 * x151 * x83,
            x115 * x147 * x83,
            x107 * x144 * x83,
            x103 * x154 * x81,
            x115 * x151 * x81,
            x107 * x147 * x81,
            x120 * x144 * x81,
            x121 * x133 * x39,
            x125 * x134 * x62,
            x123 * x136 * x62,
            x127 * x134 * x83,
            x125 * x136 * x83,
            x123 * x138 * x83,
            x131 * x134 * x81,
            x127 * x136 * x81,
            x125 * x138 * x81,
            x123 * x142 * x81,
            x156 * x39 * x91,
            x116 * x156 * x62,
            x159 * x62 * x91,
            x156 * x83 * x94,
            x116 * x159 * x83,
            x163 * x83 * x91,
            x111 * x156 * x81,
            x159 * x81 * x94,
            x116 * x163 * x81,
            x166 * x81 * x91,
            x103 * x167 * x32,
            x103 * x170 * x48,
            x115 * x167 * x48,
            x103 * x173 * x55,
            x115 * x170 * x55,
            x107 * x167 * x55,
            x103 * x175 * x79,
            x115 * x173 * x79,
            x107 * x170 * x79,
            x120 * x167 * x79,
            x134 * x144 * x32,
            x134 * x147 * x48,
            x136 * x144 * x48,
            x134 * x151 * x55,
            x136 * x147 * x55,
            x138 * x144 * x55,
            x134 * x154 * x79,
            x136 * x151 * x79,
            x138 * x147 * x79,
            x142 * x144 * x79,
            x123 * x156 * x32,
            x125 * x156 * x48,
            x123 * x159 * x48,
            x127 * x156 * x55,
            x125 * x159 * x55,
            x123 * x163 * x55,
            x131 * x156 * x79,
            x127 * x159 * x79,
            x125 * x163 * x79,
            x123 * x166 * x79,
            x176 * x32 * x91,
            x116 * x176 * x48,
            x179 * x48 * x91,
            x176 * x55 * x94,
            x116 * x179 * x55,
            x182 * x55 * x91,
            x111 * x176 * x79,
            x179 * x79 * x94,
            x116 * x182 * x79,
            x184 * x79 * x91,
            x103 * x185 * x30,
            x103 * x186 * x28,
            x115 * x185 * x28,
            x103 * x16 * x187,
            x115 * x16 * x186,
            x107 * x16 * x185,
            x10 * x189,
            x187 * x190 * x90,
            x107 * x11 * x186,
            x11 * x120 * x185,
            x134 * x167 * x30,
            x134 * x170 * x28,
            x136 * x167 * x28,
            x134 * x16 * x173,
            x136 * x16 * x170,
            x138 * x16 * x167,
            x132 * x175 * x190,
            x11 * x136 * x173,
            x11 * x138 * x170,
            x11 * x142 * x167,
            x144 * x156 * x30,
            x147 * x156 * x28,
            x144 * x159 * x28,
            x151 * x156 * x16,
            x147 * x159 * x16,
            x144 * x16 * x163,
            x11 * x154 * x156,
            x11 * x151 * x159,
            x11 * x147 * x163,
            x11 * x144 * x166,
            x123 * x176 * x30,
            x125 * x176 * x28,
            x123 * x179 * x28,
            x127 * x16 * x176,
            x125 * x16 * x179,
            x123 * x16 * x182,
            x11 * x131 * x176,
            x11 * x127 * x179,
            x11 * x125 * x182,
            x121 * x184 * x192,
            x193 * x30 * x91,
            x116 * x193 * x28,
            x194 * x28 * x91,
            x16 * x193 * x94,
            x116 * x16 * x194,
            x16 * x195 * x91,
            x11 * x111 * x193,
            x11 * x194 * x94,
            x192 * x195 * x71,
            x10 * x196,
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
            x103 * x212 * x65,
            x103 * x214 * x87,
            x115 * x212 * x87,
            x101 * x103 * x216,
            x101 * x115 * x214,
            x101 * x107 * x212,
            x103 * x113 * x220,
            x113 * x115 * x216,
            x107 * x113 * x214,
            x113 * x120 * x212,
            x133 * x197 * x65,
            x134 * x200 * x87,
            x136 * x202 * x87,
            x101 * x134 * x205,
            x101 * x136 * x200,
            x101 * x138 * x202,
            x113 * x134 * x209,
            x113 * x136 * x205,
            x113 * x138 * x200,
            x113 * x142 * x202,
            x103 * x223 * x43,
            x103 * x226 * x60,
            x115 * x223 * x60,
            x103 * x230 * x77,
            x115 * x226 * x77,
            x107 * x223 * x77,
            x103 * x112 * x233,
            x112 * x115 * x230,
            x107 * x112 * x226,
            x112 * x120 * x223,
            x134 * x212 * x43,
            x134 * x214 * x60,
            x136 * x212 * x60,
            x134 * x216 * x77,
            x136 * x214 * x77,
            x138 * x212 * x77,
            x112 * x134 * x220,
            x112 * x136 * x216,
            x112 * x138 * x214,
            x112 * x142 * x212,
            x156 * x202 * x43,
            x156 * x200 * x60,
            x159 * x202 * x60,
            x156 * x205 * x77,
            x159 * x200 * x77,
            x163 * x202 * x77,
            x112 * x156 * x209,
            x112 * x159 * x205,
            x112 * x163 * x200,
            x112 * x166 * x202,
            x103 * x234 * x37,
            x103 * x237 * x46,
            x115 * x234 * x46,
            x103 * x240 * x75,
            x115 * x237 * x75,
            x107 * x234 * x75,
            x2 * x243,
            x240 * x244 * x90,
            x107 * x237 * x73,
            x120 * x234 * x73,
            x134 * x223 * x37,
            x134 * x226 * x46,
            x136 * x223 * x46,
            x134 * x230 * x75,
            x136 * x226 * x75,
            x138 * x223 * x75,
            x132 * x233 * x244,
            x136 * x230 * x73,
            x138 * x226 * x73,
            x142 * x223 * x73,
            x156 * x212 * x37,
            x156 * x214 * x46,
            x159 * x212 * x46,
            x156 * x216 * x75,
            x159 * x214 * x75,
            x163 * x212 * x75,
            x156 * x220 * x73,
            x159 * x216 * x73,
            x163 * x214 * x73,
            x166 * x212 * x73,
            x176 * x202 * x37,
            x176 * x200 * x46,
            x179 * x202 * x46,
            x176 * x205 * x75,
            x179 * x200 * x75,
            x182 * x202 * x75,
            x176 * x209 * x73,
            x179 * x205 * x73,
            x182 * x200 * x73,
            x184 * x2 * x245,
            x103 * x23 * x246,
            x103 * x21 * x247,
            x115 * x21 * x246,
            x248 * x4,
            x247 * x249 * x90,
            x107 * x246 * x9,
            x188
            * (x121 * x242 + x3 * (x175 + 3 * x231 + 3 * x232 + 3 * x238 + 3 * x239)),
            x248 * x90,
            x107 * x247 * x8,
            x120 * x246 * x8,
            x134 * x23 * x234,
            x134 * x21 * x237,
            x136 * x21 * x234,
            x132 * x240 * x249,
            x136 * x237 * x9,
            x138 * x234 * x9,
            x132 * x243,
            x136 * x240 * x8,
            x138 * x237 * x8,
            x142 * x234 * x8,
            x156 * x223 * x23,
            x156 * x21 * x226,
            x159 * x21 * x223,
            x156 * x230 * x9,
            x159 * x226 * x9,
            x163 * x223 * x9,
            x156 * x233 * x8,
            x159 * x230 * x8,
            x163 * x226 * x8,
            x166 * x223 * x8,
            x176 * x212 * x23,
            x176 * x21 * x214,
            x179 * x21 * x212,
            x176 * x216 * x9,
            x179 * x214 * x9,
            x182 * x212 * x9,
            x176 * x220 * x8,
            x179 * x216 * x8,
            x182 * x214 * x8,
            x184 * x212 * x8,
            x193 * x202 * x23,
            x193 * x200 * x21,
            x194 * x202 * x21,
            x193 * x205 * x9,
            x194 * x200 * x9,
            x195 * x245 * x4,
            x193 * x209 * x8,
            x194 * x205 * x8,
            x195 * x200 * x8,
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
            x121 * x251 * x65,
            x125 * x254 * x87,
            x123 * x253 * x87,
            x101 * x127 * x254,
            x101 * x125 * x253,
            x101 * x123 * x257,
            x113 * x131 * x254,
            x113 * x127 * x253,
            x113 * x125 * x257,
            x113 * x123 * x260,
            x262 * x65 * x91,
            x116 * x262 * x87,
            x264 * x87 * x91,
            x101 * x262 * x94,
            x101 * x116 * x264,
            x101 * x266 * x91,
            x111 * x113 * x262,
            x113 * x264 * x94,
            x113 * x116 * x266,
            x113 * x270 * x91,
            x144 * x254 * x43,
            x147 * x254 * x60,
            x144 * x253 * x60,
            x151 * x254 * x77,
            x147 * x253 * x77,
            x144 * x257 * x77,
            x112 * x154 * x254,
            x112 * x151 * x253,
            x112 * x147 * x257,
            x112 * x144 * x260,
            x123 * x262 * x43,
            x125 * x262 * x60,
            x123 * x264 * x60,
            x127 * x262 * x77,
            x125 * x264 * x77,
            x123 * x266 * x77,
            x112 * x131 * x262,
            x112 * x127 * x264,
            x112 * x125 * x266,
            x112 * x123 * x270,
            x273 * x43 * x91,
            x116 * x273 * x60,
            x276 * x60 * x91,
            x273 * x77 * x94,
            x116 * x276 * x77,
            x280 * x77 * x91,
            x111 * x112 * x273,
            x112 * x276 * x94,
            x112 * x116 * x280,
            x112 * x283 * x91,
            x167 * x254 * x37,
            x170 * x254 * x46,
            x167 * x253 * x46,
            x173 * x254 * x75,
            x170 * x253 * x75,
            x167 * x257 * x75,
            x175 * x244 * x250,
            x173 * x253 * x73,
            x170 * x257 * x73,
            x167 * x260 * x73,
            x144 * x262 * x37,
            x147 * x262 * x46,
            x144 * x264 * x46,
            x151 * x262 * x75,
            x147 * x264 * x75,
            x144 * x266 * x75,
            x154 * x262 * x73,
            x151 * x264 * x73,
            x147 * x266 * x73,
            x144 * x270 * x73,
            x123 * x273 * x37,
            x125 * x273 * x46,
            x123 * x276 * x46,
            x127 * x273 * x75,
            x125 * x276 * x75,
            x123 * x280 * x75,
            x131 * x273 * x73,
            x127 * x276 * x73,
            x125 * x280 * x73,
            x121 * x283 * x284,
            x285 * x37 * x91,
            x116 * x285 * x46,
            x288 * x46 * x91,
            x285 * x75 * x94,
            x116 * x288 * x75,
            x291 * x75 * x91,
            x111 * x285 * x73,
            x288 * x73 * x94,
            x284 * x291 * x71,
            x2 * x294,
            x185 * x23 * x254,
            x186 * x21 * x254,
            x185 * x21 * x253,
            x187 * x249 * x250,
            x186 * x253 * x9,
            x185 * x257 * x9,
            x189 * x250,
            x187 * x253 * x8,
            x186 * x257 * x8,
            x185 * x260 * x8,
            x167 * x23 * x262,
            x170 * x21 * x262,
            x167 * x21 * x264,
            x173 * x262 * x9,
            x170 * x264 * x9,
            x167 * x266 * x9,
            x175 * x262 * x8,
            x173 * x264 * x8,
            x170 * x266 * x8,
            x167 * x270 * x8,
            x144 * x23 * x273,
            x147 * x21 * x273,
            x144 * x21 * x276,
            x151 * x273 * x9,
            x147 * x276 * x9,
            x144 * x280 * x9,
            x154 * x273 * x8,
            x151 * x276 * x8,
            x147 * x280 * x8,
            x144 * x283 * x8,
            x123 * x23 * x285,
            x125 * x21 * x285,
            x123 * x21 * x288,
            x127 * x285 * x9,
            x125 * x288 * x9,
            x121 * x291 * x295,
            x131 * x285 * x8,
            x127 * x288 * x8,
            x125 * x291 * x8,
            x121 * x294,
            x23 * x296 * x91,
            x116 * x21 * x296,
            x21 * x297 * x91,
            x296 * x9 * x94,
            x295 * x297 * x71,
            x298 * x4,
            x111 * x296 * x8,
            x297 * x8 * x94,
            x298 * x71,
            x191
            * (x132 * x293 + x3 * (x184 + 3 * x281 + 3 * x282 + 3 * x289 + 3 * x290)),
        ]
    )


def dipole3d_44(a, A, b, B, C):
    """Cartesian 3D (gg) dipole moment integrals.
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
    x12 = x5 ** 2 * x9
    x13 = x3 * x9
    x14 = 3 * x13
    x15 = x12 + x14
    x16 = x3 * (2 * x11 + x15)
    x17 = 4 * x16
    x18 = x4 * x9
    x19 = x3 * (x10 + x18)
    x20 = x11 + x13
    x21 = x20 * x5
    x22 = x19 + x21
    x23 = x22 * x5
    x24 = x3 * (3 * x12 + x14)
    x25 = x13 * x5
    x26 = 2 * x25
    x27 = x12 + x13
    x28 = x27 * x5
    x29 = x26 + x28
    x30 = x29 * x5
    x31 = x24 + x30
    x32 = x3 * (x17 + 4 * x23 + x31)
    x33 = 3 * x19
    x34 = x3 * (3 * x21 + x29 + x33)
    x35 = x16 + x23
    x36 = x35 * x5
    x37 = x34 + x36
    x38 = x2 * x37
    x39 = x32 + x38
    x40 = x2 * x39
    x41 = x2 * x35
    x42 = 8 * x25
    x43 = x3 * (4 * x28 + x42)
    x44 = x2 * x31
    x45 = x43 + x44
    x46 = x3 * (5 * x34 + x36 + 4 * x41 + x45)
    x47 = x40 + x46
    x48 = x2 * x22
    x49 = x2 * x29
    x50 = x24 + x49
    x51 = x3 * (x17 + x23 + 3 * x48 + x50)
    x52 = x34 + x41
    x53 = x2 * x52
    x54 = x3 * (5 * x24 + x30 + 4 * x49)
    x55 = x2 * x45
    x56 = x54 + x55
    x57 = x2 * x47 + x3 * (2 * x32 + 2 * x38 + 4 * x51 + 4 * x53 + x56)
    x58 = x2 * x27
    x59 = x3 * (x28 + x42 + 3 * x58)
    x60 = x2 * x50
    x61 = x59 + x60
    x62 = x2 * x20
    x63 = 2 * x62
    x64 = x26 + x58
    x65 = x3 * (x21 + x33 + x63 + x64)
    x66 = x16 + x48
    x67 = x2 * x66
    x68 = 3 * x65 + 3 * x67
    x69 = x3 * (2 * x34 + 2 * x41 + x61 + x68)
    x70 = x51 + x53
    x71 = x2 * x70
    x72 = x2 * x56 + x3 * (2 * x43 + 2 * x44 + 4 * x59 + 4 * x60)
    x73 = numpy.exp(-x6 * (A[1] - B[1]) ** 2)
    x74 = numpy.exp(-x6 * (A[2] - B[2]) ** 2)
    x75 = numpy.pi * x0 * x74
    x76 = x73 * x75
    x77 = -x0 * (a * A[1] + b * B[1])
    x78 = -x77 - B[1]
    x79 = x69 + x71
    x80 = x10 * x2
    x81 = x18 * x2
    x82 = x3 * (x11 + x14 + x80 + x81)
    x83 = x19 + x62
    x84 = x2 * x83
    x85 = 2 * x80
    x86 = x3 * (x15 + x85)
    x87 = x2 * x64
    x88 = x86 + x87
    x89 = x3 * (2 * x16 + 2 * x48 + 2 * x82 + 2 * x84 + x88)
    x90 = x65 + x67
    x91 = x2 * x90
    x92 = x2 * x61
    x93 = 3 * x86 + 3 * x87
    x94 = x3 * (2 * x24 + 2 * x49 + x93)
    x95 = x92 + x94
    x96 = x76 * (x2 * x79 + x3 * (3 * x51 + 3 * x53 + 3 * x89 + 3 * x91 + x95))
    x97 = -x0 * (a * A[2] + b * B[2])
    x98 = -x97 - B[2]
    x99 = x73 * x8
    x100 = x3 * x99
    x101 = x78 ** 2 * x99
    x102 = x100 + x101
    x103 = x89 + x91
    x104 = x2 * x9
    x105 = x3 * (x10 + x104)
    x106 = x13 + x80
    x107 = x106 * x2
    x108 = x105 + x107
    x109 = x3 * (x104 + x18)
    x110 = x13 + x81
    x111 = x110 * x2
    x112 = x109 + x111
    x113 = x3 * (x108 + x112 + 2 * x19 + x63)
    x114 = x82 + x84
    x115 = x114 * x2
    x116 = x2 * x88
    x117 = x3 * (2 * x105 + 2 * x107 + 4 * x25 + 2 * x58)
    x118 = x116 + x117
    x119 = x103 * x2 + x3 * (2 * x113 + 2 * x115 + x118 + x68)
    x120 = x74 * x8
    x121 = x76 * x98
    x122 = x120 * x3
    x123 = x120 * x98 ** 2
    x124 = x122 + x123
    x125 = x100 * x78
    x126 = 2 * x125
    x127 = x102 * x78
    x128 = x126 + x127
    x129 = x113 + x115
    x130 = x2 ** 2 * x9
    x131 = x130 + x14
    x132 = x112 * x2 + x3 * (x131 + 2 * x81)
    x133 = x3 * (x131 + x85)
    x134 = x108 * x2
    x135 = x133 + x134
    x136 = x129 * x2 + x3 * (x132 + x135 + 3 * x82 + 3 * x84)
    x137 = x120 * x98
    x138 = x78 * x99
    x139 = x122 * x98
    x140 = 2 * x139
    x141 = x124 * x98
    x142 = x140 + x141
    x143 = 3 * x100
    x144 = x3 * (3 * x101 + x143)
    x145 = x128 * x78
    x146 = x144 + x145
    x147 = x13 + x130
    x148 = 2 * x13 * x2 + x147 * x2
    x149 = x132 * x2 + x3 * (3 * x109 + 3 * x111 + x148)
    x150 = 3 * x122
    x151 = x3 * (3 * x123 + x150)
    x152 = x142 * x98
    x153 = x151 + x152
    x154 = -x77 - A[1]
    x155 = x57 * x76
    x156 = x154 * x99
    x157 = x156 * x78
    x158 = x100 + x157
    x159 = x102 * x154
    x160 = x126 + x159
    x161 = x128 * x154
    x162 = x144 + x161
    x163 = 8 * x125
    x164 = x3 * (4 * x127 + x163)
    x165 = x146 * x154
    x166 = x164 + x165
    x167 = -x97 - A[2]
    x168 = x167 * x76
    x169 = x120 * x167
    x170 = x169 * x98
    x171 = x122 + x170
    x172 = x124 * x167
    x173 = x140 + x172
    x174 = x142 * x167
    x175 = x151 + x174
    x176 = 8 * x139
    x177 = x3 * (4 * x141 + x176)
    x178 = x153 * x167
    x179 = x177 + x178
    x180 = x154 ** 2 * x99
    x181 = x100 + x180
    x182 = x3 * (x138 + x156)
    x183 = x154 * x158
    x184 = x182 + x183
    x185 = x143 + 2 * x157
    x186 = x3 * (x101 + x185)
    x187 = x154 * x160
    x188 = x186 + x187
    x189 = x3 * (x127 + 3 * x159 + x163)
    x190 = x154 * x162
    x191 = x189 + x190
    x192 = x3 * (5 * x144 + x145 + 4 * x161)
    x193 = x154 * x166
    x194 = x192 + x193
    x195 = x120 * x167 ** 2
    x196 = x122 + x195
    x197 = x3 * (x137 + x169)
    x198 = x167 * x171
    x199 = x197 + x198
    x200 = x150 + 2 * x170
    x201 = x3 * (x123 + x200)
    x202 = x167 * x173
    x203 = x201 + x202
    x204 = x3 * (x141 + 3 * x172 + x176)
    x205 = x167 * x175
    x206 = x204 + x205
    x207 = x3 * (5 * x151 + x152 + 4 * x174)
    x208 = x167 * x179
    x209 = x207 + x208
    x210 = 2 * x100 * x154 + x154 * x181
    x211 = x3 * (x180 + x185)
    x212 = x154 * x184
    x213 = x211 + x212
    x214 = x154 * x188
    x215 = x3 * (4 * x125 + 2 * x159 + 2 * x182 + 2 * x183)
    x216 = x214 + x215
    x217 = x154 * x191
    x218 = 3 * x186 + 3 * x187
    x219 = x3 * (2 * x144 + 2 * x161 + x218)
    x220 = x217 + x219
    x221 = x154 * x194 + x3 * (2 * x164 + 2 * x165 + 4 * x189 + 4 * x190)
    x222 = 2 * x122 * x167 + x167 * x196
    x223 = x3 * (x195 + x200)
    x224 = x167 * x199
    x225 = x223 + x224
    x226 = x167 * x203
    x227 = x3 * (4 * x139 + 2 * x172 + 2 * x197 + 2 * x198)
    x228 = x226 + x227
    x229 = x167 * x206
    x230 = 3 * x201 + 3 * x202
    x231 = x3 * (2 * x151 + 2 * x174 + x230)
    x232 = x229 + x231
    x233 = x167 * x209 + x3 * (2 * x177 + 2 * x178 + 4 * x204 + 4 * x205)
    x234 = x154 * x210 + x3 * (x143 + 3 * x180)
    x235 = x154 * x213 + x3 * (3 * x182 + 3 * x183 + x210)
    x236 = x154 * x216 + x3 * (2 * x211 + 2 * x212 + x218)
    x237 = x154 * x220 + x3 * (3 * x189 + 3 * x190 + 3 * x214 + 3 * x215)
    x238 = x7 * x75
    x239 = x238 * (x154 * x221 + x3 * (3 * x192 + 3 * x193 + 4 * x217 + 4 * x219))
    x240 = x238 * x4
    x241 = numpy.pi * x0 * x7 * x73
    x242 = x241 * x4
    x243 = x167 * x222 + x3 * (x150 + 3 * x195)
    x244 = x167 * x225 + x3 * (3 * x197 + 3 * x198 + x222)
    x245 = x167 * x228 + x3 * (2 * x223 + 2 * x224 + x230)
    x246 = x167 * x232 + x3 * (3 * x204 + 3 * x205 + 3 * x226 + 3 * x227)
    x247 = x241 * (x167 * x233 + x3 * (3 * x207 + 3 * x208 + 4 * x229 + 4 * x231))
    x248 = -x77 - C[1]
    x249 = x76 * (x2 * x72 + x3 * (3 * x54 + 3 * x55 + 4 * x92 + 4 * x94))
    x250 = x138 * x248
    x251 = x100 + x250
    x252 = x2 * x95 + x3 * (3 * x116 + 3 * x117 + 3 * x59 + 3 * x60)
    x253 = x248 * x99
    x254 = x3 * (x138 + x253)
    x255 = x251 * x78
    x256 = x254 + x255
    x257 = x118 * x2 + x3 * (2 * x133 + 2 * x134 + x93)
    x258 = x3 * (x101 + x143 + 2 * x250)
    x259 = x256 * x78
    x260 = x258 + x259
    x261 = x135 * x2 + x3 * (3 * x105 + 3 * x107 + x148)
    x262 = 3 * x254
    x263 = x3 * (x128 + 3 * x255 + x262)
    x264 = x260 * x78
    x265 = x263 + x264
    x266 = x148 * x2 + x3 * (3 * x130 + x14)
    x267 = x156 * x248
    x268 = x100 + x267
    x269 = x154 * x251
    x270 = x254 + x269
    x271 = x154 * x256
    x272 = x258 + x271
    x273 = x154 * x260
    x274 = x263 + x273
    x275 = 4 * x258
    x276 = x3 * (x146 + 4 * x259 + x275)
    x277 = x154 * x265
    x278 = x276 + x277
    x279 = x3 * (x156 + x253)
    x280 = x154 * x268
    x281 = x279 + x280
    x282 = x3 * (x143 + x157 + x250 + x267)
    x283 = x154 * x270
    x284 = x282 + x283
    x285 = x154 * x272
    x286 = 2 * x269
    x287 = x3 * (x160 + x255 + x262 + x286)
    x288 = x285 + x287
    x289 = x154 * x274
    x290 = x3 * (x162 + x259 + 3 * x271 + x275)
    x291 = x289 + x290
    x292 = x154 * x278
    x293 = x3 * (x166 + 5 * x263 + x264 + 4 * x273)
    x294 = x292 + x293
    x295 = x154 * x281 + x3 * (x143 + x180 + 2 * x267)
    x296 = x154 * x284
    x297 = x3 * (x184 + 2 * x254 + x281 + x286)
    x298 = x296 + x297
    x299 = x154 * x288
    x300 = x3 * (x188 + 2 * x258 + 2 * x271 + 2 * x282 + 2 * x283)
    x301 = x299 + x300
    x302 = x154 * x291
    x303 = 3 * x285 + 3 * x287
    x304 = x3 * (x191 + 2 * x263 + 2 * x273 + x303)
    x305 = x302 + x304
    x306 = x154 * x294 + x3 * (x194 + 2 * x276 + 2 * x277 + 4 * x289 + 4 * x290)
    x307 = x238 * x306
    x308 = x2 * x238
    x309 = x241 * x248
    x310 = x154 * x295 + x3 * (x210 + 3 * x279 + 3 * x280)
    x311 = x154 * x298 + x3 * (x213 + 3 * x282 + 3 * x283 + x295)
    x312 = x154 * x301 + x3 * (x216 + 2 * x296 + 2 * x297 + x303)
    x313 = x238 * (x154 * x305 + x3 * (x220 + 3 * x289 + 3 * x290 + 3 * x299 + 3 * x300))
    x314 = x238 * x5
    x315 = -x97 - C[2]
    x316 = x315 * x76
    x317 = x137 * x315
    x318 = x122 + x317
    x319 = x120 * x315
    x320 = x3 * (x137 + x319)
    x321 = x318 * x98
    x322 = x320 + x321
    x323 = x3 * (x123 + x150 + 2 * x317)
    x324 = x322 * x98
    x325 = x323 + x324
    x326 = 3 * x320
    x327 = x3 * (x142 + 3 * x321 + x326)
    x328 = x325 * x98
    x329 = x327 + x328
    x330 = x169 * x315
    x331 = x122 + x330
    x332 = x167 * x318
    x333 = x320 + x332
    x334 = x167 * x322
    x335 = x323 + x334
    x336 = x167 * x325
    x337 = x327 + x336
    x338 = 4 * x323
    x339 = x3 * (x153 + 4 * x324 + x338)
    x340 = x167 * x329
    x341 = x339 + x340
    x342 = x3 * (x169 + x319)
    x343 = x167 * x331
    x344 = x342 + x343
    x345 = x3 * (x150 + x170 + x317 + x330)
    x346 = x167 * x333
    x347 = x345 + x346
    x348 = x167 * x335
    x349 = 2 * x332
    x350 = x3 * (x173 + x321 + x326 + x349)
    x351 = x348 + x350
    x352 = x167 * x337
    x353 = x3 * (x175 + x324 + 3 * x334 + x338)
    x354 = x352 + x353
    x355 = x167 * x341
    x356 = x3 * (x179 + 5 * x327 + x328 + 4 * x336)
    x357 = x355 + x356
    x358 = x2 * x241
    x359 = x167 * x344 + x3 * (x150 + x195 + 2 * x330)
    x360 = x167 * x347
    x361 = x3 * (x199 + 2 * x320 + x344 + x349)
    x362 = x360 + x361
    x363 = x167 * x351
    x364 = x3 * (x203 + 2 * x323 + 2 * x334 + 2 * x345 + 2 * x346)
    x365 = x363 + x364
    x366 = x167 * x354
    x367 = 3 * x348 + 3 * x350
    x368 = x3 * (x206 + 2 * x327 + 2 * x336 + x367)
    x369 = x366 + x368
    x370 = x167 * x357 + x3 * (x209 + 2 * x339 + 2 * x340 + 4 * x352 + 4 * x353)
    x371 = x241 * x370
    x372 = x241 * x5
    x373 = x167 * x359 + x3 * (x222 + 3 * x342 + 3 * x343)
    x374 = x167 * x362 + x3 * (x225 + 3 * x345 + 3 * x346 + x359)
    x375 = x167 * x365 + x3 * (x228 + 2 * x360 + 2 * x361 + x367)
    x376 = x241 * (x167 * x369 + x3 * (x232 + 3 * x352 + 3 * x353 + 3 * x363 + 3 * x364))

    # 675 item(s)
    return numpy.array(
        [
            x76 * (x2 * x57 + x3 * (3 * x40 + 3 * x46 + 4 * x69 + 4 * x71 + x72)),
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
            x120 * x158 * x79,
            x121 * x154 * x79,
            x103 * x120 * x160,
            x103 * x137 * x158,
            x103 * x124 * x156,
            x120 * x129 * x162,
            x129 * x137 * x160,
            x124 * x129 * x158,
            x129 * x142 * x156,
            x120 * x132 * x166,
            x132 * x137 * x162,
            x124 * x132 * x160,
            x132 * x142 * x158,
            x132 * x153 * x156,
            x155 * x167,
            x168 * x78 * x79,
            x171 * x79 * x99,
            x102 * x103 * x169,
            x103 * x138 * x171,
            x103 * x173 * x99,
            x128 * x129 * x169,
            x102 * x129 * x171,
            x129 * x138 * x173,
            x129 * x175 * x99,
            x132 * x146 * x169,
            x128 * x132 * x171,
            x102 * x132 * x173,
            x132 * x138 * x175,
            x132 * x179 * x99,
            x120 * x181 * x47,
            x120 * x184 * x70,
            x137 * x181 * x70,
            x120 * x188 * x90,
            x137 * x184 * x90,
            x124 * x181 * x90,
            x114 * x120 * x191,
            x114 * x137 * x188,
            x114 * x124 * x184,
            x114 * x142 * x181,
            x112 * x120 * x194,
            x112 * x137 * x191,
            x112 * x124 * x188,
            x112 * x142 * x184,
            x112 * x153 * x181,
            x154 * x168 * x47,
            x158 * x169 * x70,
            x156 * x171 * x70,
            x160 * x169 * x90,
            x158 * x171 * x90,
            x156 * x173 * x90,
            x114 * x162 * x169,
            x114 * x160 * x171,
            x114 * x158 * x173,
            x114 * x156 * x175,
            x112 * x166 * x169,
            x112 * x162 * x171,
            x112 * x160 * x173,
            x112 * x158 * x175,
            x112 * x156 * x179,
            x196 * x47 * x99,
            x138 * x196 * x70,
            x199 * x70 * x99,
            x102 * x196 * x90,
            x138 * x199 * x90,
            x203 * x90 * x99,
            x114 * x128 * x196,
            x102 * x114 * x199,
            x114 * x138 * x203,
            x114 * x206 * x99,
            x112 * x146 * x196,
            x112 * x128 * x199,
            x102 * x112 * x203,
            x112 * x138 * x206,
            x112 * x209 * x99,
            x120 * x210 * x39,
            x120 * x213 * x52,
            x137 * x210 * x52,
            x120 * x216 * x66,
            x137 * x213 * x66,
            x124 * x210 * x66,
            x120 * x220 * x83,
            x137 * x216 * x83,
            x124 * x213 * x83,
            x142 * x210 * x83,
            x110 * x120 * x221,
            x110 * x137 * x220,
            x110 * x124 * x216,
            x110 * x142 * x213,
            x110 * x153 * x210,
            x169 * x181 * x39,
            x169 * x184 * x52,
            x171 * x181 * x52,
            x169 * x188 * x66,
            x171 * x184 * x66,
            x173 * x181 * x66,
            x169 * x191 * x83,
            x171 * x188 * x83,
            x173 * x184 * x83,
            x175 * x181 * x83,
            x110 * x169 * x194,
            x110 * x171 * x191,
            x110 * x173 * x188,
            x110 * x175 * x184,
            x110 * x179 * x181,
            x156 * x196 * x39,
            x158 * x196 * x52,
            x156 * x199 * x52,
            x160 * x196 * x66,
            x158 * x199 * x66,
            x156 * x203 * x66,
            x162 * x196 * x83,
            x160 * x199 * x83,
            x158 * x203 * x83,
            x156 * x206 * x83,
            x110 * x166 * x196,
            x110 * x162 * x199,
            x110 * x160 * x203,
            x110 * x158 * x206,
            x110 * x156 * x209,
            x222 * x39 * x99,
            x138 * x222 * x52,
            x225 * x52 * x99,
            x102 * x222 * x66,
            x138 * x225 * x66,
            x228 * x66 * x99,
            x128 * x222 * x83,
            x102 * x225 * x83,
            x138 * x228 * x83,
            x232 * x83 * x99,
            x110 * x146 * x222,
            x110 * x128 * x225,
            x102 * x110 * x228,
            x110 * x138 * x232,
            x110 * x233 * x99,
            x120 * x234 * x37,
            x120 * x235 * x35,
            x137 * x234 * x35,
            x120 * x22 * x236,
            x137 * x22 * x235,
            x124 * x22 * x234,
            x120 * x20 * x237,
            x137 * x20 * x236,
            x124 * x20 * x235,
            x142 * x20 * x234,
            x239 * x4,
            x237 * x240 * x98,
            x124 * x18 * x236,
            x142 * x18 * x235,
            x153 * x18 * x234,
            x169 * x210 * x37,
            x169 * x213 * x35,
            x171 * x210 * x35,
            x169 * x216 * x22,
            x171 * x213 * x22,
            x173 * x210 * x22,
            x169 * x20 * x220,
            x171 * x20 * x216,
            x173 * x20 * x213,
            x175 * x20 * x210,
            x167 * x221 * x240,
            x171 * x18 * x220,
            x173 * x18 * x216,
            x175 * x18 * x213,
            x179 * x18 * x210,
            x181 * x196 * x37,
            x184 * x196 * x35,
            x181 * x199 * x35,
            x188 * x196 * x22,
            x184 * x199 * x22,
            x181 * x203 * x22,
            x191 * x196 * x20,
            x188 * x199 * x20,
            x184 * x20 * x203,
            x181 * x20 * x206,
            x18 * x194 * x196,
            x18 * x191 * x199,
            x18 * x188 * x203,
            x18 * x184 * x206,
            x18 * x181 * x209,
            x156 * x222 * x37,
            x158 * x222 * x35,
            x156 * x225 * x35,
            x160 * x22 * x222,
            x158 * x22 * x225,
            x156 * x22 * x228,
            x162 * x20 * x222,
            x160 * x20 * x225,
            x158 * x20 * x228,
            x156 * x20 * x232,
            x166 * x18 * x222,
            x162 * x18 * x225,
            x160 * x18 * x228,
            x158 * x18 * x232,
            x154 * x233 * x242,
            x243 * x37 * x99,
            x138 * x243 * x35,
            x244 * x35 * x99,
            x102 * x22 * x243,
            x138 * x22 * x244,
            x22 * x245 * x99,
            x128 * x20 * x243,
            x102 * x20 * x244,
            x138 * x20 * x245,
            x20 * x246 * x99,
            x146 * x18 * x243,
            x128 * x18 * x244,
            x102 * x18 * x245,
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
            x120 * x268 * x72,
            x120 * x270 * x95,
            x137 * x268 * x95,
            x118 * x120 * x272,
            x118 * x137 * x270,
            x118 * x124 * x268,
            x120 * x135 * x274,
            x135 * x137 * x272,
            x124 * x135 * x270,
            x135 * x142 * x268,
            x120 * x148 * x278,
            x137 * x148 * x274,
            x124 * x148 * x272,
            x142 * x148 * x270,
            x148 * x153 * x268,
            x168 * x248 * x72,
            x169 * x251 * x95,
            x171 * x253 * x95,
            x118 * x169 * x256,
            x118 * x171 * x251,
            x118 * x173 * x253,
            x135 * x169 * x260,
            x135 * x171 * x256,
            x135 * x173 * x251,
            x135 * x175 * x253,
            x148 * x169 * x265,
            x148 * x171 * x260,
            x148 * x173 * x256,
            x148 * x175 * x251,
            x148 * x179 * x253,
            x120 * x281 * x56,
            x120 * x284 * x61,
            x137 * x281 * x61,
            x120 * x288 * x88,
            x137 * x284 * x88,
            x124 * x281 * x88,
            x108 * x120 * x291,
            x108 * x137 * x288,
            x108 * x124 * x284,
            x108 * x142 * x281,
            x120 * x147 * x294,
            x137 * x147 * x291,
            x124 * x147 * x288,
            x142 * x147 * x284,
            x147 * x153 * x281,
            x169 * x268 * x56,
            x169 * x270 * x61,
            x171 * x268 * x61,
            x169 * x272 * x88,
            x171 * x270 * x88,
            x173 * x268 * x88,
            x108 * x169 * x274,
            x108 * x171 * x272,
            x108 * x173 * x270,
            x108 * x175 * x268,
            x147 * x169 * x278,
            x147 * x171 * x274,
            x147 * x173 * x272,
            x147 * x175 * x270,
            x147 * x179 * x268,
            x196 * x253 * x56,
            x196 * x251 * x61,
            x199 * x253 * x61,
            x196 * x256 * x88,
            x199 * x251 * x88,
            x203 * x253 * x88,
            x108 * x196 * x260,
            x108 * x199 * x256,
            x108 * x203 * x251,
            x108 * x206 * x253,
            x147 * x196 * x265,
            x147 * x199 * x260,
            x147 * x203 * x256,
            x147 * x206 * x251,
            x147 * x209 * x253,
            x120 * x295 * x45,
            x120 * x298 * x50,
            x137 * x295 * x50,
            x120 * x301 * x64,
            x137 * x298 * x64,
            x124 * x295 * x64,
            x106 * x120 * x305,
            x106 * x137 * x301,
            x106 * x124 * x298,
            x106 * x142 * x295,
            x2 * x307,
            x305 * x308 * x98,
            x104 * x124 * x301,
            x104 * x142 * x298,
            x104 * x153 * x295,
            x169 * x281 * x45,
            x169 * x284 * x50,
            x171 * x281 * x50,
            x169 * x288 * x64,
            x171 * x284 * x64,
            x173 * x281 * x64,
            x106 * x169 * x291,
            x106 * x171 * x288,
            x106 * x173 * x284,
            x106 * x175 * x281,
            x167 * x294 * x308,
            x104 * x171 * x291,
            x104 * x173 * x288,
            x104 * x175 * x284,
            x104 * x179 * x281,
            x196 * x268 * x45,
            x196 * x270 * x50,
            x199 * x268 * x50,
            x196 * x272 * x64,
            x199 * x270 * x64,
            x203 * x268 * x64,
            x106 * x196 * x274,
            x106 * x199 * x272,
            x106 * x203 * x270,
            x106 * x206 * x268,
            x104 * x196 * x278,
            x104 * x199 * x274,
            x104 * x203 * x272,
            x104 * x206 * x270,
            x104 * x209 * x268,
            x222 * x253 * x45,
            x222 * x251 * x50,
            x225 * x253 * x50,
            x222 * x256 * x64,
            x225 * x251 * x64,
            x228 * x253 * x64,
            x106 * x222 * x260,
            x106 * x225 * x256,
            x106 * x228 * x251,
            x106 * x232 * x253,
            x104 * x222 * x265,
            x104 * x225 * x260,
            x104 * x228 * x256,
            x104 * x232 * x251,
            x2 * x233 * x309,
            x120 * x31 * x310,
            x120 * x29 * x311,
            x137 * x29 * x310,
            x120 * x27 * x312,
            x137 * x27 * x311,
            x124 * x27 * x310,
            x313 * x5,
            x312 * x314 * x98,
            x10 * x124 * x311,
            x10 * x142 * x310,
            x238
            * (x154 * x306 + x3 * (x221 + 3 * x292 + 3 * x293 + 4 * x302 + 4 * x304)),
            x313 * x98,
            x124 * x312 * x9,
            x142 * x311 * x9,
            x153 * x310 * x9,
            x169 * x295 * x31,
            x169 * x29 * x298,
            x171 * x29 * x295,
            x169 * x27 * x301,
            x171 * x27 * x298,
            x173 * x27 * x295,
            x167 * x305 * x314,
            x10 * x171 * x301,
            x10 * x173 * x298,
            x10 * x175 * x295,
            x167 * x307,
            x171 * x305 * x9,
            x173 * x301 * x9,
            x175 * x298 * x9,
            x179 * x295 * x9,
            x196 * x281 * x31,
            x196 * x284 * x29,
            x199 * x281 * x29,
            x196 * x27 * x288,
            x199 * x27 * x284,
            x203 * x27 * x281,
            x10 * x196 * x291,
            x10 * x199 * x288,
            x10 * x203 * x284,
            x10 * x206 * x281,
            x196 * x294 * x9,
            x199 * x291 * x9,
            x203 * x288 * x9,
            x206 * x284 * x9,
            x209 * x281 * x9,
            x222 * x268 * x31,
            x222 * x270 * x29,
            x225 * x268 * x29,
            x222 * x27 * x272,
            x225 * x27 * x270,
            x228 * x268 * x27,
            x10 * x222 * x274,
            x10 * x225 * x272,
            x10 * x228 * x270,
            x10 * x232 * x268,
            x222 * x278 * x9,
            x225 * x274 * x9,
            x228 * x272 * x9,
            x232 * x270 * x9,
            x233 * x268 * x9,
            x243 * x253 * x31,
            x243 * x251 * x29,
            x244 * x253 * x29,
            x243 * x256 * x27,
            x244 * x251 * x27,
            x245 * x253 * x27,
            x10 * x243 * x260,
            x10 * x244 * x256,
            x10 * x245 * x251,
            x246 * x309 * x5,
            x243 * x265 * x9,
            x244 * x260 * x9,
            x245 * x256 * x9,
            x246 * x251 * x9,
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
            x154 * x316 * x72,
            x158 * x319 * x95,
            x156 * x318 * x95,
            x118 * x160 * x319,
            x118 * x158 * x318,
            x118 * x156 * x322,
            x135 * x162 * x319,
            x135 * x160 * x318,
            x135 * x158 * x322,
            x135 * x156 * x325,
            x148 * x166 * x319,
            x148 * x162 * x318,
            x148 * x160 * x322,
            x148 * x158 * x325,
            x148 * x156 * x329,
            x331 * x72 * x99,
            x138 * x331 * x95,
            x333 * x95 * x99,
            x102 * x118 * x331,
            x118 * x138 * x333,
            x118 * x335 * x99,
            x128 * x135 * x331,
            x102 * x135 * x333,
            x135 * x138 * x335,
            x135 * x337 * x99,
            x146 * x148 * x331,
            x128 * x148 * x333,
            x102 * x148 * x335,
            x138 * x148 * x337,
            x148 * x341 * x99,
            x181 * x319 * x56,
            x184 * x319 * x61,
            x181 * x318 * x61,
            x188 * x319 * x88,
            x184 * x318 * x88,
            x181 * x322 * x88,
            x108 * x191 * x319,
            x108 * x188 * x318,
            x108 * x184 * x322,
            x108 * x181 * x325,
            x147 * x194 * x319,
            x147 * x191 * x318,
            x147 * x188 * x322,
            x147 * x184 * x325,
            x147 * x181 * x329,
            x156 * x331 * x56,
            x158 * x331 * x61,
            x156 * x333 * x61,
            x160 * x331 * x88,
            x158 * x333 * x88,
            x156 * x335 * x88,
            x108 * x162 * x331,
            x108 * x160 * x333,
            x108 * x158 * x335,
            x108 * x156 * x337,
            x147 * x166 * x331,
            x147 * x162 * x333,
            x147 * x160 * x335,
            x147 * x158 * x337,
            x147 * x156 * x341,
            x344 * x56 * x99,
            x138 * x344 * x61,
            x347 * x61 * x99,
            x102 * x344 * x88,
            x138 * x347 * x88,
            x351 * x88 * x99,
            x108 * x128 * x344,
            x102 * x108 * x347,
            x108 * x138 * x351,
            x108 * x354 * x99,
            x146 * x147 * x344,
            x128 * x147 * x347,
            x102 * x147 * x351,
            x138 * x147 * x354,
            x147 * x357 * x99,
            x210 * x319 * x45,
            x213 * x319 * x50,
            x210 * x318 * x50,
            x216 * x319 * x64,
            x213 * x318 * x64,
            x210 * x322 * x64,
            x106 * x220 * x319,
            x106 * x216 * x318,
            x106 * x213 * x322,
            x106 * x210 * x325,
            x221 * x308 * x315,
            x104 * x220 * x318,
            x104 * x216 * x322,
            x104 * x213 * x325,
            x104 * x210 * x329,
            x181 * x331 * x45,
            x184 * x331 * x50,
            x181 * x333 * x50,
            x188 * x331 * x64,
            x184 * x333 * x64,
            x181 * x335 * x64,
            x106 * x191 * x331,
            x106 * x188 * x333,
            x106 * x184 * x335,
            x106 * x181 * x337,
            x104 * x194 * x331,
            x104 * x191 * x333,
            x104 * x188 * x335,
            x104 * x184 * x337,
            x104 * x181 * x341,
            x156 * x344 * x45,
            x158 * x344 * x50,
            x156 * x347 * x50,
            x160 * x344 * x64,
            x158 * x347 * x64,
            x156 * x351 * x64,
            x106 * x162 * x344,
            x106 * x160 * x347,
            x106 * x158 * x351,
            x106 * x156 * x354,
            x104 * x166 * x344,
            x104 * x162 * x347,
            x104 * x160 * x351,
            x104 * x158 * x354,
            x154 * x357 * x358,
            x359 * x45 * x99,
            x138 * x359 * x50,
            x362 * x50 * x99,
            x102 * x359 * x64,
            x138 * x362 * x64,
            x365 * x64 * x99,
            x106 * x128 * x359,
            x102 * x106 * x362,
            x106 * x138 * x365,
            x106 * x369 * x99,
            x104 * x146 * x359,
            x104 * x128 * x362,
            x102 * x104 * x365,
            x358 * x369 * x78,
            x2 * x371,
            x234 * x31 * x319,
            x235 * x29 * x319,
            x234 * x29 * x318,
            x236 * x27 * x319,
            x235 * x27 * x318,
            x234 * x27 * x322,
            x237 * x314 * x315,
            x10 * x236 * x318,
            x10 * x235 * x322,
            x10 * x234 * x325,
            x239 * x315,
            x237 * x318 * x9,
            x236 * x322 * x9,
            x235 * x325 * x9,
            x234 * x329 * x9,
            x210 * x31 * x331,
            x213 * x29 * x331,
            x210 * x29 * x333,
            x216 * x27 * x331,
            x213 * x27 * x333,
            x210 * x27 * x335,
            x10 * x220 * x331,
            x10 * x216 * x333,
            x10 * x213 * x335,
            x10 * x210 * x337,
            x221 * x331 * x9,
            x220 * x333 * x9,
            x216 * x335 * x9,
            x213 * x337 * x9,
            x210 * x341 * x9,
            x181 * x31 * x344,
            x184 * x29 * x344,
            x181 * x29 * x347,
            x188 * x27 * x344,
            x184 * x27 * x347,
            x181 * x27 * x351,
            x10 * x191 * x344,
            x10 * x188 * x347,
            x10 * x184 * x351,
            x10 * x181 * x354,
            x194 * x344 * x9,
            x191 * x347 * x9,
            x188 * x351 * x9,
            x184 * x354 * x9,
            x181 * x357 * x9,
            x156 * x31 * x359,
            x158 * x29 * x359,
            x156 * x29 * x362,
            x160 * x27 * x359,
            x158 * x27 * x362,
            x156 * x27 * x365,
            x10 * x162 * x359,
            x10 * x160 * x362,
            x10 * x158 * x365,
            x154 * x369 * x372,
            x166 * x359 * x9,
            x162 * x362 * x9,
            x160 * x365 * x9,
            x158 * x369 * x9,
            x154 * x371,
            x31 * x373 * x99,
            x138 * x29 * x373,
            x29 * x374 * x99,
            x102 * x27 * x373,
            x138 * x27 * x374,
            x27 * x375 * x99,
            x10 * x128 * x373,
            x10 * x102 * x374,
            x372 * x375 * x78,
            x376 * x5,
            x146 * x373 * x9,
            x128 * x374 * x9,
            x102 * x375 * x9,
            x376 * x78,
            x241
            * (x167 * x370 + x3 * (x233 + 3 * x355 + 3 * x356 + 4 * x366 + 4 * x368)),
        ]
    )
