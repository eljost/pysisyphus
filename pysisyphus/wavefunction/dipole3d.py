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
    """Cartesian 3D (ss) dipole moment integral.
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
    S = numpy.array(
        [
            x2 * (x0 * (a * A[0] + b * B[0]) - C[0]),
            x2 * (x0 * (a * A[1] + b * B[1]) - C[1]),
            x2 * (x0 * (a * A[2] + b * B[2]) - C[2]),
        ]
    )
    return S


def dipole3d_01(a, A, b, B, C):
    """Cartesian 3D (sp) dipole moment integral.
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
    x8 = -x7 - C[0]
    x9 = x3 * (-x7 - B[0])
    x10 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x11 = numpy.pi * x0 * x10
    x12 = numpy.pi ** (3 / 2) * x4
    x13 = x12 * x9
    x14 = -x0 * (a * A[1] + b * B[1])
    x15 = x2 * (-x14 - C[1])
    x16 = x0 * x10 * x15
    x17 = -x0 * (a * A[2] + b * B[2])
    x18 = x10 * (-x17 - C[2])
    x19 = x0 * x2
    x20 = -x14 - B[1]
    x21 = x19 * x3
    x22 = x12 * x20 * x21
    x23 = x10 * x8
    x24 = -x17 - B[2]
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


def dipole3d_02(a, A, b, B, C):
    """Cartesian 3D (sd) dipole moment integral.
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
    x19 = -x18 - C[1]
    x20 = x17 * (x11 * x5**2 + x12)
    x21 = -x1 * (a * A[2] + b * B[2])
    x22 = -x21 - C[2]
    x23 = -x18 - B[1]
    x24 = x13 * x17
    x25 = x0 * x3
    x26 = x14 * x25
    x27 = x14 * x3
    x28 = x23 * x27
    x29 = x19 * x28 + x26
    x30 = numpy.pi ** (3 / 2)
    x31 = x1 * x14
    x32 = x15 * x2 * x30 * x31 * x8
    x33 = -x21 - B[2]
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


def dipole3d_03(a, A, b, B, C):
    """Cartesian 3D (sf) dipole moment integral.
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
    x9 = x5 * x8**2
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
    x20 = -x19 - C[1]
    x21 = x6 + x9
    x22 = x18 * (2 * x0 * x11 + x21 * x8)
    x23 = -x1 * (a * A[2] + b * B[2])
    x24 = -x23 - C[2]
    x25 = -x19 - B[1]
    x26 = x14 * x18
    x27 = x0 * x4
    x28 = x15 * x27
    x29 = x15 * x4
    x30 = x25 * x29
    x31 = x20 * x30
    x32 = x28 + x31
    x33 = x16 * x4
    x34 = x18 * x21
    x35 = -x23 - B[2]
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


def dipole3d_04(a, A, b, B, C):
    """Cartesian 3D (sg) dipole moment integral.
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
    x27 = -x26 - C[1]
    x28 = x25 * (x0 * (3 * x16 + x19) + x18 * x3)
    x29 = -x1 * (a * A[2] + b * B[2])
    x30 = -x29 - C[2]
    x31 = -x26 - B[1]
    x32 = x21 * x25
    x33 = x0 * x6
    x34 = x22 * x33
    x35 = x22 * x6
    x36 = x31 * x35
    x37 = x27 * x36
    x38 = x34 + x37
    x39 = x23 * x6
    x40 = x18 * x25
    x41 = -x29 - B[2]
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


def dipole3d_10(a, A, b, B, C):
    """Cartesian 3D (ps) dipole moment integral.
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
    x8 = -x7 - C[0]
    x9 = x3 * (-x7 - A[0])
    x10 = numpy.exp(-x1 * (A[2] - B[2]) ** 2)
    x11 = numpy.pi * x0 * x10
    x12 = numpy.pi ** (3 / 2) * x4
    x13 = x12 * x9
    x14 = -x0 * (a * A[1] + b * B[1])
    x15 = x2 * (-x14 - C[1])
    x16 = x0 * x10 * x15
    x17 = -x0 * (a * A[2] + b * B[2])
    x18 = x10 * (-x17 - C[2])
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


def dipole3d_11(a, A, b, B, C):
    """Cartesian 3D (pp) dipole moment integral.
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
    x20 = -x19 - C[1]
    x21 = x18 * (x12 * x8 + x13)
    x22 = -x1 * (a * A[2] + b * B[2])
    x23 = -x22 - C[2]
    x24 = -x19 - B[1]
    x25 = x18 * (x11 * x12 + x13)
    x26 = x0 * x3
    x27 = x15 * x26
    x28 = x15 * x3
    x29 = x24 * x28
    x30 = x20 * x29 + x27
    x31 = x17 * x6
    x32 = x30 * x31
    x33 = x1 * x15
    x34 = numpy.pi ** (3 / 2) * x16 * x2 * x33
    x35 = x23 * x34
    x36 = x12 * x6
    x37 = -x22 - B[2]
    x38 = x34 * x37
    x39 = x16 * x26
    x40 = x16 * x3
    x41 = x37 * x40
    x42 = x23 * x41 + x39
    x43 = numpy.pi * x33
    x44 = x43 * x6
    x45 = x42 * x44
    x46 = -x19 - A[1]
    x47 = x14 * x18
    x48 = x20 * x28
    x49 = x27 + x46 * x48
    x50 = x31 * (x27 + x29 * x46)
    x51 = x6 * x9
    x52 = -x22 - A[2]
    x53 = x34 * x52
    x54 = x23 * x40
    x55 = x39 + x52 * x54
    x56 = x44 * (x39 + x41 * x52)

    # 27 item(s)
    S = numpy.array(
        [
            x18 * (x0 * (x11 + x8) + x12 * x14),
            x20 * x21,
            x21 * x23,
            x24 * x25,
            x12 * x32,
            x24 * x35 * x36,
            x25 * x37,
            x20 * x36 * x38,
            x12 * x45,
            x46 * x47,
            x17 * x49 * x7,
            x35 * x46 * x7,
            x50 * x9,
            x31 * (x0 * (x29 + x48) + x30 * x46),
            x23 * x50,
            x38 * x46 * x51,
            x31 * x37 * x49,
            x45 * x46,
            x47 * x52,
            x20 * x53 * x7,
            x43 * x55 * x7,
            x24 * x51 * x53,
            x32 * x52,
            x24 * x44 * x55,
            x56 * x9,
            x20 * x56,
            x44 * (x0 * (x41 + x54) + x42 * x52),
        ]
    )
    return S


def dipole3d_12(a, A, b, B, C):
    """Cartesian 3D (pd) dipole moment integral.
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
    x9 = x5 * x8**2
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
    x23 = -x22 - C[1]
    x24 = x6 + x9
    x25 = x21 * (2 * x0 * x11 + x13 * x24)
    x26 = -x1 * (a * A[2] + b * B[2])
    x27 = -x26 - C[2]
    x28 = -x22 - B[1]
    x29 = x21 * (x13 * x16 + x15)
    x30 = x0 * x4
    x31 = x18 * x30
    x32 = x18 * x4
    x33 = x28 * x32
    x34 = x23 * x33
    x35 = x31 + x34
    x36 = x11 * x13 + x6
    x37 = x19 * x4
    x38 = x21 * x36
    x39 = -x26 - B[2]
    x40 = x19 * x30
    x41 = x37 * x39
    x42 = x27 * x41
    x43 = x40 + x42
    x44 = x28**2 * x32
    x45 = x31 + x44
    x46 = x13 * x14 + x6
    x47 = x23 * x32
    x48 = x0 * (x33 + x47)
    x49 = x28 * x35 + x48
    x50 = x20 * x3
    x51 = x49 * x50
    x52 = x13 * x50
    x53 = x21 * x39
    x54 = numpy.pi * x1 * x18 * x3
    x55 = x13 * x54
    x56 = x37 * x39**2
    x57 = x40 + x56
    x58 = x27 * x37
    x59 = x0 * (x41 + x58)
    x60 = x39 * x43 + x59
    x61 = x54 * x60
    x62 = -x22 - A[1]
    x63 = x17 * x21
    x64 = x31 + x47 * x62
    x65 = x21 * x24
    x66 = x31 + x33 * x62
    x67 = x50 * (x35 * x62 + x48)
    x68 = x50 * x8
    x69 = x54 * x62
    x70 = x50 * (2 * x28 * x31 + x45 * x62)
    x71 = x10 * x50
    x72 = -x26 - A[2]
    x73 = x40 + x58 * x72
    x74 = x54 * x8
    x75 = x40 + x41 * x72
    x76 = x54 * (x43 * x72 + x59)
    x77 = x54 * (2 * x39 * x40 + x57 * x72)

    # 54 item(s)
    S = numpy.array(
        [
            x21 * (x0 * (2 * x12 + 3 * x6 + x9) + x13 * x17),
            x23 * x25,
            x25 * x27,
            x28 * x29,
            x35 * x36 * x37,
            x27 * x28 * x38,
            x29 * x39,
            x23 * x38 * x39,
            x32 * x36 * x43,
            x37 * x45 * x46,
            x13 * x51,
            x27 * x45 * x52,
            x28 * x46 * x53,
            x35 * x39 * x52,
            x28 * x43 * x55,
            x32 * x46 * x57,
            x23 * x55 * x57,
            x13 * x61,
            x62 * x63,
            x24 * x37 * x64,
            x27 * x62 * x65,
            x16 * x37 * x66,
            x67 * x8,
            x27 * x66 * x68,
            x16 * x53 * x62,
            x39 * x64 * x68,
            x43 * x69 * x8,
            x10 * x70,
            x50 * (x0 * (3 * x31 + 2 * x34 + x44) + x49 * x62),
            x27 * x70,
            x39 * x66 * x71,
            x39 * x67,
            x43 * x5 * x66,
            x10 * x57 * x69,
            x5 * x57 * x64,
            x61 * x62,
            x63 * x72,
            x23 * x65 * x72,
            x24 * x32 * x73,
            x16 * x21 * x28 * x72,
            x35 * x68 * x72,
            x28 * x73 * x74,
            x16 * x32 * x75,
            x23 * x74 * x75,
            x76 * x8,
            x45 * x71 * x72,
            x51 * x72,
            x45 * x5 * x73,
            x10 * x28 * x54 * x75,
            x35 * x5 * x75,
            x28 * x76,
            x10 * x77,
            x23 * x77,
            x54 * (x0 * (3 * x40 + 2 * x42 + x56) + x60 * x72),
        ]
    )
    return S


def dipole3d_13(a, A, b, B, C):
    """Cartesian 3D (pf) dipole moment integral.
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
    x17 = x3**2 * x7
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
    x30 = -x29 - C[1]
    x31 = x28 * (x0 * (3 * x17 + x21) + x19 * x20)
    x32 = -x1 * (a * A[2] + b * B[2])
    x33 = -x32 - C[2]
    x34 = -x29 - B[1]
    x35 = x28 * (x20 * x23 + x22)
    x36 = x16 + x18 * x20
    x37 = x0 * x6
    x38 = x25 * x37
    x39 = x25 * x6
    x40 = x34 * x39
    x41 = x30 * x40
    x42 = x38 + x41
    x43 = x26 * x6
    x44 = x42 * x43
    x45 = x28 * x36
    x46 = -x32 - B[2]
    x47 = x26 * x37
    x48 = x43 * x46
    x49 = x33 * x48
    x50 = x47 + x49
    x51 = x39 * x50
    x52 = x11 + x14 * x20
    x53 = x34**2 * x39
    x54 = x38 + x53
    x55 = x43 * x54
    x56 = x30 * x39
    x57 = x0 * (x40 + x56)
    x58 = x34 * x42
    x59 = x57 + x58
    x60 = x12 + x20 * x8
    x61 = x33 * x43
    x62 = x28 * x46
    x63 = x43 * x46**2
    x64 = x47 + x63
    x65 = x39 * x64
    x66 = x0 * (x48 + x61)
    x67 = x46 * x50
    x68 = x66 + x67
    x69 = 2 * x34 * x38
    x70 = x34 * x54 + x69
    x71 = x10 * x20 + x12
    x72 = 3 * x38
    x73 = x0 * (2 * x41 + x53 + x72)
    x74 = x34 * x59 + x73
    x75 = x27 * x5
    x76 = x74 * x75
    x77 = x20 * x75
    x78 = x50 * x7
    x79 = x64 * x7
    x80 = numpy.pi * x1 * x25 * x5
    x81 = x20 * x80
    x82 = 2 * x46 * x47
    x83 = x46 * x64 + x82
    x84 = 3 * x47
    x85 = x0 * (2 * x49 + x63 + x84)
    x86 = x46 * x68 + x85
    x87 = x80 * x86
    x88 = -x29 - A[1]
    x89 = x24 * x28
    x90 = x38 + x56 * x88
    x91 = x19 * x28
    x92 = x38 + x40 * x88
    x93 = x42 * x88 + x57
    x94 = x54 * x88 + x69
    x95 = x75 * (x59 * x88 + x73)
    x96 = x3 * x75
    x97 = x80 * x88
    x98 = x75 * (x0 * (3 * x53 + x72) + x70 * x88)
    x99 = x75 * x9
    x100 = -x32 - A[2]
    x101 = x100 * x61 + x47
    x102 = x100 * x48 + x47
    x103 = x100 * x50 + x66
    x104 = x3 * x80
    x105 = x100 * x64 + x82
    x106 = x80 * (x100 * x68 + x85)
    x107 = x80 * (x0 * (3 * x63 + x84) + x100 * x83)

    # 90 item(s)
    S = numpy.array(
        [
            x28 * (x0 * (3 * x11 + 3 * x15 + x19) + x20 * x24),
            x30 * x31,
            x31 * x33,
            x34 * x35,
            x36 * x44,
            x33 * x34 * x45,
            x35 * x46,
            x30 * x45 * x46,
            x36 * x51,
            x52 * x55,
            x43 * x59 * x60,
            x54 * x60 * x61,
            x34 * x52 * x62,
            x42 * x48 * x60,
            x40 * x50 * x60,
            x52 * x65,
            x56 * x60 * x64,
            x39 * x60 * x68,
            x43 * x70 * x71,
            x20 * x76,
            x33 * x70 * x77,
            x48 * x54 * x71,
            x46 * x59 * x77,
            x20 * x54 * x78,
            x40 * x64 * x71,
            x20 * x42 * x79,
            x34 * x68 * x81,
            x39 * x71 * x83,
            x30 * x81 * x83,
            x20 * x87,
            x88 * x89,
            x19 * x43 * x90,
            x33 * x88 * x91,
            x23 * x43 * x92,
            x18 * x43 * x93,
            x18 * x61 * x92,
            x23 * x62 * x88,
            x18 * x48 * x90,
            x18 * x51 * x88,
            x14 * x43 * x94,
            x3 * x95,
            x33 * x94 * x96,
            x14 * x48 * x92,
            x46 * x93 * x96,
            x50 * x8 * x92,
            x14 * x65 * x88,
            x64 * x8 * x90,
            x3 * x68 * x97,
            x9 * x98,
            x75 * (x0 * (3 * x57 + 3 * x58 + x70) + x74 * x88),
            x33 * x98,
            x46 * x94 * x99,
            x46 * x95,
            x78 * x94,
            x10 * x64 * x92,
            x79 * x93,
            x68 * x7 * x92,
            x83 * x9 * x97,
            x7 * x83 * x90,
            x87 * x88,
            x100 * x89,
            x100 * x30 * x91,
            x101 * x19 * x39,
            x100 * x23 * x28 * x34,
            x100 * x18 * x44,
            x101 * x18 * x40,
            x102 * x23 * x39,
            x102 * x18 * x56,
            x103 * x18 * x39,
            x100 * x14 * x55,
            x100 * x59 * x96,
            x101 * x54 * x8,
            x102 * x14 * x40,
            x102 * x42 * x8,
            x103 * x104 * x34,
            x105 * x14 * x39,
            x104 * x105 * x30,
            x106 * x3,
            x100 * x70 * x99,
            x100 * x76,
            x101 * x7 * x70,
            x10 * x102 * x54,
            x102 * x59 * x7,
            x103 * x54 * x7,
            x105 * x34 * x80 * x9,
            x105 * x42 * x7,
            x106 * x34,
            x107 * x9,
            x107 * x30,
            x80 * (x0 * (3 * x66 + 3 * x67 + x83) + x100 * x86),
        ]
    )
    return S


def dipole3d_14(a, A, b, B, C):
    """Cartesian 3D (pg) dipole moment integral.
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
    x10 = x5 * x9**2
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
    x37 = -x36 - C[1]
    x38 = x35 * (x0 * (8 * x22 + 4 * x25) + x27 * x28)
    x39 = -x1 * (a * A[2] + b * B[2])
    x40 = -x39 - C[2]
    x41 = -x36 - B[1]
    x42 = x35 * (x28 * x30 + x29)
    x43 = x21 + x26 * x28
    x44 = x0 * x4
    x45 = x32 * x44
    x46 = x32 * x4
    x47 = x41 * x46
    x48 = x37 * x47
    x49 = x45 + x48
    x50 = x33 * x4
    x51 = x49 * x50
    x52 = x35 * x43
    x53 = -x39 - B[2]
    x54 = x33 * x44
    x55 = x50 * x53
    x56 = x40 * x55
    x57 = x54 + x56
    x58 = x46 * x57
    x59 = x14 + x19 * x28
    x60 = x41**2 * x46
    x61 = x45 + x60
    x62 = x50 * x61
    x63 = x23 + x24 * x28
    x64 = x37 * x46
    x65 = x0 * (x47 + x64)
    x66 = x41 * x49
    x67 = x65 + x66
    x68 = x50 * x67
    x69 = x40 * x50
    x70 = x35 * x53
    x71 = x50 * x53**2
    x72 = x54 + x71
    x73 = x46 * x72
    x74 = x0 * (x55 + x69)
    x75 = x53 * x57
    x76 = x74 + x75
    x77 = x46 * x76
    x78 = x16 + x17 * x28
    x79 = x41 * x45
    x80 = 2 * x79
    x81 = x41 * x61
    x82 = x80 + x81
    x83 = x50 * x82
    x84 = 3 * x45
    x85 = x0 * (2 * x48 + x60 + x84)
    x86 = x41 * x67
    x87 = x85 + x86
    x88 = x12 * x28 + x6
    x89 = x53 * x54
    x90 = 2 * x89
    x91 = x53 * x72
    x92 = x90 + x91
    x93 = x46 * x92
    x94 = 3 * x54
    x95 = x0 * (2 * x56 + x71 + x94)
    x96 = x53 * x76
    x97 = x95 + x96
    x98 = x0 * (3 * x60 + x84)
    x99 = x41 * x82 + x98
    x100 = x15 * x28 + x6
    x101 = x0 * (3 * x65 + 3 * x66 + x82)
    x102 = x101 + x41 * x87
    x103 = x3 * x34
    x104 = x102 * x103
    x105 = x103 * x28
    x106 = x5 * x57
    x107 = x5 * x72
    x108 = x5 * x76
    x109 = x5 * x92
    x110 = numpy.pi * x1 * x3 * x32
    x111 = x110 * x28
    x112 = x0 * (3 * x71 + x94)
    x113 = x112 + x53 * x92
    x114 = x0 * (3 * x74 + 3 * x75 + x92)
    x115 = x114 + x53 * x97
    x116 = x110 * x115
    x117 = -x36 - A[1]
    x118 = x31 * x35
    x119 = x117 * x64 + x45
    x120 = x27 * x35
    x121 = x117 * x47 + x45
    x122 = x117 * x49 + x65
    x123 = x117 * x61 + x80
    x124 = x117 * x67 + x85
    x125 = x117 * x82 + x98
    x126 = x103 * (x101 + x117 * x87)
    x127 = x103 * x9
    x128 = x110 * x117
    x129 = x103 * (x0 * (8 * x79 + 4 * x81) + x117 * x99)
    x130 = x103 * x11
    x131 = -x39 - A[2]
    x132 = x131 * x69 + x54
    x133 = x131 * x55 + x54
    x134 = x131 * x57 + x74
    x135 = x131 * x72 + x90
    x136 = x131 * x76 + x95
    x137 = x110 * x9
    x138 = x112 + x131 * x92
    x139 = x110 * (x114 + x131 * x97)
    x140 = x110 * (x0 * (8 * x89 + 4 * x91) + x113 * x131)

    # 135 item(s)
    S = numpy.array(
        [
            x35 * (x0 * (4 * x14 + 4 * x20 + x27) + x28 * x31),
            x37 * x38,
            x38 * x40,
            x41 * x42,
            x43 * x51,
            x40 * x41 * x52,
            x42 * x53,
            x37 * x52 * x53,
            x43 * x58,
            x59 * x62,
            x63 * x68,
            x61 * x63 * x69,
            x41 * x59 * x70,
            x49 * x55 * x63,
            x47 * x57 * x63,
            x59 * x73,
            x63 * x64 * x72,
            x63 * x77,
            x78 * x83,
            x50 * x87 * x88,
            x69 * x82 * x88,
            x55 * x61 * x78,
            x55 * x67 * x88,
            x57 * x61 * x88,
            x47 * x72 * x78,
            x49 * x72 * x88,
            x47 * x76 * x88,
            x78 * x93,
            x64 * x88 * x92,
            x46 * x88 * x97,
            x100 * x50 * x99,
            x104 * x28,
            x105 * x40 * x99,
            x100 * x55 * x82,
            x105 * x53 * x87,
            x106 * x28 * x82,
            x100 * x61 * x72,
            x107 * x28 * x67,
            x108 * x28 * x61,
            x100 * x47 * x92,
            x109 * x28 * x49,
            x111 * x41 * x97,
            x100 * x113 * x46,
            x111 * x113 * x37,
            x116 * x28,
            x117 * x118,
            x119 * x27 * x50,
            x117 * x120 * x40,
            x121 * x30 * x50,
            x122 * x26 * x50,
            x121 * x26 * x69,
            x117 * x30 * x70,
            x119 * x26 * x55,
            x117 * x26 * x58,
            x123 * x19 * x50,
            x124 * x24 * x50,
            x123 * x24 * x69,
            x121 * x19 * x55,
            x122 * x24 * x55,
            x121 * x24 * x57,
            x117 * x19 * x73,
            x119 * x24 * x72,
            x117 * x24 * x77,
            x125 * x17 * x50,
            x126 * x9,
            x125 * x127 * x40,
            x123 * x17 * x55,
            x124 * x127 * x53,
            x12 * x123 * x57,
            x121 * x17 * x72,
            x12 * x122 * x72,
            x12 * x121 * x76,
            x117 * x17 * x93,
            x119 * x12 * x92,
            x128 * x9 * x97,
            x11 * x129,
            x103 * (x0 * (4 * x85 + 4 * x86 + x99) + x102 * x117),
            x129 * x40,
            x125 * x130 * x53,
            x126 * x53,
            x106 * x125,
            x123 * x15 * x72,
            x107 * x124,
            x108 * x123,
            x121 * x15 * x92,
            x109 * x122,
            x121 * x5 * x97,
            x11 * x113 * x128,
            x113 * x119 * x5,
            x116 * x117,
            x118 * x131,
            x120 * x131 * x37,
            x132 * x27 * x46,
            x131 * x30 * x35 * x41,
            x131 * x26 * x51,
            x132 * x26 * x47,
            x133 * x30 * x46,
            x133 * x26 * x64,
            x134 * x26 * x46,
            x131 * x19 * x62,
            x131 * x24 * x68,
            x132 * x24 * x61,
            x133 * x19 * x47,
            x133 * x24 * x49,
            x134 * x24 * x47,
            x135 * x19 * x46,
            x135 * x24 * x64,
            x136 * x24 * x46,
            x131 * x17 * x83,
            x127 * x131 * x87,
            x12 * x132 * x82,
            x133 * x17 * x61,
            x12 * x133 * x67,
            x12 * x134 * x61,
            x135 * x17 * x47,
            x12 * x135 * x49,
            x136 * x137 * x41,
            x138 * x17 * x46,
            x137 * x138 * x37,
            x139 * x9,
            x130 * x131 * x99,
            x104 * x131,
            x132 * x5 * x99,
            x133 * x15 * x82,
            x133 * x5 * x87,
            x134 * x5 * x82,
            x135 * x15 * x61,
            x135 * x5 * x67,
            x136 * x5 * x61,
            x11 * x110 * x138 * x41,
            x138 * x49 * x5,
            x139 * x41,
            x11 * x140,
            x140 * x37,
            x110 * (x0 * (x113 + 4 * x95 + 4 * x96) + x115 * x131),
        ]
    )
    return S


def dipole3d_20(a, A, b, B, C):
    """Cartesian 3D (ds) dipole moment integral.
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
    x19 = -x18 - C[1]
    x20 = x17 * (x11 * x5**2 + x12)
    x21 = -x1 * (a * A[2] + b * B[2])
    x22 = -x21 - C[2]
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


def dipole3d_21(a, A, b, B, C):
    """Cartesian 3D (dp) dipole moment integral.
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
    x23 = -x22 - C[1]
    x24 = x5 * x8
    x25 = x11 + x6
    x26 = x21 * (x0 * (x10 + x24) + x25 * x8)
    x27 = -x1 * (a * A[2] + b * B[2])
    x28 = -x27 - C[2]
    x29 = -x22 - B[1]
    x30 = x14 + x6
    x31 = x21 * (x0 * (x13 + x24) + x30 * x8)
    x32 = x0 * x4
    x33 = x18 * x32
    x34 = x18 * x4
    x35 = x29 * x34
    x36 = x23 * x35
    x37 = x33 + x36
    x38 = x5 * x8**2 + x6
    x39 = x19 * x4
    x40 = x21 * x38
    x41 = -x27 - B[2]
    x42 = x19 * x32
    x43 = x39 * x41
    x44 = x28 * x43
    x45 = x42 + x44
    x46 = -x22 - A[1]
    x47 = x17 * x21
    x48 = x23 * x34
    x49 = x46 * x48
    x50 = x33 + x49
    x51 = x21 * x46
    x52 = x35 * x46
    x53 = x33 + x52
    x54 = x0 * (x35 + x48) + x37 * x46
    x55 = x20 * x3
    x56 = x54 * x55
    x57 = x55 * x8
    x58 = x41 * x55
    x59 = numpy.pi * x1 * x18
    x60 = x3 * x59
    x61 = x60 * x8
    x62 = -x27 - A[2]
    x63 = x21 * x62
    x64 = x28 * x39
    x65 = x62 * x64
    x66 = x42 + x65
    x67 = x29 * x60
    x68 = x43 * x62
    x69 = x42 + x68
    x70 = x0 * (x43 + x64) + x45 * x62
    x71 = x60 * x70
    x72 = x33 + x34 * x46**2
    x73 = x34 * x46
    x74 = x0 * (x48 + x73) + x46 * x50
    x75 = x20 * x9
    x76 = x55 * (x0 * (x35 + x73) + x46 * x53)
    x77 = x59 * x9
    x78 = x39 * x62**2 + x42
    x79 = x39 * x62
    x80 = x0 * (x64 + x79) + x62 * x66
    x81 = x60 * (x0 * (x43 + x79) + x62 * x69)

    # 54 item(s)
    S = numpy.array(
        [
            x21 * (x0 * (x11 + x14 + x15 + 3 * x6) + x17 * x8),
            x23 * x26,
            x26 * x28,
            x29 * x31,
            x37 * x38 * x39,
            x28 * x29 * x40,
            x31 * x41,
            x23 * x40 * x41,
            x34 * x38 * x45,
            x46 * x47,
            x25 * x39 * x50,
            x25 * x28 * x51,
            x30 * x39 * x53,
            x56 * x8,
            x28 * x53 * x57,
            x30 * x41 * x51,
            x50 * x58 * x8,
            x45 * x46 * x61,
            x47 * x62,
            x23 * x25 * x63,
            x25 * x34 * x66,
            x29 * x30 * x63,
            x37 * x57 * x62,
            x66 * x67 * x8,
            x30 * x34 * x69,
            x23 * x61 * x69,
            x71 * x8,
            x16 * x39 * x72,
            x74 * x75,
            x28 * x72 * x75,
            x12 * x76,
            x55 * (x0 * (3 * x33 + x36 + x49 + x52) + x46 * x54),
            x28 * x76,
            x12 * x58 * x72,
            x58 * x74,
            x45 * x5 * x72,
            x16 * x51 * x62,
            x50 * x62 * x75,
            x46 * x66 * x77,
            x12 * x53 * x55 * x62,
            x56 * x62,
            x5 * x53 * x66,
            x12 * x46 * x60 * x69,
            x5 * x50 * x69,
            x46 * x71,
            x16 * x34 * x78,
            x23 * x77 * x78,
            x77 * x80,
            x12 * x67 * x78,
            x37 * x5 * x78,
            x67 * x80,
            x12 * x81,
            x23 * x81,
            x60 * (x0 * (3 * x42 + x44 + x65 + x68) + x62 * x70),
        ]
    )
    return S


def dipole3d_22(a, A, b, B, C):
    """Cartesian 3D (dd) dipole moment integral.
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
    x12 = x5**2 * x9
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
    x30 = -x29 - C[1]
    x31 = x10 * x2
    x32 = x28 * (x2 * x24 + x3 * (x15 + 2 * x31))
    x33 = -x0 * (a * A[2] + b * B[2])
    x34 = -x33 - C[2]
    x35 = -x29 - B[1]
    x36 = x16 * x2
    x37 = x17 + x22
    x38 = x28 * (x2 * x37 + x3 * (x11 + x14 + x31 + x36))
    x39 = x3 * x8
    x40 = x25 * x39
    x41 = x25 * x8
    x42 = x35 * x41
    x43 = x30 * x42
    x44 = x40 + x43
    x45 = x2 * x9
    x46 = x13 + x31
    x47 = x2 * x46 + x3 * (x10 + x45)
    x48 = x26 * x8
    x49 = x28 * x47
    x50 = -x33 - B[2]
    x51 = x26 * x39
    x52 = x48 * x50
    x53 = x34 * x52
    x54 = x51 + x53
    x55 = x35**2 * x41
    x56 = x40 + x55
    x57 = x13 + x36
    x58 = x2 * x57 + x3 * (x16 + x45)
    x59 = x30 * x41
    x60 = x3 * (x42 + x59)
    x61 = x35 * x44
    x62 = x60 + x61
    x63 = x13 + x2**2 * x9
    x64 = x34 * x48
    x65 = x28 * x50
    x66 = x48 * x50**2
    x67 = x51 + x66
    x68 = x3 * (x52 + x64)
    x69 = x50 * x54
    x70 = x68 + x69
    x71 = -x29 - A[1]
    x72 = x21 * x28
    x73 = x59 * x71
    x74 = x40 + x73
    x75 = x24 * x28
    x76 = x42 * x71
    x77 = x40 + x76
    x78 = x44 * x71
    x79 = x60 + x78
    x80 = x41 * x71
    x81 = 2 * x35 * x40 + x56 * x71
    x82 = 3 * x40
    x83 = x55 + x82
    x84 = x3 * (2 * x43 + x83) + x62 * x71
    x85 = x27 * x7
    x86 = x84 * x85
    x87 = x2 * x85
    x88 = numpy.pi * x0 * x25 * x7
    x89 = x2 * x88
    x90 = -x33 - A[2]
    x91 = x64 * x90
    x92 = x51 + x91
    x93 = x28 * x90
    x94 = x48 * x90
    x95 = x52 * x90
    x96 = x51 + x95
    x97 = x54 * x90
    x98 = x68 + x97
    x99 = 2 * x50 * x51 + x67 * x90
    x100 = 3 * x51
    x101 = x100 + x66
    x102 = x3 * (x101 + 2 * x53) + x70 * x90
    x103 = x102 * x88
    x104 = x40 + x41 * x71**2
    x105 = x3 * (x59 + x80) + x71 * x74
    x106 = x3 * (x42 + x80) + x71 * x77
    x107 = x85 * (x3 * (x43 + x73 + x76 + x82) + x71 * x79)
    x108 = x5 * x85
    x109 = x85 * (x3 * (2 * x76 + x83) + x71 * x81)
    x110 = x4 * x85
    x111 = x71 * x88
    x112 = x48 * x90**2 + x51
    x113 = x3 * (x64 + x94) + x90 * x92
    x114 = x5 * x88
    x115 = x3 * (x52 + x94) + x90 * x96
    x116 = x88 * (x3 * (x100 + x53 + x91 + x95) + x90 * x98)
    x117 = x88 * (x3 * (x101 + 2 * x95) + x90 * x99)

    # 108 item(s)
    S = numpy.array(
        [
            x28 * (x2 * x21 + x3 * (3 * x17 + x19 + 2 * x22 + x24)),
            x30 * x32,
            x32 * x34,
            x35 * x38,
            x44 * x47 * x48,
            x34 * x35 * x49,
            x38 * x50,
            x30 * x49 * x50,
            x41 * x47 * x54,
            x48 * x56 * x58,
            x48 * x62 * x63,
            x56 * x63 * x64,
            x35 * x58 * x65,
            x44 * x52 * x63,
            x42 * x54 * x63,
            x41 * x58 * x67,
            x59 * x63 * x67,
            x41 * x63 * x70,
            x71 * x72,
            x24 * x48 * x74,
            x34 * x71 * x75,
            x37 * x48 * x77,
            x46 * x48 * x79,
            x46 * x64 * x77,
            x37 * x65 * x71,
            x46 * x52 * x74,
            x46 * x54 * x80,
            x48 * x57 * x81,
            x2 * x86,
            x34 * x81 * x87,
            x52 * x57 * x77,
            x50 * x79 * x87,
            x45 * x54 * x77,
            x57 * x67 * x80,
            x45 * x67 * x74,
            x70 * x71 * x89,
            x72 * x90,
            x30 * x75 * x90,
            x24 * x41 * x92,
            x35 * x37 * x93,
            x44 * x46 * x94,
            x42 * x46 * x92,
            x37 * x41 * x96,
            x46 * x59 * x96,
            x41 * x46 * x98,
            x56 * x57 * x94,
            x62 * x87 * x90,
            x45 * x56 * x92,
            x42 * x57 * x96,
            x44 * x45 * x96,
            x35 * x89 * x98,
            x41 * x57 * x99,
            x30 * x89 * x99,
            x103 * x2,
            x104 * x20 * x48,
            x105 * x23 * x48,
            x104 * x23 * x64,
            x106 * x18 * x48,
            x107 * x5,
            x106 * x108 * x34,
            x104 * x18 * x52,
            x105 * x108 * x50,
            x10 * x104 * x54,
            x109 * x4,
            x85 * (x3 * (3 * x60 + x61 + 2 * x78 + x81) + x71 * x84),
            x109 * x34,
            x106 * x110 * x50,
            x107 * x50,
            x106 * x54 * x9,
            x104 * x16 * x67,
            x105 * x67 * x9,
            x104 * x70 * x9,
            x20 * x71 * x93,
            x23 * x74 * x94,
            x23 * x80 * x92,
            x18 * x77 * x94,
            x108 * x79 * x90,
            x10 * x77 * x92,
            x18 * x80 * x96,
            x10 * x74 * x96,
            x111 * x5 * x98,
            x110 * x81 * x90,
            x86 * x90,
            x81 * x9 * x92,
            x16 * x77 * x96,
            x79 * x9 * x96,
            x77 * x9 * x98,
            x111 * x4 * x99,
            x74 * x9 * x99,
            x103 * x71,
            x112 * x20 * x41,
            x112 * x23 * x59,
            x113 * x23 * x41,
            x112 * x18 * x42,
            x10 * x112 * x44,
            x113 * x114 * x35,
            x115 * x18 * x41,
            x114 * x115 * x30,
            x116 * x5,
            x112 * x16 * x56,
            x112 * x62 * x9,
            x113 * x56 * x9,
            x115 * x35 * x4 * x88,
            x115 * x44 * x9,
            x116 * x35,
            x117 * x4,
            x117 * x30,
            x88 * (x102 * x90 + x3 * (3 * x68 + x69 + 2 * x97 + x99)),
        ]
    )
    return S


def dipole3d_23(a, A, b, B, C):
    """Cartesian 3D (df) dipole moment integral.
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
    x20 = x4**2 * x8
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
    x38 = -x37 - C[1]
    x39 = x2 * x21
    x40 = x36 * (x2 * x32 + x3 * (8 * x18 + x22 + 3 * x39))
    x41 = -x0 * (a * A[2] + b * B[2])
    x42 = -x41 - C[2]
    x43 = -x37 - B[1]
    x44 = x26 + x31
    x45 = x16 * x2
    x46 = x19 + x39
    x47 = x36 * (x2 * x44 + x3 * (x13 + x17 + 2 * x45 + x46))
    x48 = x3 * x7
    x49 = x33 * x48
    x50 = x33 * x7
    x51 = x43 * x50
    x52 = x38 * x51
    x53 = x49 + x52
    x54 = x2 * x9
    x55 = x2 * x46 + x3 * (x25 + 2 * x54)
    x56 = x34 * x7
    x57 = x36 * x55
    x58 = -x41 - B[2]
    x59 = x34 * x48
    x60 = x56 * x58
    x61 = x42 * x60
    x62 = x59 + x61
    x63 = x43**2 * x50
    x64 = x49 + x63
    x65 = x11 * x2
    x66 = x12 + x45
    x67 = x2 * x66 + x3 * (x15 + x24 + x54 + x65)
    x68 = x38 * x50
    x69 = x3 * (x51 + x68)
    x70 = x43 * x53
    x71 = x69 + x70
    x72 = x2 * x8
    x73 = x14 + x54
    x74 = x2 * x73 + x3 * (x72 + x9)
    x75 = x42 * x56
    x76 = x36 * x58
    x77 = x56 * x58**2
    x78 = x59 + x77
    x79 = x3 * (x60 + x75)
    x80 = x58 * x62
    x81 = x79 + x80
    x82 = x43 * x49
    x83 = 2 * x82
    x84 = x43 * x64
    x85 = x83 + x84
    x86 = x14 + x65
    x87 = x2 * x86 + x3 * (x11 + x72)
    x88 = 3 * x49
    x89 = x63 + x88
    x90 = x3 * (2 * x52 + x89)
    x91 = x43 * x71
    x92 = x90 + x91
    x93 = x14 + x2**2 * x8
    x94 = x58 * x59
    x95 = 2 * x94
    x96 = x58 * x78
    x97 = x95 + x96
    x98 = 3 * x59
    x99 = x77 + x98
    x100 = x3 * (2 * x61 + x99)
    x101 = x58 * x81
    x102 = x100 + x101
    x103 = -x37 - A[1]
    x104 = x30 * x36
    x105 = x103 * x68
    x106 = x105 + x49
    x107 = x32 * x36
    x108 = x103 * x51
    x109 = x108 + x49
    x110 = x103 * x53
    x111 = x110 + x69
    x112 = x103 * x50
    x113 = x103 * x64
    x114 = x113 + x83
    x115 = x103 * x71
    x116 = x115 + x90
    x117 = x103 * x85 + x3 * (3 * x63 + x88)
    x118 = 3 * x69
    x119 = x103 * x92 + x3 * (x118 + 3 * x70 + x85)
    x120 = x35 * x6
    x121 = x119 * x120
    x122 = x120 * x2
    x123 = numpy.pi * x0 * x33 * x6
    x124 = x123 * x2
    x125 = -x41 - A[2]
    x126 = x125 * x75
    x127 = x126 + x59
    x128 = x125 * x36
    x129 = x125 * x56
    x130 = x125 * x60
    x131 = x130 + x59
    x132 = x125 * x62
    x133 = x132 + x79
    x134 = x125 * x78
    x135 = x134 + x95
    x136 = x125 * x81
    x137 = x100 + x136
    x138 = x125 * x97 + x3 * (3 * x77 + x98)
    x139 = 3 * x79
    x140 = x102 * x125 + x3 * (x139 + 3 * x80 + x97)
    x141 = x123 * x140
    x142 = x103**2 * x50 + x49
    x143 = x103 * x106 + x3 * (x112 + x68)
    x144 = x103 * x109 + x3 * (x112 + x51)
    x145 = x103 * x111 + x3 * (x105 + x108 + x52 + x88)
    x146 = x103 * x114 + x3 * (2 * x108 + x89)
    x147 = x120 * (x103 * x116 + x3 * (2 * x110 + x114 + x118 + x70))
    x148 = x120 * x4
    x149 = x120 * (x103 * x117 + x3 * (3 * x113 + 8 * x82 + x84))
    x150 = x10 * x120
    x151 = x103 * x123
    x152 = x125**2 * x56 + x59
    x153 = x125 * x127 + x3 * (x129 + x75)
    x154 = x125 * x131 + x3 * (x129 + x60)
    x155 = x125 * x133 + x3 * (x126 + x130 + x61 + x98)
    x156 = x123 * x4
    x157 = x125 * x135 + x3 * (2 * x130 + x99)
    x158 = x123 * (x125 * x137 + x3 * (2 * x132 + x135 + x139 + x80))
    x159 = x123 * (x125 * x138 + x3 * (3 * x134 + 8 * x94 + x96))

    # 180 item(s)
    S = numpy.array(
        [
            x36 * (x2 * x30 + x3 * (4 * x26 + x28 + 3 * x31 + x32)),
            x38 * x40,
            x40 * x42,
            x43 * x47,
            x53 * x55 * x56,
            x42 * x43 * x57,
            x47 * x58,
            x38 * x57 * x58,
            x50 * x55 * x62,
            x56 * x64 * x67,
            x56 * x71 * x74,
            x64 * x74 * x75,
            x43 * x67 * x76,
            x53 * x60 * x74,
            x51 * x62 * x74,
            x50 * x67 * x78,
            x68 * x74 * x78,
            x50 * x74 * x81,
            x56 * x85 * x87,
            x56 * x92 * x93,
            x75 * x85 * x93,
            x60 * x64 * x87,
            x60 * x71 * x93,
            x62 * x64 * x93,
            x51 * x78 * x87,
            x53 * x78 * x93,
            x51 * x81 * x93,
            x50 * x87 * x97,
            x68 * x93 * x97,
            x102 * x50 * x93,
            x103 * x104,
            x106 * x32 * x56,
            x103 * x107 * x42,
            x109 * x44 * x56,
            x111 * x46 * x56,
            x109 * x46 * x75,
            x103 * x44 * x76,
            x106 * x46 * x60,
            x112 * x46 * x62,
            x114 * x56 * x66,
            x116 * x56 * x73,
            x114 * x73 * x75,
            x109 * x60 * x66,
            x111 * x60 * x73,
            x109 * x62 * x73,
            x112 * x66 * x78,
            x106 * x73 * x78,
            x112 * x73 * x81,
            x117 * x56 * x86,
            x121 * x2,
            x117 * x122 * x42,
            x114 * x60 * x86,
            x116 * x122 * x58,
            x114 * x62 * x72,
            x109 * x78 * x86,
            x111 * x72 * x78,
            x109 * x72 * x81,
            x112 * x86 * x97,
            x106 * x72 * x97,
            x102 * x103 * x124,
            x104 * x125,
            x107 * x125 * x38,
            x127 * x32 * x50,
            x128 * x43 * x44,
            x129 * x46 * x53,
            x127 * x46 * x51,
            x131 * x44 * x50,
            x131 * x46 * x68,
            x133 * x46 * x50,
            x129 * x64 * x66,
            x129 * x71 * x73,
            x127 * x64 * x73,
            x131 * x51 * x66,
            x131 * x53 * x73,
            x133 * x51 * x73,
            x135 * x50 * x66,
            x135 * x68 * x73,
            x137 * x50 * x73,
            x129 * x85 * x86,
            x122 * x125 * x92,
            x127 * x72 * x85,
            x131 * x64 * x86,
            x131 * x71 * x72,
            x133 * x64 * x72,
            x135 * x51 * x86,
            x135 * x53 * x72,
            x124 * x137 * x43,
            x138 * x50 * x86,
            x124 * x138 * x38,
            x141 * x2,
            x142 * x29 * x56,
            x143 * x23 * x56,
            x142 * x23 * x75,
            x144 * x27 * x56,
            x145 * x21 * x56,
            x144 * x21 * x75,
            x142 * x27 * x60,
            x143 * x21 * x60,
            x142 * x21 * x62,
            x146 * x16 * x56,
            x147 * x4,
            x146 * x148 * x42,
            x144 * x16 * x60,
            x145 * x148 * x58,
            x144 * x62 * x9,
            x142 * x16 * x78,
            x143 * x78 * x9,
            x142 * x81 * x9,
            x10 * x149,
            x120 * (x103 * x119 + x3 * (3 * x115 + x117 + 4 * x90 + x91)),
            x149 * x42,
            x146 * x150 * x58,
            x147 * x58,
            x146 * x62 * x8,
            x11 * x144 * x78,
            x145 * x78 * x8,
            x144 * x8 * x81,
            x11 * x142 * x97,
            x143 * x8 * x97,
            x102 * x142 * x8,
            x103 * x128 * x29,
            x106 * x129 * x23,
            x112 * x127 * x23,
            x109 * x129 * x27,
            x111 * x129 * x21,
            x109 * x127 * x21,
            x112 * x131 * x27,
            x106 * x131 * x21,
            x112 * x133 * x21,
            x114 * x129 * x16,
            x116 * x125 * x148,
            x114 * x127 * x9,
            x109 * x131 * x16,
            x111 * x131 * x9,
            x109 * x133 * x9,
            x112 * x135 * x16,
            x106 * x135 * x9,
            x137 * x151 * x4,
            x117 * x125 * x150,
            x121 * x125,
            x117 * x127 * x8,
            x11 * x114 * x131,
            x116 * x131 * x8,
            x114 * x133 * x8,
            x109 * x11 * x135,
            x111 * x135 * x8,
            x109 * x137 * x8,
            x10 * x138 * x151,
            x106 * x138 * x8,
            x103 * x141,
            x152 * x29 * x50,
            x152 * x23 * x68,
            x153 * x23 * x50,
            x152 * x27 * x51,
            x152 * x21 * x53,
            x153 * x21 * x51,
            x154 * x27 * x50,
            x154 * x21 * x68,
            x155 * x21 * x50,
            x152 * x16 * x64,
            x152 * x71 * x9,
            x153 * x64 * x9,
            x154 * x16 * x51,
            x154 * x53 * x9,
            x155 * x156 * x43,
            x157 * x16 * x50,
            x156 * x157 * x38,
            x158 * x4,
            x11 * x152 * x85,
            x152 * x8 * x92,
            x153 * x8 * x85,
            x11 * x154 * x64,
            x154 * x71 * x8,
            x155 * x64 * x8,
            x10 * x123 * x157 * x43,
            x157 * x53 * x8,
            x158 * x43,
            x10 * x159,
            x159 * x38,
            x123 * (x125 * x140 + x3 * (4 * x100 + x101 + 3 * x136 + x138)),
        ]
    )
    return S


def dipole3d_24(a, A, b, B, C):
    """Cartesian 3D (dg) dipole moment integral.
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
    x12 = x5**2 * x9
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
    x46 = -x45 - C[1]
    x47 = x2 * x29
    x48 = x44 * (x2 * x40 + x3 * (5 * x24 + x30 + 4 * x47))
    x49 = -x0 * (a * A[2] + b * B[2])
    x50 = -x49 - C[2]
    x51 = -x45 - B[1]
    x52 = x33 + x38
    x53 = x2 * x22
    x54 = x24 + x47
    x55 = x44 * (x2 * x52 + x3 * (x17 + x23 + 3 * x53 + x54))
    x56 = x3 * x8
    x57 = x41 * x56
    x58 = x41 * x8
    x59 = x51 * x58
    x60 = x46 * x59
    x61 = x57 + x60
    x62 = x2 * x27
    x63 = x2 * x54 + x3 * (x28 + x39 + 3 * x62)
    x64 = x42 * x8
    x65 = x44 * x63
    x66 = -x49 - B[2]
    x67 = x42 * x56
    x68 = x64 * x66
    x69 = x50 * x68
    x70 = x67 + x69
    x71 = x51**2 * x58
    x72 = x57 + x71
    x73 = x16 + x53
    x74 = x2 * x20
    x75 = x26 + x62
    x76 = x2 * x73 + x3 * (x21 + x32 + 2 * x74 + x75)
    x77 = x46 * x58
    x78 = x3 * (x59 + x77)
    x79 = x51 * x61
    x80 = x78 + x79
    x81 = x10 * x2
    x82 = x2 * x75 + x3 * (x15 + 2 * x81)
    x83 = x50 * x64
    x84 = x44 * x66
    x85 = x64 * x66**2
    x86 = x67 + x85
    x87 = x3 * (x68 + x83)
    x88 = x66 * x70
    x89 = x87 + x88
    x90 = x51 * x57
    x91 = 2 * x90
    x92 = x51 * x72
    x93 = x91 + x92
    x94 = x18 * x2
    x95 = x19 + x74
    x96 = x2 * x95 + x3 * (x11 + x14 + x81 + x94)
    x97 = 3 * x57
    x98 = x71 + x97
    x99 = x3 * (2 * x60 + x98)
    x100 = x51 * x80
    x101 = x100 + x99
    x102 = x2 * x9
    x103 = x13 + x81
    x104 = x103 * x2 + x3 * (x10 + x102)
    x105 = x66 * x67
    x106 = 2 * x105
    x107 = x66 * x86
    x108 = x106 + x107
    x109 = 3 * x67
    x110 = x109 + x85
    x111 = x3 * (x110 + 2 * x69)
    x112 = x66 * x89
    x113 = x111 + x112
    x114 = x3 * (3 * x71 + x97)
    x115 = x51 * x93
    x116 = x114 + x115
    x117 = x13 + x94
    x118 = x117 * x2 + x3 * (x102 + x18)
    x119 = 3 * x78
    x120 = x3 * (x119 + 3 * x79 + x93)
    x121 = x101 * x51
    x122 = x120 + x121
    x123 = x13 + x2**2 * x9
    x124 = x3 * (x109 + 3 * x85)
    x125 = x108 * x66
    x126 = x124 + x125
    x127 = 3 * x87
    x128 = x3 * (x108 + x127 + 3 * x88)
    x129 = x113 * x66
    x130 = x128 + x129
    x131 = -x45 - A[1]
    x132 = x37 * x44
    x133 = x131 * x77
    x134 = x133 + x57
    x135 = x40 * x44
    x136 = x131 * x59
    x137 = x136 + x57
    x138 = x131 * x61
    x139 = x138 + x78
    x140 = x131 * x58
    x141 = x131 * x72
    x142 = x141 + x91
    x143 = x131 * x80
    x144 = x143 + x99
    x145 = x131 * x93
    x146 = x114 + x145
    x147 = x101 * x131
    x148 = x120 + x147
    x149 = 8 * x90
    x150 = x116 * x131 + x3 * (x149 + 4 * x92)
    x151 = 4 * x99
    x152 = x122 * x131 + x3 * (4 * x100 + x116 + x151)
    x153 = x43 * x7
    x154 = x152 * x153
    x155 = x153 * x2
    x156 = numpy.pi * x0 * x41 * x7
    x157 = x156 * x2
    x158 = -x49 - A[2]
    x159 = x158 * x83
    x160 = x159 + x67
    x161 = x158 * x44
    x162 = x158 * x64
    x163 = x158 * x68
    x164 = x163 + x67
    x165 = x158 * x70
    x166 = x165 + x87
    x167 = x158 * x86
    x168 = x106 + x167
    x169 = x158 * x89
    x170 = x111 + x169
    x171 = x108 * x158
    x172 = x124 + x171
    x173 = x113 * x158
    x174 = x128 + x173
    x175 = 8 * x105
    x176 = x126 * x158 + x3 * (4 * x107 + x175)
    x177 = 4 * x111
    x178 = x130 * x158 + x3 * (4 * x112 + x126 + x177)
    x179 = x156 * x178
    x180 = x131**2 * x58 + x57
    x181 = x131 * x134 + x3 * (x140 + x77)
    x182 = x131 * x137 + x3 * (x140 + x59)
    x183 = x131 * x139 + x3 * (x133 + x136 + x60 + x97)
    x184 = x131 * x142 + x3 * (2 * x136 + x98)
    x185 = x131 * x144 + x3 * (x119 + 2 * x138 + x142 + x79)
    x186 = x131 * x146 + x3 * (3 * x141 + x149 + x92)
    x187 = x153 * (x131 * x148 + x3 * (x100 + 3 * x143 + x146 + x151))
    x188 = x153 * x5
    x189 = x153 * (x131 * x150 + x3 * (5 * x114 + x115 + 4 * x145))
    x190 = x153 * x4
    x191 = x131 * x156
    x192 = x158**2 * x64 + x67
    x193 = x158 * x160 + x3 * (x162 + x83)
    x194 = x158 * x164 + x3 * (x162 + x68)
    x195 = x158 * x166 + x3 * (x109 + x159 + x163 + x69)
    x196 = x158 * x168 + x3 * (x110 + 2 * x163)
    x197 = x158 * x170 + x3 * (x127 + 2 * x165 + x168 + x88)
    x198 = x156 * x5
    x199 = x158 * x172 + x3 * (x107 + 3 * x167 + x175)
    x200 = x156 * (x158 * x174 + x3 * (x112 + 3 * x169 + x172 + x177))
    x201 = x156 * (x158 * x176 + x3 * (5 * x124 + x125 + 4 * x171))

    # 270 item(s)
    S = numpy.array(
        [
            x44 * (x2 * x37 + x3 * (5 * x33 + x35 + 4 * x38 + x40)),
            x46 * x48,
            x48 * x50,
            x51 * x55,
            x61 * x63 * x64,
            x50 * x51 * x65,
            x55 * x66,
            x46 * x65 * x66,
            x58 * x63 * x70,
            x64 * x72 * x76,
            x64 * x80 * x82,
            x72 * x82 * x83,
            x51 * x76 * x84,
            x61 * x68 * x82,
            x59 * x70 * x82,
            x58 * x76 * x86,
            x77 * x82 * x86,
            x58 * x82 * x89,
            x64 * x93 * x96,
            x101 * x104 * x64,
            x104 * x83 * x93,
            x68 * x72 * x96,
            x104 * x68 * x80,
            x104 * x70 * x72,
            x59 * x86 * x96,
            x104 * x61 * x86,
            x104 * x59 * x89,
            x108 * x58 * x96,
            x104 * x108 * x77,
            x104 * x113 * x58,
            x116 * x118 * x64,
            x122 * x123 * x64,
            x116 * x123 * x83,
            x118 * x68 * x93,
            x101 * x123 * x68,
            x123 * x70 * x93,
            x118 * x72 * x86,
            x123 * x80 * x86,
            x123 * x72 * x89,
            x108 * x118 * x59,
            x108 * x123 * x61,
            x113 * x123 * x59,
            x118 * x126 * x58,
            x123 * x126 * x77,
            x123 * x130 * x58,
            x131 * x132,
            x134 * x40 * x64,
            x131 * x135 * x50,
            x137 * x52 * x64,
            x139 * x54 * x64,
            x137 * x54 * x83,
            x131 * x52 * x84,
            x134 * x54 * x68,
            x140 * x54 * x70,
            x142 * x64 * x73,
            x144 * x64 * x75,
            x142 * x75 * x83,
            x137 * x68 * x73,
            x139 * x68 * x75,
            x137 * x70 * x75,
            x140 * x73 * x86,
            x134 * x75 * x86,
            x140 * x75 * x89,
            x146 * x64 * x95,
            x103 * x148 * x64,
            x103 * x146 * x83,
            x142 * x68 * x95,
            x103 * x144 * x68,
            x103 * x142 * x70,
            x137 * x86 * x95,
            x103 * x139 * x86,
            x103 * x137 * x89,
            x108 * x140 * x95,
            x103 * x108 * x134,
            x103 * x113 * x140,
            x117 * x150 * x64,
            x154 * x2,
            x150 * x155 * x50,
            x117 * x146 * x68,
            x148 * x155 * x66,
            x102 * x146 * x70,
            x117 * x142 * x86,
            x102 * x144 * x86,
            x102 * x142 * x89,
            x108 * x117 * x137,
            x102 * x108 * x139,
            x102 * x113 * x137,
            x117 * x126 * x140,
            x102 * x126 * x134,
            x130 * x131 * x157,
            x132 * x158,
            x135 * x158 * x46,
            x160 * x40 * x58,
            x161 * x51 * x52,
            x162 * x54 * x61,
            x160 * x54 * x59,
            x164 * x52 * x58,
            x164 * x54 * x77,
            x166 * x54 * x58,
            x162 * x72 * x73,
            x162 * x75 * x80,
            x160 * x72 * x75,
            x164 * x59 * x73,
            x164 * x61 * x75,
            x166 * x59 * x75,
            x168 * x58 * x73,
            x168 * x75 * x77,
            x170 * x58 * x75,
            x162 * x93 * x95,
            x101 * x103 * x162,
            x103 * x160 * x93,
            x164 * x72 * x95,
            x103 * x164 * x80,
            x103 * x166 * x72,
            x168 * x59 * x95,
            x103 * x168 * x61,
            x103 * x170 * x59,
            x172 * x58 * x95,
            x103 * x172 * x77,
            x103 * x174 * x58,
            x116 * x117 * x162,
            x122 * x155 * x158,
            x102 * x116 * x160,
            x117 * x164 * x93,
            x101 * x102 * x164,
            x102 * x166 * x93,
            x117 * x168 * x72,
            x102 * x168 * x80,
            x102 * x170 * x72,
            x117 * x172 * x59,
            x102 * x172 * x61,
            x157 * x174 * x51,
            x117 * x176 * x58,
            x157 * x176 * x46,
            x179 * x2,
            x180 * x36 * x64,
            x181 * x31 * x64,
            x180 * x31 * x83,
            x182 * x34 * x64,
            x183 * x29 * x64,
            x182 * x29 * x83,
            x180 * x34 * x68,
            x181 * x29 * x68,
            x180 * x29 * x70,
            x184 * x22 * x64,
            x185 * x27 * x64,
            x184 * x27 * x83,
            x182 * x22 * x68,
            x183 * x27 * x68,
            x182 * x27 * x70,
            x180 * x22 * x86,
            x181 * x27 * x86,
            x180 * x27 * x89,
            x186 * x20 * x64,
            x187 * x5,
            x186 * x188 * x50,
            x184 * x20 * x68,
            x185 * x188 * x66,
            x10 * x184 * x70,
            x182 * x20 * x86,
            x10 * x183 * x86,
            x10 * x182 * x89,
            x108 * x180 * x20,
            x10 * x108 * x181,
            x10 * x113 * x180,
            x189 * x4,
            x153 * (x131 * x152 + x3 * (5 * x120 + x121 + 4 * x147 + x150)),
            x189 * x50,
            x186 * x190 * x66,
            x187 * x66,
            x186 * x70 * x9,
            x18 * x184 * x86,
            x185 * x86 * x9,
            x184 * x89 * x9,
            x108 * x18 * x182,
            x108 * x183 * x9,
            x113 * x182 * x9,
            x126 * x18 * x180,
            x126 * x181 * x9,
            x130 * x180 * x9,
            x131 * x161 * x36,
            x134 * x162 * x31,
            x140 * x160 * x31,
            x137 * x162 * x34,
            x139 * x162 * x29,
            x137 * x160 * x29,
            x140 * x164 * x34,
            x134 * x164 * x29,
            x140 * x166 * x29,
            x142 * x162 * x22,
            x144 * x162 * x27,
            x142 * x160 * x27,
            x137 * x164 * x22,
            x139 * x164 * x27,
            x137 * x166 * x27,
            x140 * x168 * x22,
            x134 * x168 * x27,
            x140 * x170 * x27,
            x146 * x162 * x20,
            x148 * x158 * x188,
            x10 * x146 * x160,
            x142 * x164 * x20,
            x10 * x144 * x164,
            x10 * x142 * x166,
            x137 * x168 * x20,
            x10 * x139 * x168,
            x10 * x137 * x170,
            x140 * x172 * x20,
            x10 * x134 * x172,
            x174 * x191 * x5,
            x150 * x158 * x190,
            x154 * x158,
            x150 * x160 * x9,
            x146 * x164 * x18,
            x148 * x164 * x9,
            x146 * x166 * x9,
            x142 * x168 * x18,
            x144 * x168 * x9,
            x142 * x170 * x9,
            x137 * x172 * x18,
            x139 * x172 * x9,
            x137 * x174 * x9,
            x176 * x191 * x4,
            x134 * x176 * x9,
            x131 * x179,
            x192 * x36 * x58,
            x192 * x31 * x77,
            x193 * x31 * x58,
            x192 * x34 * x59,
            x192 * x29 * x61,
            x193 * x29 * x59,
            x194 * x34 * x58,
            x194 * x29 * x77,
            x195 * x29 * x58,
            x192 * x22 * x72,
            x192 * x27 * x80,
            x193 * x27 * x72,
            x194 * x22 * x59,
            x194 * x27 * x61,
            x195 * x27 * x59,
            x196 * x22 * x58,
            x196 * x27 * x77,
            x197 * x27 * x58,
            x192 * x20 * x93,
            x10 * x101 * x192,
            x10 * x193 * x93,
            x194 * x20 * x72,
            x10 * x194 * x80,
            x10 * x195 * x72,
            x196 * x20 * x59,
            x10 * x196 * x61,
            x197 * x198 * x51,
            x199 * x20 * x58,
            x198 * x199 * x46,
            x200 * x5,
            x116 * x18 * x192,
            x122 * x192 * x9,
            x116 * x193 * x9,
            x18 * x194 * x93,
            x101 * x194 * x9,
            x195 * x9 * x93,
            x18 * x196 * x72,
            x196 * x80 * x9,
            x197 * x72 * x9,
            x156 * x199 * x4 * x51,
            x199 * x61 * x9,
            x200 * x51,
            x201 * x4,
            x201 * x46,
            x156 * (x158 * x178 + x3 * (5 * x128 + x129 + 4 * x173 + x176)),
        ]
    )
    return S


def dipole3d_30(a, A, b, B, C):
    """Cartesian 3D (fs) dipole moment integral.
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
    x9 = x5 * x8**2
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
    x20 = -x19 - C[1]
    x21 = x6 + x9
    x22 = x18 * (2 * x0 * x11 + x21 * x8)
    x23 = -x1 * (a * A[2] + b * B[2])
    x24 = -x23 - C[2]
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


def dipole3d_31(a, A, b, B, C):
    """Cartesian 3D (fp) dipole moment integral.
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
    x24 = x13 + x8
    x25 = x2 * x24 + x3 * (x12 + x16)
    x26 = x15 + x8
    x27 = x2 * x26 + x3 * (x12 + x18)
    x28 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x29 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x30 = numpy.pi * x0 * x29
    x31 = x28 * x30
    x32 = -x0 * (a * A[1] + b * B[1])
    x33 = -x32 - C[1]
    x34 = x2**2 * x7
    x35 = x34 + x9
    x36 = x31 * (x2 * x25 + x3 * (2 * x13 + x35))
    x37 = -x0 * (a * A[2] + b * B[2])
    x38 = -x37 - C[2]
    x39 = -x32 - B[1]
    x40 = x31 * (x2 * x27 + x3 * (2 * x15 + x35))
    x41 = x3 * x6
    x42 = x28 * x41
    x43 = x28 * x6
    x44 = x39 * x43
    x45 = x33 * x44
    x46 = x42 + x45
    x47 = x34 + x8
    x48 = 2 * x12 * x3 + x2 * x47
    x49 = x29 * x6
    x50 = x31 * x48
    x51 = -x37 - B[2]
    x52 = x29 * x41
    x53 = x49 * x51
    x54 = x38 * x53
    x55 = x52 + x54
    x56 = -x32 - A[1]
    x57 = x23 * x31
    x58 = x33 * x43
    x59 = x56 * x58
    x60 = x42 + x59
    x61 = x31 * x56
    x62 = x44 * x56
    x63 = x42 + x62
    x64 = x3 * (x44 + x58)
    x65 = x46 * x56
    x66 = x64 + x65
    x67 = x38 * x49
    x68 = x43 * x56
    x69 = -x37 - A[2]
    x70 = x31 * x69
    x71 = x67 * x69
    x72 = x52 + x71
    x73 = x49 * x69
    x74 = x53 * x69
    x75 = x52 + x74
    x76 = x3 * (x53 + x67)
    x77 = x55 * x69
    x78 = x76 + x77
    x79 = x43 * x56**2
    x80 = x42 + x79
    x81 = x3 * (x58 + x68) + x56 * x60
    x82 = x3 * (x44 + x68) + x56 * x63
    x83 = 3 * x42
    x84 = x3 * (x45 + x59 + x62 + x83) + x56 * x66
    x85 = x11 * x30
    x86 = numpy.pi * x0 * x28
    x87 = x11 * x86
    x88 = x49 * x69**2
    x89 = x52 + x88
    x90 = x3 * (x67 + x73) + x69 * x72
    x91 = x3 * (x53 + x73) + x69 * x75
    x92 = 3 * x52
    x93 = x3 * (x54 + x71 + x74 + x92) + x69 * x78
    x94 = 2 * x42 * x56 + x56 * x80
    x95 = x79 + x83
    x96 = x30 * x5
    x97 = x96 * (x3 * (2 * x59 + x95) + x56 * x81)
    x98 = x94 * x96
    x99 = x96 * (x3 * (2 * x62 + x95) + x56 * x82)
    x100 = x69 * x96
    x101 = x5 * x86
    x102 = x101 * x56
    x103 = 2 * x52 * x69 + x69 * x89
    x104 = x101 * x103
    x105 = x88 + x92
    x106 = x101 * (x3 * (x105 + 2 * x71) + x69 * x90)
    x107 = x101 * (x3 * (x105 + 2 * x74) + x69 * x91)

    # 90 item(s)
    S = numpy.array(
        [
            x31 * (x2 * x23 + x3 * (2 * x19 + 2 * x21 + x25 + x27)),
            x33 * x36,
            x36 * x38,
            x39 * x40,
            x46 * x48 * x49,
            x38 * x39 * x50,
            x40 * x51,
            x33 * x50 * x51,
            x43 * x48 * x55,
            x56 * x57,
            x25 * x49 * x60,
            x25 * x38 * x61,
            x27 * x49 * x63,
            x47 * x49 * x66,
            x47 * x63 * x67,
            x27 * x51 * x61,
            x47 * x53 * x60,
            x47 * x55 * x68,
            x57 * x69,
            x25 * x33 * x70,
            x25 * x43 * x72,
            x27 * x39 * x70,
            x46 * x47 * x73,
            x44 * x47 * x72,
            x27 * x43 * x75,
            x47 * x58 * x75,
            x43 * x47 * x78,
            x22 * x49 * x80,
            x24 * x49 * x81,
            x24 * x67 * x80,
            x26 * x49 * x82,
            x84 * x85,
            x38 * x82 * x85,
            x26 * x53 * x80,
            x51 * x81 * x85,
            x12 * x55 * x80,
            x22 * x61 * x69,
            x24 * x60 * x73,
            x24 * x68 * x72,
            x26 * x63 * x73,
            x66 * x69 * x85,
            x12 * x63 * x72,
            x26 * x68 * x75,
            x12 * x60 * x75,
            x56 * x78 * x87,
            x22 * x43 * x89,
            x24 * x58 * x89,
            x24 * x43 * x90,
            x26 * x44 * x89,
            x12 * x46 * x89,
            x39 * x87 * x90,
            x26 * x43 * x91,
            x33 * x87 * x91,
            x87 * x93,
            x20 * x49 * x94,
            x10 * x97,
            x10 * x38 * x98,
            x14 * x99,
            x96 * (x3 * (2 * x64 + 2 * x65 + x81 + x82) + x56 * x84),
            x38 * x99,
            x14 * x51 * x98,
            x51 * x97,
            x55 * x7 * x94,
            x20 * x73 * x80,
            x10 * x100 * x81,
            x16 * x72 * x80,
            x100 * x14 * x82,
            x100 * x84,
            x7 * x72 * x82,
            x18 * x75 * x80,
            x7 * x75 * x81,
            x7 * x78 * x80,
            x20 * x68 * x89,
            x16 * x60 * x89,
            x10 * x102 * x90,
            x18 * x63 * x89,
            x66 * x7 * x89,
            x63 * x7 * x90,
            x102 * x14 * x91,
            x60 * x7 * x91,
            x102 * x93,
            x103 * x20 * x43,
            x10 * x104 * x33,
            x10 * x106,
            x104 * x14 * x39,
            x103 * x46 * x7,
            x106 * x39,
            x107 * x14,
            x107 * x33,
            x101 * (x3 * (2 * x76 + 2 * x77 + x90 + x91) + x69 * x93),
        ]
    )
    return S


def dipole3d_32(a, A, b, B, C):
    """Cartesian 3D (fd) dipole moment integral.
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
    x12 = x5**2 * x9
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
    x43 = -x42 - C[1]
    x44 = x2 * x9
    x45 = x3 * (x10 + x44)
    x46 = x13 + x31
    x47 = x2 * x46
    x48 = x41 * (x2 * x37 + x3 * (4 * x26 + 2 * x28 + 2 * x45 + 2 * x47))
    x49 = -x0 * (a * A[2] + b * B[2])
    x50 = -x49 - C[2]
    x51 = -x42 - B[1]
    x52 = x33 + x35
    x53 = x45 + x47
    x54 = x13 + x32
    x55 = x2 * x54 + x3 * (x17 + x44)
    x56 = x41 * (x2 * x52 + x3 * (2 * x18 + x25 + x53 + x55))
    x57 = x3 * x8
    x58 = x38 * x57
    x59 = x38 * x8
    x60 = x51 * x59
    x61 = x43 * x60
    x62 = x58 + x61
    x63 = x2**2 * x9
    x64 = x14 + x63
    x65 = x2 * x53 + x3 * (x36 + x64)
    x66 = x39 * x8
    x67 = x41 * x65
    x68 = -x49 - B[2]
    x69 = x39 * x57
    x70 = x66 * x68
    x71 = x50 * x70
    x72 = x69 + x71
    x73 = x51**2 * x59
    x74 = x58 + x73
    x75 = x2 * x55 + x3 * (2 * x32 + x64)
    x76 = x43 * x59
    x77 = x3 * (x60 + x76)
    x78 = x51 * x62
    x79 = x77 + x78
    x80 = x13 + x63
    x81 = 2 * x13 * x2 + x2 * x80
    x82 = x50 * x66
    x83 = x41 * x68
    x84 = x66 * x68**2
    x85 = x69 + x84
    x86 = x3 * (x70 + x82)
    x87 = x68 * x72
    x88 = x86 + x87
    x89 = -x42 - A[1]
    x90 = x30 * x41
    x91 = x76 * x89
    x92 = x58 + x91
    x93 = x37 * x41
    x94 = x60 * x89
    x95 = x58 + x94
    x96 = x62 * x89
    x97 = x77 + x96
    x98 = x59 * x89
    x99 = 2 * x58
    x100 = x74 * x89
    x101 = x100 + x51 * x99
    x102 = 3 * x58
    x103 = x102 + x73
    x104 = x3 * (x103 + 2 * x61)
    x105 = x79 * x89
    x106 = x104 + x105
    x107 = -x49 - A[2]
    x108 = x107 * x82
    x109 = x108 + x69
    x110 = x107 * x41
    x111 = x107 * x66
    x112 = x107 * x70
    x113 = x112 + x69
    x114 = x107 * x72
    x115 = x114 + x86
    x116 = 2 * x69
    x117 = x107 * x85
    x118 = x116 * x68 + x117
    x119 = 3 * x69
    x120 = x119 + x84
    x121 = x3 * (x120 + 2 * x71)
    x122 = x107 * x88
    x123 = x121 + x122
    x124 = x59 * x89**2
    x125 = x124 + x58
    x126 = x3 * (x76 + x98) + x89 * x92
    x127 = x3 * (x60 + x98)
    x128 = x89 * x95
    x129 = x127 + x128
    x130 = x3 * (x102 + x61 + x91 + x94)
    x131 = x89 * x97
    x132 = x130 + x131
    x133 = 2 * x94
    x134 = x101 * x89 + x3 * (x103 + x133)
    x135 = 2 * x96
    x136 = x106 * x89 + x3 * (x101 + x135 + 3 * x77 + x78)
    x137 = x40 * x7
    x138 = x136 * x137
    x139 = x137 * x2
    x140 = numpy.pi * x0 * x38 * x7
    x141 = x140 * x2
    x142 = x107**2 * x66
    x143 = x142 + x69
    x144 = x107 * x109 + x3 * (x111 + x82)
    x145 = x3 * (x111 + x70)
    x146 = x107 * x113
    x147 = x145 + x146
    x148 = x3 * (x108 + x112 + x119 + x71)
    x149 = x107 * x115
    x150 = x148 + x149
    x151 = 2 * x112
    x152 = x107 * x118 + x3 * (x120 + x151)
    x153 = 2 * x114
    x154 = x107 * x123 + x3 * (x118 + x153 + 3 * x86 + x87)
    x155 = x140 * x154
    x156 = x125 * x89 + x89 * x99
    x157 = x102 + x124
    x158 = x126 * x89 + x3 * (x157 + 2 * x91)
    x159 = x129 * x89 + x3 * (x133 + x157)
    x160 = x137 * (x132 * x89 + x3 * (x126 + x129 + x135 + 2 * x77))
    x161 = x137 * x5
    x162 = x137 * (x134 * x89 + x3 * (2 * x100 + 2 * x127 + 2 * x128 + 4 * x51 * x58))
    x163 = x137 * x4
    x164 = x140 * x89
    x165 = x107 * x116 + x107 * x143
    x166 = x119 + x142
    x167 = x107 * x144 + x3 * (2 * x108 + x166)
    x168 = x140 * x5
    x169 = x107 * x147 + x3 * (x151 + x166)
    x170 = x140 * (x107 * x150 + x3 * (x144 + x147 + x153 + 2 * x86))
    x171 = x140 * (x107 * x152 + x3 * (2 * x117 + 2 * x145 + 2 * x146 + 4 * x68 * x69))

    # 180 item(s)
    S = numpy.array(
        [
            x41 * (x2 * x30 + x3 * (2 * x16 + 2 * x22 + 2 * x33 + 2 * x35 + x37)),
            x43 * x48,
            x48 * x50,
            x51 * x56,
            x62 * x65 * x66,
            x50 * x51 * x67,
            x56 * x68,
            x43 * x67 * x68,
            x59 * x65 * x72,
            x66 * x74 * x75,
            x66 * x79 * x81,
            x74 * x81 * x82,
            x51 * x75 * x83,
            x62 * x70 * x81,
            x60 * x72 * x81,
            x59 * x75 * x85,
            x76 * x81 * x85,
            x59 * x81 * x88,
            x89 * x90,
            x37 * x66 * x92,
            x50 * x89 * x93,
            x52 * x66 * x95,
            x53 * x66 * x97,
            x53 * x82 * x95,
            x52 * x83 * x89,
            x53 * x70 * x92,
            x53 * x72 * x98,
            x101 * x55 * x66,
            x106 * x66 * x80,
            x101 * x80 * x82,
            x55 * x70 * x95,
            x70 * x80 * x97,
            x72 * x80 * x95,
            x55 * x85 * x98,
            x80 * x85 * x92,
            x80 * x88 * x98,
            x107 * x90,
            x107 * x43 * x93,
            x109 * x37 * x59,
            x110 * x51 * x52,
            x111 * x53 * x62,
            x109 * x53 * x60,
            x113 * x52 * x59,
            x113 * x53 * x76,
            x115 * x53 * x59,
            x111 * x55 * x74,
            x111 * x79 * x80,
            x109 * x74 * x80,
            x113 * x55 * x60,
            x113 * x62 * x80,
            x115 * x60 * x80,
            x118 * x55 * x59,
            x118 * x76 * x80,
            x123 * x59 * x80,
            x125 * x23 * x66,
            x126 * x29 * x66,
            x125 * x29 * x82,
            x129 * x34 * x66,
            x132 * x46 * x66,
            x129 * x46 * x82,
            x125 * x34 * x70,
            x126 * x46 * x70,
            x125 * x46 * x72,
            x134 * x54 * x66,
            x138 * x2,
            x134 * x139 * x50,
            x129 * x54 * x70,
            x132 * x139 * x68,
            x129 * x44 * x72,
            x125 * x54 * x85,
            x126 * x44 * x85,
            x125 * x44 * x88,
            x110 * x23 * x89,
            x111 * x29 * x92,
            x109 * x29 * x98,
            x111 * x34 * x95,
            x111 * x46 * x97,
            x109 * x46 * x95,
            x113 * x34 * x98,
            x113 * x46 * x92,
            x115 * x46 * x98,
            x101 * x111 * x54,
            x106 * x107 * x139,
            x101 * x109 * x44,
            x113 * x54 * x95,
            x113 * x44 * x97,
            x115 * x44 * x95,
            x118 * x54 * x98,
            x118 * x44 * x92,
            x123 * x141 * x89,
            x143 * x23 * x59,
            x143 * x29 * x76,
            x144 * x29 * x59,
            x143 * x34 * x60,
            x143 * x46 * x62,
            x144 * x46 * x60,
            x147 * x34 * x59,
            x147 * x46 * x76,
            x150 * x46 * x59,
            x143 * x54 * x74,
            x143 * x44 * x79,
            x144 * x44 * x74,
            x147 * x54 * x60,
            x147 * x44 * x62,
            x141 * x150 * x51,
            x152 * x54 * x59,
            x141 * x152 * x43,
            x155 * x2,
            x156 * x21 * x66,
            x158 * x27 * x66,
            x156 * x27 * x82,
            x159 * x19 * x66,
            x160 * x5,
            x159 * x161 * x50,
            x156 * x19 * x70,
            x158 * x161 * x68,
            x10 * x156 * x72,
            x162 * x4,
            x137 * (x136 * x89 + x3 * (2 * x104 + 2 * x105 + 2 * x130 + 2 * x131 + x134)),
            x162 * x50,
            x159 * x163 * x68,
            x160 * x68,
            x159 * x72 * x9,
            x156 * x17 * x85,
            x158 * x85 * x9,
            x156 * x88 * x9,
            x111 * x125 * x21,
            x111 * x126 * x27,
            x109 * x125 * x27,
            x111 * x129 * x19,
            x107 * x132 * x161,
            x10 * x109 * x129,
            x113 * x125 * x19,
            x10 * x113 * x126,
            x10 * x115 * x125,
            x107 * x134 * x163,
            x107 * x138,
            x109 * x134 * x9,
            x113 * x129 * x17,
            x113 * x132 * x9,
            x115 * x129 * x9,
            x118 * x125 * x17,
            x118 * x126 * x9,
            x123 * x125 * x9,
            x143 * x21 * x98,
            x143 * x27 * x92,
            x144 * x27 * x98,
            x143 * x19 * x95,
            x10 * x143 * x97,
            x10 * x144 * x95,
            x147 * x19 * x98,
            x10 * x147 * x92,
            x150 * x164 * x5,
            x101 * x143 * x17,
            x106 * x143 * x9,
            x101 * x144 * x9,
            x147 * x17 * x95,
            x147 * x9 * x97,
            x150 * x9 * x95,
            x152 * x164 * x4,
            x152 * x9 * x92,
            x155 * x89,
            x165 * x21 * x59,
            x165 * x27 * x76,
            x167 * x27 * x59,
            x165 * x19 * x60,
            x10 * x165 * x62,
            x167 * x168 * x51,
            x169 * x19 * x59,
            x168 * x169 * x43,
            x170 * x5,
            x165 * x17 * x74,
            x165 * x79 * x9,
            x167 * x74 * x9,
            x140 * x169 * x4 * x51,
            x169 * x62 * x9,
            x170 * x51,
            x171 * x4,
            x171 * x43,
            x140
            * (x107 * x154 + x3 * (2 * x121 + 2 * x122 + 2 * x148 + 2 * x149 + x152)),
        ]
    )
    return S


def dipole3d_33(a, A, b, B, C):
    """Cartesian 3D (ff) dipole moment integral.
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
    x20 = x4**2 * x8
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
    x51 = -x50 - C[1]
    x52 = x2 * x9
    x53 = 2 * x52
    x54 = x3 * (x26 + x53)
    x55 = x2 * x41
    x56 = x49 * (x2 * x45 + x3 * (2 * x34 + 2 * x35 + 3 * x54 + 3 * x55))
    x57 = -x0 * (a * A[2] + b * B[2])
    x58 = -x57 - C[2]
    x59 = -x50 - B[1]
    x60 = x42 + x44
    x61 = x11 * x2
    x62 = x3 * (x15 + x25 + x52 + x61)
    x63 = x12 + x38
    x64 = x2 * x63
    x65 = x54 + x55
    x66 = x49 * (x2 * x60 + x3 * (2 * x27 + 2 * x33 + 2 * x62 + 2 * x64 + x65))
    x67 = x3 * x7
    x68 = x46 * x67
    x69 = x46 * x7
    x70 = x59 * x69
    x71 = x51 * x70
    x72 = x68 + x71
    x73 = x2 * x8
    x74 = x3 * (x73 + x9)
    x75 = x14 + x52
    x76 = x2 * x75
    x77 = x2 * x65 + x3 * (4 * x18 + 2 * x40 + 2 * x74 + 2 * x76)
    x78 = x47 * x7
    x79 = x49 * x77
    x80 = -x57 - B[2]
    x81 = x47 * x67
    x82 = x78 * x80
    x83 = x58 * x82
    x84 = x81 + x83
    x85 = x59**2 * x69
    x86 = x68 + x85
    x87 = x62 + x64
    x88 = x74 + x76
    x89 = x14 + x61
    x90 = x2 * x89 + x3 * (x11 + x73)
    x91 = x2 * x87 + x3 * (2 * x12 + x39 + x88 + x90)
    x92 = x51 * x69
    x93 = x3 * (x70 + x92)
    x94 = x59 * x72
    x95 = x93 + x94
    x96 = x2**2 * x8
    x97 = x25 + x96
    x98 = x2 * x88 + x3 * (x53 + x97)
    x99 = x58 * x78
    x100 = x49 * x80
    x101 = x78 * x80**2
    x102 = x101 + x81
    x103 = x3 * (x82 + x99)
    x104 = x80 * x84
    x105 = x103 + x104
    x106 = x59 * x68
    x107 = 2 * x106
    x108 = x59 * x86
    x109 = x107 + x108
    x110 = x2 * x90 + x3 * (2 * x61 + x97)
    x111 = 3 * x68
    x112 = x111 + x85
    x113 = x3 * (x112 + 2 * x71)
    x114 = x59 * x95
    x115 = x113 + x114
    x116 = x14 + x96
    x117 = x116 * x2 + 2 * x14 * x2
    x118 = x80 * x81
    x119 = 2 * x118
    x120 = x102 * x80
    x121 = x119 + x120
    x122 = 3 * x81
    x123 = x101 + x122
    x124 = x3 * (x123 + 2 * x83)
    x125 = x105 * x80
    x126 = x124 + x125
    x127 = -x50 - A[1]
    x128 = x37 * x49
    x129 = x127 * x92
    x130 = x129 + x68
    x131 = x45 * x49
    x132 = x127 * x70
    x133 = x132 + x68
    x134 = x127 * x72
    x135 = x134 + x93
    x136 = x127 * x69
    x137 = x127 * x86
    x138 = x107 + x137
    x139 = x127 * x95
    x140 = x113 + x139
    x141 = x3 * (x111 + 3 * x85)
    x142 = x109 * x127
    x143 = x141 + x142
    x144 = 3 * x93
    x145 = x3 * (x109 + x144 + 3 * x94)
    x146 = x115 * x127
    x147 = x145 + x146
    x148 = -x57 - A[2]
    x149 = x148 * x99
    x150 = x149 + x81
    x151 = x148 * x49
    x152 = x148 * x78
    x153 = x148 * x82
    x154 = x153 + x81
    x155 = x148 * x84
    x156 = x103 + x155
    x157 = x102 * x148
    x158 = x119 + x157
    x159 = x105 * x148
    x160 = x124 + x159
    x161 = x3 * (3 * x101 + x122)
    x162 = x121 * x148
    x163 = x161 + x162
    x164 = 3 * x103
    x165 = x3 * (3 * x104 + x121 + x164)
    x166 = x126 * x148
    x167 = x165 + x166
    x168 = x127**2 * x69
    x169 = x168 + x68
    x170 = x127 * x130 + x3 * (x136 + x92)
    x171 = x3 * (x136 + x70)
    x172 = x127 * x133
    x173 = x171 + x172
    x174 = x3 * (x111 + x129 + x132 + x71)
    x175 = x127 * x135
    x176 = x174 + x175
    x177 = 2 * x132
    x178 = x3 * (x112 + x177)
    x179 = x127 * x138
    x180 = x178 + x179
    x181 = x127 * x140
    x182 = 2 * x134
    x183 = x3 * (x138 + x144 + x182 + x94)
    x184 = x181 + x183
    x185 = x127 * x143 + x3 * (8 * x106 + x108 + 3 * x137)
    x186 = x127 * x147 + x3 * (4 * x113 + x114 + 3 * x139 + x143)
    x187 = x48 * x6
    x188 = x186 * x187
    x189 = x187 * x2
    x190 = numpy.pi * x0 * x46 * x6
    x191 = x190 * x2
    x192 = x148**2 * x78
    x193 = x192 + x81
    x194 = x148 * x150 + x3 * (x152 + x99)
    x195 = x3 * (x152 + x82)
    x196 = x148 * x154
    x197 = x195 + x196
    x198 = x3 * (x122 + x149 + x153 + x83)
    x199 = x148 * x156
    x200 = x198 + x199
    x201 = 2 * x153
    x202 = x3 * (x123 + x201)
    x203 = x148 * x158
    x204 = x202 + x203
    x205 = x148 * x160
    x206 = 2 * x155
    x207 = x3 * (x104 + x158 + x164 + x206)
    x208 = x205 + x207
    x209 = x148 * x163 + x3 * (8 * x118 + x120 + 3 * x157)
    x210 = x148 * x167 + x3 * (4 * x124 + x125 + 3 * x159 + x163)
    x211 = x190 * x210
    x212 = x127 * x169 + 2 * x127 * x68
    x213 = x111 + x168
    x214 = x127 * x170 + x3 * (2 * x129 + x213)
    x215 = x127 * x173 + x3 * (x177 + x213)
    x216 = x127 * x176 + x3 * (x170 + x173 + x182 + 2 * x93)
    x217 = x127 * x180 + x3 * (4 * x106 + 2 * x137 + 2 * x171 + 2 * x172)
    x218 = x187 * (x127 * x184 + x3 * (2 * x113 + 2 * x139 + 2 * x174 + 2 * x175 + x180))
    x219 = x187 * x4
    x220 = x187 * (x127 * x185 + x3 * (2 * x141 + 2 * x142 + 3 * x178 + 3 * x179))
    x221 = x10 * x187
    x222 = x127 * x190
    x223 = x148 * x193 + 2 * x148 * x81
    x224 = x122 + x192
    x225 = x148 * x194 + x3 * (2 * x149 + x224)
    x226 = x148 * x197 + x3 * (x201 + x224)
    x227 = x148 * x200 + x3 * (2 * x103 + x194 + x197 + x206)
    x228 = x190 * x4
    x229 = x148 * x204 + x3 * (4 * x118 + 2 * x157 + 2 * x195 + 2 * x196)
    x230 = x190 * (x148 * x208 + x3 * (2 * x124 + 2 * x159 + 2 * x198 + 2 * x199 + x204))
    x231 = x190 * (x148 * x209 + x3 * (2 * x161 + 2 * x162 + 3 * x202 + 3 * x203))

    # 300 item(s)
    S = numpy.array(
        [
            x49 * (x2 * x37 + x3 * (2 * x24 + 2 * x31 + 3 * x42 + 3 * x44 + x45)),
            x51 * x56,
            x56 * x58,
            x59 * x66,
            x72 * x77 * x78,
            x58 * x59 * x79,
            x66 * x80,
            x51 * x79 * x80,
            x69 * x77 * x84,
            x78 * x86 * x91,
            x78 * x95 * x98,
            x86 * x98 * x99,
            x100 * x59 * x91,
            x72 * x82 * x98,
            x70 * x84 * x98,
            x102 * x69 * x91,
            x102 * x92 * x98,
            x105 * x69 * x98,
            x109 * x110 * x78,
            x115 * x117 * x78,
            x109 * x117 * x99,
            x110 * x82 * x86,
            x117 * x82 * x95,
            x117 * x84 * x86,
            x102 * x110 * x70,
            x102 * x117 * x72,
            x105 * x117 * x70,
            x110 * x121 * x69,
            x117 * x121 * x92,
            x117 * x126 * x69,
            x127 * x128,
            x130 * x45 * x78,
            x127 * x131 * x58,
            x133 * x60 * x78,
            x135 * x65 * x78,
            x133 * x65 * x99,
            x100 * x127 * x60,
            x130 * x65 * x82,
            x136 * x65 * x84,
            x138 * x78 * x87,
            x140 * x78 * x88,
            x138 * x88 * x99,
            x133 * x82 * x87,
            x135 * x82 * x88,
            x133 * x84 * x88,
            x102 * x136 * x87,
            x102 * x130 * x88,
            x105 * x136 * x88,
            x143 * x78 * x90,
            x116 * x147 * x78,
            x116 * x143 * x99,
            x138 * x82 * x90,
            x116 * x140 * x82,
            x116 * x138 * x84,
            x102 * x133 * x90,
            x102 * x116 * x135,
            x105 * x116 * x133,
            x121 * x136 * x90,
            x116 * x121 * x130,
            x116 * x126 * x136,
            x128 * x148,
            x131 * x148 * x51,
            x150 * x45 * x69,
            x151 * x59 * x60,
            x152 * x65 * x72,
            x150 * x65 * x70,
            x154 * x60 * x69,
            x154 * x65 * x92,
            x156 * x65 * x69,
            x152 * x86 * x87,
            x152 * x88 * x95,
            x150 * x86 * x88,
            x154 * x70 * x87,
            x154 * x72 * x88,
            x156 * x70 * x88,
            x158 * x69 * x87,
            x158 * x88 * x92,
            x160 * x69 * x88,
            x109 * x152 * x90,
            x115 * x116 * x152,
            x109 * x116 * x150,
            x154 * x86 * x90,
            x116 * x154 * x95,
            x116 * x156 * x86,
            x158 * x70 * x90,
            x116 * x158 * x72,
            x116 * x160 * x70,
            x163 * x69 * x90,
            x116 * x163 * x92,
            x116 * x167 * x69,
            x169 * x32 * x78,
            x170 * x36 * x78,
            x169 * x36 * x99,
            x173 * x43 * x78,
            x176 * x41 * x78,
            x173 * x41 * x99,
            x169 * x43 * x82,
            x170 * x41 * x82,
            x169 * x41 * x84,
            x180 * x63 * x78,
            x184 * x75 * x78,
            x180 * x75 * x99,
            x173 * x63 * x82,
            x176 * x75 * x82,
            x173 * x75 * x84,
            x102 * x169 * x63,
            x102 * x170 * x75,
            x105 * x169 * x75,
            x185 * x78 * x89,
            x188 * x2,
            x185 * x189 * x58,
            x180 * x82 * x89,
            x184 * x189 * x80,
            x180 * x73 * x84,
            x102 * x173 * x89,
            x102 * x176 * x73,
            x105 * x173 * x73,
            x121 * x169 * x89,
            x121 * x170 * x73,
            x126 * x169 * x73,
            x127 * x151 * x32,
            x130 * x152 * x36,
            x136 * x150 * x36,
            x133 * x152 * x43,
            x135 * x152 * x41,
            x133 * x150 * x41,
            x136 * x154 * x43,
            x130 * x154 * x41,
            x136 * x156 * x41,
            x138 * x152 * x63,
            x140 * x152 * x75,
            x138 * x150 * x75,
            x133 * x154 * x63,
            x135 * x154 * x75,
            x133 * x156 * x75,
            x136 * x158 * x63,
            x130 * x158 * x75,
            x136 * x160 * x75,
            x143 * x152 * x89,
            x147 * x148 * x189,
            x143 * x150 * x73,
            x138 * x154 * x89,
            x140 * x154 * x73,
            x138 * x156 * x73,
            x133 * x158 * x89,
            x135 * x158 * x73,
            x133 * x160 * x73,
            x136 * x163 * x89,
            x130 * x163 * x73,
            x127 * x167 * x191,
            x193 * x32 * x69,
            x193 * x36 * x92,
            x194 * x36 * x69,
            x193 * x43 * x70,
            x193 * x41 * x72,
            x194 * x41 * x70,
            x197 * x43 * x69,
            x197 * x41 * x92,
            x200 * x41 * x69,
            x193 * x63 * x86,
            x193 * x75 * x95,
            x194 * x75 * x86,
            x197 * x63 * x70,
            x197 * x72 * x75,
            x200 * x70 * x75,
            x204 * x63 * x69,
            x204 * x75 * x92,
            x208 * x69 * x75,
            x109 * x193 * x89,
            x115 * x193 * x73,
            x109 * x194 * x73,
            x197 * x86 * x89,
            x197 * x73 * x95,
            x200 * x73 * x86,
            x204 * x70 * x89,
            x204 * x72 * x73,
            x191 * x208 * x59,
            x209 * x69 * x89,
            x191 * x209 * x51,
            x2 * x211,
            x212 * x30 * x78,
            x214 * x23 * x78,
            x212 * x23 * x99,
            x215 * x28 * x78,
            x21 * x216 * x78,
            x21 * x215 * x99,
            x212 * x28 * x82,
            x21 * x214 * x82,
            x21 * x212 * x84,
            x16 * x217 * x78,
            x218 * x4,
            x217 * x219 * x58,
            x16 * x215 * x82,
            x216 * x219 * x80,
            x215 * x84 * x9,
            x102 * x16 * x212,
            x102 * x214 * x9,
            x105 * x212 * x9,
            x10 * x220,
            x187
            * (x127 * x186 + x3 * (2 * x145 + 2 * x146 + 3 * x181 + 3 * x183 + x185)),
            x220 * x58,
            x217 * x221 * x80,
            x218 * x80,
            x217 * x8 * x84,
            x102 * x11 * x215,
            x102 * x216 * x8,
            x105 * x215 * x8,
            x11 * x121 * x212,
            x121 * x214 * x8,
            x126 * x212 * x8,
            x152 * x169 * x30,
            x152 * x170 * x23,
            x150 * x169 * x23,
            x152 * x173 * x28,
            x152 * x176 * x21,
            x150 * x173 * x21,
            x154 * x169 * x28,
            x154 * x170 * x21,
            x156 * x169 * x21,
            x152 * x16 * x180,
            x148 * x184 * x219,
            x150 * x180 * x9,
            x154 * x16 * x173,
            x154 * x176 * x9,
            x156 * x173 * x9,
            x158 * x16 * x169,
            x158 * x170 * x9,
            x160 * x169 * x9,
            x148 * x185 * x221,
            x148 * x188,
            x150 * x185 * x8,
            x11 * x154 * x180,
            x154 * x184 * x8,
            x156 * x180 * x8,
            x11 * x158 * x173,
            x158 * x176 * x8,
            x160 * x173 * x8,
            x11 * x163 * x169,
            x163 * x170 * x8,
            x167 * x169 * x8,
            x136 * x193 * x30,
            x130 * x193 * x23,
            x136 * x194 * x23,
            x133 * x193 * x28,
            x135 * x193 * x21,
            x133 * x194 * x21,
            x136 * x197 * x28,
            x130 * x197 * x21,
            x136 * x200 * x21,
            x138 * x16 * x193,
            x140 * x193 * x9,
            x138 * x194 * x9,
            x133 * x16 * x197,
            x135 * x197 * x9,
            x133 * x200 * x9,
            x136 * x16 * x204,
            x130 * x204 * x9,
            x208 * x222 * x4,
            x11 * x143 * x193,
            x147 * x193 * x8,
            x143 * x194 * x8,
            x11 * x138 * x197,
            x140 * x197 * x8,
            x138 * x200 * x8,
            x11 * x133 * x204,
            x135 * x204 * x8,
            x133 * x208 * x8,
            x10 * x209 * x222,
            x130 * x209 * x8,
            x127 * x211,
            x223 * x30 * x69,
            x223 * x23 * x92,
            x225 * x23 * x69,
            x223 * x28 * x70,
            x21 * x223 * x72,
            x21 * x225 * x70,
            x226 * x28 * x69,
            x21 * x226 * x92,
            x21 * x227 * x69,
            x16 * x223 * x86,
            x223 * x9 * x95,
            x225 * x86 * x9,
            x16 * x226 * x70,
            x226 * x72 * x9,
            x227 * x228 * x59,
            x16 * x229 * x69,
            x228 * x229 * x51,
            x230 * x4,
            x109 * x11 * x223,
            x115 * x223 * x8,
            x109 * x225 * x8,
            x11 * x226 * x86,
            x226 * x8 * x95,
            x227 * x8 * x86,
            x10 * x190 * x229 * x59,
            x229 * x72 * x8,
            x230 * x59,
            x10 * x231,
            x231 * x51,
            x190
            * (x148 * x210 + x3 * (2 * x165 + 2 * x166 + 3 * x205 + 3 * x207 + x209)),
        ]
    )
    return S


def dipole3d_34(a, A, b, B, C):
    """Cartesian 3D (fg) dipole moment integral.
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
    x12 = x5**2 * x9
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
    x58 = -x57 - C[1]
    x59 = x2 * x27
    x60 = x3 * (x28 + x41 + 3 * x59)
    x61 = x2 * x48
    x62 = x56 * (x2 * x52 + x3 * (2 * x42 + 2 * x43 + 4 * x60 + 4 * x61))
    x63 = -x0 * (a * A[2] + b * B[2])
    x64 = -x63 - C[2]
    x65 = -x57 - B[1]
    x66 = x49 + x51
    x67 = x2 * x20
    x68 = 2 * x67
    x69 = x26 + x59
    x70 = x3 * (x21 + x33 + x68 + x69)
    x71 = x16 + x46
    x72 = x2 * x71
    x73 = x60 + x61
    x74 = x56 * (x2 * x66 + x3 * (2 * x34 + 2 * x40 + 3 * x70 + 3 * x72 + x73))
    x75 = x3 * x8
    x76 = x53 * x75
    x77 = x53 * x8
    x78 = x65 * x77
    x79 = x58 * x78
    x80 = x76 + x79
    x81 = x10 * x2
    x82 = 2 * x81
    x83 = x3 * (x15 + x82)
    x84 = x2 * x69
    x85 = x2 * x73 + x3 * (2 * x24 + 2 * x47 + 3 * x83 + 3 * x84)
    x86 = x54 * x8
    x87 = x56 * x85
    x88 = -x63 - B[2]
    x89 = x54 * x75
    x90 = x86 * x88
    x91 = x64 * x90
    x92 = x89 + x91
    x93 = x65**2 * x77
    x94 = x76 + x93
    x95 = x70 + x72
    x96 = x18 * x2
    x97 = x3 * (x11 + x14 + x81 + x96)
    x98 = x19 + x67
    x99 = x2 * x98
    x100 = x83 + x84
    x101 = x2 * x95 + x3 * (x100 + 2 * x16 + 2 * x46 + 2 * x97 + 2 * x99)
    x102 = x58 * x77
    x103 = x3 * (x102 + x78)
    x104 = x65 * x80
    x105 = x103 + x104
    x106 = x2 * x9
    x107 = x3 * (x10 + x106)
    x108 = x13 + x81
    x109 = x108 * x2
    x110 = x100 * x2 + x3 * (2 * x107 + 2 * x109 + 4 * x25 + 2 * x59)
    x111 = x64 * x86
    x112 = x56 * x88
    x113 = x86 * x88**2
    x114 = x113 + x89
    x115 = x3 * (x111 + x90)
    x116 = x88 * x92
    x117 = x115 + x116
    x118 = x65 * x76
    x119 = 2 * x118
    x120 = x65 * x94
    x121 = x119 + x120
    x122 = x97 + x99
    x123 = x107 + x109
    x124 = x13 + x96
    x125 = x124 * x2 + x3 * (x106 + x18)
    x126 = x122 * x2 + x3 * (x123 + x125 + 2 * x19 + x68)
    x127 = 3 * x76
    x128 = x127 + x93
    x129 = x3 * (x128 + 2 * x79)
    x130 = x105 * x65
    x131 = x129 + x130
    x132 = x2**2 * x9
    x133 = x132 + x14
    x134 = x123 * x2 + x3 * (x133 + x82)
    x135 = x88 * x89
    x136 = 2 * x135
    x137 = x114 * x88
    x138 = x136 + x137
    x139 = 3 * x89
    x140 = x113 + x139
    x141 = x3 * (x140 + 2 * x91)
    x142 = x117 * x88
    x143 = x141 + x142
    x144 = x3 * (x127 + 3 * x93)
    x145 = x121 * x65
    x146 = x144 + x145
    x147 = x125 * x2 + x3 * (x133 + 2 * x96)
    x148 = 3 * x103
    x149 = x3 * (3 * x104 + x121 + x148)
    x150 = x131 * x65
    x151 = x149 + x150
    x152 = x13 + x132
    x153 = 2 * x13 * x2 + x152 * x2
    x154 = x3 * (3 * x113 + x139)
    x155 = x138 * x88
    x156 = x154 + x155
    x157 = 3 * x115
    x158 = x3 * (3 * x116 + x138 + x157)
    x159 = x143 * x88
    x160 = x158 + x159
    x161 = -x57 - A[1]
    x162 = x45 * x56
    x163 = x102 * x161
    x164 = x163 + x76
    x165 = x52 * x56
    x166 = x161 * x78
    x167 = x166 + x76
    x168 = x161 * x80
    x169 = x103 + x168
    x170 = x161 * x77
    x171 = x161 * x94
    x172 = x119 + x171
    x173 = x105 * x161
    x174 = x129 + x173
    x175 = x121 * x161
    x176 = x144 + x175
    x177 = x131 * x161
    x178 = x149 + x177
    x179 = 8 * x118
    x180 = x3 * (4 * x120 + x179)
    x181 = x146 * x161
    x182 = x180 + x181
    x183 = 4 * x129
    x184 = x3 * (4 * x130 + x146 + x183)
    x185 = x151 * x161
    x186 = x184 + x185
    x187 = -x63 - A[2]
    x188 = x111 * x187
    x189 = x188 + x89
    x190 = x187 * x56
    x191 = x187 * x86
    x192 = x187 * x90
    x193 = x192 + x89
    x194 = x187 * x92
    x195 = x115 + x194
    x196 = x114 * x187
    x197 = x136 + x196
    x198 = x117 * x187
    x199 = x141 + x198
    x200 = x138 * x187
    x201 = x154 + x200
    x202 = x143 * x187
    x203 = x158 + x202
    x204 = 8 * x135
    x205 = x3 * (4 * x137 + x204)
    x206 = x156 * x187
    x207 = x205 + x206
    x208 = 4 * x141
    x209 = x3 * (4 * x142 + x156 + x208)
    x210 = x160 * x187
    x211 = x209 + x210
    x212 = x161**2 * x77
    x213 = x212 + x76
    x214 = x161 * x164 + x3 * (x102 + x170)
    x215 = x3 * (x170 + x78)
    x216 = x161 * x167
    x217 = x215 + x216
    x218 = x3 * (x127 + x163 + x166 + x79)
    x219 = x161 * x169
    x220 = x218 + x219
    x221 = 2 * x166
    x222 = x3 * (x128 + x221)
    x223 = x161 * x172
    x224 = x222 + x223
    x225 = x161 * x174
    x226 = 2 * x168
    x227 = x3 * (x104 + x148 + x172 + x226)
    x228 = x225 + x227
    x229 = x3 * (x120 + 3 * x171 + x179)
    x230 = x161 * x176
    x231 = x229 + x230
    x232 = x161 * x178
    x233 = x3 * (x130 + 3 * x173 + x176 + x183)
    x234 = x232 + x233
    x235 = x161 * x182 + x3 * (5 * x144 + x145 + 4 * x175)
    x236 = x161 * x186 + x3 * (5 * x149 + x150 + 4 * x177 + x182)
    x237 = x55 * x7
    x238 = x236 * x237
    x239 = x2 * x237
    x240 = numpy.pi * x0 * x53 * x7
    x241 = x2 * x240
    x242 = x187**2 * x86
    x243 = x242 + x89
    x244 = x187 * x189 + x3 * (x111 + x191)
    x245 = x3 * (x191 + x90)
    x246 = x187 * x193
    x247 = x245 + x246
    x248 = x3 * (x139 + x188 + x192 + x91)
    x249 = x187 * x195
    x250 = x248 + x249
    x251 = 2 * x192
    x252 = x3 * (x140 + x251)
    x253 = x187 * x197
    x254 = x252 + x253
    x255 = x187 * x199
    x256 = 2 * x194
    x257 = x3 * (x116 + x157 + x197 + x256)
    x258 = x255 + x257
    x259 = x3 * (x137 + 3 * x196 + x204)
    x260 = x187 * x201
    x261 = x259 + x260
    x262 = x187 * x203
    x263 = x3 * (x142 + 3 * x198 + x201 + x208)
    x264 = x262 + x263
    x265 = x187 * x207 + x3 * (5 * x154 + x155 + 4 * x200)
    x266 = x187 * x211 + x3 * (5 * x158 + x159 + 4 * x202 + x207)
    x267 = x240 * x266
    x268 = x161 * x213 + 2 * x161 * x76
    x269 = x127 + x212
    x270 = x161 * x214 + x3 * (2 * x163 + x269)
    x271 = x161 * x217 + x3 * (x221 + x269)
    x272 = x161 * x220 + x3 * (2 * x103 + x214 + x217 + x226)
    x273 = x161 * x224 + x3 * (4 * x118 + 2 * x171 + 2 * x215 + 2 * x216)
    x274 = x161 * x228 + x3 * (2 * x129 + 2 * x173 + 2 * x218 + 2 * x219 + x224)
    x275 = x161 * x231 + x3 * (2 * x144 + 2 * x175 + 3 * x222 + 3 * x223)
    x276 = x237 * (x161 * x234 + x3 * (2 * x149 + 2 * x177 + 3 * x225 + 3 * x227 + x231))
    x277 = x237 * x5
    x278 = x237 * (x161 * x235 + x3 * (2 * x180 + 2 * x181 + 4 * x229 + 4 * x230))
    x279 = x237 * x4
    x280 = x161 * x240
    x281 = x187 * x243 + 2 * x187 * x89
    x282 = x139 + x242
    x283 = x187 * x244 + x3 * (2 * x188 + x282)
    x284 = x187 * x247 + x3 * (x251 + x282)
    x285 = x187 * x250 + x3 * (2 * x115 + x244 + x247 + x256)
    x286 = x187 * x254 + x3 * (4 * x135 + 2 * x196 + 2 * x245 + 2 * x246)
    x287 = x187 * x258 + x3 * (2 * x141 + 2 * x198 + 2 * x248 + 2 * x249 + x254)
    x288 = x240 * x5
    x289 = x187 * x261 + x3 * (2 * x154 + 2 * x200 + 3 * x252 + 3 * x253)
    x290 = x240 * (x187 * x264 + x3 * (2 * x158 + 2 * x202 + 3 * x255 + 3 * x257 + x261))
    x291 = x240 * (x187 * x265 + x3 * (2 * x205 + 2 * x206 + 4 * x259 + 4 * x260))

    # 450 item(s)
    S = numpy.array(
        [
            x56 * (x2 * x45 + x3 * (2 * x32 + 2 * x38 + 4 * x49 + 4 * x51 + x52)),
            x58 * x62,
            x62 * x64,
            x65 * x74,
            x80 * x85 * x86,
            x64 * x65 * x87,
            x74 * x88,
            x58 * x87 * x88,
            x77 * x85 * x92,
            x101 * x86 * x94,
            x105 * x110 * x86,
            x110 * x111 * x94,
            x101 * x112 * x65,
            x110 * x80 * x90,
            x110 * x78 * x92,
            x101 * x114 * x77,
            x102 * x110 * x114,
            x110 * x117 * x77,
            x121 * x126 * x86,
            x131 * x134 * x86,
            x111 * x121 * x134,
            x126 * x90 * x94,
            x105 * x134 * x90,
            x134 * x92 * x94,
            x114 * x126 * x78,
            x114 * x134 * x80,
            x117 * x134 * x78,
            x126 * x138 * x77,
            x102 * x134 * x138,
            x134 * x143 * x77,
            x146 * x147 * x86,
            x151 * x153 * x86,
            x111 * x146 * x153,
            x121 * x147 * x90,
            x131 * x153 * x90,
            x121 * x153 * x92,
            x114 * x147 * x94,
            x105 * x114 * x153,
            x117 * x153 * x94,
            x138 * x147 * x78,
            x138 * x153 * x80,
            x143 * x153 * x78,
            x147 * x156 * x77,
            x102 * x153 * x156,
            x153 * x160 * x77,
            x161 * x162,
            x164 * x52 * x86,
            x161 * x165 * x64,
            x167 * x66 * x86,
            x169 * x73 * x86,
            x111 * x167 * x73,
            x112 * x161 * x66,
            x164 * x73 * x90,
            x170 * x73 * x92,
            x172 * x86 * x95,
            x100 * x174 * x86,
            x100 * x111 * x172,
            x167 * x90 * x95,
            x100 * x169 * x90,
            x100 * x167 * x92,
            x114 * x170 * x95,
            x100 * x114 * x164,
            x100 * x117 * x170,
            x122 * x176 * x86,
            x123 * x178 * x86,
            x111 * x123 * x176,
            x122 * x172 * x90,
            x123 * x174 * x90,
            x123 * x172 * x92,
            x114 * x122 * x167,
            x114 * x123 * x169,
            x117 * x123 * x167,
            x122 * x138 * x170,
            x123 * x138 * x164,
            x123 * x143 * x170,
            x125 * x182 * x86,
            x152 * x186 * x86,
            x111 * x152 * x182,
            x125 * x176 * x90,
            x152 * x178 * x90,
            x152 * x176 * x92,
            x114 * x125 * x172,
            x114 * x152 * x174,
            x117 * x152 * x172,
            x125 * x138 * x167,
            x138 * x152 * x169,
            x143 * x152 * x167,
            x125 * x156 * x170,
            x152 * x156 * x164,
            x152 * x160 * x170,
            x162 * x187,
            x165 * x187 * x58,
            x189 * x52 * x77,
            x190 * x65 * x66,
            x191 * x73 * x80,
            x189 * x73 * x78,
            x193 * x66 * x77,
            x102 * x193 * x73,
            x195 * x73 * x77,
            x191 * x94 * x95,
            x100 * x105 * x191,
            x100 * x189 * x94,
            x193 * x78 * x95,
            x100 * x193 * x80,
            x100 * x195 * x78,
            x197 * x77 * x95,
            x100 * x102 * x197,
            x100 * x199 * x77,
            x121 * x122 * x191,
            x123 * x131 * x191,
            x121 * x123 * x189,
            x122 * x193 * x94,
            x105 * x123 * x193,
            x123 * x195 * x94,
            x122 * x197 * x78,
            x123 * x197 * x80,
            x123 * x199 * x78,
            x122 * x201 * x77,
            x102 * x123 * x201,
            x123 * x203 * x77,
            x125 * x146 * x191,
            x151 * x152 * x191,
            x146 * x152 * x189,
            x121 * x125 * x193,
            x131 * x152 * x193,
            x121 * x152 * x195,
            x125 * x197 * x94,
            x105 * x152 * x197,
            x152 * x199 * x94,
            x125 * x201 * x78,
            x152 * x201 * x80,
            x152 * x203 * x78,
            x125 * x207 * x77,
            x102 * x152 * x207,
            x152 * x211 * x77,
            x213 * x39 * x86,
            x214 * x44 * x86,
            x111 * x213 * x44,
            x217 * x50 * x86,
            x220 * x48 * x86,
            x111 * x217 * x48,
            x213 * x50 * x90,
            x214 * x48 * x90,
            x213 * x48 * x92,
            x224 * x71 * x86,
            x228 * x69 * x86,
            x111 * x224 * x69,
            x217 * x71 * x90,
            x220 * x69 * x90,
            x217 * x69 * x92,
            x114 * x213 * x71,
            x114 * x214 * x69,
            x117 * x213 * x69,
            x231 * x86 * x98,
            x108 * x234 * x86,
            x108 * x111 * x231,
            x224 * x90 * x98,
            x108 * x228 * x90,
            x108 * x224 * x92,
            x114 * x217 * x98,
            x108 * x114 * x220,
            x108 * x117 * x217,
            x138 * x213 * x98,
            x108 * x138 * x214,
            x108 * x143 * x213,
            x124 * x235 * x86,
            x2 * x238,
            x235 * x239 * x64,
            x124 * x231 * x90,
            x234 * x239 * x88,
            x106 * x231 * x92,
            x114 * x124 * x224,
            x106 * x114 * x228,
            x106 * x117 * x224,
            x124 * x138 * x217,
            x106 * x138 * x220,
            x106 * x143 * x217,
            x124 * x156 * x213,
            x106 * x156 * x214,
            x106 * x160 * x213,
            x161 * x190 * x39,
            x164 * x191 * x44,
            x170 * x189 * x44,
            x167 * x191 * x50,
            x169 * x191 * x48,
            x167 * x189 * x48,
            x170 * x193 * x50,
            x164 * x193 * x48,
            x170 * x195 * x48,
            x172 * x191 * x71,
            x174 * x191 * x69,
            x172 * x189 * x69,
            x167 * x193 * x71,
            x169 * x193 * x69,
            x167 * x195 * x69,
            x170 * x197 * x71,
            x164 * x197 * x69,
            x170 * x199 * x69,
            x176 * x191 * x98,
            x108 * x178 * x191,
            x108 * x176 * x189,
            x172 * x193 * x98,
            x108 * x174 * x193,
            x108 * x172 * x195,
            x167 * x197 * x98,
            x108 * x169 * x197,
            x108 * x167 * x199,
            x170 * x201 * x98,
            x108 * x164 * x201,
            x108 * x170 * x203,
            x124 * x182 * x191,
            x186 * x187 * x239,
            x106 * x182 * x189,
            x124 * x176 * x193,
            x106 * x178 * x193,
            x106 * x176 * x195,
            x124 * x172 * x197,
            x106 * x174 * x197,
            x106 * x172 * x199,
            x124 * x167 * x201,
            x106 * x169 * x201,
            x106 * x167 * x203,
            x124 * x170 * x207,
            x106 * x164 * x207,
            x161 * x211 * x241,
            x243 * x39 * x77,
            x102 * x243 * x44,
            x244 * x44 * x77,
            x243 * x50 * x78,
            x243 * x48 * x80,
            x244 * x48 * x78,
            x247 * x50 * x77,
            x102 * x247 * x48,
            x250 * x48 * x77,
            x243 * x71 * x94,
            x105 * x243 * x69,
            x244 * x69 * x94,
            x247 * x71 * x78,
            x247 * x69 * x80,
            x250 * x69 * x78,
            x254 * x71 * x77,
            x102 * x254 * x69,
            x258 * x69 * x77,
            x121 * x243 * x98,
            x108 * x131 * x243,
            x108 * x121 * x244,
            x247 * x94 * x98,
            x105 * x108 * x247,
            x108 * x250 * x94,
            x254 * x78 * x98,
            x108 * x254 * x80,
            x108 * x258 * x78,
            x261 * x77 * x98,
            x102 * x108 * x261,
            x108 * x264 * x77,
            x124 * x146 * x243,
            x106 * x151 * x243,
            x106 * x146 * x244,
            x121 * x124 * x247,
            x106 * x131 * x247,
            x106 * x121 * x250,
            x124 * x254 * x94,
            x105 * x106 * x254,
            x106 * x258 * x94,
            x124 * x261 * x78,
            x106 * x261 * x80,
            x241 * x264 * x65,
            x124 * x265 * x77,
            x241 * x265 * x58,
            x2 * x267,
            x268 * x37 * x86,
            x270 * x31 * x86,
            x111 * x268 * x31,
            x271 * x35 * x86,
            x272 * x29 * x86,
            x111 * x271 * x29,
            x268 * x35 * x90,
            x270 * x29 * x90,
            x268 * x29 * x92,
            x22 * x273 * x86,
            x27 * x274 * x86,
            x111 * x27 * x273,
            x22 * x271 * x90,
            x27 * x272 * x90,
            x27 * x271 * x92,
            x114 * x22 * x268,
            x114 * x27 * x270,
            x117 * x268 * x27,
            x20 * x275 * x86,
            x276 * x5,
            x275 * x277 * x64,
            x20 * x273 * x90,
            x274 * x277 * x88,
            x10 * x273 * x92,
            x114 * x20 * x271,
            x10 * x114 * x272,
            x10 * x117 * x271,
            x138 * x20 * x268,
            x10 * x138 * x270,
            x10 * x143 * x268,
            x278 * x4,
            x237
            * (x161 * x236 + x3 * (2 * x184 + 2 * x185 + 4 * x232 + 4 * x233 + x235)),
            x278 * x64,
            x275 * x279 * x88,
            x276 * x88,
            x275 * x9 * x92,
            x114 * x18 * x273,
            x114 * x274 * x9,
            x117 * x273 * x9,
            x138 * x18 * x271,
            x138 * x272 * x9,
            x143 * x271 * x9,
            x156 * x18 * x268,
            x156 * x270 * x9,
            x160 * x268 * x9,
            x191 * x213 * x37,
            x191 * x214 * x31,
            x189 * x213 * x31,
            x191 * x217 * x35,
            x191 * x220 * x29,
            x189 * x217 * x29,
            x193 * x213 * x35,
            x193 * x214 * x29,
            x195 * x213 * x29,
            x191 * x22 * x224,
            x191 * x228 * x27,
            x189 * x224 * x27,
            x193 * x217 * x22,
            x193 * x220 * x27,
            x195 * x217 * x27,
            x197 * x213 * x22,
            x197 * x214 * x27,
            x199 * x213 * x27,
            x191 * x20 * x231,
            x187 * x234 * x277,
            x10 * x189 * x231,
            x193 * x20 * x224,
            x10 * x193 * x228,
            x10 * x195 * x224,
            x197 * x20 * x217,
            x10 * x197 * x220,
            x10 * x199 * x217,
            x20 * x201 * x213,
            x10 * x201 * x214,
            x10 * x203 * x213,
            x187 * x235 * x279,
            x187 * x238,
            x189 * x235 * x9,
            x18 * x193 * x231,
            x193 * x234 * x9,
            x195 * x231 * x9,
            x18 * x197 * x224,
            x197 * x228 * x9,
            x199 * x224 * x9,
            x18 * x201 * x217,
            x201 * x220 * x9,
            x203 * x217 * x9,
            x18 * x207 * x213,
            x207 * x214 * x9,
            x211 * x213 * x9,
            x170 * x243 * x37,
            x164 * x243 * x31,
            x170 * x244 * x31,
            x167 * x243 * x35,
            x169 * x243 * x29,
            x167 * x244 * x29,
            x170 * x247 * x35,
            x164 * x247 * x29,
            x170 * x250 * x29,
            x172 * x22 * x243,
            x174 * x243 * x27,
            x172 * x244 * x27,
            x167 * x22 * x247,
            x169 * x247 * x27,
            x167 * x250 * x27,
            x170 * x22 * x254,
            x164 * x254 * x27,
            x170 * x258 * x27,
            x176 * x20 * x243,
            x10 * x178 * x243,
            x10 * x176 * x244,
            x172 * x20 * x247,
            x10 * x174 * x247,
            x10 * x172 * x250,
            x167 * x20 * x254,
            x10 * x169 * x254,
            x10 * x167 * x258,
            x170 * x20 * x261,
            x10 * x164 * x261,
            x264 * x280 * x5,
            x18 * x182 * x243,
            x186 * x243 * x9,
            x182 * x244 * x9,
            x176 * x18 * x247,
            x178 * x247 * x9,
            x176 * x250 * x9,
            x172 * x18 * x254,
            x174 * x254 * x9,
            x172 * x258 * x9,
            x167 * x18 * x261,
            x169 * x261 * x9,
            x167 * x264 * x9,
            x265 * x280 * x4,
            x164 * x265 * x9,
            x161 * x267,
            x281 * x37 * x77,
            x102 * x281 * x31,
            x283 * x31 * x77,
            x281 * x35 * x78,
            x281 * x29 * x80,
            x283 * x29 * x78,
            x284 * x35 * x77,
            x102 * x284 * x29,
            x285 * x29 * x77,
            x22 * x281 * x94,
            x105 * x27 * x281,
            x27 * x283 * x94,
            x22 * x284 * x78,
            x27 * x284 * x80,
            x27 * x285 * x78,
            x22 * x286 * x77,
            x102 * x27 * x286,
            x27 * x287 * x77,
            x121 * x20 * x281,
            x10 * x131 * x281,
            x10 * x121 * x283,
            x20 * x284 * x94,
            x10 * x105 * x284,
            x10 * x285 * x94,
            x20 * x286 * x78,
            x10 * x286 * x80,
            x287 * x288 * x65,
            x20 * x289 * x77,
            x288 * x289 * x58,
            x290 * x5,
            x146 * x18 * x281,
            x151 * x281 * x9,
            x146 * x283 * x9,
            x121 * x18 * x284,
            x131 * x284 * x9,
            x121 * x285 * x9,
            x18 * x286 * x94,
            x105 * x286 * x9,
            x287 * x9 * x94,
            x240 * x289 * x4 * x65,
            x289 * x80 * x9,
            x290 * x65,
            x291 * x4,
            x291 * x58,
            x240
            * (x187 * x266 + x3 * (2 * x209 + 2 * x210 + 4 * x262 + 4 * x263 + x265)),
        ]
    )
    return S


def dipole3d_40(a, A, b, B, C):
    """Cartesian 3D (gs) dipole moment integral.
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
    x27 = -x26 - C[1]
    x28 = x25 * (x0 * (3 * x16 + x19) + x18 * x3)
    x29 = -x1 * (a * A[2] + b * B[2])
    x30 = -x29 - C[2]
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


def dipole3d_41(a, A, b, B, C):
    """Cartesian 3D (gp) dipole moment integral.
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
    x34 = x2**2 * x7
    x35 = x34 + x9
    x36 = x2 * x28 + x3 * (2 * x12 + x35)
    x37 = x2 * x32 + x3 * (2 * x14 + x35)
    x38 = numpy.exp(-x4 * (A[1] - B[1]) ** 2)
    x39 = numpy.exp(-x4 * (A[2] - B[2]) ** 2)
    x40 = numpy.pi * x0 * x39
    x41 = x38 * x40
    x42 = -x0 * (a * A[1] + b * B[1])
    x43 = -x42 - C[1]
    x44 = x34 + x8
    x45 = 2 * x11 * x3 + x2 * x44
    x46 = x41 * (x2 * x36 + x3 * (3 * x25 + 3 * x27 + x45))
    x47 = -x0 * (a * A[2] + b * B[2])
    x48 = -x47 - C[2]
    x49 = -x42 - B[1]
    x50 = x41 * (x2 * x37 + x3 * (3 * x29 + 3 * x31 + x45))
    x51 = x3 * x6
    x52 = x38 * x51
    x53 = x38 * x6
    x54 = x49 * x53
    x55 = x43 * x54
    x56 = x52 + x55
    x57 = x2 * x45 + x3 * (3 * x34 + x9)
    x58 = x39 * x6
    x59 = x41 * x57
    x60 = -x47 - B[2]
    x61 = x39 * x51
    x62 = x58 * x60
    x63 = x48 * x62
    x64 = x61 + x63
    x65 = -x42 - A[1]
    x66 = x33 * x41
    x67 = x43 * x53
    x68 = x65 * x67
    x69 = x52 + x68
    x70 = x41 * x65
    x71 = x54 * x65
    x72 = x52 + x71
    x73 = x3 * (x54 + x67)
    x74 = x56 * x65
    x75 = x73 + x74
    x76 = x48 * x58
    x77 = x53 * x65
    x78 = -x47 - A[2]
    x79 = x41 * x78
    x80 = x76 * x78
    x81 = x61 + x80
    x82 = x58 * x78
    x83 = x62 * x78
    x84 = x61 + x83
    x85 = x3 * (x62 + x76)
    x86 = x64 * x78
    x87 = x85 + x86
    x88 = x53 * x65**2
    x89 = x52 + x88
    x90 = x3 * (x67 + x77)
    x91 = x65 * x69
    x92 = x90 + x91
    x93 = x3 * (x54 + x77)
    x94 = x65 * x72
    x95 = x93 + x94
    x96 = 3 * x52
    x97 = x3 * (x55 + x68 + x71 + x96)
    x98 = x65 * x75
    x99 = x97 + x98
    x100 = x58 * x78**2
    x101 = x100 + x61
    x102 = x3 * (x76 + x82)
    x103 = x78 * x81
    x104 = x102 + x103
    x105 = x3 * (x62 + x82)
    x106 = x78 * x84
    x107 = x105 + x106
    x108 = 3 * x61
    x109 = x3 * (x108 + x63 + x80 + x83)
    x110 = x78 * x87
    x111 = x109 + x110
    x112 = 2 * x52 * x65 + x65 * x89
    x113 = x88 + x96
    x114 = x3 * (x113 + 2 * x68) + x65 * x92
    x115 = x3 * (x113 + 2 * x71) + x65 * x95
    x116 = x3 * (2 * x73 + 2 * x74 + x92 + x95) + x65 * x99
    x117 = x40 * x5
    x118 = x116 * x117
    x119 = x117 * x2
    x120 = numpy.pi * x0 * x38 * x5
    x121 = x120 * x2
    x122 = x101 * x78 + 2 * x61 * x78
    x123 = x100 + x108
    x124 = x104 * x78 + x3 * (x123 + 2 * x80)
    x125 = x107 * x78 + x3 * (x123 + 2 * x83)
    x126 = x111 * x78 + x3 * (x104 + x107 + 2 * x85 + 2 * x86)
    x127 = x120 * x126
    x128 = x112 * x65 + x3 * (3 * x88 + x96)
    x129 = x117 * (x114 * x65 + x3 * (x112 + 3 * x90 + 3 * x91))
    x130 = x117 * x128
    x131 = x117 * (x115 * x65 + x3 * (x112 + 3 * x93 + 3 * x94))
    x132 = x117 * x78
    x133 = x120 * x65
    x134 = x122 * x78 + x3 * (3 * x100 + x108)
    x135 = x120 * x134
    x136 = x120 * (x124 * x78 + x3 * (3 * x102 + 3 * x103 + x122))
    x137 = x120 * (x125 * x78 + x3 * (3 * x105 + 3 * x106 + x122))

    # 135 item(s)
    S = numpy.array(
        [
            x41 * (x2 * x33 + x3 * (3 * x17 + 3 * x23 + x36 + x37)),
            x43 * x46,
            x46 * x48,
            x49 * x50,
            x56 * x57 * x58,
            x48 * x49 * x59,
            x50 * x60,
            x43 * x59 * x60,
            x53 * x57 * x64,
            x65 * x66,
            x36 * x58 * x69,
            x36 * x48 * x70,
            x37 * x58 * x72,
            x45 * x58 * x75,
            x45 * x72 * x76,
            x37 * x60 * x70,
            x45 * x62 * x69,
            x45 * x64 * x77,
            x66 * x78,
            x36 * x43 * x79,
            x36 * x53 * x81,
            x37 * x49 * x79,
            x45 * x56 * x82,
            x45 * x54 * x81,
            x37 * x53 * x84,
            x45 * x67 * x84,
            x45 * x53 * x87,
            x24 * x58 * x89,
            x28 * x58 * x92,
            x28 * x76 * x89,
            x32 * x58 * x95,
            x44 * x58 * x99,
            x44 * x76 * x95,
            x32 * x62 * x89,
            x44 * x62 * x92,
            x44 * x64 * x89,
            x24 * x70 * x78,
            x28 * x69 * x82,
            x28 * x77 * x81,
            x32 * x72 * x82,
            x44 * x75 * x82,
            x44 * x72 * x81,
            x32 * x77 * x84,
            x44 * x69 * x84,
            x44 * x77 * x87,
            x101 * x24 * x53,
            x101 * x28 * x67,
            x104 * x28 * x53,
            x101 * x32 * x54,
            x101 * x44 * x56,
            x104 * x44 * x54,
            x107 * x32 * x53,
            x107 * x44 * x67,
            x111 * x44 * x53,
            x112 * x22 * x58,
            x114 * x26 * x58,
            x112 * x26 * x76,
            x115 * x30 * x58,
            x118 * x2,
            x115 * x119 * x48,
            x112 * x30 * x62,
            x114 * x119 * x60,
            x11 * x112 * x64,
            x22 * x82 * x89,
            x26 * x82 * x92,
            x26 * x81 * x89,
            x30 * x82 * x95,
            x119 * x78 * x99,
            x11 * x81 * x95,
            x30 * x84 * x89,
            x11 * x84 * x92,
            x11 * x87 * x89,
            x101 * x22 * x77,
            x101 * x26 * x69,
            x104 * x26 * x77,
            x101 * x30 * x72,
            x101 * x11 * x75,
            x104 * x11 * x72,
            x107 * x30 * x77,
            x107 * x11 * x69,
            x111 * x121 * x65,
            x122 * x22 * x53,
            x122 * x26 * x67,
            x124 * x26 * x53,
            x122 * x30 * x54,
            x11 * x122 * x56,
            x121 * x124 * x49,
            x125 * x30 * x53,
            x121 * x125 * x43,
            x127 * x2,
            x128 * x20 * x58,
            x10 * x129,
            x10 * x130 * x48,
            x13 * x131,
            x117 * (x116 * x65 + x3 * (x114 + x115 + 3 * x97 + 3 * x98)),
            x131 * x48,
            x13 * x130 * x60,
            x129 * x60,
            x128 * x64 * x7,
            x112 * x20 * x82,
            x10 * x114 * x132,
            x112 * x15 * x81,
            x115 * x13 * x132,
            x118 * x78,
            x115 * x7 * x81,
            x112 * x18 * x84,
            x114 * x7 * x84,
            x112 * x7 * x87,
            x101 * x20 * x89,
            x101 * x15 * x92,
            x104 * x15 * x89,
            x101 * x18 * x95,
            x101 * x7 * x99,
            x104 * x7 * x95,
            x107 * x18 * x89,
            x107 * x7 * x92,
            x111 * x7 * x89,
            x122 * x20 * x77,
            x122 * x15 * x69,
            x10 * x124 * x133,
            x122 * x18 * x72,
            x122 * x7 * x75,
            x124 * x7 * x72,
            x125 * x13 * x133,
            x125 * x69 * x7,
            x127 * x65,
            x134 * x20 * x53,
            x10 * x135 * x43,
            x10 * x136,
            x13 * x135 * x49,
            x134 * x56 * x7,
            x136 * x49,
            x13 * x137,
            x137 * x43,
            x120 * (x126 * x78 + x3 * (3 * x109 + 3 * x110 + x124 + x125)),
        ]
    )
    return S


def dipole3d_42(a, A, b, B, C):
    """Cartesian 3D (gd) dipole moment integral.
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
    x12 = x5**2 * x9
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
    x62 = -x61 - C[1]
    x63 = x2**2 * x9
    x64 = x14 + x63
    x65 = x3 * (x40 + x64)
    x66 = x2 * x48
    x67 = x60 * (x2 * x56 + x3 * (3 * x41 + 3 * x42 + 2 * x65 + 2 * x66))
    x68 = -x0 * (a * A[2] + b * B[2])
    x69 = -x68 - C[2]
    x70 = -x61 - B[1]
    x71 = x53 + x55
    x72 = x65 + x66
    x73 = x2 * x52 + x3 * (2 * x36 + x64)
    x74 = x60 * (x2 * x71 + x3 * (3 * x37 + 3 * x39 + x72 + x73))
    x75 = x3 * x8
    x76 = x57 * x75
    x77 = x57 * x8
    x78 = x70 * x77
    x79 = x62 * x78
    x80 = x76 + x79
    x81 = x13 + x63
    x82 = x2 * x81 + 2 * x3 * x34
    x83 = x2 * x72 + x3 * (3 * x45 + 3 * x47 + x82)
    x84 = x58 * x8
    x85 = x60 * x83
    x86 = -x68 - B[2]
    x87 = x58 * x75
    x88 = x84 * x86
    x89 = x69 * x88
    x90 = x87 + x89
    x91 = x70**2 * x77
    x92 = x76 + x91
    x93 = x2 * x73 + x3 * (3 * x49 + 3 * x51 + x82)
    x94 = x62 * x77
    x95 = x3 * (x78 + x94)
    x96 = x70 * x80
    x97 = x95 + x96
    x98 = x2 * x82 + x3 * (x14 + 3 * x63)
    x99 = x69 * x84
    x100 = x60 * x86
    x101 = x84 * x86**2
    x102 = x101 + x87
    x103 = x3 * (x88 + x99)
    x104 = x86 * x90
    x105 = x103 + x104
    x106 = -x61 - A[1]
    x107 = x44 * x60
    x108 = x106 * x94
    x109 = x108 + x76
    x110 = x56 * x60
    x111 = x106 * x78
    x112 = x111 + x76
    x113 = x106 * x80
    x114 = x113 + x95
    x115 = x106 * x77
    x116 = 2 * x76
    x117 = x106 * x92
    x118 = x116 * x70 + x117
    x119 = 3 * x76
    x120 = x119 + x91
    x121 = x3 * (x120 + 2 * x79)
    x122 = x106 * x97
    x123 = x121 + x122
    x124 = -x68 - A[2]
    x125 = x124 * x99
    x126 = x125 + x87
    x127 = x124 * x60
    x128 = x124 * x84
    x129 = x124 * x88
    x130 = x129 + x87
    x131 = x124 * x90
    x132 = x103 + x131
    x133 = 2 * x87
    x134 = x102 * x124
    x135 = x133 * x86 + x134
    x136 = 3 * x87
    x137 = x101 + x136
    x138 = x3 * (x137 + 2 * x89)
    x139 = x105 * x124
    x140 = x138 + x139
    x141 = x106**2 * x77
    x142 = x141 + x76
    x143 = x3 * (x115 + x94)
    x144 = x106 * x109
    x145 = x143 + x144
    x146 = x3 * (x115 + x78)
    x147 = x106 * x112
    x148 = x146 + x147
    x149 = x3 * (x108 + x111 + x119 + x79)
    x150 = x106 * x114
    x151 = x149 + x150
    x152 = 2 * x111
    x153 = x3 * (x120 + x152)
    x154 = x106 * x118
    x155 = x153 + x154
    x156 = x106 * x123
    x157 = 2 * x113
    x158 = x3 * (x118 + x157 + 3 * x95 + x96)
    x159 = x156 + x158
    x160 = x124**2 * x84
    x161 = x160 + x87
    x162 = x3 * (x128 + x99)
    x163 = x124 * x126
    x164 = x162 + x163
    x165 = x3 * (x128 + x88)
    x166 = x124 * x130
    x167 = x165 + x166
    x168 = x3 * (x125 + x129 + x136 + x89)
    x169 = x124 * x132
    x170 = x168 + x169
    x171 = 2 * x129
    x172 = x3 * (x137 + x171)
    x173 = x124 * x135
    x174 = x172 + x173
    x175 = x124 * x140
    x176 = 2 * x131
    x177 = x3 * (3 * x103 + x104 + x135 + x176)
    x178 = x175 + x177
    x179 = x106 * x116 + x106 * x142
    x180 = x119 + x141
    x181 = x106 * x145 + x3 * (2 * x108 + x180)
    x182 = x3 * (x152 + x180)
    x183 = x106 * x148
    x184 = x182 + x183
    x185 = x106 * x151
    x186 = x3 * (x145 + x148 + x157 + 2 * x95)
    x187 = x185 + x186
    x188 = x106 * x155 + x3 * (2 * x117 + 2 * x146 + 2 * x147 + 4 * x70 * x76)
    x189 = x106 * x159 + x3 * (2 * x121 + 2 * x122 + 2 * x149 + 2 * x150 + x155)
    x190 = x33 * x59
    x191 = numpy.pi * x0 * x57
    x192 = x191 * x33
    x193 = x124 * x133 + x124 * x161
    x194 = x136 + x160
    x195 = x124 * x164 + x3 * (2 * x125 + x194)
    x196 = x3 * (x171 + x194)
    x197 = x124 * x167
    x198 = x196 + x197
    x199 = x124 * x170
    x200 = x3 * (2 * x103 + x164 + x167 + x176)
    x201 = x199 + x200
    x202 = x124 * x174 + x3 * (2 * x134 + 2 * x165 + 2 * x166 + 4 * x86 * x87)
    x203 = x124 * x178 + x3 * (2 * x138 + 2 * x139 + 2 * x168 + 2 * x169 + x174)
    x204 = x106 * x179 + x3 * (x119 + 3 * x141)
    x205 = x106 * x181 + x3 * (3 * x143 + 3 * x144 + x179)
    x206 = x106 * x184 + x3 * (3 * x146 + 3 * x147 + x179)
    x207 = x59 * x7
    x208 = x207 * (x106 * x187 + x3 * (3 * x149 + 3 * x150 + x181 + x184))
    x209 = x207 * x5
    x210 = x207 * (x106 * x188 + x3 * (3 * x153 + 3 * x154 + 2 * x182 + 2 * x183))
    x211 = x124 * x207
    x212 = x191 * x7
    x213 = x106 * x212
    x214 = x124 * x193 + x3 * (x136 + 3 * x160)
    x215 = x124 * x195 + x3 * (3 * x162 + 3 * x163 + x193)
    x216 = x212 * x5
    x217 = x124 * x198 + x3 * (3 * x165 + 3 * x166 + x193)
    x218 = x212 * (x124 * x201 + x3 * (3 * x168 + 3 * x169 + x195 + x198))
    x219 = x212 * (x124 * x202 + x3 * (3 * x172 + 3 * x173 + 2 * x196 + 2 * x197))

    # 270 item(s)
    S = numpy.array(
        [
            x60 * (x2 * x44 + x3 * (3 * x24 + 3 * x31 + 2 * x53 + 2 * x55 + x56)),
            x62 * x67,
            x67 * x69,
            x70 * x74,
            x80 * x83 * x84,
            x69 * x70 * x85,
            x74 * x86,
            x62 * x85 * x86,
            x77 * x83 * x90,
            x84 * x92 * x93,
            x84 * x97 * x98,
            x92 * x98 * x99,
            x100 * x70 * x93,
            x80 * x88 * x98,
            x78 * x90 * x98,
            x102 * x77 * x93,
            x102 * x94 * x98,
            x105 * x77 * x98,
            x106 * x107,
            x109 * x56 * x84,
            x106 * x110 * x69,
            x112 * x71 * x84,
            x114 * x72 * x84,
            x112 * x72 * x99,
            x100 * x106 * x71,
            x109 * x72 * x88,
            x115 * x72 * x90,
            x118 * x73 * x84,
            x123 * x82 * x84,
            x118 * x82 * x99,
            x112 * x73 * x88,
            x114 * x82 * x88,
            x112 * x82 * x90,
            x102 * x115 * x73,
            x102 * x109 * x82,
            x105 * x115 * x82,
            x107 * x124,
            x110 * x124 * x62,
            x126 * x56 * x77,
            x127 * x70 * x71,
            x128 * x72 * x80,
            x126 * x72 * x78,
            x130 * x71 * x77,
            x130 * x72 * x94,
            x132 * x72 * x77,
            x128 * x73 * x92,
            x128 * x82 * x97,
            x126 * x82 * x92,
            x130 * x73 * x78,
            x130 * x80 * x82,
            x132 * x78 * x82,
            x135 * x73 * x77,
            x135 * x82 * x94,
            x140 * x77 * x82,
            x142 * x32 * x84,
            x145 * x43 * x84,
            x142 * x43 * x99,
            x148 * x54 * x84,
            x151 * x48 * x84,
            x148 * x48 * x99,
            x142 * x54 * x88,
            x145 * x48 * x88,
            x142 * x48 * x90,
            x155 * x52 * x84,
            x159 * x81 * x84,
            x155 * x81 * x99,
            x148 * x52 * x88,
            x151 * x81 * x88,
            x148 * x81 * x90,
            x102 * x142 * x52,
            x102 * x145 * x81,
            x105 * x142 * x81,
            x106 * x127 * x32,
            x109 * x128 * x43,
            x115 * x126 * x43,
            x112 * x128 * x54,
            x114 * x128 * x48,
            x112 * x126 * x48,
            x115 * x130 * x54,
            x109 * x130 * x48,
            x115 * x132 * x48,
            x118 * x128 * x52,
            x123 * x128 * x81,
            x118 * x126 * x81,
            x112 * x130 * x52,
            x114 * x130 * x81,
            x112 * x132 * x81,
            x115 * x135 * x52,
            x109 * x135 * x81,
            x115 * x140 * x81,
            x161 * x32 * x77,
            x161 * x43 * x94,
            x164 * x43 * x77,
            x161 * x54 * x78,
            x161 * x48 * x80,
            x164 * x48 * x78,
            x167 * x54 * x77,
            x167 * x48 * x94,
            x170 * x48 * x77,
            x161 * x52 * x92,
            x161 * x81 * x97,
            x164 * x81 * x92,
            x167 * x52 * x78,
            x167 * x80 * x81,
            x170 * x78 * x81,
            x174 * x52 * x77,
            x174 * x81 * x94,
            x178 * x77 * x81,
            x179 * x23 * x84,
            x181 * x30 * x84,
            x179 * x30 * x99,
            x184 * x38 * x84,
            x187 * x46 * x84,
            x184 * x46 * x99,
            x179 * x38 * x88,
            x181 * x46 * x88,
            x179 * x46 * x90,
            x188 * x50 * x84,
            x189 * x190,
            x188 * x190 * x69,
            x184 * x50 * x88,
            x187 * x190 * x86,
            x184 * x34 * x90,
            x102 * x179 * x50,
            x102 * x181 * x34,
            x105 * x179 * x34,
            x128 * x142 * x23,
            x128 * x145 * x30,
            x126 * x142 * x30,
            x128 * x148 * x38,
            x128 * x151 * x46,
            x126 * x148 * x46,
            x130 * x142 * x38,
            x130 * x145 * x46,
            x132 * x142 * x46,
            x128 * x155 * x50,
            x124 * x159 * x190,
            x126 * x155 * x34,
            x130 * x148 * x50,
            x130 * x151 * x34,
            x132 * x148 * x34,
            x135 * x142 * x50,
            x135 * x145 * x34,
            x140 * x142 * x34,
            x115 * x161 * x23,
            x109 * x161 * x30,
            x115 * x164 * x30,
            x112 * x161 * x38,
            x114 * x161 * x46,
            x112 * x164 * x46,
            x115 * x167 * x38,
            x109 * x167 * x46,
            x115 * x170 * x46,
            x118 * x161 * x50,
            x123 * x161 * x34,
            x118 * x164 * x34,
            x112 * x167 * x50,
            x114 * x167 * x34,
            x112 * x170 * x34,
            x115 * x174 * x50,
            x109 * x174 * x34,
            x106 * x178 * x192,
            x193 * x23 * x77,
            x193 * x30 * x94,
            x195 * x30 * x77,
            x193 * x38 * x78,
            x193 * x46 * x80,
            x195 * x46 * x78,
            x198 * x38 * x77,
            x198 * x46 * x94,
            x201 * x46 * x77,
            x193 * x50 * x92,
            x193 * x34 * x97,
            x195 * x34 * x92,
            x198 * x50 * x78,
            x198 * x34 * x80,
            x192 * x201 * x70,
            x202 * x50 * x77,
            x192 * x202 * x62,
            x192 * x203,
            x204 * x21 * x84,
            x205 * x28 * x84,
            x204 * x28 * x99,
            x19 * x206 * x84,
            x208 * x5,
            x206 * x209 * x69,
            x19 * x204 * x88,
            x205 * x209 * x86,
            x10 * x204 * x90,
            x210 * x4,
            x207
            * (x106 * x189 + x3 * (3 * x156 + 3 * x158 + 2 * x185 + 2 * x186 + x188)),
            x210 * x69,
            x206 * x207 * x4 * x86,
            x208 * x86,
            x206 * x9 * x90,
            x102 * x17 * x204,
            x102 * x205 * x9,
            x105 * x204 * x9,
            x128 * x179 * x21,
            x128 * x181 * x28,
            x126 * x179 * x28,
            x128 * x184 * x19,
            x187 * x211 * x5,
            x10 * x126 * x184,
            x130 * x179 * x19,
            x10 * x130 * x181,
            x10 * x132 * x179,
            x188 * x211 * x4,
            x189 * x211,
            x126 * x188 * x9,
            x130 * x17 * x184,
            x130 * x187 * x9,
            x132 * x184 * x9,
            x135 * x17 * x179,
            x135 * x181 * x9,
            x140 * x179 * x9,
            x142 * x161 * x21,
            x145 * x161 * x28,
            x142 * x164 * x28,
            x148 * x161 * x19,
            x10 * x151 * x161,
            x10 * x148 * x164,
            x142 * x167 * x19,
            x10 * x145 * x167,
            x10 * x142 * x170,
            x155 * x161 * x17,
            x159 * x161 * x9,
            x155 * x164 * x9,
            x148 * x167 * x17,
            x151 * x167 * x9,
            x148 * x170 * x9,
            x142 * x17 * x174,
            x145 * x174 * x9,
            x142 * x178 * x9,
            x115 * x193 * x21,
            x109 * x193 * x28,
            x115 * x195 * x28,
            x112 * x19 * x193,
            x10 * x114 * x193,
            x10 * x112 * x195,
            x115 * x19 * x198,
            x10 * x109 * x198,
            x201 * x213 * x5,
            x118 * x17 * x193,
            x123 * x193 * x9,
            x118 * x195 * x9,
            x112 * x17 * x198,
            x114 * x198 * x9,
            x112 * x201 * x9,
            x202 * x213 * x4,
            x109 * x202 * x9,
            x203 * x213,
            x21 * x214 * x77,
            x214 * x28 * x94,
            x215 * x28 * x77,
            x19 * x214 * x78,
            x10 * x214 * x80,
            x215 * x216 * x70,
            x19 * x217 * x77,
            x216 * x217 * x62,
            x218 * x5,
            x17 * x214 * x92,
            x214 * x9 * x97,
            x215 * x9 * x92,
            x212 * x217 * x4 * x70,
            x217 * x80 * x9,
            x218 * x70,
            x219 * x4,
            x219 * x62,
            x212
            * (x124 * x203 + x3 * (3 * x175 + 3 * x177 + 2 * x199 + 2 * x200 + x202)),
        ]
    )
    return S


def dipole3d_43(a, A, b, B, C):
    """Cartesian 3D (gf) dipole moment integral.
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
    x20 = x4**2 * x8
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
    x71 = -x70 - C[1]
    x72 = x2 * x60
    x73 = x2 * x8
    x74 = x3 * (x73 + x9)
    x75 = x14 + x52
    x76 = x2 * x75
    x77 = x3 * (4 * x18 + 2 * x40 + 2 * x74 + 2 * x76)
    x78 = x69 * (x2 * x65 + x3 * (3 * x41 + 3 * x42 + 3 * x72 + 3 * x77))
    x79 = -x0 * (a * A[2] + b * B[2])
    x80 = -x79 - C[2]
    x81 = -x70 - B[1]
    x82 = x61 + x63
    x83 = x74 + x76
    x84 = x3 * (x11 + x73)
    x85 = x14 + x53
    x86 = x2 * x85
    x87 = x84 + x86
    x88 = x3 * (2 * x12 + x45 + x83 + x87)
    x89 = x54 + x56
    x90 = x2 * x89
    x91 = x72 + x77
    x92 = x69 * (x2 * x82 + x3 * (x50 + 2 * x88 + 2 * x90 + x91))
    x93 = x3 * x7
    x94 = x66 * x93
    x95 = x66 * x7
    x96 = x81 * x95
    x97 = x71 * x96
    x98 = x94 + x97
    x99 = x2**2 * x8
    x100 = x25 + x99
    x101 = x3 * (x100 + x57)
    x102 = x2 * x83
    x103 = x2 * x91 + x3 * (2 * x101 + 2 * x102 + x64)
    x104 = x67 * x7
    x105 = x103 * x69
    x106 = -x79 - B[2]
    x107 = x67 * x93
    x108 = x104 * x106
    x109 = x108 * x80
    x110 = x107 + x109
    x111 = x81**2 * x95
    x112 = x111 + x94
    x113 = x88 + x90
    x114 = x101 + x102
    x115 = x2 * x87 + x3 * (x100 + 2 * x53)
    x116 = x113 * x2 + x3 * (x114 + x115 + 3 * x54 + 3 * x56)
    x117 = x71 * x95
    x118 = x3 * (x117 + x96)
    x119 = x81 * x98
    x120 = x118 + x119
    x121 = x14 + x99
    x122 = x121 * x2 + 2 * x14 * x2
    x123 = x114 * x2 + x3 * (x122 + 3 * x74 + 3 * x76)
    x124 = x104 * x80
    x125 = x106 * x69
    x126 = x104 * x106**2
    x127 = x107 + x126
    x128 = x3 * (x108 + x124)
    x129 = x106 * x110
    x130 = x128 + x129
    x131 = x81 * x94
    x132 = 2 * x131
    x133 = x112 * x81
    x134 = x132 + x133
    x135 = x115 * x2 + x3 * (x122 + 3 * x84 + 3 * x86)
    x136 = 3 * x94
    x137 = x111 + x136
    x138 = x3 * (x137 + 2 * x97)
    x139 = x120 * x81
    x140 = x138 + x139
    x141 = x122 * x2 + x3 * (x25 + 3 * x99)
    x142 = x106 * x107
    x143 = 2 * x142
    x144 = x106 * x127
    x145 = x143 + x144
    x146 = 3 * x107
    x147 = x126 + x146
    x148 = x3 * (2 * x109 + x147)
    x149 = x106 * x130
    x150 = x148 + x149
    x151 = -x70 - A[1]
    x152 = x51 * x69
    x153 = x117 * x151
    x154 = x153 + x94
    x155 = x65 * x69
    x156 = x151 * x96
    x157 = x156 + x94
    x158 = x151 * x98
    x159 = x118 + x158
    x160 = x151 * x95
    x161 = x112 * x151
    x162 = x132 + x161
    x163 = x120 * x151
    x164 = x138 + x163
    x165 = x3 * (3 * x111 + x136)
    x166 = x134 * x151
    x167 = x165 + x166
    x168 = 3 * x118
    x169 = x3 * (3 * x119 + x134 + x168)
    x170 = x140 * x151
    x171 = x169 + x170
    x172 = -x79 - A[2]
    x173 = x124 * x172
    x174 = x107 + x173
    x175 = x172 * x69
    x176 = x104 * x172
    x177 = x108 * x172
    x178 = x107 + x177
    x179 = x110 * x172
    x180 = x128 + x179
    x181 = x127 * x172
    x182 = x143 + x181
    x183 = x130 * x172
    x184 = x148 + x183
    x185 = x3 * (3 * x126 + x146)
    x186 = x145 * x172
    x187 = x185 + x186
    x188 = 3 * x128
    x189 = x3 * (3 * x129 + x145 + x188)
    x190 = x150 * x172
    x191 = x189 + x190
    x192 = x151**2 * x95
    x193 = x192 + x94
    x194 = x3 * (x117 + x160)
    x195 = x151 * x154
    x196 = x194 + x195
    x197 = x3 * (x160 + x96)
    x198 = x151 * x157
    x199 = x197 + x198
    x200 = x3 * (x136 + x153 + x156 + x97)
    x201 = x151 * x159
    x202 = x200 + x201
    x203 = 2 * x156
    x204 = x3 * (x137 + x203)
    x205 = x151 * x162
    x206 = x204 + x205
    x207 = x151 * x164
    x208 = 2 * x158
    x209 = x3 * (x119 + x162 + x168 + x208)
    x210 = x207 + x209
    x211 = x3 * (8 * x131 + x133 + 3 * x161)
    x212 = x151 * x167
    x213 = x211 + x212
    x214 = x151 * x171
    x215 = x3 * (4 * x138 + x139 + 3 * x163 + x167)
    x216 = x214 + x215
    x217 = x104 * x172**2
    x218 = x107 + x217
    x219 = x3 * (x124 + x176)
    x220 = x172 * x174
    x221 = x219 + x220
    x222 = x3 * (x108 + x176)
    x223 = x172 * x178
    x224 = x222 + x223
    x225 = x3 * (x109 + x146 + x173 + x177)
    x226 = x172 * x180
    x227 = x225 + x226
    x228 = 2 * x177
    x229 = x3 * (x147 + x228)
    x230 = x172 * x182
    x231 = x229 + x230
    x232 = x172 * x184
    x233 = 2 * x179
    x234 = x3 * (x129 + x182 + x188 + x233)
    x235 = x232 + x234
    x236 = x3 * (8 * x142 + x144 + 3 * x181)
    x237 = x172 * x187
    x238 = x236 + x237
    x239 = x172 * x191
    x240 = x3 * (4 * x148 + x149 + 3 * x183 + x187)
    x241 = x239 + x240
    x242 = x151 * x193 + 2 * x151 * x94
    x243 = x136 + x192
    x244 = x151 * x196 + x3 * (2 * x153 + x243)
    x245 = x3 * (x203 + x243)
    x246 = x151 * x199
    x247 = x245 + x246
    x248 = x151 * x202
    x249 = x3 * (2 * x118 + x196 + x199 + x208)
    x250 = x248 + x249
    x251 = x151 * x206
    x252 = x3 * (4 * x131 + 2 * x161 + 2 * x197 + 2 * x198)
    x253 = x251 + x252
    x254 = x151 * x210
    x255 = x3 * (2 * x138 + 2 * x163 + 2 * x200 + 2 * x201 + x206)
    x256 = x254 + x255
    x257 = 3 * x204 + 3 * x205
    x258 = x151 * x213 + x3 * (2 * x165 + 2 * x166 + x257)
    x259 = 3 * x207 + 3 * x209
    x260 = x151 * x216 + x3 * (2 * x169 + 2 * x170 + x213 + x259)
    x261 = x6 * x68
    x262 = x260 * x261
    x263 = x2 * x261
    x264 = numpy.pi * x0 * x6 * x66
    x265 = x2 * x264
    x266 = 2 * x107 * x172 + x172 * x218
    x267 = x146 + x217
    x268 = x172 * x221 + x3 * (2 * x173 + x267)
    x269 = x3 * (x228 + x267)
    x270 = x172 * x224
    x271 = x269 + x270
    x272 = x172 * x227
    x273 = x3 * (2 * x128 + x221 + x224 + x233)
    x274 = x272 + x273
    x275 = x172 * x231
    x276 = x3 * (4 * x142 + 2 * x181 + 2 * x222 + 2 * x223)
    x277 = x275 + x276
    x278 = x172 * x235
    x279 = x3 * (2 * x148 + 2 * x183 + 2 * x225 + 2 * x226 + x231)
    x280 = x278 + x279
    x281 = 3 * x229 + 3 * x230
    x282 = x172 * x238 + x3 * (2 * x185 + 2 * x186 + x281)
    x283 = 3 * x232 + 3 * x234
    x284 = x172 * x241 + x3 * (2 * x189 + 2 * x190 + x238 + x283)
    x285 = x264 * x284
    x286 = x151 * x242 + x3 * (x136 + 3 * x192)
    x287 = x151 * x244 + x3 * (3 * x194 + 3 * x195 + x242)
    x288 = x151 * x247 + x3 * (3 * x197 + 3 * x198 + x242)
    x289 = x151 * x250 + x3 * (3 * x200 + 3 * x201 + x244 + x247)
    x290 = x151 * x253 + x3 * (2 * x245 + 2 * x246 + x257)
    x291 = x261 * (x151 * x256 + x3 * (2 * x248 + 2 * x249 + x253 + x259))
    x292 = x261 * x4
    x293 = x261 * (x151 * x258 + x3 * (3 * x211 + 3 * x212 + 3 * x251 + 3 * x252))
    x294 = x10 * x261
    x295 = x151 * x264
    x296 = x172 * x266 + x3 * (x146 + 3 * x217)
    x297 = x172 * x268 + x3 * (3 * x219 + 3 * x220 + x266)
    x298 = x172 * x271 + x3 * (3 * x222 + 3 * x223 + x266)
    x299 = x172 * x274 + x3 * (3 * x225 + 3 * x226 + x268 + x271)
    x300 = x264 * x4
    x301 = x172 * x277 + x3 * (2 * x269 + 2 * x270 + x281)
    x302 = x264 * (x172 * x280 + x3 * (2 * x272 + 2 * x273 + x277 + x283))
    x303 = x264 * (x172 * x282 + x3 * (3 * x236 + 3 * x237 + 3 * x275 + 3 * x276))

    # 450 item(s)
    S = numpy.array(
        [
            x69 * (x2 * x51 + x3 * (3 * x33 + 3 * x38 + 3 * x61 + 3 * x63 + x65)),
            x71 * x78,
            x78 * x80,
            x81 * x92,
            x103 * x104 * x98,
            x105 * x80 * x81,
            x106 * x92,
            x105 * x106 * x71,
            x103 * x110 * x95,
            x104 * x112 * x116,
            x104 * x120 * x123,
            x112 * x123 * x124,
            x116 * x125 * x81,
            x108 * x123 * x98,
            x110 * x123 * x96,
            x116 * x127 * x95,
            x117 * x123 * x127,
            x123 * x130 * x95,
            x104 * x134 * x135,
            x104 * x140 * x141,
            x124 * x134 * x141,
            x108 * x112 * x135,
            x108 * x120 * x141,
            x110 * x112 * x141,
            x127 * x135 * x96,
            x127 * x141 * x98,
            x130 * x141 * x96,
            x135 * x145 * x95,
            x117 * x141 * x145,
            x141 * x150 * x95,
            x151 * x152,
            x104 * x154 * x65,
            x151 * x155 * x80,
            x104 * x157 * x82,
            x104 * x159 * x91,
            x124 * x157 * x91,
            x125 * x151 * x82,
            x108 * x154 * x91,
            x110 * x160 * x91,
            x104 * x113 * x162,
            x104 * x114 * x164,
            x114 * x124 * x162,
            x108 * x113 * x157,
            x108 * x114 * x159,
            x110 * x114 * x157,
            x113 * x127 * x160,
            x114 * x127 * x154,
            x114 * x130 * x160,
            x104 * x115 * x167,
            x104 * x122 * x171,
            x122 * x124 * x167,
            x108 * x115 * x162,
            x108 * x122 * x164,
            x110 * x122 * x162,
            x115 * x127 * x157,
            x122 * x127 * x159,
            x122 * x130 * x157,
            x115 * x145 * x160,
            x122 * x145 * x154,
            x122 * x150 * x160,
            x152 * x172,
            x155 * x172 * x71,
            x174 * x65 * x95,
            x175 * x81 * x82,
            x176 * x91 * x98,
            x174 * x91 * x96,
            x178 * x82 * x95,
            x117 * x178 * x91,
            x180 * x91 * x95,
            x112 * x113 * x176,
            x114 * x120 * x176,
            x112 * x114 * x174,
            x113 * x178 * x96,
            x114 * x178 * x98,
            x114 * x180 * x96,
            x113 * x182 * x95,
            x114 * x117 * x182,
            x114 * x184 * x95,
            x115 * x134 * x176,
            x122 * x140 * x176,
            x122 * x134 * x174,
            x112 * x115 * x178,
            x120 * x122 * x178,
            x112 * x122 * x180,
            x115 * x182 * x96,
            x122 * x182 * x98,
            x122 * x184 * x96,
            x115 * x187 * x95,
            x117 * x122 * x187,
            x122 * x191 * x95,
            x104 * x193 * x39,
            x104 * x196 * x43,
            x124 * x193 * x43,
            x104 * x199 * x62,
            x104 * x202 * x60,
            x124 * x199 * x60,
            x108 * x193 * x62,
            x108 * x196 * x60,
            x110 * x193 * x60,
            x104 * x206 * x89,
            x104 * x210 * x83,
            x124 * x206 * x83,
            x108 * x199 * x89,
            x108 * x202 * x83,
            x110 * x199 * x83,
            x127 * x193 * x89,
            x127 * x196 * x83,
            x130 * x193 * x83,
            x104 * x213 * x87,
            x104 * x121 * x216,
            x121 * x124 * x213,
            x108 * x206 * x87,
            x108 * x121 * x210,
            x110 * x121 * x206,
            x127 * x199 * x87,
            x121 * x127 * x202,
            x121 * x130 * x199,
            x145 * x193 * x87,
            x121 * x145 * x196,
            x121 * x150 * x193,
            x151 * x175 * x39,
            x154 * x176 * x43,
            x160 * x174 * x43,
            x157 * x176 * x62,
            x159 * x176 * x60,
            x157 * x174 * x60,
            x160 * x178 * x62,
            x154 * x178 * x60,
            x160 * x180 * x60,
            x162 * x176 * x89,
            x164 * x176 * x83,
            x162 * x174 * x83,
            x157 * x178 * x89,
            x159 * x178 * x83,
            x157 * x180 * x83,
            x160 * x182 * x89,
            x154 * x182 * x83,
            x160 * x184 * x83,
            x167 * x176 * x87,
            x121 * x171 * x176,
            x121 * x167 * x174,
            x162 * x178 * x87,
            x121 * x164 * x178,
            x121 * x162 * x180,
            x157 * x182 * x87,
            x121 * x159 * x182,
            x121 * x157 * x184,
            x160 * x187 * x87,
            x121 * x154 * x187,
            x121 * x160 * x191,
            x218 * x39 * x95,
            x117 * x218 * x43,
            x221 * x43 * x95,
            x218 * x62 * x96,
            x218 * x60 * x98,
            x221 * x60 * x96,
            x224 * x62 * x95,
            x117 * x224 * x60,
            x227 * x60 * x95,
            x112 * x218 * x89,
            x120 * x218 * x83,
            x112 * x221 * x83,
            x224 * x89 * x96,
            x224 * x83 * x98,
            x227 * x83 * x96,
            x231 * x89 * x95,
            x117 * x231 * x83,
            x235 * x83 * x95,
            x134 * x218 * x87,
            x121 * x140 * x218,
            x121 * x134 * x221,
            x112 * x224 * x87,
            x120 * x121 * x224,
            x112 * x121 * x227,
            x231 * x87 * x96,
            x121 * x231 * x98,
            x121 * x235 * x96,
            x238 * x87 * x95,
            x117 * x121 * x238,
            x121 * x241 * x95,
            x104 * x242 * x32,
            x104 * x244 * x37,
            x124 * x242 * x37,
            x104 * x247 * x48,
            x104 * x250 * x46,
            x124 * x247 * x46,
            x108 * x242 * x48,
            x108 * x244 * x46,
            x110 * x242 * x46,
            x104 * x253 * x55,
            x104 * x256 * x75,
            x124 * x253 * x75,
            x108 * x247 * x55,
            x108 * x250 * x75,
            x110 * x247 * x75,
            x127 * x242 * x55,
            x127 * x244 * x75,
            x130 * x242 * x75,
            x104 * x258 * x85,
            x2 * x262,
            x258 * x263 * x80,
            x108 * x253 * x85,
            x106 * x256 * x263,
            x110 * x253 * x73,
            x127 * x247 * x85,
            x127 * x250 * x73,
            x130 * x247 * x73,
            x145 * x242 * x85,
            x145 * x244 * x73,
            x150 * x242 * x73,
            x176 * x193 * x32,
            x176 * x196 * x37,
            x174 * x193 * x37,
            x176 * x199 * x48,
            x176 * x202 * x46,
            x174 * x199 * x46,
            x178 * x193 * x48,
            x178 * x196 * x46,
            x180 * x193 * x46,
            x176 * x206 * x55,
            x176 * x210 * x75,
            x174 * x206 * x75,
            x178 * x199 * x55,
            x178 * x202 * x75,
            x180 * x199 * x75,
            x182 * x193 * x55,
            x182 * x196 * x75,
            x184 * x193 * x75,
            x176 * x213 * x85,
            x172 * x216 * x263,
            x174 * x213 * x73,
            x178 * x206 * x85,
            x178 * x210 * x73,
            x180 * x206 * x73,
            x182 * x199 * x85,
            x182 * x202 * x73,
            x184 * x199 * x73,
            x187 * x193 * x85,
            x187 * x196 * x73,
            x191 * x193 * x73,
            x160 * x218 * x32,
            x154 * x218 * x37,
            x160 * x221 * x37,
            x157 * x218 * x48,
            x159 * x218 * x46,
            x157 * x221 * x46,
            x160 * x224 * x48,
            x154 * x224 * x46,
            x160 * x227 * x46,
            x162 * x218 * x55,
            x164 * x218 * x75,
            x162 * x221 * x75,
            x157 * x224 * x55,
            x159 * x224 * x75,
            x157 * x227 * x75,
            x160 * x231 * x55,
            x154 * x231 * x75,
            x160 * x235 * x75,
            x167 * x218 * x85,
            x171 * x218 * x73,
            x167 * x221 * x73,
            x162 * x224 * x85,
            x164 * x224 * x73,
            x162 * x227 * x73,
            x157 * x231 * x85,
            x159 * x231 * x73,
            x157 * x235 * x73,
            x160 * x238 * x85,
            x154 * x238 * x73,
            x151 * x241 * x265,
            x266 * x32 * x95,
            x117 * x266 * x37,
            x268 * x37 * x95,
            x266 * x48 * x96,
            x266 * x46 * x98,
            x268 * x46 * x96,
            x271 * x48 * x95,
            x117 * x271 * x46,
            x274 * x46 * x95,
            x112 * x266 * x55,
            x120 * x266 * x75,
            x112 * x268 * x75,
            x271 * x55 * x96,
            x271 * x75 * x98,
            x274 * x75 * x96,
            x277 * x55 * x95,
            x117 * x277 * x75,
            x280 * x75 * x95,
            x134 * x266 * x85,
            x140 * x266 * x73,
            x134 * x268 * x73,
            x112 * x271 * x85,
            x120 * x271 * x73,
            x112 * x274 * x73,
            x277 * x85 * x96,
            x277 * x73 * x98,
            x265 * x280 * x81,
            x282 * x85 * x95,
            x265 * x282 * x71,
            x2 * x285,
            x104 * x286 * x30,
            x104 * x23 * x287,
            x124 * x23 * x286,
            x104 * x28 * x288,
            x104 * x21 * x289,
            x124 * x21 * x288,
            x108 * x28 * x286,
            x108 * x21 * x287,
            x110 * x21 * x286,
            x104 * x16 * x290,
            x291 * x4,
            x290 * x292 * x80,
            x108 * x16 * x288,
            x106 * x289 * x292,
            x110 * x288 * x9,
            x127 * x16 * x286,
            x127 * x287 * x9,
            x130 * x286 * x9,
            x10 * x293,
            x261
            * (x151 * x260 + x3 * (3 * x214 + 3 * x215 + 3 * x254 + 3 * x255 + x258)),
            x293 * x80,
            x106 * x290 * x294,
            x106 * x291,
            x110 * x290 * x8,
            x11 * x127 * x288,
            x127 * x289 * x8,
            x130 * x288 * x8,
            x11 * x145 * x286,
            x145 * x287 * x8,
            x150 * x286 * x8,
            x176 * x242 * x30,
            x176 * x23 * x244,
            x174 * x23 * x242,
            x176 * x247 * x28,
            x176 * x21 * x250,
            x174 * x21 * x247,
            x178 * x242 * x28,
            x178 * x21 * x244,
            x180 * x21 * x242,
            x16 * x176 * x253,
            x172 * x256 * x292,
            x174 * x253 * x9,
            x16 * x178 * x247,
            x178 * x250 * x9,
            x180 * x247 * x9,
            x16 * x182 * x242,
            x182 * x244 * x9,
            x184 * x242 * x9,
            x172 * x258 * x294,
            x172 * x262,
            x174 * x258 * x8,
            x11 * x178 * x253,
            x178 * x256 * x8,
            x180 * x253 * x8,
            x11 * x182 * x247,
            x182 * x250 * x8,
            x184 * x247 * x8,
            x11 * x187 * x242,
            x187 * x244 * x8,
            x191 * x242 * x8,
            x193 * x218 * x30,
            x196 * x218 * x23,
            x193 * x221 * x23,
            x199 * x218 * x28,
            x202 * x21 * x218,
            x199 * x21 * x221,
            x193 * x224 * x28,
            x196 * x21 * x224,
            x193 * x21 * x227,
            x16 * x206 * x218,
            x210 * x218 * x9,
            x206 * x221 * x9,
            x16 * x199 * x224,
            x202 * x224 * x9,
            x199 * x227 * x9,
            x16 * x193 * x231,
            x196 * x231 * x9,
            x193 * x235 * x9,
            x11 * x213 * x218,
            x216 * x218 * x8,
            x213 * x221 * x8,
            x11 * x206 * x224,
            x210 * x224 * x8,
            x206 * x227 * x8,
            x11 * x199 * x231,
            x202 * x231 * x8,
            x199 * x235 * x8,
            x11 * x193 * x238,
            x196 * x238 * x8,
            x193 * x241 * x8,
            x160 * x266 * x30,
            x154 * x23 * x266,
            x160 * x23 * x268,
            x157 * x266 * x28,
            x159 * x21 * x266,
            x157 * x21 * x268,
            x160 * x271 * x28,
            x154 * x21 * x271,
            x160 * x21 * x274,
            x16 * x162 * x266,
            x164 * x266 * x9,
            x162 * x268 * x9,
            x157 * x16 * x271,
            x159 * x271 * x9,
            x157 * x274 * x9,
            x16 * x160 * x277,
            x154 * x277 * x9,
            x280 * x295 * x4,
            x11 * x167 * x266,
            x171 * x266 * x8,
            x167 * x268 * x8,
            x11 * x162 * x271,
            x164 * x271 * x8,
            x162 * x274 * x8,
            x11 * x157 * x277,
            x159 * x277 * x8,
            x157 * x280 * x8,
            x10 * x282 * x295,
            x154 * x282 * x8,
            x151 * x285,
            x296 * x30 * x95,
            x117 * x23 * x296,
            x23 * x297 * x95,
            x28 * x296 * x96,
            x21 * x296 * x98,
            x21 * x297 * x96,
            x28 * x298 * x95,
            x117 * x21 * x298,
            x21 * x299 * x95,
            x112 * x16 * x296,
            x120 * x296 * x9,
            x112 * x297 * x9,
            x16 * x298 * x96,
            x298 * x9 * x98,
            x299 * x300 * x81,
            x16 * x301 * x95,
            x300 * x301 * x71,
            x302 * x4,
            x11 * x134 * x296,
            x140 * x296 * x8,
            x134 * x297 * x8,
            x11 * x112 * x298,
            x120 * x298 * x8,
            x112 * x299 * x8,
            x10 * x264 * x301 * x81,
            x301 * x8 * x98,
            x302 * x81,
            x10 * x303,
            x303 * x71,
            x264
            * (x172 * x284 + x3 * (3 * x239 + 3 * x240 + 3 * x278 + 3 * x279 + x282)),
        ]
    )
    return S


def dipole3d_44(a, A, b, B, C):
    """Cartesian 3D (gg) dipole moment integral.
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
    x12 = x5**2 * x9
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
    x78 = -x77 - C[1]
    x79 = x2 * x61
    x80 = x10 * x2
    x81 = 2 * x80
    x82 = x3 * (x15 + x81)
    x83 = x2 * x64
    x84 = 3 * x82 + 3 * x83
    x85 = x3 * (2 * x24 + 2 * x49 + x84)
    x86 = x76 * (x2 * x72 + x3 * (3 * x54 + 3 * x55 + 4 * x79 + 4 * x85))
    x87 = -x0 * (a * A[2] + b * B[2])
    x88 = -x87 - C[2]
    x89 = -x77 - B[1]
    x90 = x69 + x71
    x91 = x18 * x2
    x92 = x3 * (x11 + x14 + x80 + x91)
    x93 = x19 + x62
    x94 = x2 * x93
    x95 = x82 + x83
    x96 = x3 * (2 * x16 + 2 * x48 + 2 * x92 + 2 * x94 + x95)
    x97 = x65 + x67
    x98 = x2 * x97
    x99 = x79 + x85
    x100 = x76 * (x2 * x90 + x3 * (3 * x51 + 3 * x53 + 3 * x96 + 3 * x98 + x99))
    x101 = x3 * x8
    x102 = x101 * x73
    x103 = x73 * x8
    x104 = x103 * x89
    x105 = x104 * x78
    x106 = x102 + x105
    x107 = x2 * x95
    x108 = x2 * x9
    x109 = x3 * (x10 + x108)
    x110 = x13 + x80
    x111 = x110 * x2
    x112 = x3 * (2 * x109 + 2 * x111 + 4 * x25 + 2 * x58)
    x113 = x2 * x99 + x3 * (3 * x107 + 3 * x112 + 3 * x59 + 3 * x60)
    x114 = x74 * x8
    x115 = x113 * x76
    x116 = -x87 - B[2]
    x117 = x101 * x74
    x118 = x114 * x116
    x119 = x118 * x88
    x120 = x117 + x119
    x121 = x103 * x89**2
    x122 = x102 + x121
    x123 = x96 + x98
    x124 = x109 + x111
    x125 = x3 * (x108 + x18)
    x126 = x13 + x91
    x127 = x126 * x2
    x128 = x125 + x127
    x129 = x3 * (x124 + x128 + 2 * x19 + x63)
    x130 = x92 + x94
    x131 = x130 * x2
    x132 = x107 + x112
    x133 = x123 * x2 + x3 * (2 * x129 + 2 * x131 + x132 + x68)
    x134 = x103 * x78
    x135 = x3 * (x104 + x134)
    x136 = x106 * x89
    x137 = x135 + x136
    x138 = x2**2 * x9
    x139 = x138 + x14
    x140 = x3 * (x139 + x81)
    x141 = x124 * x2
    x142 = x132 * x2 + x3 * (2 * x140 + 2 * x141 + x84)
    x143 = x114 * x88
    x144 = x116 * x76
    x145 = x114 * x116**2
    x146 = x117 + x145
    x147 = x3 * (x118 + x143)
    x148 = x116 * x120
    x149 = x147 + x148
    x150 = x102 * x89
    x151 = 2 * x150
    x152 = x122 * x89
    x153 = x151 + x152
    x154 = x129 + x131
    x155 = x140 + x141
    x156 = x128 * x2 + x3 * (x139 + 2 * x91)
    x157 = x154 * x2 + x3 * (x155 + x156 + 3 * x92 + 3 * x94)
    x158 = 3 * x102
    x159 = x121 + x158
    x160 = x3 * (2 * x105 + x159)
    x161 = x137 * x89
    x162 = x160 + x161
    x163 = x13 + x138
    x164 = 2 * x13 * x2 + x163 * x2
    x165 = x155 * x2 + x3 * (3 * x109 + 3 * x111 + x164)
    x166 = x116 * x117
    x167 = 2 * x166
    x168 = x116 * x146
    x169 = x167 + x168
    x170 = 3 * x117
    x171 = x145 + x170
    x172 = x3 * (2 * x119 + x171)
    x173 = x116 * x149
    x174 = x172 + x173
    x175 = x3 * (3 * x121 + x158)
    x176 = x153 * x89
    x177 = x175 + x176
    x178 = x156 * x2 + x3 * (3 * x125 + 3 * x127 + x164)
    x179 = 3 * x135
    x180 = x3 * (3 * x136 + x153 + x179)
    x181 = x162 * x89
    x182 = x180 + x181
    x183 = x164 * x2 + x3 * (3 * x138 + x14)
    x184 = x3 * (3 * x145 + x170)
    x185 = x116 * x169
    x186 = x184 + x185
    x187 = 3 * x147
    x188 = x3 * (3 * x148 + x169 + x187)
    x189 = x116 * x174
    x190 = x188 + x189
    x191 = -x77 - A[1]
    x192 = x57 * x76
    x193 = x134 * x191
    x194 = x102 + x193
    x195 = x72 * x76
    x196 = x104 * x191
    x197 = x102 + x196
    x198 = x106 * x191
    x199 = x135 + x198
    x200 = x103 * x191
    x201 = x122 * x191
    x202 = x151 + x201
    x203 = x137 * x191
    x204 = x160 + x203
    x205 = x153 * x191
    x206 = x175 + x205
    x207 = x162 * x191
    x208 = x180 + x207
    x209 = 8 * x150
    x210 = x3 * (4 * x152 + x209)
    x211 = x177 * x191
    x212 = x210 + x211
    x213 = 4 * x160
    x214 = x3 * (4 * x161 + x177 + x213)
    x215 = x182 * x191
    x216 = x214 + x215
    x217 = -x87 - A[2]
    x218 = x143 * x217
    x219 = x117 + x218
    x220 = x217 * x76
    x221 = x114 * x217
    x222 = x118 * x217
    x223 = x117 + x222
    x224 = x120 * x217
    x225 = x147 + x224
    x226 = x146 * x217
    x227 = x167 + x226
    x228 = x149 * x217
    x229 = x172 + x228
    x230 = x169 * x217
    x231 = x184 + x230
    x232 = x174 * x217
    x233 = x188 + x232
    x234 = 8 * x166
    x235 = x3 * (4 * x168 + x234)
    x236 = x186 * x217
    x237 = x235 + x236
    x238 = 4 * x172
    x239 = x3 * (4 * x173 + x186 + x238)
    x240 = x190 * x217
    x241 = x239 + x240
    x242 = x103 * x191**2
    x243 = x102 + x242
    x244 = x3 * (x134 + x200)
    x245 = x191 * x194
    x246 = x244 + x245
    x247 = x3 * (x104 + x200)
    x248 = x191 * x197
    x249 = x247 + x248
    x250 = x3 * (x105 + x158 + x193 + x196)
    x251 = x191 * x199
    x252 = x250 + x251
    x253 = 2 * x196
    x254 = x3 * (x159 + x253)
    x255 = x191 * x202
    x256 = x254 + x255
    x257 = x191 * x204
    x258 = 2 * x198
    x259 = x3 * (x136 + x179 + x202 + x258)
    x260 = x257 + x259
    x261 = x3 * (x152 + 3 * x201 + x209)
    x262 = x191 * x206
    x263 = x261 + x262
    x264 = x191 * x208
    x265 = x3 * (x161 + 3 * x203 + x206 + x213)
    x266 = x264 + x265
    x267 = x3 * (5 * x175 + x176 + 4 * x205)
    x268 = x191 * x212
    x269 = x267 + x268
    x270 = x191 * x216
    x271 = x3 * (5 * x180 + x181 + 4 * x207 + x212)
    x272 = x270 + x271
    x273 = x114 * x217**2
    x274 = x117 + x273
    x275 = x3 * (x143 + x221)
    x276 = x217 * x219
    x277 = x275 + x276
    x278 = x3 * (x118 + x221)
    x279 = x217 * x223
    x280 = x278 + x279
    x281 = x3 * (x119 + x170 + x218 + x222)
    x282 = x217 * x225
    x283 = x281 + x282
    x284 = 2 * x222
    x285 = x3 * (x171 + x284)
    x286 = x217 * x227
    x287 = x285 + x286
    x288 = x217 * x229
    x289 = 2 * x224
    x290 = x3 * (x148 + x187 + x227 + x289)
    x291 = x288 + x290
    x292 = x3 * (x168 + 3 * x226 + x234)
    x293 = x217 * x231
    x294 = x292 + x293
    x295 = x217 * x233
    x296 = x3 * (x173 + 3 * x228 + x231 + x238)
    x297 = x295 + x296
    x298 = x3 * (5 * x184 + x185 + 4 * x230)
    x299 = x217 * x237
    x300 = x298 + x299
    x301 = x217 * x241
    x302 = x3 * (5 * x188 + x189 + 4 * x232 + x237)
    x303 = x301 + x302
    x304 = 2 * x102 * x191 + x191 * x243
    x305 = x158 + x242
    x306 = x191 * x246 + x3 * (2 * x193 + x305)
    x307 = x3 * (x253 + x305)
    x308 = x191 * x249
    x309 = x307 + x308
    x310 = x191 * x252
    x311 = x3 * (2 * x135 + x246 + x249 + x258)
    x312 = x310 + x311
    x313 = x191 * x256
    x314 = x3 * (4 * x150 + 2 * x201 + 2 * x247 + 2 * x248)
    x315 = x313 + x314
    x316 = x191 * x260
    x317 = x3 * (2 * x160 + 2 * x203 + 2 * x250 + 2 * x251 + x256)
    x318 = x316 + x317
    x319 = x191 * x263
    x320 = 3 * x254 + 3 * x255
    x321 = x3 * (2 * x175 + 2 * x205 + x320)
    x322 = x319 + x321
    x323 = x191 * x266
    x324 = 3 * x257 + 3 * x259
    x325 = x3 * (2 * x180 + 2 * x207 + x263 + x324)
    x326 = x323 + x325
    x327 = x191 * x269 + x3 * (2 * x210 + 2 * x211 + 4 * x261 + 4 * x262)
    x328 = x191 * x272 + x3 * (2 * x214 + 2 * x215 + 4 * x264 + 4 * x265 + x269)
    x329 = x7 * x75
    x330 = x328 * x329
    x331 = x2 * x329
    x332 = numpy.pi * x0 * x7 * x73
    x333 = x2 * x332
    x334 = 2 * x117 * x217 + x217 * x274
    x335 = x170 + x273
    x336 = x217 * x277 + x3 * (2 * x218 + x335)
    x337 = x3 * (x284 + x335)
    x338 = x217 * x280
    x339 = x337 + x338
    x340 = x217 * x283
    x341 = x3 * (2 * x147 + x277 + x280 + x289)
    x342 = x340 + x341
    x343 = x217 * x287
    x344 = x3 * (4 * x166 + 2 * x226 + 2 * x278 + 2 * x279)
    x345 = x343 + x344
    x346 = x217 * x291
    x347 = x3 * (2 * x172 + 2 * x228 + 2 * x281 + 2 * x282 + x287)
    x348 = x346 + x347
    x349 = x217 * x294
    x350 = 3 * x285 + 3 * x286
    x351 = x3 * (2 * x184 + 2 * x230 + x350)
    x352 = x349 + x351
    x353 = x217 * x297
    x354 = 3 * x288 + 3 * x290
    x355 = x3 * (2 * x188 + 2 * x232 + x294 + x354)
    x356 = x353 + x355
    x357 = x217 * x300 + x3 * (2 * x235 + 2 * x236 + 4 * x292 + 4 * x293)
    x358 = x217 * x303 + x3 * (2 * x239 + 2 * x240 + 4 * x295 + 4 * x296 + x300)
    x359 = x332 * x358
    x360 = x191 * x304 + x3 * (x158 + 3 * x242)
    x361 = x191 * x306 + x3 * (3 * x244 + 3 * x245 + x304)
    x362 = x191 * x309 + x3 * (3 * x247 + 3 * x248 + x304)
    x363 = x191 * x312 + x3 * (3 * x250 + 3 * x251 + x306 + x309)
    x364 = x191 * x315 + x3 * (2 * x307 + 2 * x308 + x320)
    x365 = x191 * x318 + x3 * (2 * x310 + 2 * x311 + x315 + x324)
    x366 = x191 * x322 + x3 * (3 * x261 + 3 * x262 + 3 * x313 + 3 * x314)
    x367 = x329 * (x191 * x326 + x3 * (3 * x264 + 3 * x265 + 3 * x316 + 3 * x317 + x322))
    x368 = x329 * x5
    x369 = x329 * (x191 * x327 + x3 * (3 * x267 + 3 * x268 + 4 * x319 + 4 * x321))
    x370 = x329 * x4
    x371 = x191 * x332
    x372 = x217 * x334 + x3 * (x170 + 3 * x273)
    x373 = x217 * x336 + x3 * (3 * x275 + 3 * x276 + x334)
    x374 = x217 * x339 + x3 * (3 * x278 + 3 * x279 + x334)
    x375 = x217 * x342 + x3 * (3 * x281 + 3 * x282 + x336 + x339)
    x376 = x217 * x345 + x3 * (2 * x337 + 2 * x338 + x350)
    x377 = x217 * x348 + x3 * (2 * x340 + 2 * x341 + x345 + x354)
    x378 = x332 * x5
    x379 = x217 * x352 + x3 * (3 * x292 + 3 * x293 + 3 * x343 + 3 * x344)
    x380 = x332 * (x217 * x356 + x3 * (3 * x295 + 3 * x296 + 3 * x346 + 3 * x347 + x352))
    x381 = x332 * (x217 * x357 + x3 * (3 * x298 + 3 * x299 + 4 * x349 + 4 * x351))

    # 675 item(s)
    S = numpy.array(
        [
            x76 * (x2 * x57 + x3 * (3 * x40 + 3 * x46 + 4 * x69 + 4 * x71 + x72)),
            x78 * x86,
            x86 * x88,
            x100 * x89,
            x106 * x113 * x114,
            x115 * x88 * x89,
            x100 * x116,
            x115 * x116 * x78,
            x103 * x113 * x120,
            x114 * x122 * x133,
            x114 * x137 * x142,
            x122 * x142 * x143,
            x133 * x144 * x89,
            x106 * x118 * x142,
            x104 * x120 * x142,
            x103 * x133 * x146,
            x134 * x142 * x146,
            x103 * x142 * x149,
            x114 * x153 * x157,
            x114 * x162 * x165,
            x143 * x153 * x165,
            x118 * x122 * x157,
            x118 * x137 * x165,
            x120 * x122 * x165,
            x104 * x146 * x157,
            x106 * x146 * x165,
            x104 * x149 * x165,
            x103 * x157 * x169,
            x134 * x165 * x169,
            x103 * x165 * x174,
            x114 * x177 * x178,
            x114 * x182 * x183,
            x143 * x177 * x183,
            x118 * x153 * x178,
            x118 * x162 * x183,
            x120 * x153 * x183,
            x122 * x146 * x178,
            x137 * x146 * x183,
            x122 * x149 * x183,
            x104 * x169 * x178,
            x106 * x169 * x183,
            x104 * x174 * x183,
            x103 * x178 * x186,
            x134 * x183 * x186,
            x103 * x183 * x190,
            x191 * x192,
            x114 * x194 * x72,
            x191 * x195 * x88,
            x114 * x197 * x90,
            x114 * x199 * x99,
            x143 * x197 * x99,
            x144 * x191 * x90,
            x118 * x194 * x99,
            x120 * x200 * x99,
            x114 * x123 * x202,
            x114 * x132 * x204,
            x132 * x143 * x202,
            x118 * x123 * x197,
            x118 * x132 * x199,
            x120 * x132 * x197,
            x123 * x146 * x200,
            x132 * x146 * x194,
            x132 * x149 * x200,
            x114 * x154 * x206,
            x114 * x155 * x208,
            x143 * x155 * x206,
            x118 * x154 * x202,
            x118 * x155 * x204,
            x120 * x155 * x202,
            x146 * x154 * x197,
            x146 * x155 * x199,
            x149 * x155 * x197,
            x154 * x169 * x200,
            x155 * x169 * x194,
            x155 * x174 * x200,
            x114 * x156 * x212,
            x114 * x164 * x216,
            x143 * x164 * x212,
            x118 * x156 * x206,
            x118 * x164 * x208,
            x120 * x164 * x206,
            x146 * x156 * x202,
            x146 * x164 * x204,
            x149 * x164 * x202,
            x156 * x169 * x197,
            x164 * x169 * x199,
            x164 * x174 * x197,
            x156 * x186 * x200,
            x164 * x186 * x194,
            x164 * x190 * x200,
            x192 * x217,
            x195 * x217 * x78,
            x103 * x219 * x72,
            x220 * x89 * x90,
            x106 * x221 * x99,
            x104 * x219 * x99,
            x103 * x223 * x90,
            x134 * x223 * x99,
            x103 * x225 * x99,
            x122 * x123 * x221,
            x132 * x137 * x221,
            x122 * x132 * x219,
            x104 * x123 * x223,
            x106 * x132 * x223,
            x104 * x132 * x225,
            x103 * x123 * x227,
            x132 * x134 * x227,
            x103 * x132 * x229,
            x153 * x154 * x221,
            x155 * x162 * x221,
            x153 * x155 * x219,
            x122 * x154 * x223,
            x137 * x155 * x223,
            x122 * x155 * x225,
            x104 * x154 * x227,
            x106 * x155 * x227,
            x104 * x155 * x229,
            x103 * x154 * x231,
            x134 * x155 * x231,
            x103 * x155 * x233,
            x156 * x177 * x221,
            x164 * x182 * x221,
            x164 * x177 * x219,
            x153 * x156 * x223,
            x162 * x164 * x223,
            x153 * x164 * x225,
            x122 * x156 * x227,
            x137 * x164 * x227,
            x122 * x164 * x229,
            x104 * x156 * x231,
            x106 * x164 * x231,
            x104 * x164 * x233,
            x103 * x156 * x237,
            x134 * x164 * x237,
            x103 * x164 * x241,
            x114 * x243 * x47,
            x114 * x246 * x56,
            x143 * x243 * x56,
            x114 * x249 * x70,
            x114 * x252 * x61,
            x143 * x249 * x61,
            x118 * x243 * x70,
            x118 * x246 * x61,
            x120 * x243 * x61,
            x114 * x256 * x97,
            x114 * x260 * x95,
            x143 * x256 * x95,
            x118 * x249 * x97,
            x118 * x252 * x95,
            x120 * x249 * x95,
            x146 * x243 * x97,
            x146 * x246 * x95,
            x149 * x243 * x95,
            x114 * x130 * x263,
            x114 * x124 * x266,
            x124 * x143 * x263,
            x118 * x130 * x256,
            x118 * x124 * x260,
            x120 * x124 * x256,
            x130 * x146 * x249,
            x124 * x146 * x252,
            x124 * x149 * x249,
            x130 * x169 * x243,
            x124 * x169 * x246,
            x124 * x174 * x243,
            x114 * x128 * x269,
            x114 * x163 * x272,
            x143 * x163 * x269,
            x118 * x128 * x263,
            x118 * x163 * x266,
            x120 * x163 * x263,
            x128 * x146 * x256,
            x146 * x163 * x260,
            x149 * x163 * x256,
            x128 * x169 * x249,
            x163 * x169 * x252,
            x163 * x174 * x249,
            x128 * x186 * x243,
            x163 * x186 * x246,
            x163 * x190 * x243,
            x191 * x220 * x47,
            x194 * x221 * x56,
            x200 * x219 * x56,
            x197 * x221 * x70,
            x199 * x221 * x61,
            x197 * x219 * x61,
            x200 * x223 * x70,
            x194 * x223 * x61,
            x200 * x225 * x61,
            x202 * x221 * x97,
            x204 * x221 * x95,
            x202 * x219 * x95,
            x197 * x223 * x97,
            x199 * x223 * x95,
            x197 * x225 * x95,
            x200 * x227 * x97,
            x194 * x227 * x95,
            x200 * x229 * x95,
            x130 * x206 * x221,
            x124 * x208 * x221,
            x124 * x206 * x219,
            x130 * x202 * x223,
            x124 * x204 * x223,
            x124 * x202 * x225,
            x130 * x197 * x227,
            x124 * x199 * x227,
            x124 * x197 * x229,
            x130 * x200 * x231,
            x124 * x194 * x231,
            x124 * x200 * x233,
            x128 * x212 * x221,
            x163 * x216 * x221,
            x163 * x212 * x219,
            x128 * x206 * x223,
            x163 * x208 * x223,
            x163 * x206 * x225,
            x128 * x202 * x227,
            x163 * x204 * x227,
            x163 * x202 * x229,
            x128 * x197 * x231,
            x163 * x199 * x231,
            x163 * x197 * x233,
            x128 * x200 * x237,
            x163 * x194 * x237,
            x163 * x200 * x241,
            x103 * x274 * x47,
            x134 * x274 * x56,
            x103 * x277 * x56,
            x104 * x274 * x70,
            x106 * x274 * x61,
            x104 * x277 * x61,
            x103 * x280 * x70,
            x134 * x280 * x61,
            x103 * x283 * x61,
            x122 * x274 * x97,
            x137 * x274 * x95,
            x122 * x277 * x95,
            x104 * x280 * x97,
            x106 * x280 * x95,
            x104 * x283 * x95,
            x103 * x287 * x97,
            x134 * x287 * x95,
            x103 * x291 * x95,
            x130 * x153 * x274,
            x124 * x162 * x274,
            x124 * x153 * x277,
            x122 * x130 * x280,
            x124 * x137 * x280,
            x122 * x124 * x283,
            x104 * x130 * x287,
            x106 * x124 * x287,
            x104 * x124 * x291,
            x103 * x130 * x294,
            x124 * x134 * x294,
            x103 * x124 * x297,
            x128 * x177 * x274,
            x163 * x182 * x274,
            x163 * x177 * x277,
            x128 * x153 * x280,
            x162 * x163 * x280,
            x153 * x163 * x283,
            x122 * x128 * x287,
            x137 * x163 * x287,
            x122 * x163 * x291,
            x104 * x128 * x294,
            x106 * x163 * x294,
            x104 * x163 * x297,
            x103 * x128 * x300,
            x134 * x163 * x300,
            x103 * x163 * x303,
            x114 * x304 * x39,
            x114 * x306 * x45,
            x143 * x304 * x45,
            x114 * x309 * x52,
            x114 * x312 * x50,
            x143 * x309 * x50,
            x118 * x304 * x52,
            x118 * x306 * x50,
            x120 * x304 * x50,
            x114 * x315 * x66,
            x114 * x318 * x64,
            x143 * x315 * x64,
            x118 * x309 * x66,
            x118 * x312 * x64,
            x120 * x309 * x64,
            x146 * x304 * x66,
            x146 * x306 * x64,
            x149 * x304 * x64,
            x114 * x322 * x93,
            x110 * x114 * x326,
            x110 * x143 * x322,
            x118 * x315 * x93,
            x110 * x118 * x318,
            x110 * x120 * x315,
            x146 * x309 * x93,
            x110 * x146 * x312,
            x110 * x149 * x309,
            x169 * x304 * x93,
            x110 * x169 * x306,
            x110 * x174 * x304,
            x114 * x126 * x327,
            x2 * x330,
            x327 * x331 * x88,
            x118 * x126 * x322,
            x116 * x326 * x331,
            x108 * x120 * x322,
            x126 * x146 * x315,
            x108 * x146 * x318,
            x108 * x149 * x315,
            x126 * x169 * x309,
            x108 * x169 * x312,
            x108 * x174 * x309,
            x126 * x186 * x304,
            x108 * x186 * x306,
            x108 * x190 * x304,
            x221 * x243 * x39,
            x221 * x246 * x45,
            x219 * x243 * x45,
            x221 * x249 * x52,
            x221 * x252 * x50,
            x219 * x249 * x50,
            x223 * x243 * x52,
            x223 * x246 * x50,
            x225 * x243 * x50,
            x221 * x256 * x66,
            x221 * x260 * x64,
            x219 * x256 * x64,
            x223 * x249 * x66,
            x223 * x252 * x64,
            x225 * x249 * x64,
            x227 * x243 * x66,
            x227 * x246 * x64,
            x229 * x243 * x64,
            x221 * x263 * x93,
            x110 * x221 * x266,
            x110 * x219 * x263,
            x223 * x256 * x93,
            x110 * x223 * x260,
            x110 * x225 * x256,
            x227 * x249 * x93,
            x110 * x227 * x252,
            x110 * x229 * x249,
            x231 * x243 * x93,
            x110 * x231 * x246,
            x110 * x233 * x243,
            x126 * x221 * x269,
            x217 * x272 * x331,
            x108 * x219 * x269,
            x126 * x223 * x263,
            x108 * x223 * x266,
            x108 * x225 * x263,
            x126 * x227 * x256,
            x108 * x227 * x260,
            x108 * x229 * x256,
            x126 * x231 * x249,
            x108 * x231 * x252,
            x108 * x233 * x249,
            x126 * x237 * x243,
            x108 * x237 * x246,
            x108 * x241 * x243,
            x200 * x274 * x39,
            x194 * x274 * x45,
            x200 * x277 * x45,
            x197 * x274 * x52,
            x199 * x274 * x50,
            x197 * x277 * x50,
            x200 * x280 * x52,
            x194 * x280 * x50,
            x200 * x283 * x50,
            x202 * x274 * x66,
            x204 * x274 * x64,
            x202 * x277 * x64,
            x197 * x280 * x66,
            x199 * x280 * x64,
            x197 * x283 * x64,
            x200 * x287 * x66,
            x194 * x287 * x64,
            x200 * x291 * x64,
            x206 * x274 * x93,
            x110 * x208 * x274,
            x110 * x206 * x277,
            x202 * x280 * x93,
            x110 * x204 * x280,
            x110 * x202 * x283,
            x197 * x287 * x93,
            x110 * x199 * x287,
            x110 * x197 * x291,
            x200 * x294 * x93,
            x110 * x194 * x294,
            x110 * x200 * x297,
            x126 * x212 * x274,
            x108 * x216 * x274,
            x108 * x212 * x277,
            x126 * x206 * x280,
            x108 * x208 * x280,
            x108 * x206 * x283,
            x126 * x202 * x287,
            x108 * x204 * x287,
            x108 * x202 * x291,
            x126 * x197 * x294,
            x108 * x199 * x294,
            x108 * x197 * x297,
            x126 * x200 * x300,
            x108 * x194 * x300,
            x191 * x303 * x333,
            x103 * x334 * x39,
            x134 * x334 * x45,
            x103 * x336 * x45,
            x104 * x334 * x52,
            x106 * x334 * x50,
            x104 * x336 * x50,
            x103 * x339 * x52,
            x134 * x339 * x50,
            x103 * x342 * x50,
            x122 * x334 * x66,
            x137 * x334 * x64,
            x122 * x336 * x64,
            x104 * x339 * x66,
            x106 * x339 * x64,
            x104 * x342 * x64,
            x103 * x345 * x66,
            x134 * x345 * x64,
            x103 * x348 * x64,
            x153 * x334 * x93,
            x110 * x162 * x334,
            x110 * x153 * x336,
            x122 * x339 * x93,
            x110 * x137 * x339,
            x110 * x122 * x342,
            x104 * x345 * x93,
            x106 * x110 * x345,
            x104 * x110 * x348,
            x103 * x352 * x93,
            x110 * x134 * x352,
            x103 * x110 * x356,
            x126 * x177 * x334,
            x108 * x182 * x334,
            x108 * x177 * x336,
            x126 * x153 * x339,
            x108 * x162 * x339,
            x108 * x153 * x342,
            x122 * x126 * x345,
            x108 * x137 * x345,
            x108 * x122 * x348,
            x104 * x126 * x352,
            x106 * x108 * x352,
            x333 * x356 * x89,
            x103 * x126 * x357,
            x333 * x357 * x78,
            x2 * x359,
            x114 * x360 * x37,
            x114 * x31 * x361,
            x143 * x31 * x360,
            x114 * x35 * x362,
            x114 * x29 * x363,
            x143 * x29 * x362,
            x118 * x35 * x360,
            x118 * x29 * x361,
            x120 * x29 * x360,
            x114 * x22 * x364,
            x114 * x27 * x365,
            x143 * x27 * x364,
            x118 * x22 * x362,
            x118 * x27 * x363,
            x120 * x27 * x362,
            x146 * x22 * x360,
            x146 * x27 * x361,
            x149 * x27 * x360,
            x114 * x20 * x366,
            x367 * x5,
            x366 * x368 * x88,
            x118 * x20 * x364,
            x116 * x365 * x368,
            x10 * x120 * x364,
            x146 * x20 * x362,
            x10 * x146 * x363,
            x10 * x149 * x362,
            x169 * x20 * x360,
            x10 * x169 * x361,
            x10 * x174 * x360,
            x369 * x4,
            x329
            * (x191 * x328 + x3 * (3 * x270 + 3 * x271 + 4 * x323 + 4 * x325 + x327)),
            x369 * x88,
            x116 * x366 * x370,
            x116 * x367,
            x120 * x366 * x9,
            x146 * x18 * x364,
            x146 * x365 * x9,
            x149 * x364 * x9,
            x169 * x18 * x362,
            x169 * x363 * x9,
            x174 * x362 * x9,
            x18 * x186 * x360,
            x186 * x361 * x9,
            x190 * x360 * x9,
            x221 * x304 * x37,
            x221 * x306 * x31,
            x219 * x304 * x31,
            x221 * x309 * x35,
            x221 * x29 * x312,
            x219 * x29 * x309,
            x223 * x304 * x35,
            x223 * x29 * x306,
            x225 * x29 * x304,
            x22 * x221 * x315,
            x221 * x27 * x318,
            x219 * x27 * x315,
            x22 * x223 * x309,
            x223 * x27 * x312,
            x225 * x27 * x309,
            x22 * x227 * x304,
            x227 * x27 * x306,
            x229 * x27 * x304,
            x20 * x221 * x322,
            x217 * x326 * x368,
            x10 * x219 * x322,
            x20 * x223 * x315,
            x10 * x223 * x318,
            x10 * x225 * x315,
            x20 * x227 * x309,
            x10 * x227 * x312,
            x10 * x229 * x309,
            x20 * x231 * x304,
            x10 * x231 * x306,
            x10 * x233 * x304,
            x217 * x327 * x370,
            x217 * x330,
            x219 * x327 * x9,
            x18 * x223 * x322,
            x223 * x326 * x9,
            x225 * x322 * x9,
            x18 * x227 * x315,
            x227 * x318 * x9,
            x229 * x315 * x9,
            x18 * x231 * x309,
            x231 * x312 * x9,
            x233 * x309 * x9,
            x18 * x237 * x304,
            x237 * x306 * x9,
            x241 * x304 * x9,
            x243 * x274 * x37,
            x246 * x274 * x31,
            x243 * x277 * x31,
            x249 * x274 * x35,
            x252 * x274 * x29,
            x249 * x277 * x29,
            x243 * x280 * x35,
            x246 * x280 * x29,
            x243 * x283 * x29,
            x22 * x256 * x274,
            x260 * x27 * x274,
            x256 * x27 * x277,
            x22 * x249 * x280,
            x252 * x27 * x280,
            x249 * x27 * x283,
            x22 * x243 * x287,
            x246 * x27 * x287,
            x243 * x27 * x291,
            x20 * x263 * x274,
            x10 * x266 * x274,
            x10 * x263 * x277,
            x20 * x256 * x280,
            x10 * x260 * x280,
            x10 * x256 * x283,
            x20 * x249 * x287,
            x10 * x252 * x287,
            x10 * x249 * x291,
            x20 * x243 * x294,
            x10 * x246 * x294,
            x10 * x243 * x297,
            x18 * x269 * x274,
            x272 * x274 * x9,
            x269 * x277 * x9,
            x18 * x263 * x280,
            x266 * x280 * x9,
            x263 * x283 * x9,
            x18 * x256 * x287,
            x260 * x287 * x9,
            x256 * x291 * x9,
            x18 * x249 * x294,
            x252 * x294 * x9,
            x249 * x297 * x9,
            x18 * x243 * x300,
            x246 * x300 * x9,
            x243 * x303 * x9,
            x200 * x334 * x37,
            x194 * x31 * x334,
            x200 * x31 * x336,
            x197 * x334 * x35,
            x199 * x29 * x334,
            x197 * x29 * x336,
            x200 * x339 * x35,
            x194 * x29 * x339,
            x200 * x29 * x342,
            x202 * x22 * x334,
            x204 * x27 * x334,
            x202 * x27 * x336,
            x197 * x22 * x339,
            x199 * x27 * x339,
            x197 * x27 * x342,
            x200 * x22 * x345,
            x194 * x27 * x345,
            x200 * x27 * x348,
            x20 * x206 * x334,
            x10 * x208 * x334,
            x10 * x206 * x336,
            x20 * x202 * x339,
            x10 * x204 * x339,
            x10 * x202 * x342,
            x197 * x20 * x345,
            x10 * x199 * x345,
            x10 * x197 * x348,
            x20 * x200 * x352,
            x10 * x194 * x352,
            x356 * x371 * x5,
            x18 * x212 * x334,
            x216 * x334 * x9,
            x212 * x336 * x9,
            x18 * x206 * x339,
            x208 * x339 * x9,
            x206 * x342 * x9,
            x18 * x202 * x345,
            x204 * x345 * x9,
            x202 * x348 * x9,
            x18 * x197 * x352,
            x199 * x352 * x9,
            x197 * x356 * x9,
            x357 * x371 * x4,
            x194 * x357 * x9,
            x191 * x359,
            x103 * x37 * x372,
            x134 * x31 * x372,
            x103 * x31 * x373,
            x104 * x35 * x372,
            x106 * x29 * x372,
            x104 * x29 * x373,
            x103 * x35 * x374,
            x134 * x29 * x374,
            x103 * x29 * x375,
            x122 * x22 * x372,
            x137 * x27 * x372,
            x122 * x27 * x373,
            x104 * x22 * x374,
            x106 * x27 * x374,
            x104 * x27 * x375,
            x103 * x22 * x376,
            x134 * x27 * x376,
            x103 * x27 * x377,
            x153 * x20 * x372,
            x10 * x162 * x372,
            x10 * x153 * x373,
            x122 * x20 * x374,
            x10 * x137 * x374,
            x10 * x122 * x375,
            x104 * x20 * x376,
            x10 * x106 * x376,
            x377 * x378 * x89,
            x103 * x20 * x379,
            x378 * x379 * x78,
            x380 * x5,
            x177 * x18 * x372,
            x182 * x372 * x9,
            x177 * x373 * x9,
            x153 * x18 * x374,
            x162 * x374 * x9,
            x153 * x375 * x9,
            x122 * x18 * x376,
            x137 * x376 * x9,
            x122 * x377 * x9,
            x332 * x379 * x4 * x89,
            x106 * x379 * x9,
            x380 * x89,
            x381 * x4,
            x381 * x78,
            x332
            * (x217 * x358 + x3 * (3 * x301 + 3 * x302 + 4 * x353 + 4 * x355 + x357)),
        ]
    )
    return S
