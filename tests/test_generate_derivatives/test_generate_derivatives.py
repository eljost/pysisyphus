from pysisyphus.intcoords.generate_derivatives import generate_wilson


def test_generate_derivatives():
    # generate = ("bond", "bend", "linear_bend", "dihedral")
    generate = ("bond", "bend", "linear_bend")
    out_fn = False

    derivs = generate_wilson(generate, out_fn)

    # bond, bend, linear_bend = 
    bond, bend, linear_bend = derivs

    assert len(bond.d1) == 6
    assert len(bond.d2) == 36

    assert len(bend.d1) == 9
    assert len(bend.d2) == 81

    assert len(linear_bend.d1) == 9
    assert len(linear_bend.d2) == 81
