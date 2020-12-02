from collections import namedtuple


Molecule = namedtuple("Molecule", "name fn charge mult density")

MOLECULES = (
    # Densities in g/ml at 20°C
    Molecule(
        "water",
        "water.xyz",
        0,
        1,
        0.9982,
    ),
    Molecule(
        "acetone",
        "acetone.xyz",
        0,
        1,
        0.790,
    ),
    Molecule(
        "acetonitrile",
        "acetonitrile.xyz",
        0,
        1,
        0.7822,
    ),
    Molecule(
        "dichloromethane",
        "dichloromethane.xyz",
        0,
        1,
        1.322,
    ),
    Molecule(
        "chloroform",
        "chloroform.xyz",
        0,
        1,
        1.489,
    ),
)
