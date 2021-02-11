from collections import namedtuple


Molecule = namedtuple("Molecule", "name fn charge mult density resname")

MOLECULES = (
    # Densities in g/ml at 20°C
    Molecule(
        "water",
        "water.xyz",
        0,
        1,
        0.9982,
        "H2O",
    ),
    Molecule(
        "acetone",
        "acetone.xyz",
        0,
        1,
        0.790,
        "ACO",
    ),
    Molecule(
        "acetonitrile",
        "acetonitrile.xyz",
        0,
        1,
        0.7822,
        "ACN",
    ),
    Molecule(
        "dichloromethane",
        "dichloromethane.xyz",
        0,
        1,
        1.322,
        "DCM",
    ),
    Molecule(
        "chloroform",
        "chloroform.xyz",
        0,
        1,
        1.489,
        "CLF",
    ),
    Molecule(
        "tetrahydrofuran",
        "tetrahydrofuran.xyz",
        0,
        1,
        0.8876,
        "THF",
    ),
    Molecule(
        "14dioxane",
        "14dioxane.xyz",
        0,
        1,
        1.03,
        "DOX",
    ),
)
