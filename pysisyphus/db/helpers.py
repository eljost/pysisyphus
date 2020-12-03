import itertools as it
from pysisyphus.db import MOLECULES, LEVELS, LEVEL_DIR


def full_name(molecule_name, level_name):
    return f"{molecule_name}_{level_name}"


def get_path(molecule_name, level_name):
    return LEVEL_DIR / level_name / f"{molecule_name}.xyz"


def get_molecules_levels():
    # return list(
    molecules_levels = list()
    for level, molecule in it.product(LEVELS, MOLECULES):
        level_name, *_ = level
        molecule_name, *_ = molecule
        molecules_levels.append(full_name(molecule_name, level_name))
    return molecules_levels


if __name__ == "__main__":
    from pprint import pprint
    pprint(get_molecules_levels())
    print(get_path("water", "gfn0_xtb"))
