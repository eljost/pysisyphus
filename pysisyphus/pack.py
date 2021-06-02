import argparse
from math import pi as PI, ceil
from pathlib import Path
import sys

from pysisyphus.constants import AMU2KG
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import get_input, highlight_text
from pysisyphus.io.pdb import geom_to_pdb_str
from pysisyphus.wrapper.packmol import make_input, call_packmol

from pysisyphus.db import LEVELS, MOLECULES
from pysisyphus.db.helpers import get_path as db_get_path


AMU2G = AMU2KG * 1e3
CM2ANG = 1e-8


def parse_args(args):
    parser = argparse.ArgumentParser()

    # Solvent
    solvent_group = parser.add_mutually_exclusive_group(required=True)
    solvent_group.add_argument("--solv", help="Filename of solvent geometry.")
    solvent_group.add_argument(
        "--db", action="store_true", help="Choose from internal database."
    )

    parser.add_argument(
        "--solv_num", type=int, help="Number of solvent molecules to pack."
    )
    parser.add_argument("--solv_dens", type=float, help="Solvent density in g/cm³.")

    parser.add_argument(
        "--output", default="output.pdb", help="Filename of packed molecules."
    )

    # Solute
    parser.add_argument("--solute", default=None, help="Filename of solute geometry.")
    parser.add_argument(
        "--solute_num", type=int, default=1, help="Number of solute molecules to pack."
    )

    return parser.parse_args(args)


def as_pdb(fn):
    if not str(fn).endswith(".pdb"):
        geom = geom_loader(fn)
        pdb_str = geom_to_pdb_str(geom)
        cwd = Path(".")
        pdb_fn = cwd / Path(fn).with_suffix(".pdb").name
        with open(pdb_fn, "w") as handle:
            handle.write(pdb_str)
        print(f"Converted '{fn}' to PDB format ('{pdb_fn}')")
    return pdb_fn


def sphere_radius_from_volume(volume):
    radius = (3 / 4 * volume / PI) ** (1 / 3)
    return radius


def volume_for_density(molecule_num, mol_mass, density):
    # Convert density from g/cm³ to amu/Å³
    density_au = density / AMU2G * CM2ANG ** 3
    # The molar mass in g/mol is numerically equal to the value in AMU (dalton)
    # so we can use it as it is.
    total_mass = mol_mass * molecule_num
    # Volume in Å
    volume = total_mass / density_au

    return volume


def print_info(title, geom):
    print(title)
    print(f"\t{geom}")
    print(f"\tMolar mass: {geom.total_mass:.2f} g mol⁻¹")
    print()


def run():
    args = parse_args(sys.argv[1:])

    solute_fn = args.solute
    if solute_fn:
        solute = geom_loader(solute_fn)
        solute_num = args.solute_num
        solute_mass = solute.total_mass
        print_info("Solute", solute)
    else:
        solute = None
        solute_mass = 0.0

    solv_fn = args.solv
    if solv_fn:
        solv_dens = args.solv_dens
    # Load from internal db
    else:
        print(highlight_text("Interactive solvent selection"))
        level = get_input(LEVELS, "Level of theory", lbl_func=lambda lvl: lvl[0])
        print()
        molecule = get_input(MOLECULES, "Molecule", lbl_func=lambda mol: mol.name)
        print()

        solv_fn = db_get_path(molecule.name, level[0])
        solv_dens = molecule.density

    solv = geom_loader(solv_fn)
    solv_num = args.solv_num
    solv_mass = solv.total_mass
    print_info("Solvent", solv)

    solute_solv_mass = solute_mass + solv_num * solv_mass
    print(f"Total mass of solute(s) and solvent(s): {solute_solv_mass:.2f} amu")
    print()

    # Solvent volume
    solv_vol = volume_for_density(solv_num, solv_mass, solv_dens)
    print(f"Solvent volume: {solv_vol:>10.2f} Å³")

    # Solute volume; Use the solvent density for this calculation
    if solute:
        solute_vol = volume_for_density(solute_num, solute_mass, solv_dens)
        print(f" Solute volume: {solute_vol:>10.2f} Å³")
    else:
        solute_vol = 0.0

    total_vol = solv_vol + solute_vol
    print(f"  Total volume: {total_vol:>10.2f} Å³")
    print()

    radius = sphere_radius_from_volume(total_vol)
    print(f"     Sphere radius: {radius:>8.2f} Å")
    cradius = ceil(radius)
    print(f"Using ceil(radius): {cradius:>8.2f} Å")
    print()

    # Create solute/solvent PDBs if needed

    inp_kwargs = {
        "output_fn": args.output,
        "solvent_fn": as_pdb(solv_fn),
        "solvent_num": solv_num,
        "sphere_radius": cradius,
    }
    if solute:
        inp_kwargs.update({"solute_fn": as_pdb(solute_fn), "solute_num": solute_num})

    inp = make_input(**inp_kwargs)
    inp_fn = "packmol.inp"
    with open(inp_fn, "w") as handle:
        handle.write(inp)
    print(f"Wrote packmol input to '{inp_fn}'")

    proc = call_packmol(inp)

    log_fn = "packmol.log"
    with open(log_fn, "w") as handle:
        handle.write(proc.stdout)
    print(f"Wrote packmol ouput to '{log_fn}'")
    print()

    return_ = proc.returncode
    if return_ != 0:
        print(proc.stdout)
    else:
        print("packmol returned successfully!")

    return return_
