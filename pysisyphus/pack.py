import argparse
from math import pi as PI, ceil
from pathlib import Path
import sys

from pysisyphus.constants import AMU2KG
from pysisyphus.helpers import geom_loader
from pysisyphus.io.pdb import geom_to_pdb_str
from pysisyphus.wrapper.packmol import make_input, call_packmol


AMU2G = AMU2KG * 1e3
CM2ANG = 1e-8


def parse_args(args):
    parser = argparse.ArgumentParser()

    # Solute
    parser.add_argument("solute",
        help="Filename of solute geometry."
    )
    parser.add_argument("--solute_num", type=int, default=1,
        help="Number of solute molecules to pack."
    )

    # Solvent
    parser.add_argument("--solv",
        help="Filename of solvent geometry."
    )
    parser.add_argument("--solv_num", type=int,
        help="Number of solvent molecules to pack."
    )
    parser.add_argument("--solv_dens", type=float,
        help="Solvent density in g/cm³."
    )

    parser.add_argument("--output", default="output.pdb",
        help="Filename of packed molecules."
    )

    return parser.parse_args(args)


def as_pdb(fn):
    if not fn.endswith(".pdb"):
        geom = geom_loader(fn)
        pdb_str = geom_to_pdb_str(geom)
        pdb_fn = str(Path(fn).with_suffix(".pdb"))
        with open(pdb_fn, "w") as handle:
            handle.write(pdb_str)
        print(f"Converted '{fn}' to PDB format ('{pdb_fn}')")
        fn = pdb_fn
    return fn


def sphere_radius_from_volume(volume):
    radius = (3/4 * volume / PI)**(1/3)
    return radius


def volume_for_density(molecule_num, mol_mass, density):
    # Convert density from g/cm³ to amu/Å³
    density_au = density / AMU2G * CM2ANG**3
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
    solute = geom_loader(solute_fn)
    solute_num = args.solute_num
    solute_mass = solute.total_mass

    solv_fn = args.solv
    solv = geom_loader(solv_fn)
    solv_dens = args.solv_dens
    solv_num = args.solv_num
    solv_mass = solv.total_mass

    print_info("Solute", solute)
    print_info("Solvent", solv)

    solute_solv_mass = solute_mass + solv_mass
    print(f"Total mass of solute(s) and solvent(s): {solute_solv_mass:.2f} amu")
    print()

    # Solvent volume
    solv_vol = volume_for_density(solv_num, solv_mass, solv_dens)
    print(f"Solvent volume: {solv_vol:>10.2f} Å³")

    # Solute volume; Use the solvent density for this calculation
    solute_vol = volume_for_density(solute_num, solute_mass, solv_dens)
    print(f" Solute volume: {solute_vol:>10.2f} Å³")

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
        "solute_fn": as_pdb(solute_fn),
        "solute_num": solute_num,
        "solvent_fn": as_pdb(solv_fn),
        "solvent_num": solv_num,
        "sphere_radius": cradius,
    }

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
