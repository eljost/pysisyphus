#!/usr/bin/env python

import argparse
import sys

import h5py
import physical_validation
import physical_validation as pv

from pysisyphus.constants import AU2KJPERMOL, BOHR2ANG

BOHR2NM = BOHR2ANG / 10


def pysis_to_pv(h5_fn="md.h5", h5_group="run", remove_com_v=True):

    with h5py.File(h5_fn, "r") as handle:
        group = handle[h5_group]

        dt = group.attrs["dt"]
        T = group.attrs["T_target"]
        masses = group.attrs["masses"]
        atoms = group.attrs["atoms"]
        tot_ene = group["energy_tot"][:] * AU2KJPERMOL
        kin_ene = group["energy_kin"][:] * AU2KJPERMOL
        pot_ene = group["energy_pot"][:] * AU2KJPERMOL
        # Trajectory data
        position = group["cart_coords"][:] * BOHR2NM
        velocity = group["velocity"][:] * BOHR2NM

    natoms = len(atoms)
    ndof_reduction_tra = 3 if remove_com_v else 0

    system = pv.data.SystemData(
        natoms=natoms,
        nconstraints=0,
        ndof_reduction_tra=ndof_reduction_tra,
        ndof_reduction_rot=0,
        mass=masses,
    )

    units = pv.data.UnitData(
        kb=8.314462435405199e-3,  # R in kJ/(mol K)
        energy_str="kJ/mol",
        energy_conversion=1.0,
        length_str="nm",
        length_conversion=1.0,
        # volume_str='nm^3',
        volume_conversion=1.0,
        temperature_str="K",
        temperature_conversion=1.0,
        # pressure_str='bar',
        pressure_conversion=1.0,
        time_str="fs",
        time_conversion=1.0,
    )

    ensemble = pv.data.EnsembleData(
        ensemble="NVT",
        natoms=natoms,
        # volume=3.01125**3,
        temperature=T,
    )

    observables = pv.data.ObservableData(
        kinetic_energy=kin_ene,
        potential_energy=pot_ene,
        total_energy=tot_ene,
        # TODO: add conserved quantitiy
        # constant_of_motion=,
    )

    trajectory = pv.data.TrajectoryData(
        position=position,
        velocity=velocity,
    )

    res = pv.data.SimulationData(
        units=units,
        ensemble=ensemble,
        system=system,
        observables=observables,
        trajectory=trajectory,
        dt=dt,
    )

    return res


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("--groups", type=str, nargs="*", default=("run",))

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    groups = args.groups
    results = [pysis_to_pv(h5_group=grp) for grp in groups]

    for grp, res in zip(groups, results):
        print(f"### Validating kinetic energy for {grp}")
        _ = pv.kinetic_energy.distribution(
            res, verbosity=2, strict=True, filename=f"{grp}_kin_test"
        )

    # Needs conserved quantity
    # if len(groups) > 1:
    # print("### Validating integrator energy")
    # i = pv.integrator.convergence(results, verbose=True, filename=f"{grp}_integrator")
    # print(i)

    # Singular Matrix
    # if res2:
    # print("### Validating ensemble")
    # quantiles = pv.ensemble.check(
    # # res, quiet=False, screen=False, filename=f"{grp}_ensemble_test"
    # res, res2, screen=False, filename=f"{grp}_ensemble_test"
    # )
    # print(quantiles)


if __name__ == "__main__":
    run()
