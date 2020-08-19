import numpy as np

from pysisyphus.constants import FORCE2ACC
from pysisyphus.dynamics.helpers import energy_forces_getter_closure, remove_com_velocity


def rattle_closure(geom, constraints, dt, tol=1e-3, max_cycles=25,
                   energy_forces_getter=None, remove_com_v=True):
    # Inverse atom masses
    masses = geom.masses
    inv_masses = 1 / masses
    inv_masses_rep = np.repeat(inv_masses, 3)

    dt2 = 0.5 * dt
    dt2_sq = dt2 * dt
    tol_sq = tol**2

    constraints = np.array(constraints, dtype=int)
    inv_masses_sums = inv_masses[constraints[:, 0]] + inv_masses[constraints[:, 1]]
    coords3d = geom.coords3d

    def get_bond_vecs(coords3d):
        return coords3d[constraints[:, 0]] - coords3d[constraints[:, 1]]

    # Original bond lengths
    lengths = np.linalg.norm(get_bond_vecs(coords3d), axis=1)
    lengths_sq = lengths**2

    if energy_forces_getter is None:
        energy_forces_getter = energy_forces_getter_closure(geom)
        # def forces_getter(coords):
            # return geom.get_energy_and_forces_at(coords)["forces"]

    def rattle(coords, velocities, forces):
        acceleration = forces * inv_masses_rep * FORCE2ACC
        # Bond vectors at time t
        bond_vecs = get_bond_vecs(coords.reshape(-1, 3))

        # Unconstrained position update
        coords3d_updated = (coords + dt*velocities + dt2_sq * acceleration).reshape(-1, 3)
        # Unconstrained velocity update
        velocities3d_updated = (velocities + dt2 * acceleration).reshape(-1, 3)
        if remove_com_v:
            velocities3d_updated = remove_com_velocity(velocities3d_updated, masses)

        # First part of RATTLE algorithm.
        #   Yields updated positions for t+dt and half-updated velocities for t+dt/2
        for i in range(max_cycles):
            corrected = False
            for j, (atom_1, atom_2) in enumerate(constraints):
                bond_vec_updated = coords3d_updated[atom_1] - coords3d_updated[atom_2]
                bond_length_sq = bond_vec_updated.dot(bond_vec_updated)
                ref_length_sq = lengths_sq[j]
                diff_sq = ref_length_sq - bond_length_sq

                # We have to correct coordinates and velocities if the deviations is
                # above the tolerance.
                if abs(diff_sq) <= (ref_length_sq * tol_sq):
                    continue

                corrected = True
                # Otherwise try to satisfy the constraint by calculating the
                # approximate lagrange multiplier 'g'.
                dot = bond_vec_updated.dot(bond_vecs[j])
                # TODO: test for constraint failure
                g = diff_sq / (2 * inv_masses_sums[j] * dot)
                
                # Update positions to satify constraint
                atom_1_factor = g * bond_vecs[j] * inv_masses[atom_1]
                atom_2_factor = g * bond_vecs[j] * inv_masses[atom_2]
                coords3d_updated[atom_1] += atom_1_factor
                coords3d_updated[atom_2] -= atom_2_factor
                # Update velocities to satify constraint
                velocities3d_updated[atom_1] += atom_1_factor / dt
                velocities3d_updated[atom_2] -= atom_2_factor / dt

            # Stop macro-iterations when no correction was done
            if not corrected:
                # print(f"RATTLE_1 finished after {i+1} macro cycle(s)!")
                break
        else:
            raise Exception("First part of RATTLE did not converge!")

        # Calculate energy (E_pot) and forces at updated coordinates
        energy_new, forces_new = energy_forces_getter(coords3d_updated.flatten())
        velocities3d_updated += dt2 * (forces_new * inv_masses_rep * FORCE2ACC).reshape(-1, 3)
        if remove_com_v:
            velocities3d_updated = remove_com_velocity(velocities3d_updated, masses)
        # The coordinates aren't updated anymore in the second part, so we can
        # calculate the bond vectors here once and store them.
        bond_vecs = get_bond_vecs(coords3d_updated)

        # Second part of RATTLE algorithm
        #   Update velocities to full time step t+dt
        for i in range(max_cycles):
            corrected = False
            for j, (atom_1, atom_2) in enumerate(constraints):
                velocity_diff = (velocities3d_updated[atom_1]
                                 - velocities3d_updated[atom_2])
                dot = bond_vecs[j].dot(velocity_diff)
                # Approximate Lagrange multiplier
                k = dot / (inv_masses_sums[j] * lengths_sq[j])

                if abs(k) <= tol_sq:
                    continue

                corrected = True
                # Update velocities to satify constraint
                velocities3d_updated[atom_1] -= k * bond_vecs[j] * inv_masses[atom_1]
                velocities3d_updated[atom_2] += k * bond_vecs[j] * inv_masses[atom_2]

            # Stop macro-iterations when no correction was done
            if not corrected:
                # print(f"RATTLE_2 finished after {i+1} macro cycle(s)!")
                break
        else:
            raise Exception("Second part of RATTLE did not converge!")

        return coords3d_updated.flatten(), velocities3d_updated.flatten(), energy_new, forces_new

    return rattle
