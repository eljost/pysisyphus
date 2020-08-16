import numpy as np


def lincs_closure(geom, constraints, order=4):
    """
    Drop conn_nums and keep conns_ in a list of list, instead of an array.
    We could do everything with
        for i, conn_constr in enumerate(conns)
            pass
    instead of the pseudo-code way as given in the paper.
    """
    # Number of constraints
    K = len(constraints)
    # Inverse atom masses
    inv_masses = 1 / geom.masses

    # List of constraint sets
    constraint_sets = [set(c) for c in constraints]
    constraints = np.array(constraints, dtype=int)
    # Indices of the atoms that make up the bond constraints
    atom_1, atom_2 = constraints.T

    # Original bond lengths
    lengths = np.linalg.norm(geom.coords3d[atom_1] - geom.coords3d[atom_2], axis=1)

    # Determine number of constraints that are connected to every constraint
    conn_nums = np.zeros(K, dtype=int)
    conns_ = [list() for _ in range(K)]
    for i, con_1 in enumerate(constraint_sets):
        for j, con_2 in enumerate(constraint_sets):
            if i == j:
                continue

            if con_1 & con_2:
                conn_nums[i] += 1
                conns_[i].append(j)
    # Maximum number of connected constraints
    c_max = conn_nums.max()
    S_diag = 1 / np.sqrt(inv_masses[atom_1] + inv_masses[atom_2])

    # Array containing the indices of the connected constraints for every constraint
    conns = np.zeros((K, c_max), dtype=int)
    for i, c in enumerate(conns_):
        conns[i, :len(c)] = c

    coefs = np.zeros((K, c_max))
    signs = {
        True: -1,
        False: 1,
    }
    for i, con_1 in enumerate(constraint_sets):
        for j in range(c_max):
            conns_ij = conns[i, j]
            sign = signs[(atom_1[i] == atom_1[conns_ij]) or (atom_2[i] == atom_2[conns_ij])]
            # Connected atom
            c = tuple(con_1 & constraint_sets[conns_ij])
            assert len(c) == 1
            c = int(c[0])
            coefs[i, j] = sign * inv_masses[c] * S_diag[i] * S_diag[conns_ij]

    def solve(rhs, sol, new_coords3d, A, B):
        w = 1
        for _ in range(order):
            for i, _ in enumerate(constraints):
                rhs[w, i] = 0.
                for j in range(conn_nums[i]):
                    # rhs[w, i] += A[i, j] * rhs[2-w, conns_[i][j]]
                    rhs[w, i] += A[i, j] * rhs[2-w, conns[i, j]]
                sol[i] += rhs[w, i]
            w = 2 - w

        for i, (a1, a2) in enumerate(constraints):
            new_coords3d[a1] -= inv_masses[a1] * S_diag[i] * sol[i] * B[i].sum()
            new_coords3d[a2] += inv_masses[a2] * S_diag[i] * sol[i] * B[i].sum()

    def lincs(prev_coords3d, new_coords3d):
        # Calculate constraint directions
        B = prev_coords3d[atom_1] - prev_coords3d[atom_2]
        B /= np.linalg.norm(B, axis=1)[:, None]
        print("B", B)

        A = np.zeros((K, c_max))
        rhs = np.zeros((2, K))
        sol = np.zeros(K)

        for i, (a1, a2) in enumerate(constraints):
            print(i, a1, a2)
            for j, k in enumerate(conns[i]):
                A[i, j] = coefs[i, j] * (B[i] * B[k]).sum()
            
            rhs[0, i] = S_diag[i] * ((B[i] * (new_coords3d[a1] - new_coords3d[a2])).sum() - lengths[i])
            sol[i] = rhs[0, i]
        solve(rhs, sol, new_coords3d, A, B)

        # Correction for rotational lenghthening
        for i, (a1, a2) in enumerate(constraints):
            length_i = lengths[i]
            p = (2 * length_i**2 - (-(new_coords3d[a1] - new_coords3d[a2])**2).sum())**0.5
            rhs[0, i] = S_diag[i] * (length_i - p)
            sol[i] = rhs[0, i]

        solve(rhs, sol, new_coords3d, A, B)
        # import pdb; pdb.set_trace()
        return new_coords3d
    return lincs
