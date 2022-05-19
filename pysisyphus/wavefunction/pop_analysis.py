import numpy as np
import scipy as sp


def mulliken(P, S, atom_num, nuc_charges, ao_centers):
    def mulliken_atom_pops(P, S):
        mo_populations = np.einsum("ij,ji->i", P, S)
        atom_populations = np.zeros(atom_num)
        for i, center in enumerate(ao_centers):
            atom_populations[center] += mo_populations[i]
        return atom_populations

    if P.ndim == 3:
        atom_populations_a = mulliken_atom_pops(P[0], S)
        atom_populations_b = mulliken_atom_pops(P[1], S)
    else:
        atom_populations_a = mulliken_atom_pops(P, S) / 2
        atom_populations_b = atom_populations_a

    charges = nuc_charges - atom_populations_a - atom_populations_b
    return charges


def mulliken_from_wf(wf):
    return mulliken(
        P=wf.P,
        S=wf.S,
        atom_num=wf.atom_num,
        nuc_charges=wf.nuc_charges,
        ao_centers=wf.ao_centers,
    )
