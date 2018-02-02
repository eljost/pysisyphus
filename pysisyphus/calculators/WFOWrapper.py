#!/usr/bin/env python3

import itertools

import numpy as np
from pyscf import gto

from pysisyphus.config import Config


class WFOWrapper:
    def __init__(self, occ_mos, virt_mos):
        #self.base_cmd = Config["wfoverlap"]["cmd"]
        #self.atoms = atoms
        self.occ_mos = occ_mos
        self.virt_mos = virt_mos
        self.base_det_str = "d"*self.occ_mos + "e"*self.virt_mos

        # Molecule coordinates
        self.coords_list = list()
        # MO coefficients
        self.mos_list = list()
        self.ci_coeffs_list = list()
        self.mo_inds_list = list()
        self.from_set_list = list()
        self.to_set_list = list()

    def build_mole(self, coords):
        coords3d = coords.reshape(-1, 3)
        return [[atom, c] for atom, c in zip(self.atoms, coords3d)]

    def overlap(self, a1, a2, basis="sto-3g"):
        def prepare(atom):
            mol = gto.Mole()
            mol.atom = atom
            mol.basis = basis
            # Charge or spin aren't needed for overlap integrals
            mol.build()
            return mol
        mol1 = prepare(a1)
        mol2 = prepare(a2)
        ao_ovlp = gto.mole.intor_cross("int1e_ovlp_sph", mol1, mol2)
        return ao_ovlp


    def make_det_string(self, inds):
        """Return spin adapted strings."""
        from_mo, to_mo = inds
        # Until now the first virtual MO (to_mo) has index 0. To subsitute
        # the base_str at the correct index we have to increase all to_mo
        # indices by the number off occupied MO.
        to_mo += self.occ_mos
        # Make string for excitation of an alpha electron
        ab = list(self.base_det_str)
        ab[from_mo] = "b"
        ab[to_mo] = "a"
        ab_str = "".join(ab)
        # Make string for excitation of an beta electron
        ba = list(self.base_det_str)
        ba[from_mo] = "a"
        ba[to_mo] = "b"
        ba_str = "".join(ba)
        return ab_str, ba_str

    def generate_all_dets(self, occ_set1, virt_set1, occ_set2, virt_set2):
        """Generate all possible single excitation determinant strings
        from union(occ_mos) to union(virt_mos)."""
        # Unite the respective sets of both calculations
        occ_set = occ_set1 | occ_set2
        virt_set = virt_set1 | virt_set2
        # Genrate all possible excitations (combinations) from the occupied
        # MO set to (and) the virtual MO set.
        all_inds = [(om, vm) for om, vm
                    in itertools.product(occ_set, virt_set)]
        det_strings = [self.make_det_string(inds) for inds in all_inds]
        return all_inds, det_strings

    def set_from_nested_list(self, nested):
        return set([i for i in itertools.chain(*nested)])

    def store_iteration(self, atoms, coords, mos, ci_coeffs, mo_inds):
        from_mos, to_mos = zip(*mo_inds)
        from_set = self.set_from_nested_list(from_mos)
        to_set = self.set_from_nested_list(to_mos)

        self.atoms = atoms
        self.mos_list = mos
        self.coords_list.append(coords)
        self.ci_coeffs_list.append(ci_coeffs)
        self.mo_inds_list.append(mo_inds)
        self.from_set_list.append(from_set)
        self.to_set_list.append(to_set)

    def get_iteration(self, ind):
        return (self.mos_list[ind], self.coords_list[ind],
                self.ci_coeffs_list[ind], self.mo_inds_list[ind],
                self.from_set_list[ind], self.to_set_list[ind])

    def track(self):
        if len(self.ci_coeffs_list) < 2:
            print("can't track yet!")
            return None
        mos1, coords1, cic1, moi1, fs1, ts1 = self.get_iteration(-2)
        mos2, coords2, cic2, moi2, fs2, ts2 = self.get_iteration(-1)
        a1 = self.build_mole(coords1)
        a2 = self.build_mole(coords2)
        ao_ovlp = self.overlap(a1, a2, basis="sto-3g")
        all_inds, det_strings = self.generate_all_dets(fs1, ts1, fs2, ts2)
        import pdb; pdb.set_trace()
        return None



if __name__ == "__main__":
    oc1 = np.arange(5)
    oc2 = oc1.copy()
    vs1 = np.arange(2)
    vs2 = vs1.copy()
    occ_mos = len(oc1)
    virt_mos = len(vs1)
    wfow = WFOWrapper(occ_mos, virt_mos)
    all_inds, det_strings = wfow.generate_all_dets(oc1, vs1, oc2, vs2)
    import pdb; pdb.set_trace()
