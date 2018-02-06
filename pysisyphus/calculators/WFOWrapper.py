#!/usr/bin/env python3

import itertools
import logging
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
try:
    from pyscf import gto
except ImportError:
    print("Couldn't import pyscf!")
    #raise Exception
import pyparsing as pp

from pysisyphus.config import Config


CIOVL="""mix_aoovl=ao_ovl
a_mo=mos.1
b_mo=mos.2
a_det=dets.1
b_det=dets.2
a_mo_read=2
b_mo_read=2"""


class WFOWrapper:
    def __init__(self, occ_mos, virt_mos):
        self.base_cmd = Config["wfoverlap"]["cmd"]
        self.occ_mos = occ_mos
        self.virt_mos = virt_mos
        self.mos = self.occ_mos + self.virt_mos
        self.base_det_str = "d"*self.occ_mos + "e"*self.virt_mos
        self.fmt = "{: .10f}"

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
            mol.unit = "bohr"
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

    def make_full_dets_list(self, all_inds, det_strings, ci_coeffs):
        dets_list = list()
        for inds, det_string in zip(all_inds, det_strings):
            ab, ba = det_string
            from_mo, to_mo = inds
            per_state =  ci_coeffs[:,from_mo,to_mo]
            # See 10.1063/1.3000012 Eq. (5) and 10.1021/acs.jpclett.7b01479 SI
            per_state *= 1/2**0.5
            as_str = lambda arr: " ".join([self.fmt.format(cic)
                                           for cic in arr])
            ps_str = as_str(per_state)
            mps_str = as_str(-per_state)
            dets_list.append(f"{ab}\t{ps_str}")
            dets_list.append(f"{ba}\t{mps_str}")
        return dets_list

    def set_from_nested_list(self, nested):
        return set([i for i in itertools.chain(*nested)])

    def store_iteration(self, atoms, coords, mos, ci_coeffs, mo_inds):
        from_mos, to_mos = zip(*mo_inds)
        from_set = self.set_from_nested_list(from_mos)
        to_set = self.set_from_nested_list(to_mos)

        self.atoms = atoms
        self.mos_list.append(mos)
        self.coords_list.append(coords)
        self.ci_coeffs_list.append(ci_coeffs)
        self.mo_inds_list.append(mo_inds)
        self.from_set_list.append(from_set)
        self.to_set_list.append(to_set)

    def get_iteration(self, ind):
        return (self.mos_list[ind], self.coords_list[ind],
                self.ci_coeffs_list[ind], self.mo_inds_list[ind],
                self.from_set_list[ind], self.to_set_list[ind])

    def make_dets_header(self, cic, dets_list):
        return f"{len(cic)} {self.mos} {len(dets_list)}"

    def parse_wfoverlap(self, text):
        """Returns overlap matrix."""
        header = pp.Literal("Overlap matrix <PsiA_i|PsiB_j>")
        float_ = pp.Word(pp.nums+"-.")
        psi_bra = pp.Literal("<Psi") + pp.Word(pp.alphas) \
                  + pp.Word(pp.nums) + pp.Literal("|")
        psi_ket = pp.Literal("|Psi") + pp.Word(pp.alphas) \
                  + pp.Word(pp.nums) + pp.Literal(">")
        matrix_line = pp.Suppress(psi_bra) + pp.OneOrMore(float_)

        parser = pp.SkipTo(header, include=True) \
                 + psi_ket + psi_ket \
                 + pp.OneOrMore(matrix_line).setResultsName("overlap")

        result = parser.parseString(text)
        return np.array(list(result["overlap"]), dtype=np.float)

    def track(self, old_root=None):
        mos1, coords1, cic1, moi1, fs1, ts1 = self.get_iteration(-2)
        mos2, coords2, cic2, moi2, fs2, ts2 = self.get_iteration(-1)
        # Create a fake array for the ground state where all CI coefficients
        # are zero and add it.
        gs_cic = np.zeros_like(cic1[0])
        cic1_with_gs = np.concatenate((gs_cic[None,:,:], cic1))
        cic2_with_gs = np.concatenate((gs_cic[None,:,:], cic2))
        a1 = self.build_mole(coords1)
        a2 = self.build_mole(coords2)
        ao_ovlp = self.overlap(a1, a2, basis="sto-3g")
        ao_header = "{} {}".format(*ao_ovlp.shape)
        logging.warning("!using sto3g! for ao overlaps")

        all_inds, det_strings = self.generate_all_dets(fs1, ts1, fs2, ts2)
        # Prepare line for ground state
        gs_coeffs = np.zeros(len(cic1_with_gs))
        # Ground state is 100% HF configuration
        gs_coeffs[0] = 1
        gs_coeffs_str = " ".join([self.fmt.format(c)
                                  for c in gs_coeffs])
        gs_line = f"{self.base_det_str}\t{gs_coeffs_str}"
        dets1 = [gs_line] + self.make_full_dets_list(all_inds,
                                                     det_strings,
                                                     cic1_with_gs)
        dets2 = [gs_line] + self.make_full_dets_list(all_inds,
                                                     det_strings,
                                                     cic2_with_gs)
        header1 = self.make_dets_header(cic1_with_gs, dets1)
        header2 = self.make_dets_header(cic2_with_gs, dets2)

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            mos1_path = shutil.copy(mos1, tmp_path / "mos.1")
            mos2_path = shutil.copy(mos2, tmp_path / "mos.2")
            dets1_path = tmp_path / "dets.1"
            with open(dets1_path, "w") as handle:
                handle.write(header1+"\n"+"\n".join(dets1))
            dets2_path = tmp_path / "dets.2"
            with open(dets2_path, "w") as handle:
                handle.write(header2+"\n"+"\n".join(dets2))
            ao_ovl_path = tmp_path / "ao_ovl"
            np.savetxt(ao_ovl_path, ao_ovlp, fmt="%+22.15E", header=ao_header,
                       comments="")
            ciovl_fn = "ciovl.in"
            with open(tmp_path / ciovl_fn, "w") as handle:
                handle.write(CIOVL)
            cmd = f"{self.base_cmd} -f {ciovl_fn}".split()
            result = subprocess.Popen(cmd, cwd=tmp_path,
                                      stdout=subprocess.PIPE)
            result.wait()
            stdout = result.stdout.read().decode("utf-8")
        print(stdout)
        overlap_matrix = self.parse_wfoverlap(stdout)
        overlap_matrix = overlap_matrix.reshape(-1, len(cic2_with_gs))
        print(overlap_matrix)

        if old_root:
            return overlap_matrix[old_root]
        else:
            return overlap_matrix
