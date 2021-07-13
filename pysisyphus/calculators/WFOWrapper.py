from collections import OrderedDict
import itertools
import logging
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
import pyparsing as pp

from pysisyphus.config import Config
from pysisyphus.helpers_pure import chunks


CIOVL="""mix_aoovl=ao_ovl
a_mo=mos.1
b_mo=mos.2
ncore={ncore}
a_det=dets.1
b_det=dets.2
a_mo_read=2
b_mo_read=2
"""

CIOVL_NO_SAO="""ao_read=-1
same_aos=.true.
a_mo=mos.1
b_mo=mos.2
ncore={ncore}
a_det=dets.1
b_det=dets.2
a_mo_read=2
b_mo_read=2"""


class WFOWrapper:
    logger = logging.getLogger("wfoverlap")
    matrix_types = OrderedDict((
        ("ovlp", "Overlap matrix"),
        ("renorm", "Renormalized overlap matrix"),
        ("ortho", "Orthonormalized overlap matrix")
    ))

    def __init__(self, occ_mo_num, virt_mo_num, conf_thresh=1e-4,
                 calc_number=0, out_dir="./", wfow_mem=8000,
                 ncore=0, debug=False):
        try:
            self.base_cmd = Config["wfoverlap"]["cmd"]
        except KeyError:
            self.log("WFOverlap cmd not found in ~/.pysisyphusrc!")
        # Should correspond to the attribute of the parent calculator
        self.calc_number = calc_number
        self.conf_thresh = conf_thresh
        self.out_dir = Path(out_dir).resolve()
        self.wfow_mem = int(wfow_mem)
        self.ncore = int(ncore)
        self.debug = debug

        self.name = f"WFOWrapper_{self.calc_number}"
        self.log(f"Using -m {self.wfow_mem} for wfoverlap.")

        self.mo_inds_list = list()
        self.from_set_list = list()
        self.to_set_list = list()
        self.turbo_mos_list = list()

        self.occ_mo_num = int(occ_mo_num)
        self.virt_mo_num = int(virt_mo_num)
        self.mo_num = self.occ_mo_num + self.virt_mo_num
        self.base_det_str = "d"*self.occ_mo_num + "e"*self.virt_mo_num
        self.fmt = "{: .10f}"

        self.iter_counter = 0

    def log(self, message):
        self.logger.debug(f"{self.name}, " + message)

    def fake_turbo_mos(self, mo_coeffs):
        """Create a mos file suitable for TURBOMOLE input. All MO eigenvalues
        are set to 0.0. There is also a little deviation in the formatting
        (see turbo_fmt()) but it works ..."""

        def turbo_fmt(num):
            """Not quite the real TURBOMOLE format, but it works ...
            In TURBOMOLE the first character is always 0 for positive doubles
            and - for negative doubles."""
            return f"{num:+20.13E}".replace("E", "D")

        base = "$scfmo    scfconv=7  format(4d20.14)\n# from pysisyphus\n" \
               "{mo_strings}\n$end"

        # WFOverlap expects the string eigenvalue starting at 16, so we have
        mo_str = "{mo_index:>6d}  a      eigenvalue=-.00000000000000D+00   " \
                 "nsaos={nsaos}\n{joined}"
        nsaos = mo_coeffs.shape[0]

        mo_strings = list()
        for mo_index, mo in enumerate(mo_coeffs, 1):
            in_turbo_fmt = [turbo_fmt(c) for c in mo]
            # Combine into chunks of four
            lines = ["".join(chnk) for chnk in chunks(in_turbo_fmt, 4)]
            # Join the lines
            joined = "\n".join(lines)
            mo_strings.append(mo_str.format(mo_index=mo_index, nsaos=nsaos,
                                            joined=joined))
        return base.format(mo_strings="\n".join(mo_strings))

    def ci_coeffs_above_thresh(self, ci_coeffs, thresh=None):
        # Drop unimportant configurations, that are configurations
        # having low weights in all states under consideration.
        if thresh is None:
            thresh = self.conf_thresh
        mo_inds = np.where(np.abs(ci_coeffs) >= thresh)
        return mo_inds

    def make_det_string(self, inds):
        """Return spin adapted strings."""
        from_mo, to_mo = inds
        # Until now the first virtual MO (to_mo) has index 0. To subsitute
        # the base_str at the correct index we have to increase all to_mo
        # indices by the number off occupied MO.
        to_mo += self.occ_mo_num
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
            if (np.abs(per_state) < self.conf_thresh).all():
                continue
            # A singlet determinant can be formed in two ways:
            # (up down) (up down) (up down) ...
            # or
            # (down up) (down up) (down up) ...
            # We take this into account by expanding the singlet determinants
            # and using a proper normalization constant.
            # See 10.1063/1.3000012 Eq. (5) and 10.1021/acs.jpclett.7b01479 SI
            # and "Principles of Molecular Photochemistry: An Introduction",
            # Section 2.27 Vector model of Two Coupled Electron Spins, p. 91-92
            per_state *= 1/(2**0.5)
            as_str = lambda arr: " ".join([self.fmt.format(cic)
                                           for cic in arr])
            ps_str = as_str(per_state)
            mps_str = as_str(-per_state)
            dets_list.append(f"{ab}\t{ps_str}")
            dets_list.append(f"{ba}\t{mps_str}")
        return dets_list

    def set_from_nested_list(self, nested):
        return set([i for i in itertools.chain(*nested)])

    def make_dets_header(self, cic, dets_list):
        return f"{len(cic)} {self.mo_num} {len(dets_list)}"

    def parse_wfoverlap_out(self, text, type_="ortho"):
        """Returns overlap matrix."""
        header_str = self.matrix_types[type_] + " <PsiA_i|PsiB_j>"
        header = pp.Literal(header_str)
        float_ = pp.Word(pp.nums+"-.")
        psi_bra = pp.Literal("<Psi") + pp.Word(pp.alphas) \
                  + pp.Word(pp.nums) + pp.Literal("|")
        psi_ket = pp.Literal("|Psi") + pp.Word(pp.alphas) \
                  + pp.Word(pp.nums) + pp.Literal(">")
        matrix_line = pp.Suppress(psi_bra) + pp.OneOrMore(float_)

        # I really don't know why this is needed but otherwise I can't parse
        # overlap calculations with the true AO overlap matrix, even though
        # the files appear completely similar regarding printing of the matrices.
        # WTF. WTF!
        text = text.replace("\n", " ")
        parser = pp.SkipTo(header, include=True) \
                 + pp.OneOrMore(psi_ket) \
                 + pp.OneOrMore(matrix_line).setResultsName("overlap")

        result = parser.parseString(text)

        return np.array(list(result["overlap"]), dtype=np.float)

    def get_from_to_sets(self, ci_coeffs):
        all_mo_inds = [self.ci_coeffs_above_thresh(per_state)
                       for per_state in ci_coeffs]

        from_mos, to_mos = zip(*all_mo_inds)
        from_set = self.set_from_nested_list(from_mos)
        to_set = self.set_from_nested_list(to_mos)

        return from_set, to_set

    def get_gs_line(self, ci_coeffs_with_gs):
        gs_coeffs = np.zeros(len(ci_coeffs_with_gs))
        # Ground state is 100% HF configuration
        gs_coeffs[0] = 1
        gs_coeffs_str = " ".join([self.fmt.format(c)
                                  for c in gs_coeffs])
        gs_line = f"{self.base_det_str}\t{gs_coeffs_str}"
        return gs_line

    def wf_overlap(self, cycle1, cycle2, ao_ovlp=None):
        mos1, cic1 = cycle1
        mos2, cic2 = cycle2

        fs1, ts1 = self.get_from_to_sets(cic1)
        fs2, ts2 = self.get_from_to_sets(cic2)

        # Create a fake array for the ground state where all CI coefficients
        # are zero and add it.
        gs_cic = np.zeros_like(cic1[0])
        cic1_with_gs = np.concatenate((gs_cic[None,:,:], cic1))
        cic2_with_gs = np.concatenate((gs_cic[None,:,:], cic2))

        all_inds, det_strings = self.generate_all_dets(fs1, ts1, fs2, ts2)
        # Prepare lines for ground state
        gs_line1 = self.get_gs_line(cic1_with_gs)
        gs_line2 = self.get_gs_line(cic2_with_gs)
        dets1 = [gs_line1] + self.make_full_dets_list(all_inds,
                                                     det_strings,
                                                     cic1_with_gs)
        dets2 = [gs_line2] + self.make_full_dets_list(all_inds,
                                                     det_strings,
                                                     cic2_with_gs)
        header1 = self.make_dets_header(cic1_with_gs, dets1)
        header2 = self.make_dets_header(cic2_with_gs, dets2)

        backup_path = self.out_dir / f"wfo_{self.calc_number}.{self.iter_counter:03d}"
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            self.log(f"Calculation in {tmp_dir}")
            # Write fake TURBOMOLE mo files
            for i, mos in enumerate((mos1, mos2), 1):
                turbo_mos = self.fake_turbo_mos(mos)
                with open(tmp_path / f"mos.{i}", "w") as handle:
                    handle.write(turbo_mos)
            dets1_path = tmp_path / "dets.1"
            with open(dets1_path, "w") as handle:
                handle.write(header1+"\n"+"\n".join(dets1))
            dets2_path = tmp_path / "dets.2"
            with open(dets2_path, "w") as handle:
                handle.write(header2+"\n"+"\n".join(dets2))

            # Decide wether to use a double molecule overlap matrix or
            # (approximately) reconstruct the ao_ovlp matrix from the MO
            # coefficients.
            if ao_ovlp is None:
                ciovl_in = CIOVL_NO_SAO
                self.log("Got no ao_ovl-matrix. Using ao_read=-1 and "
                         "same_aos=.true. to reconstruct the AO-overlap matrix!")
            else:
                ciovl_in = CIOVL
                ao_header = "{} {}".format(*ao_ovlp.shape)
                ao_ovl_path = tmp_path / "ao_ovl"
                np.savetxt(ao_ovl_path, ao_ovlp, fmt="%22.15E", header=ao_header,
                           comments="")

            ciovl_in_rendered = ciovl_in.format(ncore=self.ncore)
            ciovl_fn = "ciovl.in"
            with open(tmp_path / ciovl_fn, "w") as handle:
                handle.write(ciovl_in_rendered)

            # Create a backup of the whole temporary directory
            try:
                shutil.rmtree(backup_path)
            except FileNotFoundError:
                pass
            shutil.copytree(tmp_dir, backup_path)

            # Currently, debug==True crashes the subsequent parsing
            debug_str = "--debug" if self.debug else ""
            cmd = f"{self.base_cmd} -m {self.wfow_mem} -f {ciovl_fn} {debug_str}".split()
            result = subprocess.Popen(cmd, cwd=tmp_path,
                                      stdout=subprocess.PIPE)
            result.wait()
            stdout = result.stdout.read().decode("utf-8")
            self.iter_counter += 1
        if "differs significantly" in stdout:
            self.log("WARNING: Orthogonalized matrix differs significantly "
                     "from original matrix! There is probably mixing with "
		     "external states.")

        wfo_log_fn = self.out_dir / f"wfo_{self.calc_number}.{self.iter_counter:03d}.out"
        with open(wfo_log_fn, "w") as handle:
            handle.write(stdout)
        # Also copy the WFO-output to the input backup
        shutil.copy(wfo_log_fn, backup_path)

        matrices = [self.parse_wfoverlap_out(stdout, type_=key)
                    for key in self.matrix_types.keys()]

        reshaped_mats = [mat.reshape(-1, len(cic2_with_gs))
                         for mat in matrices]
        for key, mat in zip(self.matrix_types.keys(), reshaped_mats):
            mat_fn = backup_path / f"{key}_mat.dat"
            np.savetxt(mat_fn, mat)

        return reshaped_mats

    def __str__(self):
        return self.name
