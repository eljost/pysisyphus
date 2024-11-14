from collections import OrderedDict
import itertools
import logging
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Literal

import numpy as np
import pyparsing as pp

from pysisyphus.config import Config
from pysisyphus.helpers_pure import chunks


# Tempate for cases w/ explicit overlap matrix
CIOVL = """mix_aoovl=ao_ovl
a_mo=mos.1
b_mo=mos.2
ncore={ncore}
a_det=dets.1
b_det=dets.2
a_mo_read=2
b_mo_read=2
"""

# Tempate for cases w/o overlap matrix. wfoverlap will reconstruct it.
CIOVL_NO_SAO = """ao_read=-1
same_aos=.true.
a_mo=mos.1
b_mo=mos.2
ncore={ncore}
a_det=dets.1
b_det=dets.2
a_mo_read=2
b_mo_read=2"""


# Dictionary that mimics a creation operator.
# The key comprises two entries (curr. MO, spin of electron to be added)
# Will return None for invalid applications
# 2: doubly occupeid, a: singly occuped w/ spin alpha, b: singly occupied w/ spin beta
ADD_TO_MO = {
    # Add alpha
    ("a", "a"): None,
    ("b", "a"): "2",
    ("e", "a"): "a",
    ("2", "a"): None,
    # Add beta
    ("a", "b"): "2",
    ("b", "b"): None,
    ("e", "b"): "b",
    ("2", "b"): None,
}

# Dictionary that mimics an annihilation operator.
# The key comprises two entries (curr. MO, spin of electron to be deleted)
# See comment on ADD_TO_MO.
DEL_FROM_MO = {
    # Delete alpha
    ("a", "a"): "e",
    ("b", "a"): None,
    ("e", "a"): None,
    ("2", "a"): "b",
    # Delete beta
    ("a", "b"): None,
    ("b", "b"): "e",
    ("e", "b"): None,
    ("2", "b"): "a",
}


class WFOWrapper:
    logger = logging.getLogger("wfoverlap")
    matrix_types = OrderedDict(
        (
            ("ovlp", "Overlap matrix"),
            ("renorm", "Renormalized overlap matrix"),
            ("ortho", "Orthonormalized overlap matrix"),
        )
    )

    def __init__(
        self,
        occa: int,
        virta: int,
        occb: int,
        virtb: int,
        conf_thresh: float = 1e-3,
        calc_number: int = 0,
        out_dir: str | Path = "./",
        wfow_mem: int = 8000,
        ncore: int = 0,
        debug: bool = False,
    ):
        try:
            self.base_cmd = Config["wfoverlap"]["cmd"]
        except KeyError:
            self.log("WFOverlap cmd not found in ~/.pysisyphusrc!")
        # Should correspond to the attribute of the parent calculator
        self.calc_number = calc_number
        self.name = f"WFOWrapper_{self.calc_number}"

        self.conf_thresh = conf_thresh
        self.out_dir = Path(out_dir).resolve()
        self.wfow_mem = int(wfow_mem)
        self.ncore = int(ncore)
        self.debug = debug

        self.log(f"Using -m {self.wfow_mem} for wfoverlap.")

        self.mo_inds_list = list()
        self.from_set_list = list()
        self.to_set_list = list()
        self.turbo_mos_list = list()

        self.occa = int(occa)
        self.virta = int(virta)
        self.occb = int(occb)
        self.virtb = int(virtb)
        self.nmosa = self.occa + self.virta
        self.nmosb = self.occb + self.virtb
        self.nmos = self.nmosa + self.nmosb

        # We begin with alpha MOs, followed by beta MOs
        self.base_det_str = (
            # Alpha
            "a" * self.occa
            + "e" * self.virta
            # Beta
            + "b" * self.occb
            + "e" * self.virtb
        )

        self.fmt = "{: .10f}"
        self.iter_counter = 0

    @property
    def conf_thresh(self):
        return self._conf_thresh

    @conf_thresh.setter
    def conf_thresh(self, conf_thresh):
        self._conf_thresh = conf_thresh
        self.log(f"Set CI-coeff threshold to {self.conf_thresh:.4e}")

    def log(self, message):
        self.logger.debug(f"{self.name}, " + message)

    def fake_turbo_mos(self, mo_coeffs):
        """Create a mos file suitable for TURBOMOLE input. All MO eigenvalues
        are set to 0.0. There is also a little deviation in the formatting
        (see turbo_fmt()) but it works ...

        Supplied MO coefficients are expected in the shape (naos, nmos), that
        is MOs must be supplied in columns!
        """

        def turbo_fmt(num):
            """Not quite the real TURBOMOLE format, but it works ...
            In TURBOMOLE the first character is always 0 for positive doubles
            and - for negative doubles."""
            return f"{num:+20.13E}".replace("E", "D")

        base = (
            "$scfmo    scfconv=7  format(4d20.14)\n# from pysisyphus\n"
            "{mo_strings}\n$end"
        )

        # WFOverlap expects the string eigenvalue starting at 16, so we have
        mo_str = (
            "{mo_index:>6d}  a      eigenvalue=-.00000000000000D+00   "
            "nsaos={nsaos}\n{joined}"
        )
        nsaos, _ = mo_coeffs.shape

        mo_strings = list()
        # Loop over MOs; columns of 'mo_coeffs'
        for mo_index, mo in enumerate(mo_coeffs.T, 1):
            in_turbo_fmt = [turbo_fmt(c) for c in mo]
            # Combine into chunks of four
            lines = ["".join(chnk) for chnk in chunks(in_turbo_fmt, 4)]
            # Join the lines
            joined = "\n".join(lines)
            mo_strings.append(
                mo_str.format(mo_index=mo_index, nsaos=nsaos, joined=joined)
            )
        return base.format(mo_strings="\n".join(mo_strings))

    def ci_coeffs_above_thresh(self, ci_coeffs, thresh=None):
        # Drop unimportant configurations, that are configurations
        # having low weights in all states under consideration.
        if thresh is None:
            thresh = self.conf_thresh
        mo_inds = np.where(np.abs(ci_coeffs) >= thresh)
        return mo_inds

    def make_det_string(self, inds, exc_of: Literal["a", "b"]):
        """Return spin adapted strings."""
        from_mo, to_mo = inds
        # Until now the first virtual MO (to_mo) has index 0. To subsitute
        # the base_str at the correct index we have to increase all to_mo
        # indices by the number off occupied MO.
        det = list(self.base_det_str)
        # Excitation of an alpha electron
        if exc_of == "a":
            to_mo += self.occa
        # Excitation of a beta electron
        elif exc_of == "b":
            from_mo += self.nmosa
            to_mo += self.nmosa + self.occb
        from_key = (det[from_mo], exc_of)
        det[from_mo] = DEL_FROM_MO[from_key]
        add_key = (det[to_mo], exc_of)
        det[to_mo] = ADD_TO_MO[add_key]
        det = "".join(det)
        return det

    def generate_all_dets(
        self, occ_set1, virt_set1, occ_set2, virt_set2, exc_of: Literal["a", "b"]
    ):
        """Generate all possible single excitation determinant strings
        from union(occ_mos) to union(virt_mos)."""
        # Unite the respective sets of both calculations
        occ_set = occ_set1 | occ_set2
        virt_set = virt_set1 | virt_set2
        # Genrate all possible excitations (combinations) from the occupied
        # MO set to (and) the virtual MO set.
        all_inds = [(om, vm) for om, vm in itertools.product(occ_set, virt_set)]
        det_strings = [self.make_det_string(inds, exc_of) for inds in all_inds]
        return all_inds, det_strings

    """
    def make_full_dets_list(self, all_inds, det_strings, ci_coeffs):
        dets_list = list()
        for inds, det_string in zip(all_inds, det_strings):
            ab, ba = det_string
            from_mo, to_mo = inds
            per_state = ci_coeffs[:, from_mo, to_mo]
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
            per_state *= 1 / (2**0.5)
            as_str = lambda arr: " ".join([self.fmt.format(cic) for cic in arr])
            ps_str = as_str(per_state)
            mps_str = as_str(-per_state)
            dets_list.append(f"{ab}\t{ps_str}")
            dets_list.append(f"{ba}\t{mps_str}")
        return dets_list
    """

    def make_full_dets_list_unr(self, all_inds, det_strings, ci_coeffs):
        dets_list = list()
        for inds, det in zip(all_inds, det_strings):
            from_mo, to_mo = inds
            # CI-coefficients for a given MO pair for all states
            per_state = ci_coeffs[:, from_mo, to_mo]
            # Skip this determinant when the coefficients are very small for all states
            if (np.abs(per_state) < self.conf_thresh).all():
                continue
            coeff_str = " ".join([self.fmt.format(cic) for cic in per_state])
            dets_list.append(f"{det}\t{coeff_str}")
        return dets_list

    def set_from_nested_list(self, nested):
        return set([i for i in itertools.chain(*nested)])

    def make_dets_header(self, nstates, nmos, ndets):
        return f"{nstates} {nmos} {ndets}"

    def parse_wfoverlap_out(self, text, type_="ortho"):
        """Returns overlap matrix."""
        header_str = self.matrix_types[type_] + " <PsiA_i|PsiB_j>"
        header = pp.Literal(header_str)
        float_ = pp.Word(pp.nums + "-.")
        psi_bra = (
            pp.Literal("<Psi") + pp.Word(pp.alphas) + pp.Word(pp.nums) + pp.Literal("|")
        )
        psi_ket = (
            pp.Literal("|Psi") + pp.Word(pp.alphas) + pp.Word(pp.nums) + pp.Literal(">")
        )
        matrix_line = pp.Suppress(psi_bra) + pp.OneOrMore(float_)

        # I really don't know why this is needed but otherwise I can't parse
        # overlap calculations with the true AO overlap matrix, even though
        # the files appear completely similar regarding printing of the matrices.
        # WTF. WTF!
        text = text.replace("\n", " ")
        parser = (
            pp.SkipTo(header, include=True)
            + pp.OneOrMore(psi_ket)
            + pp.OneOrMore(matrix_line).setResultsName("overlap")
        )

        result = parser.parseString(text)

        return np.array(list(result["overlap"]), dtype=np.float64)

    def get_from_to_sets(self, ci_coeffs):
        all_mo_inds = [
            self.ci_coeffs_above_thresh(per_state) for per_state in ci_coeffs
        ]

        from_mos, to_mos = zip(*all_mo_inds)
        from_set = self.set_from_nested_list(from_mos)
        to_set = self.set_from_nested_list(to_mos)

        return from_set, to_set

    def get_gs_line(self, ci_coeffs_with_gs):
        gs_coeffs = np.zeros(len(ci_coeffs_with_gs))
        # Ground state is 100% HF configuration
        gs_coeffs[0] = 1
        gs_coeffs_str = " ".join([self.fmt.format(c) for c in gs_coeffs])
        gs_line = f"{self.base_det_str}\t{gs_coeffs_str}"
        return gs_line

    def wf_overlap(self, mos1, cic1a, cic1b, mos2, cic2a, cic2b, ao_ovlp=None):
        assert mos1.shape[1] == mos2.shape[1]
        # Add ground states as 'states' where all CI-coefficients are 0.0.
        # Appropriate determinant strings are constructed later
        cic1a_gs = np.zeros_like(cic1a[0])
        cic1b_gs = np.zeros_like(cic1b[0])
        cic2a_gs = np.zeros_like(cic2a[0])
        cic2b_gs = np.zeros_like(cic2b[0])
        cic1a = np.concatenate((cic1a_gs[None, :, :], cic1a))
        cic1b = np.concatenate((cic1b_gs[None, :, :], cic1b))
        cic2a = np.concatenate((cic2a_gs[None, :, :], cic2a))
        cic2b = np.concatenate((cic2b_gs[None, :, :], cic2b))

        nstates1 = len(cic1a)
        nstates2 = len(cic2a)

        # Determine which orbital pairs (from, to) are involved/important
        #
        # Reference
        fs1a, ts1a = self.get_from_to_sets(cic1a)
        fs1b, ts1b = self.get_from_to_sets(cic1b)
        # Current
        fs2a, ts2a = self.get_from_to_sets(cic2a)
        fs2b, ts2b = self.get_from_to_sets(cic2b)

        # Create a list of strings containing all possible determinants
        all_indsa, det_stringsa = self.generate_all_dets(fs1a, ts1a, fs2a, ts2a, "a")
        all_indsb, det_stringsb = self.generate_all_dets(fs1b, ts1b, fs2b, ts2b, "b")

        # Create GS determinant and determinants for alpha and beta excitations
        gs_line1 = self.get_gs_line(cic1a)
        dets1a = self.make_full_dets_list_unr(all_indsa, det_stringsa, cic1a)
        dets1b = self.make_full_dets_list_unr(all_indsb, det_stringsb, cic1b)
        dets1 = [gs_line1] + dets1a + dets1b

        # Create GS determinant and determinants for alpha and beta excitations
        gs_line2 = self.get_gs_line(cic2a)
        dets2a = self.make_full_dets_list_unr(all_indsa, det_stringsa, cic2a)
        dets2b = self.make_full_dets_list_unr(all_indsb, det_stringsb, cic2b)
        dets2 = [gs_line2] + dets2a + dets2b

        header1 = self.make_dets_header(len(cic1a), self.nmos, len(dets1))
        header2 = self.make_dets_header(len(cic2a), self.nmos, len(dets2))

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
                handle.write(header1 + "\n" + "\n".join(dets1))
            dets2_path = tmp_path / "dets.2"
            with open(dets2_path, "w") as handle:
                handle.write(header2 + "\n" + "\n".join(dets2))

            # Decide wether to use a double molecule overlap matrix or
            # (approximately) reconstruct the ao_ovlp matrix from the MO
            # coefficients.
            if ao_ovlp is None:
                ciovl_in = CIOVL_NO_SAO
                self.log(
                    "Got no ao_ovl-matrix. Using ao_read=-1 and "
                    "same_aos=.true. to reconstruct the AO-overlap matrix!"
                )
            else:
                ciovl_in = CIOVL
                ao_header = "{} {}".format(*ao_ovlp.shape)
                ao_ovl_path = tmp_path / "ao_ovl"
                np.savetxt(
                    ao_ovl_path, ao_ovlp, fmt="%22.15E", header=ao_header, comments=""
                )

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
            cmd = (
                f"{self.base_cmd} -m {self.wfow_mem} -f {ciovl_fn} {debug_str}".split()
            )
            result = subprocess.Popen(cmd, cwd=tmp_path, stdout=subprocess.PIPE)
            result.wait()
            stdout = result.stdout.read().decode("utf-8")
            self.iter_counter += 1
        if "differs significantly" in stdout:
            self.log(
                "WARNING: Orthogonalized matrix differs significantly "
                "from original matrix! There is probably mixing with "
                "external states."
            )

        wfo_log_fn = (
            self.out_dir / f"wfo_{self.calc_number}.{self.iter_counter:03d}.out"
        )
        with open(wfo_log_fn, "w") as handle:
            handle.write(stdout)
        # Also copy the WFO-output to the input backup
        shutil.copy(wfo_log_fn, backup_path)

        matrices = [
            self.parse_wfoverlap_out(stdout, type_=key)
            for key in self.matrix_types.keys()
        ]

        reshaped_mats = [mat.reshape(nstates1, nstates2) for mat in matrices]
        for key, mat in zip(self.matrix_types.keys(), reshaped_mats):
            mat_fn = backup_path / f"{key}_mat.dat"
            np.savetxt(mat_fn, mat)

        return reshaped_mats

    def __str__(self):
        return self.name
