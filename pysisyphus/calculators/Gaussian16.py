from collections import namedtuple
import io
from pathlib import Path
import re
import shutil
import subprocess
import textwrap

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.OverlapCalculator import OverlapCalculator
from pysisyphus.constants import AU2EV, BOHR2ANG
from pysisyphus.helpers_pure import file_or_str


class Gaussian16(OverlapCalculator):

    conf_key = "gaussian16"

    def __init__(
        self,
        route,
        gbs="",
        gen="",
        keep_chk=False,
        stable="",
        fchk=None,
        **kwargs,
    ):
        super().__init__(**kwargs)

        self.route = route.lower()
        invalid_keywords = ("symmetry", "nosymm", "force", "opt", "freq", "irc")
        invalid_kws_given = [kw for kw in invalid_keywords if kw in self.route]
        assert (
            not invalid_kws_given
        ), f"Invalid keywords given in route: {invalid_kws_given}. Please remove them."
        self.gbs = gbs
        assert "@" not in gbs, "Give only the path to the .gbs file, " "without the @!"
        self.gen = gen
        self.keep_chk = keep_chk
        self.stable = stable
        self.fchk = fchk

        keywords = {
            kw: option
            for kw, option in [self.parse_keyword(kw) for kw in self.route.split()]
        }
        exc_keyword = [key for key in "td tda cis".split() if key in keywords]
        self.root = None
        self.nstates = None
        if exc_keyword:
            self.exc_key = exc_keyword[0]
            exc_dict = keywords[self.exc_key]
            self.nstates = int(exc_dict["nstates"])
            try:
                self.root = int(exc_dict["root"])
            except KeyError:
                self.root = 1
                self.log("No explicit root was specified! Using root=1 as default!")
            # Collect remaining options if specified
            self.exc_args = {
                k: v for k, v in exc_dict.items() if k not in ("nstates", "root")
            }
            # Delete exc keyword, as we build it later on
            self.route = re.sub(r"((?:td|cis|tda).+?(:?\s|$))", "", self.route)

        self.to_keep = ("com", "fchk", "log", "dump_635r", "input.xyz")
        if self.keep_chk:
            self.to_keep += (".chk",)

        self.fn_base = "gaussian16"
        self.inp_fn = f"{self.fn_base}.com"
        self.out_fn = f"{self.fn_base}.log"
        self.chk_fn = f"{self.fn_base}.chk"
        self.dump_base_fn = f"{self.fn_base}_rwfdump"

        self.gaussian_input = """
        %nproc={pal}
        %mem={mem}MB
        {chk_link0}
        {add_link0}
        #P {calc_type} {route} {exc}
        # Symmetry=None {reuse_data}
        
        title

        {charge} {mult}
        {coords}

        {gen}{gbs}



        """
        self.gaussian_input = textwrap.dedent(self.gaussian_input)

        self.parser_funcs = {
            "force": self.parse_force,
            "hessian": self.parse_hessian,
            "stable": self.parse_stable,
            "noparse": lambda path: None,
            "double_mol": self.parse_double_mol,
        }

        self.base_cmd = self.get_cmd("cmd")
        self.formchk_cmd = self.get_cmd("formchk")
        self.unfchk_cmd = self.get_cmd("unfchk")
        self.rwfdump_cmd = self.get_cmd("rwfdump")
        self.log(
            f"Using commands: g16={self.base_cmd}, formchk={self.formchk_cmd}, "
            f"unfchk={self.unfchk_cmd}, rwfdump={self.rwfdump_cmd}"
        )

    def make_exc_str(self):
        # Ground state calculation
        if not self.root:
            return ""
        root = f"root={self.root}"
        nstates = f"nstates={self.nstates}"
        pair2str = lambda k, v: f"{k}" + (f"={v}" if v else "")
        arg_str = ",".join([pair2str(k, v) for k, v in self.exc_args.items()])
        exc_str = f"{self.exc_key}=({root},{nstates},{arg_str})"
        return exc_str

    def reuse_data(self, path):
        # Nothing to reuse if no fchk or chk present
        if (self.fchk is None) and not hasattr(self, "chk"):
            return ""
        new_chk = path / self.chk_fn
        prev_fchk = new_chk.with_suffix(".fchk")
        shutil.copy(self.fchk, prev_fchk)
        cmd = f"{self.unfchk_cmd} {prev_fchk}".split()
        subprocess.run(cmd, stdout=subprocess.PIPE, cwd=path)
        self.log(f"Using MO guess from '{self.fchk}'.")

        reuse_str = "guess=read"
        # Also try to reuse information of previous TD calculation
        # if self.nstates and hasattr(self, "chk"):
        # DISABLED for now, as gaussian also resorts the states
        # internally, which then messes up our internal numbering
        # of states.
        # shutil.copy(self.chk, new_chk)
        # reuse_str += " td=read"
        # self.log("Using td=read")

        return reuse_str

    def make_gbs_str(self):
        if self.gbs:
            return f"@{self.gbs}"
        else:
            return ""

    def prepare_input(
        self, atoms, coords, calc_type, did_stable=False, point_charges=None
    ):
        path = self.prepare_path(use_in_run=True)
        xyz_str = self.prepare_xyz_string(atoms, coords)
        with open(path / "input.xyz", "w") as handle:
            handle.write(xyz_str)

        coords = self.prepare_coords(atoms, coords)
        kwargs = {
            "pal": self.pal,
            "mem": self.pal * self.mem,
            "chk_link0": f"%chk={self.chk_fn}",
            "add_link0": "",
            "route": self.route,
            # td/tda/cis(...)
            "exc": self.make_exc_str(),
            # guess=read td=read
            "reuse_data": self.reuse_data(path),
            "calc_type": calc_type,
            "charge": self.charge,
            "mult": self.mult,
            "coords": coords,
            "gbs": self.make_gbs_str(),
            "gen": self.gen,
        }
        if calc_type == "double_mol":
            update = {
                "chk_link0": "",
                "add_link0": "%KJob L302 1",
                "route": self.route + " iop(3/33=1) geom=notest",
                "calc_type": "",
                "reuse_data": "",
                "exc": "",
            }
            kwargs.update(update)
        # stable and td/tda can't be used together, so we omit the
        # td/tda string here
        if "stable" in calc_type.lower():
            kwargs["exc"] = ""
        # after a stability analysis we already got a converged wavefuntion
        # and we dont want to 'optimize' it any further
        if did_stable:
            kwargs["route"] += " scf=maxcycles=0"
        if point_charges is not None:
            kwargs["route"] += " charge"
            point_charges = point_charges.copy()
            point_charges[:, :3] *= BOHR2ANG
            with io.StringIO() as handle:
                np.savetxt(handle, point_charges, fmt="%16.10f")
                pc_str = handle.getvalue()
            # Append point charges to coords
            kwargs["coords"] += "\n\n" + pc_str
        inp = self.gaussian_input.format(**kwargs)
        return inp

    def make_fchk(self, path):
        cmd = f"{self.formchk_cmd} {self.chk_fn}".split()
        subprocess.run(cmd, stdout=subprocess.PIPE, cwd=path)
        self.log("Created .fchk")

    def run_rwfdump(self, path, rwf_index, chk_path=None):
        if chk_path is None:
            chk_path = path / self.chk_fn
        dump_fn = path / f"{self.dump_base_fn}_dump_{rwf_index}"
        cmd = f"rwfdump {chk_path} {dump_fn} {rwf_index}".split()
        subprocess.run(cmd)
        self.log(f"Dumped {rwf_index} from .chk.")
        return dump_fn

    def run_after(self, path):
        chk_path = path / self.chk_fn
        if not chk_path.exists():
            self.log("No .chk file found.")
            return
        # Create the .fchk file so we can keep it and parse it later on.
        self.make_fchk(path)
        if self.track:
            self.run_rwfdump(path, "635r")
            self.nmos, self.roots = self.parse_log(path / self.out_fn)

    def parse_keyword(self, text):
        word = pp.Word(pp.alphanums + "-" + "/")

        keyword = word.setResultsName("keyword")
        equals = pp.Literal("=")
        option = pp.Group(
            word
            + pp.Suppress(pp.Optional(equals))
            + pp.Optional(word, default="")
            + pp.Suppress(pp.Optional(","))
        )
        options = (
            pp.Suppress(pp.Optional(equals))
            + pp.Suppress(pp.Optional("("))
            + pp.OneOrMore(option)
            + pp.Suppress(pp.Optional(")"))
        ).setResultsName("options")

        parser = keyword + pp.Optional(options, default=[])

        result = parser.parseString(text)
        as_dict = result.asDict()
        kw = as_dict["keyword"]
        opt_dict = {key: value for key, value in as_dict["options"]}
        return kw, opt_dict

    @staticmethod
    def parse_fchk(fchk_path, keys):
        with open(fchk_path) as handle:
            text = handle.read()

        def to_float(s, loc, toks):
            return float(toks[0])

        # Matches -4.10693837E-16 and 1.60184209E-15
        float_ = pp.Word(pp.nums + "+-E.").setParseAction(to_float)
        # Start with Empty so we can progessively build the
        # parser for all keys.
        parser = pp.Empty()

        def parser_for_key(key):
            return pp.Group(
                pp.Suppress(pp.SkipTo(key))
                + key
                + pp.restOfLine
                + pp.ZeroOrMore(float_)
            )

        for key in keys:
            parser += parser_for_key(key)
        results = parser.parseString(text)
        results_dict = {}
        for key, res in zip(keys, results):
            # This handles scalar entries like
            # Total Energy  [...] R     -9.219940072302333E+01
            if len(res) == 2:
                results_dict[key] = float(res[-1].split()[-1])
            # This handles matrices like
            # Cartesian Gradient [...] R   N=           9 \
            # [Matrix entries]
            if len(res) > 2:
                results_dict[key] = np.array(res[2:])
        return results_dict

    def get_energy(self, atoms, coords, **prepare_kwargs):
        results = self.get_forces(atoms, coords, **prepare_kwargs)
        del results["forces"]
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        did_stable = False
        if self.stable:
            is_stable = self.run_stable(atoms, coords)
            self.log(f"Wavefunction is now stable: {is_stable}")
            did_stable = True
        inp = self.prepare_input(
            atoms, coords, "force", did_stable=did_stable, **prepare_kwargs
        )
        kwargs = {
            "calc": "force",
        }
        results = self.run(inp, **kwargs)
        if self.track:
            self.store_overlap_data(atoms, coords)
            if self.track_root():
                # Redo the calculation with the updated root
                results = self.get_forces(atoms, coords)
        return results

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        inp = self.prepare_input(atoms, coords, "freq", **prepare_kwargs)
        kwargs = {
            "calc": "hessian",
        }
        results = self.run(inp, **kwargs)
        return results

    def run_stable(self, atoms, coords, **prepare_kwargs):
        inp = self.prepare_input(atoms, coords, self.stable, **prepare_kwargs)
        self.log(f"Running stability analysis with {self.stable}")
        kwargs = {
            "calc": "stable",
        }
        results = self.run(inp, **kwargs)
        return results

    def parse_stable(self, path):
        log_path = path / self.out_fn
        with open(log_path) as handle:
            text = handle.read()
        instab_re = "wavefunction.*instability.*"
        instab_mobj = re.search(instab_re, text)
        if instab_mobj:
            self.log("Found instability!")
        stable_re = "wavefunction is stable"
        mobj = re.search(stable_re, text)
        is_stable = bool(mobj)
        return is_stable

    def run_calculation(self, atoms, coords, **prepare_kwargs):
        inp = self.prepare_input(atoms, coords, "", **prepare_kwargs)
        kwargs = {
            "calc": "noparse",
        }
        results = self.run(inp, **kwargs)
        if self.track:
            self.store_overlap_data(atoms, coords)
            self.track_root()
            self.log(
                "This track_root() call is a bit superfluous as the "
                "as the result is ignored :)"
            )
        return results

    def run_double_mol_calculation(self, atoms, coords1, coords2):
        self.log("Running double molecule calculation")
        double_atoms = atoms + atoms
        double_coords = np.hstack((coords1, coords2))
        inp = self.prepare_input(double_atoms, double_coords, "double_mol")
        kwargs = {
            "calc": "double_mol",
            "keep": False,
            "inc_counter": False,
            "run_after": False,
        }
        double_mol_ovlps = self.run(inp, **kwargs)
        return double_mol_ovlps

    def parse_tddft(self, path):
        with open(path / self.out_fn) as handle:
            text = handle.read()
        td_re = r"Excited State\s*\d+:\s*[\.\w\?-]+\s*([\d\.-]+?)\s*eV"
        matches = re.findall(td_re, text)
        assert len(matches) == self.nstates
        # Excitation energies in eV
        exc_energies = np.array(matches, dtype=np.float)
        # Convert to Hartree
        exc_energies /= AU2EV
        return exc_energies

    @file_or_str(".log", method=True)
    def parse_log(self, text):
        def parse(text, regex, func):
            mobj = re.search(regex, text)
            return func(mobj[1])

        # Depending on wether we did the calculation with td=read or not
        # roots will be at a different value. Without reading the CI coeffs
        # from the checkpoint Gaussian will calculate four times as much roots
        # as requested in the first iterations of the calculation. This will
        # lead to a much higher number of expected number of CI-coefficients
        # when parsing the 635r dump later on.
        roots_re = r"Root\s+(\d+)"
        roots = np.array(re.findall(roots_re, text), dtype=int).max()

        # NBasis=    16 NAE=    12 NBE=    12 NFC=     6 NFV=     0
        basis_re = "NBasis=(.+)NAE=(.+)NBE=(.+)NFC=(.+)NFV=(.+)"
        basis_mobj = re.search(basis_re, text)
        basis_funcs, alpha, beta, _, _ = [int(n) for n in basis_mobj.groups()]
        a_occ = alpha
        b_occ = beta
        a_vir = basis_funcs - a_occ
        b_vir = basis_funcs - b_occ
        restricted = alpha == beta
        act_re = "NROrb=(.*)NOA=(.*)NOB=(.*)NVA=(.*)NVB=(.*)"
        act_mobj = re.search(act_re, text)
        _, a_act, b_act, _, _ = [int(n) for n in act_mobj.groups()]
        a_core = a_occ - a_act
        b_core = b_occ - b_act

        NMOs = namedtuple(
            "NMOs",
            "a_core, a_occ a_act a_vir " "b_core b_occ b_act b_vir " "restricted",
        )

        nmos = NMOs(
            a_core, a_occ, a_act, a_vir, b_core, b_occ, b_act, b_vir, restricted
        )
        self.log(str(nmos))
        return nmos, roots

    def parse_635r_dump(self, dump_path, roots, nmos):
        self.log(f"Parsing 635r dump '{dump_path}'")
        with open(dump_path) as handle:
            text = handle.read()
        regex = r"read left to right\):\s*(.+)"
        mobj = re.search(regex, text, re.DOTALL)
        arr_str = mobj[1].replace("D", "E")
        # Drop the first 12 items as they are always 0
        tmp = np.array(arr_str.split()[12:], dtype=np.float64)

        # the core electrons are frozen in TDDFT/TDA
        expected = (nmos.a_act * nmos.a_vir + nmos.b_act * nmos.b_vir) * roots * 2
        self.log(
            f"Expecting {expected} ci coefficients in {dump_path}. "
            f"There are {tmp.size} items (including eigenvalues)."
        )
        coeffs = tmp[:expected]
        # 1. dim: X+Y, X-Y -> 2
        # 2. dim: roots -> variable, usually higher than Nstates as gaussian
        #         calculates more roots in the first iteration of TDDFT
        # 3. dim: alpha, beta -> 2
        # The remainder depends on the number of alpha and beta electrons.
        # The vectors are given for active orbs. x virt. obs that is the
        # core orbitals are neglected.
        if nmos.restricted:
            XpY, XmY = coeffs.reshape(2, roots, 2, -1)
            X = 0.5 * (XpY + XmY)
            # Within a TDA calculation XpY and XmY are the same and
            # Y will only contain zeros.
            Y = XpY - X
        # In the unrestricted case we can't use 2 for the 3. dimension anymore
        # (same number of alpha and beta electrons). The sizes for the respective
        # spins would be (a_act*a_vir) and (b_act*b_vir).
        else:
            raise Exception("Unrestricted not supported yet!")
        # Shapes of X and Y
        # 1. dim: roots (more than initially requested) -> variable
        # 2. dim: alpha, beta                           -> 2
        # 3. dim: act. x virt.                          -> variable

        # xpy_fn = f"{dump_fn}_XpY"
        # np.save(xpy_fn, XpY)
        # xmy_fn = f"{dump_fn}_XmY"
        # np.save(xmy_fn, XmY)
        # x_fn = f"{dump_fn}_X"
        # np.save(x_fn, X)

        # As we only handle restricted calculations for now the alpha
        # and beta part are identical. So we only keep the alpha part.
        if self.nmos.restricted:
            # Drop beta part and restrict to the requested states
            X = X[: self.nstates, 0, :]
            Y = Y[: self.nstates, 0, :]

        X_full = np.zeros((self.nstates, nmos.a_occ, nmos.a_vir))
        X_full[:, nmos.a_core :] = X.reshape(-1, nmos.a_act, nmos.a_vir)
        Y_full = np.zeros((self.nstates, nmos.a_occ, nmos.a_vir))
        Y_full[:, nmos.a_core :] = Y.reshape(-1, nmos.a_act, nmos.a_vir)

        return X_full, Y_full

    def prepare_overlap_data(self, path):
        # Parse X eigenvector from 635r dump
        X, Y = self.parse_635r_dump(self.dump_635r, self.roots, self.nmos)

        # From http://gaussian.com/cis/, Examples tab, Normalization:
        #
        # 'For closed shell calculations, the sum of the squares of the
        # expansion coefficients is normalized to total 1/2 (as the beta
        # coefficients are not shown).'
        #
        # Right now we only deal with restricted calculatios, so alpha == beta
        # and we ignore beta. So we are lacking a factor of sqrt(2). Another
        # option would be to normalize all states to 1.
        # ci_coeffs *= 2**0.5

        # Parse mo coefficients from .fchk file and write a 'fake' turbomole
        # mos file.
        keys = (
            "SCF Energy",
            "Alpha Orbital Energies",
            "Alpha MO coefficients",
            "ETran state values",
        )
        scf_key, mo_energies_key, mo_key, exc_key = keys
        fchk_dict = self.parse_fchk(self.fchk, keys=keys)
        mo_coeffs = fchk_dict[mo_key]
        mo_energies = fchk_dict[mo_energies_key]
        mo_coeffs = mo_coeffs.reshape(-1, mo_energies.size)

        gs_energy = fchk_dict[scf_key]
        exc_data = fchk_dict[exc_key].reshape(-1, 16)
        exc_energies = exc_data[:, 0]
        all_energies = np.zeros(len(exc_energies) + 1)
        all_energies[0] = gs_energy
        all_energies[1:] += exc_energies
        return mo_coeffs, X, Y, all_energies

    def parse_force(self, path):
        results = {}
        keys = ("SCF Energy", "Total Energy", "Cartesian Gradient")
        fchk_path = Path(path) / f"{self.fn_base}.fchk"
        fchk_dict = self.parse_fchk(fchk_path, keys)
        results["energy"] = fchk_dict["SCF Energy"]
        results["forces"] = -fchk_dict["Cartesian Gradient"]

        if self.nstates:
            # This sets the proper excited state energy in the
            # results dict and also stores all energies.
            exc_energies = self.parse_tddft(path)
            # G16 root input is 1 based, so we substract 1 to get
            # the right index here.
            root_exc_en = exc_energies[self.root - 1]
            gs_energy = fchk_dict["SCF Energy"]
            # Add excitation energy to ground state energy.
            results["energy"] += root_exc_en
            # Create a new array including the ground state energy
            # to save the energies of all states.
            all_ens = np.full(len(exc_energies) + 1, gs_energy)
            all_ens[1:] += exc_energies
            results["all_energies"] = all_ens

        return results

    def parse_hessian(self, path):
        keys = (
            "Total Energy",
            "Cartesian Gradient",
            "Cartesian Force Constants",
        )
        fchk_path = Path(path) / f"{self.fn_base}.fchk"
        fchk_dict = self.parse_fchk(fchk_path, keys=keys)
        en_key, grad_key, hess_key = keys
        grad = fchk_dict[grad_key]
        tril_hess = fchk_dict[hess_key]
        tril_inds = np.tril_indices(grad.size)
        full_hessian = np.zeros((grad.size, grad.size))
        full_hessian[tril_inds] = tril_hess
        triu_inds = np.triu_indices(grad.size, k=1)
        full_hessian[triu_inds] = full_hessian.T[triu_inds]
        results = {
            "energy": fchk_dict[en_key],
            "forces": -grad,
            "hessian": full_hessian,
        }
        return results

    def parse_double_mol(self, path, out_fn=None):
        def repl_double(s, loc, toks):
            return toks[0].replace("D", "E")

        if out_fn is None:
            out_fn = self.out_fn

        with open(path / out_fn) as handle:
            text = handle.read()
        # Number of basis functions in the double molecule
        nbas = int(re.search(r"NBasis =\s*(\d+)", text)[1])
        assert nbas % 2 == 0
        # Gaussian prints columns of a triangular matrix including
        # the diagonal
        backup_white = pp.ParserElement.DEFAULT_WHITE_CHARS
        pp.ParserElement.setDefaultWhitespaceChars(" \t")
        int_ = pp.Suppress(pp.Word(pp.nums))

        float_ = pp.Word(pp.nums + ".D+-").setParseAction(repl_double)
        nl = pp.Suppress(pp.Literal("\n"))
        header = pp.OneOrMore(int_) + nl
        line = int_ + pp.OneOrMore(float_) + nl

        block = pp.Group(header + pp.OneOrMore(~header + line))

        parser = pp.Suppress(
            pp.SkipTo("*** Overlap *** \n", include=True)
        ) + pp.OneOrMore(block).setResultsName("blocks")

        result = parser.parseFile(path / out_fn)
        pp.ParserElement.setDefaultWhitespaceChars(backup_white)

        # The full double molecule overlap matrix (square)
        full_mat = np.zeros((nbas, nbas))

        def get_block_inds(block_ind, nbas, cols=5):
            """Returns the indices as required to assign the matrix
            printed by Gaussian. Gaussian prints a lower triangle
            matrix in blocks of five columns each."""
            start_row = block_ind * cols
            start_col = start_row
            inds = list()
            for row in range(start_row, nbas):
                col = start_col
                while (col <= row) and (col < start_row + cols):
                    inds.append((row, col))
                    col += 1
            return inds

        for i, block in enumerate(result["blocks"]):
            rows, cols = zip(*get_block_inds(i, nbas))
            full_mat[rows, cols] = block.asList()

        nbas_single = nbas // 2
        """The whole matrix consists of four blocks:
            Original overlaps of molecule 1
                b1 = full_mat[:nbas_single, :nbas_single]
            Zero, as we only get the lower triangle
                b2 = full_mat[:nbas_single, nbas_single:]
            Double molecule overlaps between basis functions of
            molecule 1 and molecule 2
                b3 = full_mat[nbas_single:, :nbas_single]
            Original overlaps of molecule 2
                b4 = full_mat[nbas_single:, nbas_single:]
        """
        double_mol_ovlp = full_mat[nbas_single:, :nbas_single]
        return double_mol_ovlp

    def parse_charges(self, path=None):
        if path is None and self.fchk is not None:
            fchk_path = self.fchk
        elif path is not None and path.endswith(".fchk"):
            fchk_path = path
        elif path is not None:
            fchk_path = Path(path) / f"{self.fn_base}.fchk"

        keys = ("Mulliken Charges",)
        fchk_dict = self.parse_fchk(fchk_path, keys=keys)
        charges = np.array(fchk_dict["Mulliken Charges"])

        return charges

    def keep(self, path):
        kept_fns = super().keep(path)
        if self.keep_chk:
            self.chk = kept_fns[".chk"]
        try:
            self.fchk = kept_fns["fchk"]
            self.mwfn_wf = self.fchk
        except KeyError:
            self.log("No .fchk file found!")
            return
        if self.track:
            self.dump_635r = kept_fns["dump_635r"]

    def get_chkfiles(self):
        return {
            "fchk": self.fchk,
        }

    def set_chkfiles(self, chkfiles):
        try:
            fchk = chkfiles["fchk"]
            self.fchk = fchk
            self.log(f"Set chkfile '{fchk}' as fchk.")
        except KeyError:
            self.log("Found no fchk information in chkfiles!")

    def __str__(self):
        return "Gaussian16 calculator"
