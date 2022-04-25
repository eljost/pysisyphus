import glob
import os
import re
import struct
import shutil
import subprocess

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.OverlapCalculator import OverlapCalculator
from pysisyphus.constants import BOHR2ANG, ANG2BOHR
from pysisyphus.helpers_pure import file_or_str


def make_sym_mat(table_block):
    mat_size = int(table_block[1])
    # Orca prints blocks of 5 columns
    arr = np.array(table_block[2:], dtype=float)
    assert arr.size == mat_size ** 2
    block_size = 5 * mat_size
    cbs = [
        arr[i * block_size : (i + 1) * block_size].reshape(mat_size, -1)
        for i in range(arr.size // block_size + 1)
    ]
    return np.concatenate(cbs, axis=1)


def save_orca_pc_file(point_charges, pc_fn, hardness=None):
    point_charges = point_charges.copy()
    # ORCA excepcts point charge positions in Angstrom
    point_charges[:, :3] *= BOHR2ANG

    # ORCA also expects the ordering <q> <x> <y> <z>, so we have to resort.
    shape = point_charges.shape
    if hardness is not None:
        shape = shape[0], shape[1] + 1
    point_charges_orca = np.zeros_like(point_charges)
    point_charges_orca = np.zeros(shape)
    point_charges_orca[:, 0] = point_charges[:, 3]
    point_charges_orca[:, 1:4] = point_charges[:, :3]

    if hardness:
        point_charges_orca[:, 4] = hardness
    np.savetxt(
        pc_fn,
        point_charges_orca,
        fmt="%16.10f",
        header=str(len(point_charges)),
        comments="",
    )


class ORCA(OverlapCalculator):

    conf_key = "orca"

    def __init__(
        self,
        keywords,
        blocks="",
        gbw=None,
        do_stable=False,
        numfreq=False,
        **kwargs,
    ):
        """ORCA calculator.

        Wrapper for creating ORCA input files for energy, gradient
        and Hessian calculations. The PAL and memory inputs must not
        be given in the keywords and/or blocks, as they are handled
        by the 'pal' and 'memory' arguments.

        Parameters
        ----------
        keywords : str
            Keyword line, as normally given in ORCA, excluding the
            leading "!".
        blocks : str, optional
            ORCA block input(s), e.g. for TD-DFT calculations (%tddft ... end).
            As the blocks start with a leading "%", wrapping the input in quotes
            ("") is required, otherwise the parsing will fail.
        gbw : str, optional
            Path to an input gbw file, which will be used as initial guess
            for the first calculation. Will be overriden later, with the
            path to the gbw file of a previous calculation.
        do_stable: bool, optional
            Run stability analysis until a stable wavefunction is obtained,
            before every calculation.
        numfreq : boo, optional
            Use numerical frequencies instead of analytical ones.
        mem : int
            Mememory per core in MB.
        """
        super().__init__(**kwargs)

        self.keywords = keywords.lower()
        self.blocks = blocks.lower()
        self.gbw = gbw
        self.do_stable = bool(do_stable)
        self.freq_keyword = "numfreq" if numfreq else "freq"

        assert ("pal" not in keywords) and ("nprocs" not in blocks), (
            "PALn/nprocs not " "allowed! Use 'pal: n' in the 'calc' section instead."
        )
        assert "maxcore" not in blocks, (
            "maxcore not allowed! " "Use 'mem: n' in the 'calc' section instead!"
        )

        self.to_keep = (
            "inp",
            "out:orca.out",
            "gbw",
            "engrad",
            "hessian",
            "cis",
            "molden:orca.molden",
            "hess",
            "pcgrad",
        )
        self.do_tddft = False
        if "tddft" in self.blocks:
            self.do_tddft = True
            try:
                self.root = int(re.search(r"iroot\s*(\d+)", self.blocks).group(1))
            except AttributeError:
                self.log("Doing TDA/TDDFT calculation without gradient.")
        self.triplets = bool(re.search(r"triplets\s+true", self.blocks))
        self.inp_fn = "orca.inp"
        self.out_fn = "orca.out"

        self.orca_input = """!{keywords} {calc_type}
        {moinp}

        %pal nprocs {pal} end
        %maxcore {mem}

        {blocks}
        {pointcharges}

        *xyz {charge} {mult}
        {coords}
        *
        """

        self.parser_funcs = {
            "energy": self.parse_energy,
            "grad": self.parse_engrad,
            "hessian": self.parse_hessian,
            "noparse": lambda path: None,
            "stable": self.parse_stable,
        }

        self.base_cmd = self.get_cmd()

    def reattach(self, last_calc_cycle):
        # Use the latest .gbw
        gbw = self.make_fn("gbw", last_calc_cycle)
        self.log(f"Restarted. using {gbw}")

    def get_moinp_str(self, gbw):
        moinp_str = ""
        if gbw:
            moinp_str = f"""!moread
            %moinp "{os.path.abspath(gbw)}" """
        return moinp_str

    def prepare_input(
        self, atoms, coords, calc_type, point_charges=None, do_stable=False
    ):
        coords = self.prepare_coords(atoms, coords)
        if self.gbw:
            self.log(f"Using {self.gbw}")
        else:
            self.log("Using initial guess provided by ORCA")
        if calc_type == "noparse":
            calc_type = ""

        pc_str = ""
        if point_charges is not None:
            pc_fn = self.make_fn("pointcharges_inp.pc")
            save_orca_pc_file(point_charges, pc_fn)
            pc_str = f'%pointcharges "{pc_fn}"'
        stable_block = "\n%scf stabperform true hftyp uhf end" if do_stable else ""
        blocks = self.get_block_str() + stable_block

        inp = self.orca_input.format(
            keywords=self.keywords,
            calc_type=calc_type,
            moinp=self.get_moinp_str(self.gbw),
            pal=self.pal,
            mem=self.mem,
            blocks=blocks,
            pointcharges=pc_str,
            coords=coords,
            charge=self.charge,
            mult=self.mult,
        )
        return inp

    def get_block_str(self):
        block_str = self.blocks
        # Use the correct root if we track it
        if self.track:
            block_str = re.sub(r"iroot\s+(\d+)", f"iroot {self.root}", self.blocks)
            self.log(f"Using iroot '{self.root}' for excited state gradient.")
        return block_str

    def get_stable_wavefunction(self, atoms, coords):
        self.log("Trying to get a stable wavefunction")

        stable = False
        max_cycles = 10
        for i in range(max_cycles):
            inp = self.prepare_input(atoms, coords, calc_type="", do_stable=True)
            stable = self.run(inp, calc="stable")
            self.log(f"{i:02d} stable: {stable}")

            if stable:
                self.log(f"Found stable wavefunction in cycle {i}!")
                break
        else:
            raise Exception(
                "Could not find stable wavefunction in {max_cycles}! " "Aborting."
            )

    def parse_stable(self, path):
        with open(path / "orca.out") as handle:
            text = handle.read()

        stable_re = re.compile("Stability Analysis indicates a stable")
        stable = bool(stable_re.search(text))

        unstable_re = re.compile("Stability Analysis indicates an UNSTABLE")
        unstable = bool(unstable_re.search(text))

        stable = stable and not unstable

        return stable

    def store_and_track(self, results, func, atoms, coords, **prepare_kwargs):
        if self.track:
            self.store_overlap_data(atoms, coords)
            if self.track_root():
                # Redo the calculation with the updated root
                results = func(atoms, coords, **prepare_kwargs)
        return results

    def get_energy(self, atoms, coords, **prepare_kwargs):
        calc_type = ""

        if self.do_stable:
            self.get_stable_wavefunction(atoms, coords)

        inp = self.prepare_input(atoms, coords, calc_type, **prepare_kwargs)
        results = self.run(inp, calc="energy")
        results = self.store_and_track(
            results, self.get_energy, atoms, coords, **prepare_kwargs
        )
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        if self.do_stable:
            self.get_stable_wavefunction(atoms, coords)

        calc_type = "engrad"
        inp = self.prepare_input(atoms, coords, calc_type, **prepare_kwargs)
        kwargs = {
            "calc": "grad",
        }
        results = self.run(inp, **kwargs)
        results = self.store_and_track(
            results, self.get_forces, atoms, coords, **prepare_kwargs
        )
        return results

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        calc_type = self.freq_keyword

        if self.do_stable:
            self.get_stable_wavefunction(atoms, coords)

        inp = self.prepare_input(atoms, coords, calc_type, **prepare_kwargs)
        results = self.run(inp, calc="hessian")
        # results = self.store_and_track(
        # results, self.get_hessian, atoms, coords, **prepare_kwargs
        # )
        return results

    def run_calculation(self, atoms, coords, **prepare_kwargs):
        """Basically some kind of dummy method that can be called
        to execute ORCA with the stored cmd of this calculator."""
        inp = self.prepare_input(atoms, coords, "noparse", **prepare_kwargs)
        kwargs = {
            "calc": "energy",
        }
        results = self.run(inp, **kwargs)
        if self.track:
            self.store_overlap_data(atoms, coords)
        return results

    def run_after(self, path):
        # Create .molden file when CDDs are requested
        if self.cdds:
            cmd = "orca_2mkl orca -molden".split()
            proc = subprocess.Popen(
                cmd,
                cwd=path,
                universal_newlines=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            proc.wait()
            shutil.copy(path / "orca.molden.input", path / "orca.molden")

    @staticmethod
    @file_or_str(".hess", method=False)
    def parse_hess_file(text):
        integer = pp.Word(pp.nums)
        float_ = pp.Word(pp.nums + ".-")
        plus = pp.Literal("+")
        minus = pp.Literal("-")
        E = pp.Literal("E")
        scientific = pp.Combine(float_ + E + pp.Or([plus, minus]) + integer)

        table_header_line = pp.Suppress(integer + pp.restOfLine)
        scientific_line = pp.Suppress(integer) + pp.OneOrMore(scientific)
        scientific_block = table_header_line + pp.OneOrMore(scientific_line)
        float_line = pp.Suppress(integer) + float_
        comment_line = pp.Literal("#") + pp.restOfLine
        mass_xyz_line = pp.Group(
            pp.Word(pp.alphas) + float_ + pp.Group(pp.OneOrMore(float_))
        )

        block_name = pp.Word(pp.alphas + "$_")
        block_length = integer

        block_int = block_name + block_length
        block_float = block_name + float_
        block_table = block_name + integer + pp.OneOrMore(scientific_block)
        block_table_two_int = (
            block_name + integer + pp.Suppress(integer) + pp.OneOrMore(scientific_block)
        )
        block_float_table = block_name + integer + pp.OneOrMore(float_line)
        block_atoms = block_name + integer + pp.OneOrMore(mass_xyz_line)

        act_atom = block_int.setResultsName("act_atom")
        act_coord = block_int.setResultsName("act_coord")
        act_energy = block_float.setResultsName("act_energy")
        hessian = block_table.setResultsName("hessian")
        vib_freqs = block_float_table.setResultsName("vib_freqs")
        normal_modes = block_table_two_int.setResultsName("normal_modes")
        atoms = block_atoms.setResultsName("atoms")

        parser = (
            block_name
            + act_atom
            + act_coord
            + act_energy
            + hessian
            + vib_freqs
            + normal_modes
            + pp.OneOrMore(comment_line)
            + atoms
        )
        parsed = parser.parseString(text)
        return parsed

    def parse_hessian(self, path):
        hessian_fn = glob.glob(os.path.join(path, "*.hess"))
        assert len(hessian_fn) == 1
        hessian_fn = hessian_fn[0]
        if not hessian_fn:
            raise Exception("ORCA calculation failed.")

        parsed = ORCA.parse_hess_file(hessian_fn)

        # logging.warning("Hacky orca energy parsing in orca hessian calculation!")
        orca_log_fn = os.path.join(path, self.out_fn)
        with open(orca_log_fn) as handle:
            log_text = handle.read()

        energy_re = r"FINAL SINGLE POINT ENERGY\s*([-\.\d]+)"
        energy_mobj = re.search(energy_re, log_text)
        energy = float(energy_mobj.groups()[0])

        results = {
            "energy": energy,
            "hessian": make_sym_mat(parsed["hessian"]),
        }

        return results

    def parse_energy(self, path):
        log_fn = glob.glob(os.path.join(path / "orca.out"))
        if not log_fn:
            raise Exception("ORCA calculation failed.")

        assert len(log_fn) == 1
        log_fn = log_fn[0]
        with open(log_fn) as handle:
            text = handle.read()
        mobj = re.search(r"FINAL SINGLE POINT ENERGY\s+([\d\-\.]+)", text)
        energy = float(mobj[1])
        return {"energy": energy}

    def parse_engrad(self, path):
        results = {}
        engrad_fn = glob.glob(os.path.join(path, "*.engrad"))
        if not engrad_fn:
            raise Exception("ORCA calculation failed.")

        assert len(engrad_fn) == 1
        engrad_fn = engrad_fn[0]
        with open(engrad_fn) as handle:
            engrad = handle.read()
        engrad = re.findall(r"([\d\-\.]+)", engrad)
        atoms = int(engrad.pop(0))
        energy = float(engrad.pop(0))
        force = -np.array(engrad[: 3 * atoms], dtype=float)
        results["energy"] = energy
        results["forces"] = force

        return results

    def parse_cis(self, cis):
        """
        Read binary CI vector file from ORCA.
            Adapted from TheoDORE 1.7.1, Authors: S. Mai, F. Plasser
            https://sourceforge.net/p/theodore-qc
        """
        cis_handle = open(cis, "rb")
        self.log(f"Parsing CI vectors from {cis_handle}")

        # the header consists of 9 4-byte integers, the first 5
        # of which give useful info.
        nvec = struct.unpack("i", cis_handle.read(4))[0]
        # header array contains:
        # [0] index of first alpha occ,  is equal to number of frozen alphas
        # [1] index of last  alpha occ
        # [2] index of first alpha virt
        # [3] index of last  alpha virt, header[3]+1 is equal to number of bfs
        # [4] index of first beta  occ,  for restricted equal to -1
        # [5] index of last  beta  occ,  for restricted equal to -1
        # [6] index of first beta  virt, for restricted equal to -1
        # [7] index of last  beta  virt, for restricted equal to -1
        header = [struct.unpack("i", cis_handle.read(4))[0] for i in range(8)]

        # Assert that all flags regarding unrestricted calculations are -1
        if any([flag != -1 for flag in header[4:8]]):
            raise Exception("parse_cis, no support for unrestricted MOs")

        nfrzc = header[0]
        nocc = header[1] + 1
        nact = nocc - nfrzc
        nmo = header[3] + 1
        nvir = nmo - header[2]
        lenci = nact * nvir
        self.log(f"nmo = {nmo}, nocc = {nocc}, nact = {nact}, nvir = {nvir}")

        # Loop over states. For non-TDA order is: X+Y of 1, X-Y of 1,
        # X+Y of 2, X-Y of 2, ...
        prev_root = -1
        prev_mult = 1
        iroot_triplets = 0

        # Flags that may later be set to True
        triplets = False
        tda = False
        Xs = list()
        Ys = list()

        for ivec in range(nvec):
            # header of each vector
            # contains 6 4-byte ints, then 1 8-byte double, then 8 byte unknown
            nele, d1, mult, d2, iroot, d3 = struct.unpack("iiiiii", cis_handle.read(24))

            # Will evaluate True only once when triplets were requested.
            if prev_mult != mult:
                triplets = True
                prev_root = -1

            # When we encounter the second "state" we can decide if it is a TDA
            # calculation (without Y-vector).
            if (ivec == 1) and (iroot == prev_root + 1):
                tda = True

            if triplets:
                iroot = iroot_triplets

            ene, d3 = struct.unpack("dd", cis_handle.read(16))
            self.log(f"ivec={ivec}, nele={nele}, mult={mult}, iroot={iroot}")
            # Then come nact * nvirt 8-byte doubles with the coefficients
            coeffs = struct.unpack(lenci * "d", cis_handle.read(lenci * 8))
            coeffs = np.array(coeffs).reshape(-1, nvir)
            # create full array, i.e nocc x nvirt
            coeffs_full = np.zeros((nocc, nvir))
            coeffs_full[nfrzc:] = coeffs

            # In this case, we have a non-TDA state, where Y is present!
            # We can recover the original X and Y by first computing X as
            #   X = (X+Y + X-Y) / 2
            # and then
            #   Y = X+Y - X
            if prev_root == iroot:
                X_plus_Y = Xs[-1]
                X_minus_Y = coeffs_full
                X = 0.5 * (X_plus_Y + X_minus_Y)
                Y = X_plus_Y - X
                Xs[-1] = X
                Ys[-1] = Y
            else:
                Xs.append(coeffs_full)
                Ys.append(np.zeros_like(coeffs_full))

            # Somehow ORCA stops to update iroot correctly after the singlet states.
            if (mult == 3) and (tda or (ivec % 2) == 1):
                iroot_triplets += 1

            prev_root = iroot
            prev_mult = mult
        cis_handle.close()

        Xs = np.array(Xs)
        Ys = np.array(Ys)

        # Only return triplet states if present
        if triplets:
            assert (len(Xs) % 2) == 0
            states = len(Xs) // 2
            Xs = Xs[states:]
            Ys = Ys[states:]
        return Xs, Ys

    def parse_gbw(self, gbw_fn):
        """Adapted from
        https://orcaforum.kofo.mpg.de/viewtopic.php?f=8&t=3299

        The first 5 long int values represent pointers into the file:

        Pointer @+0:  Internal ORCA data structures
        Pointer @+8:  Geometry
        Pointer @+16: BasisSet
        Pointer @+24: Orbitals
        Pointer @+32: ECP data
        """

        with open(gbw_fn, "rb") as handle:
            handle.seek(24)
            offset = struct.unpack("<q", handle.read(8))[0]
            handle.seek(offset)
            operators = struct.unpack("<i", handle.read(4))[0]
            dimension = struct.unpack("<i", handle.read(4))[0]

            # print('Offset: {}'.format(offset))
            # print('Number of Operators: {}'.format(operators))
            # print('Basis Dimension: {}'.format(dimension))

            coeffs_fmt = "<" + dimension ** 2 * "d"

            assert operators == 1, "Unrestricted case is not implemented!"

            for i in range(operators):
                # print('\nOperator: {}'.format(i))
                coeffs = struct.unpack(coeffs_fmt, handle.read(8 * dimension ** 2))
                occupations = struct.iter_unpack("<d", handle.read(8 * dimension))
                energies = struct.iter_unpack("<d", handle.read(8 * dimension))
                irreps = struct.iter_unpack("<i", handle.read(4 * dimension))
                cores = struct.iter_unpack("<i", handle.read(4 * dimension))

                coeffs = np.array(coeffs).reshape(-1, dimension).T
                energies = np.array([en[0] for en in energies])

                # print('Coefficients')
                # for coef in coefficients:
                # print('{:16.12f}'.format(*coef))
                # print('Occupations')
                # for occupation in occupations:
                # print('{:16.12f}'.format(*occupation))
                # print('Energies')
                # for energy in energies:
                # print('{:16.12f}'.format(*energy))
                # print('Irreps')
                # for irrep in irreps:
                # print('{}'.format(*irrep))
                # print('Core')
                # for core in cores:
                # print('{}'.format(*core))
            return coeffs, energies

    @staticmethod
    def set_mo_coeffs_in_gbw(in_gbw_fn, out_gbw_fn, mo_coeffs):
        """See self.parse_gbw."""

        with open(in_gbw_fn, "rb") as handle:
            handle.seek(24)
            offset = struct.unpack("<q", handle.read(8))[0]
            handle.seek(offset)
            operators = struct.unpack("<i", handle.read(4))[0]
            dimension = struct.unpack("<i", handle.read(4))[0]
            assert operators == 1, "Unrestricted case is not implemented!"

            handle.seek(0)
            gbw_bytes = handle.read()

        tot_offset = offset + 4 + 4
        start = gbw_bytes[:tot_offset]
        end = gbw_bytes[tot_offset + 8 * dimension ** 2 :]
        # Construct new gbw content by replacing the MO coefficients in the middle
        mod_gbw_bytes = start + mo_coeffs.T.tobytes() + end

        with open(out_gbw_fn, "wb") as handle:
            handle.write(mod_gbw_bytes)

    def parse_all_energies(self, text=None, triplets=None):
        if triplets is None:
            triplets = self.triplets

        if text is None:
            with open(self.out) as handle:
                text = handle.read()

        energy_re = r"FINAL SINGLE POINT ENERGY\s*([-\.\d]+)"
        energy_mobj = re.search(energy_re, text)
        gs_energy = float(energy_mobj.groups()[0])
        all_energies = [gs_energy]

        if self.do_tddft:
            scf_re = re.compile(r"E\(SCF\)\s+=\s*([\d\-\.]+) Eh")
            scf_mobj = scf_re.search(text)
            scf_en = float(scf_mobj.group(1))
            gs_energy = scf_en
            tddft_re = re.compile(r"STATE\s*(\d+):\s*E=\s*([\d\.]+)\s*au")
            states, exc_ens = zip(
                *[(int(state), float(en)) for state, en in tddft_re.findall(text)]
            )
            if triplets:
                roots = len(states) // 2
                exc_ens = exc_ens[-roots:]
                states = states[-roots:]
            assert len(exc_ens) == len(set(states))
            all_energies = np.full(1 + len(exc_ens), gs_energy)
            all_energies[1:] += exc_ens
        all_energies = np.array(all_energies)
        return all_energies

    @staticmethod
    @file_or_str(".out", method=False)
    def parse_atoms_coords(text):
        ac_re = re.compile(
            r"CARTESIAN COORDINATES \(ANGSTROEM\)\s+\-{33}(.+?)\s+\-{28}", re.DOTALL
        )
        mobj = ac_re.search(text)
        atoms_coords = mobj.group(1).strip().split()
        # atoms, *coords = np.array(atoms_coords).reshape(-1, 4).T
        atoms_coords = np.array(atoms_coords).reshape(-1, 4)
        atoms = tuple(atoms_coords[:, 0])
        coords = atoms_coords[:, 1:].astype(float).flatten() * ANG2BOHR
        return atoms, coords

    @staticmethod
    @file_or_str(".out", method=False)
    def parse_engrad_info(text):
        soi_re = re.compile(r"State of interest\s+\.{3}\s+(\d+)")
        try:
            root = soi_re.search(text).group(1)
            root = int(root)
        except AttributeError:
            root = None
        triplets = bool(re.search(r"triplets\s+true", text))
        return root, triplets

    def parse_mo_numbers(self, out_fn):
        with open(out_fn) as handle:
            text = handle.read()
        electron_re = r"NEL\s*....\s*(\d+)"
        electrons = int(re.search(electron_re, text)[1])
        assert electrons % 2 == 0, "unrestricted is not yet supported!"
        occ_num = int(electrons / 2)

        mo_re = r"Dim\s*....\s*(\d+)"
        mo_num = int(re.search(mo_re, text)[1])
        virt_num = mo_num - occ_num
        self.log(
            f"found {electrons} electrons, {mo_num} MOs, with "
            f"{occ_num} occupied and {virt_num} virtual."
        )
        return occ_num, virt_num

    def set_mo_coeffs(self, mo_coeffs=None, gbw=None):
        if mo_coeffs is not None:
            self.mo_coeffs = mo_coeffs
            return
        if not gbw and self.gbw:
            gbw = self.gbw
        else:
            raise Exception("Got no .gbw file to parse!")
        self.log(f"Setting MO coefficients from {gbw}.")
        self.mo_coeffs, _ = self.parse_gbw(self.gbw)

    def prepare_overlap_data(self, path):
        # Parse eigenvectors from tda/tddft calculation
        X, Y = self.parse_cis(self.cis)
        # Parse mo coefficients from gbw file and write a 'fake' turbomole
        # mos file.
        mo_coeffs, _ = self.parse_gbw(self.gbw)
        all_energies = self.parse_all_energies()
        return mo_coeffs, X, Y, all_energies

    def keep(self, path):
        kept_fns = super().keep(path)
        self.gbw = kept_fns["gbw"]
        self.out = kept_fns["out"]
        if self.do_tddft:
            self.cis = kept_fns["cis"]
        try:
            self.mwfn_wf = kept_fns["molden"]
        except KeyError:
            self.log("Didn't set 'mwfn_wf'. No .molden file in kept_fns.")

    def get_chkfiles(self):
        return {
            "gbw": self.gbw,
        }

    def set_chkfiles(self, chkfiles):
        try:
            gbw = chkfiles["gbw"]
            self.gbw = gbw
            self.log(f"Set chkfile '{gbw}' as gbw.")
        except KeyError:
            self.log("Found no gbw information in chkfiles!")

    @file_or_str(".out", method=True)
    def check_termination(self, text):
        term_re = re.compile(r"\*{4}ORCA TERMINATED NORMALLY\*{4}")
        mobj = term_re.search(text)
        return bool(mobj)

    def clean_tmp(self, path):
        tmp_fns = path.glob("*.tmp")
        for tmp in tmp_fns:
            os.remove(tmp)
            self.log(f"Removed '{tmp}'")
        # try:
        # os.remove(path / "orca.gbw")
        # except FileNotFoundError:
        # pass

    def __str__(self):
        return f"ORCA({self.name})"
