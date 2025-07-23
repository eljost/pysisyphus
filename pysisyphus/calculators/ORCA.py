import glob
import io
from math import sqrt
from pathlib import Path
import os
import re
import struct
import shutil
import warnings

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.MOCoeffs import MOCoeffs
from pysisyphus.calculators.OverlapCalculator import OverlapCalculator
from pysisyphus.constants import BOHR2ANG, ANG2BOHR
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.wavefunction import norm_ci_coeffs, Wavefunction


def make_sym_mat(table_block):
    mat_size = int(table_block[1])
    # Orca prints blocks of 5 columns
    arr = np.array(table_block[2:], dtype=float)
    assert arr.size == mat_size**2
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


def parse_orca_gbw(gbw_fn):
    """Adapted from
    https://orcaforum.kofo.mpg.de/viewtopic.php?f=8&t=3299

    The first 5 long int values represent pointers into the file:

    Pointer @+0:  Internal ORCA data structures
    Pointer @+8:  Geometry
    Pointer @+16: BasisSet
    Pointer @+24: Orbitals
    Pointer @+32: ECP data
    """

    warnings.warn(
        "'parse_orca_gbw()' is deprecated. Please use 'parse_orca_gbw_new()'.",
        warnings.DeprecationWarning,
    )

    with open(gbw_fn, "rb") as handle:
        handle.seek(24)
        offset = struct.unpack("<q", handle.read(8))[0]
        handle.seek(offset)
        operators = struct.unpack("<i", handle.read(4))[0]
        dimension = struct.unpack("<i", handle.read(4))[0]

        coeffs_fmt = "<" + dimension**2 * "d"
        assert operators == 1, "Unrestricted case is not implemented!"

        for i in range(operators):
            # print('\nOperator: {}'.format(i))
            coeffs = struct.unpack(coeffs_fmt, handle.read(8 * dimension**2))
            occupations = struct.iter_unpack("<d", handle.read(8 * dimension))
            energies = struct.iter_unpack("<d", handle.read(8 * dimension))
            irreps = struct.iter_unpack("<i", handle.read(4 * dimension))
            cores = struct.iter_unpack("<i", handle.read(4 * dimension))

            coeffs = np.array(coeffs).reshape(-1, dimension)
            energies = np.array([en[0] for en in energies])

        # MOs are returned in columns
        return coeffs, energies


def parse_orca_gbw_new(gbw_fn: str) -> MOCoeffs:
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

        coeffs_fmt = "<" + dimension**2 * "d"
        assert operators in (1, 2)

        keys = {
            0: "a",
            1: "b",
        }
        kwargs = {}
        for i in range(operators):
            key = keys[i]
            coeffs = struct.unpack(coeffs_fmt, handle.read(8 * dimension**2))
            occupations = struct.iter_unpack("<d", handle.read(8 * dimension))
            energies = struct.iter_unpack("<d", handle.read(8 * dimension))
            irreps = struct.iter_unpack("<i", handle.read(4 * dimension))
            cores = struct.iter_unpack("<i", handle.read(4 * dimension))

            coeffs = np.array(coeffs).reshape(-1, dimension)
            occupations = np.array([occ[0] for occ in occupations])
            energies = np.array([en[0] for en in energies])
            kwargs[f"C{key}"] = coeffs
            kwargs[f"ens{key}"] = energies
            kwargs[f"occs{key}"] = occupations

        mo_coeffs = MOCoeffs(**kwargs)
        return mo_coeffs


def set_mo_coeffs_in_gbw(
    gbw_in,
    gbw_out,
    alpha_mo_coeffs=None,
    beta_mo_coeffs=None,
    alpha_energies=None,
    alpha_occs=None,
    beta_energies=None,
    beta_occs=None,
):
    """MOs are expected to be in columns."""

    with open(gbw_in, "rb") as handle:
        handle.seek(24)
        offset = struct.unpack("<q", handle.read(8))[0]
        handle.seek(offset)
        operators = struct.unpack("<i", handle.read(4))[0]
        if beta_mo_coeffs is not None:
            assert operators == 2
        # Number of MOs
        dimension = struct.unpack("<i", handle.read(4))[0]

        handle.seek(0)
        gbw_bytes = handle.read()

    # offset + operators + dimension
    tot_offset = offset + 4 + 4
    start = gbw_bytes[:tot_offset]

    # Update alpha MO coefficients
    _4dim = 4 * dimension
    _8dim = 2 * _4dim
    _8dim2 = _8dim * dimension
    coeffs_size = _8dim2
    occs_size = energies_size = _8dim
    irreps_size = cores_size = _4dim
    block_size = coeffs_size + occs_size + energies_size + irreps_size + cores_size

    occ_start = coeffs_size
    ens_start = occ_start + occs_size
    ens_end = ens_start + energies_size

    def update_block(start, mo_coeffs, occs, energies):
        block = gbw_bytes[start : start + block_size]

        # Read existing data from block or convert provided arguments to bytes.

        # Occupation numbers
        if occs is None:
            occs = block[occ_start:ens_start]
        else:
            occs = np.array(occs, dtype=float).tobytes()
        # MO energies
        if energies is None:
            energies = block[ens_start:ens_end]
        else:
            energies = np.array(energies, dtype=float).tobytes()
        # MO coefficients
        if mo_coeffs is None:
            mo_coeffs = block[:coeffs_size]
        else:
            mo_coeffs = np.array(mo_coeffs, dtype=float).tobytes()

        # Reassemble block
        block = mo_coeffs + occs + energies + block[coeffs_size + _8dim + _8dim :]
        return block

    # The alpha block is always present in the .gbw ...
    alpha_block = gbw_bytes[tot_offset : tot_offset + block_size]
    alpha_block = update_block(tot_offset, alpha_mo_coeffs, alpha_occs, alpha_energies)

    # The beta block is only present for operators == 2.
    if operators == 2:
        beta_block = update_block(
            tot_offset + block_size, beta_mo_coeffs, beta_occs, beta_energies
        )
    else:
        beta_block = bytes()

    # Reassemble modified .gbw file
    mod_gbw_bytes = start + alpha_block + beta_block
    # and write the bytes.
    with open(gbw_out, "wb") as handle:
        handle.write(mod_gbw_bytes)


def parse_orca_cis(cis_fn):
    """
    Read binary CI vector file from ORCA.
        Loosly based on TheoDORE 1.7.1, Authors: S. Mai, F. Plasser
        https://sourceforge.net/p/theodore-qc
    """
    cis_handle = open(cis_fn, "rb")
    # self.log(f"Parsing CI vectors from {cis_handle}")

    # The header consists of 9 4-byte integers, the first 5 of which give useful info.
    nvec = struct.unpack("i", cis_handle.read(4))[0]
    # [0] index of first alpha occ,  is equal to number of frozen alphas
    # [1] index of last  alpha occ
    # [2] index of first alpha virt
    # [3] index of last  alpha virt, header[3]+1 is equal to number of bfs
    # [4] index of first beta  occ,  for restricted equal to -1
    # [5] index of last  beta  occ,  for restricted equal to -1
    # [6] index of first beta  virt, for restricted equal to -1
    # [7] index of last  beta  virt, for restricted equal to -1
    header = [struct.unpack("i", cis_handle.read(4))[0] for i in range(8)]

    def parse_header(header):
        first_occ, last_occ, first_virt, last_virt = header
        frozen = first_occ
        occupied = last_occ + 1
        active = occupied - frozen
        mo_num = last_virt + 1
        virtual = mo_num - first_virt
        return frozen, active, occupied, virtual

    a_frozen, a_active, a_occupied, a_virtual = parse_header(header[:4])
    b_header = parse_header(header[4:])
    unrestricted = all([bh != -1 for bh in (b_header)])
    b_frozen, b_active, b_occupied, b_virtual = b_header
    a_lenci = a_active * a_virtual
    b_lenci = b_active * b_virtual
    a_nel = a_frozen + a_active
    b_nel = b_frozen + b_active
    if not unrestricted:
        b_nel = a_nel
    # expect_mult = a_nel - b_nel + 1

    # Loop over states. For non-TDA order is: X+Y of 1, X-Y of 1,
    # X+Y of 2, X-Y of 2, ...
    prev_root = -1
    prev_mult = None
    iroot_triplets = 0

    # Flags that may later be set to True
    triplets = False
    spin_flip = False
    tda = False
    Xs_a = list()
    Ys_a = list()
    Xs_b = list()
    Ys_b = list()

    def parse_coeffs(lenci, frozen, occupied, virtual):
        coeffs = struct.unpack(lenci * "d", cis_handle.read(lenci * 8))
        coeffs = np.array(coeffs).reshape(-1, virtual)
        # create full array, i.e nocc x nvirt
        coeffs_full = np.zeros((occupied, virtual))
        coeffs_full[frozen:] = coeffs
        return coeffs_full

    def handle_X_Y(root_updated, Xs, Ys, coeffs):
        # When the root was not incremented compared to the previous root we have
        # just parsed X-Y (and parsed X+Y before.)
        #
        # We can recover the separate X and Y vectors by first computing X as
        #   X = (X + Y + X - Y) / 2
        # and then
        #   Y = X + Y - X
        if root_updated:
            X_plus_Y = Xs[-1]  # Parsed in previous cycle
            X_minus_Y = coeffs  # Parsed in current cycle
            X = 0.5 * (X_plus_Y + X_minus_Y)
            Y = X_plus_Y - X
            # Update the X and Y vectors that were already saved with their correct values.
            Xs[-1] = X
            Ys[-1] = Y
        # When the root was incremented we either have a TDA calculation without Y or
        # we parsed X-Y in the previous cycle and now moved to a new root.
        else:
            Xs.append(coeffs)
            Ys.append(np.zeros_like(coeffs))

    for ivec in range(nvec):
        # Header of each vector, contains 6 4-byte ints.
        ncoeffs, _, mult, _, iroot, _ = struct.unpack("iiiiii", cis_handle.read(24))

        # Check if we deal with a spin-flip calculation. There, excitations are from
        # α_activate -> β_virtual.
        if unrestricted and ncoeffs == (a_active * b_virtual):
            unrestricted = False  # Don't expect β_active -> β_virtual.
            spin_flip = True
            a_lenci = ncoeffs
            a_virtual = b_virtual
            warnings.warn(
                "Spin-flip calculation detected. Pysisyphus can parse it, but "
                "the transition density matrix is not yet handled properly by the "
                "OverlapCalculator-class or the Wavefunction-class!"
            )

        if prev_mult is None:
            prev_mult = mult

        # 2 x 8 bytes unknown?!
        ene, _ = struct.unpack("dd", cis_handle.read(16))

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

        root_updated = prev_root == iroot

        # self.log(f"{ivec=}, {nele=}, {mult=}, {iroot=}, {root_updated=}")
        # Then come nact * nvirt 8-byte doubles with the coefficients
        coeffs_a = parse_coeffs(a_lenci, a_frozen, a_occupied, a_virtual)
        handle_X_Y(root_updated, Xs_a, Ys_a, coeffs_a)
        if unrestricted:
            coeffs_b = parse_coeffs(b_lenci, b_frozen, b_occupied, b_virtual)
            handle_X_Y(root_updated, Xs_b, Ys_b, coeffs_b)

        # Somehow ORCA stops to update iroot correctly after the singlet states.
        if (mult == 3) and (tda or (ivec % 2) == 1):
            iroot_triplets += 1

        prev_root = iroot
        prev_mult = mult
    # Verify, that we are at the EOF. We request 1 byte, but we only get 0.
    assert len(cis_handle.read(1)) == 0
    cis_handle.close()

    # Convert everything to numpy arrays.
    Xs_a, Ys_a, Xs_b, Ys_b = [np.array(_) for _ in (Xs_a, Ys_a, Xs_b, Ys_b)]

    def handle_triplets(Xs, Ys):
        assert (len(Xs) % 2) == 0
        states = len(Xs) // 2
        Xs = Xs[states:]
        Ys = Ys[states:]
        return Xs, Ys

    # Only return triplet states if present
    if triplets:
        Xs_a, Ys_a = handle_triplets(Xs_a, Ys_a)
        assert len(Xs_b) == 0
        assert len(Ys_b) == 0

    # Beta part will be empty
    if not unrestricted:
        assert len(Xs_b) == 0
        assert len(Ys_b) == 0
        Xs_b = np.zeros_like(Xs_a)
        Ys_b = np.zeros_like(Xs_b)

    return Xs_a, Ys_a, Xs_b, Ys_b


@file_or_str(".log", ".out")
def parse_orca_all_energies(text, triplets=False, do_tddft=False):
    energy_re = r"FINAL SINGLE POINT ENERGY\s*([-\.\d]+)"
    energy_mobj = re.search(energy_re, text)
    gs_energy = float(energy_mobj.groups()[0])
    all_energies = [gs_energy]

    if do_tddft:
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


def get_name(text: bytes):
    """Return string that comes before first \x00 character & offset."""
    until = text.find(b"\x00")
    return text[:until].decode(), until


@file_or_str(".densities", mode="rb")
def parse_orca_densities(text: bytes):
    handle = io.BytesIO(text)

    # Determine file size
    handle.seek(0, 2)
    file_size = handle.tell()
    handle.seek(0, 0)

    offset, _ = struct.unpack(
        "ii", handle.read(8)
    )  # Don't know about the second integer
    # print("offset", offset)
    dens_size = offset - 8
    assert dens_size % 8 == 0
    dens_floats = dens_size // 8
    # print(f"Expecting {dens_floats} density doubles")
    densities = struct.unpack("d" * dens_floats, handle.read(dens_size))
    ndens = struct.unpack("i", handle.read(4))[0]

    # Block of 512 bytes meta data. I don't really know what is contained in there.
    meta = handle.read(512)
    base_name, _ = get_name(meta)

    # Now multiple 512 byte blocks for each density follow
    dens_names = list()
    for i in range(ndens):
        dens_meta = handle.read(512)
        dens_name, _ = get_name(dens_meta)
        dens_names.append(dens_name)
        # don't know about the first item, 2nd items seems to 0
        _, _, nao1, nao2 = struct.unpack("iiii", handle.read(16))
        assert nao1 == nao2
        _ = struct.unpack("b", handle.read(1))[0]  # 0 byte
        assert _ == 0
    dens_exts = [Path(dens_name).suffix[1:] for dens_name in dens_names]

    # Verify that we parsed to whole file
    assert file_size - handle.tell() == 0
    handle.close()

    # Construct density matrices
    assert dens_floats % ndens == 0
    nao = int(sqrt(dens_floats // ndens))
    dens_shape = (nao, nao)
    densities = np.array(densities).reshape(ndens, *dens_shape)

    dens_dict = {dens_ext: dens for dens_ext, dens in zip(dens_exts, densities)}
    # This check could be removed but I'll keep if for now, so I only deal with
    # known densities.
    # scfp : HF/DFT electronic density
    # scfr : HF/DFT spin density
    # cisp : TDA/TD-DFT/CIS electronic density
    # cisr : TDA/TD-DFT/CIS spin density
    assert set(dens_dict) <= {"scfp", "scfr", "cisp", "cisr"}

    return dens_dict


def get_exc_ens_fosc(wf_fn, cis_fn, log_fn):
    wf = Wavefunction.from_file(wf_fn)
    Xa, Ya, Xb, Yb = parse_orca_cis(cis_fn)
    all_energies = parse_orca_all_energies(log_fn, do_tddft=True)
    Xa, Ya = norm_ci_coeffs(Xa, Ya)
    exc_ens = all_energies[1:] - all_energies[0]
    tdens = wf.get_transition_dipole_moment(Xa + Ya)
    warnings.warn("Only alpha TDM is currently taken into account!")
    fosc = wf.oscillator_strength(exc_ens, tdens)
    return exc_ens, fosc


class ORCA(OverlapCalculator):
    conf_key = "orca"
    _set_plans = (
        "gbw",
        "out",
        "cis",
        "densities",
        ("molden", "mwfn_wf"),
        "json",
    )

    def __init__(
        self,
        keywords,
        blocks="",
        gbw=None,
        do_stable=False,
        numfreq=False,
        json_dump=None,
        wavefunction_dump=True,
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
        numfreq : bool, optional
            Use numerical frequencies instead of analytical ones.
        json_dump : bool, optional, deprecated
            Use 'wavefunction_dump' instead.
        wavefunction_dump : bool, optional
            Whether to dump the wavefunction to JSON via orca_2json. The JSON can become
            very large in calculations comprising many basis functions.
        """
        super().__init__(**kwargs)

        self.keywords = keywords.lower()
        self.blocks = blocks.lower()
        self.gbw = gbw
        self.do_stable = bool(do_stable)
        self.freq_keyword = "numfreq" if numfreq else "freq"
        if json_dump is not None:
            warnings.warn(
                "Use of 'json_dump' is deprecated! Use 'wavefunction_dump' instead!",
                DeprecationWarning,
            )
            wavefunction_dump = json_dump
        self.wavefunction_dump = bool(wavefunction_dump)

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
            "densities:orca.densities",
            "json",
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
        results["all_energies"] = self.parse_all_energies()
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
        results = self.store_and_track(
            results, self.get_hessian, atoms, coords, **prepare_kwargs
        )
        return results

    def get_stored_wavefunction(self, **kwargs):
        return self.load_wavefunction_from_file(self.json, **kwargs)

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
            cmd = "orca_2mkl orca -molden"
            self.popen(cmd, cwd=path)
            shutil.copy(path / "orca.molden.input", path / "orca.molden")

        if self.wavefunction_dump:
            # Will silently fail with ECPs
            cmd = "orca_2json orca"
            proc = self.popen(cmd, cwd=path)
            if (ret := proc.returncode) != 0:
                self.log(f"orca_2json call failed with return-code {ret}!")

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

    @staticmethod
    def parse_cis(cis):
        """Simple wrapper of external function.

        Currently, only returns Xα and Yα.
        """
        return parse_orca_cis(cis)[:2]

    @staticmethod
    def parse_gbw(gbw_fn):
        return parse_orca_gbw(gbw_fn)

    @staticmethod
    def set_mo_coeffs_in_gbw(in_gbw_fn, out_gbw_fn, mo_coeffs):
        """See self.parse_gbw."""
        warnings.warn(
            "Previously, this method expected MO coefficients in rows! "
            "Did you update your code?!"
        )
        return set_mo_coeffs_in_gbw(in_gbw_fn, out_gbw_fn, mo_coeffs)

    def parse_all_energies(self, text=None, triplets=None):
        if text is None:
            with open(self.out) as handle:
                text = handle.read()

        if triplets is None:
            triplets = self.triplets

        return parse_orca_all_energies(text, triplets, self.do_tddft)

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
        C, _ = self.parse_gbw(self.gbw)
        all_energies = self.parse_all_energies()
        return C, X, Y, all_energies

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
