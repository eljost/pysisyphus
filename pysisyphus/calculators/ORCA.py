#!/usr/bin/env python3

import glob
import logging
import os
from pathlib import Path
import re
import struct
import subprocess

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.OverlapCalculator import OverlapCalculator
from pysisyphus.calculators.WFOWrapper import WFOWrapper


def make_sym_mat(table_block):
    mat_size = int(table_block[1])
    # Orca prints blocks of 5 columns
    arr = np.array(table_block[2:], dtype=np.float)
    assert(arr.size == mat_size**2)
    block_size = 5*mat_size
    cbs = [arr[i*block_size:(i+1)*block_size].reshape(mat_size, -1)
           for i in range(arr.size // block_size + 1)
    ]
    return np.concatenate(cbs, axis=1)


class ORCA(OverlapCalculator):

    conf_key = "orca"

    def __init__(self, keywords, gbw="", blocks="", **kwargs):
        super().__init__(**kwargs)

        self.keywords = keywords
        # Only call when we are not restarting
        if not ("last_calc_cycle" in kwargs):
            self.set_moinp_str(gbw)
        self.blocks = blocks

        assert (("pal" not in keywords.lower())
                and ("nprocs" not in blocks.lower())), "PALn/nprocs not " \
                "allowed! Use 'pal: n' in the 'calc' section instead."

        self.to_keep = ("inp", "out:orca.out", "gbw", "engrad", "hessian", "cis")
        self.do_tddft = False
        if "tddft" in self.blocks:
            self.do_tddft = True
            try:
                self.root = int(re.search("iroot\s*(\d+)", self.blocks).group(1))
            except AttributeError:
                self.log("Doing TDA/TDDFT calculation without gradient.")
        if self.track:
            # Defer the creation of the WFOWrapper object until the after the
            # calculation as then the number of occupied/virtual MOs are
            # available.
            self.wfow = None
        self.inp_fn = "orca.inp"
        self.out_fn = "orca.out"

        self.orca_input="""!{keywords} {calc_type}
        {moinp}

        %pal nprocs {pal} end

        {blocks}

        *xyz {charge} {mult}
        {coords}
        *
        """

        self.parser_funcs = {
            "grad": self.parse_engrad,
            "hessian": self.parse_hessian,
            "noparse": lambda path: None,
        }

        self.base_cmd = self.get_cmd("cmd")

    def reattach(self, last_calc_cycle):
        # Use the latest .gbw
        gbw = self.make_fn("gbw", last_calc_cycle)
        self.log(f"restarted. using {gbw}")
        self.set_moinp_str(gbw)

    def set_moinp_str(self, gbw):
        if not gbw:
            self.moinp = ""
            self.gbw = ""
        else:
            self.moinp = f"""!moread
            %moinp "{gbw}" """
            self.gbw = gbw

    def prepare_input(self, atoms, coords, calc_type):
        coords = self.prepare_coords(atoms, coords)
        if self.gbw:
            self.log(f"using {self.gbw}")
        else:
            self.log("using initial guess provided by ORCA")
        if calc_type == "noparse":
            calc_type = ""
        inp = self.orca_input.format(
                                keywords=self.keywords,
                                calc_type=calc_type,
                                moinp=self.moinp,
                                pal=self.pal,
                                blocks=self.get_block_str(),
                                coords=coords,
                                charge=self.charge,
                                mult=self.mult,
        )
        return inp

    def get_block_str(self):
        block_str = ""
        if self.track:
            block_str = re.sub("iroot\s+(\d+)", f"iroot {self.root}", self.blocks)
            self.log(f"Using iroot '{self.root}' for excited state gradient.")
        return block_str

    def get_energy(self, atoms, coords):
        results = self.get_forces(atoms, coords)
        del results["forces"]
        return results

    def get_forces(self, atoms, coords):
        calc_type = "engrad"
        inp = self.prepare_input(atoms, coords, calc_type)
        kwargs = {
            "calc": "grad",
        }
        results = self.run(inp, **kwargs)
        if self.track:
            if self.track_root(atoms, coords):
                # Redo the calculation with the updated root
                results = self.get_forces(atoms, coords)
        return results

    def get_hessian(self, atoms, coords):
        calc_type = "freq"
        inp = self.prepare_input(atoms, coords, calc_type)
        results = self.run(inp, calc="hessian")
        return results

    def run_calculation(self, atoms, coords):
        """Basically some kind of dummy method that can be called
        to execute ORCA with the stored cmd of this calculator."""
        inp = self.prepare_input(atoms, coords, "noparse")
        kwargs = {
                "calc": "noparse",
        }
        results = self.run(inp, **kwargs)
        if self.track:
            self.store_overlap_data(atoms, coords)
        return results

    def parse_hessian(self, path):
        results = {}
        hessian_fn = glob.glob(os.path.join(path, "*.hess"))
        assert(len(hessian_fn) == 1)
        hessian_fn = hessian_fn[0]
        if not hessian_fn:
            raise Exception("ORCA calculation failed.")
        with open(hessian_fn) as handle:
            text = handle.read()

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
        mass_xyz_line = (pp.Word(pp.alphas) + float_ +
                         pp.Group(pp.OneOrMore(float_))
        )

        block_name = pp.Word(pp.alphas + "$_")
        block_length = integer

        block_int = block_name + block_length
        block_float = block_name + float_
        block_table = block_name + integer + pp.OneOrMore(scientific_block)
        block_table_two_int = (block_name + integer + pp.Suppress(integer)
                               + pp.OneOrMore(scientific_block)
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

        parser = (block_name + act_atom + act_coord + act_energy
                  + hessian + vib_freqs + normal_modes
                  + pp.OneOrMore(comment_line) + atoms
        )
        parsed = parser.parseString(text)
        results["hessian"] = make_sym_mat(parsed["hessian"])

        logging.warning("Hacky orca energy parsing in orca hessian calculation!")
        orca_log_fn = os.path.join(path, self.out_fn)
        with open(orca_log_fn) as handle:
            log_text = handle.read()

        energy_re = "FINAL SINGLE POINT ENERGY\s*([-\.\d]+)"
        energy_mobj = re.search(energy_re, log_text)
        energy = float(energy_mobj.groups()[0])
        results["energy"] = energy

        if self.do_tddft:
            """FIXME: Store the right energy etc. similar to
            parse_engrad."""
            raise Exception("Proper handling of TDDFT and hessian "
                            " is not yet implemented.")
        return results



    def parse_engrad(self, path):
        results = {}
        engrad_fn = glob.glob(os.path.join(path, "*.engrad"))
        if not engrad_fn:
            raise Exception("ORCA calculation failed.")
            self.print_out_fn(path)

        assert(len(engrad_fn) == 1)
        engrad_fn = engrad_fn[0]
        with open(engrad_fn) as handle:
            engrad = handle.read()
        engrad = re.findall("([\d\-\.]+)", engrad)
        atoms = int(engrad.pop(0))
        energy = float(engrad.pop(0))
        force = -np.array(engrad[:3*atoms], dtype=np.float)
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
        nvec  = struct.unpack('i', cis_handle.read(4))[0]
        # header array contains:
        # [0] index of first alpha occ,  is equal to number of frozen alphas
        # [1] index of last  alpha occ
        # [2] index of first alpha virt
        # [3] index of last  alpha virt, header[3]+1 is equal to number of bfs
        # [4] index of first beta  occ,  for restricted equal to -1
        # [5] index of last  beta  occ,  for restricted equal to -1
        # [6] index of first beta  virt, for restricted equal to -1
        # [7] index of last  beta  virt, for restricted equal to -1
        header = [struct.unpack('i', cis_handle.read(4))[0]
                  for i in range(8)]

        if any([flag != -1 for flag in header[4:8]]):
            raise Exception("parse_cis, no support for unrestricted MOs")

        nfrzc = header[0]
        nocc = header[1] + 1
        nact = nocc - nfrzc
        nmo  = header[3] + 1
        nvir = nmo - header[2]
        lenci = nact * nvir
        self.log(f"nmo = {nmo}, nocc = {nocc}, nact = {nact}, nvir = {nvir}")

        # Loop over states. For non-TDA order is: X+Y of 1, X-Y of 1,
        # X+Y of 2, X-Y of 2, ...
        prevroot = -1
        istate = 0
        coeffs = list()
        for ivec in range(nvec):
            # header of each vector
            # contains 6 4-byte ints, then 1 8-byte double, then 8 byte unknown
            nele,d1,mult,d2,iroot,d3 = struct.unpack('iiiiii', cis_handle.read(24))
            ene,d3 = struct.unpack('dd', cis_handle.read(16))
            self.log(f"nele = {nele}, mult = {mult}, iroot = {iroot}")
            # then comes nact * nvirt 8-byte doubles with the coefficients
            coeff = struct.unpack(lenci*'d', cis_handle.read(lenci*8))
            coeff = np.array(coeff).reshape(-1, nvir)
            # create full array, i.e nocc x nvirt
            coeff_full = np.zeros((nocc, nvir))
            coeff_full[nfrzc:] = coeff

            # in this case, we have a non-TDA state!
            # and we need to compute (prevvector+currentvector)/2 = X vector
            if prevroot == iroot:
                self.log('Constructing X-vector of RPA state')
                x_plus_y = coeffs[-1]
                x_minus_y = coeff_full
                x = 0.5*(x_plus_y + x_minus_y)
                coeffs[-1] = x
            else:
                coeffs.append(coeff_full)

            prevroot=iroot
        cis_handle.close()
        return np.array(coeffs)

    def parse_gbw(self, gbw_fn):
        """Adapted from
        https://orcaforum.cec.mpg.de/viewtopic.php?f=8&t=3299&hilit=perl&start=20#p13273

        The first 5 long int values represent pointers into the file:

        Pointer @+0:  Internal ORCA data structures
        Pointer @+8:  Geometry
        Pointer @+16: BasisSet
        Pointer @+24: Orbitals
        Pointer @+32: ECP data
        """

        with open(gbw_fn, 'rb') as handle:
            handle.seek(24)
            offset = struct.unpack('<q', handle.read(8))[0]
            handle.seek(offset)
            operators = struct.unpack('<i', handle.read(4))[0]
            dimension = struct.unpack('<i', handle.read(4))[0]

            # print('Offset: {}'.format(offset))
            # print('Number of Operators: {}'.format(operators))
            # print('Basis Dimension: {}'.format(dimension))

            coeffs_fmt = "<" + dimension**2 * "d"

            assert operators == 1, "Unrestricted case is not implemented!"

            for i in range(operators):
                #print('\nOperator: {}'.format(i))
                coeffs = struct.unpack(coeffs_fmt, handle.read(8*dimension**2))
                occupations = struct.iter_unpack('<d', handle.read(8*dimension))
                energies = struct.iter_unpack('<d', handle.read(8*dimension))
                irreps = struct.iter_unpack('<i', handle.read(4*dimension))
                cores = struct.iter_unpack('<i', handle.read(4*dimension))

                coeffs = np.array(coeffs).reshape(-1, dimension).T

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
            return coeffs

    def parse_all_energies(self):
        with open(self.out) as handle:
            text = handle.read()

        energy_re = "FINAL SINGLE POINT ENERGY\s*([-\.\d]+)"
        energy_mobj = re.search(energy_re, text)
        gs_energy = float(energy_mobj.groups()[0])
        all_energies = [gs_energy]

        if self.do_tddft:
            tddft_re = re.compile("STATE\s*\d+:\s*E=\s*([\d\.]+)\s*au")
            excitation_ens = [float(en) for en in tddft_re.findall(text)]
            all_energies += (np.array(excitation_ens) + gs_energy).tolist()
        all_energies = np.array(all_energies)
        return all_energies

    def parse_mo_numbers(self, out_fn):
        with open(out_fn) as handle:
            text = handle.read()
        electron_re = "NEL\s*....\s*(\d+)"
        electrons = int(re.search(electron_re, text)[1])
        assert(electrons % 2 == 0), "unrestricted is not yet supported!"
        occ_num = int(electrons / 2)

        mo_re = "Dim\s*....\s*(\d+)"
        mo_num = int(re.search(mo_re, text)[1])
        virt_num = mo_num - occ_num
        self.log(f"found {electrons} electrons, {mo_num} MOs, with "
                 f"{occ_num} occupied and {virt_num} virtual.")
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
        self.mo_coeffs = self.parse_gbw(self.gbw)

    def set_ci_coeffs(self, ci_coeffs=None, cis=None):
        if ci_coeffs is not None:
            self.ci_coeffs = ci_coeffs
            return
        if not cis and self.cis:
            cis = self.cis
        else:
            raise Exception("Got no .cis file to parse!")
        self.log(f"Setting CI coefficients from {cis}.")
        self.ci_coeffs = self.parse_cis(cis)

    def prepare_overlap_data(self):
        # Create the WFOWrapper object if it is not already there
        if self.wfow == None:
            occ_num, virt_num = self.parse_mo_numbers(self.out)
            self.wfow = WFOWrapper(occ_num, virt_num, calc_number=self.calc_number,
                                   basis=None, charge=None, out_dir=self.out_dir)
        # Parse eigenvectors from tda/tddft calculation
        ci_coeffs = self.parse_cis(self.cis)
        # Parse mo coefficients from gbw file and write a 'fake' turbomole
        # mos file.
        mo_coeffs = self.parse_gbw(self.gbw)
        fake_mos_fn = Path(self.make_fn("mos"))
        if not fake_mos_fn.exists:
            fake_mos_str = self.wfow.fake_turbo_mos(mo_coeffs)
            with open(fake_mos_fn, "w") as handle:
                handle.write(fake_mos_str)
        else:
            self.log("Skipping creation of MOs in TURBOMOLE format, as the file "
                     f"'{fake_mos_fn}' already exists.")

        all_energies = self.parse_all_energies()
        return fake_mos_fn, mo_coeffs, ci_coeffs, all_energies

    def keep(self, path):
        kept_fns = super().keep(path)
        self.set_moinp_str(kept_fns["gbw"])
        self.out = kept_fns["out"]
        if self.do_tddft:
            self.cis = kept_fns["cis"]

    def __str__(self):
        return f"ORCA({self.name})"


if __name__ == "__main__":
    from pysisyphus.helpers import geom_from_library
    geom = geom_from_library("dieniminium_cation_s1_opt.xyz")
    keywords = "BP86 def2-SV(P)"
    blocks = "tddft iroot 1 end"
    charge = 1
    mult = 1
    orca = ORCA(keywords, blocks, charge=charge, mult=mult)
    """
    geom.set_calculator(orca)
    forces = geom.forces
    print(forces)
    """
    res = orca.parse_engrad("/scratch/test/pysis_orca/neu")
    orca.set_moinp_str("")
    print(orca.moinp)
    orca.set_moinp_str("path/to/gbw")
    print(orca.moinp)
