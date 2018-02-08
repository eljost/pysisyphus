#!/usr/bin/env python3

import glob
import logging
import os
from pathlib import Path
import re
import subprocess

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.config import Config

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


class ORCA(Calculator):

    def __init__(self, keywords, gbw="", blocks="", **kwargs):
        super(ORCA, self).__init__(**kwargs)

        self.keywords = keywords
        # Only call when we are not restarting
        if not ("last_calc_cycle" in kwargs):
            self.set_moinp_str(gbw)
        self.blocks = blocks

        self.to_keep = ("out", "gbw", "engrad", "hessian")
        self.do_tddft = False
        if "tddft" in self.blocks:
            self.do_tddft = True
            self.iroot = int(re.search("iroot\s*(\d+)", self.blocks).group(1))
        self.inp_fn = "orca.inp"
        self.out_fn = "orca.out"

        self.orca_input="""!{keywords} {calc_type}
        {moinp}

        {blocks}

        *xyz {charge} {mult}
        {coords}
        *
        """

        self.parser_funcs = {
            "grad": self.parse_engrad,
            "hessian": self.parse_hessian,
        }

        self.base_cmd = Config["orca"]["cmd"]

    def reattach(self, last_calc_cycle):
        # Use the latest .gbw
        gbw = self.make_fn("gbw", last_calc_cycle, True)
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
        inp = self.orca_input.format(
                                keywords=self.keywords,
                                calc_type=calc_type,
                                moinp=self.moinp,
                                blocks=self.blocks,
                                coords=coords,
                                charge=self.charge,
                                mult=self.mult,
        )
        return inp

    def get_energy(self, atoms, coords):
        logging.info("orca, energy_calculation!")
        logging.warning("orca energy not implemented properly!")
        logging.warning("Called energy, exiting!")
        import sys; sys.exit()

    def get_forces(self, atoms, coords):
        calc_type = "engrad"
        inp = self.prepare_input(atoms, coords, calc_type)
        results = self.run(inp, calc="grad")
        return results

    def get_hessian(self, atoms, coords):
        calc_type = "freq"
        inp = self.prepare_input(atoms, coords, calc_type)
        results = self.run(inp, calc="hessian")
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

        if self.do_tddft:
            # This sets the proper excited state energy in the
            # results dict and also stores all energies.
            excitation_ens = self.parse_tddft(path)
            # ORCA iroot input is 1 based, so we substract 1 to get
            # the right index here.
            iroot_exc_en = excitation_ens[self.iroot-1]
            gs_energy = results["energy"]
            # Add excitation energy to ground state energy.
            results["energy"] += iroot_exc_en
            all_ens = np.full(len(excitation_ens)+1, gs_energy)
            all_ens[1:] += excitation_ens
            results["tddft_energies"] = all_ens

        return results

    def parse_tddft(self, path):
        results = {}
        orca_out = Path(path) / self.out_fn
        with open(orca_out) as handle:
            text = handle.read()
        tddft_re = re.compile("STATE\s*\d+:\s*E=\s*([\d\.]+)\s*au")
        excitation_ens = [float(en) for en in tddft_re.findall(text)]
        return excitation_ens


    def keep(self, path):
        kept_fns = super().keep(path)
        self.set_moinp_str(kept_fns["gbw"])

    def __str__(self):
        return "ORCA calculator"


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
