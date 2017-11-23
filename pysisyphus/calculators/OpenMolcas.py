#!/usr/bin/env python3

import os
import re

import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.config import Config
from pysisyphus.constants import BOHR2ANG

from pysisyphus.xyzloader import make_xyz_str


class OpenMolcas(Calculator):

    def __init__(self, basis, inporb, roots, rlxroot, **kwargs):
        super(OpenMolcas, self).__init__(**kwargs)

        self.basis = basis
        self.inporb = inporb
        self.roots = roots
        self.rlxroot = rlxroot

        self.inp_fn = "openmolcas.in"
        self.out_fn = "openmolcas.out"
        self.float_regex = "([\d\.\-]+)"

        self.openmolcas_input = """
        >> copy {inporb}  $Project.RasOrb
        &gateway
         coord
          {xyz_str}
         basis
          {basis}
         group
          nosym
        ricd
 
        &seward
         doanalytical

        &rasscf
         charge
          {charge}
         spin
          {mult}
         fileorb
          $Project.RasOrb
         thrs
          1.0e-6,1.0e-2,1.0e-2
         ciroot
          {roots} {roots} 1
         rlxroot
          {rlxroot}

        &alaska
         pnew
        """

        self.parser_funcs = {
            "grad": self.parse_gradient,
        }

        self.base_cmd = Config["openmolcas"]["cmd"]

    def prepare_coords(self, atoms, coords):
        coords = coords * BOHR2ANG
        return make_xyz_str(atoms, coords.reshape((-1, 3)))

    def prepare_input(self, atoms, coords):
        xyz_str = self.prepare_coords(atoms, coords)
        inp = self.openmolcas_input.format(
                                        inporb=self.inporb,
                                        xyz_str=xyz_str,
                                        basis=self.basis,
                                        charge=self.charge,
                                        mult=self.mult,
                                        roots=self.roots,
                                        rlxroot=self.rlxroot,
        )
        return inp

    def get_forces(self, atoms, coords):
        #self.logger.debug(f"Using inporb: {self.inporb}")
        inp = self.prepare_input(atoms, coords)
        add_args = ("-clean", "-oe", self.out_fn)
        env = os.environ.copy()
        env["MOLCAS_PROJECT"] = f"{self.name}_{self.counter}"
        results = self.run(inp, calc="grad", add_args=add_args, env=env)
        return results

    def keep(self, path):
        kept_fns = super().keep(path, ("RasOrb", "out"))
        self.inporb = kept_fns["RasOrb"]


    def parse_energies(self, text):
        # Energy of root for which gradient was computed
        energy_regex = "RASSCF state energy =\s*" + self.float_regex
        energy = float(re.search(energy_regex, text).groups()[0])

        # All state average energies
        root_re = "RASSCF root number.+Total energy.+?" + self.float_regex
        matches = re.findall(root_re, text)
        sa_energies = np.array(matches, dtype=np.float)

        return energy, sa_energies

    def parse_gradient(self, path):
        results = {}
        gradient_fn = os.path.join(path, self.out_fn)
        with open(gradient_fn) as handle:
            text = handle.read()

        # Search for the block containing the gradient table
        regex = "Molecular gradients(.+?)--- Stop Module:\s*alaska"
        floats = [self.float_regex for i in range(3)]
        line_regex = "([A-Z\d]+)\s*" + "\s*".join(floats)

        mobj = re.search(regex, text, re.DOTALL)
        gradient = list()
        for line in mobj.groups()[0].split("\n"):
            # Now look for the lines containing the gradient
            mobj = re.match(line_regex, line.strip())
            if not mobj:
                continue
            # Discard first column (atom+number)
            gradient.append(mobj.groups()[1:])
        gradient = np.array(gradient, dtype=np.float).flatten()

        energy, sa_energies = self.parse_energies(text)

        results["energy"] = energy
        results["sa_energies"] = sa_energies
        results["forces"] = -gradient

        return results

    def __str__(self):
        return "OpenMolcas calculator"


if __name__ == "__main__":
    from pysisyphus.helpers import geom_from_library
    fileorb = "/scratch/test/ommin/excrp.es_opt.RasOrb"
    basis = "6-31G*"
    roots = 2
    rlxroot = 2
    om = OpenMolcas(basis, fileorb, roots, rlxroot)
    geom = geom_from_library("dieniminium_cation_s1_opt.xyz")
    geom.set_calculator(om)
    #print(geom.forces)
    om.parse_gradient("/scratch/test/satest")
