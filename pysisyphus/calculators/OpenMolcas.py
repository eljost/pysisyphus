import os
from pathlib import Path
import re

import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG
from pysisyphus.xyzloader import make_xyz_str


class OpenMolcas(Calculator):

    conf_key = "openmolcas"

    def __init__(self, basis, inporb, roots, mdrlxroot,
                 supsym=None, track=True, **kwargs):
        super(OpenMolcas, self).__init__(**kwargs)

        assert self.pal == 1, (
            "RI SA-CASSCF analytical gradients do not work correctly in "
            "parallel (yet). Consider using pal=1 instead of the current "
           f"pal={self.pal}!"
        )

        self.basis = basis
        self.inporb = Path(inporb).absolute()
        self.roots = roots
        self.mdrlxroot = mdrlxroot
        self.supsym = self.build_supsym_str(supsym)
        self.track = track

        self.to_keep = ("RasOrb", "out", "in", "JobIph", "rasscf.molden")
        self.jobiph = ""

        self.inp_fn = "openmolcas.in"
        self.out_fn = "openmolcas.out"
        self.float_regex = r"([\d\.\-E]+)"

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
         mdrlxroot
          {mdrlxroot}
         {supsym}

        >> copy $Project.JobIph $CurrDir/$Project.JobIph

        {rassi}

        &alaska
         pnew
        """

        self.parser_funcs = {
            "grad": self.parse_gradient,
        }

        self.base_cmd = self.get_cmd()

    def reattach(self, last_calc_cycle):
        self.inporb = self.make_fn("RasOrb", last_calc_cycle)
        self.jobiph = self.make_fn("JobIph", last_calc_cycle)
        self.log(f"restarted. using {self.inporb}, {self.jobiph}")

    def build_supsym_str(self, supsym):
        """Can handle only one subgroup for now."""
        if not supsym:
            return ""

        num_orbitals = len(supsym.split())
        return f"supsym\n1\n{num_orbitals} {supsym};"

    def build_rassi_str(self):
        # In the first iteration self.jobiph isn't set yet.
        if (not self.track) or (self.calc_counter == 0):
            return ""
        else:
            # JOB001 corresponds to the current iteration,
            # JOB002 to the previous iteration.
            return f"""
            >> copy $Project.JobIph JOB001
            >> copy {self.jobiph} JOB002
            &rassi
             track
            """

    def get_pal_env(self):
        env_copy = os.environ.copy()
        env_copy["MOLCAS_NPROCS"] = str(self.pal)

        return env_copy

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
                                        mdrlxroot=self.mdrlxroot,
                                        supsym=self.supsym,
                                        rassi=self.build_rassi_str(),
        )
        return inp

    def get_energy(self, atoms, coords):
        results = self.get_forces(atoms, coords)
        del results["forces"]
        return results

    def get_forces(self, atoms, coords):
        self.log(f"using inporb: {self.inporb}")
        inp = self.prepare_input(atoms, coords)
        add_args = ("-clean", "-oe", self.out_fn)
        env_copy = self.get_pal_env()
        env_copy["MOLCAS_PROJECT"] = f"{self.name}_{self.calc_counter}"
        kwargs = {
            "calc": "grad",
            "add_args": add_args,
            "env": env_copy,
        }
        results = self.run(inp, **kwargs)
        return results

    def keep(self, path):
        kept_fns = super().keep(path)
        self.inporb = kept_fns["rasorb"]
        # Keep references to the .JobIph file to be used in
        # &rassi to track our root in a state average
        # calculation.
        self.jobiph = kept_fns["jobiph"]
        self.log(f"current JobIph is {self.jobiph}")

    def parse_energies(self, text):
        # Energy of root for which gradient was computed
        energy_regex = r"RASSCF state energy =\s*" + self.float_regex
        energy = float(re.search(energy_regex, text).groups()[0])

        # All state average energies
        root_re = "RASSCF root number.+Total energy.+?" + self.float_regex
        matches = re.findall(root_re, text)
        sa_energies = np.array(matches, dtype=float)

        return energy, sa_energies

    def parse_gradient(self, path):
        results = {}
        gradient_fn = os.path.join(path, self.out_fn)
        with open(gradient_fn) as handle:
            text = handle.read()

        # Search for the block containing the gradient table
        regex = r"Molecular gradients(.+?)--- Stop Module:\s*alaska"
        floats = [self.float_regex for i in range(3)]
        line_regex = r"([A-Z\d]+)\s*" + r"\s*".join(floats)

        mobj = re.search(regex, text, re.DOTALL)
        gradient = list()
        for line in mobj.groups()[0].split("\n"):
            # Now look for the lines containing the gradient
            mobj = re.match(line_regex, line.strip())
            if not mobj:
                continue
            # Discard first column (atom+number)
            gradient.append(mobj.groups()[1:])
        gradient = np.array(gradient, dtype=float).flatten()

        if self.track and self.calc_counter > 0:
            self.parse_rassi_track(path)

        energy, sa_energies = self.parse_energies(text)

        results["energy"] = energy
        # results["sa_energies"] = sa_energies
        results["forces"] = -gradient

        return results

    def parse_rassi_track(self, path):
        gradient_fn = path / self.out_fn
        with open(gradient_fn) as handle:
            text = handle.read()
        track_re = r"Initial root:\s*(\d+)\s*Overlaps with current " \
                   r"states:(.+)New root:\s*(\d+)"
        #overlap_re = "OVERLAP MATRIX FOR THE ORIGINAL STATES:(.+?)##"
        mobj = re.search(track_re, text, re.DOTALL)

        initial_root, overlaps, new_root = mobj.groups()
        overlaps = np.array(overlaps.strip().split(), dtype=float).reshape(-1, 2)
        # Filters for overlaps > 10% (0.1**2 ~ 0.31622)
        thresh = 0.1
        inds = np.where(np.abs(overlaps[:,1]) > thresh**0.5)
        ov_perc_str = ", ".join([f"{nr:.0f}: {ov**2:.2%}"
                                 for nr, ov in overlaps[inds]])
        self.log(f"Overlaps between previous root {initial_root} and "
                 f"new roots bigger {thresh:.0%}:  {ov_perc_str}. Will "
                 f"use root {new_root} for the following gradient calculation.")
        if new_root != initial_root:
            self.log("Found a root flip!")
        self.mdrlxroot = new_root

    def __str__(self):
        return "OpenMolcas calculator"
