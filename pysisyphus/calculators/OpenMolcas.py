import os
from pathlib import Path
import re
import warnings

import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG
from pysisyphus.xyzloader import make_xyz_str


class OpenMolcas(Calculator):
    conf_key = "openmolcas"
    _set_plans = (
        ("rasscf.h5", "inporb"),
        "jobiph",
    )

    def __init__(
        self,
        basis,
        inporb,
        rasscf=None,
        gateway=None,
        mcpdft=None,
        rassi=None,
        nprocs=1,
        track=True,
        omp_var="OMP_NUM_THREADS",
        **kwargs,
    ):
        super(OpenMolcas, self).__init__(**kwargs)

        env = os.environ
        if "MOLCAS" not in env:
            warnings.warn(
                "$MOLCAS environment variable is not set! Couldn't find it in the "
                "environment!"
            )

        self.basis = basis
        inporb = Path(inporb).absolute()
        self.inporb = inporb
        if rasscf is None:
            rasscf = {}
        self.rasscf = rasscf
        if gateway is None:
            gateway = {}
        self.gateway = gateway
        if mcpdft is None:
            mcpdft = {}
        self.mcpdft = mcpdft
        if rassi is None:
            rassi = {}
        self.rassi = rassi

        self.nprocs = nprocs
        self.omp_var = omp_var
        assert self.nprocs == 1, (
            "RI SA-CASSCF analytical gradients do not work correctly in "
            "parallel (yet). Consider using nprocs=1 instead of the current "
            f"nprocs={self.nprocs}!"
        )
        self.track = track

        self.to_keep = (
            "RasOrb",
            "out",
            "in",
            "JobIph",
            "rasscf.molden",
            "rassi.h5",
            "rasscf.h5",
        )
        self.jobiph = ""

        self.inp_fn = "openmolcas.in"
        self.out_fn = "openmolcas.out"
        self.float_regex = r"([\d\.\-E]+)"

        self.openmolcas_input = """
        >> copy {inporb}  $Project.rasscf.{inporb_ext}
        &gateway
         coord
          {xyz_str}
         basis
          {basis}
         group
          nosym
         {gateway_kwargs}

        &seward

        &rasscf
         charge
          {charge}
         spin
          {mult}
         fileorb
          $Project.rasscf.{inporb_ext}
         {rasscf_kwargs}

        {mcpdft}

        >> copy $Project.JobIph $CurrDir/$Project.JobIph

        {rassi}

        {alaska}
        """

        self.parser_funcs = {
            "energy": self.parse_energy,
            "grad": self.parse_gradient,
        }

        self.base_cmd = self.get_cmd()

    def reattach(self, last_calc_cycle):
        self.inporb = self.make_fn("RasOrb", last_calc_cycle)
        self.jobiph = self.make_fn("JobIph", last_calc_cycle)
        self.log(f"restarted. using {self.inporb}, {self.jobiph}")

    def build_str_from_dict(self, dct):
        strs = list()
        for key, val in dct.items():
            if val is None:
                val = ""
            strs.append(f"{key}\n{val}")
        return "\n".join(strs)

    def build_gateway_str(self):
        return self.build_str_from_dict(self.gateway)

    def build_rasscf_str(self):
        return self.build_str_from_dict(self.rasscf)

    def build_rassi_str(self):
        # In the first iteration self.jobiph isn't set yet.
        if self.rassi:
            rassi_keywords = self.build_str_from_dict(self.rassi)
            rassi_str = f"""
            &rassi
             {rassi_keywords}
            """
        elif (not self.track) or (self.calc_counter == 0):
            rassi_str = ""
        else:
            # JOB001 corresponds to the current iteration,
            # JOB002 to the previous iteration.
            rassi_str = f"""
            >> copy $Project.JobIph JOB001
            >> copy {self.jobiph} JOB002
            &rassi
             track
            """
        return rassi_str

    def build_mcpdft_str(self):
        if not self.mcpdft:
            return ""

        mcpdft_kwargs = self.build_str_from_dict(self.mcpdft)
        return f"&mcpdft\n{mcpdft_kwargs}"

    def get_pal_env(self):
        env_copy = os.environ.copy()
        env_copy["MOLCAS_NPROCS"] = str(self.nprocs)
        env_copy[self.omp_var] = str(self.pal)
        # Memory is per process
        env_copy["MOLCAS_MEM"] = str(self.mem * self.pal)

        return env_copy

    def prepare_coords(self, atoms, coords):
        coords = coords * BOHR2ANG
        return make_xyz_str(atoms, coords.reshape((-1, 3)))

    def prepare_input(self, atoms, coords, calc_type):
        self.log(f"using inporb: {self.inporb}")
        xyz_str = self.prepare_coords(atoms, coords)
        alaska_str = "&alaska\npnew" if calc_type == "grad" else ""
        inporb_ext = Path(self.inporb).suffix[1:]  # Drop dot
        inp = self.openmolcas_input.format(
            inporb=self.inporb,
            inporb_ext=inporb_ext,
            xyz_str=xyz_str,
            basis=self.basis,
            charge=self.charge,
            mult=self.mult,
            gateway_kwargs=self.build_gateway_str(),
            rasscf_kwargs=self.build_rasscf_str(),
            mcpdft=self.build_mcpdft_str(),
            rassi=self.build_rassi_str(),
            alaska=alaska_str,
        )
        return inp

    def run_calculation(self, atoms, coords, calc_type="energy"):
        path = self.prepare_path(use_in_run=True)
        inp = self.prepare_input(atoms, coords, calc_type)
        add_args = ("-clean", "-oe", self.out_fn)
        env_copy = self.get_pal_env()
        env_copy["MOLCAS_PROJECT"] = "openmolcas"
        env_copy["MOLCAS_WORKDIR"] = str(path)
        kwargs = {
            "calc": calc_type,
            "add_args": add_args,
            "env": env_copy,
        }
        results = self.run(inp, **kwargs)
        return results

    def get_energy(self, atoms, coords):
        return self.run_calculation(atoms, coords, "energy")

    def get_forces(self, atoms, coords):
        return self.run_calculation(atoms, coords, "grad")

    def get_root(self):
        ras = self.rasscf
        return int(ras.get("mdrlxroot", ras.get("rlxroot", 1)))

    def parse_energies(self, text):
        if self.mcpdft:
            root_re = r"PDFT Root\s*\d+\s*Total energy:\s*" + self.float_regex
            matches = re.findall(root_re, text)
            root_energies = np.array(matches, dtype=float)
            root = self.get_root()
            energy = root_energies[root - 1]
        else:
            # All state average energies
            root_re = "RASSCF root number.+Total energy.+?" + self.float_regex
            matches = re.findall(root_re, text)
            root_energies = np.array(matches, dtype=float)

            # Energy of root for which gradient was computed
            energy_regex = r"RASSCF state energy =\s*" + self.float_regex
            try:
                energy = float(re.search(energy_regex, text).groups()[0])
            # Fallback to first energy, when no gradient was computed
            except AttributeError:
                energy = root_energies[0]
        return energy, root_energies

    def parse_energy(self, path):
        with open(path / self.out_fn) as handle:
            text = handle.read()
        energy, all_energies = self.parse_energies(text)
        return {
            "energy": energy,
            "all_energies": all_energies,
        }

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

        energy, _ = self.parse_energies(text)

        results["energy"] = energy
        # results["sa_energies"] = sa_energies
        results["forces"] = -gradient

        return results

    def parse_rassi_track(self, path):
        gradient_fn = path / self.out_fn
        with open(gradient_fn) as handle:
            text = handle.read()
        track_re = (
            r"Initial root:\s*(\d+)\s*Overlaps with current "
            r"states:(.+)New root:\s*(\d+)"
        )
        # overlap_re = "OVERLAP MATRIX FOR THE ORIGINAL STATES:(.+?)##"
        mobj = re.search(track_re, text, re.DOTALL)

        initial_root, overlaps, new_root = mobj.groups()
        overlaps = np.array(overlaps.strip().split(), dtype=float).reshape(-1, 2)
        # Filters for overlaps > 10% (0.1**2 ~ 0.31622)
        thresh = 0.1
        inds = np.where(np.abs(overlaps[:, 1]) > thresh**0.5)
        ov_perc_str = ", ".join([f"{nr:.0f}: {ov**2:.2%}" for nr, ov in overlaps[inds]])
        self.log(
            f"Overlaps between previous root {initial_root} and "
            f"new roots bigger {thresh:.0%}:  {ov_perc_str}. Will "
            f"use root {new_root} for the following gradient calculation."
        )
        if new_root != initial_root:
            self.log("Found a root flip!")
        self.mdrlxroot = new_root

    def get_chkfiles(self):
        return {
            "inporb": self.inporb,
        }

    def set_chkfiles(self, chkfiles):
        try:
            inporb = chkfiles["inporb"]
            self.inporb = inporb
            self.log(f"Set chkfile '{inporb}' as inporb.")
        except KeyError:
            self.log("Found no inporb information in chkfiles!")

    def __str__(self):
        return "OpenMolcas calculator"
