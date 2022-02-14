from math import ceil
from pathlib import Path
import re
import shutil
import textwrap

import jinja2
import numpy as np

from pysisyphus.calculators.OverlapCalculator import OverlapCalculator
from pysisyphus.constants import BOHR2ANG, AU2EV


def parse_mo(eigvec):
    coeffs = list()
    # Skip first line
    for line in eigvec.split("\n")[1:]:
        line = line.strip()
        if line == "":
            continue
        line = line.split()
        if len(line) == 5:
            line = line[-3:]
        coeffs.append(line[1])
    return coeffs


def parse_xplusy(text):
    lines = text.split("\n")
    size, states = [int(i) for i in lines.pop(0).split()]
    # Add one for header line
    block_size = ceil(size / 6) + 1
    assert len(lines) == (states * block_size)

    xpys = list()
    for i in range(states):
        block = lines[i * block_size : (i + 1) * block_size]
        _, *rest = block
        xpy = np.array([line.split() for line in rest], dtype=float)
        xpys.append(xpy)
    return size, states, np.array(xpys)


class DFTBp(OverlapCalculator):

    conf_key = "dftbp"
    max_ang_moms = {
        "mio-ext": {
            "H": "s",
            "C": "p",
            "N": "p",
            "O": "p",
        },
    }

    def __init__(self, parameter, *args, slakos=None, root=None, **kwargs):
        super().__init__(*args, **kwargs)

        assert self.mult == 1, "Open-shell not yet supported!"
        self.parameter = parameter
        if slakos is None:
            slakos = self.get_cmd("slakos")
        self.slakos_prefix = str(slakos)
        assert (Path(self.slakos_prefix) / self.parameter).exists(), (
            f"Expected '{self.parameter}' sub-directory in '{self.slakos_prefix}' "
            "but could not find it!"
        )
        self.root = root

        self.base_cmd = self.get_cmd()
        self.gen_geom_fn = "geometry.gen"
        self.inp_fn = "dftb_in.hsd"
        self.out_fn = "dftb.out"
        self.to_keep = (
            "dftb_in.hsd",
            "detailed.out",
            self.gen_geom_fn,
            self.out_fn,
            "XplusY.DAT",
            "EXC.DAT",
        )
        self.parser_funcs = {
            "energy": self.parse_energy,
            "forces": self.parse_forces,
        }

        self.dftb_tpl = jinja2.Template(
            textwrap.dedent(
                """
            Geometry = GenFormat {
              <<< "{{ gen_geom_fn }}"
            }
            Hamiltonian = DFTB {
              Scc = Yes
              Charge = {{ charge }}
              SpinPolarisation = {{ spinpol }}
              SlaterKosterFiles = Type2FileNames {
                  Prefix = "{{ slakos_prefix }}/{{ parameter }}/"
                  Separator = "-"
                  Suffix = ".skf"
                }
              MaxAngularMomentum {
               {%- for atom, ang_mom in max_ang_moms %}
                {{ atom }} = "{{ ang_mom }}"
               {%- endfor %}
              }
            }

            ExcitedState {
             {{ excited_state_str }}
            }

            Analysis {
             {% for anal in analysis %}
              {{ anal }}
             {%- endfor %}
            }

            ParserOptions {
              ParserVersion = 8
            }
        """
            ).strip()
        )

    @staticmethod
    def get_gen_str(atoms, coords):
        gen_fmt_tpl = jinja2.Template(
            textwrap.dedent(
                """
            {{ atom_num }} C
            {% for atom in atom_types %}{{ atom }} {% endfor %}
            {% for type_, x, y, z in types_xyz %}
             {{ loop.index }} {{type_}}  {{x}} {{y}} {{z}}
            {%- endfor %}
        """
            ).strip()
        )
        atom_num = len(atoms)
        unique_atoms = tuple(set(atoms))
        atom_types = {atom: i for i, atom in enumerate(unique_atoms, 1)}
        c3d = coords.reshape(-1, 3) * BOHR2ANG
        types_xyz = [(atom_types[atom], x, y, z) for atom, (x, y, z) in zip(atoms, c3d)]
        gen_str = gen_fmt_tpl.render(
            atom_num=atom_num,
            atom_types=unique_atoms,
            types_xyz=types_xyz,
        )
        return gen_str

    @staticmethod
    def get_excited_state_str(root, forces=False):
        if root is None:
            return ""

        casida_tpl = jinja2.Template(
            textwrap.dedent(
                """
            Casida {
                NrOfExcitations = {{ nstates }}
                Symmetry = Singlet
                StateOfInterest = {{ root }}
                WriteXplusY = Yes
                {{ es_forces }}
            }
            """
            )
        )
        es_forces = "ExcitedStateForces = Yes" if forces else ""
        es_str = casida_tpl.render(
            nstates=root + 5,
            root=root,
            es_forces=es_forces,
        )
        return es_str

    def prepare_input(self, atoms, coords, calc_type):
        path = self.prepare_path(use_in_run=True)
        gen_str = self.get_gen_str(atoms, coords)
        with open(path / self.gen_geom_fn, "w") as handle:
            handle.write(gen_str)
        analysis = list()
        if calc_type == "forces":
            analysis.append("CalculateForces = Yes")
        if self.root:
            analysis.extend(("WriteEigenvectors = Yes", "EigenvectorsAsText = Yes"))
        ang_moms = self.max_ang_moms[self.parameter]
        max_ang_moms = [(atom, ang_moms[atom]) for atom in set(atoms)]

        # spinpol = (
        # "{}"
        # if (self.mult == 1)
        # else f"Colinear {{ UnpairedElectrons = {self.mult-1} }}"
        # )
        spinpol = "{}"
        es_forces = calc_type == "forces"

        inp = self.dftb_tpl.render(
            gen_geom_fn=self.gen_geom_fn,
            charge=self.charge,
            spinpol=spinpol,
            slakos_prefix=self.slakos_prefix,
            parameter=self.parameter,
            max_ang_moms=max_ang_moms,
            excited_state_str=self.get_excited_state_str(self.root, es_forces),
            analysis=analysis,
        )
        return inp, path

    def get_energy(self, atoms, coords, **prepare_kwargs):
        inp, path = self.prepare_input(atoms, coords, "energy")
        results = self.run(inp, "energy")
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        inp, path = self.prepare_input(atoms, coords, "forces")
        run_kwargs = {
            "calc": "forces",
            "hold": self.track,
        }
        results = self.run(inp, **run_kwargs)
        if self.track:
            self.calc_counter += 1
            self.store_overlap_data(atoms, coords, path)
            if self.track_root():
                # Redo the calculation with the updated root
                results = self.get_forces(atoms, coords)
        try:
            shutil.rmtree(path)
        except FileNotFoundError:
            self.log(f"'{path}' has already been deleted!")
        return results

    def run_calculation(self, atoms, coords, **prepare_kwargs):
        inp, path = self.prepare_input(atoms, coords, "energy")
        run_kwargs = {
            "calc": "energy",
            "hold": self.track,
        }
        results = self.run(inp, **run_kwargs)
        if self.track:
            self.calc_counter += 1
            self.store_overlap_data(atoms, coords, path)
        return results

    def parse_total_energy(self, text):
        energy_re = re.compile(r"Total energy:\s*([-\d\.]+)\s*H")
        exc_energy_re = re.compile(r"Excitation Energy:\s*([\-\.\d]+)\s*H")
        energy = float(energy_re.search(text)[1])
        exc_mobj = exc_energy_re.search(text)
        if exc_mobj:
            energy += float(exc_mobj[1])
        return energy

    def parse_energy(self, path):
        detailed = path / "detailed.out"
        with open(detailed) as handle:
            text = handle.read()
        results = {
            "energy": self.parse_total_energy(text),
        }
        return results

    @staticmethod
    def parse_exc_dat(text):
        exc_re = re.compile("=+(.+)", re.DOTALL)
        mobj = exc_re.search(text)
        exc_lines = mobj[1].strip().split("\n")
        exc_ens = (
            np.array([line.strip().split()[0] for line in exc_lines], dtype=float)
            / AU2EV
        )
        return exc_ens

    def parse_forces(self, path):
        forces_re = re.compile("Total Forces(.+)Maximal derivative", re.DOTALL)
        detailed = path / "detailed.out"
        with open(detailed) as handle:
            text = handle.read()
        mobj = forces_re.search(text)
        forces = np.array(mobj[1].strip().split(), dtype=float).reshape(-1, 4)[:, 1:]
        results = {
            "energy": self.parse_total_energy(text),
            "forces": forces.flatten(),
        }
        return results

    def prepare_overlap_data(self, path):
        #
        # Excitation energies
        #
        with open(path / "detailed.out") as handle:
            detailed = handle.read()
        gs_energy = self.parse_total_energy(detailed)
        with open(path / "EXC.DAT") as handle:
            exc_dat = handle.read()
        exc_ens = self.parse_exc_dat(exc_dat)
        all_energies = np.full(len(exc_ens) + 1, gs_energy)
        all_energies[1:] += exc_ens

        #
        # MO coefficients
        #
        with open(path / "eigenvec.out") as handle:
            eigenvecs = handle.read().strip()
        eigenvecs = eigenvecs.split("Eigenvector")[1:]

        mo_coeffs = np.array([parse_mo(eigvec) for eigvec in eigenvecs], dtype=float)
        assert mo_coeffs.shape[0] == mo_coeffs.shape[1]

        #
        # CI coefficients
        #
        electron_re = re.compile(r"Nr. of electrons \(up\):\s*([\d\.]+)")
        mobj = electron_re.search(detailed)
        electrons = int(float(mobj[1]))
        assert electrons % 2 == 0
        occ = electrons // 2
        mo_num = mo_coeffs.shape[0]
        vir = mo_num - occ
        # X+Y
        with open(path / "XplusY.DAT") as handle:
            xpy_text = handle.read().strip()
        size, states, xpy = parse_xplusy(xpy_text)
        ci_coeffs = xpy.reshape(states, occ, vir)
        assert size == occ * vir
        return mo_coeffs, ci_coeffs, all_energies
