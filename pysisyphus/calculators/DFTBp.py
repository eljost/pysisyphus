import itertools as it
from math import ceil
from pathlib import Path
import re
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
        xpy = list()
        for line in block[1:]:
            xpy.extend(line.split())
        xpy = np.array(xpy, dtype=float)
        xpys.append(xpy)
    return size, states, np.array(xpys)


class DFTBp(OverlapCalculator):
    conf_key = "dftbp"
    _set_plans = (
        "out",
        "exc_dat",
        "eigenvec",
        "xplusy_dat",
    )

    max_ang_moms = {
        "mio-ext": {
            "H": "s",
            "C": "p",
            "N": "p",
            "O": "p",
        },
        "3ob": {
            "Br": "d",
            "C": "p",
            "Ca": "p",
            "Cl": "d",
            "F": "p",
            "H": "s",
            "I": "d",
            "K": "p",
            "Mg": "p",
            "N": "p",
            "Na": "p",
            "O": "p",
            "P": "d",
            "S": "d",
            "Zn": "d",
        },
    }
    hubbard_derivs = {
        "3ob": {
            "Br": -0.0573,
            "Mg": -0.0200,
            "C": -0.1492,
            "N": -0.1535,
            "Ca": -0.0340,
            "Na": -0.0454,
            "Cl": -0.0697,
            "O": -0.1575,
            "F": -0.1623,
            "P": -0.1400,
            "H": -0.1857,
            "S": -0.1100,
            "I": -0.0433,
            "Zn": -0.0300,
            "K": -0.0339,
        },
    }

    def __init__(self, parameter, *args, slakos=None, **kwargs):
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
        self.base_cmd = self.get_cmd()
        self.gen_geom_fn = "geometry.gen"
        self.inp_fn = "dftb_in.hsd"
        self.out_fn = "dftb.out"
        self.to_keep = (
            "dftb_in.hsd",
            "out:detailed.out",
            self.gen_geom_fn,
            self.out_fn,
            "xplusy_dat:XplusY.DAT",
            "exc_dat:EXC.DAT",
            "eigenvec:eigenvec.out",
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
              {% if parameter == "3ob" %}
                ThirdOrderFull = Yes
                HubbardDerivs {
                {% for atom, hd in hubbard_derivs -%}
                 {{ atom }} = {{ hd }}
                {% endfor %}
                }
                HCorrection = Damping {
                 Exponent = 4.00
                }
              {% endif %}
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
    def get_excited_state_str(track, root, nroots, forces=False):
        if root is None and (track == False):
            return ""

        casida_tpl = jinja2.Template(
            textwrap.dedent(
                """
            Casida {
                NrOfExcitations = {{ nstates }}
                Symmetry = Singlet
                {% if root %}StateOfInterest = {{ root }}{% endif %}
                WriteXplusY = Yes
                {{ es_forces }}
            }
            """
            )
        )
        es_forces = "ExcitedStateForces = Yes" if forces else ""
        es_str = casida_tpl.render(
            nstates=nroots if nroots else root + 5,
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
        if self.track or self.root:
            analysis.extend(("WriteEigenvectors = Yes", "EigenvectorsAsText = Yes"))
        ang_moms = self.max_ang_moms[self.parameter]
        unique_atoms = set(atoms)
        max_ang_moms = [(atom, ang_moms[atom]) for atom in unique_atoms]

        # spinpol = (
        # "{}"
        # if (self.mult == 1)
        # else f"Colinear {{ UnpairedElectrons = {self.mult-1} }}"
        # )
        spinpol = "{}"
        es_forces = calc_type == "forces"
        try:
            param_hubbard_derivs = self.hubbard_derivs[self.parameter]
            hubbard_derivs = [
                (atom, param_hubbard_derivs[atom]) for atom in unique_atoms
            ]
        except KeyError:
            hubbard_derivs = list()

        inp = self.dftb_tpl.render(
            gen_geom_fn=self.gen_geom_fn,
            charge=self.charge,
            spinpol=spinpol,
            slakos_prefix=self.slakos_prefix,
            parameter=self.parameter,
            max_ang_moms=max_ang_moms,
            hubbard_derivs=hubbard_derivs,
            excited_state_str=self.get_excited_state_str(
                self.track, self.root, self.nroots, es_forces
            ),
            analysis=analysis,
        )
        return inp, path

    def store_and_track(self, results, func, atoms, coords, **prepare_kwargs):
        if self.track:
            self.store_overlap_data(atoms, coords)
            if self.track_root():
                # Redo the calculation with the updated root
                results = func(atoms, coords, **prepare_kwargs)
        results["all_energies"] = self.parse_all_energies()
        return results

    def get_energy(self, atoms, coords, **prepare_kwargs):
        inp, path = self.prepare_input(atoms, coords, "energy")
        results = self.run(inp, "energy")
        results = self.store_and_track(
            results, self.get_energy, atoms, coords, **prepare_kwargs
        )
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        inp, path = self.prepare_input(atoms, coords, "forces")
        run_kwargs = {
            "calc": "forces",
            "hold": self.track,
        }
        results = self.run(inp, **run_kwargs)
        results = self.store_and_track(
            results, self.get_forces, atoms, coords, **prepare_kwargs
        )
        # if self.track:
        # self.calc_counter += 1
        # self.store_overlap_data(atoms, coords, path)
        # if self.track_root():
        # # Redo the calculation with the updated root
        # results = self.get_forces(atoms, coords)
        # try:
        # shutil.rmtree(path)
        # except FileNotFoundError:
        # self.log(f"'{path}' has already been deleted!")
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

    def parse_all_energies(self, out_fn=None, exc_dat=None):
        if out_fn is None:
            out_fn = self.out
        if exc_dat is None:
            try:
                exc_dat = self.exc_dat
            except AttributeError:
                exc_dat = None

        with open(out_fn) as handle:
            detailed = handle.read()

        gs_energy = self.parse_total_energy(detailed)
        if (exc_dat is not None) and exc_dat.exists():
            with open(exc_dat) as handle:
                exc_dat = handle.read()
            exc_ens = self.parse_exc_dat(exc_dat)
            all_energies = np.full(len(exc_ens) + 1, gs_energy)
            all_energies[1:] += exc_ens
        else:
            all_energies = np.array((gs_energy,))
        return all_energies

    def prepare_overlap_data(self, path):
        #
        # Excitation energies
        #
        # import pdb; pdb.set_trace()  # fmt: skip
        # with open(path / "detailed.out") as handle:
        with open(self.out) as handle:
            detailed = handle.read()
        all_energies = self.parse_all_energies(path)

        #
        # MO coefficients
        #
        # with open(path / "eigenvec.out") as handle:
        with open(self.eigenvec) as handle:
            eigenvecs = handle.read().strip()
        eigenvecs = eigenvecs.split("Eigenvector")[1:]

        C = np.array([parse_mo(eigvec) for eigvec in eigenvecs], dtype=float).T
        assert C.shape[0] == C.shape[1]

        #
        # CI coefficients
        #
        electron_re = re.compile(r"Nr. of electrons \(up\):\s*([\d\.]+)")
        mobj = electron_re.search(detailed)
        electrons = int(float(mobj[1]))
        assert electrons % 2 == 0
        occ = electrons // 2
        mo_num = C.shape[0]
        vir = mo_num - occ
        # X+Y
        with open(self.xplusy_dat) as handle:
            xpy_text = handle.read().strip()
        size, states, xpy = parse_xplusy(xpy_text)
        # As of DFTB+ 22.2, we can only get X+Y, but not X-Y, so we can't reconstruct
        # the true X and Y vectors.
        XpY = xpy.reshape(states, occ, vir)
        X = XpY
        Y = np.zeros_like(XpY)
        assert size == occ * vir
        return C, X, Y, all_energies
