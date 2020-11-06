import re
import textwrap

import jinja2
import numpy as np

from pysisyphus.calculators.OverlapCalculator import OverlapCalculator
from pysisyphus.constants import BOHR2ANG


class DFTBp(OverlapCalculator):

    conf_key = "dftbp"
    max_ang_moms = {
        "mio-ext": {
            "H": "s",
            "O": "p",
        },
    }

    def __init__(self, parameter, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.parameter = parameter

        self.inp_fn = "dftb_in.hsd"
        self.parser_funcs = {
            "forces": self.parse_forces,
        }

        self.base_cmd = self.get_cmd("cmd")

        self.gen_geom_fn = "geometry.gen"
        self.slako_prefix = self.get_cmd("slako")
        self.dftb_tpl = jinja2.Template(
            textwrap.dedent(
                """
            Geometry = GenFormat {
              <<< "{{ gen_geom_fn }}"
            }
            Hamiltonian = DFTB {
              Scc = Yes
              SlaterKosterFiles = Type2FileNames {
                  Prefix = "{{ slako_prefix }}/{{ parameter }}/"
                  Separator = "-"
                  Suffix = ".skf"
                }
              MaxAngularMomentum {
               {%- for atom, ang_mom in max_ang_moms %}
                {{ atom }} = "{{ ang_mom }}"
               {%- endfor %}
              }
            }
            Analysis {
             {% for anal in analysis -%}
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
            textwrap.dedent("""
            {{ atom_num }} C
            {% for atom in atom_types %}{{ atom }} {% endfor %}
            {% for type_, x, y, z in types_xyz %}
             {{ loop.index }} {{type_}}  {{x}} {{y}} {{z}}
            {%- endfor %}
        """).strip())
        atom_num = len(atoms)
        unique_atoms = tuple(set(atoms))
        atom_types = {atom: i for i, atom in enumerate(unique_atoms, 1)}
        c3d = coords.reshape(-1, 3) * BOHR2ANG
        types_xyz = [
            (atom_types[atom], x, y, z) for atom, (x, y, z) in zip(atoms, c3d)
        ]
        gen_str = gen_fmt_tpl.render(
            atom_num=atom_num,
            atom_types=unique_atoms,
            types_xyz=types_xyz,
        )
        return gen_str

    def prepare_input(self, atoms, coords, calc_type):
        path = self.prepare_path(use_in_run=True)
        gen_str = self.get_gen_str(atoms, coords)
        with open(path / self.gen_geom_fn, "w") as handle:
            handle.write(gen_str)
        analysis = list()
        if calc_type == "forces":
            analysis.append("CalculateForces = Yes")
        ang_moms = self.max_ang_moms[self.parameter]
        max_ang_moms = [(atom, ang_moms[atom]) for atom in set(atoms)]

        inp = self.dftb_tpl.render(
              gen_geom_fn=self.gen_geom_fn,
              slako_prefix=self.slako_prefix,
              parameter=self.parameter,
              max_ang_moms=max_ang_moms,
              analysis=analysis,
        )
        return inp

    def get_forces(self, atoms, coords, prepare_kwargs=None):
        inp = self.prepare_input(atoms, coords, "forces")
        results = self.run(inp, "forces")
        return results

    def parse_total_energy(self, text):
        energy_re = re.compile("Total energy:\s*([-\d\.]+)\s*H")
        return float(energy_re.search(text)[1])

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
