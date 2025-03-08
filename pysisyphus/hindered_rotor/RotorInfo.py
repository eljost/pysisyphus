import dataclasses
import textwrap
from typing import Literal, Optional

import jinja2
import numpy as np

from pysisyphus.constants import AMU2AU, AU2EV
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.hindered_rotor import fragment as hr_fragment, inertmom
from pysisyphus import xyzloader


@dataclasses.dataclass
class Fragment:
    key: Literal["left", "right"]
    natoms: int
    tot_mass: float
    indices: list[int]
    imom: float

    def __post_init__(self):
        self.tot_mass_au = self.tot_mass * AMU2AU
        self.B_const = 1 / (2.0 * self.imom)  # in au
        self.B_const_mueV = self.B_const * AU2EV * 1e6  # in μeV


ROTOR_REPORT_TPL = jinja2.Template(
    textwrap.dedent(
        """
    {{ title }}

    Whole Geometry
    --------------
                Atoms: {{ natoms }}
           Total Mass: {{ "{0:.6f}".format(tot_mass) }} amu ({{ tot_mass_au }} au)
      Torsion indices: {{ indices }}
    {% for frag in fragments %}
    {{ frag.key.capitalize() }} Fragment
    ---------------
                Atoms: {{ frag.natoms }}
        Fragment Mass: {{ "{0:.6f}".format(frag.tot_mass) }} amu ({{ "{0:.6f}".format(frag.tot_mass_au) }} au)
              Indices: {{ frag.indices }}
    Moment of inertia: {{ "{0:.6}".format(frag.imom) }} au
      Rot. constant B: {{ "{0:.2f}".format(frag.B_const_mueV) }} μeV
    {% endfor %}
    Moments of interia determined by I(m={{m}}, n={{n}}) approximation. """
    )
)


@dataclasses.dataclass
class RotorInfo:
    atoms: tuple[str, ...]
    coords3d: np.ndarray
    # Masses in amu
    masses: np.ndarray
    indices_left: list[int]
    indices_right: list[int]
    indices: list[int]
    bond: list[int]
    # Moments of inertia in atomic units NOT in atomic mass units!
    imom_left: float
    imom_right: float
    m: int
    n: int

    @staticmethod
    def from_torsion(geom, torsion_indices: list[int], **kwargs):
        return prepare_rotor_info(geom, torsion_indices, **kwargs)

    @property
    def natoms(self):
        return len(self.atoms)

    @property
    def tot_mass(self):
        return self.masses.sum()

    def get_fragment(self, key: Literal["left", "right"]):
        if key == "left":
            indices = self.indices_left
            imom = self.imom_left
        elif key == "right":
            indices = self.indices_right
            imom = self.imom_right
        tot_mass = self.masses[indices].sum()
        return Fragment(
            key=key,
            natoms=len(indices),
            tot_mass=tot_mass,
            indices=indices,
            imom=imom,
        )

    @property
    def fragment_left(self):
        return self.get_fragment("left")

    @property
    def fragment_right(self):
        return self.get_fragment("right")

    def render_report(self) -> str:
        title = highlight_text("Rotor info")
        tot_mass = self.tot_mass
        fragments = (self.fragment_left, self.fragment_right)
        rendered = ROTOR_REPORT_TPL.render(
            title=title,
            natoms=self.natoms,
            tot_mass=tot_mass,
            tot_mass_au=tot_mass * AMU2AU,
            indices=self.indices,
            fragments=fragments,
            m=self.m,
            n=self.n,
        )
        return rendered

    def as_xyzs(self):
        c3d_left = self.coords3d.copy()[self.indices_left]
        c3d_right = self.coords3d.copy()[self.indices_right]
        c3d_blocked = np.concatenate((c3d_left, c3d_right), axis=0)

        atoms_left = [self.atoms[i] for i in self.indices_left]
        atoms_left_rest = ["x"] * len(self.indices_right)
        atoms_left = atoms_left + atoms_left_rest
        xyz_left = xyzloader.make_xyz_str_au(
            atoms_left, c3d_blocked, comment="Fragment left"
        )

        atoms_right = [self.atoms[i] for i in self.indices_right]
        atoms_right_rest = ["x"] * len(self.indices_left)
        atoms_right = atoms_right_rest + atoms_right
        xyz_right = xyzloader.make_xyz_str_au(
            atoms_right, c3d_blocked, comment="Fragment right"
        )
        return xyz_left, xyz_right

    def dump_trj(self, fn):
        xyzs = self.as_xyzs()
        trj = "\n".join(xyzs)
        with open(fn, "w") as handle:
            handle.write(trj)


def prepare_rotor_info(
    geom: Geometry,
    torsion_indices: list[int],
    rotor_indices: Optional[list[int]] = None,
    m: Literal[1, 2, 3] = 2,
    n: Literal[1, 2, 3] = 3,
) -> RotorInfo:
    bond = torsion_indices[1:3]
    # If no indices for the left fragment were given we try to determine them.
    if rotor_indices is None:
        print("Using automatic fragmentation to determine rotor indices.")
        indices_left, indices_right = hr_fragment.fragment_geom(geom, bond)
    # ... otherwise we put all atoms not present in the left fragment into the right
    # fragment.
    else:
        print("Using user supplied indices for rotor.")
        natoms = len(geom.atoms)
        indices_left = rotor_indices
        indices_right = [i for i in range(natoms) if i not in rotor_indices]

    imom_left, imom_right = inertmom.get_top_moment_of_inertia(
        geom.coords3d, geom.masses, indices_left, bond, m=m, n=n
    )
    imom_left *= AMU2AU
    imom_right *= AMU2AU
    rot_info = RotorInfo(
        atoms=geom.atoms,
        coords3d=geom.coords3d.copy(),
        masses=geom.masses.copy(),
        indices_left=indices_left,
        indices_right=indices_right,
        indices=torsion_indices,
        bond=bond,
        imom_left=imom_left,
        imom_right=imom_right,
        m=m,
        n=n,
    )
    return rot_info
