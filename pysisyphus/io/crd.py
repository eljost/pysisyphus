import re

import jinja2
import numpy as np

from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import file_or_str


@file_or_str(".crd")
def parse_crd(crd):
    lines = [l.strip() for l in crd.strip().split("\n") if not l.startswith("*")]
    expect_line = lines.pop(0)
    expect_num = int(expect_line.split()[0])
    coords3d = np.zeros((expect_num, 3), dtype=float)
    atoms = list()
    re_ = re.compile(r"\d+")
    for i, line in enumerate(lines):
        _, _, res, atom, x, y, z, *_ = line.split()
        coords3d[i] = (x, y, z)
        atom = re_.sub("", atom)  # Delte number
        atoms.append(atom)
    assert len(atoms) == expect_num
    coords3d *= ANG2BOHR
    return atoms, coords3d.flatten()


def geom_from_crd(fn, **kwargs):
    atoms, coords = parse_crd(fn)
    return Geometry(atoms, coords, **kwargs)


# I10   I10 2X A8 2X A8       3F20.10     2X A8 2X A8 F20.10
CRD_TPL = """* Created by pysisyphus
*
{{ i10(atomnum) }}  EXT
{% for atom, x, y, z in atoms_xyz -%}
 {{ i10(loop.index) }}{{ i10(resno) }}  {{ a8(res) }}  {{ a8(atom) }}{{ f20(x) }}{{ f20(y) }}{{ f20(z) }}  {{ a8(segid) }}  {{ a8(resid) }}{{ f20(0.0) }}
{% endfor %}
"""


def i10(int_):
    return f"{int_:10d}"


def f20(float_):
    return f"{float_:>20.10f}"


def a8(str_):
    return f"{str(str_):<8s}"


def atoms_coords_to_crd_str(
    atoms,
    coords,
    resno=1,
    res="UNL1",
    segid=None,
    resid=1,
    ref_atoms=None,
    del_atoms=None,
):
    if segid is None:
        segid = res
    if del_atoms is None:
        del_atoms = []

    env = jinja2.Environment(loader=jinja2.BaseLoader)
    env.globals.update(i10=i10, f20=f20, a8=a8)
    tpl = env.from_string(CRD_TPL)

    if ref_atoms is not None:
        atoms_for_names = ref_atoms
    else:
        atoms_for_names = atoms
    counter = dict()
    atom_names = list()
    for atom in atoms_for_names:
        counter.setdefault(atom, 1)
        atom_names.append(f"{atom}{counter[atom]}")
        counter[atom] += 1
    atom_names = [an for i, an in enumerate(atom_names) if i not in del_atoms]

    coords3d = coords.reshape(-1, 3) / ANG2BOHR
    atoms_xyz = list(zip(atom_names, *coords3d.T))
    rendered = tpl.render(
        atomnum=len(atoms_xyz),
        atoms_xyz=atoms_xyz,
        resno=resno,
        res=res,
        resid=resid,
        segid=segid,
    )
    return rendered


def geom_to_crd_str(geom, **kwargs):
    return atoms_coords_to_crd_str(geom.atoms, geom.cart_coords, **kwargs)
