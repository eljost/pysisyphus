import numpy as np

try:
    from openbabel import openbabel as ob
    from openbabel import pybel
except ModuleNotFoundError:
    pass

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG, AU2KJPERMOL, AU2KCALMOL


def _getpluginnames(ptype):
    plugins = ob.vectorString()
    ob.OBPlugin.ListAsVector(ptype, None, plugins)
    return [p.split()[0].lower() for p in plugins if p.strip()]


def _getplugins(findplugin, names):
    return dict([(x, findplugin(x)) for x in names if findplugin(x)])


class OBabel(Calculator):
    conv_dict = {
        "kj/mol": AU2KJPERMOL,
        "kcal/mol": AU2KCALMOL,
    }

    def __init__(self, ff="gaff", mol=None, **kwargs):
        super().__init__(**kwargs)

        if self.charge != 0:
            print("'charge' keyword is currently ignored!")

        self.ff_name = ff
        avail_ffs = _getplugins(
            ob.OBForceField.FindType, _getpluginnames("forcefields")
        )
        self.ff = avail_ffs[self.ff_name]
        self.mol = mol

        unit = self.ff.GetUnit()
        an_grad = self.ff.HasAnalyticalGradients()
        self.log(f"Energy unit: {unit}, analytical gradient: {an_grad}")
        # Used to convert energy from force field units to atomic units
        self.conv_fac = self.conv_dict[unit.lower()]

        self.is_setup = False

    def setup(self, atoms, coords):
        xyz = self.prepare_xyz_string(atoms, coords)
        if self.mol is None:
            self.mol = pybel.readstring("xyz", xyz)
        else:
            self.mol.OBMol.SetCoordinates(ob.double_array(coords*BOHR2ANG))

        if not self.is_setup:
            self.is_setup = self.ff.Setup(self.mol.OBMol)
            assert self.is_setup
        else:
            self.ff.SetCoordinates(self.mol.OBMol)
        return self.mol

    def get_energy(self, atoms, coords, **prepare_kwargs):
        self.setup(atoms, coords)
        energy = self.ff.Energy() / self.conv_fac

        results = {"energy": energy}
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        mol = self.setup(atoms, coords)
        ff = self.ff

        energy = ff.Energy(True) / self.conv_fac
        grad = [ff.GetGradient(atom.OBAtom) for atom in mol.atoms]
        grad = np.array([(g.GetX(), g.GetY(), g.GetZ()) for g in grad])
        # WTH?`It seems openbabel does not return the gradient, but the force...
        # So we don't multiply with -1 here...
        forces = grad.flatten() * BOHR2ANG / self.conv_fac

        results = {"energy": energy, "forces": forces}
        return results
