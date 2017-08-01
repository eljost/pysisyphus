import numpy as np

from pysisyphus.calculators.Calculator import Calculator

class AnaPot(Calculator):

    def __init__(self): 
        super(AnaPot, self).__init__()

    def get_energy(self, atoms, coords):
        x, y, z = coords
        energy = 4 + 4.5*x - 4*y + x**2 + 2*y**2-2*x*y + x**4 - 2*x**2*y
        return {"energy": energy}

    def get_forces(self, atoms, coords):
        x, y, z = coords
        forces = (4.5 + 2*x -2*y + 4*x**3 - 4*x*y,
                 (-4 + 4*y - 2*x - 2*x**2),
                  0
        )
        forces = -np.array(forces)
        results = self.get_energy(atoms, coords)
        results["forces"] = forces
        return results

    """
    def get_hessian(self, atoms, coords):
        x, y, z = coords
        self._hessian = ((12*x**2 + 2 - 4*y, -4*x-2),
                         (-4*x-2, 4)
        )
    """

    def __str__(self):
        return "AnaPot calculator"
