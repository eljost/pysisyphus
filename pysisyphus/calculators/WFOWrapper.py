from pyscf import gto

from pysisyphus.config import Config

class WFOWrapper:
    def __init__(self):
        self.base_cmd = Config["wfoverlap"]["cmd"]

    def build_atom(atoms, coords):

    def overlap(atoms, coords1, coords2, basis, mos1, mos2):
        atom1 = [
            ["O", (0.23316300,     0.26027600,     0.49195600)],
            ["H", (1.11389400,     1.55866300,    -0.47754700)],
            ["H", (1.47215800,    -1.10158000,     0.60031200)],
        ]
        atom2 = [
            ["O", (0.24316300,     0.27027600,     0.48195600)],
            ["H", (1.10389400,     1.56866300,    -0.46754700)],
            ["H", (1.48215800,    -1.11158000,     0.61031200)],
        ]

        def prepare(atom):
            mol = gto.Mole()
            mol.atom = atom
            mol.basis = "sto-3g"
            mol.charge = 0
            mol.spin = 0
            mol.build()
            return mol
        mol1 = prepare(atom1)
        mol2 = prepare(atom2)
        ao_ovlp = gto.mole.intor_cross("int1e_ovlp_sph", mol1, mol2)
        print(ao_ovlp)
        ref = np.loadtxt("mix_ovl.new")
        print()
        print(ref)
