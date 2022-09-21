import numpy as np

try:
    from pyscf import gto, scf
    from pyscf.dft import numint
    from pyscf.tools import cubegen
except ModuleNotFoundError:
    pass

from pysisyphus.wavefunction import Shells
from pysisyphus.testing import using


@using("pyscf")
def test_pyscf_shells():
    mol = gto.M(atom="H 0 0 0; H 0 0 1", basis="def2qzvpp")
    shells = Shells.from_pyscf_mol(mol)
    S_ref = mol.intor("int1e_ovlp_sph")
    S = shells.S_sph
    np.set_printoptions(suppress=True, precision=4, linewidth=240)
    np.testing.assert_allclose(S, S_ref, atol=1e-12)


@using("pyscf")
def test_grid_density():
    mol = gto.M(atom="H 0 0 0; H 0 0 1", basis="def2qzvpp")
    mf = scf.RHF(mol)
    mf.scf()
    dm = mf.make_rdm1()
    nx = ny = nz = 8
    cc = cubegen.Cube(mol, nx, ny, nz, resolution=None, margin=3.0)

    GTOval = "GTOval"
    coords3d = cc.get_coords()
    ao = mol.eval_gto(GTOval, coords3d)
    rho_ref = numint.eval_rho(mol, ao, dm)

    shells = Shells.from_pyscf_mol(mol)
    vals = shells.eval(coords3d, spherical=True)
    np.testing.assert_allclose(vals, ao, atol=1e-10)
    rho = np.einsum("uv,iu,iv->i", dm, vals, vals)
    np.testing.assert_allclose(rho, rho_ref)


@using("pyscf")
def test_quadrupole_ints():
    def _charge_center(mol):
        charges = mol.atom_charges()
        coords = mol.atom_coords()
        return np.einsum("z,zr->r", charges, coords) / charges.sum()

    mol = gto.M(verbose=0, atom="H 0 0 0; H 0 0 1.0;", basis="ccpvtz")

    nao = mol.nao

    origin = _charge_center(mol)
    with mol.with_common_orig(origin):
        quad = mol.intor("int1e_rr").reshape(3, 3, nao, nao)
    shells = Shells.from_pyscf_mol(mol)
    pysis_quad = shells.get_quadrupole_ints_sph(origin)
    np.testing.assert_allclose(pysis_quad, quad, atol=1e-14)
