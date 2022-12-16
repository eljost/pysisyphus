import numpy as np

try:
    from pyscf import gto, scf, df
    from pyscf.dft import numint
    from pyscf.tools import cubegen
except ModuleNotFoundError:
    pass
import pytest

from pysisyphus.wavefunction import Shells
from pysisyphus.testing import using


@using("pyscf")
def test_pyscf_shells():
    mol = gto.M(atom="H 0 0 0; H 0 0 1", basis="def2qzvpp", parse_arg=False)
    shells = Shells.from_pyscf_mol(mol)
    S_ref = mol.intor("int1e_ovlp_sph")
    S = shells.S_sph
    np.set_printoptions(suppress=True, precision=4, linewidth=240)
    np.testing.assert_allclose(S, S_ref, atol=1e-12)


@using("pyscf")
def test_grid_density():
    mol = gto.M(atom="H 0 0 0; H 0 0 1", basis="def2qzvpp", parse_arg=False)
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

    mol = gto.M(verbose=0, atom="H 0 0 0; H 0 0 1.0;", basis="ccpvtz", parse_arg=False)

    nao = mol.nao

    origin = _charge_center(mol)
    with mol.with_common_orig(origin):
        quad = mol.intor("int1e_rr").reshape(3, 3, nao, nao)
    shells = Shells.from_pyscf_mol(mol)
    pysis_quad = shells.get_quadrupole_ints_sph(origin)
    np.testing.assert_allclose(pysis_quad, quad, atol=1e-14)

    # Only diagonal elements
    pysis_dquad = shells.get_diag_quadrupole_ints_sph(origin)
    np.testing.assert_allclose(
        pysis_dquad, pysis_quad[[0, 1, 2], [0, 1, 2]], atol=1e-14
    )


@using("pyscf")
@pytest.fixture
def mol_auxmol():

    mol = gto.Mole()
    mol.atom = "He 0 0 0"
    basis = {
        "He": gto.basis.parse(
            """
    He    S
         1.0 1.0
         2.0 1.0
    He    P
         1.0 1.0
         2.0 1.0
    """
        )
    }
    mol.basis = basis
    mol.build(parse_arg=False)
    auxbasis = {
        "He": gto.basis.parse(
            """
    He    S
         1.0 1.0
         4.0 1.0
         5.0 1.0
    He    D
         4.0 1.0
         5.0 1.0
    """
        )
    }
    auxmol = df.make_auxmol(mol, auxbasis)
    return (mol, auxmol)


def test_2c2e(mol_auxmol):
    _, auxmol = mol_auxmol
    int2c = auxmol.intor("int2c2e_cart")
    # Cartesian d-Orbitals are not normalized in pyscf
    S_aux = auxmol.intor("int1e_ovlp_cart")
    N_aux = 1 / np.diag(S_aux) ** 0.5
    NaNa = N_aux[:, None] * N_aux[None, :]

    aux_shells = Shells.from_pyscf_mol(auxmol)
    integrals = aux_shells.get_2c2el_ints_cart()

    np.testing.assert_allclose(integrals, int2c * NaNa)


def test_3c2e(mol_auxmol):
    mol, auxmol = mol_auxmol
    int3c = df.incore.aux_e2(mol, auxmol, "int3c2e_sph")

    shells = Shells.from_pyscf_mol(mol)
    aux_shells = Shells.from_pyscf_mol(auxmol)
    integrals = shells.get_3c2el_ints_sph(aux_shells)

    np.testing.assert_allclose(integrals, int3c, atol=1e-12)
