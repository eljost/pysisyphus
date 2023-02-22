import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.wavefunction.Basis import (
    basis_from_orca_str,
    basis_from_pyscf_str,
    shells_with_basis,
)


def test_orca_bas():
    basis_str = """
    %basis
     newgto 1
     s 2
     1 0.7 0.5
     1 0.8 0.6
     p 2
     1 0.7 0.5
     1 0.8 0.6
     d 2
     1 0.7 0.5
     1 0.8 0.6
     d 3
     1 0.7 0.5
     1 0.8 0.6
     1 0.6 0.4
     end
    end
    """
    bas = basis_from_orca_str(basis_str)
    shells = bas[1]["electron_shells"]
    assert len(shells) == 4


def test_pyscf_bas():
    basis_str = """
    # Comment1
    He    S
         13.6267000              0.1752300
          1.9993500              0.8934830
          0.3829930              0.0000000
    He    S
         13.6267000              0.0000000
          1.9993500              0.0000000
          0.3829930              1.0000000
    """
    bas = basis_from_pyscf_str(basis_str)
    shells = bas[2]["electron_shells"]
    assert len(shells) == 2


@pytest.mark.parametrize(
    "fn, ref_len",
    (
        ("631G_CH.json", 13),  # SP-shells
        ("ccpvdz_CH.json", 18),  # General contraction
    ),
)
def test_json_basis(fn, ref_len, this_dir):
    geom = geom_loader("lib:methane.xyz")
    name = this_dir / "bases" / fn
    shells = shells_with_basis(geom, name=name)
    assert len(shells) == ref_len
