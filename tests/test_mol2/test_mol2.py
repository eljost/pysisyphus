import textwrap

import pytest

from pysisyphus.io import geom_from_mol2


@pytest.fixture
def benzene_mol2():
    return textwrap.dedent(
        """
# Name: benzene
# Creating user name: ChemicBook
# Creation time: Sat Mar 20 00:18:30 2021

@<TRIPOS>MOLECULE
benzene
 12 12 0 0 0
SMALL
NO_CHARGES
****
Any comment about the molecule goes here. No need for a pound sign here

@<TRIPOS>ATOM
      1 C          -0.7600    1.1691   -0.0005 C.ar    1  BENZENE       0.000
      2 C           0.6329    1.2447   -0.0012 C.ar    1  BENZENE       0.000
      3 C           1.3947    0.0765    0.0004 C.ar    1  BENZENE       0.000
      4 C           0.7641   -1.1677    0.0027 C.ar    1  BENZENE       0.000
      5 C          -0.6288   -1.2432    0.0001 C.ar    1  BENZENE       0.000
      6 C          -1.3907   -0.0751   -0.0015 C.ar    1  BENZENE       0.000
      7 H          -1.3536    2.0792    0.0005 H       1  BENZENE       0.000
      8 H           1.1243    2.2140   -0.0028 H       1  BENZENE       0.000
      9 H           2.4799    0.1355   -0.0000 H       1  BENZENE       0.000
     10 H           1.3576   -2.0778    0.0063 H       1  BENZENE       0.000
     11 H          -1.1202   -2.2126   -0.0005 H       1  BENZENE       0.000
     12 H          -2.4759   -0.1340   -0.0035 H       1  BENZENE       0.000
@<TRIPOS>BOND
     1     1     2   ar
     2     2     3   ar
     3     3     4   ar
     4     4     5   ar
     5     5     6   ar
     6     1     6   ar
     7     1     7    1
     8     2     8    1
     9     3     9    1
    10     4    10    1
    11     5    11    1
    12     6    12    1
"""
    )


def test_benzene_mol2(benzene_mol2):
    # fn = "test.mol2"
    # with open(fn) as handle:
        # text = handle.read()
    geom = geom_from_mol2(benzene_mol2)
    assert len(geom.atoms) == 12
    # geom.jmol()
