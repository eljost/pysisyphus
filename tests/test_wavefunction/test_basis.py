from pysisyphus.wavefunction.Basis import basis_from_orca_str, basis_from_pyscf_str


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
