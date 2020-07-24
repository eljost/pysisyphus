let nixpkgs = import ./nixpkgs.nix;
in  with nixpkgs; python37Packages.callPackage ./pysisyphus.nix {
      orca = null;
      turbomole = null;
      gaussian = null;
      jmol = null;
      multiwfn = null;
      xtb = null;
      openmolcas = null;
      pyscf = null;
      psi4 = null;
      mopac = null;
    }
