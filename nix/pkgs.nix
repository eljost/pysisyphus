let
  sources = import ./sources.nix;
  qchemOverlay = import sources.NixOS-QChem;
  pysisOverlay = import ./overlay.nix;
  postOverlay = self: super:
    { # Example to change MPI/BLAS+LAPACK providers. Uncomment to switch to MKL and MVAPICH, i.e.
      /*
      blas = super.blas.override { blasProvider = super.mkl; };
      lapack = super.lapack.override { lapackProvider = super.mkl; };
      mpi = super.mvapich;
      */
    };
  nixpkgs = import sources.nixpkgs {
    overlays = [ qchemOverlay pysisOverlay postOverlay ];
    allowUnfree = false; # Change to true if you want to use proprietary software
    qchem-config = {
      optAVX = false; # Change to true if your CPU supports at least AVX2 (Haswell upwards)
      useCuda = false;
    };
  };

in nixpkgs
