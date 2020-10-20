let nixpkgs = import ./nixpkgs.nix;
    config = import ./nixwithchemistry/config.nix;
in  with nixpkgs; python3Packages.callPackage ./pysisyphus.nix {
      orca = if config.orca then orca else null;
      turbomole = if config.turbomole then turbomole else null;
      gaussian = if config.gaussian then gaussian else null;
      gamess-us = if config.gamess-us then gamess-us else null;
      cfour = if config.cfour then cfour else null;
      mopac = null;
    }
