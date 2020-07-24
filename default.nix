let nixpkgs = import ./nixpkgs.nix;
in  with nixpkgs; python37Packages.callPackage ./pysisyphus.nix {}
