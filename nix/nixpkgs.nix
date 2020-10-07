let
  # nixos-20.09 at 07.10.2020
  nixpkgs = import (builtins.fetchGit {
    url = "https://github.com/nixos/nixpkgs";
    name = "nixos-20.09";
    rev = "7badbf18c45b7490d893452beb8950d966327831";
    ref = "refs/heads/nixos-20.09";
  }) { overlays = [NixWithChemistry]; };

  NixWithChemistry = import ./nixwithchemistry/default.nix;

in nixpkgs
