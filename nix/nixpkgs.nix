let
  # nixos-20.09 at 19.10.2020
  nixpkgs = import (builtins.fetchGit {
    url = "https://github.com/nixos/nixpkgs";
    name = "nixos-20.09";
    rev = "d105075a1fd870b1d1617a6008cb38b443e65433";
    ref = "refs/heads/nixos-20.09";
  }) { overlays = [NixWithChemistry]; };

  NixWithChemistry = import ./nixwithchemistry/default.nix;

in nixpkgs
