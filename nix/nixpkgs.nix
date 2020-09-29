let
  # nixos-20.03 at 21.07.2020
  nixpkgs = import (builtins.fetchGit {
    url = "https://github.com/nixos/nixpkgs-channels/";
    name = "nixos-20.09";
    rev = "a9226f2b3a52fcbbc5587d2fa030729e714f40fe";
    ref = "refs/heads/nixos-20.09";
  }) { overlays = [NixWithChemistry]; };

  NixWithChemistry = import ./nixwithchemistry/default.nix;

in nixpkgs
