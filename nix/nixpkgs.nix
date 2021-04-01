let
  # nixos-20.09 at 19.10.2020
  nixpkgs = import (builtins.fetchGit {
    url = "https://github.com/nixos/nixpkgs";
    name = "nixos-20.09";
    rev = "9cea2bf89b5cbe90d933b2a0b3018692342657b4";
    ref = "refs/heads/nixos-20.09";
  }) { overlays = [NixWithChemistry]; };

  NixWithChemistry = import ./nixwithchemistry/default.nix;

in nixpkgs
