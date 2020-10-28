let
  # nixos-20.09 at 19.10.2020
  nixpkgs = import (builtins.fetchGit {
    url = "https://github.com/nixos/nixpkgs";
    name = "nixos-20.09";
    rev = "13d0c311e3ae923a00f734b43fd1d35b47d8943a";
    ref = "refs/heads/nixos-20.09";
  }) { overlays = [NixWithChemistry]; };

  NixWithChemistry = import ./nixwithchemistry/default.nix;

in nixpkgs
