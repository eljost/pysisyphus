let
  # nixos-20.03 at 21.07.2020
  nixpkgs = import (builtins.fetchGit {
    url = "https://github.com/nixos/nixpkgs-channels/";
    name = "nixos-20.03-9ea61f7";
    rev = "9ea61f7bc4454734ffbff73c9b6173420fe3147b";
    ref = "refs/heads/nixos-20.03";
  }) { overlays = [NixWithChemistry]; };

  NixWithChemistry = import ./nixwithchemistry/default.nix;

in nixpkgs
