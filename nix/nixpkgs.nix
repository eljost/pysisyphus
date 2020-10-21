let
  # nixos-20.09 at 19.10.2020
  nixpkgs = import (builtins.fetchGit {
    url = "https://github.com/nixos/nixpkgs";
    name = "nixos-20.09";
    rev = "51aaa3fa1b69559456f9bd4968bd5b179a784f67";
    ref = "refs/heads/nixos-20.09";
  }) { overlays = [NixWithChemistry]; };

  /*
  Containts the fixes for fsspec python package already.
  */
  /*
  nixpkgs = import (builtins.fetchGit {
    url = "https://github.com/nixos/nixpkgs";
    name = "release-20.09";
    rev = "4188bf40dd247f7af3d6105f1c7097a4478ad376";
    ref = "refs/heads/release-20.09";
  }) { overlays = [NixWithChemistry]; };
  */

  NixWithChemistry = import ./nixwithchemistry/default.nix;

in nixpkgs
