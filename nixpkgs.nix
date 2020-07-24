let
  # nixos-20.03 at 21.07.2020
  nixpkgs = import (builtins.fetchGit {
    url = "https://github.com/nixos/nixpkgs-channels/";
    name = "nixos-20.03-9ea61f7";
    rev = "9ea61f7bc4454734ffbff73c9b6173420fe3147b";
    ref = "refs/heads/nixos-20.03";
  }) { overlays = [NixWithChemistry]; };

  NixWithChemistry =
    let
      repoPath = builtins.fetchGit {
        url = "https://gitlab.com/theoretical-chemistry-jena/nixwithchemistry.git";
        name = "NixWithChemistry";
        rev = "6e47f912176b70122cc337f14838e59418af867e";
        ref = "refs/heads/master";
      };
    in import "${repoPath}/default.nix";

in nixpkgs
