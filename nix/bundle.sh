#! /usr/bin/env nix-shell
#! nix-shell -i bash -p nix-bundle -I nixpkgs=https://github.com/NixOS/nixpkgs/archive/71326cd12ddfa0fac40fdb451fcba7dad763c56e.tar.gz

export NIX_PATH=nixpkgs=https://github.com/NixOS/nixpkgs/archive/71326cd12ddfa0fac40fdb451fcba7dad763c56e.tar.gz

for i in pysis pysisfilter pysispack pysisplot pysisthermo pysistrj; do
  nix-bundle '(import ./default.nix { })' /bin/$i > $i
done
