# Flake compatibility definition for the legacy package set as if it had never
# seen flakes
(import ./default.nix).legacyPackages."${builtins.currentSystem}"
