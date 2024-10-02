{
  description = "Python suite for optimization of stationary points on ground- and excited states PES and determination of reaction paths";

  inputs = {
    qchem.url = "github:nix-qchem/nixos-qchem/master";

    flake-compat = {
      url = "github:edolstra/flake-compat";
      flake = false;
    };

    flake-utils.url = "github:numtide/flake-utils";

    nixBundlers.url = "github:NixOS/bundlers/master";

    pysisyphus-addons.url = "github:eljost/pysisyphus-addons";

    sympleints.url = "github:eljost/sympleints";
  };

  nixConfig = {
    # Custom prompt in nix develop shell
    bash-prompt = ''\[\e[0;1;38;5;215m\]pysisyphus\[\e[0;1m\]:\[\e[0;1;38;5;75m\]\w\[\e[0;1m\]$ \[\e[0m\]'';
    extra-subtituters = [ "https://pysisyphus.cachix.org" ];
  };

  outputs = { self, qchem, sympleints, pysisyphus-addons, flake-utils, nixBundlers, ... }:
    flake-utils.lib.eachSystem [ "x86_64-linux" ]
      (system:
        let
          pkgs = import qchem.inputs.nixpkgs {
            inherit system;
            overlays = [
              sympleints.overlays.default
              pysisyphus-addons.overlays.default
              qchem.overlays.default
              (import ./nix/overlay.nix)
            ];
            config = {
              allowUnfree = true;
              qchem-config = {
                optAVX = true;
                allowEnv = false;
              };
            };
          };

          toSingularityImage = drv: pkgs.singularity-tools.buildImage {
            name = drv.name;
            contents = with pkgs; [
              drv
              bashInteractive
              coreutils
              findutils
              gnused
              which
            ];
            diskSize = 40000;
            memSize = 2000;
          };

        in
        {
          packages = {
            default = self.packages."${system}".pysisyphus;

            pysisyphus = pkgs.pysisyphus.override { };

            pysisyphusOrca = self.packages."${system}".pysisyphus.override { enableOrca = true; };

            pysisyphusTurbomole = self.packages."${system}".pysisyphus.override { enableTurbomole = true; };

            pysisyphusFull = self.packages."${system}".pysisyphus.override {
              enableOrca = true;
              enableTurbomole = true;
              enableCfour = true;
              enableMolpro = true;
              enableGamess = true;
            };
          };

          devShells.default =
            let
              pythonEnv = pkgs.python3.withPackages (p: with p;
                [ pip ]
                ++ p.pysisyphus.nativeBuildInputs
                ++ p.pysisyphus.buildInputs
                ++ p.pysisyphus.propagatedBuildInputs
              );
            in
            pkgs.mkShell {
              buildInputs = [
                pkgs.nixpkgs-fmt
                pkgs.black
                pythonEnv
              ];

              shellHook = ''
                export PYSISRC=${pkgs.python3.pkgs.pysisyphus.passthru.pysisrc}
              '';
            };

          bundlers = {
            inherit toSingularityImage;
            inherit (nixBundlers.bundlers."${system}") toArx toDEB toRPM toDockerImage;
            default = self.bundlers."${system}".toArx;
          };

          formatter = pkgs.nixpkgs-fmt;

          checks = {
            inherit (self.packages."${system}") pysisyphus;
          };
        }) // {
      overlays.default = import ./nix/overlay.nix;

      hydraJobs = self.checks // {
        x86_64-linux = {
          pysisyphusTurbomole = self.packages.x86_64-linux.pysisyphusTurbomole;
          pysisyphusOrca = self.packages.x86_64-linux.pysisyphusOrca;
        };
      };
    };
}
