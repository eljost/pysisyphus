{
  description = "Python suite for optimization of stationary points on ground- and excited states PES and determination of reaction paths";

  inputs = {
    nixpkgs.url = "nixpkgs/nixpkgs-unstable";

    qchem = {
      url = "github:markuskowa/nixos-qchem/master";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    flake-compat = {
      url = "github:edolstra/flake-compat";
      flake = false;
    };

    flake-utils.url = "github:numtide/flake-utils";
  };

  nixConfig = {
    # Custom prompt in nix develop shell
    bash-prompt = ''\[\e[0;1;38;5;215m\]pysisyphus\[\e[0;1m\]:\[\e[0;1;38;5;75m\]\w\[\e[0;1m\]$ \[\e[0m\]'';
    extra-subtituters = [ "https://pysisyphus.cachix.org" ];
  };

  outputs = { self, nixpkgs, qchem, flake-utils, ... }:
    flake-utils.lib.eachSystem [ "x86_64-linux" ]
      (system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ qchem.overlays.default (import ./nix/overlay.nix) ];
            config = {
              allowUnfree = true;
              qchem-config = {
                optAVX = true;
                allowEnv = false;
              };
            };
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

            pysisyphusDocker = pkgs.dockerTools.buildImage {
              name = self.packages."${system}".pysisyphus.pname;
              tag = self.packages."${system}".pysisyphus.version;

              fromImageName = null;
              fromImageTag = null;

              copyToRoot = pkgs.buildEnv {
                name = "image-root";
                paths = [
                  # Necessary for interactive usage
                  pkgs.bashInteractive
                  pkgs.coreutils

                  # Computational chemistry software
                  self.packages."${system}".pysisyphus
                ];
                pathsToLink = [ "/bin" ];
              };
            };

            pysisyphusSingularity = pkgs.singularity-tools.buildImage {
              name = self.packages."${system}".pysisyphus.name;
              contents = with pkgs; [ bashInteractive coreutils self.packages."${system}".pysisyphus ];

              # Size of the virtual disk used for building the container.
              # This is NOT the final disk size.
              diskSize = 10000;

              # Memory for the virtualisation environment that BUILDS the image.
              # This is not a runtime parameter of the image.
              memSize = 6000;
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

        }) // {
      overlays.default = import ./nix/overlay.nix;
    };
}
