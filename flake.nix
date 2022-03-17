{
  description =
    "Python suite for optimization of stationary points on ground- and excited states PES and determination of reaction paths";

  inputs = {
    nixpkgs.url = "nixpkgs/nixpkgs-unstable";

    qchem.url = "github:markuskowa/nixos-qchem/master";

    flake-compat = {
      url = "github:edolstra/flake-compat";
      flake = false;
    };

    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, qchem, flake-utils, ... }:
    flake-utils.lib.eachSystem [ "x86_64-linux" ] (system:

      let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ qchem.overlay (import ./nix/overlay.nix) ];
          config = {
            allowUnfree = true;
            qchem-config = {
              optAVX = true;
              allowEnv = true;
            };
          };
        };
      in {

        legacyPackages = pkgs;

        packages = {
          default = self.packages."${system}".pysisyphus;
          pysisyphus = pkgs.python3.pkgs.pysisyphus;

          pysisyphusLib = self.packages."${system}".pysisyphus;

          pysisyphusApp = pkgs.pysisyphus;

          pysisyphusOrca = pkgs.python3.pkgs.pysisyphus.override { enableOrca = true; };
          pysisyphusTurbomole = pkgs.python3.pkgs.pysisyphus.override { enableTurbomole = true; };
          pysisyphusFull = pkgs.python3.pkgs.pysisyphus.override {
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

            contents = [
              # Necessary for interactive usage
              pkgs.bashInteractive
              pkgs.coreutils

              # Computational chemistry software
              self.packages."${system}".pysisyphus
            ];
          };

          pysisyphusSingularity =
            pkgs.runCommand "${self.packages."${system}".pysisyphus.pname}.sif"
            { } ''
              cp -r /etc etc
              touch etc/resolv.conf
              ${pkgs.bubblewrap}/bin/bwrap \
                --ro-bind /nix /nix \
                --ro-bind etc /etc \
                --ro-bind ${pkgs.bash}/bin /usr/local/bin \
                --bind . /out \
                --dev-bind /dev /dev \
                --proc /proc \
                --uid 1000 \
                ${pkgs.singularity}/bin/singularity build /out/${
                  self.packages."${system}".pysisyphus.pname
                }.sif docker-archive://${
                  self.packages."${system}".pysisyphusDocker
                }
              cp ${self.packages."${system}".pysisyphus.pname}.sif $out
            '';
        };

        defaultPackage = self.packages."${system}".pysisyphusApp;

        defaultApp = {
          type = "app";
          program = "${self.packages."${system}".pysisyphusApp}/bin/pysis";
        };

        devShells.default = let
          pythonEnv = pkgs.python3.withPackages (p:
            with p;
            [ pip ] ++ p.pysisyphus.nativeBuildInputs ++ p.pysisyphus.buildInputs
            ++ p.pysisyphus.propagatedBuildInputs);
        in pkgs.mkShell {
          buildInputs = [ pkgs.black pythonEnv ];

          shellHook = ''
            export PYSISRC=${pkgs.python3.pkgs.pysisyphus.passthru.pysisrc}
          '';
        };

        devShell = self.devShells."${system}".default;

      }) // {
        overlays.default = import ./nix/overlay.nix;
        overlay = self.overlays.default;
      };
}
