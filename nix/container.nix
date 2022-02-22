let
  pkgs = import ./pkgs.nix;
  pysisyphus = pkgs.python3.pkgs.pysisyphus.overrideAttrs (_: {
    doCheck = false;
    doInstallCheck = false;
  });

  dockerContainer = pkgs.dockerTools.buildImage {
    name = pysisyphus.pname;
    tag = pysisyphus.version;

    fromImageName = null;
    fromImageTag = null;

    contents = [
      # Necessary for interactive usage
      pkgs.bashInteractive
      pkgs.coreutils

      # Computational chemistry software
      pysisyphus
    ];
  };

  singularityContainer = pkgs.runCommand "${pysisyphus.pname}.sif" {} ''
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
      ${pkgs.singularity}/bin/singularity build /out/${pysisyphus.pname}.sif docker-archive://${dockerContainer}
    cp ${pysisyphus.pname}.sif $out
  '';

in { inherit dockerContainer singularityContainer; }
