final: prev: {
  pysisyphus = prev.python3.pkgs.toPythonApplication final.python3.pkgs.pysisyphus;

  qchem = prev.qchem // {
    orca = prev.qchem.orca.overrideAttrs (oldAttrs: {
      version = "5.0.4";

      src = prev.requireFile {
        name = "orca_5_0_4_linux_x86-64_shared_openmpi411.tar.xz";
        sha256 = "sha256-xOpa6mDae8sYprcEJgkgb76yp2XG+pWMVonUULWIsDY=";
        message = "Please add orca_5_0_4_linux_x86-64_shared_openmpi411.tar.xz to the nix store";
      };

      installPhase = ''
        mkdir -p $out/bin $out/lib $out/share/doc/orca

        cp autoci_* $out/bin
        cp orca_* $out/bin
        cp orca $out/bin
        cp otool_* $out/bin

        cp -r ORCACompoundMethods $out/bin/.

        cp *.so.5 $out/lib/.

        cp *.pdf $out/share/doc/orca

        wrapProgram $out/bin/orca --prefix PATH : '${final.openmpi}/bin:${final.openssh}/bin'

        ln -s ${final.qchem.xtb}/bin/xtb $out/bin/otool_xtb
      '';
    });
  };

  python3 = prev.python3.override (old: {
    packageOverrides = prev.lib.composeExtensions (old.packageOverrides or (_: _: { })) (pfinal: pprev: {
      pysisyphus = pfinal.callPackage ./pysisyphus.nix {
        inherit (final.qchem)
          multiwfn
          xtb
          molcas
          psi4
          wfoverlap
          nwchem
          orca # Requires ORCA 5.0 series
          turbomole
          gaussian
          cfour
          molpro
          gamess-us;
      };
    });
  });
}
