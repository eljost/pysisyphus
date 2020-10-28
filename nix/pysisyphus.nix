{ fetchFromGitHub, buildPythonApplication, lib, writeTextFile, writeScript, makeWrapper
# Python dependencies
, autograd
, dask
, distributed
, h5py
, jinja2
, matplotlib
, numpy
, natsort
, pytest
, pyyaml
, rmsd
, scipy
, sympy
, bash
, orca ? orca # or null
, turbomole ? turbomole # or null
, gaussian ? gaussian # or null
, jmol ? jmol # or null
, multiwfn ? multiwfn # or null
, xtb ? xtb # or null
, openmolcas ? openmolcas # or null
, pyscf ? pyscf # or null
, psi4 ? psi4 # or null
, mopac ? null
, wfoverlap ? wfoverlap # or null
, nwchem ? nwchem # or null
, gamess-us ? gamess-us # or null
, cfour ? cfour  # or null
, molpro ? molpro # or null
, qcengine
, ase
}:
let
  psi4Wrapper = writeScript "psi4.sh" ''
    #!${bash}/bin/bash
    ${psi4}/bin/psi4 -o stdout $1
  '';
  pysisrc =
    let
      gaussian16Conf = {
        cmd = "${gaussian}/bin/g16";
        formchk_cmd = "${gaussian}/bin/formchk";
        unfchk_cmd = "${gaussian}/bin/unfchk";
      };
      openmolcasConf = {
        cmd = "${openmolcas}/bin/pymolcas";
      };
      orcaConf = {
        cmd = "${orca}/bin/orca";
      };
      mopacConf = {
      };
      psi4Conf = {
        cmd = "${toString psi4Wrapper}";
      };
      xtbConf = {
        cmd = "${xtb}/bin/xtb";
      };
      wfoverlapConf = {
        cmd = "${wfoverlap}/bin/wfoverlap.x";
      };
      multiwfnConf = {
        cmd = "${multiwfn}/bin/Multiwfn";
      };
      jmolConf = {
        cmd = "${jmol}/bin/jmol";
      };
      text = with lib.lists; lib.generators.toINI {} (builtins.listToAttrs ([]
        ++ optional (gaussian != null) { name = "gaussian16"; value = gaussian16Conf; }
        ++ optional (openmolcas != null) { name = "openmolcas"; value = openmolcasConf; }
        ++ optional (orca != null) { name = "orca"; value = orcaConf; }
        ++ optional (mopac != null) { name = "mopac"; value = mopacConf; }
        ++ optional (psi4 != null) { name = "psi4"; value = psi4Conf; }
        ++ optional (xtb != null) { name = "xtb"; value = xtbConf; }
        ++ optional (wfoverlap != null) { name = "wfoverlap"; value = wfoverlapConf; }
        ++ optional (multiwfn != null) { name = "multiwfn"; value = multiwfnConf; }
        ++ optional (jmol != null) { name = "jmol"; value = jmolConf; }
      ));
    in
      writeTextFile {
        inherit text;
        name = "pysisrc";
      };

  pname = "pysisyphus";
  version = "dev";
  description = "Python suite for optimization of stationary points on ground- and excited states PES and determination of reaction paths";
  homepage = "https://github.com/eljost/pysisyphus";

  pkgConfig = writeTextFile {
    name = "${pname}.pc";
    text =  ''
      prefix=@out@
      exec_prefix=''${prefix}
      libdir=''${exec_prefix}/lib
      sharedlibdir=''${libdir}
      includedir=''${prefix}/include

      Name: ${pname}
      Description: Python suite for geometry optimisations
      Version: ${version}

      Requires:
      Libs:
      Cflags:
      URL: ${homepage}
    '';
  };

in
  buildPythonApplication rec {
    inherit pname version;

    nativeBuildInputs = [
      makeWrapper
    ];

    propagatedBuildInputs = with lib.lists; [
        autograd
        dask
        distributed
        h5py
        jinja2
        matplotlib
        numpy
        natsort
        pytest
        pyyaml
        rmsd
        scipy
        sympy
        qcengine
        ase
    ] ++ optional (orca != null) orca
      ++ optional (turbomole != null) turbomole
      ++ optional (gaussian != null) gaussian
      ++ optional (jmol != null) jmol
      ++ optional (multiwfn != null) multiwfn
      ++ optional (xtb != null) xtb
      ++ optional (openmolcas != null) openmolcas
      ++ optional (pyscf != null) pyscf
      ++ optional (psi4 != null) psi4
      ++ optional (mopac != null) mopac
      ++ optional (wfoverlap != null) wfoverlap
      ++ optional (nwchem != null) nwchem
      ++ optional (gamess-us != null) gamess-us
      ++ optional (cfour != null) cfour
      ++ optional (molpro != null) molpro
    ;

    src = builtins.path {
      name = "pysisyphus";
      path = ./..;
    };

    doCheck = false;

    binSearchPath = lib.makeSearchPath "bin" ([ ]
      ++ lib.optional (orca != null) orca
      ++ lib.optional (turbomole != null) turbomole
      ++ lib.optional (gaussian != null) gaussian
      ++ lib.optional (jmol != null) jmol
      ++ lib.optional (multiwfn != null) multiwfn
      ++ lib.optional (xtb != null) xtb
      ++ lib.optional (openmolcas != null) openmolcas
      ++ lib.optional (pyscf != null) pyscf
      ++ lib.optional (psi4 != null) psi4
      ++ lib.optional (mopac != null) mopac
      ++ lib.optional (wfoverlap != null) wfoverlap
      ++ lib.optional (nwchem != null) nwchem
      ++ lib.optional (gamess-us != null) gamess-us
      ++ lib.optional (cfour != null) cfour
      ++ lib.optional (molpro != null) molpro
    );

    postInstall = ''
      mkdir -p $out/share/pysisyphus
      cp ${pysisrc} $out/share/pysisyphus/pysisrc

      for exe in $out/bin/*; do
        wrapProgram $exe \
          --prefix PATH : ${binSearchPath} \
          --set-default "PYSISRC" "$out/share/pysisyphus/pysisrc"
      done

      mkdir -p $out/lib/pkgconfig
      substitute ${pkgConfig} $out/lib/pkgconfig/${pkgConfig.name} \
        --subst-var "out"
    '';

    meta = with lib; {
      inherit description homepage;
      license = licenses.gpl3;
      platforms = platforms.unix;
    };
  }
