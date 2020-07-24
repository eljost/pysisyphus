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
, orca ? null
, turbomole ? null
, gaussian ? null
, jmol ? null
, multiwfn ? null
, xtb ? null
, openmolcas ? null
, pyscf ? null
, psi4 ? null
, mopac ? null
, wfoverlap ? null
}:
let
  psi4Wrapper = writeScript {
    name = "psi4.sh";
    text = ''
      ${psi4}/bin/psi4 -o stdout $1
    '';
  };
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
        cmd = "${psi4Wrapper}";
      };
      xtbConf = {
        cmd = "${xtb}/bin/xtb";
      };
      wfoverlapConf = {
        cmd = "${wfoverlap}/bin/wfoverlap";
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
in
  buildPythonApplication rec {
    pname = "pysisyphus";
    version = "dev";

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
    ]
    ++ optional (orca != null) orca
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
    ;

    src = builtins.path {
      name = "pysisyphus";
      path = ./.;
    };

    doCheck = false;

    postInstall = ''
      mkdir -p $out/share/pysisyphus
      cp ${pysisrc} $out/share/pysisyphus/pysisrc

      for exe in $out/bin/*; do
        wrapProgram $exe \
          --set-default "PYSISRC" "$out/share/pysisyphus/pysisrc"
      done;
    '';

    meta = with lib; {
      description = "Python suite for optimization of stationary points on ground- and excited states PES and determination of reaction paths";
      license = licenses.gpl3;
      homepage = "https://github.com/eljost/pysisyphus";
      platforms = platforms.unix;
    };
  }
