{ fetchPypi, fetchFromGitHub, buildPythonPackage, lib, writeTextFile, writeScript, makeWrapper,
 # Python dependencies
 autograd, dask, distributed, h5py, jinja2, matplotlib, numpy, natsort, pytest, pyyaml, rmsd, scipy,
 sympy, scikitlearn, qcengine, ase, xtb-python, openbabel-bindings,
 # Runtime dependencies
 bash, jmol, multiwfn, xtb, openmolcas, pyscf, psi4, wfoverlap, nwchem, orca ? null ,
 turbomole ? null, gaussian ? null, gamess-us ? null, cfour ? null, molpro ? null,
 # Test dependencies
 openssh,
 # Configuration
 fullTest ? false
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
      psi4Conf = {
        cmd = "${psi4Wrapper}";
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
      text = with lib.lists; lib.generators.toINI {} (builtins.listToAttrs ([
        { name = "openmolcas"; value = openmolcasConf; }
        { name = "psi4"; value = psi4Conf; }
        { name = "wfoverlap"; value = wfoverlapConf; }
        { name = "multiwfn"; value = multiwfnConf; }
        { name = "jmol"; value = jmolConf; }
        { name = "xtb"; value = xtbConf; }
      ] ++ optional (gaussian != null) { name = "gaussian16"; value = gaussian16Conf; }
        ++ optional (orca != null) { name = "orca"; value = orcaConf; }
        ++ optional (jmol != null) { name = "jmol"; value = jmolConf; }
      ));
    in
      writeTextFile {
        inherit text;
        name = "pysisrc";
      };

in
  buildPythonPackage rec {
    pname = "pysisyphus";
    version = "55ce1636358c184c08aabf072a2e82d5b0eea21d";

    nativeBuildInputs = [ makeWrapper ];

    propagatedBuildInputs = with lib; [
      autograd
      dask
      distributed
      h5py
      jinja2
      matplotlib
      numpy
      natsort
      pyyaml
      rmsd
      scipy
      sympy
      scikitlearn
      qcengine
      ase
      xtb-python
      openbabel-bindings
      # Syscalls
      jmol
      multiwfn
      xtb
      openmolcas
      pyscf
      psi4
      wfoverlap
      nwchem
    ] ++ lists.optional (orca != null) orca
      ++ lists.optional (turbomole != null) turbomole
      ++ lists.optional (gaussian != null) gaussian
      ++ lists.optional (cfour != null) cfour
      ++ lists.optional (molpro != null) molpro
    ;

    src = fetchFromGitHub {
      owner = "eljost";
      repo = pname;
      rev = "55ce1636358c184c08aabf072a2e82d5b0eea21d";
      sha256 = "1b9knz6dms918irr84nqyn99m6ly0pqs0c5p9ygcihwjsd9rq40s";
    };

    doCheck = true;

    checkInputs = [
      pytest
      openssh
    ];

    checkPhase = ''
      export PYSISRC=${pysisrc}
      export PATH=$PATH:${binSearchPath}
      export OMPI_MCA_rmaps_base_oversubscribe=1

      ${if fullTest then "pytest -v tests" else "pytest -v --pyargs pysisyphus.tests"}
    '';

    binSearchPath = with lib; makeSearchPath "bin" ([
      jmol
      multiwfn
      xtb
      openmolcas
      pyscf
      psi4
      wfoverlap
      nwchem
    ] ++ lists.optional (orca != null) orca
      ++ lists.optional (turbomole != null) turbomole
      ++ lists.optional (gaussian != null) gaussian
      ++ lists.optional (cfour != null) cfour
      ++ lists.optional (molpro != null) molpro
    );

    postInstall = ''
      mkdir -p $out/share/pysisyphus
      cp ${pysisrc} $out/share/pysisyphus/pysisrc
      for exe in $out/bin/*; do
        wrapProgram $exe \
          --prefix PATH : ${binSearchPath} \
          --set-default "PYSISRC" "$out/share/pysisyphus/pysisrc"
      done
    '';

    meta = with lib; {
      description = "Python suite for optimization of stationary points on ground- and excited states PES and determination of reaction paths";
      homepage = "https://github.com/eljost/pysisyphus";
      license = licenses.gpl3;
      platforms = platforms.linux;
    };
  }
