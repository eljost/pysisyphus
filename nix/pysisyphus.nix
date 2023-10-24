{ buildPythonPackage
, lib
, writeTextFile
, writeScript
, makeWrapper
, mpiCheckPhaseHook
, pytestCheckHook
, nix-gitignore
  # Python dependencies
, setuptools-scm
, autograd
, dask
, distributed
, h5py
, jinja2
, matplotlib
, numpy
, natsort
, pyyaml
, rmsd
, scipy
, sympy
, scikit-learn
, Fabric
, psutils
, qcengine
, ase
, xtb-python
, openbabel-bindings
, pyscf
  # Runtime dependencies
, runtimeShell
, jmol
, enableJmol ? true
, multiwfn
, enableMultiwfn ? true
, xtb
, enableXtb ? true
, molcas
, enableMolcas ? true
, psi4
, enablePsi4 ? true
, wfoverlap
, enableWfoverlap ? true
, nwchem
, enableNwchem ? true
, orca
, enableOrca ? false
, turbomole
, enableTurbomole ? false
, gaussian
, enableGaussian ? false
, cfour
, enableCfour ? false
, molpro
, enableMolpro ? false
, gamess-us
, enableGamess ? false
  # Test dependencies
, openssh
}:
let
  psi4Wrapper = writeScript "psi4.sh" ''
    #!${runtimeShell}
    ${psi4}/bin/psi4 -o stdout $1
  '';
  pysisrc =
    let
      gaussian16Conf = {
        cmd = "${gaussian}/bin/g16";
        formchk = "${gaussian}/bin/formchk";
        unfchk = "${gaussian}/bin/unfchk";
        rwfdump = "${gaussian}/bin/rwfdump";
      };
      text = lib.generators.toINI { } (builtins.listToAttrs ([ ]
        ++ lib.optional enableMolcas { name = "openmolcas"; value.cmd = "${molcas}/bin/pymolcas"; }
        ++ lib.optional enablePsi4 { name = "psi4"; value.cmd = "${psi4Wrapper}"; }
        ++ lib.optional enableWfoverlap { name = "wfoverlap"; value.cmd = "${wfoverlap}/bin/wfoverlap.x"; }
        ++ lib.optional enableMultiwfn { name = "mwfn"; value.cmd = "${multiwfn}/bin/Multiwfn"; }
        ++ lib.optional enableJmol { name = "jmol"; value.cmd = "${jmol}/bin/jmol"; }
        ++ lib.optional enableXtb { name = "xtb"; value.cmd = "${xtb}/bin/xtb"; }
        ++ lib.optional enableGaussian { name = "gaussian16"; value = gaussian16Conf; }
        ++ lib.optional enableOrca { name = "orca"; value.cmd = "${orca}/bin/orca"; }
        ++ lib.optional enableGamess { name = "gamess"; value.cmd = "${gamess-us}/bin/rungms"; }
      ));
    in
    writeTextFile {
      inherit text;
      name = "pysisrc";
    };

  binSearchPath = lib.makeSearchPath "bin" ([ ]
    ++ lib.optional enableJmol jmol
    ++ lib.optional enableMultiwfn multiwfn
    ++ lib.optional enableXtb xtb
    ++ lib.optional enableMolcas molcas
    ++ lib.optional enablePsi4 psi4
    ++ lib.optional enableWfoverlap wfoverlap
    ++ lib.optional enableNwchem nwchem
    ++ lib.optional enableOrca orca
    ++ lib.optional enableTurbomole turbomole
    ++ lib.optional enableGaussian gaussian
    ++ lib.optional enableCfour cfour
    ++ lib.optional enableMolpro molpro
    ++ lib.optional enableGamess gamess-us
  );

in
buildPythonPackage rec {
  pname = "pysisyphus";
  version = "0.8.0b0";

  nativeBuildInputs = [
    makeWrapper
    setuptools-scm
    mpiCheckPhaseHook
  ];

  propagatedBuildInputs = [
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
    scikit-learn
    Fabric
    psutils
    qcengine
    ase
    openbabel-bindings
    openssh
    pyscf
  ] # Syscalls
  ++ lib.optional enableXtb xtb-python
  ++ lib.optional enableXtb xtb
  ++ lib.optional enableJmol jmol
  ++ lib.optional enableMultiwfn multiwfn
  ++ lib.optional enableMolcas molcas
  ++ lib.optional enablePsi4 psi4
  ++ lib.optional enableWfoverlap wfoverlap
  ++ lib.optional enableNwchem nwchem
  ++ lib.optional enableOrca orca
  ++ lib.optional enableTurbomole turbomole
  ++ lib.optional enableGaussian gaussian
  ++ lib.optional enableCfour cfour
  ++ lib.optional enableMolpro molpro
  ++ lib.optional enableGamess gamess-us
  ;

  src = nix-gitignore.gitignoreSource [ ] ../.;

  pyproject = true;

  preBuild = "export SETUPTOOLS_SCM_PRETEND_VERSION=${version}";

  checkInputs = [ openssh pytestCheckHook ];

  preCheck = ''
    export PYSISRC=${pysisrc}
    export PATH=$PATH:${binSearchPath}
    export XTBPATH=${xtb}/share/xtb
  '';

  pytestFlagsArray =
    let
      faultyPyscf = [
        "ohch3f_anion"
        "afir_hessian"
        "opt_restart"
        "dimer_hcn"
        "numhess"
        "backtransform_hessian"
        "hcn_iso_gs2"
        "water_fd_hessian"
        "gradient"
        "hcn_neb_dimer_irc"
        "test_ci_opt[RFOptimizer-opt_kwargs1--78.2487951]" # Broken as of PySCF >= 2.3 as a DFT functional definition was changed
        "test_composite_oniom[lib:subst_effect/toluene_minus_H_b3lypG_631g.xyz-2-high_inds1--270.824806805671]" # Slight numerical deviations
      ];
    in
    [
      "-v"
      "--show-capture=no"
      "--durations=0"
      "-m 'not benchmark and not skip_ci'"
      "-k '${lib.strings.concatMapStringsSep " and " (s: "not ${s}") faultyPyscf}'"
      "tests"
    ];

  pythonImportsCheck = [ "pysisyphus" ];

  postInstall = ''
    mkdir -p $out/share/pysisyphus
    cp ${pysisrc} $out/share/pysisyphus/pysisrc
    for exe in $out/bin/*; do
      wrapProgram $exe \
        ${if binSearchPath == "" then "" else "--prefix PATH : ${binSearchPath}"} \
        --set-default PYSISRC $out/share/pysisyphus/pysisrc \
        --set SCRATCH "./"
    done
  '';

  passthru = {
    inherit
      pysisrc
      enableXtb
      enableJmol
      enableMultiwfn
      enableMolcas
      enablePsi4
      enableWfoverlap
      enableNwchem
      enableOrca
      enableTurbomole
      enableGaussian
      enableCfour
      enableMolpro
      enableGamess
      ;
  };

  meta = with lib; {
    description = "Python suite for optimization of stationary points on ground- and excited states PES and determination of reaction paths";
    homepage = "https://github.com/eljost/pysisyphus";
    license = licenses.gpl3Plus;
    platforms = platforms.linux;
    maintainers = [ maintainers.sheepforce ];
    mainProgram = "pysis";
  };
}
