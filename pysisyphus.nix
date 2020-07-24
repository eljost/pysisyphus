{ fetchFromGitHub, buildPythonApplication, lib
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
}:
buildPythonApplication rec {
    pname = "pysisyphus";
    version = "dev";

    propagatedBuildInputs = [
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
      ];

    src = builtins.path {
      name = "pysisyphus";
      path = ./.;
    };

    doCheck = false;

    meta = with lib; {
      description = "Python suite for optimization of stationary points on ground- and excited states PES and determination of reaction paths";
      license = licenses.gpl3;
      homepage = "https://github.com/eljost/pysisyphus";
      platforms = platforms.unix;
    };
}
