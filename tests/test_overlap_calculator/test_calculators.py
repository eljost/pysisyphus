import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.calculators import DFTBp, Gaussian16, ORCA, Turbomole
from pysisyphus.testing import using


_ROOT_REF_ENERGIES = {
    # -74.9607 au is the GS
    # -74.4893 au is the S1, not the GS
    (ORCA, None): -74.9607,
    (ORCA, 2): -74.3984,
    (PySCF, None): -74.9607,
    (PySCF, 2): -74.3984,  # ???
    (Turbomole, None): -74.9607,
    (Turbomole, 2): -74.39837,
    (Gaussian16, None): -74.9607,
    (Gaussian16, 2): -74.39837,
    (DFTBp, None): -4.077751,
    (DFTBp, 2): -3.33313,
}


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs",
    (
        pytest.param(
            ORCA,
            {
                "keywords": f"hf sto-3g",
                "blocks": "%tddft tda false nroots 3 end",
            },
            marks=(using("orca")),
        ),
        pytest.param(
            PySCF,
            {
                "basis": "sto3g",
                "method": "tdhf",
                "nroots": 3,
            },
            marks=(using("pyscf")),
        ),
        pytest.param(
            DFTBp,
            {
                "parameter": "mio-ext",
                "nroots": 4,
            },
            marks=(using("dftbp")),
        ),
        pytest.param(
            Turbomole,
            {
                "simple_input": {
                    "basis": "sto-3g hondo",
                    "scfinstab": "rpas",
                    "soes": {
                        "a": 3,
                    },
                }
            },
            marks=(using("turbomole")),
        ),
        pytest.param(
            Gaussian16,
            {"route": "hf sto-3g td=(nstates=3)"},
            marks=(using("gaussian16")),
        ),
    ),
)
@pytest.mark.parametrize("root", (None, 2))
def test_root(calc_cls, calc_kwargs, root):
    ref_energy = _ROOT_REF_ENERGIES[(calc_cls, root)]
    calc_kwargs.update(
        {
            "root": root,
            "base_name": f"root_{root}",
            "track": True,
        }
    )
    geom = calc_cls.geom_from_fn("lib:h2o.xyz", **calc_kwargs)
    energy = geom.energy
    assert energy == pytest.approx(ref_energy)
