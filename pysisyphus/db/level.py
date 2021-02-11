from pysisyphus.calculators import XTB, ORCA


LEVELS = (
    (
        "gfn0_xtb",
        XTB,
        {
            "gfn": 0,
            "quiet": True,
        },
    ),
    (
        "gfn1_xtb",
        XTB,
        {
            "gfn": 1,
            "quiet": True,
        },
    ),
    (
        "gfn2_xtb",
        XTB,
        {
            "gfn": 2,
            "quiet": True,
        },
    ),
    (
        "b973c_orca",
        ORCA,
        {
            "keywords": "b97-3c Grid4 FinalGrid5 tightscf",
            "pal": 4,
        },
    ),
)
