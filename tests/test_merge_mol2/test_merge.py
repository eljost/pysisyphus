import filecmp
import tempfile

import pytest

from pysisyphus.drivers.merge_mol2 import merge_mol2_geoms


@pytest.mark.parametrize(
    "fn1, fn2, bonds, del1, del2, ref",
    (
        (
            "phalloidin-F_linker.mol2",
            "rhodamin_amid.mol2",
            ((106, 28),),
            (103, 104, 105, 107, 125),
            (64,),
            "rhodamin_merged.mol2",
        ),
        (
            "phalloidin-F_linker.mol2",
            "alexa-488_amid.mol2",
            ((106, 11),),
            (103, 104, 105, 107, 125),
            (49,),
            "alexa_merged.mol2",
        ),
        (
            "dy-490_amid.mol2",
            "linker.mol2",
            ((9, 16),),
            (63,),
            (1, 19, 22, 23, 24),
            "dy_linker.mol2",
        ),
        (
            "dy_linker.mol2",
            "phallo.mol2",
            ((63, 36),),
            (82,),
            (43, 44),
            "dy_linker_phallo.mol2",
        ),
    ),
)
def test_merge_mol2_geoms(fn1, fn2, bonds, del1, del2, ref, this_dir):
    merged = merge_mol2_geoms(this_dir / fn1, this_dir / fn2, bonds, del1, del2)

    with tempfile.NamedTemporaryFile(mode="w") as handle:
        handle.write(merged)
        handle.flush()
        assert filecmp.cmp(handle.name, this_dir / ref)
