import pytest

from pysisyphus.helpers_pure import increment_fn



@pytest.mark.parametrize(
        "org_fn, suffix, ref_fn", (
            ("opt", "rebuilt", "opt_rebuilt_000"),
            ("opt_rebuilt_000", "rebuilt", "opt_rebuilt_001"),
            ("opt", None, "opt_000"),
            ("opt_000", None, "opt_001"),
            )
)
def test_increment_fn(org_fn, suffix, ref_fn):
    incr_fn = increment_fn(org_fn, suffix=suffix)
    assert incr_fn == ref_fn
