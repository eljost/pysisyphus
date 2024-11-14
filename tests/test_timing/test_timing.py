import pytest
from pysisyphus import timing


@pytest.mark.parametrize(
    "dur, ref",
    (
        (60, "1m 0s"),
        (60**2, "1h 0m 0s"),
        (60**2 * 24, "1d 0h 0m 0s"),
        (17 * 60**2 * 24, "17d 0h 0m 0s"),
        (59.12345, "59s"),
    ),
)
def test_timing(dur, ref):
    rendered = timing.render(dur)
    # print(rendered)
    assert rendered == ref
