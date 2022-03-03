import pytest

from pysisyphus.helpers_pure import full_expand



@pytest.mark.parametrize(
    "to_expand, expected", (
    ("0..10", list(range(10))),
    ("0..3,6,8,9..12", [0, 1, 2, 6, 8, 9, 10, 11]),
    (list(range(10)), list(range(10))),
    (range(10), list(range(10))),
    (1, [1, ]),
))
def test_full_expand(to_expand, expected):
    expanded = full_expand(to_expand)
    assert expanded == expected
