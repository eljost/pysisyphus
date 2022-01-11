from collections import namedtuple

import psutil
import pytest

from pysisyphus.helpers_pure import check_mem
from pysisyphus.init_logging import init_logging
from pysisyphus.calculators.Calculator import Calculator


init_logging()


def get_mock_mv(available_mb):
    #                          kbytes bytes
    available = available_mb * 1024 * 1024
    vm = namedtuple("vm", ("available",))

    def mock_virtual_memory():
        return vm(available=available)

    return mock_virtual_memory


@pytest.mark.parametrize(
    "pal",
    range(1, 4),
)
@pytest.mark.parametrize(
    "mem",
    (1000, 8000),
)
def test_check_mem(pal, mem, monkeypatch):
    available_mb = 10_000
    mock_virtual_memory = get_mock_mv(available_mb)
    monkeypatch.setattr(psutil, "virtual_memory", mock_virtual_memory)

    avail_frac = 0.5
    mb_corr = check_mem(mem, pal, avail_frac=avail_frac)
    assert mb_corr == int(min(mem, available_mb * avail_frac / pal))


def test_check_mem_calc(monkeypatch):
    available_mb = 10_000
    mock_virtual_memory = get_mock_mv(available_mb)
    monkeypatch.setattr(psutil, "virtual_memory", mock_virtual_memory)
    calc = Calculator(pal=10, mem=20_000)
    assert calc.mem == 850
