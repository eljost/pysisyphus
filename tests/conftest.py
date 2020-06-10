from collections import OrderedDict
from pathlib import Path

import pytest

"""
    Adapted from the pytest-harvest plugin.
"""


FixtureStore = OrderedDict()


class ResultsBag(dict):
    
    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]

    def __hash__(self):
        return id(self)


@pytest.fixture(scope="session")
def results_bag(request):
    bag = ResultsBag()
    test_name = request._pyfuncitem.name
    FixtureStore.setdefault("results_bag", dict())[test_name] = bag
    return bag


@pytest.fixture(scope="session")
def fixture_store(request):
    return FixtureStore


@pytest.fixture(scope="module")
def this_dir(request):
    path = Path(request.fspath)
    return path.parent
