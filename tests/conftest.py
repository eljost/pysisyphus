from collections import OrderedDict
from pathlib import Path

import pytest

"""
    Adapted from the pytest-harvest plugin.
"""


FixtureStore = OrderedDict((("results_bag", dict()),))


class ResultsBag(dict):
    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]

    def __hash__(self):
        return id(self)


@pytest.fixture(scope="function")
def results_bag(request):
    bag = ResultsBag()
    test_name = request._pyfuncitem.name
    FixtureStore["results_bag"][test_name] = bag
    return bag


@pytest.fixture(scope="session")
def fixture_store(request):
    return FixtureStore


@pytest.fixture(scope="module")
def this_dir(request):
    path = Path(request.fspath)
    return path.parent


# def pytest_collection_modifyitems(session, config, items):
# modified = list()
# synths = dict()
# for i, item in enumerate(items):
# name = item.name
# if "synthesis" in name:
# base_name = name.split("_synthesis")[0]
# synths[base_name] = name
# else:
# modified.append(item)
# import pdb; pdb.set_trace()
# print("Hallo")
