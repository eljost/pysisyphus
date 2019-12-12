#!/usr/bin/env python3

from pathlib import Path

import pytest
from distutils import dir_util
import os

@pytest.fixture
def data_dir(tmpdir, request):
    """As described in https://stackoverflow.com/questions/29627341"""
    print("request", request)
    filename = request.module.__file__
    print("filename", filename)
    test_dir, _ = os.path.splitext(filename)
    print("test_dir", test_dir)

    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, str(tmpdir))

    return Path(tmpdir)
