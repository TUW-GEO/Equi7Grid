from pathlib import Path

import pytest


@pytest.fixture
def input_dir():
    return Path(__file__).parent / "prj_accuracy_test"