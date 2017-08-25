#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from equi7grid.skeleton import fib

__author__ = "Bernhard Bauer-Marschallinger"
__copyright__ = "Bernhard Bauer-Marschallinger"
__license__ = "none"


def test_fib():
    assert fib(1) == 1
    assert fib(2) == 1
    assert fib(7) == 13
    with pytest.raises(AssertionError):
        fib(-10)
