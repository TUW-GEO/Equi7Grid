# Copyright (c) 2026, TU Wien
# Licensed under the MIT License. See LICENSE file.

from importlib.util import find_spec

KM = 1e3  # kilometers in meters
MAX_SAMPLING = 10 * KM  # maximum sampling in meters

WARP_INSTALLED = None not in [find_spec(pkg) for pkg in ["rasterio", "scipy"]]
