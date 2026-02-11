import sys

import pytest


@pytest.fixture(autouse=True)
def warp_installed(request):  # noqa: ANN001
    if request.node.get_closest_marker("warp_installed"):
        vis_mods_installed = "rasterio" in sys.modules and "scipy" in sys.modules
        if not vis_mods_installed:
            pytest.skip()
