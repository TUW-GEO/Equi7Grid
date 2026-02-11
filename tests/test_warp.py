import os
import tempfile
from pathlib import Path

import numpy as np
import pyproj
import pytest
from approvaltests.approvals import verify_file
from approvaltests.namer import NamerFactory

from equi7grid import get_standard_equi7grid
from equi7grid._const import WARP_INSTALLED
from equi7grid._core import Equi7Grid
from equi7grid.warp import resample_to_equi7_tiles

if WARP_INSTALLED:
    import rasterio as rio
    from rasterio.transform import from_bounds
    from rasterio.warp import Resampling


@pytest.fixture
def output_dirpath(tmp_path: Path):
    return tmp_path


@pytest.fixture(scope="module")
def e7grid():
    return get_standard_equi7grid({"T6": 1000})


@pytest.fixture(scope="module")
def e7grid_eu():
    return get_standard_equi7grid({"T6": 1000}, continent_order=["EU"])


@pytest.fixture(scope="module")
def e7grid_eu_coarse():
    return get_standard_equi7grid({"T6": 10_000}, continent_order=["EU"])


@pytest.fixture(scope="module")
def lonlat_img_continuous():
    if not WARP_INSTALLED:
        return None
    tmp_dirpath = Path(tempfile.mkdtemp())
    tmp_filepath = tmp_dirpath / "lonlat_img_continuous.tif"
    extent = (50, 60.7, 51.2, 63.2)
    sampling = 0.1
    nodata = -9999
    dtype = "int16"
    compress = "zstd"

    width, height = (
        int((extent[2] - extent[0]) / sampling),
        int((extent[3] - extent[1]) / sampling),
    )
    w_ar, h_ar = np.meshgrid(np.arange(width), np.arange(height))
    ar = w_ar**2 + h_ar**2
    ar[int(height / 2), int(width / 2)] = nodata

    transform = from_bounds(*extent, width, height)
    kwargs = {}
    kwargs.update(
        {
            "crs": pyproj.CRS.from_epsg(4326),
            "transform": transform,
            "width": width,
            "height": height,
            "nodata": nodata,
            "dtype": dtype,
            "compress": compress,
            "count": 1,
        }
    )

    with rio.open(tmp_filepath, "w", **kwargs) as rd:
        rd.write(ar.astype(np.int16), 1)

    return tmp_filepath


@pytest.fixture(scope="module")
def webmercator_img_discrete():
    if not WARP_INSTALLED:
        return None
    tmp_dirpath = Path(tempfile.mkdtemp())
    tmp_filepath = tmp_dirpath / "webmercator_img_discrete.tif"
    extent = (4_900_000, 7_810_000, 6_180_000, 8_750_000)
    sampling = 10_000
    nodata = 255
    dtype = "uint8"
    compress = "zstd"

    width, height = (
        int((extent[2] - extent[0]) / sampling),
        int((extent[3] - extent[1]) / sampling),
    )
    w_ar, h_ar = np.meshgrid(np.arange(width), np.arange(height))
    ar = w_ar**2 + h_ar**2
    ar[ar < (width + height)] = 0
    ar[(ar >= (width + height)) & (ar < (width**2 + height))] = 1
    ar[ar >= (width**2 + height)] = 2
    ar[int(height / 4), int(width / 4)] = nodata

    transform = from_bounds(*extent, width, height)
    kwargs = {}
    kwargs.update(
        {
            "crs": pyproj.CRS.from_epsg(3857),
            "transform": transform,
            "width": width,
            "height": height,
            "nodata": nodata,
            "dtype": dtype,
            "compress": compress,
            "count": 1,
        }
    )

    with rio.open(tmp_filepath, "w", **kwargs) as rd:
        rd.write(ar.astype(np.uint8), 1)

    return tmp_filepath


@pytest.mark.skipif(
    os.name == "nt", reason="CI Windows has troubles creating directories"
)
@pytest.mark.warp_installed
def test_reprojecting_lonlat_continuous(
    lonlat_img_continuous: Path, e7grid: Equi7Grid, output_dirpath: Path
):
    tile_filepaths = resample_to_equi7_tiles(
        lonlat_img_continuous,
        e7grid,
        output_dirpath,
        tiling_id="T6",
        accurate_boundary=True,
        compress_type="zstd",
    )

    n_expected_tiles = 4
    assert len(tile_filepaths) == n_expected_tiles

    verify_file(
        (
            output_dirpath
            / "EQUI7_AS/E018N072T6/lonlat_img_continuous_AS_E018N072T6.tif"
        ).as_posix(),
        options=NamerFactory.with_parameters("AS_E018N072T6"),
    )
    verify_file(
        (
            output_dirpath
            / "EQUI7_AS/E018N066T6/lonlat_img_continuous_AS_E018N066T6.tif"
        ).as_posix(),
        options=NamerFactory.with_parameters("AS_E018N066T6"),
    )
    verify_file(
        (
            output_dirpath
            / "EQUI7_EU/E066N030T6/lonlat_img_continuous_EU_E066N030T6.tif"
        ).as_posix(),
        options=NamerFactory.with_parameters("EU_E066N030T6"),
    )
    verify_file(
        (
            output_dirpath
            / "EQUI7_EU/E072N030T6/lonlat_img_continuous_EU_E072N030T6.tif"
        ).as_posix(),
        options=NamerFactory.with_parameters("EU_E072N030T6"),
    )


@pytest.mark.skipif(
    os.name == "nt", reason="CI Windows has troubles creating directories"
)
@pytest.mark.warp_installed
def test_reprojecting_lonlat_continuous_eu(
    lonlat_img_continuous: Path, e7grid_eu: Equi7Grid, output_dirpath: Path
):
    tile_filepaths = resample_to_equi7_tiles(
        lonlat_img_continuous,
        e7grid_eu,
        output_dirpath,
        tiling_id="T6",
        accurate_boundary=True,
        compress_type="zstd",
    )

    n_expected_tiles = 2
    assert len(tile_filepaths) == n_expected_tiles

    verify_file(
        (
            output_dirpath
            / "EQUI7_EU/E066N030T6/lonlat_img_continuous_EU_E066N030T6.tif"
        ).as_posix(),
        options=NamerFactory.with_parameters("EU_E066N030T6"),
    )
    verify_file(
        (
            output_dirpath
            / "EQUI7_EU/E072N030T6/lonlat_img_continuous_EU_E072N030T6.tif"
        ).as_posix(),
        options=NamerFactory.with_parameters("EU_E072N030T6"),
    )


@pytest.mark.skipif(
    os.name == "nt", reason="CI Windows has troubles creating directories"
)
@pytest.mark.warp_installed
def test_reprojecting_webmercator_discrete_eu(
    webmercator_img_discrete: Path, e7grid_eu_coarse: Equi7Grid, output_dirpath: Path
):
    tile_filepaths = resample_to_equi7_tiles(
        webmercator_img_discrete,
        e7grid_eu_coarse,
        output_dirpath,
        tiling_id="T6",
        accurate_boundary=True,
        compress_type="zstd",
        resampling_type=Resampling.nearest,
    )

    n_expected_tiles = 4
    assert len(tile_filepaths) == n_expected_tiles

    verify_file(
        (
            output_dirpath
            / "EQUI7_EU/E066N030T6/webmercator_img_discrete_EU_E066N030T6.tif"
        ).as_posix(),
        options=NamerFactory.with_parameters("EU_E066N030T6"),
    )
    verify_file(
        (
            output_dirpath
            / "EQUI7_EU/E072N030T6/webmercator_img_discrete_EU_E072N030T6.tif"
        ).as_posix(),
        options=NamerFactory.with_parameters("EU_E072N030T6"),
    )
    verify_file(
        (
            output_dirpath
            / "EQUI7_EU/E066N024T6/webmercator_img_discrete_EU_E066N024T6.tif"
        ).as_posix(),
        options=NamerFactory.with_parameters("EU_E066N024T6"),
    )
    verify_file(
        (
            output_dirpath
            / "EQUI7_EU/E072N024T6/webmercator_img_discrete_EU_E072N024T6.tif"
        ).as_posix(),
        options=NamerFactory.with_parameters("EU_E072N024T6"),
    )


if __name__ == "__main__":
    pass
