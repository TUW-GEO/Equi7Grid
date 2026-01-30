# Copyright (c) 2026, TU Wien
# Licensed under the MIT License. See LICENSE file.

"""Code for resampling image raster data to the Equi7Grid.

Yields tiled raster images organised in folders for
a) in the continental zones + grid sampling
b) the subgrid tiling

The module should help to easily bring your raster images to the Equi7Grid
spatial reference.
"""

import multiprocessing as mp
import tempfile
from collections.abc import Callable, Mapping
from functools import partial
from pathlib import Path
from typing import Any, cast

import numpy as np
import shapely
from pytileproj.projgeom import ProjGeom
from shapely.geometry import shape

from equi7grid._const import WARP_INSTALLED
from equi7grid._core import Equi7Grid, Equi7Tile
from equi7grid._types import Extent

if WARP_INSTALLED:
    import rasterio
    import rasterio.features
    from rasterio.warp import Affine, Resampling, reproject
    from scipy import ndimage


def requires_warp(f: Callable) -> Callable:
    """Check if the warp extension is installed."""

    def wrapper(*args: Any, **kwargs: Any) -> Any:  # noqa: ANN401
        if not WARP_INSTALLED:
            err_msg = "It is required to install the 'warp' extension."
            raise ImportError(err_msg)
        return f(*args, **kwargs)

    return wrapper


@requires_warp
def pixel_to_world_coords(tf: "Affine", pixel_coords: np.ndarray) -> np.ndarray:
    """Convert pixel to world system coordinates."""
    sa, sb, sc, sd, se, sf, _, _, _ = tf
    world_coords = np.zeros_like(pixel_coords, dtype=np.float64)
    world_coords[:, 0] = pixel_coords[:, 0] * sa + pixel_coords[:, 1] * sb + sc
    world_coords[:, 1] = pixel_coords[:, 0] * sd + pixel_coords[:, 1] * se + sf

    return world_coords


@requires_warp
def get_raster_boundary(filepath: Path) -> ProjGeom:
    """Get accurate boundary of a raster file."""
    qlook_size = 400.0

    with (
        tempfile.NamedTemporaryFile(suffix=filepath.suffix) as tmp_filepath,
        rasterio.open(filepath) as src,
    ):
        src_transform = src.transform
        dst_transform = src_transform * Affine.scale(qlook_size)
        dst_width = src.width / qlook_size
        dst_height = src.height / qlook_size

        kwargs = src.meta.copy()
        kwargs.update(
            {
                "crs": src.crs,
                "transform": dst_transform,
                "width": dst_width,
                "height": dst_height,
            }
        )

        with rasterio.open(tmp_filepath, "w", **kwargs) as dst:
            reproject(
                source=rasterio.band(src, 1),
                destination=rasterio.band(dst, 1),
                src_transform=src_transform,
                src_crs=src.crs,
                dst_transform=dst_transform,
                dst_crs=src.crs,
                resampling=Resampling.nearest,
            )

            data = dst.read(1)
            mask = data != dst.nodata

            # morphologic dilation
            iterations = 3
            struct = ndimage.generate_binary_structure(2, 2)
            new_mask = ndimage.binary_dilation(
                mask, iterations=iterations, structure=struct
            )
            roi_polys = []
            for geom in rasterio.features.shapes(new_mask.astype(np.uint8)):
                if geom[1] == 1:
                    roi_poly = shapely.transform(
                        shape(geom[0]),
                        lambda x: pixel_to_world_coords(dst_transform, x),
                    )
                    roi_polys.append(roi_poly)

            if len(roi_polys) == 1:
                roi_poly = roi_polys[0]
            else:
                roi_poly = shapely.MultiPolygon(roi_polys)

            return ProjGeom(geom=roi_poly, crs=src.crs)


@requires_warp
def get_raster_extent(filepath: Path) -> ProjGeom:
    """Get extent of a raster file."""
    with rasterio.open(filepath) as src:
        min_x, min_y, max_x, max_y = src.bounds
        extent_poly = shapely.Polygon(
            [(min_x, min_y), (max_x, min_y), (max_x, max_y), (min_x, max_y)]
        )
        extent_crs = src.crs

    return ProjGeom(geom=extent_poly, crs=extent_crs)


def get_default_e7_filename(filepath: Path, ftilename: str) -> str:
    """Get default Equi7Grid filename."""
    return f"{filepath.stem}_{ftilename}{filepath.suffix}"


@requires_warp
def resample_tile(  # noqa: PLR0913
    e7tile: Equi7Tile,
    filepath: Path,
    output_dirpath: Path,
    *,
    band: int = 1,
    image_nodata: float | None = None,
    resampling_type: "Resampling" = Resampling.bilinear,
    compress_type: str = "LZW",
    naming_traffo: Callable | None = None,
    tile_nodata: float | None = None,
    tile_dtype: np.dtype | None = None,
    tile_scale: float | None = None,
    tile_offset: float | None = None,
    tile_blocksize: int | None = None,
    tif_is_tiled: bool = True,
    overwrite: bool = False,
    create_e7_folder: bool = True,
) -> Path:
    """Resample a portion of an image to an Equi7Grid tile."""
    ftilename = cast("str", e7tile.name)

    if create_e7_folder:
        grid_foldername = f"EQUI7_{ftilename[0:6]}"
        tile_dirpath = output_dirpath / grid_foldername / ftilename[7:]
        tile_dirpath.mkdir(exist_ok=True, parents=True)
    else:
        tile_dirpath = output_dirpath

    if naming_traffo is None:
        tile_filename = get_default_e7_filename(filepath, ftilename)
    else:
        tile_filename = naming_traffo(filepath, ftilename)

    tile_filepath = tile_dirpath / tile_filename

    with rasterio.open(filepath, nodata=image_nodata) as src:
        transform = Affine.from_gdal(*e7tile.geotrans)
        blockxsize, blockysize = src.block_shapes[band - 1]
        kwargs = src.meta.copy()
        kwargs.update(
            {
                "overwrite": overwrite,
                "crs": e7tile.pyproj_crs,
                "transform": transform,
                "width": e7tile.width,
                "height": e7tile.height,
                "nodata": tile_nodata or src.nodatavals[band - 1],
                "dtype": tile_dtype or src.dtypes[band - 1],
                "compress": compress_type,
                "tiled": tif_is_tiled,
                "blockxsize": tile_blocksize or blockxsize,
                "blockysize": tile_blocksize or blockysize,
            }
        )
        with rasterio.open(tile_filepath, "w", **kwargs) as dst:
            reproject(
                source=rasterio.band(src, band),
                destination=rasterio.band(dst, band),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=dst.transform,
                dst_crs=dst.crs,
                resampling=resampling_type,
            )

            dst._set_all_scales([tile_scale or src.scales[band - 1]])  # noqa: SLF001
            dst._set_all_offsets([tile_offset or src.offsets[band - 1]])  # noqa: SLF001
            dst.update_tags(**src.tags())
            dst.update_stats()

    return tile_filepath


def get_overlapping_tiles(  # noqa: PLR0913
    e7grid: Equi7Grid,
    *,
    tiling_id: str | int = 0,
    xy_bbox_map: Mapping[str, Extent] | None = None,
    geog_bbox: Extent | None = None,
    roi_geom: ProjGeom | None = None,
    filepath: Path | None = None,
    cover_land: bool = False,
    accurate_boundary: bool = False,
) -> list[Equi7Tile]:
    """Get all overlapping Equi7Grid tiles."""
    e7tiles = []
    if xy_bbox_map is not None:
        for continent, xy_bbox in xy_bbox_map.items():
            e7tiles_cont = e7grid[continent].get_tiles_in_bbox(
                xy_bbox, tiling_id=tiling_id
            )
            e7tiles.extend(list(e7tiles_cont))
    elif geog_bbox is not None:
        e7tiles = e7grid.get_tiles_in_geog_bbox(
            geog_bbox, tiling_id, cover_land=cover_land
        )
    elif roi_geom is not None:
        e7tiles = e7grid.get_tiles_in_geom(roi_geom, tiling_id, cover_land=cover_land)
    elif filepath is not None:
        if accurate_boundary:
            roi_geom = get_raster_boundary(filepath)
        else:
            roi_geom = get_raster_extent(filepath)
        e7tiles = e7grid.get_tiles_in_geom(roi_geom, tiling_id, cover_land=cover_land)
    else:
        err_msg = (
            "One of the arguments must be set: 'xy_bbox_map', 'geog_bbox', "
            "'roi_geom', or 'filepath'."
        )
        raise ValueError(err_msg)

    return list(e7tiles)


def ftilenames_to_tiles(ftilenames: list[str], e7grid: Equi7Grid) -> list[Equi7Tile]:
    """Convert Equi7Grid tilenames to tiles."""
    return [e7grid.get_tile_from_name(ftilename) for ftilename in ftilenames]


@requires_warp
def resample_to_equi7_tiles(  # noqa: PLR0913
    filepath: Path,
    e7grid: Equi7Grid,
    output_dirpath: Path,
    *,
    tiling_id: str | int = 1,
    xy_bbox_map: Mapping[str, Extent] | None = None,
    geog_bbox: Extent | None = None,
    roi_geom: ProjGeom | None = None,
    cover_land: bool = False,
    accurate_boundary: bool = False,
    ftilenames: list[str] | None = None,
    band: int = 1,
    image_nodata: float | None = None,
    resampling_type: "Resampling" = Resampling.bilinear,
    compress_type: str = "LZW",
    naming_traffo: Callable | None = None,
    tile_nodata: float | None = None,
    tile_dtype: np.dtype | None = None,
    tile_scale: float | None = None,
    tile_offset: float | None = None,
    tile_blocksize: int | None = None,
    tif_is_tiled: bool = True,
    overwrite: bool = False,
    create_e7_folder: bool = True,
    n_tasks: int = 1,
) -> list[Path]:
    """Resample an image to Equi7Grid tiles."""
    if ftilenames is not None:
        e7tiles = ftilenames_to_tiles(ftilenames, e7grid)
    else:
        e7tiles = list(
            get_overlapping_tiles(
                e7grid,
                tiling_id=tiling_id,
                xy_bbox_map=xy_bbox_map,
                geog_bbox=geog_bbox,
                roi_geom=roi_geom,
                filepath=filepath,
                cover_land=cover_land,
                accurate_boundary=accurate_boundary,
            )
        )

    resample_kwargs = {
        "filepath": filepath,
        "output_dirpath": output_dirpath,
        "band": band,
        "image_nodata": image_nodata,
        "resampling_type": resampling_type,
        "compress_type": compress_type,
        "naming_traffo": naming_traffo,
        "tile_nodata": tile_nodata,
        "tile_dtype": tile_dtype,
        "tile_scale": tile_scale,
        "tile_offset": tile_offset,
        "tile_blocksize": tile_blocksize,
        "tif_is_tiled": tif_is_tiled,
        "overwrite": overwrite,
        "create_e7_folder": create_e7_folder,
    }

    if n_tasks == 1:
        tile_filepaths = []
        for e7tile in e7tiles:
            tile_filepath = resample_tile(e7tile, **resample_kwargs)  # ty: ignore[invalid-argument-type]
            tile_filepaths.append(tile_filepath)
    else:
        n_tiles = len(e7tiles)
        num_cpu = min(n_tasks, n_tiles)
        with mp.Pool(processes=num_cpu) as pool:
            func = partial(resample_tile, **resample_kwargs)  # ty: ignore[invalid-argument-type]
            tile_filepaths = pool.map(func, e7tiles)
            pool.close()
            pool.join()

        tile_filepaths = [x for x in tile_filepaths if x]

    return tile_filepaths


if __name__ == "__main__":
    pass
