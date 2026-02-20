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
from typing import Any, Union, cast

import numpy as np
import shapely
from pytileproj import ProjGeom
from shapely.geometry import shape

from equi7grid._const import WARP_INSTALLED
from equi7grid._types import Extent
from equi7grid.core import Equi7Grid, Equi7Tile

if WARP_INSTALLED:
    import rasterio
    import rasterio.features
    from rasterio.warp import Affine, Resampling, reproject
    from scipy import ndimage


def _requires_warp(f: Callable) -> Callable:
    """Check if the warp extension is installed."""

    def wrapper(*args: Any, **kwargs: Any) -> Any:  # noqa: ANN401
        if not WARP_INSTALLED:
            err_msg = "It is required to install the 'warp' extension."
            raise ImportError(err_msg)
        return f(*args, **kwargs)

    return wrapper


@_requires_warp
def pixel_to_world_coords(tf: "Affine", pixel_coords: np.ndarray) -> np.ndarray:
    """Convert pixel to world system coordinates.

    Parameters
    ----------
    tf: Affine
        Affine object representing affine transformation parameters.
    pixel_coords: np.ndarray
        2D array with shape (n, 2) containing pixel coordinates:
            - first column are pixel column coordinates
            - second column are pixel row coordinates

    Returns
    -------
    np.ndarray
        World system coordinates array with shape (n, 2):
            - first column are X coordinates
            - second column are Y coordinates

    """
    sa, sb, sc, sd, se, sf, _, _, _ = tf
    world_coords = np.zeros_like(pixel_coords, dtype=np.float64)
    world_coords[:, 0] = pixel_coords[:, 0] * sa + pixel_coords[:, 1] * sb + sc
    world_coords[:, 1] = pixel_coords[:, 0] * sd + pixel_coords[:, 1] * se + sf

    return world_coords


@_requires_warp
def get_raster_boundary(filepath: Path) -> ProjGeom:
    """Get accurate boundary of a raster file.

    Parameters
    ----------
    filepath: Path
        Geospatial filepath to retrieve raster boundary from.

    Returns
    -------
    ProjGeom
        Accurate raster boundary represented by ProjGeom instance
        (shapely polygon + pyproj.CRS).

    """
    qlook_size = 400.0

    with (
        tempfile.NamedTemporaryFile(suffix=filepath.suffix) as tmp_filepath,
        rasterio.open(filepath) as src,
    ):
        src_transform = src.transform
        dst_width = src.width / qlook_size
        dst_height = src.height / qlook_size
        scale = [qlook_size, qlook_size]
        if dst_width < 1:
            dst_width = src.width
            scale[0] = 1
        if dst_height < 1:
            dst_height = src.height
            scale[1] = 1
        dst_transform = src_transform * Affine.scale(*scale)

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


@_requires_warp
def get_raster_extent(filepath: Path) -> ProjGeom:
    """Get extent of a raster file.

    Parameters
    ----------
    filepath: Path
        Geospatial filepath to retrieve raster extent from.

    Returns
    -------
    ProjGeom
        Raster extent represented by ProjGeom instance
        (shapely polygon + pyproj.CRS).

    """
    with rasterio.open(filepath) as src:
        min_x, min_y, max_x, max_y = src.bounds
        extent_poly = shapely.Polygon(
            [(min_x, min_y), (max_x, min_y), (max_x, max_y), (min_x, max_y)]
        )
        extent_crs = src.crs

    return ProjGeom(geom=extent_poly, crs=extent_crs)


def get_default_e7_filename(filepath: Path, ftilename: str) -> str:
    """Get default Equi7Grid filename.

    Parameters
    ----------
    filepath: Path
        Full system path to geospatial image file.
    ftilename: str
        Full Equi7Grid tilename.

    Returns
    -------
    str
        Equi7Grid standard filename.

    """
    return f"{filepath.stem}_{ftilename}{filepath.suffix}"


@_requires_warp
def resample_tile(  # noqa: PLR0913
    e7tile: Equi7Tile,
    filepath: Path,
    output_dirpath: Path,
    *,
    band: int = 1,
    image_nodata: float | None = None,
    resampling_type: Union["Resampling", None] = None,
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
    """Resample a portion of an image to an Equi7Grid tile.

    Parameters
    ----------
    e7tile: Equi7Tile
        Equi7Tile instance to resample and slice to.
    filepath: Path
        Full system path to geospatial image file.
    output_dirpath: Path
        Full system path to the output directory.
    band: int, optional
        Band number (defaults to 1).
    image_nodata: float | None, optional
        Nodata value of geospatial image.
        Defaults to nodata value provided in the file metadata.
    resampling_type: Resampling | None, optional
        Resampling method provided by rasterio. Defaults to neareast.
    compress_type: str, optional
        Compression type (defaults to "LZW").
    naming_traffo: Callable | None, optional
        Callable/function to define naming convention.
        It expects two input arguments:
            - `filepath`: the full system path to geospatial image file
            - `ftilename`: the full tilename.
        Defaults to standard naming: "[filename]_[full tilename].[suffix]"
    tile_nodata: float | None, optional
        Nodata value of the tiled image. Defaults to nodata value of the input image.
    tile_dtype: np.dtype | None, optional
        Data type of the tiled image. Defaults to data type of the input image.
    tile_scale: float | None, optional
        Scale factor of the tiled image. Defaults to scale factor of the input image.
    tile_offset: float | None, optional
        Offset of the tiled image. Defaults to offset of the input image.
    tile_blocksize: int | None, optional
        Block size of the tiled image. Defaults to block size of the input image.
    tif_is_tiled: bool, optional
        Flag to define if the output tile image should be tiled w.r.t. the GeoTIFF
        format (defined by the block size) or not (defaults to True).
    overwrite: bool, optional
        Flag if data should be overwritten or not (defaults to false).
    create_e7_folder: bool, optional
        Flag if a root Equi7Grid folder ("EQUI7_[continent]") should be
        created or not (defaults to true).

    Returns
    -------
    Path
        Full system path to the image tile in the Equi7Grid.

    """
    ftilename = cast("str", e7tile.name)

    if create_e7_folder:
        grid_foldername = f"EQUI7_{ftilename[0:2]}"
        tile_dirpath = output_dirpath / grid_foldername / ftilename[3:]
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
        if (blockxsize % 16) != 0:
            blockxsize = None
        if (blockysize % 16) != 0:
            blockysize = None
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
        resampling_type = resampling_type or Resampling.nearest
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
    """Get all overlapping Equi7Grid tiles.

    Parameters
    ----------
    e7grid: Equi7Grid
        Equi7Grid instance to get the tiles from.
    tiling_id: int | str, optional
        Tiling level or name.
        Defaults to 0, the first tiling level.
    xy_bbox_map: Mapping[str, Extent] | None, optional
        Defines projected bounding boxes in the Equi7Grid
        to consider for reprojection. The keys are the continents,
        and the values the spatial extent/bounding box.
        Defaults to image extent.
    geog_bbox: Extent | None, optional
        Defines geographic bounding box to consider for reprojection.
        Defaults to image extent.
    roi_geom: ProjGeom | None, optional
        Defines geospatial geometry to consider only a certain region
        for reprojection. Defaults to image extent.
    filepath: Path | None, optional
        Full system path to geospatial image file.
    cover_land: bool, optional
        Flag defining if only tiles covering land should be considered
        (defaults to false).
    accurate_boundary: bool, optional
        Flag to define if the accurate raster boundary
        (excluding nodata values) should be used.
        If false, then the image extent will be used (default).

    """
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
    """Convert Equi7Grid tilenames to tiles.

    Parameters
    ----------
    ftilenames: list[str]
        List of full Equi7Grid tilenames.
    e7grid: Equi7Grid
        Equi7Grid instance, where the given tilenames are part of.

    Returns
    -------
    list[Equi7Tile]
        List of Equi7Grid tiles.

    """
    return [e7grid.get_tile_from_name(ftilename) for ftilename in ftilenames]


@_requires_warp
def resample_to_equi7_tiles(  # noqa: PLR0913
    filepath: Path,
    e7grid: Equi7Grid,
    output_dirpath: Path,
    *,
    tiling_id: str | int = 0,
    xy_bbox_map: Mapping[str, Extent] | None = None,
    geog_bbox: Extent | None = None,
    roi_geom: ProjGeom | None = None,
    cover_land: bool = False,
    accurate_boundary: bool = False,
    ftilenames: list[str] | None = None,
    band: int = 1,
    image_nodata: float | None = None,
    resampling_type: Union["Resampling", None] = None,
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
    """Resample an image to Equi7Grid tiles.

    Parameters
    ----------
    filepath: Path
        Full system path to geospatial image file.
    e7grid: Equi7Grid
        Equi7Grid instance to resample and slice to.
    output_dirpath: Path
        Full system path to the output directory.
    tiling_id: int | str, optional
        Tiling level or name.
        Defaults to 0, the first tiling level.
    xy_bbox_map: Mapping[str, Extent] | None, optional
        Defines projected bounding boxes in the Equi7Grid to consider
        for reprojection. The keys are the continents, and the values
        the spatial extent/bounding box. Defaults to image extent.
    geog_bbox: Extent | None, optional
        Defines geographic bounding box to consider for reprojection.
        Defaults to image extent.
    roi_geom: ProjGeom | None, optional
        Defines geospatial geometry to consider only a certain region for
        reprojection.Defaults to image extent.
    cover_land: bool, optional
        Flag defining if only tiles covering land should be considered
        (defaults to false).
    accurate_boundary: bool, optional
        Flag to define if the accurate raster boundary
        (excluding nodata values) should be used.
        If false, then the image extent will be used (default).
    ftilenames: list[str] | None, optional
        List of full tilenames to resample to.
        Defaults to all tiles intersecting with the image.
    band: int, optional
        Band number (defaults to 1).
    image_nodata: float | None, optional
        Nodata value of geospatial image.
        Defaults to nodata value provided in the file metadata.
    resampling_type: Resampling | None, optional
        Resampling method provided by rasterio. Defaults to neareast.
    compress_type: str, optional
        Compression type (defaults to "LZW").
    naming_traffo: Callable | None, optional
        Callable/function to define naming convention.
        It expects two input arguments:
            - `filepath`: the full system path to geospatial image file
            - `ftilename`: the full tilename.
        Defaults to standard naming: "[filename]_[full tilename].[suffix]"
    tile_nodata: float | None, optional
        Nodata value of the tiled image. Defaults to nodata value of the input image.
    tile_dtype: np.dtype | None, optional
        Data type of the tiled image. Defaults to data type of the input image.
    tile_scale: float | None, optional
        Scale factor of the tiled image. Defaults to scale factor of the input image.
    tile_offset: float | None, optional
        Offset of the tiled image. Defaults to offset of the input image.
    tile_blocksize: int | None, optional
        Block size of the tiled image. Defaults to block size of the input image.
    tif_is_tiled: bool, optional
        Flag to define if the output tile image should be tiled w.r.t. the GeoTIFF
        format (defined by the block size) or not (defaults to True).
    overwrite: bool, optional
        Flag if data should be overwritten or not (defaults to false).
    create_e7_folder: bool, optional
        Flag if a root Equi7Grid folder ("EQUI7_[continent]") should be created
        or not (defaults to true).
    n_tasks: int, optional
        Number of parallel processes to resample multiple tiles in parallel
        (defaults to 1).

    Returns
    -------
    list[Path]
        List of full system paths to the image tiles in the Equi7Grid.

    """
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
    resampling_type = resampling_type or Resampling.nearest

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
