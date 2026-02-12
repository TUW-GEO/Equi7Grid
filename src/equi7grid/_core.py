# Copyright (c) 2026, TU Wien
# Licensed under the MIT License. See LICENSE file.

"""Core module defining the Equi7Grid classes."""

import warnings
from collections.abc import Generator, Mapping
from pathlib import Path
from typing import Any

import numpy as np
import shapely
from geographiclib.geodesic import Geodesic
from morecantile.models import Tile as RegularTile
from pytileproj import (
    GeogGeom,
    ProjGeom,
    ProjSystemDefinition,
    RasterTile,
    RegularGrid,
    RegularProjTilingSystem,
    RegularTilingDefinition,
    TileOutOfZoneError,
)
from pytileproj._const import DEF_SEG_LEN_DEG
from pytileproj._types import SamplingFloatOrMap
from pytileproj.projgeom import (
    convert_any_to_geog_geom,
    transform_geometry,
)

from equi7grid._types import Extent, T_co
from equi7grid.create_grids import get_standard_tilings, get_system_definitions

Equi7TileGenerator = Generator["Equi7Tile", None, None]


class Equi7Tile(RasterTile[Any]):
    """Defines a tile in the Equi7Grid."""

    covers_land: bool


class Equi7TilingSystem(RegularProjTilingSystem):
    """Defines a tiling system for each Equi7Grid continent."""

    land_zone: ProjGeom

    def _extent_covers_land(self, extent: Extent) -> bool:
        """Evaluate if the given extent covers land masses.

        Parameters
        ----------
        extent: Extent
            A given spatial extent (x_min, y_min, x_max, y_max).

        Returns
        -------
        bool:
            True if the extent intersects land masses.

        """
        ll_x, ll_y, ur_x, ur_y = extent
        poly = shapely.Polygon(((ll_x, ll_y), (ll_x, ur_y), (ur_x, ur_y), (ur_x, ll_y)))

        return shapely.intersects(poly, self.land_zone.geom)

    def _tile_to_raster_tile(
        self, tile: RegularTile, name: str | None = None
    ) -> Equi7Tile:
        """Create an Equi7 tile object from a given regular tile.

        Parameters
        ----------
        tile: RegularTile
            Regular tile object from TMS.
        name: str | None, optional
            Tilename.

        Returns
        -------
        Equi7Tile
            Equi7 tile object.

        """
        extent = self._tms.xy_bounds(tile)
        sampling = self[tile.z].sampling
        covers_land = self._extent_covers_land(extent)
        return Equi7Tile.from_extent(
            extent,
            self.pyproj_crs,
            sampling,
            sampling,
            name=name,
            covers_land=covers_land,
        )

    def _tile_to_partial_name(self, tile: RegularTile) -> str:
        """Create a tilename from a given regular tile object.

        The generated tilename is unique within a tiling system.

        Parameters
        ----------
        tile: RegularTile
            Regular tile object from TMS.

        Returns
        -------
        str
            Partial Equi7 tilename (unique at tiling system level).

        """
        e7_tile = self._tile_to_raster_tile(tile)
        ll_x, ll_y, _, _ = e7_tile.outer_boundary_extent
        tile_x, tile_y = int(ll_x / 1e5), int(ll_y / 1e5)
        y_label = "S" if tile_y < 0 else "N"
        if tile_y < 0:
            pass

        return f"E{tile_x:03}{y_label}{abs(tile_y):03}T{self.tilings[tile.z].name[1:]}"

    def _tile_to_full_name(self, tile: RegularTile) -> str:
        """Create a tilename from a given regular tile object.

        The generated tilename is unique within the grid.

        Parameters
        ----------
        tile: RegularTile
            Regular tile object from TMS.

        Returns
        -------
        str
            Full Equi7 tilename (unique at grid level).

        """
        tilename = self._tile_to_partial_name(tile)
        return f"{self.name}_{tilename}"

    def _tile_to_name(self, tile: RegularTile) -> str:
        """Create a full Equi7 tilename.

        The generated tilename is unique within
        the whole grid.

        Parameters
        ----------
        tile: RegularTile
            Regular tile object from TMS.

        Returns
        -------
        str
            Equi7 tilename.

        """
        return self._tile_to_full_name(tile)

    def _name_to_tile(self, tilename: str) -> RegularTile:
        """Get regular tile from a tilename.

        Parameters
        ----------
        tilename: str
            Equi7 tile name.

        Returns
        -------
        RegularTile
            Regular tile object in the TMS corresponding
            to the given tile name.

        """
        tiling_id = tilename[-2:]
        tiling_level = self.tiling_id_to_level(tiling_id)

        y_label = tilename[4]
        y_mirror_fac = 1 if y_label == "N" else -1
        x, y = float(tilename[1:4]) * 1e5, float(tilename[5:8]) * 1e5 * y_mirror_fac

        tilesize = self[tiling_level].tile_shape[0]
        sampling = self[tiling_level].sampling
        tile = RasterTile.from_extent(
            (x, y, x + tilesize, y + tilesize), self.pyproj_crs, sampling, sampling
        )
        if not self._tile_in_zone(tile):
            raise TileOutOfZoneError(tile)

        return self._tms._tile(x + tilesize / 2.0, y + tilesize / 2.0, tiling_level)  # noqa: SLF001

    def get_tile_from_name(self, ftilename: str) -> Equi7Tile:
        """Get Equi7 tile from a full tilename.

        Parameters
        ----------
        ftilename: str
            Full Equi7 tile name.

        Returns
        -------
        Equi7Tile
            Equi7 tile corresponding to the given tile name.

        """
        tile = self._name_to_tile(ftilename.split("_")[1])
        e7_tile = self._tile_to_raster_tile(tile, name=ftilename)
        if not self._tile_in_zone(e7_tile):
            raise TileOutOfZoneError(e7_tile)

        return e7_tile

    def get_tiles_in_geog_bbox(
        self,
        bbox: tuple[float, float, float, float],
        tiling_id: int | str,
        *,
        cover_land: bool = False,
    ) -> Equi7TileGenerator:
        """Get all Equi7 tiles intersecting with the geographic bounding box.

        Parameters
        ----------
        bbox: tuple[float, float, float, float]
            Bounding box (x_min, y_min, x_max, y_max) for selecting tiles.
        tiling_id: int | str
            Tiling level or name.
            Defaults to the first tiling level.
        cover_land: bool, optional
            True if only tiles which cover land should be returned.
            Defaults to false.

        Returns
        -------
        Equi7TileGenerator
            Yields Equi7 tile after tile, which intersects with the given
            bounding box.

        """
        for tile in super().get_tiles_in_geog_bbox(bbox, tiling_id=tiling_id):
            if cover_land and not tile.covers_land:
                continue

            yield tile

    def get_tiles_in_geom(
        self, proj_geom: ProjGeom, tiling_id: int | str, *, cover_land: bool = False
    ) -> Equi7TileGenerator:
        """Get all Equi7 tiles intersecting with the geographic bounding box.

        Parameters
        ----------
        proj_geom : ProjGeom
            Projected geometry representing the region of interest.
        tiling_id: int | str
            Tiling level or name.
            Defaults to the first tiling level.
        cover_land: bool, optional
            True if only tiles which cover land should be returned.
            Defaults to false.

        Returns
        -------
        Equi7TileGenerator
            Yields Equi7 tile after tile, which intersects with the given
            bounding box.

        """
        for tile in super().get_tiles_in_geom(proj_geom, tiling_id=tiling_id):
            if cover_land and not tile.covers_land:
                continue

            yield tile

    def calc_length_distortion(
        self, x: float | np.ndarray, y: float | np.ndarray
    ) -> float | np.ndarray:
        """Compute local maximum length distortion k.

        k equals the local areal distortion (as always h=1 for the Azimuthal
        Equidistant projection). Uses planar coordinates thus being much faster,
        than 'calc_length_distortion_on_ellipsoid'.

        Parameters
        ----------
        x: float | np.ndarray
            World system coordinates in X direction.
        y: float | np.ndarray
            World system coordinates in Y direction.

        Returns
        -------
        k: float | np.ndarray
            Local max length distortion = local areal distortion.

        """
        ellaxis = Geodesic.WGS84.a

        # ignore loss of information during CRS to PROJ4 conversion
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            proj4_dict = self.pyproj_crs.to_dict()

        fe, fn = proj4_dict["x_0"], proj4_dict["y_0"]
        dists = np.sqrt((np.array(x) - fe) ** 2 + (np.array(y) - fn) ** 2)

        return dists / ellaxis / np.sin(dists / ellaxis)

    def get_children_from_name(self, tilename: str) -> Equi7TileGenerator:
        """Get all child tiles (next higher tiling level).

        Parameters
        ----------
        tilename: str
            Tilename.

        Returns
        -------
        Equi7TileGenerator
            Yields all tile children as Equi7 tiles.

        """
        yield from super().get_children_from_name(tilename.split("_")[1])

    def get_parent_from_name(self, tilename: str) -> Equi7Tile:
        """Get parent tile (next lower tiling level).

        Parameters
        ----------
        tilename: str
            Tilename.

        Returns
        -------
        Equi7Tile
            Parent Equi7 tile.

        """
        return super().get_parent_from_name(tilename.split("_")[1])


class Equi7Grid(RegularGrid[T_co]):
    """Defines Equi7Grid with all sub-grid."""

    AF: Equi7TilingSystem | None = None
    AN: Equi7TilingSystem | None = None
    AS: Equi7TilingSystem | None = None
    EU: Equi7TilingSystem | None = None
    NA: Equi7TilingSystem | None = None
    OC: Equi7TilingSystem | None = None
    SA: Equi7TilingSystem | None = None

    _rpts_cls = Equi7TilingSystem

    @staticmethod
    def _create_rpts_from_def(
        proj_def: ProjSystemDefinition,
        sampling: SamplingFloatOrMap,
        tiling_defs: dict[int, RegularTilingDefinition],
    ) -> Equi7TilingSystem:
        """Create regular projected tiling system from Equi7 grid definitions.

        Create a regular, projected tiling system instance from given Equi7
        tiling system definitions and a grid sampling.

        Parameters
        ----------
        proj_def: ProjSystemDefinition
            Projection system definition (stores name, CRS, extent,
            and axis orientation).
        sampling: float | int | Dict[int | str, float | int]
            Grid sampling/pixel size specified as a single value or a dictionary with
            tiling IDs as keys and samplings as values.
        tiling_defs: Dict[int, RegularTilingDefinition]
            Tiling definition (stores name/tiling level and tile size).

        Returns
        -------
        Equi7TilingSystem
            Equi7 tiling system instance.

        """
        land_zone_geog = convert_any_to_geog_geom(
            Path(__file__).parent / "data" / "land.parquet"
        )
        land_proj_zone_geog = shapely.intersection(
            proj_def.proj_zone_geog.geom, land_zone_geog.geom
        )
        land_zone = transform_geometry(
            GeogGeom(geom=land_proj_zone_geog), proj_def.crs, segment=DEF_SEG_LEN_DEG
        )
        land_zone.geom = shapely.buffer(land_zone.geom, 0)
        return Equi7TilingSystem.from_sampling(
            sampling, proj_def, tiling_defs, land_zone=land_zone
        )

    def calc_length_distortion_on_ellipsoid(self, lon: float, lat: float) -> float:
        """Compute local maximum length distortion k on the ellipsoid.

        k equals the local areal distortion (as always h=1 for the Azimuthal
        Equidistant projection)

        Parameters
        ----------
        lon: float
            Longitude.
        lat: float
            Latitude.

        Returns
        -------
        k: float
            Local max length distortion = local areal distortion.

        """
        rpts = self.get_systems_from_lonlat(lon, lat)[0]

        # ignore loss of information during CRS to PROJ4 conversion
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            proj4_dict = rpts.pyproj_crs.to_dict()

        lon_0, lat_0 = proj4_dict["lon_0"], proj4_dict["lat_0"]
        # get spherical distance and azimuth between projection
        # centre and point of interest
        geod = Geodesic.WGS84
        gi = geod.Inverse(lat_0, lon_0, lat, lon)
        c1 = gi["s12"]

        # apply equation for distortion in direction
        # perpendicular to the radius, k:
        return c1 / geod.a / np.sin(c1 / geod.a)

    def get_tiles_in_geog_bbox(
        self,
        bbox: tuple[float, float, float, float],
        tiling_id: int | str,
        *,
        cover_land: bool = False,
    ) -> Equi7TileGenerator:
        """Get all Equi7 tiles intersecting with the geographic bounding box.

        Parameters
        ----------
        bbox: tuple[float, float, float, float]
            Bounding box (x_min, y_min, x_max, y_max) for selecting tiles.
        tiling_id: int | str
            Tiling level or name.
            Defaults to the first tiling level.
        cover_land: bool, optional
            True if only tiles which cover land should be returned.
            Defaults to false.

        Returns
        -------
        Equi7TileGenerator
            Yields Equi7 tile after tile, which intersects with the given
            bounding box.

        """
        continents = self.system_order or []
        for continent in continents:
            e7_tiling_sys = self[continent]
            yield from e7_tiling_sys.get_tiles_in_geog_bbox(
                bbox, tiling_id, cover_land=cover_land
            )

    def get_tiles_in_geom(
        self, proj_geom: ProjGeom, tiling_id: int | str, *, cover_land: bool = False
    ) -> Equi7TileGenerator:
        """Get all Equi7 tiles intersecting with the geographic bounding box.

        Parameters
        ----------
        proj_geom : ProjGeom
            Projected geometry representing the region of interest.
        tiling_id: int | str
            Tiling level or name.
            Defaults to the first tiling level.
        cover_land: bool, optional
            True if only tiles which cover land should be returned.
            Defaults to false.

        Returns
        -------
        Equi7TileGenerator
            Yields Equi7 tile after tile, which intersects with the given
            bounding box.

        """
        continents = self.system_order or []
        for continent in continents:
            e7_tiling_sys = self[continent]
            yield from e7_tiling_sys.get_tiles_in_geom(
                proj_geom, tiling_id, cover_land=cover_land
            )

    def get_tile_from_name(self, ftilename: str) -> Equi7Tile:
        """Get Equi7 tile from a full tilename.

        Parameters
        ----------
        ftilename: str
            Full Equi7 tile name.

        Returns
        -------
        Equi7Tile
            Equi7 tile corresponding to the given tile name.

        """
        continent = ftilename[:2]
        return self[continent].get_tile_from_name(ftilename)


def get_user_equi7grid(
    sampling: SamplingFloatOrMap,
    tiling_defs: Mapping[int, RegularTilingDefinition],
    *,
    buffered: bool = False,
    continent_order: list[str] | None = None,
) -> Equi7Grid:
    """Get user-defined Equi7Grid definition.

    Parameters
    ----------
    sampling: float | int | Dict[int | str, float | int]
            Grid sampling/pixel size specified as a single value or a dictionary with
            tiling IDs as keys and samplings as values.
    tiling_defs: Dict[int, RegularTilingDefinition]
            Tiling definition (stores name/tiling level and tile size).
    buffered: bool, optional
        If this flag is set to true, then the buffered projection zone
        will be used (defaults to false).
    continent_order: list[str] | None, optional
        Defines the usage and order of the continents.

    Returns
    -------
    Equi7Grid
        Equi7Grid instance.

    """
    proj_defs = get_system_definitions(buffered=buffered)
    return Equi7Grid.from_sampling(
        sampling, proj_defs, tiling_defs, system_order=continent_order
    )


def get_standard_equi7grid(
    sampling: SamplingFloatOrMap,
    *,
    buffered: bool = False,
    continent_order: list[str] | None = None,
) -> Equi7Grid:
    """Get standard Equi7Grid definition.

    Parameters
    ----------
    sampling: float | int | Dict[int | str, float | int]
            Grid sampling/pixel size specified as a single value or a dictionary with
            tiling IDs as keys and samplings as values.
    buffered: bool, optional
        If this flag is set to true, then the buffered projection zone
        will be used (defaults to false).
    continent_order: list[str] | None, optional
        Defines the usage and order of the continents.

    Returns
    -------
    Equi7Grid
        Equi7Grid instance.

    """
    return get_user_equi7grid(
        sampling,
        get_standard_tilings(),
        buffered=buffered,
        continent_order=continent_order,
    )


if __name__ == "__main__":
    pass
