# Copyright (c) 2026, TU Wien
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# The views and conclusions contained in the software and documentation are
# those of the authors and should not be interpreted as representing official
# policies, either expressed or implied, of the FreeBSD Project.

"""Core module defining the Equi7Grid classes."""

from collections.abc import Generator
from pathlib import Path

import shapely
from morecantile.models import Tile as RegularTile
from pytileproj._errors import TileOutOfZoneError
from pytileproj.grid import RegularGrid
from pytileproj.projgeom import (
    ProjGeom,
    convert_any_to_geog_geom,
    transform_geom_to_geog,
)
from pytileproj.tile import RasterTile
from pytileproj.tiling_system import (
    ProjSystemDefinition,
    RegularProjTilingSystem,
    RegularTilingDefinition,
    _tiling_access,
)

from equi7grid._const import KM, MAX_SAMPLING
from equi7grid.create_grids import get_system_definitions

Extent = tuple[int | float, int | float, int | float, int | float]
Equi7TileGenerator = Generator["Equi7Tile", "Equi7Tile", "Equi7Tile"]


class Equi7Tile(RasterTile):
    """Defines a tile in the Equi7Grid."""

    covers_land: bool


class Equi7TilingSystem(RegularProjTilingSystem):
    """Defines a tiling system for each Equi7Grid continent."""

    land_zone_geog: ProjGeom

    def _get_sampling_label(self, tiling_level: int) -> str:
        """Get the sampling label for a tiling.

        Parameters
        ----------
        tiling_level: int
            Tiling level or zoom.

        Returns
        -------
        str
            Sampling label (string representation of the sampling).

        """
        sampling = self[tiling_level].sampling

        if sampling < KM:
            sampling_str = f"{sampling}:03"
        elif (sampling >= KM) and (sampling < MAX_SAMPLING):
            sampling_str = f"{sampling}:.1f".replace(".", "K")
        else:
            err_msg = "Sampling labelling is not supported."
            raise ValueError(err_msg)

        return sampling_str

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
        e7_poly = ProjGeom(geom=poly, crs=self.pyproj_crs)
        geog_poly = transform_geom_to_geog(e7_poly)

        return shapely.intersects(geog_poly.geom, self.land_zone_geog.geom)

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

    def _tile_to_name(self, tile: RegularTile) -> str:
        """Create a tilename from a given regular tile object.

        Parameters
        ----------
        tile: RegularTile
            Regular tile object from TMS.

        Returns
        -------
        str
            Equi7 tilename.

        """
        e7_tile = self._tile_to_raster_tile(tile)
        ll_x, ll_y, _, _ = e7_tile.outer_boundary_extent
        tile_x, tile_y = int(ll_x / 1e5), int(ll_y / 1e5)
        y_label = "S" if tile_y < 0 else "N"
        if tile_y < 0:
            pass

        return f"E{tile_x:03}{y_label}{abs(tile_y):03}T{self.tilings[tile.z].name[1:]}"

    def _get_ftilename(self, tile: RegularTile) -> str:
        """Create full Equi7 tile name.

        Parameters
        ----------
        tile: RegularTile
            Regular tile object from TMS.

        Returns
        -------
        str
            Full Equi7 tile name.

        """
        tilename = self._tile_to_name(tile)
        return f"{self.name}{self._get_sampling_label(tile.z)}M_{tilename}"

    def _get_tile(self, tilename: str) -> RegularTile:
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
        tiling_level = int(tilename.split("T")[-1])
        y_label = tilename[4]
        y_mirror_fac = 1 if y_label == "N" else -1
        x, y = float(tilename[1:4]) * 1e5, float(tilename[5:8]) * 1e5 * y_mirror_fac

        return self._tms._tile(x, y, tiling_level)  # noqa: SLF001

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
        tile = self._get_tile(ftilename.split("_")[1])
        e7_tile = self._tile_to_raster_tile(tile, name=ftilename)
        if not self._tile_in_zone(e7_tile):
            raise TileOutOfZoneError(e7_tile)

        return e7_tile

    @_tiling_access
    def get_tiles_in_lonlat_bbox(
        self, bbox: tuple[float, float, float, float], tiling_id: int | str
    ) -> Equi7TileGenerator:
        """Get all Equi7 tiles intersecting with the geographic bounding box.

        Parameters
        ----------
        bbox: tuple[float, float, float, float]
            Bounding box (x_min, y_min, x_max, y_max) for selecting tiles.
        tiling_id: int | str
            Tiling level or name.
            Defaults to the first tiling level.

        Returns
        -------
        Equi7TileGenerator
            Yields Equi7 tile after tile, which intersects with the given
            bounding box.

        """
        for tile in self._get_tiles_in_geog_bbox(bbox, tiling_id):
            ftilename = self._get_ftilename(tile)
            e7_tile = self._tile_to_raster_tile(tile, name=ftilename)
            if not self._tile_in_zone(e7_tile):
                continue

            yield e7_tile


class Equi7Grid(RegularGrid):
    """Defines Equi7Grid with all sub-grid."""

    _rpts_cls = Equi7TilingSystem

    @staticmethod
    def _create_rpts_from_def(
        proj_def: ProjSystemDefinition,
        sampling: float | dict[int, float | int],
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
        return Equi7TilingSystem.from_sampling(
            sampling, proj_def, tiling_defs, land_zone_geog=land_zone_geog
        )


def get_standard_equi7grid(sampling: float | dict[int | str, float | int]) -> Equi7Grid:
    """Get standard Equi7Grid definition.

    Parameters
    ----------
    sampling: float | int | Dict[int | str, float | int]
            Grid sampling/pixel size specified as a single value or a dictionary with
            tiling IDs as keys and samplings as values.

    Returns
    -------
    Equi7Grid
        Equi7Grid instance.

    """
    json_path = Path(__file__).parent / "data" / "grids" / "equi7standard.json"
    return Equi7Grid.from_grid_def(json_path, sampling)


def get_user_equi7grid(
    sampling: float | dict[int | str, float | int],
    tiling_defs: dict[int, RegularTilingDefinition],
) -> Equi7Grid:
    """Get user-defined Equi7Grid definition.

    Parameters
    ----------
    sampling: float | int | Dict[int | str, float | int]
            Grid sampling/pixel size specified as a single value or a dictionary with
            tiling IDs as keys and samplings as values.
    tiling_defs: Dict[int, RegularTilingDefinition]
            Tiling definition (stores name/tiling level and tile size).

    Returns
    -------
    Equi7Grid
        Equi7Grid instance.

    """
    proj_defs = get_system_definitions()
    return Equi7Grid.from_sampling(sampling, proj_defs, tiling_defs)


if __name__ == "__main__":
    pass
