# Copyright (c) 2026, TU Wien
# Licensed under the MIT License. See LICENSE file.

"""Helper module to create Equi7Grid definition objects."""

from pathlib import Path
from typing import Literal

from pytileproj import ProjSystemDefinition, RegularTilingDefinition
from pytileproj.grid import write_grid_def

from equi7grid._types import T_co


class Equi7ProjSystemDefinition(ProjSystemDefinition[T_co]):
    """Defines a projection system for an Equi7 continent."""

    axis_orientation: tuple[Literal["E"], Literal["S"]] = ("E", "S")


def get_system_definitions(
    *,
    buffered: bool = False,
) -> dict[str, Equi7ProjSystemDefinition]:
    """Get the 7 Equi7 projection system definitions.

    Parameters
    ----------
    buffered: bool, optional
        If this flag is set to true, then the buffered
        projection zone will be used (defaults to false).

    Returns
    -------
    dict[str, Equi7ProjSystemDefinition]

    """
    suffix = "_bfrd" if buffered else ""

    zone_path = Path(__file__).parent / "data" / "zones"
    e7af = Equi7ProjSystemDefinition(
        name="AF",
        crs=27701,
        min_xy=(0, 0),
        proj_zone_geog=zone_path / f"af_zone{suffix}.parquet",
    )
    e7an = Equi7ProjSystemDefinition(
        name="AN",
        crs=27702,
        min_xy=(0, 0),
        proj_zone_geog=zone_path / f"an_zone{suffix}.parquet",
    )
    e7as = Equi7ProjSystemDefinition(
        name="AS",
        crs=27703,
        min_xy=(0, -1_800_000),
        proj_zone_geog=zone_path / f"as_zone{suffix}.parquet",
    )
    e7eu = Equi7ProjSystemDefinition(
        name="EU",
        crs=27704,
        min_xy=(0, -600_000),
        proj_zone_geog=zone_path / f"eu_zone{suffix}.parquet",
    )
    e7na = Equi7ProjSystemDefinition(
        name="NA",
        crs=27705,
        min_xy=(0, 0),
        proj_zone_geog=zone_path / f"na_zone{suffix}.parquet",
    )
    e7oc = Equi7ProjSystemDefinition(
        name="OC",
        crs=27706,
        min_xy=(0, 0),
        proj_zone_geog=zone_path / f"oc_zone{suffix}.parquet",
    )
    e7sa = Equi7ProjSystemDefinition(
        name="SA",
        crs=27707,
        min_xy=(0, 0),
        proj_zone_geog=zone_path / f"sa_zone{suffix}.parquet",
    )
    return {
        e7af.name: e7af,
        e7an.name: e7an,
        e7as.name: e7as,
        e7eu.name: e7eu,
        e7oc.name: e7oc,
        e7na.name: e7na,
        e7sa.name: e7sa,
    }


def get_standard_tilings() -> dict[int, RegularTilingDefinition]:
    """Get standard Equi7Grid tiling defintions."""
    return {
        1: RegularTilingDefinition(name="T6", tile_shape=600_000),
        2: RegularTilingDefinition(name="T3", tile_shape=300_000),
        3: RegularTilingDefinition(name="T1", tile_shape=100_000),
    }


def create_standard_equi7() -> None:
    """Create standard Equi7Grid defintion on disk."""
    tiling_defs = get_standard_tilings()
    proj_defs = get_system_definitions()
    json_path = Path(__file__).parent / "data" / "grids" / "equi7standard.json"
    write_grid_def(json_path, proj_defs, tiling_defs)


if __name__ == "__main__":
    create_standard_equi7()
