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

"""Helper module to create Equi7Grid definition objects."""

from pathlib import Path
from typing import Generic, Literal

from pytileproj.grid import write_grid_def
from pytileproj.tiling_system import ProjSystemDefinition, RegularTilingDefinition

from equi7grid._types import T_co


class Equi7ProjSystemDefinition(ProjSystemDefinition[T_co], Generic[T_co]):
    """Defines a projection system for an Equi7 continent."""

    axis_orientation: tuple[Literal["E"], Literal["S"]] = ("E", "S")


def get_system_definitions() -> dict[str, Equi7ProjSystemDefinition]:
    """Get the 7 Equi7 projection system definitions."""
    zone_path = Path(__file__).parent / "data" / "zones"
    e7af = Equi7ProjSystemDefinition(
        name="AF",
        crs=27701,
        min_xy=(0, 0),
        proj_zone_geog=zone_path / "af_zone.parquet",
    )
    e7an = Equi7ProjSystemDefinition(
        name="AN",
        crs=27702,
        min_xy=(0, 0),
        proj_zone_geog=zone_path / "an_zone.parquet",
    )
    e7as = Equi7ProjSystemDefinition(
        name="AS",
        crs=27703,
        min_xy=(0, -1_800_000),
        proj_zone_geog=zone_path / "as_zone.parquet",
    )
    e7eu = Equi7ProjSystemDefinition(
        name="EU",
        crs=27704,
        min_xy=(0, -600_000),
        proj_zone_geog=zone_path / "eu_zone.parquet",
    )
    e7na = Equi7ProjSystemDefinition(
        name="NA",
        crs=27705,
        min_xy=(0, 0),
        proj_zone_geog=zone_path / "na_zone.parquet",
    )
    e7oc = Equi7ProjSystemDefinition(
        name="OC",
        crs=27706,
        min_xy=(0, 0),
        proj_zone_geog=zone_path / "oc_zone.parquet",
    )
    e7sa = Equi7ProjSystemDefinition(
        name="SA",
        crs=27707,
        min_xy=(0, 0),
        proj_zone_geog=zone_path / "sa_zone.parquet",
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
