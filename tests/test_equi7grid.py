# Copyright (c) 2026, TU Wien
# Licensed under the MIT License. See LICENSE file.

"""
Tests for the Equi7Grid.
"""

import numpy.testing as nptest
import pytest
import shapely
from pyproj import CRS, Transformer
from pytileproj import GeogGeom, GeomOutOfZoneError, TileOutOfZoneError

from equi7grid import get_standard_equi7grid
from equi7grid._core import Equi7Grid, Equi7TileGenerator


@pytest.fixture(scope="module")
def e7grid():
    return get_standard_equi7grid(500)


@pytest.fixture(scope="module")
def poly_siberia_alaska() -> GeogGeom:
    points = [
        (177.6545884597184, 67.05574774066811),
        (179.0195867605756, 65.33232820668778),
        (198.4723636216472 - 360, 66.06909015550372),
        (198.7828129097253 - 360, 68.14247939909886),
    ]
    return GeogGeom(geom=shapely.Polygon(points))


@pytest.fixture(scope="module")
def poly_spitzbergen() -> GeogGeom:
    points = [
        (8.391827331539572, 77.35762113396143),
        (25.43098663332705, 75.61353436967198),
        (40.50119498304080, 79.73786853853339),
        (16.87007957357446, 81.59290885863483),
    ]
    return GeogGeom(geom=shapely.Polygon(points))


def assert_tiles(tiles: Equi7TileGenerator, ref_tilenames: list[str]):
    tilenames = [str(tile.name) for tile in tiles]
    assert sorted(tilenames) == sorted(ref_tilenames)


def test_xylonlat_doubles(e7grid: Equi7Grid):
    x = 5138743.127891
    y = 1307029.157093
    lon_should, lat_should = 15.1, 45.3
    geog_coord = e7grid.EU.xy_to_lonlat(x, y)
    nptest.assert_allclose(lon_should, geog_coord.x)
    nptest.assert_allclose(lat_should, geog_coord.y)


def test_lonlatxy_doubles_get_rpts(e7grid: Equi7Grid):
    x_should = 5138743.127891
    y_should = 1307029.157093
    lon, lat = 15.1, 45.3
    rpts = e7grid.get_systems_from_lonlat(lon, lat)[0]
    assert rpts.name == "EU"
    proj_coord = rpts.lonlat_to_xy(lon, lat)
    nptest.assert_allclose(x_should, proj_coord.x)
    nptest.assert_allclose(y_should, proj_coord.y)


def test_lonlatxy_antimeridian(e7grid: Equi7Grid):
    # Far-east Siberia
    x_should = 7048122.707876
    y_should = 9238361.594967
    lon, lat = -178.5, 67.75
    rpts = e7grid.get_systems_from_lonlat(lon, lat)[0]
    assert rpts.name == "AS"
    e7_coord = rpts.lonlat_to_xy(lon, lat)
    nptest.assert_allclose(x_should, e7_coord.x)
    nptest.assert_allclose(y_should, e7_coord.y)

    # Far-west Aleutian Islands
    x_should = 3887311.532849
    y_should = 7756934.345841
    lon, lat = -178.0, 51.75
    rpts = e7grid.get_systems_from_lonlat(lon, lat)[0]
    assert rpts.name == "NA"
    e7_coord = rpts.lonlat_to_xy(lon, lat)
    nptest.assert_allclose(x_should, e7_coord.x)
    nptest.assert_allclose(y_should, e7_coord.y)

    # Far-east Aleutian Islands
    x_should = 3865149.386282
    y_should = 8432250.89933
    lon, lat = 173.0, 53.0
    rpts = e7grid.get_systems_from_lonlat(lon, lat)[0]
    assert rpts.name == "NA"
    e7_coord = rpts.lonlat_to_xy(lon, lat)
    nptest.assert_allclose(x_should, e7_coord.x)
    nptest.assert_allclose(y_should, e7_coord.y)


def test_calc_length_distortion_on_ellipsoid(e7grid: Equi7Grid):
    lon = -29.1
    lat = 24
    k = e7grid.calc_length_distortion_on_ellipsoid(lon, lat)

    k_should = 1.1432971727108836
    nptest.assert_almost_equal(k, k_should, decimal=6)


def test_calc_length_distortion(e7grid: Equi7Grid):
    x = 470916.85506916
    y = 8289708.44085851
    k = e7grid.AF.calc_length_distortion(x, y)

    k_should = 1.1432971727108836
    nptest.assert_almost_equal(k, k_should, decimal=6)


def test_rc2xy(e7grid: Equi7Grid):
    x_should = 3166500
    y_should = 5178000
    tile = e7grid.EU.get_tile_from_xy(3245631, 5146545, tiling_id="T6")
    x, y = tile.rc2xy(444, 333)
    nptest.assert_allclose(x_should, x)
    nptest.assert_allclose(y_should, y)

    # pixel center coordinates
    x_should = 3166750
    y_should = 5177750
    tile = e7grid.EU.get_tile_from_xy(x=3245631, y=5146545, tiling_id="T6")
    x, y = tile.rc2xy(444, 333, px_origin="c")
    nptest.assert_allclose(x_should, x)
    nptest.assert_allclose(y_should, y)

    # lowerleft case
    x_should = 4800000
    y_should = 1200000
    tile = e7grid.EU.get_tile_from_xy(x=4800123, y=1200123, tiling_id="T6")
    x, y = tile.rc2xy(tile.n_rows - 1, 0, px_origin="ll")
    nptest.assert_allclose(x_should, x)
    nptest.assert_allclose(y_should, y)


def test_xy2rc(e7grid: Equi7Grid):
    c_should = 333
    r_should = 444
    tile = e7grid.EU.get_tile_from_xy(x=3245631, y=5146545)
    r, c = tile.xy2rc(3166500, 5178000)
    nptest.assert_allclose(c_should, c)
    nptest.assert_allclose(r_should, r)


def test_lonlat2rc_in_tile(e7grid: Equi7Grid):
    lon, lat = 18.507, 44.571
    tile = e7grid.EU.get_tile_from_lonlat(lon, lat, tiling_id="T6")
    e7_coord = e7grid.EU.lonlat_to_xy(lon, lat)
    r, c = tile.xy2rc(e7_coord.x, e7_coord.y)
    tile_should = "EU_E048N012T6"
    c_should = 1199
    r_should = tile.n_rows - 1
    nptest.assert_equal(r, r_should)
    nptest.assert_equal(c, c_should)
    nptest.assert_equal(tile.name, tile_should)


def test_proj4_reprojection_accuracy():
    aeqd_wkt = (
        'PROJCS["Azimuthal_Equidistant",'
        'GEOGCS["GCS_WGS_1984",'
        'DATUM["D_WGS_1984",'
        'SPHEROID["WGS_1984",6378137.0,298.257223563]],'
        'PRIMEM["Greenwich",0.0],'
        'UNIT["Degree",0.0174532925199433]],'
        'PROJECTION["Azimuthal_Equidistant"],'
        'PARAMETER["false_easting",5621452.01998],'
        'PARAMETER["false_northing",5990638.42298],'
        'PARAMETER["longitude_of_center",21.5],'
        'PARAMETER["latitude_of_center",8.5],UNIT["Meter",1.0]]'
    )
    aeqd_crs = CRS.from_wkt(aeqd_wkt)

    # test locations in Africa
    points = [
        (30.306273, -31.627336),
        (-43.880131, -14.589038),
        (-35.261658, 79.423313),
        (10.457987, 23.456413),
    ]

    for pt in points:
        # from lat/lon to aeqd
        tf = Transformer.from_crs(("epsg", "4326"), aeqd_crs)
        aeqd_x, aeqd_y = tf.transform(pt[0], pt[1])

        # from aeqd to lat/lon
        tf = Transformer.from_crs(aeqd_crs, ("epsg", "4326"))
        lon, lat = tf.transform(aeqd_x, aeqd_y)

        nptest.assert_allclose(pt[0], lon)
        nptest.assert_allclose(pt[1], lat)


def test_decode_tilename(e7grid: Equi7Grid):
    tile = e7grid.get_tile_from_name("EU_E042N006T6")
    sampling = 500
    assert tile.x_pixel_size == sampling
    assert tile.outer_boundary_corners[0] == (4200000, 600000)

    try:
        tile = e7grid.get_tile_from_name("EU_E242N006T6")
        raise AssertionError
    except TileOutOfZoneError:
        assert True


def test_find_overlapping_tilenames(e7grid: Equi7Grid):
    tiles_should = [
        "EU_E042N006T3",
        "EU_E042N009T3",
        "EU_E045N006T3",
        "EU_E045N009T3",
    ]
    tiles = e7grid.EU.get_children_from_name("EU_E042N006T6")
    assert_tiles(tiles, tiles_should)

    tile_should = "EU_E039N009T3"
    tile = e7grid.EU.get_parent_from_name("EU_E041N011T1")
    assert tile.name == tile_should

    tile_should = "EU_E042N012T6"
    tile = e7grid.EU.get_parent_from_name("EU_E044N015T1")
    tile = e7grid.EU.get_parent_from_name(str(tile.name))
    assert tile.name == tile_should


def test_search_tiles_geog_bbox(e7grid: Equi7Grid):
    tiles = e7grid.get_tiles_in_geog_bbox(
        bbox=(0, 30, 10, 40), tiling_id="T6", cover_land=True
    )
    tiles_should = [
        "EU_E036N006T6",
        "EU_E042N000T6",
        "EU_E042N006T6",
        "AF_E030N084T6",
        "AF_E030N090T6",
        "AF_E036N084T6",
        "AF_E036N090T6",
        "AF_E042N084T6",
        "AF_E042N090T6",
    ]
    assert_tiles(tiles, tiles_should)


# noqa: TODO: bbm (validate if tiles are as desired)
def test_find_all_tiles_with_global_bbox(e7grid: Equi7Grid):
    tiles_all = e7grid.get_tiles_in_geog_bbox(
        bbox=(-179.9, -89.9, 179.9, 89.9), tiling_id="T6", cover_land=True
    )
    n_t6_tiles_global = 864
    assert len(list(tiles_all)) == n_t6_tiles_global  # 854 was the old number


def test_search_tiles_geog_extent_poles(e7grid: Equi7Grid):
    tiles = e7grid.get_tiles_in_geog_bbox(bbox=(-170, 88, 150.0, 90), tiling_id="T6")
    tiles_should = [
        "NA_E078N084T6",
        "NA_E078N090T6",
        "NA_E084N084T6",
        "NA_E084N090T6",
    ]
    assert_tiles(tiles, tiles_should)

    tiles = e7grid.get_tiles_in_geog_bbox(bbox=(-170, -90, 150.0, -89), tiling_id="T6")
    tiles_should = ["AN_E036N030T6"]
    assert_tiles(tiles, tiles_should)


def test_search_tiles_geog_extent_antimeridian(e7grid: Equi7Grid):
    tiles = e7grid.get_tiles_in_geog_bbox(bbox=(179, 66, -150, 67.85), tiling_id="T6")

    tiles_should = [
        "AS_E066N090T6",
        "AS_E066N096T6",
        "AS_E072N090T6",
        "AS_E072N096T6",
        "NA_E054N072T6",
        "NA_E054N078T6",
        "NA_E060N072T6",
    ]

    assert_tiles(tiles, tiles_should)


def test_search_tiles_spitzbergen(e7grid: Equi7Grid, poly_spitzbergen: GeogGeom):
    tiles_should = [
        "EU_E054N042T6",
        "EU_E054N048T6",
        "EU_E060N042T6",
        "EU_E060N048T6",
    ]
    tiles = e7grid.get_tiles_in_geom(poly_spitzbergen, tiling_id="T6")

    assert_tiles(tiles, tiles_should)


def test_search_tiles_siberia_antimeridian(
    e7grid: Equi7Grid, poly_siberia_alaska: GeogGeom
):
    tiles_should = [
        "AS_E066N090T6",
        "AS_E066N096T6",
        "AS_E072N090T6",
        "AS_E072N096T6",
        "NA_E054N072T6",
        "NA_E054N078T6",
        "NA_E060N078T6",
    ]
    tiles = e7grid.get_tiles_in_geom(
        poly_siberia_alaska, tiling_id="T6", cover_land=False
    )

    assert_tiles(tiles, tiles_should)


def test_search_tiles_siberia_antimeridian_land(
    e7grid: Equi7Grid, poly_siberia_alaska: GeogGeom
):
    tiles_should = [
        "AS_E066N090T6",
        "AS_E072N090T6",
        "NA_E054N072T6",
        "NA_E054N078T6",
        "NA_E060N078T6",
    ]
    tiles = e7grid.get_tiles_in_geom(
        poly_siberia_alaska, tiling_id="T6", cover_land=True
    )

    assert_tiles(tiles, tiles_should)


def test_identify_tiles_overlapping_xybbox(e7grid: Equi7Grid):
    tiles_should = [
        "EU_E048N006T6",
        "EU_E054N006T6",
        "EU_E060N006T6",
        "EU_E048N012T6",
        "EU_E054N012T6",
        "EU_E060N012T6",
    ]

    tiles = e7grid.EU.get_tiles_in_bbox(
        [5138743, 1111111, 6200015, 1534657], tiling_id="T6"
    )

    assert_tiles(tiles, tiles_should)


def test_lonlat_to_xy_continental_transition(e7grid: Equi7Grid):
    lon, lat = 51.12, 61.71
    proj_coords = e7grid.lonlat_to_xy(lon, lat)
    assert len(proj_coords) == 1
    assert next(iter(proj_coords.keys())) == "AS"


def test_lonlat_to_xy_continental_transition_bfrd():
    lon, lat = 51.12, 61.71
    e7grid_bfrd = get_standard_equi7grid(500, buffered=True)
    proj_coords = e7grid_bfrd.lonlat_to_xy(lon, lat)
    n_continents = 2
    assert len(proj_coords) == n_continents
    assert list(proj_coords.keys()) == ["AS", "EU"]


def test_lonlat_to_xy_system_order():
    lon, lat = 51.12, 61.71
    e7grid_eu = get_standard_equi7grid(500, continent_order=["EU"])
    try:
        _ = e7grid_eu.lonlat_to_xy(lon, lat)
        raise AssertionError
    except GeomOutOfZoneError:
        assert True


def test_lonlat_to_xy_system_order_bfrd():
    lon, lat = 51.12, 61.71
    e7grid_bfrd = get_standard_equi7grid(
        500, buffered=True, continent_order=["EU", "AS"]
    )
    proj_coords = e7grid_bfrd.lonlat_to_xy(lon, lat)
    n_continents = 2
    assert len(proj_coords) == n_continents
    assert list(proj_coords.keys()) == ["EU", "AS"]


if __name__ == "__main__":
    pass
