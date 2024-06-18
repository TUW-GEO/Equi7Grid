# Copyright (c) 2022, TU Wien, Department of Geodesy and Geoinformation (GEO).
# All rights reserved.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL VIENNA UNIVERSITY OF TECHNOLOGY, DEPARTMENT OF
# GEODESY AND GEOINFORMATION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
"""
Tests for the Equi7Grid.
"""

import unittest
import numpy as np
import numpy.testing as nptest
from pyproj import Transformer
from pyproj import CRS

from pytileproj.geometry import setup_test_geom_spitzbergen
from pytileproj.geometry import setup_geom_kamchatka
from pytileproj.geometry import setup_test_geom_siberia_antimeridian_180plus
from pytileproj.geometry import setup_test_geom_siberia_alaska

from equi7grid.equi7grid import Equi7Grid


### for testing at BBM machine
# # gdal 2
# import os
# os.environ["GDAL_DATA"] = r"C:\ProgramData\OSGeoW\share\gdal"
# os.environ["GDAL_DRIVER_PAT"] = r"C:\ProgramData\OSGeoW\bin\gdalplugins"

# # gdal 3
# import os
# os.environ["GDAL_DATA"] = r"C:\Program Files\GDAL\gdal-data"
# os.environ["GDAL_DRIVER_PAT"] = r"C:\Program Files\GDAL\gdalplugins"


class TestEqui7Grid(unittest.TestCase):

    def test_xy2lonlat_doubles(self):
        """
        Tests xy to lonlat projection using double numbers.
        """
        e7 = Equi7Grid(500)
        x = 5138743.127891
        y = 1307029.157093
        lon_should, lat_should = 15.1, 45.3
        lon, lat = e7.EU.xy2lonlat(x, y)
        nptest.assert_allclose(lon_should, lon)
        nptest.assert_allclose(lat_should, lat)

    def test_xy2lonlat_numpy_array(self):
        """
        Tests xy to lonlat projection using numpy arrays.
        """
        e7 = Equi7Grid(500)
        x = np.array([5138743.127891])
        y = np.array([1307029.157093])
        lon_should, lat_should = 15.1, 45.3
        lon, lat = e7.EU.xy2lonlat(x, y)
        nptest.assert_allclose(lon_should, lon)
        nptest.assert_allclose(lat_should, lat)

    def test_lonlat2xy_doubles(self):
        """
        Tests lonlat to xy projection using double numbers.
        """
        e7 = Equi7Grid(500)
        x_should = 5138743.127891
        y_should = 1307029.157093
        lon, lat = 15.1, 45.3
        sgrid_id, x, y = e7.lonlat2xy(lon, lat)
        assert sgrid_id == 'EU'
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)

    def test_lonlat2xy_antimeridian(self):
        """
        Tests lonlat to xy projection for locations close to the antimeridian.
        """
        # Far-east Siberia
        e7 = Equi7Grid(500)
        x_should = 7048122.707876
        y_should = 9238361.594967
        lon, lat = -178.5, 67.75
        sgrid_id, x, y = e7.lonlat2xy(lon, lat)
        assert sgrid_id == 'AS'
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)

        # Far-west Aleutian Islands
        x_should = 3887311.532849
        y_should = 7756934.345841
        lon, lat = -178.0, 51.75
        sgrid_id, x, y = e7.lonlat2xy(lon, lat)
        assert sgrid_id == 'NA'
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)

        # Far-east Aleutian Islands
        x_should = 3865149.386282
        y_should = 8432250.89933
        lon, lat = 173.0, 53.0
        sgrid_id, x, y = e7.lonlat2xy(lon, lat)
        assert sgrid_id == 'NA'
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)

    def test_lonlat2xy_numpy_array(self):
        """
        Tests lonlat to xy projection using numpy arrays.
        """
        e7 = Equi7Grid(500)
        x_should = np.array([5138743.127891, 5138743.127891])
        y_should = np.array([1307029.157093, 1307029.157093])
        lon = np.array([15.1, 15.1])
        lat = np.array([45.3, 45.3])
        sgrid_id, x, y = e7.lonlat2xy(lon, lat, 'EU')
        nptest.assert_array_equal(sgrid_id, np.array(['EU', 'EU']))
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)

    def test_lonlat2xy_numpy_array_no_sgrid(self):
        """
        Tests lonlat to xy projection using numpy arrays.
        """
        e7 = Equi7Grid(500)
        x_should = np.array([5138743.127891, 5138743.127891])
        y_should = np.array([1307029.157093, 1307029.157093])
        lon = np.array([15.1, 15.1])
        lat = np.array([45.3, 45.3])
        sgrid_id, x, y = e7.lonlat2xy(lon, lat)
        nptest.assert_array_equal(sgrid_id, np.array(['EU', 'EU']))
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)

    def test_calc_length_distortion_on_ellipsoid(self):
        """
        Tests calculation of length distortion with spherical solution
        """
        lon = -29.1
        lat = 24

        e7 = Equi7Grid(500)
        k = e7.calc_length_distortion_on_ellipsoid(lon, lat)

        k_should = 1.1432971727108836
        nptest.assert_almost_equal(k, k_should, decimal=6)

    def test_calc_length_distortion(self):
        """
        Tests calculation of length distortion with planar solution
        """
        x = 470916.85506916
        y = 8289708.44085851

        e7 = Equi7Grid(500)
        k = e7.AF.calc_length_distortion(x, y)

        k_should = 1.1432971727108836
        nptest.assert_almost_equal(k, k_should, decimal=6)

    def test_ij2xy(self):
        """
        Tests tile indices to xy coordination in the subgrid projection
        """
        # upperleft case (regular numpy array indexing)
        e7 = Equi7Grid(500)
        x_should = 3166500
        y_should = 5178000
        tile = e7.EU.tilesys.create_tile(x=3245631, y=5146545)
        x, y = tile.ij2xy(333, 444)
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)

        # pixel center coordinates
        x_should = 3166750
        y_should = 5177750
        tile = e7.EU.tilesys.create_tile(x=3245631, y=5146545)
        x, y = tile.ij2xy(333, 444, offset='center')
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)

        # lowerleft case
        x_should = 4800000
        y_should = 1200000
        tile = e7.EU.tilesys.create_tile(x=4800123, y=1200123)
        x, y = tile.ij2xy(0, 0, lowerleft=True, offset='ll')
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)

    def test_xy2ij(self):
        """
        Tests xy to tile array indices.
        """
        e7 = Equi7Grid(500)
        column_should = 333
        row_should = 444
        tile = e7.EU.tilesys.create_tile(x=3245631, y=5146545)
        column, row = tile.xy2ij(3166500, 5178000)
        nptest.assert_allclose(column_should, column)
        nptest.assert_allclose(row_should, row)

    def test_lonlat2ij_in_tile(self):
        """
        Tests the identification of column n rows indices in a Equi7Grid's tile

        """
        e7 = Equi7Grid(500)
        column_should = 1199
        row_should = 0
        tile_should = 'EU500M_E048N012T6'
        tilename, i, j = e7.lonlat2ij_in_tile(18.507, 44.571, lowerleft=True)
        nptest.assert_equal(i, column_should)
        nptest.assert_equal(j, row_should)
        nptest.assert_equal(tilename, tile_should)

    def test_lonlat2ij_in_tile(self):
        """
        Tests the tile name with option tile_names_in_m as True or False

        """
        e7 = Equi7Grid(3000, tile_names_in_m=True)
        column_should = 199
        row_should = 0
        tile_should = 'EU3000M_E048N012T6'
        tilename, i, j = e7.lonlat2ij_in_tile(18.507, 44.571, lowerleft=True)
        nptest.assert_equal(i, column_should)
        nptest.assert_equal(j, row_should)
        nptest.assert_equal(tilename, tile_should)

        e7 = Equi7Grid(3000, tile_names_in_m=False)
        column_should = 199
        row_should = 0
        tile_should = 'EU3K0M_E048N012T6'
        tilename, i, j = e7.lonlat2ij_in_tile(18.507, 44.571, lowerleft=True)
        nptest.assert_equal(i, column_should)
        nptest.assert_equal(j, row_should)
        nptest.assert_equal(tilename, tile_should)

    def test_proj4_reprojection_accuracy(self):
        """
        Tests the proj4 reproject accuracy by forward and backward reprojection.
        """
        # for subgrid Africa
        aeqd_wkt = ('PROJCS["Azimuthal_Equidistant",'
                    'GEOGCS["GCS_WGS_1984",'
                    'DATUM["D_WGS_1984",'
                    'SPHEROID["WGS_1984",6378137.0,298.257223563]],'
                    'PRIMEM["Greenwich",0.0],'
                    'UNIT["Degree",0.0174532925199433]],'
                    'PROJECTION["Azimuthal_Equidistant"],'
                    'PARAMETER["false_easting",5621452.01998],'
                    'PARAMETER["false_northing",5990638.42298],'
                    'PARAMETER["longitude_of_center",21.5],'
                    'PARAMETER["latitude_of_center",8.5],UNIT["Meter",1.0]]')
        aeqd_crs = CRS.from_wkt(aeqd_wkt)

        # test locations in Africa
        points = [(30.306273, -31.627336), (-43.880131, -14.589038),
                  (-35.261658, 79.423313), (10.457987, 23.456413)]

        for i, pt in enumerate(points):
            # from lat/lon to aeqd
            tf = Transformer.from_crs(('epsg', '4326'), aeqd_crs)
            aeqd_x, aeqd_y = tf.transform(pt[0], pt[1])

            # from aeqd to lat/lon
            tf = Transformer.from_crs(aeqd_crs, ('epsg', '4326'))
            lon, lat = tf.transform(aeqd_x, aeqd_y)

            nptest.assert_allclose(pt[0], lon)
            nptest.assert_allclose(pt[1], lat)

    def test_decode_tilename(self):
        """
        Tests the decoding of tilenames.
        """
        e7_500 = Equi7Grid(500)
        e7_10 = Equi7Grid(10)

        assert e7_500.EU.tilesys.decode_tilename('EU500M_E042N006T6') == \
               ('EU', 500, 600000, 4200000, 600000, 'T6')

        assert e7_10.OC.tilesys.decode_tilename('OC010M_E085N091T1') == \
               ('OC', 10, 100000, 8500000, 9100000, 'T1')

        assert e7_500.EU.tilesys.decode_tilename('E042N006T6') == \
               ('EU', 500, 600000, 4200000, 600000, 'T6')

        with nptest.assert_raises(ValueError) as excinfo:
            e7_10.EU.tilesys.decode_tilename('E042N006T6')
        assert str(excinfo.exception).startswith(
            '"tilename" is not properly defined!')

    def test_find_overlapping_tilenames(self):
        """
        Tests search for tiles which share the same extent_m but
        with different resolution and tilecode.
        """
        e7_500 = Equi7Grid(500)
        e7_10 = Equi7Grid(10)

        tiles1_should = ['EU025M_E042N006T3', 'EU025M_E042N009T3',
                         'EU025M_E045N006T3', 'EU025M_E045N009T3']
        tiles1 = e7_500.EU.tilesys.get_congruent_tiles_from_tilename(
            'EU500M_E042N006T6', target_sampling=25)
        assert sorted(tiles1) == sorted(tiles1_should)

        tiles2_should = [
            'E042N006T3', 'E042N009T3', 'E045N006T3', 'E045N009T3'
        ]
        tiles2 = e7_500.EU.tilesys.get_congruent_tiles_from_tilename(
            'E042N006T6', target_tiletype='T3')
        assert sorted(tiles2) == sorted(tiles2_should)

        tiles3_should = ['EU500M_E042N012T6']

        tiles3 = e7_10.EU.tilesys.get_congruent_tiles_from_tilename(
            'E044N015T1', target_sampling=500)
        assert sorted(tiles3) == sorted(tiles3_should)

        tiles4_should = ['E039N009T3']
        tiles4 = e7_10.EU.tilesys.get_congruent_tiles_from_tilename(
            'E041N011T1', target_tiletype='T3')
        assert sorted(tiles4) == sorted(tiles4_should)

    def test_search_tiles_lon_lat_extent_by_points(self):
        """
        Tests searching for tiles with input of lon lat points
        """
        e7 = Equi7Grid(500)
        tiles = e7.search_tiles_in_roi(points=[(10, 40), (5, 50),
                                               (-90.9, -1.2), (-175.2, 66)],
                                       coverland=True)

        desired_tiles = ['EU500M_E042N006T6', 'EU500M_E042N018T6',
                         'AS500M_E072N090T6', 'SA500M_E036N066T6']

        assert sorted(tiles) == sorted(desired_tiles)

    def test_search_tiles_lon_lat_extent(self):
        # begin-snippet: search-tiles-in-lon-lat-roi
        tiles = Equi7Grid(500).search_tiles_in_roi(bbox=[(0, 30), (10, 40)], coverland=True)
        assert sorted(tiles) == sorted([
            'EU500M_E036N006T6', 'EU500M_E042N000T6', 'EU500M_E042N006T6',
            'AF500M_E030N084T6', 'AF500M_E030N090T6', 'AF500M_E036N084T6',
            'AF500M_E036N090T6', 'AF500M_E042N084T6', 'AF500M_E042N090T6'])
        # end-snippet

    def test_find_all_tiles_in_global_lon_lat_extent(self):
        tiles_all = Equi7Grid(500).search_tiles_in_roi(bbox=[(-179.9, -89.9),
                                                             (179.9, 89.9)], coverland=True)
        assert (len(tiles_all)) == 854

    def test_search_tiles_lon_lat_extent_poles(self):
        e7 = Equi7Grid(500)

        tiles = e7.search_tiles_in_roi(bbox=[(-170, 88), (150.0, 90)])
        desired_tiles = ['NA500M_E078N084T6', 'NA500M_E078N090T6',
                         'NA500M_E084N084T6', 'NA500M_E084N090T6']
        assert sorted(tiles) == sorted(desired_tiles)

        tiles = e7.search_tiles_in_roi(bbox=[(-170, -90), (150.0, -89)])
        desired_tiles = ['AN500M_E036N030T6']
        assert tiles == desired_tiles

    def test_search_tiles_lon_lat_extent_antimeridian(self):
        """
        test antimeridian crossing bounding box
        &
        is only correct when segmentation in the function is done properly
            (as the curvature at this high latitude region is significant)

        """
        e7 = Equi7Grid(500)

        # test longitude values larger than 180 degrees
        tiles1 = e7.search_tiles_in_roi(bbox=[(179, 66), (210, 67.85)])

        desired_tiles = ['AS500M_E066N090T6', 'AS500M_E066N096T6', 'AS500M_E072N090T6',
                         'AS500M_E072N096T6', 'NA500M_E054N072T6', 'NA500M_E054N078T6',
                         'NA500M_E060N072T6']

        assert sorted(tiles1) == sorted(desired_tiles)

        # test longitudes that are more eastern, but lower in value
        tiles2 = e7.search_tiles_in_roi(bbox=[(179, 66), (-150, 67.85)])
        assert sorted(tiles2) == sorted(desired_tiles)

    def test_search_tiles_spitzbergen(self):
        """
        Tests the tile searching over Spitzbergen in the polar zone; ROI defined
        by a 4-corner polygon over high latitudes (is much curved on the globe).
        """

        grid = Equi7Grid(500)

        spitzbergen_geom = setup_test_geom_spitzbergen()
        spitzbergen_geom_tiles = sorted(['EU500M_E054N042T6', 'EU500M_E054N048T6',
                                         'EU500M_E060N042T6', 'EU500M_E060N048T6'])
        tiles = sorted(
            grid.search_tiles_in_roi(spitzbergen_geom, coverland=False))

        assert sorted(tiles) == sorted(spitzbergen_geom_tiles)

    def test_search_tiles_siberia_antimeridian(self):
        """
        Tests the tile searching over Siberia and Alaska in the polar zone; ROI defined
        by a 4-corner polygon over high latitudes (is much curved on the globe).

        Comprises:
            - do not return tiles when not intersecting the zone
            - interpret correctly longitudes higher than 180 degrees
        """

        grid = Equi7Grid(500)

        geom_siberia_tiles = sorted(['AS500M_E066N090T6', 'AS500M_E072N090T6'])
        poly_siberia_antim_180plus = setup_test_geom_siberia_antimeridian_180plus()
        tiles = sorted(grid.search_tiles_in_roi(poly_siberia_antim_180plus, coverland=False))

        assert sorted(tiles) == sorted(geom_siberia_tiles)

        geom_siberia_alaska_tiles = sorted([
            'AS500M_E066N090T6', 'AS500M_E072N090T6', 'AS500M_E072N096T6',
            'NA500M_E054N072T6', 'NA500M_E054N078T6', 'NA500M_E060N078T6'])

        poly_siberia_alaska = setup_test_geom_siberia_alaska()
        tiles = sorted(
            grid.search_tiles_in_roi(poly_siberia_alaska, coverland=True))

        assert sorted(tiles) == sorted(geom_siberia_alaska_tiles)

    def test_search_tiles_kamchatka(self):
        """
        Tests the tile searching over Kamchatka in far east Sibiria;

        This test is especially nice, as it contains also a tile that covers both,
        the ROI and the continental zone, but the intersection of the tile and
        the ROI is outside of the zone.

        Furthermore, it also covers Equi7Grid subgrids that consist of a multipolygon,
        as they overspan the 180deg/dateline.
        """

        grid = Equi7Grid(500)

        kamchatka_geom = setup_geom_kamchatka()
        kamchatka_geom_tiles = sorted([
            'AS500M_E072N078T6', 'AS500M_E078N078T6', 'AS500M_E078N084T6',
            'NA500M_E036N078T6', 'NA500M_E036N084T6', 'NA500M_E042N078T6',
            'NA500M_E042N084T6'
        ])

        tiles = grid.search_tiles_in_roi(kamchatka_geom, coverland=False)

        assert sorted(tiles) == sorted(kamchatka_geom_tiles)

    def test_identify_tiles_overlapping_xybbox(self):
        """
        Tests identification of tiles covering a bounding box
        given in equi7 coordinats
        """

        e7_500 = Equi7Grid(500)
        e7_10 = Equi7Grid(10)

        tiles1_should = [
            'EU500M_E048N006T6', 'EU500M_E054N006T6', 'EU500M_E060N006T6',
            'EU500M_E048N012T6', 'EU500M_E054N012T6', 'EU500M_E060N012T6'
        ]

        tiles2_should = [
            'EU010M_E051N011T1', 'EU010M_E052N011T1', 'EU010M_E051N012T1',
            'EU010M_E052N012T1'
        ]

        tiles1 = e7_500.EU.tilesys.identify_tiles_overlapping_xybbox(
            [5138743, 1111111, 6200015, 1534657])

        tiles2 = e7_10.EU.tilesys.identify_tiles_overlapping_xybbox(
            [5138743, 1111111, 5299999, 1234657])

        assert sorted(tiles1) == sorted(tiles1_should)
        assert sorted(tiles2) == sorted(tiles2_should)

    def test_create_tiles_overlapping_xybbox(self):
        e7_500 = Equi7Grid(500)

        tiles = e7_500.EU.tilesys.create_tiles_overlapping_xybbox(
            [5138743, 1111111, 6200015, 1777777])

        nptest.assert_equal(tiles.shape, (2, 3))

        nptest.assert_equal([x.name for x in tiles.flatten()], [
            'EU500M_E048N012T6', 'EU500M_E054N012T6', 'EU500M_E060N012T6',
            'EU500M_E048N006T6', 'EU500M_E054N006T6', 'EU500M_E060N006T6'
        ])

        nptest.assert_equal(tiles[0, 0].active_subset_px, (677, 0, 1200, 1155))
        nptest.assert_equal(tiles[1, 0].active_subset_px,
                            (677, 1022, 1200, 1200))
        nptest.assert_equal(tiles[0, 2].active_subset_px, (0, 0, 400, 1155))

    def test_get_covering_tiles(self):
        """
        Tests the search for co-locating tiles of other type.
        """

        e7g_40 = Equi7Grid(40)
        e7g_10 = Equi7Grid(10)

        fine_tiles = [
            'EU010M_E005N058T1', 'EU010M_E005N059T1', 'EU010M_E005N060T1',
            'EU010M_E005N061T1'
        ]

        target_tiletype = e7g_40.get_tiletype()
        target_sampling = e7g_40.core.sampling

        # invoke the results as tile name in short form
        coarse_tiles_shortform = e7g_10.EU.tilesys.collect_congruent_tiles(
            fine_tiles, target_tiletype=target_tiletype)

        # invoke the results as tile name in long form
        coarse_tiles_longform = e7g_10.EU.tilesys.collect_congruent_tiles(
            fine_tiles, target_sampling=target_sampling)

        assert sorted(coarse_tiles_shortform) == ['E003N057T3', 'E003N060T3']
        assert sorted(coarse_tiles_longform) == [
            'EU040M_E003N057T3', 'EU040M_E003N060T3'
        ]

    def test_tile_get_extent_geometry_geog(self):
        """
        Tests the geometry functions of the tile object.

        """
        e7g = Equi7Grid(500)

        tile_up_north = e7g.create_tile('EU500M_E054N054T6')

        # more precise would be: (-35.42781, 81.57133, 55.67043, 87.77046) with decimal=5
        nptest.assert_almost_equal(tile_up_north.bbox_geog,
                                   (-35.43, 81.57, 55.67, 87.77),
                                   decimal=2)

        tile_antimeridian = e7g.create_tile('NA500M_E042N078T6')

        nptest.assert_almost_equal(tile_antimeridian.bbox_geog,
                                   (174.93808, 54.32175, -172.33224, 60.55437),
                                   decimal=5)


if __name__ == '__main__':
    unittest.main()
