# Copyright (c) 2016,Vienna University of Technology,
# Department of Geodesy and Geoinformation
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

'''
Tests for the equi7 grid class.
'''
import os

from equi7grid.equi7grid import Equi7Grid

import numpy as np
import numpy.testing as nptest

def test_ij2xy():
    """
    Test xy to lonlat projection using double numbers.
    """
    e7 = Equi7Grid(500)
    x_should = 3166500
    y_should = 5178000
    tile = e7.EU.tilesys.create_tile(x=3245631, y=5146545)
    x, y = tile.ij2xy(333, 444)
    nptest.assert_allclose(x_should, x)
    nptest.assert_allclose(y_should, y)

def test_xy2ij():
    """
    Test xy to lonlat projection using double numbers.
    """
    e7 = Equi7Grid(500)
    column_should = 333
    row_should = 444
    tile = e7.EU.tilesys.create_tile(x=3245631, y=5146545)
    column, row = tile.xy2ij(3166500, 5178000)
    nptest.assert_allclose(column_should, column)
    nptest.assert_allclose(row_should, row)


def test_equi7xy2lonlat_doubles():
    """
    Test xy to lonlat projection using double numbers.
    """
    e7 = Equi7Grid(500)
    x = 5138743.127891
    y = 1307029.157093
    lon_should, lat_should = 15.1, 45.3
    lon, lat = e7.EU.xy2lonlat(x,y)
    nptest.assert_allclose(lon_should, lon)
    nptest.assert_allclose(lat_should, lat)


def test_equi7xy2lonlat_numpy_array():
    """
    Test xy to lonlat projection using numpy arrays.
    """
    e7 = Equi7Grid(500)
    x = np.array([5138743.127891])
    y = np.array([1307029.157093])
    lon_should, lat_should = 15.1, 45.3
    lon, lat = e7.EU.xy2lonlat(x,y)
    nptest.assert_allclose(lon_should, lon)
    nptest.assert_allclose(lat_should, lat)


def test_lonlat2equi7xy_doubles():
    """
    Test lonlat to xy projection using double numbers.
    """
    e7 = Equi7Grid(500)
    x_should = 5138743.127891
    y_should = 1307029.157093
    lon, lat = 15.1, 45.3
    sgrid_id, x, y = e7.lonlat2xy(lon, lat)
    assert sgrid_id == 'EU'
    nptest.assert_allclose(x_should, x)
    nptest.assert_allclose(y_should, y)


def test_lonlat2equi7xy_numpy_array():
    """
    Test lonlat to xy projection using numpy arrays.
    """
    e7 = Equi7Grid(500)
    x_should = np.array([5138743.127891,
                         5138743.127891])
    y_should = np.array([1307029.157093,
                         1307029.157093])
    lon = np.array([15.1, 15.1])
    lat = np.array([45.3, 45.3])
    sgrid_id, x, y = e7.lonlat2xy(lon, lat, 'EU')
    nptest.assert_array_equal(sgrid_id, np.array(['EU', 'EU']))
    nptest.assert_allclose(x_should, x)
    nptest.assert_allclose(y_should, y)


def test_lonlat2equi7xy_numpy_array_no_sgrid():
    """
    Test lonlat to xy projection using numpy arrays.
    """
    e7 = Equi7Grid(500)
    x_should = np.array([5138743.127891,
                         5138743.127891])
    y_should = np.array([1307029.157093,
                         1307029.157093])
    lon = np.array([15.1, 15.1])
    lat = np.array([45.3, 45.3])
    sgrid_id, x, y = e7.lonlat2xy(lon, lat)
    nptest.assert_array_equal(sgrid_id, np.array(['EU', 'EU']))
    nptest.assert_allclose(x_should, x)
    nptest.assert_allclose(y_should, y)


def test_proj4_reprojection_accuracy():

    # test the proj4 reproject accuracy by forward and backward reprojection
    from osgeo import osr
    import pyproj

    geo_sr = osr.SpatialReference()
    geo_sr.SetWellKnownGeogCS("EPSG:4326")
    # Africa
    aeqd_proj = ('PROJCS["Azimuthal_Equidistant",'
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
    aeqd_sr = osr.SpatialReference()
    aeqd_sr.ImportFromWkt(aeqd_proj)
    p_grid = pyproj.Proj(aeqd_sr.ExportToProj4())
    p_geo = pyproj.Proj(geo_sr.ExportToProj4())

    # test locations in Africa
    points = [(-31.627336, 30.306273),
              (-14.589038, -43.880131),
              (79.423313, -35.261658),
              (23.456413, 10.457987)]

    for i, pt in enumerate(points):
        # from lat/lon to aeqd
        aeqd_x, aeqd_y = pyproj.transform(p_geo, p_grid, pt[0], pt[1])
        # from aeqd to lat/lon
        lon, lat = pyproj.transform(p_grid, p_geo, aeqd_x, aeqd_y)
        # print info

        print("testing location {}:".format(i))
        _info = "   ({:f},{:f}) -> ({:f},{:f}) -> ({:f},{:f})"
        print(info.format(pt[0], pt[1], aeqd_x, aeqd_y, lon, lat))
        print("    difference: ({:f},{:f})".format(lon - pt[0], lat - pt[1]))
        nptest.assert_allclose(pt[0], lon)
        nptest.assert_allclose(pt[1], lat)


def test_search_tile_500_lon_lat_extent():
    """
    Test searching of tile with input of lon lat extent
    """
    e7 = Equi7Grid(500)
    tiles = e7.search_tiles_in_geo_roi(extent=[10, 40, 20, 50],
                                       coverland=True)

    assert tiles == ['EU500M_E042N006T6',
                     'EU500M_E042N012T6',
                     'EU500M_E048N006T6',
                     'EU500M_E048N012T6',
                     'EU500M_E048N018T6',
                     'EU500M_E054N006T6',
                     'EU500M_E054N012T6',
                     'EU500M_E054N018T6',
                     'AF500M_E042N090T6']