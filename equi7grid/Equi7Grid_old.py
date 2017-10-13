# Copyright (c) 2017, Vienna University of Technology (TU Wien), Department of
# Geodesy and Geoinformation (GEO).
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


'''
Created on November 05, 2015

Equi7Grid

@author: Bernhard Bauer-Marschallinger, bbm@geo.tuwien.ac.at

Note:

    Geos support of Gdal library should be enabled for the accurate spatial
    operation, otherwise, the search overlapped tiles will not be completely
    accurate. Equi7Grid also depends the grids' zone shape_px files of both aeqd
    and wgs84.

    terminology:
    grid =  grid name (Equi7)
    res= grid resolution/grid spacing (5,10,40,75,500 meters)
    subgrid= subgrid full name (e.g. "AF500M")
    subgrid_id= sub-grid ini (e.g. "AF")
    tile= tile name ("E036N090T6")
    tile_code= tile code representing width of tile (T1:100km, T3:300km,
                T6:600km)
    ftile= full name of tile/subgrid+tile ("AF500M_E036N090T6")


'''


import os
import itertools
import pickle
import click
from command_line import cli

import pyproj
import numpy as np
from osgeo import ogr
from osgeo import osr

import platform
import gdalport

def _load_equi7grid_data(module_path):
    # load the data, raise the error if failed to load equi7grid.dat
    equi7_data = None
    system = platform.system()
    fname = os.path.join(os.path.dirname(module_path), "data", "equi7grid_{0}.dat".format(system))
    with open(fname, "rb") as f:
        equi7_data = pickle.load(f)
    return equi7_data


class Equi7Grid(object):
    """
    Equi7 Grid

    Parameters
    ----------
    res : float
        The tile resolution
    """

    # static attribute
    _static_equi7_data = _load_equi7grid_data(__file__)
    # sub grid IDs
    _static_subgrid_ids = ["NA", "EU", "AS", "SA", "AF", "OC", "AN"]
    # supported tile widths(resolution)
    _static_tilecodes = ["T6", "T3", "T1"]
    # supported grid spacing (resolution)
    _static_res = [1000, 800, 750, 600, 500, 400, 300, 250, 200,
                   150, 125, 100, 96, 80, 75, 64, 60, 50, 48, 40,
                   32, 30, 25, 24, 20, 16, 10, 8, 5, 4, 2, 1]

    def __init__(self, res):
        """
        construct Equi7 grid system.

        """
        self.res = int(res)
        # initializing
        self._initialize()

    def _initialize(self):
        """
        initialization
        """

        # check if the equi7grid.data have been loaded successfully
        if Equi7Grid._static_equi7_data is None:
            self.res = None
            raise ValueError("cannot load Equi7Grid ancillary data!")

        self.tilecode, self.tile_size_m = self.link_res_2_tile(self.res, get_span=True)

        # keep a reference to the _static_data
        self._equi7_data = Equi7Grid._static_equi7_data


    @property
    def subgrids(self):
        return Equi7Grid._static_subgrid_ids

    def is_coverland(self, ftile):
        """check if tile covers land

        Parameters
        ----------
        ftile : string
        full tile name
        """
        return Equi7Tile(ftile).covers_land()

    def get_subgrid_zone_geom(self, subgrid_id):
        """return sub-grid extent_m geometry

        Parameters
        ----------
        subgrid_id : string
            sub-grid id string e.g. EU for Europe

        Return
        ------
        OGRGeomtery
            a geometry representing the extent_m of given sub-grid

        """
        geom = ogr.CreateGeometryFromWkt(self._equi7_data[subgrid_id]["extent"])
        geo_sr = osr.SpatialReference()
        geo_sr.SetWellKnownGeogCS("EPSG:4326")
        geom.AssignSpatialReference(geo_sr)
        return geom

    def get_subgrid_projection(self, subgrid_id):
        """return sub-grid spatial reference in wkt format

        Parameters
        ----------
        subgrid_id : string
            sub-grid id string e.g. EU for Europe

        Return
        ------
        string
            wkt string representing sub-grid spatial reference

        """
        return self._equi7_data[subgrid_id]["project"]

    def get_subgrid_tiles(self, subgrid_id, tilecode=None):
        """return all the tiles in the given sub grid"""
        if tilecode is None:
            tilecode = self.tilecode

        return list(self._equi7_data[subgrid_id]["coverland"][tilecode])

    def identify_subgrid(self, geom):
        """return the overlapped grid ids."""
        subgrid_ids = [subgrid_id for subgrid_id in Equi7Grid._static_subgrid_ids
                     if geom.Intersects(self.get_subgrid_zone_geom(subgrid_id))]
        return subgrid_ids

    def identify_tile_per_xy(self, subgrid_id, location):
        """Return the tile name."""
        east = (int(location[0])
                / self.tile_size_m) * self.tile_size_m / 100000
        north = (int(location[1])
                 / self.tile_size_m) * self.tile_size_m / 100000
        return "{}{:03d}M_E{:03d}N{:03d}{}".format(subgrid_id, self.res,
                                                   east, north, self.tilecode)

    @staticmethod
    def get_tile_extent_m(ftile):
        """return the extent_m of the tile in the terms of [minX,minY,maxX,maxY]
        """
        return Equi7Tile(ftile).extent_m

    @staticmethod
    def get_tile_geotransform(ftile):
        """
        Return the GDAL geotransform list

        Parameters
        ----------
        ftile : string
            full tile name e.g. EU075M_E048N012T6

        Returns
        -------
        list
            a list contain the geotransfrom elements

        """
        return Equi7Tile(ftile).geotransform()

    def get_tile_geotags(self, ftile):
        """
        Return the geotags for given tile used as geoinformation for Geotiff
        """
        tile = Equi7Tile(ftile)
        geotags = {'geotransform': tile.geotransform(),
                   'spatialreference': tile.projection()}
        return geotags

    def lonlat2equi7xy(self, lon, lat, subgrid_id=None):
        """ convert (lon, lat) to Equi7 (subgrid_id, x, y)

        Parameters
        ----------
        lon: number or numpy.ndarray
            longitude
        lat: number or numpy.ndarray
            latitude
        subgrid_id: string, optional
            subgrid_id if known, otherwise it will be found from the coordinates.
            Without a known subgrid_id this method is much slower since each
            point has to be handled separately.

        Return
        ------
        subgrid_id: string
            subgrid id
        x: like longitude
            projected x coordinate
        y: like latitude
            projected y coordinate

        Raises
        ------
        TypeError
            if lon, lat are numpy arrays but no subgrid_id is given
        """
        if subgrid_id is None:

            vfunc = np.vectorize(self._lonlat2xy_no_subgrid)
            return vfunc(lon, lat)
        else:
            return self._lonlat2xy_subgrid(lon, lat, subgrid_id=subgrid_id)

    def _lonlat2xy_no_subgrid(self, lon, lat):
        """ convert (lon, lat) to Equi7 (subgrid_id, x, y)
        without knowledge of subgrid

        Parameters
        ----------
        lon: number
            longitude
        lat: number
            latitude

        Return
        ------
        subgrid_id: string
            subgrid id
        x: number
            projected x coordinate
        y: number
            projected y coordinate
        """

        # create point geometry
        geo_sr = osr.SpatialReference()
        geo_sr.SetWellKnownGeogCS("EPSG:4326")
        geom_pt = ogr.Geometry(ogr.wkbPoint)
        geom_pt.AddPoint(lon, lat)
        geom_pt.AssignSpatialReference(geo_sr)
        # any valid location should be in only one grid
        subgrid_id = self.identify_subgrid(geom_pt)[0]
        # convert to Equi7
        grid_sr = osr.SpatialReference()
        grid_sr.ImportFromWkt(self.get_subgrid_projection(subgrid_id))
        tx = osr.CoordinateTransformation(geo_sr, grid_sr)
        x, y, _ = tx.TransformPoint(lon, lat)

        return np.full_like(x, subgrid_id, dtype=object), x, y

    def _lonlat2xy_subgrid(self, lon, lat, subgrid_id):
        """ convert (lon, lat) to Equi7 (subgrid_id, x, y)
        with knowledge of subgrid id

        Parameters
        ----------
        lon: number or numpy.ndarray
            longitude
        lat: number or numpy.ndarray
            latitude
        subgrid_id: string
            subgrid id

        Return
        ------
        subgrid_id: string
            subgrid id
        x: like longitude
            projected x coordinate
        y: like latitude
            projected y coordinate
        """

        grid_sr = osr.SpatialReference()
        grid_sr.ImportFromWkt(self.get_subgrid_projection(subgrid_id))
        p_grid = pyproj.Proj(grid_sr.ExportToProj4())
        p_geo = pyproj.Proj(init="EPSG:4326")
        x, y = pyproj.transform(p_geo, p_grid, lon, lat)

        return subgrid_id, x, y

    def equi7xy2lonlat(self, subgrid_id, x, y):
        """ convert Equi7 (subgrid_id, x, y) to (lon, lat)

        Parameters
        ----------
        subgrid_id : string
            the sub-grid id (represent continent) of the coordination,
            should be one of ["NA", "EU", "AS", "SA", "AF", "OC", "AN"]
        x : number or numpy.ndarray
            x coordination
        y : number or numpy.ndarray
            y coordination
        Return
        ------
        lon: number of numpy.ndarray
            longitudes in EPSG 4326
        lat: number of numpy.ndarray
            latitudes in EPSG 4326
        """
        grid_sr = osr.SpatialReference()
        grid_sr.ImportFromWkt(self.get_subgrid_projection(subgrid_id))
        p_grid = pyproj.Proj(grid_sr.ExportToProj4())
        p_geo = pyproj.Proj(init="EPSG:4326")
        lon, lat = pyproj.transform(p_grid, p_geo, x, y)
        return lon, lat

    def equi7xy2en(self, subgrid_id, x, y, projstring=None, epsg='4326'):
        '''
        todo: make new method and replace it here. keep the name here?
        e.g.
        equi7xy2en(self, subgrid_id, x, y, projstring=None, epsg='4326')
            xy2uv(subgrid_id, x, y, projstring=None, epsg='4326')
        '''

        """ convert Equi7 (subgrid_id, x, y) to (easting, northing)
        in output projection

        Parameters
        ----------
        subgrid_id : string
            the sub-grid id (represent continent) of the coordination,
            should be one of ["NA", "EU", "AS", "SA", "AF", "OC", "AN"]
        x : number or numpy.ndarray
            x coordination
        y : number or numpy.ndarray
            y coordination
        Return
        ------
        lon: number of numpy.ndarray
            easting in output projection
        lat: number of numpy.ndarray
            northing in output projection
        """
        grid_sr = osr.SpatialReference()
        grid_sr.ImportFromWkt(self.get_subgrid_projection(subgrid_id))
        p_grid = pyproj.Proj(grid_sr.ExportToProj4())
        if projstring is None:
            projection = ':'.join(['EPSG', epsg])
            p_geo = pyproj.Proj(init=projection)
        else:
            projection = projstring
            p_geo = pyproj.Proj(projection)
        eastings, northings = pyproj.transform(p_grid, p_geo, x, y)
        return eastings, northings

    @staticmethod
    def decode_tilename(ftile):
        """keep for compatibility"""
        # example ftile "AF500M_E060N036T6"
        # check the validation of ftile
        subgrid_id = ftile[0:2]
        res = int(ftile[2:5])
        tile_size_m = int(ftile[-1]) * 100000
        east = int(ftile[8:11]) * 100000
        north = int(ftile[12:15]) * 100000
        tile_code = ftile[-2:]

        if int(res) not in Equi7Grid._static_res:
            raise ValueError("The grid resolution is not supported!")

        equi7_data = Equi7Grid._static_equi7_data
        if ftile[7:] not in equi7_data[subgrid_id]["coverland"][tile_code]:
            msg = "Given tile name was not found in land masses tile list!"
            raise ValueError(msg)

        return (subgrid_id, res, tile_size_m, east, north, tile_code)

    @staticmethod
    def get_index(dst_ftile, src_ftile, get_px_counts=False):
        """
        Return the index for oversammpling src ftile to dst ftile.

        Parameters
        ----------
        dst_ftile : string
            dst full tile name e.g. AF075M_E066N030T6
        src_ftile : string
            src full tile name e.g. AF500M_E066N030T6
        get_px_counts : bool
            keyword for giving as second return output the number of
            fine pixels in individual coarse pixels

        Return
        ------
        index : numpy array
            The index array with the same size_px as the dst tilename
        px_counts : numpy array
            The number number of fine pixels per coarse pixel
        """

        dtile = Equi7Tile(dst_ftile)
        stile = Equi7Tile(src_ftile)

        # check if dst tile is a sub tile of src tile
        if dtile.llx < stile.llx or dtile.lly < stile.lly \
           or dtile.llx + dtile.size_m > stile.llx + stile.size_m \
           or dtile.lly + dtile.size_m > stile.lly + stile.size_m:
            raise ValueError("dst tile should be a sub tile of src tile!")

        index_pattern = {"075T6-500T6": (7, 6, 7),
                         "040T3-500T6": (12, 13),
                         "010T1-500T6": (50,),
                         "040T3-075T6": (2, 2, 2, 1, 2, 2, 2, 2),
                         "010T1-075T6": (7, 8),
                         "010T1-040T3": (4,),
                         "010T1-020T1": (4,)
                         }

        index_id = "{:03d}{}-{:03d}{}".format(dtile.res, dtile.tilecode,
                                              stile.res, stile.tilecode)
        if index_id in index_pattern:
            pattern = index_pattern[index_id]
            pattern_sum = int(np.sum(pattern))
        else:
            raise ValueError("Unsupported indexing!")

        # create template
        pattern_tmpl = list()
        for i in range(len(pattern)):
            pattern_tmpl.extend([i] * pattern[i])
        # k number of patterns that dtile takes
        k = dtile.size_px / pattern_sum + 2
        idx = np.tile(pattern_tmpl, k)
        corr = np.repeat(np.arange(k) * len(pattern), pattern_sum)
        idx += corr

        # make x index
        xoff = (dtile.llx - stile.llx) / dtile.res
        # x_n skip number of patterns
        # x_m start index of dtile in the remaining patterns
        x_n, x_m = xoff / pattern_sum, xoff % pattern_sum
        x_idx = idx + (x_n * len(pattern))
        # shift idx to the correct start point
        x_idx = x_idx[x_m:x_m + dtile.size_px]

        # make y index
        yoff = (dtile.lly + dtile.size_m - stile.lly - stile.size_m) / -dtile.res
        y_n, y_m = yoff / pattern_sum, yoff % pattern_sum
        y_idx = idx + (y_n * len(pattern))
        # shift idx to the correct start point
        y_idx = y_idx[y_m:y_m + dtile.size_px]

        # create index array
        index = np.empty((dtile.size_px, dtile.size_px), dtype=np.uint32)
        for i, v in enumerate(y_idx):
            index[i, :] = x_idx + v * stile.size_px

        if get_px_counts:
            n_pixels = (np.unique(x_idx,return_counts=True)[1]).astype(np.uint16)
            return index, n_pixels
        else:
            return index

    @staticmethod
    def link_res_2_tile(res, get_span=False):
        res = int(res)
        tile_code = None
        tile_size_m = None
        # allowing sampling of [1000, 800, 750, 600, 500, 400, 300, 250, 200, 150, 125, 100, 96, 80, 75, 64] metres
        if ((res in range(64, 1001)) and (600000 % res == 0)):
            tile_code = "T6"
            tile_size_m = 600000
        # allowing sampling of [60, 50, 48, 40, 32, 30, 25, 24, 20] metres
        elif ((res in range(20, 61)) and (300000 % res == 0)):
            tile_code = "T3"
            tile_size_m = 300000
        # allowing sampling of [16, 10, 8, 5, 4, 2, 1] metres
        elif ((res in range(1, 17)) and (100000 % res == 0)):
            tile_code = "T1"
            tile_size_m = 100000
        else:
            msg = "Error: Given resolution %d is not supported!" % res
            msg += " Supported resolutions: {}".format(
                str(Equi7Grid._static_res))
            raise ValueError(msg)

        if get_span == True:
            result = (tile_code, tile_size_m)
        else:
            result = tile_code

        return result

    @staticmethod
    def find_overlapped_tiles(ftile, res):
        """ This function will return the corresponding tile of the
        given tile in the given resolution grid system.

        Parameters
        ----------
        ftile : string
            full tile name e.g. AF075M_E066N030T6
        res : int
            the resolution of the grid system

        Return
        ------
        list
            list of found tiles.
        """
        return Equi7Tile(ftile).find_family_tiles(res=res)

    def identify_tiles_per_bbox(self, subgrid_id, bbox):
        """Light-weight routine that returns
           the name of tiles intersecting the bounding box.

        Parameters
        ----------
        subgrid_id : string
            sub-grid id string e.g. EU for Europe
        extent_m : list
            list of equi7-coordinates limiting the bounding box.
            scheme: [xmin, ymin, xmax, ymax]

        Return
        ------
        tilenames : list
            list of tile names intersecting the bounding box,

        """
        xmin, ymin, xmax, ymax = bbox

        if (xmin >= xmax) or (ymin >= ymax):
            raise ValueError("Check order of coordinates of bbox! "
                             "Scheme: [xmin, ymin, xmax, ymax]")

        tcode = self.tilecode
        trs = self.res

        sp = self.tile_size_m
        sp_multi = sp / 100000

        x_anchors_100km = range(xmin / sp * sp_multi, xmax / sp * sp_multi + 1, sp_multi)
        y_anchors_100km = range(ymin / sp * sp_multi, ymax / sp * sp_multi + 1, sp_multi)

        tx, ty = np.meshgrid(x_anchors_100km, y_anchors_100km)
        tx = tx.flatten()
        ty = ty.flatten()

        tilenames = []
        for i, _ in enumerate(tx):
            tilenames.append("{}{:03d}M_E{:03d}N{:03d}{}".format(subgrid_id, trs, tx[i], ty[i], tcode))

        return tilenames

    def get_tiles_per_bbox(self, subgrid_id, bbox):
        """Light-weight routine that returns
           the name of tiles intersecting the bounding box.

        Parameters
        ----------
        subgrid_id : string
            sub-grid id string e.g. EU for Europe
        bbox : list
            list of equi7-coordinates limiting the bounding box.
            scheme: [xmin, ymin, xmax, ymax]

        Return
        ------
        e7tiles : list
            list of Equi7Tiles() intersecting the bounding box,
            with .subset() not exceeding the bounding box.

        """
        tilenames = self.identify_tiles_per_bbox(subgrid_id, bbox)

        e7tiles = list()

        for t in tilenames:

            tile = Equi7Tile(t)

            le, te, re, be = tile.active_subset_px
            extent = tile.extent_m

            #left_edge
            if extent[0] <= bbox[0]:
                le = (bbox[0] - extent[0]) / tile.res
            #top_edge
            if extent[1] <= bbox[1]:
                te = (bbox[1] - extent[1]) / tile.res
            #right_edge
            if extent[2] > bbox[2]:
                re = (bbox[2] - extent[2] + tile.size_m) / tile.res
            #bottom_edge
            if extent[3] > bbox[3]:
                be = (bbox[3] - extent[3] + tile.size_m) / tile.res

            tile.active_subset_px = le, te, re, be

            e7tiles.append(tile)

        return e7tiles

    def search_tiles(self,
                     geom_area=None,
                     extent_m=None,
                     epsg=4326,
                     subgrid_ids=None,
                     coverland=False,
                     gdal_path=None):

        """
        Search the tiles which are intersected by the poly_roi area.

        Parameters
        ----------
        geom_area : geometry
            a polygon or multipolygon geometery object representing the ROI
        extent_m : list
            It is a polygon representing the rectangle region of interest
            in the format of [xmin, ymin, xmax, ymax].
        epsg : str
            EPSG CODE defining the spatial reference system, in which
            the geometry or extent_m is given. Default is LatLon (EPSG:4326)
        subgrid_ids : list
            grid ID to specified which continents you want to search. Default
            value is None for searching all continents.

        Returns
        -------
        list
            return a list of  the overlapped tiles' name.
            If not found, return empty list.
        """
        # check input grids
        if subgrid_ids is None:
            subgrid_ids = Equi7Grid._static_subgrid_ids
        else:
            subgrid_ids = [x.upper() for x in subgrid_ids]
            if set(subgrid_ids).issubset(set(Equi7Grid._static_subgrid_ids)):
                subgrid_ids = list(subgrid_ids)
            else:
                raise ValueError("Invalid agrument: grid must one of [ %s ]." %
                                 " ".join(Equi7Grid._static_subgrid_ids))

        if not geom_area and not extent_m:
            print("Error: either geom or extent_m should be given as the ROI.")
            return list()

        # obtain the geometry of ROI
        if not geom_area:
            geom_area = gdalport.extent2polygon(extent_m)
            geom_sr = osr.SpatialReference()
            geom_sr.ImportFromEPSG(epsg)
            geom_area.AssignSpatialReference(geom_sr)

        # load lat-lon spatial reference as the default
        geo_sr = osr.SpatialReference()
        geo_sr.ImportFromEPSG(4326)

        geom_sr = geom_area.GetSpatialReference()
        if geom_sr is None:
            geom_area.AssignSpatialReference(geo_sr)
        elif not geom_sr.IsSame(geo_sr):
            geom_area.TransformTo(geo_sr)

        # intersect the given grid ids and the overlapped ids
        overlapped_grids = self.identify_subgrid(geom_area)
        subgrid_ids = list(set(subgrid_ids) & set(overlapped_grids))

        # finding tiles
        overlapped_tiles = list()
        for subgrid_id in subgrid_ids:
            overlapped_tiles.extend(
                self.__search_subgrid_tiles(geom_area, subgrid_id, coverland))
        return overlapped_tiles

    def __search_subgrid_tiles(self, area_geometry, subgrid_id, coverland):
        """
        Search the tiles which are overlapping with the given grid.

        Parameters
        ----------
        area_geometry : geometry
            It is a polygon geometry representing the region of interest.
        subgrid_id : string
            sub grid ID to specified which continent you want to search.
            Default value is None for searching all continents.

        Returns
        -------
        list
            Return a list of  the overlapped tiles' name.
            If not found, return empty list.
        """
        # get the intersection of the area of interest and grid zone
        intersect = area_geometry.Intersection(self.get_subgrid_zone_geom(subgrid_id))
        if not intersect:
            return list()
        # The spatial reference need to be set again after intersection
        intersect.AssignSpatialReference(area_geometry.GetSpatialReference())
        # transform the area of interest to the grid coordinate system
        grid_sr = osr.SpatialReference()
        grid_sr.ImportFromWkt(self.get_subgrid_projection(subgrid_id))
        intersect.TransformTo(grid_sr)

        # get envelope of the Geometry and cal the bounding tile of the
        envelope = intersect.GetEnvelope()
        x_min = int(envelope[0]) / self.tile_size_m * self.tile_size_m
        x_max = (int(envelope[1]) / self.tile_size_m + 1) * self.tile_size_m
        y_min = int(envelope[2]) / self.tile_size_m * self.tile_size_m
        y_max = (int(envelope[3]) / self.tile_size_m + 1) * self.tile_size_m

        # make sure x_min and y_min greater or equal 0
        x_min = 0 if x_min < 0 else x_min
        y_min = 0 if y_min < 0 else y_min

        # get overlapped tiles
        overlapped_tiles = list()
        for x, y in itertools.product(range(x_min, x_max, self.tile_size_m),
                                      range(y_min, y_max, self.tile_size_m)):
            geom_tile = gdalport.extent2polygon((x, y, x + self.tile_size_m,
                                                 y + self.tile_size_m))
            if geom_tile.Intersects(intersect):
                ftile = self.identify_tile_per_xy(subgrid_id, [x, y])
                if not coverland or self.is_coverland(ftile):
                    overlapped_tiles.append(ftile)

        return overlapped_tiles

    def resample_to_geotiff(self, image, output_dir, gdal_path=None, subgrid_ids=None,
                            accurate_boundary=True, e7_folder=True, ftiles=None,
                            roi=None, coverland=True, outshortname=None,
                            withtilenameprefix=False, withtilenamesuffix=True,
                            compress=True, compresstype="LZW", resampling_type="near",
                            overwrite=False, image_nodata=None, tile_nodata=None,
                            tiledtiff=True, blocksize=512):
        """Resample the image to tiles.

        Parameters
        ----------
        image : string
            Image file path.
        ouput_dir : string
            Output directory path.
        e7_folder: boolean
            if set (default), the output data will be stored in equi7 folder structure
        gdal_path : string
            Gdal utilities location.
        subgrid_ids : list
            Only resample to the specified continents,
            default is to resample to all continents.
        roi : OGRGeometry
            Region of interest.
            The roi will beignored if ftiles keyword is provided
        ftiles : list of tile names
            full name of tiles to which data should be resampled
        coverland : bool
            Only resample to the tiles that cover land.
        outshortname : string
            The short name will be main part of the output tile name.
        withtilenameprefix : bool
            Prepend the tile name in the tile file name.
        withtilenamesuffix : bool
            Append the tile name in the tile file name.
        compress : bool
            The output tiles is compressed or not.
        resampling_type : string
            Resampling method.
        overwrite : bool
            Overwrite the old tile or not.
        image_nodata : double
            The nodata value of input images.
        tile_nodata : double
            The nodata value of tile.
        tiledtiff : bool
            Set to yes for tiled geotiff output
        blocksize: integer
            sets tile width, in x and y direction

        """
        # get the output shortname if not provided
        if not outshortname:
            outshortname = os.path.splitext(os.path.basename(image))[0]

        # find overlapping tiles
        if ftiles is None:
            if roi is None:
                if accurate_boundary:
                    try:
                        roi_geom = retrieve_raster_boundary(image,
                                                            gdal_path=gdal_path,
                                                            nodata=image_nodata)
                    except Exception as e:
                        print("retrieve_raster_boundary failed:", str(e))
                        roi_geom = None
                else:
                    roi_geom = None
                if roi_geom:
                    ftiles = self.search_tiles(geom_area=roi_geom,
                                               subgrid_ids=subgrid_ids,
                                               coverland=coverland,
                                               gdal_path=gdal_path)
                else:
                    ds = gdalport.open_image(image)
                    img_extent = ds.get_extent()
                    geo_extent = gdalport.extent2polygon(img_extent)
                    geo_sr = osr.SpatialReference()
                    geo_sr.ImportFromWkt(ds.projection())
                    geo_extent.AssignSpatialReference(geo_sr)
                    ftiles = self.search_tiles(geom_area=geo_extent,
                                               subgrid_ids=subgrid_ids,
                                               coverland=coverland,
                                               gdal_path=gdal_path)
            else:
                ftiles = self.search_tiles(geom_area=roi,
                                           subgrid_ids=subgrid_ids,
                                           coverland=coverland,
                                           gdal_path=gdal_path)
        else:
            if type(ftiles) != list:
                ftiles = [ftiles]

        # keep the full path of the output(resampled) files
        dst_file_names = []

        # resample for each tile sequentially
        for ftile in ftiles:
            # create grid folder
            if e7_folder:
                grid_folder = "EQUI7_{}".format(ftile[0:6])
                tile_path = os.path.join(output_dir, grid_folder, ftile[7:])
                if not os.path.exists(tile_path):
                    platform.makedirs(tile_path)
            else:
                tile_path = output_dir

            # make output filename
            outbasename = outshortname
            if withtilenameprefix:
                outbasename = "_".join((ftile, outbasename))
            if withtilenamesuffix:
                outbasename = "_".join((outbasename, ftile))
            filename = os.path.join(tile_path, "".join((outbasename, ".tif")))

            # using gdalwarp to resample
            extent_m = self.get_tile_extent_m(ftile)
            tile_project = self.get_subgrid_projection(ftile[0:2])

            # prepare options for gdalwarp
            options = {'-t_srs': tile_project, '-of': 'GTiff',
                       '-r': resampling_type,
                       '-te': " ".join(map(str, extent_m)),
                       '-tr': "{} -{}".format(self.res, self.res)}

            options["-co"] = list()
            if compress:
                options["-co"].append("COMPRESS={0}".format(compresstype))
            if image_nodata != None:
                options["-srcnodata"] = image_nodata
            if tile_nodata != None:
                options["-dstnodata"] = tile_nodata
            if overwrite:
                options["-overwrite"] = " "
            if tiledtiff:
                options["-co"].append("TILED=YES")
                options["-co"].append("BLOCKXSIZE={0}".format(blocksize))
                options["-co"].append("BLOCKYSIZE={0}".format(blocksize))

            # call gdalwarp for resampling
            succeed, _ = gdalport.call_gdal_util('gdalwarp', src_files=image,
                                                 dst_file=filename, gdal_path=gdal_path,
                                                 options=options)

            # full path to the output(resampled) files
            if succeed:
                dst_file_names.extend([filename])

        return dst_file_names


class Equi7Tile(object):

    """ Equi7 Tile class

    A tile in the Equi7 grid system.
    """

    def __init__(self, ftile):
        """ fetch information from full tile name
        """
        self.ftile = None
        if not Equi7Tile.is_valid(ftile):
            raise ValueError("invalid full tile name!")
        self.ftile = ftile.upper()
        self.subgrid = self.ftile[0:2]
        self.res = int(self.ftile[2:5])
        self.tilecode = self.ftile[-2:]
        self.llx = int(self.ftile[8:11]) * 100000
        self.lly = int(self.ftile[12:15]) * 100000
        self.size_m = int(self.ftile[-1]) * 100000
        self.size_px = int(self.size_m / self.res)
        self._subset_px = (0, 0, self.size_px, self.size_px)


    @property
    def fullname(self):
        return self.ftile


    @property
    def shortname(self):
        return self.ftile[7:]


    @property
    def extent_m(self):
        """return the extent_m of the tile in the terms of [minX,minY,maxX,maxY]
        """
        return [self.llx, self.lly,
                self.llx + self.size_m, self.lly + self.size_m]

    @property
    def shape_px(self):
        return (self.size_px, self.size_px)


    @property
    def active_subset_px(self):
        """
        holds indices of the active_subset_px-of-interest
        :return: active_subset_px-of-interest
        """
        return self._subset_px

    @active_subset_px.setter
    def active_subset_px(self, limits):
        """
        changes the indices of the active_subset_px-of-interest,
        mostly to a smaller extent, for efficient reading

        limits : tuple
            the limits of subsets as (xmin, ymin, xmax, ymax).

        """

        string = ['xmin', 'ymin', 'xmax', 'ymax']
        if len(limits) != 4:
            raise ValueError('Limits are not properly set!')

        _max = self.size_px

        for l, limit in enumerate(limits):
            if (limit < 0) or (limit > _max):
                raise ValueError('{} is out of bounds!'.format(string[l]))

        xmin, ymin, xmax, ymax = limits

        if xmin >= xmax:
            raise ValueError('xmin >= xmax!')
        if ymin >= ymax:
            raise ValueError('ymin >= ymax!')

        self._subset_px = limits


    def projection(self):
        return Equi7Grid._static_equi7_data[self.subgrid]["project"]

    def geotransform(self):
        """
        Return the GDAL geotransform list

        Parameters
        ----------
        ftile : string
            full tile name e.g. EU075M_E048N012T6

        Returns
        -------
        list
            a list contain the geotransfrom elements

        """
        geot = [self.llx, self.res, 0,
                self.lly + self.size_m, 0, -self.res]
        return geot


    def px_2_coord(self, x, y):
        """
        Returns the Equi7 coordinates of a tile pixel

        Parameters
        ----------
        x : number
            pixel row number
        y : number
            pixel collumn number


        Returns
        -------
        tuple
            Equi7 coordinates of a tile pixel

        """
        gt = self.geotransform()

        xgeo = gt[0] + x * gt[1] + y * gt[2]
        ygeo = gt[3] + x * gt[4] + y * gt[5]

        return (xgeo, ygeo)

    def get_tile_geotags(self):
        """
        Return geotags for given tile used as geoinformation for Geotiff
        """
        geotags = {'geotransform': self.geotransform(),
                   'spatialreference': self.projection()}

        return geotags

    def find_family_tiles(self, res=None, target_tilecode=None):
        """find the family tiles which share the same extent_m but
        with different resolution and tilecode

        Parameters
        ----------
        ftile : string
            full tile name e.g. AF075M_E066N030T6
        res : int
            the resolution of the grid system
        tilecode : string
            tile code string

        Returns
        -------
        list
            list of found tiles.

        Notes
        -----
        Either the res or tile code should be given.
        But if both are given, the res will be used.
        """
        if res is not None:
            target_tilecode = Equi7Grid.link_res_2_tile(res)
        elif target_tilecode is not None:
            target_tilecode = target_tilecode
        else:
            raise ValueError("either res or tilecode should be given!")

        # found family tiles
        family_tiles = list()

        if target_tilecode >= self.tilecode:
            t_span = int(target_tilecode[-1]) * 100000
            t_east = (self.llx / t_span) * t_span / 100000
            t_north = (self.lly / t_span) * t_span / 100000
            name = "E{:03d}N{:03d}{}".format(t_east, t_north, target_tilecode)
            family_tiles.append(name)
        else:
            sub_span = int(target_tilecode[-1]) * 100000
            n = int(self.size_m / sub_span)
            for x, y in itertools.product(range(n), range(n)):
                s_east = (self.llx + x * sub_span) / 100000
                s_north = (self.lly + y * sub_span) / 100000
                name = "E{:03d}N{:03d}{}".format(s_east, s_north, target_tilecode)
                family_tiles.append(name)
        return family_tiles

    def covers_land(self):
        """check if tile covers land"""
        land_tiles = Equi7Grid._static_equi7_data[self.subgrid]["coverland"]
        return self.shortname in land_tiles[self.tilecode]

    @staticmethod
    def create(equi7xy=None, lonlat=None):
        raise NotImplementedError()

    @staticmethod
    def is_valid(ftile):
        """check if ftile is a valid tile name"""

        ftile = ftile.upper()
        # check the constant
        if len(ftile) != 17:
            return False
        if ftile[5:8] != "M_E" or ftile[11] != "N":
            return False
        # check variables
        if ftile[0:2] not in Equi7Grid._static_subgrid_ids:
            return False
        if int(ftile[2:5]) not in Equi7Grid._static_res:
            return False
        _east = int(ftile[8:11])
        _north = int(ftile[12:15])
        if _east < 0 or _north < 0:
            return False
        if ftile[-2:] not in Equi7Grid._static_tilecodes:
            return False
        return True


@click.command()
@click.argument('resolution')
@click.option('--subgrid_ids', '-s',
              default=None,
              multiple=True,
              type=click.Choice(['NA', 'EU', 'AS',
                                 'SA', 'AF', 'OC', 'AN']),
              help=('Equi7 subgrid id, '
                    'Default are all subgrids'))
@click.option('--extent_m', '-e',
              default=(-180.0, -90.0, 180.0, 90.0),
              type=float,
              nargs=4,
              help='extent_m to search tiles in. lonmin latmin lonmax latmax')
@click.option('--coverland', '-c',
              type=bool,
              default=True,
              help='Only search tile that cover land area.')

def find_tiles(resolution, subgrid_ids, extent_m, coverland):
    """
    Find Tiles for Equi7 grid of resolution.
    """
    eq7 = Equi7Grid(resolution)
    if len(subgrid_ids) == 0:
        subgrid_ids = None
    for tile in eq7.search_tiles(subgrid_ids=subgrid_ids,
                                 extent_m=extent_m,
                                 coverland=coverland):
        click.echo(tile)


# add a new subgroup for equi7 related commands
@click.group()
def equi7():
    pass

equi7.add_command(find_tiles)
cli.add_command(equi7)
