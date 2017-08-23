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
Created on March 1, 2017

Code for Tiled Projection Systems.

@author: Bernhard Bauer-Marschallinger, bbm@geo.tuwien.ac.at

'''

import abc

import numpy as np
import dask.array as da

from osgeo import osr

import geometry as geometry



class TPSCoreProperty(object):
    """
    Class holding information needed at every level of `TiledProjectionSystem`,
    the alltime-valid "core properties".
    With this, core parameters are everywhere accessible per same name.

    Parameters
    ----------
    epsg : integer
        sdfsd
    """
    def __init__(self, tag, projection, res, tiletype, tile_xsize_m, tile_ysize_m):
        self.tag = tag
        self.projection = projection
        self.res = res
        self.tiletype = tiletype
        self.tile_xsize_m = tile_xsize_m
        self.tile_ysize_m = tile_ysize_m


class TPSProjection():
    """
    Projection class holding and translating the definitions of a projection when initialising.

    Parameters
    ----------
    epsg : integer
        The EPSG-code of the spatial reference. As from http://www.epsg-registry.org
        Not all reference do have a EPSG code.
    proj4 : string
        The proj4-string defining the spatial reference.
    wkt : string
        The wkt-string (well-know-text) defining the spatial reference.

    """
    def __init__(self, epsg=None, proj4=None, wkt=None):

        checker = {epsg, proj4, wkt}
        checker.discard(None)
        if len(checker) == 0:
            raise ValueError('Projection is not defined!')

        if len(checker) != 1:
            raise ValueError('Projection is defined ambiguously!')

        spref = osr.SpatialReference()

        if epsg is not None:
            spref.ImportFromEPSG(epsg)
            self.osr_spref = spref
            self.proj4 = spref.ExportToProj4()
            self.wkt = spref.ExportToWkt()
            self.epsg = epsg

        if proj4 is not None:
            spref.ImportFromProj4(proj4)
            self.osr_spref = spref
            self.proj4 = proj4
            self.wkt = spref.ExportToWkt()
            self.epsg = self.extract_epsg(self.wkt)

        if wkt is not None:
            spref.ImportFromWkt(wkt)
            self.osr_spref = spref
            self.proj4 = spref.ExportToProj4()
            self.wkt = wkt
            self.epsg = self.extract_epsg(self.wkt)

    def extract_epsg(self, wkt):
        """
        Checks if the WKT contains an EPSG code for the spatial reference, a
        and returns it, if found.

        Parameters
        ----------
        wkt : string
            The wkt-string (well-know-text) defining the spatial reference.

        Return
        ------
        epsg : integer, None
            the EPSG code of the spatial reference (if found). Else: None

        """
        pos_last_code = wkt.rfind('EPSG')
        pos_end = len(wkt)
        if pos_end - pos_last_code < 16:
            epsg = int(wkt[pos_last_code+7:pos_last_code+11])
        else:
            epsg = None

        return epsg


def dummy(self, thing):
    """
    Description

    Parameters
    ----------
    thing : type
        more words

    Return
    ------
    thing : type
        more words
    """

    return thing

class TiledProjectionSystem(object):

    __metaclass__ = abc.ABCMeta

    def __init__(self, res, nametag='TPS'):

        tiletype = self.get_tiletype(res)
        tile_xsize_m, tile_ysize_m = self.get_tilesize(res)

        self.core = TPSCoreProperty(nametag, None, res, tiletype, tile_xsize_m, tile_ysize_m)

        self.subgrids = self.define_subgrids()
        pass

    def __getattr__(self, item):
        '''
        short link for items of subgrids and core
        '''
        if item in self.subgrids:
            return self.subgrids[item]
        elif item in self.core.__dict__:
            return self.core.__dict__[item]
        else:
            return self.__dict__[item]


    @abc.abstractmethod
    def define_subgrids(self):
        pass


    def locate_geometry_in_subgrids(self, geom):
        """
        finds overlapping subgrids of given geometry.
        """
        #covering_subgrid = dict()
        covering_subgrid = list()
        for x in self.subgrids.keys():
            if geom.Intersects(self.subgrids.get(x).polygon_geog):
                #covering_subgrid[x] = self.subgrids.get(x)
                covering_subgrid.append(x)
        return covering_subgrid


    def lonlat2xy(self, lat, lon, subgrid=None):
        '''
        converts latitude and longitude coordinates to TPS grid coordinates

        :param lat:
        :param lon:
        :param subgrid:
        :return:
        '''

        # TODO use pyproj.transform for performance when having many points
        if subgrid is None:
            vfunc = np.vectorize(self._lonlat2xy)
            return vfunc(lat, lon)
        else:
            vfunc = np.vectorize(self._lonlat2xy_subgrid)
            return vfunc(lat, lon, subgrid)


    def _lonlat2xy(self, lat, lon):
        """
        finds overlapping subgrids of given geometry and computes the projected coordinates
        """
        # create point geometry
        lonlatprojection = TPSProjection(epsg=4326)
        point_geom = geometry.create_point_geom(lon, lat, lonlatprojection)

        # search for co-locating subgrid
        subgrid = self.locate_geometry_in_subgrids(point_geom)[0]

        x, y, = geometry.uv2xy(lonlatprojection.osr_spref,
                      self.subgrids[subgrid].core.projection.osr_spref,
                      lon,
                      lat)

        return np.full_like(x, subgrid, dtype=(np.str, len(subgrid))), x, y


    def _lonlat2xy_subgrid(self, lat, lon, subgrid):
        '''
        computes the projected coordinates in given subgrid

        :param lat:
        :param lon:
        :return:
        '''

        # set up spatial references
        lonlatprojection = TPSProjection(epsg=4326)

        x, y, = geometry.uv2xy(lonlatprojection.osr_spref,
                      self.subgrids[subgrid].core.projection.osr_spref,
                      lon,
                      lat)

        return np.full_like(x, subgrid, dtype=(np.str, len(subgrid))), x, y

    @abc.abstractmethod
    def create_tile(self, name):
        pass

    @abc.abstractmethod
    def get_tiletype(self, res):
        pass

    @abc.abstractmethod
    def get_tilesize(self, res):
        pass

class TiledProjection(object):
    """
    Class holding the projection and tiling definition of a tiled projection space.

    Parameters
    ----------
    Projection : Projection()
        A Projection object defining the spatial reference.
    tile_definition: TilingSystem()
        A TilingSystem object defining the tiling system.
        If None, the whole space is one single tile.
    """

    __metaclass__ = abc.ABCMeta

    staticdata = None

    def __init__(self, core, polygon_geog=None, tilingsystem=None):

        self.core = core

        if polygon_geog is None:
            polygon_geog = GlobalTile(self.core.projection).polygon()
        self.polygon_geog = polygon_geog
        self.polygon_proj= geometry.transform_geometry(polygon_geog, self.core.projection)
        self.bbox_geog = geometry.get_geom_boundaries(self.polygon_geog, rounding=self.core.res/1000000.0)
        self.bbox_proj = geometry.get_geom_boundaries(self.polygon_proj, rounding=self.core.res)

        if tilingsystem is None:
            tilingsystem = GlobalTile(self.core.projection)
        self.tilesys = tilingsystem


    def __getattr__(self, item):
        '''
        short link for items of core
        '''
        if item in self.core.__dict__:
            return self.core.__dict__[item]
        else:
            return self.__dict__[item]


    @abc.abstractmethod
    def get_polygon(self):
        return None


    def get_bbox_geog(self):
        bbox = self.polygon_geog.GetEnvelope()
        return bbox


    def xy2lonlat(self, x, y):
        vfunc = np.vectorize(self._xy2lonlat)
        return vfunc(x, y)

    def _xy2lonlat(self, x, y):
        '''
        convert projected coordinates to WGS84 longitude and latitude
        :param x:
        :param y:
        :return:
        '''

        # set up spatial references
        lonlatprojection = TPSProjection(epsg=4326)

        return geometry.uv2xy(self.core.projection.osr_spref, lonlatprojection.osr_spref, x, y)



class TilingSystem(object):
    """
    Class defining the tiling system and providing methods for queries and handling.

    Parameters (BBM: init(stuff))
    ----------
    projection : :py:class:`Projection`
        A Projection object defining the spatial reference.
    tile_definition: TilingSystem
        A TilingSystem object defining the tiling system.
        If None, the whole space is one single tile.

    Attributes (BBM: .stuff that needs to be explained)
    ----------
    extent_geog:
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, core, polygon_proj, x0, y0):

        self.core = core
        self.x0 = x0
        self.y0 = y0
        self.xstep = self.core.tile_xsize_m
        self.ystep = self.core.tile_ysize_m
        self.polygon_proj = polygon_proj
        self.bbox_proj = geometry.get_geom_boundaries(self.polygon_proj, rounding=self.core.res)


    def __getattr__(self, item):
        '''
        short link for items of core
        '''
        if item in self.core.__dict__:
            return self.core.__dict__[item]
        else:
            return self.__dict__[item]

    @abc.abstractmethod
    def create_tile(self, name=None, x=None, y=None):
        return

    def round_xy2lowerleft(self, x0, y0):
        llx = x0 / self.core.tile_xsize_m * self.core.tile_xsize_m
        lly = y0 / self.core.tile_ysize_m * self.core.tile_ysize_m
        return llx, lly

    @abc.abstractmethod
    def point2tilename(self, x, y):
        return

    @abc.abstractmethod
    def _encode_tilename(self, x, y):
        return

    @abc.abstractmethod
    def decode_tilename(self, name):
        a = None
        return a

    @abc.abstractmethod
    def identify_tiles_overlapping_bbox(self, bbox):
        """Light-weight routine that returns
           the name of tiles intersecting the bounding box.

        Parameters
        ----------
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

        tiletype = self.core.tiletype
        tres = self.core.res
        tilenames = list()
        return tilenames


    def create_tiles_overlapping_bbox(self, bbox):
        """Light-weight routine that returns
           the name of tiles intersecting the bounding box.

        Parameters
        ----------
        bbox : list
            list of equi7-coordinates limiting the bounding box.
            scheme: [xmin, ymin, xmax, ymax]

        Return
        ------
        tiles : list
            list of Equi7Tiles() intersecting the bounding box,
            with .subset() not exceeding the bounding box.

        """
        tilenames = self.identify_tiles_overlapping_bbox(bbox)
        tiles = list()

        for t in tilenames:

            tile = self.create_tile(name=t)
            le, te, re, be = tile.active_subset_px
            extent = tile.limits_m()

            # left_edge
            if extent[0] <= bbox[0]:
                le = (bbox[0] - extent[0]) / tile.core.res
            # top_edge
            if extent[1] <= bbox[1]:
                te = (bbox[1] - extent[1]) / tile.core.res
            # right_edge
            if extent[2] > bbox[2]:
                re = (bbox[2] - extent[2] + self.core.tile_xsize_m) / tile.core.res
            # bottom_edge
            if extent[3] > bbox[3]:
                be = (bbox[3] - extent[3] + self.core.tile_ysize_m) / tile.core.res

            tile.active_subset_px = le, te, re, be
            tiles.append(tile)

        return tiles


    def create_daskarray_overlapping_bbox(self, bbox):

        tiles = self.create_tiles_per_bbox(bbox)
        tilenames = [x.name for x in tiles]

        x_anchors = set([t.llx for t in tiles])
        y_anchors = set([t.lly for t in tiles])

        box = np.zeros(((bbox[2]-bbox[0])/self.core.res, (bbox[3]-bbox[1])/self.core.res))

        d = da.from_array(box, chunks=1000)
        for i in y_anchors:

            pass
        x = np.arange(1200 ** 2).reshape((1200, 1200))
        d = da.from_array(x, chunks=120)
        g = da.ghost.ghost(d, depth={0: 1, 1: 1}, boundary={0: da.concatenate, 1: da.concatenate})

        m = g.map_blocks(func)

        y = da.ghost.trim_internal(m, {0: 1, 1: 1})



def func(block):
    return np.mean(block)

class Tile(object):
    """
    Class defining a tile and providing methods for handling.

    Parameters
    ----------
    projection : :py:class:`Projection`
        A Projection object defining the spatial reference.

    Attributes (BBM: .stuff that needs to be explained)
    ----------
    extent_geog:
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, core, name, x0, y0):
        self.core = core
        self.name = name
        self.typename = core.tiletype
        self.llx = x0
        self.lly = y0
        self.x_size_px = self.core.tile_xsize_m / self.core.res
        self.y_size_px = self.core.tile_ysize_m / self.core.res
        self._subset_px = (0, 0, self.x_size_px, self.y_size_px)

    def __getattr__(self, item):
        '''
        short link for items of core
        '''
        if item in self.core.__dict__:
            return self.core.__dict__[item]
        else:
            return self.__dict__[item]


    def shape_px(self):
        """
        :returns the shape of the pixel array
        """
        return (self.x_size_px, self.y_size_px)

    def limits_m(self):
        """
        :returns the limits of the tile in the terms of (xmin, ymin, xmax, ymax)
        """
        return (self.llx, self.lly,
                self.llx + self.core.tile_xsize_m, self.lly + self.core.tile_ysize_m)

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

        _max = [self.x_size_px, self.y_size_px, self.x_size_px, self.y_size_px]

        for l, limit in enumerate(limits):
            if (limit < 0) or (limit > _max):
                raise ValueError('{} is out of bounds!'.format(string[l]))

        xmin, ymin, xmax, ymax = limits

        if xmin >= xmax:
            raise ValueError('xmin >= xmax!')
        if ymin >= ymax:
            raise ValueError('ymin >= ymax!')

        self._subset_px = limits

    def geotransform(self):
        """
        :returns the GDAL geotransform list

        Parameters
        ----------
        ftile : string
            full tile name e.g. EU075M_E048N012

        Returns
        -------
        list
            a list contain the geotransfrom elements

        """
        geot = [self.llx, self.res, 0,
                self.lly + self.core.tile_ysize_m, 0, -self.res]

        return geot

    def ij2xy(self, i, j):
        """
        Returns the projected coordinates of a tile pixel in the TilingSystem

        Parameters
        ----------
        i : number
            pixel row number
        j : number
            pixel collumn number

        Returns
        -------
        x : number
            x coordinate in the projection
        y : number
            y coordinate in the projection
        """

        gt = self.geotransform()

        x = gt[0] + i * gt[1] + j * gt[2]
        y = gt[3] + i * gt[4] + j * gt[5]

        if self.core.res <= 1.0:
            precision = len(str(int(1.0/self.core.res)))+1
            return round(x, precision), round(y, precision)
        else:
            return x, y

    def xy2ij(self, x, y):
        """
        returns the column and row number (i, j) of a projection coordinate (x, y)

        parameters
        ----------
        x : number
            x coordinate in the projection
        y : number
            y coordinate in the projection

        returns
        -------
        i : number
            pixel row number
        j : number
            pixel column number
        """

        gt = self.geotransform()

        #TODO: check if 1) round-to-nearest-integer or 2) round-down-to-integer
        i = int(round(-1.0*(gt[2] * gt[3] - gt[0] * gt[5] + gt[5] * x - gt[2] * y) /
                      (gt[2] * gt[4] - gt[1] * gt[5])))
        j = int(round(-1.0*(-1 * gt[1] * gt[3] + gt[0] * gt[4] - gt[4] * x + gt[1] * y) /
                      (gt[2] * gt[4] - gt[1] * gt[5])))

        return i, j

    def get_geotags(self):
        """
        Return geotags for given tile used as geoinformation for GDAL
        """
        geotags = {'geotransform': self.geotransform(),
                   'spatialreference': self.core.projection.wkt}

        return geotags

class GlobalTile(Tile):

    __metaclass__ = abc.ABCMeta

    def __init__(self, core, name):
        super(GlobalTile, self).__init__(core, name, 0, 0)




