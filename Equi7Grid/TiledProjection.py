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

import os
import abc
import copy
import pyproj

import numpy as np

from osgeo import osr
from osgeo import ogr


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

    def __init__(self, res, nametag='TPS'):

        tiletype, \
        tile_xsize_m, \
        tile_ysize_m = self.link_res_2_tilesize(res, get_size=True)

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

    #TODO: check how to force this correctly
    @abc.abstractmethod
    def define_subgrids(self):
        pass

    @abc.abstractmethod
    def link_res_2_tilesize(self):
        pass

    @abc.abstractmethod
    def latlon2xy(self):
        pass

    @abc.abstractmethod
    def xy2latlon(self):
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

    staticdata = None

    def __init__(self, core, polygon_geog=None, tilingsystem=None):

        self.core = core

        if polygon_geog is None:
            polygon_geog = GlobalTile(self.core.projection).polygon()
        self.polygon_geog = polygon_geog
        self.polygon_proj= transform_geometry(polygon_geog, self.core.projection)
        self.bbox_geog = get_geom_boundaries(self.polygon_geog, rounding=self.core.res/1000000.0)
        self.bbox_proj = get_geom_boundaries(self.polygon_proj, rounding=self.core.res)

        if tilingsystem is None:
            tilingsystem = GlobalTile(self.core.projection)
        self.tilingsystem = tilingsystem

    def __getattr__(self, item):
        '''
        short link for items of core
        '''
        if item in self.core.__dict__:
            return self.core.__dict__[item]
        else:
            return self.__dict__[item]

    @abc.abstractmethod
    def get_polygon(self, staticdata=staticdata):
        return None

    @abc.abstractmethod
    def get_bbox_geog(self):
        bbox = self.polygon_geog.GetEnvelope()
        return bbox

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
    def __init__(self, core, polygon_proj, x0, y0):

        self.core = core
        self.x0 = x0
        self.y0 = y0
        self.xstep = self.core.tile_xsize_m
        self.ystep = self.core.tile_ysize_m
        self.polygon_proj = polygon_proj
        self.bbox_proj = get_geom_boundaries(self.polygon_proj, rounding=self.core.res)


    def __getattr__(self, item):
        '''
        short link for items of core
        '''
        if item in self.core.__dict__:
            return self.core.__dict__[item]
        else:
            return self.__dict__[item]


    def get_tile(self, x,y):
        return Tile(self.core.projection, 'xybounds')

    @abc.abstractmethod
    def get_tile_name(self, x, y):
        return


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
    def __init__(self, core, name, x0, y0):
        self.name = name
        self.typename = self.get_type_name()
        self.core = core
        self.llx = x0
        self.lly = y0
        self.x_size_m = x0 + self.core.tile_xsize_m
        self.y_size_m = y0 + self.core.tile_ysize_m
        self.x_size_px = self.x_size_m / self.core.res
        self.y_size_px = self.y_size_m / self.core.res

    def __getattr__(self, item):
        '''
        short link for items of core
        '''
        if item in self.core.__dict__:
            return self.core.__dict__[item]
        else:
            return self.__dict__[item]


    @abc.abstractmethod
    def get_type_name(self):
        """
        :returns the tile type name
        """
        return

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
                self.llx + self.x_size_m, self.lly + self.y_size_m)

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
            full tile name e.g. EU075M_E048N012T6

        Returns
        -------
        list
            a list contain the geotransfrom elements

        """
        geot = [self.llx, self.res, 0,
                self.lly + self.y_size_m, 0, -self.res]

        return geot

    def ij_2_xy(self, i, j):
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
            x coordinate in the TilingSystem
        y : number
            y coordinate in the TilingSystem
        """

        gt = self.geotransform()

        x = gt[0] + i * gt[1] + j * gt[2]
        y = gt[3] + i * gt[4] + j * gt[5]

        return x, y

    def get_geotags(self):
        """
        Return geotags for given tile used as geoinformation for GDAL
        """
        geotags = {'geotransform': self.geotransform(),
                   'spatialreference': self.projection}

        return geotags

class GlobalTile(object):

    def __init__(self, projection, name):
        self.core = projection
        self.name = name

    @abc.abstractmethod
    def polygon(self):
        return


def create_wkt_geometry(geometry_wkt, epsg=4326):
    """
    return extent geometry

    Parameters
    ----------
    geometry_wkt : string
        WKT text containing points of geometry (e.g. polygon)
    epsg : int
        EPSG code of spatial reference of the points.

    Return
    ------
    OGRGeomtery
        a geometry representing the extent_m of given sub-grid

    """
    geom = ogr.CreateGeometryFromWkt(geometry_wkt)
    geo_sr = osr.SpatialReference()
    geo_sr.SetWellKnownGeogCS("EPSG:{}".format(str(epsg)))
    geom.AssignSpatialReference(geo_sr)
    return geom

def transform_geometry(geometry, projection):
    """
    return extent geometry

    Parameters
    ----------
    geometry_wkt : string
        WKT text containing points of geometry (e.g. polygon)
    epsg : int
        EPSG code of spatial reference of the points.

    Return
    ------
    OGRGeomtery
        a geometry representing the extent_m of given sub-grid

    """

    out_srs = projection.osr_spref
    geometry_out = geometry.Clone()
    geometry_out.TransformTo(out_srs)
    geometry = None
    return geometry_out


def get_geom_boundaries(geometry, rounding=1.0):
    limits = geometry.GetEnvelope()
    limits = [int(x / rounding) * rounding for x in limits]
    return limits

'''
class TiledProjectedLocation(object):

    #Spatial information of a location in a
    tiled projection system


    def __init__(self,
                 grid=None,
                 u = 0.0,
                 v = 0.0):

        self.grid = grid
        self.tile = grid.get_tile(u, v)
        self.geoscoords = (u, v)
        self.projcoords = self.geog2proj(u, v)
        self.tilecoords = self.geog2tile(u, v)


    def geog2proj(self, u, v):
        x, y = self.grid.latlon2xy
        return x, y

    def proj2geog(self, u, v):
        x, y = self.grid.xy2latlon
        return x, y

    def proj2tile(self, u, v):
        x, y = self.grid.xy2ij
        return x, y
    def tile2proj(self, u, v):
        x, y = self.grid.ij2xy
        return x, y

    def geog2tile(self, u, v):
        a, b = self.geog2proj(u, v)
        x, y = self.proj2tile(a, b)
        return x, y
    def tile2geog(self, u, v):
        a, b = self.tile2proj(u, v)
        x, y = self.proj2geog(a, b)
        return x, y
'''
