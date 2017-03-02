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
import pyproj

import numpy as np

from osgeo import osr
from osgeo import ogr

class Projection():
    """
    Projection class holding and translating the definitions of a projection when initialising.

    Parameters
    ----------
    epsg : integer
        The EPSG-code of the sptaial reference. As from http://www.epsg-registry.org
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

    def __init__(self, res, subgrids):

        self.res = res
        self.tiletype, \
        self.tile_xsize_m, \
        self.tile_ysize_m = self.link_res_2_tilesize(self.res, get_size=True)

        self.subgrids = subgrids
        pass

    def __getattr__(self, item):
        if item in self.subgrids:
            return self.subgrids[item]
        else:
            return self.__dict__[item]

    @abc.abstractmethod
    def link_res_2_tilesize(self):
        pass

    @abc.abstractmethod
    def latlon2xy(self):
        pass

    @abc.abstractmethod
    def xy2latlon(self):
        pass




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

def transform_geometry(geometry, Projection):
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

    out_srs = Projection.osr_spref
    geometry.TransformTo(out_srs)

    return geometry

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

    def __init__(self, Projection, TilingSystem=None):
        self.projection = Projection
        if TilingSystem is None:
            TilingSystem = GlobalTile(Projection)
        self.tilingsystem = TilingSystem


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
    def __init__(self, projection, geog_polygon, res,  x0, y0, xstep, ystep):

        self.projection = Projection
        self.res = res
        self.x0 = x0
        self.y0 = y0
        self.xstep = xstep
        self.ystep = ystep
        self.polygon_geog = geog_polygon
        self.polygon_proj = transform_geometry(geog_polygon, projection)
        self.bbox_geog = self.get_boundaries(self.polygon_geog, rounding=0.001)
        self.bbox_proj = self.get_boundaries(self.polygon_proj, rounding=self.res)

    def get_boundaries(self, geometry, rounding=1):
        limits = self.extent_proj.GetEnvelope()
        limits = [int(x / rounding) * rounding for x in limits]
        return limits

    def get_tile(self, x,y):
        return Tile(self.projection, 'xybounds')

    @abc.abstractmethod
    def get_tile_name(self, x, y):
        return


class Tile(object):

    def __init__(self, Projection, name, limits):
        self.name = name

class GlobalTile(object):

    def __init__(self, Projection, name, limits):
        self.name = name


class Equi7TilingSystem(TilingSystem):

    pass


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

class Equi7Tile(Tile):

    pass