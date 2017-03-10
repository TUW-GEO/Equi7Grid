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
import fractions
import math
from geometry import *


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
        point_geom = create_point_geom(lon, lat, lonlatprojection)

        # search for co-locating subgrid
        subgrid = self.locate_geometry_in_subgrids(point_geom)[0]

        x, y, = uv2xy(lonlatprojection.osr_spref,
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

        x, y, = uv2xy(lonlatprojection.osr_spref,
                      self.subgrids[subgrid].core.projection.osr_spref,
                      lon,
                      lat)

        return np.full_like(x, subgrid, dtype=(np.str, len(subgrid))), x, y


    @abc.abstractmethod
    def link_res_2_tilesize(self):
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
        self.polygon_proj= transform_geometry(polygon_geog, self.core.projection)
        self.bbox_geog = get_geom_boundaries(self.polygon_geog, rounding=self.core.res/1000000.0)
        self.bbox_proj = get_geom_boundaries(self.polygon_proj, rounding=self.core.res)

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

        return uv2xy(self.core.projection.osr_spref, lonlatprojection.osr_spref, x, y)

def sigma0_resample_operations(array):

    ind_exclude = (array <= -2000) | (array >= -500)
    array[ind_exclude] = -9999
    array_en = array * 0.01
    array_en[array == -9999] = -9999

    return array

def downsample_via_pixel_indices(array, res_f, res_c, bbox, function=sigma0_resample_operations):
    """ downsample with pixel averaging and consecutive filtering.

    This method perform the masking, averaging, and resampling
    (and consecutive gaussian filtering)

    """

    fine_id, n_pixels_x, n_pixels_y = translate_indices(res_f, res_c, bbox, get_px_counts=True)
    #cum_pixels = np.cumsum(n_pixels).astype(np.uint16)

    ratio = 1.0 * res_c / res_f

    # retrieve the indexes of the course pixels
    if ratio.is_integer():
        course_id = (fine_id[0::int(ratio), 0::int(ratio)]).copy()
    else:
        # if ratio is float: np.unique is too slow for T1->T6
        # therefore changes in values of fine_idx are detected by subtraction fine_idx with fine_idx shifted by
        # one index. With this resulting mask course_idx can be gained.
        temp_idx = np.append(fine_id[1:, 0] - fine_id[:-1, 0], 1)
        temp_idy = np.append(fine_id[0, 1:] - fine_id[0, :-1], 1)
        temp_id_2d = np.outer(temp_idx, temp_idy)
        course_id = fine_id[temp_id_2d.astype(bool)]
        course_id = course_id.reshape(int(len(temp_idx)/ratio), int(len(temp_idy)/ratio))
        temp_idx = 0

    # calculate the  number of pixels (in A tile) overlapping
    # within each pixel of C tile
    pix_num = np.outer(n_pixels_y, n_pixels_x)
    # minimum number of valid fine pixels in coarse pixels
    limit = (pix_num / 100.0 * 99.0).astype(np.uint16)

    if function is not None:
        array = function(array)

    pass



def fast_mask_counting(mask, course_shape, pattern):
    """
    counting the number of masked fine pixels in individual coarse pixels.
    uses jit from numba for speeding up the loops.

    Parameters
    ----------
    mask : numpy array
        logical array: "1" at index of valid fine pixels,
                       "0" at index of non-valid fine pixels
    course_shape : tuple
        (dim_x, dim_y) shape of the target coarse pixel array
    pattern : numpy array
        array, telling for each coarse-scale position the highest fine pixel index

    Returns
    -------
    result: numpy array
        array holding the counts of non-valid (masked) fine pixels per coarse pixel
    """

    # input check
    error_msg = 'fast_mask_counting: input is not valid!'
    types = [np.ndarray, np.array, np.memmap, np.flatiter]
    if (type(mask) not in types) or (type(pattern) not in types) or \
       (type(course_shape) != tuple) or (mask.dtype != 'int8') or (pattern.dtype != 'uint16'):
        raise TypeError(error_msg)

    # fast computation using just in time (jit)
    @jit()
    def run_fast_mask_counting(mask, course_shape, pattern):
        m, n = mask.shape
        result = np.zeros(course_shape, dtype=np.int32)
        x = 0
        for i in range(m):
            if i == pattern[x]:
                x += 1
            y = 0
            for j in range(n):
                if j == pattern[y]:
                    y += 1
                if mask[i, j] != 1:
                    result[x, y] += 1
        return result

    result = run_fast_mask_counting(mask, course_shape, pattern)
    return result



def translate_indices(res_f, res_c, bbox, get_px_counts=False):

    # resolution of grid
    res_f = res_f
    # target resolution
    res_c = res_c
    if res_f >= res_c:
        raise ValueError('"res_c (={}m) must be larger than '
                         'source grid resolution (={}m)'.format(res_c, res_f))

    xsize_m = bbox[2] - bbox[0]
    ysize_m = bbox[3] - bbox[1]
    xsize_f = xsize_m / float(res_f)
    ysize_f = ysize_m / float(res_f)
    xsize_c = xsize_m / float(res_c)
    ysize_c = ysize_m / float(res_c)

    pattern_f = calc_pixel_index_pattern(res_f, res_c)

    pattern_length_f = sum(pattern_f)
    pattern_length_c = len(pattern_f)

    if (xsize_m % (res_c*pattern_length_c) != 0) or \
       (ysize_m % (res_c*pattern_length_c) != 0):
        raise ValueError('"bbox" must have width and height '
                         'dividable by {}!'.format(res_c*pattern_length_c))

    # create template
    pattern_tmpl = list()
    for i in range(len(pattern_f)):
        pattern_tmpl.extend([i] * pattern_f[i])
    # kx, ky:  number of patterns that bbox spans in x and y direction
    kx = int(xsize_f) / pattern_length_f
    ky = int(ysize_f) / pattern_length_f
    idx = np.tile(pattern_tmpl, kx)
    idy = np.tile(pattern_tmpl, ky)
    idx += np.repeat(np.arange(kx) * len(pattern_f), pattern_length_f)
    idy += np.repeat(np.arange(ky) * len(pattern_f), pattern_length_f)

    # create index array
    # TODO: transpose the index? flip the array upside-down?
    index = np.zeros((ysize_f, xsize_f), dtype=np.uint32)
    for i, v in enumerate(idy):
        index[i, :] = idx + v * xsize_c

    # gets a vector holding the fine pixel per coarse pixel
    if get_px_counts:
        n_pixels_x = (np.unique(idx, return_counts=True)[1]).astype(np.uint16)
        n_pixels_y = (np.unique(idy, return_counts=True)[1]).astype(np.uint16)
        return index, n_pixels_x, n_pixels_y
    else:
        return index


def calc_pixel_index_pattern(res1, res2):
    '''
    finds best  representation of overlapping pixels lattices through optimal rounding.
    is periodic and aims for most symmetry-

    :param res1:
    :param res2:
    :return:
    '''


    # get fraction of the two resolutions
    frac = fractions.Fraction('{}/{}'.format(res1, res2))
    numerator = frac._numerator
    denominator = frac._denominator

    # pattern showing digitized relation of pixel sizes.
    pattern = list()

    # algorithm finding the integers of the pixel lattice relation
    n_items = 0
    excess_sum = 0.0
    remainder = denominator
    for i in range(numerator):
        r = remainder / float(numerator - n_items)
        item = math.ceil(r)
        excess = item - r
        excess_sum += excess
        if excess_sum >= 0.5:
            item -= 1
            excess_sum -= 1.0
        n_items += 1
        remainder -= item
        pattern.append(int(item))

    return pattern

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
        self.bbox_proj = get_geom_boundaries(self.polygon_proj, rounding=self.core.res)


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
    def encode_tilename(self, x, y):
        return

    @abc.abstractmethod
    def decode_tilename(self, name):
        a = None
        return a

    @abc.abstractmethod
    def identify_tiles_per_bbox(self, bbox):
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


    def create_tiles_per_bbox(self, bbox):
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
        tilenames = self.identify_tiles_per_bbox(bbox)
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


    def create_daskarray_per_bbox(self, bbox):

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
                self.llx + self.core.tile_xsize_m, self.lly + self.core.tile_xsize_m)

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

    __metaclass__ = abc.ABCMeta

    def __init__(self, projection, name):
        self.projection = projection
        self.name = name

    @abc.abstractmethod
    def polygon(self):
        return


