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

Code for the Equi7 Grid.

@author: Bernhard Bauer-Marschallinger, bbm@geo.tuwien.ac.at

'''

import os
import pickle
import copy
import itertools

import numpy as np

from TiledProjection import TiledProjectionSystem
from TiledProjection import TiledProjection
from TiledProjection import TPSProjection
from TiledProjection import TilingSystem
from TiledProjection import Tile
from TiledProjection import create_wkt_geometry

def _load_static_data(module_path):
    # load the data, raise the error if failed to load equi7grid.dat
    equi7_data = None
    fname = os.path.join(os.path.dirname(module_path), "data", "equi7grid.dat")
    with open(fname, "rb") as f:
        equi7_data = pickle.load(f)
    return equi7_data


class Equi7Grid(TiledProjectionSystem):
    """
    Equi7 Grid

    Parameters
    ----------
    res : float
        The tile resolution
    """

    # static attribute
    _static_equi7_data = _load_static_data(__file__)
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

        # check if the equi7grid.data have been loaded successfully
        if Equi7Grid._static_equi7_data is None:
            raise ValueError("cannot load Equi7Grid ancillary data!")
        # check if res is allowed
        if res not in Equi7Grid._static_res:
            raise ValueError("Resolution {}m is not supported!".format(res))

        # initializing
        super(Equi7Grid, self).__init__(res, nametag='Equi7')
        self.core.projection = 'multiple'

    def define_subgrids(self):
        subgrids = dict()
        for sg in self._static_subgrid_ids:
            subgrids[sg] = Equi7Subgrid(self.core, sg)
        return subgrids

    @staticmethod
    def link_res_2_tilesize(res, get_size=False):
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

        if get_size == True:
            result = (tile_code, tile_size_m, tile_size_m)
        else:
            result = tile_code

        return result


    def find_overlapping_tiles(self, tilename, target_tiletype):
        """
        find the family tiles which share the same extent_m but
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

        # found family tiles

        # allow also long-form of tilename
        if len(tilename) == 17:
            tilename = tilename[7:]

        family_tiles = list()

        if target_tiletype >= self.core.tiletype:
            t_span = int(target_tiletype[-1]) * 100000
            t_east = (self.llx / t_span) * t_span / 100000
            t_north = (self.lly / t_span) * t_span / 100000
            name = "E{:03d}N{:03d}{}".format(t_east, t_north, target_tiletype)
            family_tiles.append(name)
        else:
            sub_span = int(target_tiletype[-1]) * 100000
            n = int(self.core.tile_ysize_m / sub_span)
            for x, y in itertools.product(range(n), range(n)):
                s_east = (self.llx + x * sub_span) / 100000
                s_north = (self.lly + y * sub_span) / 100000
                name = "E{:03d}N{:03d}{}".format(s_east, s_north, target_tiletype)
                family_tiles.append(name)
        return family_tiles


class Equi7Subgrid(TiledProjection):

    def __init__(self, core, continent):

        data = Equi7Grid._static_equi7_data[continent]
        _core = copy.copy(core)
        _core.tag = continent
        _core.projection = TPSProjection(wkt=data['project'])

        self.core = _core
        self.polygon_geog = create_wkt_geometry(data['extent'])
        self.tilesys = Equi7TilingSystem(self.core, self.polygon_geog)

        super(Equi7Subgrid, self).__init__(self.core, self.polygon_geog, self.tilesys)


    def get_polygon(self):
        pass

class Equi7TilingSystem(TilingSystem):
    """
    Equi7 tiling system class, providing methods for queries and handling.

    A tile in the Equi7 core system.
    """



    def __init__(self, core, polygon):

        super(Equi7TilingSystem, self).__init__(core, polygon, 0, 0)

        self.msg1 = '"tilename" is not properly defined! Examples: ' \
                    '"{0}{1:03d}M_E012N036{2}" ' \
                    'or "E012N036{2}"'.format(self.core.tag, self.core.res, self.core.tiletype)
        self.msg2 = 'East and North coordinates of lower-left-pixel ' \
                    'must be multiples of {}00km!'.format(self.core.tile_ysize_m / 100000)

    def create_tile(self, name=None, x=None, y=None):

        if x is not None and y is not None and name is None:
            llx, lly = self.round_xy2lowerleft(x, y)
        elif name is not None and x is None and y is None:
            llx, lly = self.tilename2lowerleft(name)
        else:
            raise AttributeError('"name" or "x"&"y" must be defined!')

        # get name of tile (assures long-form of tilename, even if short-form is given)
        name = self.encode_tilename(llx, lly)
        # set True if land in the tile
        covers_land = self.check_land_coverage(tilename=name[7:])

        return Equi7Tile(self.core, name, llx, lly, covers_land=covers_land)


    def point2tilename(self, x0, y0):
        llx, lly = self.round_xy2lowerleft(x0, y0)
        return self.encode_tilename(llx, lly)


    def tilename2lowerleft(self, name):
        _, _, _, llx, lly, _ = self.decode_tilename(name)
        return llx, lly


    def encode_tilename(self, llx, lly):

        # gives long-form of tilename (e.g. "EU500M_E012N018T6")
        tilename = "{}{:03d}M_E{:03d}N{:03d}{}".format(self.core.tag, self.core.res,
                                                        llx / 100000, lly / 100000,
                                                        self.core.tiletype)
        return tilename


    def check_tilename(self, name):

        check = False
        self.decode_tilename(name)
        check = True
        return check

    def decode_tilename(self, name):

        tf = self.core.tile_ysize_m / 100000

        # allow short-form of tilename (e.g. "E012N018T6")
        if len(name) == 10:
            tile_size_m = int(name[-1]) * 100000
            if tile_size_m != self.core.tile_xsize_m:
                raise ValueError(self.msg1)
            llx = int(name[1:4])
            if llx % tf:
                raise ValueError(self.msg2)
            lly = int(name[5:8])
            if lly % tf:
                raise ValueError(self.msg2)
            tile_code = name[-2:]
            if tile_code != self.core.tiletype:
                raise ValueError(self.msg1)
            subgrid_id = self.core.tag
            res = self.core.res

        # allow long-form of tilename (e.g. "EU500M_E012N018T6")
        elif len(name) == 17:
            subgrid_id = name[0:2]
            if subgrid_id != self.core.tag:
                raise ValueError(self.msg1)
            res = int(name[2:5])
            if res != self.core.res:
                raise ValueError(self.msg1)
            tile_size_m = int(name[-1]) * 100000
            if tile_size_m != self.core.tile_xsize_m:
                raise ValueError(self.msg1)
            llx = int(name[8:11])
            if llx % tf:
                raise ValueError(self.msg2)
            lly = int(name[12:15])
            if lly % tf:
                raise ValueError(self.msg2)
            tile_code = name[-2:]
            if tile_code != self.core.tiletype:
                raise ValueError(self.msg1)

        # wrong length
        else:
            raise ValueError(self.msg1)

        return subgrid_id, res, tile_size_m, llx*100000, lly*100000, tile_code


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

        tsize = self.core.tile_xsize_m
        factor = tsize / 100000

        llxs = range(xmin / tsize * factor, xmax / tsize * factor + 1, factor)
        llys = range(ymin / tsize * factor, ymax / tsize * factor + 1, factor)
        tx, ty = np.meshgrid(llxs, llys)
        tx = tx.flatten()
        ty = ty.flatten()

        tilenames = list()
        for i, _ in enumerate(tx):
            tilenames.append(self.encode_tilename(tx[i]*100000, ty[i]*100000))
        return tilenames


    def check_land_coverage(self, tilename=None, all_tiles=False):
        """
        check if a tile covers land
        """

        land_tiles = Equi7Grid._static_equi7_data[self.core.tag]["coverland"][self.core.tiletype]
        if all_tiles:
            return list(land_tiles)
        if tilename is not None:
            if self.check_tilename(tilename):
                return tilename in land_tiles


class Equi7Tile(Tile):
    """
    Equi7 Tile class

    A tile in the Equi7 grid system.
    """

    def __init__(self, core, name, x0, y0, covers_land):
        super(Equi7Tile, self).__init__(core, name, x0, y0)
        self.covers_land = covers_land

    @property
    def shortname(self):
        return self.name[7:]

