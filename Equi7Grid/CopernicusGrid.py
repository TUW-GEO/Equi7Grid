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

from geometry import extent2polygon



class CopernicusGrid(TiledProjectionSystem):
    """
    Equi7 Grid

    Parameters
    ----------
    res : float
        The tile resolution
    """

    # sub grid IDs
    _static_subgrid_ids = ["GLOBAL"]
    # supported tile widths(resolution)
    _static_tilecodes = ["T10", "GLOBAL"]
    # supported grid spacing (resolution)
    _static_res = [1.0/112]
    # set True if one single global tile is needed
    _global_tile = True

    def __init__(self, res):
        """
        construct Copernicus grid system.

        """

        # check if res is allowed
        if res not in CopernicusGrid._static_res:
            raise ValueError("Resolution {}m is not supported!".format(res))

        # initializing
        super(CopernicusGrid, self).__init__(res, nametag='Equi7')
        self.core.projection = 'Latlon'

    def define_subgrids(self):
        subgrids = dict()
        for sg in self._static_subgrid_ids:
            subgrids[sg] = CopernicusSubgrid(self.core, sg)
        return subgrids

    def get_tiletype(self, res, global_tile=_global_tile):
        res = res
        tile_code = None
        if global_tile:
            tile_code = "GLOBAL"
        elif (res == 1.0/112) and (10.0 % res <= 0.000000001):
            tile_code = "T10"
        else:
            msg = "Error: Given resolution %d is not supported!" % res
            msg += " Supported resolutions: {}".format(
                str(CopernicusGrid._static_res))
            raise ValueError(msg)

        return tile_code

    def get_tilesize(self, res):
        xsize = {'T10': 10.0, 'GLOBAL': 360.0}[self.get_tiletype(res)]
        ysize = {'T10': 10.0, 'GLOBAL': 180.0}[self.get_tiletype(res)]
        return xsize, ysize


class CopernicusSubgrid(TiledProjection):

    def __init__(self, core, continent):

        _core = copy.copy(core)
        _core.tag = continent
        _core.projection = TPSProjection(epsg=4326)

        self.core = _core
        self.polygon_geog = extent2polygon([-180,-90,180,90], wkt=self.core.projection.wkt)
        self.tilesys = CopernicusTilingSystem(self.core, self.polygon_geog)

        super(CopernicusSubgrid, self).__init__(self.core, self.polygon_geog, self.tilesys)


    def get_polygon(self):
        pass

class CopernicusTilingSystem(TilingSystem):
    """
    Equi7 tiling system class, providing methods for queries and handling.

    A tile in the Equi7 core system.
    """



    def __init__(self, core, polygon):

        super(CopernicusTilingSystem, self).__init__(core, polygon, 0, 0)


    def create_tile(self, name='GLOBAL', x=None, y=None):
        return CopernicusTile(self.core, name)


    def point2tilename(self, x0, y0):
        return


    def _encode_tilename(self, llx, lly):
        return

    def tilename2lowerleft(self, name):
        return

    def check_tilename(self, name):
        check = False
        self.decode_tilename(name)
        check = True
        return check

    def decode_tilename(self, name):
        pass

    def identify_tiles_overlapping_bbox(self, bbox):
        return



class CopernicusTile(Tile):
    """
    Equi7 Tile class

    A tile in the Equi7 grid system.
    """

    def __init__(self, core, name):
        super(CopernicusTile, self).__init__(core, name, -180, -90)


