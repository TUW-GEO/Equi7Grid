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
import platform
import pickle
import copy

from TiledProjection import TiledProjectionSystem
from TiledProjection import TiledProjection
from TiledProjection import TPSProjection
from TiledProjection import TilingSystem
from TiledProjection import Tile
from TiledProjection import create_wkt_geometry

def _load_static_data(module_path):
    # load the data, raise the error if failed to load equi7grid.dat
    equi7_data = None
    system = platform.system()
    fname = os.path.join(os.path.dirname(module_path), "data", "equi7grid_{0}.dat".format(system))
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
            self.res = None
            raise ValueError("cannot load Equi7Grid ancillary data!")

        # check if res is allowed
        if res not in Equi7Grid._static_res:
            self.res = None
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


    def latlon2xy(self):
        pass


    def xy2latlon(self):
        pass

class Equi7Subgrid(TiledProjection):

    def __init__(self, core, continent):

        data = Equi7Grid._static_equi7_data[continent]
        _core = copy.copy(core)
        _core.tag = continent
        _core.projection = TPSProjection(wkt=data['project'])

        self.core = _core
        self.polygon_geog = create_wkt_geometry(data['extent'])
        self.tilingsystem = Equi7TilingSystem(self.core, self.polygon_geog)

        super(Equi7Subgrid, self).__init__(self.core, self.polygon_geog, self.tilingsystem)


class Equi7TilingSystem(TilingSystem):
    """
    Equi7 tiling system class, providing methods for queries and handling.

    A tile in the Equi7 core system.
    """

    def __init__(self, core, polygon):

        super(Equi7TilingSystem, self).__init__(core, polygon, 0, 0)

    def create_tile(self, x0, y0):

        name = self.create_tile_name(self.core.tag, x0/10000, y0/10000)
        return Equi7Tile(self.core, name, x0/ 10000 * 10000, y0 / 10000 * 10000)

    def create_tile_name(self, continent, x0, y0):

        if (x0 % int(self.core.tiletype[1])) != 0:
            raise ValueError('x0 coordinate does not fit tile type!')
        if (y0 % int(self.core.tiletype[1])) != 0:
            raise ValueError('y0 coordinate does not fit tile type!')

        name = "{}{:03d}M_E{:03d}N{:03d}{}".format(continent, self.core.res,
                                                   x0, y0, self.core.tiletype)
        return name

    def ask_tile_cover_land(self):
        """
        check if a tile covers land
        """
        land_tiles = Equi7Grid._static_equi7_data[Equi7Grid.subgrid]["coverland"]
        return self.shortname in land_tiles[self.tilecode]

class Equi7Tile(Tile):
    """
    Equi7 Tile class

    A tile in the Equi7 grid system.
    """

    def __init__(self, core, name, x0, y0):
        super(Equi7Tile, self).__init__(core, name, x0, y0)

    @property
    def shortname(self):
        return self.name[7:]