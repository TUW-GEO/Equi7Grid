# Copyright (c) 2022, TU Wien, Department of Geodesy and Geoinformation (GEO).
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
"""
Code for the Equi7 Grid.
"""

import os
import pickle
import copy
import itertools

import numpy as np

from pytileproj.base import TiledProjectionSystem
from pytileproj.base import TiledProjection
from pytileproj.base import TPSProjection
from pytileproj.base import TilingSystem
from pytileproj.base import Tile
from pytileproj.geometry import create_geometry_from_wkt
from geographiclib.geodesic import Geodesic


def _load_static_data(module_path):
    """
    load the data, raise the error if failed to load equi7grid.dat

    Parameters
    ----------
    module_path : string
        mainpath of the Equi7Grid module

    Returns
    -------
    equi7data : dict
        dictionary containing for each subgrid...
            a) the multipolygon 'zone_extent'
            b) the WKT-string 'projection'
            c) the sets for T6/T3/T1-tiles covering land 'coverland'
            d) the Equi7Grid version 'version'

    """
    equi7_data = None
    fname = os.path.join(os.path.dirname(module_path), "data", "equi7grid.dat")
    with open(fname, "rb") as f:
        equi7_data = pickle.load(f)
    return equi7_data


class Equi7Grid(TiledProjectionSystem):
    """
    Equi7Grid class object, inheriting TiledProjectionSystem() from pytileproj.

    Attributes
    ----------
    _static_data  : dict
        dictionary containing for each subgrid...
            a) the multipolygon 'zone_extent'
            b) the WKT-string 'projection'
            c) the sets for T6/T3/T1-tiles covering land 'coverland'
            d) the Equi7Grid version 'version'
    _static_subgrid_ids : list of strings
        lists the acronyms of the 7 (continental) subgrids
    _static_tilecodes : list of strings
        lists the 3 tile acronyms
    _static_sampling : list of int
        lists all allowed grid samplings
    """

    # static attribute
    _static_data = _load_static_data(__file__)
    # sub grid IDs
    _static_subgrid_ids = list(_static_data.keys())
    # supported tile widths (linked to the grid sampling)
    _static_tilecodes = ["T6", "T3", "T1"]
    # supported grid spacing ( = the pixel sampling)
    _static_sampling = [
        6000, 3000, 1000, 800, 750, 600, 500, 400, 300, 250, 200, 150, 125, 100, 96, 80,
        75, 64, 60, 50, 48, 40, 32, 30, 25, 24, 20, 16, 10, 8, 5, 4, 2, 1
    ]

    def __init__(self, sampling, tile_names_in_m=False):
        """
        Initialises an Equi7Grid class for a specified sampling.

        Parameters
        ----------
        sampling : int
            the grid sampling = size of pixels; in metres.
        tile_names_in_m : bool, optional
            controls whether the tile names are in metres or in km when > 1000m

        """

        self._tile_names_in_m = tile_names_in_m

        # check if the equi7grid.data have been loaded successfully
        if Equi7Grid._static_data is None:
            raise ValueError("cannot load Equi7Grid ancillary data!")
        # check if sampling is allowed
        if sampling not in Equi7Grid._static_sampling:
            raise ValueError("Sampling {}m is not supported!".format(sampling))

        # initializing
        super(Equi7Grid, self).__init__(sampling, tag='Equi7')
        self.core.projection = 'multiple'
        self.core.tile_names_in_m = tile_names_in_m

    @staticmethod
    def encode_sampling(sampling, tile_names_in_m=False):
        """
        provides a string representing the sampling (e.g. for the tilenames)

        Parameters
        ----------
        sampling : int
            the grid sampling = size of pixels; in metres.

        Returns
        -------
        sampling_str : str
            string representing the sampling
        tile_names_in_m : bool, optional
            controls whether the tile names are in metres or in km when > 1000m
        """
        if tile_names_in_m:
            sampling_str = str(sampling).rjust(3, '0')
        else:
            if sampling <= 999:
                sampling_str = str(sampling).rjust(3, '0')
            if sampling >= 1000:
                sampling_str = "".join((str(sampling / 1000.0)[0], 'K', str(sampling / 1000.0)[2]))
        return sampling_str

    @staticmethod
    def decode_sampling(sampling_str, tile_names_in_m=False):
        """
        converts the string representing the sampling (e.g. from the tilenames)
        to an integer value in metres

        Parameters
        ----------
        sampling_str : str
            string representing the sampling
        tile_names_in_m : bool, optional
            controls whether the tile names are in metres or in km when > 1000m

        Returns
        -------
        sampling : int
            the grid sampling = size of pixels; in metres.
        """
        if tile_names_in_m:
            sampling = int(sampling_str)
        else:
            if len(sampling_str) != 3:
                raise ValueError('Resolution is badly defined!')
            if sampling_str[1] == 'K':
                sampling = int(sampling_str[0]) * 1000 + int(sampling_str[2]) * 100
            else:
                sampling = int(sampling_str)
        return sampling

    def define_subgrids(self):
        """
        Builds the grid's subgrids from a static file.

        Returns
        -------
        subgrids : dict of Equi7Subgrid
            dict of all subgrids of the grid
        """
        # Set the tile names in metres or in km
        # Reinforce copying into self.core, as it's initialized in TiledProjectionSystem not in Equi7Grid so can't
        # take this custom property at init time. This is a bit of a hack but its needed in Equi7Subgrid()
        self.core.tile_names_in_m = self._tile_names_in_m

        subgrids = dict()
        for sg in self._static_subgrid_ids:
            subgrids[sg] = Equi7Subgrid(self.core, sg)
        return subgrids

    def get_tiletype(self, sampling=None):
        """
        Returns the tilecode defined for the grid's sampling

        Parameters
        ----------
        sampling : int, optional
            the grid sampling = size of pixels; in metres.

        Returns
        -------
        tilecode : str
            tilecode (related the tile size of the grid)
        """

        # get the tile code of the grid instance
        if sampling is None:
            return self._get_tiletype()

        sampling = int(sampling)

        # allowing sampling of [6000, 3000, 1000, 800, 750, 600, 500, 400, 300, 250, 200,
        # 150, 125, 100, 96, 80, 75, 64] metres
        if ((sampling in range(64, 6001)) and (600000 % sampling == 0)):
            tilecode = "T6"
        # allowing sampling of [60, 50, 48, 40, 32, 30, 25, 24, 20] metres
        elif ((sampling in range(20, 61)) and (300000 % sampling == 0)):
            tilecode = "T3"
        # allowing sampling of [16, 10, 8, 5, 4, 2, 1] metres
        elif ((sampling in range(1, 17)) and (100000 % sampling == 0)):
            tilecode = "T1"
        else:
            msg = "Error: Given resolution %d is not supported!" % sampling
            msg += " Supported resolutions: {}".format(str(Equi7Grid._static_sampling))
            raise ValueError(msg)

        return tilecode

    def get_tilesize(self, sampling):
        """
        Return the tile size in metres defined for the grid's sampling

        Parameters
        ----------
        sampling : int
            the grid sampling = size of pixels; in metres.

        Returns
        -------
        xsize, ysize : int
            tile size in x and y direction defined for the grid's sampling

        """
        xsize = {
            'T6': 600000,
            'T3': 300000,
            'T1': 100000
        }[self.get_tiletype(sampling)]
        ysize = {
            'T6': 600000,
            'T3': 300000,
            'T1': 100000
        }[self.get_tiletype(sampling)]
        return xsize, ysize

    def create_tile(self, name):
        """
        shortcut to create_tile, returning a Equi7Tile object

        Parameters
        ----------
        name : str
            name of the tile; e.g EU500M_E012N018T6

        Returns
        -------
        Equi7Tile
            object containing info of the specified tile

        """
        return self.subgrids[name[0:2]].tilesys.create_tile(name)

    def lonlat2ij_in_tile(self, lon, lat, lowerleft=False):
        """
        finds the tile and the pixel indices of a given point in lon-lat-space.
        pixel indices comprise the column and row number (i, j)

        columns go from left to right (easting)
        rows go either
            top to bottom (lowerleft=False)
            bottom to top (lowerleft=True)

        Parameters
        ----------
        lon : number
            longitude coordinate
        lat : number
            latitude coordinate)
        lowerleft : bool, optional
            should the row numbering start at the bottom?
            If yes, it returns lowerleft indices.

        Returns
        -------
        tilename : str
            long form of the tilename containing the lon-lat position
        i : integer
            pixel column number; starts with 0
        j : integer
            pixel row number; starts with 0

        """
        # get the xy-coordinates
        subgrid, x, y = self._lonlat2xy(lon, lat)

        tilename, i, j = self.subgrids[str(subgrid)].tilesys.xy2ij_in_tile(
            x, y, lowerleft=lowerleft)

        return tilename, i, j

    def calc_length_distortion_on_ellipsoid(self, lon, lat):
        """
        function returing the local maximum length distortion k,
        which equals the local areal distortion (as always h=1 for the Azimuthal Equidistant projection)

        Parameters
        ----------
        lon : number
            longitude coordinate
        lat : number
            latitude coordinate)

        Returns
        -------
        k : float
            value of local max lenght distortion = local areal distortion

        """

        # get the subgrid
        sg, _, _ = self.lonlat2xy(lon, lat)

        lon0 = self.subgrids[str(sg)].core.projection.osr_spref.GetProjParm('central_meridian')
        lat0 = self.subgrids[str(sg)].core.projection.osr_spref.GetProjParm('latitude_of_origin')

        # get spherical distance and azimuth between projection centre and point of interest
        geod = Geodesic.WGS84
        gi = geod.Inverse(lat0, lon0, lat, lon)
        c1 = gi['s12']
        az1 = gi['azi1']

        # apply equation for distortion in direction perpendicular to the radius, k:
        # k = c/geod.a / np.sin(c/geod.a)
        k = c1 / geod.a / np.sin(c1 / geod.a)

        return k


class Equi7Subgrid(TiledProjection):
    """
    Equi7Subgrid class object, inheriting TiledProjection() from pytileproj.

    """

    def __init__(self, core, continent):
        """
        Initialises an Equi7Subgrid class for a specified continent.

        Parameters
        ----------
        core : TPSCoreProperty
            defines core parameters of the (sub-) grid
        continent : str
            acronym of the continent, e.g. 'EU' or 'SA'.
        """

        # load WKT string and extent shape
        data = Equi7Grid._static_data[continent]

        _core = copy.copy(core)
        _core.tag = continent
        _core.projection = TPSProjection(wkt=data['wkt'])

        # holds core parameters of the (sub-) grid
        self.core = _core

        # holds name of the subgrid
        self.name = ''.join(('EQUI7_', continent, Equi7Grid.encode_sampling(core.sampling, core.tile_names_in_m), 'M'))

        # holds the extent of the subgrid in the lonlat-space
        self.polygon_geog = create_geometry_from_wkt(data['zone_extent'], epsg=4326)

        # defines the tilingsystem of the subgrid
        self.tilesys = Equi7TilingSystem(self.core, self.polygon_geog)

        super(Equi7Subgrid, self).__init__(self.core, self.polygon_geog, self.tilesys)

    def calc_length_distortion(self, x, y):
        """
        function returing the local maximum length distortion k,
        which equals the local areal distortion (as always h=1 for the Azimuthal Equidistant projection)

        uses the planar coordinates and is much faster, and allows multiple input values

        Parameters
        ----------
        x : number or list of numbers
            projected x coordinate(s) in metres
        y : number or list of numbers
            projected y coordinate(s) in metres


        Returns
        -------
        k : number or np.array
            value of local max lenght distortion = local areal distortion

        """

        # get the major axis of the used Earth ellipsoid
        ellaxis = Geodesic.WGS84.a

        # get the centre of the subgrid's projection
        fe = self.core.projection.osr_spref.GetProjParm('false_easting')
        fn = self.core.projection.osr_spref.GetProjParm('false_northing')

        # create the distances to the projection centre
        dists = np.sqrt((np.array(x) - fe)**2 + (np.array(y) - fn)**2)

        # apply equation for distortion in direction perpendicular to the radius, k:
        # k = c/geod.a / np.sin(c/geod.a)
        #
        # is it just about the distance to the centre (c), and as are equally long
        # on the ellipsoid and on the projected plane (the core of of AEQD!)
        k = dists / ellaxis / np.sin(dists / ellaxis)

        return k


class Equi7TilingSystem(TilingSystem):
    """
    Equi7TilingSystem class, inheriting TilingSystem() from pytileproj.
    provides methods for queries and handling.
    """

    def __init__(self, core, polygon_geog):
        """
        Initialises an Equi7TilingSystem class for a specified continent.

        Parameters
        ----------
        core : TPSCoreProperty
            defines core parameters of the (sub-) grid
        polygon_geog : OGRGeometry
            geometry defining the extent/outline of the subgrid
        """

        super(Equi7TilingSystem, self).__init__(core, polygon_geog, 0, 0)

        self.msg1 = '"tilename" is not properly defined! Examples: ' \
                    '"{0}{1:03d}M_E012N036{2}" or "E012N036{2}"'.format(
                        self.core.tag, self.core.sampling, self.core.tiletype)
        self.msg2 = 'East and North coordinates of lower-left-pixel ' \
                    'must be multiples of {}00km!'.format(
                        self.core.tile_ysize_m // 100000)
        self.msg3 = 'Tilecode must be one of T6, T3, T1!'

    def create_tile(self, name=None, x=None, y=None):
        """
        Returns a Equi7Tile object

        Parameters
        ----------
        name : str
            name of the tile; e.g EU500M_E012N018T6 or E012N018T6
        x : int
            x (right) coordinate of a pixel located in the desired tile
            must to given together with y
        y : int
            y (up) coordinate of a pixel located in the desired tile
            must to given together with x

        Returns
        -------
        Equi7Tile
            object containing info of the specified tile

        Notes
        -----
        either name, or x and y, must be given.
        """

        # use the x and y coordinates for specifing the tile
        if x is not None and y is not None and name is None:
            llx, lly = self.round_xy2lowerleft(x, y)
        # use the tile name for specifing the tile
        elif name is not None and x is None and y is None:
            llx, lly = self.tilename2lowerleft(name)
        else:
            raise AttributeError('"name" or "x"&"y" must be defined!')

        # get name of tile (assures long-form of tilename, even if short-form
        # is given)
        name = self._encode_tilename(llx, lly)
        # set True if land in the tile
        covers_land = self.check_tile_covers_land(tilename=name)

        return Equi7Tile(self.core, name, llx, lly, covers_land=covers_land)

    def point2tilename(self, x, y, shortform=False):
        """
        Returns the name string of an Equi7Tile in which the point,
        defined by x and y coordinates (in metres), is located.

        Parameters
        ----------
        x : int
            x (right) coordinate in the desired tile
            must to given together with y
        y : int
            y (up) coordinate in the desired tile
            must to given together with x
        shortform : Boolean
            option for giving back the shortform of the tile name
            e.g. 'E012N018T6'.

        Returns
        -------
        str
            the tilename in longform e.g. 'EU500M_E012N018T6'
            or in shortform e.g. 'E012N018T6'.

        """
        llx, lly = self.round_xy2lowerleft(x, y)
        return self._encode_tilename(llx, lly, shortform=shortform)

    def encode_tilename(self, llx, lly, sampling, tilecode, shortform=False):
        """
        Encodes a tilename

        Parameters
        ----------
        llx : int
            Lower-left x coordinate.
        lly : int
            Lower-left y coordinate.
        sampling : int
            the grid sampling = size of pixels; in metres.
        tilecode : str
            tilecode
        shortform : boolean, optional
            return shortform of tilename (default: False).

        Returns
        -------
        str
            the tilename in longform e.g. 'EU500M_E012N018T6'
            or in shortform e.g. 'E012N018T6'.
        """

        # gives long-form of tilename (e.g. "EU500M_E012N018T6")
        tilename = "{}{}M_E{:03d}N{:03d}{}".format(
                        self.core.tag,
                        Equi7Grid.encode_sampling(sampling, self.core.tile_names_in_m),
                        int(llx) // 100000,
                        int(lly) // 100000,
                        tilecode)

        if shortform:
            tilename = self.tilename2short(tilename)

        return tilename

    def _encode_tilename(self, llx, lly, shortform=False):
        """
        Encodes a tilename defined by the lower-left coordinates of the tile,
        using inherent information

        Parameters
        ----------
        llx : int
            lower-left x coordinate.
        lly : int
            lower-left y coordinate.
        shortform : boolean, optional
            return shortform of tilename (default: False).

        Returns
        -------
        str
            the tilename in longform e.g. 'EU500M_E012N018T6'
            or in shortform e.g. 'E012N018T6'.

        """
        return self.encode_tilename(llx,
                                    lly,
                                    self.core.sampling,
                                    self.core.tiletype,
                                    shortform=shortform)

    def tilename2short(self, tilename):
        """
        Converts a tilename in longform to shortform
        e.g. 'EU500M_E012N018T6' --> 'E012N018T6'

        Parameters
        ----------
        tilename : str
            longform of the tilename

        Returns
        -------
        str
            shortform of the tilename

        """
        self.check_tilename(tilename)
        if len(tilename) == 17:
            tilename = tilename[7:]
        return tilename

    def tilename2lowerleft(self, tilename):
        """
        Return the lower-left coordinates of the tile

        Parameters
        ----------
        tilename : str
            the tilename in longform e.g. 'EU500M_E012N018T6'
            or in shortform e.g. 'E012N018T6'.

        Returns
        -------
        llx, lly: int
            lower-left coordinates of the tile
        """
        _, _, _, llx, lly, _ = self.decode_tilename(tilename)
        return llx, lly

    def check_tilename(self, tilename):
        """
        checks if the given tilename is valid

        Parameters
        ----------
        tilename : str
            the tilename in longform e.g. 'EU500M_E012N018T6'
            or in shortform e.g. 'E012N018T6'.

        Returns
        -------
        Boolean
            result of the check
        """

        check = False
        self.decode_tilename(tilename)
        check = True
        return check

    def decode_tilename(self, tilename):
        """
        Returns the information assigned to the tilename

        Parameters
        ----------
        tilename : str
            the tilename in longform e.g. 'EU500M_E012N018T6'
            or in shortform e.g. 'E012N018T6'.

        Returns
        -------
        subgrid_id : str
            ID acronym of the subgrid, e.g. 'EU'
        sampling : int
            the grid sampling = size of pixels; in metres.
        tile_size_m : int
            extent/size of the tile; in metres.
        llx : int
            lower-left x coordinate.
        lly : int
            lower-left y coordinate.
        tilecode : str
            tilecode (related the tile size of the grid)
        """
        tf = self.core.tile_ysize_m // 100000

        # allow short-form of tilename (e.g. "E012N018T6")
        if len(tilename) == 10:
            tile_size_m = int(tilename[-1]) * 100000
            if tile_size_m != self.core.tile_xsize_m:
                raise ValueError(self.msg1)
            llx = int(tilename[1:4])
            if llx % tf:
                raise ValueError(self.msg2)
            lly = int(tilename[5:8])
            if lly % tf:
                raise ValueError(self.msg2)
            tilecode = tilename[-2:]
            if tilecode != self.core.tiletype:
                raise ValueError(self.msg1)
            subgrid_id = self.core.tag
            sampling = self.core.sampling

        # allow long-form of tilename (e.g. "EU500M_E012N018T6")
        # tile_names_in_m True or False
        else:
            subgrid_id = tilename[0:2]
            if subgrid_id != self.core.tag:
                raise ValueError(self.msg1)
            tilename_sampling = tilename[2:].split('M')
            sampling = Equi7Grid.decode_sampling(tilename_sampling[0], self.core.tile_names_in_m)
            tilename_remaining = tilename.split('_')[1]
            if sampling != self.core.sampling:
                raise ValueError(self.msg1)
            tile_size_m = int(tilename_remaining[-1]) * 100000
            if tile_size_m != self.core.tile_xsize_m:
                raise ValueError(self.msg1)
            llx = int(tilename_remaining[1:4])
            if llx % tf:
                raise ValueError(self.msg2)
            lly = int(tilename_remaining[5:8])
            if lly % tf:
                raise ValueError(self.msg2)
            tilecode = tilename_remaining[-2:]
            if tilecode != self.core.tiletype:
                raise ValueError(self.msg1)


        return subgrid_id, sampling, tile_size_m, llx * 100000, lly * 100000, tilecode

    def get_congruent_tiles_from_tilename(self,
                                          tilename,
                                          target_sampling=None,
                                          target_tiletype=None):
        """
        finds the "family tiles", which share a congruent or partial overlap,
        but with different sampling and tilecode

        Parameters
        ----------
        tilename : str
            the tilename in longform e.g. 'EU500M_E012N018T6'
            or in shortform e.g. 'E012N018T6'.
        target_sampling : int
            the sampling of the target grid system
        target_tiletype : string
            tilecode string

        Returns
        -------
        list
            list of found tiles
            for smaller tiles: tiles contained in tile with 'tilename'
            for larger tiles: the tile overlap the with 'tilename'

        Notes
        -----
        Either the sampling or tilecode should be given.
        But if both are given, the sampling will be used.
        """

        # return tilenames in shortform or longform?
        if target_sampling is None:
            shortform = True
        else:
            shortform = False

        # check search input
        if target_sampling is not None and target_tiletype is None:
            sampling = target_sampling
        if target_sampling is None and target_tiletype is not None:
            if target_tiletype == 'T1':
                sampling = 10
            elif target_tiletype == 'T3':
                sampling = 20
            elif target_tiletype == 'T6':
                sampling = 500
            else:
                raise ValueError(self.msg3)
        if target_sampling is not None and target_tiletype is not None:
            sampling = target_sampling

        # grid of the searched tiles
        target_grid = Equi7Grid(sampling=sampling)
        target_tiletype = target_grid.core.tiletype
        target_tilesize = target_grid.core.tile_xsize_m

        # features of the input tile(name)
        _, src_sampling, src_tile_size_m, src_llx, src_lly, src_tiletype = \
            self.decode_tilename(tilename)

        # overlapping tiles
        family_tiles = list()

        # for larger tiles
        if target_tiletype >= src_tiletype:
            t_east = (src_llx // target_tilesize) * target_tilesize
            t_north = (src_lly // target_tilesize) * target_tilesize
            name = target_grid.subgrids[self.core.tag].\
                                tilesys.encode_tilename(t_east, t_north,
                                                        sampling,
                                                        target_tiletype,
                                                        shortform=shortform)
            family_tiles.append(name)

        # for smaller tiles
        else:
            n = int(src_tile_size_m // target_tilesize)
            for x, y in itertools.product(range(n), range(n)):
                s_east = (src_llx + x * target_tilesize)
                s_north = (src_lly + y * target_tilesize)
                name = target_grid.subgrids[self.core.tag].\
                                tilesys.encode_tilename(s_east, s_north,
                                                        sampling,
                                                        target_tiletype,
                                                        shortform=shortform)
                family_tiles.append(name)

        return family_tiles

    def check_tile_covers_land(self, tilename=None):
        """
        checks if a tile covers land

        Parameters
        ----------
        tilename : str
            the tilename in longform e.g. 'EU500M_E012N018T6'

        Returns
        -------
        Boolean
        """
        land_tiles = self.list_tiles_covering_land()
        if self.check_tilename(tilename):
            tilename = self.tilename2short(tilename)
            return tilename in land_tiles

    def list_tiles_covering_land(self):
        """
        Returns a list of all tiles in the subgrid covering land

        Returns
        -------
        list
            list containing land tiles
        """

        land_tiles = Equi7Grid._static_data[self.core.tag]["coverland"][
            self.core.tiletype]
        return list(land_tiles)


class Equi7Tile(Tile):
    """
    The Equi7Tile class, inheriting Tile() from pytileproj.

    A tile in the Equi7Grid system, holding characteristics of the tile.
    """

    def __init__(self, core, name, xll, yll, covers_land):
        super(Equi7Tile, self).__init__(core, name, xll, yll)
        self.covers_land = covers_land

    @property
    def shortname(self):
        return self.name[7:]
