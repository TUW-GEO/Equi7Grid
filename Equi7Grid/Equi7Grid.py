# Copyright (c) 2015, Vienna University of Technology (TU Wien), Department of
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
    accurate. Equi7Grid also depends the grids' zone shape files of both aeqd
    and wgs84.

    terminology:
    grid =  grid name (Equi7)
    res= grid resolution/grid spacing (5,10,40,75,500 meters)
    sgrid= subgrid full name (e.g. "AF500M")
    sgrid_id= sub-grid ini (e.g. "AF")
    tile= tile name ("E036N090T6")
    tile_code= tile code representing width of tile (T1:100km, T3:300km,
                T6:600km)
    ftile= full name of tile/sgrid+tile ("AF500M_E036N090T6")


'''

import os
import itertools
import pickle
import numpy as np
from osgeo import ogr
from osgeo import osr
from Equi7Grid import gdalport
from Equi7Grid import parallel


# refactor


def _load_equi7grid_data(module_path):
	# load the data, raise the error if failed to load equi7grid.dat
	equi7_data = None
	fname = os.path.join(os.path.dirname(module_path), "equi7grid.dat")
	with open(fname, "rb") as f:
		equi7_data = pickle.load(f)
	return equi7_data

class Equi7Grid(object):
	"""Equi7 Grid
    """

	# static attribute
	_static_equi7_data = _load_equi7grid_data(__file__)
	# sub grid IDs
	_static_sgrid_ids = ["NA", "EU", "AS", "SA", "AF", "OC", "AN"]
	# supported grid spacing (resolution)
	_static_res = [1000, 800, 750, 600, 500, 400, 300, 250, 200, 150, 125, 100, 96, 80, 75, 64, 60, 50, 48, 40, 32,
				   30, 25, 24, 20, 16, 10, 8, 5, 4, 2, 1]
	# supported tile widths(resolution)
	_static_tilecodes = ["T6", "T3", "T1"]

	def __init__(self, res):
		""" construct Equi7 grid system.

        Parameters
        ----------
        res : float
            The tile resolution
        """
		self.res = int(res)
		# initializing
		self._initialize()

	def _initialize(self):
		""" initialization"""

		# check if the equi7grid.data have been loaded successfully
		if Equi7Grid._static_equi7_data is None:
			self.res = None
			raise ValueError("cannot load Equi7Grid ancillary data!")

		# allowing sampling of [1000, 800, 750, 600, 500, 400, 300, 250, 200, 150, 125, 100, 96, 80, 75, 64] metres
		if ((self.res in range(64, 1000)) and (600000 % self.res == 0)):
			self._tilecode = "T6"
			self._tile_xspan, self._tile_yspan = (600000, 600000)
		# allowing sampling of [60, 50, 48, 40, 32, 30, 25, 24, 20] metres
		elif ((self.res in range(20, 60)) and (300000 % self.res == 0)):
			self._tilecode = "T3"
			self._tile_xspan, self._tile_yspan = (300000, 300000)
		# allowing sampling of [16, 10, 8, 5, 4, 2, 1] metres
		elif ((self.res in range(1, 16)) and (100000 % self.res == 0)):
			self._tilecode = "T1"
			self._tile_xspan, self._tile_yspan = (100000, 100000)
		else:
			invalid_res = self.res
			self.res = None
			msg = "Unsupported resolution {:d}!".format(invalid_res)
			msg += "Supported resolution: {}".format(
				str(Equi7Grid._static_res))
			raise ValueError(msg)

		# keep a reference to the _static_equi7_data
		self._equi7_data = Equi7Grid._static_equi7_data

	@property
	def span(self):
		return self._tile_xspan

	@property
	def tile_xspan(self):
		return self._tile_xspan

	@property
	def tile_yspan(self):
		return self._tile_yspan

	@property
	def sgrids(self):
		return Equi7Grid._static_sgrid_ids

	@property
	def tilecode(self):
		return self._tilecode

	def is_coverland(self, ftile):
		"""check if tile covers land

        Parameters
        ----------
        ftile : string
            full tile name
        """
		return Equi7Tile(ftile).covers_land()

	def get_sgrid_zone_geom(self, sgrid_id):
		"""return sub-grid extent geometry

        Parameters
        ----------
        sgrid_id : string
            sub-grid id string e.g. EU for Europe

        Return
        ------
        OGRGeomtery
            a geometry representing the extent of given sub-grid

        """
		return ogr.CreateGeometryFromWkt(self._equi7_data[sgrid_id]["extent"])

	def get_sgrid_projection(self, sgrid_id):
		"""return sub-grid spatial reference in wkt format

        Parameters
        ----------
        sgrid_id : string
            sub-grid id string e.g. EU for Europe

        Return
        ------
        string
            wkt string representing sub-grid spatial reference

        """
		return self._equi7_data[sgrid_id]["project"]

	def get_sgrid_tiles(self, sgrid_id, tilecode=None):
		"""return all the tiles in the given sub grid"""
		if tilecode is None:
			tilecode = self.tilecode

		return list(self._equi7_data[sgrid_id]["tiles"][tilecode])

	def identify_sgrid(self, geom):
		"""return the overlapped grid ids."""
		sgrid_ids = [sgrid_id for sgrid_id in Equi7Grid._static_sgrid_ids
					 if geom.Intersects(self.get_sgrid_zone_geom(sgrid_id))]
		return sgrid_ids

	def identfy_tile(self, sgrid_id, location):
		"""Return the tile name."""
		east = (int(location[0]) / self._tile_xspan) * self._tile_xspan / 100000
		north = (int(location[1]) / self._tile_yspan) * self._tile_yspan / 100000
		return "{}{:03d}M_E{:03d}N{:03d}{}".format(sgrid_id, self.res,
												   east, north, self.tilecode)

	@staticmethod
	def get_tile_extent(ftile):
		"""return the extent of the tile in the terms of [minX,minY,maxX,maxY]
        """
		return Equi7Tile(ftile).extent

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

	def lonlat2equi7xy(self, lon, lat):
		""" convert (lon, lat) to Equi7 (sgrid_id, x, y)

        Return
        ------
        tuple
            (sgrid_id, x, y) representing the sub-grid id and coordinates
        """
		# create point geometry
		geo_sr = osr.SpatialReference()
		geo_sr.SetWellKnownGeogCS("EPSG:4326")
		geom_pt = ogr.Geometry(ogr.wkbPoint)
		geom_pt.AddPoint(lon, lat)
		geom_pt.AssignSpatialReference(geo_sr)
		# any valid location should be in only one grid
		sgrid_id = self.identify_sgrid(geom_pt)[0]
		# convert to Equi7
		grid_sr = osr.SpatialReference()
		grid_sr.ImportFromWkt(self.get_sgrid_projection(sgrid_id))
		tx = osr.CoordinateTransformation(geo_sr, grid_sr)
		x, y, _ = tx.TransformPoint(lon, lat)

		return (sgrid_id, x, y)

	def equi7xy2lonlat(self, sgrid_id, x, y):
		""" convert Equi7 (sgrid_id, x, y) to (lon, lat)

        Parameters
        ----------
        sgrid_id : string
            the sub-grid id (represent continent) of the coordination,
            should be one of ["NA", "EU", "AS", "SA", "AF", "OC", "AN"]
        x : number
            x coordination
        y : number
            y coordination
        Return
        ------
        tuple
            (lat, lon)
        """
		grid_sr = osr.SpatialReference()
		grid_sr.ImportFromWkt(self.get_sgrid_projection(sgrid_id))
		geo_sr = osr.SpatialReference()
		geo_sr.SetWellKnownGeogCS("EPSG:4326")
		tx = osr.CoordinateTransformation(grid_sr, geo_sr)
		lon, lat, _ = tx.TransformPoint(x, y)
		return (lon, lat)

	@staticmethod
	def decode_tilename(ftile):
		"""keep for compatibility"""
		# example ftile "AF500M_E060N036T6"
		# check the validation of ftile
		sgrid_id = ftile[0:2]
		res = int(ftile[2:5])
		span = int(ftile[-1]) * 100000
		east = int(ftile[8:11]) * 100000
		north = int(ftile[12:15]) * 100000
		tile_code = ftile[-2:]

		if int(res) not in Equi7Grid._static_res:
			raise ValueError("The grid resolution is not supported!")

		equi7_data = Equi7Grid._static_equi7_data
		if ftile[7:] not in equi7_data[sgrid_id]["coverland"][tile_code]:
			msg = "Given tile name was not found in land masses tile list!"
			raise ValueError(msg)

		return (sgrid_id, res, span, east, north, tile_code)

	@staticmethod
	def get_index(dst_ftile, src_ftile):
		"""
        Return the index for oversammpling src ftile to dst ftile.

        Parameters
        ----------
        dst_ftile : string
            dst full tile name e.g. AF075M_E066N030T6
        src_ftile : string
            src full tile name e.g. AF500M_E066N030T6

        Return
        ------
        numpy array
            The index array with the same size as the dst tilename

        """

		dtile = Equi7Tile(dst_ftile)
		stile = Equi7Tile(src_ftile)

		# check if dst tile is a sub tile of src tile
		if dtile.llx < stile.llx or dtile.lly < stile.lly \
				or dtile.llx + dtile.span > stile.llx + stile.span \
				or dtile.lly + dtile.span > stile.lly + stile.span:
			raise ValueError("dst tile should be a sub tile of src tile!")

		index_pattern = {"075T6-500T6": (7, 6, 7),
						 "040T3-500T6": (12, 13),
						 "010T1-500T6": (50,),
						 "040T3-075T6": (2, 2, 2, 1, 2, 2, 2, 2),
						 "010T1-075T6": (7, 8),
						 "010T1-040T3": (4,)
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
		k = dtile.size / pattern_sum + 2
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
		x_idx = x_idx[x_m:x_m + dtile.size]

		# make y index
		yoff = (dtile.lly + dtile.span - stile.lly - stile.span) / -dtile.res
		y_n, y_m = yoff / pattern_sum, yoff % pattern_sum
		y_idx = idx + (y_n * len(pattern))
		# shift idx to the correct start point
		y_idx = y_idx[y_m:y_m + dtile.size]

		# create index array
		index = np.empty((dtile.size, dtile.size), dtype=np.uint32)
		for i, v in enumerate(y_idx):
			index[i, :] = x_idx + v * stile.size
		return index

	@staticmethod
	def get_tile_code(res):
		res = int(res)
		tile_code = None
		if res in [500, 75]:
			tile_code = "T6"
		elif res in [40]:
			tile_code = "T3"
		elif res in [10, 5]:
			tile_code = "T1"
		else:
			msg = "Error: Given resolution %d is not supported!" % res
			raise ValueError(msg)

		return tile_code

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

	def search_tiles(self, geom_area=None, extent=None,
					 sgrid_ids=None, coverland=False):
		"""
        Search the tiles which are intersected by the poly_roi area.

        Parameters
        ----------
        geom : geometry
            a polygon or multipolygon geometery object representing the ROI
        extent : list
            It is a polygon representing the rectangle region of interest
            in the foramt of [lonmin, latmin, lonmax, latmax].

        sgrid_ids : list
            grid ID to specified which continents you want to search. Default
            value is None for searching all continents.

        Returns
        -------
        list
            return a list of  the overlapped tiles' name.
            If not found, return empty list.
        """
		# check input grids
		if sgrid_ids is None:
			sgrid_ids = Equi7Grid._static_sgrid_ids
		elif set(sgrid_ids).issubset(set(Equi7Grid._static_sgrid_ids)):
			sgrid_ids = list(sgrid_ids)
		else:
			raise ValueError("Invalid agrument: grid must one of [ %s ]." %
							 " ".join(Equi7Grid._static_sgrid_ids))

		if not geom_area and not extent:
			print "Error: either geom or extent should be given as the ROI."
			return list()

		# obtain the geometry of ROI
		geo_sr = osr.SpatialReference()
		geo_sr.SetWellKnownGeogCS("EPSG:4326")
		if not geom_area:
			geom_area = gdalport.extent2polygon(extent)
			geom_area.AssignSpatialReference(geo_sr)

		geom_sr = geom_area.GetSpatialReference()
		if geom_sr is None:
			geom_area.AssignSpatialReference(geo_sr)
		elif not geom_sr.IsSame(geo_sr):
			geom_area.TransformTo(geo_sr)

		# intersect the given grid ids and the overlapped ids
		overlapped_grids = self.identify_sgrid(geom_area)
		sgrid_ids = list(set(sgrid_ids) & set(overlapped_grids))

		# finding tiles
		overlapped_tiles = list()
		for sgrid_id in sgrid_ids:
			overlapped_tiles.extend(
				self.__search_sgrid_tiles(geom_area, sgrid_id, coverland))
		return overlapped_tiles

	def __search_sgrid_tiles(self, geom, sgrid_id, coverland):
		"""
        Search the tiles which are overlapping with the given grid.

        Parameters
        ----------
        area_geometry : geometry
            It is a polygon geometry representing the region of interest.
        sgrid_id : string
            sub grid ID to specified which continent you want to search.
            Default value is None for searching all continents.

        Returns
        -------
        list
            Return a list of  the overlapped tiles' name.
            If not found, return empty list.
        """
		# get the intersection of the area of interest and grid zone
		intersect = geom.Intersection(self.get_sgrid_zone_geom(sgrid_id))
		if not intersect:
			return list()
		# The spatial reference need to be set again after intersection
		intersect.AssignSpatialReference(geom.GetSpatialReference())
		# transform the area of interest to the grid coordinate system
		grid_sr = osr.SpatialReference()
		grid_sr.ImportFromWkt(self.get_sgrid_projection(sgrid_id))
		intersect.TransformTo(grid_sr)

		# get envelope of the Geometry and cal the bounding tile of the
		envelope = intersect.GetEnvelope()
		x_min = int(envelope[0]) / self._tile_xspan * self._tile_xspan
		x_max = (int(envelope[1]) / self._tile_xspan + 1) * self._tile_xspan
		y_min = int(envelope[2]) / self._tile_yspan * self._tile_yspan
		y_max = (int(envelope[3]) / self._tile_yspan + 1) * self._tile_yspan

		# make sure x_min and y_min greater or equal 0
		x_min = 0 if x_min < 0 else x_min
		y_min = 0 if y_min < 0 else y_min

		# get overlapped tiles
		overlapped_tiles = list()
		for x, y in itertools.product(xrange(x_min, x_max, self._tile_xspan),
									  xrange(y_min, y_max, self._tile_yspan)):
			geom_tile = gdalport.extent2polygon((x, y, x + self._tile_xspan,
												 y + self._tile_yspan))
			if geom_tile.Intersects(intersect):
				ftile = self.identfy_tile(sgrid_id, [x, y])
				if not coverland or self.is_coverland(ftile):
					overlapped_tiles.append(ftile)

		return overlapped_tiles

	def resample_to_geotiff(self, image, output_dir, gdal_path=None, sgrid_ids=None,
							e7_folder=True, ftiles=None,
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
        sgrid_ids : list
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
				try:
					roi_geom = gdalport.retrieve_raster_boundary(image,
																 gdal_path=gdal_path,
																 nodata=image_nodata)
				except Exception as e:
					print "retrieve_raster_boundary failed:", str(e)
					roi_geom = None
				if roi_geom:
					ftiles = self.search_tiles(geom_area=roi_geom,
											   sgrid_ids=sgrid_ids,
											   coverland=coverland)
				else:
					ds = gdalport.open_image(image)
					img_extent = ds.get_extent()
					ftiles = self.search_tiles(extent=img_extent,
											   sgrid_ids=sgrid_ids,
											   coverland=coverland)
			else:
				ftiles = self.search_tiles(geom_area=roi,
										   sgrid_ids=sgrid_ids,
										   coverland=coverland)
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
					parallel.makedirs(tile_path)
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
			tile_extent = self.get_tile_extent(ftile)
			tile_project = self.get_sgrid_projection(ftile[0:2])

			# prepare options for gdalwarp
			options = {'-t_srs': tile_project, '-of': 'GTiff',
					   '-r': resampling_type,
					   '-te': " ".join(map(str, tile_extent)),
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
		self._sgrid = self.ftile[0:2]
		self.res = int(self.ftile[2:5])
		self.llx = int(self.ftile[8:11]) * 100000
		self.lly = int(self.ftile[12:15]) * 100000
		self._xspan = int(self.ftile[-1]) * 100000
		self._yspan = self._xspan
		self.tilecode = self.ftile[-2:]
		self._xsize = int(self._xspan / self.res)
		self._ysize = int(self._yspan / self.res)

	@property
	def fullname(self):
		return self.ftile

	@property
	def shortname(self):
		return self.ftile[7:]

	@property
	def sgrid(self):
		return self._sgrid

	@property
	def span(self):
		return self._xspan

	@property
	def size(self):
		return self._xsize

	@property
	def shape(self):
		return (self._ysize, self._xsize)

	@property
	def extent(self):
		"""return the extent of the tile in the terms of [minX,minY,maxX,maxY]
        """
		return [self.llx, self.lly,
				self.llx + self._xspan, self.lly + self._yspan]

	def projection(self):
		return Equi7Grid._static_equi7_data[self._sgrid]["project"]

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
				self.lly + self._yspan, 0, -self.res]
		return geot

	def get_tile_geotags(self):
		"""
        Return geotags for given tile used as geoinformation for Geotiff
        """
		geotags = {'geotransform': self.geotransform(),
				   'spatialreference': self.projection()}

		return geotags

	def find_family_tiles(self, res=None, tilecode=None):
		"""find the family tiles which share the same tile extent but
        with different resoluion and tilecode

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
			tilecode = Equi7Grid.get_tile_code(res)
		elif tilecode is not None:
			tilecode = tilecode
		else:
			raise ValueError("either res or tilecode should be given!")

		# found family tiles
		family_tiles = list()

		if tilecode >= self.tilecode:
			t_span = int(tilecode[-1]) * 100000
			t_east = (self.llx / t_span) * t_span / 100000
			t_north = (self.lly / t_span) * t_span / 100000
			name = "E{:03d}N{:03d}{}".format(t_east, t_north, tilecode)
			family_tiles.append(name)
		else:
			sub_span = int(tilecode[-1]) * 100000
			n = int(self.span / sub_span)
			for x, y in itertools.product(range(n), range(n)):
				s_east = (self.llx + x * sub_span) / 100000
				s_north = (self.lly + y * sub_span) / 100000
				name = "E{:03d}N{:03d}{}".format(s_east, s_north, tilecode)
				family_tiles.append(name)
		return family_tiles

	def covers_land(self):
		"""check if tile covers land"""
		land_tiles = Equi7Grid._static_equi7_data[self._sgrid]["coverland"]
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
		if ftile[0:2] not in Equi7Grid._static_sgrid_ids:
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
