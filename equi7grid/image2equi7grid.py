# Copyright (c) 2021, Vienna University of Technology (TU Wien), Department of
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
Created on 2018 July 19

Code for resampling image raster data to the Equi7Grid, yielding tiled raster
images organised in folders for
a) in the continental zones + grid sampling
b) the subgrid tiling

The module should help to easily bring your raster images to the Equi7Grid
spatial reference, using the Equi7TilingSystem() for the file tiling.

@author: Bernhard Bauer-Marschallinger, bbm@geo.tuwien.ac.at
'''

import os
import subprocess
from datetime import datetime
import numpy as np
from scipy import ndimage
from osgeo import gdal
from osgeo import osr
from osgeo import ogr

try:
    from equi7grid.equi7grid import Equi7Grid
except:
    from equi7grid import Equi7Grid

# dict for transfer the datatype and resample type
gdal_datatype = {"uint8": gdal.GDT_Byte,
                 "int16": gdal.GDT_Int16,
                 "int32": gdal.GDT_Int32,
                 "uint16": gdal.GDT_UInt16,
                 "uint32": gdal.GDT_UInt32,
                 "float32": gdal.GDT_Float32,
                 "float64": gdal.GDT_Float64,
                 "complex64": gdal.GDT_CFloat32,
                 "complex128": gdal.GDT_CFloat64
                 }


gdal_resample_type = {"nearst": gdal.GRA_NearestNeighbour,
                      "bilinear": gdal.GRA_Bilinear,
                      "cubic": gdal.GRA_Cubic,
                      "cubicspline": gdal.GRA_CubicSpline,
                      "lanczos": gdal.GRA_Lanczos,
                      "average": gdal.GRA_Average,
                      "mode": gdal.GRA_Mode
                      }


class GdalImage:

    """A sample class to access a image with GDAL library"""

    def __init__(self, gdaldataset, filepath):
        self.dataset = gdaldataset
        self.filepath = filepath

    def close(self):
        """close the dataset"""
        self.dataset = None

    def read_band(self, band_idx, subset=None):
        """Read data from given band.

        Parameters
        ----------
        band_idx : int
            The band index starting from 1.
        subset : list or tuple
            The subset should be in pixels, like this (xmin, ymin, xmax, ymax).

        Returns
        ------
        2darray
            the 2d array including data reading from given band

        """
        if band_idx < 1 or band_idx > self.dataset.RasterCount:
            raise IndexError("band index is out of range")
        band = self.dataset.GetRasterBand(band_idx)
        if subset is None:
            data = band.ReadAsArray(0, 0, band.XSize, band.YSize)
        else:
            data = band.ReadAsArray(subset[0], subset[1], subset[2], subset[3])
        return data

    def XSize(self):
        """get the width of the image"""
        return self.dataset.RasterXSize if self.dataset else None

    def YSize(self):
        """get the height of the image"""
        return self.dataset.RasterYSize if self.dataset else None

    @property
    def shape(self):
        return (self.YSize(), self.XSize())

    def get_band_nodata(self, band_idx=1):
        """get band nodata value.

        Parameters
        ----------
        band_idx : int
            the band index starting from 1

        Returns
        -------
        nodata
            nodata value if it's available, otherwise it will return None

        """
        if band_idx < 1 or band_idx > self.dataset.RasterCount:
            raise IndexError("band index is out of range")
        return self.dataset.GetRasterBand(band_idx).GetNoDataValue()


    def get_raster_nodata(self):
        """get the nodata value for all bands in a list
        this is compatible with the write_image's nodata parameter

        Returns
        -------
        nodata
            no data values in a list if it's available, otherwise it will return None
        """

        nodata = list()
        for i in xrange(0, self.dataset.RasterCount):
            nodata.append(self.dataset.GetRasterBand(i + 1).GetNoDataValue())

        return nodata if len(nodata) >= 0 and not all(d is None for d in nodata) else None

    def read_all_band(self):
        """read the data of all the bands"""
        m = np.full((self.band_count(), self.YSize(), self.XSize()), 0.0)

        for bandIdx in xrange(self.band_count()):
            m[bandIdx] = self.read_band(bandIdx + 1)

        return m

    def get_band_dtype(self, band_idx=1):
        """get the data type of given band"""
        if band_idx < 1 or band_idx > self.dataset.RasterCount:
            raise IndexError("band index is out of range")
        return self.dataset.GetRasterBand(band_idx).DataType

    def geotransform(self):
        """get the geotransform data"""
        return self.dataset.GetGeoTransform() if self.dataset else None

    def projection(self):
        """get the projection string in wkt format"""
        return self.dataset.GetProjection() if self.dataset else None

    def colormap(self, band_idx=1):
        """get the colormap of given band"""
        if band_idx < 1 or band_idx > self.dataset.RasterCount:
            raise IndexError("band index is out of range")
        ct = self.dataset.GetRasterBand(band_idx).GetColorTable()
        if ct is None:
            return None

        colormap = []
        for i in xrange(ct.GetCount()):
            colormap.append(ct.GetColorEntry(i))

        return colormap

    def band_count(self):
        """get the band count"""
        return self.dataset.RasterCount if self.dataset else None

    def get_extent(self):
        """get the extent of the image as (xmin, ymin, xmax, ymax)."""
        geot = self.geotransform()
        return (geot[0], geot[3] + self.YSize() * geot[5],
                geot[0] + self.XSize() * geot[1], geot[3])

    def pixel2coords(self, x, y):
        """Returns global coordinates from pixel x, y coords"""
        xoff, a, b, yoff, d, e = self.geotransform()

        xp = a * x + b * y + xoff
        yp = d * x + e * y + yoff
        return(xp, yp)

    def coords2pixel(self, l1, l2):
        gt = self.geotransform()
        col = int((l1 - gt[0]) / gt[1])
        row = int((l2 - gt[3]) / gt[5])
        if col < 0 or col >= self.XSize() or row < 0 or row >= self.YSize():
            return None
        return [row, col]

    def inside(self, l1, l2):
        """Checks if a pixel is in this image"""
        x, y = self.coords2pixel(l1, l2)

        return x >= 0 and x < self.XSize() and y >= 0 and y < self.YSize()


def image2equi7grid(e7grid, image, output_dir, gdal_path=None, inband=None,
                    subgrid_ids=None, accurate_boundary=True, e7_folder=True,
                    ftiles=None, coverland=True, roi=None, naming_convention=None,
                    compress_type="LZW", resampling_type="bilinear",
                    subfolder=None, overwrite=False, data_type=None,
                    image_nodata=None, tile_nodata=None,
                    scale=None, offset=None, tiled=True, blocksize=512):

    """
    Resample an raster image to tiled images in the Equi7Grid.
    Currently only GTiff is supported as output format.

    Parameters
    ----------
    e7grid : Equi7Grid
        Equi7Grid object defining output projection and sampling.
    image : str
        Image file path.
    output_dir : str
        Output directory path.
    gdal_path : str
        Path to GDAL utilities location.
    inband : str, optional
        Name of the band for multi-band inputs, e.g. NETCDF.
    subgrid_ids : list, optional
        Only resample to the specified continents.
        Default is to resample to all continents.
    accurate_boundary : bool, optional
        If true (default), the accurate raster boundary of the input image is retrieved.
        If false, only the extent is used.
    e7_folder : bool, optional
        If true (default), the output data will be stored in the Equi7Grid folder structure, i.e. "subgrid/tilename".
    ftiles : list, optional
        List of full name of tiles to which data should be resampled. If it is not set, all tiles are used.
    coverland : bool, optional
        If true (default), only tiles oovering land should be resampled.
    roi : ogr.Geometry, optional
        Region of interest defined by a ogr.Geometry.
        `roi` will be ignored if `ftiles` is provided.
    naming_convention : SmartFilename, optional
        If provided, the keys "grid_name" and "tile_name" of the file naming convention are filled with the
        corresponding values of the resampled tile. The file naming convention is then used to specify the file name
        of the resampled files.
    compress_type : str, optional
        GeoTIFF compression type. Defaults to "LZW".
    resampling_type : str, optional
        GDAL resampling method. Defaults to "bilinear".
    subfolder : str, optional
        If it's set, it creates the sub-folder within the tile folder where the
        image will be resampled to.
    overwrite : bool, optional
        Overwrite the old tile or not (defaults to false).
    data_type : str, optional
        Force the output image bands to have a specific data type supported by the driver,
        which may be one of the following: Byte, UInt16, Int16, UInt32, Int32, Float32, Float64,
        CInt16, CInt32, CFloat32 or CFloat64.
    image_nodata : float, optional
        The no data value of the input image. If it is not set, the no data value is automatically taken
        from the input image (default).
    tile_nodata : float, optional
        The no data value of the tile. If it is not set, the no data value is automatically taken
        from the input image (default).
    scale : float, optional
        Scale factor that should be applied to the pixel values.
    offset : float, optional
        Offset value that should be applied to the pixel values.
    tiled : bool, optional
        If true, tiled blocks are used for the GeoTIFF file (default).
    blocksize : integer, optional
        Sets the GeoTIFF block size (defaults to 512).
        If `tiled=True`, the block size is used for both X and Y dimensions.
        If `tiled=False`, then stripped GeoTIFF files are created, where `blocksize`
        specifies the block size in Y direction.

    """

    # find overlapping tiles
    if ftiles is None:
        if roi is None:
            if accurate_boundary:
                try:
                    geo_extent = retrieve_raster_boundary(image,
                                                        gdal_path=gdal_path,
                                                        nodata=image_nodata)
                except Exception as e:
                    print("retrieve_raster_boundary failed:", str(e))
                    geo_extent = None
            else:
                geo_extent = None
            if geo_extent:
                ftiles = e7grid.search_tiles_in_roi(roi_geometry=geo_extent,
                                                    subgrid_ids=subgrid_ids,
                                                    coverland=coverland)
            else:
                ds = open_image(image)
                img_extent = ds.get_extent()
                bbox = (img_extent[0:2], img_extent[2:4])
                img_spref = osr.SpatialReference()
                img_spref.ImportFromWkt(ds.projection())
                ftiles = e7grid.search_tiles_in_roi(bbox=bbox,
                                                    subgrid_ids=subgrid_ids,
                                                    osr_spref=img_spref,
                                                    coverland=coverland)
        else:
            ftiles = e7grid.search_tiles_in_roi(roi_geometry=roi,
                                                subgrid_ids=subgrid_ids,
                                                coverland=coverland)
    else:
        if type(ftiles) != list:
            ftiles = [ftiles]

    # keep the full path of the output(resampled) files
    dst_file_names = []

    # resample for each tile sequentially
    for ftile in ftiles:
        e7tile = e7grid.create_tile(ftile)
        if coverland:  # skip tiles not covering land
            if not e7tile.covers_land:
                continue
        # create grid folder
        if e7_folder:
            grid_folder = "EQUI7_{}".format(ftile[0:6])
            tile_path = os.path.join(output_dir, grid_folder, ftile[7:])
            if not os.path.exists(tile_path):
                os.makedirs(tile_path)
        else:
            tile_path = output_dir

        # make output filename
        if naming_convention is None:
            out_filename = os.path.splitext(os.path.basename(image))[0]
            out_filename = "_".join((out_filename, ftile + ".tif"))
        else:
            try:
                naming_convention['grid_name'] = ftile.split('_')[0]
                naming_convention['tile_name'] = ftile.split('_')[1]
            except KeyError:
                err_msg = "File naming convention does not contain 'grid_name' or 'tile_name'."
                raise KeyError(err_msg)
            out_filename = str(naming_convention)
        if subfolder:
            out_filepath = os.path.join(tile_path, subfolder, out_filename)
        else:
            out_filepath = os.path.join(tile_path, out_filename)

        # using gdalwarp to resample
        bbox = e7grid.get_tile_bbox_proj(ftile)
        tile_project = '"{}"'.format(e7grid.subgrids[ftile[0:2]].core.projection.proj4)

        # prepare options for gdalwarp
        options = {'-t_srs': tile_project,
                   '-of': 'GTiff',
                   '-r': resampling_type,
                   '-te': " ".join(map(str, bbox)),
                   '-tr': "{} -{}".format(e7grid.core.sampling,
                                          e7grid.core.sampling)}

        options["-co"] = list()
        if compress_type is not None:
            options["-co"].append("COMPRESS={0}".format(compress_type))
        if image_nodata != None:
            options["-srcnodata"] = image_nodata
        if tile_nodata != None:
            options["-dstnodata"] = tile_nodata
        if data_type != None:
            options["-ot"] = data_type
            options["-wt"] = data_type # test if this is what we want
        if overwrite:
            options["-overwrite"] = " "
        if tiled:  # tiled, square blocks
            options["-co"].append("TILED=YES")
            options["-co"].append("BLOCKXSIZE={0}".format(blocksize))
            options["-co"].append("BLOCKYSIZE={0}".format(blocksize))
        else:  # stripped blocks
            blockxsize = e7grid.core.tile_xsize_m // e7grid.core.sampling
            blockysize = blocksize
            options["-co"].append("TILED=NO")
            options["-co"].append("BLOCKXSIZE={0}".format(blockxsize))
            options["-co"].append("BLOCKYSIZE={0}".format(blockysize))

        # call gdalwarp for resampling
        succeed, _ = call_gdal_util('gdalwarp', src_files=image, src_band=inband,
                                    dst_file=out_filepath, gdal_path=gdal_path,
                                    options=options)

        # full path to the output(resampled) files
        if succeed:
            dst_file_names.append(out_filepath)

        if scale is not None and offset is not None:
            # prepare options for gdal_translate
            options = {'-a_scale': scale,
                       '-a_offset': offset}
            options["-co"] = list()
            if tile_nodata != None:
                options["-a_nodata "] = tile_nodata
            if compress_type is not None:
                options["-co"].append("COMPRESS={0}".format(compress_type))
            if blocksize is not None:
                options["-co"].append("TILED=YES")
                options["-co"].append("BLOCKXSIZE={0}".format(blocksize))
                options["-co"].append("BLOCKYSIZE={0}".format(blocksize))

            succeed, _ = call_gdal_util('gdal_translate', src_files=out_filepath,
                                        dst_file=out_filepath, gdal_path=gdal_path,
                                        options=options)

    return dst_file_names


def open_image(filename):
    """ open an image file

    Parameters
    ----------
    filename : string. full path string of input file

    Returns
    -------
    GdalImage
        GdalImage object if successful, otherwise None

    Raise
    ------
    IOError
        if fail to open the image file

    """

    dataset = gdal.Open(filename, gdal.GA_ReadOnly)
    if dataset is None:
        raise IOError("cannot open %s" % filename)

    return GdalImage(dataset, filename)


def call_gdal_util(util_name, gdal_path=None, src_files=None, src_band=None,
                   dst_file=None, options={}):

    """call gdal utility to run the operation.
        http://www.gdal.org/gdal_utilities.html

    Parameters
    ----------
    util_name : string
        pre-defined name of the utility
        (e.g. "gdal_translate": convert raster data between different formats,
        potentially performing some operations like subsettings, resampling,
        and rescaling pixels in the process.)
    src_files : string or list of strings
        The source dataset name(s). It can be either file name(s),
        URL of data source or subdataset name for multi-dataset files.
    dst_file : string
        The destination file name.
    gdal_path : string
        It the path where your gdal installed. If gpt command can not found by
        the os automatically, you should give this path.
    options : dict
        A dictionary of options. You can find all options at
        http://www.gdal.org/gdal_utilities.html

    """
    # define specific options
    _opt_2b_in_quote = ["-mo", "-co"]

    # get the gdal installed path if it is set in system environmental variable
    if not gdal_path:
        gdal_path = _find_gdal_path()
    if not gdal_path:
        raise OSError("gdal utility not found in system environment!")

    # prepare the command string
    cmd = []
    gdal_cmd = os.path.join(gdal_path, util_name) if gdal_path else util_name
    # put gdal_cmd in double quotation
    cmd.append('"%s"' % gdal_cmd)

    for k, v in iter(options.items()):
        if k in _opt_2b_in_quote:
            if (k == "-mo" or k == "-co") and isinstance(v, (tuple, list)):
                for i in range(len(v)):
                    cmd.append(" ".join((k, '"%s"' % v[i])))
            else:
                cmd.append(" ".join((k, '"%s"' % v)))
        else:
            if v is not None:
                cmd.append(k)
            # if hasattr(v, "__iter__"):
            #    cmd.append(' '.join(map(str, v)))
            # else:
            #    cmd.append(str(v))
                cmd.append(str(v))

    # add source files and destination file (in double quotation)
    dst_file_ori = None
    # switch for multiple source files (e.g. for gdal_merge.py)
    if isinstance(src_files, list):
        src_files_str = " ".join(src_files)
    # switch for single source files
    else:
        src_files_str = '"%s"' % src_files
        # NETCDF input case
        if src_files.endswith('.nc'):
            src_files = 'NETCDF:{}:{}'.format(src_files, src_band)
        # create an interim existing file
        if src_files == dst_file:
            fileparts = os.path.splitext(dst_file)
            dst_file_tmp = fileparts[0] + '_temp' + fileparts[1]
            dst_file_ori = dst_file
            dst_file = dst_file_tmp

    # build the final call
    cmd.append(src_files_str)
    cmd.append('"%s"' % dst_file)

    # create the directory if not exists
    if dst_file is not None:
        if not os.path.exists(os.path.dirname(dst_file)):
            os.makedirs(os.path.dirname(dst_file))

    # check for success
    output = subprocess.check_output(" ".join(cmd), shell=True, cwd=gdal_path)
    succeed = _analyse_gdal_output(output)

    # restore old filename
    if succeed and dst_file_ori is not None:
        os.remove(dst_file_ori)
        os.rename(dst_file, dst_file_ori)

    return succeed, output


def _find_gdal_path():
    """find the gdal installed path from the system enviroment variables."""
    evn_name = 'GDAL_DATA'
    # print os.environ[gdal_env_var_name]
    return os.environ[evn_name].replace('\\gdal-data', '') if evn_name in os.environ else None


def _analyse_gdal_output(output):
    """analyse the output from gpt to find if it executes successfully."""

    # return false if "Error" is found.
    if b'error' in output.lower():
        return False
    # return true if "100%" is found.
    elif b'100 - done' in output.lower():
        return True
    # otherwise return false.
    else:
        return False


def retrieve_raster_boundary(infile, gdal_path=None, nodata=None,
                             resize_factor=None):
    """
    Retrieve the boundary of raster image excluding the nodata

    Parameters
    ----------
    infile : filepath
        input raster image file
    gdal_path : string
        gdal utilities path
    nodata : number
        value of no data / background value
    resize_factor : number
        factor for downsampling infile for the quicklook

    Return
    ------
        a multipolygon geometry object

    Note
    ----
    The accurate boundary is retrieved via the following steps:
        1. create the quicklook (an overview)
        2. using morphological operaton to fill the samll gaps in the quicklook
        3. using gdalpolgonize function to get the geometry

    """
    # create overview
    cur_time = datetime.utcnow().strftime("%Y%m%d%H%M%S")
    infile_bname = os.path.splitext(os.path.basename(infile))[0]
    qlook = os.path.join(os.path.dirname(infile),
                         "temp_{}_{}.tif".format(infile_bname, cur_time))

    img_ds = open_image(infile)
    width, height = img_ds.XSize(), img_ds.YSize()
    img_geot = img_ds.geotransform()
    img_extent = img_ds.get_extent()
    img_ds = None
    qlook_size = 400.0 if resize_factor is None else resize_factor
    factor = int(np.round(max(width / qlook_size, height / qlook_size)))
    x_res = img_geot[1] * factor
    y_res = img_geot[5] * factor

    options = {'-of': 'GTiff',
               '-r': 'near',
               '-te': ' '.join(map(str, img_extent)),
               '-tr': '{} -{}'.format(x_res, np.abs(y_res)),
               '-co': 'COMPRESS=LZW',
               '-srcnodata': nodata,
               '-dstnodata': nodata
               }

    succeed, _ = call_gdal_util('gdalwarp', src_files=infile,
                                dst_file=qlook, gdal_path=gdal_path,
                                options=options)
    if not succeed:
        return None

    # make binary image
    qlook_dst = open_image(qlook)
    src_arr = qlook_dst.read_band(1)
    mask = (src_arr != nodata)
    src_arr = None

    # morphologic dilation
    pixels = 3
    struct = ndimage.generate_binary_structure(2, 2)
    struct = ndimage.morphology.iterate_structure(struct, pixels)
    new_mask = ndimage.binary_dilation(mask, structure=struct)
    src_arr = np.zeros_like(mask, dtype=np.uint8)
    src_arr[new_mask] = 1

    # make a memory dataset
    mem_drv = gdal.GetDriverByName('MEM')
    binary_dst = mem_drv.Create('', src_arr.shape[1], src_arr.shape[0], 1,
                                gdal.GDT_Byte)
    binary_dst.SetGeoTransform(qlook_dst.geotransform())
    binary_dst.SetProjection(qlook_dst.projection())
    binary_band = binary_dst.GetRasterBand(1)
    binary_band.WriteArray(src_arr)
    binary_band.SetNoDataValue(0)
    maskband = binary_band.GetMaskBand()

    # Polygonize
    driver = ogr.GetDriverByName('Memory')
    dst_ds = driver.CreateDataSource("")

    srs = None
    if qlook_dst.projection():
        srs = osr.SpatialReference()
        srs.ImportFromWkt(qlook_dst.projection())

    dst_layer = dst_ds.CreateLayer("out", srs=srs)
    fd = ogr.FieldDefn('DN', ogr.OFTInteger)
    dst_layer.CreateField(fd)
    dst_field = 0

    geom = None
    if gdal.CE_None == gdal.Polygonize(binary_band, maskband,
                                       dst_layer, dst_field, callback=None):
        # get polygons from dataset
        multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
        for feature in dst_ds.GetLayerByIndex():
            polygon = feature.GetGeometryRef()

            # When the polygon has hole inside, only the out line ring will 
            # be used as extent polygon
            if polygon.GetGeometryCount() > 1:
                ring = polygon.GetGeometryRef(0)
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(ring)
                polygon = poly

            multipolygon.AddGeometry(polygon)
        multipolygon.AssignSpatialReference(dst_ds.GetLayer(0).GetSpatialRef())

        geom = multipolygon

    # clean tmp file
    dst_ds.Destroy()
    binary_dst = None
    qlook_dst = None
    # remove temp qlook
    os.remove(qlook)

    return geom


def equi72lonlat(e7grid, image, output_dir, gdal_path=None, subgrid_ids=None,
                 accurate_boundary=True, e7_folder=True, ftiles=None,
                 roi=None, outshortname=None,
                 withtilenameprefix=False, withtilenamesuffix=True,
                 compress=True, compresstype="LZW",
                 resampling_type="bilinear",
                 overwrite=False, image_nodata=None, tile_nodata=None,
                 tiledtiff=True, blocksize=512):

    ftile = 'EU010M_E073N032T1'
    tile_path= r'D:\temp\out'
    outshortname = 'aaaaa'
    # make output filename
    outbasename = outshortname
    if withtilenameprefix:
        outbasename = "_".join((ftile, outbasename))
    if withtilenamesuffix:
        outbasename = "_".join((outbasename, ftile))
    filename = os.path.join(tile_path, "".join((outbasename, ".tif")))

    # using gdalwarp to resample

    res500 = 0.0000892857142857142857142857142857
    extent_m = e7grid.get_tile_bbox_geog(ftile)

    # prepare options for gdalwarp
    options = {'-t_srs': 'EPSG:4326', '-of': 'GTiff',
               '-r': 'bilinear', '': '-overwrite',
               '-te': " ".join(map(str, extent_m)),
               '-tr': "{} -{}".format(res500, res500)}

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
    succeed, _ = call_gdal_util('gdalwarp', src_files=image,
                                dst_file=filename, gdal_path=gdal_path,
                                options=options)
