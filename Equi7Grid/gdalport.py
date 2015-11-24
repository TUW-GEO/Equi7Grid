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

A sample class and Some handy functions for using gdal library to read/write
data

@author: Bernhard Bauer-Marschallinger, bbm@geo.tuwien.ac.at

'''


import os
import subprocess
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from datetime import datetime
import numpy as np
from scipy import ndimage  # for extract raster boundary


class GdalImage:

    """A sample class to access a image with GDAL library"""

    def __init__(self, gdaldataset):
        self.dataset = gdaldataset

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
            The subset should be in pixles, like this (xmin, ymin, xmax, ymax).

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

    def get_raster_nodata_value(self, band_idx=1):
        """get the nodatave.

        Parameters
        ----------
        band_idx : int
            the band index starting from 1

        Returns
        -------
        nodata
            no data value if it's available, otherwise it will return None

        """
        if band_idx < 1 or band_idx > self.dataset.RasterCount:
            raise IndexError("band index is out of range")
        return self.dataset.GetRasterBand(band_idx).GetNoDataValue()

    def read_all_band(self):
        """read the data of all the bands"""
        raise NotImplementedError(
            "sorry, read_all_band has not been implemented!")

    def get_raster_dtype(self, band_idx=1):
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

    def band_count(self):
        """get the band count"""
        return self.dataset.RasterCount if self.dataset else None

    def get_extent(self):
        """get the extent of the image as (xmin, ymin, xmax, ymax)."""
        geot = self.geotransform()
        return (geot[0], geot[3] + self.YSize() * geot[5],
                geot[0] + self.XSize() * geot[1], geot[3])


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

    return GdalImage(dataset)


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


def get_gdal_datatype(datatype):
    """get gdal data type from datatype

    Parameters
    ----------
    datatype :  string.
            data type string in python such as "uint8", "float32" and so forth.

    Returns
    -------
    gdal_datatype
            gdal data type
    """
    datatype = datatype.lower()
    return gdal_datatype[datatype] if datatype in gdal_datatype else None


def get_gdal_resampling_type(resampling_method):
    """get gdal resample type from resampling method string"""
    mthd = resampling_method.lower()
    return gdal_resample_type[mthd] if mthd in gdal_resample_type else None


def create_dataset(filename, datatype, dims, frmt="GTiff", geotransform=None,
                   projection=None, option=None):
    """create GDAL dataset.

    Parameters
    ----------
    filename : string
        full path of output filename
    datatype : string
        data type string like numpy's dtpye
    dims : tuple
        Dimension of the dataset in the format of (bands, XSize, YSize)
    frmt :  string
        The format of output image should be a string that gdal supported
    geotransform : array like
        It contains six geotransform parameters
    projection : string
        projection definition string

    Returns
    -------
    GDAL dataset

    Raise:
    ------
    IOError
        if fail to obtain driver with specific format or to create the output
        dataset

    """
    driver = gdal.GetDriverByName(frmt)
    gdaldatatype = get_gdal_datatype(datatype)
    if driver is None:
        raise IOError("cannot get driver of {}".format(frmt))
    band_count, xsize, ysize = dims
    if option is None:
        out_ds = driver.Create(
            filename, xsize, ysize, band_count, gdaldatatype)
    else:
        out_ds = driver.Create(filename, xsize, ysize, band_count,
                               gdaldatatype, option)
    if out_ds is None:
        raise IOError("cannot create file of {}".format(filename))
    if geotransform is not None:
        out_ds.SetGeoTransform(geotransform)
    if projection is not None:
        out_ds.SetProjection(projection)
    return out_ds


def write_image(image, filename, frmt="GTiff", nodata=None,
                geotransform=None, projection=None, option=None):
    """output image into filename with specific format

    Parameters
    ----------
    image : array like
        two dimension array containing data that will be stored
    filename : string
        full path of output filename
    frmt :  string
        the format of output image should be a string that gdal supported.
    nodata : list, optional
        a list contian the nodata values of each band
    geotransform : array like
        contain six geotransform parameters
    projection : string
        projection definition string

    Returns
    -------
    dataset
        the gdal dataset object of output dataset

    Raise
    -----
        IOError
            if IO error happens
        ValueError
            if some invalid parameters are given

    """
    # to make sure dim of image is 2 or 3
    if image is None or image.ndim < 2 or image.ndim > 3:
        raise ValueError(
            "The image is None or it's dimension isn't in two or three.")
    dims = (1, image.shape[1], image.shape[0]) if image.ndim == 2 \
        else (image.shape[2], image.shape[1], image.shape[0])
    # create dataset
    ds = create_dataset(filename, str(image.dtype), dims, frmt, geotransform,
                        projection, option)
    # write data
    if image.ndim == 2:
        ds.GetRasterBand(1).WriteArray(image, 0, 0)
    else:
        for i in xrange(ds.RasterCount):
            ds.GetRasterBand(i + 1).WriteArray(image[:, :, i], 0, 0)
    # set nodata for each band
    if nodata is not None:
        for i in xrange(ds.RasterCount):
            ds.GetRasterBand(i + 1).SetNoDataValue(nodata[i])

    return ds


def write_geometry(geom, fname, format="shapefile"):
    """ write an geometry to a vector file.

    parameters
    ----------
    geom : Geometry
        geometry object
    fname : string
        full path of the output file name
    format : string
        format name. currently only shape file is supported
    """

    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(fname)
    srs = geom.GetSpatialReference()
    
    dst_layer = dst_ds.CreateLayer("out", srs=srs)
    fd = ogr.FieldDefn('DN', ogr.OFTInteger)
    dst_layer.CreateField(fd)
    dst_field = 0
    
    feature = ogr.Feature(dst_layer.GetLayerDefn())
    feature.SetField("DN", 1)
    feature.SetGeometry(geom)
    dst_layer.CreateFeature(feature)
    feature.Destroy()
    # clean tmp file
    dst_ds.Destroy()
    return 


def extent2polygon(extent):
    """create a polygon geometry from extent."""
    area = [(extent[0], extent[1]), (extent[2], extent[1]),
            (extent[2], extent[3]), (extent[0], extent[3])]
    edge = ogr.Geometry(ogr.wkbLinearRing)
    [edge.AddPoint(x, y) for x, y in area]
    edge.CloseRings()
    geom_area = ogr.Geometry(ogr.wkbPolygon)
    geom_area.AddGeometry(edge)
    return geom_area


def call_gdal_util(util_name, gdal_path=None, src_files=None, dst_file=None,
                   options={}):
    """call gdal utility to run the operation.
        http://www.gdal.org/gdal_utilities.html

    Parameters
    ----------
    util_name : string
        pre-defined name of the utility
        (e.g. "gdal_translate": convert raster data between different formats,
        potentially performing some operations like subsettings, resampling,
        and rescaling pixels in the process.)
    src_files : string
        The source dataset name. It can be either file name,
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

    for k, v in options.iteritems():
        if k in _opt_2b_in_quote:
            if (k == "-mo" or k == "-co") and isinstance(v, (tuple, list)):
                for i in range(len(v)):
                    cmd.append(" ".join((k, '"%s"' % v[i])))
            else:
                cmd.append(" ".join((k, '"%s"' % v)))
        else:
            cmd.append(k)
            if hasattr(v, "__iter__"):
                cmd.append(' '.join(map(str, v)))
            else:
                cmd.append(str(v))

    # add source files and destination file (in double quotation)
    if hasattr(src_files,  "__iter__"):
        src_files_str = " ".join(src_files)
    else:
        src_files_str = '"%s"' % src_files
    cmd.append(src_files_str)
    cmd.append('"%s"' % dst_file)

    output = subprocess.check_output(" ".join(cmd), shell=True, cwd=gdal_path)
    succeed = _analyse_gdal_output(output)
    return succeed, output


def _find_gdal_path():
    """find the gdal installed path from the system enviroment variables."""
    evn_name = "gdal_util_HOME"
    # print os.environ[gdal_env_var_name]
    return os.environ[evn_name] if evn_name in os.environ else None


def _analyse_gdal_output(output):
    """analyse the output from gpt to find if it executes successfully."""

    # return false if "Error" is found.
    if 'error' in output.lower():
        return False
    # return true if "100%" is found.
    elif '100 - done' in output.lower():
        return True
    # otherwise return false.
    else:
        return False


def retrieve_raster_boundary(infile, gdal_path=None, nodata=None):
    """return a geometry of multipolygon representing the accurate boundary.

    Parameters
    ----------
    infile : filepath
        input raster image file
    gdal_path : string
        fullpath of the gdal utilities


    Return
    ------
        a multipolygon geometry object
    Note
    ----
    retrieve the accurate boundary ploygon in the following step:
        1. create the qlook if not provided.
        2. using morphological operaton to fill the samll gaps in the qlook
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
    qlook_size = 400.0
    factor = int(np.round(max(width/qlook_size, height/qlook_size)))
    x_res = img_geot[1] * factor
    y_res = img_geot[5] * factor
    
    options = {'-of': 'GTiff',
               '-r': "near",
               '-te': " ".join(map(str, img_extent)),
               '-tr': "{} -{}".format(x_res, np.abs(y_res)),
               '-co': "COMPRESS=LZW",
               '-srcnodata' : nodata,
               '-dstnodata' : nodata
               }


#    options = {'-of': 'GTiff', '-co': 'COMPRESS=LZW',
#               '-outsize': ("2%", "2%"), '-ot': 'Byte',
#               '-scale': (-9999, 9999, 0, 255)}

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

    # src_arr[mask] = 1

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

#    outfile = os.path.join(os.path.dirname(infile),
#             "{}_{}.shp".format(
#                 os.path.splitext(os.path.basename(infile))[0],
#                 datetime.now().strftime("%Y%m%d_%H%M%S")))
#    drv = ogr.GetDriverByName("ESRI Shapefile")
#    dst_ds = drv.CreateDataSource(outfile)

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
        for feature in dst_ds.GetLayer():
            multipolygon.AddGeometry(feature.GetGeometryRef())
        multipolygon.AssignSpatialReference(dst_ds.GetLayer(0).GetSpatialRef())
        geom = multipolygon

    # clean tmp file
    dst_ds.Destroy()
    binary_dst = None
    qlook_dst = None
    # remove temp qlook
    os.remove(qlook)

    return geom

