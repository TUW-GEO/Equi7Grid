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
Created on March 7, 2017

Code for osgeo geometry operations.

@author: Bernhard Bauer-Marschallinger, bbm@geo.tuwien.ac.at

'''






from osgeo import ogr
from osgeo import osr



def uv2xy(src_ref, dst_ref, u, v):
    # tranform the point
    tx = osr.CoordinateTransformation(src_ref, dst_ref)
    x, y, _ = tx.TransformPoint(u, v)
    return x, y

def create_multipoint_geom(u, v, projection):
    geog_spref = projection.osr_spref
    point_geom = ogr.Geometry(ogr.wkbMultiPoint)
    point_geom.AssignSpatialReference(geog_spref)
    for p, _ in enumerate(u):
        point = ogr.Geometry(ogr.wkbPoint)
        point.SetPoint(0, u[p], v[p])
        point_geom.AddGeometry(point)


    return point_geom

def create_point_geom(u, v, projection):
    geog_spref = projection.osr_spref
    point_geom = ogr.Geometry(ogr.wkbPoint)
    point_geom.AddPoint(u, v)
    point_geom.AssignSpatialReference(geog_spref)

    return point_geom


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


def open_geometry(fname, format="shapefile"):
    '''
    opens a geometry from a vector file.

    fname : string
        full path of the output file name
    format : string
        format name. currently only shape file is supported
    def __getitem__(key):
        return geoms[key]
    '''
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.Open(fname, 0)
    feature = ds.GetLayer(0).GetFeature(0)
    geom = feature.GetGeometryRef()

    out = geom.Clone()
    ds, feature, geom, = None, None, None
    return out

def write_geometry(geom, fname, format="shapefile"):
    """ write a geometry to a vector file.

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
    #dst_field = 0

    feature = ogr.Feature(dst_layer.GetLayerDefn())
    feature.SetField("DN", 1)
    feature.SetGeometry(geom)
    dst_layer.CreateFeature(feature)
    feature.Destroy()
    # clean tmp file
    dst_ds.Destroy()
    return

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

def extent2polygon(extent, wkt=None):
    """create a polygon geometry from extent.

    extent : list
        extent in terms of [xmin, ymin, xmax, ymax]
    wkt : string
        project string in well known text format

    """
    area = [(extent[0], extent[1]), ((extent[0] + extent[2])/2, extent[1]), (extent[2], extent[1]),
            (extent[2], (extent[1] + extent[3])/2),
            (extent[2], extent[3]), ((extent[0] + extent[2])/2, extent[3]), (extent[0], extent[3]),
            (extent[0], (extent[1] + extent[3])/2)]
    edge = ogr.Geometry(ogr.wkbLinearRing)
    [edge.AddPoint(x, y) for x, y in area]
    edge.CloseRings()
    geom_area = ogr.Geometry(ogr.wkbPolygon)
    geom_area.AddGeometry(edge)
    if wkt:
        geo_sr = osr.SpatialReference()
        geo_sr.ImportFromWkt(wkt)
        geom_area.AssignSpatialReference(geo_sr)
    return geom_area