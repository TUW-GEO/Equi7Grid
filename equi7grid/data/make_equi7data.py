#! /usr/bin/env python 
# Copyright (c) 2018, Vienna University of Technology (TU Wien), Department of
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
Created on July 10, 2018

make equi7grid.dat file for Equi7Grid class

@author: Senmao Cao, Senmao.Cao@geo.tuwien.ac.at
'''


import os
import argparse
import pickle
from osgeo import ogr, osr


def make_equi7data(outpath, version="V13"):
    """ Make the equi7grid.dat file

    Parameters
    ----------
    outpath : string
        output file directory path.
    
    Returns
    -------
    int
        0 if succeeded, otherwise error code

    Notes
    -----
    equi7grid.dat is a dictionary including neccessary information required by
    Equi7Grid.py class in the following structure.
    { "AF": { "projection": "projection in wkt format",
              "zone_extent": "subgrid zone geometry in wkt format",
              "coverland": {"T6": set(), # A set includes all tiles covering land
                            "T3": set(),
                            "T1": set(),
                            }
            }
        "AN" : ...
    
    }

    """
    
    outfile = os.path.join(outpath, "equi7grid.dat")
    if os.path.exists(outfile):
        print "Error: File Already Exist!"
        return -1
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    module_path = os.path.dirname(os.path.abspath(__file__))
    grids_path = os.path.join(os.path.dirname(module_path), "grids")
    
    subgrids = ["AF", "AN", "AS", "EU", "NA", "OC", "SA"] 
    tilecodes = ["T1", "T3", "T6"]

    equi7_data = dict()

    for subgrid in subgrids:
        subgrid_data = dict()

        zone_fpath = os.path.join(grids_path, subgrid, "GEOG",
                                  "EQUI7_{}_{}_GEOG_ZONE.shp".format(version, subgrid))
        zone_extent = load_zone_extent(zone_fpath)
        subgrid_data["zone_extent"] = zone_extent.ExportToWkt()
       
        subgrid_data["coverland"] = dict()
        for tilecode in tilecodes:
            tilepath = os.path.join(grids_path, subgrid, "GEOG",
                                    "EQUI7_{}_{}_GEOG_TILE_{}.shp".format(version, subgrid, tilecode))
            tiles_coversland = load_coverland_tiles(tilepath)
            subgrid_data["coverland"][tilecode] = set(tiles_coversland)
        
        # Use spatial reference of T1 tile as the subgrid spatial reference
        sr_path = os.path.join(grids_path, subgrid, "PROJ", "EQUI7_{}_{}_PROJ_TILE_T1.shp".format(version, subgrid))
        sr_wkt = load_spatial_reference(sr_path)
        subgrid_data["projection"] = sr_wkt

        equi7_data[subgrid] = subgrid_data 
    
    # Serialize equi7 data by pickle with protocal=2
    with open(outfile, "wb") as f:
        pickle.dump(equi7_data, f, protocol=2)

    return 0


def load_zone_extent(zone_fpath):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.Open(zone_fpath, update=False)
    layer = ds.GetLayer()
    num_features = layer.GetFeatureCount()
    if num_features < 0:
        raise ValueError("No feature found in {}".format(zone_fpath))
    elif num_features == 1:
        f = layer.GetFeature(0)
        geom = f.GetGeometryRef().Clone()
    else:
        geom = ogr.Geometry(ogr.wkbMultiPolygon)
        for f in layer:
            geom.AddGeometry(f.GetGeometryRef().Clone())
    return geom
    

def load_coverland_tiles(tile_fpath):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.Open(tile_fpath, update=False)
    layer = ds.GetLayer()
    num_features = layer.GetFeatureCount()
    if num_features < 0:
        raise ValueError("No features found in {}".format(tile_fpath))
    
    interval = 100000
    tiles_coversland = list()
    for f in layer:
        coversland = f.GetField("COVERSLAND")
        if coversland:
            extent = int(f.GetField("EXTENT"))
            east = int(f.GetField("EASTINGLL"))
            north = int(f.GetField("NORTHINGLL"))
            tilename = "E{:03d}N{:03d}T{}".format(east/interval, north/interval, extent/interval)
            tiles_coversland.append(tilename)

    return tiles_coversland
    

def load_spatial_reference(fpath):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.Open(fpath, update=False)
    sr_wkt = ds.GetLayer().GetSpatialRef().ExportToWkt()
    return sr_wkt


def main():
    parser = argparse.ArgumentParser(description='Make Equi7Grid Data File')
    parser.add_argument("outpath", help="output folder")
    parser.add_argument("-v", "--version", dest="version", nargs=1, metavar="", help="Equi7 Grids Version. Default is V13.")
    args = parser.parse_args()
    
    outpath = os.path.abspath(args.outpath)
    version = args.version[0] if args.version else "V13"
    return make_equi7data(outpath, version)


if __name__ == "__main__":
    import sys
    sys.argv.append("C:\code\TPS\Equi7Grid\equi7grid\data")
    main()

