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
Checks the accuracy of the projection from gdal/proj4
"""

import os
import tempfile
import numpy as np
from osgeo import osr
from osgeo import gdal

from equi7grid.image2equi7grid import call_gdal_util
from equi7grid.image2equi7grid import open_image
from equi7grid.image2equi7grid import _find_gdal_path
from equi7grid.equi7grid import Equi7Grid


def check_gdal_aeqd_accuracy(quiet=False,
                             check_gdalwarp=True,
                             lib_dir=None,
                             gdal_path=None):

    print("\nTesting GDAL used via python-bindings.")
    print("Version: {}".format(gdal.__version__))

    # check if the osgeo with GDAL is accurate
    # test point in lat/lon projection
    points = [(-31.627336, 30.306273), (-14.589038, -43.880131),
              (79.423313, -35.261658)]
    geo_sr = osr.SpatialReference()
    geo_sr.SetWellKnownGeogCS("EPSG:4326")

    aeqd_proj = (
        'PROJCS["Azimuthal_Equidistant",GEOGCS["GCS_WGS_1984",'
        'DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],'
        'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],'
        'PROJECTION["Azimuthal_Equidistant"],'
        'PARAMETER["false_easting",5621452.01998],'
        'PARAMETER["false_northing",5990638.42298],'
        'PARAMETER["longitude_of_center",21.5],'
        'PARAMETER["latitude_of_center",8.5],UNIT["Meter",1.0]]')
    aeqd_sr = osr.SpatialReference()
    aeqd_sr.ImportFromWkt(aeqd_proj)

    test_results = []
    for pt in points:
        # from lat/lon to aeqd
        forward_tx = osr.CoordinateTransformation(aeqd_sr, geo_sr)
        (aeqd_x, aeqd_y, _) = forward_tx.TransformPoint(pt[0], pt[1])
        # from aeqd to lat/lon
        backward_tx = osr.CoordinateTransformation(geo_sr, aeqd_sr)
        (geo_x, geo_y, _) = backward_tx.TransformPoint(aeqd_x, aeqd_y)

        if abs(geo_x - pt[0]) < 1e-6 and abs(geo_y - pt[1]) < 1e-6:
            test_results.append(True)
        else:
            test_results.append(False)

    is_osgeo_accurate = all(test_results)
    if not is_osgeo_accurate:
        error_message = (
            "Error: The proj4 library used by GDAL-python-bindings is "
            "not accurately reprojecting between WGS84-lon/lat "
            "and Equi7Grid! Please updated GDAL, or "
            "contact GEO TUWien for the modified proj4 library.")
        raise RuntimeError(error_message)
    else:
        if not quiet:
            print(
                "Success: The proj4 library used by GDAL-python-bindings is "
                "accurately reprojecting between WGS84-lon/lat and Equi7Grid!")

    # check the gdalwarp is accurate or not
    if check_gdalwarp:

        _, gdal_cmd_version = call_gdal_util("gdalinfo",
                                             gdal_path=gdal_path,
                                             options={'': '--version'})
        pos_str = str(gdal_cmd_version).find('GDAL')
        version = str(gdal_cmd_version)[pos_str + 5:pos_str + 10]

        if not gdal_path:
            gdal_path = _find_gdal_path()
        print("\nTesting GDAL used via cmd-line at path {}".format(gdal_path))
        print("Version: {}".format(version))

        # 10m in degrees
        sampling_10m = 0.0000892857142857142857142857142857

        test_data_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), r"prj_accuracy_test")
        input_image = os.path.join(test_data_dir,
                                   r"lake_in_russia_equi7grid.tif")
        reprojected_image = tempfile.NamedTemporaryFile(suffix=".tif",
                                                        delete=False).name
        re_reprojected_image = tempfile.NamedTemporaryFile(suffix=".tif",
                                                           delete=False).name

        expected_image = os.path.join(test_data_dir,
                                      r"lake_in_russia_lonlat.tif")

        opt = {
            "-co": ["COMPRESS=LZW"],
            "-t_srs": '"' + geo_sr.ExportToProj4() + '"',
            '-r': "bilinear",
            '-tr': "{} -{}".format(sampling_10m, sampling_10m),
            "-dstnodata": -9999,
            "-overwrite": " "
        }

        succeed, _ = call_gdal_util("gdalwarp",
                                    src_files=input_image,
                                    dst_file=reprojected_image,
                                    gdal_path=gdal_path,
                                    options=opt)

        if not succeed:
            raise RuntimeError("Error: gdalwrap is not working!")

        # check the reproject accuracy
        reprojected_data = open_image(reprojected_image).read_band(1)

        expected_data = open_image(expected_image).read_band(1)
        diff_data = np.abs(expected_data - reprojected_data)
        if diff_data.sum() > 100 or diff_data.max() > 10:
            raise RuntimeError(
                "Error: gdalwarp cmd-line-tool is not accurately reprojecting images to lonlat!"
            )

        # reproject back to inital Equi7Grid format
        e7g = Equi7Grid(10)
        opt = {
            "-co": ["COMPRESS=LZW"],
            "-t_srs": '"' + e7g.AS.core.projection.proj4 + '"',
            '-r': "bilinear",
            '-tr': "{} -{}".format(10, 10),
            '-te': "2116900 6970860 2125490 6978790",
            "-srcnodata": -9999,
            "-dstnodata": -9999,
            "-overwrite": " "
        }

        succeed, _ = call_gdal_util("gdalwarp",
                                    src_files=reprojected_image,
                                    dst_file=re_reprojected_image,
                                    gdal_path=gdal_path,
                                    options=opt)

        # check the reproject accuracy
        re_reprojected_data = open_image(re_reprojected_image).read_band(1)
        indnan = np.where(re_reprojected_data != -9999)

        input_data = open_image(input_image).read_band(1)
        diff_data = input_data[indnan] - re_reprojected_data[indnan]
        if np.abs(diff_data.mean()) > 0.001 or np.abs(
                diff_data).sum() > 500000 or np.abs(diff_data).max() > 25:
            raise RuntimeError(
                "Error: gdalwarp cmd-line-tool is not accurately reprojecting images back to Equi7Grid!"
            )

        os.remove(reprojected_image)
        os.remove(re_reprojected_image)

        print(
            "Success: gdalwarp cmd-line-tool is accurately (re-) projecting images!"
        )


if __name__ == "__main__":
    check_gdal_aeqd_accuracy()
