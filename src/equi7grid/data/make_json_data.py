# Copyright (c) 2023, TU Wien, Department of Geodesy and Geoinformation (GEO).
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




"""

from equi7grid.equi7grid import Equi7Grid
from geospade.raster import RegularMosaicGeometry
from geospade.raster import Tile
from geospade.crs import SpatialRef


def run1():

    eg = Equi7Grid(500)

    for sg in ['EU']:#eg.subgrids:

        # all tiles covering the full AF BB
        x1, x2, y1, y2 = list(eg.subgrids[sg].get_bbox_proj())
        at_bb = eg.subgrids[sg].tilesys.identify_tiles_overlapping_xybbox([x1, y1, x2, y2])
        at_bb.sort()
        # 1248 for AF

        at_sh = eg.subgrids[sg].search_tiles_over_geometry(eg.subgrids[sg].polygon_geog, coverland=False)
        at_sh.sort()
        # 1004 for AF

        at_cl = eg.subgrids[sg].search_tiles_over_geometry(eg.subgrids[sg].polygon_geog)
        at_cl.sort()
        # 498 for AF

        tile_objects = []
        for tile in at_bb:

            t = eg.subgrids[sg].tilesys.create_tile(tile)

            n_rows, n_cols = t.shape_px()
            geotransform = t.geotransform()
            name = t.shortname
            sref = SpatialRef(t.core.projection.wkt)
            md = {'covers_land': tile in at_cl}

            gst = Tile(n_rows, n_cols, sref, geotrans=geotransform, name=name, active=tile in at_sh, metadata=md)

            #gst.to_json(r"D:\Arbeit\atasks\202307_equi7grid_update\geojson\tests.json")
            tile_objects.append(gst)

        mosaic_name = "{} subgrid".format(sg)
        mosaic_geom = RegularMosaicGeometry.from_tile_list(tile_objects, name=mosaic_name, boundary=eg.subgrids[sg].polygon_geog)
        mosaic_geom.to_json(r"D:\Arbeit\atasks\202307_equi7grid_update\geojson\{}.json".format(sg))

        pass


def run2():

    jsonfile = r"D:\Arbeit\atasks\202307_equi7grid_update\geojson\EU.json"

    mosaicgeom = RegularMosaicGeometry.from_json(jsonfile)

    mosaicgeom.plot(label_tiles=True, show=True, plot_boundary=False)
    pass

if __name__ == '__main__':
    run1()



