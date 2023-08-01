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
import numpy as np
from equi7grid.equi7grid import Equi7Grid
from geospade.raster import RegularMosaicGeometry
from geospade.raster import Tile
from geospade.crs import SpatialRef


def run1():

    eg = Equi7Grid(500)

    for sg in ['EU']:#eg.subgrids:

        boundary = eg.subgrids[sg].polygon_proj

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

            mosaic_topo = 'OUTER'
            within = t.get_extent_geometry_proj().Within(boundary)
            intersects = t.get_extent_geometry_proj().Intersects(boundary)
            if within:
                mosaic_topo = 'INNER'
            elif intersects:
                mosaic_topo = 'BOUNDARY'

            gst = Tile(n_rows, n_cols, sref, geotrans=geotransform, name=name, mosaic_topology=mosaic_topo,
                       active=tile in at_sh, px_origin='ll',
                       metadata=md)

            #gst.to_json(r"D:\Arbeit\atasks\202307_equi7grid_update\geojson\tests.json")
            tile_objects.append(gst)

        mosaic_name = "{} subgrid".format(sg)


        mosaic_geom = RegularMosaicGeometry.from_tile_list(tile_objects, name=mosaic_name, boundary=boundary)
        mosaic_geom.to_json(r"D:\Arbeit\atasks\202307_equi7grid_update\geojson\{}_5.json".format(sg))
        pass

def run3():

    sampling = 500

    eg = Equi7Grid(sampling)

    tilesize = eg.core.tile_xsize_m
    tilecode = eg.get_tiletype()
    sample_size = int(tilesize / sampling)


    for sg in ['EU']:  # eg.subgrids:

        sref = SpatialRef(eg.subgrids[sg].core.projection.wkt)

        bd = eg.subgrids[sg].polygon_proj
        bbox = tuple([np.round(x, decimals=5) for x in bd.GetEnvelope()])
        n_rows = int(bbox[3] // tilesize + 1)
        n_cols = int(bbox[1] // tilesize + 1)

        gt = (0, sampling, 0, 0, 0, sampling * -1)

        mosaic_name = "{} subgrid".format(sg)



        mg = RegularMosaicGeometry.from_rectangular_definition(n_rows, n_cols, tilesize, tilesize, sref,
                                    geotrans=gt, tile_class=Tile, tile_kwargs=None,
                                    name_frmt="E{:03d}N{:03d}" + tilecode, boundary=bd, name=mosaic_name, description="")

        mg.to_json(r"D:\Arbeit\atasks\202307_equi7grid_update\geojson\{}_6.json".format(sg))
        mg.plot(label_tiles=True, show=True)
        pass

def run2():
    import matplotlib.pyplot as plt
    jsonfile = r"D:\Arbeit\atasks\202307_equi7grid_update\geojson\EU_5.json"

    mosaicgeom = RegularMosaicGeometry.from_json(jsonfile)


    #mosaicgeom.poi2tile()
    tiles = mosaicgeom.get_neighbouring_tiles('E066N030T6', active_only=True)
    #plt.imshow(tiles['E066N036T6'].mask)
    #plt.show()
    mosaicgeom.plot(label_tiles=True, show=True)
    RegularMosaicGeometry.from_rectangular_definition()
    pass

if __name__ == '__main__':
    run1()



