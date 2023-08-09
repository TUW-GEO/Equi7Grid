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

from osgeo import ogr

class Equi7Tile(Tile):

    def __init__(self, n_rows, n_cols, sref,
                 geotrans=(0, 1, 0, 0, 0, -1),
                 mosaic_topology="INNER",
                 active=True,
                 metadata=None,
                 name=None,
                 description="",
                 px_origin="ul",
                 parent=None):

        super().__init__(n_rows, n_cols, sref,
                         geotrans=geotrans,
                         mosaic_topology=mosaic_topology,
                         active=active,
                         metadata=metadata,
                         name=name,
                         description=description,
                         px_origin=px_origin,
                         parent=parent)

    def to_definition(self) -> dict:
        tile_dict = super().to_definition()
        del tile_dict['spatial_reference']
        del tile_dict['pixel_origin']
        del tile_dict['number_of_rows']
        del tile_dict['number_of_columns']

        return tile_dict

class Equi7Subgrid(RegularMosaicGeometry):

    def __init__(self, tiles, boundary=None, adjacency_matrix=None, name="",
                 description="",
                 check_consistency=True,
                 **kwargs):

        super().__init__(tiles, boundary=boundary, adjacency_matrix=adjacency_matrix, name=name,
                         description=description,
                         check_consistency=check_consistency, **kwargs)

    def to_definition(self) -> dict:
        subgrid_dict_temp = super().to_definition()
        ref_tile = self.tiles[0]
        subgrid_dict_temp['spatial_reference'] = ref_tile.sref.to_proj4_dict()
        subgrid_dict_temp['spatial_reference_wkt'] = ref_tile.sref.wkt
        subgrid_dict_temp['pixel_origin'] = ref_tile.px_origin
        subgrid_dict_temp['number_of_rows'] = ref_tile.n_rows
        subgrid_dict_temp['number_of_columns'] = ref_tile.n_cols

        grid_header_atts = ['name', 'type', 'tile_class', 'number_of_rows', 'number_of_columns',
                            'pixel_origin', 'description',
                            'spatial_reference', 'spatial_reference_wkt']

        subgrid_dict = {k: None for k in grid_header_atts}
        subgrid_dict.update(subgrid_dict_temp)
        return subgrid_dict

    @classmethod
    def from_definition(cls, definition, check_consistency=True) -> "Equi7Subgrid":
        """
        Creates a mosaic geometry from a human-readable mosaic definition, which is a
        dictionary containing the following elements:
            - 'type'
            - 'boundary'
            - 'name'
            - 'description'
            - 'adjacency_matrix'
            - 'tile_class'
            - 'tiles'

        The expected values can be taken from the `MosaicGeometry`, `Tile`, and `RasterGeometry` constructor docs.

        Parameters
        ----------
        definition : dict
            Human-readable definition of a mosaic.
        check_consistency : bool, optional
            If True, the tiles are checked for consistency, i.e. to be non-overlapping (defaults to True).

        Returns
        -------
        geospade.raster.MosaicGeometry

        """
        mosaic_type = definition['type']
        if mosaic_type != cls._type:
            err_msg = "Mosaic type of definition '{}' does not match expected mosaic type '{}'".format(mosaic_type, cls._type)
            raise ValueError(err_msg)
        mosaic_boundary = ogr.CreateGeometryFromWkt(definition['boundary'])
        mosaic_name = definition['name']
        description = definition['description']
        adjacency_matrix = None if definition['adjacency_matrix'] is None else np.array(definition['adjacency_matrix'])

        tiles = []
        tile_class_name = definition['tile_class']
        tile_class = globals().get(tile_class_name)
        if tile_class is None:
            err_msg = "Tile class '{}' must be imported.".format(tile_class_name)
            raise ImportError(err_msg)

        for key in definition['tiles'].keys():
            tile_def = definition['tiles'][key]
            tile_def['spatial_reference'] = definition['spatial_reference']
            tile_def['pixel_origin'] = definition['pixel_origin']
            tile_def['number_of_rows'] = definition['number_of_rows']
            tile_def['number_of_columns'] = definition['number_of_columns']
            tiles.append(tile_class.from_definition(tile_def))

        mosaic_boundary.AssignSpatialReference(tiles[0].sref.osr_sref)
        return cls.from_tile_list(tiles, boundary=mosaic_boundary, adjacency_matrix=adjacency_matrix, name=mosaic_name, description=description,
                                  check_consistency=check_consistency)


    def from_rectangular_definition(cls, my_keywords) -> "Equi7Subgrid":

        pass

        return cls.from_tile_list(tile_objects, name=mosaic_name, boundary=boundary)

def run1():

    eg = Equi7Grid(10)

    for sg in ['AS']:#eg.subgrids:

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

            gst = Equi7Tile(n_rows, n_cols, sref, geotrans=geotransform, name=name, mosaic_topology=mosaic_topo,
                       active=tile in at_sh, px_origin='ul',
                       metadata=md)

            #gst.to_json(r"D:\Arbeit\atasks\202307_equi7grid_update\geojson\tests.json")
            tile_objects.append(gst)

        mosaic_name = "{} subgrid".format(sg)


        mosaic_geom = Equi7Subgrid.from_tile_list(tile_objects, name=mosaic_name, boundary=boundary)
        mosaic_geom.to_json(r"D:\Arbeit\atasks\202307_equi7grid_update\geojson\{}_13.json".format(sg))
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

def run4():

    a = Equi7Subgrid.from_json(r"D:\Arbeit\atasks\202307_equi7grid_update\geojson\EU_13.json")

    a.plot(show=True, label_tiles=True)

    pass

if __name__ == '__main__':
    run1()



