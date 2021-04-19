import numpy as np
import os
from datetime import datetime
from veranda.io.geotiff import GeoTiffFile
from equi7grid.equi7grid import Equi7Grid
from geographiclib.geodesic import Geodesic
from geopathfinder.naming_conventions.yeoda_naming import YeodaFilename


def run(outpath, sampling):

    # initialise objects
    e7grid = Equi7Grid(sampling)
    ellaxis = Geodesic.WGS84.a

    # loop over the continents
    for sg in e7grid.subgrids.keys():

        subgrid = e7grid.subgrids[sg]

        sg_str = sg + str(sampling).zfill(3) + 'M'
        sg_path = os.path.join(outpath, 'EQUI7_' + sg_str)
        if not os.path.exists(sg_path):
            os.makedirs(sg_path)

        # get the projected projection centre
        fe = subgrid.core.projection.osr_spref.GetProjParm('false_easting')
        fn = subgrid.core.projection.osr_spref.GetProjParm('false_northing')

        # loop over the tiles that cover land
        tiles = subgrid.tilesys.list_tiles_covering_land()
        for t in tiles:

            t_path = os.path.join(sg_path, t)
            if not os.path.exists(t_path):
                os.makedirs(t_path)

            yfields = {'var_name': 'K-DISTORTION',
                       'datetime_1': '',
                       'datetime_2': '',
                       'band': '',
                       'extra_field': '',
                       'tile_name': t,
                       'grid_name': sg_str,
                       'version_run_id': 'V01R01',
                       'sensor_field': 'EQUI7GRIDV14'}
            filename = YeodaFilename(yfields)

            tile = subgrid.tilesys.create_tile(t)
            gt = tile.geotransform()

            # create all pixels in projected coordinates
            xv, yv = np.meshgrid(np.arange(gt[0] + gt[1]/2, gt[0] + tile.core.tile_xsize_m, gt[1]),
                                 np.arange(gt[3] + gt[5]/2, gt[3] - tile.core.tile_ysize_m, gt[5]))

            # create the distances to the projection centre
            dists = np.sqrt((xv - fe)**2 + (yv - fn)**2)

            # apply equation for distortion in direction perpendicular to the radius, k:
            # k = c/geod.a / np.sin(c/geod.a)
            #
            # is it just about the distance to the centre (c), and as they should be equally long
            # on the ellipsoid and on the projected plane (the core of of AEQD!), this should be fine.
            ks = dists / ellaxis / np.sin(dists / ellaxis)


            md = {}
            md['creation_date'] = datetime.now()
            md['scale_factor'] = 1.0/10000
            md['add_offset'] = 0
            md['nodata'] = -9999

            with GeoTiffFile(os.path.join(t_path, filename._build_fn()), mode='w', count=1,
                             gdal_opt={'COMPRESS': 'ZSTD'},
                             geotransform=gt, spatialref=tile.core.projection.wkt) as gt_file:
                gt_file.write(np.round(ks*10000).astype(np.int16), band=1,
                              scale_factor=md['scale_factor'],
                              add_offset=md['add_offset'],
                              nodata=md['nodata'])


if __name__ == '__main__':

    outpath = r'D:\Arbeit\atasks\202104_equi7grid_distortion_plots_2'
    sampling = 250

    run(outpath, sampling)