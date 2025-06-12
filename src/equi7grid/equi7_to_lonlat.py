import os
import glob
from collections import OrderedDict

import numpy as np
from osgeo import gdal
from equi7grid import Equi7Grid
from image2equi7grid import call_gdal_util, open_image



def chunk_calc(method, nodata, input_files):
    """
    Function to chunkwise perform pixelwise mean/min/max calculation of rasters.
    Overwrites first input file with result.

    Parameters
    ----------
    method : str
        'mean', 'min' or 'max' are supported.
    nodata : int/float
        nodata mask value will be ignored for calculations.
    input_files : list
        list of files to be used as input
        (Must have identical extent/geotransform/proj etc., but no check is done)
    """

    datasets = [gdal.Open(input_files[0], gdal.GA_Update)] + [gdal.Open(path) for path in input_files[1:]]
    bands = [ds.GetRasterBand(1) for ds in datasets]

    x_size, y_size = datasets[0].RasterXSize, datasets[0].RasterYSize
    block_size_x, block_size_y = bands[0].GetBlockSize()

    for y in range(0, y_size, block_size_y):
        rows = min(block_size_y, y_size - y)

        for x in range(0, x_size, block_size_x):
            cols = min(block_size_x, x_size - x)

            arrays = [band.ReadAsArray(x, y, cols, rows) for band in bands]
            stack = np.stack(arrays)
            if nodata is not None:
                stack = np.ma.masked_equal(stack, nodata)
            if method=='min':
                stack_result = stack.min(axis=0)
            elif method=='max':
                stack_result = stack.max(axis=0)
            else:
                stack_result = stack.mean(axis=0)
            bands[0].WriteArray(stack_result, x, y)

    for band in bands:
        band.FlushCache()
        band = None
    for ds in datasets:
        ds = None



def equi7_to_lonlat(roi,
                    pixelsize,
                    input_folder_path,
                    input_file_path,
                    full_output_path,
                    gdal_path=None,
                    pixelsize_deg=None,
                    input_nodata=None,
                    output_nodata=None,
                    aggregation_method='max',
                    num_threads=None,
                    resampling_type='near',
                    data_type=None,
                    compress=True,
                    compress_type="LZW",
                    tiled=True,
                    blocksize=512,
                    bigtiff=True,
                    ):
    """
    Function to generate a GeoTIFF-file in lonlat-projection(epsg:4326) with a
    user-defined extent and resolution from tiled Equi7Grid files.

    Parameters
    ----------
    roi : tuple
        region of interest in [deg] as (lon_min, lat_min, lon_max, lat_max).
    pixelsize : int/float
        pixelsize of Equi7 input data in [m].
    input_folder_path : str
        regex string pattern of input folder structure. Use SUBGRID and TILE as placeholder:
        e.g.: ../Sentinel-1_CSAR/IWGRDH/parameters/datasets/par/B0104/EQUI7_EU010M/E014N023T1/tmensig38/...
        should be: ../Sentinel-1_CSAR/IWGRDH/parameters/datasets/par/B0104/EQUI7_SUBGRID/TILE/tmensig38/...
    input_file_path : regex string pattern of input file structure. Same rule as input_folder_path
        e.g..: M20170527_20171229_TMENSIG38_S1-IWGRDH1VH-_---_B0104_EU010M_E014N023T1.tif
        should be: M*_TMENSIG38_S2-IWGRDH1VH-_---_B0104_SUBGRID_TILE.tif
    full_output_path : str
        folder+file name of the output image.
    gdal_path : str, optional
        Path to the gdal installation, should be set incase not found by the os
    pixelsize_deg : int/float
        The output pixelsize in [deg]. If None the pixelsize will be automatically estimated by gdal.
    input_nodata : int/float, optional
        nodata value to be excluded from processing.
    output_nodata : int/float, optional
        value flagged as nodata in the output image
    aggregation_method : str, optional
        Incase roi overlaps with several Equi7-subgrids, this parameter decides
        how mutliple values per geographic location are handled.
        minimum: 'min', maximum: 'max'(default), or mean: 'mean' are supported
    num_threads : int, optional
        number of threads used for gdal warping.
    resampling_type : str, optional
        gdal resampling method: bilinear, cubic, near (default), ...
    data_type : str, optional
        force output image to have a specific datatype. Byte, UInt16, Int16,
        UInt32, Int32, Float32, Float64, CInt16, CInt32, CFloat32 or CFloat64.
    compress : bool, optional
        Decision if output file is compressed. Defaults to True
    compress_type : str, optional
        Compression method of the output file
    tiled : bool, optional
        Decide if output image is blockwise tiled, else stripped
    blocksize : int, optional
        internal file blocksize of output file. The default is 512.
    bigtiff : bool, optional
        Decision if outputfile is a bigtiff format. Defaults to True

    Returns
    -------
    True, if successfull.
    """

    #dateline check
    if roi[0] > roi[2]:
        raise ValueError('Given bbox crosses antimeridian!')

    #initialise the grid
    e7grid = Equi7Grid(pixelsize)
    #find intersecting tiles with region of interest
    tiles = e7grid.search_tiles_in_roi(bbox=[(roi[0], roi[1]), (roi[2], roi[3])])
    if not tiles:
        raise ValueError(f"No intersecting tiles found for region of interest {roi}.")

    input_filepaths = OrderedDict()
    #get the file for each intersecting tile if it exists and save it grouped per subgrid (=continent-wise)
    for tile in tiles:
        sg = tile.split('_')[0]
        t = tile.split('_')[1]

        #replace wildcards by actual values
        input_path = os.path.join(input_folder_path, input_file_path)
        input_path = input_path.replace('SUBGRID', sg)
        input_path = input_path.replace('TILE', t)

        #find file matching the pattern
        file_path = glob.glob(input_path)

        #if a file was found:
        if file_path:
            if len(file_path)>1:
                print(f'Multiple files found for search path: {input_path}. Only first will be used.')
                file_path = sorted(file_path)[0]

            #create new dict entry incase not yet existing for sg
            if sg not in input_filepaths:
                input_filepaths[sg] = []

            #put filepath into correct list
            input_filepaths[sg].extend(file_path)
        else: print(f'No file found for search path: {input_path}')

    #if not a single input file was found
    if not input_filepaths:
        raise FileNotFoundError("No input files were found for any intersecting tile. Aborting script.")

    #sort groups by their size, since later parameters might be estimated by the first filegroup so make this the one with most input files
    sorted_input_filepaths = sorted(input_filepaths.items(), key=lambda item: len(item[1]), reverse=True)

    #incase there are only files of one subgrid, no zonal overlap hast to be dealt with --> warping result of first subgrid will already be the final result
    onezone_flag = len(sorted_input_filepaths) == 1


    #loop through subgrids
    warping_results = []
    for sg, filepaths in sorted_input_filepaths:
        #build vrt for all files of a subgrid (no complex mosaicing needed, since tiles are co-registered and non-overlapping)
        vrt_path = full_output_path.replace('.tif', f'_{sg}_mosaic.vrt')
        vrt_options = {'xRes': "{}".format(e7grid.core.sampling),
                       'yRes': "{}".format(e7grid.core.sampling)}
        vrt = gdal.BuildVRT(vrt_path,
                            filepaths,
                            options=gdal.BuildVRTOptions(**vrt_options)
                            )
        vrt=None

        # warp the vrt to the target extent/crs
        options = {
            '-t_srs': 'EPSG:4326',
            '-of': 'GTiff',
            '-te': " ".join(map(str, roi)),
            '-r': resampling_type,
            '-multi': ""
            }

        options["-co"] = []
        if compress:
            options["-co"].append("COMPRESS={0}".format(compress_type))
        if input_nodata != None:
            options["-srcnodata"] = input_nodata
        if output_nodata != None:
            options["-dstnodata"] = output_nodata
        if data_type != None:
            options["-ot"] = data_type
        if tiled:
            options["-co"].append("TILED=YES")
            options["-co"].append("BLOCKXSIZE={0}".format(blocksize))
            options["-co"].append("BLOCKYSIZE={0}".format(blocksize))
        else:
            blockxsize = e7grid.core.tile_xsize_m // e7grid.core.sampling
            blockysize = blocksize
            options["-co"].append("TILED=NO")
            options["-co"].append("BLOCKXSIZE={0}".format(blockxsize))
            options["-co"].append("BLOCKYSIZE={0}".format(blockysize))
        if bigtiff:
            options["-co"].append("BIGTIFF=YES")
        if num_threads!= None:
            options["-wo"] = "\"NUM_THREADS={}\"".format(num_threads)


        #grab geotransform from first subgrid iteration if several subgrids involved
        if warping_results:
            ds = open_image(warping_results[0])
            geotransform = ds.geotransform()
            output_bounds = [geotransform[0],
                             geotransform[3]+geotransform[5]*ds.YSize(),
                             geotransform[0]+geotransform[1]*ds.XSize(),
                             geotransform[3]]
            ds.close()
            pixelsize_deg = (geotransform[1], geotransform[5])
            options['-tr'] =  "{} {}".format(pixelsize_deg[0], pixelsize_deg[1])
            options['-te'] = " ".join(map(str, output_bounds))

        #if there was no iteration before, but given px, then use it
        elif pixelsize_deg!=None:
            options['-tr'] =  "{} -{}".format(pixelsize_deg, pixelsize_deg)

        # if only 1 subgrid, warping result is already final file
        if onezone_flag:
            warp_path = full_output_path
        else: warp_path = full_output_path.replace('.tif', f'_{sg}.tif')

        succeed, _ = call_gdal_util('gdalwarp',
                                    src_files=vrt_path,
                                    dst_file=warp_path,
                                    gdal_path=gdal_path,
                                    options=options)

        #remove the temporary vrt file
        os.remove(vrt_path)

        #if only 1 subgrid -- done!
        if onezone_flag:
            return True
        #else also warp other subgrids and merge afterwards
        else:
            warping_results.append(warp_path)

    #Perform chunkwise file aggregation if data is fetched from several Equi7 subgrids
    chunk_calc(aggregation_method, output_nodata, warping_results)
    #Rename first file since aggregation result overwrite first input file
    os.rename(warping_results[0], full_output_path)

    #remove temp files
    for wr in warping_results[1:]:
        os.remove(wr)

    return True

if __name__ == '__main__':
    print('Running Test')
    # TEST: like test1, but with data on R: --> windows test suitable
    #       byte esa-worldcover dataset
    equi7_to_lonlat(roi=(9, 46, 18, 50),
                    pixelsize=20,
                    input_folder_path = r'R:\Datapool\ESA_WorldCover\02_processed\2021_v200\EQUI7_SUBGRID\TILE',
                    input_file_path = 'ESA-WorldCover-2021-V200_SUBGRID_TILE.tif',
                    full_output_path = r'D:\Arbeit\ctemp\EQUI7_to_LATLON2\WorldCover_46-50LON_9-18LAT_bbm.tif',
                    gdal_path = r'C:\Python\Miniconda310_22.11.1\envs\equi7grid\Library\bin',
                    pixelsize_deg=None,
                    input_nodata = 0,
                    output_nodata = 0,
                    aggregation_method=None,
                    resampling_type='near',
                    data_type=None,
                    compress=True,
                    compress_type='LZW',
                    blocksize=512,
                    bigtiff=True,
                    num_threads=1)
    print('Test done!')
