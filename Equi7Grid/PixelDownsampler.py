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
Created on March 14, 2017

Code for downsampling an array-based image to a smaller array-based image.

Is fast through intelligent pixel averaging.
Gives physical meaningful results through proper filtering.

@author: Bernhard Bauer-Marschallinger, bbm@geo.tuwien.ac.at

'''


import fractions
import math
import copy

import numpy as np
from scipy import ndimage
from numba import jit
from astropy.convolution import convolve
from kernel import Kernel


class PixelDownsampler(object):

    def __init__(self, spacing_fine, spacing_coarse, bbox):

        if spacing_fine >= spacing_coarse:
            raise ValueError('"res_c (={}m) must be larger than '
                             'source grid resolution (={}m!)'.format(spacing_fine, spacing_coarse))

        self.spacing_fine = spacing_fine
        self.spacing_coarse = spacing_coarse
        self.target_bbox = bbox

        pixelmap_fine, pixelmap_coarse, pixel_count, n_cum_x, n_cum_y= self.translate_pixelmaps()
        self.pixelmap_fine = pixelmap_fine
        self.pixelmap_coarse = pixelmap_coarse
        self.pixel_count = pixel_count
        self.n_cum_x = n_cum_x
        self.n_cum_y = n_cum_y

        pass

    def calc_pixel_index_pattern(self):
        '''
        finds best  representation of overlapping pixels lattices through optimal rounding.
        is periodic and aims for most symmetry-


        :return:
        '''

        # get fraction of the two resolutions
        frac = fractions.Fraction('{}/{}'.format(self.spacing_fine, self.spacing_coarse))
        numerator = frac._numerator
        denominator = frac._denominator

        # pattern showing digitized relation of pixel sizes.
        pattern = list()

        # algorithm finding the integers of the pixel lattice relation
        n_items = 0
        excess_sum = 0.0
        remainder = denominator
        for i in range(numerator):
            r = remainder / float(numerator - n_items)
            item = math.ceil(r)
            excess = item - r
            excess_sum += excess
            if excess_sum >= 0.5:
                item -= 1
                excess_sum -= 1.0
            n_items += 1
            remainder -= item
            pattern.append(int(item))

        return pattern


    def translate_pixelmaps(self):

        # resolution of grid
        res_f = self.spacing_fine
        # target resolution
        res_c = self.spacing_coarse

        ratio = 1.0 * self.spacing_coarse / self.spacing_fine

        pattern_f = self.calc_pixel_index_pattern()

        pattern_length_f = sum(pattern_f)
        pattern_length_c = len(pattern_f)

        kgV = res_c * pattern_length_c
        needed_bbox = copy.copy(self.target_bbox)
        needed_bbox[0] -= kgV
        needed_bbox[1] -= kgV
        needed_bbox[2] += kgV
        needed_bbox[3] += kgV
        self.needed_bbox = needed_bbox

        xsize_m = needed_bbox[2] - needed_bbox[0]
        ysize_m = needed_bbox[3] - needed_bbox[1]
        if (xsize_m % (kgV) != 0) or (ysize_m % (kgV) != 0):
            raise ValueError('"bbox" must have width and height '
                             'dividable by {}!'.format(kgV))

        xsize_f = int(xsize_m / float(res_f))
        ysize_f = int(ysize_m / float(res_f))
        xsize_c = int(xsize_m / float(res_c))
        ysize_c = int(ysize_m / float(res_c))

        self.needed_shape = (ysize_f, xsize_f)

        # create template
        pattern_tmpl = list()
        for i in range(len(pattern_f)):
            pattern_tmpl.extend([i] * pattern_f[i])
        # kx, ky:  number of patterns that bbox spans in x and y direction
        kx = int(xsize_f) / pattern_length_f
        ky = int(ysize_f) / pattern_length_f
        idx = np.tile(pattern_tmpl, kx)
        idy = np.tile(pattern_tmpl, ky)
        idx += np.repeat(np.arange(kx) * len(pattern_f), pattern_length_f)
        idy += np.repeat(np.arange(ky) * len(pattern_f), pattern_length_f)

        # create index array
        # TODO: transpose the index? flip the array upside-down?
        pixelmap_fine = np.zeros((ysize_f, xsize_f), dtype=np.uint32)
        for i, v in enumerate(idy):
            pixelmap_fine[i, :] = idx + v * xsize_c

        # gets vectors holding the fine pixel per coarse pixel
        n_pixels_x = (np.unique(idx, return_counts=True)[1]).astype(np.uint16)
        n_pixels_y = (np.unique(idy, return_counts=True)[1]).astype(np.uint16)
        
        cum_pixels_x = np.cumsum(n_pixels_x).astype(np.uint16)
        cum_pixels_y = np.cumsum(n_pixels_y).astype(np.uint16)
        # calculate the  number of pixels (in A tile) overlapping
        # within each pixel of C tile
        pix_count = np.outer(n_pixels_y, n_pixels_x)

        # retrieve the indexes of the course pixels
        if ratio.is_integer():
            pixelmap_coarse = (pixelmap_fine[0::int(ratio), 0::int(ratio)]).copy()
        else:
            # if ratio is float: np.unique is too slow for T1->T6
            # therefore changes in values of pixelmap_finex are detected by subtraction pixelmap_finex with pixelmap_finex shifted by
            # one index. With this resulting mask pixelmap_coarsex can be gained.
            temp_idx = np.append(pixelmap_fine[1:, 0] - pixelmap_fine[:-1, 0], 1)
            temp_idy = np.append(pixelmap_fine[0, 1:] - pixelmap_fine[0, :-1], 1)
            temp_id_2d = np.outer(temp_idx, temp_idy)
            pixelmap_coarse = pixelmap_fine[temp_id_2d.astype(bool)]
            pixelmap_coarse = pixelmap_coarse.reshape(int(len(temp_idx)/ratio), int(len(temp_idy)/ratio))

        
        return pixelmap_fine, pixelmap_coarse, pix_count, cum_pixels_x, cum_pixels_y


    def downsample_via_pixel_indices(self, array, function=None):
        """ downsample with pixel averaging and consecutive filtering.

        This method perform the masking, averaging, and resampling
        (and consecutive gaussian filtering)

        """
        needed_shape =self.pixelmap_fine.shape
        if array.shape != needed_shape:
            raise ValueError('Input "array" must have shape={}!'.format(needed_shape))


        no_data_value = -9999.0
        datatype = None

        pixelmap_fine = self.pixelmap_fine
        pixelmap_coarse = self.pixelmap_coarse

        # minimum number of valid fine pixels in coarse pixels
        limit = (self.pixel_count / 100.0 * 40.0).astype(np.uint16)

        # apply function to data array
        if function is not None:
            array = function(array)

        # use max of numpy uint32 as NaN value for index array
        nan_uint32 = np.iinfo(np.uint32).max
        # create a mask of NaN values in LUT index array
        #mask = np.logical_or(np.isnan(array), (array == no_data_value))
        mask = np.isnan(array)
        # applying mask
        pixelmap_fine[mask] = nan_uint32

        # invert and convert the mask
        mask_num = ((~mask) * 1).astype(np.int8)

        # count of valid fine pixels in coarse pixel
        masked_pixels = self.fast_mask_counting(mask_num)

        # only coarse pixels with enough valid fine pixels
        id_enough = masked_pixels <= limit
        ide = pixelmap_coarse[id_enough]

        # resample parameter (do the pixel averaging) to target grid
        array_ds = np.ndarray.ravel(np.full(pixelmap_coarse.shape, np.nan, np.float32))
        array_ds[ide] = ndimage.mean(array, labels=pixelmap_fine, index=ide)

        # encode
        array_en = array_ds.astype(np.float32).reshape(pixelmap_coarse.shape)
        invidx = array_en == no_data_value
        array_en[invidx] = np.nan

        # apply Gaussian filter
        kernel = Kernel(size=3)
        kernel.gauss_weightratio_kernel()
        gauss_kernel = kernel.array

        '''
        gauss_kernel = np.array([(0.04973561603529225, 0.12354360159590561, 0.04973561603529225),
                                 (0.12354360159590561, 0.30688312947520852, 0.12354360159590561),
                                 (0.04973561603529225, 0.12354360159590561, 0.04973561603529225)])
        '''
        array_conv = convolve(array_en, gauss_kernel, boundary='extend')
        array_conv[invidx] = no_data_value
        array_ds[:] = array_conv.flatten()[:]

        array_out = array_ds.reshape(pixelmap_coarse.shape)
        array_out[array_out == np.nan] = no_data_value
        array_out = array_out.astype(datatype)

        return array_out


    def fast_mask_counting(self, mask):
        """
        counting the number of masked fine pixels in individual coarse pixels.
        uses jit from numba for speeding up the loops.

        Parameters
        ----------
        mask : numpy array
            logical array: "1" at index of valid fine pixels,
                           "0" at index of non-valid fine pixels
        course_shape : tuple
            (dim_x, dim_y) shape of the target coarse pixel array
        pattern : numpy array
            array, telling for each coarse-scale position the highest fine pixel index

        Returns
        -------
        result: numpy array
            array holding the counts of non-valid (masked) fine pixels per coarse pixel
        """

        # type check
        course_shape = self.pixelmap_coarse.shape
        pattern_x = self.n_cum_x
        pattern_y = self.n_cum_y

        error_msg = 'fast_mask_counting: input is not valid!'
        types = [np.ndarray, np.array, np.memmap, np.flatiter]
        if (type(mask) not in types) or (type(pattern_x) not in types) or (type(pattern_y) not in types) or \
           (type(course_shape) != tuple) or (mask.dtype != 'int8') or (pattern_x.dtype != 'uint16') or \
                (pattern_y.dtype != 'uint16'):
            raise TypeError(error_msg)

        # fast computation using just in time (jit)
        @jit()
        def run_fast_mask_counting(mask, course_shape, pattern_x, pattern_y):
            m, n = mask.shape
            result = np.zeros(course_shape, dtype=np.int32)
            x = 0
            for i in range(m):
                if i == pattern_y[x]:
                    x += 1
                y = 0
                for j in range(n):
                    if j == pattern_x[y]:
                        y += 1
                    if mask[i, j] != 1:
                        result[x, y] += 1
            return result

        result = run_fast_mask_counting(mask, course_shape, pattern_x, pattern_y)
        return result


def dummy_function(array):

    ind_exclude = (array <= -2000) | (array >= -500)
    array[ind_exclude] = -9999
    array_en = array * 0.01
    array_en[array == -9999] = -9999

    return array