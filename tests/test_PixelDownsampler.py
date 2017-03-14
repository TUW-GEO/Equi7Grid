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

Tests for the PixelDownsampler class.

@author: Bernhard Bauer-Marschallinger, bbm@geo.tuwien.ac.at

'''

import numpy.testing as nptest
import numpy as np
from Equi7Grid.Equi7Grid.PixelDownsampler import PixelDownsampler


def test_downsample_via_pixel_indices():

    a_should = np.array([
       [-53.65834427, -52.49286652, -50.99286652, -49.82738876],
       [-46.66547775, -45.5       , -44.        , -42.83452225],
       [-37.66547775, -36.5       , -35.        , -33.83452225],
       [-28.66547775, -27.5       , -26.        , -24.83452225],
       [-19.66547775, -18.5       , -17.        , -15.83452225],
       [-12.67261124, -11.50713348, -10.00713348,  -8.84165573]])

    array = np.arange(9 * 6, dtype=np.float32).reshape((9, 6)) - 56
    array[1, 3] = np.nan

    ds = PixelDownsampler(200, 300, [0, 0, 1200, 1800])
    a = ds.downsample_via_pixel_indices(array)

    nptest.assert_allclose(a_should, a)


if __name__ == '__main__':
    test_downsample_via_pixel_indices()