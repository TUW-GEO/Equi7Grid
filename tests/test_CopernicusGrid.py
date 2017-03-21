# Copyright (c) 2016,Vienna University of Technology,
# Department of Geodesy and Geoinformation
# All rights reserved.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL VIENNA UNIVERSITY OF TECHNOLOGY, DEPARTMENT OF
# GEODESY AND GEOINFORMATION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

'''
Tests for the CopernicusGrid class.
'''

import os

from Equi7Grid.Equi7Grid.CopernicusGrid import CopernicusGrid

import numpy as np
import numpy.testing as nptest


def test_ij2xy():
    """
    Test xy to lonlat projection using double numbers.
    """
    cg = CopernicusGrid(1.0/112)
    lon_should = -150.2411
    lat_should = 50.3214
    globaltile = cg.GLOBAL.tilesys.create_tile()
    lon, lat = globaltile.ij2xy(3333, 4444)
    nptest.assert_allclose(lon_should, lon)
    nptest.assert_allclose(lat_should, lat)

def test_xy2ij():
    """
    Test xy to lonlat projection using double numbers.
    """
    cg = CopernicusGrid(1.0/112)
    column_should = 3333
    row_should = 4444
    globaltile = cg.GLOBAL.tilesys.create_tile()
    column, row = globaltile.xy2ij(-150.2411, 50.3214)
    nptest.assert_allclose(column_should, column)
    nptest.assert_allclose(row_should, row)

if __name__ == '__main__':
    test_ij2xy()
    test_xy2ij()





