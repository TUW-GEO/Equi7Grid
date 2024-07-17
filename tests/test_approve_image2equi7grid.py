import os

import pytest
from approvaltests.approvals import verify_file
from approvaltests.namer import NamerFactory

from equi7grid.equi7grid import Equi7Grid
from equi7grid.image2equi7grid import TileResampler, image2equi7grid


@pytest.fixture
def out_dir(tmp_path):
    return tmp_path


@pytest.mark.skipif(os.name == 'nt', reason="CI Windows has troubles creating directories")
def test_approve_imag2equi7grid(input_dir, out_dir):
    # begin-snippet: image2equi7grid-example
    input_file = input_dir / "lake_in_russia_lonlat.tif"
    image2equi7grid(Equi7Grid(100), input_file.as_posix(), out_dir.as_posix())

    assert (out_dir / "EQUI7_AS100M/E018N066T6/lake_in_russia_lonlat_AS100M_E018N066T6.tif").exists()
    assert (out_dir / "EQUI7_EU100M/E072N030T6/lake_in_russia_lonlat_EU100M_E072N030T6.tif").exists()
    # end-snippet: image2equi7grid-example

    verify_file((out_dir / "EQUI7_AS100M/E018N066T6/lake_in_russia_lonlat_AS100M_E018N066T6.tif").as_posix(),
                options=NamerFactory.with_parameters("E018N066T6"))
    verify_file((out_dir / "EQUI7_EU100M/E072N030T6/lake_in_russia_lonlat_EU100M_E072N030T6.tif").as_posix(),
                options=NamerFactory.with_parameters("E072N030T6"))


@pytest.mark.skipif(os.name == 'nt', reason="CI Windows has troubles creating directories")
def test_multiprocessing_10m(input_dir, out_dir):
    # begin-snippet: image2equi7grid-example
    input_file = input_dir / "lake_in_russia_lonlat.tif"

    tile_resampler = TileResampler(Equi7Grid(10), input_file.as_posix(), out_dir.as_posix(), parallelism=12)
    results = tile_resampler.resample_tiles()

    print(results)

    assert (out_dir / "EQUI7_AS010M/E021N069T1/lake_in_russia_lonlat_AS010M_E021N069T1.tif").exists()
    assert (out_dir / "EQUI7_EU010M/E073N032T1/lake_in_russia_lonlat_EU010M_E073N032T1.tif").exists()


@pytest.mark.skipif(os.name == 'nt', reason="CI Windows has troubles creating directories")
def test_multiprocessing_20m(input_dir, out_dir):
    # begin-snippet: image2equi7grid-example
    input_file = input_dir / "lake_in_russia_lonlat.tif"

    tile_resampler = TileResampler(Equi7Grid(20), input_file.as_posix(), out_dir.as_posix(), parallelism=12)
    results = tile_resampler.resample_tiles()

    print(results)

    assert (out_dir / "EQUI7_EU020M/E072N030T3/lake_in_russia_lonlat_EU020M_E072N030T3.tif").exists()
    assert (out_dir / "EQUI7_AS020M/E021N069T3/lake_in_russia_lonlat_AS020M_E021N069T3.tif").exists()