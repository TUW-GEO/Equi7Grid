# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

<!-- insertion marker -->

## [v1.1.1](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v1.1.1) - 2026-04-14

<small>[Compare with v1.1.0](https://github.com/TUW-GEO/Equi7Grid/compare/v1.1.0...v1.1.1)</small>

### Summary

- updated changelog


## [v1.1.0](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v1.1.0) - 2026-04-14

<small>[Compare with v1.0.0](https://github.com/TUW-GEO/Equi7Grid/compare/v1.0.0...v1.1.0)</small>

### Summary

- updated origin of Equi7 tiles to lower-left
- added global/grid-based tile export
- simplified string representation of classes
- documentation updates
- added segmentation to tile export
- fixed bug for functions accepting wrong tilenames
- replaced EPSG zones with internal ones


## [v1.0.0](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v1.0.0) - 2026-03-10

<small>[Compare with v1.0.0](https://github.com/TUW-GEO/Equi7Grid/compare/v1.0.0...v1.0.0)</small>

### Summary

- complete rebuild of package and underlying engine `pytileproj`
- new interface to generate an Equi7Grid object
- support for buffered zone boundaries
- reduced package data to land and projection zones
- complete rebuild of "image2equi7grid" to `warp` module

## [v0.2.6](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.6) - 2025-06-12

<small>[Compare with v0.2.5](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.5...v0.2.6)</small>

### Summary

- update project metadata

## [v0.2.5](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.5) - 2025-03-14

<small>[Compare with v0.2.4](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.4...v0.2.5)</small>

### Summary

- parallel processing for equi7grid tiling
- fixes dirty tag version handling

## [v0.2.4](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.4) - 2023-08-16

<small>[Compare with v0.2.3](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.3...v0.2.4)</small>

### Summary

- adds manifest to package data

## [v0.2.3](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.3) - 2023-08-09

<small>[Compare with v0.2.2](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.2...v0.2.3)</small>

### Summary

- update of installation resources
- update of readme

## [v0.2.2](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.2) - 2023-07-27

<small>[Compare with v0.2.1](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.1...v0.2.2)</small>

### Summary

- update pyproj -> 3.6 (x/y axis order)
- update numpy -> 1.2 (int)
- added /wkt/ files

## [v0.2.1](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.1) - 2023-06-21

<small>[Compare with v0.2.0](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.0...v0.2.1)</small>

### Summary

- removes deprecated np.int()

## [v0.2.0](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.0) - 2022-09-30

<small>[Compare with v0.1.0](https://github.com/TUW-GEO/Equi7Grid/compare/v0.1.0...v0.2.0)</small>

### Summary

- update project files using pyscaffold v3.3
- code formatting using yapf

## [v0.1.0](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.1.0) - 2021-07-07

<small>[Compare with v0.0.14](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.14...v0.1.0)</small>

### Summary

- support for stripped and tiled GeoTIFF blocks

## [v0.0.14](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.14) - 2021-04-19

<small>[Compare with v0.0.13](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.13...v0.0.14)</small>

### Summary

- minor updates in image2equi7
- adds calc_length_distortion functions
- adds create_distortion_files.py module

## [v0.0.13](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.13) - 2021-04-06

<small>[Compare with v0.0.12](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.12...v0.0.13)</small>

### Summary

- full update of all files \_LAND.shp
- revision of all .prj files
- newly generated equi7grid.dat
- full update of coverland lists within equi7grid.dat

## [v0.0.12](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.12) - 2021-02-09

<small>[Compare with v0.0.11](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.11...v0.0.12)</small>

### Summary

- Changed image2equi7grid interface and added geopathfinder file naming classes usage

## [v0.0.11](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.11) - 2020-10-02

<small>[Compare with v0.0.10](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.10...v0.0.11)</small>

### Summary

- add short description
- adds test for center pixel coordinates in ij2xy()

## [v0.0.10](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.10) - 2019-10-30

<small>[Compare with v0.0.9](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.9...v0.0.10)</small>

### Summary

- in image2equi7()...
  - float64 removed
  - same outputfile as input possible
  - scale/offset is now accepted

## [v0.0.9](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.9) - 2019-10-22

<small>[Compare with v0.0.8](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.8...v0.0.9)</small>

### Summary

- bugfix in antimeridian handling
- image2equi7 ingests now .nc with "bands"

## [v0.0.8](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.8) - 2019-10-10

<small>[Compare with v0.0.7](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.7...v0.0.8)</small>

### Summary

- extended search_tiles_in_roi()
- fixes issues with projection strings when using GDAL>=3.0.0

## [v0.0.7](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.7) - 2019-07-11

<small>[Compare with v0.0.6](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.6...v0.0.7)</small>

### Summary

- fixed issues on antimeridian/dateline

## [v0.0.6](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.6) - 2018-11-23

<small>[Compare with v0.0.5](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.5...v0.0.6)</small>

### Summary

- minor update on lonlat2xy functions

## [v0.0.5](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.5) - 2018-10-17

<small>[Compare with v0.0.4](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.4...v0.0.5)</small>

### Summary

- minor update on tilesystem functions

## [v0.0.4](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.4) - 2018-07-30

<small>[Compare with v0.0.3](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.3...v0.0.4)</small>

### Summary

- compatible with Python 3

## [v0.0.3](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.3) - 2018-01-31

<small>[Compare with v0.0.2](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.2...v0.0.3)</small>

### Summary

- changes due to Python 3

## [v0.0.2](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.2) - 2017-11-14

<small>[Compare with 0.0.4](https://github.com/TUW-GEO/Equi7Grid/compare/0.0.4...v0.0.2)</small>

### Summary

- second Draft

## [v0.0.1](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.1) - 2017-09-27

<small>[Compare with first commit](https://github.com/TUW-GEO/Equi7Grid/compare/28b01b331ac8a5bb901a407e0b07eab17ef4cb5c...v0.0.1)</small>

### Summary

- initial Draft
