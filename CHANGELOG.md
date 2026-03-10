# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

<!-- insertion marker -->

## [v1.0.0](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v1.0.0) - 2026-03-10

<small>[Compare with v1.0.0](https://github.com/TUW-GEO/Equi7Grid/compare/v1.0.0...v1.0.0)</small>

### Added

- added static HTML image ([f0c994d](https://github.com/TUW-GEO/Equi7Grid/commit/f0c994d24bf7be9d948cadf0a9b8773c83cf5a09) by claxn).
- added first parts of the webviewer documentation ([d2b542a](https://github.com/TUW-GEO/Equi7Grid/commit/d2b542a0e999f1b02df2d99322d7c7b629a7346e) by claxn).

### Fixed

- fix: remove artifact from jinja template ([1a402ce](https://github.com/TUW-GEO/Equi7Grid/commit/1a402cef25726b5dae46fde383f80fb88645fcb3) by npikall).
- fixed relative path ([96459e6](https://github.com/TUW-GEO/Equi7Grid/commit/96459e6acb767b60a965e23e3e23ab7f53825ed0) by claxn).
- fixed bug in docs build;finished layer manager and tile query docs ([5779cf5](https://github.com/TUW-GEO/Equi7Grid/commit/5779cf58ebc89bdff272e1bf4a7ff2cfbe1ff306) by claxn).
- fixed doc tag build; finised coordinate transformer docs ([c21d51e](https://github.com/TUW-GEO/Equi7Grid/commit/c21d51eff746f7e248d62734e24faf66b1d72076) by claxn).
- fixed typo in doc version switcher;finished home view webviewer docs ([076c330](https://github.com/TUW-GEO/Equi7Grid/commit/076c330753d65b0c88395d8bfb9b908f63196069) by claxn).

## [v1.0.0](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v1.0.0) - 2026-02-23

<small>[Compare with v0.2.6](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.6...v1.0.0)</small>

### Added

- added notes for bounding box longitude order ([7fcecc5](https://github.com/TUW-GEO/Equi7Grid/commit/7fcecc5573db1bb4f6024c9ac6c838a202b58fb7) by claudio.schein-navacchi).
- added scrollable outputs ([cc1281c](https://github.com/TUW-GEO/Equi7Grid/commit/cc1281cb0f3c3b0cb3516461aab4bdd45479249f) by claudio.schein-navacchi).
- adds case in test_search_tiles_geog_extent_poles() ([f91b48a](https://github.com/TUW-GEO/Equi7Grid/commit/f91b48a05fa33c5965ca7e6289a4e659cab328f8) by bbm).
- added approvaltests for warp module;fixed accurate polygon boundary size issue ([e89bc79](https://github.com/TUW-GEO/Equi7Grid/commit/e89bc7994eeeab3bbcc77564ad59f3b752452aaf) by claudio.schein-navacchi).
- added lonlat_to_xy on a grid level ([799fca1](https://github.com/TUW-GEO/Equi7Grid/commit/799fca1ea2a8286c6fd787fa6c5e4d60c1de0ba5) by claudio.schein-navacchi).
- added doc files;updated ruff version;updated readme ([acd5ba3](https://github.com/TUW-GEO/Equi7Grid/commit/acd5ba32d84e2ab52c05e667594b727979b1ff3e) by claudio.schein-navacchi).
- added and finished warp module ([a3d8e3d](https://github.com/TUW-GEO/Equi7Grid/commit/a3d8e3d56ffd7b98407eab4596690d54752f7afe) by claudio.schein-navacchi).

### Fixed

- fixed warp tests ([6ed51de](https://github.com/TUW-GEO/Equi7Grid/commit/6ed51deee485d722a5abfa281f723c42419dc9a1) by claudio.schein-navacchi).
- fixed test_search_tiles_geog_extent_antimeridian() ([2f65f87](https://github.com/TUW-GEO/Equi7Grid/commit/2f65f87a69f841a078e04e38f4c2ea9f25035957) by claudio.schein-navacchi).
- fixed doc install and imports;added buffered zones;deleted standard grid creation from JSON ([32ffa9f](https://github.com/TUW-GEO/Equi7Grid/commit/32ffa9fdb4c55cae90d829046999d49bd3913e13) by claudio.schein-navacchi).
- fixed wrong workflow folder ([902e431](https://github.com/TUW-GEO/Equi7Grid/commit/902e431dfef6ead98935816433159566461d4046) by claudio.schein-navacchi).
- fixed covariant tping ([616ce9d](https://github.com/TUW-GEO/Equi7Grid/commit/616ce9da500eb4704e4a38fdb819bde72d4fb555) by claudio.schein-navacchi).
- fix project description link ([f692155](https://github.com/TUW-GEO/Equi7Grid/commit/f69215563801aa05b06114eff4cc3abe4f78964a) by Sebastian Hahn).

### Changed

- change switcher package name ([3d4dd9c](https://github.com/TUW-GEO/Equi7Grid/commit/3d4dd9c69246a6a7c2f3aa1b9bc359b03b29da23) by claudio.schein-navacchi).

### Removed

- removes without replacement test_search_tiles_lon_lat_extent_by_points() ([a61e64b](https://github.com/TUW-GEO/Equi7Grid/commit/a61e64baf6c5bad1f460fa074208d3943aedd924) by bbm).
- removed sampling in full tilenames;aligned license information ([e4aade1](https://github.com/TUW-GEO/Equi7Grid/commit/e4aade1efc70dc0aa735419f45372bbf76209e42) by claudio.schein-navacchi).

## [v0.2.6](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.6) - 2025-06-12

<small>[Compare with v0.2.5](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.5...v0.2.6)</small>

### Added

- add missing package dependency ([0f1364b](https://github.com/TUW-GEO/Equi7Grid/commit/0f1364bcc2a4405310bbf86c20ed448a6276144f) by Sebastian Hahn).

### Fixed

- fix snippet syntax ([4781e75](https://github.com/TUW-GEO/Equi7Grid/commit/4781e758b5932992d2c94b5097d20615a1b77640) by Sebastian Hahn).
- fix typo ([2bb1eb8](https://github.com/TUW-GEO/Equi7Grid/commit/2bb1eb8ebe6da2a2645f955934317dd9e6907730) by Sebastian Hahn).
- fix missing snippet close ([8129eb4](https://github.com/TUW-GEO/Equi7Grid/commit/8129eb492dc7a60594c353f12741946ca3aa8e52) by Sebastian Hahn).
- fix format issue ([9325eac](https://github.com/TUW-GEO/Equi7Grid/commit/9325eacf4e9dfcad7f3db802fcdb0e6199cee1a8) by Sebastian Hahn).

## [v0.2.5](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.5) - 2025-03-14

<small>[Compare with v0.2.4](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.4...v0.2.5)</small>

### Added

- add simple test for 20m data ([a2adbad](https://github.com/TUW-GEO/Equi7Grid/commit/a2adbad27bf44975f8783a60ba1f6bad73d20e52) by Tobias Stachl).
- add simple test ([4b4805a](https://github.com/TUW-GEO/Equi7Grid/commit/4b4805a420176a1e7630ff251d3ac0046a22c821) by Tobias Stachl).
- add parallel processing of equi7grid tiling ([ab134f8](https://github.com/TUW-GEO/Equi7Grid/commit/ab134f8b84666362fb26440dd30bbe04c9b2a0ec) by Tobias Stachl).
- adds news section in README.md ([1589032](https://github.com/TUW-GEO/Equi7Grid/commit/1589032e4db7ce3bdb3ad3efe1400efd84ed01ce) by Bernhard BM).

### Fixed

- fixes \_version.py ([f44f491](https://github.com/TUW-GEO/Equi7Grid/commit/f44f491e45f4bc4fb9c303d6125125413e14e1a6) by bbm).
- fix outfile naming from naming_convention ([c2fc7b4](https://github.com/TUW-GEO/Equi7Grid/commit/c2fc7b4820b28a77c3985da44dcc7e8806ab84b1) by Tobias Stachl).
- fix scipy deprecation warning ([a4ceb11](https://github.com/TUW-GEO/Equi7Grid/commit/a4ceb11fa9097d87d8133798fb7d0daa3004af3e) by Tobias Stachl).

## [v0.2.4](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.4) - 2023-08-16

<small>[Compare with v0.2.3](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.3...v0.2.4)</small>

## [v0.2.3](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.3) - 2023-08-09

<small>[Compare with v0.2.2](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.2...v0.2.3)</small>

## [v0.2.2](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.2) - 2023-07-27

<small>[Compare with v0.2.1](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.1...v0.2.2)</small>

## [v0.2.1](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.1) - 2023-06-21

<small>[Compare with v0.2.0](https://github.com/TUW-GEO/Equi7Grid/compare/v0.2.0...v0.2.1)</small>

### Added

- add updated license: Equi7_V14_License.pdf ([9cb0811](https://github.com/TUW-GEO/Equi7Grid/commit/9cb0811e26636f25de0faf90d8e859a30aeeff2b) by Bernhard BM).

### Removed

- removes deprecated np.int() ([2cdbb15](https://github.com/TUW-GEO/Equi7Grid/commit/2cdbb15742cb5e14365a07ea6951158283911989) by bbm).

## [v0.2.0](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.2.0) - 2022-09-30

<small>[Compare with v0.1.0](https://github.com/TUW-GEO/Equi7Grid/compare/v0.1.0...v0.2.0)</small>

### Fixed

- fix warnings ([01de1b8](https://github.com/TUW-GEO/Equi7Grid/commit/01de1b8319dcb5c3ea41af335dad9c3482323546) by Sebastian Hahn).
- fix filenames ([7507b81](https://github.com/TUW-GEO/Equi7Grid/commit/7507b8194a09616e3ec9b66adea545849db2a03b) by Sebastian Hahn).

## [v0.1.0](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.1.0) - 2021-07-07

<small>[Compare with v0.0.14](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.14...v0.1.0)</small>

### Added

- added geographiclib in dependencies ([c84b69b](https://github.com/TUW-GEO/Equi7Grid/commit/c84b69bca84f19935ae808ac273bd2da29c8345e) by bbm).
- added stripped and tiled GeoTIFF creation to image2equi7grid ([8d1adda](https://github.com/TUW-GEO/Equi7Grid/commit/8d1addab64673ddc2ee048615189e0a62e5db774) by cnavacch).

## [v0.0.14](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.14) - 2021-04-19

<small>[Compare with v0.0.13](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.13...v0.0.14)</small>

### Added

- adds create_distortion_files.py module ([ac89226](https://github.com/TUW-GEO/Equi7Grid/commit/ac892268e5b7a9ea2ee070e0579f8ed36eb56a68) by bbm).
- adds calc_length_distortion functions, updates pytileproj dependency ([4fb6269](https://github.com/TUW-GEO/Equi7Grid/commit/4fb626981d29fa48360c054458610797a6814f86) by bbm).
- adds data_type keyword to image2equi7grid.py ([2e7898c](https://github.com/TUW-GEO/Equi7Grid/commit/2e7898c204a3f2bfa4e46d1373b90c605fe0e791) by bbm).

### Changed

- changelog ([460d161](https://github.com/TUW-GEO/Equi7Grid/commit/460d1615594441902ec30a90971239462bfda664) by bbm).

## [v0.0.13](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.13) - 2021-04-06

<small>[Compare with v0.0.12](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.12...v0.0.13)</small>

## [v0.0.12](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.12) - 2021-02-09

<small>[Compare with v0.0.11](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.11...v0.0.12)</small>

### Changed

- Changed image2equi7grid interface and added geopathfinder file naming classes usage ([10ce4ff](https://github.com/TUW-GEO/Equi7Grid/commit/10ce4ffcf5c8e50e45db20bfe4297b0611c7f25c) by cnavacch).

## [v0.0.11](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.11) - 2020-10-02

<small>[Compare with v0.0.10](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.10...v0.0.11)</small>

### Added

- add short description ([11e51f2](https://github.com/TUW-GEO/Equi7Grid/commit/11e51f2df1919f8139a625c88dee7124ddccef65) by Bernhard Bauer-Marschallinger).
- adds test for center pixel coordinates in ij2xy() ([ba041ba](https://github.com/TUW-GEO/Equi7Grid/commit/ba041ba64fd488e06143b5ea4c04c77b1aff4d4f) by Bernhard Bauer-Marschallinger).

## [v0.0.10](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.10) - 2019-10-30

<small>[Compare with v0.0.9](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.9...v0.0.10)</small>

## [v0.0.9](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.9) - 2019-10-22

<small>[Compare with v0.0.8](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.8...v0.0.9)</small>

### Changed

- changelog ([e29e83c](https://github.com/TUW-GEO/Equi7Grid/commit/e29e83c9d355fa1dba0d3475a458a53b94785beb) by Bernhard Bauer-Marschallinger).

## [v0.0.8](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.8) - 2019-10-10

<small>[Compare with v0.0.7](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.7...v0.0.8)</small>

### Added

- adds create_multipoint_geometry(); new case for split_polygon_by_antimeridian(); big fix in get_geometry_envelope() - to be continued ([382e203](https://github.com/TUW-GEO/Equi7Grid/commit/382e203ec66275b47f620489a3021fb2bcaa8896) by Bernhard Bauer-Marschallinger).
- adds get_tile_bbox_geog/proj(); adds envelope to Tile(); cleaning up search_tiles_over_geometry(); renaming collect_congruent_tiles() ([4abc830](https://github.com/TUW-GEO/Equi7Grid/commit/4abc83063760cff85f979f054d2c07c7df2625c0) by Bernhard Bauer-Marschallinger).

### Fixed

- fix #13 support subfolder with tile folder when resampling. ([67d357c](https://github.com/TUW-GEO/Equi7Grid/commit/67d357c0a6e94bba2054ad01f22833dea4dff1d3) by Senmao Cao).

### Changed

- changes around get_geometry_envelope() ([aa2889c](https://github.com/TUW-GEO/Equi7Grid/commit/aa2889c814093dfc103f0513df2a803dea82e862) by Bernhard Bauer-Marschallinger).

## [v0.0.7](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.7) - 2019-07-11

<small>[Compare with v0.0.6](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.6...v0.0.7)</small>

### Added

- added cut_polygon_by_antimeridian() & co - part 2 ([f8a5632](https://github.com/TUW-GEO/Equi7Grid/commit/f8a5632bd73174c6c5897febd12cb61888e5d9a1) by Bernhard Bauer-Marschallinger).
- added cut_polygon_by_antimeridian() & co - part 1 ([5700a6a](https://github.com/TUW-GEO/Equi7Grid/commit/5700a6aa6b63cd3625d1ceb173de7d8c43939354) by Bernhard Bauer-Marschallinger).
- added test_create_tiles_overlapping_xybbox() ([87b5847](https://github.com/TUW-GEO/Equi7Grid/commit/87b58477df0342c00a559129e630a79d6580a399) by Bernhard Bauer-Marschallinger).
- added functions lonlat2ij_in_tile() & xy2ij_in_tile() ([2a9ee97](https://github.com/TUW-GEO/Equi7Grid/commit/2a9ee974159ceec5762bc637516199c73d2a4de8) by Bernhard Bauer-Marschallinger).

## [v0.0.6](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.6) - 2018-11-23

<small>[Compare with v0.0.5](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.5...v0.0.6)</small>

## [v0.0.5](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.5) - 2018-10-17

<small>[Compare with v0.0.4](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.4...v0.0.5)</small>

### Added

- added get_covering_tiles(); clarified integer divisions, needs new pytileproj=0.0.8 ([d86e7ed](https://github.com/TUW-GEO/Equi7Grid/commit/d86e7ed67df1eb494c7271e58afe41852c39ecac) by Bernhard Bauer-Marschallinger).

### Removed

- removed test for py3.4 ([4bdd61b](https://github.com/TUW-GEO/Equi7Grid/commit/4bdd61b4247e178281bc96c824e61ba0df8be109) by Bernhard Bauer-Marschallinger).

## [v0.0.4](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.4) - 2018-07-30

<small>[Compare with v0.0.3](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.3...v0.0.4)</small>

### Added

- added image2equi7grid.py for heritage of resample() ([151465e](https://github.com/TUW-GEO/Equi7Grid/commit/151465e77ca8e572cade21383eaa4a6980cc0288) by Bernhard Bauer-Marschallinger).
- adds nose to test requirements ([8e155cb](https://github.com/TUW-GEO/Equi7Grid/commit/8e155cb634a405e94b081e87bf3d1fcdc2f5982b) by Bernhard Bauer-Marschallinger).
- adds a test for tile search over kamchatka, and updates other tile search tests ([1ed9630](https://github.com/TUW-GEO/Equi7Grid/commit/1ed9630194aeef35d200f84b207331947324e510) by Bernhard Bauer-Marschallinger).
- add module to make equi7grid.dat ([260202f](https://github.com/TUW-GEO/Equi7Grid/commit/260202f701108480200d2da5498122a6765a4551) by Senmao Cao).
- added test_search_tiles_spitzbergen() and small docu edits ([7b6b127](https://github.com/TUW-GEO/Equi7Grid/commit/7b6b127e8579c308d5e2b42f1ca11f79bf71a2b2) by Bernhard Bauer-Marschallinger).

### Fixed

- fix in documentation ([aea43c5](https://github.com/TUW-GEO/Equi7Grid/commit/aea43c5083e9e7a220fae2e3ada0ae3d7e56f3f6) by Bernhard Bauer-Marschallinger).

## [v0.0.3](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.3) - 2018-01-31

<small>[Compare with v0.0.2](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.2...v0.0.3)</small>

### Added

- add additional test in test_search_tile_500_lon_lat_extent, if all T6 land tiles are returned ([ba66e2d](https://github.com/TUW-GEO/Equi7Grid/commit/ba66e2d5e18b03af8c8bcc0764504c8c958a00f1) by Simon Hochstöger).

### Fixed

- Fix unused imports ([0472a77](https://github.com/TUW-GEO/Equi7Grid/commit/0472a77bfbaee4293ac40cecf44e9d0a680eac61) by Sebastian Hahn).
- Fix test, sort results ([1775a0b](https://github.com/TUW-GEO/Equi7Grid/commit/1775a0bdbcd5afe4bf32335f71bfad0bd4f842a2) by Sebastian Hahn).
- Fix tests ([3371a8a](https://github.com/TUW-GEO/Equi7Grid/commit/3371a8a3dee3ed59966784211f04b45932ce4022) by Sebastian Hahn).

### Changed

- changes due to python3 ([ff2c3ee](https://github.com/TUW-GEO/Equi7Grid/commit/ff2c3ee31beef9e760fdb574b728e2b623d75e4b) by Simon Hochstöger).

### Removed

- removed empty lines and conftest.py ([cdd8c7c](https://github.com/TUW-GEO/Equi7Grid/commit/cdd8c7ca4535d635cc62338ec59fcaaa16334fc4) by bbm).

## [v0.0.2](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.2) - 2017-11-14

<small>[Compare with 0.0.4](https://github.com/TUW-GEO/Equi7Grid/compare/0.0.4...v0.0.2)</small>

### Changed

- changed version of required package pytileproj ([c25f06d](https://github.com/TUW-GEO/Equi7Grid/commit/c25f06d34c803528e579c6b581b698259e9ef2bb) by Bernhard Bauer-Marschallinger).

## [0.0.4](https://github.com/TUW-GEO/Equi7Grid/releases/tag/0.0.4) - 2017-11-06

<small>[Compare with v0.0.1](https://github.com/TUW-GEO/Equi7Grid/compare/v0.0.1...0.0.4)</small>

### Added

- added test_geom_interset() ([2d328eb](https://github.com/TUW-GEO/Equi7Grid/commit/2d328ebf1ebb677010969140642273318480a5b6) by Bernhard Bauer-Marschallinger).

## [v0.0.1](https://github.com/TUW-GEO/Equi7Grid/releases/tag/v0.0.1) - 2017-09-27

<small>[Compare with first commit](https://github.com/TUW-GEO/Equi7Grid/compare/28b01b331ac8a5bb901a407e0b07eab17ef4cb5c...v0.0.1)</small>

### Added

- added conda_environment.yml ([abd34c5](https://github.com/TUW-GEO/Equi7Grid/commit/abd34c5e354075f19bb27c934ad29bbfd075f373) by Bernhard Bauer-Marschallinger).

### Removed

- removed ([55e370b](https://github.com/TUW-GEO/Equi7Grid/commit/55e370b8b51db289435e42f52fc77740461122ef) by Bernhard Bauer-Marschallinger).
