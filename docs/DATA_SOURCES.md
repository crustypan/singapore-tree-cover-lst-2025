# Sources and data reuse

The current snapshot was accessed on 14 November 2026 (Tokyo time). The coverage output identifies the analyzed files by SHA-256. The MODIS inventory and Landsat screening table record unsigned source URLs and scene identifiers so that the inputs can be traced after time-limited access URLs expire.

## Geography and land cover

| Input | Version and preparation | Licence |
| --- | --- | --- |
| URA Master Plan 2025 Planning Area Boundary (No Sea) | December 2025, 55 planning areas. Download the GeoJSON via the data.gov.sg dataset API. Dissolve the repaired polygons for island-wide calculations. | [Singapore Open Data Licence 1.0](https://data.gov.sg/open-data-licence) |
| ESA WorldCover | 2021 v200, 10 m nominal, EPSG:4326. Download the N00E102 tile and use class 80 for permanent water. The WorldCover non-water mask is the primary land-area denominator. | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| Esri / Impact Observatory Sentinel-2 land use and cover | 2025, version 003, 10 m, UTM 48N. Read the Singapore window from the publisher's 48N GeoTIFF by HTTP byte range and preserve its exact grid and color table. | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |

Contains information from Master Plan 2025 Planning Area Boundary (No Sea), accessed on 14 November 2026 from [URA via data.gov.sg](https://data.gov.sg/datasets/d_2cc750190544007400b2cfd5d7f53209/view) and made available under the [Singapore Open Data Licence version 1.0](https://data.gov.sg/open-data-licence).

© ESA WorldCover project 2021 / Contains modified Copernicus Sentinel data (2021) processed by the ESA WorldCover consortium. Cite Zanaga et al. (2022), ESA WorldCover 10 m 2021 v200. [Dataset DOI](https://doi.org/10.5281/zenodo.7254221) · [Data access and attribution](https://esa-worldcover.org/en/data-access)

Land-cover data from Esri and Impact Observatory, derived from Copernicus Sentinel-2 imagery. Credit the 2025 time series and Karra et al. (2021), *Global land use/land cover with Sentinel-2 and deep learning*. The Esri land-cover classes used in the study are water (1), trees (2), built area (7), clouds (10), and rangeland (11). Rangeland includes open herbaceous cover and some shrub/scrub cover. The snow and ice label occurs in a small number of tropical pixels, illustrating classification error. [Publisher release and licence](https://www.esri.com/about/newsroom/arcnews/latest-land-cover-data-release-shows-more-change-over-time) · [Service legend](https://ic.imagery1.arcgis.com/arcgis/rest/services/Sentinel2_10m_LandCover/ImageServer) · [Paper DOI](https://doi.org/10.1109/IGARSS47720.2021.9553499)

## Temperature and QA

NASA Terra MOD11A1 Collection 6.1 (Wan, Hook, and Hulley, 2021) is a daily nominal 1 km Level 3 product. The pipeline uses daytime `LST_Day_1km` and `QC_Day`, and queries NASA's CMR for all 2025 granules in h28v08. It reads native-grid windows of available cloud-optimized GeoTIFF pairs from Microsoft's public Planetary Computer data mirror. The audit retains dates without available files separately from dates with available data and zero passing pixels. NASA LP DAAC permits data reuse and redistribution. [Dataset DOI](https://doi.org/10.5067/MODIS/MOD11A1.061) · [Band and bit definitions](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD11A1) · [Mirror catalog](https://planetarycomputer.microsoft.com/dataset/modis-11A1-061)

USGS Landsat 8/9 Collection 2 Level-2 Science Product (L2SP) supplies daytime `ST_B10` and `QA_PIXEL` for candidate date screening. The selected-date study also uses `ST_CDIST` and `ST_QA` for sensitivity checks. These encode cloud distance in 0.01 km and product-reported temperature uncertainty in 0.01 K. The source windows for the selected scenes were retrieved on 15 November 2026 (Tokyo time). The effective thermal measurement is coarser than the 30 m delivery grid. The USGS Landsat data can be reused and redistributed. [Dataset DOI](https://doi.org/10.5066/P9OGBGM6) · [USGS temperature documentation](https://www.usgs.gov/landsat-missions/landsat-collection-2-surface-temperature) · [QA definitions](https://www.usgs.gov/landsat-missions/landsat-collection-2-quality-assessment-bands)
