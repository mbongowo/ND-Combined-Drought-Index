ND Combined Drought Index
=========================

Code for:

Mbuh, M., Caparo, A., & Alu, O. (2026). Mapping Drought Patterns in North Dakota, USA, with an
Experimental Landsat-Based Combined Drought Index: An Evaluation of Fourteen Years Between 2000
and 2023. Journal of Earth Observation and Geospatial Applications, 2(2), 112-153.
https://doi.org/10.65372/k588jq80

Combined Drought Index (Section 2.3.6)
  CDI_corrected.js           Google Earth Engine (Code Editor) script for the corrected CDI.
                             Exports the annual CDI at 30 m for the fourteen analyzed years.
                             The normalization limits and weights in it come from
                             cdi_corrected_run.py.
  cdi_corrected_run.py       The same computation through the Earth Engine Python API. Draws the
                             pooled sample (5,000 cells per year), derives the normalization limits
                             (2nd and 98th percentiles) and the principal-component weights, and
                             writes the annual CDI and its layers on the 837.9 m analysis grid.
  export_cdi_30m.py          Queues the 30 m exports to Google Drive from the Python run.

Analyses on the 837.9 m grid (Sections 3.3 to 3.8)
  run_corrected_pipeline.sh  Runs the scripts below in order.
  analyze_cdi_corrected.py   Tables 1 to 3 and Figures 7 to 11.
  analyze.py                 Statewide benchmark against the USDM, nClimDiv PDSI and SPI-3, and
                             state yields (Tables 5, 7 and 8), including the area-weighted check.
  county_benchmark.py        County benchmark (Table 9) and the Landsat-era test.
  variance_decomposition.py  Share of CDI variance fixed in space, common to a year, or residual.
  checks_corrected.py        Checks of Section 2.3.6, Williams' tests and the rankings of Table 6.
  fig6.py                    Figure 6.
  fig16_corrected.py         Figure 16.
  year_coverage.py           Landsat scene counts and June-August composite coverage, 2000-2023
                             (Section 2.3.5).

Other figures and tables
  MSI_moisture_stress_index_GEE.js   MSI maps of Figures 2 to 5 (Section 2.3.5).
  Table4_landcover_P_ET.py           Table 4, precipitation and evapotranspiration by land cover.
  fig1_study_area.py                 Figure 1.
  fix_fig12_13_aspect.py             Proportion correction of Figures 12 and 13.

Requirements: Python 3 with earthengine-api, numpy, pandas, scipy, rasterio, geopandas and
matplotlib; an Earth Engine account. Reference data are public (see the Data Availability
Statement): U.S. Drought Monitor, NOAA NCEI nClimDiv, USDA NASS Quick Stats, U.S. Census
Bureau cartographic boundaries.

Data
  data/annual_reference_series.csv   Annual June-August reference series for the fourteen analyzed
                                     years: USDM DSCI and D0-D2 area, nClimDiv PDSI, SPI-3,
                                     precipitation and temperature, and detrended state yield
                                     anomalies (%) of nine crops. Pass it as VALIDATION_CSV to
                                     run_corrected_pipeline.sh; the CDI columns are added from the
                                     rasters written by cdi_corrected_run.py.

  The Earth Engine scripts read the North Dakota boundary and the NDVI/LST limit images from
  projects/ee-geodrought/assets/. Other users can replace them with their own copies of the
  state boundary (U.S. Census cartographic boundary) and the limit images.

Citation
  Please cite the article above if you use this code.

License
  MIT (see LICENSE).
