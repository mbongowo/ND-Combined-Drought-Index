"""Queue the 30 m corrected-CDI exports to Google Drive (folder MSI), one task per analyzed year.

Reuses the computation of cdi_corrected_run.py (its cached limits and weights), so the 30 m
rasters match the 837.9 m analysis rasters. Same export settings as CDI_corrected.js;
masked cells are written as -9999 and flagged as NoData.

usage: python export_cdi_30m.py OUT_DIR   (OUT_DIR = the cdi_corrected_run.py output folder)
"""
import runpy
import sys

import ee

g = runpy.run_path(str(__import__("pathlib").Path(__file__).with_name("cdi_corrected_run.py")), run_name="cdi")
for y in g["YEARS"]:
    task = ee.batch.Export.image.toDrive(
        image=g["cdi_year"](y).toFloat().unmask(-9999), description=f"CDI_corrected_30m_{y}", fileNamePrefix=f"CDI_corrected_30m_{y}",
        folder="MSI", region=g["shp"].geometry(), scale=30, crs="EPSG:2266", maxPixels=1e13,
        formatOptions={"cloudOptimized": True, "noData": -9999})
    task.start()
    print(y, task.id, flush=True)
