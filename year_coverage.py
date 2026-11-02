"""Clear-pixel coverage of the June-August Landsat composite for every year 2000-2023.

For each year: number of scenes over North Dakota, and the share of state cells (837.9 m grid)
with a valid MSI value in each month and in the June-August mean. Uses the corrected
Collection 2 mask and collections of cdi_corrected_run.py.

usage: python year_coverage.py CDI_RUN_OUT_DIR OUT_CSV
"""
import runpy
import sys
from pathlib import Path

import ee
import numpy as np
import pandas as pd

g = runpy.run_path(str(Path(__file__).with_name("cdi_corrected_run.py")), run_name="cdi")
L875, shp, grid = g["L875"], g["shp"], g["grid"]
state = ee.Image.constant(1).clip(shp).mask()

def valid_share(img):
    a = ee.data.computePixels({"expression": ee.Image.cat([img.mask().And(state).unmask(0).rename("v"),
                                                           state.unmask(0).rename("s")]).toByte(),
                               "fileFormat": "NUMPY_NDARRAY", "grid": grid})
    return float(a["v"].sum() / a["s"].sum())

rows = []
for y in range(2000, 2024):
    r = dict(year=y)
    for m in (6, 7, 8):
        d = ee.Date.fromYMD(y, m, 1); c = L875.filterDate(d, d.advance(1, "month"))
        r[f"scenes_{m}"] = c.size().getInfo()
        r[f"cover_{m}"] = valid_share(c.select("SR_B5").mean()) if r[f"scenes_{m}"] else 0.0
    c = L875.filterDate(ee.Date.fromYMD(y, 6, 1), ee.Date.fromYMD(y, 9, 1))
    r["cover_jja"] = valid_share(c.select("SR_B5").mean()) if c.size().getInfo() else 0.0
    rows.append(r); print(r, flush=True)
pd.DataFrame(rows).to_csv(sys.argv[2], index=False)
