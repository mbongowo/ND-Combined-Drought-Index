"""Share of CDI variance that is fixed in space (pixel climatology), common to a year, or residual.

Computed on the 259,800 cells valid in all fourteen years (two-way decomposition without replication).
usage: python variance_decomposition.py PROJ_DIR
"""
import sys
from pathlib import Path
import numpy as np
import rasterio

PROJ = Path(sys.argv[1])
YEARS = [2000, 2001, 2002, 2003, 2005, 2006, 2007, 2009, 2011, 2012, 2017, 2018, 2021, 2023]
arrs = []
for y in YEARS:
    with rasterio.open(PROJ / "CDI_final_package/preview_rasters" / f"CDI_{y}_utm14_preview.tif") as d:
        a = d.read(1).astype(float)
        if d.nodata is not None: a[a == d.nodata] = np.nan
        a[~np.isfinite(a) | (np.abs(a) > 1e30)] = np.nan
        arrs.append(a)
stack = np.stack(arrs)
v = stack[:, np.all(np.isfinite(stack), axis=0)]          # years x cells
g = v.mean()
pix = v.mean(axis=0) - g; yr = v.mean(axis=1) - g
res = v - g - pix[None, :] - yr[:, None]
tot = ((v - g) ** 2).sum()
ss_pix = v.shape[0] * (pix ** 2).sum(); ss_yr = v.shape[1] * (yr ** 2).sum(); ss_res = (res ** 2).sum()
print(f"cells {v.shape[1]}, total SS {tot:.1f}")
print(f"pixel (fixed spatial) share {ss_pix/tot:.3f}; year share {ss_yr/tot:.3f}; residual (year-specific pattern) share {ss_res/tot:.3f}")
