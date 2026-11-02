"""Corrected CDI, run through the Earth Engine Python API.

Same computation as CDI_corrected.js:
  - Collection 2 QA_PIXEL mask on bits 1, 3, 4, 5 (dilated cloud, cloud, shadow, snow)
  - June-August monthly composites, continuous VHI, MSI = SWIR1/NIR, WB = PRISM ppt - TerraClimate AET
  - Equation 9 limits: 2nd/98th percentiles of a pooled sample (5000 cells x 14 years, 300 m)
  - PCA weights: |loadings| on PC1 of the correlation matrix of the normalized layers
  - Equation 10: CDI = [a(1 - VHI) + b MSI + g(1 - WB)] / 3

Writes the annual CDI on the 837.9 m NAD83 / UTM 14N analysis grid used by the other
analysis scripts (700 x 426 cells), plus limits.json and weights.json.

usage: python cdi_corrected_run.py OUT_DIR
"""
import json
import sys
from pathlib import Path

import ee
import numpy as np
import rasterio
from rasterio.transform import Affine

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
ee.Initialize()

YEARS = [2000, 2001, 2002, 2003, 2005, 2006, 2007, 2009, 2011, 2012, 2017, 2018, 2021, 2023]
MONTHS = [6, 7, 8]
A = "projects/ee-geodrought/assets/"
shp = ee.FeatureCollection(A + "nd_shape")
maxLST, maxNDVI = ee.Image(A + "maxLST"), ee.Image(A + "maxNDVI")
minLST, minNDVI = ee.Image(A + "minLST"), ee.Image(A + "minNDVI")

# analysis grid (same as CDI_final_package/preview_rasters)
GRID_T = Affine(837.91, 0.0, 108825.29, 0.0, -837.91, 5440059.95)
W, H = 700, 426


def mask_c2(img):
    qa = img.select("QA_PIXEL")
    bits = (1 << 1) | (1 << 3) | (1 << 4) | (1 << 5)
    return img.updateMask(qa.bitwiseAnd(bits).eq(0))


def scale(img):
    opt = img.select("SR_B.").multiply(0.0000275).add(-0.2).clip(shp)
    th = img.select("ST_B.*").multiply(0.00341802).add(149.0).clip(shp)
    return img.addBands(opt, None, True).addBands(th, None, True)


jja = ee.Filter.calendarRange(6, 8, "month")
L5 = (ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(shp).filterDate("2000-06-01", "2011-09-01")
      .select(["SR_B3", "SR_B4", "SR_B5", "ST_B6", "QA_PIXEL"]).filter(jja)
      .map(scale).map(mask_c2).select(["SR_B3", "SR_B4", "SR_B5", "ST_B6"]))
L7 = (ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(shp).filterDate("2012-06-01", "2012-09-01")
      .select(["SR_B3", "SR_B4", "SR_B5", "ST_B6", "QA_PIXEL"]).filter(jja)
      .map(scale).map(mask_c2).select(["SR_B3", "SR_B4", "SR_B5", "ST_B6"]))
L8 = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(shp)
      .select(["SR_B4", "SR_B5", "SR_B6", "ST_B10", "QA_PIXEL"]).filter(jja)
      .map(scale).map(mask_c2)
      .map(lambda i: i.rename(["SR_B3", "SR_B4", "SR_B5", "ST_B6", "QA_PIXEL"]))
      .select(["SR_B3", "SR_B4", "SR_B5", "ST_B6"]))
L875 = L5.merge(L7).merge(L8)

PRISM = ee.ImageCollection("OREGONSTATE/PRISM/ANm").filterDate("2000-01-01", "2024-01-01").select("ppt")
ET = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate("2000-01-01", "2024-01-01").select("aet")


def summer_indices(y):
    imgs = []
    for m in MONTHS:
        d = ee.Date.fromYMD(y, m, 1)
        c = L875.filterDate(d, d.advance(1, "month"))
        imgs.append(c.mean().set("n", c.size()))
    monthly = ee.ImageCollection.fromImages(imgs).filter(ee.Filter.gt("n", 0))

    def idx(img):
        ndvi = img.normalizedDifference(["SR_B4", "SR_B3"])
        lst = img.select("ST_B6")
        msi = img.select("SR_B5").divide(img.select("SR_B4")).rename("MSI")
        vci = ndvi.subtract(minNDVI).divide(maxNDVI.subtract(minNDVI))
        tci = maxLST.subtract(lst).divide(maxLST.subtract(minLST))
        vhi = vci.multiply(0.5).add(tci.multiply(0.5)).rename("VHI")
        return msi.addBands(vhi)

    return monthly.map(idx).mean().clip(shp)


def summer_wb(y):
    imgs = []
    for m in MONTHS:
        d = ee.Date.fromYMD(y, m, 1)
        p = PRISM.filterDate(d, d.advance(1, "month")).sum()
        e = ET.filterDate(d, d.advance(1, "month")).sum().multiply(0.1)
        imgs.append(p.subtract(e).rename("WB"))
    return ee.ImageCollection.fromImages(imgs).mean().clip(shp)


def layers(y):
    return summer_indices(y).addBands(summer_wb(y))


# pooled sample: 5000 random cells per year, one request per year
# (all fourteen years in one request fails with "Too many concurrent aggregations")
smp_file = OUT / "pooled_sample.csv"
if smp_file.exists():
    S = np.loadtxt(smp_file, delimiter=",", skiprows=1)
else:
    parts = []
    for y in YEARS:
        r = (layers(y).sample(region=shp.geometry(), scale=300, numPixels=5000, seed=y, geometries=False)
             .reduceColumns(ee.Reducer.toList(3), ["VHI", "MSI", "WB"]).get("list").getInfo())
        parts.append(np.column_stack([np.full(len(r), y), np.array(r, float)]))
        print("sample", y, len(r), flush=True)
    S = np.vstack(parts)
    np.savetxt(smp_file, S, delimiter=",", header="year,VHI,MSI,WB", comments="", fmt="%.6f")
X_raw = S[:, 1:4]

# Equation 9 limits: 2nd and 98th percentiles of the pooled sample
lim_file = OUT / "limits.json"
if lim_file.exists():
    limits = json.load(open(lim_file))
else:
    limits = {b: np.percentile(X_raw[:, i], [2, 98]).tolist() for i, b in enumerate(["VHI", "MSI", "WB"])}
    limits["n_sample"] = int(len(X_raw))
    json.dump(limits, open(lim_file, "w"), indent=1)
print("limits", limits, flush=True)


def normalize(img, lim):
    return img.subtract(lim[0]).divide(lim[1] - lim[0]).clamp(0, 1)


def norm_layers(y):
    L = layers(y)
    return (normalize(L.select("VHI"), limits["VHI"]).rename("VHI")
            .addBands(normalize(L.select("MSI"), limits["MSI"]).rename("MSI"))
            .addBands(normalize(L.select("WB"), limits["WB"]).rename("WB")))


# PCA weights
w_file = OUT / "weights.json"
if w_file.exists():
    weights = json.load(open(w_file))
else:
    # same pooled cells, normalized as in Equation 9
    X = np.column_stack([np.clip((X_raw[:, i] - limits[b][0]) / (limits[b][1] - limits[b][0]), 0, 1)
                         for i, b in enumerate(["VHI", "MSI", "WB"])])
    R = np.corrcoef(X, rowvar=False)
    val, vec = np.linalg.eigh(R)
    o = np.argsort(val)[::-1]; val, vec = val[o], vec[:, o]
    load = np.abs(vec[:, 0]) * np.sqrt(val[0])
    weights = dict(alpha_VHI=float(load[0]), beta_MSI=float(load[1]), gamma_WB=float(load[2]),
                   pc1_share=float(val[0] / val.sum()), eigenvalues=val.tolist(), n=len(X),
                   corr=R.tolist(), signs_pc1=np.sign(vec[:, 0]).tolist())
    json.dump(weights, open(w_file, "w"), indent=1)
print("weights", weights, flush=True)
a, b, g = weights["alpha_VHI"], weights["beta_MSI"], weights["gamma_WB"]


def cdi_year(y):
    N = norm_layers(y)
    return (ee.Image(1).subtract(N.select("VHI")).multiply(a)
            .add(N.select("MSI").multiply(b))
            .add(ee.Image(1).subtract(N.select("WB")).multiply(g))
            .divide(3).rename("CDI"))


grid = {"dimensions": {"width": W, "height": H},
        "affineTransform": {"scaleX": GRID_T.a, "shearX": 0, "translateX": GRID_T.c,
                            "shearY": 0, "scaleY": GRID_T.e, "translateY": GRID_T.f},
        "crsCode": "EPSG:26914"}


def fetch(img, path):
    if path.exists():
        return
    arr = ee.data.computePixels({"expression": img.toFloat().unmask(-9999), "fileFormat": "NUMPY_NDARRAY", "grid": grid})
    names = arr.dtype.names
    with rasterio.open(path, "w", driver="GTiff", width=W, height=H, count=len(names), dtype="float32",
                       crs="EPSG:26914", transform=GRID_T, nodata=-9999, compress="deflate") as d:
        for i, n in enumerate(names):
            d.write(arr[n].astype("float32"), i + 1)
            d.set_band_description(i + 1, n)


rd = OUT / "CDI_final_package/preview_rasters"; rd.mkdir(parents=True, exist_ok=True)
ld = OUT / "layers"; ld.mkdir(exist_ok=True)
for y in YEARS:
    fetch(cdi_year(y), rd / f"CDI_{y}_utm14_preview.tif")
    fetch(layers(y).addBands(norm_layers(y).rename(["VHI_n", "MSI_n", "WB_n"])), ld / f"layers_{y}.tif")
    print("done", y, flush=True)
