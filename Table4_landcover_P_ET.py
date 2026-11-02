"""
Table 4. June-August precipitation (P) and actual ET (AET) by land-cover class,
North Dakota, severe-drought years 2000, 2003, 2007, 2017 and 2023.

Inputs (in data_dir)
  cdl/CDL_<year>_38.tif                USDA NASS Cropland Data Layer, North Dakota (CropScape)
  prism/prism_ppt_us_25m_<yyyymm>.tif  PRISM monthly precipitation, 4 km
  terra/aet_<year>.nc                  TerraClimate v1.1 monthly AET (NCSS subset)

Sources (accessed 19 October 2026)
  https://nassgeodata.gmu.edu/axis2/services/CDLService/GetCDLFile?year=<year>&fips=38
  https://services.nacse.org/prism/data/get/us/4km/ppt/<yyyymm>
  http://thredds.northwestknowledge.net:8080/thredds/ncss/TERRACLIMATE_ALL/data/TerraClimate_aet_<year>.nc
    ?var=aet&north=49.2&west=-104.2&east=-96.4&south=45.8&accept=netcdf
    &time_start=<year>-06-01T00:00:00Z&time_end=<year>-08-31T00:00:00Z

Each CDL cell takes the P and AET of the climate cell it falls in (nearest
neighbor); area-weighted means are then computed by CDL class and by group.

python Table4_landcover_P_ET.py <data_dir> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject
from rasterio.windows import Window

DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
YEARS = [2000, 2003, 2007, 2017, 2023]
MONTHS = ["06", "07", "08"]

# CDL class codes grouped into the land-cover classes used in the figures
CROP = ([1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14] + list(range(21, 40)) + list(range(41, 62))
        + list(range(66, 73)) + list(range(74, 78)) + list(range(204, 251)) + [254])
GROUPS = {
    "Cropland": CROP,
    "Grassland": [176],
    "Forest": [63, 141, 142, 143],
    "Shrubland": [64, 152],
    "Wetland": [87, 190, 195],
    "Water": [83, 111],
    "Urban": [82, 121, 122, 123, 124],
    "Barren": [65, 131],
    "Other": [81, 88, 92, 112],  # Clouds/No Data, Nonag/Undefined, Aquaculture, Ice/Snow
}
CODE2GROUP = {c: g for g, codes in GROUPS.items() for c in codes}


def prism_jja(year):
    total, prof = None, None
    for m in MONTHS:
        with rasterio.open(DATA / "prism" / f"prism_ppt_us_25m_{year}{m}.tif") as d:
            a = d.read(1, masked=True).astype("float64")
            prof = (d.transform, d.crs)
        total = a if total is None else total + a
    return total.filled(np.nan), prof


def terra_jja(year):
    ds = xr.open_dataset(DATA / "terra" / f"aet_{year}.nc")  # xarray applies the 0.1 scale factor
    aet = ds["aet"].sel(time=ds.time.dt.month.isin([6, 7, 8]))
    assert aet.sizes["time"] == 3, f"expected 3 months of AET for {year}"
    a = aet.sum("time", skipna=False).transpose("lat", "lon")
    lat, lon = a["lat"].values, a["lon"].values
    arr = a.values.astype("float64")
    if lat[0] < lat[-1]:  # rows north to south
        arr, lat = arr[::-1], lat[::-1]
    res_x, res_y = abs(lon[1] - lon[0]), abs(lat[1] - lat[0])
    transform = from_origin(lon[0] - res_x / 2, lat[0] + res_y / 2, res_x, res_y)
    return arr, (transform, "EPSG:4326")


def process(year):
    p, (p_tr, p_crs) = prism_jja(year)
    e, (e_tr, e_crs) = terra_jja(year)
    n = 256
    sums = {k: np.zeros(n) for k in ("cells", "cells_valid", "p", "e")}
    with rasterio.open(DATA / "cdl" / f"CDL_{year}_38.tif") as d:
        cell_area = abs(d.transform.a * d.transform.e)  # m2, Albers equal area (30 m; 56 m in 2007)
        assert d.crs.to_epsg() == 5070 or "Albers" in d.crs.to_wkt(), d.crs
        for row0 in range(0, d.height, 2048):
            win = Window(0, row0, d.width, min(2048, d.height - row0))
            lc = d.read(1, window=win)
            tr = d.window_transform(win)
            pw = np.full(lc.shape, np.nan); ew = np.full(lc.shape, np.nan)
            reproject(p, pw, src_transform=p_tr, src_crs=p_crs, dst_transform=tr, dst_crs=d.crs,
                      resampling=Resampling.nearest, src_nodata=np.nan, dst_nodata=np.nan)
            reproject(e, ew, src_transform=e_tr, src_crs=e_crs, dst_transform=tr, dst_crs=d.crs,
                      resampling=Resampling.nearest, src_nodata=np.nan, dst_nodata=np.nan)
            inside = lc > 0  # 0 = outside the state
            ok = inside & np.isfinite(pw) & np.isfinite(ew)
            sums["cells"] += np.bincount(lc[inside], minlength=n)
            sums["cells_valid"] += np.bincount(lc[ok], minlength=n)
            sums["p"] += np.bincount(lc[ok], weights=pw[ok], minlength=n)
            sums["e"] += np.bincount(lc[ok], weights=ew[ok], minlength=n)
    rows = []
    for code in np.nonzero(sums["cells"])[0]:
        v = sums["cells_valid"][code]
        rows.append(dict(year=year, cdl_code=int(code), group=CODE2GROUP.get(int(code), "UNASSIGNED"),
                         area_km2=sums["cells"][code] * cell_area / 1e6,
                         valid_area_km2=v * cell_area / 1e6,
                         p_sum=sums["p"][code], e_sum=sums["e"][code], cells_valid=v))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.concat([process(y) for y in YEARS], ignore_index=True)
    df["P_mm"] = df.p_sum / df.cells_valid
    df["AET_mm"] = df.e_sum / df.cells_valid
    df["P_minus_AET_mm"] = df.P_mm - df.AET_mm
    df.drop(columns=["p_sum", "e_sum"]).to_csv(OUT / "Table4_by_CDL_code.csv", index=False)
    g = df.groupby(["year", "group"]).agg(area_km2=("area_km2", "sum"), valid_area_km2=("valid_area_km2", "sum"),
                                          p_sum=("p_sum", "sum"), e_sum=("e_sum", "sum"),
                                          cells_valid=("cells_valid", "sum")).reset_index()
    g["P_mm"] = g.p_sum / g.cells_valid
    g["AET_mm"] = g.e_sum / g.cells_valid
    g["P_minus_AET_mm"] = g.P_mm - g.AET_mm
    tot = g.groupby("year").area_km2.transform("sum")
    g["pct_of_state"] = 100 * g.area_km2 / tot
    g = g.drop(columns=["p_sum", "e_sum", "cells_valid"])
    g.to_csv(OUT / "Table4_by_group.csv", index=False)
    print(g.round(1).to_string(index=False))
    print("unassigned codes:", sorted(df.loc[df.group == "UNASSIGNED", "cdl_code"].unique()))


if __name__ == "__main__":
    main()
