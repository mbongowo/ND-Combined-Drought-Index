"""County-level benchmark of the CDI.

County mean CDI (from the 837.9 m analysis rasters) is compared, for the fourteen analyzed
years, with detrended NASS county yields and with the county June-August USDM DSCI.
Yield anomalies are percentage departures from a linear trend fitted to each county's
1990-2023 series (same detrending as the state yields in Section 2.3.7).

Three views of each relationship:
  pooled  - all county-years together (mixes between-year and between-county variation)
  within-year (spatial) - correlation across counties in each year; median over the 14 years
  within-county (temporal) - county means removed first, then pooled

Also tests whether the CDI separates by Landsat era (L5/7 2000-2012 vs L8 2017-2023).

usage: python county_benchmark.py PROJ_DIR DATA_DIR OUT_DIR
"""
import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from scipy import stats

PROJ = Path(sys.argv[1]); DATA = Path(sys.argv[2]); OUT = Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
YEARS = [2000, 2001, 2002, 2003, 2005, 2006, 2007, 2009, 2011, 2012, 2017, 2018, 2021, 2023]
CROPS = {
    "Spring wheat": "WHEAT, SPRING, (EXCL DURUM) - YIELD, MEASURED IN BU / ACRE",
    "Durum wheat": "WHEAT, SPRING, DURUM - YIELD, MEASURED IN BU / ACRE",
    "Barley": "BARLEY - YIELD, MEASURED IN BU / ACRE",
    "Corn (grain)": "CORN, GRAIN - YIELD, MEASURED IN BU / ACRE",
    "Soybeans": "SOYBEANS - YIELD, MEASURED IN BU / ACRE",
    "Alfalfa hay": "HAY, ALFALFA - YIELD, MEASURED IN TONS / ACRE",
    "Other hay": "HAY, (EXCL ALFALFA) - YIELD, MEASURED IN TONS / ACRE",
    "Sunflower": "SUNFLOWER, OIL TYPE - YIELD, MEASURED IN LB / ACRE",
    "Canola": "CANOLA - YIELD, MEASURED IN LB / ACRE",
}

# county mean CDI
cty = gpd.read_file(DATA / "counties/cb_2023_us_county_500k.shp")
cty = cty[cty.STATEFP == "38"].to_crs(26914).reset_index(drop=True)
cty["fips"] = (cty.STATEFP + cty.COUNTYFP).astype(int)
rows = []
for y in YEARS:
    with rasterio.open(PROJ / "CDI_final_package/preview_rasters" / f"CDI_{y}_utm14_preview.tif") as d:
        a = d.read(1).astype(float)
        if d.nodata is not None: a[a == d.nodata] = np.nan
        a[~np.isfinite(a) | (a < -1e30) | (a > 1e30)] = np.nan
        zones = rasterize(((g, i + 1) for i, g in enumerate(cty.geometry)), out_shape=a.shape,
                          transform=d.transform, fill=0, dtype="int32")
    for i, f in enumerate(cty.fips):
        v = a[zones == i + 1]; v = v[np.isfinite(v)]
        rows.append(dict(fips=f, year=y, cdi=v.mean() if v.size else np.nan, ncell=v.size))
cdi = pd.DataFrame(rows)

# county JJA DSCI
u = pd.read_csv(DATA / "usdm_county_dsci.csv")
u["date"] = pd.to_datetime(u.MapDate.astype(str)); u = u[u.date.dt.month.isin([6, 7, 8])]
u["year"] = u.date.dt.year
dsci = u.groupby(["FIPS", "year"]).DSCI.mean().rename("dsci").reset_index().rename(columns={"FIPS": "fips"})

# county yield anomalies
q = pd.read_csv(DATA / "nd_yields.tsv", sep="\t", dtype=str,
                usecols=["SOURCE_DESC", "SHORT_DESC", "AGG_LEVEL_DESC", "COUNTY_ANSI", "YEAR",
                         "REFERENCE_PERIOD_DESC", "DOMAIN_DESC", "VALUE"])
q = q[(q.SOURCE_DESC == "SURVEY") & (q.AGG_LEVEL_DESC == "COUNTY") & (q.REFERENCE_PERIOD_DESC == "YEAR")
      & (q.DOMAIN_DESC == "TOTAL") & q.COUNTY_ANSI.notna() & (q.COUNTY_ANSI != "")]
q["value"] = pd.to_numeric(q.VALUE.str.replace(",", ""), errors="coerce")
q["year"] = q.YEAR.astype(int); q["fips"] = 38000 + q.COUNTY_ANSI.astype(int)
q = q[q.year.between(1990, 2023) & q.value.notna()]
anoms = []
for crop, sd in CROPS.items():
    c = q[q.SHORT_DESC == sd]
    for f, g in c.groupby("fips"):
        if g.year.nunique() < 15: continue
        b = np.polyfit(g.year, g.value, 1); fit = np.polyval(b, g.year)
        anoms.append(pd.DataFrame(dict(crop=crop, fips=f, year=g.year, anom=100 * (g.value - fit) / fit)))
yld = pd.concat(anoms)

panel = cdi.merge(dsci, on=["fips", "year"], how="left")
panel.to_csv(OUT / "county_cdi_dsci.csv", index=False)

def summarize(df, x, yv):
    df = df.dropna(subset=[x, yv])
    r_pool, p_pool = stats.pearsonr(df[x], df[yv])
    per_year = [stats.pearsonr(g[x], g[yv])[0] for _, g in df.groupby("year")
                if len(g) >= 10 and g[x].std() > 0 and g[yv].std() > 0]
    dm = df.assign(**{x: df[x] - df.groupby("fips")[x].transform("mean"),
                      yv: df[yv] - df.groupby("fips")[yv].transform("mean")})
    r_wc, p_wc = stats.pearsonr(dm[x], dm[yv])
    # county-years share a year, so p-values above are optimistic; resample whole years instead
    yrs = df.year.unique(); rng = np.random.default_rng(1); boot = []
    for _ in range(2000):
        bd = pd.concat([df[df.year == yy] for yy in rng.choice(yrs, len(yrs))])
        boot.append(np.corrcoef(bd[x], bd[yv])[0, 1])
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    return dict(n=len(df), ci_lo=lo, ci_hi=hi, counties=df.fips.nunique(), r_pooled=r_pool, p_pooled=p_pool,
                r_within_year_median=float(np.median(per_year)),
                r_within_year_min=float(np.min(per_year)), r_within_year_max=float(np.max(per_year)),
                years_neg=int(np.sum(np.array(per_year) < 0)), years=len(per_year),
                r_within_county=r_wc, p_within_county=p_wc)

res = []
for crop in CROPS:
    d = panel.merge(yld[yld.crop == crop], on=["fips", "year"])
    res.append(dict(crop=crop, ref="CDI", **summarize(d, "cdi", "anom")))
    res.append(dict(crop=crop, ref="USDM DSCI", **summarize(d.dropna(subset=["cdi"]), "dsci", "anom")))
tab = pd.DataFrame(res); tab.to_csv(OUT / "county_yield_benchmark.csv", index=False)
cd = summarize(panel, "cdi", "dsci")

# Landsat era test (state level)
V = pd.read_csv(PROJ / "validation/CDI_validation_dataset.csv").set_index("year").loc[YEARS]
era = np.array([0 if y <= 2012 else 1 for y in YEARS])
l57 = V[era == 0]; l8 = V[era == 1]
t_era = stats.ttest_ind(l57.cdi_mean, l8.cdi_mean, equal_var=False)
r57 = stats.pearsonr(l57.cdi_mean, l57.usdm_dsci_jja)
cdi_dm = V.cdi_mean - pd.Series(np.where(era == 0, l57.cdi_mean.mean(), l8.cdi_mean.mean()), index=V.index)
r_dm = stats.pearsonr(cdi_dm, V.usdm_dsci_jja)
r_dm_pdsi = stats.pearsonr(cdi_dm, V.ncei_pdsi_jja)
era_out = dict(mean_l57=l57.cdi_mean.mean(), mean_l8=l8.cdi_mean.mean(), t=t_era.statistic, p=t_era.pvalue,
               r_dsci_l57_only=r57[0], p_dsci_l57_only=r57[1],
               r_dsci_era_demeaned=r_dm[0], p_dsci_era_demeaned=r_dm[1],
               r_pdsi_era_demeaned=r_dm_pdsi[0], p_pdsi_era_demeaned=r_dm_pdsi[1])

json.dump(dict(cdi_vs_county_dsci=cd, era=era_out), open(OUT / "county_benchmark_summary.json", "w"),
          indent=1, default=float)
pd.set_option("display.width", 200)
print(tab.round(3).to_string())
print("CDI vs county DSCI:", {k: round(v, 3) if isinstance(v, float) else v for k, v in cd.items()})
print("Era:", {k: round(v, 4) for k, v in era_out.items()})
