import sys, json, numpy as np, pandas as pd, rasterio, geopandas as gpd
from scipy import stats
from pathlib import Path

PROJ = Path(sys.argv[1]); DL = Path(sys.argv[2])
YEARS = [2000,2001,2002,2003,2005,2006,2007,2009,2011,2012,2017,2018,2021,2023]
out = {}

# raster grid
arrs = []
for y in YEARS:
    with rasterio.open(PROJ/"CDI_final_package/preview_rasters"/f"CDI_{y}_utm14_preview.tif") as d:
        a = d.read(1).astype(float); nd = d.nodata
        if nd is not None: a[a == nd] = np.nan
        a[~np.isfinite(a)] = np.nan
        arrs.append(a)
        meta = dict(crs=str(d.crs), res=d.res, width=d.width, height=d.height, bounds=list(d.bounds), nodata=nd)
out["grid"] = meta
stack = np.stack(arrs)
valid_each = np.isfinite(stack).sum(axis=(1,2))
shared = np.all(np.isfinite(stack), axis=0)
cell_km2 = abs(meta["res"][0]*meta["res"][1])/1e6
out["valid_each"] = dict(zip(YEARS, valid_each.tolist()))
out["shared_pixels"] = int(shared.sum()); out["shared_area_km2"] = float(shared.sum()*cell_km2)
out["union_pixels"] = int(np.any(np.isfinite(stack),axis=0).sum())
out["cell_km2"] = cell_km2
vals = stack[:, shared]
q75 = float(np.quantile(vals, 0.75))
rows = []
for i,y in enumerate(YEARS):
    v = vals[i]
    rows.append(dict(year=y, mean=v.mean(), median=np.median(v), std=v.std(), min=v.min(), max=v.max(), pct_q75=100*(v>=q75).mean()))
shared_df = pd.DataFrame(rows); out["shared_q75"] = q75
out["shared_table"] = shared_df.round(4).to_dict(orient="records")
out["shared_pooled_sd"] = float(vals.std())

# climate divisions, equal vs area weights
cd = gpd.read_file(DL/"shp/GIS.OFFICIAL_CLIM_DIVISIONS.shp")
cols = list(cd.columns); out["cd_cols"] = cols
nd = cd[cd["STATE"].astype(str).str.upper().isin(["NORTH DAKOTA","ND","38"])] if "STATE" in cd.columns else None
if nd is None or len(nd)==0:
    for c in cols:
        if cd[c].astype(str).str.contains("North Dakota", case=False).any():
            nd = cd[cd[c].astype(str).str.contains("North Dakota", case=False)]; break
nd = nd.to_crs(5070)
divcol = [c for c in nd.columns if c.upper() in ("CD_NEW","CD_2DIG","DIVISION","CLIMDIV","CD")]
out["divcols"] = divcol
nd["area_km2"] = nd.geometry.area/1e6
nd["div"] = nd[divcol[0]].astype(int) % 100 if divcol else range(1,len(nd)+1)
areas = nd.groupby("div")["area_km2"].sum()
out["div_areas"] = areas.round(0).to_dict()

def read_cd(fn):
    recs = {}
    for line in open(DL/fn):
        code = line[:10]; st, dv, el, yr = code[:2], code[2:4], code[4:6], int(code[6:10])
        if st != "32": continue
        v = [float(x) for x in line[10:].split()]
        recs[(int(dv), yr)] = v
    return recs
w = areas/areas.sum()
cdres = {}
for name, fn in [("pdsi","pdsidv.txt"),("spi3","sp03dv.txt"),("pcp","pcpndv.txt"),("tmp","tmpcdv.txt")]:
    r = read_cd(fn); eq=[]; aw=[]
    for y in YEARS:
        dv = np.array([np.mean(r[(d,y)][5:8]) for d in areas.index])
        if name=="pcp": dv = np.array([np.sum(r[(d,y)][5:8])/3 for d in areas.index])
        eq.append(dv.mean()); aw.append((dv*w.values).sum())
    cdres[name] = dict(eq=eq, aw=aw)
out["cd_series"] = {k:{"eq":np.round(v["eq"],3).tolist(),"aw":np.round(v["aw"],3).tolist()} for k,v in cdres.items()}

V = pd.read_csv(PROJ/"validation/CDI_validation_dataset.csv").set_index("year").loc[YEARS]
cdi = V["cdi_mean"].values; dsci = V["usdm_dsci_jja"].values
out["check_pdsi_eq_vs_dataset"] = float(np.corrcoef(cdres["pdsi"]["eq"], V["ncei_pdsi_jja"])[0,1])

def fisher_ci(r, n=14):
    z = np.arctanh(r); se = 1/np.sqrt(n-3)
    return np.tanh(z-1.96*se), np.tanh(z+1.96*se)
rng = np.random.default_rng(42)
def boot_ci(x, y, B=10000):
    n=len(x); idx = rng.integers(0,n,(B,n)); rs=[]
    for ii in idx:
        xx, yy = x[ii], y[ii]
        if xx.std()==0 or yy.std()==0: continue
        rs.append(np.corrcoef(xx,yy)[0,1])
    return np.percentile(rs,[2.5,97.5])
def row(x, y):
    r, p = stats.pearsonr(x, y); rho, prho = stats.spearmanr(x, y)
    lo, hi = fisher_ci(r); blo, bhi = boot_ci(np.asarray(x), np.asarray(y))
    loo = [stats.pearsonr(np.delete(x,i), np.delete(y,i))[0] for i in range(len(x))]
    return dict(r=r, p=p, lo=lo, hi=hi, blo=blo, bhi=bhi, rho=rho, prho=prho, loo_min=min(loo), loo_max=max(loo),
                loo_min_year=YEARS[int(np.argmin(loo))], loo_max_year=YEARS[int(np.argmax(loo))])
bench = {}
pairs = {
 "DSCI~CDI": (V.usdm_dsci_jja, V.cdi_mean), "DSCI~Q75": (V.usdm_dsci_jja, V.cdi_pct_ge_q75),
 "D0~CDI": (V.usdm_d0_jja_pct, V.cdi_mean), "D1~CDI": (V.usdm_d1_jja_pct, V.cdi_mean), "D2~CDI": (V.usdm_d2_jja_pct, V.cdi_mean),
 "PDSI~CDI": (V.ncei_pdsi_jja, V.cdi_mean), "PDSI~Q75": (V.ncei_pdsi_jja, V.cdi_pct_ge_q75),
 "SPI3~CDI": (V.ncei_spi3_jja, V.cdi_mean), "SPI3~Q75": (V.ncei_spi3_jja, V.cdi_pct_ge_q75),
 "PCP~CDI": (V.ncei_precip_jja_in, V.cdi_mean), "TMP~CDI": (V.ncei_tmean_jja_F, V.cdi_mean),
 "PDSI~DSCI": (V.ncei_pdsi_jja, V.usdm_dsci_jja), "SPI3~DSCI": (V.ncei_spi3_jja, V.usdm_dsci_jja),
 "PCP~DSCI": (V.ncei_precip_jja_in, V.usdm_dsci_jja), "TMP~DSCI": (V.ncei_tmean_jja_F, V.usdm_dsci_jja),
}
for k,(a,b) in pairs.items(): bench[k] = row(a.values.astype(float), b.values.astype(float))
crops = [c for c in V.columns if c.startswith("yieldanom")]
for c in crops:
    bench["Y:"+c+"~CDI"] = row(V[c].values, cdi); bench["Y:"+c+"~DSCI"] = row(V[c].values, dsci)
# area-weighted sensitivity
for nm in ["pdsi","spi3"]:
    for wt in ["eq","aw"]:
        bench[f"CD_{nm}_{wt}~CDI"] = row(np.array(cdres[nm][wt]), cdi)
        bench[f"CD_{nm}_{wt}~DSCI"] = row(np.array(cdres[nm][wt]), dsci)
out["eq_aw_corr"] = {nm: float(np.corrcoef(cdres[nm]["eq"], cdres[nm]["aw"])[0,1]) for nm in cdres}
out["eq_aw_maxdiff"] = {nm: float(np.max(np.abs(np.array(cdres[nm]["eq"])-np.array(cdres[nm]["aw"])))) for nm in cdres}
# influential-year exclusions
excl = {}
for drop in [[2021],[2011],[2009],[2011,2009],[2021,2011,2009]]:
    m = ~np.isin(YEARS, drop)
    excl[str(drop)] = {ref: dict(zip(["r","p"], stats.pearsonr(V[col].values[m], cdi[m])))
                       for ref,col in [("DSCI","usdm_dsci_jja"),("PDSI","ncei_pdsi_jja"),("SPI3","ncei_spi3_jja")]}
out["exclusions"] = excl
# Spearman from ranks with shared-footprint CDI too
out["shared_vs_published_mean_r"] = float(np.corrcoef(shared_df["mean"], cdi)[0,1])
out["shared_bench"] = {ref: row(V[col].values.astype(float), shared_df["mean"].values) for ref,col in [("DSCI","usdm_dsci_jja"),("PDSI","ncei_pdsi_jja"),("SPI3","ncei_spi3_jja")]}
out["bench"] = bench
def conv(o):
    if isinstance(o,(np.floating,)): return round(float(o),4)
    if isinstance(o,(np.integer,)): return int(o)
    if isinstance(o,np.ndarray): return o.tolist()
    return str(o)
print(json.dumps(out, default=conv, indent=1))
