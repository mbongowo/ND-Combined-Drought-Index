"""Checks and summary numbers for the corrected CDI (Section 2.3.6, 3.8, 4.1, 4.3).

usage: python checks_corrected.py ROOT
ROOT holds CDI_final_package/preview_rasters, layers/, analysis/CDI_summary_statistics.csv,
validation/CDI_validation_dataset.csv, weights.json, limits.json.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy import stats

ROOT = Path(sys.argv[1])
YEARS = [2000, 2001, 2002, 2003, 2005, 2006, 2007, 2009, 2011, 2012, 2017, 2018, 2021, 2023]
CLASS = {2002: "extreme", 2006: "extreme", 2012: "extreme", 2021: "extreme",
         2000: "severe", 2003: "severe", 2007: "severe", 2017: "severe", 2023: "severe",
         2001: "moderate", 2005: "moderate", 2009: "moderate", 2011: "moderate", 2018: "moderate"}
MID = {"extreme": 10, "severe": 25, "moderate": 35}
out = {}

w = json.load(open(ROOT / "weights.json")); lim = json.load(open(ROOT / "limits.json"))
out["weights"] = w; out["limits"] = lim
out["upper_bound"] = (w["alpha_VHI"] + w["beta_MSI"] + w["gamma_WB"]) / 3

stack = []
for y in YEARS:
    with rasterio.open(ROOT / "CDI_final_package/preview_rasters" / f"CDI_{y}_utm14_preview.tif") as d:
        a = d.read(1).astype(float); a[a == d.nodata] = np.nan; stack.append(a)
stack = np.stack(stack)
out["cell_min"] = float(np.nanmin(stack)); out["cell_max"] = float(np.nanmax(stack))
out["valid_each"] = {y: int(np.isfinite(stack[i]).sum()) for i, y in enumerate(YEARS)}
shared = np.all(np.isfinite(stack), axis=0)
out["shared_cells"] = int(shared.sum())

S = pd.read_csv(ROOT / "analysis/CDI_summary_statistics.csv").set_index("year").loc[YEARS]
V = pd.read_csv(ROOT / "validation/CDI_validation_dataset.csv").set_index("year").loc[YEARS]
mean = S.mean_cdi.values
out["between_year_sd"] = float(np.std(mean, ddof=1))
out["within_sd_range"] = [float(S.std_cdi.min()), float(S.std_cdi.max())]
out["cv_max"] = float(S.max_cdi.std(ddof=1) / S.max_cdi.mean())
out["cv_p90"] = float(S.p90_cdi.std(ddof=1) / S.p90_cdi.mean())
out["grand_mean"] = float(mean.mean())

# shared-footprint vs own-cells
sm = np.array([np.mean(stack[i][shared]) for i in range(len(YEARS))])
out["shared_mean_maxdiff"] = float(np.max(np.abs(sm - mean)))
out["shared_vs_own_r"] = float(np.corrcoef(sm, mean)[0, 1])

# sign check: CDI against the VHI class midpoint, and class means
mids = np.array([MID[CLASS[y]] for y in YEARS])
r, p = stats.pearsonr(mids, mean); out["class_midpoint_r"] = [float(r), float(p)]
out["class_means"] = {c: float(np.mean([m for y, m in zip(YEARS, mean) if CLASS[y] == c])) for c in MID}

# component sign check on the layers: annual mean of each raw layer vs annual mean CDI
comp = {}
for b, idx in [("MSI", 1), ("VHI", 2), ("WB", 3)]:   # band order of layers_YYYY.tif
    vals = []
    for y in YEARS:
        with rasterio.open(ROOT / "layers" / f"layers_{y}.tif") as d:
            x = d.read(idx).astype(float); x[x == d.nodata] = np.nan
            vals.append(np.nanmean(x))
    comp[b] = vals
    r, p = stats.pearsonr(vals, mean); out[f"cdi_vs_{b}_mean_r"] = [float(r), float(p)]
out["layer_annual_means"] = comp
r, p = stats.pearsonr(comp["VHI"], comp["MSI"]); out["VHI_MSI_annual_r"] = [float(r), float(p)]
r, p = stats.pearsonr(comp["VHI"], comp["WB"]); out["VHI_WB_annual_r"] = [float(r), float(p)]
r, p = stats.pearsonr(comp["MSI"], comp["WB"]); out["MSI_WB_annual_r"] = [float(r), float(p)]

# Landsat era
era = np.array([y <= 2012 for y in YEARS])
t = stats.ttest_ind(mean[era], mean[~era], equal_var=False)
out["era"] = dict(l57=float(mean[era].mean()), l8=float(mean[~era].mean()), t=float(t.statistic), p=float(t.pvalue))
dm = mean - np.where(era, mean[era].mean(), mean[~era].mean())
for ref, col in [("DSCI", "usdm_dsci_jja"), ("PDSI", "ncei_pdsi_jja")]:
    r, p = stats.pearsonr(dm, V[col]); out["era"][f"demeaned_{ref}"] = [float(r), float(p)]


# Williams' t for dependent correlations r(DSCI,CDI) vs r(DSCI,ref), n = 14 (Steiger 1980)
def williams(r12, r13, r23, n):
    detR = 1 - r12 ** 2 - r13 ** 2 - r23 ** 2 + 2 * r12 * r13 * r23
    rbar = (r12 + r13) / 2
    tt = (r12 - r13) * np.sqrt((n - 1) * (1 + r23) /
                               (2 * (n - 1) / (n - 3) * detR + rbar ** 2 * (1 - r23) ** 3))
    return float(tt), float(2 * stats.t.sf(abs(tt), n - 3))


ds = V.usdm_dsci_jja.values
for ref, col, sgn in [("PDSI", "ncei_pdsi_jja", -1), ("SPI3", "ncei_spi3_jja", -1)]:
    x = sgn * V[col].values  # orient so drought is positive, like the CDI
    out[f"williams_vs_{ref}"] = williams(np.corrcoef(ds, mean)[0, 1], np.corrcoef(ds, x)[0, 1],
                                         np.corrcoef(mean, x)[0, 1], len(YEARS))

# rankings
out["rank_cdi"] = [int(y) for y in np.array(YEARS)[np.argsort(-mean)]]
out["rank_dsci"] = [int(y) for y in np.array(YEARS)[np.argsort(-ds)]]
out["rank_pdsi"] = [int(y) for y in np.array(YEARS)[np.argsort(V.ncei_pdsi_jja.values)]]
r, p = stats.spearmanr(S.high_cdi_pixels_pct_global_q75, ds); out["q75share_spearman_dsci"] = [float(r), float(p)]
print(json.dumps(out, indent=1))
