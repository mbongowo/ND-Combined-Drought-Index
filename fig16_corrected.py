"""Figure 16: VHI against MSI and against WB, fourteen annual statewide means.

Uses the June-August layers of the corrected CDI (layers_YYYY.tif, band order
MSI, VHI, WB) averaged over the valid cells of each analyzed year.

usage: python fig16_corrected.py ROOT OUT_PNG
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from scipy import stats

ROOT = Path(sys.argv[1]); OUT = sys.argv[2]
YEARS = [2000, 2001, 2002, 2003, 2005, 2006, 2007, 2009, 2011, 2012, 2017, 2018, 2021, 2023]
CLASS = {2002: "Extreme", 2006: "Extreme", 2012: "Extreme", 2021: "Extreme",
         2000: "Severe", 2003: "Severe", 2007: "Severe", 2017: "Severe", 2023: "Severe",
         2001: "Moderate", 2005: "Moderate", 2009: "Moderate", 2011: "Moderate", 2018: "Moderate"}
COL = {"Extreme": "#b2182b", "Severe": "#ef8a62", "Moderate": "#4393c3"}

msi, vhi, wb = [], [], []
for y in YEARS:
    with rasterio.open(ROOT / "layers" / f"layers_{y}.tif") as d:
        a = d.read([1, 2, 3]).astype(float); a[a == d.nodata] = np.nan
    msi.append(np.nanmean(a[0])); vhi.append(np.nanmean(a[1])); wb.append(np.nanmean(a[2]))
msi, vhi, wb = map(np.array, (msi, vhi, wb))

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.4))
for ax, x, lab, tag in [(axes[0], msi, "Moisture Stress Index (MSI, June–August mean)", "(a) VHI versus MSI"),
                        (axes[1], wb, "Water balance (WB = P − AET, mm month⁻¹, June–August mean)", "(b) VHI versus WB")]:
    fit = stats.linregress(x, vhi)
    xx = np.linspace(x.min(), x.max(), 50)
    ax.plot(xx, fit.intercept + fit.slope * xx, "--", color="#555555", lw=1.3)
    for c in COL:
        m = np.array([CLASS[y] == c for y in YEARS])
        ax.scatter(x[m], vhi[m], s=42, color=COL[c], edgecolor="k", lw=0.5, label=f"{c} drought (Section 3.1)", zorder=3)
    for xi, yi, y in zip(x, vhi, YEARS):
        ax.annotate(str(y), (xi, yi), xytext=(4, 3), textcoords="offset points", fontsize=7.5)
    ax.set_title(f"{tag}:  slope = {fit.slope:.3g}, R² = {fit.rvalue**2:.2f}, p = {fit.pvalue:.2g}, n = 14".replace("-", "−"), fontsize=10, loc="left")
    ax.set_xlabel(lab); ax.set_ylabel("VHI (0–1, June–August mean)")
    ax.grid(True, color="#dddddd", lw=0.6)
    print(tag, "slope", fit.slope, "intercept", fit.intercept, "r", fit.rvalue, "R2", fit.rvalue**2, "p", fit.pvalue)
axes[0].legend(frameon=False, fontsize=8, loc="upper right")
fig.tight_layout()
fig.savefig(OUT, dpi=250, facecolor="white")
