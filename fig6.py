import sys, numpy as np, rasterio, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

PROJ = Path(sys.argv[1]); OUT = Path(sys.argv[2])
YEARS = [2000,2001,2002,2003,2005,2006,2007,2009,2011,2012,2017,2018,2021,2023]
arrs = []
for y in YEARS:
    with rasterio.open(PROJ/"CDI_final_package/preview_rasters"/f"CDI_{y}_utm14_preview.tif") as d:
        a = d.read(1).astype(float); a[(a == d.nodata) | ~np.isfinite(a)] = np.nan
        arrs.append(a); b = d.bounds
stack = np.stack(arrs)
shared = np.all(np.isfinite(stack), axis=0)
vals = stack[:, shared]
q75 = np.quantile(vals, 0.75)
vmin, vmax = np.floor(np.nanmin(vals)*20)/20, np.ceil(np.nanmax(vals)*20)/20
ext = [b.left/1000, b.right/1000, b.bottom/1000, b.top/1000]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
fig, axes = plt.subplots(4, 4, figsize=(11, 7.6))
cmap = plt.get_cmap("RdYlGn_r").copy(); cmap.set_bad("white")
for ax in axes.flat: ax.axis("off")
for i, y in enumerate(YEARS):
    ax = axes.flat[i]
    a = np.where(shared, stack[i], np.nan)
    im = ax.imshow(a, cmap=cmap, vmin=vmin, vmax=vmax, extent=ext, interpolation="nearest")
    ax.set_title(f"{y}   mean {np.nanmean(a):.3f} | {100*np.nanmean(a[shared] >= q75):.1f}% ≥ Q75", fontsize=8.5)
# scale bar + north arrow in spare panels
ax = axes.flat[14]  # same extent and aspect as the map panels, so the bar is true to scale
x0 = ext[0] + 20; y0 = ext[2] + 0.55 * (ext[3] - ext[2])
ax.plot([x0, x0 + 300], [y0, y0], color="k", lw=4); ax.text(x0 + 150, y0 + 25, "300 km", ha="center", fontsize=9)
ax.text(x0, ext[2] + 40, "NAD83 / UTM zone 14N (EPSG:26914)\nCell size 837.9 m; shared valid footprint", fontsize=7.5)
ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal")
ax = axes.flat[15]; ax.annotate("N", xy=(0.5, 0.9), xytext=(0.5, 0.35), ha="center", fontsize=12,
                                arrowprops=dict(arrowstyle="-|>", lw=2, color="k"), xycoords="axes fraction")
cax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
cb = fig.colorbar(im, cax=cax); cb.set_label("CDI value (low → high drought-positive signal)")
cb.ax.axhline(q75, color="k", lw=1.2); cb.ax.text(1.6, q75, "Q75", va="center", fontsize=8, transform=cb.ax.get_yaxis_transform())
fig.suptitle("Annual Combined Drought Index, North Dakota (fourteen analyzed years, 2000–2023)", fontsize=12, fontweight="bold")
fig.subplots_adjust(left=0.02, right=0.91, top=0.92, bottom=0.03, wspace=0.05, hspace=0.25)
fig.savefig(OUT, dpi=250)
m = axes.flat[13]; t = axes.flat[14]
print("map px/km", m.transData.transform((ext[0]+300,0))[0]-m.transData.transform((ext[0],0))[0], "bar", t.transData.transform((ext[0]+300,0))[0]-t.transData.transform((ext[0],0))[0])
print(vmin, vmax, q75)
