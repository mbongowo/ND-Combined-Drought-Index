"""Figure 1: study area with physiographic regions.

Regions follow the physiographic divisions used in the text and are drawn from the
EPA Level III/IV ecoregions of North Dakota (Bryce et al., 1996; Omernik and Griffith, 2014):
  Red River Valley  = Level III 48 (Lake Agassiz Plain)
  Drift Prairie     = Level III 46 (Northern Glaciated Plains)
  Missouri Coteau   = Level IV 42a, 42b, 42d
  Missouri Plateau  = Level IV 42c, 42i and all of Level III 43
Both panels are drawn with equal axis scaling, so North Dakota keeps its true proportions.

usage: python fig1_study_area.py DATA_DIR OUT_PNG
"""
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from shapely.geometry import box

DATA = Path(sys.argv[1]); OUT = Path(sys.argv[2])
UTM = 26914; ALBERS = 5070

states = gpd.read_file(DATA / "states/cb_2023_us_state_20m.shp")
conus = states[~states.STUSPS.isin(["AK", "HI", "PR"])].to_crs(ALBERS)
nd = states[states.STUSPS == "ND"].to_crs(UTM)

eco = gpd.read_file(DATA / "nd_eco_l4/nd_eco_l4.shp").to_crs(UTM)
def region(r):
    if r.US_L3CODE == "48": return "Red River Valley"
    if r.US_L3CODE == "46": return "Drift Prairie"
    if r.US_L4CODE in ("42a", "42b", "42d"): return "Missouri Coteau"
    return "Missouri Plateau"
eco["region"] = eco.apply(region, axis=1)
reg = eco.dissolve("region").reset_index()
reg["geometry"] = reg.geometry.intersection(nd.geometry.iloc[0])
reg["km2"] = reg.geometry.area / 1e6
print(reg[["region", "km2"]])

COLORS = {"Missouri Plateau": "#e6d3a3", "Missouri Coteau": "#c9b27c",
          "Drift Prairie": "#b9d7a8", "Red River Valley": "#8fc1d4"}
LABELS = {"Missouri Plateau": (-102.9, 46.55), "Missouri Coteau": (-100.0, 46.9),
          "Drift Prairie": (-98.9, 47.75), "Red River Valley": (-97.25, 47.35)}

ndll = nd.to_crs(4326).geometry.iloc[0]
clipbox = gpd.GeoSeries([ndll.buffer(0.05)], crs=4326)
def clip(path, names, field="name"):
    g = gpd.read_file(path)
    g = g[g[field].isin(names)]
    return gpd.clip(g, clipbox.iloc[0]).to_crs(UTM)
rivers = clip(DATA / "rivers/ne_10m_rivers_lake_centerlines.shp", ["Missouri", "Red", "James"])
rivers2 = clip(DATA / "rivers_na/ne_10m_rivers_north_america.shp",
               ["Sheyenne", "Souris", "Little Missouri", "Heart", "Cannonball", "Pembina"])
lakes = clip(DATA / "lakes/ne_10m_lakes.shp", ["Lake Sakakawea", "Lake Oahe"])
lakes2 = clip(DATA / "lakes_na/ne_10m_lakes_north_america.shp", ["Devils Lake"])

cities = gpd.GeoDataFrame(
    {"name": ["Fargo", "Bismarck", "Grand Forks", "Minot", "Dickinson", "Williston"]},
    geometry=gpd.points_from_xy([-96.789, -100.784, -97.033, -101.296, -102.789, -103.618],
                                [46.877, 46.808, 47.925, 48.233, 46.879, 48.147]), crs=4326).to_crs(UTM)
OFFS = {"Fargo": (-8, -2, "right"), "Grand Forks": (-8, 4, "right"), "Bismarck": (6, -10, "left"),
        "Minot": (6, 4, "left"), "Dickinson": (6, 4, "left"), "Williston": (6, 4, "left")}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
fig = plt.figure(figsize=(12, 5.2), dpi=300)
axa = fig.add_axes([0.01, 0.10, 0.36, 0.78])
axb = fig.add_axes([0.39, 0.10, 0.60, 0.78])

# (a) locator
conus.plot(ax=axa, color="#f3efe4", edgecolor="#7a7a7a", linewidth=0.4)
conus[conus.STUSPS == "ND"].plot(ax=axa, color="#c0392b", edgecolor="black", linewidth=0.6)
axa.set_aspect("equal"); axa.set_axis_off()
axa.set_title("(a) North Dakota within the conterminous United States\n(Albers Equal Area, EPSG:5070)",
              fontsize=10, fontweight="bold")

# (b) study area
for name, g in reg.groupby("region"):
    g.plot(ax=axb, color=COLORS[name], edgecolor="none")
reg.boundary.plot(ax=axb, color="#5a4a2a", linewidth=0.7, linestyle=(0, (4, 2)))
lakes.plot(ax=axb, color="#6fa8dc", edgecolor="none"); lakes2.plot(ax=axb, color="#6fa8dc", edgecolor="none")
rivers.plot(ax=axb, color="#2e75b6", linewidth=1.1); rivers2.plot(ax=axb, color="#2e75b6", linewidth=0.6)
nd.boundary.plot(ax=axb, color="black", linewidth=1.0)
cities.plot(ax=axb, color="black", markersize=14, zorder=5)
for _, c in cities.iterrows():
    dx, dy, ha = OFFS[c["name"]]
    axb.annotate(c["name"], (c.geometry.x, c.geometry.y), xytext=(dx, dy), textcoords="offset points",
                 ha=ha, fontsize=8, zorder=6)
lab = gpd.GeoDataFrame({"name": list(LABELS)}, geometry=gpd.points_from_xy(
    [v[0] for v in LABELS.values()], [v[1] for v in LABELS.values()]), crs=4326).to_crs(UTM)
for _, r in lab.iterrows():
    axb.text(r.geometry.x, r.geometry.y, r["name"].replace(" ", "\n", 1) if r["name"] != "Drift Prairie" else r["name"],
             ha="center", va="center", fontsize=9, fontstyle="italic", fontweight="bold", color="#3b2f14", zorder=6)

x0, y0, x1, y1 = nd.total_bounds
pad = 15000
axb.set_xlim(x0 - pad, x1 + pad); axb.set_ylim(y0 - pad, y1 + pad); axb.set_aspect("equal")
# graticule labels
import numpy as np
from pyproj import Transformer
tr = Transformer.from_crs(4326, UTM, always_xy=True)
for lon in [-104, -102, -100, -98]:
    xs, ys = tr.transform([lon] * 50, np.linspace(45.8, 49.1, 50))
    axb.plot(xs, ys, color="#999999", linewidth=0.3, linestyle=":", zorder=1)
    xb, _ = tr.transform(lon, 45.935)
    axb.text(xb, y0 - pad - 4000, f"{abs(lon)}°W", ha="center", va="top", fontsize=7.5)
for lat in [46, 47, 48, 49]:
    xs, ys = tr.transform(np.linspace(-104.2, -96.4, 50), [lat] * 50)
    axb.plot(xs, ys, color="#999999", linewidth=0.3, linestyle=":", zorder=1)
    _, yb = tr.transform(-104.05, lat)
    axb.text(x0 - pad - 4000, yb, f"{lat}°N", ha="right", va="center", fontsize=7.5)
axb.set_xticks([]); axb.set_yticks([])
for s in axb.spines.values(): s.set_linewidth(0.6)
# scale bar and north arrow
sx, sy = x0 + 10000, y0 + 8000
for i in range(4):
    axb.add_patch(plt.Rectangle((sx + i * 25000, sy), 25000, 5000, facecolor="black" if i % 2 == 0 else "white",
                                edgecolor="black", linewidth=0.6, zorder=7))
for d in (0, 50, 100):
    axb.text(sx + d * 1000, sy + 8000, f"{d}", ha="center", va="bottom", fontsize=7.5, zorder=7)
axb.text(sx + 100000 + 6000, sy + 2500, "km", va="center", fontsize=7.5, zorder=7)
axb.annotate("N", xy=(x1 - 5000, y1 - 10000), xytext=(x1 - 5000, y1 - 55000), ha="center", fontsize=11,
             fontweight="bold", arrowprops=dict(arrowstyle="-|>", color="black", lw=1.5), zorder=7)
handles = [Patch(facecolor=COLORS[k], edgecolor="#5a4a2a", label=k) for k in
           ["Missouri Plateau", "Missouri Coteau", "Drift Prairie", "Red River Valley"]]
handles += [plt.Line2D([], [], color="#2e75b6", lw=1.1, label="Major rivers and reservoirs"),
            plt.Line2D([], [], marker="o", color="black", lw=0, markersize=4, label="Cities")]
axb.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=6, fontsize=8,
           frameon=False, handlelength=1.6, columnspacing=1.2)
axb.set_title("(b) North Dakota physiographic regions (NAD83 / UTM zone 14N, EPSG:26914)",
              fontsize=10, fontweight="bold")
fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
w, h = x1 - x0, y1 - y0
print(f"ND extent {w/1000:.1f} x {h/1000:.1f} km, ratio {w/h:.3f}")
