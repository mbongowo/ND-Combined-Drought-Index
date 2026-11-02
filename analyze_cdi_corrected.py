from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import rasterio


import sys
ROOT = Path(sys.argv[1])
OUT_DIR = ROOT
PREVIEW_DIR = ROOT / "CDI_final_package/preview_rasters"
ANALYSIS_DIR = OUT_DIR / "analysis"
GRAPH_DIR = ANALYSIS_DIR / "graphs"

YEARS = np.array(
    [2000, 2001, 2002, 2003, 2005, 2006, 2007, 2009, 2011, 2012, 2017, 2018, 2021, 2023],
    dtype=float,
)

CDI_CMAP = LinearSegmentedColormap.from_list(
    "cdi_drought",
    ["#1a9850", "#91cf60", "#fee08b", "#fdae61", "#d73027", "#7f0000"],
)
MISSING_COLOR = "#f4f1ea"


def ensure_dirs() -> None:
    ANALYSIS_DIR.mkdir(exist_ok=True)
    GRAPH_DIR.mkdir(exist_ok=True)


def read_preview(year: int) -> tuple[np.ndarray, tuple[float, ...]]:
    path = PREVIEW_DIR / f"CDI_{year}_utm14_preview.tif"
    with rasterio.open(path) as dataset:
        array = dataset.read(1).astype("float64")
        nodata = dataset.nodata
        t = dataset.transform
    if nodata is not None:
        array[array == nodata] = np.nan
    array[array <= -3.0e38] = np.nan
    return array, (t.c, t.a, t.b, t.f, t.d, t.e)


def valid_values(array: np.ndarray) -> np.ndarray:
    return array[np.isfinite(array)]


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def linear_regression(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    n = float(len(x))
    x_mean = float(np.mean(x))
    y_mean = float(np.mean(y))
    dx = x - x_mean
    dy = y - y_mean
    ssx = float(np.sum(dx * dx))
    slope = float(np.sum(dx * dy) / ssx)
    intercept = y_mean - slope * x_mean
    fitted = intercept + slope * x
    residual = y - fitted
    sse = float(np.sum(residual * residual))
    sst = float(np.sum((y - y_mean) ** 2))
    r2 = 1.0 - (sse / sst)
    rmse = math.sqrt(sse / n)
    se_slope = math.sqrt((sse / (n - 2.0)) / ssx)
    t_stat = slope / se_slope
    p_approx = math.erfc(abs(t_stat) / math.sqrt(2.0))
    return {
        "n": n,
        "slope_per_year": slope,
        "intercept": intercept,
        "r_squared": r2,
        "rmse": rmse,
        "std_error_slope": se_slope,
        "t_stat_approx": t_stat,
        "p_value_normal_approx": p_approx,
    }


def pixel_slopes(stack: np.ndarray) -> np.ndarray:
    x = YEARS - np.mean(YEARS)
    valid = np.isfinite(stack)
    counts = np.sum(valid, axis=0)
    sums = np.nansum(stack, axis=0)
    y_mean = np.full(stack.shape[1:], np.nan, dtype="float64")
    y_mean[counts > 0] = sums[counts > 0] / counts[counts > 0]
    numerator = np.nansum((stack - y_mean[None, :, :]) * x[:, None, None], axis=0)
    denominator = np.sum((x[:, None, None] ** 2) * valid, axis=0)
    slopes = np.full(stack.shape[1:], np.nan, dtype="float64")
    np.divide(numerator, denominator, out=slopes, where=denominator > 0)
    slopes[counts < 3] = np.nan
    return slopes


def build_analysis() -> None:
    ensure_dirs()
    mpl.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42})

    arrays: list[np.ndarray] = []
    geotransform: tuple[float, ...] | None = None
    for year in YEARS.astype(int):
        array, transform = read_preview(year)
        arrays.append(array)
        geotransform = geotransform or transform

    stack = np.stack(arrays)
    all_values = stack[np.isfinite(stack)]
    q25, q50, q75 = np.percentile(all_values, [25, 50, 75])
    pixel_area_km2 = abs(geotransform[1] * geotransform[5]) / 1_000_000.0 if geotransform else np.nan

    summary_rows: list[dict[str, object]] = []
    distributions: list[np.ndarray] = []
    for year, array in zip(YEARS.astype(int), arrays):
        values = valid_values(array)
        distributions.append(values)
        high_pct = float(np.mean(values >= q75) * 100.0)
        summary_rows.append(
            {
                "year": year,
                "valid_pixels": int(values.size),
                "approx_area_km2": round(float(values.size * pixel_area_km2), 3),
                "mean_cdi": round(float(np.mean(values)), 6),
                "median_cdi": round(float(np.median(values)), 6),
                "std_cdi": round(float(np.std(values, ddof=1)), 6),
                "min_cdi": round(float(np.min(values)), 6),
                "p10_cdi": round(float(np.percentile(values, 10)), 6),
                "p25_cdi": round(float(np.percentile(values, 25)), 6),
                "p75_cdi": round(float(np.percentile(values, 75)), 6),
                "p90_cdi": round(float(np.percentile(values, 90)), 6),
                "max_cdi": round(float(np.max(values)), 6),
                "high_cdi_pixels_pct_global_q75": round(high_pct, 3),
            }
        )

    summary_path = ANALYSIS_DIR / "CDI_summary_statistics.csv"
    write_csv(
        summary_path,
        list(summary_rows[0].keys()),
        summary_rows,
    )

    mean_values = np.array([float(row["mean_cdi"]) for row in summary_rows])
    median_values = np.array([float(row["median_cdi"]) for row in summary_rows])
    high_pct_values = np.array([float(row["high_cdi_pixels_pct_global_q75"]) for row in summary_rows])
    regression = linear_regression(YEARS, mean_values)

    change_rows: list[dict[str, object]] = []
    for idx in range(1, len(YEARS)):
        year = int(YEARS[idx])
        previous = int(YEARS[idx - 1])
        delta = mean_values[idx] - mean_values[idx - 1]
        change_rows.append(
            {
                "from_year": previous,
                "to_year": year,
                "year_gap": year - previous,
                "mean_cdi_change": round(float(delta), 6),
                "annualized_mean_cdi_change": round(float(delta / (year - previous)), 6),
            }
        )
    change_path = ANALYSIS_DIR / "CDI_year_to_year_changes.csv"
    write_csv(change_path, list(change_rows[0].keys()), change_rows)

    flattened = [array.reshape(-1) for array in arrays]
    corr = np.full((len(YEARS), len(YEARS)), np.nan)
    for i in range(len(YEARS)):
        for j in range(len(YEARS)):
            mask = np.isfinite(flattened[i]) & np.isfinite(flattened[j])
            corr[i, j] = np.corrcoef(flattened[i][mask], flattened[j][mask])[0, 1]

    corr_path = ANALYSIS_DIR / "CDI_correlation_matrix.csv"
    corr_rows = []
    for year, row in zip(YEARS.astype(int), corr):
        corr_rows.append({"year": year, **{str(int(y)): round(float(v), 6) for y, v in zip(YEARS, row)}})
    write_csv(corr_path, ["year", *[str(int(y)) for y in YEARS]], corr_rows)

    regression_path = ANALYSIS_DIR / "CDI_regression_summary.csv"
    write_csv(
        regression_path,
        ["metric", "value"],
        [{"metric": key, "value": round(value, 10)} for key, value in regression.items()],
    )

    fig, ax = plt.subplots(figsize=(9, 5.2), dpi=220)
    ax.plot(YEARS, mean_values, color="#7f0000", marker="o", linewidth=2.3, label="Mean CDI")
    ax.plot(YEARS, median_values, color="#2166ac", marker="s", linewidth=1.8, label="Median CDI")
    fitted = regression["intercept"] + regression["slope_per_year"] * YEARS
    ax.plot(YEARS, fitted, color="#333333", linestyle="--", linewidth=1.7, label="Linear trend")
    ax.set_title("CDI Annual Mean and Median Trend")
    ax.set_xlabel("Year")
    ax.set_ylabel("CDI")
    ax.grid(True, color="#dddddd", linewidth=0.7)
    ax.legend(frameon=False)
    fig.tight_layout()
    trend_graph = GRAPH_DIR / "CDI_mean_median_trend.png"
    fig.savefig(trend_graph, facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5.2), dpi=220)
    ax.bar(YEARS.astype(int), high_pct_values, color="#d73027", edgecolor="#7f0000", linewidth=0.5)
    ax.set_title("Share of Pixels in the Global Upper CDI Quartile")
    ax.set_xlabel("Year")
    ax.set_ylabel("Pixels at or above global 75th percentile (%)")
    ax.grid(True, axis="y", color="#dddddd", linewidth=0.7)
    fig.tight_layout()
    high_graph = GRAPH_DIR / "CDI_high_value_area_by_year.png"
    fig.savefig(high_graph, facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.5, 5.4), dpi=220)
    ax.boxplot(distributions, tick_labels=[str(int(y)) for y in YEARS], patch_artist=True, showfliers=False)
    for patch in ax.patches:
        patch.set_facecolor("#fee08b")
        patch.set_edgecolor("#a6611a")
    ax.set_title("CDI Distribution by Year")
    ax.set_xlabel("Year")
    ax.set_ylabel("CDI")
    ax.tick_params(axis="x", rotation=45)
    ax.grid(True, axis="y", color="#dddddd", linewidth=0.7)
    fig.tight_layout()
    box_graph = GRAPH_DIR / "CDI_distribution_boxplots.png"
    fig.savefig(box_graph, facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 7), dpi=220)
    image = ax.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_title("Spatial Correlation Between CDI Years")
    ax.set_xticks(np.arange(len(YEARS)), [str(int(y)) for y in YEARS], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(YEARS)), [str(int(y)) for y in YEARS])
    colorbar = fig.colorbar(image, ax=ax, shrink=0.82)
    colorbar.set_label("Pearson correlation")
    fig.tight_layout()
    corr_graph = GRAPH_DIR / "CDI_correlation_heatmap.png"
    fig.savefig(corr_graph, facecolor="white")
    plt.close(fig)

    slopes = pixel_slopes(stack)
    vmax = float(np.nanpercentile(np.abs(slopes), 98))
    fig, ax = plt.subplots(figsize=(8.5, 5.5), dpi=220)
    slope_cmap = mpl.colormaps["RdBu_r"].copy()
    slope_cmap.set_bad(MISSING_COLOR)
    im = ax.imshow(slopes, cmap=slope_cmap, vmin=-vmax, vmax=vmax, interpolation="nearest")
    ax.set_title("Pixel-wise CDI Linear Trend, 2000-2023")
    ax.set_xticks([])
    ax.set_yticks([])
    colorbar = fig.colorbar(im, ax=ax, shrink=0.82)
    colorbar.set_label("CDI change per year")
    fig.tight_layout()
    slope_graph = GRAPH_DIR / "CDI_pixelwise_linear_trend.png"
    fig.savefig(slope_graph, facecolor="white")
    plt.close(fig)

    top_mean = summary_rows[int(np.argmax(mean_values))]
    low_mean = summary_rows[int(np.argmin(mean_values))]
    top_high = summary_rows[int(np.argmax(high_pct_values))]
    avg_neighbor_corr = float(np.mean([corr[i, i + 1] for i in range(len(YEARS) - 1)]))

    report_path = ANALYSIS_DIR / "CDI_interpretation_report.md"
    report_path.write_text(
        "\n".join(
            [
                "# Combined Drought Index (CDI) Statistical Interpretation",
                "",
                "## Data and Method",
                f"The analysis uses the 14 NAD83 / UTM zone 14N preview rasters in `{PREVIEW_DIR}`. "
                "All years are evaluated on the same 700 by 426 grid, using valid non-nodata pixels only.",
                "",
                "## Map Interpretation",
                "The updated CDI map uses a drought-oriented green-to-yellow-to-orange-to-red palette. "
                "Lower CDI values appear green, intermediate values appear yellow/orange, and higher CDI values appear red to dark maroon. "
                "Because a single shared color scale is used for every panel, color differences can be compared directly across years.",
                "",
                "## Main Statistical Findings",
                f"The overall CDI range across the preview rasters is {np.nanmin(stack):.4f} to {np.nanmax(stack):.4f}. "
                f"The global quartile thresholds are Q25={q25:.4f}, median={q50:.4f}, and Q75={q75:.4f}.",
                f"The highest mean CDI year is {int(top_mean['year'])} with mean CDI={float(top_mean['mean_cdi']):.4f}.",
                f"The lowest mean CDI year is {int(low_mean['year'])} with mean CDI={float(low_mean['mean_cdi']):.4f}.",
                f"The largest share of high-CDI pixels occurs in {int(top_high['year'])}, where {float(top_high['high_cdi_pixels_pct_global_q75']):.2f}% of valid pixels are at or above the global 75th percentile.",
                "",
                "## Regression Analysis",
                f"A simple linear regression of annual mean CDI against year gives a slope of {regression['slope_per_year']:.6f} CDI units per year "
                f"(about {regression['slope_per_year'] * 10.0:.5f} per decade), with R2={regression['r_squared']:.3f}.",
                f"The approximate t statistic for the slope is {regression['t_stat_approx']:.3f} "
                f"(normal-approximation p={regression['p_value_normal_approx']:.4f}).",
                "Interpretation: a positive slope indicates increasing mean CDI through time; a negative slope indicates decreasing mean CDI. "
                "R2 describes how much of the year-to-year mean CDI variation is explained by a straight-line trend.",
                "",
                "## Correlation Analysis",
                f"The mean correlation between adjacent observed years is {avg_neighbor_corr:.3f}. "
                "High positive correlations indicate that the spatial pattern of CDI is broadly persistent between years, while lower correlations indicate stronger spatial reorganization.",
                "",
                "## Output Files",
                f"- Summary statistics table: `{summary_path}`",
                f"- Year-to-year change table: `{change_path}`",
                f"- Correlation matrix: `{corr_path}`",
                f"- Regression summary: `{regression_path}`",
                f"- Mean/median trend graph: `{trend_graph}`",
                f"- High-CDI area graph: `{high_graph}`",
                f"- Distribution boxplots: `{box_graph}`",
                f"- Correlation heatmap: `{corr_graph}`",
                f"- Pixel-wise trend map: `{slope_graph}`",
            ]
        ),
        encoding="utf-8",
    )

    print(f"Wrote {summary_path}")
    print(f"Wrote {change_path}")
    print(f"Wrote {corr_path}")
    print(f"Wrote {regression_path}")
    print(f"Wrote {report_path}")
    print(f"Wrote graphs to {GRAPH_DIR}")


if __name__ == "__main__":
    build_analysis()
