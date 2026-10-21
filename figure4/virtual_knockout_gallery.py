#!/usr/bin/env python3
# virtual_knockout_gallery.py
# Render Figure 4 panels g (PCA displacement field) and h (spatial perturbation
# maps) for the Hub-TLS manuscript.
#
# Panel h uses the TLS-module readout (consistent with panels e/f/i), with
# leave-one-out de-circularized deltas and a per-gene LOO reference recomputed
# from the underlying matrix. Original TLS score is shown in paper green;
# collapse delta in paper blue; a black contour anchors the hub region on all
# four spatial maps.
#
# Design:
#   1. Hub-high contour on all four spatial maps; TLS-like basin contour on
#      the PCA field.
#   2. Paper palette on panel h; panel g keeps magma/coolwarm with the basin
#      contour added.
#   3. Liberation Sans (Arial-compatible) throughout.
#
# Outputs: ko_spatial_tls_delta_enhanced.*, ko_pca_field_enhanced.*

from pathlib import Path
import importlib.util

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA

SCRIPTS = Path("data/spatial/scripts_pseudotime_v2")
SPEC = importlib.util.spec_from_file_location("flashy47", SCRIPTS / "47_flashy_knockout_figures.py")
F47 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(F47)

FLASHY = Path("data/spatial/hub_precursor_analysis/virtual_knockout/figures/flashy")
CSV = FLASHY / "per_spot_ko_deltas.csv"
GENES = ["CCL19", "CCL21", "CXCL13"]
MODULES = ["Hub", "TLS", "Bcell", "Tcell", "Chemo_core", "HEV"]

plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.titleweight": "bold",
    "axes.labelsize": 8,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

CMAP_HUB = LinearSegmentedColormap.from_list("hubred", ["#FFFFFF", "#F2C2CC", "#C41E3A", "#7A0F24"])
CMAP_COLLAPSE = LinearSegmentedColormap.from_list("collapse", ["#2F4C7A", "#6C8EBF", "#C9D6E8", "#FFFFFF"])
BASIN_COL = "#2F4C7A"


def load():
    df = pd.read_csv(CSV)
    print(f"loaded {CSV.name}: {df.shape[0]} spots")
    return df


def hub_contour_levels(df, q=0.75):
    return np.quantile(df["orig_Hub"], q)


def gridify(df, value="orig_Hub"):
    gx = np.sort(df["array_x"].unique())
    gy = np.sort(df["array_y"].unique())
    grid = np.full((len(gy), len(gx)), np.nan)
    xi = np.searchsorted(gx, df["array_x"].values)
    yi = np.searchsorted(gy, df["array_y"].values)
    grid[yi, xi] = df[value].values
    smooth = gaussian_filter(np.nan_to_num(grid, nan=0.0), sigma=1.5)
    return gx, gy, smooth


COMBO_KEYS = ["CCL19_CCL21", "CCL19_CCL21_CXCL13"]
COMBO_LABELS = {"CCL19_CCL21": "CCL19+CCL21 KO", "CCL19_CCL21_CXCL13": "Triple KO"}
REF_KEYS = GENES + COMBO_KEYS

STATE_COLORS = {"TLS-like": "#5DA05D", "Hub-like": "#6C8EBF", "Other": "#E3E3E3"}


def compute_loo_reference(df):
    """Per-gene (and per-combo) LOO reference module scores, all six modules.
    The saved per_spot_ko_deltas_loo.csv keeps only the last gene's LOO module
    set in its orig_ columns (loop overwrite in 47), biasing CCL19/CCL21
    reference frames; recompute each key's own LOO reference from the
    unperturbed matrix (no graph needed). Cached to loo_reference_all.csv so
    style-only iterations skip the matrix."""
    cache = FLASHY / "loo_reference_all.csv"
    if cache.exists():
        print(f"loaded LOO reference from cache ({cache.name})")
        return pd.read_csv(cache)
    adata = F47.load_inputs()[0]
    assert adata.n_obs == len(df), "spot count mismatch"
    full_tls = F47.score_modules(adata)["TLS"].values
    corr = np.corrcoef(full_tls, df["orig_TLS"])[0, 1]
    print(f"alignment check: corr(rescored full TLS, csv orig_TLS) = {corr:.4f}")
    ref = {}
    for key in REF_KEYS:
        gs = key.split("_")
        scores = F47.score_modules(adata, modules=F47.make_loo_modules(gs))
        for m in scores.columns:
            ref[f"{key}__{m}"] = scores[m].values
    ref_df = pd.DataFrame(ref)
    ref_df.to_csv(cache, index=False)
    print(f"cached LOO reference -> {cache.name}")
    return ref_df


def plot_spatial(df, df_loo, ref_df):
    fig, axes = plt.subplots(1, 4, figsize=(12.8, 3.4))
    gx, gy, grid = gridify(df)
    level = hub_contour_levels(df)

    def draw_contour(ax):
        ax.contour(gx, gy, grid, levels=[level], colors="black", linewidths=1.1)

    ax = axes[0]
    scat = ax.scatter(df["array_x"], df["array_y"], c=df["orig_TLS"], s=2,
                      cmap=CMAP_HUB, rasterized=True)
    draw_contour(ax)
    ax.set_title("Original TLS score")
    plt.colorbar(scat, ax=ax, fraction=0.046, pad=0.04, label="TLS score")

    for ax, gene in zip(axes[1:], GENES):
        delta = df_loo[f"{gene}_TLS"].values - ref_df[f"{gene}__TLS"].values
        vmin = np.quantile(delta, 0.01)
        scat = ax.scatter(df["array_x"], df["array_y"], c=delta, s=2,
                          cmap=CMAP_COLLAPSE, norm=Normalize(vmin=vmin, vmax=0),
                          rasterized=True)
        draw_contour(ax)
        ax.set_title(f"Δ TLS after {gene} KO (LOO)")
        plt.colorbar(scat, ax=ax, fraction=0.046, pad=0.04, label="Δ TLS score")

    for ax in axes:
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.axis("off")
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FLASHY / f"ko_spatial_tls_delta_enhanced.{ext}",
                    dpi=600 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print("saved ko_spatial_tls_delta_enhanced.*")


def classify_states(tls, hub, tls_t, hub_t):
    s = np.full(len(tls), "Other", dtype=object)
    s[tls > tls_t] = "TLS-like"
    s[(tls <= tls_t) & (hub > hub_t)] = "Hub-like"
    return s


def plot_state_dissolution(df, df_loo, df_combo, ref_df):
    """Six-map series: state map before KO, after each single KO, and after
    double/triple KO. States use the same median-threshold rule as 47, with
    each panel's thresholds taken from its own LOO reference (median of the
    reference = 50% TLS-like by construction)."""
    gx, gy, grid = gridify(df)
    level = hub_contour_levels(df)

    panels = [("Before KO", df["orig_TLS"].values, df["orig_Hub"].values,
               np.median(df["orig_TLS"]), np.median(df["orig_Hub"]))]
    for g in GENES:
        panels.append((f"{g} KO", df_loo[f"{g}_TLS"].values, df_loo[f"{g}_Hub"].values,
                       np.median(ref_df[f"{g}__TLS"]), np.median(ref_df[f"{g}__Hub"])))
    for key in COMBO_KEYS:
        panels.append((COMBO_LABELS[key], df_combo[f"{key}_TLS"].values,
                       df_combo[f"{key}_Hub"].values,
                       np.median(ref_df[f"{key}__TLS"]), np.median(ref_df[f"{key}__Hub"])))

    n0 = None
    fig, axes = plt.subplots(1, len(panels), figsize=(3.1 * len(panels), 3.6))
    for ax, (title, tls, hub, tls_t, hub_t) in zip(axes, panels):
        state = classify_states(tls, hub, tls_t, hub_t)
        n_tls = int((state == "TLS-like").sum())
        if n0 is None:
            n0 = n_tls
        retention = 100.0 * n_tls / n0
        for st, size in [("Other", 1), ("Hub-like", 2.2), ("TLS-like", 2.6)]:
            m = state == st
            ax.scatter(df["array_x"].values[m], df["array_y"].values[m],
                       c=STATE_COLORS[st], s=size, rasterized=True, linewidths=0)
        ax.contour(gx, gy, grid, levels=[level], colors="black", linewidths=1.0)
        ax.set_title(title, fontsize=9)
        ax.text(0.5, -0.045, f"TLS-like retained: {retention:.0f}%",
                transform=ax.transAxes, ha="center", fontsize=7.5, color="#404040")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.axis("off")
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=c,
                      markeredgecolor="none", markersize=6, label=s)
               for s, c in STATE_COLORS.items()]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FLASHY / f"ko_state_dissolution_enhanced.{ext}",
                    dpi=600 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print("saved ko_state_dissolution_enhanced.*")


def plot_pca(df):
    orig = df[[f"orig_{m}" for m in MODULES]].values
    pca = PCA(n_components=2, random_state=1)
    orig_pca = pca.fit_transform(orig)

    tls_like = df["orig_TLS"] > df["orig_TLS"].median()
    basin_xy = orig_pca[tls_like.values]
    if len(basin_xy) > 5000:
        rng = np.random.default_rng(1)
        basin_xy = basin_xy[rng.choice(len(basin_xy), 5000, replace=False)]
    kde = gaussian_kde(basin_xy.T, bw_method=0.25)
    gx = np.linspace(orig_pca[:, 0].min(), orig_pca[:, 0].max(), 120)
    gy = np.linspace(orig_pca[:, 1].min(), orig_pca[:, 1].max(), 120)
    GX, GY = np.meshgrid(gx, gy)
    dens = kde(np.vstack([GX.ravel(), GY.ravel()])).reshape(GX.shape)
    dens /= dens.max()

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0))
    for ax, gene in zip(axes, GENES):
        ko = df[[f"{gene}_{m}" for m in MODULES]].values
        ko_pca = pca.transform(ko)
        delta = ko_pca - orig_pca
        step = max(1, len(orig_pca) // 300)
        idx = np.arange(0, len(orig_pca), step)

        ax.scatter(orig_pca[:, 0], orig_pca[:, 1], c=df["orig_TLS"], cmap="magma",
                   s=1, alpha=0.3, rasterized=True)
        mag = np.sqrt(delta[idx, 0] ** 2 + delta[idx, 1] ** 2)
        q = ax.quiver(orig_pca[idx, 0], orig_pca[idx, 1], delta[idx, 0], delta[idx, 1],
                      mag, cmap="coolwarm", alpha=0.85,
                      scale=8, width=0.004, headwidth=3, headlength=4)
        # white halo under the basin contour keeps it readable on the magma background
        ax.contour(GX, GY, dens, levels=[0.35, 0.7], colors=["white"],
                   linewidths=[2.4, 2.8], linestyles=["dashed", "solid"])
        ax.contour(GX, GY, dens, levels=[0.35, 0.7], colors=[BASIN_COL],
                   linewidths=[0.9, 1.4], linestyles=["dashed", "solid"])
        ax.set_title(f"{gene} KO displacement")
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2" if ax is axes[0] else None)
        if ax is axes[-1]:
            plt.colorbar(q, ax=ax, fraction=0.046, pad=0.04,
                         label="Displacement magnitude")
            ax.legend(handles=[Line2D([0], [0], color=BASIN_COL, lw=1.2,
                                      label="TLS-like basin")],
                      loc="upper right", fontsize=7, frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FLASHY / f"ko_pca_field_enhanced.{ext}",
                    dpi=600 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print("saved ko_pca_field_enhanced.*")


def main():
    df = load()
    df_loo = pd.read_csv(FLASHY / "per_spot_ko_deltas_loo.csv")
    df_combo = pd.read_csv(FLASHY / "per_spot_combo_ko_deltas.csv")
    ref_df = compute_loo_reference(df)
    plot_spatial(df, df_loo, ref_df)
    plot_pca(df)
    plot_state_dissolution(df, df_loo, df_combo, ref_df)


if __name__ == "__main__":
    main()
