#!/usr/bin/env python3
"""Fig 7 per-cancer spatial-stats panels (C1 forest + C2 stats).
Source: metabolism_target/spatial_deep_full_cohort.csv (132 samples).
Outputs: fig7_C1_forest.{svg,png,pdf}, fig7_C2_stats.{svg,png,pdf}
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "Liberation Sans", "svg.fonttype": "none", "pdf.fonttype": 42})

CSV = "data/hub_tls_paper_outline/metabolism_target/spatial_deep_full_cohort.csv"
OUTDIR = "data/hub_tls_paper_outline/figure7_assembled/panels"

FIX = {
    "GSE175540_Lung_Cancer": "Renal",
    "GSE189843_Breast_Cancer": "Ovarian",
    "GSE208253_Ovarian_Cancer": "Oral",
    "GSE211956_PDAC": "Ovarian",
    "BreastCancer_Spatial_CID": "Breast",
    "GSE213688_Breast": "Breast",
    "GSE237183_GBM": "GBM",
    "GSE238264_HCC": "HCC",
    "GSE246011_Bladder_Cancer": "Bladder",
    "GSE246011_Gastric_Cancer": "Gastric",
    "GSE251950_Gastric_Cancer": "Gastric",
    "GSE267401_Colorectal_Cancer": "Colorectal",
    "GSE278687_PDAC": "PDAC",
    "GSE308486_Prostate_Cancer": "Prostate",
}

df = pd.read_csv(CSV)
df["cancer"] = df["dataset"].map(FIX)
assert df["cancer"].notna().all()

C_COLORS = {
    "Renal": "#4E79A7", "Breast": "#F28E2B", "Ovarian": "#76B7B2",
    "Oral": "#E15759", "PDAC": "#59A14F", "GBM": "#EDC948",
    "HCC": "#B07AA1", "Bladder": "#FF9DA7", "Gastric": "#9C755F",
    "Colorectal": "#BAB0AC", "Prostate": "#86BCB6",
}
SIG_NEG, SIG_POS, NS = "#C41E3A", "#6C8EBF", "#C4C4C4"

# ---------- C1: forest (two sub-panels) ----------
rho = df.dropna(subset=["rho_mat", "p_mat"])
n_sig_neg = int(((rho["p_mat"] < 0.05) & (rho["rho_mat"] < 0)).sum())
n_sig_pos = int(((rho["p_mat"] < 0.05) & (rho["rho_mat"] > 0)).sum())
n_ns = int((rho["p_mat"] >= 0.05).sum())
n_rho = len(rho)

lor = df.dropna(subset=["lor_e_vs_sfl", "se_lor"])
# REML (DerSimonian-Laird)
y, v = lor["lor_e_vs_sfl"].values, lor["se_lor"].values ** 2
w = 1 / v
mu_f = (y * w).sum() / w.sum()
Q = float((w * (y - mu_f) ** 2).sum())
C = w.sum() - (w ** 2).sum() / w.sum()
tau2 = max(0.0, (Q - (len(y) - 1)) / C)
wr = 1 / (v + tau2)
mu = float((y * wr).sum() / wr.sum())
from scipy import stats as st
p_reml = 2 * st.norm.sf(abs(mu / np.sqrt(1 / wr.sum())))

fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.4, 4.4), dpi=300)

order_L = rho.groupby("cancer")["rho_mat"].median().sort_values(ascending=False).index.tolist()
for i, c in enumerate(order_L):
    sub = rho[rho["cancer"] == c]
    yy = np.full(len(sub), i) + np.random.default_rng(0).normal(0, 0.12, len(sub))
    cols = np.where(sub["p_mat"] >= 0.05, NS, np.where(sub["rho_mat"] < 0, SIG_NEG, SIG_POS))
    axL.scatter(sub["rho_mat"], yy, s=7, c=cols, alpha=0.85, linewidths=0)
    med = sub["rho_mat"].median()
    axL.plot([med, med], [i - 0.35, i + 0.35], color=C_COLORS[c], lw=3, solid_capstyle="round")
axL.axvline(0, ls="--", lw=0.8, color="black")
axL.set_yticks(range(len(order_L)))
axL.set_yticklabels([f"{c} (n={int((rho['cancer'] == c).sum())})" for c in order_L], fontsize=7)
axL.set_xlabel("Within-sample Spearman ρ (IDO1 vs TLS maturity)", fontsize=8)
axL.set_title(f"Within-sample IDO1–maturity correlation\n(sig-negative {n_sig_neg} / sig-positive {n_sig_pos} / n.s. {n_ns} of {n_rho} samples)", fontsize=8)
axL.tick_params(axis="x", labelsize=7)

order_R = lor.groupby("cancer")["lor_e_vs_sfl"].median().sort_values(ascending=False).index.tolist()
for i, c in enumerate(order_R):
    sub = lor[lor["cancer"] == c].sort_values("lor_e_vs_sfl")
    yy = np.arange(len(sub)) * 0.0 + i
    yy = np.linspace(i - 0.32, i + 0.32, len(sub))
    axR.errorbar(sub["lor_e_vs_sfl"], yy, xerr=1.96 * sub["se_lor"],
                 fmt="o", ms=3, color=C_COLORS[c], ecolor=C_COLORS[c], elinewidth=0.5,
                 alpha=0.8)
axR.axvline(0, ls="-", lw=0.8, color="black")
axR.axvline(mu, ls=":", lw=1.2, color="black")
axR.text(mu, len(order_R) - 0.2, f" REML pooled = {mu:.2f} (p = {p_reml:.1e})", fontsize=7, va="top")
axR.set_yticks(range(len(order_R)))
axR.set_yticklabels([f"{c} (n={int((lor['cancer'] == c).sum())})" for c in order_R], fontsize=7)
axR.set_xlabel("log OR (E-TLS vs SFL-TLS membership of IDO1+ spots)", fontsize=8)
axR.set_title(f"IDO1+ enrichment: E-TLS vs SFL-TLS (n = {len(lor)} samples)", fontsize=8)
axR.tick_params(axis="x", labelsize=7)

for ax in (axL, axR):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("svg", "png", "pdf"):
    fig.savefig(f"{OUTDIR}/fig7_C1_forest.{ext}", dpi=300)

# ---------- C2: spatial stats (three sub-panels) ----------
fig2, (a1, a2, a3) = plt.subplots(1, 3, figsize=(10.5, 3.4), dpi=300)
cancers = sorted(df["cancer"].unique())
for c in cancers:
    sub = df[df["cancer"] == c]
    a1.scatter(sub["I_mat"], sub["I_ido"], s=7, color=C_COLORS[c], alpha=0.75, label=c, linewidths=0)
lim = [min(df["I_mat"].min(), df["I_ido"].min()) - 0.05, max(df["I_mat"].max(), df["I_ido"].max()) + 0.05]
a1.plot(lim, lim, ls=":", lw=0.8, color="grey")
a1.set_xlim(lim); a1.set_ylim(lim)
a1.set_xlabel("Moran's I of TLS maturity (control)", fontsize=8)
a1.set_ylabel("Moran's I of IDO1", fontsize=8)
a1.set_title("IDO1 clusters at similar scale\nas TLS maturity itself", fontsize=8)
a1.legend(fontsize=5.5, frameon=False, ncol=2, markerscale=0.9, handletextpad=0.15, borderaxespad=0.1)

d2 = df.dropna(subset=["rho_mat", "rho_d"])
for c in cancers:
    sub = d2[d2["cancer"] == c]
    a2.scatter(sub["rho_mat"], sub["rho_d"], s=7, color=C_COLORS[c], alpha=0.75, linewidths=0)
a2.axhline(0, ls=":", lw=0.8, color="grey"); a2.axvline(0, ls=":", lw=0.8, color="grey")
a2.set_xlabel("ρ (IDO1 vs maturity), within sample", fontsize=8)
a2.set_ylabel("ρ (IDO1 vs distance to SFL-TLS)", fontsize=8)
a2.set_title("Two spatial axes per sample\n(brake quadrant: left & up)", fontsize=8)

d3 = df.dropna(subset=["rho_d"])
order3 = d3.groupby("cancer")["rho_d"].median().sort_values().index.tolist()
data3 = [d3[d3["cancer"] == c]["rho_d"].values for c in order3]
bp = a3.boxplot(data3, vert=False, patch_artist=True, showfliers=True,
                flierprops=dict(ms=2, marker="o", mfc="grey", mec="grey", alpha=0.5),
                medianprops=dict(color="black", lw=0.8), widths=0.6)
for patch, c in zip(bp["boxes"], order3):
    patch.set_facecolor(C_COLORS[c]); patch.set_alpha(0.75); patch.set_edgecolor("grey")
a3.axvline(0, ls="--", lw=0.8, color="black")
a3.set_yticklabels(order3, fontsize=7)
a3.set_xlabel("ρ (IDO1 vs distance to nearest SFL-TLS), per sample", fontsize=8)
a3.set_title("Distance-to-SFL gradient by cancer type\n(positive = IDO1 avoids mature TLS)", fontsize=8)
a3.tick_params(axis="x", labelsize=7)

for ax in (a1, a2, a3):
    ax.tick_params(axis="x", labelsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig2.tight_layout()
for ext in ("svg", "png", "pdf"):
    fig2.savefig(f"{OUTDIR}/fig7_C2_stats.{ext}", dpi=300)

print(f"C1: {n_rho} samples with rho_mat | sig-neg {n_sig_neg}, sig-pos {n_sig_pos}, ns {n_ns}")
print(f"    lor arm n={len(lor)} | REML logOR={mu:.3f}, p={p_reml:.1e}, tau2={tau2:.3f}")
print(f"C2 done. Cancers: {len(cancers)}")
