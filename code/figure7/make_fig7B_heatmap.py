#!/usr/bin/env python3
"""Generate the fig7 per-cancer IDO1 axis heatmap.

Canonical cancer labels are assigned per dataset via the DATASET_TO_CANCER
mapping below. Source: metabolism_target/hub_spot_candidate_genes_with_pt.csv
(37,230 hub spots).
Output: fig7_B_heatmap.{svg,png,pdf} + fig7_B_heatmap.csv
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "Liberation Sans",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

CSV = "data/hub_tls_paper_outline/metabolism_target/hub_spot_candidate_genes_with_pt.csv"
OUT = "data/hub_tls_paper_outline/figure7_assembled/panels/fig7_B_heatmap"

FIX = {
    "GSE175540_Lung_Cancer": "Renal",
    "GSE189843_Breast_Cancer": "Ovarian",
    "GSE208253_Ovarian_Cancer": "Oral",
    "GSE211956_PDAC": "Ovarian",
    "GSE246011_Lung_Cancer": "Lung",
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

STAGES = ["E-TLS", "PFL-TLS", "SFL-TLS"]

df = pd.read_csv(CSV)
df["cancer"] = df["dataset"].map(FIX)
assert df["cancer"].notna().all(), "unmapped dataset"
df = df[df["tls_stage"].isin(STAGES)]
df["ido1_pos"] = df["IDO1"] > 0

agg = df.groupby(["cancer", "tls_stage"]).agg(
    mean_ido1=("IDO1", "mean"), frac_pos=("ido1_pos", "mean"),
).reset_index()
mean_piv = agg.pivot(index="cancer", columns="tls_stage", values="mean_ido1")[STAGES]
frac_piv = agg.pivot(index="cancer", columns="tls_stage", values="frac_pos")[STAGES]
mean_piv = mean_piv.sort_values("E-TLS", ascending=False)
frac_piv = frac_piv.loc[mean_piv.index]

n_decline = int((mean_piv["E-TLS"] > mean_piv["SFL-TLS"]).sum())
print(f"cancer types: {len(mean_piv)}, E>SFL decline: {n_decline}")
print(mean_piv.round(3))

agg_out = agg.copy()
agg_out["tls_stage"] = pd.Categorical(agg_out["tls_stage"], STAGES, ordered=True)
agg_out.sort_values(["cancer", "tls_stage"]).to_csv(OUT + ".csv", index=False)

fig, ax = plt.subplots(figsize=(4.6, 5.6), dpi=300)
vals = mean_piv.values
im = ax.imshow(vals, cmap="magma", aspect="auto", vmin=0)

for i, c in enumerate(mean_piv.index):
    for j, s in enumerate(STAGES):
        m = vals[i, j]
        f = frac_piv.iloc[i, j]
        color = "white" if m < vals.max() * 0.55 else "black"
        ax.text(j, i, f"{m:.2f}\n({f*100:.0f}%)", ha="center", va="center",
                fontsize=7, color=color)

ax.set_xticks(range(len(STAGES)))
ax.set_xticklabels(STAGES, fontsize=8)
ax.set_yticks(range(len(mean_piv.index)))
ax.set_yticklabels(mean_piv.index, fontsize=8)
ax.tick_params(length=0)
for spine in ax.spines.values():
    spine.set_visible(False)

cb = fig.colorbar(im, ax=ax, shrink=0.55, pad=0.03)
cb.set_label("Mean IDO1", fontsize=8)
cb.ax.tick_params(labelsize=7)
cb.outline.set_visible(False)

fig.tight_layout()
fig.savefig(OUT + ".svg")
fig.savefig(OUT + ".png", dpi=300)
fig.savefig(OUT + ".pdf")
print("saved:", OUT + ".{svg,png,pdf,csv}")
