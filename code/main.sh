#!/usr/bin/env bash
# main.sh — runs all 16 panel-plotting scripts in figure order.
#
# Prerequisites:
#   1. The manuscript's source data must be placed under ./data/ preserving
#      the layout described in README.md (e.g. data/hub_tls_paper_outline/figure1_assembled/…).
#   2. R (>= 4.3) with packages: Seurat, ggplot2, patchwork, svglite, ragg,
#      cowplot, RANN, dplyr, tidyr, data.table.
#   3. Python (>= 3.10) with packages: matplotlib, pandas, numpy, scipy,
#      scikit-learn.
#
# Each script is independent and can also be run individually, e.g.:
#   Rscript code/figure1/assemble_figure1.R
#   python3 code/figure4/virtual_knockout_gallery.py

set -euo pipefail
cd "$(dirname "$0")/.."

# Figure 1
Rscript code/figure1/assemble_figure1.R
Rscript code/figure1/generate_panel_c.R

# Figure 2
Rscript code/figure2/assemble_figure2.R
Rscript code/figure2/subtype_robustness_check.R

# Figure 3
Rscript code/figure3/render_panel_a_axis_tiles.R
Rscript code/figure3/render_panel_b_strip.R
Rscript code/figure3/export_slingshot_paired_slope.R

# Figure 4
Rscript code/figure4/assemble_figure4.R
python3 code/figure4/virtual_knockout_gallery.py

# Figure 5
Rscript code/figure5/assemble_figure5.R
Rscript code/figure5/panel_g_sample_level_standalone.R

# Figure 6
Rscript code/figure6/plot_immunotherapy_bars.R

# Figure 7
python3 code/figure7/make_fig7B_heatmap.py
python3 code/figure7/make_fig7C_panels.py
python3 code/figure7/make_ido1_lineage_split_figures.py

# Figure 8
python3 code/figure8/make_brake_schematic.py

echo "All figure scripts completed."
