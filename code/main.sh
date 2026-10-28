#!/usr/bin/env bash
# Runs all 16 figure scripts in order. Each one also works standalone, e.g.:
#   Rscript code/figure1/assemble_figure1.R
#   python3 code/figure4/virtual_knockout_gallery.py
#
# Needs the manuscript's source data under ./data/ (see README for layout).

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
