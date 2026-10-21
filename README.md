# Code — Immune hubs drive directional maturation of tertiary lymphoid structures under an IDO1 metabolic brake

Panel-plotting code for the manuscript's main figures (Figures 1–8).

## Scope

Sixteen scripts, at most two per figure, covering the key panels of each main figure:

| Directory | Scripts | Panels covered |
|---|---|---|
| `figure1/` | `assemble_figure1.R`, `generate_panel_c.R` | radial gradients, proximal–distal comparison, cell-type composition, producer cells |
| `figure2/` | `assemble_figure2.R`, `subtype_robustness_check.R` | Base/Stem hub directionality; TLS-shared-gene-removal robustness check |
| `figure3/` | `render_panel_a_axis_tiles.R`, `render_panel_b_strip.R`, `export_slingshot_paired_slope.R` | axis tiles, pseudotime strip, per-sample Slingshot ordering |
| `figure4/` | `assemble_figure4.R`, `virtual_knockout_gallery.py` | virtual-knockout rankings, LOO de-circularization, PCA displacement field, spatial perturbation maps |
| `figure5/` | `assemble_figure5.R`, `panel_g_sample_level_standalone.R` | causal mediation, ligand–receptor spatial proximity, hub-proximal versus hub-distal TLS |
| `figure6/` | `plot_immunotherapy_bars.R` | Cohort-B outcome bars and hub→TLS→response flow |
| `figure7/` | `make_fig7B_heatmap.py`, `make_fig7C_panels.py`, `make_ido1_lineage_split_figures.py` | IDO1 decline along the maturation axis, within-sample correlation distributions, dual-trajectory pseudotime |
| `figure8/` | `make_brake_schematic.py` | IDO1 metabolic-brake model schematic |

This repository provides panel-level plotting code for the main figures. Preprocessing pipelines, computational-analysis scripts (spatial QC, hub/TLS calling, pseudotime inference, virtual knockout, causal mediation, deconvolution, IDO1 screening), and supplementary-figure scripts are available from the corresponding authors upon request.

## Dependencies

- **R** ≥ 4.3: Seurat, ggplot2, patchwork, svglite, ragg, cowplot, RANN, dplyr, tidyr
- **Python** ≥ 3.10: matplotlib, pandas, numpy

## Notes for reuse

- `code/main.sh` runs all 16 scripts in figure order. Each script is independent and can also be run individually.
- All input paths are **relative to a `data/` root** (e.g. `data/hub_tls_paper_outline/figure1_assembled/…`). Unpack the manuscript's source data under `data/` preserving this layout, or set your working directory accordingly.
- Figure conventions: 183 mm double-column width, Arial throughout, RGB color mode; SVGs are delivered without `textLength` attributes for Illustrator compatibility.

## Contact

Corresponding authors (see manuscript). Issues with figure rendering: open an issue with the figure number and panel letter.
