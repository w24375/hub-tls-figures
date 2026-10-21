"""Fig 8b: 'metabolic brake' model schematic v2 — fully vector, editable text (SVG)."""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

C_BRAKE = '#C41E3A'; C_BLUE = '#6C8EBF'; C_HUB = '#9467BD'; C_TLS = '#5DA05D'; C_GREY = '#666666'

fig, ax = plt.subplots(figsize=(9.2, 4.4))
ax.set_xlim(0, 10); ax.set_ylim(0, 5); ax.axis('off')

# ---------- LEFT: stalled hub ----------
ax.text(2.2, 4.7, 'Stalled immune hub', ha='center', fontsize=10.5, fontweight='bold', color=C_HUB)
ax.add_patch(Circle((1.7, 2.75), 1.45, facecolor=C_HUB, alpha=0.12, edgecolor=C_HUB, linewidth=1.2, linestyle='--', zorder=1))
ax.add_patch(Circle((1.7, 2.75), 0.58, facecolor=C_HUB, alpha=0.55, edgecolor='none', zorder=2))
ax.text(1.7, 2.75, 'CCL19\nCCL21\nCXCL13', ha='center', va='center', fontsize=6.3, color='white', zorder=3)
# mregDC + IDO1
ax.add_patch(Circle((0.72, 3.9), 0.30, facecolor=C_BRAKE, alpha=0.9, edgecolor='white', linewidth=1, zorder=4))
ax.text(0.72, 3.9, 'IDO1', ha='center', va='center', fontsize=6.5, color='white', fontweight='bold', zorder=5)
ax.text(0.72, 4.32, 'mregDC', ha='center', fontsize=7, color=C_BRAKE)
# kynurenine arrow down to Tfh/B
ax.add_patch(FancyArrowPatch((0.95, 3.62), (2.15, 1.55), arrowstyle='-|>', mutation_scale=11,
                             color=C_BRAKE, linewidth=1.4, zorder=4))
ax.text(2.3, 1.62, 'Kynurenine\n(Trp depletion)', fontsize=6.5, color=C_BRAKE, ha='left', va='center')
# suppressed Tfh/B
ax.add_patch(Circle((2.15, 1.05), 0.26, facecolor=C_GREY, alpha=0.5, edgecolor='white', zorder=4))
ax.text(2.15, 0.55, 'Tfh / B cell\n(suppressed)', ha='center', fontsize=6.5, color=C_GREY)
# faded TLS + blocked arrow
ax.add_patch(Circle((3.7, 2.75), 0.72, facecolor=C_TLS, alpha=0.06, edgecolor=C_TLS, linewidth=1.0, linestyle='--', zorder=1))
ax.add_patch(FancyArrowPatch((2.75, 2.75), (3.35, 2.75), arrowstyle='-|>', mutation_scale=11,
                             color=C_GREY, linewidth=1.6, alpha=0.5, zorder=3, linestyle='--'))
ax.plot([3.52], [2.75], marker='x', color=C_BRAKE, markersize=9, markeredgewidth=2.4, zorder=5)
ax.text(3.15, 3.3, 'maturation\nblocked', ha='center', fontsize=6.3, color=C_BRAKE)
ax.text(3.7, 1.85, 'TLS fails\nto mature', ha='center', fontsize=6.3, color=C_GREY)
ax.text(2.2, 0.06, 'Non-responder pattern: IDO1+ hub far from mature TLS',
        ha='center', fontsize=7.5, color=C_BRAKE, style='italic')

# ---------- MIDDLE: intervention arrow ----------
ax.add_patch(FancyArrowPatch((4.62, 2.75), (5.5, 2.75), arrowstyle='-|>', mutation_scale=16,
                             color='black', linewidth=2.2, zorder=5))
ax.text(5.06, 3.12, 'IDO1 blockade', ha='center', fontsize=8.5, fontweight='bold')
ax.text(5.06, 2.38, '± anti-PD-1', ha='center', fontsize=7, color=C_GREY)

# ---------- RIGHT: brake released ----------
ax.text(8.3, 4.7, 'Mature TLS (SFL)', ha='center', fontsize=10.5, fontweight='bold', color=C_TLS)
ax.add_patch(Circle((6.75, 2.75), 1.0, facecolor=C_HUB, alpha=0.05, edgecolor=C_HUB, linewidth=1.0, linestyle='--', zorder=1))
ax.text(6.75, 1.55, 'hub\n(resolving)', ha='center', fontsize=6.3, color=C_HUB)
ax.add_patch(Circle((8.6, 2.75), 1.25, facecolor=C_TLS, alpha=0.14, edgecolor=C_TLS, linewidth=1.4, zorder=1))
ax.add_patch(Circle((8.6, 2.75), 0.70, facecolor=C_TLS, alpha=0.4, edgecolor='none', zorder=2))
ax.text(8.6, 2.75, 'GC\nB cells', ha='center', va='center', fontsize=7, color='white', zorder=3)
# conversion arrow
ax.add_patch(FancyArrowPatch((7.35, 2.75), (8.0, 2.75), arrowstyle='-|>', mutation_scale=13,
                             color=C_TLS, linewidth=2.4, zorder=4))
ax.text(7.68, 3.05, 'conversion', ha='center', fontsize=6.5, color=C_TLS)
# active Tfh/B
ax.add_patch(Circle((6.35, 1.0), 0.26, facecolor=C_BLUE, alpha=0.9, edgecolor='white', zorder=4))
ax.text(6.35, 0.5, 'Tfh / B cell\n(active)', ha='center', fontsize=6.5, color=C_BLUE)
ax.add_patch(FancyArrowPatch((6.6, 1.3), (8.05, 2.3), arrowstyle='-|>', mutation_scale=10,
                             color=C_BLUE, linewidth=1.2, zorder=4))
ax.text(7.1, 1.05, 'Tfh help\nrestored', fontsize=6.5, color=C_BLUE, ha='center')
ax.text(8.3, 0.06, 'Responder pattern: IDO1+ hub adjacent to mature TLS',
        ha='center', fontsize=7.5, color=C_TLS, style='italic')

plt.savefig('data/hub_tls_paper_outline/figure_ido1_assembly/panels_svg/fig8b_brake_model.svg', bbox_inches='tight')
plt.savefig('data/hub_tls_paper_outline/figure_ido1_assembly/panels_svg/fig8b_brake_model.png', dpi=300, bbox_inches='tight')
print('schematic v2 saved')
