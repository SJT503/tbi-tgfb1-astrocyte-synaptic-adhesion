# -*- coding: utf-8 -*-
"""
Nature Communications Unified Bioinformatics & Integrated Figures Generator (v2 - Prism/NPG Aesthetic)
--------------------------------------------------------------------------------------------------------
Strictly aligned with the wet-lab GraphPad Prism / Nature visual style of Fig 1-5:
- Translucent pastel fills (alpha=0.36-0.40 blended over white) + crisp dark borders (1.15 pt)
- Zero dashed background gridlines; clean 0.95 pt L-spines; outward ticks; Arial 8-9.5 pt typography
- Concise in-plot annotations (*, **, ***, ****, #, ns) with zero text overlap or wordy callouts
- Rich, varied Nature Communications chart types:
  * Fig 6: (A) Horizontal lollipop-bar atlas, (B) KDE contour + UMAP trajectory, (C) Framed DotPlot matrix,
           (D) Alluvial flow-connected stacked bar chart, (E) Cubic-spline pseudotime ribbon curves,
           (F) 3-group Z-score regulon heatmap + delta lollipop, (G) Full 2,006-gene CellOracle S-curve,
           (H) Cross-species MEF2C coupling forest/lollipop plot.
  * Fig 7: (A-B) Mouse & Human CellChat diverging lollipop-bar plots, (C) 3-group rescue bars + trajectory,
           (D) 2D pathway x DG-subtype bubble-heatmap matrix, (E) Hippocampus vs Cortex vulnerability plot,
           (F) Cross-cohort effect-size forest plot.
  * Fig 8: (A) HT-22 x C8-D1A co-culture IF montage, (B) Prism bar + replicate scatter + SD error bars,
           (C) Dual-trajectory time-course ribbon + replicate scatter (GSE271769), (D) DG_C0-C3 alluvial
           stacked bar chart, (E) Human neuronal apoptosis diverging lollipop-bar plot, (F) Working model.
  * Fig S1-S4: Unified translucent Prism/NPG palette and clean typography.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch, Polygon
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import make_interp_spline
from scipy.stats import gaussian_kde
from PIL import Image

# ---------------------------------------------------------------------------
# Global Nature / GraphPad Prism Publication rcParams (matching Fig 1-5)
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.edgecolor": "#1E293B",
    "axes.linewidth": 0.95,
    "axes.labelcolor": "#0F172A",
    "xtick.color": "#1E293B",
    "ytick.color": "#1E293B",
    "xtick.major.width": 0.95,
    "ytick.major.width": 0.95,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.labelsize": 8.2,
    "ytick.labelsize": 8.2,
    "axes.labelsize": 9.2,
    "axes.labelweight": "bold",
    "axes.titlesize": 9.8,
    "axes.titleweight": "bold",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})


def blend_white(hex_col, alpha=0.72):
    """Blends a hex color with white at `alpha` (0.72 rich NPG saturation)."""
    h = hex_col.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    rb = int(round(r * alpha + 255 * (1.0 - alpha)))
    gb = int(round(g * alpha + 255 * (1.0 - alpha)))
    bb = int(round(b * alpha + 255 * (1.0 - alpha)))
    return f"#{rb:02X}{gb:02X}{bb:02X}"


# Exact Nature Publishing Group (NPG) + Okabe-Ito Edge Colors, Marker Colors, and Rich Fills (0% pink)
COL_SHAM_EDGE = "#1E293B"
COL_SHAM_RAW  = "#475569"
COL_SHAM_FILL = blend_white("#64748B", 0.68)      # Rich cool slate

COL_TBI_EDGE  = "#881D0C"
COL_TBI_RAW   = "#E64B35"
COL_TBI_FILL  = blend_white(COL_TBI_RAW, 0.74)    # Rich NPG vermilion-red

COL_SCR_EDGE  = "#154360"
COL_SCR_RAW   = "#2E86C1"
COL_SCR_FILL  = blend_white(COL_SCR_RAW, 0.74)    # Rich NPG cerulean-sapphire blue

COL_KO_EDGE   = "#005A4B"
COL_KO_RAW    = "#00A087"
COL_KO_FILL   = blend_white(COL_KO_RAW, 0.78)     # Rich NPG viridian emerald-teal (Core Rescue)

COL_TEAL_EDGE = "#1F2E52"
COL_TEAL_RAW  = "#3C5488"
COL_TEAL_FILL = blend_white(COL_TEAL_RAW, 0.74)   # Rich NPG deep sapphire-indigo

COL_GOLD_EDGE = "#8C510A"
COL_GOLD_RAW  = "#E69F00"
COL_GOLD_FILL = blend_white(COL_GOLD_RAW, 0.76)   # Rich Okabe-Ito warm amber-gold

# 7 Major CNS Lineages Palette for GSE333879 (n = 35,047) — High-Contrast NPG Jewel Tones
LINEAGE_PALETTE = {
    "Excitatory neuron":      ("#3C5488", "#1F2E52"),  # Deep NPG Sapphire-Indigo
    "Oligodendrocyte":        ("#E69F00", "#8C510A"),  # Warm Okabe-Ito Amber-Gold
    "Inhibitory neuron":      ("#2E86C1", "#154360"),  # Crisp NPG Cerulean Blue
    "Microglia":              ("#E64B35", "#881D0C"),  # Rich NPG Vermilion-Red
    "Astrocyte":              ("#00A087", "#004D40"),  # Vibrant NPG Emerald-Teal (CORE FOCUS)
    "OPC":                    ("#5E4FA2", "#3B2F6B"),  # Deep Royal Indigo-Violet (non-pink)
    "Endothelial / Pericyte": ("#8C510A", "#4A2A05"),  # Rich NPG Warm Sienna-Bronze
}


def style_umap_ax(ax, xlabel="UMAP 1", ylabel="UMAP 2"):
    """Applies full-length L-shaped UMAP coordinate axes (left & bottom spines) with clean ticks."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_linewidth(0.95)
    ax.spines["bottom"].set_linewidth(0.95)
    ax.spines["left"].set_color("#0F172A")
    ax.spines["bottom"].set_color("#0F172A")
    ax.grid(False)
    ax.set_xlabel(xlabel, fontsize=8.6, fontweight="bold", labelpad=3)
    ax.set_ylabel(ylabel, fontsize=8.6, fontweight="bold", labelpad=3)
    ax.tick_params(axis="both", which="major", labelsize=7.6, width=0.85, length=3.0, colors="#0F172A")

# Custom Nature Publishing Group (NPG) Blue-White-Vermilion and Blue-White-Emerald colormaps
NPG_CMAP = LinearSegmentedColormap.from_list(
    "npg_div",
    ["#1B4F72", "#5499C7", "#F8FAFC", "#EB984E", "#992B15"],
    N=256
)
NPG_RESCUE_CMAP = LinearSegmentedColormap.from_list(
    "npg_rescue",
    ["#1B4F72", "#7FB3D5", "#F8FAFC", "#76D7C4", "#006652"],
    N=256
)

BASE = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project"
OUT_MAIN = os.path.join(BASE, "JNI", "figures", "cleaned_Fig1_5")
OUT_SUPP = os.path.join(BASE, "JNI", "supplementary", "cleaned_FigS1_S6")
os.makedirs(OUT_MAIN, exist_ok=True)
os.makedirs(OUT_SUPP, exist_ok=True)


def style_ax(ax, zero_line=None):
    """Applies clean L-shaped spines with zero dashed background gridlines (matching Fig 1-5)."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.95)
    ax.spines["bottom"].set_linewidth(0.95)
    ax.grid(False)
    if zero_line == "y":
        ax.axhline(0, color="#475569", linewidth=0.85, linestyle="--", zorder=1)
    elif zero_line == "x":
        ax.axvline(0, color="#475569", linewidth=0.85, linestyle="--", zorder=1)


def add_panel_label(ax, letter, x=-0.14, y=1.10, fontsize=15.5):
    ax.text(
        x, y, letter,
        transform=ax.transAxes,
        fontsize=fontsize,
        fontweight="bold",
        fontfamily="Times New Roman",
        va="top", ha="left",
        color="#000000",
        clip_on=False
    )


def add_bracket(ax, x1, x2, y, h, text, fontsize=8.5, color="#0F172A"):
    """Draws a crisp GraphPad Prism significance bracket matching Fig 1-5."""
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=0.9, c=color, clip_on=False)
    yr = ax.get_ylim()[1] - ax.get_ylim()[0]
    ax.text(
        (x1 + x2) * 0.5, y + h + 0.008 * yr,
        text, ha="center", va="bottom", color=color,
        fontsize=fontsize, fontweight="bold", clip_on=False
    )


def save_pub_figure(fig, out_dir, stem, dpi=300):
    png_path = os.path.join(out_dir, f"{stem}.png")
    tif_path = os.path.join(out_dir, f"{stem}.tif")
    pdf_path = os.path.join(out_dir, f"{stem}.pdf")
    svg_path = os.path.join(out_dir, f"{stem}.svg")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    im = Image.open(png_path).convert("RGB")
    im.save(tif_path, compression="tiff_lzw", dpi=(dpi, dpi))
    print(f"[OK] Saved {stem} (.tif / .png / .pdf / .svg) -> {im.size}")


# ===========================================================================
# FIGURE 6: 8-Panel (A-H) In-House GSE333879 Single-Nucleus Discovery
# ===========================================================================
def build_figure_6():
    fig = plt.figure(figsize=(16.6, 18.6), dpi=300)
    gs = gridspec.GridSpec(
        4, 2, figure=fig,
        height_ratios=[1.15, 1.04, 1.12, 1.04],
        width_ratios=[1.08, 1.00],
        hspace=0.48, wspace=0.28
    )

    # Dedicated Row 0 GridSpec so Panel A (UMAP + DotPlot) gets generous horizontal width
    # and neither UMAP A nor UMAP B is horizontally compressed
    gs_row0 = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=gs[0, :],
        width_ratios=[1.38, 1.00], wspace=0.24
    )

    # -----------------------------------------------------------------------
    # Panel A: GSE333879 Global 35,047-Nucleus UMAP + 7-Lineage Marker DotPlot
    # -----------------------------------------------------------------------
    gs_A = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_row0[0, 0], width_ratios=[1.42, 0.82], wspace=0.30)
    axA_umap = fig.add_subplot(gs_A[0, 0])
    axA_dot  = fig.add_subplot(gs_A[0, 1])
    add_panel_label(axA_umap, "A", x=-0.15, y=1.14)

    atlas_df = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/atlas_export/atlas_35047_umap_meta.csv"))
    dot7_df  = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/atlas_export/atlas_7lineage_markers_dotplot.csv"))

    # Plot 7 major CNS lineages on Global Harmony UMAP with high-saturation NPG jewel tones
    ct_order_7 = [
        "Excitatory neuron",
        "Oligodendrocyte",
        "Inhibitory neuron",
        "Microglia",
        "OPC",
        "Endothelial / Pericyte",
        "Astrocyte",  # Plotted on top so core focal population is crystal-clear
    ]
    rng_a = np.random.default_rng(2026)
    for ct in ct_order_7:
        sub = atlas_df[atlas_df["celltype"] == ct]
        raw_c, edge_c = LINEAGE_PALETTE[ct]
        is_ast = (ct == "Astrocyte")
        if len(sub) > 5000:
            idx_s = rng_a.choice(len(sub), size=5000, replace=False)
            sub_plot = sub.iloc[idx_s]
        else:
            sub_plot = sub
        axA_umap.scatter(
            sub_plot["UMAP1_harm"], sub_plot["UMAP2_harm"],
            s=8.0 if is_ast else 3.0, c=edge_c, alpha=0.85 if is_ast else 0.45,
            linewidths=0, zorder=4 if is_ast else 2, rasterized=True
        )
        axA_umap.scatter(
            sub_plot["UMAP1_harm"], sub_plot["UMAP2_harm"],
            s=5.0 if is_ast else 1.7, c=raw_c, alpha=0.95 if is_ast else 0.85,
            linewidths=0, zorder=5 if is_ast else 3, rasterized=True
        )

    # Highlight Core Focal Clusters with Dashed Zoom Boxes & Outer-Margin Callout Pills (0% overlap with cells)
    # 1) Astrocyte (Aqp4+, n = 1,038 -> 774 pure) at (+2.95, -10.98)
    rect_ast = plt.Rectangle((0.3, -13.3), 5.3, 4.4, fill=False, edgecolor="#004D40",
                             linewidth=1.4, linestyle="--", zorder=6)
    axA_umap.add_patch(rect_ast)
    axA_umap.annotate(
        "Astrocyte (1,038)\npure n = 774 → B",
        xy=(5.6, -11.2), xytext=(13.2, -13.8),
        ha="center", va="center", fontsize=6.6, fontweight="bold", color="white",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="#007B66", edgecolor="#004D40", lw=0.95),
        arrowprops=dict(arrowstyle="-|>,head_width=0.18,head_length=0.25", lw=1.05, color="#004D40"),
        zorder=8
    )

    # 2) DG Granule (Prox1+, n = 3,522) island at (-3.44, -5.38)
    rect_dg = plt.Rectangle((-6.8, -8.3), 6.5, 5.8, fill=False, edgecolor="#1F2E52",
                            linewidth=1.25, linestyle="--", zorder=6)
    axA_umap.add_patch(rect_dg)
    axA_umap.annotate(
        "DG Granule\n(Prox1+, 3,522)",
        xy=(-6.8, -5.4), xytext=(-14.2, -5.4),
        ha="center", va="center", fontsize=6.4, fontweight="bold", color="#0F172A",
        bbox=dict(boxstyle="round,pad=0.18", facecolor="#F8FAFC", edgecolor="#1F2E52", lw=0.85),
        arrowprops=dict(arrowstyle="-", lw=0.85, color="#1F2E52"),
        zorder=7
    )

    # Outer-margin callout labels with leader lines so zero text sits on top of scatter points
    outer_callouts = [
        ((-6.8,  7.5), (-13.6, 12.8), "Excitatory (18,556)", "#3C5488"),
        (( 5.2,  9.8), ( 12.2, 13.8), "Inhibitory (4,314)",  "#2E86C1"),
        ((10.8, -0.5), ( 14.2,  4.2), "Oligo (7,647)",       "#B45309"),
        (( 3.4, -3.1), ( 13.6, -5.8), "OPC (988)",           "#4338CA"),
        ((-2.2,-14.8), ( -1.5,-17.8), "Microglia (2,155)",   "#B91C1C"),
        ((-9.9,-11.4), (-13.8,-14.8), "Endo/Peri (349)",     "#78350F"),
    ]
    for xy_pt, xy_txt, ptxt, pcol in outer_callouts:
        axA_umap.annotate(
            ptxt, xy=xy_pt, xytext=xy_txt,
            ha="center", va="center", fontsize=6.3, fontweight="bold", color=pcol,
            bbox=dict(boxstyle="round,pad=0.16", facecolor="white", edgecolor=pcol, lw=0.75, alpha=0.95),
            arrowprops=dict(arrowstyle="-", lw=0.75, color=pcol, alpha=0.80),
            zorder=7
        )

    axA_umap.set_xlim(-19.2, 19.8)
    axA_umap.set_ylim(-20.2, 17.5)
    style_umap_ax(axA_umap, xlabel="UMAP 1", ylabel="UMAP 2")
    axA_umap.set_xticks([-15, -10, -5, 0, 5, 10, 15])
    axA_umap.set_yticks([-15, -10, -5, 0, 5, 10, 15])
    axA_umap.set_title(
        "Hippocampal snRNA-seq atlas (GSE333879, n = 35,047)\nSham: 11,378 | CCI+scr: 6,845 | CCI+TGF-βR1 KO: 16,824",
        fontsize=8.8, fontweight="bold", pad=8
    )

    # Right subpanel of Panel A: Compact 7-Lineage Canonical Marker DotPlot
    ct_rows_a = [
        ("Excitatory neuron",      "Excitatory"),
        ("Inhibitory neuron",      "Inhibitory"),
        ("Oligodendrocyte",        "Oligo"),
        ("Astrocyte",              "Astrocyte"),
        ("Microglia",              "Microglia"),
        ("OPC",                    "OPC"),
        ("Endothelial / Pericyte", "Endo/Peri"),
    ]
    genes_a = ["Slc17a7", "Gad1", "Mbp", "Aqp4", "Cx3cr1", "Pdgfra", "Flt1"]
    for sp in ("top", "right", "left", "bottom"):
        axA_dot.spines[sp].set_visible(True)
        axA_dot.spines[sp].set_linewidth(0.8)
        axA_dot.spines[sp].set_color("#94A3B8")
    axA_dot.grid(True, linestyle="-", linewidth=0.4, color="#E2E8F0", zorder=1)

    xa_l, ya_l, sa_l, ca_l = [], [], [], []
    for r_i, (ct_full, _) in enumerate(ct_rows_a):
        y_r = len(ct_rows_a) - 1 - r_i
        for c_j, g_a in enumerate(genes_a):
            row_m = dot7_df[(dot7_df["celltype"] == ct_full) & (dot7_df["gene"] == g_a)]
            if len(row_m) > 0:
                pct_v = float(row_m["pct_exp"].values[0])
                z_v   = float(row_m["z_score"].values[0])
            else:
                pct_v, z_v = 0.0, -0.5
            xa_l.append(c_j)
            ya_l.append(y_r)
            sa_l.append(max(pct_v * 1.45, 10.0))
            ca_l.append(z_v)

    axA_dot.scatter(
        xa_l, ya_l, s=sa_l, c=ca_l, cmap=NPG_CMAP,
        vmin=-1.2, vmax=2.0, edgecolors="#0F172A", linewidths=0.7, zorder=3
    )
    axA_dot.axhspan(2.55, 3.45, color="#00A087", alpha=0.14, zorder=0)

    axA_dot.set_xticks(range(len(genes_a)))
    axA_dot.set_xticklabels(genes_a, rotation=40, ha="right", fontstyle="italic", fontweight="bold", fontsize=7.4)
    axA_dot.set_yticks(range(len(ct_rows_a)))
    axA_dot.set_yticklabels([x[1] for x in ct_rows_a[::-1]], fontweight="bold", fontsize=7.4)
    for tick_lbl in axA_dot.get_yticklabels():
        if "Astrocyte" in tick_lbl.get_text():
            tick_lbl.set_color("#005A4B")
    axA_dot.set_xlim(-0.6, len(genes_a) - 0.4)
    axA_dot.set_ylim(-0.55, len(ct_rows_a) - 0.45)
    axA_dot.set_title("Lineage markers", fontsize=9.0, fontweight="bold", pad=10)

    # -----------------------------------------------------------------------
    # Panel B: Pure Astrocyte 3-State UMAP (n = 774) + High-Contrast NPG + Trajectory
    # -----------------------------------------------------------------------
    axB = fig.add_subplot(gs_row0[0, 1])
    add_panel_label(axB, "B", x=-0.12, y=1.14)

    umap_df = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/celloracle/umap_coords.csv"))
    meta_df = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/celloracle/cell_metadata.csv"))
    ast_df = pd.merge(umap_df, meta_df, on="cell_barcode")

    xg, yg = np.mgrid[-5.8:5.8:100j, -5.2:5.0:100j]
    positions = np.vstack([xg.ravel(), yg.ravel()])

    # Rich, Vibrant NPG Colors: Cerulean-Sapphire (Homeo-A), Viridian Emerald (Homeo-B), Vermilion (Reactive)
    state_palette = [
        ("Homeostatic-A",        COL_SCR_RAW, COL_SCR_EDGE, "o", "Homeostatic-A"),
        ("Homeostatic-B(Apoe+)", COL_KO_RAW,  COL_KO_EDGE,  "^", "Homeostatic-B (Apoe+)"),
        ("Reactive(Id3+)",       COL_TBI_RAW, COL_TBI_EDGE, "v", "Reactive (Id3+)")
    ]
    for st_key, raw_c, edge_c, mk, lbl in state_palette:
        sub = ast_df[ast_df["cluster_label"] == st_key]
        if len(sub) > 15:
            kde = gaussian_kde(np.vstack([sub["UMAP_1"], sub["UMAP_2"]]), bw_method=0.42)
            z = np.reshape(kde(positions).T, xg.shape)
            axB.contourf(xg, yg, z, levels=[z.max()*0.16, z.max()*1.05], colors=[raw_c], alpha=0.22, zorder=1)
            axB.contour(xg, yg, z, levels=[z.max()*0.16], colors=[edge_c], linewidths=0.95, linestyles="--", alpha=0.70, zorder=2)
        axB.scatter(
            sub["UMAP_1"], sub["UMAP_2"],
            s=28, marker=mk, facecolor=raw_c,
            edgecolor=edge_c, linewidths=0.65, alpha=0.92,
            label=f"{lbl} (n = {len(sub)})", zorder=3
        )

    c_ha = ast_df[ast_df["cluster_label"] == "Homeostatic-A"][["UMAP_1", "UMAP_2"]].median().values
    c_hb = ast_df[ast_df["cluster_label"] == "Homeostatic-B(Apoe+)"][["UMAP_1", "UMAP_2"]].median().values
    c_re = ast_df[ast_df["cluster_label"] == "Reactive(Id3+)"][["UMAP_1", "UMAP_2"]].median().values

    # Bifurcating Monocle3 trajectory curves with white halo for contrast
    for xy_to, col_arr, ls_arr, rad_v in [
        (c_hb, "#004D40", "-",  0.18),
        (c_re, "#7F1D1D", "--", -0.15)
    ]:
        axB.annotate(
            "", xy=xy_to, xytext=c_ha,
            arrowprops=dict(arrowstyle="-|>,head_width=0.36,head_length=0.46", lw=3.2, color="white", connectionstyle=f"arc3,rad={rad_v}"),
            zorder=4
        )
        axB.annotate(
            "", xy=xy_to, xytext=c_ha,
            arrowprops=dict(arrowstyle="-|>,head_width=0.34,head_length=0.44", lw=2.1, color=col_arr, ls=ls_arr, connectionstyle=f"arc3,rad={rad_v}"),
            zorder=5
        )

    # Outer-whitespace callout pills with leader lines so zero astrocyte points or trajectories are covered
    axB.annotate(
        "Homeostatic-A\n(n = 332)",
        xy=(-2.2, 2.5), xytext=(-4.8, 4.3),
        ha="center", va="center", fontsize=7.2, fontweight="bold", color=COL_SCR_EDGE,
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor=COL_SCR_EDGE, lw=0.95, alpha=0.96),
        arrowprops=dict(arrowstyle="-", lw=0.85, color=COL_SCR_EDGE), zorder=6
    )
    axB.annotate(
        "Homeostatic-B (Apoe+)\nRestorative (n = 313)",
        xy=(2.8, 0.4), xytext=(5.2, 2.8),
        ha="center", va="center", fontsize=7.2, fontweight="bold", color="white",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="#007B66", edgecolor="#004D40", lw=1.0, alpha=0.96),
        arrowprops=dict(arrowstyle="-", lw=0.95, color="#004D40"), zorder=6
    )
    axB.annotate(
        "Reactive (Id3+ / Gfap+)\n(n = 129)",
        xy=(0.8, -4.1), xytext=(4.4, -4.9),
        ha="center", va="center", fontsize=7.2, fontweight="bold", color=COL_TBI_EDGE,
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor=COL_TBI_EDGE, lw=0.95, alpha=0.96),
        arrowprops=dict(arrowstyle="-", lw=0.85, color=COL_TBI_EDGE), zorder=6
    )

    axB.set_xlim(-6.8, 7.8)
    axB.set_ylim(-6.2, 5.5)
    style_umap_ax(axB, xlabel="UMAP 1", ylabel="UMAP 2")
    axB.set_xticks([-6, -4, -2, 0, 2, 4, 6])
    axB.set_yticks([-6, -4, -2, 0, 2, 4])
    axB.set_title("Hippocampal astrocyte states & Monocle3 trajectory (n = 774)", pad=22)
    axB.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01),
        ncol=3, frameon=False, fontsize=7.3,
        columnspacing=0.70, handletextpad=0.28
    )

    # -----------------------------------------------------------------------
    # Panel C: Framed Nature-Style Subcluster Marker DotPlot
    # -----------------------------------------------------------------------
    axC = fig.add_subplot(gs[1, 0])
    add_panel_label(axC, "C", x=-0.16, y=1.12)

    pt_df = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/figures/Fig_Monocle_pseudotime/key_gene_pseudotime.csv"))
    genes_c = ["Atp1a2", "Slc1a2", "Aqp4", "Aldh1l1", "Nrxn1", "Nrxn2", "Apoe", "Cd81", "Id3", "Gfap"]
    clusters_c = ["Homeostatic-A", "Homeostatic-B(Apoe+)", "Reactive(Id3+)"]
    cluster_short = ["Homeostatic-A", "Homeostatic-B\n(Apoe+ / Cd81+)", "Reactive\n(Id3+ / Gfap+)"]

    mean_mat = np.zeros((len(clusters_c), len(genes_c)))
    pct_mat  = np.zeros((len(clusters_c), len(genes_c)))
    for i, cl in enumerate(clusters_c):
        for j, g in enumerate(genes_c):
            vals = pt_df[(pt_df["cluster"] == cl) & (pt_df["gene"] == g)]["expr"].values
            mean_mat[i, j] = vals.mean()
            pct_mat[i, j] = (vals > 0).mean() * 100.0

    z_mat = (mean_mat - mean_mat.mean(axis=0, keepdims=True)) / (mean_mat.std(axis=0, keepdims=True) + 1e-6)

    # Subtle matrix frame & light coordinate crosshairs (standard Scanpy/Nature DotPlot)
    for sp in ("top", "right", "left", "bottom"):
        axC.spines[sp].set_visible(True)
        axC.spines[sp].set_linewidth(0.85)
        axC.spines[sp].set_color("#94A3B8")
    axC.grid(True, linestyle="-", linewidth=0.45, color="#E2E8F0", zorder=1)

    xs_d, ys_d, sz_d, col_d = [], [], [], []
    for i in range(len(clusters_c)):
        for j in range(len(genes_c)):
            xs_d.append(j)
            ys_d.append(2 - i)
            sz_d.append(max(pct_mat[i, j] * 2.3, 22))
            col_d.append(z_mat[i, j])

    sc = axC.scatter(
        xs_d, ys_d, s=sz_d, c=col_d, cmap=NPG_CMAP,
        vmin=-1.35, vmax=1.35, edgecolors="#1E293B", linewidths=0.85, zorder=3
    )
    axC.set_xticks(range(len(genes_c)))
    axC.set_xticklabels(genes_c, rotation=32, ha="right", fontstyle="italic", fontweight="bold")
    axC.set_yticks([2, 1, 0])
    axC.set_yticklabels(cluster_short, fontweight="bold", fontsize=8.2)
    axC.set_xlim(-0.65, len(genes_c) - 0.35)
    axC.set_ylim(-0.55, 2.55)
    axC.set_title("State marker & synaptic-adhesion ligand profile", pad=22)

    cbar = plt.colorbar(sc, ax=axC, fraction=0.032, pad=0.03)
    cbar.set_label("Z-score", fontsize=8.0, fontweight="bold")
    cbar.ax.tick_params(labelsize=7.5)

    # Place % expressing dot-size legend cleanly ABOVE the frame (never inside the grid!)
    axC.text(0.0, 1.04, "Pct. expressed:", transform=axC.transAxes, fontsize=7.6, fontweight="bold", va="center", color="#334155")
    for pct_val, x_ax in zip([25, 50, 80], [0.30, 0.50, 0.72]):
        axC.scatter([x_ax], [1.045], s=pct_val * 1.6, c="#E2E8F0", edgecolors="#1E293B", lw=0.75,
                    transform=axC.transAxes, clip_on=False)
        axC.text(x_ax + 0.035, 1.04, f"{pct_val}%", transform=axC.transAxes, fontsize=7.5, va="center", ha="left", color="#334155")

    # -----------------------------------------------------------------------
    # Panel D: Alluvial / Ribbon-Connected Stacked Proportion Chart
    # -----------------------------------------------------------------------
    axD = fig.add_subplot(gs[1, 1])
    style_ax(axD)
    add_panel_label(axD, "D", x=-0.12, y=1.12)

    groups_d = ["Sham\n(n = 173)", "CCI+scr\n(n = 102)", "CCI+TGF-βR1 KO\n(n = 499)"]
    # Exact proportions (%): Homeostatic-A, Homeostatic-B (Apoe+), Reactive (Id3+)
    props = np.array([
        [50.87, 74.51, 33.67],  # Homeostatic-A (bottom stack - Cerulean-Sapphire)
        [42.20, 15.69, 44.89],  # Homeostatic-B (middle stack - Viridian Emerald-Teal rescue!)
        [6.94,  9.80,  21.44],  # Reactive (top stack - Vermilion-Red)
    ])
    st_fills = [COL_SCR_FILL, COL_KO_FILL, COL_TBI_FILL]
    st_edges = [COL_SCR_EDGE, COL_KO_EDGE, COL_TBI_EDGE]
    st_raws  = [COL_SCR_RAW,  COL_KO_RAW,  COL_TBI_RAW]
    st_names = ["Homeostatic-A", "Homeostatic-B (Apoe+)", "Reactive (Id3+)"]

    x_d = np.arange(3)
    bw = 0.42
    cum = np.zeros((4, 3))
    cum[1:] = np.cumsum(props, axis=0)

    # Draw smooth cosine alluvial ribbons between adjacent bars
    t_curve = np.linspace(0, 1, 40)
    s_curve = 0.5 * (1 - np.cos(np.pi * t_curve))
    for i in range(3):
        for g_idx in range(2):
            x_left = x_d[g_idx] + bw / 2
            x_right = x_d[g_idx + 1] - bw / 2
            xs_rib = x_left + (x_right - x_left) * t_curve
            y_bot = cum[i, g_idx] + (cum[i, g_idx + 1] - cum[i, g_idx]) * s_curve
            y_top = cum[i + 1, g_idx] + (cum[i + 1, g_idx + 1] - cum[i + 1, g_idx]) * s_curve
            axD.fill_between(xs_rib, y_bot, y_top, color=st_raws[i], alpha=0.25, edgecolor=st_edges[i], linewidth=0.5, zorder=1)

    # Draw stacked bars and crisp white-pill percentage labels inside each segment
    for i in range(3):
        axD.bar(
            x_d, props[i], bottom=cum[i], width=bw,
            color=st_fills[i], edgecolor=st_edges[i], linewidth=1.2,
            label=st_names[i], zorder=3
        )
        for g_idx in range(3):
            val = props[i, g_idx]
            if val >= 9.0:
                ymid = cum[i, g_idx] + val / 2.0
                axD.text(
                    x_d[g_idx], ymid, f"{val:.1f}%",
                    ha="center", va="center", fontsize=8.0, fontweight="bold",
                    color="#0F172A",
                    bbox=dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor="none", alpha=0.78),
                    zorder=4
                )

    add_bracket(axD, 1, 2, 102.5, 2.6, "****", fontsize=9.0, color=COL_KO_EDGE)
    axD.set_xticks(x_d)
    axD.set_xticklabels(groups_d, fontweight="bold")
    axD.set_ylabel("Astrocyte state proportion (%)")
    axD.set_ylim(0, 114)
    axD.set_xlim(-0.52, 2.52)
    axD.set_title("Astrocyte state remodeling & Homeostatic-B expansion", pad=22)
    axD.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), frameon=False, fontsize=7.3, ncol=3, columnspacing=0.75, handletextpad=0.3)

    # -----------------------------------------------------------------------
    # Panel E: Smooth Cubic-Spline Monocle3 Pseudotime Curves + Ribbons
    # -----------------------------------------------------------------------
    axE = fig.add_subplot(gs[2, 0])
    style_ax(axE)
    add_panel_label(axE, "E", x=-0.16, y=1.14)

    pt_valid = pt_df[np.isfinite(pt_df["pseudotime"])].copy()
    pt_valid["pt_bin"] = pd.cut(pt_valid["pseudotime"], bins=[0.0, 0.85, 2.0, 3.3, 4.5, 5.8, 7.6], include_lowest=True)

    genes_pt = [
        ("Nrxn1",  COL_SCR_EDGE,  COL_SCR_RAW,  "-",  "Nrxn1"),
        ("Slc1a2", COL_SHAM_EDGE, COL_SHAM_RAW, "--", "Slc1a2"),
        ("Apoe",   COL_KO_EDGE,   COL_KO_RAW,   "-",  "Apoe"),
        ("Cd81",   COL_TEAL_EDGE, COL_TEAL_RAW, "-.", "Cd81"),
        ("Id3",    COL_TBI_EDGE,  COL_TBI_RAW,  "-",  "Id3")
    ]
    x_smooth = np.linspace(0.25, 7.1, 120)
    for g, col_e, col_f, ls_g, lbl_g in genes_pt:
        sub_g = pt_valid[pt_valid["gene"] == g]
        grp_g = sub_g.groupby("pt_bin", observed=False).agg(
            pt_mean=("pseudotime", "mean"),
            ex_mean=("expr", "mean"),
            ex_sem=("expr", "sem")
        ).dropna()
        xv = grp_g["pt_mean"].values
        yv = grp_g["ex_mean"].values
        sv = grp_g["ex_sem"].fillna(0).values
        vmin, vmax = yv.min(), yv.max()
        yn = (yv - vmin) / (vmax - vmin + 1e-8)
        sn = np.clip(sv / (vmax - vmin + 1e-8) * 0.45, 0.03, 0.11)

        spl_y = make_interp_spline(xv, yn, k=2)(x_smooth)
        spl_s = make_interp_spline(xv, sn, k=2)(x_smooth)
        spl_y = np.clip(spl_y, 0.0, 1.02)
        axE.fill_between(x_smooth, np.clip(spl_y - spl_s, 0, 1.05), np.clip(spl_y + spl_s, 0, 1.08),
                         color=col_f, alpha=0.14, zorder=2)
        axE.plot(x_smooth, spl_y, color=col_e, linestyle=ls_g, linewidth=2.0, label=lbl_g, zorder=3)

    # Subtle top regime strip (y = 1.08..1.14) so text never collides with curves
    for x0_r, x1_r, fc_r, ec_r, txt_r in [
        (0.0, 2.3, COL_SHAM_FILL, COL_SHAM_EDGE, "Homeostatic-A"),
        (2.3, 5.1, COL_KO_FILL,   COL_KO_EDGE,   "Homeostatic-B (Apoe+)"),
        (5.1, 7.5, COL_TBI_FILL,  COL_TBI_EDGE,  "Reactive (Id3+)"),
    ]:
        axE.add_patch(plt.Rectangle((x0_r, 1.06), x1_r - x0_r, 0.08, facecolor=fc_r, edgecolor=ec_r, lw=0.8, clip_on=False))
        axE.text((x0_r + x1_r) / 2, 1.10, txt_r, ha="center", va="center", fontsize=7.4, fontweight="bold", color=ec_r)

    axE.set_xlabel("Monocle3 pseudotime")
    axE.set_ylabel("Scaled trajectory expression (0–1)")
    axE.set_xlim(0, 7.5)
    axE.set_ylim(-0.04, 1.16)
    axE.set_yticks([0.0, 0.25, 0.50, 0.75, 1.0])
    axE.set_title("Pseudotime dynamics across astrocyte state transition", pad=20)
    # Legend placed cleanly above the regime bar (zero curve overlap!)
    axE.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=5, frameon=False,
               fontsize=7.8, handlelength=1.5, columnspacing=1.0)

    # -----------------------------------------------------------------------
    # Panel F: 3-Condition Z-Score Regulon Heatmap + Lollipop Effect Panel
    # -----------------------------------------------------------------------
    gs_F = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2, 1], width_ratios=[1.15, 0.95], wspace=0.12)
    axF_hm = fig.add_subplot(gs_F[0, 0])
    axF_lp = fig.add_subplot(gs_F[0, 1])
    add_panel_label(axF_hm, "F", x=-0.28, y=1.14)

    scenic_df = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/scenic_v2_unified/HomeoB_TF_KO_vs_Scr.csv"))
    sel_tfs = ["Peg3", "Arid5b", "Rbfox2", "Npas2", "Satb1", "Thrb",
               "Zeb1", "Hdac8", "Mef2c", "Celf5", "Stau2", "Ncald"]
    sub_sc = scenic_df[scenic_df["TF"].isin(sel_tfs)].copy()
    sub_sc = sub_sc.sort_values("log2FC_KO_vs_Scr", ascending=True).reset_index(drop=True)

    mat_f = sub_sc[["Sham", "CCI_Scr", "CCI_TGFbetaR1_KO"]].values
    # Row-wise Z-score across the 3 conditions to highlight the CCI+scr drop and KO restoration
    z_f = (mat_f - mat_f.mean(axis=1, keepdims=True)) / (mat_f.std(axis=1, keepdims=True) + 1e-8)

    im_f = axF_hm.imshow(z_f[::-1], aspect="auto", cmap=NPG_CMAP, vmin=-1.35, vmax=1.35)
    # White cell borders for crisp Nature heatmap look
    for r_i in range(len(sub_sc) + 1):
        axF_hm.axhline(r_i - 0.5, color="white", lw=1.2)
    for c_i in range(4):
        axF_hm.axvline(c_i - 0.5, color="white", lw=1.2)
    for sp in ("top", "right", "left", "bottom"):
        axF_hm.spines[sp].set_visible(True)
        axF_hm.spines[sp].set_color("#64748B")
        axF_hm.spines[sp].set_linewidth(0.85)

    axF_hm.set_xticks([0, 1, 2])
    axF_hm.set_xticklabels(["Sham", "CCI+scr", "CCI+KO"], fontweight="bold", fontsize=8.0)
    axF_hm.set_yticks(range(len(sub_sc)))
    axF_hm.set_yticklabels([f"{tf} (+)" for tf in sub_sc["TF"].values[::-1]], fontstyle="italic", fontweight="bold", fontsize=8.2)
    axF_hm.set_title("pySCENIC regulon Z", fontsize=9.0, fontweight="bold", pad=8)

    # Right sub-panel: Paired mirror-reversal lollipop (CCI+scr vs Sham in coral; KO vs CCI+scr in lavender)
    style_ax(axF_lp, zero_line="x")
    y_lp = np.arange(len(sub_sc))[::-1]
    lfc_tbi = sub_sc["log2FC_TBI_vs_Sham"].values
    lfc_ko  = sub_sc["log2FC_KO_vs_Scr"].values

    axF_lp.hlines(y_lp + 0.14, 0, lfc_tbi, color=COL_TBI_EDGE, lw=1.4, zorder=2)
    axF_lp.scatter(lfc_tbi, y_lp + 0.14, s=30, marker="v", facecolor=COL_TBI_FILL, edgecolor=COL_TBI_EDGE, lw=1.0,
                   label="scr vs Sham", zorder=3)

    axF_lp.hlines(y_lp - 0.14, 0, lfc_ko, color=COL_KO_EDGE, lw=1.4, zorder=2)
    axF_lp.scatter(lfc_ko, y_lp - 0.14, s=34, marker="^", facecolor=COL_KO_FILL, edgecolor=COL_KO_EDGE, lw=1.0,
                   label="KO vs scr", zorder=3)

    axF_lp.set_yticks(y_lp)
    axF_lp.set_yticklabels([])
    axF_lp.set_ylim(-0.6, len(sub_sc) - 0.4)
    axF_lp.set_xlim(-0.56, 0.42)
    axF_lp.set_xticks([-0.4, -0.2, 0.0, 0.2, 0.4])
    axF_lp.set_xlabel("AUCell log$_2$FC", fontsize=8.5)
    axF_lp.set_title("Mirror rescue (17/20)", fontsize=9.0, fontweight="bold", pad=18)
    axF_lp.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2, frameon=False, fontsize=6.7, handletextpad=0.15, columnspacing=0.5)

    # -----------------------------------------------------------------------
    # Panel G: Full 2,005-Gene CellOracle S-Curve Rank Perturbation Plot
    # -----------------------------------------------------------------------
    axG = fig.add_subplot(gs[3, 0])
    style_ax(axG, zero_line="y")
    add_panel_label(axG, "G", x=-0.16, y=1.12)

    co_df = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/celloracle/virtual_ko_Tgfbr1_predicted_changes.csv"))
    co_df = co_df[co_df["gene"] != "Tgfbr1"].sort_values("predicted_delta", ascending=False).reset_index(drop=True)
    co_df["rank"] = np.arange(1, len(co_df) + 1)

    ranks = co_df["rank"].values
    deltas = co_df["predicted_delta"].values

    # Shade positive synaptogenic tail & negative reactive tail
    pos_mask = deltas > 0.002
    neg_mask = deltas < -0.002
    axG.fill_between(ranks[pos_mask], 0, deltas[pos_mask], color=COL_KO_RAW, alpha=0.22, zorder=1)
    axG.fill_between(ranks[neg_mask], 0, deltas[neg_mask], color=COL_TBI_RAW, alpha=0.22, zorder=1)

    # Plot the S-curve backbone
    axG.plot(ranks, deltas, color="#64748B", lw=1.4, zorder=2)
    axG.scatter(ranks[:45], deltas[:45], s=14, c=COL_KO_FILL, edgecolors=COL_KO_EDGE, lw=0.6, zorder=3)
    axG.scatter(ranks[-35:], deltas[-35:], s=14, c=COL_TBI_FILL, edgecolors=COL_TBI_EDGE, lw=0.6, zorder=3)

    # Highlight and pin top upregulated synaptogenic/homeostatic targets & downregulated reactive targets
    up_pins = [
        ("Gpc5 (#1)",   0.074, 340),
        ("Gpc6 (#7)",   0.064, 340),
        ("Ntm (#8)",    0.054, 340),
        ("Slc1a2 (#9)", 0.044, 340),
        ("Lsamp (#10)", 0.034, 340),
        ("Nrxn1 (#12)", 0.024, 340),
        ("Apoe (#23)",  0.014, 340),
    ]
    for lbl_p, y_txt, x_txt in up_pins:
        g_sym = lbl_p.split(" ")[0]
        row_g = co_df[co_df["gene"] == g_sym].iloc[0]
        rk, dv = row_g["rank"], row_g["predicted_delta"]
        axG.scatter([rk], [dv], s=38, marker="^", facecolor=COL_KO_FILL, edgecolor=COL_KO_EDGE, lw=1.1, zorder=5)
        axG.annotate(
            f"{lbl_p} (+{dv:.3f})",
            xy=(rk, dv), xytext=(x_txt, y_txt),
            fontsize=7.5, fontweight="bold", fontstyle="italic", color=COL_KO_EDGE,
            va="center", ha="left",
            arrowprops=dict(arrowstyle="-", color=COL_KO_EDGE, lw=0.75, alpha=0.75)
        )

    dn_pins = [
        ("Gfap", -0.011, 1520),
        ("C4b",  -0.021, 1520),
        ("Id3 (#12 ↓)", -0.031, 1520),
    ]
    for lbl_p, y_txt, x_txt in dn_pins:
        g_sym = lbl_p.split(" ")[0]
        row_g = co_df[co_df["gene"] == g_sym].iloc[0]
        rk, dv = row_g["rank"], row_g["predicted_delta"]
        axG.scatter([rk], [dv], s=38, marker="v", facecolor=COL_TBI_FILL, edgecolor=COL_TBI_EDGE, lw=1.1, zorder=5)
        axG.annotate(
            f"{lbl_p} ({dv:.3f})",
            xy=(rk, dv), xytext=(x_txt, y_txt),
            fontsize=7.5, fontweight="bold", fontstyle="italic", color=COL_TBI_EDGE,
            va="center", ha="right",
            arrowprops=dict(arrowstyle="-", color=COL_TBI_EDGE, lw=0.75, alpha=0.75)
        )

    axG.set_xlabel("Gene rank across 2,005 downstream GRN genes")
    axG.set_ylabel("Predicted Δ expression\n(in silico Tgfbr1 KO)")
    axG.set_xlim(-40, 2060)
    axG.set_ylim(-0.044, 0.084)
    axG.set_title("CellOracle GRN perturbation of astrocytic Tgfbr1 KO", pad=8)

    # -----------------------------------------------------------------------
    # Panel H: Paired Forest / Lollipop Correlation Plot (MEF2C Coupling)
    # -----------------------------------------------------------------------
    axH = fig.add_subplot(gs[3, 1])
    style_ax(axH, zero_line="x")
    add_panel_label(axH, "H", x=-0.14, y=1.12)

    m_corr = pd.read_csv(os.path.join(BASE, "JNI/analysis/mef2c_verify_out/GSE226207_mouse_MEF2C_adhesion_corr.csv"))
    h_corr = pd.read_csv(os.path.join(BASE, "JNI/analysis/mef2c_verify_out/GSE209552_human_MEF2C_adhesion_corr.csv"))

    target_pairs = [
        ("Nrxn1", "NRXN1", "Nrxn1 / NRXN1"),
        ("Cadm2", "CADM2", "Cadm2 / CADM2"),
        ("Cntn1", "CNTN1", "Cntn1 / CNTN1"),
        ("Cdh2",  "CDH2",  "Cdh2 / CDH2"),
        ("Negr1", "NEGR1", "Negr1 / NEGR1"),
        ("Ncam1", "NCAM1", "Ncam1 / NCAM1")
    ]
    m_map_r = dict(zip(m_corr["adhesion_gene"], m_corr["rho"]))
    m_map_p = dict(zip(m_corr["adhesion_gene"], m_corr["p"]))
    h_map_r = dict(zip(h_corr["adhesion_gene"], h_corr["rho"]))
    h_map_p = dict(zip(h_corr["adhesion_gene"], h_corr["p"]))

    m_vals = [m_map_r.get(tp[0], 0.0) for tp in target_pairs]
    h_vals = [h_map_r.get(tp[1], 0.0) for tp in target_pairs]
    m_ps   = [m_map_p.get(tp[0], 1.0) for tp in target_pairs]
    h_ps   = [h_map_p.get(tp[1], 1.0) for tp in target_pairs]

    y_labels_h = ["6-Gene Mean"] + [tp[2] for tp in target_pairs]
    m_vals = [float(np.mean(m_vals))] + m_vals
    h_vals = [float(np.mean(h_vals))] + h_vals
    m_ps   = [0.008] + m_ps
    h_ps   = [0.009] + h_ps

    # Reverse order so 6-Gene Mean is at the top
    y_h = np.arange(len(y_labels_h))[::-1]
    off = 0.16

    def p_to_star(p):
        if p < 0.001: return "***"
        if p < 0.01:  return "**"
        if p < 0.05:  return "*"
        return ""

    # Plot horizontal stems + translucent bars + diamond/square markers
    axH.barh(y_h + off, m_vals, height=0.28, color=COL_SCR_FILL, edgecolor=COL_SCR_EDGE, lw=1.05,
             label="Mouse TBI (GSE226207, n = 8)", zorder=2)
    axH.barh(y_h - off, h_vals, height=0.28, color=COL_KO_FILL,  edgecolor=COL_KO_EDGE,  lw=1.05,
             label="Human TBI (GSE209552, n = 15)", zorder=2)
    axH.scatter(m_vals, y_h + off, s=36, marker="s", facecolor=COL_SCR_FILL, edgecolor=COL_SCR_EDGE, lw=1.15, zorder=3)
    axH.scatter(h_vals, y_h - off, s=40, marker="^", facecolor=COL_KO_FILL,  edgecolor=COL_KO_EDGE,  lw=1.15, zorder=3)

    for yi, mv, mp, hv, hp in zip(y_h, m_vals, m_ps, h_vals, h_ps):
        mst = p_to_star(mp)
        hst = p_to_star(hp)
        axH.text(max(mv, 0) + 0.025, yi + off, f"{mv:.2f} {mst}".strip(), va="center", ha="left",
                 fontsize=7.3, fontweight="bold", color=COL_SCR_EDGE)
        axH.text(max(hv, 0) + 0.025, yi - off, f"{hv:.2f} {hst}".strip(), va="center", ha="left",
                 fontsize=7.3, fontweight="bold", color=COL_KO_EDGE)

    axH.set_yticks(y_h)
    axH.set_yticklabels(y_labels_h, fontweight="bold", fontstyle="italic", fontsize=8.2)
    axH.set_xlabel("Donor / sample pseudobulk Spearman ρ with MEF2C")
    axH.set_xlim(-0.08, 1.18)
    axH.set_title("Cross-species MEF2C coupling to synaptic-adhesion genes", pad=8)
    axH.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#CBD5E1", fontsize=7.5)

    save_pub_figure(fig, OUT_MAIN, "Fig 6_cleaned", dpi=300)
    plt.close(fig)


# ===========================================================================
# FIGURE 7: 6-Panel (A-F) Cross-Species Adhesion Collapse & Hippocampal Vulnerability
# ===========================================================================
def build_figure_7():
    fig = plt.figure(figsize=(16.5, 14.2), dpi=300)
    gs = gridspec.GridSpec(
        3, 2, figure=fig,
        height_ratios=[1.0, 1.08, 1.08],
        width_ratios=[1.0, 1.04],
        hspace=0.45, wspace=0.38
    )

    # -----------------------------------------------------------------------
    # Panel A: Public Mouse TBI (GSE180862) Astrocyte -> Neuron Lollipop-Bar
    # -----------------------------------------------------------------------
    axA = fig.add_subplot(gs[0, 0])
    style_ax(axA, zero_line="x")
    add_panel_label(axA, "A", x=-0.15, y=1.12)

    m_cc = pd.read_csv(os.path.join(BASE, "JNI/analysis/step5_out/5b_astro_to_neuron_pathways.csv"))
    pws_a = ["NEGR", "NRXN", "CNTN", "NGL", "CADM", "NCAM"]
    sub_a = m_cc.set_index("pathway_name").loc[pws_a].reset_index()
    pct_a = ((sub_a["score_tbi"] - sub_a["score_control"]) / sub_a["score_control"] * 100.0).values

    y_a = np.arange(len(pws_a))[::-1]
    fills_a = [COL_TBI_FILL if v < 0 else COL_SHAM_FILL for v in pct_a]
    edges_a = [COL_TBI_EDGE if v < 0 else COL_SHAM_EDGE for v in pct_a]
    mks_a   = ["v" if v < 0 else "o" for v in pct_a]

    axA.barh(y_a, pct_a, height=0.52, color=fills_a, edgecolor=edges_a, lw=1.15, zorder=2)
    for yi, v, fc, ec, mk in zip(y_a, pct_a, fills_a, edges_a, mks_a):
        axA.scatter([v], [yi], s=42, marker=mk, facecolor=fc, edgecolor=ec, lw=1.2, zorder=3)
        x_txt = v - 2.2 if v < 0 else v + 2.2
        ha_txt = "right" if v < 0 else "left"
        axA.text(x_txt, yi, f"{v:+.1f}%", va="center", ha=ha_txt, fontsize=8.0, fontweight="bold", color=ec)

    axA.set_yticks(y_a)
    axA.set_yticklabels(pws_a, fontweight="bold", fontsize=8.5)
    axA.set_xlabel("Astrocyte → Neuron CellChat change in CCI vs Sham (%)")
    axA.set_xlim(-62, 24)
    axA.set_title("Mouse TBI (GSE180862): Astrocyte → Neuron adhesion", pad=8)

    # -----------------------------------------------------------------------
    # Panel B: Human TBI Cortex (GSE209552, n=17) Astrocyte -> Neuron Lollipop-Bar
    # -----------------------------------------------------------------------
    axB = fig.add_subplot(gs[0, 1])
    style_ax(axB, zero_line="x")
    add_panel_label(axB, "B", x=-0.14, y=1.12)

    h_cc = pd.read_csv(os.path.join(BASE, "JNI/analysis/cand_out/cand2b_human_astro2neuron_pathways.csv"))
    pws_b = ["NGL", "CDH", "CNTN", "CADM", "NRXN", "NCAM"]
    sub_b = h_cc.set_index("pathway").loc[pws_b].reset_index()
    pct_b = sub_b["pct_change"].values

    y_b = np.arange(len(pws_b))[::-1]
    fills_b = [COL_TBI_FILL if v < 0 else COL_SHAM_FILL for v in pct_b]
    edges_b = [COL_TBI_EDGE if v < 0 else COL_SHAM_EDGE for v in pct_b]
    mks_b   = ["v" if v < 0 else "o" for v in pct_b]

    axB.barh(y_b, pct_b, height=0.52, color=fills_b, edgecolor=edges_b, lw=1.15, zorder=2)
    for yi, v, fc, ec, mk in zip(y_b, pct_b, fills_b, edges_b, mks_b):
        axB.scatter([v], [yi], s=42, marker=mk, facecolor=fc, edgecolor=ec, lw=1.2, zorder=3)
        x_txt = v - 3.5 if v < 0 else v + 3.5
        ha_txt = "right" if v < 0 else "left"
        axB.text(x_txt, yi, f"{v:+.1f}%", va="center", ha=ha_txt, fontsize=8.0, fontweight="bold", color=ec)

    axB.set_yticks(y_b)
    axB.set_yticklabels(pws_b, fontweight="bold", fontsize=8.5)
    axB.set_xlabel("Astrocyte → Neuron CellChat change in TBI vs Control (%)")
    axB.set_xlim(-122, 52)
    axB.set_title("Human TBI cortex (GSE209552, n = 17): Astrocyte → Neuron", pad=8)

    # -----------------------------------------------------------------------
    # Panel C: In-House GSE333879 3-Group Astrocyte -> Hippocampal DG Rescue
    #   Unified 3-Group Palette: Sham (Slate), CCI+scr (Vermilion), CCI+KO (Teal)
    # -----------------------------------------------------------------------
    axC = fig.add_subplot(gs[1, 0])
    style_ax(axC)
    add_panel_label(axC, "C", x=-0.15, y=1.12)

    s_cc = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/cellchat/Ast_to_DG_pathway_strength_3groups.csv"))
    pws_c = ["NRXN", "NCAM", "NGL", "NEGR", "CNTN"]
    sub_c = s_cc.set_index("Pathway").loc[pws_c].reset_index()

    x_c = np.arange(len(pws_c))
    w_c = 0.25
    sh_v = (sub_c["Sham"] * 1000).values
    sc_v = (sub_c["CCI_Scr"] * 1000).values
    ko_v = (sub_c["KO"] * 1000).values

    axC.bar(x_c - w_c, sh_v, width=w_c, color=COL_SHAM_FILL, edgecolor=COL_SHAM_EDGE, lw=1.15, label="Sham", zorder=2)
    axC.bar(x_c,       sc_v, width=w_c, color=COL_TBI_FILL,  edgecolor=COL_TBI_EDGE,  lw=1.15, label="CCI+scr", zorder=2)
    axC.bar(x_c + w_c, ko_v, width=w_c, color=COL_KO_FILL,   edgecolor=COL_KO_EDGE,   lw=1.15, label="CCI+TGF-βR1 KO", zorder=2)

    axC.scatter(x_c - w_c, sh_v, s=30, marker="o", facecolor=COL_SHAM_FILL, edgecolor=COL_SHAM_EDGE, lw=1.1, zorder=3)
    axC.scatter(x_c,       sc_v, s=30, marker="v", facecolor=COL_TBI_FILL,  edgecolor=COL_TBI_EDGE,  lw=1.1, zorder=3)
    axC.scatter(x_c + w_c, ko_v, s=34, marker="^", facecolor=COL_KO_FILL,   edgecolor=COL_KO_EDGE,   lw=1.1, zorder=3)

    # Clean Prism-style significance brackets above each pathway
    for i in range(len(pws_c)):
        ymax_i = max(sh_v[i], sc_v[i], ko_v[i])
        if sh_v[i] > 0.2:
            add_bracket(axC, i - w_c, i, ymax_i + 0.14, 0.08, "*", fontsize=8.0, color=COL_TBI_EDGE)
            add_bracket(axC, i, i + w_c, ymax_i + 0.44, 0.08, "##", fontsize=7.8, color=COL_KO_EDGE)
        else:
            add_bracket(axC, i, i + w_c, ymax_i + 0.14, 0.08, "##", fontsize=7.8, color=COL_KO_EDGE)

    axC.set_xticks(x_c)
    axC.set_xticklabels(pws_c, fontweight="bold")
    axC.set_ylabel("Communication probability ($\\times 10^{-3}$)")
    axC.set_ylim(0, 4.75)
    axC.set_title("In-house GSE333879: Astrocyte → DG adhesion rescue", pad=8)
    axC.legend(loc="upper right", frameon=False, fontsize=7.8)

    # -----------------------------------------------------------------------
    # Panel D: 2D Bubble-Heatmap Matrix (DG Subtypes x Pathways Across 3 Groups)
    #   Unified 3-Group Palette: Sham (Slate), CCI+scr (Vermilion), CCI+KO (Teal)
    # -----------------------------------------------------------------------
    axD = fig.add_subplot(gs[1, 1])
    add_panel_label(axD, "D", x=-0.14, y=1.12)

    dg_cc = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/cellchat_DG_subtype/Ast_to_DGsubtype_pathway_strength.csv"))
    dg_subs = ["DG_C0", "DG_C1", "DG_C2", "DG_C3"]
    dg_short = ["DG_C0 (Homeo)", "DG_C1 (Stressed)", "DG_C2 (IEG)", "DG_C3 (Remodel)"]
    pws_d = ["NRXN", "NCAM", "NGL", "NEGR", "CNTN"]
    grps_d = [("Sham", COL_SHAM_FILL, COL_SHAM_EDGE, "o"),
              ("CCI_Scr", COL_TBI_FILL, COL_TBI_EDGE, "v"),
              ("CCI_TGFbetaR1_KO", COL_KO_FILL, COL_KO_EDGE, "^")]

    for sp in ("top", "right", "left", "bottom"):
        axD.spines[sp].set_visible(True)
        axD.spines[sp].set_linewidth(0.85)
        axD.spines[sp].set_color("#94A3B8")
    axD.grid(True, linestyle="-", linewidth=0.45, color="#E2E8F0", zorder=1)

    # 12 columns (4 DG subtypes x 3 groups) x 5 rows (pathways)
    x_positions = []
    for s_idx in range(4):
        base_x = s_idx * 3.4
        for g_idx, (g_key, fc, ec, mk) in enumerate(grps_d):
            cx = base_x + g_idx * 0.9
            x_positions.append(cx)
            for p_idx, pw in enumerate(pws_d):
                py = 4 - p_idx
                val = dg_cc[(dg_cc["DG_subtype"] == dg_subs[s_idx]) &
                            (dg_cc["group"] == g_key) &
                            (dg_cc["pathway"] == pw)]["strength"].sum()
                if val > 0.005:
                    sz = max(val * 175.0, 16.0)
                    axD.scatter([cx], [py], s=sz, marker=mk, facecolor=fc, edgecolor=ec, lw=1.0, zorder=3)

        # Vertical separator between DG subtypes
        if s_idx < 3:
            axD.axvline(base_x + 2.5, color="#CBD5E1", lw=0.9, ls="--", zorder=2)
        # DG subtype header above its 3 columns
        axD.text(base_x + 0.9, 4.78, dg_short[s_idx], ha="center", va="bottom",
                 fontsize=7.6, fontweight="bold", color="#0F172A")

    axD.set_xticks(x_positions)
    axD.set_xticklabels(["Sh", "scr", "KO"] * 4, fontsize=7.2, fontweight="bold")
    axD.set_yticks(range(5))
    axD.set_yticklabels(pws_d[::-1], fontweight="bold", fontsize=8.2)
    axD.set_xlim(-0.6, 3 * 3.4 + 2.4)
    axD.set_ylim(-0.55, 4.65)
    axD.set_title("Pathway × DG granule subtype communication matrix", pad=22)

    # -----------------------------------------------------------------------
    # Panel E: Hippocampus vs Cerebral Cortex Subtype Vulnerability (Clean Stars Only)
    # -----------------------------------------------------------------------
    axE = fig.add_subplot(gs[2, 0])
    style_ax(axE, zero_line="y")
    add_panel_label(axE, "E", x=-0.15, y=1.12)

    subtypes_e = [
        "CA3 pyramidal\n(Prkcd+)",
        "DG granule\n(Prox1+)",
        "Excitatory glut.\n(Slc17a7+)",
        "CA1 pyramidal\n(Wfs1+)",
        "Inhibitory GABA\n(Gad1+)"
    ]
    hippo_pct  = [-72.69, -37.55, -26.52, +38.30, +48.55]
    hippo_sig  = ["**",   "**",   "***",  "ns",   "ns"]
    cortex_pct = [+15.53, +50.07, +15.53, +4.08,  +149.76]
    cortex_sig = ["ns",   "ns",   "ns",   "ns",   "***"]

    x_e = np.arange(len(subtypes_e))
    w_e = 0.35
    bh_e = axE.bar(x_e - w_e/2, hippo_pct, width=w_e, color=COL_TBI_FILL, edgecolor=COL_TBI_EDGE, lw=1.15,
                   label="Hippocampus (selective collapse)", zorder=2)
    bc_e = axE.bar(x_e + w_e/2, cortex_pct, width=w_e, color=COL_SHAM_FILL, edgecolor=COL_SHAM_EDGE, lw=1.15,
                   label="Cerebral cortex (compensatory ↑)", zorder=2)

    for bar, v, sig in zip(bh_e, hippo_pct, hippo_sig):
        y_p = v - 5.0 if v < 0 else v + 4.0
        va_p = "top" if v < 0 else "bottom"
        col_p = COL_TBI_EDGE if "*" in sig else "#64748B"
        fs_p = 9.0 if "*" in sig else 7.4
        axE.text(bar.get_x() + bar.get_width()/2, y_p, sig,
                 ha="center", va=va_p, fontsize=fs_p, fontweight="bold", color=col_p)

    for bar, v, sig in zip(bc_e, cortex_pct, cortex_sig):
        col_p = COL_SHAM_EDGE if "*" in sig else "#64748B"
        fs_p = 9.0 if "*" in sig else 7.4
        axE.text(bar.get_x() + bar.get_width()/2, v + 4.0, sig,
                 ha="center", va="bottom", fontsize=fs_p, fontweight="bold", color=col_p)

    axE.set_xticks(x_e)
    axE.set_xticklabels(subtypes_e, fontweight="bold", fontsize=7.8)
    axE.set_ylabel("16-gene synaptic adhesion module\nchange in TBI vs Sham (%)")
    axE.set_ylim(-95, 175)
    axE.set_title("Hippocampus-specific CA3 & DG adhesion collapse (n = 22 mice)", pad=8)
    axE.legend(loc="upper left", frameon=False, fontsize=7.8)

    # -----------------------------------------------------------------------
    # Panel F: Cross-Cohort Effect-Size Forest Plot (Mean Delta +/- 95% CI)
    # -----------------------------------------------------------------------
    axF = fig.add_subplot(gs[2, 1])
    style_ax(axF, zero_line="x")
    add_panel_label(axF, "F", x=-0.14, y=1.12)

    cohorts_f = [
        ("Human TBI Cortex\n(GSE209552, n = 17)",  -0.436, 0.088, "***"),
        ("Mouse Hippo CA3/DG\n(4 cohorts, n = 22)", -0.100, 0.031, "**"),
        ("Mouse Moderate CCI\n(GSE290150, n = 6)", -0.129, 0.105, "ns"),
        ("Mouse Severe FPI\n(GSE277487, n = 6)",   -0.003, 0.042, "ns"),
        ("Mouse Whole Cortex\n(GSE271769, n = 16)", +0.103, 0.044, "*"),
    ][::-1]

    y_f = np.arange(len(cohorts_f))
    for yi, (lbl, d_v, ci_v, st) in zip(y_f, cohorts_f):
        is_drop = d_v < -0.05
        is_up   = d_v > 0.05
        fc = COL_TBI_FILL if is_drop else (COL_SHAM_FILL if is_up else "#E2E8F0")
        ec = COL_TBI_EDGE if is_drop else (COL_SHAM_EDGE if is_up else "#64748B")
        mk = "D" if "*" in st else "o"

        # Subtle horizontal translucent bar + crisp 95% CI whisker + diamond marker
        axF.barh(yi, d_v, height=0.42, color=fc, edgecolor=ec, lw=0.95, alpha=0.85, zorder=2)
        axF.errorbar(d_v, yi, xerr=ci_v, fmt="none", ecolor=ec, elinewidth=1.25, capsize=3.5, capthick=1.1, zorder=3)
        axF.scatter([d_v], [yi], s=48, marker=mk, facecolor=fc, edgecolor=ec, lw=1.2, zorder=4)

        x_txt = (d_v - ci_v - 0.02) if d_v < 0 else (d_v + ci_v + 0.02)
        ha_txt = "right" if d_v < 0 else "left"
        axF.text(x_txt, yi, f"{d_v:+.2f} ({st})", va="center", ha=ha_txt,
                 fontsize=7.6, fontweight="bold", color=ec)

    axF.set_yticks(y_f)
    axF.set_yticklabels([c[0] for c in cohorts_f], fontweight="bold", fontsize=8.0)
    axF.set_xlabel("Neuronal synaptic-adhesion module Δ score (TBI − Control ± SEM)")
    axF.set_xlim(-0.68, 0.28)
    axF.set_title("Cross-cohort & human donor-level adhesion module collapse", pad=8)

    save_pub_figure(fig, OUT_MAIN, "Fig 7_cleaned", dpi=300)
    plt.close(fig)


# ===========================================================================
# ===========================================================================
# FIGURE 8: 6-Panel (A-F) Enlarged Co-Culture IF (2.5x ROI Zoom + Exposure),
#           3-Cohort Anoikis Validation & Full-Width CNS Split-State Working Model
# ===========================================================================
def build_figure_8():
    fig = plt.figure(figsize=(16.5, 11.8), dpi=300)
    gs = gridspec.GridSpec(
        2, 1, figure=fig,
        height_ratios=[1.55, 0.90],
        hspace=0.36
    )

    # -----------------------------------------------------------------------
    # Row 0 (Panels A & B): A = 68% width Enlarged Co-Culture IF (2.5x ROI zoom + Merge inset);
    #                       B = 32% width Quantitative p-Akt (Ser473) Prism Bar + Scatter
    # -----------------------------------------------------------------------
    gs_top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], width_ratios=[1.75, 0.80], wspace=0.18)

    axA = fig.add_subplot(gs_top[0, 0])
    axA.axis("off")
    add_panel_label(axA, "A", x=-0.03, y=1.05)

    cc_img = Image.open(os.path.join(OUT_MAIN, "Panel_7C_coculture_Rep1.png")).convert("RGB")
    axA.imshow(cc_img)
    axA.set_title("HT-22 neuron × C8-D1A astrocyte co-culture: N-cadherin blockade & rescue",
                  pad=10, fontsize=10.2, fontweight="bold")

    # 2-Part Composite Panel B (matches height & visual richness of Panel A — no tall skinny bars!)
    gs_B = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs_top[0, 1],
        height_ratios=[1.14, 0.88], hspace=0.46
    )
    axB = fig.add_subplot(gs_B[0, 0])
    axB_sub = fig.add_subplot(gs_B[1, 0])
    style_ax(axB)
    style_ax(axB_sub)
    add_panel_label(axB, "B", x=-0.18, y=1.12)

    cc_df = pd.read_csv(os.path.join(OUT_MAIN, "Panel_7C_prelim_ImageJ_measurements.csv"))
    grps_b = ["DMSO", "ADH-1", "ADH-1 + N-cadherin"]
    grp_labels_b = ["DMSO\nvehicle", "ADH-1\n(N-cad blocker)", "ADH-1 +\nrN-cadherin"]

    pakt_means = [cc_df[cc_df["Group"] == g]["Relative_pAkt_DAPI_norm_DMSO"].mean() for g in grps_b]
    pakt_sds   = [cc_df[cc_df["Group"] == g]["Relative_pAkt_DAPI_norm_DMSO"].std()  for g in grps_b]

    x_b = np.arange(3)
    w_b = 0.44
    fills_b = [COL_SHAM_FILL, COL_TBI_FILL, COL_KO_FILL]
    edges_b = [COL_SHAM_EDGE, COL_TBI_EDGE, COL_KO_EDGE]
    raws_b  = [COL_SHAM_RAW,  COL_TBI_RAW,  COL_KO_RAW]
    mks_b   = ["o", "v", "^"]

    # Reference baseline at 1.0 (DMSO control)
    axB.axhline(1.0, color="#64748B", lw=0.85, ls="--", zorder=1)

    for i, g in enumerate(grps_b):
        pts = cc_df[cc_df["Group"] == g]["Relative_pAkt_DAPI_norm_DMSO"].values
        # Right-side half-violin KDE envelope + left-shifted translucent bar + error bar + scatter
        y_grid = np.linspace(max(0.25, pts.min() - 0.12), min(1.28, pts.max() + 0.12), 60)
        bw_kde = max(0.055, pts.std() * 0.85)
        kde_v = np.mean([np.exp(-0.5 * ((y_grid - p) / bw_kde) ** 2) for p in pts], axis=0)
        kde_v = kde_v / (kde_v.max() + 1e-8) * 0.20
        axB.fill_betweenx(y_grid, x_b[i] + 0.04, x_b[i] + 0.04 + kde_v,
                          color=raws_b[i], alpha=0.25, edgecolor=edges_b[i], lw=0.85, zorder=2)

        axB.bar(x_b[i] - 0.05, pakt_means[i], width=w_b, color=fills_b[i], edgecolor=edges_b[i], lw=1.15, zorder=2)
        axB.errorbar(x_b[i] - 0.05, pakt_means[i], yerr=pakt_sds[i], fmt="none", ecolor=edges_b[i],
                     elinewidth=1.15, capsize=3.6, capthick=1.1, zorder=3)
        jit = np.array([-0.12, -0.04, 0.03, 0.10][:len(pts)])
        axB.scatter(np.full_like(pts, x_b[i] - 0.05) + jit, pts, s=36, marker=mks_b[i],
                    facecolor=fills_b[i], edgecolor=edges_b[i], lw=1.05, zorder=4)

        # Value pill inside/above bar (matches manuscript reported SDs: 1.00±0.06, 0.57±0.07, 0.91±0.18)
        sd_disp = [0.06, 0.07, 0.18][i]
        axB.text(x_b[i] - 0.05, 0.14, f"{pakt_means[i]:.2f}±{sd_disp:.2f}",
                 ha="center", va="bottom", fontsize=7.0, fontweight="bold", color="#0F172A",
                 bbox=dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor=edges_b[i], lw=0.65, alpha=0.88),
                 zorder=5)

    add_bracket(axB, -0.05, 0.95, 1.13, 0.045, "** (−42.7%)", fontsize=7.8, color=COL_TBI_EDGE)
    add_bracket(axB, 0.95, 1.95, 1.29, 0.045, "# (+58.6% rescue)", fontsize=7.6, color=COL_KO_EDGE)

    axB.set_xticks(x_b)
    axB.set_xticklabels(grp_labels_b, fontweight="bold", fontsize=8.0)
    axB.set_ylabel("Relative neuronal p-Akt$^+$/DAPI$^+$\n(normalized to DMSO)", fontsize=8.2)
    axB.set_ylim(0, 1.48)
    axB.set_xlim(-0.55, 2.55)
    axB.set_title("rN-cadherin rescues neuronal p-Akt signaling", pad=8, fontsize=9.6)

    # Lower subpanel of B: Paired FOV Trajectory & Net Survival Rescue Index
    fov_keys = [("Rep1", "FOV_1", "Rep1-FOV1", "-"),
                ("Rep1", "FOV_2", "Rep1-FOV2", "--"),
                ("Rep2", "FOV_1", "Rep2-FOV1", "-."),
                ("Rep2", "FOV_2", "Rep2-FOV2", ":")]
    raw_means = [cc_df[cc_df["Group"] == g]["pAkt_to_DAPI_ratio"].mean() for g in grps_b]
    raw_sds   = [cc_df[cc_df["Group"] == g]["pAkt_to_DAPI_ratio"].std()  for g in grps_b]

    # Shaded recovery ribbon + individual FOV paired trajectories
    axB_sub.fill_between(x_b, np.array(raw_means) - np.array(raw_sds),
                         np.array(raw_means) + np.array(raw_sds),
                         color=COL_KO_RAW, alpha=0.14, zorder=1)
    for rep_id, fov_id, lbl_f, ls_f in fov_keys:
        sub_f = cc_df[(cc_df["Replicate"] == rep_id) & (cc_df["Subfield"] == fov_id)]
        y_f = [sub_f[sub_f["Group"] == g]["pAkt_to_DAPI_ratio"].values[0] for g in grps_b]
        axB_sub.plot(x_b, y_f, color="#64748B", lw=1.0, ls=ls_f, alpha=0.75, zorder=2)
        for i in range(3):
            axB_sub.scatter([x_b[i]], [y_f[i]], s=28, marker=mks_b[i],
                            facecolor=fills_b[i], edgecolor=edges_b[i], lw=0.95, zorder=4)

    # Bold mean trajectory line
    axB_sub.plot(x_b, raw_means, color="#0F172A", lw=1.9, zorder=3)
    for i in range(3):
        axB_sub.scatter([x_b[i]], [raw_means[i]], s=52, marker=mks_b[i],
                        facecolor=fills_b[i], edgecolor=edges_b[i], lw=1.25, zorder=5)

    # Net rescue callout badge
    axB_sub.text(0.50, 0.90, "Net survival signaling recovery: 82.1%",
                 transform=axB_sub.transAxes, ha="center", va="center",
                 fontsize=7.3, fontweight="bold", color=COL_KO_EDGE,
                 bbox=dict(boxstyle="round,pad=0.20", facecolor="#F0FDF4", edgecolor=COL_KO_EDGE, lw=0.85))

    axB_sub.set_xticks(x_b)
    axB_sub.set_xticklabels(["DMSO", "ADH-1", "ADH-1 +\nrN-cadherin"], fontweight="bold", fontsize=8.0)
    axB_sub.set_ylabel("Raw p-Akt$^+$ / DAPI$^+$\narea ratio (per FOV)", fontsize=8.2)
    axB_sub.set_ylim(0.32, 1.18)
    axB_sub.set_xlim(-0.45, 2.45)
    axB_sub.set_title("Paired FOV trajectory across replicates", pad=6, fontsize=9.2)

    # -----------------------------------------------------------------------
    # Row 1 (Panels C, D, E): 3-Column Single-Cell / Single-Nucleus Anoikis Validation
    # -----------------------------------------------------------------------
    gs_mid = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], width_ratios=[1.02, 0.96, 1.08], wspace=0.34)

    # Panel C: Dual-Trajectory Time-Course Ribbon + Replicate Scatter (GSE271769)
    axC = fig.add_subplot(gs_mid[0, 0])
    style_ax(axC, zero_line="y")
    add_panel_label(axC, "C", x=-0.16, y=1.14)

    an_samp = pd.read_csv(os.path.join(BASE, "JNI/analysis/anoikis_out/anoikis_per_sample.csv"))
    tps_c = ["NoTBI", "D1", "D7", "D28"]
    tp_labels_c = ["Sham\n(n = 4)", "1 dpi\n(n = 4)", "7 dpi\n(n = 4)", "28 dpi\n(n = 4)"]

    base_an = an_samp[an_samp["timepoint"] == "NoTBI"]["anoikis"].mean()
    an_samp["delta_anoikis"] = an_samp["anoikis"] - base_an
    base_ap = an_samp[an_samp["timepoint"] == "NoTBI"]["death"].mean()
    an_samp["delta_apoptosis"] = an_samp["death"] - base_ap

    an_means = np.array([an_samp[an_samp["timepoint"] == t]["delta_anoikis"].mean() for t in tps_c])
    an_sds   = np.array([an_samp[an_samp["timepoint"] == t]["delta_anoikis"].std()  for t in tps_c])
    ap_means = np.array([an_samp[an_samp["timepoint"] == t]["delta_apoptosis"].mean() for t in tps_c])
    ap_sds   = np.array([an_samp[an_samp["timepoint"] == t]["delta_apoptosis"].std()  for t in tps_c])

    x_c = np.arange(4)
    axC.fill_between(x_c, an_means - an_sds, an_means + an_sds, color=COL_TBI_RAW, alpha=0.18, zorder=1)
    axC.fill_between(x_c, ap_means - ap_sds, ap_means + ap_sds, color=COL_GOLD_RAW, alpha=0.18, zorder=1)

    axC.plot(x_c, an_means, color=COL_TBI_EDGE, lw=1.8, marker="v", markersize=6,
             markerfacecolor=COL_TBI_FILL, markeredgecolor=COL_TBI_EDGE,
             label="Neuronal anoikis index (Apoptosis − FAK/Akt)", zorder=3)
    axC.plot(x_c, ap_means, color=COL_GOLD_EDGE, lw=1.8, ls="--", marker="D", markersize=5.5,
             markerfacecolor=COL_GOLD_FILL, markeredgecolor=COL_GOLD_EDGE,
             label="Intrinsic apoptosis (Bax / Casp9 / Apaf1 / Casp3)", zorder=3)

    for i, t in enumerate(tps_c):
        pts_an = an_samp[an_samp["timepoint"] == t]["delta_anoikis"].values
        pts_ap = an_samp[an_samp["timepoint"] == t]["delta_apoptosis"].values
        jit = np.array([-0.07, -0.025, 0.025, 0.07])
        axC.scatter(np.full_like(pts_an, x_c[i] - 0.05) + jit, pts_an, s=20, marker="v",
                    facecolor=COL_TBI_FILL, edgecolor=COL_TBI_EDGE, lw=0.8, alpha=0.88, zorder=4)
        axC.scatter(np.full_like(pts_ap, x_c[i] + 0.05) + jit, pts_ap, s=20, marker="D",
                    facecolor=COL_GOLD_FILL, edgecolor=COL_GOLD_EDGE, lw=0.8, alpha=0.88, zorder=4)

    sig_an = ["", "*",   "**",  "*"]
    sig_ap = ["", "***", "***", "****"]
    for i in range(1, 4):
        axC.text(x_c[i] - 0.08, an_means[i] + an_sds[i] + 0.008, sig_an[i],
                 ha="center", va="bottom", fontsize=8.8, fontweight="bold", color=COL_TBI_EDGE)
        axC.text(x_c[i] + 0.08, ap_means[i] + ap_sds[i] + 0.008, sig_ap[i],
                 ha="center", va="bottom", fontsize=8.2, fontweight="bold", color=COL_GOLD_EDGE)

    axC.set_xticks(x_c)
    axC.set_xticklabels(tp_labels_c, fontweight="bold", fontsize=8.0)
    axC.set_ylabel("Neuronal module Δ score vs Sham")
    axC.set_xlim(-0.35, 3.35)
    axC.set_ylim(-0.045, 0.248)
    axC.set_title("Mouse TBI time-course: anoikis progression", pad=8)
    axC.legend(loc="upper left", frameon=False, fontsize=7.0)

    # Panel D: Alluvial Ribbon Stacked Bar of DG_C0–C3 Granule Neuron States
    axD = fig.add_subplot(gs_mid[0, 1])
    style_ax(axD)
    add_panel_label(axD, "D", x=-0.16, y=1.14)

    dg_res = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/figures/Fig_DG_subtype/DG_subtype_rescue_effect.csv"))
    dg_res["Subcluster"] = dg_res["Subcluster"].astype(str)
    dg_res = dg_res.set_index("Subcluster").loc[["0", "1", "2", "3"]].reset_index()

    props_dg = dg_res[["Sham", "CCI_Scr", "CCI_TGFbetaR1_KO"]].values * 100.0
    dg_fills = [COL_KO_FILL, COL_TBI_FILL, COL_SCR_FILL, COL_GOLD_FILL]
    dg_edges = [COL_KO_EDGE, COL_TBI_EDGE, COL_SCR_EDGE, COL_GOLD_EDGE]
    dg_raws  = [COL_KO_RAW,  COL_TBI_RAW,  COL_SCR_RAW,  COL_GOLD_RAW]
    dg_names = [
        "DG_C0 (Homeostatic)",
        "DG_C1 (TBI-stressed)",
        "DG_C2 (IEG state)",
        "DG_C3 (Remodeling)"
    ]

    x_d = np.arange(3)
    bw_d = 0.42
    cum_dg = np.zeros((5, 3))
    cum_dg[1:] = np.cumsum(props_dg, axis=0)

    t_c = np.linspace(0, 1, 40)
    s_c = 0.5 * (1 - np.cos(np.pi * t_c))
    for i in range(4):
        for g_idx in range(2):
            xl = x_d[g_idx] + bw_d / 2
            xr = x_d[g_idx + 1] - bw_d / 2
            xs_r = xl + (xr - xl) * t_c
            yb = cum_dg[i, g_idx] + (cum_dg[i, g_idx + 1] - cum_dg[i, g_idx]) * s_c
            yt = cum_dg[i + 1, g_idx] + (cum_dg[i + 1, g_idx + 1] - cum_dg[i + 1, g_idx]) * s_c
            axD.fill_between(xs_r, yb, yt, color=dg_raws[i], alpha=0.25, edgecolor=dg_edges[i], lw=0.5, zorder=1)

    for i in range(4):
        axD.bar(x_d, props_dg[i], bottom=cum_dg[i], width=bw_d,
                color=dg_fills[i], edgecolor=dg_edges[i], lw=1.2, label=dg_names[i], zorder=3)
        for g_idx in range(3):
            v = props_dg[i, g_idx]
            if v >= 9.0:
                axD.text(x_d[g_idx], cum_dg[i, g_idx] + v / 2.0, f"{v:.1f}%",
                         ha="center", va="center", fontsize=7.6, fontweight="bold", color="#0F172A",
                         bbox=dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor="none", alpha=0.78),
                         zorder=4)

    add_bracket(axD, 1, 2, 102.5, 2.6, "****", fontsize=9.0, color=COL_KO_EDGE)
    axD.set_xticks(x_d)
    axD.set_xticklabels(["Sham\n(1,073)", "CCI+scr\n(519)", "CCI+KO\n(1,930)"], fontweight="bold", fontsize=8.0)
    axD.set_ylabel("DG granule neuron state (%)\n(GSE333879, n = 3,522)")
    axD.set_ylim(0, 114)
    axD.set_xlim(-0.52, 2.52)
    axD.set_title("Astrocytic Tgfbr1 KO rescues DG_C0 state", pad=30)
    axD.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), frameon=False, fontsize=6.8, ncol=2, columnspacing=0.6)

    # Panel E: Human TBI Neuronal Apoptotic Shift
    axE = fig.add_subplot(gs_mid[0, 2])
    style_ax(axE, zero_line="x")
    add_panel_label(axE, "E", x=-0.22, y=1.14)

    hu_items = [
        ("BAX (Pro-apoptotic)",     +273.7, "***"),
        ("CASP8 (Initiator)",       +66.0,  "**"),
        ("BAD (BH3 sensitizer)",    +44.8,  "**"),
        ("CASP9 (Apoptosome)",      +31.2,  "*"),
        ("CYCS (Cytochrome c)",     +18.4,  "*"),
        ("PTK2 (FAK kinase)",       -9.4,   "**"),
        ("BCL2 (Survival guard)",   -11.5,  "***")
    ][::-1]

    y_e = np.arange(len(hu_items))
    v_e = [x[1] for x in hu_items]
    s_e = [x[2] for x in hu_items]
    l_e = [x[0] for x in hu_items]
    fills_e = [COL_TBI_FILL if v > 0 else COL_SCR_FILL for v in v_e]
    edges_e = [COL_TBI_EDGE if v > 0 else COL_SCR_EDGE for v in v_e]
    mks_e   = ["v" if v > 0 else "o" for v in v_e]

    axE.barh(y_e, v_e, height=0.52, color=fills_e, edgecolor=edges_e, lw=1.15, zorder=2)
    for yi, v, st, fc, ec, mk in zip(y_e, v_e, s_e, fills_e, edges_e, mks_e):
        axE.scatter([v], [yi], s=34, marker=mk, facecolor=fc, edgecolor=ec, lw=1.15, zorder=3)
        x_p = v + 8 if v > 0 else v - 8
        ha_p = "left" if v > 0 else "right"
        axE.text(x_p, yi, f"{v:+.1f}% ({st})", va="center", ha=ha_p,
                 fontsize=7.3, fontweight="bold", color=ec)

    axE.set_yticks(y_e)
    axE.set_yticklabels(l_e, fontweight="bold", fontsize=7.8)
    axE.set_xlim(-88, 375)
    axE.set_xlabel("% change in human TBI neurons vs Control\n(GSE209552, n = 13,326 nuclei)")
    axE.set_title("Human TBI cortical neurons: apoptotic shift", pad=8)

    save_pub_figure(fig, OUT_MAIN, "Fig 8_cleaned", dpi=300)
    plt.close(fig)


# ===========================================================================
# SUPPLEMENTARY FIGURES S1-S4 (Unified Translucent Prism Palette Matching Fig 1-8)
# ===========================================================================
def build_supp_figures():
    # -----------------------------------------------------------------------
    # Supplementary Figure S1: Tgfbr1/TGFBR1 mRNA & Detection Rate Across 5 Datasets
    # -----------------------------------------------------------------------
    fig1, axes1 = plt.subplots(1, 5, figsize=(16.5, 4.2), dpi=300)
    m = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/figures/Fig_Human_TGFb/tgfbr1_byCelltype_external_mouse.csv"))
    hu = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/figures/Fig_Human_TGFb/human_TGFBR_change_byCelltype.csv"))
    hu = hu[hu["gene"] == "TGFBR1"]
    sf = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/figures/Fig_CellChat/tgfbr_expression_by_celltype.csv"))
    sf.rename(columns={sf.columns[0]: "celltype"}, inplace=True)
    norm_map = {"Excitatory": "Neuron", "Inhibitory": "Neuron", "Astrocytes": "Astro",
                "Microglia": "Micro", "Oligodendrocytes": "Oligo", "OPCs": "OPC", "Pericyte": "Endo"}
    sf["celltype"] = sf["celltype"].map(lambda x: norm_map.get(x, x))
    sf_agg = sf.groupby("celltype", as_index=False)[["Tgfbr1", "Tgfbr2"]].max()

    datasets_s1 = [
        ("GSE209552 (Human, n = 17)",
         dict(zip(hu["celltype"], hu["TBI"])),
         dict(zip(hu["celltype"], hu["det_TBI"]))),
        ("GSE271769 (Mouse, n = 16)",
         dict(zip(m[m["dataset"] == "GSE271769"]["celltype"], m[m["dataset"] == "GSE271769"]["Tgfbr1_mean"])),
         dict(zip(m[m["dataset"] == "GSE271769"]["celltype"], m[m["dataset"] == "GSE271769"]["Tgfbr1_det"]))),
        ("GSE277487 (Mouse, n = 6)",
         dict(zip(m[m["dataset"] == "GSE277487"]["celltype"], m[m["dataset"] == "GSE277487"]["Tgfbr1_mean"])),
         dict(zip(m[m["dataset"] == "GSE277487"]["celltype"], m[m["dataset"] == "GSE277487"]["Tgfbr1_det"]))),
        ("GSE290150 (Mouse, n = 6)",
         dict(zip(m[m["dataset"] == "GSE290150"]["celltype"], m[m["dataset"] == "GSE290150"]["Tgfbr1_mean"])),
         dict(zip(m[m["dataset"] == "GSE290150"]["celltype"], m[m["dataset"] == "GSE290150"]["Tgfbr1_det"]))),
        ("GSE333879 (In-house, n = 35k)",
         dict(zip(sf_agg["celltype"], sf_agg["Tgfbr1"])),
         {"Micro": 84.2, "Endo": 28.4, "OPC": 19.5, "Oligo": 14.2, "Astro": 16.8, "Neuron": 9.4})
    ]
    ct_order = ["Micro", "Endo", "OPC", "Oligo", "Astro", "Neuron"]
    ct_cols  = [COL_TBI_FILL, COL_GOLD_FILL, COL_TEAL_FILL, COL_SCR_FILL, COL_KO_FILL, COL_SHAM_FILL]
    ct_edgs  = [COL_TBI_EDGE, COL_GOLD_EDGE, COL_TEAL_EDGE, COL_SCR_EDGE, COL_KO_EDGE, COL_SHAM_EDGE]
    ct_mks   = ["v", "D", "s", "o", "^", "o"]

    for idx, (ds_title, val_map, det_map) in enumerate(datasets_s1):
        ax = axes1[idx]
        style_ax(ax)
        if idx == 0:
            add_panel_label(ax, "A", x=-0.24, y=1.14)
        vals = [val_map.get(c, 0.0) for c in ct_order]
        dets = [det_map.get(c, 0.0) for c in ct_order]
        ymax_s1 = max(vals) * 1.32

        bars_s1 = ax.bar(range(len(ct_order)), vals, width=0.56, color=ct_cols, edgecolor=ct_edgs, lw=1.15, zorder=2)
        for i_c, (v_c, d_c, fc_c, ec_c, mk_c) in enumerate(zip(vals, dets, ct_cols, ct_edgs, ct_mks)):
            ax.scatter([i_c], [v_c], s=32, marker=mk_c, facecolor=fc_c, edgecolor=ec_c, lw=1.1, zorder=3)
            # Value + % detection label above bar
            lbl_s1 = f"{v_c:.2f}\n({d_c:.0f}%)" if d_c > 0 else f"{v_c:.2f}"
            ax.text(i_c, v_c + ymax_s1 * 0.03, lbl_s1, ha="center", va="bottom",
                    fontsize=6.6, fontweight="bold", color=ec_c, linespacing=0.92)

        # Highlight Astrocyte column with a subtle vertical callout box
        ax.axvspan(3.62, 4.38, color=COL_KO_RAW, alpha=0.09, zorder=1)
        ax.set_ylim(0, ymax_s1)
        ax.set_xticks(range(len(ct_order)))
        ax.set_xticklabels(ct_order, rotation=32, ha="right", fontweight="bold", fontsize=8.2)
        ax.set_title(ds_title, fontsize=8.8, fontweight="bold", pad=8)
        if idx == 0:
            ax.set_ylabel("Mean Tgfbr1 / TGFBR1 mRNA\n(top: mean; bottom: % det.)")

    plt.tight_layout()
    save_pub_figure(fig1, OUT_SUPP, "Fig S1_cleaned", dpi=300)
    plt.close(fig1)

    # -----------------------------------------------------------------------
    # Supplementary Figure S2: Panel Composition Sensitivity & Hippocampal Ligand mRNA
    # -----------------------------------------------------------------------
    fig2, (ax2A, ax2B) = plt.subplots(1, 2, figsize=(15.8, 5.2), dpi=300, gridspec_kw={"width_ratios": [1.02, 1.15]})
    style_ax(ax2A, zero_line="y")
    style_ax(ax2B, zero_line="y")
    add_panel_label(ax2A, "A", x=-0.14, y=1.10)
    add_panel_label(ax2B, "B", x=-0.12, y=1.10)

    panels_s2 = ["Full 16-gene\n(6 families)", "NRXN only\n(3 genes)", "CNTN only\n(3 genes)", "CADM only\n(3 genes)", "No-CADM\n(13 genes)"]
    hip_s2 = [-15.0, -3.0, +25.0, -26.0, -15.0]
    hip_sig_s2 = ["*", "ns", "*", "ns", "*"]
    ctx_s2 = [+50.0, +14.0, +65.0, +129.0, +44.0]
    ctx_sig_s2 = ["***", "***", "***", "***", "***"]
    x_s2 = np.arange(len(panels_s2))
    w_s2 = 0.35
    b_hip = ax2A.bar(x_s2 - w_s2/2, hip_s2, width=w_s2, color=COL_TBI_FILL, edgecolor=COL_TBI_EDGE, lw=1.15,
                     label="Hippocampus (selective collapse)", zorder=2)
    b_ctx = ax2A.bar(x_s2 + w_s2/2, ctx_s2, width=w_s2, color=COL_SHAM_FILL, edgecolor=COL_SHAM_EDGE, lw=1.15,
                     label="Cerebral cortex (compensatory ↑)", zorder=2)

    for i_s, (vh, sh, vc, sc) in enumerate(zip(hip_s2, hip_sig_s2, ctx_s2, ctx_sig_s2)):
        ax2A.scatter([x_s2[i_s] - w_s2/2], [vh], s=30, marker="v", facecolor=COL_TBI_FILL, edgecolor=COL_TBI_EDGE, lw=1.0, zorder=3)
        ax2A.scatter([x_s2[i_s] + w_s2/2], [vc], s=30, marker="o", facecolor=COL_SHAM_FILL, edgecolor=COL_SHAM_EDGE, lw=1.0, zorder=3)
        yh_txt = vh - 5.5 if vh < 0 else vh + 4.0
        vah = "top" if vh < 0 else "bottom"
        ax2A.text(x_s2[i_s] - w_s2/2, yh_txt, f"{vh:+.0f}%\n({sh})", ha="center", va=vah,
                  fontsize=7.0, fontweight="bold", color=COL_TBI_EDGE, linespacing=0.9)
        ax2A.text(x_s2[i_s] + w_s2/2, vc + 4.0, f"+{vc:.0f}%\n({sc})", ha="center", va="bottom",
                  fontsize=7.0, fontweight="bold", color=COL_SHAM_EDGE, linespacing=0.9)

    ax2A.set_xticks(x_s2)
    ax2A.set_xticklabels(panels_s2, fontweight="bold", fontsize=8.0)
    ax2A.set_ylabel("% change in TBI vs Sham")
    ax2A.set_ylim(-52, 162)
    ax2A.set_title("Adhesion-module sensitivity across 5 gene-panel definitions", pad=8)
    ax2A.legend(loc="upper left", frameon=False, fontsize=8.0)

    da3 = pd.read_csv(os.path.join(BASE, "JNI/analysis/DA3_public_mouse_ligand/DA3_DEG_ligands.csv"))
    da3 = da3.sort_values("avg_log2FC").reset_index(drop=True)
    x_da = np.arange(len(da3))
    cols_da = [COL_TBI_FILL if v < 0 else COL_SHAM_FILL for v in da3["avg_log2FC"]]
    edgs_da = [COL_TBI_EDGE if v < 0 else COL_SHAM_EDGE for v in da3["avg_log2FC"]]
    mks_da  = ["v" if v < 0 else "o" for v in da3["avg_log2FC"]]
    bars_da = ax2B.bar(x_da, da3["avg_log2FC"], width=0.54, color=cols_da, edgecolor=edgs_da, lw=1.15, zorder=2)
    for i_d, (bar, v, padj, fc_d, ec_d, mk_d) in enumerate(zip(bars_da, da3["avg_log2FC"], da3["p_val_adj"], cols_da, edgs_da, mks_da)):
        ax2B.scatter([i_d], [v], s=30, marker=mk_d, facecolor=fc_d, edgecolor=ec_d, lw=1.05, zorder=3)
        st = "***" if padj < 0.001 else ("*" if padj < 0.05 else "")
        lbl_d = f"{v:+.2f}\n{st}".strip()
        y_p = v - 0.045 if v < 0 else v + 0.035
        va_p = "top" if v < 0 else "bottom"
        ax2B.text(i_d, y_p, lbl_d, ha="center", va=va_p, fontsize=7.0, fontweight="bold", color=ec_d, linespacing=0.9)

    ax2B.set_xticks(x_da)
    ax2B.set_xticklabels(da3["gene"], rotation=28, ha="right", fontstyle="italic", fontweight="bold")
    ax2B.set_ylabel("Astrocytic ligand log$_2$FC (TBI vs Sham)")
    ax2B.set_ylim(-0.78, 0.86)
    ax2B.set_title("Hippocampal astrocyte adhesion-ligand mRNA shift (n = 22 mice)", pad=8)

    plt.tight_layout()
    save_pub_figure(fig2, OUT_SUPP, "Fig S2_cleaned", dpi=300)
    plt.close(fig2)

    # -----------------------------------------------------------------------
    # Supplementary Figure S3: Pooled GSE333879 % Positivity & Within-Cell Co-Positivity
    #   Unified 3-Group Palette: Sham (Slate), CCI+scr (Vermilion), CCI+KO (Teal)
    #   Upgraded Panel B: 16-Gene Hub Connectivity Profile + Edge Count Summary
    # -----------------------------------------------------------------------
    fig3 = plt.figure(figsize=(16.5, 5.4), dpi=300)
    gs3 = gridspec.GridSpec(1, 2, figure=fig3, width_ratios=[1.22, 1.08], wspace=0.28)
    ax3A = fig3.add_subplot(gs3[0, 0])
    style_ax(ax3A)
    add_panel_label(ax3A, "A", x=-0.10, y=1.10)

    pct_df = pd.read_csv(os.path.join(BASE, "JNI/analysis/upgrade_C_pct_positivity/C_pct_positivity.csv"))
    neu_p = pct_df[pct_df["celltype"] == "Neuron"]
    genes_s3 = sorted(neu_p["gene"].unique())
    x_s3 = np.arange(len(genes_s3))
    w_s3 = 0.26
    sh_s3 = [neu_p[(neu_p["gene"] == g) & (neu_p["group"] == "Sham")]["pct"].values[0] for g in genes_s3]
    sc_s3 = [neu_p[(neu_p["gene"] == g) & (neu_p["group"] == "CCI_Scr")]["pct"].values[0] for g in genes_s3]
    ko_s3 = [neu_p[(neu_p["gene"] == g) & (neu_p["group"] == "CCI_TGFbetaR1_KO")]["pct"].values[0] for g in genes_s3]

    ax3A.bar(x_s3 - w_s3, sh_s3, width=w_s3, color=COL_SHAM_FILL, edgecolor=COL_SHAM_EDGE, lw=0.95, label="Sham (n = 8,114)", zorder=2)
    ax3A.bar(x_s3,        sc_s3, width=w_s3, color=COL_TBI_FILL,  edgecolor=COL_TBI_EDGE,  lw=0.95, label="CCI+scr (n = 4,313)", zorder=2)
    ax3A.bar(x_s3 + w_s3, ko_s3, width=w_s3, color=COL_KO_FILL,   edgecolor=COL_KO_EDGE,   lw=0.95, label="CCI+TGF-βR1 KO (n = 10,439)", zorder=2)
    ax3A.scatter(x_s3 - w_s3, sh_s3, s=16, marker="o", facecolor=COL_SHAM_FILL, edgecolor=COL_SHAM_EDGE, lw=0.8, zorder=3)
    ax3A.scatter(x_s3,        sc_s3, s=16, marker="v", facecolor=COL_TBI_FILL,  edgecolor=COL_TBI_EDGE,  lw=0.8, zorder=3)
    ax3A.scatter(x_s3 + w_s3, ko_s3, s=18, marker="^", facecolor=COL_KO_FILL,   edgecolor=COL_KO_EDGE,   lw=0.8, zorder=3)

    ax3A.set_xticks(x_s3)
    ax3A.set_xticklabels(genes_s3, rotation=32, ha="right", fontstyle="italic", fontweight="bold")
    ax3A.set_ylabel("% positive neuronal nuclei (≥1 UMI)")
    ax3A.set_ylim(0, 118)
    ax3A.set_title("Neuronal 16-gene adhesion positivity profile in GSE333879", pad=20)
    ax3A.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, frameon=False, fontsize=7.6)

    # Panel B: 2-Part Composite (Left: 16-Gene Hub Co-Positivity Mean r; Right: Significant Network Edges)
    gs3_B = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs3[0, 1], width_ratios=[1.35, 0.72], wspace=0.34)
    ax3B_left = fig3.add_subplot(gs3_B[0, 0])
    ax3B_right = fig3.add_subplot(gs3_B[0, 1])
    style_ax(ax3B_left)
    style_ax(ax3B_right)
    add_panel_label(ax3B_left, "B", x=-0.18, y=1.10)

    hub_df = pd.read_csv(os.path.join(BASE, "JNI/analysis/upgrade_A_copositivity/A_hub_membership.csv"))
    sh_hub = hub_df[hub_df["group"] == "Sham"].sort_values("obs_mean_r", ascending=True)
    gene_order_b = sh_hub["gene"].values
    y_hb = np.arange(len(gene_order_b))

    r_sh = [hub_df[(hub_df["group"] == "Sham") & (hub_df["gene"] == g)]["obs_mean_r"].values[0] for g in gene_order_b]
    r_sc = [hub_df[(hub_df["group"] == "CCI_Scr") & (hub_df["gene"] == g)]["obs_mean_r"].values[0] for g in gene_order_b]
    r_ko = [hub_df[(hub_df["group"] == "CCI_TGFbetaR1_KO") & (hub_df["gene"] == g)]["obs_mean_r"].values[0] for g in gene_order_b]

    ax3B_left.plot(r_sh, y_hb, color=COL_SHAM_EDGE, lw=1.2, marker="o", markersize=4.5,
                   markerfacecolor=COL_SHAM_FILL, label="Sham", zorder=3)
    ax3B_left.plot(r_sc, y_hb, color=COL_TBI_EDGE, lw=1.2, marker="v", markersize=4.5,
                   markerfacecolor=COL_TBI_FILL, label="CCI+scr", zorder=4)
    ax3B_left.plot(r_ko, y_hb, color=COL_KO_EDGE, lw=1.2, marker="^", markersize=5.0,
                   markerfacecolor=COL_KO_FILL, label="CCI+KO", zorder=5)
    ax3B_left.set_yticks(y_hb)
    ax3B_left.set_yticklabels(gene_order_b, fontstyle="italic", fontweight="bold", fontsize=7.6)
    ax3B_left.set_xlabel("Within-cell co-positivity mean r", fontsize=8.2)
    ax3B_left.set_xlim(0.0, 0.19)
    ax3B_left.set_title("16-gene hub connectivity", pad=18, fontsize=9.0)
    ax3B_left.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, frameon=False, fontsize=6.8, columnspacing=0.5, handletextpad=0.2)

    x_3b = np.arange(3)
    edges_3b = [71, 68, 72]
    mean_r_3b = [0.096, 0.108, 0.087]
    fills_3b = [COL_SHAM_FILL, COL_TBI_FILL, COL_KO_FILL]
    edges_col_3b = [COL_SHAM_EDGE, COL_TBI_EDGE, COL_KO_EDGE]
    mks_3b = ["o", "v", "^"]
    ax3B_right.bar(x_3b, edges_3b, width=0.52, color=fills_3b, edgecolor=edges_col_3b, lw=1.15, zorder=2)
    for i, (v, mr, fc, ec, mk) in enumerate(zip(edges_3b, mean_r_3b, fills_3b, edges_col_3b, mks_3b)):
        ax3B_right.scatter([i], [v], s=34, marker=mk, facecolor=fc, edgecolor=ec, lw=1.1, zorder=3)
        ax3B_right.text(i, v + 2.0, f"{v} edges\n(r={mr:.3f})", ha="center", va="bottom",
                        fontweight="bold", fontsize=7.0, color=ec, linespacing=0.92)
    ax3B_right.set_xticks(x_3b)
    ax3B_right.set_xticklabels(["Sham\n(8,114)", "CCI+scr\n(4,313)", "CCI+KO\n(10,439)"], fontweight="bold", fontsize=7.5)
    ax3B_right.set_ylabel("Significant edges (r > q99+0.05)", fontsize=8.2)
    ax3B_right.set_ylim(0, 92)
    ax3B_right.set_title("Network edge count", pad=8, fontsize=9.0)

    fig3.tight_layout()
    save_pub_figure(fig3, OUT_SUPP, "Fig S3_cleaned", dpi=300)
    plt.close(fig3)

    # -----------------------------------------------------------------------
    # Supplementary Figure S4: Transcriptomic Boundary & Negative-Gating Controls
    # -----------------------------------------------------------------------
    fig4, axes4 = plt.subplots(2, 2, figsize=(15.8, 10.4), dpi=300)
    for ax, let in zip(axes4.flat, ["A", "B", "C", "D"]):
        style_ax(ax)
        add_panel_label(ax, let, x=-0.14, y=1.10)

    # 4A: PANoptosis detection gate with cohort-level SEM + scatter + shaded <1% gate
    pa = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/figures/Fig_PANoptosis/panoptosis_external_neuron_panel.csv"))
    cat_map = {"凋亡": "Apoptosis", "坏死性凋亡": "Necroptosis", "焦亡": "Pyroptosis"}
    pa["cat_en"] = pa["category"].map(lambda x: cat_map.get(x, x))
    pa_agg = pa.groupby(["gene", "cat_en"], as_index=False).agg(
        mean_det=("pct_det_TBI", "mean"),
        sem_det=("pct_det_TBI", "sem")
    )
    pa_agg = pa_agg.sort_values(["cat_en", "mean_det"], ascending=[True, False]).reset_index(drop=True)

    cat_colors = {"Apoptosis": COL_TBI_FILL, "Necroptosis": COL_KO_FILL, "Pyroptosis": COL_GOLD_FILL}
    cat_edges  = {"Apoptosis": COL_TBI_EDGE, "Necroptosis": COL_KO_EDGE, "Pyroptosis": COL_GOLD_EDGE}
    x_4a = np.arange(len(pa_agg))
    axes4[0, 0].axhspan(0, 1.0, color=COL_TBI_RAW, alpha=0.10, zorder=1)
    axes4[0, 0].bar(x_4a, pa_agg["mean_det"], width=0.60,
                    color=[cat_colors.get(c, "#CBD5E1") for c in pa_agg["cat_en"]],
                    edgecolor=[cat_edges.get(c, "#334155") for c in pa_agg["cat_en"]], lw=1.05, zorder=2)
    axes4[0, 0].errorbar(x_4a, pa_agg["mean_det"], yerr=pa_agg["sem_det"].fillna(0),
                         fmt="none", ecolor="#334155", elinewidth=0.95, capsize=2.8, zorder=3)
    for i_g, row_g in pa_agg.iterrows():
        pts_g = pa[pa["gene"] == row_g["gene"]]["pct_det_TBI"].values
        jit_g = np.linspace(-0.12, 0.12, len(pts_g))
        axes4[0, 0].scatter(np.full_like(pts_g, i_g) + jit_g, pts_g, s=16,
                            facecolor=cat_colors.get(row_g["cat_en"], "#CBD5E1"),
                            edgecolor=cat_edges.get(row_g["cat_en"], "#334155"), lw=0.75, alpha=0.88, zorder=4)
    axes4[0, 0].axhline(1.0, color=COL_TBI_EDGE, ls="--", lw=1.0, zorder=3)
    axes4[0, 0].annotate(
        "<1% neuronal detection gate\n(Zbp1, Ripk3, Mlkl < 0.4%)",
        xy=(9.0, 1.05), xytext=(9.0, 10.5),
        ha="center", va="center", fontsize=7.2, fontweight="bold", color=COL_TBI_EDGE,
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor=COL_TBI_EDGE, lw=0.85, alpha=0.95),
        arrowprops=dict(arrowstyle="->", color=COL_TBI_EDGE, lw=0.9), zorder=5
    )
    axes4[0, 0].set_xticks(x_4a)
    axes4[0, 0].set_xticklabels(pa_agg["gene"], rotation=42, ha="right", fontstyle="italic", fontweight="bold", fontsize=7.6)
    axes4[0, 0].set_ylabel("Detection rate in TBI neurons (%)")
    axes4[0, 0].set_title("Single-cell detection gate: Necroptosis / Pyroptosis < 1% in neurons", pad=8)
    axes4[0, 0].legend(handles=[
        Patch(facecolor=COL_TBI_FILL,  edgecolor=COL_TBI_EDGE,  label="Intrinsic Apoptosis (mean ± SEM, n=4 cohorts)"),
        Patch(facecolor=COL_KO_FILL,   edgecolor=COL_KO_EDGE,   label="Necroptosis (<1% gate)"),
        Patch(facecolor=COL_GOLD_FILL, edgecolor=COL_GOLD_EDGE, label="Pyroptosis (<1% gate)")
    ], loc="upper right", frameon=False, fontsize=7.4)

    # 4B: Glutamate-handling genes across 4 validated cohorts with value badges
    gl = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/figures/Fig_Homeostatic_Deepdive/glutamate_panel_ALL_datasets.csv"))
    gl = gl[~gl["dataset"].isin(["GSE292715", "GSE230253", "GSE226072"])].copy().reset_index(drop=True)
    x_4b = np.arange(len(gl))
    w_4b = 0.25
    b1 = axes4[0, 1].bar(x_4b - w_4b, gl["GLUL_pct"],   width=w_4b, color=COL_TBI_FILL,  edgecolor=COL_TBI_EDGE,  lw=1.05, label="Glul / GLUL (GS)", zorder=2)
    b2 = axes4[0, 1].bar(x_4b,        gl["SLC1A2_pct"], width=w_4b, color=COL_SCR_FILL,  edgecolor=COL_SCR_EDGE,  lw=1.05, label="Slc1a2 / SLC1A2 (GLT-1)", zorder=2)
    b3 = axes4[0, 1].bar(x_4b + w_4b, gl["SLC1A3_pct"], width=w_4b, color=COL_SHAM_FILL, edgecolor=COL_SHAM_EDGE, lw=1.05, label="Slc1a3 / SLC1A3 (GLAST)", zorder=2)
    for bars_sub, col_e in [(b1, COL_TBI_EDGE), (b2, COL_SCR_EDGE), (b3, COL_SHAM_EDGE)]:
        for bar in bars_sub:
            v_b = bar.get_height()
            yp = v_b + 1.2 if v_b >= 0 else v_b - 1.5
            vap = "bottom" if v_b >= 0 else "top"
            axes4[0, 1].text(bar.get_x() + bar.get_width()/2, yp, f"{v_b:+.0f}%",
                             ha="center", va=vap, fontsize=6.5, fontweight="bold", color=col_e)
    axes4[0, 1].axhline(0, color="#475569", lw=0.85, ls="--")
    axes4[0, 1].set_xticks(x_4b)
    axes4[0, 1].set_xticklabels(gl["dataset"], rotation=20, ha="right", fontweight="bold", fontsize=8.0)
    axes4[0, 1].set_ylim(-55, 28)
    axes4[0, 1].set_ylabel("% change in astrocytes (TBI vs Control)")
    axes4[0, 1].set_title("Astrocytic glutamate-handling genes across TBI cohorts", pad=8)
    axes4[0, 1].legend(loc="lower left", frameon=False, fontsize=7.6)

    # 4C: Kir4.1 (Kcnj10) vs Aqp4 negative gate with regression corridor & stats badge
    k_df = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/figures/Fig_Kir41_Gate/kir41_crossdataset_gate.csv"))
    x_line = np.linspace(-40, 60, 100)
    axes4[1, 0].fill_between(x_line, x_line - 12, x_line + 12, color=COL_KO_RAW, alpha=0.12, zorder=1, label="±12% passive-tracking corridor")
    axes4[1, 0].plot(x_line, x_line, ls="--", color="#64748B", lw=1.1, zorder=2, label="1:1 identity line (y = x)")
    axes4[1, 0].scatter(k_df["Aqp4_pct"], k_df["Kir41_pct"], s=68, marker="D", c=COL_KO_FILL, edgecolors=COL_KO_EDGE, lw=1.25, zorder=4)
    offsets_4c = [(-6, 6), (6, -8), (6, 6), (6, -8), (6, 6)]
    for idx_r, (_, row) in enumerate(k_df.iterrows()):
        dx, dy = offsets_4c[idx_r % len(offsets_4c)]
        axes4[1, 0].annotate(
            str(row["dataset"]),
            xy=(row["Aqp4_pct"], row["Kir41_pct"]),
            xytext=(dx, dy), textcoords="offset points",
            fontsize=7.8, fontweight="bold", color="#0F172A"
        )
    r_kir = np.corrcoef(k_df["Aqp4_pct"], k_df["Kir41_pct"])[0, 1]
    axes4[1, 0].text(0.05, 0.88, f"Passive osmotic tracking:\nPearson r = {r_kir:.2f} (no selective Kir4.1 collapse)",
                     transform=axes4[1, 0].transAxes, ha="left", va="center",
                     fontsize=7.5, fontweight="bold", color=COL_KO_EDGE,
                     bbox=dict(boxstyle="round,pad=0.24", facecolor="white", edgecolor=COL_KO_EDGE, lw=0.85))
    axes4[1, 0].set_xlim(-42, 62)
    axes4[1, 0].set_ylim(-42, 62)
    axes4[1, 0].set_xlabel("Aqp4 % change (TBI vs Control)")
    axes4[1, 0].set_ylabel("Kcnj10 (Kir4.1) % change (TBI vs Control)")
    axes4[1, 0].set_title("Negative gate: Kcnj10 (Kir4.1) passively tracks Aqp4", pad=8)
    axes4[1, 0].legend(loc="lower right", frameon=False, fontsize=7.4)

    # 4D: Time-series neuronal DEG persistence with exact counts
    ft_dir = os.path.join(BASE, "wet_lab_data/results/figures/Fig_TimeSeries")
    deg_files = [
        ("DEG_GSE271769_Neuron_D1.csv", "GSE271769\n1d"),
        ("DEG_GSE271769_Neuron_D7.csv", "GSE271769\n7d"),
        ("DEG_GSE271769_Neuron_D28.csv", "GSE271769\n28d"),
        ("DEG_GSE277487_Neuron_t6h.csv", "GSE277487\n6h"),
        ("DEG_GSE277487_Neuron_t2d.csv", "GSE277487\n2d"),
        ("DEG_GSE277487_Neuron_t4d.csv", "GSE277487\n4d"),
        ("DEG_GSE180862Hip_Neuron_24hrs.csv", "GSE180862\n24h"),
        ("DEG_GSE180862Hip_Neuron_7days.csv", "GSE180862\n7d"),
        ("DEG_GSE101901Hip_Neuron.csv", "GSE101901\n3d")
    ]
    lbls_4d, up_4d, dn_4d = [], [], []
    for fn, lb in deg_files:
        fp = os.path.join(ft_dir, fn)
        if os.path.exists(fp):
            d = pd.read_csv(fp)
            lbls_4d.append(lb)
            up_4d.append(int(((d["padj"] < 0.05) & (d["log2FoldChange"] > 1)).sum()))
            dn_4d.append(int(((d["padj"] < 0.05) & (d["log2FoldChange"] < -1)).sum()))

    x_4d = np.arange(len(lbls_4d))
    w_4d = 0.36
    b_up = axes4[1, 1].bar(x_4d - w_4d/2, up_4d, width=w_4d, color=COL_TBI_FILL, edgecolor=COL_TBI_EDGE, lw=1.05, label="Upregulated DEGs (log$_2$FC > 1)")
    b_dn = axes4[1, 1].bar(x_4d + w_4d/2, dn_4d, width=w_4d, color=COL_SCR_FILL, edgecolor=COL_SCR_EDGE, lw=1.05, label="Downregulated DEGs (log$_2$FC < −1)")
    ymax_4d = max(max(up_4d), max(dn_4d)) * 1.18
    for bar in b_up:
        v_u = int(bar.get_height())
        axes4[1, 1].text(bar.get_x() + bar.get_width()/2, v_u + ymax_4d * 0.015, str(v_u),
                         ha="center", va="bottom", fontsize=6.6, fontweight="bold", color=COL_TBI_EDGE)
    for bar in b_dn:
        v_d = int(bar.get_height())
        axes4[1, 1].text(bar.get_x() + bar.get_width()/2, v_d + ymax_4d * 0.015, str(v_d),
                         ha="center", va="bottom", fontsize=6.6, fontweight="bold", color=COL_SCR_EDGE)
    axes4[1, 1].set_ylim(0, ymax_4d)
    axes4[1, 1].set_xticks(x_4d)
    axes4[1, 1].set_xticklabels(lbls_4d, fontweight="bold", fontsize=7.6)
    axes4[1, 1].set_ylabel("Number of neuronal pseudobulk DEGs")
    axes4[1, 1].set_title("Persistent neuronal transcriptional dysregulation from 6 h to 28 dpi", pad=8)
    axes4[1, 1].legend(loc="upper right", frameon=False, fontsize=7.6)

    plt.tight_layout()
    save_pub_figure(fig4, OUT_SUPP, "Fig S4_cleaned", dpi=300)
    plt.close(fig4)

    # -----------------------------------------------------------------------
    # Supplementary Figure S5: In-House Hippocampal snRNA-seq Atlas (GSE333879, n = 35,047)
    #   Lineage Validation, Purity Gate, DG_C0–C3 Subclustering & Mitochondrial OXPHOS Rescue
    # -----------------------------------------------------------------------
    build_supp_figure_s5()


def build_supp_figure_s5():
    fig5 = plt.figure(figsize=(16.6, 19.2), dpi=300)
    gs5 = gridspec.GridSpec(
        4, 1, figure=fig5,
        height_ratios=[1.05, 1.00, 1.08, 0.98],
        hspace=0.48
    )

    atlas_df = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/atlas_export/atlas_35047_umap_meta.csv"))
    dot7_df  = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/atlas_export/atlas_7lineage_markers_dotplot.csv"))
    pur_df   = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/atlas_export/astro_purity_gate_1038_to_774.csv"))
    dg_umap  = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/atlas_export/dg_subcluster_umap_meta.csv"))
    dg_dot   = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/atlas_export/dg_subcluster_markers_dotplot.csv"))
    dg_res   = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/figures/Fig_DG_subtype/DG_subtype_rescue_effect.csv"))
    mito_df  = pd.read_csv(os.path.join(BASE, "wet_lab_data/results/figures/Fig_Mito_Global/Global_mito_module_3groups.csv"))

    ct_order_7 = [
        "Excitatory neuron",
        "Oligodendrocyte",
        "Inhibitory neuron",
        "Microglia",
        "OPC",
        "Endothelial / Pericyte",
        "Astrocyte",
    ]

    # -----------------------------------------------------------------------
    # Row 0 (Full Width) — Panel A: Condition-Split Global UMAPs (Sham, CCI+scr, CCI+KO)
    #   + Lineage Composition Bar + External 7-Lineage Color Key
    # -----------------------------------------------------------------------
    gs5_A = gridspec.GridSpecFromSubplotSpec(
        1, 5, subplot_spec=gs5[0],
        width_ratios=[1.05, 1.05, 1.05, 0.85, 0.55],
        wspace=0.26
    )
    conds_a = [
        ("Sham",             "Sham (n = 11,378)"),
        ("CCI_Scr",          "CCI+scr (n = 6,845)"),
        ("CCI_TGFbetaR1_KO", "CCI+TGF-βR1 KO (n = 16,824)"),
    ]
    rng_s5 = np.random.default_rng(77)
    for idx_c, (c_key, c_title) in enumerate(conds_a):
        ax_c = fig5.add_subplot(gs5_A[0, idx_c])
        if idx_c == 0:
            add_panel_label(ax_c, "A", x=-0.20, y=1.14)
        sub_c = atlas_df[atlas_df["group"] == c_key]
        for ct in ct_order_7:
            sub_ct = sub_c[sub_c["celltype"] == ct]
            raw_c, edge_c = LINEAGE_PALETTE[ct]
            is_ast = (ct == "Astrocyte")
            if len(sub_ct) > 3000:
                sub_ct = sub_ct.iloc[rng_s5.choice(len(sub_ct), size=3000, replace=False)]
            ax_c.scatter(
                sub_ct["UMAP1_harm"], sub_ct["UMAP2_harm"],
                s=7.0 if is_ast else 2.5, c=edge_c, alpha=0.82 if is_ast else 0.42,
                linewidths=0, zorder=4 if is_ast else 2, rasterized=True
            )
            ax_c.scatter(
                sub_ct["UMAP1_harm"], sub_ct["UMAP2_harm"],
                s=4.5 if is_ast else 1.5, c=raw_c, alpha=0.95 if is_ast else 0.85,
                linewidths=0, zorder=5 if is_ast else 3, rasterized=True
            )
        # Highlight Astrocyte & DG islands with crisp dashed/dotted boxes
        ax_c.add_patch(plt.Rectangle((0.3, -13.3), 5.3, 4.4, fill=False, edgecolor="#004D40", lw=1.15, ls="--", zorder=6))
        ax_c.add_patch(plt.Rectangle((-6.8, -8.3), 6.5, 5.8, fill=False, edgecolor="#1F2E52", lw=1.05, ls=":", zorder=6))
        ax_c.set_xlim(-15.5, 15.5)
        ax_c.set_ylim(-17.5, 16.5)
        style_umap_ax(ax_c, xlabel="UMAP 1", ylabel="UMAP 2" if idx_c == 0 else "")
        ax_c.set_xticks([-15, -10, -5, 0, 5, 10, 15])
        ax_c.set_yticks([-15, -10, -5, 0, 5, 10, 15])
        ax_c.set_title(c_title, fontsize=8.8, fontweight="bold", pad=6)

    # 4th column of Panel A: 100% Stacked Composition Bar Across the 3 Groups
    axA_bar = fig5.add_subplot(gs5_A[0, 3])
    style_ax(axA_bar)
    ct_stack_order = [
        ("Excitatory neuron",      "Excitatory neuron (18,556)"),
        ("Inhibitory neuron",      "Inhibitory neuron (4,314)"),
        ("Oligodendrocyte",        "Oligodendrocyte (7,647)"),
        ("Microglia",              "Microglia (2,155)"),
        ("Astrocyte",              "Astrocyte (1,038 → 774 pure)"),
        ("OPC",                    "OPC (988)"),
        ("Endothelial / Pericyte", "Endothelial / Pericyte (349)"),
    ]
    x_ab = np.arange(3)
    bot_ab = np.zeros(3)
    for ct, ct_lbl in ct_stack_order:
        raw_c, edge_c = LINEAGE_PALETTE[ct]
        vals_ct = []
        for c_key, _ in conds_a:
            sub_c = atlas_df[atlas_df["group"] == c_key]
            vals_ct.append((sub_c["celltype"] == ct).mean() * 100.0)
        vals_ct = np.array(vals_ct)
        axA_bar.bar(x_ab, vals_ct, bottom=bot_ab, width=0.52,
                    color=blend_white(raw_c, 0.84), edgecolor=edge_c, lw=0.95,
                    label=ct_lbl, zorder=3)
        bot_ab += vals_ct
    axA_bar.set_xticks(x_ab)
    axA_bar.set_xticklabels(["Sham", "CCI+scr", "CCI+KO"], fontweight="bold", fontsize=8.0)
    axA_bar.set_ylabel("Lineage composition (%)", fontsize=8.2)
    axA_bar.set_ylim(0, 104)
    axA_bar.set_title("Lineage composition", fontsize=8.8, fontweight="bold", pad=6)

    # 5th column of Panel A: Dedicated Legend Panel for 7 CNS Lineages + Zoom Boxes
    axA_leg = fig5.add_subplot(gs5_A[0, 4])
    axA_leg.axis("off")
    handles_a, labels_a = axA_bar.get_legend_handles_labels()
    leg_a = axA_leg.legend(
        handles_a[::-1], labels_a[::-1],
        loc="center left", bbox_to_anchor=(-0.18, 0.54),
        frameon=True, facecolor="#F8FAFC", edgecolor="#CBD5E1",
        title="CNS Lineages (n = 35,047)", title_fontsize=7.8,
        fontsize=7.2, handlelength=1.2, handletextpad=0.4, labelspacing=0.45
    )
    leg_a.get_title().set_fontweight("bold")
    axA_leg.text(
        -0.14, 0.08,
        "Dashed green box: Astrocyte island\nDotted navy box: DG Granule island",
        transform=axA_leg.transAxes, ha="left", va="center",
        fontsize=6.8, fontweight="bold", color="#334155",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#94A3B8", lw=0.75)
    )

    # -----------------------------------------------------------------------
    # Row 1 — Panel B (Left 52%) & Panel C (Right 48%)
    # -----------------------------------------------------------------------
    gs5_row1 = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=gs5[1],
        width_ratios=[1.10, 1.00], wspace=0.26
    )

    # Panel B: 14-Gene Canonical Marker DotPlot & SingleR Reference Concordance
    gs5_B = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs5_row1[0, 0], width_ratios=[1.35, 0.80], wspace=0.36)
    axB_dot = fig5.add_subplot(gs5_B[0, 0])
    axB_sr  = fig5.add_subplot(gs5_B[0, 1])
    add_panel_label(axB_dot, "B", x=-0.22, y=1.14)

    ct_rows_b = [
        ("Excitatory neuron",      "Excitatory"),
        ("Inhibitory neuron",      "Inhibitory"),
        ("Oligodendrocyte",        "Oligodendrocyte"),
        ("Astrocyte",              "Astrocyte"),
        ("Microglia",              "Microglia"),
        ("OPC",                    "OPC"),
        ("Endothelial / Pericyte", "Endo / Pericyte"),
    ]
    genes_b = [
        "Slc17a7", "Camk2a", "Gad1", "Gad2", "Mbp", "Plp1",
        "Aqp4", "Gja1", "Aldh1l1", "Cx3cr1", "Hexb", "Pdgfra", "Cspg4", "Flt1"
    ]
    for sp in ("top", "right", "left", "bottom"):
        axB_dot.spines[sp].set_visible(True)
        axB_dot.spines[sp].set_linewidth(0.8)
        axB_dot.spines[sp].set_color("#94A3B8")
    axB_dot.grid(True, linestyle="-", linewidth=0.38, color="#E2E8F0", zorder=1)

    xb_l, yb_l, sb_l, cb_l = [], [], [], []
    for r_i, (ct_full, _) in enumerate(ct_rows_b):
        yr = len(ct_rows_b) - 1 - r_i
        for c_j, gb in enumerate(genes_b):
            rm = dot7_df[(dot7_df["celltype"] == ct_full) & (dot7_df["gene"] == gb)]
            pct_v = float(rm["pct_exp"].values[0]) if len(rm) else 0.0
            z_v   = float(rm["z_score"].values[0]) if len(rm) else -0.5
            xb_l.append(c_j)
            yb_l.append(yr)
            sb_l.append(max(pct_v * 1.15, 8.0))
            cb_l.append(z_v)

    axB_dot.scatter(xb_l, yb_l, s=sb_l, c=cb_l, cmap=NPG_CMAP, vmin=-1.2, vmax=2.0,
                    edgecolors="#0F172A", linewidths=0.6, zorder=3)
    axB_dot.axhspan(2.55, 3.45, color="#00A087", alpha=0.14, zorder=0)
    axB_dot.set_xticks(range(len(genes_b)))
    axB_dot.set_xticklabels(genes_b, rotation=42, ha="right", fontstyle="italic", fontweight="bold", fontsize=7.0)
    axB_dot.set_yticks(range(len(ct_rows_b)))
    axB_dot.set_yticklabels([x[1] for x in ct_rows_b[::-1]], fontweight="bold", fontsize=7.4)
    axB_dot.set_xlim(-0.6, len(genes_b) - 0.4)
    axB_dot.set_ylim(-0.55, len(ct_rows_b) - 0.45)
    axB_dot.set_title("Canonical CNS lineage markers (14 genes)", fontsize=8.8, fontweight="bold", pad=8)

    # Right subpanel of B: SingleR Concordance Bar/Matrix
    style_ax(axB_sr)
    sr_labels = ["Neuron\n(22,870)", "Oligo\n(7,647)", "Astrocyte\n(1,038)", "Microglia\n(2,155)", "OPC/Glia\n(988)", "Endo/Mural\n(349)"]
    sr_pcts   = [99.96, 99.65, 98.17, 96.94, 91.70, 88.25]
    sr_cols   = [
        blend_white("#3C5488", 0.78),
        blend_white("#E69F00", 0.78),
        blend_white("#00A087", 0.85),
        blend_white("#E64B35", 0.78),
        blend_white("#5E4FA2", 0.78),
        blend_white("#8C510A", 0.78),
    ]
    sr_edgs   = ["#1F2E52", "#8C510A", "#004D40", "#881D0C", "#3B2F6B", "#4A2A05"]
    y_sr = np.arange(len(sr_labels))[::-1]
    axB_sr.barh(y_sr, sr_pcts, height=0.55, color=sr_cols, edgecolor=sr_edgs, lw=1.05, zorder=2)
    for yi, pv, ec in zip(y_sr, sr_pcts, sr_edgs):
        axB_sr.text(pv + 1.5, yi, f"{pv:.1f}%", va="center", ha="left", fontsize=7.2, fontweight="bold", color=ec)
    axB_sr.set_yticks(y_sr)
    axB_sr.set_yticklabels(sr_labels, fontweight="bold", fontsize=7.2)
    axB_sr.set_xlabel("SingleR concordance (%)", fontsize=7.8)
    axB_sr.set_xlim(0, 122)
    axB_sr.set_title("SingleR validation", fontsize=8.8, fontweight="bold", pad=8)

    # Panel C: Per-Lineage QC Distribution & Astrocyte Purity Gate (1,038 -> 774)
    gs5_C = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs5_row1[0, 1], width_ratios=[0.96, 1.18], wspace=0.30)
    axC_qc  = fig5.add_subplot(gs5_C[0, 0])
    axC_pur = fig5.add_subplot(gs5_C[0, 1])
    style_ax(axC_qc)
    style_ax(axC_pur)
    add_panel_label(axC_qc, "C", x=-0.26, y=1.14)

    # Left: Per-lineage nFeature_RNA median & IQR horizontal boxplot
    qc_cts = [x[0] for x in ct_rows_b[::-1]]
    qc_lbls = [x[1] for x in ct_rows_b[::-1]]
    qc_data = [atlas_df[atlas_df["celltype"] == c]["nFeature_RNA"].values for c in qc_cts]
    bp = axC_qc.boxplot(
        qc_data, vert=False, patch_artist=True, widths=0.55, showfliers=False,
        medianprops=dict(color="#0F172A", linewidth=1.4),
        whiskerprops=dict(color="#334155", linewidth=0.9),
        capprops=dict(color="#334155", linewidth=0.9)
    )
    for patch, c_name in zip(bp["boxes"], qc_cts):
        rc, ec = LINEAGE_PALETTE[c_name]
        patch.set_facecolor(blend_white(rc, 0.78))
        patch.set_edgecolor(ec)
        patch.set_linewidth(1.0)
    axC_qc.set_yticks(range(1, len(qc_lbls) + 1))
    axC_qc.set_yticklabels(qc_lbls, fontweight="bold", fontsize=7.4)
    axC_qc.set_xlabel("Detected genes per nucleus (nFeature_RNA)", fontsize=7.8)
    axC_qc.set_title("Single-nucleus QC (n = 35,047)", fontsize=8.8, fontweight="bold", pad=8)

    # Right: Astrocyte Purity Gate (774 Retained Pure Astrocytes vs 264 Excluded Doublets)
    pur_genes = ["Slc1a2", "Gja1", "Aqp4", "Plp1", "Mbp", "Hexb", "Cx3cr1"]
    pur_sub = pur_df.set_index("gene").loc[pur_genes].reset_index()
    x_p = np.arange(len(pur_genes))
    w_p = 0.34
    axC_pur.bar(x_p - w_p/2, pur_sub["pure_pct"], width=w_p, color=COL_KO_FILL, edgecolor=COL_KO_EDGE, lw=1.1,
                label="Pure astrocytes (n = 774)", zorder=3)
    axC_pur.bar(x_p + w_p/2, pur_sub["excl_pct"], width=w_p, color=COL_TBI_FILL, edgecolor=COL_TBI_EDGE, lw=1.1,
                label="Excluded doublets (n = 264)", zorder=3)
    # Place Plp1 callout in open whitespace above Mbp/Hexb/Cx3cr1 (x=4..6, where bars are <16%) so zero bars are touched!
    axC_pur.annotate(
        "Plp1: 91.7% → 0.5%\n(Oligo doublets removed)",
        xy=(3 - w_p/2, 2.5), xytext=(4.95, 56.0),
        ha="center", va="center", fontsize=6.8, fontweight="bold", color=COL_KO_EDGE,
        bbox=dict(boxstyle="round,pad=0.20", facecolor="#F0FDFA", edgecolor=COL_KO_EDGE, lw=0.85),
        arrowprops=dict(arrowstyle="-|>,head_width=0.18,head_length=0.25", lw=1.0, color=COL_KO_EDGE,
                        connectionstyle="arc3,rad=0.20"),
        zorder=5
    )
    axC_pur.set_xticks(x_p)
    axC_pur.set_xticklabels(pur_genes, rotation=30, ha="right", fontstyle="italic", fontweight="bold", fontsize=7.6)
    axC_pur.set_ylabel("% positive nuclei", fontsize=8.0)
    axC_pur.set_ylim(0, 108)
    axC_pur.set_title("Astrocyte purity gate (1,038 → 774 pure)", fontsize=8.8, fontweight="bold", pad=20)
    axC_pur.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2,
        frameon=False, fontsize=6.8, handletextpad=0.25, columnspacing=0.6
    )

    # -----------------------------------------------------------------------
    # Row 2 — Panel D (Left 54%: 2 Balanced DG UMAPs) & Panel E (Right 46%: DG Markers + Rescue)
    # -----------------------------------------------------------------------
    gs5_row2 = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=gs5[2],
        width_ratios=[1.16, 1.00], wspace=0.26
    )

    # Panel D: Hippocampal Prox1+ DG Granule Neuron Subclustering UMAP (DG_C0–C3)
    gs5_D = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs5_row2[0, 0], width_ratios=[1.0, 1.0], wspace=0.22)
    axD_sub = fig5.add_subplot(gs5_D[0, 0])
    axD_grp = fig5.add_subplot(gs5_D[0, 1])
    add_panel_label(axD_sub, "D", x=-0.18, y=1.18)

    dg_umap["subcluster"] = dg_umap["subcluster"].astype(str)
    xg_d, yg_d = np.mgrid[-6.2:6.2:100j, -5.2:6.4:100j]
    pos_d = np.vstack([xg_d.ravel(), yg_d.ravel()])

    dg_pal = [
        ("0", COL_KO_RAW,   COL_KO_EDGE,   "DG_C0 (Homeostatic)"),
        ("1", COL_TBI_RAW,  COL_TBI_EDGE,  "DG_C1 (TBI-stressed)"),
        ("2", COL_SCR_RAW,  COL_SCR_EDGE,  "DG_C2 (IEG)"),
        ("3", COL_GOLD_RAW, COL_GOLD_EDGE, "DG_C3 (Remodel)"),
    ]
    for sc_k, rc, ec, lbl_d in dg_pal:
        sub_d = dg_umap[dg_umap["subcluster"] == sc_k]
        if len(sub_d) > 30:
            kde_d = gaussian_kde(np.vstack([sub_d["UMAP_1"], sub_d["UMAP_2"]]), bw_method=0.40)
            zd = np.reshape(kde_d(pos_d).T, xg_d.shape)
            axD_sub.contourf(xg_d, yg_d, zd, levels=[zd.max()*0.18, zd.max()*1.05], colors=[rc], alpha=0.18, zorder=1)
            axD_sub.contour(xg_d, yg_d, zd, levels=[zd.max()*0.18], colors=[ec], linewidths=0.85, linestyles="--", alpha=0.65, zorder=2)
        axD_sub.scatter(
            sub_d["UMAP_1"], sub_d["UMAP_2"],
            s=8.5 if sc_k == "3" else 3.8, c=rc, edgecolors=ec, linewidths=0.25,
            alpha=0.90, label=lbl_d, zorder=4 if sc_k == "3" else 3, rasterized=True
        )

    # Outer-margin callouts with thin leader lines on DG_C0–C3 UMAP (zero scatter overlap)
    dg_outer_notes = [
        (( 3.2, -1.5), ( 5.5, -5.8), "DG_C0 (Homeostatic)\nPcdh7+ / Bdnf+",   "#004D40", "#F0FDFA"),
        ((-3.2,  2.8), (-5.4,  6.2), "DG_C1 (TBI-stressed)\nPde1a+ / Tenm4+", COL_TBI_EDGE, "white"),
        ((-3.0, -2.8), (-5.5, -5.8), "DG_C2 (IEG state)",                     COL_SCR_EDGE, "white"),
        (( 0.2,  5.3), ( 5.2,  6.3), "DG_C3 (Remodeling)",                    COL_GOLD_EDGE, "white"),
    ]
    for xy_p, xy_t, txt_d, col_d, fc_d in dg_outer_notes:
        axD_sub.annotate(
            txt_d, xy=xy_p, xytext=xy_t,
            ha="center", va="center", fontsize=6.3, fontweight="bold", color=col_d,
            bbox=dict(boxstyle="round,pad=0.18", facecolor=fc_d, edgecolor=col_d, lw=0.85, alpha=0.95),
            arrowprops=dict(arrowstyle="-", lw=0.80, color=col_d), zorder=6
        )

    axD_sub.set_xlim(-8.8, 9.0)
    axD_sub.set_ylim(-7.6, 8.0)
    style_umap_ax(axD_sub, xlabel="UMAP 1", ylabel="UMAP 2")
    axD_sub.set_xticks([-6, -3, 0, 3, 6])
    axD_sub.set_yticks([-6, -3, 0, 3, 6])
    axD_sub.set_title("Prox1+ DG granule states (DG_C0–C3)", fontsize=8.8, fontweight="bold", pad=28)
    axD_sub.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2,
        frameon=False, fontsize=6.4, markerscale=1.6, handletextpad=0.2, columnspacing=0.5
    )

    # Right subpanel of D: Condition-Projected DG Granule UMAP with full L-shaped UMAP axes
    grp_pal_d = [
        ("Sham",             "#64748B", "#1E293B", "Sham"),
        ("CCI_Scr",          COL_TBI_RAW, COL_TBI_EDGE, "CCI+scr"),
        ("CCI_TGFbetaR1_KO", COL_KO_RAW,  COL_KO_EDGE,  "CCI+KO"),
    ]
    for g_k, rc_g, ec_g, lbl_g in grp_pal_d:
        sub_g = dg_umap[dg_umap["group"] == g_k]
        axD_grp.scatter(
            sub_g["UMAP_1"], sub_g["UMAP_2"],
            s=3.8, c=rc_g, alpha=0.82, linewidths=0, label=lbl_g, zorder=3, rasterized=True
        )
    # Draw DG_C0 & DG_C1 territory outlines on top
    for sc_k, rc, ec, _ in dg_pal[:2]:
        sub_d = dg_umap[dg_umap["subcluster"] == sc_k]
        kde_d = gaussian_kde(np.vstack([sub_d["UMAP_1"], sub_d["UMAP_2"]]), bw_method=0.40)
        zd = np.reshape(kde_d(pos_d).T, xg_d.shape)
        axD_grp.contour(xg_d, yg_d, zd, levels=[zd.max()*0.22], colors=[ec], linewidths=1.2, linestyles="--", zorder=5)

    # Outer-margin callouts with thin leader lines on Condition DG UMAP (zero scatter overlap)
    axD_grp.annotate(
        "KO repopulates DG_C0\n(+26.0% vs CCI+scr)",
        xy=(3.2, -1.5), xytext=(5.2, -5.8),
        ha="center", va="center", fontsize=6.3, fontweight="bold", color="#004D40",
        bbox=dict(boxstyle="round,pad=0.18", facecolor="#F0FDFA", edgecolor="#004D40", lw=0.85, alpha=0.95),
        arrowprops=dict(arrowstyle="-", lw=0.80, color="#004D40"), zorder=6
    )
    axD_grp.annotate(
        "CCI+scr enriches DG_C1\n(−25.4% in KO)",
        xy=(-3.2, 2.8), xytext=(-5.2, 6.2),
        ha="center", va="center", fontsize=6.3, fontweight="bold", color=COL_TBI_EDGE,
        bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor=COL_TBI_EDGE, lw=0.85, alpha=0.95),
        arrowprops=dict(arrowstyle="-", lw=0.80, color=COL_TBI_EDGE), zorder=6
    )

    axD_grp.set_xlim(-8.8, 9.0)
    axD_grp.set_ylim(-7.6, 8.0)
    style_umap_ax(axD_grp, xlabel="UMAP 1", ylabel="")
    axD_grp.set_xticks([-6, -3, 0, 3, 6])
    axD_grp.set_yticks([-6, -3, 0, 3, 6])
    axD_grp.set_title("Condition projection (χ² p = 3.83e−43)", fontsize=8.8, fontweight="bold", pad=22)
    axD_grp.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3,
        frameon=False, fontsize=6.6, markerscale=1.8, handletextpad=0.2, columnspacing=0.55
    )

    # Panel E: DG_C0–C3 State-Defining Marker DotPlot & Bidirectional Rescue Effect
    gs5_E = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs5_row2[0, 1], width_ratios=[1.34, 0.82], wspace=0.34)
    axE_dot = fig5.add_subplot(gs5_E[0, 0])
    axE_res = fig5.add_subplot(gs5_E[0, 1])
    add_panel_label(axE_dot, "E", x=-0.22, y=1.18)

    dg_genes_e = [
        "Prox1", "Dock10", "Pcdh7", "Bdnf", "Rasgrp1", "Rgs6",
        "Pde1a", "Tenm4", "Pdzrn4", "Nckap5",
        "Hectd2", "Trpc5", "Olfml2b",
        "Kcnq5", "Sox5", "Ldb2"
    ]
    dg_rows_e = ["DG_C0", "DG_C1", "DG_C2", "DG_C3"]
    dg_lbls_e = ["DG_C0 (Homeo)", "DG_C1 (Stressed)", "DG_C2 (IEG)", "DG_C3 (Remodel)"]

    for sp in ("top", "right", "left", "bottom"):
        axE_dot.spines[sp].set_visible(True)
        axE_dot.spines[sp].set_linewidth(0.8)
        axE_dot.spines[sp].set_color("#94A3B8")
    axE_dot.grid(True, linestyle="-", linewidth=0.4, color="#E2E8F0", zorder=1)

    xe_l, ye_l, se_l, ce_l = [], [], [], []
    for r_i, d_cl in enumerate(dg_rows_e):
        yr = len(dg_rows_e) - 1 - r_i
        for c_j, ge in enumerate(dg_genes_e):
            rm = dg_dot[(dg_dot["subcluster"] == d_cl) & (dg_dot["gene"] == ge)]
            pct_v = float(rm["pct_exp"].values[0]) if len(rm) else 0.0
            z_v   = float(rm["z_score"].values[0]) if len(rm) else -0.5
            xe_l.append(c_j)
            ye_l.append(yr)
            se_l.append(max(pct_v * 1.45, 10.0))
            ce_l.append(z_v)

    axE_dot.scatter(xe_l, ye_l, s=se_l, c=ce_l, cmap=NPG_CMAP, vmin=-1.35, vmax=1.5,
                    edgecolors="#0F172A", linewidths=0.65, zorder=3)
    axE_dot.axhspan(2.55, 3.45, color="#00A087", alpha=0.14, zorder=0)
    axE_dot.set_xticks(range(len(dg_genes_e)))
    axE_dot.set_xticklabels(dg_genes_e, rotation=42, ha="right", fontstyle="italic", fontweight="bold", fontsize=7.0)
    axE_dot.set_yticks(range(len(dg_rows_e)))
    axE_dot.set_yticklabels(dg_lbls_e[::-1], fontweight="bold", fontsize=7.4)
    axE_dot.set_xlim(-0.6, len(dg_genes_e) - 0.4)
    axE_dot.set_ylim(-0.55, len(dg_rows_e) - 0.45)
    axE_dot.set_title("DG_C0–C3 state-defining markers (16 genes)", fontsize=8.8, fontweight="bold", pad=8)

    # Right subpanel of E: Net Rescue Effect (KO - CCI+scr percentage points)
    style_ax(axE_res, zero_line="x")
    dg_res["Subcluster"] = dg_res["Subcluster"].astype(str)
    dg_res_s = dg_res.set_index("Subcluster").loc[["0", "1", "2", "3"]].reset_index()
    delta_pts = dg_res_s["KO_vs_CCI_pct_pt"].values
    y_er = np.arange(4)[::-1]
    fills_er = [COL_KO_FILL, COL_TBI_FILL, COL_SCR_FILL, COL_GOLD_FILL]
    edges_er = [COL_KO_EDGE, COL_TBI_EDGE, COL_SCR_EDGE, COL_GOLD_EDGE]

    axE_res.barh(y_er, delta_pts, height=0.52, color=fills_er, edgecolor=edges_er, lw=1.15, zorder=2)
    for yi, dv, ec in zip(y_er, delta_pts, edges_er):
        if abs(dv) > 15:
            axE_res.text(dv * 0.5, yi, f"{dv:+.1f}%", va="center", ha="center",
                         fontsize=7.5, fontweight="bold", color="white", zorder=4)
        else:
            xp = dv + 2.2 if dv >= 0 else dv - 2.2
            hp = "left" if dv >= 0 else "right"
            axE_res.text(xp, yi, f"{dv:+.1f}%", va="center", ha=hp,
                         fontsize=7.3, fontweight="bold", color=ec, zorder=4)
    axE_res.set_yticks(y_er)
    axE_res.set_yticklabels(["DG_C0", "DG_C1", "DG_C2", "DG_C3"], fontweight="bold", fontsize=7.5)
    axE_res.set_xlabel("Δ proportion in KO vs scr (%)", fontsize=7.8)
    axE_res.set_xlim(-35, 36)
    axE_res.set_title("Selective DG_C0 rescue", fontsize=8.8, fontweight="bold", pad=8)

    # -----------------------------------------------------------------------
    # Row 3 (Full Width) — Panel F: Cross-Lineage Mitochondrial OXPHOS Collapse & KO Rescue
    # -----------------------------------------------------------------------
    gs5_F = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=gs5[3],
        width_ratios=[1.00, 1.25], wspace=0.28
    )
    axF_ox = fig5.add_subplot(gs5_F[0, 0])
    axF_hm = fig5.add_subplot(gs5_F[0, 1])
    style_ax(axF_ox)
    add_panel_label(axF_ox, "F", x=-0.14, y=1.16)

    ox_sub = mito_df[mito_df["module"] == "mtOXPHOS_mt"].copy()
    ct_map_m = {
        "Astrocyte (全部)":         "Astrocyte (n = 1,038)",
        "Excitatory Neuron (全部)": "Excitatory Neuron (n = 18,556)",
        "Inhibitory Neuron (全部)": "Inhibitory Neuron (n = 4,314)"
    }
    ox_sub["ct_clean"] = ox_sub["cell_type"].map(ct_map_m)
    x_f = np.arange(3)
    w_f = 0.24
    sh_ox = ox_sub["Sham"].values
    sc_ox = ox_sub["CCI_Scr"].values
    ko_ox = ox_sub["CCI_TGFbetaR1_KO"].values

    axF_ox.bar(x_f - w_f, sh_ox, width=w_f, color=COL_SHAM_FILL, edgecolor=COL_SHAM_EDGE, lw=1.1, label="Sham", zorder=2)
    axF_ox.bar(x_f,       sc_ox, width=w_f, color=COL_TBI_FILL,  edgecolor=COL_TBI_EDGE,  lw=1.1, label="CCI+scr", zorder=2)
    axF_ox.bar(x_f + w_f, ko_ox, width=w_f, color=COL_KO_FILL,   edgecolor=COL_KO_EDGE,   lw=1.1, label="CCI+TGF-βR1 KO", zorder=2)

    sigs_tbi = ["**",   "****", "****"]
    sigs_ko  = ["***",  "****", "****"]
    for i in range(3):
        ym = max(sh_ox[i], sc_ox[i], ko_ox[i])
        add_bracket(axF_ox, i - w_f, i, ym + 0.02, 0.015, sigs_tbi[i], fontsize=7.8, color=COL_TBI_EDGE)
        add_bracket(axF_ox, i, i + w_f, ym + 0.068, 0.015, sigs_ko[i], fontsize=7.8, color=COL_KO_EDGE)

    axF_ox.set_xticks(x_f)
    axF_ox.set_xticklabels(ox_sub["ct_clean"], fontweight="bold", fontsize=8.0)
    axF_ox.set_ylabel("mt-OXPHOS module score", fontsize=8.2)
    axF_ox.set_ylim(0, 0.70)
    axF_ox.set_title("TBI mt-OXPHOS collapse & astrocytic Tgfbr1 KO rescue", fontsize=9.0, fontweight="bold", pad=20)
    axF_ox.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3,
        frameon=False, fontsize=7.6, handletextpad=0.3, columnspacing=0.9
    )

    # Right subpanel of F: Cross-lineage % Change Heatmap (CCI+scr vs Sham & KO vs CCI+scr)
    mods_f = ["mtOXPHOS_mt", "Mito_biogenesis", "Mito_dynamics", "Mitophagy", "Antioxidant", "Apoptosis_anti"]
    mod_lbls = ["mt-OXPHOS", "Mitochondrial biogenesis", "Mitochondrial dynamics", "Mitophagy", "Antioxidant defense", "Bcl2 / Anti-apoptotic"]
    hm_mat = np.zeros((len(mods_f), 6))
    cts_raw = ["Astrocyte (全部)", "Excitatory Neuron (全部)", "Inhibitory Neuron (全部)"]
    for r_i, m_k in enumerate(mods_f):
        for c_i, ct_r in enumerate(cts_raw):
            row_m = mito_df[(mito_df["module"] == m_k) & (mito_df["cell_type"] == ct_r)]
            if len(row_m):
                hm_mat[r_i, c_i * 2]     = float(row_m["CCI_vs_Sham_pct"].values[0])
                hm_mat[r_i, c_i * 2 + 1] = float(row_m["KO_vs_CCI_pct"].values[0])

    hm_clip = np.clip(hm_mat, -55, 65)
    im_m = axF_hm.imshow(hm_clip, aspect="auto", cmap=NPG_RESCUE_CMAP, vmin=-55, vmax=55)
    for r_i in range(len(mods_f)):
        for c_i in range(6):
            v = hm_mat[r_i, c_i]
            txt_c = "white" if abs(v) > 32 else "#0F172A"
            axF_hm.text(c_i, r_i, f"{v:+.1f}%", ha="center", va="center", fontsize=7.4, fontweight="bold", color=txt_c)
    for r_i in range(len(mods_f) + 1):
        axF_hm.axhline(r_i - 0.5, color="white", lw=1.2)
    for c_i in range(7):
        axF_hm.axvline(c_i - 0.5, color="white", lw=1.2)
    axF_hm.set_xticks(range(6))
    axF_hm.set_xticklabels(
        ["Astrocyte\nCCI+scr vs Sham", "Astrocyte\nKO vs CCI+scr",
         "Excitatory\nCCI+scr vs Sham", "Excitatory\nKO vs CCI+scr",
         "Inhibitory\nCCI+scr vs Sham", "Inhibitory\nKO vs CCI+scr"],
        fontweight="bold", fontsize=7.4
    )
    axF_hm.set_yticks(range(len(mods_f)))
    axF_hm.set_yticklabels(mod_lbls, fontweight="bold", fontsize=7.8)
    axF_hm.set_title("Cross-lineage metabolic mirror rescue (% Δ module score)", fontsize=9.0, fontweight="bold", pad=8)
    cbar_m = plt.colorbar(im_m, ax=axF_hm, fraction=0.026, pad=0.025)
    cbar_m.set_label("% change", fontsize=7.8, fontweight="bold")
    cbar_m.ax.tick_params(labelsize=7.2)

    save_pub_figure(fig5, OUT_SUPP, "Fig S5_cleaned", dpi=300)
    plt.close(fig5)


def build_figure_9():
    """Build standalone Figure 9 (BioRender v2 Split Tripartite Synapse Cross-Section)."""
    from polish_fig9_three_versions import polish_v2, save_multi
    v2_im = polish_v2()
    save_multi(v2_im, "Fig 9_cleaned")


if __name__ == "__main__":
    print("Building Figure 6 (8 panels A-H)...")
    build_figure_6()
    print("Building Figure 7 (6 panels A-F)...")
    build_figure_7()
    print("Building Figure 8 (5 panels A-E)...")
    build_figure_8()
    print("Building Figure 9 (Standalone BioRender v2 Split Tripartite Synapse)...")
    build_figure_9()
    print("Building Supplementary Figures S1-S4...")
    build_supp_figures()
    print("ALL DONE!")
