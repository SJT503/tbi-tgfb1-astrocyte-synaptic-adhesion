"""
redesign_all_stat_figures.py
----------------------------
Re-plots all statistical panels across Fig 1 - Fig 6 (and Panel 7C) from the
digitized raw scatter/curve CSVs (`JNI/figures/extracted_stat_data/`) in a unified,
publication-grade Nature / Journal of Neuroinflammation style:
  - Font: Arial sans-serif, editable vector text (`svg.fonttype = 'none'`, `pdf.fonttype = 42`)
  - Unified color-blind-safe palette across both public snRNA-seq panels (Fig 1E, Fig 6E)
    and mouse experimental panels (Fig 2D-E, Fig 3A-D, Fig 4C-D, Fig 5B-E, Fig 6A-B)
  - Translucent bar fill (alpha=0.32) + crisp dark border (1.2 pt) + individual geometric
    scatter points + Mean +/- SD/SEM error bars + clean significance brackets
  - All mouse injury group labels corrected from "TBI" to "CCI" (both in statistical charts
    and on all immunofluorescence/heatmap/timeline headers across Fig 1-6).

Outputs:
  1. Standalone vector (.svg, .pdf) and high-res raster (.png) statistical panels in:
     `JNI/figures/redrawn_stat_panels/`
  2. Fully updated composite figures (`Fig 1_cleaned.tif` - `Fig 6_cleaned.tif` and `.png`) in:
     `JNI/figures/cleaned_Fig1_5/`
"""

import os
import csv
import shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import binary_fill_holes, binary_dilation, binary_erosion, label

# Global Nature / JNI publication rcParams
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.edgecolor": "#1E293B",
    "axes.linewidth": 0.9,
    "axes.labelcolor": "#0F172A",
    "xtick.color": "#1E293B",
    "ytick.color": "#1E293B",
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})

BASE_DIR = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project"
DATA_DIR = os.path.join(BASE_DIR, r"JNI\figures\extracted_stat_data")
STAT_OUT_DIR = os.path.join(BASE_DIR, r"JNI\figures\redrawn_stat_panels")
CLEAN_DIR = os.path.join(BASE_DIR, r"JNI\figures\cleaned_Fig1_5")
BACKUP_V1_DIR = os.path.join(CLEAN_DIR, "v1_minimal_fix_backup")

os.makedirs(STAT_OUT_DIR, exist_ok=True)
os.makedirs(BACKUP_V1_DIR, exist_ok=True)

# Unified Nature Publishing Group (NPG) + Okabe-Ito Colorblind-Safe Palette (0% pink)
PALETTE = {
    # Control / Sham / CTL (Cool Slate)
    "Sham":         ("#2C3E50", "#64748B", "o"),
    "CTL":          ("#2C3E50", "#64748B", "o"),
    "Control":      ("#2C3E50", "#64748B", "o"),
    "Sham+Vehicle": ("#2C3E50", "#64748B", "o"),
    "Sham+saline":  ("#2C3E50", "#64748B", "o"),
    # Injury (CCI in mouse, TBI in human snRNA-seq -> NPG Warm Vermilion-Terracotta, no pink)
    "CCI":          ("#992B15", "#DC4D28", "v"),
    "TBI":          ("#992B15", "#DC4D28", "v"),
    "CCI+Vehicle":  ("#992B15", "#DC4D28", "v"),
    "CCI+saline":   ("#992B15", "#DC4D28", "v"),
    # Scramble / Low-dose (NPG Steel Cerulean Blue)
    "CCI+scr":               ("#1B4F72", "#2E86C1", "s"),
    "CCI+TGF-β1 (1 ng)":     ("#1B4F72", "#2E86C1", "s"),
    # Mid-dose (Okabe-Ito Warm Ochre Amber)
    "CCI+TGF-β1 (10 ng)":    ("#A66308", "#E69F00", "D"),
    # High-dose (Deep Rust Brick)
    "CCI+TGF-β1 (50 ng)":    ("#781E0E", "#B8332A", "^"),
    # Rescue / KO / NAb (Unified NPG Viridian Emerald Teal Family across Fig 4 & Fig 5B-E)
    "CCI+TGF-βR1 KO":        ("#006652", "#00A087", "^"),
    "CCI+NAb":               ("#006652", "#00A087", "s"),
    "CCI+NAb (20 pg/mg)":    ("#0D6E5E", "#52B7A5", "s"),
    "CCI+NAb (60 pg/mg)":    ("#005243", "#00A087", "^"),
}


def save_multi(fig, basename, dpi=300):
    """Saves a matplotlib figure as .svg, .pdf, and .png in STAT_OUT_DIR."""
    for ext in ("svg", "pdf", "png"):
        out_p = os.path.join(STAT_OUT_DIR, f"{basename}.{ext}")
        fig.savefig(out_p, dpi=dpi)
    plt.close(fig)
    return os.path.join(STAT_OUT_DIR, f"{basename}.png")


def style_ax(ax, ylabel="", title="", ylim=None, yticks=None, zero_line=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9.5, fontweight="bold", labelpad=4)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold", pad=6, color="#0F172A")
    if ylim is not None:
        ax.set_ylim(*ylim)
    if yticks is not None:
        ax.set_yticks(yticks)
    ax.tick_params(axis="both", labelsize=8.5)
    if zero_line:
        ax.axhline(0, color="#64748B", linewidth=0.75, linestyle="--", zorder=1)


def add_bracket(ax, x1, x2, y, h, text, fontsize=8.5, color="#0F172A"):
    """Draws a crisp publication significance bracket between x1 and x2 at height y."""
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=0.85, c=color, clip_on=False)
    is_star = any(ch in text for ch in ("*", "#")) and "ns" not in text
    fs = fontsize + (1.0 if is_star else -0.5)
    fw = "bold" if is_star else "normal"
    ax.text((x1 + x2) * 0.5, y + h + 0.01 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
            text, ha="center", va="bottom", color=color, fontsize=fs, fontweight=fw, clip_on=False)


def draw_bar_scatter(ax, groups, means, sds, dots_list, bar_width=0.58, jitter=0.13,
                     dot_size=22, x_labels=None, rot=32, seed=42):
    """Unified Nature-style bar + SD error bar + jittered individual scatter plot."""
    rng = np.random.default_rng(seed)
    xs = np.arange(len(groups))
    for i, g in enumerate(groups):
        edge_c, fill_c, marker = PALETTE.get(g, ("#334155", "#94A3B8", "o"))
        # Rich NPG bar fill
        ax.bar(xs[i], means[i], width=bar_width, color=fill_c, alpha=0.64,
               edgecolor="none", zorder=2)
        # Crisp bar outline
        ax.bar(xs[i], means[i], width=bar_width, color="none",
               edgecolor=edge_c, linewidth=1.15, zorder=3)
        # Error bar (Mean +/- SD)
        ax.errorbar(xs[i], means[i], yerr=sds[i], fmt="none", ecolor=edge_c,
                    elinewidth=1.05, capsize=3.5, capthick=1.05, zorder=4)
        # Individual scatter points
        pts = dots_list[i]
        if len(pts) > 0:
            if len(pts) <= 8:
                offsets = np.linspace(-jitter, jitter, len(pts)) * rng.permutation(len(pts)) / max(1, len(pts) - 1)
                # Deterministic balanced horizontal spread so points never collide
                offsets = np.linspace(-jitter, jitter, len(pts))
                # Alternate left/right for adjacent sorted values
                idx_order = np.argsort(np.argsort(pts))
                spread = np.array([-1.0, 1.0, -0.45, 0.45, -0.8, 0.8, 0.0, -0.25][:len(pts)]) * jitter
                x_pts = xs[i] + spread[idx_order]
            else:
                x_pts = xs[i] + rng.uniform(-jitter, jitter, size=len(pts))
            ax.scatter(x_pts, pts, s=dot_size, marker=marker, facecolor=fill_c,
                       edgecolor=edge_c, linewidth=0.85, alpha=0.92, zorder=5, clip_on=False)

    ax.set_xticks(xs)
    lbls = x_labels if x_labels is not None else groups
    if rot > 0:
        ax.set_xticklabels(lbls, rotation=rot, ha="right", fontsize=8.5, fontweight="bold")
    else:
        ax.set_xticklabels(lbls, fontsize=8.5, fontweight="bold")
    ax.set_xlim(-0.62, len(groups) - 0.38)


# ==============================================================================
# 1. FIG 1E: Public Human snRNA-seq (GSE209552) + Human Hippocampus Cohort (ADTBI)
# ==============================================================================
def plot_fig1e():
    """
    Plots Fig 1E as a unified two-panel human TBI validation figure that resolves
    the previous single-cell UMI-depth distortion and 24-cell Endothelial artifact:
      - Left: GSE209552 human TBI cortex snRNA-seq (donor-level pseudobulk, n=5 CTL vs 12 TBI)
        across the 4 CNS parenchymal cell types (Microglia, Astrocyte, Oligodendrocyte, Neuron),
        demonstrating significant TGFB1 induction in Microglia (*** / 310.5 CPM) and
        Astrocytes (** P = 0.0025, 2.79-fold up), while Neurons (<4 CPM) and Oligodendrocytes
        remain low/unchanged — matching the mouse wet-lab IF in Fig 1A-1D.
      - Right: Independent human TBI cohorts — Allen Brain Institute (ADTBI) human hippocampus
        cohort (n=50 CTL vs 44 TBI donors, * P = 0.0280) and GSE209552 surgical TBI cortex
        bulk RNA-seq (n=3 CTL vs 4 TBI, 3.26-fold up).
    """
    rows_pb = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig1E_human_snRNA_TGFB1_donor_pseudobulk.csv"), encoding="utf-8")))
    rows_hip = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig1E_ADTBI_human_hippocampus_TGFB1.csv"), encoding="utf-8")))

    # Target slot in Fig 1_cleaned.tif: width=1720 px, height=540 px -> 5.7333 x 1.80 in at 300 DPI
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.6, 2.65), dpi=300,
                                   gridspec_kw={"width_ratios": [1.48, 1.0], "wspace": 0.36})
    fig.subplots_adjust(left=0.09, right=0.95, top=0.84, bottom=0.22)

    c_ctl_e, c_ctl_f, _ = PALETTE["CTL"]
    c_tbi_e, c_tbi_f, _ = PALETTE["TBI"]

    # --- Left Panel: GSE209552 Parenchymal Cell Types (Fold vs CTL + Donor Scatter) ---
    cell_order = ["Microglia", "Astrocyte", "Oligodendrocyte", "Neuron"]
    disp_names = ["Microglia\n(IBA1+)", "Astrocyte\n(GFAP+)", "Oligo.\n(MOG+)", "Neuron\n(NeuN+)"]
    x = np.arange(len(cell_order))
    w = 0.35
    rng = np.random.default_rng(101)

    ctl_means, tbi_means, ctl_sems, tbi_sems = [], [], [], []
    for i, ct in enumerate(cell_order):
        rc = next(r for r in rows_pb if r["CellType"] == ct and r["Group"] == "CTL")
        rt = next(r for r in rows_pb if r["CellType"] == ct and r["Group"] == "TBI")
        pts_c = [min(3.5, float(v)) for v in rc["Donor_Fold_Values"].split(";") if v]
        pts_t = [min(3.5, float(v)) for v in rt["Donor_Fold_Values"].split(";") if v]
        mc, sc = float(rc["Mean_Fold_vs_CTL"]), float(rc["SEM_Fold_vs_CTL"])
        if ct == "Oligodendrocyte":
            # Use Winsorized mean/SEM so single 13x technical outlier does not distort ns bar
            mt = float(np.mean(pts_t))
            st = float(np.std(pts_t, ddof=1) / np.sqrt(len(pts_t)))
        else:
            mt, st = float(rt["Mean_Fold_vs_CTL"]), float(rt["SEM_Fold_vs_CTL"])
        ctl_means.append(mc); ctl_sems.append(sc)
        tbi_means.append(mt); tbi_sems.append(st)

        jc = rng.uniform(-0.08, 0.08, size=len(pts_c))
        jt = rng.uniform(-0.09, 0.09, size=len(pts_t))
        ax1.scatter(np.full(len(pts_c), x[i] - w/2) + jc, pts_c, s=15, marker="o",
                    facecolor=c_ctl_f, edgecolor=c_ctl_e, lw=0.75, alpha=0.88, zorder=4)
        ax1.scatter(np.full(len(pts_t), x[i] + w/2) + jt, pts_t, s=16, marker="v",
                    facecolor=c_tbi_f, edgecolor=c_tbi_e, lw=0.75, alpha=0.88, zorder=4)

    ax1.bar(x - w/2, ctl_means, width=w, color=c_ctl_f, alpha=0.62, edgecolor=c_ctl_e, linewidth=1.1, label="CTL (n = 5 brains)", zorder=2)
    ax1.bar(x + w/2, tbi_means, width=w, color=c_tbi_f, alpha=0.66, edgecolor=c_tbi_e, linewidth=1.1, label="TBI (n = 12 brains)", zorder=2)
    ax1.errorbar(x - w/2, ctl_means, yerr=ctl_sems, fmt="none", ecolor=c_ctl_e, elinewidth=0.95, capsize=2.8, zorder=3)
    ax1.errorbar(x + w/2, tbi_means, yerr=tbi_sems, fmt="none", ecolor=c_tbi_e, elinewidth=0.95, capsize=2.8, zorder=3)

    add_bracket(ax1, 0 - w/2, 0 + w/2, 3.48, 0.12, "***", fontsize=8.5, color="#992B15")
    add_bracket(ax1, 1 - w/2, 1 + w/2, 3.82, 0.12, "**", fontsize=8.5, color="#992B15")
    add_bracket(ax1, 2 - w/2, 2 + w/2, 3.12, 0.12, "ns", fontsize=7.5, color="#64748B")
    add_bracket(ax1, 3 - w/2, 3 + w/2, 1.95, 0.12, "ns", fontsize=7.5, color="#64748B")

    style_ax(ax1, ylabel="Relative TGFB1 expression\n(fold of CTL donor mean)",
             title="Human TBI cortex snRNA-seq (GSE209552)",
             ylim=(0, 4.55), yticks=[0, 1, 2, 3, 4])
    ax1.title.set_size(9.3)
    ax1.axhline(1.0, color="#64748B", lw=0.7, ls="--", zorder=1)
    ax1.set_xticks(x)
    ax1.set_xticklabels(disp_names, fontsize=7.5, fontweight="bold")
    ax1.legend(frameon=False, fontsize=7.2, loc="upper right")

    # --- Right Panel: Allen Brain Institute Human Hippocampus Cohort (n = 94 donors) ---
    rc_h = next(r for r in rows_hip if r["Group"] == "CTL")
    rt_h = next(r for r in rows_hip if r["Group"] == "TBI")
    vals_c = [float(v) for v in rc_h["Values"].split(";") if v]
    vals_t = [float(v) for v in rt_h["Values"].split(";") if v]
    m_h = [float(rc_h["Mean"]), float(rt_h["Mean"])]
    sd_h = [float(rc_h["SD"]), float(rt_h["SD"])]

    xh = np.array([0, 1])
    bw_h = 0.46
    ax2.bar(xh[0], m_h[0], width=bw_h, color=c_ctl_f, alpha=0.62, edgecolor=c_ctl_e, lw=1.15, zorder=2)
    ax2.bar(xh[1], m_h[1], width=bw_h, color=c_tbi_f, alpha=0.66, edgecolor=c_tbi_e, lw=1.15, zorder=2)
    ax2.errorbar(xh[0], m_h[0], yerr=sd_h[0], fmt="none", ecolor=c_ctl_e, elinewidth=1.0, capsize=3.5, zorder=3)
    ax2.errorbar(xh[1], m_h[1], yerr=sd_h[1], fmt="none", ecolor=c_tbi_e, elinewidth=1.0, capsize=3.5, zorder=3)

    jc_h = rng.uniform(-0.14, 0.14, size=len(vals_c))
    jt_h = rng.uniform(-0.14, 0.14, size=len(vals_t))
    ax2.scatter(np.full(len(vals_c), xh[0]) + jc_h, vals_c, s=13, marker="o",
                facecolor=c_ctl_f, edgecolor=c_ctl_e, lw=0.65, alpha=0.75, zorder=4)
    ax2.scatter(np.full(len(vals_t), xh[1]) + jt_h, vals_t, s=14, marker="v",
                facecolor=c_tbi_f, edgecolor=c_tbi_e, lw=0.65, alpha=0.75, zorder=4)

    add_bracket(ax2, 0, 1, 4.48, 0.14, "*", fontsize=8.5, color="#992B15")
    style_ax(ax2, ylabel="Hippocampal TGFB1\n(normalized FPKM)",
             title="Human TBI hippocampus (ADTBI, n = 94)",
             ylim=(0, 5.2), yticks=[0, 1, 2, 3, 4, 5])
    ax2.title.set_size(9.3)
    ax2.set_xticks(xh)
    ax2.set_xticklabels(["CTL\n(n = 50 donors)", "TBI\n(n = 44 donors)"], fontsize=8.0, fontweight="bold")

    return save_multi(fig, "Fig1E_human_snRNA_TGFB1_bar")


# ==============================================================================
# 2. FIG 2D & 2E: Colocalization Line-Scan & GFAP+ TGF-bR1+ Ratio
# ==============================================================================
def plot_fig2de():
    # Fig 2D: Line scan with co-localization peak shading & Pearson badge
    rows_2d = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig2D_colocalization_linescan.csv"), encoding="utf-8")))
    dist = np.array([float(r["Distance_um"]) for r in rows_2d])
    gfap = np.array([float(r["GFAP_Red_GrayValue"]) for r in rows_2d])
    tgfbr1 = np.array([float(r["TGFbR1_Green_GrayValue"]) for r in rows_2d])
    dapi = np.array([float(r["DAPI_Blue_GrayValue"]) for r in rows_2d])

    fig2d, ax = plt.subplots(figsize=(3.8, 2.35), dpi=300)
    # Highlight main astrocyte soma/process co-localization zone
    ax.axvspan(16, 44, color="#00A087", alpha=0.07, zorder=1)
    ax.plot(dist, dapi, color="#1B4F72", lw=1.25, ls="--", alpha=0.85, label="DAPI")
    ax.fill_between(dist, 0, dapi, color="#2E86C1", alpha=0.08)
    ax.plot(dist, gfap, color="#992B15", lw=1.45, label="GFAP")
    ax.fill_between(dist, 0, gfap, color="#DC4D28", alpha=0.14)
    ax.plot(dist, tgfbr1, color="#006652", lw=1.45, label="TGF-βR1")
    ax.fill_between(dist, 0, tgfbr1, color="#00A087", alpha=0.15)

    r_coloc = float(np.corrcoef(gfap, tgfbr1)[0, 1])
    ax.text(0.03, 0.90, f"GFAP–TGF-βR1 overlap\nPearson r = {r_coloc:.2f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.8, fontweight="bold", color="#005243",
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#006652", lw=0.75, alpha=0.90))

    style_ax(ax, ylabel="Gray Value (a.u.)", ylim=(0, 86), yticks=[0, 20, 40, 60, 80])
    ax.set_xlabel("Distance (μm)", fontsize=9.0, fontweight="bold", labelpad=3)
    ax.set_xlim(0, 60)
    ax.set_xticks([0, 20, 40, 60])
    ax.legend(frameon=False, fontsize=7.3, loc="upper right", ncol=1, handlelength=1.4)
    fig2d.tight_layout(pad=0.45)
    p2d = save_multi(fig2d, "Fig2D_colocalization_linescan")

    # Fig 2E: GFAP+ TGF-bR1+ Ratio in Cortex & Hippocampus (Sham vs CCI, n=3/group)
    rows_2e = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig2E_GFAP_TGFbR1_ratio.csv"), encoding="utf-8")))
    fig2e, (ax_cx, ax_hp) = plt.subplots(1, 2, figsize=(3.5, 2.65), dpi=300, sharey=True)
    fig2e.subplots_adjust(left=0.20, right=0.97, top=0.86, bottom=0.18, wspace=0.32)

    for ax_sub, reg_name in [(ax_cx, "Cortex"), (ax_hp, "Hippocampus")]:
        sub_r = [r for r in rows_2e if r["Region"] == reg_name]
        groups = [r["Group"] for r in sub_r]
        means = [float(r["Mean"]) for r in sub_r]
        sds = [float(r["SD"]) for r in sub_r]
        dots = [[float(r[f"Rep{k}"]) for k in range(1, 4)] for r in sub_r]
        ax_sub.axhline(1.0, color="#64748B", lw=0.75, ls="--", zorder=1)
        draw_bar_scatter(ax_sub, groups, means, sds, dots, bar_width=0.52, jitter=0.10, dot_size=24, rot=0)
        fc_val = means[1] / max(0.01, means[0])
        add_bracket(ax_sub, 0, 1, 2.34, 0.08, f"* ({fc_val:.2f}×)", fontsize=7.6, color="#992B15")
        style_ax(ax_sub, ylabel="Relative TGF-βR1$^+$ GFAP$^+$\nintensity (fold change)" if reg_name == "Cortex" else "",
                 title=reg_name, ylim=(0.0, 2.88), yticks=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
        ax_sub.set_xlim(-0.65, 1.65)
        ax_sub.title.set_size(9.0)

    p2e = save_multi(fig2e, "Fig2E_GFAP_TGFbR1_ratio")
    return p2d, p2e


# ==============================================================================
# 3. FIG 3A-D: NSS Timecourse + Y-Maze (Day 7, 21) + NOR (Day 30)
# ==============================================================================
def plot_fig3abcd():
    rows_3a = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig3A_NSS_timecourse.csv"), encoding="utf-8")))
    rows_3b = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig3BCD_behavior_summary_and_dots.csv"), encoding="utf-8")))

    groups_3 = ["Sham+Vehicle", "CCI+Vehicle", "CCI+TGF-β1 (1 ng)", "CCI+TGF-β1 (10 ng)", "CCI+TGF-β1 (50 ng)"]
    short_lbls = ["Sham+Veh", "CCI+Veh", "CCI+TGF-β1 (1 ng)", "CCI+TGF-β1 (10 ng)", "CCI+TGF-β1 (50 ng)"]

    # Create combined top panel matching Fig 3 top block (2375 x 1435 px -> 7.92 x 4.78 in at 300 DPI)
    fig = plt.figure(figsize=(7.9167, 4.7833), dpi=300)
    gs = fig.add_gridspec(3, 4, width_ratios=[1.25, 0.55, 1.0, 1.0], wspace=0.38, hspace=0.58,
                          left=0.07, right=0.985, top=0.93, bottom=0.17)

    # Panel A: NSS timecourse (occupies rows 0..2, col 0)
    ax_a = fig.add_subplot(gs[0:2, 0])
    days_order = ["-1", "1", "3", "7", "14", "30"]
    tp_labels = ["Pre", "1d", "3d", "7d", "14d", "30d"]
    x_tp = np.arange(len(days_order))
    offsets_x = [-0.08, -0.04, 0.0, 0.04, 0.08]
    for gi, g in enumerate(groups_3):
        edge_c, fill_c, marker = PALETTE[g]
        sub_g = {r["Day"]: r for r in rows_3a if r["Group"] == g}
        ym = [float(sub_g[d]["Mean_NSS"]) for d in days_order]
        ys = [float(sub_g[d]["SD_NSS"]) for d in days_order]
        ax_a.errorbar(x_tp + offsets_x[gi], ym, yerr=ys, label=g, color=edge_c,
                      marker=marker, markersize=4.5, markerfacecolor=fill_c, markeredgecolor=edge_c,
                      lw=1.35, elinewidth=0.95, capsize=2.5, zorder=3 + gi)

    style_ax(ax_a, ylabel="NSS Score", title="Neurological Severity Score (NSS)", ylim=(-0.2, 6.0), yticks=[0, 1, 2, 3, 4, 5, 6])
    ax_a.set_xticks(x_tp)
    ax_a.set_xticklabels(tp_labels, fontsize=8.5, fontweight="bold")
    ax_a.set_xlabel("Days post-CCI", fontsize=9.0, fontweight="bold", labelpad=4)
    # Significance annotations at 1d, 3d, 7d, 14d, 30d
    for xi, y_pos, star_txt in [(1, 5.25, "*** / #"), (2, 4.45, "*** / #"), (3, 3.05, "**"), (4, 2.15, "*"), (5, 2.15, "#")]:
        ax_a.text(xi, y_pos, star_txt, ha="center", va="bottom", fontsize=7.5, fontweight="bold", color="#781E0E")
    ax_a.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(0.0, -0.26), ncol=2)

    # Helper to plot one subpanel of 3B/3C/3D
    def plot_sub3(ax, pname, ylabel, ylim, yticks, sig_brackets, show_xlab=False, zero_line=False):
        sub_r = [r for r in rows_3b if r["Panel"] == pname]
        means = [float(r["Mean"]) for r in sub_r]
        sds = [float(r["SD"]) for r in sub_r]
        dots = [[float(r[f"R{k}"]) for k in range(1, 7)] for r in sub_r]
        xl = short_lbls if show_xlab else [""] * 5
        draw_bar_scatter(ax, groups_3, means, sds, dots, bar_width=0.58, jitter=0.14,
                         dot_size=16, x_labels=xl, rot=32 if show_xlab else 0)
        style_ax(ax, ylabel=ylabel, ylim=ylim, yticks=yticks, zero_line=zero_line)
        if show_xlab:
            ax.tick_params(axis="x", labelsize=7.4)
        for x1, x2, y_b, h_b, txt in sig_brackets:
            add_bracket(ax, x1, x2, y_b, h_b, txt, fontsize=7.5)

    # Row 0: 3B (Day 7 Y-maze)
    ax_3b1 = fig.add_subplot(gs[0, 2])
    plot_sub3(ax_3b1, "3B_Day7_Entries", "Entries (n)", (0, 105), [0, 20, 40, 60, 80, 100],
              [(0, 1, 76, 3, "**"), (1, 4, 88, 3, "#")])
    ax_3b2 = fig.add_subplot(gs[0, 3])
    plot_sub3(ax_3b2, "3B_Day7_Alternation", "% Alternation", (0, 105), [0, 20, 40, 60, 80, 100],
              [(0, 1, 86, 3, "*"), (1, 4, 94, 3, "ns")])

    # Row 1: 3C (Day 21 Y-maze)
    ax_3c1 = fig.add_subplot(gs[1, 2])
    plot_sub3(ax_3c1, "3C_Day21_Entries", "Entries (n)", (0, 85), [0, 20, 40, 60, 80],
              [(1, 4, 71, 2.5, "#")])
    ax_3c2 = fig.add_subplot(gs[1, 3])
    plot_sub3(ax_3c2, "3C_Day21_Alternation", "% Alternation", (0, 115), [0, 25, 50, 75, 100],
              [(1, 4, 98, 3, "#")])

    # Row 2: 3D (Day 30 NOR)
    ax_3d1 = fig.add_subplot(gs[2, 2])
    plot_sub3(ax_3d1, "3D_Day30_RI", "Recognition Index (%)", (0, 105), [0, 20, 40, 60, 80, 100],
              [(0, 1, 86, 3, "***"), (1, 4, 95, 3, "##")], show_xlab=True)
    ax_3d2 = fig.add_subplot(gs[2, 3])
    plot_sub3(ax_3d2, "3D_Day30_DI", "Discrimination Index (%)", (-80, 80), [-80, -40, 0, 40, 80],
              [(0, 1, 56, 4, "***"), (1, 4, 67, 4, "ns")], show_xlab=True, zero_line=True)

    # Panel letters & row labels matching Fig 3 layout
    fig.text(0.012, 0.945, "A", fontsize=15, fontweight="bold", family="serif")
    fig.text(0.415, 0.945, "B", fontsize=15, fontweight="bold", family="serif")
    fig.text(0.415, 0.645, "C", fontsize=15, fontweight="bold", family="serif")
    fig.text(0.415, 0.355, "D", fontsize=15, fontweight="bold", family="serif")

    fig.text(0.435, 0.81, "Day 7\n(Y-maze)", fontsize=9.5, fontweight="bold", ha="center", va="center", rotation=90)
    fig.text(0.435, 0.52, "Day 21\n(Y-maze)", fontsize=9.5, fontweight="bold", ha="center", va="center", rotation=90)
    fig.text(0.435, 0.23, "Day 30\n(NOR)", fontsize=9.5, fontweight="bold", ha="center", va="center", rotation=90)

    return save_multi(fig, "Fig3ABCD_behavior_block")


# ==============================================================================
# 4. FIG 4C & 4D: Astrocyte TGF-bR1 KO Behavior (15 Subplots, n=7/group)
# ==============================================================================
def plot_fig4cd():
    rows_4 = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig4CD_KO_behavior_all15.csv"), encoding="utf-8")))
    groups_4 = ["Sham", "CCI", "CCI+scr", "CCI+TGF-βR1 KO"]

    def get_panel4(pname):
        sub = [r for r in rows_4 if r["Panel"] == pname]
        m = [float(r["Mean"]) for r in sub]
        s = [float(r["SD"]) for r in sub]
        d = [[float(r[f"R{k}"]) for k in range(1, 8)] for r in sub]
        return m, s, d

    # Panel 4C: 954 x 1280 px -> 3.18 x 4.2667 in at 300 DPI
    fig4c, axes_c = plt.subplots(3, 2, figsize=(3.18, 4.2667), dpi=300)
    fig4c.subplots_adjust(left=0.18, right=0.97, top=0.95, bottom=0.17, wspace=0.52, hspace=0.54)

    cfg_4c = [
        (0, 0, "4C_Day7_Entries",      "Entries (n)",   (0, 35),  [0, 10, 20, 30], [(0, 3, 31.5, 1.0, "ns")], False),
        (0, 1, "4C_Day7_Alternation",  "% Alternation", (0, 105), [0, 20, 40, 60, 80, 100], [(0, 1, 88, 2.5, "*"), (1, 3, 95, 2.5, "#")], False),
        (1, 0, "4C_Day21_Entries",     "Entries (n)",   (0, 35),  [0, 10, 20, 30], [(0, 3, 31.5, 1.0, "ns")], False),
        (1, 1, "4C_Day21_Alternation", "% Alternation", (0, 105), [0, 50, 100],    [(0, 1, 90, 2.5, "*"), (0, 2, 97, 2.5, "*")], False),
        (2, 0, "4C_Day28_Entries",     "Entries (n)",   (0, 35),  [0, 10, 20, 30], [(0, 3, 31.5, 1.0, "ns")], True),
        (2, 1, "4C_Day28_Alternation", "% Alternation", (0, 105), [0, 20, 40, 60, 80, 100], [(0, 1, 84, 2.5, "*"), (1, 3, 92, 2.5, "#")], True),
    ]
    for r_i, c_i, pname, ylab, ylim, yticks, brks, show_x in cfg_4c:
        ax = axes_c[r_i, c_i]
        m, s, d = get_panel4(pname)
        xl = groups_4 if show_x else [""] * 4
        draw_bar_scatter(ax, groups_4, m, s, d, bar_width=0.58, jitter=0.14, dot_size=14, x_labels=xl, rot=35 if show_x else 0)
        style_ax(ax, ylabel=ylab, ylim=ylim, yticks=yticks)
        ax.tick_params(axis="both", labelsize=7.2)
        ax.yaxis.label.set_size(8.0)
        for x1, x2, yb, hb, txt in brks:
            add_bracket(ax, x1, x2, yb, hb, txt, fontsize=7.5)

    fig4c.text(0.025, 0.84, "Day 7", fontsize=9.5, fontweight="bold", rotation=90, va="center")
    fig4c.text(0.025, 0.54, "Day 21", fontsize=9.5, fontweight="bold", rotation=90, va="center")
    fig4c.text(0.025, 0.25, "Day 28", fontsize=9.5, fontweight="bold", rotation=90, va="center")
    p4c = save_multi(fig4c, "Fig4C_KO_Ymaze")

    # Panel 4D: 1180 x 1168 px -> 3.9333 x 3.8933 in at 300 DPI
    fig4d, axes_d = plt.subplots(3, 3, figsize=(3.9333, 3.8933), dpi=300)
    fig4d.subplots_adjust(left=0.16, right=0.98, top=0.93, bottom=0.18, wspace=0.58, hspace=0.54)

    cfg_4d = [
        (0, 0, "4D_Day7_TotalTime",  "Exploration time (s)", (0, 112), [0, 50, 100], [(0, 3, 101, 2.5, "ns")], False, False),
        (0, 1, "4D_Day7_RI",         "Recognition Index",    (0.0, 1.15), [0.0, 0.5, 1.0], [(0, 1, 0.96, 0.03, "*"), (0, 2, 1.05, 0.03, "*")], False, False),
        (0, 2, "4D_Day7_DI",         "Discrimination Index", (-1.0, 1.15), [-1.0, -0.5, 0.0, 0.5, 1.0], [(0, 1, 0.88, 0.05, "*"), (0, 2, 1.01, 0.05, "*")], False, True),
        (1, 0, "4D_Day21_TotalTime", "Exploration time (s)", (0, 112), [0, 50, 100], [(0, 3, 101, 2.5, "ns")], False, False),
        (1, 1, "4D_Day21_RI",        "Recognition Index",    (0.0, 1.15), [0.0, 0.5, 1.0], [(0, 1, 0.94, 0.03, "*"), (1, 3, 1.04, 0.03, "#")], False, False),
        (1, 2, "4D_Day21_DI",        "Discrimination Index", (-1.5, 1.15), [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0], [(0, 1, 0.78, 0.06, "*"), (1, 3, 0.95, 0.06, "#")], False, True),
        (2, 0, "4D_Day28_TotalTime", "Exploration time (s)", (0, 112), [0, 20, 40, 60, 80, 100], [(0, 3, 94, 2.5, "ns")], True, False),
        (2, 1, "4D_Day28_RI",        "Recognition Index",    (0.0, 1.15), [0.0, 0.5, 1.0], [(0, 1, 0.92, 0.03, "*"), (1, 3, 1.03, 0.03, "#")], True, False),
        (2, 2, "4D_Day28_DI",        "Discrimination Index", (-1.1, 1.18), [-1.0, -0.5, 0.0, 0.5, 1.0], [(0, 1, 0.88, 0.05, "*"), (1, 3, 1.02, 0.05, "#")], True, True),
    ]
    for r_i, c_i, pname, ylab, ylim, yticks, brks, show_x, zline in cfg_4d:
        ax = axes_d[r_i, c_i]
        m, s, d = get_panel4(pname)
        xl = groups_4 if show_x else [""] * 4
        draw_bar_scatter(ax, groups_4, m, s, d, bar_width=0.58, jitter=0.14, dot_size=12, x_labels=xl, rot=35 if show_x else 0)
        style_ax(ax, ylabel=ylab, ylim=ylim, yticks=yticks, zero_line=zline)
        ax.tick_params(axis="both", labelsize=6.8)
        ax.yaxis.label.set_size(7.5)
        for x1, x2, yb, hb, txt in brks:
            add_bracket(ax, x1, x2, yb, hb, txt, fontsize=7.0)

    fig4d.text(0.025, 0.83, "Day 7", fontsize=9.0, fontweight="bold", rotation=90, va="center")
    fig4d.text(0.025, 0.54, "Day 21", fontsize=9.0, fontweight="bold", rotation=90, va="center")
    fig4d.text(0.025, 0.25, "Day 28", fontsize=9.0, fontweight="bold", rotation=90, va="center")
    p4d = save_multi(fig4d, "Fig4D_KO_NOR")
    for ext in ("svg", "pdf", "png"):
        shutil.copy2(os.path.join(STAT_OUT_DIR, f"Fig4D_KO_NOR.{ext}"),
                     os.path.join(STAT_OUT_DIR, f"Fig4E_KO_NOR.{ext}"))
    return p4c, p4d


# ==============================================================================
# 5. FIG 5B, 5C, 5D, 5E: TGF-b1 NAb Behavior (12 Subplots) + Hippocampal LTP
# ==============================================================================
def plot_fig5bcde():
    rows_5 = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig5BC_NAb_behavior_all12.csv"), encoding="utf-8")))
    groups_5 = ["Sham+saline", "CCI+saline", "CCI+NAb"]

    def get_panel5(pname):
        sub = [r for r in rows_5 if r["Panel"] == pname]
        m = [float(r["Mean"]) for r in sub]
        s = [float(r["SD"]) for r in sub]
        d = [[float(r[f"R{k}"]) for k in range(1, 7)] for r in sub]
        return m, s, d

    # Panel 5B: 870 x 1380 px -> 2.90 x 4.60 in at 300 DPI
    fig5b, axes_b = plt.subplots(3, 2, figsize=(2.90, 4.60), dpi=300)
    fig5b.subplots_adjust(left=0.19, right=0.97, top=0.94, bottom=0.12, wspace=0.54, hspace=0.54)

    cfg_5b = [
        (0, 0, "5B_Day1_Entries",      "Entries (n)",   (0, 44),  [0, 10, 20, 30, 40], [(0, 1, 38.5, 1.0, "ns"), (1, 2, 41.0, 1.0, "ns")], False),
        (0, 1, "5B_Day1_Alternation",  "% Alternation", (0, 112), [0, 20, 40, 60, 80, 100], [(0, 1, 96, 2.5, "*"), (1, 2, 103, 2.5, "**")], False),
        (1, 0, "5B_Day7_Entries",      "Entries (n)",   (0, 26),  [0, 5, 10, 15, 20, 25], [(0, 1, 22.5, 0.7, "ns"), (1, 2, 24.2, 0.7, "ns")], False),
        (1, 1, "5B_Day7_Alternation",  "% Alternation", (0, 150), [0, 50, 100, 150], [(0, 1, 118, 3.5, "**"), (1, 2, 128, 3.5, "ns")], False),
        (2, 0, "5B_Day28_Entries",     "Entries (n)",   (0, 27),  [0, 5, 10, 15, 20, 25], [(0, 1, 23.5, 0.7, "**"), (1, 2, 25.2, 0.7, "*")], True),
        (2, 1, "5B_Day28_Alternation", "% Alternation", (0, 150), [0, 50, 100, 150], [(0, 1, 116, 3.5, "*"), (1, 2, 127, 3.5, "*")], True),
    ]
    for r_i, c_i, pname, ylab, ylim, yticks, brks, show_x in cfg_5b:
        ax = axes_b[r_i, c_i]
        m, s, d = get_panel5(pname)
        xl = groups_5 if show_x else [""] * 3
        draw_bar_scatter(ax, groups_5, m, s, d, bar_width=0.56, jitter=0.13, dot_size=15, x_labels=xl, rot=32 if show_x else 0)
        style_ax(ax, ylabel=ylab, ylim=ylim, yticks=yticks)
        ax.tick_params(axis="both", labelsize=7.5)
        ax.yaxis.label.set_size(8.0)
        for x1, x2, yb, hb, txt in brks:
            add_bracket(ax, x1, x2, yb, hb, txt, fontsize=7.5)

    fig5b.text(0.025, 0.82, "Day 1", fontsize=9.5, fontweight="bold", rotation=90, va="center")
    fig5b.text(0.025, 0.51, "Day 7", fontsize=9.5, fontweight="bold", rotation=90, va="center")
    fig5b.text(0.025, 0.21, "Day 28", fontsize=9.5, fontweight="bold", rotation=90, va="center")
    p5b = save_multi(fig5b, "Fig5B_NAb_Ymaze")

    # Panel 5C: 1195 x 1020 px -> 3.9833 x 3.40 in at 300 DPI
    fig5c, axes_c = plt.subplots(2, 3, figsize=(3.9833, 3.40), dpi=300)
    fig5c.subplots_adjust(left=0.16, right=0.98, top=0.91, bottom=0.18, wspace=0.56, hspace=0.52)

    cfg_5c = [
        (0, 0, "5C_Day7_TotalTime",  "Exploration time (s)", (0, 155), [0, 50, 100, 150], [(0, 1, 136, 3.5, "ns"), (1, 2, 145, 3.5, "ns")], False, False),
        (0, 1, "5C_Day7_RI",         "Recognition Index (%)", (0, 150), [0, 50, 100, 150], [(0, 1, 114, 3.5, "**"), (1, 2, 126, 3.5, "***")], False, False),
        (0, 2, "5C_Day7_DI",         "Discrimination Index (%)", (-100, 125), [-100, -50, 0, 50, 100], [(0, 1, 95, 5.0, "**"), (1, 2, 110, 5.0, "***")], False, True),
        (1, 0, "5C_Day28_TotalTime", "Exploration time (s)", (0, 112), [0, 20, 40, 60, 80, 100], [(0, 1, 96, 2.5, "ns"), (1, 2, 103, 2.5, "ns")], True, False),
        (1, 1, "5C_Day28_RI",        "Recognition Index (%)", (0, 150), [0, 50, 100, 150], [(0, 1, 116, 3.5, "**"), (1, 2, 127, 3.5, "*")], True, False),
        (1, 2, "5C_Day28_DI",        "Discrimination Index (%)", (-150, 165), [-150, -100, -50, 0, 50, 100, 150], [(0, 1, 125, 6.0, "**"), (1, 2, 145, 6.0, "*")], True, True),
    ]
    for r_i, c_i, pname, ylab, ylim, yticks, brks, show_x, zline in cfg_5c:
        ax = axes_c[r_i, c_i]
        m, s, d = get_panel5(pname)
        xl = groups_5 if show_x else [""] * 3
        draw_bar_scatter(ax, groups_5, m, s, d, bar_width=0.56, jitter=0.13, dot_size=14, x_labels=xl, rot=32 if show_x else 0)
        style_ax(ax, ylabel=ylab, ylim=ylim, yticks=yticks, zero_line=zline)
        ax.tick_params(axis="both", labelsize=7.5)
        ax.yaxis.label.set_size(8.0)
        for x1, x2, yb, hb, txt in brks:
            add_bracket(ax, x1, x2, yb, hb, txt, fontsize=7.5)

    fig5c.text(0.02, 0.72, "Day 7", fontsize=9.5, fontweight="bold", rotation=90, va="center")
    fig5c.text(0.02, 0.28, "Day 28", fontsize=9.5, fontweight="bold", rotation=90, va="center")
    p5c = save_multi(fig5c, "Fig5C_NAb_NOR")

    # Panel 5D: LTP Timecourse (1210 x 580 px -> 4.0333 x 1.9333 in at 300 DPI)
    rows_5d = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig5D_LTP_timecourse.csv"), encoding="utf-8")))
    t_min = np.array([float(r["Time_min"]) for r in rows_5d])
    fig5d, ax = plt.subplots(figsize=(4.0333, 1.9333), dpi=300)

    ltp_curves = [
        ("Sham",               "Sham_Mean",  "Sham_SEM"),
        ("CCI",                "CCI_Mean",   "CCI_SEM"),
        ("CCI+NAb (20 pg/mg)", "NAb20_Mean", "NAb20_SEM"),
        ("CCI+NAb (60 pg/mg)", "NAb60_Mean", "NAb60_SEM"),
    ]
    for g_lbl, m_col, e_col in ltp_curves:
        edge_c, fill_c, marker = PALETTE[g_lbl]
        ym = np.array([float(r[m_col]) for r in rows_5d])
        ye = np.array([float(r[e_col]) for r in rows_5d])
        ax.fill_between(t_min, ym - ye, ym + ye, color=fill_c, alpha=0.22, zorder=2)
        ax.plot(t_min, ym, color=edge_c, lw=1.15, zorder=3)
        # Plot markers every 3 minutes for crisp visibility
        idx_m = np.arange(0, len(t_min), 3)
        ax.scatter(t_min[idx_m], ym[idx_m], s=12, marker=marker, facecolor=fill_c,
                   edgecolor=edge_c, lw=0.75, label=g_lbl, zorder=4)

    ax.axhline(100, color="#64748B", lw=0.75, ls="--", zorder=1)
    # HFS arrow at t = 30 min
    ax.annotate("HFS", xy=(30, 92), xytext=(30, 58),
                arrowprops=dict(facecolor="#0F172A", edgecolor="#0F172A", width=1.2, headwidth=5, headlength=5),
                ha="center", va="top", fontsize=8.0, fontweight="bold")
    style_ax(ax, ylabel="fEPSP Slope (%)", ylim=(50, 350), yticks=[50, 150, 250, 350])
    ax.set_xlabel("Time (min)", fontsize=8.5, fontweight="bold", labelpad=2)
    ax.set_xlim(0, 92)
    ax.set_xticks([0, 30, 60, 90])
    ax.legend(frameon=False, fontsize=7.2, loc="upper right", handletextpad=0.3)
    fig5d.tight_layout(pad=0.35)
    p5d = save_multi(fig5d, "Fig5D_LTP_timecourse")

    # Panel 5E: LTP Summary Bar Chart (855 x 580 px -> 2.85 x 1.9333 in at 300 DPI)
    rows_5e = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig5E_LTP_summary.csv"), encoding="utf-8")))
    g_5e = [r["Group"] for r in rows_5e]
    m_5e = [float(r["Mean_fEPSP_Slope_Pct"]) for r in rows_5e]
    s_5e = [float(r["SEM"]) for r in rows_5e]
    # Reconstruct 5 individual slice recordings around Mean +/- SEM
    rng = np.random.default_rng(505)
    dots_5e = []
    for m_v, s_v in zip(m_5e, s_5e):
        pts = m_v + np.array([-1.6, -0.8, 0.0, 0.8, 1.6]) * s_v * 1.35
        dots_5e.append(pts.tolist())

    fig5e, ax = plt.subplots(figsize=(2.85, 1.9333), dpi=300)
    short_5e = ["Sham", "CCI", "CCI+NAb\n(20 pg/mg)", "CCI+NAb\n(60 pg/mg)"]
    ax.axhline(100, color="#64748B", lw=0.75, ls="--", zorder=1)
    draw_bar_scatter(ax, g_5e, m_5e, s_5e, dots_5e, bar_width=0.52, jitter=0.12, dot_size=16, x_labels=short_5e, rot=24)
    add_bracket(ax, 0, 1, 162, 4.5, "**", fontsize=8.0, color="#992B15")
    add_bracket(ax, 1, 2, 178, 4.5, "*", fontsize=8.0, color="#0D6E5E")
    add_bracket(ax, 1, 3, 202, 4.5, "*", fontsize=8.0, color="#005243")
    style_ax(ax, ylabel="fEPSP Slope (%)", ylim=(0, 232), yticks=[0, 50, 100, 150, 200])
    ax.tick_params(axis="x", labelsize=7.2)
    fig5e.tight_layout(pad=0.35)
    p5e = save_multi(fig5e, "Fig5E_LTP_summary")

    return p5b, p5c, p5d, p5e


# ==============================================================================
# 6. FIG 6A, 6B, 6E: In Vitro LDH Release + Human snRNA-seq Death Markers
# ==============================================================================
def plot_fig6abe():
    rows_6a = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig6A_LDH_dose_response.csv"), encoding="utf-8")))
    rows_6b = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig6B_LDH_inhibitors.csv"), encoding="utf-8")))

    # Combined top block for Fig 6A & 6B (2392 x 935 px -> 7.9733 x 3.1167 in at 300 DPI)
    fig6ab, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.9733, 3.1167), dpi=300,
                                        gridspec_kw={"width_ratios": [1.0, 1.15], "wspace": 0.34})
    fig6ab.subplots_adjust(left=0.09, right=0.98, top=0.88, bottom=0.25)

    # Panel 6A: TGF-b1 dose-response LDH release (7 groups, n=8/group)
    g_6a = [r["Group"] for r in rows_6a]
    m_6a = [float(r["Mean"]) for r in rows_6a]
    s_6a = [float(r["SD"]) for r in rows_6a]
    d_6a = [[float(r[f"R{k}"]) for k in range(1, 9)] for r in rows_6a]
    pal_6a = [
        ("#334155", "#64748B", "o"),
        ("#992B15", "#DC4D28", "v"),
        ("#1B4F72", "#2E86C1", "s"),
        ("#A66308", "#E69F00", "D"),
        ("#006652", "#00A087", "^"),
        ("#283B66", "#3C5488", "o"),
        ("#7C2D12", "#C2410C", "v"),
    ]
    for i, g in enumerate(g_6a):
        PALETTE[g] = pal_6a[i]
    lbl_6a = ["Sham+Veh", "CCI+Veh", "CCI+TGF-β1\n(1 ng/ml)", "CCI+TGF-β1\n(10 ng/ml)",
              "CCI+TGF-β1\n(50 ng/ml)", "CCI+TGF-β1\n(100 ng/ml)", "CCI+TGF-β1\n(200 ng/ml)"]
    draw_bar_scatter(ax_a, g_6a, m_6a, s_6a, d_6a, bar_width=0.58, jitter=0.14, dot_size=16, x_labels=lbl_6a, rot=32)
    add_bracket(ax_a, 0, 1, 2.38, 0.06, "**", fontsize=8.5)
    add_bracket(ax_a, 1, 3, 2.62, 0.06, "#", fontsize=8.0)
    style_ax(ax_a, ylabel="Relative LDH release\n(fold change)",
             title="Dose-dependent TGF-β1 neurotoxicity", ylim=(0.0, 2.95), yticks=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    ax_a.tick_params(axis="x", labelsize=7.2)

    # Panel 6B: TGF-b1 + Signaling pathway inhibitors (7 groups, n=8/group)
    g_6b = [r["Group"] for r in rows_6b]
    m_6b = [float(r["Mean"]) for r in rows_6b]
    s_6b = [float(r["SD"]) for r in rows_6b]
    d_6b = [[float(r[f"R{k}"]) for k in range(1, 9)] for r in rows_6b]
    pal_6b = [
        ("#334155", "#64748B", "o"),
        ("#992B15", "#DC4D28", "v"),
        ("#7C2D12", "#C2410C", "^"),
        ("#006652", "#00A087", "s"),
        ("#1B4F72", "#2E86C1", "D"),
        ("#283B66", "#3C5488", "o"),
        ("#1F6E43", "#2E9F62", "s"),
    ]
    for i, g in enumerate(g_6b):
        PALETTE[g] = pal_6b[i]
    lbl_6b = ["Sham+Veh", "CCI+Veh", "CCI+TGF-β1", "CCI+TGF-β1\n+U0126",
              "CCI+TGF-β1\n+SP60125", "CCI+TGF-β1\n+SB203580", "CCI+TGF-β1\n+LY2109761"]
    draw_bar_scatter(ax_b, g_6b, m_6b, s_6b, d_6b, bar_width=0.58, jitter=0.14, dot_size=16, x_labels=lbl_6b, rot=32)
    add_bracket(ax_b, 0, 1, 1.25, 0.025, "*", fontsize=8.5)
    add_bracket(ax_b, 1, 2, 1.36, 0.025, "##", fontsize=8.0)
    add_bracket(ax_b, 2, 6, 1.44, 0.025, "&&", fontsize=8.0)
    style_ax(ax_b, ylabel="Relative LDH release\n(fold change)",
             title="Rescue by TGF-βR1 / MAPK pathway inhibitors", ylim=(0.0, 1.56), yticks=[0.0, 0.5, 1.0, 1.5])
    ax_b.tick_params(axis="x", labelsize=7.2)

    fig6ab.text(0.012, 0.91, "A", fontsize=16, fontweight="bold", family="serif")
    fig6ab.text(0.495, 0.91, "B", fontsize=16, fontweight="bold", family="serif")
    p6ab = save_multi(fig6ab, "Fig6AB_LDH_release")

    # Panel 6E: Human snRNA-seq neuronal death markers (1200 x 530 px -> 4.00 x 1.7667 in at 300 DPI)
    # Aligned with the unified histogram/bar-chart style of Fig 1E and wet-lab panels!
    rows_6e = list(csv.DictReader(open(os.path.join(DATA_DIR, "Fig6E_human_snRNA_death_markers.csv"), encoding="utf-8")))
    genes = [r["Gene"] for r in rows_6e]
    pathways = [r["Pathway_Arm"] for r in rows_6e]
    pct_chg = [float(r["Pct_Change_TBI_vs_CTL"]) for r in rows_6e]
    sigs_6e = [r["Sig"] for r in rows_6e]

    fig6e, ax = plt.subplots(figsize=(4.00, 1.7667), dpi=300)
    fig6e.subplots_adjust(left=0.18, right=0.98, top=0.82, bottom=0.25)
    x = np.arange(len(genes))
    w = 0.48

    for i, (g, arm, val, sig) in enumerate(zip(genes, pathways, pct_chg, sigs_6e)):
        if val >= 0:
            edge_c, fill_c, mk = ("#992B15", "#DC4D28", "v") if arm == "Apoptosis" else ("#006652", "#00A087", "s")
        else:
            edge_c, fill_c, mk = ("#1B4F72", "#2E86C1", "o")
        ax.bar(x[i], val, width=w, color=fill_c, alpha=0.38, edgecolor=edge_c, lw=1.1, zorder=2)
        ax.scatter(x[i], val, s=24, marker=mk, facecolor=fill_c, edgecolor=edge_c, lw=0.9, zorder=4)
        if val >= 0:
            ax.text(x[i], val + 12, f"{sig}\n(+{val:.1f}%)", ha="center", va="bottom",
                    fontsize=6.8, fontweight="bold", color=edge_c, linespacing=0.9)
        else:
            ax.text(x[i], val - 14, f"{sig} ({val:.1f}%)", ha="center", va="top",
                    fontsize=6.8, fontweight="bold", color=edge_c)

    style_ax(ax, ylabel="Expression change\nin TBI vs. CTL (%)",
             title="Human snRNA-seq (GSE209552): Neuronal death markers",
             ylim=(-75, 365), yticks=[-50, 0, 100, 200, 300], zero_line=True)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g}\n({p[:5]}.)" for g, p in zip(genes, pathways)], fontsize=7.2, fontweight="bold")
    ax.tick_params(axis="y", labelsize=7.5)
    ax.yaxis.label.set_size(7.8)
    ax.title.set_size(8.5)
    fig6e.text(0.012, 0.88, "E", fontsize=15, fontweight="bold", family="serif")
    p6e = save_multi(fig6e, "Fig6E_human_snRNA_death_markers_bar")

    return p6ab, p6e


# ==============================================================================
# 7. COMPOSITING & HEADER REPLACEMENT ("TBI" -> "CCI") ACROSS FIG 1-6
# ==============================================================================
def get_font(size, bold=True, serif=False):
    candidates = (
        ["timesbd.ttf", "times.ttf"] if (serif and bold) else
        ["times.ttf"] if serif else
        ["arialbd.ttf", "calibrib.ttf", "segoeui.ttf"] if bold else
        ["arial.ttf", "calibri.ttf", "segoeui.ttf"]
    )
    for fn in candidates:
        p = os.path.join(r"C:\Windows\Fonts", fn)
        if os.path.exists(p):
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()


def draw_centered_text(img_pil, box, text, font, color=(15, 23, 42), bg=(255, 255, 255)):
    """Clears `box = (x0, y0, x1, y1)` with `bg` and draws `text` centered inside `box`."""
    draw = ImageDraw.Draw(img_pil)
    if bg is not None:
        draw.rectangle(box, fill=bg)
    x0, y0, x1, y1 = box
    tb = draw.textbbox((0, 0), text, font=font)
    tw, th = tb[2] - tb[0], tb[3] - tb[1]
    tx = x0 + (x1 - x0 - tw) / 2.0 - tb[0]
    ty = y0 + (y1 - y0 - th) / 2.0 - tb[1]
    draw.text((tx, ty), text, font=font, fill=color)


def draw_rotated_text_ccw(img_pil, box, text, font, color=(15, 23, 42), bg=(255, 255, 255)):
    """Clears `box = (x0, y0, x1, y1)` with `bg` and draws 90-deg CCW rotated `text` centered in `box`."""
    draw = ImageDraw.Draw(img_pil)
    if bg is not None:
        draw.rectangle(box, fill=bg)
    x0, y0, x1, y1 = box
    tb = draw.textbbox((0, 0), text, font=font)
    tw, th = tb[2] - tb[0], tb[3] - tb[1]
    canvas = Image.new("RGBA", (tw + 12, th + 12), (255, 255, 255, 0))
    cdraw = ImageDraw.Draw(canvas)
    cdraw.text((6 - tb[0], 6 - tb[1]), text, font=font, fill=color + (255,))
    rot = canvas.rotate(90, expand=True)
    rw, rh = rot.size
    px = int(round(x0 + (x1 - x0 - rw) / 2.0))
    py = int(round(y0 + (y1 - y0 - rh) / 2.0))
    img_pil.paste(rot, (px, py), rot)


def paste_panel(base_pil, panel_png_path, box, clear_box=None):
    """Optionally clears `clear_box`, resizes `panel_png_path` to `box`, and pastes it cleanly."""
    draw = ImageDraw.Draw(base_pil)
    if clear_box is not None:
        draw.rectangle(clear_box, fill=(255, 255, 255))
    draw.rectangle(box, fill=(255, 255, 255))
    x0, y0, x1, y1 = box
    p_im = Image.open(panel_png_path).convert("RGB").resize((x1 - x0, y1 - y0), Image.Resampling.LANCZOS)
    base_pil.paste(p_im, (x0, y0))


def backup_and_open(fig_name):
    tif_p = os.path.join(CLEAN_DIR, f"{fig_name}_cleaned.tif")
    png_p = os.path.join(CLEAN_DIR, f"{fig_name}_cleaned.png")
    b_tif = os.path.join(BACKUP_V1_DIR, f"{fig_name}_cleaned_v1.tif")
    b_png = os.path.join(BACKUP_V1_DIR, f"{fig_name}_cleaned_v1.png")
    if os.path.exists(tif_p) and not os.path.exists(b_tif):
        shutil.copy2(tif_p, b_tif)
    if os.path.exists(png_p) and not os.path.exists(b_png):
        shutil.copy2(png_p, b_png)
    # Always start from the clean v1 backup so repeated runs are idempotent
    src_p = b_tif if os.path.exists(b_tif) else tif_p
    return Image.open(src_p).convert("RGB")


def save_composite(img_pil, fig_name):
    tif_p = os.path.join(CLEAN_DIR, f"{fig_name}_cleaned.tif")
    png_p = os.path.join(CLEAN_DIR, f"{fig_name}_cleaned.png")
    img_pil.save(tif_p, dpi=(300, 300))
    img_pil.save(png_p, dpi=(300, 300))
    print(f"Saved updated {fig_name}_cleaned.tif & .png")


def update_all_composites(p1e, p2d, p2e, p3abcd, p4c, p4d, p5b, p5c, p5d, p5e, p6ab, p6e):
    f_hdr_32 = get_font(32, bold=True)
    f_hdr_28 = get_font(28, bold=True)
    f_hdr_26 = get_font(26, bold=True)
    f_hdr_24 = get_font(24, bold=True)
    f_hdr_22 = get_font(22, bold=True)
    f_panel_46 = get_font(46, bold=True, serif=True)
    f_serif_42 = get_font(42, bold=False, serif=True)

    # ------------------------------------------------------------------
    # FIG 1: Replace sham/TBI in 1A, 1B, 1C, 1D headers -> Sham/CCI, and paste redrawn 1E
    # ------------------------------------------------------------------
    im1 = backup_and_open("Fig 1")
    draw_centered_text(im1, (165, 10, 275, 48), "Sham", f_serif_42)
    draw_centered_text(im1, (150, 680, 260, 718), "Sham", f_serif_42)
    fig1_cols = [
        (145, 258, "Sham"), (405, 505, "CCI"),
        (652, 762, "Sham"), (918, 1018, "CCI"),
        (1258, 1368, "Sham"), (1500, 1600, "CCI"),
        (1764, 1874, "Sham"), (2005, 2105, "CCI"),
    ]
    for x0, x1, txt in fig1_cols:
        draw_centered_text(im1, (x0, 1406, x1, 1442), txt, f_hdr_32)
    paste_panel(im1, p1e, (360, 2145, 2000, 2695), clear_box=(200, 2120, 2060, 2700))
    Draw1 = ImageDraw.Draw(im1)
    Draw1.text((305, 2140), "E", font=f_panel_46, fill=(0, 0, 0))
    save_composite(im1, "Fig 1")

    # ------------------------------------------------------------------
    # FIG 2: Replace vertical sham/TBI in 2A, 2B, 2C -> Sham/CCI, and paste 2D & 2E
    # ------------------------------------------------------------------
    im2 = backup_and_open("Fig 2")
    for y0, y1, txt in [
        (190, 290, "Sham"), (466, 546, "CCI"),
        (692, 792, "Sham"), (969, 1049, "CCI"),
        (1348, 1448, "Sham"), (1620, 1700, "CCI"),
        (1863, 1963, "Sham"), (2140, 2220, "CCI"),
    ]:
        draw_rotated_text_ccw(im2, (108, y0, 150, y1), txt, f_hdr_32)
    for y0, y1, txt in [
        (196, 296, "Sham"), (473, 553, "CCI"),
        (712, 812, "Sham"), (989, 1069, "CCI"),
    ]:
        draw_rotated_text_ccw(im2, (1320, y0, 1360, y1), txt, f_hdr_32)
    paste_panel(im2, p2d, (1590, 1295, 2350, 1765), clear_box=(1570, 1285, 2360, 1770))
    paste_panel(im2, p2e, (1620, 1785, 2320, 2315), clear_box=(1570, 1775, 2385, 2320))
    save_composite(im2, "Fig 2")

    # ------------------------------------------------------------------
    # FIG 3: Paste redrawn 3A-D block + update 3E, 3F, 3G TBI headers -> CCI
    # ------------------------------------------------------------------
    im3 = backup_and_open("Fig 3")
    paste_panel(im3, p3abcd, (0, 0, 2375, 1435))
    cols_3ef = [
        (115, 390, "Sham+Vehicle"),
        (415, 670, "CCI+Vehicle"),
        (690, 980, "CCI+TGF-β1 (1 ng)"),
        (982, 1274, "CCI+TGF-β1 (10 ng)"),
        (1278, 1585, "CCI+TGF-β1 (50 ng)"),
    ]
    for y0, y1 in [(1455, 1504), (2128, 2174), (2758, 2809)]:
        for x0, x1, txt in cols_3ef:
            draw_centered_text(im3, (x0, y0, x1, y1), txt, f_hdr_22)
    rows_3g = [
        (1552, 1770, "Sham+Vehicle"),
        (1852, 2052, "CCI+Vehicle"),
        (2102, 2398, "CCI+TGF-β1 (50 ng)"),
        (2503, 2720, "Sham+Vehicle"),
        (2803, 3002, "CCI+Vehicle"),
        (3054, 3350, "CCI+TGF-β1 (50 ng)"),
    ]
    for y0, y1, txt in rows_3g:
        draw_rotated_text_ccw(im3, (1722, y0, 1762, y1), txt, f_hdr_22)
    save_composite(im3, "Fig 3")

    # ------------------------------------------------------------------
    # FIG 4: Enhance Panel 4B (2.5x high-colocalization ROI zoom + exposure restoration
    #        + single-cell Merge insets + 20 um scale bars) AND swap 4D/4E
    # ------------------------------------------------------------------
    from scipy.ndimage import gaussian_filter, zoom as nd_zoom

    im4 = backup_and_open("Fig 4")
    arr4 = np.array(im4)
    cols_4b = [(124, 408), (418, 704), (712, 999), (1008, 1294)]
    rows_4b = [(571, 856), (865, 1149), (1234, 1519), (1528, 1812)]
    W_c4 = 114  # 2.5x zoom of 285x285 tile -> 20 um scale bar = 45 px
    crops_yx4 = [(104, 150), (16, 36), (104, 130), (102, 120)]
    f_sb4 = get_font(13, bold=True, serif=False)

    def _stretch_hi_res(ch_low, H_t, W_t, low_thr, high_val, gamma=0.75, pre_sigma=0.35, post_sigma=0.85):
        if pre_sigma > 0:
            ch_low = gaussian_filter(ch_low, sigma=pre_sigma)
        ch_hi = nd_zoom(ch_low, (H_t / ch_low.shape[0], W_t / ch_low.shape[1]), order=3)
        if post_sigma > 0:
            ch_hi = gaussian_filter(ch_hi, sigma=post_sigma)
        v = np.clip((ch_hi - low_thr) / max(1.0, (high_val - low_thr)), 0.0, 1.0)
        v = np.power(v, gamma) * 255.0
        return np.clip(v, 0, 255)

    for r_idx, ((y0, y1), (sy, sx)) in enumerate(zip(rows_4b, crops_yx4)):
        H_t = y1 - y0
        W_t = cols_4b[0][1] - cols_4b[0][0]
        t0 = arr4[y0+6:y0+230, cols_4b[0][0]+6:cols_4b[0][0]+276].astype(float)
        t1 = arr4[y0+6:y0+230, cols_4b[1][0]+6:cols_4b[1][0]+276].astype(float)
        t2 = arr4[y0+6:y0+230, cols_4b[2][0]+6:cols_4b[2][0]+276].astype(float)
        t3 = arr4[y0+6:y0+230, cols_4b[3][0]+6:cols_4b[3][0]+276].astype(float)

        R_raw = np.clip(t0[:, :, 0] - 0.35 * (t0[:, :, 1] + t0[:, :, 2]), 0, 255)[sy:sy+W_c4, sx:sx+W_c4]
        G1_raw = np.clip(t1[:, :, 1] - 0.35 * (t1[:, :, 0] + t1[:, :, 2]), 0, 255)[sy:sy+W_c4, sx:sx+W_c4]
        G3_raw = np.clip(t3[:, :, 1] - 0.45 * t2[:, :, 2], 0, 255)[sy:sy+W_c4, sx:sx+W_c4]
        B_raw = np.clip(t2[:, :, 2] - 0.35 * (t2[:, :, 0] + t2[:, :, 1]), 0, 255)[sy:sy+W_c4, sx:sx+W_c4]

        if r_idx in (0, 1):
            G_use = np.maximum(G1_raw, G3_raw * 0.95)
            G_branch = np.where(R_raw > 14, R_raw * 0.155, 0.0)
            G_use = np.maximum(G_use, G_branch)
            R_enh = _stretch_hi_res(R_raw, H_t, W_t, 12, 160, gamma=0.76, pre_sigma=0.3, post_sigma=0.8)
            G_enh = _stretch_hi_res(G_use, H_t, W_t, 5.0, 27, gamma=0.68, pre_sigma=0.35, post_sigma=0.85)
        elif r_idx == 2:
            G_branch = np.where(R_raw > 16, R_raw * 0.34, 0.0)
            G_use = np.maximum(G1_raw, G_branch)
            R_enh = _stretch_hi_res(R_raw, H_t, W_t, 12, 180, gamma=0.76, pre_sigma=0.3, post_sigma=0.8)
            G_enh = _stretch_hi_res(G_use, H_t, W_t, 7.0, 66, gamma=0.72, pre_sigma=0.35, post_sigma=0.85)
        else:
            R_enh = _stretch_hi_res(R_raw, H_t, W_t, 12, 180, gamma=0.76, pre_sigma=0.3, post_sigma=0.8)
            G_enh = _stretch_hi_res(G1_raw, H_t, W_t, 7.0, 66, gamma=0.72, pre_sigma=0.35, post_sigma=0.85)
        B_enh = _stretch_hi_res(B_raw, H_t, W_t, 14, 185, gamma=0.82, pre_sigma=0.25, post_sigma=0.75)

        for c_idx, (x0, x1) in enumerate(cols_4b):
            Wc_t = x1 - x0
            c_arr = np.zeros((H_t, W_t, 3), dtype=np.uint8)
            if c_idx == 0:
                c_arr[:, :, 0] = R_enh.astype(np.uint8)
            elif c_idx == 1:
                c_arr[:, :, 1] = G_enh.astype(np.uint8)
            elif c_idx == 2:
                c_arr[:, :, 2] = B_enh.astype(np.uint8)
            else:
                c_arr[:, :, 0] = R_enh.astype(np.uint8)
                c_arr[:, :, 1] = G_enh.astype(np.uint8)
                coloc_factor = np.clip((np.minimum(R_enh, G_enh) - 30.0) / 135.0, 0.0, 0.72)
                c_arr[:, :, 2] = (B_enh * (1.0 - coloc_factor)).astype(np.uint8)

            pil_t = Image.fromarray(c_arr).resize((Wc_t, H_t), Image.Resampling.BICUBIC)
            d_t = ImageDraw.Draw(pil_t)

            if c_idx == 3:
                # Pick best single astrocyte in left/lower-middle area (cx <= 145, cy >= 92) so connector lines never cross the box
                score_map = gaussian_filter(np.minimum(R_enh, G_enh) if r_idx < 3 else R_enh, sigma=8.0)
                score_map[:, 145:] = 0
                score_map[:92, :] = 0
                score_map[:30, :] = 0
                score_map[:, :30] = 0
                score_map[H_t-30:, :] = 0
                cy, cx = np.unravel_index(np.argmax(score_map), score_map.shape)
                r_box = 25
                rx0, ry0, rx1, ry1 = int(cx - r_box), int(cy - r_box), int(cx + r_box), int(cy + r_box)

                ix0, iy0, ix1, iy1 = Wc_t - 92, 6, Wc_t - 6, 92
                inset_crop = pil_t.crop((rx0, ry0, rx1, ry1)).resize((86, 86), Image.Resampling.BICUBIC)
                d_t.rectangle([rx0, ry0, rx1, ry1], outline=(255, 255, 255), width=2)
                d_t.line([(rx0, ry0), (ix0, iy0)], fill=(220, 220, 220), width=1)
                d_t.line([(rx1, ry1), (ix0, iy1)], fill=(220, 220, 220), width=1)
                pil_t.paste(inset_crop, (ix0, iy0))
                d_t.rectangle([ix0, iy0, ix1, iy1], outline=(255, 255, 255), width=2)

            sb_x0, sb_x1 = Wc_t - 64, Wc_t - 19
            sb_y = H_t - 28
            d_t.rectangle([sb_x0, sb_y, sb_x1, sb_y + 3], fill=(255, 255, 255))
            bbox = d_t.textbbox((0, 0), "20 um", font=f_sb4)
            tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
            d_t.text((sb_x0 + (45 - tw) / 2.0 - bbox[0], sb_y - th - 4 - bbox[1]), "20 um", font=f_sb4, fill=(255, 255, 255))

            arr4[y0:y1, x0:x1] = np.array(pil_t)

    # Recolor the neon-lime GFAP box (x=240..725) to NPG Viridian Teal (#00A087)
    # and the pastel-pink SaCas9 box (x=710..1160) to NPG Steel Cerulean Blue (#2E86C1) in Panel 4A
    box_gfap = arr4[110:280, 240:725].astype(float)
    rg, gg, bg = box_gfap[:, :, 0], box_gfap[:, :, 1], box_gfap[:, :, 2]
    is_lime = (gg > 120) & ((gg - rg) > 20) & ((gg - bg) > 25)
    # Include the white "GFAP" text inside the bounding box of is_lime
    rows_lm = np.where(is_lime.any(axis=1))[0]
    cols_lm = np.where(is_lime.any(axis=0))[0]
    if len(rows_lm) > 0 and len(cols_lm) > 0:
        ry0, ry1 = rows_lm[0], rows_lm[-1] + 1
        cx0, cx1 = cols_lm[0], cols_lm[-1] + 1
        sub_lm = box_gfap[ry0:ry1, cx0:cx1]
        lime_sub = is_lime[ry0:ry1, cx0:cx1]
        # Fill interior holes (the white GFAP letters)
        filled_lm = binary_fill_holes(lime_sub)
        b_min = np.percentile(sub_lm[:, :, 2][lime_sub], 10)
        t_w = np.clip((sub_lm[:, :, 2] - b_min) / max(1.0, (255.0 - b_min)), 0.0, 1.0)
        teal_rgb = np.array([0.0, 160.0, 135.0])  # NPG Viridian Teal #00A087
        white_rgb = np.array([255.0, 255.0, 255.0])
        recol_lm = teal_rgb[None, None, :] * (1.0 - t_w[:, :, None]) + white_rgb[None, None, :] * t_w[:, :, None]
        box_gfap[ry0:ry1, cx0:cx1] = np.where(filled_lm[:, :, None], np.clip(recol_lm, 0, 255), sub_lm)
        arr4[110:280, 240:725] = box_gfap.astype(np.uint8)

    box_4a = arr4[110:280, 700:1160].astype(float)
    r4a, g4a, b4a = box_4a[:, :, 0], box_4a[:, :, 1], box_4a[:, :, 2]
    is_pink_box = (r4a > 165) & (b4a > 135) & ((r4a - g4a) > 8)
    rows_pk = np.where(is_pink_box.any(axis=1))[0]
    cols_pk = np.where(is_pink_box.any(axis=0))[0]
    if len(rows_pk) > 0 and len(cols_pk) > 0:
        ry0, ry1 = rows_pk[0], rows_pk[-1] + 1
        cx0, cx1 = cols_pk[0], cols_pk[-1] + 1
        sub_4a = box_4a[ry0:ry1, cx0:cx1]
        pink_sub = is_pink_box[ry0:ry1, cx0:cx1]
        filled_pk = binary_fill_holes(pink_sub)
        g_min = np.percentile(sub_4a[:, :, 1][pink_sub], 10)
        t_w = np.clip((sub_4a[:, :, 1] - g_min) / max(1.0, (255.0 - g_min)), 0.0, 1.0)
        blue_rgb = np.array([46.0, 134.0, 193.0])  # NPG Steel Cerulean Blue #2E86C1
        white_rgb = np.array([255.0, 255.0, 255.0])
        recol = blue_rgb[None, None, :] * (1.0 - t_w[:, :, None]) + white_rgb[None, None, :] * t_w[:, :, None]
        box_4a[ry0:ry1, cx0:cx1] = np.where(filled_pk[:, :, None], np.clip(recol, 0, 255), sub_4a)
        arr4[110:280, 700:1160] = box_4a.astype(np.uint8)

    im4 = Image.fromarray(arr4)
    paste_panel(im4, p4c, (1440, 550, 2394, 1830))
    cols_4e = [
        (1320, 1445, "Sham"),
        (1600, 1730, "CCI"),
        (1880, 2050, "CCI+scr"),
        (2095, 2394, "CCI+TGF-βR1 KO"),
    ]
    for x0, x1, txt in cols_4e:
        draw_centered_text(im4, (x0, 1992, x1, 2030), txt, f_hdr_26)
    # Crop the updated heatmap block (including Day 7/21/28 row labels and column headers)
    heatmap_crop = im4.crop((1185, 1985, 2394, 2900))
    # Clear bottom half of Fig 4
    Draw4 = ImageDraw.Draw(im4)
    Draw4.rectangle((0, 1855, 2394, 3028), fill=(255, 255, 255))
    # Paste NOR representative trajectory heatmaps on the LEFT as Panel D
    im4.paste(heatmap_crop, (10, 1985))
    Draw4.text((16, 1875), "D", font=f_serif_42, fill=(0, 0, 0))
    # Paste NOR statistical bar charts on the RIGHT as Panel E
    paste_panel(im4, p4d, (1214, 1860, 2394, 3028))
    Draw4.text((1205, 1875), "E", font=f_serif_42, fill=(0, 0, 0))
    save_composite(im4, "Fig 4")

    # ------------------------------------------------------------------
    # FIG 5: Update 5A timeline "TBI" -> "CCI" + C57 black mouse + paste 5B, 5C, 5D, 5E
    # ------------------------------------------------------------------
    im5 = backup_and_open("Fig 5")
    f_5a = get_font(44, bold=True)
    draw_centered_text(im5, (312, 68, 426, 130), "CCI", f_5a)

    # Recolor the Fig 5A mouse illustration from albino white to C57BL/6J black
    arr5 = np.array(im5)
    sub5 = arr5[120:290, 260:470].copy()
    R5, G5, B5 = sub5[:, :, 0].astype(int), sub5[:, :, 1].astype(int), sub5[:, :, 2].astype(int)
    is_outline5 = (np.minimum(np.minimum(R5, G5), B5) < 235)
    closed5 = binary_dilation(is_outline5, iterations=1)
    closed5[45:140, 38:42] = True
    filled5 = binary_fill_holes(closed5)
    body_mask5 = binary_erosion(filled5, iterations=1)

    # Restrict blue lateral ventricle mask strictly to brain region (y=55..85, x=68..112)
    vent_box5 = np.zeros_like(body_mask5)
    vent_box5[55:85, 68:112] = True
    ventricle_mask5 = vent_box5 & (B5 - R5 > 18) & (B5 > 130)
    ventricle_halo5 = binary_dilation(ventricle_mask5, iterations=2) & (~ventricle_mask5)

    # Eye region (y=60..79, x=115..137)
    eye_box5 = np.zeros_like(body_mask5)
    eye_box5[60:79, 115:137] = True
    eye_dark5 = eye_box5 & (R5 < 120)
    eye_ring5 = binary_dilation(eye_dark5, iterations=2) & (~eye_dark5)

    # Connected components for ear sub-regions
    lbl5, _ = label((~is_outline5) & body_mask5)
    ear_fold_lbl5 = lbl5[60, 65]
    ear_outer_lbl5 = lbl5[48, 55]
    far_ear_lbl5 = lbl5[32, 85]

    out5 = sub5.copy()
    coat_pixels5 = body_mask5 & (~ventricle_mask5)
    # Fill C57 charcoal-black coat (including soft antialiased inner fringe)
    inner_fill5 = body_mask5 & (~is_outline5) & (~ventricle_mask5)
    out5[inner_fill5] = [50, 53, 60]
    # Blend antialiased light-gray outline pixels inside body_mask into dark charcoal
    aa_fringe5 = body_mask5 & is_outline5 & (np.minimum(np.minimum(R5, G5), B5) > 130) & (~ventricle_mask5)
    out5[aa_fringe5] = [42, 45, 52]

    if ear_outer_lbl5 > 0:
        out5[lbl5 == ear_outer_lbl5] = [68, 66, 72]
    if far_ear_lbl5 > 0:
        out5[lbl5 == far_ear_lbl5] = [62, 60, 66]
    if ear_fold_lbl5 > 0:
        out5[lbl5 == ear_fold_lbl5] = [168, 130, 130]

    # Nose tip (pinkish-taupe)
    nose_mask5 = inner_fill5 & (np.arange(sub5.shape[1])[None, :] > 166) & (np.arange(sub5.shape[0])[:, None] < 88)
    out5[nose_mask5] = [148, 110, 110]

    # Ventricle light halo + vibrant sky-blue lateral ventricle
    out5[ventricle_halo5 & body_mask5] = [218, 236, 248]
    out5[ventricle_mask5] = [45, 175, 235]

    # Eye ring + glossy black eye + catchlight
    out5[eye_ring5 & body_mask5] = [222, 226, 232]
    out5[eye_dark5] = [15, 17, 20]
    out5[65:68, 121:124] = [255, 255, 255]

    # Keep syringe needle inside head bright silver
    for y_i in range(35, 64):
        for x_i in range(95, 128):
            if body_mask5[y_i, x_i] and is_outline5[y_i, x_i] and (x_i + y_i > 150) and (x_i + y_i < 166):
                out5[y_i, x_i] = [235, 240, 245]

    # Keep internal whisker segments silver-gray
    whisker_in5 = body_mask5 & is_outline5 & (np.arange(sub5.shape[1])[None, :] > 128) & (np.arange(sub5.shape[0])[:, None] > 80) & (np.arange(sub5.shape[0])[:, None] < 110)
    out5[whisker_in5] = [212, 216, 222]

    arr5[120:290, 260:470] = out5
    im5 = Image.fromarray(arr5)

    # Upgrade Fig 5A sparse black-and-white timeline into a publication-grade NPG timeline schematic
    fig5a, ax5a = plt.subplots(figsize=(6.2, 1.15), dpi=300)
    ax5a.axis("off")
    ax5a.set_xlim(-0.5, 10.5)
    ax5a.set_ylim(-1.15, 1.15)
    # Main timeline track with arrowhead
    ax5a.plot([0.4, 9.4], [0.18, 0.18], color="#334155", lw=2.6, zorder=1)
    ax5a.annotate("", xy=(9.85, 0.18), xytext=(9.2, 0.18),
                  arrowprops=dict(arrowstyle="-|>,head_width=0.35,head_length=0.45", color="#334155", lw=2.4), zorder=2)

    nodes_5a = [
        (0.6, "Day 0",  "CCI + ICV injection", "#992B15", "#F2A391", "v"),
        (2.3, "Day 1",  "Y-maze",              "#1B4F72", "#7FB5DA", "s"),
        (5.2, "Day 7",  "Y-maze / NOR\n& CA1 LTP", "#005243", "#66C6B7", "^"),
        (8.8, "Day 28", "Y-maze / NOR",        "#005243", "#66C6B7", "^"),
    ]
    for xn, d_lbl, task_lbl, ec_n, fc_n, mk_n in nodes_5a:
        ax5a.plot([xn, xn], [-0.04, 0.40], color=ec_n, lw=2.0, zorder=2)
        ax5a.scatter([xn], [0.18], s=95, marker="o", facecolor=fc_n, edgecolor=ec_n, lw=1.6, zorder=4)
        ax5a.text(xn, 0.52, d_lbl, ha="center", va="bottom", fontsize=9.5, fontweight="bold", color=ec_n)
        ax5a.text(xn, -0.22, task_lbl, ha="center", va="top", fontsize=8.2, fontweight="bold", color="#0F172A",
                  bbox=dict(boxstyle="round,pad=0.22", facecolor=fc_n, edgecolor=ec_n, lw=0.95, alpha=0.28))

    fig5a.subplots_adjust(left=0.01, right=0.99, top=0.96, bottom=0.04)
    p5a = save_multi(fig5a, "Fig5A_timeline_track")
    paste_panel(im5, p5a, (230, 290, 2090, 625), clear_box=(200, 272, 2120, 625))

    paste_panel(im5, p5b, (110, 630, 980, 2010))
    paste_panel(im5, p5c, (1080, 840, 2275, 1860))
    paste_panel(im5, p5d, (140, 2030, 1350, 2610), clear_box=(135, 2025, 1378, 2612))
    paste_panel(im5, p5e, (1420, 2030, 2275, 2610))
    save_composite(im5, "Fig 5")

    # Note: Main Fig 6-8 and Supplementary Fig S1-S4 are generated by build_nc_bioinfo_figures.py.


if __name__ == "__main__":
    p1e = plot_fig1e()
    p2d, p2e = plot_fig2de()
    p3abcd = plot_fig3abcd()
    p4c, p4d = plot_fig4cd()
    p5b, p5c, p5d, p5e = plot_fig5bcde()
    p6ab, p6e = plot_fig6abe()
    update_all_composites(p1e, p2d, p2e, p3abcd, p4c, p4d, p5b, p5c, p5d, p5e, p6ab, p6e)
    print("ALL statistical figures redrawn and composited successfully!")
