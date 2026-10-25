#!/usr/bin/env Rscript
# style_lib_v2.R — JNI 投稿统一可视化风格 (publication-ready, colorblind-safe)
# 用法: source("style_lib_v2.R") 然后用 jni_theme(), jni_palette
suppressMessages({library(ggplot2)})

# ── PALETTES (colorblind-safe) ──────────────────────────────────────
JNI_COL <- list(
  control = "#4DAECD",   # teal-blue (sham/control)
  tbi     = "#E07B54",   # warm orange (TBI/CCI)
  ko      = "#6A9A5B",   # sage green (TGFβRI KO / rescue)
  neutral = "#9099A8",   # cool grey (intermediate/neutral)
  hippo   = "#7AAFE6",   # light blue (hippocampus)
  cortex  = "#B07AA1",   # mauve (cortex)
  human   = "#D4A24C",   # gold (human cohort)
  up      = "#C44E52",   # crimson (up-regulated)
  down    = "#4878D0"    # blue (down-regulated)
)

# canonical 3-group palette for Sham/CCI/KO
jni_palette_3 <- c(Sham = JNI_COL$control,
                   CCI_Scr = JNI_COL$tbi,
                   CCI_TGFbetaR1_KO = JNI_COL$ko)

# canonical 2-group palette for control/TBI
jni_palette_2 <- c(control = JNI_COL$control,
                   tbi = JNI_COL$tbi)

# brain region
jni_palette_region <- c(hippocampus = JNI_COL$hippo,
                        cortex = JNI_COL$cortex)

# ── THEMES ──────────────────────────────────────────────────────────
jni_theme <- function(base = 9, legend = "top") {
  theme_bw(base_size = base) +
    theme(plot.title       = element_text(face = "bold", size = base * 1.08),
          plot.subtitle    = element_text(size = base * 0.82, colour = "grey35"),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text       = element_text(face = "bold", size = base * 0.84),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.key.size  = unit(3.2, "mm"),
          legend.title     = element_blank(),
          legend.position  = legend,
          axis.title       = element_text(size = base * 0.95))
}

jni_theme_violin <- function(base = 9) {
  jni_theme(base) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

# ── COMMON LABELS ───────────────────────────────────────────────────
jni_labels <- list(
  panel_letter = function(x) labs(tag = x) + theme(plot.tag = element_text(face = "bold", size = 11)),
  star_sig = function(p) ifelse(p < 0.001, "***",
                          ifelse(p < 0.01,  "**",
                          ifelse(p < 0.05,  "*", "ns")))
)

# ── FIGURE OUTPUT HELPERS ───────────────────────────────────────────
jni_save <- function(plot, file, w = 190, h = 150, unit = "mm", dpi = 300) {
  ggsave(paste0(file, ".pdf"), plot, width = w, height = h, units = unit)
  ggsave(paste0(file, ".png"), plot, width = w, height = h, units = unit, dpi = dpi)
}

# usage example -----------------------------------------------------
# source("D:/.../JNI/scripts/style_lib_v2.R")
# p <- ggplot(df, aes(x, y, fill=group)) + geom_col() +
#      scale_fill_manual(values=jni_palette_3) +
#      labs(title="A. ...") +
#      jni_theme()
# jni_save(p, "JNI/figures/FigN_panel_A", w=95, h=70)
