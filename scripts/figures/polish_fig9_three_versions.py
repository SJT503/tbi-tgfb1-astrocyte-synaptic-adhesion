# -*- coding: utf-8 -*-
"""
Polish and export 3 distinct BioRender-style CNS Working Model versions for Figure 9:
  - Fig 9_v1_Horizontal_MultiScale.png (.tif / .pdf)
  - Fig 9_v2_Split_Tripartite_Synapse.png (.tif / .pdf)
  - Fig 9_v3_TopDown_DualArm.png (.tif / .pdf)
  - Fig 9_cleaned.png (.tif / .pdf) [Default set to v2 Split Tripartite Synapse, easily switchable]
"""
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ART_DIR = r"d:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\JNI\figures\raw_fig9_biorender"
OUT_DIR = r"d:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\JNI\figures\cleaned_Fig1_5"

V1_SRC = os.path.join(ART_DIR, "fig9_ai_biorender_candidate_1792921501782.jpg")
V2_SRC = os.path.join(ART_DIR, "fig9_v2_split_tripartite_synapse_1792921621132.jpg")
V3_SRC = os.path.join(ART_DIR, "fig9_v3_top_to_bottom_dual_arm_1792921645130.jpg")


def get_font(size, bold=False, italic=False):
    candidates = []
    if bold and italic:
        candidates = ["arialbi.ttf", "calibriz.ttf"]
    elif bold:
        candidates = ["arialbd.ttf", "calibrib.ttf", "segoeuib.ttf"]
    elif italic:
        candidates = ["ariali.ttf", "calibrii.ttf"]
    else:
        candidates = ["arial.ttf", "calibri.ttf", "segoeui.ttf"]
    for c in candidates:
        p = os.path.join(r"C:\Windows\Fonts", c)
        if os.path.exists(p):
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()


def sample_bg(arr, x1, y1, x2, y2):
    """Sample median background color around the perimeter of a box."""
    top = arr[max(0, y1-3):y1+1, x1:x2]
    bot = arr[y2:min(arr.shape[0], y2+3), x1:x2]
    ring = np.concatenate([top.reshape(-1, 3), bot.reshape(-1, 3)], axis=0)
    med = np.median(ring, axis=0).astype(np.uint8)
    return tuple(int(v) for v in med)


def erase_box(im, draw, x1, y1, x2, y2, bg=None):
    arr = np.array(im)
    if bg is None:
        bg = sample_bg(arr, x1, y1, x2, y2)
    draw.rectangle([x1, y1, x2, y2], fill=bg)
    return bg


def draw_centered_lines(draw, cx, y_top, lines, font, fill=(25, 25, 28), line_gap=4):
    y = y_top
    for ln in lines:
        bb = draw.textbbox((0, 0), ln, font=font)
        tw = bb[2] - bb[0]
        th = bb[3] - bb[1]
        draw.text((cx - tw / 2, y), ln, font=font, fill=fill)
        y += th + line_gap


def save_multi(im_2x, stem):
    png_p = os.path.join(OUT_DIR, f"{stem}.png")
    tif_p = os.path.join(OUT_DIR, f"{stem}.tif")
    pdf_p = os.path.join(OUT_DIR, f"{stem}.pdf")
    im_2x.save(png_p, dpi=(300, 300))
    im_2x.save(tif_p, dpi=(300, 300), compression="tiff_lzw")
    im_2x.save(pdf_p, "PDF", resolution=300.0)
    print(f"[OK] Saved {stem} (.png / .tif / .pdf) -> {im_2x.size}")


def polish_v1():
    """Version 1: Horizontal 3-Stage Multi-Scale BioRender Panorama."""
    im = Image.open(V1_SRC).convert("RGB")
    # Upscale 2x for crisp 300-DPI publication output (2752 x 1536)
    im2 = im.resize((im.width * 2, im.height * 2), Image.Resampling.LANCZOS)
    save_multi(im2, "Fig 9_v1_Horizontal_MultiScale")
    return im2


def polish_v2():
    """Version 2: Split Tripartite Synapse Cross-Section (Left = TBI Pathology vs Right = Therapeutic Rescue)."""
    im = Image.open(V2_SRC).convert("RGB")
    draw = ImageDraw.Draw(im)

    f_reg20 = get_font(20, bold=False)
    f_reg19 = get_font(19, bold=False)
    f_bold19 = get_font(19, bold=True)

    # 1. Replace 'TGF-beta1' -> 'TGF-β1' (x=455..575, y=82..116)
    erase_box(im, draw, 455, 82, 578, 116)
    draw_centered_lines(draw, 514, 86, ["TGF-β1"], f_reg20, fill=(20, 20, 24))

    # 2. Replace 'TGF-beta1 NAb' -> 'TGF-β1 NAb' (x=992..1158, y=72..106)
    erase_box(im, draw, 992, 72, 1158, 106)
    draw.text((995, 76), "TGF-β1 NAb", font=f_reg20, fill=(20, 20, 24))

    # 3. Replace 'Loss of expression' -> 'Loss of Gpc5 expression' (x=182..298, y=334..392)
    erase_box(im, draw, 182, 334, 298, 392)
    draw_centered_lines(draw, 240, 336, ["Loss of Gpc5", "expression"], f_reg19, fill=(20, 20, 24), line_gap=3)

    # 4. Refine 'Homeostatic Astrocyte' -> 'Homeostatic-B Astrocyte (Apoe+)' (x=1185..1348, y=275..335)
    erase_box(im, draw, 1185, 275, 1355, 338)
    draw_centered_lines(draw, 1268, 274, ["Homeostatic-B", "Astrocyte (Apoe+)"], f_reg19, fill=(20, 20, 24), line_gap=3)

    # 5. Replace 'Mitochondrial' -> 'Mitochondrial MOMP' (x=235..398, y=662..696)
    erase_box(im, draw, 235, 662, 398, 696)
    draw_centered_lines(draw, 316, 665, ["Mitochondrial MOMP"], f_reg19, fill=(20, 20, 24))

    # 6. Flip arrow between red Cyt-c dots and 'Cytochrome c' label (x=440..506, y=604..638)
    # Sample clean background just below at y=640..648
    arr_cur = np.array(im)
    bg_left_cyto = tuple(int(v) for v in np.median(arr_cur[640:648, 445:500].reshape(-1, 3), axis=0))
    draw.rectangle([440, 604, 506, 638], fill=bg_left_cyto)
    draw.line([(446, 618), (492, 618)], fill=(35, 35, 40), width=2)
    draw.polygon([(501, 618), (489, 612), (489, 624)], fill=(35, 35, 40))

    # 7. Fix right-edge clipping of 'Postsynaptic hippocampal neuron dendritic spine' (x=1170..1376, y=506..585)
    # Sample clean background inside the right-side cytoplasm at x=1120..1165, y=520..575
    bg_right_cyto = tuple(int(v) for v in np.median(arr_cur[520:575, 1120:1165].reshape(-1, 3), axis=0))
    draw.rectangle([1170, 506, 1376, 585], fill=bg_right_cyto)
    draw_centered_lines(
        draw, 1272, 508,
        ["Postsynaptic", "hippocampal neuron", "dendritic spine"],
        f_reg19, fill=(20, 20, 24), line_gap=3
    )

    im2 = im.resize((im.width * 2, im.height * 2), Image.Resampling.LANCZOS)
    save_multi(im2, "Fig 9_v2_Split_Tripartite_Synapse")
    return im2


def polish_v3():
    """Version 3: Top-to-Bottom Dual-Arm BioRender Architecture (Direct Upgrade of Figure 6F)."""
    im = Image.open(V3_SRC).convert("RGB")
    draw = ImageDraw.Draw(im)

    f_reg19 = get_font(19, bold=False)
    f_reg18 = get_font(18, bold=False)
    f_bold22 = get_font(22, bold=True)

    # 1. Fix 'TGF-βea1' -> 'TGF-β1' at top center (x=648..752, y=20..52)
    erase_box(im, draw, 648, 20, 752, 52, bg=(255, 255, 255))
    draw_centered_lines(draw, 700, 24, ["TGF-β1"], f_reg19, fill=(20, 20, 24))

    # 2. Fix 'TGF-βea1 NAb' -> 'TGF-β1 NAb' at top right (x=1068..1235, y=45..82)
    erase_box(im, draw, 1068, 45, 1235, 82, bg=(255, 255, 255))
    draw.text((1074, 51), "TGF-β1 NAb", font=get_font(21, bold=False), fill=(20, 20, 24))

    # 3. Fix 'N-cacherin/CDH2' -> 'N-cadherin / CDH2' at center synaptic cleft (x=648..752, y=382..428)
    erase_box(im, draw, 648, 382, 752, 428, bg=(255, 255, 255))
    draw_centered_lines(draw, 700, 384, ["N-cadherin /", "CDH2"], f_reg18, fill=(20, 20, 24), line_gap=2)

    # 4. Fix 'Dendritic' -> 'Dendritic spine' at left arm bottom (x=95..195, y=582..609)
    erase_box(im, draw, 95, 582, 195, 609)
    draw.text((100, 585), "Dendritic spine", font=f_reg18, fill=(20, 20, 24))

    # 5. Refine bottom left 'Synaptic Plasticity' -> 'LTP Impairment / Rescue' (x=205..448, y=620..660)
    erase_box(im, draw, 205, 620, 448, 660, bg=(255, 255, 255))
    draw_centered_lines(draw, 326, 626, ["LTP Impairment / Rescue"], f_bold22, fill=(20, 20, 24))

    # 6. Refine bottom right 'Neuronal Anoikis' -> 'Anoikis vs Survival Rescue' (x=945..1175, y=620..660)
    erase_box(im, draw, 945, 620, 1175, 660, bg=(255, 255, 255))
    draw_centered_lines(draw, 1060, 626, ["Anoikis vs Survival"], f_bold22, fill=(20, 20, 24))

    im2 = im.resize((im.width * 2, im.height * 2), Image.Resampling.LANCZOS)
    save_multi(im2, "Fig 9_v3_TopDown_DualArm")
    return im2


if __name__ == "__main__":
    v1_im = polish_v1()
    v2_im = polish_v2()
    v3_im = polish_v3()
    # Save v2 as default Fig 9_cleaned as well
    save_multi(v2_im, "Fig 9_cleaned")
