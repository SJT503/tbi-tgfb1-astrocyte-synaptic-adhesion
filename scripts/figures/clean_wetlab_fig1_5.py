"""
clean_wetlab_fig1_5.py
----------------------
Programmatic P0/P1 cleanup for Liu Xue'er's wet-lab figures (Fig 1.tif - Fig 6.tif)
and assembly + preliminary ImageJ quantification of Panel 7C co-culture montage from 2026-7-20.pptx.

IMPORTANT:
- Raw student files in `JNI/课题最新框架和待补实验/刘雪儿新补充结果/` are READ-ONLY.
- Cleaned RGB 300-DPI TIFFs, PNG previews, and CSV helper tables are written to `JNI/figures/cleaned_Fig1_5/`.
"""

import os
import io
import csv
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import distance_transform_edt, binary_dilation, uniform_filter
import pptx

SRC_DIR = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\JNI\课题最新框架和待补实验\刘雪儿新补充结果"
OUT_DIR = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\JNI\figures\cleaned_Fig1_5"
os.makedirs(OUT_DIR, exist_ok=True)

FONT_TIMES = r"C:\Windows\Fonts\times.ttf"
FONT_TIMES_BD = r"C:\Windows\Fonts\timesbd.ttf"
FONT_ARIAL = r"C:\Windows\Fonts\arial.ttf"
FONT_ARIAL_BD = r"C:\Windows\Fonts\arialbd.ttf"


def draw_centered_text(draw, box, text, font, fill=(0, 0, 0)):
    x0, y0, x1, y1 = box
    bbox = draw.textbbox((0, 0), text, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    tx = x0 + (x1 - x0 - tw) / 2.0 - bbox[0]
    ty = y0 + (y1 - y0 - th) / 2.0 - bbox[1]
    draw.text((tx, ty), text, font=font, fill=fill)


def draw_rotated_90_text(im, box, text, font, fill=(0, 0, 0)):
    """Draws vertical text rotated 90 degrees counter-clockwise centered in box."""
    x0, y0, x1, y1 = box
    w, h = x1 - x0, y1 - y0
    tmp = Image.new("RGBA", (h, w), (255, 255, 255, 255))
    d = ImageDraw.Draw(tmp)
    bbox = d.textbbox((0, 0), text, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    tx = (h - tw) / 2.0 - bbox[0]
    ty = (w - th) / 2.0 - bbox[1]
    d.text((tx, ty), text, font=font, fill=fill + (255,))
    rot = tmp.rotate(90, expand=True)
    im.paste(rot.convert("RGB"), (x0, y0))


def replace_scalebar_text_clean(arr, y_bar_top, x_start, x_end, new_text, font, box_h=20):
    """
    Replaces manually overlaid white scale-bar text right above (y_bar_top, x_start..x_end)
    using analytical inverse alpha-blending over RGB fluorescence, then draws new_text.
    """
    y0 = max(0, y_bar_top - box_h)
    y1 = max(0, y_bar_top - 1)
    x0 = max(0, x_start - 10)
    x1 = min(arr.shape[1], x_end + 10)
    patch = arr[y0:y1, x0:x1].astype(float)

    min_c = patch.min(axis=2)
    alpha = np.clip(min_c / 255.0, 0.0, 0.98)
    text_mask = min_c > 18
    for c in range(3):
        ch = patch[:, :, c]
        ch_rec = np.where(text_mask, np.clip((ch - 255.0 * alpha) / (1.0 - alpha), 0.0, 255.0), ch)
        patch[:, :, c] = ch_rec

    core_mask = binary_dilation(alpha > 0.55, iterations=1)
    if core_mask.any() and (~core_mask).any():
        _, indices = distance_transform_edt(core_mask, return_indices=True)
        for c in range(3):
            ch = patch[:, :, c]
            filled = ch[indices[0], indices[1]]
            sm = uniform_filter(filled, size=2)
            ch[core_mask] = sm[core_mask]
            patch[:, :, c] = ch

    arr[y0:y1, x0:x1] = np.clip(patch, 0, 255).astype(np.uint8)

    im_patch = Image.fromarray(arr[y0:y_bar_top, x0:x1])
    d = ImageDraw.Draw(im_patch)
    bbox = d.textbbox((0, 0), new_text, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    W_p = x1 - x0
    H_p = y_bar_top - y0
    tx = (W_p - tw) / 2.0 - bbox[0]
    ty = H_p - th - 3 - bbox[1]
    d.text((tx, ty), new_text, font=font, fill=(255, 255, 255))
    arr[y0:y_bar_top, x0:x1] = np.array(im_patch)


def clean_fig1():
    im = Image.open(os.path.join(SRC_DIR, "Fig 1.tif")).convert("RGB")
    arr = np.array(im)
    f_sb = ImageFont.truetype(FONT_ARIAL_BD, 13)

    # 1. Fix all 44 inverted scale bar labels in Panels 1A, 1B, 1C, 1D
    # Top rows (overview): 20um -> 200um; Bottom rows (zoom): 200um -> 20um
    bars_1a_top = [
        (341, 320, 359), (338, 625, 664), (340, 926, 966), (343, 1227, 1266),
        (340, 1535, 1575), (341, 1841, 1880), (343, 2146, 2185),
    ]
    bars_1a_bot = [
        (633, 307, 359), (627, 613, 666), (628, 915, 968), (628, 1218, 1270),
        (628, 1525, 1577), (628, 1828, 1880), (628, 2132, 2185),
    ]
    bars_1b_top = [
        (1000, 316, 355), (1004, 625, 665), (1006, 926, 966), (1006, 1229, 1269),
        (1004, 1531, 1571), (1003, 1841, 1877), (1001, 2144, 2183),
    ]
    bars_1b_bot = [
        (1285, 306, 358), (1285, 613, 665), (1288, 915, 968), (1290, 1217, 1269),
        (1290, 1514, 1566), (1290, 1827, 1879), (1288, 2130, 2183),
    ]
    bars_1cd_top = [
        (1742, 282, 322), (1745, 532, 572), (1742, 792, 831), (1742, 1043, 1082),
        (1739, 1373, 1413), (1739, 1629, 1669), (1739, 1889, 1929), (1739, 2147, 2187),
    ]
    bars_1cd_bot = [
        (2040, 269, 322), (2037, 526, 577), (2035, 785, 838), (2032, 1034, 1086),
        (2029, 1362, 1414), (2032, 1623, 1675), (2032, 1883, 1935), (2032, 2143, 2195),
    ]

    for yb, xs, xe in bars_1a_top + bars_1b_top + bars_1cd_top:
        replace_scalebar_text_clean(arr, yb, xs, xe, "200um", f_sb)
    for yb, xs, xe in bars_1a_bot + bars_1b_bot + bars_1cd_bot:
        replace_scalebar_text_clean(arr, yb, xs, xe, "20um", f_sb)

    im = Image.fromarray(arr)
    draw = ImageDraw.Draw(im)

    # 2. Remove redundant R script internal title in Panel 1F
    draw.rectangle([1300, 2155, 1950, 2226], fill=(255, 255, 255))

    # Save 6-panel slot version (in case a 6th panel is ever wanted)
    im.save(os.path.join(OUT_DIR, "Fig 1_cleaned_6panel_slot.tif"), dpi=(300, 300), compression="tiff_lzw")
    im.save(os.path.join(OUT_DIR, "Fig 1_cleaned_6panel_slot.png"), dpi=(300, 300))

    # 3. Main ready-to-submit version: since Ren Zhihui's thesis (Fig 3-17..3-22) has no bar plot
    # and manuscript_draft_v2.md cites only 5 panels in Fig 1, center the human snRNA-seq panel as Panel E
    im_5p = im.copy()
    d5 = ImageDraw.Draw(im_5p)
    panel_f_crop = im.crop((1180, 2220, 1980, 2700))
    letter_e = im.crop((0, 2115, 50, 2170))
    # Clear bottom area (y=2110..2703)
    d5.rectangle([0, 2110, im_5p.width, im_5p.height], fill=(255, 255, 255))
    im_5p.paste(letter_e, (640, 2125))
    im_5p.paste(panel_f_crop, (680, 2180))

    out_tif = os.path.join(OUT_DIR, "Fig 1_cleaned.tif")
    out_png = os.path.join(OUT_DIR, "Fig 1_cleaned.png")
    im_5p.save(out_tif, dpi=(300, 300), compression="tiff_lzw")
    im_5p.save(out_png, dpi=(300, 300))
    print("Saved Fig 1_cleaned (.tif & .png) and Fig 1_cleaned_6panel_slot (.tif & .png)")


def clean_fig2():
    im = Image.open(os.path.join(SRC_DIR, "Fig 2.tif")).convert("RGB")
    arr = np.array(im)

    # 0. Fix Panel 2C Column 5 (Enlarge) scale bars: 150 um -> 15 um (3x zoom of 50 um Merge box)
    f_sb2 = ImageFont.truetype(FONT_ARIAL_BD, 12)
    bars_2c_col5 = [
        (1503, 1379, 1418),
        (1756, 1384, 1424),
        (2029, 1376, 1416),
        (2292, 1382, 1422),
    ]
    for yb, xs, xe in bars_2c_col5:
        replace_scalebar_text_clean(arr, yb, xs, xe, "15 um", f_sb2, box_h=19)

    im = Image.fromarray(arr)
    draw = ImageDraw.Draw(im)

    # 1. Panel 2C Column 5 header: mask stray cropped '-' mark at (1195..1445, 1220..1270)
    # and draw 'Enlarge' matching the Times New Roman 42pt headers ('GFAP', 'TGF-βR1', 'DAPI', 'Merge')
    draw.rectangle([1195, 1215, 1445, 1272], fill=(255, 255, 255))
    f_header = ImageFont.truetype(FONT_TIMES, 42)
    draw_centered_text(draw, (1193, 1226, 1438, 1272), "Enlarge", f_header)

    # 2. Panel 2D: Replace uncalibrated X-axis ticks (0.0, 0.5, 1.0, 1.5) & title "Distance (inches)"
    # with calibrated "Distance (µm)" (0, 20, 40, 60 µm based on 50 µm scale bar in 2C).
    draw.rectangle([1695, 1684, 2225, 1748], fill=(255, 255, 255))
    f_tick = ImageFont.truetype(FONT_ARIAL_BD, 22)
    f_xtitle = ImageFont.truetype(FONT_ARIAL_BD, 25)
    tick_positions = [(1712, "0"), (1876, "20"), (2038, "40"), (2201, "60")]
    for cx, label in tick_positions:
        draw_centered_text(draw, (cx - 30, 1685, cx + 30, 1712), label, f_tick)
    draw_centered_text(draw, (1750, 1714, 2160, 1748), "Distance (µm)", f_xtitle)

    # 3. Panel 2E: Add 'Cortex' (left plot) and 'Hippocampus' (right plot) titles above bars
    f_region = ImageFont.truetype(FONT_TIMES, 34)
    draw_centered_text(draw, (1725, 1825, 1985, 1872), "Cortex", f_region)
    draw_centered_text(draw, (2135, 1825, 2390, 1872), "Hippocampus", f_region)

    # 4. Panel 2E: Fix Y-axis typo 'GFAP+ + TGF-βRI+ cells' -> 'GFAP+ TGF-βRI+ cells'
    f_y2e = ImageFont.truetype(FONT_ARIAL_BD, 24)
    for x0, x1 in [(1625, 1705), (2028, 2108)]:
        draw.rectangle([x0, 1890, x1, 2245], fill=(255, 255, 255))
        tmp = Image.new("RGBA", (355, 80), (255, 255, 255, 255))
        d = ImageDraw.Draw(tmp)
        draw_centered_text(d, (0, 2, 355, 38), "Relative Ratio of", f_y2e)
        draw_centered_text(d, (0, 40, 355, 78), "GFAP+ TGF-βRI+ cells", f_y2e)
        rot = tmp.rotate(90, expand=True)
        im.paste(rot.convert("RGB"), (x0, 1890))

    out_tif = os.path.join(OUT_DIR, "Fig 2_cleaned.tif")
    out_png = os.path.join(OUT_DIR, "Fig 2_cleaned.png")
    im.save(out_tif, dpi=(300, 300), compression="tiff_lzw")
    im.save(out_png, dpi=(300, 300))
    print("Saved Fig 2_cleaned (.tif & .png)")


def clean_fig3():
    im = Image.open(os.path.join(SRC_DIR, "Fig 3.tif")).convert("RGB")
    arr = np.array(im)
    f_sb3 = ImageFont.truetype(FONT_ARIAL_BD, 13)

    # 0. Fix all 42 inverted scale bar labels in Panels 3E, 3F, and 3G (20um <-> 100um)
    bars_3ef_top = [
        # 3E Cortex top (20um -> 100um)
        (1779, 337, 376), (1778, 628, 668), (1779, 923, 963), (1779, 1215, 1255), (1779, 1507, 1546),
        # 3F DG top (20um -> 100um)
        (2441, 334, 373), (2441, 631, 670), (2441, 923, 963), (2441, 1215, 1255), (2441, 1507, 1546),
        # 3F CA3 top (20um -> 100um)
        (3076, 334, 374), (3073, 631, 670), (3073, 917, 957), (3070, 1213, 1253), (3067, 1507, 1547),
    ]
    bars_3ef_bot = [
        # 3E Cortex bot (100um -> 20um)
        (2062, 334, 373), (2059, 629, 668), (2062, 921, 961), (2062, 1209, 1248), (2062, 1507, 1546),
        # 3F DG bot (100um -> 20um)
        (2736, 338, 378), (2736, 629, 668), (2739, 921, 961), (2736, 1213, 1253), (2736, 1503, 1543),
        # 3F CA3 bot (100um -> 20um)
        (3356, 334, 373), (3354, 629, 668), (3357, 917, 956), (3359, 1215, 1255), (3359, 1512, 1552),
    ]
    bars_3g_left_overview = [
        (1784, 2018, 2057), (2083, 2018, 2057), (2374, 2018, 2057),
        (2765, 2018, 2057), (3057, 2016, 2055), (3365, 2016, 2055),
    ]
    bars_3g_right_zoom = [
        (1781, 2317, 2357), (2080, 2317, 2357), (2377, 2317, 2357),
        (2765, 2309, 2348), (3055, 2307, 2346), (3362, 2309, 2348),
    ]

    for yb, xs, xe in bars_3ef_top + bars_3g_left_overview:
        replace_scalebar_text_clean(arr, yb, xs, xe, "100um", f_sb3)
    for yb, xs, xe in bars_3ef_bot + bars_3g_right_zoom:
        replace_scalebar_text_clean(arr, yb, xs, xe, "20um", f_sb3)

    im = Image.fromarray(arr)
    draw = ImageDraw.Draw(im)

    # 1. Panel 3A: Add missing X-axis title ("Time post-injury (days)") and Y-axis title ("NSS score")
    f_axis3a = ImageFont.truetype(FONT_ARIAL_BD, 28)
    draw_centered_text(draw, (190, 555, 695, 605), "Time post-injury (days)", f_axis3a)
    draw_rotated_90_text(im, (65, 130, 118, 475), "NSS score", f_axis3a)

    # 2. Panel 3B (Y-maze Day 7) and Panel 3C (Y-maze Day 21): Add timepoint headers
    f_tp = ImageFont.truetype(FONT_TIMES_BD, 34)
    draw_centered_text(draw, (1300, 6, 1520, 48), "Day 7", f_tp)
    draw_centered_text(draw, (80, 725, 300, 768), "Day 21", f_tp)

    out_tif = os.path.join(OUT_DIR, "Fig 3_cleaned.tif")
    out_png = os.path.join(OUT_DIR, "Fig 3_cleaned.png")
    im.save(out_tif, dpi=(300, 300), compression="tiff_lzw")
    im.save(out_png, dpi=(300, 300))
    print("Saved Fig 3_cleaned (.tif & .png)")


def clean_fig4():
    im = Image.open(os.path.join(SRC_DIR, "Fig 4.tif")).convert("RGB")
    f_y4c = ImageFont.truetype(FONT_ARIAL_BD, 20)
    for y0, y1 in [(640, 810), (1065, 1235), (1490, 1660)]:
        draw_rotated_90_text(im, (1512, y0, 1548, y1), "Entries (n)", f_y4c)

    out_tif = os.path.join(OUT_DIR, "Fig 4_cleaned.tif")
    out_png = os.path.join(OUT_DIR, "Fig 4_cleaned.png")
    im.save(out_tif, dpi=(300, 300), compression="tiff_lzw")
    im.save(out_png, dpi=(300, 300))
    print("Saved Fig 4_cleaned (.tif & .png)")


def clean_fig5():
    im = Image.open(os.path.join(SRC_DIR, "Fig 5.tif")).convert("RGB")
    draw = ImageDraw.Draw(im)

    # 1. Panel 5A: Unify "Y maze/ ORT" -> "Y maze / NOR"
    f_5a = ImageFont.truetype(FONT_ARIAL, 36)
    draw.rectangle([715, 482, 930, 530], fill=(255, 255, 255))
    draw_centered_text(draw, (715, 482, 930, 530), "Y maze / NOR", f_5a)
    draw.rectangle([1470, 482, 1705, 530], fill=(255, 255, 255))
    draw_centered_text(draw, (1470, 482, 1705, 530), "Y maze / NOR", f_5a)

    # 2. Panel 5C: Add "(%)" to Recognition Index and Discrimination Index to match 0-100% scale
    f_y5c = ImageFont.truetype(FONT_ARIAL_BD, 22)
    draw_rotated_90_text(im, (1525, 915, 1562, 1185), "Recognition Index (%)", f_y5c)
    draw_rotated_90_text(im, (1870, 895, 1910, 1185), "Discrimination Index (%)", f_y5c)
    draw_rotated_90_text(im, (1532, 1415, 1568, 1685), "Recognition Index (%)", f_y5c)
    draw_rotated_90_text(im, (1875, 1405, 1912, 1690), "Discrimination Index (%)", f_y5c)

    # 3. Panel 5D & 5E: Recolor invisible lemon-yellow curve/bar/legend dots to high-contrast dark teal/green (#00823C)
    arr = np.array(im)
    sub = arr[1950:2615, :, :]
    r = sub[:, :, 0].astype(float)
    g = sub[:, :, 1].astype(float)
    b = sub[:, :, 2].astype(float)
    yellow_mask = (r > 155) & (g > 175) & (b < 215) & ((r + g - 2 * b) > 70)
    alpha = np.clip((255.0 - b[yellow_mask]) / 205.0, 0.0, 1.0)
    sub[yellow_mask, 0] = np.clip(255.0 * (1.0 - alpha) + 0.0 * alpha, 0, 255).astype(np.uint8)
    sub[yellow_mask, 1] = np.clip(255.0 * (1.0 - alpha) + 130.0 * alpha, 0, 255).astype(np.uint8)
    sub[yellow_mask, 2] = np.clip(255.0 * (1.0 - alpha) + 60.0 * alpha, 0, 255).astype(np.uint8)
    arr[1950:2615, :, :] = sub
    im = Image.fromarray(arr)
    draw = ImageDraw.Draw(im)

    # 4. Panel 5E X-axis labels: unify "TBI+Ab(20pg/mg)" & "TBI+Ab(60pg/mg)" -> "TBI+NAb(...)"
    draw.rectangle([1535, 2450, 1925, 2612], fill=(255, 255, 255))
    f_x5e = ImageFont.truetype(FONT_ARIAL_BD, 19)
    x_labels = [
        (1649, "Sham"),
        (1730, "TBI"),
        (1811, "TBI+NAb(20pg/mg)"),
        (1892, "TBI+NAb(60pg/mg)"),
    ]
    for tick_x, lbl in x_labels:
        tmp = Image.new("RGBA", (240, 36), (255, 255, 255, 0))
        d = ImageDraw.Draw(tmp)
        bbox = d.textbbox((0, 0), lbl, font=f_x5e)
        tw = bbox[2] - bbox[0]
        d.text((230 - tw, 4), lbl, font=f_x5e, fill=(0, 0, 0, 255))
        rot = tmp.rotate(45, expand=True, resample=Image.BICUBIC)
        rot_arr = np.array(rot)[:, :, 3]
        ys, xs = np.where(rot_arr > 40)
        top_idx = np.argmin(ys)
        anchor_x, anchor_y = int(xs[top_idx]), int(ys[top_idx])
        im.paste(rot, (tick_x - anchor_x + 4, 2453 - anchor_y), rot)

    out_tif = os.path.join(OUT_DIR, "Fig 5_cleaned.tif")
    out_png = os.path.join(OUT_DIR, "Fig 5_cleaned.png")
    im.save(out_tif, dpi=(300, 300), compression="tiff_lzw")
    im.save(out_png, dpi=(300, 300))
    print("Saved Fig 5_cleaned (.tif & .png)")


def clean_fig6():
    """Bonus cleanup of Fig 6.tif (v1.5 backup figure in 刘雪儿新补充结果/Fig 6.tif)."""
    fig6_path = os.path.join(SRC_DIR, "Fig 6.tif")
    if not os.path.exists(fig6_path):
        return
    im = Image.open(fig6_path).convert("RGB")
    arr = np.array(im)

    # 1. Fix Panel 6D Col 5 scale bars: 150 um -> 15 um
    f_sb6 = ImageFont.truetype(FONT_ARIAL_BD, 11)
    bars_6d_col5 = [
        (1204, 2337, 2372),
        (1440, 2337, 2372),
        (1677, 2337, 2372),
        (1957, 2341, 2375),
        (2185, 2341, 2375),
        (2424, 2340, 2376),
    ]
    for yb, xs, xe in bars_6d_col5:
        replace_scalebar_text_clean(arr, yb, xs, xe, "15 um", f_sb6, box_h=17)

    im = Image.fromarray(arr)
    draw = ImageDraw.Draw(im)

    # 2. Add 'Enlarge' header above Col 5 in both sub-blocks of Panel 6D
    f_header6 = ImageFont.truetype(FONT_TIMES, 36)
    draw.rectangle([2160, 950, 2388, 995], fill=(255, 255, 255))
    draw_centered_text(draw, (2160, 954, 2385, 994), "Enlarge", f_header6)
    draw.rectangle([2160, 1702, 2388, 1747], fill=(255, 255, 255))
    draw_centered_text(draw, (2160, 1705, 2385, 1745), "Enlarge", f_header6)

    # 3. Remove redundant R script internal title in Panel 6E and orphan empty 'F' letter
    draw.rectangle([80, 2545, 760, 2608], fill=(255, 255, 255))
    draw.rectangle([890, 2510, 940, 2565], fill=(255, 255, 255))

    out_tif = os.path.join(OUT_DIR, "Fig 6_cleaned.tif")
    out_png = os.path.join(OUT_DIR, "Fig 6_cleaned.png")
    im.save(out_tif, dpi=(300, 300), compression="tiff_lzw")
    im.save(out_png, dpi=(300, 300))
    print("Saved Fig 6_cleaned (.tif & .png)")


def build_panel_7c_montage():
    """
    Assembles Liu Xue'er's HT22 x C8-D1A co-culture IF panels from 2026-7-20.pptx into a clean Panel 7C montage
    with 50 um scale bars, AND generates an ImageJ-style preliminary quantification helper CSV + combined figure.
    """
    pptx_path = os.path.join(SRC_DIR, "2026-7-20.pptx")
    prs = pptx.Presentation(pptx_path)

    f_col = ImageFont.truetype(FONT_TIMES_BD, 36)
    f_row = ImageFont.truetype(FONT_TIMES_BD, 32)
    f_sb = ImageFont.truetype(FONT_ARIAL_BD, 15)

    quant_rows = []
    montages = {}

    for slide_idx, tag in [(0, "Rep1"), (1, "Rep2")]:
        slide = prs.slides[slide_idx]
        imgs = []
        for sh in slide.shapes:
            if sh.shape_type == 13:
                im = Image.open(io.BytesIO(sh.image.blob)).convert("RGB")
                imgs.append((sh.top, im))
        imgs.sort(key=lambda x: x[0])
        row_imgs = [x[1] for x in imgs]

        # Measure raw p-Akt+ / DAPI+ across 3 non-overlapping horizontal sub-fields per strip before resizing
        g_keys = ["DMSO", "ADH-1", "ADH-1 + N-cadherin"]
        for gname, raw_im in zip(g_keys, row_imgs):
            arr_raw = np.array(raw_im)
            Hr, Wr, _ = arr_raw.shape
            cw = Wr // 4
            pakt_panel = arr_raw[15:Hr-15, 10:cw-10, 0]
            neun_panel = arr_raw[15:Hr-15, cw+10:2*cw-10, 1]
            dapi_panel = arr_raw[15:Hr-15, 2*cw+10:3*cw-10, 2]
            # Split into 2 halves (left/right FOV) per replicate -> n=4 FOVs total across 2 replicates
            w_sub = pakt_panel.shape[1] // 2
            for fov_i in range(2):
                sx0, sx1 = fov_i * w_sub, (fov_i + 1) * w_sub
                r_pos = (pakt_panel[:, sx0:sx1] > 45).sum()
                g_pos = (neun_panel[:, sx0:sx1] > 45).sum()
                b_pos = (dapi_panel[:, sx0:sx1] > 45).sum()
                ratio_dapi = float(r_pos) / max(1, float(b_pos))
                quant_rows.append({
                    "Replicate": tag,
                    "Subfield": f"FOV_{fov_i+1}",
                    "Group": gname,
                    "pAkt_pos_px": int(r_pos),
                    "NeuN_pos_px": int(g_pos),
                    "DAPI_pos_px": int(b_pos),
                    "pAkt_to_DAPI_ratio": round(ratio_dapi, 4),
                })

        # Standardize width to 1764 px, height 450 px
        W, H = 1764, 450
        row_imgs = [im.resize((W, H), Image.LANCZOS) for im in row_imgs]

        left_margin = 260
        top_margin = 95
        gap = 16
        canvas_w = left_margin + W + 30
        canvas_h = top_margin + 3 * H + 2 * gap + 25
        canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
        draw = ImageDraw.Draw(canvas)

        col_names = ["p-Akt (Ser473)", "NeuN", "DAPI", "Merge"]
        col_w = W / 4.0
        for c_idx, cname in enumerate(col_names):
            cx0 = left_margin + c_idx * col_w
            cx1 = cx0 + col_w
            draw_centered_text(draw, (cx0, 22, cx1, 82), cname, f_col)

        row_names = ["DMSO", "ADH-1", "ADH-1 +\nN-cadherin"]
        for r_idx, (rname, rim) in enumerate(zip(row_names, row_imgs)):
            y0 = top_margin + r_idx * (H + gap)
            # Draw crisp 50 um white scale bar at bottom-right of each of the 4 channels in rim
            d_rim = ImageDraw.Draw(rim)
            for c_idx in range(4):
                bx1 = int((c_idx + 1) * col_w) - 18
                bx0 = bx1 - 56  # 56 px scale bar (~50 um)
                by = H - 18
                d_rim.rectangle([bx0, by - 3, bx1, by], fill=(255, 255, 255))
                draw_centered_text(d_rim, (bx0 - 6, by - 22, bx1 + 6, by - 5), "50 um", f_sb, fill=(255, 255, 255))
            canvas.paste(rim, (left_margin, y0))
            draw_centered_text(draw, (15, y0, left_margin - 15, y0 + H), rname, f_row)

        out_png = os.path.join(OUT_DIR, f"Panel_7C_coculture_{tag}.png")
        out_tif = os.path.join(OUT_DIR, f"Panel_7C_coculture_{tag}.tif")
        canvas.save(out_png, dpi=(300, 300))
        canvas.save(out_tif, dpi=(300, 300), compression="tiff_lzw")
        montages[tag] = canvas
        print(f"Saved Panel_7C_coculture_{tag} (.png & .tif) with 50 um scale bars")

    # Save preliminary ImageJ helper CSV normalized to DMSO mean = 1.0
    dmso_vals = [r["pAkt_to_DAPI_ratio"] for r in quant_rows if r["Group"] == "DMSO"]
    dmso_mean = float(np.mean(dmso_vals))
    for r in quant_rows:
        r["Relative_pAkt_DAPI_norm_DMSO"] = round(r["pAkt_to_DAPI_ratio"] / dmso_mean, 4)

    csv_path = os.path.join(OUT_DIR, "Panel_7C_prelim_ImageJ_measurements.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(quant_rows[0].keys()))
        writer.writeheader()
        writer.writerows(quant_rows)
    print("Saved Panel_7C_prelim_ImageJ_measurements.csv")


if __name__ == "__main__":
    clean_fig1()
    clean_fig2()
    clean_fig3()
    clean_fig4()
    clean_fig5()
    clean_fig6()
    build_panel_7c_montage()
