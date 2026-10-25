# -*- coding: utf-8 -*-
"""Nature Communications (NC) Submission — 3 Word Tables (100% aligned with Fig. 1–9 & Supplementary Fig. S1–S4)."""
import os
import re
import shutil
from docx import Document
from docx.shared import Pt, Cm

OUT = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\JNI\tables"
OUT_SYNC = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\JNI\课题最新框架和待补实验\表格"
os.makedirs(OUT, exist_ok=True)
os.makedirs(OUT_SYNC, exist_ok=True)


def set_cell_text_formatted(cell, text_str, font_size=8):
    """Writes text into a table cell while converting ^{...} into native Word superscript runs."""
    cell.text = ""
    p = cell.paragraphs[0]
    parts = re.split(r"(\^\{[^}]+\})", str(text_str))
    for part in parts:
        if not part:
            continue
        if part.startswith("^{") and part.endswith("}"):
            r = p.add_run(part[2:-1])
            r.font.superscript = True
            r.font.size = Pt(font_size)
        else:
            r = p.add_run(part)
            r.font.size = Pt(font_size)


def make_table(fname, title, caption, headers, rows, widths=None):
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Arial"
    st.font.size = Pt(9)
    h = doc.add_paragraph()
    r = h.add_run(title)
    r.bold = True
    r.font.size = Pt(11)
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = "Table Grid"
    for i, hd in enumerate(headers):
        c = t.rows[0].cells[i]
        c.text = hd
        for p in c.paragraphs:
            for run in p.runs:
                run.bold = True
                run.font.size = Pt(8.5)
    for row in rows:
        cells = t.add_row().cells
        for i, val in enumerate(row):
            set_cell_text_formatted(cells[i], val, font_size=8)
    if widths:
        for row in t.rows:
            for i, w in enumerate(widths):
                row.cells[i].width = Cm(w)
    cap = doc.add_paragraph()
    cr = cap.add_run(caption)
    cr.font.size = Pt(8)
    cr.italic = True
    fp = os.path.join(OUT, fname)
    doc.save(fp)
    shutil.copy2(fp, os.path.join(OUT_SYNC, fname))
    print("saved & synced:", fp)


# ========== Table 1 — Experimental design & multi-cohort architecture (Main Table 1 for JNI) ==========
t1_head = [
    "Cohort / Experiment", "Model / System", "Groups", "Sample size",
    "Intervention / Design", "Time points", "Key readouts (Figures)", "Data / Source"
]
t1_rows = [
    [
        "1. Endogenous TGF-β1 time-course & receptor localization",
        "Mouse CCI (C57BL/6J)",
        "Sham; CCI 0.5, 1, 3, 7, 14, 30 d; Sham vs CCI 7 d co-localization",
        "n = 3–6 mice/group",
        "Natural post-injury time-course (3.0 m/s, 2.0 mm depth, 500 ms dwell)",
        "0.5, 1, 3, 7, 14, 30 dpi",
        "Cortical & hippocampal TGF-β1 IF; TGF-β1 + Iba1/GFAP double IF (Fig. 1A–D); TGF-βRI + NeuN/Iba1/GFAP double IF & WB (Fig. 4A–E)",
        "In-house wet-lab (SUMC 2021-134)"
    ],
    [
        "2. Exogenous rmTGF-β1 ICV dose-response",
        "Mouse CCI (C57BL/6J)",
        "Sham+Vehicle; CCI+Vehicle; CCI+rmTGF-β1 (1, 10, 50 ng / 25 g)",
        "n = 6 (NSS); n = 3–6 (IF/behavior)",
        "ICV stereotaxic injection (AP +1.5, ML −1.1, DV −2.0 mm) at injury",
        "NSS −1 to 30 dpi; Y-maze 7, 21 dpi; NOR 30 dpi; IF 30 dpi",
        "NSS; Y-maze; NOR; cortical & hippocampal CA3/DG neuronal loss (NeuN) & gliosis (Iba1/GFAP) (Fig. 2A–G)",
        "In-house wet-lab (SUMC 2021-134)"
    ],
    [
        "3. Anti-TGF-β1 neutralizing mAb (immediate & delayed LTP cohorts)",
        "Mouse CCI (C57BL/6J)",
        "Sham; CCI+Saline/Vehicle; CCI+anti-TGF-β1 mAb (20, 30, 60 pg/mg)",
        "n = 6 (behavior & Schaffer–CA1 LTP)",
        "ICV anti-TGF-β1 mAb (BioLegend #521707) immediately or at 3 dpi; Schaffer–CA1 fEPSP recorded via MultiClamp 700B / Digidata 1550B / pClamp 10.7",
        "Y-maze 1, 7, 28/30 dpi; NOR 28 dpi; hippocampal LTP at 7 dpi",
        "Y-maze & NOR rescue; gliosis suppression; Schaffer collateral–CA1 fEPSP LTP restoration (Fig. 3A–E)",
        "In-house wet-lab (SUMC 2021-134)"
    ],
    [
        "4. Astrocyte-specific Tgfbr1 knockout in hippocampal CA3",
        "Mouse CCI (C57BL/6J)",
        "Sham; CCI; CCI+scr; CCI+TGF-βR1 KO",
        "n = 8 (behavior); n = 4–6 (IF)",
        "Hippocampal CA3 stereotaxic AAV-PHP.eB-GFAP-SaCas9-sgTgfbr1 (3 µL, 1×10^{12} vg/mL; AP +3.0, ML −2.8, DV −3.0 mm) 3 wk pre-CCI",
        "Y-maze & NOR at 7, 21, 28 dpi; IF at 7 & 28 dpi",
        "HA-tag & GFAP+TGF-βRI knockout verification; Y-maze alternation & entries; NOR trajectory heatmaps & preference index (Fig. 5A–E)",
        "In-house wet-lab (SUMC 2021-134)"
    ],
    [
        "5. In-house hippocampal single-nucleus RNA-seq atlas (GSE333879)",
        "Mouse CCI, ipsilateral hippocampus",
        "Sham; CCI+scr; CCI+TGF-βR1 KO",
        "35,047 QC-passed nuclei (8 mice pooled/group)",
        "Same as Cohort 4 (AAV-GFAP-SaCas9-sgTgfbr1), harvested at 30 dpi; Seurat (v4.1.3) CCA integration & SingleR annotation",
        "30 dpi",
        "7-cell-type atlas & SingleR validation; 774 pure astrocytes (Homeo-A, Homeo-B 14× expansion, Reactive); Monocle3; pySCENIC (17/20 regulons rescued); CellOracle in silico Tgfbr1 KO (Fig. 6A–G); Astrocyte→DG CellChat rescue (Fig. 7C–D); DG_C0–C3 granule neuron state rescue (Fig. 8D)",
        "GEO: GSE333879 (GSM9776429–31; SUMC 2021-134)"
    ],
    [
        "6. Multi-cohort mouse & human single-cell/bulk cross-validation",
        "Mouse TBI (FPI/CCI/concussion) & Human surgical/post-mortem TBI",
        "Control vs TBI across 7 independent public cohorts",
        "n = 22 mice (Hippo vs Cortex); n = 17 human snRNA donors; n = 94 human bulk hippocampal samples",
        "Donor/sample-level pseudobulk DESeq2, CellChat v2, MEF2C coupling & 16-gene adhesion scoring",
        "6 h to 28 dpi (mouse); acute & chronic (human)",
        "Human TGFB1 surge (Fig. 1E); MEF2C–adhesion coupling (Fig. 6H); mouse & human CellChat collapse & hippocampal CA3/DG selective vulnerability (Fig. 7A–B, E–F); mouse time-course anoikis & human neuronal apoptosis (Fig. 8C, E; SFig. S1–S5)",
        "GEO: GSE209552, GSE180862, GSE101901, GSE271769, GSE277487, GSE290150, GSE226207; Allen ADTBI"
    ],
    [
        "7. Astrocyte–neuron direct co-culture adhesion blockade & rescue",
        "HT-22 hippocampal neurons × C8-D1A astrocytes",
        "DMSO vehicle; ADH-1; ADH-1 + recombinant N-cadherin",
        "n = 3 independent biological replicates",
        "Cyclic CHAVC N-cadherin antagonist peptide ADH-1 (0.2 mg/mL) ± recombinant N-cadherin Fc chimera (10 µg/mL, 24 h)",
        "24 h co-culture",
        "Triple IF (p-Akt Ser473 / NeuN / DAPI) & quantitative p-Akt survival signaling rescue (Fig. 8A–B)",
        "In-house wet-lab"
    ],
]
make_table(
    "Table1_experimental_design.docx",
    "Table 1. Summary of animal cohorts, in-house single-nucleus transcriptomics, multi-cohort validations, and co-culture design.",
    "CCI, controlled cortical impact; ICV, intracerebroventricular; mAb, monoclonal antibody; NSS, neurological severity score; "
    "NOR, novel object recognition; LTP, long-term potentiation; fEPSP, field excitatory postsynaptic potential; snRNA-seq, single-nucleus RNA sequencing.",
    t1_head, t1_rows, widths=[2.4, 2.0, 3.0, 1.8, 3.2, 2.0, 3.8, 2.2]
)

# ========== Table 2 (Supplementary Table 1) — Antibodies & key reagents (All verified catalog numbers filled) ==========
t2_head = ["Antibody / Reagent", "Type / Host (Clone)", "Application", "Dilution / Dose", "Vendor", "Catalog No."]
t2_rows = [
    ["Anti-TGF-β1", "Rabbit pAb", "IF", "1:500", "Servicebio", "GB11179"],
    ["Anti-TGF-β1 neutralizing mAb", "Mouse mAb (Clone 19D8)", "ICV in vivo", "20 / 30 / 60 pg/mg", "BioLegend", "521707"],
    ["Anti-TGF-βRI (ALK5)", "Rabbit pAb", "WB / IF", "WB 1:1000; IF 1:200", "Affinity Biosciences", "AF5347"],
    ["Anti-NeuN", "Mouse mAb (Clone E4M5P)", "IF", "1:2000", "Cell Signaling Technology", "94403"],
    ["Anti-GFAP", "Mouse mAb (Clone GA5)", "IF", "1:2000", "Cell Signaling Technology", "3670"],
    ["Anti-Iba1", "Rabbit pAb", "IF", "1:500", "GeneTex", "GTX100042"],
    ["Anti-HA-Tag", "Rabbit mAb (Clone C29F4)", "IF", "1:1000", "Cell Signaling Technology", "3724"],
    ["Anti-β-actin", "Mouse mAb (Clone AC-15)", "WB", "1:5000", "Sigma-Aldrich", "A5441"],
    ["Anti-Phospho-Akt (Ser473)", "Rabbit mAb (Clone D9E)", "IF (Co-culture)", "1:300", "Cell Signaling Technology", "4060"],
    ["Goat anti-Rabbit IgG (H+L), Alexa Fluor 488", "Goat pAb (Secondary)", "IF", "1:500", "Invitrogen (Thermo Fisher)", "A-11034"],
    ["Goat anti-Rabbit IgG (H+L), Alexa Fluor 594", "Goat pAb (Secondary)", "IF", "1:500", "Invitrogen (Thermo Fisher)", "A-11037"],
    ["Goat anti-Mouse IgG (H+L), Alexa Fluor 488", "Goat pAb (Secondary)", "IF", "1:500", "Invitrogen (Thermo Fisher)", "A-11029"],
    ["Goat anti-Mouse IgG (H+L), Alexa Fluor 594", "Goat pAb (Secondary)", "IF", "1:500", "Invitrogen (Thermo Fisher)", "A-11032"],
    ["Antifade Mounting Medium with DAPI", "Nuclear counterstain / mountant", "IF", "Ready-to-use", "Beyotime", "P0131"],
    ["Recombinant mouse TGF-β1", "Recombinant protein", "ICV in vivo", "1, 10, 50 ng / 25 g", "Cell Signaling Technology", "5231"],
    ["ADH-1 (Exherin, cyclic CHAVC)", "N-cadherin antagonist peptide", "In vitro co-culture", "0.2 mg/mL (24 h)", "MedChemExpress", "HY-13541"],
    ["Recombinant Mouse N-Cadherin Fc Chimera", "Recombinant protein (Fc chimera)", "In vitro co-culture rescue", "10 µg/mL (24 h)", "R&D Systems", "6626-NC"],
    ["ViraTrap AAV Purification Maxiprep Kit", "Viral purification kit", "AAV-PHP.eB preparation", "Per manufacturer protocol", "Biomiga", "BW-V1369"],
    ["One-Step TUNEL Apoptosis Kit", "Fluorescence assay kit", "IF", "Per manufacturer protocol", "Beyotime", "C1090"],
]
make_table(
    "Table2_antibodies_reagents.docx",
    "Supplementary Table 1. Primary and secondary antibodies and key experimental reagents.",
    "IF, immunofluorescence; WB, Western blot; ICV, intracerebroventricular; mAb, monoclonal antibody; pAb, polyclonal antibody.",
    t2_head, t2_rows, widths=[3.6, 2.6, 2.2, 2.6, 3.2, 2.8]
)

# ========== Table 3 (Supplementary Table 2) — Single-cell & transcriptomic datasets ==========
t3_head = [
    "Accession / Cohort", "Species & Model", "Severity", "Brain Region", "Platform",
    "Samples (Control vs TBI)", "Role in This Study (Figures)"
]
t3_rows = [
    [
        "GSE333879 (In-house)", "Mouse CCI + AAV-GFAP-sgTgfbr1", "Moderate", "Ipsilateral hippocampus", "10x Chromium 3' v3 snRNA-seq",
        "Sham (11,378), CCI+scr (6,845), CCI+TGF-βR1 KO (16,824); 35,047 nuclei total (8 mice pooled/group)",
        "Atlas & SingleR validation (Fig. 6A); 774 pure astrocytes, Homeo-B 14× expansion, Monocle3, pySCENIC, CellOracle (Fig. 6B–G); Astrocyte→DG CellChat rescue (Fig. 7C–D); DG_C0–C3 granule neuron state rescue (Fig. 8D; SFig. S1, S3)"
    ],
    [
        "GSE101901", "Mouse TBI (fluid percussion)", "Mild–Moderate", "Hippocampus", "Drop-seq scRNA-seq",
        "Sham vs TBI (integrated with GSE180862: n = 9 Sham vs 13 TBI mice)",
        "Hippocampus-specific CA3 (−72.7%) & DG (−37.6%) neuronal adhesion collapse (Fig. 7E–F); panel sensitivity & ligand mRNA (SFig. S2, S4)"
    ],
    [
        "GSE180862", "Mouse FPI", "Mild–Moderate", "Hippocampus + cerebral cortex", "Drop-seq scRNA-seq",
        "Sham vs 24 h & 7 dpi (integrated n = 22 mice across hippocampus & cortex)",
        "Public mouse Astrocyte→Neuron CellChat adhesion collapse (Fig. 7A); hippocampal vs cortical subtype contrast (Fig. 7E; SFig. S2, S4)"
    ],
    [
        "GSE271769", "Mouse CCI", "Severe", "Cerebral cortex", "10x Chromium scRNA-seq",
        "NoTBI (n = 4) vs D1 (n = 4), D7 (n = 4), D28 (n = 4) biological replicates",
        "Multi-replicate time-course neuronal survival withdrawal & anoikis index progression (Fig. 8C); cross-severity module (Fig. 7F; SFig. S1, S4)"
    ],
    [
        "GSE290150", "Mouse CCI", "Moderate", "Cerebral cortex", "10x Chromium scRNA-seq",
        "Sham (n = 3) vs TBI (n = 3)",
        "Cross-severity neuronal adhesion module (Fig. 7F); Tgfbr1 mRNA distribution & glutamate panel (SFig. S1, S4)"
    ],
    [
        "GSE277487", "Mouse CCI", "Severe", "Cerebral cortex", "10x Chromium snRNA-seq",
        "Control (n = 3) vs TBI (6 h, 2 d, 4 d)",
        "Cross-severity neuronal adhesion module (Fig. 7F); Tgfbr1 mRNA & time-series neuronal DEGs (SFig. S1, S4)"
    ],
    [
        "GSE226207", "Mouse TBI", "Moderate", "Cerebral cortex", "10x Chromium snRNA-seq",
        "Control (n = 4) vs TBI (n = 4)",
        "Mouse pseudobulk Mef2c–synaptic adhesion correlation (Fig. 6H)"
    ],
    [
        "GSE209552", "Human surgical TBI", "Severe", "Cerebral cortex", "10x Chromium snRNA-seq",
        "Control (n = 5) vs acute surgical TBI cortex (n = 12; 17 donors total, 13,326 neuronal nuclei)",
        "Human glial TGFB1 surge (Fig. 1E); MEF2C–adhesion coupling (Fig. 6H); human Astrocyte→Neuron CellChat collapse (Fig. 7B) & donor-level adhesion collapse (P = 1.77×10^{−5}, Fig. 7F); human neuronal apoptosis shift (Fig. 8E; SFig. S1)"
    ],
    [
        "Allen ADTBI Cohort", "Human post-mortem TBI", "History of TBI with LOC", "Hippocampus", "Bulk RNA-seq",
        "Control (n = 50) vs TBI (n = 44) human hippocampal samples (n = 94)",
        "Human hippocampal cohort context & bulk dilution boundary analysis (Fig. 1E; Discussion)"
    ],
]
make_table(
    "Table3_public_datasets.docx",
    "Supplementary Table 2. In-house and public transcriptomic datasets analyzed in this study.",
    "CCI, controlled cortical impact; FPI, fluid percussion injury; LOC, loss of consciousness; snRNA-seq, single-nucleus RNA sequencing; scRNA-seq, single-cell RNA sequencing.",
    t3_head, t3_rows, widths=[2.2, 2.3, 1.6, 2.2, 2.2, 3.6, 4.5]
)

print("\nAll 3 Word tables generated and synced ->", OUT, "and", OUT_SYNC)
