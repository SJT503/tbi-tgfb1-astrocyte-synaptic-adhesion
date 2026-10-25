# Glia-Derived TGF-β1 Drives Post-Traumatic Cognitive Impairment by Reprogramming Astrocyte–Neuron Synaptic Adhesion

This repository contains the complete custom R and Python analysis, figure-generation, and table-generation scripts (101 scripts) and processed summary tables (94 tables) for:

> **Sheng J†, Liu X†, Zhu L, Ren Z, Lin R, Sun Q, Li K*, Chen W*.** *"Glia-derived TGF-β1 drives post-traumatic cognitive impairment by reprogramming astrocyte–neuron synaptic adhesion"* (*Journal of Neuroinflammation*, 2026).

---

## Overview

Combining mouse controlled cortical impact (CCI), intracerebroventricular (ICV) gain- and loss-of-function dosing, hippocampal CA3 astrocyte-specific CRISPR/SaCas9 *Tgfbr1* deletion (`AAV-PHP.eB-GFAP-SaCas9-sgRNA`), Schaffer collateral–CA1 long-term potentiation (LTP) electrophysiology, in-house single-nucleus RNA sequencing of 35,047 hippocampal nuclei (**GSE333879**), cross-species validation across seven public mouse and human TBI cohorts (including **GSE209552** and **Allen ADTBI**), and direct HT-22 × C8-D1A astrocyte–neuron co-culture N-cadherin blockade (`ADH-1`) and recombinant N-cadherin (`rN-cadherin`) rescue, this study resolves how glia-derived TGF-β1 signals through astrocytic TGF-βRI (ALK5) to suppress a synaptogenic Homeostatic-B (*Apoe*⁺ / *Cd81*⁺ / *Nrxn1*⁺) astrocyte program (`MEF2C` / `Gpc5`), dismantle trans-synaptic `NRXN` and `N-cadherin (CDH2)` adhesion onto hippocampal CA3 and dentate gyrus (DG) neurons, impair CA1 LTP, and trigger neuronal anoikis.

---

## Data Availability

- **In-house mouse hippocampal snRNA-seq (Sham, CCI+scr, CCI+TGF-βR1 KO; 35,047 nuclei)**: NCBI GEO [**GSE333879**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333879) (`GSM9776429`–`GSM9776431`)
- **Public human & mouse TBI cohorts re-analyzed**:
  - **GSE209552**: Human surgical TBI cortex snRNA-seq (*n* = 17 donors; 13,326 neuronal nuclei)
  - **Allen Brain Map ADTBI**: Human post-mortem hippocampal bulk RNA-seq (*n* = 94 donors: 50 controls vs 44 TBI)
  - **GSE101901 + GSE180862**: Paired mouse hippocampal and cortical TBI scRNA-seq (*n* = 22 mice: 9 Sham vs 13 TBI)
  - **GSE271769**: Longitudinal mouse TBI scRNA-seq across Sham, 1, 7, and 28 dpi (*n* = 16 mice)
  - **GSE277487**: Severe mouse TBI snRNA-seq across 6 h, 2 d, and 4 d
  - **GSE290150**: Moderate mouse CCI scRNA-seq (*n* = 6 mice)
  - **GSE226207**: Moderate mouse TBI snRNA-seq (*n* = 8 mice for donor-level *Mef2c*–adhesion coupling)

---

## Repository Structure

| Directory | Contents |
|---|---|
| `scripts/analysis/` | **87 R & Python scripts** covering QC, ratio-based doublet removal (`774` pure hippocampal astrocytes; `3,522` *Prox1*⁺ DG granule neurons), Monocle 3 pseudotime trajectories, sequencing-depth-independent pySCENIC AUCell regulon scoring, CellOracle GRN in silico *Tgfbr1* perturbation (`2,005` genes), CellChat v2 astrocyte→neuron synaptic-adhesion inference, 16-gene synaptic-adhesion module & sensitivity panels, cross-species donor-level *Mef2c/MEF2C* coupling, and neuronal anoikis / negative-gating analyses (`100`–`106`). |
| `scripts/figures/` | **14 R & Python scripts** generating Main Figures 1–9 (`Fig 1_cleaned`–`Fig 9_cleaned`), Supplementary Figures S1–S5 (`Fig S1_cleaned`–`Fig S5_cleaned`), and Main/Supplementary Tables (`Table 1`, `Supplementary Tables S1–S2`) in a unified colorblind-safe NPG + Okabe–Ito palette. |
| `tables/` | **94 processed CSV summary tables** containing per-cell, per-cluster, and per-donor/per-mouse statistical readouts backing Figures 1–8 and Supplementary Figures S1–S5. |

---

## Software Requirements

- **R (>= 4.3)**: `Seurat` (v4.1.3 / v5), `SingleR`, `monocle3`, `CellChat` (v2), `DESeq2`, `ggplot2`, `patchwork`, `ComplexHeatmap`
- **Python (>= 3.10)**: `scanpy`, `pyscenic` (AUCell), `celloracle`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`, `python-docx`

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
