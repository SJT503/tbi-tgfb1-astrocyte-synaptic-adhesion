#!/usr/bin/env python
# ADTBI human hippocampus TBI vs control - adhesion ligand bulk RNA-seq
# Data: Allen Aging Dementia TBI project (107 donors, bulk RNA-seq, 4 brain regions)
# Subset: HIP only; TBI = ever_tbi_w_loc='Y'; control = 'N'
# Test: log2(FPKM+1) Mann-Whitney U per gene
import json, csv, sys
import numpy as np
from scipy import stats

BASE = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\JNI\analysis\ADTBI_human_hippo"
PANEL = ["NRXN1","NRXN2","NRXN3","CADM1","CADM2","CADM3","CNTN1","CNTN2","CNTN4","NEGR1","NCAM1","NCAM2"]

# --- 1. donor metadata (TBI Y/N) ---
print("[1] donor metadata")
d = json.load(open(f"{BASE}/ApiTbiDonorDetail.json", encoding='utf-8'))
donor_tbi = {}
for r in d['msg']:
    donor_tbi[int(r['donor_id'])] = r['ever_tbi_w_loc']
n_y = sum(1 for v in donor_tbi.values() if v=='Y')
n_n = sum(1 for v in donor_tbi.values() if v=='N')
print(f"  donors: {len(donor_tbi)} | TBI(Y): {n_y} | control(N): {n_n}")

# --- 2. sample metadata: HIP only ---
print("\n[2] HIP samples")
hip_samples = []   # list of (rnaseq_profile_id, donor_id, tbi_label)
with open(f"{BASE}/columns-samples.csv", encoding='utf-8') as f:
    rd = csv.DictReader(f)
    for r in rd:
        if r['structure_acronym'] == 'HIP':
            did = int(r['donor_id'])
            tbi = donor_tbi.get(did, '?')
            hip_samples.append((r['rnaseq_profile_id'], did, tbi))
print(f"  HIP samples: {len(hip_samples)}")
hip_y = sum(1 for s in hip_samples if s[2]=='Y')
hip_n = sum(1 for s in hip_samples if s[2]=='N')
print(f"  HIP TBI: {hip_y} | HIP control: {hip_n}")

# --- 3. find PANEL genes in rows-genes.csv ---
print("\n[3] gene index")
gene_idx = {}  # gene_symbol -> row_index (1-based for fpkm matrix)
with open(f"{BASE}/rows-genes.csv", encoding='utf-8') as f:
    rd = csv.DictReader(f)
    # header columns may include: 'gene_id','chromosome','gene_entrez_id','gene_symbol','gene_name'
    print(f"  header: {rd.fieldnames}")
    for i, r in enumerate(rd):
        sym = r.get('gene_symbol','').strip().upper()
        if sym in PANEL:
            gene_idx[sym] = i  # 0-based; matrix is 0-based row positionally
print(f"  PANEL genes found: {gene_idx}")

# --- 4. read FPKM matrix, extract PANEL gene rows, HIP sample columns ---
print("\n[4] read FPKM matrix (normalized)")
sample_ids = []   # column header sample IDs
panel_rows = {g: None for g in gene_idx}  # gene -> np array of FPKM (per column)
with open(f"{BASE}/fpkm_table_normalized.csv", encoding='utf-8') as f:
    header = f.readline().strip().split(',')
    sample_ids = header[1:]  # first col is gene_id
    print(f"  total columns: {len(sample_ids)}; sample IDs sample: {sample_ids[:3]}")
    # build positional map (row_index -> gene if in panel)
    target_rows = {v: k for k, v in gene_idx.items()}
    found = {}
    for i, line in enumerate(f):
        if i in target_rows:
            parts = line.strip().split(',')
            vals = np.array(parts[1:], dtype=float)
            found[target_rows[i]] = vals
        if len(found) == len(gene_idx):
            print(f"  all {len(found)} panel genes captured at row index up to {i}")
            break
    panel_rows = found

# --- 5. find HIP column positions ---
hip_col_pos = []
hip_labels = []
# strip quotes from sample IDs
sample_to_col = {sid.strip('"'): pos for pos, sid in enumerate(sample_ids)}
for prof_id, did, tbi in hip_samples:
    pos = sample_to_col.get(str(prof_id))
    if pos is None: continue
    hip_col_pos.append(pos)
    hip_labels.append(tbi)
hip_col_pos = np.array(hip_col_pos, dtype=int)
hip_labels = np.array(hip_labels)
mask_y = hip_labels == 'Y'
mask_n = hip_labels == 'N'
print(f"  HIP matched to columns: {len(hip_col_pos)} | TBI {mask_y.sum()} | control {mask_n.sum()}")

# --- 6. per-gene Mann-Whitney U on log2(FPKM+1), HIP TBI vs control ---
print("\n[6] DEG per gene (HIP TBI vs control)")
print(f"{'gene':<8}{'n_TBI':>6}{'n_ctrl':>8}{'med_TBI':>12}{'med_ctrl':>12}{'pct_change':>12}{'mwU_p':>10}")
results = []
for g in PANEL:
    if g not in panel_rows:
        results.append({'gene': g, 'note': 'not_in_matrix'})
        continue
    fpkm = panel_rows[g]
    log_fpkm = np.log2(fpkm + 1)
    vals_tbi = log_fpkm[hip_col_pos[mask_y]]
    vals_ctrl = log_fpkm[hip_col_pos[mask_n]]
    med_tbi = np.median(2**vals_tbi - 1)  # back to FPKM scale
    med_ctrl = np.median(2**vals_ctrl - 1)
    if med_ctrl < 0.001:
        pct = float('nan')
    else:
        pct = 100 * (med_tbi - med_ctrl) / med_ctrl
    try:
        u, p = stats.mannwhitneyu(vals_tbi, vals_ctrl, alternative='two-sided')
    except ValueError:
        p = float('nan')
    print(f"{g:<8}{len(vals_tbi):>6}{len(vals_ctrl):>8}{med_tbi:>12.3f}{med_ctrl:>12.3f}{pct:>11.1f}%{p:>10.3g}")
    results.append({'gene': g, 'n_TBI': len(vals_tbi), 'n_ctrl': len(vals_ctrl),
                    'median_TBI_fpkm': med_tbi, 'median_ctrl_fpkm': med_ctrl,
                    'pct_change': pct, 'mw_p': p})

# --- 7. save ---
import json as J
with open(f"{BASE}/ADTBI_hippo_adhesion_DEG_results.json", 'w', encoding='utf-8') as f:
    J.dump({'panel': PANEL, 'n_TBI': int(mask_y.sum()), 'n_ctrl': int(mask_n.sum()),
            'results': results}, f, indent=2, default=str)
import csv as C
with open(f"{BASE}/ADTBI_hippo_adhesion_DEG_results.csv", 'w', newline='', encoding='utf-8') as f:
    w = C.DictWriter(f, fieldnames=['gene','n_TBI','n_ctrl','median_TBI_fpkm','median_ctrl_fpkm','pct_change','mw_p'])
    w.writeheader()
    for r in results:
        if 'mw_p' in r:
            w.writerow(r)

print(f"\nDONE. Saved to {BASE}/ADTBI_hippo_adhesion_DEG_results.{{csv,json}}")
