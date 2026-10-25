# -*- coding: utf-8 -*-
"""虚拟敲除 step2 —— 基于 GRN 的 FOXG1 在硅敲除 (Ridge 回归 GRN + 敲除传播)
方法等价于 CellOracle 核心算法; 仅用 sklearn/scipy, 不依赖 celloracle 包。
真实数据: 整合小鼠星形胶质 19808 细胞 (GSE101901+GSE180862)。"""
import numpy as np, pandas as pd, scipy.io, scipy.sparse, os
from sklearn.linear_model import Ridge
from scipy.stats import mannwhitneyu

D = "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/vko"
P = "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"

genes = [l.strip() for l in open(os.path.join(D,"genes.txt"),encoding="utf-8")]
cells = [l.strip() for l in open(os.path.join(D,"cells.txt"),encoding="utf-8")]
X = scipy.io.mmread(os.path.join(D,"data.mtx")).tocsr().T.toarray()  # cells x genes (log-norm)
print("表达矩阵:", X.shape, "(细胞 x 基因)")
gidx = {g:i for i,g in enumerate(genes)}

tf_mm = [l.strip() for l in open(os.path.join(P,"scenic_input/databases/allTFs_mm.txt"),encoding="utf-8")]
tfs = [g for g in genes if g in set(tf_mm)]
print("可用 TF:", len(tfs))
assert "Foxg1" in genes, "Foxg1 不在矩阵中"

# 标准化
mu, sd = X.mean(0), X.std(0); sd[sd==0]=1
Z = (X - mu)/sd
Xtf = Z[:, [gidx[t] for t in tfs]]            # 预测变量: TF (标准化)
# Ridge GRN: 每个基因 ~ 全部 TF  (多输出一次解出)
ridge = Ridge(alpha=1.0)
ridge.fit(Xtf, Z)                              # coef_: nGene x nTF
coef = ridge.coef_
tf_pos = {t:i for i,t in enumerate(tfs)}
print("GRN 拟合完成: coef", coef.shape)

# FOXG1 虚拟敲除: FOXG1 表达置 0
fi = gidx["Foxg1"]
foxg1_raw = X[:, fi]
ko_std = (0 - mu[fi])/sd[fi]                   # 敲除后标准化值
d_foxg1 = ko_std - Z[:, fi]                    # 每细胞 ΔFOXG1 (标准化)
# 一步敲除位移: Δgene = coef[:,FOXG1] * ΔFOXG1
shift1 = np.outer(d_foxg1, coef[:, tf_pos["Foxg1"]])   # cells x genes
# 一步传播: 受扰动的 TF 再影响下游
dtf = shift1[:, [gidx[t] for t in tfs]]
shift2 = shift1 + dtf @ coef.T
mean_shift = shift2.mean(0)                    # 每基因平均预测位移
res = pd.DataFrame({"gene":genes, "ko_shift":mean_shift,
                    "foxg1_coef":coef[:, tf_pos["Foxg1"]]}).sort_values("ko_shift")
res.to_csv(os.path.join(D,"foxg1_ko_shift.csv"), index=False)

# 黏附程序基因
ADH = ["Nrxn1","Nrxn3","Cadm1","Cadm2","Cntn1","Cdh2"]
adh = [g for g in ADH if g in gidx]
print("\n=== 黏附程序基因 FOXG1 敲除预测位移 ===")
for g in adh:
    print(f"  {g:8s} shift={mean_shift[gidx[g]]:+.4f}  FOXG1_coef={coef[gidx[g],tf_pos['Foxg1']]:+.4f}")
adh_sh = np.array([mean_shift[gidx[g]] for g in adh])
bg = mean_shift
u,p = mannwhitneyu(adh_sh, bg, alternative="less")
print(f"\n黏附基因 vs 全基因背景: 黏附均值={adh_sh.mean():+.4f}  背景中位数={np.median(bg):+.4f}")
print(f"Mann-Whitney (黏附更低): p={p:.2e}")
print(f"黏附 6 基因下调比例: {(adh_sh<0).sum()}/6")

# 全量: 下调/上调 top
print("\n=== 全量 top10 下调 / top10 上调 ===")
print("下调:", list(res.head(10)["gene"]))
print("上调:", list(res.tail(10)["gene"]))

np.save(os.path.join(D,"mean_shift.npy"), mean_shift)
pd.DataFrame({"gene":genes}).to_csv(os.path.join(D,"genes_order.csv"),index=False)
# 每细胞扰动强度(位移L2范数), 按分组
meta = pd.read_csv(os.path.join(D,"meta.csv"))
cell_pert = np.sqrt((shift2**2).sum(1))
pd.DataFrame({"cell":cells,"condition":meta["condition"].values,
              "pert":cell_pert}).to_csv(os.path.join(D,"cell_pert.csv"),index=False)
print("\nDONE. 输出 ->", D)
