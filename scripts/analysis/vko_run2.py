# -*- coding: utf-8 -*-
"""虚拟敲除 step2b —— 海马星形胶质子集 + 去稀释 GRN
区域分层(海马 vs 皮层)检 FOXG1-黏附相关; 海马子集上跑去稀释 GRN 虚拟敲除。"""
import numpy as np, pandas as pd, scipy.io, os
from sklearn.linear_model import Ridge
from scipy.stats import spearmanr, mannwhitneyu

D = "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/vko"
P = "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
genes = [l.strip() for l in open(os.path.join(D,"genes.txt"),encoding="utf-8")]
X = scipy.io.mmread(os.path.join(D,"data.mtx")).tocsr().T.toarray()
meta = pd.read_csv(os.path.join(D,"meta.csv"))
gidx = {g:i for i,g in enumerate(genes)}
ADH = [g for g in ["Nrxn1","Nrxn3","Cadm1","Cadm2","Cntn1","Cdh2"] if g in gidx]
adh_i = [gidx[g] for g in ADH]

# 区域分层检相关性
print("=== FOXG1 与黏附模块 Spearman 相关(区域分层) ===")
fi = gidx["Foxg1"]
for reg in ["hippocampus","cortex"]:
    m = (meta["region"]==reg).values
    if m.sum()<50: print(f"  {reg}: 细胞太少({m.sum()})"); continue
    adh_mod = X[np.ix_(m, adh_i)].mean(1)
    rho,p = spearmanr(X[m, fi], adh_mod)
    print(f"  {reg}: n={m.sum()}  FOXG1-黏附 rho={rho:+.3f} (p={p:.1e})")

# 选海马子集
hip = (meta["region"]=="hippocampus").values
print(f"\n海马星形胶质子集: n={hip.sum()}")
Xh = X[hip]
# 去稀释: 预测变量取 top-150 高变 TF
tf_mm = set(l.strip() for l in open(os.path.join(P,"scenic_input/databases/allTFs_mm.txt"),encoding="utf-8"))
tfs_all = [g for g in genes if g in tf_mm]
tf_var = {t: Xh[:,gidx[t]].var() for t in tfs_all}
tfs = sorted(tf_var, key=tf_var.get, reverse=True)[:150]
if "Foxg1" not in tfs: tfs = tfs[:149]+["Foxg1"]
print(f"GRN 预测变量: top-{len(tfs)} 高变 TF (去稀释)")

mu,sd = Xh.mean(0), Xh.std(0); sd[sd==0]=1
Z = (Xh-mu)/sd
Xtf = Z[:, [gidx[t] for t in tfs]]
ridge = Ridge(alpha=1.0); ridge.fit(Xtf, Z)
coef = ridge.coef_                                  # nGene x nTF
fpos = tfs.index("Foxg1")
print("GRN 拟合完成。")

# FOXG1 敲除 + 一步传播
ko_std = (0-mu[fi])/sd[fi]
d_foxg1 = ko_std - Z[:, fi]
shift1 = np.outer(d_foxg1, coef[:, fpos])
dtf = shift1[:, [gidx[t] for t in tfs]]
shift2 = shift1 + dtf @ coef.T
ms = shift2.mean(0)

print("\n=== 海马: FOXG1 敲除黏附程序基因预测位移 ===")
for g in ADH:
    print(f"  {g:8s} shift={ms[gidx[g]]:+.4f}  FOXG1_coef={coef[gidx[g],fpos]:+.4f}")
adh_sh = np.array([ms[i] for i in adh_i])
u,pp = mannwhitneyu(adh_sh, ms, alternative="less")
print(f"\n黏附均值={adh_sh.mean():+.4f}  背景中位数={np.median(ms):+.4f}")
print(f"下调比例 {(adh_sh<0).sum()}/{len(ADH)}  Mann-Whitney(黏附更低) p={pp:.2e}")
pd.DataFrame({"gene":genes,"ko_shift":ms,"foxg1_coef":coef[:,fpos]}
            ).sort_values("ko_shift").to_csv(os.path.join(D,"foxg1_ko_hippo.csv"),index=False)
print("\nDONE")
