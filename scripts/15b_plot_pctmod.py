# -*- coding: utf-8 -*-
"""15b_plot_pctmod.py — 只读 _pctmod_cache.npz 秒级出图
(1) %Modified 关联图：x=表达 log2FC(PR8/MOCK), y=delta %Modified(PR8-MOCK)
    展示真正解耦——表达变但甲基化密度不变 / 反之。
(2) %Modified 差异：6 样本 t-test on %m6A，PR8 vs MOCK 火山图（真甲基化密度差异）。
对齐 Fig2 风格。venue 红线：描述性，不写机制因果。
"""
import sys; sys.path.insert(0, r"D:\IAV_m6A\SciData_descriptor\scripts")
import _figstyle as fs
import os, numpy as np, math
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fs.apply_style()

RESDIR=r"D:\IAV_m6A\SciData_descriptor\results"; FIGDIR=os.path.join(RESDIR,"figures")
d=np.load(os.path.join(RESDIR,"_pctmod_cache.npz"),allow_pickle=True)
ge=d["ge_l2"]; dpm=d["dpm"]; samp=d["samp"]; mockpm=d["mockpm"]; pr8pm=d["pr8pm"]

# ===== (1) decoupling scatter: expr log2FC vs delta %Modified =====
fig,ax=plt.subplots(figsize=(5.8,5.2))
TH_E=1.0; TH_M=0.05
# color: methylation-density changed (|dpm|>=TH_M) vs not
mch=np.abs(dpm)>=TH_M
ax.scatter(ge[~mch],dpm[~mch],s=3,color=fs.PALETTE["grey_pt"],alpha=0.25,rasterized=True,linewidths=0)
ax.scatter(ge[mch&(dpm>0)],dpm[mch&(dpm>0)],s=5,color=fs.PALETTE["PR8"],alpha=0.55,rasterized=True,linewidths=0,
           label=f"Higher %m6A in PR8 (n={int((mch&(dpm>0)).sum())})")
ax.scatter(ge[mch&(dpm<0)],dpm[mch&(dpm<0)],s=5,color=fs.PALETTE["MOCK"],alpha=0.55,rasterized=True,linewidths=0,
           label=f"Higher %m6A in MOCK (n={int((mch&(dpm<0)).sum())})")
ax.axhline(0,color="#999",lw=0.6); ax.axvline(0,color="#999",lw=0.6)
ax.axhline(TH_M,ls="--",c="grey",lw=0.5); ax.axhline(-TH_M,ls="--",c="grey",lw=0.5)
ax.axvline(TH_E,ls="--",c="grey",lw=0.5); ax.axvline(-TH_E,ls="--",c="grey",lw=0.5)
ax.set_xlabel("Expression log2 fold change (PR8 / MOCK)")
ax.set_ylabel("$\\Delta$ m6A methylation fraction (PR8 $-$ MOCK)")
ax.set_title("Methylation stoichiometry versus expression change", fontweight="bold", fontsize=12, pad=8)
xlim=np.percentile(np.abs(ge),99.5); ylim=np.percentile(np.abs(dpm),99.5)
ax.set_xlim(-xlim,xlim); ax.set_ylim(-ylim,ylim)
ax.legend(fontsize=7.5,markerscale=2,loc="upper right",frameon=True,edgecolor="none",framealpha=0.9)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
fig.savefig(os.path.join(FIGDIR,"Figure4_methylation_stoichiometry.pdf"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(FIGDIR,"Figure4_methylation_stoichiometry.png"),dpi=150,bbox_inches="tight")
plt.close(fig); print("SAVED Figure4 (stoichiometry decoupling)")

# ===== (2) %Modified differential volcano: t-test on 6 samples =====
mock=samp[:,:3]; pr8=samp[:,3:]
# paired? no — independent groups n=3. Welch t-test per transcript
with np.errstate(all="ignore"):
    t,p=stats.ttest_ind(pr8,mock,axis=1,equal_var=False)
delta=pr8.mean(axis=1)-mock.mean(axis=1)   # delta %Modified
valid=np.isfinite(p)
x=delta; y=-np.log10(np.clip(p,1e-300,None))
sig=valid&(np.abs(delta)>=0.05)&(p<0.05)
hiP=sig&(delta>0); hiM=sig&(delta<0)
fig2,ax2=plt.subplots(figsize=(5.6,4.9))
ax2.scatter(x[valid&~sig],y[valid&~sig],s=3,color=fs.PALETTE["grey_pt"],alpha=0.3,rasterized=True,linewidths=0)
ax2.scatter(x[hiP],y[hiP],s=6,color=fs.PALETTE["PR8"],alpha=0.6,rasterized=True,linewidths=0,label=f"Higher %m6A in PR8 (n={int(hiP.sum())})")
ax2.scatter(x[hiM],y[hiM],s=6,color=fs.PALETTE["MOCK"],alpha=0.6,rasterized=True,linewidths=0,label=f"Higher %m6A in MOCK (n={int(hiM.sum())})")
ax2.axvline(0.05,ls="--",c="grey",lw=0.6); ax2.axvline(-0.05,ls="--",c="grey",lw=0.6)
ax2.axhline(-math.log10(0.05),ls="--",c="grey",lw=0.6)
ax2.set_xlabel("$\\Delta$ m6A methylation fraction (PR8 $-$ MOCK)")
ax2.set_ylabel("$-\\log_{10}$ p-value (Welch t-test)")
ax2.set_title("Supplementary Figure 4  Differential m6A methylation fraction (%modified)", loc="left", fontweight="bold", fontsize=11, pad=6)
ax2.legend(fontsize=7.5,markerscale=2,loc="upper center")
for sp in ("top","right"): ax2.spines[sp].set_visible(False)
fig2.savefig(os.path.join(FIGDIR,"SuppFigure4_pctmod_diff.pdf"),dpi=300,bbox_inches="tight")
fig2.savefig(os.path.join(FIGDIR,"SuppFigure4_pctmod_diff.png"),dpi=150,bbox_inches="tight")
plt.close(fig2); print("SAVED SuppFig4 (pctmod differential)")

# stats dump for text/legend
with open(os.path.join(RESDIR,"_pctmod_diff_stats.txt"),"w",encoding="utf-8") as f:
    f.write(f"transcripts tested={int(valid.sum())}\n")
    f.write(f"diff %Modified (|delta|>=0.05 & p<0.05): total={int(sig.sum())} higher_PR8={int(hiP.sum())} higher_MOCK={int(hiM.sum())}\n")
    f.write(f"delta %Mod range=[{delta[valid].min():.3f},{delta[valid].max():.3f}]\n")
# QC
from PIL import Image
for n in ["SuppFigure4_pctmod_assoc","SuppFigure5_pctmod_diff"]:
    im=np.asarray(Image.open(os.path.join(FIGDIR,n+".png")).convert("RGB"),dtype=float)
    nw=(im.min(axis=2)<245); h,w,_=im.shape; b=max(2,int(min(h,w)*0.01))
    edge=np.concatenate([nw[:b,:].ravel(),nw[-b:,:].ravel(),nw[:,:b].ravel(),nw[:,-b:].ravel()]).mean()
    print(f"QC {n}: {w}x{h} nonwhite={nw.mean():.3f} edge_ink={edge:.3f}")
