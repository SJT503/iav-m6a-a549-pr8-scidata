# -*- coding: utf-8 -*-
"""14b_fig4_assoc_lncrna.py — 只读缓存秒级出图
NewFig: 表达-甲基化四象限关联(主图)
SuppFig3: lncRNA 差异 m6A 火山图(补充)
对齐 Fig2 风格 (_figstyle)。方向已在缓存里转为 PR8/MOCK。
"""
import sys; sys.path.insert(0, r"D:\IAV_m6A\SciData_descriptor\scripts")
import _figstyle as fs
import os, numpy as np, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fs.apply_style()

RESDIR=r"D:\IAV_m6A\SciData_descriptor\results"
FIGDIR=os.path.join(RESDIR,"figures")
d=np.load(os.path.join(RESDIR,"_assoc_lnc_cache.npz"),allow_pickle=True)
ge=d["ge_l2"]; m6=d["m6_l2"]; lnc=d["lnc_l2"]; lncp=d["lnc_p"]

# ============ Figure 4: expression-methylation association (4-quadrant) ============
fig,ax=plt.subplots(figsize=(5.6,5.2))
TH=1.0
both=(np.abs(ge)>=TH)&(np.abs(m6)>=TH)
rest=~both
ax.scatter(ge[rest],m6[rest],s=3,color=fs.PALETTE["grey_pt"],alpha=0.25,rasterized=True,linewidths=0)
# four quadrant colors
qa=both&(ge>0)&(m6>0)   # hyper-up
qb=both&(ge<0)&(m6>0)   # hypo-mRNA up-m6A
qc=both&(ge>0)&(m6<0)
qd=both&(ge<0)&(m6<0)   # both down
ax.scatter(ge[qa],m6[qa],s=6,color="#C44E52",alpha=0.6,rasterized=True,linewidths=0,label=f"Up expr & up m6A (n={int(qa.sum())})")
ax.scatter(ge[qd],m6[qd],s=6,color="#4C72B0",alpha=0.6,rasterized=True,linewidths=0,label=f"Down expr & down m6A (n={int(qd.sum())})")
ax.scatter(ge[qb],m6[qb],s=6,color="#55A868",alpha=0.6,rasterized=True,linewidths=0,label=f"Down expr & up m6A (n={int(qb.sum())})")
ax.scatter(ge[qc],m6[qc],s=6,color="#E1A33A",alpha=0.6,rasterized=True,linewidths=0,label=f"Up expr & down m6A (n={int(qc.sum())})")
ax.axhline(0,color="#999",lw=0.6); ax.axvline(0,color="#999",lw=0.6)
ax.axhline(TH,ls="--",c="grey",lw=0.5); ax.axhline(-TH,ls="--",c="grey",lw=0.5)
ax.axvline(TH,ls="--",c="grey",lw=0.5); ax.axvline(-TH,ls="--",c="grey",lw=0.5)
ax.set_xlabel("Expression log2 fold change (PR8 / MOCK)")
ax.set_ylabel("m6A log2 fold change (PR8 / MOCK)")
ax.set_title("Expression–m6A methylation association", fontweight="bold", fontsize=12, pad=8)
lim=max(np.percentile(np.abs(ge),99.5),np.percentile(np.abs(m6),99.5))
ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
ax.legend(fontsize=7,markerscale=1.6,loc="upper left",framealpha=0.9,frameon=True,edgecolor="none")
for sp in ("top","right"): ax.spines[sp].set_visible(False)
fig.savefig(os.path.join(FIGDIR,"Figure4_assoc.pdf"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(FIGDIR,"Figure4_assoc.png"),dpi=150,bbox_inches="tight")
plt.close(fig)
print("SAVED Figure4_assoc")

# ============ SuppFig3: lncRNA differential m6A volcano ============
x=lnc; y=-np.log10(np.clip(lncp,1e-300,None))
sig=(np.abs(x)>=1)&(lncp<0.05)
hiP=sig&(x>0); hiM=sig&(x<0)
fig2,ax2=plt.subplots(figsize=(5.4,4.8))
ax2.scatter(x[~sig],y[~sig],s=3,color=fs.PALETTE["grey_pt"],alpha=0.3,rasterized=True,linewidths=0)
ax2.scatter(x[hiP],y[hiP],s=5,color=fs.PALETTE["PR8"],alpha=0.6,rasterized=True,linewidths=0,label=f"Higher m6A in PR8 (n={int(hiP.sum())})")
ax2.scatter(x[hiM],y[hiM],s=5,color=fs.PALETTE["MOCK"],alpha=0.6,rasterized=True,linewidths=0,label=f"Higher m6A in MOCK (n={int(hiM.sum())})")
ax2.axvline(1,ls="--",c="grey",lw=0.6); ax2.axvline(-1,ls="--",c="grey",lw=0.6)
ax2.axhline(-math.log10(0.05),ls="--",c="grey",lw=0.6)
ax2.set_xlabel("lncRNA m6A log2 fold change (PR8 / MOCK)")
ax2.set_ylabel("$-\\log_{10}$ p-value")
ax2.set_title("Supplementary Figure 3  Differential m6A methylation of lncRNAs",
              loc="left",fontweight="bold",fontsize=11)
ax2.legend(fontsize=7.5,markerscale=2,loc="upper right")
for sp in ("top","right"): ax2.spines[sp].set_visible(False)
fig2.savefig(os.path.join(FIGDIR,"SuppFigure3_lncRNA.pdf"),dpi=300,bbox_inches="tight")
fig2.savefig(os.path.join(FIGDIR,"SuppFigure3_lncRNA.png"),dpi=150,bbox_inches="tight")
plt.close(fig2)
print("SAVED SuppFigure3_lncRNA")

# QC
from PIL import Image
for n in ["Figure4_assoc","SuppFigure3_lncRNA"]:
    im=np.asarray(Image.open(os.path.join(FIGDIR,n+".png")).convert("RGB"),dtype=float)
    nw=(im.min(axis=2)<245); h,w,_=im.shape
    b=max(2,int(min(h,w)*0.01))
    edge=np.concatenate([nw[:b,:].ravel(),nw[-b:,:].ravel(),nw[:,:b].ravel(),nw[:,-b:].ravel()]).mean()
    print(f"QC {n}: {w}x{h} nonwhite={nw.mean():.3f} edge_ink={edge:.3f}")
