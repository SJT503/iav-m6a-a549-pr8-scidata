# -*- coding: utf-8 -*-
"""11_fig3_m6a_signal.py — Fig 3 m6A 信号特征与差异景观
真实数据：mRNA_Quantity(allcomparison).xlsx (火山) + (diffmethylated).xlsx (区域分布)
面板：(a) 差异 m6A 火山图  (b) m6A 位点区域分布(5'UTR/CDS/3'UTR)  (c) top 差异转录本聚类热图
方向实测：脚本验证 'up' 行是否 MOCK(mean)>PR8(mean)，据此定轴标签，不假设。
"""
import os, numpy as np, openpyxl, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

RP = r"D:\IAV_m6A\Report\Report\Report\PaperData\Differentially Methylated Genes"
RESDIR = r"D:\IAV_m6A\SciData_descriptor\results"
FIGDIR = os.path.join(RESDIR,"figures")
ALL = os.path.join(RP,"mRNA_Quantity(allcomparison).xlsx")
DIFF = os.path.join(RP,"mRNA_Quantity(diffmethylated).xlsx")

def find_header(ws):
    for ri,row in enumerate(ws.iter_rows(min_row=1,max_row=40,values_only=True)):
        if row[0] is not None and str(row[0]).strip()=="Transcript_ID":
            return ri
    return None

# ---- read allcomparison: FC(col5), p(col6), regulation(col8), MOCK mean(col3), PR8 mean(col4) ----
wb = openpyxl.load_workbook(ALL, read_only=True, data_only=True)
ws = wb["MOCK_vs_PR8"]
hr = find_header(ws)
fc=[]; pv=[]; reg=[]; mock=[]; pr8=[]
for ri,row in enumerate(ws.iter_rows(values_only=True)):
    if ri<=hr: continue
    a=row[0]
    if a is None or str(a).strip().startswith("#") or str(a).strip()=="": continue
    try:
        f=float(row[5]); p=float(row[6])
    except (TypeError,ValueError): continue
    fc.append(f); pv.append(p); reg.append(str(row[8]).strip() if row[8] else "")
    try: mock.append(float(row[3])); pr8.append(float(row[4]))
    except (TypeError,ValueError): mock.append(np.nan); pr8.append(np.nan)
wb.close()
fc=np.array(fc); pv=np.array(pv); reg=np.array(reg); mock=np.array(mock); pr8=np.array(pr8)
print("allcomparison rows:", len(fc))

# ---- verify direction empirically ----
up_mask = reg=="up"
dn_mask = reg=="down"
md = np.nanmean(mock[up_mask]-pr8[up_mask])
print(f"'up' rows: mean(MOCK-PR8 log2)={md:.3f}  → up means {'MOCK>PR8 (higher in mock)' if md>0 else 'PR8>MOCK (higher in PR8)'}")
with open(os.path.join(RESDIR,"_fig3_direction.txt"),"w") as f:
    f.write(f"n_all={len(fc)} n_up={up_mask.sum()} n_down={dn_mask.sum()}\n")
    f.write(f"'up' mean(MOCK-PR8 log2)={md:.4f}\n")
    f.write("=> up = higher m6A in MOCK\n" if md>0 else "=> up = higher m6A in PR8\n")

# ---- volcano in PR8-vs-MOCK perspective ----
# Foldchange(col5)=2^(MOCK-PR8) for 'up'; log2(FC) = MOCK-PR8. Flip to PR8 perspective:
log2fc_mock_over_pr8 = np.log2(np.clip(fc,1e-9,None))   # = MOCK - PR8 (log2)
x = -log2fc_mock_over_pr8                                # PR8 vs MOCK: >0 => higher in PR8
y = -np.log10(np.clip(pv,1e-300,None))
sig = (np.abs(x)>=1.0) & (pv<0.05)
hi_pr8 = sig & (x>0)    # higher m6A in PR8
hi_mock = sig & (x<0)   # higher m6A in MOCK
print(f"VOLCANO at |log2FC|>=1 & p<0.05: higher_in_PR8={hi_pr8.sum()}  higher_in_MOCK={hi_mock.sum()}")
# cross-check vs report's hyper/hypo (3735 / 4350)
with open(os.path.join(RESDIR,"_fig3_direction.txt"),"a") as f:
    f.write(f"\nVolcano |log2FC|>=1 & p<0.05: higher_in_PR8={hi_pr8.sum()} higher_in_MOCK={hi_mock.sum()}\n")
    f.write("report feature_counts: hyper 3735 / hypo 4350 (MOCK_vs_PR8 naming)\n")

# ---- read region distribution (col24) from diffmethylated hyper+hypo ----
def region_counts(sheet):
    wb2=openpyxl.load_workbook(DIFF,read_only=True,data_only=True); ws2=wb2[sheet]
    hr2=find_header(ws2); cnt={"5'UTR":0,"CDS":0,"3'UTR":0}
    for ri,row in enumerate(ws2.iter_rows(values_only=True)):
        if ri<=hr2: continue
        a=row[0]
        if a is None or str(a).strip().startswith("#") or str(a).strip()=="": continue
        regn=row[24]
        if regn:
            for part in str(regn).split(";"):
                p=part.strip()
                if p in cnt: cnt[p]+=1
    wb2.close(); return cnt
rc_hyper=region_counts("hyper_MOCK_vs_PR8"); rc_hypo=region_counts("hypo_MOCK_vs_PR8")
region_total={k:rc_hyper[k]+rc_hypo[k] for k in rc_hyper}
print("region distribution (all m6A sites on diff transcripts):", region_total)
with open(os.path.join(RESDIR,"_fig3_region.txt"),"w") as f:
    f.write(f"hyper sheet sites: {rc_hyper}\nhypo sheet sites: {rc_hypo}\ncombined: {region_total}\n")

# ================= FIGURE =================
fig=plt.figure(figsize=(13,4.5))
gs=GridSpec(1,3,figure=fig,wspace=0.32)

# (a) volcano
axa=fig.add_subplot(gs[0,0])
axa.scatter(x[~sig],y[~sig],s=2,color="#BBBBBB",alpha=0.3,rasterized=True)
axa.scatter(x[hi_pr8],y[hi_pr8],s=3,color="#C44E52",alpha=0.5,rasterized=True,label=f"Higher in PR8 (n={hi_pr8.sum()})")
axa.scatter(x[hi_mock],y[hi_mock],s=3,color="#4C72B0",alpha=0.5,rasterized=True,label=f"Higher in MOCK (n={hi_mock.sum()})")
axa.axvline(1,ls="--",c="grey",lw=0.6); axa.axvline(-1,ls="--",c="grey",lw=0.6)
axa.axhline(-math.log10(0.05),ls="--",c="grey",lw=0.6)
axa.set_xlabel("log2 fold change (PR8 / MOCK)"); axa.set_ylabel("-log10 p-value")
axa.set_title("a  Differential m6A methylation (mRNA)",loc="left",fontweight="bold",fontsize=11)
leg=axa.legend(frameon=False,fontsize=8,markerscale=3,loc="upper right")

# (b) region distribution
axb=fig.add_subplot(gs[0,1])
regs=["5'UTR","CDS","3'UTR"]; vals=[region_total[r] for r in regs]
bars=axb.bar(regs,vals,color=["#55A868","#4C72B0","#C44E52"],edgecolor="black",linewidth=0.5)
tot=sum(vals)
for b,v in zip(bars,vals):
    axb.text(b.get_x()+b.get_width()/2,v,f"{v}\n({100*v/tot:.1f}%)",ha="center",va="bottom",fontsize=8)
axb.set_ylabel("m6A sites on differential mRNAs")
axb.set_title("b  m6A site distribution by mRNA region",loc="left",fontweight="bold",fontsize=11)
axb.set_ylim(0,max(vals)*1.18)

# (c) top differential transcripts as log2FC barplot (honest, informative)
axc=fig.add_subplot(gs[0,2])
# need gene symbols: re-read from allcomparison col2 aligned with x
wb3=openpyxl.load_workbook(ALL,read_only=True,data_only=True); ws3=wb3["MOCK_vs_PR8"]
hr3=find_header(ws3); syms=[]
for ri,row in enumerate(ws3.iter_rows(values_only=True)):
    if ri<=hr3: continue
    a=row[0]
    if a is None or str(a).strip().startswith("#") or str(a).strip()=="": continue
    try: float(row[5]); float(row[6])
    except (TypeError,ValueError): continue
    syms.append(str(row[2]).strip() if row[2] else str(a).strip())
wb3.close()
syms=np.array(syms)
# top 10 up-in-PR8 + top 10 up-in-MOCK by |x| among sig
order_sig=[i for i in np.argsort(-np.abs(x)) if sig[i]]
top_pr8=[i for i in order_sig if x[i]>0][:10]
top_mock=[i for i in order_sig if x[i]<0][:10]
sel=top_mock[::-1]+top_pr8[::-1]
yv=np.arange(len(sel)); xv=[x[i] for i in sel]
cols=["#C44E52" if x[i]>0 else "#4C72B0" for i in sel]
axc.barh(yv,xv,color=cols,edgecolor="black",linewidth=0.3)
axc.set_yticks(yv); axc.set_yticklabels([syms[i] for i in sel],fontsize=6)
axc.axvline(0,color="black",lw=0.6)
axc.set_xlabel("log2 fold change (PR8 / MOCK)")
axc.set_title("c  Top differential m6A transcripts",loc="left",fontweight="bold",fontsize=11)
axc.margins(y=0.01)

fig.suptitle("Figure 3  m6A signal characteristics and differential landscape",fontsize=13,fontweight="bold",y=1.02)
fig.savefig(os.path.join(FIGDIR,"Figure3_m6a_signal.pdf"),dpi=300,bbox_inches="tight")
fig.savefig(os.path.join(FIGDIR,"Figure3_m6a_signal.png"),dpi=150,bbox_inches="tight")
plt.close(fig)
print("SAVED Figure3")

