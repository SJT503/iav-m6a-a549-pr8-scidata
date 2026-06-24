# -*- coding: utf-8 -*-
"""_qc_cross_panel.py — 跨面板碰撞检测：panel a 右伸标签 vs panel b 左侧 GO 文字
扫 wspace,找到两面板文字不重叠的最小间距。
"""
import sys; sys.path.insert(0, r"D:\IAV_m6A\SciData_descriptor\scripts")
import _figstyle as fs
import numpy as np, openpyxl, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
fs.apply_style()

ENZ = r"D:\IAV_m6A\Report\Report\A549-PR8的m6A数据的再分析 - 2022-4-13\m6A相关的酶\m6A酶表达情况xlsx.xlsx"
GO  = r"D:\IAV_m6A\Report\Report\Report\PaperData\GO Analysis Report\hyper_MOCK_vs_PR8(mRNA_Quantity)\BP_result.txt"
GRP_COLOR={"Writer":"#7B3294","Eraser":"#D95F02","Reader":"#1B9E77"}
wb=openpyxl.load_workbook(ENZ,read_only=True,data_only=True)
classes={"WRITER":"Writer","READERS":"Reader","ERASERS":"Eraser"}; rows=[]
for sh,cls in classes.items():
    ws=wb[sh]; agg={}
    for row in ws.iter_rows(values_only=True):
        g=row[4] if len(row)>4 else None
        if not g or str(g).strip() in ("GeneSymbol","","None"): continue
        try: fge=float(row[1]); fm=float(row[2])
        except (TypeError,ValueError): continue
        agg.setdefault(str(g).strip(),[]).append((fge,fm))
    for g,v in agg.items(): rows.append((cls,g,np.mean([a for a,b in v]),np.mean([b for a,b in v])))
wb.close()
order={"Writer":0,"Eraser":1,"Reader":2}; rows.sort(key=lambda r:(order[r[0]],r[1]))
genes=[r[1] for r in rows]; mat=np.array([[np.log2(r[2]),np.log2(r[3])] for r in rows])

def _trunc(s,n):
    s=str(s).replace(" - Homo sapiens (human)","")
    return s if len(s)<=n else s[:n-1]+"…"
go=[]
with open(GO,encoding="utf-8",errors="replace") as f:
    rd=csv.reader(f,delimiter="\t"); next(rd)
    for r in rd:
        if len(r)<=8: continue
        try: p=float(r[8])
        except (ValueError,IndexError): continue
        go.append((r[1],p))
        if len(go)>=10: break

def test(wspace, wr, label_x, go_trunc):
    fig=plt.figure(figsize=(11,6.0))
    gs=GridSpec(1,2,figure=fig,wspace=wspace,width_ratios=wr)
    ax1=fig.add_subplot(gs[0,0])
    ax1.imshow(mat,aspect="auto",cmap=fs.DIV_CMAP,vmin=-1.2,vmax=1.2)
    ax1.set_yticks(range(len(genes))); ax1.set_yticklabels(genes,fontsize=8)
    cum=0; alab=[]
    for cls in ["Writer","Eraser","Reader"]:
        n=sum(1 for r in rows if r[0]==cls)
        ax1.add_patch(Rectangle((1.60,cum-0.5),0.18,n,facecolor=GRP_COLOR[cls],edgecolor="none",clip_on=False))
        t=ax1.text(label_x,cum+n/2-0.5,cls,rotation=90,va="center",ha="center",fontsize=9,fontweight="bold",clip_on=False)
        alab.append(t); cum+=n
    ax1.set_xlim(-0.5,2.35)
    ax2=fig.add_subplot(gs[0,1])
    terms=[_trunc(t,go_trunc) for t,p in go][::-1]; sc=[-np.log10(p) for t,p in go][::-1]
    ax2.barh(range(len(terms)),sc,height=0.72); ax2.set_yticks(range(len(terms)))
    blab=ax2.set_yticklabels(terms,fontsize=7.5)
    fig.canvas.draw(); r=fig.canvas.get_renderer()
    a_right=max(t.get_window_extent(renderer=r).x1 for t in alab)
    b_left=min(t.get_window_extent(renderer=r).x0 for t in blab)
    plt.close(fig)
    return a_right,b_left,b_left-a_right

print("=== actual script uses trunc=42; scan wspace + trunc ===")
print("wspace  trunc  a_label_right  b_text_left   gap(px)  verdict")
for tr in [42, 34, 28]:
    for ws in [0.55, 0.9, 1.2]:
        ar,bl,gap=test(ws,[1.0,1.45],1.96,tr)
        print(f"{ws:.2f}    {tr}     {ar:8.1f}    {bl:8.1f}    {gap:7.1f}  {'OK' if gap>10 else 'OVERLAP'}")
