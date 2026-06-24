# -*- coding: utf-8 -*-
"""_qc_overlap.py — 精确检测 SuppFig1 分组标签/侧边条与热图单元格的几何碰撞
用 matplotlib renderer 拿到每个文本/patch 的真实显示坐标(像素bbox)，
再算热图 imshow 的数据区像素范围，判断是否重叠。不靠肉眼。
"""
import sys; sys.path.insert(0, r"D:\IAV_m6A\SciData_descriptor\scripts")
import _figstyle as fs
import numpy as np, openpyxl
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
fs.apply_style()

ENZ = r"D:\IAV_m6A\Report\Report\A549-PR8的m6A数据的再分析 - 2022-4-13\m6A相关的酶\m6A酶表达情况xlsx.xlsx"
GRP_COLOR = {"Writer": "#7B3294", "Eraser": "#D95F02", "Reader": "#1B9E77"}
wb = openpyxl.load_workbook(ENZ, read_only=True, data_only=True)
classes = {"WRITER": "Writer", "READERS": "Reader", "ERASERS": "Eraser"}
rows = []
for sh, cls in classes.items():
    ws = wb[sh]; agg = {}
    for row in ws.iter_rows(values_only=True):
        g = row[4] if len(row) > 4 else None
        if not g or str(g).strip() in ("GeneSymbol","","None"): continue
        try: fge=float(row[1]); fm=float(row[2])
        except (TypeError,ValueError): continue
        agg.setdefault(str(g).strip(),[]).append((fge,fm))
    for g,v in agg.items():
        rows.append((cls,g,np.mean([a for a,b in v]),np.mean([b for a,b in v])))
wb.close()
order={"Writer":0,"Eraser":1,"Reader":2}
rows.sort(key=lambda r:(order[r[0]],r[1]))
genes=[r[1] for r in rows]
mat=np.array([[np.log2(r[2]),np.log2(r[3])] for r in rows])

fig=plt.figure(figsize=(11,6.0))
gs=GridSpec(1,2,figure=fig,wspace=0.55,width_ratios=[1.0,1.45])
ax1=fig.add_subplot(gs[0,0])
im=ax1.imshow(mat,aspect="auto",cmap=fs.DIV_CMAP,vmin=-1.2,vmax=1.2)
ax1.set_yticks(range(len(genes))); ax1.set_yticklabels(genes,fontsize=8)
# group labels (current geometry)
cum=0; label_artists=[]
for cls in ["Writer","Eraser","Reader"]:
    n=sum(1 for r in rows if r[0]==cls)
    ax1.add_patch(Rectangle((1.60,cum-0.5),0.18,n,facecolor=GRP_COLOR[cls],edgecolor="none",clip_on=False))
    t=ax1.text(1.96,cum+n/2-0.5,cls,rotation=90,va="center",ha="center",fontsize=9,fontweight="bold",clip_on=False)
    label_artists.append((cls,t))
    cum+=n
ax1.set_xlim(-0.5,2.35)

fig.canvas.draw()
r=fig.canvas.get_renderer()
# heatmap cell region in display pixels: data x from -0.5 to 1.5 (2 columns)
inv=ax1.transData
def datapx(x,y): return ax1.transData.transform((x,y))
# heatmap right edge (cells end at x=1.5) in px
cell_right_px = datapx(1.5,0)[0]
band_left_px  = datapx(1.60,0)[0]
print(f"heatmap cells right edge x-px = {cell_right_px:.1f}")
print(f"group color band left x-px    = {band_left_px:.1f}  (gap={band_left_px-cell_right_px:.1f}px, should be >0)")
for cls,t in label_artists:
    bb=t.get_window_extent(renderer=r)
    overlap = bb.x0 < cell_right_px
    print(f"  label '{cls}': text bbox x0={bb.x0:.1f} x1={bb.x1:.1f}  -> {'OVERLAP heatmap!' if overlap else 'clear of heatmap'}")
plt.close(fig)
