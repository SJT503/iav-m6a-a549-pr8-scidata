# -*- coding: utf-8 -*-
"""10_fig2_quality.py — Fig 2 数据质量与样本关系（Technical Validation 核心）
真实数据：Matrix.xlsx(40,956探针×12) + pca.txt(mRNA_Quantity)
面板：(a) PCA  (b) 样本相关热图  (c) 重复 Pearson R 散点  (d) 归一化强度箱线
数据真实性铁律：全部来自实读，禁止编造。
"""
import os, numpy as np, openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

RESDIR = r"D:\IAV_m6A\SciData_descriptor\results"
FIGDIR = r"D:\IAV_m6A\SciData_descriptor\results\figures"
os.makedirs(FIGDIR, exist_ok=True)
MATRIX = r"D:\IAV_m6A\Report\Report\Report\PaperData\Raw Data Files\Matrix.xlsx"
PCA = r"D:\IAV_m6A\Report\Report\盛江涛_H2108172_Human_m6A-mRNA&lncRNA Epi_补充分析_20220913_mk_zrb\pca\mRNA_Quantity_PR8_vs_MOCK\pca.txt"

SAMPLES = ["MOCK_1","MOCK_2","MOCK_3","PR8_1","PR8_2","PR8_3"]
COLORS = {"MOCK":"#4C72B0","PR8":"#C44E52"}

# ---- load IP normalized matrix (6 IP columns) ----
wb = openpyxl.load_workbook(MATRIX, read_only=True, data_only=True)
ws = wb["Matrix"]
rows = ws.iter_rows(values_only=True)
hdr = next(rows)
# header: probeID, MOCK_1(IP,Normalized)...PR8_3(IP,Normalized), then Sup...
ip_idx = [i for i,h in enumerate(hdr) if h and "(IP,Normalized)" in str(h)]
ip_names = [str(hdr[i]).split("(")[0] for i in ip_idx]
data = []
for r in rows:
    vals = [r[i] for i in ip_idx]
    if all(isinstance(v,(int,float)) for v in vals):
        data.append(vals)
wb.close()
M = np.array(data, dtype=float)  # probes x 6
print("IP matrix:", M.shape, "samples:", ip_names)
np.save(os.path.join(RESDIR,"_ip_matrix.npy"), M)
with open(os.path.join(RESDIR,"_ip_meta.txt"),"w") as f:
    f.write(f"shape={M.shape}\nsamples={ip_names}\n")

grp = ["MOCK" if n.startswith("MOCK") else "PR8" for n in ip_names]

# ---- load PCA (real, from provider supplementary) ----
pca = {}
with open(PCA, encoding="utf-8") as f:
    next(f)  # header
    for line in f:
        p = line.strip().split("\t")
        if len(p) >= 3 and p[0] in SAMPLES:
            pca[p[0]] = (float(p[1]), float(p[2]))
# provider reported variance proportions for mRNA_Quantity: PC1 50.3%, PC2 48.0% (pca.txt SS/ProportionVar)
PC1V, PC2V = 50.3, 48.0

# ---- correlation matrix (Pearson) ----
R = np.corrcoef(M.T)  # 6x6

# ================= FIGURE =================
fig = plt.figure(figsize=(11, 9))
gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.30)

# (a) PCA
axa = fig.add_subplot(gs[0,0])
for n in ip_names:
    x,y = pca[n]
    axa.scatter(x, y, s=120, color=COLORS[("MOCK" if n.startswith("MOCK") else "PR8")],
                edgecolor="black", linewidth=0.6, zorder=3)
    axa.annotate(n, (x,y), fontsize=7, xytext=(4,4), textcoords="offset points")
axa.set_xlabel(f"PC1 ({PC1V:.1f}%)"); axa.set_ylabel(f"PC2 ({PC2V:.1f}%)")
axa.set_title("a  PCA of m6A quantity profiles", loc="left", fontweight="bold", fontsize=11)
from matplotlib.lines import Line2D
axa.legend(handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor=COLORS[g],
           markeredgecolor='k',markersize=9,label=g) for g in ["MOCK","PR8"]],
           frameon=False, fontsize=9, loc="best")

# (b) correlation heatmap
axb = fig.add_subplot(gs[0,1])
im = axb.imshow(R, cmap="viridis", vmin=0.85, vmax=1.0)
axb.set_xticks(range(6)); axb.set_yticks(range(6))
axb.set_xticklabels(ip_names, rotation=45, ha="right", fontsize=7)
axb.set_yticklabels(ip_names, fontsize=7)
for i in range(6):
    for j in range(6):
        # viridis: dark (low) -> light/yellow (high). White text only on dark cells.
        axb.text(j,i,f"{R[i,j]:.2f}",ha="center",va="center",
                 color="white" if R[i,j]<0.93 else "black",fontsize=7)
axb.set_title("b  Sample-sample correlation (Pearson)", loc="left", fontweight="bold", fontsize=11)
fig.colorbar(im, ax=axb, fraction=0.046, pad=0.04, label="Pearson r")

# (c) within-group replicate scatter (representative: MOCK_1 vs MOCK_2, PR8_1 vs PR8_2)
axc = fig.add_subplot(gs[1,0])
pairs = [(0,1,"MOCK","MOCK_1 vs MOCK_2"),(3,4,"PR8","PR8_1 vs PR8_2")]
for xi,yi,g,lab in pairs:
    r = R[xi,yi]
    axc.scatter(M[:,xi], M[:,yi], s=2, alpha=0.12, color=COLORS[g],
                label=f"{lab} (r={r:.3f})", rasterized=True)
axc.set_xlabel("Normalized m6A signal (replicate 1)")
axc.set_ylabel("Normalized m6A signal (replicate 2)")
axc.set_title("c  Within-group reproducibility", loc="left", fontweight="bold", fontsize=11)
leg = axc.legend(frameon=False, fontsize=8, markerscale=4)
for lh in leg.legend_handles:
    lh.set_alpha(1)

# (d) per-sample normalized intensity boxplot
axd = fig.add_subplot(gs[1,1])
bp = axd.boxplot([M[:,i] for i in range(6)], positions=range(6), widths=0.6,
                 patch_artist=True, showfliers=False)
for i,box in enumerate(bp['boxes']):
    box.set(facecolor=COLORS[grp[i]], alpha=0.7)
for med in bp['medians']: med.set(color="black")
axd.set_xticks(range(6)); axd.set_xticklabels(ip_names, rotation=45, ha="right", fontsize=7)
axd.set_ylabel("Normalized m6A signal (log2)")
axd.set_title("d  Per-sample signal distribution", loc="left", fontweight="bold", fontsize=11)

fig.suptitle("Figure 2  Data quality and sample relationships", fontsize=13, fontweight="bold", y=0.98)
out_pdf = os.path.join(FIGDIR,"Figure2_quality.pdf")
out_png = os.path.join(FIGDIR,"Figure2_quality.png")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
fig.savefig(out_png, dpi=150, bbox_inches="tight")
plt.close(fig)
print("SAVED:", out_pdf)
# dump real r values for legend/text use
with open(os.path.join(RESDIR,"_fig2_stats.txt"),"w") as f:
    f.write("Pearson r matrix (IP normalized, 6 samples):\n")
    f.write("samples: "+",".join(ip_names)+"\n")
    for i in range(6):
        f.write(" ".join(f"{R[i,j]:.4f}" for j in range(6))+"\n")
    f.write(f"\nWithin-MOCK r: {R[0,1]:.4f},{R[0,2]:.4f},{R[1,2]:.4f}\n")
    f.write(f"Within-PR8 r: {R[3,4]:.4f},{R[3,5]:.4f},{R[4,5]:.4f}\n")

