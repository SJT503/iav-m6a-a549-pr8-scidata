# -*- coding: utf-8 -*-
"""13_supp_figs.py — 补充图（重修美化版，对齐 Figure 2 视觉风格）
SuppFig1: (a) m6A 调控酶 log2FC 热图(Writer/Eraser/Reader 分组) + (b) GO(BP) top10 富集
SuppFig2: KEGG pathway top10 富集 (hyper m6A)
真实数据：m6A酶表达情况xlsx.xlsx + GO/Pathway result.txt
数据真实性铁律：全部来自实读，禁止编造。read_only=True 只读原始报告。
"""
import sys; sys.path.insert(0, r"D:\IAV_m6A\SciData_descriptor\scripts")
import _figstyle as fs
import os, csv
import numpy as np, openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

fs.apply_style()

FIGDIR = r"D:\IAV_m6A\SciData_descriptor\results\figures"
os.makedirs(FIGDIR, exist_ok=True)
ENZ  = r"D:\IAV_m6A\Report\Report\A549-PR8的m6A数据的再分析 - 2022-4-13\m6A相关的酶\m6A酶表达情况xlsx.xlsx"
GO   = r"D:\IAV_m6A\Report\Report\Report\PaperData\GO Analysis Report\hyper_MOCK_vs_PR8(mRNA_Quantity)\BP_result.txt"
KEGG = r"D:\IAV_m6A\Report\Report\Report\PaperData\Pathway Analysis Report\hyper_MOCK_vs_PR8(mRNA_Quantity)\hsa_pathwayResult.txt"

# group accent colors (Dark2-style, distinct from RdBu diverging heatmap)
GRP_COLOR = {"Writer": "#7B3294", "Eraser": "#D95F02", "Reader": "#1B9E77"}

# ============================================================
# 1) Load m6A enzymes — aggregate per gene (mean of probes), keep class
# ============================================================
wb = openpyxl.load_workbook(ENZ, read_only=True, data_only=True)
classes = {"WRITER": "Writer", "READERS": "Reader", "ERASERS": "Eraser"}
rows = []
for sh, cls in classes.items():
    ws = wb[sh]
    agg = {}
    for row in ws.iter_rows(values_only=True):
        g = row[4] if len(row) > 4 else None
        if not g or str(g).strip() in ("GeneSymbol", "", "None"):
            continue
        try:
            fge = float(row[1]); fm = float(row[2])
        except (TypeError, ValueError):
            continue
        agg.setdefault(str(g).strip(), []).append((fge, fm))
    for g, v in agg.items():
        ge = np.mean([a for a, b in v]); m6 = np.mean([b for a, b in v])
        rows.append((cls, g, ge, m6))
wb.close()

order = {"Writer": 0, "Eraser": 1, "Reader": 2}
rows.sort(key=lambda r: (order[r[0]], r[1]))
genes = [r[1] for r in rows]
mat = np.array([[np.log2(r[2]), np.log2(r[3])] for r in rows])
print("enzymes:", len(rows))

# ============================================================
# 2) Enrichment table loader
# ============================================================
def _trunc(s, n):
    s = str(s)
    # strip trailing " - Homo sapiens (human)" for KEGG readability
    s = s.replace(" - Homo sapiens (human)", "")
    return s if len(s) <= n else s[:n - 1] + "…"

def read_tab(path, n=10, term_i=1, p_i=8):
    """Read tab file, keep top-n by file order; return [(term, pvalue), ...]."""
    out = []
    with open(path, encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f, delimiter="\t")
        next(rd)  # header
        for r in rd:
            if len(r) <= max(term_i, p_i):
                continue
            try:
                p = float(r[p_i])
            except (ValueError, IndexError):
                continue
            out.append((r[term_i], p))
            if len(out) >= n:
                break
    return out

go = read_tab(GO, 10, term_i=1, p_i=8)            # GO BP: Term col1, Pvalue col8
kegg = read_tab(KEGG, 10, term_i=1, p_i=3)        # KEGG: Definition col1, Fisher-Pvalue col3
print("GO top:", len(go), "| KEGG top:", len(kegg))

# ============================================================
# 3) SuppFigure 1 — (a) enzyme heatmap  (b) GO BP barh
# ============================================================
fig = plt.figure(figsize=(11, 6.0))
gs = GridSpec(1, 2, figure=fig, wspace=0.90, width_ratios=[1.0, 1.45])

# ---- (a) enzyme heatmap ----
ax1 = fig.add_subplot(gs[0, 0])
VLIM = 1.2
im = ax1.imshow(mat, aspect="auto", cmap=fs.DIV_CMAP, vmin=-VLIM, vmax=VLIM)
ax1.set_xticks([0, 1])
ax1.set_xticklabels(["Expression\nlog2FC", "m6A\nlog2FC"], fontsize=8.5)
ax1.set_yticks(range(len(genes)))
ax1.set_yticklabels(genes, fontsize=8)
ax1.tick_params(length=0)
# overlay numeric values
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        val = mat[i, j]
        ax1.text(j, i, f"{val:+.2f}", ha="center", va="center", fontsize=6.2,
                 color="white" if abs(val) > 0.7 else "#222222")
# class separators + RIGHT-side colored group sidebar.
# Colorbar is moved to the BOTTOM (horizontal) so the entire right side is free for
# the group sidebar -> no collision between group labels and the colorbar.
cum = 0
for cls in ["Writer", "Eraser", "Reader"]:
    n = sum(1 for r in rows if r[0] == cls)
    # colored band just right of the heatmap cells (cells end at x=1.5)
    ax1.add_patch(Rectangle((1.60, cum - 0.5), 0.18, n, facecolor=GRP_COLOR[cls],
                            edgecolor="none", clip_on=False))
    # group name to the right of the band, rotated, color-matched
    ax1.text(1.96, cum + n / 2 - 0.5, cls, rotation=90, va="center", ha="center",
             fontsize=9, fontweight="bold", color=GRP_COLOR[cls], clip_on=False)
    if cum > 0:
        ax1.axhline(cum - 0.5, color="black", lw=1.4)
    cum += n
ax1.set_xlim(-0.5, 2.35)
fs.panel_title(ax1, "a", "m6A regulatory enzymes")
# horizontal colorbar below the heatmap
cb = fig.colorbar(im, ax=ax1, orientation="horizontal", fraction=0.05, pad=0.10, aspect=30)
cb.set_label("log2 fold change (PR8/MOCK)", fontsize=8.5)
cb.ax.tick_params(labelsize=7.5)

# ---- (b) GO BP top10 horizontal bars ----
ax2 = fig.add_subplot(gs[0, 1])
terms = [_trunc(t, 34) for t, p in go][::-1]
scores = [-np.log10(p) for t, p in go][::-1]
y = range(len(terms))
ax2.barh(y, scores, color=fs.PALETTE["green"], edgecolor="black", linewidth=0.5, height=0.72)
ax2.set_yticks(y)
ax2.set_yticklabels(terms, fontsize=7.5)
ax2.set_xlabel("$-\\log_{10}$ p-value")
ax2.set_ylim(-0.6, len(terms) - 0.4)
for spine in ("top", "right"):
    ax2.spines[spine].set_visible(False)
# annotate values at bar ends
xmax = max(scores)
for yi, s in zip(y, scores):
    ax2.text(s + xmax * 0.012, yi, f"{s:.1f}", va="center", ha="left", fontsize=6.5, color="#333333")
ax2.set_xlim(0, xmax * 1.12)
fs.panel_title(ax2, "b", "GO biological process (hyper-m6A, top 10)")

fig.suptitle("Supplementary Figure 1  m6A regulatory enzymes and functional enrichment",
             **fs.SUPTITLE_KW)
p1 = os.path.join(FIGDIR, "SuppFigure1_enzymes_GO.pdf")
p1png = os.path.join(FIGDIR, "SuppFigure1_enzymes_GO.png")
fig.savefig(p1, dpi=300, bbox_inches="tight")
fig.savefig(p1png, dpi=150, bbox_inches="tight")
plt.close(fig)
print("SAVED SuppFig1:", p1)

# ============================================================
# 4) SuppFigure 2 — KEGG top10 horizontal bars (blue)
# ============================================================
fig2, ax = plt.subplots(figsize=(8.0, 4.6))
terms = [_trunc(t, 46) for t, p in kegg][::-1]
scores = [-np.log10(p) for t, p in kegg][::-1]
y = range(len(terms))
ax.barh(y, scores, color=fs.PALETTE["MOCK"], edgecolor="black", linewidth=0.5, height=0.72)
ax.set_yticks(y)
ax.set_yticklabels(terms, fontsize=8)
ax.set_xlabel("$-\\log_{10}$ Fisher p-value")
ax.set_ylim(-0.6, len(terms) - 0.4)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
xmax = max(scores)
for yi, s in zip(y, scores):
    ax.text(s + xmax * 0.012, yi, f"{s:.1f}", va="center", ha="left", fontsize=6.8, color="#333333")
ax.set_xlim(0, xmax * 1.12)
fs.panel_title(ax, "", "Supplementary Figure 2  KEGG pathway enrichment (hyper-m6A, top 10)", fontsize=12)

p2 = os.path.join(FIGDIR, "SuppFigure2_KEGG.pdf")
p2png = os.path.join(FIGDIR, "SuppFigure2_KEGG.png")
fig2.savefig(p2, dpi=300, bbox_inches="tight")
fig2.savefig(p2png, dpi=150, bbox_inches="tight")
plt.close(fig2)
print("SAVED SuppFig2:", p2)

# ============================================================
# 5) Programmatic self-check (PIL+numpy; Read 渲染不出图像)
# ============================================================
from PIL import Image
def qc(png, panels):
    im = np.asarray(Image.open(png).convert("RGB"), dtype=float)
    h, w, _ = im.shape
    nonwhite = (im.min(axis=2) < 245)
    overall = nonwhite.mean()
    # edge margin check: outer 1% border should be mostly white (no clipping)
    b = max(2, int(min(h, w) * 0.01))
    edge = np.concatenate([
        nonwhite[:b, :].ravel(), nonwhite[-b:, :].ravel(),
        nonwhite[:, :b].ravel(), nonwhite[:, -b:].ravel()])
    edge_ink = edge.mean()
    print(f"\nQC {os.path.basename(png)}: size={w}x{h}, overall non-white={overall:.3f}, edge ink={edge_ink:.3f}")
    for name, (r0, r1, c0, c1) in panels.items():
        sub = nonwhite[int(r0*h):int(r1*h), int(c0*w):int(c1*w)]
        frac = sub.mean()
        print(f"   panel {name}: non-white={frac:.3f}  {'OK' if frac > 0.01 else 'BLANK!!'}")
    return overall

o1 = qc(p1png, {"a(left)": (0.1, 0.95, 0.0, 0.45), "b(right)": (0.1, 0.95, 0.45, 1.0)})
o2 = qc(p2png, {"full": (0.1, 0.95, 0.0, 1.0)})
print(f"\nSuppFig2 overall non-white={o2:.3f}  {'OK (>0.1)' if o2 > 0.1 else 'FAIL <0.1 (blank bug!)'}")
print("\nALL DONE.")
