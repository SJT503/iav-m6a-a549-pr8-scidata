# -*- coding: utf-8 -*-
"""12_fig1_design.py — Fig 1 研究设计与数据流程图（示意图，无统计数据点）
重写美化版：分阶段蛇形流程，对齐 Fig2 视觉风格（共享 _figstyle 调色板/字体/suptitle）。
真实参数全部来自 anchored_facts.md，禁止改数字：
  A549(ATCC) + IAV PR8 MOI 0.5 24h, MOCK 对照, n=3/组(6 样本);
  m6A MeRIP(anti-m6A IP) -> IP=Cy5/Sup=Cy3 双色标记 + spike-in;
  Arraystar Human mRNA&lncRNA Epitranscriptomic Microarray (8x60K);
  杂交 65C 17h -> Agilent G2505C -> Feature Extraction v11.0.1.1;
  spike-in 归一化 -> 40,956 探针 x 6 -> m6A quantity(IP)+expression(Sup);
  差异分析 |FC|>=2, p<0.05, mRNA/lncRNA/ncRNA.
"""
import os, sys
sys.path.insert(0, r"D:\IAV_m6A\SciData_descriptor\scripts")
import _figstyle as fs

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import to_rgb

fs.apply_style()

FIGDIR = r"D:\IAV_m6A\SciData_descriptor\results\figures"
os.makedirs(FIGDIR, exist_ok=True)

P = fs.PALETTE
BLUE, RED, GREEN, NEUT = P["MOCK"], P["PR8"], P["green"], P["neutral"]


def tint(color, white=0.82):
    """Blend a hex color toward white (white=fraction of white)."""
    r, g, b = to_rgb(color)
    return (r + (1 - r) * white, g + (1 - g) * white, b + (1 - b) * white)


# ------------------------------------------------------------------ figure
fig, ax = plt.subplots(figsize=(10.2, 8.0))
ax.set_xlim(2.5, 92.5)
ax.set_ylim(7.5, 100.5)
ax.axis("off")
fig.subplots_adjust(left=0.015, right=0.985, top=0.915, bottom=0.02)

EDGE_LW = 0.8
ARR_COL = "#555555"


def box(x, y, w, h, text, family, emphasize=False):
    """Rounded box; family in {'blue','red','green','grey'} drives fill+edge."""
    fam = {"blue": BLUE, "red": RED, "green": GREEN, "grey": NEUT}[family]
    fc = tint(fam, 0.80 if emphasize else 0.84)
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.0,rounding_size=2.2",
        fc=fc, ec=fam, lw=EDGE_LW, zorder=3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=8.2, color="#222222", zorder=4, linespacing=1.32)


def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=13,
        lw=1.2, color=ARR_COL, shrinkA=0, shrinkB=0, zorder=2))


def stage_panel(y0, y1, label, lblcolor):
    ax.add_patch(FancyBboxPatch(
        (9, y0), 82, y1 - y0,
        boxstyle="round,pad=0.0,rounding_size=2.2",
        fc="#F5F5F5", ec="none", zorder=0))
    ax.text(5.0, (y0 + y1) / 2, label, ha="center", va="center",
            rotation=90, fontsize=9.2, fontweight="bold",
            color=lblcolor, zorder=1)


# ------------------------------------------------------------- stage panels
stage_panel(77.5, 99.0, "Biological model", "#3A3A3A")
stage_panel(55.5, 73.5, "m6A enrichment", "#3A3A3A")
stage_panel(33.5, 51.5, "Microarray & scanning", "#3A3A3A")
stage_panel(9.5, 29.5, "Data processing & records", "#3A3A3A")

# ----------------------------------------------------- Stage 1 (left->right)
box(11, 82.0, 19, 11, "A549 cells\n(ATCC, human\nalveolar epithelial)", "blue")
box(37, 88.3, 21, 8.2, "Mock control\n(n = 3)", "blue")
box(37, 79.0, 21, 8.2, "IAV PR8 infection\nMOI 0.5, 24 h  (n = 3)", "red")
box(66, 81.0, 24, 12.0,
    "Total RNA\n6 samples (MOCK 1–3,\nPR8 1–3)\nNanoDrop / Bioanalyzer QC", "grey")

arrow(30, 88.5, 37, 92.4)            # A549 -> Mock
arrow(30, 86.0, 37, 83.1)            # A549 -> PR8
arrow(58, 92.4, 66, 89.0)            # Mock -> Total RNA
arrow(58, 83.1, 66, 86.0)            # PR8  -> Total RNA

# ----------------------------------------------------- Stage 2 (right->left)
box(66, 58.0, 24, 12.0,
    "m6A MeRIP\nanti-m6A antibody IP\n+ magnetic-bead pulldown", "green")
box(37, 58.0, 24, 12.0,
    "Two-colour labelling\nIP = Cy5 / Sup = Cy3\n+ m6A spike-in controls", "green")

arrow(78, 81.0, 78, 70.0)            # Total RNA -> MeRIP
arrow(66, 64.0, 61, 64.0)            # MeRIP -> labelling

# ----------------------------------------------------- Stage 3 (left->right)
box(37, 37.0, 24, 12.0,
    "Arraystar Human mRNA &\nlncRNA Epitranscriptomic\nMicroarray (8 × 60K)", "grey")
box(66, 37.0, 24, 12.0,
    "Hybridize 65 °C, 17 h\nAgilent G2505C scan\nFeature Extraction v11.0.1.1", "grey")

arrow(49, 58.0, 49, 49.0)            # labelling -> microarray
arrow(61, 43.0, 66, 43.0)            # microarray -> scan

# ----------------------------------------------------- Stage 4 (right->left)
box(64, 13.0, 26, 13.0,
    "Spike-in normalization\n40,956 probes × 6 samples\nm6A quantity (IP) +\nexpression (Sup)", "grey")
box(33, 13.0, 28, 13.0,
    "Differential m6A & expression\nmRNA / lncRNA / ncRNA\n|FC| ≥ 2,  p < 0.05", "green", emphasize=True)

arrow(78, 37.0, 77, 26.0)            # scan -> normalization
arrow(64, 19.5, 61, 19.5)            # normalization -> differential

# ------------------------------------------------------------------ title
fig.suptitle("Figure 1  Study design and data workflow", **fs.SUPTITLE_KW)

out_pdf = os.path.join(FIGDIR, "Figure1_design.pdf")
out_png = os.path.join(FIGDIR, "Figure1_design.png")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
fig.savefig(out_png, dpi=150, bbox_inches="tight")
plt.close(fig)
print("SAVED:", out_pdf)
print("SAVED:", out_png)

# ============================================================ self-check
import numpy as np
from PIL import Image

img = Image.open(out_png).convert("RGB")
A = np.asarray(img)
H, W, _ = A.shape
nonwhite = np.any(A < 245, axis=2)          # True where pixel is not near-white
ratio = nonwhite.mean()

# four quadrants content check
hh, hw = H // 2, W // 2
quads = {
    "TL": nonwhite[:hh, :hw].mean(),
    "TR": nonwhite[:hh, hw:].mean(),
    "BL": nonwhite[hh:, :hw].mean(),
    "BR": nonwhite[hh:, hw:].mean(),
}

# min margin (fully-white border rows/cols from each edge)
def first_content(line_has):
    for i, v in enumerate(line_has):
        if v:
            return i
    return len(line_has)

rows_has = nonwhite.any(axis=1)
cols_has = nonwhite.any(axis=0)
top_m = first_content(rows_has)
bot_m = first_content(rows_has[::-1])
left_m = first_content(cols_has)
right_m = first_content(cols_has[::-1])

print("\n===== SELF-CHECK (Figure1_design.png) =====")
print(f"size              : {W} x {H} px")
print(f"non-white ratio   : {ratio*100:.2f}%")
print(f"quadrant content  : " +
      ", ".join(f"{k}={v*100:.2f}%" for k, v in quads.items()))
print(f"min edge margins  : top={top_m}px bottom={bot_m}px "
      f"left={left_m}px right={right_m}px")
all_quads_ok = all(v > 0.005 for v in quads.values())
no_clip = min(top_m, bot_m, left_m, right_m) >= 2
print(f"all 4 quadrants non-empty (>0.5%) : {all_quads_ok}")
print(f"no edge clipping (margin >= 2px)  : {no_clip}")
print(f"VERDICT: {'PASS' if (all_quads_ok and no_clip and 0.05 < ratio < 0.6) else 'REVIEW'}")
