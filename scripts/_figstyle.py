# -*- coding: utf-8 -*-
"""_figstyle.py — 共享出图样式基准（以 Figure 2 为准，所有图对齐）
用法：import 后 apply_style() 设全局 rcParams；用 PALETTE / panel_title() / SUPTITLE_KW 等辅助。
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- 调色板（焊死，与 Fig2 一致）----
PALETTE = {
    "MOCK": "#4C72B0",   # 蓝
    "PR8":  "#C44E52",   # 红
    "green":"#55A868",
    "neutral":"#8C8C8C",
    "grey_pt":"#BBBBBB",
}
SEQ_CMAP = "viridis"     # 顺序型（相关热图等）
DIV_CMAP = "RdBu_r"      # 发散型（z-score/log2FC 热图）

# ---- 全局 rcParams ----
def apply_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9.5,
        "axes.linewidth": 0.8,
        "axes.edgecolor": "#333333",
        "axes.grid": False,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,   # embed TrueType (editable text in Illustrator)
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

# ---- panel 小标题 (a/b/c, 左对齐粗体 11pt) ----
def panel_title(ax, letter, text, fontsize=11):
    ax.set_title(f"{letter}  {text}", loc="left", fontweight="bold", fontsize=fontsize, pad=6)

# ---- suptitle 统一参数 ----
SUPTITLE_KW = dict(fontsize=13, fontweight="bold", y=1.00)

# ---- colorbar 统一参数 ----
CBAR_KW = dict(fraction=0.046, pad=0.04)

# ---- 分组着色辅助 ----
def group_color(sample_name):
    return PALETTE["MOCK"] if str(sample_name).upper().startswith("MOCK") else PALETTE["PR8"]
