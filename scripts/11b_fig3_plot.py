# -*- coding: utf-8 -*-
"""11b_fig3_plot.py — read npz cache only, render Fig 3 (seconds), then programmatic QC.
No xlsx reads here. Cache produced by 11a_fig3_cache.py.
Panels: (a) volcano PR8/MOCK  (b) m6A region distribution  (c) top20 log2FC bidirectional bars.
"""
import os, math, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image
import sys; sys.path.insert(0, r"D:\IAV_m6A\SciData_descriptor\scripts")
import _figstyle as fs
fs.apply_style()

RESDIR = r"D:\IAV_m6A\SciData_descriptor\results"
FIGDIR = os.path.join(RESDIR, "figures")
NPZ = os.path.join(RESDIR, "_fig3_cache.npz")
PDF = os.path.join(FIGDIR, "Figure3_m6a_signal.pdf")
PNG = os.path.join(FIGDIR, "Figure3_m6a_signal.png")
os.makedirs(FIGDIR, exist_ok=True)

RED  = fs.PALETTE["PR8"]    # higher in PR8
BLUE = fs.PALETTE["MOCK"]   # higher in MOCK
GREEN = fs.PALETTE["green"]
GREY = fs.PALETTE["grey_pt"]

d = np.load(NPZ, allow_pickle=True)
log2fc = d["log2fc"]; nlogp = d["nlogp"]
sig = d["sig"]; hi_pr8 = d["hi_pr8"]; hi_mock = d["hi_mock"]
region_keys = list(d["region_keys"]); region_vals = d["region_vals"]
top_syms = list(d["top_syms"]); top_log2fc = d["top_log2fc"]
n_pr8 = int(d["n_pr8"]); n_mock = int(d["n_mock"])

fig = plt.figure(figsize=(13, 4.6))
gs = GridSpec(1, 3, figure=fig, wspace=0.32)

# (a) volcano
axa = fig.add_subplot(gs[0, 0])
ns = ~sig
axa.scatter(log2fc[ns], nlogp[ns], s=2.5, color=GREY, alpha=0.30, rasterized=True, linewidths=0)
axa.scatter(log2fc[hi_pr8], nlogp[hi_pr8], s=4, color=RED, alpha=0.55, rasterized=True,
            linewidths=0, label=f"Higher in PR8 (n={n_pr8:,})")
axa.scatter(log2fc[hi_mock], nlogp[hi_mock], s=4, color=BLUE, alpha=0.55, rasterized=True,
            linewidths=0, label=f"Higher in MOCK (n={n_mock:,})")
axa.axvline(1, ls="--", c="#999999", lw=0.7); axa.axvline(-1, ls="--", c="#999999", lw=0.7)
axa.axhline(-math.log10(0.05), ls="--", c="#999999", lw=0.7)
xmax = np.nanmax(np.abs(log2fc)); ymax = np.nanmax(nlogp)
axa.set_xlim(-xmax * 1.08, xmax * 1.32)
axa.set_ylim(0, ymax * 1.20)
axa.set_xlabel("log$_2$ fold change (PR8 / MOCK)")
axa.set_ylabel("$-$log$_{10}$ $p$-value")
fs.panel_title(axa, "a", "Differential m6A methylation (mRNA)")
axa.legend(loc="upper right", markerscale=2.5, handletextpad=0.3, borderpad=0.4)

# (b) region distribution
axb = fig.add_subplot(gs[0, 1])
vals = [int(v) for v in region_vals]; tot = sum(vals)
bar_colors = [GREEN, BLUE, RED]
bars = axb.bar(region_keys, vals, color=bar_colors, edgecolor="#333333", linewidth=0.6, width=0.62)
for b, v in zip(bars, vals):
    axb.text(b.get_x() + b.get_width() / 2, v, f"{v:,}\n({100*v/tot:.1f}%)",
             ha="center", va="bottom", fontsize=8.5)
axb.set_ylabel("m6A sites on differential mRNAs")
fs.panel_title(axb, "b", "m6A site distribution by mRNA region")
axb.set_ylim(0, max(vals) * 1.22)
axb.tick_params(axis="x", labelsize=9.5)
axb.spines[["top", "right"]].set_visible(False)

# (c) top differential transcripts bidirectional bars
axc = fig.add_subplot(gs[0, 2])
yv = np.arange(len(top_syms))
cols = [RED if v > 0 else BLUE for v in top_log2fc]
axc.barh(yv, top_log2fc, color=cols, edgecolor="#333333", linewidth=0.4, height=0.72)
axc.set_yticks(yv)
fsz = 6.5 if len(top_syms) > 16 else 7.5
axc.set_yticklabels(top_syms, fontsize=fsz)
axc.axvline(0, color="#333333", lw=0.7)
axc.set_xlabel("log$_2$ fold change (PR8 / MOCK)")
fs.panel_title(axc, "c", "Top differential m6A transcripts")
axc.set_ylim(-0.7, len(top_syms) - 0.3)
axc.spines[["top", "right"]].set_visible(False)
# legend for c (color meaning)
from matplotlib.patches import Patch
axc.legend(handles=[Patch(fc=RED, ec="#333333", lw=0.4, label="Higher in PR8"),
                    Patch(fc=BLUE, ec="#333333", lw=0.4, label="Higher in MOCK")],
           loc="lower right", fontsize=7, borderpad=0.4)

fig.suptitle("Figure 3  m6A signal characteristics and differential landscape", **fs.SUPTITLE_KW)
fig.savefig(PDF, dpi=300, bbox_inches="tight")
fig.savefig(PNG, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"SAVED {PDF}")
print(f"SAVED {PNG}")

# ================= programmatic QC =================
print("\n===== QC =====")
im = Image.open(PNG).convert("RGB")
arr = np.asarray(im)
H, W, _ = arr.shape
print(f"PNG size: {W} x {H} px")

# non-white pixel ratio (overall)
nonwhite = np.any(arr < 245, axis=2)
print(f"overall non-white pixel ratio: {100*nonwhite.mean():.1f}%")

# three sub-panels by horizontal thirds; check each has content
thirds = {"(a)": (0, W // 3), "(b)": (W // 3, 2 * W // 3), "(c)": (2 * W // 3, W)}
for name, (x0, x1) in thirds.items():
    sub = nonwhite[:, x0:x1]
    frac = 100 * sub.mean()
    print(f"  panel {name} non-white: {frac:.1f}%  -> {'CONTENT' if frac > 1.5 else 'BLANK!'}")

# (a) legend region: top-right quadrant of panel (a) — count strongly red/blue sig-like pixels
ax0, ax1 = 0, W // 3
lx0, lx1 = ax0 + int((ax1 - ax0) * 0.55), ax1
ly0, ly1 = 0, int(H * 0.32)
leg = arr[ly0:ly1, lx0:lx1, :]
r, g, b = leg[:, :, 0].astype(int), leg[:, :, 1].astype(int), leg[:, :, 2].astype(int)
# strong red scatter (#C44E52-ish) and strong blue scatter (#4C72B0-ish), excluding legend swatch text lines
red_pts  = (r > 150) & (g < 110) & (b < 110)
blue_pts = (b > 130) & (r < 110) & (g < 130)
red_n = int(red_pts.sum()); blue_n = int(blue_pts.sum())
leg_area = leg.shape[0] * leg.shape[1]
overlap_ratio = 100 * (red_n + blue_n) / leg_area
print(f"  (a) legend-zone colored-point pixels: red={red_n} blue={blue_n} "
      f"({overlap_ratio:.2f}% of zone)")
# legend itself contains 2 small swatches; a few hundred px expected. Heavy overlap = thousands.
verdict = "OK (legend not covering significant points)" if overlap_ratio < 1.0 else \
          "WARNING: legend zone has many colored points -> may overlap data"
print(f"  (a) legend overlap verdict: {verdict}")
print("===== QC done =====")
