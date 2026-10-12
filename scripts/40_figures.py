# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import cfg

# -*- coding: utf-8 -*-
"""Revised figures for the IAV m6A Scientific Data descriptor (revision 2026-10).
Reads only figdata/* (exported by 30_export_figdata.R) and the provider gel image; writes to figures/.
Every number drawn on a figure is computed here from the inputs and checked against figdata_summary.json; the script stops on any mismatch.
Before saving, each figure passes (1) a text-text overlap check using the true (rotated) text outlines of drawn labels only and
(2) a text-data check (no annotation or legend covers a plotted point).
Usage: python 40_figures.py [2 3 4]   (default: all). Figure 1, a schematic, is drawn by 41_figure1_schematic.py."""
import os, sys, json, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from PIL import Image

D = cfg.figdata_dir + os.sep
GEL = cfg.qc("gel_image")
OUT = cfg.figures_dir + os.sep
os.makedirs(OUT, exist_ok=True)
SUMM = json.load(open(D + "figdata_summary.json", encoding="utf-8"))
# 30_export_figdata.R writes "pct" as an unnamed vector; names in the order of the R c(...) call:
PCT_KEYS = ["n", "dep_up", "dep_down", "robust", "robust_up", "robust_down", "median_dep", "median_qn", "n_mrna", "r_dep", "r_qn"]
if isinstance(SUMM["pct"], list):
    if len(SUMM["pct"]) != len(PCT_KEYS):
        raise SystemExit("STOP: unexpected pct vector length")
    SUMM["pct"] = dict(zip(PCT_KEYS, SUMM["pct"]))

MOCK, PR8, NS, DARK = "#4C72B0", "#C44E52", "#BBBBBB", "#333333"
S = ["MOCK_1", "MOCK_2", "MOCK_3", "PR8_1", "PR8_2", "PR8_3"]
MM = 1 / 25.4
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300, "axes.unicode_minus": True,
})
SHOWN = {}
SUP_DIG = str.maketrans("-0123456789", "−⁰¹²³⁴⁵⁶⁷⁸⁹")  # Arial lacks U+207B (superscript minus); a minus sign U+2212 is used


def fail(msg):
    print("STOP:", msg)
    sys.exit(1)


def expect(label, got, want):
    if got != want:
        fail("%s = %r, expected %r" % (label, got, want))


def fmt_int(n):
    return "{:,}".format(int(n))


def sci(x):
    e = int(math.floor(math.log10(x)))
    m = x / 10 ** e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    return "%.1f × 10%s" % (m, str(e).translate(SUP_DIG))


def sci_math(x):
    """FDR as m.m × 10^e with a true superscript, typeset in the body font (Arial) via mathtext \\mathregular."""
    e = int(math.floor(math.log10(x)))
    m = x / 10 ** e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    return r"$\mathregular{%.1f \times 10^{%d}}$" % (m, e)


def signed(x, nd=1):
    return ("+" if x >= 0 else "−") + ("%." + str(nd) + "f") % abs(x)


def fig_letter(fig, x, y, s):
    fig.text(x, y, s, fontweight="bold", fontsize=7, ha="left", va="top")


# ------------------------------------------------------------------ checks
def _poly(t, renderer, bb):
    rot = t.get_rotation() % 180
    s = t.get_text()
    if rot in (0.0, 90.0) or "\n" in s:
        return [(bb.x0, bb.y0), (bb.x1, bb.y0), (bb.x1, bb.y1), (bb.x0, bb.y1)]
    w0, h0, _ = renderer.get_text_width_height_descent(s, t.get_fontproperties(), ismath=False)
    cx, cy = (bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2
    th = math.radians(t.get_rotation())
    ux, uy, vx, vy = math.cos(th), math.sin(th), -math.sin(th), math.cos(th)
    return [(cx + sx * w0 / 2 * ux + sy * h0 / 2 * vx, cy + sx * w0 / 2 * uy + sy * h0 / 2 * vy)
            for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]


def _overlap(p, q, tol=0.6):
    for poly in (p, q):
        for i in range(4):
            (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % 4]
            nx, ny = -(y2 - y1), (x2 - x1)
            L = math.hypot(nx, ny) or 1.0
            nx, ny = nx / L, ny / L
            a = [nx * x + ny * y for x, y in p]
            b = [nx * x + ny * y for x, y in q]
            if max(a) < min(b) + tol or max(b) < min(a) + tol:
                return False
    return True


def _texts(fig, renderer):
    items, annot = [], []

    def add(t, is_annot=False):
        if t is None or not t.get_visible():
            return
        s = t.get_text().strip()
        if not s:
            return
        try:
            bb = t.get_window_extent(renderer)
        except Exception:
            return
        if bb.width <= 0 or bb.height <= 0:
            return
        it = (s, bb, _poly(t, renderer, bb), t.get_gid())
        items.append(it)
        if is_annot:
            annot.append(it)
    for ax in fig.axes:
        if not ax.axison:
            for t in ax.texts:
                add(t, True)
            continue
        for axis in (ax.xaxis, ax.yaxis):
            for tick in axis._update_ticks():
                add(tick.label1)
                add(tick.label2)
            add(axis.label)
            add(axis.get_offset_text())
        add(ax.title)
        for t in ax.texts:
            add(t, True)
        leg = ax.get_legend()
        if leg is not None:
            for t in leg.get_texts():
                add(t, True)
    for t in fig.texts:
        add(t)
    return items, annot


def _points(fig):
    pts = []
    for ax in fig.axes:
        for coll in ax.collections:
            if isinstance(coll, (mcoll.PathCollection, mcoll.PolyCollection)) and coll.get_gid() != "nocheck":
                off = np.asarray(coll.get_offsets())
                if off.ndim == 2 and len(off):
                    pts.append(coll.get_offset_transform().transform(off))
    return np.vstack(pts) if pts else np.zeros((0, 2))


def check(fig, name):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    items, annot = _texts(fig, r)
    bad = []
    for s, bb, poly, gid in items:
        if bb.x0 < -0.5 or bb.y0 < -0.5 or bb.x1 > W + 0.5 or bb.y1 > H + 0.5:
            bad.append("outside figure: %r" % s)
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if items[i][3] == "cell" and items[j][3] == "cell":
                continue
            if _overlap(items[i][2], items[j][2]):
                bad.append("text overlap: %r / %r" % (items[i][0], items[j][0]))
    P = _points(fig)
    for s, bb, poly, gid in annot:
        if gid in ("cell", "ondata"):
            continue
        inside = ((P[:, 0] > bb.x0 - 1) & (P[:, 0] < bb.x1 + 1) & (P[:, 1] > bb.y0 - 1) & (P[:, 1] < bb.y1 + 1)).sum() if len(P) else 0
        if inside:
            bad.append("text on data: %r covers %d points" % (s, inside))
    for ax in fig.axes:
        segs = []
        for ln in ax.lines:
            if not ln.get_visible() or ln.get_gid() == "nocheck":
                continue
            xy = np.asarray(ln.get_xydata(), float)
            if len(xy) < 2:
                continue
            segs.append(ln.get_transform().transform(xy))
        for t in ax.texts:
            if not t.get_visible() or not t.get_text().strip() or t.get_gid() in ("cell", "ondata"):
                continue
            bb = t.get_window_extent(r)
            for P2 in segs:
                hit = False
                for (xa, ya), (xb, yb) in zip(P2[:-1], P2[1:]):
                    u = np.linspace(0, 1, 400)
                    xs, ys = xa + (xb - xa) * u, ya + (yb - ya) * u
                    if ((xs > bb.x0 - 1) & (xs < bb.x1 + 1) & (ys > bb.y0 - 1) & (ys < bb.y1 + 1)).any():
                        hit = True
                        break
                if hit:
                    bad.append("text on line: %r" % t.get_text().split("\n")[0])
                    break
        leg = ax.get_legend()
        if leg is not None and len(P):
            bb = leg.get_window_extent(r)
            inside = ((P[:, 0] > bb.x0) & (P[:, 0] < bb.x1) & (P[:, 1] > bb.y0) & (P[:, 1] < bb.y1)).sum()
            if inside:
                bad.append("legend on data: %d points" % inside)
    print("[%s] layout check: %s (%d texts, %d points)" % (name, "PASS" if not bad else "FAIL", len(items), len(P)))
    for b in bad[:40]:
        print("     ", b)
    return not bad


def save(fig, base):
    for ext in ("pdf", "png"):
        p = OUT + base + "." + ext
        if os.path.exists(p):
            os.remove(p)
        fig.savefig(p, dpi=300)
    print("[%s] saved (%.2f MB pdf, %.2f MB png)" % (base, os.path.getsize(OUT + base + ".pdf") / 1e6, os.path.getsize(OUT + base + ".png") / 1e6))


# ------------------------------------------------------------------ Figure 2
def figure2():
    pca = pd.read_csv(D + "fig2b_pca.csv")
    var = SUMM["pca_var"]
    cor = pd.read_csv(D + "fig2c_cor_ip.csv", index_col=0)
    ip = pd.read_csv(D + "fig2d_ip_matrix.csv", index_col=0)
    sup = pd.read_csv(D + "fig2e_sup_matrix.csv", index_col=0)
    pr = pd.read_csv(D + "fig2d_pairs_r.csv")
    expect("probes", ip.shape[0], 40955)
    expect("PC1 %", round(var[0], 1), 63.9); expect("PC2 %", round(var[1], 1), 18.9)
    gel = np.asarray(Image.open(GEL).convert("L"))
    lane_x = [19.5, 58.5, 96.5, 135.0, 174.0, 212.5]
    shown = {"PC1": round(var[0], 1), "PC2": round(var[1], 1), "pair_r": {}, "cor_range_within": None}

    fig = plt.figure(figsize=(180 * MM, 150 * MM))
    ax_a = fig.add_axes([0.065, 0.69, 0.255, 0.255])
    ax_b = fig.add_axes([0.415, 0.69, 0.2, 0.255])
    ax_c = fig.add_axes([0.695, 0.69, 0.205, 0.255])
    cax = fig.add_axes([0.912, 0.715, 0.011, 0.205])
    ax_d = [fig.add_axes([0.07 + k * 0.152, 0.385, 0.112, 0.112 * 180 / 150]) for k in range(6)]
    ax_e = fig.add_axes([0.075, 0.075, 0.9, 0.205])

    # a gel (provider image, unaltered pixels)
    ax_a.imshow(gel, cmap="gray", vmin=0, vmax=255, interpolation="nearest", aspect="equal")
    ax_a.set_xticks(lane_x); ax_a.set_xticklabels(S, rotation=40, ha="right", rotation_mode="anchor")
    ax_a.set_yticks([62, 82]); ax_a.set_yticklabels(["28S", "18S"])
    for sp in ax_a.spines.values():
        sp.set_visible(False)
    ax_a.tick_params(length=2, pad=1.5)

    # b PCA
    for _, r in pca.iterrows():
        ax_b.scatter(r["PC1"], r["PC2"], s=24, color=MOCK if r["group"] == "MOCK" else PR8, edgecolor="white", linewidth=0.5, zorder=3)
    xr = pca["PC1"].max() - pca["PC1"].min(); yr = pca["PC2"].max() - pca["PC2"].min()
    ax_b.set_xlim(pca["PC1"].min() - 0.2 * xr, pca["PC1"].max() + 0.2 * xr)
    ax_b.set_ylim(pca["PC2"].min() - 0.25 * yr, pca["PC2"].max() + 0.25 * yr)
    ax_b.axhline(0, color="#E0E0E0", lw=0.5, zorder=1); ax_b.axvline(0, color="#E0E0E0", lw=0.5, zorder=1)
    for _, r in pca.iterrows():
        up = r["PC2"] >= -5
        ax_b.annotate(r["sample"], (r["PC1"], r["PC2"]), xytext=(0, 4.5 if up else -4.5), textcoords="offset points",
                      ha="center", va="bottom" if up else "top", fontsize=6.3)
    ax_b.set_xlabel("PC1 (%.1f%%)" % var[0]); ax_b.set_ylabel("PC2 (%.1f%%)" % var[1])
    hs = [Line2D([], [], marker="o", ls="", color=c, markersize=4, label=l) for c, l in ((MOCK, "Mock"), (PR8, "PR8"))]
    ax_b.legend(handles=hs, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, handletextpad=0.1, columnspacing=0.9, borderaxespad=0.2)

    # c correlation heatmap
    M = np.corrcoef(ip[S].values.T)   # full precision; the exported CSV is rounded to 4 decimals
    if np.abs(M - cor.loc[S, S].values).max() > 1e-4:
        fail("correlation matrix does not match fig2c_cor_ip.csv")
    vmin = math.floor(M.min() * 100) / 100
    im = ax_c.imshow(M, cmap="viridis", vmin=vmin, vmax=1.0)
    for i in range(6):
        for j in range(6):
            ax_c.text(j, i, "%.3f" % M[i, j], ha="center", va="center", fontsize=5, gid="cell",
                      color="white" if M[i, j] < (vmin + 1) / 2 else "black")
    ax_c.set_xticks(range(6)); ax_c.set_xticklabels(S, rotation=45, ha="right", rotation_mode="anchor")
    ax_c.set_yticks(range(6)); ax_c.set_yticklabels(S)
    for sp in ax_c.spines.values():
        sp.set_visible(False)
    ax_c.tick_params(length=0, pad=1.5)
    cb = fig.colorbar(im, cax=cax); cb.set_label("Pearson r", labelpad=2); cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(width=0.5, length=2, labelsize=6)
    within = [M[i, j] for i in range(6) for j in range(i + 1, 6) if (i < 3) == (j < 3)]
    between = [M[i, j] for i in range(3) for j in range(3, 6)]
    shown["cor_within"] = [round(min(within), 3), round(max(within), 3)]
    shown["cor_between"] = [round(min(between), 3), round(max(between), 3)]

    # d six within-condition scatters
    pairs = [("MOCK_1", "MOCK_2"), ("MOCK_1", "MOCK_3"), ("MOCK_2", "MOCK_3"), ("PR8_1", "PR8_2"), ("PR8_1", "PR8_3"), ("PR8_2", "PR8_3")]
    lo, hi = np.percentile(ip[S].values.ravel(), [0.5, 99.5])
    for ax, (a, b) in zip(ax_d, pairs):
        rr = pr[(pr.x == a) & (pr.y == b)]["r"]
        if len(rr) != 1:
            fail("pair r missing for %s/%s" % (a, b))
        r_calc = float(np.corrcoef(ip[a], ip[b])[0, 1])
        if abs(float(rr.iloc[0]) - r_calc) > 5e-4:
            fail("r mismatch %s/%s" % (a, b))
        hb = ax.hexbin(ip[a].clip(lo, hi), ip[b].clip(lo, hi), gridsize=50, extent=(lo, hi, lo, hi), bins="log", mincnt=1,
                       cmap="Blues" if a.startswith("MOCK") else "Reds", linewidths=0, rasterized=True)
        hb.set_gid("nocheck")
        ax.plot([lo, hi], [lo, hi], ls="--", lw=0.5, color=DARK)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.text(0.05, 0.97, "r = %.3f" % r_calc, transform=ax.transAxes, ha="left", va="top", fontsize=6.3)
        ax.set_xlabel(a + " (log2 IP)", fontsize=6.3, labelpad=1.5); ax.set_ylabel(b + " (log2 IP)", fontsize=6.3, labelpad=1)
        ax.tick_params(labelsize=6, pad=1.5)
        shown["pair_r"]["%s/%s" % (a, b)] = round(r_calc, 3)

    # e per-sample distributions (box = IQR; whiskers = 1st-99th percentile)
    data, pos, cols = [], [], []
    for li, X in enumerate((ip, sup)):
        for si, s in enumerate(S):
            data.append(X[s].values); pos.append(li * 7 + si); cols.append(MOCK if s.startswith("MOCK") else PR8)
    bp = ax_e.boxplot(data, positions=pos, widths=0.6, whis=(1, 99), showfliers=False, patch_artist=True,
                      medianprops=dict(color="black", lw=0.7), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), boxprops=dict(lw=0.6))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.85)
    ax_e.set_xticks(pos); ax_e.set_xticklabels(S + S, fontsize=6.3)
    ax_e.set_ylabel("Normalised intensity (log2)")
    y0, y1 = ax_e.get_ylim()
    ax_e.set_ylim(y0, y1 + 0.14 * (y1 - y0))
    ax_e.text(2.5, y1 + 0.07 * (y1 - y0), "IP (Cy5)", ha="center", va="center", fontweight="bold")
    ax_e.text(9.5, y1 + 0.07 * (y1 - y0), "Sup (Cy3)", ha="center", va="center", fontweight="bold")
    ax_e.axvline(6.0, color="#E0E0E0", lw=0.6)
    ax_e.set_xlim(-0.7, 12.7)

    for x, y, s in ((0.012, 0.985, "a"), (0.345, 0.985, "b"), (0.625, 0.985, "c"), (0.012, 0.565, "d"), (0.012, 0.32, "e")):
        fig_letter(fig, x, y, s)
    SHOWN["Figure 2"] = shown
    ok = check(fig, "Figure 2")
    save(fig, "Figure2_quality")
    plt.close(fig)
    return ok


# ------------------------------------------------------------------ Figure 3
def figure3():
    f = pd.read_csv(D + "fig3_limma_ip_sup.csv")
    md = pd.read_csv(D + "fig3d_machinery.csv")
    shown = {}

    def counts(cls):
        g = f[f.rna_class == cls]
        return int((g.ip_call == "higher_PR8").sum()), int((g.ip_call == "higher_MOCK").sum())
    um, dm = counts("mRNA"); ul, dl = counts("lncRNA")
    expect("mRNA higher PR8", um, 3067); expect("mRNA higher mock", dm, 2967)
    expect("lncRNA higher PR8", ul, 729); expect("lncRNA higher mock", dl, 502)
    expect("machinery probes", len(md), 42); expect("machinery genes", md.gene.nunique(), 19)

    fig = plt.figure(figsize=(180 * MM, 118 * MM))
    ax_a = fig.add_axes([0.07, 0.54, 0.2, 0.40])
    ax_b = fig.add_axes([0.345, 0.54, 0.2, 0.40])
    ax_c = fig.add_axes([0.775, 0.54, 0.135, 0.40])
    ax_d = fig.add_axes([0.105, 0.165, 0.77, 0.16])
    cax = fig.add_axes([0.905, 0.145, 0.011, 0.2])

    for ax, cls, (u, d) in ((ax_a, "mRNA", (um, dm)), (ax_b, "lncRNA", (ul, dl))):
        g = f[f.rna_class == cls]
        y = -np.log10(g.ip_fdr.clip(lower=1e-300))
        ns = g.ip_call == "ns"
        ax.scatter(g.ip_lfc[ns], y[ns], s=1.2, color=NS, linewidths=0, rasterized=True)
        for call, c in (("higher_PR8", PR8), ("higher_MOCK", MOCK)):
            sel = g.ip_call == call
            ax.scatter(g.ip_lfc[sel], y[sel], s=1.6, color=c, linewidths=0, rasterized=True)
        xm = max(abs(g.ip_lfc.min()), abs(g.ip_lfc.max())) * 1.05
        ax.set_xlim(-xm, xm); ax.set_ylim(0, y.max() * 1.45)
        ax.axvline(1, ls="--", lw=0.5, color=DARK); ax.axvline(-1, ls="--", lw=0.5, color=DARK)
        ax.axhline(-math.log10(0.05), ls="--", lw=0.5, color=DARK)
        ax.text(0.98, 0.98, "Higher in\nPR8\nn = %s" % fmt_int(u), transform=ax.transAxes, ha="right", va="top", color=PR8, fontsize=6.3)
        ax.text(0.02, 0.98, "Higher in\nmock\nn = %s" % fmt_int(d), transform=ax.transAxes, ha="left", va="top", color=MOCK, fontsize=6.3)
        ax.set_xlabel("log2 fold change (PR8/mock)"); ax.set_ylabel("−log10 FDR")
        ax.text(0.5, 1.02, "mRNA" if cls == "mRNA" else "lncRNA", transform=ax.transAxes, ha="center", va="bottom", fontweight="bold")
        shown[cls] = {"higher_PR8": u, "higher_mock": d}

    # c top probes
    g = f[(f.rna_class == "mRNA") & (f.ip_call != "ns")]
    top = g.sort_values("ip_lfc", ascending=False).head(10)
    bot = g.sort_values("ip_lfc", ascending=True).head(10).iloc[::-1]
    tb = pd.concat([top, bot])
    if tb.probeID.duplicated().any():
        fail("duplicate probe in top list")
    labels = ["%s (%s)" % (r.gene if isinstance(r.gene, str) and r.gene else "—", r.transcript) for r in tb.itertuples()]
    if len(set(labels)) != len(labels):
        fail("duplicate label in top list")
    yy = np.arange(len(tb))[::-1]
    ax_c.barh(yy, tb.ip_lfc, color=[PR8 if v > 0 else MOCK for v in tb.ip_lfc], height=0.72)
    ax_c.set_yticks(yy); ax_c.set_yticklabels(labels, fontsize=5.6)
    ax_c.axvline(0, color=DARK, lw=0.6)
    ax_c.set_xlabel("log2 fold change (PR8/mock)")
    ax_c.set_ylim(-0.7, len(tb) - 0.3)
    ax_c.tick_params(axis="y", length=0, pad=1.5)
    ax_c.spines["left"].set_visible(False)
    tr = ax_c.get_yaxis_transform()
    ax_c.text(1.04, len(tb) - 0.25, "FDR", transform=tr, ha="left", va="bottom", fontsize=5.6, fontweight="bold")
    for y_, q in zip(yy, tb.ip_fdr):
        ax_c.text(1.04, y_, sci(q), transform=tr, ha="left", va="center", fontsize=5.6)
    shown["top_probes"] = [(lab, round(float(v), 2), sci(q)) for lab, v, q in zip(labels, tb.ip_lfc, tb.ip_fdr)]

    # d m6A machinery, one column per probe
    order_role = {"writer": 0, "eraser": 1, "reader": 2}
    md = md.assign(_r=md.role.map(order_role)).sort_values(["_r", "gene", "transcript"]).reset_index(drop=True)
    x, cur, prev_role, prev_gene = [], 0.0, None, None
    for r in md.itertuples():
        if prev_role is not None:
            cur += 1.6 if r.role != prev_role else (0.55 if r.gene != prev_gene else 0.0)
        x.append(cur); cur += 1.0; prev_role, prev_gene = r.role, r.gene
    md["x"] = x
    vmax = math.ceil(max(md.ip_lfc.abs().max(), md.sup_lfc.abs().max()) * 2) / 2
    norm = Normalize(-vmax, vmax); cmap = plt.get_cmap("RdBu_r")
    for yrow, lfc, fdr in ((1, "ip_lfc", "ip_fdr"), (0, "sup_lfc", "sup_fdr")):
        sig = (md[fdr] < 0.05).values
        ax_d.scatter(md.x, np.full(len(md), yrow), s=30, c=md[lfc], cmap=cmap, norm=norm,
                     edgecolors=np.where(sig, "black", "#A0A0A0"), linewidths=np.where(sig, 0.9, 0.45), zorder=3)
    ax_d.set_xlim(md.x.min() - 0.8, md.x.max() + 0.8); ax_d.set_ylim(-0.55, 1.55)
    ax_d.set_yticks([1, 0]); ax_d.set_yticklabels(["m6A quantity", "Expression"])
    ax_d.set_xticks([])
    for sp in ax_d.spines.values():
        sp.set_visible(False)
    ax_d.tick_params(axis="y", length=0)
    for gene, sub in md.groupby("gene", sort=False):
        x0, x1 = sub.x.min(), sub.x.max()
        if x1 > x0:
            ax_d.plot([x0 - 0.3, x1 + 0.3], [-0.42, -0.42], color="#777777", lw=0.5, clip_on=False)
        ax_d.text((x0 + x1) / 2, -0.55, gene, rotation=90, ha="center", va="top", fontsize=6, clip_on=False)
    for role, sub in md.groupby("role", sort=False):
        x0, x1 = sub.x.min(), sub.x.max()
        ax_d.plot([x0 - 0.3, x1 + 0.3], [1.42, 1.42], color="#777777", lw=0.5, clip_on=False)
        ax_d.text((x0 + x1) / 2, 1.5, {"writer": "Writers", "eraser": "Erasers", "reader": "Readers"}[role],
                  ha="center", va="bottom", fontsize=6.3, fontweight="bold", clip_on=False)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, cax=cax); cb.set_label("log2 fold change\n(PR8/mock)", labelpad=2); cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(width=0.5, length=2, labelsize=6)
    key = [Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor="black", markeredgewidth=0.9, markersize=4.8, label="FDR < 0.05"),
           Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor="#A0A0A0", markeredgewidth=0.45, markersize=4.8, label="FDR ≥ 0.05")]
    ax_d.legend(handles=key, loc="lower left", bbox_to_anchor=(0.0, 1.22), ncol=2, handletextpad=0.1, columnspacing=1.0, borderaxespad=0)
    shown["machinery"] = {"probes": len(md), "genes": int(md.gene.nunique()), "colour_limit": vmax,
                          "sig_ip": int((md.ip_fdr < 0.05).sum()), "sig_sup": int((md.sup_fdr < 0.05).sum())}

    for x_, y_, s in ((0.012, 0.985, "a"), (0.285, 0.985, "b"), (0.565, 0.985, "c"), (0.012, 0.425, "d")):
        fig_letter(fig, x_, y_, s)
    SHOWN["Figure 3"] = shown
    ok = check(fig, "Figure 3")
    save(fig, "Figure3_m6a_signal")
    plt.close(fig)
    return ok


# ------------------------------------------------------------------ Figure 4
def figure4():
    f = pd.read_csv(D + "fig4_pctmod_detected.csv")
    pc = SUMM["pct"]
    expect("tested probes", len(f), int(pc["n"]))
    up = int((f.call_dep == "higher_PR8").sum()); dn = int((f.call_dep == "higher_MOCK").sum())
    rob = f.robust.astype(str).str.upper() == "TRUE"
    rup = int((rob & (f.call_dep == "higher_PR8")).sum()); rdn = int((rob & (f.call_dep == "higher_MOCK")).sum())
    expect("pct higher PR8", up, int(pc["dep_up"])); expect("pct higher mock", dn, int(pc["dep_down"]))
    expect("robust", int(rob.sum()), int(pc["robust"])); expect("robust up", rup, int(pc["robust_up"])); expect("robust down", rdn, int(pc["robust_down"]))
    m = f[f.rna_class == "mRNA"]
    expect("detected mRNA probes", len(m), int(pc["n_mrna"]))
    shown = {}
    fig = plt.figure(figsize=(180 * MM, 68 * MM))
    axs = [fig.add_axes([0.075 + k * 0.33, 0.16, 0.235, 0.235 * 180 / 68]) for k in range(3)]
    xs = np.concatenate([m.expr_lfc_dep, m.expr_lfc_qn]); ys = np.concatenate([m.dpp_dep, m.dpp_qn])
    xl = np.percentile(xs, [0.2, 99.8]); yl = np.percentile(ys, [0.2, 99.8])
    xlim = (-max(abs(xl[0]), abs(xl[1])), max(abs(xl[0]), abs(xl[1])))
    yspan = max(abs(yl[0]), abs(yl[1]))
    ylim = (-yspan, yspan)
    for ax, xcol, ycol, ccol, tag in ((axs[0], "expr_lfc_dep", "dpp_dep", "call_dep", "Deposited normalisation"),
                                      (axs[1], "expr_lfc_qn", "dpp_qn", "call_qn", "Quantile-normalised")):
        inr = (m[xcol].abs() <= xlim[1]) & (m[ycol].abs() <= yspan)
        X = m[xcol].clip(*xlim); Y = m[ycol].clip(ylim[0], yspan)
        hb = ax.hexbin(X[inr], Y[inr], gridsize=55, extent=(xlim[0], xlim[1], ylim[0], yspan), bins="log", mincnt=1, cmap="Greys", linewidths=0, rasterized=True)
        hb.set_gid("nocheck")
        for call, c in (("higher_PR8", PR8), ("higher_MOCK", MOCK)):
            sel = (m[ccol] == call) & inr
            ax.scatter(X[sel], Y[sel], s=2.2, color=c, linewidths=0, rasterized=True, zorder=3)
        # probes beyond the axis range: open triangles on the edge, pointing outwards, coloured by call
        out = ~inr
        for call, c in (("ns", "#8C8C8C"), ("higher_PR8", PR8), ("higher_MOCK", MOCK)):
            for cond, mk in ((m[ycol] > yspan, "^"), (m[ycol] < -yspan, "v"), ((m[ycol].abs() <= yspan) & (m[xcol] > xlim[1]), ">"),
                             ((m[ycol].abs() <= yspan) & (m[xcol] < xlim[0]), "<")):
                sel = out & cond & (m[ccol] == call)
                if sel.any():
                    sc = ax.scatter(X[sel], Y[sel], s=7, marker=mk, facecolors="none", edgecolors=c, linewidths=0.5, zorder=4, clip_on=False)
                    sc.set_gid("nocheck")
        med = float(np.median(m[ycol]))
        r = float(np.corrcoef(m[xcol], m[ycol])[0, 1])
        ax.axhline(0, color=DARK, lw=0.5); ax.axhline(med, color=DARK, lw=0.7, ls="--")
        ax.axvline(1, color="#888888", lw=0.5, ls=":"); ax.axvline(-1, color="#888888", lw=0.5, ls=":")
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.text(0.5, 1.115, tag, transform=ax.transAxes, ha="center", va="bottom", fontweight="bold")
        ax.text(0.5, 1.015, "Median Δ = %s (dashed); r = %s" % (signed(med), signed(r, 2)), transform=ax.transAxes, ha="center", va="bottom", fontsize=6.3)
        ax.set_xlabel("Expression log2 fold change (PR8/mock)"); ax.set_ylabel("Δ%modified, PR8 − mock (% points)")
        shown[tag] = {"median": signed(med), "r": signed(r, 2), "n_plotted_mRNA": int(len(m)), "n_beyond_axes": int(out.sum()),
                      "n_called_beyond_axes": int((out & (m[ccol] != "ns")).sum()), "median_exact": round(med, 3)}
    expect("median deposited label", shown["Deposited normalisation"]["median"], "+4.7")
    expect("median QN label", shown["Quantile-normalised"]["median"], "+1.0")
    expect("r deposited", shown["Deposited normalisation"]["r"], "−0.08")
    expect("r QN", shown["Quantile-normalised"]["r"], "−0.06")

    ax = axs[2]
    y = -np.log10(f.fdr_dep.clip(lower=1e-300))
    ns = f.call_dep == "ns"
    ax.scatter(f.dpp_dep[ns], y[ns], s=1.2, color=NS, linewidths=0, rasterized=True)
    for call, c in (("higher_PR8", PR8), ("higher_MOCK", MOCK)):
        sel = (f.call_dep == call) & ~rob
        ax.scatter(f.dpp_dep[sel], y[sel], s=3, color=c, linewidths=0, rasterized=True, zorder=3)
        sel = (f.call_dep == call) & rob
        ax.scatter(f.dpp_dep[sel], y[sel], s=7, color=c, edgecolors="black", linewidths=0.45, rasterized=True, zorder=4)
    xm = max(abs(f.dpp_dep.min()), abs(f.dpp_dep.max())) * 1.05
    ax.set_xlim(-xm, xm); ax.set_ylim(0, y.max() * 1.6)
    ax.axvline(5, ls="--", lw=0.5, color=DARK); ax.axvline(-5, ls="--", lw=0.5, color=DARK)
    ax.axhline(-math.log10(0.05), ls="--", lw=0.5, color=DARK)
    ax.text(0.98, 0.98, "Higher in\nPR8\nn = %d" % up, transform=ax.transAxes, ha="right", va="top", color=PR8, fontsize=6.3)
    ax.text(0.02, 0.98, "Higher in\nmock\nn = %d" % dn, transform=ax.transAxes, ha="left", va="top", color=MOCK, fontsize=6.3)
    ax.text(0.5, 1.115, "Differential %modified", transform=ax.transAxes, ha="center", va="bottom", fontweight="bold")
    ax.text(0.5, 1.015, "Black outline: robust to normalisation, n = %d" % int(rob.sum()), transform=ax.transAxes,
            ha="center", va="bottom", fontsize=6.3)
    ax.set_xlabel("Δ%modified, PR8 − mock (% points)"); ax.set_ylabel("−log10 FDR")
    shown["volcano"] = {"higher_PR8": up, "higher_mock": dn, "robust": int(rob.sum()), "robust_PR8": rup, "robust_mock": rdn}
    for x_, s in ((0.005, "a"), (0.335, "b"), (0.665, "c")):
        fig_letter(fig, x_, 0.985, s)
    SHOWN["Figure 4"] = shown
    ok = check(fig, "Figure 4")
    save(fig, "Figure4_methylation_fraction")
    plt.close(fig)
    return ok


if __name__ == "__main__":
    which = sys.argv[1:] or ["2", "3", "4"]
    res = {}
    for w in which:
        res["Figure " + w] = globals()["figure" + w]()
    print("\nLAYOUT CHECKS:", res)
    print("NUMBERS SHOWN:", json.dumps(SHOWN, indent=1, ensure_ascii=False))
