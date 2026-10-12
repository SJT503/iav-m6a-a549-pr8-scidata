# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import cfg

# -*- coding: utf-8 -*-
"""Figure 1 -- study design and data workflow.

A schematic drawn entirely in code: the icons are built from matplotlib primitives and every label is text set
here, so the figure contains no raster or generated artwork and can be edited in any vector editor. The one number
it shows (the probe count) is read from figdata/fig2d_ip_matrix.csv and checked. Typeset in Arial 7 pt with 0.6 pt
outlines, as Figures 2-4; MOCK and PR8 use the colours of Figures 2-4. Coordinates are millimetres from the top-left
corner of the figure (y increases downwards).

Before saving, the layout passes (1) text/text, text/icon, icon/icon and arrow overlap checks, (2) every label and
icon lies inside its panel, (3) every icon is drawn inside its declared box. If the font that is available is wider
than Arial, the text size is reduced in 0.2 pt steps until the checks pass; the script stops if they cannot.

Needs figdata/ (stage 30). Writes figures/Figure1_design.pdf and .png. Stage 40 draws Figures 2-4.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import findfont, FontProperties
from matplotlib.patches import (FancyBboxPatch, Polygon, Circle, Ellipse, Rectangle, Arc, FancyArrowPatch,
                                PathPatch)
from matplotlib.path import Path

MM = 1 / 25.4
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
    "axes.unicode_minus": True,
})

# ---------------------------------------------------------------- style
PT_STEPS = [7.0, 6.8, 6.6, 6.4, 6.2, 6.0, 5.8, 5.6]
EDGE, LW = "#4D4D4D", 0.6
TXT, HEADC = "#262626", "#1A1A1A"
MOCK, PR8 = "#4C72B0", "#C44E52"             # as in Figures 2-4
IPC, SUPC, CAL, SPK = "#8172B2", "#55A868", "#CCB974", "#DD8452"   # IP, Sup, calibration and m6A spike-ins
RNA, AB = "#7A7A7A", "#3B4F6B"
GLASS, GLASS2 = "#F7FAFC", "#E6EBF0"
TINT = ["#FBF1F1", "#FCF7EE", "#F1F0F8", "#EDF6F2"]

# ---------------------------------------------------------------- geometry (mm)
W = 180.0
MARGIN, GUT = 1.5, 6.0
PW = [35.0, 45.0, 40.0, 39.0]
PX = [MARGIN + sum(PW[:i]) + i * GUT for i in range(4)]
assert abs(PX[3] + PW[3] + MARGIN - W) < 1e-6
HEAD_Y, BOX_TOP = 2.9, 5.4
PAD_T, PAD_B, BOT = 3.0, 3.0, 1.2
TXT_GAP, GAP_MIN, GAP_MAX = 0.9, 2.8, 7.0

SZ = dict(cells=(22.65, 7.55), flask=(18.0, 8.0), tube=(5.0, 11.5), gel=(13.0, 9.0),
          drop=(4.6, 9.9), ip=(22.0, 7.8), split=(10.6, 13.6), dyes=(8.0, 6.9),
          slide=(32.0, 12.0), scanner=(18.0, 11.0), db=(10.0, 12.0), files=(14.6, 10.8),
          matrix=(15.6, 10.0), tables=(13.0, 10.0))


def fail(msg):
    print("STOP:", msg)
    sys.exit(1)


def expect(label, got, want):
    if got != want:
        fail("%s = %r, expected %r" % (label, got, want))


def fmt_int(n):
    return "{:,}".format(int(n))


# ---------------------------------------------------------------- drawing surface
class Canvas:
    def __init__(self, w, h):
        self.w, self.h = w, h
        self.fig = plt.figure(figsize=(w * MM, h * MM))
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, w)
        self.ax.set_ylim(h, 0)
        self.ax.axis("off")
        self.cur = []                       # artists of the icon being drawn

    def patch(self, p):
        self.ax.add_patch(p)
        self.cur.append(p)
        return p

    def line(self, xs, ys, **kw):
        ln, = self.ax.plot(xs, ys, gid="nocheck", **kw)
        self.cur.append(ln)
        return ln


def P(x0, y0, pts):
    return [(x0 + a, y0 + b) for a, b in pts]


# ---------------------------------------------------------------- icons
def f_cells(c, x0, y0):                     # A549 monolayer: cobblestone epithelial cells
    cw, ch, g = 5.4, 3.6, 0.35
    for r, (off, n) in enumerate(((cw / 2 + g / 2, 3), (0.0, 4))):
        for i in range(n):
            x, y = x0 + off + i * (cw + g), y0 + r * (ch + g)
            c.patch(FancyBboxPatch((x, y), cw, ch, boxstyle="round,pad=0,rounding_size=1.0",
                                   fc="#FFFFFF", ec="#8E7471", lw=LW, zorder=4))
            c.patch(Ellipse((x + cw / 2, y + ch / 2), 1.9, 1.25, fc="#C59C98", ec="none", zorder=5))


def virion(c, cx, cy, r, col):
    for k in range(12):
        a = 2 * np.pi * k / 12
        c.line([cx + r * np.cos(a), cx + (r + 0.42) * np.cos(a)],
               [cy + r * np.sin(a), cy + (r + 0.42) * np.sin(a)],
               color=col, lw=0.7, solid_capstyle="round", zorder=7)
    c.patch(Circle((cx, cy), r, fc=col, ec="none", zorder=7))


def f_flask(c, x0, y0, col, virus=False):   # culture flask lying flat, medium tinted by condition
    c.patch(Polygon(P(x0, y0, [(18, 0), (5.2, 0), (3.3, 2.4), (3.3, 5.6), (5.2, 8), (18, 8)]),
                    closed=True, fc=GLASS, ec=EDGE, lw=LW, joinstyle="round", zorder=4))
    c.patch(Polygon(P(x0, y0, [(3.45, 5.3), (17.85, 5.3), (17.85, 7.85), (5.25, 7.85), (3.45, 5.75)]),
                    closed=True, fc=col, alpha=0.30, ec="none", zorder=5))
    c.line([x0 + 3.45, x0 + 17.85], [y0 + 5.3, y0 + 5.3], color=col, lw=0.6, zorder=5)
    c.patch(Polygon(P(x0, y0, [(3.4, 2.95), (1.7, 2.95), (1.7, 5.05), (3.4, 5.05)]), closed=False,
                    fc=GLASS, ec=EDGE, lw=LW, joinstyle="round", zorder=4))
    c.patch(FancyBboxPatch((x0 + 0.05, y0 + 2.45), 1.85, 3.1, boxstyle="round,pad=0,rounding_size=0.45",
                           fc=col, ec=EDGE, lw=LW, zorder=5))
    if virus:
        for vx, vy in ((8.6, 2.5), (12.2, 1.9), (15.5, 2.9)):
            virion(c, x0 + vx, y0 + vy, 0.8, col)


def f_tube(c, x0, y0, s=1.0, frac=None):    # microcentrifuge tube; frac = colour of the liquid
    Q = lambda pts: [(x0 + a * s, y0 + b * s) for a, b in pts]
    c.patch(Polygon(Q([(0.45, 1.9), (4.55, 1.9), (4.55, 6.8), (2.5, 11.5), (0.45, 6.8)]), closed=True,
                    fc=GLASS, ec=EDGE, lw=LW, joinstyle="round", zorder=4))
    if frac is None:
        c.patch(Polygon(Q([(1.05, 8.0), (3.95, 8.0), (2.5, 11.15)]), closed=True,
                        fc="#DCE3EC", ec="none", zorder=5))
    else:
        c.patch(Polygon(Q([(0.6, 4.2), (4.4, 4.2), (4.4, 6.85), (2.5, 11.15), (0.6, 6.85)]), closed=True,
                        fc=frac, alpha=0.75, ec="none", zorder=5))
    c.patch(Rectangle((x0, y0 + 1.1 * s), 5.0 * s, 0.8 * s, fc=GLASS2, ec=EDGE, lw=LW, zorder=5))
    c.patch(FancyBboxPatch((x0 + 0.3 * s, y0), 4.4 * s, 1.1 * s,
                           boxstyle="round,pad=0,rounding_size=%.3f" % (0.35 * s),
                           fc=GLASS2, ec=EDGE, lw=LW, zorder=5))


def f_gel(c, x0, y0):                       # denaturing gel, six lanes: 28S, 18S and low-molecular-weight bands
    c.patch(FancyBboxPatch((x0, y0), 13.0, 9.0, boxstyle="round,pad=0,rounding_size=0.6",
                           fc="#2F3542", ec=EDGE, lw=LW, zorder=4))
    for i in range(6):
        cx = x0 + 1.75 + i * 1.9
        c.patch(Rectangle((cx - 0.6, y0 + 0.75), 1.2, 0.4, fc="#191C23", ec="none", zorder=5))
        for yy, hh, a in ((2.9, 0.55, 0.95), (4.7, 0.48, 0.80), (7.3, 0.40, 0.35)):
            c.patch(FancyBboxPatch((cx - 0.6, y0 + yy), 1.2, hh, boxstyle="round,pad=0,rounding_size=0.15",
                                   fc="#F2F2F2", ec="none", alpha=a, zorder=5))


def f_drop(c, x0, y0, col):                 # pipette tip and a falling drop = a spike-in is added
    c.patch(Polygon(P(x0, y0, [(1.25, 0.0), (3.35, 0.0), (3.35, 0.9), (2.55, 4.6), (2.05, 4.6), (1.25, 0.9)]),
                    closed=True, fc=GLASS2, ec=EDGE, lw=LW, joinstyle="round", zorder=5))
    cx, cy, r, ty = x0 + 2.3, y0 + 8.1, 1.7, y0 + 5.3
    k = 0.5523
    verts = [(cx, ty),
             (cx + 0.3 * r, ty + 0.55 * r), (cx + r, cy - 0.75 * r), (cx + r, cy),
             (cx + r, cy + k * r), (cx + k * r, cy + r), (cx, cy + r),
             (cx - k * r, cy + r), (cx - r, cy + k * r), (cx - r, cy),
             (cx - r, cy - 0.75 * r), (cx - 0.3 * r, ty + 0.55 * r), (cx, ty),
             (cx, ty)]
    codes = [Path.MOVETO] + [Path.CURVE4] * 12 + [Path.CLOSEPOLY]
    c.patch(PathPatch(Path(verts, codes), fc=col, ec=EDGE, lw=LW, zorder=5))
    c.patch(Ellipse((cx - 0.55, cy - 0.45), 0.7, 1.0, angle=25, fc="#FFFFFF", ec="none", alpha=0.55, zorder=6))


def f_ip(c, x0, y0):                        # RNA carrying m6A marks, captured by an antibody on a bead
    yc, A, lam, xs = 6.5, 0.7, 5.0, 4.0
    xx = np.linspace(0.6, 21.4, 300)
    c.line(x0 + xx, y0 + yc + A * np.sin(2 * np.pi * (xx - xs) / lam), color=RNA, lw=1.0,
           solid_capstyle="round", zorder=5)
    bx, by, br = x0 + 14.0, y0 + 1.95, 1.85
    hy = y0 + 4.85
    c.line([bx, bx], [by + br, hy], color=AB, lw=1.5, solid_capstyle="round", zorder=6)
    for dx in (-1.05, 1.05):
        c.line([bx, bx + dx], [hy, y0 + 5.95], color=AB, lw=1.5, solid_capstyle="round", zorder=6)
    c.patch(Circle((bx, by), br, fc="#BDBDBD", ec=EDGE, lw=LW, zorder=6))
    c.patch(Circle((bx - 0.6, by - 0.6), 0.55, fc="#FFFFFF", ec="none", alpha=0.7, zorder=7))
    for mx in (4.0, 14.0):
        c.patch(Circle((x0 + mx, y0 + yc), 0.55, fc=IPC, ec="#FFFFFF", lw=0.4, zorder=7))


def f_split(c, x0, y0):                     # IP and Sup fractions; returns the two label anchors
    ym, yu, yl, gc, s = 7.0, 3.4, 10.6, "#8C8C8C", 0.5
    c.line([x0, x0 + 2.4], [y0 + ym, y0 + ym], color=gc, lw=0.8, solid_capstyle="round", zorder=5)
    for yt in (yu, yl):
        c.line([x0 + 2.4, x0 + 4.6, x0 + 6.0], [y0 + ym, y0 + yt, y0 + yt], color=gc, lw=0.8,
               solid_joinstyle="round", solid_capstyle="round", zorder=5)
        c.patch(FancyArrowPatch((x0 + 5.9, y0 + yt), (x0 + 7.6, y0 + yt), arrowstyle="-|>", mutation_scale=6,
                                lw=0.8, color=gc, shrinkA=0, shrinkB=0, zorder=5))
    f_tube(c, x0 + 8.0, y0 + yu - 11.5 * s / 2, s, IPC)
    f_tube(c, x0 + 8.0, y0 + yl - 11.5 * s / 2, s, SUPC)
    xr = x0 + 8.0 + 5.0 * s
    return (xr, y0 + yu), (xr, y0 + yl)


def f_dyes(c, x0, y0):                      # the IP and Sup fractions after labelling: one tube per dye
    s = 0.6
    f_tube(c, x0, y0, s, IPC)
    f_tube(c, x0 + 5.0 * s + 2.0, y0, s, SUPC)


def f_slide(c, x0, y0):                     # 8x60K slide: frosted end and eight sub-arrays
    c.patch(FancyBboxPatch((x0, y0), 32.0, 12.0, boxstyle="round,pad=0,rounding_size=0.7",
                           fc="#F7FAFD", ec=EDGE, lw=LW, zorder=4))
    c.patch(Rectangle((x0 + 0.35, y0 + 0.35), 4.6, 11.3, fc="#DADFE5", ec="none", zorder=5))
    rng = np.random.default_rng(11)
    cols = [IPC, SUPC, "#C8C8C8"]
    for r in range(2):
        for k in range(4):
            sx, sy = x0 + 6.35 + k * 6.3, y0 + 1.2 + r * 5.2
            c.patch(Rectangle((sx, sy), 5.6, 4.4, fc="#FFFFFF", ec="#9A9A9A", lw=0.4, zorder=5))
            for i in range(6):
                for j in range(4):
                    c.patch(Circle((sx + 0.55 + i * 0.9, sy + 0.62 + j * 1.05), 0.27,
                                   fc=cols[rng.integers(0, 3)], ec="none", zorder=6))


def f_scanner(c, x0, y0):
    c.patch(Polygon(P(x0, y0, [(1.2, 4.2), (3.4, 0.6), (16.8, 0.6), (16.8, 4.2)]), closed=True,
                    fc="#D5DAE0", ec=EDGE, lw=LW, joinstyle="round", zorder=4))
    c.patch(FancyBboxPatch((x0, y0 + 4.2), 18.0, 6.6, boxstyle="round,pad=0,rounding_size=0.8",
                           fc="#E7EAEE", ec=EDGE, lw=LW, zorder=5))
    c.patch(FancyBboxPatch((x0 + 2.6, y0 + 6.7), 9.6, 1.1, boxstyle="round,pad=0,rounding_size=0.3",
                           fc="#5B616B", ec="none", zorder=6))
    c.patch(Circle((x0 + 15.0, y0 + 7.25), 0.55, fc=SUPC, ec="none", zorder=6))


def f_db(c, x0, y0):                        # database cylinder
    w, eh = 10.0, 3.0
    cx, top, bot = x0 + w / 2, y0 + eh / 2, y0 + 12.0 - eh / 2
    fc, ft = "#D8E2EE", "#EEF3F8"
    c.patch(Ellipse((cx, bot), w, eh, fc=fc, ec=EDGE, lw=LW, zorder=4))
    c.patch(Rectangle((x0, top), w, bot - top, fc=fc, ec="none", zorder=4))
    c.line([x0, x0], [top, bot], color=EDGE, lw=LW, zorder=5)
    c.line([x0 + w, x0 + w], [top, bot], color=EDGE, lw=LW, zorder=5)
    for yy in (top + (bot - top) / 3, top + 2 * (bot - top) / 3):
        c.patch(Arc((cx, yy), w, eh, theta1=0, theta2=180, ec=EDGE, lw=LW, zorder=5))
    c.patch(Ellipse((cx, top), w, eh, fc=ft, ec=EDGE, lw=LW, zorder=5))


def f_doc(c, x, y, w=4.2, h=5.0):
    f = 1.1
    c.patch(Polygon([(x, y), (x + w - f, y), (x + w, y + f), (x + w, y + h), (x, y + h)], closed=True,
                    fc="#FFFFFF", ec=EDGE, lw=LW, joinstyle="round", zorder=5))
    c.patch(Polygon([(x + w - f, y), (x + w - f, y + f), (x + w, y + f)], closed=True,
                    fc="#E2E2E2", ec=EDGE, lw=0.45, joinstyle="round", zorder=6))
    for kk in range(3):
        yy = y + 2.0 + kk * 0.95
        c.line([x + 0.7, x + w - 0.7], [yy, yy], color="#B3B3B3", lw=0.5, zorder=6)


def f_files(c, x0, y0):                     # the six raw files
    for r in range(2):
        for k in range(3):
            f_doc(c, x0 + k * 5.2, y0 + r * 5.8)


def f_matrix(c, x0, y0):                    # 12 channels: 6 IP (violet) and 6 Sup (green)
    c.patch(FancyBboxPatch((x0, y0), 15.6, 10.0, boxstyle="round,pad=0,rounding_size=0.5",
                           fc="#FFFFFF", ec=EDGE, lw=LW, zorder=4))
    c.patch(Rectangle((x0 + 0.3, y0 + 0.3), 15.0, 1.3, fc="#E3E3E3", ec="none", zorder=5))
    rng = np.random.default_rng(5)
    cw, chh = 15.0 / 12, 7.6 / 6
    for k in range(12):
        for r in range(6):
            c.patch(Rectangle((x0 + 0.3 + k * cw, y0 + 2.0 + r * chh), cw - 0.12, chh - 0.12,
                              fc=IPC if k < 6 else SUPC, alpha=float(rng.uniform(0.2, 0.85)),
                              ec="none", zorder=5))
    c.line([x0 + 0.3 + 6 * cw - 0.06] * 2, [y0 + 0.3, y0 + 9.7], color=EDGE, lw=0.5, zorder=6)


def f_tables(c, x0, y0):                    # the provider tables (two sheets)
    for dx, dy, z in ((2.2, 0.0, 4.0), (0.0, 1.8, 5.0)):
        x, y = x0 + dx, y0 + dy
        c.patch(Rectangle((x, y), 10.8, 8.2, fc="#FFFFFF", ec=EDGE, lw=LW, zorder=z))
        c.patch(Rectangle((x + 0.3, y + 0.3), 10.2, 1.3, fc="#E3E3E3", ec="none", zorder=z + 0.1))
        for kk in range(1, 4):
            c.line([x + kk * 10.8 / 4] * 2, [y + 1.6, y + 7.9], color="#BDBDBD", lw=0.45, zorder=z + 0.2)
        for kk in range(1, 5):
            c.line([x + 0.3, x + 10.5], [y + 1.6 + kk * 1.26] * 2, color="#BDBDBD", lw=0.45, zorder=z + 0.2)


FN = dict(cells=f_cells, flask=f_flask, tube=f_tube, gel=f_gel, drop=f_drop, ip=f_ip, dyes=f_dyes,
          slide=f_slide, scanner=f_scanner, db=f_db, files=f_files, matrix=f_matrix, tables=f_tables)


# ---------------------------------------------------------------- content
def columns(n_probes):
    return [
        dict(head="Cell model", rows=[
            dict(icons=[("cells", {})], lines=["A549 cells", "(ATCC CCL-185)"]),
            dict(icons=[("flask", dict(col=MOCK))], lines=["Mock-infected", "n = 3"]),
            dict(icons=[("flask", dict(col=PR8, virus=True))], lines=["IAV PR8 (H1N1)", "MOI 0.5, 24 h", "n = 3"]),
            dict(icons=[], lines=["Biological replicates from", "separate cell passages", "and infection days"]),
        ]),
        dict(head="RNA and MeRIP", rows=[
            dict(icons=[("tube", {}), ("gel", {})], sep=2.6,
                 lines=["Total RNA; QC by NanoDrop", "and denaturing agarose gel"]),
            dict(icons=[("drop", dict(col=SPK))], lines=["m6A spike-in controls", "(positive and negative)"]),
            dict(icons=[("ip", {})], lines=["Anti-m6A immunoprecipitation"]),
            dict(split=True, labels=[("IP", "m6A-enriched"), ("Sup", "non-immunoprecipitated")]),
            dict(icons=[("drop", dict(col=CAL))], lines=["Calibration spike-in", "added to IP and Sup"]),
        ]),
        dict(head="Microarray", rows=[
            dict(icons=[("dyes", {})], lines=["Labelling: IP Cy5, Sup Cy3", "Co-hybridised, 65 °C, 17 h"]),
            dict(icons=[("slide", {})], lines=["Arraystar Human mRNA&lncRNA", "Epitranscriptomic Microarray",
                                               "(8×60K; Agilent-085753)"]),
            dict(icons=[("scanner", {})], lines=["Agilent G2505C scanner", "Feature Extraction v11.0.1.1"]),
        ]),
        dict(head="Data records", rows=[
            dict(icons=[("db", {})], lines=["GEO GSE336465"]),
            dict(icons=[("files", {})], lines=["6 raw Feature", "Extraction files"]),
            dict(icons=[("matrix", {})], lines=["Normalised matrix:",
                                                "%s probes × 12 channels" % fmt_int(n_probes)]),
            dict(icons=[("tables", {})], lines=["Provider profile and", "differential tables"]),
            dict(icons=[], lines=["m6A quantity = IP", "Expression = Sup",
                                  "%modified = 100 × IP/(IP + Sup)", "(linear scale)"]),
        ]),
    ]


# ---------------------------------------------------------------- build and check one layout
def build(pt, n_probes):
    LINE = pt * 1.22 / 72 * 25.4                 # one text line, mm
    COLS = columns(n_probes)

    def row_h(r):
        if r.get("split"):
            return SZ["split"][1]
        ih = max([SZ[n][1] for n, _ in r["icons"]], default=0.0)
        n = len(r["lines"])
        return ih + (((TXT_GAP if ih else 0.0) + n * LINE) if n else 0.0)

    base = [sum(row_h(r) for r in col["rows"]) + GAP_MIN * (len(col["rows"]) - 1) for col in COLS]
    inner = max(base)
    H = BOX_TOP + PAD_T + inner + PAD_B + BOT
    box_bot = H - BOT
    c = Canvas(W, H)
    ax = c.ax
    R = c.fig.canvas.get_renderer()

    def text_w(s, weight="normal"):
        t = ax.text(0, 0, s, fontsize=pt, fontweight=weight)
        e = t.get_window_extent(renderer=R)
        t.remove()
        return e.width / c.fig.dpi * 25.4

    icons, texts, heads, arrows = [], [], [], []

    def draw_icon(name, x0, y0, panel, **kw):
        c.cur = []
        ret = (f_split if name == "split" else FN[name])(c, x0, y0, **kw)
        w, h = SZ[name]
        icons.append(dict(name=name, artists=c.cur, panel=panel, decl=(x0, y0, x0 + w, y0 + h)))
        return ret

    def put(x, y, s, panel, **kw):
        t = ax.text(x, y, s, fontsize=pt, color=TXT, zorder=8, **kw)
        texts.append((t, panel))
        return t

    for k, col in enumerate(COLS):
        x0p, cx = PX[k], PX[k] + PW[k] / 2
        ax.add_patch(FancyBboxPatch((x0p, BOX_TOP), PW[k], box_bot - BOX_TOP,
                                    boxstyle="round,pad=0,rounding_size=1.6", fc=TINT[k], ec="none", zorder=1))
        heads.append((ax.text(cx, HEAD_Y, col["head"], ha="center", va="center", fontsize=pt,
                              fontweight="bold", color=HEADC, zorder=8), k))
        rows = col["rows"]
        hs = [row_h(r) for r in rows]
        gap = float(np.clip((inner - sum(hs)) / (len(rows) - 1), GAP_MIN, GAP_MAX))
        y = BOX_TOP + PAD_T + (inner - (sum(hs) + gap * (len(rows) - 1))) / 2
        for r, hh in zip(rows, hs):
            if r.get("split"):
                lw_ = max(max(text_w(a, "bold"), text_w(b)) for a, b in r["labels"])
                gw = SZ["split"][0] + 1.0 + lw_
                anchors = draw_icon("split", cx - gw / 2, y, k)
                for (axx, ayy), (a, b) in zip(anchors, r["labels"]):
                    put(axx + 1.0, ayy - 0.2, a, k, ha="left", va="bottom", fontweight="bold")
                    put(axx + 1.0, ayy + 0.2, b, k, ha="left", va="top")
            else:
                ih = max([SZ[n][1] for n, _ in r["icons"]], default=0.0)
                sep = r.get("sep", 0.0)
                tw = sum(SZ[n][0] for n, _ in r["icons"]) + sep * max(0, len(r["icons"]) - 1)
                xx = cx - tw / 2
                for n, kw in r["icons"]:
                    w, h = SZ[n]
                    draw_icon(n, xx, y + (ih - h) / 2, k, **kw)
                    xx += w + sep
                yt = y + ih + (TXT_GAP if ih else 0.0)
                for i, s in enumerate(r["lines"]):
                    put(cx, yt + i * LINE, s, k, ha="center", va="top")
            y += hh + gap

    ymid = (BOX_TOP + box_bot) / 2
    for k in range(3):
        a = FancyArrowPatch((PX[k] + PW[k] + 0.9, ymid), (PX[k + 1] - 0.9, ymid), arrowstyle="-|>",
                            mutation_scale=8, lw=0.9, color="#8C8C8C", shrinkA=0, shrinkB=0, zorder=8)
        ax.add_patch(a)
        arrows.append(a)

    # ---- checks
    c.fig.canvas.draw()
    R = c.fig.canvas.get_renderer()
    inv = ax.transData.inverted()

    def bb(a):
        e = a.get_window_extent(renderer=R)
        (x1, y1), (x2, y2) = inv.transform([(e.x0, e.y0), (e.x1, e.y1)])
        return (min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))

    def union(bs):
        return (min(b[0] for b in bs), min(b[1] for b in bs), max(b[2] for b in bs), max(b[3] for b in bs))

    def ov(a, b, tol=0.05):
        return min(a[2], b[2]) - max(a[0], b[0]) > tol and min(a[3], b[3]) - max(a[1], b[1]) > tol

    def inside(a, b, m):
        return a[0] >= b[0] + m and a[1] >= b[1] + m and a[2] <= b[2] - m and a[3] <= b[3] - m

    PANEL = [(PX[k], BOX_TOP, PX[k] + PW[k], box_bot) for k in range(4)]
    tb = [(bb(t), p, t.get_text()) for t, p in texts]
    ib = [(union([bb(a) for a in d["artists"]]), d["panel"], d["name"], d["decl"]) for d in icons]
    hb = [(bb(t), p, t.get_text()) for t, p in heads]
    ab = [bb(a) for a in arrows]
    prob = []
    for i in range(len(tb)):
        for j in range(i + 1, len(tb)):
            if ov(tb[i][0], tb[j][0]):
                prob.append("text/text: %r | %r" % (tb[i][2], tb[j][2]))
    for b, p, s in tb:
        for cc, q, n, _ in ib:
            if ov(b, cc):
                prob.append("text/icon: %r on %s" % (s, n))
        if not inside(b, PANEL[p], 0.6):
            prob.append("text outside panel %d: %r (x %.1f-%.1f, panel %.1f-%.1f)"
                        % (p + 1, s, b[0], b[2], PANEL[p][0], PANEL[p][2]))
        for a in ab:
            if ov(b, a):
                prob.append("text/arrow: %r" % s)
    for i in range(len(ib)):
        for j in range(i + 1, len(ib)):
            if ov(ib[i][0], ib[j][0]):
                prob.append("icon/icon: %s | %s" % (ib[i][2], ib[j][2]))
    for cc, q, n, d in ib:
        if not inside(cc, PANEL[q], 0.6):
            prob.append("icon outside panel %d: %s" % (q + 1, n))
        if cc[0] < d[0] - 0.4 or cc[1] < d[1] - 0.4 or cc[2] > d[2] + 0.4 or cc[3] > d[3] + 0.4:
            prob.append("icon %s drawn outside its declared box: %s vs %s"
                        % (n, tuple(round(v, 2) for v in cc), tuple(round(v, 2) for v in d)))
        for a in ab:
            if ov(cc, a):
                prob.append("icon/arrow: %s" % n)
    for b, p, s in hb:
        if b[0] < PANEL[p][0] or b[2] > PANEL[p][2] or b[1] < 0 or b[3] > BOX_TOP:
            prob.append("header outside its column: %r" % s)
    info = dict(H=H, n_text=len(tb), n_icon=len(ib), n_arrow=len(ab),
                slack=[round(PW[k] - 1.2 - max(b[2] - b[0] for b, p, s in tb if p == k), 1) for k in range(4)],
                labels=[t.get_text() for t, _ in heads] + [s for _, _, s in tb])
    return c, prob, info


def main():
    ip = pd.read_csv(cfg.figdata_dir + os.sep + "fig2d_ip_matrix.csv", index_col=0)
    n_probes = ip.shape[0]
    expect("probes", n_probes, 40955)
    font = os.path.basename(findfont(FontProperties(family=["sans-serif"])))
    print("[Figure 1] font: %s%s" % (font, "" if font.lower().startswith("arial") else
                                     "  (Arial not found: text size is reduced to fit if needed)"))
    used = None
    for pt in PT_STEPS:
        c, prob, info = build(pt, n_probes)
        if not prob:
            used = pt
            break
        print("[Figure 1] %.1f pt: %d layout problem(s)" % (pt, len(prob)))
        for p in prob[:6]:
            print("     ", p)
        plt.close(c.fig)
    if used is None:
        fail("Figure 1 layout checks fail down to %.1f pt" % PT_STEPS[-1])
    print("[Figure 1] layout check: PASS (%d labels, %d icons, %d arrows; %.1f pt; %.0f x %.1f mm)"
          % (info["n_text"], info["n_icon"], info["n_arrow"], used, W, info["H"]))
    print("[Figure 1] spare width per panel (mm): %s" % info["slack"])
    os.makedirs(cfg.figures_dir, exist_ok=True)
    out = cfg.figures_dir + os.sep + "Figure1_design"
    for ext in ("pdf", "png"):
        if os.path.exists(out + "." + ext):
            os.remove(out + "." + ext)
        c.fig.savefig(out + "." + ext, dpi=300)
    print("[Figure1_design] saved (%.2f MB pdf, %.2f MB png)"
          % (os.path.getsize(out + ".pdf") / 1e6, os.path.getsize(out + ".png") / 1e6))
    plt.close(c.fig)
    print("NUMBERS SHOWN:", json.dumps({"Figure 1": {"probes": n_probes, "text_pt": used}}, ensure_ascii=False))
    print("LABELS DRAWN:", json.dumps(info["labels"], ensure_ascii=False))


if __name__ == "__main__":
    main()
