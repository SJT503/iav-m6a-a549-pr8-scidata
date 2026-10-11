# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import cfg

# -*- coding: utf-8 -*-
"""(1) 最后一次尝试: 按 FE 质量标记取 spike-in 特征作参照, 能否精确复现 provider 因子
   (2) 生成 probeID -> GeneSymbol / transcript_type / m6A_mRNA_region 映射, 供 R 端用"""
import os, csv, math, statistics
RAW = cfg.raw_dir
CACHE = cfg.interim_dir
S = cfg.samples
PROV = {"r": [10.65155, 10.17241, 10.85019, 10.54372, 10.31455, 10.12483],
        "g": [9.70681, 9.50237, 9.77775, 9.74515, 9.63554, 9.68526]}

def parse(path):
    feats, header = [], None
    with open(path, encoding="utf-8", errors="replace", newline="") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if not p:
                continue
            if p[0] == "FEATURES":
                header = p[1:]; continue
            if p[0] == "DATA" and header:
                feats.append(dict(zip(header, p[1:])))
            elif p[0] == "TYPE":
                header = None
    return feats

rules = {
    "all_spike": lambda d, c: True,
    "IsWellAboveBG": lambda d, c: d.get(c + "IsWellAboveBG") == "1",
    "IsPosAndSignif": lambda d, c: d.get(c + "IsPosAndSignif") == "1",
    "WellAbove_and_notSat": lambda d, c: d.get(c + "IsWellAboveBG") == "1" and d.get(c + "IsSaturated") == "0",
}
res = {k: {"r": [], "g": []} for k in rules}
for s in S:
    sp = [d for d in parse(os.path.join(RAW, s + ".txt")) if (d.get("ProbeName") or "").upper().startswith("ASSPINKEIN")]
    for k, fn in rules.items():
        for c in ("r", "g"):
            v = [math.log2(float(d[c + "ProcessedSignal"])) for d in sp if fn(d, c) and float(d[c + "ProcessedSignal"]) > 0]
            res[k][c].append(statistics.mean(v) if v else float("nan"))
print("(1) flag-rule reference vs provider (max abs error over 6 arrays)")
for k in rules:
    er = max(abs(a - b) for a, b in zip(res[k]["r"], PROV["r"]))
    eg = max(abs(a - b) for a, b in zip(res[k]["g"], PROV["g"]))
    print("   %-22s IP: %s | maxerr=%.4f   Sup: %s | maxerr=%.4f" % (
        k, " ".join("%.3f" % x for x in res[k]["r"]), er, " ".join("%.3f" % x for x in res[k]["g"]), eg))

# (2) probe map
hdr = None; n = 0
with open(os.path.join(CACHE, "mrna_q_full.csv"), encoding="utf-8-sig", newline="") as f, \
     open(os.path.join(CACHE, "probe_map_mrna.csv"), "w", newline="", encoding="utf-8") as g:
    w = csv.writer(g); w.writerow(["probeID", "Transcript_ID", "GeneSymbol", "transcript_type", "m6A_mRNA_region"])
    for r in csv.reader(f):
        if hdr is None:
            if r and (r[0] or "").strip() == "Transcript_ID":
                hdr = {c: i for i, c in enumerate(r)}
            continue
        if r and r[0]:
            w.writerow([r[hdr["probeID"]], r[0], r[hdr["GeneSymbol"]], r[hdr["transcript_type"]], r[hdr["m6A_mRNA_region"]]])
            n += 1
print("\n(2) probe_map_mrna.csv rows:", n)
