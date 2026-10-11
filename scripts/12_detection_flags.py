# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import cfg

# -*- coding: utf-8 -*-
"""%modified 需要分母(Sup)可测: 统计每个探针在几个阵列上 gIsWellAboveBG==1 (Sup 检出) / rIsWellAboveBG==1 (IP 检出),
   并检查稳健集与 ISG 在 MOCK 中是否接近背景"""
import os, csv, collections, json
RAW = cfg.raw_dir
RES = cfg.results_dir
S = cfg.samples
det = {}
for k, s in enumerate(S):
    header = None
    with open(os.path.join(RAW, s + ".txt"), encoding="utf-8", errors="replace", newline="") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if not p:
                continue
            if p[0] == "FEATURES":
                header = p[1:]; continue
            if p[0] == "DATA" and header:
                d = dict(zip(header, p[1:]))
                nm = d["ProbeName"].strip()
                if nm.upper().startswith("ASH"):
                    det.setdefault(nm, [[0] * 6, [0] * 6])
                    det[nm][0][k] = int(d["gIsWellAboveBG"] == "1")
                    det[nm][1][k] = int(d["rIsWellAboveBG"] == "1")
            elif p[0] == "TYPE":
                header = None
matrix = [r["probeID"] for r in csv.DictReader(open(os.path.join(RES, "limma_pct_SI.csv"), encoding="utf-8"))]
sup_all6 = {p for p in matrix if sum(det[p][0]) == 6}
sup_mock3 = {p for p in matrix if sum(det[p][0][:3]) == 3}
print("matrix probes: %d | Sup detected in all 6 arrays: %d (%.1f%%) | Sup detected in all 3 MOCK: %d"
      % (len(matrix), len(sup_all6), 100 * len(sup_all6) / len(matrix), len(sup_mock3)))

def calls(fn):
    out = {}
    for r in csv.DictReader(open(os.path.join(RES, fn), encoding="utf-8")):
        fdr, df = float(r["FDR"]), float(r["dfrac"])
        out[r["probeID"]] = (1 if df >= 0.05 else -1 if df <= -0.05 else 0) if fdr < 0.05 else 0
    return out
cSI, cQ = calls("limma_pct_SI.csv"), calls("limma_pct_QNj.csv")
rob = {p for p in matrix if cSI[p] != 0 and cSI[p] == cQ[p]}
rob_f = rob & sup_all6
print("robust set: %d -> with Sup detected in all 6 arrays: %d (up %d / down %d)"
      % (len(rob), len(rob_f), sum(cSI[p] == 1 for p in rob_f), sum(cSI[p] == -1 for p in rob_f)))
print("robust probes lost by the filter: up %d / down %d"
      % (sum(cSI[p] == 1 for p in rob - rob_f), sum(cSI[p] == -1 for p in rob - rob_f)))
for tag, c in (("SI", cSI), ("QNj", cQ)):
    allc = {p for p in matrix if c[p] != 0}
    fc = allc & sup_all6
    print("%-3s calls: %d (up %d/down %d) -> filtered %d (up %d/down %d)" % (tag, len(allc),
          sum(c[p] == 1 for p in allc), sum(c[p] == -1 for p in allc), len(fc),
          sum(c[p] == 1 for p in fc), sum(c[p] == -1 for p in fc)))

gene = {r["probeID"]: (r["GeneSymbol"], r["direction"]) for r in csv.DictReader(open(os.path.join(RES, "pctmod_robust_set.csv"), encoding="utf-8"))}
print("\nISG / top robust probes: Sup detected (MOCK1-3 | PR8 1-3), IP detected (MOCK | PR8)")
for p, (g, d) in gene.items():
    if g in ("OASL", "GBP4", "IFIT1", "CCL5", "IFNL3", "MX2", "OAS2", "RSAD2", "BLM", "PRPF3", "TAF2"):
        print("   %-7s %-11s Sup %s | %s   IP %s | %s" % (g, d, det[p][0][:3], det[p][0][3:], det[p][1][:3], det[p][1][3:]))
json.dump({"sup_detected_all6": len(sup_all6), "robust": len(rob), "robust_filtered": len(rob_f),
           "robust_filtered_up": sum(cSI[p] == 1 for p in rob_f), "robust_filtered_down": sum(cSI[p] == -1 for p in rob_f)},
          open(os.path.join(RES, "pctmod_detection_filter.json"), "w"), indent=2)
with open(os.path.join(RES, "probe_detection_flags.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["probeID", "sup_detected_n", "ip_detected_n", "sup_detected_mock", "sup_detected_pr8"])
    for p in matrix:
        w.writerow([p, sum(det[p][0]), sum(det[p][1]), sum(det[p][0][:3]), sum(det[p][0][3:])])
