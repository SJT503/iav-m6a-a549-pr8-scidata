# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import cfg

# -*- coding: utf-8 -*-
"""(a) Matrix.xlsx 的归一化因子 == 转录本表的 provider 因子?  (matrix_norm - log2 FE signal 是否逐阵列常数)
   (b) 缓存 lncRNA / others 的 probeID->类别 映射"""
import os, csv, math, statistics, openpyxl, time
RAW = cfg.raw_dir
CACHE = cfg.interim_dir
BASE = cfg.processed_dir
S = cfg.samples

def fe_signal(s):
    out, header = {}, None
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
                r, g = float(d["rProcessedSignal"]), float(d["gProcessedSignal"])
                if nm.upper().startswith("ASH") and r > 0 and g > 0:
                    out.setdefault(nm, []).append((math.log2(r), math.log2(g)))
            elif p[0] == "TYPE":
                header = None
    return {k: (statistics.mean(x[0] for x in v), statistics.mean(x[1] for x in v)) for k, v in out.items()}

with open(os.path.join(CACHE, "matrix.csv"), encoding="utf-8-sig", newline="") as f:
    rd = csv.reader(f); mh = next(rd)
    rows = [r for r in rd if r and r[0]]
print("(a) Matrix normalized - log2(FE ProcessedSignal), per array")
for k, s in enumerate(S):
    fe = fe_signal(s)
    ci = mh.index("%s(IP,Normalized)" % s); cs = mh.index("%s(Sup,Normalized)" % s)
    dI = [float(r[ci]) - fe[r[0]][0] for r in rows if r[0] in fe]
    dS = [float(r[cs]) - fe[r[0]][1] for r in rows if r[0] in fe]
    print("   %-7s n=%d  IP: mean=%+.5f sd=%.1e   Sup: mean=%+.5f sd=%.1e" % (
        s, len(dI), statistics.mean(dI), statistics.pstdev(dI), statistics.mean(dS), statistics.pstdev(dS)))

print("\n(b) caching lncRNA / others probe maps")
for key, fn in (("lncrna", "lncRNA_Quantity(flag_pass).xlsx"), ("others", "others_Quantity(flag_pass).xlsx")):
    t0 = time.time()
    wb = openpyxl.load_workbook(os.path.join(BASE, fn), read_only=True, data_only=True)
    ws = wb.worksheets[0]; hdr = None; n = 0
    with open(os.path.join(CACHE, "probe_map_%s.csv" % key), "w", newline="", encoding="utf-8") as g:
        w = csv.writer(g); w.writerow(["probeID", "Transcript_ID", "GeneSymbol", "transcript_type"])
        for r in ws.iter_rows(values_only=True):
            if hdr is None:
                if r and str(r[0] or "").strip() == "Transcript_ID":
                    hdr = {str(c): i for i, c in enumerate(r) if c is not None}
                continue
            if r and r[0]:
                w.writerow([r[hdr["probeID"]], r[0], r[hdr.get("GeneSymbol", 2)], r[hdr.get("transcript_type", 1)]]); n += 1
    wb.close()
    print("   %-7s rows=%d (%.1fs) | probeID col=%s" % (key, n, time.time() - t0, hdr.get("probeID") if hdr else None))
