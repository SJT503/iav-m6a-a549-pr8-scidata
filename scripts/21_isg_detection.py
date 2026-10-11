# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import cfg

# -*- coding: utf-8 -*-
"""Detection status of the 12 interferon transcripts in each channel, from the FE files
   (rIsWellAboveBG = IP/Cy5; gIsWellAboveBG = Sup/Cy3), MOCK_1-3 | PR8_1-3."""
import os, csv
RAW = cfg.raw_dir
CACHE = cfg.interim_dir
S = cfg.samples
T = {"ENST00000333625": "IFNL1", "ENST00000380232": "IFNB1", "ENST00000605509": "CCL5", "ENST00000382040": "RSAD2",
     "ENST00000379389": "ISG15", "ENST00000565898": "PML", "ENST00000613087": "IFNL3", "ENST00000342315": "OAS2",
     "ENST00000202917": "OAS1", "ENST00000620239": "OASL", "ENST00000330714": "MX2", "ENST00000371811": "IFIT3"}
probe = {r["Transcript_ID"]: r["probeID"] for r in csv.DictReader(open(os.path.join(CACHE, "probe_map_mrna.csv"), encoding="utf-8")) if r["Transcript_ID"] in T}
want = set(probe.values()); flags = {p: {"r": [], "g": []} for p in want}
for s in S:
    hdr = None
    for ln in open(os.path.join(RAW, s + ".txt"), encoding="utf-8", errors="replace"):
        q = ln.rstrip("\n").split("\t")
        if q[0] == "FEATURES":
            hdr = {c: i for i, c in enumerate(q)}; continue
        if q[0] == "DATA" and hdr and q[hdr["ProbeName"]] in want:
            p = q[hdr["ProbeName"]]
            flags[p]["r"].append(int(q[hdr["rIsWellAboveBG"]])); flags[p]["g"].append(int(q[hdr["gIsWellAboveBG"]]))
        elif q[0] == "TYPE":
            hdr = None
n_ip = n_sup = 0
for t, g in T.items():
    f = flags[probe[t]]
    off_on = lambda v: v[:3] == [0, 0, 0] and v[3:] == [1, 1, 1]
    n_ip += off_on(f["r"]); n_sup += off_on(f["g"])
    print("%-6s IP %s|%s %-3s  Sup %s|%s %s" % (g, f["r"][:3], f["r"][3:], "*" if off_on(f["r"]) else "", f["g"][:3], f["g"][3:], "*" if off_on(f["g"]) else ""))
print("off in all mock & on in all PR8:  IP channel %d/12 | Sup channel %d/12" % (n_ip, n_sup))
