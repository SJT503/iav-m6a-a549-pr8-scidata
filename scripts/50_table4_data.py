# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import cfg

# -*- coding: utf-8 -*-
"""Regenerate the Table 4 data with one rule for every gene (probe with the largest increase in expression, i.e. Sup log2 FC), and
recompute every number quoted from it. Also count interferon genes among the top m6A-quantity increases (mRNA probes, as in Figure 3c)."""
import os, csv, shutil
import pandas as pd, numpy as np
F = cfg.figdata_dir + os.sep
R = cfg.results_dir + os.sep
RAW = cfg.raw_dir + os.sep
S = cfg.samples
GENES = ["IFNL1", "IFNB1", "CCL5", "RSAD2", "ISG15", "PML", "IFNL3", "OAS2", "OAS1", "OASL", "MX2", "IFIT3"]
f3 = pd.read_csv(F + "fig3_limma_ip_sup.csv")
rows = []
for g in GENES:
    sub = f3[(f3.gene == g) & (f3.rna_class == "mRNA")]
    r = sub.loc[sub.sup_lfc.idxmax()]
    rows.append(dict(gene=g, transcript_id=r.transcript, probeID=r.probeID, sup_log2fc_pr8_vs_mock=r.sup_lfc, sup_fdr=r.sup_fdr,
                     ip_log2fc_pr8_vs_mock=r.ip_lfc, ip_fdr=r.ip_fdr, n_probes=len(sub)))
t = pd.DataFrame(rows).sort_values("sup_log2fc_pr8_vs_mock", ascending=False)
dst = R + "isg_transcript_table.csv"
# compare with a previous table only if one is already present (the first run on a clean
# checkout has none); keep a copy of the previous table when it exists.
if os.path.exists(dst):
    old = pd.read_csv(dst)
    changed = sorted(set(old.transcript_id) ^ set(t.transcript_id))
    print("transcripts changed vs previous table:", changed)
    prev = R + "isg_transcript_table_previous.csv"
    if not os.path.exists(prev):
        shutil.copyfile(dst, prev)
else:
    print("no previous isg_transcript_table.csv; writing a fresh one")
t.drop(columns=["probeID", "n_probes"]).to_csv(dst, index=False)
print(t[["gene", "transcript_id", "n_probes", "sup_log2fc_pr8_vs_mock", "sup_fdr", "ip_log2fc_pr8_vs_mock", "ip_fdr"]].to_string(index=False))

top4 = t[t.gene.isin(["IFNL1", "IFNB1", "CCL5", "RSAD2"])]; rest = t[~t.gene.isin(["IFNL1", "IFNB1", "CCL5", "RSAD2"])]
print("\ntop-4 fold: %.1f - %.1f (log2 %.2f - %.2f)" % (2 ** top4.sup_log2fc_pr8_vs_mock.min(), 2 ** top4.sup_log2fc_pr8_vs_mock.max(),
      top4.sup_log2fc_pr8_vs_mock.min(), top4.sup_log2fc_pr8_vs_mock.max()))
print("other 8 fold: %.2f - %.2f" % (2 ** rest.sup_log2fc_pr8_vs_mock.min(), 2 ** rest.sup_log2fc_pr8_vs_mock.max()))
print("max Sup FDR over 12 = %.6f (%s) | max IP FDR = %.6f" % (t.sup_fdr.max(), t.loc[t.sup_fdr.idxmax(), "gene"], t.ip_fdr.max()))

want = dict(zip(t.probeID, t.gene)); det = {p: [] for p in want}
for s in S:
    hdr = None
    for ln in open(RAW + s + ".txt", encoding="utf-8", errors="replace"):
        q = ln.rstrip("\n").split("\t")
        if q[0] == "FEATURES":
            hdr = {c: i for i, c in enumerate(q)}; continue
        if q[0] == "DATA" and hdr and q[hdr["ProbeName"]] in want:
            det[q[hdr["ProbeName"]]].append(int(q[hdr["gIsWellAboveBG"]]))
        elif q[0] == "TYPE":
            hdr = None
off_on = [want[p] for p, v in det.items() if v[:3] == [0, 0, 0] and v[3:] == [1, 1, 1]]
print("Sup channel undetected in all mock and detected in all PR8: %d of 12: %s" % (len(off_on), off_on))

ISG = {"IFNB1", "IFNL1", "IFNL2", "IFNL3", "RSAD2", "OASL", "CCL5", "PML", "ISG15", "IFIT1", "IFIT2", "IFIT3", "OAS1", "OAS2", "OAS3", "MX1", "MX2",
       "IFI44", "IFI44L", "GBP4", "CXCL10", "IFI6", "IFITM1", "STAT1", "IRF7", "DDX58", "IFIH1", "HERC5"}
g = f3[(f3.rna_class == "mRNA") & (f3.ip_call == "higher_PR8")].sort_values("ip_lfc", ascending=False)
for n in (10, 20):
    top = g.head(n)
    hits = [x for x in top.gene if x in ISG]
    print("mRNA top %d m6A-quantity increases: %d interferon-related entries %s" % (n, len(hits), hits))
