# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import cfg

# -*- coding: utf-8 -*-
"""流式缓存大文件（只取前 30 列 + 关键注释），逐块写 CSV"""
import openpyxl, os, time, csv
BASE = cfg.processed_dir
CACHE = cfg.interim_dir
cfg.make_dirs()

JOBS = [
    ("mrna_q_flagpass", "mRNA_Quantity(flag_pass).xlsx", 0, 30),
    ("mrna_expr_flagpass", "mRNA_ExpressionLevel(flag_pass).xlsx", 0, 16),
]

for key, fn, sidx, ncol in JOBS:
    dst = os.path.join(CACHE, key + ".csv")
    if os.path.exists(dst):
        print("[skip]", key); continue
    t0 = time.time()
    p = os.path.join(BASE, fn)
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    ws = wb.worksheets[sidx]
    n = 0
    with open(dst, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        for row in ws.iter_rows(max_col=ncol, values_only=True):
            w.writerow(["" if c is None else c for c in row])
            n += 1
            if n % 5000 == 0:
                print("   %s: %d rows (%.0fs)" % (key, n, time.time()-t0), flush=True)
    wb.close()
    print("[ok] %-20s %d rows x %d cols (%.1fs)" % (key, n, ncol, time.time()-t0), flush=True)
print("DONE")
