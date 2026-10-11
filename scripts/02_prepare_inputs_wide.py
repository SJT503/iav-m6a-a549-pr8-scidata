# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import cfg

# -*- coding: utf-8 -*-
"""Step B: 重新缓存 flagpass 文件（50 列，覆盖 probeID 与表达层），供 ISG / 折叠 / 公式核验使用"""
import openpyxl, os, time, csv

BASE = cfg.processed_dir
CACHE = cfg.interim_dir

JOBS = [
    ("mrna_q_full", "mRNA_Quantity(flag_pass).xlsx", 0, 50),
    ("mrna_expr_full", "mRNA_ExpressionLevel(flag_pass).xlsx", 0, 50),
]

for key, fn, sidx, ncol in JOBS:
    dst = os.path.join(CACHE, key + ".csv")
    if os.path.exists(dst):
        print("[skip]", key, flush=True)
        continue
    t0 = time.time()
    wb = openpyxl.load_workbook(os.path.join(BASE, fn), read_only=True, data_only=True)
    ws = wb.worksheets[sidx]
    n = 0
    with open(dst, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        for row in ws.iter_rows(max_col=ncol, values_only=True):
            w.writerow(["" if c is None else c for c in row])
            n += 1
            if n % 5000 == 0:
                print("   %s: %d rows (%.0fs)" % (key, n, time.time() - t0), flush=True)
    wb.close()
    print("[ok] %s %d rows (%.1fs)" % (key, n, time.time() - t0), flush=True)
print("DONE")
