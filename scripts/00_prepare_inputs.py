# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import cfg

# -*- coding: utf-8 -*-
"""Stage 00: cache the nine provider xlsx tables as CSV so that later stages iterate quickly.
Input: data/raw/Matrix.xlsx and the eight tables in data/processed/ (the GEO deposit).
Output: data/interim/*.csv"""
import os, time, json
import pandas as pd

# Matrix.xlsx is deposited in raw/ with the Feature Extraction files; the provider
# result tables are in processed/. cfg.deposit() looks in both.
CACHE = cfg.interim_dir
cfg.make_dirs()

# file names come from [input_files] in config/config.ini
FILES = {"matrix": cfg.input_files["matrix"]}
FILES.update({k: v for k, v in cfg.input_files.items() if k != "fe_files" and k != "matrix"})

meta = {}
for key, fn in FILES.items():
    t0 = time.time()
    src = cfg.deposit(fn)
    dst = os.path.join(CACHE, key + ".csv")
    if os.path.exists(dst):
        print("[skip]", key, "(cached)")
        continue
    try:
        # header 可能不在第 1 行（有 # 注释行）
        raw = pd.read_excel(src, header=None, nrows=12)
        hdr = 0
        for i in range(min(12, len(raw))):
            row = [str(x) for x in raw.iloc[i].tolist() if pd.notna(x)]
            if len(row) >= 3 and "Transcript" in " ".join(row) + " " + " ".join(row):
                hdr = i; break
            if len(row) >= 3 and row[0] == "probeID":
                hdr = i; break
        df = pd.read_excel(src, header=hdr)
        df.to_csv(dst, index=False, encoding="utf-8-sig")
        meta[key] = {"file": fn, "header_row": int(hdr), "shape": list(df.shape),
                     "cols": [str(c) for c in df.columns][:30], "sec": round(time.time()-t0, 1)}
        print("[ok] %-22s shape=%s (hdr=%d, %.1fs)" % (key, df.shape, hdr, time.time()-t0))
    except Exception as e:
        meta[key] = {"file": fn, "error": "%s: %s" % (type(e).__name__, str(e)[:200])}
        print("[ERR]", key, type(e).__name__, str(e)[:160])

with open(os.path.join(CACHE, "_manifest.json"), "w", encoding="utf-8") as f:
    json.dump(meta, f, indent=2, ensure_ascii=False)
print("\nDONE")
