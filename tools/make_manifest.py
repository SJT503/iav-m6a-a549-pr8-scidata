# -*- coding: utf-8 -*-
"""Regenerate MANIFEST.md5 for the release tree (files that are distributed, excluding
generated outputs and the GEO inputs, which have their own checksum manifests).

    python tools/make_manifest.py
"""
import hashlib, os

REL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKIP_DIRS = {".git", "results", "figdata", "figures", "tables", "interim", "__pycache__",
             "raw", "processed", "qc_evidence"}
SKIP_EXT = {".pyc", ".log", ".png", ".pdf", ".docx", ".xlsx", ".jpeg"}
OUT = os.path.join(REL, "MANIFEST.md5")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


rows = []
for root, dirs, files in os.walk(REL):
    dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS)
    for name in sorted(files):
        if name == "MANIFEST.md5" or os.path.splitext(name)[1] in SKIP_EXT:
            continue
        p = os.path.join(root, name)
        rel = os.path.relpath(p, REL).replace(os.sep, "/")
        rows.append((md5(p), rel))

with open(OUT, "w", encoding="utf-8", newline="\n") as f:
    for h, rel in rows:
        f.write("%s  %s\n" % (h, rel))
print("wrote %s with %d entries" % (OUT, len(rows)))
for h, rel in rows:
    print("  %s  %s" % (h, rel))
