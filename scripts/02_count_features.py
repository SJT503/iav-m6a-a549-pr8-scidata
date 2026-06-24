# -*- coding: utf-8 -*-
"""02_count_features.py — 精确统计 DEG / 差异m6A 的真实条目数
Arraystar xlsx 前段是 # Column 注释行，需定位真实 header(含 Transcript_ID) 后计非空数据行。
数据真实性铁律：只报实际计数。
"""
import os, openpyxl
RP = r"D:\IAV_m6A\Report\Report\Report\PaperData"

def count_data_rows(path, sheet):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet]
    header_row = None
    n = 0
    for ri, row in enumerate(ws.iter_rows(values_only=True)):
        a = row[0]
        if header_row is None:
            if a is not None and str(a).strip() == "Transcript_ID":
                header_row = ri
            continue
        # after header: count rows whose col A looks like a transcript id
        if a is None:
            continue
        s = str(a).strip()
        if s and not s.startswith("#") and s.lower() != "transcript_id":
            n += 1
    wb.close()
    return header_row, n

targets = [
    ("DEG up",   os.path.join(RP, "Differentially Expressed Genes", "mRNA_ExpressionLevel(diffexpressed).xlsx"), "up_MOCK_vs_PR8"),
    ("DEG down", os.path.join(RP, "Differentially Expressed Genes", "mRNA_ExpressionLevel(diffexpressed).xlsx"), "down_MOCK_vs_PR8"),
    ("DM mRNA hyper",  os.path.join(RP, "Differentially Methylated Genes", "mRNA_Quantity(diffmethylated).xlsx"), "hyper_MOCK_vs_PR8"),
    ("DM mRNA hypo",   os.path.join(RP, "Differentially Methylated Genes", "mRNA_Quantity(diffmethylated).xlsx"), "hypo_MOCK_vs_PR8"),
    ("DM lncRNA hyper",os.path.join(RP, "Differentially Methylated Genes", "lncRNA_Quantity(diffmethylated).xlsx"), "hyper_MOCK_vs_PR8"),
    ("DM lncRNA hypo", os.path.join(RP, "Differentially Methylated Genes", "lncRNA_Quantity(diffmethylated).xlsx"), "hypo_MOCK_vs_PR8"),
    ("DM others hyper",os.path.join(RP, "Differentially Methylated Genes", "others_Quantity(diffmethylated).xlsx"), "hyper_MOCK_vs_PR8"),
    ("DM others hypo", os.path.join(RP, "Differentially Methylated Genes", "others_Quantity(diffmethylated).xlsx"), "hypo_MOCK_vs_PR8"),
]
out = []
for label, path, sheet in targets:
    hr, n = count_data_rows(path, sheet)
    out.append(f"{label:20s} header_row={hr}  data_rows={n}")

txt = "REAL FEATURE COUNTS (|FC|>=2.0 & p<0.05, MOCK_vs_PR8)\n" + "\n".join(out)
OUT = r"D:\IAV_m6A\SciData_descriptor\results\feature_counts.txt"
with open(OUT, "w", encoding="utf-8") as f:
    f.write(txt)
