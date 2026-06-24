# -*- coding: utf-8 -*-
"""04_probe_columns.py — 精确探查绘图所需的列结构与现成 PCA 数据"""
import os, openpyxl
RP = r"D:\IAV_m6A\Report\Report\Report\PaperData"
out = []

def header_and_sample(path, sheet, label):
    out.append("="*70); out.append(f"[{label}] {sheet}")
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet]
    hdr_row = None
    for ri, row in enumerate(ws.iter_rows(min_row=1, max_row=40, values_only=True)):
        a = row[0]
        if a is not None and str(a).strip() == "Transcript_ID":
            hdr_row = ri
            # print full header with col index
            for ci, v in enumerate(row):
                if v is not None and str(v).strip():
                    out.append(f"   col{ci}: {v}")
            break
    # one data row
    if hdr_row is not None:
        for ri, row in enumerate(ws.iter_rows(min_row=hdr_row+2, max_row=hdr_row+2, values_only=True)):
            out.append("   --- sample data row ---")
            for ci, v in enumerate(row):
                if v is not None and str(v).strip():
                    out.append(f"   col{ci}= {v}")
    wb.close()

header_and_sample(os.path.join(RP,"Differentially Methylated Genes","mRNA_Quantity(diffmethylated).xlsx"),"hyper_MOCK_vs_PR8","DM mRNA")
header_and_sample(os.path.join(RP,"Differentially Methylated Genes","mRNA_Quantity(allcomparison).xlsx"),"MOCK_vs_PR8","DM mRNA ALL (for volcano)")

# PCA files
out.append("="*70); out.append("[PCA files]")
PCA = r"D:\IAV_m6A\Report\Report\盛江涛_H2108172_Human_m6A-mRNA&lncRNA Epi_补充分析_20220913_mk_zrb\pca"
for sub in os.listdir(PCA):
    p = os.path.join(PCA, sub, "pca.txt")
    if os.path.exists(p):
        out.append(f"--- {sub}/pca.txt ---")
        with open(p, encoding="utf-8", errors="replace") as f:
            for i, line in enumerate(f):
                if i < 10:
                    out.append("   "+line.rstrip())
                else: break

open(r"D:\IAV_m6A\SciData_descriptor\results\plot_columns.txt","w",encoding="utf-8").write("\n".join(out))
print("done")
