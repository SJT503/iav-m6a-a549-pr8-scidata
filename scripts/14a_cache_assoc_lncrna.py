# -*- coding: utf-8 -*-
"""14a_cache_assoc_lncrna.py — 缓存四象限关联 + lncRNA 差异数据(慢提取一次)
数据真实性铁律：实读，禁编造。方向：Foldchange=MOCK/PR8，画图转 PR8/MOCK。
"""
import os, numpy as np, openpyxl
RESDIR = r"D:\IAV_m6A\SciData_descriptor\results"
ASSOC = r"D:\IAV_m6A\Report\Report\Report\GE_m6A_Association\mRNA_(ExpressionLevel_Quantity)_MOCK_vs_PR8.xlsx"
LNC   = r"D:\IAV_m6A\Report\Report\Report\PaperData\Differentially Methylated Genes\lncRNA_Quantity(allcomparison).xlsx"

def find_header(ws, key="Transcript_ID"):
    for ri,row in enumerate(ws.iter_rows(min_row=1,max_row=40,values_only=True)):
        if row[0] is not None and str(row[0]).strip()==key:
            return ri
    return None

# ---- 1) association: FC(GE) col1, FC(m6A) col2 ----
wb=openpyxl.load_workbook(ASSOC,read_only=True,data_only=True)
ws=wb["MOCK_vs_PR8(all)"]; hr=find_header(ws)
ge=[]; m6=[]; sym=[]
for ri,row in enumerate(ws.iter_rows(values_only=True)):
    if ri<=hr: continue
    a=row[0]
    if a is None or str(a).strip().startswith("#") or str(a).strip()=="": continue
    try: g=float(row[1]); m=float(row[2])
    except (TypeError,ValueError): continue
    ge.append(g); m6.append(m); sym.append(str(row[4]).strip() if row[4] else "")
wb.close()
ge=np.array(ge); m6=np.array(m6)
# to PR8/MOCK log2
ge_l2 = -np.log2(np.clip(ge,1e-9,None))
m6_l2 = -np.log2(np.clip(m6,1e-9,None))
print(f"assoc rows={len(ge)}")

# ---- 2) lncRNA differential: FC col5, p col6, regulation col8 ----
wb2=openpyxl.load_workbook(LNC,read_only=True,data_only=True)
ws2=wb2["MOCK_vs_PR8"]; hr2=find_header(ws2)
lfc=[]; lp=[]
for ri,row in enumerate(ws2.iter_rows(values_only=True)):
    if ri<=hr2: continue
    a=row[0]
    if a is None or str(a).strip().startswith("#") or str(a).strip()=="": continue
    try: f=float(row[5]); p=float(row[6])
    except (TypeError,ValueError): continue
    lfc.append(f); lp.append(p)
wb2.close()
lfc=np.array(lfc); lp=np.array(lp)
lnc_l2 = -np.log2(np.clip(lfc,1e-9,None))  # PR8/MOCK
print(f"lncRNA rows={len(lfc)}")

np.savez(os.path.join(RESDIR,"_assoc_lnc_cache.npz"),
         ge_l2=ge_l2, m6_l2=m6_l2, sym=np.array(sym,dtype=object),
         lnc_l2=lnc_l2, lnc_p=lp)
# meta for human check
sig_lnc=(np.abs(lnc_l2)>=1)&(lp<0.05)
with open(os.path.join(RESDIR,"_assoc_lnc_meta.txt"),"w",encoding="utf-8") as f:
    f.write(f"assoc rows={len(ge)}\nlncRNA rows={len(lfc)}\n")
    f.write(f"lncRNA sig(|log2FC|>=1&p<0.05)={int(sig_lnc.sum())}  up_in_PR8={int((sig_lnc&(lnc_l2>0)).sum())} up_in_MOCK={int((sig_lnc&(lnc_l2<0)).sum())}\n")
    # 4-quadrant counts (both |log2|>=1)
    q=(np.abs(ge_l2)>=1)&(np.abs(m6_l2)>=1)
    f.write(f"assoc both|log2|>=1: n={int(q.sum())}\n")
    f.write(f"  GE+ m6A+ (both up in PR8): {int((q&(ge_l2>0)&(m6_l2>0)).sum())}\n")
    f.write(f"  GE+ m6A- : {int((q&(ge_l2>0)&(m6_l2<0)).sum())}\n")
    f.write(f"  GE- m6A+ : {int((q&(ge_l2<0)&(m6_l2>0)).sum())}\n")
    f.write(f"  GE- m6A- (both down in PR8): {int((q&(ge_l2<0)&(m6_l2<0)).sum())}\n")
print("CACHED")
