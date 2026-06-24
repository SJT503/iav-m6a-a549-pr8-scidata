# -*- coding: utf-8 -*-
"""15a_cache_pctmod.py (pandas版) — 快速缓存 %Modified 数据
换策略：openpyxl 逐格读 49K×64 太慢(>4min挂住)，改 pandas 只读需要的列。
数据真实性铁律：实读，禁编造。方向 FC=MOCK/PR8 -> 画图转 PR8/MOCK。
"""
import os, numpy as np, pandas as pd
RESDIR=r"D:\IAV_m6A\SciData_descriptor\results"
ASSOC=r"D:\IAV_m6A\Report\Report\Report\GE_m6A_Association\mRNA_(ExpressionLevel_Quantity)_MOCK_vs_PR8.xlsx"
QUANT=r"D:\IAV_m6A\Report\Report\Report\PaperData\Differentially Methylated Genes\mRNA_Quantity(allcomparison).xlsx"

def load_with_header(path, sheet, key="Transcript_ID", maxscan=40):
    # find header row by scanning first column
    raw=pd.read_excel(path, sheet_name=sheet, header=None, nrows=maxscan, engine="openpyxl")
    hr=None
    for i in range(len(raw)):
        if str(raw.iloc[i,0]).strip()==key: hr=i; break
    if hr is None: raise RuntimeError(f"header not found in {sheet}")
    df=pd.read_excel(path, sheet_name=sheet, header=hr, engine="openpyxl")
    return df

print("reading association...")
a=load_with_header(ASSOC,"MOCK_vs_PR8(all)")
print("assoc cols:", list(a.columns)[:6], "rows", len(a))
# col1 = Foldchange(GE); col0 Transcript_ID; col4 GeneSymbol — use positional
a_tid=a.iloc[:,0].astype(str).str.strip()
a_gefc=pd.to_numeric(a.iloc[:,1], errors="coerce")
a_sym=a.iloc[:,4].astype(str).str.strip()
ge_map=dict(zip(a_tid,a_gefc)); sym_map=dict(zip(a_tid,a_sym))

print("reading quant...")
q=load_with_header(QUANT,"MOCK_vs_PR8")
print("quant cols n=", len(q.columns), "rows", len(q))
q_tid=q.iloc[:,0].astype(str).str.strip()
# col9/10 mean %m6A; col17-22 sample %m6A
mockpm=pd.to_numeric(q.iloc[:,9],errors="coerce")
pr8pm =pd.to_numeric(q.iloc[:,10],errors="coerce")
samp=np.column_stack([pd.to_numeric(q.iloc[:,c],errors="coerce") for c in (17,18,19,20,21,22)])

# merge on tid present in association
mask=q_tid.isin(ge_map.keys()) & mockpm.notna() & pr8pm.notna() & np.isfinite(samp).all(axis=1)
idx=np.where(mask.values)[0]
tids=q_tid.values[idx]
ge=np.array([ge_map[t] for t in tids],float)
syms=np.array([sym_map.get(t,"") for t in tids],dtype=object)
mockpm=mockpm.values[idx]; pr8pm=pr8pm.values[idx]; samp=samp[idx]
ge_l2=-np.log2(np.clip(ge,1e-9,None))
dpm=pr8pm-mockpm
print(f"merged records={len(tids)}")

np.savez(os.path.join(RESDIR,"_pctmod_cache.npz"),
         tids=tids, ge_l2=ge_l2, mockpm=mockpm, pr8pm=pr8pm, dpm=dpm, samp=samp, syms=syms)
with open(os.path.join(RESDIR,"_pctmod_meta.txt"),"w",encoding="utf-8") as f:
    f.write(f"merged records={len(tids)}\n")
    f.write(f"delta %Modified (PR8-MOCK): mean={dpm.mean():.4f} sd={dpm.std():.4f} range=[{dpm.min():.3f},{dpm.max():.3f}]\n")
    f.write(f"%Modified MOCK mean={mockpm.mean():.4f}  PR8 mean={pr8pm.mean():.4f}\n")
    de=np.abs(ge_l2)>=1
    f.write(f"\nDECOUPLING CHECK (genes |expr log2FC|>=1, n={int(de.sum())}):\n")
    f.write(f"  |delta %Mod|<0.05 (expr changed, methylation density ~unchanged): {int((de&(np.abs(dpm)<0.05)).sum())}\n")
    f.write(f"  |delta %Mod|>=0.05: {int((de&(np.abs(dpm)>=0.05)).sum())}\n")
    se=np.abs(ge_l2)<0.5
    f.write(f"  expr ~stable(|log2FC|<0.5) but |delta %Mod|>=0.1 (true methylation change w/o expr change): {int((se&(np.abs(dpm)>=0.1)).sum())}\n")
print("CACHED pctmod (pandas)")
