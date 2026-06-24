# -*- coding: utf-8 -*-
"""11a_fig3_cache.py — one-shot extract Fig 3 data from the two big xlsx -> small npz cache.
Real data only (read_only, no writes to source). Entry-guard checks against locked facts:
  allcomparison rows ~49,849; volcano higher_in_PR8=4350 / higher_in_MOCK=3735.
Columns: col0 Transcript_ID, col2 GeneSymbol, col3 MOCK(mean), col4 PR8(mean),
         col5 Foldchange=2^(MOCK-PR8), col6 Pvalue, col8 Regulation, col24 m6A_mRNA_region.
Volcano perspective PR8/MOCK: log2FC = -log2(Foldchange).
"""
import os, time, numpy as np, openpyxl

RP   = r"D:\IAV_m6A\Report\Report\Report\PaperData\Differentially Methylated Genes"
RESDIR = r"D:\IAV_m6A\SciData_descriptor\results"
ALL  = os.path.join(RP, "mRNA_Quantity(allcomparison).xlsx")
DIFF = os.path.join(RP, "mRNA_Quantity(diffmethylated).xlsx")
NPZ  = os.path.join(RESDIR, "_fig3_cache.npz")
META = os.path.join(RESDIR, "_fig3_cache_meta.txt")

def find_header(ws):
    for ri, row in enumerate(ws.iter_rows(min_row=1, max_row=40, values_only=True)):
        if row[0] is not None and str(row[0]).strip() == "Transcript_ID":
            return ri
    return None

t0 = time.time()
# ---- read allcomparison MOCK_vs_PR8 ----
wb = openpyxl.load_workbook(ALL, read_only=True, data_only=True)
ws = wb["MOCK_vs_PR8"]
hr = find_header(ws)
fc=[]; pv=[]; reg=[]; mock=[]; pr8=[]; syms=[]
for ri, row in enumerate(ws.iter_rows(values_only=True)):
    if hr is None or ri <= hr: continue
    a = row[0]
    if a is None or str(a).strip().startswith("#") or str(a).strip() == "": continue
    try:
        f = float(row[5]); p = float(row[6])
    except (TypeError, ValueError): continue
    fc.append(f); pv.append(p)
    reg.append(str(row[8]).strip() if row[8] else "")
    try: mock.append(float(row[3])); pr8.append(float(row[4]))
    except (TypeError, ValueError): mock.append(np.nan); pr8.append(np.nan)
    syms.append(str(row[2]).strip() if row[2] else str(a).strip())
wb.close()
fc=np.array(fc); pv=np.array(pv); reg=np.array(reg)
mock=np.array(mock); pr8=np.array(pr8); syms=np.array(syms, dtype=object)
n_all = len(fc)

# ---- direction (empirical) ----
up_mask = reg == "up"; dn_mask = reg == "down"
md = float(np.nanmean(mock[up_mask] - pr8[up_mask]))
up_is_mock = md > 0  # 'up' => higher m6A in MOCK

# ---- volcano in PR8/MOCK perspective ----
log2fc = -np.log2(np.clip(fc, 1e-9, None))   # >0 => higher in PR8
nlogp  = -np.log10(np.clip(pv, 1e-300, None))
sig = (np.abs(log2fc) >= 1.0) & (pv < 0.05)
hi_pr8  = sig & (log2fc > 0)
hi_mock = sig & (log2fc < 0)
n_pr8 = int(hi_pr8.sum()); n_mock = int(hi_mock.sum())

# ---- top20 transcripts: 10 higher-in-PR8 + 10 higher-in-MOCK by |log2fc| among sig ----
order_sig = [i for i in np.argsort(-np.abs(log2fc)) if sig[i]]
top_pr8  = [i for i in order_sig if log2fc[i] > 0][:10]
top_mock = [i for i in order_sig if log2fc[i] < 0][:10]
sel = top_mock[::-1] + top_pr8[::-1]
top_syms = np.array([syms[i] for i in sel], dtype=object)
top_log2fc = np.array([log2fc[i] for i in sel], dtype=float)

# ---- region distribution col24 from diffmethylated hyper+hypo ----
def region_counts(sheet):
    wb2 = openpyxl.load_workbook(DIFF, read_only=True, data_only=True); ws2 = wb2[sheet]
    hr2 = find_header(ws2); cnt = {"5'UTR":0, "CDS":0, "3'UTR":0}
    for ri, row in enumerate(ws2.iter_rows(values_only=True)):
        if hr2 is None or ri <= hr2: continue
        a = row[0]
        if a is None or str(a).strip().startswith("#") or str(a).strip() == "": continue
        regn = row[24]
        if regn:
            for part in str(regn).split(";"):
                pp = part.strip()
                if pp in cnt: cnt[pp] += 1
    wb2.close(); return cnt
rc_hyper = region_counts("hyper_MOCK_vs_PR8")
rc_hypo  = region_counts("hypo_MOCK_vs_PR8")
region_total = {k: rc_hyper[k] + rc_hypo[k] for k in rc_hyper}
region_keys = ["5'UTR", "CDS", "3'UTR"]
region_vals = np.array([region_total[k] for k in region_keys], dtype=int)

# ---- entry guard against locked facts ----
warns = []
if not (49849 - 50 <= n_all <= 49849 + 50):
    warns.append(f"WARNING: allcomparison rows={n_all} outside 49849+-50 -> wrong columns/header?")
if abs(n_pr8 - 4350) > 50:
    warns.append(f"WARNING: higher_in_PR8={n_pr8} deviates >50 from 4350 -> column/direction read error?")
if abs(n_mock - 3735) > 50:
    warns.append(f"WARNING: higher_in_MOCK={n_mock} deviates >50 from 3735 -> column/direction read error?")
for w in warns: print(w)
if not warns: print("ENTRY GUARD: all locked-fact checks passed.")

# ---- save cache ----
np.savez_compressed(
    NPZ,
    log2fc=log2fc.astype(np.float32), nlogp=nlogp.astype(np.float32),
    sig=sig, hi_pr8=hi_pr8, hi_mock=hi_mock,
    region_keys=np.array(region_keys, dtype=object), region_vals=region_vals,
    top_syms=top_syms, top_log2fc=top_log2fc,
    n_all=n_all, n_pr8=n_pr8, n_mock=n_mock,
)
with open(META, "w", encoding="utf-8") as f:
    f.write("Fig3 cache meta (human-readable verification)\n")
    f.write(f"allcomparison parsed rows : {n_all}  (locked 49,849)\n")
    f.write(f"Regulation up/down        : up={int(up_mask.sum())} down={int(dn_mask.sum())}  (locked up=18,749 down=31,100)\n")
    f.write(f"'up' mean(MOCK-PR8 log2)  : {md:.4f} => up = higher m6A in {'MOCK' if up_is_mock else 'PR8'}\n")
    f.write(f"Volcano |log2FC|>=1 & p<0.05: higher_in_PR8={n_pr8} (locked 4350)  higher_in_MOCK={n_mock} (locked 3735)\n")
    f.write(f"Region hyper sheet        : {rc_hyper}\n")
    f.write(f"Region hypo sheet         : {rc_hypo}\n")
    f.write(f"Region combined           : {region_total}  (locked 5'UTR 12,159 / CDS 136,007 / 3'UTR 64,744)\n")
    tot = int(region_vals.sum())
    f.write(f"Region pct                : " + " / ".join(f"{k} {100*v/tot:.1f}%" for k,v in zip(region_keys, region_vals)) + "\n")
    f.write(f"top20 transcripts (sym, log2FC PR8/MOCK):\n")
    for s, v in zip(top_syms, top_log2fc):
        f.write(f"    {s}\t{v:+.3f}\n")
    if warns:
        f.write("\n".join(warns) + "\n")
    else:
        f.write("ENTRY GUARD: all locked-fact checks passed.\n")

print(f"region combined: {region_total}")
print(f"Regulation up={int(up_mask.sum())} down={int(dn_mask.sum())}")
print(f"SAVED cache -> {NPZ}  ({time.time()-t0:.1f}s)")
print(f"META -> {META}")
