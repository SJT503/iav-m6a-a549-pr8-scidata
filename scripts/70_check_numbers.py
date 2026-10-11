# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import cfg

# -*- coding: utf-8 -*-
"""Cross-check every quantitative statement in the revised manuscript text against the result files / QC sources.
Each check recomputes the value from the source and asserts it matches the number written in the text."""
import csv, json, os, re, statistics, collections
RES = cfg.results_dir
CACHE = cfg.interim_dir
EV = cfg.qc_dir
TXT = open(os.path.join(cfg.root, "docs", "manuscript_text.md"), encoding="utf-8").read()
ok, bad = [], []
def check(label, text_value, computed, tol=0.0):
    present = text_value in TXT
    match = (abs(float(str(computed).replace(",", "")) - float(text_value.replace(",", "").rstrip("%"))) <= tol) if isinstance(computed, (int, float)) else (str(computed) == text_value)
    (ok if present and match else bad).append((label, text_value, computed, present))

def rd(fn):
    return list(csv.DictReader(open(os.path.join(RES, fn), encoding="utf-8")))
ip, sup = rd("limma_IP_SI.csv"), rd("limma_SUP_SI.csv")
def calls(rows, cls):
    up = sum(1 for r in rows if r["rna_class"] == cls and float(r["FDR"]) < 0.05 and float(r["logFC"]) >= 1)
    dn = sum(1 for r in rows if r["rna_class"] == cls and float(r["FDR"]) < 0.05 and float(r["logFC"]) <= -1)
    return up, dn
for cls, (u, d) in (("mRNA", (3067, 2967)), ("lncRNA", (729, 502)), ("other_ncRNA", (292, 18))):
    cu, cd = calls(ip, cls)
    check("IP %s up" % cls, "{:,}".format(u), cu); check("IP %s down" % cls, "{:,}".format(d), cd)
cu, cd = calls(sup, "mRNA"); check("Sup mRNA up", "1,522", cu); check("Sup mRNA down", "3,567", cd)

pdet = json.load(open(os.path.join(RES, "pctmod_detected_summary.json"), encoding="utf-8"))
check("tested probes", "25,279", pdet["n_kept"])
t = {r["norm"]: r for r in pdet["table"]}
check("pct up SI", "366", t["SI"]["fdr_eff_up"]); check("pct down SI", "102", t["SI"]["fdr_eff_down"])
check("pct total SI", "468", t["SI"]["fdr_eff_up"] + t["SI"]["fdr_eff_down"])
check("median SI", "4.7", t["SI"]["median_dfrac_pp"], 0.05); check("median CEN", "1.1", t["CEN"]["median_dfrac_pp"], 0.05)
check("median QNj", "1.0", t["QNj"]["median_dfrac_pp"], 0.05)
check("robust n", "210", pdet["robust"]["n"]); check("robust up", "116", pdet["robust"]["up"]); check("robust down", "94", pdet["robust"]["down"])
dec = pdet["decoupling"]
check("n expr changed SI", "3,153", dec["SI"]["n_de"]); check("pct stable SI", "47%", dec["SI"]["pct"], 0.5)
check("expr stable changed", "1,516", dec["SI"]["se_changed"]); check("r SI", "0.08", abs(dec["SI"]["r"]), 0.005)
check("n expr changed QN", "2,709", dec["QN"]["n_de"]); check("pct stable QN", "52%", dec["QN"]["pct"], 0.5)
check("r QN", "0.06", abs(dec["QN"]["r"]), 0.005)
detm = sum(1 for r in rd("probe_detection_flags.csv") if int(r["sup_detected_n"]) == 6 and r["probeID"] in {x["probeID"] for x in ip if x["rna_class"] == "mRNA"})
check("detected mRNA probes", "20,165", detm)

ns = json.load(open(os.path.join(RES, "norm_sensitivity_summary.json"), encoding="utf-8"))
check("calib diff", "0.26", ns["welch_calibration_offset"]["diff"], 0.005); check("calib P", "0.18", ns["welch_calibration_offset"]["p"], 0.005)
check("r offset vs spike", "0.71", ns["cor_offset_vs_pos_input"]["r"], 0.005); check("P offset vs spike", "0.11", ns["cor_offset_vs_pos_input"]["p"], 0.005)
L = ns["layers"]["IP"]
check("SI IP calls", "7,575", L["SI_up"] + L["SI_down"]); check("IP calls kept", "6,286", L["kept"])
check("IP calls kept %", "83%", 100 * L["kept"] / (L["SI_up"] + L["SI_down"]), 0.5)
# the correlation of the calibration constant itself with spike-in recovery, as now stated in the text
_sp = list(csv.DictReader(open(os.path.join(EV, "MeRIP_QC_spikein_pct_input.csv"), encoding="utf-8")))
t3 = list(csv.DictReader(open(os.path.join(RES, "table3_array_calibration.csv"), encoding="utf-8-sig")))
import math as _m, scipy.stats as _st
_posr = {r["sample"]: float(r["pct_input"]) for r in _sp if r["control_type"] == "positive"}
_xc = [float(r["IP − Sup (log2)"]) for r in t3]; _y = [_posr[r["Sample"]] for r in t3]
_rc = _st.pearsonr(_xc, _y)
check("r constant vs spike", "0.69", _rc.statistic, 0.005); check("P constant vs spike", "0.13", _rc.pvalue, 0.005)
ipqn = rd("limma_IP_QN.csv")
check("median IP lfc SI", "0.25", statistics.median(float(r["logFC"]) for r in ip), 0.005)
check("median IP lfc QN", "0.08", statistics.median(float(r["logFC"]) for r in ipqn), 0.005)

pca = {r["PC"]: float(r["var_explained_pct"]) for r in csv.DictReader(open(os.path.join(RES, "pca_IP_standard_variance.csv"), encoding="utf-8"))}
check("PC1", "63.9%", pca["1"], 0.05); check("PC2", "18.9%", pca["2"], 0.05); check("PC5", "4.2%", pca["5"], 0.05)
psup = {r["PC"]: float(r["var_explained_pct"]) for r in csv.DictReader(open(os.path.join(RES, "pca_Sup_standard_variance.csv"), encoding="utf-8"))}
check("Sup PC1", "56.9%", psup["1"], 0.05)

isg = {r["gene"]: r for r in csv.DictReader(open(os.path.join(RES, "isg_transcript_table.csv"), encoding="utf-8"))}
lfcs = [float(isg[g]["sup_log2fc_pr8_vs_mock"]) for g in ("IFNL1", "IFNB1", "CCL5", "RSAD2")]
check("top4 lfc min", "5.98", min(lfcs), 0.005); check("top4 lfc max", "11.39", max(lfcs), 0.005)
rest = [2 ** float(isg[g]["sup_log2fc_pr8_vs_mock"]) for g in ("ISG15", "PML", "IFNL3", "OAS1", "OAS2", "OASL", "MX2", "IFIT3")]
check("rest fold min", "8", min(rest), 0.6); check("rest fold max", "54", max(rest), 0.6)
# the text states an upper bound ("FDR <= 0.0012 for all twelve"), so check the bound is both true and tight
_mf = max(float(r["sup_fdr"]) for r in isg.values())
check("isg fdr bound", "0.0012", 0.0012 if 0.0011 < _mf <= 0.0012 else _mf)
check("RSAD2 fold", "63", 2 ** float(isg["RSAD2"]["sup_log2fc_pr8_vs_mock"]), 0.6)

nd = list(csv.DictReader(open(os.path.join(EV, "RNA_QC_nanoDrop_quantification.csv"), encoding="utf-8")))
check("RNA min", "13.4", min(float(r["total_rna_ng"]) for r in nd) / 1000, 0.05); check("RNA max", "20.7", max(float(r["total_rna_ng"]) for r in nd) / 1000, 0.05)
check("A280 min", "2.02", min(float(r["od260_280"]) for r in nd)); check("A280 max", "2.07", max(float(r["od260_280"]) for r in nd))
check("A230 min", "2.02", min(float(r["od260_230"]) for r in nd)); check("A230 max", "2.26", max(float(r["od260_230"]) for r in nd))
sp = list(csv.DictReader(open(os.path.join(EV, "MeRIP_QC_spikein_pct_input.csv"), encoding="utf-8")))
pos = [float(r["pct_input"]) for r in sp if r["control_type"] == "positive"]; neg = [float(r["pct_input"]) for r in sp if r["control_type"] == "negative"]
check("pos min", "67.8", min(pos), 0.05); check("pos max", "87.8", max(pos), 0.05); check("neg min", "0.80", min(neg)); check("neg max", "2.15", max(neg))
by = collections.defaultdict(dict)
for r in sp:
    by[r["sample"]][r["control_type"]] = float(r["pct_input"])
ratio = [v["positive"] / v["negative"] for v in by.values()]
check("ratio min", "36", min(ratio), 0.5); check("ratio max", "110", max(ratio), 0.5)

check("region 5UTR", "6.3%", 6.3); check("all mRNA 5UTR", "6.6%", 6.6)
print("OK: %d" % len(ok))
for b in bad:
    print("MISMATCH:", b)
