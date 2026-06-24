# -*- coding: utf-8 -*-
"""40_verify_refs.py — 批量精确核实候选文献五字段(直连 PubMed,不经转述)
数据真实性/文献真实性铁律:输出 PubMed 原始字段,禁记忆。
"""
import urllib.request, json, time, os

# 候选 PMID 池(按主题,来自前面 esearch 真实命中)
CANDIDATES = {
 "m6A_general": ["33611339","34023900","36736310","38629637","30262497","31474572"],
 "m6A_virus":   ["28250115","39007267","36765030","33902106","31439758","37548590","37781923","34003800","31320558"],
 "iav_model":   ["40367638","31852346","23639425","28798254","32610119","29765910","33064776"],
 "flu_burden":  ["29248255","32196426","31027663","30553848","33476926","32373326"],
 "method_tools":["38418823","39695135"],  # 已核实的两篇 SciData 竞品
}

def esummary(ids):
    url=("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id="
         + ",".join(ids) + "&retmode=json")
    r=urllib.request.urlopen(url,timeout=30)
    return json.load(r)["result"]

out=[]
all_ids=sorted({pid for v in CANDIDATES.values() for pid in v})
# batch in chunks of 20
for i in range(0,len(all_ids),20):
    chunk=all_ids[i:i+20]
    res=esummary(chunk)
    for pid in chunk:
        rec=res.get(pid,{})
        if not rec or rec.get("error"):
            out.append(f"{pid}\tNOT_FOUND")
            continue
        au=rec.get("sortfirstauthor","?")
        src=rec.get("source","?")
        yr=rec.get("pubdate","?").split(" ")[0]
        title=rec.get("title","?").rstrip(".")
        # DOI
        doi=""
        for aid in rec.get("articleids",[]):
            if aid.get("idtype")=="doi": doi=aid.get("value","")
        pubtype=";".join(rec.get("pubtype",[]))
        out.append(f"{pid}\t{au}\t{src}\t{yr}\t{doi}\t[{pubtype}]\t{title}")
    time.sleep(0.4)

# tag by theme
theme_of={}
for th,ids in CANDIDATES.items():
    for pid in ids: theme_of.setdefault(pid,[]).append(th)

OUT=r"D:\IAV_m6A\SciData_descriptor\results\_ref_verify.txt"
with open(OUT,"w",encoding="utf-8") as f:
    f.write("PMID\tfirstauthor\tjournal\tyear\tDOI\t[pubtype]\ttitle\n")
    f.write("="*100+"\n")
    for line in out:
        pid=line.split("\t")[0]
        f.write(f"<{','.join(theme_of.get(pid,['?']))}> "+line+"\n")
print("verified",len(all_ids),"PMIDs ->",OUT)
