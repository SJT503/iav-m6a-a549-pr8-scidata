# -*- coding: utf-8 -*-
"""41_verify_refs2.py — 第二批候选文献精确核实(直连 PubMed)"""
import urllib.request, json, time

NEW = {
 "merip_method": ["26121403","26253969","33669207","31548708","29131848","36793308"],
 "iav_m6a_mech": ["35320754","35637949","38194973","34706873","35273176","38001065","34995800"],
 "arraystar":    ["35741720","36175708"],
 "topgo":        ["16199517","19910308"],  # topGO/Alexa 经典(手填,下面会核实)
 "m6a_immune":   ["37734910","37903470","36632449","40177773"],
}
def esummary(ids):
    url=("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id="
         +",".join(ids)+"&retmode=json")
    return json.load(urllib.request.urlopen(url,timeout=30))["result"]

theme_of={}
for th,ids in NEW.items():
    for pid in ids: theme_of.setdefault(pid,[]).append(th)
all_ids=sorted(theme_of)
out=[]
for i in range(0,len(all_ids),20):
    res=esummary(all_ids[i:i+20])
    for pid in all_ids[i:i+20]:
        rec=res.get(pid,{})
        if not rec or rec.get("error"):
            out.append((pid,"NOT_FOUND","","","","","")); continue
        doi=""
        for aid in rec.get("articleids",[]):
            if aid.get("idtype")=="doi": doi=aid.get("value","")
        out.append((pid,rec.get("sortfirstauthor","?"),rec.get("source","?"),
                    rec.get("pubdate","?").split(" ")[0],doi,
                    ";".join(rec.get("pubtype",[])),rec.get("title","?").rstrip(".")))
    time.sleep(0.4)
with open(r"D:\IAV_m6A\SciData_descriptor\results\_ref_verify2.txt","w",encoding="utf-8") as f:
    for pid,au,src,yr,doi,pt,ti in out:
        f.write(f"<{','.join(theme_of.get(pid,['?']))}> {pid}\t{au}\t{src}\t{yr}\t{doi}\t[{pt}]\t{ti}\n")
print("verified",len(all_ids))
