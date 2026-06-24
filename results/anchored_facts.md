# IAV_m6A 真实数据锚点 (anchored_facts.md)

> **数据真实性铁律**：本文件所有数字均来自 `scripts/01_extract_anchors.py` / `02_count_features.py` / `03_enzymes.py` 对
> `D:\IAV_m6A\Report\` 商业报告文件的**实际读取**（2026-06-22 抽取）。
> descriptor 正文引用任何数字必须能追溯到本文件对应行。禁止凭记忆补全。
> 原始抽取产物见 `results/anchor_extraction_raw.txt` / `anchor_xlsx_dims.txt` / `feature_counts.txt` / `m6a_enzymes.txt`。

---

## 1. 项目与样本（来源：主报告 docx TABLE 0/1）

| 项 | 真实值 |
|---|---|
| 项目标题 | m6A-mRNA&lncRNA Epitranscriptomic Microarray Service for 6 Human samples |
| 项目编号 | **H2108172** |
| 委托日期 | Nov 18, 2021 |
| 委托人 / 单位 | 盛江涛 / 汕头大学医学院 |
| 服务商 | **Aksomics Inc.（Aksomics Services）** ⚠️ 旧 CLAUDE.md 推测的"康测"不准，以报告落款为准 |
| 物种 / 样本类型 | Human / Cell（A549 人肺泡上皮细胞，感染模型见下） |
| 样本总数 | **6**（MOCK ×3 + PR8 ×3，均有生物学重复 n=3/组） |
| 样本名 | PR8_1, PR8_2, PR8_3, MOCK_1, MOCK_2, MOCK_3 |
| RNA QC | 6 样本全部 pass / Status OK |

⚠️ 注：报告本身只写 "Human / Cell"，未在交付 docx 中明写 "A549" 与 "PR8 株感染" 细节。
湿实验参数来自用户（盛老师）口径（2026-06-22 确认）：

| 湿实验参数 | 真实值（用户确认） |
|---|---|
| 细胞系 | **A549 人肺泡上皮细胞，来源 ATCC** |
| 病毒株 | Influenza A virus **PR8（A/Puerto Rico/8/1934, H1N1）** |
| **MOI** | **0.5** |
| **感染时长** | **24 h**（收样时间点） |
| 对照 | MOCK（未感染） |

仍待核实（投稿前补）：A549 具体 ATCC 货号/批次、培养基与血清、病毒来源/扩增方式、感染时的吸附条件、
是否含 TPCK-trypsin、RNA 提取试剂。这些可在 Methods 细胞培养/病毒感染段补全。

---

## 2. 平台与实验流程（来源：主报告 docx P13/P15/P193-P217）

| 项 | 真实值 |
|---|---|
| 芯片平台 | **Arraystar Human mRNA&lncRNA Epitranscriptomic Microarray (8×60K, Arraystar)** |
| 原理 | MeRIP（m6A 抗体免疫沉淀）+ 双色芯片杂交（IP=Cy5 / Sup=Cy3） |
| m6A 抗体 | anti-N6-methyladenosine rabbit polyclonal（**Synaptic Systems, 202003**），2 µg/样本 |
| 磁珠 | Dynabeads™ M-280 Sheep Anti-Rabbit IgG（**Invitrogen, 11203D**） |
| RNase inhibitor | Enzymatics, Y9240L |
| 投入量 | 3–5 µg total RNA + m6A spike-in control |
| 标记试剂盒 | **Arraystar Super RNA Labeling Kit（AL-SE-005）** |
| 纯化 | RNeasy Mini Kit（QIAGEN, 74105） |
| 杂交 | 65 °C 17 h，Agilent Hybridization Oven |
| 扫描仪 | **Agilent Scanner G2505C**（双色通道） |
| 图像/原始数据软件 | **Agilent Feature Extraction software v11.0.1.1** |
| RNA 定量/质检仪 | NanoDrop ND-1000 + Bioanalyzer 2100 / Mops 电泳 |

---

## 3. 数据处理与归一化（来源：docx P15/P216-P217）

- 归一化：IP(Cy5) 与 Sup(Cy3) 原始强度按 **log2-scaled Spike-in RNA 强度均值**归一化。
- 探针保留标准：在 6 样本中**至少 3 个**具 Present(P) 或 Marginal(M) QC flag 的探针保留用于 "m6A quantity" 分析。
- **m6A quantity** = 基于 IP(Cy5) 归一化强度的甲基化量。
- **%Modified（% m6A）** = 基于 IP 与 Sup 归一化强度算出的甲基化百分比。
- 差异判定阈值：**|FC| ≥ 2.0 且 p < 0.05**（unpaired t-test）；另含 FDR 列。
- 聚类：R 软件层次聚类；GO：topGO（Fisher exact，p≤0.05）；Pathway：Fisher exact test。

---

## 4. 数据规模（来源：Matrix.xlsx 表头 + flag_pass dims）

| 项 | 真实值 | 来源 |
|---|---|---|
| 原始探针矩阵 | **40,956 探针 × 12 列**（6 IP + 6 Sup，normalized） | Raw Data Files\Matrix.xlsx 实读 R0 表头 + 40,956 行 |
| 原始 FE 文本 | 6 个 `{MOCK,PR8}_{1,2,3}.txt`，各 **62,986 行**，Agilent Feature Extraction 三段式（TYPE/FEPARAMS/DATA），含全部探针特征强度 | Raw Data Files\ 实读 |
| 扫描日期 / 仪器条码 | **2021-11-11**，Agilent G2505C（仪器号 US10450393），扫描协议 GE2_1100_Jul11 | MOCK_1.txt FEPARAMS 实读 |
| 基因组 build | **hg19 / GRCh37 (Feb 2009)** | .txt Grid_GenomicBuild |
| 芯片 design ID | **085753_D_F_20180926**（Grid 384×164） | .txt 实读 |
| 探针 ID 前缀 | ASHGV4… / ASHG19AP1B…（Arraystar 探针命名） | Matrix.xlsx R1-R3 |
| flag_pass mRNA | ~32,0xx 行（含 ~30 行 # 注释头，**不作精确条目数**） | mRNA_Quantity(flag_pass).xlsx dims 32074 |
| flag_pass lncRNA | ~6,9xx 行（同上含注释头） | lncRNA_Quantity(flag_pass).xlsx dims 6933 |
| flag_pass others | ~2,0xx 行（pre-/pri-miRNA, snoRNA, snRNA） | others_Quantity(flag_pass).xlsx dims 2056 |

⚠️ flag_pass 行数含 Arraystar 注释头（每文件 ~30 行 `# Column X`），精确探针数须减注释头后重数；
正文如需"通过 QC 的探针总数"应另写脚本精确计数，勿直接用 dims。**当前唯一干净可引用的总数 = 40,956 探针（原始矩阵）。**

---

## 5. 差异分析真实条目数（精确计数，来源：02_count_features.py，MOCK_vs_PR8，|FC|≥2 & p<0.05）

| 类别 | hyper / up | hypo / down |
|---|---|---|
| **差异表达 mRNA (DEG)** | up **3,275** | down **2,271** |
| **差异 m6A 甲基化 mRNA** | hyper **3,735** | hypo **4,350** |
| **差异 m6A 甲基化 lncRNA** | hyper **586** | hypo **795** |
| **差异 m6A 其他 ncRNA** | hyper **17** | hypo **275** |

（计数方法：定位含 "Transcript_ID" 的真实表头行后，统计 col A 为转录本 ID 的非空非注释行。）

---

## 6. m6A 调控酶变化（来源：03_enzymes.py，完整清单，FC = PR8/MOCK 方向待核对）

> ⚠️ 全量输出（筛选铁律 Step1）。同一基因多行=多探针/多转录本。FC(GE)=表达倍数，FC(m6A)=甲基化倍数。

**Writers（甲基转移酶）**：METTL3, METTL14, WTAP, RBM15, RBM15B, KIAA1429(VIRMA)
- 突出：**METTL14** FC(GE)≈1.99, FC(m6A)≈1.98（表达与甲基化同步升高）；WTAP 部分探针 FC 低至 0.24–0.43。

**Readers（识别蛋白）**：YTHDF1, YTHDF2, YTHDF3, HNRNPC, HNRNPA2B1, EIF3A, EIF3B, ELAVL1, SRSF2, DGCR8
- 突出：**HNRNPA2B1** FC(GE)≈2.2–2.34, FC(m6A)≈2.41–2.49（升高）；**YTHDF3** 部分探针 FC(m6A)≈2.06–2.23（升高）；
  **YTHDF2** FC(GE)≈0.55–0.73, FC(m6A)≈0.55–0.60（降低）；ELAVL1 FC(m6A)≈0.45（降低）。

**Erasers（去甲基化酶）**：FTO, ALKBH5
- 突出：**FTO** FC(GE)≈0.51, FC(m6A)≈0.41（显著降低）；ALKBH5 变化小（FC≈0.91–1.18）。

⚠️ 推荐候选（筛选铁律 Step2，每个均"需 PubMed 核实"再写进正文）：METTL14↑、FTO↓、YTHDF2↓、HNRNPA2B1↑。
文献支持一律走 PubMed 实检，未核实不写作者/年份。

---

## 7. 已交付的现成分析模块（A 档·可直接整理成图表）

- QC 三件套：RNA_QC.pdf / MeRIP_QC.pdf / Labeling_Report.pdf（docx P229-P245）
- 层次聚类热图 8 张（全 m6A mRNA / lncRNA + 差异 mRNA / lncRNA…，GeneralData\Heatmap Plot）
- 火山图 / 散点图（MOCK_vs_PR8 mRNA ExpressionLevel）
- 表达-甲基化四象限关联（GE_m6A_Association）
- GO（BP/CC/MF）+ KEGG Pathway 富集（hyper/hypo 分别）
- 6 张芯片扫描原图（Array_Images，MOCK×3 + PR8×3）
- 2022-09 补充分析（PCA 等，盛江涛_H2108172_…补充分析_20220913 目录）

---

## 8. 关键缺口（投稿前必补，禁止编造）

1. ~~**湿实验参数**~~ ✅ 已确认（2026-06-22，用户）：A549 来源 **ATCC**，IAV **PR8**，**MOI=0.5**，**感染 24 h**。
   仅余次要细节（ATCC 货号、培养基、病毒扩增方式、是否加 TPCK-trypsin）待投稿前补，不阻塞起草。
2. ~~**原始数据是否在手**~~ ✅ 已确认在文件夹：`PaperData\Raw Data Files\{MOCK,PR8}_{1,2,3}.txt`（6 个 Agilent FE 原始强度文件）
   + `Matrix.xlsx`（40,956 探针归一化矩阵）+ `芯片数据上传说明-CN.ppt`（GEO 上传指引）。
3. **GEO accession**：未投，Data Records 先占位（GSExxxxxx）。原始 .txt + 矩阵齐全，可按 ppt 指引上传。
4. **平台 GPL 号**：Arraystar Human mRNA&lncRNA Epitranscriptomic Microarray (8×60K) 的 GEO 平台 GPL 编号待查/确认（投稿前）。
5. **服务商共享/发表政策**：docx P254 "数据仅保留一年" → 但原始 .txt 已在本地留存，公开无障碍；Aksomics 署名/致谢口径确认。
6. **生物学重复优势**：n=3/组 → 与 TBI 那篇 n=1 不同，**可做统计检验 + PCA 分群**（补充分析目录已有 PCA），Technical Validation 呈现。
