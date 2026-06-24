# HANDOFF —— IAV m6A 表观转录组芯片 · Scientific Data 数据描述符

> **给新会话的话**：这份文件是本子项目的完整起点。新对话开窗后第一步：读完本 HANDOFF + `results/anchored_facts.md`，
> 输出"故事身份卡"，再开始干活。本项目与 TBI 系列项目内容严格隔离，数据只用 `D:\IAV_m6A\Report\`。

---

## 0 · 一句话定义

把 **A549 人肺泡上皮细胞感染流感病毒 PR8 株 vs MOCK 的 Arraystar m6A-mRNA & lncRNA 表观转录组芯片数据
（MeRIP + 双色芯片，MOCK n=3 / PR8 n=3）** 做成一篇 **Scientific Data 数据描述符（Data Descriptor）**
——发表的是数据集本身，不是机制结论。

数据物理位置：`D:\IAV_m6A\Report\`（商业服务商 Aksomics Inc. 交付，项目号 H2108172，2021-11 委托，2022 交付）。
本子项目产出全部放 `D:\IAV_m6A\SciData_descriptor\`。

---

## 1 · 为什么投 Scientific Data（venue 论证）

- Data Descriptor 发表的是"数据集"本身：不需要假设检验、不需要机制结论、不需要 novelty，
  只需证明①数据真实可靠 ②技术质量可评估 ③可被他人复用。
- 本数据是一套**完整、QC 通过、有生物学重复（n=3/组）**的人源细胞 m6A 表观转录组芯片，
  覆盖 mRNA + lncRNA + 其他 ncRNA 的甲基化与表达双层信息，适合作为流感病毒-宿主 m6A 互作的复用资源。
- **与已发表 TBI snRNA-seq descriptor 子项目同构**（同一研究者、同一 venue、同一五节结构），
  可直接复用其投稿包结构、风格地图方法、三锁对齐流程。

### 与 TBI descriptor 的关键差异（适配点，焊死）

| 维度 | TBI snRNA-seq descriptor | 本 IAV m6A 芯片 descriptor |
|---|---|---|
| 数据类型 | 单核测序（有 FASTQ） | **双色荧光芯片（无 FASTQ，有 .txt 原始强度 + 探针矩阵）** |
| 重复 | n=1 pool/组（无生物学重复，是硬伤） | **n=3/组（有生物学重复，是优势）** |
| 统计 | 禁止组间定量检验 | **可做组内相关/PCA/t-test 描述性呈现**（但仍按 descriptor 定位，不下机制结论） |
| Data Records | GEO + FASTQ + matrix + Seurat rds | GEO（芯片平台 GPL）+ 原始 .txt + 探针强度矩阵 + 处理后 xlsx |
| Technical Validation 五轴 | 测序/比对质量、QC、相关、doublet/ambient、注释 | **杂交/标记 QC、MeRIP 富集效率（spike-in qPCR）、Spike-in 归一化、组内重复性/PCA、m6A motif/region 分布** |
| 纯度门禁 | 细胞类型标志物比值 | N/A（细胞系，无细胞类型混杂）；改为**探针 QC flag + spike-in 富集**门禁 |

---

## 2 · 数据家底（精确，全部来自 results/anchored_facts.md 实读）

见 `results/anchored_facts.md` 全文。核心：
- 6 样本（MOCK_1-3 / PR8_1-3），全部 RNA QC pass。
- 平台 Arraystar Human mRNA&lncRNA Epitranscriptomic Microarray (8×60K)；扫描 Agilent G2505C；Feature Extraction v11.0.1.1。
- 原始探针矩阵 40,956 探针 × 12（6 IP + 6 Sup，normalized）。
- 差异（|FC|≥2 & p<0.05，MOCK_vs_PR8）：DEG up 3,275 / down 2,271；差异 m6A mRNA hyper 3,735 / hypo 4,350；
  lncRNA hyper 586 / hypo 795；其他 ncRNA hyper 17 / hypo 275。
- m6A 酶：Writer METTL14↑；Eraser FTO↓；Reader YTHDF2↓、HNRNPA2B1↑、YTHDF3↑（FC 见 anchored_facts §6）。

---

## 3 · 数据隔离铁律（焊死）

1. 本子项目**只用** `D:\IAV_m6A\Report\` 下的自测芯片数据。
2. 公共数据集（GEO 其他 m6A/IAV 数据）**只能**作 Technical Validation 的外部一致性参照，不混入被描述的数据集本身。
3. 产出放 `D:\IAV_m6A\SciData_descriptor\`，不污染 TBI 系列项目。
4. 任何读数据的脚本，读取后立即加 QC 门禁（探针数=40,956 / 样本=6 / 组别 MOCK,PR8 各 3），不符 `stop()`。

---

## 4 · 文章结构（Scientific Data 固定五节，无 Results/Discussion）

Background & Summary · Methods · Data Records · Technical Validation · Usage Notes
（风格地图见 `refs/template_style_dissection.md`，由模板检索子 agent 拆解真实 SciData 芯片/m6A descriptor 得出。）

---

## 5 · 投稿前硬性要求

1. **GEO 上传**：原始 .txt + 探针强度矩阵 + 样本 metadata，拿到 GSE accession（descriptor 无 accession 不能投）。
   ⚠️ docx P254："数据仅保留一年" → 先确认原始 .txt 是否仍在手 / 能否公开。
2. **代码仓库**：QC→差异→富集→可视化全流程脚本上 GitHub + Zenodo DOI。
3. **APC**：Scientific Data 纯 OA，投稿前核实当前 APC（约 USD 2000+，查 nature.com/sdata）。
4. **湿实验参数核实**：A549 来源、IAV PR8 株、MOI、感染时长、收样点——商业报告不含，**向盛老师要**。
5. **服务商政策**：确认 Aksomics 数据共享/发表政策。

---

## 6 · 当前进度

- **Phase**：精准科研 Phase 1（项目骨架）→ 数据锚点已完成。
- **Step**：治理三件套建立中 → 待模板风格地图 → 起草五节正文。
- **GEO**：未投（用户确认）→ Data Records 用 accession 占位，投稿前回填。
- **上次返工原因**：无（新子项目起点 2026-06-22）。

*建立 2026-06-22 · 复用 TBI SciData_descriptor 子项目标准*
