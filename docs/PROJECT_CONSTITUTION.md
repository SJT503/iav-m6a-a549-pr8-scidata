# PROJECT_CONSTITUTION — IAV m6A 芯片 · Scientific Data 数据描述符

> 生信分析门禁铁律落地文件（芯片版）。任何读数据的脚本必须经入口检验。
> 建立 2026-06-22。

---

## 第一章 · 唯一数据源对照表（焊死）

本子项目**只读** `D:\IAV_m6A\Report\` 商业交付数据，以**只读**方式引用，不复制污染、不修改原始文件。

| 用途 | 定稿文件（唯一） | 规模 | 状态 |
|------|-----------------|------|------|
| 原始探针强度矩阵（被描述主数据） | `PaperData\Raw Data Files\Matrix.xlsx` | 40,956 探针 × 12（6 IP + 6 Sup） | ✅ 定稿 |
| 原始强度 .txt | `PaperData\Raw Data Files\*.txt` | 6 样本 | ⚠️ 需确认在手（docx P254：仅保留一年） |
| m6A quantity 谱（flag_pass） | `GeneralData\Methylation Profiling\{mRNA,lncRNA,others}_Quantity(flag_pass).xlsx` | mRNA/lncRNA/others | ✅ 定稿 |
| 表达谱（flag_pass） | `GeneralData\Expression Profiling\mRNA_ExpressionLevel(flag_pass).xlsx` | mRNA | ✅ 定稿 |
| 差异 m6A | `PaperData\Differentially Methylated Genes\*_Quantity(diffmethylated).xlsx` | hyper/hypo | ✅ 定稿 |
| 差异表达 | `PaperData\Differentially Expressed Genes\mRNA_ExpressionLevel(diffexpressed).xlsx` | up/down | ✅ 定稿 |
| 主报告 | `Report\Human m6A-mRNA&lncRNA Epitranscriptomic Microarray Report.docx` | — | ✅ 方法/平台权威来源 |
| 补充分析（PCA 等） | `盛江涛_H2108172_…补充分析_20220913_mk_zrb\` | — | ✅ 定稿 |

数据物理路径前缀（从本目录出发）：`..\Report\Report\Report\`（注意报告解压后是三层 Report 嵌套）。

## 第二章 · 期望参数焊死（来源：results/anchored_facts.md 实读）

| 参数 | 焊死值 | 来源 |
|------|--------|------|
| 样本总数 | 6（MOCK ×3 + PR8 ×3） | docx TABLE 1 |
| 组别 | MOCK / PR8，各 n=3 | docx TABLE 1 |
| 平台 | Arraystar Human mRNA&lncRNA Epitranscriptomic Microarray (8×60K) | docx P13 |
| 扫描仪 / 软件 | Agilent G2505C / Feature Extraction v11.0.1.1 | docx P13/P15 |
| 原始探针数 | 40,956 | Matrix.xlsx 实读 |
| 差异阈值 | \|FC\| ≥ 2.0 & p < 0.05（unpaired t-test） | docx P99/P101 |
| DEG | up 3,275 / down 2,271 | 02_count_features.py |
| 差异 m6A mRNA | hyper 3,735 / hypo 4,350 | 02_count_features.py |
| 差异 m6A lncRNA | hyper 586 / hypo 795 | 02_count_features.py |
| 差异 m6A 其他 ncRNA | hyper 17 / hypo 275 | 02_count_features.py |

下游脚本禁止"重算"出与上表不符的数字而不报错；不一致先停下查原因。

## 第三章 · 入口检验（每个读数据的脚本必含·芯片版）

读取后立即执行，任一失败 `stop()`：

1. **样本数** ✓ —— 6 列样本（或 12 列 IP+Sup），不符报错
2. **组别完整** ✓ —— MOCK / PR8 各 3，组名正确
3. **探针数** ✓ —— 原始矩阵 40,956（处理后子集另记）
4. **QC flag 门禁** ✓ —— 用 spike-in 归一化后强度 / P-M flag，验证 6 样本分布无异常离群
5. **数据来源** ✓ —— 路径在 `..\Report\` 下，非跨项目（禁止误读任何 TBI 项目文件）

⚠️ 芯片无"细胞类型纯度"概念（A549 单一细胞系），故不套用 snRNA-seq 的标志物比值门禁；
改为 **spike-in 富集 + 探针 QC flag** 作为质量门禁。

## 第四章 · venue 红线（descriptor 专属）

- descriptor 定位 = 提供数据 + 评估质量 + 支持复用，**不下"PR8 通过 m6A 调控 X 致病"的机制结论**。
- n=3 自带统计力 → 差异分析（DEG/差异m6A/酶变化/富集）可呈现，但**作为"数据可复用性与内部一致性证据"**，
  不作机制论证；措辞用描述性（"the dataset captures…", "differential signals are provided for reuse"），
  不写因果（"PR8 infection drives…via m6A"）。
- 保持描述性也为后续机制论文（若有）保住发表资格（Nature 政策：descriptor 不阻碍后续研究论文，
  只要后者超越描述性分析）。

## 第五章 · 宪法修订流程（5 步）

改数据源 / 阈值 / 启用其他文件必须：提案 → 评估（影响范围）→ 用户审批 → 执行 → 验证（入口检验复跑）。
修订记录追加到本文件末尾。

---

## 修订记录

- 2026-06-22 建立。复用 TBI SciData_descriptor 子项目宪法结构，按芯片 + n=3 特性适配。
