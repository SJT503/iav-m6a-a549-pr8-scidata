# STORY_FRAMEWORK — IAV m6A 芯片 · Scientific Data 数据描述符

> Data Descriptor 没有"假说→证据链"，它的"故事"是对**数据集本身**的三条断言。
> 跨日会话恢复时第一句先汇报本框架进度，再问当前要做什么。

---

## 核心问题（descriptor 要回答的）

把 **A549 感染 IAV PR8 株 vs MOCK 的 Arraystar m6A-mRNA & lncRNA 表观转录组芯片数据（6 样本，n=3/组）**
做成一篇 **Scientific Data Data Descriptor**——发表的是数据集本身，不是机制结论。

## 主线断言（descriptor 的三条 claim）

| # | 断言 | 用什么证明 |
|---|------|-----------|
| C1 | **数据真实可靠** | RNA QC（NanoDrop/Bioanalyzer 全 pass）、MeRIP 富集效率（spike-in qPCR）、标记 QC、6 样本来源可追溯（项目号 H2108172） |
| C2 | **技术质量可评估** | Spike-in 归一化、探针 QC flag 保留标准（P/M ≥3/6）、组内重复性（n=3 相关/PCA 分群）、m6A motif/region 分布、%Modified 分布 |
| C3 | **可被他人复用** | 标准 GEO 芯片文件结构（原始 .txt + 探针矩阵 + 处理后 xlsx）、mRNA+lncRNA+ncRNA 双层信息、可复现 pipeline 代码 |

## 不可妥协核心（焊死·违反即推翻全稿）

1. **venue = Scientific Data**：保持描述性。虽然 n=3 可做 t-test（数据自带），但 descriptor 定位是
   "提供数据 + 评估质量 + 支持复用"，**不下流感-m6A 机制结论**，差异分析作为"数据可复用性证据"呈现，不作机制论证。
2. **数据真实性**：所有数字必须来自 `results/anchored_facts.md`（实读商业报告）。禁止凭记忆补全任何
   样本数/平台型号/探针数/差异条目数/酶 FC。
3. **数据隔离**：只用 `D:\IAV_m6A\Report\` 自测数据。公共数据集仅作 Technical Validation 外部参照。
4. **文献真实性**：引用任何文献前 PubMed 核实五字段，未核实不引用。
5. **缺口诚实**：A549/PR8 湿实验参数商业报告不含 → 必须向盛老师核实补全，不得编造 MOI/感染时长/细胞来源。

## 可解释空间（允许讨论/调整）

- 具体补哪些 Technical Validation 项、用哪个公共 m6A/IAV 数据集做外部一致性
- Background & Summary 写多长、Methods 细分到什么粒度
- 图表数量与主/补分配（服从"讲清数据"最低需要）
- 是否把 m6A 酶变化、GO/KEGG 富集纳入正文 Technical Validation（作为复用价值演示）还是降为补充

## 当前进度

- **Phase**：精准科研 Phase 1（课题设计 / 项目骨架）
- **Step**：数据锚点✓ → 治理三件套建立中 → 待模板风格地图 → 起草五节正文
- **GEO**：未投 → Data Records 占位，投稿前回填
- **上次返工原因**：无（新子项目起点 2026-06-22）

## 文章结构（Scientific Data 固定五节，无 Results/Discussion）

Background & Summary · Methods · Data Records · Technical Validation · Usage Notes
（详见 HANDOFF.md 第 4 节 + refs/template_style_dissection.md）
