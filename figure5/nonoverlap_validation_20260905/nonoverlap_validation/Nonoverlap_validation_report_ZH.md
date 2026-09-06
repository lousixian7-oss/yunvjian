# GSE16134 无患者重叠子集：锁定模型探索性评估

2026-09-05。主分析在 66 个保留样本的原始 CEL 文件上单独完成 RMA，再按既定外部队列内均值／样本标准差标准化规则应用锁定三基因 logistic 模型。未重新选基因、拟合预测模型系数、调阈值或应用校准修正。

## 主分析结果

| 指标 | 结果 |
|---|---|
| 样本与患者 | 66 个部位，30 位患者；60 PD、6 NC |
| AUC | 0.986（95% CI 0.950–1.000） |
| Brier score | 0.205（95% CI 0.131–0.286） |
| 校准截距（理想值 0） | 4.011（95% CI 3.159–5.001） |
| 校准斜率（理想值 1） | 1.727；bootstrap 不稳定，见下文 |
| 灵敏度（固定阈值 0.5） | 0.717（95% CI 0.600–0.833） |
| 特异度（固定阈值 0.5） | 1.000（6/6；精确二项 95% CI 0.541–1.000） |
| 平均预测概率 | 0.578 |
| 实际病变部位比例 | 0.909 |
| 逐个排除患者后的 AUC 范围 | 0.983–1.000 |

6 个 NC 均判对会使普通非参数 bootstrap 特异度区间退化为 1–1，不能理解为确定无误。此处对 6 个不同患者各贡献的一个 NC 部位，另给精确二项区间；bootstrap 原始结果仍完整保留。

校准斜率 bootstrap 仅 1289/2000 次得到收敛估计（其余包含 3 次单类别重抽样及 708 次不收敛拟合）。已保存的收敛结果百分位范围为 1.250–4.430，但因大量不收敛，不把该范围作为可靠的常规置信区间。

## 预处理对照（预先保留全部结果）

| 路径 | AUC | Brier | 校准截距 | 校准斜率 |
|---|---|---|---|---|
| 主分析：66 芯片单独 RMA + 66 样本标准化 | 0.986 | 0.205 | 4.011 | 1.727 |
| 对照：原 310 芯片 RMA + 66 样本标准化 | 0.992 | 0.208 | 3.967 | 2.102 |
| 对照：直接筛出原 310 样本的旧预测 | 0.989 | 0.122 | 2.826 | 2.529 |

这些处理路径是敏感性对照；主分析事先确定为仅用 66 个非重叠样本的原始 CEL 预处理，不根据 AUC 高低选择。旧 310 样本整体 AUC 不能用作这里的独立验证结果。

## 解释与限制

这是发现重叠后进行的事后敏感性评估；66 个样本此前包含在完整 GSE16134 评估中，因此不能将其描述为从未查看过结果的全新留出队列。

- 此子集与开发集没有已识别的患者重叠，但仍来自同一来源研究／机构，应称为患者不重叠的评估子集。它不是来自另一机构的新外部队列。
- 6 个 NC 样本是 6 位牙周炎患者的健康取材部位，不能写成 6 位无牙周炎健康受试者；所有 30 位患者均贡献病变部位。模型在此评估的是部位分类。
- 区分能力、概率校准和临床用途需分别解释。病变部位占比约 91%，此为选定样本构成；Brier 和校准不能直接解释为目标临床人群的诊断概率准确性。
- 仅按本子集患病部位比例计算的常数预测 Brier 为 0.083，仅作为描述性参照，不是外部训练的竞争模型。
- 决策曲线保留为探索性结果；没有经验证的目标人群患病率及临床决策情境，不据此声称临床净获益。
- 在检查的 0.05–0.95（步长 0.01）阈值中，模型同时优于 treat-all 和 treat-none 的网格点数为 8/91；完整数值见 decision_curve.csv。
- 置信区间来自 2,000 次按患者有放回抽样，患者内所有取材部位一起抽取。预测和预处理参数固定，反映对当前预测的条件评估不确定性；不包含开发模型或预处理重新估计的不确定性。单一类别重抽样不能计算 AUC，被跳过并保留缺失值，各指标有效重复数见 primary_metrics_95CI.csv。
- 外部队列内标准化使用整批 66 个样本的无标签表达分布，因而是批次适配的评估，尚不能证明对单个新患者直接部署时的性能。
- 逐个排除患者是固定预测下的影响诊断，未重新训练或重新标准化。

## 英文结果段落（供后续修订）

After excluding participants represented in the development data, GSE16134 contributed 66 tissue samples from 30 participants, including 60 affected and six unaffected sites. The retained arrays were independently preprocessed using RMA, and the locked AGT/CXCR4/FOS model was applied using the previously specified cohort-wise standardization rule without refitting its coefficients. The AUC was 0.986 (95% participant-cluster bootstrap CI, 0.950–1.000), and the Brier score was 0.205. The calibration intercept and slope were 4.011 and 1.727, respectively. This participant-disjoint subset originated from the same source study; the small number of unaffected sites and cohort-wise preprocessing limit interpretation as independent clinical validation.

## 复现与核验

- `../prepare_nonoverlap.py`：选取 66 个 CEL 并记录哈希，提取历史对照所需三基因数据。
- `../nonoverlap_validation.R`：完整 RMA、探针注释／均值聚合、锁定模型预测、患者 bootstrap、影响诊断和绘图。
- `predictions_66.csv`：逐样本三种处理路径预测；`transport_parameters.csv`：训练及外部标准化参数。
- `three_gene_probe_mapping.csv`：复用原规则：有效 SYMBOL、每 PROBEID 保留第一条，再按基因均值聚合。
- `sessionInfo.txt`：R 与软件包版本；`CEL_manifest.json`：原始文件哈希。
- Python 独立复算主分析概率、两两比较 AUC 及 Brier，与 R 输出一致（概率最大差异 < 1e-10）。
- 本次未覆盖正式论文、Figure、Table S6 或 TRIPOD；其中有关原 310 个样本独立验证的表述仍需后续统一修改。
