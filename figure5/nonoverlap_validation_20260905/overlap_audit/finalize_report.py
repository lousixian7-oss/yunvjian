from pathlib import Path
import csv,json,hashlib
p=Path(__file__).resolve().parent
s=json.loads((p/'audit_summary.json').read_text())
x=list(csv.DictReader((p/'GSE10334_GSE16134_crosswalk.csv').open(encoding='utf-8-sig')))
d=list(csv.DictReader((p/'GSE16134_sample_disposition.csv').open(encoding='utf-8-sig')))
keep=list(csv.DictReader((p/'GSE16134_nonoverlapping_patient_samples.csv').open(encoding='utf-8-sig')))
assert len(x)==247 and len(d)==310 and len(keep)==66
assert len({r['GSE16134_sample'] for r in x if r['GSE16134_sample']})==244
assert all(r['development_patient_key']=='GSE10334_P'+r['patient'] for r in x if r['used_in_current_development']=='True')
assert not ({r['patient'] for r in keep}&{r['patient'] for r in x})
assert sum(r['group']=='NC' for r in keep)==6
inputs=[Path('G:/1Yunvjian/periodontitis基因/GEO/数据处理/GSE10334_series_matrix.txt.gz'),Path('G:/1Yunvjian/periodontitis基因/GEO/外部数据/GSE16134_analysis/GSE16134_series_matrix.txt.gz'),Path('G:/1Yunvjian/1A玉女煎/code/coredata/bulk/patient_split.csv')]
manifest=[]
for f in inputs:
 h=hashlib.sha256()
 with f.open('rb') as handle:
  for block in iter(lambda:handle.read(1024*1024),b''):h.update(block)
 manifest.append({'path':str(f),'bytes':f.stat().st_size,'sha256':h.hexdigest()})
(p/'input_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
report='''# GSE10334 与 GSE16134 患者／样本重叠核对

核对日期：2026-09-05。范围：原始 GEO 系列矩阵、当前 patient_split.csv、原论文，以及三个疑点样本的原始 CEL 文件。未修改模型、预测概率、论文、图、表或 TRIPOD 文件。

## 已核实结果

| 项目 | 核对结果 |
|---|---:|
| 当前开发集全部样本 | 172 |
| 其中 GSE10334 样本 | 154（训练 108、内部测试 46） |
| 这 154 个样本覆盖的 GSE10334 患者 | 90 |
| GSE16134 全部样本／患者 | 310／120 |
| 来自上述同一批 90 位患者的 GSE16134 样本 | 244（181 PD、63 NC） |
| 当前开发样本中有标题与表达谱双重支持的对应样本 | 152（训练 106、内部测试 46） |
| 按患者排除后的候选评估子集 | 66 个样本／30 位患者 |
| 候选子集构成 | 60 PD、6 NC；30 位患者贡献 PD 部位，其中 6 位也贡献 NC 部位 |

NC 在此表示牙周炎患者的临床健康取材部位，不表示 6 位没有牙周炎的健康受试者。

## 对应关系如何验证

1. 从两个 GEO 原始 series matrix 读取 GSM、患者编号、取材部位和病变／健康标签，与当前开发集的实际样本清单交叉核对。GSM 不同不能证明生物样本独立。
2. 对全部 54,675 个共同探针计算跨队列表达谱 Pearson 相关，并将每个 GSE10334 样本与全部 310 个 GSE16134 样本比较。
3. 244 个原队列样本的表达谱最近邻同时满足患者、取材部位、表型相同，相关系数范围 0.9978448914–0.9987776525。这构成强对应证据；两次处理的表达数值并非完全相同，不能称为逐字节相同的数据，也没有对这 244 对原始 CEL 逐一做哈希核验。
4. 当前 GSE10334 开发集患者键与原始 GEO 标题中的患者编号全部一致，覆盖全部 90 位源队列患者。排除范围因此必须涵盖这些患者的所有 GSE16134 样本，而不限于重复的 152 个样本或训练集重复样本。

## 三个未建立对应的样本

| GSE10334 样本 | 源患者／部位 | 当前开发集中是否使用 | 处理 |
|---|---|---|---|
| GSM261185 | 患者 42，部位 4，PD | 否 | 不把表达最近邻当成同一样本 |
| GSM261186 | 患者 43，部位 1，PD | 是，训练集 | 同上 |
| GSM261293 | 患者 82，部位 3，NC | 是，训练集 | 同上 |

这些样本的最近邻相关仅为 0.9839–0.9884，且部位、表型或患者编号不一致，不能建立可靠样本身份对应。下载这三个 GEO 原始 CEL.gz、解压后计算 SHA-256，并与本地 GSE16134_RAW.tar 中全部 310 个 CEL 比较，均未发现完全相同的原始文件。哈希不匹配本身不排除文件重写等可能；因此对应表将其标为 unresolved_no_identity_claim，而不是擅自改患者编号或分组。三位患者仍在共享患者排除范围内。

## 原论文与实际 GEO 清单的差异

Papapanou 等 2009 年论文 Discussion 明确说明，120 位患者／310 个芯片的研究包含此前报告的 90 位患者／247 个芯片。实际本地 GEO 元数据及表达谱支持 244 个对应样本；上述三个样本未能建立对应。论文的总数陈述不能替代实际样本审核，也不能直接将 310−247=63 作为可用评估样本数。具体差异原因没有证据确定，本核对不推测。

原文：[Papapanou et al., BMC Microbiology 2009](https://link.springer.com/article/10.1186/1471-2180-9-221)。
原始数据库：[GSE10334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10334)、[GSE16134](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16134)。

## 对论文结论的影响与下一步

- 原 GSE16134 全部 310 个样本的 AUC 等结果不能继续作为独立外部验证的证据；仅修改措辞不能消除重叠对性能估计的影响。
- 已输出 66 个候选样本清单，用于之后按预先固定的分析方案重新评估同一个锁定三基因模型。此次未计算这个子集的 AUC、校准或净获益，不能沿用全部 310 个样本的指标。
- 该子集按可用元数据与开发集无患者重叠，但仍来自同一来源研究／机构，应描述为来自同一来源研究的患者不重叠评估子集，不宜宣称为新机构、前瞻性或完全独立采集的验证队列。
- 只有 6 个 NC 部位，性能和校准不确定性可能较大。后续评估应以患者为聚类单位计算不确定性，公开所有结果，不按结果好坏选择样本。
- 后续还应明确预处理参数的估计范围：目前外部矩阵／标准化使用了完整队列信息。若重算子集，必须说明保留原处理还是在子集上重新估计，并保持锁定模型系数；不能简单筛选旧预测概率后称为完整独立验证。若要严格隔离原始处理，应在非重叠样本的原始 CEL 上重做相应预处理。这不是本次身份核对已完成的分析。

可立即采用、不会把未完成分析写成已完成的英文表述：

> GSE16134 overlaps with the GSE10334 development data at the participant level and was therefore not considered an independent external validation cohort. An audit of GEO metadata identified 244 GSE16134 samples from the 90 participants represented in the development data. The remaining 66 samples from 30 participants constitute a candidate participant-disjoint evaluation subset; performance in this subset has not yet been reassessed.

## 文件与复现

- `GSE10334_GSE16134_crosswalk.csv`：247 行源样本对应／疑点记录。未解决样本的对应 GSM 留空，最近邻只用于审核。
- `GSE16134_sample_disposition.csv`：310 个样本的患者重叠标志及排除建议。
- `GSE16134_nonoverlapping_patient_samples.csv`：66 个候选样本。
- `audit_summary.json`：计数。
- `raw_CEL_exception_check.json`：三个源 CEL 哈希和匹配结果。
- `input_manifest.json`：元数据与开发集划分输入的文件哈希。
- `audit_overlap.py`：元数据与表达谱核对脚本；matrix_cache.npz 为本次矩阵缓存，输入变化时需删除该缓存或重建后运行。
- `check_raw_exceptions.py`：三个例外原始文件下载及与完整外部原始归档的哈希比较。

结论边界：本报告核对的是队列身份与重叠，未对模型开发全流程进行新的泄漏审计，也未证明剩余子集满足其他全部独立验证条件。
'''
(p/'Overlap_audit_report_ZH.md').write_text(report,encoding='utf-8')
print('Verified tables and wrote audit report and input hashes.')
