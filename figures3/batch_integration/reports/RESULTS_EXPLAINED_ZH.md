# 批次整合补充分析：结果解读与修改建议

## 一句话结论

Harmony 对批次混合产生了幅度很大的改善，而且重新建图聚类后原有细胞结构高度稳定；但残余批次结构没有完全消失。因此，论文应写成“批次混合大幅改善并保留生物学结构，同时仍存在细胞类型特异的残余效应”，不应写成“批次效应被完全消除”。

## 1. 批次混合到底改善了多少

正式 kBET（kBET 0.99.6，按细胞类型分层、每个 dataset-celltype 最多抽取 500 个细胞、30 个近邻、100 次重复）的总体结果如下：

| 表示空间 | 加权 kBET 拒绝率 | 怎样理解 |
|---|---:|---|
| PCA，整合前 | 0.940 | 约 94% 的局部邻域与应有的数据集构成不一致，批次结构很强 |
| 原 Harmony | 0.603 | 相比整合前明显下降，但仍存在残余结构 |
| 显式重跑 Harmony | 0.571 | 三者中最好；相比 PCA 绝对下降 0.370、相对下降约 39.3% |

kBET 的随机背景期望约为 0.044。拒绝率越低越好；0.571 仍高于随机背景，所以它证明的是“明显改善”，不是“完全消除”。kBET 对大样本和很小的局部偏差非常敏感，因此不能单独用它否定整个整合结果，需要同时看 iLISI、近邻混合、silhouette、细胞组成和生物学结构保存。

按细胞类型看，显式 Harmony 后：

- 混合相对较好：Neutrophils 0.338、Mast cells 0.380、Cycling cells 0.435。
- 中等残余：Endothelial cells 0.500、Epithelial cells 0.547、Fibroblasts 0.549、pDCs 0.568、Plasma cells 0.602、B cells 0.644。
- 残余最明显：Macrophages 0.753、Pericytes 0.834。
- Cycling cells、Mast cells 和 pDCs 至少有一个数据集的局部期望细胞数小于 5，结果受组成不均衡限制，应谨慎解释。

## 2. 其他批次指标怎么看

| 指标 | 整合前 | 显式 Harmony | 越大/越小越好 | 当前含义 |
|---|---:|---:|---|---|
| 30-NN 跨数据集邻居比例 | 0.268 | 0.509 | 越大越好 | 一个细胞周围来自其他数据集的邻居明显增多 |
| 跨数据集比例/组成期望值 | 0.424 | 0.797 | 越接近 1 越好 | 已达到按细胞类型组成所允许理想混合的约 80% |
| iLISI | 1.532 | 2.212 | 越大越好；3 个数据集时上限接近 3 | 局部邻域的数据集多样性明显提高 |
| dataset silhouette | 0.046 | -0.034 | 越接近 0 或略低越好 | 数据集标签不再形成清楚分隔的几何团块 |

这些指标与正式 kBET 的方向一致，因此“改善”不是 UMAP 视觉印象造成的。

## 3. 重新去批次以后，原聚类有没有被推翻

在 `harmony_explicit` 上重新运行 FindNeighbors 和 FindClusters（resolution = 0.3）后，仍得到 17 个 cluster，与原结果比较：

| 指标 | 结果 | 通俗解释 |
|---|---:|---|
| ARI | 0.952 | 1 表示完全相同；0.952 表明新旧细胞分群几乎一致 |
| NMI | 0.943 | 两套聚类包含的信息高度一致 |
| 新 cluster 映射回原 cluster 的细胞一致率 | 96.3% | 约 96 个/100 个细胞仍归入对应的原 cluster |
| 新 cluster 的细胞类型纯度 | 96.3% | 新聚类没有把主要细胞类型大范围混在一起 |
| 样本层面 cluster 构成相关性中位数 | 0.999 | 各样本的群体构成基本不变 |
| 样本层面总变差距离中位数 | 0.0208 | 新旧样本构成平均只相差约 2.1% |

四条预设稳定性规则全部通过。因此，显式重跑 Harmony 没有推翻原细胞类型结构和主要生物学结论。

## 4. 是否需要把所有单细胞分析重做

不建议把全部内容无条件重做。Harmony 修改的是降维坐标和邻居图，不修改原始 counts；UMAP 只负责展示，不是差异表达的输入。

推荐采用“原 Harmony 主分析 + 显式 Harmony 敏感性验证”的策略：

1. 主文或补充材料替换/增加整合前 PCA UMAP、整合后 Harmony UMAP、按 sample 着色 UMAP和按 cell type 着色 UMAP。
2. 补充正式 kBET、iLISI、30-NN 混合比例、dataset silhouette 和细胞类型分层结果。
3. 补充新旧聚类 ARI/NMI、重叠热图和样本构成稳定性，说明主要结构对 Harmony 重跑稳健。
4. 原 celltype 注释及依赖这些稳定注释的表达可视化、评分和多数下游结论可以保留；不需要因为 UMAP 更新而全部重跑。
5. 对 Pericytes 和 Macrophages 的疾病差异结论应做 dataset/sample 分层敏感性分析，避免残余批次驱动结论。

只有在论文决定把 `clusters_explicit` 作为全新的主聚类定义时，才需要重新检查 marker、细胞注释以及所有 cluster-dependent 分析。当前稳定性结果并不要求这样做。

## 5. NC 与 PD 差异分析仍需单独处理

批次整合良好并不能替代差异分析中的实验单位控制。细胞不是独立生物学重复；如果目前 NC vs PD 使用逐细胞 FindMarkers，仍建议改成 sample-level pseudobulk，并采用类似 `~ dataset + group` 的设计。若坚持细胞层面模型，应至少把 sample 作为随机效应、dataset 作为固定或随机效应，并说明模型假设。

这一修改与是否重新运行 Harmony 无关，是为了避免伪重复和 dataset/group 混杂。

## 6. 可用于回复审稿人的表述

### 中文

感谢审稿人指出数据集着色 UMAP 中可见的局部聚集。我们补充了整合前后对比及多项定量评估。按细胞类型分层的正式 kBET 分析显示，加权拒绝率由整合前 PCA 空间的 0.940 降至显式 Harmony 空间的 0.571；同时 iLISI 由 1.532 升至 2.212，跨数据集 30-NN 混合比例由 0.268 升至 0.509，dataset silhouette 由 0.046 降至 -0.034。这些结果表明 Harmony 大幅改善了批次混合，但仍存在细胞类型特异的残余数据集结构，主要见于周细胞和巨噬细胞。进一步在显式 Harmony 空间重新构建邻居图和聚类后，新旧聚类的 ARI 为 0.952，细胞映射一致率为 96.3%，样本层面 cluster 构成相关性中位数为 0.999，说明主要生物学结构对整合重跑稳健。我们已在修订稿中增加整合前后 UMAP、定量批次指标、细胞组成表和聚类稳定性分析，并在疾病差异分析中控制 dataset/sample 层面的影响。

### English

We thank the reviewer for noting the local dataset enrichment visible in the dataset-colored UMAP. We therefore added before-versus-after integration visualizations and quantitative, cell-type-stratified diagnostics. The weighted rejection rate from the official kBET implementation decreased from 0.940 in the pre-integration PCA space to 0.571 in the explicit Harmony space. In parallel, iLISI increased from 1.532 to 2.212, the cross-dataset 30-nearest-neighbor fraction increased from 0.268 to 0.509, and the dataset silhouette decreased from 0.046 to -0.034. These results demonstrate a substantial improvement in batch mixing, while also indicating cell-type-specific residual dataset structure, most notably in pericytes and macrophages. Reconstructing the neighbor graph and clusters in the explicit Harmony space yielded an adjusted Rand index of 0.952, 96.3% cell-level cluster concordance after majority mapping, and a median sample-level cluster-composition correlation of 0.999, supporting preservation of the major biological structure. We have added the integration diagnostics, dataset-composition tables, and clustering-stability analysis and have controlled dataset/sample-level effects in disease-associated comparisons.

## 7. 不能过度声称的内容

- 不写“批次效应被完全消除”。
- 不把 UMAP 上同色聚集全部解释为技术批次；其中相当一部分来自数据集间真实细胞类型组成差异。
- 不把 kBET 的改善称为传统意义上的统计学“显著差异”，除非另行定义并检验前后差值；更稳妥的用语是“substantially improved”。
- 不用细胞数量代替样本数量进行 NC vs PD 推断。

