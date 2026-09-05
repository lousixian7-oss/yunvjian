library(clusterProfiler)#GO富集分析、KEGG通路富集分析
library(org.Hs.eg.db)#基因注释数据库
library(enrichplot)
library(ggplot2)
library(GOplot)
library(dplyr)

#读入数
genes_df <- read.csv("核心基因集.csv")  

# 将基因符号转换为ENTREZ ID
entrezIDs <- bitr(genes_df$gene, fromType = "SYMBOL", 
               toType = c("ENTREZID", "SYMBOL"),
               OrgDb = org.Hs.eg.db) # `OrgDb`参数指定使用的人类基因组注释数据库
# gene_entrez
# # 创建手动补全的行
# mt_fix <- data.frame(
#   SYMBOL = c("MT-CO1", "MT-CO2"),
#   ENTREZID = c("4512", "4513")
# )

# # 合并到原始注释结果中（去重以防重复）
# entrezIDs <- rbind(entrezIDs, mt_fix) %>% distinct()

#使用entrezIDs 
gene<- entrezIDs$ENTREZID
##GO富集分析
go<- enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 1,ont="all",readable =T)
write.table(go,file="GO.txt",sep="\t",quote=F,row.names = F) #
library(data.table)
go <- fread("GO.txt", header = TRUE, sep = "\t", 
            stringsAsFactors = FALSE, 
            fill = TRUE,           # 自动填充缺失列
            blank.lines.skip = TRUE,
            quote = "")            # 忽略引号
# 假设 go 是 enrichResult 对象
library(ggplot2)
library(clusterProfiler)
selected_ids <- c(
  "GO:0005583",  # fibrillar collagen trimer
  "GO:0098643",  # fibrillar collagen
  "GO:0031012",  # extracellular matrix
  "GO:0005788",  # endoplasmic reticulum lumen
  "GO:0030020"   # extracellular matrix structural constituent conferring tensile strength
)
# 从 go 对象提取数据框
go_df <- as.data.frame(go)

# 筛选选中的 ID
go_sel <- go_df[go_df$ID %in% selected_ids, ]

# 计算 GeneRatio 数值（如果列为字符 "5/29"）
go_sel$GeneRatio_num <- sapply(go_sel$GeneRatio, function(x) eval(parse(text = x)))

# 按 ONTOLOGY 和 p.adjust 排序
go_sel <- go_sel %>%
  arrange(ONTOLOGY, p.adjust) %>%
  group_by(ONTOLOGY) %>%
  mutate(rank = row_number()) %>%
  ungroup()
library(scales)
# 绘图
p <- ggplot(go_sel, aes(x = GeneRatio_num, y = reorder(Description, -p.adjust))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "#F8766D", high = "#00BFC4", name = "p.adjust",
                       
                       labels = scientific_format(digits = 2)   # 使用 scientific_format
  ) +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
  scale_size_continuous(range = c(3, 8), name = "Count") +

  labs(x = "GeneRatio", y = NULL, title = "GO terms (AGT virtual overexpression)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, lineheight = 0.8),
        strip.text.y = element_text(angle = 0))

pdf("GO-气泡agtOE.pdf", width =10, height = 4)
print(p)
dev.off()


library(clusterProfiler)
library(dplyr)

# 读取两个 GMT 文件
gmt_legacy <- read.gmt("c2.cp.kegg_legacy.v2026.1.Hs.entrez.gmt")
gmt_medicus <- read.gmt("c2.cp.kegg_medicus.v2026.1.Hs.entrez.gmt")

# 合并并去重（使用正确的列名 `term` 和 `gene`）
gmt_combined <- bind_rows(gmt_legacy, gmt_medicus) %>%
  distinct(term, gene)   # 注意这里用 term，不是 gs_name

# 进行富集分析（gene 应为 Entrez ID 向量）
kk <- enricher(gene = gene,
               TERM2GENE = gmt_combined,
               pvalueCutoff = 0.1,        # 放宽 p 值
               qvalueCutoff = 1,        # 适度控制 FDR
               pAdjustMethod = "fdr",
               minGSSize = 1,
               maxGSSize = 1000)

write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)    
# 1. 提取 KEGG 结果数据框（如果还没有）
# 如果 kegg_res 是 enrichResult 对象
library(data.table)
kk <- fread("KEGG.txt", header = TRUE, sep = "\t", 
            stringsAsFactors = FALSE, 
            fill = TRUE,           # 自动填充缺失列
            blank.lines.skip = TRUE,
            quote = "")            # 忽略引号
# 假设 go 是 enrichResult 对象
kegg_df <- as.data.frame(kk)

id_mapping <- data.frame(
  ID = c(
    "KEGG_ECM_RECEPTOR_INTERACTION",
    "KEGG_FOCAL_ADHESION",
    "KEGG_MEDICUS_REFERENCE_WNT5A_ROR_SIGNALING_PATHWAY",
    "KEGG_MEDICUS_REFERENCE_WNT_SIGNALING_MODULATION_WNT_ACYLATION",
    "KEGG_RENIN_ANGIOTENSIN_SYSTEM"
  ),
  hsaID = c(
    "hsa04512",
    "hsa04510",
    "medicus",
    "medicus",
    "hsa04614"
  ),
  StandardName = c(
    "ECM-receptor interaction",
    "Focal adhesion",
    "WNT5A-ROR signaling pathway",
    "WNT signaling modulation",
    "Renin-angiotensin system"
  ),
  stringsAsFactors = FALSE
)
# 2. 筛选通路（只保留存在于结果中的）
selected_ids <- id_mapping$ID
kegg_sel <- kegg_df[kegg_df$ID %in% selected_ids, ]

# 3. 合并映射，生成新 Description（格式：hsaID | StandardName）
kegg_sel <- kegg_sel %>%
  left_join(id_mapping, by = "ID") %>%
  mutate(Description = paste0(hsaID, " | ", StandardName))

# 4. 计算 GeneRatio 数值
kegg_sel$GeneRatio_num <- sapply(kegg_sel$GeneRatio, function(x) eval(parse(text = x)))

# 5. 按 GeneRatio 降序排列（顶部最大）
kegg_sel <- kegg_sel %>%
  arrange(desc(GeneRatio_num)) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

# 5. 绘图（无分面）
p <- ggplot(kegg_sel, aes(x = GeneRatio_num, 
                          y = reorder(Description, -pvalue))) +
  geom_point(aes(size = Count, color = pvalue)) +
  scale_color_gradient(low = "#F8766D", high = "#00BFC4", name = "p.adjust") +
  scale_size_continuous(range = c(3, 8), name = "Count") +
  labs(x = "GeneRatio", y = NULL, title = "KEGG pathways (AGT virtual overexpression)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, lineheight = 0.9),
        plot.title = element_text(hjust = 0.5))

# 6. 保存
pdf("KEGG-气泡图agtoe.pdf", width = 8, height = 4)
print(p)
dev.off()

