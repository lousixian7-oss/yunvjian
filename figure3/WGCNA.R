library("WGCNA")

datExpr0 = read.csv(
  "去除批次效应后表达矩阵.csv",
  row.names=1,
  check.names=FALSE
)

datExpr0 <- t(datExpr0)

datExpr0 <- as.data.frame(datExpr0)

dim(datExpr0)
keepGenes <- apply(
  datExpr0,
  2,
  function(x){
    sum(x > 1) >= 0.5*nrow(datExpr0)
  }
)
datExpr0 <- datExpr0[,keepGenes]

nGenes_keep <- 5000

geneVar <- apply(datExpr0,2,var)

topGenes <- names(
  sort(geneVar,decreasing=TRUE)[1:nGenes_keep]
)

datExpr <- datExpr0[,topGenes]

#主要看缺失值
gsg = goodSamplesGenes(datExpr0, verbose = 3)

gsg$allOK 
# 返回TRUE则继续
if (!gsg$allOK){
  # 把含有缺失值的基因或样本打印出来
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # 去掉那些缺失值
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


##样本过滤

sampleTree = hclust(dist(datExpr0), method = "average")

pdf(file = "1、聚类树.pdf", width = 10, height = 8)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


#如果有异常值就需要去掉，根据聚类图自己设置cutHeight 参数的值

##单独的图
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) +
  #想用哪里切，就把“h = 110”和“cutHeight = 110”中换成你的cutoff
  abline(h = 95, col = "red")

clust = cutreeStatic(sampleTree, cutHeight = 95, minSize = 10)
keepSamples = (clust==1)    #1是保留的
datExpr = datExpr[keepSamples, ]
datTraits <- datTraits[rownames(datExpr),]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)
dim(datTraits)
#没有异常样本就不需要去除
# datExpr = datExpr0
#输入表型信息

#要求必须是数值型，要么像年龄那样的数字，要么搞成0，1，或者是1，2，3等。
##如果上一步删除了某个样本，这一步一定记得也要抱=把样本删除！！！

datTraits = read.csv("分类（新）.csv", row.names = 1)

sampleTree2 = hclust(dist(datExpr), method = "average")
# 各个样本的表现: 白色表示低，红色为高，灰色为缺失
traitColors = numbers2colors(datTraits, signed = FALSE)
# 把样本聚类和表型绘制在一起
pdf(file = "2、样本表现.pdf", width = 10, height = 8)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

##下面是正式的WGCNA分析

##软阈值的筛选
##设置一系列软阈值，范围是1-30之间
powers = c(1:10, seq(from = 12, to=30, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#这个结果就是推荐的软阈值
sft$powerEstimat


cex1 = 0.9 #一般是0.9.不能低于0.85
pdf(file = "3.软阈值的选择.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (Power)",
     ylab="Scale Free Topology Model Fit,Signed R^2",type="n",
     main = paste("Scale Independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (Power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#进一步构建网络

power = sft$powerEstimat

##需要时间

net = blockwiseModules(
  datExpr,
  
  TOMType="signed",
  
  minModuleSize=100,
  
  deepSplit=2,
  
  mergeCutHeight=0.35,
  
  numericLabels=TRUE,
  
  pamRespectsDendro=FALSE
)


#此处展示得到了多少模块，每个模块里面有多少基因。
table(net$colors)
##如果结果不合理，可以适当调整参数进行修改

#具体细节可参考：https://zhuanlan.zhihu.com/p/34697561



mergedColors = labels2colors(net$colors)
pdf(file = "4.聚类趋势.pdf", width = 6, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#每种颜色都代表着一种基因模块，即里面的基因功能是相似的
#灰色没意义

#保存每个模块对应的基因
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
gm = data.frame(net$colors)
gm$color = moduleColors
head(gm)

genes = split(rownames(gm),gm$color)
save(genes,file = "genes.Rdata")


#模块与表型的相关性

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
# 自动寻找与PD相关性最高模块

target_trait <- "PD"

target_module <- rownames(
  moduleTraitCor
)[which.max(
  abs(moduleTraitCor[,target_trait])
)]

target_module
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#热图
pdf(file = "5.相关性热图.pdf", width = 6, height = 10)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed (50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()


#相关系数最好大于0.8，实在没有，0.6或0.7也行，再小就不合适了。


#把gene module输出到文件
module_list <- split(
  colnames(datExpr),
  moduleColors
)

for(i in names(module_list)){
  
  expr_module <- datExpr[,module_list[[i]]]
  
  write.csv(
    expr_module,
    paste0(i,"_expression.csv")
  )
}

#GS与MM
#GS代表模块里的每个基因与形状的相关性
#MM代表每个基因和所在模块之间的相关性，表示是否与模块的趋势一致

filter_hub <- function(module_name,
                       gs_cutoff=0.2,
                       mm_cutoff=0.7){
  
  genes <- names(moduleColors)[moduleColors == module_name]
  
  GS <- geneTraitSignificance[genes,1]
  
  module_index <- match(module_name, modNames)
  
  MM <- geneModuleMembership[
    genes,
    module_index
  ]
  
  
  result <- data.frame(
    Gene=genes,
    GS=as.numeric(GS),
    MM=as.numeric(MM)
  )
  
  
  result <- result[
    abs(result$GS)>gs_cutoff &
      abs(result$MM)>mm_cutoff,
  ]
  
  
  result <- result[
    order(-abs(result$GS)),
  ]
  
  return(result)
}

turquoise_hub <- filter_hub(
  "turquoise",
  gs_cutoff=0.2,
  mm_cutoff=0.7
)


brown_hub <- filter_hub(
  "blue",
  gs_cutoff=0.2,
  mm_cutoff=0.7
)
# 输出

write.csv(
  turquoise_hub,
  "turquoise_hub.csv",
  row.names=FALSE
)


write.csv(
  blue_hub,
  "blue_hub.csv",
  row.names=FALSE
)


#第几列的表型是最关心的，下面的i就设置为几。
#与关心的表型相关性最高的模块赋值给下面的module。
traitData = datTraits 
i = 2 #替换自己想要的表型列
module ="blue"##替换自己的样本颜色#重新计算MM

MEs0 = moduleEigengenes(
  datExpr,
  moduleColors
)$eigengenes

MEs = orderMEs(MEs0)


geneModuleMembership = as.data.frame(
  cor(datExpr,
      MEs,
      use="p")
)


MMPvalue = as.data.frame(
  corPvalueStudent(
    as.matrix(geneModuleMembership),
    nSamples
  )
)


names(geneModuleMembership) =
  paste(
    "MM.",
    names(MEs),
    sep=""
  )
geneTraitSignificance =
  as.data.frame(
    cor(datExpr,
        instrait,
        use="p")
  )

module="blue"


column = match(
  paste0("MM.ME",module),
  names(geneModuleMembership)
)


moduleGenes = moduleColors==module

pdf(
  "6.模块相关性散点图.pdf",
  width=5,
  height=5
)
verboseScatterplot(
  abs(geneModuleMembership[moduleGenes,column]),
  abs(geneTraitSignificance[moduleGenes,1]),
  
  xlab=paste(
    "Module Membership in",
    module,
    "Module"
  ),
  
  ylab="Gene Significance",
  
  main=paste(
    "Module Membership vs Gene Significance\n",
    module,
    "module"
  ),
  
  cex.main=1.2,
  cex.lab=1.2,
  cex.axis=1.2,
  col=module
)
dev.off()

module="turquoise"


column = match(
  paste0("MM.ME",module),
  names(geneModuleMembership)
)


moduleGenes = moduleColors==module

pdf(
  "6.turquoise.pdf",
  width=5,
  height=5
)
verboseScatterplot(
  abs(geneModuleMembership[moduleGenes,column]),
  abs(geneTraitSignificance[moduleGenes,1]),
  
  xlab=paste(
    "Module Membership in",
    module,
    "Module"
  ),
  
  ylab="Gene Significance",
  
  main=paste(
    "Module Membership vs Gene Significance\n",
    module,
    "module"
  ),
  
  cex.main=1.2,
  cex.lab=1.2,
  cex.axis=1.2,
  col=module
)
dev.off()


#也可以通过这里找到GS和MM都大的基因
##大家有兴趣可以自己看一下，一般没见过文章里用这个数据的

f = data.frame(GS = abs(geneModuleMembership[moduleGenes, column]),
               MM = abs(geneTraitSignificance[moduleGenes, 1]))
rownames(f) = rownames(gm[moduleGenes,])
head(f)


##用基因相关性热图的方式展示加权网络

#每行每列代表一个基因。一般取400个基因画就够啦，不然电脑运行不了

nSelect = 400
set.seed(10)
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6)

select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]

# 再计算基因之间的距离树(对于基因的子集，需要重新聚类)
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
library(gplots)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
pdf(file = "7.加权网络.pdf", width = 4, height = 4)
plotDiss = selectTOM^7
diag(plotDiss) = NA #将对角线设成NA，在图形中显示为白色的点，更清晰显示趋势
TOMplot(plotDiss, selectTree, selectColors, col=myheatcol,main = "Network Heatmap Plot")
dev.off()


#模块与表型的相关性
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(cbind(MEs, instrait))

pdf(file = "8、模块和表型的相关性.pdf", width = 12, height = 8)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(4,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()

#也可以把上面的图分开来画

pdf(file = "9、模块聚类树.pdf", width = 10, height = 8)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()

pdf(file = "10、表型相关热图.pdf",  width = 10, height = 8)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(4,5,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)


#========================
# 读取127基因
#========================

gene127 <- read.csv(
  "pd&yunajian交集基因.csv",
  stringsAsFactors = FALSE
)
gene127 <- as.character(
  gene127$Intersection_Genes
)

wgcna_genes <- colnames(datExpr)

gene127_wgcna <- intersect(
  gene127,
  wgcna_genes
)

length(gene127_wgcna)

head(gene127_wgcna)

names(moduleColors) <- colnames(datExpr)

gene_module <- data.frame(
  Gene = gene127_wgcna,
  Module = moduleColors[gene127_wgcna]
)


table(gene_module$Module)

library(ggplot2)

p <- ggplot(
  gene_module,
  aes(x=Module)
)+
  geom_bar()+
  theme_classic()+
  labs(
    x="WGCNA Module",
    y="Number of overlapping genes"
  )

ggsave(
  "127gene_module_distribution.pdf",
  p,
  width=5,
  height=4
)

gene_module[
  gene_module$Gene %in% c(
    "BRCA1",
    "AGT",
    "C3",
    "CYP3A4",
    "CXCR4",
    "FOS"
  ),
]

genes6_wgcna <- genes6[
  genes6 %in% colnames(datExpr)
]


genes6_wgcna

ME_cor <- cor(
  datExpr[,genes6_wgcna],
  MEs,
  use="p"
)


ME_cor

library(pheatmap)


pdf(
  "Diagnostic_genes_WGCNA_module_correlation.pdf",
  width=5,
  height=4
)


pheatmap(
  ME_cor,
  cluster_rows=FALSE,
  cluster_cols=FALSE,
  display_numbers=TRUE,
  fontsize=10
)


dev.off()