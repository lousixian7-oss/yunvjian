

###加载R包
library(tidyverse)
library(GEOquery)
library(tinyarray)
library(AnnoProbe)
library(limma)
library(stringr)
library(writexl)# 导出数据为xlsx格式
library(data.table)


############################1.GEO数据下载##########################


############################2.读取临床信息##########################
####1.读取数据：GSE23586
###1.读取表达矩阵：
setwd("E:\\GEO\\3.GEO数据下载与ID转换\\输入数据")
geoID = 'GSE23586'
gset = getGEO(geoID, destdir=".", 
              AnnotGPL = F, getGPL = F)#后面的是注释文件，为F的就不下载了
gset
setwd("E:\\GEO\\3.GEO数据下载与ID转换\\输出数据")


###2.通过pData函数获取分组信息（获取临床数据）
#临床信息数据：
phe=pData(gset[[1]])
table(phe$title) 
#分为两组：
group_list <- ifelse(str_detect(phe$title, "Healthy_gingiva_CP"), "NC",
                     "PD")#str_detect:字符串中的每个元素都重叠的话返回一个逻辑向量TRUE，这样就不需要把disease state:或者Disease State:列出
table(group_list)
#因子型
group_list = factor(group_list,
                    levels = c("NC","PD"))  
table(group_list)
phe$group = group_list#增加一列分组信息
#保存临床信息
write.table(phe,file=paste0(geoID,".clinic.txt"),sep="\t",quote=F)#quote=F是文件的内容不加上双引号
write.csv(phe, file=paste0(geoID,".clinic.csv"),quote=F)



############################3.提取表达矩阵##########################
#提取表达矩阵
exp <- exprs(gset[[1]])
dim(exp)#看一下dat这个矩阵的维度
exp <- as.data.frame(exp)


############################4.探针ID转换##########################
index = gset[[1]]@annotation#查看平台注释
index
ids <- AnnoProbe::idmap(index)

#对ids的列名进行重命名
colnames(ids) <- c("probe_id", "symbol") #重命名
write.table(ids,file=paste0(geoID,".ids.txt"),sep="\t",quote=F)#quote=F是文件的内容不加上双引号
write_xlsx(ids, paste0(geoID,".ids.xlsx"))



############################5.基因ID注释##########################
#合并探针和基因名
exp = as.data.frame(exp)
exp$probe_id <- rownames(exp)
exp1 <- dplyr::inner_join(ids,exp,by=c("probe_id"="probe_id"))

#去除重复基因
exp2 <- exp1[,-1]
exp <- as.data.frame(limma::avereps(exp2[,-1],ID = exp2$symbol))



############################6.数据标准化##########################
###以下为limma包对于微矩阵数据的标准化
####看数据是否经过归一化处理，若处理过就不用矫正,直接进行下一步id转换
pdf(file=paste0(geoID,".boxplot-数据标准化前.pdf"), width=10,height=8)
par(mar = c(5, 5, 2, 1), oma = c(2, 2, 2, 2))
boxplot(exp,outline=FALSE, notch=F,col=group_list, las=2)#画图
#标准化之前的counts数据——数据量大且分散
dev.off()
range(exp)#查看数据值范围
write.table(exp,file=paste0(geoID,".数据标准化前.txt"),sep="\t",quote=F)
write.csv(exp, file=paste0(geoID,".数据标准化前.csv"),quote=F)


##（避免多次进行normalizeBetweenArrays）
exp <- log2(exp+1)#矫正数据
exp=normalizeBetweenArrays(exp)#校正：数据均值已经接近，不需要再次标准化
pdf(file=paste0(geoID,".boxplot-数据标准化后.pdf"), width=10,height=8)
par(mar = c(5, 5, 2, 1), oma = c(2, 2, 2, 2))
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)#画图
#标准化之后的vst数据——数据量集中
dev.off()
range(exp)#查看数据值范围

#保存表达量矩阵及分组
save(exp,phe,group_list,file = paste0(geoID,"分组及表达矩阵.Rdata"))
write.table(exp,file=paste0(geoID,".数据标准化后.txt"),sep="\t",,quote=F)
write.csv(exp, file=paste0(geoID,".数据标准化后.csv"),quote=F)

#数据备份：
exp_GSE23586 = exp
phe_GSE23586 = phe


##################处理第二个数据##############
####1.读取数据：GSE156993
###1.读取表达矩阵：
setwd("E:\\GEO\\3.GEO数据下载与ID转换\\输入数据")
geoID = 'GSE156993'
gset = getGEO(geoID, destdir=".", 
              AnnotGPL = F, getGPL = F)#后面的是注释文件，为F的就不下载了
gset
setwd("E:\\GEO\\3.GEO数据下载与ID转换\\输出数据")


###2.通过pData函数获取分组信息（获取临床数据）
#临床信息数据：
phe=pData(gset[[1]])
table(phe$title) 
phe = phe[-c(1:18),]
#分为两组：
group_list <- ifelse(str_detect(phe$title, "P"), "PD",
                     "NC")#str_detect:字符串中的每个元素都重叠的话返回一个逻辑向量TRUE，这样就不需要把disease state:或者Disease State:列出
table(group_list)
#因子型
group_list = factor(group_list,
                    levels = c("NC","PD"))  
table(group_list)
phe$group = group_list#增加一列分组信息
#保存临床信息
write.table(phe,file=paste0(geoID,".clinic.txt"),sep="\t",quote=F)#quote=F是文件的内容不加上双引号
write.csv(phe, file=paste0(geoID,".clinic.csv"),quote=F)



############################3.提取表达矩阵##########################
#提取表达矩阵
exp <- exprs(gset[[1]])
dim(exp)#看一下dat这个矩阵的维度
exp <- as.data.frame(exp)


############################4.探针ID转换##########################
index = gset[[1]]@annotation#查看平台注释
index
ids <- AnnoProbe::idmap(index)

#对ids的列名进行重命名
colnames(ids) <- c("probe_id", "symbol") #重命名
write.table(ids,file=paste0(geoID,".ids.txt"),sep="\t",quote=F)#quote=F是文件的内容不加上双引号
write_xlsx(ids, paste0(geoID,".ids.xlsx"))



############################5.基因ID注释##########################
#去除已经删除临床数据的样本
exp <- exp[,phe$geo_accession]
#合并探针和基因名
exp = as.data.frame(exp)
exp$probe_id <- rownames(exp)
exp1 <- dplyr::inner_join(ids,exp,by=c("probe_id"="probe_id"))

#去除重复基因
exp2 <- exp1[,-1]
exp <- as.data.frame(limma::avereps(exp2[,-1],ID = exp2$symbol))


############################6.数据标准化##########################
###以下为limma包对于微矩阵数据的标准化
####看数据是否经过归一化处理，若处理过就不用矫正,直接进行下一步id转换
pdf(file=paste0(geoID,".boxplot-数据标准化前.pdf"), width=10,height=8)
par(mar = c(5, 5, 2, 1), oma = c(2, 2, 2, 2))
boxplot(exp,outline=FALSE, notch=F,col=group_list, las=2)#画图
#标准化之前的counts数据——数据量大且分散
dev.off()
range(exp)#查看数据值范围
write.table(exp,file=paste0(geoID,".数据标准化前.txt"),sep="\t",quote=F)
write.csv(exp, file=paste0(geoID,".数据标准化前.csv"),quote=F)


##（避免多次进行normalizeBetweenArrays）
#exp <- log2(exp+1)#矫正数据
exp=normalizeBetweenArrays(exp)#校正：数据均值已经接近，不需要再次标准化
pdf(file=paste0(geoID,".boxplot-数据标准化后.pdf"), width=10,height=8)
par(mar = c(5, 5, 2, 1), oma = c(2, 2, 2, 2))
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)#画图
#标准化之后的vst数据——数据量集中
dev.off()
range(exp)#查看数据值范围

#保存表达量矩阵及分组
save(exp,phe,group_list,file = paste0(geoID,"分组及表达矩阵.Rdata"))
write.table(exp,file=paste0(geoID,".数据标准化后.txt"),sep="\t",,quote=F)
write.csv(exp, file=paste0(geoID,".数据标准化后.csv"),quote=F)

#数据备份：
exp_GSE156993 = exp
phe_GSE156993 = phe


##################处理第三个数据##############
####1.读取数据：GSE6751
###1.读取表达矩阵：
setwd("E:\\GEO\\3.GEO数据下载与ID转换\\输入数据")
geoID = 'GSE6751'
gset = getGEO(geoID, destdir=".", 
              AnnotGPL = F, getGPL = F)#后面的是注释文件，为F的就不下载了
gset
setwd("E:\\GEO\\3.GEO数据下载与ID转换\\输出数据")


###2.通过pData函数获取分组信息（获取临床数据）
#临床信息数据：
phe=pData(gset[[1]])
table(phe$title) 

phe$title
phe = phe[-c(2,3,4,6,7,9,10,11,13,14,15,17,18,19,21,22,23,25,26,27,29,30,31,33,34,35,37,38,39,41,42,43,45,46,47,49,50,51,53,54,55,57,58,59),]

#分为两组：
group_list <- ifelse(str_detect(phe$title, "Patient"), "PD",
                     "NC")#str_detect:字符串中的每个元素都重叠的话返回一个逻辑向量TRUE，这样就不需要把disease state:或者Disease State:列出
table(group_list)
#因子型
group_list = factor(group_list,
                    levels = c("NC","PD"))  
table(group_list)
phe$group = group_list#增加一列分组信息
#保存临床信息
write.table(phe,file=paste0(geoID,".clinic.txt"),sep="\t",quote=F)#quote=F是文件的内容不加上双引号
write.csv(phe, file=paste0(geoID,".clinic.csv"),quote=F)



############################3.提取表达矩阵##########################
#提取表达矩阵
exp <- exprs(gset[[1]])
dim(exp)#看一下dat这个矩阵的维度
exp <- as.data.frame(exp)


############################4.探针ID转换##########################
index = gset[[1]]@annotation#查看平台注释
index
ids <- AnnoProbe::idmap(index)

#对ids的列名进行重命名
colnames(ids) <- c("probe_id", "symbol") #重命名
write.table(ids,file=paste0(geoID,".ids.txt"),sep="\t",quote=F)#quote=F是文件的内容不加上双引号
write_xlsx(ids, paste0(geoID,".ids.xlsx"))



############################5.基因ID注释##########################
#去除已经删除临床数据的样本
exp <- exp[,phe$geo_accession]
#合并探针和基因名
exp = as.data.frame(exp)
exp$probe_id <- rownames(exp)
exp1 <- dplyr::inner_join(ids,exp,by=c("probe_id"="probe_id"))

#去除重复基因
exp2 <- exp1[,-1]
exp <- as.data.frame(limma::avereps(exp2[,-1],ID = exp2$symbol))


############################6.数据标准化##########################
###以下为limma包对于微矩阵数据的标准化
####看数据是否经过归一化处理，若处理过就不用矫正,直接进行下一步id转换
pdf(file=paste0(geoID,".boxplot-数据标准化前.pdf"), width=10,height=8)
par(mar = c(5, 5, 2, 1), oma = c(2, 2, 2, 2))
boxplot(exp,outline=FALSE, notch=F,col=group_list, las=2)#画图
#标准化之前的counts数据——数据量大且分散
dev.off()
range(exp)#查看数据值范围
write.table(exp,file=paste0(geoID,".数据标准化前.txt"),sep="\t",quote=F)
write.csv(exp, file=paste0(geoID,".数据标准化前.csv"),quote=F)


##（避免多次进行normalizeBetweenArrays）
#exp <- log2(exp+1)#矫正数据
exp=normalizeBetweenArrays(exp)#校正：数据均值已经接近，不需要再次标准化
pdf(file=paste0(geoID,".boxplot-数据标准化后.pdf"), width=10,height=8)
par(mar = c(5, 5, 2, 1), oma = c(2, 2, 2, 2))
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)#画图
#标准化之后的vst数据——数据量集中
dev.off()
range(exp)#查看数据值范围

#保存表达量矩阵及分组
save(exp,phe,group_list,file = paste0(geoID,"分组及表达矩阵.Rdata"))
write.table(exp,file=paste0(geoID,".数据标准化后.txt"),sep="\t",,quote=F)
write.csv(exp, file=paste0(geoID,".数据标准化后.csv"),quote=F)

#数据备份：
exp_GSE6751=exp
phe_GSE6751 = phe

##################处理第四个数据##############
####1.读取数据：GSE10334
###1.读取表达矩阵：

geoID = 'GSE10334'
gset = getGEO(geoID, destdir=".", 
              AnnotGPL = F, getGPL = F)#后面的是注释文件，为F的就不下载了
gset

###2.通过pData函数获取分组信息（获取临床数据）
#临床信息数据：
phe=pData(gset[[1]])
table(phe$title) 

library(dplyr)

# 1. 添加原始行号（数字向量）
phe$orig_row <- 1:nrow(phe)

# 2. 提取患者编号并清理空格
phe$patient <- gsub(".*patient ([0-9]+).*", "\\1", phe$title)
phe$patient <- trimws(phe$patient)

# 3. 标记类型
phe$type <- ifelse(grepl("Unaffected", phe$title, ignore.case = TRUE), "NC", "PD")

# 4. 按患者和类型分组，每组取第一个样本
selected <- phe %>%
  group_by(patient, type) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(as.numeric(patient))

# 5. 提取原始行号（现在完全是数字）
selected_rows <- sort(selected$orig_row)   # 无需 as.numeric

# 6. 用这些行索引筛选 phe
phe <- phe[selected_rows, ]

# 7. 更新分组信息
group_list <- phe$type

# 8. 验证
table(group_list)   # 应相等
View(phe)           # 检查是否每个患者只有1个PD和1个NC
phe$group = group_list#增加一列分组信息
#保存临床信息
write.table(phe,file=paste0(geoID,".clinic.txt"),sep="\t",quote=F)#quote=F是文件的内容不加上双引号
write.csv(phe, file=paste0(geoID,".clinic.csv"),quote=F)



############################3.提取表达矩阵##########################
#提取表达矩阵
exp <- exprs(gset[[1]])
dim(exp)#看一下dat这个矩阵的维度
exp <- as.data.frame(exp)


############################4.探针ID转换##########################
index = gset[[1]]@annotation#查看平台注释
index
ids <- AnnoProbe::idmap(index)

#对ids的列名进行重命名
colnames(ids) <- c("probe_id", "symbol") #重命名
write.table(ids,file=paste0(geoID,".ids.txt"),sep="\t",quote=F)#quote=F是文件的内容不加上双引号
write_xlsx(ids, paste0(geoID,".ids.xlsx"))



############################5.基因ID注释##########################
#去除已经删除临床数据的样本
exp <- exp[,phe$geo_accession]
#合并探针和基因名
exp = as.data.frame(exp)
exp$probe_id <- rownames(exp)
exp1 <- dplyr::inner_join(ids,exp,by=c("probe_id"="probe_id"))

#去除重复基因
exp2 <- exp1[,-1]
exp <- as.data.frame(limma::avereps(exp2[,-1],ID = exp2$symbol))


############################6.数据标准化##########################
###以下为limma包对于微矩阵数据的标准化
####看数据是否经过归一化处理，若处理过就不用矫正,直接进行下一步id转换
pdf(file=paste0(geoID,".boxplot-数据标准化前.pdf"), width=10,height=8)
par(mar = c(5, 5, 2, 1), oma = c(2, 2, 2, 2))
box_col <- ifelse(group_list == "NC", "#00AFBB", "#E7B800")
boxplot(exp,outline=FALSE, notch=F,col=box_col, las=2)
#标准化之前的counts数据——数据量大且分散
dev.off()
range(exp)#查看数据值范围
write.table(exp,file=paste0(geoID,".数据标准化前.txt"),sep="\t",quote=F)
write.csv(exp, file=paste0(geoID,".数据标准化前.csv"),quote=F)


##（避免多次进行normalizeBetweenArrays）
#exp <- log2(exp+1)#矫正数据
exp=normalizeBetweenArrays(exp)#校正：数据均值已经接近，不需要再次标准化
pdf(file=paste0(geoID,".boxplot-数据标准化后.pdf"), width=10,height=8)
par(mar = c(5, 5, 2, 1), oma = c(2, 2, 2, 2))
box_col <- ifelse(group_list == "NC", "#00AFBB", "#E7B800")
boxplot(exp,outline=FALSE, notch=F,col=box_col, las=2)
#标准化之后的vst数据——数据量集中
dev.off()
range(exp)#查看数据值范围

#保存表达量矩阵及分组
save(exp,phe,group_list,file = paste0(geoID,"分组及表达矩阵.Rdata"))
write.table(exp,file=paste0(geoID,".数据标准化后.txt"),sep="\t",,quote=F)
write.csv(exp, file=paste0(geoID,".数据标准化后.csv"),quote=F)

#数据备份：
exp_GSE10334=exp
phe_GSE10334 = phe
#######################7.芯片数据合并#####################
###表达矩阵合并
exp_GSE23586 <- as.data.frame(exp_GSE23586)
exp_GSE23586 <- read.csv("GSE23586.数据标准化后.csv",row.names = 1)
exp_GSE156993 <- as.data.frame(exp_GSE156993)
exp_GSE23586$Gene_name <- rownames(exp_GSE23586)
exp_GSE156993$Gene_name <- rownames(exp_GSE156993)
exp_GSE156993 <- read.csv("GSE156993.数据标准化后.csv",row.names = 1)
exp_GSE6751 <- as.data.frame(exp_GSE6751)
exp_GSE6751$Gene_name <- rownames(exp_GSE6751)
exp_GSE6751 <- read.csv("GSE6751.数据标准化后.csv",row.names = 1)
exp_GSE10334 <- read.csv("GSE10334.数据标准化后.csv",row.names = 1)
exp_GSE10334 <- as.data.frame(exp_GSE10334)
exp_GSE10334$Gene_name <- rownames(exp_GSE10334)
#提取交集基因
merge_eset <- inner_join(exp_GSE23586, exp_GSE156993, by = "Gene_name") #suffix = c(".data1", ".data2")：当两个数据框中存在相同名称的列时，inner_join 会在这些列名后添加后缀。这里指定了 .data1 和 .data2 作为后缀。
# merge_eset <- inner_join(merge_eset, exp_GSE6751, by = "Gene_name") 
merge_eset <- inner_join(merge_eset, exp_GSE10334, by = "Gene_name") 
#Gene_name作为行名
n=which(colnames(merge_eset) == "Gene_name") #找到Gene_name列的所在位置
exp <- as.data.frame(limma::avereps(merge_eset[,-n],ID = merge_eset$Gene_name))
dim(exp) 
##保存
write.table(exp,file="合并后的数据.txt",sep="\t",quote=F)
write.csv(exp, file="合并后的数据.csv",quote=F)


###临床信息合并
phe_GSE23586 <- phe_GSE23586[,c(2,39)]
phe_GSE156993 <- phe_GSE156993[,c(2,44)]
# phe_GSE6751 <- phe_GSE6751[,c(2,39)]
phe_GSE10334 <- phe_GSE10334[,c(2,39)]
#使用rbind()函数进行按行合并
group <- as.data.frame(rbind(phe_GSE23586,phe_GSE156993,phe_GSE10334)) 
#验证group与表达矩阵的列名是否一致
identical(group$geo_accession,colnames(exp))

##设置分组信息
group_list <- group$group
group_list = factor(group_list,
                    levels = c("NC","PD"))  
group_list


###5.去批次效应
#5.1去批次前处理
nrow(phe_GSE23586)
nrow(phe_GSE156993)
nrow(phe_GSE6751)
nrow(phe_GSE10334)
#制作batchType文件
batchType=c(rep(1,nrow(phe_GSE23586)),
            rep(2,nrow(phe_GSE156993)),
            # rep(3,nrow(phe_GSE6751)),
            rep(3,nrow(phe_GSE10334))
)  

#制作mod文件
group_list
mod = model.matrix(~group_list)
mod

#5.2绘制去批次前的箱线图和PCA图，以便后续进行对比
pdf(file = "combat前-boxplot.pdf",width = 10,height = 10)
par(mar = c(5, 5, 2, 1), oma = c(2, 2, 2, 2))
boxplot(exp,outline=F, notch=T,col=group_list, las=2)#绘制去批次前箱线图
dev.off()
#绘制去批次前PCA图
library(FactoMineR)
library(factoextra)
exp <- read.csv(
  "合并后的数据.csv",
  row.names=1,
  check.names=FALSE
)
dat.pca <- PCA(as.data.frame(t(exp)), graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                         col.ind = group_list,
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, 
                         legend.title = "Groups")
pca_plot
ggsave(plot = pca_plot,filename ="combat前-prenormal_PCA.pdf",width = 8,height = 8)
dat.pca.before <- PCA(t(exp), graph = FALSE) # t()把样本转成行
# 以校正前按批次着色为例
p1 <- fviz_pca_ind(dat.pca.before,
                   geom.ind = "point",
                   col.ind = as.factor(batchType),
                   palette = c("#FF6B6B", "#E7B800", "#00AFBB", "#96CEB4"),
                   addEllipses = TRUE,      # 画椭圆，清晰圈出4个数据集的范围
                   ellipse.level = 0.95,
                   pointshape = 19,         # 实心圆点
                   pointsize = 2,           # 点调小一点避免完全覆盖
                   alpha.ind = 0.7,         # 增加透明度！重叠越多颜色越深
                   legend.title = "Dataset")
print(p1)


#5.3进行去批次处理
library(sva)
library(limma)
#sva包ComBat函数去除批次:
exp_ComBat=ComBat(dat=exp, batch=batchType, #使用ComBat法去批次
                  mod=mod, par.prior=TRUE)
range(exp_ComBat)

#绘制去批次后PCA图
pdf(file = "combat后-boxplot.pdf",width = 10,height = 10)
par(mar = c(5, 5, 2, 1), oma = c(2, 2, 2, 2))
boxplot(exp_ComBat,outline=F, notch=T,col=group_list, las=2)#绘制去批次后箱线图
dev.off()
exp_ComBat=normalizeBetweenArrays(exp_ComBat)
boxplot(exp_ComBat,outline=F, notch=T,col=group_list, las=2)#绘制去批次后箱线图
dat.pca2 <- PCA(as.data.frame(t(exp_ComBat)), graph = FALSE)
pca_plot2 <- fviz_pca_ind(dat.pca2,
                          geom.ind = "point",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                          col.ind = group_list,
                          palette = c("#00AFBB", "#E7B800"),
                          addEllipses = TRUE, 
                          legend.title = "Groups")
pca_plot2
ggsave(plot = pca_plot2,filename ="combat后_PCA.pdf",width = 8,height = 8)
# 以校正前按批次着色为例
# 基于你已经校正好的 exp_ComBat 数据
dat.pca.after_batch <- PCA(as.data.frame(t(exp_ComBat)), graph = FALSE)

pca_after_batch <- fviz_pca_ind(at.pca.after_batch,
                                geom.ind = "point",
                                col.ind = as.factor(batchType),  # 按4个数据集着色
                                palette = c("#FF6B6B", "#E7B800", "#00AFBB", "#96CEB4"),
                                addEllipses = TRUE,
                                ellipse.level = 0.95,
                                pointshape = 19,
                                pointsize = 2,
                                alpha.ind = 0.7,
                                legend.title = "Dataset") +
  ggtitle("PCA After ComBat")

print(pca_after_batch)
ggsave(plot = pca_after_batch, filename = "combat后_PCA_by_batch.pdf", width = 8, height = 8)
#保存
write.table(exp_ComBat,file="去除批次效应后表达矩阵.txt",sep="\t",quote=F)
write.csv(exp_ComBat, file="去除批次效应后表达矩阵.csv",quote=F)

#保存矩阵数据
save(exp_ComBat,group,group_list,file = "去除批次效应后表达矩阵.Rdata")




