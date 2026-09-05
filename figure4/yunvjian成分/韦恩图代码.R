# 加载必要的包
library(ggvenn)
library(eulerr)
library(scales)
library(dplyr)  
library(data.table)
library(readxl)
# 读取数据
data <- read_xls("yunvjian.xls")
library(dplyr)

# 将所有 NA 替换为空字符串
data[is.na(data)] <- ""


# 提取两列基因作为集合
comgene <- list(
  'Mai Dong' = unique(data$maidong),
  'Niu Xi' = unique(data$niuxi),
  'Shu Di Huang  ' = unique(data$shudihuang),
  'Zhi Mu' = unique(data$zhimu),
  "Shi Gao"=unique(data$shigao)
)
library (VennDiagram)  
venn.diagram(x=comgene,
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.5, #透明度
             
             lwd=1,lty=1,col=c("#6699CC", "#9966CC","#669966","#EB7E60","#F5AA61"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fontfamily = "Times New Roma",  # 字体
             
             fill=c("#6699CC", "#9966CC","#669966","#EB7E60","#F5AA61"), # 填充色 配色https://www.58pic.com/
             
             category.names = c('Mai Dong',
                                'Niu Xi' ,
                                'Shu Di Huang    ' ,
                                'Zhi Mu' ,
                                "Shi Gao"),#标签名
             
             cat.dist = c(0.2, 0.2, 0.2, 0.2, 0.2), # 标签距离圆圈的远近
             
             cat.pos = c(0, -10, 240, 120, 20), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 0.5, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col=c("#6699CC", "#9966CC","#669966","#EB7E60","#F5AA61") ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             cat.fontfamily = "Times New Roma",  # 字体
             
             xrotation = 2,  # 1 2 3 旋转确定大打头数据集
             
             filename='5组.png',# 文件保存
             
             output=TRUE,
             
             imagetype="png",  # 类型（tiff png svg）
             
             resolution = 400,  # 分辨率
             
             compression = "lzw",  # 压缩算法
             
             ext.text = T, # 增加指示线和标签
             
             ext.percent = c(0.1,0.1,0.1), # 出现指示线和标签的条件
             
             ext.dist = c(0.01,0.01),
             
             ext.length = 0.8)
# 计算交集基因
com_genes <- union(comgene$`Mai Dong`, comgene$`Niu Xi`)
com_genes <- union(com_genes, comgene$`Shu Di Huang`)
com_genes <- union(com_genes,comgene$`Zhi Mu`)
#保存为 CSV
{
write.csv(data.frame(Intersection_Genes = com_genes),
          "并集基因.csv", row.names = FALSE)
}

#普通韦恩图
p1 <- ggvenn(
  comgene,
  show_percentage = F,
  show_elements = FALSE,
  label_sep = ",",
  digits = 2,
  stroke_color = NA,  # 去除边框
  fill_color = c("#6699CC", "#9966CC","#669966","#EB7E60"), 
  set_name_color = c("#6699CC", "#9966CC","#669966","#EB7E60"),
  text_color = "black",  # 调整文字颜色
  text_size = 5          # 增大文字
)
print(p1)

# 绘制比例图
p2 <- euler(comgene)
plot(
  p2,
  labels = list(col = "black", font = 1, cex = 1.0),  # 调整字体大小和颜色
  edges = NULL,  # 去除边框
  fills = c("#C2E0F7", "#A4CBA8","#FFD39B","#FF9999"),  
  alpha = 0.5,
  quantities = list(cex = 1.0, col = 'black')  # 设置数量文字更突出
)

