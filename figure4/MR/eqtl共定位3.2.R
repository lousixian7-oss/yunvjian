library(data.table)
library(coloc)
library(TwoSampleMR)

GWASname <- "periodontal.out.zip"  ####注意需要改 最后汇总的列名

library(vroom)
a<-vroom(file.choose())
a<-read.table(file="eqtlgen-lite.txt",header=T)

# c<-unique(b$gene)
# write.csv(c,file="gene.csv")
c<-read.csv("gene.csv")
write.table(c,file="gene.txt",quote = FALSE,row.names = FALSE, col.names = FALSE)
c<-read.table("gene.txt")

d<-a[a$GENE==AANAT]

for(i in c){
  if(i %in% a){
  d<-a[a$GENE==i]
  write.table(b,file=paste(i,".txt", sep =""))
  }
}
# 读取 GWAS 数据
data <- fread(GWASname, header = TRUE) ##gwas数据
head(data)
data = data[!duplicated(data$SNP),]
data <- data[,c('SNP','chr','pos','effect_allele','other_allele','beta','se','eaf','pval','nearest_genes')]
colnames(data) <- c('SNP','CHR','BP','A1','A2','BETA','SE','FREQ','P','GENE')

workdir <- getwd()  # 获取当前工作目录的路径
workdir <- paste0(getwd(), "/eqtl")  
workdir  # 打印新的工作目录路径

files <- list.files(path = workdir, pattern = ".txt", full.names = TRUE) ####（原始数据）

result_coloc <- list()

gene_chr <- fread("基因起始-37.txt")###蛋白共定位的染色体和位置，我已经做好

for (filename in files) {
  # 从文件中读取数据
  eqtl <- fread(filename)
  # eqtl <- fread(files[1])
  # head(eqtl)
  eqtl <- eqtl[,c('SNP','Chr','BP','A1','A2','b','SE','Freq','p','GENE')]#######需要根据你暴露数据列名修改
  colnames(eqtl) <- c('SNP','CHR','BP','A1','A2','BETA','SE','FREQ','P','GENE')
  
   #  files[1]
     #filename="ACTG1.txt" #举例文件名
  #去除beta值不为数字或无限的行
  eqtl <- eqtl[!is.na(as.numeric(eqtl$BETA)) & is.finite(as.numeric(eqtl$BETA)), ]
   

  
  # 使用strsplit函数提取第一个_之前的字符
  gene <- tools::file_path_sans_ext(basename(filename))
  #gene <- eqtl$GENE[1]
  gene_chr_subset <- gene_chr[toupper(gene_chr$GENE) == toupper(gene), "CHR"]
  gene_chr_subset <- as.numeric(gene_chr_subset)
  #gene_chr_subset <- paste0("chr", gene_chr_subset)
  gene_start <- gene_chr[toupper(gene_chr$GENE) == toupper(gene), "START"]
  gene_end <- gene_chr[toupper(gene_chr$GENE) == toupper(gene), "END"]
  gene_start <- as.numeric(gene_start)
  gene_end <- as.numeric(gene_end)
  #要删除rsids列中不是以"rs"开头的行
  eqtl <- eqtl[eqtl$SNP %like% "^rs", ]

  

  my_eqtl <- eqtl[eqtl$BP > (gene_start - 100000) & eqtl$BP < (gene_end + 100000), ]
  my_eqtl <- my_eqtl[!duplicated(my_eqtl$SNP), ]
  
  merged_data <- merge(my_eqtl, data, by = "SNP") 
  
  my_eqtl <- merged_data[,c('SNP','FREQ.y','BETA.x','SE.x','P.x')]
  colnames(my_eqtl) <- c('snp','MAF','beta','SE','pvalues')
  
  gwas <- merged_data[,c('SNP','P.y','BETA.y','SE.y')]
  colnames(gwas) <- c('snp','pvalues','beta','se')
  gwas$varbeta <- (gwas$se)^2
  
  ##计算MAF值
  my_eqtl$MAF <- ifelse(my_eqtl$MAF < 0.5, my_eqtl$MAF, 1 - my_eqtl$MAF)

  # 将 P 值为 0 的行替换为 1e-310 否则0值报错

  my_eqtl$pvalues[my_eqtl$pvalues == 0] <- 1e-310
  
  
  my_eqtl <- as.list(my_eqtl)
  my_eqtl[['type']] <- 'quant'
  my_eqtl[['N']] <- 31684
  
  gwas <- as.list(gwas)
  gwas[['type']] <- 'cc'
  gwas[['N']] <- 406876
  
  result <- coloc.abf(dataset1 = gwas, dataset2 = my_eqtl)

  result_summary_df <- as.data.frame(result$summary)
  result_summary <- as.data.frame(t(result_summary_df))
  result_summary$Gene <- toupper(gene)
  result_summary$GWAS <- GWASname
  df <- result_summary[, c("GWAS", "Gene",  setdiff(names(result_summary), c("GWAS", "Gene" )))]
  result_coloc[[filename]] <- df
}

# 汇总所有共定位结果
final_result <- rbindlist(result_coloc, fill = TRUE)

# 将最后的汇总结果写入文件
fwrite(final_result, file = "result_coloc2.csv")
###老铁 完事了！