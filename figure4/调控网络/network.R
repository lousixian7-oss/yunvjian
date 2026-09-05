a<-read.csv("net_network.csv")
gene<-read.csv("pd&yunajian交集基因.csv")
b=a[a$X2 %in% gene$Intersection_Genes,]
write.csv(b,file="net_network1")
library(openxlsx)
a=read.csv("Maidong.csv",header=T)
b=read.xlsx("Niuxi.xlsx")
c=read.xlsx("Shudihuang.xlsx")
d=read.xlsx("Zhimu.xlsx")
e<-read.csv("type.csv")

for (i in 17:32){
  if (e$node[i] %in% a$X1){
    e[i,3]=1
    e[i,2]="maidong"
  }else if(e$node[i] %in% b$X1){
    e[i,4]=1
    e[i,2]="niuxi"
  }else if(e$node[i] %in% c$X1){
    e[i,6]=1
    e[i,2]="shudihuang"
  }else if(e$node[i] %in% d$X1){
    e[i,5]=1
    e[i,2]="zhimu"
  }
}


for (i in 17:32){
  if (e$node[i] %in% b$X1){
    e[i,4]=1
    e[i,2]="niuxi"
  }}
  
  for (i in 17:32){
    if (e$node[i] %in% c$X1){
      e[i,6]=1
      e[i,2]="shudihuang"
    }}
  
for (i in 17:32){
  if (e$node[i] %in% d$X1){
    e[i,5]=1
    e[i,2]="zhimu"
  }}
write.csv(e,file="type.csv")
