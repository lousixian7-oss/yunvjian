setwd("E:/Yunvjian/yunvjian")

library(openxlsx)
a=read.csv("Maidong.csv",header=T)
b=read.xlsx("Niuxi.xlsx")
c=read.xlsx("Shudihuang.xlsx")
d=read.xlsx("Zhimu.xlsx")
e=rbind(a,b,c,d)
write.csv(e,file="net_network.csv")
