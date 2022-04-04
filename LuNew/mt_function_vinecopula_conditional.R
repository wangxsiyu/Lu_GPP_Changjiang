###2021.12.3###
##藤Copula求条件概率##
##安装加载packages##
library(VineCopula)
library(CDVineCopulaConditional)

get_prob_cond <- function(group,gppcond,spaeicond,heatcond){
##建模数据##
  ff<-matrix(0,nrow(group),3)
  ff[,1]<- group[[1]][,1]
  ff[,2]<- group[[1]][,2]
  ff[,3]<- group[[1]][,3]
  
  
##藤Copula建模##
  cvine<-CDVineCondFit(ff,Nx=2,c(1,2,3,4,5,9),rotation=F,treecrit="AIC",type="CVine",selectioncrit="AIC")
  
##条件概率分布和联合概率##
#CVine#
  d<-dim(cvine$Matrix)[1]
  cond1<-ff[,cvine$Matrix[(d+1)-1,(d+1)-1]]
  cond2<-ff[,cvine$Matrix[(d+1)-2,(d+1)-2]]
  condition<-cbind(cond1,cond2)#构造2个条件分布矩阵#
  f1<-CDVineCondSim(cvine,condition)
  #f1[,1]#条件概率分布# 公式10
  
#假设u1是0.6，u2是0.7，u3是0.8，求联合分布概率#
  fff<-matrix(0,1,3)
  fff[1,1]<-gppcond
  fff[1,2]<-spaeicond
  fff[1,3]<-heatcond
  prob = RVinePIT(fff,cvine)[3]#联合概率# 公式11
  #RVinePDF(fff,cvine)#联合概率密度#
  return(prob)
}
