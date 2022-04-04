setwd("/Users/mengtian/Desktop/sy_analysis/mt_analysis/function/")
source('../../function/tools_wang.R')
source('../../function/function_wang.R')
source("../001copula/mt_function_vinecopula_conditional.R")

######  1. load data
datadir = '../../../check_datamap/3source_matrix_data/'
data = loadraw(datadir,
               c("AVSDgpp","spei","heatindex")) 

te1 = data$AVSDgpp
te2 = data$spei
te3 = data$heatindex
commonyr = intersect(intersect(te1$yr,te2$yr),te3$yr)
idx1 = which(te1$yr %in% commonyr)
idx2 = which(te2$yr %in% commonyr)
idx3 = which(te3$yr %in% commonyr)

######  2. calculate probability of conditional vine copula
gppperc = 0.3
spaeicond = 
heatcond = 

nx = dim(te1$data[[1]])[1]
ny = dim(te1$data[[1]])[2]
for( xi in 1:nx ){
  for( yi in 1:ny ){
    #############   1 单变量的边缘分布
    load(paste0("../data/",xi,"_",yi,"_1d_marginal_distribution_6-9.RData")) ### group
    
    #############   2 三维Vine Copula函数拟合
    library(VineCopula)
    library(CDVineCopulaConditional)
    
    
    for( mi in 1:length(group) ){ ## month6-9
      emp = group[[mi]]
      gppcond = quantile(emp$Gpp,gppperc)
      
      #####  vine Copula函数拟合和计算条件概率
      prob[[mi]][xi,yi] = get_prob_cond(emp,gppcond,spaeicond,heatcond)
    }
  }
}



# RVM = RVineStructureSelect(emp,family = c(1:10,13,16,18,20,23,24,26:30),indeptest = T,progress = TRUE,se = TRUE,method = 'itau',rotations = TRUE)
# summary(RVM)
# contour(RVM)
# RVineTreePlot(RVM)

###### tau test
# n = nrow(group[[1]])
# tau_thres = 1.96/sqrt(9*n*(n-1)/(2*(2*n+5)))

