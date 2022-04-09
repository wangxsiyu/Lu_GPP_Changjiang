library(WangTools)
library(LuGPP)
library(stringr)
library(devtools)
reload(pkgload::inst("WangTools"))
clc()
hs = c(0.95, 0.9, 0.8, 0.5)
ds = c(0.05, 0.1, 0.2, 0.5)
file = './copula_temporal_plotdata_tmaxppet_detrend.RData'
savename = str_replace(file,'_plotdata','')
savename = str_replace(savename, 'copula_', '')
savename = str_replace(savename, '.RData', '')
load(file)

source('./function_copula_MT.R')
for (i in 1:4){
  plt_cpT_veg_month(cpls[[i]]$prob, sprintf("%s_ratio%.2f.png", savename, cpls[[i]]$ratio),hs, ds)
}