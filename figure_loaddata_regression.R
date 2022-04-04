library(WangTools)
library(LuGPP)
library(pracma)
library(devtools)
reload(pkgload::inst("LuGPP"))
clc()

formatname <- function(d){
  nms = names(d)
  nms = nms[nms!=""]
  d = d[names(d) %in% nms]
  if (all(grepl('~', nms))){
    nms = arrayfun(function(x)strsplit(x,' ~')[[1]][1],nms)
  }
  if (all(grepl('_', nms))){
    nms = arrayfun(function(x)strsplit(x,'_')[[1]][2],nms)
  }
  names(d) = nms
  return(d)
}

dir_cor = "./output/parcor/"
d_cor = W_load(dir_cor)
d_cor = formatname(d_cor)

dir_reg = "./output/regression/"
d_reg = W_load(dir_reg)
d_reg = formatname(d_reg)

dir_lag = "./output/reglag/"
d_lag = W_load(dir_lag)
d_lag = formatname(d_lag)

