library(WangTools)
library(LuGPP)
library(pracma)
library(devtools)
reload(pkgload::inst("LuGPP"))

clc()
dir_cor = "./output/parcor/"
d = W_load(dir_cor)
nms = names(d)
nms = nms[nms!=""]
d = d[names(d) %in% nms]

nms = arrayfun(function(x)strsplit(x,' ~')[[1]][1],nms)
nms = arrayfun(function(x)strsplit(x,'_')[[1]][2],nms)
names(d) = nms
nms = setdiff(nms, c("NIRvgpp","RSgpp","MODISgpp","ECLUEgpp"))
# nms = setdiff(nms, c("AVgpp","SDgpp","PCA4gpp"))
ng = length(nms)
month = 1:12#5:9
nm = length(month)

# correlation map
file = file.path('./figures',c("cor_tmax.png","cor_ppet.png"))
file2 = file.path('./figures',c("summary_cor_tmax.png","summary_cor_ppet.png"))
for (vi in 1:2){
  pvalHD = corHD = matrix(list(), ng, nm)
  for (gi in 1:ng){
    for(mi in 1:nm){
      corHD[[gi, mi]] = d[[nms[gi]]]$output_reg$coef_cor[[vi,month[mi]]]
      pvalHD[[gi, mi]] = d[[nms[gi]]]$output_reg$pval_cor[[vi,month[mi]]]
    }
  }
  plt_rasters(corHD[,5:9], pvalHD[,5:9], list(legend = F,breaks = seq(-1,1,0.1), posinset = c(-.1,.3,.1,.3)), nmsgpp = nms, savedir = file[vi])
  te = computes_veg_month(corHD)
  plt_veg_months(te, nmsgpp = nms, savedir = file2[vi])
}

# partial correlation map
file = file.path('./figures',c("parcor_tmax.png","parcor_ppet.png"))
file2 = file.path('./figures',c("summary_parcor_tmax.png","summary_parcor_ppet.png"))
for (vi in 1:2){
  pvalHD = corHD = matrix(list(), ng, nm)
  for (gi in 1:ng){
    for(mi in 1:nm){
      corHD[[gi, mi]] = d[[nms[gi]]]$output_reg$coef_partialcor[[vi,month[mi]]]
      pvalHD[[gi, mi]] = d[[nms[gi]]]$output_reg$pval_partialcor[[vi,month[mi]]]
    }
  }
  plt_rasters(corHD[,5:9], pvalHD[,5:9], list(breaks = seq(-1,1,0.1)),nmsgpp = nms, savedir = file[vi])
  te = computes_veg_month(corHD)
  plt_veg_months(te, nmsgpp = nms, savedir = file2[vi])
}




dir_cor = "./output/regression/"
d = W_load(dir_cor)
nms = names(d)
nms = nms[nms!=""]
d = d[names(d) %in% nms]

nms = arrayfun(function(x)strsplit(x,' ~')[[1]][1],nms)
# nms = arrayfun(function(x)strsplit(x,'_')[[1]][2],nms)
names(d) = nms
nms = setdiff(nms, c("NIRvgpp","RSgpp","MODISgpp","ECLUEgpp"))
ng = length(nms)

# standardized coef map
file = file.path('./figures',c("regstd_tmax.png","regstd_ppet.png"))
file2 = file.path('./figures',c("summary_regstd_tmax.png","summary_regstd_ppet.png"))
for (vi in 1:2){
  corHD = pvalHD = matrix(list(), ng, nm)
  for (gi in 1:ng){
    for(mi in 1:nm){
      corHD[[gi, mi]] = d[[nms[gi]]]$output_reg$coef_standreg[[vi,month[mi]]]
      pvalHD[[gi, mi]] = d[[nms[gi]]]$output_reg$pval_standreg[[vi,month[mi]]]
    }
  }
  plt_rasters(corHD[,5:9], pvalHD[,5:9], list(breaks = seq(-1,1,0.1)),nmsgpp = nms, savedir = file[vi])
  te = computes_veg_month(corHD)
  plt_veg_months(te, nmsgpp = nms, savedir = file2[vi])
}

# relative weight rate map
file = file.path('./figures',c("rws_tmax.png","rws_ppet.png"))
file2 = file.path('./figures',c("summary_rws_tmax.png","summary_rws_ppet.png"))
for (vi in 1:2){
  corHD = pvalHD = matrix(list(), ng, nm)
  for (gi in 1:ng){
    for(mi in 1:nm){
      corHD[[gi, mi]] = d[[nms[gi]]]$output_reg$rws[[vi,month[mi]]]
      # pvalHD[[gi, mi]] = d[[nms[gi]]]$output_reg$rws[[vi,month[mi]]]
    }
  }
  plt_rasters(corHD[,5:9], NULL, list(breaks = seq(0,100,1)),nmsgpp = nms, savedir = file[vi])
  te = computes_veg_month(corHD)
  plt_veg_months(te, nmsgpp = nms, savedir = file2[vi])
}


# # relative weight rate map
# file = file.path('./figures',c("rw_tmax.png","rw_ppet.png"))
# for (vi in 1:2){
#   corHD = matrix(list(), ng, nm)
#   for (gi in 1:ng){
#     for(mi in 1:nm){
#       corHD[[gi, mi]] = d[[nms[gi]]]$output_reg$rw[[vi,month[mi]]]
#     }
#   }
#   plt_rasters(corHD, NULL, list(breaks = seq(0,0.5,0.01)),nmsgpp = nms, savedir = file[vi])
# }

# total variance
file = file.path('./figures',c("r2_reg.png"))
file2 = file.path('./figures',c("summary_r2_reg.png"))
corHD = matrix(list(), ng, nm)
for (gi in 1:ng){
  for(mi in 1:nm){
    corHD[[gi, mi]] = d[[nms[gi]]]$output_reg$rs[[month[mi]]]
  }
}
plt_rasters(corHD, NULL, list(breaks = seq(0,1,0.01)),nmsgpp = nms, savedir = file)
te = computes_veg_month(corHD)
plt_veg_months(te, nmsgpp = nms, savedir = file2)


dir_cor = "./output/reglag/"
d = W_load(dir_cor)
nms = names(d)
nms = nms[nms!=""]
d = d[names(d) %in% nms]

nms = arrayfun(function(x)strsplit(x,' ~')[[1]][1],nms)
nms = arrayfun(function(x)strsplit(x,'_')[[1]][3],nms)
names(d) = nms
nms = setdiff(nms, c("NIRvgpp","RSgpp","MODISgpp","ECLUEgpp"))
ng = length(nms)


# total variance
file = file.path('./figures',c("LAG_r2_reg.png"))
R = P = matrix(list(), ng, nm)
for (gi in 1:ng){
  for(mi in 1:nm){
    R[[gi, mi]] = d[[nms[gi]]]$output_reg$rs[[month[mi]]]
    P[[gi, mi]] = d[[nms[gi]]]$output_reg$ps[[month[mi]]]
  }
}
plt_rasters(R, P, list(breaks = seq(0,1,0.01)),nmsgpp = nms, savedir = file)
te = computes_veg_month(R)
plt_veg_months(te, nmsgpp = nms, savedir = file)


gppv = c("AVgpp","SDgpp")
vi =2
ddd = 6
R = P = matrix(list(), 6, nm)
for (gi in 1:6){
  for(mi in 1:nm){
    R[[gi, mi]] = d[[gppv[vi]]]$output_reg$coef_standreg[[gi+ddd,month[mi]]]
    P[[gi, mi]] = d[[gppv[vi]]]$output_reg$pval_standreg[[gi+ddd,month[mi]]]
  }
}
te = computes_veg_month(R)

plt_veg_months(te, nmsgpp = nms, savedir = file.path('./figures/','lag_tmax.png'))










gppv = c("AVgpp","SDgpp")
file = file.path('./figures',paste('LAG_TMAX_', gppv,".png", sep = ""))
for (vi in 1:2){
  R = P = matrix(list(), 6, nm)
  for (gi in 1:6){
    for(mi in 1:nm){
      R[[gi, mi]] = d[[gppv[vi]]]$output_reg$coef_standreg[[gi,month[mi]]]
      P[[gi, mi]] = d[[gppv[vi]]]$output_reg$pval_standreg[[gi,month[mi]]]
    }
  }
  plt_rasters(R, P, list(breaks = seq(0,1,0.01)),nmsgpp = nms, savedir = file[vi])
}

gppv = c("AVgpp","SDgpp")
file = file.path('./figures',paste('LAG_PPET_', gppv,".png", sep = ""))
for (vi in 1:2){
  R = P = matrix(list(), 6, nm)
  for (gi in 1:6){
    for(mi in 1:nm){
      R[[gi, mi]] = d[[gppv[vi]]]$output_reg$coef_standreg[[gi+6,month[mi]]]
      P[[gi, mi]] = d[[gppv[vi]]]$output_reg$pval_standreg[[gi+6,month[mi]]]
    }
  }
  plt_rasters(R, P, list(breaks = seq(0,1,0.01)),nmsgpp = nms, savedir = file[vi])
}
