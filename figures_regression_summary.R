# library(WangTools)
# library(LuGPP)
# library(pracma)
# library(devtools)
# reload(pkgload::inst("LuGPP"))
#
# clc()
# dir_cor = "./output/parcor/"
# d = W_load(dir_cor)
# nms = names(d)
# nms = nms[nms!=""]
# d = d[names(d) %in% nms]
#
# nms = arrayfun(function(x)strsplit(x,' ~')[[1]][1],nms)
# nms = arrayfun(function(x)strsplit(x,'_')[[1]][2],nms)
# names(d) = nms
# nms = setdiff(nms, c("NIRvgpp","RSgpp","MODISgpp","ECLUEgpp"))
# ng = length(nms)
# month = 1:12#5:9
# nm = length(month)
#
# # correlation map
# file = file.path('./figures',c("summary_cor_tmax.pdf","summary_cor_ppet.pdf"))
# for (vi in 1:2){
#   pvalHD = corHD = matrix(list(), ng, nm)
#   for (gi in 1:ng){
#     for(mi in 1:nm){
#       corHD[[gi, mi]] = d[[nms[gi]]]$output_reg$coef_cor[[vi,month[mi]]]
#       pvalHD[[gi, mi]] = d[[nms[gi]]]$output_reg$pval_cor[[vi,month[mi]]]
#     }
#   }
#   te = computes_veg_month(corHD)
#   plt_veg_months(te, nmsgpp = nms, savedir = file[vi])
#
# # partial correlation map
# file = file.path('./figures',c("parcor_tmax.pdf","parcor_ppet.pdf"))
# for (vi in 1:2){
#   pvalHD = corHD = matrix(list(), ng, nm)
#   for (gi in 1:ng){
#     for(mi in 1:nm){
#       corHD[[gi, mi]] = d[[nms[gi]]]$output_reg$coef_partialcor[[vi,month[mi]]]
#       pvalHD[[gi, mi]] = d[[nms[gi]]]$output_reg$pval_partialcor[[vi,month[mi]]]
#     }
#   }
#   # plt_rasters(corHD, pvalHD, list(breaks = seq(-1,1,0.1)),nmsgpp = nms, savedir = file[vi])
# }
#
#
#
#
# dir_cor = "./output/regression/"
# d = W_load(dir_cor)
# nms = names(d)
# nms = nms[nms!=""]
# d = d[names(d) %in% nms]
#
# nms = arrayfun(function(x)strsplit(x,' ~')[[1]][1],nms)
# # nms = arrayfun(function(x)strsplit(x,'_')[[1]][2],nms)
# names(d) = nms
# nms = setdiff(nms, c("NIRvgpp","RSgpp","MODISgpp","ECLUEgpp"))
# ng = length(nms)
#
# # standardized coef map
# file = file.path('./figures',c("regstd_tmax.pdf","regstd_ppet.pdf"))
# for (vi in 1:2){
#   corHD = pvalHD = matrix(list(), ng, nm)
#   for (gi in 1:ng){
#     for(mi in 1:nm){
#       corHD[[gi, mi]] = d[[nms[gi]]]$output_reg$coef_standreg[[vi,month[mi]]]
#       pvalHD[[gi, mi]] = d[[nms[gi]]]$output_reg$pval_standreg[[vi,month[mi]]]
#     }
#   }
#   plt_rasters(corHD, pvalHD, list(breaks = seq(-1,1,0.1)),nmsgpp = nms, savedir = file[vi])
# }
#
# # relative weight rate map
# file = file.path('./figures',c("rws_tmax.pdf","rws_ppet.pdf"))
# for (vi in 1:2){
#   corHD = pvalHD = matrix(list(), ng, nm)
#   for (gi in 1:ng){
#     for(mi in 1:nm){
#       corHD[[gi, mi]] = d[[nms[gi]]]$output_reg$rws[[vi,month[mi]]]
#       # pvalHD[[gi, mi]] = d[[nms[gi]]]$output_reg$rws[[vi,month[mi]]]
#     }
#   }
#   plt_rasters(corHD, NULL, list(breaks = seq(0,100,1)),nmsgpp = nms, savedir = file[vi])
# }
#
#
# # relative weight rate map
# file = file.path('./figures',c("rw_tmax.pdf","rw_ppet.pdf"))
# for (vi in 1:2){
#   corHD = matrix(list(), ng, nm)
#   for (gi in 1:ng){
#     for(mi in 1:nm){
#       corHD[[gi, mi]] = d[[nms[gi]]]$output_reg$rw[[vi,month[mi]]]
#     }
#   }
#   plt_rasters(corHD, NULL, list(breaks = seq(0,0.5,0.01)),nmsgpp = nms, savedir = file[vi])
# }
#
# # total variance
# file = file.path('./figures',c("r2_reg.pdf"))
# corHD = matrix(list(), ng, nm)
# for (gi in 1:ng){
#   for(mi in 1:nm){
#     corHD[[gi, mi]] = d[[nms[gi]]]$output_reg$rs[[month[mi]]]
#   }
# }
# plt_rasters(corHD, NULL, list(breaks = seq(0,1,0.01)),nmsgpp = nms, savedir = file)
#
#
#
# dir_cor = "./output/reglag/"
# d = W_load(dir_cor)
# nms = names(d)
# nms = nms[nms!=""]
# d = d[names(d) %in% nms]
#
# nms = arrayfun(function(x)strsplit(x,' ~')[[1]][1],nms)
# nms = arrayfun(function(x)strsplit(x,'_')[[1]][3],nms)
# names(d) = nms
# nms = setdiff(nms, c("NIRvgpp","RSgpp","MODISgpp","ECLUEgpp"))
# ng = length(nms)
#
#
# # total variance
# file = file.path('./figures',c("LAG_r2_reg.pdf"))
# R = P = matrix(list(), ng, nm)
# for (gi in 1:ng){
#   for(mi in 1:nm){
#     R[[gi, mi]] = d[[nms[gi]]]$output_reg$rs[[month[mi]]]
#     P[[gi, mi]] = d[[nms[gi]]]$output_reg$ps[[month[mi]]]
#   }
# }
# plt_rasters(R, P, list(breaks = seq(0,1,0.01)),nmsgpp = nms, savedir = file)
#
#
# gppv = c("AVgpp","SDgpp")
# file = file.path('./figures',paste('LAG_TMAX_', gppv,".pdf", sep = ""))
# for (vi in 1:2){
#   R = P = matrix(list(), 6, nm)
#   for (gi in 1:6){
#     for(mi in 1:nm){
#       R[[gi, mi]] = d[[gppv[vi]]]$output_reg$coef_standreg[[gi,month[mi]]]
#       P[[gi, mi]] = d[[gppv[vi]]]$output_reg$pval_standreg[[gi,month[mi]]]
#     }
#   }
#   plt_rasters(R, P, list(breaks = seq(0,1,0.01)),nmsgpp = nms, savedir = file[vi])
# }
#
# gppv = c("AVgpp","SDgpp")
# file = file.path('./figures',paste('LAG_PPET_', gppv,".pdf", sep = ""))
# for (vi in 1:2){
#   R = P = matrix(list(), 6, nm)
#   for (gi in 1:6){
#     for(mi in 1:nm){
#       R[[gi, mi]] = d[[gppv[vi]]]$output_reg$coef_standreg[[gi+6,month[mi]]]
#       P[[gi, mi]] = d[[gppv[vi]]]$output_reg$pval_standreg[[gi+6,month[mi]]]
#     }
#   }
#   plt_rasters(R, P, list(breaks = seq(0,1,0.01)),nmsgpp = nms, savedir = file[vi])
# }
