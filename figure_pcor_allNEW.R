library(WangTools)
library(LuGPP)
library(pracma)
library(devtools)
reload(pkgload::inst("LuGPP"))
reload(pkgload::inst("WangTools"))
clc()

formatname <- function(d){
  nms = names(d)
  nms = nms[nms!=""]
  d = d[names(d) %in% nms]
  if (all(grepl('~', nms))){
    nms = arrayfun(function(x)strsplit(x,' ~')[[1]][1],nms)
  }
  if (all(grepl('_', nms))){
    nms = arrayfun(function(x)strsplit(x,'_')[[1]][3],nms)
  }
  names(d) = nms
  return(d)
}
dir_cor = "./output/reglagPcor/"
d_cor = W_load(dir_cor)
d_cor = formatname(d_cor)

lagveg2veglag <- function(te){
  out = matrix(list(),1,6)
  nn = length(te)
  for(i in 1:6){
    out[[i]] = matrix(NA,nn, 12)
    for (j in 1:nn){
      out[[i]][j,] = te[[j]][i,]
    }
  }
  return(out)
}


nms = names(d_cor)
nms = nms[!nms %in% c("AVgpp","ECLUEgpp")]
nms = nms[c(5,4,2,6,7,3,1)]


nms_veg = c("forest","open forest","shrub land","paddy field","dry land","grassland")
cols = c("#23FFDC", "#00008F", "#800000", "#FF4A00", "#005AFF", "#ECFF13")

# library(devEMF)
# emf('./figures/FIGURE_all_lag.emf')
dev.off()
{
  png('./figures/FIGURE_all_PCOR_new.png',width = 1600, height = 900)
  tout = matrix(list(),2,6)
  for (i in 1:6){
    tout[[1,i]] = matrix(NA, 7,12)
    tout[[2,i]] = matrix(NA, 7,12)
  }
  for (j in 1:7){
    d = d_cor[[nms[j]]]$output_reg

    tet = computes_veg_month(d$parcor_predH[1,])
    tet2 = computes_veg_month(d$parcor_predD[1,])
    tt = tt2 = matrix(NA, 6, 12)
    for (i in 1: 12){
      tt[,i] = tet[[i]]
      tt2[,i] = tet2[[i]]
    }
    for (i in 1:6){
      tout[[1,i]][j,] = tt[i,]
      tout[[2,i]][j,] = tt2[i,]
    }
  }

  W_figure(2,6,c(0,0.2,0.2,0))
  par(xpd = NA)
  strtype = c("Tmax","P-Pet")
  for (k in 1:2){
    if (k == 1){
      ylm = c(-0.1,1)
    }    else{
      ylm = c(-0.1,1)
    }
    for (i in 1:6){
      if (i == 1){
        ylb = sprintf("%s\nPartial Correlation", strtype[k])
      }else {
        ylb = "Partial Correlation"
      }
      par(new = F)
      boxplot(tout[[k,i]][2:7,], cex.main = 4, cex.axis = 2, cex.lab = 3,
              xlab = "Month", ylab =ylb, main = nms_veg[i], ylim = ylm, xlim = c(0.5, 12.5))
      par(new = T)
      plot(1:12,tout[[k,i]][1,], type = 'l',lty = 1, col = "blue", lwd = 3,cex.main = 4, cex.axis = 2, cex.lab = 3,
           xlab = "Month", ylab =ylb, main = nms_veg[i], ylim = ylm, xlim = c(0.5, 12.5))

    }
  }

  graphics.off()
}

