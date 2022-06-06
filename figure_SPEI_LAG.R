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
    nms = arrayfun(function(x)strsplit(x,'_')[[1]][2],nms)
  }
  names(d) = nms
  return(d)
}
dir_lag = "./output/lagSPEI/"
d_lag = W_load(dir_lag)
d_lag = formatname(d_lag)

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


nms = names(d_lag)
nms = nms[!nms %in% c("AVgpp","ECLUEgpp")]
nms = nms[c(5,4,2,6,7,3,1)]


nms_veg = c("forest","open forest","shrub land","paddy field","dry land","grassland")
cols = c("#23FFDC", "#00008F", "#800000", "#FF4A00", "#005AFF", "#ECFF13")

###
{
  png('./figures/FIGURE_all_lagSPEI.png',width = 1600, height = 1920)
  W_figure(7,6,c(0,0.2,0,0))
  par(xpd = NA)

  for (j in 1:7){
    d = d_lag[[nms[j]]]$output_reg




    tet = computes_veg_month(d$coef_partialcor[,])
    tet = lagveg2veglag(tet)
    tet = computes_veg_month(d$coef_partialcor[,])
    tet = lagveg2veglag(tet)
    od = order(d$lag)

    for (i in 1:6){
      # par(mar = c(3,3,5,2) + 0.1,xpd = NA)
      tt = tet[[i]]
      tt = t(tt)
      tt = tt[5:9,]
      if (i == 1){
        ylb = sprintf("%s\nPartial Correlation", nms[j])
      }
      else
      {
        ylb = "";
      }
      boxplot(tt[,od], names = NA,ylim = c(-0.2, 0.4),
              xlab = "Month (lag)", ylab = ylb,
              main = nms_veg[i], col = cols[i], cex.main = 4, cex.axis = 2, cex.lab = 3)
      axis(side = 1,at = 1:length(od),labels = d$lag[od], cex.axis = 1.1, cex.lab = 2.25)
    }

  }
  graphics.off()
}


