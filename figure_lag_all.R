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
dir_lag = "./output/reglag/"
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
nms = nms[!nms %in% "AVgpp"]
nms = nms[c(5,4,2,6,7,3,1)]


nms_veg = c("forest","open forest","shrub land","paddy field","dry land","grassland")
cols = c("#23FFDC", "#00008F", "#800000", "#FF4A00", "#005AFF", "#ECFF13")

# library(devEMF)
# emf('./figures/FIGURE_all_lag.emf')
dev.off()
{
  png('./figures/FIGURE_all_lag_ppet.png',width = 1600, height = 1920)
  W_figure(7,6,c(0,0.2,0,0))
  par(xpd = NA)

  for (j in 1:7){
    d = d_lag[[nms[j]]]$output_reg

    tet = computes_veg_month(d$coef_standreg[1:6,])
    tep = computes_veg_month(d$coef_standreg[(1:6)+6,])

    tet = lagveg2veglag(tet)
    tep = lagveg2veglag(tep)
    for (i in 1:6){
      tt = t(tet[[i]])
      tp = t(tep[[i]])
      if (i == 1){
        ylb = sprintf("%s\nStandard coefficients", nms[j])
      }
      else
      {
        ylb = "";
      }
      boxplot(tp[5:9,], at = 1:6, names = 0:5,ylim = c(-0.2, 0.5), xlab = "Month (lag)", ylab = ylb,
              main = nms_veg[i], col = cols[i], cex.main = 4, cex.axis = 2, cex.lab = 3)
    }

  }
  graphics.off()
}
{
  png('./figures/FIGURE_all_lag_tmax.png',width = 1600, height = 1920)
  W_figure(7,6,c(0,0.2,0,0))
  par(xpd = NA)

  for (j in 1:7){
    d = d_lag[[nms[j]]]$output_reg

    tet = computes_veg_month(d$coef_standreg[1:6,])
    tep = computes_veg_month(d$coef_standreg[(1:6)+6,])

    tet = lagveg2veglag(tet)
    tep = lagveg2veglag(tep)
    for (i in 1:6){
      tt = t(tet[[i]])
      tp = t(tep[[i]])
      if (i == 1){
        ylb = sprintf("%s\nStandard coefficients", nms[j])
      }
      else
      {
        ylb = "";
      }
      boxplot(tt[5:9,], at = 1:6, names = 0:5,ylim = c(-0.4, 0.8), xlab = "Month (lag)", ylab = ylb,
              main = nms_veg[i], col = cols[i], cex.main = 4, cex.axis = 2, cex.lab = 3)
    }

  }
  graphics.off()
}






{
  png('./figures/FIGURE_all_regcor.png',width = 1600, height = 800)

  tout = matrix(list(),2,6)
  for (i in 1:6){
    tout[[1,i]] = matrix(NA, 7,12)
    tout[[2,i]] = matrix(NA, 7,12)
  }
  for (j in 1:7){
    d = d_cor[[nms[j]]]$output_reg

    tet = computes_veg_month(d$coef_partialcor[1,])
    tet2 = computes_veg_month(d$coef_partialcor[2,])
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

  W_figure(2,6,c(0,0,0.2,0))
  par(xpd = NA)
  strtype = c("Tmax","P-Pet")
  for (k in 1:2){
    if (k == 1){
      ylm = c(0,1)
    }    else{
      ylm = c(-0.5, 0.5)
    }
    for (i in 1:6){
      if (i == 1){
        ylb = sprintf("%s\nPartial Correlation", strtype[k])
      }else {
        ylb = "Partial Correlation"
      }
      par(new = F)
      for (j in 1:7){
        plot(1:12, tout[[k,i]][j,], type = 'l', col = cols[j], cex.main = 2, cex.axis = 1.2, cex.lab = 1.5,
             xlab = "Month", ylab =ylb, main = nms_veg[i], ylim = ylm)
        par(new = T)
      }
    }
  }

  graphics.off()
}

