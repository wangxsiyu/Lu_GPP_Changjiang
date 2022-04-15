library(WangTools)
library(LuGPP)
library(stringr)
library(devtools)
reload(pkgload::inst("WangTools"))
clc()

load('./copula_temporal_CHWIspei.RData')
load('./copula_temporal_CHWIspei_detrend.RData')
library(pracma)
tpd = arrayfun(function(x){
  if (is.null(x)){
    x = ""
  } else{
    x = x$distnames[x$idbest]
  }
  return(x)
}, dry1)
tph = arrayfun(function(x){
  if (is.null(x)){
    x = ""
  } else{
    x = x$distnames[x$idbest]
  }
  return(x)
}, CHWI1)
tpv = arrayfun(function(x){
  if (is.null(x)){
    x = ""
  } else{
    x = x$distnames[x$idbest]
  }
  return(x)
}, veg1)

tpd = as.factor(c(tpd[tpd!= ""]))
tph = as.factor(c(tph[tph!= ""]))
tpv = as.factor(c(tpv[tpv!= ""]))

{
  WD = 1920
  HT = 700
  jpeg("./figs/pieplt.jpg", width = WD, height = HT)
  W_figure(1,3)
  # plot(vmar$fits[[idv]])
  par(cex = 1.5)
  par(mar =  c(8.1, 7.1, 7.1, 2.1))
  pie(table(tpv), main = "GPP")
  # dev.off()
  # jpeg("./figs/qqplot_h.jpg", width = WD, height = HT)
  par(mar =  c(8.1, 7.1, 7.1, 2.1))

  pie(table(tph), main = "CHWI")
  par(mar =  c(8.1, 7.1, 7.1, 2.1))
  # pie(table(tpd), main = "SPEI")

  dev.off()
}

hs = c(0.95, 0.9, 0.8, 0.5)
ds = c(0.05, 0.1, 0.2, 0.5)
file = './copula_temporal_plotdata_tmaxppet_detrend.RData'
file = './copula_temporal_plotdata_CHWIspei_detrend.RData'
savename = str_replace(file,'_plotdata','')
savename = str_replace(savename, 'copula_', '')
savename = str_replace(savename, '.RData', '')
load(file)

source('./function_copula_MT.R')
for (i in 1:4){
  # plt_cpT_veg_month(cpls[[i]]$prob, sprintf("./figs/%s_ratio%.2f.png", savename, cpls[[i]]$ratio),hs, ds)
  plt_cpT_veg_month(cpls[[i]]$prob, sprintf("./figs/%s_ratio%.2f.png", savename, cpls[[i]]$ratio),hs, ds,set = 2)
}


ratios = arrayfun(function(x)x$ratio, cpls)



cplt = rbind(cpls[[1]]$prob[4,],cpls[[2]]$prob[4,],cpls[[3]]$prob[4,],cpls[[4]]$prob[4,])
cplt = t(cplt)
plt_cpT_veg_month(cplt, sprintf("./figs/%s_D_H0.png", savename),ratios,ds,set = 3)

lgd = c("CHWI = 0.5\nSPEI = 0.05","CHWI = 0.5\nSPEI = 0.1",
        "CHWI = 0.5\nSPEI = 0.2","CHWI = 0.5\nSPEI = 0.5")
plt_final(cplt, lgd, "./figs/final_D_H0.png")





cplt = cbind(cpls[[1]]$prob[,4],cpls[[2]]$prob[,4],cpls[[3]]$prob[,4],cpls[[4]]$prob[,4])
plt_cpT_veg_month(cplt, sprintf("./figs/%s_H_D0.png", savename),ratios,hs,set = 4)

lgd = c("CHWI = 0.95\nSPEI = 0.5","CHWI = 0.9\nSPEI = 0.5",
        "CHWI = 0.8\nSPEI = 0.5","CHWI = 0.5\nSPEI = 0.5")
plt_final(cplt, lgd, "./figs/final_H_D0.png")





cplt = matrix(list(),4,4)
for (i in 1:4){
  cplt[[1, i]] = cpls[[i]]$prob[[1,1]]
  cplt[[2, i]] = cpls[[i]]$prob[[4,1]]
  cplt[[3, i]] = cpls[[i]]$prob[[1,4]]
  cplt[[4, i]] = cpls[[i]]$prob[[4,4]]
}
tds =t(matrix(c(0.95,0.05,0.5,0.05,0.95,0.5,0.5,0.5),2,))
plt_cpT_veg_month(cplt, sprintf("./figs/%s_HD.png", savename),ratios,tds,set = 5)

lgd = c("CHWI = 0.95\nSPEI = 0.05","CHWI = 0.5\nSPEI = 0.05",
        "CHWI = 0.95\nSPEI = 0.5","CHWI = 0.5\nSPEI = 0.5")
cols = c("red","gold","pink", "green")
plt_final(cplt, lgd, "./figs/final_HD.png",cols)

plt_final <-function(cplt, lgd, fname = NULL,    cols = c("red","pink", "gold","green")){
  nms_veg = c("Forest","Open Forest","Shrub Land","Paddy Field","Dry Land","Grassland")
  png(filename = fname,width = 1600, height = 1920, units = "px",
      bg = "transparent",  res = NA)
  par(bg = "#ffffff")
  W_figure(6,4, marg = c(0.2, 0.5, 0.12, 1),w = c(1,1,1,1), h = c(1,1,1,1,1,1))
  for (i in 1:6){
    tidx = vegmap == i
    tidx = which(tidx == 1)
    avs = c(0.1,0.2,0.3,0.5)
    sqs = seq(0.5,1,0.01)
    for (j in 1:4){
      tmp = matrix(NA,length(tidx),4)
      tav = matrix(NA, 4, length(sqs))
       for (k in 1:4){
          tmp[,k] = cplt[[k, j]][tidx]
          tav[k,] = quantile(tmp[,k], sqs, na.rm = T)
       }
      par(new = F, cex = 1, xpd = NA)
      par(mar = c(1,1,1,1)*1.5)
      for (k in 1:4){
            if (j == 1 && k == 4){
               ylb = sprintf("%s", nms_veg[i])
            } else{
              ylb  = ""
            }
           if (i == 6 && k == 4){
              xlb = "quantile"
           }else{
             xlb = ""
           }
           if (i == 1 && k == 4){
              man = sprintf("p(GPP < %.2f)\n", avs[j])
           }else{
             man= ""
           }
         plot(sqs, tav[k,], type = 'l', col = cols[k], xlim = c(0.5,1), ylim = c(0,1),
              cex.axis = 1.5, cex.lab = 2, cex.main = 2, main = man,xlab = xlb, ylab = ylb,
              lwd = 5)
         par(new = T)
      }
      if (j == 4 && i == 3){
        par(xpd = NA)
        legend(x = 1, y = 0.5, legend = c(lgd), bty = "n", lty = c(1,1,1,1),
               col = cols, lwd = 5, cex = 2, y.intersp = 2)
      }
      # tav = colMeans(tmp > avs[j], na.rm=T)
      # tav = arrayfun(function(x)mean(tmp[tmp[,x] > avs[j],x], na.rm=T), 1:4)
     # pmin = arrayfun(function(x)min(x[tidx], na.rm = T), cplt)
     # pmax = arrayfun(function(x)max(x[tidx], na.rm = T), cplt)
     # pav = arrayfun(function(x)mean(x[tidx], na.rm = T), cplt)
    }
  }
  dev.off()
}



cplt = matrix(list(),4,3)
for (i in 1:4){
  cplt[[i,1]] = cpls[[1]]$prob[[i,4]]
  cplt[[i,2]] = cpls[[1]]$prob[[4,i]]
}
cplt[[1, 3]] = cpls[[1]]$prob[[1,1]]
cplt[[2, 3]] = cpls[[1]]$prob[[4,1]]
cplt[[3, 3]] = cpls[[1]]$prob[[1,4]]
cplt[[4, 3]] = cpls[[1]]$prob[[4,4]]
lgd = matrix(list(),1,3)
lgd[[1]] = c("CHWI = 0.5\nSPEI = 0.05","CHWI = 0.5\nSPEI = 0.1",
        "CHWI = 0.5\nSPEI = 0.2","CHWI = 0.5\nSPEI = 0.5")
lgd[[2]] = c("CHWI = 0.95\nSPEI = 0.5","CHWI = 0.9\nSPEI = 0.5",
        "CHWI = 0.8\nSPEI = 0.5","CHWI = 0.5\nSPEI = 0.5")
lgd[[3]] = c("CHWI = 0.95\nSPEI = 0.05","CHWI = 0.5\nSPEI = 0.05",
        "CHWI = 0.95\nSPEI = 0.5","CHWI = 0.5\nSPEI = 0.5")




cplt = matrix(list(),4,1)
for (i in 1:1){
  cplt[[1, i]] = cpls[[i]]$prob[[1,1]]
  cplt[[2, i]] = cpls[[i]]$prob[[4,1]]
  cplt[[3, i]] = cpls[[i]]$prob[[1,4]]
  cplt[[4, i]] = cpls[[i]]$prob[[4,4]]
}

lgd = c("Compound","Drought",
        "Heat","Control")

plt_final2(cplt, lgd, fname = 'finalfinal.png')
plt_final2 <-function(cplt, lgd,fname = NULL){
  nms_veg = c("Forest","Open Forest","Shrub Land","Paddy Field","Dry Land","Grassland")
  png(filename = fname,width = 1920, height = 1280, units = "px",
      bg = "transparent",  res = NA)
  par(bg = "#ffffff")
  # mat = matrix(c(1,4,7,10,13,16,0,0,0,0,0,0,2,5,8,11,14,17,0,0,0,0,0,0,3,6,9,12,15,18,0,0,0,0,0,0),6,)
  W_figure(2,3,marg = c(0.1,0.1,0.1,0),w = c(1,1,1), h = c(1,1))
  for (i in 1:6){
    tidx = vegmap == i
    tidx = which(tidx == 1)
    avs = c(0.1,0.2,0.3,0.5)
    sqs = seq(0.5,1,0.01)
    cols = c("red","gold","pink","green")
    for (j in 1:1){
      tmp = matrix(NA,length(tidx),length(cols))
      tav = yyy = matrix(NA, length(cols), length(sqs))
      for (k in 1:length(cols)){
        tmp[,k] = cplt[[k, j]][tidx]
        tav[k,] = quantile(tmp[,k], sqs, na.rm = T)
      }
      xxx = 1 - sqs
      for (k in 1:length(cols)){
        yyy[k,] = tav[k,] #(tav[k,] - tav[length(cols),])/ tav[length(cols),]
      }
      par(new = F, cex = 1, xpd = NA)
      par(mar = c(4,5,3,1))
      for (k in 1:length(cols)){
        if (k == length(cols)){
          man = sprintf("%s", nms_veg[i])
          xlb = ""
          ylb = ""
        } else{
          xlb = ""
          ylb = ""
          man = ""
        }
        plot(xxx, yyy[k,], type = 'l', col = cols[k], xlim = c(0,0.5), ylim = c(0,1),
             cex.axis = 2, cex.lab = 3, cex.main = 3, main = man, xlab = xlb, ylab = ylb,
             lwd = 5)
        par(new = T)
      }
        par(xpd = NA)
        if (i == 1)
        legend('topright', legend = c(lgd), bty = "n", lty = c(1,1,1),
               col = cols, lwd = 3, cex = 3)
      # tav = colMeans(tmp > avs[j], na.rm=T)
      # tav = arrayfun(function(x)mean(tmp[tmp[,x] > avs[j],x], na.rm=T), 1:4)
      # pmin = arrayfun(function(x)min(x[tidx], na.rm = T), cplt)
      # pmax = arrayfun(function(x)max(x[tidx], na.rm = T), cplt)
      # pav = arrayfun(function(x)mean(x[tidx], na.rm = T), cplt)
    }
  }
  dev.off()
}




plt_final3 <-function(cplt, lgd,fname = NULL){
  nms_veg = c("Forest","Open Forest","Shrub Land","Paddy Field","Dry Land","Grassland")
  png(filename = fname,width = 1920, height = 1280, units = "px",
      bg = "transparent",  res = NA)
  par(bg = "#ffffff")
  # mat = matrix(c(1,4,7,10,13,16,0,0,0,0,0,0,2,5,8,11,14,17,0,0,0,0,0,0,3,6,9,12,15,18,0,0,0,0,0,0),6,)
  W_figure(2,3,marg = c(0.1,0.1,0.1,0),w = c(1,1,1), h = c(1,1))
  for (i in 1:6){
    tidx = vegmap == i
    tidx = which(tidx == 1)
    avs = c(0.1,0.2,0.3,0.5)
    sqs = seq(0.5,1,0.01)
    cols = c("red","gold","pink","green")
    for (j in 1:1){
      tmp = matrix(NA,length(tidx),length(cols))
      tav = yyy = matrix(NA, length(cols), length(sqs))
      for (k in 1:length(cols)){
        tmp[,k] = cplt[[k, j]][tidx]
        tav[k,] = quantile(tmp[,k], sqs, na.rm = T)
      }
      xxx = 1 - sqs
      yyy = c(0,0,0)
      for (k in 1:3){
        yyy[k] = xxx[min(which(tav[k,] > 2* tav[4,]))] #(tav[k,] - tav[length(cols),])/ tav[length(cols),]
      }
      par(new = F, cex = 1, xpd = NA)
      par(mar = c(4,5,3,1))
          man = sprintf("%s", nms_veg[i])
          xlb = ""
          ylb = ""
      barplot(yyy, type = 'l', col = cols[k], xlim = c(0,4), ylim = c(0,1),
              cex.axis = 2, cex.lab = 3, cex.main = 3, main = man, xlab = xlb, ylab = ylb,
              lwd = 5)
      par(new = T)

      par(xpd = NA)
      if (i == 1)
        legend('topright', legend = c(lgd), bty = "n", lty = c(1,1,1),
               col = cols, lwd = 3, cex = 3)
      # tav = colMeans(tmp > avs[j], na.rm=T)
      # tav = arrayfun(function(x)mean(tmp[tmp[,x] > avs[j],x], na.rm=T), 1:4)
      # pmin = arrayfun(function(x)min(x[tidx], na.rm = T), cplt)
      # pmax = arrayfun(function(x)max(x[tidx], na.rm = T), cplt)
      # pav = arrayfun(function(x)mean(x[tidx], na.rm = T), cplt)
    }
  }
  dev.off()
}
plt_final3(cplt, lgd, fname = 'finalfinal2.png')
