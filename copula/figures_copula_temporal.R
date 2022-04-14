library(WangTools)
library(LuGPP)
library(stringr)
library(devtools)
reload(pkgload::inst("WangTools"))
clc()

load('./copula_temporal_tmaxppet_detrend.RData')



hs = c(0.95, 0.9, 0.8, 0.5)
ds = c(0.05, 0.1, 0.2, 0.5)
file = './copula_temporal_plotdata_tmaxppet_detrend.RData'
file = './copula_temporal_plotdata_heatspei_detrend.RData'
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

lgd = c("HEAT = 0.5, SPEI = 0.05","HEAT = 0.5, SPEI = 0.1",
        "HEAT = 0.5, SPEI = 0.2","HEAT = 0.5, SPEI = 0.5")
plt_final(cplt, lgd, "./figs/final_D_H0.png")





cplt = cbind(cpls[[1]]$prob[,4],cpls[[2]]$prob[,4],cpls[[3]]$prob[,4],cpls[[4]]$prob[,4])
plt_cpT_veg_month(cplt, sprintf("./figs/%s_H_D0.png", savename),ratios,hs,set = 4)

lgd = c("HEAT = 0.95, SPEI = 0.5","HEAT = 0.9, SPEI = 0.5",
        "HEAT = 0.8, SPEI = 0.5","HEAT = 0.5, SPEI = 0.5")
plt_final(cplt, lgd, "./figs/final_H_D0.png")





cplt = matrix(list(),4,4)
for (i in 1:4){
  cplt[[1, i]] = cpls[[i]]$prob[[1,1]]
  cplt[[2, i]] = cpls[[i]]$prob[[4,1]]
  cplt[[3, i]] = cpls[[i]]$prob[[1,4]]
  cplt[[4, i]] = cpls[[i]]$prob[[4,4]]
}
tds =t(matrix(c(0.95,0.05,0.95,0.5,0.5,0.95,0.5,0.5),2,))
plt_cpT_veg_month(cplt, sprintf("./figs/%s_HD.png", savename),ratios,tds,set = 5)

lgd = c("HEAT = 0.95, SPEI = 0.05","HEAT = 0.95, SPEI = 0.5",
        "HEAT = 0.5, SPEI = 0.05","HEAT = 0.5, SPEI = 0.5")
plt_final(cplt, lgd, "./figs/final_HD.png")

plt_final <-function(cplt, lgd, fname = NULL){
  nms_veg = c("Forest","Open Forest","Shrub Land","Paddy Field","Dry Land","Grassland")
  png(filename = fname,width = 1280, height = 1920, units = "px",
      bg = "transparent",  res = NA)
  par(bg = "#ffffff")
  W_figure(6,4, marg = c(0.2, 0.2, 0.12, 0.05),w = c(1,1,1,1), h = c(1,1,1,1,1,1))
  for (i in 1:6){
    tidx = vegmap == i
    tidx = which(tidx == 1)
    avs = c(0.1,0.2,0.3,0.5)
    sqs = seq(0.5,1,0.01)
    cols = c("red1","gold1","yellow2","green")
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
      legend(x = 0.5, y = 1, legend = c(lgd), bty = "n", lty = c(1,1,1,1), col = cols, lwd = 5)
      # tav = colMeans(tmp > avs[j], na.rm=T)
      # tav = arrayfun(function(x)mean(tmp[tmp[,x] > avs[j],x], na.rm=T), 1:4)
     # pmin = arrayfun(function(x)min(x[tidx], na.rm = T), cplt)
     # pmax = arrayfun(function(x)max(x[tidx], na.rm = T), cplt)
     # pav = arrayfun(function(x)mean(x[tidx], na.rm = T), cplt)
    }
  }
  dev.off()
}
