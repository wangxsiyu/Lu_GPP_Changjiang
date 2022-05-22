library(WangTools)
library(LuGPP)
library(pracma)
library(ppcor)
library(rwa)
library(devtools)
reload(pkgload::inst("LuGPP"))
filehead = 'lagSPEI_'
option = c(1)
clc()
datadir = './data/gpp/'
data = loadraw(datadir)
var = loadraw('./data/','tmax')
spei = loadraw('./data/spei/')
nspei = length(spei)
lag = names(spei)
lag = strsplit(lag, 'spei')
lag = arrayfun(function(x){as.numeric(x[2])}, lag)
veg = names(data)
veg = veg[veg%in%c("Pgpp","SIFgpp","VPMgpp","Musyqgpp")]
months = 1:12

for (vegi in 1:length(veg)) {
  tfile = sprintf("%s%s.RData", filehead, veg[vegi])
  tfilename = file.path("./output/lagSPEI",tfile)

  output_reg = list()
  output_reg$coef_partialcor = output_reg$pval_partialcor = matrix(list(), length(spei), 12)
  for (si in 1:length(spei)){
    tveg = c(var, data[veg[vegi]])
    tveg$spei = spei[[si]]
    tveg = merge_vars(tveg)
    tstr = paste(veg[vegi], "~ tmax + spei")
    for (mi in 1:length(months)){

      print(sprintf('spei %d/%d, mt%d/%d', si, length(spei), mi, length(months)))
      mt = months[mi]
      tvegm = get_month(tveg, mt)
      nx = dim(tvegm[[veg[vegi]]][[1]])[1]
      ny = dim(tvegm[[veg[vegi]]][[1]])[2]
      output_reg$coef_partialcor[[si, mi]] = matrix(NA, nx, ny)
      output_reg$pval_partialcor[[si, mi]] = matrix(NA, nx, ny)
      nt = length(tvegm$yyyymm)
      for (xi in 1:nx){
        if (xi %% 10 == 0)
        {        print(sprintf('%d/%d', xi, nx))}
        for (yi in 1:ny){
          te = matrix(NA, nt , 3)
          te[,1] = arrayfun(function(x)x[xi, yi], tvegm[[veg[vegi]]])
          te[,2] = arrayfun(function(x)x[xi, yi], tvegm$tmax)
          te[,3] = arrayfun(function(x)x[xi, yi], tvegm$spei)
          te = as.data.frame(te)
          names(te) = c(veg[vegi], 'tmax','spei')

          idxall =colMeans(is.na(t(te))) ==0
          if ((sum(idxall) < 10) || (sd(te[,1], na.rm = T)==0)){

          }
          else{
            te = te[idxall,]
            if (1 %in% option){
              # partial correlation
              ttt <- try({
                tc = pcor.test(te[,1], te$spei, te$tmax,'pearson')
                output_reg$coef_partialcor[[si, mi]][xi,yi] = tc$estimate
                output_reg$pval_partialcor[[si, mi]][xi,yi] = tc$p.value
              }, silent = T)
            }
          }
        }
      }
    }
  }
  output_reg$vegname = veg[vegi]
  output_reg$lag = lag
  save(output_reg, file = tfilename)
}



