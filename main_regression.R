library(WangTools)
library(LuGPP)
library(pracma)
library(ppcor)
library(rwa)
filehead = 'parcor_'
option = c(1,2)
clc()
datadir = './data/gpp/'
data = loadraw(datadir)
vars = loadraw('./data/',c('tmax','ppet'))
## veg vs var
veg = names(data)
var = names(vars)
nvar = length(var)
## get the data for the variable first
dvar = get_common(vars)
months = 1:12
for (vegi in 1:length(veg)) {
  tveg = c(dvar, data[vegi])
  tveg = merge_vars(tveg)
  tstr = paste(veg[vegi], "~",paste(var, collapse = "+"))
  tfile = sprintf("%s%s.RData", filehead, tstr)
  tfilename = file.path("./output/regression/",tfile)
  # tfilemark = file.path("./output/regression/",sprintf("%s_mark.RData", tstr))
  if (file.exists(tfilename)){
    next;
  }
  # else {
  #   save(tfilemark, file = tfilemark)
  # }
  output_reg = list()
  output_reg$rs = matrix(list(),1,12)
  output_reg$rw = output_reg$rws = output_reg$sn = matrix(list(), nvar, 12)
  output_reg$coef_cor = output_reg$pval_cor = matrix(list(), nvar, 12)
  output_reg$coef_partialcor = output_reg$pval_partialcor = matrix(list(), nvar, 12)
  output_reg$coef_standreg = output_reg$pval_standreg = matrix(list(), nvar, 12)
  for (mi in 1:length(months)){
    print(sprintf('veg %d/%d, mt%d/%d', vegi, length(veg), mi, length(months)))
    mt = months[mi]
    tvegm = get_month(tveg, mt)
    nx = dim(tvegm[[veg[vegi]]][[1]])[1]
    ny = dim(tvegm[[veg[vegi]]][[1]])[2]
    nt = length(tvegm$yyyymm)
    output_reg$rs[[mi]] = matrix(NA, nx, ny)
    for (vi in 1:nvar){
      output_reg$rw[[vi, mi]] = matrix(NA, nx, ny)
      output_reg$rws[[vi, mi]] = matrix(NA, nx, ny)
      output_reg$sn[[vi, mi]] = matrix(NA, nx, ny)
      output_reg$coef_cor[[vi, mi]] = matrix(NA, nx, ny)
      output_reg$pval_cor[[vi, mi]] = matrix(NA, nx, ny)
      output_reg$coef_partialcor[[vi, mi]] = matrix(NA, nx, ny)
      output_reg$pval_partialcor[[vi, mi]] = matrix(NA, nx, ny)
      output_reg$coef_standreg[[vi, mi]] = matrix(NA, nx, ny)
      output_reg$pval_standreg[[vi, mi]] = matrix(NA, nx, ny)
    }
    for (xi in 1:nx){
      if (xi %% 10 == 0)
      {        print(sprintf('%d/%d', xi, nx))}
      for (yi in 1:ny){
        te = matrix(NA, nt , nvar + 1)
        te[,1] = arrayfun(function(x)x[xi, yi], tvegm[[veg[vegi]]])
        for (vi in 1:nvar){
          te[,vi+1] = arrayfun(function(x)x[xi, yi], tvegm[[var[vi]]])
        }
        te = as.data.frame(te)
        names(te) = c(veg[vegi], var)

        idxall =colMeans(is.na(t(te))) ==0
        if ((sum(idxall) < 10) || (sd(te[,1], na.rm = T)==0)){

        }
        else{
          te = te[idxall,]
          if (1 %in% option){
            # correlation
            for (vi in 1:nvar){
              tc = cor.test(te[,1], te[,vi+1])
              output_reg$coef_cor[[vi, mi]][xi,yi] = as.numeric(tc$estimate)
              output_reg$pval_cor[[vi, mi]][xi,yi] = tc$p.value
            }
          }
          if (2 %in% option){
            # partial correlation
            for (vi in 1:nvar){
              ttt <- try({tc = pcor.test(te[,1], te[,vi+1], te[,setdiff(1:nvar, vi)+1],'pearson')
              output_reg$coef_partialcor[[vi, mi]][xi,yi] = tc$estimate
              output_reg$pval_partialcor[[vi, mi]][xi,yi] = tc$p.value
              }, silent = T)
            }
          }
          # standardize
          
          t0 = as.data.frame(scale(te))
          
          if (3 %in% option){
            # relative weight
            tlm = rwa(t0, veg[vegi], var)
            output_reg$rs[[mi]][xi, yi] = tlm$rsquare
            for (vi in 1:nvar){
              output_reg$rw[[vi, mi]][xi,yi] = tlm$result$Raw.RelWeight[vi]
              output_reg$rws[[vi, mi]][xi,yi] = tlm$result$Rescaled.RelWeight[vi]
              output_reg$sn[[vi, mi]][xi,yi] = tlm$result$Sign[vi]
            }
          }
          
          if (4 %in% option){
            # standardized regression weight
            tlm = lm(tstr, t0)
            tlm = summary(tlm)
            for (vi in 1:nvar){
              output_reg$coef_standreg[[vi, mi]][xi,yi] = as.numeric(tlm$coefficients[vi+1,1])
              output_reg$pval_standreg[[vi, mi]][xi,yi] = as.numeric(tlm$coefficients[vi+1,4])
            }
          }
        }
      }
    }
  }
  output_reg$vegname = veg[vegi]
  save(output_reg, file = tfilename)
}



