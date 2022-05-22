{
library(WangTools)
library(LuGPP)
library(pracma)
library(ppcor)
library(rwa)
filehead = 'reglagNEW_'
clc()
datadir = './data/gpp/'
data = loadraw(datadir)
vars = loadraw('./data/',c('tmax','ppet'))
## veg vs var
veg = names(data)
var = names(vars)
nvar = length(var)
## get the data for the variable first


hd = list()
nlag = 5
for (i in 0:nlag){
  hd[[sprintf('h_%d', i)]] = get_lag(vars$tmax, i)
  hd[[sprintf('d_%d', i)]] = get_lag(vars$ppet, i)
}
hnms = arrayfun(function(x){sprintf('h_%d',x)}, 0:nlag)
dnms = arrayfun(function(x){sprintf('d_%d',x)}, 0:nlag)
hdnms = c(hnms, dnms)
strhd = paste(hdnms, collapse = "+")

hd = get_common(hd)
months = 1:12
}

veg = c("Pgpp","SDgpp","SIFgpp","VPMgpp","Musyqgpp","LUEgpp")
for (vegi in 1:length(veg)) {
  tveg = hd
  tveg$veg = data[[veg[vegi]]]
  tveg = merge_vars(tveg)
  tstr = paste("veg", "~", strhd)
  tfile = sprintf("%slag%d_%s.RData", filehead, nlag, veg[vegi])
  tfilename = file.path("./output/regression/",tfile)
  if (file.exists(tfilename)){
    next;
  }
  output_reg = list()
  output_reg$rs = output_reg$ps = output_reg$parcor_predH = output_reg$parcor_predD = matrix(list(),1,12)
  output_reg$coef_standreg = output_reg$pval_standreg = matrix(list(), nlag*2+2, 12)
  for (mi in 1:length(months)){
    print(sprintf('veg %d/%d, mt%d/%d', vegi, length(veg), mi, length(months)))
    mt = months[mi]
    tvegm = get_month(tveg, mt)
    nx = dim(tvegm$veg[[1]])[1]
    ny = dim(tvegm$veg[[1]])[2]
    nt = length(tvegm$yyyymm)
    output_reg$rs[[mi]] = output_reg$ps[[mi]] = matrix(NA, nx, ny)
    output_reg$parcor_predH[[mi]] = matrix(NA, nx, ny)
    output_reg$parcor_predD[[mi]] = matrix(NA, nx, ny)
    for (vi in 1:(2*(nlag+1))){
      output_reg$coef_standreg[[vi, mi]] = matrix(NA, nx, ny)
      output_reg$pval_standreg[[vi, mi]] = matrix(NA, nx, ny)
    }
    for (xi in 1:nx){
      if (xi %% 10 == 0)
      {        print(sprintf('%d/%d', xi, nx))}
      for (yi in 1:ny){
        te = matrix(NA, nt , nlag*2 + 3)
        tnms = names(tveg)
        tnms = setdiff(tnms, c("yyyymm","yr","mt"))
        for (vi in 1:length(tnms)){
          te[,vi] = arrayfun(function(x)x[xi, yi], tvegm[[tnms[vi]]])
        }
        te = as.data.frame(te)
        names(te) = tnms

        idxall =colMeans(is.na(t(te))) ==0
        if ((sum(idxall) <= 10) || (sd(te$veg, na.rm = T)==0)){

        }
        else{
          te = te[idxall,]

          # standardize
          t0 = as.data.frame(scale(te))


          # standardized regression weight
          tlm = lm(tstr, t0)
          tlm = summary(tlm)

          t0nms = names(t0)
          t0nms = t0nms[!t0nms %in% c("veg")]
          cfs = tlm$coefficients[,1]
          cfs = arrayfun(function(x){cfs[names(cfs) %in% x]}, t0nms)

          y0 = t0[,dim(t0)[2]]
          tempx = as.matrix(t0[,1:dim(t0)[2]-1] )
          tempw = matrix(cfs,length(cfs),1)

          ttid = seq(1,12,2)
          predH = tempx[,ttid] %*% tempw[ttid]
          ttid = seq(2,12,2)
          predD = tempx[,ttid] %*% tempw[ttid]

          output_reg$parcor_predH[[mi]][xi,yi] = cor(predH, y0)
          output_reg$parcor_predD[[mi]][xi,yi] = cor(predD, y0)

          for (vi in 1:(dim(tlm$coefficients)[1]-1) ){
            output_reg$coef_standreg[[vi, mi]][xi,yi] = as.numeric(tlm$coefficients[vi+1,1])
            output_reg$pval_standreg[[vi, mi]][xi,yi] = as.numeric(tlm$coefficients[vi+1,4])
          }

          output_reg$rs[[mi]][xi, yi] = tlm$r.squared
          output_reg$ps[[mi]][xi, yi] = pf(tlm$fstatistic[1],tlm$fstatistic[2],tlm$fstatistic[3],lower.tail=F)

        }
      }
    }
  }
  output_reg$vegname = veg[vegi]
  save(output_reg, file = tfilename)
}





######## regression - lag
# outputdir = '../output/reglagnew/'
# library(pracma)
# for (xxi in 1:length(xx)){
#   for (loopj in 1:length(gpp)){
#     # print(sprintf('%d', loopj))
#     pair = c(names(xx)[xxi], names(gpp)[loopj])
#     td0 = list()
#     td0$x0 = getlag(xx[[xxi]], 0)
#     td0$x1 = getlag(xx[[xxi]], 1)
#     td0$x2 = getlag(xx[[xxi]], 2)
#     td0$x3 = getlag(xx[[xxi]], 3)
#     td0$x4 = getlag(xx[[xxi]], 4)
#     td0$x5 = getlag(xx[[xxi]], 5)
#     td0$y = gpp[[loopj]]
# for (mt in 1:12){
#   tic()
#   if (file.exists(file.path(outputdir, sprintf("Reg5_%s_%s_mt%d.RData", pair[1], pair[2], mt)))){
#     print(sprintf('skipping month %d/12', mt));
#     next
#   }
#   print(sprintf('calculating month %d, %d/12', loopj, mt))
#   td = get_months(td0, mt)
#   td = get_common(td)
#   # calculate correlations
#   nx = dim(td[[1]]$data[[1]])[1]
#   ny = dim(td[[1]]$data[[1]])[2]
#   nd = length(td)
#   rs = ps = matrix(list(), 1,nd-1)
#   Rsq = pval = matrix(NA, nx, ny)
#   for (i in 1:(nd-1)){
#     rs[[i]] = ps[[i]] = matrix(NA, nx, ny)
#   }
#   for (xi in 1:nx){
#     if (xi %% 20 == 1){
#       print(sprintf('xy loop: %d/%d', xi, nx))
#     }
#     for (yi in 1:ny){
#       dd = matrix(NA, length(td[[1]]$yr),nd)
#       for (di in 1:nd){
#         tttt = arrayfun(function(x)x[xi, yi], td[[di]]$data)
#         dd[,di] = c(scale(tttt))
#       }
      # dd = as.data.frame(dd)
      # names(dd) = names(td)
      # if (sum(colSums(t(is.na(dd))) == 0) > 10){
        # treg = lm('y~x0+x1+x2+x3+x4+x5', dd)
        # sreg = summary(treg)
        # for (i in 1:(nd-1)){
        #   rs[[i]][xi, yi] = sreg$coefficients[i+1,1]
        #   ps[[i]][xi, yi] = sreg$coefficients[i+1,4]
        # }
        # Rsq[xi, yi] = sreg$r.squared
        # pval[xi, yi] = pf(sreg$fstatistic[1],sreg$fstatistic[2],sreg$fstatistic[3],lower.tail=F)
  #     }
  #   }
  # }
  # out = NULL
  # out$pval = pval
  # out$R2 = Rsq
  # out$ps = ps
  # out$rs = rs
  # save(out, file = file.path(outputdir, sprintf("Reg5_%s_%s_mt%d.RData", pair[1], pair[2], mt)))
  # toc()
# }
# }
# }





