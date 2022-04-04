library(WangTools)
library(LuGPP)
library(pracma)
library(rwa)
clc()
datadir = './data/gpp/'
data = loadraw(datadir)
# get spatial average

vars = loadraw('./data/',c("tmax","ppet"))
mts = 5:9
nms = names(data)
map = matrix(list(), length(data), length(mts))
for (di in 1:length(data)){
  for (mi in 1:length(mts)){
    td = get_month(data[[nms[di]]], mts[mi])
    map[[di, mi]] = list_mean(td$data)
  }
}

xp = matrix(list(), 2, length(mts))
for (mi in 1:length(mts)){
  td = get_month(vars$tmax, mts[mi])
  xp[[1, mi]] = list_mean(td$data)
  td = get_month(vars$ppet, mts[mi])
  xp[[2, mi]] = list_mean(td$data)
}

gppi = which(nms == "SDgpp")

## sensitivity analysis
out = matrix(list(), 6, length(mts))
for (mi in 1:length(mts)){
  for (vegi in 1:6){
    print(sprintf("%d,%d", mi, vegi))
      tout = list()
      idx = vegmap == vegi
      te_gpp = map[[gppi, mi]][idx]
      th_gpp = xp[[1, mi]][idx]
      td_gpp = xp[[2,mi]][idx]
      tdd = data.frame(gpp = te_gpp, tmax = th_gpp, ppet = td_gpp)
      idx = which(colMeans(t(is.na(tdd))) == 0)
      tdd = tdd[idx,]
      tdd$t0 = scale(tdd$tmax)
      tdd$p0 = scale(tdd$ppet)
      tdd$g0 = scale(tdd$gpp)
      avt = mean(tdd$tmax)
      sdt = sd(tdd$tmax)
      avp = mean(tdd$ppet)
      sdp = sd(tdd$ppet)
      # compute sensitivity, do regression for the closest 100 points near a particular combination of values
      rg_tmax = seq(4,35,0.5) #c(min(tdd$tmax), max(tdd$tmax))
      rg_ppet = seq(-160,360,10)
      nntt = length(rg_tmax)
      nnpp = length(rg_ppet)
      # rw_tmax = rw_ppet =
      tout$coef_tmax = tout$coef_ppet = tout$c0 =  tout$p_tmax = tout$p_ppet = matrix(NA, nntt, nnpp)
      for (xi in 1:(length(rg_tmax)-0)){
        for (yi in 1:(length(rg_ppet)-0)){
           # idx = tdd$tmax >= rg_tmax[xi] & tdd$tmax < rg_tmax[xi+1] &
           #   tdd$ppet >= rg_ppet[yi] & tdd$ppet < rg_ppet[yi+1]
            ttt = (rg_tmax[xi] - avt)/sdt
            ppp = (rg_ppet[yi] - avp)/sdp
            dis = (tdd$t0 - ttt)^2 + (tdd$p0 - ppp)^2
            oodd = order(dis)
            thres = ((10)/sdp)^2 + ((0.5)/sdt)^2
            oodd = oodd[dis[oodd] < thres]
            if (length(oodd) < 30){
               if (length(oodd) > 10){
                  temp = tdd[oodd,]
               } else{
                  temp = NULL
               }
            }else{
               oodd = oodd[1:30]
               temp = tdd[oodd,]
            }
            if (!is.null(temp)){
              # tlm = rwa(temp, 'g0', c('t0','p0'))
              # output_reg$rs[[mi]][xi, yi] = tlm$rsquare
              # for (vi in 1:nvar){
              #   output_reg$rw[[vi, mi]][xi,yi] = tlm$result$Raw.RelWeight[vi]
              #   output_reg$rws[[vi, mi]][xi,yi] = tlm$result$Rescaled.RelWeight[vi]
              #   output_reg$sn[[vi, mi]][xi,yi] = tlm$result$Sign[vi]
              # }
              tlm = lm('g0~t0 + p0', temp)
              tlm = summary(tlm)
                tout$coef_tmax[xi,yi] = as.numeric(tlm$coefficients[2,1])
                tout$p_tmax[xi,yi] = as.numeric(tlm$coefficients[2,4])
                tout$coef_ppet[xi,yi] = as.numeric(tlm$coefficients[3,1])
                tout$p_ppet[xi,yi] = as.numeric(tlm$coefficients[3,4])
                tout$c0[xi,yi] = as.numeric(tlm$coefficients[1,1])
            }
        }
      }
      out[[vegi, mi]] = tout
  }
}
