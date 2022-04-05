library(WangTools)
library(LuGPP)
library(pracma)
clc()
datadir = '../data/gpp/'
data = loadraw(datadir,"SDgpp")
vars = loadraw('../data/',c("tmax","ppet"))
# get spatial average
mts = 5:9
nmts = length(mts)
nms = names(data)
map = matrix(list(), length(data),nmts)
for (di in 1:length(data)){
  for (mi in 1:length(mts)){
    td = get_month(data[[nms[di]]], mts[mi])
    map[[di, mi]] = list_mean(td$data)
  }
}
xp = matrix(list(), 2,nmts)
for (mi in 1:length(mts)){
  td = get_month(vars$tmax, mts[mi])
  xp[[1, mi]] = list_mean(td$data)
  td = get_month(vars$ppet, mts[mi])
  xp[[2, mi]] = list_mean(td$data)
}
####################### compute copula
source('./function_copula_MT.R')
library(VineCopula)
library(CDVineCopulaConditional)
cvine = d0 = d1 = veg1 = heat1 = dry1 = matrix(list(), 6,nmts)
for (vegi in 1:6){
  for (mi in 1:length(mts)){
    print(sprintf('%d,%d', vegi, mi))
    idx = vegmap == vegi
    tveg = map[[1,mi]][idx]
    theat = xp[[1, mi]][idx]
    tdry = xp[[2, mi]][idx]
    td = data.frame(gpp = tveg, H = theat, D = tdry)
    td = td[colMeans(t(is.na(td))) == 0,]
    d0[[vegi, mi]] = td
    ## compute marginal
    veg1[[vegi, mi]] = MT_copula_marginal(td$gpp,1)
    heat1[[vegi, mi]] = MT_copula_marginal(td$H,1)
    dry1[[vegi, mi]] = MT_copula_marginal(td$D, 1)# c(1,3,4,5,6))
    d1[[vegi, mi]] = data.frame(gpp = veg1[[vegi, mi]]$x, H = heat1[[vegi, mi]]$x, D = dry1[[vegi, mi]]$x)
    ## fit cvine
    ff = as.matrix(d1[[vegi, mi]])
    cvine[[vegi, mi]]<-CDVineCondFit(ff,Nx=2,c(1,2,3,4,5,9),rotation=F,treecrit="AIC",type="CVine",selectioncrit="AIC")
  }
}
save(d0, d1, dry1, veg1, heat1, cvine, file = './copula_spatial.RData')