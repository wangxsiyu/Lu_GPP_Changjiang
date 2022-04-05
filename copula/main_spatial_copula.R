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


ratio = c(0.05, 0.1, 0.2, 0.5)
cpls = matrix(list(),1,length(ratio))
for (ri in 1:length(ratio)){
  tratio = ratio[ri]
  cpls[[ri]]$ratio = tratio
  ### compute conditional prob
  library(akima) 
  library(fields) 
  sss = seq(0,1,0.01)
  nsss = length(sss)
  msh = meshgrid(sss)
  fff<-matrix(0,nsss*nsss,4)
  fff[,2] = matrix(msh$X, ,1)
  fff[,3] = matrix(msh$Y, ,1)
  fff[,1] = tratio
  cpl = matrix(list(), 6,nmts)
  for (vegi in 1:6){
    for (mi in 1:length(mts)){
      print(sprintf('%d, %d,%d', ri, vegi, mi))
      te = RVinePIT(fff[,1:3],cvine[[vegi, mi]])[,1]
      tid = heat1[[vegi,mi]]$idbest
      thhh = getinvmarginal(fff[,2], heat1[[vegi, mi]]$pars[[tid]], heat1[[vegi, mi]]$distnames[tid])
      tid = dry1[[vegi,mi]]$idbest
      tddd = getinvmarginal(fff[,3], dry1[[vegi, mi]]$pars[[tid]], dry1[[vegi, mi]]$distnames[tid])
      tf = data.frame(gpp = te, H = thhh, D = tddd)
      tf = tf[colMeans(t(is.infinite(as.matrix(tf))))==0,]
      cpl[[vegi, mi]] = interp(tf$H, tf$D, tf$gpp, nx = 500, ny = 500)
    }
  }
  cpls[[ri]] = cpl
}
save(cpls, file = './copula_spatial_plotdata.RData')
## plot


ratio = c(0.05, 0.1, 0.2, 0.5)
for (ri in 1:length(ratio)){
  nms_veg = c("forest","open forest","shrub land","paddy field","dry land","grassland")
  png(filename = sprintf('./cpl_spatial_%.2f.jpg', ratio[ri]),width = 1920, height = 1080, units = "px",
    bg = "transparent",  res = NA)
library(fields)
# library(RColorBrewer)
set.panel()
ind <- split.screen(c(6,nmts))
for(ii in 1:length(ind)){
  vegi = ceiling(ii /nmts)
  mi = mod0(ii, nmts)
  print(sprintf('%d,%d', vegi, mi))
  if (vegi == 1){
    tmain = sprintf("month = %d", mts[mi])
  } else {
    tmain = "";
  }
  if (mi == 1){
    tylb = nms_veg[vegi]
  } else {
    tylb = "p-pet"
  }
  Fdata <- cpls[[ri]][[vegi, mi]]
  screen(ind[ii])
  image.plot(Fdata$x,Fdata$y,Fdata$z,main=tmain,cex.main=3,xlab="tmax",ylab= tylb,
             axes=F,col =tim.colors(20),axis.args = list(cex.axis=2),
             legend.args=list( text=""),smallplot=c(0.85,0.9,0.1,0.9)
  )
}
close.screen(all=TRUE);
dev.off();
}



################### 
library(akima) 
library(fields) 
sss0 = seq(0,1,0.05)
nsss0 = length(sss0)
msh0 = meshgrid(sss0)
sss = seq(0,1,0.01)
nsss = length(sss)
msh = meshgrid(sss)
fff<-matrix(0,nsss0*nsss0,4)
fff[,2] = matrix(msh0$X, ,1)
fff[,3] = matrix(msh0$Y, ,1)
cpl2 = matrix(list(), 6,nmts)
for (vegi in 1:6){
  for (mi in 1:length(mts)){
    print(sprintf('%d,%d', vegi, mi))
    f4 = matrix(0,dim(fff)[1],1)
    for (xi in 1:dim(fff)[1]){
      tfff = matrix(0, nsss, 4)
      tfff[,2] = fff[xi, 2]
      tfff[,3] = fff[xi, 3]
      tfff[,1] = sss
      te = RVinePDF(tfff[,1:3],cvine[[vegi, mi]])
      f4[xi] = mean(sss[max(te) == te])
    }
    tid = veg1[[vegi,mi]]$idbest
    ff4 = getinvmarginal(f4, veg1[[vegi, mi]]$pars[[tid]], veg1[[vegi, mi]]$distnames[tid])
    tid = heat1[[vegi,mi]]$idbest
    thhh = getinvmarginal(fff[,2], heat1[[vegi, mi]]$pars[[tid]], heat1[[vegi, mi]]$distnames[tid])
    tid = dry1[[vegi,mi]]$idbest
    tddd = getinvmarginal(fff[,3], dry1[[vegi, mi]]$pars[[tid]], dry1[[vegi, mi]]$distnames[tid])
    tf = data.frame(gpp = ff4, H = thhh, D = tddd)
    tf = tf[colMeans(t(is.infinite(as.matrix(tf))))==0,]
    cpl2[[vegi, mi]] = interp(tf$H, tf$D, tf$gpp, nx = 500, ny = 500)
  }
}

{
nms_veg = c("forest","open forest","shrub land","paddy field","dry land","grassland")
png(filename = sprintf('./cpl_spatial_gpp.jpg'),width = 1920, height = 1080, units = "px",
    bg = "transparent",  res = NA)
library(fields)
# library(RColorBrewer)
set.panel()
ind <- split.screen(c(6,nmts))
for(ii in 1:length(ind)){
  vegi = ceiling(ii /nmts)
  mi = mod0(ii, nmts)
  print(sprintf('%d,%d', vegi, mi))
  if (vegi == 1){
    tmain = sprintf("month = %d", mts[mi])
  } else {
    tmain = "";
  }
  if (mi == 1){
    tylb = nms_veg[vegi]
  } else {
    tylb = "p-pet"
  }
  Fdata <- cpl2[[vegi, mi]]
  screen(ind[ii])
  image.plot(Fdata$x,Fdata$y,Fdata$z,main=tmain,cex.main=3,xlab="tmax",ylab= tylb,
             axes=T,col =rev(tim.colors(20)),axis.args = list(cex.axis=2),
             legend.args=list( text=""),smallplot=c(0.85,0.9,0.1,0.9)
  )
}
close.screen(all=TRUE);
dev.off();
}
