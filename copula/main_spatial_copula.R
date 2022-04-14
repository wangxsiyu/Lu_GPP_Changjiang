library(WangTools)
library(LuGPP)
library(pracma)
clc()
source('./function_copula_MT.R')
source('./function_spatial.R')
library(VineCopula)
library(CDVineCopulaConditional)
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
####################### compute marginal
getall <- function(lst){
    out = c()
    for (i in 1:length(lst)){
      tout = lst[[i]]
      tout = tout[!is.na(tout)]
      out = c(out, tout)
    }
    return(out)
}
# veg = getall(map)
# h = getall(xp[1,])
# d = getall(xp[2,])
# vmar = MT_copula_marginal(veg)
# hmar = MT_copula_marginal(h)
# dmar = MT_copula_marginal(d, c(1,5,6))
# idv = vmar$idbest
# print(vmar$distnames[idv])
# idh = hmar$idbest
# print(hmar$distnames[idh])
# idd = dmar$idbest
# print(dmar$distnames[idd])
# 
# # plot marginal Q-Q plot
# {
# WD = 1920
# HT = 700
# jpeg("./figs/qqplot.jpg", width = WD, height = HT)
# W_figure(1,3)
# # plot(vmar$fits[[idv]])
# par(cex = 1.5)
# par(mar =  c(8.1, 7.1, 7.1, 2.1))
# qqplot(ecdf(veg)(veg), vmar$x, xlab = "Empirical quantile", ylab = "Theoretical quantile", main = "GPP",
#        cex.axis = 1, cex.lab = 1.5, cex.main = 2)
# # dev.off()
# # jpeg("./figs/qqplot_h.jpg", width = WD, height = HT)
# par(mar =  c(8.1, 7.1, 7.1, 2.1))
# qqplot(ecdf(h)(h), hmar$x, xlab = "Empirical quantile", ylab = "Theoretical quantile", main = "Tmax",
#        cex.axis = 1, cex.lab = 1.5, cex.main = 2)
# # plot(hmar$fits[[idh]])
# # dev.off()
# # jpeg("./figs/qqplot_d.jpg", width = WD, height = HT)
# par(mar =  c(8.1, 7.1, 7.1, 2.1))
# qqplot(ecdf(d)(d), dmar$x, xlab = "Empirical quantile", ylab = "Theoretical quantile", main = "P-PET",
#        cex.axis = 1, cex.lab = 1.5, cex.main = 2)
# # plot(dmar$fits[[idd]])
# dev.off()
# }

#### get marginals
# vmar_data = convert2marginal(map[1,], vmar$pars[[vmar$idbest]], vmar$distnames[vmar$idbest])
# hmar_data = convert2marginal(xp[1,], hmar$pars[[hmar$idbest]], hmar$distnames[hmar$idbest])
# dmar_data = convert2marginal(xp[2,], dmar$pars[[dmar$idbest]], dmar$distnames[dmar$idbest])
#### compute copula by vegmap
# savename = './copula_spatial_0413.RData'
# compute_spatial(vmar_data, hmar_data, dmar_data, savename)
# load(savename)
savename = './copula_spatial_xyz.RData'
compute_xyz_spatial(map[1,], xp, mts, vmar$distnames[vmar$idbest], hmar$distnames[hmar$idbest], dmar$distnames[dmar$idbest], savename)
load(savename)
savename = './copula_spatial_xyz_emp.RData'
compute_xyz_spatial_emp(map[1,], xp, mts, savename)
load(savename)

{
  WD = 1920
  HT = 700
  jpeg("./figs/qqplot.jpg", width = WD, height = HT)
  W_figure(1,3)
  # plot(vmar$fits[[idv]])
  par(cex = 1.5)
  par(mar =  c(8.1, 7.1, 7.1, 2.1))
  veg = h = d = vmarx = hmarx = dmarx = c()
  for (vegi in 1:6){
    for (mi in 1:5){
      veg = c(veg, d0[[vegi,mi]]$gpp)
      vmarx = c(vmarx, d1[[vegi,mi]]$gpp)
      h = c(h, d0[[vegi,mi]]$H)
      hmarx = c(hmarx, d1[[vegi,mi]]$H)
      d = c(d, d0[[vegi,mi]]$D)
      dmarx = c(dmarx, d1[[vegi,mi]]$D)
    }
  }
  qqplot(vmarx, ecdf(veg)(veg), ylab = "Empirical quantile", xlab = "Theoretical quantile", main = "GPP",
         cex.axis = 1, cex.lab = 1.5, cex.main = 2)
  # dev.off()
  # jpeg("./figs/qqplot_h.jpg", width = WD, height = HT)
  par(mar =  c(8.1, 7.1, 7.1, 2.1))
  qqplot(hmarx, ecdf(h)(h), ylab = "Empirical quantile", xlab = "Theoretical quantile", main = "Tmax",
         cex.axis = 1, cex.lab = 1.5, cex.main = 2)
  # plot(hmar$fits[[idh]])
  # dev.off()
  # jpeg("./figs/qqplot_d.jpg", width = WD, height = HT)
  par(mar =  c(8.1, 7.1, 7.1, 2.1))
  qqplot(dmarx, ecdf(d)(d), ylab = "Empirical quantile", xlab = "Theoretical quantile", main = "P-PET",
         cex.axis = 1, cex.lab = 1.5, cex.main = 2)
  # plot(dmar$fits[[idd]])
  dev.off()
}


### get 2-D
# plt_3d_copula(d0[[1,4]])
library(lattice)
tv = c(arrayfun(function(x)x$distnames[x$idbest], veg1))
ttt = histogram(as.factor(tv))
th = c(arrayfun(function(x)x$distnames[x$idbest], heat1))
plot(as.factor(th))
td = c(arrayfun(function(x)x$distnames[x$idbest], dry1))
plot(as.factor(td))
# fit cvine
savename = './copula_spatial_cvine_emp.RData'
compute_spatial_old(d1, savename)
load(savename)
# simulate
savename = './copula_spatial_simu_emp.RData'
simulate_spatial(d1, cvine, savename)
load(savename)
### plot 2-D (heat vs drought)
library(MASS)
vegi = 1; mi = 1
te2 <- kde2d(d2[[vegi, mi]]$H, d2[[vegi, mi]]$gpp, n = 100)
image(te2)
te1 <- kde2d(d1[[vegi, mi]]$H, d1[[vegi, mi]]$gpp, n = 100)
image(te1)

### sim vs obs
# [vegi, mi]



#### compute MLE
library(akima) 
library(fields) 
# sss0 = seq(0,1,0.05)
# nsss0 = length(sss0)
# msh0 = meshgrid(sss0)
# sss = seq(0,1,0.01)
# nsss = length(sss)
# msh = meshgrid(sss)
# fff<-matrix(0,nsss0*nsss0,4)
# fff[,2] = matrix(msh0$X, ,1)
# fff[,3] = matrix(msh0$Y, ,1)
cpl2 = matrix(list(), 6,nmts)
for (vegi in 1:6){
  for (mi in 1:length(mts)){
    print(sprintf('%d,%d', vegi, mi))
    
    hs = linspace(quantile(heat1[[vegi,mi]]$raw, 0.05), quantile(heat1[[vegi,mi]]$raw, 0.95), 100)
    ds = linspace(quantile(dry1[[vegi,mi]]$raw, 0.05), quantile(dry1[[vegi,mi]]$raw, 0.95), 100)
    h20 = convert2marginal(list(hs), heat1[[vegi,mi]]$pars[[heat1[[vegi,mi]]$idbest]], heat1[[vegi,mi]]$distnames[heat1[[vegi,mi]]$idbest])[[1]]
    d20 = convert2marginal(list(ds), dry1[[vegi,mi]]$pars[[dry1[[vegi,mi]]$idbest]], dry1[[vegi,mi]]$distnames[dry1[[vegi,mi]]$idbest])[[1]]
    h2 = unique(h20)
    d2 = unique(d20)
    hs = hs[arrayfun(function(x)min(which(x == h20)), h2)]
    ds = ds[arrayfun(function(x)min(which(x == d20)), d2)]
    # h2 = seq(0,1,0.01)
    # d2 = seq(0,1,0.01)
    # hs = h2
    # ds = d2
    nh = length(h2)
    nd = length(d2)
    msh = meshgrid(h2, d2)
    fff<-matrix(0,nh * nd,2)
    fff[,1] = matrix(msh$X, ,1)
    fff[,2] = matrix(msh$Y, ,1)
    # f4 = matrix(0,dim(fff)[1],1)
    # for (xi in 1:dim(fff)[1]){
    #   tfff = matrix(0, nsss, 4)
    #   tfff[,2] = fff[xi, 2]
    #   tfff[,3] = fff[xi, 3]
    #   tfff[,1] = sss
    #   te = RVinePDF(tfff[,1:3],cvine[[vegi, mi]])
    #   f4[xi] = mean(sss[max(te) == te])
    # }
    sss0 = linspace(quantile(veg1[[vegi, mi]]$raw, 0), quantile(veg1[[vegi, mi]]$raw, 1), 100)
    sss = convert2marginal(sss0, veg1[[vegi, mi]]$pars[[veg1[[vegi, mi]]$idbest]], veg1[[vegi, mi]]$distnames[veg1[[vegi, mi]]$idbest])
    nsss = length(sss)
    f4 = matrix(0,dim(fff)[1],1)
    for (xi in 1:dim(fff)[1]){
      tfff = matrix(0, nsss, 4)
      tfff[,2] = fff[xi, 1]
      tfff[,3] = fff[xi, 2]
      tfff[,1] = sss
      te = RVinePDF(tfff[,1:3],cvine[[vegi, mi]])
      f4[xi] = mean(sss0[max(te) == te])
    }
    ipt = cbind(fff, f4)
    ipt = as.data.frame(ipt)
    names(ipt) = c("x","y","z")
    tout = list()
    tout$x = hs
    tout$y = ds
    tout$z = acast(ipt, x~y, value.var="z")
    cpl2[[vegi, mi]] = tout
    
    
    
    
    # tid = veg1[[vegi,mi]]$idbest
    # ff4 = getinvmarginal(f4, veg1[[vegi, mi]]$pars[[tid]], veg1[[vegi, mi]]$distnames[tid])
    # tid = heat1[[vegi,mi]]$idbest
    # thhh = getinvmarginal(fff[,2], heat1[[vegi, mi]]$pars[[tid]], heat1[[vegi, mi]]$distnames[tid])
    # tid = dry1[[vegi,mi]]$idbest
    # tddd = getinvmarginal(fff[,3], dry1[[vegi, mi]]$pars[[tid]], dry1[[vegi, mi]]$distnames[tid])
    # tf = data.frame(gpp = ff4, H = thhh, D = tddd)
    # tf = tf[colMeans(t(is.infinite(as.matrix(tf))))==0,]
    # cpl2[[vegi, mi]] = interp(tf$H, tf$D, tf$gpp, nx = 500, ny = 500)
  }
}

# ##temp
# tmp = matrix(list(),6,5)
# for (i in 1:6){
#   for (j in 1:5){
#     tmp[[i,j]] = cpl2[[i,j]]$z
#   }
# }
{
  plt_cp_veg_month(cpl2,sprintf('./figs/cpl_spatial_gpp.jpg'), list(xlm = c(5,35), ylm = c(-150, 300), zlm = c(0,350), col = rev(tim.colors(100))))
}




#### compute ratio
library(fitdistrplus)
library(extRemes)
library(extraDistr)
library(raster)
# hs = seq(quantile(h, 0.01), quantile(h, 0.99), 0.05)
# ds = seq(quantile(d, 0.01), quantile(d, 0.99), 0.5)
# h2 = convert2marginal(list(hs), hmar$pars[[hmar$idbest]], hmar$distnames[hmar$idbest])[[1]]
# d2 = convert2marginal(list(ds), dmar$pars[[dmar$idbest]], dmar$distnames[dmar$idbest])[[1]]
# nh = length(h2)
# nd = length(d2)
# msh = meshgrid(h2, d2)
# fff<-matrix(0,nh * nd,2)
# fff[,1] = matrix(msh$X, ,1)
# fff[,2] = matrix(msh$Y, ,1)
ratio = c(0.05, 0.1,0.3, 0.5)
cpls = matrix(list(),1,length(ratio))
for (ri in 1:length(ratio)){
  tratio = ratio[ri]
  cpls[[ri]]$ratio = tratio
  cpl = matrix(list(), 6,nmts)
  for (vegi in 1:6){
    for (mi in 1:length(mts)){
      print(sprintf('%d, %d,%d', ri, vegi, mi))
      
      hs = linspace(quantile(heat1[[vegi,mi]]$raw, 0), quantile(heat1[[vegi,mi]]$raw, 1), 200)
      ds = linspace(quantile(dry1[[vegi,mi]]$raw, 0), quantile(dry1[[vegi,mi]]$raw, 1), 200)
      h20 = convert2marginal(list(hs), heat1[[vegi,mi]]$pars[[heat1[[vegi,mi]]$idbest]], heat1[[vegi,mi]]$distnames[heat1[[vegi,mi]]$idbest])[[1]]
      d20 = convert2marginal(list(ds), dry1[[vegi,mi]]$pars[[dry1[[vegi,mi]]$idbest]], dry1[[vegi,mi]]$distnames[dry1[[vegi,mi]]$idbest])[[1]]
      h2 = unique(h20)
      d2 = unique(d20)
      hs = hs[arrayfun(function(x)min(which(x == h20)), h2)]
      ds = ds[arrayfun(function(x)min(which(x == d20)), d2)]
      # h2 = seq(0,1,0.01)
      # d2 = seq(0,1,0.01)
      # hs = h2
      # ds = d2
      nh = length(h2)
      nd = length(d2)
      msh = meshgrid(h2, d2)
      fff<-matrix(0,nh * nd,2)
      fff[,1] = matrix(msh$X, ,1)
      fff[,2] = matrix(msh$Y, ,1)
      
      tfff = cbind(tratio, fff)
      
      te = RVinePIT(tfff,cvine[[vegi, mi]])[,1]
      # cpl[[vegi, mi]] = interp(fff[,1], ds, te)
      ipt = cbind(fff, te)
      ipt = as.data.frame(ipt)
      names(ipt) = c("x","y","z")
      library(reshape2)
      tout = list()
      tout$x = hs
      tout$y = ds
      tout$z = acast(ipt, x~y, value.var="z")
      cpl[[vegi, mi]] = tout;
      
      # tout = rasterFromXYZ(ipt)
      # tid = heat1[[vegi,mi]]$idbest
      # thhh = getinvmarginal(fff[,2], heat1[[vegi, mi]]$pars[[tid]], heat1[[vegi, mi]]$distnames[tid])
      # tid = dry1[[vegi,mi]]$idbest
      # tddd = getinvmarginal(fff[,3], dry1[[vegi, mi]]$pars[[tid]], dry1[[vegi, mi]]$distnames[tid])
      # tf = data.frame(gpp = te, H = thhh, D = tddd)
      # tf = tf[colMeans(t(is.infinite(as.matrix(tf))))==0,]
      # cpl[[vegi, mi]] = interp(tf$H, tf$D, tf$gpp, nx = 500, ny = 500)
    }
  }
  cpls[[ri]] = cpl
}
# save(cpls, file = './copula_spatial_plotdata.RData')
# load('./copula_spatial_plotdata.RData')
## plot
library(fields)
for (ri in 1:length(ratio))
{
   plt_cp_veg_month(cpls[[ri]], sprintf('./figs/cpl_spatial_emp_%.2f.jpg', ratio[ri]), ratio[ri],
                    params = list(xlm = NULL, ylm = NULL, zlm = c(0,1),  col= hcl.colors(12, "YlOrRd", rev = TRUE)))
}



