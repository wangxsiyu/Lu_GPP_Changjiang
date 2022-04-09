library(WangTools)
library(LuGPP)
library(pracma)
clc()
datadir = '../data/yearly/'
datadir = '../data/yearly2/'
data = loadraw(datadir)
idx = data$spei$yr %in% data$AVSDgpp$yr
data$spei$data = data$spei$data[idx]
data$spei$yr = data$spei$yr[idx]
####################### compute copula
source('./function_copula_MT.R')
library(VineCopula)
library(CDVineCopulaConditional)
nx = dim(data$AVSDgpp$data[[1]])[1]
ny = dim(data$AVSDgpp$data[[1]])[2]


cvine = d1 = veg1 = heat1 = dry1 = matrix(list(), nx, ny)

for (xi in 1:nx){
  print(sprintf('%d/%d', xi, nx))
  for (yi in 1:ny){
    gpp = arrayfun(function(x)x[xi,yi], data$AVSDgpp$data)
    h = arrayfun(function(x)x[xi,yi], data$heatindex$data)
    d = arrayfun(function(x)x[xi,yi], data$spei$data)
    td = data.frame(gpp = detrend(c(gpp)), H = detrend(c(h)), D = detrend(c(d)))
    td = td[colMeans(t(is.na(td))) == 0,]
    if (dim(td)[1] >= 30){
      ## compute marginal
      veg1[[xi, yi]] = MT_copula_marginal(td$gpp)
      heat1[[xi, yi]] = MT_copula_marginal(td$H)
      dry1[[xi, yi]] = MT_copula_marginal(td$D, c(1,5,6))
      d1[[xi, yi]] = data.frame(gpp = veg1[[xi, yi]]$x, H = heat1[[xi, yi]]$x, D = dry1[[xi, yi]]$x)
      ## fit cvine
      ff = as.matrix(d1[[xi, yi]])
      invisible(capture.output(cvine[[xi, yi]]<-CDVineCondFit(ff,Nx=2,c(1,2,3,4,5,9),rotation=F,treecrit="AIC",type="CVine",selectioncrit="AIC")))
    } else{
      
    }
  }
}
save(d1, dry1, veg1, heat1, cvine, file = './copula_temporal_heatspei_detrend.RData')


library(akima) 
library(fields) 
ratio = c(0.1, 0.2, 0.3, 0.5)
hs = c(0.95, 0.9, 0.8, 0.5)
ds = c(0.05, 0.1, 0.2, 0.5)
cpls = matrix(list(),1,length(ratio))
for (ri in 1:length(ratio)){
  cpls[[ri]]$prob = matrix(list(), length(hs), length(ds))
  tratio = ratio[ri]
  cpls[[ri]]$ratio = tratio
  ### compute conditional prob
  for (hi in 1:length(hs)){
    for (di in 1:length(ds)){
      cpls[[ri]]$prob[[hi, di]] = matrix(NA, nx, ny)
      Hval = hs[hi]
      Dval = ds[di]
      for (xi in 1:nx){
        print(sprintf('%d,%d,%d, %d/%d', ri, hi, di, xi, nx))
        for (yi in 1:ny){
          if (!is.null(cvine[[xi,yi]])){
            fff = matrix(NA, 1,3)
            fff[1] = tratio
            fff[2] = Hval
            fff[3] = Dval
            cpls[[ri]]$prob[[hi, di]][xi,yi]= RVinePIT(fff,cvine[[xi, yi]])[1]
          }
        }
      }
    }
  }
}
save(cpls, file = './copula_temporal_plotdata_heatspei_detrend.RData')
## plot


