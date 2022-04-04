library(WangTools)
library(LuGPP)
library(pracma)
clc()
datadir = './data/gpp/'
veg = list.files(datadir)
veg = tools::file_path_sans_ext(veg)
veg = setdiff(veg, c("AVgpp", "SDgpp"))
veg = setdiff(veg, c("NIRvgpp","RSgpp","ECLUEgpp","MODISgpp"))
data = loadraw(datadir, veg)
data = data[names(data)%in%veg]
ng = length(data)
## get the data for the variable first
yr = c()
for (i in 1:ng){
  yr = c(yr, data[[i]]$yyyymm)
}
yr = sort(unique(yr))

gpp0 = list()
gpp0$yyyymm = yr
gpp0$mt = yr %% 100
gpp0$yr = floor(yr/100)
gpp0$data = matrix(list(), 1, length(yr))
nyr = length(gpp0$data)
sz = dim(data[[1]]$data[[1]])
for (i in 1:length(gpp0$data)){
  gpp0$data[[i]] = matrix(NA, sz[1],sz[2])
}

# AV/SD
av = sd = matrix(list(),1,ng)
for (i in 1:ng){
  av[[i]] = list_mean(data[[i]]$data)
  sd[[i]] = list_sd(data[[i]]$data)
}
avall = list_mean(av)
sdall = list_mean(sd)

# d1 = d2 = data
# for (i in 1:ng){
#   for (j in 1:length(data[[i]]$data)){
#      d1[[i]]$data[[j]] = (d1[[i]]$data[[j]] - av[[i]])
#   }
# }

gpp1 = gpp2 = gpp0
for (xi in 1:sz[1]){
  print(sprintf("%d", xi))
  for (yi in 1:sz[2]){
     tt1 = tt2 = matrix(NA, ng, nyr)
     for (gi in 1:ng){
       x = data[[gi]]
       idx = yr %in% x$yyyymm
       tt1[gi, idx] = arrayfun(function(t)t[xi,yi], x$data) - av[[gi]][xi,yi]
       tt2[gi, idx] = tt1[gi, idx]/sd[[gi]][xi,yi]
     }
     tt1 = colMeans(tt1, na.rm = T)
     tt2 = colMeans(tt2, na.rm = T)
     for (i in 1:nyr){
        gpp1$data[[i]][xi,yi] = tt1[i] + avall[xi,yi]
        gpp2$data[[i]][xi,yi] = tt2[i]*sdall[xi,yi] + avall[xi,yi]
     }
  }
}

gpp = gpp1
save(gpp, file = './data/gpp/AVgpp.RData')
gpp = gpp2
save(gpp, file = './data/gpp/SDgpp.RData')

################################### PCA
library(FactoMineR)
library(factoextra)

dd = data[1:4]
ng2 = 4
dd = get_common(dd)
gpp = gpp0
pcev = matrix(NA, sz[1], sz[2])
for (xi in 1:sz[1]){
  print(sprintf("%d", xi))
  for (yi in 1:sz[2]){
    tt = matrix(NA, ng2, nyr)
    for (gi in 1:ng2){
      x = dd[[gi]]
      idx = yr %in% x$yyyymm
      tt[gi, idx] = arrayfun(function(t)t[xi,yi], x$data)
    }
    tt = t(tt)
    tt = tt[colMeans(is.na(t(tt))) == 0,]
    if (dim(tt)[1] > 10){
      tsd = arrayfun(function(t)sd(tt[,t]), 1:ng2)
      te = prcomp(tt[, tsd > 0], scale = T)
      tte = summary(te)
      pcev[xi,yi] = tte$importance[2,1]
      pc1 = te$x[,1]
      pc1 = c(scale(pc1))

      pc1 = matrix(pc1,ncol = 1) %*% matrix(tsd, nrow = 1)
      pc1 = c(colMeans(t(pc1), na.rm = T))
      for (i in 1:length(pc1)){
        gpp$data[[which(idx)[i]]][xi,yi] = pc1[i] + mean(te$center)
      }
    }
  }
}
save(gpp, file = './data/gpp/PCA4gpp.RData')


plt_raster(pcev, list(legend = T, breaks = seq(0.45,0.95,0.02)))
save(pcev, file = './output/pca4var.RData')
