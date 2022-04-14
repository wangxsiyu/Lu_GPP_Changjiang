MT_copula_marginal <- function(te, idx_dist = NULL, isjitter = T){
  library(fitdistrplus)
  library(extRemes)
  library(extraDistr)
  distnames = c("norm","weibull","lnorm","gamma","logis","gev")
  if (is.null(idx_dist)){  
    ndist = length(distnames)
    idx_dist = 1:ndist#matrix(1, 1, ndist)
  }
  if (is.numeric(idx_dist)){
    distnames = distnames[(1:length(distnames) %in% idx_dist)]
  } else {
    distnames = idx_dist
  }
  # if (any(te < 0) && "weibull" %in% distnames){
  #   disp("removing weibull, negative values exist")
  #   distnames = setdiff(distnames, "weibull")
  # }
  ndist = length(distnames)
  
  pars = matrix(list(), 1, ndist)
  # set initial parsameter
  for (i in 1:ndist){
    pars[[i]] = getparsinit(distnames[i])
  }
  fits = matrix(list(), 1, ndist)
  ## fit model
  for (i in 1:ndist){
    fits[[i]] = getfitdist(te, distnames[i])
    pars[[i]] = getparsfit(fits[[i]], distnames[i])
  }
    
  # extreme value distribution -- http://tougao.ecoagri.ac.cn/html/zgstny/2019/12/2019-1206.htm#outline_anchor_8
  # plot(fit)
  
  ## Goodness-of-fit statistics
  # ks test
  ks<-matrix(0,1,ndist)
  for (i in 1:ndist){
    if (isjitter){
       tte = jitter(te)
    } else{
       tte = te
    }
    tname = paste("p",distnames[i], sep = "")
    if (tname %in% "pgev"){
      ks[i] = ks.test(x = tte, y = tname, pars[[i]][1],pars[[i]][2],pars[[i]][3])[["p.value"]]
    } else {
      ks[i] = ks.test(x = tte, y = tname, pars[[i]][1],pars[[i]][2])[["p.value"]]
    }
  }
  ## aic
  aic_re<-matrix(0,1, ndist)
  for (i in 1:ndist) {
    if (distnames[i] == "gev"){
      aic_re[i] = summary(fits[[i]], silent = T)$AIC
    } else {
      aic_re[i] = gofstat(fits[[i]], fitnames = c(distnames[i]))$aic
    }
  }
  
  ## marginal count
  marginaldist<- which.min(aic_re)
  
  ##  2. calculte marginal distribution
  ## "norm","weibull","lognorm","gamma","logistic","gev"
  Xmarginal = getmarginal(te, pars[[marginaldist]], distnames[marginaldist])
  out = list()
  out$x = Xmarginal
  out$raw = te
  out$fits = fits
  out$pars = pars
  out$idbest = marginaldist
  out$distnames = distnames
  return(out)
}

getparsinit<- function(name){
  if (name %in% c("gev")){
    out = matrix(0, 3, 1)
  } else {
    out = matrix(0,2,1)
  }
  return(out)
}
getparsfit <- function(fits, name){
  if (name == "gev"){
    pars = fits$results$par
  } else {
    pars = fits$estimate
  }
  return(pars)
}
getfitdist<- function(te, name){
  if (name == "gev"){
    fits = fevd(te)
  } else if (name == "weibull"){
    fits = fitdist(te, name, method = "mle")
  } else if (name == "gamma"){
    fits = fitdist(te, name, method = "mme")
  } else {
    fits = fitdist(te, name)
  }
  return(fits)
}
getmarginal<- function(te, mpar, name){
  func = get(paste("p", name, sep = ""))
  if (name == "gev"){
    OUT = func(te, mpar[1], mpar[2], mpar[3])
  } else {
    OUT = func(te, mpar[1], mpar[2])
  }
  return(OUT)
}

getinvmarginal<- function(te, mpar, name){
  func = get(paste("q", name, sep = ""))
  if (name == "gev"){
    OUT = func(te, mpar[1], mpar[2], mpar[3])
  } else {
    OUT = func(te, mpar[1], mpar[2])
  }
  return(OUT)
}
library(fields)#tim.colors(100)
plt_cp_veg_month <- function(xxx, fname, ratio, params = list(xlm = NULL, ylm = NULL, zlm = c(0,1), col= hcl.colors(12, "YlOrRd", rev = TRUE))){
  nms_veg = c("Forest","Open Forest","Shrub Land","Paddy Field","Dry Land","Grassland")
  png(filename = fname,width = 1280, height = 1200, units = "px",
      bg = "transparent",  res = NA)
  par(bg = "#ffffff")
  mat <-t(matrix(1:30,5,6)) 
  mat = cbind(mat, c(31,31,0,0,0,0))
  nf <- W_figure(mat, marg = c(0.15, 0.22, 0.15, 0), w = c(1,1,1,1,1,0.3), h = c(1,1,1,1,1,1))
  for(ii in 1:30){
    par(mar = c(3.5 ,4.5, 2, 0.5), xpd = NA)
    vegi = ceiling(ii /nmts)
    mi = mod0(ii, nmts)
    print(sprintf('%d,%d', vegi, mi))
    if (vegi == 1){
      tmain = sprintf("month = %d\n", mts[mi])
    } else {
      tmain = "";
    }
    if (mi == 1){
      tylb = paste(nms_veg[vegi], "p-pet",sep ='\n')
    } else {
      tylb = "p-pet"
    }
    Fdata <- xxx[[vegi, mi]]
    if (!is.null(params$xlm)){
       xlm = params$xlm
    } else {
      xlm = c(min(Fdata$x, na.rm = T), max(Fdata$x, na.rm = T))
    }
    if (!is.null(params$ylm)){
      ylm = params$ylm
    } else {
      ylm = c(min(Fdata$y, na.rm = T), max(Fdata$y, na.rm = T))
    }
    # plot(1:10)
    pppp <- par('usr')
    do.call('clip', as.list(pppp))
    image(Fdata$x,Fdata$y,Fdata$z,zlim = params$zlm,col =params$col,
          main=tmain,xlab="tmax", ylab= tylb, axes=T,
          cex.main = 3, cex.lab = 3, cex.axis = 2, useRaster = FALSE, xlim = xlm, ylim = ylm)
    
    # legend.args=list( text=""),smallplot=c(0.85,0.9,0.1,0.9)
  }
  par(mar = c(2 ,1, 2, 5), pty = "m", err = -1,cex.axis = 2)
  breaks = linspace(0,1,length(params$col)+1) * max(params$zlm)
  ix <- 1:2
  iy <- breaks
  nBreaks <- length(breaks)
  midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
  iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
  image(ix,iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "",col =params$col,breaks = breaks, las = 2, cex.lab = 3)
  axis.args <- c(list(side = 4, mgp = c(3, 1, 0), 
                      las = 2, 
                      at = seq(0,1,0.1)* max(params$zlm), cex.axis = 2))
  do.call(axis,axis.args)
  text(2,par('usr')[3] - 0.3, labels = sprintf("P(GPP < %.f%%)", ratio * 100),
       srt = 270, cex = 3)
  # axis.args <- c(list(side = 1, mgp = c(3, 1, 0), 
  #                     las = 2, cex.lab = 3))
  # do.call(axis,axis.args)
  dev.off();
}


plt_cpT_veg_month <- function(xxx, fname, hs, ds,set = 1,params = list(zlm = c(0,1), col= hcl.colors(12, "YlOrRd", rev = TRUE))){
  library(pracma)
  library(fields)
  nms_veg = c("forest","open forest","shrub land","paddy field","dry land","grassland")
  png(filename = fname,width = 1920, height = 1024, units = "px",
      bg = "transparent",  res = NA)
  par(bg = "#ffffff")
  mat <-t(matrix(1:16,4,4)) 
  mat = cbind(mat, c(17,0,0,0))
  nf = W_figure(mat, marg = c(0.15, 0.1, 0.1, 0), w = c(1,1,1,1,0.2))
  # set.panel()
  # ind <- split.screen(c(6,nmts))
  breaks = linspace(0,1,length(params$col)+1) * max(params$zlm)
  for(ii in 1:16){
    par(mar = c(2.5 ,3.5, 1.5, 0.5), mgp = c(2,0.8,0), xpd = NA)
    hi = ceiling(ii /4)
    di = mod0(ii, 4)
    print(sprintf('%d,%d', hi, di))
    if (hi == 1){
      if (set == 1)
       tmain = sprintf("P-PET = %.2f",ds[di])
      else if (set == 2)
        tmain = sprintf("SPEI = %.2f",ds[di])
      else if (set == 3 || set == 4 || set == 5)
        tmain = sprintf("p(GPP < %.2f)",hs[di])
    } else {
      tmain = "";
    }
    if (di==1){
      if (set == 1)
        tylb = sprintf("Tmax = %.2f", hs[hi])
      else if (set == 2)
        tylb = sprintf("HEAT = %.2f", hs[hi])
      else if (set == 3)     
        tylb = sprintf("SPEI = %.2f\nHEAT = 0.5",ds[hi])
      else if (set == 4)     
        tylb = sprintf("HEAT = %.2f\nSPEI = 0.5",ds[hi])
      else if (set == 5)     
        tylb = sprintf("HEAT = %.2f\nSPEI = %.2f",ds[hi,1],ds[hi,2])

        
      
    } else {
      tylb = ""
    }
    Fdata <- xxx[[hi, di]]
    # plot(1:10)
    # tmap<- as.matrix(Fdata) #t(map);
    # lon <- c(90.55001, 121.85)
    # lat <- c(24.55, 35.75)
    # Fmap <-raster(tmap, #[nrow(tmap):1,],
    #               xmn=min(lon), xmx=max(lon),
    #               ymn=min(lat), ymx=max(lat),
    #               crs=CRS("+proj=longlat +datum=WGS84"))
    # plot(Fmap, col = params$col, breaks = breaks, main = "",xlab ="", ylab = "",axes = F)
    image(t(Fdata[dim(Fdata)[1]:1,]),zlim = params$zlm,col =params$col,
          main=tmain, cex.main=3,xlab="", ylab= tylb, axes=F,
          cex.main = 3, cex.lab = 3, cex.axis = 2, useRaster = TRUE)
    # legend.args=list( text=""),smallplot=c(0.85,0.9,0.1,0.9)
  }
  par(mar = c(3.5 ,1.5, 1.5, 3.5), pty = "m", err = -1, xpd = T)
  ix <- 1:2
  iy <- breaks
  nBreaks <- length(breaks)
  midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
  iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
  image(ix,iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "",col =params$col,breaks = breaks)
  axis.args <- c(list(side = 4, mgp = c(3, 1, 0), 
                      las = 2, 
                      at = seq(0,1,0.1)* max(params$zlm), cex.axis = 2))
  do.call(axis,axis.args)
  dev.off();
}
