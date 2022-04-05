MT_copula_marginal <- function(te, idx_dist = NULL, isjitter = T){
  library(fitdistrplus)
  library(extRemes)
  library(extraDistr)
  distnames = c("norm","weibull","lnorm","gamma","logis","gev")
  if (is.null(idx_dist)){  
    idx_dist = matrix(1, 1, ndist)
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



plt_cp_veg_month <- function(xxx, fname, params = list(zlm = c(0,1), col= tim.colors(100))){
  png(filename = fname,width = 1080, height = 1080, units = "px",
      bg = "transparent",  res = NA)
  par(bg = "#ffffff")
  mat <-t(matrix(1:30,5,6)) 
  mat = cbind(mat, c(31,0,0,0,0,0))
  # mat = cbind(c(1,6,11,16,21,26), mat)
  # mat = rbind(mat, c(26,27,28,29,30,0))
  nf <- layout(mat, widths = c(1,1,1,1,1,0.2), heights = c(1,1,1,1,1,1))
  layout.show(nf)
  # set.panel()
  # ind <- split.screen(c(6,nmts))
  for(ii in 1:30){
    par(mar = c(3.5 ,4.5, 1.5, 0.5), mgp = c(2,0.8,0))
    vegi = ceiling(ii /nmts)
    mi = mod0(ii, nmts)
    print(sprintf('%d,%d', vegi, mi))
    if (vegi == 1){
      tmain = sprintf("month = %d", mts[mi])
    } else {
      tmain = "";
    }
    if (mi == 1){
      tylb = paste(nms_veg[vegi], "p-pet",sep ='\n')
    } else {
      tylb = "p-pet"
    }
    Fdata <- xxx[[vegi, mi]]
    # plot(1:10)
    image(Fdata$x,Fdata$y,Fdata$z,zlim = params$zlm,col =params$col,
          main=tmain, cex.main=3,xlab="tmax", ylab= tylb, axes=T,
          cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2)
    # legend.args=list( text=""),smallplot=c(0.85,0.9,0.1,0.9)
  }
  par(mar = c(3.5 ,0.5, 1.5, 2.5), pty = "m", err = -1)
  breaks = linspace(0,1,length(params$col)+1) * max(params$zlm)
  ix <- 1:2
  iy <- breaks
  nBreaks <- length(breaks)
  midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
  iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
  image(ix,iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "",col =params$col,breaks = breaks)
  axis.args <- c(list(side = 4, mgp = c(3, 1, 0), 
                      las = 2, 
                      at = seq(0,1,0.1)* max(params$zlm)))
  do.call(axis,axis.args)
  dev.off();
}
