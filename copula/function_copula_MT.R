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
