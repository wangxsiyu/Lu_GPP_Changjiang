convert2marginal <- function(xs, pars, distname){
    for (i in 1:length(xs)){
      if(distname == "ecdf"){
        xs[[i]] = ecdf(pars)(xs[[i]])
      }else{
        xs[[i]] = getmarginal(xs[[i]], pars, distname)
      }
    }
  return(xs)
}

compute_spatial <- function(map, xp, yp, savename){
### compute conditional prob
library(akima) 
library(fields) 
################ compute marginal
nmts = length(map)
cvine = d1 = veg1 = heat1 = dry1 = matrix(list(), 6,nmts)
for (vegi in 1:6){
  for (mi in 1:length(mts)){
    print(sprintf('%d,%d', vegi, mi))
    idx = vegmap == vegi
    tveg = map[[mi]][idx]
    theat = xp[[mi]][idx]
    tdry = yp[[mi]][idx]
    td = data.frame(gpp = tveg, H = theat, D = tdry)
    td = td[colMeans(t(is.na(td))) == 0,]
    d1[[vegi, mi]] = td
    ## compute marginal
    # veg1[[vegi, mi]] = MT_copula_marginal(td$gpp)
    # heat1[[vegi, mi]] = MT_copula_marginal(td$H)
    # dry1[[vegi, mi]] = MT_copula_marginal(td$D, c(1,5,6))
    # d1[[vegi, mi]] = data.frame(gpp = veg1[[vegi, mi]]$x, H = heat1[[vegi, mi]]$x, D = dry1[[vegi, mi]]$x)
    ## fit cvine
    ff = as.matrix(d1[[vegi, mi]])
    cvine[[vegi, mi]]<-CDVineCondFit(ff,Nx=2,c(1,2,3,4,5,9),rotation=F,treecrit="AIC",type="CVine",selectioncrit="AIC")
  }
}
save(cvine, file = savename)
}




compute_xyz_spatial_emp <- function(map, xp, mts, savename){
  ### compute conditional prob
  library(akima)
  library(fields)
  ################ compute marginal
  nmts = length(mts)
  d0 = d1 = veg1 = heat1 = dry1 = matrix(list(), 6,nmts)
  for (vegi in 1:6){
    for (mi in 1:length(mts)){
      print(sprintf('%d,%d', vegi, mi))
      idx = vegmap == vegi
      tveg = map[[mi]][idx]
      theat = xp[[1, mi]][idx]
      tdry = xp[[2, mi]][idx]
      td = data.frame(gpp = tveg, H = theat, D = tdry)
      td = td[colMeans(t(is.na(td))) == 0,]
      d0[[vegi, mi]] = td
      ## compute marginal
      veg1[[vegi, mi]]$raw = ecdf(td$gpp)
      veg1[[vegi, mi]]$x = ecdf(td$gpp)(td$gpp)
      veg1[[vegi, mi]]$distnames = "ecdf"
      veg1[[vegi, mi]]$idbest = 1
      veg1[[vegi, mi]]$pars = list(td$gpp)
      heat1[[vegi, mi]]$raw = ecdf(td$H)
      heat1[[vegi, mi]]$x = ecdf(td$H)(td$H)
      heat1[[vegi, mi]]$distnames = "ecdf"
      heat1[[vegi, mi]]$idbest = 1
      heat1[[vegi, mi]]$pars = list(td$H)
      dry1[[vegi, mi]]$raw = ecdf(td$D)
      dry1[[vegi, mi]]$x = ecdf(td$D)(td$D)
      dry1[[vegi, mi]]$distnames = "ecdf"
      dry1[[vegi, mi]]$idbest = 1
      dry1[[vegi, mi]]$pars = list(td$D)
      d1[[vegi, mi]] = data.frame(gpp = veg1[[vegi, mi]]$x, H = heat1[[vegi, mi]]$x, D = dry1[[vegi, mi]]$x)
    }
  }
  save(d0, d1, dry1, veg1, heat1, file = savename)
}
## old codes
compute_xyz_spatial <- function(map, xp, mts, nmgpp, nmh, nmd, savename){
  ### compute conditional prob
  library(akima)
  library(fields)
  ################ compute marginal
  nmts = length(mts)
  d0 = d1 = veg1 = heat1 = dry1 = matrix(list(), 6,nmts)
  for (vegi in 1:6){
    for (mi in 1:length(mts)){
      print(sprintf('%d,%d', vegi, mi))
      idx = vegmap == vegi
      tveg = map[[mi]][idx]
      theat = xp[[1, mi]][idx]
      tdry = xp[[2, mi]][idx]
      td = data.frame(gpp = tveg, H = theat, D = tdry)
      td = td[colMeans(t(is.na(td))) == 0,]
      d0[[vegi, mi]] = td
      ## compute marginal
      veg1[[vegi, mi]] = MT_copula_marginal(td$gpp)
      heat1[[vegi, mi]] = MT_copula_marginal(td$H)
      dry1[[vegi, mi]] = MT_copula_marginal(td$D, c(1,5,6))
      d1[[vegi, mi]] = data.frame(gpp = veg1[[vegi, mi]]$x, H = heat1[[vegi, mi]]$x, D = dry1[[vegi, mi]]$x)
    }
  }
  out = list()
  out$d0 = d0
  out$d1 = d1
  out$dry1 = dry1
  out$heat1 = heat1
  out$veg1 = veg1
  save(d0, d1, dry1, veg1, heat1, file = savename)
}

compute_spatial_old <- function(d1, savename){
  ### compute conditional prob
  library(akima)
  library(fields)
  ################ compute marginal
  nmts = dim(d1)[2]
  cvine = matrix(list(), 6,nmts)
  for (vegi in 1:6){
    for (mi in 1:nmts){
      print(sprintf('%d,%d', vegi, mi))
       ## fit cvine
      ff = as.matrix(d1[[vegi, mi]])
      cvine[[vegi, mi]]<-CDVineCondFit(ff,Nx=2,c(1,2,3,4,5,9),rotation=F,treecrit="AIC",type="CVine",selectioncrit="AIC")
    }
  }
  save(cvine, file = savename)
}

simulate_spatial <- function(d1, cvine, savename){
    nmts = dim(d1)[2]
    d2 = matrix(list(), 6,nmts)
    for (vegi in 1:6){
      for (mi in 1:nmts){
        print(sprintf('%d,%d', vegi, mi))
        d<-dim(cvine[[vegi, mi]]$Matrix)[1]
        ff = d1[[vegi, mi]]
        cond1<-ff[,cvine[[vegi, mi]]$Matrix[(d+1)-1,(d+1)-1]]
        cond2<-ff[,cvine[[vegi, mi]]$Matrix[(d+1)-2,(d+1)-2]]
        condition<-cbind(cond1,cond2)# create conditional distribution matrix#
        f1<-CDVineCondSim(cvine[[vegi, mi]],condition)
        f1 = data.frame(f1)
        names(f1) = names(ff)
        d2[[vegi, mi]] = f1
      }
    }
    save(d2, file = savename)
    return(d2)
}
# 
plt_3d_copula <- function(td1){
  library(plotly)
  fig <- plot_ly(td1, x = ~H, y = ~D, z = ~gpp,
                 marker = list(color = ~gpp, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE))
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'heat'),
                                     yaxis = list(title = 'drought'),
                                     zaxis = list(title = 'gpp')),
                        annotations = list(
                          x = 1.13,
                          y = 1.05,
                          text = 'GPP',
                          xref = 'paper',
                          yref = 'paper',
                          showarrow = FALSE
                        ))
  fig
}

