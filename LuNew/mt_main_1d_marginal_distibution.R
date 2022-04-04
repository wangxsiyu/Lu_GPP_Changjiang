rm(list = ls())
#library(robustloggamma)
library(fitdistrplus)
library(extRemes)
library(extraDistr)
## load raw data
setwd("/Users/mengtian/Desktop/sy_analysis/mt_analysis/function/")
source('../../function/tools_wang.R')
source('../../function/function_wang.R')
datadir = '../../../check_datamap/3source_matrix_data/'
data = loadraw(datadir, names = c("AVSDgpp","heatindex","spei"))

te1 = data[["AVSDgpp"]]
te2 = data$spei
te3 = data[["heatindex"]]
commonyr = intersect(intersect(te1$yr,te2$yr),te3$yr)
idx1 = which(te1$yr %in% commonyr)
idx2 = which(te2$yr %in% commonyr)
idx3 = which(te3$yr %in% commonyr)

## get pixel data
nx = dim(te1$data[[1]])[1]
ny = dim(te1$data[[1]])[2]
month = 6:9  ## select interested month
nmon = length(month)
nyear = length(commonyr)
marginal_gpp <- matrix(0,6,nmon) ## 6distribution * nmon
marginal_heat <- matrix(0,2,nmon)
rownames(marginal_gpp)<-c("norm","weibull","lognorm","gamma","logistic","gev")
rownames(marginal_heat) <- c("norm","gev")
colnames(marginal_gpp) = colnames(marginal_gpp) = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")[month]


### 1. marginal distribution count
for(xi in 1:nx){
  if (xi %% 100 == 1){
    print(sprintf("xi is:%d/%d", xi, nx))
  }
  for(yi in 1:ny){
    gpp = spaei = heatindex = matrix(NA,length(commonyr),12)
    rownames(gpp) = rownames(spaei) = rownames(heatindex) = commonyr
    for( mi in 1:12 ){
      gpp[,mi] = arrayfun(function(x) {te1$data[[x]][[xi,yi]]},intersect(idx1,which(te1$mt == mi)))
      spaei[,mi] = arrayfun(function(x) {te2$data[[x]][[xi,yi]]}, intersect(idx2,which(te2$mt == mi)))
      heatindex[,mi] = arrayfun(function(x) {te3$data[[x]][[xi,yi]]}, intersect(idx3,which(te3$mt == mi)))
    }
    ndim = length(commonyr)*12
    if(length(which(is.na(gpp)))<ndim & length(which(is.na(spaei)))<ndim & length(which(is.na(heatindex)))<ndim){
      
      ### 1. select marginal distribution -- get distribution parameters
      ## normal / 
      ## weibull / lognormal / gamma /  ---- not include  0 or NA
      ## logistic /  ---- not include   NA
      
      SPAEI_norm_par<-matrix(0,2,nmon)
      
      Heat_norm_par<-matrix(0,2,nmon)
      Heat_gev_par<-matrix(0,3,nmon)
      
      GPP_norm_par<-matrix(0,2,nmon)
      GPP_weibull_par<-matrix(0,2,nmon)
      GPP_lognorm_par<-matrix(0,2,nmon)
      GPP_gamma_par<-matrix(0,2,nmon)
      GPP_logis_par<-matrix(0,2,nmon)
      GPP_gev_par<-matrix(0,3,nmon)
      
      spaei_fit1 = heat_fit1 = heat_fit2 = gpp_fit1 = gpp_fit2 = gpp_fit3 = gpp_fit4 = gpp_fit5 = gpp_fit6 = matrix(list(),1,nmon)
      for (i in 1:nmon) {
        te = spaei[,month[i]]
        idx = which(is.na(te))
        if(length(idx)>0){
          te = te[-idx]
        }
        ## fit model
        spaei_fit1[[i]] = fitdist(te,"norm")
        
        ## spaei - normal
        SPAEI_norm_par[,i]<-fitdist(te,"norm")$estimate
      }
      
      for (i in 1:nmon) {
        te = heatindex[,month[i]]
        idx = which(is.na(te))
        if(length(idx)>0){
          te = te[-idx]
        }
        ## fit model
        #heat_fit1[[i]] = fitdist(te,"norm")
        heat_fit2[[i]] = fevd(te)
        
        ## heatindex - normal / gev
        #Heat_norm_par[,i]<-fitdist(te,"norm")$estimate
        Heat_gev_par[,i]<-fevd(te)$results$par
      }
      
      for (i in 1:nmon) {
        te =  gpp[,month[i]]
        idx = which( te == 0 | is.na(te))
        te[idx] = mean(te, na.rm = T)
        idx = which(is.na(te))
        if(length(idx)>0){
          te = te[-idx]
        }
        
        ## fit model
        gpp_fit1[[i]] = fitdist(te,"norm")
        gpp_fit2[[i]] = fitdist(te,"weibull")
        gpp_fit3[[i]] = fitdist(te,"lnorm")
        gpp_fit4[[i]] = fitdist(te,"gamma")
        gpp_fit5[[i]] = fitdist(te,"logis")
        gpp_fit6[[i]] = fevd(te)
        
        ## gpp - normal / weibull / lognormal / gamma / logistic / gev
        GPP_norm_par[,i]<-fitdist(te,"norm")$estimate
        GPP_weibull_par[,i]<-fitdist(te,"weibull", method="mle")$estimate
        GPP_lognorm_par[,i]<-fitdist(te,"lnorm")$estimate
        GPP_gamma_par[,i]<-fitdist(te,"gamma", method="mme")$estimate
        GPP_logis_par[,i]<-fitdist(te,"logis")$estimate
        GPP_gev_par[,i]<-fevd(te)$results$par
      }
      
      # ## wue - extreme value distribution -- http://tougao.ecoagri.ac.cn/html/zgstny/2019/12/2019-1206.htm#outline_anchor_8
      # library(extRemes)
      # fit = fevd(gpp[,month[i]])
      # par = fit$results$par
      # plot(fit)
      
      #########  2. Goodness-of-fit statistics
      ##### ks test
      ks<-matrix(0,8,nmon)
      for (i in 1:nmon) {
        ks[1,i]<-ks.test(x=jitter(spaei[,month[i]]),y="pnorm",SPAEI_norm_par[1,i],SPAEI_norm_par[2,i])[["p.value"]]
        
        he = heatindex[,month[i]]
        ks[2,i]<-ks.test(x=jitter(he),y="pgev",Heat_gev_par[1,i],Heat_gev_par[2,i],Heat_gev_par[3,i])[["p.value"]]
        
        ge = gpp[,month[i]]
        ks[3,i]<-ks.test(x=jitter(ge),y="pnorm",GPP_norm_par[1,i],GPP_norm_par[2,i])[["p.value"]]
        ks[4,i]<-ks.test(x=jitter(ge),y="pweibull",GPP_weibull_par[1,i],GPP_weibull_par[2,i])[["p.value"]]
        ks[5,i]<-ks.test(x=jitter(ge),y="plnorm",GPP_lognorm_par[1,i],GPP_lognorm_par[2,i])[["p.value"]]
        ks[6,i]<-ks.test(x=jitter(ge),y="pgamma",GPP_gamma_par[1,i],GPP_gamma_par[2,i])[["p.value"]]
        ks[7,i]<-ks.test(x=jitter(ge),y="plogis",GPP_logis_par[1,i],GPP_logis_par[2,i])[["p.value"]]
        ks[8,i]<-ks.test(x=jitter(ge),y="pgev",GPP_gev_par[1,i],GPP_gev_par[2,i],GPP_gev_par[3,i])[["p.value"]]
        
      }
      colnames(ks)<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")[month]
      rownames(ks)<-c("norm",
                      "gev",
                      "norm","weibull","lognorm","gamma","logistic","gev")#,"robust gamma")
      
      ##### aic
      aic_re<-matrix(0,8,nmon)
      for (i in 1:nmon) {
        aic_re[1,i] = gofstat(spaei_fit1[[i]],fitnames=c("norm"))$aic
        
#        aic_re[2,i] = gofstat(heat_fit1[[i]],fitnames=c("norm"))$aic
        aic_re[2,i] = summary(heat_fit2[[i]],silent=TRUE)$AIC
      
        aic_re[3,i] = gofstat(gpp_fit1[[i]],fitnames=c("norm"))$aic
        aic_re[4,i] = gofstat(gpp_fit2[[i]],fitnames=c("weibull"))$aic
        aic_re[5,i] = gofstat(gpp_fit3[[i]],fitnames=c("lnorm"))$aic
        aic_re[6,i] = gofstat(gpp_fit4[[i]],fitnames=c("gamma"))$aic
        aic_re[7,i] = gofstat(gpp_fit5[[i]],fitnames=c("logis"))$aic
        aic_re[8,i] = summary(gpp_fit6[[i]],silent=TRUE)$AIC
      }

      
      ##### marginal count
      marginaldist<-matrix(0,3,nmon)
      for (i in 1:nmon) {
        marginaldist[1,i]<-which.max(ks[1,i])
        marginaldist[2,i]<-which.min(aic_re[2,i])
        marginaldist[3,i]<-which.max(ks[3:8,i])
      }
      rownames(marginaldist)<-c("SPAEI","Heatindex","GPP")
      
      # ###### check numbers of different theoretical marginal distributions for each pixel
      # for( i in 1:nmon ){
      #   id = which((1:6 %in% marginaldist[3,i]) == T) 
      #   marginal_gpp[id,i] = marginal_gpp[id,i]+1
      #   
      #   id2 = which((1:6 %in% marginaldist[2,i]) == T) 
      #   marginal_heat[id2,i] = marginal_heat[id2,i]+1
      # }
      
      ######  1. get parameter
      SPAEImarginal<-matrix(0,nyear,nmon)
 
      Heatmarginal<-matrix(0,nyear,nmon)
      
      GPPmarginal<-matrix(0,nyear,nmon)
      
      
      #####  2. calculte marginal distribution
      ## "norm","weibull","lognorm","gamma","logistic","gev"
      for (i in 1:nyear) {
        ## spaei
        for (j in 1:nmon) {
          se = spaei[i,month[j]]
          SPAEImarginal[i,j]<-pnorm(se,SPAEI_norm_par[1,j],SPAEI_norm_par[2,j])
        }
        ## heat
        for (j in 1:nmon) {
          he = heatindex[i,month[j]]
          Heatmarginal[i,j]<-pgev(he,Heat_gev_par[1,j],Heat_gev_par[2,j],Heat_gev_par[3,j])
        }
        
        ## gpp "norm","weibull","lognorm","gamma","logistic","gev"
        for (j in 1:nmon) {
          ge = gpp[i,month[j]]
          if(marginaldist[2,j]==1){
            GPPmarginal[i,j]<-pnorm(ge,GPP_norm_par[1,j],GPP_norm_par[2,j])
          }else if(marginaldist[2,j]==2){
            GPPmarginal[i,j]<-pweibull(ge,GPP_weibull_par[1,j],GPP_weibull_par[2,j])
          }else if(marginaldist[2,j]==3){
            GPPmarginal[i,j]<-plnorm(ge,GPP_lognorm_par[1,j],GPP_lognorm_par[2,j])
          }else if(marginaldist[2,j]==4){
            GPPmarginal[i,j]<-pgamma(ge,GPP_gamma_par[1,j],GPP_gamma_par[2,j])
          }else if(marginaldist[2,j]==5){
            GPPmarginal[i,j]<-plogis(ge,GPP_logis_par[1,j],GPP_logis_par[2,j])
          }else{
            GPPmarginal[i,j]<-pgev(ge,GPP_gev_par[1,j],GPP_gev_par[2,j],GPP_gev_par[3,j])
          }
        }
      }
      
      group = matrix(list(),1,nmon)
      for( i in 1:nmon ){
        group[[i]] = data.frame(Gpp = GPPmarginal[,i],Spaei = SPAEImarginal[,i],Heatindex = Heatmarginal[,i])
      }
      names(group)<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")[month]
      tfilename = paste0("../001copula/data/",xi,"_",yi,"_1d_marginal_par_6-9.RData")
      #if (!file.exists(tfilename)){
      save(marginaldist,SPAEI_norm_par,Heat_gev_par,
           GPP_norm_par,
           GPP_weibull_par,
           GPP_lognorm_par,
           GPP_gamma_par,
           GPP_logis_par,
           GPP_gev_par,ks,aic_re,
           file = tfilename)
      save(group,file = paste0("../001copula/data/",xi,"_",yi,"_1d_marginal_distribution_6-9.RData"))
      #}

      
    }
  }
}     
      
      
      

