d = d_lag
nms = c("GLASSgpp", "LUEgpp", "Musyqgpp", "Pgpp", "SIFgpp", "VPMgpp")
cols = c("#23FFDC", "#00008F", "#800000", "#FF4A00", "#005AFF", "#ECFF13")

nms_veg = c("forest","open forest","shrub land","paddy field","dry land","grassland")
d = d$SDgpp$output_reg
# significance of the models/explained variance
avr = list_mean(d$rs[8])

library(devEMF)
emf('./figures/explainedvar.emf')
plt_raster(avr, insetmap = NULL,
           list(legend = T, clim = c(0,1), breaks = seq(0,1,0.1),
                titlename = "month = 8, R^2"))
graphics.off()
avp = list_mean(d$ps[8])
emf('./figures/explainedvar_p.emf')
plt_raster(avp, insetmap = NULL,
           list(legend = F, clim = c(0,0.05), breaks = c(0,0.001, 0.01, 0.05),
                titlename = "month = 8, p < 0.05", colors = c("red","orange","yellow","white")))
graphics.off()

par(mfrow = c(1,2))
# emf('./figures/explainedvar_bymonth.emf')
te = computes_veg_month(d$rs)[[1]]
boxplot(te, xlab = "month", ylab = "R^2")
# graphics.off()
# emf('./figures/explainedvar_byveg.emf')
boxplot(t(te[,5:9]), ylab = "R^2", xaxt = "n", col = cols)
axis(side = 1, at = 1:6, label = nms_veg, srt = 45, las = 2)
# graphics.off()
# emf('./figures/explainedvar_bymonth.emf')

te = computes_veg_month(d$ps, 0.05)[[1]]
boxplot(te, xlab = "month", ylab = "percent significant")

boxplot(t(te[,5:9]), ylab = "percent significant", xaxt = "n", col = cols)
axis(side = 1, at = 1:6, label = nms_veg, srt = 45, las = 2)

# R vs Lag by vegetation

lagveg2veglag <- function(te){
  out = matrix(list(),1,6)
  nn = length(te)
  for(i in 1:6){
    out[[i]] = matrix(NA,nn, 12)
    for (j in 1:nn){
      out[[i]][j,] = te[[j]][i,]
    }
  }
  return(out)
}

tet = computes_veg_month(d$coef_standreg[1:6,])
tep = computes_veg_month(d$coef_standreg[(1:6)+6,])
# par(mfrow = c(2,3), mai = c(1, .1, 0.2, 0.2),
#     oma = c(3,3,0,0) + 0.1,
    # mar = c(2,2,1,1) + 0.1)
par(mfrow = c(2,3))

for (i in 1:6){
  tt = tet[[i]]
  tp = tep[[i]]
  boxplot(tp[,], ylim = c(-0.1, 0.7), xlab = "lag", ylab = "coef",
          main = sprintf('lag = %d', i-1), col = cols[i])
}


# 
# emf('./figures/lag_tmax.emf')
tet = lagveg2veglag(tet)
tep = lagveg2veglag(tep)
# par(mfrow = c(2,3), mai = c(1, .1, 0.2, 0.2),
    # oma = c(3,3,0,) + 0.1,
    # mar = c(2,2,1,1) + 0.1)
par(new = F)
par(mfrow = c(2,3))
for (i in 1:6){
  tt = t(tet[[i]])
  tp = t(tep[[i]])
  boxplot(tp[5:9,], ylim = c(-0.1, 0.5), xlab = "lag", ylab = "coef",
          main = nms_veg[i], col = cols[i])
}
# graphics.off()

###
load('./output/lagSPEI/lagSPEI_SDgpp.RData')
d = output_reg

tet = computes_veg_month(d$coef_partialcor[,])
tet = lagveg2veglag(tet)
par(mfrow = c(2,3), mai = c(1, .1, 0.2, 0.2),
    oma = c(3,3,0,0) + 0.1,
    mar = c(2,2,1,1) + 0.1)
od = order(d$lag)
for (i in 1:6){
  tt = tet[[i]]
  tt = t(tt)
  tt = tt[5:9,]
  boxplot(tt[,od],labels = d$lag[od], names = d$lag[od], ylim = c(-0.1, 0.7), xlab = "lag", ylab = "coef",
          main = nms_veg[i], col = cols[i])
}
