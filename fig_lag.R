nms = c("GLASSgpp", "LUEgpp", "Musyqgpp", "Pgpp", "SIFgpp", "VPMgpp")
cols = c("#23FFDC", "#00008F", "#800000", "#FF4A00", "#005AFF", "#ECFF13")

nms_veg = c("Forest","Open Forest","Shrub Land","Paddy Field","Dry Land","Grassland")


d = d_lag
d = d$SDgpp$output_reg
# significance of the models/explained variance
# avr = list_mean(d$rs[8])

#library(devEMF)
#emf('./figures/explainedvar.emf')
#plt_raster(avr, insetmap = NULL,
           # list(legend = T, clim = c(0,1), breaks = seq(0,1,0.1),
           #      titlename = "month = 8, R^2"))
#graphics.off()
#avp = list_mean(d$ps[8])
#emf('./figures/explainedvar_p.emf')
#plt_raster(avp, insetmap = NULL,
#           list(legend = F, clim = c(0,0.05), breaks = c(0,0.001, 0.01, 0.05),
#                titlename = "month = 8, p < 0.05", colors = c("red","orange","yellow","white")))
#graphics.off()

# par(mfrow = c(1,2))
# # emf('./figures/explainedvar_bymonth.emf')
# te = computes_veg_month(d$rs)[[1]]
# boxplot(te, xlab = "month", ylab = "R^2")
# # graphics.off()
# # emf('./figures/explainedvar_byveg.emf')
# boxplot(t(te[,5:9]), ylab = "R^2", xaxt = "n", col = cols)
# axis(side = 1, at = 1:6, label = nms_veg, srt = 45, las = 2)
# graphics.off()
# emf('./figures/explainedvar_bymonth.emf')

# te = computes_veg_month(d$ps, 0.05)[[1]]
# boxplot(te, xlab = "month", ylab = "percent significant")
#
# boxplot(t(te[,5:9]), ylab = "percent significant", xaxt = "n", col = cols)
# axis(side = 1, at = 1:6, label = nms_veg, srt = 45, las = 2)

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
# par(mfrow = c(2,3))
#
# for (i in 1:6){
#   tt = tet[[i]]
#   tp = tep[[i]]
#   boxplot(tt[,], ylim = c(-0.1, 0.7), xlab = "lag", ylab = "coef",
#           main = sprintf('lag = %d', i-1), col = cols[i])
# }
#
#
#
# emf('./figures/lag_tmax.emf')
tet = lagveg2veglag(tet)
tep = lagveg2veglag(tep)
# par(mfrow = c(2,3), mai = c(1, .1, 0.2, 0.2),
    # oma = c(3,3,0,) + 0.1,
    # mar = c(2,2,1,1) + 0.1)
W_figure(2,3)
for (i in 1:6){
  tt = t(tet[[i]])
  tp = t(tep[[i]])
  boxplot(tp[5:9,], at = 1:6, names = 0:5,ylim = c(-0.2, 0.5), xlab = "Month (lag)", ylab = "Standard coefficients",
          main = nms_veg[i], col = cols[i], cex.main = 2, cex.axis = 1.2, cex.lab = 1.5)
}
# graphics.off()


pp = d$coef_standreg[(1:6)+6,]
rmax = mmax = matrix(list(),1,12)
for (mi in 1:12){
  rmax[[mi]] = pp[[1, mi]]
  mmax[[mi]] = rmax[[mi]] * 0
  for (i in 2:6){
     tid = pp[[i,mi]] < rmax[[mi]]
     tid = which(tid)
     mmax[[mi]][tid] = i-1
     rmax[[mi]][tid] = pp[[i,mi]][tid]
  }
}

tet = computes_veg_month(rmax[1,])
tet2 = computes_veg_month(mmax[1,])
tt = tt2 = matrix(NA, 6, 12)
for (i in 1: 12){
  tt[,i] = tet[[i]]
  tt2[,i] = tet2[[i]]
}
out =out2 = matrix(NA, 141702, 12)
for (i in 1:12){
  out[,i] = c(rmax[[i]])
  out2[,i] = c(mmax[[i]])
}
a = W_av(out)
a2 = W_av(out2)
dev.off()

{
W_figure(2,2)
W_barplot(a,xlab = "Month", ylab ="Peak Partial correlation", main ="GPP vs P-PET by Month",
          ylim = c(0,0.7), names.arg = 1:12, cex.axis = 1.2, cex.lab = 1.5)
plt <- boxplot(t(tt[,5:9]), ylim = c(0, 0.7), xlab = "", ylab = "Peak Partial correlation",
               main = "GPP vs P-PET by Vegetation", col = cols, xaxt = "n", cex.lab = 1.5, cex.axis = 1.2)
text(1:6, par("usr")[3], labels = nms_veg, srt = 45, adj = c(1.1,1.1),
     xpd = TRUE, cex=1.2)
W_barplot(a2,xlab = "Month", ylab ="Peak lag", main ="Peak lag by Month",
          ylim = c(0,5), names.arg = 1:12, cex.axis = 1.2, cex.lab = 1.5)
plt <- boxplot(t(tt2[,5:9]), ylim = c(0, 5), xlab = "", ylab = "Peak lag",
               main = "Peak lag by Vegetation", col = cols, xaxt = "n", cex.lab = 1.5, cex.axis = 1.2)
text(1:6, par("usr")[3], labels = nms_veg, srt = 45, adj = c(1.1,1.1),
     xpd = TRUE, cex=1.2)
}


###
load('./output/lagSPEI/lagSPEI_SDgpp.RData')
d = output_reg

tet = computes_veg_month(d$coef_partialcor[,])
tet = lagveg2veglag(tet)
tet = computes_veg_month(d$coef_partialcor[,])
tet = lagveg2veglag(tet)
W_figure(2,3)
od = order(d$lag)
for (i in 1:6){
  par(mar = c(3,3,5,2) + 0.1,xpd = NA)
  tt = tet[[i]]
  tt = t(tt)
  tt = tt[5:9,]
  boxplot(tt[,od], names = NA,ylim = c(-0.2, 0.4),
          xlab = "Month (lag)", ylab = "Partial correlation",
          main = nms_veg[i], col = cols[i], cex.main = 2* 1.5, cex.axis = 1.2* 1.5, cex.lab = 1.5 * 1.5)
  axis(side = 1,at = 1:length(od),labels = d$lag[od], cex.axis = 1.1, cex.lab = 2.25)
}

