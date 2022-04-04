d = d_cor
nms = c("GLASSgpp", "LUEgpp", "Musyqgpp", "Pgpp", "SIFgpp", "VPMgpp")
cols = c("#23FFDC", "#00008F", "#800000", "#FF4A00", "#005AFF", "#ECFF13")

nms_veg = c("forest","open forest","shrub land","paddy field","dry land","grassland")
par(mfrow = c(2,2), oma = c(5,4,1,1))
emf('./figures/parcor_tmax.emf')
avr = list_mean(d$SDgpp$output_reg$coef_partialcor[1,8])
plt_raster(avr, insetmap = NULL,
           list(legend = T, clim = c(-1,1), breaks = seq(-1,1,0.1),
                titlename = "example month: partial correlation of tmax in August"))
graphics.off()
emf('./figures/parcor_p_tmax.emf')
avp = list_mean(d$SDgpp$output_reg$pval_partialcor[1,8])
plt_raster(avp, insetmap = NULL,
           list(legend = F, clim = c(0,0.05), breaks = c(0,0.001, 0.01, 0.05),
                titlename = "month = 8, p < 0.05", colors = c("red","orange","yellow","white")))
graphics.off()


tet = computes_veg_month(d$SDgpp$output_reg$coef_partialcor[1,])
tt = matrix(NA, 6, 12)
for (i in 1: 12){
  tt[,i] = tet[[i]]
}
boxplot(t(tt[,5:9]), ylim = c(0.1, 0.6), xlab = "", ylab = "partial correlation",
        main = "tmax", col = cols, xaxt = "n")
axis(side = 1, at = 1:6, label = nms_veg, srt = 45, las = 2)

out = matrix(NA, 141702, 12)
for (i in 1:12){
  out[,i] = c(d$SDgpp$output_reg$coef_partialcor[[1,i]])
}

boxplot(out, outline = F)
a = colMeans(out, na.rm = T)

barplot(a,xlab = "month", ylab ="partial correlation", main ="tmax", ylim = c(0,0.5))

emf('./figures/parcor_ppet.emf')

avr = list_mean(d$SDgpp$output_reg$coef_partialcor[2,8])
plt_raster(avr, insetmap = NULL,
           list(legend = T, clim = c(-1,1), breaks = seq(-1,1,0.1),
                titlename = "month = 8, par cor"))
graphics.off()
emf('./figures/parcor_p_ppet.emf')

avp = list_mean(d$SDgpp$output_reg$pval_partialcor[2,8])
plt_raster(avp, insetmap = NULL,
           list(legend = F, clim = c(0,0.05), breaks = c(0,0.001, 0.01, 0.05),
                titlename = "month = 8, p < 0.05", colors = c("red","orange","yellow","white")))
graphics.off()



tet = computes_veg_month(d$SDgpp$output_reg$coef_partialcor[2,])
tt = matrix(NA, 6, 12)
for (i in 1: 12){
  tt[,i] = tet[[i]]
}
boxplot(t(tt[,5:9]), ylim = c(-0.1, 0.2), xlab = "", ylab = "partial correlation",
        main = "P-PET", col = cols, xaxt = "n")
axis(side = 1, at = 1:6, label = nms_veg, srt = 45, las = 2)

out = matrix(NA, 141702, 12)
for (i in 1:12){
  out[,i] = c(d$SDgpp$output_reg$coef_partialcor[[2,i]])
}

boxplot(out, outline = F)
a = colMeans(out, na.rm = T)

barplot(a,xlab = "month", ylab ="partial correlation", main ="P-PET", ylim = c(0,0.5))
