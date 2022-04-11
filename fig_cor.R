d = d_cor
nms = c("GLASSgpp", "LUEgpp", "Musyqgpp", "Pgpp", "SIFgpp", "VPMgpp")
cols = c("#23FFDC", "#00008F", "#800000", "#FF4A00", "#005AFF", "#ECFF13")

nms_veg = c("Forest","Open Forest","Shrub Land","Paddy Field","Dry Land","Grassland")
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

out = matrix(NA, 141702, 12)
for (i in 1:12){
  out[,i] = c(d$SDgpp$output_reg$coef_partialcor[[1,i]])
}

boxplot(out, outline = F)
a = W_av(out)
# b = W_se(out)
dev.off()
W_figure(1,2)
W_barplot(a,xlab = "Month", ylab ="Partial correlation", main ="GPP vs Tmax by Month",
          ylim = c(0,0.7), names.arg = 1:12, cex.axis = 1.2, cex.lab = 1.5)
plt <- boxplot(t(tt[,5:9]), ylim = c(0, 0.7), xlab = "", ylab = "Partial correlation",
        main = "GPP vs Tmax by Vegetation", col = cols, xaxt = "n", cex.lab = 1.5, cex.axis = 1.2)
# nms_veg = c("Forest","Open\nForest","Shrub\nLand","Paddy\nField","Dry\nLand","Grassland")
text(1:6, par("usr")[3], labels = nms_veg, srt = 45, adj = c(1.1,1.1),
     xpd = TRUE, cex=1.2)

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

out = matrix(NA, 141702, 12)
for (i in 1:12){
  out[,i] = c(d$SDgpp$output_reg$coef_partialcor[[2,i]])
}

boxplot(out, outline = F)
a = colMeans(out, na.rm = T)
W_figure(1,2)
W_barplot(a,xlab = "Month", ylab ="Partial correlation", main ="GPP vs P-PET by Month",
          ylim = c(-0.2, 0.5), names.arg = 1:12, cex.axis = 1.2, cex.lab = 1.5)
plt <- boxplot(t(tt[,5:9]), ylim = c(-0.2, 0.5), xlab = "", ylab = "Partial correlation",
               main = "GPP vs P-PET by Vegetation", col = cols, xaxt = "n", cex.lab = 1.5, cex.axis = 1.2)
text(1:6, par("usr")[3], labels = nms_veg, srt = 45, adj = c(1.1,1.1),
     xpd = TRUE, cex=1.2)
