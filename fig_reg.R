d = d_reg
nms = c("GLASSgpp", "LUEgpp", "Musyqgpp", "Pgpp", "SIFgpp", "VPMgpp")
cols = c("#23FFDC", "#00008F", "#800000", "#FF4A00", "#005AFF", "#ECFF13")

nms_veg = c("forest","open forest","shrub land","paddy field","dry land","grassland")
par(mfrow = c(2,2), oma = c(5,4,1,1))

emf('./figures/relativeweight.emf')
avr = list_mean(d$SDgpp$output_reg$rws[1,8])
plt_raster(avr, insetmap = NULL,
           list(legend = T, clim = c(0,100), breaks = seq(0,100,10),
                titlename = "example month: relative importance in August"))
graphics.off()
emf('./figures/relativeweight_p.emf')

avr = list_mean(d$SDgpp$output_reg$rws[2,8])
plt_raster(avr, insetmap = NULL,
           list(legend = T, clim = c(0,100), breaks = seq(0,100,10),
                titlename = "month = 8, par cor"))
graphics.off()
avr = list_mean(d$SDgpp$output_reg$rws[1,1])
plt_raster(avr, insetmap = NULL,
           list(legend = T, clim = c(0,100), breaks = seq(0,100,10),
                titlename = "month = 1, par cor"))
avr = list_mean(d$SDgpp$output_reg$rws[2,1])
plt_raster(avr, insetmap = NULL,
           list(legend = T, clim = c(0,100), breaks = seq(0,100,10),
                titlename = "month = 1, par cor"))


out = matrix(NA, 141702, 12)
for (i in 1:12){
   out[,i] = c(d$SDgpp$output_reg$rws[[1,i]])
}

boxplot(out)

tet = computes_veg_month(d$SDgpp$output_reg$rws[1,])
tt = matrix(NA, 6, 12)
for (i in 1: 12){
  tt[,i] = tet[[i]]
}
a = colMeans(out, na.rm = T)
b = 100-a
plot(1:12, a/100, type = 'l', col = "red", lwd = 2, ylim  = c(0,1), ylab = "", xlab = "month")
par(new= T)
plot(1:12, b/100, type = "l", col = "blue", lwd = 2, ylim = c(0,1), xlab = "month", ylab = "relative weight")
# par(mfrow = c(1,1))
# boxplot(tt/100, ylim = c(-0.1, 1.1), xlab = "month", ylab = "relative weight",
#         main = "tmax", col = cols[i])


boxplot(t(tt[,5:9])/100, ylim = c(-0.1, 1.1), xlab = "", ylab = "relative weight",
        main = "red:tmax, blue:p-pet", col = "red", xaxt = "n")
par(new = T)
boxplot(1-t(tt[,5:9])/100, ylim = c(-0.1, 1.1), xlab = "", ylab = "relative weight",
        main = "", col = "blue", xaxt = "n")
axis(side = 1, at = 1:6, label = nms_veg, srt = 45, las = 2)
