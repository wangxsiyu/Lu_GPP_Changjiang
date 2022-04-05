
## simulated data
nsim = 18000
td = c()
td$H = rnorm(nsim, mean = 1, sd = 1)
td$D = rnorm(nsim, mean = 2, sd = 1)
td$gpp = jitter(-td$H + td$D)
td = as.data.frame(td)


d<-dim(cvine$Matrix)[1]
cond1<-ff[,cvine$Matrix[(d+1)-1,(d+1)-1]]
cond2<-ff[,cvine$Matrix[(d+1)-2,(d+1)-2]]
condition<-cbind(cond1,cond2)# create conditional distribution matrix#
f1<-CDVineCondSim(cvine,condition)
f1 = data.frame(f1)
names(f1) = names(td1)

par(mfrow = c(1,3))
plot(f1$gpp, f1$H)
plot(f1$gpp, f1$D)
plot(f1$H, f1$D)
# library(plotly)
# fig <- plot_ly(td1, x = ~H, y = ~D, z = ~gpp,
#                marker = list(color = ~gpp, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE))
# fig <- fig %>% add_markers()
# fig <- fig %>% layout(scene = list(xaxis = list(title = 'heat'),
#                                    yaxis = list(title = 'drought'),
#                                    zaxis = list(title = 'gpp')),
#                       annotations = list(
#                         x = 1.13,
#                         y = 1.05,
#                         text = 'GPP',
#                         xref = 'paper',
#                         yref = 'paper',
#                         showarrow = FALSE
#                       ))
# fig



# {
#   par(mfrow=c(1,1))
#   fff<-matrix(0,nsss,4)
#   fff[,1] = sss
#   fff[,2] = 0.1
#   fff[,3] = 0.9
#   fff[,4] = RVinePDF(fff[,1:3],cvine)
#   plot(sss, fff[,4], type = 'l')
# }
# 
# 
# {
#   par(mfrow=c(1,3))
#   fff<-matrix(0,nsss,4)
#   fff[,1] = sss
#   fff[,2] = 0.1
#   fff[,3] = 0.9
#   ttt = RVinePIT(fff[,1:3],cvine)
#   plot(sss, ttt[,1], type = 'l')
#   plot(sss, ttt[,2], type = 'l')
#   plot(sss, ttt[,3], type = 'l')
# }

library(stars)
library(sf)
# sfObj <- st_as_sf(tf, coords=1:2, crs=st_crs(4326))
# plot(sfObj, axes = T)