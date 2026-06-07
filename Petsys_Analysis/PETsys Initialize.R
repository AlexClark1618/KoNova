# distance in cm
# time in ps

if( !require( fs ) )install.packages( "fs" )
library( fs )

speed.of.light <- 3.00e8 * 100 / 1e12 # in cm/ps

delta.z <- 11.75 * 2.54 # in cm from x layer on top to x layer on bottom.  Same for y.  Measured 4/29/2026 see Photos
channel.spacing <- 33 / 2 / 10 # in cm == L/2
layer.height <- 17 / 10 # in cm == h

cluster.gap <- 30000 # minimum gap between clusters in ps = 30 ns

min.hits.per.cluster <- 4 # >=

layer.numbers <- 1:4 # assumes these are ordered from bottom to top

# MPPC mapping
mapping <- scan( paste0( program.root, "/Mapping.csv" ), skip = 1, sep = ",", what = list( bts.label = "", con1 = 0, or = 0, con2 = 0, samtec = 0, J2 = 0, ChNumber = 0, SiPMx = 0, SiPMy = 0 ) )
sort( mapping$ChNumber )
# hist( mapping$ChNumber, breaks = seq( 65 ) - 2 )
# plot( mapping$SiPMx, mapping$SiPMy, type = "n", xlim = c( 0, 9 ), ylim = c( 0, 9 ), main = "SiPM Channel Mapping\nFront View, Channel/Hamamatsu Label", xlab = "SiPMx", ylab = "SiPMy" )
# text( mapping$SiPMx, mapping$SiPMy, paste0( mapping$ChNumber, "/", mapping$bts.label ), cex = 0.75 )
# points( 0.25, 8, pch = 20 )
# points( 0.25, 1, pch = 20 )
# points( 8.75, 1, pch = 20 )

# Bar mapping
ch.2.bar.mat <- matrix( 0, nrow = 64, ncol = length( layer.numbers ) )

for ( layer.number in layer.numbers ) {
  mat <- matrix( scan( file = paste0( program.root, "/Layer Maps/bar_ch_map_layer", layer.number, ".csv"  ), skip = 1, sep = "," ), ncol = 2, byrow = TRUE )
  # mat[,1] = channel number
  # mat[,2] = bar number
  ch.2.bar.mat[,layer.number] <- mat[,2] # take a channel, add +1, then read off bar number
}

layer.min.channel.number <- seq( layer.numbers ) * 0 
layer.max.channel.number <- seq( layer.numbers ) * 0 

layer.min.channel.number[1] <- 0
layer.max.channel.number[1] <- 200

layer.min.channel.number[2] <- 200
layer.max.channel.number[2] <- 500

layer.min.channel.number[3] <- 500
layer.max.channel.number[3] <- 700

layer.min.channel.number[4] <- 700
layer.max.channel.number[4] <- 1e6

layer.subtract.this <- c( 64, 320, 576, 832 )

layer.qdc.offsets <- c( 1.31481355062237, 1.35048687884656, 0.868484977921783, 1.36229120996146 )

# Post stats analysis parameters.  These can be changed after get hit stats and get events stats have been run.

delta.x.add <- 0.01 # in cm
delta.y.add <- 0.35 # in cm

delta.x.stretch <- 1.003
delta.y.stretch <- 0.9965


