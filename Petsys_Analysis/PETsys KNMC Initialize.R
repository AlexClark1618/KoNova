# distance in cm
# time in ps

speed.of.light <- 3.00e8 * 100 / 1e12 # in cm/ps

delta.z <- 30.0 # in cm from x layer on top to x layer on bottom.  Same for y.
channel.spacing <- 33 / 2 / 10 # in cm

cluster.gap <- 15000 # minimum gap between clusters in ps = 15 ns

min.hits.per.cluster <- 4 # >=

mapping <- scan( "Mapping KNMC.csv", skip = 1, sep = ",", what = list( bts.label = "", con1 = 0, or = 0, con2 = 0, samtec = 0, J2 = 0, ChNumber = 0, SiPMx = 0, SiPMy = 0 ) )
sort( mapping$ChNumber )
hist( mapping$ChNumber, breaks = seq( 65 ) - 2 )
plot( mapping$SiPMx, mapping$SiPMy, type = "n", xlim = c( 0, 9 ), ylim = c( 0, 9 ), main = "SiPM Channel Mapping\nFront View, Channel/Hamamatsu Label", xlab = "SiPMx", ylab = "SiPMy" )
text( mapping$SiPMx, mapping$SiPMy, paste0( mapping$ChNumber, "/", mapping$bts.label ), cex = 0.75 )
points( 0.25, 8, pch = 20 )
points( 0.25, 1, pch = 20 )
points( 8.75, 1, pch = 20 )

# Tells you which bar belongs to which channel
# bar.map[ch, layer]
bar.map <- matrix( -1e6, nrow = 64, ncol = 4 )

# Reduced channels go from 0 to 63
# Array indices go from 1 to 64

for ( j in 1:4 ) {
  bar <- 1
  for ( i in 0:63 ) {
    bar.map[i + 1,j] <- bar
    bar <- bar + 1
  }
}
