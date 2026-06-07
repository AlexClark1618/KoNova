# distance in cm
# time in ps

speed.of.light <- 3.00e8 * 100 / 1e12 # in cm/ps

delta.z <- 13.6 # in cm from x layer on top to x layer on bottom.  Same for y.
channel.spacing <- 33 / 2 / 10 # in cm

cluster.gap <- 15000 # minimum gap between clusters in ps = 15 ns

min.hits.per.cluster <- 4 # >=

mapping <- scan( "Mapping.csv", skip = 1, sep = ",", what = list( bts.label = "", con1 = 0, or = 0, con2 = 0, samtec = 0, J2 = 0, ChNumber = 0, SiPMx = 0, SiPMy = 0 ) )
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

# Layer #1, X Top, Use after 12/27/2025
bar.map[ 0+1,1] <- 8
bar.map[10+1,1] <- 7
bar.map[ 5+1,1] <- 6
bar.map[14+1,1] <- 5
bar.map[20+1,1] <- 4
bar.map[23+1,1] <- 3
bar.map[17+1,1] <- 2
bar.map[18+1,1] <- 1

# Layer #2, Y Top
bar.map[48+1,2] <- 8
bar.map[47+1,2] <- 7
bar.map[44+1,2] <- 6
bar.map[42+1,2] <- 5
bar.map[51+1,2] <- 4
bar.map[54+1,2] <- 3
bar.map[55+1,2] <- 2
bar.map[60+1,2] <- 1

# Layer #3, X Bottom
bar.map[ 0+1,3] <- 8
bar.map[10+1,3] <- 7
bar.map[ 5+1,3] <- 6
bar.map[14+1,3] <- 5
bar.map[20+1,3] <- 4
bar.map[23+1,3] <- 3
bar.map[17+1,3] <- 2
bar.map[18+1,3] <- 1

# Layer #4, Y Bottom
bar.map[60+1,4] <- 1
bar.map[55+1,4] <- 2
bar.map[54+1,4] <- 3
bar.map[51+1,4] <- 4
bar.map[42+1,4] <- 5
bar.map[44+1,4] <- 6
bar.map[47+1,4] <- 7
bar.map[48+1,4] <- 8
