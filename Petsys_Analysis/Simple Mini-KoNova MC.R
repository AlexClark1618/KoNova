# Simple MC
# Initialize R

setwd( "/Users/ifft/Desktop/Muography/PETsys Analysis" )

rm( list = ls() )

number.muons <- 1000000
L.detector <- 14.85 # in cm
L.start <- L.detector # in cm
exposure.time <- number.muons / L.start^2 / 0.01099556
z1 <- 0
z2 <- 2 * 1.7 + 10 + 2 * 1.7

# z1 <- 0
# z2 <- 0.1

zenith <- seq( from = 0, to = 85, by = 0.1 ) * pi / 180 # 90 degrees produces a warning
this.sum <- cumsum( cos( zenith )^2 * sin( zenith ) )
ran.sum <- runif( number.muons, min = 0, max = max( this.sum ) )
ran.zeniths <- approx( this.sum, zenith, ran.sum )$y
hist( ran.zeniths * 180 / pi )
h <- hist( ran.zeniths * 180 / pi, plot = FALSE )
lines( zenith * 180 / pi, cos( zenith )^2 * sin( zenith ) * max( h$counts )  / ( cos( 0.606 )^2 * sin( 0.606 ) ) )

ran.azimuths <- runif( number.muons, min = 0, max = 2 * pi )

ran.x <- runif( number.muons, min = -L.start / 2, max = L.start / 2 )
ran.y <- runif( number.muons, min = -L.start / 2, max = L.start / 2 )
#plot( ran.x, ran.y )

x1 <- ran.x + z1 * tan( ran.zeniths ) * cos( ran.azimuths )
y1 <- ran.y + z1 * tan( ran.zeniths ) * sin( ran.azimuths )

x2 <- ran.x + z2 * tan( ran.zeniths ) * cos( ran.azimuths )
y2 <- ran.y + z2 * tan( ran.zeniths ) * sin( ran.azimuths )

stdcut <- x1 >= -L.detector / 2 & x1 <= L.detector / 2
stdcut <- stdcut & y1 >= -L.detector / 2 & y1 <= L.detector / 2
stdcut <- stdcut & x2 >= -L.detector / 2 & x2 <= L.detector / 2
stdcut <- stdcut & y2 >= -L.detector / 2 & y2 <= L.detector / 2

length( ran.zeniths[stdcut] )
exposure.time
mc.rate <- length( ran.zeniths[stdcut] ) / exposure.time
mc.rate
sqrt( length( ran.zeniths[stdcut] ) ) / exposure.time
hist( ran.zeniths[stdcut] * 180 / pi )

integral.rate <- 0.01099556 * L.detector^2 # for no gap
integral.rate
integral.rate <- 0.0003357554 # for 1000 cm gap
integral.rate
integral.rate <- ( 70 / 100^2 ) * ( L.detector^2 / (z2 - z1)^2 ) * L.detector^2
integral.rate
integral.rate <- 0.6974 # for 10 cm gap
integral.rate

integral.rate / mc.rate

# plot( x1[stdcut], y1[stdcut], pch = 20, cex = 0.1 )
# 
# plot( x2[stdcut], y2[stdcut], pch = 20, cex = 0.1 )
hist( x2[stdcut] )
hist( y2[stdcut] )



