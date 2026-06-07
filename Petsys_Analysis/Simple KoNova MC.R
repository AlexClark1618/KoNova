# Simple MC
# Initialize R

setwd( "/Users/ifft/Desktop/Muography/PETsys Analysis" )

rm( list = ls() )

number.muons <- 10000000
spacing <- 3.3 / 2
L.detector <- ( spacing * 63 ) + 2 * spacing # maximum length over which a muon can touch a detector in cm
L.start <- L.detector # in cm
exposure.time <- number.muons / L.start^2 / 0.01099556
z1 <- 0
# z2 <- 15.25 * 2.54 - 2 * 0.7375 # from bottom of bottom scintillator to top of top scintillator
# z2
z2 <- 11.75 * 2.54 # in cm from x layer on top to x layer on bottom.  Same for y.  Measured 4/29/2026 see Photos


# integral.rate <- 0.01099556 * L.detector^2 # for no gap
# integral.rate
# integral.rate <- 0.0003357554 # for 1000 cm gap
# integral.rate
# integral.rate <- ( 70 / 100^2 ) * ( L.detector^2 / (z2 - z1)^2 ) * L.detector^2
# integral.rate
# integral.rate <- 0.6974 # for 10 cm gap
# integral.rate
# 
# integral.rate / mc.rate

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
y2 <- ran.y + z2 * tan( ran.zeniths ) * sin( ran.azimuths ) + 2

stdcut <- x1 >= -L.detector / 2 & x1 <= L.detector / 2
stdcut <- stdcut & y1 >= -L.detector / 2 & y1 <= L.detector / 2
stdcut <- stdcut & x2 >= -L.detector / 2 & x2 <= L.detector / 2
stdcut <- stdcut & y2 >= -L.detector / 2 & y2 <= L.detector / 2
length( ran.zeniths[stdcut] )

r.max <- ( L.detector - 2 * spacing ) / 2
r1 <- sqrt( x1^2 + y1^2 )
r2 <- sqrt( x2^2 + y2^2 )
r.cut <- r1 < r.max & r2 < r.max
stdcut <- stdcut & r.cut
length( ran.zeniths[stdcut] )

exposure.time
mc.rate <- length( ran.zeniths[stdcut] ) / exposure.time
mc.rate
sqrt( length( ran.zeniths[stdcut] ) ) / exposure.time

main.label <- paste0( "Theoretical Zeniths\nRate = ", signif( mc.rate, 4 ), " Hz" )
hist( ran.zeniths[stdcut] * 180 / pi, nclass = 100, main = main.label )
h <- hist( ran.zeniths[stdcut] * 180 / pi, breaks = seq( from = 0, to = 90, by = 1 ), plot = FALSE )
write( t( cbind( h$mids, h$density ) ), file = "KoNova Zenith Response.csv", ncolumns = 2, sep = "," )
sum( h$density )


hist( ran.azimuths[stdcut] * 180 / pi, nclass = 50 )
h <- hist( ran.azimuths[stdcut] * 180 / pi, nclass = 50, plot = FALSE )
max( h$counts ) / min( h$counts )


