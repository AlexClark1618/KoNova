# Make plots directory
graphics.off()

this.plot.directory <- paste0( program.root, plot.root, "/", analysis.indicator )
if ( dir.exists( this.plot.directory ) ) {
  system( paste0( "rm -r ", shQuote( this.plot.directory ) ) )
}
dir.create( this.plot.directory )


main.label <- paste0( "Histogram of Event Delta Ts\n", analysis.indicator )
hist( event.delta.t[muon.stdcut] / 1000, xlab = "Event delta.t (ns)", main = main.label, nclass = 100 )
pdf( paste0( this.plot.directory, "/", analysis.indicator, "-DeltaT.pdf" ) )
hist( event.delta.t[muon.stdcut] / 1000, xlab = "Event delta.t (ns)", main = main.label, nclass = 100 )
dev.off()

main.label <- paste0( "Histogram of Event QDCs\n", analysis.indicator )
hist( event.qdc[muon.stdcut], xlab = "Event QDCs", main = main.label, nclass = 100 )
pdf( paste0( this.plot.directory, "/", analysis.indicator, "-DeltaT.pdf" ) )
hist( event.qdc[muon.stdcut], xlab = "Event QDCs", main = main.label, nclass = 100 )
dev.off()


# Layer statistics

# Hits

par( mfrow = c( 2, 2 ), oma = c( 0, 0, 3, 0 ) )
for ( layer in layer.numbers ) {
  this.hits <- layer.hits[muon.stdcut,layer]
  hist( this.hits[this.hits < 6], main = paste0( "Layer ", layer ) )
}
mtext( analysis.indicator, side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
par( mfrow = c( 1, 1 ), oma = c( 0, 0, 0, 0 ) )

pdf( paste0( this.plot.directory, "/", analysis.indicator, "-LayerHits.pdf" ) )
par( mfrow = c( 2, 2 ), oma = c( 0, 0, 3, 0 ) )
for ( layer in layer.numbers ) {
  this.hits <- layer.hits[muon.stdcut,layer]
  hist( this.hits[this.hits < 6], main = paste0( "Layer ", layer ) )
}
mtext( analysis.indicator, side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
par( mfrow = c( 1, 1 ), oma = c( 0, 0, 0, 0 ) )
dev.off()

# QDC

max.qdc <- 7.5
par( mfrow = c( 2, 2 ), oma = c( 0, 0, 3, 0 ) )
for ( layer in layer.numbers ) {
  x <- layer.qdc[muon.stdcut,layer]
  x <- x[x > 0 & x < max.qdc]
  hist( x, xlab = "QDC after cuts", main = paste0( "Layer = ", layer, "\nmean = ", round( mean( x ), 2 ) ), xlim = c( 0, max.qdc ), nclass = 50 )
}
mtext( analysis.indicator, side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
par( mfrow = c( 1, 1 ), oma = c( 0, 0, 0, 0 ) )

pdf( paste0( this.plot.directory, "/", analysis.indicator, "-LayerQDC.pdf" ) )
max.qdc <- 7.5
par( mfrow = c( 2, 2 ), oma = c( 0, 0, 3, 0 ) )
for ( layer in layer.numbers ) {
  x <- layer.qdc[muon.stdcut,layer]
  x <- x[x > 0 & x < max.qdc]
  hist( x, xlab = "QDC after cuts", main = paste0( "Layer = ", layer, "\nmean = ", round( mean( x ), 2 ) ), xlim = c( 0, max.qdc ), nclass = 50 )
}
mtext( analysis.indicator, side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
par( mfrow = c( 1, 1 ), oma = c( 0, 0, 0, 0 ) )
dev.off()

# XY

main.label <- "Bottom x Layer"
hist( x1[muon.stdcut & x1 > -1 & x1 < 120], main = main.label, xlab = "x1 in cm", nclass = 400 )
main.label <- "Bottom x Layer (Zoom)"
hist( x1[muon.stdcut & x1 > -1 & x1 < 25], main = main.label, xlab = "x1 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )
hist( x1[muon.stdcut & x1 > 80 & x1 < 107], main = main.label, xlab = "x1 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

main.label <- "Bottom y Layer"
hist( y1[muon.stdcut & y1 > -1 & y1 < 120], main = main.label, xlab = "y1 in cm", nclass = 400 )
main.label <- "Bottom y Layer (Zoom)"
hist( y1[muon.stdcut & y1 > -1 & y1 < 25], main = main.label, xlab = "y1 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )
hist( y1[muon.stdcut & y1 > 80 & y1 < 107], main = main.label, xlab = "y1 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

main.label <- "Top x Layer"
hist( x2[muon.stdcut & x2 > -1 & x2 < 120], main = main.label, xlab = "x2 in cm", nclass = 400 )
main.label <- "Top x Layer (Zoom)"
hist( x2[muon.stdcut & x2 > -1 & x2 < 25], main = main.label, xlab = "x2 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )
hist( x2[muon.stdcut & x2 > 80 & x2 < 107], main = main.label, xlab = "x2 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

main.label <- "Top y Layer"
hist( y2[muon.stdcut & y2 > -1 & y2 < 120], main = main.label, xlab = "y2 in cm", nclass = 400 )
main.label <- "Top y Layer (Zoom)"
hist( y2[muon.stdcut & y2 > -1 & y2 < 25], main = main.label, xlab = "y2 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )
hist( y2[muon.stdcut & y2 > 80 & y2 < 107], main = main.label, xlab = "y2 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

pdf( paste0( this.plot.directory, "/", analysis.indicator, "-XY.pdf" ) )
main.label <- "Bottom x Layer"
hist( x1[muon.stdcut & x1 > -1 & x1 < 120], main = main.label, xlab = "x1 in cm", nclass = 400 )
main.label <- "Bottom x Layer (Zoom)"
hist( x1[muon.stdcut & x1 > -1 & x1 < 25], main = main.label, xlab = "x1 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )
hist( x1[muon.stdcut & x1 > 80 & x1 < 107], main = main.label, xlab = "x1 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

main.label <- "Bottom y Layer"
hist( y1[muon.stdcut & y1 > -1 & y1 < 120], main = main.label, xlab = "y1 in cm", nclass = 400 )
main.label <- "Bottom y Layer (Zoom)"
hist( y1[muon.stdcut & y1 > -1 & y1 < 25], main = main.label, xlab = "y1 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )
hist( y1[muon.stdcut & y1 > 80 & y1 < 107], main = main.label, xlab = "y1 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

main.label <- "Top x Layer"
hist( x2[muon.stdcut & x2 > -1 & x2 < 120], main = main.label, xlab = "x2 in cm", nclass = 400 )
main.label <- "Top x Layer (Zoom)"
hist( x2[muon.stdcut & x2 > -1 & x2 < 25], main = main.label, xlab = "x2 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )
hist( x2[muon.stdcut & x2 > 80 & x2 < 107], main = main.label, xlab = "x2 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

main.label <- "Top y Layer"
hist( y2[muon.stdcut & y2 > -1 & y2 < 120], main = main.label, xlab = "y2 in cm", nclass = 400 )
main.label <- "Top y Layer (Zoom)"
hist( y2[muon.stdcut & y2 > -1 & y2 < 25], main = main.label, xlab = "y2 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )
hist( y2[muon.stdcut & y2 > 80 & y2 < 107], main = main.label, xlab = "y2 in cm", nclass = 100 )
abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )
dev.off()

# Get intercepts
delta.break <- 0.1
intercepts <- NULL
par( mfrow = c( 2, 2 ), oma = c(0, 0, 3, 0) )
for ( layer in layer.numbers ) {
  x <- 1 / cos( zenith[muon.stdcut] )
  y <- layer.qdc[muon.stdcut, layer]
  # sample.indices <- sample( seq( x ), 1000 )
  # plot( x[sample.indices], y[sample.indices], pch = 20, cex = 0.1 )
  
  mean.y <- NULL
  breaks = seq( from = 1.0, to = 2.0, by = delta.break )
  for ( this.break in breaks ) {
    stdcut <- x > this.break & x < this.break + delta.break & y > 0 & y < 15
    mean.y <- c( mean.y, mean( y[stdcut] ) )
    #hist( y[stdcut], nclass = 100, main = paste0( mean( y[stdcut] ) ) )
  }
  ls.fit <- lsfit( breaks + 0.05, mean.y )
  plot( breaks + 0.05, mean.y, xlab = "mean( 1 / cos( zenith ) )", ylab = "Layer QDC", main = paste0( "Layer ", layer, "\nIntercept = ", signif( ls.fit$coefficients[1], 4 ) ) )
  abline( ls.fit )
  intercepts <- c( intercepts, ls.fit$coefficients[1] )
}
mtext( analysis.indicator, side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
par( mfrow = c( 1, 1 ), oma = c(0, 0, 0, 0) )
cbind( layer.numbers, intercepts )

paste0( "Put this in Initialize file" )
paste0( "layer.qdc.offsets <- c( ", intercepts[1], ", ", intercepts[2], ", ", intercepts[3], ", ", intercepts[4], " )" )

pdf( paste0( this.plot.directory, "/", analysis.indicator, "-Intercepts.pdf" ) )
par( mfrow = c( 2, 2 ), oma = c(0, 0, 3, 0) )
for ( layer in layer.numbers ) {
  x <- 1 / cos( zenith[muon.stdcut] )
  y <- layer.qdc[muon.stdcut, layer]
  
  mean.y <- NULL
  breaks = seq( from = 1.0, to = 2.0, by = delta.break )
  for ( this.break in breaks ) {
    stdcut <- x > this.break & x < this.break + delta.break & y > 0 & y < 15
    mean.y <- c( mean.y, mean( y[stdcut] ) )
  }
  ls.fit <- lsfit( breaks + 0.05, mean.y )
  plot( breaks + 0.05, mean.y, xlab = "mean( 1 / cos( zenith ) )", ylab = "Layer QDC", main = paste0( "Layer ", layer, "\nIntercept = ", signif( ls.fit$coefficients[1], 4 ) ) )
  abline( ls.fit )
  intercepts <- c( intercepts, ls.fit$coefficients[1] )
}
mtext( analysis.indicator, side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
par( mfrow = c( 1, 1 ), oma = c(0, 0, 0, 0) )
dev.off()

# Get inter-fiber frequencies

par( mfrow = c( 2, 2 ), oma = c( 0, 0, 3, 0 ) )

hist( x1[muon.stdcut & x1 > -1 & x1 < 120] %% channel.spacing, main = main.label, xlab = "x1 in cm", nclass = 100 )
hist( y1[muon.stdcut & y1 > -1 & y1 < 120] %% channel.spacing, main = main.label, xlab = "y1 in cm", nclass = 100 )
hist( x2[muon.stdcut & x2 > -1 & x2 < 120] %% channel.spacing, main = main.label, xlab = "x2 in cm", nclass = 100 )
hist( y2[muon.stdcut & y2 > -1 & y2 < 120] %% channel.spacing, main = main.label, xlab = "y2 in cm", nclass = 100 )

mtext( analysis.indicator, side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
par( mfrow = c( 1, 1 ), oma = c( 0, 0, 0, 0 ) )

pdf( paste0( this.plot.directory, "/", analysis.indicator, "-InterFiber.pdf" ) )
par( mfrow = c( 2, 2 ), oma = c( 0, 0, 3, 0 ) )

hist( x1[muon.stdcut & x1 > -1 & x1 < 120] %% channel.spacing, main = main.label, xlab = "x1 in cm", nclass = 100 )
hist( y1[muon.stdcut & y1 > -1 & y1 < 120] %% channel.spacing, main = main.label, xlab = "y1 in cm", nclass = 100 )
hist( x2[muon.stdcut & x2 > -1 & x2 < 120] %% channel.spacing, main = main.label, xlab = "x2 in cm", nclass = 100 )
hist( y2[muon.stdcut & y2 > -1 & y2 < 120] %% channel.spacing, main = main.label, xlab = "y2 in cm", nclass = 100 )

mtext( analysis.indicator, side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
par( mfrow = c( 1, 1 ), oma = c( 0, 0, 0, 0 ) )
dev.off()

# This is for KoNova
# plot( x1, y1, main = "Top Detector", xlab = "x (cm)", ylab = "y (cm)", xlim = c( -3, 106 ), ylim = c( -3, 106 ), type = "n" )
# stdcut <- hits.cut & adjacent.cut
# points( xydetector.x[stdcut,1], xydetector.y[stdcut,1], col = "brown", pch = 20, cex = 0.1 )
# points( xydetector.x[muon.stdcut,1], xydetector.y[muon.stdcut,1], col = "blue", pch = 20, cex = 0.1 )
# 
# hist( xydetector.x[muon.stdcut,1], nclass = 100 )
# hist( xydetector.y[muon.stdcut,1], nclass = 100 )
# 
# plot( xydetector.x[,2], xydetector.y[,2], main = "Top Detector", xlab = "x (cm)", ylab = "y (cm)", xlim = c( -3, 106 ), ylim = c( -3, 106 ), type = "n" )
# stdcut <- hits.cut & adjacent.cut
# points( xydetector.x[stdcut,2], xydetector.y[stdcut,1], col = "brown", pch = 20, cex = 0.1 )
# points( xydetector.x[muon.stdcut,2], xydetector.y[muon.stdcut,1], col = "blue", pch = 20, cex = 0.1 )
# 
# hist( xydetector.x[muon.stdcut,2], nclass = 100 )
# hist( xydetector.y[muon.stdcut,2], nclass = 100 )

mat <- matrix( scan( "./konova_v1_zenith_distribution.csv", sep = ",", skip = 1 ), byrow = TRUE, ncol = 2 )
# mat <- matrix( scan( "./KoNova Zenith Response.csv", sep = ",", skip = 1 ), byrow = TRUE, ncol = 2 )
forward.zenith <- mat[,1] * pi / 180
forward.zenith.flux <- mat[,2]

these.breaks <- seq( from = 0, to = max( zenith[muon.stdcut & zenith > 0] * 180 / pi ) + 3, by = 1 )
main.label <- paste0( "Histogram of Zenith\n", analysis.indicator )
hist( zenith[muon.stdcut & zenith >0] * 180 / pi, nclass = 200, main = main.label, xlab = "Zenith (degrees)" )
h <- hist( zenith * 180 / pi, nclass = 200, plot = FALSE )
lines( forward.zenith * 180 / pi, forward.zenith.flux * length( zenith[muon.stdcut] ) / sum( forward.zenith.flux ) / ( length( h$counts ) / length( forward.zenith ) ) / 0.9, col = "green" )

pdf( paste0( this.plot.directory, "/", analysis.indicator, "-Zenith.pdf" ) )
hist( zenith[muon.stdcut & zenith >0] * 180 / pi, nclass = 200, main = main.label, xlab = "Zenith (degrees)" )
lines( forward.zenith * 180 / pi, forward.zenith.flux * length( zenith[muon.stdcut] ) / sum( forward.zenith.flux ) / ( length( h$counts ) / length( forward.zenith ) ) / 0.9, col = "green" )
dev.off()


mat <- matrix( scan( "./konova_v1_azimuth_distribution.csv", sep = ",", skip = 1 ), byrow = TRUE, ncol = 2 )
forward.azimuth <- mat[,1] * pi / 180
forward.azimuth.flux <- mat[,2]

these.breaks <- seq( from = 0, to = 360, by = 5 )
main.label <- paste0( "Histogram of Azimuth\n", analysis.indicator )
hist( azimuth[muon.stdcut] * 180 / pi, nclass = 200, main = main.label, xlab = "Azimuth (degrees) East = 0 degrees" )
h <- hist( azimuth[muon.stdcut] * 180 / pi, nclass = 200, plot = FALSE )
#abline( v = c( 0, 90, 180, 270, 360 ), col = "blue" )
# lines( forward.azimuth * 180 / pi, forward.azimuth.flux * length( azimuth[muon.stdcut] ) / sum( forward.azimuth.flux ) / ( length( h$counts ) / length( forward.azimuth ) ), col = "green" )

pdf( paste0( this.plot.directory, "/", analysis.indicator, "-Azimuth.pdf" ) )
hist( azimuth[muon.stdcut] * 180 / pi, nclass = 200, main = main.label, xlab = "Azimuth (degrees) East = 0 degrees" )
# lines( forward.azimuth * 180 / pi, forward.azimuth.flux * length( azimuth[muon.stdcut] ) / sum( forward.azimuth.flux ) / ( length( h$counts ) / length( forward.azimuth ) ), col = "green" )
dev.off()

summary.file <- paste0( this.plot.directory, "/", analysis.indicator, ", Summary.txt" )

this.line <- paste0( "Analysis indicator = ", analysis.indicator )
write( this.line, summary.file, append = FALSE )

this.line <- paste0( "Runs included = ", unique( run.names ) )
write( this.line, summary.file, append = TRUE )

this.line <- paste0( "Number of cycles = ", length( filenames ) )
write( this.line, summary.file, append = TRUE )

this.line <- paste0( "Number of events = ", length( event.hits[muon.stdcut] ) )
write( this.line, summary.file, append = TRUE )

this.line <- paste0( "Run time = ", signif( total.run.time / 60, 3 ), " min = ", signif( total.run.time / 3600, 3 ), " hours = ", signif( total.run.time / ( 24 * 3600 ), 3 ), " days" )
write( this.line, summary.file, append = TRUE )

this.line <- paste0( "Average hit rate = ", signif( length( event.hits[muon.stdcut] ) / ( total.run.time ), 4 ), " +/= ", signif( sqrt( length( event.hits[muon.stdcut] ) ) / ( total.run.time ), 1 ), " Hz" )
write( this.line, summary.file, append = TRUE )

system( paste0( "cat ", shQuote( summary.file ) ) )

# Background rate
# min.hits.per.cluster * ( 3000 )^min.hits.per.cluster * ( 15e-9 )^( min.hits.per.cluster - 1 )

# Look at MPPC distribution lots of junk

# total.hits <- 0
# number.noise <- 0
# avg.coinc.number <- 
#   par( mfrow = c( 2, 2 ) )
# for ( this.layer.number in layer.numbers ) {
#   stdcut <- coinc.layer == this.layer.number
#   stdcut <- stdcut & coinc.number.of.hits >= min.hits
#   avg.coinc.number <- length( coinc.number[stdcut] ) / 64
#   plot( 0, 0, xlim = c( 0, 9 ), ylim = c( 0, 9 ), xlab = "MPPC x", ylab = "MPPC y", main = paste0( "Layer = ", this.layer.number ), type = "n" )
#   for ( i in 1:8 ) {
#     for ( j in 1:8 ) {
#       text( i, j, signif( length( coinc.number[stdcut & coinc.MPPC.x == i & coinc.MPPC.y == j] ) / avg.coinc.number, 3 ), cex = 0.5 )
#       #rect(xleft = i - 0.5, ybottom = j - 0.5, xright = i + 0.5, ytop = j + 0.5, col = gray( length( coinc.number[stdcut & coinc.MPPC.x == i & coinc.MPPC.y == j] ) / avg.coinc.number / 2 ), border = "white", lwd = 1 )
#       total.hits <- total.hits + length( coinc.number[stdcut & coinc.MPPC.x == i & coinc.MPPC.y == j] )
#       if ( i > 1 & i < 8 ) {
#         number.noise <- number.noise + length( coinc.number[stdcut & coinc.MPPC.x == i & coinc.MPPC.y == j] )
#       }
#     }
#   }
# }
# par( mfrow = c( 1, 1 ) )
# number.noise
# total.events <- length( coinc.layer )
# total.hits
# number.noise / total.hits

