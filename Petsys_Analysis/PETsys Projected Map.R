
# Get projected XY statistics

depth <- 20 # in m

R <- depth / cos( zenith[muon.stdcut] )

x <- as.double( R * sin( zenith[muon.stdcut] ) * cos( azimuth[muon.stdcut] ) )
y <- as.double( R * sin( zenith[muon.stdcut] ) * sin( azimuth[muon.stdcut] ) )

# Asymmetries

meter.cutoff <- 0
wall.position <- 0
main.descriptor <- "Cliff (3 in)"
this.xlab <- "<- West (m) East ->"
this.ylab <- "<- South (m) North ->"

measurement( length( x[x > 0 & y > 0] ), sqrt( length( x[x > 0 & y > 0] ) ) )
measurement( length( x[x < 0 & y > 0] ), sqrt( length( x[x < 0 & y > 0] ) ) )
measurement( length( x[x < 0 & y < 0] ), sqrt( length( x[x < 0 & y < 0] ) ) )
measurement( length( x[x > 0 & y < 0] ), sqrt( length( x[x > 0 & y < 0] ) ) )

ratio <- measurement( length( x[x > meter.cutoff] ), sqrt( length( x[x > meter.cutoff] ) ) ) / measurement( length( x[x < -meter.cutoff] ), sqrt( length( x[x < -meter.cutoff] ) ) )
ratio
hist.title <- paste0( analysis.indicator, ", Histogram of X\nRatio = ", signif( ratio$val, 4 ), " +/- ", signif( ratio$err, 2 ) )
hist( x[abs(x)<10], main = hist.title, xlab = this.xlab, nclass = 30, xlim = c( -10, 10 ) )
abline( v = wall.position, col = "red" )
abline( v = c( meter.cutoff, -meter.cutoff ), col = "red", lty = 2 )

# jpeg( paste0( plot.root, "/", "HistX-", analysis.indicator, ".jpeg" ) )
# hist( x[abs(x)<10], main = hist.title, xlab = this.xlab, nclass = 30, xlim = c( -10, 10 ) )
# abline( v = wall.position, col = "red" )
# abline( v = c( meter.cutoff, -meter.cutoff ), col = "red", lty = 2 )
# dev.off()

ratio <- measurement( length( y[y > meter.cutoff] ), sqrt( length( y[y > meter.cutoff] ) ) ) / measurement( length( y[y < -meter.cutoff] ), sqrt( length( y[y < -meter.cutoff] ) ) )
ratio
hist.title <- paste0( analysis.indicator, ", Histogram of Y\nRatio = ", signif( ratio$val, 4 ), " +/- ", signif( ratio$err, 2 ) )
hist( y[abs(y)<10], main = hist.title, xlab = this.ylab, nclass = 30, xlim = c( -10, 10 ) )
#abline( v = wall.position, col = "red" )
abline( v = c( meter.cutoff, -meter.cutoff ), col = "red", lty = 2 )

# jpeg( paste0( plot.root, "/", "HistY-", analysis.indicator, ".jpeg" ) )
# hist( y[abs(y)<10], main = hist.title, xlab = this.ylab, nclass = 30, xlim = c( -10, 10 ) )
# #abline( v = wall.position, col = "red" )
# abline( v = c( meter.cutoff, -meter.cutoff ), col = "red", lty = 2 )
# dev.off()


# Make heat map

L <- 40 # side length of the map in m
number.pixels.per.side <- 81
delta.L <- L / number.pixels.per.side

x.pixel <- seq( from = -L/2 + delta.L/2, to = L/2 - delta.L/2, length.out = number.pixels.per.side )
y.pixel <- x.pixel

data <- expand.grid(X=x.pixel, Y=y.pixel)

numbers <- data$X * 0
for ( i in seq( data$X ) ) {
  stdcut <- x >= data$X[i] - delta.L/2 & x < data$X[i] + delta.L/2
  stdcut <- stdcut & y >= data$Y[i] - delta.L/2 & y < data$Y[i] + delta.L/2
  numbers[i] <- length( x[stdcut] )
}
data$Z <- numbers

main.label <- paste0( analysis.indicator, "\n", main.descriptor )
image.plot( unique( data$X ), unique( data$Y ), matrix( data$Z, ncol = length( unique( data$X ) ), nrow = length( unique( data$Y ) ) ), col = heat.colors(100), xlab = this.xlab, ylab = this.ylab, main = main.label )

# jpeg( paste0( plot.root, "/", "HeatMap-", analysis.indicator, ".jpeg" ) )
# image.plot( unique( data$X ), unique( data$Y ), matrix( data$Z, ncol = length( unique( data$X ) ), nrow = length( unique( data$Y ) ) ), col = heat.colors(100), xlab = this.xlab, ylab = this.ylab, main = main.label )
# dev.off()

plot( data$X, data$Y, type = "n", xlab = this.xlab, ylab = this.ylab, main = paste0( analysis.indicator, "\n", main.descriptor, " Counts per Pixel in ", round( total.run.time/3600, 1 ), " hrs" ) )
abline( v = wall.position, col = "red" )
abline( v = meter.cutoff, col = "red", lty = 2 )
abline( v = -meter.cutoff, col = "red", lty = 2 )
text( data$X, data$Y, paste( data$Z ), cex = 0.25 )

# jpeg( paste0( plot.root, "/", "NumbersMap-", analysis.indicator, ".jpeg" ) )
# plot( data$X, data$Y, type = "n", xlab = this.xlab, ylab = this.ylab, main = paste0( analysis.indicator, "\n", main.descriptor, " Counts per Pixel in ", round( total.run.time/3600, 1 ), " hrs" ) )
# abline( v = wall.position, col = "red" )
# abline( v = meter.cutoff, col = "red", lty = 2 )
# abline( v = -meter.cutoff, col = "red", lty = 2 )
# text( data$X, data$Y, paste( data$Z ) )
# dev.off()


hist( coinc.qdc[coinc.layer == 1 & coinc.reduced.channel == 1], nclass = 100 )
length( coinc.qdc[coinc.layer == 1 & coinc.reduced.channel == 1] )
median( coinc.qdc[coinc.layer == 1 & coinc.reduced.channel == 1] )
hist( coinc.qdc[coinc.layer == 1 & coinc.reduced.channel == 2], nclass = 100)
length( coinc.qdc[coinc.layer == 1 & coinc.reduced.channel == 2] )
median( coinc.qdc[coinc.layer == 1 & coinc.reduced.channel == 2] )
hist( coinc.qdc[coinc.layer == 1 & coinc.reduced.channel == 3], nclass = 100)
length( coinc.qdc[coinc.layer == 1 & coinc.reduced.channel == 3] )
median( coinc.qdc[coinc.layer == 1 & coinc.reduced.channel == 3] )


