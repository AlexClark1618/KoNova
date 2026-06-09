# Use for PETsys data
# This file has code for analysis of _single.dat files
# Initialize R

setwd( "/Users/ifft/Desktop/Muography/PETsys Analysis" )

rm( list = ls() )

program.root <- "/Users/ifft/Desktop/Muography/PETsys Analysis"
data.root <- "PETsys Data"
plot.root <- "PETsys Plots"
stats.root <- "PETsys Stats"

source( paste0( program.root, "/PETsys Initialize.R" ) )

# Load one single file

filename <- "KNVA-20251115-01-00001"

source( paste0( program.root, "/PETsys Load Single Data.R" ) )

# Get coincidence file

source( paste0( program.root, "/PETsys Cluster.R" ) )
source( paste0( program.root, "/PETsys Get Coincidence File.R" ) )

# Load one coincidence file

source( paste0( program.root, "/PETsys Load Coinc File.R" ) )

par( mfrow = c( 2, 2 ) )
for ( this.layer.number in 1:4 ) {
  stdcut <- coinc.layer == this.layer.number
  plot( 0, 0, xlim = c( 0, 9 ), ylim = c( 0, 9 ), xlab = "MPPC x", ylab = "MPPC y", main = paste0( "Layer = ", this.layer.number ), type = "n" )
  for ( i in 1:8 ) {
    for ( j in 1:8 ) {
      text( i, j, length( coinc.number[stdcut & coinc.MPPC.x == i & coinc.MPPC.y == j] ) )
    }
  }
}
par( mfrow = c( 1, 1 ) )

# Get event statistics

percent.print <- 0.01
event.hits <- unique( coinc.number ) * 0 + 1e6 # overall number of hits, 4 or greater because of above
event.delta.t <- unique( coinc.number ) * 0 + 1e6 # max delta t for all hits
event.qdc <- unique( coinc.number ) * 0 + 1e6 # total qdc
layer.hits <- matrix( 1e6, nrow = length( unique( coinc.number ) ), ncol = 4 )
layer.qdc <- matrix( 1e6, nrow = length( unique( coinc.number ) ), ncol = 4 )
xydetector.x <- matrix( 1e6, nrow = length( unique( coinc.number ) ), ncol = 2 )
xydetector.y <- matrix( 1e6, nrow = length( unique( coinc.number ) ), ncol = 2 )
for ( this.coinc.number in unique( coinc.number ) ) {
  if ( this.coinc.number / length( coinc.number ) > percent.print ) {
    print( paste0( percent.print * 100, "% done" ) )
    percent.print <- 2 * percent.print
  }
  coinc.stdcut <- coinc.number == this.coinc.number
  
  # event statistics
  event.hits[this.coinc.number] <- length( coinc.time[coinc.stdcut] )
  event.delta.t[this.coinc.number] <- diff( range( coinc.time[coinc.stdcut] ) )
  event.qdc[this.coinc.number] <- sum( coinc.qdc[coinc.stdcut] )
  
  # layer statistics
  for ( layer in 1:4 ) {
    layer.stdcut <- coinc.layer == layer
    layer.hits[this.coinc.number, layer] <- length( coinc.time[coinc.stdcut & layer.stdcut] )
    layer.qdc[this.coinc.number, layer] <- sum( coinc.qdc[coinc.stdcut & layer.stdcut] )
  }
  
  # xydetector statistics
  layer.stdcut <- coinc.layer == 1 # x analysis on top layer
  xydetector.x[this.coinc.number, 1] <- weighted.mean( coinc.MPPC.y[coinc.stdcut & layer.stdcut], w = coinc.qdc[coinc.stdcut & layer.stdcut] ) * channel.spacing
  layer.stdcut <- coinc.layer == 2 # y analysis on top layer
  xydetector.y[this.coinc.number, 1] <- weighted.mean( coinc.MPPC.y[coinc.stdcut & layer.stdcut], w = coinc.qdc[coinc.stdcut & layer.stdcut] ) * channel.spacing
  layer.stdcut <- coinc.layer == 3 # x analysis on bottom layer
  xydetector.x[this.coinc.number, 2] <- weighted.mean( coinc.MPPC.y[coinc.stdcut & layer.stdcut], w = coinc.qdc[coinc.stdcut & layer.stdcut] ) * channel.spacing
  layer.stdcut <- coinc.layer == 4 # y analysis on bottom layer
  xydetector.y[this.coinc.number, 2] <- weighted.mean( coinc.MPPC.y[coinc.stdcut & layer.stdcut], w = coinc.qdc[coinc.stdcut & layer.stdcut] ) * channel.spacing
}

hist( layer.qdc[,1] )
hist( layer.qdc[,2] )
hist( layer.qdc[,3] )
hist( layer.qdc[,4] )

# Graphs for notebook

hist( event.hits, breaks = seq( from = 3.5, to = 24.5, by = 1 ) )

hist( event.delta.t, xlab = "Delta t (ps)" )

hist( event.qdc )

layer.hits

plot( xydetector.x[,1], xydetector.y[,1] )
plot( xydetector.x[,2], xydetector.y[,2] )

# Cuts one per event

hits.layer.cut <- 2 # minimum hits per layer
hits.cut <- layer.hits[,1] >= hits.layer.cut & layer.hits[,2] >= hits.layer.cut & layer.hits[,3] >= hits.layer.cut & layer.hits[,4] >= hits.layer.cut

adjacent.cut <- seq( event.hits ) > -1 # All TRUE to start.  False if on any layer there are non-adjacent hits
for ( i in seq( event.hits ) ) {
  coinc.stdcut <- coinc.number == i
  for ( this.layer.number in 1:4 ) {
    layer.stdcut <- coinc.layer == this.layer.number # x analysis on top layer
    if( length( coinc.MPPC.y[coinc.stdcut & layer.stdcut] ) >= 2 ) {
      if( max( diff( sort( coinc.MPPC.y[coinc.stdcut & layer.stdcut] ) ) ) >= 2 ) {
        print( paste0( "This good event number = ", i ) )
        print( paste0( "Layer #", this.layer.number ) )
        print( coinc.MPPC.y[coinc.stdcut & layer.stdcut] )
        adjacent.cut[i] <- FALSE
      }
    }
  }
}

length( event.hits )
stdcut <- hits.cut
length( event.hits[stdcut] )
stdcut <- hits.cut & adjacent.cut
length( event.hits[stdcut] )


# Graphs for notebook

length( event.hits[stdcut] )

layer.hits[stdcut,]

hist( event.hits[stdcut], breaks = seq( from = 3.5, to = 10.5, by = 1 ) )

hist( event.delta.t[stdcut], xlab = "Event Delta t (ps)" )

hist( event.qdc[stdcut], xlab = "QDC", main = paste0( "QDC distribution\nmean = ", round( mean( event.qdc[stdcut] ), 1 ) ) )


plot( xydetector.x[stdcut,1], xydetector.y[stdcut,1], main = "Top Layer", xlab = "x (cm)", ylab = "y (cm)" )
plot( xydetector.x[stdcut,2], xydetector.y[stdcut,2], main = "Bottom Layer", xlab = "x (cm)", ylab = "y (cm)" )


plot( xydetector.x[stdcut,1], xydetector.y[stdcut,1], main = "Top Layer", xlab = "x (cm)", ylab = "y (cm)", type = "n" )
text( xydetector.x[stdcut,1], xydetector.y[stdcut,1], seq( event.hits[stdcut] ) )
plot( xydetector.x[stdcut,2], xydetector.y[stdcut,2], main = "Bottom Layer", xlab = "x (cm)", ylab = "y (cm)", type = "n" )
text( xydetector.x[stdcut,2], xydetector.y[stdcut,2], seq( event.hits[stdcut] ) )

# layer analysis with muons

for ( i in 1:4 ) {
  print( paste0( "Layer #", i, " QDC = ", round( mean( layer.qdc[stdcut,i] ), 2 ), " +/- ", round( sqrt( var( layer.qdc[stdcut,i] ) / length( layer.qdc[,i] ) ), 2 ) ) )
}



# Using the PETsys generated coinc file

source( paste0( program.root, "/PETsys Load PETsys Coinc Data.R" ) )

# Permanent change from original file

increasing.time.order <- order( coinc.time )
coinc.time <- coinc.time[increasing.time.order]
coinc.qdc <- coinc.qdc[increasing.time.order]
coinc.ch <- coinc.ch[increasing.time.order]
range( diff( coinc.time ) ) 
hist( diff( coinc.time )[diff( coinc.time ) < 10000] )

source( paste0( program.root, "/PETsys Cluster Coinc Data.R" ) )









hist( event.delta.t[stdcut] )
range( event.delta.t[stdcut] )
length( event.delta.t[stdcut] )

hist( event.qdc.0, nclass = 100 )
abline( v = qdc.layer.cut, col = "red" )
length( event.qdc.0[event.qdc.0 > qdc.layer.cut] )

hist( event.qdc.1, nclass = 100 )
abline( v = qdc.layer.cut, col = "red" )
length( event.qdc.0[event.qdc.0 > qdc.layer.cut] )

plot( event.qdc.0, event.qdc.1, xlim = c( 0, 20 ), ylim = c( 0, 20 ), pch = 20, cex = 0.1 )
abline( v = qdc.layer.cut, col = "red" )
abline( h = qdc.layer.cut, col = "red" )
abline( a = qdc.layer.cut * 2, b = -1, col = "red" )

stdcut <- event.qdc.total > qdc.layer.cut * 2
main.label <- paste0( "Weighted Muon Locations, XY-Detector #1\n", filename )
plot( event.x[stdcut], event.y[stdcut], xlab = "x (mm)", ylab = "y (mm)", main = main.label )
abline( v = 0:8 * channel.spacing, col = "gray" )
abline( h = 0:8 * channel.spacing, col = "gray" )
length( event.qdc.total[stdcut])

stdcut <- event.hits.0 >= hits.layer.cut & event.hits.1 >= hits.layer.cut
stdcut <- stdcut & event.qdc.total > qdc.layer.cut * 2
stdcut <- stdcut & event.delta.t < delta.t.cut
main.label <- paste0( "Weighted Muon Locations, XY-Detector #1\n", filename, "\n2+2, energy and timing cuts, N = ", length( event.hits.total[stdcut] ) )
plot( event.x[stdcut], event.y[stdcut], xlab = "x (mm)", ylab = "y (mm)", main = main.label )
abline( v = 0:8 * channel.spacing, col = "gray" )
abline( h = 0:8 * channel.spacing, col = "gray" )
pdf( file = paste0( "XY plot - ", filename, ".pdf" ) )
plot( event.x[stdcut], event.y[stdcut], xlab = "x (mm)", ylab = "y (mm)", main = main.label )
abline( v = 0:8 * channel.spacing, col = "gray" )
abline( h = 0:8 * channel.spacing, col = "gray" )
dev.off()




length( coinc.qdc[coinc.channel < 200] )
length( coinc.qdc[coinc.channel > 200] )

# hist( coinc.qdc[coinc.channel < 200], breaks = seq( from = -15, to = 20, by = 0.05 ) )
# hist( coinc.qdc[coinc.channel > 200], breaks = seq( from = -15, to = 20, by = 0.05 ) )

stdcut <- coinc.channel < 200
hist( coinc.channel[stdcut], breaks = seq( from = min( coinc.channel[stdcut] ), to = max( coinc.channel[stdcut] ), by = 1 ) )
range( coinc.channel[stdcut] )

stdcut <- coinc.channel > 200 & coinc.channel < 400
hist( coinc.channel[stdcut], breaks = seq( from = min( coinc.channel[stdcut] ), to = max( coinc.channel[stdcut] ), by = 1 ) )
range( coinc.channel[stdcut] )

stdcut <- coinc.channel > 400 & coinc.channel < 700
hist( coinc.channel[stdcut], breaks = seq( from = min( coinc.channel[stdcut] ), to = max( coinc.channel[stdcut] ), by = 1 ) )
range( coinc.channel[stdcut] )

stdcut <- coinc.channel > 700
hist( coinc.channel[stdcut], breaks = seq( from = min( coinc.channel[stdcut] ), to = max( coinc.channel[stdcut] ), by = 1 ) )
range( coinc.channel[stdcut] )

plot( coinc.channel[coinc.channel < 200], coinc.qdc[coinc.channel < 200], pch = 20, cex = 0.1 )
plot( coinc.channel[coinc.channel > 200], coinc.qdc[coinc.channel > 200], pch = 20, cex = 0.1 )

hist( diff( coinc.time )[ diff( coinc.time ) < 100000 ], nclass = 200 )
# seems to be a +/- 10 or 11 ns coincidence window

coinc.bts.label <- array( "", dim = length( coinc.channel ) )
coinc.layer <- seq( coinc.channel ) * 0 - 1e6
coinc.SiPMx <- seq( coinc.channel ) * 0 - 1e6
coinc.SiPMy <- seq( coinc.channel ) * 0 - 1e6
coinc.event.number <- seq( coinc.channel ) * 0
for ( i in seq( coinc.channel ) ) {
  if ( i == 1 ) {
    this.event.number <- 1
  } else {
    if ( abs( coinc.time[i] - coinc.time[i-1] ) > gap ) {
      this.event.number <- this.event.number + 1
    }
  }
  layer.number <- floor( coinc.channel[i] / 64 )
  coinc.layer[i] <- ( layer.number - 1 ) / 2 / 2 # first layer is 0, C style, check for more than 2 layers!
  reduced.channel <- ( coinc.channel[i] / 64 - floor( coinc.channel[i] / 64 ) ) * 64
  coinc.bts.label[i] <- mapping$bts.label[reduced.channel == mapping$ChNumber]
  coinc.SiPMx[i] <- mapping$SiPMx[reduced.channel == mapping$ChNumber] - 1 # 0 to 7
  coinc.SiPMy[i] <- mapping$SiPMy[reduced.channel == mapping$ChNumber] - 1 # 0 to 7
  coinc.event.number[i] <- this.event.number
}
length( coinc.event.number )
range( coinc.event.number )
hist( coinc.event.number, seq( from = 0.5, to = max( coinc.event.number ) + 0.5, by = 1 ) )
h <- hist( coinc.event.number, seq( from = 0.5, to = max( coinc.event.number ) + 0.5, by = 1 ), plot = FALSE )
cbind( h$mids, h$counts )

length( h$mids[h$counts > 6] )
for ( this.coinc.event.number in h$mids[h$counts > 6] ) {
  stdcut <- coinc.event.number == this.coinc.event.number
  print( paste0( "Number of hits = ", length( unique( coinc.channel[stdcut] ) ) ) )
  print( paste0( "Delta t = ", diff( range( coinc.time[stdcut] ) ) / 1000, " ns" ) )
  # print( unique( coinc.channel[stdcut] ) )
  for ( this.layer in 0:3 ) {
    stdcut <- coinc.event.number == this.coinc.event.number & coinc.layer == this.layer
    print( paste0( "Layer ", this.layer ) )
    print( unique( coinc.bts.label[stdcut] ) )
  }
}




# 08 started 4:11 PM

# Study other cuts

hist( event.qdc.total[stdcut], nclass = 100 )
abline( v = qdc.layer.cut * 2, col = "red" )

hist( event.delta.t[stdcut], nclass = 100 )




# Load one file

filename <- "KNVA-20251022-02-00001"

source( paste0( program.root, "/PETsys Load Data.R" ) )

sort( unique( ch ) )
plot( sort( unique( ch ) ), type = "n" )
lines( sort( unique( ch ) ) )
hist( ch, nclass = 200 )
hist( ch[ch < 200], breaks = seq( from = 63, to = 128, by = 1 ) )
hist( ch[ch > 200], breaks = seq( from = 320, to = 384, by = 1 ) )

length( ch[ch < 200] )
length( ch[ch > 200] )
length( ch[ch > 200] ) / length( ch[ch < 200] )

hist( QDC[ch < 200], breaks = seq( from = -1.5, to = 10, by = 0.05 ) )
hist( QDC[ch > 200], breaks = seq( from = -1.5, to = 10, by = 0.05 ) )

plot( ch, QDC, type = "n" )
points( ch, QDC, pch = 20, cex = 0.1 )


# Look at one file

this.event.number <- 0

this.event.number <- this.event.number + 1
stdcut <- ms.number < max.ms.number & ns.number < max.ns.number & RF == 0 & event.number == this.event.number & ID != 112
these.ns.times <- ns.time[stdcut] - min( ns.time[stdcut] )
ID[stdcut]
sort( these.ns.times )
cbind( ID[stdcut], these.ns.times )

hist( these.ns.times, breaks = seq( from = min( these.ns.times ) - 200, to = max( these.ns.times ) + 200, by = 200 ) )

plot( these.ns.times )


# Run "GPS Select Filenames.R"


# Analyze multiple files

max.ms.number <- 7 * 24 * 3600 * 1000
max.ns.number <- 1e6
exclude.these.macs <- c( 1e6 )
cycle.times <- NULL
cycle.numbers <- NULL
coincident.dectors.hits <- NULL
IDs.hit <- NULL
for ( filename in filenames ) {
  print( paste0( "Analyzing ", filename ) )
  source( paste0( program.root, "/GPS Load Data.R" ) )
  
  number.stdcut <- ms.number < max.ms.number & ns.number < max.ns.number
  mac.stdcut <- ms.number * 0 > -1 # all TRUE
  for ( exclude.this.mac in exclude.these.macs ) {
    mac.stdcut <- mac.stdcut & ID != exclude.this.mac
  }

  # diff.ms.number <- diff( ms.number )
  # hist( abs( diff.ms.number )[abs( diff.ms.number ) < 10 * 1000], nclass = 1000 )
  # length( abs( diff.ms.number )[abs( diff.ms.number ) < 10 * 1000] )
  # sort( abs( diff.ms.number )[abs( diff.ms.number ) < 10 * 1000] )
  # ten.second.stdcut <- c( TRUE, abs( diff.ms.number ) < 10 * 1000 )
  # ten.second.stdcut <- ten.second.stdcut & c( abs( diff.ms.number ) < 10 * 1000, TRUE ) # cut either up or down
  # 
  # hist( ms.number[ten.second.stdcut & number.stdcut], nclass = 1000 )
  
  diff.mean.cut <- 1e7 # runs not longer than about 5 hours
  diff.mean <- ms.number - median( ms.number[number.stdcut & mac.stdcut] )
  hist( diff.mean[number.stdcut & mac.stdcut & abs( diff.mean ) < diff.mean.cut], nclass = 100 )
  abline( v = c( -diff.mean.cut, diff.mean.cut ), col = "red" )
  
  this.cycle.number <- length( ms.number[number.stdcut & mac.stdcut] )
  cycle.numbers <- c( cycle.numbers, this.cycle.number )
  
  this.cycle.time <- diff( range( diff.mean[number.stdcut & mac.stdcut & abs( diff.mean ) < diff.mean.cut] ) ) / 1000 # in s
  cycle.times <- c( cycle.times, this.cycle.time ) # in s
  
  this.file <- paste0( program.root, stats.root, "/", filename, "_run_stats.txt" )
  write( c( "Number of Events", "Cycle Time (s)", "Cycle Rate (Hz)" ), this.file, ncolumns = 3, sep = ", ", append = FALSE )
  write( c( this.cycle.number, this.cycle.time, this.cycle.number / this.cycle.time ), this.file, ncolumns = 3, sep = ", ", append = TRUE )

  # number.stdcut <- ms.number < max.ms.number & ns.number < max.ns.number
  # diff.mean <- ms.number - mean( ms.number[number.stdcut] )
  # hist( diff.mean[number.stdcut & abs( diff.mean ) < 1e7], nclass = 1000 ) # no longer than 28 hours
  # cycle.time <- diff( range( diff.mean[abs( diff.mean ) < 1e7] ) ) / 1000 # in s
  # print( length( ms.number ) )
  # print( length( ms.number[number.stdcut] ) )
  # print( cycle.time )
  
  # LR analysis from data taken with the borehole detector
  
  # length( ms.number )
  # length( ms.number[number.stdcut] )
  # length( ms.number ) / length( ms.number[ms.number < max.ms.number  & ns.number < max.ns.number] )
  # length( ms.number[number.stdcut & RF == 0] )
  print( paste0( "Number of borehole requests = ", length( ms.number[number.stdcut & RF == 0 & ID == 48] ) ) )
  print( paste0( "This cycle time = ", this.cycle.time, " s" ) )
  print( paste0( "Borehole rate = ", signif( length( ms.number[number.stdcut & RF == 0 & ID == 48] ) / this.cycle.time, 4 ), " Hz" ) )
  live.time <- length( ms.number[number.stdcut & RF == 0 & ID == 48] ) * 0.002
  live.time
  
  this.file <- paste0( program.root, plot.root, "/", filename, "_LR_plots.pdf" )
  pdf( file = this.file )
  for ( this.ID in sort( unique( ID[ID != 48] ) ) ) {
    print( paste0( "ID = ", this.ID ) )
    print( paste0( "Overall rate ", round( ( length( ms.number[number.stdcut & RF == 0 & ID == this.ID] ) ) / live.time ), " Hz" ) )
    L.rate <- ( length( ms.number[number.stdcut & RF == 0 & Ch == 0 & ID == this.ID] ) ) / live.time
    print( paste0( "L rate ", round( L.rate ), " Hz" ) )
    R.rate <- ( length( ms.number[number.stdcut & RF == 0 & Ch == 1 & ID == this.ID] ) ) / live.time
    print( paste0( "R rate ", round( R.rate ), " Hz" ) )
    print( "" )
    
    ns.time.L <- ns.time[number.stdcut & RF == 0 & Ch == 0 & ID == this.ID]
    ns.time.R <- ns.time[number.stdcut & RF == 0 & Ch == 1 & ID == this.ID]
    
    
    these.indexs <- seq( ns.time )[number.stdcut & RF == 0 & Ch == 0 & ID == this.ID]
    ToT.L <- ns.time[these.indexs+1] - ns.time[these.indexs] # in ns
    if ( length( ToT.L[abs( ToT.L ) < 500] ) > 0 ) hist( ToT.L[abs( ToT.L ) < 500], xlab = "ToT in ns", main = paste0( "Left (", this.ID, ") ToT\nL rate = ", round( L.rate ), " Hz" ) )
    these.indexs <- seq( ns.time )[number.stdcut & RF == 0 & Ch == 1 & ID == this.ID]
    ToT.R <- ns.time[these.indexs+1] - ns.time[these.indexs] # in ns
    if ( length( ToT.R[abs( ToT.R ) < 500] ) > 0 ) hist( ToT.R[abs( ToT.R ) < 500], xlab = "ToT in ns", main = paste0( "Right (", this.ID, ") ToT\nR rate = ", round( R.rate ), " Hz" ) )
    
    if ( length( ns.time.L ) > 0 ) plot( ( ns.time.L - min( ns.time.L ) ) / 1e9, ylab = "Event.time (s)", pch = 20, cex = 0.1, main = paste0( "Event times\n", this.ID, "L" ) )
    if ( length( ns.time.L ) > 1 ) plot( diff( ns.time.L ) / 1e9, ylab = "Difference in time (s)", pch = 20, cex = 0.1, main = paste0( "Diff times\n", this.ID, "L" ) )
    if ( length( ns.time.L ) > 1 ) hist( diff( ns.time.L ) / 1e9, main = paste0( "Diff times\n", this.ID, "L" ), nclass = 100 )
    if ( length( ns.time.R ) > 0 ) plot( ( ns.time.R - min( ns.time.R ) ) / 1e9, ylab = "Event.time (s)", pch = 20, cex = 0.1, main = paste0( "Event times\n", this.ID, "R" ) )
    if ( length( ns.time.R ) > 1 ) plot( diff( ns.time.R ) / 1e9, ylab = "Difference in time (s)", pch = 20, cex = 0.1, main = paste0( "Diff times\n", this.ID, "R" ) )
    if ( length( ns.time.R ) > 1 ) hist( diff( ns.time.R ) / 1e9, main = paste0( "Diff times\n", this.ID, "R" ), nclass = 100 )
    
    time.diff <- NULL
    
    min.ToT <- NULL
    for ( i in 1:length( ns.time.L ) ) {
      this.ns.time.L <- ns.time.L[i]
      for ( j in 1:length( ns.time.R ) ) {
        this.ns.time.R <- ns.time.R[j]
        if ( abs( this.ns.time.L - this.ns.time.R ) < 100 / 2 ) {
          # print( i )
          # print( paste0( "Found one!  L - R = ", ( this.ns.time.L - this.ns.time.R ), " ns" ) )
          time.diff <- c( time.diff, this.ns.time.L - this.ns.time.R ) # L - R
          # print( ToT.L[i] )
          # print( ToT.R[j] )
          ToTs <- c( ToT.L[i], ToT.R[j] )
          if ( length( min( ToTs ) ) == 1 ) {
            min.ToT <- c( min.ToT, min( ToTs ) )
          } else {
            min.ToT <- c( min.Tot, ToTs[1] )
          }
        }
      }
    }
    
    # Coincidence rate
    length( time.diff )
    live.time
    LR.rate <- length( time.diff ) / live.time
    
    stdcut <- abs( time.diff ) < 1000 & min.ToT > 0
    if ( length( time.diff[stdcut] ) > 0 ) {
      hist( time.diff[stdcut], nclass = 100, xlab = "L - R (ns)", main = paste0( "Histogram of L - R (", this.ID, ") times within 100 ns\nmean = ", signif( mean( time.diff[stdcut] ), 3 ), ", RMS = ", signif( sqrt( var( time.diff[stdcut] ) ), 3 ), ", LR Rate = ", round( LR.rate ), " Hz" ) ) 
      abline( v = c( -10, 10 ), col = "red" )
    }
    
    # plot( min.ToT[stdcut], time.diff[stdcut], pch = 20, cex = 0.5, main = "Time diff vs ToT", xlab = "Minimum ToT (ns)", ylab = "Time Diff (ns)" )
    # 
    # diff.sorted.times <- diff( sort( ns.time[ms.number < max.ms.number & RF == 0 & ( ID == 16 | ID == 200 )] ) )
    # sort( diff.sorted.times )
    # hist( diff.sorted.times[diff.sorted.times < 200e-6], nclass = 50 )
    
  }
  dev.off()
  
  # Here are Alex's notes
  # 
  # Req Code; ID; RF; Cal; Ch, W#; t_ow mil; t_ow submil; Event #
  # You can ignore req code
  # ID is the last byte of the esp mac address (20 for the borehole in this case and 200 for the array)
  # RF (Rise =0, Fall = 1)
  # Ignore cal
  # Ch is the gps channel (the encoding of the borehole is being sent to channel 1
  # Week number
  # Time of week in milliseconds
  # Time of week nanoseconds
  # Event # is something I add on (all proper borehole signals should have array data with a matching event #)
  # Note: The time stamp is not UNIX its from Jan 6 1980 (if I remember correctly)
  
  # Get borehole request event.numbers
  
  stdcut <- ID == 48 & number.stdcut & mac.stdcut & RF == 0
  
  unique.bh.request.numbers <- unique( event.number[stdcut] )
  length( unique.bh.request.numbers )
  
  # Find time diffs within 1 micros
  
  # min.time.diff <- unique.bh.request.numbers * 0 + 1e6
  print( filename )
  percent.print <- 0.01
  for ( i in 1:length( unique.bh.request.numbers ) ) {
    if ( i / length( unique.bh.request.numbers ) > percent.print ) { 
      print( paste0( percent.print * 100, "% done" ) )
      percent.print <- percent.print + percent.print
    }
    this.event.number <- unique.bh.request.numbers[i]
    stdcut <- ID != 48 & event.number == this.event.number & RF == 0 # AS
    as.ms.times <- ms.number[stdcut]
    as.ns.times <- ns.number[stdcut]
    as.IDs <- ID[stdcut]
    stdcut <- ID == 48 & event.number == this.event.number & RF == 0 # BH
    bh.ms.time <- ms.number[stdcut]
    bh.ns.time <- ns.number[stdcut]
    if ( length( as.ns.times ) > 0 ) {
      time.diffs <- ( as.ms.times - bh.ms.time ) * 1e6 + ( as.ns.times - bh.ns.time ) # in ns
      time.diff.stdcut <- abs( time.diffs ) < 500 # within a micros
      if ( length( time.diffs[time.diff.stdcut] ) > 0 ) {
        print( i )
        print( paste0( "Found one! Event Number = ", this.event.number, ", Length = ", length( time.diffs[time.diff.stdcut] ) + 1, ", Detectors hit = ", length( unique( as.IDs[time.diff.stdcut] ) ) + 1 ) )
        print( paste0( "IDs = ", c( 48, as.IDs[time.diff.stdcut] ) ) )
        print( paste0( "time diffs (AS - BH) = ", time.diffs[time.diff.stdcut], " ns" ) )
        print( "" )
        
        coincident.dectors.hits <- c( coincident.dectors.hits, length( unique( as.IDs[time.diff.stdcut] ) ) + 1)
        IDs.hit <- c( IDs.hit, as.IDs[time.diff.stdcut] )
      }
    }
  }
  
  # 
  # hist( min.time.diff[ min.time.diff > -50 & min.time.diff < 50], nclass = 50, xlab = "Min time difference (AS - BH) in ns", main = "Minimum Time Difference - 100 NS range" )
  # hist( min.time.diff[ min.time.diff > -500 & min.time.diff < 500], nclass = 50, xlab = "Min time difference (AS - BH) in ns", main = "Minimum Time Difference - 1 micros range" )
  # hist( min.time.diff[ min.time.diff > -5000 & min.time.diff < 5000], nclass = 50, xlab = "Min time difference (AS - BH) in ns", main = "Minimum Time Difference - 10 micros range" )
  # hist( min.time.diff[ min.time.diff > -50000 & min.time.diff < 50000], xlab = "Min time difference (AS - BH) in ns", main = "Minimum Time Difference - 100 micros range" )
  # hist( min.time.diff[ min.time.diff > -500000 & min.time.diff < 500000], xlab = "Min time difference (AS - BH) in ns", main = "Minimum Time Difference - 1 ms range" )
  # hist( min.time.diff[ min.time.diff > -100e5 & min.time.diff < 100e5], xlab = "Min time difference (AS - BH) in ns", main = "Minimum Time Difference - 20 ms range" )
  
}

cycle.rates <- cycle.numbers / cycle.times # Hz

plot( cycle.rates, ylab = "Rate (Hz) After Cuts", main = analysis.indicator )

plot( cycle.numbers, ylab = "Number of Selected Events", main = analysis.indicator )
abline( h = 0 )

plot( cycle.times / 3600, ylab = "Cycle Time (hrs)", main = analysis.indicator )
abline( h = 0 )

cbind( filenames, cycle.numbers, cycle.times, cycle.rates )

sum( cycle.times ) / 3600

hist( coincident.dectors.hits, breaks = seq( from = 1, to = 11, by = 1 ) - 0.5, xlab = "Number of detectors (including borehole) hit within 1 micros", main = paste0( "GPS DAQ ", analysis.indicator, "\n>4 = ", length( coincident.dectors.hits[coincident.dectors.hits>4] ), " in ", round( sum( cycle.times ) / 3600 ), " hours" ) )

hist( IDs.hit )

