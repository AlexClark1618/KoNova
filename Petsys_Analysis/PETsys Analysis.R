# Use for PETsys data
# This file uses PETsys _coinc.dat compact text files
# Initialize R
# Units are cm and ps

setwd( "C:/Users/alexc/Desktop/Petsys_Analysis")

rm( list = ls() )

program.root <- getwd()

data.root <- "/PETsys Data"
stats.root <- "/PETsys Stats"
plot.root <- "/PETsys Plots"

# source( "/users/ifft/documents/Error Analysis/measurement 2.0.r" )

# Run "GPS Select Filenames.R"
source( './PETsys Select Filenames.R' )

# Run "GPS Get Stats Next.R" here in terminals
source( './PETsys Get Stats Next.R' )

# Get collated data

if ( Sys.info()["sysname"] == "Darwin" ) { # MacOS
  
  unique.run.names <- unique( substr( filenames, start = 1, stop = 16 ) )
  
  for ( run.name in unique.run.names ) {
    system( paste0( "cd \'", program.root, "/", stats.root, "\'; cat ", run.name, "*hitstats.dat  | grep -v Event > Collated-hitstats.dat" ) )
    system( paste0( "cd \'", program.root, "/", stats.root, "\'; cat ", run.name, "*eventstats.dat  | grep -v Event > Collated-eventstats.dat" ) )
  }
  filename <- "Collated"
  source( paste0( program.root, "/PETsys Load Hitstats File.R" ) )
  system( paste0( "cd \'", program.root, "/", stats.root, "\'; rm Collated-hitstats.dat" ) )
  source( paste0( program.root, "/PETsys Load Eventstats File.R" ) )
  system( paste0( "cd \'", program.root, "/", stats.root, "\'; rm Collated-eventstats.dat" ) )
  
} else if ( Sys.info()["sysname"] == "Windows" ) { # PC
  
  this.file <- paste0( program.root, stats.root, "/", "Collated", "-hitstats.dat" )
  write( c( "Event Number", "Time (ps)", "QDC", "Channel", "Layer", "Reduced Channel", "MPPC X", "MPPC Y", "Number of hits" ), this.file, ncolumns = 9, append = FALSE, sep = ", " )
  
  this.file <- paste0( program.root, stats.root, "/", "Collated", "-eventstats.dat" )
  write( c( "Event Number", "Hits", "Delta t", "QDC", "Layer 1 Hits", "Layer 1 QDC", "Layer 2 Hits", "Layer 2 QDC", "Layer 3 Hits", "Layer 3 QDC", "Layer 4 Hits", "Layer 4 QDC", "x1 (cm)", "y1 (cm)", "x2 (cm)", "y2 (cm)", "Adjacent Cut", "Edge Cut" ), this.file, ncolumns = 18, append = FALSE, sep = ", " )
  
  for ( file.i in 1:length( filenames ) ) {
    filename <- filenames[file.i]
    
    source( './PETsys Load Hitstats File.R' )
    this.file <- paste0( program.root, stats.root, "/", "Collated", "-hitstats.dat" )
    write( t( cbind( coinc.number, as.character( coinc.time ), coinc.qdc, coinc.channel, coinc.layer, coinc.reduced.channel, coinc.MPPC.x, coinc.MPPC.y, coinc.number.of.hits ) ), this.file, ncolumns = 9, append = TRUE, sep = ", " )
    
    source( './PETsys Load Eventstats File.R' )
    this.file <- paste0( program.root, stats.root, "/", "Collated", "-eventstats.dat" )
    write( t( cbind( event.coinc.number, event.hits, event.delta.t, event.qdc, layer.hits[, 1], layer.qdc[, 1], layer.hits[, 2], layer.qdc[, 2], layer.hits[, 3], layer.qdc[, 3], layer.hits[, 4], layer.qdc[, 4], xydetector.x[, 1], xydetector.y[, 1], xydetector.x[, 2], xydetector.y[, 2], adjacent.cut, edge.cut ) ), this.file, ncolumns = 18, append = TRUE, sep = ", " )
  }
  
  filename <- "Collated"
  source( paste0( program.root, "/PETsys Load Hitstats File.R" ) )
  #system( paste0( "cd \'", program.root, "/", stats.root, "\'; rm Collated-hitstats.dat" ) )
  source( paste0( program.root, "/PETsys Load Eventstats File.R" ) )
  #system( paste0( "cd \'", program.root, "/", stats.root, "\'; rm Collated-eventstats.dat" ) )
  
} else {
  print( "Cannot find operating system!" )
  break
}


# Get cycle stats

cycle.event.numbers <- NULL
cycle.run.times <- NULL
for ( filename in filenames ) {
  # string <- system( paste0( "tail -n 1 \'", program.root, "/", stats.root, "/", filename, "-hitstats.dat\'" ), intern = TRUE )
  
  last.line <- tail( readLines( paste0( program.root, "/", stats.root, "/", filename, "-hitstats.dat" ) ), n = 1 )
  
  this.event.number <- as.integer( unlist( strsplit( last.line, "," ) )[1] )
  cycle.event.numbers <- c( cycle.event.numbers, this.event.number )

  this.run.time <- as.double( unlist( strsplit( last.line, "," ) )[2] ) / 1e12 # in seconds
  cycle.run.times <- c( cycle.run.times, this.run.time )
}
total.run.time <- sum( cycle.run.times )
cycle.rates <- cycle.event.numbers / cycle.run.times


# For MC files
if ( substr( filename, 1, 4 ) == "KNMC" ) {
  coinc.qdc <- coinc.qdc / 160 # to bring the qdc to roughly within the right range
}

# Initialize the cuts.  This assumes all of the initialization files are the same as the first
run.name <- substr( filenames[1], 1, 16 )
this.file <- paste0( program.root, stats.root, "/", run.name, "-ini.R" )
source( this.file )

# Apply final cuts

min.hits <- 4
hits.layer.min <- 1
hits.cut <- event.hits >= min.hits
hits.cut <- hits.cut & layer.hits[,1] >= hits.layer.min & layer.hits[,2] >= hits.layer.min & layer.hits[,3] >= hits.layer.min & layer.hits[,4] >= hits.layer.min


# Working on coordinates

x1 <- xydetector.x[,1]
y1 <- xydetector.y[,1]
x2 <- xydetector.x[,2]
y2 <- xydetector.y[,2]

x1[is.na(x1)] <- 1e6
y1[is.na(y1)] <- 1e6
x2[is.na(x2)] <- 1e6
y2[is.na(y2)] <- 1e6

delta.x <- x2 - x1
delta.y <- y2 - y1

#zenith <- atan( sqrt( delta.x^2 + delta.y^2 ) / delta.z )
zenith <- atan( sqrt( (delta.x * delta.x.stretch + delta.x.add)^2 + (delta.y * delta.y.stretch + delta.y.add)^2 ) / delta.z )
#azimuth <- atan2( delta.y, delta.x )
azimuth <- atan2( delta.y * delta.y.stretch + delta.y.add, delta.x * delta.x.stretch + delta.x.add )
azimuth[azimuth<0] <- 2 * pi + azimuth[azimuth<0] # For Bill

x.max.max <- 107.25
position.stdcut <- x1 > -channel.spacing / 2 & x1 < x.max.max
position.stdcut <- position.stdcut & y1 > -channel.spacing / 2 & y1 < x.max.max
position.stdcut <- position.stdcut & x2 > -channel.spacing / 2 & x2 < x.max.max
position.stdcut <- position.stdcut & y2 > -channel.spacing / 2 & y2 < x.max.max

radius.max <- 51 # in cm
radius.1 <- sqrt( ( x1 - (105.6/2) )^2 + ( y1 - (105.6/2) )^2 )
radius.2 <- sqrt( ( x2 - (105.6/2) )^2 + ( y2 - (105.6/2) )^2 )
radius.cut <- radius.1 < radius.max & radius.2 < radius.max

length( event.hits )
stdcut <- hits.cut
length( event.hits[stdcut] )
stdcut <- hits.cut & adjacent.cut
length( event.hits[stdcut] )
stdcut <- hits.cut & adjacent.cut & radius.cut
length( event.hits[stdcut] )
stdcut <- hits.cut & adjacent.cut & radius.cut & position.stdcut
length( event.hits[stdcut] )
muon.stdcut <- stdcut # save this
muon.stdcut[is.na(muon.stdcut)] <- FALSE
length( event.hits[stdcut] )

source( paste0( program.root, "/PETsys Make Summary Plots.R" ) )

