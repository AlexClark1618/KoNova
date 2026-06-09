# Assumes compact text format

con <- file( paste0( program.root, data.root, "/", filename, "_coinc.dat" ) ) 
open( con )

line <- readLines( con, n = 1, warn = FALSE )
results.list <- as.character( unlist( strsplit( line, split = "\t" ) ) )
results.list <- results.list[results.list != ""]

total.hits <- 0
while ( length( line ) > 0 ) {
  this.hits <- as.integer( results.list[1] ) + as.integer( results.list[2] )
  total.hits <- total.hits + this.hits
  for ( i in 1:this.hits ) {
    line <- readLines( con, n = 1, warn = FALSE )
  }
  
  line <- readLines( con, n = 1, warn = FALSE )
  if ( length( line ) > 0 ) {
    results.list <- as.character( unlist( strsplit( line, split = "\t" ) ) )
    results.list <- results.list[results.list != ""]
  }
}

close( con )

con <- file( paste0( program.root, data.root, "/", filename, "_coinc.dat" ) ) 
open( con )

line <- readLines( con, n = 1, warn = FALSE )
results.list <- as.character( unlist( strsplit( line, split = "\t" ) ) )
results.list <- results.list[results.list != ""]

petsys.number <- 1
hit.number <- 1
petsys.numbers <- seq( total.hits ) * 0
petsys.time <- seq( total.hits ) * 0 - 1e6
petsys.qdc <- seq( total.hits ) * 0 - 1e6
petsys.ch <- seq( total.hits ) * 0 - 1e6
petsys.number.of.hits <- seq( total.hits ) * 0 - 1e6
ch8187.count <- 0
ch8784.count <- 0
while ( length( line ) > 0 ) {
  number.of.hits <- as.integer( results.list[1] ) + as.integer( results.list[2] )
  these.channels <- NULL
  for ( i in 1:number.of.hits ) {
    petsys.numbers[hit.number] <- petsys.number
    
    line <- readLines( con, n = 1, warn = FALSE )
    results.list <- as.character( unlist( strsplit( line, split = "\t" ) ) )
    results.list <- results.list[results.list != ""]
    
    petsys.time[hit.number] <- as.double( results.list[1] )
    petsys.qdc[hit.number] <- as.double( results.list[2] )
    petsys.ch[hit.number] <- as.integer( results.list[3] )
    these.channels <- c( these.channels, as.integer( results.list[3] ) )
    petsys.number.of.hits[hit.number] <- number.of.hits
    
    hit.number <- hit.number + 1
  }
  if ( min( abs( these.channels - 87 ) ) == 0 ) {
    if ( min( abs( these.channels - 81 ) ) == 0 ) {
      ch8187.count <- ch8187.count + 1
    }
    if ( min( abs( these.channels - 84 ) ) == 0 ) {
      ch8784.count <- ch8784.count + 1
    }
  }
  
  line <- readLines( con, n = 1, warn = FALSE )
  results.list <- as.character( unlist( strsplit( line, split = "\t" ) ) )
  results.list <- results.list[results.list != ""]
  petsys.number <- petsys.number + 1
}

close( con )

# print( paste0( "ch8187.count = ", ch8187.count ) )
# print( paste0( "ch8784.count = ", ch8784.count ) )

hist( ch.2.bar.mat[petsys.ch[petsys.ch < 200] - layer.subtract.this[1] + 1, 1], seq( from = 0.5, to = 64.5, by = 1 ) )

# Permanent change from original file

increasing.time.order <- order( petsys.time )
petsys.numbers <- petsys.numbers[increasing.time.order]
petsys.time <- petsys.time[increasing.time.order]
petsys.qdc <- petsys.qdc[increasing.time.order]
petsys.ch <- petsys.ch[increasing.time.order]
range( diff( petsys.time ) ) 
hist( diff( petsys.time )[diff( petsys.time ) < 10000] )

