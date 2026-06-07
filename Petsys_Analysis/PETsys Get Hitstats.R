# Get and save coincidences.  4 or more hardwired

print( "Getting Hit Stats" )

event.number.hits <- diff( cluster.start.indices )
length( event.number.hits )
length( cluster.start.times ) # lost 2
length( event.number.hits[event.number.hits > 3] )
hist( event.number.hits[event.number.hits > 3] )

stdcut <- event.number.hits >= min.hits.per.cluster # forget anything with less than min.hits.per.cluster
good.indices <- seq( event.number.hits )[stdcut]
print( paste0( length( good.indices ), " coincidences" ) )

# Remove duplicates and write out a coincidence file

this.file <- paste0( program.root, stats.root, "/", filename, "-hitstats.dat" )
write( c( "Event Number", "Time (ps)", "QDC", "Channel", "Layer", "Reduced Channel", "MPPC X", "MPPC Y" ), this.file, ncolumns = 8, append = FALSE, sep = ", " )

event.number <- 1
old.petsys.time <- petsys.time[good.indices[1]]
old.petsys.qdc <- petsys.qdc[good.indices[1]]
old.petsys.ch <- petsys.ch[good.indices[1]]

percent.print <- 0.01
for ( i in good.indices ) {
  if ( i / max( good.indices ) > percent.print ) {
    print( paste0( percent.print * 100, "% done" ) )
    percent.print <- 2 * percent.print
  }
  for ( j in cluster.start.indices[i]:(cluster.start.indices[i+1]-1) ) {
    # Only write it to file if it is a non-duplicate
    if ( !( petsys.time[j] == old.petsys.time & petsys.qdc[j] == old.petsys.qdc & petsys.ch[j] == old.petsys.ch ) ) {
      
      this.layer.number <- 0
      for ( layer.number in layer.numbers ) {
        if ( petsys.ch[j] > layer.min.channel.number[layer.number] & petsys.ch[j] < layer.max.channel.number[layer.number] ) {
          this.layer.number <- layer.number
        }
      }
      
      layer <- this.layer.number # first layer is 1 check for more than 2 layers!
      reduced.channel <- petsys.ch[j] - layer.subtract.this[layer]
      # reduced.channel <- ( petsys.ch[j] / 64 - floor( petsys.ch[j] / 64 ) ) * 64
      # bts.label <- mapping$bts.label[reduced.channel == mapping$ChNumber] # If you ever need it
      
      MPPC.x <- mapping$SiPMx[reduced.channel == mapping$ChNumber] # 1 to 8 # Location on MPPC
      MPPC.y <- mapping$SiPMy[reduced.channel == mapping$ChNumber] # 1 to 8 # Location on MPPC
      
      write( c( event.number, as.character( petsys.time[j] ), petsys.qdc[j], petsys.ch[j], layer, reduced.channel, MPPC.x, MPPC.y, petsys.number.of.hits[j] ), this.file, ncolumns = 9, append = TRUE, sep = ", " )
    }
    old.petsys.time <- petsys.time[j]
    old.petsys.qdc <- petsys.qdc[j]
    old.petsys.ch <- petsys.ch[j]
  }
  event.number <- event.number + 1
}
