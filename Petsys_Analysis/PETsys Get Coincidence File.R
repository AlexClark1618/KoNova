# Get and save coincidences.  4 or more hardwired

print( "Getting Coincidence File" )

event.number.hits <- diff( cluster.start.indices )
length( event.number.hits )
length( cluster.start.times ) # lost 2
length( event.number.hits[event.number.hits > 3] )
hist( event.number.hits[event.number.hits > 3] )

stdcut <- event.number.hits > 3 # forget anything with less than 4 hits.  15k speed improvement!
good.indices <- seq( event.number.hits )[stdcut]
print( paste0( length( good.indices ), " coincidences" ) )

# Write out a coincidence file

this.file <- paste0( program.root, stats.root, "/", filename, "-coinc.dat" )
write( c( "Event Number", "Time (ps)", "QDC", "Channel", "Layer", "Reduced Channel", "MPPC X", "MPPC Y" ), this.file, ncolumns = 8, append = FALSE, sep = ", " )

event.number <- 1
for ( i in good.indices ) {
  for ( j in cluster.start.indices[i]:(cluster.start.indices[i+1]-1) ) {
    layer <- ( floor( single.ch[j] / 64 ) - 1 ) / 2 / 2 + 1 # first layer is 1 check for more than 2 layers!
    reduced.channel <- ( single.ch[j] / 64 - floor( single.ch[j] / 64 ) ) * 64
    # bts.label <- mapping$bts.label[reduced.channel == mapping$ChNumber] # If you ever need it
    MPPC.x <- mapping$SiPMx[reduced.channel == mapping$ChNumber] # 1 to 8 # Location on MPPC
    MPPC.y <- mapping$SiPMy[reduced.channel == mapping$ChNumber] # 1 to 8 # Location on MPPC
    
    write( c( event.number, as.character( single.time[j] ), single.qdc[j], single.ch[j], layer, reduced.channel, MPPC.x, MPPC.y ), this.file, ncolumns = 8, append = TRUE, sep = ", " )
  }
  event.number <- event.number + 1
}
