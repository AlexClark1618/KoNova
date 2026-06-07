# Load file

print( "Loading File" )

mat <- matrix( scan( file = paste0( program.root, data.root, "/", filename, "_single.dat" ) ), ncol = 3, byrow = TRUE )
single.time <- mat[,1]
single.qdc <- mat[,2]
single.ch <- mat[,3]
range( diff( single.time ) ) 
hist( diff( single.time )[diff( single.time ) < 10000] ) # a tiny fraction are out of order
rm( mat )

# Permanently changed from original file

increasing.time.order <- order( single.time )
single.time <- single.time[increasing.time.order]
single.qdc <- single.qdc[increasing.time.order]
single.ch <- single.ch[increasing.time.order]
range( diff( single.time ) ) 
hist( diff( single.time )[diff( single.time ) < 10000] )





