# Load saved coincidences

print( "Loading Hitstats File" )

mat <- matrix( scan( file = paste0( program.root, stats.root, "/", filename, "-hitstats.dat" ), skip = 1, sep = "," ), ncol = 9, byrow = TRUE )
coinc.number <- mat[,1]
coinc.time <- mat[,2]
coinc.qdc <- mat[,3]
coinc.channel <- mat[,4]
coinc.layer <- mat[,5]
coinc.reduced.channel <- mat[,6]
coinc.MPPC.x <- mat[,7]
coinc.MPPC.y <- mat[,8]
coinc.number.of.hits <- mat[,9]

#print( paste0( max( coinc.number ), " coincidences" ) )

run.time <- ( coinc.time[length(coinc.time)] - coinc.time[1] ) / 1e12 # in s
