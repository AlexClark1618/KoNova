# Load saved coincidences

print( "Loading Eventstats File" )

mat <- matrix( scan( file = paste0( program.root, stats.root, "/", filename, "-eventstats.dat" ), skip = 1, sep = "," ), ncol = 18, byrow = TRUE )
event.coinc.number <- mat[,1]
event.hits <- mat[,2]
event.delta.t <- mat[,3]
event.qdc <- mat[,4]
layer.hits <- mat[,c(5,7,9,11)]
layer.qdc <- mat[,c(6,8,10,12)]
xydetector.x <- mat[,c(13,15)]
xydetector.y <- mat[,c(14,16)]
adjacent.cut <- mat[,17] > 0.5
edge.cut <- mat[,18] > 0.5

