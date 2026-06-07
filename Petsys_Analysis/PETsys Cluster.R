# Cluster data
# Sorts and clusters data based on cluster.gap

print( "Clustering Data" )

stdcut <- c( TRUE, diff( single.time ) > cluster.gap )
cluster.start.indices <- seq( single.time )[stdcut] # start of a new cluster
cluster.start.times <- single.time[stdcut] # start of cluster in ps
cluster.start.times <- c( cluster.start.times, max( single.time ) + 1 )
number.of.clusters <- length( cluster.start.times ) # number of clusters
length( single.time )
number.of.clusters
length( single.time ) - number.of.clusters
