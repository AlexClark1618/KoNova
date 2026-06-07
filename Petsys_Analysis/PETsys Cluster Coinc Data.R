# Cluster data from PETsys coincidence file
# Sorts and clusters data based on cluster.gap

print( "Clustering Coincidence Data" )

stdcut <- c( TRUE, diff( petsys.time ) > cluster.gap )
cluster.start.indices <- seq( petsys.time )[stdcut] # start of a new cluster
cluster.start.times <- petsys.time[stdcut] # start of cluster in ps
cluster.start.times <- c( cluster.start.times, max( petsys.time ) + 1 )
number.of.clusters <- length( cluster.start.times ) # number of clusters
print( paste0( "Number of hits = ", length( petsys.time ) ) )
print( paste0( "Number of clusters = ", number.of.clusters ) )
length( petsys.time ) - number.of.clusters
