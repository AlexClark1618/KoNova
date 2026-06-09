# Just in case

this.file <- paste0( program.root, data.root, "/", filename, "_single.dat" )

# Req Code; ID; RF; Cal; Ch, Locked, W#; t_ow mil; t_ow submil; Event #
mat <- matrix( scan( this.file, skip = 1, sep = "\t" ), ncol = 3, byrow = TRUE )
ps <- mat[,1]
QDC <- mat[,2]
ch <- mat[,3]
