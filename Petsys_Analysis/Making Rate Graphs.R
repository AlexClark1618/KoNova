# Use for PETsys data
# This file uses PETsys _coinc.dat compact text files
# Initialize R

setwd( "/Users/ifft/Desktop/Muography/PETsys Analysis" )

rm( list = ls() )

program.root <- getwd()

data.root <- "/PETsys Data"
stats.root <- "/PETsys Stats"
plot.root <- "/PETsys Plots"

source( "/users/ifft/documents/Error Analysis/measurement 2.0.r" )

mat <- matrix( scan( file = "PETsys Coincidence Data.csv", sep = ",", skip = 2 ), byrow = TRUE, ncol = 11 )
exp.time <- mat[,1]
vth.e <- mat[,2]
file.size.kb <- mat[,5]
number.accepted.2222 <- mat[,6]
number.accepted.1111 <- mat[,9]

plot( vth.e, file.size.kb, xlab = "Vth_e", ylab = "File size (kB)", main = "File Size vs Vth_e", col = "red", log = "y" )

rate.2222 <- measurement( number.accepted.2222, sqrt( number.accepted.2222 ) ) / exp.time
plot( measurement( vth.e, 0 ), rate.2222, ylim = c( 0, 0.5 ), xlab = "Vth_e", ylab = "2-2-2-2 Rate (Hz)", main = "2-2-2-2 Rate vs Vth_e" )
abline( h = 0.44, col = "red" )

rate.1111 <- measurement( number.accepted.1111, sqrt( number.accepted.1111 ) ) / exp.time
plot( measurement( vth.e, 0 ), rate.1111, ylim = c( 0, 1.5 ), xlab = "Vth_e", ylab = "1-1-1-1 Rate (Hz)", main = "1-1-1-1 Rate vs Vth_e" )
abline( h = 1.02, col = "red" )

plot( measurement( vth.e, 0 ), rate.1111, type = "n", xlab = "Vth_e", ylab = "Rate (Hz)", main = "Rates vs Vth_e" )
points( measurement( vth.e, 0 ), rate.1111, col = "blue" )
points( measurement( vth.e, 0 ), rate.2222, col = "green" )
legend( "topright", c( "1111", "2222" ), text.col = c( "blue", "green" ) )


mat <- matrix( scan( file = "PETsys Singles Data.csv", sep = ",", skip = 2 ), byrow = TRUE, ncol = 6 )
exp.time <- mat[,1]
vth.e <- mat[,2]
number.accepted <- mat[,5]

rate.singles <- number.accepted / exp.time
plot( vth.e, rate.singles / 1e6, xlab = "Vth_e", ylab = "Rate (MHz)", main = "Singles Rate vs Vth_e", col = "brown" )

