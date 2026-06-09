# Use for PETsys data
# Initialize R

setwd( "/Users/ifft/Desktop/Muography/PETsys Analysis" )

rm( list = ls() )

program.root <- "/Users/ifft/Desktop/Muography/PETsys Analysis"
data.root <- "PETsys Data"
plot.root <- "PETsys Plots"
stats.root <- "PETsys Stats"

source( paste0( program.root, "/PETsys Initialize.R" ) )

# Load one single file

# filename.bkg <- "KNVA-20251203-22-00001" # background with thres = 14
filename.src <- "KNVA-20251203-23-00001" # source with thres = 14

# filename <- filename.bkg
# source( paste0( program.root, "/PETsys Load Single Data.R" ) )
# # hist( single.ch, breaks = seq( from = -0.5, to = 1024.5, by = 1 ), main = filename )
# # h.bkg <- hist( single.ch, breaks = seq( from = -0.5, to = 1024.5, by = 1 ), plot = FALSE )
# # hist( single.qdc[abs(single.qdc)<5], nclass = 1000 )
# single.ch.bkg <- single.ch
# single.qdc.bkg <- single.qdc
# single.time.bkg <- single.time

filename <- filename.src
source( paste0( program.root, "/PETsys Load Single Data.R" ) )
# hist( single.ch, breaks = seq( from = -0.5, to = 1024.5, by = 1 ), main = filename )
# h.src <- hist( single.ch, breaks = seq( from = -0.5, to = 1024.5, by = 1 ), plot = FALSE )
# hist( single.qdc[abs(single.qdc)<5], nclass = 1000 )
single.ch.src <- single.ch
single.qdc.src <- single.qdc
single.time.src <- single.time

# cbind( h.src$counts - h.bkg$counts )
# plot( h.src$counts - h.bkg$counts )
# 
# this.channel <- 119
# delta.bin <- 0.1
# filename <- filename.src
# source( paste0( program.root, "/PETsys Load Single Data.R" ) )
# hist( single.qdc[single.ch == this.channel], breaks = seq( from = -1, to = 20, by = delta.bin ) )
# h.src <- hist( single.qdc[single.ch == this.channel], breaks = seq( from = -1, to = 20, by = delta.bin ), plot = FALSE )
# filename <- filename.bkg
# source( paste0( program.root, "/PETsys Load Single Data.R" ) )
# hist( single.qdc[single.ch == this.channel], breaks = seq( from = -1, to = 20, by = delta.bin ) )
# h.bkg <- hist( single.qdc[single.ch == this.channel], breaks = seq( from = -1, to = 20, by = delta.bin ), plot = FALSE )
# cbind( h.src$counts - h.bkg$counts )
# plot( h.src$mids, h.src$counts - h.bkg$counts, type = "n", xlim = c( 0, 5 ) )
# lines( h.src$mids, h.src$counts - h.bkg$counts )
# 


# Motion analysis

delta.t <- 2 * 1e12 # in ps
start.ts <- seq( from = min( single.time.src ), to = max( single.time.src ), by = delta.t )

this.start.t <- start.ts[1]
stdcut <- single.time.src >= this.start.t & single.time.src < this.start.t + delta.t
stdcut <- stdcut & single.ch.src >= 64 & single.ch.src <= 127
hist( single.ch.src[stdcut] - 64, breaks = seq( from = -0.5, to = 63.5, by = 1 ), main = filename )
h.bkg <- hist( single.ch.src[stdcut] - 64, seq( from = -0.5, to = 63.5, by = 1 ), plot = FALSE ) # assuming the source starts away from the detector

channel.numbers <- c( 48, 47, 44, 42, 51, 54, 55, 60 )
y.offset <- -200
delta.y <- -60
for ( this.start.t in start.ts[2:(length(start.ts)-1)] ) {
    stdcut <- single.time.src >= this.start.t & single.time.src < this.start.t + delta.t
    stdcut <- stdcut & single.ch.src >= 64 & single.ch.src <= 127
    h.src <- hist( single.ch.src[stdcut] - 64, breaks = seq( from = -0.5, to = 63.5, by = 1 ), plot = FALSE )
    plot( h.src$mids, h.src$counts - h.bkg$counts, ylim = c( -700, 2000 ), main = filename )
    this.y <- y.offset
    for ( this.bar in 1:8 ) {
      abline( v = channel.numbers[this.bar], col = "gray" )
      text( channel.numbers[this.bar], this.y, paste0( "Bar ", this.bar ), cex = 0.5 )
      this.y <- this.y + delta.y
    }
}

