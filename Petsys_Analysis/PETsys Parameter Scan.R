# Parameter scan
# Run this after PETsys Analysis loads all of the data

delta.x.add <- seq( from = -0.03, to = 0.03, by = 0.005 )
length( delta.x.add )
delta.y.add <- seq( from = 0.32, to = 0.41, by = 0.005 )
length( delta.y.add )

delta.x.stretch <- seq( from = 1.006, to = 1.012, by = 0.0005 )
length( delta.x.stretch )
delta.y.stretch <- seq( from = 1.000, to = 1.007, by = 0.0005 )
length( delta.y.stretch )

number.trials <- length( delta.x.add ) * length( delta.y.add ) * length( delta.x.stretch ) * length( delta.y.stretch )
number.trials

results.mat <- matrix( 0, nrow = number.trials, ncol = 5 )

these.breaks <- seq( from = 0, to = 360, by = 20 )
zenith.cut <- 50 * pi / 180

trial.number <- 1
min.reduced.chi.squared <- 1e6
for ( this.delta.x.add in delta.x.add ) {
  for ( this.delta.y.add in delta.y.add ) {
    for ( this.delta.x.stretch in delta.x.stretch ) {
      for ( this.delta.y.stretch in delta.y.stretch ) {
        
        azimuth <- atan2( delta.y * this.delta.y.stretch + this.delta.y.add, delta.x * this.delta.x.stretch + this.delta.x.add )
        azimuth[azimuth<0] <- 2 * pi + azimuth[azimuth<0] # For Bill
        
        #main.label <- paste0( "Histogram of Azimuth\n", this.delta.y.add )
        #hist( azimuth[muon.stdcut & zenith < zenith.cut] * 180 / pi, breaks = these.breaks, main = main.label, xlab = "Azimuth (degrees) East = 0 degrees" )
        h <- hist( azimuth[muon.stdcut & zenith < zenith.cut] * 180 / pi, breaks = these.breaks, plot = FALSE )
        #abline( h = mean( h$counts ), col = "blue" )
        
        this.reduced.chi.squared <- sum( ( h$counts - mean( h$counts ) )^2 / h$counts ) / ( length( h$counts ) - 1 )
        results.mat[trial.number,] <- c( this.delta.x.add, this.delta.y.add, this.delta.x.stretch, this.delta.y.stretch, this.reduced.chi.squared )
        if ( this.reduced.chi.squared < min.reduced.chi.squared ) {
          min.reduced.chi.squared <- this.reduced.chi.squared
          print( results.mat[trial.number,] )
        }
        trial.number <- trial.number + 1
      }
    }
  }
}
min.trial.number <- seq( results.mat[,5] )[results.mat[,5] == min( results.mat[,5] )]
results.mat[min.trial.number,]

print( paste0( "Min delta.x.add = ", results.mat[min.trial.number,1], " cm" ) )
print( paste0( "Min delta.y.add = ", results.mat[min.trial.number,2], " cm" ) )
print( paste0( "Min delta.x.stretch = ", results.mat[min.trial.number,3] ) )
print( paste0( "Min delta.y.stretch = ", results.mat[min.trial.number,4] ) )
print( paste0( "Min reduced chi-squared = ", results.mat[min.trial.number,5] ) )

stdcut <- results.mat[,2] == results.mat[min.trial.number,2]
stdcut <- stdcut & results.mat[,3] == results.mat[min.trial.number,3]
stdcut <- stdcut & results.mat[,4] == results.mat[min.trial.number,4]
plot( results.mat[stdcut,1], results.mat[stdcut,5], xlab = "delta.x.add (cm)", ylab = "reduced.chi.squared" )

stdcut <- results.mat[,1] == results.mat[min.trial.number,1]
stdcut <- stdcut & results.mat[,3] == results.mat[min.trial.number,3]
stdcut <- stdcut & results.mat[,4] == results.mat[min.trial.number,4]
plot( results.mat[stdcut,2], results.mat[stdcut,5], xlab = "delta.y.add (cm)", ylab = "reduced.chi.squared" )

stdcut <- results.mat[,1] == results.mat[min.trial.number,1]
stdcut <- stdcut & results.mat[,2] == results.mat[min.trial.number,2]
stdcut <- stdcut & results.mat[,4] == results.mat[min.trial.number,4]
plot( results.mat[stdcut,3], results.mat[stdcut,5], xlab = "delta.x.stretch (cm)", ylab = "reduced.chi.squared" )

stdcut <- results.mat[,1] == results.mat[min.trial.number,1]
stdcut <- stdcut & results.mat[,2] == results.mat[min.trial.number,2]
stdcut <- stdcut & results.mat[,3] == results.mat[min.trial.number,3]
plot( results.mat[stdcut,4], results.mat[stdcut,5], xlab = "delta.y.stretch (cm)", ylab = "reduced.chi.squared" )


this.delta.x.add <- results.mat[min.trial.number,1]
this.delta.y.add <- results.mat[min.trial.number,2]
this.delta.x.stretch <- results.mat[min.trial.number,3]
this.delta.y.stretch <- results.mat[min.trial.number,4]
this.reduced.chi.squared <- results.mat[min.trial.number,5]

azimuth <- atan2( delta.y * this.delta.y.stretch + this.delta.y.add, delta.x * this.delta.x.stretch + this.delta.x.add )
azimuth[azimuth<0] <- 2 * pi + azimuth[azimuth<0] # For Bill

main.label <- paste0( "Histogram of Azimuth\ndx add = ", this.delta.x.add * 10, " mm, dy add = ",this.delta.y.add * 10, " mm\ndx stretch = ", this.delta.x.stretch, ", dy stretch = ",this.delta.y.stretch, "\nReduced chi-squared = ", signif( this.reduced.chi.squared, 3 ) )
hist( azimuth[muon.stdcut & zenith < zenith.cut] * 180 / pi, breaks = these.breaks, main = main.label, xlab = "Azimuth (degrees) East = 0 degrees" )
h <- hist( azimuth[muon.stdcut & zenith < zenith.cut] * 180 / pi, breaks = these.breaks, plot = FALSE )
abline( h = mean( h$counts ), col = "blue" )



this.delta.x.add <- 0.00
this.delta.y.add <- 0.00
this.delta.x.stretch <- 1.000
this.delta.y.stretch <- 1.000

azimuth <- atan2( delta.y * this.delta.y.stretch + this.delta.y.add, delta.x * this.delta.x.stretch + this.delta.x.add )
azimuth[azimuth<0] <- 2 * pi + azimuth[azimuth<0] # For Bill

main.label <- paste0( "Histogram of Azimuth\ndx add = ", this.delta.x.add * 10, " mm, dy add = ",this.delta.y.add * 10, " mm\ndx stretch = ", this.delta.x.stretch, ", dy stretch = ",this.delta.y.stretch )
hist( azimuth[muon.stdcut & zenith < zenith.cut] * 180 / pi, breaks = these.breaks, main = main.label, xlab = "Azimuth (degrees) East = 0 degrees" )
h <- hist( azimuth[muon.stdcut & zenith < zenith.cut] * 180 / pi, breaks = these.breaks, plot = FALSE )
abline( h = mean( h$counts ), col = "blue" )
