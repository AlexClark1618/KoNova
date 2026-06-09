# Use for PETsys data
# Initialize R

setwd( "/Users/ifft/Desktop/Muography/PETsys Analysis" )

rm( list = ls() )

vertical.flux <- 70 / 100^2 # in number per cm^2 per s per sr
number.of.boxes <- 50
# L.x <- 14.85 # in cm for mini-KoNova
# L.y <- 14.85
L.x <- ( 64 / 2 ) * ( 33 / 10 ) - 1 * ( 33 / 10 ) # in cm for KoNova
L.y <- ( 64 / 2 ) * ( 33 / 10 ) - 1 * ( 33 / 10 )
box.size.x <- L.x / number.of.boxes
box.size.y <- L.y / number.of.boxes
box.area <- box.size.x * box.size.y
x.bb <- seq( from = box.size.x / 2, to = L.x - box.size.x / 2, by = box.size.x )
y.bb <- seq( from = box.size.y / 2, to = L.y - box.size.y / 2, by = box.size.y )
#delta.z <- 1.7 * 2 + 1000 + 1.7 * 2 # in cm top of top layer to bottom of bottom layer for mini-KoNova.  Enter gap
delta.z <- 30.0 # in cm top of top layer to bottom of bottom layer

total.rate <- 0.0 # in number per s
for ( this.x.bb in x.bb ) {
  for ( this.y.bb in y.bb ) {
    for ( this.x.tt in x.bb ) {
      for ( this.y.tt in y.bb ) {
        this.zenith <- atan( sqrt( (this.x.tt - this.x.bb)^2 + (this.y.tt - this.y.bb)^2 ) / delta.z )
        this.sr <- box.area * cos( this.zenith ) / ( (this.x.tt - this.x.bb)^2 + (this.y.tt - this.y.bb)^2 + delta.z^2 )
        total.rate <- total.rate + vertical.flux * ( cos( this.zenith ) )^2 * this.sr * box.area * cos( this.zenith )
      }
    }
  }
}
total.rate


vertical.flux <- 70 / 100^2 # in number per cm^2 per s per sr
delta.zenith <- 0.1 * pi / 180
zeniths <- seq( from = 0, to = 90 * pi / 180, by = delta.zenith )
dN.dzenith <- 2 * pi * vertical.flux * cos( zeniths )^2 * cos( zeniths ) * sin( zeniths )
plot( zeniths * 180 / pi, dN.dzenith )
sum( dN.dzenith ) * delta.zenith

