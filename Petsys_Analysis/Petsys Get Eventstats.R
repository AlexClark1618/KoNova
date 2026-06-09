# Get event statistics

print( "Getting Event Stats" )

event.coinc.number <- unique( coinc.number ) * 0 - 1e6 # Link back to coinc data
event.hits <- unique( coinc.number ) * 0 + 1e6 # overall number of hits, 4 or greater because of above
event.delta.t <- unique( coinc.number ) * 0 + 1e6 # max delta t for all hits
event.qdc <- unique( coinc.number ) * 0 + 1e6 # total qdc
layer.hits <- matrix( -1e6, nrow = length( unique( coinc.number ) ), ncol = 4 )
layer.qdc <- matrix( -1e6, nrow = length( unique( coinc.number ) ), ncol = 4 )
# layer.min.qdc <- matrix( -1e6, nrow = length( unique( coinc.number ) ), ncol = 4 )
# layer.max.qdc <- matrix( -1e6, nrow = length( unique( coinc.number ) ), ncol = 4 )
xydetector.x <- matrix( -1e6, nrow = length( unique( coinc.number ) ), ncol = 2 )
xydetector.y <- matrix( -1e6, nrow = length( unique( coinc.number ) ), ncol = 2 )
adjacent.cut <- seq( event.hits ) > -1 # All TRUE to start.  False if on any layer there are non-adjacent hits
edge.cut <- seq( event.hits ) > -1 # All TRUE to start.  False if on any layer there are non-adjacent hits

this.file <- paste0( program.root, stats.root, "/", filename, "-eventstats.dat" )
write( c( "Event Number", "Hits", "Delta t", "QDC", "Layer 1 Hits", "Layer 1 QDC", "Layer 2 Hits", "Layer 2 QDC", "Layer 3 Hits", "Layer 3 QDC", "Layer 4 Hits", "Layer 4 QDC", "x1 (cm)", "y1 (cm)", "x2 (cm)", "y2 (cm)", "Adjacent Cut", "Edge Cut" ), this.file, ncolumns = 18, append = FALSE, sep = ", " )

percent.print <- 0.01
# count.12 <- 0
# count.23 <- 0
for ( this.coinc.number in unique( coinc.number ) ) {
  if ( this.coinc.number / length( unique( coinc.number ) ) > percent.print ) {
    print( paste0( percent.print * 100, "% done" ) )
    percent.print <- 2 * percent.print
  }

  #coinc.stdcut <- coinc.number == this.coinc.number
  coinc.seq <- seq( coinc.number )[coinc.number == this.coinc.number]
  
  # event statistics
  event.coinc.number <- this.coinc.number
  event.hits[this.coinc.number] <- length( coinc.time[coinc.seq] )
  event.delta.t[this.coinc.number] <- diff( range( coinc.time[coinc.seq] ) )
  event.qdc[this.coinc.number] <- sum( coinc.qdc[coinc.seq] )
  
  # layer statistics, bottom to top
  for ( layer in layer.numbers ) {
    # this.stdcut <- coinc.layer == layer & coinc.stdcut
    this.seq <- coinc.seq[coinc.layer[coinc.seq]==layer]
    if ( length( this.seq ) > 0 ) {
      layer.hits[this.coinc.number, layer] <- length( this.seq )
      layer.qdc[this.coinc.number, layer] <- sum( coinc.qdc[this.seq] )
      # layer.min.qdc[this.coinc.number, layer] <- min( coinc.qdc[coinc.stdcut & layer.stdcut] )
      # layer.max.qdc[this.coinc.number, layer] <- max( coinc.qdc[coinc.stdcut & layer.stdcut] )
    }
  }
  
  # xydetector statistics
  # in cm
  # this analysis takes 1 bar hits
  # reduced channels goes from 0 to 63 but bar.map goes from 1 to 64!
  
  x1 <- -1e6
  y1 <- -1e6
  x2 <- -1e6
  y2 <- -1e6
  
  this.layer <- 1
  this.seq <- coinc.seq[coinc.layer[coinc.seq]==this.layer]
  these.bars <- ch.2.bar.mat[coinc.reduced.channel[this.seq] + 1,this.layer]
  if ( length( these.bars ) > 0 ) {
    x1 <- weighted.mean( these.bars - 1, w = coinc.qdc[this.seq] ) * channel.spacing # this is crude and the default
    xydetector.x[this.coinc.number, 1] <- x1
  }
  
  this.layer <- 2
  this.seq <- coinc.seq[coinc.layer[coinc.seq]==this.layer]
  these.bars <- ch.2.bar.mat[coinc.reduced.channel[this.seq] + 1,this.layer]
  if ( length( these.bars ) > 0 ) {
    y1 <- weighted.mean( these.bars - 1, w = coinc.qdc[this.seq] ) * channel.spacing
    xydetector.y[this.coinc.number, 1] <- y1
  }
  
  this.layer <- 3
  this.seq <- coinc.seq[coinc.layer[coinc.seq]==this.layer]
  these.bars <- ch.2.bar.mat[coinc.reduced.channel[this.seq] + 1,this.layer]
  if ( length( these.bars ) > 0 ) {
    x2 <- weighted.mean( these.bars - 1, w = coinc.qdc[this.seq] ) * channel.spacing
    xydetector.x[this.coinc.number, 2] <- x2
  }
  
  this.layer <- 4
  this.seq <- coinc.seq[coinc.layer[coinc.seq]==this.layer]
  these.bars <- ch.2.bar.mat[coinc.reduced.channel[this.seq] + 1,this.layer]
  if ( length( these.bars ) > 0 ) {
    y2 <- weighted.mean( these.bars - 1, w = coinc.qdc[this.seq] ) * channel.spacing
    xydetector.y[this.coinc.number, 2] <- y2
  }
  
  # Do one more iteration for 2 bars for now
  
  if ( x1 > -1e5 & y1 > -1e5 & x2 > -1e5 & y2 > -1e5 & !is.na( x1 ) & !is.na( y1 ) & !is.na( x2 ) & !is.na( y2 ) ) {

    delta.x <- x2 - x1
    delta.y <- y2 - y1
    zenith <- atan( sqrt( delta.x^2 + delta.y^2 ) / delta.z )
    azimuth <- atan2( delta.y, delta.x )

    this.layer <- 1
    this.seq <- coinc.seq[coinc.layer[coinc.seq]==this.layer]
    these.bars <- ch.2.bar.mat[coinc.reduced.channel[this.seq] + 1,this.layer]
    these.qdcs <- coinc.qdc[this.seq]
    if ( length( these.bars ) == 1 ) {
      x1 <- ( these.bars[1] - 1 + (  runif( 1 ) - 0.5 ) ) * channel.spacing
      xydetector.x[this.coinc.number, 1] <- x1
    } else if ( length( these.bars ) == 2 ) {
      tan.eta <- tan( zenith ) * cos( azimuth )
      increasing.position.order <- order( these.bars )
      these.bars <- these.bars[increasing.position.order]
      these.qdcs <- these.qdcs[increasing.position.order]
      if ( these.bars[1] %% 2 == 0 ) {  # First bar is even
        UD <- -1
      } else {
        UD <- 1
      }
      # this.offset <- layer.qdc.offsets[layer]
      # QR <- ( these.qdcs[2] - this.offset ) / ( sum( these.qdcs ) - 2 * this.offset )
      QR <- runif( 1 )
      
      x1 <- ( these.bars[1] - 1 ) * channel.spacing + UD * ( layer.height / 2.0 ) * tan.eta + QR * ( channel.spacing - UD * layer.height * tan.eta )
      # print( paste0( "Old x1 = ", round( xydetector.x[this.coinc.number, 1], 1 ), ", New x1 = ", round( x1, 1 ) ) )
      xydetector.x[this.coinc.number, 1] <- x1
    }

    this.layer <- 2
    this.seq <- coinc.seq[coinc.layer[coinc.seq]==this.layer]
    these.bars <- ch.2.bar.mat[coinc.reduced.channel[this.seq] + 1,this.layer]
    these.qdcs <- coinc.qdc[this.seq]
    if ( length( these.bars ) == 1 ) {
      y1 <- ( these.bars[1] - 1 + ( runif( 1 ) - 0.5 ) ) * channel.spacing
      xydetector.y[this.coinc.number, 1] <- y1
    } else if ( length( these.bars ) == 2 ) {
      tan.eta <- tan( zenith ) * sin( azimuth )
      increasing.position.order <- order( these.bars )
      these.bars <- these.bars[increasing.position.order]
      these.qdcs <- these.qdcs[increasing.position.order]
      if ( these.bars[1] %% 2 == 0 ) {  # First bar is even
        UD <- 1
      } else {
        UD <- -1
      }
      # this.offset <- layer.qdc.offsets[layer]
      # QR <- ( these.qdcs[2] - this.offset ) / ( sum( these.qdcs ) - 2 * this.offset )
      QR <- runif( 1 )
      
      y1 <- ( these.bars[1] - 1 ) * channel.spacing + UD * ( layer.height / 2.0 ) * tan.eta + QR * ( channel.spacing - UD * layer.height * tan.eta )
      # print( paste0( "Old y1 = ", round( xydetector.y[this.coinc.number, 1], 1 ), ", New y1 = ", round( y1, 1 ) ) )
      xydetector.y[this.coinc.number, 1] <- y1
    }

    this.layer <- 3
    this.seq <- coinc.seq[coinc.layer[coinc.seq]==this.layer]
    these.bars <- ch.2.bar.mat[coinc.reduced.channel[this.seq] + 1,this.layer]
    these.qdcs <- coinc.qdc[this.seq]
    if ( length( these.bars ) == 1 ) {
      x2 <- ( these.bars[1] - 1 + (  runif( 1 ) - 0.5 ) ) * channel.spacing
      xydetector.x[this.coinc.number, 2] <- x2
    } else if ( length( these.bars ) == 2 ) {
      tan.eta <- tan( zenith ) * cos( azimuth )
      increasing.position.order <- order( these.bars )
      these.bars <- these.bars[increasing.position.order]
      these.qdcs <- these.qdcs[increasing.position.order]
      if ( these.bars[1] %% 2 == 0 ) {  # First bar is even
        UD <- -1
      } else {
        UD <- 1
      }
      # this.offset <- layer.qdc.offsets[layer]
      # QR <- ( these.qdcs[2] - this.offset ) / ( sum( these.qdcs ) - 2 * this.offset )
      QR <- runif( 1 )
      
      x2 <- ( these.bars[1] - 1 ) * channel.spacing + UD * ( layer.height / 2.0 ) * tan.eta + QR * ( channel.spacing - UD * layer.height * tan.eta )
      # print( paste0( "Old x2 = ", round( xydetector.x[this.coinc.number, 2], 1 ), ", New x2 = ", round( x2, 1 ) ) )
      xydetector.x[this.coinc.number, 2] <- x2
    }

    this.layer <- 4
    this.seq <- coinc.seq[coinc.layer[coinc.seq]==this.layer]
    these.bars <- ch.2.bar.mat[coinc.reduced.channel[this.seq] + 1,this.layer]
    these.qdcs <- coinc.qdc[this.seq]
    if ( length( these.bars ) == 1 ) {
      y2 <- ( these.bars[1] - 1 + (  runif( 1 ) - 0.5 ) ) * channel.spacing
      xydetector.y[this.coinc.number, 2] <- y2
    } else if ( length( these.bars ) == 2 ) {
      tan.eta <- tan( zenith ) * sin( azimuth )
      increasing.position.order <- order( these.bars )
      these.bars <- these.bars[increasing.position.order]
      these.qdcs <- these.qdcs[increasing.position.order]
      if ( these.bars[1] %% 2 == 0 ) {  # First bar is even
        UD <- 1
      } else {
        UD <- -1
      }
      # this.offset <- layer.qdc.offsets[layer]
      # QR <- ( these.qdcs[2] - this.offset ) / ( sum( these.qdcs ) - 2 * this.offset )
      QR <- runif( 1 )
      
      y2 <- ( these.bars[1] - 1 ) * channel.spacing + UD * ( layer.height / 2.0 ) * tan.eta + QR * ( channel.spacing - UD * layer.height * tan.eta )
      # print( paste0( "Old y2 = ", round( xydetector.y[this.coinc.number, 2], 1 ), ", New y2 = ", round( y2, 1 ) ) )
      xydetector.y[this.coinc.number, 2] <- y2
    }
  }

  this.min.time <- min( coinc.time[coinc.seq] )
  
  for ( this.layer.number in layer.numbers ) {
    this.seq <- coinc.seq[coinc.layer[coinc.seq]==this.layer.number]
    
    these.bars <- ch.2.bar.mat[coinc.reduced.channel[this.seq] + 1, this.layer.number]
    if( length( these.bars ) >= 2 ) {
      if( max( diff( sort( these.bars ) ) ) >= 2 ) {
        adjacent.cut[this.coinc.number] <- FALSE
      }
    }
    
    if( length( these.bars ) > 0 ) {
      if ( min( these.bars ) == 1 ) {
        edge.cut[this.coinc.number] <- FALSE
      }
      
      if ( max( these.bars ) == 64 ) {
        edge.cut[this.coinc.number] <- FALSE
      }
    }
  }
  
  write( c( event.coinc.number, event.hits[this.coinc.number], event.delta.t[this.coinc.number], event.qdc[this.coinc.number], layer.hits[this.coinc.number, 1], layer.qdc[this.coinc.number, 1], layer.hits[this.coinc.number, 2], layer.qdc[this.coinc.number, 2], layer.hits[this.coinc.number, 3], layer.qdc[this.coinc.number, 3], layer.hits[this.coinc.number, 4], layer.qdc[this.coinc.number, 4], xydetector.x[this.coinc.number, 1], xydetector.y[this.coinc.number, 1], xydetector.x[this.coinc.number, 2], xydetector.y[this.coinc.number, 2], adjacent.cut[this.coinc.number], edge.cut[this.coinc.number] ), this.file, ncolumns = 18, append = TRUE, sep = ", " )
  
}


