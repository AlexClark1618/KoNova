# Close old graphics devices
graphics.off()

# Make plots directory

this.plot.directory <- paste0( program.root, plot.root, "/", analysis.indicator )

if ( dir.exists( this.plot.directory ) ) {
  #system( paste0( "rm -r ", shQuote( this.plot.directory ) ) )
  dir_delete( this.plot.directory )
}

dir.create( this.plot.directory )

# =========================================================
# Event Delta T
# =========================================================

main.label <- paste0( "Histogram of Event Delta Ts\n", analysis.indicator )

pdf(
  paste0( this.plot.directory, "/", analysis.indicator, "-DeltaT.pdf" ),
  width = 10,
  height = 8
)

hist(
  event.delta.t[muon.stdcut] / 1000,
  xlab = "Event delta.t (ns)",
  main = main.label,
  nclass = 100
)

dev.off()

# =========================================================
# Event QDC
# =========================================================

main.label <- paste0( "Histogram of Event QDCs\n", analysis.indicator )

pdf(
  paste0( this.plot.directory, "/", analysis.indicator, "-EventQDC.pdf" ),
  width = 10,
  height = 8
)

hist(
  event.qdc[muon.stdcut],
  xlab = "Event QDCs",
  main = main.label,
  nclass = 100
)

dev.off()

# =========================================================
# Layer statistics
# Hits
# =========================================================

pdf(
  paste0( this.plot.directory, "/", analysis.indicator, "-LayerHits.pdf" ),
  width = 12,
  height = 10
)

par(
  mfrow = c( 2, 2 ),
  mar = c( 3, 3, 2, 1 ),
  oma = c( 0, 0, 2, 0 )
)

for ( layer in layer.numbers ) {
  
  this.hits <- layer.hits[muon.stdcut,layer]
  
  hist(
    this.hits[this.hits < 6],
    main = paste0( "Layer ", layer )
  )
}

mtext(
  analysis.indicator,
  side = 3,
  line = 0.5,
  outer = TRUE,
  cex = 1.5,
  font = 2
)

par( mfrow = c( 1, 1 ), oma = c( 0, 0, 0, 0 ) )

dev.off()

# =========================================================
# QDC
# =========================================================

max.qdc <- 7.5

pdf(
  paste0( this.plot.directory, "/", analysis.indicator, "-LayerQDC.pdf" ),
  width = 12,
  height = 10
)

par(
  mfrow = c( 2, 2 ),
  mar = c( 3, 3, 2, 1 ),
  oma = c( 0, 0, 2, 0 )
)

for ( layer in layer.numbers ) {
  
  x <- layer.qdc[muon.stdcut,layer]
  x <- x[x > 0 & x < max.qdc]
  
  hist(
    x,
    xlab = "QDC after cuts",
    main = paste0(
      "Layer = ",
      layer,
      "\nmean = ",
      round( mean( x ), 2 )
    ),
    xlim = c( 0, max.qdc ),
    nclass = 50
  )
}

mtext(
  analysis.indicator,
  side = 3,
  line = 0.5,
  outer = TRUE,
  cex = 1.5,
  font = 2
)

par( mfrow = c( 1, 1 ), oma = c( 0, 0, 0, 0 ) )

dev.off()

# =========================================================
# XY
# =========================================================

pdf(
  paste0( this.plot.directory, "/", analysis.indicator, "-XY.pdf" ),
  width = 12,
  height = 10
)

main.label <- "Bottom x Layer"

hist(
  x1[muon.stdcut & x1 > -1 & x1 < 120],
  main = main.label,
  xlab = "x1 in cm",
  nclass = 400
)

main.label <- "Bottom x Layer (Zoom)"

hist(
  x1[muon.stdcut & x1 > -1 & x1 < 25],
  main = main.label,
  xlab = "x1 in cm",
  nclass = 100
)

abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

hist(
  x1[muon.stdcut & x1 > 80 & x1 < 107],
  main = main.label,
  xlab = "x1 in cm",
  nclass = 100
)

abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

main.label <- "Bottom y Layer"

hist(
  y1[muon.stdcut & y1 > -1 & y1 < 120],
  main = main.label,
  xlab = "y1 in cm",
  nclass = 400
)

main.label <- "Bottom y Layer (Zoom)"

hist(
  y1[muon.stdcut & y1 > -1 & y1 < 25],
  main = main.label,
  xlab = "y1 in cm",
  nclass = 100
)

abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

hist(
  y1[muon.stdcut & y1 > 80 & y1 < 107],
  main = main.label,
  xlab = "y1 in cm",
  nclass = 100
)

abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

main.label <- "Top x Layer"

hist(
  x2[muon.stdcut & x2 > -1 & x2 < 120],
  main = main.label,
  xlab = "x2 in cm",
  nclass = 400
)

main.label <- "Top x Layer (Zoom)"

hist(
  x2[muon.stdcut & x2 > -1 & x2 < 25],
  main = main.label,
  xlab = "x2 in cm",
  nclass = 100
)

abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

hist(
  x2[muon.stdcut & x2 > 80 & x2 < 107],
  main = main.label,
  xlab = "x2 in cm",
  nclass = 100
)

abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

main.label <- "Top y Layer"

hist(
  y2[muon.stdcut & y2 > -1 & y2 < 120],
  main = main.label,
  xlab = "y2 in cm",
  nclass = 400
)

main.label <- "Top y Layer (Zoom)"

hist(
  y2[muon.stdcut & y2 > -1 & y2 < 25],
  main = main.label,
  xlab = "y2 in cm",
  nclass = 100
)

abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

hist(
  y2[muon.stdcut & y2 > 80 & y2 < 107],
  main = main.label,
  xlab = "y2 in cm",
  nclass = 100
)

abline( v = ( 1:64 - 1 ) * channel.spacing, col = "blue" )

dev.off()

# =========================================================
# Get intercepts
# =========================================================

delta.break <- 0.1
intercepts <- NULL

pdf(
  paste0( this.plot.directory, "/", analysis.indicator, "-Intercepts.pdf" ),
  width = 12,
  height = 10
)

par(
  mfrow = c( 2, 2 ),
  mar = c( 3, 3, 2, 1 ),
  oma = c( 0, 0, 2, 0 )
)

for ( layer in layer.numbers ) {
  
  x <- 1 / cos( zenith[muon.stdcut] )
  y <- layer.qdc[muon.stdcut, layer]
  
  mean.y <- NULL
  
  breaks = seq( from = 1.0, to = 2.0, by = delta.break )
  
  for ( this.break in breaks ) {
    
    stdcut <- x > this.break &
      x < this.break + delta.break &
      y > 0 &
      y < 15
    
    mean.y <- c( mean.y, mean( y[stdcut] ) )
  }
  
  ls.fit <- lsfit( breaks + 0.05, mean.y )
  
  plot(
    breaks + 0.05,
    mean.y,
    xlab = "mean( 1 / cos( zenith ) )",
    ylab = "Layer QDC",
    main = paste0(
      "Layer ",
      layer,
      "\nIntercept = ",
      signif( ls.fit$coefficients[1], 4 )
    )
  )
  
  abline( ls.fit )
  
  intercepts <- c(
    intercepts,
    ls.fit$coefficients[1]
  )
}

mtext(
  analysis.indicator,
  side = 3,
  line = 0.5,
  outer = TRUE,
  cex = 1.5,
  font = 2
)

par( mfrow = c( 1, 1 ), oma = c(0, 0, 0, 0) )

dev.off()

# =========================================================
# Get inter-fiber frequencies
# =========================================================

pdf(
  paste0( this.plot.directory, "/", analysis.indicator, "-InterFiber.pdf" ),
  width = 12,
  height = 10
)

par(
  mfrow = c( 2, 2 ),
  mar = c( 3, 3, 2, 1 ),
  oma = c( 0, 0, 2, 0 )
)

hist(
  x1[muon.stdcut & x1 > -1 & x1 < 120] %% channel.spacing,
  main = main.label,
  xlab = "x1 in cm",
  nclass = 100
)

hist(
  y1[muon.stdcut & y1 > -1 & y1 < 120] %% channel.spacing,
  main = main.label,
  xlab = "y1 in cm",
  nclass = 100
)

hist(
  x2[muon.stdcut & x2 > -1 & x2 < 120] %% channel.spacing,
  main = main.label,
  xlab = "x2 in cm",
  nclass = 100
)

hist(
  y2[muon.stdcut & y2 > -1 & y2 < 120] %% channel.spacing,
  main = main.label,
  xlab = "y2 in cm",
  nclass = 100
)

mtext(
  analysis.indicator,
  side = 3,
  line = 0.5,
  outer = TRUE,
  cex = 1.5,
  font = 2
)

par( mfrow = c( 1, 1 ), oma = c( 0, 0, 0, 0 ) )

dev.off()