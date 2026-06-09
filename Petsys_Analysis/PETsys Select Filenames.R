# Get filenames and orders them to analyze

filenames <- NULL

# search.string <- "KNVA-20260425-02" # Lab ~1 hour run
# search.string <- "KNVA-20260425-03" # Lab ~1 day run redone
# search.string <- "KNVA-20260501-01" # Cliff run ~1.5 of data (subset)
# search.string <- "KNVA-20260502-03" # Cliff run ~15 minutes
#search.string <- "KNVA-20260502-04" # Cliff run ~1 hour
# search.string <- "KNVA-20260503-06" # Testing with old, default threshold setting 20, 20, 15
# search.string <- "KNVA-20260503-03" # old threshold at 40
# search.string <- "KNVA-20260503-07" # Testing with old, default threshold setting 20, 20, 10
# search.string <- "KNVA-20260503-08" # Testing with old, default threshold setting 20, 15, 10
# search.string <- "KNVA-20260503-09" # Back to 50
# search.string <- "KNVA-20260504-02" # Using normal thresholding 20, 20, 15
# search.string <- "KNVA-20260505-01" # Cliff run ~15 min (subset)
# search.string <- "KNVA-20260505-03" # Parking lot run ~5 min (subset)
# search.string <- "KNVA-20260506-05" # Parking lot Rotated 90 100 s cycles
# search.string <- "KNVA-20260506-06" # Parking lot Rotated 90 600 s cycles
# search.string <- "KNVA-20260509-01" # Parking lot Rotated 90 100 s cycles with 50 0 0 thresholds
# search.string <- "KNVA-20260509-02" # Parking lot Rotated 90 600 s cycles with 50 0 0 thresholds
# search.string <- "KNVA-20260511-07" # Blue Sky Run ~2 hours
# search.string <- "KNVA-20260511-09" # Hammer throw field, 600 s cycles
 search.string <- "KNVA-20260514-01" # Loading dock
analysis.indicator <- "Practice Run"

temp <- list.files( path = paste0( program.root, data.root ), pattern = "KNVA" )

filenames <- c( filenames, temp[grep( search.string, temp )] )
filenames

filenames <- substr( filenames, start = 1, stop = 22 )

run.names <- substr( filenames, start = 1, stop = 16 )

run.roots <- substr( filenames, start = 1, stop = 16 )
run.roots <- unique( run.roots )
run.roots
length( run.roots )

print( filenames )

rm( temp )
