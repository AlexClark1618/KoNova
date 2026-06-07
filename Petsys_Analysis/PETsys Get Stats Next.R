# Initialize Needs updating

# cd "/Users/ifft/Desktop/Muography/PETsys Analysis"
# mkdir Run1
# cd Run1
# RScript '/Users/ifft/Desktop/Muography/PETsys Analysis/PETsys Get Stats Next.R'

rm( list = ls() )

program.root <- "C:/Users/alexc/Desktop/Petsys_Analysis"

data.root <- "/PETsys Data"
stats.root <- "/PETsys Stats"
plot.root <- "/PETsys Plots"

#source( "/users/ifft/documents/Error Analysis/measurement 2.0.r" )

# Run "GPS Select Filenames.R

source( paste0( program.root, '/PETsys Select Filenames.R' ) )

initialize.file <- paste0( program.root, "/", "PETsys Initialize.R" )
for ( file.i in seq( filenames ) ) {
  
  filename <- filenames[file.i]
  
  run.name <- substr( filename, 1, 16 )
  this.file <- paste0( program.root, stats.root, "/", run.name, "-ini.R" )
  if ( file.exists( this.file ) ) {
    print( "Using existing Initialize.R file")
    source( this.file )
  } else {
    print( "Making a new Initialize.R file")
    if ( substr( filename, 1, 4 ) == "KNMC" ) {
      source.file <- paste0( program.root, "/PETsys KNMC Initialize.R" )
    } else {
      source.file <- paste0( program.root, "/PETsys Initialize.R" )
    }
    source( source.file )
    file.copy( source.file, this.file )
  }
  
  # Set up LR directory if needed
  # LR.directory.name <- paste0( this.run.root, " LR directory" )
  # LR.directory.path <- paste0( program.root, plot.root, "/", LR.directory.name )
  # if ( !dir.exists( LR.directory.path ) ) dir.create( LR.directory.path )
  
  if ( file.exists( paste0( program.root, stats.root, "/", filenames[file.i], "-hitstats.dat" ) ) ) {
    print( paste0( "Skipping ", filenames[file.i], ", hitstats file already exists." ) )
  } else {
    
    # Establish this file for this Run Number.  This will be re-written later.  Just a placeholder
    this.file <- paste0( program.root, stats.root, "/", filename, "-hitstats.dat" )
    write( c( "Event Number", "Time (ps)", "QDC", "Channel", "Layer", "Reduced Channel", "MPPC X", "MPPC Y" ), this.file, ncolumns = 8, append = FALSE, sep = ", " )
    
    # Using the PETsys generated coinc file
    
    source( paste0( program.root, "/PETsys Load PETsys Coinc Data.R" ) )
    
    # coinc.81.87 <- 0
    # coinc.84.87 <- 0
    # for ( this.petsys.number in unique( petsys.numbers ) ) {
    #   these.channels <- petsys.ch[petsys.numbers==this.petsys.number]
    #   if ( min( abs( these.channels - 87 ) ) == 0 ) {
    #     if ( min( abs( these.channels - 81 ) ) == 0 ) {
    #       coinc.81.87 <- coinc.81.87 + 1
    #     }
    #     if ( min( abs( these.channels - 84 ) ) == 0 ) {
    #       coinc.84.87 <- coinc.84.87 + 1
    #     }
    #     # print( c( coinc.81.87, coinc.84.87 ) )
    #   }
    # }
    
    # Set time
    if ( substr( filename, 1, 4 ) == "KNMC" ) {
      this.file <- paste0( program.root, data.root, "/", filename, "_runinfo.txt" )
      source( this.file )
      exposure.time <- num_total_muons  / ((x0_max - x0_min)*(y0_max - y0_min)) / 0.0109
    } else {
      exposure.time <- max( petsys.time ) # in ps
    }
    
    source( paste0( program.root, "/PETsys Cluster Coinc Data.R" ) )
    
    # Get hitstats and eventstats files, this will skip if the hitstats file exists
    
    this.file <- paste0( program.root, stats.root, "/", filename, "-hitstats.dat" )
    print( paste0( "Generating ", filename, "-hitstats.dat" ) )
    source( paste0( program.root, "/PETsys Get Hitstats.R" ) )
    source( paste0( program.root, "/PETsys Load Hitstats File.R" ) ) # All this does is change names to coinc.XXX
    source( paste0( program.root, "/PETsys Get Eventstats.R" ) )
    
  }
}

