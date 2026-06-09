# R equivalent for reading the packed binary format
# A few important notes:
#   The C structs are marked __packed__, so there is no padding between fields.
# This code assumes:
#   little-endian binary layout
# long long = 8-byte signed integer
# float = IEEE 754 4-byte float
# int = 4-byte signed integer
# R’s readBin(..., size=8, what="integer") may return doubles for large 64-bit integers because base R does not have native 64-bit integer storage. If exact 64-bit handling matters, use the bit64 package.

# library(readBin) # Preloaded I guess

# Define file path
con <- file("/Users/ifft/Desktop/Muography/PETsys Analysis/PETsys Data/KNVA-20260507-03-00084_coinc.ldat", "rb")

# Size of packed C struct Event:
# long long = 8 bytes
# float     = 4 bytes
# int       = 4 bytes
# total     = 16 bytes

read_event <- function(con) {
  list(
    time = readBin(con, what = "integer", size = 8, n = 1,
                   signed = TRUE, endian = "big"),
    e    = readBin(con, what = "numeric", size = 4, n = 1,
                   endian = "big"),
    id   = readBin(con, what = "integer", size = 4, n = 1,
                   signed = TRUE, endian = "big")
  )
}

read_event_array <- function(con, n) {
  if (n == 0) {
    return(data.frame(
      time = numeric(0),
      e = numeric(0),
      id = integer(0)
    ))
  }
  
  data.frame(
    time = readBin(con, what = "integer", size = 8, n = n,
                   signed = TRUE, endian = "big"),
    e    = readBin(con, what = "numeric", size = 4, n = n,
                   endian = "big"),
    id   = readBin(con, what = "integer", size = 4, n = n,
                   signed = TRUE, endian = "big")
  )
}

results <- list()
i <- 1

#repeat {
  
  # Stop if end of file reached
  if (isIncomplete(con)) break
  
  # Read CoincidenceGroupHeader
  # uint8_t nHits1;
  # uint8_t nHits2;
  
  hdr <- tryCatch({
    list(
      nHits1 = readBin(con, what = "integer", size = 1,
                       signed = FALSE, n = 1),
      nHits2 = readBin(con, what = "integer", size = 1,
                       signed = FALSE, n = 1)
    )
  }, error = function(e) NULL)
  
  if (is.null(hdr)) break
  
  # Read event arrays
  coinc1 <- read_event_array(con, hdr$nHits1)
  coinc2 <- read_event_array(con, hdr$nHits2)
  
  results[[i]] <- list(
    header = hdr,
    coinc1 = coinc1,
    coinc2 = coinc2
  )
  
  i <- i + 1
#}

close(con)

install.packages("bit64")
library(bit64)
large_int <- as.integer64(coinc1$time[1])
print(coinc1$time[1], scientific = TRUE)
# [1] 9223372036854775807


large_int <- as.integer64("9223372036854775807")
print(large_int)
# [1] 9223372036854775807
