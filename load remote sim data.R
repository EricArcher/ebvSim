rm(list = ls())
library(googledrive)

drive.id <- "185soK9wFB9G7NS4bcAA-y0FxGOxPpvHQ"
x <- drive_ls(as_id(drive.id))



temp <- tempfile(fileext = ".csv")
dl <- drive_download(
  as_id("1AiZda_1-2nwrxI8fLD0Y6e5rTg7aocv0"), path = temp, overwrite = TRUE)
out <- unzip(temp, exdir = tempdir())
bank <- read.csv(out[14], sep = ";")