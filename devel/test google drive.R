rm(list = ls())
library(ebvSim)

# get run labels from google drive
run.labels <- availRuns(google.drive.id)
run.labels

# download a run from google drive
dl.dir <- downloadRun(run.labels[1], google.drive.id, "dl.run")
dl.dir
dir(dl.dir)