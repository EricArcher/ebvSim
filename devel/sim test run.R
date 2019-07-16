# run next two lines if you don't have the latest devel version of strataG
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("ericarcher/strataG", ref = "redevel2018", dependencies = TRUE)

# you will also need the latest version of rmetasim from GitHub.
# uncomment the next two lines and run to install it
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("stranda/rmetasim", dependencies = TRUE)

rm(list = ls())
library(ebvSim)
library(strataG)

label <- "ebvSim.snps_test_1"
num.rep <- 3
google.drive.id <- "1TGI2TVFnOAx0ib1GdBaL80Pwq-7ruqG1"

# set scenarios
scenarios <- makeScenarios(
  num.pops = c(1, 5),
  Ne = c(100, 1000),
  num.samples = c(10, 100),
  mig.rate = 1e-5,
  mig.type = "island",
  dvgnc.time = 200
)

# run fastsimcoal2
out.dir <- runEBVsim(
  label = label,
  scenarios = scenarios,
  genetics = fscSettingsGenetics(fscBlock_snp(1, 1e-4), num.chrom = 1000),
  num.rep = num.rep,
  google.drive.id = NULL, #google.drive.id,
  run.rmetasim = TRUE
)

source("analysis test.R")
