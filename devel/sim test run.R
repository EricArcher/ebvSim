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

label <- "ebvSim.snps_google_drive_1"
num.rep <- 10
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

# specify SNP loci
genetics <- fscSettingsGenetics(fscBlock_snp(1, 1e-4), num.chrom = 1000)

# run fastsimcoal2
out.dir <- runEBVsim(
  label = label,
  scenarios = scenarios,
  genetics = genetics,
  num.rep = num.rep,
  google.drive.id = google.drive.id,
  run.rmetasim = TRUE
)

save.image(paste0(label, "_ws.rdata"))
