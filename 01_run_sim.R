rm(list = ls())

# run next two lines if you don't have the latest devel version of strataG
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("ericarcher/strataG", ref = "redevel2018", dependencies = TRUE)
stopifnot(require(tidyverse))
stopifnot(require(strataG))
source("misc_funcs.R")


label <- "ebvSim.snps2"
google.drive.id <- "1TGI2TVFnOAx0ib1GdBaL80Pwq-7ruqG1"

# set scenarios
scenarios <- expand.grid(
  num.pops = c(1, 5),
  Ne = c(100, 1000),
  num.samples = c(10, 100),
  mig.rate = 1e-5,
  mig.type = "island",
  dvgnc.time = 200,
  stringsAsFactors = FALSE
) %>% 
  mutate(scenario = 1:n()) %>% 
  select(scenario, everything())

# specify SNP loci
genetics <- fscSettingsGenetics(fscBlock_snp(1, 1e-4), num.chrom = 1000)

# run fastsimcoal2
out.dir <- ebvSim(
  scenarios = scenarios, 
  genetics = genetics, 
  label = label, 
  num.sim = 3,
  google.drive.id = google.drive.id
)
save.image("ebvSim ws.rdata")

# get run labels from google drive
run.labels <- availRuns(google.drive.id)
run.labels

# download a run from google drive
dl.dir <- downloadRun(run.labels[1], google.drive.id, "dl.run")
dl.dir
dir(dl.dir)


# run rmetasim to establish linkage disequilibrium
# num.rms.gens <- 5
#af <- alleleFreqs(sim.snps.g, by.strata = T, type = "prop")
# rl <- loadLandscape(sc, af, num.rms.gens)
# for(i in 1:num.rms.gens) {
#   rl <- rmetasim::landscape.simulate(rl, 1)
#   rl <- killExcess(rl, sc$ne)
# }
# rl.g <- landscape2gtypes(rl)

