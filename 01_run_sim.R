rm(list = ls())

# run next two lines if you don't have the latest devel version of strataG
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("ericarcher/strataG", ref = "redevel2018", dependencies = TRUE)

library(strataG)
source("misc_funcs.R")

scenarios <- expand.grid(
  num.pops = c(1, 5),
  Ne = c(100, 1000),
  num.samples = c(10, 100),
  mig.rate = 1e-5,
  mig.type = "island",
  dvgnc.time = 200,
  stringsAsFactors = FALSE
)
scenarios$scenario <- 1:nrow(scenarios)

genetics <- fscSettingsGenetics(fscBlock_snp(1, 1e-4), num.chrom = 1000)

# run fastsimcoal2
ebvSim(scenarios, genetics, label = "ebvSim.snps", num.sim = 3)

save.image("ebvSim ws.rdata")





# convert SNPs to gtypes
# sim.snps.g <- df2gtypes(sim.snps, ploidy = 2)

# run rmetasim to establish linkage disequilibrium
# num.rms.gens <- 5
#af <- alleleFreqs(sim.snps.g, by.strata = T, type = "prop")
# rl <- loadLandscape(sc, af, num.rms.gens)
# for(i in 1:num.rms.gens) {
#   rl <- rmetasim::landscape.simulate(rl, 1)
#   rl <- killExcess(rl, sc$ne)
# }
# rl.g <- landscape2gtypes(rl)

