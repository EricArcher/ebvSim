rm(list = ls())
library(ebvSim)
library(strataG)

scenarios <- makeScenarios(
  num.pops = 3,
  Ne = 100,
  num.samples = 10,
  mig.rate = 1e-5,
  mig.type = "island",
  dvgnc.time = 100, 
  rmetasim.ngen = 10
)

output <- runEBVsim(
  label = "ebvSim.snps_test",
  scenarios = scenarios,
  genetics = fscSettingsGenetics(fscBlock_snp(1, 1e-4), num.chrom = 1000),
  ploidy = 2,
  num.rep = 10,
  num.cores = 4
)