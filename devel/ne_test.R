rm(list = ls())
library(ebvSim)
library(strataG)
library(tidyverse)

label <- "ne_test2"
num.rep <- 5

# set scenarios
scenarios <- makeScenarios(
  num.pops = 5,
  Ne = c(100, 1000),
  num.samples = c(10, 100, 1000),
  mig.rate = 1e-5,
  mig.type = "island",
  dvgnc.time = 100,
  rmetasim.ngen = 10
)

runEBVsim(
  label = label,
  scenarios = scenarios,
  genetics = fscSettingsGenetics(fscBlock_snp(1, 1e-4), num.chrom = 1000),
  num.rep = num.rep,
  google.drive.id = NULL
)
ne <- analyzeReps("ldNe", label)
ne

left_join(
  scenarios, 
  ne %>% 
    group_by(scenario) %>% 
    summarize(median.ne = median(Ne[!is.infinite(Ne)])),
  by = "scenario"
)

save.image("ne test.rdata")