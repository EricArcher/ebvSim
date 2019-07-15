rm(list = ls())
library(tidyverse)
library(ebvSim)

label <- "ebvSim.snps_freq_test"
num.reps <- 3

scenarios <- loadScenarios(label)

replicates <- expand.grid(scenario = 1:nrow(scenarios), replicate = 1:num.reps)
n <- nrow(replicates)

ne.est <- do.call(
  rbind,
  lapply(1:n, function(i) {
    cat(i, "/", n, "\n")
    sc <- replicates$scenario[i]
    rep <- replicates$replicate[i]
    loadGenotypes(label, sc, rep) %>% 
      calc_ldNe() %>% 
      mutate(
        scenario = sc,
        replicate = rep
      )
        
  })
)
        