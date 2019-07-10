rm(list = ls())
library(tidyverse)
library(rmetasim)
source("misc_funcs.R")

label <- "ebvSim.snps2"
scenario <- 1
replicate <- 2
num.gens <- 10
scenarios <- loadScenarios(label)
sc <- scenarios[scenario, ]

g <- loadLandscape(label, sc, 2, num.gens) %>% 
  rmetasim::landscape.make.genind() %>% 
  genind2gtypes()
