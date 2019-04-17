rm(list = ls())
library(strataG)
library(tidyverse)
load("ebvSim ws.rdata")

sim.snps <- fscReadArp(snp.p, drop.mono = TRUE)
sim.snps.g <- df2gtypes(sim.snps, ploidy = 2)

source("02.heterozygosity.R")
source("02.ldNe.R")
source("02.g2.R")
source("02.froh.R")

cat(format(Sys.time()), "Analyses done!\n")
save.image("ebv analysis ws.rdata")