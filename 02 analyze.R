rm(list = ls())
library(strataG)
load("ebvSim ws.rdata")

sim.snps <- fscReadArp(snp.p)
sim.snps.g <- df2gtypes(sim.snps, ploidy = 2)

source("02.locus.summary.R")
source("02.ldNe.R")
source("02.g2.R")
source("02.froh.R")

save.image("ebv analysis ws.rdata")