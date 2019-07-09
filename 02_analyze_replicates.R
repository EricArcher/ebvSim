rm(list = ls())
library(strataG)
library(tidyverse)
load("ebvSim ws.rdata")

source("02a.heterozygosity.R")
source("02b.ldNe.R")
source("02c.g2.R")
source("02d.froh.R")

cat(format(Sys.time()), "Analyses done!\n")
save.image("ebv analysis ws.rdata")