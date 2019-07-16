rm(list = ls())
library(tidyverse)
library(ebvSim)

label <- "ebvSim.snps_google_drive_7"

ne <- analyzeReps("ldNe", label)
g2 <- analyzeReps("g2", label)
het <- analyzeReps("het", label)
froh <- analyzeReps("froh", label)

all.ebvs <- ne %>% 
  left_join(g2, by = c("scenario", "replicate", "stratum")) %>% 
  left_join(het, by = c("scenario", "replicate", "stratum")) %>% 
  left_join(froh, by = c("scenario", "replicate", "stratum"))
