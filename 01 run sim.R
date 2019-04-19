rm(list = ls())
# run next two lines if you don't have the latest devel version of strataG
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("ericarcher/strataG", ref = "redevel2018", dependencies = TRUE)
library(strataG)
source("misc funcs.R")

# run fastsimcoal2
snp.p <- ebvSim(
  fscSettingsGenetics(fscBlock_snp(1, 1e-5), num.chrom = 1000),
  label = "ebvSim.snps",
  num.sim = 1
)

# example DNA and microsat blocks
# dna.msat.p <- ebvSim(
#   fscSettingsGenetics(fscBlock_dna(1000, 1e-6), fscBlock_microsat(15, 5e-4)),
#   label = "ebvSim.dna.msat",
#   num.sim = 1
# )

# convert SNPs to gtypes
sim.snps <- fscReadArp(snp.p, drop.mono = TRUE)
sim.snps.g <- df2gtypes(sim.snps, ploidy = 2)

# run rmetasim to establish linkage disequilibrium
num.rms.gens <- 5
af <- alleleFreqs(sim.snps.g, by.strata = T)
rl <- loadLandscape(sc, af, num.rms.gens)
for(i in 1:num.rms.gens) {
  rl <- rmetasim::landscape.simulate(rl, 1)
  rl <- killExcess(rl, sc$ne)
}
rl.g <- landscape2gtypes(rl)

save.image("ebvSim ws.rdata")