rm(list = ls())
library(strataG)
library(ebvSim)
library(tidyverse)
library(rmetasim)
source("landscape funcs.R")

Ne <- 1000
sample.size <- 10
num.pops <- 10
mig.rate <- 1e-5
mig.type <- "island"
mig.mat <- makeMigMat(1e-5, num.pops, "island")

rmetasim.ngen <- 10

demes <- lapply(1:num.pops, function(i) {
  fscDeme(deme.size = Ne, sample.size = sample.size)
})
demes$ploidy <- 2

p <- strataG::fscWrite(
  demes = do.call(fscSettingsDemes, demes),
  migration = fscSettingsMigration(mig.mat),
  events = makeEventSettings(20000, num.pops),
  genetics = fscSettingsGenetics(fscBlock_snp(1, 1e-5), num.chrom = 1000),
  label = "rmetasim.test"
) %>% 
  strataG::fscRun(num.sims = 3)

gen.data <- fscReadArp(p, sim = c(1, 1))

freqs <- calcFreqs(gen.data)

npop <- length(freqs$pop)
localS <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2)
localR <- matrix(c(0, 0, 1.2, 0), nrow = 2, ncol = 2)
localM <- matrix(c(0, 0, 0, 1.2), nrow = 2, ncol = 2)
S <- M <- matrix(0, nrow = npop * 2, ncol = npop * 2)
diag(S) <- diag(M) <- 1
R <- rmetasim::landscape.mig.matrix(
  h = nrow(mig.mat), 
  s = 2, 
  mig.model = "custom", 
  R.custom = mig.mat
)$R

Rland <- rmetasim::landscape.new.empty() %>% 
  rmetasim::landscape.new.intparam(
    h = num.pops, s = 2, cg = 0, ce = 0, totgen = rmetasim.ngen + 1
  ) %>% 
  rmetasim::landscape.new.switchparam() %>% 
  rmetasim::landscape.new.floatparam() %>% 
  rmetasim::landscape.new.local.demo(localS, localR, localM) %>% 
  rmetasim::landscape.new.epoch(R = R, carry = rep(Ne, npop)) 

for(i in 1:length(freqs$global)) {
  Rland <- rmetasim::landscape.new.locus(
    Rland, type = 2, ploidy = 2, mutationrate = 0,
    transmission = 0, numalleles = 2, allelesize = 1,
    frequencies = freqs$global[[i]], states = names(freqs$global[[i]])
  )
}

Rland <- rmetasim::landscape.new.individuals(Rland, rep(c(Ne, 0), npop))
Rland <- rmetasim::landscape.setpopfreq(Rland, freqs$pop)

for(i in 1:rmetasim.ngen) {
  Rland <- rmetasim::landscape.simulate(Rland, numit = 1)
  Rland <- killExcess(Rland, Ne)
}


  rland.gi <- rmetasim::landscape.make.genind(Rland)
  rland.g <- strataG::genind2gtypes(rland.gi)
  g.ne <- ldNe(rland.g)
  print(g.ne)
mean(g.ne$Ne)



