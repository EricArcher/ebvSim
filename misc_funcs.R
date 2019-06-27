makeMigSettings <- function(mig.rate, num.pops, 
                        type = c("island", "stepping.stone")) {
  type <- match.arg(type)
  switch(
    type,      
    island = {
      m <- mig.rate / (num.pops - 1)
      mat <- matrix(rep(m, num.pops ^ 2), nrow = num.pops)
      diag(mat) <- 1 - mig.rate
      fscSettingsMigration(mat)
    },
    stepping.stone = {
      mat <- matrix(0, nrow = num.pops, ncol = num.pops)
      m <- mig.rate / 2
      for (k in 1:(num.pops - 1)) {
        mat[k, k + 1] <- mat[k + 1, k] <- m
      }
      mat[1, num.pops] <- mat[num.pops, 1] <- m
      diag(mat) <- 1 - mig.rate
      fscSettingsMigration(mat)
    }
  )
}


makeEventSettings <- function(dvgnc.time, num.pops) {
  if(num.pops == 1) return(NULL)
  pop.pairs <- t(combn(num.pops, 2) - 1)
  pop.pairs <- pop.pairs[pop.pairs[, 1] == 0, , drop = FALSE]
  do.call(
    fscSettingsEvents, 
    lapply(
      1:nrow(pop.pairs),
      function(i) fscEvent(dvgnc.time, pop.pairs[i, 2], pop.pairs[i, 1])
    )
  )
}


ebvSim <- function(scenarios, genetics, label, num.sim, ploidy = 2) {
  out.folder <- paste0(label, ".scenario.replicates")
  if(!dir.exists(out.folder)) dir.create(out.folder)
  fname <- paste0(label, "_scenarios.csv")
  write.csv(scenarios, file = file.path(out.folder, fname), row.names = FALSE)
  
  for(sc.i in 1:nrow(scenarios)) {
    sc <- scenarios[sc.i, ]
    
    deme.list <- lapply(1:sc$num.pops, function(i) {
      fscDeme(deme.size = sc$Ne, sample.size = sc$num.samples)
    })
    deme.list$ploidy <- ploidy
    
    p <- fscWrite(
      demes = do.call(fscSettingsDemes, deme.list),
      migration = if(sc$num.pops > 1) {
        makeMigSettings(sc$mig.rate, sc$num.pops, sc$mig.type) 
      } else NULL,
      events = makeEventSettings(sc$dvgnc.time, sc$num.pops),
      genetics = genetics,
      label = label
    )
    
    p <- fscRun(p, num.sim = num.sim)
    
    for(sim.i in 1:num.sim) {
      gen.data <- fscReadArp(p, sim = c(1, sim.i), drop.mono = TRUE)
      fname <- paste0(
        label, "_", 
        "scenario_", sc$scenario, 
        "_replicate_", sim.i, 
        ".csv"
      )
      write.csv(gen.data, file = file.path(out.folder, fname), row.names = FALSE)
    }
  }
  
  invisible(out.folder)
}


killExcess <- function(Rland, n) {
  library(magrittr)
  to.kill <- tapply(
    1:nrow(Rland$individuals), 
    Rland$individuals[, 1], 
    function(i) if(length(i) > n) sample(i, length(i) - n) else NULL
  ) %>% 
    unname() %>% 
    unlist()
  if(length(to.kill) > 0) Rland$individuals <- Rland$individuals[-to.kill, ]
  Rland
}


loadLandscape <- function(sc, AlleleFreqs, num.gens) {
  library(magrittr)
  library(rmetasim)
  
  localS <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2)
  localR <- matrix(c(0, 0, 1.2, 0), nrow = 2, ncol = 2)
  localM <- matrix(c(0, 0, 0, 1.2), nrow = 2, ncol = 2)
  S <- M <- matrix(0, nrow = sc$num.pops * 2, ncol = sc$num.pops * 2)
  diag(S) <- diag(M) <- 1
  R <- rmetasim::landscape.mig.matrix(
    h = nrow(sc$mig.mat), s = 2, mig.model = "custom", R.custom = sc$mig.mat
  )$R
  
  Rland <- rmetasim::landscape.new.empty() %>% 
    rmetasim::landscape.new.intparam(
      h = sc$num.pops, s = 2, cg = 0, ce = 0, totgen = num.gens + 1
    ) %>% 
    rmetasim::landscape.new.switchparam() %>% 
    rmetasim::landscape.new.floatparam() %>% 
    rmetasim::landscape.new.local.demo(localS, localR, localM) %>% 
    rmetasim::landscape.new.epoch(R = R, carry = rep(sc$ne, sc$num.pops))
  
  # just make loci that have the correct type and ploidy
  for(i in 1:length(AlleleFreqs)) {
    Rland <- rmetasim::landscape.new.locus(
      Rland, type = 2, ploidy = 2, mutationrate = 0,
      transmission = 0, numalleles = 2, allelesize = 1,
      frequencies = AlleleFreqs[[i]], states = rownames(AlleleFreqs[[i]])
    )
  }

  rmetasim::landscape.new.individuals(Rland, rep(c(sc$ne, 0), sc$num.pops))
}