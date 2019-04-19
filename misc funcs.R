ebvSim <- function(genetics, label, num.sim, ploidy = 2) {
  p <- fscWrite(
    demes = fscSettingsDemes(
      fscDeme(500, 100),
      fscDeme(500, 100),
      fscDeme(500, 100), 
      ploidy = ploidy
    ),
    migration = fscSettingsMigration(
      matrix(c(0, 0.5, 0.005, 0.5, 0, 0.0005, 0.005, 0.0005, 0), ncol = 3)
    ),
    events = fscSettingsEvents(fscEvent(2000, 2, 1), fscEvent(2000, 1, 0)),
    genetics = genetics,
    label = label
  )
  fscRun(p, num.sim = num.sim)
}


killExcess <- function(rl, n) {
  to.kill <- tapply(1:nrow(rl$individuals), rl$individuals[, 1], function(i) {
    if (length(i) > n) {
      sample(i, length(i) - n)
    } else NULL
  })
  to.kill <- unlist(unname(to.kill))
  if (length(to.kill) > 0) rl$individuals <- rl$individuals[-to.kill, ]
  rl
}


loadLandscape <- function(sc, AlleleFreqs, num.gens) {
  rl <- rmetasim::landscape.new.intparam(
    rmetasim::landscape.new.empty(), h = sc$num.pops, 
    s = 2, cg = 0, ce = 0, totgen = num.gens + 1
  )
  
  rl <- rmetasim::landscape.new.floatparam(rmetasim::landscape.new.switchparam(rl))
  
  localS <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2)
  localR <- matrix(c(0, 0, 1.2, 0), nrow = 2, ncol = 2)
  localM <- matrix(c(0, 0, 0, 1.2), nrow = 2, ncol = 2)
  rl <- rmetasim::landscape.new.local.demo(rl, localS, localR, localM)
  
  S <- M <- matrix(0, nrow = sc$num.pops * 2, ncol = sc$num.pops * 2)
  diag(S) <- diag(M) <- 1
  R <- rmetasim::landscape.mig.matrix(
    h = nrow(sc$mig.mat), s = 2, mig.model = "custom", R.custom = sc$mig.mat
  )$R
  rl <- rmetasim::landscape.new.epoch(rl, R = R, carry = rep(sc$ne, sc$num.pops))
  
  #just make loci that have the correct type and ploidy
  for (i in 1:length(AlleleFreqs)) {
    rl <- rmetasim::landscape.new.locus(
      rl, type = 2, ploidy = 2, mutationrate = 0, #sc$mut.rate,
      transmission = 0, numalleles = 2, states = NULL
    )
  }
  
  landscape.new.ind.genos(rl, rep(c(sc$ne, 0), sc$num.pops), AlleleFreqs)
}