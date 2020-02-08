#' @noRd
#' 
.runFscSim <- function(rep.i, params) {
  sc.i <- params$rep.df$sc[rep.i]
  sc <- params$scenarios[sc.i, ]
  
  deme.list <- lapply(1:sc$num.pops, function(i) {
    strataG::fscDeme(deme.size = sc$Ne, sample.size = sc$num.samples)
  })
  deme.list$ploidy <- sc$ploidy
  
  genetics <- switch(
    as.character(sc$marker.type),
    snp = strataG::fscSettingsGenetics(
      strataG::fscBlock_snp(1, sc$mut.rate), 
      num.chrom = sc$num.loci
    )
  )
  
  label <- paste0(
    params$label, ".sc_", sc$scenario, ".rep_", params$rep.df$rep[rep.i]
  )
  strataG::fscWrite(
    demes = do.call(strataG::fscSettingsDemes, deme.list),
    migration = if(sc$num.pops > 1) {
      mig.mat <- makeMigMat(sc$mig.rate, sc$num.pops, sc$mig.type) 
      strataG::fscSettingsMigration(mig.mat)
    } else NULL,
    events = makeEventSettings(sc$dvgnc.time, sc$num.pops),
    genetics = genetics,
    label = label,
    use.wd = params$use.wd
  ) %>% 
    strataG::fscRun(num.cores = 1, exec = params$fsc.exec)
}

#' @noRd
#' 
.setupScRland <- function(sc.num, params) {  
  if(params$scenarios$rmetasim.ngen[sc.num] == 0) return(NA)
  sc <- params$scenarios[sc.num, ]
  localS <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2)
  localR <- matrix(c(0, 0, 1, 0), nrow = 2, ncol = 2)
  localM <- matrix(c(0, 0, 0, 1), nrow = 2, ncol = 2)
  S <- M <- matrix(0, nrow = sc$num.pops * 2, ncol = sc$num.pops * 2)
  diag(S) <- diag(M) <- 1
  R <- if(sc$num.pops == 1) NULL else {
    rmetasim::landscape.mig.matrix(
      h = sc$num.pops, s = 2, mig.model = "custom", 
      R.custom = makeMigMat(sc$mig.rate, sc$num.pops, sc$mig.type)
    )$R
  }
  
  rmetasim::landscape.new.empty() %>% 
    rmetasim::landscape.new.intparam(
      h = sc$num.pops, s = 2, cg = 0, ce = 0, totgen = sc$rmetasim.ngen + 1
    ) %>% 
    rmetasim::landscape.new.switchparam() %>% 
    rmetasim::landscape.new.floatparam() %>% 
    rmetasim::landscape.new.local.demo(localS, localR, localM) %>% 
    rmetasim::landscape.new.epoch(R = R, carry = rep(sc$Ne, sc$num.pops)) 
}

#' @noRd
#' 
.runRmetasim <- function(freqs, Rland, sc) {
  for(i in 1:length(freqs$global)) {
    Rland <- rmetasim::landscape.new.locus(
      Rland, type = 2, ploidy = sc$ploidy, mutationrate = 0,
      transmission = 0, numalleles = 2, allelesize = 1,
      frequencies = freqs$global[[i]], states = names(freqs$global[[i]])
    )
  }
  
  Rland %>% 
    rmetasim::landscape.new.individuals(rep(c(sc$Ne, 0), sc$num.pops)) %>% 
    rmetasim::landscape.setpopfreq(freqs$pop) %>% 
    rmetasim::landscape.simulate(numit = sc$rmetasim.ngen)
}

#' @noRd
#' 
.killExcess <- function(rl, n) {
  to.keep <- tapply(1:nrow(rl$individuals), rl$individuals[, 1], function(i) {
    if(length(i) > n) sample(i, n) else i
  })
  to.keep <- sort(unlist(unname(to.keep)))
  rl$individuals <- rl$individuals[to.keep, , drop = FALSE]
  rl
}
