#' @title Run fastsimcoal
#' @description Run a fastsimcoal simulation given scenario parameters
#' 
#' @param label label of simulation run.
#' @param sc a one row data frame of scenario parameters.
#' @param genetics specification of genetic data to be simulated as output from 
#'   \code{\link[strataG]{fscSettingsGenetics}}.
#' @param ploidy ploidy of genetic data.
#' @param num.rep number of replicates to run for each scenario.
#' @param use.wd use working directory for fastsimcoal files?
#' @param fsc.exec name of fastsimcoal executable.
#' 
#' @return parameter output from `strataG::fscRun()`
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
runFscSim <- function(label, sc, genetics, ploidy, num.rep, use.wd, fsc.exec) {
  deme.list <- lapply(1:sc$num.pops, function(i) {
    strataG::fscDeme(deme.size = sc$Ne, sample.size = sc$num.samples)
  })
  deme.list$ploidy <- ploidy
  
  strataG::fscWrite(
    demes = do.call(strataG::fscSettingsDemes, deme.list),
    migration = if(sc$num.pops > 1) {
      mig.mat <- makeMigMat(sc$mig.rate, sc$num.pops, sc$mig.type) 
      strataG::fscSettingsMigration(mig.mat)
    } else NULL,
    events = makeEventSettings(sc$dvgnc.time, sc$num.pops),
    genetics = genetics,
    label = label,
    use.wd = use.wd
  ) %>% 
    strataG::fscRun(num.sims = num.rep, num.cores = 1, exec = fsc.exec)
}
