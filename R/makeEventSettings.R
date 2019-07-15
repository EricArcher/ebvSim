#' @title Make Event Settings
#' @description Create event input for strataG::fscRun.
#' 
#' @param dvgnc.time number of generations since divergence.
#' @param num.pops number of populations.
#' 
#' @return fscEvents object output from strataG::fscSettingsEvents.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'  
makeEventSettings <- function(dvgnc.time, num.pops) {
  if(num.pops == 1) return(NULL)
  pop.pairs <- t(utils::combn(num.pops, 2) - 1)
  pop.pairs <- pop.pairs[pop.pairs[, 1] == 0, , drop = FALSE]
  do.call(
    strataG::fscSettingsEvents, 
    lapply(
      1:nrow(pop.pairs),
      function(i) {
        strataG::fscEvent(dvgnc.time, pop.pairs[i, 2], pop.pairs[i, 1])
      }
    )
  )
}
