#' @title Make Scenario Parameters
#' @description Create data frame of scenario parameters.
#' 
#' @param num.pops number of populations.
#' @param Ne effective population size.
#' @param num.samples number of samples to simulate. Must be <= Ne.
#' @param mig.rate migration rate specified as proportion of population 
#'   migrating per generaton.
#' @param mig.type of migration matrix structure. "island" = rate between all 
#'   populations is the same. "stepping.stone" = migration only occurs between 
#'   neighboring populations. Metapopulation is ring shaped, not a linear chain.
#' @param dvgnc.time number of generations since divergence.
#' 
#' @note Each parameter can be a single value or vector of values. Output will 
#'   be all unique combinations of values.
#' 
#' @return data frame of scenario specifications. 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'  
makeScenarios <- function(num.pops, Ne, num.samples, mig.rate, 
                          mig.type = c("island", "stepping.stone"),
                          dvgnc.time) {
  expand.grid(
    num.pops = num.pops,
    Ne = Ne,
    num.samples = num.samples,
    mig.rate = mig.rate,
    mig.type = match.arg(mig.type),
    dvgnc.time = dvgnc.time,
    stringsAsFactors = FALSE
  ) %>% 
    dplyr::mutate(scenario = 1:dplyr::n()) %>% 
    dplyr::select(.data$scenario, dplyr::everything())
}
