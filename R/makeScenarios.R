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
#' @param marker.type type of marker to simulate ("snp").
#' @param mut.rate mutation rate of marker in units of # of mutations per 
#'   generation.
#' @param num.loci number of independent loci.
#' @param ploidy ploidy of marker. Set to 2 for diploid.
#' @param rmetasim.ngen number of generations to run Rmetasim for. Set to 
#'   \code{0} to skip Rmetasim.
#' 
#' @note Each parameter can be a single value or vector of values. Output will 
#'   be all unique combinations of values, with the following exceptions. 
#'   Scenarios wth \code{num.samples} > \code{Ne} are deleted. Scenarios with
#'   \code{num.pops = 1} will have \code{mig.rate} set to NA. 
#' 
#' @return data frame of scenario specifications. 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'  
makeScenarios <- function(num.pops, Ne, num.samples, mig.rate, mig.type,
                          dvgnc.time, marker.type, mut.rate, num.loci, 
                          ploidy, rmetasim.ngen) {
  expand.grid(
    num.pops = num.pops,
    Ne = Ne,
    num.samples = num.samples,
    mig.rate = mig.rate,
    mig.type = tolower(mig.type),
    dvgnc.time = dvgnc.time,
    marker.type = tolower(marker.type),
    mut.rate = mut.rate,
    num.loci = num.loci,
    ploidy = ploidy,
    rmetasim.ngen = rmetasim.ngen,
    stringsAsFactors = FALSE
  ) %>% 
    dplyr::mutate(
      mig.rate = ifelse(.data$num.pops == 1, NA, .data$mig.rate),
      num.samples = ifelse(
        is.na(.data$num.samples), 
        NA, 
        ifelse(.data$Ne <= .data$num.samples, NA, .data$num.samples)
      )
    ) %>% 
    dplyr::filter(!duplicated(.data)) %>% 
    dplyr::mutate(scenario = 1:dplyr::n()) %>% 
    dplyr::select(.data$scenario, dplyr::everything())
}
