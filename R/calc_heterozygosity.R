#' @title Calculate Heterozygosity
#' @description Calculate observed and expected heterozygosity
#' 
#' @param x data frame of genotypes as read from simulation replicate .csv file
#' 
#' @return data frame of mean, median and 95% CI of observed and expected 
#'   heterozygosity for each population
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'
calc_heterozygosity <- function(x) {
  sim.snps.g <- strataG::df2gtypes(x, ploidy = 2) 
  
  dplyr::left_join(
    strataG::heterozygosity(sim.snps.g, by.strata = TRUE, type = "observed"),
    strataG::heterozygosity(sim.snps.g, by.strata = TRUE, type = "expected"),
    by = c("stratum", "locus")
  ) %>% 
    dplyr::group_by(.data$stratum) %>% 
    dplyr::summarize(
      obs.het.mean = mean(.data$obsvd.het),
      obs.het.sd = stats::sd(.data$obsvd.het),
      obs.het.median = stats::median(.data$obsvd.het),
      obs.het.lci = unname(stats::quantile(.data$obsvd.het, 0.025)),
      obs.het.uci = unname(stats::quantile(.data$obsvd.het, 0.975)),
      exp.het.mean = mean(.data$exptd.het),
      exp.het.sd = stats::sd(.data$exptd.het),
      exp.het.median = stats::median(.data$exptd.het),
      exp.het.lci = unname(stats::quantile(.data$exptd.het, 0.025)),
      exp.het.uci = unname(stats::quantile(.data$exptd.het, 0.975))
    ) %>% 
    dplyr::ungroup() %>% 
    as.data.frame(stringsAsFactors = FALSE)
}