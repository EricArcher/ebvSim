#' @title Calculate inbreeding (g2)
#' @description Calculate inbreeding (g2)
#' 
#' @param x data frame of genotypes as read from simulation replicate .csv file.
#' @param nboot number of bootstrap replicates to estimate confidence interval.
#' 
#' @return data frame of g2 and confidence intervals for each population 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'
calc_g2 <- function(x, nboot = 100) {
  split(x, x$stratum) %>% 
    purrr::map(function(x.st) inbreedR::convert_raw(x.st[, -(1:2)])) %>% 
    #purrr::list_modify(global = inbreedR::convert_raw(x[, -(1:2)])) %>% 
    purrr::map(function(x.st) {
      res <- inbreedR::g2_snps(x.st, nboot = nboot, verbose = FALSE)
      data.frame(
        g2 = res$g2, 
        g2.lci = unname(res$CI_boot[1]), 
        g2.uci = unname(res$CI_boot[2])
      )
    }) %>% 
    dplyr::bind_rows() %>% 
    tibble::rownames_to_column("stratum")
}