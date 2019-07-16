#' @title Analyze Replicates
#' @description Conduct analysis of all replicates for a given EBV
#'   
#' @param analysis character vector of EBV analysis to run (`ldNe`, `g2`, 
#'   `het`, or `froh`).   
#' @param label label of simulation run.
#' 
#' @return data frame of summary metrics
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'  
analyzeReps <- function(analysis, label) {
  analysis.func <- switch(
    analysis,
    ldNe = calc_ldNe,
    g2 = calc_g2,
    het = calc_heterozygosity,
    froh = calc_froh
  )
  num.rep <- scenarios <- NULL
  load(paste0(label, "_ws.rdata"))
  replicates <- expand.grid(scenario = 1:nrow(scenarios), replicate = 1:num.rep)
  n <- nrow(replicates)
  purrr::map(1:n, function(i) {
    cat(i, "/", n, "\n")
    sc <- replicates$scenario[i]
    rep <- replicates$replicate[i]
    loadGenotypes(label, sc, rep) %>% 
      analysis.func() %>% 
      dplyr::mutate(
        scenario = sc,
        replicate = rep
      ) %>% 
      dplyr::select(.data$scenario, .data$replicate, .data$stratum, dplyr::everything())
    
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::arrange(.data$scenario, .data$replicate, .data$stratum)
}