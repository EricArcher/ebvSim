#' @title Analyze Replicates
#' @description Conduct analysis of all replicates for a given EBV
#'   
#' @param analysis character vector of EBV analysis to run (`ldNe`, `g2`, 
#'   `het`, or `froh`).   
#' @param params parameter list output from \code{runEBVsim()}.
#' @param num.cores number of cores to use to split replicate analyses among.
#' 
#' @return data frame of summary metrics
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'  
analyzeReps <- function(analysis, params, num.cores = 1) {
    
  params$replicates <- expand.grid(
    scenario = 1:nrow(params$scenarios), 
    replicate = 1:params$num.rep
  )
  params$analysis.func <- switch(
    analysis,
    ldNe = calc_ldNe,
    g2 = calc_g2,
    het = calc_heterozygosity,
    froh = calc_froh
  )
  
  n <- nrow(params$replicates)
  cat(format(Sys.time()), "Starting", analysis, "analysis of", n, "replicates...\n")
  
  result <- if(num.cores == 1) {
    lapply(1:n, function(i, p = params) {
      cat(i, "/", n, "\n")
      .repAnalysis(i, p)
    })
  } else {
    cl <- strataG:::.setupClusters(num.cores)
    tryCatch({
      parallel::clusterEvalQ(cl, require(ebvSim))
      parallel::clusterExport(cl, "params", environment())
      parallel::parLapply(cl, 1:n, .repAnalysis, p = params)
    }, finally = parallel::stopCluster(cl))
  } 

  cat(format(Sys.time()), analysis, "analysis complete!\n")
  
  result %>% 
    dplyr::bind_rows() %>% 
    dplyr::arrange(.data$scenario, .data$replicate, .data$stratum)
}

#' @noRd
#' 
.repAnalysis <- function(i, p) {
  sc <- p$replicates$scenario[i]
  rep <- p$replicates$replicate[i]
  results <- loadGenotypes(p$label, sc, rep) %>% 
    p$analysis.func() 
  if(!is.null(results)) {
    results %>% 
      dplyr::mutate(scenario = sc, replicate = rep) %>% 
      dplyr::select(.data$scenario, .data$replicate, .data$stratum, dplyr::everything())
  } else NULL
}
