#' @title Analyze Replicates
#' @description Conduct analysis of all replicates for a given EBV
#'
#' @param folder folder containing replicate .csv files.
#' @param scenarios vector of scenario numbers to analyze.
#' @param analyses named list of EBV analyses to run.   
#' @param num.cores number of cores to use to split replicate analyses among.
#' 
#' @return list of data frames of summary metrics
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'  
analyzeReps <- function(folder, scenarios = NULL, 
                        analyses = NULL, num.cores = 1) {
  if(is.null(analyses)) {
    analyses <- list(
      ldNe = calc_ldNe,
      g2 = calc_g2,
      het = calc_heterozygosity,
      froh = calc_froh
    )
  }
  if(is.null(names(analyses))) names(analyses) <- 1:length(analyses)
  sc.df <- repCSVfiles(folder) 
  if(is.null(scenarios)) scenarios <- unique(sc.df$scenario)
  sc.df <- sc.df[sc.df$scenario %in% scenarios, ]
  n <- nrow(sc.df)
  
  cat(format(Sys.time()), "Starting analysis of", n, "replicates...\n")
  
  files.written <- if(num.cores == 1) {
    sapply(1:n, .repAnalysis, rep.df = sc.df, analyses = analyses)
  } else {
    cl <- strataG:::.setupClusters(num.cores)
    tryCatch({
      parallel::clusterEvalQ(cl, require(ebvSim))
      parallel::clusterEvalQ(cl, require(strataG))
      parallel::clusterExport(cl, c("sc.df", "analyses"), environment())
      parallel::parSapplyLB(
        cl, 1:n, .repAnalysis, rep.df = sc.df, analyses = analyses
      )
    }, finally = parallel::stopCluster(cl))
  } 
  
  cat(format(Sys.time()), "Analysis complete!\n")
  
  sc.df$files.written <- files.written
  sc.df
}

#' @noRd
#' 
.repAnalysis <- function(i, rep.df, analyses) {
  sc <- rep.df$scenario[i]
  rep <- rep.df$replicate[i]
  
  out.folder <- "analysis.results"
  if(!dir.exists(out.folder)) dir.create(out.folder)
  out.file <- file.path(
    out.folder,
    paste0("scenario.", sc, "_replicate.", rep, "_analysis.results.csv")
  )
  if(file.exists(out.file)) return(FALSE)
                     
  cat(format(Sys.time()), "Scenario", sc, "/ Replicate", rep, "\n")
  
  f.df <- data.table::fread(file = rep.df$fname[i])
  results <- sapply(names(analyses), function(x) {
    cat(format(Sys.time()), "...", x, "\n")
    suppressMessages(results <- analyses[[x]](f.df))
    if(is.null(results)) return(NULL)
    results %>% 
      dplyr::mutate(scenario = sc, replicate = rep) %>% 
      dplyr::select(
        .data$scenario, .data$replicate, .data$stratum, dplyr::everything()
      )
  }, simplify = FALSE) 
  
  results <- results[!sapply(results, is.null)]
  if(length(results) > 1) {
    results <- purrr::reduce(
      results,
      dplyr::full_join, by = c("scenario", "replicate", "stratum")
    )
  }
  
  utils::write.csv(results, file = out.file, row.names = FALSE)
  return(TRUE)
}
