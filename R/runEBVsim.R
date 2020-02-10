#' @title Run EBV Simulations
#' @description Run a set of scenarios and replicates to generate genotypes.
#' 
#' @param label label used for output files and folders for the run.
#' @param scenarios data.frame where each row provides different parameters for 
#'   each scenario. See \code{Details} for more information on its structure.
#' @param num.rep number of replicates to run for each scenario.
#' @param google.drive.id id of root Google Drive to save run to.
#' @param delete.fsc.files delete fastsimcoal files when the run is complete?
#' @param use.wd use the working directory to save fastsimcoal files? If 
#'   \code{FALSE} (default), a temporary folder is used, which is deleted when 
#'   the R session is closed.
#' @param num.cores number of cores to use. Each replciate is allocated to a 
#'   separate core. If set to 1, then progress is output to the console.
#' @param fsc.exec name of fastsimcoal executable.
#' 
#' @details The \code{scenarios} data.frame must have the following columns:
#' \tabular{ll}{
#'   \code{num.pops} \tab number of populations.\cr
#'   \code{Ne} \tab effective population size of each population.\cr
#'   \code{num.samples} \tab number of samples to use.\cr
#'   \code{mig.rate} \tab migration rate.\cr
#'   \code{mig.type} \tab type of migration matrix: "island" or "stepping.stone".\cr
#'   \code{dvgnc.time} \tab divergence time of populations.\cr
#'   \code{rmetasim.ngen} \tab number of generations to run \code{rmetasim} for.\cr
#'  }
#' 
#' @return vector of folders where scenario data and replicates were written.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'  
runEBVsim <- function(label, scenarios, num.rep,
                      google.drive.id = NULL, delete.fsc.files = TRUE,
                      use.wd = FALSE, num.cores = 1, fsc.exec = "fsc26") {
  
  if(!all(c("num.pops", "Ne", "num.samples", "mig.rate", "mut.rate",
            "num.loci", "ploidy", "rmetasim.ngen", "mig.type",
            "marker.type") %in% colnames(scenarios))) {
    stop("'scenarios' is missing some required columns")
  }
  for(x in colnames(scenarios)) {
    if(x %in% c("mig.type", "marker.type")) {
      scenarios[[x]] <- tolower(as.character(scenarios[[x]]))
    }
    if(x %in% c("num.pops", "Ne", "num.samples", "mig.rate", "mut.rate", 
                "num.loci", "ploidy", "rmetasim.ngen")) {
      scenarios[[x]] <- as.numeric(scenarios[[x]])
    }
  }
  
  params <- list(
    label = label,
    scenarios = scenarios,
    num.rep = num.rep,
    google.drive.id = google.drive.id,
    delete.fsc.files = delete.fsc.files,
    use.wd = use.wd,
    num.cores = num.cores,
    fsc.exec = fsc.exec,
    folders = writeScenarios(label, scenarios, google.drive.id)
  )
  
  params$rep.id <- if(!is.null(params$folders$google)) {
    googledrive::as_id(params$folders$google) 
  } else NULL
  
  params$rep.df <- expand.grid( 
    rep = 1:params$num.rep, 
    sc = 1:nrow(params$scenarios),
    stringsAsFactors = FALSE
  )
  
  params$Rland <- lapply(1:nrow(params$scenarios), .setupScRland, params = params)
  
  sc.rep.vec <- 1:nrow(params$rep.df)
  
  start.time <- Sys.time()
  cat(format(start.time), "Running", nrow(params$scenarios), "scenarios...\n")
  params$scenario.runs <- if(num.cores == 1) {  
    tryCatch(lapply(sc.rep.vec, .runWithLabel, params = params))
  } else {
    cl <- strataG:::.setupClusters(num.cores)
    tryCatch({
      parallel::clusterEvalQ(cl, require(ebvSim))
      parallel::clusterExport(cl, "params", environment())
      parallel::parLapplyLB(cl, sc.rep.vec, .runScRep, params = params)
    }, finally = parallel::stopCluster(cl))
  }
  cat(format(Sys.time()), "Run complete!!\n")
  cat("Total time:", format(round(difftime(Sys.time(), start.time), 2)), "\n")
  
  save(params, file = paste0(params$label, "_params.rdata"))
  invisible(params)
}

#' @noRd
#' 
.runWithLabel <- function(rep.i, params) {
  cat(paste0(
    format(Sys.time()), 
    " ---- Scenario ", params$scenarios$scenario[params$rep.df$sc[rep.i]], 
    ", Replicate ", params$rep.df$rep[rep.i],
    " (", round(100 * rep.i / nrow(params$rep.df)), "%) ----\n"
  ))
  .runScRep(rep.i, params)
}

#' @noRd
#' 
.runScRep <- function(rep.i, params) {
  sc.num <- params$rep.df$sc[rep.i]
  rep.num <- params$rep.df$rep[rep.i]
  
  tryCatch(
    {
      p <- .runFscSim(rep.i, params)
      gen.data <- strataG::fscReadArp(p)
      sc <- params$scenarios[sc.num, ]
      if(sc$rmetasim.ngen > 0) {
        cat(format(Sys.time()), "running rmetasim...\n")
        gen.data <- gen.data %>% 
          calcFreqs() %>% 
          .runRmetasim(Rland = params$Rland[[sc.num]], sc = sc) %>% 
          strataG::landscape2df()
      }
      
      if(!is.na(sc$num.samples)) {
        to.keep <- tapply(1:nrow(gen.data), gen.data$strata, function(i) {
          if(length(i) <= sc$num.samples) i else sample(i, sc$num.samples)
        })
        gen.data <- gen.data[unlist(to.keep), ]
      }
      gen.data$id <- 1:nrow(gen.data)
      
      fname <- repFname(params$label, sc$scenario, rep.num)
      out.name <- file.path(params$folders$out, fname)
      utils::write.csv(gen.data, file = out.name, row.names = FALSE)
      if(!is.null(params$rep.id)) {
        googledrive::drive_upload(
          out.name, 
          path = params$rep.id, 
          name = fname, 
          verbose = FALSE
        )
      }
      
      if(params$delete.fsc.files) strataG::fscCleanup(p$label, p$folder)
      list(fsc.p = p, file = out.name)
    },
    error = function(e) {
      stop(
        format(Sys.time()),
        " Scenario ", params$scenarios$scenario[sc.num],
        ", Replicate ", rep.num, 
        ": ", e
      )
    }
  )
}