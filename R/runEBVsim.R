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
#' @param num.cores number of cores to use. Each scenario is allocated to a 
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
  
  num.sc <- nrow(params$scenarios)
  params$scenario.runs <- if(num.cores == 1) {  
    lapply(1:num.sc, function(sc.i) {  
      cat(paste0(
        format(Sys.time()), 
        " ---- Scenario ", 
        params$scenarios$scenario[sc.i], 
        " (", sc.i, " / ", num.sc , 
        ") ----\n"
      ))
      .runScenario(sc.i, params = params, quiet = FALSE)
    })
  } else {
    cl <- strataG:::.setupClusters(num.cores)
    tryCatch({
      parallel::clusterEvalQ(cl, require(ebvSim))
      parallel::clusterExport(cl, "params", environment())
      parallel::parLapply(cl, 1:num.sc, .runScenario, params = params)
    }, finally = parallel::stopCluster(cl))
  }
  
  save(params, file = paste0(params$label, "_params.rdata"))
  invisible(params)
}

#' @noRd
#' 
.runScenario <- function(sc.i, params, quiet = TRUE) {
  p <- .runFscSim(params, sc.i)
  
  files <- sapply(1:params$num.rep, function(sim.i) {
    if(!quiet) {
      cat(paste0(
        format(Sys.time()),
        " *** Replicate ", sim.i, " / ", params$num.rep, "***\n"
      ))
    }
    gen.data <- strataG::fscReadArp(p, sim = c(1, sim.i))
    if(params$scenarios$rmetasim.ngen[sc.i] > 0) {
      gen.data <- gen.data %>% 
        calcFreqs() %>% 
        .runRmetasim(params$scenarios[sc.i, ]) %>% 
        strataG::landscape2gtypes() %>% 
        strataG::as.data.frame()
    }
    
    fname <- repFname(params$label, params$scenarios$scenario[sc.i], sim.i)
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
    out.name
  })
  
  if(params$delete.fsc.files) strataG::fscCleanup(p$label, p$folder)
  list(fsc.p = p, files = files)
}