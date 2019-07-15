#' @title Replicate File Name
#' @description Return file name containing genotypes from a simulation 
#'   replicate.
#'   
#' @param label label of simulation run.
#' @param scenario scenario number.
#' @param replicate replicate number within the scenario.
#' 
#' @return for `repFname`, the name of file that should contain genotypes,
#'   for `repFolder`, the name of the folder that should contain simulation 
#'   replicate files.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'  
repFname <- function(label, scenario, replicate) {
  paste0(label, "_scenario.", scenario, "_replicate.", replicate, ".csv")
}

#' @rdname repFname
#' @export
#' 
repFolder <- function(label) paste0(label, "_scenario.replicates")

#' @rdname repFname
#' @export
#' 
scenariosFname <- function(label) paste0(label, "_scenarios.csv")