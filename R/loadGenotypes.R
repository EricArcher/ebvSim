#' @title Load Genotypes
#' @description Load a data frame of genotypes based on run label, scenario and 
#'   replicate number, and folder.
#' 
#' @param label label of simulation run to download.
#' @param scenario scenario number to retrieve.
#' @param replicate replicate number within the scenario to retrieve.
#' @param folder folder containing the folder of replicates. Default is current 
#'   working directory.
#' 
#' @return data frame of genotypes
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'  
loadGenotypes <- function(label, scenario, replicate, folder = getwd()) {
  fname <- file.path(
    folder,
    repFolder(label), 
    repFname(label, scenario, replicate)
  )
  if(!file.exists(fname)) stop("'", fname, "' cannot be found.")
  utils::read.csv(fname, colClasses = "character", stringsAsFactors = FALSE)
}
