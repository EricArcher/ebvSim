#' @title Load Scenarios
#' @description Load a data frame of scenario specifications for a simulation 
#'   run.
#' 
#' @param label label of simulation run.
#' @param folder folder containing the scenario data file. Default is current 
#'   working directory.
#' 
#' @return data frame of scenario specifications.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'  
loadScenarios <- function(label, folder = getwd()) {
  fname <- scenariosFname(label)
  utils::read.csv(file.path(folder, fname), stringsAsFactors = FALSE)
}
