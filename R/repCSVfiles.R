#' @title Get replicate .csv file names
#' @description Return data frame of replicate file names, and extract
#' scenario and replicate numbers
#'   
#' @param folder path of folder to search in.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'  
repCSVfiles <- function(folder) {
  fnames <- dir(
    folder, 
    pattern = "_scenario[///.][[:digit:]]+_replicate[///.][[:digit:]]+[///.]csv$", 
    full.names = TRUE, 
    recursive = TRUE
  )
  
  data.frame(fname = fnames, stringsAsFactors = FALSE) %>% 
    dplyr::mutate(
      scenario = stringr::str_extract(.data$fname, "_scenario[///.][[:digit:]]+_replicate"),
      scenario = as.numeric(stringr::str_extract(.data$scenario, "[[:digit:]]+")),
      replicate = stringr::str_extract(.data$fname, "_replicate[///.][[:digit:]]+[///.]csv$"),
      replicate = as.numeric(stringr::str_extract(.data$replicate, "[[:digit:]]+"))
    ) %>% 
    dplyr::arrange(.data$scenario, .data$replicate)}