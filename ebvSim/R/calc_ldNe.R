#' @title Calculate ldNe
#' @description Calculate ldNe for data frame of genotypes
#' 
#' @param x data frame of genotypes as read from simulation replicate .csv file
#' 
#' @return estimated Ne for each stratum in file with lower and upper CIs
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'
calc_ldNe <- function(x) {
  x.ne <- x %>% 
    strataG::df2gtypes(ploidy = 2) %>% 
    strataG::ldNe()
  
  if(!is.null(x.ne)) {
    x.ne %>% 
      dplyr::select(.data$stratum, .data$Ne, .data$param.lci, .data$param.uci) %>% 
      stats::setNames(c("stratum", "Ne", "Ne.lci", "Ne.uci"))
  } else NULL
}