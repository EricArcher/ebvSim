#' @title Available Runs
#' @description Return run labels available on a Google Drive.
#' 
#' @param drive.id ID of Google Drive where runs are stored
#' 
#' @return character vector of labels of runs available for download on the
#'   Google Drive requested in 'drive.id'
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'    
availRuns <- function(drive.id) {
  contents <- drive.id %>% 
    googledrive::as_id() %>% 
    googledrive::drive_ls()
  
  folder.ext <- "_scenario.replicates"
  i <- grepl(folder.ext, contents$name) & googledrive::is_folder(contents) %>% 
    which()
  
  gsub(folder.ext, "", contents$name[i])
}
