#' @title Download a Simulation Run
#' @description Download scenario specification and genotype csv files from 
#'   Google Drive.
#' 
#' @param label label of simulation run to download.
#' @param google.drive.id id of root Google Drive to download run from.
#' @param out.folder local folder to download run to. If `NULL`, a temporary 
#'   folder is created and used.
#' 
#' @return path of folder where run was downloaded to.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'    
downloadRun <- function(label, google.drive.id, out.folder = NULL) {
  if(is.null(out.folder)) out.folder <- file.path(tempdir(), label)
  if(!dir.exists(out.folder)) dir.create(out.folder)
  contents <- googledrive::drive_ls(googledrive::as_id(google.drive.id))
  i <- grep(label, contents$name)
  
  # download csv
  csv <- grep(".csv$", contents$name[i]) 
  csv.fname <- contents$name[csv]
  cat(csv.fname, "\n")
  googledrive::drive_download(
    googledrive::as_id(contents$id[csv]), 
    path = file.path(out.folder, csv.fname), 
    overwrite = TRUE, 
    verbose = FALSE
  )
  
  # download scenario replicates
  folder <- grepl("_scenario.replicates", contents$name[i]) &
    googledrive::is_folder(contents[i, ])
  folder <- which(folder)
  rep.folder <- file.path(out.folder, contents$name[folder])
  if(!dir.exists(rep.folder)) dir.create(rep.folder)
  fnames <- googledrive::drive_ls(googledrive::as_id(contents$id[folder]))
  for(f in 1:nrow(fnames)) {
    rep.fname <- fnames$name[f]
    cat(f, "/", nrow(fnames), " : ", rep.fname, "\n", sep = "")
    googledrive::drive_download(
      googledrive::as_id(fnames$id[f]), 
      path = file.path(rep.folder, rep.fname), 
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  
  out.folder
}
