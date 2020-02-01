#' @title Calculate Runs of Homozygosity
#' @description Calculate Runs of Homozygosity
#' 
#' @param x data frame of genotypes as read from simulation replicate .csv file.
#' 
#' @return data frame of ROH summary 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
calc_froh <- function(x) {
  # Create 'hierfstat' input file by converting MAC genotypes to integer alleles
  snp.hf <- x %>% 
    strataG::df2gtypes(ploidy = 2) %>% 
    strataG::as.data.frame(coded = T)
  for(i in 3:ncol(snp.hf)) snp.hf[, i] <- c(11, 12, 22)[snp.hf[, i] + 1]
  ids <- snp.hf$id
  snp.hf$id <- NULL
  
  # Write .ped and .map files
  
  #ped.label <- file.path(tempdir(), "testhierf")
  ped.label <- tempfile()
  hierfstat::write.ped(snp.hf, ilab = ids, fname = ped.label)
  
  # Change separator of .ped file to space for detectRUNS
  froh.ped <- paste0(ped.label, "FROH.ped")
  utils::write.table(
    utils::read.table(paste0(ped.label, ".ped"), header = FALSE),
    froh.ped, sep = " ", col.names = FALSE, row.names = FALSE
  )
  
  # detectRUNS needs base-pair coordinates to run (column 4 in a .map file, equal
  # to 0 for all our simulated SNPs at this point) here, I arbitrarily gave a
  # base-pair coordinate to each variant ranging from 1 to the maximum number of
  # variants
  mapfile <- utils::read.table(paste0(ped.label, ".map"), header = FALSE)
  mapfile[1:nrow(mapfile), 4] <- 1:nrow(mapfile)
  froh.map <- paste0(ped.label, "FROH.map")
  utils::write.table(mapfile, froh.map, sep=" ", col.names = FALSE, row.names = FALSE)
  
  # Detect runs using the consecutive method (Marras et al 2015 An Gen) and 
  # calculate Froh for each individual
  froh.res <- detectRUNS::consecutiveRUNS.run(
    froh.ped, froh.map, ROHet = FALSE, maxOppRun = 0, maxMissRun = 0, 
    minSNP = 15, minLengthBps = 10, maxGap = 10^6
  ) %>% 
    dplyr::mutate(group = gsub("\"", "", .data$group), id <- gsub("\"", "", id)) %>% 
    detectRUNS::Froh_inbreeding(froh.map, genome_wide = TRUE)
  
  # Delete files created
  file.pattern <- paste0("^", ped.label, "[[:alnum:]]*(.map$)|(.ped$)")
  file.remove(dir(pattern = file.pattern))
  
  # Calculte mean, standard deviation, median, and 95% CI of Froh for each 
  # population
  tapply(froh.res[, 4], froh.res$group, function(x) {
    data.frame(
      froh.mean = mean(x), 
      froh.sd = stats::sd(x), 
      froh.median = stats::median(x), 
      froh.lci = unname(stats::quantile(x, 0.025)), 
      froh.uci = unname(stats::quantile(x, 0.975))
    )
  }) %>% 
    dplyr::bind_rows() %>% 
    tibble::rownames_to_column("stratum")
}