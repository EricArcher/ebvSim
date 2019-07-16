#' @title Calculate Allele Frequencies
#' @description Calculate overall and by-strata allele frequencies for input to 
#'   rmetasim.
#' 
#' @param snps data frame of SNP genotypes with first column called 'id', 
#'   second column called 'stratum' and every two columns thereafter being 
#'   SNP loci.
#' 
#' @return list of overall ('global') and by-strata ('pop') allele frequencies
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'    
calcFreqs <- function(snps) {
  mac.df <- snps %>% 
    strataG::df2gtypes(ploidy = 2) %>% 
    strataG::as.data.frame(coded = T) %>% 
    dplyr::select(-.data$id) %>% 
    tidyr::gather("locus", "mac", -.data$stratum) %>% 
    dplyr::mutate(
      stratum = as.numeric(factor(.data$stratum)),
      locus = as.numeric(factor(.data$locus))
    )
  
  list(
    global = lapply(split(mac.df, mac.df$locus), function(loc.df) {
      alleleProp(loc.df$mac)
    }),
    pop = lapply(split(mac.df, mac.df$stratum), function(st.df) {
      lapply(split(st.df, st.df$locus), function(loc.df) alleleProp(loc.df$mac))
    })
  )
}


#' @rdname calcFreqs
#' 
#' @param mac a vector of numerics giving the count of the minor allele, where 
#'   0 = homozygous for reference allele, 1 = heterozygous, 2 = homozygous for
#'   minor allele
#' 
#' @export
#'    
alleleProp <- function(mac) {
  maf <- mean(mac == 0) + (mean(mac == 1) / 2)
  c('1' = maf, '2' = 1 - maf)
}