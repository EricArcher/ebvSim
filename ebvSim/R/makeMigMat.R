#' @title Make Migration Matrix
#' @description Create island or stepping stone migration matrix given 
#'   migration rate and number of populations.
#' 
#' @param mig.rate migration rate specified as proportion of population 
#'   migrating per generaton.
#' @param num.pops number of populations.
#' @param type of migration matrix structure. "island" = rate between all 
#'   populations is the same. "stepping.stone" = migration only occurs between 
#'   neighboring populations. Metapopulation is ring shaped, not a linear chain.
#' 
#' @return matrix of pairwise migration rates.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'  
makeMigMat <- function(mig.rate, num.pops, 
                       type = c("island", "stepping.stone")) {
  type <- match.arg(type)
  mig.mat <- switch(
    type,      
    island = {
      m <- mig.rate / (num.pops - 1)
      matrix(rep(m, num.pops ^ 2), nrow = num.pops)
    },
    stepping.stone = {
      mat <- matrix(0, nrow = num.pops, ncol = num.pops)
      m <- mig.rate / 2
      for (k in 1:(num.pops - 1)) {
        mat[k, k + 1] <- mat[k + 1, k] <- m
      }
      mat[1, num.pops] <- mat[num.pops, 1] <- m
      mat
    }
  )
  diag(mig.mat) <- 1 - mig.rate
  mig.mat
}
