#' ebvSim
#' 
#' @docType package
#' @name ebvSim-package
#' @aliases ebvSim
#' @title Essential Biodiversity Variables Genetic Simulations
#' @keywords package
#' @importFrom rlang .data
#' @importFrom magrittr %>%
NULL

#' @docType data
#' @name trial.1
#' @title Baseline Scenarios - Trial 1
#' @description A data.frame of 96 scenarios for first trial of baseline 
#'   EBV simulations. All scenarios generate 1000 unlinked SNP loci that 
#'   have gone through 10 generations of random mating in the rmetasim 
#'   forward simulator. The mutation rate is set at 1e-4 substitutions 
#'   per locus per generation. Multiple populations are connected in 
#'   an island arrangement and diverged 1000 generations in the past. 
#'   Indivdiual scenarios are then all permutations of:
#'   \tabular{ll}{
#'     \code{number of populations} \tab 1, 2, 10, 50\cr
#'     \code{effective population size} \tab 20, 50, 100, 250, 500, 1000\cr
#'     \code{migration rate} \tab 1e-4, 1e-3, 1e-2, 0.1, 0.25\cr
#'   }
#' @usage data(trial.1)
#' @format data.frame
#' @keywords datasets
NULL