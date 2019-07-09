library(inbreedR)

cat(format(Sys.time()), "Computing g2...\n")

# create list of genotypes in inbreedR input format
# one element for each population and a final element of all genotypes
pop.geno <- split(sim.snps, sim.snps$deme)
geno.inbrd <- sapply(
  pop.geno, function(x) convert_raw(x[, -(1:2)]), simplify = FALSE
)
geno.inbrd$global <- convert_raw(sim.snps[, -(1:2)])

#calculating g2 for each population (+ global) using snp data. Can be also
#calculated with SSRs if necessary (I didn't tried yet) here, 10 bootstraps are
#performed for CI calculation. Also, a p-value can be computed using
#permutations both parameters will impact computation time
#Create a matrix for extracting g2 estimates, lower (LCI) and upper (UCI)
#confidence intervals for each population and for the global dataset
g2 <- t(sapply(geno.inbrd, function(x) {
  res <- g2_snps(x, nboot = 10)
  c(g2 = res$g2, LCI = res$CI_boot[1], UCI = res$CI_boot[2])
})) %>% 
  as.data.frame() %>% 
  rownames_to_column("stratum")
