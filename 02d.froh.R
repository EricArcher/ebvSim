library(detectRUNS)
library(hierfstat)

cat(format(Sys.time()), "Computing Froh...\n")

#' Create 'hierfstat' input file by converting 0:2 genotypes to integer alleles
snp.hf <- fscReadArp(snp.p, coded.snps = TRUE)
for(i in 3:ncol(snp.hf)) snp.hf[, i] <- c(11, 12, 22)[snp.hf[, i] + 1]
ids <- snp.hf$id
snp.hf$id <- NULL

#' Write .ped and .map files
ped.label <- "testhierf"
write.ped(snp.hf, ilab = ids, fname = ped.label)

#' Change separator of .ped file to space for detectRUNS
froh.ped <- paste0(ped.label, "FROH.ped")
write.table(
  read.table(paste0(ped.label, ".ped"), header = FALSE),
  froh.ped, sep = " ", col.names = FALSE, row.names = FALSE
)

#' detectRUNS needs base-pair coordinates to run (column 4 in a .map file, equal
#' to 0 for all our simulated SNPs at this point) here, I arbitrarily gave a
#' base-pair coordinate to each variant ranging from 1 to the maximum number of
#' variants
mapfile <- read.table(paste0(ped.label, ".map"), header = FALSE)
mapfile[1:nrow(mapfile), 4] <- 1:nrow(mapfile)
froh.map <- paste0(ped.label, "FROH.map")
write.table(mapfile, froh.map, sep=" ", col.names = FALSE, row.names = FALSE)

#' Detect runs using the consecutive method (Marras et al 2015 An Gen) and 
#' calculate Froh for each individual
froh.res <- consecutiveRUNS.run(
  froh.ped, froh.map, ROHet = FALSE, maxOppRun = 0, maxMissRun = 0, 
  minSNP = 15, minLengthBps = 10, maxGap = 10^6
) %>% 
  mutate(group = gsub("\"", "", group), id <- gsub("\"", "", id)) %>% 
  Froh_inbreeding(froh.map, genome_wide = TRUE)

#' Calculte mean, standard deviation, median, and 95% CI of Froh for each 
#' population
froh.smry <- do.call(rbind, tapply(froh.res[, 4], froh.res$group, function(x) {
  c(
    mean = mean(x), sd = sd(x), median = median(x), 
    LCI = unname(quantile(x, 0.025)), UCI = unname(quantile(x, 0.975))
  )
})) %>% 
  as.data.frame() %>% 
  rownames_to_column("stratum")

#' Delete files created
file.pattern <- paste0("^", ped.label, "[[:alnum:]]*(.map$)|(.ped$)")
file.remove(dir(pattern = file.pattern))