library(detectRUNS)
##use 'strataG' and 'hierfstat' packages to create Plink-like .ped .map files
##for calculating Froh with 'detectRUNS' package
snp.hf <- genind2hierfstat(gtypes2genind(sim.snps.g, type = "codom"))
ped.label <- "testhierfped"
write.ped(snp.hierfstat, ilab = getIndNames(sim.snps.g), fname = ped.label)

#the following lines must be run to change the separator of the .ped file
#generated with 'write.ped' to a space. 'detectRUNS' will not work without this
#step. apparently this is not necessary for the .map file
pedgood <- read.table(paste0(ped.label, ".ped"), header = FALSE)
froh.ped <- "testpedFROH.ped"
write.table(pedgood, froh.ped, sep=" ", col.names = FALSE, row.names = FALSE)

#detectRUNS needs base-pair coordinates to run (column 4 in a .map file, equal
#to 0 for all our simulated SNPs at this point) here, I arbitrarily gave a
#base-pair coordinate to each variant ranging from 1 to the maximum number of
#variants
mapfile <- read.table(paste0(ped.label, ".map"), header = FALSE)
mapfile[1:nrow(mapfile), 4] <- 1:nrow(mapfile)
froh.map <- "testmapFROH.map"
write.table(mapfile, froh.map, sep=" ", col.names = FALSE, row.names = FALSE)

#detect runs using the consecutive method (Marras et al 2015 An Gen). An
#alternative method can be also calculated with 'detectRUNS' I didn't really
#checked the validity of the maxOppRun, maxMissRun etc... parameters. I just
#focused on 'runability'
runs <- consecutiveRUNS.run(
  froh.ped, 
  froh.mat,
  ROHet = FALSE, 
  maxOppRun = 0,
  maxMissRun = 0, 
  minSNP = 15, 
  minLengthBps = 10, 
  maxGap = 10^6
)

#calculate Froh for all individuals
froh.res <- Froh_inbreeding(runs, froh.map, genome_wide = TRUE)

#create a matrix for extracting mean and sd Froh estimates for each population
froh.smry <- tapply(res[, 4], res$group, function(x) {
  c(mean = mean(x), sd = sd(x))
})