######## Works on R 3.5.3 , Ubuntu v 18.10 running on a Virtual Machine (VMware workstation player 15)
######## 6 cores (2.8Ghz), 22 Gb Ram

library(strataG) #genetic data simulation (SNP, SSR, DNA sequence) with fastsimcoal, LDNe calculation, intermediate I/O file conversion step for calculating Froh with detectRuns
library(radiator) #intermediate I/O file conversion step for calculating g2 with inbreedR
library(vcfR) #intermediate I/O file conversion step for calculating g2 with inbreedR
library(reshape) #intermediate I/O file conversion step for calculating g2 with inbreedR
library(inbreedR) #g2 statistic calculation
library(hierfstat) #intermediate I/O file conversion step for calculating Froh with detectRuns
library(detectRUNS) #Froh statistic calculation

start_time=Sys.time()

# Set fastsimcoal parameters # Population information: 3 populations with Ne = 500, drawing 100 samples from each.
pop.info <- fscPopInfo(pop.size = rep(500, 3), sample.size = rep(100, 3))

# Migration rates among the 3 populations
mig.rates <- matrix(c(0, 0.5, 0.005, 0.5, 0, 0.0005, 0.005, 0.0005, 0), ncol = 3)

# Define historical events in which populations diverged 2000 generations in past
hist.ev <- fscHistEv(num.gen = c(2000, 2000), source.deme = c(2, 1),sink.deme = c(1, 0), prop.migrants = 1)

# snp, dna and ssr params
snp.params <- fscLocusParams(locus.type = "snp", num.loci = 1000, mut.rate = 0.001, recomb.rate=0,ploidy=2)
dna.params <- fscLocusParams(locus.type = "dna", sequence.length = 1000, mut.rate = 0.0000001, recomb.rate=0, transition.rate=0.66, ploidy=2)
ssr.params <- fscLocusParams(locus.type = "msat", num.loci = 15, mut.rate = 0.0005, recomb.rate=0,ploidy=2)

# Run simulations using fastsimcoal 2.6 (fsc26 exec file must be in the working directory)
sim.snp <- fastsimcoal(pop.info, snp.params, mig.rates, hist.ev, exec=paste(getwd(),"/fsc26",sep=""))
sim.dna <- fastsimcoal(pop.info, dna.params, mig.rates, hist.ev, exec=paste(getwd(),"/fsc26",sep=""))
sim.ssr <- fastsimcoal(pop.info, ssr.params, mig.rates, hist.ev, exec=paste(getwd(),"/fsc26",sep=""))

# Create a reference table to store model parameter values and calculated summary statistics
refTable=matrix(nrow=1,ncol=8)
colnames(refTable)=c("pop_size_sample_1","pop_size_sample_2","pop_size_sample_3","sampl_size_sample_1","sampl_size_sample_2","sampl_size_sample_3","t_histevent_1","t_histevent_2")

# fill the table with model parameter values
refTable[1,1:3]=pop.info[1:3,1]
refTable[1,4:6]=pop.info[1:3,2]
refTable[1,7:8]=hist.ev[1:2,1]

#Calculate Linkage Disequilibrium calculated for the SNP dataset (strataG package)
ldsnp=ldNe(sim.snp)

#Create a matrix for extracting LDNe estimates, lower (LCI) and upper (UCI) confidence intervals for each population
ldne=matrix(nrow=1,ncol=9)
colnames(ldne)=c("LDNE_sample_1","LDNE_LCI_sample_1","LDNE_UCI_sample_1","LDNE_sample_2","LDNE_LCI_sample_2","LDNE_UCI_sample_2","LDNE_sample_3","LDNE_LCI_sample_3","LDNE_UCI_sample_3")
ldne[1,]=c(ldsnp[1,5:7],ldsnp[2,5:7],ldsnp[3,5:7])

#Add LDNe results to the Reference Table
refTable=cbind(refTable,ldne)

##Use package 'radiator' to convert the snp dataset from gtypes format to vcf format
genomic_converter(sim.snp,output="vcf",filename="testvcf")

#Copy the vcf format of interest to the working directory, then remove the directory created by 'radiator'
#this may be helpful when performing many simulations, since radiator will create a new directory for each simulation

rad_cmd1=paste("cp ",getwd(), "/*radiator*/testvcf.vcf"," ", getwd(),sep="")
rad_cmd2=paste("rm -r *radiator*")
system(rad_cmd1,intern=F)
system(rad_cmd2,intern=F)

#need packages 'vcfR', 'reshape' and 'inbreedR' for converting the raw snp data on vcf format to the input format of 'inbreedR'

#read vcf using 'vcfR'
vcf <- read.vcfR("testvcf.vcf", verbose = FALSE )

#extract genotypes
gt <- extract.gt(vcf,IDtoRowNames=F)

#transpose and data.frame
gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)

#NA handling (not necessary for simulated data)
#gt[gt == "."] <- NA

#split columns using do.call from 'reshape' package
snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))

#convert to inbreed input. Here, I created three subdatasets (one per pop) + a global dataset to calculate g2 at both local & global scales
testinbreedR_pop1 <- inbreedR::convert_raw(snp_geno[1:100,])
testinbreedR_pop2 <- inbreedR::convert_raw(snp_geno[101:200,])
testinbreedR_pop3 <- inbreedR::convert_raw(snp_geno[201:300,])
testinbreedR_global <- inbreedR::convert_raw(snp_geno)

#calculating g2 for each population (+ global) using snp data. Can be also calculated with SSRs if necessary (I didn't tried yet)
#here, 10 bootstraps are performed for CI calculation. Also, a p-value can be computed using permutations
#both parameters will impact computation time
g2pop1=g2_snps(testinbreedR_pop1,nboot=10)
g2pop2=g2_snps(testinbreedR_pop2,nboot=10)
g2pop3=g2_snps(testinbreedR_pop3,nboot=10)
g2global=g2_snps(testinbreedR_global,nboot=10)

#Create a matrix for extracting g2 estimates, lower (LCI) and upper (UCI) confidence intervals for each population and for the global dataset
g2res=matrix(nrow=1,ncol=12)
colnames(g2res)=c("g2_sample_1","g2_LCI_sample_1","g2_UCI_sample_1","g2_sample_2","g2_LCI_sample_2","g2_UCI_sample_2","g2_sample_3","g2_LCI_sample_3","g2_UCI_sample_3","g2_global","g2_LCI_global","g2_UCI_global")
g2res[,1:3]=c(g2pop1$g2,g2pop1$CI_boot[1],g2pop1$CI_boot[2])
g2res[,4:6]=c(g2pop2$g2,g2pop2$CI_boot[1],g2pop2$CI_boot[2])
g2res[,7:9]=c(g2pop3$g2,g2pop3$CI_boot[1],g2pop3$CI_boot[2])
g2res[,10:12]=c(g2global$g2,g2global$CI_boot[1],g2global$CI_boot[2])

#Add g2 stats to the reference table
refTable=cbind(refTable,g2res)

##use 'strataG' and 'hierfstat' packages to create Plink-like .ped .map files for calculating Froh with 'detectRUNS' package

snp.genind=gtypes2genind(sim.snp, type = "codom")
snp.hierfstat=genind2hierfstat(snp.genind,pop=c(rep(1,100),rep(2,100),rep(3,100)))
write.ped(snp.hierfstat, ilab=seq(1,300,1),fname="testhierfped")

#the following lines must be run to change the separator of the .ped file generated with 'write.ped' to a space. 'detectRUNS' will not work without this step.
#apparently this is not necessary for the .map file
pedgood=read.table("testhierfped.ped",h=F)
write.table(pedgood,"testpedFROH.ped",sep=" ",col.names=F,row.names=F)

#detectRUNS needs base-pair coordinates to run (column 4 in a .map file, equal to 0 for all our simulated SNPs at this point)
#here, I arbitrarily gave a base-pair coordinate to each variant ranging from 1 to the maximum number of variants
mapfile=read.table("testhierfped.map",h=F)
mapfile[1:nrow(mapfile),4]=seq(1,nrow(mapfile),1)
write.table(mapfile,"testmapFROH.map",sep=" ",col.names=F,row.names=F)

#detect runs using the consecutive method (Marras et al 2015 An Gen). An alternative method can be also calculated with 'detectRUNS'
#I didn't really checked the validity of the maxOppRun, maxMissRun etc... parameters. I just focused on 'runability' 
FROH_ped_input=paste(getwd(),"/testpedFROH.ped",sep="")
FROH_map_input=paste(getwd(),"/testmapFROH.map",sep="")
runs=consecutiveRUNS.run(FROH_ped_input, FROH_map_input, ROHet = FALSE, maxOppRun = 0,maxMissRun = 0, minSNP = 15, minLengthBps = 10, maxGap = 10^6)

#calculate Froh for all individuals
res=Froh_inbreeding(runs, FROH_map_input, genome_wide = TRUE)

#create a matrix for extracting mean and sd Froh estimates for each population
matFROH=matrix(nrow=1,ncol=6)
colnames(matFROH)=c("mean_FROH_sample_1","sd_FROH_sample_1","mean_FROH_sample_2","sd_FROH_sample_2","mean_FROH_sample_3","sd_FROH_sample_3")

matFROH[,1]=mean(res[which(res$group==1),4])
matFROH[,2]=sd(res[which(res$group==1),4])
matFROH[,3]=mean(res[which(res$group==2),4])
matFROH[,4]=sd(res[which(res$group==2),4])
matFROH[,5]=mean(res[which(res$group==3),4])
matFROH[,6]=sd(res[which(res$group==3),4])

#add Froh stats to the reference table
refTable=cbind(refTable,matFROH)

#remove all intermediate files from the working directory
rm_interm_files=paste("rm test*")
system(rm_interm_files,intern=F)

#remove all R objects excepting the refTable
rm(list=setdiff(ls(),c("refTable","start_time")))

write.table(refTable,"GEOBON_TESTMODEL_REFTABLE.txt", dec=".",sep="\t")

end_time=Sys.time()
runtime=end_time-start_time
gc()


