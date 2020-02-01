# ebvSim

### Description

*ebvSim* is a package for simulating SNP data and testing the performance of various genetic diversity metrics. The package will simulate replicate data for multiple scenarios. It is a wrapper that the coalescent simulator *fastsimcoal2* through the *strataG* package, and then runs a few generations of a forward simulator, *rmetasim*, initialized with allele frequencies from the genotypes output by *fastsimcoal2*.

***

### Installation

1\. Download and install  [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2) so that it can be executed on the command line from any location. This requires it to be in a folder that is somewhere in the execution PATH. The test is to try to execute fastsimcoal (`fsc26`) from any folder. Here are sites with guidance for setting the path for different operating systems:
[Windows](https://www.java.com/en/download/help/path.xml), 
[Mac OSX](http://osxdaily.com/2014/08/14/add-new-path-to-path-command-line/), or 
[LINUX](http://www.wikihow.com/Change-the-Path-Variable-in-Linux). On Mac OS or LINUX/UNIX systems, executables can be placed in the _/usr/local/bin_ folder which is usually a default in the PATH.

2\. Install the `ebvSim` package from GitHub with: 

```r
devtools::install_github("ericarcher/ebvSim", dependencies = TRUE)
```
This *should* also install *strataG* and *rmetasim* from their GitHub repositories. If these packages are not available after the *ebvSim* installation, install them from:

```r
devtools::install_github("ericarcher/strataG", ref = "redevel2018", dependencies = TRUE)
devtools::install_github("stranda/rmetasim", dependencies = TRUE)
```

3\. By default, *rmetasim* can only simulate a maximum of 1001 loci. If this needs to be increased, it can be done so by changing this constant and recompiling. Instructions for this can be found [here](https://thierrygosselin.github.io/grur/articles/rad_genomics_computer_setup.html#rmetasim).

***

### Create simulation scenarios

A data frame of scenarios must first be created. The data frame should have the following columns:

* __num.pops__: the number of populations.
* __Ne__: the effective population size.
* __num.samples__:	the number of samples to simulate. Must be <= Ne.
* __mig.rate__: the migration rate specified as proportion of population migrating per generaton.
* __mig.type__: the type of of migration matrix structure. "island" = rate between all populations is the same. "stepping.stone" = migration only occurs between neighboring populations. Metapopulation is ring shaped, not a linear chain.
* __dvgnc.time__: the number of generations since divergence of the populations.
* __rmetasim.ngen__: the number of generations to run Rmetasim for. Set to 0 to skip Rmetasim.

This data frame can be created in an R script or read from an external file. If you want to create a data frame of all combinatons of vectors of parameters, you can use the function `makeScenarios()`:

```r
library(ebvSim)

# creates six scenarios of all listed combinations of number of populations and effective populations size
scenarios <- makeScenarios(
  num.pops = c(1, 3, 5),
  Ne = c(10, 100),
  num.samples = 10,
  mig.rate = 1e-5,
  mig.type = "island",
  dvgnc.time = 100, 
  rmetasim.ngen = 10
)
```

***

### Run the simulations

With the scenarios specified, the simulation can be run with the `runEBVsim()` function:

```r
library(ebvSim)

# a single scenario
scenarios <- data.frame(
  num.pops = 2,
  Ne = 100,
  num.samples = 10,
  mig.rate = 1e-5,
  mig.type = "island",
  dvgnc.time = 100, 
  rmetasim.ngen = 10
)

genetics <- strataG::fscSettingsGenetics(
  strataG::fscBlock_snp(
    sequence.length = 1, mut.rate = 1e-4
  ), 
  num.chrom = 1000
)

params <- runEBVsim(
  label = "ebvSim.snps_test",
  scenarios = scenarios,
  genetics = genetics,
  ploidy = 2,
  num.rep = 10,
  num.cores = 4
)
```

The parameters for `runEBVsim()` are:

* __label__: text to use to label this set of of scenarios.
* __scenarios__: the scenarios data frame previously created.
* __genetics__: specification for the genetic marker to be simulated. This should be the output of the `fscSettingsGenetics()` function from the _strataG_ package. This function takes the output of a function like `fscBlock_snp()` to simulate SNP loci.
* __ploidy__: the ploidy of the markers to be simulated. Set to 2 for diploid.
* __num.rep__: the number of replicates to simulate for each scenario.
* __num.cores__: the number of cores to use. replicates for each scenario are assigned to one core. If this is set to a value greater than 1, progress notifications will not be printed on the console.

When the simulations are complete, there will be three new items in the working directory:

* \<label\>_scenario.replicates: a folder containig .csv files of genotypes for each scenario replicate.
* \<label\>_scenarios.csv: a .csv file of the scenario specifications.
* \<label\>_params.rdata: an R workspace file containing a a list called `params` that contains the parameters used to to run the scenarios and an element called `$scenario.runs` which is a list of the fastsimcoal output (`$fsc.p`) and genotype files output for each replicate (`$files`). This file will be written when the simulation is complete.

The `runEBVsim()` function also invisibly returns the same summary list contained in \<label\>_params.rdata. 

***

### Analyze the replicates

Sets of scenarios can be analyzed with the `analyzeReps()` function. The `params` list saved or returned by `runEBVsim()` needs to be provided which can be found in the saved workspace file. Note that the folder containing the scenario replicate genotype .csv files must also be in the current working directory.

Each analysis is run separately, by providing the name of the desired analysis as the first argument to `analyzeReps()`. Current options are `ldNe` (linkage disequilibrium effective population size), `g2` (inbreeding g2), `het` (observed heterozygosity), `froh` (runs of homozygosity).

The parameter `num.cores` can be specified to spread analyses of replicates among multiple cores.

```r
# load params into the workspace
load("ebvSim.snps_test_params.rdata")

# calculate effective population size
ne <- analyzeReps("ldNe", params)

# calculate inbreeding, heterozygosity, and runs of homozygosity using four cores
g2 <- analyzeReps("g2", params, num.cores = 4)
het <- analyzeReps("het", params, 4)
froh <- analyzeReps("froh", params, 4)
```

All results from each `analyzeReps()` analysis can then be combined by joining the data frames by the __scenario__, __replicate__, and __stratum__ columns.

```r
library(tidyverse)
join.by <- c("scenario", "replicate", "stratum")
all.ebvs <- ne %>% 
  left_join(g2, by = join.by) %>% 
  left_join(het, by = join.by) %>% 
  left_join(froh, by = join.by)
  
# save results to workspace file
save.image(paste0(params$label, "_analysis results.rdata"))
```

***

### Contact

* submit suggestions and bug-reports: <https://github.com/ericarcher/ebvSim/issues>
* send a pull request: <https://github.com/ericarcher/ebvSim/>
* e-mail: <eric.archer@noaa.gov>