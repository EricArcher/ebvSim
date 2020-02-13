# ebvSim

### Description

*ebvSim* is a package for simulating SNP data and testing the performance of various genetic diversity metrics. The package will simulate replicate data for multiple scenarios. It is a wrapper for the coalescent simulator *fastsimcoal2* run through the *strataG* package, which then runs a few generations of a forward simulator, *rmetasim*, initialized with allele frequencies from the genotypes output by *fastsimcoal2*.

***

### Installation

1\. Download and install  [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2) so that it can be executed on the command line from any location. This requires it to be in a folder that is somewhere in the execution PATH. The test is to try to execute fastsimcoal (`fsc26`) from any folder. Here are sites with guidance for setting the path for different operating systems:
[Windows](https://www.java.com/en/download/help/path.xml), 
[Mac OSX](http://osxdaily.com/2014/08/14/add-new-path-to-path-command-line/), or 
[LINUX](http://www.wikihow.com/Change-the-Path-Variable-in-Linux). On Mac OS or LINUX/UNIX systems, executables can be placed in the _/usr/local/bin_ folder which is usually a default in the PATH.

2\. Install the `ebvSim` package from GitHub with: 

```r
devtools::install_github("ericarcher/ebvSim", dependencies = TRUE, force = TRUE)
```
This *should* also install *strataG* and *rmetasim* from their GitHub repositories. If these packages are not available after the *ebvSim* installation, install them from:

```r
devtools::install_github("ericarcher/strataG", dependencies = TRUE, force = TRUE)
devtools::install_github("stranda/rmetasim", dependencies = TRUE, force = TRUE)
```

3\. By default, *rmetasim* can only simulate a maximum of 1001 loci. If this needs to be increased, it can be done so by changing this constant and recompiling. Instructions for this can be found [here](https://thierrygosselin.github.io/grur/articles/rad_genomics_computer_setup.html#rmetasim).

***

### Select simulation scenarios
The code is designed to run multiple replicates of a set of demographic scenarios. The scenarios are defined by rows in data frames that are included in the package. The first set of scenarios is called `trial.1`, and has the following columns: 

* __num.pops__: the number of populations.
* __Ne__: the effective population size.
* __num.samples__:	the number of samples to simulate. Must be <= Ne. If NA, then num.samples = Ne.
* __mig.rate__: the migration rate specified as proportion of population migrating per generaton.
* __mig.type__: the type of of migration matrix structure. "island" = rate between all populations is the same. "stepping.stone" = migration only occurs between neighboring populations. Metapopulation is ring shaped, not a linear chain.
* __dvgnc.time__: the number of generations since divergence of the populations.
* __marker.type__: type of marker to simulate. At this point, only "snp" is available.
* __mut.rate__: mutation rate (# of mutations per generation) of the markers to simulate.
* __num.loci__: number of independent loci to simulate.
* __ploidy__: the ploidy of the markers to be simulated. Set to 2 for diploid.
* __rmetasim.ngen__: the number of generations to run Rmetasim for. Set to 0 to skip Rmetasim.

For the first set of trials, we would like people to sign up to run __100 replicates__ of each scenario on the google spreadsheet located [here](https://docs.google.com/spreadsheets/d/1-o7dPFz9l8w2Eh0sox7CV-s6d3NP1XWy5-vCp7kaFEA/edit?usp=sharing).

***

### Run the simulations

To run a specific set of scenarios, the simulation can be run with the `runEBVsim()` function:

```r
rm(list = ls())
library(ebvSim)
data(trial.1)

# create a vector of specific scenarios
i <- c(1, 5, 10, 12, 20)

runEBVsim(
  label = "EIA_trial.1_1",
  scenarios = trial.1[i, ],
  num.rep = 10,
  num.cores = 4
)
```

The parameters for `runEBVsim()` are:

* __label__: text to use to label this set of of scenarios. The format should be "initials_trial.name_attempt.num". attempt.num can be a number you use to separate different attempts if you're doing either multiple trials or multiple runs of the same set of scenarios on different systems.
* __scenarios__: the data frame of scenarios to run.
* __num.rep__: the number of replicates to simulate for each scenario.
* __num.cores__: the number of cores to use. Replicates for all scenarios are load balanced amongst the number of cores selected. That is, as cores are freed, the next replicates will be allocated to those cores, so that cores are always working until there are fewer than num.cores replicates left. If this parameter is set to a value greater than 1, progress notifications will not be printed on the console.

__NOTE:__ Some scenarios can take a lot of memory and a long time to run. Resource use scales proportional to the number of individuals (Ne * num.pops). On many systems, multiple cores will share the available memory. Thus, increasing the number of cores to use has the potential to exponentially increase memory usage and cause system crashes. If several large scenarios are in the group being run, it is suggested to use fewer cores and let the simulation run longer. These choices will be system dependant and it suggested that a few test runs be done with a small number of replicates to ensure that crashes do not occur.

***

### Upload simulation results

When the simulations are complete, there will be three new items in the working directory:

* \<label\>_scenario.replicates: a folder containig .csv files of genotypes for each scenario replicate.
* \<label\>_scenarios.csv: a .csv file of the scenario specifications.
* \<label\>_params.rdata: an R workspace file containing a a list called `params` that contains the parameters used to to run the scenarios and an element called `$run.smry` which is a data frame of the replicates, their start and stop times, the total run time, and the filename of genotypes created.

The `runEBVsim()` function also invisibly returns the same summary list contained in \<label\>_params.rdata. 

When the run is complete, compress the folder of results along with the "_params.rdata" file and upload them to the Google Drive folder [here](https://drive.google.com/open?id=1TGI2TVFnOAx0ib1GdBaL80Pwq-7ruqG1). Name the compressed file "\<label\>_results.tar.gz" (or .zip, or whatever compression algorithm you use).

***

### Contact

* submit suggestions and bug-reports: <https://github.com/ericarcher/ebvSim/issues>
* send a pull request: <https://github.com/ericarcher/ebvSim/>
* e-mail: <eric.archer@noaa.gov>
