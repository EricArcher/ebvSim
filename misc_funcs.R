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


makeScenarios <- function(num.pops, Ne, num.samples, mig.rate, 
                          mig.type = c("island", "stepping.stone"),
                          dvgnc.time) {
  stopifnot(require(tidyverse))
  
  expand.grid(
    num.pops = num.pops,
    Ne = Ne,
    num.samples = num.samples,
    mig.rate = mig.rate,
    mig.type = match.arg(mig.type),
    dvgnc.time = dvgnc.time,
    stringsAsFactors = FALSE
  ) %>% 
    mutate(scenario = 1:n()) %>% 
    select(scenario, everything())
}

makeEventSettings <- function(dvgnc.time, num.pops) {
  stopifnot(require(strataG))
  if(num.pops == 1) return(NULL)
  pop.pairs <- t(combn(num.pops, 2) - 1)
  pop.pairs <- pop.pairs[pop.pairs[, 1] == 0, , drop = FALSE]
  do.call(
    fscSettingsEvents, 
    lapply(
      1:nrow(pop.pairs),
      function(i) fscEvent(dvgnc.time, pop.pairs[i, 2], pop.pairs[i, 1])
    )
  )
}


scenariosFname <- function(label) paste0(label, "_scenarios.csv")

repFolder <- function(label) paste0(label, "_scenario.replicates")

repFname <- function(label, scenario, replicate) {
  paste0(label, "_scenario.", scenario, "_replicate.", replicate, ".csv")
}


writeScenarios <- function(label, scenarios, google.drive.id = NULL) {
  fname <- scenariosFname(label)
  out.folder <- repFolder(label)
  google.rep.dir <- NULL
  if(!is.null(google.drive.id)) {
    label.pattern = paste0(fname, "|", out.folder)
    id <- as_id(google.drive.id)
    drive.files <- drive_ls(id, pattern = label.pattern, verbose = FALSE)
    if(nrow(drive.files) > 0) {
      stop(
        "The scenario file and/or replicate folder for '", label, 
        "' already exists on this Google Drive. ",
        "Either delete them or change the label for the run."
      )
    } else {
      google.rep.dir <- drive_mkdir(out.folder, parent = id, verbose = FALSE)
    }
  }
  
  if(!dir.exists(out.folder)) dir.create(out.folder)
  write.csv(scenarios, file = fname, row.names = FALSE)
  if(!is.null(google.drive.id)) {
    drive_upload(fname, path = as_id(google.drive.id), verbose = FALSE)
  }
  
  list(out = out.folder, google = google.rep.dir)
}


runFscSim <- function(sc, ploidy, num.rep) {
  stopifnot(require(strataG))
  stopifnot(require(magrittr))
  deme.list <- lapply(1:sc$num.pops, function(i) {
    fscDeme(deme.size = sc$Ne, sample.size = sc$num.samples)
  })
  deme.list$ploidy <- ploidy
  
  p <- fscWrite(
    demes = do.call(fscSettingsDemes, deme.list),
    migration = if(sc$num.pops > 1) {
      mig.mat <- makeMigMat(sc$mig.rate, sc$num.pops, sc$mig.type) 
      fscSettingsMigration(mig.mat)
    } else NULL,
    events = makeEventSettings(sc$dvgnc.time, sc$num.pops),
    genetics = genetics,
    label = label
  ) %>% 
    fscRun(num.sims = num.rep)
}

macFreqs <- function(mac) {
  maf <- mean(mac == 0) + (mean(mac == 1) / 2)
  c('1' = maf, '2' = 1 - maf)
}


calcFreqs <- function(snps) {
  stopifnot(require(magrittr))
  stopifnot(require(strataG))
  
  mac.df <- snps %>% 
    df2gtypes(ploidy = 2) %>% 
    as.data.frame(coded = T) %>% 
    select(-id) %>% 
    gather(locus, mac, -stratum) %>% 
    mutate(
      stratum = as.numeric(factor(stratum)),
      locus = as.numeric(factor(locus))
    )
  
  list(
    global = lapply(split(mac.df, mac.df$locus), function(loc.df) {
      macFreqs(loc.df$mac)
    }),
    pop = lapply(split(mac.df, mac.df$stratum), function(st.df) {
      lapply(split(st.df, st.df$locus), function(loc.df) macFreqs(loc.df$mac))
    })
  )
}


loadLandscape <- function(freqs, sc) {
  stopifnot(require(magrittr))
  stopifnot(require(rmetasim))
  
  cat(format(Sys.time()), "running rmetasim\n")
  
  localS <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2)
  localR <- matrix(c(0, 0, 1.2, 0), nrow = 2, ncol = 2)
  localM <- matrix(c(0, 0, 0, 1.2), nrow = 2, ncol = 2)
  S <- M <- matrix(0, nrow = sc$num.pops * 2, ncol = sc$num.pops * 2)
  diag(S) <- diag(M) <- 1
  R <- rmetasim::landscape.mig.matrix(
    h = nrow(sc$mig.mat), 
    s = 2, 
    mig.model = "custom", 
    R.custom = makeMigMat(sc$mig.rate, sc$num.pops, sc$mig.type)
  )$R
  
  Rland <- rmetasim::landscape.new.empty() %>% 
    rmetasim::landscape.new.intparam(
      h = sc$num.pops, s = 2, cg = 0, ce = 0, totgen = 3
    ) %>% 
    rmetasim::landscape.new.switchparam() %>% 
    rmetasim::landscape.new.floatparam() %>% 
    rmetasim::landscape.new.local.demo(localS, localR, localM) %>% 
    rmetasim::landscape.new.epoch(R = R, carry = rep(sc$Ne, sc$num.pops)) 
  
  for(i in 1:length(freqs$global)) {
    Rland <- rmetasim::landscape.new.locus(
      Rland, type = 2, ploidy = 2, mutationrate = 0,
      transmission = 0, numalleles = 2, allelesize = 1,
      frequencies = freqs$global[[i]], states = names(freqs$global[[i]])
    )
  }
  
  Rland %>% 
    rmetasim::landscape.new.individuals(rep(c(sc$Ne, 0), sc$num.pops)) %>% 
    rmetasim::landscape.setpopfreq(freqs$pop) %>% 
    rmetasim::landscape.simulate(numit = 2)
}


ebvSim <- function(scenarios, genetics, label, num.rep, ploidy = 2, 
                   google.drive.id = NULL, run.rmetasim = TRUE) {
  stopifnot(require(strataG))
  stopifnot(require(magrittr))
  stopifnot(require(googledrive))
  
  folders <- writeScenarios(label, scenarios, google.drive.id)
  rep.id <- if(!is.null(folders$google)) as_id(folders$google) else NULL
               
  for(sc.i in 1:nrow(scenarios)) {
    cat(format(Sys.time()), "---- Scenario", sc.i, "----\n")
    
    sc <- scenarios[sc.i, ]
    p <- runFscSim(sc, ploidy, num.rep)
    
    for(sim.i in 1:num.rep) {
      gen.data <- fscReadArp(p, sim = c(1, sim.i), drop.mono = TRUE)
      if(run.rmetasim) {
        gen.data <- gen.data %>% 
          calcFreqs() %>% 
          loadLandscape(sc) %>% 
          rmetasim::landscape.make.genind() %>% 
          genind2gtypes() %>% 
          as.data.frame()
      }
      
      fname <- repFname(label, sc$scenario, sim.i)
      out.name <- file.path(folders$out, fname)
      write.csv(gen.data, file = out.name, row.names = FALSE)
      if(!is.null(rep.id)) {
        drive_upload(out.name, path = rep.id, name = fname, verbose = FALSE)
      }
    }
  }
  
  invisible(folders)
}


loadScenarios <- function(label, folder = getwd()) {
  fname <- scenariosFname(label)
  read.csv(file.path(folder, fname), stringsAsFactors = FALSE)
}


loadGenotypes <- function(label, scenario, replicate, folder = getwd()) {
  rep.folder <- repFolder(label)
  fname <- file.path(rep.folder, repFname(label, scenario, replicate))
  read.csv(
    file.path(folder, fname), 
    colClasses = "character", 
    stringsAsFactors = FALSE
  )
}


availRuns <- function(drive.id) {
  stopifnot(require(googledrive))
  
  contents <- drive_ls(as_id(drive.id))
  i <- grepl("_scenario.replicates", contents$name) &
    is_folder(contents)
  i <- which(i)
  gsub("_scenario.replicates", "", contents$name[i])
}


downloadRun <- function(label, drive.id, out.folder = NULL) {
  stopifnot(require(googledrive))
  
  if(is.null(out.folder)) out.folder <- file.path(tempdir(), label)
  if(!dir.exists(out.folder)) dir.create(out.folder)
  contents <- drive_ls(as_id(drive.id))
  i <- grep(label, contents$name)
  
  # download csv
  csv <- grep(".csv$", contents$name[i]) 
  csv.fname <- contents$name[csv]
  cat(csv.fname, "\n")
  drive_download(
    as_id(contents$id[csv]), 
    path = file.path(out.folder, csv.fname), 
    overwrite = TRUE, 
    verbose = FALSE
  )
  
  # download scenario replicates
  folder <- grepl("_scenario.replicates", contents$name[i]) &
    is_folder(contents[i, ])
  folder <- which(folder)
  rep.folder <- file.path(out.folder, contents$name[folder])
  if(!dir.exists(rep.folder)) dir.create(rep.folder)
  fnames <- drive_ls(as_id(contents$id[folder]))
  for(f in 1:nrow(fnames)) {
    rep.fname <- fnames$name[f]
    cat(f, "/", nrow(fnames), " : ", rep.fname, "\n", sep = "")
    drive_download(
      as_id(fnames$id[f]), 
      path = file.path(rep.folder, rep.fname), 
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  
  out.folder
}


# killExcess <- function(Rland, n) {
#   library(magrittr)
#   to.kill <- tapply(
#     1:nrow(Rland$individuals), 
#     Rland$individuals[, 1], 
#     function(i) if(length(i) > n) sample(i, length(i) - n) else NULL
#   ) %>% 
#     unname() %>% 
#     unlist()
#   if(length(to.kill) > 0) Rland$individuals <- Rland$individuals[-to.kill, ]
#   Rland
# }
