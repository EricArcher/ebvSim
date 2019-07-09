makeMigSettings <- function(mig.rate, num.pops, 
                        type = c("island", "stepping.stone")) {
  stopifnot(require(strataG))
  
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
  fscSettingsMigration(mig.mat)
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


ebvSim <- function(scenarios, genetics, label, num.sim, ploidy = 2, 
                   google.drive.id = NULL) {
  stopifnot(require(strataG))
  stopifnot(require(googledrive))
  
  fname <- paste0(label, "_scenarios.csv")
  out.folder <- paste0(label, "_scenario.replicates")
  
  google.rep.dir <- NULL
  if(!is.null(google.drive.id)) {
    label.pattern = paste0(fname, "|", out.folder)
    drive.files <- drive_ls(
      as_id(google.drive.id), 
      pattern = label.pattern, 
      verbose = FALSE
    )
    if(nrow(drive.files) > 0) {
      stop(
        "the scenario file and/or replicate folder for '", label, 
        "' already exists on this Google Drive.\n",
        "Delete them before running the sim with this label."
      )
    } else {
      google.rep.dir <- drive_mkdir(
        out.folder, 
        parent = as_id(google.drive.id), 
        verbose = FALSE
      )
    }
  }
    
  if(!dir.exists(out.folder)) dir.create(out.folder)
  write.csv(scenarios, file = fname, row.names = FALSE)
  if(!is.null(google.drive.id)) {
    drive_upload(fname, path = as_id(google.drive.id), verbose = FALSE)
  }
  
  for(sc.i in 1:nrow(scenarios)) {
    sc <- scenarios[sc.i, ]
    
    deme.list <- lapply(1:sc$num.pops, function(i) {
      fscDeme(deme.size = sc$Ne, sample.size = sc$num.samples)
    })
    deme.list$ploidy <- ploidy
    
    p <- fscWrite(
      demes = do.call(fscSettingsDemes, deme.list),
      migration = if(sc$num.pops > 1) {
        makeMigSettings(sc$mig.rate, sc$num.pops, sc$mig.type) 
      } else NULL,
      events = makeEventSettings(sc$dvgnc.time, sc$num.pops),
      genetics = genetics,
      label = label
    )
    
    p <- fscRun(p, num.sim = num.sim)
    
    for(sim.i in 1:num.sim) {
      gen.data <- fscReadArp(p, sim = c(1, sim.i), drop.mono = TRUE)
      fname <- paste0(
        label, "_", 
        "scenario.", sc$scenario, 
        "_replicate.", sim.i, 
        ".csv"
      )
      
      out.name <- file.path(out.folder, fname)
      write.csv(gen.data, file = out.name, row.names = FALSE)
      if(!is.null(google.rep.dir)) {
        drive_upload(
          out.name, 
          path = as_id(google.rep.dir), 
          name = fname,
          verbose = FALSE
        )
      }
    }
  }
  
  invisible(out.folder)
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

killExcess <- function(Rland, n) {
  library(magrittr)
  to.kill <- tapply(
    1:nrow(Rland$individuals), 
    Rland$individuals[, 1], 
    function(i) if(length(i) > n) sample(i, length(i) - n) else NULL
  ) %>% 
    unname() %>% 
    unlist()
  if(length(to.kill) > 0) Rland$individuals <- Rland$individuals[-to.kill, ]
  Rland
}


loadLandscape <- function(sc, AlleleFreqs, num.gens) {
  library(magrittr)
  library(rmetasim)
  
  localS <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2)
  localR <- matrix(c(0, 0, 1.2, 0), nrow = 2, ncol = 2)
  localM <- matrix(c(0, 0, 0, 1.2), nrow = 2, ncol = 2)
  S <- M <- matrix(0, nrow = sc$num.pops * 2, ncol = sc$num.pops * 2)
  diag(S) <- diag(M) <- 1
  R <- rmetasim::landscape.mig.matrix(
    h = nrow(sc$mig.mat), s = 2, mig.model = "custom", R.custom = sc$mig.mat
  )$R
  
  Rland <- rmetasim::landscape.new.empty() %>% 
    rmetasim::landscape.new.intparam(
      h = sc$num.pops, s = 2, cg = 0, ce = 0, totgen = num.gens + 1
    ) %>% 
    rmetasim::landscape.new.switchparam() %>% 
    rmetasim::landscape.new.floatparam() %>% 
    rmetasim::landscape.new.local.demo(localS, localR, localM) %>% 
    rmetasim::landscape.new.epoch(R = R, carry = rep(sc$ne, sc$num.pops))
  
  # just make loci that have the correct type and ploidy
  for(i in 1:length(AlleleFreqs)) {
    Rland <- rmetasim::landscape.new.locus(
      Rland, type = 2, ploidy = 2, mutationrate = 0,
      transmission = 0, numalleles = 2, allelesize = 1,
      frequencies = AlleleFreqs[[i]], states = rownames(AlleleFreqs[[i]])
    )
  }

  rmetasim::landscape.new.individuals(Rland, rep(c(sc$ne, 0), sc$num.pops))
}