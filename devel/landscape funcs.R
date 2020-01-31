killExcess <- function(rl, n) {
  to.keep <- tapply(1:nrow(rl$individuals), rl$individuals[, 1], function(i) {
    if(length(i) > n) sample(i, n) else i
  })
  to.keep <- sort(unlist(unname(to.keep)))
  rl$individuals <- rl$individuals[to.keep, , drop = FALSE]
  rl
}

landscape.setallelefreq.2 <- function (rland, af = NULL, states = TRUE) {
  if (is.null(af)) {
    cat("specify at least some allele frequencies\n")
  } else {
    locposition <- as.character(landscape.locusvec(rland))
    possLoc <- as.character(1:length(rland$loci))
    possStgs <- as.character(0:((rland$intparam$habitats * rland$intparam$stages) - 1))
    democol <- landscape.democol()
    if (length(which(!(names(af) %in% possStgs))) > 0) {
      stop(paste("some stages in af do not exist in landscape"))
    }
    else {
      for (s in names(af)) {
        s.inds <- rland$individuals[, 1] == as.numeric(s)
        inds <- rland$individuals[s.inds, ]
        n <- dim(inds)[1]
        if(!all(names(af[[s]]) %in% possLoc)) {
          stop("a locus name is not found in the landscape")
        }
        for (l in names(af[[s]])) {
          av <- af[[s]][[l]]
          if (sum(av) > 1) stop("allele freq vector sums to > 1")
          aindex <- if (states) {
            ai <- landscape.locus.states(rland, as.numeric(l), do.check = F)
            ai$aindex[ai$state %in% names(av)]
          } else names(av)
          lcols <- which(locposition == l) + democol
          inds[, lcols] <- as.numeric(
            sample(aindex, n * length(lcols), replace = T, prob = av)
          )
        }
        rland$individuals[s.inds, ] <- inds
      }
    }
  }
  rland
}


landscape.setpopfreq.2 <- function (rland, af = NULL, states = TRUE) {
  posspops <- sort(unique(landscape.populations(rland)))
  s = rland$intparam$stages
  h = rland$intparam$habitats
  newaf <- lapply(names(af), function(p) {
    pnum <- as.numeric(p)
    stages <- ((pnum - 1) * s):((pnum - 1) * s + s - 1)
    lapply(stages, function(s) stats::setNames(list(af[[p]]), s))
  })
  newaf <- unlist(newaf, recursive = F)
  landscape.setallelefreq.2(rland, newaf, states)
}

landscape.ind.freq.2 <- function (Rland, include.states = TRUE) {
  l <- Rland
  ploidy <- landscape.ploidy(l)
  aml <- vector("list", length(ploidy))
  for (loc in 1:length(aml)) {
    genos <- landscape.locus(l, loc)[, -(1:landscape.democol())]
    loc.ploidy <- ploidy[loc]
    if (l$loci[[loc]]$type != 253) {
      lst <- landscape.locus.states(l, loc)
      names(lst$state) <- lst$aindex
      if (loc.ploidy == 2) {
        genos[, 1] <- unname(lst$state[as.character(genos[, 1])])
        genos[, 2] <- unname(lst$state[as.character(genos[, 2])])
      } else {
        genos <- unname(lst$state[as.character(genos)])
      }
    }
    unique.genos <- sort(as.character(unique(as.vector(genos))))
    aml[[loc]] <- sapply(unique.genos, function(x) {
      if (loc.ploidy == 2) {
        (as.character(genos[, 1]) == x) + (as.character(genos[, 2]) == x)
      } else {
        as.character(genos) == x
      }
    }) / loc.ploidy
  }
  do.call(cbind, aml)
}


landscape.make.genind.2 <- function (Rland) {
  tab <- landscape.ind.freq.2(Rland) * 2
  dimnames(tab) <- list(
    rownames = 1:dim(tab)[1], 
    colnames = landscape.freq.locnames(Rland)
  )
  gi <- adegenet::genind(tab, pop = as.factor(landscape.populations(Rland)), ploidy = 2)
  gi[, loc = which(landscape.ploidy(Rland) > 1)]
}