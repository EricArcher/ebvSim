rm(list = ls())
# run next two lines if you don't have the latest devel version of strataG
if(!require(devtools)) install.packages("devtools")
install_github("ericarcher/strataG", ref = "redevel2018", dependencies = TRUE)
library(strataG)

ebvSim <- function(genetics, label, num.sim, ploidy = 2) {
  p <- fscWrite(
    demes = fscSettingsDemes(
      fscDeme(500, 100),
      fscDeme(500, 100),
      fscDeme(500, 100), 
      ploidy = ploidy
    ),
    migration = fscSettingsMigration(
      matrix(c(0, 0.5, 0.005, 0.5, 0, 0.0005, 0.005, 0.0005, 0), ncol = 3)
    ),
    events = fscSettingsEvents(fscEvent(2000, 2, 1), fscEvent(2000, 1, 0)),
    genetics = genetics,
    label = label
  )
  fscRun(p, num.sim = num.sim)
}

snp.p <- ebvSim(
  fscSettingsGenetics(fscBlock_snp(10, 1e-3), num.chrom = 1000),
  label = "ebvSim.snps",
  num.sim = 1
)

dna.msat.p <- ebvSim(
  fscSettingsGenetics(fscBlock_dna(1000, 1e-6), fscBlock_microsat(15, 5e-4)),
  label = "ebvSim.dna.msat",
  num.sim = 1
)

save.image("ebvSim ws.rdata")