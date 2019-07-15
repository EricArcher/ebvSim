cat(format(Sys.time()), "Computing ldNe...\n")

snps.ne <- ldNe(sim.snps.g) %>% 
  select(stratum, Ne, param.lci, param.uci)
