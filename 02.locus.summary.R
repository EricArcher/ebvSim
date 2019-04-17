snp.het <- left_join(
  heterozygosity(sim.snps.g, by.strata = TRUE, type = "observed"),
  heterozygosity(sim.snps.g, by.strata = TRUE, type = "expected"),
  by = c("stratum", "locus")
) %>% 
  group_by(stratum) %>% 
  summarize(
    mean.obs = mean(obsvd.het),
    median.obs = median(obsvd.het),
    LCI.obs = unname(quantile(obsvd.het, 0.025)),
    UCI.obs = unname(quantile(obsvd.het, 0.975)),
    mean.exp = mean(exptd.het),
    median.exp = median(exptd.het),
    LCI.exp = unname(quantile(exptd.het, 0.025)),
    UCI.exp = unname(quantile(exptd.het, 0.975))
  ) %>% 
  ungroup() %>% 
  as.data.frame(stringsAsFactors = FALSE)
    
  