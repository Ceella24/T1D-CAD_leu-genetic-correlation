library(MRInstruments)
library(TwoSampleMR)
library(dplyr)
library(xlsx)
library(RadialMR)

#ao <- available_outcomes()

mr_func <- function(name){
  
  require(MRInstruments)
  require(TwoSampleMR)
  require(dplyr)
  require(xlsx)
  require(RadialMR)
  
  
  exposure_data <- read_exposure_data(filename = paste0("~/Downloads/eQTL/010724/eqtlGen/sumstats/", name, "_eqtlgen_eqtl.tsv"), 
                                      clump = F, sep = "\t", snp_col = "SNP", beta_col = "Beta", se_col = "SE",
                                      effect_allele_col = "AssessedAllele", eaf_col = "AlleleB_all",
                                      other_allele_col = "OtherAllele", pval_col = "Pvalue", chr_col = "SNPChr", 
                                      pos_col = "SNPPos", gene_col = "GeneSymbol", samplesize_col = "NrSamples", 
                                      )
  exposure_data <- clump_data(dat = exposure_data, clump_r2 = 0.1)


  # #t1d
  # filepath <-paste0("/Users/priscillasaarah/Downloads/eQTL/010724/eqtlGen/clumped_r2=0.1/ivw/", name, "/T1D.xlsx")
  # 
  # outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = "ebi-a-GCST90014023")
  # H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
  # 
  # ivw_rad <- ivw_radial(H_data)
  # #egger_rad <- egger_radial(H_data)
  # 
  # H_data <- subset(H_data, !(SNP %in% ivw_rad$outliers$SNP)) #(| SNP %in% egger_rad$outliers$SNP))
  # 
  # #mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
  # mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
  # odds_ratio <- generate_odds_ratios(mr_results)
  # plei <- mr_pleiotropy_test(H_data)
  # heter <- mr_heterogeneity(H_data)
  # 
  # write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
  # write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
  # write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
  # write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  # # write.csv(H_data, file = filepath1, row.names = F)

  # #eos
  # filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/eqtlGen/clumped_r2=0.1/ivw/", name, "/EOS.xlsx")
  # #filepath1 <- paste0("/Users/priscillasaarah/Downloads/eQTL/010324/", name, "/EOS.csv")
  # 
  # outcome_data <- read_outcome_data(filename = "~/Downloads/eQTL/123123/sum_stats/eos_Kachuri.csv", snps = exposure_data$SNP,
  #                                   sep = ",", snp_col = "rsid", beta_col = "beta", se_col = "se",
  #                                   effect_allele_col = "effect_allele", other_allele_col = "other_allele",
  #                                   pval_col = "P_value", chr_col = "Chr")
  # 
  # H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
  # 
  # 
  # ivw_rad <- ivw_radial(H_data)
  # #egger_rad <- egger_radial(H_data)
  # 
  # H_data <- subset(H_data, !(SNP %in% ivw_rad$outliers$SNP)) #(| SNP %in% egger_rad$outliers$SNP))
  # 
  # mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
  # #mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
  # odds_ratio <- generate_odds_ratios(mr_results)
  # plei <- mr_pleiotropy_test(H_data)
  # heter <- mr_heterogeneity(H_data)
  # 
  # write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
  # write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
  # write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
  # write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  # #write.csv(H_data, file = filepath1, row.names = F)

  #leu
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/eqtlGen/clumped_r2=0.1/ivw/", name, "/LEU.xlsx")
  #filepath1 <- paste0("/Users/priscillasaarah/Downloads/eQTL/010324/", name, "/LEU.csv")

  outcome_data <- read_outcome_data(filename = "~/Downloads/eQTL/123123/sum_stats/leu_Kachuri.csv",
                                    snps = exposure_data$SNP, sep = ",", snp_col = "rsid", beta_col = "beta", se_col = "se",
                                    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                                    pval_col = "P_value", chr_col = "Chr")

  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)


  ivw_rad <- ivw_radial(H_data)
  #egger_rad <- egger_radial(H_data)

  H_data <- subset(H_data, !(SNP %in% ivw_rad$outliers$SNP)) #(| SNP %in% egger_rad$outliers$SNP))

  mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
  #mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
  odds_ratio <- generate_odds_ratios(mr_results)
  plei <- mr_pleiotropy_test(H_data)
  heter <- mr_heterogeneity(H_data)

  write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
  write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
  write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
  write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  # write.csv(H_data, file = filepath1, row.names = F)

  #lym
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/eqtlGen/clumped_r2=0.1/ivw/", name, "/LYM.xlsx")
  #filepath1 <- paste0("/Users/priscillasaarah/Downloads/eQTL/010324/", name, "/LYM.csv")

  outcome_data <- read_outcome_data(filename = "~/Downloads/eQTL/123123/sum_stats/lym_Kachuri.csv",
                                    snps = exposure_data$SNP, sep = ",", snp_col = "rsid", beta_col = "beta", se_col = "se",
                                    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                                    pval_col = "P_value", chr_col = "Chr")

  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)

  ivw_rad <- ivw_radial(H_data)
  #egger_rad <- egger_radial(H_data)

  H_data <- subset(H_data, !(SNP %in% ivw_rad$outliers$SNP)) #(| SNP %in% egger_rad$outliers$SNP))

  mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
  #mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
  odds_ratio <- generate_odds_ratios(mr_results)
  plei <- mr_pleiotropy_test(H_data)
  heter <- mr_heterogeneity(H_data)

  write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
  write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
  write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
  write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  # write.csv(H_data, file = filepath1, row.names = F)

  #mono
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/eqtlGen/clumped_r2=0.1/ivw/", name, "/MONO.xlsx")
  #filepath1 <- paste0("/Users/priscillasaarah/Downloads/eQTL/010324/", name, "/MONO.csv")

  outcome_data <- read_outcome_data(filename = "~/Downloads/eQTL/123123/sum_stats/mono_Kachuri.csv",
                                    snps = exposure_data$SNP, sep = ",", snp_col = "rsid", beta_col = "beta", se_col = "se",
                                    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                                    pval_col = "P_value", chr_col = "Chr")

  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)

  ivw_rad <- ivw_radial(H_data)
  #egger_rad <- egger_radial(H_data)

  H_data <- subset(H_data, !(SNP %in% ivw_rad$outliers$SNP)) #(| SNP %in% egger_rad$outliers$SNP))

  mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
  #mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
  odds_ratio <- generate_odds_ratios(mr_results)
  plei <- mr_pleiotropy_test(H_data)
  heter <- mr_heterogeneity(H_data)

  write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
  write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
  write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
  write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  #write.csv(H_data, file = filepath1, row.names = F)

  #baso
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/eqtlGen/clumped_r2=0.1/ivw/", name, "/BASO.xlsx")
  #filepath1 <- paste0("/Users/priscillasaarah/Downloads/eQTL/010324/", name, "/BASO.csv")

  outcome_data <- read_outcome_data(filename = "~/Downloads/eQTL/123123/sum_stats/baso_Kachuri.csv",
                                    snps = exposure_data$SNP, sep = ",", snp_col = "rsid", beta_col = "beta", se_col = "se",
                                    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                                    pval_col = "P_value", chr_col = "Chr")

  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)

  ivw_rad <- ivw_radial(H_data)
  #egger_rad <- egger_radial(H_data)

  H_data <- subset(H_data, !(SNP %in% ivw_rad$outliers$SNP)) #(| SNP %in% egger_rad$outliers$SNP))

  mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
  #mr_results <- mr(H_data, method_list =c("mr_wald_ratio"))
  odds_ratio <- generate_odds_ratios(mr_results)
  plei <- mr_pleiotropy_test(H_data)
  heter <- mr_heterogeneity(H_data)

  write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
  write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
  write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
  write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  #write.csv(H_data, file = filepath1, row.names = F)

  #neu
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/eqtlGen/unclumped/ivw/", name, "/NEU.xlsx")
  #filepath1 <- paste0("/Users/priscillasaarah/Downloads/eQTL/010324/", name, "/NEU.csv")

  outcome_data <- read_outcome_data(filename = "~/Downloads/eQTL/123123/sum_stats/neu_Kachuri.csv",
                                    snps = exposure_data$SNP, sep = ",", snp_col = "rsid", beta_col = "beta", se_col = "se",
                                    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                                    pval_col = "P_value", chr_col = "Chr")

  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)

  ivw_rad <- ivw_radial(H_data)
  #egger_rad <- egger_radial(H_data)

  H_data <- subset(H_data, !(SNP %in% ivw_rad$outliers$SNP))  #(| SNP %in% egger_rad$outliers$SNP))

  mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
  #mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
  odds_ratio <- generate_odds_ratios(mr_results)
  plei <- mr_pleiotropy_test(H_data)
  heter <- mr_heterogeneity(H_data)

  write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
  write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
  write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
  write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  #write.csv(H_data, file = filepath1, row.names = F)

  #cad
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/eqtlGen/clumped_r2=0.1/ivw/", name, "/CAD.xlsx")
  #filepath1 <- paste0("/Users/priscillasaarah/Downloads/eQTL/010324/", name, "/CAD.csv")

  outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = "ebi-a-GCST005195")

  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)

  ivw_rad <- ivw_radial(H_data)
  #egger_rad <- egger_radial(H_data)

  H_data <- subset(H_data, !(SNP %in% ivw_rad$outliers$SNP)) #(| SNP %in% egger_rad$outliers$SNP))

  mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
  #mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
  odds_ratio <- generate_odds_ratios(mr_results)
  plei <- mr_pleiotropy_test(H_data)
  heter <- mr_heterogeneity(H_data)

  write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
  write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
  write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
  write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  #write.csv(H_data, file = filepath1, row.names = F)
}

gene_name <- c("ZBTB9", "MORF4L1")
# "SMAD3", "BAK1", "CTSH", "SUOX", "SH2B3", "BACH2", "CLEC16A", "IL2RA", "PTPN22", "IFIH1", "IKZF4",
for (i in 1:length(gene_name)) {
  mr_func(gene_name[i])
}

