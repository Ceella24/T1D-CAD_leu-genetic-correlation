library(MRInstruments)
library(TwoSampleMR)
library(dplyr)
library(xlsx)

ao <- available_outcomes()   

ifih1_only <- read_exposure_data(
                  filename = "/Users/priscillasaarah/Downloads/ifih1_only_snps.tsv",
                  sep = "\t",
                  snp_col = "variant_id",
                  beta_col = "beta",
                  se_col = "standard_error",
                  effect_allele_col = "effect_allele",
                  other_allele_col = "other_allele",
                  eaf_col = "effect_allele_frequency",
                  pval_col = "p_value",
                  chr_col = "chromosome",
                  pos_col = "base_pair_location"
 )

#method_list_5 <- c("mr_ivw", "mr_ivw_fe", "mr_ivw_mre", "mr_simple_median", "mr_weighted_median","mr_weighted_mode","mr_egger_regression", "mr_two_sample_ml")


mr_func <- function(outcome_code, shorthand) {
  
  require(MRInstruments)
  require(TwoSampleMR)
  require(xlsx)
  require(dplyr)
  
  filepath <- paste0("/Users/priscillasaarah/Downloads/ifih1_PS_Nov14/mr_results_all/", shorthand, ".xlsx")
  
  outcome_data <- extract_outcome_data(snps = ifih1_only$SNP, outcomes = outcome_code)
  H_data <- harmonise_data(exposure_dat = ifih1_only, outcome_dat = outcome_data)
  
  mr_results <- mr(H_data, method_list = c("mr_penalised_weighted_median", "mr_ivw", "mr_ivw_mre", "mr_ivw_fe", "mr_simple_mode", "mr_weighted_mode",
                                           "mr_simple_median", "mr_weighted_median", "mr_egger_regression", "mr_two_sample_ml", "mr_ivw_radial"))
  odds_ratio <- generate_odds_ratios(mr_results)
  plei <- mr_pleiotropy_test(H_data)
  heter <- mr_heterogeneity(H_data)
  
  write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
  write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
  write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
  write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
}

outcome_code <- c("ebi-a-GCST004618", "ebi-a-GCST004606", "ebi-a-GCST004627", "ebi-a-GCST004625", "ebi-a-GCST004629", "ebi-a-GCST004610")
shorthand <- c("BASO", "EOS", "LYM", "MONO", "NEU", "LEU" )

for (x in 1:length(shorthand)){
  mr_func(outcome_code[x], shorthand[x])
}