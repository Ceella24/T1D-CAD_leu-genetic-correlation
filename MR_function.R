library(MRInstruments)
library(TwoSampleMR)
library(dplyr)
library(xlsx)

#gwas_list <- c("ebi-a-GCST005195", "ebi-a-GCST90014023", "ebi-a-GCST90002407", "ebi-a-GCST90002379", "ebi-a-GCST90002381", "ebi-a-GCST90002316", "ebi-a-GCST90002340", "ebi-a-GCST90002398")
#shorthand <- c("CAD", "T1D", "leu", "baso", "eos", "lym", "mono", "neu")

mr_func <- function(abbrev_expo, abbrev_out, expo_code, out_code){
  
  require(MRInstruments)
  require(TwoSampleMR)
  require(dplyr)
  require(xlsx)
  
  ao <- available_outcomes()
  
  filepath <- paste0("/Users/priscillasaarah/Downloads/MR_Oct/reverse_causality/", abbrev_expo, "/", abbrev_out, ".xlsx")
  #filepath1 <- paste0("/Users/priscillasaarah/Downloads/MR_Nov23/forward_causality/CAD/", shorthand, ".csv")
  
  expo_gwas <- extract_instruments(outcomes = expo_code)
  #exposure_data <- format_data(exposure_gwas)
  #exposure_data <-clump_data(expo_gwas, clump_r2 = 0.001)
  
  outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = out_code)
  
  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
  
  mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
  odds_ratio <- generate_odds_ratios(mr_results)
  plei <- mr_pleiotropy_test(H_data)
  heter <- mr_heterogeneity(H_data)
  #press <- run_mr_presso(H_data, NbDistribution = 1500)
  
  write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
  write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
  write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
  write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
}

#forward causality
  out_code <- c("ebi-a-GCST005195", "ebi-a-GCST90014023")
  expo_abrev <- c("LEU", "BASO", "EOS", "LYM", "MONO", "NEU")
  out_abrev <- c("CAD", "T1D")
  expo_code <- c("ebi-a-GCST90002407", "ebi-a-GCST90002379", "ebi-a-GCST90002381", "ebi-a-GCST90002316", "ebi-a-GCST90002340", "ebi-a-GCST90002398")
    
  #T1D as exposure
  #CAD as exposure
  for (x in 1:length(expo_code)){
   for (y in 1:length(out_abrev)){
     mr_func(expo_abrev[x], out_abrev[y], expo_code[x], out_code[y])
   }
  }
   
    
