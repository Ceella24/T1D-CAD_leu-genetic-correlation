library(MRInstruments)
library(TwoSampleMR)
library(dplyr)
library(xlsx)
library(RadialMR)

#ao <- available_outcomes()
#ivw_rad <- ivw_radial(H_data)
#egger_rad <- egger_radial(H_data)

#H_data <- subset(H_data, !(SNP %in% ivw_rad$outliers$SNP)) #(| SNP %in% egger_rad$outliers$SNP))

mr_func <- function(name, small_out_name, big_out_name){
  
  require(MRInstruments)
  require(TwoSampleMR)
  require(dplyr)
  require(xlsx)
  require(RadialMR)
  
  count = 0
  
  exposure_data <- read_exposure_data(filename = paste0("~/Downloads/eQTL/010724/GTEx/sumstats/", name, "_GTEx_eqtl.tsv"), 
                                      clump = F, sep = "\t", snp_col = "SNP.Id", beta_col = "slope", se_col = "slope_se",
                                      effect_allele_col = "effect_allele", eaf_col = "maf",
                                      other_allele_col = "other_allele", pval_col = "pval_nominal", 
                                      gene_col = "Gene.Symbol", phenotype_col = "Gene.Symbol", id_col = "Tissue"
  )
  exposure_data <- clump_data(dat = exposure_data, clump_r2 = 0.1)
  
  #t1d
  count = 0
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/GTEx/clumped_r2=0.1/nothing_removed/", name, "/T1D.xlsx")

  outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = "ebi-a-GCST90014023")

  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)

  for (i in H_data$mr_keep){
    if (i == TRUE){
    count = count + 1
    }
  }
  if (count == 1){
    mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
    odds_ratio <- generate_odds_ratios(mr_results)

    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)

  }else if (count > 1){

    mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
    odds_ratio <- generate_odds_ratios(mr_results)
    plei <- mr_pleiotropy_test(H_data)
    heter <- mr_heterogeneity(H_data)

    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
    write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  }

  #cad
  count = 0
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/GTEx/clumped_r2=0.1/nothing_removed/", name, "/CAD.xlsx")

  outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = "ebi-a-GCST005195")
  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)

  for (i in H_data$mr_keep){
    if (i == TRUE){
      count = count + 1
    }
  }
  if (count == 1){
    mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
    odds_ratio <- generate_odds_ratios(mr_results)

    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)

  }else if (count > 1){

    mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
    odds_ratio <- generate_odds_ratios(mr_results)
    plei <- mr_pleiotropy_test(H_data)
    heter <- mr_heterogeneity(H_data)

    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
    write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  }
   
  #eos  
  count = 0
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/GTEx/clumped_r2=0.1/nothing_removed/", name, "/EOS.xlsx")
  
  outcome_data <- read_outcome_data(filename = "~/Downloads/Kachuri_sumstats/eos_Kachuri.tsv", snps = exposure_data$SNP,
                                    sep = "\t", snp_col = "rsid", beta_col = "beta", se_col = "se",
                                    effect_allele_col = "risk_allele", other_allele_col = "other_allele",
                                    pval_col = "P_value", chr_col = "Chr")
  outcome_data$outcome = "eos"

  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
  
  for (i in H_data$mr_keep){
    if (i == TRUE){
    count = count + 1
    }
  }
  if (count == 1){
     mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
    odds_ratio <- generate_odds_ratios(mr_results)
    
    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
    
  }else if (count > 1){
    
    mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
    odds_ratio <- generate_odds_ratios(mr_results)
    plei <- mr_pleiotropy_test(H_data)
    heter <- mr_heterogeneity(H_data)
    
    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
    write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  }
  
  #leu
  count = 0
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/GTEx/clumped_r2=0.1/nothing_removed/", name, "/LEU.xlsx")
  
  outcome_data <- read_outcome_data(filename = "~/Downloads/Kachuri_sumstats/leu_Kachuri.tsv",
                                    snps = exposure_data$SNP, sep = "\t", snp_col = "rsid", beta_col = "beta", se_col = "se",
                                    effect_allele_col = "risk_allele", other_allele_col = "other_allele",
                                    pval_col = "P_value", chr_col = "Chr")
  outcome_data$outcome = "leu"
 
  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
  
  for (i in H_data$mr_keep){
    if (i == TRUE){
      count = count + 1
    }
  }
  if (count == 1){
    mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
    odds_ratio <- generate_odds_ratios(mr_results)
    
    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
    
  }else if (count > 1){
    
    mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
    odds_ratio <- generate_odds_ratios(mr_results)
    plei <- mr_pleiotropy_test(H_data)
    heter <- mr_heterogeneity(H_data)
    
    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
    write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  }
  
  #lym
  count = 0
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/GTEx/clumped_r2=0.1/nothing_removed/", name, "/LYM.xlsx")
  
  outcome_data <- read_outcome_data(filename = "~/Downloads/Kachuri_sumstats/lym_Kachuri.tsv",
                                    snps = exposure_data$SNP, sep = "\t", snp_col = "rsid", beta_col = "beta", se_col = "se",
                                    effect_allele_col = "risk_allele", other_allele_col = "other_allele",
                                    pval_col = "P_value", chr_col = "Chr")
  outcome_data$outcome = "lym"
  
  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
  
  for (i in H_data$mr_keep){
    if (i == TRUE){
      count = count + 1
    }
  }
  if (count == 1){
    mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
    odds_ratio <- generate_odds_ratios(mr_results)
    
    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
    
  }else if (count > 1){
    
    mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
    odds_ratio <- generate_odds_ratios(mr_results)
    plei <- mr_pleiotropy_test(H_data)
    heter <- mr_heterogeneity(H_data)
    
    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
    write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  }
  
  
  #mono
  count = 0
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/GTEx/clumped_r2=0.1/nothing_removed/", name, "/MONO.xlsx")
  
  outcome_data <- read_outcome_data(filename = "~/Downloads/Kachuri_sumstats/mono_Kachuri.tsv",
                                    snps = exposure_data$SNP, sep = "\t", snp_col = "rsid", beta_col = "beta", se_col = "se",
                                    effect_allele_col = "risk_allele", other_allele_col = "other_allele",
                                    pval_col = "P_value", chr_col = "Chr")
  outcome_data$outcome = "mono"
  
  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
  
  for (i in H_data$mr_keep){
    if (i == TRUE){
      count = count + 1
    }
  }
  if (count == 1){
    mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
    odds_ratio <- generate_odds_ratios(mr_results)
    
    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
    
  }else if (count > 1){
    
    mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
    odds_ratio <- generate_odds_ratios(mr_results)
    plei <- mr_pleiotropy_test(H_data)
    heter <- mr_heterogeneity(H_data)
    
    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
    write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  }
  
  #baso
  count = 0
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/GTEx/clumped_r2=0.1/nothing_removed/", name, "/BASO.xlsx")
  
  outcome_data <- read_outcome_data(filename = "~/Downloads/Kachuri_sumstats/baso_Kachuri.tsv",
                                    snps = exposure_data$SNP, sep = "\t", snp_col = "rsid", beta_col = "beta", se_col = "se",
                                    effect_allele_col = "risk_allele", other_allele_col = "other_allele",
                                    pval_col = "P_value", chr_col = "Chr")
  outcome_data$outcome = "baso"
  
  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
  
  for (i in H_data$mr_keep){
    if (i == TRUE){
      count = count + 1
    }
  }
  if (count == 1){
    mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
    odds_ratio <- generate_odds_ratios(mr_results)
    
    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
    
  }else if (count > 1){
    
    mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
    odds_ratio <- generate_odds_ratios(mr_results)
    plei <- mr_pleiotropy_test(H_data)
    heter <- mr_heterogeneity(H_data)
    
    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
    write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  }
  
  #neu
  count = 0
  filepath <- paste0("/Users/priscillasaarah/Downloads/eQTL/010724/GTEx/clumped_r2=0.1/nothing_removed/", name, "/NEU.xlsx")
  
  outcome_data <- read_outcome_data(filename = "~/Downloads/Kachuri_sumstats/neu_Kachuri.tsv",
                                    snps = exposure_data$SNP, sep = "\t", snp_col = "rsid", beta_col = "beta", se_col = "se",
                                    effect_allele_col = "risk_allele", other_allele_col = "other_allele",
                                    pval_col = "P_value", chr_col = "Chr")
  outcome_data$outcome = "neu"
  
  H_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
  
  for (i in H_data$mr_keep){
    if (i == TRUE){
      count = count + 1
    }
  }
  if (count == 1){
    mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))
    odds_ratio <- generate_odds_ratios(mr_results)
    
    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
    
  }else if (count > 1){
    
    mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
    odds_ratio <- generate_odds_ratios(mr_results)
    plei <- mr_pleiotropy_test(H_data)
    heter <- mr_heterogeneity(H_data)
    
    write.xlsx(odds_ratio, file = filepath, sheetName = "odds_ratios", row.names = FALSE)
    write.xlsx(heter, file = filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
    write.xlsx(plei, file = filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
    write.xlsx(H_data, file = filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
  }

}
#"IL2RA_sp", "MORF4L1_LCL", "SMAD3_subadi", "CTSH_wb", "CXCL12_pan", 
# "IKZF4_LCL", "CTRB1_pan", "BAK1_pan", "BAK1_wb", "CTRB2_pan", 
# "SH2B3_wb", "SUOX_pan", "SUOX_wb","HORMAD2_subadi", "PTPN22_wb", 
# "TCP11_sp", "BACH2_pan", "PTPN22_pan", "SH2B3_pan", 
# "IKZF4_pan", "CTSH_LCL", "MORF4L1_wb","CLEC16A_wb", "IL2RA_pan", 
#"CTSH_pan", "IFIH1_pan", "SMAD3_sp",
gene <- c("SH2B3_LCL", "MORF4L1_pan", "IFIH1_wb", "TCP11_LCL", 
          "IFIH1_increase_wb", "CFDP1_wb", "CFDP1_pan", "CFDP1_LCL", "CFDP1_increase_LCL", 
          "SMAD3_increase_sp", "BACH2_increase_pan", "SH2B3_increase_LCL")
#  
for (i in 1:length(gene)){
  mr_func(gene[i])
}
