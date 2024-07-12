library(MRInstruments)
library(TwoSampleMR)
library(dplyr)
library(xlsx)
library(RadialMR)
 

#ivw and ivw_egger function implementations 
het_rem <- function(gene_name, out_name){  
  require(MRInstruments)
  require(TwoSampleMR)
  require(dplyr)
  require(xlsx)
  require(RadialMR)
  
  ivw_filepath <-paste0("~/Downloads/eQTL/010724/GTEx/clumped_r2=0.1/ivw/", gene_name, "/", out_name, ".xlsx")
  ivw_egger_filepath <-paste0("~/Downloads/eQTL/010724/GTEx/clumped_r2=0.1/ivw_egger/", gene_name, "/", out_name, ".xlsx")
  
  H_data_path <- paste0("~/Downloads/eQTL/010724/GTEx/H_data_clum_csv/", gene_name, "_", out_name, "_H_data.csv")
  read_excel <- paste0("~/Downloads/eQTL/010724/GTEx/clumped_r2=0.1/nothing_removed/", gene_name, "/", out_name, ".xlsx")
  
  H <- read.xlsx(file = read_excel, sheetName = "H_data")
  
  write.csv(H, H_data_path, row.names = F)
  H_data <- read.csv(file = H_data_path, header = T)
  
  ivw_rad <- ivw_radial(H_data)
  egger_rad <- egger_radial(H_data)
  count = 0
  print(c(gene_name, out_name))
  for (n in ivw_rad$data$Outliers){
    if (n == "Variant"){
      next
    }
    H_data <- subset(H_data, !(SNP %in% ivw_rad$outliers$SNP))
    
    for (i in H_data$mr_keep){
      if (i == TRUE){
        count = count + 1
      }
    }
    if (count == 1){
      mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))  
      odds_ratio <- generate_odds_ratios(mr_results)
      
      write.xlsx(odds_ratio, file = ivw_filepath, sheetName = "odds_ratios", row.names = FALSE)
      write.xlsx(H_data, file = ivw_filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
      
    }else if (count > 1){
      
      mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
      odds_ratio <- generate_odds_ratios(mr_results)
      plei <- mr_pleiotropy_test(H_data)
      heter <- mr_heterogeneity(H_data)
      
      write.xlsx(odds_ratio, file = ivw_filepath, sheetName = "odds_ratios", row.names = FALSE)
      write.xlsx(heter, file = ivw_filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
      write.xlsx(plei, file = ivw_filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
      write.xlsx(H_data, file = ivw_filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
    }
    
    count = 0
    for (m in egger_rad$data$Outliers){
      if (m == "Variant"){
        next
      }

      H_data <- subset(H_data, !(SNP %in% ivw_rad$outliers$SNP | SNP %in% egger_rad$outliers$SNP))
      
      for (i in H_data$mr_keep){
        if (i == TRUE){
          count = count + 1
        }
      }
      if (count == 1){
        mr_results <- mr(H_data, method_list = c("mr_wald_ratio"))  
        odds_ratio <- generate_odds_ratios(mr_results)
        
        write.xlsx(odds_ratio, file = ivw_egger_filepath, sheetName = "odds_ratios", row.names = FALSE)
        write.xlsx(H_data, file = ivw_egger_filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
        
      }else if (count > 1){
        
        mr_results <- mr(H_data, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"))
        odds_ratio <- generate_odds_ratios(mr_results)
        plei <- mr_pleiotropy_test(H_data)
        heter <- mr_heterogeneity(H_data)
        
        write.xlsx(odds_ratio, file = ivw_egger_filepath, sheetName = "odds_ratios", row.names = FALSE)
        write.xlsx(heter, file = ivw_egger_filepath, sheetName = "heterogeneity_test", append = TRUE, row.names = FALSE)
        write.xlsx(plei, file = ivw_egger_filepath, sheetName = "pleiotropy_test", append = TRUE, row.names = FALSE)
        write.xlsx(H_data, file = ivw_egger_filepath, sheetName = "H_data", append = TRUE, row.names = FALSE)
      }
      break
    }
    break
  }
}
 
outcomes <- c("T1D", "CAD", "BASO", "EOS", "LEU", "NEU", "LYM", "MONO") 
   
genes <- c("CFDP1_wb", "CFDP1_pan", "CFDP1_LCL", "CFDP1_increase_LCL", "PTPN22_pan", "IL2RA_pan",
           "IFIH1_increase_wb", "SH2B3_increase_wb", "SH2B3_LCL", "SH2B3_increase_LCL", 
           "SMAD3_subadi","MORF4L1_wb", "SH2B3_pan", "MORF4L1_pan", "IKZF4_pan", "SMAD3_sp","IKZF4_LCL", 
           "MORF4L1_LCL",  "TCP11_LCL", "BAK1_pan", "CLEC16A_wb", "IFIH1_pan", "HORMAD2_subadi",  
           "SUOX_pan", "BACH2_pan", "SH2B3_wb", "CTSH_pan","IFIH1_wb", "PTPN22_wb", 
           "TCP11_sp", "SUOX_wb", "SMAD3_increase_sp", "BACH2_increase_pan", "CFDP1_wb",
           "CFDP1_pan", "CFDP1_LCL", "CFDP1_increase_LCL", "BAK1_wb", "CTSH_wb", "CXCL12_pan", 
           "CTSH_LCL", "IL2RA_sp", "CTRB1_pan", "CTRB2_pan", "IKZF4_LCL", "TCP11 sp", "SMAD3")

#main function implementation
for (i in 1:length(genes)){
  for (j in 1:length(outcomes)){
    het_rem(genes[i], outcomes[j])
  }
}

