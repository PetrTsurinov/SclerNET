#!/usr/bin/env Rscript
library(limma)
library(annotate)
library(stringr)
library(biomaRt)
library(data.table)
library(dplyr)
library(filesstrings)

args = commandArgs(trailingOnly=TRUE)
print(args)

work_dir <- args[1]
output_path <- args[2]
chip_type <- args[3]

#define type of microarray platform
read_raw_data <- function(work_dir, type = c('affy', 'oligo')) {
  setwd(work_dir)
  
  control_files <- as.data.frame(dir(path='control', pattern = '.CEL'))
  colnames(control_files) <- 'filename'
  control_files$filepath <- paste0(work_dir, "/control/", control_files$filename)
  control_files$type <- 'control'
  
  admission_files <- as.data.frame(dir(path='admission', pattern = '.CEL'))
  colnames(admission_files) <- 'filename'
  admission_files$filepath <- paste0(work_dir, "/admission/", admission_files$filename)
  admission_files$type <- 'admission'
  
  all_files <- rbind(control_files, admission_files)
  for (file in all_files$filepath) {
    file.move(file, work_dir)
  }
  
  if (type == 'affy') {
    library(affy)
   
    write.table(all_files, 'targets.txt', row.names = F, col.names = T, quote = F, sep="\t")
    targets <- readTargets('targets.txt', sep="\t", row.names="filename")
    raw_data <- ReadAffy(filenames = targets$fileName)
  }
  
  if (type == 'oligo') {
    library(oligo)
    
    celFiles <- list.celfiles(work_dir, full.names=TRUE)
    raw_data <- read.celfiles(celFiles)
  }
  
  design <- cbind(control = 1, level = all_files$type == "admission")
  return(c(raw_data, as.data.frame(design)))
}

array_process <- function(raw_data, design, output_path){
  #create and normalise expression set object
  eset <- rma(raw_data)
  
  fit <- lmFit(eset, design) ##create linear model and get DEGs via design matrix
  fit <- eBayes(fit, trend = F, robust = F) ##normalise
  all_BH <- topTable(fit, coef=2, number = Inf, adjust.method = 'BH') ##create normal dataframe, p.adj default
  filtered_BH <- all_BH[all_BH$adj.P.Val < 0.05,]
  
  print('Amount of gene, p-adj < 0.05:')
  print(nrow(filtered_BH))
  
  ##create annotation (AGI code) to Affy ID
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  #symbols_table <- getBM(attributes = c('affy_hugene_1_0_st_v1', 'hgnc_symbol', 'chromosome_name',
  #                                      'start_position', 'end_position', 'band'),
  #                       filters = 'affy_hugene_1_0_st_v1', 
  #                       values = rownames(filtered_BH), 
  #                       mart = ensembl)
  load('../data/symbols_table.RData')
  tab <- table(symbols_table$affy_hugene_1_0_st_v1)
  
  #find unique ID, multi-ID remove 
  unique_ids <- names(tab)[tab == 1]
  unique_ids_symbols_table <- symbols_table[symbols_table$affy_hugene_1_0_st_v1 %in% unique_ids, ]
  filtered_BH <- setDT(filtered_BH, keep.rownames = TRUE)[]
  colnames(filtered_BH)[1] <- "affy_id"
  filtered_BH$affy_id <- as.numeric(as.character(filtered_BH$affy_id))
  joined_BH <- merge(filtered_BH, unique_ids_symbols_table, by.x=c("affy_id"), by.y=c("affy_hugene_1_0_st_v1"))
  joined_BH <- joined_BH[joined_BH$hgnc_symbol != ""]
  filtered_joined_BH <- joined_BH %>%
    group_by(hgnc_symbol) %>%
    filter(AveExpr == max(AveExpr))
  filtered_joined_BH$hgnc_symbol <- toupper(filtered_joined_BH$hgnc_symbol)
  df_filt <- filtered_joined_BH[(filtered_joined_BH$logFC > log2(1.5) | filtered_joined_BH$logFC < -log2(1.5)), ]
  
  write.csv(filtered_joined_BH, output_path, row.names = F, quote = F)
  write.csv(df_filt, gsub("\\.csv", "_filtered\\.csv", output_path), row.names = F, quote = F)
  return(filtered_joined_BH)
}

raw_data <- read_raw_data(work_dir, chip_type)
result_df <- array_process(raw_data[[1]], 
                           design = as.data.frame(cbind(raw_data[[2]], raw_data[[3]])), 
                           output_path = output_path)