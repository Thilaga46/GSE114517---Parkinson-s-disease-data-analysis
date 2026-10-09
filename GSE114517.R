#script for differential gene expression analysis (GSE114517)

setwd("~/GSE114517 Data analysis/GSE114517---Parkinson-s-disease-data-analysis")

#load libraries
library(dplyr) 
library(tidyverse) 
library(GEOquery)

#select .txt files in the directory
file_list <- list.files(pattern = "\\.txt$")

#combine target files into a single file using the common denominator (here GeneID is common to all files)
combined_counts <- file_list %>%
  map(~ read_tsv(.x, col_names = c("GeneID", .x), show_col_types = FALSE)) %>%
  reduce(full_join, by = "GeneID")

#filter and replace the missing values with 0
combined_counts <- combined_counts %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  filter(!str_starts(GeneID, "__"))

write_csv(combined_counts, "GSE114517_SN_combined_counts.csv")
cat("Success! Real matrix dimensions are now:", dim(combined_counts), "\n")

#head(combined_counts)

options(download.file.method = "libcurl")

#get metadata in tabular format instead of downloading from NCBI and avoid multiple steps of manipulation
gse <- getGEO(GEO = "GSE114517", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)

geo_metadata <- pData(phenoData(gse[[1]]))

#select desired columns in the metadata and rename them for further analysis, which is to join with the manipulated data file
geo_metadata_clean <- geo_metadata %>%
  select(geo_accession, characteristics_ch1) %>%
  dplyr::rename(condition = characteristics_ch1) %>%
  mutate(condition = gsub("subject status: ", "", condition)) %>% 
  mutate(condition = case_when(
    str_detect(condition, "PD") ~ "PD",
    str_detect(condition, "Control") ~ "Control",
    TRUE ~ condition
  ))

#head(geo_metadata_clean)

#extract only the sample names from the merged data matrix, while intentionally leaving out the "GeneID" column
sample_cols <- colnames(combined_counts)[-1]
sample_ids <- str_extract(sample_cols, "GSM[0-9]+")

#head(sample_ids)

metadata <- tibble(
  sampleID = sample_cols,
  geo_accession = sample_ids) %>%
  left_join(geo_metadata_clean, by = "geo_accession")

#head(metadata)

#quality control check -> counting how many missing values exist in the condition column of the metadata table 
#(0 missing values indicate successful matching of the samples "control" or "condition" with the metadata)
sum(is.na(metadata$condition))

#reshape metadata table to fit the strict formatting rules required by the Bioconductor DESeq2 package
metadata$Condition <- factor(metadata$condition, levels = c("Control", "PD"))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sampleID

all(rownames(metadata) == colnames(combined_counts)[-1])

library(DESeq2)

#create matrix
count_matrix <- combined_counts %>%
  column_to_rownames("GeneID") %>%
  as.matrix() %>%
  round()

#feed the pre-made matrix directly into DESeq2 
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ Condition)

#calculate normalization factors
dds <- estimateSizeFactors(dds)

#extract the normalized counts matrix
normalized_counts <- counts(dds, normalized = TRUE)

normalized_df <- as.data.frame(normalized_counts) %>%
  rownames_to_column("GeneID")

write_csv(normalized_df, "GSE114517_SN_normalized_counts.csv")
cat("Success! Saved normalized matrix with dimensions:", dim(normalized_df), "\n")

#metadata file with gene list
mapping_file <- read_csv("GSE114517_Diff_PDvsCont_OG.csv")

#select rows with gene list in the metadata
gene_map <- mapping_file %>%
  select(Row.names, gene_name) %>%
  distinct() 

#map metadata file with gene list to the manipulated and normalized dataframes
normalized_with_gene_names <- normalized_df %>%
  left_join(gene_map, by = c("GeneID" = "Row.names")) %>%
  select(GeneID, gene_name, everything())

write_csv(normalized_with_gene_names, "GSE114517_SN_normalized_with_genenames.csv")
cat("Mapping complete! New dimensions:", dim(normalized_with_gene_names), "\n")

#reshaping data 
dat.long <- normalized_with_gene_names %>%
  select(-GeneID) %>%
  pivot_longer(cols = -gene_name, names_to = "sampleID", values_to = "Normalized")

#head(dat.long)

metadata_slim <- metadata %>%
  select(sampleID, Condition)

#join dataframes
dat_joined <- dat.long %>%
  left_join(metadata_slim, by = "sampleID")

sum(is.na(dat_joined$Condition))

#select, filter, and group genes of interest
genes_of_interest <- c('CYP46A1', 'KEAP1', 'NFE2L2', 'CUL3', 'GSK3B', 'SOD1', 'HO1', 'NQO1', 'HMOX1', 'GCLC', 'SNCA', 'TH', 'ESR1')

dat_filtered <- dat_joined %>%
  filter(gene_name %in% genes_of_interest)

dat_filtered %>%
  filter(gene_name %in% c('CYP46A1', 'KEAP1', 'NFE2L2', 'CUL3', 'GSK3B', 'SOD1', 'HO1', 'NQO1', 'HMOX1', 'GCLC', 'SNCA', 'TH', 'ESR1')) %>%
  group_by(gene_name, Condition) %>%
  summarize(mean_Normalized = mean(Normalized),
            median_Normalized = median(Normalized),
            .groups = 'drop')

library(ggplot2)

#data visualization by plotting using Heatmap
#pdf("heatmap_GSE114517.pdf", width = 10, height = 8) --> save the results in pdf format
dat_filtered %>%
  filter(gene_name %in% c('CYP46A1', 'KEAP1', 'NFE2L2', 'CUL3', 'GSK3B', 'SOD1', 'HO1', 'NQO1', 'HMOX1', 'GCLC', 'SNCA', 'TH', 'ESR1')) %>%
  ggplot(., aes(x = sampleID, y = gene_name, fill = Normalized)) + #check the labels properly from the data to plot them correctly
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'darkgreen') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~ Condition, scales = "free_x", space = "free_x") #to group "control" and "PD" for better visualization



