library(tidyverse)
cols<- c("chrom","pos","ref","alt","qual","vaf","depth","freq")
a <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/COLO829.COLO829_BL/snv_gnomad_af_anno.tsv", col_names = cols)
tp <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL_passed/tp.vcf",col_names = c("chrom","pos","id","ref","alt")) |>
  dplyr::select(chrom, pos, ref, alt) |> add_column(type="TP")

fp <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL_passed/fp.vcf",col_names = c("chrom","pos","id","ref","alt")) |>
  dplyr::select(chrom, pos, ref, alt) |> add_column(type="FP")

data <- a |> left_join(rbind(tp,fp), by=c("chrom","pos","ref","alt"))
data$log_depth <- log(data$depth + 1)
data_scaled <- scale(data[, c("vaf", "qual", "log_depth", "freq")])

pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
pca_data <- as.data.frame(pca_result$x)

pca_data <- pca_data %>% 
  dplyr::select(PC1, PC2,PC3)

# Combine PCA data with the original data (optional)
pca_data <- cbind(pca_data, data)

ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = type), size = 3, alpha=.5) +  # You can color points by any variable
  labs(color="",
       x = "PC1",
       y = "PC2") + scale_color_manual(values=mypal2) + 
  theme_minimal() 

ggplot(pca_data, aes(x = PC1, y = PC3)) +
  geom_point(aes(color = type), size = 3, alpha=.5) +  # Change VAF to any other variable you want to use for coloring
  labs(title = "PCA: PC1 vs PC3",
       x = "Principal Component 1",
       y = "Principal Component 3") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


a
tp <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL_passed/tp.vcf",col_names = c("chrom","pos","id","ref","alt")) |>
  dplyr::select(chrom, pos, ref, alt) |> add_column(type="TP")

fp <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL_passed/fp.vcf",col_names = c("chrom","pos","id","ref","alt")) |>
  dplyr::select(chrom, pos, ref, alt) |> add_column(type="FP")




library(Rtsne)

# Perform t-SNE
duplicates <- duplicated(data_scaled)
data_scaled_unique <- data_scaled[!duplicates, ]
original_data_unique <- data[!duplicates, ]

tsne_result <- Rtsne(data_scaled_unique, dims = 2, perplexity = 30)

# Create a data frame with t-SNE results
tsne_data <- as.data.frame(tsne_result$Y)
colnames(tsne_data) <- c("Dim1", "Dim2")

# Combine t-SNE results with the original data (optional)
tsne_data <- cbind(tsne_data, original_data_unique)

# Plot t-SNE results
ggplot(tsne_data, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = type), size = 3) +
  labs(title = "t-SNE of Variants",
       x = "Dimension 1",
       y = "Dimension 2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
