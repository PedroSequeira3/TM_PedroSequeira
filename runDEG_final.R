library(NOISeq)
library(edgeR)
library(DESeq)
library(dplyr)
library(tidyverse)
library(tibble)
library(VennDiagram)
library(ggplot2)

# Function to read and process STAR --quantMode output files
STAR_read_and_preprocess <- function(filepath, count_col) {
  data <- read.table(filepath, header = FALSE, sep = "\t", skip = 3, col.names = c("Gene", "U", "F", "R"))
  data <- data[, c("Gene", count_col)]
  return(data)
}

#Select STAR --quantMode alignment files (.out.tab) 
ControlSample_Replicate1 <-  #Replicate 1 for control conditions
ControlSample_Replicate2 <-  #Replicate 2 for control conditions
Sample_Replicate1 <-         #Replicate 1 for stress conditions
Sample_Replicate2 <-         #Replicate 2 for stress conditions
  
  
#Data preparation
Sample1_1 <- STAR_read_and_preprocess(ControlSample_Replicate1, "U")
Sample1_2 <- STAR_read_and_preprocess(ControlSample_Replicate2, "U")
Sample2_1 <- STAR_read_and_preprocess(Sample_Replicate1, "U")
Sample2_2 <- STAR_read_and_preprocess(Sample_Replicate2, "U")

Sample1 <- merge(Sample1_1, Sample1_2, by = "Gene", suffixes = c(".CELTA1", ".CELTA2"))
Sample2 <- merge(Sample2_1, Sample2_2, by = "Gene", suffixes = c(".TE1", ".TE2"))
Data <- merge(Sample1, Sample2, by = "Gene")

rownames(Data) <- Data$Gene
Data <- Data[-1, ]
Data_HC <- Data %>% #Remove header
  filter(!grepl("LC", Gene))

Data_HC <- Data_HC[, -1]

colnames(Data_HC) <- c("CELTA1", "CELTA2", "TE1", "TE2")

#Number of genes were expressed in the samples
genes_expressed <- sum(rowSums(Data_HC) > 0)

#DESeq data preparation
Sample1_Sum <- Data_HC$CELTA1 + Data_HC$CELTA2
Sample2_Sum <- Data_HC$TE1 + Data_HC$TE2
Genes <- rownames_to_column(Data_HC, var="rowname")
Genes <- Genes[,1]
Data_HC2 <- data.frame(Genes, Sample1 = Sample1_Sum, Sample2 = Sample2_Sum)

rownames(Data_HC2) <- Data_HC2$Gene
Data_HC2 <- Data_HC2[, -1]

#NOISeq DEG
myfactors = data.frame(Group = c("XCELTA", "XCELTA", "TE", "TE"))

Data_NOI_HC <- readData(data = Data_HC, factors = myfactors)
myTMM_HC = tmm(assayData(Data_NOI_HC)$exprs, long = 1000, lc = 0)

mynoiseq_HC = noiseq(Data_NOI_HC, k = 0.5, norm = "tmm", factor = "Group", pnr = 0.2, nss = 5, v = 0.01, lc = 1, replicates = "technical")
mynoiseq_results_HC <- mynoiseq_HC@results[[1]]

upregulated_NOISeq_HC = degenes(mynoiseq_HC, q = 0.95, M = "up") 
downregulated_NOISeq_HC = degenes(mynoiseq_HC, q = 0.95, M = "down")
de_NOISeq_HC = degenes(mynoiseq_HC, q = 0.95)

mynoiseq_de_HC <- mynoiseq_results_HC[mynoiseq_results_HC[, 6] != 0, ]

mynoiseq_de_genes_HC <- rownames(mynoiseq_de_HC)

de_genes_NOISeq_HC <- rbind(downregulated_NOISeq_HC, upregulated_NOISeq_HC)

#EdgeR DEG
Group <- factor(c("CELTA", "CELTA", "TE", "TE"))
edgeR_HC <- DGEList(counts = Data_HC, genes = rownames(Data_HC), group = Group)

keep <- filterByExpr(edgeR_HC)
edgeR_HC <- edgeR_HC[keep, , keep.lib.sizes = FALSE]

edgeR_HC <- calcNormFactors(edgeR_HC, method = "TMM")

edgeR_HC <- estimateDisp(edgeR_HC)

bcv <- 0.01 
edgeR_et_HC <- exactTest(edgeR_HC, dispersion = bcv^2)

genes_edgeR_HC <- topTags(edgeR_et_HC, n=Inf, p.value = 0.05)$table

upregulated_edgeR_HC <- genes_edgeR_HC[genes_edgeR_HC$logFC > 2, ]
downregulated_edgeR_HC <- genes_edgeR_HC[genes_edgeR_HC$logFC < -2, ]

de_genes_edgeR_HC <- rbind(downregulated_edgeR_HC, upregulated_edgeR_HC)

#DESeq DEG
Group <- factor(c("CELTA", "TE"))
cds_HC <- newCountDataSet( Data_HC2, Group )

cds_HC <- estimateSizeFactors( cds_HC )
size_factors_HC <- sizeFactors(cds_HC)

cds_HC <- estimateDispersions(cds_HC, method = "blind", sharingMode="fit-only")
res_HC <- nbinomTest( cds_HC, "CELTA", "TE" )
res4_HC <- res_HC[is.na(res_HC$log2FoldChange), ] #Remove NaN log2FoldChange
res5_HC <- res_HC[is.finite(res_HC$log2FoldChange), ] #Remove NaN & Inf log2FoldChange

DESeq_de_genes_HC <- res5_HC[, c(1,6,8)]

downregulated_DESeq_HC <- DESeq_de_genes_HC[DESeq_de_genes_HC$log2FoldChange < -2 & DESeq_de_genes_HC$padj < 0.05, 1:3]
upregulated_DESeq_HC <- DESeq_de_genes_HC[DESeq_de_genes_HC$log2FoldChange > 2 & DESeq_de_genes_HC$padj < 0.05, 1:3]

de_genes_DESeq_HC <- rbind(downregulated_DESeq_HC, upregulated_DESeq_HC)

#Intersection of DEGs
common_upregulated_genes_HC <- Reduce(intersect, list(rownames(upregulated_edgeR_HC), rownames(upregulated_NOISeq_HC), upregulated_DESeq_HC$id))
common_downregulated_genes_HC <- Reduce(intersect, list(rownames(downregulated_edgeR_HC), rownames(downregulated_NOISeq_HC), downregulated_DESeq_HC$id))
common_de_genes_HC <- Reduce(intersect, list(rownames(de_genes_edgeR_HC), rownames(de_genes_NOISeq_HC), de_genes_DESeq_HC$id))

#Functional Annotation
file_path <- "~/Desktop/Tese/Genome/iwgsc_refseqv2.1_functional_annotation.csv"
annotation_data <- read.csv(file_path, header = FALSE, stringsAsFactors = FALSE)
colnames(annotation_data) <- c("Element", "Source", "Category", "Identifier", "Description")

annotation_list <- list()
for (element in common_de_genes_HC) {
  filtered_data <- annotation_data[annotation_data$Element == element, ]
  human_readable <- filtered_data$Source == "Human readable description"
  filtered_data[human_readable, c("Identifier", "Description")] <- filtered_data[human_readable, c("Description", "Identifier")]
  filtered_data[human_readable, "Identifier"] <- ""
  annotation_list[[element]] <- filtered_data
}
annotation_results <- do.call(rbind, annotation_list)
row.names(annotation_results) <- NULL

#Retrieve the lines that correspond to Gene Ontology Source
GO_annotations <- annotation_results[annotation_results[, 2] == "Gene Ontology", ]
clean_GO <- GO_annotations[,c(1,3,4,5)]
colnames(clean_GO) <- c("Gene", "Category","GO_ID", "GO_Description")

#Retrieve the lines that correspond to Blast-Hit-Accession Source
BLAST_annotations <- annotation_results[annotation_results[, 2] == "Blast-Hit-Accession", ]
clean_BLAST <- BLAST_annotations[,c(1,4)]
colnames(clean_BLAST) <- c("Gene", "Accession_ID")

#Retrieve the lines that correspond to Human readable description Source
HUMAN_annotations <- annotation_results[annotation_results[, 2] == "Human readable description", ]
clean_HUMAN <- HUMAN_annotations[,c(1,5)]
colnames(clean_HUMAN) <- c("Gene", "Human_Description")

combined_annotations <- rbind(GO_annotations, BLAST_annotations, HUMAN_annotations)
combined_annotations <- combined_annotations[order(combined_annotations[[1]]), ]

#Combine multiples blast-hit-accession entries for the same gene
clean_BLAST2 <- clean_BLAST %>%
  group_by(Gene) %>%
  summarise(Accession_ID = paste(Accession_ID, collapse = ", "), .groups = "drop")

#Merge the GO entries with the blast entries
GO_BLAST <- clean_GO %>%
  left_join(clean_BLAST2, by = "Gene")

#Combine multiples human readable description entries for the same gene
clean_HUMAN2 <- clean_HUMAN %>%
  group_by(Gene) %>%
  summarise(Human_Description = paste(Human_Description, collapse = ", "), .groups = "drop")

#Merge the GO+blast entries with the human entries
GO_BLAST_HUMAN <- GO_BLAST %>%
  left_join(clean_HUMAN2, by = "Gene")
GO_BLAST_HUMAN <- GO_BLAST_HUMAN[order(GO_BLAST_HUMAN[[1]]), ]

edgeR_logFC <- data_frame(rownames(de_genes_edgeR_HC), de_genes_edgeR_HC$logFC)
colnames(edgeR_logFC) <- c("Gene", "edgeRlogFC")
NOISeq_logFC <- data_frame(rownames(de_genes_NOISeq_HC), de_genes_NOISeq_HC$M)
colnames(NOISeq_logFC) <- c("Gene", "NOISeqlogFC")
DESeq_logFC <- data.frame(de_genes_DESeq_HC$id, de_genes_DESeq_HC$log2FoldChange)
colnames(DESeq_logFC) <- c("Gene", "DESeqlogFC")

all_logFC <- merge(edgeR_logFC, NOISeq_logFC, by = "Gene")
all_logFC <- merge(all_logFC, DESeq_logFC, by = "Gene")

#Check if for the same gene all logFC are positive/negative
logFC_check <- all(apply(all_logFC[, 2:4], 1, function(row) all(row > 0) | all(row < 0)))

GO_BLAST_HUMAN_logFC <- merge(GO_BLAST_HUMAN, all_logFC, by = "Gene")
GO_BLAST_HUMAN_logFC <- GO_BLAST_HUMAN_logFC[, c(1,6, 3, 2, 4, 5, 7, 8, 9)]

BP_ID_DESC <- GO_BLAST_HUMAN_logFC[GO_BLAST_HUMAN_logFC[,4] == "BP", c(1,3,5)]
colnames(BP_ID_DESC) <- c("Gene_ID","GO_ID", "Description")

MF_ID_DESC <- GO_BLAST_HUMAN_logFC[GO_BLAST_HUMAN_logFC[,4] == "MF", c(1,3,5)]
colnames(MF_ID_DESC) <- c("Gene_ID","GO_ID", "Description")

CC_ID_DESC <- GO_BLAST_HUMAN_logFC[GO_BLAST_HUMAN_logFC[,4] == "CC", c(1,3,5)]
colnames(CC_ID_DESC) <- c("Gene_ID","GO_ID", "Description")

FINAL <- GO_BLAST_HUMAN_logFC[,c(1,7,8,9,2,6)]
colnames(FINAL) <- c("Gene_ID", "edgeR_logFC","NOISeq_logFC", "DESeq_logFC", "Description", "Accession_ID")


#Export csv files for Biological Analysis
write.csv(FINAL, file = "~/Desktop/TE_CSvsTE_FS_Genes.csv", row.names = FALSE)
write.csv(BP_ID_DESC, file = "~/Desktop/TE_CSvsTE_FS_BiologicalProcess.csv", row.names = FALSE)
write.csv(MF_ID_DESC, file = "~/Desktop/TE_CSvsTE_FS_MolecularFunction.csv", row.names = FALSE)
write.csv(CC_ID_DESC, file = "~/Desktop/TE_CSvsTE_FS_CelularComponent.csv", row.names = FALSE)


#Statistic check
Stats <- data.frame(nrow(Data), nrow(Data_HC), sum(rowSums(Data_HC) > 0),logFC_check, nrow(de_genes_edgeR_HC), nrow(de_genes_NOISeq_HC), nrow(de_genes_DESeq_HC), length(common_de_genes_HC), sum(rowSums(abs(all_logFC[, 2:4]) < 3) == 3), sum(rowSums(abs(all_logFC[, 2:4]) < 4) == 3)   )
colnames(Stats) <- c("All genes","All non LC genes","Genes expressed","logFC check","edgeR DE genes","NOISeq DE genes","DESeq DE genes","Intersection DE genes","|logFC| < 3"," |logFC| < 4")
#View(Stats)

#Plots
values <- FINAL[[2]] 
df <- data.frame(values = values)

#Density Plot
ggplot(data.frame(values), aes(x = values)) +
  geom_density(data = subset(df, values < 0), fill = "red", alpha = 0.5) + 
  geom_density(data = subset(df, values >= 0), fill = "green", alpha = 0.5) +
  scale_x_continuous(breaks = c(seq(floor(min(values)), -2, by = 1), 2, seq(3, ceiling(max(values)), by = 1))) +
  theme_minimal() +
  labs(title = "log2FC Density", x = "Values", y = "Density")
  
print(((sum(values >= 2 & values <= 3) / length(values)) * 100)+((sum(values >= -3 & values <= -2) / length(values)) * 100))

#Volcano Plot using edgeR pvalue 
volcano_data <- FINAL %>%
  distinct(Gene_ID, .keep_all = TRUE) %>%  
  inner_join(de_genes_edgeR_HC[, c("genes", "PValue")], by = c("Gene_ID" = "genes")) %>% 
  mutate(
    neg_log10_pval = -log10(ifelse(PValue == 0, 1e-323, PValue)),  
    Color = ifelse(edgeR_logFC > 0, "Positive", "Negative")  
  ) %>%
  select(Gene_ID, edgeR_logFC, neg_log10_pval, Color)

ggplot(volcano_data, aes(x = edgeR_logFC, y = neg_log10_pval, color = Color)) +
  geom_point(size = 0.8, alpha = 0.6) +  
  theme_minimal() +
  scale_color_manual(values = c("Negative" = "red", "Positive" = "green")) +  
  labs(title = "logFC p-value Volcano Plot", x = "log Fold Change (logFC)", y = "-log10(p-value)") + 
  theme(legend.position = "none") +
  ylim(0, 325)

#Intersection Plots
grid.newpage()
venn.plot <- venn.diagram(
  x = list(edgeR = de_genes_edgeR_HC[[1]], DESeq = de_genes_DESeq_HC[[1]], NOISeq = rownames(de_genes_NOISeq_HC)),
  fill = c("blue", "red", "green"),
  col = "black", # Black border lines
  filename = NULL,
  lty = "solid", # Solid line type
  lwd = 2, # Line width
  alpha = 0, # Fully transparent (no fill)
  cex = 1.5, # Font size for counts
  cat.cex = 1.2, # Font size for category names
  cat.pos = c(0,0,180) # Automatic category placement
)
grid::grid.draw(venn.plot)

grid.newpage()
venn.plot <- venn.diagram(
  x = list(
    edgeR = rownames(downregulated_edgeR_HC),   # full set, not [1]
    DESeq = downregulated_DESeq_HC[,1],   # full set
    NOISeq = rownames(downregulated_NOISeq_HC)  # full set
  ),
  fill = c("blue", "red", "green"), 
  col = "black",
  filename = NULL,
  euler.d = FALSE,
  scaled = FALSE,
  lty = "solid",
  lwd = 2,
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(0,0,180)
)

grid::grid.draw(venn.plot)

grid.newpage()
venn.plot <- venn.diagram(
  x = list(
    edgeR = rownames(upregulated_edgeR_HC),   # full set, not [1]
    DESeq = upregulated_DESeq_HC[,1],   # full set
    NOISeq = rownames(upregulated_NOISeq_HC)  # full set
  ),
  fill = c("blue", "red", "green"), 
  col = "black",
  filename = NULL,
  euler.d = FALSE,
  scaled = FALSE,
  lty = "solid",
  lwd = 2,
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(0,0,180)
)

grid::grid.draw(venn.plot)