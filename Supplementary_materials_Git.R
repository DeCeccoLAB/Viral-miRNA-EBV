######### Supplementary materials ###########
# Run this script AFTER the "main article" script

# In this script you will find how we get to the panel in the supplementary materials and the star-Plot in Fig2.

library(edgeR)
library(reshape2)
library(ggplot2)
library(dplyr)

setwd("path/to/Data")

tab <- read.table(file = "EBV_miRNA_table.tsv", header = TRUE, row.names = 1)
m_cluster <- read.csv("metadata.csv", header = TRUE, sep = ";", row.names = 1)
m_cluster_filtered <- m_cluster[m_cluster$Cluster != "Cl4", , drop = FALSE]
rownames(m_cluster_filtered)[rownames(m_cluster_filtered) == "OG60-OE78"] <- "OE78" 
common <- intersect(rownames(m_cluster_filtered), colnames(tab))

m_cluster_filtered <- m_cluster_filtered[common, , drop = FALSE]

tab <- tab[, common, drop = FALSE]
data_batch <- read.table(file = "DATA_batch")
rownames(data_batch)[rownames(data_batch) == "OG60"] <- "OE78"
data_batch <- data_batch[common, ,drop = FALSE]

tab <- ComBat_seq(
  counts = as.matrix(tab),
  batch = data_batch$DATA)

dge <- DGEList(counts = tab)

all_mirna <- rownames(tab)

get_core <- function(x) {
  if (grepl("BHRF1-", x)) {
    sub(".*(BHRF1-[0-9]).*", "\\1", x)
  } else if (grepl("BART", x)) {
    sub(".*(BART[0-9]+).*", "\\1", x)
  } else {
    NA
  }
}

mirna_core <- sapply(all_mirna, get_core)

cluster_BHRF_core  <- c("BHRF1-1", "BHRF1-2", "BHRF1-3")

cluster_1_core <- c("BART1", "BART3", "BART4", "BART5",
                    "BART6", "BART15", "BART16", "BART17")

cluster_2_core <- c("BART7", "BART8", "BART9", "BART10",
                    "BART11", "BART12", "BART13", "BART14",
                    "BART18", "BART19", "BART20", "BART22")

cluster_BART2_core <- c("BART2")  

get_by_core <- function(core_vec) {
  unlist(lapply(core_vec, function(core) {
    hits <- all_mirna[mirna_core == core]
    hits[order(grepl("5p", hits), decreasing = TRUE)]
  }))
}

mirna_order <- c(
  get_by_core(cluster_BHRF_core),
  get_by_core(cluster_1_core),
  get_by_core(cluster_2_core),
  get_by_core(cluster_BART2_core),
  all_mirna[!(all_mirna %in% get_by_core(cluster_BHRF_core) |
                all_mirna %in% get_by_core(cluster_1_core) |
                all_mirna %in% get_by_core(cluster_2_core) |
                all_mirna %in% get_by_core(cluster_BART2_core))]
)


assign_cluster <- function(core) {
  if (is.na(core)) return("Unclustered")
  if (core %in% cluster_BHRF_core)   return("BHRF")
  if (core %in% cluster_1_core)      return("Cluster 1")
  if (core %in% cluster_2_core)      return("Cluster 2")
  if (core %in% cluster_BART2_core)  return("BART2")
  return("Unclustered")
}

plot_global_miRNA_from_dge <- function(dge_object) {
  
  library(edgeR)
  library(reshape2)
  library(ggplot2)
  library(dplyr)
  
  logcpm <- cpm(dge_object, log = TRUE)
  
  df <- melt(logcpm)
  colnames(df) <- c("miRNA", "Sample", "logCPM")
  
  df$miRNA  <- factor(df$miRNA, levels = mirna_order)
  df$core   <- mirna_core[as.character(df$miRNA)]
  df$Cluster <- sapply(df$core, assign_cluster)
  
  
  df$Cluster <- factor(df$Cluster, levels = c("BHRF", "Cluster 1", "Cluster 2", "BART2", "Unclustered"))
  
  cluster_colors <- c(
    "BHRF"       = "#8A2BE2",
    "Cluster 1"  = "#9370DB",
    "Cluster 2"  = "#BA55D3",
    "BART2"      = "#DDA0DD",
    "Unclustered"= "grey70"
  )
  
  ggplot(df, aes(x = miRNA, y = logCPM, fill = Cluster)) +
    geom_boxplot(outlier.size = 0.7) +
    scale_fill_manual(
      values = cluster_colors,
      breaks = c("BHRF", "Cluster 1", "Cluster 2", "BART2", "Unclustered"),
      drop = FALSE
    ) +
    facet_grid(. ~ Cluster, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 14) +
    labs(
      x = "miRNA",
      y = "logCPM",
      title = "Expression profiles of EBV miRNAs (logCPM)"
    ) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title  = element_text(size = 16, face = "bold")
    )
}

plot_global_miRNA_from_dge(dge)

### supp fig 2
tab <- read.table(file = "Tissue_ebv_arm.tsv", header = TRUE, row.names = 1)
m_cluster <- read.csv("m_cluster.csv", header = TRUE, sep = ";", row.names = 1)
m_cluster_filtered <- m_cluster[m_cluster$Cluster != "Cl4", , drop = FALSE]

rownames(m_cluster_filtered)[rownames(m_cluster_filtered) == "OG60-OE78"] <- "OE78" 
common <- intersect(rownames(m_cluster_filtered), colnames(tab))

m_cluster_filtered <- m_cluster_filtered[common, , drop = FALSE]

tab <- tab[, common, drop = FALSE]
data_batch <- read.table(file = "DATA_batch")
rownames(data_batch)[rownames(data_batch) == "OG60"] <- "OE78"
data_batch <- data_batch[common, ,drop = FALSE]

tab <- ComBat_seq(
  counts = as.matrix(tab),
  batch = data_batch$DATA)

metadata_edgeR <- as.data.frame(m_cluster_filtered)

samples_tab_edgeR <- colnames(tab)
group_edgeR <- metadata_edgeR$Cluster[match(samples_tab_edgeR, rownames(metadata_edgeR))]
group_edgeR <- factor(group_edgeR)

dge <- DGEList(counts = tab, group = group_edgeR)
keep <- filterByExpr(dge)
tab_filtrata <- tab[keep, ]
dge <- DGEList(counts = tab_filtrata, group = group_edgeR)
logcpm_before <- cpm(dge, log=TRUE)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log=TRUE)
design <- model.matrix(~ group_edgeR)
outliers <- c("OE92", "LO76", "OE75", "LO99")
tab_complete_out <- tab[, !colnames(tab) %in% outliers]
names(group_edgeR) <- colnames(tab)

group_edgeRO <- group_edgeR[!names(group_edgeR) %in% outliers]


dge1 <- DGEList(counts = tab_complete_out, group = group_edgeRO)
keep1 <- filterByExpr(dge1)
tab_filtrata1 <- tab_complete_out[keep1, ]
dim(tab_filtrata1) 
dge1 <- DGEList(counts = tab_filtrata1, group = group_edgeRO)
logcpm_before1 <- cpm(dge1, log=TRUE)
dge1 <- calcNormFactors(dge1, method = "TMM")
logcpm1 <- cpm(dge1, log=TRUE)
boxplot(logcpm_before1)
boxplot(logcpm1) 
design1 <- model.matrix(~ group_edgeRO)
dge1 <- estimateDisp(dge1, design1) 
fit1 <- glmQLFit(dge1, design1)
qlf.2vs1 <- glmQLFTest(fit1, coef=2) 
deg.2vs1 <- topTags(qlf.2vs1, n=20000, adjust.method = "BH", sort.by = "PValue")$table 
deg.2vs1_rawP <- deg.2vs1[deg.2vs1$PValue < 0.05, ]
up.genes.2vs1_rawP <- row.names(deg.2vs1_rawP[deg.2vs1_rawP$logFC > 0, ])
down.genes.2vs1_rawP <- row.names(deg.2vs1_rawP[deg.2vs1_rawP$logFC < 0, ])


qlf.3vs1 <- glmQLFTest(fit1, coef=3)
deg.3vs1 <- topTags(qlf.3vs1, n=20000, adjust.method = "BH", sort.by = "PValue")$table 
deg.3vs1_rawP <- deg.3vs1[deg.3vs1$PValue < 0.05, ]
up.genes.3vs1_rawP <- row.names(deg.3vs1_rawP[deg.3vs1_rawP$logFC > 0, ])
down.genes.3vs1_rawP <- row.names(deg.3vs1_rawP[deg.3vs1_rawP$logFC < 0, ])

qlf.3vs2 <- glmQLFTest(fit1, contrast=c(0,-1,1)) 
deg.3vs2 <- topTags(qlf.3vs2, n=20000, adjust.method = "BH", sort.by = "PValue")$table 
deg.3vs2_rawP <- deg.3vs2[deg.3vs2$PValue < 0.05, ]
up.genes.3vs2_rawP <- row.names(deg.3vs2_rawP[deg.3vs2_rawP$logFC > 0, ])
down.genes.3vs2_rawP <- row.names(deg.3vs2_rawP[deg.3vs2_rawP$logFC < 0, ])
selected_miRNA <- unique(c(up.genes.3vs2_rawP, down.genes.3vs2_rawP, up.genes.2vs1_rawP, down.genes.2vs1_rawP, down.genes.3vs1_rawP, up.genes.3vs1_rawP))


all_miRNA <- unique(c(
  rownames(deg.2vs1),   
  rownames(deg.3vs1),
  rownames(deg.3vs2)
))


pval_table <- data.frame(
  miRNA = all_miRNA,
  p_2vs1 = NA,
  p_3vs1 = NA,
  p_3vs2 = NA
)

pval_table$p_2vs1[match(rownames(deg.2vs1), pval_table$miRNA)] <- deg.2vs1$PValue
pval_table$p_3vs1[match(rownames(deg.3vs1), pval_table$miRNA)] <- deg.3vs1$PValue
pval_table$p_3vs2[match(rownames(deg.3vs2), pval_table$miRNA)] <- deg.3vs2$PValue

plot_boxplot_all_clusters <- function(logcpm_matrix, group_vector, mirna_list, pval_table) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  
  sub <- logcpm_matrix[mirna_list, , drop = FALSE]
  df <- melt(sub)
  colnames(df) <- c("miRNA", "Sample", "logCPM")
  
  df <- df %>%
    left_join(
      data.frame(Sample = names(group_vector),
                 Cluster = group_vector,
                 stringsAsFactors = FALSE),
      by = "Sample"
    ) %>%
    mutate(Cluster = factor(Cluster, levels = c("Cl1", "Cl2", "Cl3")))
  
  if (any(is.na(df$Cluster))) {
    stop("Errore: alcuni Sample non trovano corrispondenza in group_vector.")
  }
  
  cluster_labels <- c(
    Cl1 = "Immune-active cluster",
    Cl2 = "Defence response cluster",
    Cl3 = "Proliferation cluster"
  )
  
  cluster_colors <- annotation_colors$Cluster
  
  pval_long <- pval_table %>%
    filter(miRNA %in% mirna_list) %>%
    pivot_longer(cols = starts_with("p_"),
                 names_to = "comparison",
                 values_to = "p") %>%
    mutate(
      group1 = case_when(
        comparison == "p_2vs1" ~ "Cl1",
        comparison == "p_3vs1" ~ "Cl1",
        comparison == "p_3vs2" ~ "Cl2"
      ),
      group2 = case_when(
        comparison == "p_2vs1" ~ "Cl2",
        comparison == "p_3vs1" ~ "Cl3",
        comparison == "p_3vs2" ~ "Cl3"
      )
    )
  
  # 5) Calcolo y.position per ogni miRNA
  ymax_df <- df %>%
    group_by(miRNA) %>%
    summarise(ymax = max(logCPM))
  
  pval_long <- pval_long %>%
    left_join(ymax_df, by = "miRNA") %>%
    mutate(
      p = signif(p, 3),
      y.position = case_when(
        comparison == "p_2vs1" ~ ymax + 0.5,
        comparison == "p_3vs1" ~ ymax + 1.0,
        comparison == "p_3vs2" ~ ymax + 1.5
      ),
      x = NA,
      fill = NA
    )
  
  
  # 6) Plot finale
  ggplot(df, aes(x = Cluster, y = logCPM, fill = Cluster)) +
    geom_boxplot(outlier.size = 0.8, notch = TRUE) +
    facet_wrap(~ miRNA, scales = "free_y") +
    scale_fill_manual(values = cluster_colors, labels = cluster_labels) +
    scale_x_discrete(labels = cluster_labels) +
    stat_pvalue_manual(
      pval_long,
      label = "p",
      tip.length = 0.01,
      size = 3,
      inherit.aes = FALSE
    ) +
    theme_bw(base_size = 14) +
    labs(x = "Cluster", y = "logCPM") +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
}


plot_boxplot_all_clusters(
  logcpm_matrix = logcpm1,
  group_vector = group_edgeRO,
  mirna_list = selected_miRNA,
  pval_table = pval_table
)


### starplot
plot_STAR_miRNA <- function(logcpm_matrix, group_vector, mirna_list, cluster_labels, cluster_colors) {
  
  library(ggplot2)
  library(dplyr)
  
  # 1) miRNA subset
  sub <- logcpm_matrix[mirna_list, , drop = FALSE]
  
  # 2) PCA
  pca <- prcomp(t(sub), scale. = TRUE)
  
  df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Cluster = group_vector[rownames(pca$x)]
  )
  
  # 3) cluster name
  df$Cluster_label <- factor(cluster_labels[df$Cluster],
                             levels = cluster_labels)
  
  # 4) Centroids calculation
  centroids <- df %>%
    group_by(Cluster_label) %>%
    summarise(
      PC1 = mean(PC1),
      PC2 = mean(PC2)
    )
  
  # 5) Merge
  df <- df %>%
    left_join(centroids, by = "Cluster_label", suffix = c("", "_centroid"))
  
  # 6) STAR PLOT
  ggplot(df) +
    geom_segment(aes(x = PC1_centroid, y = PC2_centroid,
                     xend = PC1, yend = PC2,
                     color = Cluster_label),
                 alpha = 0.5, linewidth = 0.8) +
    geom_point(aes(x = PC1, y = PC2, color = Cluster_label),
               size = 3, alpha = 0.9) +
    geom_point(data = centroids,
               aes(x = PC1, y = PC2, color = Cluster_label),
               size = 6, shape = 21, fill = "white", stroke = 1.8) +
    geom_text(data = centroids,
              aes(x = PC1, y = PC2, label = Cluster_label),
              vjust = -1, fontface = "bold", size = 4) +
    theme_bw(base_size = 14) +
    scale_color_manual(values = setNames(cluster_colors, cluster_labels)) +
    labs(
      x = "PC1",
      y = "PC2",
      color = "Cluster",
      title = "Star plot of selected miRNAs (logCPM)"
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
}

plot_STAR_miRNA(
  logcpm_matrix = logcpm1,
  group_vector = group_edgeRO,
  mirna_list = selected_miRNA,
  cluster_labels = cluster_labels,
  cluster_colors = cluster_colors
)

