library(edgeR) #4.4.1
library(sva) #3.54.0
library(pheatmap) #1.0.12
library(org.Hs.eg.db) #3.20.0
library(clusterProfiler) #4.14.6
library(patchwork) #1.3.0
library(ggplot2) #3.5.1
library(readxl) #1.4.3
library(ggpubr) #0.6.0.999
library(igraph) #2.1.2
library(ggraph) #2.2.2
library(ReactomePA) #1.50.0
library(fgsea) #1.32.2
library(immunedeconv) #2.1.4
library(stringr) #1.5.1
library(patchwork) #1.3.0
library(tidyr) #1.3.1
library(dplyr) #2.5.0
library(immunedeconv) #2.1.4 !!! this package can be very hard to download after 2024, if you are NOT interested in the ESTIMATE part, just skip it, otherwise I would suggest to download it using the package "remote" or using the latest version of Rtools

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

metadata_edgeR <- as.data.frame(m_cluster_filtered)

samples_tab_edgeR <- colnames(tab)
group_edgeR <- metadata_edgeR$Cluster[match(samples_tab_edgeR, rownames(metadata_edgeR))]
group_edgeR <- factor(group_edgeR)

dge <- DGEList(counts = tab, group = group_edgeR)
keep <- filterByExpr(dge)
tab_filtrata <- tab[keep, ]

dim(tab_filtrata)
dge <- DGEList(counts = tab_filtrata, group = group_edgeR)
logcpm_before <- cpm(dge, log=TRUE)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log=TRUE)
design <- model.matrix(~ group_edgeR)
design
plotMDS(logcpm)

outliers <- c("OE92", "LO76", "OE75", "LO99")
tab_complete_out <- tab[, !colnames(tab) %in% outliers]
names(group_edgeR) <- colnames(tab)

group_edgeRO <- group_edgeR[!names(group_edgeR) %in% outliers]
table(group_edgeRO)  #20 33 34
group_edgeRCl1 <- factor(ifelse(group_edgeRO == "Cl1", "Cl1", "Other"))
table(group_edgeRCl1)

dge1 <- DGEList(counts = tab_complete_out, group = group_edgeRCl1)
keep1 <- filterByExpr(dge1)
tab_filtrata1 <- tab_complete_out[keep1, ]
dim(tab_filtrata1) #37 87
dge1 <- DGEList(counts = tab_filtrata1, group = group_edgeRCl1)
logcpm_before1 <- cpm(dge1, log=TRUE)
dge1 <- calcNormFactors(dge1, method = "TMM")
logcpm1 <- cpm(dge1, log=TRUE)
design1 <- model.matrix(~ group_edgeRCl1)
plotMDS(logcpm1, labels = group_edgeRCl1)

dge1 <- estimateDisp(dge1, design1) 
fit1 <- glmQLFit(dge1, design1) 
qlf.Cl1 <- glmQLFTest(fit1, coef=2) 
deg.Cl1 <- topTags(qlf.Cl1, n=20000, adjust.method = "BH", sort.by = "PValue")$table 
deg.Cl1_rawP <- deg.Cl1[deg.Cl1$PValue < 0.05, ]
up.genes.Cl1_rawP   <- rownames(deg.Cl1_rawP[deg.Cl1_rawP$logFC < 0, ])
down.genes.Cl1_rawP <- rownames(deg.Cl1_rawP[deg.Cl1_rawP$logFC > 0, ])
print(c(up.genes.Cl1_rawP, deg.Cl1$PValue[deg.Cl1$FDR < 0.05]))
print(c(up.genes.Cl1_rawP, deg.Cl1$PValue[rownames(deg.Cl1) == "ebv-miR-BART17-5p"]))


#
miRNACl1 <- "ebv-miR-BART17-5p"

df <- data.frame(
  expr = as.numeric(logcpm1[miRNACl1, ]),
  group = factor(group_edgeRCl1, levels = c("Cl1", "Other"))
)
df$group <- factor(df$group,
                   levels = c("Cl1", "Other"),
                   labels = c("Immune-active", "Other"))


ggboxplot(df, x = "group", y = "expr",
          fill = "group",
          palette = c("Immune-active" = "indianred1",
                      "Other" = "cornsilk1")) +
  
  annotate("text",
           x = 1.5,   
           y = max(df$expr) * 1.1,
           label = "FDR = 0.029",
           size = 5,
           fontface = "bold") +
  
  labs(
    title = "Expression of ebv-miR-BART17-5p",
    x = "",
    y = "logCPM"
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )


group_edgeRCl3 <- factor(ifelse(group_edgeRO == "Cl3", "Cl3", "Other"))
table(group_edgeRCl3)

dge1 <- DGEList(counts = tab_complete_out, group = group_edgeRCl3)
keep1 <- filterByExpr(dge1)
tab_filtrata1 <- tab_complete_out[keep1, ]
dim(tab_filtrata1) #37 87
dge1 <- DGEList(counts = tab_filtrata1, group = group_edgeRCl3)
logcpm_before1 <- cpm(dge1, log=TRUE)
dge1 <- calcNormFactors(dge1, method = "TMM")
logcpm1 <- cpm(dge1, log=TRUE)
design1 <- model.matrix(~ group_edgeRCl3)
plotMDS(logcpm1, labels = group_edgeRCl3)

dge1 <- estimateDisp(dge1, design1)
fit1 <- glmQLFit(dge1, design1) 
qlf.Cl3 <- glmQLFTest(fit1, coef=2) 
deg.Cl3 <- topTags(qlf.Cl3, n=20000, adjust.method = "BH", sort.by = "PValue")$table 
deg.Cl3_rawP <- deg.Cl3[deg.Cl3$PValue < 0.05, ]
up.genes.Cl3_rawP   <- rownames(deg.Cl3_rawP[deg.Cl3_rawP$logFC < 0, ])
down.genes.Cl3_rawP <- rownames(deg.Cl3_rawP[deg.Cl3_rawP$logFC > 0, ])
print(c(up.genes.Cl3_rawP, deg.Cl3$PValue[rownames(deg.Cl3) == "ebv-miR-BART3-5p"], deg.Cl3$PValue[rownames(deg.Cl3) == "ebv-miR-BART15"]))


miRNACl1 <- "ebv-miR-BART3-5p"

df <- data.frame(
  expr = as.numeric(logcpm1[miRNACl1, ]),
  group = factor(group_edgeRCl3, levels = c("Cl3", "Other"))
)
df$group <- factor(df$group,
                   levels = c("Cl3", "Other"),
                   labels = c("Proliferation", "Other"))


ggboxplot(df, x = "group", y = "expr",
          fill = "group",
          palette = c("Proliferation" = "palegreen1",
                      "Other" = "cornsilk1")) +
  
  
  annotate("text",
           x = 1.5,   
           y = max(df$expr) * 1.1,
           label = "FDR = 0.0098",
           size = 5,
           fontface = "bold") +
  
  labs(
    title = "Expression of ebv-miR-BART3-5p",
    x = "",
    y = "logCPM"
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )


miRNACl1 <- "ebv-miR-BART15"

df <- data.frame(
  expr = as.numeric(logcpm1[miRNACl1, ]),
  group = factor(group_edgeRCl3, levels = c("Cl3", "Other"))
)
df$group <- factor(df$group,
                   levels = c("Cl3", "Other"),
                   labels = c("Proliferation", "Other"))


ggboxplot(df, x = "group", y = "expr",
          fill = "group",
          palette = c("Proliferation" = "palegreen1",
                      "Other" = "cornsilk1")) +
  
  
  annotate("text",
           x = 1.5,   
           y = max(df$expr) * 1.1,
           label = "FDR = 0.027",
           size = 5,
           fontface = "bold") +
  
  labs(
    title = "Expression of ebv-miR-BART15",
    x = "",
    y = "logCPM"
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )


# ESTIMATE
xcelldf <- read.table("matrix2_xcell.txt", header = TRUE, row.names = 1)
colnames(xcelldf)[colnames(xcelldf) == "OG60"] <- "OE78"

xcelldf1 <- xcelldf[,colnames(xcelldf) %in% colnames(tab_complete_out), drop = FALSE]
dim(xcelldf1)

metadata_ge <- metadata_edgeR[rownames(metadata_edgeR) %in% colnames(xcelldf1), , drop = FALSE]

metadata_ge <- metadata_ge[match(colnames(xcelldf1), rownames(metadata_ge)), , drop = FALSE]

data_batch <- data_batch[colnames(xcelldf1), , drop = FALSE]

xcelldf1_cs <- ComBat_seq(
  counts = as.matrix(xcelldf1),
  batch = data_batch$DATA)

dim(xcelldf1_cs)
dgege <- DGEList(counts = xcelldf1_cs, group = group_edgeRO)

genes <- rownames(dgege)

map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = genes,
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)

dgege$genes <- data.frame(
  ENSEMBL = genes,
  SYMBOL = map$SYMBOL[match(genes, map$ENSEMBL)]
)


keep <- !is.na(dgege$genes$SYMBOL)
dgege <- dgege[keep, , keep.lib.sizes = FALSE]


dgege <- calcNormFactors(dgege)
logcpm <- cpm(dgege, log = TRUE, prior.count = 1)


mat <- as.data.frame(logcpm)
mat$SYMBOL <- dgege$genes$SYMBOL

mat_collapsed <- mat %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), sum)) %>%
  as.data.frame()

rownames(mat_collapsed) <- mat_collapsed$SYMBOL
mat_collapsed$SYMBOL <- NULL


res <- immunedeconv::deconvolute(mat_collapsed, "estimate")


purity <- res %>%
  dplyr::filter(cell_type == "tumor purity") %>%
  tidyr::pivot_longer(-cell_type, names_to = "sample", values_to = "purity") %>%
  dplyr::select(-cell_type)

meta <- data.frame(
  sample = names(group_edgeRO),
  cluster = group_edgeRO
)

df <- merge(purity, meta, by = "sample")
df$cluster <- recode(df$cluster,
                     "Cl1" = "Immune-active",
                     "Cl2" = "Defence response",
                     "Cl3" = "Proliferation"
)

df$cluster <- factor(df$cluster,
                     levels = c("Immune-active", "Defence response", "Proliferation")
)




comparisons <- list(
  c("Immune-active", "Defence response"),
  c("Immune-active", "Proliferation"),
  c("Defence response", "Proliferation")
)

ggplot(df, aes(x = cluster, y = purity, fill = cluster)) +
  
  # violin
  geom_violin(
    trim = FALSE,
    color = "black",
    alpha = 0.8
  ) +
  
  # boxplot
  geom_boxplot(
    width = 0.12,
    fill = "white",
    outlier.shape = NA,
    color = "black"
  ) +
  
  # test
  stat_compare_means(
    method = "kruskal.test",
    label.y = 0.95
  ) +
  
  # pairwise
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    step.increase = 0.08
  ) +
  
  # custom colors
  scale_fill_manual(values = c(
    "Immune-active" = "indianred1",
    "Defence response" = "deepskyblue",
    "Proliferation" = "palegreen1"
  )) +
  
  labs(
    title = "Tumor purity across clusters",
    x = "",
    y = "Tumor purity"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      face = "bold",
      size = 12
    ),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )


#### de genes
dgege <- DGEList(counts = xcelldf1_cs, group = group_edgeRCl1)
keepge <- filterByExpr(dgege)
ge_fil <- xcelldf1_cs[keepge, ]
dim(xcelldf1_cs) #56145 87
dim(ge_fil)        #13882 87
dgege <- DGEList(counts = ge_fil, group = group_edgeRCl1)
logcpm_beforege <- cpm(dgege, log=TRUE)
dgege <- calcNormFactors(dgege, method = "TMM")
logcpmge <- cpm(dgege, log=TRUE)
designge <- model.matrix(~ group_edgeRCl1)
dgege <- estimateDisp(dgege, designge)
fitge <- glmQLFit(dgege, designge)
qlfge.Cl1 <- glmQLFTest(fitge, coef=2) 
degge.Cl1 <- topTags(qlfge.Cl1, n=20000, adjust.method = "BH", sort.by = "PValue")$table 
degge.Cl1_rawP <- degge.Cl1[degge.Cl1$FDR < 0.05, ]
up.genesge.Cl1_rawP <- row.names(degge.Cl1_rawP[degge.Cl1_rawP$logFC < -1.0, ])
down.genesge.Cl1_rawP <- row.names(degge.Cl1_rawP[degge.Cl1_rawP$logFC > 1.0, ])


## de geni Cl3
dgege <- DGEList(counts = xcelldf1_cs, group = group_edgeRCl3)
keepge <- filterByExpr(dgege)
ge_fil <- xcelldf1_cs[keepge, ]
dgege <- DGEList(counts = ge_fil, group = group_edgeRCl3)
logcpm_beforege <- cpm(dgege, log=TRUE)
dgege <- calcNormFactors(dgege, method = "TMM")
logcpmge <- cpm(dgege, log=TRUE)
designge <- model.matrix(~ group_edgeRCl3)

dgege <- estimateDisp(dgege, designge)
fitge <- glmQLFit(dgege, designge)
qlfge.Cl3 <- glmQLFTest(fitge, coef=2) 
degge.Cl3 <- topTags(qlfge.Cl3, n=20000, adjust.method = "BH", sort.by = "PValue")$table 
degge.Cl3_rawP <- degge.Cl3[degge.Cl3$FDR < 0.05, ]
up.genesge.Cl3_rawP <- row.names(degge.Cl3_rawP[degge.Cl3_rawP$logFC < -1.0, ])
down.genesge.Cl3_rawP <- row.names(degge.Cl3_rawP[degge.Cl3_rawP$logFC > 1.0, ])
up.genes.Cl1_rawP <- tolower(up.genes.Cl1_rawP) 
up.genes.Cl3_rawP <- tolower(up.genes.Cl3_rawP)

## Cl1
selected_miRNA <- unique(c(up.genes.Cl1_rawP, up.genes.Cl3_rawP))
mirbase <- read.table("interactions_human-ebv.microT.mirbase.txt", header = TRUE)
selected_miRNA_lower <- tolower(selected_miRNA)
mirbase$mirna_lower <- tolower(mirbase$mirna)

pattern <- paste0("^(", paste(selected_miRNA_lower, collapse = "|"), ")(-3p|-5p)?$")

matches <- mirbase$mirna[grepl(pattern, mirbase$mirna_lower)]

mirbase_filtered <- mirbase[grepl(pattern, mirbase$mirna_lower), ]

mirbase_filtered70 <- mirbase_filtered[mirbase_filtered$interaction_score > 0.70,, drop = FALSE]
mirbasebasrt21_70 <- mirbase_filtered70[mirbase_filtered70$mirna_lower %in% up.genes.Cl1_rawP, ,drop = FALSE]
mirbasebasrt21_70genes <- as.vector(mirbasebasrt21_70$ensembl_gene_id)

genes21 <- intersect(down.genesge.Cl1_rawP, mirbasebasrt21_70genes)
mirbasebasrt21_70genessymbols1 <- unname(mapIds(org.Hs.eg.db, keys = genes21, keytype = "ENSEMBL", column="SYMBOL"))


interactions21 <- mirbase_filtered70[
  mirbase_filtered70$ensembl_gene_id %in% genes21 &
    mirbase_filtered70$mirna_lower %in% up.genes.Cl1_rawP,
  c("mirna", "ensembl_gene_id")
]
interactions21$symbol <- mapIds(org.Hs.eg.db,
                                keys = interactions21$ensembl_gene_id,
                                keytype = "ENSEMBL",
                                column = "SYMBOL")
edges <- interactions21[, c("mirna", "symbol")]
colnames(edges) <- c("from", "to")
nodes <- data.frame(
  name = unique(c(edges$from, edges$to)),
  type = ifelse(unique(c(edges$from, edges$to)) %in% edges$from, "miRNA", "gene")
)
g21 <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
ggraph(g21, layout = "fr") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(aes(color = type), size = 4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c(miRNA = "red", gene = "skyblue2")) +
  theme_void()

##Cl3
entrez21 <- mapIds(org.Hs.eg.db,
                   keys = mirbasebasrt21_70genessymbols1,
                   keytype = "SYMBOL",
                   column = "ENTREZID",
                   multiVals = "first") |> na.omit()

kegg21 <- enrichKEGG(gene = entrez21,
                     organism = "hsa",
                     pvalueCutoff = 0.1)

kegg_df <- as.data.frame(kegg21@result)

ggplot(kegg_df[1:14, ], aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_col(fill = "skyblue3", color = "black") +
  coord_flip() +
  labs(title = "KEGG enrichment",
       x = "Pathway",
       y = "-log10(p-value)") +
  theme_minimal(base_size = 18)

#react
react21 <- enrichPathway(gene = entrez21,
                         organism = "human",
                         pvalueCutoff = 0.1)

react_df <- as.data.frame(react21@result)

ggplot(react_df[1:14, ],  #package stringr required for a good quality/readability of the plot
       aes(x = reorder(str_wrap(Description, width = 50), -pvalue),
           y = -log10(pvalue))) +
  geom_col(fill = "lightsalmon1", color = "black") +
  coord_flip() +
  labs(title = "Reactome enrichment",
       x = "Pathway",
       y = "-log10(p-value)") +
  theme_minimal(base_size = 18)

### hallmark
pathways.hallmark <- gmtPathways("h.all.v2025.1.Hs.symbols.gmt")
hallmark_t2g <- stack(pathways.hallmark)
colnames(hallmark_t2g) <- c("gene", "term")  # gene = SYMBOL, term = hallmark name
hallmark_t2g <- hallmark_t2g[, c("term", "gene")]


e21_hallmark <- enricher(
  gene = mirbasebasrt21_70genessymbols1,
  TERM2GENE = hallmark_t2g,
  pvalueCutoff = 0.1
)
res21_hallmark <- as.data.frame(e21_hallmark@result)
res21_hallmark_sig <- res21_hallmark[res21_hallmark$pvalue < 0.1, ]
res21_hallmark_sig$Description_clean <- gsub("^HALLMARK_", "", res21_hallmark_sig$Description)
ggplot(res21_hallmark_sig,
       aes(x = reorder(Description_clean, -pvalue),
           y = -log10(pvalue))) +
  geom_col(fill = "lightpink1", color = "black") +
  coord_flip() +
  labs(title = "Hallmark enrichment",
       x = "Hallmark",
       y = "-log10(pvalue)") +
  theme_minimal(base_size = 18)


### Cl3 network
mirbasebasrt21_70 <- mirbase_filtered70[mirbase_filtered70$mirna_lower %in% up.genes.Cl3_rawP, ,drop = FALSE]
mirbasebasrt21_70genes <- as.vector(mirbasebasrt21_70$ensembl_gene_id)

genes21 <- intersect(down.genesge.Cl3_rawP, mirbasebasrt21_70genes)
mirbasebasrt21_70genessymbols1 <- unname(mapIds(org.Hs.eg.db, keys = genes21, keytype = "ENSEMBL", column="SYMBOL"))


interactions21 <- mirbase_filtered70[
  mirbase_filtered70$ensembl_gene_id %in% genes21 &
    mirbase_filtered70$mirna_lower %in% up.genes.Cl3_rawP,
  c("mirna", "ensembl_gene_id")
]
interactions21$symbol <- mapIds(org.Hs.eg.db,
                                keys = interactions21$ensembl_gene_id,
                                keytype = "ENSEMBL",
                                column = "SYMBOL")
edges <- interactions21[, c("mirna", "symbol")]
colnames(edges) <- c("from", "to")
nodes <- data.frame(
  name = unique(c(edges$from, edges$to)),
  type = ifelse(unique(c(edges$from, edges$to)) %in% edges$from, "miRNA", "gene")
)
g21 <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
ggraph(g21, layout = "fr") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(aes(color = type), size = 4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c(miRNA = "red", gene = "skyblue2")) +
  theme_void()

##kegg immune-active
entrez21 <- mapIds(org.Hs.eg.db,
                   keys = mirbasebasrt21_70genessymbols1,
                   keytype = "SYMBOL",
                   column = "ENTREZID",
                   multiVals = "first") |> na.omit()

kegg21 <- enrichKEGG(gene = entrez21,
                     organism = "hsa",
                     pvalueCutoff = 0.1)

kegg_df <- as.data.frame(kegg21@result)

ggplot(kegg_df, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_col(fill = "skyblue3", color = "black") +
  coord_flip() +
  labs(title = "KEGG enrichment",
       x = "Pathway",
       y = "-log10(p-value)") +
  theme_minimal(base_size = 18)

#react immune active
react21 <- enrichPathway(gene = entrez21,
                         organism = "human",
                         pvalueCutoff = 0.1)

react_df <- as.data.frame(react21@result)

ggplot(react_df[1:14, ],
       aes(x = reorder(str_wrap(Description, width = 50), -pvalue),
           y = -log10(pvalue))) +
  geom_col(fill = "lightsalmon1", color = "black") +
  coord_flip() +
  labs(title = "Reactome enrichment",
       x = "Pathway",
       y = "-log10(p-value)") +
  theme_minimal(base_size = 18)


##hallmark proliferation
e21_hallmark <- enricher(
  gene = mirbasebasrt21_70genessymbols1,
  TERM2GENE = hallmark_t2g,
  pvalueCutoff = 0.1
)
res21_hallmark <- as.data.frame(e21_hallmark@result)
res21_hallmark_sig <- res21_hallmark[res21_hallmark$pvalue < 0.1, ]
res21_hallmark_sig$Description_clean <- gsub("^HALLMARK_", "", res21_hallmark_sig$Description)
ggplot(res21_hallmark_sig,
       aes(x = reorder(Description_clean, -pvalue),
           y = -log10(pvalue))) +
  geom_col(fill = "lightpink1", color = "black") +
  coord_flip() +
  labs(title = "Hallmark enrichment",
       x = "Hallmark",
       y = "-log10(pvalue)") +
  theme_minimal(base_size = 18)



### Viral load EBV
load_table <- read_xlsx("viral_load.xlsx")
load_table$Campione[load_table$Campione == "OG60"] <- "OE78"
load_table <- load_table[load_table$Campione %in% rownames(data_batch),]
load_table$Cluster <- metadata_ge[load_table$Campione, "Cluster"]
load_table$Cluster <- dplyr::recode(
  load_table$Cluster,
  "Cl1" = "Immune-active",
  "Cl2" = "Defence response",
  "Cl3" = "Proliferation"
)
load_table$Cluster <- factor(load_table$Cluster, levels = c("Immune-active", "Defence response", "Proliferation"))

ggboxplot(load_table, x = "Cluster", y = "EBV_miRNA_load",
          fill = "Cluster",
          palette = c("Immune-active" = "indianred1",
                      "Defence response" = "deepskyblue",
                      "Proliferation" = "palegreen1")) +
  
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(
                       c("Immune-active", "Defence response"),
                       c("Immune-active", "Proliferation"),
                       c("Defence response", "Proliferation")
                     ),
                     label = "p.format") +
  
  ggtitle("EBV miRNA load") +
  
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.position = "right"
  )
