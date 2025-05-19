# Detailed analysis of GPR89 and FRMPD2 forebrain cells

# Expression pattern of FRMDP2 and GPR89 in zebrafish forebrain.
```{r}
# Heatmap with Z-score across tissues:
meta <- all.data.combined.selected.detailed.forebrain@meta.data
HSDs <- c("frmpd2", "GPR89B")
expr <- all.data.combined.selected.detailed.forebrain@assays$RNA$data[HSDs,] %>% as.matrix()
heatmap_df <- bind_cols(data.frame(clusters=meta$hc_subset_euclidean_forebrain_10_renamed), t(expr)) %>% group_by(clusters) %>% summarize_all(sum) %>% as.data.frame()
rownames(heatmap_df) <- heatmap_df$clusters
heatmap_df <- heatmap_df[,-1] %>% as.matrix()

pdf("final_analysis/plots/HSDs_expression_heatmap_forebrain_GPR89_FRMPD2.pdf", height=8, width=8)
pheatmap::pheatmap(heatmap_df,cluster_rows = FALSE, cluster_cols = TRUE, scale = "column")
dev.off()

```

# Re-cluster to focus on telencephalon.
```{r}
# Grab the cells of the groups of interest.
all.data.combined.selected.detailed <- subset(all.data.combined.selected.filtered, subset= genotype %in% c("gpr89KO", "GPR89B",
                                                                                                           "frmpd2KO", "FRMPD2B",
                                                                                                           "scrambled", "GFP")) # 31,925 cells

# Grab forebrain cells.
all.data.combined.selected.detailed.forebrain <- subset(all.data.combined.selected.detailed, subset= specific.celltypes.final %in% c("forebrain", "thalamus","pallium", "hypothalamus","olfactory-organ")) # 10,040 cells

# Get distance between all cells from corrected-PC embeddings:
d3 <- dist(all.data.combined.selected.detailed.forebrain@reductions$harmony@cell.embeddings, method = "euclidean")

# Perform the hierarchical clustering using the euclidean distance:
h_euclidean3 <- fastcluster::hclust(d3, method = "ward.D2") # use the hclust from 'fastcluster' to deal with large datasets
saveRDS(h_euclidean3, file="final_analysis/saved_datasets/h_euclidean3.rds")

# Cut the tree at a high K:
all.data.combined.selected.detailed.forebrain$hc_subset_euclidean_forebrain_10 <- cutree(h_euclidean3, k= 10)
all.data.combined.selected.detailed.forebrain$hc_subset_euclidean_forebrain_10 <- as.factor(all.data.combined.selected.detailed.forebrain$hc_subset_euclidean_forebrain_10)

# Obtain gene markers per cluster:
Idents(all.data.combined.selected.detailed.forebrain) <- all.data.combined.selected.detailed.forebrain$hc_subset_euclidean_forebrain_10
all.data.combined.selected.detailed.forebrain <- PrepSCTFindMarkers(all.data.combined.selected.detailed.forebrain, assay="SCT", verbose = TRUE)

gene.markers.selected.k10 <- FindAllMarkers(object = all.data.combined.selected.detailed.forebrain, 
                                            assay="SCT", 
                                            logfc.threshold = 0.10, 
                                            test.use = "MAST", 
                                            only.pos = TRUE, 
                                            min.pct = 0.10, 
                                            random.seed = 1234, 
                                            verbose = TRUE,
                                            recorrect_umi=FALSE)

write.csv(gene.markers.selected.k10, file="final_analysis/data/gene.markers.selected.k10.csv")
write.csv(as.data.frame(table(all.data.combined.selected.detailed.forebrain$hc_subset_euclidean_forebrain_10)), file="final_analysis/data/gene.markers.selected.k10.clustercounts.csv")

# Run UMAP:
all.data.combined.selected.detailed.forebrain <- RunUMAP(
  all.data.combined.selected.detailed.forebrain,
  reduction.name = 'UMAP',
  dims= 1:50,
  reduction = 'harmony',
  assay="SCT",
  seed.use = 1234)

# Plot UMAP with clusters (basic):
Idents(all.data.combined.selected.detailed.forebrain) <- all.data.combined.selected.detailed.forebrain$hc_subset_euclidean_forebrain_10
pdf('final_analysis/plots/UMAP_subset_euclidean_detailed_forebrain_k10.pdf', height=10, width=10)
DimPlot_scCustom(seurat_object = all.data.combined.selected.detailed.forebrain)
dev.off()

# Plot UMAP with clusters (with labels for publication):

UMAP_centers_cluster <- tibble(
    UMAP_1 = as.data.frame(all.data.combined.selected.detailed.forebrain@reductions$UMAP@cell.embeddings)$UMAP_1,
    UMAP_2 = as.data.frame(all.data.combined.selected.detailed.forebrain@reductions$UMAP@cell.embeddings)$UMAP_2,
    cluster = all.data.combined.selected.detailed.forebrain@meta.data$hc_subset_euclidean_forebrain_10
  ) %>%
  group_by(cluster) %>%
  summarize(x = median(UMAP_1), y = median(UMAP_2))

plot_umap_forebrain <- bind_cols(all.data.combined.selected.detailed.forebrain@meta.data, as.data.frame(all.data.combined.selected.detailed.forebrain@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = hc_subset_euclidean_forebrain_10)) +
  geom_point(size = 0.05) +
  geom_label(
    data = UMAP_centers_cluster,
    mapping = aes(x, y, label = cluster),
    size = 2.5,
    fill = 'white',
    color = 'black',
    fontface = 'bold',
    alpha = 0.5,
    label.size = 0,
    show.legend = FALSE
  ) +
  theme_bw() +
  scale_color_manual(
    values = custom_colors$discrete,
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed()

ggsave(
  'final_analysis/plots/UMAP_subset_euclidean_detailed_forebrain_k10_labels.pdf',
  plot_umap_forebrain,
  height = 6,
  width = 6
)

# Save current object.
saveRDS(all.data.combined.selected.detailed.forebrain, file="final_analysis/saved_datasets/all.data.combined.selected.detailed.forebrain.rds")


# Re-name clusters:
forebrain_clusters <- c("ventral_tel", "dorsal_tel", "inner_dien",
                        "dorsal_tel", "thalamus", "hypothalamus",
                        "habenula", "TBD", "epiphysis", "olfactory")
names(forebrain_clusters) <- levels(all.data.combined.selected.detailed.forebrain$hc_subset_euclidean_forebrain_10)

Idents(all.data.combined.selected.detailed.forebrain) <- "hc_subset_euclidean_forebrain_10"
all.data.combined.selected.detailed.forebrain <- RenameIdents(all.data.combined.selected.detailed.forebrain, forebrain_clusters)
all.data.combined.selected.detailed.forebrain$hc_subset_euclidean_forebrain_10_renamed <- Idents(all.data.combined.selected.detailed.forebrain)

# Remove cell types not of interest:
all.data.combined.selected.detailed.forebrain <- subset(all.data.combined.selected.detailed.forebrain, hc_subset_euclidean_forebrain_10_renamed %in% c("ventral_tel", "dorsal_tel", "inner_dien", "thalamus", "hypothalamus", "habenula", "epiphysis", "olfactory")) # 9,866
all.data.combined.selected.detailed.forebrain$hc_subset_euclidean_forebrain_10_renamed <- droplevels(all.data.combined.selected.detailed.forebrain$hc_subset_euclidean_forebrain_10_renamed)

# Order clusters for plotting:
all.data.combined.selected.detailed.forebrain$hc_subset_euclidean_forebrain_10_renamed <- factor(all.data.combined.selected.detailed.forebrain$hc_subset_euclidean_forebrain_10_renamed, levels=c("dorsal_tel", "ventral_tel", "olfactory","epiphysis", "habenula", "thalamus", "inner_dien", "hypothalamus"))

# Run dimensionality reductions:
all.data.combined.selected.detailed.forebrain <- RunTSNE(
  all.data.combined.selected.detailed.forebrain,
  reduction.name = 'TSNE',
  dims= 1:50,
  reduction = 'harmony',
  assay="SCT",
  seed.use = 1234)

# Plot UMAP/tSNE with clusters (with labels for publication):

TSNE_centers_cluster <- tibble(
    tSNE_1 = as.data.frame(all.data.combined.selected.detailed.forebrain@reductions$TSNE@cell.embeddings)$tSNE_1,
    tSNE_2 = as.data.frame(all.data.combined.selected.detailed.forebrain@reductions$TSNE@cell.embeddings)$tSNE_2,
    cluster = all.data.combined.selected.detailed.forebrain@meta.data$hc_subset_euclidean_forebrain_10_renamed
  ) %>%
  group_by(cluster) %>%
  summarize(x = median(tSNE_1), y = median(tSNE_2))

plot_forebrain_tsne_final <- bind_cols(all.data.combined.selected.detailed.forebrain@meta.data, as.data.frame(all.data.combined.selected.detailed.forebrain@reductions$TSNE@cell.embeddings)) %>%
  ggplot(aes(tSNE_1, tSNE_2, color = hc_subset_euclidean_forebrain_10_renamed)) +
  geom_point(size = 0.05) +
  geom_label(
    data = TSNE_centers_cluster,
    mapping = aes(x, y, label = cluster),
    size = 2.5,
    fill = 'white',
    color = 'black',
    fontface = 'bold',
    alpha = 0.5,
    label.size = 0,
    show.legend = FALSE
  ) +
  theme_bw() +
  scale_color_manual(
    values = custom_colors$discrete,
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed()

ggsave(
  'final_analysis/plots/TSNE_subset_euclidean_detailed_forebrain_k10_labels.pdf',
  plot_forebrain_tsne_final,
  height = 6,
  width = 6
)

```

# DGE cell-type-specific:
```{r}
# Grab cells from telencephalon
all.data.combined.selected.detailed.telencephalon <- subset(all.data.combined.selected.detailed.forebrain.sampled, subset= hc_subset_euclidean_forebrain_10_renamed %in% c("dorsal_tel", "ventral_tel"))
all.data.combined.selected.detailed.telencephalon$hc_subset_euclidean_forebrain_10_renamed <- droplevels(all.data.combined.selected.detailed.telencephalon$hc_subset_euclidean_forebrain_10_renamed)

# Subsample to balance numbers of cells across genotypes.
all.data.combined.selected.detailed.telencephalon.sample <- balanceSamples(object = all.data.combined.selected.detailed.telencephalon, group = "genotype")

# Pseudobulk the counts based on donor-condition-celltype
all.data.combined.selected.detailed.telencephalon.sampled.pseudo <- AggregateExpression(all.data.combined.selected.detailed.forebrain.sample, assays = "RNA", return.seurat = T, group.by = c("genotype", "sample"))

# DGE:
tel.deseq2 <- DESeqDataSetFromMatrix(countData = all.data.combined.selected.detailed.telencephalon.sampled.pseudo@assays$RNA$counts+1,colData= all.data.combined.selected.detailed.telencephalon.sampled.pseudo@meta.data, design= ~genotype)
tel.deseq2 <- estimateSizeFactors(tel.deseq2)
tel.deseq2 <- estimateDispersionsGeneEst(tel.deseq2)
dispersions(tel.deseq2) <- mcols(tel.deseq2)$dispGeneEst
tel.deseq2 <- nbinomWaldTest(tel.deseq2)

res.tel.gpr89KO <- as.data.frame(results(tel.deseq2, contrast = c("genotype", "gpr89KO", "scrambled")))
res.tel.GPR89B <- as.data.frame(results(tel.deseq2, contrast = c("genotype", "GPR89B", "GFP")))

res.tel.frmpd2KO <- as.data.frame(results(tel.deseq2, contrast = c("genotype", "frmpd2KO", "scrambled")))
res.tel.FRMPD2B <- as.data.frame(results(tel.deseq2, contrast = c("genotype", "FRMPD2B", "GFP")))

# Genes in common for FRMPD2:
res.tel.frmpd2KO$direction <- ifelse(res.tel.frmpd2KO$log2FoldChange < 0, "down",
                                     ifelse(res.tel.frmpd2KO$log2FoldChange > 0, "up", "none"))
res.tel.FRMPD2B$direction <- ifelse(res.tel.FRMPD2B$log2FoldChange < 0, "down",
                                     ifelse(res.tel.FRMPD2B$log2FoldChange > 0, "up", "none"))

merged_FRMPD2 <- merge(res.tel.frmpd2KO, res.tel.FRMPD2B, by="row.names")

merged_FRMPD2$group <- ifelse(merged_FRMPD2$direction.x == "up" & merged_FRMPD2$direction.y == "up", "same_up", ifelse(merged_FRMPD2$direction.x == "down" & merged_FRMPD2$direction.y == "down", "same_down", "other"))

write.csv(merged_FRMPD2[merged_FRMPD2$group == "same_up",], file="final_analysis/data/Forebrain_FRMPD2_DEGs_same_up.csv")

write.csv(merged_FRMPD2[merged_FRMPD2$group == "same_down",], file="final_analysis/data/Forebrain_FRMPD2_DEGs_same_down.csv")

p_FRMPD2 <- ggplot() +
  geom_point(data= merged_FRMPD2[merged_FRMPD2$group == "same_up",], aes(x= log2FoldChange.x, y= log2FoldChange.y), col= "red") +
  geom_point(data= merged_FRMPD2[merged_FRMPD2$group == "same_down",], aes(x= log2FoldChange.x, y= log2FoldChange.y), col= "blue")  +
  theme_bw() + xlab("log2(FC) frmpd2 KO") + ylab("log2(FC) FRMPD2B")

# Distribution of genes:
# Up FRMPD2B, up in frmpd2 KO: 9355 genes (31.2%)
# Down FRMPD2B, down in frmpd2 KO: 15123 genes (50.5%)
# Total genes: 29945 genes

ggsave("final_analysis/plots/Correlation_telencephalon_FRMPD2.pdf", p_FRMPD2, height = 8, width= 8)

keep_background <- rowSums(all.data.combined.selected.detailed.telencephalon$RNA$counts) >= 1
background_genes <- rownames(all.data.combined.selected.detailed.telencephalon$RNA$counts)[keep_background]


GO_common_up_FRMPD2 <- enrichGO(gene      = merged_FRMPD2$Row.names[merged_FRMPD2$group == "same_up"],
                        OrgDb         = org.Dr.eg.db,
                        keyType       = "SYMBOL",
                        universe      = background_genes,
                        ont           = "BP",
                        pvalueCutoff  = 0.1,
                        pAdjustMethod = "BH",
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)
write.csv(GO_common_up_FRMPD2@result, file="final_analysis/plots/Correlation_telencephalon_FRMPD2_barplot_same_up_GO.csv")

GO_common_down_FRMPD2 <- enrichGO(gene      = merged_FRMPD2$Row.names[merged_FRMPD2$group == "same_down"],
                        OrgDb         = org.Dr.eg.db,
                        keyType       = "SYMBOL",
                        universe      = background_genes,
                        ont           = "BP",
                        pvalueCutoff  = 0.1,
                        pAdjustMethod = "BH",
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)
write.csv(GO_common_down_FRMPD2@result, file="final_analysis/plots/Correlation_telencephalon_FRMPD2_barplot_same_down_GO.csv")

pdf("final_analysis/plots/Correlation_telencephalon_FRMPD2_barplot_same_up_GO.pdf")
mutate(GO_common_up_FRMPD2, logq = -log(qvalue, base=10)) %>% 
    barplot(x="logq", showCategory= 15, ordered=T, color="qvalue")
dev.off()

pdf("final_analysis/plots/Correlation_telencephalon_FRMPD2_barplot_down_GO.pdf")
mutate(GO_common_down_FRMPD2, logq = -log(qvalue, base=10)) %>% 
    barplot(x="logq", showCategory= 15, ordered=T, color="qvalue")
dev.off()

write.csv(GO_common_up_FRMPD2@result, file="final_analysis/data/Correlation_telencephalon_FRMPD2_barplot_up_GO.csv")

write.csv(GO_common_down_FRMPD2@result, file="final_analysis/data/Correlation_telencephalon_FRMPD2_barplot_down_GO.csv")

# Genes opposite for GPR89:

res.tel.gpr89KO$direction <- ifelse(res.tel.gpr89KO$log2FoldChange < 0, "down",
                                     ifelse(res.tel.gpr89KO$log2FoldChange > 0, "up", "none"))
res.tel.GPR89B$direction <- ifelse(res.tel.GPR89B$log2FoldChange < 0, "down",
                                     ifelse(res.tel.GPR89B$log2FoldChange > 0, "up", "none"))

merged_GPR89 <- merge(res.tel.gpr89KO, res.tel.GPR89B, by="row.names")

merged_GPR89$group <- ifelse(merged_GPR89$direction.x == "up" & merged_GPR89$direction.y == "down", "opposite_low", ifelse(merged_GPR89$direction.x == "down" & merged_GPR89$direction.y == "up", "opposite_high", "other"))

write.csv(merged_GPR89[merged_GPR89$group == "opposite_low",], file="final_analysis/data/Forebrain_GPR89_DEGs_low.csv")

write.csv(merged_GPR89[merged_GPR89$group == "opposite_high",], file="final_analysis/data/Forebrain_GPR89_DEGs_high.csv")

p_GPR89 <- ggplot() +
  geom_point(data= merged_GPR89[merged_GPR89$group=="opposite_high",], aes(x= log2FoldChange.x, y= log2FoldChange.y), col= "red") +
  geom_point(data= merged_GPR89[merged_GPR89$group=="opposite_low",], aes(x= log2FoldChange.x, y= log2FoldChange.y), col= "blue") +
  xlab("log2(FC) gpr89 KO") + ylab("log2(FC) GPR89B") +
  theme_bw()

# Distribution of genes:
# Up GPR89B, down in gpr89 KO: 2530 genes (8.4%)
# Down GPR89B, up in gpr89 KO: 5865 genes (19.6%)
# Total genes: 29945 genes

ggsave("final_analysis/plots/Correlation_telencephalon_GPR89.pdf", p_GPR89, height = 8, width= 8)

GO_high_GPR89 <- enrichGO(gene      = merged_GPR89$Row.names[merged_GPR89$group=="opposite_high"],
                        OrgDb         = org.Dr.eg.db,
                        universe= background_genes,
                        keyType       = "SYMBOL",
                        ont           = "BP",
                        pvalueCutoff  = 0.1,
                        pAdjustMethod = "BH",
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)
write.csv(GO_high_GPR89@result, file="final_analysis/plots/Correlation_telencephalon_GPR89_barplot_opposite_high_GO.csv")

GO_low_GPR89 <- enrichGO(gene      = merged_GPR89$Row.names[merged_GPR89$group=="opposite_low"],
                        OrgDb         = org.Dr.eg.db,
                        keyType       = "SYMBOL",
                        ont           = "BP",
                        pvalueCutoff  = 0.1,
                        pAdjustMethod = "BH",
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)
write.csv(GO_low_GPR89@result, file="final_analysis/plots/Correlation_telencephalon_GPR89_barplot_opposite_low_GO.csv")


pdf("final_analysis/plots/Correlation_telencephalon_GPR89_barplot_opposite_high_GO.pdf")
mutate(GO_high_GPR89, logq = -log(qvalue, base=10)) %>% 
    barplot(x="logq", showCategory= 15, ordered=T, color="qvalue")
dev.off()

pdf("final_analysis/plots/Correlation_telencephalon_GPR89_barplot_opposite_low_GO.pdf")
mutate(GO_low_GPR89, logq = -log(qvalue, base=10)) %>% 
    barplot(x="logq", showCategory= 15, ordered=T, color="qvalue")
dev.off()

write.csv(GO_high_GPR89@result, file="final_analysis/data/Correlation_telencephalon_GPR89_barplot_high_GO.csv")

write.csv(GO_low_GPR89@result, file="final_analysis/data/Correlation_telencephalon_GPR89_barplot_low_GO.csv")

```

# Cell cycle and progenitor markers.
```{r}
progenitors <- c("sox19a", "sox2", "rpl5a", "npm1a", "s100b", "dla")
neuron <- c("elavl3", "elavl4", "tubb5")

# Progenitors:

counts_progenitors <- Matrix::Matrix(all.data.combined.selected.detailed.forebrain@assays$RNA@counts[progenitors, ])
presence_prog <- as.integer(colSums(counts_progenitors) >= 1)
all.data.combined.selected.detailed.forebrain$presence_prog <- presence_prog

all.data.combined.selected.detailed.forebrain$genotype <- relevel(all.data.combined.selected.detailed.forebrain$genotype, ref="GFP")

all.data.combined.selected.detailed.forebrain$presence_prog <- as.factor(all.data.combined.selected.detailed.forebrain$presence_prog)
all.data.combined.selected.detailed.forebrain$presence_prog <- relevel(all.data.combined.selected.detailed.forebrain$presence_prog, ref="1")

glm_prog <- glm(presence_prog ~ genotype, data = all.data.combined.selected.detailed.forebrain@meta.data, family = binomial)
summary(glm_prog)

#                  Estimate Std. Error z value Pr(>|z|)
#(Intercept)        1.21698    0.05071  23.998  < 2e-16 ***
#genotypeGPR89B     0.29293    0.11525   2.542   0.0110 *
#genotypegpr89KO   -0.28439    0.07041  -4.039 5.37e-05 ***
#genotypescrambled -0.11735    0.06270  -1.871   0.0613 .

coef_prog <- summary(glm_prog)$coef
z_critical <- qnorm(0.975)
coef_df_neuron <- data.frame(
  Genotype = rownames(coef_prog),
  Estimate = coef_prog[, 1],
  OR= exp(coef_prog[,1]),
  Lower_CI = exp(coef_prog[, 1]) * exp(-z_critical * coef_prog[, 2]),
  Upper_CI = exp(coef_prog[, 1]) * exp(z_critical * coef_prog[, 2])
)
coef_df_neuron <- coef_df_neuron[coef_df_neuron$Genotype != "(Intercept)",]

#                           Genotype   Estimate        OR  Lower_CI  Upper_CI
#genotypeGPR89B       genotypeGPR89B  0.2929253 1.3403426 1.0693466 1.6800151
#genotypegpr89KO     genotypegpr89KO -0.2843876 0.7524749 0.6554767 0.8638271
#genotypescrambled genotypescrambled -0.1173509 0.8892731 0.7864316 1.0055632

p_prog <- ggplot(coef_df_neuron, aes(x = OR, y = Genotype)) +
  geom_point() +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height=0) +
  labs(x = "OR", y = "") +
  theme_bw() + geom_vline(xintercept = 1)

ggsave("final_analysis/plots/Forestplot_progenitors_OR_GPR89.pdf", p_prog)

```

# Estimation of progenitor cells.
```{r}
progenitors <- c("sox19a", "sox2", "rpl5a", "npm1a", "s100b", "dla")
neuron <- c("elavl3", "elavl4", "tubb5")

## Estimate the % of expression of each data set per cell:
all.data.combined.selected.detailed.telencephalon.sample[["percent.prog"]] <- PercentageFeatureSet(all.data.combined.selected.detailed.telencephalon.sample, features = progenitors, assay = 'RNA')

## Tests:

# Progenitors:
counts_neuro <- Matrix::Matrix(all.data.combined.selected.detailed.telencephalon@assays$RNA@counts[neuron, ])
presence_neuro <- as.integer(colSums(counts_neuro) >= 1)
all.data.combined.selected.detailed.telencephalon$presence_neuro <- presence_neuro

table(all.data.combined.selected.detailed.telencephalon$sample, all.data.combined.selected.detailed.telencephalon$presence_neuro)

# Are there changes in the % of committed neurons across genotypes?

neuro_abundance <- table(all.data.combined.selected.detailed.telencephalon$presence_neuro, all.data.combined.selected.detailed.telencephalon$genotype)

df_GPR98B <- neuro_abundance[, c("GFP", "GPR89B")]
chisq.test(df_GPR98B) # X-squared = 69.787, df = 1, p-value < 2.2e-16

df_gpr89KO <- neuro_abundance[, c("scrambled", "gpr89KO")]
chisq.test(df_gpr89KO) # X-squared = 6.5174, df = 1, p-value = 0.01068

```

# Obtain DEGs for SRGAP2 and ARHGAP11 models in forebrain cells.
```{r}
# Grab the cells of the groups of interest.
all.data.combined.selected.detailed2 <- subset(all.data.combined.selected.filtered, subset= genotype %in% c("srgap2KO", "SRGAP2C","ARHGAP11B","scrambled", "GFP")) # 24,288 cells

# Grab forebrain cells.
all.data.combined.selected.detailed2.forebrain <- subset(all.data.combined.selected.detailed2, subset= specific.celltypes.final %in% c("forebrain", "thalamus","pallium", "hypothalamus","olfactory-organ")) # 8,067 cells

# Subsample to balance numbers of cells across genotypes.
all.data.combined.selected.detailed2.forebrain.sample <- balanceSamples(object = all.data.combined.selected.detailed2.forebrain, group = "genotype")

# Pseudobulk the counts based on donor-condition-celltype
all.data.combined.selected.detailed2.forebrain.sample.pseudo <- AggregateExpression(all.data.combined.selected.detailed2.forebrain.sample, assays = "RNA", return.seurat = T, group.by = c("genotype", "sample"))

# DGE:
forebrain.deseq2 <- DESeqDataSetFromMatrix(countData = all.data.combined.selected.detailed2.forebrain.sample.pseudo@assays$RNA$counts+1,colData= all.data.combined.selected.detailed2.forebrain.sample.pseudo@meta.data, design= ~genotype)
forebrain.deseq2 <- estimateSizeFactors(forebrain.deseq2)
forebrain.deseq2 <- estimateDispersionsGeneEst(forebrain.deseq2)
dispersions(forebrain.deseq2) <- mcols(forebrain.deseq2)$dispGeneEst
forebrain.deseq2 <- nbinomWaldTest(forebrain.deseq2)

res.SRGAP2C <- as.data.frame(results(forebrain.deseq2, contrast = c("genotype", "SRGAP2C", "GFP")))
res.srgap2KO <- as.data.frame(results(forebrain.deseq2, contrast = c("genotype", "srgap2KO", "scrambled")))

res.tel.ARHGAP11B <- as.data.frame(results(forebrain.deseq2, contrast = c("genotype", "ARHGAP11B", "GFP")))

## Genes upregulated SRGAP2C/srgap2 KO:

res.SRGAP2C$direction <- ifelse(res.SRGAP2C$log2FoldChange < 0, "down",
                                     ifelse(res.SRGAP2C$log2FoldChange > 0, "up", "none"))
res.srgap2KO$direction <- ifelse(res.srgap2KO$log2FoldChange < 0, "down",
                                     ifelse(res.srgap2KO$log2FoldChange > 0, "up", "none"))

merged_SRGAP2 <- merge(res.srgap2KO, res.SRGAP2C, by="row.names")

merged_SRGAP2$group <- ifelse(merged_SRGAP2$direction.x == "up" & merged_SRGAP2$direction.y == "up", "same_up", ifelse(merged_SRGAP2$direction.x == "down" & merged_SRGAP2$direction.y == "down", "same_down", "other"))

write.csv(merged_SRGAP2[merged_SRGAP2$group == "same_up",], file="final_analysis/data/Forebrain_SRGAP2_DEGs_same_up.csv")

write.csv(merged_FRMPD2[merged_FRMPD2$group == "same_down",], file="final_analysis/data/Forebrain_SRGAP2_DEGs_same_down.csv")

# Distribution of genes:
# Up SRGAP2C, up in srgap2 KO: 9064 genes (30.3%)
# Down SRGAP2C, down in srgap2 KO: 15210 genes (50.8%)
# Total genes: 29945 genes

keep_background <- rowSums(all.data.combined.selected.detailed2.forebrain$RNA$counts) >= 1
background_genes <- rownames(all.data.combined.selected.detailed2.forebrain$RNA$counts)[keep_background]


GO_common_up_SRGAP2 <- enrichGO(gene      = merged_SRGAP2$Row.names[merged_SRGAP2$group == "same_up"],
                        OrgDb         = org.Dr.eg.db,
                        keyType       = "SYMBOL",
                        universe      = background_genes,
                        ont           = "BP",
                        pvalueCutoff  = 0.1,
                        pAdjustMethod = "BH",
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)

pdf("final_analysis/plots/Correlation_forebrain_SRGAP2_barplot_up_GO.pdf")
mutate(GO_common_up_SRGAP2, logq = -log(qvalue, base=10)) %>% 
    barplot(x="logq", showCategory= 15, ordered=T, color="qvalue")
dev.off()

write.csv(GO_common_up_SRGAP2@result, file="final_analysis/data/Correlation_forebrain_SRGAP2_barplot_up_GO.csv")

## Upregulated in ARHGAP11B:

res.tel.ARHGAP11B$direction <- ifelse(res.tel.ARHGAP11B$log2FoldChange < 0, "down",
                                     ifelse(res.tel.ARHGAP11B$log2FoldChange > 0, "up", "none"))
res.tel.ARHGAP11B <- res.tel.ARHGAP11B[complete.cases(res.tel.ARHGAP11B),]
res.tel.ARHGAP11B$Genes <- rownames(res.tel.ARHGAP11B)

GO_up_ARHGAP11B <- enrichGO(gene      = res.tel.ARHGAP11B$Genes[res.tel.ARHGAP11B$log2FoldChange > 0 & res.tel.ARHGAP11B$padj < 0.05],
                        OrgDb         = org.Dr.eg.db,
                        keyType       = "SYMBOL",
                        universe      = background_genes,
                        ont           = "BP",
                        pvalueCutoff  = 0.1,
                        pAdjustMethod = "BH",
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)

pdf("final_analysis/plots/Correlation_forebrain_ARHGAP11B_barplot_up_GO.pdf")
mutate(GO_up_ARHGAP11B, logq = -log(qvalue, base=10)) %>% 
    barplot(x="logq", showCategory= 15, ordered=T, color="qvalue")
dev.off()

write.csv(GO_up_ARHGAP11B@result, file="final_analysis/data/Correlation_forebrain_ARHGAP11B_barplot_up_GO.csv")


```

