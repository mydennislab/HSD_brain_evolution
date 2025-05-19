# R version: 4.2.2
# BiocManager version: 3.16

# Data handling
#install.packages("tidyverse")
#install.packages("data.table")
library(tidyverse)
library(data.table)

# Visualization
#install.packages("gridExtra")
#install.packages("gplots")
#install.packages("devtools")
#library(devtools)
#devtools::install_github("kevinblighe/CorLevelPlot")
#install.packages("RColorBrewer")
library(gridExtra)
library(gplots)
library(CorLevelPlot)
library(RColorBrewer)

# Network analysis
# install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz", repos=NULL, type="source")
#BiocManager::install("DESeq2")
#install.packages("WGCNA")
# install.packages("psych")
# install.packages("igraph")
library(DESeq2)
library(WGCNA)
library(psych)
library(igraph)

# install.packages("venneuler")
# BiocManager::install("limma")
# BiocManager::install("edgeR")
# install.packages("qvalue")
# install_github("jtlovell/limmaDE2")
library(limmaDE2)

# Gene ontology over-representation
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(GO.db) # GO.db_3.16.0
library("org.Hs.eg.db", character.only = TRUE)

sessionInfo()

dir.create("plots")
dir.create("results")

c("#9e0142", "#d53e4f", "#f46d43","#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd","#5e4fa2")

# Expression analysis of SD-98 genes in BrainSpan dataset ----------------------

## 1. Loading files ------------------------------------------------------------

# Reading counts
data <- read.table(file = 'counts/allSamples_brainspan_counts.tsv', sep = '\t', header = TRUE)
data_long <- data %>% pivot_longer(-gene, names_to = "run", values_to = "counts")
rm(data)

# Reading metadata
metadata <- read.table(file = "metadata/BrainSpan_metadata.txt", sep="\t", fill=TRUE, strip.white=TRUE, quote="", header=TRUE)
data_annot <- left_join(data_long, metadata, by = c("run" = "run"))
rm(data_long)

# Removing mitochondrial genes
mt <- read.table("annotations/Mt_genes.txt", col.names = "gene_id") %>% pull()
data_nomt <- data_annot %>% filter(!gene %in% mt)
rm(data_annot)

# Selecting Geschwind samples only
data_brainspan <- data_nomt %>% 
  filter(geschwind == "yes") %>%
  filter(developmental_stage %in% c("Early prenatal", "Early mid-prenatal", "Late mid-prenatal")) %>%
  filter(body_site %in% c("DFC", "VFC", "MFC", "OFC"))
rm(data_nomt)

## 2. Data pre-processing ------------------------------------------------------

### 2.1. Removing genes with no expression -------------------------------------

# Geschwind filtered genes that were present at an RPKM of 1 in 80% of the 
# samples from at least one neocortical region at one major temporal epoch.
# We performed a similar filter using TPMs instead.

genes_to_include <- fread("brainspan_wgcna_genes.txt", header = FALSE) %>% pull(V1)
length(genes_to_include) # 17,408 genes

data_de <- data_brainspan %>%
  filter(gene %in% genes_to_include) %>%
  arrange(order) %>%
  mutate(run = factor(run, levels = unique(run))) %>%
  dplyr::select(gene, run, counts) %>%
  pivot_wider(names_from = run, values_from = counts) %>%
  column_to_rownames(var = "gene")

data_brainspan %>%
  dplyr::select(run, body_site, developmental_stage) %>%
  distinct() %>%
  ggplot(aes(body_site)) + 
  geom_bar(width = 0.8) +
  facet_wrap(~developmental_stage, ncol = 3) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
ggsave("plots/BrainSpan_metadata.pdf", width = 6, height = 3)

# converting the first column to row names
data_de_mat <- as.matrix(data_de)

### 2.2 Looking for outliers ---------------------------------------------------

# we first used the goodSamplesGenes function to remove flagged samples/genes
gsg <- goodSamplesGenes(t(data_de_mat))
summary(gsg)
gsg$allOK # TRUE
data_de_mat <- data_de_mat[gsg$goodGenes == TRUE,]
dim(data_de_mat) # 17,388 genes remaining

# we then identified outliers with hierarchical clustering and PCA

# hierarchical clustering
htree <- hclust(dist(t(data_de_mat)), method = "average")

data_de_mat.meta <-
  data.frame(run = colnames(data_de_mat)) %>%
  left_join(metadata) %>%
  mutate(developmental_stage = factor(developmental_stage, ordered = TRUE, 
                                      levels = c("Early prenatal", "Early mid-prenatal", 
                                                 "Late mid-prenatal"))) %>%
  arrange(developmental_stage) %>%
  dplyr::select(run, developmental_stage, body_site, age, order)

traitColors <- labels2colors(data_de_mat.meta$developmental_stage, 
                             colorSeq = c("#eccbae", "#cb7a5c","#5785c1"))

pdf("plots/samples_dendrogram.pdf", width = 20, height = 10)
plotDendroAndColors(htree, traitColors,
                    groupLabels = names(htree),
                    dendroLabels = colnames(data_de_mat),
                    cex.dendroLabels = 0.5,
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# PCA by developmental stage and age
pca <- prcomp(t(data_de_mat))
pca.data <- as.data.frame(pca$x)
pca.var <- pca$sdev
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

mycolors <- brewer.pal(8, "Spectral")

ggplot(pca.data, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = data_de_mat.meta$developmental_stage, color = data_de_mat.meta$age), size = 2) +
  #geom_text(label = rownames(pca.data), size = 3) +
  labs(x = paste0('PC1:', pca.var.percent[1], '%'),
       y = paste0('PC2:', pca.var.percent[2], '%')) +
  #scale_shape_manual(values = c(21,22,23,24,25,3,4,8)) +
  #scale_color_manual(values = c("#eccbae", "#cb7a5c","#5785c1")) +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom")
ggsave("plots/samples_pca.pdf", width = 4, height = 4)

# PCA by developmental stage and site
mycolors <- brewer.pal(8, "Spectral")

ggplot(pca.data, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = data_de_mat.meta$developmental_stage, color = data_de_mat.meta$body_site), size = 2) +
  #geom_text(label = rownames(pca.data), size = 3) +
  labs(x = paste0('PC1:', pca.var.percent[1], '%'),
       y = paste0('PC2:', pca.var.percent[2], '%')) +
  #scale_shape_manual(values = c(21,22,23,24,25,3,4,8)) +
  #scale_color_manual(values = c("#eccbae", "#cb7a5c","#5785c1")) +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom")
ggsave("plots/samples_pca2.pdf", width = 4, height = 4)

### 2.3 Sample filtering -------------------------------------------------------

# Removing sample outliers 
# Example: data_de_mat_select <- data_de_mat[, !(colnames(data_de_mat) %in% c("SRR3583472", "SRR3583765", "SRR3584065", "SRR3584071"))]
# We did not remove any samples from this analysis
data_de_mat_select <- data_de_mat

# generating a traits file for remaining samples
traits <- data.frame(run = colnames(data_de_mat_select)) %>%
  left_join(metadata) %>%
  dplyr::select(run, developmental_stage, body_site, age, order)

write.table(traits, "results/traits.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

### 2.4 Variance stabilizing transformation ------------------------------------

# converting counts into DESeq2 object
#dds <- DESeqDataSetFromMatrix(countData = round(data_de_mat_select),
#                              colData = traits,
#                              design = body_site ~ developmental_stage)

dds <- DESeqDataSetFromMatrix(countData = round(data_de_mat_select),
                              colData = traits,
                              design = ~1)

# running VST
ddsDE <- DESeq(dds)
dds.vst <- assay(varianceStabilizingTransformation(ddsDE, blind = FALSE))

norm.counts <- t(dds.vst)
write.table(cbind(run = rownames(norm.counts), norm.counts), "results/normalized_counts.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

## 3. WGCNA analysis -----------------------------------------------------------

### 3.1 Selecting power parameter ----------------------------------------------

# trying different power parameters and evaluating network topology
powers <- c(c(1:10), seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(
  norm.counts,
  powerVector = powers,
  networkType = "signed",
  verbose = 5)

sft.data <- sft$fitIndices

# Visualizing topology metrics
pdf("plots/sft_power.pdf")
p1 <- ggplot(sft.data, aes(x = Power, y = SFT.R.sq, label = Power)) +
  geom_point() + 
  geom_text(nudge_y = 0.1) + 
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x = "Power", y = "Scale free topology model fit, unsigned R") + 
  theme_classic()

p2 <- ggplot(sft.data, aes(x = Power, y = mean.k., label = Power)) +
  geom_point() + 
  geom_text(nudge_y = 0.1) + 
  labs(x = "Power", y = "Mean connectivity") + 
  theme_classic()

grid.arrange(p1, p2, nrow = 2)
dev.off()

# observing the graph, we aimed for an R2>=0.8 and a minimum mean connectivity
soft_power <- 24

### 3.2 Network construction and module detection ------------------------------

# Setting default correlation function to WGCNA's function
temp_cor <- cor
cor <- WGCNA::cor

## Geschwind et al cut the tree using cutreeHybrid after generating the TOM
## Key parameters are corType="bicor", networkType="signed"

bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 20000, # depends on memory available
                          
                          # adjacency function options
                          power = soft_power,
                          networkType = "signed",
                          
                          # topological overlap options
                          TOMType = "signed", # TOM = topological overlap matrix
                          
                          # tree cut options
                          deepSplit = 4,
                          detectCutHeight = 0.995, 
                          minModuleSize = 30,
                          
                          # module merging options
                          mergeCutHeight = 0.25, # mergeCloseModules: Merging modules whose distance is less than value
                          
                          # other options
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

# Reassigning cor function to R
cor <- temp_cor

# Saving/loading network
#saveRDS(bwnet, file = "bwnet.rds")
bwnet <- readRDS(file = "bwnet.rds")

traits <- fread("results/traits.tsv") %>%
  mutate(developmental_stage = factor(developmental_stage, ordered = TRUE, 
                                      levels = c("Early prenatal", "Early mid-prenatal", 
                                                 "Late mid-prenatal", "Late prenatal", 
                                                 "Early infancy", "Late infancy", 
                                                 "Early childhood", "Late childhood", 
                                                 "Adolescence", "Adulthood"))) %>%
  arrange(developmental_stage)

norm.counts <- as.matrix(fread("results/normalized_counts.tsv") %>% column_to_rownames("run"))

# Counting the number of genes assigned to each module
table(bwnet$colors)
dim(table(bwnet$colors))
median(table(bwnet$colors))
mean(table(bwnet$colors))
max(table(bwnet$colors))
min(table(bwnet$colors))

as.data.frame(table(bwnet$colors)) %>%
  ggplot(aes(y = Var1, x = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = sort(unique(bwnet$colors))) +
  #scale_y_discrete(expand = c(0, 0)) +
  theme_classic() +
  xlab("No. of genes") +
  theme(legend.position="none",
        axis.title.y = element_blank())
ggsave("plots/module_gene_number.pdf", width = 4, height = 6)

# Saving gene assignment to a file
df.modules <- data.frame(module = bwnet$colors) %>% rownames_to_column("gene")
fwrite(df.modules, "results/gene_module_assignment.tsv", sep = "\t")

# Plotting the gene dendrogram and the assigned modules underneath
unmergedColors <- bwnet$unmergedColors[!is.na(bwnet$blocks)]
mergedColors <- bwnet$colors[!is.na(bwnet$blocks)]

pdf("plots/genes_dendrogram.pdf", height = 5, width = 8)
plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(unmergedColors, mergedColors),
                    c("unmerge","merge"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)
dev.off()

# Note: the grey module contains all genes that weren't clustered into a module

### 3.3 Module comparison with eigengenes --------------------------------------

# Obtaining module eigengenes from network and ordering per developmental stage
module_eigengenes <- bwnet$MEs
colnames(module_eigengenes) <- sub("^ME", "B-", colnames(module_eigengenes))
  
# Calculating the module's eigengenes dissimilarity (distance) and clustering
module_eigengenes.diss = 1-cor(module_eigengenes)
module_eigengenes.tree = hclust(as.dist(module_eigengenes.diss), method = "average")

pdf("plots/modules_corr_tree.pdf", width = 10, height = 5)
plot(module_eigengenes.tree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# abline(h=0.6, col = "red")
dev.off()

# Obtaining the co-expression heatmap between eigengenes
cormat <- round(cor(module_eigengenes),2)

pdf("plots/modules_corr_heatmap.pdf", width = 10, height = 5)
heatmap.2(cormat, density = "none", trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# Visualizing eigengenes - barplots
traits <- traits %>% arrange(order) %>% mutate(order = factor(order, level = unique(order)))
traits %>% dplyr::select(age, order) %>% distinct()

n1 <- length(unique(subset(traits, developmental_stage == "Early prenatal")$order))
n2 <- length(unique(subset(traits, developmental_stage == "Early mid-prenatal")$order))
n3 <- length(unique(subset(traits, developmental_stage == "Late mid-prenatal")$order))

blues  <- brewer.pal(3, "Blues")[1:n1]   # 8 distinct blues
greens <- brewer.pal(3, "Greens")[1:n2]  # 8 distinct greens
reds   <- brewer.pal(3, "Reds")[1:n3]    # 8 distinct reds

stage_colors <- c(blues, greens, reds)
names(stage_colors) <- unique(traits$order)
run_colors <- stage_colors[traits$order]

colors <- sort(gsub("B-","", names(module_eigengenes)))
module_eigengenes %>%
  rownames_to_column("run") %>%
  mutate(run = factor(run, ordered = T, levels = run)) %>%
  pivot_longer(!run, names_to = "module", values_to = "value") %>%
  ggplot(aes(x = run, y = value, fill = module)) + 
    geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
    scale_fill_manual(values = colors) +
    facet_wrap(~module, ncol = 9) + # , scales = "free"
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = run_colors)) +
    guides(fill = "none")
ggsave("plots/modules_eigengenes.pdf", height = 10, width = 20)

# Visualizing eigengenes - heatmap v1
module_eigengenes %>%
  rownames_to_column("run") %>%
  mutate(run = factor(run, ordered = T, levels = run)) %>%
  #mutate(run = factor(run, ordered = T, levels = order)) %>%
  pivot_longer(-run, names_to = "name", values_to = "value") %>%
  mutate(name = gsub("ME", "", name)) %>%
  ggplot(., aes(x=run, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module Eigengene \"expression\"", y = "Modules")
ggsave("plots/modules_eigengenes_heatmap1.pdf", height = 10, width = 10)

# Visualizing eigengenes - heatmap v2
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100))
traitColors <- labels2colors(traits$developmental_stage, 
                             colorSeq = c("#abdda4", "#fee08b", "#d53e4f", "#9e0142", "#e6f598",  "#f46d43"))

pdf('plots/modules_eigengenes_heatmap2.pdf', height = 10, width = 10)
heatmap.2(t(as.matrix(module_eigengenes)), 
          col = col, dendrogram = "row", Colv = F,
          ColSideColors = traitColors,
          trace = "none", density.info = "none")
dev.off()

### 3.4 Relating modules to condition ------------------------------------------

# Correlation with brain site and developmental stage

# Defining the number of genes and samples
nGenes = ncol(norm.counts)
nSamples = nrow(norm.counts)

# Converting traits to integers
traits <- traits %>% mutate(site_dev = paste0(developmental_stage, "-", body_site))

traits$developmental_stage <- factor(traits$developmental_stage)
binarized <- binarizeCategoricalColumns(traits$developmental_stage, 
                                        dropFirstLevelVsAll = FALSE,
                                        dropUninformative = FALSE,
                                        includePairwise = FALSE,
                                        includeLevelVsAll = TRUE,
                                        minCount = 1)
bin_traits <- cbind(traits, binarized) %>%
  column_to_rownames(var = "run") %>%
  dplyr::select(!c(developmental_stage, body_site, site_dev, order, age))

# Calculating correlation between eigengenes expression and traits
module.trait.corr <- cor(module_eigengenes, bin_traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

heatmap.data <- merge(module_eigengenes, bin_traits,  by = 'row.names', all = TRUE)
heatmap.data <- column_to_rownames(heatmap.data, var = "Row.names")

new_col_order <- c(sort(names(heatmap.data)[1:23]), c("data.Late mid-prenatal.vs.all", "data.Early mid-prenatal.vs.all", "data.Early prenatal.vs.all"))
heatmap.data.reordered <- heatmap.data[, new_col_order]

names(heatmap.data.reordered)
pdf("plots/module_devsite_relationship.pdf", width = 20, height = 3.5)
CorLevelPlot(heatmap.data.reordered,
             y = names(heatmap.data.reordered)[24:26],
             x = names(heatmap.data.reordered)[1:23],
             col = c("blue", "skyblue", "white", "pink", "red"),
             rotLabX = 45,
             cexCorval = 0.8,
             signifSymbols = c("*", ""),
             signifCutpoints = c(0, 0.05, 1))
dev.off()

### 3.5 Obtaining module membership and gene significance ----------------------

# Assigning genes on each module to variable

module_genes <- data.frame(module = bwnet$colors) %>% rownames_to_column(var = "gene")

# We can obtain a measure of how all genes are associated with a module by
# calculating their correlation to the module eigengene. This is known as
# module membership or intramodular connectivity.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pval <- corPvalueStudent(module.membership.measure, nSamples)
write.table(module.membership.measure, "results/module.membership.measure.tsv", 
            sep = "\t", row.names = T, quote = F, col.names = NA)

as.data.frame(module.membership.measure) %>%
  rownames_to_column("module") %>% 
  pivot_longer(!module, names_to = "gene", values_to = "membership") %>% 
  fwrite("results/module.membership.measure.long.tsv", sep = "\t")

# We can also obtain the association of a gene to a certain trait calculating
# gene significance, which is the correlation between expression data and 
# conditions.

gene.signif.corr <- cor(norm.counts, bin_traits, use = 'p')
gene.signif.corr.pvals <- corPvalueStudent(gene.signif.corr, nSamples)
fwrite(gene.signif.corr, "results/gene.signif.corr.tsv", sep = "\t", row.names = TRUE)

# Distribution of module membership across modules

df.module.1 <- 
  module.membership.measure %>%
  t() %>% as.data.frame() %>% rownames_to_column(var = "gene") %>%
  left_join(df.modules, by = "gene") %>% 
  pivot_longer(!c(gene, module), names_to = "MM_module", values_to = "MM_value") %>% 
  filter(paste0("B-", module) == MM_module) %>%
  dplyr::select(gene, MM_module, MM_value) %>%
  rename("assigned_module" = MM_module, "assigned_value" = MM_value)

df.module.1 %>% 
  filter(assigned_module %in% c("B-turquoise", "B-blue", "B-brown", "B-yellow", "B-grey")) %>%
  ggplot(aes(assigned_value, color = assigned_module)) + 
  geom_density() +
  scale_color_manual(values = sort(c("turquoise", "blue", "brown", "yellow", "grey"))) +
  theme_classic()
ggsave("plots/module_membership_distribution_selected.pdf", width = 10, height = 5)

df.module.1 %>% 
  group_by(assigned_module) %>%
  summarize(median(assigned_value)) %>%
  as.data.frame()

df.module.2 <- module.membership.measure %>%
  t() %>% as.data.frame() %>% rownames_to_column(var = "gene") %>% 
  pivot_longer(!c(gene), names_to = "MM_module", values_to = "MM_value") %>% 
  group_by(gene) %>%
  mutate(MM_max = max(MM_value)) %>%
  filter(MM_max == MM_value) %>%
  dplyr::select(gene, MM_module, MM_value) %>%
  rename("max_module" = MM_module, "max_value" = MM_value)

df.module.3 <- module.membership.measure %>%
  t() %>% as.data.frame() %>% rownames_to_column(var = "gene") %>% 
  pivot_longer(!c(gene), names_to = "MM_module", values_to = "MM_value") %>% 
  group_by(gene) %>%
  mutate(MM_min = min(MM_value)) %>%
  filter(MM_min == MM_value) %>%
  dplyr::select(gene, MM_module, MM_value) %>%
  rename("min_module" = MM_module, "min_value" = MM_value)

df.module.4 <- full_join(df.module.1, df.module.2, by = "gene") %>% full_join(df.module.3, by = "gene")
fwrite(df.module.4, "results/module.membership.results.tsv", sep = "\t")

## 4. Enrichment analyses ------------------------------------------------------

### 4.1 GO enrichment analysis -------------------------------------------------

modules_colors <- names(table(bwnet$colors))

for(module in modules_colors){
  print(module)
  
  selected_color <- module
  selected_eigengene <- paste0("ME", module)
  
  geneid2ensembl <- fread("annotations/geneid2ensembl.tsv", sep = "\t", header = FALSE,
                          col.names = c("ensembl", "geneid"))
  
  module_interest.genes <- names(bwnet$colors[bwnet$colors == selected_color])
  module_interest.genes.ensembl <- geneid2ensembl[geneid2ensembl$geneid %in% module_interest.genes,]$ensembl
  
  module_interest.ego <- enrichGO(gene          = module_interest.genes.ensembl,
                                  OrgDb         = org.Hs.eg.db,
                                  keyType       = 'ENSEMBL',
                                  ont           = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.1)
  
  if(!is.null(module_interest.ego)){
    significant_gos <- length(module_interest.ego@result$p.adjust[module_interest.ego@result$p.adjust<0.05])
    
    if(significant_gos > 0){
      i <- dotplot(module_interest.ego, showCategory = 10, title = module) +
        scale_color_gradient(low = "blue", 
                             high = "grey",
                             limits = c(0, 0.05))
      
      pdf(paste0("plots/ego_", selected_color, ".pdf"), height = 4, width = 6)
      print(i)
      dev.off()
    }
  }
  
}

### 4.2 ASD-related genes enrichment -------------------------------------------

mssng <- fread("annotations/MSSNG_genes.tsv", sep = "\t", header = FALSE, col.names = c("symbol", "gene_id", "dataset"))
module_genes.mssng <- left_join(module_genes, mssng, by=c("gene" = "gene_id"))
module_genes.mssng.total <- module_genes.mssng %>% group_by(module) %>% summarize(total = n())
module_genes.mssng.mssng <- module_genes.mssng %>% drop_na() %>% group_by(module) %>% summarize(mssng = n())

module_genes.mssng.counts <- 
  Reduce(function(x, y) merge(x, y, all=TRUE, by = "module"),
         list(module_genes.mssng.total, module_genes.mssng.mssng)) %>%
  mutate(mssng_perc = mssng/total*100)

nTotalGenes <- sum(module_genes.mssng.counts$total)
nMssngGenes <- sum(module_genes.mssng.counts$mssng, na.rm=TRUE)

hypergeometric.pvals <- c()
for (row in 1:nrow(module_genes.mssng.counts)) { 
  total_module <- module_genes.mssng.counts[row, 2]
  mssng_module <- module_genes.mssng.counts[row, 3]
  hypergeometric.pvals[row] <- phyper(mssng_module, 
                                      nMssngGenes, 
                                      nTotalGenes-nMssngGenes, 
                                      total_module, 
                                      lower.tail= FALSE)
}
module_genes.mssng.counts$pvals <- round(hypergeometric.pvals, 5)
module_genes.mssng.counts %>% filter(pvals <= 0.05)
module_genes.mssng.counts %>% fwrite("results/ASD_enrichment.tsv", sep = "\t")

### 4.3 SD-98 genes enrichment -------------------------------------------------

# Identifying modules enriched in SD-98 genes

sd98 <- fread("annotations/SD98_genes.tsv", sep = "\t", header = FALSE,
              col.names = c("name", "ensembl", "gene_id", "biotype"))
rownames(data_de_mat)[rownames(data_de_mat) %in% sd98$gene_id] %>% length() # 448

module_genes.annot <- left_join(module_genes, sd98, by=c("gene" = "gene_id"))
module_genes.annot.total <- module_genes.annot %>% group_by(module) %>% summarize(total = n())
module_genes.annot.sd98 <- module_genes.annot %>% drop_na() %>% group_by(module) %>% summarize(sd98 = n())
module_genes.annot.sd98.protein <-  module_genes.annot %>% drop_na() %>% 
  filter(biotype == "protein_coding") %>% group_by(module) %>% summarize(sd98_protein = n())

module_genes.annot.counts <- 
  Reduce(function(x, y) merge(x, y, all=TRUE, by = "module"),
         list(module_genes.annot.total, module_genes.annot.sd98, module_genes.annot.sd98.protein)) %>%
  mutate(sd98_perc = sd98/total*100)

nTotalGenes <- sum(module_genes.annot.counts$total)
nSD98Genes <- sum(module_genes.annot.counts$sd98, na.rm=TRUE)

hypergeometric.pvals <- c()
for (row in 1:nrow(module_genes.annot.counts)) { 
  total_module <- module_genes.annot.counts[row, 2]
  sd98_module <- module_genes.annot.counts[row, 3]
  hypergeometric.pvals[row] <- phyper(sd98_module, nSD98Genes, nTotalGenes-nSD98Genes, total_module, lower.tail= FALSE)
}
module_genes.annot.counts$pvals <- round(hypergeometric.pvals, 5)
module_genes.annot.counts %>% filter(pvals <= 0.05)
module_genes.annot.counts %>% fwrite("results/SD98_enrichment.tsv", sep = "\t")

### 4.4 Human-specific duplicated genes ----------------------------------------

humdups <- fread("annotations/human_specific_sd98_genes.tsv", sep = "\t", header = FALSE,
              col.names = c("name", "gene_id", "type", "cn"))
rownames(data_de_mat)[rownames(data_de_mat) %in% humdups$gene_id] %>% length() # 226

module_genes.annot <- left_join(module_genes, humdups, by=c("gene" = "gene_id"))
module_genes.annot.total <- module_genes.annot %>% group_by(module) %>% summarize(total = n())
module_genes.annot.group <- module_genes.annot %>% drop_na() %>% group_by(module) %>% summarize(group = n())

module_genes.annot.counts <- 
  Reduce(function(x, y) merge(x, y, all=TRUE, by = "module"),
         list(module_genes.annot.total, module_genes.annot.group)) %>%
  mutate(group_perc = group/total*100)

nTotalGenes <- sum(module_genes.annot.counts$total) # 17388
nGroupGenes <- sum(module_genes.annot.counts$group, na.rm=TRUE) # 226

hypergeometric.pvals <- c()
for (row in 1:nrow(module_genes.annot.counts)) { 
  total_module <- module_genes.annot.counts[row, 2]
  group_module <- module_genes.annot.counts[row, 3]
  hypergeometric.pvals[row] <- phyper(group_module, nGroupGenes, nTotalGenes-nGroupGenes, total_module, lower.tail= FALSE)
}
module_genes.annot.counts$pvals <- round(hypergeometric.pvals, 5)
module_genes.annot.counts %>% filter(pvals <= 0.05)
module_genes.annot.counts %>% fwrite("results/HumDups_enrichment.tsv", sep = "\t")

### 4.5 Human-specific CN-constrianed duplicated genes -------------------------

humdups <- fread("annotations/human_specific_sd98_genes.tsv", sep = "\t", header = FALSE,
                 col.names = c("name", "gene_id", "type", "cn")) %>%
  filter(cn != "Polymorphic")

module_genes.annot <- left_join(module_genes, humdups, by=c("gene" = "gene_id"))
module_genes.annot.total <- module_genes.annot %>% group_by(module) %>% summarize(total = n())
module_genes.annot.group <- module_genes.annot %>% drop_na() %>% group_by(module) %>% summarize(group = n())

module_genes.annot.counts <- 
  Reduce(function(x, y) merge(x, y, all=TRUE, by = "module"),
         list(module_genes.annot.total, module_genes.annot.group)) %>%
  mutate(group_perc = group/total*100)

nTotalGenes <- sum(module_genes.annot.counts$total) # 17388
nGroupGenes <- sum(module_genes.annot.counts$group, na.rm=TRUE) # 226

hypergeometric.pvals <- c()
for (row in 1:nrow(module_genes.annot.counts)) { 
  total_module <- module_genes.annot.counts[row, 2]
  group_module <- module_genes.annot.counts[row, 3]
  hypergeometric.pvals[row] <- phyper(group_module, nGroupGenes, nTotalGenes-nGroupGenes, total_module, lower.tail= FALSE)
}
module_genes.annot.counts$pvals <- round(hypergeometric.pvals, 5)
module_genes.annot.counts %>% filter(pvals <= 0.05)
module_genes.annot.counts %>% fwrite("results/HumDups_fixed_enrichment.tsv", sep = "\t")

### 4.6 Cell markers enrichment ------------------------------------------------

humdups <- fread("annotations/human_specific_sd98_genes.tsv", sep = "\t", header = FALSE,
                 col.names = c("name", "gene_id", "type", "cn")) %>%
  filter(cn != "Polymorphic")

module_genes.markers <- left_join(module_genes, cell_makers, by=c("gene" = "MarkerID"))
module_genes.markers.total <- module_genes.markers %>% group_by(module) %>% summarize(total = n())
module_genes.markers.markerTotal <- module_genes.markers %>% group_by(`Cell-type`) %>% summarize(markersTotal = n()) %>% drop_na()
module_genes.markers.markersModule <- module_genes.markers %>% drop_na() %>% group_by(module, `Cell-type`) %>% summarize(markers = n()) 

module_genes.markers.counts <- 
  left_join(module_genes.markers.markersModule, module_genes.markers.markerTotal, by = "Cell-type") %>%
  left_join(module_genes.markers.total, by = "module")
  
nTotalGenes <- length(module_genes.markers$gene)

hypergeometric.pvals <- c()
for (row in 1:nrow(module_genes.markers.counts)) { 
  total_module <- module_genes.markers.counts[[row, 5]]
  markers_module <- module_genes.markers.counts[[row, 3]]
  nMarkersGenes <- module_genes.markers.counts[[row, 4]]
  
  hypergeometric.pvals[row] <- phyper(markers_module, 
                                      nMarkersGenes, 
                                      nTotalGenes-nMarkersGenes, 
                                      total_module, 
                                      lower.tail= FALSE)
}
module_genes.markers.counts$pvals <- round(hypergeometric.pvals, 5)
module_genes.markers.counts %>% filter(pvals <= 0.05)
module_genes.markers.counts %>% fwrite("results/CelMarker_enrichment.tsv", sep = "\t")

### 4.7 Cell markers enrichment PanglaoDB --------------------------------------

cell_makers <- fread("annotations/PanglaoDB.tsv", sep = "\t", header = TRUE) %>% separate_rows(cell_type, sep = ";")

module_genes.markers <- left_join(module_genes, cell_makers, by=c("gene" = "gene_id"))
module_genes.markers.total <- module_genes.markers %>% group_by(module) %>% summarize(total = n())
module_genes.markers.markerTotal <- module_genes.markers %>% group_by(cell_type) %>% summarize(markersTotal = n()) %>% drop_na()
module_genes.markers.markersModule <- module_genes.markers %>% drop_na() %>% group_by(module, cell_type) %>% summarize(markers = n()) 

module_genes.markers.counts <- 
  left_join(module_genes.markers.markersModule, module_genes.markers.markerTotal, by = "cell_type") %>%
  left_join(module_genes.markers.total, by = "module")

nTotalGenes <- length(module_genes.markers$gene)

hypergeometric.pvals <- c()
for (row in 1:nrow(module_genes.markers.counts)) { 
  total_module <- module_genes.markers.counts[[row, 5]]
  markers_module <- module_genes.markers.counts[[row, 3]]
  nMarkersGenes <- module_genes.markers.counts[[row, 4]]
  
  hypergeometric.pvals[row] <- phyper(markers_module, 
                                      nMarkersGenes, 
                                      nTotalGenes-nMarkersGenes, 
                                      total_module, 
                                      lower.tail= FALSE)
}
module_genes.markers.counts$pvals <- round(hypergeometric.pvals, 5)
module_genes.markers.counts %>% filter(pvals <= 0.05)
module_genes.markers.counts %>% fwrite("results/CelMarkerPanglaoDB_enrichment.tsv", sep = "\t")
