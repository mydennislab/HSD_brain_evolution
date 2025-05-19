# R version: 4.2.2
# BiocManager version: 3.16

# Data handling
#install.packages("tidyverse")
#install.packages("data.table")
library(tidyverse)
library(data.table)

# Visualization
# install.packages("gridExtra")
# install.packages("gplots")
# install.packages("devtools")
# library(devtools)
# devtools::install_github("kevinblighe/CorLevelPlot")
# install.packages("RColorBrewer")
library(gridExtra)
library(gplots)
library(CorLevelPlot)
library(RColorBrewer)

# Network analysis
# install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz", repos=NULL, type="source")
# BiocManager::install("DESeq2")
# BiocManager::install("impute")
# BiocManager::install("preprocessCore")
# install.packages("WGCNA")
# install.packages("psych")
# install.packages("igraph")
library(DESeq2)
library(WGCNA)
library(psych)
library(igraph)

# Gene ontology over-representation
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(GO.db) # GO.db_3.16.0
library("org.Hs.eg.db", character.only = TRUE)

# BiocManager::install("limma")
# BiocManager::install("edgeR")
# install.packages("venneuler")
# install.packages("qvalue")
# devtools::install_github("jtlovell/limmaDE2")
library(limmaDE2)

sessionInfo()

# WGCNA expression analysis of SD-98 genes in CORTECON dataset -----------------

# Resources used in this analysis:
# - https://www.youtube.com/watch?v=gYE59uEMXT4&ab_channel=Bioinformagician
# - https://www.youtube.com/watch?v=mzXIxjPr_Mc&t=514s&ab_channel=Bioinformagician
# - https://bioinformaticsworkbook.org/tutorials/wgcna.html
# - https://rstudio-pubs-static.s3.amazonaws.com/687551_ed469310d8ea4652991a2e850b0018de.html
# 
# The CORTECON dataset (https://pubmed.ncbi.nlm.nih.gov/24991954/) includes 
# transcriptomic data from human cerebral cortex at different developmental time points.
# 
# In this analysis, we sought to identify genes co-expressed during human brain 
# development building an unsupervised co-expression network with the WGCNA 
# package. In particular, we aimed to answer the following questions:
# 1. Are there any modules with enrichment of duplicated genes? If so, what are
#    GO terms for those modules
# 2. Can we pull out a module that has clear connections with proliferation or 
#    autism (based on known genes) and then generate a nice network image of 
#    this/these, highlighting duplicated genes that would be worth following up?

# The overall pipeline involves the following steps:
# - Quantification of expression using Salmon/tximport. 
# - Importantly, to perform between sample normalization we need to use the raw 
#   counts instead of TPM.
# - Between sample normalization with DESeq2's variance stabilizing transformation.
# - Building network using WGCNA.
# - Identification of modules enriched in duplicated genes and GO 
#   over-representation ("guilt by association")
# - Using markers to identify modules associated with neural proliferation and 
#   autism, and look for SD-98 genes in these modules. 
# - Generate network plot.

## 1. Loading files ------------------------------------------------------------

# Reading counts
data <- read.table(file = 'counts/allSamples_vandeLeemput_counts.tsv', sep = '\t', header = TRUE)
data_pivot <- data %>% 
  pivot_longer(-gene, names_to = "run", values_to = "counts") %>%
  mutate(run = str_remove(run, "_1"))

# Reading metadata
metadata <- read.table(file = 'metadata/CORTECON_metadata.tsv', sep = '\t', header = TRUE) %>%
  dplyr::select(Run, RUN_ID, Stage) %>% 
  rename("run" = Run, "dataset" = RUN_ID, "condition" = Stage)
data_annot <- left_join(data_pivot, metadata, by = c("run" = "run"))

metadata %>%
  mutate(condition = factor(condition, 
                            levels = c("Pluripotency", 
                                       "Neural Differentation", 
                                       "Cortical Specification", 
                                       "Deep Layer Formation", 
                                       "Upper Layer Formation"))) %>%
  filter(dataset == "Cortecon") %>%
  ggplot(aes(y=condition)) + 
  geom_bar(width = 0.75, fill = c("#0fa09d", "#cfb277", "#e45025", "#6d8646", "#c18748")) +
  theme_classic()
ggsave("plots/cortecon_metadata.pdf", height = 3, width = 4)

# Filtering for main CORTECON
data_cortecon <- data_annot %>% filter(dataset == "Cortecon")

## 2. Data pre-processing ------------------------------------------------------

# According to WGCNA's FAQ website, (1) it is recommended to remove features 
# whose counts are consistently low and (2) counts should be normalize using 
# variance-stabilizing transformation. The authors indicate that as long as the
# samples were processed in the same way, the use of FPKM, RPKM, or TPM will not
# impact the network. They also recommend to use ComBat for batch effect removal.
# Finally, they usually check quantile scatterplots to make sure there are no 
# systematic shifts between samples. 
# Source: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html

### 2.1 Removing outliers ------------------------------------------------------

# According to a thread inthe Bioconductor forum, outliers samples indicate strong
# drivers of expressions, which will results in very large modules. Thus, it is
# recommended to remove them.

# transforming data to a matrix
data_de <- data_cortecon %>%
  dplyr::select(gene, run, counts) %>%
  pivot_wider(names_from = run, values_from = counts) %>%
  column_to_rownames(var = "gene")

# converting the first column to row names
data_de_mat <- as.matrix(data_de)
dim(data_de)

# First, we used the goodSamplesGenes function. The documentation states:
# > This function checks data for missing entries, entries with weights below a
# threshold, and zero-variance genes, and returns a list of samples and genes 
# that pass criteria on maximum number of missing or low weight values. If 
# necessary, the filtering is iterated.

gsg <- goodSamplesGenes(t(data_de_mat))
summary(gsg)
gsg$allOK
data_de_mat <- data_de_mat[gsg$goodGenes == TRUE,]
dim(data_de_mat)

# We observed that while all samples passed this test, around ~20,000 genes did
# not, so we removed them from the data set, going from 62216 to 47160 genes.

# Then, we identified outliers with hierarchical clustering and PCA.

# Hierarchical clustering
htree <- hclust(dist(t(data_de_mat)), method = "average")

data_de_mat.meta <- data.frame(run = colnames(data_de_mat)) %>%
  left_join(metadata) %>%
  dplyr::select(run, condition)

traitColors <- labels2colors(sort(data_de_mat.meta$condition), 
                             colorSeq = c("#0fa09d", "#cfb277", "#e45025", "#6d8646", "#c18748"))

pdf("plots/samples_dendrogram.pdf", width = 15, height = 10)
plotDendroAndColors(htree, traitColors,
                    groupLabels = names(htree), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# PCA
pca <- prcomp(t(data_de_mat))
pca.data <- as.data.frame(pca$x)
pca.var <- pca$sdev
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

ggplot(pca.data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = factor(data_de_mat.meta$condition, 
                                levels = c("Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation"))
                 ), 
             size = 2) +
  geom_text(label = rownames(pca.data), size = 3) +
  scale_color_manual(values = c("#0fa09d", "#cfb277", "#e45025", "#6d8646", "#c18748")) +
  labs(x = paste0('PC1:', pca.var.percent[1], '%'),
       y = paste0('PC2:', pca.var.percent[2], '%')) +
  theme_bw() +
  labs(fill = "", color = "")
ggsave("plots/samples_pca.pdf", width = 6, height = 4)

# PCA and hierarchical clustering show that samples SRR1238515 and SRR1238516 
# associated with Upper Layer Formation, seem to be outliers. Therefore we opted
# for removing them from the network construction.

samples.to.be.excluded <- c('SRR1238515','SRR1238516')
data_de_mat.subset <- data_de_mat[, !(colnames(data_de_mat) %in% samples.to.be.excluded)]

dim(data_de_mat) # 23
dim(data_de_mat.subset) # 21

# Generating a metadata file associating every run with a condition *after* filtering outliers.
traits <- data.frame(run = colnames(data_de_mat.subset)) %>%
  left_join(metadata) %>%
  dplyr::select(run, condition)

write.table(traits, "results/traits.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

### 2.2 Variance stabilizing transformation ------------------------------------

# Then, we performed variance stabilizing transformation (VST) with DESeq2 
# package. The goal of variance-stabilizing transformation is to overcome the 
# effect of having different number of sequences in different samples. From the
# Bioconductor forum: "The variance stabilizing transformations are very 
# different from TPM and RPKM. These latter normalizations allow for comparison
# of values across genes, because they are proportional to original counts of 
# transcripts. However, you will see that they are not variance stabilizing. 
# Distances between samples will be highly weighted by contributions from gene 
# with highest TPM." Regarding blind parameter, in a Bioconductor forum, they
# suggest "blind=TRUE, when the transformation will be used to potentially 
# identify outliers or when doing QC. However, especially when there are many
# differences in counts which can be accounted for by terms in the design, 
# blind=FALSE is a better choice for visualization and downstream analysis. 
# And it's not so much of a concern to use blind=FALSE, as the only information 
# used from the steps involving the design is the experiment-wide trend of 
# dispersion over mean, not the individual gene-wise estimates. So, use the full
# design for the VST. In order to subtract effects of covariates, you can apply 
# limma's removeBatchEffect() to the transformed count matrix." Importantly, 
# "the VST transforms to a log2-like count which stabilizes the variance, and 
# additionally accounts for library size differences of the samples. The VST 
# does nothing to correct for different genes having different length."
# To use the VST function, we rounded the counts. This should not imply a 
# significant data loss, as counts are in other order of magnitude.

# Converting counts into DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(data_de_mat.subset),
                              colData = traits,
                              design = ~ condition)

# Considering WGCNA authors' recommendations, we removed features whose counts 
# are consistently low, meaning features with counts<10 in 90% samples.
dds90 <- dds[rowSums(counts(dds) >= 10) >= 19,]
nrow(dds90) 

# We ended up with 15,697 genes with over 10 counts in 90% of the samples.

# Running VST
dds90 <- DESeq(dds90)

# method 1
dds.vst.1 <- assay(varianceStabilizingTransformation(dds90)) # blind=TRUE is default here
# method 2
dds.vst.2 <- assay(varianceStabilizingTransformation(dds90, blind = FALSE))
# method 3
dds.vst.3 <- getVarianceStabilizedData(dds90) # bypasses blind=TRUE
# There is also wrapper function called vst()

# Methods 2 and 3 give the same results.
identical(dds.vst.1, dds.vst.3) # FALSE
identical(dds.vst.2, dds.vst.3) # TRUE

# We move forward with blind=FALSE option and transpose for downstream analyses
norm.counts <- t(dds.vst.2)
write.table(norm.counts, "results/normalized_counts.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

## 3. WGCNA analysis -----------------------------------------------------------

### 3.1 Selecting power parameter ----------------------------------------------

# The goal is to build a scale-free topology network, which is assumed to better 
# represent biological networks. A scale-free network is a network where the 
# distribution of node connections follows a power law distribution. This means 
# that most nodes have few connections with a few nodes acting as hubs with 
# multiple connections.

# To build a scale-free network from expression data, we need to apply a 
# correction factor. From WGCNA tutorial, "this factor enhances the differences 
# between strong and weak correlations. Weak values become closer to zero. 
# This power value must produce a graph similar to a scale-free network."

# We chose to perform a unsigned network to obtain any correlation between
# genes: "In signed networks, these highly connected hub genes may up-regulate 
# adjacent genes since they are positively correlated with them, while in 
# unsigned networks they may activate or repress their neighboring genes."

# Trying different power parameters and evaluating network topology
powers <- c(c(1:10), seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(
  norm.counts,
  powerVector = powers,
  networkType = "signed",
  verbose = 5)

sft.data <- sft$fitIndices

# Visualizing topology metrics
pdf("plots/soft_power.pdf")
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

# There are two metrics to decide the right power that generate a scale-free
# topology: SFT.R.sq. and mean connectivity (mean.k). We want to choose the 
# power that gives the max R square and the min mean connectivity.

# Observing the graph, we aimed for an R2>=0.8 and a minimum mean connectivity
soft_power <- 24

### 3.2 Network construction and module detection ------------------------------

# Setting default correlation function to WGCNA's function
temp_cor <- cor
cor <- WGCNA::cor

# To tune the right number of modules we constructed networks with different
# values of mergeCutHeight and deepSplit:
# - deepSplit ranges from 0 to 4, increase for more modules.
# - mergeCutHeight around 0.25-0.3 for large data sets. Reduce for more modules.

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
               mergeCutHeight = 0.15, # mergeCloseModules: Merging modules whose distance is less than value

               # other options
               numericLabels = FALSE,
               randomSeed = 1234,
               verbose = 3)

# Reassigning cor function to R
cor <- temp_cor

# Saving/loading network

# saveRDS(bwnet, file = "bwnet.rds")
bwnet <- readRDS(file = "bwnet.rds")
traits <- read.table("results/traits.tsv", sep = "\t", header = TRUE)
norm.counts <- read.table("results/normalized_counts.tsv", header = TRUE, sep = "\t")

# Counting the number of genes assigned to each module
table(bwnet$colors)
dim(table(bwnet$colors))
median(table(bwnet$colors))
mean(table(bwnet$colors))
min(table(bwnet$colors))
max(table(bwnet$colors))

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
pdf("plots/genes_dendrogram.pdf", height = 5, width = 8)
plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerge","merge"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)
dev.off()
# Note: the turquoise module contains all genes that weren't clustered into a module

### 3.3 Module comparison with eigengenes --------------------------------------

# The eigengene is the first PC of the module gene expression. It can be
# interpreted as the representative gene of the module.

# Obtaining module eigengenes from network
module_eigengenes <- bwnet$MEs
colnames(module_eigengenes) <- sub("^ME", "C-", colnames(module_eigengenes))

# Calculating the module's eigengenes dissimilarity (distance) and clustering
module_eigengenes.diss = 1-cor(module_eigengenes)
module_eigengenes.tree = hclust(as.dist(module_eigengenes.diss), method = "average")

pdf("plots/modules_corr_tree.pdf")
op <- par(mar = c(5, 4, 4, 8))
plot(as.dendrogram(module_eigengenes.tree),
     horiz = TRUE, xlab = "", sub = "", cex   = 0.5)
# abline(h=0.6, col = "red")
dev.off()

par(op)

# Obtaining the co-expression heatmap between eigengenes
cormat <- round(cor(module_eigengenes),2)

pdf("plots/modules_corr_heatmap.pdf")
heatmap.2(cormat, density = "none", trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# Visualizing eigengenes

# barplots
colors <- sort(gsub("C-","", names(module_eigengenes)))
module_eigengenes %>%
  rownames_to_column("run") %>%
  pivot_longer(!run, names_to = "module", values_to = "value") %>%
  ggplot(aes(x = run, y = value, fill = module)) + 
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = colors) +
    facet_wrap(~module, ncol = 8) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill = "none") 
ggsave("plots/modules_eigengenes.pdf", height = 15, width = 20)

# heatmap v1
module_eigengenes %>%
  rownames_to_column("run") %>%
  pivot_longer(-run, names_to = "name", values_to = "value") %>%
  mutate(name = gsub("ME", "", name)) %>%
  ggplot(., aes(x=run, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "#fde725",
    mid = "#21918c",
    high = "#440154",
    midpoint = 0) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module Eigengene \"expression\"", y = "Modules")
ggsave("plots/modules_eigengenes_heatmap1.pdf")

# heatmap v2
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100))
traitColors <- labels2colors(traits$condition,colorSeq = c("#e41a1c", "#984ea3", "#377eb8", "#4daf4a", "#ff7f00"))
pdf('plots/modules_eigengenes_heatmap2.pdf')
heatmap.2(t(as.matrix(module_eigengenes)), 
          col = col, dendrogram = "row", Colv = F,
          ColSideColors = traitColors,
          trace = "none", density.info = "none")
dev.off()

### 3.4 Relating modules to condition ------------------------------------------

# Defining the number of genes and samples
nGenes = ncol(norm.counts)
nSamples = nrow(norm.counts)

# Converting traits to integers
traits$condition <- factor(traits$condition, 
                           levels = c("Pluripotency", 
                                      "Neural Differentation",
                                      "Cortical Specification",
                                      "Deep Layer Formation",
                                      "Upper Layer Formation"))
binarized <- binarizeCategoricalColumns(traits$condition, 
                                        dropFirstLevelVsAll = FALSE,
                                        dropUninformative = FALSE,
                                        includePairwise = FALSE,
                                        includeLevelVsAll = TRUE,
                                        minCount = 1)
bin_traits <- cbind(traits, binarized) %>%
  column_to_rownames(var = "run") %>%
  dplyr::select(-condition)

# Calculating correlation between eigengenes expression and traits
module.trait.corr <- cor(module_eigengenes, bin_traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

heatmap.data <- merge(module_eigengenes, bin_traits,  by = 'row.names', all = TRUE)
heatmap.data <- column_to_rownames(heatmap.data, var = "Row.names")

new_col_order <- c(sort(names(heatmap.data)[1:37]), names(heatmap.data)[38:42])
heatmap.data.reordered <- heatmap.data[, new_col_order]

names(heatmap.data.reordered)
pdf("plots/module_trait_relationship.pdf", width = 20, height = 4)
CorLevelPlot(heatmap.data.reordered,
             y = names(heatmap.data.reordered)[38:42],
             x = names(heatmap.data.reordered)[1:37],
             col = c("blue", "skyblue", "white", "pink", "red"),
             rotLabX = 45,
             cexCorval = 0.8,
             signifSymbols = c("*", ""),
             signifCutpoints = c(0, 0.05, 1))
dev.off()

### 3.5 Obtaining module membership and gene significance ----------------------

# Assigning genes on each module to variable

module_genes <- data.frame(module = bwnet$colors) %>%
  rownames_to_column(var = "gene")

# We can obtain a measure of how all genes are associated with a module by
# calculating their correlation to the module eigengene. This is known as
# module membership or intramodular connectivity.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pval <- corPvalueStudent(module.membership.measure, nSamples)
write.table(module.membership.measure, "results/module.membership.measure.tsv", 
            sep = "\t", row.names = T, quote = F, col.names = NA)

module.membership.measure.long <- as.data.frame(module.membership.measure) %>%
  rownames_to_column("module") %>% 
  pivot_longer(!module, names_to = "gene", values_to = "membership")

fwrite(module.membership.measure.long, "results/module.membership.measure.long.tsv", sep = "\t")

# We can also obtain the association of a gene to a certain trait calculating
# gene significance, which is the correlation between expression data and 
# conditions.

gene.signif.corr <- cor(norm.counts, bin_traits, use = 'p')
gene.signif.corr.pvals <- corPvalueStudent(gene.signif.corr, nSamples)
write.table(gene.signif.corr.pvals, "results/gene.signif.corr.tsv", 
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Distribution of module membership across modules

df.modules.2 <- 
  module.membership.measure %>%
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "gene") %>%
  left_join(df.modules, by = "gene") %>%
  pivot_longer(!c(gene, module), names_to = "MM_module", values_to = "MM_value") %>%
  filter(paste0("C-", module) == MM_module) %>%
  dplyr::select(gene, MM_module, MM_value) %>%
  rename("assigned_module" = MM_module, "assigned_value" = MM_value)

df.modules.2 %>% 
  filter(assigned_module %in% c("C-turquoise", "C-blue", "C-green", "C-yellow", "C-red", "C-grey")) %>%
  ggplot(aes(assigned_value, color = assigned_module)) + 
  geom_density() +
  scale_color_manual(values = sort(c("turquoise", "blue", "green", "yellow", "red", "grey"))) +
  theme_classic()
ggsave("plots/module_membership_distribution_selected.pdf", width = 10, height = 5)

df.modules.2 %>% 
  group_by(assigned_module) %>%
  summarize(median(assigned_value)) %>%
  as.data.frame()

## 4. Enrichment analyses ------------------------------------------------------

### 4.1 GO enrichment analysis -------------------------------------------------

modules_colors <- names(table(bwnet$colors))
modules_gos <- data.frame()

for(module in modules_colors){
    selected_color <- module
    selected_eigengene <- paste0("C-", module)
    
    geneid2ensembl <- fread("metadata/geneid2ensembl.tsv", sep = "\t", header = FALSE,
                            col.names = c("ensembl", "geneid"))
    
    module_interest.genes <- names(bwnet$colors[bwnet$colors == selected_color])
    module_interest.genes.ensembl <- geneid2ensembl[geneid2ensembl$geneid %in% module_interest.genes,]$ensembl
    
    module_interest.ego <- enrichGO(gene          = module_interest.genes.ensembl,
                                    OrgDb         = org.Hs.eg.db,
                                    keyType       = 'ENSEMBL',
                                    ont           = "BP",
                                    pAdjustMethod = "fdr",
                                    pvalueCutoff  = 0.05,
                                    qvalueCutoff  = 0.05)
    
    #module_interest.ego.filter <- gofilter(module_interest.ego, level = 6)
    module_interest.ego.filter <- module_interest.ego
    
    modules_go <-
      module_interest.ego.filter@result %>% 
      filter(p.adjust < 0.05) %>%
      mutate(module = selected_color)
    modules_gos <- rbind(modules_gos, modules_go)
    
    significant_gos <- length(module_interest.ego.filter@result$p.adjust[module_interest.ego.filter@result$p.adjust<0.05])
    
    if(significant_gos > 0){
        options(enrichplot.colours = c("blue","lightgray")) 
        i <- dotplot(module_interest.ego.filter, showCategory = 10, size = 1) + 
             scale_color_gradient(limits = c(0, 0.05))
        
        pdf(paste0("plots/ego_", selected_color, ".pdf"), height = 4, width = 5)
        print(i)
        i + scale_fill_gradientn(name = "adjusted p-values",
                                 colours = c("blue", "white"),
                                 limits= c(0, 0.05), 
                                 breaks=c(0, 0.025, 0.05) )
        dev.off()
    }
}

fwrite(modules_gos, "results/module_gos.tsv", sep = "\t", row.names = TRUE)

### 4.2 Obtaining top hub in cluster -------------------------------------------

# We can use the WGCNA function chooseTopHubInEachModule to obtain the hub of
# each module, i.e., the gene with the highest connectivity.
hub_genes <- 
  chooseTopHubInEachModule(
  norm.counts,
  bwnet$colors,
  power = 24, 
  type = "signed")

# Alternatively, we can generate a network and obtain degree centrality for
# each gene in the network. We did this analysis for specific modules. See
# later sections.

# Note: Besides connectivity, "hub" genes are also found as module membership,
# meaning the correlation between the gene expression and a module's eigengene.
# Often a correlation |r| > 0.8 is used to define relevant genes within a module.

### 4.3 ASD-related genes in clusters ------------------------------------------

mssng <-
  fread("metadata/MSSNG_genes.tsv", sep = "\t", header = FALSE,
        col.names = c("symbol", "gene_id", "dataset"))

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
module_genes.annot.counts %>% fwrite("results/ASD_enrichment.tsv", sep = "\t")

### 4.4 SD-98 genes module analysis --------------------------------------------

# Identifying modules enriched in SD-98 genes

sd98 <- fread("metadata/SD98_genes.tsv", sep = "\t", header = FALSE,
              col.names = c("name", "ensembl", "gene_id", "biotype"))

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

# Plotting of number of SD-98 genes per module
module_genes.annot.counts %>%
  mutate(label = paste0(round(sd98_perc, 2), "%", "; ", round(pvals,2))) %>%
  mutate(total = -total) %>%
  dplyr::select(module, total, sd98, label) %>%
  pivot_longer(!c(module, label), values_to = "number", names_to = "type") %>%
  ggplot(aes(x = module, y = number, fill = type)) + 
    geom_bar(stat="identity", position="identity") +
    geom_text(aes(y = 800, label = label)) +
    scale_fill_manual(values = c("skyblue", "lightgrey")) +
    scale_y_continuous(limits = c(-4000, 1000)) +
    coord_flip() + theme_minimal()
ggsave("plots/sd98_per_module.pdf")

### 4.4 Human duplicated genes module analysis ---------------------------------

humdups <- fread("metadata/human_specific_sd98_genes.tsv", sep = "\t", header = FALSE,
                 col.names = c("name", "gene_id", "type", "cn"))

module_genes.annot <- left_join(module_genes, humdups, by=c("gene" = "gene_id"))
module_genes.annot.total <- module_genes.annot %>% group_by(module) %>% summarize(total = n())
module_genes.annot.group <- module_genes.annot %>% drop_na() %>% group_by(module) %>% summarize(group = n())

module_genes.annot.counts <- 
  Reduce(function(x, y) merge(x, y, all=TRUE, by = "module"),
         list(module_genes.annot.total, module_genes.annot.group)) %>%
  mutate(group_perc = group/total*100)

nTotalGenes <- sum(module_genes.annot.counts$total) # 15695
nGroupGenes <- sum(module_genes.annot.counts$group, na.rm=TRUE) # 200

hypergeometric.pvals <- c()
for (row in 1:nrow(module_genes.annot.counts)) { 
  total_module <- module_genes.annot.counts[row, 2]
  group_module <- module_genes.annot.counts[row, 3]
  hypergeometric.pvals[row] <- phyper(group_module, nGroupGenes, nTotalGenes-nGroupGenes, total_module, lower.tail= FALSE)
}
module_genes.annot.counts$pvals <- round(hypergeometric.pvals, 5)
module_genes.annot.counts %>% filter(pvals <= 0.05)
module_genes.annot.counts %>% fwrite("results/HumDups_enrichment.tsv", sep = "\t")

### 4.5 Human duplicated CN-constrained genes module analysis ------------------

humdups_fixed <- fread("metadata/human_specific_sd98_genes.tsv", sep = "\t", header = FALSE,
                 col.names = c("name", "gene_id", "type", "cn")) %>%
  filter(cn != "Polymorphic")

module_genes.annot <- left_join(module_genes, humdups_fixed, by=c("gene" = "gene_id"))
module_genes.annot.total <- module_genes.annot %>% group_by(module) %>% summarize(total = n())
module_genes.annot.group <- module_genes.annot %>% drop_na() %>% group_by(module) %>% summarize(group = n())

module_genes.annot.counts <- 
  Reduce(function(x, y) merge(x, y, all=TRUE, by = "module"),
         list(module_genes.annot.total, module_genes.annot.group)) %>%
  mutate(group_perc = group/total*100)

nTotalGenes <- sum(module_genes.annot.counts$total) # 15695
nGroupGenes <- sum(module_genes.annot.counts$group, na.rm=TRUE) # 200

hypergeometric.pvals <- c()
for (row in 1:nrow(module_genes.annot.counts)) { 
  total_module <- module_genes.annot.counts[row, 2]
  group_module <- module_genes.annot.counts[row, 3]
  hypergeometric.pvals[row] <- phyper(group_module, nGroupGenes, nTotalGenes-nGroupGenes, total_module, lower.tail= FALSE)
}
module_genes.annot.counts$pvals <- round(hypergeometric.pvals, 5)
module_genes.annot.counts %>% filter(pvals <= 0.05)
module_genes.annot.counts %>% fwrite("results/HumDups_fixed_enrichment.tsv", sep = "\t")

### 4.6 Cell marker enrichment -------------------------------------------------
cell_makers <- fread("metadata/cell_markers.tsv", sep = "\t", header = TRUE)

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

### 4.7 Genomic hotspots enrichment --------------------------------------------

hotspots <-
  fread("metadata/genomic_disorders_genes.tsv", sep = "\t", header = FALSE,
        col.names = c("symbol", "gene_id", "hotspot"))

#hotspots <- hotspots %>% filter(grepl("15q25.2", hotspot))

module_genes.hotspots <- 
  left_join(module_genes, hotspots, by=c("gene" = "gene_id")) %>%
  group_by(gene, module, symbol) %>%
  summarize(hostpots = paste(hotspot, collapse = ";"))
  
module_genes.hotspots.total <- module_genes.hotspots %>% group_by(module) %>% summarize(total = n())
module_genes.hotspots.hotspots <- module_genes.hotspots %>% drop_na() %>% group_by(module) %>% summarize(hotspot = n())

module_genes.hotspots.counts <- 
  Reduce(function(x, y) merge(x, y, all=TRUE, by = "module"),
         list(module_genes.hotspots.total, module_genes.hotspots.hotspots)) %>%
  mutate(hotspot_perc = hotspot/total*100)

nTotalGenes <- sum(module_genes.hotspots.counts$total)
nHotspotsGenes <- sum(module_genes.hotspots.counts$hotspot, na.rm=TRUE)

hypergeometric.pvals <- c()
for (row in 1:nrow(module_genes.hotspots.counts)) { 
  total_module <- module_genes.hotspots.counts[row, 2]
  hotspots_module <- module_genes.hotspots.counts[row, 3]
  hypergeometric.pvals[row] <- phyper(hotspots_module, 
                                      nHotspotsGenes, 
                                      nTotalGenes-nHotspotsGenes, 
                                      total_module, 
                                      lower.tail= FALSE)
}
module_genes.hotspots.counts$pvals <- round(hypergeometric.pvals, 5)
module_genes.hotspots.counts %>% filter(pvals <= 0.05)
module_genes.markers.counts %>% fwrite("results/Hotspots_enrichment.tsv", sep = "\t")

## 5. Visualizing modules of interest ------------------------------------------

selected_color <- "yellow"

### 5.0 Visualizing expression patterns ----------------------------------------

# Looking at the expression patterns of genes within modules
# Example for black module
selected_module_genes <- names(bwnet$colors[bwnet$colors == selected_color])
selected_module_expression <- norm.counts[,colnames(norm.counts) %in% selected_module_genes]
pdf(paste0(selected_color, "_heatmap.pdf"))
heatmap.2(t(selected_module_expression), trace = "none", density = "none", 
          col = bluered(20), scale = "row", dendrogram = "row", Colv = FALSE)
dev.off()

### 5.1 Network using adjacency ------------------------------------------------

# select genes with high module membership
module_interest.genes <- names(bwnet$colors[bwnet$colors == selected_color])

module_interest.genes.highMM <-
  module.membership.measure.long %>%
  filter(gene %in% module_interest.genes) %>%
  filter(module == paste0("ME",selected_color)) %>%
  filter(membership > 0.5) %>%
  pull(gene)

# select genes within relevant categories, including:
# 1. Human duplicated genes
# 2. SD98 genes
# 3. MSSNG genes
# 4. Genomic disorders

human_specific <- fread("metadata/human_specific_sd98_genes.tsv", sep = "\t", header = FALSE,
                        col.names = c("symbol", "gene_id", "status"))
sd98 <- fread("metadata/SD98_genes.tsv", sep = "\t", header = FALSE,
              col.names = c("name", "ensembl", "gene_id", "biotype"))
mssng <- fread("metadata/MSSNG_genes.tsv", sep = "\t", header = FALSE, 
               col.names = c("symbol", "gene_id", "dataset"))
hotspots <- fread("metadata/genomic_disorders_genes.tsv", sep = "\t", header = FALSE, 
                  col.names = c("symbol", "gene_id", "hotspot"))

module_interest.genes.highMM <-
  module_interest.genes.highMM[module_interest.genes.highMM %in% c(sd98$gene_id, 
                                                                   mssng$gene_id)]

module_interest.expression <- norm.counts[,colnames(norm.counts) %in% module_interest.genes.highMM]

# define soft threshold for adjacency matrix
powers <- c(c(1:10), seq(from = 12, to=70, by=2))
sft <- pickSoftThreshold(module_interest.expression, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red")
abline(h=0.8, col="blue", lty=2)

adjacency <- adjacency(module_interest.expression, power=18, type="signed")
#adjacency[adjacency < 0] = 0
#adjacency[adjacency > 1] = 1
TOM <- TOMsimilarity(adjacency, TOMType="signed") # number of shared neighbors
row.names(TOM) = colnames(module_interest.expression)
colnames(TOM) = colnames(module_interest.expression)
adj <- TOM
#adj[adj > 0.1] = 1
#adj[adj != 1] = 0
network <- graph_from_adjacency_matrix(adj, weighted=TRUE, mode="undirected")
network <- igraph::simplify(network, remove.multiple=TRUE, remove.loops=TRUE)  # removes self-loops

max(E(network)$weight)
min(E(network)$weight)
network <- delete_edges(network, E(network)[which(E(network)$weight<0.1)]) # remove edges below pearson r
network <- delete_vertices(network, degree(network)==0) # remove unconnected nodes

vcol <- V(network)$name
vcol[vcol %in% human_specific$gene_id] <- "#f7423e"
vcol[vcol %in% sd98$gene_id] <- "#eab5b5"
vcol[vcol %in% mssng$gene_id] <- "#fcd262"
#vcol[vcol %in% hotspots$gene_id] <- "white"
vcol[grepl("CHM13", vcol)] <- selected_color
vcol[grepl("LOFF", vcol)] <- selected_color
V(network)$color <- vcol

vname <- V(network)$name
vname1 <- sd98[match(V(network)$name, sd98$gene_id), name]
vname2 <- mssng[match(V(network)$name, mssng$gene_id), symbol]
#vname3 <- hotspots[match(V(network)$name, hotspots$gene_id), symbol]
vname <- pmin(vname1, vname2, na.rm = TRUE)

V(network)$frame.color <- "black"

vwidth <- V(network)$name
vwidth[vwidth %in% hotspots$gene_id] <- 1.5
vwidth[grepl("CHM13", vwidth)] <- 0.5
vwidth[grepl("LOFF", vwidth)] <- 0.5
V(network)$frame.width <- vwidth

l <- layout_with_fr(network, weights = E(network)$weight)
l <- norm_coords(l, ymin=-0.6, ymax=0.6, xmin=-1, xmax=1)

pdf(paste0(selected_color, ".pdf"), width=12, height=12)
plot(network,
     layout=l,
     rescale=FALSE,
     edge.width=E(network)$weight,
     vertex.label.color="black",
     vertex.label.cex=0.8, 
     vertex.label.dist=0.1,
     vertex.label.family="Helvetica", 
     vertex.label.color="black", 
     vertex.label=vname, 
     vertex.size=4*log10(degree(network, mode="all")),
     edge.color="grey")
dev.off()

### 5.2 Network using Cytoscape RCy3 -------------------------------------------

library(RCy3) # & open Cytoscape
cytoscapePing() # connect with Cytoscape
cytoscapeVersionInfo() # check if Cytoscape connection worked

# create network in Cytoscape
createNetworkFromIgraph(network)

as_edgelist(network, names=T)
as_adjacency_matrix(network, attr="weight")
as_data_frame(network, what="edges")         # data frame describing edges
as_data_frame(network, what="vertices")      # data frame describing nodes

### 5.3 Network using wgcna2igraph ---------------------------------------------
simplify <- igraph::simplify

# does not preserve weights
net <- wgcna2igraph(bwnet, 
                    datExpr = norm.counts,
                    kME.threshold = 0.8, 
                    adjacency.threshold = 0.1,
                    adj.power = 24,
                    modules2plot = c("saddlebrown"), 
                    colors2plot = c("saddlebrown"))

net <- simplify(net, remove.multiple = F, remove.loops = T)
deg <- degree(net)


pdf("saddlebrown2.pdf", width = 20, height = 10)
plot(net, 
     layout=layout_with_fr,
     vertex.color="orange",
     vertex.label.family="Arial", 
     vertex.label.color="black", 
     vertex.frame.color="#555555",
     #vertex.label=V(graph), 
     vertex.label.color="black",
     vertex.label.cex=0.5,
     #vertex.size=2,
     vertex.size=log10(deg),
     edge.arrow.size=0,
     edge.width=1,
     edge.color="grey")
dev.off()

graph_hub_score <- hub_score(graph)$vector %>% sort(decreasing = T)

### 5.4 Using TOM --------------------------------------------------------
soft_power <- 18
TOM <- TOMsimilarityFromExpr(module_interest.expression,
                             power = soft_power)

row.names(TOM) = colnames(module_interest.expression)
colnames(TOM) = colnames(module_interest.expression)

edge_list <- data.frame(TOM) %>%
  mutate(gene1 = row.names(.)) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, weight = value) %>%
  unique() %>%
  subset(!(gene1==gene2))

# Network with igraph
module_interest.network <- graph_from_data_frame(d = edge_list, directed = F)
edge.attributes(module_interest.network)$weight

# Simplify the adjacency object
module_interest.network <- igraph::simplify(module_interest.network, remove.multiple=TRUE, remove.loops=TRUE)

# Remove edges below absolute Pearson correlation t
module_interest.network <- delete_edges(
  module_interest.network, 
  E(module_interest.network)[which(E(module_interest.network)$weight<0.1)])

# Remove any vertices remaining that have no edges
module_interest.network <- delete_vertices(module_interest.network, degree(module_interest.network)==0)

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
#module_interest.network <- igraph::mst(module_interest.network, algorithm="prim") 

# Network in igraph
deg <- degree(module_interest.network, mode="all")
plot(module_interest.network, layout = layout.fruchterman.reingold, edge.arrow.mode=0, 
     vertex.label.cex = 0.5, vertex.size = log10(deg))  # vertex.label = NA

# Network for Cytoscape
fwrite(as_data_frame(module_interest.network), paste0(selected_color, "_network_filter.tsv"), sep = "\t")

## 5. Other --------------------------------------------------------------------

### 5.1 Data for Gerald --------------------------------------------------------

# SRGAP2, HYDIN, ARHGAP11, PTPN20, NCF1, GPR89, PDZK1, NPY4R, FAM72, FRMPD2
pHSDs <- sd98 %>% 
  filter(name %in% c("SRGAP2", "SRGAP2B", "SRGAP2C", "HYDIN", "HYDIN2", 
                     "ARHGAP11A",  "ARHGAP11B", "PTPN20", "PTPN20CP", "NCF1", 
                     "NCF1B", "NCF1C", "GPR89A", "GPR89B", "PDZK1", "PDZK1P1", 
                     "NPY4R", "NPY4R2", "FAM72A", "FAM72B", "FAM72D",   
                     "FRMPD2", "FRMPD2B"))

norm.counts.pHSDs <- norm.counts[,colnames(norm.counts) %in% pHSDs$gene_id]
colnames(norm.counts.pHSDs) <- pHSDs$name[pHSDs$gene_id %in% colnames(norm.counts.pHSDs)]

norm.counts.pHSDs.corr <- corr.test(norm.counts.pHSDs, method="pearson", ci=F)

norm.counts.pHSDs.corr$p[lower.tri(norm.counts.pHSDs.corr$p, diag=-TRUE)] = NA # correlation p-values
Pval.adj <- as.data.frame(as.table(norm.counts.pHSDs.corr$p))

norm.counts.pHSDs.corr$r [lower.tri(norm.counts.pHSDs.corr$r, diag=TRUE)] = NA # correlation coeff
Correlation <- as.data.frame(as.table(norm.counts.pHSDs.corr$r))

Cor.table <- na.exclude(cbind( Correlation, Pval.adj))[,c(1,2,3,6)]
colnames(Cor.table) <- c("gene1","gene2","cor","p.adj")

Cor.table.filt <- Cor.table[(abs(Cor.table[,3])>0.8 & Cor.table[,4] <0.05 ),]

fwrite(Cor.table, "Cor.table.CORTECON.tsv")
fwrite(Cor.table.filt, "Cor.table.filt.CORTECON.tsv")

# ##############################################################################

###### Module expression ######

pdf(paste0(selected_color,"_expression.pdf"))
heatmap.2(t(module_interest.expression), trace = "none", density = "none", 
          col = bluered(20), scale = "row", dendrogram = "row", Colv = FALSE)
dev.off()

p1 <- as.data.frame(module_interest.eigengenes) %>%
  as.data.frame() %>%
  rownames_to_column(var = "run") %>%
  pivot_longer(-run, names_to = "gene", values_to = "expression") %>%
  ggplot(aes(x = run, y = expression, group = gene)) +
  geom_line(alpha = 0.5) + theme_bw()

p2 <- module_interest.expression %>%
  as.data.frame() %>%
  rownames_to_column(var = "run") %>%
  pivot_longer(-run, names_to = "gene", values_to = "expression") %>%
  #mutate(label = ifelse(gene == "CHM13_G0018308", "y", "n")) %>%
  ggplot(aes(x = run, y = expression, group = gene)) +
  geom_line(alpha = 0.2) + scale_color_manual(values = c("grey","red")) + theme_bw()

pdf(paste0(selected_color,"_expression_line.pdf"))
grid.arrange(p1, p2)
dev.off()

###### Module membership and gene significance ######

# Most significant gene memberships
module_interest.mm.pval <- module.membership.measure.pval[rownames(module.membership.measure.pval) %in% selected_eigengene,]
module_interest.mm.pval.signif <- module_interest.mm[module_interest.mm.pval <= 0.05]
module_interest.mm.pval.signif.top <- sort(module_interest.mm.pval.signif) %>% head(100) %>% names()
module_genes.annot %>% drop_na() %>% filter(gene %in% module_interest.mm.pval.signif.top)
# None of these were SD-98 genes.

# Module membership vs gene significance
module_interest.mm <- module.membership.measure[rownames(module.membership.measure) %in% selected_eigengene,]

gene_significance.pluripotency <- as.data.frame(gene.signif.corr)$`data.Pluripotency.vs.all`
gene_significance.neural <- as.data.frame(gene.signif.corr)$`data.Neural Differentation.vs.all`
gene_significance.cortical <- as.data.frame(gene.signif.corr)$`data.Cortical Specification.vs.all`
gene_significance.deep <- as.data.frame(gene.signif.corr)$`data.Deep Layer Formation.vs.all`
gene_significance.upper <- as.data.frame(gene.signif.corr)$`data.Upper Layer Formation.vs.all`

module_interest.mm.scatterplot <- as.data.frame(
  cbind(module_interest.mm,
        gene_significance.pluripotency,
        gene_significance.neural,
        gene_significance.cortical,
        gene_significance.deep,
        gene_significance.upper))
rownames(module_interest.mm.scatterplot) <- names(module_interest.mm)

module_interest.mm.scatterplot %>% 
  pivot_longer(-module_interest.mm, names_to = "stage", values_to = "gene_significance") %>%
  ggplot(aes(x = module_interest.mm, y = gene_significance)) +
  geom_point(alpha = 0.2, color = selected_color) + 
  facet_wrap(~stage) +
  theme_bw() +
  geom_hline(yintercept = 0.8, color = "grey") + geom_hline(yintercept = -0.8, color = "grey") + 
  geom_vline(xintercept = 0.8, color = "grey") + geom_vline(xintercept = -0.8, color = "grey")
ggsave(paste0(selected_color,"module_interest_MMvsGS.pdf"))

fwrite(module_interest.mm.scatterplot, paste0(selected_color,"_membership.tsv"), sep = "\t", row.names = TRUE)

# Genes high high module membership and high gene significance are likely
# drivers of the correlation between the module and traits.
selected_stage <- gene_significance.pluripotency

module_interest.mm.scatterplot.top <- module_interest.mm.scatterplot %>% 
  filter(selected_stage >= 0.8 | selected_stage <= -0.8) %>%
  filter(module_interest.mm >= 0.8 | module_interest.mm <= -0.8) %>% 
  rownames()

module_genes.annot %>% drop_na() %>% filter(gene %in% module_interest.mm.scatterplot.top)

# Amplify or decrease the width of the edges
edgeweights <- E(module_interest.network)$weight * 2.0
module_interest.adjacency <- adjacency(module_interest.expression, power = soft_power)
module_interest.adjacency.s <- simplify(module_interest.adjacency, remove.multiple=TRUE, remove.loops=TRUE)


module_interest.network2 <- simplify(module_interest.network, remove.multiple = F, remove.loops = T)

# degree
deg <- degree(module_interest.network, mode="all")
deg.dist <- degree_distribution(module_interest.network , cumulative=T, mode="all")
pdf("test.pdf")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")
dev.off()

# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "module_intereset_edgelist.tsv",
            delim = "\t")
