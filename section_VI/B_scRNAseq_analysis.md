# Set up

## Libraries
```{r}
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(MAST)
library(harmony)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(DoubletFinder)
library(SingleCellExperiment)
library(clusterProfiler)
library(org.Dr.eg.db)
library(DESeq2)
library(edgeR)
library(viridis)
library(scCustomize)
library(DescTools)
library(dendextend)
library(dunn.test)
library(glmGamPoi)
library(hooke)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(scMiko)
library(dittoSeq)
library(BiocParallel)
library(pheatmap)

# not needed at the moment:
library(UpSetR)
library(condiments)
library(intrinsicDimension)
library(hdWGCNA)
library(WGCNA)
library(tradeSeq)
library(slingshot)
```

## Custom color palettes and functions:
```{r}
# Source of palettes: https://flatuicolors.com/

custom_colors <- list()

colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)

colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)

colors_indian <- c(
  '#FEA47F', '#25CCF7', '#EAB543', '#55E6C1', '#CAD3C8',
  '#F97F51', '#1B9CFC', '#F8EFBA', '#58B19F', '#2C3A47',
  '#B33771', '#3B3B98', '#FD7272', '#9AECDB', '#D6A2E8',
  '#6D214F', '#182C61', '#FC427B', '#BDC581', '#82589F'
)

colors_turkish <- c(
  '#cd84f1', '#ffcccc', '#ff4d4d', '#ffaf40', '#fffa65',
  '#c56cf0', '#ffb8b8', '#ff3838', '#ff9f1a', '#fff200',
  '#32ff7e', '#7efff5', '#18dcff', '#7d5fff', '#4b4b4b',
  '#3ae374', '#67e6dc', '#17c0eb', '#7158e2', '#3d3d3d'
)

custom_colors$discrete <- c(colors_dutch, colors_spanish, colors_indian, colors_turkish)

custom_colors$cell_cycle <- setNames(
  c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d'),
  c('G1',      'S',       'G2M',     '-')
)


combine_exon_intron <- function (df_gene, gene_count) 
{
    gene_count_exon = gene_count[df_gene$exon_intron == "exon", 
        ]
    gene_count_intron = gene_count[df_gene$exon_intron == "intron", 
        ]
    if (nrow(gene_count_exon) == nrow(gene_count_intron)) {
        gene_count_combine = gene_count_exon + gene_count_intron
    }
    else {
        gene_count_combine = gene_count_exon[-nrow(gene_count_exon), 
            ] + gene_count_intron
        gene_count_combine = rbind(gene_count_combine, gene_count_exon[nrow(gene_count_exon), 
            ])
    }
    return(gene_count_combine)
}

sciRNAseq_gene_count_summary <- function (gene_count_folder) {
    gene_matrix = paste(gene_count_folder, "/count.MM", sep = "")
    df_gene = paste(gene_count_folder, "/gene_name_annotate.txt", 
        sep = "")
    df_cell = paste(gene_count_folder, "/cell_annotate.txt", 
        sep = "")
    df_report = paste(gene_count_folder, "/report.MM", sep = "")
    report_annotate = paste(gene_count_folder, "/report_annotate.txt", 
        sep = "")
    df_gene = read.csv(df_gene, header = F)
    df_cell = read.csv(df_cell, header = F)
    gene_matrix = read.csv(gene_matrix, header = F)
    colnames(df_gene) = c("gene_id", "gene_type", "exon_intron", 
        "gene_name", "index")
    colnames(df_cell) = c("sample", "index")
    rownames(df_gene) = df_gene$gene_id
    rownames(df_cell) = df_cell$cell_name
    gene_count = sparseMatrix(i = gene_matrix$V1, j = gene_matrix$V2, 
        x = gene_matrix$V3)
    df_gene = df_gene[1:nrow(gene_count), ]
    rownames(gene_count) = df_gene$gene_name
    colnames(gene_count) = df_cell$cell_name
    gene_count = combine_exon_intron(df_gene, gene_count)
    df_gene = df_gene %>% filter(exon_intron == "exon")
    reportMM = read.csv(df_report, header = F)
    df_report = sparseMatrix(i = reportMM$V1, j = reportMM$V2, 
        x = reportMM$V3)
    df_report = as.matrix(t(df_report))
    df_report_annotate = read.csv(report_annotate, header = F)
    colnames(df_report) = df_report_annotate$V2
    df_report = data.frame(df_report)
    df_report["index"] = as.numeric(rownames(df_report))
    df_cell_combine = inner_join(df_cell, df_report, by = "index")
    df_cell_combine["all_exon"] = df_cell_combine$X.Perfect.intersect.exon.match + 
        df_cell_combine$X.Nearest.intersect.exon.match + df_cell_combine$X.Perfect.combine.exon.match + 
        df_cell_combine$X.Nearest.combine.exon.match
    df_cell_combine["all_intron"] = df_cell_combine$X.Perfect.intersect.gene.match + 
        df_cell_combine$X.Nearest.intersect.gene.match + df_cell_combine$X.Perfect.combine.gene.match + 
        df_cell_combine$X.Nearest.combine.gene.match
    df_cell_combine["all_reads"] = df_cell_combine$all_exon + 
        df_cell_combine$all_intron + df_cell_combine$X.No.match
    df_cell_combine["unmatched_rate"] = df_cell_combine$X.No.match/df_cell_combine$all_reads
    df_cell = df_cell_combine %>% dplyr::select(sample, unmatched_rate)
    df_cell$UMI_count = df_cell_combine$all_exon + df_cell_combine$all_intron
    df_gene = df_gene %>% dplyr::select(gene_id, gene_type, gene_name)
    return(list(df_cell, df_gene, gene_count))
}

PrctCellExpringGene <- function(object, genes, group.by = "all"){
    if(group.by == "all"){
        prct = unlist(lapply(genes,calc_helper, object=object))
        result = data.frame(Markers = genes, Cell_proportion = prct)
        return(result)
    }

    else{        
        list = SplitObject(object, group.by)
        factors = names(list)

        results = lapply(list, PrctCellExpringGene, genes=genes)
        for(i in 1:length(factors)){
        results[[i]]$Feature = factors[i]
        }
        combined = do.call("rbind", results)
        return(combined)
    }
}

calc_helper <- function(object,genes){
    counts = object[['RNA']]@counts
    ncells = ncol(counts)
    if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
    }else{return(NA)}
}

`%notin%` <- Negate(`%in%`)

```

# General processing:

## Extract data per folder:
```{r}
# Humanized libraries:
result_hum1 = sciRNAseq_gene_count_summary("sciRNAseq_hum1_newTx/gene_count")
df_cell_hum1 = result_hum1[[1]]
df_gene_hum1 = result_hum1[[2]]
gene_count_hum1 = result_hum1[[3]]

result_hum2 = sciRNAseq_gene_count_summary("sciRNAseq_hum2_newTx/gene_count")
df_cell_hum2 = result_hum2[[1]]
df_gene_hum2 = result_hum2[[2]]
gene_count_hum2 = result_hum2[[3]]

# Knockout library:
result_KO = sciRNAseq_gene_count_summary("sciRNAseq_KO_newTx/gene_count")
df_cell_KO = result_KO[[1]]
df_gene_KO = result_KO[[2]]
gene_count_KO = result_KO[[3]]

```

## Assign cells to their original cell wells based on their RT barcodes and add information on sample, genotype, and processing batch.
```{r}
plateinfo_KO <- read.csv("PlateDescription_knockout.csv", head=T)
plateinfo_hum2 <- read.csv("PlateDescription_humanized_lib2.csv", head=T)
plateinfo_hum1 <- read.csv("PlateDescription_humanized_lib1.csv", head=T)
# note: plateinfo_hum includes two RTbarcodes that are labeled "empty". this is because an initial check found 9 cells with no assigned RTbarcode from the original list, so included these at the end to avoid conflicts in ordering matrices. Will remove these 9 cells during QC anyway.

df_cell_hum1_processed <- df_cell_hum1 %>% separate(sample,c('well','fullbarcode'))
df_cell_hum1_processed$RTbarcode <- substr(df_cell_hum1_processed$fullbarcode,(nchar(df_cell_hum1_processed$fullbarcode)+1)-10,nchar(df_cell_hum1_processed$fullbarcode)) # the last 10 bases of the fullbarcode are the RT barcode
df_cell_hum1_processed$order <- as.numeric(rownames(df_cell_hum1_processed)) # need to keep the original order from the sci-RNA-seq pipeline because cells are in this order in the matrix
cellwell_hum1 <- merge(df_cell_hum1_processed, plateinfo_hum1, by="RTbarcode") # merge the 2 datasets by their RT values to match RT barcode and initial well
cellwell_hum1 <- cellwell_hum1[order(cellwell_hum1$order, decreasing = F),] # order the spreadsheet based on their original position
cellwell_hum1$sample_code <- paste(cellwell_hum1$sample, cellwell_hum1$fullbarcode, sep=".")
rownames(cellwell_hum1) <- paste(cellwell_hum1$well.x, cellwell_hum1$fullbarcode, sep=".")
rownames(cellwell_hum1) <- make.unique(cellwell_hum1$sample_code, sep="-dup")

df_cell_hum2_processed <- df_cell_hum2 %>% separate(sample,c('well','fullbarcode'))
df_cell_hum2_processed$RTbarcode <- substr(df_cell_hum2_processed$fullbarcode,(nchar(df_cell_hum2_processed$fullbarcode)+1)-10,nchar(df_cell_hum2_processed$fullbarcode))
df_cell_hum2_processed$order <- as.numeric(rownames(df_cell_hum2_processed))
cellwell_hum2 <- merge(df_cell_hum2_processed, plateinfo_hum2, by="RTbarcode")
cellwell_hum2 <- cellwell_hum2[order(cellwell_hum2$order, decreasing = F),]
cellwell_hum2$sample_code <- paste(cellwell_hum2$sample, cellwell_hum2$fullbarcode, sep=".")
rownames(cellwell_hum2) <- paste(cellwell_hum2$well.x, cellwell_hum2$fullbarcode, sep=".")
rownames(cellwell_hum2) <- make.unique(cellwell_hum2$sample_code, sep="-dup")

df_cell_KO_processed <- df_cell_KO %>%
  separate(sample,c('well','fullbarcode'))
df_cell_KO_processed$RTbarcode <- substr(df_cell_KO_processed$fullbarcode,(nchar(df_cell_KO_processed$fullbarcode)+1)-10,nchar(df_cell_KO_processed$fullbarcode))
df_cell_KO_processed$order <- rownames(df_cell_KO_processed)
df_cell_KO_processed$order <- as.numeric(df_cell_KO_processed$order)
cellwell_KO <- merge(df_cell_KO_processed, plateinfo_KO, by="RTbarcode")
cellwell_KO <- cellwell_KO[order(cellwell_KO$order, decreasing = F),]
cellwell_KO$sample_code <- paste(cellwell_KO$sample, cellwell_KO$fullbarcode, sep=".")
rownames(cellwell_KO) <- paste(cellwell_KO$well.x, cellwell_KO$fullbarcode, sep=".")
rownames(cellwell_KO) <- make.unique(cellwell_KO$sample_code, sep="-dup")

```

## Make Seurat objects per library:
```{r}
# Make unique cell and gene names, add "dup" to duplicated names:
df_gene_hum1$gene_name <- make.unique(df_gene_hum1$gene_name, sep="-dup")
df_gene_hum2$gene_name <- make.unique(df_gene_hum2$gene_name, sep="-dup")
df_gene_KO$gene_name <- make.unique(df_gene_KO$gene_name, sep="-dup")

# Set up column and rownames of the expression matrix for each sample:
colnames(gene_count_hum1) <- rownames(cellwell_hum1)
rownames(gene_count_hum1) <- df_gene_hum1$gene_name

colnames(gene_count_hum2) <- rownames(cellwell_hum2)
rownames(gene_count_hum2) <- df_gene_hum2$gene_name

colnames(gene_count_KO) <- rownames(cellwell_KO)
rownames(gene_count_KO) <- df_gene_KO$gene_name

# Create Seurat objects and add a "library" variable:
scobj_hum1 <- CreateSeuratObject(counts= gene_count_hum1, min.features = 1, min.cells = 1, project = "hum1", meta.data = cellwell_hum1) # Initial: 10,780 cells
scobj_hum2 <- CreateSeuratObject(counts= gene_count_hum2, min.features = 1, min.cells = 1, project = "hum2", meta.data = cellwell_hum2) # Initial: 82,985  cells
scobj_KO <- CreateSeuratObject(counts= gene_count_KO, min.features = 1, min.cells = 1, project = "KO", meta.data = cellwell_KO) # Initial: 114,540 cells
```

## Merge Seurat objects.
```{r}
scobj_all <- merge(scobj_hum1, y= c(scobj_hum2, scobj_KO), project="combined") # 208,305 cells
```

# Quality control and filtering

##  Filtering.
```{r}
# Calculate percentage of mitochondrial genes:
scobj_all[["percent.mt"]] <- PercentageFeatureSet(scobj_all, pattern = "^mt-")

# Calculate percentage of ribosomal genes:
scobj_all[["percent.ribo"]] <- PercentageFeatureSet(scobj_all, pattern = "^rb[sl]")

# Find doublets using DoubletFinder:
scobj_all <- NormalizeData(scobj_all)
scobj_all <- FindVariableFeatures(scobj_all, selection.method = "vst", nfeatures = 2000)
scobj_all <- ScaleData(scobj_all)
scobj_all <- RunPCA(scobj_all)
sweep <- paramSweep(scobj_all, PCs = 1:10)
sweep.stats <- summarizeSweep(sweep)
bcmvn <- find.pK(sweep.stats) # selected pK= 0.28
nExp_poi <- round(0.04*nrow(scobj_all@meta.data))  # ~4% doublet rate expected based on split-seq paper (https://www.biorxiv.org/content/10.1101/2022.08.27.505512v1.full) and sci-RNA-seq3 https://www.nature.com/articles/s41586-019-0969-x#Sec8 
scobj_all <- doubletFinder(scobj_all, PCs = 1:10, pN = 0.25, pK=0.28, nExp= nExp_poi)

# Remove samples not to be included in this project:
scobj_all$genotype <- as.factor(scobj_all$genotype)
genotypes <- levels(scobj_all$genotype)
genotypes_to_keep <- subset(genotypes, genotypes %notin% c('YTHDF2', 'ythdf2KO', 'NCF1A', 'NCF1C', 'ncf1KO', "empty"))
scobj_all <- subset(scobj_all, subset= genotype %in% genotypes_to_keep) # 189,301 initial cells
saveRDS(scobj_all, file="final_analysis/saved_datasets/scobj_all.rds")

# Stats UMI:
mean(scobj_all$UMI_count) # 1215.195
sd(scobj_all$UMI_count) # 5341.986
lower_threshold= 150
upper_threshold= mean(scobj_all$UMI_count) + 4*sd(scobj_all$UMI_count)

# Filtering data:
# Following sciRNAseq experiments from Trapnell lab: "After the single-cell gene-count matrix was generated, lower  unique molecular identifier (UMI) thresholds were determined for each experiment  (from 100–250), followed by removal of cells with UMIs greater than four standard deviations from the mean." source: https://www.nature.com/articles/s41586-023-06720-2#MOESM5

scobj_all_filt <- subset(scobj_all, subset = 
                           UMI_count > lower_threshold & 
                           UMI_count < upper_threshold & 
                           percent.mt < 5 & 
                           percent.ribo < 5 &
                           DF.classifications_0.25_0.28_8332 == "Singlet")

# Drop empty levels from the metadata:
scobj_all_filt$genotype <- droplevels(scobj_all_filt$genotype)
table(scobj_all_filt$orig.ident)
# hum1  hum2    KO 
# 8814 51185 94479 

# Plots post-filtering:
QCplots <- VlnPlot(scobj_all_filt, pt.size = 0.05, features = c("nFeature_RNA", "nCount_RNA", "UMI_count","percent.mt", "percent.ribo"), ncol=5)

ggsave("final_analysis/plots/QCplot_general_plots.pdf",
       QCplots,
       height = 8, width = 12)

QCplots2 <- ggplot(scobj_all_filt@meta.data, aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
  geom_point(size = 0.5) +
  scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
  theme_bw() +
  scale_color_viridis(
    name = '% mito',
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
  )

ggsave("final_analysis/plots/QCplot_transcript_UMIs.pdf",
       QCplots2,
       height = 8, width = 12)

# Correlation nFeature, nCounts
cor.test(scobj_all_filt$nCount_RNA, scobj_all_filt$UMI_count) # rho= 0.8738963, p-value < 2.2e-16
cor.test(scobj_all_filt$UMI_count, scobj_all_filt$nFeature_RNA) # rho= 0.956256, p-value < 2.2e-16

# Average number of UMIs per cell:
mean(scobj_all_filt$UMI_count) # 948.4134 cells

# Average number of cells per sample:
mean(table(scobj_all_filt$sample)) # 1391.694 cells

# Average number of cells per genotype:
mean(table(scobj_all_filt$genotype)) # 6179.12 cells

# Save RDS files:
saveRDS(scobj_all_filt, file="final_analysis/saved_datasets/scobj_all_filt.rds")
```

## Figure for sample descriptions.
```{r}
# Plot number of cells per genotype:
p1 <- scobj_all_filt@meta.data %>%
  group_by(genotype) %>%
         summarise(n = n()) %>%
         drop_na() %>%
  ggplot() +
  geom_bar(aes(x = genotype,
               y = n,
               fill = genotype),
           stat = "identity",
           color = "black",
           linewidth = 0.2) +
  theme_classic() +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values = custom_colors$discrete) +
  coord_flip()

p2 <- scobj_all_filt@meta.data %>%
  group_by(genotype, replicate) %>%
         summarise(n = n()) %>%
         drop_na() %>% ggplot() +
  geom_boxplot(aes(x = genotype,
                   y = log10(n),
                   fill = genotype),
               color = "black",
               size = 0.1,
               outlier.size = 0.2,
               outlier.stroke = 0) +
  theme_classic() +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values = custom_colors$discrete) +
  coord_flip()
  
ggsave("final_analysis/plots/QCplot_barplot_samples_description.pdf",
       p1 + p2,
       height = 8, width = 12)
```

# Integration and hierarchical clustering.

## Integrate data from different sequencing libraries using harmony.
```{r}
# Make sure variables are factors:
scobj_all_filt$orig.ident <- as.factor(scobj_all_filt$orig.ident) # sequencing library
scobj_all_filt$batch <- as.factor(scobj_all_filt$batch) # injection/preparation batch
scobj_all_filt$replicate <- as.factor(scobj_all_filt$replicate) # biological replicates

# Split dataset by sequencing library:
all.obj.list <- SplitObject(scobj_all_filt, split.by="orig.ident")

# Normalize each dataset using SCTransform:
all.obj.list <- lapply(X = all.obj.list, FUN = SCTransform, assay = "RNA", method = "glmGamPoi", new.assay.name = "SCT", variable.features.n = 3000, vars.to.regress = c("percent.mt", "percent.ribo"), seed.use = 1234, vst.flavor = "v2", verbose= FALSE)

# Select features for integration:
features <- SelectIntegrationFeatures(all.obj.list, nfeatures= 3000)

# Merge normalized datasets:
scobj_all_filt_norm <- Merge_Seurat_List(all.obj.list, add.cell.ids = NULL, merge.data = TRUE, project = "combined")
DefaultAssay(scobj_all_filt_norm) <- "SCT"

# Manually set variable features of merged Seurat object:
VariableFeatures(scobj_all_filt_norm) <- features

# Calculate PCs using manually set variable features:
scobj_all_filt_norm <- RunPCA(scobj_all_filt_norm, assay = "SCT", npcs = 50)

# Integration with Harmony with co-variates:
png(filename= "final_analysis/plots/integration_convergence.png")
all.data.combined <- RunHarmony(scobj_all_filt_norm, 
				group.by.vars = c("batch", "replicate"), 
				reduction = "pca", assay.use = "SCT", reduction.save = "harmony",
				kmeans_init_nstart=100, kmeans_init_iter_max=5000,
				plot_convergence=T, lambda=NULL)
dev.off()

# Visualize integration:
p <- DimPlot(all.data.combined, reduction = "harmony", group.by='genotype', raster=T)
ggsave("final_analysis/plots/PCA_integrated_genotype.pdf",
       p,
       height = 5, width = 6)

# Format co-variates in the metadata:
all.data.combined$orig.ident <- as.factor(all.data.combined$orig.ident)
all.data.combined$sample <- as.factor(all.data.combined$sample)
all.data.combined$batch <- as.factor(all.data.combined$batch)
all.data.combined$genotype <- as.factor(all.data.combined$genotype)

# Save RDS files:
saveRDS(all.data.combined, file="final_analysis/saved_datasets/all.data.combined.rds")
```

## Perform preliminary cells clustering and UMAP reductions.
```{r}
# Find neighbors:
all.data.combined <- FindNeighbors(all.data.combined, 
                                   reduction = "harmony")

all.data.combined <- RunUMAP(
  all.data.combined,
  reduction.name = 'UMAP',
  dims= 1:50,
  reduction = 'harmony',
  assay="SCT",
  seed.use = 1234)

# Plot UMAPs:
plot_umap_by_nCount <- bind_cols(all.data.combined@meta.data, as.data.frame(all.data.combined@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = nCount_RNA)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_viridis(
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
    labels = scales::comma,
  ) +
  labs(color = 'Transcripts') +
  coord_fixed()

plot_umap_by_genotype <- bind_cols(all.data.combined@meta.data, as.data.frame(all.data.combined@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = genotype)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(values = custom_colors$discrete) +
  labs(color = 'Genotype') +
  coord_fixed()

ggsave(
  'final_analysis/plots/UMAP_initial_combined.pdf',
  plot_umap_by_nCount + plot_umap_by_genotype,
  height = 8,
  width = 10
)

# Save RDS files:
saveRDS(all.data.combined, file="final_analysis/saved_datasets/all.data.combined.rds")
```

## UMAP per genotype to identify obvious abundance differences.
```{r}
p <- bind_cols(all.data.combined@meta.data, as.data.frame(all.data.combined@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = genotype)) +
  geom_point(size = 0.1, show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = custom_colors$discrete) +
  coord_fixed() +
  facet_wrap(~genotype)

ggsave(
  'final_analysis/plots/UMAPs_per_genotype.pdf',
  p,
  height = 8,
  width = 8
)
```

## Hierarchical clustering

```{r}
# Get distance between all cells from corrected-PC embeddings:
d <- dist(all.data.combined@reductions$harmony@cell.embeddings, method = "euclidean")

# Perform the hierarchical clustering using the euclidean distance:
h_euclidean <- fastcluster::hclust(d, method = "ward.D2") # use the hclust from 'fastcluster' to deal with large datasets
saveRDS(h_euclidean, file="final_analysis/saved_datasets/h_euclidean.rds")

# Cut the tree at a high K:
all.data.combined$hc_euclidean_50 <- cutree(h_euclidean, k= 50)
all.data.combined$hc_euclidean_50 <- as.factor(all.data.combined$hc_euclidean_50)

# Plot the clustering tree:
dend <- as.dendrogram(h_euclidean)
saveRDS(dend, file="final_analysis/saved_datasets/dend.rds")

# Plot UMAP with clusters:
Idents(all.data.combined) <- all.data.combined$hc_euclidean_50
pdf('final_analysis/plots/UMAP_euclidean_k50.pdf', height=10, width=10)
DimPlot_scCustom(seurat_object = all.data.combined)
dev.off()

# Obtain gene markers per cluster:
all.data.combined <- PrepSCTFindMarkers(all.data.combined, assay="SCT", verbose = TRUE)
gene.markers.k50 <- FindAllMarkers(object = all.data.combined, 
                                   assay="SCT", 
                                   logfc.threshold = 0.10, 
                                   test.use = "MAST", 
                                   only.pos = TRUE, 
                                   min.pct = 0.15,
                                   min.diff.pct = 0.10,
                                   verbose = TRUE, 
                                   random.seed = 1234)
write.csv(gene.markers.k50, file="final_analysis/data/gene.markers.k50.csv")

```

## Plotting known gene markers.
```{r}
# Plot the expression of known marker genes to narrow in on cell types:
markers.ei <- c("slc17a6a", "slc17a6b", "gad1b", "gad2")
markers.telencephalon <- c("ascl1a", "eomesa", "pvalb5", "calb2a", "epcam")
markers.diencephalon <- c("rbp4l", "zic1", "zic3", "cbln1", "tac1", "pou4f2")
markers.mesencephalon <- c("tal1", "gata3", "six3a", "pax7a", "pax7b")
markers.rhombencephalon <- c("nfia")
markers.radialglia <- c( "mfge8a", "s100b")
markers.microglia <- c("mpeg1.1", "apoeb", "apoc1", "havcr1")
markers.neuralcrest <- c("sox10")
markers.oligodendrocytes <- c("mbpa", "olig1", "olig2")
markers.cornea <- c("matn4", "col1a2")
markers.lens <- c("crybb1", "cryba4", "crygn2")
markers.RGC <- c( "pax6a", "robo2", "isl2b", "rbpms2b")
markers.horizontal <- c("ompa", "rprmb", "isl1")
markers.amacrine <- c("ptf1a", "hes2.2", "atoh7")
markers.bipolar <- c("gnb3a", "vsx1", "crx")
markers.PR <- c("opn1mw1", "pde6h", "rho")
markers.RPE <- c("pmela")
markers.muellerglia <- c("pax10", "tkta", "glula", "glulb", "gfap")
markers.progenitors <- c("sox2", "sox19a")

DefaultAssay(object = all.data.combined) <- "RNA"
pdf("final_analysis/plots/UMAP_euclidean_k50_ei.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.ei)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_telencephalon.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.telencephalon)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_diencencephalon.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.diencephalon)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_mesencephalon.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.mesencephalon)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_rhombencephalon.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.rhombencephalon)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_radialglia.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.radialglia)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_microglia.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.microglia)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_neuralcrest.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.neuralcrest)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_oligodendrocytes.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.oligodendrocytes)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_cornea.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.cornea)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_lens.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.lens)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_RGC.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.RGC)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_horizontal.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.horizontal)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_amacrine.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.amacrine)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_bipolar.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.bipolar)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_PR.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.PR)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_RPE.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.RPE)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_muellerglia.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.muellerglia)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_ei.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.ei)
dev.off()
pdf("final_analysis/plots/UMAP_euclidean_k50_neuralprogenitors.pdf")
FeaturePlot_scCustom(seurat_object = all.data.combined, features = markers.progenitors)
dev.off()

pdf("final_analysis/plots/Violin_euclidean_k50_ei.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.ei, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_telencephalon.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.telencephalon, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_diencencephalon.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.diencephalon, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_mesencephalon.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.mesencephalon, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_rhombencephalon.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.rhombencephalon, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_radialglia.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.radialglia, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_microglia.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.microglia, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_neuralcrest.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.neuralcrest, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_oligodendrocytes.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.oligodendrocytes, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_cornea.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.cornea, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_lens.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.lens, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_RGC.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.RGC, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_horizontal.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.horizontal, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_amacrine.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.amacrine, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_bipolar.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.bipolar, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_PR.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.PR, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_RPE.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.RPE, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_mullerglia.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.mullerglia, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_ei.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.ei, pt.size = 0)
dev.off()
pdf("final_analysis/plots/Violin_euclidean_k50_neuralprogenitors.pdf", height = 10, width= 10)
VlnPlot_scCustom(seurat_object = all.data.combined, features = markers.progenitors, pt.size = 0)
dev.off()

DefaultAssay(object = all.data.combined) <- "SCT"

```

## Use the Raj et al 2018 reference data to identify common clusters and further confidently define cell types.
```{r}
## Load reference dataset from Raj et al 3dpf:
reference_obj <- readRDS("reference_datasets/Rajetal2018_3dpf_updated.rds")
reference_obj$res.5 <- as.factor(reference_obj$res.5)
reference_obj <- SCTransform(reference_obj, vars.to.regress = "pt.mito") # need to have both objects with same normalization type
DefaultAssay(reference_obj) <- "SCT"
reference_obj <- RunPCA(reference_obj, npcs = 50, verbose = FALSE) # need to re-run PCA with same number of dimensions

# Load marker genes (resolution used by Raj et al = 5):
reference_markers <- read.csv("reference_datasets/Rajetal2018_3dpf_markers.csv", head=T)

# Add marker IDs (broad and narrow) to the reference cells:
broad.celltypes <- reference_markers$clusterID_broad
narrow.celltypes <- reference_markers$clusterID_narrow
names(broad.celltypes) <- levels(reference_obj$res.5)
names(narrow.celltypes) <- levels(reference_obj$res.5)

Idents(reference_obj) <- "res.5"
reference_obj <- RenameIdents(reference_obj, broad.celltypes)
reference_obj$broad.celltypes <- Idents(reference_obj)

Idents(reference_obj) <- "res.5"
reference_obj <- RenameIdents(reference_obj, narrow.celltypes)
reference_obj$narrow.celltypes <- Idents(reference_obj)

# Transfer labels from reference to current sciRNAseq data:

ref.anchors <- FindTransferAnchors(
  reference_obj,
  all.data.combined,
  normalization.method = "SCT",
  recompute.residuals = TRUE,
  reduction = "pcaproject",
  reference.reduction = "pca",
  project.query = TRUE,
  scale = TRUE,
  npcs = 50,
  l2.norm = TRUE,
  dims = 1:50,
  k.anchor = 5,
  k.score = 30,
  max.features = 200,
  nn.method = "annoy",
  n.trees = 50,
  eps = 0,
  approx.pca = TRUE,
  mapping.score.k = NULL,
  verbose = TRUE
)

ref.labels <- TransferData(
  anchorset = ref.anchors,
  refdata= reference_obj$narrow.celltypes,
  reference = reference_obj,
  query = all.data.combined,
  query.assay = "SCT",
  weight.reduction = "pcaproject",
  l2.norm = FALSE,
  k.weight = 50,
  sd.weight = 1,
  eps = 0,
  n.trees = 50,
  verbose = TRUE,
  slot = "data",
  prediction.assay = FALSE,
  only.weights = FALSE,
  store.weights = TRUE
)

```

# Renaming cell cluster identities and retain cell types of interest.
```{r}
# After processing the spreadsheets with gene markers the cluster identities at K50 are:
markers.processed <- read.csv("final_analysis/data/gene.markers.k50_summary.csv", head=T)
major.celltypes <- markers.processed$major_celltype
names(major.celltypes) <- levels(all.data.combined$hc_euclidean_50)

# Save new classifications to metadata following order of current identities:
Idents(all.data.combined) <- "hc_euclidean_50"
all.data.combined <- RenameIdents(all.data.combined, major.celltypes)
all.data.combined$major.celltypes <- Idents(all.data.combined)

# Select clusters of interest
all.data.combined.selected <- subset(all.data.combined, major.celltypes %in% c("neuronal", "eye", "glial")) # 116,180 cells

# Average number of cells per sample/genotype:
mean(table(all.data.combined.selected$sample)) # 1046.667 cells
sd(table(all.data.combined.selected$sample)) # 1019.343 cells

mean(table(all.data.combined.selected$genotype)) # 4647.2 cells
sd(table(all.data.combined.selected$genotype)) # 3672.988 cells

saveRDS(all.data.combined.selected, file="final_analysis/saved_datasets/all.data.combined.selected.rds")
```

# Re-cluster and re-do UMAP reduction on subsetted cells.
```{r}
# Get distance between all cells from corrected-PC embeddings:
d2 <- dist(all.data.combined.selected@reductions$harmony@cell.embeddings, method = "euclidean")

# Perform the hierarchical clustering using the euclidean distance:
h_euclidean2 <- fastcluster::hclust(d2, method = "ward.D2") # use the hclust from 'fastcluster' to deal with large datasets
saveRDS(h_euclidean2, file="final_analysis/saved_datasets/h_euclidean2.rds")

# Cut the tree at a high K:
all.data.combined.selected$hc_subset_euclidean_50 <- cutree(h_euclidean2, k= 50)
all.data.combined.selected$hc_subset_euclidean_50 <- as.factor(all.data.combined.selected$hc_subset_euclidean_50)

# Obtain gene markers per cluster:
Idents(all.data.combined.selected) <- all.data.combined.selected$hc_subset_euclidean_50
all.data.combined.selected <- PrepSCTFindMarkers(all.data.combined.selected, assay="SCT", verbose = TRUE)

gene.markers.selected.k50 <- FindAllMarkers(object = all.data.combined.selected, 
                                            assay="SCT", 
                                            logfc.threshold = 0.10, 
                                            test.use = "MAST", 
                                            only.pos = TRUE, 
                                            min.pct = 0.15, 
                                            random.seed = 1234, 
                                            verbose = TRUE)

write.csv(gene.markers.selected.k50, file="final_analysis/data/gene.markers.selected.k50.csv")
write.csv(as.data.frame(table(all.data.combined.selected$hc_subset_euclidean_50)), file="final_analysis/data/gene.markers.selected.k50.clustercounts.csv")

# Run UMAP:
all.data.combined.selected <- RunUMAP(
  all.data.combined.selected,
  reduction.name = 'UMAP',
  dims= 1:50,
  reduction = 'harmony',
  assay="SCT",
  seed.use = 1234)

# Plot UMAP with clusters (basic):
Idents(all.data.combined.selected) <- all.data.combined.selected$hc_subset_euclidean_50
pdf('final_analysis/plots/UMAP_subset_euclidean_k50.pdf', height=10, width=10)
DimPlot_scCustom(seurat_object = all.data.combined.selected)
dev.off()

# Plot UMAP with clusters (with labels for publication):

UMAP_centers_cluster <- tibble(
    UMAP_1 = as.data.frame(all.data.combined.selected@reductions$UMAP@cell.embeddings)$UMAP_1,
    UMAP_2 = as.data.frame(all.data.combined.selected@reductions$UMAP@cell.embeddings)$UMAP_2,
    cluster = all.data.combined.selected@meta.data$hc_subset_euclidean_50
  ) %>%
  group_by(cluster) %>%
  summarize(x = median(UMAP_1), y = median(UMAP_2))

plot_umap_selected <- bind_cols(all.data.combined.selected@meta.data, as.data.frame(all.data.combined.selected@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = hc_subset_euclidean_50)) +
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
  'final_analysis/plots/UMAP_subset_euclidean_k50_labels.pdf',
  plot_umap_selected,
  height = 6,
  width = 6
)

# Save current object.
saveRDS(all.data.combined.selected, file="final_analysis/saved_datasets/all.data.combined.selected.rds")

## Select the 'neurons' cell types and do further clustering of these since groups are still very broad:

# Load processed cluser ID based on Daniocell/ZFIN databases:
markers.subset.k50.processed <- read.csv("final_analysis/data/gene.markers.selected.k50_summary_temp.csv")
major.celltypes <- markers.subset.k50.processed$major_celltype
names(major.celltypes) <- levels(all.data.combined.selected$hc_euclidean_50)

# Save new classifications to metadata following order of current identities:
Idents(all.data.combined.selected) <- "hc_subset_euclidean_50"
all.data.combined.selected <- RenameIdents(all.data.combined.selected, major.celltypes)
all.data.combined.selected$major.celltypes.processed <- Idents(all.data.combined.selected)

# Select clusters of interest
all.data.combined.selected.neuronal <- subset(all.data.combined.selected, major.celltypes.processed %in% c("neurons")) # 65,013 cells
all.data.combined.selected.neuronal$hc_subset_euclidean_50 <- droplevels(all.data.combined.selected.neuronal$hc_subset_euclidean_50) # remove empty non-neuronal clusters

# Obtain gene markers per neuronal cluster:
Idents(all.data.combined.selected.neuronal) <- all.data.combined.selected.neuronal$hc_subset_euclidean_50
gene.markers.selected.k50.neuronal <- FindAllMarkers(object = all.data.combined.selected.neuronal, 
                                            assay="SCT", 
                                            logfc.threshold = 0.10, 
                                            test.use = "MAST", 
                                            only.pos = TRUE, 
                                            random.seed = 1234, 
                                            verbose = TRUE,
                                            recorrect_umi= FALSE)

write.csv(gene.markers.selected.k50.neuronal, file="final_analysis/data/gene.markers.selected.k50.neuronal.csv")

# Plot known markers from Zhang et al 2021
DefaultAssay(object = all.data.combined.selected) <- "RNA"

pdf("final_analysis/plots/UMAP_subset_euclidean_k50_neuronal_forebrain.pdf", height = 10, width= 10)
FeaturePlot_scCustom(seurat_object = all.data.combined.selected, features = c("foxg1a", "gnrh3", "egr1", "mcl1a", "junbb", "cbx7a"), pt.size = 0)
dev.off()

pdf("final_analysis/plots/UMAP_subset_euclidean_k50_neuronal_pallium.pdf", height = 10, width= 10)
FeaturePlot_scCustom(seurat_object = all.data.combined.selected, features = c("dlx5a"), pt.size = 0)
dev.off()

pdf("final_analysis/plots/UMAP_subset_euclidean_k50_neuronal_subpallium.pdf", height = 10, width= 10)
FeaturePlot_scCustom(seurat_object = all.data.combined.selected, features = c("eomesa"), pt.size = 0)
dev.off()

pdf("final_analysis/plots/UMAP_subset_euclidean_k50_neuronal_diencephalon.pdf", height = 10, width= 10)
FeaturePlot_scCustom(seurat_object = all.data.combined.selected, features = c("pitx2"), pt.size = 0)
dev.off()

pdf("final_analysis/plots/UMAP_subset_euclidean_k50_neuronal_thalamus.pdf", height = 10, width= 10)
FeaturePlot_scCustom(seurat_object = all.data.combined.selected, features = c("barhl2"), pt.size = 0)
dev.off()

pdf("final_analysis/plots/UMAP_subset_euclidean_k50_neuronal_optictectum.pdf", height = 10, width= 10)
FeaturePlot_scCustom(seurat_object = all.data.combined.selected, features = c("tal1", "en2a"), pt.size = 0)
dev.off()

pdf("final_analysis/plots/UMAP_subset_euclidean_k50_neuronal_hindbrain.pdf", height = 10, width= 10)
FeaturePlot_scCustom(seurat_object = all.data.combined.selected, features = c("hoxa3a", "hoxb3a"), pt.size = 0)
dev.off()

pdf("final_analysis/plots/UMAP_subset_euclidean_k50_neuronal_glut.pdf", height = 10, width= 10)
FeaturePlot_scCustom(seurat_object = all.data.combined.selected, features = c("slc17a6b"), pt.size = 0)
dev.off()

pdf("final_analysis/plots/UMAP_subset_euclidean_k50_neuronal_gaba.pdf", height = 10, width= 10)
FeaturePlot_scCustom(seurat_object = all.data.combined.selected, features = c("gad1b", "gad2"), pt.size = 0)
dev.off()

pdf("final_analysis/plots/UMAP_subset_euclidean_k50_neuronal_glyc.pdf", height = 10, width= 10)
FeaturePlot_scCustom(seurat_object = all.data.combined.selected, features = c("slc6a5"), pt.size = 0)
dev.off()

DefaultAssay(object = all.data.combined.selected) <- "SCT"


## Final re-naming of clusters:

# Load processed cluser ID based on Daniocell/ZFIN databases:
markers.subset.k50.processed.final <- read.csv("final_analysis/data/gene.markers.selected.k50_summary_final.csv")
major.celltypes <- markers.subset.k50.processed.final$major_celltype
broad.celltypes <- markers.subset.k50.processed.final$broad_celltype
specific.celltypes <- markers.subset.k50.processed.final$specific_celltype

names(major.celltypes) <- levels(all.data.combined.selected$hc_euclidean_50)
names(broad.celltypes) <- levels(all.data.combined.selected$hc_euclidean_50)
names(specific.celltypes) <- levels(all.data.combined.selected$hc_euclidean_50)

# Save new classifications to metadata following order of current identities:
Idents(all.data.combined.selected) <- "hc_subset_euclidean_50"
all.data.combined.selected <- RenameIdents(all.data.combined.selected, major.celltypes)
all.data.combined.selected$major.celltypes.final <- Idents(all.data.combined.selected)

Idents(all.data.combined.selected) <- "hc_subset_euclidean_50"
all.data.combined.selected <- RenameIdents(all.data.combined.selected, broad.celltypes)
all.data.combined.selected$broad.celltypes.final <- Idents(all.data.combined.selected)

Idents(all.data.combined.selected) <- "hc_subset_euclidean_50"
all.data.combined.selected <- RenameIdents(all.data.combined.selected, specific.celltypes)
all.data.combined.selected$specific.celltypes.final <- Idents(all.data.combined.selected)

# Remove cell types not of interest:

all.data.combined.selected.filtered <- subset(all.data.combined.selected, major.celltypes.final %in% c("neurons", "glia", "olfactory", "retinal_neural", "retinal_noneural")) # 95,555
all.data.combined.selected.filtered$specific.celltypes.final <- droplevels(all.data.combined.selected.filtered$specific.celltypes.final)


# Order clusters for plotting:
all.data.combined.selected.filtered$specific.celltypes.final <- factor(all.data.combined.selected.filtered$specific.celltypes.final, levels=c("forebrain", "pallium", "hypothalamus", "thalamus", "olfactory_organ", "midbrain", "hindbrain", "spinal_cord","RGC", "bipolar_cells", "amacrine_horizontal", "photoreceptors_precursors", "cones", "rods", "radial_glia", "mueller_glia", "microglia"))

# Get counts per final cluster:
write.csv(as.data.frame(table(all.data.combined.selected.filtered$specific.celltypes.final)), file="final_analysis/data/gene.markers.selected.specific.celltypes.counts.csv")

p_cluster_counts <- ggplot(data= as.data.frame(table(all.data.combined.selected.filtered$specific.celltypes.final)), 
                           aes(x = Var1, y = 1, size = Freq)) +
  geom_point(alpha = 0.8) +
  theme_minimal()
ggsave("final_analysis/plots/balloonplot_cluster_counts.pdf",
       p_cluster_counts, width= 6, height= 4)

Idents(all.data.combined.selected.filtered) <- "specific.celltypes.final"

# Run dimensionality reductions:
all.data.combined.selected.filtered <- RunUMAP(
  all.data.combined.selected.filtered,
  reduction.name = 'UMAP',
  dims= 1:50,
  reduction = 'harmony',
  assay="SCT",
  seed.use = 1234)

all.data.combined.selected.filtered <- RunTSNE(
  all.data.combined.selected.filtered,
  reduction.name = 'TSNE',
  dims= 1:50,
  reduction = 'harmony',
  assay="SCT",
  seed.use = 1234)

# Plot UMAP/tSNE with clusters (with labels for publication):

UMAP_centers_cluster <- tibble(
    UMAP_1 = as.data.frame(all.data.combined.selected.filtered@reductions$UMAP@cell.embeddings)$UMAP_1,
    UMAP_2 = as.data.frame(all.data.combined.selected.filtered@reductions$UMAP@cell.embeddings)$UMAP_2,
    cluster = all.data.combined.selected.filtered@meta.data$specific.celltypes.final
  ) %>%
  group_by(cluster) %>%
  summarize(x = median(UMAP_1), y = median(UMAP_2))

plot_umap_final <- bind_cols(all.data.combined.selected.filtered@meta.data, as.data.frame(all.data.combined.selected.filtered@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = specific.celltypes.final)) +
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
  'final_analysis/plots/UMAP_subset_euclidean_final_labels.pdf',
  plot_umap_final,
  height = 6,
  width = 6
)

TSNE_centers_cluster <- tibble(
    tSNE_1 = as.data.frame(all.data.combined.selected.filtered@reductions$TSNE@cell.embeddings)$tSNE_1,
    tSNE_2 = as.data.frame(all.data.combined.selected.filtered@reductions$TSNE@cell.embeddings)$tSNE_2,
    cluster = all.data.combined.selected.filtered@meta.data$specific.celltypes.final
  ) %>%
  group_by(cluster) %>%
  summarize(x = median(tSNE_1), y = median(tSNE_2))

plot_tsne_selected <- bind_cols(all.data.combined.selected.filtered@meta.data, as.data.frame(all.data.combined.selected.filtered@reductions$TSNE@cell.embeddings)) %>%
  ggplot(aes(tSNE_1, tSNE_2, color = specific.celltypes.final)) +
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
  'final_analysis/plots/TSNE_subset_euclidean_final_labels.pdf',
  plot_tsne_selected,
  height = 6,
  width = 6
)

# Save final object:
saveRDS(all.data.combined.selected.filtered, file="final_analysis/saved_datasets/all.data.combined.selected.filtered.rds")
```


# Cell cycle analyses
```{r}
# Load cell cycle markers:
ccmarkers <- read.csv("cell_cycle_markers.csv", head=T, stringsAsFactors = T)

# Infer cell cycle status:
all.data.combined.selected.filtered <- CellCycleScoring(
  all.data.combined.selected.filtered,
  assay = 'SCT',
  s.features = ccmarkers$Zebrafish_gene[ccmarkers$phase=="S"],
  g2m.features = ccmarkers$Zebrafish_gene[ccmarkers$phase=="G2/M"]
)

all.data.combined.selected.filtered@meta.data$cell_cycle_seurat <- all.data.combined.selected.filtered@meta.data$Phase
all.data.combined.selected.filtered@meta.data$Phase <- NULL
all.data.combined.selected.filtered@meta.data$cell_cycle_seurat <- factor(all.data.combined.selected.filtered@meta.data$cell_cycle_seurat, levels = c('G1', 'S', 'G2M'))

# Plot:
pdf('final_analysis/plots/UMAP_selected_cell_cycle.pdf', height=5, width=5)
Meta_Highlight_Plot(seurat_object = all.data.combined.selected.filtered, meta_data_column = "cell_cycle_seurat",
    meta_data_highlight = c("G1", "S", "G2M"), highlight_color = c("#45aaf2", "#f1c40f", "#e74c3c"),
    background_color = "lightgray")
dev.off()
```

# Expression pattern of the HSD genes.
```{r}
# Heatmap with Z-score across tissues:
meta <- all.data.combined.selected.filtered@meta.data
HSDs <- c("srgap2", "arhgap11a", "GPR89B", "npy4r", "FAM72B", "frmpd2", "hydin", "pdzk1", "ptpn20")
expr <- all.data.combined.selected.filtered@assays$RNA$data[HSDs,] %>% as.matrix()
heatmap_df <- bind_cols(data.frame(clusters=meta$specific.celltypes.final), t(expr)) %>% group_by(clusters) %>% summarize_all(mean) %>% as.data.frame()
rownames(heatmap_df) <- heatmap_df$clusters
heatmap_df <- heatmap_df[,-1] %>% as.matrix()

pdf("final_analysis/plots/HSDs_expression_heatmap.pdf", height=8, width=8)
pheatmap::pheatmap(heatmap_df,cluster_rows = FALSE, cluster_cols = TRUE, scale = "column")
dev.off()

```

# Correlation in fold-change across samples
```{r}
# Subsample to balance numbers of cells across genotypes.
all.data.combined.selected.cor <- balanceSamples(object = all.data.combined.selected.filtered, group = "genotype")

# Merge control genotypes to have a common control across groups:
all.data.combined.selected.cor$genotype_grouped <- fct_recode(all.data.combined.selected.cor$genotype,
                                                     control= "scrambled",
                                                     control= "GFP",
                                                     control= "controlMO")

# Pseudobulk the counts based on donor-condition-celltype (recommended via Seurat tutorial) https://satijalab.org/seurat/articles/de_vignette
all.data.combined.selected.cor.pseudo <- AggregateExpression(all.data.combined.selected.cor, assays = "SCT", return.seurat = T, group.by = c("genotype", "sample", "specific.celltypes.final"))

## Correlations between samples:
data.deseq2.all <- DESeqDataSetFromMatrix(countData = all.data.combined.selected.cor.pseudo@assays$SCT$counts+1,
                                          colData= all.data.combined.selected.cor.pseudo@meta.data,
                                          design= ~genotype)
# Perform DGE:
data.deseq2.all <- estimateSizeFactors(data.deseq2.all)
data.deseq2.all <- estimateDispersionsGeneEst(data.deseq2.all)
dispersions(data.deseq2.all) <- mcols(data.deseq2.all)$dispGeneEst
data.deseq2.all <- nbinomWaldTest(data.deseq2.all)

# Get correlation:
res.KO <- as.data.frame(results(data.deseq2.all, contrast = c("genotype", "gpr89KO", "controlMO"), parallel=T))
res.hum <- as.data.frame(results(data.deseq2.all, contrast = c("genotype", "GPR89B", "GFP"), parallel=T))

# Correlation:
merged.df <- merge(res.KO, res.hum, by= "row.names") # all genes
merged.df <- merged.df[complete.cases(merged.df),]
correlation <- cor.test(merged.df$log2FoldChange.x, merged.df$log2FoldChange.y, method= "spearman")

# Plot:
rho <- correlation$estimate
n <- nrow(merged.df)
se_rho <- sqrt((1 - rho^2) / (n - 2))

correlation_line <- function(x, rho, x_mean, y_mean) {
  return(rho * (x - x_mean) + y_mean)
}
# Calculate mean of fold changes
x_mean <- mean(merged.df$log2FoldChange.x)
y_mean <- mean(merged.df$log2FoldChange.y)
# Create a data frame for the line
line_data <- data.frame(
  x = seq(min(merged.df$log2FoldChange.x), max(merged.df$log2FoldChange.x), length.out = 100)
)
line_data$y <- correlation_line(line_data$x, rho, x_mean, y_mean)
# Create data frames for the upper and lower lines
upper_line_data <- data.frame(
  x = line_data$x,
  y = correlation_line(line_data$x, rho + se_rho, x_mean, y_mean)
)

lower_line_data <- data.frame(
  x = line_data$x,
  y = correlation_line(line_data$x, rho - se_rho, x_mean, y_mean)
)

# Make plot:
p <- ggplot(data= merged.df, aes(x= log2FoldChange.x, y= log2FoldChange.y)) +
  geom_point(col= "gray", size= 0.2, alpha= 0.5)  +
  geom_smooth(method= "lm", se=T, col="blue") +
geom_vline(xintercept = 0, lty= 2) + 
  geom_hline(yintercept = 0, lty=2) +
geom_line(data = line_data, aes(x = x, y = y), color = "black") +
  geom_line(data = upper_line_data, aes(x = x, y = y), color = "black", linetype = "dashed") +
  geom_line(data = lower_line_data, aes(x = x, y = y), color = "black", linetype = "dashed") +  xlab("log2(FC) KO") + ylab("log2(FC) Hum")  +
  theme_bw()

ggsave("final_analysis/plots/DEGs_correlation_SRGAP2.pdf",
       p,
       height = 2,
       width= 2)

####### Get cell-type-specific DEGs:

# Create functions to use:
perform_GO_enrichment <- function(deg_list, ontology = "BP") {
  ego <- enrichGO(gene          = deg_list,
                  OrgDb         = org.Dr.eg.db,
                  keyType       = 'SYMBOL',
                  ont           = ontology,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.1,
                  qvalueCutoff  = 0.1)
  return(ego)
}

collect_GO_terms <- function(deg_list, genotype_name) {
  ego <- perform_GO_enrichment(deg_list)
  if (!is.null(ego)) {
    ego_filtered <- gofilter(ego, level=4)
    go_terms <- ego_filtered@result
    
    if (nrow(go_terms) == 0) {
      return(NULL)
    }
      go_terms$genotype <- genotype_name
    
    return(go_terms)
  } else {
    return(NULL)
  }
}

# Subsample to balance numbers of cells across genotypes
all.data.combined.selected.cor <- balanceSamples(object = all.data.combined.selected.filtered, group = "genotype")

# Pseudobulk the counts based on donor-condition-celltype
all.data.combined.selected.cor.pseudo <- AggregateExpression(
  all.data.combined.selected.cor,
  assays = "SCT",
  return.seurat = TRUE,
  group.by = c("genotype", "sample", "specific.celltypes.final")
)

# List of genotypes and contrasts
genotypes <- c("srgap2KO", "SRGAP2C", "ARHGAP11B", "gpr89KO", "GPR89B", "npy4rKO", "NPY4R","fam72KO", "FAM72B", "frmpd2KO", "FRMPD2B", "hydinKO", "pdzk1KO", "PDZK1P1","ptpn20KO", "PTPN20CP")

contrasts <- list(c("orig.ident", "srgap2KO", "scrambled"),
                  c("orig.ident", "SRGAP2C", "GFP"),
                  c("orig.ident", "ARHGAP11B", "GFP"),
                  c("orig.ident", "gpr89KO", "scrambled"),
                  c("orig.ident", "GPR89B", "GFP"),
                  c("orig.ident", "npy4rKO", "scrambled"),
                  c("orig.ident", "NPY4R", "GFP"),
                  c("orig.ident", "fam72KO", "scrambled"),
                  c("orig.ident", "FAM72B", "GFP"),
                  c("orig.ident", "frmpd2KO", "scrambled"),
                  c("orig.ident", "FRMPD2B", "GFP"),
                  c("orig.ident", "hydinKO", "scrambled"),
                  c("orig.ident", "pdzk1KO", "scrambled"),
                  c("orig.ident", "PDZK1P1", "GFP"),
                  c("orig.ident", "ptpn20KO", "scrambled"),
                  c("orig.ident", "PTPN20CP", "GFP"))

# Forebrain:

# Extract rows that correspond to 'forebrain'
forebrain_rows <- grep("forebrain", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'forebrain' data
forebrain.data <- subset(all.data.combined.selected.cor.pseudo, cells = forebrain_rows)

# Create DESeq2 dataset
forebrain.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(forebrain.data@assays$SCT@counts + 1),
  colData = forebrain.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
forebrain.deseq2 <- estimateSizeFactors(forebrain.deseq2)
forebrain.deseq2 <- estimateDispersionsGeneEst(forebrain.deseq2)
dispersions(forebrain.deseq2) <- mcols(forebrain.deseq2)$dispGeneEst
forebrain.deseq2 <- nbinomWaldTest(forebrain.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(forebrain.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/forebrain_DEGs_complete_results.csv", row.names = TRUE)

# Collect significant DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(forebrain.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res),]
  rownames(res[res$padj < 0.1 & abs(res$log2FoldChange) > 0,])
})

# Perform GO enrichment and collect GO terms
go_terms_list <- lapply(seq_along(genotypes), function(i) {
  genotype <- genotypes[i]
  deg_list <- deg_lists[[i]]
  collect_GO_terms(deg_list, genotype)
})

# Combine all GO terms into one data frame and get top 5 per genotype:
all_go_terms <- bind_rows(go_terms_list)
top_go_terms <- all_go_terms %>%
  group_by(genotype) %>%  
  arrange(qvalue) %>%       
  slice_min(qvalue, n = 5)

# Aggregate duplicated GO terms by taking the minimum p-value for each genotype
go_aggregated <- top_go_terms %>%
  group_by(Description, genotype) %>%
  summarise(qvalue= qvalue)

# Spread the data into a wide format (GO terms x genotypes)
go_matrix_wide <- go_aggregated %>%
  spread(key = genotype, value = qvalue, fill = NA)

# Convert to a matrix for heatmap
go_matrix_qvalues <- as.matrix(go_matrix_wide[,-1])  # Exclude the Description column
rownames(go_matrix_qvalues) <- go_matrix_wide$Description

# Save as a CSV file:
write.csv(go_matrix_qvalues, file="final_analysis_new/forebrain_GOs_FC_matrix.csv")

# Plot the heatmap

# Apply log2 transformation to q-values
go_matrix_log2 <- -log2(go_matrix_qvalues)
go_matrix_log2[is.infinite(go_matrix_log2)] <- 0
go_matrix_log2[is.na(go_matrix_log2)] <- 0
go_matrix_log2[(go_matrix_log2) >= 12] <- 12 # for visual purposes

# Custom row and column labels
row_labels <- rownames(go_matrix_log2)
col_labels <- colnames(go_matrix_log2)

# Plot the heatmap without clustering
pdf("final_analysis_new/plots/forebrain_GOs_qvalue_heatmap.pdf", width=8, height= 10)
pheatmap(go_matrix_log2,
         color = colorRampPalette(c("white", "red"))(10),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         labels_row = row_labels,
         labels_col = col_labels)
dev.off()

pdf("final_analysis_new/plots/forebrain_GOs_qvalue_heatmap_colsclustered.pdf", width=8, height= 10)
pheatmap(go_matrix_log2,
         color = colorRampPalette(c("white", "red"))(10),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         labels_row = row_labels,
         labels_col = col_labels)
dev.off()

# Pallium:

# Extract rows that correspond to 'pallium'
pallium_rows <- grep("pallium", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'pallium' data
pallium.data <- subset(all.data.combined.selected.cor.pseudo, cells = pallium_rows)

# Create DESeq2 dataset
pallium.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(pallium.data@assays$SCT@counts + 1),
  colData = pallium.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
pallium.deseq2 <- estimateSizeFactors(pallium.deseq2)
pallium.deseq2 <- estimateDispersionsGeneEst(pallium.deseq2)
dispersions(pallium.deseq2) <- mcols(pallium.deseq2)$dispGeneEst
pallium.deseq2 <- nbinomWaldTest(pallium.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(pallium.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/pallium_DEGs_complete_results.csv", row.names = TRUE)

# Hypothalamus:

# Extract rows that correspond to 'hypothalamus'
hypothalamus_rows <- grep("hypothalamus", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'hypothalamus' data
hypothalamus.data <- subset(all.data.combined.selected.cor.pseudo, cells = hypothalamus_rows)

# Create DESeq2 dataset
hypothalamus.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(hypothalamus.data@assays$SCT@counts + 1),
  colData = hypothalamus.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
hypothalamus.deseq2 <- estimateSizeFactors(hypothalamus.deseq2)
hypothalamus.deseq2 <- estimateDispersionsGeneEst(hypothalamus.deseq2)
dispersions(hypothalamus.deseq2) <- mcols(hypothalamus.deseq2)$dispGeneEst
hypothalamus.deseq2 <- nbinomWaldTest(hypothalamus.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(hypothalamus.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/hypothalamus_DEGs_complete_results.csv", row.names = TRUE)

# Thalamus:

# Extract rows that correspond to 'thalamus'
thalamus_rows <- grep("thalamus", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'thalamus' data
thalamus.data <- subset(all.data.combined.selected.cor.pseudo, cells = thalamus_rows)

# Create DESeq2 dataset
thalamus.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(thalamus.data@assays$SCT@counts + 1),
  colData = thalamus.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
thalamus.deseq2 <- estimateSizeFactors(thalamus.deseq2)
thalamus.deseq2 <- estimateDispersionsGeneEst(thalamus.deseq2)
dispersions(thalamus.deseq2) <- mcols(thalamus.deseq2)$dispGeneEst
thalamus.deseq2 <- nbinomWaldTest(thalamus.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(thalamus.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/thalamus_DEGs_complete_results.csv", row.names = TRUE)

# Midbrain:

# Extract rows that correspond to 'midbrain'
midbrain_rows <- grep("midbrain", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'forebrain' data
midbrain.data <- subset(all.data.combined.selected.cor.pseudo, cells = midbrain_rows)

# Create DESeq2 dataset
midbrain.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(midbrain.data@assays$SCT@counts + 1),
  colData = midbrain.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
midbrain.deseq2 <- estimateSizeFactors(midbrain.deseq2)
midbrain.deseq2 <- estimateDispersionsGeneEst(midbrain.deseq2)
dispersions(midbrain.deseq2) <- mcols(midbrain.deseq2)$dispGeneEst
midbrain.deseq2 <- nbinomWaldTest(midbrain.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(midbrain.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/midbrain_DEGs_complete_results.csv", row.names = TRUE)

# Collect significant DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(midbrain.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res),]
  rownames(res[res$padj < 0.1 & abs(res$log2FoldChange) > 0,])
})

# Perform GO enrichment and collect GO terms
go_terms_list <- lapply(seq_along(genotypes), function(i) {
  genotype <- genotypes[i]
  deg_list <- deg_lists[[i]]
  collect_GO_terms(deg_list, genotype)
})

# Combine all GO terms into one data frame
all_go_terms <- bind_rows(go_terms_list)

all_go_terms <- all_go_terms %>%
  mutate(
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    BgRatio_num = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    FoldChange = GeneRatio_num / BgRatio_num
  )

# Function to extract top GO terms for a genotype
top_go_terms <- all_go_terms %>%
  group_by(genotype) %>%  
  arrange(qvalue) %>%       
  slice_min(qvalue, n = 5)

# Aggregate duplicated GO terms by taking the minimum p-value for each genotype
go_aggregated <- top_go_terms %>%
  group_by(Description, genotype) %>%
  summarise(qvalue= qvalue)

# Spread the data into a wide format (GO terms x genotypes)
go_matrix_wide <- go_aggregated %>%
  spread(key = genotype, value = qvalue, fill = NA)

# Convert to a matrix for heatmap
go_matrix_qvalues <- as.matrix(go_matrix_wide[,-1])  # Exclude the Description column
rownames(go_matrix_qvalues) <- go_matrix_wide$Description

# Save as a CSV file:
write.csv(go_matrix_qvalues, file="final_analysis_new/midbrain_GOs_qvalue_matrix.csv")

# Plot the heatmap

# Apply log2 transformation to q-values
go_matrix_log2 <- -log2(go_matrix_qvalues)
go_matrix_log2[is.infinite(go_matrix_log2)] <- 0
go_matrix_log2[is.na(go_matrix_log2)] <- 0
go_matrix_log2[(go_matrix_log2) >= 12] <- 12 # for visual purposes

# Custom row and column labels
row_labels <- rownames(go_matrix_log2)
col_labels <- colnames(go_matrix_log2)

# Plot the heatmap without clustering
pdf("final_analysis_new/plots/midbrain_GOs_qvalue_heatmap.pdf", width=8, height= 10)
pheatmap(go_matrix_log2,
         color = colorRampPalette(c("white", "red"))(10),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         labels_row = row_labels,
         labels_col = col_labels)
dev.off()

pdf("final_analysis_new/plots/midbrain_GOs_qvalue_heatmap_colsclustered.pdf", width=8, height= 10)
pheatmap(go_matrix_log2,
         color = colorRampPalette(c("white", "red"))(10),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         labels_row = row_labels,
         labels_col = col_labels)
dev.off()

# Hindbrain:

# Extract rows that correspond to 'hindbrain'
hindbrain_rows <- grep("hindbrain", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'hindbrain' data
hindbrain.data <- subset(all.data.combined.selected.cor.pseudo, cells = hindbrain_rows)

# Create DESeq2 dataset
hindbrain.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(hindbrain.data@assays$SCT@counts + 1),
  colData = hindbrain.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
hindbrain.deseq2 <- estimateSizeFactors(hindbrain.deseq2)
hindbrain.deseq2 <- estimateDispersionsGeneEst(hindbrain.deseq2)
dispersions(hindbrain.deseq2) <- mcols(hindbrain.deseq2)$dispGeneEst
hindbrain.deseq2 <- nbinomWaldTest(hindbrain.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(hindbrain.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/hindbrain_DEGs_complete_results.csv", row.names = TRUE)

# spinal-cord:

# Extract rows that correspond to 'spinal-cord'
spinalcord_rows <- grep("spinal_cord", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'spinal-cord' data
spinalcord.data <- subset(all.data.combined.selected.cor.pseudo, cells = spinalcord_rows)

# Create DESeq2 dataset
spinalcord.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(spinalcord.data@assays$SCT@counts + 1),
  colData = spinalcord.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
spinalcord.deseq2 <- estimateSizeFactors(spinalcord.deseq2)
spinalcord.deseq2 <- estimateDispersionsGeneEst(spinalcord.deseq2)
dispersions(spinalcord.deseq2) <- mcols(spinalcord.deseq2)$dispGeneEst
spinalcord.deseq2 <- nbinomWaldTest(spinalcord.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(spinalcord.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/spinalcord_DEGs_complete_results.csv", row.names = TRUE)

# RGC:
# Extract rows that correspond to 'RGC'
RGC_rows <- grep("RGC", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'spinal-cord' data
RGC.data <- subset(all.data.combined.selected.cor.pseudo, cells = RGC_rows)

# Create DESeq2 dataset
RGC.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(RGC.data@assays$SCT@counts + 1),
  colData = RGC.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
RGC.deseq2 <- estimateSizeFactors(RGC.deseq2)
RGC.deseq2 <- estimateDispersionsGeneEst(RGC.deseq2)
dispersions(RGC.deseq2) <- mcols(RGC.deseq2)$dispGeneEst
RGC.deseq2 <- nbinomWaldTest(RGC.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(RGC.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/RGC_DEGs_complete_results.csv", row.names = TRUE)

# bipolar:
# Extract rows that correspond to 'bipolar'
bipolar_rows <- grep("bipolar", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'bipolar' data
bipolar.data <- subset(all.data.combined.selected.cor.pseudo, cells = bipolar_rows)

# Create DESeq2 dataset
bipolar.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(bipolar.data@assays$SCT@counts + 1),
  colData = bipolar.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
bipolar.deseq2 <- estimateSizeFactors(bipolar.deseq2)
bipolar.deseq2 <- estimateDispersionsGeneEst(bipolar.deseq2)
dispersions(bipolar.deseq2) <- mcols(bipolar.deseq2)$dispGeneEst
bipolar.deseq2 <- nbinomWaldTest(bipolar.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(bipolar.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/bipolar_DEGs_complete_results.csv", row.names = TRUE)

# amacrine-horizontal:
# Extract rows that correspond to 'amacrine-horizontal'
amacrine_horizontal_rows <- grep("amacrine_horizontal", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'amacrine_horizontal' data
amacrine_horizontal.data <- subset(all.data.combined.selected.cor.pseudo, cells = amacrine_horizontal_rows)

# Create DESeq2 dataset
amacrine_horizontal.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(amacrine_horizontal.data@assays$SCT@counts + 1),
  colData = amacrine_horizontal.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
amacrine_horizontal.deseq2 <- estimateSizeFactors(amacrine_horizontal.deseq2)
amacrine_horizontal.deseq2 <- estimateDispersionsGeneEst(amacrine_horizontal.deseq2)
dispersions(amacrine_horizontal.deseq2) <- mcols(amacrine_horizontal.deseq2)$dispGeneEst
amacrine_horizontal.deseq2 <- nbinomWaldTest(amacrine_horizontal.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(amacrine_horizontal.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/amacrine_horizontal_DEGs_complete_results.csv", row.names = TRUE)

# photoreceptors-precursors:

# Extract rows that correspond to 'photoreceptors-precursors'
photoreceptors_precursors_rows <- grep("photoreceptors_precursors", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'amacrine_horizontal' data
photoreceptor_precursors.data <- subset(all.data.combined.selected.cor.pseudo, cells = photoreceptors_precursors_rows)

# Create DESeq2 dataset
photoreceptor_precursors.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(photoreceptor_precursors.data@assays$SCT@counts + 1),
  colData = photoreceptor_precursors.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
photoreceptor_precursors.deseq2 <- estimateSizeFactors(photoreceptor_precursors.deseq2)
photoreceptor_precursors.deseq2 <- estimateDispersionsGeneEst(photoreceptor_precursors.deseq2)
dispersions(photoreceptor_precursors.deseq2) <- mcols(amacrine_horizontal.deseq2)$dispGeneEst
photoreceptor_precursors.deseq2 <- nbinomWaldTest(photoreceptor_precursors.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(photoreceptor_precursors.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/photoreceptors_precursors_DEGs_complete_results.csv", row.names = TRUE)

# cones:
# Extract rows that correspond to 'cones'
cones_rows <- grep("cones", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'cones' data
cones.data <- subset(all.data.combined.selected.cor.pseudo, cells = cones_rows)

# Create DESeq2 dataset
cones.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(cones.data@assays$SCT@counts + 1),
  colData = cones.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
cones.deseq2 <- estimateSizeFactors(cones.deseq2)
cones.deseq2 <- estimateDispersionsGeneEst(cones.deseq2)
dispersions(cones.deseq2) <- mcols(cones.deseq2)$dispGeneEst
cones.deseq2 <- nbinomWaldTest(cones.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(cones.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/cones_DEGs_complete_results.csv", row.names = TRUE)

# rods:
# Extract rows that correspond to 'rods'
rods_rows <- grep("rods", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'rods' data
rods.data <- subset(all.data.combined.selected.cor.pseudo, cells = rods_rows)

# Create DESeq2 dataset
rods.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(rods.data@assays$SCT@counts + 1),
  colData = rods.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
rods.deseq2 <- estimateSizeFactors(rods.deseq2)
rods.deseq2 <- estimateDispersionsGeneEst(rods.deseq2)
dispersions(rods.deseq2) <- mcols(rods.deseq2)$dispGeneEst
rods.deseq2 <- nbinomWaldTest(rods.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(rods.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/rods_DEGs_complete_results.csv", row.names = TRUE)

# radial-glia:
# Extract rows that correspond to 'radial-glia'
radialglia_rows <- grep("radial_glia", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'radial-glia' data
radialglia.data <- subset(all.data.combined.selected.cor.pseudo, cells = radialglia_rows)

# Create DESeq2 dataset
radialglia.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(radialglia.data@assays$SCT@counts + 1),
  colData = radialglia.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
radialglia.deseq2 <- estimateSizeFactors(radialglia.deseq2)
radialglia.deseq2 <- estimateDispersionsGeneEst(radialglia.deseq2)
dispersions(radialglia.deseq2) <- mcols(radialglia.deseq2)$dispGeneEst
radialglia.deseq2 <- nbinomWaldTest(radialglia.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(radialglia.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/radialglia_DEGs_complete_results.csv", row.names = TRUE)

# mueller-glia:
# Extract rows that correspond to 'mueller-glia'
muellerglia_rows <- grep("mueller_glia", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'mueller-glia' data
muellerglia.data <- subset(all.data.combined.selected.cor.pseudo, cells = muellerglia_rows)

# Create DESeq2 dataset
muellerglia.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(muellerglia.data@assays$SCT@counts + 1),
  colData = muellerglia.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
muellerglia.deseq2 <- estimateSizeFactors(muellerglia.deseq2)
muellerglia.deseq2 <- estimateDispersionsGeneEst(muellerglia.deseq2)
dispersions(muellerglia.deseq2) <- mcols(muellerglia.deseq2)$dispGeneEst
muellerglia.deseq2 <- nbinomWaldTest(muellerglia.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(muellerglia.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/muellerglia_DEGs_complete_results.csv", row.names = TRUE)

# microglia:
# Extract rows that correspond to 'mueller-glia'
microglia_rows <- grep("microglia", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'mueller-glia' data
microglia.data <- subset(all.data.combined.selected.cor.pseudo, cells = microglia_rows)

# Create DESeq2 dataset
microglia.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(microglia.data@assays$SCT@counts + 1),
  colData = microglia.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
microglia.deseq2 <- estimateSizeFactors(microglia.deseq2)
microglia.deseq2 <- estimateDispersionsGeneEst(microglia.deseq2)
dispersions(microglia.deseq2) <- mcols(microglia.deseq2)$dispGeneEst
microglia.deseq2 <- nbinomWaldTest(microglia.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(microglia.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/microglia_DEGs_complete_results.csv", row.names = TRUE)

# olfactory-organ:
# Extract rows that correspond to 'olfactory-organ'
olfactory_rows <- grep("olfactory_organ", rownames(all.data.combined.selected.cor.pseudo@meta.data), value = TRUE)

# Subset the Seurat object to include only 'olfactory-organ' data
olfactory.data <- subset(all.data.combined.selected.cor.pseudo, cells = olfactory_rows)

# Create DESeq2 dataset
olfactory.deseq2 <- DESeqDataSetFromMatrix(
  countData = round(olfactory.data@assays$SCT@counts + 1),
  colData = olfactory.data@meta.data,
  design = ~orig.ident
)

# Run DESeq2 analysis
olfactory.deseq2 <- estimateSizeFactors(olfactory.deseq2)
olfactory.deseq2 <- estimateDispersionsGeneEst(olfactory.deseq2)
dispersions(olfactory.deseq2) <- mcols(olfactory.deseq2)$dispGeneEst
olfactory.deseq2 <- nbinomWaldTest(olfactory.deseq2)

# Collect all DEGs:
deg_lists <- lapply(contrasts, function(contrast) {
  res <- as.data.frame(results(olfactory.deseq2, contrast = contrast, parallel = TRUE))
  res <- res[complete.cases(res), ]
  res$group <- paste(contrast, collapse = "_vs_")
  return(res)
})
all_deg_results <- do.call(rbind, deg_lists)
write.csv(all_deg_results, file = "final_analysis_new/olfactory_DEGs_complete_results.csv", row.names = TRUE)

```

# Make heatmap across genotypes to highlight genes of interest:
```{r}
library(RColorBrewer)

# Heatmap to show total counts of DEGs:
DEGs_up <- read.csv("DEGs_counts_up.csv")
DEGs_down <- read.csv("DEGs_counts_down.csv")
rownames(DEGs_up) <- DEGs_up$Cell.type
rownames(DEGs_down) <- DEGs_down$Cell.type
DEGs_up <- DEGs_up[,-1]
DEGs_down <- DEGs_down[,-1]

DEGs_up <- as.matrix(DEGs_up)
DEGs_down <- as.matrix(DEGs_down)

DEGs_up <- scale(DEGs_up, center= T, scale= T)
DEGs_down <- scale(DEGs_down, center= T, scale= T)

pdf("DEGs_counts_up.pdf", height = 5, width=5)
x <- pheatmap(DEGs_up, 
         scale="column", 
         cluster_rows=FALSE, 
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "red"))(100),
         breaks = c(0,1,2,3,4),
         na_col= "white")
dev.off()

pdf("DEGs_counts_down.pdf", height = 5, width=5)
pheatmap(DEGs_down, 
         scale="column", 
         cluster_rows=FALSE, 
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "blue"))(100),
         breaks = c(0,1,2,3,4),
         na_col = "white")
dev.off()

# Make zoomed in heatmap for regions of interest:

# Forebrain:
res.list <- list(
  res.forebrain.srgap2KO,
  res.forebrain.SRGAP2C,
  res.forebrain.gpr89KO,
  res.forebrain.GPR89B,
  res.forebrain.fam72KO,
  res.forebrain.FAM72B,
  res.forebrain.frmpd2KO,
  res.forebrain.FRMPD2B,
  res.forebrain.hydinKO,
  res.forebrain.pdzk1KO,
  res.forebrain.PDZK1P1
)

# Function to filter genes
filter_genes <- function(res) {
  res$gene <- rownames(res)
  DEGs <- res[res$padj < 0.2, "gene"]
  return(DEGs)
}

DEGs_genes_list <- lapply(res.list, filter_genes)
common_rownames <- Reduce(intersect, DEGs_genes_list)
cleaned_common_rownames <- common_rownames[!grepl("^NA\\.", common_rownames)]
cleaned_common_rownames <- cleaned_common_rownames[complete.cases(cleaned_common_rownames)]

# Heatmap of log2FC:
log2fc_matrix <- data.frame(row.names = cleaned_common_rownames)
for (i in 1:length(res.list)) {
  df <- res.list[[i]]
  df <- df[cleaned_common_rownames, , drop=FALSE]
  log2fc_matrix <- cbind(log2fc_matrix, df$log2FoldChange)
}
write.csv(log2fc_matrix, file="final_analysis/data/Forebrain_log2fc.csv")

# Set column names
colnames(log2fc_matrix) <- c("srgap2KO", "SRGAP2C", "gpr89KO", "GPR89B", "fam72KO", "FAM72B", "frmpd2KO", "FRMPD2B","hydinKO", "pdzk1KO", "PDZK1P1")

# Convert to matrix
log2fc_matrix <- log2fc_matrix[complete.cases(log2fc_matrix),] # 96 genes
log2fc_matrix <- as.matrix(log2fc_matrix)
log2fc_matrix2 <- scale(log2fc_matrix, center = T, scale= T)

pdf("final_analysis/plots/Forebrain_heatmap.pdf", height = 4, width= 4)
pheatmap(log2fc_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         display_numbers = FALSE,
         drop_levels = TRUE,
         border_color = NA,
         color = colorRampPalette(c("blue", "white", "red"))(100), fontsize_row = 0.1
         )
dev.off()

# Check GO terms of top genes:
background_genes_to_keep <- apply(forebrain.data@assays$SCT$counts, 1, function(row) any(row != 0))
background_genes <-rownames(count_matrix[background_genes_to_keep,])

forebrain_ego <- enrichGO(gene= cleaned_common_rownames,
                          universe= background_genes,
                          OrgDb= org.Dr.eg.db,
                          keyType = "SYMBOL",
                          ont= "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.1,
                          qvalueCutoff = 0.1,
                          readable = T)
forebrain_ego <- simplify(forebrain_ego, cutoff= 0.5)
pdf("final_analysis/plots/Forebrain_GO_dotplot.pdf", width = 5, height= 5)
dotplot(forebrain_ego)
dev.off()
write.csv(forebrain_ego@result, file="final_analysis/plots/Forebrain_GO.csv")

# Network:
pdf("final_analysis/plots/Forebrain_GO_network.pdf")
cnetplot(forebrain_ego)
dev.off()

# Midbrain:
res.list <- list(
  res.midbrain.srgap2KO,
  res.midbrain.SRGAP2C,
  res.midbrain.fam72KO,
  res.midbrain.FAM72B,
  res.midbrain.frmpd2KO,
  res.midbrain.FRMPD2B,
  res.midbrain.hydinKO
)

# Function to filter genes
filter_genes <- function(res) {
  res$gene <- rownames(res)
  DEGs <- res[res$padj < 0.2, "gene"]
  return(DEGs)
}

DEGs_genes_list <- lapply(res.list, filter_genes)
common_rownames <- Reduce(intersect, DEGs_genes_list)
cleaned_common_rownames <- common_rownames[!grepl("^NA\\.", common_rownames)]
cleaned_common_rownames <- cleaned_common_rownames[complete.cases(cleaned_common_rownames)] # 211 genes

# Heatmap of log2FC:
log2fc_matrix <- data.frame(row.names = cleaned_common_rownames)
for (i in 1:length(res.list)) {
  df <- res.list[[i]]
  df <- df[cleaned_common_rownames, , drop=FALSE]
  log2fc_matrix <- cbind(log2fc_matrix, df$log2FoldChange)
}
write.csv(log2fc_matrix, file="final_analysis/data/FMidbrain_log2fc.csv")

# Set column names
colnames(log2fc_matrix) <- c("srgap2KO", "SRGAP2C", "fam72KO", "FAM72B", "frmpd2KO", "FRMPD2B","hydinKO")

# Convert to matrix
log2fc_matrix <- log2fc_matrix[complete.cases(log2fc_matrix),] # 96 genes
log2fc_matrix <- as.matrix(log2fc_matrix)
log2fc_matrix2 <- scale(log2fc_matrix, center = T, scale= T)

pdf("final_analysis/plots/Midbrain_heatmap.pdf", height = 4, width= 4)
pheatmap(log2fc_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         display_numbers = FALSE,
         drop_levels = TRUE,
         border_color = NA,
         color = colorRampPalette(c("blue", "white", "red"))(100), fontsize_row = 0.1
         )
dev.off()

# Check GO terms of top genes:
background_genes_to_keep <- apply(midbrain.data@assays$SCT$counts, 1, function(row) any(row != 0))
background_genes <-rownames(count_matrix[background_genes_to_keep,])

midbrain_ego <- enrichGO(gene= cleaned_common_rownames,
                          universe= background_genes,
                          OrgDb= org.Dr.eg.db,
                          keyType = "SYMBOL",
                          ont= "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.1,
                          qvalueCutoff = 0.1,
                          readable = T)
midbrain_ego <- simplify(midbrain_ego, cutoff= 0.5)
pdf("final_analysis/plots/Midbrain_GO_dotplot.pdf", width = 5, height= 5)
dotplot(midbrain_ego)
dev.off()
write.csv(midbrain_ego@result, file="final_analysis/plots/Midbrain_GO.csv")

# Network:
pdf("final_analysis/plots/Midbrain_GO_network.pdf")
cnetplot(midbrain_ego)
dev.off()

```

# Get the number of expressed genes per cell type.
```{r}
# Extract metadata and expression data
metadata <- all.data.combined.selected.filtered@meta.data
expression_data <- all.data.combined.selected.filtered@assays$SCT@counts

# Initialize a list to store results:
results_list <- list()

# Celltype levels:
cell_types <- unique(metadata$specific.celltypes.final)

for (cell_type in cell_types) {
  cells_of_type <- which(metadata$specific.celltypes.final == cell_type)
  
  expression_subset <- expression_data[, cells_of_type, drop = FALSE]
  
  num_expressed_genes_per_cell <- Matrix::rowSums(expression_subset)
  total_expressed_genes <- length(num_expressed_genes_per_cell[num_expressed_genes_per_cell>1])
    results_list[[cell_type]] <- total_expressed_genes
}

# Convert to a data frame
results_df <- data.frame(
  celltype = names(results_list),
  total_expressed_genes = unlist(results_list)
)

# Plot:

results_df$celltype <- factor(results_df$celltype, levels=c("forebrain","pallium", "hypothalamus","thalamus", "olfactory_organ", "midbrain", "hindbrain", "spinal_cord", "RGC", "bipolar_cells", "amacrine_horizontal", "photoreceptors_precursors", "cones","rods","radial_glia", "mueller_glia", "myeloid"))

p_gene_counts <- ggplot(data= results_df, 
                           aes(x = celltype, y = 1, size = total_expressed_genes)) +
  geom_point(alpha = 0.8) +
  theme_minimal()

ggsave("final_analysis/plots/balloonplot_cluster_number_expressed_genes.pdf",
       p_gene_counts, width= 6, height= 4)

```
